module molecular_mod

 ! written by tjh
 ! 21/11/06

 ! made my own by dar
 ! 01/X/07

   use kind_mod
   use constants_mod
   use utils_mod
   use messages_mod
   use timing
!   use grid_mod
   use gridio_mod
   use math_mod
   use datacube_mod
   use nrutil
   use nrtype
   use parallel_mod
   use vtk_mod
!   use mkl_vsl_type
!   use mkl_vsl

   implicit none

   integer :: mintrans, minlevel, maxlevel, maxtrans, nlevels
   integer :: grand_iter, nray
   Logical :: molebench

!   include '/Library/Frameworks/Intel_MKL.framework/Headers/mkl_vml.fi'

 ! Define data structures used within code - MOLECULETYPE holds parameters unique to a particular molecule
 ! Telescope holds data about particular telescopes

 ! Much of this code is our interpretation of the 2000 paper by F. v.d. Tak and M. Hogerheijde - and the 
 ! 2002 paper by G-J. v. Zadelhoff et al

   type MOLECULETYPE

      character(len=10) :: molecule
      real :: molecularweight
      real(double) :: abundance
      integer :: nLevels
      real(double), pointer :: energy(:)
      real(double), pointer :: g(:)
      real(double), pointer :: j(:)
      integer :: nTrans
      real(double), allocatable :: einsteinA(:)
      real(double), allocatable :: einsteinBlu(:)
      real(double), allocatable :: einsteinBul(:)
!      real(double), pointer:: einsteinBul(:)
!      real(double), pointer:: einsteinBlu(:)
!      real(double), pointer:: einsteinA(:)
      real(double), pointer :: transfreq(:)
      integer, pointer :: itransUpper(:)
      integer, pointer :: itransLower(:)
      real(double), pointer :: Eu(:)
      integer, pointer :: iCollUpper(:,:)
      integer, pointer :: iCollLower(:,:)
      integer :: nCollPart
      character(len=20), pointer :: collBetween(:)
      integer, pointer :: nCollTrans(:)
      integer, pointer :: nCollTemps(:)
      real(double), pointer :: collTemps(:,:)
      real(double), pointer :: collRates(:,:,:)
   end type MOLECULETYPE
   
 contains
   ! Read in molecular parameters from file - note: abundance hard-coded here

   subroutine readMolecule(thisMolecule, molFilename)
     use input_variables, only : molAbundance
     type(MOLECULETYPE) :: thisMolecule
     character(len=*) :: molFilename
     character(len=80) :: junk
     character(len=200):: dataDirectory, filename
     integer :: i, j, iLow, iUp, iPart
     real(double) :: a, freq, eu, c(20)
     logical :: preprocess1 = .true.
     logical :: preprocess2 = .true.
     integer :: maxnCollTrans, maxnCollTemps

     !thisMolecule%abundance = tiny(thisMolecule%abundance)
     thisMolecule%abundance = molAbundance ! fixed at benchmark value here

     do while(preprocess2)
        if(.not. preprocess1) preprocess2 = .false.

        call unixGetenv("TORUS_DATA", dataDirectory)
        filename = trim(dataDirectory)//"/"//molfilename
        
        open(30, file=filename, status="old", form="formatted")

        read(30,*) junk
        read(30,'(a)') thisMolecule%molecule
        if( preprocess2) call writeInfo("Reading data for molecule: "//trim(thisMolecule%molecule),IMPORTANT)

        read(30,*) junk
        read(30,*) thisMolecule%molecularWeight
        
        read(30,*) junk
        read(30,*) thisMolecule%nLevels
        
        allocate(thisMolecule%energy(1:thisMolecule%nLevels))
        allocate(thisMolecule%g(1:thisMolecule%nLevels))
        allocate(thisMolecule%j(1:thisMolecule%nLevels))

        read(30,*) junk
        do i = 1, thisMolecule%nLevels
           read(30,*) j, thisMolecule%energy(i), thisMolecule%g(i), thisMolecule%j(i)
           thisMolecule%energy(i) = thisMolecule%energy(i) / 8065.541  ! convert from per cm to ev - e/hc (CGS)
           
        enddo
        
        read(30,*) junk
        read(30,*) thisMolecule%nTrans
                
        allocate(thisMolecule%einsteinA(1:thisMolecule%nTrans))
        allocate(thisMolecule%einsteinBul(1:thisMolecule%nTrans))
        allocate(thisMolecule%einsteinBlu(1:thisMolecule%nTrans))
        allocate(thisMolecule%transfreq(1:thisMolecule%nTrans))
        allocate(thisMolecule%itransUpper(1:thisMolecule%nTrans))
        allocate(thisMolecule%itransLower(1:thisMolecule%nTrans))
        allocate(thisMolecule%eu(1:thisMolecule%nTrans))

        read(30,*) junk
        do i = 1, thisMolecule%nTrans
           read(30,*) j, iUp, iLow, a, freq, eu
           
           thisMolecule%einsteinA(i) = a
           thisMolecule%transfreq(i) = freq*1.d9
           thisMolecule%eu(i) = eu
           thisMolecule%itransUpper(i) = iUp
           thisMolecule%itransLower(i) = iLow
           
           thisMolecule%einsteinBul(i) = a * (cspeed**2)/(2.d0*hcgs*(freq*1.d9)**3) ! transform Aul -> Bul
           
           thisMolecule%einsteinBlu(i) = thisMolecule%einsteinBul(i) &
                * thisMolecule%g(iUp)/thisMolecule%g(iLow)

        enddo

        read(30,*) junk
        read(30,*) thisMolecule%nCollPart
        
        if(.not. preprocess2) then
           maxnCollTrans = maxval(thisMolecule%nCollTrans(:))
           maxnCollTemps = maxval(thisMolecule%nCollTemps(:))
           deallocate(thisMolecule%nCollTrans, thisMolecule%nCollTemps)
        else
           maxnCollTrans = 5000
           maxnCollTemps = 20
        endif

        allocate(thisMolecule%nCollTrans(1:thisMolecule%nCollPart))
        allocate(thisMolecule%nCollTemps(1:thisMolecule%nCollPart))
        allocate(thisMolecule%collTemps(1:maxnCollTrans, 1:maxnCollTemps))
        allocate(thisMolecule%collBetween(1:thisMolecule%nCollPart))

        allocate(thisMolecule%collRates(1:thisMolecule%nCollPart, maxnCollTrans, &
             maxnCollTemps))
        allocate(thisMolecule%iCollUpper(thisMolecule%nCollPart, maxnCollTrans))
        allocate(thisMolecule%iCollLower(thisMolecule%nCollPart, maxnCollTrans))

        do iPart = 1, thisMolecule%nCollPart

           read(30,*) junk
           read(30,*) thisMolecule%collBetween(iPart)

           read(30,*) junk
           read(30,*) thisMolecule%nCollTrans(iPart)

           read(30,*) junk
           read(30,*) thisMolecule%nCollTemps(iPart)

           read(30,*) junk
           read(30,*) thisMolecule%collTemps(iPart,1:thisMolecule%nCollTemps(ipart))

           thisMolecule%collRates = 0.d0
           read(30,*) junk

           do j = 1, thisMolecule%nCollTrans(iPart)

              read(30,*) i, thisMolecule%iCollUpper(ipart,j), &
                   thisMolecule%iCollLower(ipart,j), c(1:thisMolecule%nCollTemps(iPart))
              thisMolecule%collRates(iPart, j, 1:thisMolecule%nCollTemps(iPart)) = c(1:thisMolecule%nCollTemps(iPart))

           enddo
        enddo
        close(30)

        if(preprocess1) then
           deallocate(thisMolecule%energy, thisMolecule%g, thisMolecule%j, thisMolecule%einsteinA, thisMolecule%einsteinBlu, &
                      thisMolecule%einsteinBul, thisMolecule%transfreq, thisMolecule%itransUpper, thisMolecule%itransLower, &
                      thisMolecule%Eu, thisMolecule%iCollUpper, thisMolecule%iCollLower, thisMolecule%collBetween, &
                      thisMolecule%collTemps, thisMolecule%collRates)
        endif

        preprocess1 = .false.
     enddo

     call writeInfo("Done.", IMPORTANT)
   end subroutine readMolecule

 ! Same routine for benchmark molecule (slightly different file format)
   subroutine readBenchmarkMolecule(thisMolecule, molFilename)
     use input_variables, only : molAbundance
     type(MOLECULETYPE) :: thisMolecule
     character(len=*) :: molFilename
     character(len=200):: dataDirectory, filename
     character(len=80) :: junk
     integer :: i, iPart

 !    thisMolecule%abundance = tiny(thisMolecule%abundance)
     thisMolecule%abundance = molAbundance!hcoabundance ! fixed at benchmark value here

     call unixGetenv("TORUS_DATA", dataDirectory)
     filename = trim(dataDirectory)//"/"//molfilename
     
     call writeInfo("Opening file: "//trim(filename),TRIVIAL)
     open(30, file=filename, status="old", form="formatted")

     read(30,'(a)') thisMolecule%molecule
     call writeInfo("Reading data for molecule: "//trim(thisMolecule%molecule),IMPORTANT)

     read(30,*) thisMolecule%molecularWeight

     read(30,*) thisMolecule%nLevels, thisMolecule%nTrans

     allocate(thisMolecule%energy(1:thisMolecule%nLevels))
     allocate(thisMolecule%g(1:thisMolecule%nLevels))
     allocate(thisMolecule%j(1:thisMolecule%nLevels))

     allocate(thisMolecule%einsteinA(1:thisMolecule%nTrans))
     allocate(thisMolecule%einsteinBul(1:thisMolecule%nTrans))
     allocate(thisMolecule%einsteinBlu(1:thisMolecule%nTrans))
     allocate(thisMolecule%transfreq(1:thisMolecule%nTrans))
     allocate(thisMolecule%itransUpper(1:thisMolecule%nTrans))
     allocate(thisMolecule%itransLower(1:thisMolecule%nTrans))
     allocate(thisMolecule%eu(1:thisMolecule%nTrans))

     read(30,'(7f11.7)') thisMolecule%energy(1:thisMolecule%nLevels)

     thisMolecule%energy = thisMolecule%energy / 8065.541  ! per cm to ev

     read(30,*) thisMolecule%g(1:thisMolecule%nLevels)
     thisMolecule%j(1:thisMolecule%nLevels) = (thisMolecule%g(1:thisMolecule%nLevels) - 1.d0)*0.5d0

     read(30,*) thisMolecule%itransUpper(1:thisMolecule%nTrans)
     read(30,*) thisMolecule%itransLower(1:thisMolecule%nTrans)

     read(30,*) thisMolecule%einsteinA(1:thisMolecule%nTrans)

     do i = 1, thisMolecule%nTrans
        thisMolecule%transFreq(i) = ((thisMolecule%energy(thisMolecule%iTransUpper(i)) - &
             thisMolecule%energy(thisMolecule%iTransLower(i)))/ergToEV)/hcgs
        thisMolecule%einsteinBul(i) = thisMolecule%einsteinA(i) * (cspeed**2)/(2.d0*hcgs*thisMolecule%transFreq(i)**3)
        thisMolecule%einsteinBlu(i) = thisMolecule%einsteinBul(i) &
             * thisMolecule%g(thisMolecule%iTransUpper(i))/thisMolecule%g(thisMolecule%iTransLower(i))
     enddo

     read(30,*) junk

     thisMolecule%nCollPart = 1 ! only one collision partner considered in benchmark
     iPart = 1 ! follows from ^

     allocate(thisMolecule%nCollTrans(1:thisMolecule%nCollPart))
     allocate(thisMolecule%nCollTemps(1:thisMolecule%nCollPart))
     allocate(thisMolecule%collTemps(1:thisMolecule%nCollPart, 1:4))

     read(30,*) thisMolecule%nCollTrans(1), thisMolecule%nCollTemps(1),thisMolecule%collTemps(1,1:4)

     allocate(thisMolecule%iCollUpper(1:thisMolecule%nCollPart, 1:thisMolecule%nCollTrans(iPart)))
     allocate(thisMolecule%iCollLower(1:thisMolecule%nCollPart, 1:thisMolecule%nCollTrans(iPart)))
     allocate(thisMolecule%collRates(1:thisMolecule%nCollPart, &
          1:thisMolecule%nCollTrans(iPart), 1:thisMolecule%nCollTemps(ipart)))

     thisMolecule%collrates = -1.d0

     read(30,*) thisMolecule%iCollUpper(1, 1:thisMolecule%nCollTrans(1))
     read(30,*) thisMolecule%iCollLower(1, 1:thisMolecule%nCollTrans(1))

     do i = 1, thisMolecule%nCollTemps(1)
        read(30,*) thisMolecule%collRates(iPart,  1:thisMolecule%nCollTrans(ipart), i)
     enddo

     close(30)
     call writeInfo("Done.", IMPORTANT)
   end subroutine readBenchmarkMolecule

 ! This subroutine allocates memory to the octal for storing molecular level data *if* none currently exists then
 ! stores the non-zero local emission coefficient   
   recursive subroutine  allocateMolecularLevels(grid, thisOctal, thisMolecule, restart, isinlte)

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i
     logical, optional :: restart, isinlte
     logical :: dorestart, lte

     if(.not. present(restart)) then
        dorestart = .false.
     else
        dorestart = restart
     endif

     if(.not. present(isinlte)) then
        lte = .false.
     else
        lte = isinlte
     endif

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren, 1
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call allocateMolecularLevels(grid, child, thisMolecule, restart, isinlte)
                 exit
              end if
           end do
        else

           if(dorestart) then

              allocate(thisOctal%newmolecularLevel(1:thisOctal%maxChildren,1:maxlevel))
              thisOctal%newmolecularLevel = 1.d-20
              allocate(thisOctal%oldmolecularLevel(1:thisOctal%maxChildren,1:maxlevel))
              thisOctal%oldmolecularLevel = 1.d-20

              if (.not.associated(thisOctal%bnu)) then
                 allocate(thisOctal%bnu(1:thisOctal%maxChildren, 1:maxlevel))
              endif
              thisOctal%bnu(subcell,i) = bnu(thisMolecule%transFreq(i), dble(thisOctal%temperature(subcell)))

           else

              if (.not. associated(thisOctal%newmolecularLevel)) then
                 allocate(thisOctal%newmolecularLevel(1:thisOctal%maxChildren,1:maxlevel))
              endif
              thisOctal%newmolecularLevel = 1.d-20

              if (.not.associated(thisOctal%oldmolecularLevel)) then
                 allocate(thisOctal%oldmolecularLevel(1:thisOctal%maxChildren,1:maxlevel))
              endif
              thisOctal%oldmolecularLevel = 1.d-20

              if (.not.associated(thisOctal%molecularLevel)) then
                 allocate(thisOctal%molecularLevel(1:thisOctal%maxChildren,1:maxlevel))
              endif
              thisOctal%molecularLevel = 1.d-20

              if(lte) then           
                 call LTEpops(thisMolecule, dble(thisOctal%temperature(subcell)), &
                              thisOctal%molecularLevel(subcell,1:maxlevel))
              else      
                 call LTEpops(thisMolecule, tcbr, thisOctal%molecularLevel(subcell,1:maxlevel))
              endif

              if((grid%geometry .eq. "h2obench1") .or. (grid%geometry .eq. "h2obench2")) then
                 thisOctal%molecularLevel(subcell,1) = 1.0
                 thisOctal%molecularLevel(subcell,2) = 0.0000
              endif

              if (.not.associated(thisOctal%bnu)) then
                 allocate(thisOctal%bnu(thisOctal%maxChildren, 1:maxtrans))
              endif

              do i = 1, maxtrans
                 thisOctal%bnu(subcell,i) = bnu(thisMolecule%transFreq(i), dble(thisOctal%temperature(subcell)))
              enddo

              if (.not.associated(thisOctal%jnu)) then
                 allocate(thisOctal%jnu(1:thisOctal%maxChildren, 1:maxtrans))
              endif
              thisOctal%oldmolecularLevel = 1.d-20

              if(lte) then
                 do i = 1, maxtrans
                    thisOctal%jnu(subcell,i) = bnu(thisMolecule%transFreq(i), dble(thisoctal%temperature(subcell)))
 !                if(thisOctal%temperature(subcell) .gt. Tcbr) then
!                    thisOctal%jnu(subcell,i) = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell) * &
!                                               thisOctal%molecularLevel(subcell,i) * thisMolecule%einsteinA(i) * hCgsOverfourPi * &
!                                              thisOctal%microturb(subcell) / sqrtPi
!                    write(*,*) thisoctal%jnu(subcell,i)
                  enddo
              else
                 do i = 1, maxtrans
                    thisOctal%jnu(subcell,i) = bnu(thisMolecule%transFreq(i), tcbr)
                 enddo
              endif
              thisOctal%oldmolecularLevel = 1.d-20
           endif

           if (.not.associated(thisOctal%molAbundance)) then
              allocate(thisOctal%molAbundance(1:thisOctal%maxChildren))
              thisOctal%molAbundance(subcell) = thisMolecule%abundance
           endif

          if(thisOctal%temperaturedust(subcell) .eq. thisOctal%temperaturegas(subcell)) then
             thisOctal%temperaturegas(subcell) = thisOctal%temperature(subcell)
             thisOctal%temperaturedust(subcell) = thisOctal%temperature(subcell)
          endif

          thisOctal%molmicroturb(subcell) = 1.d0 / thisOctal%microturb(subcell)
!          if (grid%geometry .eq. "molcluster") CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
       endif
    enddo

!    if (.not. usedust) then
       grid%oneKappaAbs = 1e-30
       grid%oneKappaSca = 1e-30
!    endif

   end subroutine allocateMolecularLevels

 ! Does a lot of work - do more rays whilst problem not converged -            
   subroutine molecularLoop(grid, thisMolecule)

     use input_variables, only : blockhandout, tolerance, debug, geometry, lucyfilenamein, openlucy,&
          usedust, rinner, router, readmol, amr1d, amr2d, amr3d, plotlevels
     use messages_mod, only : myRankIsZero
#ifdef MPI
     include 'mpif.h'
#endif

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(VECTOR) :: position, direction
     integer :: nOctal, iOctal, subcell
     real(double), allocatable :: ds(:), phi(:), i0(:,:)

     type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
     type(OCTAL), pointer :: thisOctal
     integer, parameter :: maxIter = 1000, maxRay = 1000000
     logical :: popsConverged, gridConverged, gridConvergedTest
     character(len=200) :: message
     integer :: iRay, iTrans, iter, i 
     integer :: iStage
     real(double), allocatable :: oldpops(:)
     real(double) :: fac
#ifdef MPI
     ! For MPI implementations
     integer       ::   my_rank        ! my processor rank
     integer       ::   np             ! The number of processes
     integer       ::   ierr           ! error flag
     integer       ::   tag = 0 
     logical       ::   rankComplete
     integer, dimension(:), allocatable :: octalsBelongRank
     real(double), allocatable :: tArrayd(:),tempArrayd(:) 
#endif
      integer :: nVoxels
      integer :: ioctal_beg, ioctal_end
      logical :: fixedRays
      integer :: isize
      integer, allocatable :: iseed(:)
      real(double) :: maxFracChange 
      real(double), allocatable :: maxFracChangePerLevel(:), avgFracChange(:,:)
      integer, allocatable :: convergenceCounter(:,:)
 !     real :: r1(2) !random numbers
      real(double) ::  maxavgfracChange, maxRMSfracChange
      integer :: maxavgtrans(1),maxRMStrans(1)

      character(len=30) :: filename

!      real(double) :: convtestarray(200,200,16)
      real(double) :: r(100)

      logical :: restart = .false.
      logical :: firsttime = .true.
      
      integer :: ntrans , lst !least significant transition

!      integer, allocatable, save :: ioctalArray(:)
      
      real :: tol
      real(double) :: dummy, tau, tauarray(40) = -1.d0    
      logical :: quasi = .true.

 ! blockhandout must be off for fixed ray case, otherwise setting the
 ! seed is not enough to ensure the same directions are done for
 ! each cell every iteration

      nlevels = thisMolecule%nlevels
      if(grid%geometry .eq. 'molebench') molebench = .true.

      call countVoxels(grid%octreeRoot,nOctal,nVoxels)

!     finds max temperature and works out highest inhabited level at a given fraction. Highest necessary level is 4 above that
      call findmaxlevel(grid, grid%octreeroot, thisMolecule, maxlevel, nlevels, nVoxels)

      maxlevel = max(maxlevel,minlevel)
      maxtrans = maxlevel - 1
      
      write(message, *) "Maximum Interesting Level", maxlevel
      call writeinfo(message, FORINFO)

      blockHandout = .false. 

#ifdef MPI
     ! FOR MPI IMPLEMENTATION=======================================================
     !  Get my process rank # 
     call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

     ! Find the total # of precessor being used in this run
     call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#endif

      minlevel = maxlevel
      mintrans = maxtrans

 ! Write Column headings for file containing level populations and convergence data
     grand_iter = 0

     if(writeoutput) then
18      format(a,tr3,8(a,tr6),a,6(tr1,a))
        open(138,file="fracChanges.dat",status="unknown",form="formatted")
        write(138,18) "nRays","J=0","J=1","J=2","J=3","J=4","J=5","J=6","J=7","Fixed","1% Conv","2.5% Conv","5% Conv","Max","Avg"
        close(138)

        open(139,file="avgChange.dat",status="replace",form="formatted")
        open(140,file="avgRMSChange.dat",status="replace",form="formatted")
        open(141,file="tau.dat",status="replace",form="formatted")
        close(139)
        close(140)
        close(141)
     endif

      if(openlucy) then

         if(writeoutput) then
            write(message,*) "Reading in lucy temp files: ",lucyfilenamein
            call writeinfo(message,FORINFO)
         endif

         call freeGrid(grid)
         call readAmrGrid(lucyfilenamein,.false.,grid)

         call writeinfo("Successfully read in Lucy Grid file", FORINFO)
         call writeinfo("Plotting Temperatures...",TRIVIAL)
      endif

      allocate(octalArray(grid%nOctals))
      call getOctalArray(grid%octreeRoot,octalArray, nOctal)

      if(restart) then

         call writeinfo("Reading in previous grid", FORINFO)
         call freeGrid(grid)
         call readAMRgrid("molecular_tmp.grid",.false.,grid)
         call allocateMolecularLevels(grid, grid%octreeRoot, thisMolecule, .true., .true.)
         call writeinfo("Done!", TRIVIAL)

      else
         
         call writeinfo("Allocating and initialising molecular levels", FORINFO)
         call AllocateMolecularLevels(grid, grid%octreeRoot, thisMolecule, .false., .true.)
         call writeinfo("Done!", TRIVIAL)
         
         if(myrankiszero) call writeAMRgrid("molecular_lte.grid",.false.,grid)
      endif

      nRay = 1 ! number of rays used to establish estimate of jnu and pops

      if((amr1d .or. amr2d) .and. writeoutput) call dumpresults(grid, thisMolecule)!, convtestarray) ! find radial pops on final grid     

      if(usedust) then 
         call continuumIntensityAlongRay(vector(1d10,0.d0,0.d0),vector(-1.d0,1d-20,1d-20), grid, 1e4, 0.d0, dummy, &
              tau, .true.)
         if(writeoutput) write(*,*) "TAU", tau
      endif

      nRay = 100 ! number of rays used to establish estimate of jnu and pops

      allocate(oldPops(1:maxtrans))

      call init_random_seed()

     !       if (myRankIsZero) &
 !    call writeAmrGrid("molecular_tmp.grid",.false.,grid)

      call random_seed(size=iSize)
      allocate(iSeed(1:iSize))
      call random_seed(get=iSeed)
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      call MPI_BCAST(iSeed, iSize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
      
      do iStage = 1, 2 ! fixed rays to reduce variance between cells or 2)  random rays to ensure sufficient spatial sampling
         
         if (iStage == 1) then
            fixedRays = .true.
            nRay = 100
         else
            fixedRays = .false.
            nRay = 100
         endif

         gridConvergedTest = .false.
         gridConverged = .false.

         do while (.not.gridConverged)

            grand_iter = grand_iter + 1

             call taualongray(VECTOR(1.d-10,1.d-10,1.d-10), VECTOR(1.d0, -1.d-20, -1.d-20), &
	     grid, thisMolecule, 0.d0, tauarray(1:maxtrans))

!            write(*,*) "tauarray", tauarray

            do i = size(tauarray), 1, -1
               mintrans = i + 2
               minlevel = mintrans + 1
!               write(*,*) i, tauarray(i)
               if(tauarray(i) .gt. 0.01) exit
            enddo
            
            write(message, *) "Minimum important level", minlevel
            call writeinfo(message, TRIVIAL)

            write(message,*) "Iteration ",grand_iter
            call writeinfo(message, FORINFO)

            if(writeoutput) then
               write(message,*) "Done ",nray," rays"
               call tune(6, message)  ! start a stopwatch
            endif

            if(Writeoutput .and. plotlevels .and. .not.(amr1d)) then
               write(filename, *) "./plots/data_",grand_iter,".vtk"
               call  writeVtkFile(grid, filename)
            endif

            allocate(ds(1:nRay))
            allocate(phi(1:nRay))
            allocate(i0(1:maxtrans, 1:nRay))

            if (fixedRays) then
               call random_seed(put=iseed)   ! same seed for fixed rays
            else
               call init_random_seed()
            endif
            
            ! default loop indicies
            ioctal_beg = 1
            ioctal_end = SIZE(octalArray)         

#ifdef MPI

     ! we will use an array to store the rank of the process
     !   which will calculate each octal's variables
     allocate(octalsBelongRank(size(octalArray)))

     if (my_rank == 0) then
        print *, ' '
        print *, 'molecular  computed by ', np-1, ' processors.'
        print *, ' '
        call mpiBlockHandout(np,octalsBelongRank,blockDivFactor=1,tag=tag,&
                             maxBlockSize=10,setDebug=.false.)

     endif
     ! ============================================================================

  if (my_rank /= 0) then
   blockLoop: do     
  call mpiGetBlock(my_rank,iOctal_beg,iOctal_end,rankComplete,tag,setDebug=.false.)
    if (rankComplete) exit blockLoop 
#endif

!    if(firsttime) then
!       allocate(iOctalArray(ioctal_end - ioctal_beg + 1))
!       call decideOctalOrder(iOctal_beg, iOctal_end, iOctalArray, shouldinterlace = .true.)
!       firsttime = .false.
!    endif


 ! iterate over all octals, all rays, solving the system self-consistently
           if(writeoutput) then
              write(message,*) "Getray"
              call tune(6, message)  ! start a stopwatch
           endif

           do iOctal = ioctal_beg, ioctal_end

!              if (debug .and. writeoutput) then
!                 write(message,*) iOctal,ioctal_beg,ioctal_end,ioctalArray(ioctal)
!                 call writeInfo(message,TRIVIAL)
!              endif

 !             thisOctal => octalArray(iOctalArray(iOctal))%content
             thisOctal => octalArray(ioctal)%content

              do subcell = 1, thisOctal%maxChildren

                 if (.not.thisOctal%hasChild(subcell)) then
!                    if(fixedrays) call sobseq(r1,-1)                     

                    do iRay = 1, nRay


                       call getRay(grid, thisOctal, subcell, position, direction, &
                            ds(iRay), phi(iRay), i0(1:maxtrans,iRay), &
                            thisMolecule,fixedrays) ! does the hard work - populates i0 etc

                    enddo

                    iter = 0
                    popsConverged = .false.

                    thisOctal%oldMolecularLevel(subcell,:) = thisOctal%molecularLevel(subcell,:) ! Store previous level so can use updated one in future
                    thisOctal%newMolecularLevel(subcell,:) = thisOctal%molecularLevel(subcell,:) ! Store previous level so can use it for testing convergence at end

                    do while (.not. popsConverged)
                       iter = iter + 1
                       oldpops = thisOctal%newmolecularLevel(subcell,1:maxlevel) ! retain old pops before calculating new one

                       do iTrans = 1, maxtrans

                             call calculateJbar(grid, thisOctal, subcell, thisMolecule, ds(1:nRay), &
                                  phi(1:nRay), i0(iTrans,1:nRay), iTrans, thisOctal%jnu(subcell,iTrans), &
                                  thisOctal%newMolecularLevel(subcell,1:maxlevel)) ! calculate updated Jbar

                       enddo

                       call solveLevels(thisOctal%newMolecularLevel(subcell,1:maxlevel), &
                            thisOctal%jnu(subcell,1:maxtrans), dble(thisOctal%temperature(subcell)), &
                            thisMolecule, thisOctal%nh2(subcell))

                       fac = abs(maxval((thisOctal%newMolecularLevel(subcell,1:minlevel) - oldpops(1:minlevel)) &
                            / oldpops(1:minlevel))) ! convergence criterion ! 6 or 8?
                         
                       if (fac < 1.d-6) then
                          popsConverged = .true.
                       endif

                       if (iter == maxIter) then
                          popsConverged = .true.
                          call writeWarning("Maximum number of iterations reached in pop solver")
                       endif
                    enddo

!                    thisOctal%oldmolecularlevel(1:maxlevel,subcell) = thisOctal%molecularlevel(1:maxlevel,subcell) ! Store previous iteration levels for convergence
!                    thisOctal%molecularlevel(1:maxlevel,subcell) = thisOctal%newmolecularlevel(1:maxlevel,subcell) ! Store new levels as good current ones. Should improve convergence?

                 endif
              enddo
           end do

           if(writeoutput) then
              write(message,*) "Getray"
              call tune(6, message)  ! start a stopwatch
           endif
           
#ifdef MPI
           if (.not.blockHandout) exit blockloop
        end do blockLoop
     end if ! (my_rank /= 0)


      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        if(my_rank == 0) write(*,*) "Updating MPI grids"

      ! have to send out the 'octalsBelongRank' array
      call MPI_BCAST(octalsBelongRank,SIZE(octalsBelongRank),  &
                     MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

      call countVoxels(grid%octreeRoot,nOctal,nVoxels)
      allocate(tArrayd(1:nVoxels))
      allocate(tempArrayd(1:nVoxels))
      tArrayd = 0.d0
      tempArrayd = 0.d0
      do i = 1, maxlevel
        tArrayd = 0.d0
        call packMoleLevel(octalArray, nVoxels, tArrayd, octalsBelongRank, i)
        call MPI_ALLREDUCE(tArrayd,tempArrayd,nVoxels,MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,ierr)
        tArrayd = tempArrayd
        call unpackMoleLevel(octalArray, nVoxels, tArrayd, octalsBelongRank, i)
     enddo
      deallocate(tArrayd, tempArrayd)

      if(my_rank == 0) write(*,*) "Done updating"
#endif

      ! Calculate convergence towards solution (per level)
      
      maxRMSFracChange = -1.d30

!      call solveAllPops(grid, grid%octreeroot, thisMolecule)
      call writeinfo("Calculating convergence data", FORINFO)
      call calculateConvergenceData(grid, nvoxels, fixedrays, maxRMSFracChange)

 ! If you think it's converged then test with the same number of rays to make sure

        if(fixedrays) then 
           tol = tolerance * 0.5
        else
           tol = tolerance
        endif

        if(maxRMSFracChange < (tol ** 2) * real(nvoxels)) then
           if (gridConvergedTest) gridConverged = .true.
!           gridConverged = .true.
           gridConvergedTest = .true.
        else
           gridConvergedTest = .false.
           gridConverged = .false.
        endif

!        if (writeoutput) then
!           write(*,*) maxRMSFracChange, tol ** 2 * real(nvoxels)
!           write(*,*) gridconvergedtest, gridconverged, "a"
!        endif

        if(myrankiszero) then
           call writeAmrGrid("molecular_tmp.grid",.false.,grid)
           call writeinfo("written grid", TRIVIAL)
	endif

        call writeinfo("Dumping results", FORINFO)
        if((amr1d .or. amr2d) .and. writeoutput .and. (.not. gridconvergedtest)) call dumpresults(grid, thisMolecule)!, convtestarray) ! find radial pops on final grid     

        if(writeoutput) then
           write(message,*) "Done ",nray," rays"
           call tune(6, message)  ! stop a stopwatch

        endif
        
!        write(*,*) "fixedrays", fixedrays
!        write(*,*) "gridconverged", gridconverged
!        write(*,*) "gridconvergedtest", gridconvergedtest
!        write(*,*) "molebench", molebench

!	if(writeoutput) write(*,*) gridconvergedtest, gridconverged, fixedrays, molebench

        if (.not.gridConvergedTest) then
           if (.not.gridConverged) then               
              if (.not.fixedRays) nRay = nRay * 2 !double number of rays if convergence criterion not met and not using fixed rays - revise!!! Can get away with estimation?
              write(message,*) "Trying ",nRay," Rays"
              call writeInfo(message,FORINFO)
              if (molebench .and. (nray .gt. 400000)) then 
                 if(writeoutput) write(*,*) "Molebench Test Ended, Exiting..."
                 gridconverged = .true.
              endif
 
           endif
        else
           if(.not. gridconverged) call writeinfo("Doing all rays again", FORINFO)
        endif
        
        if (nRay > maxRay) then
           nRay = maxRay  ! stop when it's not practical to do more rays
           call writeWarning("Maximum number of rays exceeded - capping")
        endif

!        if(gridconvergedtest) gridconverged = .true.
!        write(*,*) "gridconverged", gridconverged

#ifdef MPI
        deallocate(octalsBelongRank)
#endif
        deallocate(ds, phi, i0)
     enddo

!     if(writeoutput) then
!        call writeinfo("Writing convergence profile", TRIVIAL)
!     endif

!     if(writeoutput .and. molebench .and. amr1d) then
!        open(141,file="convergenceprofile.dat",status="unknown",form="formatted")

!        do i=1,100
!           r(i) = log10(rinner) + dble(i-1)/dble(100-1)*((log10(router) - log10(rinner))) ! log(radius) -> radius
!        enddo
        
!        write(141,'(tr3,100(es8.3e1,tr1))') 10.**r

!        do i = 1,grand_iter
!           write(141,'(i2,tr1,100(f8.6,tr1))') i, convtestarray(i,:,1)
!        enddo

!        close(141)
!     endif
  enddo
  
  close(33)

  call writeinfo("mole loop done.", FORINFO)
  stop
end subroutine molecularLoop

   subroutine getRay(grid, fromOctal, fromSubcell, position, direction, ds, phi, i0, thisMolecule, fixedrays)

     use input_variables, only : useDust, amrgridsize

     type(MOLECULETYPE) :: thisMolecule
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal, startOctal, fromOctal

     integer, parameter :: maxSamplePoints = 100

     logical, optional :: fixedrays
     logical :: stage1
     integer :: fromSubcell
     integer :: subcell
!     real(double) :: ds, phi, i0(:), r, di0
     real(double) :: ds, phi, i0(:), r, di0(maxtrans)
     integer :: iTrans,nTrans
     type(VECTOR) :: position, direction, currentPosition, thisPosition, thisVel, rayvel
     type(VECTOR) :: startVel, endVel, endPosition

     real(double) :: alphanu(maxtrans,2), alpha(maxtrans), alphatemp(maxtrans), snu(maxtrans), jnu(maxtrans)
     integer :: iLower(maxtrans) , iUpper(maxtrans)
     real(double) :: dv, deltaV
     integer :: i
     real(double) :: dist, dds, tval
     integer :: nTau
     real(double) :: nLower(maxtrans)
!     real(double),pointer :: nLower(:)
     real(double) :: nUpper(maxtrans)
!     real(double),pointer :: nupper(:)
!     real(double) :: dTau, etaline(maxtrans), kappaAbs 
     real(double) :: dTau(maxtrans), kappaAbs, localradiationfield(maxtrans), attenuation(maxtrans)
!     real(double), allocatable :: tau(:)
     real(double) :: tau(maxtrans)
     real(double) :: dvAcrossCell, projVel, endprojVel
     real(double) :: CellEdgeInterval, phiProfVal
!     real :: lambda
     real(double), allocatable, save :: lambda(:), OneArray(:)
     integer :: ilambda
     integer :: iCount

     real(double), save :: r1(5) !! quasi random number generators
     real(double), save :: OneOverNTauArray(maxsamplePoints)

     logical,save :: firsttime = .true.
     logical :: quasi = .false.

     logical :: realdust
     real(double),save :: BnuBckGrnd(128)
     real(double) :: balance(maxtrans), spontaneous(maxtrans), snugas(maxtrans)
     real(double) :: nMol
     integer :: done(maxtrans)
 !    real(double) :: z,

     realdust = .true.

     if(present(fixedrays)) then
        stage1 = fixedrays
     else
        stage1 = .false.
     endif

     nTrans = thisMolecule%nTrans

!     allocate(tau(1:maxTrans))

     position = randomPositionInCell(fromOctal, fromsubcell)

     icount = 1
     thisOctal => fromOctal!grid%octreeRoot
     subcell = fromSubcell

     rayVel = Velocity(position, grid, thisOctal, subcell)

 !!! quasi random code!!!

     if(firsttime) then

        allocate(lambda(maxtrans))
        allocate(OneArray(maxtrans))
        OneArray = 1.0d0
       
        do itrans = 1, maxtrans
           BnuBckGrnd(itrans) = Bnu(thisMolecule%transfreq(itrans), Tcbr)
!          Bnudust(itrans) = Bnu(thisMolecule%transfreq(itrans), thisOctal%dusttemperature(subcell))
        enddo

       do ntau = 2,maxsamplepoints
          OneOverNtauArray(ntau) = 1.d0 / (dble(nTau) - 1.d0)
       enddo

       lambda(:) = (cspeed / thisMolecule%transfreq(1:maxtrans)) * 1d8 ! cm to Angstroms

       firsttime =.false.
       if (quasi) call sobseq(r1,-1)
    endif

     if(quasi) then     
        call sobseq(r1)    
        direction = specificUnitVector(r1(1),r1(2))
 !       r = (r1(3)+r1(4)+r1(5)) / 3.
        deltaV = 4.3 * thisOctal%microturb(subcell) * (r1(3) - 0.5d0) ! random frequency near line centre
     else
        call random_number(r)

        direction = randomUnitVector() ! Put this line back if you want to go back to pseudorandom

!       randompointingrid = randompointinRegion(amrgridsize * 0.5)
!       randompointingrid = 1d7 * randomUnitVECTOR()
!       call normalize(randompointingrid)
!       unitposition = position
!       call normalize(unitposition)
!       write(*,*) "unitposition",unitposition
!       write(*,*) "randompointingrid",randompointingrid
!       direction = UnitVECTORBetweenTwoPoints(randompointingrid,position)
!       write(*,*) "direction", direction

        deltaV = 4.3 * thisOctal%microturb(subcell) * (r - 0.5d0) ! random frequency near line spectrum peak. 

!        r = GaussianComplementRand()
!        deltaV = thisOctal%microturb(subcell) * r ! random frequency near line spectrum peak. 

!        write(55,*) r

     endif

!!! SCIENCE LINES

!     deltaV = deltaV + (rayVel .dot. direction) ! transform to take account of grid velocity
!     projVel = deltaV - (rayVel .dot. direction) ! transform back to velocity relative to local flow

!!! COMPUTATIONALLY QUICKER TO DO THIS

     projvel = deltaV
     deltaV =  deltaV + (rayVel .dot. direction)

     call distanceToCellBoundary(grid, position, direction, ds, sOctal=thisOctal)

     currentPosition = position + ds * direction

     if (inOctal(grid%octreeRoot, currentPosition)) then ! check that we're still on the grid
        thisVel = velocity(currentposition, grid) 
        endprojVel = deltaV - (thisVel.dot.direction)
        phi = phiProf((endprojVel+projvel)*0.5d0, thisOctal%molmicroturb(subcell))
     else
        phi = phiProf(projVel, thisOctal%molmicroturb(subcell)) ! if fell off grid then assume phi is unchanged
     endif

!     write(*,*) "phi", phi
!     currentPosition = position !This line is wrong because the source function does not apply locally
     ds = ds * 1.d10 ! convert from cm to torus units

     i0 = 0.d0
     tau = 0.d0

     iUpper(:) = thisMolecule%iTransUpper(:)
     iLower(:) = thisMolecule%iTransLower(:)

     do while(inOctal(grid%octreeRoot, currentPosition))

        call findSubcellLocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

        if(useDust) then             
           do itrans = 1,maxtrans
              if(realdust) then
                 call locate(grid%lamArray, size(grid%lamArray), real(lambda(itrans)), ilambda)
                 call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = real(lambda(itrans)), kappaAbs = kappaAbs)
              else
!             kappaAbs = thisOctal%rho(subcell) * 0.1 * thisMolecule%transfreq(itrans) * 1e-12 * 1d10 !multiplied by density !cm -> torus
                 kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10 !multiplied by density !cm -> torus
              endif
              alphanu(itrans,2) = kappaAbs * 1.d-10 !torus units -> cm !1d-10 is right    !kappa already multiplied by density in amr_mod
           enddo
        endif

        nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)

!        nlower(:) = 
!        nUpper(:) = thisOctal%molecularLevel(iUpper(:),subcell) * nMol
        balance(:) = (hcgsOverFourPi * nmol) * (thisOctal%molecularLevel(subcell,ilower(:)) * thisMolecule%einsteinBlu(:) - &
                     thisOctal%molecularLevel(subcell,iupper(:)) * thisMolecule%einsteinBul(:))

        spontaneous(:) = (hCgsOverfourPi * nmol) * thisMolecule%einsteinA(:) * thisOctal%molecularLevel(subcell,iupper(:))

        snu(:) = spontaneous(:) / balance(:) ! Source function  -only true if no dust else replaced by gas test


!        alphaTemp(:) = hCgsOverFourPi * balance(:) ! Equation 8
          
        startVel = velocity(currentposition, grid) 
        endPosition = currentPosition + tval * direction
        endVel = velocity(endposition, grid)

        dvAcrossCell = ((startVel - endvel) .dot. direction) ! start - end
        dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))

        nTau = min(max(2, nint(dvAcrossCell * 5.d0)), maxSamplePoints) ! selects dVacrossCell as being between 0.1 and 10 (else nTau bounded by 2 and 200)

        CellEdgeInterval = OneOverNtauArray(ntau)
        dds = tval * cellEdgeInterval

        dist = 0.d0
        done = 0

        do i = 2, nTau

           dist = dist + dds ! increment along distance counter
              
           thisPosition = currentPosition + dist * direction

           thisVel = velocity(thisPosition, grid)
           dv = deltaV - (thisVel .dot. direction)
           PhiProfVal = phiProf(dv, thisOctal%molmicroturb(subcell))
           
           if(usedust) then
              alphanu(1:maxtrans,1) = phiprofval * balance(1:maxtrans) ! Equation 8
              alpha(1:maxtrans) = alphanu(1:maxtrans,1) + alphanu(1:maxtrans,2)
              jnu(1:maxtrans) = spontaneous(1:maxtrans) * phiProfVal + &
                   alphanu(1:maxtrans,2) * thisOctal%bnu(subcell,1:maxtrans) ! jnu, emission coefficient - equation 7
              
              do itrans = 1,maxtrans
                 if(alpha(itrans) .ne. 0.d0) then 
                    snu(itrans) = jnu(itrans) / alpha(itrans)
                 else
                    snu(itrans) = tiny(snu)
                    alpha(itrans) = tiny(alpha)
                 endif
              enddo
           else
              alpha(:) = balance(:) * phiprofVal 
           endif

           dTau(:) = alpha(:) * dds * 1.d10 ! dds is interval width & optical depth, dTau = alphanu*dds - between eqs (3) and (4)
           attenuation(:) = exp(-tau(:))
!           call vdexp(maxtrans, -tau, attenuation)
           localradiationfield(:) = exp(-dtau(:)) 
           localradiationfield(:) = OneArray(:) - localradiationfield(:)
           di0(:) = attenuation(:) * localradiationfield(:) * snu(:) ! 2nd term is local radiation field from this cell 
!     write(*,*) "test",attenuation(1),localradiationfield(1), snu(1), balance(1), spontaneous(1), nlower(1), nupper(1)
!     write(*,*) "test2",currentposition,thisoctal%temperature(subcell),thisoctal%nh2(subcell),thisoctal%molabundance(subcell)
           do itrans = 1, maxtrans

              if(di0(itrans) .gt. 1d-6 * i0(itrans)) then
!                 if(tau(itrans) .lt. 10.d0) then
                 i0(itrans) = i0(itrans) + di0(itrans) ! summed radiation intensity from line integral 
                 tau(itrans) = tau(itrans) + dtau(itrans) ! contribution to optical depth from this line integral
              else
                 done(itrans) = 1
                 i0(itrans) = i0(itrans) + di0(itrans) ! summed radiation intensity from line integral 
                 tau(itrans) = tau(itrans) + dtau(itrans) ! contribution to optical depth from this line integral
!                 if(sum(done(1:mintrans)) .eq. mintrans) then
                 if(sum(done) .ge. mintrans) then
                    goto 118 
                 endif
              endif
              
           enddo
           
        enddo
        
        currentPosition = currentPosition + (tval + 1.d-3*grid%halfSmallestSubcell) * direction ! FUDGE - make sure that new position is in new cell
     enddo
     
118  continue
     
     i0(:) = i0(:) + BnuBckGrnd(:) * attenuation(:)

   end subroutine getRay

   subroutine calculateJbar(grid, thisOctal, subcell, thisMolecule, ds, phi, i0, iTrans, jbar, nPops)

     use input_variables, only : useDust
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(MOLECULETYPE) :: thisMolecule

     real(double) :: ds(:), phi(:), i0(:), nPops(:)
     real(double) :: phids(nray), tauArray(nray), opticaldepthArray(nray), jbarinternalArray(nray), jbarExternalArray(nray)
     integer :: iTrans
     real(double) :: jbar
     integer :: iRay
     real(double) :: nLower, nUpper, nMol
     real(double) :: jBarInternal, jBarExternal
     real(double) :: alphanu(2), jnu, jnuDust, etaline, etalineBase, alphanuBase, kappaAbs, alpha
     integer :: iUpper, iLower
     real(double) :: tau, opticaldepth, snu, sumPhi
     real :: lambda
     integer :: ilambda, i
     logical :: realdust = .true.
     
     jBarExternal = 1.d-60
     jBarInternal = 1.d-60

     jnu = 1.d-60

     iUpper = thisMolecule%iTransUpper(iTrans)
     iLower = thisMolecule%iTransLower(iTrans)
     ! commented elsewhere in the code
     sumPhi = 0.d0

     nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
     etaLineBase = thisMolecule%einsteinA(iTrans)! * thisMolecule%transFreq(iTrans)

     nLower = nPops(iLower) * nMol
     nUpper = nPops(iUpper) * nMol

     etaLine = etaLineBase * nUpper

     alphanuBase =(nLower * thisMolecule%einsteinBlu(iTrans) - &
                  nUpper * thisMolecule%einsteinBul(iTrans))

     if(useDust) then
        lambda = cspeed / thisMolecule%transfreq(iTrans) * 1d8
        call locate(grid%lamArray, size(grid%lamArray), lambda, ilambda)

        if(realdust) then
           call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = lambda, kappaAbs = kappaAbs)
        else
           kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10 !multiplied by density !cm -> torus
        endif

        alphanu(2) = kappaAbs * 1.d-10 !* thisOctal%rho(subcell) - already done
        jnuDust = alphanu(2) * bnu(thisMolecule%transfreq(iTrans), dble(thisOctal%temperaturedust(subcell)))

        do iRay = 1, nRay
           
           alphanu(1) = alphanuBase * phi(iray) * hCgsOverFourPi!/thisMolecule%transFreq(iTrans)
           alpha = alphanu(1) + alphanu(2)

           jnu = (etaLine) * phi(iRay) * hCgsOverFourPi!/thisMolecule%transFreq(iTrans)
           jnu = jnu + jnuDust
           
           if (alpha /= 0.d0) then
              snu = jnu/alpha
           else
              snu = tiny(snu)
           endif

           tau = alpha * ds(iray)
           opticaldepth = exp(-tau)

           jBarExternal = jBarExternal + i0(iray) * opticaldepth * phi(iRay)
           jBarInternal = jBarInternal + snu * (1.d0 - opticaldepth) * phi(iRay)

           sumPhi = sumPhi + phi(iRay)
        enddo
     else

        alpha = alphanuBase
        jnu = etaLine

        if (alpha /= 0.d0) then
           snu = jnu/alpha
        else
           snu = tiny(snu)
        endif

        phids(1:nRay) = phi(1:nray) * ds(1:nray)
        tauArray(1:nRay) = alpha * hCgsOverFourPi * phids(1:nRay)
        opticaldepthArray(1:nRay) = exp(-1.d0 * tauArray(1:nRay))

        jBarExternalArray(1:nRay) = i0(1:nRay) * opticaldepthArray(1:nRay) * phi(1:nRay)
        jBarInternalArray(1:nRay) = snu * (1.d0 - opticaldepthArray(1:nRay)) * phi(1:nRay)

!        write(*,*) "tau",opticaldepthArray
!        write(*,*) "phi",phi
!        write(*,*) "INT",jbarinternalarray
!        write(*,*) "ext",jbarinternalarray

        sumPhi = SUM(phi(:))
        jbar = (sum(jBarExternalArray(1:nRay)) + sum(jBarInternalArray(1:nRay)))/sumPhi

     endif

     if(usedust)then
        sumPhi = SUM(phi(:))
        jbar = (jBarExternal + jBarInternal)/sumPhi
     endif

   
   end subroutine calculateJbar

 ! solves rate equation in matrix format - equation 10
   subroutine solveLevels(nPops, jnu,  temperature, thisMolecule, nh2)
     real(double) :: nPops(:)
     real(double) :: temperature
     real(double) :: jnu(:)
     real(double) :: nh2
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: matrixA(maxlevel+1,maxlevel+1), matrixB(maxlevel+1), collMatrix(maxlevel+1,maxlevel+1), cTot(maxlevel)
     real(double) :: boltzFac
     integer :: nLevels
     integer :: i, j
     integer :: itrans, l, k, iPart
     real(double) :: collEx, colldeEx

     character(len=80) :: message

     matrixA = 1.d-60 ! Initialise rates to negligible to avoid divisions by zero ! used to be 1d-10
     matrixB = 1.d-60 ! Solution vector - all components (except last) => equilibrium ! used to be 1d-10

     matrixB(maxLevel+1) = 1.d0 ! Sum over all levels = 1 - Conservation constraint

 ! This do loop calculates the contribution to transition rates of each level from every other level. NB emission is +ve here
     do iTrans = 1, maxtrans

        k = thisMolecule%iTransUpper(iTrans)
        l = thisMolecule%iTransLower(iTrans)

        matrixA(k,k) = matrixA(k,k) + thisMolecule%einsteinBul(iTrans) * jnu(iTrans) + thisMolecule%einsteinA(iTrans)
        matrixA(l,l) = matrixA(l,l) + thisMolecule%einsteinBlu(iTrans) * jnu(iTrans)
        matrixA(k,l) = matrixA(k,l) - thisMolecule%einsteinBlu(iTrans) * jnu(iTrans)
        matrixA(l,k) = matrixA(l,k) - thisMolecule%einsteinBul(iTrans) * jnu(iTrans) - thisMolecule%einsteinA(iTrans)

     enddo

     collMatrix = 1.d-60 ! used to be 1.d-10

 ! Calculate contribution from collisions - loop over all collision partners and all molecular levels

     do iPart = 1, thisMolecule%nCollPart
        do iTrans = 1, maxtrans !thisMolecule%nCollTrans(iPart)

           k = thisMolecule%iCollUpper(iPart, iTrans)
           l = thisMolecule%iCollLower(iPart, iTrans)

           boltzFac = exp(-abs(thisMolecule%energy(k)-thisMolecule%energy(l)) / (kev*temperature))
           colldeEx = collRate(thisMolecule, temperature, iPart, iTrans) * nh2
           collEx = colldeEx * boltzFac * thisMolecule%g(k) / thisMolecule%g(l)

           collMatrix(l, k) = collMatrix(l, k) + collEx
           collMatrix(k, l) = collMatrix(k, l) + colldeEx

        enddo
     enddo

     cTot = 1d-60

     do k = 1, maxlevel
        do l = 1, maxlevel
           cTot(k) = cTot(k) + collMatrix(k,l) ! sum over all collisional rates out of each level
        enddo
     enddo

     do i = 1, maxlevel
           matrixA(i,i) = matrixA(i,i) + cTot(i)
        do j = 1, maxlevel
           if (i .ne. j) then
              matrixA(i,j) = matrixA(i,j) - collMatrix(j, i)
           endif
        enddo
     enddo

     matrixA(maxlevel+1,1:maxlevel+1) = 1.d0 ! sum of all level populations
     matrixA(1:maxlevel+1,maxlevel+1) = 1.d-60 ! fix highest population to small non-zero value
     
!     do i = 1, thisMolecule%nlevels
!        do j = 1, thisMolecule%nlevels
!           if( i .eq. j ) then
!              matrixA(i,j) = matrixA(i,j)
!           else
!              matrixA(i,j) = 1.d-10
!           endif
!        enddo
!     enddo

 ! finished creating equation 10, now solve it to find new level populations
   call luSlv(matrixA, matrixB, maxlevel+1)

   matrixB = abs(matrixB) ! stops negative level populations causing problems

   if(.not.isnan(matrixB(1))) then
      nPops(1:maxLevel) = matrixB(1:maxLevel)
   else
      call LTEpops(thisMolecule, temperature, npops(1:maxlevel))

      write(message, *) "bad",temperature,jnu(1)
      call writeinfo(message, IMPORTANT)
!      write(*,*) npops
!      stop
   endif
   
 end subroutine solveLevels


 ! Calculate collision rates between partners for given temperature
   real(double) function collRate(thisMolecule, temperature, iPart, iTrans)
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: temperature, r
     integer :: iTrans, k, iPart

     collRate = 0.d0

!     call locate(thisMolecule%collTemps(iPart,1:thisMolecule%nCollTemps(iPart)), temperature, k)
     call locate(thisMolecule%collTemps(iPart,1:thisMolecule%nCollTemps(iPart)), &
          thisMolecule%nCollTemps(iPart), temperature, k)


     r = (temperature - thisMolecule%collTemps(iPart,k)) / &
          (thisMolecule%collTemps(iPart,k+1) - thisMolecule%collTemps(iPart, k))

     collRate = collRate + thisMolecule%collRates(iPart, iTrans, k) + &
          r * ( thisMolecule%collRates(iPart, iTrans, k+1) -  thisMolecule%collRates(iPart, iTrans, k))

   end function collRate

 ! this subroutine calculates the maximum fractional change in the first 6 energy levels
   recursive subroutine  swapPops(thisOctal, maxFracChangePerLevel, avgFracChange, counter, &
                          iter, nVoxels,fixedrays)

     use input_variables, only : tolerance

     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i
     real(double) :: maxFracChangePerLevel(:), maxFracChange, avgFracChange(:,:)
     real(double) :: newFracChangePerLevel(minlevel), temp(minlevel)
     integer :: counter(:,:),j

     integer :: iter
     integer :: nVoxels

     logical :: fixedrays
     integer :: ntrans

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call swapPops(child, maxFracChangePerLevel, avgFracChange, counter, iter, &
                      nVoxels, fixedrays)
                 exit
              end if
           end do
        else
           
!           thisOctal%molecularLevel(:,subcell) = thisOctal%oldmolecularLevel(:,subcell) 

           maxFracChange = MAXVAL(maxFracChangePerLevel(1:minlevel))

           newFracChangePerLevel = abs(((thisOctal%newMolecularLevel(subcell,1:minlevel) - &
                thisOctal%molecularLevel(subcell,1:minlevel)) / &
                thisOctal%newmolecularLevel(subcell,1:minlevel)))

           temp = newFracChangePerLevel(1:minlevel)

           do j=1, mintrans
                 if(newFracChangePerLevel(j) < 0.5 * tolerance) counter(4,j) = counter(4,j) + 1 ! used for itransdone
                 if(newFracChangePerLevel(j) < tolerance) counter(1,j) = counter(1,j) + 1
                 if(newFracChangePerLevel(j) < 2. * tolerance) counter(2,j) = counter(2,j) + 1
                 if(newFracChangePerLevel(j) < 5. * tolerance) counter(3,j) = counter(3,j) + 1
           enddo

           if(maxval(temp) .lt. 5. * tolerance) then
              counter(3,mintrans + 1) = counter(3,mintrans + 1) + 1
              if(maxval(temp) .lt. 2. * tolerance) then
                 counter(2,mintrans + 1) = counter(2,mintrans + 1) + 1
                 if(maxval(temp) .lt. tolerance) then
                    counter(1,mintrans + 1) = counter(1,mintrans + 1) + 1
                 endif
              endif
           endif

              avgFracChange(:,1) = avgFracChange(:,1) + temp
              avgFracChange(:,2) = avgFracChange(:,2) + temp**2

           If (maxval(temp) > maxFracChange) then
              maxFracChangePerLevel = newFracChangePerLevel ! update maxFracChange if fractional change is great => not converged yet
           endif

          thisOctal%molecularLevel(subcell,:) = thisOctal%newmolecularLevel(subcell,:)

!           diff = thisOctal%molecularLevel(:,subcell) - &
!                  thisOctal%newmolecularLevel(:,subcell)

 !          if(fixedrays .and. modulus(SubcellCentre(thisOctal,subcell)) .lt. 1e7 .and. modulus(SubcellCentre(thisOctal,subcell)) .gt. 3e6) then
!           if(.not. fixedrays) then
!              if (maxval(temp) .gt. tolerance * 6.) then
!                 thisOctal%molecularLevel(:,subcell) = &
!                      thisOctal%newmolecularLevel(:,subcell)! - (0.0)*diff! (enhanced) update molecular levels
!              else
!                 thisOctal%molecularLevel(:,subcell) = &
!                      thisOctal%newmolecularLevel(:,subcell)! - (0.0)*diff! (damped) update molecular levels - reduced from 0.2 to quicken convrgnce
!              endif

!              if(maxval(temp) .lt. tolerance * 2.)then
!                 thisOctal%molecularLevel(:,subcell) = &
!                      thisOctal%newmolecularLevel(:,subcell)! - 0.0 * diff! (damped) update molecular levels - reduced from 0.2 to quicken convrgnce
!              endif

!           else

!              if (maxval(temp) .gt. tolerance * 10.) then
!                 thisOctal%molecularLevel(:,subcell) = &
!                      thisOctal%newmolecularLevel(:,subcell)! - (0.0)*diff! (enhanced) update molecular levels
!              else
!                 thisOctal%molecularLevel(:,subcell) = &
!                      thisOctal%newmolecularLevel(:,subcell)! - (0.0)*diff! (damped) update molecular levels - reduced from 0.2 to quicken convrgnce
!              endif

!              if(maxval(temp) .lt. tolerance * 2.)then
!                 thisOctal%molecularLevel(:,subcell) = &
!                      thisOctal%newmolecularLevel(:,subcell)! - 0.0 * diff! (damped) update molecular levels - reduced from 0.2 to quicken convrgnce
!              endif
!           endif

        endif
     enddo

   end subroutine swapPops

!!!!READ ROUTINES!!!!

 subroutine calculateMoleculeSpectrum(grid, thisMolecule)
   use input_variables, only : itrans, nSubpixels, inc, usedust

#ifdef MPI
       include 'mpif.h'
#endif

   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   type(VECTOR) :: unitvec, posvec, centrevec, viewvec
   type(DATACUBE) ::  cube

   character (len=80) :: filename

   real(double) :: mean(6) = 0.d0
   integer :: icount = 0


#ifdef MPI
     ! For MPI implementations
     integer       ::   my_rank        ! my processor rank
     integer       ::   np             ! The number of processes
     integer       ::   ierr           ! error flag

     ! FOR MPI IMPLEMENTATION=======================================================
     !  Get my process rank # 
     call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

     ! Find the total # of precessor being used in this run
     call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#endif

!     if(grid%geometry .eq. "iras04158") call plotdiscvalues(grid, thisMolecule)

!     call readAMRgrid("molecular_tmp.grid",.false.,grid)

    !if(usedust .and. grid%geometry .eq. "iras04158") call testOpticalDepth(grid, thisMolecule)
!call testOpticalDepth(grid, thisMolecule)
    
     if(writeoutput) call writeinfo('Setting observer parameters', TRIVIAL)
     call setObserverVectors(inc, viewvec, unitvec, posvec, centreVec, grid%octreeroot%subcellsize)
     if(writeoutput) call writeinfo('Creating Image', TRIVIAL)
!     write(*,*)"posvec",posvec,"inc", inc
     Call createimage(cube, grid, unitvec, posvec, thismolecule, itrans, nSubpixels) ! create an image using the supplied parameters (also calculate spectra)
     
     call cubeIntensityToFlux(cube, thismolecule, itrans) ! convert intensity (ergcm-2sr-1Hz-1) to flux (Wm-2Hz-1) (so that per pixel flux is correct)

     !Commented out lines are for removing background intensity - will want to do this one day

     !   call TranslateCubeIntensity(cube,1.d0*Tcbr) ! 
     !   call cubeIntensityToFlux(cube, thismolecule, itrans) ! convert intensity (ergcm-2sr-1Hz-1) to flux (Wm-2Hz-1) (so that per pixel flux is correct)
     
     call createFluxSpectra(cube, thismolecule, itrans)

     if(writeoutput) call writedatacube(cube, 'MolRT.fits')

     write(filename,'(a,a,i1,a,i1,a)') trim(cube%telescope%label),'fluxcubeJ', &
          thisMolecule%itransUpper(itrans)-1,'-',thisMolecule%itransLower(itrans)-1,'.ps/vcps'
     stop

    end subroutine calculateMoleculeSpectrum

   function makeImageGrid(cube, unitvec, posvec, grid, thisMolecule, iTrans, deltaV, nsubpixels) result (imagegrid)

     use input_variables, only : npixels, imageside

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     integer :: itrans
     type(VECTOR) :: unitvec, viewVec, posvec
     type(VECTOR) :: imagebasis(2), pixelcorner
     integer :: nsubpixels, subpixels
     real :: imagegrid(npixels,npixels)
     integer :: ipixels, jpixels
     real(double) :: pixelside
     real(double) :: deltaV
     type(datacube) :: cube
     integer :: index(2)

     logical, save :: firsttime = .true.

     if(firsttime) then

        imagebasis(2) = unitVec .cross. VECTOR(1d-20,1d-20,1d0) ! gridvector
        imagebasis(1) = imagebasis(2) .cross. unitVec ! gridvector perp

        !    write(*,*) unitvec, "HELP"

        call normalize(imagebasis(1))
        call normalize(imagebasis(2))
        !   call normalize(posvec)

        pixelside = imageside / dble(npixels)
        imagebasis(1) = imagebasis(1) * pixelside ! rescale basis vectors so that stepsize is simplified 
        imagebasis(2) = imagebasis(2) * pixelside
        viewvec = (-1.d0) * unitvec ! look *towards* origin from posvec
        
!        firsttime = .false.
        endif

        write(*,*) "viewvec", viewvec

        pixelcorner = posvec - (dble(npixels)/2.d0)*(imagebasis(2) - imagebasis(1)) + imagebasis(1) ! fudge at the end to 'fit in' with the do loop

        if (nsubpixels .gt. 0) then 
           subpixels = nsubpixels
        else
           subpixels = 0
        endif

        do jpixels = 1, npixels ! raster over image
           if (jpixels .eq. 1) then 
              pixelcorner = pixelcorner - imagebasis(1)
           else
              pixelcorner = pixelcorner - (imagebasis(1) + dble(npixels)*imagebasis(2))
           endif

           do ipixels = 1, npixels
              index = (/ipixels,jpixels/)
              imagegrid(ipixels,jpixels) = newPixelIntensity(cube,pixelside,viewvec,pixelcorner,imagebasis,grid,thisMolecule,&
                   iTrans,deltaV, subpixels,index)
              pixelcorner = pixelcorner + imagebasis(2)
           enddo
        enddo

 end function makeImageGrid

 !!! Calculates the intensity for a square pixel of arbitrary size, position, orientation

 function PixelIntensity   (cube,pixelside,viewvec,pixelcorner,&
                           imagebasis,grid,thisMolecule,iTrans,deltaV,subpixels,index)&
                           result(totalpixelintensity)

   use input_variables, only : tolerance
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     integer :: itrans
     type(VECTOR) :: viewVec
     real(double) :: i0, opticaldepth
     real(double) :: totalPixelIntensity, oldTotalPixelIntensity
     type(VECTOR) :: imagebasis(2), pixelbasis(2), pixelcorner, newposvec

     integer :: nsubpixels, subpixels
     real(double) :: OneOverSubPixelsSquared
     integer :: i,j
     integer :: index(2)
     integer, parameter :: maxSubPixels = 32
 !    real(double),allocatable :: subpixelgrid(:,:)
     real(double) :: pixelside,subpixelsize !,subPixelSolidAngle

     logical :: converged, romberg
     real(double) :: deltaV
     type(DATACUBE) :: cube

     converged = .false. ! failed flag

     nsubpixels = subpixels ! dummy variable
     cube%converged(index(1),index(2),cube%nv) = 1
     oldtotalPixelIntensity = 1.d-30

     if(subpixels .eq. 0) then 
        romberg= .true.
        subpixels = 1
     endif

     do while((.not. converged))

 !    allocate(subpixelgrid(subpixels,subpixels)) ! square array of pixels that make up the larger pixel 

     cube%nsubpixels(index(1),index(2),cube%nv) = subpixels
     OneOverSubpixelsSquared = 1.d0/dble(subpixels)**2

     pixelbasis(1) = imagebasis(1) / dble(subpixels)
     pixelbasis(2) = imagebasis(2) / dble(subpixels)
     newposvec = pixelcorner + (pixelbasis(1) + pixelbasis(2)) / 2.d0 ! Code takes position passed from makeImageGrid and puts it on the true pixel corner
     subpixelsize = (pixelside/dble(subpixels))**2 ! Area the pixel covers
 !    subpixelSolidAngle = subpixelSize/(fourpi*cube%obsdistance**2) ! Fraction of solid area at distance

 !    subpixelgrid = 0.d0 
     totalpixelintensity = 0.d0
 !    romberg = .true. ! always on at the moment
 !    if(subpixels .ne. nsubpixels) romberg = .true.

     do j = 1,subpixels ! this whole loop rasters across the imagegrid calculating the intensities at each position
        if (j .eq. 1) then 
           newposvec = newposvec - pixelbasis(1)
        else
           newposvec = newposvec - (pixelbasis(1) + dble(subpixels)*pixelbasis(2))
        endif

        do i = 1,subpixels
           call intensityalongray(newposvec,viewvec,grid,thisMolecule,itrans,deltaV,i0,opticaldepth)
 !          subpixelgrid(i,j) = i0
           totalPixelIntensity = totalPixelIntensity + i0!+ subpixelgrid(i,j) 
           newposvec = newposvec + pixelbasis(2)
        enddo
     enddo

 !    subpixelgrid = subpixelgrid / dble(subpixels)**2
     totalPixelIntensity = totalPixelIntensity * OneOverSubpixelsSquared

     if(romberg) then
        if(abs((totalPixelIntensity - oldtotalPixelIntensity)/TotalPixelIntensity) > tolerance) then
           oldTotalPixelIntensity = TotalPixelIntensity
           subpixels = subpixels * 2

        else

           converged = .true.
 !          TotalPixelIntensity = (dble(2**(log(dble(subpixels))/log(2.d0)))*TotalPixelIntensity - OldTotalPixelIntensity) &
 !                               /(dble(2**(log(dble(subpixels))/log(2.d0)))-1.d0) ! Richardson Extrapolation

           TotalPixelIntensity = (dble(subpixels)**2 * TotalPixelIntensity - OldTotalPixelIntensity) / (dble(subpixels)**2 -1.d0) 
           cube%nsubpixels(index(1),index(2),cube%nv) = subpixels
           subpixels = 0
        endif
     else
        converged = .true.
     endif

 !    deallocate(subpixelgrid)
  enddo

 end function PixelIntensity

 function newPixelIntensity(cube,pixelside,viewvec,pixelcorner,&
      imagebasis,grid,thisMolecule,iTrans,deltaV,subpixels,index)&
      result(pixelintensity)
   
   use input_variables, only : tolerance, nsubpixels
   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   integer :: itrans
   type(VECTOR) :: viewVec
   real(double) :: i0, opticaldepth

   type(VECTOR) :: imagebasis(2), pixelbasis(2), pixelcorner, rayposition
   
   integer :: subpixels, minrays
   integer :: i, iray
   integer :: index(2)

   real(double) :: pixelside
   real(double) :: avgIntensityNew, avgIntensityOld
   real(double) :: varIntensityNew, varIntensityOld
   real(double) :: PixelIntensity
   real(double) :: rtemp(2)
   real(double), save ::  r(10000,2)

   logical :: converged
   real(double) :: deltaV
   type(DATACUBE) :: cube
   logical, save :: firsttime = .true.
   
   if(firsttime) then
      
      call sobseq(rtemp,-1)
      
      do i = 1, 10000
         call sobseq(rtemp)
         r(i,:) = dble(rtemp)
      enddo
      
      firsttime = .false.
   endif

   pixelbasis(1) = imagebasis(1)
   pixelbasis(2) = imagebasis(2)
  
   avgIntensityOld = 0.
   varIntensityOld = 0.
   
   converged = .false. ! failed flag
     
   if(subpixels .ne. 0) then 
      converged = .true.
      minrays = nsubpixels
      cube%converged(index(1),index(2),cube%nv) = 1
   else
      minrays = -1
   endif
   
   iray = 1
   
   do while((.not. converged) .or. iray .le. minrays)  

      rayposition = pixelcorner + r(iray,2) * pixelbasis(2) - r(iray,1) * pixelbasis(1) ! random position in pixel
      
!      if(index(1) .eq. 12 .and. index(2) .eq. 14) then
!         call intensityalongray(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0,opticaldepth, .true.)
!         write(*,*) "ray",i0, opticaldepth
!      else

     call intensityalongray(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0,opticaldepth)
!      call lteintensityalongray(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0,opticaldepth)

!      endif

      avgIntensityNew = ((iray - 1) * avgIntensityOld + i0) / dble(iray)
      varIntensityNew = ((iray - 1) * varIntensityOld + ((i0 - avgIntensityNew) * (i0 - avgIntensityOld))) / dble(iray)
      
      if(varIntensityNew .lt. iray * (tolerance* avgIntensityNew)**2 .and. iray .gt. 1) then
         converged = .true.
         PixelIntensity = avgIntensityNew
!         write(*,*) "nrays = ", iray, "for",index
      !   write(*,*) "var   = ",varIntensityNew,varIntensityOld
      !   write(*,*) "avg   = ",avgIntensityNew,avgIntensityOld
         cube%nsubpixels(index(1),index(2),cube%nv) = iray
         cube%converged(index(1),index(2),cube%nv) = 1
         iray = iray + 1
      elseif(iray .gt. 10000) then
         PixelIntensity = avgIntensityNew
         converged = .true.
         cube%nsubpixels(index(1),index(2),cube%nv) = iray
         cube%converged(index(1),index(2),cube%nv) = 0
      else
!      write(*,'(a,es15.9,tr3,es15.9,tr3,es16.9)') "badvar   = ",varIntensityNew, iray * (tolerance*0.001 * avgIntensityNew)**2,varIntensityNew - iray * (tolerance*0.001 * avgIntensityNew)**2
         avgIntensityOld = avgIntensityNew
         varIntensityOld = varIntensityNew
         PixelIntensity = avgIntensityNew
         cube%nsubpixels(index(1),index(2),cube%nv) = iray
         iray = iray + 1
      endif

   enddo
   
 end function NewPixelIntensity

 !!! This subroutine takes the parameters supplied to it and makes an image by calling more subroutines 

   subroutine createimage(cube, grid, unitvec, posvec, thisMolecule, iTrans, nSubpixels)

     use input_variables, only : gridDistance, beamsize, npixels, nv, imageside, maxVel
#ifdef MPI
     include 'mpif.h'
#endif
     type(TELESCOPE) :: mytelescope
     type(MOLECULETYPE) :: thisMolecule
     type(GRIDTYPE) :: grid
     type(DATACUBE) :: cube
     type(VECTOR) :: unitvec, posvec
     real(double) :: minVel
     real(double) :: deltaV
     integer :: iTrans
     integer :: i !, j, k
 !    real(double) :: xval, yval, r
 !    integer :: nMonte, imonte
     integer :: iv
     integer :: nsubpixels
     character(len=200) :: message
     real(double) :: intensitysum, fluxsum, ddv!, dummy
     real(double), save :: background
     real(double), allocatable :: weightedfluxmap(:,:)
     real(double) :: weightedfluxsum, weightedflux
     real(double), allocatable :: fineweightedfluxmap(:,:)
     real(double) :: fineweightedfluxsum, fineweightedflux
!     real(double), pointer :: weight(npixels,npixels) => null()
     real(double), allocatable :: weight(:,:)

#ifdef MPI
     ! For MPI implementations
     integer       ::   my_rank        ! my processor rank
     integer       ::   np             ! The number of processes
     integer       ::   ierr           ! error flag
     integer :: ix1, ix2, n
     real(double), allocatable :: tempArray(:), tempArray2(:)

     ! FOR MPI IMPLEMENTATION=======================================================
     !  Get my process rank # 
     call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

     ! Find the total # of precessor being used in this run
     call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#endif

 ! SHOULD HAVE CASE(TELESCOPE) STATEMENT HERE

     mytelescope%label = 'JCMT'
     mytelescope%diameter = 15.d2 ! diameter in cm
     mytelescope%beamsize = beamsize

     minVel = (-1.d0) * maxVel

     call writeinfo("initcube",TRIVIAL)

     if(nv .eq. 0) then
        call initCube(cube, npixels, npixels, 200, mytelescope) ! Make cube
     else
        call initCube(cube, npixels, npixels, nv, mytelescope) ! Make cube
     endif

     cube%obsDistance = gridDistance * 1d10!(in cm) Additional information that will be useful
     write(message,*) "Observer Distance: ",gridDistance/pctocm, " pc"
     call writeinfo(message, FORINFO) 
     call addSpatialAxes(cube, -imageside/2.d0, imageside/2.d0, -imageside/2.d0, imageside/2.d0)

     if(nv .ne. 0) then 
        call addvelocityAxis(cube, minVel, maxVel) ! velocities in km/s from +ve (redder, away) to -ve (bluer,towards)
     else
        cube%vAxis(1) = minVel
     endif

#ifdef MPI
     ix1 = (my_rank) * (cube%nx / (np)) + 1
     ix2 = (my_rank+1) * (cube%nx / (np))
     if (my_rank == (np-1)) ix2 = cube%nx
#endif

     deltaV = minVel * 1.e5/cspeed_sgl

     if(nv .ne. 0) then

        do iv = 1,nv

           deltaV = (cube%vAxis(iv)*1.e5/cSpeed_sgl) ! velocities in fraction of c

           if(writeoutput) then
              write(message,*) "Done ",iv," velocity"
              call tune(6, message)  ! start a stopwatch
           endif

           if(iv .eq. 1) then
              call writeinfo("Filling Octal parameters for first time",TRIVIAL)
              call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule, deltaV,.true.)
           else
              call writeinfo("Filling Octal parameters again",TRIVIAL)
              call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule, deltaV,.false.)
           endif

           cube%intensity(:,:,iv) = makeImageGrid(cube,unitvec,posvec,grid,thisMolecule,itrans,deltaV,nsubpixels) ! 

           if(writeoutput) then
              call tune(6, message)  ! stop a stopwatch
           endif

           cube%converged(:,:,iv) = cube%converged(:,:,cube%nv) ! pretty dodgy hack to save making more variables. Convergence criterion unsound too
           cube%nsubpixels(:,:,iv) = cube%nsubpixels(:,:,cube%nv)! as above - used to record number of subpixels used to sample grid

           if(writeoutput) write(*,*) dble(sum(cube%nsubpixels(:,:,iv))) / dble(npixels**2) , maxval(cube%nsubpixels(:,:,iv))
           intensitysum = sum(cube%intensity(:,:,iv)) / cube%nx**2
           fluxsum = intensitytoflux(intensitysum, dble(imageside), dble(gridDistance), thisMolecule)

           if(iv .eq. 1) then
              background = cube%intensity(1,1,1)
              write(*,*) "background",background
              background = intensitytoflux(background, dble(imageside), dble(gridDistance), thisMolecule)
              write(*,*) "Background Flux: ", background
 !             call GaussianWeighting(cube, cube%nx, beamsize)!, NormalizeArea = .true.)
              
              
!              call fineGaussianWeighting(cube, cube%nx, beamsize, weight)!, NormalizeArea = .true.) 

 !          open(77, file="weight.dat",status="unknown",form="formatted")
 !          open(78, file="fineweight.dat",status="unknown",form="formatted")

 !          do i = 1,npixels
 !             write(77,'(50(f9.5,1x))') cube%Weight(:,i)
 !             write(78,'(50(f9.5,1x))') Weight(:,i)
 !          enddo

 !          close(77)
 !          close(78)
           endif

 !          weightedfluxmap  = cube%intensity(:,:,iv) * cube%weight(:,:)
 !          weightedfluxsum = sum(weightedfluxmap) / cube%nx**2
 !          weightedflux = intensitytoflux(weightedfluxsum, dble(imageside), dble(gridDistance))

!           fineweightedfluxmap  = cube%intensity(:,:,iv) * weight(:,:)
!           fineweightedfluxsum = sum(fineweightedfluxmap) / cube%nx**2
!           fineweightedflux = intensitytoflux(fineweightedfluxsum, dble(imageside), dble(gridDistance), thisMolecule)

           write(message,'(a,es11.4e1,tr3,a,f8.4,tr3,a,es12.4,a,es12.4,es12.4,es12.4)') &
                "DELTAV(v/c):",deltaV," V (km/s):",real(cube%vAxis(iv)), "Average Intensity:",intensitysum, &
                " FLUX: ", fluxsum, (fluxsum / thisMolecule%transfreq(itrans)) * 1e26, (fluxsum - background) &
                / thisMolecule%transfreq(itrans) * 1e26  
           open(10, file="tempfile.dat",status="unknown",form="formatted",position="append")
           write(10,'(es11.4e1,f8.4,es12.4,es12.4,es12.4,es12.4)') &
                real(cube%vAxis(iv)), deltaV, intensitysum, &
                fluxsum / thisMolecule%transfreq(itrans)* 1e26, (fluxsum - background) / thisMolecule%transfreq(itrans) * 1e26 
           close(10)
           call writeinfo(message,FORINFO)
        end do

     else
        iv = 1
        do while (deltaV*cspeed*1.d-5 .lt. maxVel)

           ddV = (maxVel-minVel)/100.
           cube%vAxis(iv) = deltaV * cSpeed*1.d-5

           if(iv .gt. 4) then 
              cube%vAxis(iv) = cube%vAxis(iv-1) + ddV * nextStep(cube, iv) ! NEXTSTEP
              deltaV = cube%vAxis(iv) * 1d5 / cSpeed
           else
              cube%vAxis(iv) = minvel + 4.d0 * ddV * iv          
              deltaV = cube%vAxis(iv) * 1d5 / cSpeed
           endif

           cube%intensity(:,:,iv) = makeImageGrid(cube,unitvec,posvec,grid,thisMolecule,itrans,deltaV,nsubpixels) ! 

           cube%converged(:,:,iv) = cube%converged(:,:,200)
           cube%nsubpixels(:,:,iv) = cube%nsubpixels(:,:,200)

           intensitysum = sum(cube%intensity(:,:,iv)) / cube%nx**2
           fluxsum = intensitytoflux(intensitysum, dble(imageside), dble(gridDistance), thisMolecule)

           if(iv .eq. 1) then

           call GaussianWeighting(cube, cube%nx, beamsize, NormalizeArea = .true.) 
           call fineGaussianWeighting(cube, cube%nx, beamsize, weight, NormalizeArea = .true.) 

           open(77, file="weight.dat",status="unknown",form="formatted")
           open(78, file="fineweight.dat",status="unknown",form="formatted")

           do i = 1,npixels
           write(77,'(50(f9.5,1x))') cube%Weight(:,i)
           write(78,'(50(f9.5,1x))') Weight(:,i)
           enddo

           close(77)
           close(78)
           endif

           weightedfluxmap  = cube%intensity(:,:,iv) * cube%weight(:,:)
           weightedfluxsum = sum(weightedfluxmap) / cube%nx**2
           weightedflux = intensitytoflux(weightedfluxsum, dble(imageside), dble(gridDistance), thisMolecule)

           fineweightedfluxmap  = cube%intensity(:,:,iv) * weight(:,:)
           fineweightedfluxsum = sum(fineweightedfluxmap) / cube%nx**2
           fineweightedflux = intensitytoflux(fineweightedfluxsum, dble(imageside), dble(gridDistance), thisMolecule)

           write(message,'(a,es11.4e1,tr3,a,f8.4,tr3,a,es12.4,a,es12.4,es12.4,es12.4)') &
                "DELTAV(v/c):",real(deltaV)," V (km/s):",real(cube%vAxis(iv)), "Average Intensity:",intensitysum, &
                " nu*FLUX: ", fluxsum * thisMolecule%transfreq(itrans), weightedflux * thisMolecule%transfreq(itrans), &
                fineweightedflux * thisMolecule%transfreq(itrans)
           open(10, file="tempfile.dat",status="unknown",form="formatted",position="append")
           write(10,'(es11.4e1,f8.4,es12.4,es12.4,es12.4,es12.4)') &
                real(cube%vAxis(iv)), real(deltaV), intensitysum, &
                fluxsum * thisMolecule%transfreq(itrans)!, weightedflux * thisMolecule%transfreq(itrans), &  
                !fineweightedflux * thisMolecule%transfreq(itrans)
           close(10)
           call writeinfo(message,FORINFO)

           iv = iv + 1
        enddo
     endif

 !    cube%intensity(:,:,1:cube%nv) = cube%intensity(:,:,1:cube%nv) - cube%intensity(1,1,1) ! TJH background subtract

 !    do i = ix1, ix2
 !       Write(*,*) i
 !	write(*,*) "You are here"
 !       do j = 1, cube%ny
 !
 !          do iMonte = 1, nMonte
 !             if (nMonte > 1) then
 !                call random_number(r)
 !                xVal = cube%xAxis(i) + (r-0.5d0)*(cube%xAxis(2)-cube%xAxis(1))
 !                call random_number(r)
 !!                yVal = cube%yAxis(j) + (r-0.5d0)*(cube%yAxis(2)-cube%yAxis(1))
 !             else
 !                xVal = cube%xAxis(i)
 !                yVal = cube%yAxis(j)
 !             endif
 !             rayPos =  (xval * xProj) + (yval * yProj)
 !             raypos = rayPos + ((-1.d0*distance) * viewVec)
 !             do k = 1, cube%nv
 !                deltaV = cube%vAxis(k)*1.d5/cSpeed
 !                cube%intensity(i,j,k) = intensityAlongRay(rayPos, viewVec, grid, thisMolecule, iTrans, deltaV)
 !             enddo
 !          enddo
 !          cube%intensity(i,j,1:cube%nv) = cube%intensity(i,j,1:cube%nv) / dble(nMonte)
 !       enddo
 !    enddo

#ifdef MPI
     n = (cube%nx*cube%ny*cube%nv)
     allocate(tempArray(1:n), tempArray2(1:n))
     tempArray = reshape(cube%intensity, (/  n /))

      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

        call MPI_ALLREDUCE(tempArray,tempArray2,n,MPI_DOUBLE_PRECISION,&
            MPI_SUM,MPI_COMM_WORLD,ierr)

     cube%intensity = reshape(tempArray2, (/ cube%nx, cube%ny, cube%nv /))
     deallocate(tempArray, tempArray2)
#endif

   end subroutine createimage

 subroutine createFluxSpectra(cube,thismolecule,itrans)

   use input_variables, only : gridDistance, molAbundance
   character(len=80) :: filename
   type(datacube) :: cube
   type(moleculetype) :: thisMolecule
   real(double),allocatable :: fluxspec(:), inspec(:) !, &
!                               convolvedfluxspec(:), brightnessspec(:), antspec(:), weightedspec(:)
   real(double) :: fac, dA, fac2
   integer :: i
   integer :: itrans

!   allocate(weight(1:cube%nx,cube%ny))

!   call fineGaussianWeighting(cube, cube%nx, beamsize, weight)!, NormalizeArea = .true.)

   allocate(fluxspec(1:cube%nv))
!   allocate(convolvedfluxspec(1:cube%nv))
   allocate(inspec(1:cube%nv))
!   allocate(antspec(1:cube%nv))
!   allocate(weightedspec(1:cube%nv))
!   allocate(brightnessspec(1:cube%nv))

   fac =  pi*(cube%telescope%Diameter/2.d0)**2 / (2.d0*kerg) ! Effective Telescope Area / 2k - Flux -> Antenna Temp

   dA = ((cube%xaxis(2) - cube%xaxis(1)) / (gridDistance/1.d10)) ** 2

   fac2 = (1.d0* (cspeed**2 / thisMolecule%transFreq(iTrans)**2) &
                    / (pi * (cube%telescope%diameter / 2.d0)**2 )) / dA

   if(writeoutput) then
      write(filename,'(a,a,i1,a,i1,a,f4.1,a,f4.1,a,es8.2e2,a)') trim(cube%telescope%label),'lineprofilesJ', &
        thisMolecule%itransUpper(itrans)-1,'-',thisMolecule%itransLower(itrans)-1,&
        'beam',cube%telescope%beamsize,'diam',cube%telescope%diameter/100.,'abund',molAbundance,'.dat'

   open(30,file=filename,status="unknown",form="formatted")
   write(30,'(a,tr5,a,tr6,a,tr9,a,tr8,a)') "V (km/s)", "Intensity Wm-2sr-1Hz-1", "Flux Wm-2","Flux in Jy", "Convolved Flux"

 !  continuumflux = cube%flux(1,1,1)

!   convolvedcube = cube

!   call convolveCube(convolvedcube, dble(beamsize))

   do i = 1, cube%nv
      fluxspec(i) = SUM(cube%flux(1:cube%nx,1:cube%ny,i))
 !     convolvedfluxspec(i) = SUM(convolvedcube%flux(1:cube%nx,1:cube%ny,i))
      inspec(i) = SUM(cube%intensity(1:cube%nx,1:cube%ny,i))
!      cube%flux(:,:,i)=cube%flux(:,:,i)*weight(:,:)
!      weightedspec(i) = SUM(cube%flux(1:cube%nx,1:cube%ny,i))
!      antspec(i) = SUM(cube%flux(1:cube%nx,1:cube%ny,i)) * fac
!      brightnessspec(i) = fluxspec(i) / fac2
            
      write(30,'(f6.3,tr8,es13.4,tr6,es13.4,tr8,f11.5,tr8,es13.4)') real(cube%vAxis(i)), inspec(i)/(1e3*(cube%nx**2)), &
                  fluxspec(i), fluxspec(i) / thismolecule%transfreq(itrans)* 1e26 !Jy is not multiplied by frequency!
   enddo
   
   close(30)
endif

   deallocate(inspec, fluxspec)
!   deallocate(convolvedfluxspec, antspec, weightedspec, brightnessspec)

 end subroutine createFluxSpectra

 subroutine cubeIntensityToFlux(cube,thisMolecule,itrans)

   use input_variables, only : gridDistance, npixels, nv
   type(moleculetype) :: thismolecule
   type(datacube) :: cube
   integer :: ipixel,jpixel,iv, itrans
   real(double) :: dx

   dx = cube%xAxis(2) - cube%xAxis(1) ! pixelwidth in torus units
   dx = (dx / (gridDistance*1e-10))**2 ! not in steradians (* 2 pi)

   !write(*,*) "x2", cube%xAxis(2) , "x1", cube%xAxis(1),  "DX", dx
   if(writeoutput) write(*,*) "dx", cube%xAxis(2) - cube%xAxis(1), "griddistance", griddistance 
 !  cube%flux(1:npixels,1:npixels,dx**2) = cube%intensity(1:npixels,1:npixels,1:nv) * (dx**2) ! Flux through solid angle covered by one pixel in ergs/s
   do ipixel = 1,npixels
      do jpixel = 1,npixels
         do iv = 1,nv
            cube%flux(ipixel,jpixel,iv) = cube%intensity(ipixel,jpixel,iv) * dx * 1e-3 * thismolecule%transfreq(itrans)
 !           write(*,*) cube%flux(ipixel,jpixel,iv)
         enddo
      enddo
   enddo

   write(*,*) "DX", dx,  "gridist", griddistance
 !  cube%flux = cube%flux * 1d23
 end subroutine cubeIntensityToFlux

 real(double) pure function IntensityToFlux(intensity,dx,distance,thismolecule) result (flux)

   use input_variables, only : itrans
   real(double), intent(in) :: dx, intensity, distance
   type(moleculetype), intent(in) :: thismolecule

   flux = intensity * ((dx * 1.d10 / distance)**2) * 1e-3 * thisMolecule%transfreq(itrans) ! Flux through solid angle covered by one pixel in ergs/s

 end function IntensityToFlux

 subroutine GaussianWeighting(cube,npixels,FWHM,NormalizeArea)

   use input_variables, only : gridDistance
   type(datacube) :: cube
   integer :: npixels, i, j
 ! real(double),pointer :: PixelWeight(:,:)
 ! type(VECTOR), allocatable :: SubcellCentreArray(:), GridArray(:), PixelPositionArray(:,:)
   real(double) :: ArcSecsToRadians = pi/6.48d5 ! 6.48d5 = 180*60*60
   real :: FWHM
   real(double) :: beamsize,beamsizeInRadians,sigma2,Int,rr
 !  character(len=80) :: message
   logical,optional :: NormalizeArea

 ! allocate(PixelPositionArray(npixels,npixels))

   !! Calculate Gaussian Weighting 

   beamsize = FWHM
   beamsizeInRadians = beamsize * ArcSecsToRadians !(gridDistance/1d10)*(beamsize/(6.48d5/pi)) 
 ! beamSizeInRadians = 1.22d0*(cspeed/thisMolecule%transFreq(itrans))/15.d2 / 2.d0
   sigma2 = (beamSizeInRadians/2.35d0)**2

 !  write(message,*) "Beam size (in radians): ",beamSizeinRadians
 !  call writeinfo(message,TRIVIAL) 

   Int = 0.d0
   cube%Weight = 0.d0

   do i=1,npixels ! Where are all the pixels (top left corners) - what are all their weights
      do j=1,npixels
         rr = (cube%xAxis(i)/(gridDistance*1d-10))**2 + & 
              (cube%yAxis(j)/(gridDistance*1d-10))**2

         cube%Weight(i,j) = exp(-0.5d0*(rr/sigma2))

         Int = Int + cube%Weight(i,j)/dble(npixels)**2

      enddo
   enddo

   if(present(NormalizeArea) .and. NormalizeArea) then
      cube%Weight = cube%Weight / Int ! Normalize to npixels**2 (because implicit uniform instrument function is 1 in each pixel not 1/npixels**2 
   else
      cube%Weight = cube%Weight / maxval(cube%weight(:,:))
   endif
 end subroutine GaussianWeighting

 subroutine fineGaussianWeighting(cube,npixels,FWHM,weight,NormalizeArea)
   use input_variables, only : gridDistance
   type(datacube) :: cube
   integer :: npixels, i, j,k,l
 ! real(double),pointer :: PixelWeight(:,:)
 ! type(VECTOR), allocatable :: SubcellCentreArray(:), GridArray(:), PixelPositionArray(:,:)
   real(double) :: ArcSecsToRadians = pi/6.48d5 ! 6.48d5 = 180*60*60
   real :: FWHM
   real(double) :: beamsize,beamsizeInRadians,sigma2,Int,rr,deltax, deltay,x,y,f,weight(:,:)
   character(len=80) :: message
   logical,optional :: NormalizeArea

   !! Calculate Gaussian Weighting 

   beamsize = FWHM
   beamsizeInRadians = beamsize * ArcSecsToRadians !(gridDistance/1d10)*(beamsize/(6.48d5/pi)) 
 ! beamSizeInRadians = 1.22d0*(cspeed/thisMolecule%transFreq(itrans))/15.d2 / 2.d0
   sigma2 = (beamSizeInRadians/2.35d0)**2

   write(message,*) "Beam size (in radians): ",beamSizeinRadians
   call writeinfo(message,TRIVIAL) 

   Int = 0.d0
   Weight = 0.d0

   deltax = cube%xAxis(3) - cube%xAxis(2)
   deltay = cube%yAxis(3) - cube%yAxis(2)

   do i=1,npixels ! Where are all the pixels (top left corners) - what are all their weights
      do j=1,npixels

         do k = 1,19,2
            do l = 1,19,2
               x = cube%xAxis(i) + real(k-10)/20. * deltax  
               y = cube%yAxis(j) + real(l-10)/20. * deltay

               rr = (x/(gridDistance*1d-10))**2 + (y/(gridDistance*1d-10))**2
               f = exp(-0.5d0*(rr/sigma2))
               Weight(i,j) = Weight(i,j) + f*0.01

            enddo
         enddo

         Int = Int + Weight(i,j)/(dble(npixels)**2)

      enddo
   enddo

   if(present(NormalizeArea) .and. NormalizeArea) then
      Weight = Weight / Int ! Normalize to npixels**2 (because implicit uniform instrument function is 1 in each pixel not 1/npixels**2 
   else
      Weight = Weight / maxval(weight(:,:))
   endif

 end subroutine fineGaussianWeighting

!sub intensity
   subroutine intensityAlongRay(position, direction, grid, thisMolecule, iTrans, deltaV,i0,tau,tautest)

     use input_variables, only : useDust
     type(VECTOR) :: position, direction, dsvector
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(VECTOR) :: currentPosition, thisPosition, endPosition
     type(VECTOR) :: startVel, endVel, thisVel, veldiff
     real(double) :: alphanu1, alphanu2, jnu, snu
     real(double) :: alpha
     real(double)  :: dv, deltaV, dVacrossCell
     integer :: i, icount
     real(double) :: distArray(2000), tval
     integer :: nTau
     real(double) ::  OneOvernTauMinusOne, ds
     real(double) :: nLower, nUpper, balance
     real(double) :: dTau, etaline, dustjnu
     real(double), intent(out), optional :: tau
     real(double),save :: BnuBckGrnd

     real(double) :: phiProfVal

     logical,save :: firsttime = .true.
     logical, optional :: tautest
     logical :: dotautest
     logical :: lowvelgrad = .false.
     logical :: realdust = .true.

     if(present(tautest)) then
        dotautest = tautest
     else
        dotautest = .false.
     endif

     if(firsttime .or. dotautest) then
        BnuBckGrnd = Bnu(thisMolecule%transfreq(itrans), Tcbr)
        if(.not. dotautest) firsttime = .false.
     endif

     if(inOctal(grid%octreeRoot, Position)) then
!        if( dotautest) write(*,*) "inside"
        disttogrid = 0.
     else
        distToGrid = distanceToGridFromOutside(grid, position, direction) 
!        if( dotautest) write(*,*) "outside", disttogrid
     endif

     if (distToGrid > 1.e29) then
        write(*,*) "ray does not intersect grid",position,direction
        i0 = 1.d-60
        tau = 1.d-60
        goto 666
     endif

     currentPosition = position + (distToGrid + 5.d-4*grid%halfSmallestSubcell) * direction
    if(grid%geometry .eq. 'iras04158') currentPosition = position + (distToGrid + 2500.d0*Grid%halfsmallestsubcell) * direction
!     if(grid%geometry .eq. 'iras04158') currentPosition = position + (distToGrid * (1.+1d-5)) * direction

     i0 = 0.d0
     tau = 0.d0

     thisOctal => grid%octreeRoot
     icount = 0

        if(.not. inOctal(grid%octreeRoot, currentPosition) .and. icount .eq. 0) stop

        do while(inOctal(grid%octreeRoot, currentPosition))
           icount = icount + 1

        call findSubcelllocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

        nMol = thisOctal%molcellparam(subcell,1)
        nLower = thisOctal%molcellparam(subcell,2)
        nUpper = thisOctal%molcellparam(subcell,3)
        balance = thisOctal%molcellparam(subcell,4)
        etaline = thisOctal%molcellparam(subcell,5)
        alphanu1 = thisOctal%molcellparam(subcell,6)
        alphanu2 = thisOctal%molcellparam(subcell,7)
        dustjnu = thisOctal%molcellparam(subcell,8)

        thisPosition = currentPosition

        if(grid%geometry .eq. 'iras04158') then
           startVel = keplerianVelocity(currentposition, grid)
           endPosition = currentPosition + tval * direction
           endVel = keplerianVelocity(endposition, grid)
        else
           startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell)
           endPosition = currentPosition + tval * direction
           endVel = amrGridVelocity(grid%octreeRoot, endPosition, startOctal = thisOctal, actualSubcell = subcell)
        endif

        Veldiff = endVel - startVel

        dvAcrossCell = (veldiff.dot.direction)
        dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))

        nTau = min(max(2, nint(dvAcrossCell * 5.d0)), 100) ! ensure good resolution / 5 chosen as its the magic number!

        distArray(1) = 0.d0
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)

        ds = tval * OneOvernTauMinusOne

        dsvector = ds * direction

        do i = 2, nTau 

           thisPosition = thisPosition + dsvector
           thisvel = startvel + real(i-1) * OneOvernTauMinusOne * Veldiff
           dv = (thisVel .dot. direction) - deltaV

           phiProfval = phiProf(dv, thisOctal%molmicroturb(subcell))
           alphanu1 = thisOctal%molcellparam(subcell, 6) * phiprofval

           alpha = alphanu1 + alphanu2
           dTau = alpha * ds * 1.d10

           jnu = etaLine * phiProfVal

           if(useDust) jnu = jnu + dustjnu

           if (alpha .ne. 0.d0) then

              snu = jnu/alpha
           else

              snu = tiny(snu)
              i0 = i0 + tiny(i0)
           endif

           if(tau .lt. 500. .and. dtau .lt. 500.) then
              i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
           else
              i0 = i0
           endif

           tau = tau + dtau

        enddo

        currentPosition = currentPosition + (tval + 1d-3 * grid%halfSmallestSubcell) * direction
     enddo

 666 continue

           if(tau .lt. 500.) then
              i0 = i0 + bnuBckGrnd * exp(-tau) ! from far side
           else
              i0 = i0
           endif
   end subroutine intensityAlongRay

   subroutine tauAlongRay(position, direction, grid, thisMolecule, deltaV, tau)

     use input_variables, only : useDust
     type(VECTOR) :: position, direction
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: nMol
     type(OCTAL), pointer :: thisOctal
     integer :: subcell

     type(VECTOR) :: currentPosition, thisPosition, endPosition
     type(VECTOR) :: startVel, endVel, thisVel, veldiff
     real(double) :: alphanu(maxtrans, 2), alpha(maxtrans)
     real(double)  :: dv, deltaV, dVacrossCell
     integer :: i
     real(double) :: dist, tval
     integer :: nTau
     real(double) :: dds, celledgeinterval
     real(double) :: nLower(maxtrans), nUpper(maxtrans), balance(maxtrans), alphatemp(maxtrans)
     integer :: iLower(maxtrans), iUpper(maxtrans)
     real(double) :: dTau(maxtrans)
     real(double), intent(out) :: tau(maxtrans)

     real(double) :: phiProfVal

     currentPosition = position
     
     tau = 0.d0
     thisOctal => grid%octreeRoot

     iUpper(:)  = thisMolecule%iTransUpper(:)
     iLower(:)  = thisMolecule%iTransLower(:)
     
     do while(inOctal(grid%octreeRoot, currentPosition))
        
        call findSubcelllocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

        nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
        nLower(:)  = thisOctal%molecularLevel(subcell,iLower(:)) * nMol
        nUpper(:)  = thisOctal%molecularLevel(subcell,iUpper(:)) * nMol
        balance(:) = (nLower(:) * thisMolecule%einsteinBlu(:) - &
                     nUpper(:) * thisMolecule%einsteinBul(:))
        alphaTemp(:) = hCgsOverFourPi * balance(:) ! Equation 8
!write(*,*) "nmol", nmol
!write(*,*) "nlower", nlower
!write(*,*) "nupper", nupper
!write(*,*) "balance", balance

        thisPosition = currentPosition
        startVel = Velocity(currentposition, grid)
!        write(*,*) "startvel", startvel
        endPosition = currentPosition + tval * direction
        endVel = Velocity(endposition, grid)
!        write(*,*) "endvel", endvel
        Veldiff = endVel - startVel
!        write(*,*) "veldiff", veldiff
        dvAcrossCell = (veldiff .dot. direction)
!        write(*,*) "dvacrosscell", dvacrosscell
        dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))

        nTau = min(max(2, nint(dvAcrossCell * 5.d0)), 100) ! ensure good resolution / 5 chosen as its the magic number!

        CellEdgeInterval = 1.d0 / (dble(ntau) - 1.0)
        dds = tval * cellEdgeInterval

        dist = 0.d0

        do i = 2, nTau

           dist = dist + dds ! increment along distance counter
              
           thisPosition = currentPosition + dist * direction

           thisVel = velocity(thisPosition, grid)
           dv = deltaV - (thisVel .dot. direction)
           PhiProfVal = phiProf(dv, thisOctal%molmicroturb(subcell))
           
           alphanu(1:maxtrans,1) = phiprofval * alphaTemp(1:maxtrans) ! Equation 8
           alpha(1:maxtrans) = alphanu(1:maxtrans,1) !+ alphanu(1:maxtrans,2)              
!        write(*,*) "alpha", alpha
!        write(*,*) dds
           dTau(:) = alpha(:) * dds * 1.d10 ! dds is interval width & optical depth, dTau = alphanu*dds - between eqs (3) and (4)
              
           tau(:) = tau(:) + dtau(:) ! contribution to optical depth from this line integral
           
        enddo
        
        currentPosition = currentPosition + (tval + 1.d-3*grid%halfSmallestSubcell) * direction ! FUDGE - make sure that new position is in new cell
     enddo

   end subroutine tauAlongRay

   subroutine lteintensityAlongRay(position, direction, grid, thisMolecule, iTrans, deltaV,i0,tau,tautest)

     use input_variables, only : useDust
     type(VECTOR) :: position, direction, dsvector
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0
     type(OCTAL), pointer :: thisOctal, startOctal
     integer :: subcell
     type(VECTOR) :: currentPosition, thisPosition, thisVel
     type(VECTOR) :: rayVel, startVel, endVel, endPosition
     real(double) :: alphanu(2), jnu, snu
     real(double) :: alpha
     integer, save :: iLower , iUpper
     real(double) :: dv, deltaV, dvAcrossCell
     integer :: i, icount
     real(double) :: distArray(200), tval
     integer :: nTau
     real(double) ::  OneOvernTauMinusOne, ds
     real(double) :: nLower, nUpper, balance, phiProfVal
     real(double) :: dTau, etaline, kappaAbs
     real(double), intent(out), optional :: tau
     real(double),save :: BnuBckGrnd
     real :: lambda
     integer :: ilambda

     logical,save :: firsttime = .true.
     logical, optional :: tautest
     logical :: dotautest
     logical :: lowvelgrad = .false.
     logical :: realdust = .true.

     if(present(tautest)) then
        dotautest = tautest
     else
        dotautest = .false.
     endif

     if(firsttime .or. dotautest) then
        BnuBckGrnd = Bnu(thisMolecule%transfreq(itrans), Tcbr)
        iUpper = thisMolecule%iTransUpper(iTrans)
        iLower = thisMolecule%iTransLower(iTrans)
 !       write(*,*) iupper,ilower
        if(.not. dotautest) firsttime = .false.
     endif

     if(useDust) then
        lambda = (cspeed_sgl * 1e8) / thisMolecule%transfreq(iTrans) ! cm to Angstroms
        lambda = lambda * (1.d0 - deltaV) ! when dv +ve wavelength gets shorter!
        call locate(grid%lamArray, size(grid%lamArray), lambda, ilambda)
!        write(*,*) lambda,ilambda
     endif

     if(inOctal(grid%octreeRoot, Position)) then
        disttogrid = 0.
     else
        distToGrid = distanceToGridFromOutside(grid, position, direction) 
     endif

     if (distToGrid > 1.e29) then
        write(*,*) "ray does not intersect grid",position,direction
        i0 = 1.d-60
        tau = 1.d-60
        goto 666
     endif

    currentPosition = position + (distToGrid + 5000.d0*grid%halfSmallestSubcell) * direction
 !    currentPosition = position + (distToGrid + 5.d-2*grid%halfSmallestSubcell) * direction
     i0 = 0.d0
     tau = 0.d0
     rayVel = VECTOR(0.d0, 0.d0, 0.d0)

     thisOctal => grid%octreeRoot
     icount = 0

     if(.not. inOctal(grid%octreeRoot, currentPosition) .and. icount .eq. 0) stop

     do while(inOctal(grid%octreeRoot, currentPosition))
        icount = icount + 1

        call findSubcelllocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

!        if(.not. thisOctal%done) then

!           thisOctal%molmicroturb(subcell) = 1.d0 / thisOctal%microturb(subcell)

!           allocate(thisOctal%molcellparam(subcell,8))
!           thisOctal%molcellparam(subcell,1) = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
           
!           nMol = thisOctal%molcellparam(subcell,1)

!           thisOctal%molcellparam(subcell,2) = thisOctal%molecularLevel(subcell,iLower) * nMol
!           thisOctal%molcellparam(subcell,3) = thisOctal%molecularLevel(subcell,iUpper) * nMol
           
!           nLower = thisOctal%molcellparam(subcell,2)
!           nUpper = thisOctal%molcellparam(subcell,3)
!           thisOctal%molcellparam(subcell,4) = nLower * thisMolecule%einsteinBlu(iTrans) - nUpper * thisMolecule%einsteinBul(iTrans)

!           etaLine = hCgsOverFourPi * thisMolecule%einsteinA(iTrans)
!           thisOctal%molcellparam(subcell,5) = etaLine * nUpper
        
!           thisOctal%done = .true.

!           thisOctal%molcellparam(subcell,6) = nTau
!           thisOctal%molcellparam(subcell,7) = OneOvernTauMinusOne
!        endif

        nMol = thisOctal%molcellparam(subcell,1)
        nLower = thisOctal%molcellparam(subcell,2)
        nUpper = thisOctal%molcellparam(subcell,3)
        balance = thisOctal%molcellparam(subcell,4)
        etaline = thisOctal%molcellparam(subcell,5)
       
!        write(*,*) thisOctal%molmicroturb(subcell), thisOctal%microturb(subcell)


 !       write(*,*) dvAcrossCell, OneOvernTauMinusOne, startvel, endvel, thisOctal%microturb(subcell)                                                                                                                                                                                

        if(usedust) then
!           if(realdust) then
              call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = lambda, kappaAbs = kappaAbs)
!           else
!             kappaAbs = thisOctal%rho(subcell) * 0.1 * thisMolecule%transfreq(itrans) * 1e-12 * 1d20 !multiplied by density !cm -> torus
!              kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10 !multiplied by density !cm -> torus
!           endif
           alphanu(2) = kappaAbs * 1.d-10 ! * thisOctal%rho(subcell)
        else
           alphanu(2) = 0.d0
        endif

        ds = tval * OneOvernTauMinusOne
        dsvector = ds * direction
        thisPosition = currentPosition

        startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell)
        endPosition = currentPosition + tval * direction
        endVel = amrGridVelocity(grid%octreeRoot, endPosition)

        dvAcrossCell = ((startVel - endvel).dot.direction)

!        write(*,*) "before",dvAcrossCell

        dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))

!        write(*,*) "afeter",dvAcrossCell

        nTau = min(max(2, nint(dvAcrossCell * 20.d0)), 100) ! ensure good resolution / selects dVacrossCell as being between 0.1 and 10 (else nTau bounded by 2 and 200)

        distArray(1) = 0.d0
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)

        do i = 2, nTau 

           startOctal => thisOctal
           thisPosition = thisPosition + dsvector

!           if(lowvelgrad) then 
              thisvel = startvel + (i-1) * OneOvernTauMinusOne * (endVel - startVel)
!           else
!              thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell)
!           endif

!              write(*,*) modulus(altthisvel-thisvel),thisvel%x-altthisvel%x,thisvel%y-altthisvel%y,thisvel%z-altthisvel%z

           dv = (thisVel .dot. direction) - deltaV

           phiProfval = phiProf(dv, thisOctal%molmicroturb(subcell))

!           write(*,*) phiProfVal

           alphanu(1) = hCgsOverfourPi * phiprofval
           alphanu(1) = alphanu(1) * balance

           alpha = alphanu(1) + alphanu(2)
           dTau = alpha * ds * 1.d10

           jnu = etaLine * phiProfVal

           if(useDust) jnu = jnu + alphanu(2) * thisOctal%bnu(subcell,itrans)

           if (alpha .ne. 0.d0) then

 !             snu = jnu/alpha
              snu = thisOctal%bnu(subcell,itrans)
!              if((i0 .eq. 0.) .or. tau .lt. 10. .or. exp(-tau) * (1.d0-exp(-dtau))*snu .gt. i0*(1d-5)) then
                 i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
!              else
!                 write(*,*) "i0", i0, "tau", tau, "skipping"
!                 goto 666
!              endif
           else
              snu = tiny(snu)
              i0 = i0 + tiny(i0)
           endif
           tau = tau + dtau
        enddo
        
        currentPosition = currentPosition + (tval+1d-3*grid%halfSmallestSubcell) * direction

!        if((grid%geometry .eq. 'IRAS04158') .and. (modulus(currentposition) < 0.5*rinner)) then
!           currentPosition = currentPosition + (rinner/(VECTOR(1d-20,1d-20,1.d0) .dot. direction)) * direction
!        endif
        
!        if(debug) write(*,*) i0,tau
     enddo

 666 continue

     i0 = i0 + bnuBckGrnd * exp(-tau) ! from far side                                                                                                                                                                                                                                
!     if(debug) write(*,*) i0,tau
 !    i0 = i0 - bnu(thisMolecule%transFreq(iTrans), Tcbr)                                                                                                                                                                                                                            

 !    if (debug) write(*,*) "bound condition",i0,bnu(thisMolecule%transFreq(iTrans), Tcbr) * exp(-tau)                                                                                                                                                                               

     ! convert to brightness T                                                                                                                                                                                                                                                       

 !    i0 = i0 * (cSpeed**2 / (2.d0 * thisMolecule%transFreq(iTrans)**2 * kerg))                                                                                                                                                                                                      

   end subroutine lteintensityAlongRay

   subroutine continuumIntensityAlongRay(position, direction, grid, wavelength, deltaV, i0, tau, tautest)

     type(VECTOR) :: position, direction
     type(GRIDTYPE) :: grid
     real(double) :: disttoGrid
     real(double), intent(out) :: i0
     type(OCTAL), pointer :: thisOctal, startOctal
     integer :: subcell
     type(VECTOR) :: currentPosition, thisPosition, thisVel
     type(VECTOR) :: startVel, endVel, endPosition
     real(double) :: alphanu(2), jnu, snu
     real(double) :: alpha
     real(double) :: dv, deltaV
     integer :: i, icount
     real(double) :: distArray(20000), tval
     integer :: nTau
     real(double) ::  OneOvernTauMinusOne
     real(double) :: dTau, kappaAbs, kappaSca
     real(double), intent(out), optional :: tau

     real(double) :: dvAcrossCell
     real :: lambda, wavelength
     integer :: ilambda

     logical, optional :: tautest
     logical :: dotautest
     logical :: lowvelgrad = .false.

     if(present(tautest)) then
        dotautest = tautest
     else
        dotautest = .false.
     endif

     lambda = wavelength * (1.d0 - deltaV) ! when dv +ve wavelength gets shorter!
!     write(*,*) "lambda", lambda
     
     call locate(grid%lamArray, size(grid%lamArray), lambda, ilambda)
     
     if(inOctal(grid%octreeRoot, Position)) then
        disttogrid = 0.
!        write(*,*) "inside"
     else
        distToGrid = distanceToGridFromOutside(grid, position, direction) 
       ! write(*,*) "outside", position
       ! write(*,*) "outside", direction
       ! write(*,*) "disttogrid", disttogrid
     endif

     if (distToGrid > 1.e29) then
        write(*,*) "ray does not intersect grid",position,direction
        i0 = 1.d-60
        tau = 1.d-60
        goto 666
     endif

    currentPosition = position + (distToGrid + 5.d-2*grid%halfSmallestSubcell) * direction
    if(grid%geometry .eq. 'iras04158') currentPosition = position + (distToGrid + 2500.d0*Grid%halfsmallestsubcell) * direction
     !   write(*,*) "outside", position
     !   write(*,*) "outside", direction
     !   write(*,*) "disttogrid", disttogrid
 !    currentPosition = position + (distToGrid + 5.d-2*grid%halfSmallestSubcell) * direction
     i0 = 0.d0
     tau = 0.d0
    
     thisOctal => grid%octreeRoot
     icount = 0

     do while(inOctal(grid%octreeRoot, currentPosition))
        icount = icount + 1

        call findSubcelllocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

        startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell)
        endPosition = currentPosition + tval * direction
        endVel = amrGridVelocity(grid%octreeRoot, endPosition)

 !       dvAcrossCell = ((startVel - rayVel).dot.direction) - ((endVel - rayVel).dot.direction)
         dvAcrossCell = ((startVel - endvel).dot.direction)
        dvAcrossCell = abs(dvAcrossCell / thisOctal%microturb(subcell))

        !nTau = min(5, nint(dvAcrossCell*10.d0)) !!! CHECK                                                                                                                                                                                                                           
        nTau = min(max(2, nint(dvAcrossCell * 20.d0)), 200) ! ensure good resolution / selects dVacrossCell as being between 0.1 and 10 (else nTau bounded by 2 and 200)                                                                                                             
!        ntau = 100

        distArray(1) = 0.d0
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)
 !       write(*,*) dvAcrossCell, OneOvernTauMinusOne, startvel, endvel, thisOctal%microturb(subcell)                                                                                                                                                                                
        do i = 2, nTau

           distArray(i) = tval * dble(i-1) * OneOvernTauMinusOne

           startOctal => thisOctal
           thisPosition = currentPosition + distArray(i)*direction

        if(grid%geometry .eq. 'iras04158') then
           thisVel = keplerianVelocity(thisPosition, grid)
        else
           thisvel = startvel + (i-1) * OneOvernTauMinusOne * (endvel - startvel)
        endif

           dv = (thisVel .dot. direction) - deltaV

           call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = lambda, kappaAbs = kappaAbs, kappaSca = kappaSca)

           alphanu(2) = (kappaSca + kappaAbs) * 1.d-10  !* thisOctal%rho(subcell)

           alpha = alphanu(2)
           dTau = alpha * distArray(2) * 1.d10

           jnu = alphanu(2) * bnu(cspeed/lambda, dble(thisOctal%temperaturedust(subcell)))

           if (alpha .ne. 0.d0) then

              snu = jnu/alpha
 !             snu = thisOctal%bnu(subcell,itrans)
              i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu

           else
              snu = tiny(snu)
              i0 = i0 + tiny(i0)
           endif
           tau = tau + dtau
!           if(writeoutput) write(*,*) "i",i,"tau at",thisposition%x,"is",tau
        enddo

        currentPosition = currentPosition + (tval+1.d-6*grid%halfSmallestSubcell) * direction

     enddo

 666 continue

     i0 = i0 + bnu(cspeed/lambda, dble(tcbr))  * exp(-tau) ! from far side                                                                                                                                                                                                                                

   end subroutine continuumIntensityAlongRay

   recursive subroutine calculateOctalParams(grid, thisOctal, thisMolecule, deltaV, firsttime)

     use input_variables, only : useDust, iTrans

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     type(VECTOR)  :: rvec
     integer :: subcell, i
     integer :: iupper, ilower
     integer, save :: ilambda = 1
     real(double) :: nmol, nlower, nupper, etaline, kappaAbs, deltaV
     real :: lambda
     logical,save :: onetime = .true.
     logical :: firsttime
     character(len=10) :: out

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren, 1
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call calculateOctalParams(grid, child, thisMolecule, deltaV, firsttime)
                 exit
              end if
           end do
        else
          if(firsttime) then

             if(grid%geometry .eq. "iras04158") then
                 thisOctal%microturb(subcell) = max(3d-7,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / &
                      (28.0 * amu)) + 0.3**2 ) / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence

                 thisOctal%velocity(subcell) = keplerianVelocity(subcellcentre(thisOctal,subcell), grid)
                 CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)

                 out = 'abundance'
                 rvec = subcellCentre(thisOctal,subcell)
                 thisOctal%molabundance(subcell) = readparameterfrom2dmap(rvec,out,.true.)
!                 thisOctal%molabundance(subcell) = 1d-4
!                 thisOctal%nh2(subcell) = thisOctal%nh2(subcell) * 1.05
!                 thisOctal%rho(subcell) = thisOctal%rho(subcell) * 1.05
              endif
!              else
!                 thisOctal%microturb(subcell) = max(1d-8,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / &
!                      (29.0 * amu)) + 0.3**2) / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence
!                 thisOctal%velocity(subcell) = keplerianVelocity(subcellcentre(thisOctal,subcell), grid)
!                 CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
!              endif

              thisOctal%molmicroturb(subcell) = 1.d0 / thisOctal%microturb(subcell)

              if (.not.associated(thisOctal%molcellparam)) then
                 allocate(thisOctal%molcellparam(1:thisOctal%maxChildren,8))
              endif

              thisOctal%molcellparam(subcell,1) = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
              
              nMol = thisOctal%molcellparam(subcell,1)

              iUpper = thisMolecule%iTransUpper(iTrans)
              iLower = thisMolecule%iTransLower(iTrans)

              thisOctal%molcellparam(subcell,2) = thisOctal%molecularLevel(subcell,iLower) * nMol
              thisOctal%molcellparam(subcell,3) = thisOctal%molecularLevel(subcell,iUpper) * nMol
              
              nLower = thisOctal%molcellparam(subcell,2)
              nUpper = thisOctal%molcellparam(subcell,3)
              thisOctal%molcellparam(subcell,4) = nLower * thisMolecule%einsteinBlu(iTrans) &
                                                - nUpper * thisMolecule%einsteinBul(iTrans)
              
              etaLine = hCgsOverFourPi * thisMolecule%einsteinA(iTrans)
              
              thisOctal%molcellparam(subcell,5) = etaLine * nUpper
              thisOctal%molcellparam(subcell,6) = hCgsOverFourPi * thisOctal%molcellparam(subcell,4)! balance
           endif

           if(usedust) then

              lambda = (cspeed_sgl * 1e8) / thisMolecule%transfreq(iTrans) ! cm to Angstroms
              lambda = lambda * (1.d0 - deltaV) ! when dv +ve wavelength gets shorter
              if(onetime) then
                 call locate(grid%lamArray, size(grid%lamArray), lambda, ilambda)
                 onetime = .false.
              endif
              call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = lambda, kappaAbs = kappaAbs)
              thisOctal%molcellparam(subcell,7) = kappaAbs * 1.d-10 ! * thisOctal%rho(subcell)
              thisOctal%molcellparam(subcell,8) = thisOctal%molcellparam(subcell,7) * thisOctal%bnu(subcell,itrans)
           endif
        endif

     end do

   end subroutine calculateOctalParams


 ! Equation 9 - calculating the doppler broadening due to a turbulent velocity field (b - constant locally)
 ! NB b here = b*nu0/c in the paper ! actually since April 2008 b = 1/old b to speed up code

   real(double) pure function phiProf(dv, b) result (phi)
     real(double), intent(in):: dv, b
     real(double) :: fac

     fac = min((dv * b)**2,200.d0) ! avoid floating point underflows
     phi = exp(-fac) * (b * OneOversqrtPi)
   end function phiProf

    ! sample level populations at logarithmically spaced annuli
   subroutine dumpResults(grid, thisMolecule)!, convtestarray)
     use input_variables, only : rinner, router,amr2d
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: r, ang
     real(double), save :: logrinner, logrouter, OneOverNradiusMinusOne, OneOverNangle
     integer :: i
     real(double) :: pops(minlevel), fracChange(minlevel), tauarray(40)!, convtestarray(:,:,:), tauarray(40)
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     integer :: j
     real(double) :: x, z
     type(VECTOR) :: posvec
     integer :: iter
     integer, save :: nradius, nangle
     character(len=30) :: resultfile, resultfile2
     logical, save :: firsttime = .true.

     iter = grand_iter
!     write(*,*) iter, grand_iter,nray

     write(resultfile,'(a,I7.7)') "results.", nRay
     write(resultfile2,'(a,I2)') "fracChangeGraph.", iter
     
     open(31,file=resultfile,status="unknown",form="formatted")

     open(32,file=resultfile2,status="unknown",form="formatted")

!     if( .not. amr2d) &
          call taualongray(VECTOR(1.d-10,1.d-10,1.d-10), VECTOR(1.d0, -1.d-20, -1.d-20),&
	  grid, thisMolecule, 0.d0, tauarray(1:maxtrans))

     if(firsttime) then
        open(33,file="tau.dat",status="unknown",form="formatted")	

        nradius = 100! number of radii used to sample distribution
        nangle = 360 ! number of rays used to sample distribution
        
        logrinner = log10(rinner)
        logrouter = log10(router)
        OneOverNradiusMinusOne = 1.d0/dble(nradius - 1)
        OneOverNangle = 1.d0 / dble(nangle)
        
        firsttime = .false.
     endif

     thisOctal => grid%octreeroot

     do i = 1, nradius
        r = logrinner + dble(i - 1) * OneOverNradiusMinusOne * (logrouter - logrinner)
        r = 10.d0**r
        pops = 0.d0
        fracChange = 0.d0

        do j = 1 , nangle
           ang = dble(j + 0.5d0) * OneOverNangle - OneOverNangle
           ang = ang * twopi 
           z = r * cos(ang)
           x = r * sin(ang)
           posVec = VECTOR(x, 0.d0, z) ! get put in octal on grid
           call findSubcellLocal(posVec, thisOctal,subcell) 
           pops = pops + thisOctal%molecularLevel(subcell,1:minlevel) ! interested in first 8 levels
           fracChange = fracChange + abs((thisOctal%molecularLevel(subcell,1:minlevel) - &
                        thisOctal%oldmolecularLevel(subcell,1:minlevel)) / thisOctal%molecularlevel(subcell,1:minlevel))           
        enddo

        pops = pops * OneOverNangle ! normalised level population at the 20 random positions 
        fracChange = fracChange * OneOverNangle

!        if(iter .gt. 0) convtestarray(iter,i,1:minlevel-1) = fracChange(1:minlevel-1)

        write(31,'(es11.5e2,4x,14(es14.6e2,tr2))') r*1.e10, pops(1:minlevel)
        write(32,'(es11.5e2,4x,14(f7.5,tr2))') r*1.e10, fracChange(1:minlevel)
     enddo
     write(33,'(i2,tr3,12(f11.6,tr3))') grand_iter, tauarray(1:mintrans)
     
     close(31)
     close(32)
   end subroutine dumpResults

   recursive subroutine  findmaxlevel(grid, thisOctal, thisMolecule, maxinterestinglevel, nlevels, nVoxel)

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i, nlevels
     real(double) :: mollevels(nlevels)
     integer :: maxinterestinglevel, nVoxel
     integer, save :: icount = 0
     real(double) :: maxtemp = 0.d0
     character(len = 80) :: message
         
     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call findmaxlevel(grid, child, thisMolecule, maxinterestinglevel, nlevels, nVoxel)
                 exit
              end if
           end do
        else
!           write(*,*) "temp", thisoctal%temperature(subcell),subcell

           maxtemp = max(maxtemp, thisOctal%temperature(subcell))     
           icount = icount + 1
        endif
     end do

     if(icount .eq. nVoxel) then

        write(message, *)  "Maximum Global Temperature", maxtemp
        call writeinfo(message, FORINFO)

        call LTEpops(thisMolecule, maxtemp, mollevels(1:thisMolecule%nlevels))
                
        do i = thisMolecule%nlevels, 1, -1
           maxinterestinglevel = i
!           write(*,*) i, mollevels(i)
           if(mollevels(i) .gt. 1d-8) exit
        enddo
        
        maxinterestinglevel = maxinterestinglevel
     endif
     
   End subroutine findmaxlevel

   real(double) pure function DustModel1(freq) result (kappa) !! GGTAU paper
      
     real(double), intent(in) :: freq
     real(double) :: lambda

     lambda = cspeed * 1d-2 / freq

     if(lambda .gt. 250e-6) then
        kappa = 1.6e8 * lambda**2
     else
        kappa = 4.8159e5 * lambda**1.3
     endif
      
   end function DustModel1

   real function nextStep(cube, ivplusone) result (step)

     type(DATACUBE) :: cube
     integer :: iv,i,j,ivplusone
     real(double), allocatable :: x(:), y(:), y2(:), intensitysum(:)
     real(double) :: xquad(3),yquad(3),xa(3)
     real(double) :: sumx(4),sumxy(4)
     real(double) :: a(3,3),b(3)
     real(double) :: fac
     real(double), save :: alty2 !, ydash
     real(double), save :: oldstep = 4.d0

     iv = ivplusone - 1

     allocate(x(iv))
     allocate(y(iv))
     allocate(y2(iv))
     allocate(intensitysum(iv))

     do i = 1,iv
        intensitysum(i) = SUM(cube%intensity(1:cube%nx,1:cube%ny,i))
     enddo

     fac = intensitysum(1)
     x = cube%vAxis(1:iv)
     y =  intensitysum / fac 

 !    ydash = ((y(iv) - y(iv-1)) / (x(iv) - x(iv-1))) + (y2(iv-1) * (x(iv) - x(iv-1)))

 !    write(*,*) "x",x
 !    write(*,*) "y",y

     call spline(x,y,iv,1.d31,1.d31,y2)

 !    write(*,*) y2

     xquad = x(iv-3:iv-1)
     yquad = y2(iv-3:iv-1)
     xa = xquad

     do i=1,4
        sumx(i) = sum(xa)
        sumxy(i) = sum(xa * yquad)

        xa = xa*xquad
     enddo

     do i = 1,3
       do j = 1,3
          if(6-i-j .ne. 0) then 
             a(i,j) = sumx(6-i-j)
          else
             a(i,j) = 3.d0
          endif
       enddo; enddo

       b = (/sumxy(2),sumxy(1),sum(yquad)/)
       call luSlv(a,b,3)

         write(*,*) b(1)*x(iv)**2+b(2)*x(iv)+b(3)

         alty2 = b(1)*x(iv)**2+b(2)*x(iv)+b(3)

          step = (oldstep + min(max(25.d0/abs(alty2),0.5),4.d0))/2.d0
          oldstep = step

 !      write(*,*) "x",xquad,"y2",yquad
 !      write(*,*) "quad", alty2, min(max(25.d0/abs(alty2),0.5),4.d0), oldstep, (oldstep + min(max(25.d0/abs(alty2),0.5),4.d0))/2.d0


   end function nextStep

#ifdef MPI 
       subroutine packMoleLevel(octalArray, nTemps, tArray, octalsBelongRank, iLevel)
     include 'mpif.h'
         type(OCTALWRAPPER) :: octalArray(:)
         integer :: octalsBelongRank(:)
         integer :: nTemps
         real(double) :: tArray(:)
         integer :: iOctal, iSubcell, my_rank, ierr
         integer :: iLevel
         type(OCTAL), pointer :: thisOctal

        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
        !
        ! Update the edens values of grid computed by all processors.
        !
        nTemps = 0
        do iOctal = 1, SIZE(octalArray)

           thisOctal => octalArray(iOctal)%content

           do iSubcell = 1, thisOctal%maxChildren
               if (.not.thisOctal%hasChild(iSubcell)) then
                  nTemps = nTemps + 1
                  if (octalsBelongRank(iOctal) == my_rank) then
                    tArray(nTemps) = thisOctal%newMolecularLevel(isubcell,iLevel)
                  else 
                    tArray(nTemps) = 0.d0
                  endif
               endif
           end do
        end do
      end subroutine packMoleLevel

      subroutine unpackMoleLevel(octalArray, nTemps, tArray, octalsBelongRank, iLevel)
        include 'mpif.h'
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
        integer :: nTemps
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr
        integer :: iLevel
        type(OCTAL), pointer :: thisOctal
        
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
        !
        ! Update the edens values of grid computed by all processors.
        !
        nTemps = 0
        do iOctal = 1, SIZE(octalArray)
           
           thisOctal => octalArray(iOctal)%content
           
           do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 thisOctal%newMolecularLevel(isubcell,iLevel) = tArray(nTemps) 
              endif
           end do
        end do
      end subroutine unpackMoleLevel
#endif
      
      
      
      SUBROUTINE sobseq(x,init)
        USE nrtype; USE nrutil, ONLY : nrerror
        IMPLICIT NONE
        REAL(double), DIMENSION(:), INTENT(OUT) :: x
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: init
        INTEGER(I4B), PARAMETER :: MAXBIT=30,MAXDIM=6
        REAL(SP), SAVE :: fac
        INTEGER(I4B) :: i,im,ipp,j,k,l
        INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE:: iu
        INTEGER(I4B), SAVE :: in
        INTEGER(I4B), DIMENSION(MAXDIM), SAVE :: ip,ix,mdeg
        INTEGER(I4B), DIMENSION(MAXDIM*MAXBIT), SAVE :: iv
        DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
        DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
        if (present(init)) then
           ix=0
           in=0
           if (iv(1) /= 1) RETURN
           fac=1.0_sp/2.0_sp**MAXBIT
           allocate(iu(MAXDIM,MAXBIT))
           iu=reshape(iv,shape(iu))
           do k=1,MAXDIM
              do j=1,mdeg(k)
                 iu(k,j)=iu(k,j)*2**(MAXBIT-j)
              end do
              do j=mdeg(k)+1,MAXBIT
                 ipp=ip(k)
                 i=iu(k,j-mdeg(k))
                 i=ieor(i,i/2**mdeg(k))
                 do l=mdeg(k)-1,1,-1
                    if (btest(ipp,0)) i=ieor(i,iu(k,j-l))
                    ipp=ipp/2
                 end do
                 iu(k,j)=i
              end do
           end do
           iv=reshape(iu,shape(iv))
           deallocate(iu)
        else
           im=in
           do j=1,MAXBIT 
              if (.not. btest(im,0)) exit
              im=im/2
           end do
           if (j > MAXBIT) call nrerror('MAXBIT too small in sobseq')
           im=(j-1)*MAXDIM
           j=min(size(x),MAXDIM)
           ix(1:j)=ieor(ix(1:j),iv(1+im:j+im))
           x(1:j)=ix(1:j)*fac
           in=in+1
        end if
      END SUBROUTINE sobseq

      subroutine setObserverVectors(inc, viewvec, unitvec, posvec, centreVec, gridsize)
        
        type(VECTOR) :: unitvec, viewvec, posvec, centreVec
        real(double) :: farAway, gridsize
        real :: inc
        
        farAway = 500.0 * gridsize
        
        
        if(inc .ge. 0. .and. inc .le. 360.) then
           unitvec = VECTOR(sin(inc*degtorad),1.d-10,cos(inc*degtorad))
        else
           unitvec = randomunitVector()
        endif
        
        posvec = faraway * unitvec
        centrevec = VECTOR(0d7,0d7,0d7)
        viewvec = centrevec - posvec
        call normalize(viewvec)
        unitvec = (-1.d0) * viewvec
        
      end subroutine setObserverVectors
      
      recursive subroutine  findtempdiff(grid, thisOctal, thisMolecule, mean, icount)
        use input_variables, only : rinner, router
        type(GRIDTYPE) :: grid
        type(MOLECULETYPE) :: thisMolecule
        type(octal), pointer   :: thisOctal
        type(octal), pointer  :: child 
        integer :: subcell, i
        integer :: icount
        real(double) :: mean(6),r
        type(VECTOR) :: rvec
        character(len=10) :: out
        
        do subcell = 1, thisOctal%maxChildren
           if (thisOctal%hasChild(subcell)) then
              ! find the child
              do i = 1, thisOctal%nChildren, 1
                 if (thisOctal%indexChild(i) == subcell) then
                    child => thisOctal%child(i)
                    call findtempdiff(grid, child, thisMolecule, mean, icount)
                    exit
                 end if
              end do
           else
              rvec = subcellcentre(thisOctal, subcell)
              r = sqrt(rvec%x**2 + rvec%y**2)
              if(icount .ge.  1) then
                 if(r .lt. router .and. r .gt. rinner) then
                    if(icount .eq. 2) then 
                       out='td'
                       thisOctal%temperature(subcell) = readparameterfrom2dmap(rvec,out,.false.)
                       out = 'nh2'           
                       thisOctal%rho(subcell) = readparameterfrom2dmap(rvec,out,.true.) * 2. * mHydrogen
                       out = 'abundance'
                       thisOctal%molAbundance(subcell) = readparameterfrom2dmap(rvec,out,.true.)
                    else
                       out = 'abundance'
                       thisOctal%molAbundance(subcell) = readparameterfrom2dmap(rvec,out,.true.)
                    endif
                 else
                    thisOctal%temperature(subcell) = tcbr
                    thisOctal%rho(subcell) = 1d-20
                    thisOctal%molAbundance(subcell) = 1d-20
                 endif
              else
                 
                 out='td'
                 thisOctal%temperaturegas(subcell) = readparameterfrom2dmap(rvec,out,.false.)
                 thisOctal%temperaturedust(subcell) = thisOctal%temperature(subcell)
                 thisOctal%temperature(subcell) = thisOctal%temperaturegas(subcell) / thisOctal%temperaturedust(subcell)
                 
!           mean(1) = mean(1) + thisOctal%temperaturegas(subcell)
!           mean(2) = mean(2) + thisOctal%temperaturedust(subcell)
!           mean(3) = mean(3) + thisOctal%temperature(subcell)
                 out = 'nh2'
                 
                 thisOctal%temperaturegas(subcell) = readparameterfrom2dmap(rvec,out,.true.)
                 
!           mean(4) = mean(4) + thisOctal%rho(subcell)
       
                 thisOctal%nh2(subcell) = (thisOctal%temperaturegas(subcell) * 2. *mhydrogen) / thisOctal%rho(subcell)
                 thisOctal%rho(subcell) = thisOctal%nh2(subcell)
!           mean(5) = mean(5) + thisOctal%temperaturegas(subcell)* 2. *mhydrogen
!           mean(6) = mean(6) + thisOctal%nh2(subcell)
       
!           icount = icount+1

                 out = 'abundance'
                 thisOctal%molAbundance(subcell) = readparameterfrom2dmap(rvec,out,.false.)

!           thisOctal%velocity(subcell) = keplerianvelocity(rvec,grid)
!           CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
!           write(*,*) thisOctal%nh2(subcell)

              end if
           endif
        end do

      end subroutine findtempdiff

!!!subroutine testOpticalDepth
      subroutine testOpticalDepth(grid,thisMolecule)

        use input_variables, only : lamstart, lamend, rinner, router, debug
        type(GRIDTYPE) :: grid
        type(MOLECULETYPE) :: thisMolecule
        type(octal), pointer   :: thisOctal
        type(VECTOR) :: currentposition(3), posvec, viewvec, unitvec, centrevec
        integer :: subcell, i, itrans
        integer :: ilamb, nlamb
        real(double) :: xmidplane, gridsize
        real :: lamb
        real(double) :: tau, dummy, kappaAbs, kappaSca, i0
        character(len=50) :: message
  
  call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule,0.d0,.true.)
  
  write(message,*) "Angular dependence"
  call writeinfo(message, FORINFO)
  
  gridsize = grid%octreeroot%subcellsize
  nlamb = 500

  do i = 1, 90

     call setObserverVectors(real(i), viewvec, unitvec, posvec, centreVec, gridsize)
     call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
     
     write(message,*) i, tau
     call writeinfo(message, FORINFO)
  enddo
  
  write(message,*) "Tau from above"
  call writeinfo(message, FORINFO)
  
  do i = 1, 100
     
     unitvec = VECTOR(1.d-20,1.d-20,1.d0)
     posvec = VECTOR(real(i) * grid%octreeroot%subcellsize * 0.01,0d7,2.d0 * grid%octreeroot%subcellsize)
     viewvec = (-1.d0) * unitvec
     
     call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
     
     write(message,*) i, tau
     call writeinfo(message, FORINFO)
  enddo
  
  write(message,*) "Tau from the side"
  call writeinfo(message, FORINFO)
  
  do i = 1, 100
     unitvec = VECTOR(1.d0,0.d0,0.d0)
     posvec = VECTOR(2.d0 * grid%octreeroot%subcellsize,0.d0, real(i) * grid%octreeroot%subcellsize * 0.01)
     viewvec = (-1.d0) * unitvec
     call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
     
     write(message,*) i, tau
     call writeinfo(message, FORINFO)
  enddo

call intensityAlongRay(VECTOR(0.d0,0.d0,0.d0),VECTOR(1d-20,1d-20,1.d0), grid, thisMolecule, 1, 0.d0, dummy, tau, .true.)

xmidplane = rinner + (router-rinner) * 0.5

if(writeoutput) write(*,*) "Midplane tau!"
open(50, file = 'dustcheck.dat', status="unknown", form = "formatted") 
currentposition(1) = VECTOR(xmidplane,0.d0,0.d0)

call findSubcellTD(currentPosition(1), grid%octreeroot, thisOctal, subcell)

do i = 0, nlamb
   
   lamb = (10.0**(real(i) * log10(lamend/lamstart) / 500.0)) * lamstart
   
   call continuumIntensityAlongRay(VECTOR(1.d-10,1d-10,1d-10),VECTOR(1.0,1d-20,1.d-20), grid, lamb, 0.d0, dummy, &
        tau, .true.)
   call locate(grid%lamArray, size(grid%lamArray), lamb, ilamb)
   call returnKappa(grid, thisOctal, subcell, ilambda = ilamb, lambda = lamb, kappaAbs = kappaAbs, kappaSca = kappaSca)
   write(50, '(f8.4,tr3,f10.4,tr3,f10.4,tr3f10.4)') lamb *1e-4, tau, kappaAbs*1e-10 / thisOctal%rho(subcell),&
        (kappaAbs+kappaSca)*1e-10 / thisOctal%rho(subcell)
enddo

close(50)

if(debug) then
   
   currentposition(1) =  VECTOR(xmidplane,0.d0,0.d0)
   currentposition(2) =  VECTOR(xmidplane,0.d0,xmidplane)
   currentposition(3) = VECTOR(xmidplane/sqrt(2.),xmidplane/sqrt(2.),2. * xmidplane)
   
   do i = 1,3
      
      call findSubcellTD(currentPosition(i), grid%octreeroot, thisOctal, subcell)
      currentposition(i) = subcellcentre(thisOctal, subcell)
      
      write(message,*) "currentposition", currentposition(i)
      call writeinfo(message, FORINFO)
      write (message,*) "r", sqrt(currentposition(i)%x**2+currentposition(i)%y**2)
      call writeinfo(message, FORINFO)
      write(message,*) "cell velocity",modulus(thisOctal%velocity(subcell)) * cspeed / 1d5
      call writeinfo(message, FORINFO)
      write(message,*) "calc V ",modulus(keplerianVelocity(currentposition(i),grid)) * cspeed / 1d5
      call writeinfo(message, FORINFO)
      
   enddo
   
endif

call continuumIntensityAlongRay(VECTOR(-1.d10,1.d-10,1.d-10),VECTOR(1.d0,-1d-10,-1d-10), grid, 1e4,&
     0.d0, dummy, tau, .true.)
write(message,*) "TAU @ 1micron", tau
call writeinfo(message, FORINFO)

end subroutine testOpticalDepth

subroutine plotdiscValues(grid, thisMolecule)

  type(GRIDTYPE) :: grid
  type(MOLECULETYPE) :: thisMolecule
  real(double) :: mean(6)
  
  call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule,0.d0,.true.)
    
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 2)
  
  call readAMRgrid("molecular_tmp.grid",.false.,grid)
  
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 0)
  
  call readAMRgrid("molecular_tmp.grid",.false.,grid)
  
  ! Set everything back to the way it was?
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 1)

end subroutine plotdiscValues

     function velocity(position, grid, startOctal, subcell) RESULT(out)

       implicit none
       type(VECTOR) :: out
       type(VECTOR), intent(in) :: position
       type(gridtype), intent(in) :: grid
       type(octal), pointer, optional :: startOctal
       integer, optional :: subcell

!    select case (grid%geometry)
       if(molebench)          Out = Molebenchvelocity(Position, grid)
!       case("molebench")

          
!       case("iras04158","shakara")
!          out = keplerianVelocity(position, grid)

!       case default
!          if(present(startOctal) .and. present(subcell)) then
!             out = amrGridVelocity(grid%octreeRoot, position, startOctal = startOctal, actualSubcell = subcell)
!          else
!             out = amrGridVelocity(grid%octreeRoot, position)
!          endif
!    end select

  end function velocity

  subroutine DecideOctalOrder(ioctal_beg, ioctal_end, iOctalArray, shouldinterlace)
    
    integer :: ioctal_beg, ioctal_end
    integer :: iOctalArray(ioctal_end - ioctal_beg + 1)
    logical, optional :: shouldinterlace
    logical :: interlace
    integer :: i
    
    if(present(shouldinterlace)) then
       interlace = shouldinterlace
    else
       interlace = .false.
    endif
    
    if(interlace) then
       
       do i = 1, (ioctal_end - iOctal_beg)/2 + 1
          iOctalArray(i) = 2 * (iOctal_beg/2) + 1 + 2 * (i - 1)
       enddo
       
       do i = (ioctal_end - iOctal_beg)/2 + 2, iOctal_end - iOctal_beg + 1
          iOctalArray(i) = 2 * ((iOctal_beg + 1)/2) + 2 * (i - ((ioctal_end - iOctal_beg)/2 + 2))
       enddo
       
    else
       
       do i = 1, size(iOctalArray)
          iOctalArray(i) = i
       enddo
       
    endif
  end subroutine DecideOctalOrder
  
  subroutine calculateConvergenceData(grid, nvoxels, fixedrays, maxRMSFracChange)
    
    use input_variables, only : tolerance

    integer :: nvoxels
    logical :: fixedrays

    type(GRIDTYPE) :: grid
    real(double) :: maxFracChangePerLevel(minlevel)
    integer :: ConvergenceCounter(4, minlevel) ! 3 different levels of convergence
    real(double) :: avgFracChange(minlevel - 1, 2) ! don't want the uppermost level to count for convergence
    real(double) :: maxavgFracChange, maxRMSFracChange, maxFracChange
    integer :: maxavgtrans(1), maxRMStrans(1)
    integer :: i

    character(len=160) :: message
    
10006 format(i6,tr3,6(f8.5,1x), 27x, 3x,l1,tr3,3(f5.3,tr3),f7.5,tr3,f7.5)
10007 format(i6,tr3,7(f8.5,1x), 18x, 3x,l1,tr3,3(f5.3,tr3),f7.5,tr3,f7.5)
10008 format(i6,tr3,8(f8.5,1x),  9x, 3x,l1,tr3,3(f5.3,tr3),f7.5,tr3,f7.5)
10009 format(i6,tr3,9(f8.5,1x),      3x,l1,tr3,3(f5.3,tr3),f7.5,tr3,f7.5)
10010 format(i6,tr3,9(f8.5,1x),      3x,l1,tr3,3(f5.3,tr3),f7.5,tr3,f7.5)
10011 format(i6,tr3,9(f8.5,1x),      3x,l1,tr3,3(f5.3,tr3),f7.5,tr3,f7.5)
10012 format(i6,tr3,9(f8.5,1x),      3x,l1,tr3,3(f5.3,tr3),f7.5,tr3,f7.5)
    
    maxFracChangePerLevel = -1.d30
    ConvergenceCounter    = 0
    avgFracChange = 0.d0

    call swapPops(grid%octreeRoot, maxFracChangePerLevel, avgFracChange, &
         convergenceCounter, grand_iter, nVoxels, fixedrays) ! compares level populations between this and previous levels 

    maxavgFracChange = maxval(avgFracChange(:,1))
    maxRMSFracChange = maxval(avgFracChange(:,2))

    maxavgtrans = maxloc(avgFracChange(:,1)) - 1 ! array index -> level
    maxRMStrans = maxloc(avgFracChange(:,2)) - 1

    maxFracChange = MAXVAL(maxFracChangePerLevel(1:mintrans)) ! Largest change of any level < 6 in any voxel
  
    write(message,'(a,1x,f11.7,1x,a,1x,f5.3,1x,a,2x,l1,1x,a,1x,i6,1x)') &
         "Maximum fractional change this iteration ", maxFracChange, "tolerance", tolerance, "fixed rays", fixedrays, &
         "nray", nray
    call writeInfo(message,FORINFO)

    write(message,'(a,f11.7,a,i1)') "  Average fractional change this iteration  ", &
         maxavgFracChange/real(nVoxels)," in level ",maxavgTrans
    call writeInfo(message,FORINFO)
    write(message,'(a,f11.7,a,i1)') "  RMS fractional change this iteration      ", &
         sqrt(maxRMSFracChange/real(nVoxels))," in level ",maxRMSTrans
        call writeInfo(message,FORINFO)
    write(message,'(a,f11.7)') "  Std Dev                                   ", &
         sqrt(maxRMSFracChange/real(nVoxels)-(maxavgFracChange/real(nVoxels))**2)
    call writeInfo(message,FORINFO)

    
    do i=3,1,-1
       if(i .lt. 3) write(message,'(a,f6.4,a, 14(f6.4,2x))') "Individual levels converged @ ",tolerance * i," | ",&
            real(convergenceCounter(i,1:minlevel))/real(nVoxels) ! Fraction of first 8 levels converged at 1,2,5 * tolerance respectively
       if(i .eq. 3) write(message,'(a,f6.4,a, 14(f6.4,2x))') "Individual levels converged @ ",tolerance * (i+2)," | ", &
            real(convergenceCounter(i,1:minlevel))/real(nVoxels) ! Fraction of first 8 levels converged at 1,2,5 * tolerance respectively
       call writeInfo(message,FORINFO)         
    enddo
  
    call writeinfo("",FORINFO)
    write(message,'(a,f6.4,a,12(f6.4,2x))') "Individual levels converged @ ",tolerance * 0.5," | ", &
         real(convergenceCounter(4,1:minlevel))/real(nVoxels) ! Fraction of first 8 levels converged at 1,2,5 * tolerance respectively
    call writeInfo(message,FORINFO)
  
    if (writeoutput) then
       open(138,file="fracChanges.dat",position="append",status="unknown")
  
    select case(minlevel)
    case(6)
       write(138,10006) nray, maxFracChangePerLevel(1:minlevel), fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels),&
            maxFracChange, maxavgFracChange/real(nVoxels)
    case(7)
       write(138,10007) nray, maxFracChangePerLevel(1:minlevel), fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels),&
            maxFracChange, maxavgFracChange/real(nVoxels)
    case(8)
       write(138,10008) nray, maxFracChangePerLevel(1:minlevel), fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels),&
            maxFracChange, maxavgFracChange/real(nVoxels)
    case(9)
       write(138,10009) nray, maxFracChangePerLevel(1:minlevel), fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels),&
            maxFracChange, maxavgFracChange/real(nVoxels)
    case(10)
       write(138,10010) nray, maxFracChangePerLevel(1:9), fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels),  &
            maxFracChange, maxavgFracChange/real(nVoxels)
    case(11)
       write(138,10011) nray, maxFracChangePerLevel(1:9), fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels),  &
            maxFracChange, maxavgFracChange/real(nVoxels)
    case(12)
       write(138,10012) nray, maxFracChangePerLevel(1:9), fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels),  &
            maxFracChange, maxavgFracChange/real(nVoxels)
    case default
       write(138,*) nray, maxFracChangePerLevel(1:9), fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels),  &
            maxFracChange, maxavgFracChange/real(nVoxels)
    end select
       
    close(138)
  
    open(139, file="avgChange.dat", position="append", status="unknown")
20  format(i2,tr3,i6,tr3,12(f7.5,1x))
    write(139,20) grand_iter, nray, avgFracChange(1:minlevel-1,1)/real(nvoxels)
    
    close(139)
    
    open(140, file="avgRMSchange.dat", position="append", status="unknown")
    write(140,20) grand_iter, nray, sqrt(avgFracChange(1:minlevel-1,2)/real(nvoxels))
       
    close(140)
 endif

end subroutine calculateConvergenceData
     
  real(double) function bNum(nu,T)
    
    real(double) :: fac1, fac2, fac3, nu
    real :: T
    real(double), parameter :: twoHOverCSquared = 1.4745d-43 
    
    fac1 = twoHOverCSquared * nu**3
    fac3 =  (hCgs*nu)/ (kErg * T) 
    if (fac3 > 100.d0) then
       fac2 = 0.d0
    else
       fac2 = 1.d0/(exp(fac3) - 1.d0)
    endif
    
    bNum = fac1 * fac2
  end function bNum
  
   recursive subroutine  solveAllPops(grid, thisOctal, thisMolecule)
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren ! What's this?
              if (thisOctal%indexChild(i) .eq. subcell) then
                 child => thisOctal%child(i)
                 call solveAllpops(grid, child, thisMolecule)
                 exit
              end if
           end do
        else ! once octal has no more children, solve for current parameter set
           call solveLevels(thisOctal%molecularLevel(subcell,1:maxlevel), &
                thisOctal%jnu(subcell,1:maxtrans),  &
                dble(thisOctal%temperature(subcell)), thisMolecule, thisOctal%nh2(subcell))
        endif
     enddo
   end subroutine solveAllPops

 subroutine LTEpops(thisMolecule, temperature, levelpops)

   type(MOLECULETYPE) :: thisMolecule
   real(double) :: temperature, nsum
   real(double) :: levelpops(:), fac
   integer :: i

   levelpops(1) = 1.d0
   nsum = 1.d0

   do i=2,size(levelpops)
      fac = abs(thisMolecule%energy(i)-thisMolecule%energy(i-1)) / (kev*temperature)
      levelpops(i) = max(1d-30,levelpops(i-1) * thisMolecule%g(i)/thisMolecule%g(i-1) * exp(-fac))
      nsum = nsum + levelpops(i)
   enddo

   levelpops = levelpops / nsum

 end subroutine LTEpops


end module molecular_mod
    
