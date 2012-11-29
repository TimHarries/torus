#ifdef MOLECULAR
module molecular_mod 
 ! written by tjh
 ! 21/11/06

 ! made my own by dar
 ! 01/X/07

   use kind_mod
   use constants_mod
   use utils_mod
   use messages_mod
   use timing, only: tune
   use amr_mod
   use random_mod
   use octal_mod
   use math_mod
   use datacube_mod, only: DATACUBE, TELESCOPE, initCube, npixels, &
        addspatialaxes, addvelocityaxis, convertspatialaxes
#ifdef USECFITSIO
   use datacube_mod, only : writeDataCube
#endif
   use parallel_mod, only : torus_mpi_barrier
   use gridio_mod, only: readamrgrid, writeamrgrid
   use atom_mod, only: bnu
   use vtk_mod, only: writeVtkFile
   use mpi_global_mod
#ifdef USEMKL
   use mkl95_lapack
   use mkl95_precision
#endif
#ifdef USEIEEEISNAN
   use ieee_arithmetic, isnan=> ieee_is_nan
#endif
   implicit none

   interface LTEpops
      module procedure LTEpopsReal
      module procedure LTEpopsDouble
   end interface

   integer, parameter :: maxray = 1048577
   integer :: mintrans, minlevel, maxlevel, maxtrans, nlevels, ntrans
   integer :: grand_iter, ngcounter
   integer :: nray
   Logical :: molebench, molcluster, chrisdisc, ggtau, hhobench,agbstar
   real(double), allocatable :: lamarray(:), lambda(:)
   integer :: ilambda

   integer(bigInt) :: iseedpublic
   real(double) :: r1(5) !! quasi random number generators
   real(double) :: r 
   integer :: accstep = 5
   integer :: accstepgrand = 4
   integer :: deptharray(50) = 0

   integer :: iupper(200), iLower(200)

   character(len=20) :: molgridfilename, molgridltefilename

   logical :: allCellsConverged, gridConverged, gridConvergedTest

   real(double) :: vexact, vlin, vquad, vquaddiff, vlindiff, vlindiffcounter
   real(double) :: vlindiffnormcounter, vquaddiffcounter, vquaddiffnormcounter = 0.d0


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
      real(double), pointer :: einsteinA(:)
      real(double), pointer :: einsteinBlu(:)
      real(double), pointer :: einsteinBul(:)
      real(double), pointer :: transfreq(:)
      integer, pointer :: itransUpper(:)
      integer, pointer :: itransLower(:)
      real(double), pointer :: Eu(:)
      integer, pointer :: iCollUpper(:,:)
      integer, pointer :: iCollLower(:,:)
      integer :: nCollPart
      character(len=10),pointer :: collPartnerName(:)
      integer, pointer :: iCollPartner(:)
      integer, pointer :: nCollTrans(:)
      integer, pointer :: nCollTemps(:)
      real(double), pointer :: collTemps(:,:)
      real(double), pointer :: collRates(:,:,:)
   end type MOLECULETYPE

   type(MOLECULETYPE) :: globalMolecule

 contains
   ! Read in molecular parameters from file - note: abundance hard-coded here

   subroutine readMolecule(thisMolecule, molFilename)
     use unix_mod, only: unixGetenv
     use inputs_mod, only : molAbundance
     type(MOLECULETYPE) :: thisMolecule
     character(len=*) :: molFilename
     character(len=80) :: junk
     character(len=200):: dataDirectory, filename
     integer :: i, j, iLow, iUp, iPart
     real(double) :: a, freq, eu, c(20)
     logical :: preprocess1 = .true.
     logical :: preprocess2 = .true.
     integer :: maxnCollTrans, maxnCollTemps

     if (trim(molFilename) == "hco_benchmark.mol") then
        call readBenchmarkMolecule(thisMolecule, molFilename)
        goto 666
     endif

     thisMolecule%abundance = molAbundance ! fixed at benchmark value here

       do while(preprocess2)
        if(.not. preprocess1) preprocess2 = .false.

        call unixGetenv("TORUS_DATA", dataDirectory)
        filename = trim(dataDirectory)//"/"//molfilename
        
        open(30, file=filename, status="old", form="formatted")

        read(30,*) junk
        read(30,'(a)') thisMolecule%molecule
        if( preprocess2) call writeInfo("Reading data for molecule: "//trim(thisMolecule%molecule),IMPORTANT)

        write(molgridfilename,*) trim(thismolecule%molecule),"_grid.grid"

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
                
        call allocateAttribute(thisMolecule%einsteinA, thisMolecule%nTrans)
        call allocateAttribute(thisMolecule%einsteinBul, thisMolecule%nTrans)
        call allocateAttribute(thisMolecule%einsteinBlu, thisMolecule%nTrans)
        call allocateAttribute(thisMolecule%transfreq, thisMolecule%nTrans)
        call allocateAttribute(thisMolecule%itransupper, thisMolecule%nTrans)
        call allocateAttribute(thisMolecule%itranslower, thisMolecule%nTrans)
        call allocateAttribute(thisMolecule%eu, thisMolecule%nTrans)

!        allocate(thisMolecule%einsteinBul(1:thisMolecule%nTrans))
!        allocate(thisMolecule%einsteinBlu(1:thisMolecule%nTrans))
!        allocate(thisMolecule%transfreq(1:thisMolecule%nTrans))
!        allocate(thisMolecule%itransUpper(1:thisMolecule%nTrans))
!        allocate(thisMolecule%itransLower(1:thisMolecule%nTrans))
!        allocate(thisMolecule%eu(1:thisMolecule%nTrans))

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
        allocate(thisMolecule%collPartnerName(1:thisMolecule%nCollPart))
        allocate(thisMolecule%iCollPartner(1:thisMolecule%nCollPart))
        allocate(thisMolecule%collRates(1:thisMolecule%nCollPart, maxnCollTrans, maxnCollTemps))
        allocate(thisMolecule%iCollUpper(thisMolecule%nCollPart, maxnCollTrans))
        allocate(thisMolecule%iCollLower(thisMolecule%nCollPart, maxnCollTrans))

        thisMolecule%collRates = 0.d0

        do iPart = 1, thisMolecule%nCollPart

           read(30,*) junk
           read(30,*) junk

           call parseCollisionPartner(junk, iPart, thisMolecule)

           read(30,*) junk
           read(30,*) thisMolecule%nCollTrans(iPart)

           read(30,*) junk
           read(30,*) thisMolecule%nCollTemps(iPart)

           read(30,*) junk
           read(30,*) thisMolecule%collTemps(iPart,1:thisMolecule%nCollTemps(ipart))

           read(30,*) junk

           do j = 1, thisMolecule%nCollTrans(iPart)

              read(30,*) i, thisMolecule%iCollUpper(ipart,j), &
                   thisMolecule%iCollLower(ipart,j), c(1:thisMolecule%nCollTemps(iPart))
              thisMolecule%collRates(iPart, j, 1:thisMolecule%nCollTemps(iPart)) = c(1:thisMolecule%nCollTemps(iPart))

           enddo
        enddo
        close(30)

        if(preprocess1) then
           deallocate(thisMolecule%energy, thisMolecule%g, thisMolecule%j, &
                      thisMolecule%einsteinA, thisMolecule%einsteinBlu, thisMolecule%einsteinBul, &
                      thisMolecule%transfreq, thisMolecule%itransUpper, thisMolecule%itransLower, &
                      thisMolecule%Eu, thisMolecule%iCollUpper, thisMolecule%iCollLower, &
                      thisMolecule%collPartnerName, thisMolecule%iCollPartner, &
                      thisMolecule%collTemps, thisMolecule%collRates)
        endif
        preprocess1 = .false.
     enddo

     call writeInfo("Done.", IMPORTANT)
666 continue
   end subroutine readMolecule

   subroutine parseCollisionPartner(cstring, iPart, thisMolecule)
     type(MOLECULETYPE) :: thisMolecule
     integer :: iPart
     character(len=*) :: cString
     
     read(cString,'(i1,a)') thisMolecule%iCollPartner(iPart),thisMolecule%collPartnerName(iPart)

   end subroutine parseCollisionPartner

 ! Same routine for benchmark molecule (slightly different file format)
   subroutine readBenchmarkMolecule(thisMolecule, molFilename)
     use unix_mod, only: unixGetenv
     use inputs_mod, only : molAbundance, zeroGhosts
     type(MOLECULETYPE) :: thisMolecule
     character(len=*) :: molFilename
     character(len=200):: dataDirectory, filename
     character(len=80) :: junk
     integer :: i, iPart

     thisMolecule%abundance = molAbundance ! fixed at benchmark value here

     call unixGetenv("TORUS_DATA", dataDirectory)
     filename = trim(dataDirectory)//"/"//molfilename
     call writeInfo("Opening file: "//trim(filename),TRIVIAL)
     open(30, file=filename, status="old", form="formatted")

     read(30,'(a)') thisMolecule%molecule
     call writeInfo("Reading data for molecule: "//trim(thisMolecule%molecule),IMPORTANT)

     write(molgridfilename,'(a)') "HCO+_grid.grid"

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
     allocate(thisMolecule%iCollPartner(1:1), thisMolecule%collPartnerName(1:1))
     thisMolecule%iCollPartner(iPart) = 1
     thisMolecule%collPartnerName(iPart) = "H2"

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
   recursive subroutine  allocateMolecularLevels(grid, thisOctal, thisMolecule)

     use grid_mod, only: freeGrid
     use inputs_mod, only : vturb, restart, isinLTE, &
          addnewmoldata, setmaxlevel, doCOchemistry, x_d, x_0, removeHotMolecular, &
          molAbundance, usedust, getdepartcoeffs, constantAbundance, photoionPhysics, zeroghosts

!         plotlevels
!     type(VECTOR) :: pos
!    REAL :: T1, T2, TSTART

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, ichild, isubcell
     integer :: i
     integer :: nOctal, nVoxels
     logical, save :: firsttime1 = .true.
     logical, save :: firsttime2 = .true.
     logical :: inlte
     logical, save :: firstAbundanceWarning = .true.
     logical, save :: photoionmessage=.true.
     real(double) :: nupper(200), nlower(200)
     real(double) :: levelpops(200), alphanubase(200), nmol

     !$OMP THREADPRIVATE( firstTime1, firstTime2, firstAbundanceWarning )

     inlte = isinlte

!Code to handle a restart. Read in grid from file 
     if(firsttime1 .and. restart) then
        call writeinfo("Reading in previous grid", FORINFO)
        call freeGrid(grid)
        call readAMRgrid(molgridfilename,.false.,grid)
        call writeinfo("Done", FORINFO)
        firsttime1 = .false.
     elseif(addnewmoldata .and. firsttime1) then ! read intermediate grid
!- grid is now read in elsewhere THAW
!        call readAMRgrid("notmolecular.grid",.false.,grid)
        restart = .false.
        firsttime1 = .false.
     endif

! Recursively fill all the subcells with molecular level data
! To save space when deailing with large numbers of octals determine maxlevels first.
     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do ichild= 1, thisOctal%nChildren
              if (thisOctal%indexChild(ichild) == subcell) then
                 child => thisOctal%child(ichild)
                 call allocateMolecularLevels(grid, child, thisMolecule)
                 exit
              end if
           end do
        else
! MPI stuff
           if (grid%splitOverMPI.and.(.not.octalOnThread(thisOctal, subcell, myrankGlobal))) cycle
! Allocate space for n_H2 and initialise, if not already done
           if (.not.associated(thisOctal%nh2)) then 
              allocate(thisOctal%nh2(1:thisOctal%maxChildren))
!- There was a bug here, it would only assign values to one subcell per octal THAW
!              thisOctal%nh2(subcell) = thisOctal%rho(subcell) / (2.d0*mHydrogen)
              thisOctal%nh2 = thisOctal%rho / (2.d0*mHydrogen)
           end if
! Temporary graph making code for molcluster. V1 code should ultimately be removed 
!           if(plotlevels .and. molcluster) then
!              if(firsttime2 .and. molcluster) then
!                 pos = clusterparameter(VECTOR(0.d0,0.d0,0.d0), thisoctal, subcell, isdone = .true.)
!                 call new_read_sph_data(sphdatafilename)
!                 call CPU_TIME(tstart)
!              endif
              
!              if(debug) then
!                 CALL CPU_TIME(T1)
!                 pos = subcellcentre(thisoctal, subcell)
!                 if(addnewmoldata) thisoctal%velocity(subcell) = molclustervelocity(pos, grid)
!                 thisoctal%linearvelocity(subcell) = &
!                      amrGridVelocity(grid%octreeRoot, pos, startOctal = thisOctal, &
!                      actualSubcell = subcell, linearinterp = .true.)
!                 thisoctal%quadvelocity(subcell) = &
!                      amrGridVelocity(grid%octreeRoot, pos, startOctal = thisOctal, &
!                      actualSubcell = subcell, linearinterp = .false.)
!                 counter = counter + 1
!                 CALL CPU_TIME(T2)

!                 write(111, *) counter, t2 - t1, t2 - tstart
!                 vexact = modulus(thisoctal%velocity(subcell)) * cspeed/1d5
                
!                 if(thisoctal%ndepth .gt. 3 .and. vexact .lt. 10.d0) then
!                    vlin = modulus(thisoctal%linearvelocity(subcell)) * cspeed/1d5
!                    vquad = modulus(thisoctal%quadvelocity(subcell)) * cspeed/1d5
!                    vquaddiff = modulus(thisoctal%velocity(subcell) - &
!                                        thisoctal%quadvelocity(subcell)) * cspeed/1d5
!                    vlindiff = modulus(thisoctal%velocity(subcell) - &
!                                       thisoctal%linearvelocity(subcell)) * cspeed/1d5
!                    vquaddiffcounter = vquaddiffcounter + vquaddiff
!                    vlindiffcounter = vlindiffcounter + vlindiff
!                    vquaddiffnormcounter = vquaddiffnormcounter + vquaddiff / vexact
!                    vlindiffnormcounter = vlindiffnormcounter + vlindiff / vexact
                 
!                    write(113,'(i2,tr2,3(f10.7,tr2),4(f10.3,tr2))') &
!                         thisoctal%ndepth, vexact, vquaddiff, vlindiff, &
!                         vquaddiffcounter, vlindiffcounter, vquaddiffnormcounter, vlindiffnormcounter 
!                 endif
!              endif
!           endif

! maxlevel - the greatest level that will be converged (but not necessarily counted for convergence)
           if(restart) then
! if restart then use previous maxlevel              
              if(firsttime2) then
                 if(setmaxlevel .eq. 0) then
                    maxlevel = size(thisOctal%molecularLevel(:,1))
                 else
                    maxlevel = setmaxlevel
                 endif
! TODO V2:
! maxtrans currently assumes linear molecule so maxlevel - 1 is J=maxlevel - maxlevel-1 (not nec the case)
                 maxtrans = maxlevel - 1
                 firsttime2 = .false.
              endif

              if(usedust) then
                 if (.not.associated(thisOctal%bnu)) then
                    allocate(thisOctal%bnu(1:maxlevel, 1:thisOctal%maxChildren))
                    
                    do isubcell = 1, thisoctal%maxchildren
                       do i = 1, maxtrans
                          thisOctal%bnu(i,isubcell) = bnu(thisMolecule%transFreq(i), &
                               dble(thisOctal%temperature(isubcell)))
                       enddo
                    enddo
                    
                 endif
              endif
! First 5 LTE level populations stored in departcoeff to measure departure from LTE
              if(getdepartcoeffs) then
                 if (.not.associated(thisOctal%departcoeff)) then
                    allocate(thisOctal%departcoeff(1:5,1:thisOctal%maxChildren))
                    
                    do isubcell = 1, thisoctal%maxchildren
                       call LTEpops(thisMolecule, dble(thisOctal%temperature(isubcell)), &
                       dble(thisOctal%departcoeff(1:5,isubcell)))
                       thisoctal%departcoeff(1:5,isubcell) = real(1.d0 / thisoctal%departcoeff(1:5,isubcell))
                    enddo
                 endif
              endif
! Implement CO drop model - See drundle thesis for details
              if(doCOchemistry) then
                 do isubcell = 1, thisoctal%maxchildren
                    if(thisOctal%nh2(isubcell) .gt. 3e4 .and. &
                       thisOctal%temperature(isubcell) .lt. 30.) &
                       thisOctal%molAbundance(isubcell) = x_D ! reduced fraction
                 enddo
              endif
 
           else ! from scratch, i.e. no restart part-converged grid
! Find maxlevel - the greatest level that will be converged (but not necessarily counted for convergence)
              if(firsttime2) then 
                 call countVoxels(grid%octreeRoot,nOctal,nVoxels)
                 if(setmaxlevel .eq. 0) then 
                    call findmaxlevel(grid, grid%octreeroot, thisMolecule, &
                         maxlevel, thismolecule%nlevels, nVoxels, lte = .true.)
                 else
                    maxlevel = setmaxlevel
                 endif
! TODO V2:
! maxtrans currently assumes linear molecule so maxlevel - 1 is J=maxlevel - maxlevel-1 (not nec the case)
                 maxtrans = maxlevel - 1
                 firsttime2 = .false.
              endif
! fill microturbulent velocity array with constant (in TORUS V1) 
! this should change to some function at some point
! set up microturbulence - this shouldn't be done here but is at the moment. It can easily be moved into a function.
! this just catches stuff if it's been allocated and not set.

              if (.not.associated(thisOctal%microturb)) then
                 allocate(thisOctal%microturb(1:thisOctal%maxChildren))
                 thisOctal%microTurb = 0.d0
              endif
              if(thisoctal%microturb(subcell) .le. 1d-20) then
!                 if(.not. molebench) 
                 thisOctal%microturb(subcell) = max(1d-7, &
                      sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / &
                      (thisMolecule%molecularWeight * amu)) + vturb**2 ) / (cspeed * 1d-5))
! 1d-10 is conversion from kerg -> k km^2.g.s^-2.K^-1 (10^-7 (erg->J) * 10^-6 (m^2-km^2) * 10^3 (kg->g))
! *2 because using 1/e definition of thermal line width
! molebench has vturb already
              endif

! Fill cells with molecular abundance data
              if (.not.associated(thisOctal%molAbundance)) &
                   allocate(thisOctal%molAbundance(1:thisOctal%maxChildren))

              if (constantAbundance) then
                 thisOctal%molAbundance(subcell) = molAbundance
              elseif(doCOchemistry) then
                 if(thisOctal%nh2(subcell) .gt. 3e4 .and. thisOctal%temperature(subcell) .lt. 30.) then
                    thisOctal%molAbundance(subcell) = x_D ! depleted abundance
                 else
                    thisOctal%molAbundance(subcell) = x_0 ! normal abundance
                 endif
              else
                 if (firstAbundanceWarning) then
                    call writeinfo("constantAbundance off and doCOchemistry off")
                    call writeinfo("Implement a new molecular abundance elsewhere")
                    call writeinfo("because molabundance hasn't been changed here at all")
                    call writeinfo("... which might of course be what you wanted :)")
                    firstAbundanceWarning = .false.
                 endif
              endif
              if (photoionPhysics.or.removeHotMolecular .and. photoionmessage) then
                 call writeInfo("Applying photoion chemistry")
                 call photoionChemistry(grid, thisOctal, subcell)
                 photoionmessage = .false.
              endif


#ifdef MPI
#ifdef HYDRO
              if(zeroghosts) then
                 if(thisOctal%ghostcell(subcell)) then
                    thisOctal%molAbundance(subcell) = 1.d-30
                 end if
              end if
#endif
#endif

! Fill cells with LTE data to determine departure coefficents if required
              if(getdepartcoeffs) then
                 if (.not.associated(thisOctal%departcoeff)) &
                      allocate(thisOctal%departcoeff(1:5,1:thisOctal%maxChildren))
                                 
                 call LTEpops(thisMolecule, dble(thisOctal%temperature(subcell)), levelpops(1:5))
                 do i = 1, 5
                    thisoctal%departcoeff(i,subcell) = real(max(levelpops(i),1d-30))
                 enddo
                 thisoctal%departcoeff(1:5,subcell) = real(1.d0 / thisoctal%departcoeff(1:5,subcell))
              endif
! Fill cells with molecular level populations - LTE or small
              if (.not. associated(thisOctal%molecularLevel)) &
                   allocate(thisOctal%molecularLevel(1:maxlevel,1:thisOctal%maxChildren))

                 if((grid%geometry .eq. "h2obench1") .or. (grid%geometry .eq. "h2obench2")) then
                    thisOctal%molecularLevel(1,subcell) = 1.d0-1d-10
                    thisOctal%molecularLevel(2,subcell) = 1d-10
                 else
                    if(inlte) then
                       call LTEpops(thisMolecule, dble(thisOctal%temperature(subcell)), &
                            levelpops(1:maxlevel))
                       thisOctal%molecularLevel(1:maxlevel,subcell) = levelpops(1:maxlevel)
                    elseif(maxlevel .gt. 3) then
                       thisOctal%molecularLevel(1:2,subcell) = 0.5d0               
                       thisOctal%molecularLevel(3:maxlevel,subcell) = 1.d-10
                    else
                       thisOctal%molecularLevel(1:maxlevel,subcell) = 1.d-10
                    endif
                 endif

              if (.not.associated(thisOctal%bnu)) &
                 allocate(thisOctal%bnu(1:maxtrans, thisOctal%maxChildren))
              do i = 1, maxtrans
                 thisOctal%bnu(i,subcell) = bnu(thisMolecule%transFreq(i), dble(thisOctal%temperature(subcell)))
              enddo

! Get jnu from this cell
              if (.not.associated(thisOctal%jnu)) &
                 allocate(thisOctal%jnu(1:maxtrans,1:thisOctal%maxChildren))
              
              nmol = thisoctal%nh2(subcell) * thisoctal%molabundance(subcell)

              nlower(1:maxtrans) = thisOctal%molecularLevel(iLower(1:maxtrans),subcell) * nMol
              nupper(1:maxtrans) = thisOctal%molecularLevel(iUpper(1:maxtrans),subcell) * nMol
              alphanuBase(1:maxtrans) = nLower(1:maxtrans) * thisMolecule%einsteinBlu(1:maxtrans) - &
                   nUpper(1:maxtrans) * thisMolecule%einsteinBul(1:maxtrans)

              thisOctal%jnu(1:maxtrans, subcell) = &
                   thisoctal%bnu(1:maxtrans,subcell) * alphanuBase(1:maxtrans)
              
! deallocate bnu if not using dust. It takes up a lot of space.
              if((.not. usedust).and.(.not.isinLTE)) deallocate(thisoctal%bnu)
              
           endif ! if(restart) - these things are common to all allocations

! If gas and dust temperature are different then deallocate dust and gas temperature and use common temp
           if(associated(thisoctal%temperaturedust) .and. associated(thisoctal%temperaturegas)) then
              if(thisOctal%temperaturedust(subcell) .eq. thisOctal%temperaturegas(subcell)) then
                 deallocate(thisoctal%temperaturedust)
                 deallocate(thisoctal%temperaturegas)
              endif
           endif

! molmicroturb = 1/microturb which is used far more commonly. For speed
           if (.not.associated(thisOctal%molmicroturb)) then
              allocate(thisOctal%molmicroturb(1:thisOctal%maxChildren))
           endif
           thisOctal%molmicroturb(subcell) = 1.d0 / thisOctal%microturb(subcell)

           deptharray(thisoctal%ndepth + 1) = deptharray(thisoctal%ndepth + 1) + 1
        endif ! if haschild
        
     enddo ! do over all children

   end subroutine allocateMolecularLevels

   recursive subroutine  allocateOther(grid, thisOctal)

     use inputs_mod, only : gettau

     type(GRIDTYPE) :: grid
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, ichild
     real(double) :: temparray(100,8)
     integer :: temparray_int(100,8)
     integer :: h

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do ichild= 1, thisOctal%nChildren, 1
              if (thisOctal%indexChild(ichild) == subcell) then
                 child => thisOctal%child(ichild)
                 call allocateOther(grid, child)
                 exit
              end if
           end do
        else



! new/old/oldest are all needed for ng acceleration

           if (grid%splitOverMPI.and.(.not.octalOnThread(thisOctal, subcell, myrankGlobal))) cycle

           if(associated(thisOctal%newmolecularlevel)) then
              temparray(1:maxlevel,1:thisoctal%maxchildren) = &
                   thisoctal%newmolecularlevel(1:maxlevel,1:thisoctal%maxchildren)
              deallocate(thisoctal%newmolecularlevel)
              allocate(thisoctal%newmolecularlevel(maxlevel,1:thisoctal%maxchildren))
              thisoctal%newmolecularlevel(1:maxlevel,1:thisoctal%maxchildren) = &
                   temparray(1:maxlevel,1:thisoctal%maxchildren)
           else
              allocate(thisoctal%newmolecularlevel(maxlevel,1:thisoctal%maxchildren))
              thisoctal%newmolecularlevel(1:maxlevel,1:thisOctal%maxChildren) = &
                   thisoctal%molecularlevel(1:maxlevel,1:thisOctal%maxChildren)
           endif

           if(associated(thisOctal%oldmolecularlevel)) then
              h = size(thisoctal%oldmolecularlevel(:,1))
              temparray(1:h,1:thisoctal%maxchildren) = thisoctal%oldmolecularlevel(1:h,1:thisoctal%maxchildren)
              deallocate(thisoctal%oldmolecularlevel)
              allocate(thisoctal%oldmolecularlevel(1:minlevel,1:thisoctal%maxchildren))
              thisoctal%oldmolecularlevel(1:minlevel,1:thisoctal%maxchildren) = &
                   temparray(1:minlevel,1:thisoctal%maxchildren)
           else
              allocate(thisoctal%oldmolecularlevel(minlevel,1:thisoctal%maxchildren))
              thisoctal%oldmolecularlevel(1:minlevel,1:thisOctal%maxChildren) = &
                   thisoctal%molecularlevel(1:minlevel,1:thisOctal%maxChildren)
           endif

           if(associated(thisOctal%oldestmolecularlevel)) then
              h = size(thisoctal%oldestmolecularlevel(:,1))
              temparray(1:h,1:thisoctal%maxchildren) = thisoctal%oldestmolecularlevel(1:h,1:thisoctal%maxchildren)
              deallocate(thisoctal%oldestmolecularlevel)
              allocate(thisoctal%oldestmolecularlevel(1:minlevel,1:thisoctal%maxchildren))
              thisoctal%oldestmolecularlevel(1:minlevel,1:thisoctal%maxchildren) = &
                   temparray(1:minlevel,1:thisoctal%maxchildren)
           else
              allocate(thisoctal%oldestmolecularlevel(minlevel,1:thisoctal%maxchildren))
              thisoctal%oldestmolecularlevel(1:minlevel,1:thisOctal%maxChildren) = &
                   thisoctal%molecularlevel(1:minlevel,1:thisOctal%maxChildren)
           endif

           if(associated(thisOctal%levelconvergence)) then
              h = size(thisoctal%levelconvergence(:,1))
              temparray_int(1:h,1:thisoctal%maxchildren) = thisoctal%levelconvergence(1:h,1:thisoctal%maxchildren)
              deallocate(thisoctal%levelconvergence)
              allocate(thisoctal%levelconvergence(minlevel,1:thisoctal%maxchildren))
              thisoctal%levelconvergence(1:minlevel,1:thisoctal%maxchildren) = &
                   temparray_int(1:minlevel,1:thisoctal%maxchildren)
           else
              allocate(thisoctal%levelconvergence(minlevel,1:thisoctal%maxchildren))
           endif

           if(gettau) then
              if(associated(thisOctal%tau)) then
                 h = size(thisoctal%tau(:,1))
                 temparray(1:h,1:thisoctal%maxchildren) = thisoctal%tau(1:h,1:thisoctal%maxchildren)
                 deallocate(thisoctal%tau)
                 allocate(thisoctal%tau(minlevel,1:thisoctal%maxchildren))
                 thisoctal%tau(1:minlevel,1:thisoctal%maxchildren) = temparray(1:minlevel,1:thisoctal%maxchildren)
              else
                 allocate(thisoctal%tau(minlevel,1:thisoctal%maxchildren))
              endif
           endif

           if(associated(thisOctal%nsplit)) then
              continue
           else
              allocate(thisoctal%nsplit(1:thisoctal%maxchildren))
              thisoctal%nsplit(subcell) = 2
           endif

           if(associated(thisOctal%convergence)) then
              continue
           else
              allocate(thisoctal%convergence(1:thisoctal%maxchildren))
           endif
                     
        endif

     enddo
     
   end subroutine allocateOther

   recursive subroutine  deallocateUnused(grid, thisOctal, upperlev, everything)

     use inputs_mod, only : densitysubsample, debug, lowmemory
 
     type(GRIDTYPE) :: grid
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, ichild
     real(double) :: temparray(100,8)
     integer, optional :: upperlev
     integer :: upperlevel
     logical, optional :: everything
     logical :: deallocateeverything


     if(present(everything)) then
        deallocateeverything = everything
     else
        deallocateeverything = .false.
     endif

     if(present(upperlev)) then
        upperlevel = upperlev
     else
        upperlevel = 0
     endif

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do ichild= 1, thisOctal%nChildren, 1
              if (thisOctal%indexChild(ichild) == subcell) then
                 child => thisOctal%child(ichild)
                 call allocateOther(grid, child)
                 exit
              end if
           end do
        else

           if (grid%splitOverMPI.and.(.not.octalOnThread(thisOctal, subcell, myrankGlobal))) cycle

           if(deallocateeverything) then
              if(.not. (densitysubsample .and. associated(thisoctal%cornerrho))) deallocate(thisoctal%cornerrho)
              if(associated(thisoctal%microturb)) deallocate(thisoctal%microturb)
              if(.not. debug) deallocate(thisoctal%newmolecularlevel)

              if(.not. lowmemory) then
                 deallocate(thisoctal%molecularlevel)
              endif
           else
              if(lowmemory) then
                 temparray(1:upperlevel,1:thisoctal%maxchildren) = &
                      thisoctal%molecularlevel(1:upperlevel,1:thisoctal%maxchildren)
                 deallocate(thisoctal%molecularlevel)
                 allocate(thisoctal%molecularlevel(1:upperlevel,1:thisoctal%maxchildren))
                 thisoctal%molecularlevel(1:upperlevel,1:thisoctal%maxchildren) = &
                      temparray(1:upperlevel,1:thisoctal%maxchildren)
              endif
              if(debug) then
                 temparray(1:5,1:thisoctal%maxchildren) = thisoctal%newmolecularlevel(1:5,1:thisoctal%maxchildren)
                 deallocate(thisoctal%newmolecularlevel)
                 allocate(thisoctal%newmolecularlevel(1:5,1:thisoctal%maxchildren))
                 thisoctal%newmolecularlevel(1:5,1:thisoctal%maxchildren) = &
                      temparray(1:upperlevel,1:thisoctal%maxchildren)
              endif
           endif
        endif
     enddo
   end subroutine deallocateUnused

 ! Does a lot of work - do more rays whilst problem not converged -            
   subroutine molecularLoop(grid, thisMolecule)

     use inputs_mod, only : blockhandout, tolerance, &
          usedust, amr1d, amr3d, plotlevels,  &
          debug, restart, isinlte, quasi, dongstep, initnray, outputconvergence, dotune, &
          zeroGhosts, forceIniRay
     use messages_mod, only : myRankIsZero
     use dust_mod
     use parallel_mod
     use vtk_mod
#ifdef PHOTOION
     use photoion_utils_mod, only : quicksublimate
#endif
#ifdef MPI
     use mpi
#ifdef HYDRO
     use hydrodynamics_mod, only : setupedges, setupghosts
#endif
#endif

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(VECTOR) :: position, direction
     integer :: nOctal, iOctal, subcell
     character(len=80) :: mpiFilename
     type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
     type(OCTAL), pointer :: thisOctal
     integer, parameter :: maxIter = 50
     logical :: popsConverged, renewinputrays
     character(len=200) :: message
     integer :: iRay, iter, i 
     integer :: iStage
     real(double), allocatable :: oldpops1(:), oldpops2(:), oldpops3(:), oldpops4(:)
     real(double) :: fac
     real(double), allocatable :: ds(:), phi(:), i0(:,:), i0temp(:), dirWeight(:)

#ifdef MPI
     ! For MPI implementations
     integer       ::   my_rank        ! my processor rank
     integer       ::   np             ! The number of processes
     integer       ::   m
     integer       ::   n_rmdr
     integer       ::   ierr           ! error flag
     integer       ::   j
     real(double), allocatable :: tArrayd(:,:),tempArrayd(:,:) 
     real(double), allocatable :: jArrayd(:),tempjArrayd(:) 
#endif

     integer :: nVoxels
     integer :: ioctal_beg, ioctal_end
     logical :: fixedRays
     real(double) :: maxRMSfracChange
     
     character(len=30) :: filename
     
     real(double) :: tauarray(60) = -1.d0    
     logical :: warn = .true.
     integer :: warncount, warncount_all
     real(double) :: error(50)
     real(double) :: diffmax
     integer :: maxerrorloc, maxlocerror(1)
     integer :: mintransold
     
     real(double) :: nh, nHe, ne, nProtons
     integer :: dummy, ijunk
     logical :: ljunk
     logical :: juststarted = .true.
     logical :: ng
     integer :: status
     integer(bigint) :: fixedRaySeed

     real(double) :: collmatrix(50,50), ctot(50)

     logical :: warned_neg_dtau

     call writeinfo("molecular_mod 20100428.1400",TRIVIAL)

#ifdef PHOTOION
     if(usedust) then
        call quickSublimate(grid%octreeRoot, 0.01) ! do dust sublimation  
     end if
#endif

     write(mpiFilename,'(a, i4.4, a)') "quickDump.vtk"                                                                          
     call writeVtkFile(grid, mpiFilename, &                                                                                     
          valueTypeString=(/"rho          ", "dust1        " , "temperature  "/))

!     print *, "WRITTEN FILE"                         

     
! logicals are quicker to access than strings 
     if(grid%geometry .eq. 'molebench') molebench = .true.
     if(grid%geometry .eq. 'molcluster') molcluster = .true.
     if(grid%geometry .eq. 'iras04158') chrisdisc = .true.
     if(grid%geometry .eq. 'ggtau') ggtau = .true.
     if(grid%geometry .eq. 'agbstar') agbstar = .true.
     if(grid%geometry .eq. 'h2obench2') hhobench = .true.
!     usedust = .false.
     debug=.false.

! get pairs of radiative transitions stored in thismolecule. Re-assign for readability     
     iUpper(1:thismolecule%ntrans) = thisMolecule%iTransUpper(1:thismolecule%ntrans)
     iLower(1:thismolecule%ntrans) = thisMolecule%iTransLower(1:thismolecule%ntrans)

     nlevels = thisMolecule%nlevels     

! Let user know that Ng acceleration is turned on. inputs_mod doesn't say because it's optional
     ng = dongstep
     write(message,*) "Use Ng Acceleration: " , ng
     call writeinfo(message,TRIVIAL)

! blockhandout must be off for fixed ray case, otherwise setting the
! seed is not enough to ensure the same directions are done for
! each cell every iteration
      
     blockHandout = .false. 

#ifdef MPI
     ! FOR MPI IMPLEMENTATION=======================================================
     !  Get my process rank # 
     call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

     ! Find the total # of precessor being used in this run
     call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
#endif

! Write Column headings for file containing level populations and convergence data
     grand_iter = 0
     
     if(writeoutput .and. outputconvergence) then
18      format(a,tr3,8(a,tr8),a,6(tr2,a))
        open(138,file="fracChanges.dat",status="unknown",form="formatted")
        write(138,18) "nRays","J=0","J=1","J=2","J=3","J=4","J=5","J=6","J=7",&
                      "Fixed","5* Conv","2* Conv","1* Conv","Max","Avg"
        close(138)

        open(139,file="avgChange.dat",status="replace",form="formatted")
        open(140,file="avgRMSChange.dat",status="replace",form="formatted")
        open(141,file="tau.dat",status="replace",form="formatted")
        open(142,file="criticaldensities.dat",status="replace",form="formatted")
        open(999,file="tempcheck.dat",status="replace",form="formatted")
        close(999)
        close(139)
        close(140)
        close(141)
     endif

!! TORUS V1 code - if interfacing with lucy_mod, read in lucy temp file and start from here.
!! This will ultimately be removed     
!     if(openlucy) then
!        if(writeoutput) then
!           write(message,*) "Reading in lucy temp files: ",lucyfilenamein
!           call writeinfo(message,FORINFO)
!        endif
        
!        call freeGrid(grid)
!        call readAmrGrid(lucyfilenamein,.false.,grid)
        
!        call writeinfo("Successfully read in Lucy Grid file", FORINFO)
!        call writeinfo("Plotting Temperatures...",TRIVIAL)
!     endif

#ifdef MPI
#ifdef HYDRO
     if(zeroGhosts) then
        call setupedges(grid%octreeRoot, grid)
        call setupghosts(grid%octreeRoot, grid)
     end if
#endif
#endif


! IMPORTANT - Here's where molecular levels and other critical data get allocated.     
     call writeinfo("Allocating and initialising molecular levels", FORINFO)
     call allocateMolecularLevels(grid, grid%octreeRoot, thisMolecule)

! Useful temporary file containing the depth of each array - filled in allocate molecularlevels
     if(debug) then
        do i = 1, 50
           write(1003,*) i, deptharray(i)
        enddo
     endif

! Count number of subcells. Fill octalarray with pointers to all octals
     allocate(octalArray(grid%nOctals))
     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     dummy = 0
     call getOctalArray(grid%octreeRoot,octalArray, dummy)

! Write grid  with LTE populations if initialised to LTE     
     if(isinlte .and. .not. restart) then
        write(molgridltefilename,*) trim(thismolecule%molecule),"_lte.grid"

        call writeAMRgrid(molgridltefilename,.false.,grid)
!       goto 666
     endif

! Write maximum interesting level as determined by molecularlevel
     write(message, *) "Maximum Interesting Level", maxlevel
     call writeinfo(message, TRIVIAL)
! minlevel used for determining convergence and allocating less important variables (allocateother)
      minlevel = min(10, maxlevel-2)
      mintrans = minlevel - 1
! allocateother allocates new/old/oldest molecular levels for ng acceleration
! also tau and levelconvergence and convergence etc.
      call allocateother(grid, grid%octreeroot)
! set up lambda array if using dust
      if(usedust) then
         call allocateMemoryForDust(grid%octreeRoot)
         allocate(lamarray(size(grid%lamarray)))
         allocate(lambda(maxtrans))
         lamarray = grid%lamarray
         lambda(1:maxtrans) = cspeed / thisMolecule%transfreq(1:maxtrans) * 1d8
         call locate(lamArray, size(lamArray), lambda(maxtrans/2), ilambda)
      endif

! dump results if in LTE
      nRay = 1
      if(writeoutput .and. (.not. restart) .and. isinlte .and. outputconvergence) then
         call writeinfo("Writing LTE levels", TRIVIAL)
         call dumpresults(grid, thisMolecule)!, convtestarray) ! find radial pops on final grid     
      endif
! set number of rays used to estimate jnu and determine level pops 
      nRay = initnray
! initialise random seed for AMC method
! this line could be altered to give a fixed seed if re-producibility were required.
      call randomNumberGenerator(randomSeed=.true.)
      if(.not. restart) then 
         call randomNumberGenerator(randomSeed=.true.)
      else ! get previous random seed if restarting
         open(95, file="restart.dat",status="unknown",form="formatted")
         read(95,*) ljunk, ijunk, ijunk
         read(95,*) iseedpublic
         close(95)
      endif

      call randomNumberGenerator(getIseed = fixedRaySeed)

! This is the loop that controls everything
! 1) istage = 1, fixed rays to reduce variance between cells or 
! 2) istage = 2, random rays to ensure sufficient spatial/frequency sampling
      do iStage = 1, 2
         
         if (iStage .eq. 1) then
            fixedRays = .true.
            nRay = initnray
         else
            fixedRays = .false.
            nRay = max(initnray, nray)
! blockhandout will allow better load balancing if octals are divided unevenly
            blockhandout = .true.
         endif
! If restarting, get fixedrays, number of rays and iteration from previous run
         if(restart .and. juststarted) then 
            open(95, file="restart.dat",status="unknown",form="formatted")
            read(95, *) fixedrays, nray, grand_iter
            close(95)
            ngcounter = 0
            juststarted = .false.

            renewInputRays = .true.
            if(renewInputRays) then
               nRay = 1000
            end if
            if(fixedrays) cycle
         endif

         gridConvergedTest = .false.
         gridConverged = .false.


! while grid not converged iteratively find level populations
         do while (.not. gridConverged)
! new iteration
            grand_iter = grand_iter + 1
! controls Ng Acceleration, pronounced do-ng-step, not dongstep
            if(ng) ngcounter = ngcounter + 1
            mintransold = mintrans
! determine minlevel/mintrans
            if(.not. amr3d) then
               call taualongray(VECTOR(10.d0,10.d0,10.d0), VECTOR(0.577350269, 0.577350269, 0.577350269), &
                    grid, thisMolecule, 0.d0, tauarray(1:maxtrans))

               do i = size(tauarray), 1, -1
                  minlevel = min(i + 2,maxlevel - 2)
                  mintrans = minlevel - 1

!                  if(writeoutput .and. debug) write(99,*) i, mintrans, tauarray(i)
                  if(tauarray(i) .gt. 0.01) exit
               enddo

            else !minlevel = maxlevel-2 (upto maxlevel = 16)then  = maxlevel/2 
               minlevel = min(max(7,Maxlevel / 2),maxlevel - 2)
               mintrans = minlevel - 1
            endif

            minLevel = 6 !TJH !!!!!!!!!!!!!!!

            if(grid%geometry .eq. 'agbstar' .or. &
                 grid%geometry .eq. 'h2obench1' .or. &
                 grid%geometry .eq. 'h2obench2') then
               minlevel = maxlevel
               mintrans = maxtrans
            endif
! if mintrans has to change then need to change allocations for grid
            if((mintransold .ne. mintrans) .or. (grand_iter .eq. 0)) then
               call writeinfo("(Re)allocating memory", TRIVIAL)
               call allocateother(grid, grid%octreeroot)
               ngcounter = 0
            endif
            
            write(message,'(a,i3)') "Iteration ",grand_iter
            call writeinfo(message, FORINFO)

            write(message, '(a,i2)') "Minimum important level ", minlevel
            call writeinfo(message, TRIVIAL)
! Time from now until all cells level pops recalculated
            if(writeoutput .and. dotune) then
               write(message,*) "Done ",nray," rays"
               call tune(6, message)  ! start a stopwatch
            endif

! Use random seed determined earlier to generate random numbers or
! Initialise quasi-random number generator
            if (fixedRays) then
               call randomNumberGenerator(putIseed=fixedRaySeed)
               if(quasi) call sobseq(r1, -1)
            else
               call randomNumberGenerator(randomSeed=.true.)
            endif

!            call torus_mpi_barrier
!            do i = 0, nThreadsGlobal-1
!            call torus_mpi_barrier
!               if (i==myrankGlobal) then
!                  call randomNumberGenerator(getDouble=r)
!                  write(*,*) myrankGlobal,r
!               endif
!            call torus_mpi_barrier
!            enddo


! default loop indicies for single processor otherwise do MPI stuff
            ioctal_beg = 1
            ioctal_end = SIZE(octalArray)         

            warn = .true.
            warncount = 0

#ifdef MPI

! we will use an array to store the rank of the process
! which will calculate each octal's variables

            ! Set the range of index for octal loop used later.     
            np = nThreadsGlobal
            n_rmdr = MOD(SIZE(octalArray),np)
            m = SIZE(octalArray)/np
            
            if (myRankGlobal .lt. n_rmdr ) then
               ioctal_beg = (m+1)*myRankGlobal + 1
               ioctal_end = ioctal_beg + m
            else
               ioctal_beg = m*myRankGlobal + 1 + n_rmdr
               ioctal_end = ioctal_beg + m - 1
            end if
            

#endif

! Make sure all MPI processes and OpenMP threads have the correct random seed
            if (fixedRays .and. TorusOpenmp) then
               call randomNumberGenerator(putIseed=fixedRaySeed)
               call randomNumberGenerator(syncIseed=.true.)
               call test_same_hybrid
            endif

! iterate over all octals, all rays, solving the system self-consistently

                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP PRIVATE(iOctal, thisOctal, iray, warned_neg_dtau, dirweight) &
		!$OMP PRIVATE(nh, ne, nhe, nprotons, collMatrix, subcell) &
		!$OMP PRIVATE(ctot, ds, phi, direction, i0temp, i0, position, iter, popsconverged) &
		!$OMP PRIVATE(oldpops1, oldpops2, oldpops3, oldpops4, error, maxlocerror, maxerrorloc, fac) &
		!$OMP SHARED(grid, iOctal_beg, iOctal_end, octalArray, maxlevel, thisMolecule, nray,maxtrans) &
		!$OMP SHARED(fixedRays, debug, ng, minlevel, mintrans,  Warncount, accstep, fixedRaySeed, myRankGlobal) &
                !$OMP SHARED(forceIniRay)
            

! Allocate the main working arrays for i0, ds and phi
            allocate(ds(1:maxray))
            allocate(phi(1:maxray))
            allocate(dirWeight(1:maxray))
            allocate(i0temp(1:maxtrans))
            allocate(i0(1:maxray,1:maxtrans))
! set-up temporary arrays for ngstep
            allocate(oldPops1(1:maxlevel), oldPops2(1:maxlevel), oldPops3(1:maxlevel), oldPops4(1:maxlevel))

            warned_neg_dtau = .false. 

!#ifdef MPI
!            if(myrankglobal /= 0) then
!#endif

            !$OMP DO SCHEDULE(static)

    do iOctal = ioctal_beg, ioctal_end
!       if (debug .and. writeoutput) then
!          write(message,*) iOctal,ioctal_beg,ioctal_end
!          call writeInfo(message,TRIVIAL)
!       endif
! point thisoctal at the corresponding grid octal
       thisOctal => octalArray(ioctal)%content
! over all subcells in this octal
       do subcell = 1, thisOctal%maxChildren
! at the tip of the branch find the level populations in all subcells...
          if (.not.thisOctal%hasChild(subcell)) then
! Actually, first thing is to get collision matrix because this doesn't require knowing npops
             nh = 0.d0
             nHe = 0.d0
             ne = 0.d0
             nProtons = 0.d0

             call getCollMatrix(thisOctal%nh2(subcell), dble(thisOctal%temperature(subcell)), &
                  thismolecule, nh, nHe, ne, nprotons, collmatrix(1:maxlevel,1:maxlevel), &
                  ctot(1:maxlevel))

! First populate ds, phi and i0 so that they can be passed on to calculatejbar             

! thaw - the first ray is optionally forced towards a star
! this is done outside of the main getray loop to reduce slow down.
!             print *, "getting forced ray"
             if(forceIniRay) then
                call getRay(grid, thisOctal, subcell, position, direction, &
                     ds(1), phi(1), i0temp(1:maxtrans), &
                     thisMolecule,dirweight(1), fixedrays, warned_neg_dtau,tostar=.true.) ! does the hard work - populates i0 etc
                i0(1,1:maxtrans) = i0temp(1:maxtrans)
             else
                
              call getRay(grid, thisOctal, subcell, position, direction, &
                     ds(1), phi(1), i0temp(1:maxtrans), &
                     thisMolecule,dirweight(1), fixedrays, warned_neg_dtau, tostar=.false.) ! does the hard work - populates i0 etc
                i0(1,1:maxtrans) = i0temp(1:maxtrans)
             end if
!            print *, "getting main rays"

! main getray loop
             do iRay = 2, nRay
! does the hard work - populates i0 etc
                call getRay(grid, thisOctal, subcell, position, direction, &
                     ds(iRay), phi(iRay), i0temp(1:maxtrans), &
                     thisMolecule,dirweight(iRay), fixedrays, warned_neg_dtau, tostar=.false.) 
                i0(iray,1:maxtrans) = i0temp(1:maxtrans)
             enddo

!             print *, "got main rays"

             if(debug) where(isnan(i0)) i0 = 0.d0
! set iteration within subcell to 0
             iter = 0
             popsConverged = .false.
! Iterate between calculate jbar and solvelevels until converged
             do while (.not. popsConverged)
                iter = iter + 1
! update levels so that current become old, old -> older etc.
! If not using ng then just store current levels for comparison
                if(ng) then
                   oldpops1(1:minlevel) = oldpops2(1:minlevel)
                   oldpops2(1:minlevel) = oldpops3(1:minlevel)
                   oldpops3(1:minlevel) = thisOctal%newmolecularLevel(1:minlevel,subcell)
                else
                   oldpops3(1:minlevel) = thisOctal%newmolecularLevel(1:minlevel,subcell)
                endif
! calculate the average radiation field, jnu in this cell given ds, phi and i0.
! This feeds into and affects solvelevels which in turns affect calculatejbar

                call calculateJbar(nray, grid, thisOctal, subcell, thisMolecule, ds(1:nRay), & 
                     phi(1:nRay), i0(1:nray,1:maxtrans), thisOctal%newMolecularLevel(1:maxlevel,subcell), &
                     thisOctal%jnu(1:maxtrans,subcell), dirweight) ! calculate updated Jbar

!                write(*,*) "a",thisOctal%newMolecularLevel(1:3,subcell)
                if(debug) where(isnan(thisOctal%jnu(1:maxtrans,subcell))) thisOctal%jnu(1:maxtrans,subcell) = 0.d0
! use updated jnu and collision matrix to determine updated level populations
!                print *, "solving levels"
                call solveLevels(thisOctal%newMolecularLevel(1:maxlevel,subcell), &
                     thisOctal%jnu(1:maxtrans,subcell), dble(thisOctal%temperature(subcell)), &
                     thisMolecule, thisOctal%nh2(subcell), collmatrix(1:maxlevel,1:maxlevel), &
                     ctot(1:maxlevel))
                where(isnan(thisOctal%newMolecularLevel(1:maxlevel,subcell)))
                   thisOctal%newMolecularLevel(1:maxlevel,subcell) = oldpops3
                   endwhere
! accelerate convergence using Ng acceleration on the previous 3 iterations and current.
! Oldpops3 is the previous set
! Because the uppermost levels have the potential to be noisy, we don't use them in the process (length = minlevel -2)
! but they do get modified by the acceleration (1:minlevel)
                if(ng) then
                   if(mod(iter, accstep) .eq. 0) then
                      oldpops4(1:minlevel) = thisOctal%newmolecularLevel(1:minlevel,subcell)
                      if(minlevel .eq. 2) then
                         call ngStep(thisOctal%newmolecularLevel(1:minlevel,subcell), &
                              oldpops1(1:minlevel), oldpops2(1:minlevel), &
                              oldpops3(1:minlevel), oldpops4(1:minlevel), &
                              length = 2)
!!                      else
!                         call ngStep(thisOctal%newMolecularLevel(1:minLevel, subcell), &
!                              oldpops1(1:minlevel), oldpops2(1:minlevel), &
!                              oldpops3(1:minlevel), oldpops4(1:minlevel), &
!                              abs(1.d-60 + (1.d0 / thisoctal%jnu(1:minlevel,subcell))), &
!                              length = max(2,minlevel-1))
                      endif
                   endif
                endif
! Quantify the error - nb. April 2010 - this is recently changed. It will be important that higher levels are converged
! even if they are responsible for only a small fraction of the input to lower levels.
! Oldpops3 is the previous set
                error(1:minlevel) = abs((thisOctal%newMolecularLevel(1:minlevel,subcell) - oldpops3(1:minlevel)) &
                                       / thisOctal%newMolecularLevel(1:minlevel,subcell))
! Old style error
!                error(1:minlevel) = abs((thisOctal%newMolecularLevel(1:minlevel,subcell) - oldpops3(1:minlevel)))

                maxlocerror = maxloc(error(1:minlevel))
                maxerrorloc = maxlocerror(1)
! 3 parameters shoehorned into 1 variable (for space).
                thisoctal%convergence(subcell) = real(100.0 * iter + &
                                                 real(maxerrorloc-1) + &
                                                 min(0.99_db,maxval(error(1:minlevel))) )

! (R)MeanSquare error over minlevel-1 levels. RMS must be less than 1e-10 (MS < 1e-20) 
                fac = sum(error(1:max(2,minlevel-1))**2) ! convergence criterion


!Old style error - absoutle max value.
!                fac = abs(maxval(error(1:max(2,minlevel-1)))) ! convergence criterion
               if (fac < 1.d-12 .or. (iter .eq. maxiter)) then

!                if (fac < 1.d-20 .or. (iter .eq. maxiter)) then
                   popsConverged = .true.
                   if(iter .eq. maxiter) then
                      warncount = warncount + 1
                   endif
! Not entirely sure why you would want to do this. Probably diagnostic. Hence it's commented out.
!                   if(gettau) thisOctal%tau(1:mintrans,subcell) = &
!                        calculatetau(grid, thisoctal, subcell, thismolecule, phi(1:nRay), ds(1:nRay))
                endif
             enddo ! while not converged
          endif ! if no child
       enddo ! all subcells
    enddo ! all octals
 !$OMP END DO

!#ifdef MPI
! end if
!#endif


 deallocate(ds, phi, i0, i0temp, oldpops1, oldpops2, oldpops3, oldpops4, dirweight)
 !$OMP BARRIER
 !$OMP END PARALLEL
           
#ifdef MPI


      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      if(my_rank == 0) write(*,*) "Updating MPI grids"

      call MPI_ALLREDUCE(warncount,warncount_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

      call countVoxels(grid%octreeRoot,nOctal,nVoxels)

      allocate(tArrayd(3,1:nVoxels))
      allocate(tempArrayd(3,1:nVoxels))
      tArrayd = 0.d0
      tempArrayd = 0.d0
      do i = 1, maxlevel
        tArrayd = 0.d0
        call packMoleLevel(octalArray, nVoxels, tArrayd, i,ioctal_beg,ioctal_end)
        do j = 1, 3
           call MPI_ALLREDUCE(tArrayd(j,1:nVoxels),tempArrayd(j,1:nVoxels),nVoxels,MPI_DOUBLE_PRECISION,&
                MPI_SUM,MPI_COMM_WORLD,ierr)
        enddo
        tArrayd = tempArrayd
        call unpackMoleLevel(octalArray, tArrayd, i)
     enddo
     deallocate(tArrayd, tempArrayd)

      allocate(jArrayd(1:nVoxels))
      allocate(tempjArrayd(1:nVoxels))
      jArrayd = 0.d0
      tempjArrayd = 0.d0
      do i = 1, maxtrans
        jArrayd = 0.d0
        call packjnutrans(octalArray, nVoxels, jArrayd, i,ioctal_beg,ioctal_end)
        call MPI_ALLREDUCE(jArrayd(1:nVoxels),tempjArrayd(1:nVoxels),nVoxels,MPI_DOUBLE_PRECISION,&
             MPI_SUM,MPI_COMM_WORLD,ierr)
        jArrayd = tempjArrayd
        call unpackjnutrans(octalArray, jArrayd, i)
     enddo
     deallocate(jArrayd, tempjArrayd)





      if(my_rank == 0) write(*,*) "Done updating"
#else
      warncount_all = warncount
#endif

      ! Calculate convergence towards solution (per level)

      maxRMSFracChange = -1.d30

      write(message,'(a,i6)') "Number of unconverged cells ", warncount_all
      call writeinfo("", FORINFO)

      call writeinfo(message, FORINFO)
      call updateLevels(grid, nvoxels, fixedrays, maxRMSFracChange)

 ! If you think it's converged then test with the same number of rays to make sure

        if(maxRMSFracChange < (tolerance ** 2) * real(nvoxels) .or. allCellsConverged) then
           if (gridConvergedTest) gridConverged = .true.
           gridConvergedTest = .true.
        else
           gridConvergedTest = .false.
           gridConverged = .false.
        endif

        if(.not. fixedrays) then
           call writeAmrGrid(molgridfilename,.false.,grid)
           call writeinfo("",FORINFO)  
           if(writeoutput) then
              open(95, file="restart.dat",status="unknown",form="formatted")
              write(95, '(l1,1x,i7,1x,i3)') fixedrays, nray, grand_iter
              write(95, *) iseedpublic
              close(95)
           endif
        endif
        
        if(writeoutput .and. outputconvergence) then
           call writeinfo("Dumping results", FORINFO)
           call dumpresults(grid, thisMolecule) ! find radial pops on final grid     
        endif

        if(myrankiszero .and. plotlevels .and. .not.(amr1d)) then
           write(filename, '(a,i3.3,a)') "./plots/data",grand_iter,".vtk"
           write(message, *) "Wrote VTK file to", filename
           call writeinfo(message, TRIVIAL)
           call writeVtkFile(grid, filename, "vtk.txt")
           call writeinfo("",FORINFO)  
           call writeinfo("",FORINFO)  
        endif
        
        if(writeoutput .and. dotune) then
           write(message,*) "Done ",nray," rays"
           call tune(6, message)  ! stop a stopwatch
        endif
        
        if(molebench .and. outputconvergence) then
           call torus_mpi_barrier
           call compare_molbench(diffmax)
           write(message,*) "Maximum difference is ", diffmax
           open(96, file="status.dat")
           read(96,*) status
           close(96)
           
           if(status .eq. 1) then 
              call writeinfo("Convergence triggered early!",TRIVIAL)
              gridConverged = .true.
              gridConvergedTest = .true.
           else
              call writeinfo("Convergence test not passed",TRIVIAL)
              call writeinfo(message,TRIVIAL)
           endif
        endif
! Double number of rays if convergence criterion not met and not using fixed rays
!THAW - changed this logical so that doubling occurs for non-fixed rays
!        if (.not.gridConvergedTest) then
        if (.not.gridConvergedTest) then
           if (.not.gridConverged) then               
              if (.not.fixedRays .and. nray < 5000) then
!                 if(mod(ngcounter, accstepgrand) .eq. 0 .and. ngcounter .ne. 0) then
                    nRay = nRay * 2 
!                 endif
              endif
              write(message,'(a,i7,a)') "Now trying ",nRay," Rays"
              call writeInfo(message,FORINFO)
              call writeinfo("", FORINFO)
              if (molebench .and. (nray .gt. maxray)) then 
                 call writeinfo("Molebench Test Ended, Exiting...", TRIVIAL)
                 gridconverged = .true.
              endif
              
           endif
        else
           if(.not. gridconverged) call writeinfo("Doing all rays again to ensure convergence", FORINFO)
        endif

        if (nRay > maxRay) then
           nRay = maxRay  ! stop when it's not practical to do more rays
           call writeWarning("Maximum number of rays exceeded - capping")
        endif
     enddo
  enddo
  
  close(33)
!  666 continue
end subroutine molecularLoop

!!! find the radiation incident at a point in a cell from a pencil beam along a particular direction
   subroutine getRay(grid, fromOctal, fromSubcell, position, direction, ds, phi, i0, thisMolecule, dirweight, fixedrays, &
        warned_neg_dtau, tostar)
!   subroutine getRay(grid, fromOctal, fromSubcell, position, direction, thisMolecule, fixedrays)

     use inputs_mod, only : useDust, realdust, quasi, sourcePos, sourceRadius, sourceTeff


     type(MOLECULETYPE), intent(in) :: thisMolecule
     type(GRIDTYPE), intent(in) :: grid
     type(OCTAL), pointer :: thisOctal, fromOctal
     integer :: fromSubcell
     logical, optional :: fixedrays
     logical, optional :: tostar
     logical :: stage1
     real(double) :: di0(maxtrans), dirWeight
     real(double), intent(out) :: ds, phi, i0(:)
	 	 
     integer, parameter :: maxSamplePoints = 100
	 
     logical :: antithetic

     logical, save :: firsttime = .true.
     logical, save :: conj = .false.
     type(VECTOR), save :: possave, dirsave
     real(double), save :: rsave, s
	 	 	 
     real(double), allocatable, save :: OneArray(:)
     real(double), save :: OneOverNTauArray(maxsamplePoints)
     real(double), save :: BnuBckGrnd(128)
	 
     integer :: subcell

     integer :: nTau     
     integer :: ilambda, iTrans, iTau

     real(double) :: alphanu(maxtrans,2), alpha(maxtrans),  snu(maxtrans), jnu(maxtrans)
     real(double) :: tau(maxtrans), dTau(maxtrans), localradiationfield(maxtrans), attenuation(maxtrans), kappaAbs
     real(double) :: balance(maxtrans), spontaneous(maxtrans)
	 	  
     type(VECTOR) :: position, currentPosition, thisPosition, endposition, direction, halfstep
     type(VECTOR) :: startVel, endVel, thisVel, rayvel

     real(double) :: dv, deltaV
     real(double) :: dist, dds, tval, dvAcrossCell, projVel, endprojVel, CellEdgeInterval, phiProfVal
     real(double) :: r, nmol 

     logical :: warned_neg_dtau

!$OMP THREADPRIVATE (firstTime, conj, possave, dirsave, rsave, s, oneArray, oneOVerNTauArray, BnuBckGrnd)
     antithetic = .false.


!Set/initialise runtime parameters
     if(firsttime) then
        if(.not. allocated(onearray)) then
           allocate(OneArray(maxtrans))
        end if
        OneArray = 1.0d0
       
        do itrans = 1, maxtrans
           BnuBckGrnd(itrans) = Bnu(thisMolecule%transfreq(itrans), Tcbr)
        enddo

        do ntau = 2, maxsamplepoints
           OneOverNtauArray(ntau) = 1.d0 / (dble(nTau) - 1.d0)
        enddo
        
        firsttime =.false.
    endif

    if(tostar) then       
        do itrans = 1, maxtrans
           BnuBckGrnd(itrans) = Bnu(thisMolecule%transfreq(itrans), sourceTeff(1))
        enddo        
        firsttime=.true.
    end if


!Set up antithetic variates to improve convergence
    if(antithetic) then
       
       if(conj) then
          position = possave
       else
          position = randomPositionInCell(fromOctal, fromsubcell)
          possave = position
       endif
    else
       position = randomPositionInCell(fromOctal, fromsubcell)
    endif

    if(present(fixedrays)) then
       stage1 = fixedrays
    else
        stage1 = .false.
     endif

!Choose position in cell & determine velocity at that position
     position = randomPositionInCell(fromOctal, fromsubcell)
     thisOctal => fromOctal
     subcell = fromSubcell
     rayVel = Velocity(position, grid, thisOctal, subcell)


!Generate random direction and frequency using quasi/pseudo random number generators
     dirWeight = 1.d0
     if(quasi) then     
        call sobseq(r1)    
        direction = specificUnitVector(r1(1),r1(2))
        call normalize(direction)
        deltaV = 4.3 * thisOctal%microturb(subcell) * (r1(3) - 0.5d0) ! random frequency near line centre
     else
        if(antithetic) then

           if(tostar) then
              call randomNumberGenerator(getDouble=r)
              direction%x = sourcePos(1)%x - position%x
              direction%x = sourcePos(1)%y - position%y
              direction%x = sourcePos(1)%z - position%z
              call normalize(direction)
              deltaV = 4.3 * thisOctal%microturb(subcell) * (r - 0.5d0) ! random frequency near line spectrum peak.
              dirWeight = stellarRayWeight(sourcePos(1), position, sourceRadius(1))
!           if(tostar .and. conj) then
!              r = 1. - rsave
!              direction = (-1.d0) * dirsave
!              deltaV = 2.15 * thisOctal%microturb(subcell) * r ! random frequency near line spectrum peak. 
!              dirWeight = stellarRayWeight(sourcePos(1), position, sourceRadius(1))
!              call randomNumberGenerator(getdouble=s)
!              if(s > 0.5) deltaV = deltaV * (-1.d0)
!
!           else if (tostar .and. .not. conj) then
!              call randomNumberGenerator(getDouble=r)
!              rsave = r
!              direction = sourcePos(1)              
!              dirsave = direction
!              deltaV = 2.15 * thisOctal%microturb(subcell) * r ! random frequency near line spectrum peak. 
!              dirWeight = stellarRayWeight(sourcePos(1), position, sourceRadius(1))
!              call randomNumberGenerator(getdouble=s)
!              if(s > 0.5) deltaV = deltaV * (-1.d0)

           elseif(conj) then
              r = 1. - rsave
              direction = (-1.d0) * dirsave
  
              deltaV = 2.15 * thisOctal%microturb(subcell) * r ! random frequency near line spectrum peak.
              call randomNumberGenerator(getdouble=s)
               if(s > 0.5) deltaV = deltaV * (-1.d0)
           else
              call randomNumberGenerator(getDouble=r)
              rsave = r
              direction = randomUnitVector()
              dirsave = direction

              deltaV = 2.15 * thisOctal%microturb(subcell) * r
              
              call randomNumberGenerator(getdouble=s)
              if(s > 0.5) deltaV = deltaV * (-1.d0)
           endif
        else ! pseudorandom
           !thaw - force a ray to a source
           if(tostar) then
              call randomNumberGenerator(getDouble=r)
              direction%x = sourcePos(1)%x - position%x
              direction%y = sourcePos(1)%y - position%y
              direction%z = sourcePos(1)%z - position%z
              call normalize(direction)
              deltaV = 4.3 * thisOctal%microturb(subcell) * (r - 0.5d0) ! random frequency near line spectrum peak. 
              dirWeight = stellarRayWeight(sourcePos(1), position, sourceRadius(1))
           else
              call randomNumberGenerator(getDouble=r)
              direction = randomUnitVector() ! Put this line back if you want to go back to pseudorandom
              deltaV = 4.3 * thisOctal%microturb(subcell) * (r - 0.5d0) ! random frequency near line spectrum peak. 
           end if
        endif
        conj = .not. conj
     endif

!!! SCIENCE LINES

!     deltaV = deltaV + (rayVel .dot. direction) ! transform to take account of grid velocity
!     projVel = deltaV - (rayVel .dot. direction) ! transform back to velocity relative to local flow

!!! COMPUTATIONALLY QUICKER TO DO THIS
!     print *, "direction ", direction 
     projvel = deltaV
     deltaV =  deltaV + (rayVel .dot. direction)

     call distanceToCellBoundary(grid, position, direction, ds, sOctal=thisOctal)

! Get to cell boundary
     currentPosition = position + ds * direction

! Get weighting for ray (see sumphi in calculatejbar)
     if (inOctal(grid%octreeRoot, currentPosition)) then ! check that we're still on the grid
        thisVel = velocity(currentposition, grid) 
        endprojVel = deltaV - (thisVel.dot.direction)
        phi = phiProf((endprojVel+projvel)*0.5d0, thisOctal%molmicroturb(subcell))
     else
        phi = phiProf(projVel, thisOctal%molmicroturb(subcell)) ! if fell off grid then assume phi is unchanged
     endif

! Initialise more variables
     ds = ds * 1.d10 ! convert from torus units to cm for use in calculatejbar
     i0 = 0.d0
     tau = 0.d0

! Follow long characteristic of ray until edge of grid
     do while(inOctal(grid%octreeRoot, currentPosition))

! Find distance to edge of cell
        call findSubcellLocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, soctal = thisoctal)
		
! Determine how many velocity samples should be taken per grid cell.
! If fixedrays then always check.
! If ever greater than 2 then always calculate in this cell for rest of calculation.
        if(fixedrays .or. thisOctal%nsplit(subcell) .gt. 2) then
           startVel = velocity(currentposition, grid) 
           endPosition = currentPosition + tval * direction
           
           if(inOctal(grid%octreeRoot, endposition)) then
              endVel = velocity(endposition, grid)
           else
              endVel = startvel
           endif
       
           dvAcrossCell = ((startVel - endvel) .dot. direction) ! start - end
           dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))
           
           nTau = min(max(2, nint(dvAcrossCell * 5.d0)), maxSamplePoints) 

           !THAW - TEMPORARY, USE ABOVE
!           nTau = maxSamplePoints
           ! selects dVacrossCell as being between 0.1 and 10 (else nTau bounded by 2 and 200)

           if(ntau .gt. 2) thisOctal%nsplit(subcell) = ntau
        else
           ntau = 2
        endif

        CellEdgeInterval = OneOverNtauArray(ntau) ! if ntau = 2 then the segment is taken in one chunk
        dds = tval * cellEdgeInterval ! determine line segment length
        halfstep = dds * 0.5 * direction ! find point halfway between start and end to take representative velocity
        dist = 0.d0
					
! Calculate absorption from (various types of) dust 		
        if(useDust) then     
           if(realdust) then
              do itrans = 1, maxtrans               
                 call locate(grid%lamArray, size(grid%lamArray), real(lambda(itrans)), ilambda)
                 call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = real(lambda(itrans)), &
                      kappaAbs = kappaAbs)
                 alphanu(itrans,2) = kappaAbs * 1.d-10 !torus units -> cm !1d-10 is right 
             enddo
           else
              if(grid%geometry .eq. 'agbstar') then
                 do itrans = 1, maxtrans
                    kappaAbs = 0.1 * thisMolecule%transfreq(itrans) * 1e-12 * 1d10 !multiplied by density !cm -> torus
                    kappaAbs = kappaAbs * thisOctal%rho(subcell)
                    alphanu(itrans,2) = kappaAbs * 1.d-10 !torus units -> cm !1d-10 is right
                 enddo
              else
                 do itrans = 1, maxtrans
                    kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10
                    kappaAbs = kappaAbs * thisOctal%rho(subcell)
                    !multiplied by density !cm -> torus
                    alphanu(itrans,2) = kappaAbs * 1.d-10 !torus units -> cm !1d-10 is right
                 enddo
              endif
           endif
!           do itrans = 1, maxtrans
!              alphanu(itrans,2) = kappaAbs * 1.d-10 !torus units -> cm !1d-10 is right
!              !kappa already multiplied by density in amr_mod
!           end do
        endif

! Calculate number density of molecules
        nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
! Absorption by gas per length
        balance(1:maxtrans) = (hcgsOverFourPi * nmol) * &
             (thisOctal%molecularLevel(ilower(1:maxtrans),subcell) * thisMolecule%einsteinBlu(1:maxtrans) - &
             thisOctal%molecularLevel(iupper(1:maxtrans),subcell) * thisMolecule%einsteinBul(1:maxtrans))
! Emission by gas per length        
!        if(.not. thisOctal%ghostcell(subcell)) then
        spontaneous(1:maxtrans) = (hCgsOverfourPi * nmol) * &
             thisMolecule%einsteinA(1:maxtrans) * thisOctal%molecularLevel(iupper(1:maxtrans),subcell)
!        else
!           spontaneous(1:maxtrans) = 0.d0
!        end if
! Source function 
        where (balance /= 0.d0)
           snu(1:maxtrans) = spontaneous(1:maxtrans) / balance(1:maxtrans)
        end where
! Calculate total emission and absorption over entire cell (sum over all line segments)
        do itau = 2, nTau
! Get next position	
           dist = dist + dds
           thisPosition = currentPosition + dist * direction
! Get velocity at midpoint   
           thisVel = velocity(thisPosition-halfstep, grid)
! Get velocity difference and weighting
           dv = deltaV - (thisVel .dot. direction)
           PhiProfVal = phiProf(dv, thisOctal%molmicroturb(subcell))
! Get weighted gas absorption (+dust if req). If usedust then update jnu and snu            
           if(usedust) then
              alphanu(1:maxtrans,1) = phiprofval * balance(1:maxtrans) 
              alpha(1:maxtrans) = alphanu(1:maxtrans,1) + alphanu(1:maxtrans,2)
              jnu(1:maxtrans) = spontaneous(1:maxtrans) * phiProfVal + &
                   alphanu(1:maxtrans,2) * thisOctal%bnu(1:maxtrans,subcell) 
              
              do itrans = 1,maxtrans
                 if(alpha(itrans) .ne. 0.d0) then 
                    snu(itrans) = jnu(itrans) / alpha(itrans)
                 else
                    snu(itrans) = tiny(snu)
                    alpha(itrans) = tiny(alpha)
                 endif
              enddo
           else
              alpha(1:maxtrans) = balance(1:maxtrans) * phiprofVal 

!THAW - this stuff is unnecessary because it is handled prior to this do loop
!              jnu(1:maxtrans) = spontaneous(1:maxtrans) * phiProfVal
!              
!              do itrans = 1,maxtrans
!                 if(alpha(itrans) .ne. 0.d0) then 
!                    snu(itrans) = jnu(itrans) / alpha(itrans)
!                 else
!                    snu(itrans) = tiny(snu)
!                    alpha(itrans) = tiny(alpha)
!                 endif
!              enddo

           endif

! Calculate optical depth and associated attenuation factor
           dTau = alpha * dds * 1.d10
           attenuation = exp(-tau)
! Intensity along ray owing to line segment
           do iTrans = 1, maxtrans
              if (dtau(itrans) < -1.d0 .and. .not. warned_neg_dtau ) then
                 write(*,*) "dtau ERROR: ",dtau(itrans),alpha(itrans),dds,phiprofval
                 warned_neg_dtau = .true. 
              endif
           enddo
           localradiationfield = exp(-dtau) 
           localradiationfield = OneArray - localradiationfield
           di0 = localradiationfield * snu
! Add contribution to whole line           
           i0 = i0 + attenuation * di0
           tau = tau + dtau
        enddo
! Update position into next cell making sure that new position is definitely in new cell
        currentPosition = currentPosition + (tval + 1.d-3*grid%halfSmallestSubcell) * direction 
     enddo

! Add attenuated CMB to intensity at point.     
     if((grid%geometry .eq. "h2obench1") .or. (grid%geometry .eq. "h2obench2")) then
        continue
     else
        i0(1:maxtrans) = i0(1:maxtrans) + BnuBckGrnd(1:maxtrans) * attenuation(1:maxtrans)
     endif
   end subroutine getRay


!thaw - return the directional weight associated with a ray towards a star
   
   real(double) function stellarRayWeight(sourcePos, position, rstar)
     type(VECTOR) :: sourcePos
     type(VECTOR) :: position
     real(double) :: distance, rstar, theta, omega_star
     
     distance = (sourcePos%x-position%x)**2 + (sourcePos%y-position%y)**2 + (sourcePos%z-position%z)**2
     distance = sqrt(distance)
     
     theta = tan((rstar*rsol)/distance)
     
!solid angle subtended by star
     omega_star = 2.d0*pi*(1.d0-cos(theta))
     
     stellarRayWeight = omega_star/fourpi

   end function stellarRayWeight


! This subroutine calculates the average local radiation field in a cell, Jbar. So the SE can be worked out
   subroutine calculateJbar(nr, grid, thisOctal, subcell, thisMolecule, tempds, tempphi, i0, nPops, jbar, dirweight)

     use inputs_mod, only : useDust, realdust, debug

     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(MOLECULETYPE) :: thisMolecule
     integer :: nr
     real(double), intent(in) :: tempds(:), tempphi(:), i0(:,:), nPops(:), dirweight(:)
     real(double), intent(out) :: jbar(:)
     real(double), allocatable :: phids(:), temptauArray(:), opticaldepthArray(:), otp(:), &
          jbarinternalArray(:), jbarExternalArray(:)
     integer :: iTrans

     real(double) :: nLower(maxtrans), nUpper(maxtrans), nMol, etaline(maxtrans), snugas(maxtrans), &
          kappaAbs, alphanuBase(maxtrans), jnuDust(maxtrans)
     real(double), allocatable :: alphanu(:,:),  alpha(:), &
          jnu(:), snu(:)

     real(double) :: sumPhi

     allocate(phids(1:nray))
     allocate(tempTauArray(1:nRay))
     allocate(opticalDepthArray(1:nray))
     allocate(otp(1:nray))
     allocate(snu(1:nray)) !THAW
     allocate(jbarInternalArray(1:nray))
     allocate(jbarExternalArray(1:nray))
     allocate(alphanu(1:nray,1:2))
     allocate(alpha(1:nRay))
     allocate(jnu(1:nray))

     
! Get total number density of molecules in each level
     nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
     nLower(1:maxtrans) = nPops(iLower(1:maxtrans)) * nMol
     nUpper(1:maxtrans) = nPops(iUpper(1:maxtrans)) * nMol


! Determine gas emission and absorption and source function
     etaLine(1:maxtrans) = thisMolecule%einsteinA(1:maxtrans) * nUpper(1:maxtrans)
     alphanuBase(1:maxtrans) =(nLower(1:maxtrans) * thisMolecule%einsteinBlu(1:maxtrans) - &
                  nUpper(1:maxtrans) * thisMolecule%einsteinBul(1:maxtrans))
     where(alphanubase .ne. 0) 
        snugas(1:maxtrans) = etaline(1:maxtrans) / alphanuBase(1:maxtrans)
     endwhere

! Calculate phids and sum(phi) as required by Jbar equation (for integral)
     phids(1:nr) = tempphi(1:nr) * tempds(1:nr)
     sumPhi = sum(tempphi(1:nray))

     alphanu(:,2) = 0.d0
     jnudust(:) = 0.d0

     if(useDust) then
        do itrans = 1, maxtrans
           call locate(grid%lamArray, size(grid%lamArray), real(lambda(itrans)), ilambda)

           if(realdust) then
              call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = real(lambda(itrans)), &
                               kappaAbs = kappaAbs)
           elseif(grid%geometry .eq. 'agbstar') then
              kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10
              !multiplied by density !cm -> torus
           endif

           alphanu(itrans,2) = kappaAbs * 1.d-10

           if(associated(thisoctal%temperaturedust)) then
              jnuDust(itrans) = alphanu(itrans,2) * &
                                bnu(thisMolecule%transfreq(iTrans), dble(thisOctal%temperaturedust(subcell)))
           else
              jnuDust(itrans) = alphanu(itrans,2) * &
                                bnu(thisMolecule%transfreq(iTrans), dble(thisOctal%temperature(subcell)))
           endif

           
           jnu(1:nray) = etaLine(itrans) * tempphi(1:nRay) * hCgsOverFourPi!/thisMolecule%transFreq(iTrans)

           alphanu(1:nray,1) = alphanuBase(itrans) * tempphi(1:nray) * hCgsOverFourPi!/thisMolecule%transFreq(iTrans)
           alpha(1:nray) = alphanu(1:nray,1) + alphanu(itrans,2)

           if(alpha(itrans) .ne. 0) then
              snu(1:nray) = (jnu(1:nray) + jnudust(itrans)) / alpha(1:nray)
           else
              snu = tiny(snu)
           endif

           temptauArray(1:nRay) = alpha(1:nray) * tempds(1:nRay)
           opticaldepthArray(1:nRay) = exp(-1.d0 * temptauArray(1:nRay))
           otp(1:nray) = opticaldepthArray(1:nRay) * tempphi(1:nRay) ! intermediate stage (weighted opt depth)
! external jbar due to getray (i0)        
           jBarExternalArray(1:nRay) = i0(1:nray,itrans) * otp(1:nray)
! internal jbar inside cell
           jBarInternalArray(1:nRay) = snu(1:nray) * (tempphi(1:nray) - otp(1:nray))
! jbar weighted average of int + ext        
!thaw - multpilied by a direction weighting (should be 1 by default unless forcing towards stars)
           jbar(itrans) = (sum(jBarExternalArray(1:nRay)*dirWeight(1:nRay)+jBarInternalArray(1:nRay)*dirWeight(1:nRay))) / sumPhi
        enddo
     else
! vecorised jbar calculation
        do itrans = 1, maxtrans
           alpha(1:nray) = alphanuBase(itrans) * hCgsOverFourPi!/thisMolecule%transFreq(iTrans)

! calculate optical depth within cell based on average over all rays       
           temptauArray(1:nRay) = alpha(1:nray) * phids(1:nRay)
!write(*,*) "tta", temptauarray(1:nray), sum(temptauarray(1:nray))
           opticaldepthArray(1:nRay) = exp(-1.d0 * temptauArray(1:nRay))
           otp(1:nray) = opticaldepthArray(1:nRay) * tempphi(1:nRay) ! intermediate stage (weighted opt depth)
! external jbar due to getray (i0)        
           jBarExternalArray(1:nRay) = i0(1:nray,itrans) * otp(1:nray) 
! internal jbar inside cell
           jBarInternalArray(1:nRay) = snugas(itrans) * (tempphi(1:nray) - otp(1:nray))
! jbar weighted average of int + ext        
           jbar(itrans) = (sum(jBarExternalArray(1:nRay)*dirWeight(1:nRay)+jBarInternalArray(1:nRay)*dirWeight(1:nRay))) / sumPhi
        enddo
     endif
     if(debug) where(isnan(jbar)) jbar = 0.d0
   
     deallocate(phids)
     deallocate(tempTauArray)
     deallocate(opticalDepthArray)
     deallocate(otp)
     deallocate(jbarInternalArray)
     deallocate(jbarExternalArray)
     deallocate(alphanu)
     deallocate(jnu)


   end subroutine calculateJbar

 ! solves rate equation in matrix format
   subroutine solveLevels(nPops, jnu,  temperature, thisMolecule, nh2, collmatrix, ctot)
     
     use inputs_mod, only : debug
     
     real(double), intent(inout) :: nPops(:)
     real(double), intent(in) :: temperature
     real(double), intent(in) :: jnu(:)
     type(MOLECULETYPE), intent(in) :: thisMolecule

     real(double) :: nh2
     real(double) :: blujnu(maxtrans), buljnu(maxtrans)
     real(double) :: matrixA(maxlevel+1,maxlevel+1), matrixB(maxlevel+1,1)
     real(double) :: matrixAsave(maxlevel+1,maxlevel+1), matrixBsave(maxlevel+1,1)
     real(double) :: matrixArad(maxlevel+1,maxlevel+1)
     real(double), intent(in) :: collMatrix(:,:), ctot(:)

     integer :: l, k, itrans

     logical :: dummy
     logical :: luslvOK

     matrixA = 1.d-60 ! Initialise rates to negligible to avoid divisions by zero
     matrixB = 0.d0 ! Solution vector - all components (except last) => equilibrium

! Sum over all levels = 1 - Conservation constraint
     matrixB(maxLevel+1,1) = 1.d0 

! intermediate step
     Buljnu(1:maxtrans) = thisMolecule%einsteinBul(1:maxtrans) * jnu(1:maxtrans)
     Blujnu(1:maxtrans) = thisMolecule%einsteinBlu(1:maxtrans) * jnu(1:maxtrans)

 ! This do loop calculates the contribution to transition rates of each level from every other level. 
 ! NB emission is +ve here

     do iTrans = 1, maxtrans
        k = iUpper(itrans)
        l = iLower(itrans)

! total emission (+) from upper level into (-) all lower levels (stored on-diagonal) 
        matrixA(k,k) = matrixA(k,k) + buljnu(itrans) + thisMolecule%einsteinA(iTrans)

! total emission (+) from upper level, k into (-) lower level, l (stored off-diagonal)
!(-ve emission stored in l,k)
        matrixA(l,k) = matrixA(l,k) - buljnu(itrans) - thisMolecule%einsteinA(iTrans)

! stimulated emission (+) from lower level, l into all upper levels
        matrixA(l,l) = matrixA(l,l) + blujnu(itrans)

! stimulated emission (+) from lower level, l into upper level, k
!(-ve emission stored in k,l)
        matrixA(k,l) = matrixA(k,l) - blujnu(itrans)
     enddo
!write(*,*) "matrix",matrixA(:,:)
     if(debug) matrixArad = matrixA

! Add contributions from getcollmatrix
     matrixA(1:maxlevel,1:maxlevel) = matrixA(1:maxlevel,1:maxlevel) - collmatrix(1:maxlevel,1:maxlevel)
!write(*,*) "matrixN",matrixA(:,:)
     matrixA(maxlevel+1,1:maxlevel+1) = 1.d0 ! sum of all level populations
     matrixA(1:maxlevel+1,maxlevel+1) = 1.d-60 ! fix highest population to small non-zero value
!write(*,*) "matrixNN",matrixA(:,:)
!     do i = 1, maxLevel
!        write(*,'(100(1pe9.2))') matrixA(i,1:maxLevel)
!     enddo
     
     if(debug) then
        matrixAsave = matrixA
        matrixBsave = matrixB
     endif

! finished creating equation 10, now solve it to find new level populations using lu solver
#ifdef USEMKL
     call gesv(matrixA, matrixB)
#else
     if(maxlevel .gt. 2) then
        call luSlv(matrixA, matrixB(:,1),luslvOK)
     else
        call GAUSSJ(matrixA, maxlevel+1, maxlevel+1,matrixB(:,1), maxlevel+1, maxlevel+1,dummy)
     endif
#endif
!write(*,*) "matrixS",matrixA(:,:)
!write(*,*) "matrixSS",matrixB

! final step - if level population candidates have a problem then fix by 
! 1st) repeating solution over minlevels 
!     -these are the important levels to get right - remaining levels are set equal to 1-npops
! 2nd) setting npops = LTEpops if collisionally dominated
! 3rd) setting npops = radiative solution if radiatively dominated


     if ( (debug .and. any(isnan(matrixB))) .or. (.not.luslvOK) ) then
        matrixAsave(minlevel+1,1:minlevel+1) = 1.d0 ! sum of all level populations
        matrixAsave(1:minlevel+1,minlevel+1) = 1.d-60 ! fix highest population to small non-zero value
        matrixBsave = 1.d-60 ! Solution vector - all components (except last) => equilibrium ! used to be 1d-10
        matrixBsave(minLevel+1,1) = 1.d0 ! Sum over all levels = 1 - Conservation constraint

#ifdef USEMKL
        call gesv(matrixAsave(1:minlevel+1,1:minlevel+1), matrixBsave(1:minlevel+1,1))
#else
        call luSlv(matrixAsave(1:minlevel+1,1:minlevel+1), matrixBsave(1:minlevel+1,1),luslvOK)
#endif

        if ( (debug .and. .not. any(isnan(matrixBsave))) .or. luslvOK ) then
           nPops(1:minlevel) = matrixBsave(1:minlevel,1)
           npops(minlevel+1:maxlevel) = abs(1.d0 - sum(matrixBsave(1:minlevel,1)))
           npops = abs(npops / sum(npops))
           write(66,*) "error fixed 1 ", npops(minlevel), npops(maxlevel)
        else
           if(nh2 * 2.d0 * mhydrogen .gt. maxval(thismolecule%einsteinA(1:maxlevel) / ctot(1:maxlevel))) then
              call LTEpops(thisMolecule, temperature, npops)
              write(66,*) "error fixed 2 ", npops(minlevel), npops(maxlevel)
           else
              matrixA = matrixArad
              matrixArad(maxlevel+1,1:maxlevel+1) = 1.d0 ! sum of all level populations
              matrixArad(1:maxlevel+1,maxlevel+1) = 1.d-60 ! fix highest population to small non-zero value
              
              matrixB = 1.d-60 ! Solution vector - all components (except last) => equilibrium ! used to be 1d-10
              matrixB(maxLevel+1,1) = 1.d0 ! Sum over all levels = 1 - Conservation constraint

#ifdef USEMKL
              call gesv(matrixArad, matrixB)
#else
              call luSlv(matrixArad, matrixB(:,1),luslvOK)
#endif        
              if ( (debug .and. .not. any(isnan(matrixB))) .or. luslvOK ) then
                 matrixB = abs(matrixB) ! stops negative level populations causing problems
                 write(66,*) "error fixed 3 ", npops(minlevel), npops(maxlevel), &
                      maxval(thismolecule%einsteinA(1:maxlevel) / ctot(1:maxlevel)), nh2 * 2.d0 * mhydrogen
              else
                 write(66,*) "error still 4 ", npops(minlevel), npops(maxlevel)
              endif
! level populations out
              nPops(1:maxlevel) = matrixB(1:maxLevel,1)
           endif
        endif

     else ! standard case i.e. no problem with the first call to luSlv 
        matrixB = abs(matrixB) ! stops negative level populations causing problems
! level populations out
        nPops(1:maxlevel) = matrixB(1:maxLevel,1)
     endif

   end subroutine solveLevels

 ! Calculate collision rates between partners for given temperature
   real(double) function collRate(thisMolecule, temperature, iPart, iTrans)
     type(MOLECULETYPE) :: thisMolecule
     real(double), intent(in) :: temperature
     integer, intent(in) :: iTrans, iPart

     real(double) :: r
     integer :: k
if(thismolecule%ncolltemps(1) .gt. 1) then
! linearly interpolate or extrapolate (downwards)
     if(temperature .gt. thisMolecule%collTemps(iPart,1)) then

        call locate(thisMolecule%collTemps(iPart,1:thisMolecule%nCollTemps(iPart)), &
             thisMolecule%nCollTemps(iPart), temperature, k)

        r = (temperature - thisMolecule%collTemps(iPart,k)) / &
             (thisMolecule%collTemps(iPart,k+1) - thisMolecule%collTemps(iPart, k))

        collRate = thisMolecule%collRates(iPart, iTrans, k) + &
             r * ( thisMolecule%collRates(iPart, iTrans, k+1) - thisMolecule%collRates(iPart, iTrans, k))

     else
        
        r = (temperature - thisMolecule%collTemps(iPart,1)) / &
             (thisMolecule%collTemps(iPart,2) - thisMolecule%collTemps(iPart, 1))

!        collRate = max(thisMolecule%collRates(iPart, iTrans, 1), thisMolecule%collRates(iPart, iTrans, 1) + &
!             r * (thisMolecule%collRates(iPart, iTrans, 2) - thisMolecule%collRates(iPart, iTrans, 1)))

        collRate = thisMolecule%collRates(iPart, iTrans, 1) + &
             r * ( thisMolecule%collRates(iPart, iTrans, 2) - thisMolecule%collRates(iPart, iTrans, 1))

        if(collrate .lt. 0.d0) collrate = thisMolecule%collRates(iPart, iTrans, 1)
     endif
else
   collrate = thisMolecule%collRates(1, 1, 1)
endif
   end function collRate

   function collPartnerDensity(thisMolecule, ipart, nh2, nH, nHe, ne, nProtons) result(nx)
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: nx
     integer :: iPart
     real(double) :: nh2, ne, nh, nHe, nProtons

! note we assume an ortho:para ratio of 5.e-5
! http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2004A%26A...418.1035W&db_key=AST

     select case (thisMolecule%iCollPartner(iPart))
        case(1) ! all H2
           nx = nH2
        case(2) ! para-H2
           nx = nH2 * (1.d0-5.d-5)
        case(3) ! ortho-H2
           nx = nH2 * 5.d-5
        case(4) ! electrons
           nx = ne
           call writeWarning("collPartnerDensity: Ne not yet implemented as a collision partner")
        case(5) ! H atom
           nx = nH
           call writeWarning("collPartnerDensity: N(H) not yet implemented as a collision partner")
        case(6) ! He atom
           nx = nHe
           call writeWarning("collPartnerDensity: N(He) not yet implemented as a collision partner")
        case(7) ! H+
           nx = nProtons
           call writeWarning("collPartnerDensity: N(H+) not yet implemented as a collision partner")
        case DEFAULT
           call writeFatal("collRate: collision partner type not recognised")
     end select

   end function collPartnerDensity

   subroutine getCollMatrix(nh2, temperature, thismolecule, nh, nhe, ne, nprotons, collT, ctot)

     type(moleculetype) :: thismolecule

     integer :: iPart, iTrans
     integer :: i,l,k
     
     real(double) :: colldeEx, collEx, boltzFac
     real(double) :: nh2, temperature, nh, nhe, ne, nprotons
     real(double) :: OneOverkT
! collT is what goes out and gets subtracted from matrixA
     real(double) :: collMatrix(maxlevel, maxlevel), collT(:,:), ctot(:)

     collMatrix = 1.d-60
     cTot = 1d-60

     oneOverkT = 1.d0 / (kev * temperature)

     do iPart = 1, thisMolecule%nCollPart
        do iTrans = 1, thisMolecule%nCollTrans(iPart)

! no need to calculate contributions from/to levels gt maxlevel
           l = thisMolecule%iCollLower(iPart, iTrans)
           k = thisMolecule%iCollUpper(iPart, iTrans)
           if(l .gt. maxlevel .or. k .gt. maxlevel) cycle
!           if(l .gt. maxlevel) cycle

! boltzman factor connects collisional excitation and deexcitation (along with stat weighting)          
           boltzFac = exp(-abs(thisMolecule%energy(k)-thisMolecule%energy(l)) * OneOverKT)
! deexcition rate calculated from rates and coll partner density
           colldeEx = collRate(thisMolecule, temperature, iPart, iTrans) * &
                collPartnerDensity(thisMolecule, ipart, nh2, nH, nHe, ne, nProtons)
           collEx = colldeEx * boltzFac * thisMolecule%g(k) / thisMolecule%g(l)
!write(*,*) "PANIC", colldeex, collex
! emission from lower levels to upper levels 
           collMatrix(l, k) = collMatrix(l, k) + collEx
! emission from upper to lower
           collMatrix(k, l) = collMatrix(k, l) + colldeEx
        enddo
     enddo

! sum over all collisional rates out of each level (ctot). Add to radiative
     collT = transpose(collmatrix)
     cTot(:) = cTot(:) + sum(collT(:,:), dim = 1)
!the minus can be used to compensate for double counting - old code below.
     forall(i=1:maxlevel) 
        collT(i,i) = collT(i,i) - ctot(i)
     end forall

!     ctot = 0.d0
! this code used to be in solvelevels hence the matrixA
!    do k = 1, maxlevel
!       do l = 1, maxlevel
!          cTot(k) = cTot(k) + collMatrix(k,l)
!       enddo
!    enddo

!     do k = 1, maxlevel
!        do l = 1, maxlevel
!           if (k .ne. l) then
!! subtract deexcitations from upper to lower
!              matrixA(k, l) = matrixA(k, l) - collMatrix(l, k)
!           endif
!        enddo
!        matrixA(k,k) = matrixA(k,k) + cTot(k)
!     enddo

   end subroutine getCollMatrix

! This subroutine calculates the maximum fractional change in the first 6 energy levels
! 
! The compiler directive prevents O3 level optimisation with ifort due to a run time 
! segmentation fault (seen with ifort 12.1.0) D. Acreman 21/11/12
   recursive subroutine  swapPops(thisOctal, maxFracChangePerLevel, avgFracChange, counter, &
                          fixedrays)
!DEC$ OPTIMIZE:2
     use inputs_mod, only : tolerance, dongstep

     type(octal), pointer   :: thisOctal
     real(double) :: maxFracChangePerLevel(:),  avgFracChange(:,:)
     integer :: counter(:,:)

     logical :: fixedrays

     type(octal), pointer  :: child 
     integer :: subcell, i, j

     real(double) :: oldpops1(maxlevel), oldpops2(maxlevel), oldpops3(maxlevel), oldpops4(maxlevel)
   
     real(double) :: newFracChangePerLevel(minlevel), temp(minlevel-1), maxFracChange

     logical :: ng
     ng = dongstep

! traverse the grid
     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call swapPops(child, maxFracChangePerLevel, avgFracChange, counter, &
                      fixedrays)
                 exit
              end if
           end do
        else
! if doing ng Acceleration step the do it here
           if(ng) then
! If fixedrays then only do every accstepgrand steps
              if((mod(ngcounter, accstepgrand) .eq. 0 .and. ngcounter .ne. 0) &
!                   .or. nray .eq. 2048*initnray & 
!                   .or. nray .eq. 128 *initnray &
!                   .or. nray .eq. 8   *initnray) then
                   ) then
! Put old levels into temporary arrays
                 oldpops1(1:minlevel) = thisOctal%oldestmolecularLevel(1:minlevel,subcell)
                 oldpops2(1:minlevel) = thisOctal%oldmolecularLevel(1:minlevel,subcell)
                 oldpops3(1:minlevel) = thisOctal%molecularLevel(1:minlevel,subcell)
                 oldpops4(1:minlevel) = thisOctal%newmolecularLevel(1:minlevel,subcell)
! Calculate ng accelerated levels upto minlevel - 2 (avoids numerics noise)
                 if(minlevel .eq. 2) then
                    call ngStep( &
                    thisOctal%newmolecularLevel(1:minlevel,subcell),  &
                         oldpops1(1:minlevel), oldpops2(1:minlevel), &
                         oldpops3(1:minlevel), oldpops4(1:minlevel), &
                         length = 2)
!                 elseif(fixedrays) then
!                    call ngStep(&
!                         thisOctal%newmolecularLevel(1:minlevel,subcell), &
!                         oldpops1(1:minlevel), oldpops2(1:minlevel), &
!                         oldpops3(1:minlevel), oldpops4(1:minlevel), &
!                         abs(1.d-60 + (1.d0 / thisoctal%jnu(1:minlevel,subcell))), &
!                         length = max(2,minlevel-1), doubleweight = .false.)
!                 else
!                    call ngStep(&
!                         thisOctal%newmolecularLevel(1:minlevel,subcell), &
!                         oldpops1(1:minlevel), oldpops2(1:minlevel), &
!                         oldpops3(1:minlevel), oldpops4(1:minlevel), &
!                         abs(1.d-60 + (1.d0 / thisoctal%jnu(1:minlevel,subcell))), &
!                         length = max(2,minlevel-1), doubleweight = .false.)
                 endif
              endif
           endif
           
! determine updated fractional change in each level upto minlevels.
           newFracChangePerLevel = abs(((thisOctal%newMolecularLevel(1:minlevel,subcell) - &
                thisOctal%molecularLevel(1:minlevel,subcell)) / &
                (thisOctal%newmolecularLevel(1:minlevel,subcell) + 1d-30)))

           where ( newFracChangePerLevel(1:minlevel) <= 0.0 ) newFracChangePerLevel = 1.0e-30_db
! Cleverly store convergence in an integer array to reduce the size. Get 65536 log-spaced
! levels between 10^-9 and 10^1 which is enough.
! See why below!
           thisoctal%levelconvergence(1:minlevel,subcell) = &
                int((max(min(log10(newFracChangePerLevel),1.d0),-9.d0) + 4.0) * 6553.6)
! Maximum change in any level in this subcell
           maxFracChange = MAXVAL(maxFracChangePerLevel(1:minlevel))
! For display purposes output all levels to see how things are converging at different levels
           temp = newFracChangePerLevel(1:minlevel-1)
! Here we count whether a level is converged at a particular tolerance in this cell so that the level
! of convergence over the entire grid can be ascertained
           do j=1, minlevel-1
                 if(newFracChangePerLevel(j) < 5.0 * tolerance) counter(1,j) = counter(1,j) + 1
                 if(newFracChangePerLevel(j) < 2.0 * tolerance) counter(2,j) = counter(2,j) + 1
                 if(newFracChangePerLevel(j) < 1.0 * tolerance) counter(3,j) = counter(3,j) + 1
                 if(newFracChangePerLevel(j) < 0.5 * tolerance) counter(4,j) = counter(4,j) + 1
           enddo
! This counts whether a cell can be called converged in its entirety (all levels)
           if(maxval(temp) .lt. 5. * tolerance) then
              counter(1,minlevel) = counter(1,minlevel) + 1
              if(maxval(temp) .lt. 2. * tolerance) then
                 counter(2,minlevel) = counter(2,minlevel) + 1
                 if(maxval(temp) .lt. tolerance) then
                    counter(3,minlevel) = counter(3,minlevel) + 1
! Last case is special. If ?? cells satisfy this criterion then the grid may also be considered
! to be converged even if the RMS is too high. This condition stops badly behaved cells from stopping the
! grid from being called converged
                 endif
              endif
           endif
! avgFracChange stores the mean average (1)/RMS average change in the grid from the previous iteration.
! It stores it as a running total which is later divded by nvoxels for the average over the grid.
           avgFracChange(1:max(2,minlevel-1),1) = avgFracChange(1:max(2,minlevel-1),1) + temp
           avgFracChange(1:max(2,minlevel-1),2) = avgFracChange(1:max(2,minlevel-1),2) + temp**2
! update maxFracChange if fractional change is great => not converged yet
           If (maxval(temp) > maxFracChange) then
              maxFracChangePerLevel(1:max(2,minlevel-1)) = newFracChangePerLevel(1:max(2,minlevel-1)) 
           endif

! Store previous levels so can use updated one in future for NgStep
           thisOctal%oldestMolecularLevel(1:minlevel,subcell) = thisOctal%oldmolecularLevel(1:minlevel,subcell)
           thisOctal%oldmolecularLevel(1:minlevel,subcell) = thisOctal%molecularLevel(1:minlevel,subcell) 
! Update working level populations 
           thisOctal%molecularLevel(1:maxlevel,subcell) = thisOctal%newmolecularLevel(1:maxlevel,subcell)

        endif
     enddo

   end subroutine swapPops

!!! READ ROUTINES !!!

subroutine calculateMoleculeSpectrum(grid, thisMolecule, dataCubeFilename, inputViewVec)

   use inputs_mod, only : itrans, nSubpixels, observerpos, rgbCube, &
        gridDistance, imageside
#ifdef USECFITSIO
   use fits_utils_mod
#endif
#ifdef MPI
   use mpi
#endif

#ifdef USECFITSIO
   integer :: status
#endif

   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   type(VECTOR) :: unitvec, observerVec, centrevec, viewvec, imagebasis(2)
   type(DATACUBE) ::  cube
   type(VECTOR), optional :: inputViewVec
   character(len=*), optional :: dataCubeFilename
   character (len=80) :: filename, message
   real(double) :: pixelWidth

   type(OCTAL), pointer :: thisoctal
   integer :: subcell

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

     if (grid%splitOverMPI) then
        if (myrankGlobal == 0) then
           writeoutput = .true.
        else
           writeoutput = .false.
        endif
     endif

     if (grid%SplitOverMPI) then
        if (myrankGlobal /= 0) then
           call deallocateUnused(grid,grid%octreeroot,everything = .true.)
           call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
           call intensityAlongRayServer(grid, thisMolecule)
           goto 666
        endif
     endif
#endif

     unitVec = VECTOR(0.d0, 0.d0, 0.d0)
     observerVec = VECTOR(0.d0, 0.d0, 0.d0)
     viewVec = VECTOR(0.d0, 0.d0, 0.d0)
     centreVec = VECTOR(0.d0, 0.d0, 0.d0)
     cube%label = " "

!Date of last major change (ISO)
     call writeinfo("molecular_mod 20091210.1608",IMPORTANT)

     call writeinfo('Setting observer parameters', TRIVIAL)
     if (.not.PRESENT(inputViewVec)) then
        call setObserverVectors(viewvec, observerVec, imagebasis)
     else
        viewVec = inputViewVec
        observerVec = dble(-griddistance*1e-10) * viewVec ! This is the EXACT position of the observer in space
        imagebasis(1) = viewvec .cross. VECTOR(0.d0, 0.d0, 1.d0)
        imagebasis(2) = imagebasis(1) .cross. viewvec
        call normalize(imagebasis(1))
        call normalize(imagebasis(2))
        pixelwidth = imageside / dble(npixels)
        imagebasis(1) = imagebasis(1) * pixelwidth ! rescale basis vectors so that natural stepsize is 1 pixelwidth 
        imagebasis(2) = imagebasis(2) * pixelwidth
     endif
     write(message,'(a,3(2x,es10.3))') "Observer Position : ", observervec
     call writeinfo(message, FORINFO)

     write(message,'(a,3(2x,es10.3))') "View Vector       : ", viewvec
     call writeinfo(message, FORINFO)

     write(message,'(a,3(2x,f10.6))')  "Image X axis      : ", imagebasis(1) / modulus(imagebasis(1))
     call writeinfo(message, FORINFO)

     write(message,'(a,3(2x,f10.6))')  "Image Y axis      : ", imagebasis(2) / modulus(imagebasis(2))
     call writeinfo(message, FORINFO)

     write(message,'(a,es10.3)')       "Pixel length      : ", modulus(imagebasis(1))
     call writeinfo(message, FORINFO)

     open(91, file="imageparams.dat",status="unknown",form="formatted")

     write(91,'(a,3(2x,es10.3))') "Observer Position : ", observervec
     write(91,'(a,3(2x,es10.3))') "View Vector       : ", viewvec
     write(91,'(a,3(2x,f10.6))')  "Image X axis      : ", imagebasis(1) / modulus(imagebasis(1))
     write(91,'(a,3(2x,f10.6))')  "Image Y axis      : ", imagebasis(2) / modulus(imagebasis(2))
     write(91,'(a,es10.3)')       "Pixel length      : ", modulus(imagebasis(1))

     close(91)

     call findSubcellTD(VECTOR(1.,1.,1.), grid%octreeroot, thisOctal, subcell)
     if(writeoutput) then
        if (.not.grid%splitOverMPI) then
           write(message, *) 'Recovering unused memory', size(thisoctal%molecularlevel(:,1)) - &
                (thismolecule%itransupper(itrans) + 1), 'levels'
           call writeinfo(message, TRIVIAL)
        endif
     endif

     call deallocateUnused(grid,grid%octreeroot,thismolecule%itransupper(itrans)+1)

! create an image using the supplied parameters (also calculate spectra)
     call writeinfo('Creating Image', TRIVIAL)
! create an image using the supplied parameters (also calculate spectra)
     if (.not.rgbCube) then
        call createimage(cube, grid, viewvec, observerVec, thismolecule, itrans, nSubpixels, imagebasis)
     else
        call createimage(cube, grid, viewvec, observerVec, thismolecule, itrans, nSubpixels, imagebasis,revVel=.true.) 
     endif

! convert intensity (ergcm-2sr-1Hz-1) to flux (Wm-2Hz-1) (so that per pixel flux is correct)     
!    call writeinfo('Converting Intensity to Flux', TRIVIAL)
!     call cubeIntensityToFlux(cube, thismolecule, itrans)

! Commented out lines are for removing background intensity - will want to do this one day

!   call TranslateCubeIntensity(cube,1.d0*Tcbr) ! 
!   call cubeIntensityToFlux(cube, thismolecule, itrans)

     if(observerpos .gt. 0) then
        if(observerpos .lt. 10) then
           write(filename, '(a,i1,i1,a,i1,a)') trim(thisMolecule%molecule),thismolecule%itransupper(itrans)-1, &
                thismolecule%itranslower(itrans)-1,'_0',observerpos,'.fits'
        else
           write(filename, '(a,i1,i1,a,i2,a)') trim(thisMolecule%molecule),thismolecule%itransupper(itrans)-1, &
                thismolecule%itranslower(itrans)-1,'_',observerpos,'.fits'
        endif
     else
        write(filename, *) 'MolRT.fits' ! can be changed to a variable name someday
        if (PRESENT(dataCubeFilename)) then
           filename = dataCubeFilename
        endif
     endif


#ifdef MPI
   if (grid%SplitOverMPI) then
         call shutdownServers()
   endif
#endif
   
#ifdef USECFITSIO
   if (writeoutput) then
      status = 0
      call writeinfo('Deleting previous Fits file', TRIVIAL)
      call deleteFitsFile (filename, status)
      call writeinfo('Writing Cube to Fits file: '//trim(filename), TRIVIAL)
      call writedatacube(cube, filename)
   endif
#endif
!     call createFluxSpectra(cube, thismolecule, itrans)

!#ifdef MPI
!     call torus_mpi_barrier
!     if(molcluster) then
!        call MPI_FINALIZE(ierr)
!     endif
!#endif

#ifdef MPI
666  continue
     if (grid%splitOverMPI) then
        if (myrankGlobal == 1) then
           writeoutput = .true.
        else
           writeoutput = .false.
        endif
     endif
#endif

    end subroutine calculateMoleculeSpectrum

!!!!READ ROUTINES!!!!

    subroutine makeImageGrid(grid, thisMolecule, iTrans, deltaV, nsubpixels, &
                             ObserverVec, viewvec, imagebasis, imagegrid, ix1, ix2)

      type(GRIDTYPE), intent(IN) :: grid
      type(MOLECULETYPE), intent(IN) :: thisMolecule
      integer, intent(IN) :: itrans
      real(double), intent(IN) :: deltaV
      integer, intent(IN) :: nsubpixels
      type(VECTOR), intent(IN) :: viewvec, ObserverVec, imagebasis(:)
      real, intent(OUT) :: imagegrid(:,:,:)
      integer, intent(in) :: ix1, ix2

      real(double) :: dnpixels ! npixels as a double, save conversion
      real(double) :: mydnpixels ! number of x pixels in this slice
      type(VECTOR) :: pixelcorner, thisPixelcorner
      integer :: subpixels
      integer :: ipixels, jpixels

      dnpixels = dble(npixels) 
      mydnpixels = dble(ix2 - ix1 + 1) 
! pixelcorner initialised to TOPLEFT      
      pixelcorner = ObserverVec - (dnpixels * 0.5d0)*(imagebasis(1) - imagebasis(2)) + imagebasis(2)

! For MPI runs move the pixel to the correct x slice, for non-MPI case ix1=1
      pixelcorner = pixelcorner + real( (ix1 - 1), kind=db) * imagebasis(1)

      if (nsubpixels .gt. 0) then ! if nsubpixels = 0 then use adaptive subpixel sampling
         subpixels = nsubpixels
      else
         subpixels = 0
      endif

      do jpixels = 1, npixels ! raster over image
         pixelcorner = pixelcorner - imagebasis(2) 
!$OMP PARALLEL default(shared), private(ipixels, thisPixelcorner)
!$OMP DO
         do ipixels = ix1, ix2
            thisPixelCorner = pixelcorner + real((ipixels-ix1),db) * imagebasis(1)
            imagegrid(ipixels,jpixels,:) = real(PixelIntensity(viewvec,thisPixelcorner,imagebasis,grid,thisMolecule,&
                 iTrans,deltaV, subpixels))
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo

    end subroutine makeImageGrid

 !!! Calculates the intensity for a square pixel of arbitrary size, position, orientation
    function PixelIntensity(viewvec,pixelcorner,&
      imagebasis,grid,thisMolecule,iTrans,deltaV,subpixels) result(out)
   
   use inputs_mod, only : tolerance, lineimage, lamline, maxrhocalc,isinlte

   integer, parameter :: maxsubpixels = 100
   
   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   integer :: itrans
   type(VECTOR) :: viewVec
   real(double) :: i0, opticaldepth

   type(VECTOR) :: imagebasis(:), pixelbasis(2), pixelcorner, rayposition
   
   integer :: subpixels, minrays
   integer :: i, iray

   real(double) :: avgIntensityNew, avgIntensityOld
   real(double) :: varIntensityNew, varIntensityOld
   real(double) :: avgNColNew, avgNColOld
   real(double) :: rtemp(2)
   real(double), save ::  r(maxsubpixels,2)
   real(double) :: nCol

   logical :: converged
   real(double) :: deltaV
   real(double) :: out(3)
   logical, save :: firsttime = .true.
   real(double) :: rhomax, i0max

!   logical :: dosink
   
!   real(double) :: sink

!$OMP THREADPRIVATE (firstTime, r)
   if(firsttime) then
      
      call sobseq(rtemp,-1)
      
      do i = 1, maxsubpixels
         call sobseq(rtemp)
         r(i,:) = dble(rtemp)
      enddo
      
      firsttime = .false.
   endif

   pixelbasis(1) = imagebasis(1)
   pixelbasis(2) = imagebasis(2)
  
   avgIntensityOld = 0.
   ncol = 0.
   opticaldepth = 0.
   varIntensityOld = 0.
   avgNColOld      = 0.0
!   sink = 0.d0
   
   converged = .false. ! failed flag
     
   if(subpixels .ne. 0) then 
      converged = .true.
      minrays = subpixels
   else
      minrays = -1
   endif
   
   iray = 1
   
   do while((.not. converged) .or. (iray .le. minrays))  

      rayposition = pixelcorner + r(iray,1) * pixelbasis(1) - r(iray,2) * pixelbasis(2) ! random position in pixel

      if(maxrhocalc) then
         call intensityalongray(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0, &
                                tau = opticaldepth, rhomax = rhomax)
         i0max = i0
         call intensityalongray(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0, &
                                tau = opticaldepth, rhomax = rhomax, i0max = i0max, nCol=nCol)
      else
         if(lineimage) then
!            if(present(sink)) then
            if (.not.grid%splitOverMPI) then
               if(isinlte) then
                  call lteintensityalongray2(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0, &
                                             tau = opticaldepth, ncol = ncol)
               else
                  call intensityalongray2(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0, &
                                          tau = opticaldepth, ncol = ncol)
!                  call intensityalongray(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0, &
!                                          tau = opticaldepth, ncol = ncol)
               endif
            else
#ifdef MPI
               call intensityalongRaySplitOverMpI(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0)
#endif
            endif
!            else
!               call intensityalongray2(rayposition,viewvec,grid,thisMolecule,itrans,deltaV,i0, &
!                                       tau = opticaldepth, ncol = ncol)
!            endif
            else
               call continuumintensityalongray(rayposition,viewvec,grid,lamline,i0,tau = opticaldepth)
            endif
         endif

!      if (debug .and. isnan(i0)) then
!         write(*,*) "Got nan", opticaldepth, rayposition
!         i0 = 0.d0
!      endif

      avgIntensityNew = ((iray - 1) * avgIntensityOld + i0) / dble(iray)
      varIntensityNew = ((iray - 1) * varIntensityOld + ((i0 - avgIntensityNew) * (i0 - avgIntensityOld))) / dble(iray)
      avgNColNew      = ((iray - 1) * avgNColOld + nCol) / dble(iray)

      if((varIntensityNew .lt. iray * (tolerance* avgIntensityNew)**2) .and. (iray .ge. 3)) then
         converged = .true.
         out(1) = avgIntensityNew
         out(2) = opticaldepth ! not averaged, just last value
!         out(2) = varintensitynew ! not averaged, just last value
         out(3) = avgNColNew

         iray = iray + 1
         exit
      elseif(iray .ge. maxsubpixels) then
         out(1) = avgIntensityNew
!         out(2) = varintensitynew ! not averaged, just last value
         out(2) = opticaldepth ! not averaged, just last value
         out(3) = avgNColNew
         converged = .false.
         exit
      else
         avgIntensityOld = avgIntensityNew
         varIntensityOld = varIntensityNew
         if(maxrhocalc) then
            out(2) = avgIntensityNew
            out(1) = rhomax / (2.d0 * mhydrogen) ! density to nh2 ! primary hdu
         else
            out(1) = avgIntensityNew
            out(2) = opticaldepth ! not averaged, just last value
         endif
         out(3) = avgNColNew

         iray = iray + 1
      endif

!      if(present(dosink)) then

   enddo

!      if(sink .eq. 1.) then
!         out(3) = sink
!      else
!         out(3) = 0.
!      endif

 end function PixelIntensity

 !!! This subroutine takes the parameters supplied to it and makes an image by calling more subroutines 

   subroutine createimage(cube, grid, viewvec, observerVec, thisMolecule, iTrans, nSubpixels, imagebasis, revVel)

     use inputs_mod, only : gridDistance, beamsize, nv, imageside, &
          maxVel, usedust, lineimage, lamline, plotlevels, debug, wanttau, dotune, h21cm
#ifdef USECFITSIO
    use inputs_mod, only : writetempfits
    use fits_utils_mod
#endif
#ifdef MPI
     use mpi_global_mod, only: nThreadsGlobal
     use mpi
#endif
     type(TELESCOPE) :: mytelescope
     type(MOLECULETYPE) :: thisMolecule
     type(GRIDTYPE) :: grid
     type(DATACUBE) :: cube
     type(VECTOR) :: viewvec, observerVec, imagebasis(2)
     real(double) :: minVel
     real(double) :: deltaV
     integer :: iTrans
     integer :: iv
     integer :: nsubpixels
     integer :: ix1, ix2
     character(len=200) :: message
     real(double) :: intensitysum, fluxsum, linefreq !,ddv
     real(double), save :: background
     real, allocatable :: temp(:,:,:)
     character(len=80) :: filename
#ifdef USECFITSIO
     integer :: status
#endif
     logical, optional, intent(in) :: revVel
     logical :: doRevVel

#ifdef MPI
     ! For MPI implementations
     integer :: ierr, n           ! error flag
     real, allocatable :: tempArray(:), tempArray2(:)
#endif

!$OMP THREADPRIVATE( background )



 ! SHOULD HAVE CASE(TELESCOPE) STATEMENT HERE

     if(.not. lineimage) then
        nv = 1
        deltaV = 0.d0
        maxvel = 0.d0
        write(message, *) "Continuum @ ", lamline / 1e4, " um" ! ang to micron
        call writeinfo(message, TRIVIAL)
     endif

     mytelescope%label = 'JCMT'
     mytelescope%diameter = 15.d2 ! diameter in cm
     mytelescope%beamsize = beamsize

     minVel = (-1.d0) * maxVel

     call writeinfo("Initialising datacube",TRIVIAL)

     if(nv .eq. 0) then
        call initCube(cube, 200, mytelescope, wantTau=wantTau) ! Make cube
     else
        call initCube(cube, nv, mytelescope, wantTau=wantTau) ! Make cube
     endif

     cube%obsDistance = gridDistance * 1d10!(in cm) Additional information that will be useful
     if (gridDistance/pctocm < 1000.0) then 
        write(message,'(a,f10.3,a)') "Observer Distance        : ",gridDistance/pctocm, " pc"
     else
        write(message,'(a,f10.3,a)') "Observer Distance        : ",gridDistance/kpctocm, " kpc"
     endif
     call writeinfo(message, TRIVIAL) 
     write(message,'(a,1pe12.3,a)') "Finest grid resolution   : ",grid%halfsmallestsubcell*2d10/autocm, " AU"
     call writeinfo(message, TRIVIAL) 
     call addSpatialAxes(cube, -imageside/2.d0, imageside/2.d0, -imageside/2.d0, imageside/2.d0, gridDistance)

     if ( present(revVel) ) then 
        doRevVel = revVel
     else
        doRevVel = .false.
     end if

     if(nv .ne. 0) then 
        if ( doRevVel ) then 
           call addvelocityAxis(cube, maxVel, minVel)
        else
! velocities in km/s from +ve (redder, away) to -ve (bluer,towards)
           call addvelocityAxis(cube, minVel, maxVel)
        end if
     else
        cube%vAxis(1) = minVel
     endif

! Divide up the image along the x axis for MPI case, otherwise work on the whole image
#ifdef MPI
     if (.not.grid%splitOverMPI) then
        ix1 = (myRankGlobal)   * (cube%nx / (nThreadsGlobal)) + 1
        ix2 = (myRankGlobal+1) * (cube%nx / (nThreadsGlobal))
        if (myRankGlobal == (nThreadsGlobal-1)) ix2 = cube%nx
        n = (npixels*npixels)
        allocate(tempArray(1:n), tempArray2(1:n))
     else
        ix1 = 1
        ix2 = npixels
     endif
#else
     ix1 = 1
     ix2 = npixels
#endif

     deltaV = minVel * 1.e5/cspeed_sgl
     
     allocate(temp(npixels,npixels,3))

     if(nv .ne. 0) then
 
        do iv = 1,nv

           deltaV = (cube%vAxis(iv)*1.e5/cSpeed_sgl) ! velocities in fraction of c

           if(writeoutput .and. dotune) then
              write(message,*) "Done ",iv," velocity"
              call tune(6, message)  ! start a stopwatch
           endif

           if(iv .eq. 1) then
              call writeinfo("Filling Octal parameters for first time",TRIVIAL)
              if (.not. h21cm ) call deallocateUnused(grid,grid%octreeroot,everything = .true.)
              call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
              call writeinfo("Done.",TRIVIAL)
           endif

           if(usedust) call adddusttoOctalParams(grid, grid%OctreeRoot, thisMolecule)

           temp = 0.d0

           call makeImageGrid(grid, thisMolecule, itrans, deltaV, nSubpixels, &
                              ObserverVec, viewvec, imagebasis, temp, ix1, ix2)

           if(plotlevels .and. myrankiszero .and. debug) then
              write(filename, '(a)') "MolRT.vtk"
              write(message, *) "Wrote 3d VTK file to ", filename
              call writeinfo(message, TRIVIAL)
              call writeVtkFile(grid, filename, "imagevtk.txt")
           endif

! Put results from makeImageGrid into data cube, performing MPI communication if required
#ifdef MPI
           if (.not.grid%splitOverMPI) then
              ! Communicate intensity
              tempArray = reshape(temp(:,:,1), (/ n /))
              call MPI_ALLREDUCE(tempArray,tempArray2,n,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 
              temp(:,:,1) = reshape(tempArray2, (/ npixels, npixels /))
              
              ! Communicate tau
              tempArray = reshape(temp(:,:,2), (/ n /))
              call MPI_ALLREDUCE(tempArray,tempArray2,n,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 
              temp(:,:,2) = reshape(tempArray2, (/ npixels, npixels /))
              
              ! Communicate column density
              tempArray = reshape(temp(:,:,3), (/ n /))
              call MPI_ALLREDUCE(tempArray,tempArray2,n,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr) 
              temp(:,:,3) = reshape(tempArray2, (/ npixels, npixels /))
           endif
#endif           

           cube%intensity(:,:,iv) = temp(:,:,1)
           if(wanttau) cube%tau(:,:,iv) = temp(:,:,2)
           cube%nCol(:,:)   = temp(:,:,3)

           if(writeoutput .and. dotune) then
              write(message,*) "Done ",iv," velocity"
              call tune(6, message)  ! stop a stopwatch
           endif

           intensitysum = sum(temp(:,:,1)) / dble(npixels**2)
           fluxsum = intensitytoflux(intensitysum, dble(imageside), dble(gridDistance), thisMolecule)

           if(iv .eq. 1) then
              if(lineimage) then
                 linefreq = thisMolecule%transfreq(itrans)
              else
                 linefreq = cspeed / (lamline * 1d-8)
              endif

              background = Bnu(linefreq, Tcbr)
              write(message, *) "Background Intensity: ",background
              call writeinfo(message, TRIVIAL)
              background = intensitytoflux(background, dble(imageside), dble(gridDistance), thisMolecule)
              write(message, *) "Background Flux: ",background
              call writeinfo(message, TRIVIAL)
              write(message, *) ""
              call writeinfo(message, TRIVIAL)
              
!             call fineGaussianWeighting(cube, cube%nx, beamsize, weight)!, NormalizeArea = .true.) 

           endif

           write(message,'(a,es11.4e1,tr3,a,f10.4,tr3,a,es12.4,a,es12.4,es12.4,es12.4)') &
                "DELTAV(v/c):",deltaV," V (km/s):",real(cube%vAxis(iv)), "Average Intensity:",intensitysum, &
                " FLUX: ", fluxsum, (fluxsum / linefreq) * 1e26, (fluxsum - background) &
                / linefreq * 1d26  
           open(10, file="tempfile.dat",status="unknown",form="formatted",position="append")
           write(10,'(es11.4e1,f8.4,es12.4,es12.4,es12.4,es12.4)') &
                real(cube%vAxis(iv)), deltaV, intensitysum, &
                fluxsum / linefreq * 1e26, (fluxsum - background) / linefreq * 1e26 
           close(10)
           call writeinfo(message,FORINFO)

#ifdef USECFITSIO
           if(writeoutput .and. writetempfits) then 
              status = 0
              write(filename,'(a,i3,a)') "MolRTtemp",iv,".fits"
              call deleteFitsFile (filename, status)
              call writedatacube(cube, filename)
              if(iv .gt. 1) then
                 write(filename,'(a,i3,a)') "MolRTtemp",iv-1,".fits"
                 call deleteFitsFile (filename, status)
              endif
           endif
#endif
        end do

     endif

#ifdef MPI
     if (.not.grid%splitOverMPI) deallocate(tempArray, tempArray2)
#endif

   end subroutine createimage

! Flux through solid angle covered by one pixel in ergs/s
 real(double) pure function IntensityToFlux(intensity,dx,distance,thismolecule) result (flux)

   use inputs_mod, only : itrans
   real(double), intent(in) :: dx, intensity, distance
   type(moleculetype), intent(in) :: thismolecule

   flux = intensity * 1d20 * (dx / distance)**2 * 1d-3 * thisMolecule%transfreq(itrans) 

 end function IntensityToFlux

 subroutine tauAlongRay(position, direction, grid, thisMolecule, deltaV, tau)

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
   
   iUpper(:)  = thisMolecule%iTransUpper(1:maxtrans)
   iLower(:)  = thisMolecule%iTransLower(1:maxtrans)

!   write(*,*) "iupper, ilower ",iupper,ilower
      
   do while(inOctal(grid%octreeRoot, currentPosition))
      
      call findSubcelllocal(currentPosition, thisOctal, subcell)
      call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)
      
      nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
!      write(*,*) "nMol ",nmol
      nLower(:)  = thisOctal%molecularLevel(iLower(:),subcell) * nMol
      nUpper(:)  = thisOctal%molecularLevel(iUpper(:),subcell) * nMol

!      write(*,*) "nupper, nlower ",nupper,nlower

      balance(:) = (nLower(:) * thisMolecule%einsteinBlu(1:maxtrans) - &
           nUpper(:) * thisMolecule%einsteinBul(1:maxtrans))
!      write(*,*) "balance ", balance   
      alphaTemp(:) = hCgsOverFourPi * balance(:) ! Equation 8
!      write(*,*) "alphatemp ",alphatemp
      thisPosition = currentPosition
      startVel = Velocity(currentposition, grid)
      
      endPosition = currentPosition + tval * direction
      endVel = Velocity(endposition, grid)
      
      Veldiff = endVel - startVel
      
      dvAcrossCell = (veldiff .dot. direction)
      dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))
      
      nTau = min(max(2, nint(dvAcrossCell * 5.d0)), 100) ! ensure good resolution / 5 chosen as its the magic number!
!      write(*,*) "dv ntau ",dvacrosscell,ntau
      CellEdgeInterval = 1.d0 / (dble(ntau) - 1.0)
      dds = tval * cellEdgeInterval
!      write(*,*) "dds ",dds
      dist = 0.d0
      
      do i = 2, nTau
         
         dist = dist + dds ! increment along distance counter
         
         thisPosition = currentPosition + dist * direction
         
         thisVel = velocity(thisPosition, grid)
!         write(*,*) "thisvel ",thisvel
         dv = deltaV - (thisVel .dot. direction)
!         write(*,*) "dv ",dv
         PhiProfVal = phiProf(dv, thisOctal%molmicroturb(subcell))
!         write(*,*) "phiprofval ",phiprofval
!         print *, "ECHO"
         alphanu(1:maxtrans,1) = phiprofval * alphaTemp(1:maxtrans) ! Equation 8
!         write(*,*) "alphanu ",alphanu(1:maxtrans,1)
!         print *, "GOLF"
         alpha(1:maxtrans) = alphanu(1:maxtrans,1) !+ alphanu(1:maxtrans,2) !!!!!!!!!!!!!!!!!!!!!!!
!         write(*,*) "alpha ",alpha
         dTau(:) = alpha(:) * dds * 1.d10 
!         write(*,*) "dtau ",dtau
! dds is interval width & optical depth, dTau = alphanu*dds - between eqs (3) and (4)
         
         tau(:) = tau(:) + dtau(:) ! contribution to optical depth from this line integral
!         write(*,*) "tau ",tau
         
      enddo
      
      currentPosition = currentPosition + (tval + 1.d-3*grid%halfSmallestSubcell) * direction
! FUDGE - make sure that new position is in new cell
   enddo

 end subroutine tauAlongRay

 subroutine intensityAlongRay(position, direction, grid, thisMolecule, iTrans, deltaV,i0, &
                              tau,tautest,rhomax, i0max, nCol, observerVelocity)

   use inputs_mod, only : useDust, h21cm, densitysubsample
     type(VECTOR) :: position, direction, dsvector
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0
     real(double), optional, intent(out) :: nCol
     type(VECTOR), optional, intent(in) ::  observerVelocity
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(VECTOR) :: currentPosition, thisPosition, endPosition
     type(VECTOR) :: startVel, endVel, thisVel, veldiff
     real(double) :: alphanu1, alphanu2, jnu, snu
     real(double) :: alpha
     real(double)  :: dv, deltaV, dVacrossCell
     integer :: i, icount
     real(double) :: tval
     integer :: nTau
     real(double) ::  OneOvernTauMinusOne, ds

     real(double) :: dTau, etaline, dustjnu
     real(double), intent(out) :: tau
     real(double),save :: BnuBckGrnd


     real(double) :: phiProfVal
     real         :: sigma_thermal
     real(double) :: deps, origdeps

     logical,save :: firsttime = .true.
     logical,save :: havebeenWarned = .false.
     character(len = 80) :: message
     logical, optional :: tautest
     logical :: dotautest

     real(double) :: dI, dIovercell, attenuateddIovercell, dtauovercell, opticaldepth, n
     real(double), optional, intent(out) :: rhomax
     real(double), optional, intent(in) :: i0max

!$OMP THREADPRIVATE( firstTime, haveBeenWarned, BnuBckGrnd)

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
        disttogrid = 0.
     else
        distToGrid = distanceToGridFromOutside(grid, position, direction) 
     endif

     if (distToGrid > 1.e29) then
        if (.not. haveBeenWarned) then
           call writewarning("ray does not intersect grid")
           write(message, *) position
           call writeinfo(message, FORINFO)
           write(message, *) direction
           call writeinfo(message, FORINFO)
           havebeenWarned = .true.

        endif
        
        i0 = 1.d-60
        tau = 1.d-60
        goto 666
     endif
     
     deps = 5.d-4 * grid%halfSmallestSubcell ! small value to add on to distance from grid to ensure we get onto grid
     origdeps = deps

     currentPosition = position + (distToGrid + deps) * direction

     i0 = 0.d0
     tau = 0.d0
     if (present(nCol)) nCol = 0.d0

     if(present(rhomax)) rhomax = 0.d0

     thisOctal => grid%octreeRoot
     icount = 0

     do while((.not. inOctal(grid%octreeRoot, currentPosition)) .and. icount .eq. 0 .and. deps .lt. 1d30)
        distToGrid = distanceToGridFromOutside(grid, currentposition, direction) 
!        deps = 10.d0 * deps
        currentPosition = currentposition + (distToGrid + deps) * direction
     enddo
   
     do while(inOctal(grid%octreeRoot, currentPosition))

        ! thisOctal%molcellparam(1,subcell) = nmol
        ! thisOctal%molcellparam(2,subcell) = nLower
        ! thisOctal%molcellparam(3,subcell) = nUpper
        ! thisOctal%molcellparam(4,subcell) = (nLower * Blu - Nupper * Bul)
        ! thisOctal%molcellparam(5,subcell) = hcgs/4pi * Aul *nUpper
        ! thisOctal%molcellparam(6,subcell) = hcgs/4pi * (nLower * Blu - Nupper * Bul)
        ! thisOctal%molcellparam(7,subcell) = kappaAbs / 1d10
        ! thisOctal%molcellparam(8,subcell) = (kappaAbs / 1d10) * Bnu

        icount = icount + 1
        call findSubcelllocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

        if(present(rhomax)) then
           rhomax = max(rhomax, thisoctal%rho(subcell))
        endif

        dtauovercell = 0.d0
        dIovercell = 0.d0
        attenuateddIovercell = 0.d0

!Get nmol from somewhere

        if(densitysubsample) then
           if ( h21cm) then
              nmol = interpolated_Density(currentposition, grid) / (thisOctal%rho(subcell))
           else
              nmol = thisoctal%molabundance(subcell) * (interpolated_Density(currentposition, grid) / &
                   (2.d0 * mhydrogen))
           end if
        else
           if ( h21cm ) then
              nMol = 1.0
           else
              nMol = thisOctal%molcellparam(1,subcell)
           endif
        end if

        etaline = nmol * thisOctal%molcellparam(5,subcell)

!Calculate absorption from dust if reqd.

        if(usedust) then
!THAW
!           alphanu2 = nmol * thisOctal%molcellparam(7,subcell)
!           dustjnu = nmol * thisOctal%molcellparam(8,subcell
           alphanu2 = thisOctal%molcellparam(7,subcell)
           dustjnu = thisOctal%molcellparam(8,subcell)
        else
           alphanu2 =  0.0
        endif

        thisPosition = currentPosition

!Get velocity at start and end of ray inside cell. Assume linear gradient.

        startVel = Velocity(currentPosition, grid, startoctal = thisoctal, subcell = subcell)
        endPosition = currentPosition + tval * direction

        endVel = Velocity(endPosition, grid, startoctal = thisoctal, subcell = subcell)

!Directional velocity gradient
        Veldiff = endVel - startVel
        dvAcrossCell = (veldiff.dot.direction)

        if ( h21cm ) then 
           ! Calculate line width in cm/s.
           sigma_thermal = real(sqrt (  (kErg * thisOctal%temperature(subcell)) / mHydrogen))
           ! Convert to Torus units (v/c)
           sigma_thermal = real(sigma_thermal / cspeed)
           dvAcrossCell = abs(dvAcrossCell / sigma_thermal)
        else
           dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))
        end if
        
!Determine number of line segments based on velocity gradient to ensure 
!adequate resolution of velocity field
! ensure good resolution / 5 chosen as its the magic number!
        if(densitysubsample) then ! should replace 5 with maxdensity/mindensity * fac
           nTau = min(max(10, nint(dvAcrossCell * 5.d0)), 1000) 
        else
           nTau = min(max(2, nint(dvAcrossCell * 5.d0)), 100)
        endif
      
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)
      
        ds = tval * OneOvernTauMinusOne
        
!Calculate column density
!Factor of 1.d10 is to convert ds to cm 
        if (present(nCol)) then 
           nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * tval * 1.d10
        end if

!        if(present(sink)) then
!           if(thisOctal%rho(subcell) .gt. 1d-14) sink = 1.
!        endif

        dsvector = ds * direction

!Repeat intensity calculation of all line segments
        do i = 2, nTau 

           thisPosition = thisPosition + dsvector

!Determine local velocity.
!Check haven't fallen off the edge of the grid. If so move way back (effectively assume 0 velocity gradient
           if(.not. inoctal(grid%octreeroot, thisposition)) thisPosition = thisPosition - 0.99d0 * dsvector
! get back in grid
           thisVel = Velocity(thisPosition, grid, startoctal = thisoctal, subcell = subcell)

           if ( present(observerVelocity) ) then 
              thisVel = thisVel - observerVelocity
           end if

           dv = (thisVel .dot. direction) - deltaV

!Calculate position in line to determine emission/absorption characteristics
           if ( h21cm ) then 
              phiprofval = gauss (sigma_thermal, real(dv) ) / thisMolecule%transfreq(1)
           else  
              phiProfval = phiProf(dv, thisOctal%molmicroturb(subcell))
           end if

!Determine jnu and alphanu           
           if(densitysubsample .and. .not. h21cm ) then
              nmol = thisoctal%molabundance(subcell) * (interpolated_Density(thisposition, grid) / &
                     (2.d0 * mhydrogen))

              if(usedust) then
!                 alphanu2 = nmol * thisOctal%molcellparam(7,subcell)
!                 dustjnu =  nmol * thisOctal%molcellparam(8,subcell)
                 alphanu2 = thisOctal%molcellparam(7,subcell)
                 dustjnu = thisOctal%molcellparam(8,subcell)
              else
                 alphanu2 = 0.0
              endif

              etaline = nmol * thisOctal%molcellparam(5,subcell)
              alphanu1 = nmol * thisOctal%molcellparam(6,subcell) * phiprofval

           else if (densitysubsample .and. h21cm) then
              nmol     = interpolated_Density(thisposition, grid) / (thisOctal%rho(subcell))
              etaline  = nmol * thisOctal%molcellparam(5,subcell)
              alphanu1 = nmol * thisOctal%molcellparam(6,subcell) * phiprofval

           else
              alphanu1 = (nmol * thisOctal%molcellparam(6,subcell)) * phiprofval
           endif


           alpha = alphanu1 + alphanu2
!Contribution to tau by line segment. dtauovercell is tracked per subcell.
           dTau = alpha * ds * 1.d10
           dtauovercell = dtauovercell + dtau

!Calculate emissivity weighted by phi
           jnu = etaLine * phiProfVal

           

           if(useDust) jnu = jnu + dustjnu

!Determine local source function
           if (alpha .ne. 0.d0) then
              snu = jnu/alpha
           else
              snu = tiny(snu)
              i0 = i0 + tiny(i0)
           endif

!Calculate and store contribution of dI to integrated intensity
!dIovercell is that which comes from the cell 
!attenuateddIovercell will be the true amount contributed to total LOS intensity

           opticaldepth = exp(-tau)
           dI = (1.d0-exp(-dtau))*snu
           dIovercell = dIovercell + dI

           tau = tau + dtau

           dI = opticaldepth * dI
           attenuateddIovercell = attenuateddIovercell + dI

!Check still worth adding dI in optically thick case. If so add dI to I
!This should really exit the entire ray but doesn't so is commented out at the moment
!           if(dI .gt. i0 * 1d-10 .and. tau .lt. 100) then

              i0 = i0 + dI
!           else
!              i0 = i0
!           endif
              
           enddo

!i0max is a parameter for calculating the depth to which information is being received
!probably strongly linked to tau?! This method requires 2 passes. One to determine i0max
!the second to determine when 99% of total intensity has been reached.
        if(present(i0max) .and. i0 .gt. 0.99d0 * i0max) then 
           exit
        endif

!These lines store various parameters of interest regarding the contribution of a cell

!None of this is normalised by path length
!I presume this is all integrated over velocity as well at the moment so it would need to be taken into account
!depending on what one wanted to do with the data.

!This is a counter as they're all stored as running averages. I hope this is initialised to 0 somewhere!
        n = thisoctal%newmolecularlevel(4,subcell)

!Store average dtau inside cell
        thisoctal%newmolecularlevel(5,subcell) = &
        (n * thisoctal%newmolecularlevel(5,subcell) + dtauovercell) / (n + 1.d0)

!Store average intensity contribution inside cell
!        thisoctal%newmolecularlevel(1,subcell) = 
!        (n * thisoctal%newmolecularlevel(1,subcell) + dIovercell) / (n + 1.d0)
        thisoctal%newmolecularlevel(1,subcell) = thisoctal%newmolecularlevel(1,subcell) +  attenuateddiovercell

!Store average attenuated intensity contribution inside cell
        thisoctal%newmolecularlevel(2,subcell) = &
             (n * thisoctal%newmolecularlevel(2,subcell) + attenuateddIoverCell) / (n + 1.d0)

!Store average i0 to date. It felt like this should be useful. I cannot say why...
        thisoctal%newmolecularlevel(3,subcell) = &
             (n * thisoctal%newmolecularlevel(3,subcell) + i0) / (n + 1.d0)

!Increment counter
        thisoctal%newmolecularlevel(4,subcell) = thisoctal%newmolecularlevel(4,subcell) + 1.d0

        currentPosition = currentPosition + (tval + origdeps) * direction
     enddo
         
666  continue

     if(present(rhomax) .and. .not. present(i0max)) then
        return
     else
        i0 = i0 + bnuBckGrnd * exp(-tau) ! from far side
     endif
     
   end subroutine intensityAlongRay

   subroutine continuumIntensityAlongRay(position, direction, grid, lambda, i0, tau, tautest)

     use inputs_mod, only : densitysubsample

     type(VECTOR) :: position, direction
     real(double) :: disttoGrid
     type(GRIDTYPE) :: Grid
     real(double), intent(out) :: i0
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(VECTOR) :: currentPosition
     real(double) :: tval, ds
     real :: lambda
     real(double) :: fac
     real(double) :: alpha, jnu, snu
     real(double),save :: bnubackground
     real(double) :: dTau
     real(double), intent(out), optional :: tau
     real(double) :: deps, origdeps
     logical,save :: havebeenWarned = .false.
     character(len = 80) :: message
     integer :: icount
     logical, optional :: tautest
     logical :: dotautest

     logical,save :: firsttime = .true.

     real(double) :: dI
     !$OMP THREADPRIVATE( firstTime, haveBeenWarned, bnuBackground)
     if(firsttime) then
        bnubackground = bnu(cspeed / (lambda * 1d-8), dble(tcbr))
        firsttime = .false.
     endif

     if(present(tautest)) then
        dotautest = tautest
     else
        dotautest = .false.
     endif
     
     if(inOctal(grid%octreeRoot, Position)) then
        disttogrid = 0.
     else
        distToGrid = distanceToGridFromOutside(grid, position, direction) 
     endif

     if (distToGrid > 1.e29) then
        if (.not. haveBeenWarned) then
           call writewarning("ray does not intersect grid")
           write(message, *) position
           call writeinfo(message, FORINFO)
           write(message, *) direction
           call writeinfo(message, FORINFO)
           havebeenWarned = .true.
        endif
        
        i0 = 1.d-60
        tau = 1.d-60
        goto 666
     endif
     
     deps = 5.d-4 * grid%halfSmallestSubcell ! small value to add on to distance from grid to ensure we get onto grid
     origdeps = deps

     currentPosition = position + (distToGrid + deps) * direction

     i0 = 0.d0
     tau = 0.d0

     thisOctal => grid%octreeRoot
     icount = 0

     do while((.not. inOctal(grid%octreeRoot, currentPosition)) .and. icount .eq. 0 .and. deps .lt. 1d30)
        deps = 10.d0 * deps
        currentPosition = position + (distToGrid + deps) * direction
     enddo
     
     do while(inOctal(grid%octreeRoot, currentPosition))
        icount = icount + 1
        
        call findSubcelllocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)
        
        if(densitysubsample) then
           fac = (thisoctal%molabundance(subcell) * interpolated_Density(currentposition, grid)) / &
                 ((thisOctal%molcellparam(1,subcell))*(2.d0 * mhydrogen))
        else
           fac = 1.d0
        endif

        alpha = thisOctal%molcellparam(7,subcell)
        jnu = thisOctal%molcellparam(8,subcell)

        if (alpha /= 0.d0) then
           snu = jnu/alpha
        else
           snu = 0.
        endif

        ds = tval
           
        dTau = fac * alpha * ds * 1.d10
        
        dI = exp(-tau) * (1.d0-exp(-dtau))*snu
        tau = tau + dtau
        
        i0 = i0 + dI

        currentPosition = currentPosition + (tval + origdeps) * direction
     enddo

666     continue

     i0 = i0 + bnubackground * exp(-tau) ! from far side

   end subroutine continuumintensityAlongRay

   recursive subroutine calculateOctalParams(grid, thisOctal, thisMolecule)

     use inputs_mod, only : iTrans, h21cm, lowmemory, doCOchemistry, x_D, noturb

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i
     integer :: iupper, ilower
     real(double) :: nmol, nlower, nupper, etaline

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call calculateOctalParams(grid, child, thisMolecule)
                 exit
              end if
           end do
        else

           if (grid%splitOverMPI.and.(.not.octalOnThread(thisOctal, subcell, myrankGlobal))) cycle
           
!           if(grid%geometry .eq. "iras04158") then
              
!              thisOctal%velocity(subcell) = keplerianVelocity(subcellcentre(thisOctal,subcell), grid)
!              CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
              
!              out = 'abundance'
!              rvec = subcellCentre(thisOctal,subcell)
!              thisOctal%molabundance(subcell) = readparameterfrom2dmap(rvec,out,.true.)
!                 thisOctal%molabundance(subcell) = 1d-4
!                 thisOctal%nh2(subcell) = thisOctal%nh2(subcell) * 1.05
!                 thisOctal%rho(subcell) = thisOctal%rho(subcell) * 1.05
!           endif             

           if ( h21cm ) then 

              if (.not.associated(thisOctal%molcellparam)) then
                 allocate(thisOctal%molcellparam(1:thisOctal%maxChildren,8))
              endif

              thisOctal%molcellparam(1,subcell) = 0.0
              thisOctal%molcellparam(2,subcell) = 0.0
              thisOctal%molcellparam(3,subcell) = 0.0
              thisOctal%molcellparam(4,subcell) = 0.0
              thisOctal%molcellparam(5,subcell) = thisOctal%etaLine(subcell)
              thisOctal%molcellparam(6,subcell) = thisOctal%chiLine(subcell)
              thisOctal%molcellparam(7,subcell) = 0.0
              thisOctal%molcellparam(8,subcell) = 0.0 

              allocate(thisOctal%molmicroturb(1:thisOctal%maxChildren))
              thisOctal%molmicroturb(:) = 1.0

           else

              if (.not.associated(thisOctal%molmicroturb)) then
                 allocate(thisOctal%molmicroturb(1:thisOctal%maxChildren))
              endif

              if(.not. lowmemory) then
                 if (.not.associated(thisOctal%molcellparam)) then
                    allocate(thisOctal%molcellparam(8,1:thisOctal%maxChildren))
                 endif
              endif


                 ! POST-PROCESS MOLCLUSTER
                 if(molcluster) then
!                    if(.not. thisoctal%inflow(isubcell)) then
                    if(thisoctal%rho(subcell) .le. 1e-36) then
                    
                       thisoctal%rho(subcell) = 0.d0
                       thisoctal%nh2(subcell) = 0.d0
                       thisoctal%cornervelocity(:) = VECTOR(0.d0,0.d0,0.d0)
                       if(associated(thisoctal%cornerrho)) thisoctal%cornerrho(:) = 0.d0
                    endif
                 endif
                 
                 if(.not. lowmemory) then
                    if(noturb) then
                       thisOctal%microturb(subcell) = max(2d-7,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / &
                         (thisMolecule%molecularWeight * amu))) / (cspeed * 1d-5))
                    endif

                    thisOctal%molmicroturb(subcell) = 1.d0 / thisOctal%microturb(subcell)
                    
                    if(doCOchemistry) then
                       if(thisOctal%nh2(subcell) .gt. 3e4 .and. & !drop fraction
                            thisOctal%temperature(subcell) .lt. 30.) thisOctal%molAbundance(subcell) = x_D
                    endif

                    thisOctal%molcellparam(1,subcell) = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
                    nMol = thisOctal%molcellparam(1,subcell)
                    
                    iUpper = thisMolecule%iTransUpper(iTrans)
                    iLower = thisMolecule%iTransLower(iTrans)
                    
                    thisOctal%molcellparam(2,subcell) = thisOctal%molecularLevel(iLower,subcell)! * nMol
                    thisOctal%molcellparam(3,subcell) = thisOctal%molecularLevel(iUpper,subcell)! * nMol
                    
                    nLower = thisOctal%molcellparam(2,subcell)
                    nUpper = thisOctal%molcellparam(3,subcell)
                    thisOctal%molcellparam(4,subcell) = nLower * thisMolecule%einsteinBlu(iTrans) &
                         - nUpper * thisMolecule%einsteinBul(iTrans)
                    
                    etaLine = hCgsOverFourPi * thisMolecule%einsteinA(iTrans)
                    
                    thisOctal%molcellparam(5,subcell) = etaLine * nUpper
                    thisOctal%molcellparam(6,subcell) = hCgsOverFourPi * thisOctal%molcellparam(4,subcell)! balance

                 endif
              endif
           endif
        
        if( .not. associated(thisOctal%newmolecularlevel) ) then
           allocate(thisoctal%newmolecularlevel(5, thisoctal%maxchildren))
           thisOctal%newmolecularlevel(:,:) = 0.0d0
        end if
     end do
     
   end subroutine calculateOctalParams

   recursive subroutine addDustToOctalParams(grid, thisOctal, thisMolecule)

     use inputs_mod, only : iTrans, lineimage, lamline

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i
     integer, save :: ilambda = 1
     real(double) :: kappaAbs
     real :: lambda
     !$OMP THREADPRIVATE (ilambda)
     kappaAbs = 0.d0

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call addDustToOctalParams(grid, child, thisMolecule)
                 exit
              end if
           end do
        else

           if(lineimage) then
              lambda = real((cspeed_sgl * 1e8) / thisMolecule%transfreq(iTrans)) ! cm to Angstroms
           else
              lambda = lamline
           endif

!           if(onetime) then
!              allocate(xArray(1:nLambda))
!              logLamStart = log10(lamstart)
!              logLamEnd   = log10(lamend)
!              do i = 1, nLambda
!                 xArray(i) = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
!                 xArray(i) = 10.**xArray(i)
!              enddo
!
!              call createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, SIZE(miePhaseGlobal, 3))
!              deallocate(miephase)
!              deallocate(xarray)
!              onetime = .false.
!           endif
           call locate(grid%lamArray, size(grid%lamArray), lambda, ilambda)
           call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = lambda, kappaAbs = kappaAbs)
           thisOctal%molcellparam(7,subcell) = kappaAbs * 1.d-10
           if(associated(thisoctal%temperaturedust)) then
              thisOctal%molcellparam(8,subcell) = thisOctal%molcellparam(7,subcell) * bnu(cspeed/(lambda * 1d-8), &
                                                  dble(thisOctal%temperaturedust(subcell)))
           else
              thisOctal%molcellparam(8,subcell) = thisOctal%molcellparam(7,subcell) * bnu(cspeed/(lambda * 1d-8), &
                                                  dble(thisOctal%temperature(subcell)))
           endif
        endif
     enddo
     
   end subroutine addDustToOctalParams

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

     use inputs_mod, only : rinner, router, amr1d, getdepartcoeffs, outputconvergence
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: r, ang
     real(double), save :: logrinner, logrouter, OneOverNradiusMinusOne, OneOverNangle
     integer :: i
     real(double) :: pops(60), dc(5), fracChange(60), tauarray(60)!, convtestarray(:,:,:), tauarray(40)
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     integer :: j
     real(double) :: x, z
     type(VECTOR) :: posVec
     integer :: iter
     integer, save :: nradius, nangle
     character(len=30) :: resultfile, resultfile2
     logical, save :: firsttime = .true.

     real(double) :: molarray(minlevel,3),rarray(3),molout(minlevel)
     real(double) :: rminus, rcen, rplus, tval
     type(VECTOR) :: dir,pplus,pminus
     integer :: ilevel!,itrans

     !$OMP THREADPRIVATE (firstTime, logrinner, logrouter, oneovernradiusminusone, oneOverNangle, nradius, nangle)
     iter = grand_iter
     if (.not.writeoutput) goto 666

     write(resultfile,'(a,I7.7)') "results.", nRay
     write(resultfile2,'(a,I2.2)') "fracChangeGraph.", iter

     open(31,file=resultfile,status="unknown",form="formatted")
     open(310,file='results.dat',status="unknown",form="formatted")
     open(32,file=resultfile2,status="unknown",form="formatted")
     if(getdepartcoeffs) open(34,file='departcoeffs.dat',status="unknown",form="formatted")

     call taualongray(VECTOR(1.d-10,1.d-10,1.d-10), VECTOR(1.d0, -1.d-20, -1.d-20),&
          grid, thisMolecule, 0.d0, tauarray(1:maxtrans))

     if(firsttime) then

        open(33,file="tau.dat",status="unknown",form="formatted")

        nradius = 100! number of radii used to sample distribution

        if(amr1d) then
           nangle = 1
        else
           nangle = 360 ! number of rays used to sample distribution
        endif
        
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
        dc = 0.d0

        do j = 1 , nangle
           ang = dble(j-1) * OneOverNangle
           ang = ang * twopi 
           z = r * cos(ang)
           x = r * sin(ang)
           posVec = VECTOR(x-1.d-20, 0.d0, z-1.d-20) ! get put in octal on grid
           dir = VECTOR(x/r,0,z/r)

           if(inoctal(grid%octreeroot, posvec)) then
              call findSubcellLocal(posVec, thisOctal,subcell) 
              rcen = modulus(subcellcentre(thisoctal,subcell))
              molarray(:,2) = thisOctal%molecularLevel(1:minlevel,subcell)

              call distanceToCellBoundary(grid, posVec, dir, tVal, sOctal=thisOctal)
              pplus = posvec + (tval*1.01) *dir
              call findSubcellLocal(pplus, thisOctal,subcell) 
              rplus = modulus(subcellcentre(thisoctal,subcell))
              molarray(:,3) = thisOctal%molecularLevel(1:minlevel,subcell)

              call distanceToCellBoundary(grid, posVec, (-1.d0)*dir, tVal, sOctal=thisOctal)
              pminus = posvec + tval*(-1.01d0)*dir
              call findSubcellLocal(pminus, thisOctal,subcell) 
              rminus = modulus(subcellcentre(thisoctal,subcell))
              molarray(:,1) = thisOctal%molecularLevel(1:minlevel,subcell)

!              rarray = (/rminus*0.9999d0,rcen,rplus*1.0001d0/)
              rarray = (/rminus,rcen,rplus/)


              do ilevel=1,minlevel
                 if(molarray(ilevel,1) .ne. molarray(ilevel,2) .and. &
                      molarray(ilevel,2) .ne. molarray(ilevel,3) .and. &
                      molarray(ilevel,1) .ne. molarray(ilevel,3) .and. &
                    modulus(pplus) < router .and. modulus(pminus) > rinner) then
                    molout(ilevel) = general_quadint(rarray,molarray(ilevel,:),r)
                 else
                    molout(ilevel) = molarray(ilevel,2)
                 endif
              enddo
             
              pops(1:minlevel) = pops(1:minlevel) + molout(1:minlevel)

              if(getdepartcoeffs) &
                   dc = dc + (thisOctal%molecularLevel(1:5,subcell) * thisOctal%departcoeff(1:5,subcell))
              fracChange(1:minlevel) = fracChange(1:minlevel) + abs(((thisOctal%molecularLevel(1:minlevel,subcell) - &
                   thisOctal%oldmolecularLevel(1:minlevel,subcell)) / thisOctal%molecularlevel(1:minlevel,subcell)))
           else
              call torus_abort("something went wrong in dumpresults")
           endif
        enddo

        pops = pops * OneOverNangle ! normalised level population at the 20 random positions 
        fracChange = fracChange * OneOverNangle
        dc = dc * OneOverNangle
        if(outputconvergence) then
! results.dat
           write(31,'(es12.5e2,4x,14(es14.6e2,tr2))') r*1.e10, pops(1:min(14,maxlevel))
! results.niter
           write(310,'(es12.5e2,4x,14(es14.6e2,tr2))') r*1.e10, pops(1:min(14,maxlevel))
! fracChangeGraph
           write(32,'(es12.5e2,4x,14(f8.5,tr2))') r*1.e10, fracChange(1:min(14,minlevel))
! departcoeffs.dat
           if(getdepartcoeffs) write(34,'(es12.5e2,4x,5(es14.6e2,tr2))') r*1.e10, dc(1:5)
        endif

     enddo

!     write(33,'(i2,tr3,12(f11.6,tr3))') grand_iter, tauarray(1:min(12,mintrans))
     
     close(31)
     close(310)
     close(32)
     if(getdepartcoeffs) close(34)
666 continue
   end subroutine dumpResults

!  Simply find 1st level populated at greater than 10^-8 in LTE
! Also only works for linear molecules
   recursive subroutine  findmaxlevel(grid, thisOctal, thisMolecule, maxinterestinglevel, nlevels, nVoxel, lte)
     use inputs_mod, only: setmaxleveltemp

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i, nlevels
     real(double) :: mollevels(nlevels)
     integer :: maxinterestinglevel, nVoxel
     integer, save :: icount = 0
     real(double) :: maxtemp = 0.d0, minTemp = 1.d30
     character(len = 80) :: message
     logical :: lte
     logical, save :: done = .false.

     !$OMP THREADPRIVATE (done, icount)
     mollevels = 0.d0
         
     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call findmaxlevel(grid, child, thisMolecule, maxinterestinglevel, nlevels, nVoxel, lte)
                 exit
              end if
           end do
        else

           if (lte) then
              maxtemp = max(real(maxtemp), thisOctal%temperature(subcell))     
              mintemp = min(real(mintemp), thisOctal%temperature(subcell))     
           else
              mollevels(:) = max(mollevels(:), thisoctal%molecularlevel(1:nlevels,subcell))
           endif

           icount = icount + 1
        endif
     end do

     if(icount .eq. nVoxel .and. .not. done) then

        write(message, *)  "Maximum Global Temperature", maxtemp
        call writeinfo(message, FORINFO)
        write(message, *)  "Minumum Global Temperature", mintemp
        call writeinfo(message, FORINFO)

! Apply the maximum temperature specified in the parameters file.
        if (maxtemp > setmaxleveltemp ) then
           write(message,*) "Resetting maximum temperature to ", setmaxleveltemp
           call writeinfo(message, FORINFO)
           maxtemp = setmaxleveltemp
        endif

        if (lte) call LTEpops(thisMolecule, maxtemp, mollevels(1:thisMolecule%nlevels))
                
        do i = thisMolecule%nlevels, 1, -1
           maxinterestinglevel = i
           if(mollevels(i) .gt. 1d-8) exit
        enddo
        done = .true.
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

#ifdef MPI 
       subroutine packMoleLevel(octalArray, nTemps, tArray, iLevel,ioctal_beg,ioctal_end)
!         use inputs_mod,only : gettau

         integer :: ioctal_beg, ioctal_end
         type(OCTALWRAPPER) :: octalArray(:)
         integer :: nTemps
         real(double) :: tArray(:,:)
         integer :: iOctal, iSubcell
         integer :: iLevel
         type(OCTAL), pointer :: thisOctal

        !
        ! Update the edens values of grid computed by all processors.
        !
        nTemps = 0
        do iOctal = 1, size(octalArray)

           thisOctal => octalArray(iOctal)%content

           do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 
                 if ((iOctal >= iOctal_Beg).and.(ioctal <= ioctal_end)) then
                    tArray(1,ntemps) = thisOctal%newMolecularLevel(iLevel,isubcell)
       
!                    if(gettau.and.(iLevel .le. mintrans)) then
!                       tArray(2,ntemps) = thisOctal%tau(ilevel,isubcell)
!                    else
!                       tArray(2,ntemps) = 0.d0
!                    endif
                    
                    if(iLevel .eq. maxlevel) then
                       tArray(3,ntemps) = thisOctal%convergence(isubcell)
                    else
                       tArray(3,ntemps) = 0.d0
                    endif

                 endif
              endif
           end do
        end do
      end subroutine packMoleLevel

       subroutine packJnuTrans(octalArray, nTemps, tArray, iTrans, ioctal_beg,ioctal_end)

         integer :: ioctal_beg, ioctal_end
         type(OCTALWRAPPER) :: octalArray(:)
         integer :: nTemps
         real(double) :: tArray(:)
         integer :: iOctal, iSubcell
         integer :: iTrans
         type(OCTAL), pointer :: thisOctal

        !
        ! Update the edens values of grid computed by all processors.
        !
        nTemps = 0
        do iOctal = 1, size(octalArray)

           thisOctal => octalArray(iOctal)%content

           do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 nTemps = nTemps + 1
                 
                 if ((iOctal >= iOctal_Beg).and.(ioctal <= ioctal_end)) then
                    tArray(ntemps) = thisOctal%jnu(iTrans,isubcell)
       
                 endif
              endif
           end do
        end do
      end subroutine packJnuTrans

      subroutine unpackMoleLevel(octalArray, tArray, iLevel)
!        use inputs_mod, only : gettau

        type(OCTALWRAPPER) :: octalArray(:)
        real(double) :: tArray(:,:)
        integer :: iOctal, iSubcell
        integer :: iLevel, ntemp
        type(OCTAL), pointer :: thisOctal
        
        !
        ! Update the edens values of grid computed by all processors.
        !
        ntemp = 0
        do iOctal = 1, SIZE(octalArray)
           
           thisOctal => octalArray(iOctal)%content
           
           do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 ntemp = ntemp + 1
                 thisOctal%newMolecularLevel(iLevel,isubcell) = tArray(1,ntemp) 
!                 if(gettau.and.(iLevel .le. mintrans)) then
!                    thisOctal%tau(ilevel, isubcell) = tArray(2,ntemp)
!                 endif
                 if(iLevel .eq. maxlevel) then
                    thisOctal%convergence(isubcell) = real(tArray(3,ntemp))
                 endif
              endif
           end do
        end do
      end subroutine unpackMoleLevel

      subroutine unpackJnuTrans(octalArray, tArray, iTrans)

        type(OCTALWRAPPER) :: octalArray(:)
        real(double) :: tArray(:)
        integer :: iOctal, iSubcell
        integer :: iTrans, ntemp
        type(OCTAL), pointer :: thisOctal
        
        !
        ! Update the edens values of grid computed by all processors.
        !
        ntemp = 0
        do iOctal = 1, SIZE(octalArray)
           
           thisOctal => octalArray(iOctal)%content
           
           do iSubcell = 1, thisOctal%maxChildren
              if (.not.thisOctal%hasChild(iSubcell)) then
                 ntemp = ntemp + 1
                 thisOctal%jnu(itrans,isubcell) = tArray(ntemp) 
              endif
           end do
        end do
      end subroutine unpackJnuTrans

#endif

      subroutine setObserverVectors(viewvec, observerVec, imagebasis)

        use inputs_mod, only : imageside
        use inputs_mod, only : griddistance 
        use inputs_mod, only : observerpos ! line number of observer position in observerfile 
        use inputs_mod, only : centrevecX, centrevecY, centrevecZ
        use inputs_mod, only : rotateViewAboutX, rotateviewAboutY, rotateviewAboutZ

        logical :: paraxial 
        real(double) :: pixelwidth, theta
        integer :: iread
        type(VECTOR) :: centreVec, unitvec, viewvecprime, rotationaxis

        type(VECTOR), intent(OUT) :: viewvec, observerVec, imagebasis(2)

        paraxial = .true.
        CentreVec%x = centreVecX ! 
        CentreVec%y = centreVecY ! These coordinates will be at the centre of the projected image
        CentreVec%z = centreVecZ !         

        if(observerpos .gt. 0) then
           open(51, file = 'observerpos.dat', status="old", form = "formatted") 
           do iread = 1, observerpos
              read(51,*) observervec%x, observervec%y, observervec%z
           enddo
           viewvec = (-1.d0) * observervec
           imagebasis(1) = viewvec .cross. zhat
           imagebasis(2) = imagebasis(1) .cross. viewvec
        else
! Define a vector that points *from* a specific observer that may be rotated wlog
           viewVec = VECTOR(0.d0,1.d0,0.d0) 
           imagebasis(1) = VECTOR(1.d0,0.d0,0.d0) ! Define a horizontal basis vector for the image
           imagebasis(2) = VECTOR(0.d0,0.d0,1.d0) ! Define a vertical basis vector for the image
! This line controls inclination in a 2D geometry (e.g. IRAS04158, shakara) - ANTI-CLOCKWISE ROTATION           
           viewvec = rotateX(viewvec, -rotateViewAboutX * degtorad) 
! Image bases are necessarily orthonormal and linearly independent
           imagebasis(1) = rotateX(imagebasis(1), -rotateViewAboutX * degtorad) 
           imagebasis(2) = rotateX(imagebasis(2), -rotateViewAboutX * degtorad)
           
           viewvec = rotateY(viewvec, -rotateViewAboutY * degtorad) 
           imagebasis(1) = rotateY(imagebasis(1), -rotateViewAboutY * degtorad)
           imagebasis(2) = rotateY(imagebasis(2), -rotateViewAboutY * degtorad)
! This line varies 'longitude'. This should only affect 3D geometries
           viewvec = rotateZ(viewvec, -rotateViewAboutZ * degtorad)
! Image bases are necessarily orthonormal and linearly independent
           imagebasis(1) = rotateZ(imagebasis(1), -rotateViewAboutZ * degtorad)
           imagebasis(2) = rotateZ(imagebasis(2), -rotateViewAboutZ * degtorad)
        endif

        call normalize(viewvec) 
        call normalize(imagebasis(1)) ! These lines should be unnecessary. Remove them if you feel lucky...
        call normalize(imagebasis(2))

        pixelwidth = imageside / dble(npixels)
        imagebasis(1) = imagebasis(1) * pixelwidth ! rescale basis vectors so that natural stepsize is 1 pixelwidth 
        imagebasis(2) = imagebasis(2) * pixelwidth
! This is a vector parallel to the line of sight *from* the observer to the grid centre
        unitVec = (-1.d0) * viewvec 
        observerVec = dble(griddistance*1e-10) * unitvec ! This is the EXACT position of the observer in space

        if(paraxial) then
           ! Paraxial case 
           
           ! viewvec et al remain unchanged.
           ! observerVec is translated such that the viewvec from the observerVec passes through CentreVec
           
           ObserverVec = ObserverVec + CentreVec        
        else

           ! Near-field case

           viewVecprime = CentreVec - ObserverVec
           call normalize(viewVecPrime)

           theta = acos(viewvec .dot. viewVecprime)
           rotationaxis = viewvec .cross. viewVecprime
           call normalize(rotationaxis)

           viewvec = arbitraryRotate(viewvec,theta,rotationaxis)
           imagebasis(1) = arbitraryRotate(imagebasis(1),theta,rotationaxis)
           imagebasis(2) = arbitraryRotate(imagebasis(2),theta,rotationaxis)

        endif

      end subroutine setObserverVectors
      
!!$!!!subroutine testOpticalDepth
!!$      subroutine testOpticalDepth(grid,thisMolecule)
!!$
!!$        use inputs_mod, only : lamstart, lamend, rinner, router, debug, useDust
!!$        type(GRIDTYPE) :: grid
!!$        type(MOLECULETYPE) :: thisMolecule
!!$        type(octal), pointer   :: thisOctal
!!$        type(VECTOR) :: currentposition(3), posvec, viewvec, unitvec, centrevec
!!$        integer :: subcell, i, itrans
!!$        integer :: ilamb, nlamb
!!$        real(double) :: xmidplane, gridsize
!!$        real :: lamb
!!$        real(double) :: tau, dummy, kappaAbs, kappaSca, i0
!!$        character(len=50) :: message
!!$        itrans = 0
!!$        centreVec= VECTOR(0.d0, 0.d0, 0.d0)
!!$        kappaAbs = 0.d0; kappaSca = 0.d0
!!$
!!$        call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
!!$        call addDustToOctalParams(grid, grid%OctreeRoot, thisMolecule)
!!$  
!!$        write(message,*) "Angular dependence"
!!$        call writeinfo(message, FORINFO)
!!$        
!!$        gridsize = grid%octreeroot%subcellsize
!!$        nlamb = 500
!!$        
!!$        do i = 1, 90
!!$           call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
!!$     
!!$           write(message,*) i, tau
!!$           call writeinfo(message, FORINFO)
!!$        enddo
!!$  
!!$        write(message,*) "Tau from above"
!!$        call writeinfo(message, FORINFO)
!!$        
!!$        do i = 1, 100
!!$     
!!$           unitvec = VECTOR(1.d-20,1.d-20,1.d0)
!!$           posvec = VECTOR(real(i) * grid%octreeroot%subcellsize * 0.01,0d7,2.d0 * grid%octreeroot%subcellsize)
!!$           viewvec = (-1.d0) * unitvec
!!$           
!!$           call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
!!$           
!!$           write(message,*) i, tau
!!$           call writeinfo(message, FORINFO)
!!$        enddo
!!$        
!!$        write(message,*) "Tau from the side"
!!$        call writeinfo(message, FORINFO)
!!$        
!!$        do i = 1, 100
!!$           unitvec = VECTOR(1.d0,0.d0,0.d0)
!!$           posvec = VECTOR(2.d0 * grid%octreeroot%subcellsize,0.d0, real(i) * grid%octreeroot%subcellsize * 0.01)
!!$           viewvec = (-1.d0) * unitvec
!!$           call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
!!$           
!!$           write(message,*) i, tau
!!$           call writeinfo(message, FORINFO)
!!$        enddo
!!$        
!!$        call intensityAlongRay(VECTOR(0.d0,0.d0,0.d0),VECTOR(1d-20,1d-20,1.d0), grid, thisMolecule, 1, 0.d0, &
!!$                               dummy, tau, .true.)
!!$        
!!$        xmidplane = rinner + (router-rinner) * 0.5
!!$  
!!$        if(writeoutput) write(*,*) "Midplane tau!"
!!$        if (useDust) then
!!$           open(50, file = 'dustcheck.dat', status="unknown", form = "formatted") 
!!$           currentposition(1) = VECTOR(xmidplane,0.d0,0.d0)
!!$           
!!$           call findSubcellTD(currentPosition(1), grid%octreeroot, thisOctal, subcell)
!!$           
!!$           do i = 0, nlamb
!!$              
!!$              lamb = real((10.0**(dble(i) * log10(lamend/lamstart) / 500.0)) * lamstart)
!!$              
!!$              call continuumIntensityAlongRay(VECTOR(1.d-10,1d-10,1d-10),VECTOR(1.0,1d-20,1.d-20), &
!!$                                              grid, lamb, dummy, tau, .true.)
!!$              call locate(grid%lamArray, size(grid%lamArray), real(lamb), ilamb)
!!$              call returnKappa(grid, thisOctal, subcell, ilambda = ilamb, lambda = real(lamb),&
!!$                   kappaAbs = kappaAbs, kappaSca = kappaSca)
!!$              write(50, '(f8.4,tr3,f10.4,tr3,f10.4,tr3,f10.4)') lamb * 1d-4, tau, &
!!$                   kappaAbs*1e-10 / thisOctal%rho(subcell), (kappaAbs+kappaSca)*1e-10 / thisOctal%rho(subcell)
!!$           enddo
!!$           
!!$           close(50)
!!$        endif
!!$        
!!$        if(debug) then
!!$           
!!$           currentposition(1) =  VECTOR(xmidplane,0.d0,0.d0)
!!$           currentposition(2) =  VECTOR(xmidplane,0.d0,xmidplane)
!!$           currentposition(3) = VECTOR(xmidplane/sqrt(2.),xmidplane/sqrt(2.),2. * xmidplane)
!!$           
!!$           do i = 1,3
!!$              
!!$              call findSubcellTD(currentPosition(i), grid%octreeroot, thisOctal, subcell)
!!$              currentposition(i) = subcellcentre(thisOctal, subcell)
!!$              
!!$              write(message,*) "currentposition", currentposition(i)
!!$              call writeinfo(message, FORINFO)
!!$              write (message,*) "r", sqrt(currentposition(i)%x**2+currentposition(i)%y**2)
!!$              call writeinfo(message, FORINFO)
!!$              write(message,*) "cell velocity",modulus(thisOctal%velocity(subcell)) * cspeed / 1d5
!!$              call writeinfo(message, FORINFO)
!!$              write(message,*) "calc V ",modulus(keplerianVelocity(currentposition(i))) * cspeed / 1d5
!!$              call writeinfo(message, FORINFO)
!!$              
!!$           enddo
!!$   
!!$        endif
!!$        
!!$        call continuumIntensityAlongRay(VECTOR(-1.d10,1.d-10,1.d-10),VECTOR(1.d0,-1d-10,-1d-10), &
!!$                                        grid, 1e4, dummy, tau, .true.)
!!$        write(message,*) "TAU @ 1micron", tau
!!$        call writeinfo(message, FORINFO)
!!$        
!!$      end subroutine testOpticalDepth

      subroutine plotdiscValues(grid, thisMolecule)
        
        type(GRIDTYPE) :: grid
        type(MOLECULETYPE) :: thisMolecule
  real(double) :: mean(6)
  mean = 0.d0
  call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
  call addDustToOctalParams(grid, grid%OctreeRoot, thisMolecule)
  
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 2)
  
  call readAMRgrid("molecular_tmp.grid",.false.,grid)
  
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 0)
  
  call readAMRgrid("molecular_tmp.grid",.false.,grid)
  
  ! Set everything back to the way it was?
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 1)

end subroutine plotdiscValues

     type(VECTOR) function velocity(position, grid, startOctal, subcell) RESULT(out)
       use amr_utils_mod, only : returnVelocityVector2
       implicit none

       type(VECTOR), intent(in) :: position
       type(gridtype), intent(in) :: grid
       type(octal), pointer, optional :: startOctal
!       type(octal), pointer :: thisoctal
       integer, optional :: subcell

       type(VECTOR), save :: oldposition, oldout = VECTOR(-8.8d88,-8.8d88,-8.8d88)
       integer, save :: savecounter = 0
!$OMP THREADPRIVATE (oldposition, oldout, savecounter)

       if(oldposition .eq. position) then
          savecounter = savecounter + 1
          out = oldout
          return
       endif

       if(molebench) then
          Out = Molebenchvelocity(Position) 
       elseif(molcluster) then
          Out = amrGridVelocity(grid%octreeRoot, position, startOctal = startOctal, &
                                actualSubcell = subcell, linearinterp = .false.)
       elseif(chrisdisc) then
          Out = keplerianVelocity(Position)
       elseif(hhobench) then
          Out = WaterBenchmarkVelocity(Position)
       elseif(agbstar) then
          Out = AGBstarVelocity(Position)
       elseif(ggtau) then
          Out = ggtauvelocity(Position)
       else
          Out = amrGridVelocity(grid%octreeRoot, position, startOctal = startOctal, &
                                actualSubcell = subcell, linearinterp = .false.)
       endif

!       print *, "out", out
       oldout = out
       oldposition = position

  end function velocity


! Calculate the density based on functional form or values stored on AMR grid.
  real(double) function interpolated_Density(position, grid) RESULT(out)

!    use inputs_mod, only : sphdatafilename    
    implicit none
    type(VECTOR), intent(in) :: position
!    type(VECTOR), save :: oldposition
!    real(double) :: oldout = -8.8d88
!    type(VECTOR):: outv
    type(gridtype), intent(in) :: grid
    integer, save :: savecounter = 0
    logical, save :: firsttime = .true.

    !$OMP THREADPRIVATE (savecounter, firstTime)

! whilst all my code is commented out, stop the warnings    
!    subcell = subcell
    firsttime = .false.
!    startoctal => null()
    savecounter = 0

!    if(oldposition .eq. position) then
!       savecounter = savecounter + 1
!       out = oldout
!       return
!    endif
    
    if(molebench) then
       Out = Molebenchdensity(Position) 
    elseif(ggtau) then
       Out = GGtaudensity(Position) 
    elseif(molcluster) then
!       if(present(startoctal)) then
!          if(startoctal%ndepth .lt. 6) then
!             if(firsttime) then 
!                call new_read_sph_data(sphdatafilename)
!                firsttime = .false.
!             endif
!             outv = clusterparameter(position, startoctal, subcell, theparam = 2)
!             out = clusterdensity(position, grid)
!             out = outv%x
!          endif
!       else
          Out = amrGriddensity(position, grid, linearinterp = .true.)
!       endif
    else
       Out = amrGridDensity(position, grid, linearinterp = .true.)
    endif
    
!    oldout = out
!    oldposition = position
    
  end function interpolated_Density

  function calculatetau(grid, thisoctal, subcell, thismolecule, phi, ds) result(tauavg)

    use inputs_mod, only : usedust, realdust
    TYPE(gridtype) :: grid
    TYPE(OCTAL), pointer :: thisoctal
    type(MOLECULETYPE) :: thismolecule
    integer :: subcell
    real(double) :: phi(:), ds(:)
    real(double), pointer :: npops(:)
    real(double) :: nMol, alphanubase, nlower, nupper, alphagas(nray), alphadust, alpha(nray), temptauarray(nray)
    real(double) :: kappaabs
    integer :: iupper, ilower, itrans
    
    real(double) :: tauavg(mintrans)

    npops => thisOctal%newmolecularlevel(1:maxlevel,subcell)
 
    alphadust = 0.d0
    nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
        
    do itrans = 1, mintrans
       ! assume change in lambda is small... if this isn't true then perhaps molecular_mod isn't for you?!
       if(useDust) then
          if(realdust) then
             call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = real(lambda(itrans)), &
                              kappaAbs = kappaAbs)
          else
             kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10
!multiplied by density !cm -> torus
          endif
          alphadust = kappaAbs * 1.d-10
       endif
       
       iUpper = thisMolecule%iTransUpper(iTrans)
       iLower = thisMolecule%iTransLower(iTrans)
       
       nLower = nPops(iLower) * nMol
       nUpper = nPops(iUpper) * nMol
       
       alphanuBase =(nLower * thisMolecule%einsteinBlu(iTrans) - &
            nUpper * thisMolecule%einsteinBul(iTrans))
       
       alphagas(1:nray) = alphanuBase * hCgsOverFourPi * phi(1:nray) 
       alpha(1:nray) = alphagas(1:nray) + alphadust 
     
       temptauArray(1:nRay) = alpha(1:nray) * ds(1:nRay)
       
       tauavg(itrans) = sum(temptauarray) / dble(nray)
    enddo

  end function calculatetau ! solves rate equation in matrix format - equation 10
  
  subroutine updateLevels(grid, nvoxels, fixedrays, maxRMSFracChange)
    
    use inputs_mod, only : tolerance, dongstep, outputconvergence

    integer :: nvoxels
    logical :: fixedrays

    type(GRIDTYPE) :: grid
    real(double) :: maxFracChangePerLevel(minlevel)
    integer :: ConvergenceCounter(4, minlevel) ! 3 different levels of convergence
    real(double) :: avgFracChange(max(2,minlevel-1), 2) ! don't want the uppermost level to count for convergence
    real(double) :: maxavgFracChange, maxRMSFracChange, maxFracChange
    integer :: maxavglevel(1), maxRMSlevel(1)

    character(len=160) :: message

    logical :: ng
    
10006 format(i6,tr2,6(1pe12.3,1x), 27x, 3x,l1,tr3,3(f7.3,tr3),f7.3,tr3,f8.3)
10007 format(i6,tr2,7(1pe12.3,1x), 18x, 3x,l1,tr3,3(f7.3,tr3),f7.3,tr3,f8.3)
10008 format(i6,tr2,8(1pe12.3,1x),  9x, 3x,l1,tr3,3(f7.3,tr3),f7.3,tr3,f8.3)
10009 format(i6,tr2,9(1pe12.3,1x),      3x,l1,tr3,3(f7.3,tr3),f7.3,tr3,f8.3)
10010 format(i6,tr2,9(1pe12.3,1x),      3x,l1,tr3,3(f7.3,tr3),f7.3,tr3,f8.3)
10011 format(i6,tr2,9(1pe12.3,1x),      3x,l1,tr3,3(f7.3,tr3),f7.3,tr3,f8.3)
10012 format(i6,tr2,9(1pe12.3,1x),      3x,l1,tr3,3(f7.3,tr3),f7.3,tr3,f8.3)
    
! initialise to big negative
    maxFracChangePerLevel = -1.d30
    ConvergenceCounter    = 0
    avgFracChange         = 0.d0

! Compares level populations between this and previous iterations
! This actually swaps the level populations as well as returning convergence data
    call swapPops(grid%octreeRoot, maxFracChangePerLevel, avgFracChange, &
         convergenceCounter, fixedrays) 

!! Advise that ng acceleration has occured during this iteration
!! Ng acceleration occurs after 4 iterations when using random rays by which time initray has grown by 2^4 times
    ng = dongstep
!    if(ng .and. mod(ngcounter, accstepgrand) .eq. 0 .and. ngcounter .ne. 0) call writeinfo('Ng Step', TRIVIAL)
! maximum of change in each level averaged over entire grid
    maxavgFracChange = maxval(avgFracChange(1:max(2,minlevel-2),1))
    maxRMSFracChange = maxval(avgFracChange(1:max(2,minlevel-2),2))
! level that it happened in
    maxavglevel = maxloc(avgFracChange(1:max(2,minlevel-2),1)) - 1 ! array index -> level
    maxRMSlevel = maxloc(avgFracChange(1:max(2,minlevel-2),2)) - 1
! Largest change of any level in any subcell
    maxFracChange = MAXVAL(maxFracChangePerLevel(1:max(2,minlevel-2))) 

! Output text to std out re convergence  
    write(message,'(a,1x,f11.7,1x,a,1x,f6.3,1x,a,2x,l1,1x,a,1x,i6,1x,a,1x,i6)') &
         "Maximum fractional change this iteration ", maxFracChange, "tolerance", tolerance, "fixed rays", fixedrays, &
         "nray", nray
    call writeInfo(message,FORINFO)

    write(message,'(a,f11.7,a,i1)') "Average fractional change this iteration  ", &
         maxavgFracChange/real(nVoxels)," in level ",maxavglevel
    call writeInfo(message,FORINFO)
    write(message,'(a,f11.7,a,i1)') "RMS fractional change this iteration      ", &
         sqrt(maxRMSFracChange/real(nVoxels))," in level ",maxRMSlevel
        call writeInfo(message,FORINFO)
    write(message,'(a,f11.7)') "Std Dev                                   ", &
         sqrt(maxRMSFracChange/real(nVoxels)-(maxavgFracChange/real(nVoxels))**2)
    call writeInfo(message,FORINFO)
    call writeinfo("",FORINFO)  

! Fraction of first minlevel-1 levels converged at 5,2,1,0.5 * tolerance respectively
    write(message,'(a,f7.4,a, 14(f7.4,2x))') "Individual levels converged @ ",tolerance * 5," | ",&
            real(convergenceCounter(1,1:min(14,minlevel-1)))/real(nVoxels) 
       call writeInfo(message,FORINFO)         
    write(message,'(a,f7.4,a, 14(f7.4,2x))') "Individual levels converged @ ",tolerance * 2," | ",&
            real(convergenceCounter(2,1:min(14,minlevel-1)))/real(nVoxels) 
       call writeInfo(message,FORINFO)         
    write(message,'(a,f7.4,a, 14(f7.4,2x))') "Individual levels converged @ ",tolerance * 1," | ",&
            real(convergenceCounter(3,1:min(14,minlevel-1)))/real(nVoxels) 
       call writeInfo(message,FORINFO)         
    write(message,'(a,f7.4,a, 14(f7.4,2x))') "Individual levels converged @ ",tolerance * 0.5," | ",&
            real(convergenceCounter(4,1:min(14,minlevel-1)))/real(nVoxels) 
       call writeInfo(message,FORINFO)         
       call writeinfo("",FORINFO)
    
!99.9% of all cells are converged at 0.5 * tolerance. - this is fine
    allCellsConverged = .false.
    if(all(convergenceCounter(1,1:minlevel) > (1.-1e-3) * real(nvoxels))) then 
       allCellsConverged = .true.

       write(message,'(a,6(tr1,f8.5),tr1,i7)') "Cells converged to tolerence", &
            real(convergenceCounter(1,1:min(6,minlevel))) / real(nvoxels), nvoxels
       call writeinfo(message, FORINFO)
    endif

! Write loads of stuff to various convergence files
! Write fracChanges.dat which contains data for each level
    if (outputconvergence .and. writeoutput) then
       open(138,file="fracChanges.dat",position="append",status="unknown")
  
       select case(minlevel)
       case(6)
          write(138,10006) nray, maxFracChangePerLevel(1:minlevel), &
               fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels), &
               maxFracChange, maxavgFracChange/real(nVoxels)
       case(7)
          write(138,10007) nray, maxFracChangePerLevel(1:minlevel), &
               fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels), &
               maxFracChange, maxavgFracChange/real(nVoxels)
       case(8)
          write(138,10008) nray, maxFracChangePerLevel(1:minlevel), &
               fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels), &
               maxFracChange, maxavgFracChange/real(nVoxels)
       case(9)
          write(138,10009) nray, maxFracChangePerLevel(1:minlevel), &
               fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels), &
               maxFracChange, maxavgFracChange/real(nVoxels)
       case(10)
          write(138,10010) nray, maxFracChangePerLevel(1:9), &
               fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels), &
               maxFracChange, maxavgFracChange/real(nVoxels)
       case(11)
          write(138,10011) nray, maxFracChangePerLevel(1:9), &
               fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels), &
               maxFracChange, maxavgFracChange/real(nVoxels)
       case(12)
          write(138,10012) nray, maxFracChangePerLevel(1:9), &
               fixedrays, real(convergenceCounter(1:3,minlevel))/real(nVoxels), &
               maxFracChange, maxavgFracChange/real(nVoxels)
       case default
       end select
       if(minlevel .gt. 12) then
          write(138,10012) nray, maxFracChangePerLevel(1:9), fixedrays, & 
               real(convergenceCounter(1:3,minlevel))/real(nVoxels), maxFracChange, maxavgFracChange/real(nVoxels)
       endif
       
       close(138)
  
       open(139, file="avgChange.dat", position="append", status="unknown")
20     format(i2,tr3,i6,tr3,12(f8.5,1x))
       write(139,20) grand_iter, nray, avgFracChange(1:min(12,minlevel-1),1)/real(nvoxels)
    
       close(139)
       
       open(140, file="avgRMSchange.dat", position="append", status="unknown")
       write(140,20) grand_iter, nray, sqrt(avgFracChange(1:min(12,minlevel-1),2)/real(nvoxels))
       
       close(140)
    endif

end subroutine updateLevels

SUBROUTINE sobseq(x,init)
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

  !$OMP THREADPRIVATE (fac, in, iv, ip, ix, mdeg)

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
     if (j > MAXBIT) call torus_abort('MAXBIT too small in sobseq')
     im=(j-1)*MAXDIM
     j=min(size(x),MAXDIM)
     ix(1:j)=ieor(ix(1:j),iv(1+im:j+im))
     x(1:j)=ix(1:j)*fac
     in=in+1
  end if
END SUBROUTINE sobseq
      
 subroutine LTEpopsReal(thisMolecule, temperature, levelpops)

   type(MOLECULETYPE) :: thisMolecule
   real(double) :: temperature, nsum
   real :: levelpops(:)
   real(double) :: fac
   integer :: i

   levelpops(1) = 1.d0
   nsum = 1.d0

   do i=2,size(levelpops)
      fac = abs(thisMolecule%energy(i)-thisMolecule%energy(i-1)) / (kev*temperature)
      levelpops(i) = real(max(1d-30,levelpops(i-1) * thisMolecule%g(i)/thisMolecule%g(i-1) * exp(-fac)))
      nsum = nsum + levelpops(i)
   enddo

   levelpops = real(levelpops / nsum)

 end subroutine LTEpopsReal

 subroutine LTEpopsDouble(thisMolecule, temperature, levelpops)

   type(MOLECULETYPE) :: thisMolecule
   real(double) :: temperature, nsum
   real(double) :: levelpops(:), fac
   real(double) :: templevelpops(200)
   integer :: i

   templevelpops(1) = 1.d0
   nsum = 1.d0

   do i=2, thismolecule%nlevels
      fac = abs(thisMolecule%energy(i)-thisMolecule%energy(i-1)) / (kev*temperature)
      templevelpops(i) = max(1d-30,templevelpops(i-1) * thisMolecule%g(i)/thisMolecule%g(i-1) * exp(-fac))
      nsum = nsum + templevelpops(i)
   enddo

   levelpops = templevelpops(1:size(levelpops)) / nsum

 end subroutine LTEpopsDouble

! 21cm line stuff by Dave Acreman
!-----------------------------------------------------------------------------------------------------------

 subroutine make_h21cm_image(grid)
   
   use inputs_mod, only : nsubpixels, itrans, lineImage, maxRhoCalc
   use inputs_mod, only : useDust, isInLte, lowmemory
   use h21cm_mod, only : h21cm_lambda

#ifdef USECFITSIO
   use datacube_mod, only : datacubeFilename
#endif
   use datacube_mod, only : dumpCubeToSpectrum, convertIntensityToBrightnessTemperature

   implicit none

   TYPE(gridtype), intent(in) :: grid
   type(VECTOR)    :: observervec, viewvec, imagebasis(2)
   type(DATACUBE) ::  cube
   type(MOLECULETYPE) :: thisMolecule

   itrans           = 1
   lineImage        = .true.
   maxRhoCalc       = .false. 
   usedust          = .false. 
   isInLte          = .false.
   lowmemory        = .false.

! molecular weight is used for column density calculation
   thisMolecule%molecularWeight = mHydrogen / amu

! Set up 21cm line
   allocate( thisMolecule%transfreq(1) )
   thisMolecule%transfreq(1) = cSpeed / h21cm_lambda

   call writeinfo('Generating H 21cm image', TRIVIAL)
   call setObserverVectors(viewvec, observerVec, imagebasis)

   call createimage(cube, grid, viewvec, observerVec, thismolecule, itrans, nSubpixels, imagebasis)

   call convertSpatialAxes(cube,'kpc')

   call convertIntensityToBrightnessTemperature(cube, h21cm_lambda)

#if USECFITSIO
   if(writeoutput) then 

      call writeinfo("Writing intensity to intensity_"//trim(dataCubeFileName), TRIVIAL)
      call writedatacube(cube, "intensity_"//trim(dataCubeFileName), write_Intensity=.true., &
           write_ipos=.false., write_ineg=.false., write_Tau=.false., write_nCol=.false., write_axes=.false.)

      call writeinfo("Writing column density to nCol_"//trim(dataCubeFileName), TRIVIAL)
      call writedatacube(cube, "nCol_"//trim(dataCubeFileName), write_Intensity=.false., &
           write_ipos=.false., write_ineg=.false., write_Tau=.false., write_nCol=.true., write_axes=.false.)

      call writeinfo("Writing emission weighted velocity to WV_"//trim(dataCubeFileName), TRIVIAL)
      call writedatacube(cube, "WV_"//trim(dataCubeFileName), write_Intensity=.false., write_ipos=.false., &
           write_ineg=.false., write_Tau=.false., write_nCol=.false., write_axes=.false., write_WV=.true.)

   endif
#else
   call writeWarning("TORUS was built without FITS. Cubes will not be written.")
#endif

! Write out global line profile to ASCII file
   call dumpCubeToSpectrum(cube, "globalLineProfile.dat")

 end subroutine make_h21cm_image

!-----------------------------------------------------------------------------------------------------------

subroutine compare_molbench(diffmax)

!  use inputs_mod, only : tolerance

  implicit none 

  character(len=*), parameter :: model_file="results.dat"
  character(len=*), parameter :: bench_file="moltest.dat"
  character(len=*), parameter :: test_file="check.dat"
  character(len=*), parameter :: test2_file="tempcheck.dat"
  character(len=*), parameter :: status_file="status.dat"
  character(len=*), parameter :: diff_file="diff.dat"

! Maximum allowable fractional difference 
  real(double) :: max_diff,diffmax

! Total number of J columns
  integer, parameter :: ncols=8

! No. of columns to check
  integer, parameter :: ncheck=6

  real(double) :: model_R, model_J(ncols), model_Rarray(100)
  real(double) :: bench_R, bench_J(ncols)
  real(double) :: diff(ncols,200)

  integer :: diffmaxloc(2)!, diffmaxr(1)
  integer :: nlines, status, i
! Status = 1 is a pass, 0 is a numerical fail, -1 is an incomplete file
!  max_diff = 2.d0 * tolerance * sqrt(89.d0) ! sqrt(Nvoxels)
  max_diff = 0.05
  diff = -999.
  diffmax = -1.

  open (unit=60, file=bench_file, status='unknown')
  open (unit=61, file=model_file, status='old')
  open (unit=62, file=test2_file, position = "append")
  open (unit=63, file=test_file, position = "append")
  open (unit=64, file=status_file)
  open (unit=65, file=diff_file)
  nlines=0

  do
     nlines = nlines + 1

     read(60, *, iostat=status) bench_R, bench_J(:)
     read(61, *, iostat=status) model_R, model_J(:)


     if (status /= 0 ) then
!        write(*,*) "Reached end of ", model_file, "before end of ", bench_file
!        write(*,*) "TORUS: Test failed"
        status = -1
        nlines = nlines - 1
     end if

     model_Rarray(nlines) = model_R

     diff(1:ncheck,nlines) = abs(model_J(1:ncheck) - bench_J(1:ncheck)) / bench_J(1:ncheck) 
!     write(65,'(7(tr2,es12.5))') model_R, diff(1:ncheck,nlines)

     if(status .eq. -1) exit
     
  end do
  status = 1 ! If you're here then you may have passed

!  write(*,*) "Maximum difference = ",diffmax
!  write(*,*) "Radius = ", model_Rarray(diffmaxloc(2))
!  write(*,*) "Level ", diffmaxloc(1) - 1

  do i = 1,nlines
     if (model_Rarray(i) .gt. 1.5e16) then
        diffmax = maxval(diff(1:ncheck,:))
        diffmaxloc = maxloc(diff(1:ncheck,:))
     endif

     if(any(diff(1:ncheck,i) .gt. max_diff)) status = 0
  enddo

  if(status .ne. 1) then
     write(62,*) diffmax, model_Rarray(diffmaxloc(2)), diffmaxloc(1) - 1
  else
     write(63,*) diffmax, model_Rarray(diffmaxloc(2)), diffmaxloc(1) - 1
  endif

  write(64,'(i2)') status

  close(60)
  close(61)
  close(62)
  close(63)
  close(64)
  close(65)

end subroutine compare_molbench
!-----------------------------------------------------------------------------------------------------------
 recursive subroutine  findtempdiff(grid, thisOctal, thisMolecule, mean, icount)
   use inputs_mod, only : rinner, router
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
                  thisOctal%temperature(subcell) = real(readparameterfrom2dmap(rvec,out,.false.))
                  out = 'nh2'           
                  thisOctal%rho(subcell) = readparameterfrom2dmap(rvec,out,.true.) * 2. * mHydrogen
                  out = 'abundance'
                  thisOctal%molAbundance(subcell) = real(readparameterfrom2dmap(rvec,out,.true.))
               else
                  out = 'abundance'
                  thisOctal%molAbundance(subcell) = real(readparameterfrom2dmap(rvec,out,.true.))
               endif
            else
               thisOctal%temperature(subcell) = tcbr
               thisOctal%rho(subcell) = 1d-20
               thisOctal%molAbundance(subcell) = 1d-20
            endif
         else
            
            out='td'
            thisOctal%temperaturegas(subcell) = real(readparameterfrom2dmap(rvec,out,.false.))
            thisOctal%temperaturedust(subcell) = thisOctal%temperature(subcell)
            thisOctal%temperature(subcell) = thisOctal%temperaturegas(subcell) / thisOctal%temperaturedust(subcell)
            
            out = 'nh2'
            
            thisOctal%temperaturegas(subcell) = real(readparameterfrom2dmap(rvec,out,.true.))
            
            thisOctal%nh2(subcell) = (thisOctal%temperaturegas(subcell) * 2. *mhydrogen) / thisOctal%rho(subcell)
            thisOctal%rho(subcell) = thisOctal%nh2(subcell)
            
            out = 'abundance'
            thisOctal%molAbundance(subcell) = real(readparameterfrom2dmap(rvec,out,.false.))
         end if
      endif
   end do
   
 end subroutine findtempdiff
 
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

!! Not called. Commented out by DA, 24/10/12
!!$ real function nextStep(cube, ivplusone) result (step)
!!$   
!!$   type(DATACUBE) :: cube
!!$   integer :: iv,i,j,ivplusone
!!$   real(double), allocatable :: x(:), y(:), y2(:), intensitysum(:)
!!$   real(double) :: xquad(3),yquad(3),xa(3)
!!$   real(double) :: sumx(4),sumxy(4)
!!$   real(double) :: a(3,3),b(3)
!!$   real(double) :: fac
!!$   real(double), save :: alty2
!!$   real(double), save :: oldstep = 4.d0
!!$   !$OMP THREADPRIVATE (alty2, oldstep)
!!$
!!$   iv = ivplusone - 1
!!$   
!!$   allocate(x(iv))
!!$   allocate(y(iv))
!!$   allocate(y2(iv))
!!$   allocate(intensitysum(iv))
!!$   
!!$   do i = 1,iv
!!$      intensitysum(i) = SUM(cube%intensity(1:cube%nx,1:cube%ny,i))
!!$   enddo
!!$   
!!$   fac = intensitysum(1)
!!$   x = cube%vAxis(1:iv)
!!$   y =  intensitysum / fac 
!!$   
!!$   call spline(x,y,iv,1.d31,1.d31,y2)
!!$   
!!$   xquad = x(iv-3:iv-1)
!!$   yquad = y2(iv-3:iv-1)
!!$   xa = xquad
!!$   
!!$   do i=1,4
!!$      sumx(i) = sum(xa)
!!$      sumxy(i) = sum(xa * yquad)
!!$      
!!$      xa = xa*xquad
!!$   enddo
!!$   
!!$   do i = 1,3
!!$      do j = 1,3
!!$         if(6-i-j .ne. 0) then 
!!$            a(i,j) = sumx(6-i-j)
!!$         else
!!$            a(i,j) = 3.d0
!!$         endif
!!$      enddo; enddo
!!$      
!!$      b = (/sumxy(2),sumxy(1),sum(yquad)/)
!!$      call luSlv(a,b)
!!$
!!$      write(*,*) b(1)*x(iv)**2+b(2)*x(iv)+b(3)
!!$      
!!$      alty2 = b(1)*x(iv)**2+b(2)*x(iv)+b(3)
!!$      
!!$      step = real((oldstep + min(max(25.d0/abs(alty2),0.5_db),4.d0))/2.d0)
!!$      oldstep = step
!!$      
!!$    end function nextStep

subroutine lteintensityAlongRay2(position, direction, grid, thisMolecule, iTrans, deltaV,i0, &
     tau,tautest,rhomax, i0max, nCol, observerVelocity, startI0, startTau, lengthOfRay)

     use inputs_mod, only : useDust, h21cm, densitysubsample, lowmemory
     type(VECTOR) :: position, direction, dsvector
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double), optional :: starti0, startTau, lengthofRay
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0
!     real(double), optional, intent(out) :: nCol,ncol2
     real(double), optional, intent(out) :: nCol
!     integer, optional, intent(in) :: tunable
     type(VECTOR), optional, intent(in) ::  observerVelocity
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(VECTOR) :: currentPosition, thisPosition, endPosition
     type(VECTOR) :: startVel, endVel, thisVel, veldiff
     real(double) :: alphanu1, alphanu2, snu, balance
     real(double) :: alpha
     real(double)  :: dv, deltaV, dVacrossCell
     integer :: i, icount
     real(double) :: tval
     integer :: nTau
     real(double) ::  OneOvernTauMinusOne, ds

     real(double) :: dTau, etaline, dustjnu
     real(double), intent(out), optional :: tau
     real(double), save :: BnuBckGrnd

     real(double) :: phiProfVal
     real         :: sigma_thermal
     real(double) :: deps, origdeps

     logical,save :: firsttime = .true.
     logical,save :: havebeenWarned = .false.
     character(len = 80) :: message
     logical, optional :: tautest
     logical :: dotautest

     real(double) :: dI, dIovercell, attenuateddIovercell, dtauovercell, opticaldepth
     real(double), optional, intent(out) :: rhomax
     real(double), optional, intent(in) :: i0max

     integer :: iupper, ilower
     real(double) :: nlower, nupper
     
!     real(double) :: dpos

     real(double) :: rayLength, maxLengthofRay

!     type(VECTOR) :: curpos

  !$OMP THREADPRIVATE (firstTime, haveBeenwarned, bnuBckGrnd)

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
        disttogrid = 0.
!        call writeinfo("inside",TRIVIAL)
!        write(*,*) position
     else
        distToGrid = distanceToGridFromOutside(grid, position, direction) 
     endif

     if (distToGrid > 1.e29) then
        if (.not. haveBeenWarned) then
           call writewarning("ray does not intersect grid")
           write(message, *) position
           call writeinfo(message, FORINFO)
           write(message, *) direction
           call writeinfo(message, FORINFO)
           havebeenWarned = .true.

        endif
        
        i0 = 1.d-60
        tau = 1.d-60
        goto 666
     endif
     
     deps = 1.d-2 * grid%halfSmallestSubcell ! small value to add on to distance from grid to ensure we get onto grid
     origdeps = deps

     currentPosition = position + (distToGrid + deps) * direction

     i0 = 0.d0
     tau = 0.d0

     maxLengthOfRay = 1.d30
     if (PRESENT(lengthofRay)) maxLengthOfRay = lengthofRay

     if (PRESENT(starti0)) i0 = starti0

     if (PRESENT(startTau)) tau = startTau

     if (present(nCol)) nCol = 0.d0

     if(present(rhomax)) rhomax = 0.d0

     thisOctal => grid%octreeRoot
     icount = 0

     do while((.not. inOctal(grid%octreeRoot, currentPosition)) .and. icount .eq. 0 .and. deps .lt. 1d30)
        deps = 10.d0 * deps
        currentPosition = position + (distToGrid + deps) * direction
     enddo
   
     rayLength = 0.d0
     do while(inOctal(grid%octreeRoot, currentPosition).and.(rayLength < maxLengthOfRay))

        icount = icount + 1

        call findSubcelllocal(currentPosition, thisOctal, subcell)
        if (grid%splitOverMpi.and.(.not.octalonThread(thisOctal,subcell,myrankGlobal))) then
           exit!!!!
        endif
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal, sSubcell = subcell)
        if(tval .le. 0.d0) then
           write(*,*) "This is the end. distanceToCellBoundary is the problem..",tval,direction
           stop
        endif
                  
        if(lowmemory) then                             
           iUpper = thisMolecule%iTransUpper(iTrans)
           iLower = thisMolecule%iTransLower(iTrans)
        
           nlower = thisOctal%molecularLevel(iLower,subcell)! * nMol
           nupper = thisOctal%molecularLevel(iUpper,subcell)! * nMol
                    
           balance = nLower * thisMolecule%einsteinBlu(iTrans) &
                - nUpper * thisMolecule%einsteinBul(iTrans)
        
           etaLine = hCgsOverFourPi * thisMolecule%einsteinA(iTrans)
                    
           etaline = etaLine * nUpper
        endif

        if(present(rhomax)) then
           rhomax = max(rhomax, thisoctal%rho(subcell))
        endif

        dtauovercell = 0.d0
        dIovercell = 0.d0
        attenuateddIovercell = 0.d0

        if(densitysubsample) then
           if ( h21cm) then
              nmol = interpolated_Density(currentposition, grid) / (thisOctal%rho(subcell))
           else
              nmol = thisoctal%molabundance(subcell) * &
                   (interpolated_Density(currentposition, grid) / (2.d0 * mhydrogen))
           end if
        else
           if ( h21cm ) then
              nMol = 1.0
           endif
        end if

        if(lowmemory) then
           nmol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
           etaline = nmol * etaline
        else
           if (.not.associated(thisOctal%molCellParam)) then
              write(*,*) myrankGlobal," bug ", currentPosition, &
                   inOctal(thisOctal,currentPosition), thisOctal%mpiThread(subcell), subcell
              stop
           endif

           nmol = thisoctal%molcellparam(1,subcell)
           etaline = nmol * thisOctal%molcellparam(5,subcell)
        endif

        if(usedust) then
!           alphanu2 = nmol * thisOctal%molcellparam(7,subcell)
!           dustjnu = nmol * thisOctal%molcellparam(8,subcell)
           alphanu2 = thisOctal%molcellparam(7,subcell)
           dustjnu = thisOctal%molcellparam(8,subcell)
        else
           alphanu2 =  0.0
        endif

        thisPosition = currentPosition

        startVel = Velocity(currentPosition, grid, startoctal = thisoctal, subcell = subcell)       
        endposition = currentposition + (1.d0 - 1d-10) * tval * direction

        endVel = Velocity(endPosition, grid, startoctal = thisoctal, subcell = subcell)

        Veldiff = endVel - startVel

        dvAcrossCell = (veldiff.dot.direction)

        if ( h21cm ) then 
           ! Calculate line width in cm/s.
           sigma_thermal = real(sqrt (  (kErg * thisOctal%temperature(subcell)) / mHydrogen))
           ! Convert to Torus units (v/c)
           sigma_thermal = real(sigma_thermal / cspeed)
           dvAcrossCell = abs(dvAcrossCell / sigma_thermal)
        else
           dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))
        end if

        if(densitysubsample) then ! should replace 5 with maxdensity/mindensity * fac
           nTau = min(max(5, nint(dvAcrossCell * 5.d0)), 100) 
        else
           nTau = min(max(5, nint(dvAcrossCell * 5.d0)), 100)
        endif
        
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)
        
        ds = (1.d0 - 1d-10) * tval * OneOvernTauMinusOne
        
! Calculate column density
! Factor of 1.d10 is to convert ds to cm 
        if (present(nCol)) then 
           nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * tval * 1.d10
        end if

        dsvector = ds * direction
        do i = 2, nTau 
           
           thisPosition = thisPosition + dsvector
! THIS LINE CANNOT BE REMOVED
           if(.not. inoctal(grid%octreeroot, thisposition)) thisPosition = thisPosition - 0.99d0 * dsvector
           thisVel = Velocity(thisPosition, grid, startoctal = thisoctal, subcell = subcell)

           if ( present(observerVelocity) ) then 
              thisVel = thisVel - observerVelocity
           end if

           dv = (thisVel .dot. direction) - deltaV
         
           phiProfval = phiProf(dv, thisOctal%molmicroturb(subcell))

           if(densitysubsample .and. .not. h21cm ) then
              nmol = thisoctal%molabundance(subcell) * &
                   (interpolated_Density(thisposition, grid) / (2.d0 * mhydrogen))

              if(lowmemory) then
                 alphanu1 = nmol * hcgsOverFourPi * balance * phiprofval
                 if(usedust) then
                    call writeinfo("Must turn off lowmemory for dust")
                    stop
                 else
                    alphanu2 = 0.0
                 endif
              else
                 alphanu1 = nmol * thisOctal%molcellparam(6,subcell) * phiprofval                 
              endif

           else if (densitysubsample .and. h21cm) then
              nmol     = interpolated_Density(thisposition, grid) / (thisOctal%rho(subcell))
              alphanu1 = nmol * thisOctal%molcellparam(6,subcell) * phiprofval
           else ! no densitysubsample
              if(lowmemory) then
                 alphanu1 = nmol * hcgsOverFourPi * balance * phiprofval
              else
                 alphanu1 = (nmol * thisOctal%molcellparam(6,subcell)) * phiprofval
              endif
           endif

           alpha = alphanu1 + alphanu2
           dTau = alpha * ds * 1.d10

           snu = thisoctal%bnu(itrans,subcell)

           opticaldepth = exp(-tau)
           dI = (1.d0-exp(-dtau))*snu
           dIovercell = dIovercell + dI
           tau = tau + dtau

           dI = opticaldepth * dI
           i0 = i0 + dI
              
        enddo

        
           currentPosition = currentPosition + (tval + origdeps) * direction

           rayLength = rayLength + tVal + origDeps

!        if(writeoutput.and.(inoctal(thisoctal,currentposition)) &
!             .and. (whichsubcell(thisoctal, currentposition) .eq. subcell)) then
!           write(99,*) "fail ",grid%halfSmallestSubcell,tval,origdeps, &
!           sqrt(currentposition%x**2+currentPosition%y**2),currentposition%z
!           call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal, sSubcell = subcell)
!           write(99,*) "tval ",tval
!
!           curpos = currentposition
!!           dpos = min(thisoctal%subcellsize * 0.1, 8. * grid%halfsmallestsubcell)
!           dpos = 16.d0 * grid%halfsmallestsubcell
!
!           do while(dpos .gt. 1e-1 * grid%halfsmallestsubcell)
!              dpos = dpos * 0.5d0
!              if((inoctal(thisoctal,curpos)) .and. (whichsubcell(thisoctal, curpos) .eq. subcell)) then
!                 curpos = curpos + dpos * direction
!              else
!                 curpos = curpos - dpos * direction
!              endif
!           enddo
!
!           if(inoctal(thisoctal,curpos) .and. whichsubcell(thisoctal, curpos) .eq. subcell) then
!              curpos = curpos + 2.d0 * dpos * direction
!           endif
!           currentposition = curpos
!
!        endif

        enddo
         
666  continue

     if(present(rhomax) .and. .not. present(i0max)) then
        return
     else
        if (.not.inOctal(grid%octreeRoot, currentPosition)) &
             i0 = i0 + bnuBckGrnd * exp(-tau) ! from far side
     endif
   end subroutine lteintensityAlongRay2

 subroutine createFluxSpectra(cube,thismolecule,itrans)

   use inputs_mod, only : gridDistance, molAbundance
   character(len=80) :: filename
   type(datacube) :: cube
   type(moleculetype) :: thisMolecule
   real(double),allocatable :: fluxspec(:), inspec(:)
   real(double) :: fac, dA, fac2
   integer :: itrans

   allocate(fluxspec(1:cube%nv))

   allocate(inspec(1:cube%nv))

   fac =  pi*(cube%telescope%Diameter/2.d0)**2 / (2.d0*kerg) ! Effective Telescope Area / 2k - Flux -> Antenna Temp

   dA = ((cube%xaxis(2) - cube%xaxis(1)) / (gridDistance/1.d10)) ** 2

   fac2 = (1.d0* (cspeed**2 / thisMolecule%transFreq(iTrans)**2) &
                    / (pi * (cube%telescope%diameter / 2.d0)**2 )) / dA

   if(writeoutput) then
      write(filename,'(a,a,i1,a,i1,a,f4.1,a,f4.1,a,es9.2e2,a)') trim(cube%telescope%label),'lineprofilesJ', &
        thisMolecule%itransUpper(itrans)-1,'-',thisMolecule%itransLower(itrans)-1,&
        'beam',cube%telescope%beamsize,'diam',cube%telescope%diameter/100.,'abund',molAbundance,'.dat'

   open(30,file=filename,status="unknown",form="formatted")
   write(30,'(a,tr5,a,tr6,a,tr9,a,tr8,a)') "V (km/s)", "Intensity Wm-2sr-1Hz-1", &
                                           "Flux Wm-2","Flux in Jy", "Convolved Flux"

  ! this loop killed to render module inactive. It is currently deprecated and possibly won't make a comeback

!   do i = 1, cube%nv
!      fluxspec(i) = SUM(cube%flux(1:cube%nx,1:cube%ny,i)) 
! This line deprecated as don't store intensity and flux together now

!      inspec(i) = SUM(cube%intensity(1:cube%nx,1:cube%ny,i))
            
!      write(30,'(f6.3,tr8,es13.4,tr6,es13.4,tr8,f11.5,tr8,es13.4)') real(cube%vAxis(i)), &
!inspec(i)/(1e3*(cube%nx**2)), fluxspec(i), fluxspec(i) / thismolecule%transfreq(itrans)* 1e26 
!Jy is not multiplied by frequency!
!   enddo
   
   close(30)
endif

   deallocate(inspec, fluxspec)

 end subroutine createFluxSpectra

 subroutine cubeIntensityToFlux(cube,thisMolecule,itrans,doreverse)

   use inputs_mod, only : gridDistance, nv
   type(moleculetype) :: thismolecule
   type(datacube) :: cube
   integer :: ipixel,jpixel,iv, itrans
   real(double) :: dx
   logical, optional :: doreverse
   logical :: reverse

   if(present(doreverse)) then
      reverse = doreverse
   else
      reverse = .false.
   endif

   dx = cube%xAxis(2) - cube%xAxis(1) ! pixelwidth in torus units
   dx = (dx / (gridDistance*1e-10))**2 ! not in steradians (* 2 pi)

   do ipixel = 1,npixels
      do jpixel = 1,npixels
         do iv = 1,nv
            ! converting intensity to flux which is stored in intensity
            if(reverse) then
               cube%intensity(ipixel,jpixel,iv) = real(cube%intensity(ipixel,jpixel,iv) &
                                                  / (dx * 1e-3 * thismolecule%transfreq(itrans)))
            else
               cube%intensity(ipixel,jpixel,iv) = real(cube%intensity(ipixel,jpixel,iv) &
                                                  * dx * 1e-3 * thismolecule%transfreq(itrans))
            endif
         enddo
      enddo
   enddo

 end subroutine cubeIntensityToFlux

 subroutine GaussianWeighting(cube,npix,FWHM,NormalizeArea)

   use inputs_mod, only : gridDistance
   type(datacube) :: cube
   integer :: npix, i, j
   real :: FWHM
   real(double) :: beamsize,beamsizeInRadians,sigma2,Int,rr
   logical,optional :: NormalizeArea

   !! Calculate Gaussian Weighting 

   beamsize = FWHM
   beamsizeInRadians = beamsize * ArcSecsToRadians !(gridDistance/1d10)*(beamsize/(6.48d5/pi)) 
 ! beamSizeInRadians = 1.22d0*(cspeed/thisMolecule%transFreq(itrans))/15.d2 / 2.d0
   sigma2 = (beamSizeInRadians/2.35d0)**2

 !  write(message,*) "Beam size (in radians): ",beamSizeinRadians
 !  call writeinfo(message,TRIVIAL) 

   Int = 0.d0
   cube%Weight = 0.d0

   do i=1,npix ! Where are all the pixels (top left corners) - what are all their weights
      do j=1,npix
         rr = (cube%xAxis(i)/(gridDistance*1d-10))**2 + & 
              (cube%yAxis(j)/(gridDistance*1d-10))**2

         cube%Weight(i,j) = exp(-0.5d0*(rr/sigma2))

         Int = Int + cube%Weight(i,j)/dble(npix)**2

      enddo
   enddo
! Normalize to npix**2 (because implicit uniform instrument function is 1 in each pixel not 1/npix**2 
   if(present(NormalizeArea) .and. NormalizeArea) then
      cube%Weight = cube%Weight / Int 
   else
      cube%Weight = cube%Weight / maxval(cube%weight(:,:))
   endif
 end subroutine GaussianWeighting

 subroutine fineGaussianWeighting(cube,npix,FWHM,weight,NormalizeArea)
   use inputs_mod, only : gridDistance
   type(datacube) :: cube
   integer :: npix, i, j,k,l
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

   do i=1,npix ! Where are all the pixels (top left corners) - what are all their weights
      do j=1,npix

         do k = 1,19,2
            do l = 1,19,2
               x = cube%xAxis(i) + real(k-10)/20. * deltax  
               y = cube%yAxis(j) + real(l-10)/20. * deltay

               rr = (x/(gridDistance*1d-10))**2 + (y/(gridDistance*1d-10))**2
               f = exp(-0.5d0*(rr/sigma2))
               Weight(i,j) = Weight(i,j) + f*0.01

            enddo
         enddo

         Int = Int + Weight(i,j)/(dble(npix)**2)

      enddo
   enddo
! Normalize to npix**2 (because implicit uniform instrument function is 1 in each pixel not 1/npix**2 
   if(present(NormalizeArea) .and. NormalizeArea) then
      Weight = Weight / Int 
   else
      Weight = Weight / maxval(weight(:,:))
   endif

 end subroutine fineGaussianWeighting

subroutine intensityAlongRay2(position, direction, grid, thisMolecule, iTrans, deltaV,i0, &
     tau,tautest,rhomax, i0max, nCol, observerVelocity)

     use inputs_mod, only : useDust, h21cm, densitysubsample, lowmemory
     type(VECTOR) :: position, direction, dsvector
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0
!     real(double), optional, intent(out) :: nCol,ncol2
     real(double), optional, intent(out) :: nCol
!     integer, optional, intent(in) :: tunable
     type(VECTOR), optional, intent(in) ::  observerVelocity
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(VECTOR) :: currentPosition, thisPosition, endPosition
     type(VECTOR) :: startVel, endVel, thisVel, veldiff
     real(double) :: alphanu1, alphanu2, jnu, snu, balance
     real(double) :: alpha
     real(double)  :: dv, deltaV, dVacrossCell
     integer :: i, icount
     real(double) :: tval
     integer :: nTau
     real(double) ::  OneOvernTauMinusOne, ds

     real(double) :: dTau, etaline, dustjnu
     real(double), intent(out), optional :: tau
     real(double), save :: BnuBckGrnd

     real(double) :: phiProfVal
     real         :: sigma_thermal
     real(double) :: deps, origdeps

     logical,save :: firsttime = .true.
     logical,save :: havebeenWarned = .false.
     character(len = 80) :: message
     logical, optional :: tautest
     logical :: dotautest

     real(double) :: dI, dIovercell, attenuateddIovercell, dtauovercell, opticaldepth, n
     real(double), optional, intent(out) :: rhomax
     real(double), optional, intent(in) :: i0max

     integer :: iupper, ilower
     real(double) :: nlower, nupper
  !$OMP THREADPRIVATE (firstTime, haveBeenwarned, bnuBckGrnd)

!     logical :: lineimage
     
!     real(double) :: dpos
!     type(VECTOR) :: curpos

!     real(double) :: dscounter, tvalcounter
!     real(double) :: lowtau, midtau, hitau,vhitau,xhitau
!     real(double) :: lowi0, midi0, hii0,vhii0,xhii0

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
        disttogrid = 0.
!        call writeinfo("inside",TRIVIAL)
!        write(*,*) position
     else
        distToGrid = distanceToGridFromOutside(grid, position, direction) 
     endif

     if (distToGrid > 1.e29) then
        if (.not. haveBeenWarned) then
           call writewarning("ray does not intersect grid")
           write(message, *) position
           call writeinfo(message, FORINFO)
           write(message, *) direction
           call writeinfo(message, FORINFO)
           havebeenWarned = .true.

        endif
        
        i0 = 1.d-60
        tau = 1.d-60
        goto 666
     endif
     
     deps = 1.d-2 * grid%halfSmallestSubcell ! small value to add on to distance from grid to ensure we get onto grid
     origdeps = deps

     currentPosition = position + (distToGrid + deps) * direction

!     if(inoctal(grid%octreeroot,currentposition)) then
!        call writeinfo("inside",TRIVIAL)
!     else
!        call writeinfo("notinside",TRIVIAL)
!     endif

     i0 = 0.d0
     tau = 0.d0

     if (present(nCol)) nCol = 0.d0
!     if (present(nCol2)) nCol2 = 0.d0

     if(present(rhomax)) rhomax = 0.d0

     thisOctal => grid%octreeRoot
     icount = 0

     do while((.not. inOctal(grid%octreeRoot, currentPosition)) .and. icount .eq. 0 .and. deps .lt. 1d30)
        deps = 10.d0 * deps
        currentPosition = position + (distToGrid + deps) * direction
     enddo
   
     do while(inOctal(grid%octreeRoot, currentPosition))

!        write(*,*) "i0", i0
        icount = icount + 1
!        call writeinfo("check again",TRIVIAL)           
        call findSubcelllocal(currentPosition, thisOctal, subcell)
!        call writeinfo("passed",TRIVIAL)           
!        call findSubcellTD(currentPosition, grid%octreeroot, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal, sSubcell = subcell)
!        call writeinfo("passed again",TRIVIAL)
        if(tval .le. 0.d0) then
           write(*,*) "This is the end. distanceToCellBoundary is the problem.",tval,direction
           stop
        endif

!        currentposition = currentPosition + tval*0.1 * direction
                  
        if(lowmemory) then                             
           iUpper = thisMolecule%iTransUpper(iTrans)
           iLower = thisMolecule%iTransLower(iTrans)
        
           nlower = thisOctal%molecularLevel(iLower,subcell)! * nMol
           nupper = thisOctal%molecularLevel(iUpper,subcell)! * nMol
                    
           balance = nLower * thisMolecule%einsteinBlu(iTrans) &
                - nUpper * thisMolecule%einsteinBul(iTrans)
        
           etaLine = hCgsOverFourPi * thisMolecule%einsteinA(iTrans)
                    
           etaline = etaLine * nUpper
        endif

        if(present(rhomax)) then
           rhomax = max(rhomax, thisoctal%rho(subcell))
        endif

        dtauovercell = 0.d0
        dIovercell = 0.d0
        attenuateddIovercell = 0.d0

        if(densitysubsample) then
           if ( h21cm) then
              nmol = interpolated_Density(currentposition, grid) / (thisOctal%rho(subcell))
           else
              nmol = thisoctal%molabundance(subcell) * &
                   (interpolated_Density(currentposition, grid) / &
                   (2.d0 * mhydrogen))
           end if
        end if

        if ( h21cm ) then
           nMol = 1.0
           etaline = thisOctal%molcellparam(5,subcell)
        else if(lowmemory) then
           nmol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)
           etaline = nmol * etaline
        else
           nmol = thisoctal%molcellparam(1,subcell)
           etaline = nmol * thisOctal%molcellparam(5,subcell)
        endif

        if(usedust) then
!           alphanu2 = nmol * thisOctal%molcellparam(7,subcell)
!           dustjnu = nmol * thisOctal%molcellparam(8,subcell)
           alphanu2 = thisOctal%molcellparam(7,subcell)
           dustjnu = thisOctal%molcellparam(8,subcell)
        else
           alphanu2 =  0.0
        endif

        thisPosition = currentPosition

        startVel = Velocity(currentPosition, grid, startoctal = thisoctal, subcell = subcell)       
        endposition = currentposition + (1.d0 - 1d-10) * tval * direction

!        if(whichsubcell(thisoctal, endposition) .eq. subcell) then
!           endVel = Velocity(endPosition, grid, startoctal = thisoctal, subcell = subcell)
!        else

!           if(.not. (whichsubcell(thisoctal, endposition) .eq. subcell)) write (*,*) "break"
           endVel = Velocity(endPosition, grid, startoctal = thisoctal, subcell = subcell)
!        endif

        Veldiff = endVel - startVel

        dvAcrossCell = (veldiff.dot.direction)

        if ( h21cm ) then 
           ! Calculate line width in cm/s.
           sigma_thermal = real(sqrt (  (kErg * thisOctal%temperature(subcell)) / mHydrogen))
           ! Convert to Torus units (v/c)
           sigma_thermal = real(sigma_thermal / cspeed)
           dvAcrossCell = abs(dvAcrossCell / sigma_thermal)
        else
           dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))
        end if

        if(densitysubsample) then ! should replace 5 with maxdensity/mindensity * fac
           nTau = min(max(5, nint(dvAcrossCell * 5.d0)), 100) 
        else
!           nTau = min(max(2, nint(dvAcrossCell * 5.d0)), 1000) 
!           nTau = min(max(tunable, nint(dvAcrossCell * 5.d0)), 10) 
           nTau = min(max(5, nint(dvAcrossCell * 5.d0)), 1000)
        endif
        
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)
        
        ds = (1.d0 - 1d-10) * tval * OneOvernTauMinusOne
        
! Calculate column density
! Factor of 1.d10 is to convert ds to cm 
        if (present(nCol)) then 
!           nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * ds * 1.d10
           nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * tval * 1.d10
!           ncol2 = ncol2 + (interpolated_Density(thisposition, grid, thisoctal, subcell) / &
!                           (thisMolecule%molecularWeight * amu)) * tval * 1.d10
!           tvalcounter = tvalcounter + tval
        end if

!        if(present(sink)) then
!           if(thisOctal%rho(subcell) .gt. 1d-14) sink = 1.
!        endif

        dsvector = ds * direction

        do i = 2, nTau 

!           if (present(nCol2)) then 
!           nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * ds * 1.d10
!              nCol = nCol + (thisOctal%rho(subcell) / (thisMolecule%molecularWeight * amu) ) * ds * 1.d10
!              ncol2 = ncol2 + (interpolated_Density(thisposition, grid, thisoctal, subcell) / &
!                              (thisMolecule%molecularWeight * amu)) * ds * 1.d10
!              dscounter = dscounter + ds
!           end if
           
           thisPosition = thisPosition + dsvector
           if(.not. inoctal(grid%octreeroot, thisposition)) thisPosition = thisPosition - 0.99d0 * dsvector
! THIS LINE CANNOT BE REMOVED
           thisVel = Velocity(thisPosition, grid, startoctal = thisoctal, subcell = subcell)

           if ( present(observerVelocity) ) then 
              thisVel = thisVel - observerVelocity
           end if

           dv = (thisVel .dot. direction) - deltaV

           if ( h21cm ) then 
              phiprofval = gauss (sigma_thermal, real(dv) ) / thisMolecule%transfreq(1)
           else  
              phiProfval = phiProf(dv, thisOctal%molmicroturb(subcell))
           end if
           
           if(densitysubsample .and. .not. h21cm ) then
              nmol = thisoctal%molabundance(subcell) * &
              (interpolated_Density(thisposition, grid) / (2.d0 * mhydrogen))

              if(lowmemory) then
                 etaline = nmol * etaline
                 alphanu1 = nmol * hcgsOverFourPi * balance * phiprofval
                 if(usedust) then
                    call writeinfo("Must turn off lowmemory for dust")
                    stop
                 else
                    alphanu2 = 0.0
                 endif
              else
                 etaline = nmol * thisOctal%molcellparam(5,subcell)
                 alphanu1 = nmol * thisOctal%molcellparam(6,subcell) * phiprofval
                 
                 if(usedust) then
                    alphanu2 = thisOctal%molcellparam(7,subcell)
                    dustjnu = thisOctal%molcellparam(8,subcell)
!                    alphanu2 = nmol * thisOctal%molcellparam(7,subcell)
!                    dustjnu =  nmol * thisOctal%molcellparam(8,subcell)
                 else
                    alphanu2 = 0.0
                 endif
              endif

           else if (densitysubsample .and. h21cm) then
              nmol     = interpolated_Density(thisposition, grid) / (thisOctal%rho(subcell))
              etaline  = nmol * thisOctal%molcellparam(5,subcell)
              alphanu1 = nmol * thisOctal%molcellparam(6,subcell) * phiprofval
           else ! no densitysubsample
              if(lowmemory) then
                 alphanu1 = nmol * hcgsOverFourPi * balance * phiprofval
              else
                 alphanu1 = (nmol * thisOctal%molcellparam(6,subcell)) * phiprofval
              endif
           endif

           alpha = alphanu1 + alphanu2
           dTau = alpha * ds * 1.d10

!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e3 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e4) &
!                lowtau = lowtau + dtau
!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e4 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e5) &
!                midtau = midtau + dtau
!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e5 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e6) &
!                hitau = hitau + dtau
!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e6 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e7) &
!                vhitau = vhitau + dtau
!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e7 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e10) &
!                xhitau = xhitau + dtau
         
           jnu = etaLine * phiProfVal

           if(useDust) jnu = jnu + dustjnu

           if (alpha .ne. 0.d0) then
              snu = jnu/alpha
!           elseif (alpha .lt. 0.d0) then
!              write(*,*) "no"
           else
              snu = tiny(snu)
              i0 = i0 + tiny(i0)
           endif

           opticaldepth = exp(-tau)
           dI = (1.d0-exp(-dtau))*snu
           dIovercell = dIovercell + dI
           tau = tau + dtau

           dI = opticaldepth * dI

!           if(dI .gt. i0 * 1d-10 .and. tau .lt. 100) then
              i0 = i0 + dI

!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e3 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e4) &
!                lowi0 = lowi0 + dI
!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e4 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e5) &
!                midi0 = midi0 + dI
!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e5 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e6) &
!                hii0 = hii0 + dI
!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e6 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e7) &
!                vhii0 = vhii0 + dI
!           if(thisoctal%rho(subcell) / (amu * 2.) .gt. 1e7 .and. thisoctal%rho(subcell) / (amu * 2.) .lt. 1e10) &
!                xhii0 = xhii0 + dI

!           else
!              i0 = i0
!           endif
           enddo

           if(present(i0max) .and. i0 .gt. 0.99d0 * i0max) then 
              exit
           endif

           if(present(i0max)) then

              attenuateddIovercell = attenuateddIovercell + dI
              dtauovercell = dtauovercell + dtau
              n = thisoctal%newmolecularlevel(4,subcell)
! average intensity in this cell              
              thisoctal%newmolecularlevel(5,subcell) = &
                   (n * thisoctal%newmolecularlevel(5,subcell) + dtauovercell) / (n + 1.d0) 
! average intensity in this cell
              thisoctal%newmolecularlevel(1,subcell) = &
                   (n * thisoctal%newmolecularlevel(1,subcell) + dIovercell) / (n + 1.d0) 
! average intensity in this cell (attenuated)              
              thisoctal%newmolecularlevel(2,subcell) = &
                   (n * thisoctal%newmolecularlevel(2,subcell) + attenuateddIoverCell) / (n + 1.d0)
! average attenuated intensity at this cells edge 
              thisoctal%newmolecularlevel(3,subcell) = &
                   (n * thisoctal%newmolecularlevel(3,subcell) + i0) / (n + 1.d0) 
              thisoctal%newmolecularlevel(4,subcell) = thisoctal%newmolecularlevel(4,subcell) + 1.d0
           endif
        
           currentPosition = currentPosition + (tval + origdeps) * direction
        
     enddo
         
666  continue

     if(present(rhomax) .and. .not. present(i0max)) then
        return
     else
        i0 = i0 + bnuBckGrnd * exp(-tau) ! from far side
     endif
     
!     write(*,*) tvalcounter, dscounter
!     write(*,*) lowtau,midtau,hitau,vhitau,xhitau
!     write(*,*) lowi0,midi0,hii0,vhii0,xhii0,i0

   end subroutine intensityAlongRay2

! this routine writes a file of intensity, density, cellsize, cellvolume etc

    subroutine dumpIntensityContributions(grid, thisMolecule) 
     use mpi_global_mod, only:  myRankGlobal
     use parallel_mod
      use inputs_mod, only : itrans
      type(MOLECULETYPE) :: thisMolecule
      type(GRIDTYPE) :: grid
      type(VECTOR) :: viewVec, rayStart
      integer :: nx, ny, nv
      real(double) :: x, y, z, v, tau
      real(double) :: vmax, vmin, nuStart, nuEnd, deltaNu
      real(double) :: itot, i0, icheck
      integer :: i, j ,k, iv1, iv2
      logical, save :: firstTime=.true.
      !$OMP THREADPRIVATE(firstTime)

      viewVec = VECTOR(0.d0, 1.d0, 0.d0)
      nx = 2000
      ny = 2000
      nv = 20

      vmin = -2.d0*(1.d5)/cspeed
      vmax =  2.d0*(1.d5)/cspeed
      iv1 = 1
      iv2 = nv
      

#ifdef MPI
    iv1 = (myrankglobal) * (nv / (nThreadsGlobal)) + 1
    iv2 = (myrankglobal+1) * (nv / (nThreadsGlobal))
    if (myrankglobal == (nThreadsGlobal-1)) iv2 = nv
#endif



      nuStart = thisMolecule%transfreq(itrans) * (1.d0+vmin)
      nuEnd = thisMolecule%transfreq(itrans) * (1.d0+vmax)
      deltaNu = nuEnd - nuStart
      itot = 0.d0
      do i = 1, nx
         if (myrankglobal == 0) write(*,*) "i ",i
         do j = 1, ny 
            do k = iv1, iv2
               x = -grid%octreeRoot%subcellSize + (2.d0*grid%octreeRoot%subcellSize)*dble(i-1)/dble(nx-1)
               z = -grid%octreeRoot%subcellSize + (2.d0*grid%octreeRoot%subcellSize)*dble(j-1)/dble(ny-1)
               y = -grid%octreeRoot%subcellSize
               v = vMin + (vMax-vMin)*dble(k-1)/dble(nv-1)

               if(firstTime) then
                  call writeinfo("Filling Octal parameters for first time",TRIVIAL)
                  call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule)
                  firstTime = .false.
               endif

               rayStart = VECTOR(x, y, z)
               call intensityAlongRay(rayStart, viewvec, grid, thisMolecule, iTrans, v, i0, tau)
            enddo
         enddo
      enddo
      itot = 0.d0
      call sumDi(grid%octreeRoot, deltaNu, itot)
      if (myrankGlobal == 0) then
         open(33, file="contribs.dat", status="unknown",form="formatted")
         icheck = 0.d0
         call writeContributions(grid%octreeRoot, deltaNu, itot, icheck)
         write(*,*) "sanity check ",icheck/itot, icheck,itot
         close(33)
      endif
      call torus_mpi_barrier()
    end subroutine dumpIntensityContributions


 recursive subroutine sumDi(thisOctal, deltaNu, itot)
#ifdef MPI
   use mpi
   real(double) :: temp
   integer :: ierr
#endif
   type(octal), pointer   :: thisOctal
   type(octal), pointer  :: child 
   integer :: subcell, i
   real(double) :: deltaNu, itot
   do subcell = 1, thisOctal%maxChildren
      if (thisOctal%hasChild(subcell)) then
         ! find the child
         do i = 1, thisOctal%nChildren, 1
            if (thisOctal%indexChild(i) == subcell) then
               child => thisOctal%child(i)
               call sumdi(child, deltaNu, itot)
               exit
            end if
         end do
      else
#ifdef MPI
         temp = 0.d0
         call MPI_REDUCE(thisOctal%newMolecularlevel(1,subcell),temp,1,MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
         thisOctal%newMolecularLevel(1,subcell) = temp
#endif
         itot = itot + thisOctal%newMolecularlevel(1, subcell)
      endif
   end do
   
 end subroutine sumDi

 recursive subroutine writeContributions(thisOctal, deltaNu, itot, icheck)

   type(octal), pointer   :: thisOctal
   type(octal), pointer  :: child 
   integer :: subcell, i
   real(double) :: deltaNu, itot,icheck
   do subcell = 1, thisOctal%maxChildren
      if (thisOctal%hasChild(subcell)) then
         ! find the child
         do i = 1, thisOctal%nChildren, 1
            if (thisOctal%indexChild(i) == subcell) then
               child => thisOctal%child(i)
               call writeContributions(child, deltaNu, itot, icheck)
               exit
            end if
         end do
      else
         if (thisOctal%rho(subcell) > 1.d-30) then
            write(33, '(1p,5e12.3)') thisOctal%rho(subcell)/(2.d0*mhydrogen), &
                 thisOctal%newMolecularlevel(1, subcell), cellVolume(thisOctal,subcell), &
                 thisOctal%rho(subcell)*cellVolume(thisOctal,subcell)

            icheck = icheck + thisOctal%newMolecularlevel(1, subcell)
         endif

      endif
   end do
   
 end subroutine writeContributions

 subroutine  photoionChemistry(grid, thisOctal, subcell)
   use inputs_mod, only : molAbundance
   type(GRIDTYPE) :: grid
   type(octal), pointer   :: thisOctal
   integer :: subcell
   
   if (grid%splitOverMPI.and.(.not.octalOnThread(thisOctal, subcell, myrankGlobal))) goto 666
         
   thisOctal%molAbundance(subcell) = molAbundance
   if (thisOctal%temperature(subcell) > 100.d0) then
      thisOctal%molAbundance(subcell) = tiny(thisOctal%molAbundance)
   endif
666  continue
   end subroutine photoionChemistry
 
#ifdef MPI

 subroutine intensityAlongRaySplitOverMPI(position, direction, grid, thisMolecule, iTrans, deltaV,i0)
   use mpi
   type(VECTOR) :: position, direction
   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   real(double) :: disttoGrid
   integer :: itrans
   real(double), intent(out) :: i0
   type(OCTAL), pointer :: thisOctal
   integer :: subcell
   type(VECTOR) :: currentPosition
   real(double)  :: deltaV
   real(double) :: deps, origdeps
   
   logical,save :: havebeenWarned = .false.
   character(len = 80) :: message
   
!     real(double) :: dscounter, tvalcounter
!     real(double) :: lowtau, midtau, hitau,vhitau,xhitau
!     real(double) :: lowi0, midi0, hii0,vhii0,xhii0
   real(double) :: tempStorage(2)
   real(double) :: loc(11)
   integer :: status(MPI_STATUS_SIZE)
   integer, parameter :: tag = 50
   integer :: ierr
   integer :: nStorage
   real(double) :: lengthOfRay, tau
   integer :: iThread
   integer :: i

   !$OMP THREADPRIVATE (haveBeenwarned)
   i = thisMolecule%nTrans


   if(inOctal(grid%octreeRoot, Position)) then
      disttogrid = 0.
!     call writeinfo("inside",TRIVIAL)
!     write(*,*) position
   else
      distToGrid = distanceToGridFromOutside(grid, position, direction) 
   endif
   
   if (distToGrid > 1.e29) then
      if (.not. haveBeenWarned) then
         call writewarning("ray does not intersect grid")
         write(message, *) position
         call writeinfo(message, FORINFO)
         write(message, *) direction
         call writeinfo(message, FORINFO)
         havebeenWarned = .true.
         
      endif
      
      i0 = 1.d-60
      goto 666
   endif
   
   deps = 5.d-4 * grid%halfSmallestSubcell ! small value to add on to distance from grid to ensure we get onto grid
   origdeps = deps
   
   currentPosition = position + (distToGrid + deps) * direction
   
!     if(inoctal(grid%octreeroot,currentposition)) then
!        call writeinfo("inside",TRIVIAL)
!     else
!        call writeinfo("notinside",TRIVIAL)
!     endif

   i0 = 0.d0
   tau = 0.d0
   nstorage = 2
   
   thisOctal => grid%octreeRoot
   do while (inOctal(grid%octreeRoot, currentPosition))
      call findSubcellLocal(currentPosition, thisOctal, subcell)
      call distanceToCellBoundary(grid, currentPosition, direction, lengthOfRay, &
           sOctal=thisOctal, sSubcell = subcell)
      iThread = thisOctal%mpiThread(subcell)
      loc(1) = currentPosition%x
      loc(2) = currentPosition%y
      loc(3) = currentPosition%z
      loc(4) = direction%x
      loc(5) = direction%y
      loc(6) = direction%z
      loc(7) = lengthOfRay
      loc(8) = i0
      loc(9) = deltaV
      loc(10) = dble(itrans)
      loc(11) = tau
      
!        write(*,*) myrankGlobal, " sending message to ", ithread
      call MPI_SEND(loc, 11, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
!        write(*,*) myrankGlobal, " message sent awaiting recv"
      call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, status, ierr)
!        write(*,*) myrankGlobal, " message received from ",ithread

      i0 = tempStorage(1)
      tau = tempStorage(2)
      
      currentPosition = currentPosition + (lengthofRay+1.d-3*grid%halfSmallestSubcell)*direction
   enddo
666 continue
 end subroutine intensityAlongRaySplitOverMPI

 subroutine intensityAlongRayServer(grid, thisMolecule)
   use mpi
   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   type(VECTOR) :: startPosition, direction
   real(double) :: startIntensity, deltaV
   logical :: stillServing
   integer :: nStorage
   real(double) :: tempStorage(2)
   real(double) :: loc(11)
   integer :: status(MPI_STATUS_SIZE)
   integer, parameter :: tag = 50
   integer :: ierr
   integer :: itrans
   real(double) :: lengthOfray, i0, tau, starttau
   
   stillServing = .true.
   do while (stillServing)
!        write(*,*) myrankGlobal, " waiting for send from thread 0"
      call MPI_RECV(loc, 11, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, status, ierr)
!        write(*,*) myrankGlobal, " received a message from thread 0"
      startposition%x = loc(1)
      startposition%y = loc(2)
      startposition%z = loc(3)
      direction%x = loc(4)
      direction%y = loc(5)
      direction%z = loc(6)
      lengthOfRay = loc(7)
      startIntensity = loc(8)
      deltaV = loc(9)
      itrans = int(loc(10))
      starttau = loc(11)
      
      if (startposition%x > 1.d29) then
         stillServing= .false.
         write(*,*) myrankGlobal, " received server shutdown signal"
      else
         
         call lteintensityAlongRay2(startposition, direction, grid, thisMolecule, iTrans, deltaV, i0, tau=tau, &
              startI0=startIntensity, startTau=starttau, lengthOfRay=lengthofRay)
         tempStorage(1) = i0
         tempStorage(2) = tau
         nStorage = 2
         call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, ierr)
      endif
   enddo
 end subroutine intensityAlongRayServer

 subroutine shutdownServers()
   use mpi
   integer :: iThread
   real(double) :: loc(10)
   integer, parameter :: tag = 50
   integer :: ierr
   
   do iThread = 1, nThreadsGlobal-1
      if (iThread /= myrankGlobal) then
         loc = 0.d0
         loc(1) = 1.d30
         call MPI_SEND(loc, 10, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
      endif
   enddo
 end subroutine shutdownServers
 
#endif

end module molecular_mod
#endif

