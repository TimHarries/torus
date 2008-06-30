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
   use grid_mod
   use math_mod
   use datacube_mod
   use nrutil
   use nrtype
   use parallel_mod
 !  use mkl_vsl_type
 !  use mkl_vsl

   implicit none

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

              allocate(thisOctal%newmolecularLevel(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
              thisOctal%newmolecularLevel = 1.d-20
              allocate(thisOctal%oldmolecularLevel(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
              thisOctal%oldmolecularLevel = 1.d-20

              if (.not.associated(thisOctal%bnu)) then
                 allocate(thisOctal%bnu(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
              endif
              thisOctal%bnu(subcell,i) = bnu(thisMolecule%transFreq(i), dble(thisOctal%temperature(subcell)))

           else
              if (.not. associated(thisOctal%newmolecularLevel)) then
                 allocate(thisOctal%newmolecularLevel(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
              endif
              thisOctal%newmolecularLevel = 1.d-20

              if (.not.associated(thisOctal%oldmolecularLevel)) then
                 allocate(thisOctal%oldmolecularLevel(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
              endif
              thisOctal%oldmolecularLevel = 1.d-20

              if (.not.associated(thisOctal%molecularLevel)) then
                 allocate(thisOctal%molecularLevel(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
              endif

              if(lte) then           
                 call LTEpops(thisMolecule, dble(thisOctal%temperature(subcell)), &
                              thisOctal%molecularLevel(subcell,1:thisMolecule%nlevels))
              else      
                 call LTEpops(thisMolecule, tcbr, thisOctal%molecularLevel(subcell,1:thisMolecule%nlevels))

                 if((grid%geometry .eq. "h2obench1") .or. (grid%geometry .eq. "h2obench2")) then
                    thisOctal%molecularLevel(subcell,1) = 1.0
                    thisOctal%molecularLevel(subcell,2) = 0.0000
                 endif
              endif

              if (.not.associated(thisOctal%bnu)) then
                 allocate(thisOctal%bnu(1:thisOctal%maxChildren, 1:thisMolecule%nTrans))
              endif

              do i = 1, thisMolecule%nTrans
                 thisOctal%bnu(subcell,i) = bnu(thisMolecule%transFreq(i), dble(thisOctal%temperature(subcell)))
              enddo

              if (.not.associated(thisOctal%jnu)) then
                 allocate(thisOctal%jnu(1:thisOctal%maxChildren, 1:thisMolecule%ntrans))
              endif

              if(lte) then
                 do i = 1, thisMolecule%ntrans
                    thisOctal%jnu(subcell,i) = thisOctal%bnu(subcell, i)

 !                if(thisOctal%temperature(subcell) .gt. Tcbr) then
 !                   thisOctal%jnu(subcell,i) = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell) * &
 !                                              thisOctal%molecularLevel(subcell,i) * thisMolecule%einsteinA(i) * hCgsOverfourPi * &
 !                                              thisOctal%microturb(subcell) / sqrtPi

                 enddo
              else
                 do i = 1, thisMolecule%ntrans
                    thisOctal%jnu(subcell,:) = bnu(thisMolecule%transFreq(i), tcbr)
                 enddo
              endif


           endif

           if (.not.associated(thisOctal%molAbundance)) then
              allocate(thisOctal%molAbundance(1:thisOctal%maxChildren))
              thisOctal%molAbundance(subcell) = thisMolecule%abundance
           endif

          if(thisOctal%temperaturedust(subcell) .eq. thisOctal%temperaturegas(subcell)) then
             thisOctal%temperaturegas(subcell) = thisOctal%temperature(subcell)
             thisOctal%temperaturedust(subcell) = thisOctal%temperature(subcell)
          endif

          thisOctal%molmicroturb(subcell) = 1. / thisOctal%microturb(subcell)

       endif
    enddo
   end subroutine allocateMolecularLevels

 ! solves rate equation in matrix format - equation 10
   subroutine solveLevels(nPops, jnu,  temperature, thisMolecule, nh2)
     real(double) :: nPops(:)
     real(double) :: temperature
     real(double) :: jnu(:)
     real(double) :: nh2
     type(MOLECULETYPE) :: thisMolecule
     real(double), allocatable :: matrixA(:,:), matrixB(:), collMatrix(:,:), cTot(:)
     real(double) :: boltzFac
     integer :: nLevels
     integer :: i, j
     integer :: itrans, l, k, iPart
     real(double) :: collEx, colldeEx

     nLevels = thisMolecule%nLevels

     allocate(matrixA(1:nLevels+1,1:nLevels+1))
     allocate(matrixB(1:nLevels+1))

     matrixA = 1.d-10 ! Initialise rates to negligible to avoid divisions by zero
     matrixB = 1.d-10 ! Solution vector - all components (except last) => equilibrium

     matrixB(nLevels+1) = 1.d0 ! Sum over all levels = 1 - Conservation constraint
     allocate(collMatrix(1:nLevels+1, 1:nLevels+1))

 ! This do loop calculates the contribution to transition rates of each level from every other level. NB emission is +ve here
     do iTrans = 1, thisMolecule%nTrans

        k = thisMolecule%iTransUpper(iTrans)
        l = thisMolecule%iTransLower(iTrans)

        matrixA(k,k) = matrixA(k,k) + thisMolecule%einsteinBul(iTrans) * jnu(iTrans) + thisMolecule%einsteinA(iTrans)
        matrixA(l,l) = matrixA(l,l) + thisMolecule%einsteinBlu(iTrans) * jnu(iTrans)
        matrixA(k,l) = matrixA(k,l) - thisMolecule%einsteinBlu(iTrans) * jnu(iTrans)
        matrixA(l,k) = matrixA(l,k) - thisMolecule%einsteinBul(iTrans) * jnu(iTrans) - thisMolecule%einsteinA(iTrans)

     enddo

     collMatrix = 1.d-10

 ! Calculate contribution from collisions - loop over all collision partners and all molecular levels

     do iPart = 1, thisMolecule%nCollPart
        do iTrans = 1, thisMolecule%nCollTrans(iPart)

           k = thisMolecule%iCollUpper(iPart, iTrans)
           l = thisMolecule%iCollLower(iPart, iTrans)

           boltzFac =  exp(-abs(thisMolecule%energy(k)-thisMolecule%energy(l)) / (kev*temperature))
           colldeEx = collRate(thisMolecule, temperature, iPart, iTrans) * nh2! * (nh2 * thisOctal%molAbundance(subcell))!!!! What's happened here?
           collEx = colldeEx * boltzFac * thisMolecule%g(k)/thisMolecule%g(l)

           collMatrix(l, k) = collMatrix(l, k) + collEx
           collMatrix(k, l) = collMatrix(k, l) + colldeEx

        enddo
     enddo

     allocate(cTot(1:thisMolecule%nLevels))
     cTot = 0.d0

     do k = 1, thisMolecule%nLevels
        do l = 1, thisMolecule%nLevels
           cTot(k) = cTot(k) + collMatrix(k,l) ! sum over all collisional rates out of each level
        enddo
     enddo

     do i = 1, thisMolecule%nLevels
           matrixA(i,i) = matrixA(i,i) + cTot(i)
        do j = 1, thisMolecule%nLevels
           if (i /= j) then
              matrixA(i,j) = matrixA(i,j) - collMatrix(j, i)
           endif
        enddo
     enddo

     matrixA(nLevels+1,1:nLevels+1) = 1.d0 ! sum of all level populations
     matrixA(1:nLevels+1,nLevels+1) = 1.d-20 ! fix highest population to small non-zero value
     
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
   call luSlv(matrixA, matrixB, nLevels+1)

  if(.not.isnan(matrixB(1))) then
      nPops(1:nLevels) = matrixB(1:nLevels)
    else
       nPops(1:nLevels) = matrixB(1:nLevels)
       Write(*,*) "bad",temperature,npops(1),npops(2),npops(3),npops(4)     
       call LTEpops(thisMolecule, temperature, npops)
       if(writeoutput) write(*,*) "BADNESS AVOIDED"
      Write(*,*) "good",temperature,npops(1),npops(2),npops(3),npops(4)     
    endif

   deallocate(matrixA, matrixB, collMatrix, cTot)

 end subroutine solveLevels

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

 ! Calculate collision rates between partners for given temperature
   real(double) function collRate(thisMolecule, temperature, iPart, iTrans)
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: temperature, r
     integer :: iTrans, k, iPart

     collRate = 0.d0

     call locate(thisMolecule%collTemps(iPart,1:thisMolecule%nCollTemps(iPart)), &
          thisMolecule%nCollTemps(iPart), temperature, k)

     r = (temperature - thisMolecule%collTemps(iPart,k)) / &
          (thisMolecule%collTemps(iPart,k+1) - thisMolecule%collTemps(iPart, k))

     collRate = collRate + thisMolecule%collRates(iPart, iTrans, k) + &
          r * ( thisMolecule%collRates(iPart, iTrans, k+1) -  thisMolecule%collRates(iPart, iTrans, k))

   end function collRate

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
           call solveLevels(thisOctal%molecularLevel(subcell,1:thisMolecule%nLevels), &
                thisOctal%jnu(subcell,1:thisMolecule%nTrans),  &
                dble(thisOctal%temperature(subcell)), thisMolecule, thisOctal%nh2(subcell))
        endif
     enddo
   end subroutine solveAllPops

 ! this subroutine calculates the maximum fractional change in the first 6 energy levels
   recursive subroutine  swapPops(thisOctal, maxFracChangePerLevel, avgFracChange, counter, &
                          iter, nVoxels,itransdone,fixedrays)

     use input_variables, only : tolerance

     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i
     real(double) :: maxFracChangePerLevel(:), printmaxFracChange, globalmaxFracChange, avgFracChange(:,:)
     real(double), allocatable :: diff(:), newFracChangePerLevel(:), temp(:)
     integer :: counter(:,:),j


     integer :: iter
     integer,save :: l
     integer :: nVoxels

     logical :: itransDone(:),fixedrays
     integer :: ntrans

     ntrans = size(maxFracChangePerLevel)

     allocate(newFracChangePerLevel(ntrans))
     allocate(temp(ntrans))

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call swapPops(child, maxFracChangePerLevel, avgFracChange, counter, iter, &
                      nVoxels, itransdone, fixedrays)
                 exit
              end if
           end do
        else

          globalMaxFracChange = MAXVAL(maxFracChangePerLevel(1:ntrans)) 
          printmaxFracChange = MAXVAL(maxFracChangePerLevel(1:6),mask = (.not. itransdone(1:6)))

          newFracChangePerLevel = abs(((thisOctal%newMolecularLevel(subcell,1:ntrans) - &
               thisOctal%molecularLevel(subcell,1:ntrans)) / &
               thisOctal%newmolecularLevel(subcell,1:ntrans)))

           temp = newFracChangePerLevel(1:6)

           do j=1, ntrans
              if(.not. itransdone(j)) then
                 if(newFracChangePerLevel(j) < 0.5 * tolerance) counter(4,j) = counter(4,j) + 1 ! used for itransdone
                 if(newFracChangePerLevel(j) < tolerance) counter(1,j) = counter(1,j) + 1
                 if(newFracChangePerLevel(j) < 2. * tolerance) counter(2,j) = counter(2,j) + 1
                 if(newFracChangePerLevel(j) < 5. * tolerance) counter(3,j) = counter(3,j) + 1
              else
                 counter(:,j) = counter(:,j) + 1
              endif
           enddo

              if(maxval(temp) < 5. * tolerance) then
                 counter(3,ntrans + 1) = counter(3,ntrans + 1) + 1
                 if(maxval(temp) < 2. * tolerance) then
                    counter(2,ntrans + 1) = counter(2,ntrans + 1) + 1
                    if(maxval(temp) < tolerance) then
                       counter(1,ntrans + 1) = counter(1,ntrans + 1) + 1
                    endif
                 endif
              endif

              avgFracChange(:,1) = avgFracChange(:,1) + temp
              avgFracChange(:,2) = avgFracChange(:,2) + temp**2

           If (maxval(temp, mask = .not. itransdone(1:6)) > printmaxFracChange) then
              maxFracChangePerLevel = newFracChangePerLevel ! update maxFracChange if fractional change is great => not converged yet
           endif

           allocate(diff(size(thisOctal%molecularlevel(subcell,:))))

           diff = thisOctal%molecularLevel(subcell,:) - &
                  thisOctal%newmolecularLevel(subcell,:)

 !          if(fixedrays .and. modulus(SubcellCentre(thisOctal,subcell)) .lt. 1e7 .and. modulus(SubcellCentre(thisOctal,subcell)) .gt. 3e6) then
           if(.not. fixedrays) then
              if (maxval(temp, mask = (.not. itransdone(1:6))) .gt. tolerance * 6.) then
                 thisOctal%molecularLevel(subcell,:) = &
                      thisOctal%newmolecularLevel(subcell,:) - (0.0)*diff! (enhanced) update molecular levels
              else
                 thisOctal%molecularLevel(subcell,:) = &
                      thisOctal%newmolecularLevel(subcell,:) - (0.0)*diff! (damped) update molecular levels - reduced from 0.2 to quicken convrgnce
              endif

              if(maxval(temp, mask = (.not. itransdone(1:6))) .lt. tolerance * 2.)then
                 thisOctal%molecularLevel(subcell,:) = &
                      thisOctal%newmolecularLevel(subcell,:) - 0.0 * diff! (damped) update molecular levels - reduced from 0.2 to quicken convrgnce
              endif

           else

              if (maxval(temp, mask = (.not. itransdone(1:6))) .gt. tolerance * 10.) then
                 thisOctal%molecularLevel(subcell,:) = &
                      thisOctal%newmolecularLevel(subcell,:) - (0.0)*diff! (enhanced) update molecular levels
              else
                 thisOctal%molecularLevel(subcell,:) = &
                      thisOctal%newmolecularLevel(subcell,:) - (0.0)*diff! (damped) update molecular levels - reduced from 0.2 to quicken convrgnce
              endif

              if(maxval(temp, mask = (.not. itransdone(1:6))) .lt. tolerance * 2.)then
                 thisOctal%molecularLevel(subcell,:) = &
                      thisOctal%newmolecularLevel(subcell,:) + 0.2 * diff! (damped) update molecular levels - reduced from 0.2 to quicken convrgnce
              endif
           endif

           deallocate(diff)

           l = 1 + mod(l,nVoxels)

        endif
     enddo
     deallocate(newFracChangePerLevel,temp)
   end subroutine swapPops

   recursive subroutine calculateOctalParams(grid, thisOctal, thisMolecule, deltaV, firsttime)

     use input_variables, only : useDust, iTrans

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i
     integer :: iupper, ilower
     integer, save :: ilambda = 1
     real(double) :: nmol, nlower, nupper, etaline, kappaAbs, deltaV
     real :: lambda
     logical,save :: onetime = .true.
     logical :: firsttime

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
                 thisOctal%microturb(subcell) = max(1d-8,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / &
                      (28.0 * amu)) + 0.3**2) / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence

                 thisOctal%velocity(subcell) = keplerianVelocity(subcellcentre(thisOctal,subcell), grid)
                 CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
!              else
!                 thisOctal%microturb(subcell) = max(1d-8,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / &
!                      (29.0 * amu)) + 0.3**2) / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence
!                 thisOctal%velocity(subcell) = keplerianVelocity(subcellcentre(thisOctal,subcell), grid)
!                 CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
              endif

              thisOctal%molmicroturb(subcell) = 1.d0/thisOctal%microturb(subcell)

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

     fac = min((dv * b)**2,500.d0) ! avoid floating point underflows
     phi = exp(-fac) * (b * OneOversqrtPi)

   end function phiProf

 ! Make a 'ray' with a random position - then use it to update intensity and opacity
   subroutine getRay(grid, fromOctal, fromSubcell, position, direction, ds, phi, i0, thisMolecule,fixedrays,itransdone)

     use input_variables, only : useDust

     type(MOLECULETYPE) :: thisMolecule
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal, startOctal, fromOctal

     logical, optional :: fixedrays
     logical :: stage1
     integer :: fromSubcell
     integer :: subcell
     real(double) :: ds, phi, i0(:), r
     integer :: iTrans,nTrans
     type(OCTALVECTOR) :: position, direction, currentPosition, thisPosition, thisVel, rayvel
     type(OCTALVECTOR) :: startVel, endVel, endPosition
     real(double) :: alphanu(2), alpha, alphatemp, snu, jnu
     integer :: iLower(50) , iUpper(50)
     real(double) :: dv, deltaV
     integer :: i
     real(double) :: distArray(200), tval
     integer :: nTau
     real(double) :: nLower(50), nUpper(50)
     real(double) :: dTau, etaline, kappaAbs 
     real(double), allocatable :: tau(:)
     real(double) :: dvAcrossCell, projVel, endprojVel
     real(double) :: CellEdgeInterval, phiProfVal
     real :: lambda
     integer :: ilambda
     integer :: iCount

     real,save :: r1(5) !! quasi random number generators
     logical,save :: firsttime = .true.
     logical :: quasi = .false.

     logical :: itransdone(200)
     logical :: realdust
     real(double),save :: BnuBckGrnd(200)
     real(double) :: balance(200), spontaneous(200) 
     real(double) :: nMol
     integer, parameter :: maxSamplePoints = 100
 !    real(double) :: z,t

     realdust = .false.

     if(present(fixedrays)) then
        stage1 = fixedrays
     else
        stage1 = .false.
     endif

     nTrans = thisMolecule%nTrans

     allocate(tau(1:nTrans))

     position = randomPositionInCell(fromOctal, fromsubcell)

     icount = 1
     thisOctal => fromOctal!grid%octreeRoot
     subcell = fromSubcell

     rayVel = Velocity(position, grid, thisOctal, subcell)

 !!! quasi random code!!!

     if(firsttime) then
 
        do itrans = 1, thisMolecule%nTrans
           BnuBckGrnd(itrans) = Bnu(thisMolecule%transfreq(itrans), Tcbr)
!          Bnudust(itrans) = Bnu(thisMolecule%transfreq(itrans), thisOctal%dusttemperature(subcell))
        enddo

        firsttime =.false.
        if (quasi) call sobseq(r1,-1)
     endif

     if(quasi) then     
        call sobseq(r1)    
        direction = specificUnitVector(r1(1),r1(2))
 !       r = (r1(3)+r1(4)+r1(5)) / 3.
        deltaV = 4.3 * thisOctal%microturb(subcell) * (r1(3) - 0.5d0) ! random frequency near line spectru$
     else
        call random_number(r)
 !       r = gasdev() + 0.5
        direction = randomUnitVector() ! Put this line back if you want to go back to pseudorandom
        deltaV = 4.3 * thisOctal%microturb(subcell) * (r - 0.5d0) ! random frequency near line spectrum peak. 

     endif

     deltaV = deltaV + (rayVel .dot. direction) ! transform to take account of grid velocity

     projVel = deltaV - (rayVel .dot. direction) ! transform back to velocity relative to local flow

     call distanceToCellBoundary(grid, position, direction, ds, sOctal=thisOctal)

     currentPosition = position + ds * direction

     if (inOctal(grid%octreeRoot, currentPosition)) then ! check that we're still on the grid
        thisVel = velocity(currentposition, grid) 
        endprojVel = deltaV - (thisVel.dot.direction)
        phi = phiProf((endprojVel+projvel)*0.5d0, thisOctal%molmicroturb(subcell))
     else
        phi = phiProf(projVel, thisOctal%molmicroturb(subcell)) ! if fell off grid then assume phi is unchanged
     endif

     currentPosition = position !This line is wrong because the source function does not apply locally
     ds = ds * 1.d10 ! convert from cm to torus units

     i0 = 0.d0
     tau = 0.d0

     do while(inOctal(grid%octreeRoot, currentPosition))

        call findSubcellLocal(currentPosition, thisOctal, subcell)
        call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

        startVel = velocity(currentposition, grid) 
        endPosition = currentPosition + tval * direction
        endVel = velocity(endposition, grid)

        dvAcrossCell = ((startVel - endvel).dot.direction) ! start - end
        dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))

        nTau = min(max(2, nint(dvAcrossCell * 20.d0)), maxSamplePoints) ! selects dVacrossCell as being between 0.1 and 10 (else nTau bounded by 2 and 200)
        distArray(1) = 0.d0

        CellEdgeInterval = 1.d0 / (nTau - 1.d0)

        nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)

        do itrans = 1, nTrans
           if(itransdone(itrans)) cycle

           iUpper(itrans) = thisMolecule%iTransUpper(iTrans)
           iLower(itrans) = thisMolecule%iTransLower(iTrans)

           nLower(itrans) = thisOctal%molecularLevel(subcell,iLower(itrans)) * nMol
           nUpper(itrans) = thisOctal%molecularLevel(subcell,iUpper(itrans)) * nMol
           balance(itrans)= (nLower(itrans) * thisMolecule%einsteinBlu(iTrans) - &
                            nUpper(itrans) * thisMolecule%einsteinBul(iTrans))
           spontaneous(itrans) = thisMolecule%einsteinA(iTrans) * nUpper(itrans)

        enddo

        do i = 2, nTau    

           distArray(i) = tval * dble(i-1) * cellEdgeInterval

           startOctal => thisOctal
           thisPosition = currentPosition + distArray(i)*direction
           thisVel = velocity(thisPosition, grid)
           dv = deltaV - (thisVel .dot. direction)

           PhiProfVal = phiProf(dv, thisOctal%molmicroturb(subcell))
           alphatemp = hCgsOverfourPi * phiProfVal

           do iTrans = 1, nTrans
           if(itransdone(itrans)) cycle

              alphanu(1) = alphatemp             
              alphanu(1) = alphanu(1) * balance(itrans) ! Equation 8

              if(useDust) then             
                 if(realdust) then
                    lambda = cspeed_sgl / thisMolecule%transfreq(iTrans) * 1e8 ! cm to Angstroms
                    call locate(grid%lamArray, size(grid%lamArray), lambda, ilambda)
                    call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = lambda, kappaAbs = kappaAbs)
                 else
!                    kappaAbs = thisOctal%rho(subcell) * 0.1 * thisMolecule%transfreq(itrans) * 1e-12 * 1d10 !multiplied by density !cm -> torus
                     kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10 !multiplied by density !cm -> torus
                 endif

                 alphanu(2) = kappaAbs * 1.d-10 !torus units -> cm     !kappa already multiplied by density in amr_mod
              else
                 alphanu(2) = 0.d0
              endif

              alpha = alphanu(1) + alphanu(2)
              dTau = alpha * distArray(2) * 1.d10 ! distarray(2) is interval width & optical depth, dTau = alphanu*ds - between eqs (3) and (4)
              etaLine = hCgsOverfourPi * spontaneous(itrans) 

              if(usedust) then
                 jnu = etaLine * phiProfVal + &
                      alphanu(2) * bnu(thisMolecule%transfreq(iTrans), dble(thisOctal%temperaturedust(subcell))) ! jnu, emission coefficient - equation 7
              else
                 jnu = etaLine * phiProfVal
              endif

              if (alpha .ne. 0.d0) then
                 snu = jnu/alpha ! Source function 
                 i0(iTrans) = i0(iTrans) +  exp(-tau(iTrans)) * (1.d0-exp(-dtau))*snu ! summed radiation intensity from line integral + 2nd term is local radiation field 
              else
                 snu = tiny(snu)
                 i0 = tiny(i0)
              endif

              tau(iTrans) = tau(iTrans) + dtau ! contribution to optical depth from this line integral

           enddo
        enddo
        currentPosition = currentPosition + (tval + 1.d-3*grid%halfSmallestSubcell) * direction ! FUDGE - make sure that new position is in new cell
     enddo

     do iTrans = 1, thisMolecule%nTrans
        i0(iTrans) = i0(iTrans) + BnuBckGrnd(itrans) * exp(-tau(iTrans))
     enddo

     deallocate(tau)
   end subroutine getRay

 subroutine calculateMoleculeSpectrum(grid, thisMolecule)
   use input_variables, only : itrans, nSubpixels, inc, usedust

#ifdef MPI
       include 'mpif.h'
#endif

   type(GRIDTYPE) :: grid
   type(MOLECULETYPE) :: thisMolecule
   type(OCTALVECTOR) :: unitvec, posvec, centrevec, viewvec
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

     if(grid%geometry .eq. "iras04158") call plotdiscvalues(grid, thisMolecule)

     if(usedust) call testOpticalDepth(grid, thisMolecule)

     call setObserverVectors(inc, viewvec, unitvec, posvec, centreVec, grid%octreeroot%subcellsize)

     call createimage(cube, grid, unitvec, posvec, thismolecule, itrans, nSubpixels) ! create an image using the supplied parameters (also calculate spectra)
     
     call cubeIntensityToFlux(cube, thismolecule, itrans) ! convert intensity (ergcm-2sr-1Hz-1) to flux (Wm-2Hz-1) (so that per pixel flux is correct)

     !Commented out lines are for removing background intensity - will want to do this one day

     !   call TranslateCubeIntensity(cube,1.d0*Tcbr) ! 
     !   call cubeIntensityToFlux(cube, thismolecule, itrans) ! convert intensity (ergcm-2sr-1Hz-1) to flux (Wm-2Hz-1) (so that per pixel flux is correct)
     
     call createFluxSpectra(cube, thismolecule, itrans)

     if(writeoutput) call writedatacube(cube, 'MolRT.fits')

     write(filename,'(a,a,i1,a,i1,a)') trim(cube%telescope%label),'fluxcubeJ', &
          thisMolecule%itransUpper(itrans)-1,'-',thisMolecule%itransLower(itrans)-1,'.ps/vcps'
     call plotDataCube(cube, filename, plotflux=.true.)

     !  call plotDataCube(cube, 'subpixels.ps/vcps',withSpec=.False.,GotPixels=.true.) ! plot the image 
     stop

    end subroutine calculateMoleculeSpectrum

    ! sample level populations at logarithmically spaced annuli
   subroutine dumpResults(grid, thisMolecule, nRay, iter, convtestarray)
     use input_variables, only : rinner, router
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     integer :: nr
     real(double) :: r, ang
     integer :: i
     real(double) :: pops(8), fracChange(8), convtestarray(:,:,:)
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     integer :: j
     real(double) :: x, z
     type(OCTALVECTOR) :: posvec
     integer :: nRay, iter
     character(len=30) :: resultfile, resultfile2
     integer :: minlevels

     minlevels = min(thismolecule%nlevels, 8)

     write(resultfile,'(a,I7.7)') "results.", nRay
     write(resultfile2,'(a,I2)') "fracChangeGraph.", iter

     nr = 200 ! number of rays used to sample distribution
     open(31,file=resultfile,status="unknown",form="formatted")
     open(32,file=resultfile2,status="unknown",form="formatted")

     do i = 1, nr
        r = log10(rinner) + dble(i-1)/dble(nr-1)*(log10(router) - log10(rinner))
        r = 10.d0**r ! log(radius) -> radius
        pops = 0.d0
        fracChange = 0.d0
        do j = 1 , nr
 !          call random_number(ang)
           ang = (j - 1.) / nr
           ang = ang * pi
           z = r*cos(ang)
           x = r*sin(ang)
           posVec = OCTALVECTOR(x, 0.d0, z) ! get put in octal on grid
           thisOctal => grid%octreeroot
           call findSubcellLocal(posVec, thisOctal,subcell) 
           pops = pops + thisOctal%molecularLevel(subcell,1:minlevels) ! interested in first 8 levels
           fracChange = fracChange + abs((thisOctal%molecularLevel(subcell,1:minlevels) - &
                        thisOctal%oldmolecularLevel(subcell,1:minlevels)) / thisOctal%molecularlevel(subcell,1:minlevels))
        enddo
        pops = pops / real(nr) ! normalised level population at the 20 random positions 
        fracChange = fracChange / real(nr)

        if(iter .gt. 0) convtestarray(iter,i,:) = fracChange
        write(31,'(es11.5e2,4x,8(es12.6e2,tr2))') real(r*1.d10), real(pops(1:minlevels)) 
        write(32,'(es11.5e2,4x,8(f7.5,tr2))') real(r*1.d10), real(fracChange(1:minlevels))
     enddo
     close(31)
     close(32)
   end subroutine dumpResults

 ! Does a lot of work - do more rays whilst problem not converged -            
   subroutine molecularLoop(grid, thisMolecule)

     use input_variables, only : blockhandout, tolerance, debug, geometry, lucyfilenamein, openlucy,&
          usedust, rinner, router, readmol
     use messages_mod, only : myRankIsZero
#ifdef MPI
     include 'mpif.h'
#endif

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     type(OCTALVECTOR) :: position, direction
     integer :: nOctal, iOctal, subcell
     real(double), allocatable :: ds(:), phi(:), i0(:,:)
     integer :: nRay !, tempray
     integer :: previousnRay = 100
     type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
     type(OCTAL), pointer :: thisOctal
     integer, parameter :: maxIter = 1000, maxRay = 10000000
     logical :: popsConverged, gridConverged, gridConvergedTest
     character(len=200) :: message
     integer :: iRay, iTrans, iter,i, grand_iter 
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

      character(len=40) :: filename

      real(double) :: convtestarray(200,200,8)
      real(double) :: r(100)
      logical, save :: itransdone(200)
      real(double) :: fixValMax = 1.
      real(double) :: fixValMin = 1e-4

      logical :: restart = .false.
      logical :: plotlevels = .true.
      
      integer :: ntrans , lst !least significant transition
      
      real :: tol
      real(double) :: dummy, tau, tauarray(8)    
      integer :: mintrans ! for dumpresults

 !        real(double) :: r,rarray(4100000),sarray(4100000),tarray(4100000),phiprofval,deltav,dv
 !        real :: uarray(4100000)
 !        type(VSL_STREAM_STATE) :: stream
 !        integer :: status

 ! blockhandout must be off for fixed ray case, otherwise setting the
 ! seed is not enough to ensure the same directions are done for
 ! each cell every iteration

      mintrans = min(thismolecule%ntrans, 8) ! converge the first 8 levels
      allocate(avgfracChange(mintrans,2))

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

 18   format(a,tr3,8(a,tr5),a,6(tr1,a))
      open(138,file="fracChanges.dat",status="unknown",form="formatted")
      write(138,18) "nRays","J=0","J=1","J=2","J=3","J=4","J=5","J=6","J=7","Fixed","1% Conv","2.5% Conv","5% Conv","Max","Avg"
      close(138)

      open(139,file="avgChange.dat",status="replace",form="formatted")
      open(140,file="avgRMSChange.dat",status="replace",form="formatted")
      open(142,file="tau.dat",status="replace",form="formatted")
      open(143,file="tauhalf.dat",status="replace",form="formatted")
      close(139)
      close(140)
      close(142)
      close(143)

 !     if(writeoutput) then
 !        write(*,*) "Lucy turned on?",lucyradiativeeq
 !        write(*,*) "Geometry? ",geometry
 !     endif

      if(openlucy) then

         if(writeoutput) then
            write(message,*) "Reading in lucy temp files: ",lucyfilenamein
            call writeinfo(message,FORINFO)
         endif

         if(writeoutput) write(*,*) "READ"
         call readAmrGrid(lucyfilenamein,.false.,grid)

        if(writeoutput) write(*,*) "OK"
         if(writeoutput) then
            write(message,*) "Done! Plotting"
            call writeinfo(message,TRIVIAL)
            call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
                 "lucytemps.ps/vcps", .true., .false., &
                 width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.) 
         endif

      endif

      allocate(octalArray(grid%nOctals))
      call getOctalArray(grid%octreeRoot,octalArray, nOctal)

      if(restart) then
         if(writeoutput) then
            write(message,*) "Reading in previous grid"
            call writeinfo(message,TRIVIAL)
         endif

         call readAMRgrid("molecular_tmp.grid",.false.,grid)

         call allocateMolecularLevels(grid, grid%octreeRoot, thisMolecule, .true., .true.)
      else
         if(writeoutput) then
            write(message,*) "Allocating and initialising molecular levels"
            call writeinfo(message,TRIVIAL)
         endif
         call AllocateMolecularLevels(grid, grid%octreeRoot, thisMolecule, .false., .true.)
         if(writeoutput) then
            write(message,*) "Done!"
            call writeinfo(message,TRIVIAL)
         endif

         if(myrankiszero) call writeAMRgrid("molecular_lte.grid",.false.,grid)
      endif

      if(writeoutput) then

         do i=1,4

            write(filename,'(a,i1,a)') "./J/J=",i-1,"lte.ps/vcps"
            call plot_AMR_values(grid, "J", "x-z", real(grid%octreeRoot%centre%y), &
                 filename, .true., .false., &
                 width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false., &
                 ilam = i, fixValMin=fixValMin, fixValMax=fixValMax)
         enddo
         
         call plot_AMR_values(grid, "molAbundance", "x-z", real(grid%octreeRoot%centre%y), &
              "molAbundance.ps/vcps", .true., .false., &
              width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.,&
              fixValMin=1d-10, fixValMax=1d-3)          

         call plot_AMR_values(grid, "molAbundance", "x-z", real(grid%octreeRoot%centre%y), &
              "molAbundanceZoom.ps/vcps", .true., .false., &
              width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.,&
              fixValMin=1d-13, fixValMax=1d-3, boxfac = 0.001)          

      endif

      if(writeoutput) then

         write(*,*) "Dumping subcell parameters"
         i = 0
         open(11,file="celltest.dat",status='unknown',form='formatted')

         do iOctal = 1, SIZE(octalArray)
            thisOctal => octalArray(iOctal)%content
            do subcell = 1, thisOctal%maxChildren
               if (.not.thisOctal%hasChild(subcell)) then

                  write(11,'(6(a,es10.4,tr1))') &
                       "temp",thisOctal%temperature(subcell),&
                       "tempgas",thisOctal%temperaturegas(subcell),&
                       "tempdust",thisOctal%temperaturedust(subcell),&
                       "microturb",thisOctal%microturb(subcell),&
                       "nh2",thisOctal%nh2(subcell)

               endif
            enddo
         enddo

         close(11)        

      endif

      if(.not. grid%octreeroot%threed) call dumpresults(grid, thisMolecule, 0, grand_iter, convtestarray) ! find radial pops on final grid     
     
      if(usedust) then 
!         call continuumIntensityAlongRay(octalvector(1.d0,0.d0,0.d0),octalvector(1.0,1d-20,1d-20), grid, 1e-4, 0.d0, dummy, &
!              tau, .true.)
         call continuumIntensityAlongRay(octalvector(1d10,0.d0,0.d0),octalvector(-1.d0,1d-20,1d-20), grid, 1e4, 0.d0, dummy, &
              tau, .true.)
         if(writeoutput) write(*,*) "TAU", tau
      endif

      if(Writeoutput) write(*,*) "COUNTING VOXELS"
      call countVoxels(grid%octreeRoot,nOctal,nVoxels)

      nRay = 100 ! number of rays used to establish estimate of jnu and pops

      allocate(oldPops(1:thisMolecule%nLevels))

      call random_seed

     !       if (myRankIsZero) &
 !    call writeAmrGrid("molecular_tmp.grid",.false.,grid)

     call random_seed(size=iSize)
     allocate(iSeed(1:iSize))
     call random_seed(get=iSeed)
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
      call MPI_BCAST(iSeed, iSize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif


 !---------------------
 !TESTING!
 !    write(*,*) "PHIPROF TEST"
 !    status = vslnewstream(stream, VSL_BRNG_MCG31, 1 )

 !    call tune(6, 'gaussRANDs')  ! start a stopwatchtune
 !    status = vsrnggaussian(VSL_METHOD_SGAUSSIAN_ICDF, stream, 4100000, uarray, 0., 1.)
 !    call tune(6, 'gaussRANDs')  ! stop a stopwatchtune

 !    call tune(6, 'uniformRAND')  ! start a stopwatchtune
 !    status = vsrnguniform(VSL_METHOD_SUNIFORM_STD, stream, 4100000, uarray, 0., 1.)
 !    call tune(6, 'uniformRAND')  ! stop a stopwatchtune

 !    call tune(6, 'gaussRAND')  ! start a stopwatchtune
 !    status = vdrnggaussian(VSL_METHOD_DGAUSSIAN_ICDF, stream, 4100000, rarray, 0.d0, 1.d0)
 !    call tune(6, 'gaussRAND')  ! stop a stopwatchtune

 !    call tune(6, 'gaussRAND')  ! start a stopwatchtune
 !    status = vdrnggaussian(VSL_METHOD_DGAUSSIAN_BOXMULLER, stream, 4100000, rarray, 0.d0, 1.d0)
 !    call tune(6, 'gaussRAND')  ! stop a stopwatchtune


 !    call tune(6, 'RAND')  ! start a stopwatchtune
 !    do i = 1, 4100000
 !       call random_number(r)
 !       sarray(i) = r
 !    enddo
 !    call tune(6, 'RAND')  ! stop a stopwatchtune

 !    call tune(6, 'numrecRAND')  ! stop a stopwatchtune
 !    do i = 1, 4100000
 !       tarray(i) = gasdev()
 !    enddo
 !    call tune(6, 'numrecRAND')  ! stop a stopwatchtune

 !    open(7,file="phiprofs.dat",status="unknown",form="formatted")
 !    do i = 1, 10000


 !       deltaV = tarray(i)
 !       deltaV = 5e-7 * deltaV
 !       dv = -1.d0 * deltaV
 !       phiProfval = phiProf(dv, 5d-7)

 !       write(7,*) dv,phiprofval

 !    enddo

 !close(7)
 !stop
 !---------------------

     do iStage = 1, 2 ! fixed rays to reduce variance between cells or 2)  random rays to ensure sufficient spatial sampling

        if (iStage == 1) then
           fixedRays = .true.
           nRay = 100
        else
           fixedRays = .false.
           nRay = 100
           gridConvergedTest = .false.
        endif

        gridConverged = .false.

        do while (.not.gridConverged)

         grand_iter = grand_iter + 1

         if(writeoutput) then 
            write(*,*) "Iteration ",grand_iter
            write(message,*) "Done ",nray," rays"
            call tune(6, message)  ! start a stopwatch
         endif

         if(Writeoutput .and. plotlevels) then

            do i=1,4
               if(grand_iter .lt. 10) then
                  write(filename,'(a,i1,a,i1,a)') "./J/J=",i-1,"-Iter0",grand_iter,".ps/vcps"
               else
                  write(filename,'(a,i1,a,i2,a)') "./J/J=",i-1,"-Iter",grand_iter,".ps/vcps"
               endif

               call plot_AMR_values(grid, "J", "x-z", real(grid%octreeRoot%centre%y), &
                    filename, .true., .false., &
                    width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false., &
                    ilam = i, boxfac = 0.7, fixValMin=fixValMin, fixValMax=fixValMax)
            enddo

            do i=1,4
               if(grand_iter .lt. 10) then
                  write(filename,'(a,i1,a,i1,a)') "./J/J=",i-1,"-Iter0",grand_iter,"zoom.ps/vcps"
               else
                  write(filename,'(a,i1,a,i2,a)') "./J/J=",i-1,"-Iter",grand_iter,"zoom.ps/vcps"
               endif

               call plot_AMR_values(grid, "J", "x-z", real(grid%octreeRoot%centre%y), &
                    filename, .true., .false., &
                    width_3rd_dim=real(grid%octreeRoot%subcellsize), show_value_3rd_dim=.false., &
                    ilam = i, boxfac = 0.001, fixValMin=fixValMin, fixValMax=fixValMax)
            enddo

            do i=1,4
               if(grand_iter .lt. 10) then
                  write(filename,'(a,i1,a,i1,a)') "./J/J=",i-1,"-Iter0",grand_iter,"zoomlittle.ps/vcps"
               else
                  write(filename,'(a,i1,a,i2,a)') "./J/J=",i-1,"-Iter",grand_iter,"zoomlittle.ps/vcps"
               endif

               call plot_AMR_values(grid, "J", "x-z", real(grid%octreeRoot%centre%y), &
                    filename, .true., .false., &
                    width_3rd_dim=real(grid%octreeRoot%subcellsize) ,  show_value_3rd_dim=.false., &
                    ilam = i, boxfac = 0.009, fixValMin=fixValMin, fixValMax=fixValMax)

            enddo

         endif

   allocate(ds(1:nRay))
   allocate(phi(1:nRay))
   allocate(i0(1:thisMolecule%nTrans, 1:nRay))

   if (fixedRays) then
      call random_seed(put=iseed)   ! same seed for fixed rays

   else
      call random_seed
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

 ! iterate over all octals, all rays, solving the system self-consistently
           if(writeoutput) then
              write(message,*) "Getray"
              call tune(6, message)  ! start a stopwatch
           endif

           do iOctal = ioctal_beg, ioctal_end

              if (debug .and. writeoutput) then
                 write(message,*) iOctal,ioctal_beg,ioctal_end
                 call writeInfo(message,TRIVIAL)
              endif

              thisOctal => octalArray(iOctal)%content
              do subcell = 1, thisOctal%maxChildren

                 if (.not.thisOctal%hasChild(subcell)) then
 !                   if(fixedrays) call sobseq(r1,-1)                      

                    do iRay = 1, nRay
                       itransdone = .false.
                       
                       call getRay(grid, thisOctal, subcell, position, direction, &
                            ds(iRay), phi(iRay), i0(1:thisMolecule%nTrans,iRay), &
                            thisMolecule,fixedrays,itransdone) ! does the hard work - populates i0 etc

                    enddo

                    iter = 0
                    popsConverged = .false.

                       thisOctal%newMolecularLevel(subcell,:) = thisOctal%molecularLevel(subcell,:)
                       thisOctal%oldMolecularLevel(subcell,:) = thisOctal%molecularLevel(subcell,:)

                    do while (.not.popsConverged)
                       iter = iter + 1

                       oldpops = thisOctal%newmolecularLevel(subcell,1:thisMolecule%nLevels) ! retain old pops before calculating new one

                       do iTrans = 1, thisMolecule%nTrans!/2

                          if(.not. itransdone(itrans)) then
                             call calculateJbar(grid, thisOctal, subcell, thisMolecule, nRay, ds(1:nRay), &
                                  phi(1:nRay), i0(iTrans,1:nRay), iTrans, thisOctal%jnu(subcell,iTrans), &
                                  thisOctal%newMolecularLevel(subcell,1:thisMolecule%nLevels)) ! calculate updated Jbar
                          endif

                       enddo

                       call solveLevels(thisOctal%newMolecularLevel(subcell,1:thisMolecule%nLevels), &
                            thisOctal%jnu(subcell,1:thisMolecule%nTrans),  &
                            dble(thisOctal%temperature(subcell)), thisMolecule, thisOctal%nh2(subcell))

                       fac = abs(maxval((thisOctal%newMolecularLevel(subcell,1:mintrans) - oldpops(1:mintrans)) &
                             / oldpops(1:mintrans))) ! convergence criterion ! 6 or 8?
                       if (fac < 1.d-6) popsConverged = .true.

                       if (iter == maxIter) then
                          popsConverged = .true.
                          call writeWarning("Maximum number of iterations reached in pop solver")
                       endif
                    enddo
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
      do i = 1, thisMolecule%nLevels
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
      nTrans = min(thisMolecule%nTrans,12)
      lst = min(thisMolecule%nTrans,8) ! least significant transition
      
      allocate(maxFracChangePerLevel(lst))
      allocate(ConvergenceCounter(4,nTrans+1))
      maxFracChangePerLevel = -1.d30 ! negative to ensure comparison first time round

      convergenceCounter = 0
      avgFracChange = 0.d0

      call swapPops(grid%octreeRoot, maxFracChangePerLevel, avgFracChange, &
           convergenceCounter, grand_iter, nVoxels, iTransdone, fixedrays) ! compares level populations between this and previous levels 

      maxavgFracChange = maxval(avgFracChange(:,1))
      maxRMSFracChange = maxval(avgFracChange(:,2))

      maxavgtrans = maxloc(avgFracChange(:,1)) - 1
      maxRMStrans = maxloc(avgFracChange(:,2)) - 1

      maxFracChange = MAXVAL(maxFracChangePerLevel(1:lst)) ! Largest change of any level < 6 in any voxel

      if(writeoutput) then
         write(message,'(a,1x,f9.5,1x,a,1x,f5.3,1x,a,2x,l1,1x,a,1x,i6,1x)') &
              "Maximum fractional change this iteration ", maxFracChange, "tolerance", tolerance, "fixed rays", fixedrays, &
              "nray", nray
         write(*,'(a,f9.5,a,i1)') "  Average fractional change this iteration  ", &
              maxavgFracChange/real(nVoxels)," in level ",maxavgTrans
         write(*,'(a,f9.5,a,i1)') "  RMS fractional change this iteration      ", &
              sqrt(maxRMSFracChange/real(nVoxels))," in level ",maxRMSTrans
         write(*,'(a,f9.5)') "  Std Dev                                   ", &
              sqrt(maxRMSFracChange/real(nVoxels)-(maxavgFracChange/real(nVoxels))**2)

         call writeInfo(message,FORINFO)

         do i=3,1,-1
            if(i .lt. 3) write(message,'(a,f5.3,a, 13(f5.3,2x))') "Individual levels converged @ ",tolerance * i," | ",&
                 real(convergenceCounter(i,1:12))/real(nVoxels) ! Fraction of first 8 levels converged at 1,2,5 * tolerance respectively
            if(i .eq. 3) write(message,'(a,f5.3,a, 13(f5.3,2x))') "Individual levels converged @ ",tolerance * (i+2)," | ", &
                 real(convergenceCounter(i,1:12))/real(nVoxels) ! Fraction of first 8 levels converged at 1,2,5 * tolerance respectively
            call writeInfo(message,FORINFO)         
         enddo

         write(*,*) ""
         write(message,'(a,f5.3,a, 13(f5.3,2x))') "Individual levels converged @ ",tolerance * 0.5," | ", &
              real(convergenceCounter(4,1:12))/real(nVoxels) ! Fraction of first 8 levels converged at 1,2,5 * tolerance respectively
         call writeInfo(message,FORINFO)

         open(138,file="fracChanges.dat",position="append",status="unknown")
 19      format(i6,tr3,8(f7.5,1x),3x,l1,tr3,3(f5.3,tr3),f7.5,tr3,f7.5)
         write(138,19) nray, maxFracChangePerLevel(1:lst), fixedrays, real(convergenceCounter(1:3,ntrans+1))/real(nVoxels),  &
              maxFracChange, maxavgFracChange/real(nVoxels)

         close(138)

         open(139,file="avgChange.dat",position="append",status="unknown")
 20      format(i2,tr3,i6,tr3,8(f7.5,1x))
         write(139,20) grand_iter, nray, avgFracChange(1:lst,1)/real(nvoxels)

         close(139)

         open(140,file="avgRMSchange.dat",position="append",status="unknown")
         write(140,20) grand_iter, nray, sqrt(avgFracChange(1:lst,2)/real(nvoxels))

         close(140)

      endif

        deallocate(maxFracChangePerLevel)
        deallocate(ConvergenceCounter)

 ! If you think it's converged then test with the same number of rays to make sure

        if(fixedrays) then 
           tol = tolerance * 0.5
        else
           tol = tolerance
        endif

           if(sqrt(maxRMSFracChange/real(nVoxels)) < tol) then
              if (gridConvergedTest) gridConverged = .true.
              gridConvergedTest = .true.
           else
              gridConvergedTest = .false.
              gridConverged = .false.
           endif

           if (myRankIsZero) then
              call writeAmrGrid("molecular_tmp.grid",.false.,grid)
              
              write(*,*) "Dumping results"

              if(.not. grid%octreeroot%threed) call dumpresults(grid, thisMolecule, nRay, grand_iter, convtestarray) ! find radial pops on final grid     

              if(debug .and. readmol) then
                 open(142,file="tau.dat",status="unknown",position="append")
                 open(143,file="tauhalf.dat",status="unknown",position="append")

                 do itrans = 1, lst
                    call intensityAlongRay(octalvector(1d-20,1d-20,1d10),octalvector(-1d-20,-1d-20,-1.d0), grid, thisMolecule, &
                         iTrans, 0.d0, dummy, tau, .true.)
                    tauarray(itrans) = tau
                 enddo

                 write(142,'(i2,tr3,8(f9.4,tr3))') grand_iter, tauarray

                 do itrans = 1, lst
                    call intensityAlongRay(octalvector(0.d0,0.d0,0.d0),octalvector(1d-20,1d-20,1.d0), grid, thisMolecule, iTrans, &
                         0.d0, dummy, tau, .true.)
                    tauarray(itrans) = tau
                 enddo

                 write(143,'(i2,tr3,8(f9.4,tr3))') grand_iter, tauarray
                 close(142)
                 close(143)
              endif
           endif

#ifdef MPI
           deallocate(octalsBelongRank)
#endif
           deallocate(ds, phi, i0)

           if(writeoutput) then
              write(message,*) "Done ",nray," rays"
              call tune(6, message)  ! stop a stopwatch
           endif

           if (.not.gridConvergedTest) then
              if (.not.gridConverged) then               
                 if (.not.fixedRays) nRay = nRay * 2 !double number of rays if convergence criterion not met and not using fixed rays - revise!!! Can get away with estimation?
!                 if (.not.fixedRays) then 
!                    tempRay = nRay + previousnRay
!                    previousnRay = nRay
!                    nRay = tempray
!                 endif
                 write(message,*) "Trying ",nRay," Rays"
                 call writeInfo(message,FORINFO)
                 if (grid%geometry .eq. 'molebench' .and. nray .gt. 500000) then 
                    if(writeoutput) write(*,*) "Molebench Test Ended, Exiting..."
                    gridconverged = .true.
                 endif
                 
              endif
           else
              if(writeoutput .and. .not. gridconverged) write(*,*) "Doing all rays again"
              itransdone = .false.
              nRay = nRay
           endif

           if (nRay > maxRay) then
              nRay = maxRay  ! stop when it's not practical to do more rays
              call writeWarning("Maximum number of rays exceeded - capping")
           endif

        enddo

        if(writeoutput .and. (geometry .eq. 'molebench')) then
           open(141,file="convergenceprofile.dat",status="unknown",form="formatted")

           do i=1,100
              r(i) = log10(rinner) + dble(i-1)/dble(100-1)*((log10(router) - log10(rinner))) ! log(radius) -> radius
           enddo


           write(141,'(tr3,100(es8.3e1,tr1))') 10.**r
           do i = 1,grand_iter
              write(141,'(i2,tr1,100(f8.6,tr1))') i, convtestarray(i,:,1)
           enddo
           close(141)
        endif
     enddo

     write(*,*) "mole loop done."
     stop
   end subroutine molecularLoop

   subroutine calculateJbar(grid, thisOctal, subcell, thisMolecule, nRay, ds, phi, i0, iTrans, jbar, nPops)

     use input_variables, only : useDust
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(MOLECULETYPE) :: thisMolecule
     integer :: nRay
     real(double) :: ds(:), phi(:), i0(:), nPops(:)
     integer :: iTrans
     real(double) :: jbar
     integer :: iRay
     real(double) :: nLower, nUpper, nMol
     real(double) :: jBarInternal, jBarExternal
     real(double) :: alphanu(2), jnu, etaline, kappaAbs, alpha
     integer :: iUpper, iLower
     real(double) :: tau, opticaldepth, snu, sumPhi
     real :: lambda
     integer :: ilambda
     logical :: realdust = .false.
     
     jBarExternal = 0.d0
     jBarInternal = 0.d0

     jnu = 0.d0

     if(useDust) then
        lambda = cspeed_sgl / thisMolecule%transfreq(iTrans) * 1d8
        call locate(grid%lamArray, size(grid%lamArray), lambda, ilambda)
     endif

     iUpper = thisMolecule%iTransUpper(iTrans)
     iLower = thisMolecule%iTransLower(iTrans)
     ! commented elsewhere in the code
     sumPhi = 0.d0

     nMol = thisOctal%molAbundance(subcell) * thisOctal%nh2(subcell)

        
     do iRay = 1, nRay

        alphanu(1) = hCgsOverfourPi!*thisMolecule%transFreq(iTrans) 
        nLower = nPops(iLower) * nMol
        nUpper = nPops(iUpper) * nMol

        etaLine = hCgsOverFourPi * thisMolecule%einsteinA(iTrans)! * thisMolecule%transFreq(iTrans)
        etaLine = etaLine * nUpper

        alphanu(1) = alphanu(1) * (nLower * thisMolecule%einsteinBlu(iTrans) - &
             nUpper * thisMolecule%einsteinBul(iTrans)) * phi(iray)!/thisMolecule%transFreq(iTrans)

        jnu = (etaLine) * phi(iRay)!/thisMolecule%transFreq(iTrans)
 
        if(useDust) then
           if(realdust) then
              call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = lambda, kappaAbs = kappaAbs)
           else
              kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10 !multiplied by density !cm -> torus
           endif

           alphanu(2) = kappaAbs * 1.d-10 !* thisOctal%rho(subcell) - already done
           jnu = jnu + alphanu(2) * bnu(thisMolecule%transfreq(iTrans), dble(thisOctal%temperaturedust(subcell)))
 !       write(*,*) "alpha dust",alpha(2),"bnu", bnu(thisMolecule%transfreq(iTrans), dble(thisOctal%temperature(subcell)))
        else
           alphanu(2) = 0.d0
        endif

        alpha = alphanu(1) + alphanu(2)

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

 !    jbar = (jBarExternal + jBarInternal)/dble(nRay)
      jbar = (jBarExternal + jBarInternal)/sumPhi

    end subroutine calculateJbar

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
!sub intensity
   subroutine intensityAlongRay(position, direction, grid, thisMolecule, iTrans, deltaV,i0,tau,tautest)

     use input_variables, only : useDust
     type(OCTALVECTOR) :: position, direction, dsvector
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0
     type(OCTAL), pointer :: thisOctal
     integer :: subcell
     type(OCTALVECTOR) :: currentPosition, thisPosition, endPosition
     type(OCTALVECTOR) :: startVel, endVel, thisVel, veldiff
     real(double) :: alphanu(2), jnu, snu
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
 !       write(*,*) iupper,ilower
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
        alphanu(1) = thisOctal%molcellparam(subcell,6)
        alphanu(2) = thisOctal%molcellparam(subcell,7)
        dustjnu = thisOctal%molcellparam(subcell,8)

 !       write(*,*) dvAcrossCell, OneOvernTauMinusOne, startvel, endvel, thisOctal%microturb(subcell)                                                                                                                                                                                

!        if(usedust) then
!           if(realdust) then
!              call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, lambda = lambda, kappaAbs = kappaAbs)
!           else
!             kappaAbs = thisOctal%rho(subcell) * 0.1 * thisMolecule%transfreq(itrans) * 1e-12 * 1d20 !multiplied by density !cm -> torus
!              kappaAbs = thisOctal%rho(subcell) * DustModel1(thisMolecule%transfreq(itrans)) * 1d10 !multiplied by density !cm -> torus
!           endif
           
!        else
!           alphanu(2) = 0.d0
!        endif

        ds = tval * OneOvernTauMinusOne
        dsvector = ds * direction
        thisPosition = currentPosition

        if(grid%geometry .eq. 'IRAS04158') then
           startVel = keplerianVelocity(currentposition, grid)
           endPosition = currentPosition + tval * direction
           endVel = keplerianVelocity(currentposition, grid)
        else
           startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell)
           endPosition = currentPosition + tval * direction
           endVel = amrGridVelocity(grid%octreeRoot, endPosition, startOctal = thisOctal, actualSubcell = subcell)
        endif

        Veldiff = endVel - startVel

        dvAcrossCell = (veldiff.dot.direction)
        dvAcrossCell = abs(dvAcrossCell * thisOctal%molmicroturb(subcell))

        nTau = min(max(2, nint(dvAcrossCell * 5.d0)), 2000) ! ensure good resolution / 5 chosen as its the magic number!
        
!        if(ntau .gt. 19999) write(*,*) ntau, dvacrosscell, modulus(thisposition) * 1e10 / autocm

        distArray(1) = 0.d0
        OneOvernTauMinusOne = 1.d0/(nTau - 1.d0)

        do i = 2, nTau 

           thisPosition = thisPosition + dsvector

!           if(lowvelgrad) then 
 !        if(grid%geometry .eq. 'IRAS04158') then
 !           thisVel = keplerianVelocity(thisPosition, grid)
 !        else
           thisvel = startvel + (i-1) * OneOvernTauMinusOne * Veldiff
 !        endif
!           else
!              startOctal => thisOctal
!              thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell)
!           endif

!              write(*,*) modulus(altthisvel-thisvel),thisvel%x-altthisvel%x,thisvel%y-altthisvel%y,thisvel%z-altthisvel%z

           dv = (thisVel .dot. direction) - deltaV

           phiProfval = phiProf(dv, thisOctal%molmicroturb(subcell))
           alphanu(1) = alphanu(1) * phiprofval

           alpha = alphanu(1) + alphanu(2)
           dTau = alpha * ds * 1.d10

           jnu = etaLine * phiProfVal

           if(useDust) jnu = jnu + dustjnu

           if (alpha .ne. 0.d0) then

              snu = jnu/alpha
 !             snu = thisOctal%bnu(subcell,itrans)
!              if((i0 .eq. 0.) .or. tau .lt. 10. .or. exp(-tau) * (1.d0-exp(-dtau))*snu .gt. i0*(1d-5)) then
                 i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu

!                 if(isnan(i0)) then
!                    write(*,*) icount, tau, dtau, alpha, jnu, ds, alphanu(1), alphanu(2), phiprofval, thisOctal%microturb(subcell), thisOctal%molmicroturb(subcell), dv, (dv*thisOctal%molmicroturb(subcell))**2
!                    stop
!                    endif
!                 if(dotautest .and. i .eq. ntau) write(*,*) modulus(subcellcentre(thisOctal, subcell)), i0, &
!                                                            tau, dtau, thisOctal%temperature(subcell)
!              else
!                 write(*,*) "i0", i0, "tau", tau, "skipping"
!                 goto 666
!              endif
           else
!              if(dotautest .and. i .eq. ntau) write(*,*) modulus(subcellcentre(thisOctal, subcell)), i0, &
!                   tau, alpha, snu, jnu

              snu = tiny(snu)
              i0 = i0 + tiny(i0)
           endif
           tau = tau + dtau
        enddo
        
        currentPosition = currentPosition + (tval + 1d-3 * grid%halfSmallestSubcell) * direction

!        if((grid%geometry .eq. 'IRAS04158') .and. (modulus(currentposition) < 0.5*rinner)) then
!           currentPosition = currentPosition + (rinner/(OCTALVECTOR(1d-20,1d-20,1.d0) .dot. direction)) * direction
!        endif
        
!        if(debug) write(*,*) i0,tau
     enddo

 666 continue

     i0 = i0 + bnuBckGrnd * exp(-tau) ! from far side                                                                                                                                                                                                                                
 !    if(debug) write(*,*) i0,tau
 !    i0 = i0 - bnu(thisMolecule%transFreq(iTrans), Tcbr)                                                                                                                                                                                                                            

 !    if (debug) write(*,*) "bound condition",i0,bnu(thisMolecule%transFreq(iTrans), Tcbr) * exp(-tau)                                                                                                                                                                               

     ! convert to brightness T                                                                                                                                                                                                                                                       

 !    i0 = i0 * (cSpeed**2 / (2.d0 * thisMolecule%transFreq(iTrans)**2 * kerg))                                                                                                                                                                                                      

   end subroutine intensityAlongRay

   subroutine lteintensityAlongRay(position, direction, grid, thisMolecule, iTrans, deltaV,i0,tau,tautest)

     use input_variables, only : useDust
     type(OCTALVECTOR) :: position, direction, dsvector
     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     real(double) :: disttoGrid
     integer :: itrans
     real(double) :: nMol
     real(double), intent(out) :: i0
     type(OCTAL), pointer :: thisOctal, startOctal
     integer :: subcell
     type(OCTALVECTOR) :: currentPosition, thisPosition, thisVel
     type(OCTALVECTOR) :: rayVel, startVel, endVel, endPosition
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
     rayVel = OCTALVECTOR(0.d0, 0.d0, 0.d0)

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
!           currentPosition = currentPosition + (rinner/(OCTALVECTOR(1d-20,1d-20,1.d0) .dot. direction)) * direction
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

     type(OCTALVECTOR) :: position, direction
     type(GRIDTYPE) :: grid
     real(double) :: disttoGrid
     real(double), intent(out) :: i0
     type(OCTAL), pointer :: thisOctal, startOctal
     integer :: subcell
     type(OCTALVECTOR) :: currentPosition, thisPosition, thisVel
     type(OCTALVECTOR) :: startVel, endVel, endPosition
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
        write(*,*) "inside"
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

        if(grid%geometry .eq. 'IRAS04158') then
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

   function makeImageGrid(cube, unitvec, posvec, grid, thisMolecule, iTrans, deltaV, nsubpixels) result (imagegrid)

     use input_variables, only : npixels, imageside

     type(GRIDTYPE) :: grid
     type(MOLECULETYPE) :: thisMolecule
     integer :: itrans
     type(OCTALVECTOR) :: unitvec, viewVec, posvec
     type(OCTALVECTOR) :: imagebasis(2), pixelcorner
     integer :: nsubpixels, subpixels
     real :: imagegrid(npixels,npixels)
     integer :: ipixels, jpixels
     real(double) :: pixelside
     real(double) :: deltaV
     type(datacube) :: cube
     integer :: index(2)

     logical, save :: firsttime = .true.

     if(firsttime) then

        imagebasis(2) = unitVec .cross. OCTALVECTOR(1d-20,1d-20,1d0) ! gridvector
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
     type(OCTALVECTOR) :: viewVec
     real(double) :: i0, opticaldepth
     real(double) :: totalPixelIntensity, oldTotalPixelIntensity
     type(OCTALVECTOR) :: imagebasis(2), pixelbasis(2), pixelcorner, newposvec

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
   type(OCTALVECTOR) :: viewVec
   real(double) :: i0, opticaldepth

   type(OCTALVECTOR) :: imagebasis(2), pixelbasis(2), pixelcorner, rayposition
   
   integer :: subpixels, minrays
   integer :: i, iray
   integer :: index(2)

   real(double) :: pixelside
   real(double) :: avgIntensityNew, avgIntensityOld
   real(double) :: varIntensityNew, varIntensityOld
   real(double) :: PixelIntensity
   real :: rtemp(2)
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
     type(OCTALVECTOR) :: unitvec, posvec
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
     real(double) :: weightedfluxmap(npixels,npixels)
     real(double) :: weightedfluxsum, weightedflux
     real(double) :: fineweightedfluxmap(npixels,npixels)
     real(double) :: fineweightedfluxsum, fineweightedflux, weight(npixels,npixels)

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
              call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule, deltaV,.true.)
           else
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
 ! type(OCTALVECTOR), allocatable :: SubcellCentreArray(:), GridArray(:), PixelPositionArray(:,:)
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
 ! type(OCTALVECTOR), allocatable :: SubcellCentreArray(:), GridArray(:), PixelPositionArray(:,:)
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
                    tArray(nTemps) = thisOctal%newMolecularLevel(isubcell, iLevel)
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
                  thisOctal%newMolecularLevel(isubcell, iLevel) = tArray(nTemps) 
               endif
           end do
        end do
      end subroutine unpackMoleLevel
#endif

         SUBROUTINE sobseq(x,init)
         USE nrtype; USE nrutil, ONLY : nrerror
         IMPLICIT NONE
         REAL(SP), DIMENSION(:), INTENT(OUT) :: x
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

       type(OCTALVECTOR) :: unitvec, viewvec, posvec, centreVec
       real(double) :: farAway, gridsize
       real :: inc

       farAway = 500.0 * gridsize

       if(inc .le. 0. .or. inc .ge. 360.) then
          unitvec = OCTALVECTOR(sin(inc*degtorad),1.d-10,cos(inc*degtorad))
       else
          unitvec = randomunitVector()
       endif

       posvec = faraway * unitvec
       centrevec = OCTALVECTOR(0d7,0d7,0d7)
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
     type(OCTALVECTOR) :: rvec
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
  type(OCTALVECTOR) :: currentposition(3), posvec, viewvec, unitvec
  integer :: subcell, i, itrans
  integer :: ilamb, nlamb
  real(double) :: xmidplane
  real :: lamb
  real(double) :: tau, dummy, kappaAbs, kappaSca, i0
  character(len=50) :: message
  
  call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule,0.d0,.true.)
  
  write(message,*) "Angular dependence"
  call writeinfo(message, FORINFO)
  
  do i = 1, 90
     call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
     
     write(message,*) i, tau
     call writeinfo(message, FORINFO)
  enddo
  
  write(message,*) "Tau from above"
  call writeinfo(message, FORINFO)
  
  do i = 1, 100
     
     unitvec = OCTALVECTOR(1.d-20,1.d-20,1.d0)
     posvec = OCTALVECTOR(real(i) * grid%octreeroot%subcellsize * 0.01,0d7,2.d0 * grid%octreeroot%subcellsize)
     viewvec = (-1.d0) * unitvec
     
     call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
     
     write(message,*) i, tau
     call writeinfo(message, FORINFO)
  enddo
  
  write(message,*) "Tau from the side"
  call writeinfo(message, FORINFO)
  
  do i = 1, 100
     unitvec = OCTALVECTOR(1.d0,0.d0,0.d0)
     posvec = OCTALVECTOR(2.d0 * grid%octreeroot%subcellsize,0.d0, real(i) * grid%octreeroot%subcellsize * 0.01)
     viewvec = (-1.d0) * unitvec
     call intensityAlongRay(posvec, viewvec, grid, thisMolecule, iTrans, 0.d0, i0, tau)
     
     write(message,*) i, tau
     call writeinfo(message, FORINFO)
  enddo

call intensityAlongRay(octalvector(0.d0,0.d0,0.d0),octalvector(1d-20,1d-20,1.d0), grid, thisMolecule, 1, 0.d0, dummy, tau, .true.)

xmidplane = rinner + (router-rinner) * 0.5

if(writeoutput) write(*,*) "Midplane tau!"
open(50, file = 'dustcheck.dat', status="unknown", form = "formatted") 
currentposition(1) = octalvector(xmidplane,0.d0,0.d0)

call findSubcellTD(currentPosition(1), grid%octreeroot, thisOctal, subcell)

do i = 0, nlamb
   
   lamb = (10.0**(real(i) * log10(lamend/lamstart) / 500.0)) * lamstart
   
   call continuumIntensityAlongRay(octalvector(1.d-10,1d-10,1d-10),octalvector(1.0,1d-20,1.d-20), grid, lamb, 0.d0, dummy, &
        tau, .true.)
   call locate(grid%lamArray, size(grid%lamArray), lamb, ilamb)
   call returnKappa(grid, thisOctal, subcell, ilambda = ilamb, lambda = lamb, kappaAbs = kappaAbs, kappaSca = kappaSca)
   write(50, '(f8.4,tr3,f10.4,tr3,f10.4,tr3f10.4)') lamb *1e-4, tau, kappaAbs*1e-10 / thisOctal%rho(subcell),&
        (kappaAbs+kappaSca)*1e-10 / thisOctal%rho(subcell)
enddo

close(50)

if(debug) then
   
   currentposition(1) =  OCTALVECTOR(xmidplane,0.d0,0.d0)
   currentposition(2) =  OCTALVECTOR(xmidplane,0.d0,xmidplane)
   currentposition(3) = OCTALVECTOR(xmidplane/sqrt(2.),xmidplane/sqrt(2.),2. * xmidplane)
   
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

call continuumIntensityAlongRay(octalvector(-1.d10,1.d-10,1.d-10),octalvector(1.d0,-1d-10,-1d-10), grid, 1e4,&
     0.d0, dummy, tau, .true.)
write(message,*) "TAU @ 1micron", tau
call writeinfo(message, FORINFO)

end subroutine testOpticalDepth

subroutine plotdiscValues(grid, thisMolecule)
  

  type(GRIDTYPE) :: grid
  type(MOLECULETYPE) :: thisMolecule
  real(double) :: mean(6)
  
  call calculateOctalParams(grid, grid%OctreeRoot, thisMolecule,0.d0,.true.)
  
  call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
       "disctemperature.ps/vcps", .true.,.false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize), fixValMin=1.d0, fixValMax=500.d0, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "rho", "x-z", real(grid%octreeRoot%centre%y), &
       "discrho.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , fixValMin=1.d-20, fixValMax=1.d-10, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "Vy", "x-z", real(grid%octreeRoot%centre%y), &
       "discVelocity.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize), show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
       "disctemperatureZoom.ps/vcps", .true.,.false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize), fixValMin=1.d0, fixValMax=500.d0,&
       show_value_3rd_dim=.false., boxfac = 0.005) 
  call plot_AMR_values(grid, "rho", "x-z", real(grid%octreeRoot%centre%y), &
       "discrhoZoom.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , fixValMin=1.d-20, fixValMax=1.d-10,&
       show_value_3rd_dim=.false., boxfac = 0.005) 
  call plot_AMR_values(grid, "Vy", "x-z", real(grid%octreeRoot%centre%y), &
       "discVelocityZoom.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.,boxfac = 0.005) 
  call plot_AMR_values(grid, "molAbundance", "x-z", real(grid%octreeRoot%centre%y), &
       "discmolAbundanceWrong.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.,&
       fixValMin=1d-13, fixValMax=1d-4)          
  call plot_AMR_values(grid, "molAbundance", "x-z", real(grid%octreeRoot%centre%y), &
       "discmolAbundanceZoomWrong.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.,&
       fixValMin=1d-13, fixValMax=1d-4, boxfac = 0.005)          
  
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 2)
  
  call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
       "Cdisctemperature.ps/vcps", .true.,.false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) ,  fixValMin=1.d0, fixValMax=500.d0, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "rho", "x-z", real(grid%octreeRoot%centre%y), &
       "Cdiscrho.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , fixValMin=1.d-20, fixValMax=1.d-10, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "Vy", "x-z", real(grid%octreeRoot%centre%y), &
       "CdiscVelocity.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
       "CdisctemperatureZoom.ps/vcps", .true.,.false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) ,  fixValMin=1.d0, fixValMax=500.d0,&
       show_value_3rd_dim=.false., boxfac = 0.005) 
  call plot_AMR_values(grid, "rho", "x-z", real(grid%octreeRoot%centre%y), &
       "CdiscrhoZoom.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , fixValMin=1.d-20, fixValMax=1.d-10,&
       show_value_3rd_dim=.false., boxfac = 0.005) 
  call plot_AMR_values(grid, "Vy", "x-z", real(grid%octreeRoot%centre%y), &
       "CdiscVelocityZoom.ps/vcps", .false., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , fixValMin=0.d0, fixValMax=50.d0,&
       show_value_3rd_dim=.false.,boxfac = 0.005) 
  call plot_AMR_values(grid, "molAbundance", "x-z", real(grid%octreeRoot%centre%y), &
       "discmolAbundance.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.,&
       fixValMin=1d-13, fixValMax=1d-4)          
  call plot_AMR_values(grid, "molAbundance", "x-z", real(grid%octreeRoot%centre%y), &
       "discmolAbundanceZoom.ps/vcps", .true., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , show_value_3rd_dim=.false.,&
       fixValMin=1d-13, fixValMax=1d-4, boxfac = 0.005)          
           
  call readAMRgrid("molecular_tmp.grid",.false.,grid)
  
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 0)
  
  call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
       "tempdiff02.ps/vcps", .false.,.false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) ,  fixValMin=0.d0, fixValMax=2.d0, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
       "tempdiff0911.ps/vcps", .false.,.false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) ,  fixValMin=0.9d0, fixValMax=1.1d0, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "temperature", "x-z", real(grid%octreeRoot%centre%y), &
       "tempdiff099101.ps/vcps", .false.,.false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) ,  fixValMin=0.99d0, fixValMax=1.01d0, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "rho", "x-z", real(grid%octreeRoot%centre%y), &
       "rhodiff02.ps/vcps", .false., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , fixValMin=0.d0, fixValMax=2.d0, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "rho", "x-z", real(grid%octreeRoot%centre%y), &
       "rhodiff0911.ps/vcps", .false., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , fixValMin=0.9d0, fixValMax=1.1d0, show_value_3rd_dim=.false.) 
  call plot_AMR_values(grid, "rho", "x-z", real(grid%octreeRoot%centre%y), &
       "rhodiff099101.ps/vcps", .false., .false., &
       width_3rd_dim=real(grid%octreeRoot%subcellsize) , fixValMin=0.99d0, fixValMax=1.01d0, show_value_3rd_dim=.false.) 
  call readAMRgrid("molecular_tmp.grid",.false.,grid)
  
  ! Set everything back to the way it was?
  call findtempdiff(grid, grid%OctreeRoot, thisMolecule, mean, 1)

end subroutine plotdiscValues

     function velocity(position, grid, startOctal, subcell) RESULT(out)

       implicit none
       type(VECTOR) :: out
       type(octalVector), intent(in) :: position
       type(gridtype), intent(in) :: grid
       type(octal), pointer, optional :: startOctal
       integer, optional :: subcell

    select case (grid%geometry)

       case("molebench")
          out = molebenchVelocity(position, grid)
          
       case("IRAS04158","shakara")
          out = keplerianVelocity(position, grid)

       case default
          if(present(startOctal) .and. present(subcell)) then
             out = amrGridVelocity(grid%octreeRoot, position, startOctal = startOctal, actualSubcell = subcell)
          else
             out = amrGridVelocity(grid%octreeRoot, position)
          endif
    end select

  end function velocity

end module molecular_mod
