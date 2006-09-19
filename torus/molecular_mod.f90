module molecular_mod

! written by tjh


  use kind_mod
  use constants_mod
  use utils_mod
  use messages_mod
  use grid_mod
  use math_mod

  implicit none

  type MOLECULETYPE
     character(len=10) :: molecule
     real :: molecularweight
     real(double) :: abundance
     integer :: nLevels
     real(double), pointer :: energy(:)
     real(double), pointer :: g(:)
     real(double), pointer :: j(:)
     integer :: nTrans
     real(double), pointer :: einsteinA(:,:)
     real(double), pointer :: einsteinB(:,:)
     real(double), pointer :: transfreq(:,:)
     real(double), pointer :: Eu(:,:)
     integer :: nCollPart
     character(len=20), pointer :: collBetween(:)
     integer, pointer :: nCollTrans(:)
     integer, pointer :: nCollTemps(:)
     real(double), pointer :: collTemps(:,:)
     real(double), pointer :: collRates(:,:,:,:)
  end type MOLECULETYPE

contains

  subroutine readMolecule(thisMolecule, molFilename)
    type(MOLECULETYPE) :: thisMolecule
    character(len=*) :: molFilename
    character(len=80) :: junk, message
    character(len=200):: dataDirectory, filename
    integer :: i, j, iLow, iUp, iPart
    real(double) :: a, freq, eu, c(20)

    thisMolecule%abundance = 1.d-9 ! fixed at benchmark value here

    call unixGetenv("TORUS_DATA", dataDirectory, i)
    filename = trim(dataDirectory)//"/"//molfilename

    open(30, file=filename, status="old", form="formatted")

    read(30,*) junk
    read(30,'(a)') thisMolecule%molecule

    call writeInfo("Reading data for molecule: "//trim(thisMolecule%molecule),IMPORTANT)
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

       thisMolecule%energy(i) = thisMolecule%energy(i) / 8065.541  ! per cm to ev
       
    enddo

    read(30,*) junk
    read(30,*) thisMolecule%nTrans

    allocate(thisMolecule%einsteinA(1:thisMolecule%nLevels,1:thisMolecule%nLevels))
    allocate(thisMolecule%einsteinB(1:thisMolecule%nLevels,1:thisMolecule%nLevels))
    allocate(thisMolecule%transfreq(1:thisMolecule%nLevels,1:thisMolecule%nLevels))
    allocate(thisMolecule%eu(1:thisMolecule%nLevels,1:thisMolecule%nLevels))
    read(30,*) junk
    do i = 1, thisMolecule%nTrans
       read(30,*) j, iUp, iLow, a, freq, eu

       thisMolecule%einsteinA(iUp, iLow) = a
       thisMolecule%transfreq(iUp, iLow) = freq*1.d9
       thisMolecule%eu(iUp, iLow) = eu

       thisMolecule%einsteinB(iLow, iUp) = (thisMolecule%g(iUp)/thisMolecule%g(iLow)) * a * &
            (cspeed**2)/(2.d0*hcgs*(freq*1.d9)**3)
       thisMolecule%einsteinB(iUp, iLow) = thisMolecule%einsteinB(iLow, iUp) &
            * thisMolecule%g(iLow)/thisMolecule%g(iUp)

    enddo

    read(30,*) junk
    read(30,*) thisMolecule%nCollPart
    
    allocate(thisMolecule%nCollTrans(1:thisMolecule%nCollPart))
    allocate(thisMolecule%nCollTemps(1:thisMolecule%nCollPart))
    allocate(thisMolecule%collTemps(1:thisMolecule%nCollPart, 1:20))
    allocate(thisMolecule%collBetween(1:thisMolecule%nCollPart))

    
    do iPart = 1, thisMolecule%nCollPart

       read(30,*) junk
       read(30,*) thisMolecule%collBetween(iPart)
       
       read(30,*) junk
       read(30,*) thisMolecule%nCollTrans(iPart)

       read(30,*) junk
       read(30,*) thisMolecule%nCollTemps(iPart)


       read(30,*) junk
       read(30,*) thisMolecule%collTemps(iPart,1:thisMolecule%nCollTemps(ipart))

       allocate(thisMolecule%collRates(1:thisMolecule%nCollPart, &
            1:thisMolecule%nLevels, 1:thisMolecule%nLevels, 1:thisMolecule%nCollTemps(ipart)))

       read(30,*) junk
       do j = 1, thisMolecule%nCollTrans(iPart)
          read(30,*) i, iUp, iLow, c(1:thisMolecule%nCollTemps(iPart))
          thisMolecule%collRates(iPart, iUp, iLow, 1:thisMolecule%nCollTemps(iPart)) = c(1:thisMolecule%nCollTemps(iPart))
       enddo
    enddo
    close(30)
    call writeInfo("Done.", IMPORTANT)
  end subroutine readMolecule


  subroutine solveLevels(nPops, jnu,  temperature, thisMolecule, nh2)
    real(double) :: nPops(:)
    real(double) :: temperature
    real(double) :: jnu(:)
    real(double) :: nh2
    type(MOLECULETYPE) :: thisMolecule
    real(double), allocatable :: matrixA(:,:), matrixB(:)
    real(double) :: arateji, boltzFac
    integer :: nLevels
    integer :: iLower, iUpper, iLevel, i, j
    real(double) :: collEx, colldeEx

    nLevels = thisMolecule%nLevels

    allocate(matrixA(1:nLevels,1:nLevels))
    allocate(matrixB(1:nLevels))

    matrixA = 1.d-30
    matrixB = 0.d0

    matrixA(1,:) = 1.d0
    matrixB(1) = 1.d0

    
    do i = 2, nLevels

       ! ways to drop to lower levels (i=upper, j =lower)

       if (i > 1) then
          do j = 1, i-1
             colldeEx = collRate(thisMolecule, temperature, i , j) * nh2 * (nh2 * thisMolecule%abundance)

             matrixA(i,i) = matrixA(i,i) - thisMolecule%einsteinA(i,j) ! spont emiss
             matrixA(i,i) = matrixA(i,i) - thisMolecule%einsteinB(i,j)*jnu(i) ! stim emiss
             matrixA(i,i) = matrixA(i,i) - colldeEx
          enddo
       endif

       ! ways to lose upwards from this level (i=lower,j=upper)
       if (i < nLevels) then
          do j = i+1, nLevels
             boltzFac =  exp(-abs(thisMolecule%energy(i)-thisMolecule%energy(j)) / (kev*temperature))
             colldeEx = collRate(thisMolecule, temperature, j , i) * nh2 * (nh2 * thisMolecule%abundance)
             collEx = colldeEx * boltzFac

             matrixA(i,i) = matrixA(i,i) - thisMolecule%einsteinB(i,j)*jnu(j) ! stim abs
             matrixA(i,i) = matrixA(i,i) - collEx ! collisional ex
          enddo
       endif

       ! ways to gain from lower levels (i=upper,j=lower)

       if (i > 1) then
          do j = 1, i-1
             boltzFac =  exp(-abs(thisMolecule%energy(i)-thisMolecule%energy(j)) / (kev*temperature))
             colldeEx = collRate(thisMolecule, temperature, i , j) * nh2 * (nh2 * thisMolecule%abundance)
             collEx = colldeEx * boltzFac

             matrixA(i,j) = matrixA(i,j) + thisMolecule%einsteinB(i,j)*jnu(i) ! stim abs
             matrixA(i,j) = matrixA(i,j) + collEx
          enddo
       endif

       ! ways to gain from higher levels (i=lower, j=upper)

       if (i < nLevels) then
          do j = i+1, nLevels
             boltzFac =  exp(-abs(thisMolecule%energy(i)-thisMolecule%energy(j)) / (kev*temperature))
             colldeEx = collRate(thisMolecule, temperature, j , i) * nh2 * (nh2 * thisMolecule%abundance)
             collEx = colldeEx * boltzFac

             matrixA(i,j) = matrixA(i,j) + thisMolecule%einsteinA(j,i) ! spont emiss
             matrixA(i,j) = matrixA(i,j) + thisMolecule%einsteinB(i,j)*jnu(j) ! stim emiss
             matrixA(i,j) = matrixA(i,j) + colldeEx
          enddo
       endif
    enddo


  call luSlv(matrixA, matrixB, nLevels)

  matrixB(1:nLevels) = matrixB(1:nLevels) / SUM(matrixB(1:nLevels))
  nPops(1:nLevels) = matrixB(1:nLevels)
  
  deallocate(matrixA, matrixB)

  end subroutine solveLevels


  real(double) function collRate(thisMolecule, temperature, i , j)
    type(MOLECULETYPE) :: thisMolecule
    real(double) :: temperature, r
    integer :: i, j, k, iPart

    collRate = 0.d0

    do iPart = 1, thisMolecule%nCollPart
       
       call locate(thisMolecule%collTemps(iPart,1:thisMolecule%nCollTemps(iPart)), &
            thisMolecule%nCollTemps(iPart), temperature, k)

       r = (temperature - thisMolecule%collTemps(iPart,k)) / &
            (thisMolecule%collTemps(iPart,k+1) - thisMolecule%collTemps(iPart, k))

       collRate = collRate + thisMolecule%collRates(iPart, i , j, k) + &
            r * ( thisMolecule%collRates(iPart, i , j, k+1) -  thisMolecule%collRates(iPart, i , j, k))
    enddo
  end function collRate


  recursive subroutine  solveAllPops(grid, thisOctal, thisMolecule)
    type(GRIDTYPE) :: grid
    type(MOLECULETYPE) :: thisMolecule
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iUpper, iLower
    real, parameter :: Tcbr = 2.728
    real(double) :: etaLine
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call solveAllpops(grid, child, thisMolecule)
                exit
             end if
          end do
       else
          call solveLevels(thisOctal%molecularLevel(subcell,1:thisMolecule%nLevels), &
               thisOctal%jnu(subcell,1:thisMolecule%nLevels),  &
               dble(thisOctal%temperature(subcell)), thisMolecule, thisOctal%nh2(subcell))
       endif
    enddo
  end subroutine solveAllPops

  recursive subroutine  allocateMolecularLevels(grid, thisOctal, thisMolecule)
    type(GRIDTYPE) :: grid
    type(MOLECULETYPE) :: thisMolecule
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iUpper, iLower
    real(double) :: etaLine
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call allocateMolecularLevels(grid, child, thisMolecule)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%molecularLevel)) then
             allocate(thisOctal%molecularLevel(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
          endif
          thisOctal%molecularLevel = 1.d-30

          if (.not.associated(thisOctal%jnu)) then
             allocate(thisOctal%jnu(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
          endif
          thisOctal%jnu = 1.d-30

          if (.not.associated(thisOctal%jnugrid)) then
             allocate(thisOctal%jnugrid(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
          endif
          thisOctal%jnu = 1.d-30

       endif
    enddo
  end subroutine allocateMolecularLevels




  function phiProf(dv, b) result (phi)
    real(double) :: dv, b
    real(double) :: fac, phi
    phi = 1.d0 / (b * sqrt(pi))
    fac = dv**2 / b**2
    phi = phi * exp(-fac)
  end function phiProf




  subroutine getRay(grid, position, direction, ds, phi, i0, iTransUpper, thisMolecule)
    type(MOLECULETYPE) :: thisMolecule
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, startOctal
    integer :: subcell
    real(double) :: ds, phi, i0, r
    integer :: iTransUpper
    type(OCTALVECTOR) :: position, direction, currentPosition, thisPosition, thisVel
    type(OCTALVECTOR) :: rayVel
    real(double) :: alphanu, snu, jnu
    integer :: iLower , iUpper
    real(double) :: dv, deltaV
    integer :: i
    real(double) :: distArray(200), tval
    integer :: nTau
    real(double) :: nLower, nUpper
    real(double) :: dTau, etaline, didtau, tau
    real(double), parameter :: Tcbr = 2.782d0
    real(double) :: intensityIntegral

    thisOctal => grid%octreeRoot
    call findSubcellLocal(position, thisOctal, subcell)


    position = randomPositionInCell(thisOctal, subcell)
    direction = randomUnitVector()
    rayVel = amrGridVelocity(grid%octreeRoot, position, startOctal = thisOctal, actualSubcell = subcell)
    call random_number(r)
    deltaV = 5.d0 * thisOctal%microturb(subcell) * (2.d0*r-1.d0)

    iUpper = iTransUpper
    iLower = iTransUpper - 1

    call distanceToCellBoundary(grid, position, direction, ds, sOctal=thisOctal)

    currentPosition = position + ds * direction

    ds = ds * 1.d10

    phi = phiProf(deltaV, thisOctal%microturb(subcell))



    i0 = 0.d0
    intensityIntegral = 0.0
    tau = 0.d0



    do while(inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition, thisOctal, subcell)
       call distanceToCellBoundary(grid, currentPosition, direction, tVal, sOctal=thisOctal)

       ntau = 10

       do i = 2, nTau
          
          distArray(i) = tval * dble(i-1)/dble(nTau-1)

          startOctal => thisOctal
          thisPosition = currentPosition + distArray(i)*direction
          thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell) 
          thisVel= thisVel - rayVel


          dv = (thisVel .dot. direction) + deltaV

          alphanu = (hCgs*thisMolecule%transFreq(iUpper, iLower)/fourPi) * &
               phiProf(dv, thisOctal%microturb(subcell))/thisMolecule%transFreq(iUpper,iLower)

          nLower = thisOctal%molecularLevel(subcell,iLower) * thisMolecule%abundance * thisOctal%nh2(subcell)
          nUpper = thisOctal%molecularLevel(subcell,iUpper) * thisMolecule%abundance * thisOctal%nh2(subcell)

          alphanu = alphanu * (nLower * thisMolecule%einsteinB(iLower, iUpper) - &
               nUpper * thisMolecule%einsteinB(iUpper, iLower))

          dTau = alphaNu *  (distArray(i)-distArray(i-1)) * 1.d10

          etaLine = hCgs * thisMolecule%einsteinA(iUpper, iLower) * thisMolecule%transFreq(iUpper, iLower)
          etaLine = etaLine * thisOctal%nh2(subcell) * thisMolecule%abundance * thisOctal%molecularLevel(subcell, iUpper)
          jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell))/thisMolecule%transFreq(iUpper,iLower)

          if (alphanu > 0.d0) then
             snu = jnu/alphanu
          else
             snu = tiny(snu)
          endif

          dIdtau = -intensityIntegral + snu

          intensityIntegral = intensityIntegral + didtau * dtau
          tau = tau + dtau
          i0 = i0 + intensityIntegral * exp(-tau)
       enddo
       i0 = i0 + bnu(thisMolecule%transFreq(iUpper, iLower), Tcbr) * exp(-tau)
       currentPosition = currentPosition + distArray(ntau) * direction
    enddo
  end subroutine getRay



  subroutine dumpResults(grid, thisMolecule)
    type(GRIDTYPE) :: grid
    type(MOLECULETYPE) :: thisMolecule
    integer :: nr
    real(double) :: r, ang
    integer :: i
    real(double) :: pops(10)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: j
    real(double) :: x, z
    type(OCTALVECTOR) :: posvec
    real(double) :: router, rinner
    
    rinner = 1.e6
    router  = 4.6e7
    
    nr = 20
    open(31,file="results.dat",status="unknown",form="formatted")
    do i = 1, nr
       r = log10(rinner) + dble(i-1)/dble(nr-1)*(log10(router) - log10(rinner))
       r = 10.d0**r
       pops = 0.d0
       do j = 1 , 20
          call random_number(ang)
          ang = ang * pi
          z = r*cos(ang)
          x = r*sin(ang)
          posVec = OCTALVECTOR(x, 0.d0, z)
          thisOctal => grid%octreeroot
          call findSubcellLocal(posVec, thisOctal,subcell)
          pops = pops + thisOctal%molecularLevel(subcell,1:10)
       enddo
       pops = pops / 20.d0
       write(31,*) real(r*1.d10), real(pops(1:6))
    enddo
    close(31)
  end subroutine dumpResults

          
  subroutine molecularLoop(grid, thisMolecule)
    type(GRIDTYPE) :: grid
    type(MOLECULETYPE) :: thisMolecule
    type(OCTALVECTOR) :: position, direction
    integer :: iTransUpper
    integer :: nOctal, iOctal, subcell
    real(double), allocatable :: ds(:,:), phi(:,:), i0(:,:)
    integer :: nRay
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    type(OCTAL), pointer :: thisOctal
    integer, parameter :: maxIter = 10
    logical :: converged
    integer :: iRay, iUpper, iter
    real(double), allocatable :: oldpops(:)
    real(double) :: fac

    call allocateMolecularLevels(grid, grid%octreeRoot, thisMolecule)


    nRay = 10
    allocate(ds(1:thisMolecule%nLevels, 1:nRay))
    allocate(phi(1:thisMolecule%nLevels, 1:nRay))
    allocate(i0(1:thisMolecule%nLevels, 1:nRay))

    allocate(oldPops(1:thisMolecule%nLevels))
    

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)

    do while (.true.)
       do iOctal = 1, grid%nOctals
          thisOctal => octalArray(iOctal)%content
          do subcell = 1, thisOctal%maxChildren
             
             if (.not.thisOctal%hasChild(subcell)) then
                
                do iUpper = 2, thisMolecule%nLevels
                   do iRay = 1, nRay
                      call getRay(grid, position, direction, &
                           ds(iUpper,iRay), phi(iUpper,iRay), i0(iUpper,iRay), iUpper, thisMolecule)
                   enddo
                enddo
                iter = 0
                converged = .false.
                do while (.not.converged)
                   iter = iter + 1
                   oldpops = thisOctal%molecularLevel(subcell,1:thisMolecule%nLevels)
                   do iUpper = 2, thisMolecule%nLevels
                      call calculateJbar(thisOctal, subcell, thisMolecule, nRay, ds(iUpper,1:nRay), &
                           phi(iUpper,1:nRay), i0(iUpper,1:nRay), iUpper, thisOctal%jnu(subcell,iUpper))
                   enddo
                   call solveLevels(thisOctal%molecularLevel(subcell,1:thisMolecule%nLevels), &
                        thisOctal%jnu(subcell,1:thisMolecule%nLevels),  &
                        dble(thisOctal%temperature(subcell)), thisMolecule, thisOctal%nh2(subcell))
                   fac = maxval((thisOctal%molecularLevel(subcell,1:thisMolecule%nLevels) - oldpops)/oldpops)
                   if (fac < 1.d-6) converged = .true.
                   if (iter == maxIter) converged = .true.
                enddo
             endif
          enddo
       end do
       write(*,*) "Dumping results"
       call dumpresults(grid, thisMolecule)
    enddo

  end subroutine molecularLoop


  subroutine calculateJbar(thisOctal, subcell, thisMolecule, nRay, ds, phi, i0, iTransUpper, jbar)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(MOLECULETYPE) :: thisMolecule
    integer :: nRay
    real(double) :: ds(:), phi(:), i0(:)
    integer :: iTransUpper
    real(double) :: jbar
    integer :: iRay
    real(double) :: nLower, nUpper
    real(double) :: jBarInternal, jBarExternal
    real(double) :: alphanu, jnu, etaline
    integer :: iUpper, iLower
    real(double) :: tau, snu

    jBarExternal = 0.d0
    jBarInternal = 0.d0

    iUpper = iTransUpper
    iLower = iUpper - 1

    do iRay = 1, nRay
       alphanu = (hCgs*thisMolecule%transFreq(iUpper, iLower)/fourPi)
       nLower = thisOctal%molecularLevel(subcell,iLower) * thisMolecule%abundance * thisOctal%nh2(subcell)
       nUpper = thisOctal%molecularLevel(subcell,iUpper) * thisMolecule%abundance * thisOctal%nh2(subcell)
       
       alphanu = alphanu * (nLower * thisMolecule%einsteinB(iLower, iUpper) - &
            nUpper * thisMolecule%einsteinB(iUpper, iLower)) * phi(iray)/thisMolecule%transFreq(iUpper,iLower)
       tau = alphaNu * ds(iray)
       

       etaLine = hCgs * thisMolecule%einsteinA(iUpper, iLower) * thisMolecule%transFreq(iUpper, iLower)
       etaLine = etaLine * thisOctal%nh2(subcell) * thisMolecule%abundance * thisOctal%molecularLevel(subcell, iUpper)
       jnu = (etaLine/fourPi) * phi(iRay)/thisMolecule%transFreq(iUpper,iLower)
       
       if (alphanu > 0.d0) then
          snu = jnu/alphanu
       else
          snu = tiny(snu)
       endif

       jBarExternal = jBarExternal + i0(iray) * exp(-tau)
       jBarInternal = jBarInternal + snu * (1.d0 - exp(-tau))

    enddo
    
    jbar = (jBarExternal + jBarExternal)/dble(nRay)

  end subroutine calculateJbar



end module molecular_mod
