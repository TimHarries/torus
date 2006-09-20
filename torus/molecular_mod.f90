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
     real(double), pointer :: einsteinA(:)
     real(double), pointer :: einsteinBlu(:)
     real(double), pointer :: einsteinBul(:)
     real(double), pointer :: transfreq(:)
     real(double), pointer :: itransUpper(:)
     real(double), pointer :: itransLower(:)
     real(double), pointer :: Eu(:)
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

       thisMolecule%einsteinBlu(i) = (thisMolecule%g(iUp)/thisMolecule%g(iLow)) * a * &
            (cspeed**2)/(2.d0*hcgs*(freq*1.d9)**3)
       thisMolecule%einsteinBul(i) = thisMolecule%einsteinBlu(i) &
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
    integer :: itrans, l, k
    real(double) :: collEx, colldeEx

    nLevels = thisMolecule%nLevels

    allocate(matrixA(1:nLevels+1,1:nLevels+1))
    allocate(matrixB(1:nLevels+1))

    matrixA = 1.d-30
    matrixB = 0.d0

    matrixA(nLevels+1,1:nLevels+1) = 1.d0
    matrixA(1:nLevels+1,nLevels+1) = 1.d-30

    matrixB(nLevels+1) = 1.d0

    do iTrans = 1, thisMolecule%nTrans

       k = thisMolecule%iTransUpper(iTrans)
       l = thisMolecule%iTransLower(iTrans)
    

       matrixA(k,k) = matrixA(k,k) + thisMolecule%einsteinBul(iTrans) * jnu(iTrans) + thisMolecule%einsteinA(iTrans)
       matrixA(l,l) = matrixA(l,l) + thisMolecule%einsteinBlu(iTrans) * jnu(iTrans)
       matrixA(k,l) = matrixA(k,l) - thisMolecule%einsteinBlu(iTrans) * jnu(iTrans)
       matrixA(l,k) = matrixA(l,k) - thisMolecule%einsteinBul(iTrans) * jnu(iTrans) - thisMolecule%einsteinA(iTrans)

    enddo

!    do k = 2, nLevels
!       do l = 1, k - 1 
!          boltzFac =  exp(-abs(thisMolecule%energy(k)-thisMolecule%energy(l)) / (kev*temperature))
!          colldeEx = collRate(thisMolecule, temperature, k , l) * nh2 * (nh2 * thisMolecule%abundance)
!          collEx = colldeEx * boltzFac
!         matrixA(k,k) = matrixA(k,k) + colldeex
!          matrixA(l,l) = matrixA(l,l) + collex
!          matrixA(k,l) = matrixA(k,l) - collex
!          matrixA(l,k) = matrixA(l,k) - collex -  colldeex
!       enddo
!    enddo


  call luSlv(matrixA, matrixB, nLevels+1)



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
               thisOctal%jnu(subcell,1:thisMolecule%nTrans),  &
               dble(thisOctal%temperature(subcell)), thisMolecule, thisOctal%nh2(subcell))
       endif
    enddo
  end subroutine solveAllPops

  recursive subroutine  swapPops(thisOctal)
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
                call swapPops(child)
                exit
             end if
          end do
       else

          thisOctal%molecularLevel(subcell,:) = thisOctal%newmolecularLevel(subcell,:)
       endif
    enddo
  end subroutine swapPops

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

          if (.not.associated(thisOctal%newmolecularLevel)) then
             allocate(thisOctal%newmolecularLevel(1:thisOctal%maxChildren, 1:thisMolecule%nLevels))
          endif
          thisOctal%newmolecularLevel = 1.d-30

          if (.not.associated(thisOctal%jnu)) then
             allocate(thisOctal%jnu(1:thisOctal%maxChildren, 1:thisMolecule%nTrans))
          endif
          do i = 1, thisMolecule%nTrans
             thisOctal%jnu(subcell,i) = bnu(thisMolecule%transFreq(i), dble(thisOctal%temperature(subcell)))
          enddo

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




  subroutine getRay(grid, fromOctal, fromSubcell, position, direction, ds, phi, i0, iTrans, thisMolecule, tau)
    type(MOLECULETYPE) :: thisMolecule
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, startOctal, fromOctal
    integer :: fromSubcell
    integer :: subcell
    real(double) :: ds, phi, i0, r
    integer :: iTrans
    type(OCTALVECTOR) :: position, direction, currentPosition, thisPosition, thisVel
    type(OCTALVECTOR) :: rayVel, startVel, endVel, endPosition
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
    real(double) :: dvAcrossCell

    position = randomPositionInCell(fromOctal, fromsubcell)


    thisOctal => grid%octreeRoot
    call findSubcellLocal(position, thisOctal, subcell)


    direction = randomUnitVector()
    rayVel = amrGridVelocity(grid%octreeRoot, position, startOctal = thisOctal, actualSubcell = subcell)
    call random_number(r)
    deltaV = 0.d0 * thisOctal%microturb(subcell) * (2.d0*r-1.d0)!!!!!!!!!!!!!!!!!!

    iUpper = thisMolecule%iTransUpper(iTrans)
    iLower = thisMolecule%iTransLower(iTrans)

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

       startVel = amrGridVelocity(grid%octreeRoot, currentPosition, startOctal = thisOctal, actualSubcell = subcell) 
       endPosition = currentPosition + tval * direction
       endVel = amrGridVelocity(grid%octreeRoot, endPosition)

       dvAcrossCell = ((startVel - rayVel).dot.direction) - ((endVel - rayVel).dot.direction)
       dvAcrossCell = abs(dvAcrossCell / thisOctal%microturb(subcell))

       if (dvAcrossCell < 0.1) then
          ntau = 2
       else
          ntau = 5
       endif

       distArray(1) = 0.d0
       do i = 2, nTau
          
          distArray(i) = tval * dble(i-1)/dble(nTau-1)

          startOctal => thisOctal
          thisPosition = currentPosition + distArray(i)*direction
          thisVel = amrGridVelocity(grid%octreeRoot, thisPosition, startOctal = startOctal, actualSubcell = subcell) 
          thisVel= thisVel - rayVel


          dv = (thisVel .dot. direction) + deltaV

          alphanu = (hCgs*thisMolecule%transFreq(iTrans)/fourPi) * &
               phiProf(dv, thisOctal%microturb(subcell))/thisMolecule%transFreq(iTrans)

          nLower = thisOctal%molecularLevel(subcell,iLower) * thisMolecule%abundance * thisOctal%nh2(subcell)
          nUpper = thisOctal%molecularLevel(subcell,iUpper) * thisMolecule%abundance * thisOctal%nh2(subcell)

          alphanu = alphanu * (nLower * thisMolecule%einsteinBlu(iTrans) - &
               nUpper * thisMolecule%einsteinBul(iTrans))

          dTau = alphaNu *  (distArray(i)-distArray(i-1)) * 1.d10

          etaLine = hCgs * thisMolecule%einsteinA(iTrans) * thisMolecule%transFreq(iTrans)
          etaLine = etaLine * thisOctal%nh2(subcell) * thisMolecule%abundance * thisOctal%molecularLevel(subcell, iUpper)
          jnu = (etaLine/fourPi) * phiProf(dv, thisOctal%microturb(subcell))/thisMolecule%transFreq(iTrans)


          if (alphanu /= 0.d0) then
             snu = jnu/alphanu
          else
             snu = tiny(snu)
          endif

          i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
          tau = tau + dtau
       enddo
       currentPosition = currentPosition + distArray(ntau) * direction
    enddo
    i0 = i0 + bnu(thisMolecule%transFreq(iTrans), Tcbr) * exp(-tau)
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
    real(double) :: router, rinner, jnu
    
    rinner = 1.e6
    router  = 4.6e7
    
    nr = 20
    open(31,file="results.dat",status="unknown",form="formatted")
    do i = 1, nr
       r = log10(rinner) + dble(i-1)/dble(nr-1)*(log10(router) - log10(rinner))
       r = 10.d0**r
       pops = 0.d0
       jnu = 0.d0
       do j = 1 , 20
          call random_number(ang)
          ang = ang * pi
          z = r*cos(ang)
          x = r*sin(ang)
          posVec = OCTALVECTOR(x, 0.d0, z)
          thisOctal => grid%octreeroot
          call findSubcellLocal(posVec, thisOctal,subcell)
          pops = pops + thisOctal%molecularLevel(subcell,1:10)
          jnu = jnu + thisOctal%jnu(subcell,1)
       enddo
       pops = pops / 20.d0
       jnu = jnu / 20.d0
       write(31,*) real(r*1.d10), real(pops(1:6)),jnu
    enddo
    close(31)
  end subroutine dumpResults

          
  subroutine molecularLoop(grid, thisMolecule)
    type(GRIDTYPE) :: grid
    type(MOLECULETYPE) :: thisMolecule
    type(OCTALVECTOR) :: position, direction
    integer :: nOctal, iOctal, subcell
    real(double), allocatable :: ds(:,:), phi(:,:), i0(:,:), tau(:)
    integer :: nRay
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    type(OCTAL), pointer :: thisOctal
    integer, parameter :: maxIter = 10
    logical :: converged
    integer :: iRay, iTrans, iter
    real(double), allocatable :: oldpops(:)
    real(double) :: fac

    call allocateMolecularLevels(grid, grid%octreeRoot, thisMolecule)


    call solveAllPops(grid, grid%octreeRoot, thisMolecule)
    write(*,*) "Dumping results",log10(bnu(thisMolecule%transFreq(1),2.782d0))
    call dumpresults(grid, thisMolecule)

    nRay = 10
    allocate(ds(1:thisMolecule%nLevels, 1:nRay))
    allocate(phi(1:thisMolecule%nLevels, 1:nRay))
    allocate(i0(1:thisMolecule%nLevels, 1:nRay))
    allocate(tau(1:nRay))

    allocate(oldPops(1:thisMolecule%nLevels))
    

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)

    do while (.true.)
       do iOctal = 1, grid%nOctals
          write(*,*) iOctal
          thisOctal => octalArray(iOctal)%content
          do subcell = 1, thisOctal%maxChildren
             
             if (.not.thisOctal%hasChild(subcell)) then
                
                do iTrans = 1, thisMolecule%nTrans
                   do iRay = 1, nRay
                      call getRay(grid, thisOCtal, subcell, position, direction, &
                           ds(iTrans,iRay), phi(iTrans,iRay), i0(iTrans,iRay), iTrans, thisMolecule, tau(iray))
                   enddo
!                   if (iTrans == 1) write(*,*) "tau",sum(tau(1:nRay))/dble(nRay)
                enddo
                iter = 0
                converged = .false.
                do while (.not.converged)
                   iter = iter + 1
                   oldpops = thisOctal%molecularLevel(subcell,1:thisMolecule%nLevels)
                   do iTrans = 1, thisMolecule%nTrans
                      call calculateJbar(thisOctal, subcell, thisMolecule, nRay, ds(iTrans,1:nRay), &
                           phi(iTrans,1:nRay), i0(iTrans,1:nRay), iTrans, thisOctal%jnu(subcell,iTrans))
                   enddo
                   call solveLevels(thisOctal%newmolecularLevel(subcell,1:thisMolecule%nLevels), &
                        thisOctal%jnu(subcell,1:thisMolecule%nTrans),  &
                        dble(thisOctal%temperature(subcell)), thisMolecule, thisOctal%nh2(subcell))
                   fac = maxval((thisOctal%newmolecularLevel(subcell,1:thisMolecule%nLevels) - oldpops)/oldpops)
                   if (fac < 1.d-6) converged = .true.
                   if (iter == maxIter) converged = .true.
                enddo
             endif
          enddo
       end do
       call swapPops(grid%octreeRoot)
       write(*,*) "Dumping results"
       call dumpresults(grid, thisMolecule)
    enddo

  end subroutine molecularLoop


  subroutine calculateJbar(thisOctal, subcell, thisMolecule, nRay, ds, phi, i0, iTrans, jbar)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(MOLECULETYPE) :: thisMolecule
    integer :: nRay
    real(double) :: ds(:), phi(:), i0(:)
    integer :: iTrans
    real(double) :: jbar
    integer :: iRay
    real(double) :: nLower, nUpper
    real(double) :: jBarInternal, jBarExternal
    real(double) :: alphanu, jnu, etaline
    integer :: iUpper, iLower
    real(double) :: tau, snu

    jBarExternal = 0.d0
    jBarInternal = 0.d0

    iUpper = thisMolecule%iTransUpper(iTrans)
    iLower = thisMolecule%iTransLower(iTrans)

    do iRay = 1, nRay
       alphanu = (hCgs*thisMolecule%transFreq(iTrans)/fourPi)
       nLower = thisOctal%molecularLevel(subcell,iLower) * thisMolecule%abundance * thisOctal%nh2(subcell)
       nUpper = thisOctal%molecularLevel(subcell,iUpper) * thisMolecule%abundance * thisOctal%nh2(subcell)
       
       alphanu = alphanu * (nLower * thisMolecule%einsteinBlu(iTrans) - &
            nUpper * thisMolecule%einsteinBul(iTrans)) * phi(iray)/thisMolecule%transFreq(iTrans)
       tau = alphaNu * ds(iray)
       

       etaLine = hCgs * thisMolecule%einsteinA(iTrans) * thisMolecule%transFreq(iTrans)
       etaLine = etaLine * thisOctal%nh2(subcell) * thisMolecule%abundance * thisOctal%molecularLevel(subcell, iUpper)
       jnu = (etaLine/fourPi) * phi(iRay)/thisMolecule%transFreq(iTrans)
       
       if (alphanu /= 0.d0) then
          snu = jnu/alphanu
       else
          snu = tiny(snu)
       endif

       jBarExternal = jBarExternal + i0(iray) * exp(-tau)
       jBarInternal = jBarInternal + snu * (1.d0 - exp(-tau))

    enddo
    
    jbar = (jBarExternal + jBarInternal)/dble(nRay)

  end subroutine calculateJbar



end module molecular_mod
