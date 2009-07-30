module timedep_mod

  ! time dependent RT module

  use kind_mod
  use constants_mod
  use vector_mod
  use utils_mod
  use grid_mod
  use amr_mod
  use source_mod
  use vtk_mod
  use lucy_mod
  use math_mod

  implicit none

  type PHOTONSTACK
     integer :: nStack
     type(VECTOR) :: position(1000000)
     type(VECTOR) :: direction(1000000)
     real(double) :: eps(1000000)
     real(double) :: freq(1000000)
     logical :: fromsource(1000000)
     logical :: fromGas(1000000)
  end type PHOTONSTACK


  private

  public :: runtimeDependentRT

contains

  subroutine runTimeDependentRT(grid, source, nSource, nLambda, xArray)
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nLambda
    real :: xArray(:)
    integer :: iter
    real(double) :: deltaT
    integer :: nMonte
    type(PHOTONSTACK) :: stack
    character(len=80) :: vtkFilename

    nMonte = 1000
    iter = 0
    deltaT = 1.d-25
    do while (.true.)
       iter = iter + 1
       call  timeDependentRTStep(grid, nsource, source, stack, deltaT, nMonte, xArray, nLambda)
       write(vtkFilename, '(a,i4.4,a)') "output",iter,".vtk"
       call writeVTKFile(grid, vtkfilename)
    enddo
  end subroutine runTimeDependentRT

  subroutine timeDependentRTStep(grid, nsource, source, oldstack, deltaT, nMonte, lamArray, nLambda)
    type(GRIDTYPE) :: grid
    integer :: nSource
    integer :: nLambda
    integer :: nFreq
    real :: lamArray(:)
    real(double) :: deltaT
    type(SOURCETYPE) :: source(:)
    type(PHOTONSTACK), intent(inout) :: oldStack
    type(PHOTONSTACK) :: currentStack
    logical :: photonFromSource, photonFromGas
    integer :: nMonte, iMonte
    integer :: i
    real(double) :: freqArray(2000), dnu(2000), kAbsArray(2000), kAbsArray2(2000)
    real(double) :: prob(2000), luminosity, sourceLuminosity
    integer :: nFromMatter, nPhotons
    real(double) :: fracSource, chanceSource, chanceGas, weightSource, weightGas
    type(VECTOR) :: rVec, uHat, rHat
    real(double) :: eps, r, freq,  photonTime, tDble
    real(double) :: wavelength
    real(double) :: distToBoundary, tauToBoundary, distanceToEvent, timeToBoundary
    real(double) :: photonSpeed
    real(double) :: tau
    integer :: iLambda
    real :: treal
    type(SOURCETYPE) :: thisSource
    type(OCTAL), pointer :: thisOctal
    integer :: subcell, iSource, iLam
    logical :: absorbed, scattered, outOfTime, finished, ok
    real(double) :: kappaAbs, kappaSca, albedo
    real(double) :: newDeltaT, fac
    real(double) :: totalLineEmission, totalContEmission
    real :: lambda0
    lambda0 = 0.

    photonSpeed = cSpeed
    nFreq = nLambda
    do i = 1, nFreq
       freqArray(nFreq-i+1) = cSpeed / (lamArray(i)*1.e-8)
    enddo

    do i = 2, nFreq-1
       dnu(i) = 0.5*((freqArray(i+1)+freqArray(i))-(freqArray(i)+freqArray(i-1)))
    enddo

    dnu(1) = freqArray(2)-freqArray(1)
    dnu(nFreq) = freqArray(nFreq)-freqArray(nFreq-1)



    call zeroDistanceGrids(grid%octreeRoot)
    kabsArray = 0.d0
    call calculateEtaCont(grid, grid%octreeRoot, nFreq, freqArray, dnu, lamarray, nLambda, kabsArray)

    call computeProbDist(grid, totalLineEmission, totalContEmission, lambda0, .false.)
    luminosity = 0.d0
    call calculateGasEmissivity(Grid%octreeRoot, luminosity)
    write(*,*) "emissivity ",luminosity
    write(*,*) "time step ",deltaT

    currentStack%nStack = 0


    sourceLuminosity = SUM(source(1:nSource)%luminosity)

    nFromMatter = nMonte
    nPhotons = oldStack%nStack + nFromMatter
    fracSource = 0.8d0

    chanceSource = sourceLuminosity / (sourceLuminosity + luminosity)
    chanceGas = 1.d0-chanceSource
    weightSource = chanceSource / fracSource
    weightGas = (1.d0-chanceSource) / (1.d0-fracSource)
    do iMonte = 1, nPhotons
       write(*,*) "imonte ", imonte, nPhotons 
       if (oldStack%nStack > 0)  then
          call getPhotonFromStack(oldStack, rVec, uHat, eps, freq, &
               photonFromSource, photonFromGas)
       else
          call random_number(r)

          if (r < fracSource) then

             call randomSource(source, nSource, iSource)
             thisSource = source(iSource)
             call getPhotonPositionDirection(thisSource, rVec, uHat, rHat,grid)

             call getWavelength(thisSource%spectrum, wavelength)
             freq = cSpeed/(wavelength / 1.e8)
             photonTime = 0.d0
             eps = weightSource * sourceluminosity * deltaT / dble(nFromMatter)
             photonFromSource = .true.
             photonFromGas = .false.
          else

             call random_number(r)
             call locateContProbAMR(r,thisOctal,subcell)
             rVec = randomPositionInCell(thisOctal, subcell)
             call amrGridValues(grid%octreeRoot, rVec, startOctal=thisOctal, &
                  actualSubcell=subcell, temperature=treal,grid=grid, kappaAbsArray=kAbsArray)
             tdble = dble(tReal)
             do i = 1, nLambda
                iLam = nfreq - i + 1
                kAbsArray2(i) = kabsArray(ilam)
             enddo
             prob(1) = 0.d0
             do i = 2, nFreq
                prob(i) = prob(i-1) + kAbsArray2(i)*bnu(freqArray(i), tdble)*dnu(i)
             enddo
             prob = prob / prob(nFreq)
             call random_number(r)
             call locate(prob, nFreq, r, i)
             freq = freqArray(i)
             wavelength = (cSpeed/freq)*1.d8
             eps = weightGas * luminosity * deltaT / dble(nFromMatter)
             photonFromGas = .true.
             photonFromSource = .false.
             photonTime = 0.d0
             uHat = randomUnitVector()
          endif

       endif
       
       ilambda = findIlambda(real(wavelength), lamArray, nLambda, ok)

       photonTime = 0.d0

       absorbed = .false.
       scattered = .false.
       finished = .false.
       outOfTime = .false.
       thisOctal => grid%octreeRoot
       subcell = 1

       do while (.not.finished)
          call findSubcellLocal(rVec, thisOctal, subcell)
          call distanceToCellBoundary(grid, rVec, uHat, distToBoundary, sOctal=thisOctal)
          call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, kappaAbs = kappaAbs, kappaSca = kappaSca)

          timeToBoundary = distToBoundary / photonSpeed

          tauToBoundary = distToBoundary * (kappaAbs + kappaSca)

          albedo  = kappaSca / (kappaAbs + kappaSca)

          call random_number(r)
          tau = -log(1.d0 - r)

          if (tau > tauToBoundary) then
             distanceToEvent = distToBoundary
             absorbed = .false.
          else
             distanceToEvent = distToBoundary * tau/tauToBoundary
             call random_number(r)
             if (r > albedo) then
                absorbed = .true.
                finished = .true.
             else
                absorbed = .false.
                finished = .false.
                scattered = .true.
             endif
          endif

          if (absorbed.or.scattered) then
             if ((photonTime + distanceToEvent/photonSpeed) > deltaT) then
                absorbed = .false.
                scattered = .false.
                distanceToEvent = (deltaT - photonTime) * photonSpeed
                outOfTime = .true.
                finished = .true.
             endif
          else
             if ((photonTime + distanceToEvent/photonSpeed) > deltaT) then
                distanceToEvent = (deltaT - photonTime) * photonSpeed
                outOfTime = .true.
                finished = .true.
             endif
          endif

          photonTime = photonTime + distanceToEvent / photonSpeed


          call updateDistanceGrids(thisOctal, subcell, kappaAbs, eps, distanceToEvent, photonFromGas, &
               photonFromSource)

          rVec = rVec + (distanceToEvent + 1.d-3*grid%halfSmallestSubcell) * uHat

          if (.not.inOctal(grid%octreeRoot, rVec)) then
             finished = .true.
          endif
          if (scattered) then
             uhat = randomUnitVector()  ! ISOTROPIC SCATTERING ONLY SO FAR
             scattered = .false.
             photonFromGas = .true.
          endif
       end do
       if (outOfTime) then
          call addPhotonToStack(currentStack, rVec, uHat, photonTime, &
               eps, photonFromSource, photonFromGas)
       endif
    end do

    call calculateADot(grid%octreeRoot, deltaT)
    call calculateEnergyDensity(grid%octreeRoot, deltaT, photonSpeed)


    newDeltaT = 1.d30
    call calculateNewDeltaT(grid%octreeRoot, newDeltaT)
    fac = 1.d-3
    newDeltaT = newDeltaT  * fac


    call updateUDens(grid%octreeRoot, deltaT)
    call calculateTemperatureFromUdens(grid%octreeRoot)

    oldStack = currentStack
    currentStack%nStack = 0
  end subroutine timeDependentRTStep

  subroutine updateDistanceGrids(thisOctal, subcell, kappaAbs, eps, distanceToEvent, photonFromGas, &
       photonFromSource)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: kappaAbs
    real(double) :: eps
    real(double) :: distanceToEvent
    logical :: photonFromGas
    logical :: photonFromSource


    thisOctal%distanceGridAdot(subcell) = thisOctal%distanceGridAdot(subcell) + &
         distanceToEvent * kappaAbs* eps
    if (photonFromSource) then
       thisOctal%distanceGridPhotonFromSource(subcell) = thisOctal%distanceGridPhotonFromSource(subcell) + &
            distanceToEvent * eps
    endif
    if (photonFromGas) then
       thisOctal%distanceGridPhotonFromGas(subcell) = thisOctal%distanceGridPhotonFromGas(subcell) + &
            distanceToEvent * eps
    endif

  end subroutine updateDistanceGrids

  recursive subroutine updateUdens(thisOctal, deltaT)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: deltaT
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call updateUdens(child, deltaT)
                exit
             end if
          end do
       else
          thisOctal%udens(subcell) = thisOctal%udens(subcell) + deltaT * &
               (thisOctal%adot(subcell) - thisOctal%etaCont(subcell))
       endif
    enddo
  end subroutine updateUdens

  recursive subroutine calculateGasEmissivity(thisOctal, ems)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: ems, v
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateGasEmissivity(child, ems)
                exit
             end if
          end do
       else
          v = cellVolume(thisOctal, subcell) * 1.d30
          ems = ems + thisOctal%etaCont(subcell) * v
       endif
    enddo
  end subroutine calculateGasEmissivity

  recursive subroutine calculateTemperatureFromUdens(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateTemperatureFromUdens(child)
                exit
             end if
          end do
       else
          thisOctal%temperature(subcell) = &
               temperatureFunc(thisOctal%uDens(subcell), thisOctal%rho(subcell), 5.d0/3.d0)
       endif
    enddo
  end subroutine calculateTemperatureFromUdens

  recursive subroutine calculateNewDeltaT(thisOctal, deltaT)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: deltaT
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateNewDeltaT(child, deltaT)
                exit
             end if
          end do
       else
          deltaT = min(deltaT, thisOctal%uDens(subcell)/thisOctal%etaCont(subcell))
       endif
    enddo
  end subroutine calculateNewDeltaT

  recursive subroutine calculateAdot(thisOctal, deltaT)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: deltaT, v
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateAdot(child, deltaT)
                exit
             end if
          end do
       else
          v = cellVolume(thisOctal, subcell) * 1.d30
          thisOctal%aDot(subcell) = (1.d0 / v) * thisOctal%distancegridAdot(subcell) / deltaT
       endif
    enddo
  end subroutine calculateAdot

  recursive subroutine calculateEnergyDensity(thisOctal, deltaT, photonSpeed)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: deltaT, v, photonSpeed
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateEnergyDensity(child, deltaT, photonSpeed)
                exit
             end if
          end do
       else
          v = cellVolume(thisOctal, subcell) * 1.d30
          thisOctal%photonenergyDensityFromGas(subcell) = (1.d0/photonSpeed) * (1.d0 / v) * &
               thisOctal%distancegridphotonFromGas(subcell) / deltaT
          thisOctal%photonenergyDensityFromSource(subcell) =  (1.d0/photonSpeed) * (1.d0 / v) * &
               thisOctal%distancegridphotonFromSource(subcell) / deltaT
       endif
    enddo
  end subroutine calculateEnergyDensity

  recursive subroutine zeroDistanceGrids(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroDistanceGrids(child)
                exit
             end if
          end do
       else
          thisOctal%distancegridAdot(subcell) = 0.d0
          thisOctal%distancegridphotonFromGas(subcell) = 0.d0
          thisOctal%distancegridphotonFromSource(subcell) = 0.d0
       endif
    enddo
  end subroutine zeroDistanceGrids



  subroutine timeDependentRTtest()
    real(double) :: currentTime, endTime
    integer, parameter :: nx = 3
    integer, parameter :: nProb = nx + 1
    real(double) :: etaCont(nx)
    real(double) :: aDot(nx)
    real(double) :: rho(nx)
    real(double) :: kappa(nx)
    real(double) :: temperature(nx)
    real(double) :: distanceGridAdot(nx)
    real(double) :: distanceGridPhotonFromGas(nx)
    real(double) :: distanceGridPhotonFromSource(nx)
    logical :: photonFromSource, photonFromGas
    real(double) :: prob(nProb)
!    real(double) :: photonProb(nProb)
!    real(double) :: photonEnergy
    real(double) :: xCen(nx)
    real(double) :: xArray(nprob)
    real(double) :: uDens(nx)
    real(double) :: photonArray(nx)
    real(double) :: photonEnergyDensityFromSource(nx)
    real(double) :: photonEnergyDensityFromGas(nx)
    real(double) :: uDensAnalytical(nx)
    real(double) :: photonDensAnalytical(nx)
    real(double) :: tauBox
    real(double) :: luminosity
    real(double) :: photonTime
    logical :: absorbed, finished
    real(double) :: distToBoundary, timeToBoundary, tauToBoundary, r
    real(double) :: distanceToEvent, deltaT
    real(double) :: tau, dx, k
    real(double) :: xSize
    real(double) :: xPhoton
    integer, parameter :: maxStack = 10000000

    integer :: currentNStack, oldnStack
    real(double) :: currenttimeStack(maxStack)
    type(VECTOR) :: currentpositionStack(maxStack)
    type(VECTOR) :: currentdirectionStack(maxStack)
    real(double) :: currentepsStack(maxStack)
    logical :: currentFromsourceStack(maxStack)
    logical :: currentFromGasStack(maxStack)

    real(double) :: oldtimeStack(maxStack)
    type(VECTOR) :: oldpositionStack(maxStack)
    type(VECTOR) :: olddirectionStack(maxStack)
    real(double) :: oldepsStack(maxStack)
    logical :: oldfromsourceStack(maxStack), oldFromGasStack(maxStack)

    real(double) :: epsOverDeltaT, newDeltaT, tot, albedo, diffusionTime
!    real(double) ::t
    logical :: outOFtime, scattered
    type(VECTOR) :: rVec, uHat
    integer :: i, iMonte, nMonte, iPos,iter, nPhotons
    character(len=30) :: outfile
    real(double) :: photonSpeed,t 
    real(double) :: sourcePeriod, sourceLuminosity
    real(double) :: fracSource, weightSource, chanceSource, weightGas, chanceGas
    real(double) :: fac, meanFreePath, probNewPhoton
    integer :: nFromMatter

    fracSource = 0.0d0
    udensAnalytical  = 0.d0
    udens = 1.d-30
    photonSpeed = cSpeed
    nMonte = 1000000
    tauBox = 1.d3
    deltaT = 1.d0
    call random_seed

    albedo = 0.d0
    rho = 1.d-10
!    kappa = 4.e-1
    xSize = 1.d0 !autocm
    sourcePeriod = 100.d0
    deltaT = 1.d-30
    diffusionTime = 1.d30
    sourceLuminosity = 0.d0
    rho = 1.d0
    kappa = 1.d0

    uDens = 1.d-10
    udens = 1.d17
    dx = xSize/dble(nx)
    xCen(1) = dx/2.d0
    do i = 2, nx
       xCen(i) = xCen(i-1) + dx 
    enddo
    do i = 1, nprob
       xArray(i) = xSize * dble(i-1)/dble(nProb-1)
    enddo
!    do i = 1, nx
!       if (xCen(i) > 0.5d0*xSize) then
!          kappa(i) = 1.d7 / (rho(i)*xSize)
!       endif
!    enddo

!    do i = 1, nx
!!       uDens(i) = max(1.d1,cos(2.d0*xCen(i))*1.d7)
!       if (xCen(i) < 0.5d0*xSize) then
!          uDens(i) = 1.e14
!       else
!          uDens(i) = 1.d-10
!       endif
!       temperature(i) = temperatureFunc(uDens(i),  rho(i), 5.d0/3.d0)
!    enddo
!    do i = 1, nx
!       t = temperatureFunc(uDens(i), rho(i), 5.d0/3.d0)
!       write(*,*) "temperature ",t
!!       photonEnergyDensity(i) = max(1.d1,cos(2.d0*xCen(i))*1.d7)!(fourPi * (stefanBoltz/pi) * t**4 / cSpeed) 
!       photonEnergyDensity(i) = (fourPi * (stefanBoltz/pi) * t**4 / photonSpeed) 
!       photonEnergyDensity(i) = (8.d0*pi**5*(kerg*t)**4)/(15.d0 * (hCgs * photonSpeed)**3)
!    enddo
!    deltaT = 1.d-30

!    call calculateProb(photonProb, xArray, xCen, nProb, photonEnergyDensity)
!
!    oldNstack = nMonte
!    photonEnergy = SUM(photonEnergyDensity*dx)/dble(oldNStack)
!    photonArray = 0.d0
!    do i = 1,oldnStack
!       call getX(xArray, photonprob, nProb, xPhoton)
!       oldpositionStack(i)%x = xPhoton
!       oldPositionStack(i)%y = 0.d0
!       oldPositionStack(i)%z = 0.d0
!       olddirectionStack(i) = randomUnitVector2()
!       oldtimestack(i) = 0.d0
!       oldEpsStack(i) = photonEnergy 
!       call findArrayIndex(xCen, nx, xPhoton, iPos)
!       photonArray(iPos) = photonArray(iPos) + oldEpsStack(i)
!    enddo
!    do i = 1 ,nx
!       write(*,*) "udens " , udens(i), "rad ",photonEnergyDensity(i), &
!            " photon ",Photonarray(i)/dx
!    enddo
!
!    kappa = tauBox / (rho*xSize)
!    write(*,*) "total tau ",SUM(rho*kappa*dx)
!
    uDensAnalytical = uDens
    photonDensAnalytical = 0.d0
!    udensAnalytical = photonEnergyDensity
!    kappa = 4.d-1

!    photonSpeed = kappa(1)*rho(1)
    k = photonSpeed/(kappa(1)*rho(1))
!    write(*,*) "Diffusion coefficent ",k
!    diffusionTime=  0.3d0 * (dx**2/(2.d0*k))
    deltaT = 1.d-26!diffusionTime
    currentTime = 0.d0
    endTime = 1.d0
    iter = 0
    currentnStack = 0
    oldNStack = 0
    open(69,file="time.dat", status="unknown",form="formatted")
    do while (currentTime < 1.d30)
       iter = iter + 1
       meanFreePath = 1.d0/(kappa(1)*rho(1))

       distanceGridAdot = 0.d0
       distanceGridPhotonFromSource = 0.d0
       distanceGridPhotonFromGas = 0.d0
       do i = 1, nx
          temperature(i) = temperatureFunc(uDens(i),  rho(i), 5.d0/3.d0)
!          write(*,*) temperature(i)
       enddo
       etaCont = fourPi * (stefanBoltz/pi) * temperature**4 * kappa * rho 
       call calculateProb(prob, xArray, xCen, nProb, etaCont)
       luminosity = SUM(etaCont*dx)
       write(*,*) "emissivity ",luminosity
       write(*,*) "time step ",deltaT
!       write(*,*) "diffusion timestep ", 0.3d0*(dx**2/(2.d0*k))

       write(*,*) "current Time", currentTime
       if (oldNStack == 0) then
          write(*,*) "energy monte ",SUM(udens*dx)
       else
          write(*,*) "energy monte ",SUM(udens*dx)+SUM(oldepsStack(1:oldnStack))
          write(*,*) "number of photons ", oldNStack
          write(*,*) "percent in photons ", 100.d0 * &
               SUM(oldepsStack(1:oldnStack))/(SUM(udens*dx)+SUM(oldepsStack(1:oldnstack)))
       endif
       write(*,*) "energy analy ",SUM(udensanalytical*dx)
       write(*,*) "light distance ", deltaT*photonSpeed
       write(*,*) "mean free path ", 1.d0/(kappa(1)*rho(1))
       write(*,*) "diffusion speed (speed of light units)", ((1.d0/(kappa(1)*rho(1)))/deltaT)/photonSpeed


       nFromMatter = max(1, nMonte-oldNStack)
       nPhotons = oldNStack + nFromMatter
       tot = 0.d0
!       sourceLuminosity = 1.d10 * 0.5d0*(1.d0+cos(twoPi*currentTime/sourcePeriod))
       write(*,*) " source luminosity ", sourceluminosity

       chanceSource = sourceLuminosity / (sourceLuminosity + luminosity)
       chanceGas = 1.d0-chanceSource
       weightSource = chanceSource / fracSource
       weightGas = (1.d0-chanceSource) / (1.d0-fracSource)
       do iMonte = 1, nPhotons
!          write(*,*) "imonte ", imonte, nPhotons 
          if (oldnStack > 0)  then
             call getPhotonFromStacktest(oldnStack, rVec, uHat, photonTime, epsOverDeltaT, &
                  photonFromSource, photonFromGas, &
                  oldpositionStack, olddirectionStack, oldtimeStack, oldepsStack, oldFromSourceStack, &
                  oldFromGasStack)
          else
             call random_number(r)

             if (r < fracSource) then
                rVec = VECTOR(0.d0, 0.d0, 0.d0)
                uHat = VECTOR(1.d0, 0.d0, 0.d0)
                photonTime = 0.d0
                epsOverDeltaT = weightSource * sourceluminosity * deltaT / dble(nFromMatter)
                photonFromSource = .true.
                photonFromGas = .false.
             else
                call getX(xArray, prob, nProb, xPhoton)
                call findArrayIndex(xCen, nx, xPhoton, iPos)
                epsOverDeltaT =  weightGas * luminosity * deltaT / dble(nFromMatter)
                photonFromGas = .true.
                photonFromSource = .false.
                rVec%x = xPhoton
                rVec%y = 0.d0
                rVec%z = 0.d0
                photonTime = 0.d0
                uHat = randomUnitVector2()
                tot = tot + epsOverDeltaT
             endif
          endif
          photonTime = 0.d0

          absorbed = .false.
          scattered = .false.
          finished = .false.
          outOfTime = .false.
          do while (.not.finished)


             call findArrayIndex(xCen, nx, rVec%x, iPos)
             if (uHat%x > 0.d0) then
                distToBoundary = ((xCen(iPos)+dx/2.d0) - rVec%x)/uHat%x
             else
                distToBoundary = abs( (((xCen(iPos)-dx/2.d0) - rVec%x )/uHat%x ))
             endif

             timeToBoundary = distToBoundary / photonSpeed

             tauToBoundary = distToBoundary * kappa(iPos) * rho(iPos)

!             write(*,*) "tautoboundary ", tautoboundary, iPos, rVec%x, xCen(iPos)-dx/2.d0, xCen(iPos)+dx/2.d0
             call random_number(r)
             tau = -log(1.d0 - r)


             if (tau > tauToBoundary) then
                distanceToEvent = distToBoundary
                absorbed = .false.
             else
                distanceToEvent = distToBoundary * tau/tauToBoundary
                call random_number(r)
                if (r > albedo) then
                   absorbed = .true.
                   finished = .true.
                else
                   absorbed = .false.
                   finished = .false.
                   scattered = .true.
                endif
             endif

             if (absorbed.or.scattered) then
                if ((photonTime + distanceToEvent/photonSpeed) > deltaT) then
                   absorbed = .false.
                   scattered = .false.
                   distanceToEvent = (deltaT - photonTime) * photonSpeed
                   outOfTime = .true.
                   finished = .true.
                endif
             else
                if ((photonTime + distanceToEvent/photonSpeed) > deltaT) then
                   distanceToEvent = (deltaT - photonTime) * photonSpeed
                   outOfTime = .true.
                   finished = .true.
                endif
             endif

             photonTime = photonTime + distanceToEvent / photonSpeed

             distanceGridAdot(iPos) = distanceGridAdot(iPos) + &
                  distanceToEvent * (1.d0-albedo) * kappa(iPos) * rho(iPos) * epsOverDeltaT
             if (photonFromSource) then
                distanceGridPhotonFromSource(iPos) = distanceGridPhotonFromSource(iPos) + &
                     distanceToEvent * epsOverDeltaT
             endif
             if (photonFromGas) then
                distanceGridPhotonFromGas(iPos) = distanceGridPhotonFromGas(iPos) + &
                     distanceToEvent * epsOverDeltaT
             endif
             rVec = rVec + (distanceToEvent + 1.d-7*dx) * uHat
             if (rVec%x < (xCen(1)-dx/2.d0)) then
                rVec%x =  xCen(1)-dx/2.d0+1.d-7*dx
                uHat%x = -1.d0 * uHat%x
!                finished = .true.
             endif
             if (rVec%x > xCen(nx)+dx/2.d0) then
                rVec%x = xCen(nx)+dx/2.d0-1.d-7*dx
                uHat%x = -1.d0 * uHat%x
!                finished = .true.
             endif
             if (scattered) then
                uhat = randomUnitVector2()
                scattered = .false.
                photonFromGas = .true.
             endif
          end do
          if (outOfTime) then
             call addPhotonToStacktest(currentnStack, rVec, uHat, photonTime, &
                  epsOverDeltaT, photonFromSource, photonFromGas, &
                  currentpositionStack, currentdirectionStack, &
                  currenttimeStack, currentepsStack, currentFromSourceStack, currentFromGasStack)
             photonArray(iPos) = photonArray(iPos) + epsOverDeltaT
          endif
       end do
       write(*,*) "Sum of energy in photon array ",SUM(photonArray(1:nx))
       call calculateADottest(distanceGridAdot, aDot, xCen, dx, nx, deltaT)
       call calculatePhotonEnergyDensity(distanceGridPhotonFromGas, &
            photonEnergyDensityFromGas, xCen, dx, nx, deltaT, photonSpeed)
       call calculatePhotonEnergyDensity(distanceGridPhotonFromSource, &
            photonEnergyDensityFromSource, xCen, dx, nx, deltaT, photonSpeed)

!       call solveNewUdens(xCen, uDensAnalytical, nx, k, deltaT)


       newDeltaT = 1.d30
       fac = 1.d-3
       do i = 1, nx
             newDeltaT = min(newDeltaT, fac*uDens(i)/etacont(i))
!          newDeltaT = min(newDeltaT, fac * uDens(i) / abs(adot(i) - etaCont(i)))
       enddo

       deltaT = min(newDeltaT,diffusionTime)



       do i = 1, nx
          uDens(i) = uDens(i) + (aDot(i) - etaCont(i)) * deltaT
          if (i == 1) write(*,*) "udens ",udens(i), " adot ",adot(1), " eta ",etacont(i)
          uDens(i) = max (0.d0, uDens(i))
          temperature(i) = temperatureFunc(uDens(i), rho(i),  5.d0/3.d0)
       enddo

       do i = 1 , 1
          t = temperatureFunc(uDensAnalytical(i), rho(i), 5.d0/3.d0)
          udensAnalytical(i) = udensAnalytical(i) + &
               deltaT * (photonSpeed*kappa(i)*rho(i)*photonDensAnalytical(i) &
               - fourPi*kappa(i)*rho(i)*(stefanBoltz/pi)*t**4)
          write(*,*) "before anal ",photonDensAnalytical(i)
          photondensAnalytical(i) = photondensAnalytical(i) - &
               deltaT * (photonSpeed*kappa(i)*rho(i)*photonDensAnalytical(i) &
               - fourPi*kappa(i)*rho(i)*(stefanBoltz/pi)*t**4)
          write(*,*) "after anal ",photonDensAnalytical(i)
       enddo

       probNewPhoton = 0.d0!

       currentTime = currentTime + deltaT

          write(outfile,'(a,i6.6,a)') "udens",iter,".dat"
          write(*,*) iter
          open(21, file=outfile, status="unknown", form="formatted")
          write(21, '(a,1pe12.5)') "# current time ",currentTime
          do i = 1, nx
             write(21, *) xCen(i), uDens(i),  photonEnergyDensityFromSource(i), &
                  photonEnergyDensityFromGas(i)
          enddo
          close(21)
          write(69,*) currentTime, SUM(udens*dx),SUM(photonEnergyDensityFromGas*dx),currentNStack, &
               udensAnalytical(1)*xSize,photonDensAnalytical(1)*xSize
          write(*,*) "TOTAL ENERGY: ",SUM(udens*dx)+SUM(photonEnergyDensityFromGas*dx)

       oldnStack = currentnStack
       oldPositionStack(1:currentnStack) = currentPositionStack(1:currentnStack)
       oldDirectionStack(1:currentnStack) = currentDirectionStack(1:currentnStack)
       oldTimeStack(1:currentnStack) = currentTimeStack(1:currentnStack)
       oldEpsStack(1:currentnStack) = currentEpsStack(1:currentnStack)
       oldFromSourceStack(1:currentnStack) = currentFromSourcestack(1:currentNstack)
       oldFromGasStack(1:currentnStack) = currentFromGasstack(1:currentNstack)
       currentNStack = 0
    enddo
    close(69)
  end subroutine timeDependentRTtest

  subroutine addPhotonToStacktest(nStack, rVec, uVec, photonTime, &
       epsOverDeltaT, photonFromSource, photonFromGas, &
       positionStack, directionStack, timeStack, epsStack, &
       fromSourceStack, FromGasStack)
    integer :: nStack
    real(double) :: photonTime, epsOverDeltaT
    logical :: photonFromSource, photonFromGas
    type(VECTOR) :: rVec, uVec
    type(VECTOR), intent(inout) :: positionStack(:), directionStack(:)
    real(double), intent(inout) :: timeStack(:), epsStack(:)
    logical, intent(inout) :: fromSourceStack(:), fromGasStack(:)
    
    nStack = nStack + 1
!    write(*,*) "adding photon to stack ", nStack
    if (nStack > SIZE(positionStack)) then
       write(*,*) "nStack too large for array ",nStack
       stop
    endif
    positionStack(nStack) = rVec
    directionStack(nStack) = uVec
    timeStack(nStack) = photonTime
    epsStack(nStack) = epsOverDeltaT
    fromSourceStack(nStack) = photonFromSource
    fromGasStack(nStack) = photonFromGas
  end subroutine addPhotonToStacktest

  subroutine addPhotonToStack(stack, rVec, uVec, &
       eps, freq, photonFromSource, photonFromGas)
    type(PHOTONSTACK) :: stack
    real(double) :: eps, freq
    logical :: photonFromSource, photonFromGas
    type(VECTOR) :: rVec, uVec
    
    stack%nStack = stack%nStack + 1
    if (stack%nStack > SIZE(stack%position)) then
       write(*,*) "nStack too large for array ",stack%nStack
       stop
    endif
    stack%position(Stack%nStack) = rVec
    stack%direction(Stack%nStack) = uVec
    stack%eps(Stack%nStack) = eps
    stack%freq(Stack%nStack) = freq
    stack%fromSource(Stack%nStack) = photonFromSource
    stack%fromGas(Stack%nStack) = photonFromGas
  end subroutine addPhotonToStack

  subroutine getPhotonFromStacktest(nStack, rVec, uHat, photonTime, epsOverDeltaT, photonFromSource, photonFromGas, &
                  positionStack, directionStack, timeStack, epsStack, fromSourceStack, fromGasStack)
    integer :: nStack
    type(VECTOR) :: rVec
    type(VECTOR) :: uHat
    real(double) :: photonTime
    real(double) :: epsOverDeltaT
    logical :: photonFromSource, photonFromGas
    type(VECTOR) :: positionStack(:)
    type(VECTOR) :: directionStack(:)
    real(double) :: timeStack(:)
    real(double) :: epsStack(:)
    logical :: fromSourceStack(:), fromGasStack(:)

!    write(*,*) "getting photon from stack ",nStack
    rVec = positionStack(nStack)
    uHat = directionStack(nStack)
    photonTime = timeStack(nStack)
    epsOverDeltaT = epsStack(nStack)
    photonFromSource = fromSourceStack(nStack)
    photonFromGas = fromGasStack(nStack)
    nStack = nStack - 1
  end subroutine getPhotonFromStacktest

  subroutine getPhotonFromStack(stack, rVec, uHat, eps, freq, photonFromSource, photonFromGas)
    type(PHOTONSTACK) :: stack
    type(VECTOR) :: rVec
    type(VECTOR) :: uHat
    real(double) :: eps, freq
    logical :: photonFromSource, photonFromGas
    
    if (stack%nStack == 0) then
       write(*,*) "BUG: getPhotonFromStack called with empty stack"
       stop
    endif
    rVec = stack%position(stack%nStack)
    uHat = stack%direction(stack%nStack)
    eps = stack%eps(stack%nStack)
    freq = stack%freq(stack%nStack)
    photonFromSource = stack%fromSource(stack%nStack)
    photonFromGas = stack%fromGas(stack%nStack)
    stack%nStack = stack%nStack - 1
  end subroutine getPhotonFromStack


  subroutine calculateADottest(distanceGrid, aDot, xCen, dx, nx, deltaT)
    real(double) :: distanceGrid(:)
    real(double), intent(out) :: aDot(:)
    real(double) :: xCen(:), dx, deltaT
    integer :: nx, i


    do i = 1, nx
       aDot(i) = (1.d0 / dx) * distancegrid(i) / deltaT
    enddo
  end subroutine calculateADottest

  subroutine calculatePhotonEnergyDensity(distanceGrid, photonEnergyDensity, xCen, dx, nx, deltaT, photonSpeed)
    real(double) :: distanceGrid(:), photonSpeed
    real(double), intent(out) :: photonEnergyDensity(:)
    real(double) :: xCen(:), dx, deltaT
    integer :: nx, i


    do i = 1, nx
       PhotonEnergyDensity(i) = (1.d0 / dx) * (1.d0/photonSpeed) * distancegrid(i) / deltaT
    enddo
  end subroutine calculatePhotonEnergyDensity

  subroutine calculateProb(prob,  xArray, xCen, nProb, etaCont)
    real(double), intent(inout) :: xArray(:), prob(:)
    real(double) :: xCen(:), etaCont(:), totalEmission, dx
    integer :: i, nProb

    dx = xCen(2)-xCen(1)
    prob(1) = 0.d0
    totalEmission = 0.d0
    totalEmission = SUM(etaCont*dx)
    prob(1) = 0.d0
    do i = 2, nProb
       prob(i) = prob(i-1) + etaCont(i-1) * dx
    enddo
    prob(1:nProb) = prob(1:nProb) / prob(nProb)


!    open(30,file="prob.dat",status="unknown",form="formatted")
!    do i = 1, nProb
!       write(30,*) xArray(i), prob(i)
!    enddo
!    close(30)
!    stop
  end subroutine calculateProb
       
  subroutine getX(xArray, prob, nProb, x)
    real(double) :: xArray(:), prob(:)
    real(double), intent(out) :: x
    integer :: nProb
    real(double) :: r
    integer :: i

    call random_number(r)
    call locate(prob, nprob, r, i)
    x = xArray(i) + (xArray(i+1)-xArray(i))*(r - prob(i))/(prob(i+1)-prob(i))
  end subroutine getX

  subroutine findArrayIndex(xCen, nx, x, iPos)
    real(double) :: xCen(:), x, dx
    integer :: nx
    integer, intent(out) :: iPos

    dx = xCen(2) - xCen(1)
    if (x < xCen(1)-dx/2.d0) then
       write(*,*) "x value outside array ",x,xcen(1)-dx/2.d0
       stop
    endif
    if (x > xCen(nx)+dx/2.d0) then
       write(*,*) "x value outside array ",x
       stop
    endif
    if (x > xCen(nx)) then
       iPos = nx
       goto 666
    endif
    if (x < xCen(1)) then
       ipos = 1
       goto 666
    endif
    
     call locate(xCen, nx, x, iPos)
     if (x > xCen(iPos)+dx/2.d0) then
        iPos = iPos + 1
     endif
666 continue
  end subroutine findArrayIndex

  subroutine solveNewUdens(xArray, u, nx, k, deltaT)
    real(double) :: xArray(:), u(:), k, deltaT,dx
    real(double) :: uxx, uPrime(1000)
    integer :: nx, i
    dx = xArray(2)-xArray(1)
    
    do i = 2, nx-1
       uxx = (u(i+1) - 2.d0* u(i) + u(i-1)) / dx**2
       uprime(i) = u(i) + k * uxx * deltaT
    enddo
    uprime(1) = uprime(2)
    uprime(nx) = uprime(nx-1)

    u = uprime
  end subroutine solveNewUdens

  real(double) function temperatureFunc(u, rho,  gamma) result (t)
    real(double) :: u, gamma, rho
    real(double), parameter :: mu = 1.d0


    t = (gamma-1.d0)*mu*u / (rGas * rho)
!    t = (gamma-1.d0)*mu*u / (rGas)

  end function temperatureFunc

  type(VECTOR) function randomUnitVector2() result (v)
    real(double) :: r
    call random_number(r)
    if (r < 0.5d0) then
       v%x = -1.d0
    else
       v%x = 1.d0
    endif
    v%y = 0.d0
    v%z = 0.d0
  end function randomUnitVector2




end module timeDep_mod
