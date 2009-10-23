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
  use lucy_mod, only : setbiasontau
  use math_mod
  use gridio_mod

  implicit none


  private

  public :: runtimeDependentRT, timeDependentRTtest

  type STACKTYPE
     real(double) :: freq(10000000)
     real(double) :: eps(1000000)
     type(VECTOR) :: direction(10000000)
     type(VECTOR) :: position(10000000)
     logical :: photonFromSource(1000000)
     logical :: photonFromGas(10000000)
     logical :: beenScattered(10000000)
     logical :: radiativeEquPhoton(10000000)
  end type STACKTYPE

contains

  subroutine runTimeDependentRT(grid, source, nSource, nLambda, xArray)
    use input_variables, only : teff
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nLambda
    integer :: nStack
    real :: xArray(:)
    integer :: iter, itime
    type(STACKTYPE) :: oldStack
    real(double) :: deltaT, currentTime, deltaTmax, deltaTmin, varyUntilTime, startVaryTime
    real(double) :: newDeltaT
    integer :: nMonte
    character(len=80) :: vtkFilename
    integer :: iDump
    real(double) :: tDump, nextDumpTime, dumpFromNow
    real(double) :: photonsPerSecond
    real(double) :: observerDistance 
    logical :: dumpNow
    character(len=80) :: oldStackFilename, currentStackFilename
    real(double) :: luminosityPeriod
    logical :: varyingSource
    real(double) :: photonsPerStep
    real(double) :: gridCrossingTime
    integer, parameter :: nSedWavelength = 100, nTime = 1000
    real(double) :: sedTime(nTime)
    integer :: i
    real(double) :: sedWavelength(nSedWavelength)
    real(double) :: outputFlux(nSedWavelength, nTime)
    real(double) :: outputFluxScat(nSedWavelength, nTime)
    real(double) :: sedFlux(nSedWavelength, nTime)
    real(double) :: sedFluxScat(nSedWavelength, nTime)
    real(double) :: sedFluxStep(nSedWavelength, nTime)
    real(double) :: sedFluxScatStep(nSedWavelength, nTime)
    real(double) :: w1, w2
    real(double) :: sourceLuminosity, accretionLuminosity, tAcc, frac, accretionArea
    real(double) :: inc, endTime
    type(VECTOR) :: observerDirection, observerposition
    logical :: lastTime, ok
    logical :: seedRun

    seedRun = .false.
    inc = 60.d0 * degToRad
    observerDistance = 140.d0 * pcToCm
    observerDirection = VECTOR(sin(inc), 0.d0, cos(inc))
    observerPosition = observerDistance * observerDirection


    oldStackFilename = "stack1.dat"
    currentStackFilename  = "stack2.dat"
#ifdef MPI
    write(oldStackFilename, '(a,i3.3,a)') "stack1_",myrankGlobal,".dat"
    write(currentStackFilename, '(a,i3.3,a)') "stack2_",myrankGlobal,".dat"
    call init_random_seed()
#endif


    dumpFromNow = 0.d0
    sedFlux = 0.d0
    sedFluxScat = 0.d0
    w1 = 1200.d0
    w2 = 2.d7
    do i = 1, nSedWavelength
       sedWavelength(i) = log10(w1) + (log10(w2)-log10(w1))*dble(i-1)/dble(nSedWavelength-1)
    enddo
    sedWavelength = 10.d0**sedWavelength


    nStack = 0
    tDump = 1.d9
    deltaTmin = 1.d9
    Deltatmax = 1.d9
    gridCrossingTime = grid%octreeRoot%subcellSize*1.d10/cSpeed
    nMonte = 1000000

    photonsPerStep = 1000000 ! dble(nMonte) / (gridCrossingTime / deltaTmin)

    photonsPerSecond = photonsPerStep / deltaTmin

    deltaTmin = 1.d9

    iter = 0
    deltaT = DeltaTmin

    call zeroTemperature(grid%octreeRoot)

    varyingSource = .true.
    varyUntilTime = 1.d30
    endTime = 1.d30
    startVaryTime = 0.d0
!    call clearDust(grid%octreeRoot)


    currentTime = 0.d0
    iDump = 0
    nextDumpTime = currentTime + tDump

!    call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
!    call readAMRgrid("output.safe", .false., grid)
!    call allocateMemoryForTimeDep(grid%octreeRoot)
    if (varyingSource) then
       call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
       call readAMRgrid("output.safe", .false., grid)
       call allocateMemoryForTimeDep(grid%octreeRoot)
       call calculateUdensFromTemperature(grid%octreeRoot)
       i = findIlambda(1.e5, xArray, nLambda, ok)
       call setBiasOnTau(grid, i)
       luminosityPeriod = 3600.d0
       varyUntilTime = 10.d0 * luminosityPeriod
       startVaryTime = 0.d0
       seedRun = .true.
      
       deltaTMax = (varyUntilTime-startVaryTime)/dble(nTime)
       deltaTMin = deltaTmax
       tDump = deltaTmax
       dumpFromNow =  observerDistance/cSpeed
       endTime  = varyUntilTime
       nextDumpTime = dumpFromNow
       deltaT = deltaTmax
       do i = 1, nTime
          sedTime(i) = dumpFromNow + varyUntilTime * dble(i-1)/dble(nTime-1)
       enddo

    endif



    lastTime = .false.

    do while (currentTime < varyUntilTime)
       iter = iter + 1
       if (myrankGlobal == 0) then
!          write(vtkFilename, '(a,i4.4,a)') "iter",iter,".vtk"
!          call writeVtkFile(grid, vtkfilename, &
!               valueTypeString=(/"rho        ", "temperature", "edens_g    ","edens_s    ", "bias     "/))
!          write(vtkFilename, '(a,i4.4,a)') "temperature",iter,".dat"
!          call writeValues(vtkFilename, grid, currentTime)
!
!          write(vtkFilename, '(a,i4.4,a)') "iter",iter,".grid"
!          call writeAMRgrid(vtkFilename, .false., grid)

       endif



       if (varyingSource.and.(currentTime  < varyUntilTime)) then

          if (currentTime >= startVaryTime) then
             
             sourceLuminosity = fourPi * stefanBoltz * (source(1)%radius * 1.d10)**2 * teff**4
             if (.not.seedRun) then
                accretionLuminosity = sourceLuminosity * (1.d0+cos(pi+twoPi*currentTime/luminosityPeriod))/2.d0
             else
                accretionLuminosity = 0.d0
             endif
             accretionArea = 5.d-2 * fourPi * (source(1)%radius * 1.d10)**2
             tAcc  = (accretionLuminosity / (stefanBoltz * accretionArea))**0.25d0
             frac = 5.d-2
             write(*,*) "accretion temp ", tacc
             call fillSpectrumBB(source(1)%spectrum, dble(teff), 1200.d0, 2.d7, 200)
             call addToSpectrumBB(source(1)%spectrum, tAcc, frac)
             call normalizedSpectrum(source(1)%spectrum)
             source(1)%luminosity = sourceLuminosity + accretionLuminosity


             photonsPerStep = 1000000 ! dble(nMonte) / (gridCrossingTime / deltaTmin)
!             photonsPerSecond = photonsPerStep / deltaT
          else
             photonsPerStep = 1000000 ! dble(nMonte) / (gridCrossingTime / deltaTmin)
!             photonsPerSecond = photonsPerStep / deltaT

          endif
       else
          photonsPerStep = max(10000000, nint(5000000.d0 * min(1.d0, deltaT/gridCrossingTime)))
          Deltatmax = 1.d12
          deltaTmin = 1.d9
       endif

       dumpNow = .false.
       if ((currentTime+deltaT) >= nextDumpTime) then
          write(*,*) "deltaT trimmed ",deltaT, currentTime, nextDumpTime
          deltaT = nextDumpTime - currentTime
          dumpNow = .true.
       endif


       if (varyingSource.and.seedRun) then
          photonsPerStep = 10000000
          deltaT = 1.e10
       endif


       nMonte = photonsPerStep ! min(1000000,photonsPerSecond * deltaT)

       sedFluxStep  = 0.d0
       sedFluxScatStep = 0.d0
       if (myrankGlobal == 0) write(*,*) iter," Calling RT with timestep of ", deltaT, currentTime,nMonte
       call  timeDependentRTStep(grid, oldStack, nStack, oldStackFilename, currentStackFilename, &
            nsource, source, deltaT, newDeltaT, deltaTmax, deltaTmin, &
            nMonte, xArray, nLambda, varyingSource, currentTime, &
            nSedWavelength, nTime, sedWavelength, sedTime, sedFluxStep, sedFluxScatStep, dumpFromNow, &
            observerPosition, observerDirection, lastTime, seedRun)

       sedFlux  = sedFlux + sedFluxStep
       sedFluxScat = sedFluxScat + sedFluxScatStep

       currentTime = currentTime + deltaT

       if (seedRun) then
          seedRun = .false.
          currentTime = 0.d0


          do iTime = 1, nTime-1
             sedFlux(1:nSedWavelength, itime) = (sedFlux(1:nSedWavelength, itime) / deltaT) &
                  * (sedTime(itime+1)-sedTime(itime))
             sedFluxScat(1:nSedWavelength, itime) = (sedFluxScat(1:nSedWavelength, itime) / deltaT) &
                  * (sedTime(itime+1)-sedTime(itime))
          enddo


       endif

       deltaT = newDeltaT


       if (myrankGlobal == 0) then
             do i = 1, nSedWavelength-1
                outputFlux(i,1:nTime) = sedFlux(i,1:nTime) / &
                     (sedWavelength(i+1)-sedWavelength(i))
             enddo

             do i = 1, nSedWavelength-1
                outputFluxScat(i,1:nTime) = sedFluxScat(i,1:nTime) / &
                     (sedWavelength(i+1)-sedWavelength(i))
             enddo

          do itime = 1, nTime-1
             write(vtkFilename, '(a,i4.4,a)') "sed",itime,".dat"
             open(32, file=vtkFilename, status="unknown", form="formatted")
             write(32,'(a,1pe13.5,a)') "# ",dble(itime-1)*deltaT, " seconds"
             do i = 1, nSedWavelength-1
                write(32,'(1p2e14.5)') 0.5d0*(sedWavelength(i)+sedWavelength(i+1)), &
                     outputFlux(i, itime)/(sedTime(itime+1)-sedTime(itime))/observerDistance**2
             enddo
             close(32)
             
             write(vtkFilename, '(a,i4.4,a)') "scat",itime,".dat"
             open(32, file=vtkFilename, status="unknown", form="formatted")
             write(32,'(a,1pe13.5,a)') "# ",dble(itime-1)*deltaT, " seconds"
             
             do i = 1, nSedWavelength-1
                write(32,'(1p2e14.5)') 0.5d0*(sedWavelength(i)+sedWavelength(i+1)), &
                     outputFluxScat(i, iter)/(sedTime(itime+1)-sedTime(itime))/observerDistance**2
             enddo
             close(32)
          enddo
       endif

       if (dumpNow) then
          nextDumpTime = nextDumpTime + tDump
          iDump = idump + 1
!          if (myrankGlobal == 0) then
!             write(vtkFilename, '(a,i4.4,a)') "output",idump,".vtk"
!             call writeVtkFile(grid, vtkfilename, &
!                  valueTypeString=(/"rho        ", "temperature", "edens_g    "/))
!             
!             write(vtkFilename, '(a,i4.4,a)') "output",idump,".grid"
!             call writeAMRgrid(vtkFilename, .false., grid)
!          endif


       endif
    enddo
  end subroutine runTimeDependentRT

  subroutine timeDependentRTStep(grid, oldStack, oldStacknStack, oldStackFilename, currentStackFilename, &
       nsource, source, deltaT, newDeltaT, deltaTmax, deltaTmin, nMonte, lamArray, nLambda, varyingSource, currentTime, &
       nSedWavelength, nTime, sedWavelength, sedTime, sedFlux, sedFluxScat, dumpFromNow, observerPosition, observerDirection, &
       lastTime, seedRun)
    type(GRIDTYPE) :: grid
    integer :: nSource
    logical :: varyingSource
    integer :: nLambda
    integer :: nFreq
    logical :: lastTime
    integer :: nSedWavelength, nTime
    real(double) :: sedWavelength(:)
    real(double) :: sedTime(:)
    real(double) :: sedFlux(nSedWavelength,  nTime)
    real(double) :: sedFluxScat(nSedWavelength,  nTime)
    real(double) :: currentTime, dumpFromNow
    integer :: oldStacknStack
    integer :: currentStacknStack
    character(len=*) :: oldStackFilename, currentStackFilename
    character(len=80) :: tempFilename
    type(STACKTYPE), intent(inout) :: oldStack
    type(STACKTYPE) :: currentStack
    real :: lamArray(:)
    real(double) :: deltaT, deltaTmax, deltaTmin
    type(SOURCETYPE) :: source(:)
    logical :: photonFromSource, photonFromGas
    integer :: nMonte, iMonte, nFromGas, nFromSource
    integer :: i
    real(double) :: freqArray(2000), dnu(2000), kAbsArray(2000), kAbsArray2(2000)
    real(double) :: prob(2000), luminosity, sourceLuminosity
    integer :: nFromMatter, nPhotons
    real(double) :: fracSource, chanceSource, chanceGas, weightSource, weightGas
    type(VECTOR) :: rVec, uHat, rHat
    type(VECTOR) :: observerDirection, observerPosition
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
    integer :: oldStackUnit, currentStackUnit
    logical :: absorbed, scattered, outOfTime, finished, ok
    real(double) :: kappaAbs, kappaSca, albedo
    real(double), intent(out) :: newDeltaT
    real(double) :: fac
    real(double) :: totalLineEmission, totalContEmission, t1, t2, t3, t4, t5, t6, checkLum, t7
    real(double) :: checkLumSource
    real(double) :: photonTagTime
    real(double) :: firstObserverTime, timeToObserver, lastObserverTime
    real :: lambda0
    logical :: useFileForStack
    logical :: beenScattered
    logical :: radiativeEquPhoton, seedRun
    integer :: nEscaped
    integer :: nFromMatterThisThread
#ifdef MPI
    include 'mpif.h'
    integer :: ierr
    real(double) :: tempDouble
    real(double), allocatable :: tempDoubleArray(:)
#endif
    firstObserverTime = sedTime(1)
    lastObserverTime = sedTime(nTime)

    useFileForStack = .true.

    lambda0 = 0.

    oldStackUnit = 41
    currentStackUnit = 42

    if (useFileForStack) then
       open(oldStackUnit, file=oldStackFilename, form="unformatted")
       open(currentStackUnit, file=currentStackFilename, form="unformatted")
    endif

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
    call calculateEtaContLocal(grid, grid%octreeRoot)
!    call calculateEtaContBias(grid, grid%octreeRoot)
    call computeProbDist(grid, totalLineEmission, totalContEmission, lambda0, .false.)
    luminosity = 0.d0
    call calculateGasEmissivity(Grid%octreeRoot, luminosity)
    currentStacknStack =  0


    sourceLuminosity = SUM(source(1:nSource)%luminosity)

    nFromMatter = nMonte
    nFromMatterThisThread = nFromMatter
#ifdef MPI
    nFromMatterThisThread = dble(nFromMatter) / dble(nThreadsGlobal)
#endif


    nPhotons = oldStacknStack + nFromMatterThisThread

    fracSource = 0.1d0

    Chancesource = sourceLuminosity / (sourceLuminosity + luminosity)
    chanceGas = 1.d0-chanceSource
    weightSource = chanceSource / fracSource
    weightGas = (1.d0-chanceSource) / (1.d0-fracSource)


    weightGas = 1.d0
    weightSource = 1.d0
    fracSource =  sourceLuminosity / (sourceLuminosity + luminosity)

    nFromSource = fracSource * nFromMatter
    nFromGas = nFromMatter - nFromSource

    if (myrankGlobal == 0) then
       write(*,*) "doing loop with ",nphotons, " photons"
       write(*,*) "source lum ",sourceLuminosity
       write(*,*) "matter lum ",luminosity
       write(*,*) "frac source ",fracSource
    endif

    checkLum = 0.d0
    checkLumSource = 0.d0
    nEscaped = 0

    do iMonte = 1, nPhotons
!       write(*,*) "imonte ", imonte, nPhotons 
       if (oldStacknStack > 0)  then
          call getPhotonFromStack(useFileForStack, oldStack, oldStackUnit, &
               oldStacknStack, rVec, uHat, eps, freq, &
               photonFromSource, photonFromGas, beenScattered, radiativeEquPhoton)
             wavelength = (cSpeed/freq)*1.d8
          photonTime = 0.d0
       else
          call random_number(r)

          if (r < fracSource) then

             radiativeEquPhoton = .false.
             call randomSource(source, nSource, iSource)
             thisSource = source(iSource)
             call getPhotonPositionDirection(thisSource, rVec, uHat, rHat,grid)

             call getWavelength(thisSource%spectrum, wavelength)
             freq = cSpeed/(wavelength / 1.e8)
             photonTime = 0.d0
             eps = weightSource * sourceluminosity * deltaT / dble(nFromSource)
!             write(*,*) "star eps " , eps
             photonFromSource = .true.
             photonFromGas = .false.
             checkLumSource = checkLumSource + eps

             if (.not.seedRun) then
                call random_number(r)
                photonTime = r*deltaT
             endif

             ilambda = findIlambda(real(wavelength), lamArray, nLambda, ok)
             if (varyingSource.and.(.not.seedRun)) then
                call tauAlongPath(ilambda, grid, rVec, observerDirection, tau)
                timeToObserver = distanceToObserver(rVec, observerPosition, observerDirection)/cspeed
                photonTagTime = currentTime + photonTime + timeToObserver
                call timeBinPhoton(photonTagTime, wavelength, eps, tau, nSedWavelength, nTime, &
                  sedTime, sedWavelength, sedFlux, sedFluxScat, beenScattered)
             endif



          else

             call random_number(r)
             thisOctal => grid%octreeRoot
             call locateContProbAMR(r,thisOctal,subcell)
             rVec = randomPositionInCell(thisOctal, subcell)
             radiativeEquPhoton = .false.
             if (modulus(rVec)*1.d10 > cSpeed * currentTime) then
                radiativeEquPhoton = .true.
             endif
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
             eps = (1.d0/thisOctal%biasCont3D(subcell)) * weightGas * luminosity * deltaT / dble(nFromGas)
!             write(*,*) "gas eps ",eps, luminosity, weightgas, thisOctal%biasCont3d(subcell)
!             write(*,*) "Wavelength ", wavelength/1.e4
             photonFromGas = .true.
             photonFromSource = .false.
             photonTime = 0.d0
             uHat = randomUnitVector()

             if (.not.seedRun) then
                call random_number(r)
                photonTime = r*deltaT
             endif

             ilambda = findIlambda(real(wavelength), lamArray, nLambda, ok)
             if (varyingSource) then
                call tauAlongPath(ilambda, grid, rVec, observerDirection, tau)
                timeToObserver = distanceToObserver(rVec, observerPosition, observerDirection)/cspeed
                
                photonTagTime = currentTime + photonTime + timeToObserver

                if (seedRun) then
                   if (photonTagTime > lastObserverTime) then
                      do i = 1, nTime
                         call timeBinPhoton(sedTime(i), wavelength, eps, tau, nSedWavelength, nTime, &
                              sedTime, sedWavelength, sedFlux, sedFluxScat, beenScattered)
                      enddo
                   endif
                   if (photonTagTime < firstObserverTime) then
                      do i = 1, nTime
                         call timeBinPhoton(sedTime(i), wavelength, eps, tau, nSedWavelength, nTime, &
                              sedTime, sedWavelength, sedFlux, sedFluxScat, beenScattered)
                      enddo
                   endif
                endif
                if (.not.seedRun) then
                   call timeBinPhoton(photonTagTime, wavelength, eps, tau, nSedWavelength, nTime, &
                        sedTime, sedWavelength, sedFlux, sedFluxScat, beenScattered)
                endif
             endif
                
          endif
          checkLum = checkLum + eps
          beenScattered = .false.
       endif
       ilambda = findIlambda(real(wavelength), lamArray, nLambda, ok)
       call random_number(r)

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

          timeToBoundary = (distToBoundary*1.d10) / photonSpeed

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
                beenScattered = .true.
                if (varyingSource) then
                   call tauAlongPath(ilambda, grid, rVec, observerDirection, tau)
                   timeToObserver = distanceToObserver(rVec, observerPosition, observerDirection)/cspeed
                   photonTagTime = currentTime + photonTime + timeToObserver
                   if (seedRun) then
                      if (photonTagTime < firstObserverTime) then
                         do i = 1, nTime
                            call timeBinPhoton(sedTime(i), wavelength, eps, tau, nSedWavelength, nTime, &
                                 sedTime, sedWavelength, sedFlux, sedFluxScat, beenScattered)
                         enddo
                      endif
                      if (photonTagTime > lastObserverTime) then
                         do i = 1, nTime
                            call timeBinPhoton(sedTime(i), wavelength, eps, tau, nSedWavelength, nTime, &
                                 sedTime, sedWavelength, sedFlux, sedFluxScat, beenScattered)
                         enddo
                      endif
                   endif
                   if (.not.seedRun) then
                      call timeBinPhoton(photonTagTime, wavelength, eps, tau, nSedWavelength, nTime, &
                           sedTime, sedWavelength, sedFlux, sedFluxScat, beenScattered)
                   endif

                endif

             endif
          endif

          if (absorbed.or.scattered) then
             if ((photonTime + (distanceToEvent*1.d10)/photonSpeed) > deltaT) then
                absorbed = .false.
                scattered = .false.
                distanceToEvent = ((deltaT - photonTime) * photonSpeed)/1.d10
                outOfTime = .true.
                finished = .true.
             endif
          else
             if ((photonTime + (distanceToEvent*1.d10)/photonSpeed) > deltaT) then
                distanceToEvent = ((deltaT - photonTime) * photonSpeed)/1.d10
                outOfTime = .true.
                finished = .true.
             endif
          endif

          photonTime = photonTime + (distanceToEvent*1.d10) / photonSpeed

          call updateDistanceGrids(thisOctal, subcell, kappaAbs, eps, distanceToEvent, photonFromGas, &
               photonFromSource)

          rVec = rVec + (distanceToEvent + 1.d-3*grid%halfSmallestSubcell) * uHat
          

          if (.not.inOctal(grid%octreeRoot, rVec)) then
             finished = .true.
             outOfTime = .false.
             nEscaped = nEscaped + 1
          endif
          if (scattered) then
             uhat = randomUnitVector()  ! ISOTROPIC SCATTERING ONLY SO FAR
             scattered = .false.
             photonFromGas = .true.
             photonFromSource = .false.
          endif
       end do
       if (outOfTime) then
          call addPhotonToStack(useFileForStack, currentStack, currentStackUnit, &
               currentStackNStack, rVec, uHat, &
               eps, freq, photonFromSource, photonFromGas, beenScattered, radiativeEquPhoton)
       endif
    end do

#ifdef MPI
    write(*,*) myrankGlobal, " updating grid"
    call updateGridMPI(grid)

    tempdouble = 0.d0
     call MPI_ALLREDUCE(checkLum,tempDouble,1,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
     checkLum = tempDouble

    tempdouble = 0.d0
     call MPI_ALLREDUCE(checkLumSource,tempDouble,1,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
     checkLumSource = tempDouble

     allocate(tempDoubleArray(1:nSedWavelength))

     do i = 1, nTime
        tempDoubleArray = 0.d0
        call MPI_ALLREDUCE(sedFlux(1:nSedWavelength, i),tempDoubleArray,&
             nSedWavelength,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
        sedFlux(1:nSedWavelength,i) = tempDoubleArray
     enddo

     do i = 1, nTime
        tempDoubleArray = 0.d0
        call MPI_ALLREDUCE(sedFluxScat(1:nSedWavelength, i),tempDoubleArray,&
             nSedWavelength,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
        sedFluxScat(1:nSedWavelength,i) = tempDoubleArray
     enddo

     deallocate(tempDoubleArray)



#endif



     if (myrankGlobal == 0) then
        write(*,*) "Sanity check for luminosity ", checkLum/deltaT, luminosity+sourceLuminosity
        write(*,*) "Sanity check for source luminosity ", checkLumSource/deltaT, sourceLuminosity
     endif

    if (varyingSource) then
       call calculateADot(grid%octreeRoot, deltaT, currentTime)
    else
       call calculateADot(grid%octreeRoot, deltaT, 1.d30)
    endif
       


    call calculateEnergyDensity(grid%octreeRoot, deltaT, photonSpeed)

    if (.not.seedRun) then
       call updateUDens(grid%octreeRoot, deltaT, grid)
    endif

    call calculateTemperatureFromUdens(grid%octreeRoot)
    call calculateEtaContLocal(grid, grid%octreeRoot)

    fac = 0.3d0


    newDeltaT = 1.d30
    call calculateNewDeltaT(grid, grid%octreeRoot, newDeltaT, t1, t2, t3, t4, t5, t6, t7)
    if (myrankGlobal == 0) then
       write(*,*) "Temperature for time controlling cell is ",t1
       write(*,*) "Udens for time controlling cell is ",t2
       write(*,*) "adot for time controlling cell is ",t3
       write(*,*) "eta for time controlling cell is ",t4
       write(*,*) "rho for time controlling cell is ",t5
       write(*,*) "cooling time for controlling cell is ",fac*t2/abs(t4-t3)
       write(*,*) "equilibrium time for controlling cell is ",t6
       write(*,*) "new delta T ",newDeltaT*fac
       write(*,*) "chiline ",t7
    endif
    newDeltaT = newDeltaT  * fac

    newdeltaT = min(newDeltaT, deltaTMax)
    newdeltaT = max(newdeltaT, deltaTmin)



    oldStack  = currentStack
    if (useFileForStack) then
       close(oldStackUnit)
       close(currentStackUnit)
       tempFilename = oldStackFilename
       oldStackFilename = currentStackFilename
       currentStackFilename = tempFilename
    endif
    oldStacknStack = currentStacknStack
    currentStacknStack = 0
  end subroutine timeDependentRTStep



  subroutine writeValues(outputFilename, grid, currentTime)
    type(GRIDTYPE) :: grid
    real(double) :: currentTime, tVal
    character(len=*) :: outputFilename
    type(VECTOR) :: rVec, uHat, cVec
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    open(21, file=outputFilename, status="unknown", form="formatted")
    write(21, *) "# Current time ",currentTime
    rVec = VECTOR(0.d0, 0.d0, 0.d0)
    uHat = VECTOR(1.d0, 0.d0, 0.d0)
    thisOctal => grid%octreeRoot
    call findSubcellLocal(rVec, thisOctal, subcell)
    do while(inOctal(grid%octreeRoot, rVec))
       call findSubcellLocal(rVec, thisOctal, subcell)
       cVec = subcellCentre(thisOctal, subcell)
       write(21, *) modulus(cVec), thisOctal%temperature(subcell)
       call distanceToCellBoundary(grid, rVec, uHat, tVal, sOctal=thisOctal)
       rVec = rVec + (tVal + 1.d-3*grid%halfSmallestSubcell) * uHat
    enddo
    close (21)
  end subroutine writeValues
       

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

  recursive subroutine updateUdens(thisOctal, deltaT, grid)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: deltaT, equilibriumTime,deltaUdens, newUDens
    real(double) :: fac, oldUdens, t
    real :: kappaP
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call updateUdens(child, deltaT, grid)
                exit
             end if
          end do
       else
          thisOctal%chiLine(subcell) = 0.d0
          rVec = subcellCentre(thisOctal, subcell)
          equilibriumTime = 1.d30
          if (thisOctal%aDot(subcell) > 0.d0) then
             newudens = 1.d30
             oldUdens = thisOctal%udens(subcell)
             t = thisOctal%temperature(subcell)
             fac = 1.d30
             thisOctal%temperature(subcell) = 100.
             call returnKappa(grid, thisOctal, subcell, kappap=kappap)
             do while (fac > 1.d-5)
                thisOctal%temperature(subcell) = &
                     max(0.d0,(thisOctal%aDot(subcell) / (4.d0 * stefanBoltz * kappaP))**0.25d0) ! in radiative equilibrium
                newUdens =  uDensFunc(dble(thisOctal%temperature(subcell)), thisOctal%rho(subcell),  7.d0/5.d0, 2.33d0)
                call returnKappa(grid, thisOctal, subcell, kappap=kappap)
                fac = abs(newUdens - oldUdens)/newUdens
                oldUdens = newUdens
             enddo
             thisOctal%temperature(subcell) = t
             deltaUdens = abs(newUdens - thisOctal%uDens(subcell))
             if (thisOctal%adot(subcell) > thisOctal%etaCont(subcell)) then
                equilibriumTime = deltaUdens/thisOctal%aDot(subcell)
             else
                equilibriumTime = deltaUdens/thisOctal%etaCont(subcell)
             endif
          endif
          if ((equilibriumTime < deltaT).and.(thisOctal%adot(subcell) > 0.d0)) then
!            write(*,*) "cell in eq ", equilibriumTime, thisOctal%adot(subcell), teq
!            oldUdens = thisOctal%uDens(subcell)
!            fac = 1.d30
!            do while (fac > 1.d-5)
!               thisOctal%temperature(subcell) = max(0.d0,(thisOctal%aDot(subcell) / (4.d0 * stefanBoltz * kappaP))**0.25d0) ! in radiative equilibrium
!               thisOctal%uDens(subcell) = uDensFunc(dble(thisOctal%temperature(subcell)), thisOctal%rho(subcell),  7.d0/5.d0, 2.33d0)
!               call returnKappa(grid, thisOctal, subcell, kappap=kappap)
!               fac = abs(oldUdens - thisOctal%udens(subcell))/thisOctal%udens(subcell)
!               oldUdens = thisOctal%uDens(subcell)
!            enddo
             thisOctal%uDens(subcell) = newUdens
             thisOctal%etaCont(subcell) = thisOctal%aDot(subcell)
             thisOctal%chiLine(subcell) = 1.d10
          else
             thisOctal%udens(subcell) = thisOctal%udens(subcell) + deltaT * & ! otherwise update udens
                  (thisOctal%adot(subcell) - thisOctal%etaCont(subcell))
             if (thisOctal%udens(subcell) < 0.d0) then
                write(*,*) "neg udens ",thisOctal%adot(subcell), thisOctal%etaCont(subcell), deltaT, thisOctal%temperature(subcell)
                thisOctal%udens(subcell) = 0.d0
             endif
          endif
       endif
    Enddo
  end subroutine updateUdens

  recursive subroutine calculateEtaContLocal(grid, thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real :: kappap
    integer :: subcell, i

    Do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateEtaContLocal(grid, child)
                exit
             end if
          end do
       else
          call returnKappa(grid, thisOctal, subcell, kappap=kappap)
          thisOctal%etaCont(subcell) = fourPi * &
               kappap * (stefanBoltz/pi) * thisOctal%temperature(subcell)**4
       endif
    enddo
  end subroutine calculateEtaContLocal

  recursive subroutine calculateEtaContBias(grid, thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateEtaContBias(grid, child)
                exit
             end if
          end do
       else
          if ((thisOctal%rho(subcell) > 1.d-30).and.(thisOctal%etaCont(subcell) > 0.d0)) then
             thisOctal%biasCont3D(subcell) = 1.d0/sqrt(thisOctal%etaCont(subcell))
          else
             thisOctal%biasCont3D(subcell) = 1.d0
          endif
          thisOctal%biasCont3D(subcell) = 1.d0
          
       endif
    enddo
  end subroutine calculateEtaContBias

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
               max(0.d0,temperatureFunc(thisOctal%uDens(subcell), thisOctal%rho(subcell), 7.d0/5.d0, 2.33d0))
       endif
    enddo
  end subroutine calculateTemperatureFromUdens

  recursive subroutine calculateUdensFromTemperature(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateUdensFromTemperature(child)
                exit
             end if
          end do
       else
          thisOctal%udens(subcell) = &
               uDensFunc(dble(thisOctal%temperature(subcell)), thisOctal%rho(subcell), 7.d0/5.d0, 2.33d0)
       endif
    enddo
  end subroutine calculateUdensFromTemperature

  recursive subroutine calculateNewDeltaT(grid, thisOctal, deltaT, temp, udens, adot, eta, rho, equilibriumTime, chiline)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(VECTOR) :: rVec
    real(double),intent(inout) :: deltaT, temp, udens, adot, eta, rho, equilibriumTime, chiline
    integer :: subcell, i
    real(double) :: t

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateNewDeltaT(grid, child, deltaT, temp, udens, adot, eta, rho, equilibriumtime, chiline)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)

          if (thisOctal%etaCont(subcell) >= thisOctal%aDot(subcell)) then ! cooling
             t = thisOctal%uDens(subcell)/abs(thisOctal%aDot(subcell)-thisOctal%etaCont(subcell))
             if ((t < deltaT).and.(thisOctal%temperature(subcell) > 0.d0)) then
                deltaT = t 
                temp = thisOctal%temperature(subcell)
                udens = thisOctal%udens(subcell)
                adot = thisOctal%adot(subcell)
                eta = thisOctal%etaCont(subcell)
                rho = thisOctal%rho(subcell)
                chiline = 0.d0
             endif
          endif

!          if ((thisOctal%etaCont(subcell) < thisOctal%aDot(subcell)).and.(thisOctal%chiLine(subcell)==0.d0)) then ! heating
!             newudens = 1.d30
!             oldUdens = thisOctal%udens(subcell)
!             t = thisOctal%temperature(subcell)
!             fac = 1.d30
!             thisOctal%temperature(subcell) = 100.
!             call returnKappa(grid, thisOctal, subcell, kappap=kappap)
!             do while (fac > 1.d-5)
!                thisOctal%temperature(subcell) = &
!                     max(0.d0,(thisOctal%aDot(subcell) / (4.d0 * stefanBoltz * kappaP))**0.25d0) ! in radiative equilibrium
!                newUdens =  uDensFunc(dble(thisOctal%temperature(subcell)), thisOctal%rho(subcell),  7.d0/5.d0, 2.33d0)
!                call returnKappa(grid, thisOctal, subcell, kappap=kappap)
!                fac = abs(newUdens - oldUdens)/newUdens
!                oldUdens = newUdens
!             enddo
!             thisOctal%temperature(subcell) = t
!             deltaUdens = abs(newUdens - thisOctal%uDens(subcell))
!             t = deltaUdens/abs(thisOctal%aDot(subcell)-thisOctal%etaCont(subcell))
!             if (t < deltaT) then
!!             write(*,*) "time ",t,currentTemp, newTemp, newUdens, thisOctal%udens(subcell), thisOctal%adot(subcell)
!                deltaT = t 
!                temp = thisOctal%temperature(subcell)
!                udens = thisOctal%udens(subcell)
!                adot = thisOctal%adot(subcell)
!                eta = thisOctal%etaCont(subcell)
!                rho = thisOctal%rho(subcell)
!                equilibriumTime = t
!                chiline = thisOctal%chiLine(subcell)
!             endif
!          endif
       endif
    enddo
  end subroutine calculateNewDeltaT

  recursive subroutine calculateAdot(thisOctal, deltaT, currentTime)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: deltaT, v, currentTime
    integer :: subcell, i
    type(VECTOR) :: rVec

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateAdot(child, deltaT, currentTime)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          if (modulus(rVec) < cSpeed * (currentTime + deltaT) /1.d10) then ! can source photons have reached this cell?
             v = cellVolume(thisOctal, subcell) * 1.d30
             thisOctal%aDot(subcell) = max(0.d0,(1.d0 / v) * thisOctal%distancegridAdot(subcell) / deltaT)
          else
             thisOctal%adot(subcell) = thisOctal%etaCont(subcell) ! cell is radiative equilibrium
          endif
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

  recursive subroutine clearDust(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call clearDust(child)
                exit
             end if
          end do
       else
          if (thisOctal%rho(subcell) < 1.d-30) then
             thisOctal%dustTypeFraction(subcell,:) = 1.d-5
          endif
       endif
    enddo
  end subroutine clearDust

  recursive subroutine zeroTemperature(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroTemperature(child)
                exit
             end if
          end do
       else
          thisOctal%temperature(subcell) = 0.
       endif
    enddo
  end subroutine zeroTemperature



  subroutine timeDependentRTtest()
    real(double) :: currentTime, endTime
    integer, parameter :: nx = 101

    integer, parameter :: nProb = nx + 1
    real(double) :: etaCont(nx)
    real(double) :: aDot(nx)
    real(double) :: oldaDot(nx)
    real(double) :: oldetaCont(nx)
    real(double) :: rho(nx)
    real(double) :: kappa(nx)
    real(double) :: temperature(nx)
    real(double) :: distanceGridAdot(nx)
    real(double) :: distanceGridPhotonFromGas(nx)
    real(double) :: distanceGridPhotonFromSource(nx)
    logical :: photonFromSource, photonFromGas
    real(double) :: prob(nProb)
    real(double) :: bias(nProb)
    real(double) :: photonProb(nProb)
    real(double) :: photonEnergy
    real(double) :: xCen(nx)
    real(double) :: xArray(nprob)
    real(double) :: uDens(nx)
    real(double) :: photonArray(nx)
    real(double) :: photonEnergyDensityFromSource(nx)
    real(double) :: photonEnergyDensityFromGas(nx)
    real(double) :: photonEnergyDensity(nx)
    real(double) :: oldphotonEnergyDensity(nx)
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
    logical :: outOFtime, scattered, dumpnow
    logical :: timeBoundary, spaceBoundary, reflecting
    type(VECTOR) :: rVec, uHat
    integer :: i, iMonte, nMonte, iPos,iter, nPhotons
    character(len=30) :: outfile
    real(double) :: photonSpeed
    real(double) :: sourcePeriod, sourceLuminosity
    real(double) :: fracSource
    real(double) :: fac, meanFreePath, probNewPhoton, teq, deltaudens, newudens
    real(double) :: tdump, timeofnextdump, photonBias, Be
    real(double) :: udens_n_plus_1
    integer :: nFromMatter, nFromGas, nFromStar

    reflecting = .true.
    fracSource = 0.1d0
    udensAnalytical  = 0.d0
    udens = 1.d-30
    photonSpeed = cSpeed
    nMonte = 10000000
    tauBox = 1.d2
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
!       temperature(i) = temperatureFunc(uDens(i),  rho(i), 5.d0/3.d0, 0.6d0)
!    enddo
!    do i = 1, nx
!       t = temperatureFunc(uDens(i), rho(i), 5.d0/3.d0, 0.6d0)
!       write(*,*) "temperature ",t
!!       photonEnergyDensity(i) = max(1.d1,cos(2.d0*xCen(i))*1.d7)!(fourPi * (stefanBoltz/pi) * t**4 / cSpeed) 
!       photonEnergyDensity(i) = (fourPi * (stefanBoltz/pi) * t**4 / photonSpeed) 
!       photonEnergyDensity(i) = (8.d0*pi**5*(kerg*t)**4)/(15.d0 * (hCgs * photonSpeed)**3)
!    enddo
!    deltaT = 1.d-30

! test1

    photonEnergyDensity = 1.d12
    udens = 1.d10
    rho = 1.d-7
    kappa = 4.d-1
    oldNstack = nMonte

! test2

!    photonEnergyDensity = 0.d0
!    udens = 1.d8
!    rho = 1.d100
!    kappa = 4.d6
!    oldnStack = 0

! test3

!    teq = 1000.d0
!    udens = 1.d-20
!    photonEnergyDensity = 1.d-20
!    rho = 1.d-10
!    kappa = 1.d13
!    nmonte = 1000000
!    oldnStack =  0.4* nMonte
!    do i = 1, nx
!       fac = exp(-(0.5-xArray(i))**2 / 0.1d0**2)
!       photonEnergyDensity(i) = 1.d10 * fac
!       Teq = ((photonEnergyDensity(i)*cSpeed/fourPi)*(pi/stefanBoltz))**0.25d0
!       udens(i) = uDensFunc(Teq, rho(i),  5.d0/3.d0, 0.6d0)
!       photonEnergyDensityFromGas(i) = (fourPi/cSpeed) * (stefanBoltz/pi) * teq**4
!    enddo
!    udensAnalytical = udens
!    photonDensAnalytical = photonEnergyDensity
    

! test4 deltaT function (heat kernel)

! for pure scattering -works fine
    udens = 1.d-20
    rho = 1.d-10
    albedo = 1.d0
    nMonte = 0
    oldnStack =  1000000
    kappa = 1.e12

    photonEnergyDensity = 1.d-100
    photonEnergyDensity(nx/2+1) = 1d10/(xArray(2)-xArray(1))
    photonEnergyDensityFromGas = photonEnergyDensity
    photonDensAnalytical = photonEnergyDensity
!    do i = 1 , nx
!       Teq = ((photonEnergyDensity(i)*cSpeed/fourPi)*(pi/stefanBoltz))**0.25d0
!       udens(i) = uDensFunc(Teq, rho(i),  5.d0/3.d0, 0.6d0)
!    enddo
    udensAnalytical = uDens
    
! test5 Sinusoid source into cold wall

    udens = 1.d-20
    albedo = 0.d0
    kappa = 1.d-20
    rho = 2.d-5
    where (xCen >= 0.5) kappa = 1.e6
    nMonte = 40000
    oldnStack = 0
    reflecting = .false.


    call calculateProb(photonProb, bias, xArray, xCen, nProb, photonEnergyDensity)

    write(*,*) "SUM oldepsStack ",SUM(oldEpsStack(1:oldNStack))

    k =  photonspeed/(kappa(1)*rho(1)) 

    write(*,*) "Diffusion coefficent ",k
    diffusionTime=  0.3d0* (dx**2/(2.d0*k))


    deltaT = diffusionTime
    
    diffusionTime = 1.d30
    deltaT = 1.d-13
    tdump = 1.d-13
    timeofnextdump = tDump

    photonEnergy =  SUM(photonEnergyDensity*dx)/dble(nMonte)
    photonArray = 0.d0
    do i = 1,oldnStack
       call getX(xArray, photonprob, bias, nProb, xPhoton, photonBias)
       oldpositionStack(i)%x = xPhoton
       oldPositionStack(i)%y = 0.d0
       oldPositionStack(i)%z = 0.d0

       oldFromGasStack(i) = .true.
       oldFromSourceStack(i) = .false.
       olddirectionStack(i) = randomUnitVector2()!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       oldtimestack(i) = 0.d0
       oldEpsStack(i) = photonEnergy 
       call findArrayIndex(xCen, nx, xPhoton, iPos)
       photonArray(iPos) = photonArray(iPos) + oldEpsStack(i)
    enddo



    currentTime = 0.d0
    endTime = 1.d0
    iter = 0
    currentnStack = 0




    open(69,file="time.dat", status="unknown",form="formatted")

    do while (currentTime < 1.e30)
       meanFreePath = 1.d0/(kappa(1)*rho(1))


!       if ((currentTime  + deltaT) > timeofNextDump) then
!          deltaT = timeofNExtDump - currentTime
!          timeofnextDump = timeofNextDump + tDump
          dumpnow = .true.
!       endif

       distanceGridAdot = 0.d0
       distanceGridPhotonFromSource = 0.d0
       distanceGridPhotonFromGas = 0.d0
       do i = 1, nx
          temperature(i) = temperatureFunc(uDens(i),  rho(i), 5.d0/3.d0, 0.6d0)
!          write(*,*) temperature(i)
       enddo
       etaCont = fourPi * (stefanBoltz/pi) * temperature**4 * kappa * rho 
       call calculateProb(prob, bias, xArray, xCen, nProb, etaCont)
       luminosity = SUM(etaCont*dx)
       write(*,*) "emissivity ",luminosity
       write(*,*) "time step ",deltaT
       write(*,*) "diffusion timestep ", 0.3d0*(dx**2/(2.d0*k))
       write(*,*) "c chi deltat ", cSpeed * kappa(1) * rho(1) * deltaT
       write(*,*) "current Time", currentTime
       if (oldNStack == 0) then
          write(*,*) "energy monte ",SUM(udens*dx)+SUM(photonEnergyDensityFromGas*dx)
       else
          write(*,*) "energy monte ",SUM(udens*dx)+SUM(oldepsStack(1:oldnStack))
          write(*,*) "number of photons ", oldNStack
          write(*,*) "percent in photons ", 100.d0 * &
               SUM(oldepsStack(1:oldnStack))/(SUM(udens*dx)+SUM(oldepsStack(1:oldnstack)))
          write(*,*) "sum oldepsStack ",SUM(oldEpsStack(1:oldnStack))
          write(*,*) "sum udens ", SUM(udens*dx)
       endif
       write(*,*) "energy analy ",SUM(udensanalytical*dx)+SUM(photonDensAnalytical*dx)
       write(*,*) "light distance ", deltaT*photonSpeed
       write(*,*) "mean free path ", 1.d0/(kappa(1)*rho(1))
       write(*,*) "diffusion speed (speed of light units)", ((1.d0/(kappa(1)*rho(1)))/deltaT)/photonSpeed


       nFromMatter = nmonte
!       nFromMatter = max(1000, nMonte-oldNStack)
       nPhotons = oldNStack + nFromMatter
       tot = 0.d0
       sourcePeriod = 1.e-11
       sourceLuminosity = 1.d20 * 0.5d0*(1.d0+cos(twoPi*currentTime/sourcePeriod))
       write(*,*) " source luminosity ", sourceluminosity

       nFromStar = nint(fracSource*dble(nFromMatter))
       nFromGas = nFromMatter - nFromStar
       write(*,*) "nfromstar ",nFromStar, "nFromgas ",nFromgas
       adot = 0.d0
       photonArray = 0.d0
       do iMonte = 1, nPhotons
!          write(*,*) "imonte ", imonte, nPhotons 
          if (oldnStack > 0)  then
             call getPhotonFromStacktest(oldnStack, rVec, uHat, photonTime, epsOverDeltaT, &
                  photonFromSource, photonFromGas, &
                  oldpositionStack, olddirectionStack, oldtimeStack, oldepsStack, oldFromSourceStack, &
                  oldFromGasStack)
             photonTime = 0.d0
          else
             call random_number(r)

             if (r < fracSource) then
                rVec = VECTOR(1.d-18, 0.d0, 0.d0)
                uHat = VECTOR(1.d0, 0.d0, 0.d0)
                photonTime = 0.d0
                epsOverDeltaT = sourceluminosity * deltaT / dble(nFromStar)
                photonFromSource = .true.
                photonFromGas = .false.
                call random_number(r)
                photonTime = (r-0.5d0)*deltaT
             else
                call getX(xArray, prob, bias, nProb, xPhoton, photonBias)
!                write(*,*) xPhoton,photonBias
                call findArrayIndex(xCen, nx, xPhoton, iPos)
                epsOverDeltaT =  photonBias * luminosity * deltaT / dble(nFromGas)
                photonFromGas = .true.
                photonFromSource = .false.
                rVec%x = xPhoton
                rVec%y = 0.d0
                rVec%z = 0.d0
                photonTime = 0.d0
                uHat = randomUnitVector2()!!!!!!!!!!!!!!!!!
                tot = tot + epsOverDeltaT
             endif
             photonTime = 0.d0
          endif

          absorbed = .false.
          scattered = .false.
          finished = .false.
          outOfTime = .false.
          do while (.not.finished)

             timeBoundary = .false.
             spaceBoundary = .false.

             call findArrayIndex(xCen, nx, rVec%x, iPos)
             if (uHat%x > 0.d0) then
                distToBoundary = ((xCen(iPos)+dx/2.d0) - rVec%x)/uHat%x
             else
                distToBoundary = abs( (((xCen(iPos)-dx/2.d0) - rVec%x )/uHat%x ))
             endif

             timeToBoundary = distToBoundary / photonSpeed
             if ((timeToBoundary+photonTime) > deltaT) then
                timeToBoundary = deltaT - photonTime
                distToBoundary = timeToBoundary * photonSpeed
                timeBoundary = .true.
             endif


             tauToBoundary = distToBoundary * kappa(iPos) * rho(iPos)

             call random_number(r)
             tau = -log(1.d0 - r)
!             write(*,*) "tau ", tau, " tautoboundary ", &
!                  tautoboundary, iPos, rVec%x, xCen(iPos)-dx/2.d0, xCen(iPos)+dx/2.d0


             if (tau > tauToBoundary) then
                distanceToEvent = distToBoundary
                absorbed = .false.
                if (.not.timeBoundary) spaceBoundary = .true.
             else
                distanceToEvent = distToBoundary * tau/tauToBoundary
                timeBoundary = .false.
                spaceBoundary = .true.
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

             if (timeBoundary) then
                outOfTime = .true.
                finished = .true.
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
                if (reflecting) then
                   uHat%x = -1.d0 * uHat%x
                else
                   finished = .true.
                endif
             endif
             if (rVec%x > xCen(nx)+dx/2.d0) then
                rVec%x = xCen(nx)+dx/2.d0-1.d-7*dx
                if (reflecting) then
                   uHat%x = -1.d0 * uHat%x
                else
                   finished = .true.
                endif
             endif
             if (scattered) then
                uhat = randomUnitVector2()!!!!!!!!!!
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
          oldaDot = adot
          oldetaCont = etaCont
          oldPhotonEnergyDensity = photonEnergyDensity
          call calculateADottest(distanceGridAdot, aDot, xCen, dx, nx, deltaT)
          call calculatePhotonEnergyDensity(distanceGridPhotonFromGas, &
               photonEnergyDensityFromGas, xCen, dx, nx, deltaT, photonSpeed)
          call calculatePhotonEnergyDensity(distanceGridPhotonFromSource, &
               photonEnergyDensityFromSource, xCen, dx, nx, deltaT, photonSpeed)
          photonEnergyDensity = photonEnergyDensityFromSource + photonEnergyDensityFromGas
          call solveNewUdens(xCen, kappa, rho, photonDensAnalytical, uDensAnalytical, nx, k, deltaT, currentTime)

          write(*,*) "Doing step with deltat= ",deltaT
          do i = 1, nx

             call quarticSub(udens_n_plus_1, photonEnergyDensity(i), udens(i), oldphotonEnergyDensity(i), &
                  kappa(i), rho(i), albedo, 0.6d0, 5.d0/3.d0, deltaT)
             udens(i) = udens_n_plus_1
!             photonEnergyDensity(i) = E_n_plus_1

             Be = cSpeed*photonEnergyDensity(i) / fourPi
             Teq = ((Be * pi)/stefanBoltz)**0.25d0
             newUdens =  uDensFunc(Teq, rho(i),  5.d0/3.d0, 0.6d0)

!             if (udens(i) > 0.d0) write(*,*) "percent from eq " , 100.d0*(newuDens-udens(i))/udens(i)
!             newUdens = udensFunc(Teq,rho(i), 5.d0/3.d0, 0.6d0)
!!             Teq = max(0.d0,(aDot(i) / (4.d0 * stefanBoltz * kappa(i) *rho(i)))**0.25d0) ! in radiative equilibrium
!             newUdens =  uDensFunc(Teq, rho(i),  5.d0/3.d0, 0.6d0)
!             deltaUdens = abs(newUdens - uDens(i))
!             equilibriumTime =   deltaUDens / abs(adot(i)-etacont(i))
!             if (equilibriumTime < deltaT) then
!                udens(i) = newUdens
!             Be = cSpeed*photonEnergyDensity(i) / fourPi
!                Teq = ((Be * pi)/stefanBoltz)**0.25d0
!                write(*,*) "eq ",equilibriumTime, adot(i),etacont(i), udens(i)
!                udens(i) = udensFunc(Teq,rho(i), 5.d0/3.d0, 0.6d0)
!             else
!                write(*,*) "adding ",oldAdot(i),adot(i)
!                uDens(i) = uDens(i) + (adot(i) - etaCont(i)) * deltaT
!                photonEnergyDensityFromGas(i) = photonEnergyDensityFromGas(i) + (etaCont(i) - adot(i)) * deltaT
!                if (i == 1) write(*,*) "udens ",udens(i), " adot ",adot(1), " eta ",etacont(i), " dt ",deltaT
!                uDens(i) = max (0.d0, uDens(i))
!             endif
             temperature(i) = temperatureFunc(uDens(i), rho(i),  5.d0/3.d0, 0.6d0)
          enddo




          currentTime = currentTime + deltaT

          newDeltaT = 1.d30
          fac = 0.1d0
          do i = 1, nx
             if ((uDens(i) > 0.d0).and.(abs(adot(i) - etaCont(i)) /= 0.d0)) then
                if (adot(i) > 0.d0) then
                   newDeltaT = min(newDeltaT, fac*uDens(i)/abs(adot(i)-etacont(i)))
                   Teq = max(0.d0,(aDot(i) / (4.d0 * stefanBoltz * kappa(i) *rho(i)))**0.25d0) ! in radiative equilibrium
                   newUdens =  uDensFunc(Teq, rho(i),  5.d0/3.d0, 0.6d0)
                   deltaUdens = abs(newUdens - uDens(i))
                   newDeltaT = min(newDeltaT, fac * deltaUDens / abs(adot(i) - etaCont(i)))
                endif
             endif
          enddo
          write(*,*) "thermo time ", newDeltaT

!          if (deltaT > 1.d29) deltaT = diffusionTime
!          deltaT = min(newDeltaT,diffusionTime)

!	  deltaT = diffusionTime 


       probNewPhoton = 0.d0!

       if (dumpNow) then
          iter = iter + 1
          write(*,*) "dumping at ", currentTime, iter
          write(outfile,'(a,i6.6,a)') "udens",iter,".dat"
          write(*,*) iter
          open(21, file=outfile, status="unknown", form="formatted")
          write(21, '(a,1pe12.5)') "# current time ",currentTime
          do i = 1, nx
             write(21, *) xCen(i), uDens(i),  photonEnergyDensityFromSource(i), &
                  photonEnergyDensityFromGas(i), photondensAnalytical(i), udensAnalytical(i)
          enddo
          close(21)
          write(69,*) currentTime, SUM(udens*dx),SUM(photonEnergyDensityFromGas*dx), &
          SUM(udens*dx)+SUM(photonEnergyDensityFromGas*dx), &
               currentNStack, &
               udensAnalytical(1)*xSize,photonDensAnalytical(1)*xSize
          write(*,*) "TOTAL ENERGY: ",SUM(udens*dx)+SUM(photonEnergyDensity*dx)
          write(*,*) "STACK ",SUM( currentEpsStack(1:currentnStack))
          dumpnow = .false.
       endif
          

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

  subroutine addPhotonToStack(useFileForStack, stack,stackUnit, nStack, rVec, uVec, &
       eps, freq, photonFromSource, photonFromGas, beenScattered, radiativeEquPhoton)
    type(STACKTYPE), intent(inout) :: stack
    logical :: useFileForStack
    integer :: nStack, stackUnit
    real(double) :: eps, freq
    logical :: photonFromSource, photonFromGas
    logical :: beenScattered, radiativeEquPhoton
    type(VECTOR) :: rVec, uVec
    
    nStack = nStack + 1

    if (useFileForStack) then
       write(stackUnit) rVec, uVec, eps, freq, photonFromSource, photonFromGas, beenScattered, radiativeEquPhoton
    else
       stack%position(nStack) = rVec
       stack%direction(nStack) = uVec
       stack%eps(nstack) = eps
       stack%freq(nstack) = freq
       stack%photonFromSource(nstack) = photonFromSource
       stack%photonFromGas(nstack) = photonFromGas
       stack%beenScattered(nStack) = beenScattered
       stack%radiativeEquPhoton(nStack) = radiativeEquPhoton
    endif

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

  subroutine getPhotonFromStack(useFileForStack, stack, &
       stackUnit, stacknStack, rVec, uHat, eps, freq, photonFromSource, photonFromGas, beenScattered, &
       radiativeEquPhoton)
    integer :: stackUnit, stacknStack
    logical :: useFileForStack
    type(STACKTYPE) :: stack
    type(VECTOR) :: rVec
    type(VECTOR) :: uHat
    real(double) :: eps, freq
    logical :: photonFromSource, photonFromGas, beenScattered, radiativeEquPhoton
    
    if (stacknStack == 0) then
       write(*,*) "BUG: getPhotonFromStack called with empty stack"
       stop
    endif
    if (useFileForStack) then
       read(stackUnit) rVec, uHat, eps, freq, photonFromSource, photonFromGas, beenScattered,radiativeEquPhoton
    else
       rVec = stack%position(stacknStack)
       uHat = stack%direction(stacknStack)
       eps = stack%eps(stacknStack)
       freq = stack%freq(stacknStack)
       photonFromSource = stack%photonFromSource(stackNStack)
       photonFromGas = stack%photonFromGas(stackNStack)
       beenScattered = stack%beenScattered(stackNStack)
       radiativeEquPhoton = stack%radiativeEquPhoton(stackNStack)
    endif
    
    stacknStack = stacknStack - 1
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

  subroutine calculateProb(prob, bias,  xArray, xCen, nProb, etaCont)
    real(double), intent(inout) :: xArray(:), prob(:), bias(:)
    real(double) :: xCen(:), etaCont(:), totalEmission, dx, tot
    integer :: i, nProb
    dx = xCen(2)-xCen(1)
    prob(1) = 0.d0
    totalEmission = 0.d0
    totalEmission = SUM(etaCont*dx)
    bias(1) = 1.d0
    do i = 2, nProb
       bias(i) = max(1.d-3,min(1.d2,MAXVAL(etacont)/etaCont(i-1)))
    enddo
    bias = 1.d0

    tot = 0.d0
    prob(1) = 0.d0
    do i = 2, nProb
       prob(i) = prob(i-1) + etaCont(i-1) * bias(i) * dx
       tot = tot + etaCont(i-1) * dx
    enddo


    bias = tot * bias / prob(nProb)
 

   prob = prob / prob(nProb)
 
!    open(30,file="prob.dat",status="unknown",form="formatted")
!    do i = 1, nProb
!       write(30,*) xArray(i), prob(i)
!    enddo
!    close(30)
!    stop
  end subroutine calculateProb
       
  subroutine getX(xArray, prob, bias, nProb, x, photonBias)
    real(double) :: xArray(:), prob(:), bias(:), photonBias
    real(double), intent(out) :: x
    integer :: nProb
    real(double) :: r
    integer :: i

    call random_number(r)
    call locate(prob, nprob, r, i)
    x = xArray(i) + (xArray(i+1)-xArray(i))*(r - prob(i))/(prob(i+1)-prob(i))
    photonBias = bias(i) + (bias(i+1)-bias(i))*(r - prob(i))/(prob(i+1)-prob(i))
    photonBias = 1.d0/photonBias
!    photonBias = 1.d0
  end subroutine getX

  subroutine findArrayIndex(xCen, nx, x, iPos)
    real(double) :: xCen(:), x, dx
    integer :: nx
    integer, intent(out) :: iPos

    dx = xCen(2) - xCen(1)
    if (x < xCen(1)-dx/2.d0) then
       write(*,*) "x value outside array ",x,xcen(1)-dx/2.d0
       ipos = 1
       goto 666
    endif
    if (x > xCen(nx)+dx/2.d0) then
       write(*,*) "x value outside array ",x
       ipos = nx
       goto 666
    endif
    if (x > xCen(nx)) then
       iPos = nx
       goto 666
    endif
    if (x < xCen(1)) then
       ipos = 1
       goto 666
    endif
    
!     call locate(xCen, nx, x, iPos)
     iPos = nint( (x - xCen(1)) / dx) + 1
     if (x > xCen(iPos)+dx/2.d0) then
        iPos = iPos + 1
     endif
666 continue
  end subroutine findArrayIndex

  subroutine solveNewUdens(xArray, kappa, rho, photonDensAnalytical, udensAnalytical, nx, k, deltaT, currentTime)
    real(double) :: alpha, kappa(:), rho(:)
    real(double) :: xArray(:), photonDensAnalytical(:), udensAnalytical(:), diffusionTime
    real(double) :: k, deltaT,dx, currentTime
    real(double) :: uPrime(10000), a(10000), b(10000),c(10000)
    integer :: nx, i
    dx = xArray(2)-xArray(1)
    diffusionTime=  0.3d0* (dx**2/(2.d0*k))


    alpha = (k * deltaT)/(dx**2)


    do i = 2, nx - 1
       a(i) = -alpha
       b(i) = 1.d0 + 2.d0 * alpha
       c(i) = -alpha
    enddo
    a(1) = 0.d0
    b(1) = 1.d0
    c(1) = 0.d0
    a(nx) = 0.d0
    b(nx) = 1.d0
    c(nx) = 0.0d0

    call tridiag(a, b, c, photonDensAnalytical, uprime, nx)
    uprime(1) = uprime(2)
    uprime(nx) = uprime(nx-1)

    photonDensAnalytical = uprime ! diffusion terms


!    do i = 1,nx
!       photonDensAnalytical(i) = &
!            1.d10/(sqrt(fourPi*k*currentTime)) * exp(-(xArray(i)-0.5)**2/(4.d0*k*currentTime))
!    enddo

!       do i = 1 , nx ! matter interaction terms
!          t = temperatureFunc(uDensAnalytical(i), rho(i), 5.d0/3.d0, 0.6d0)
!          udensAnalytical(i) = udensAnalytical(i) + &
!               deltaT * (cSpeed*kappa(i)*rho(i)*photonDensAnalytical(i) &
!               - fourPi*kappa(i)*rho(i)*(stefanBoltz/pi)*t**4)
!          photondensAnalytical(i) = photondensAnalytical(i) - &
!               deltaT * (cSpeed*kappa(i)*rho(i)*photonDensAnalytical(i) &
!               - fourPi*kappa(i)*rho(i)*(stefanBoltz/pi)*t**4)
!       enddo



  end subroutine solveNewUdens

  real(double) function temperatureFunc(u, rho,  gamma, mu) result (t)
    real(double) :: u, gamma, rho,mu

    t = (gamma-1.d0)*mu*u / (rGas * rho)

  end function temperatureFunc

  real(double) function uDensFunc(t, rho,  gamma, mu) result (u)
    real(double) :: t, gamma, rho, mu
!    real(double), parameter :: mu = 2.33d0
!    real(double), parameter :: mu = 0.6d0

    u = t * rgas * rho / ((gamma-1.d0)*mu)

  end function uDensFunc

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

  subroutine timeBinPhoton(photonTagTime, wavelength, eps, tau, nSedWavelength, nTime, &
               sedTime, sedWavelength, sedFlux, sedFluxScat, beenScattered)
    real(double) :: photonTagTime
    real(double) :: wavelength
    real(double) :: eps
    real(double) :: tau
    integer :: nSedWavelength, nTime
    real(double) :: sedTime(:), sedWavelength(:)
    real(double) :: sedFlux(nSedWavelength, nTime)
    real(double) :: sedFluxScat(nSedWavelength, nTime)
    logical :: beenScattered
    integer :: iwavelength, iTime

    if ((photonTagTime >= sedTime(1)).and.(photonTagTime <= sedTime(nTime))) then
       call locate(sedWavelength, nSedWavelength, wavelength, iWavelength)
       call locate(sedTime, nTime, photonTagTime, iTime)
       sedFlux(iWavelength, iTime) = sedFlux(iWavelength, iTime) + eps * exp(-tau)/fourPi
       if (beenScattered) then
          sedFluxScat(iWavelength, iTime) = sedFluxScat(iWavelength, iTime) + eps * exp(-tau)/fourPi
       endif
    endif
  end subroutine timeBinPhoton

  recursive subroutine packvalues(thisOctal,nIndex,&
       distanceGridAdot, distanceGridPhotonFromSource, distanceGridPhotonFromGas)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: distanceGridAdot(:)
  real(double) :: distanceGridPhotonFromSource(:)
  real(double) :: distanceGridPhotonFromGas(:)
  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call packvalues(child,nIndex, distanceGridAdot, distanceGridPhotonFromSource, &
                     distanceGridPhotonFromGas)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          distanceGridAdot(nIndex) = thisOctal%distanceGridAdot(subcell)
          distanceGridPhotonFromSource(nIndex) = thisOctal%distanceGridPhotonFromSource(subcell)
          distanceGridPhotonFromGas(nIndex) = thisOctal%distanceGridPhotonFromGas(subcell)

       endif
    enddo
  end subroutine packvalues

  recursive subroutine unpackvalues(thisOctal,nIndex,distanceGridAdot, distanceGridPhotonFromSource, &
       distanceGridPhotonFromGas)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: distanceGridAdot(:)
  real(double) :: distanceGridPhotonFromSource(:)
  real(double) :: distanceGridPhotonFromGas(:)
  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unpackvalues(child,nIndex,distanceGridAdot, &
                     distanceGridPhotonFromSource, distanceGridPhotonFromGas)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          thisOctal%distanceGridAdot(subcell) = distanceGridAdot(nIndex) 
          thisOctal%distanceGridPhotonFromSource(subcell) = distanceGridPhotonFromSource(nIndex)
          thisOctal%distanceGridPhotonFromGas(subcell) = distanceGridPhotonFromGas(nIndex) 
       endif
    enddo
  end subroutine unpackvalues

#ifdef MPI

  subroutine updateGridMPI(grid)
    implicit none
    include 'mpif.h'
    type(gridtype) :: grid
    integer :: nOctals, nVoxels
    real(double), allocatable :: distanceGridAdot(:)
    real(double), allocatable :: distanceGridPhotonFromSource(:)
    real(double), allocatable :: distanceGridPhotonFromGas(:)
    real(double), allocatable :: tempDoubleArray(:)
    integer :: ierr, nIndex

    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    nOctals = 0
    nVoxels = 0
    call countVoxels(grid%octreeRoot,nOctals,nVoxels)
    allocate(distanceGridAdot(1:nVoxels))
    allocate(distanceGridPhotonFromSource(1:nVoxels))
    allocate(distanceGridPhotonFromGas(1:nVoxels))


    nIndex = 0
    call packValues(grid%octreeRoot,nIndex, &
         distanceGridAdot, distanceGridPhotonFromSource, distanceGridPhotonFromGas)

    allocate(tempDoubleArray(nVoxels))


    tempDoubleArray = 0.d0
    call MPI_ALLREDUCE(distanceGridAdot,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    distanceGridAdot = tempDoubleArray 

    tempDoubleArray = 0.d0
    call MPI_ALLREDUCE(distanceGridPhotonFromSource,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    distanceGridPhotonFromSource = tempDoubleArray 

    tempDoubleArray = 0.d0
    call MPI_ALLREDUCE(distanceGridPhotonFromGas,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
    distanceGridPhotonFromGas = tempDoubleArray 

    deallocate(tempDoubleArray)
     
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    
    nIndex = 0
    call unpackValues(grid%octreeRoot,nIndex, &
         distanceGridAdot, distanceGridPhotonFromSource, distanceGridPhotonFromGas)

    deallocate(distanceGridAdot, distanceGridPhotonFromSource,DistanceGridPhotonFromGas)

  end subroutine updateGridMPI
#endif

  recursive subroutine allocateMemoryForTimeDep(thisOctal)
    use octal_mod, only: allocateattribute
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call allocateMemoryForTimeDep(child)
                exit
             end if
          end do
       else
          call allocateAttribute(thisOctal%uDens, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%aDot, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%distancegridaDot, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%distanceGridPhotonFromGas, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%distanceGridPhotonFromSource, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%photonEnergyDensityFromGas, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%photonEnergyDensityFromSource, thisOctal%maxChildren)
          call allocateAttribute(thisOctal%probDistCont, thisOctal%maxChildren)
       endif
    enddo
  end subroutine allocateMemoryForTimeDep


  function distanceToObserver(point, observerPosition, observerDirection) result (distance)
    real(double) :: distance
    type(VECTOR), intent(in) :: point, observerPosition, observerDirection

    distance = (observerDirection .dot. (observerPosition - (1.d10*point)))
  end function distanceToObserver

  subroutine TRIDIAG(A,B,C,R,U,N)
      integer, parameter ::  NMAX=10000
      integer :: n, j
      real(double) ::  GAM(NMAX),A(N),B(N),C(N),R(N),U(N)
      real(double) :: bet
      IF(B(1).EQ.0.)PAUSE
      BET=B(1)
      U(1)=R(1)/BET
      DO J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF(BET.EQ.0.)PAUSE
        U(J)=(R(J)-A(J)*U(J-1))/BET
     enddo
      DO  J=N-1,1,-1
         U(J)=U(J)-GAM(J+1)*U(J+1)
      enddo
    end subroutine TRIDIAG
    
    subroutine quarticSub(udens_n_plus_1, E_n_plus_1, udens_n, E_n, kappa, rho, albedo, &
         mu, gamma, deltaT)
      real(double) :: udens_n_plus_1, E_n_plus_1, udens_n, E_n, kappa, rho, mu, gamma, DeltaT
      real(double) :: a1, a2, albedo
      integer :: nIter
      logical :: converged
      real(double) :: x1, x2, xm, y1, y2, ym
!      a1 = 4.d0 * kappa * rho * stefanBoltz * ( (mu*(gamma-1.d0))/(rGas * rho) )**4 * deltaT
!      a2 = cSpeed * kappa * rho * deltaT

      a1 = 4.d0 * (1.d0-albedo) * kappa * rho * stefanBoltz * ( (mu*(gamma-1.d0))/(rGas * rho) )**4 * deltaT
      a2 = -cSpeed * (1.d0-albedo) * kappa * rho * deltaT * E_n_plus_1

      nIter = 0
      converged = .false.
      
      x1 = 0.d0
      x2 = 1.d15
      
      do while(.not.converged)
         xm = 0.5d0*(x1+x2)
         y1 = quarticFunc(x1, a1, a2, udens_n, E_n)
         y2 = quarticFunc(x2, a1, a2, udens_n, E_n)
         ym = quarticFunc(xm, a1, a2, udens_n, E_n)
                         
         if (y1*ym < 0.d0) then
            x1 = x1
            x2 = xm
         else if (y2*ym < 0.d0) then
            x1 = xm
            x2 = x2
         else
            converged = .true.
            xm = 0.d0
            x1 = 0.d0
            x2 = 0.d0
         endif
                         
         if (abs((x1-x2)/x1) .le. 1.d-5) then
            converged = .true.
         endif
         niter = niter + 1
      enddo
      udens_n_plus_1 = 0.5d0*(x1+x2)



    end subroutine quarticSub

    function quarticFunc(x, a1, a2, udens_n, E_n) result (y)
      real(double) :: x, a1, a2, E_n, y, udens_n
!      y = x**4 + (1.d0+a2)/a1 * x - ((1.d0+a2)*udens_n + a2 * E_n) / a1

      y = x - udens_n + a1 * x**4 + a2
    end function quarticFunc


end module timeDep_mod
