!
!             T   O   R   U  S
!
! torus is a general purpose monte-carlo radiative transfer code used to
! produce polarized spectra of hot and cool stellar winds
!

! written by tjh

! v1.0 on 16/09/99

! raman scattering stuff added 3/3/2000

! OMP parallelization calls added 1/7/2001

! TJH: 19/7/02  Jorick's spotty star stuff added

program torus

  use constants_mod          ! physical constants
  use vector_mod             ! vector math
  use photon_mod             ! photon manipulation
  use grid_mod               ! opacity grid
  use phasematrix_mod        ! phase matrices and stokes vectors
  use math_mod               ! misc maths subroutines
  use blob_mod               ! clumps initialization and movement
  use distortion_mod         ! for distorting opacity grid
  use image_mod              ! stokes image
  use utils_mod
  use stateq_mod
  use inputs_mod
  use TTauri_mod
  use unix_mod
  use path_integral
  use puls_mod
  use input_variables         ! variables filled by inputs subroutine

  implicit none


  integer, parameter :: nOuterLoop = 10
  integer :: iOuterLoop
  integer :: nScat
  integer :: cpuTime, startTime


  ! variables for the grid

  type(GRIDTYPE) :: grid

  integer :: iPhase

  real :: totLineEmission
  real :: totContinuumEmission
  real :: totCoreContinuumEmission
  real :: totCoreContinuumEmission1
  real :: totCoreContinuumEmission2
  real :: totWindContinuumEmission
  real :: probLinePhoton 
  real :: weightContPhoton, weightLinePhoton
  real :: chanceLine, chanceContinuum

  
  ! optical depth variables

  integer, parameter :: maxTau = 100000, maxLambda = 200
  integer :: nTau
  real, allocatable :: contTau(:,:)


!  real :: tauExt(maxTau)
!  real :: tauAbs(maxTau)
!  real :: tauSca(maxTau)
!  real :: lambda(maxTau)
  real, allocatable :: tauExt(:)
  real, allocatable :: tauAbs(:)
  real, allocatable :: tauSca(:)
  real, allocatable :: lambda(:)


  real, allocatable :: kappaAbs(:), kappaSca(:), kappaExt(:)
  real :: dlambda, thisTau

  ! variables to do with dust

  real :: xMin, xMax
  integer, parameter :: nXmie = 20, nMuMie = 20
  type(PHASEMATRIX),allocatable :: miePhase(:, :)
  real :: particleMass, abundance

  ! torus images

  type(IMAGETYPE) :: obsImage, o6image
  type(PVIMAGETYPE), allocatable :: pvimage(:)
  real :: imageSize
  integer :: iSlit
  type(VECTOR) :: slitPosition

  ! intrinsic profile variables

  integer, parameter :: maxIntPro = 1000
  integer :: nIntPro, iSpline
  real :: lamIntPro(maxIntPro), intPro(maxIntPro)
  real, allocatable :: splineArray(:)

  real, allocatable :: lamArray(:)
  real, allocatable :: statArray(:)
  real, allocatable :: sourceSpectrum(:)
  real, allocatable :: sourceSpectrum2(:)

  logical :: resonanceLine


  ! variables for clumped wind models
  

  integer, parameter :: maxBlobs = 10000
  integer :: nCurrent
  type(BLOBTYPE), allocatable :: blobs(:)
  real, parameter :: blobTime = 1000.
  real :: timeEnd = 24.*60.*60.
  real :: timeStart = 0.
  real :: dTime, thisTime

  ! filenames

  character(len=80) :: filename, specFile

  ! photons

  type(PHOTON) :: thisPhoton
  type(PHOTON) :: outPhoton
  type(PHOTON) :: obsPhoton
  type(PHOTON) :: tempPhoton

  ! vectors

  type(VECTOR), parameter :: zAxis = VECTOR(0.,0.,1.)
  type(VECTOR), parameter :: yAxis = VECTOR(0.,1.,0.)
  type(VECTOR), parameter :: xAxis = VECTOR(1.,0.,0.)
  type(VECTOR) :: viewVec, outVec, thisVec, originalViewVec
  type(VECTOR) :: rotationAxis, normToRotation
  type(VECTOR) :: zeroVec, tempVec
  type(VECTOR) :: rHat


  ! output arrays

  integer :: iLambda
  type(STOKESVECTOR), allocatable :: yArray(:)
  type(STOKESVECTOR), allocatable :: errorArray(:,:)
  real, allocatable :: xArray(:)
  real, allocatable :: contWeightArray(:)
  type(STOKESVECTOR) :: tot

  ! model flags

  logical :: escaped, absorbed
  logical :: rotateView
  logical :: tiltView
  logical :: flatSpec
  logical :: ok
  logical :: greyContinuum
  logical :: hitCore
  logical :: firstPlot

  ! model parameters

  real :: vel
  real :: nuStart, nuEnd

  ! single dust blob parameters (WR137 type model)

  real :: meanDustParticleMass, getMeanMass


  real :: foreground = 0., background = 0.  ! plotting intensities
  real :: mu
  real :: t1, t2, t3
  real :: junk1


  ! raman scattering model parameters

  logical :: thruStar
  type(VECTOR) :: hotSourcePosition, coolStarPosition
  type(VECTOR) :: ramanSourceVelocity

  ! O VI spectrum stuff

  character(len=80) :: o6filename
  integer, parameter :: no6pts = 100
  real, parameter :: o6start = 1031.5, o6end=1032.8
  real :: o6xarray(no6pts), o6yarray(no6pts)


  ! misc

  real :: rotateDirection
  real :: meanr_line = 0., meanr_cont = 0.
  real :: wtot_line =0., wtot_cont = 0.
  real :: meanr0_line = 0., meanr0_cont = 0.
  real :: wtot0_line =0., wtot0_cont = 0.
  real :: junk
  character(len=80) :: plotFile

  real :: thisChi, thisSca, albedo
  logical :: normalizeSpectrum
  integer :: currentScat
  logical :: redRegion
  real :: r
  logical :: contWindPhoton
  real :: ramanWeight
  real :: thisVel
  real :: fac1, fac2, fac3
  real :: vRay, vOverCsqr
  real :: directionalWeight
  real :: weight
  real :: r1, r2
  real :: t
  integer :: i, j, k, n
  integer :: i1,i2,i3
  integer :: nTot  
  real :: observedLambda
  real :: thisLam
  real :: nu
  real :: fac
  real :: escProb
  logical :: contPhoton
  integer :: nContPhotons
  real :: phi
  real :: deltaLambda
  real :: rStar
  integer :: findilambda

  real :: Laccretion, Taccretion, fAccretion, sAccretion
  real :: theta1, theta2, chanceHotRing

  ! Spot stuff
  
  real :: chanceSpot                     ! chance of spot
  logical :: spotPhoton                  ! photon from spot?

  ! binary parameters

  type(VECTOR) :: starPos1, starPos2
  integer :: nVec
  type(VECTOR), allocatable :: distortionVec(:)


  ! initialize

  useNdf = .true.
  movie = .false.
  thinLine = .false.
  rotateView = .false.
  flatspec = .false.
  greyContinuum = .false.
  secondSource = .false.
  firstplot = .true.
  doRaman = .false.
  enhance = .false.
  rotateDirection = -1.

  contFluxFile = "none"
  intProFilename = "none"

  inputKappaSca = 0.
  inputKappaAbs = 0.

  zeroVec = VECTOR(0.,0.,0.)

  starPos1 = VECTOR(-400.,0.,0.)
  starPos2 = VECTOR(+400.,0.,0.)


  ! hardwired stuff

  abundance = 1.e-8
  particleMass = 15.9994*mHydrogen + 47.867*mHydrogen

  do i = 1, no6pts
     o6xArray(i) = o6start + (o6end-o6start)*real(i-1)/real(no6pts-1)
     o6yarray(i) = 1.d-10
  enddo

  ! get the model parameters

  call  inputs() ! variables are passed using the input_variables module


!  secondSourcePosition = zeroVec

  if (nLambda > maxLambda) then
     write(*,'(a,i3.3,a)') "nlambda is greater than the maximum value (",&
          maxLambda,")"
     stop
  endif

  if (geometry.eq."raman") then
     hotSourcePosition = 0.75*secondSourcePosition
     coolStarPosition = (-0.25)*secondSourcePosition
     secondSourcePosition = hotSourcePosition
     write(*,*) "hot",hotSourcePosition
  endif

  if (doRaman) screened = .true.

  if (.not.inputOK) goto 666

  if (dopvimage) then
     allocate(pvimage(1:nSlit))
  endif

  if (mie) then
     meanDustParticleMass = getMeanMass(amin, amax, qdist)
     write(*,*) "mean dust particle mass: ",meanDustParticleMass
  endif

  probLinePhoton = 1. - probContPhoton
  inclination = inclination*degToRad
  scale = scale * rSol

  rStar = 100.*rSol

  ! if we are computing spiral models then we rotate the view

  if ((nPhase /= 1).and.(distortionType .eq. "spiral")) rotateview = .true.

  ! we also rotate for the test model

  if (distortionType .eq. "test") rotateView = .true.

  ! we rotate wind-wind collision models


  if (distortionType .eq. "windwind") rotateView = .true.


  ! rotate for colliding winds too

  if (trim(geometry) .eq. "collide") rotateView = .true.

  if (geometry == "binary") rotateView = .true.

  ! rotate for raman

  if (geometry == "raman") rotateView = .true.


  if (geometry == "planet") rotateView = .true.

  if (geometry == "betacep") rotateView = .true.

  if (geometry == "donati") rotateView = .true.

  if ((geometry == "ttauri").and.(.not.enhance)) rotateView = .true.

  if (doRaman) then
        rotateView = .true.
        rotateDirection = 1.
  endif

  if (geometry == "rolf") rotateView = .true.

  if (geometry == "disk") rotateView = .true.

  ! switches for line emission

  if (lineEmission) then
     flatspec = .true.
     greyContinuum = .true.
  endif

  if (geometry == "planet") then
     flatspec = .true.
     greyContinuum = .true.
  endif

  if (geometry == "rolf") then
!     flatspec = .true.
     greyContinuum = .true.
  endif
     


  ! the observer's viewing direction

  rotationAxis = VECTOR(0., 0., 1.)

  normToRotation = rotationAxis .cross. yAxis
  call normalize(normToRotation)

  originalViewVec%x = 0.
  originalViewVec%y = -sin(inclination)
  originalViewVec%z = -cos(inclination)
  viewVec = originalViewvec


  ! allocate the grid - this might crash out through memory problems

  if (gridcoords /= "polar") then
     if ((.not.doRaman).and.(geometry /= "binary"))  then
        grid = initCartesianGrid(nx, ny, nz, nLambda, lamStart, lamEnd, &
             greyContinuum, ok)
        else
        grid = initCartesianGrid(nx, ny, nz, 1, lamStart, lamEnd, &
             greyContinuum, ok)
     endif
     if (.not.ok) goto 666
  else   if (plezModelOn) then
     grid = plezModel("model.plez", lamStart, lamEnd, nLambda, kfac)
  else
     grid = initPolarGrid(nr, nmu, nphi, nlambda, lamStart, &
                          lamEnd, greyContinuum, ok)
     if (.not.ok) goto 666
  endif


  ! a few dynamic arrays

  allocate(miePhase(1:nLambda,1:nMumie))
  allocate(kappaExt(1:nLambda))
  allocate(kappaAbs(1:nLambda))
  allocate(kappaSca(1:nLambda))
  allocate(lamArray(1:nLambda))
  allocate(statArray(1:nLambda))
  allocate(sourceSpectrum(1:nLambda))
  allocate(sourceSpectrum2(1:nLambda))

  grid%doRaman = doRaman

  statArray = 0.
  sourceSpectrum = 1.
  sourceSpectrum2 = 1.

  ! set up the wavelength array

  do i = 1, nLambda
     lamArray(i) = lamStart+(lamEnd-lamStart)*real(i-1)/real(nLambda-1)
  enddo

  if (plezModelOn) then
     call readTioCrossSection(lamArray, nLambda, kappaExt, kappaAbs, kappaSca)
  endif


  ! allocate the output arrays

  allocate(xArray(1:nLambda))
  allocate(yArray(1:nLambda))
  allocate(errorArray(1:nOuterLoop,1:nLambda))

  deltaLambda = (lamEnd - lamStart) / real(nLambda)

  xArray(1) = lamStart + deltaLambda/2.
  yArray(1)%i = 0.
  yArray(1)%q = 0.
  yArray(1)%u = 0.
  yArray(1)%v = 0.
  do i = 2, nLambda
     xArray(i) = xArray(i-1) + deltaLambda
     yArray(i)%i = 0.
     yArray(i)%q = 0.
     yArray(i)%u = 0.
     yArray(i)%v = 0.
  enddo


  errorArray(1:nOuterLoop,1:nLambda) = STOKESVECTOR(0.,0.,0.,0.)

  if (doRaman) then
     yarray(1:nLambda)%i = 1.e-20
     errorArray(1:nOuterLoop,1:nLambda)%i = 1.e-20
  endif



  ! fill up the grid with the appropriate opacities

  if (.not.plezModelOn) then
     select case(geometry)
     case("torus")
        call fillGridTorus(grid, rho, rTorus, rOuter)
     case("sphere")
        call fillGridSpheriod(grid, rho, radius, kFac)
     case("flared")
        call fillGridFlaredDisk(grid, meanDustParticleMass)
     case("ellipse")
        call fillGridEllipse(Grid,rho, rMin, rMaj)
     case("disk")
        call fillGridDisk(grid, rho, rCore, rInner, rOuter, height, mCore, diskTemp)
     case("star")
        call fillGridStar(grid, radius, mdot, vel, kfac, scale)
     case("spiral")
        call fillGridSpiral(grid, radius, mdot, vel, scale)
     case("shell")
        call fillGridShell(grid, radius, shellFrac, rho, kfac)
     case("stateq")
        call fillGridStateq(grid, opacityDataFile, kfac, scaleDensity)
     case("bipolar")
        call fillGridBipolar(grid, rho,  30.)
     case("collide")
        call fillGridCollide(grid, rho, momRatio, binarySep, mie, meanDustParticleMass, logMassLossRate)
     case("dustblob")
        call fillGridDustBlob(grid, dustBlobdistance, rho, phiDustBlob, &
             xDustBlobSize, yDustBlobSize, zDustBlobSize)
     case("raman")
        call  fillGridRaman(grid, 10.*rStar, mdot, vterm, rStar, coolStarPosition, beta)
     case("binary")
        call fillGridBinary(grid, radius1, radius2, mass1, mass2, &
        temp1, temp2, vNought1, vNought2, vTerm1, vTerm2, beta1, beta2, period, &
        mdot1, mdot2, deflectionAngle, shockWidth, shockFac)
	
        call initGridStateq(grid, contFluxFile, contFluxFile2, &
	     popFilename, readPops, writePops, lte, nLower, nUpper)
	
     case("rolf")
        call fillGridRolf(grid, mdot, vterm, o6width)
     case("wr137")
        call fillGridWR137(grid, rCore, mDot, vTerm, beta, temp1)
     case("planet")
        call fillGridPlanet(grid)
     case("hourglass")
        call fillGridHourglass(grid)
     case("ttauri")
           call fillGridMagneticAccretion(grid,contfluxfile, popFileName, &
                readPops, writePops, lte,  lamLine, Laccretion, Taccretion, sAccretion, &
                curtains, dipoleOffset, nLower, nUpper, theta1, theta2)
     case("ttwind")
        call fillGridTTauriWind(grid,contfluxfile, popFileName, &
             readPops, writePops, lte, nLower, nUpper)
     case("betacep")
        call fillGridBetaCep(grid)
     case("donati")
        resonanceLine = .false.
        if (lamLine < 2000.) resonanceLine = .true.
        call fillGridDonati2(grid, resonanceLine, misc)
        if (.not.resonanceLine) then
           call initGridStateq(grid, contFluxFile, contFluxFile2, popFilename, &
                readPops, writePops, lte, nLower, nUpper)
        endif
        grid%kappaAbs = 1.e-30

     case("puls")
        call fillGridPuls(grid, mDot, rcore, tEff, v0, vterm, beta)
!        call initGridStateq(grid, contFluxFile, contFluxFile2, popFilename, &
!                            readPops, writePops, lte, nLower, nUpper)
     case("wind")
	call fillGridWind(grid, mDot, rStar, tEff, v0, vterm, beta, &
       lte, contFluxFile, writePops, readPops, popFilename, nLower, nUpper)
	call initgridstateq(grid, contfluxFile, " ", popFileName, &
	     readPops, writePops, lte, nLower, nUpper)
	grid%etaCont = 1.e-30
	grid%kappaAbs = 1.e-30


	
     case("resonance")
        call fillGridResonance(grid, rCore, mDot, vTerm, beta, temp1)

     case DEFAULT
        write(*,*) "! Unrecognised grid geometry: ",trim(geometry)
        goto 666
     end select
  endif

  
  if (mie) then
     call fillGridMie(grid, scale, aMin, aMax, qDist, grainType)
  endif

  if (fillTio) then
     call fillGridTio(grid, scale)
  endif

  if (fillRayleighOpacity) then
     call fillGridRayleigh(grid,scale)
  endif

  if (fillThomson) then
     call fillGridThomson(grid)
  endif

  ! the distortion types

  select case(distortionType)
  case("spiral")
     call distortGridSpiral(grid, vRot, nSpiral)
  case("rotation")
     call distortRotation(grid, vRot)
     write(*,'(a,f5.1,a)') "Grid distorted by a rotational velocity of ", &
                            vRot/1.e5," km/s"
  case("test")
     call distortGridTest(grid)
!  case("raman")
!     call distortRaman(grid)
  case("raman")
     call  distortStrom(grid, secondSourcePosition, .true., .true., 0.5*rStar, coolStarPosition, ramanDist)

     case("wrdisk")
        call distortWRdisk(grid)
     case("windwind")
        call distortWindCollision(grid, momRatio, binarySep)

  end select
  
  if (enhance) then
     nVec = 40
     allocate(distortionVec(nVec))
     timeStart = 0.
     timeEnd = 12.*3600.
     dTime = (timeEnd - timeStart)/real(nPhase)
     call initInfallEnhancement(distortionVec, nVec)
     do i = 1, nStartPhase-1
	call infallEnhancment(grid, distortionVec, nVec, dTime, .false.)
     enddo
     
  endif

  ! The source spectrum is normally a black body

  if (.not.grid%lineEmission) then
     do i = 1, nLambda
        nu = cSpeed / (xArray(i)*angstromToCm)
        fac= log(2.)+log(hConst) + 3.*log(nu)
        sourceSpectrum(i) = exp(fac)/(exp(hConst*nu/(kConst*Teff))-1.)
     enddo
  endif



  ! set up a random seed

  call random_seed()


  outVec = (-1.)* originalViewVec


  ! compute the mie scattering phase matrices if necessary

  xMin =  2.*pi*grainSize/(lamEnd*angstromToCm)
  xMax =  2.*pi*grainSize/(lamStart*angstromToCm)

  if (mie) then
     write(*,'(a)') "Computing Mie phase grid..."
     do i = 1, nLambda
        do j = 1, nMumie
           mu = 2.*real(j-1)/real(nMumie-1)-1.
           thislam = lamStart+(lamEnd-lamStart)*real(i-1)/real(nLambda-1)
           call mieDistPhaseMatrix(aMin, aMax, qDist, thislam, &
                mu, miePhase(i,j), grainType)

        enddo
     enddo
     write(*,'(a)') "Completed."
  endif


  ! initialize the blobs if required

  if (nBlobs > 0) then

     allocate(blobs(1:maxBlobs))
     if (freshBlobs) then
        do i = 1 , maxBlobs
           blobs(i)%inUse = .false.
        enddo
        write(*,'(a)') "Running blobs for five days..."
 
        ! now we run a few of days worth of blobs

        t1 = 0.
        t2 = 120. * 60. * 60.
        dTime = (t2 - t1) / 100.
        do i = 1, 100
           call addNewBlobs(grid, maxBlobs, blobs, blobTime, dTime, &
                            nCurrent, blobContrast)
           call moveBlobs(maxBlobs, blobs, 0., dTime, grid)
           write(*,*) i,nCurrent
        enddo

        j = 0 
        do i = 1, maxBlobs
           if (blobs(i)%inUse) j = j + 1
        enddo

        write(*,'(a,i3)') "Number of blobs after 5 days: ",j


        dTime = 0.
        if (nPhase /= 1) then
           dTime = (timeEnd - timeStart) / real(nPhase-1)
        endif

        open(50,file="files.lis",status="unknown",form="formatted")

        write(*,'(a)') "Running blobs and writing configuration files"
        do i = 1, nPhase

           thisTime = timeStart + (timeEnd-timeStart)*real(i-1)/real(nPhase-1)
           write(specFile,'(a,i3.3,a)') trim(outfile),i,".dat"
           write(50,*) specFile(1:30),thisTime

           call addNewBlobs(grid, maxBlobs, blobs, blobTime, dTime, nCurrent, blobContrast)
           call moveBlobs(maxBlobs, blobs, 0., dTime, grid)
           write(filename,"(a,i3.3,a)") "run",i,".blob"
           call writeBlobs(filename, maxBlobs, blobs)

        enddo

        close(50)

     endif
  endif


  if (plotVelocity) then
     call plotVelocityVectors(grid, device)
  endif




  if (nPhase /= 1) then
     dTime = (timeEnd - timeStart) / real(nPhase-1)
  endif

  if (doRaman) then
     ramanSourceVelocity = (ramVel/cSpeed)*VECTOR(0.,-1., 0.)
     write(*,'(a,f7.1,a)') "Raman source speed: ",modulus(ramanSourceVelocity)*cSpeed/1.e5, " km/s"
  endif

  if (geometry == "ttauri") then
     rotationAxis = grid%diskNormal
     write(*,*) "rotation axis",rotationAxis
  endif

  if (geometry == "donati") then
     rotationAxis = rotateX(rotationAxis, -dipoleOffset)
  endif

  if (distortionType == "spiral") then
     rotationAxis = rotateX(rotationAxis, dipoleOffset)
  endif


  normToRotation = rotationAxis .cross. yAxis
  call normalize(normToRotation)

  tempVec  = arbitraryRotate(rotationAxis, inclination, normToRotation)
  originalViewVec =  (-1.)*tempVec


  if (modulus(thermalElectronVelocity(10000.)) == 0.) then
     write(*,'(a)') "THERMAL ELECTRON BROADENING IS OFF!!!!!!!!!!!!!!!!!!!!!!!!!!"
  endif

! From here we do multiple runs if required !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  tiltView = .false.

  do iPhase = nStartPhase, nEndPhase
          

     call unixTimes(cpuTime, startTime)

     viewVec = originalViewVec
     outVec = (-1.)*viewVec
     thisVec = viewVec

     ! we rotate the view by an appropriate amount


     if (rotateView.and.(nPhase /= 1)) then
        write(*,'(a)') "Rotating view..."
        phi = -rotateDirection * twoPi * real(iPhase-1)/real(nPhase-1)
        viewVec =  arbitraryRotate(thisVec, phi, rotationAxis)
        outVec = (-1.) * viewVec
       	write(*,'(a,f5.2,a,f5.2,a,f5.2,a)') "View vector: (",viewVec%x,",", &
	     viewVec%y,",",viewVec%z,")"
     endif

     if (phaseOffset /= 0.) then
        write(*,'(a,f5.2)') "Rotating by phase offset: ",phaseOffset
        phi = twoPi * phaseOffset
        viewVec =  arbitraryRotate(viewVec, phi, rotationAxis)
        outVec = (-1.) * viewVec
     endif


     if (tiltView.and.(nPhase /= 1)) then
        phi = -twoPi * real(iPhase-1)/real(nPhase-1)
        viewVec =  arbitraryRotate(thisVec, phi, xAxis)
        outVec = (-1.) * viewVec
     endif



     ! zero the output arrays

     xArray(1) = lamStart + deltaLambda/2.
     yArray(1)%i = 0.
     yArray(1)%q = 0.
     yArray(1)%u = 0.
     yArray(1)%v = 0.
     do i = 2, nLambda
        xArray(i) = xArray(i-1) + deltaLambda
        yArray(i)%i = 0.
        yArray(i)%q = 0.
        yArray(i)%u = 0.
        yArray(i)%v = 0.
     enddo

     if (doRaman) then
        yArray(1:nLambda)%i = 1.e-20
     endif

     if (stokesimage) then
        if (grid%cartesian) then
           imageSize = max(abs(grid%xAxis(1)),abs(grid%yAxis(1)),abs(grid%zAxis(1)))
           obsImage = initImage(50, imageSize, vmin, vmax)
           if (doRaman) then
              o6image = initImage(50, imageSize, vmin, vmax)
           endif
        else
           select case (geometry)
              case("disk")
                 obsImage = initImage(50, grid%rAxis(grid%nr), vMin, vMax)
              case("flared")
                 obsImage = initImage(50, 4.*grid%rAxis(1), vMin, vMax)
              case DEFAULT
                 obsImage = initImage(50, min(5.*grid%rAxis(1),grid%rAxis(grid%nr)), vMin, vMax)
           end select
        endif
     endif


     if (dopvImage) then

        if (nSlit  > 1) then
           rHat = slitPosition2 - slitPosition1
           do iSlit = 1, nSlit
              slitPosition = slitPosition1 + (real(iSlit-1)/real(nSlit-1))*(slitPosition2 - slitPosition1)
              pvImage(iSlit) = initPVimage(nv, vMin, vMax, np, -slitLength/2., slitLength/2., &
                   slitPosition, slitPA, slitWidth, slitLength)
           enddo
        else
           slitPosition = slitPosition1
           pvImage(1) = initPVimage(200, vMin, vMax, 200, -slitLength/2., slitLength/2., &
                slitPosition, slitPA, slitWidth, slitLength)
        endif
           
     endif


     ! refill the grids 
     

     if (.not.plezModelOn) then
        select case(geometry)
        case("torus")
!           call fillGridTorus(grid, rho, rTorus, rOuter)
        case("sphere")
           call fillGridSpheriod(grid, rho, radius, kFac)
        case("ellipse")
           call fillGridEllipse(Grid,rho,  rMin, rMaj)
        case("flared")
!           call fillGridFlaredDisk(grid)
        case("disk")
!           call fillGridDisk(grid, rho, rCore, rInner, rOuter, height, mCore, diskTemp)
        case("star")
           call fillGridStar(grid, radius, mdot, vel, kfac, scale)
        case("spiral")
           call fillGridSpiral(grid, radius, mdot, vel,  scale)
        case("shell")
           call fillGridShell(grid, radius, shellFrac, rho, kfac)
        case("stateq")
           call fillGridStateq(grid, opacityDataFile, kfac, scaleDensity)
        case("bipolar")
           call fillGridBipolar(grid, rho, 30.)
        case("collide")
           call fillGridCollide(grid, rho, momRatio, binarySep, mie, meanDustParticleMass, logMassLossRate)
        case("dustblob")
           call fillGridDustBlob(grid, dustBlobdistance, rho, phiDustBlob, &
                                 xDustBlobSize, yDustBlobSize, zDustBlobSize)
        case("raman")
           call  fillGridRaman(grid, 10.*rStar, mdot, vterm, rStar, coolStarPosition, beta)
        case("binary")
!           call fillGridBinary(grid, opacityDataFile, opacityDataFile2, momRatio, mass1, mass2, period, shockwidth, shockFac)
        case("rolf")
!           call fillGridRolf(grid, mdot, vterm)
        case("wr137")
           call fillGridWR137(grid, rCore, mDot, vTerm, beta, temp1)
        case("planet")
!           call fillGridPlanet(grid)
        case("hourglass")
!           call fillGridHourglass(grid)

        case("ttauri")
           if (enhance) then
	      call fillGridMagneticAccretion(grid,contfluxfile, popFileName, &
		   readPops, writePops, lte,  lamLine, Laccretion, Taccretion, sAccretion, &
		   curtains, dipoleOffset, nLower, nUpper, theta1, theta2)
	      call infallEnhancment(grid, distortionVec, nVec, dTime, .true.)
	   endif
        !        call fillGridMagneticAccretion(grid)
     case("ttwind")
        
     case("betacep")
        !        call fillGridBetaCep(grid)
     case("donati")
!        call fillGridDonati(grid, resonanceLine)
     case("puls")
	!        call fillGridPuls(grid, mDot, rStar, tEff, v0, vterm, beta)
	case("wind")
	   
!	call fillGridWind(grid, mDot, rStar, tEff, v0, vterm, beta, &
!       lte, contFluxFile, writePops, readPops, popFilename, nLower, nUpper)

     case("resonance")
!        call fillGridResonance(grid, rCore, mDot, vTerm, beta, temp)


	
     case DEFAULT
        write(*,*) "! Unrecognised grid geometry: ",trim(geometry)
        goto 666
     end select
  endif


!     if (mie) then
!        call fillGridMie(grid, scale, aMin, aMax, qDist, grainType)
!     endif

     if (fillTio) then
        call fillGridTio(grid, scale)
     endif

     if (fillRayleighOpacity) then
        call fillGridRayleigh(grid,scale)
     endif

     if (fillThomson) then
        call fillGridThomson(grid)
     endif


     select case(distortionType)
     case("spiral")
!        call distortGridSpiral(grid, vRot, nSpiral)
     case("rotation")
!        call distortRotation(grid, vRot)
!        write(*,'(a,f5.1,a)') "Grid distorted by a rotational velocity of ",vRot/1.e5," km/s"
     case("test")
        call distortGridTest(grid)
!     case("raman")
!        call distortRaman(grid)
      case("raman")
         call  distortStrom(grid, secondSourcePosition, .true., .true., 0.5*rStar, coolStarPosition, ramanDist)
     case("wrdisk")
        call distortWRdisk(grid)

     case("windwind")
        call distortWindCollision(grid, momRatio, binarySep)

     end select

     ! we read in the blob configuration from a file

     if (nBlobs > 0) then
        write(filename,"(a,i3.3,a)") "run",iPhase,".blob"
        call readBlobs(filename, maxBlobs, blobs)
        call distortGridWithBlobs(grid, maxBlobs, blobs)
     endif

     ! if we are producing a movie then plot this phase

     if (movie) then
        if (device(2:3) == "xs") then
           call plotGrid(grid, viewVec,  opaqueCore, &
                device,contTau, foreground, background, coolStarPosition, firstPlot)
        endif  
        write(filename,"(a,i3.3,a)") "frame",iPhase,".gif/gif"
        call plotGrid(grid, viewVec,  opaqueCore, &
             filename, contTau, foreground, background, coolStarPosition, firstPlot)
     endif

     if (geometry .eq. "rolf") then
        secondSourcePosition = VECTOR(grid%xAxis(grid%nx/2), &
             grid%yAxis(grid%ny/2), &
             grid%zAxis(grid%nz/2))
     endif
       


     ! chuck out some useful information to the user

     write(*,*) " "
     write(*,'(a)') "Some basic model parameters"
     write(*,'(a)') "---------------------------"
     write(*,*) " "

     allocate(lambda(1:maxTau))
     allocate(tauExt(1:maxTau))
     allocate(tauAbs(1:maxTau))
     allocate(tauSca(1:maxTau))
     allocate(contTau(1:maxTau,1:maxLambda))
     allocate(contWeightArray(1:maxTau))
     


     write(*,'(a,f7.1,a)') "Cross-sections at ",lamStart, " angstroms"
     write(*,'(a)') "------------------------------------------"
     write(*,*) " "

     call integratePath(lamArray(1),  lamLine, VECTOR(1.,1.,1.), zeroVec, &
                       VECTOR(1.,0.,0.), grid, lambda, tauExt, tauAbs, &
                       tauSca, maxTau, nTau, opaqueCore, escProb, .false. , &
                       lamStart, lamEnd, nLambda, contTau, hitCore, thinLine, .false., rStar, &
                       coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)

     write(*,'(a,1pe10.3)') "Optical depth in x-axis from centre: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth in x-axis from centre: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth in x-axis from centre: ",tauSca(ntau)
     write(*,*) " "
     
     call integratePath(lamArray(1), lamLine, VECTOR(1.,1.,1.), zeroVec, &
                       VECTOR(0.,1.,0.), grid, lambda, tauExt, tauAbs, &
                       tauSca, maxTau, nTau, opaqueCore, escProb, .false. , &
                       lamStart, lamEnd, nLambda, contTau, hitCore, thinLine, .false., rStar, &
                       coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)


     write(*,'(a,1pe10.3)') "Optical depth in y-axis from centre: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth in y-axis from centre: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth in y-axis from centre: ",tauSca(ntau)
     write(*,*) " "
     
     call integratePath(lamArray(1),  lamLine, VECTOR(1.,1.,1.), zeroVec, &
                       VECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
                       tauSca, maxTau, nTau, opaqueCore, escProb, .false. , &
                       lamStart, lamEnd, nLambda, contTau, hitCore, thinLine, .false., rStar, &
                       coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)


     write(*,'(a,1pe10.3)') "Optical depth in z-axis from centre: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth in z-axis from centre: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth in z-axis from centre: ",tauSca(ntau)
     write(*,*) " "



     call integratePath(lamLine,  lamLine, VECTOR(1.,1.,1.), &
          zeroVec, outVec, grid, lambda, &
          tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
          .false., lamStart, lamEnd, nLambda, contTau, &
          hitCore, thinLine,.false., rStar, coolStarPosition, 1., &
          .false., 0, 0, 0., 0., 0., junk)

     write(*,'(a,1pe10.3)') "Optical depth to observer: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth to observer: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth to observer: ",tauSca(ntau)
     write(*,*) " "


     weightLinePhoton = 0.
     weightContPhoton = 1.



     if (geometry == "hourglass") then
        call computeProbDist(grid, totLineEmission, &
             totWindContinuumEmission,lamline, .false.)
        weightLinePhoton = 1.
        weightContPhoton = 0.
        probLinePhoton = 1.
        probContPhoton = 0.
     endif

     if (doRaman) then
        call computeProbDist(grid, totLineEmission, &
             totWindContinuumEmission,lamline, .false.)
        write(*,*) "Total Raman Line Emission: ",totLineEmission
        weightLinePhoton = 1.
        weightContPhoton = 0.
        probLinePhoton = 1.
        probContPhoton = 0.
     endif


     if (grid%lineEmission) then

        ! integrate the line and continuum emission and read the 
        ! intrinsic profile

        if (coreEmissionLine) then
           totLineEmission = 0.
           totWindContinuumEmission = 0.
        else
           if (.not.grid%resonanceLine) then
              call computeProbDist(grid, totLineEmission, &
                   totWindContinuumEmission,lamline, useBias)

              ! convert from per steradian
              
              totLineEmission = totLineEmission * fourPi
              totWindContinuumEmission = totwindContinuumEmission * fourPi
           endif
        endif

        do i = 1, grid%nMu
           write(*,*) i,acos(grid%muAxis(i))*radtodeg, &
                             grid%muProbDistLine(grid%nr/2,i)
        enddo

        if (geometry == "donati") totWindContinuumEmission = 0.

        if (geometry == "binary") totWindContinuumEmission = 0.

        if (grid%resonanceLine) totWindContinuumEmission = 0.

        if (geometry == "ttauri") then
           totWindContinuumEmission = 0.
           write(*,'(a)') "! Wind continuum emission switched off."
        endif

        if (lineOff) then
           totLineEmission = 0.
        endif

        nu = cSpeed / (lamLine * angstromtocm)
        write(*,*) "Line emission: ",totLineEmission


        select case(geometry)
           case("binary")
              call contread(contFluxFile, nu, totCoreContinuumEmission1)
              call contread(contFluxFile2, nu, totCoreContinuumEmission2)
           case("puls")
              totCoreContinuumEmission = pi * blackBody(0.77 * tEff, lamLine)
           case DEFAULT
              call contread(contFluxFile, nu, totCoreContinuumEmission)
        end select


        ! reads the intrinsic core absorption profile for core-halo models


        if (trim(intProFilename) /= "none") then
           call rdIntPro(intProFilename, intPro, lamIntPro, nIntPro)
           allocate(splineArray(1:nIntPro))
           call spline(lamIntPro, intPro, nIntPro, 0., 0., splineArray)
           do iSpline = 1 , nLambda
              if ((xArray(iSpline) < lamIntPro(1)) .or.  &
                   (xArray(iSpline) > lamIntPro(nIntPro))) then
                 sourceSpectrum(iSpline) = 1.
              else
                 call splint(lamIntPro, intPro, splineArray, nIntPro, &
                      xArray(iSpline), sourceSpectrum(iSpline))
              endif
           enddo
           deallocate(splineArray)
        endif

        if (geometry == "binary") then
           call rdIntPro(intProFilename2, intPro, lamIntPro, nIntPro)
           allocate(splineArray(1:nIntPro))
           call spline(lamIntPro, intPro, nIntPro, 0., 0., splineArray)
           do iSpline = 1 , nLambda
              if ((xArray(iSpline) < lamIntPro(1)) .or.  &
                   (xArray(iSpline) > lamIntPro(nIntPro))) then
                 sourceSpectrum2(iSpline) = 1.
              else
                 call splint(lamIntPro, intPro, splineArray, nIntPro, &
                      xArray(iSpline), sourceSpectrum2(iSpline))
              endif
           enddo
           deallocate(splineArray)
        endif

        ! the core can emit as a Gaussian emission line

        if (coreEmissionLine) then
           call computeCoreEmissionProfile(xArray, sourceSpectrum, nLambda, &
                lamLine, velWidthCoreEmissionLine, relIntCoreEmissionLine)
        endif


        nuStart = cSpeed / (lamStart * angstromtocm)
        nuEnd = cSpeed / (lamEnd * angstromtocm)

        totWindContinuumEmission = totWindContinuumEmission * (nuStart - nuEnd)

        write(*,*) "Wind cont emission: ",totWindContinuumEmission


! factor of four pi taken out from below - input continuum files are expected
! to be in fluxes (erg/s/cm^2/hz) not Hnu's

        if (geometry /= "binary") then
           totCoreContinuumEmission = (totCoreContinuumEmission * &
                (nuStart - nuEnd))*fourPi*grid%rCore**2
        else
           totCoreContinuumEmission1 = totCoreContinuumEmission1 * (nuStart-nuEnd) * fourPi * grid%rStar1**2
           totCoreContinuumEmission2 = totCoreContinuumEmission2 * (nuStart-nuEnd) * fourPi * grid%rStar2**2
           totCoreContinuumEmission = totCoreContinuumEmission1 + totCoreContinuumEmission2
           grid%lumRatio = totCoreContinuumEmission1 / totCoreContinuumEmission
           write(*,*) "Binary luminosity ratio (p/s): ",totCoreContinuumEmission1/totCoreContinuumEmission2
        endif
           
        if (geometry == "ttauri") then
              fAccretion =  blackBody(Taccretion, lamLine)*(nuStart-nuEnd)*sAccretion
              write(*,*) "laccretion/core",fAccretion/totCoreContinuumEmission
              totCoreContinuumEmission = totCoreContinuumEmission + &
                   fAccretion
	      chanceHotRing = fAccretion/totCoreContinuumEmission
        endif
        
        chanceSpot = 0.
        if ((geometry == "disk").and.(nSpot > 0)) then
           chanceSpot = fSpot * blackBody(tSpot, 6562.8) / &
                ((1.-fSpot)*blackBody(tEff, 6562.8) + &
                (fSpot * blackBody(tSpot, 6562.8)))
           write(*,'(a,f5.3)') "Spot chance at 6563A: ",chanceSpot
        endif




        write(*,*) "Core Continuum Emission: ",totCorecontinuumEmission


        totContinuumEmission = totCoreContinuumEmission + totWindContinuumEmission

        write(*,*) "Continuum emission: ", totContinuumEmission

        if ((totContinuumEmission + totLineEmission) /= 0.) then
           chanceLine = totLineEmission/(totContinuumEmission + totLineEmission)
           chanceContinuum = totContinuumEmission / &
                (totContinuumEmission + totLineEmission)
           grid%chanceWindOverTotalContinuum = totWindContinuumEmission &
                / max(1.e-30,totContinuumEmission)
        else
           chanceLine =0.
           chanceContinuum = 1.
           grid%chanceWindOverTotalContinuum = 0.
        endif



        write(*,*) "Chance continuum emission in wind: ", grid%chanceWindOverTotalContinuum

        ! set up the line and continuum weights

        if ((probLinePhoton /= 0.).and.(probContphoton /= 0.)) then

           weightLinePhoton = chanceLine / probLinePhoton
           weightContPhoton = chanceContinuum / probContPhoton / real(nLambda)
        else
           if (probLinePhoton == 0.) then
              weightLinePhoton = 0.
              weightContPhoton = 1.
           else
              weightLinePhoton = 1.
              weightContPhoton = 0.
           endif
        endif




        write(*,*) "Line photon weight: ", weightLinePhoton
        write(*,*) "Continuum photon weight: ", weightContPhoton

        write(*,*) " "
        write(*,*) "Line photon prob: ", probLinePhoton
        write(*,*) "Continuum photon prob: ", probContPhoton


        write(*,*) "chance line",chanceline


        ! compute the probability distributions


     endif

     deallocate(lambda)
     deallocate(tauSca)
     deallocate(tauExt)
     deallocate(tauAbs)
     deallocate(contTau)
     deallocate(contWeightArray)




     nContPhotons = nint(probContPhoton * real(nPhotons) / real(nOuterLoop))




     write(*,*) " "
     write(*,'(a)') "Run-time messages"
     write(*,'(a)') "-----------------"
     write(*,*) " "



     tot%i = 0.
     tot%q = 0.
     tot%u = 0.
     tot%v = 0.


     k = nPhotons/10

     ntot = 0

     ! now we loop 10 times over a tenth of the photons - this will be used
     ! to help compute the errors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call unixTimes(cpuTime, startTime)


     if (grid%geometry == "donati") then
        junk1 = 0.
        rHat = VECTOR(0.,0.,0.)
        n = 0
        do i = 1, 100000
           call initPhoton(thisPhoton, grid, nLambda, xArray, sourceSpectrum, &
                lineEmission, lamLine, weightLinePhoton, &
                weightContPhoton, contPhoton, flatspec, vRot, pencilBeam, &
                secondSource, secondSourcePosition, lumRatio, &
                ramanSourceVelocity, vo6, contWindPhoton, directionalWeight, useBias, theta1, theta2, &
		chanceHotRing, &
                nSpot, chanceSpot, thetaSpot, phiSpot, fSpot, spotPhoton)
	  
           if (thisPhoton%linePhoton) then
              junk1 = junk1 + thisPhoton%position%z
              rHat = rHat + thisPhoton%velocity
              n = n + 1
           endif
           call getIndices(grid, thisPhoton%position, i1, i2, i3, t1, t2, t3)
        enddo
        write(*,*) "Average z position",junk1/real(n)
        write(*,*) "Average velocity",((cSpeed/1.e5)*(1./real(n)))*rHat
     endif



     do iOuterLoop = 1, nOuterLoop

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i, contPhoton, contWindPhoton, r, nScat) &
!$OMP PRIVATE(thisPhoton, directionalWeight) &
!$OMP PRIVATE( ilambda, sourceSpectrum, ok) &
!$OMP PRIVATE(hitCore, junk, thisLam, j, weight, thisVel) &
!$OMP PRIVATE(i1, i2, i3, t1, t2, t3, vray, vovercsqr, fac, observedLambda) &
!$OMP PRIVATE(t, rHat, islit, fac1, fac2, fac3, obsPhoton, r1, r2, thisTau) &
!$OMP PRIVATE(escaped, currentScat, absorbed, dlambda, thisChi, thisSca) &
!$OMP PRIVATE(albedo, tempPhoton, redRegion, thrustar, ramanWeight) &
!$OMP PRIVATE(outPhoton) &
!$OMP PRIVATE(nTau, escProb) &
!$OMP PRIVATE(lambda, tauExt, tauSca, tauAbs, contTau, contWeightArray) &

!$OMP SHARED(grid) &

!$OMP SHARED(meanr_Cont, wtot_cont,meanr_line,wtot_line, ntot) &
!$OMP SHARED(nContPhotons, nPhotons, lineEmission, lamLine, nLambda) &
!$OMP SHARED(weightLinePhoton, flatSpec, vRot, secondSource, secondSourcePosition) &
!$OMP SHARED(lumRatio, ramanSourceVelocity, vO6, doRaman, Xarray) &
!$OMP SHARED(weightContPhoton, useBias, pencilBeam ,outVec)&
!$OMP SHARED(opaqueCore, lamStart, lamEnd, thinLine, rStar, coolStarPosition) &
!$OMP SHARED(viewVec, o6xArray,o6yArray, rotationAxis, o6image, screened) &
!$OMP SHARED(xArray, yArray, statArray, stokesImage, obsImage, doPvimage) &
!$OMP SHARED(nSlit, pvimage, gridDistance, meanr0_line, wtot0_line) &
!$OMP SHARED(sourceSpectrum2, meanr0_cont,wtot0_cont, maxScat, mie) &
!$OMP SHARED(miePhase, zeroVec) &

        do i = 1, nPhotons/nOuterLoop

           allocate(lambda(1:maxTau))
           allocate(tauExt(1:maxTau))
           allocate(tauAbs(1:maxTau))
           allocate(tauSca(1:maxTau))
           allocate(contTau(1:maxTau,1:maxLambda))
           allocate(contWeightArray(1:maxTau))



           ! continuum or line photon?

           if (i <= nContPhotons) then
              contPhoton = .true.
           else
              contPhoton = .false.
           endif

           contWindPhoton = .false.

           ! there is a chance that this continuum photon is produced in the wind rather
           ! than at the core...

           if (lineEmission .and. contPhoton) then
              call random_number(r)
              if (r < grid%chanceWindOverTotalContinuum) then
                 contWindPhoton = .true.
              endif
           endif

           nScat = 0

           ! initialize the photon
           
           contWeightArray(1:maxLambda) = 1.

           select case(grid%geometry)
              case("planet")
                 call initPlanetPhoton(thisPhoton, grid, lamLine)
                 directionalWeight = 1.
                 contWindPhoton = .true.
              case DEFAULT
                 call initPhoton(thisPhoton, grid, nLambda, xArray, sourceSpectrum, &
                      lineEmission, lamLine, weightLinePhoton, &
                      weightContPhoton, contPhoton, flatspec, vRot, pencilBeam, &
                      secondSource, secondSourcePosition, lumRatio, &
                      ramanSourceVelocity, vo6, contWindPhoton, directionalWeight, useBias, &
		      theta1, theta2, chanceHotRing,  &
                      nSpot, chanceSpot, thetaSpot, phiSpot, fSpot, spotPhoton)

                 if (thisPhoton%resonanceLine) then
                    r1 = real(i)/real(nPhotons/nOuterLoop)
                    thisPhoton%lambda = xArray(1) + r1*(xArray(nLambda)-xArray(1))
                 endif

           end select


           if (thisPhoton%contPhoton) then

              meanr_cont = meanr_cont + modulus(thisPhoton%position)*thisPhoton%stokes%i
              wtot_cont = wtot_cont + thisPhoton%stokes%i
           else
              meanr_line = meanr_line + modulus(thisPhoton%position)*thisPhoton%stokes%i
              wtot_line = wtot_line + thisPhoton%stokes%i
           endif


           iLambda = findIlambda(thisPhoton%lambda, xArray, nLambda, ok)

           ! now we fire the photon direct to the observer

           if (doRaman) then
              
              call integratePath(thisPhoton%lambda, lamLine, &
                   thisPhoton%velocity, &
                   thisPhoton%position, outVec, grid, &
                   lambda, tauExt, tauAbs, tauSca, maxTau , nTau, opaqueCore, &
                   escProb, thisPhoton%contPhoton, lamStart, lamEnd, &
                   nLambda, contTau, hitCore, thinLine,.false., rStar,&
                   coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)
              thisLam = thisPhoton%lambda + (thisPhoton%velocity .dot. viewVec) * 1031.928
              j = findIlambda(thisLam, o6xArray, no6pts, ok)
              if (ok) then
                 weight = oneOnFourPi * exp(-tauExt(nTau))
                 o6yArray(j) = o6yArray(j) + weight
                 thisVel = (thisLam-lamLine)/lamLine
                 call addPhotonToImage(viewVec, rotationAxis,o6Image, thisPhoton, thisVel, weight)

              endif
           endif

           if (.not. screened) then


              ! find optical depths to observer



              call integratePath(thisPhoton%lambda, lamLine, &
                   thisPhoton%velocity, &
                   thisPhoton%position, outVec, grid, &
                   lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, &
                   escProb, thisPhoton%contPhoton, lamStart, lamEnd, &
                   nLambda, contTau, hitCore, thinLine,.false., rStar,&
                   coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)


              call getIndices(grid,thisPhoton%position,i1,i2,i3,t1,t2,t3)



              if (thisPhoton%resonanceLine) escProb = 1.

              if (thisPhoton%linePhoton) then
                 vray = -(thisPhoton%velocity .dot. outVec)
                 vovercsqr = modulus(thisPhoton%velocity)**2
                 fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                 observedlambda = thisPhoton%lambda / fac
                 weight = oneOnFourPi * exp(-tauExt(nTau)) * escProb
              endif

              if (thisPhoton%contPhoton.and.(.not.contWindPhoton)) then
                 t = modulus(thisPhoton%position)
                 if (t /= 0.) then
                    rHat = thisPhoton%position / t
                 else
                    rHat = thisPhoton%direction
                 endif
                 t = rHat .dot. outVec

                 if (lineEmission) then
                    weight = abs(t)*exp(-tauExt(nTau)) / pi
                 else
                    weight = oneOnFourPi*abs(t)*exp(-tauExt(nTau))
                 endif
                 observedlambda = thisPhoton%lambda
              endif


              if (.not.flatSpec.and.(.not.hitCore)) then
                 weight = oneOnFourPi*exp(-tauExt(nTau))
                 iLambda = findIlambda(observedlambda, xArray, nLambda, ok)

                 if (ok) then
                    yArray(iLambda) = yArray(iLambda) + &
                         (thisPhoton%stokes * weight)
                    statArray(iLambda) = statArray(iLambda) + 1.
                 endif
                 if (stokesImage) then
                    thisVel = 0.
                    call addPhotonToImage(viewVec, rotationAxis, obsImage, thisPhoton, thisVel, weight)
                 endif
                 if (doPVimage) then
                    do iSlit = 1, nSlit
                       call addPhotontoPVimage(pvImage(iSlit), thisPhoton, viewVec, rotationAxis, thisVel, &
                                 weight, gridDistance)
                    enddo
                 endif
              endif
              
              if (flatSpec.and.(.not.hitCore)) then
                 if (thisPhoton%linePhoton) then
                    fac2 = 1.
!                    if (thisPHoton%resonanceline) fac2 = directionalWeight
                    iLambda = findIlambda(observedlambda, xArray, nLambda, ok)
                    if (ok) then
                       yArray(iLambda) = yArray(iLambda) + &
                            (thisPhoton%stokes * (fac2 * oneOnFourPi * escProb * exp(-tauExt(nTau))))
                       statArray(iLambda) = statArray(iLambda) + 1.
                    endif

                    weight = oneOnFourPi * escProb * exp(-tauExt(nTau))*fac2

                    meanr0_line = meanr0_line + modulus(thisPhoton%position) * &
                         (thisPhoton%stokes%i * weight)
                    wtot0_line = wtot0_line + (thisPhoton%stokes%i * weight)
                    

                    thisVel = (observedLambda-lamLine)/lamLine
                    if (stokesImage) then
                       call addPhotonToImage(viewVec, rotationAxis, obsImage, thisPhoton, thisVel, weight)
                    endif
                    if (dopvImage) then
                       do iSlit = 1, nSlit
                          call addPhotontoPVimage(pvImage(iSlit), thisPhoton, viewVec,  rotationAxis,thisVel, &
                            weight, gridDistance)
                       enddo
                    endif

                 else
                    i1 = 0
                    do iLambda = 1, nLambda
                       fac1 = 1.
                       if (.not.contWindPhoton) then
                          fac1 = 2.*abs(thisPhoton%originalNormal.dot.outVec)/twoPi
                       else
                          fac1 = oneOnfourPi
                       endif

                       vray = -(thisPhoton%velocity .dot. outVec)
                       vovercsqr = modulus(thisPhoton%velocity)**2
                       fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                       observedlambda = thisPhoton%lambda / fac


                       thisLam = (lamLine-observedlambda) + xArray(iLambda)

                       call hunt(xArray,nLambda,thisLam,i1)
                       if (i1 ==0 ) i1 = 1
                       if (grid%geometry /= "binary") then
                          fac2 = sourceSpectrum(i1)

                          if ((grid%geometry == "disk").and.(.not.spotPhoton).and.(.not.photLine)) fac2 = 1.


                       else
                          if (thisPhoton%fromStar1) then
                             fac2 = sourceSpectrum(i1)
                          else
                             fac2 = sourceSpectrum2(i1)
                          endif
                       endif
                       fac3 = contTau(nTau,i1)
                       if (thinLine) fac3 = 0.
                       weight = (fac1 * exp(-(tauExt(ntau)+fac3)))*fac2
                       
                       yArray(iLambda) = yArray(iLambda) + &
                                           (thisPhoton%stokes * weight)
                       statArray(iLambda) = statArray(iLambda) + 1.
                       meanr0_cont = meanr0_cont + modulus(thisPhoton%position) * thisPhoton%stokes%i*weight
                       wtot0_cont = wtot0_cont + thisPhoton%stokes%i*weight
                       thisVel = (observedLambda-lamLine)/lamLine
                       if (stokesImage) then
                          call addPhotonToImage(viewVec,  rotationAxis,obsImage, thisPhoton, thisVel, weight)
                       endif
                       if (dopvImage) then
                          do iSlit = 1, nSlit
                             call addPhotontoPVimage(pvImage(islit), obsPhoton, viewVec,  rotationAxis, thisVel, &
                               weight, gridDistance)
                          enddo
                       endif

                    enddo
                 endif
              endif

           endif


           call integratePath(thisPhoton%lambda, lamLine, &
                thisPhoton%velocity, &
                thisPhoton%position, &
                thisPhoton%direction, grid, &
                lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, &
                escProb, thisPhoton%contPhoton, lamStart, lamEnd, &
                nLambda, contTau, hitCore, thinLine,.false., rStar, &
                coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)



           call random_number(r1)
           fac = (1.d0-exp(-tauExt(nTau)))
           thisTau = min(-log(1.d0-r1*fac),dble(tauExt(nTau)))

           if (thisPhoton%contPhoton) escProb = 1.

           if (thisPhoton%resonanceLine) escProb = 1.


           thisPhoton%stokes = thisPhoton%stokes * (fac * escProb)


           currentScat = 0

           escaped = .false.


           absorbed = .false.
           if (doRaman) hitCore = .false.


! scattering loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

           do while (.not.escaped .and. .not.absorbed .and. &
                     (nScat < maxScat)   )


              nScat = nScat + 1

              currentScat = currentScat  + 1

              ! find position of next interaction

              call locate(tauExt, nTau, thisTau, j)


              t = 0.
              if ((tauExt(j+1) - tauExt(j)) /= 0.) then
                 t = (thisTau - tauExt(j)) / (tauExt(j+1) - tauExt(j))
              endif
              dlambda = lambda(j) + (lambda(j+1)-lambda(j))*t


              thisPhoton%position = thisPhoton%position + dlambda*thisPhoton%direction

              contWeightArray(1:nLambda) = contWeightArray(1:nLambda) * &
                   exp(-(contTau(j,1:nLambda) +t*(contTau(j+1,1:nLambda)-contTau(j,1:nLambda))))


              call getIndices(grid, thisPhoton%position, i1, i2, i3, t1, t2, t3)
              thisChi = interpGridKappaAbs(grid, i1, i2, i3, 1, t1, t2, t3)
              thisSca = interpGridKappaSca(grid, i1, i2, i3, 1, t1, t2, t3)
              
              if ((thisChi+thisSca) > 0.) then
                 albedo = thisSca / (thisChi + thisSca)
              else
                 albedo = 0.
              endif

              if (grid%resonanceLine) albedo = 1.

              if (hitCore) then
                 if (grid%geometry /= "binary") then
                    r = modulus(thisPhoton%position)/grid%rCore
                    if (r < 1.01) then
                       absorbed = .true.
                       albedo = 0.
                    endif
                 else
                    r1 = modulus(thisPhoton%position - grid%starpos1)/grid%rStar1
                    r2 = modulus(thisPhoton%position - grid%starpos2)/grid%rStar2
                    r = min (r1, r2)
                    if (r < 1.01) then
                       albedo = 0.5
                    endif
                 endif
              endif
                 

              thisPhoton%stokes = thisPhoton%stokes * albedo
                   

              if (thisPhoton%stokes%i < reallySmall) then
                 absorbed = .true.
!                 write(*,*) "! Small photon weight",thisPhoton%stokes%i
              endif


              ! towards observer

              if (.not.absorbed) then

                 tempPhoton = thisPhoton
                 


                 call scatterPhoton(grid, tempPhoton, outVec, obsPhoton, mie, &
                       miePhase, nLambda, nMumie, lamStart, lamEnd)


                 ! the o6 photon might get scattered towards the observer by a rayleigh scattering

                 if (doRaman) then
                    call integratePath(obsPhoton%lambda, lamLine, obsPhoton%velocity, &
                         obsPhoton%position, obsPhoton%direction, grid, &
                         lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, &
                         escProb, obsPhoton%contPhoton, lamStart, lamEnd, nLambda, contTau, hitCore, &
                         thinLine, .false., rStar, coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)
                    vray = -(obsPhoton%velocity .dot. outVec)
                    vovercsqr = modulus(thisPhoton%velocity)**2
                    fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                    observedlambda = obsPhoton%lambda / fac
                    
                    weight = oneOnFourPi*exp(-tauExt(nTau)) * 34./(6.6+34.)

                    j = findIlambda(observedLambda, o6xArray, no6pts, ok)
                    if (ok) then
                       o6yArray(j) = o6yArray(j) + weight
                    endif
                 endif




                 redRegion = .false.

                 if (doRaman) redRegion = .true.
                 
                 call integratePath(obsPhoton%lambda, lamLine, &
                      obsPhoton%velocity, &
                      obsPhoton%position, obsPhoton%direction, grid, &
                      lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, &
                      escProb, obsPhoton%contPhoton, lamStart, lamEnd, &
                      nLambda, contTau, hitCore, &
                      thinLine, redRegion, rStar, coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)




                 if (obsPhoton%linePhoton) then

                    
                    vray = -(obsPhoton%velocity .dot. outVec)
                    vovercsqr = modulus(obsPhoton%velocity)**2
                    fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                    if (.not.doRaman) then
                       observedlambda = obsPhoton%lambda / fac
                    else
                       observedlambda = obsPhoton%redlambda / fac
                    endif


                 else
                    observedLambda = obsPhoton%lambda
                 endif

                 thruStar = .false.

                 if (.not.flatSpec.and.(.not.hitCore).and.(.not.thruStar)) then

                    ramanWeight = 1.
                    if (doRaman) ramanWeight = 6.6 / (34.+6.6)

                    weight = oneOnFourPi*exp(-tauExt(nTau)) * ramanWeight

                    iLambda = findIlambda(observedlambda, xArray, nLambda, ok)
                      

                    if (ok) then
                       yArray(iLambda) = yArray(iLambda) + obsPhoton%stokes*weight
                       statArray(iLambda) = statArray(iLambda) + 1.
                    endif
                    thisVel = (1./(1./observedLambda  + 1./1215.67))/1031.928-1.

                    if (stokesImage) then
                       call addPhotonToImage(viewVec,  rotationAxis, obsImage, obsPhoton, thisVel, weight)
                    endif
                    if (dopvImage) then
                       do iSlit = 1, nSlit
                          call addPhotontoPVimage(pvImage(iSlit), obsPhoton, viewVec,  rotationAxis,thisVel, &
                            weight, gridDistance)
                       enddo
                    endif

                    if (doRaman) then
                       thisPhoton%stokes = thisPhoton%stokes*(1.-ramanWeight)
                    endif
                 endif



                 if (flatSpec.and.(.not.hitCore)) then



                    if (obsPhoton%linePhoton) then
                       iLambda = findIlambda(observedlambda, xArray, nLambda, ok)
                       weight = oneOnFourPi*exp(-tauExt(nTau))


                       if (obsPhoton%resonanceLine) weight = weight * directionalWeight

                       if (ok) then
                          yArray(iLambda) = yArray(iLambda) + obsPhoton%stokes*weight
                          statArray(iLambda) = statArray(iLambda) + 1.
                       endif

                       thisVel = observedLambda
                       thisVel = (observedLambda-lamLine)/lamLine
                       if (stokesImage) then
                          call addPhotonToImage(viewVec,  rotationAxis, obsImage, obsPhoton, thisVel, weight)
                       endif
                       if (dopvImage) then
                          do iSlit = 1 , nSlit
                             call addPhotontoPVimage(pvImage(iSlit), obsPhoton, viewVec,  rotationAxis, thisVel, &
                               weight, gridDistance)
                          enddo
                       endif


                    else


                       vray = -(obsPhoton%velocity .dot. outVec)
                       vovercsqr = modulus(obsPhoton%velocity)**2
                       fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                       observedlambda = obsPhoton%lambda / fac

                       do iLambda = 1,nLambda

                          thisLam = (lamLine-observedLambda) + xArray(iLambda)

                          i1 = findILambda(thisLam, xArray, nLambda, ok)
                          if (grid%geometry /= "binary") then
                             fac2 = sourceSpectrum(i1)
                             if ((grid%geometry == "disk").and.(.not.spotPhoton).and.(.not.photLine)) fac2 = 1.
                          else
                             if (obsPhoton%fromStar1) then
                                fac2 = sourceSpectrum(i1)
                             else
                                fac2 = sourceSpectrum2(i1)
                             endif
                          endif

                          fac3 = contTau(nTau, i1)



                          if (thinLine) fac3 = 0.


                          weight = exp(-(tauExt(nTau)+fac3)) * oneOnFourpi & 
                          * fac2 * directionalWeight * contWeightArray(i1)

                          

                          if (ok) then
                             yArray(iLambda) = yArray(iLambda) + &
                                  (obsPhoton%stokes*weight)
                          else 
                             yarray(ilambda) = yArray(ilambda) + &
                                  obsPhoton%stokes*(oneOnFourpi* directionalweight * exp(-tauExt(ntau)))
                          endif


                          thisVel = (observedlambda-lamLine)/lamLine
                          if (stokesImage) then
                             call addPhotonToImage(viewVec,  rotationAxis, obsImage, obsPhoton, thisVel, weight)
                          endif
                          if (dopvImage) then
                             do iSlit = 1, nSlit
                                call addPhotontoPVimage(pvImage(iSlit), obsPhoton, viewVec,  rotationAxis, thisVel, &
                                  weight, gridDistance)
                             enddo
                          endif


                          statArray(iLambda) = statArray(iLambda) + 1.
                       enddo
                    endif
                 endif


                 call scatterPhoton(grid,thisPhoton, zeroVec, outPhoton, mie, &
                       miePhase, nLambda, nMuMie, lamStart, lamEnd)
                 thisPhoton = outPhoton


                 call integratePath(thisPhoton%lambda, lamLine, thisPhoton%velocity, &
                      thisPhoton%position, &
                      thisPhoton%direction, grid, &
                      lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, &
                      escProb, thisPhoton%contPhoton, lamStart, lamEnd, nLambda, contTau, hitCore, &
                      thinLine, .false.,rStar, coolStarPosition, 1., .false., 0, 0, 0., 0., 0., junk)




                 if (maxScat == 1) then
                    call random_number(r1)
                    thistau = -log(max(1.e-10,(1. - r1)))
                    if (thistau .gt. tauExt(nTau)) then
                       escaped = .true.
                    else
                       nscat = 0
                    endif
                 else
                    call random_number(r1)
                    fac = 1.d0 - exp(-tauExt(nTau))
                    thisTau = -log(1.d0-r1*fac)
                    thisPhoton%stokes = thisPhoton%stokes * fac
                 endif


              endif


           enddo
           nTot = nTot + nScat

           deallocate(lambda)
           deallocate(tauSca)
           deallocate(tauExt)
           deallocate(tauAbs)
           deallocate(contTau)
           deallocate(contWeightArray)
        enddo

!$OMP END PARALLEL DO

        write(*,'(i8,a,f7.3)') iOuterLoop*nPhotons/nOuterLoop," photons done"

        errorArray(iOuterLoop,1:nLambda) = yArray(1:nLambda)

     enddo

! fudge !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     yArray(1) = yArray(2)
     yArray(nLambda) = yArray(nLambda-3)
     yArray(nLambda-1) = yArray(nLambda-3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     write(*,*) " "
     write(*,'(a)') "Model summary"
     write(*,'(a)') "-------------"
     write(*,*) " "

     call systemInfo(startTime,nPhotons)


     if (.not.grid%cartesian.and.(grid%rCore /= 0.)) then
        if (wtot_line /= 0.) write(*,*) "Mean radius of line formation",meanr_line/wtot_line/grid%rCore
        if (wtot0_line /= 0.) write(*,*) "Mean radius of line zero",meanr0_line/wtot0_line/grid%rCore
        if (wtot_cont /= 0.) write(*,*) "Mean radius of cont formation",meanr_cont/wtot_cont/grid%rCore
        if (wtot0_cont /=0.) write(*,*) "Mean radius of cont zero",meanr0_cont/wtot0_cont/grid%rCore
     endif
     normalizeSpectrum = .false.
     if (.not.doRaman) then
        normalizeSpectrum = .true.
     endif

     if (mie) normalizeSpectrum = .false.


     if (nPhase == 1) then
        call  writeSpectrum(outFile,  nLambda, xArray, yArray,  errorArray, nOuterLoop, normalizeSpectrum, useNdf)
     else
        write(specFile,'(a,i3.3)') trim(outfile),iPhase
        call  writeSpectrum(specFile,  nLambda, xArray, yArray,  errorArray, nOuterLoop, normalizeSpectrum, useNdf)
        if (doRaman) then
           write(o6filename,'(a,a,i3.3,a)') trim(outfile),"_o6_",iPhase,".dat"
           open(20,file=o6filename,status="unknown",form="formatted")
           do i = 1, no6pts
              t1 = 1./(o6xArray(2)-o6xArray(1))
              write(20,*) o6xarray(i),o6yarray(i)*t1,1.e-20,1.e-30,1.e-20,1.e-30
           enddo
           close(20)
        endif
     endif


     if (stokesimage) then
        write(specFile,'(a,a,i3.3)') trim(outfile),"_image",iPhase
        call writeImage(obsImage, specfile)
        if (doRaman) then
        write(specFile,'(a,a,i3.3)') trim(outfile),"_o6image",iPhase
           call writeImage(o6Image, specfile)
        endif
     endif

     if (doPvimage) then
        do iSlit = 1, nSlit
           write(specFile,'(a,a,i3.3,a,i2.2)') trim(outfile),"_pvimage",iPhase,"_slit_",iSlit

           call smoothPVimage(pvImage(iSlit), vfwhm/2.35, pfwhm/2.35)

           call writePVimage(pvImage(iSlit), specfile, vSys)
           if (device /= "/xs") then
              write(plotfile,'(a,i3.3,a,i2.2,a,a)') "slitplot",iPhase,"_slit_",iSlit,".gif",trim(device)
           else
              plotfile = device
           endif
           call plotSlitOnImage(obsImage, PVimage(iSlit), plotfile, gridDistance)
        enddo
     endif

     if (stokesImage) call freeImage(obsImage)

     if (doPVimage) then
        do iSlit = 1, nSlit
           call freePVimage(pvImage(iSlit))
        enddo
     endif

  enddo




666 continue

!  call freeGrid(grid)

end program torus





subroutine writeSpectrum(outFile,  nLambda, xArray, yArray,  errorArray, nOuterLoop, &
     normalizeSpectrum, useNdf)

  use phasematrix_mod

  implicit none
  integer :: nLambda
  character(len=*) :: outFile
  real :: xArray(nLambda)
  logical :: useNdf
  integer :: nOuterLoop
  logical :: normalizeSpectrum
  type(STOKESVECTOR) :: yArray(nLambda), errorArray(nOuterloop,nLambda)
  type(STOKESVECTOR),pointer :: ytmpArray(:)
  real, allocatable :: meanQ(:), meanU(:), sigQ(:), sigU(:)
  real, allocatable :: stokes_i(:), stokes_q(:), stokes_qv(:)
  real, allocatable :: stokes_u(:), stokes_uv(:)
  real :: x
  integer :: i

  allocate(ytmpArray(1:nLambda))

  allocate(meanQ(1:nLambda))
  allocate(meanU(1:nLambda))
  allocate(sigQ(1:nLambda))
  allocate(sigU(1:nLambda))

  allocate(stokes_i(1:nLambda))
  allocate(stokes_q(1:nLambda))
  allocate(stokes_qv(1:nLambda))
  allocate(stokes_u(1:nLambda))
  allocate(stokes_uv(1:nLambda))

!  x = SUM(yArray(1:min(10,nLambda))%i)/real(min(10,nLambda))

!  x = 1./x

!  x = 

  if (normalizeSpectrum) then
     if (yArray(1)%i /= 0.) then
        x = 1.d0/yArray(1)%i
     else
        x  = 1.d0
     endif
  else
     x = 1.d0/(xArray(2)-xArray(1))
  endif


  do i = 1, nLambda
     ytmpArray(i) = yArray(i) * x
  enddo

  where(errorArray%i /= 0.)
     errorArray%q = errorArray%q / errorArray%i
     errorArray%u = errorArray%u / errorArray%i
  end where

  do i = 1, nLambda
     meanQ(i) = sum(errorArray(1:nOuterLoop,i)%q) / real(nOuterLoop)
     meanU(i) = sum(errorArray(1:nOuterLoop,i)%u) / real(nOuterLoop)
  enddo

  do i = 1, nLambda
     sigQ(i) = sqrt(sum((errorArray(1:nOuterLoop,i)%q-meanQ(i))**2)/real(nOuterLoop-1))
     sigU(i) = sqrt(sum((errorArray(1:nOuterLoop,i)%u-meanU(i))**2)/real(nOuterLoop-1))
  enddo

  stokes_i = ytmpArray%i
  stokes_q = ytmpArray%q
  stokes_u = ytmpArray%u
  stokes_qv = (ytmpArray%i * sigQ)**2
  stokes_uv = (ytmpArray%i * sigU)**2

  if (useNdf) then
     call wrtsp(nLambda,stokes_i,stokes_q,stokes_qv,stokes_u, &
          stokes_uv,xArray,outFile)
  else
     open(20,file=trim(outFile)//".dat",status="unknown",form="formatted")
     do i = 1, nLambda
        write(20,*) xArray(i),stokes_i(i), stokes_q(i), stokes_qv(i), &
             stokes_u(i), stokes_uv(i)
     enddo
     close(20)
  endif

end subroutine writeSpectrum


integer function findIlambda(lambda, xArray, nLambda, ok)
  use utils_mod
  implicit none
  integer :: nlambda, i
  real :: lambda
  real :: xArray(nLambda), dx
  logical :: ok

  ok = .true.

  dx = xArray(2) - xArray(1)
  if (lambda < (xArray(1))) then
     findiLambda = 1
     ok = .false.
     goto 666
  endif
  if (lambda > (xArray(nLambda))) then
     findiLambda = nLambda
     ok = .false.
     goto 666
  endif

  call locate(xArray, nLambda, lambda, i)
  if (lambda > (xArray(i)+dx/2.)) then
     i= i+1
  endif
  findilambda = min(i,nLambda)
666 continue
end function findilambda



