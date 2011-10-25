
module phaseloop_mod

  public :: do_phaseloop

  private :: choose_view

CONTAINS

subroutine do_phaseloop(grid, alreadyDoneInfall, meanDustParticleMass, rstar, vel,  &
     coolstarposition, Laccretion, Taccretion, fAccretion, sAccretion, corecontinuumflux, &
     starsurface, sigmaAbs0, sigmaSca0, ttauri_disc, distortionVec, nvec,       &
     infallParticleMass, maxBlobs, flatspec, maxTau, &
     miePhase, nsource, source, blobs, nmumie, dTime,overrideFilename, overrideInclinations)

  use kind_mod
  use inputs_mod 
  use phasematrix_mod
  use disc_class
  use image_mod, only: IMAGETYPE, PVIMAGETYPE, initImage, initPVImage, addPhotonToImage, addphotontopvimage, &
       createlucyimage, freeimage, freepvimage, smoothpvimage, writefalsecolourppm
#ifdef USECFITSIO
  use image_mod, only : writeFitsImage
#endif
  use source_mod, only: SOURCETYPE
  use photon_mod, only: PHOTON, initPhoton, scatterPhoton, initplanetphoton
  use math_mod, only: thermalElectronVelocity, interpGridScalar2
  use romanova_class, only: romanova
  use surface_mod, only: SURFACETYPE, buildSphere, createsurface, createttaurisurface2, createttaurisurface, &
       emptysurface, testsurface
  use filter_set_class, only: filter_set, get_nfilter, get_filter_name, get_filter_set_name, FWHM_filters, &
       lambda_eff_filters, info_filter_set, make_filter_set
  use grid_mod, only: fillgridbipolar, fillgridcollide, fillgriddustblob, fillgridellipse, fillgridraman, &
       fillgridshell,fillgridspheriod, fillgridspiral, fillgridstar, fillgridstateq, fillgridwr137, getIndices 
  use amr_mod, only: tauAlongPath2, findsubcelllocal, findsubcelltd, amrupdategrid, countVoxels, amrGridValues, tauAlongPathFast, returnKappa
  use path_integral, only: integratePath, test_optical_depth
  use stateq_mod, only: amrStateq
  use math_mod, only: interpGridKappaAbs, interpGridKappaSca, computecoreemissionprofile, computeprobdist 
  use mpi_global_mod, only: myRankGlobal
#ifdef MPI
  use mpi
  use mpi_global_mod, only: nThreadsGlobal
  use random_mod
  use grid_mod, only: freeGrid
  use gridio_mod, only: readamrgrid
#endif
  use TTauri_mod, only: fillGridMagneticAccretion, infallenhancment
  use blob_mod, only: blobtype, distortgridwithblobs, readblobs
  use lucy_mod, only: calccontinuumemissivitylucy, calccontinuumemissivitylucymono, setbiasontau
  use timing, only: tune
  use formal_solutions, only: compute_obs_line_flux
  use distortion_mod, only: distortgridtest, distortstrom, distortwindcollision, distortwrdisk
  use gridtype_mod, only: GRIDTYPE       
  use gridio_mod, only: writeamrgrid
  use parallel_mod, only: torus_mpi_barrier
  use utils_mod, only: locate, hunt, findIlambda, blackBody, spline, splint
  use dust_mod, only: createDustCrossSectionPhaseMatrix, stripDustAway
  use source_mod, only: sumSourceLuminosityMonochromatic, sumSourceLuminosity, randomSource
  use random_mod
  use sed_mod, only: getNumSedInc, getSedInc
  use physics_mod, only : setupdust
  implicit none

! Arguments
  type(GRIDTYPE) :: grid
  character(len=*), optional :: overrideFilename
  real, optional, intent(in):: overrideInclinations(:)
  logical :: alreadyDoneInfall
  real :: meanDustParticleMass
  real :: rstar
  real :: vel
  type(VECTOR) :: coolStarPosition
  real(double) :: Laccretion, finalTau
  real :: Taccretion, fAccretion, sAccretion
  real(double) :: corecontinuumflux
  type(SURFACETYPE) :: starSurface
  type(romanova) :: romData ! parameters and data for romanova geometry
  real  :: sigmaAbs0, sigmaSca0  ! cross section at the line centre
  type(alpha_disc)  :: ttauri_disc       ! parameters for ttauri disc
  integer :: nvec
  type(VECTOR) :: distortionVec(nvec)
  real :: infallParticleMass         ! for T Tauri infall models
  integer, intent(in) :: maxBlobs
  logical, intent(in) :: flatSpec
  real :: inclination
  integer, intent(in) :: maxTau
  type(PHASEMATRIX),pointer :: miePhase(:,:, :)
  integer :: nsource
  type(SOURCETYPE) :: source(:)
  type(BLOBTYPE) :: blobs(:)
  integer, intent(in) :: nmumie
  real :: dtime

! local parameters
  integer, parameter :: maxscat = 1000000  ! Maximum number of scatterings for individual photon

! local variables 
  real(double) :: objectDistance
  real :: problinephoton
  type(VECTOR) :: originalViewVec
  type(VECTOR) :: rotationAxis
  integer :: nOuterLoop = 10
  integer :: ilambda
  logical :: lineResAbs    ! T if you want to include absorption
  !                        !   of line by line resonance zones.
  real :: dTau
  integer :: nTot  
  real :: imageSize
  real(double) :: fac
  real :: escProb
  real :: nuStart, nuEnd
  real(oct) :: t1, t2, t3
  integer(kind=bigint) :: nContPhotons
  type(VECTOR) :: rHat
  type(VECTOR) :: slitPosition
  real(double) :: tau_bnd
  real :: dlambda, thisTau
  real(double) ::thisTauDble
  real :: thislam
  integer :: currentScat
  logical :: redRegion
  logical :: contWindPhoton
  logical :: contPhoton
  real :: thisVel
  real :: obs_weight, tau_tmp, exp_minus_tau
  real :: nu
  integer           :: nOctals       ! number of octals in grid
  integer           :: nVoxels       ! number of unique voxels in grid
  integer(kind=bigInt) :: nInnerLoop
  type(PHOTON) :: testPhoton

  real(double) :: totLineEmission
  real(double) :: totContinuumEmission
  real(double) :: totCoreContinuumEmission
  real(double) :: totCoreContinuumEmission1
  real(double) :: totCoreContinuumEmission2
  real(double) :: totWindContinuumEmission

  type(VECTOR), parameter :: xAxis = VECTOR(1.,0.,0.)
  real :: wtot0_line =0., wtot0_cont = 0.
  real :: meanr0_line = 0., meanr0_cont = 0.
!  real :: meanr_line = 0., meanr_cont = 0.

  real(double) :: albedo
  real :: phi

  real, allocatable,save :: tauExt(:)
  real, allocatable,save :: tauAbs(:)
  real, allocatable,save :: tauSca(:)
  real, allocatable,save :: linePhotonAlbedo(:)
  real, allocatable,save :: lambda(:)

  real, allocatable,save :: contTau(:,:)

  real, allocatable :: splineArray(:)

  real, allocatable :: sourceSpectrum(:)
  real, allocatable :: sourceSpectrum2(:)

  type(IMAGETYPE) :: o6image(1)
  type(IMAGETYPE), allocatable :: obsImageSet(:)
  type(PVIMAGETYPE), allocatable :: pvimage(:)

  real :: r, r1, r2
  real(oct) :: t
  integer :: ntau, j
  real(double), allocatable        :: flux(:)

  type(STOKESVECTOR) :: yArray(nLambda)
  type(STOKESVECTOR) :: oldyArray(nLambda)
  type(STOKESVECTOR) :: yArrayStellarDirect(nLambda)
  type(STOKESVECTOR) :: yArrayStellarScattered(nLambda)
  type(STOKESVECTOR) :: yArrayThermalDirect(nLambda)
  type(STOKESVECTOR) :: yArrayThermalScattered(nLambda)
  type(STOKESVECTOR), allocatable :: varianceArray(:), errorArray(:,:)
  logical :: rotateView
  logical :: tiltView
  real :: rotateDirection

  integer :: i1, i2, i3

  integer :: iLambdaPhoton
  integer :: iOuterLoop
  real(double) :: lCore
  real :: weightDust=1.0, weightPhoto=1.0
  integer :: iInclination, nInclination
  integer :: iPhase

  real :: dopShift  ! velocity shift in terms of thermal line width

  real :: weightContPhoton, weightLinePhoton
  real :: chanceLine, chanceContinuum

  real, allocatable,save :: contWeightArray(:)

  real(double) :: thisChi, thisSca

  real(double) :: energyPerPhoton

  ! Spot stuff
  real :: chanceSpot                     ! chance of spot
  logical :: spotPhoton                  ! photon from spot?

  logical :: velocitySpace

  real :: chanceDust = 0.
  real(double) :: totDustContinuumEmission, totEnvelopeEmission

  integer :: nFromEnv
  integer(kind=bigInt) :: iInner_beg, iInner_end ! beginning and end of the innerPhotonLoop index.

  integer(kind=bigint) :: i
  integer :: iSlit, istep, ispline 
  type(VECTOR) :: viewVec, outVec, thisVec
  type(VECTOR) :: amrGridCentre

  integer :: toofewsamples   ! number of errors from integratePathAMR
  integer :: boundaryprobs   ! number of errors from integratePathAMR
  integer :: negativeOpacity ! number of errors from integratePathAMR

  character(len=80) :: tempChar
  character(len=80) :: phasePopFilename
  character(len=80) :: originalOutFile, filename
  character(len=80) :: specfile, obsfluxfile

  logical :: ok
  type(OCTAL), pointer :: sourceOctal, currentOctal, tempOctal
  integer :: sourceSubcell, currentSubcell, tempSubcell

#ifdef MPI
  ! For MPI implementations =====================================================
  real, dimension(:), allocatable :: tempRealArray
  real, dimension(:), allocatable :: tempRealArray2
  real(double), dimension(:), allocatable :: tempDoubleArray
  real(double), dimension(:), allocatable :: tempDoubleArray2
!  integer, dimension(:), allocatable :: photonBelongsRank
!  integer, parameter :: tag = 0
!  logical :: rankComplete
   integer ::   tempInt, ierr     
#endif

  ! O VI spectrum stuff

  character(len=80) :: o6filename
  integer, parameter :: no6pts = 100
  real, parameter :: o6start = 1031.5, o6end=1032.8
  real :: o6xarray(no6pts), o6yarray(no6pts)

  ! torus images

  integer           :: nImageLocal  ! number of images in obsImageSet
  type(filter_set)  :: filters ! a set of filters used for imaging
  character(LEN=30) :: name_filter
  real(double) :: lambda_eff  ! Effective wavelength of a filter[A]
  real(double) :: bandwidth   ! Band width of a filter[A]
  real(double) :: weightsource
  real :: probContPhoton

  ! intrinsic profile variables

  integer, parameter :: maxIntPro = 1000
  integer :: nIntPro
  real :: lamIntPro(maxIntPro), intPro(maxIntPro)
  character(len=80) :: message
#ifdef USECFITSIO
  character(len=80) :: header
#endif
  ! raman scattering model parameters
  type(VECTOR) :: ramanSourceVelocity

  real :: probDust

! Variables formerly in inputs_mod but not set in Torus V2 ----------------
  real :: ramVel=0.0
  character(len=20) :: ramanDist          ! raman distortion type
  logical :: coreEmissionLine=.false.
  real :: relIntCoreEmissionLine       ! previously: 'relIntCoreEmission'
  real :: velWidthCoreEmissionLine     ! previously: 'velWidthCoreEmission'
  real :: phaseOffset
  character(len=80) :: contFluxFile2
  real :: temp1
  character(len=80) :: popFilename
  character(len=80) :: distortionType=" "
  logical :: doIntensivePeelOff=.false.
  logical :: forceFirstScat=.false.
  integer :: nLower, nUpper
  logical :: starOff=.false. 
  real :: logMassLossRate
  real :: slitPA, slitWidth, slitLength
  real :: vfwhm=1.0, pfwhm=1.0
  type(VECTOR) :: slitPosition1=VECTOR(1.0,1.0,1.0), slitPosition2=VECTOR(1.0,1.0,1.0)
  integer :: nSlit=1, np
  real :: usePhotonWavelength
  logical :: forcedWavelength=.false.
  character(len=80) :: opacityDataFile
  logical :: fillTio=.false.
  logical :: plezModelOn=.false.
  logical :: fillThomson=.false.
  logical :: screened=.false.
  logical :: VoigtProf=.false.
  logical :: photLine=.false.            ! photospheric line production
  logical :: useInterp=.false.
#ifdef MPI
  logical :: readFileFormatted  ! whether 'grid' input  file is formatted
#endif
  logical :: writeFileFormatted ! whether 'grid' output file is formatted
  logical :: sphericityTest=.false.        ! sphericity test
  logical :: fillRayleighOpacity =.false.  ! previously: 'fillRayleigh'
  logical :: doPVimage=.false.             ! previously: 'pvimage'
  logical :: noPhaseUpdate=.false.  ! disable updating AMR grid at each phase
!-------------------------------------------------------------------------------

!$OMP THREADPRIVATE(lambda, tauExt, tauSca, tauAbs, contTau, contWeightArray)

  probDust = 0.1
  phaseTime = 0.0
  phaseOffset = 0.0
  dopshift = 0.
  ok = .true.
  sourcesubcell = 0; spotPhoton = .false.
  call define_rotation_axis

  biasPhiDirection = -1.d0

  objectDistance = griddistance * pctocm

  if ( grid%geometry == "cmfgen" ) then 
     probContPhoton = 0.2
  else
     probContPhoton = 1.0
  endif

  probLinePhoton = 1. - probContPhoton

  if (grid%doRaman) then
     ramanSourceVelocity = (ramVel/cSpeed)*VECTOR(0.,-1., 0.)
     if (writeoutput) write(*,'(a,f7.1,a)') "Raman source speed: ",modulus(ramanSourceVelocity)*cSpeed/1.e5, " km/s"
  endif

! thermalElectronVelocity is used by scatterPhoton 
  if (modulus(thermalElectronVelocity(10000.)) == 0.) then
     if (writeoutput) write(*,'(a)') "THERMAL ELECTRON BROADENING IS OFF!!!!!!!!!!!!!!!!!!!!!!!!!!"
  endif

  lineResAbs = .false.  

  if (lineEmission) then
     velocitySpace = .true.
  else
     velocitySpace = .false.
  endif

  ! hardwired stuff
  do i = 1, no6pts
     o6xArray(i) = o6start + (o6end-o6start)*real(i-1)/real(no6pts-1)
     o6yarray(i) = 1.e-10
  enddo

  if (mie.or.dustPhysics) nOuterLoop = nLambda

  allocate(errorArray(nOuterLoop,1:nLambda))
  allocate(varianceArray(1:nLambda))
  allocate(sourceSpectrum(1:nLambda))
  allocate(sourceSpectrum2(1:nLambda))

  sourceSpectrum = 1.
  sourceSpectrum2 = 1.

  if (dopvimage) then
     allocate(pvimage(1:nSlit))
  endif

  amrGridCentre = VECTOR(amrGridCentreX, amrGridCentreY, amrGridCentreZ)

  ! Choose whether to rotate the view
  call choose_view( grid%geometry,   nPhase,          distortionType, doRaman, & ! Intent in
                    rotateview, rotateDirection, tiltView)                  ! Intent out


  ! prepare the filters for images here.
  if (stokesimage) then        
     ! bulid a filter to be used for imaging
     ! -- using a routine in filter_set_class
     call make_filter_set(filters, filter_set_name)
     ! write filter info to standard output and to a file
     if (myRankIsZero) then
        call info_filter_set(filters, "*") 
        call info_filter_set(filters, "info_filter_set.dat")
     endif
  end if


  if (grid%geometry(1:6) == "ttauri" .or. grid%geometry(1:9) == "luc_cir3d" .or. &
       grid%geometry == "cmfgen" .or. grid%geometry == "romanova") then
     call emptySurface(starSurface)
  end if

! From here we do multiple runs if required !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  phaseLoop: do iPhase = nStartPhase, nEndPhase
     
     grid%timeNow = phaseTime * real(iPhase-1)

     viewVec = originalViewVec
     outVec = (-1.d0)*viewVec
     thisVec = viewVec


     ! we rotate the view by an appropriate amount


     if (rotateView.and.(nPhase /= 1)) then
        if (writeoutput) write(*,'(a)') "Rotating view..."
        phi = real(-rotateDirection * twoPi * real(iPhase-1)/real(nPhase-1))
        viewVec =  arbitraryRotate(thisVec, dble(phi), rotationAxis)
        outVec = (-1.d0) * viewVec
        if (writeoutput) write(*,'(a,f5.2,a,f5.2,a,f5.2,a)') "View vector: (",viewVec%x,",", &
             viewVec%y,",",viewVec%z,")"
     endif

     if (phaseOffset /= 0.) then
        if (writeoutput) write(*,'(a,f5.2)') "Rotating by phase offset: ",phaseOffset
        phi = real(twoPi * phaseOffset)
        viewVec =  arbitraryRotate(viewVec, dble(phi), rotationAxis)
        outVec = (-1.d0) * viewVec
     endif


     if (tiltView.and.(nPhase /= 1)) then
        phi = real(-twoPi * real(iPhase-1)/real(nPhase-1))
        viewVec =  arbitraryRotate(thisVec, dble(phi), xAxis)
        outVec = (-1.d0) * viewVec
     endif


     ! refill the grids 
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!              
     if (.not.plezModelOn .and. .not. gridUsesAMR) then
        select case(grid%geometry)
        case("torus")
!           call fillGridTorus(grid, rho, rTorus, rOuter)
        case("sphere")
           call fillGridSpheriod(grid, rho, radius, kFac)
        case("ellipse")
           call fillGridEllipse(Grid,rho,  rMin, rMaj, rInner, teff)
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

        !        call fillGridMagneticAccretion(grid)
        case("ttwind")
        
        case("betacep")
           !        call fillGridBetaCep(grid)
        case("donati")
!           call fillGridDonati(grid, resonanceLine)
        case("puls")
        !        call fillGridPuls(grid, mDot, rStar, tEff, v0, vterm, beta, xfac)
           case("wind")
           
        case ("ttauri")
            call fillGridMagneticAccretion(grid,contfluxfile, popFileName, &
                     readPops, writePops, lte,  lamLine, Laccretion, Taccretion, sAccretion, &
                     curtains, dipoleOffset, nLower, nUpper)
                     
!       call fillGridWind(grid, mDot, rStar, tEff, v0, vterm, beta, &
!       lte, contFluxFile, writePops, readPops, popFilename, nLower, nUpper)

        case("resonance")
!           call fillGridResonance(grid, rCore, mDot, vTerm, beta, temp)

        case("wr104")
           continue

        case("cluster")
           ! do nothing
           continue
           
        case DEFAULT
           if (writeoutput) write(*,*) "! Unrecognised grid geometry: ",trim(grid%geometry)
           return
        end select
     endif


     if (grid%geometry(1:6) == "ttauri" .or. grid%geometry(1:9) == "luc_cir3d" .or. &
         grid%geometry == "romanova") then      
        ! Nu must be set again here since it is not assigned when the population/grid 
        ! file is read from a file!  (RK changed here.)
        nu = real(cSpeed / (lamLine * angstromtocm))
        call contread(contFluxFile, nu, coreContinuumFlux)
        call buildSphere(grid%starPos1, dble(grid%rCore), starSurface, 400, dble(teff), &
             source(1)%spectrum)
        if (grid%geometry == "ttauri") then
           call createTTauriSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion) 
        elseif (grid%geometry == "romanova") then
           call createTTauriSurface2(starSurface,  romData, nu, coreContinuumFlux,fAccretion) 
        else
           call createSurface(starSurface, nu, coreContinuumFlux,fAccretion)            
        end if           
        call testSurface(starSurface)

       if (grid%adaptive) then
         !call infallEnhancment(grid, distortionVec, nVec, nPhi, dTime, .true., &

         !                    infallParticleMass, alreadyDoneInfall)
          if (.not. noPhaseUpdate .and. iPhase /= nStartPhase) then
#ifdef MPI
 if (myRankGlobal == 0) then
#endif
             if (writeoutput) write(*,*) "Modifying grid for new phase"
             call amrUpdateGrid(limitScalar,limitScalar2,grid) 
             call countVoxels(grid%octreeRoot,nOctals,nVoxels)
             if (writeoutput) write(*,*) "Adaptive grid contains: ",nOctals," octals"
             if (writeoutput) write(*,*) "                      : ",nVoxels," unique voxels"
             grid%nOctals = nOctals

#ifdef MPI
    if (writePhasePops) then
        write(tempChar,'(i3.3)') iPhase
        phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
        call writeAMRgrid(phasePopFilename,writeFileFormatted,grid)
        print *,'Process ',myRankGlobal,' finished writing phase pop grid...' 
        print *,'Process ',myRankGlobal,' signalling grid has become available...' 
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)  ! XXX --> sync with XXX below
 else ! myRankGlobal /=0
    print *,'Process ',myRankGlobal,' waiting for grid to become available...' 
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)  ! XXX --> synch with XXX above
    if (writePhasePops) then  ! then file should should be avialable
       write(tempChar,'(i3.3)') iPhase
       phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
       call readAMRgrid(phasePopFilename,readFileFormatted,grid)
     end if
 end if
 call MPI_BARRIER(MPI_COMM_WORLD, ierr)    ! sync everybody here
#endif

             if (writeoutput) write(*,*) "Recalculating statistical equilibrium after changing grid" 
             call amrStateq(grid, lte, nLower, nUpper,  &
                   starSurface, recalcPrevious=.true.)
             call torus_mpi_barrier('returned from amrStatEq. Waiting to sync...')

             if (ttau_disc_on) then
                ! amrStateq will have messed up the disc, so we reset those cells
                call finish_grid(grid%octreeroot, grid, ttauri_disc, 1.0, &
                     sigmaAbs0, sigmaSca0, meanDustParticleMass)
             end if

             if (myRankIsZero .and. writePhasePops) then
                write(tempChar,'(i3.3)') iPhase
                phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
                call writeAMRgrid(phasePopFilename,writeFileFormatted,grid)
             end if

             call torus_mpi_barrier
             
          end if ! (.not. noPhaseUpdate .and. iPhase /= nStartPhase) 

       else
         call infallEnhancment(grid, distortionVec, nVec, nPhi, dTime, .true., &
                             infallParticleMass, alreadyDoneInfall)
      end if !(adaptive)

   end if

     if (.not. gridUsesAMR) then


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
!           call distortGridSpiral(grid, vRot, nSpiral)
        case("rotation")
!           call distortRotation(grid, vRot)
!           write(*,'(a,f5.1,a)') "Grid distorted by a rotational velocity of ",vRot/1.e5," km/s"
        case("test")
           call distortGridTest(grid)
!        case("raman")
!           call distortRaman(grid)
         case("raman")
            call  distortStrom(grid, secondSourcePosition, .true., .true., 0.5*rStar, coolStarPosition, ramanDist)
        case("wrdisk")
           call distortWRdisk(grid)
   
        case("windwind")
           call distortWindCollision(grid, momRatio, binarySep)
   
        end select
   
        ! we read in the blob configuration from a file

     end if ! (.not. gridUsesAMR)

     if (nBlobs > 0) then
        write(filename,"(a,i3.3,a)") "run",iPhase,".blob"
        call readBlobs(filename, maxBlobs, blobs, .false.)
        call distortGridWithBlobs(grid, maxBlobs, blobs)
     endif

     

     if (grid%geometry .eq. "rolf") then
        if (gridUsesAMR) then
           !      is this correct vvvvvvvvvvvvvvvvvvvvvv ?
           secondSourcePosition = grid%octreeRoot%centre
        else
           secondSourcePosition = VECTOR(grid%xAxis(grid%nx/2), &
             grid%yAxis(grid%ny/2), &
             grid%zAxis(grid%nz/2))
        end if
     endif

!     call  createDustCrossSectionPhaseMatrix(grid, grid%lamArray, nLambda, miePhase, nMuMie)
!     call setupDust(grid, grid%lamArray, nLambda, miePhase, nMumie, fileStart="sed")
     if (writeoutput) then
        write(*,*) "nlambda ",nlambda
        write(*,*) "grid%nlambda ",grid%nLambda
        write(*,*) "size(grid%lamArray) ",size(grid%lamarray)
        write(*,*) "grid%lamArray ",grid%lamArray
     endif

     if (noScattering) then
        if (writeoutput) write(*,*) "! WARNING: Scattering opacity turned off in model"
        grid%oneKappaSca(1:nDustType,1:nLambda) = TINY(grid%oneKappaSca)
     endif


     allocate(lambda(1:maxTau))
     allocate(tauExt(1:maxTau))
     allocate(tauAbs(1:maxTau))
     allocate(tauSca(1:maxTau))
     allocate(linePhotonalbedo(1:maxTau))
     if (grid%lineEmission) allocate(contTau(1:maxTau,1:nLambda))
     allocate(contWeightArray(1:nLambda))


!

     if (grid%geometry == "wr104") then
        write(*,*) "stripping dust!!!!!!!!!!!"
        call stripDustAway(grid%octreeRoot, 1.d-20, 4.d4)
     endif


     call writeInfo(" ", TRIVIAL)
     call writeInfo("Some basic model parameters",TRIVIAL)
     call writeInfo("---------------------------",TRIVIAL)
     call writeInfo(" ", TRIVIAL)

     call writeFormatted("(a,f7.1)","Inclination: ",real(radtodeg*acos(outVec%z)),TRIVIAL)
     ! THE FOLLOWING STATEMENT IS DANGEROUS. ACCORDING TO INPUT_MOD.F90 
     ! NOT ALL THE GEOMETRY HAS RCORE VALUES, AND SAME UNITS!
     ! THIS SHOULD BE DONE IN INITAMRGRID ROUTINE AS SOME GEOMETRY HAS
     ! DONE SO.
     if (grid%geometry /= "cmfgen") then
        grid%rStar1 = rcore/1.e10
     end if
     !
     ! Performs various optical depth tests here...
     !

     if ( (.not. hydrodynamics) .and. (.not. formalsol)) then
        call test_optical_depth(gridUsesAMR, VoigtProf, &
             amrGridCentre, sphericityTest,  &
             outVec, lambdatau,  lambdatau, grid, thin_disc_on, opaqueCore,  &
             thinLine, lineResAbs, nUpper, nLower, useinterp, grid%Rstar1, coolStarPosition, maxTau, nSource, source)
      end if


      call torus_mpi_barrier

!      write(*,*) "calling scattering test"
!      call testscatterPhoton(grid, miePhase, nDustType, nLambda, grid%lamArray, nMuMie)


     weightLinePhoton = 0.
     weightContPhoton = 1.
     weightPhoto = 1.

     if (mie .and. (.not. useDust)) then

        call calcContinuumEmissivityLucy(grid, grid%octreeRoot , nlambda, grid%lamArray)

        call computeProbDist(grid, totLineEmission, &
             totDustContinuumEmission,lamline, .false.)
        totDustContinuumEmission = totdustContinuumEmission 
        lcore = grid%lCore
        if (nSource > 0) then
           lCore = sumSourceLuminosity(source, nsource, grid%lamArray(1), grid%lamArray(nLambda))
        endif

        totEnvelopeEmission = totDustContinuumEmission
        chanceDust = real(totDustContinuumEmission/(totDustContinuumEmission+lCore/1.e30))
!        if (writeoutput) write(*,*) "totdustemission",totdustcontinuumemission
!        if (writeoutput) write(*,'(a,f7.2)') "Chance of continuum emission from dust: ",chanceDust

        weightDust = chanceDust / probDust
        weightPhoto = (1. - chanceDust) / (1. - probDust)

!        if (writeoutput) write(*,*) "WeightDust",weightDust
!        if (writeoutput) write(*,*) "WeightPhoto",weightPhoto
!        if (writeoutput) write(*,*) "core + envelope luminosity",lCore+totEnvelopeEmission*1.d30
        energyPerPhoton =  ((lCore + totEnvelopeEmission*1.d30) / dble(nPhotons))/1.d20
!        if (writeoutput) write(*,*) "Energy per photon: ", energyPerPhoton

     endif

     



     if (grid%geometry == "hourglass") then
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
        if (writeoutput) write(*,*) "Total Raman Line Emission: ",totLineEmission
        weightLinePhoton = 1.
        weightContPhoton = 0.
        probLinePhoton = 1.
        probContPhoton = 0.
     endif



     if (lineEmission) then

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


!        write(*,'(a)') "mu axis"
!        do i = 1, grid%nMu
!           write(*,*) i,acos(grid%muAxis(i))*radtodeg, &
!                             grid%muProbDistLine(grid%nr/2,i)
!        enddo
!        write(*,'(a)') "Phi axis"
!        do i = 1, grid%nPhi
!           write(*,*) i,grid%phiAxis(i)*radtodeg, grid%phiProbDistLine(grid%nr/2,grid%nmu/2,i)
!        enddo



        if (grid%geometry == "donati") totWindContinuumEmission = 0.

        if (grid%geometry == "binary") totWindContinuumEmission = 0.

        if (grid%resonanceLine) totWindContinuumEmission = 0.



                    
!              ! New photon position 
!              thisPhoton%position = thisPhoton%position + real(dlambda,kind=oct)*thisPhoton%direction
!              ! adjusting the photon weights 
!              if ((.not. mie) .and. (.not. thisPhoton%linePhoton)) then
!                 if (j < nTau) then
!                    contWeightArray(1:nLambda) = contWeightArray(1:nLambda) *  &
!                         EXP(-(contTau(j,1:nLambda) + t*(contTau(j+1,1:nLambda)-contTau(j,1:nLambda))) )
!                 else
!                    contWeightArray(1:nLambda) = contWeightArray(1:nLambda)*EXP(-(contTau(nTau,1:nLambda)))
!                 end if
!              end if
!
!
!              if (flatspec) then
!                 iLambda = 1
!              else
!                 iLambda = findIlambda(thisPhoton%lambda, grid%lamArray, nLambda, ok)
!              endif
!
!              if (grid%adaptive) then
!                 positionOc = thisPhoton%position
!                 call amrGridValues(grid%octreeRoot, positionOc, grid=grid, iLambda=iLambda, &
!                      kappaAbs = thisChi, kappaSca = thisSca)
!              else
!                 call getIndices(grid, thisPhoton%position, i1, i2, i3, t1, t2, t3)
!                 if (.not.grid%oneKappa) then
!                    if (.not.flatspec) then
!                       thisChi = interpGridKappaAbs(grid, i1, i2, i3, iLambda, t1, t2, t3)
!                       thisSca = interpGridKappaSca(grid, i1, i2, i3, iLambda, t1, t2, t3)
!                    else
!                       thisChi = interpGridKappaAbs(grid, i1, i2, i3, 1, t1, t2, t3)
!                       thisSca = interpGridKappaSca(grid, i1, i2, i3, 1, t1, t2, t3)
!                    endif
!                 else
!                    r = interpGridScalar2(grid%rho,grid%na1,grid%na2,grid%na3,i1,i2,i3,t1,t2, t3)
!                    thisChi = grid%oneKappaAbs(1,iLambda) * r
!                    thisSca = grid%oneKappaSca(1,iLambda) * r
!                 endif
!              end if
!              
!              if (contPhoton) then
!                 if ((thisChi+thisSca) >= 0.) then
!                    albedo = thisSca / (thisChi + thisSca)
!                 else
!                    albedo = 0.
!                    write(*,*) "Error:: thisChi+thisSca < 0 in torusMain."
!                    !                 stop
!                 endif
!              else
!                 albedo = linePhotonAlbedo(j)
!              endif

        !if (grid%geometry == "ttauri" .or. grid%geometry == "windtest") then
        if (grid%geometry == "windtest") then
           totWindContinuumEmission = 0.
           if (writeoutput) write(*,'(a)') "! Wind continuum emission switched off."
        endif

        if (lineOff) then
           totLineEmission = 0.
        endif

        nu = real(cSpeed / (lamLine * angstromtocm))
        if (writeoutput) write(*,'(a,e12.3)') "Line emission: ",totLineEmission
        totCoreContinuumEmission = 0.d0
        select case(grid%geometry)
           case("binary")
              call contread(contFluxFile, nu, totCoreContinuumEmission1)
              call contread(contFluxFile2, nu, totCoreContinuumEmission2)
           case("puls")
              totCoreContinuumEmission = pi * blackBody(0.77 * tEff, lamLine)
!           case DEFAULT
!              call contread(contFluxFile, nu, coreContinuumFlux)
!              totCoreContinuumEmission = coreContinuumFlux
        end select


        ! reads the intrinsic core absorption profile for core-halo models
        ! Note: The profile should be normaised (to contiuum) and the x axis shoudl be in Angstrome.

        if (trim(intProFilename) /= "none") then
           call rdIntPro(intProFilename, intPro, lamIntPro, nIntPro)
           allocate(splineArray(1:nIntPro))
           call spline(lamIntPro, intPro, nIntPro, 0., 0., splineArray)
           do iSpline = 1 , nLambda
              if ((grid%lamArray(iSpline) < lamIntPro(1)) .or.  &
                   (grid%lamArray(iSpline) > lamIntPro(nIntPro))) then
                 sourceSpectrum(iSpline) = 1.
              else
                 call splint(lamIntPro, intPro, splineArray, nIntPro, &
                      grid%lamArray(iSpline), sourceSpectrum(iSpline))
              endif
           enddo
           deallocate(splineArray)
        endif

        if (grid%geometry == "binary") then
           call rdIntPro(intProFilename2, intPro, lamIntPro, nIntPro)
           allocate(splineArray(1:nIntPro))
           call spline(lamIntPro, intPro, nIntPro, 0., 0., splineArray)
           do iSpline = 1 , nLambda
              if ((grid%lamArray(iSpline) < lamIntPro(1)) .or.  &
                   (grid%lamArray(iSpline) > lamIntPro(nIntPro))) then
                 sourceSpectrum2(iSpline) = 1.
              else
                 call splint(lamIntPro, intPro, splineArray, nIntPro, &
                      grid%lamArray(iSpline), sourceSpectrum2(iSpline))
              endif
           enddo
           deallocate(splineArray)
        endif

        ! the core can emit as a Gaussian emission line

        if (coreEmissionLine) then
           call computeCoreEmissionProfile(grid%lamArray, sourceSpectrum, nLambda, &
                lamLine, velWidthCoreEmissionLine, relIntCoreEmissionLine)
        endif


        nuStart = real(cSpeed / (grid%lamArray(1) * angstromtocm))
        nuEnd = real(cSpeed / (grid%lamArray(grid%nLambda) * angstromtocm))

        totWindContinuumEmission = totWindContinuumEmission * (nuStart - nuEnd)
        if (writeoutput) write(*,'(a,e12.3)') "Wind cont emission: ",totWindContinuumEmission
        if (writeoutput) write(*,*) grid%lamArray(1), &
             grid%lamArray(grid%nlambda)


! factor of four pi taken out from below - input continuum files are expected
! to be in fluxes (erg/s/cm^2/hz) not Hnu's

        if (grid%geometry /= "binary") then
           totCoreContinuumEmission = (totCoreContinuumEmission * &
                (nuStart - nuEnd))*fourPi*grid%rCore**2
        else
           totCoreContinuumEmission1 = totCoreContinuumEmission1 * (nuStart-nuEnd) * fourPi * grid%rStar1**2
           totCoreContinuumEmission2 = totCoreContinuumEmission2 * (nuStart-nuEnd) * fourPi * grid%rStar2**2
           totCoreContinuumEmission = totCoreContinuumEmission1 + totCoreContinuumEmission2
           grid%lumRatio = real(totCoreContinuumEmission1 / totCoreContinuumEmission)

           if (writeoutput) write(*,*) "Binary luminosity ratio (p/s): ",totCoreContinuumEmission1/totCoreContinuumEmission2
        endif
           
        if (grid%geometry == "ttauri") then
           write (*,'(a,e12.3)') 'T Tauri star: continuum emission: ',totCoreContinuumEmission
           fAccretion = fAccretion * (nuStart-nuEnd)
           write (*,'(a,e12.3)') 'T Tauri accretion: continuum emission: ',fAccretion
           write (*,'(a,f12.3)') 'Accretion continuum / stellar continuum: ',fAccretion/totCoreContinuumEmission
           totCoreContinuumEmission = totCoreContinuumEmission + fAccretion
           write (*,'(a,e12.3)') 'T Tauri total core continuum emission: ',totCoreContinuumEmission
        endif

        if (grid%geometry == "luc_cir3d".or. grid%geometry == "romanova") then
           write (*,'(a,e12.3)') 'Star: continuum emission: ',totCoreContinuumEmission
           fAccretion = fAccretion * (nuStart-nuEnd)
           write (*,'(a,e12.3)') 'Accretion: continuum emission: ',fAccretion
           if (fAccretion /= 0.d0) then
              write (*,'(a,f12.3)') 'Accretion continuum / stellar continuum: ',fAccretion/totCoreContinuumEmission
              totCoreContinuumEmission = totCoreContinuumEmission + fAccretion
           endif
        endif


        
        chanceSpot = 0.
        if ((grid%geometry == "disk").and.(nSpot > 0)) then
           chanceSpot = fSpot * blackBody(tSpot, 6562.8) / &
                ((1.-fSpot)*blackBody(tEff, 6562.8) + &
                (fSpot * blackBody(tSpot, 6562.8)))
           if (writeoutput) write(*,'(a,f6.3)') "Spot chance at 6563A: ",chanceSpot
        endif




        if (writeoutput) write(*,'(a,e12.3)') "Total Core Continuum Emission: ",totCorecontinuumEmission


        totContinuumEmission = totCoreContinuumEmission + totWindContinuumEmission

        if (writeoutput) write(*,'(a,e12.3)') "Total Continuum emission: ", totContinuumEmission


        if ((totContinuumEmission + totLineEmission) /= 0.) then
           chanceLine = real(totLineEmission/(totContinuumEmission + totLineEmission))
           chanceContinuum = real(totContinuumEmission / &
                (totContinuumEmission + totLineEmission))
           grid%chanceWindOverTotalContinuum = real(totWindContinuumEmission &
                / max(1.d-30,totContinuumEmission))
        else
           chanceLine =0.
           chanceContinuum = 1.
           grid%chanceWindOverTotalContinuum = 0.
        endif



        if (writeoutput) write(*,'(a,e12.3)') "Chance continuum emission in wind: ", grid%chanceWindOverTotalContinuum

        write(*,*) "nlambda ",nlambda
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

        if (writeoutput) then
           write(*,*) "Line photon weight: ", weightLinePhoton
           write(*,*) "Continuum photon weight: ", weightContPhoton
           
           write(*,*) " "
           write(*,*) "Line photon prob: ", probLinePhoton
           write(*,*) "Continuum photon prob: ", probContPhoton
        endif

           energyPerPhoton =  (totLineEmission + totContinuumEmission) / dble(nPhotons)

        if (writeoutput) then
           write(*,*) "Energy per photon: ", energyPerPhoton
           write(*,*) "chance line",chanceline
        endif


        ! compute the probability distributions


     endif
     if (.not. formalsol) then
        deallocate(lambda)
        deallocate(tauSca)
        deallocate(tauExt)
        deallocate(tauAbs)
        if (allocated(contTau)) deallocate(contTau)
        deallocate(linePhotonAlbedo)
        deallocate(contWeightArray)
     end if

     nContPhotons = nint((probContPhoton * real(nPhotons) / real(nOuterLoop)))
     if (writeoutput) write(*,*) "Number of continuum photons: ",nContPhotons


     if (formalSol) then
        if (myrankglobal == 0) then
           write(*,*) "calling create lucy image"
           call createLucyImage(grid, viewVec, 1.e4, grid%lamArray, nLambda, source, nSource)
        endif
        call torus_mpi_barrier
        stop
     endif


     !
     ! here we may loop over different inclinations
     !

! Use the override inclination values if supplied, otherwise use values from sed_mod
     if ( present(overrideInclinations) ) then
        nInclination = size(overrideInclinations)
     else
        nInclination = getNumSedInc()
     endif

  incLoop: do iInclination = 1, nInclination, 1 
       ! NB This loop is not indented in the code!
     
     if (nInclination >= 1) then  

        if ( present(overrideInclinations) ) then
           inclination = overrideInclinations(iInclination)
        else
           inclination = getSedInc(iInclination)
        end if
        inclination = max(inclination,1.e-4)

        if (writeoutput) then
           write(message,*) " "
           call writeInfo(message, TRIVIAL)
           write(message,*) "Inclination = ",inclination*radToDeg,' degrees'
           call writeInfo(message, TRIVIAL)
           write(message,*) " "
           call writeInfo(message, TRIVIAL)
        end if

       if (iPhase == nStartPhase .and. iInclination == 1) originalOutFile = outFile
         
       write(tempChar,'(i3.3)') NINT(inclination*radToDeg)
       outFile = trim(originalOutFile)//'_inc'//TRIM(tempChar)
       
       viewVec%x = 0.
       viewVec%y = -sin(inclination)
       viewVec%z = -cos(inclination)
       outVec = (-1.d0)*viewVec
       thisVec = viewVec

       if (rotateView.and.(nPhase /= 1)) then
         write(*,'(a)') "Rotating view..."
         phi = real(-rotateDirection * twoPi * real(iPhase-1)/real(nPhase-1))
         viewVec =  arbitraryRotate(thisVec, dble(phi), rotationAxis)
         outVec = (-1.d0) * viewVec
         write(*,'(a,f5.2,a,f5.2,a,f5.2,a)') "View vector: (",viewVec%x,",", &
         viewVec%y,",",viewVec%z,")"
       endif

       if (phaseOffset /= 0.) then
         write(*,'(a,f5.2)') "Rotating by phase offset: ",phaseOffset
         phi = real(twoPi * phaseOffset)
         viewVec =  arbitraryRotate(viewVec, dble(phi), rotationAxis)
         outVec = (-1.d0) * viewVec
       endif

       if (tiltView.and.(nPhase /= 1)) then
         phi = real(-twoPi * real(iPhase-1)/real(nPhase-1))
         viewVec =  arbitraryRotate(thisVec, dble(phi), xAxis)
         outVec = (-1.d0) * viewVec
       endif

    end if

     write(message,'(a,f6.3,a,f6.3,a,f6.3,a)') "Viewing vector: (",viewVec%x,",",viewVec%y,",", viewVec%z,")"
     call writeInfo(message, TRIVIAL)

     ! zero the output arrays

     yArray(1)%i = 0.
     yArray(1)%q = 0.
     yArray(1)%u = 0.
     yArray(1)%v = 0.
     do i = 2, nLambda
        yArray(i)%i = 0.
        yArray(i)%q = 0.
        yArray(i)%u = 0.
        yArray(i)%v = 0.
     enddo
     oldyArray = STOKESVECTOR(0., 0., 0., 0.)

     yArrayStellarDirect(:) = STOKESVECTOR(0., 0., 0., 0.)
     yArrayThermalDirect(:) = STOKESVECTOR(0., 0., 0., 0.)
     yArrayStellarScattered(:) = STOKESVECTOR(0., 0., 0., 0.)
     yArrayThermalScattered(:) = STOKESVECTOR(0., 0., 0., 0.)

     if (doRaman) then
        yArray(1:nLambda)%i = 1.e-20
     endif
     
     if (stokesimage .and. .not. formalsol) then
!! THIS SECTION MOVED OUTSIDE OF PHASE LOOP SINCE
!!  THIS SHOULD BE DONE ONLY ONCE!
!!        ! bulid a filter to be used for imaging
!!        ! -- using a routine in filter_set_class
!!        call make_filter_set(filters, filter_set_name)
!!
!!        ! write filter info to standard output and to a file
!! if (myRankIsZero) then
!!        call info_filter_set(filters, "*") 
!!        call info_filter_set(filters, "info_filter_set.dat")
!!     endif 
!!
        ! number of images = number of filters
        nImageLocal = get_nfilter(filters)
        ! Allocate the image array
        if (allocated(obsImageSet)) deallocate(obsImageSet)
        allocate(obsImageSet(nImageLocal))
        

        ! Initializing the images ...

        if (setImageSize == 0.) then
           if (grid%adaptive) then
              if (amr2d) then
                 imageSize = real(4.*grid%octreeRoot%subcellSize)
              else
                 imageSize = real(2.*grid%octreeRoot%subcellSize          )
              endif
           else if (grid%cartesian) then
              imageSize = grid%xAxis(grid%nx)-grid%xAxis(1)
           else
              imagesize = 2.*grid%rAxis(grid%nr)
           endif
        else
           imageSize = setImageSize / 1.e10
        endif

        if (imageInArcsec.and.(setImageSize /= 0.)) then
           imagesize = real(objectdistance * (setImageSize / 3600.) * degtorad / 1.e10)
        endif

!     if (molecular) then
!        if (writemol) call  molecularLoop(grid, co)
!        call calculateMoleculeSpectrum(grid, co)
!        call createDataCube(cube, grid, VECTOR(0.d0, 1.d0, 0.d0), co, 1)
!        if (myRankIsZero) call plotDataCube(cube, 'cube.ps/vcps')
!        stop
!     endif


        ! Initializing the images ...
        do i = 1, nImageLocal           
           if (grid%cartesian) then
              obsImageSet(i) = initImage(npixels, npixels, imageSize, imageSize, vmin, vmax)
              if (doRaman) then
                 o6image(1) = initImage(npixels, npixels, imageSize, imageSize, vmin, vmax)
              endif
           else if (grid%adaptive) then
              obsImageSet(i) = initImage(npixels, npixels, imageSize, imageSize, vmin, vmax)
           else   
              select case (grid%geometry)
              case("disk")
                 obsImageSet(i) = initImage(npixels, npixels, 2.*grid%rAxis(grid%nr), 2.*grid%rAxis(grid%nr),vMin, vMax)
              case("flared")
                 obsImageSet(i) = initImage(npixels, npixels, 8.*grid%rAxis(1),8.*grid%rAxis(1), vMin, vMax)
              case DEFAULT
                 obsImageSet(i) = initImage(npixels,npixels,  min(5.*grid%rAxis(1),grid%rAxis(grid%nr)), &
                      min(5.*grid%rAxis(1),grid%rAxis(grid%nr)), vMin, vMax)
              end select
           endif
        end do
        
     endif





     if (dopvImage) then

        if (nSlit  > 1) then
           rHat = slitPosition2 - slitPosition1
           do iSlit = 1, nSlit
              slitPosition = slitPosition1 + (dble(iSlit-1)/dble(nSlit-1))*(slitPosition2 - slitPosition1)
              pvImage(iSlit) = initPVimage(nv, vMin, vMax, np, -slitLength/2., slitLength/2., &
                   slitPosition, slitPA, slitWidth, slitLength)
           enddo
        else
           slitPosition = slitPosition1
           pvImage(1) = initPVimage(200, vMin, vMax, 200, -slitLength/2., slitLength/2., &
                slitPosition, slitPA, slitWidth, slitLength)
        endif
           
     endif

     
     write(message,*) " "
     call writeInfo(message, TRIVIAL)
     write(message,'(a)') "Run-time messages"
     call writeInfo(message, TRIVIAL)
     write(message,'(a)') "-----------------"
     call writeInfo(message, TRIVIAL)
     write(message,*) " "
     call writeInfo(message, TRIVIAL)
     ntot = 0

     ! now we loop 10 times over a tenth of the photons - this will be used
     ! to help compute the errors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




     !====================================================================================
     ! Perform formal integration
     !===================================================================================
     if ((grid%geometry == "ttauri" .or. grid%geometry == "romanova") &
          .and. formalsol) then
        if (doTuning) call tune(6, "One Formal Sol") ! start a stopwatch  
        write(*,*) "Started formal integration of flux...."
        if (.not. ALLOCATED(flux)) ALLOCATE(flux(grid%nlambda))
        if (nInclination >= 1) then
           write(tempChar,'(i3.3)') NINT(inclination*radToDeg)
           obsfluxfile = trim(originalOutFile)//'_inc'//TRIM(tempChar)
        else
           obsfluxfile = trim(OutFile)
        end if
        if (nphase > 1) then
           write(tempChar,'(i3.3)') iPhase
           obsfluxfile = trim(obsfluxfile)//TRIM(tempChar)
        end if

	! As of 04-aug-2006, this routine does not include the
        ! scattering flux!  Should be modified later.  (RK)
        call compute_obs_line_flux(lamline, REAL(mHydrogen), DBLE(grid%rstar1), &
             dble(TTauriRouter/1.0e10), dble(amrGridSize)/2.001d0/SQRT(2.0d0), &
             starSurface, &
             form_nr_wind, form_nr_acc, form_nr_core, form_nphi,  &
             outVec, objectDistance, grid, sampleFreq, opaqueCore,  &
             flux, grid%lamArray, grid%nlambda, obsfluxfile, &
             thin_disc_on, (ttau_disc_on .and. .not. ttau_turn_off_disc), &
             (ttau_jet_on .and. .not. ttau_turn_off_jet), ttau_discwind_on, npixels, &
             (ttau_disc_on .and. .not. ttau_turn_off_disc), ttauri_disc, thinLine, lineOff, &
             do_pos_disp) 
        write(*,*) "...Finished  formal integration of flux...."
        write(*,*) " "
        if (doTuning) call tune(6, "One Formal Sol") ! stop a stopwatch
        ! jumps to the end of inclination loops  (this should be replaced with while loop later)
        goto 777  
     end if

     iLambda = findIlambda(1.e5, grid%lamArray, nLambda, ok)
     if (doTuning) call tune(6,"Calculate bias on tau")
     if (mie)     call setBiasOnTau(grid, iLambda)
     if (doTuning) call tune(6,"Calculate bias on tau")


     !  These should be zero-ed for each viewing angle!
     tooFewSamples = 0 
     boundaryProbs = 0 
     negativeOpacity = 0 

     if (doTuning) call tune(6, "All Photon Loops")  ! Start a stopwatch
     
     call randomNumberGenerator(randomSeed=.true.)

     weightSource = 1.d0
     if (nSource > 0) &
     call randomSource(source, nSource, j, weightSource,grid%lamArray, nLambda, initialize=.true.)  

     if (mie.or.photoionization.or.lineEmission) then
        nInnerLoop = nPhotons / nOuterLoop
! Trap invalid ninnerloop. Can arise if nphotons < nouterloop
        if (ninnerloop < 1 ) then
           call writeFatal("phaseloop_mod: ninnerloop <1")
        end if
     endif

#ifdef MPI
     call test_random_across_threads()
#endif

#ifdef _OPENMP

#else
           allocate(lambda(1:maxTau))
           allocate(tauExt(1:maxTau))
           allocate(tauAbs(1:maxTau))
           allocate(tauSca(1:maxTau))
           allocate(linePhotonalbedo(1:maxTau))
           allocate(contTau(1:maxTau,1:nLambda)) 
           allocate(contWeightArray(1:nLambda))
#endif

     outerPhotonLoop: do iOuterLoop = 1, nOuterLoop

        if (mie) then

           iLambdaPhoton = iOuterLoop

!           if (rGapOuter > 0.d0) then
!              uhat = VECTOR(1.d0, 0.d0, 0.d0)
!              uHat = rotateZ(uHat, 1.d-4)
!              call tauAlongPath(ilambdaPhoton, grid, VECTOR(1.01d0*rGapInner,0.d0,0.d0), uHat , tau=thistaudble, &
!                   stopatdistance=dble(rGapOuter))
!              write(message,*) "Tau to gap at ",grid%lamArray(ilambdaPhoton),"angstroms: ",thistaudble
!              call writeInfo(message, TRIVIAL)
!           endif


           call calcContinuumEmissivityLucyMono(grid, grid%octreeRoot, grid%lamArray, &
                grid%lamArray(ilambdaPhoton), iLambdaPhoton)
           
!           if (doTuning) call tune(6,"Calculate bias on tau")
!           call setBiasOnTau(grid, iLambdaPhoton)
!           if (doTuning) call tune(6,"Calculate bias on tau")

           call computeProbDist(grid, totLineEmission, &
                totDustContinuumEmission,lamline, .false.)


           totDustContinuumEmission = totdustContinuumEmission 
           lcore = grid%lCore
           if (nSource > 0) then              
              if (.not.starOff) then
                 lCore = sumSourceLuminosityMonochromatic(grid, source, nsource, dble(grid%lamArray(iLambdaPhoton)))
                 if (writeoutput) write(*,*) "Core luminosity is: ",lcore, " erg/s/A ", lcore/(fourpi * objectDistance**2)
              else
                 lcore = tiny(lcore)
              endif
           endif

           totEnvelopeEmission = totDustContinuumEmission
           chanceDust = real(totDustContinuumEmission/(totDustContinuumEmission+lCore/1.e30))


!           if (writeoutput) write(*,*) "totdustemission",totdustcontinuumemission
!           if (writeoutput) write(*,*) "totcontemission",lcore/1.e30
!           if (writeoutput) write(*,'(a,f9.7)') "Chance of continuum emission from dust: ",chanceDust
           




           probDust = chanceDust
           weightDust = 1.
           weightPhoto = 1.

           if (chanceDust > 0.99) then
              probDust = 0.9
              weightDust = chanceDust / probDust
              weightPhoto = (1. - chanceDust) / (1. - probDust)
           endif


           energyPerPhoton =  (totDustContinuumEmission*1.d30 + lCore)/1.d20/dble(nInnerLoop)
!           if (writeoutput) write(*,*) "WeightDust",weightDust
!           if (writeoutput) write(*,*) "WeightPhoto",weightPhoto
!           if (writeoutput) write(*,*) "core + envelope luminosity",lCore+totEnvelopeEmission*1.d30
!           if (writeoutput) write(*,*) "Energy per photon: ", energyPerPhoton
!           if (writeoutput) write(*,*) "nInnerloop ",nInnerloop

        endif


        if (doTuning) call tune(6, "One Outer Photon Loop") ! Start a stop watch

        call do_one_outer_photon_loop

        if (doTuning) call tune(6, "One Outer Photon Loop") ! Stop a stop watch        

!        yArray(1:nLambda) = STOKESVECTOR(0.,0.,0.,0.)
        do i = 1, nLambda
           errorArray(iOuterLoop,1:nLambda) = yArray(i) - oldyArray(i)
        enddo
         oldyArray = yArray  
  


     end do outerPhotonLoop ! outer photon loop

#ifdef _OPENMP

#else
   deallocate(lambda)
   deallocate(tauSca)
   deallocate(tauExt)
   deallocate(tauAbs)
   deallocate(linePhotonalbedo)
   deallocate(contTau)
   deallocate(contWeightArray)
#endif


     if (doTuning) call tune(6, "All Photon Loops")  ! Stop a stopwatch

#ifdef MPI
     tempInt = 0
     call MPI_REDUCE(tooFewSamples,tempInt,1,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,ierr)
     tooFewSamples = tempInt

     call MPI_REDUCE(boundaryProbs,tempInt,1,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,ierr)
     boundaryProbs = tempInt

     call MPI_REDUCE(negativeOpacity,tempInt,1,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,ierr)
     negativeOpacity = tempInt

     call MPI_REDUCE(nTot,tempInt,1,MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD,ierr)
     nTot = tempInt

 if (stokesimage) then
   do i = 1, nImageLocal
     allocate(tempRealArray(SIZE(obsImageSet(i)%pixel)))
     allocate(tempRealArray2(SIZE(obsImageSet(i)%pixel)))
     allocate(tempDoubleArray(SIZE(obsImageSet(i)%pixel)))
     allocate(tempDoubleArray2(SIZE(obsImageSet(i)%pixel)))
     tempRealArray = 0.0
     tempRealArray2 = 0.0
     tempDoubleArray = 0.0_db
     tempDoubleArray2 = 0.0_db

     tempDoubleArray = reshape(obsImageSet(i)%pixel%i,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     obsImageSet(i)%pixel%i = reshape(tempDoubleArray2,SHAPE(obsImageSet(i)%pixel%i))

     tempDoubleArray = reshape(obsImageSet(i)%pixel%q,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     obsImageSet(i)%pixel%q = reshape(tempDoubleArray2,SHAPE(obsImageSet(i)%pixel%q))

     tempDoubleArray = reshape(obsImageSet(i)%pixel%u,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     obsImageSet(i)%pixel%u = reshape(tempDoubleArray2,SHAPE(obsImageSet(i)%pixel%u))

     tempDoubleArray = reshape(obsImageSet(i)%pixel%v,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     obsImageSet(i)%pixel%v = reshape(tempDoubleArray2,SHAPE(obsImageSet(i)%pixel%v))


     tempRealArray = reshape(obsImageSet(i)%vel,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     obsImageSet(i)%vel = reshape(tempRealArray2,SHAPE(obsImageSet(i)%vel))

     tempRealArray = reshape(obsImageSet(i)%totWeight,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     obsImageSet(i)%totWeight = reshape(tempRealArray2,SHAPE(obsImageSet(i)%totWeight))

     deallocate(tempRealArray)
     deallocate(tempRealArray2)
     deallocate(tempDoubleArray)
     deallocate(tempDoubleArray2)
   end do
     if (doRaman) then
        print *, 'MPI o6Image not implemented!'
        stop
     endif
  endif ! (stokesimage)

 if (doPvimage) then
   do iSlit = 1, nSlit
     allocate(tempDoubleArray(SIZE(pvImage(i)%pixel)))
     allocate(tempDoubleArray2(SIZE(pvImage(i)%pixel)))

     tempDoubleArray = reshape(pvImage(i)%pixel,(/SIZE(tempDoubleArray)/))
     tempDoubleArray2 = 0.0_db
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     pvImage(i)%pixel = reshape(tempDoubleArray2,SHAPE(pvImage(i)%pixel))

     deallocate(tempDoubleArray)
     deallocate(tempDoubleArray2)

     allocate(tempRealArray(SIZE(pvImage(i)%vAxis)))
     tempRealArray = 0.0
     call MPI_REDUCE(pvImage(i)%vAxis,tempRealArray,SIZE(pvImage(i)%vAxis),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     pvImage(i)%vAxis = tempRealArray
     deallocate(tempRealArray)

     allocate(tempRealArray(SIZE(pvImage(i)%pAxis)))
     tempRealArray = 0.0
     call MPI_REDUCE(pvImage(i)%pAxis,tempRealArray,SIZE(pvImage(i)%pAxis),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     pvImage(i)%pAxis = tempRealArray
     deallocate(tempRealArray)

   end do ! iSlit 
endif ! (doPvimage)
#endif

 if (myRankIsZero) then 

     write(message,*) " "
     call writeInfo(message, TRIVIAL)
     write(message,'(a)') "Model summary"
     call writeInfo(message, TRIVIAL)
     write(message,'(a)') "-------------"
     call writeInfo(message, TRIVIAL)
     write(message,*) " "
     call writeInfo(message, TRIVIAL)

     
     if (grid%adaptive) then
        write(message,*)  tooFewSamples, ' rays had 2 or less samples.'
        call writeInfo(message, TRIVIAL)
        write(message,*) BoundaryProbs, ' rays had numerical problems with octal boundaries.'
        call writeInfo(message, TRIVIAL)
        write(message,*) negativeOpacity, ' rays had problems with negative opacity values.'
        call writeInfo(message, TRIVIAL)
     end if
     write(message,*) " "
     call writeInfo(message, TRIVIAL)
     write(message,*) "Average # of scattering per photon:: ", real(nTot)/real(nPhotons)
     call writeInfo(message, TRIVIAL)
     write(message,*) " " 
     call writeInfo(message, TRIVIAL)

 end if 

!     if (.not.grid%cartesian.and.(grid%rCore /= 0.)) then
!        if (wtot_line /= 0.) write(*,*) "Mean radius of line formation",meanr_line/wtot_line/grid%rCore
!        if (wtot0_line /= 0.) write(*,*) "Mean radius of line zero",meanr0_line/wtot0_line/grid%rCore
!        if (wtot_cont /= 0.) write(*,*) "Mean radius of cont formation",meanr_cont/wtot_cont/grid%rCore
!        if (wtot0_cont /=0.) write(*,*) "Mean radius of cont zero",meanr0_cont/wtot0_cont/grid%rCore
!     endif

 varianceArray = STOKESVECTOR(0.d0, 0.d0, 0.d0, 0.d0)
! do i = 1, nLambda
!    do j = 1, nOuterLoop
!    varianceArray(i)%i = varianceArray(i)%i + (errorArray(j,i)%i - yArray(i)%i/dble(nOuterLoop))**2
!    varianceArray(i)%q = varianceArray(i)%q + (errorArray(j,i)%q - yArray(i)%q/dble(nOuterLoop))**2
!    varianceArray(i)%u = varianceArray(i)%u + (errorArray(j,i)%u - yArray(i)%u/dble(nOuterLoop))**2
!    varianceArray(i)%v = varianceArray(i)%v + (errorArray(j,i)%v - yArray(i)%v/dble(nOuterLoop))**2
! enddo
!enddo
 if (myRankIsZero) then 
    if (PRESENT(overrideFilename)) outfile = overrideFilename
    if (nLambda > 1) then
       if (nPhase == 1) then
          
          call writeSpectrum(outFile,  nLambda, grid%lamArray, yArray, varianceArray,&
               .false., objectDistance, .false., lamLine)
          
          specFile = trim(outfile)//"_stellar_direct"
          call writeSpectrum(specFile,  nLambda, grid%lamArray, yArrayStellarDirect, varianceArray,&
            .false., objectDistance, .false., lamLine)
          
          specFile = trim(outfile)//"_stellar_scattered"
          call writeSpectrum(specFile,  nLambda, grid%lamArray, yArrayStellarScattered, varianceArray,&
               .false., objectDistance, .false., lamLine)
          
          specFile = trim(outfile)//"_thermal_direct"
       call writeSpectrum(specFile,  nLambda, grid%lamArray, yArrayThermalDirect, varianceArray,&
            .false., objectDistance, .false., lamLine)

       specFile = trim(outfile)//"_thermal_scattered"
       call writeSpectrum(specFile,  nLambda, grid%lamArray, yArrayThermalScattered, varianceArray,&
            .false., objectDistance, .false., lamLine)
          
       
       if (velocitySpace) then
          specFile = trim(outfile)//"_v"
          call writeSpectrum(specFile,  nLambda, grid%lamArray, yArray, varianceArray,&
               .true., objectDistance, velocitySpace, lamLine)
       endif
       
    else
       write(tempChar,'(i3.3)') iPhase
       specFile = trim(outfile)//trim(tempChar)
       
        call writeSpectrum(specFile,  nLambda, grid%lamArray, yArray, varianceArray, &
             .false., objectDistance, velocitySpace, lamLine)
        
        if (velocitySpace) then
           tempChar = trim(specFile)//"_v"
           call writeSpectrum(tempChar,  nLambda, grid%lamArray, yArray, varianceArray,&
                .true., objectDistance, velocitySpace, lamLine)
        endif
     endif
        
        if (doRaman) then
           write(tempChar,'(i3.3)') iPhase
           o6filename = trim(outfile)//"_o6_"//trim(tempChar)//".dat"
           open(20,file=o6filename,status="unknown",form="formatted")
           do i = 1, no6pts
              t1 = 1./(o6xArray(2)-o6xArray(1))
              write(20,*) o6xarray(i),o6yarray(i)*t1,1.e-20,1.e-30,1.e-20,1.e-30
           enddo
           close(20)
        endif
     endif

     if (stokesimage) then
        do i1 = 1, nImageLocal
           name_filter = get_filter_name(filters, i1)
           bandwidth = 0.5*FWHM_filters(filters, i1)  ! 1/2 of FWHM  [A]
           lambda_eff = lambda_eff_filters(filters, i1) ! Effective wavelength of filter in [A]   
!           write(specFile,'(a,a,a,i3.3)') trim(outfile),"_"//trim(name_filter),"_image",iPhase
!           call writeImage(obsImageSet(i1), specfile, objectDistance, imageInArcsec, lambda_eff, bandwidth)
           write(specFile,'(a,a,a,i3.3,a)') trim(outfile),"_"//trim(name_filter),"_image",iPhase,".fits"
           if (torusVersion(2:2) == "2") specfile = originaloutfile


#ifdef USECFITSIO
           call writeFitsImage(obsImageSet(i1), trim(specfile), objectDistance, "intensity")
           if (polarizationImages) then
              header = specfile(1:index(specfile,".fits")-1)
              write(specFile,'(a,a)') trim(header)//"_pol.fits"
              call writeFitsImage(obsImageSet(i1), trim(specfile), objectDistance, "pol")
              write(specFile,'(a,a)') trim(header)//"_q.fits"
              call writeFitsImage(obsImageSet(i1), trim(specfile), objectDistance, "stokesq")
              write(specFile,'(a,a)') trim(header)//"_u.fits"
              call writeFitsImage(obsImageSet(i1), trim(specfile), objectDistance, "stokesu")
           endif
#endif

        end do
        if (doRaman) then
           write(specFile,'(a,a,i3.3)') trim(outfile),"_o6image",iPhase
           !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           ! Check the the effectieve wavelength and the 
           ! bandwith of O6 image here later and replace 5000 and 1.0d0 below!
           !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!           call writeImage(o6Image(1), specfile, objectDistance, imageInArcsec, 5000d0, 1.0d0)
        endif
        if (get_filter_set_name(filters) == "pn") then
           write(specFile,'(a,a,i3.3,a)') trim(outfile),"_image",iPhase,".ppm"
           call writeFalseColourPPM(trim(specfile), obsImageSet)
        endif

     endif

     if (doPvimage) then
        do iSlit = 1, nSlit
           write(specFile,'(a,a,i3.3,a,i2.2)') trim(outfile),"_pvimage",iPhase,"_slit_",iSlit

           call smoothPVimage(pvImage(iSlit), vfwhm/2.35, pfwhm/2.35)

        enddo
     endif
  end if ! (myRankIsZero)

     if (stokesImage) then
        do i = 1, nImageLocal
           call freeImage(obsImageSet(i))
        end do
     end if

     if (doPVimage) then
        do iSlit = 1, nSlit
           call freePVimage(pvImage(iSlit))
        enddo
     endif

777  continue
     end do incLoop ! end of multiple inclination loop
     
     if (associated(starSurface%element)) then
       call emptySurface(starSurface)
     end if

#ifdef MPI
! No need to free the grid if there is only one trip on phaseloop
 if (myRankGlobal /= 0 .and. .not.noPhaseUpdate .and. nStartPhase /= nEndPhase ) call freeGrid(grid)
#endif
  call torus_mpi_barrier !('waiting inside end of phase loop...')
  enddo phaseLoop

  deallocate(sourceSpectrum)
  deallocate(sourceSpectrum2)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine do_one_outer_photon_loop

    implicit none

    logical :: hitcore
    logical :: thrustar
    logical :: escaped, absorbed
    logical :: photonFromEnvelope

    type(VECTOR) :: rHatinStar
    integer :: nScat
    real :: directionalWeight
    real :: junk
    real :: fac1, fac2, fac3
    real(double) :: vRay, vOverCsqr
    real :: observedLambda
    real :: ramanWeight
    type(PHOTON) :: obsPhoton
    type(VECTOR) :: positionOc    ! photon position position

  ! photons
    type(PHOTON) :: outPhoton
    type(PHOTON) :: tempPhoton
    type(PHOTON) :: thisPhoton
    integer :: intPathError

    ! default inner loop indices
    iInner_beg = 1
    iInner_end = nInnerLoop

#ifdef MPI
  !====================================================================================
  ! Splitting the innerPhoton loop for multiple processors.
 ! if (myRankGlobal == 0) then
 !    print *, ' '
 !    print *, 'innerPhotonLoop computed by ', nThreadsGlobal, ' processors.'
 !    print *, ' '
 ! endif
  
    iInner_beg = myRankGlobal * (nInnerLoop/nThreadsGlobal) + 1
    iInner_end = (myrankGlobal+1) * (nInnerLoop/nThreadsGlobal) 
!  write(*,*) "rank ", myrankglobal, " doing ", iInner_beg, " to " , iinner_end


!  ! No need to use some processors if there are more processors
!  ! than the number of photons....
!  if (myRankGlobal > nPhotons/nOuterLoop - 1)  return
!    
!  if (myRankGlobal == 0) then
!     ! we will use an array to store the rank of the process
!     !   which will calculate each photon
!     allocate(photonBelongsRank(nPhotons/nOuterLoop))
!    
!     call mpiBlockHandout(nThreadsGlobal,photonBelongsRank,blockDivFactor=40,tag=tag,&
!                          setDebug=.false.)
!     deallocate(photonBelongsRank) ! we don't really need this here. 
!  end if
!  !====================================================================================
!
!    
!    
!  if (myRankGlobal /= 0) then
!    mpiBlockLoop: do  
!      call mpiGetBlock(myRankGlobal,iInner_beg, iInner_end,rankComplete,tag,setDebug=.false.)  
!      if (rankComplete) exit mpiBlockLoop  
#endif


!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(i, contPhoton, contWindPhoton, r, nScat) &
!$OMP PRIVATE(thisPhoton, directionalWeight) &
!$OMP PRIVATE(ilambda, ok) &
!$OMP PRIVATE(hitCore, junk, thisLam, j, obs_weight, thisVel) &
!$OMP PRIVATE(i1, i2, i3, t1, t2, t3, vray, vovercsqr, fac, observedLambda) &
!$OMP PRIVATE(t, rHat, islit, fac1, fac2, fac3, obsPhoton, r1, r2, thisTau) &
!$OMP PRIVATE(escaped, currentScat, absorbed, dlambda, thisChi, thisSca) &
!$OMP PRIVATE(albedo, tempPhoton, redRegion, thrustar, ramanWeight) &
!$OMP PRIVATE(outPhoton,intPathError, nTau, escProb, spotPhoton) &
!$OMP PRIVATE(rHatinStar, positionOc, linePhotonalbedo, dopShift, lineResAbs, tau_bnd) &
!$OMP PRIVATE(photonfromEnvelope, sourceOctal, nFromEnv) &
!$OMP PRIVATE(sourceSubcell, tau_tmp, exp_minus_tau) &
!$OMP PRIVATE(testPhoton, dtau, currentOctal, currentSubcell) &
!$OMP PRIVATE(tempOctal, tempsubcell, thistaudble, finaltau ) &
!$OMP SHARED(iLambdaPhoton, maxTau, nOuterLoop, pointSource, doIntensivePeelOff, nMuMie) &
!$OMP SHARED(grid, nContPhotons, nPhotons, lineEmission, lamLine, nLambda) &
!$OMP SHARED(weightLinePhoton, weightSource, flatSpec, secondSource, secondSourcePosition) &
!$OMP SHARED(ramanSourceVelocity, doRaman) &
!$OMP SHARED(weightContPhoton, outVec)&
!$OMP SHARED(opaqueCore, lamStart, lamEnd, thinLine, coolStarPosition) &
!$OMP SHARED(viewVec, o6xArray, rotationAxis, o6image, screened) &
!$OMP SHARED(stokesImage, obsImageSet, doPvimage) &
!$OMP SHARED(nSlit, pvimage, gridDistance, meanr0_line, wtot0_line) &
!$OMP SHARED(sourceSpectrum, sourceSpectrum2, meanr0_cont,wtot0_cont, mie, noScattering) &
!$OMP SHARED(miePhase, nSpot, chanceSpot, fSpot, gridUsesAMR) &
!$OMP SHARED(useInterp, photLine ) &
!$OMP SHARED(probDust, WeightDust, WeightPhoto, source, nsource) &
!$OMP SHARED(energyPerPhoton, filters, nUpper, nLower, nImageLocal) &
!$OMP SHARED(iInner_beg, iInner_end) &
!$OMP SHARED(curtains, starSurface, VoigtProf, nDustType, ttauri_disc, ttau_disc_on) &
!$OMP SHARED(forcedWavelength, usePhotonWavelength, thin_disc_on, forceFirstScat, fastIntegrate) &
!$OMP SHARED(o6yArray, yArray, yArrayStellarScattered, yArrayStellarDirect, yArrayThermalScattered, yArrayThermalDirect) &
!$OMP REDUCTION(+: ntot,tooFewSamples, boundaryProbs, negativeOpacity)


   call returnKappa(grid, grid%OctreeRoot, 1, reset_kappa=.true.)

    rhatinStar = VECTOR(0.d0, 0.d0, 0.d0)
    hitCore = .false.
    obsPhoton%lambda = 0.
    outPhoton%lambda = 0.
    intPathError = 0
    photonFromEnvelope = .true.
    nFromEnv = 0
!$OMP DO SCHEDULE(DYNAMIC,10)
        innerPhotonLoop: do i = iInner_beg, iInner_end

!           write(*,*) omp_get_thread_num(), i
#ifdef MPI
 !  if (MOD(i,nThreadsGlobal) /= myRankGlobal) cycle innerPhotonLoop
#endif
           ! The following six arrays must be allocated and deallocated for each 
           ! innerPhotonLoop to make the program work with OpenMP! (RK)

#ifdef _OPENMP
           allocate(lambda(1:maxTau))
           allocate(tauExt(1:maxTau))
           allocate(tauAbs(1:maxTau))
           allocate(tauSca(1:maxTau))
           allocate(linePhotonalbedo(1:maxTau))
           allocate(contTau(1:maxTau,1:nLambda)) 
           allocate(contWeightArray(1:nLambda))
#endif




#ifdef MPI
 !  if (MOD(i,nThreadsGlobal) /= myRankGlobal) goto 999
#endif

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
              call randomNumberGenerator(getReal=r)
              if (r < grid%chanceWindOverTotalContinuum) then
                 contWindPhoton = .true.
              endif
           endif

           nScat = 0

           ! initialize the photon

           
           contWeightArray = 1.

           select case(grid%geometry)
              case("planet")
                 call initPlanetPhoton(thisPhoton, grid, lamLine)
                 directionalWeight = 1.
                 contWindPhoton = .true.
              case DEFAULT
                 call initPhoton(thisPhoton, grid, nLambda, grid%lamArray, sourceSpectrum, &
                      lamLine, weightLinePhoton, &
                      weightContPhoton, contPhoton, flatspec, &
                      secondSource, secondSourcePosition, &
                      ramanSourceVelocity, contWindPhoton, directionalWeight, &
                      chanceSpot, fSpot, spotPhoton,  probDust, weightDust, weightPhoto, weightSource, &
                      source, nSource, rHatinStar, energyPerPhoton, filters, mie,&
                      starSurface, forcedWavelength, usePhotonWavelength, iLambdaPhoton,VoigtProf, &
                      photonfromEnvelope, dopShift=dopShift, sourceOctal=sourceOctal, sourcesubcell = sourceSubcell)
!                 write(*,*) "r, weight, i ", modulus(thisPhoton%position)/rcore,thisPhoton%weight, weightContPhoton, &
!                      Thisphoton%stokes%i
!                 if (.not.inOctal(sourceOctal, thisPhoton%position)) then
!                    write(*,*) "bug initializing photon"
!                 endif
                 if (thisPhoton%resonanceLine) then
                    r1 = real(i)/real(nPhotons/nOuterLoop)
                    thisPhoton%lambda = grid%lamArray(1) + r1*(grid%lamArray(nLambda)-grid%lamArray(1))
                 endif
           end select
           if (photonFromEnvelope) then
              nFromEnv = nFromEnv + 1
           endif

           if (thisPhoton%thermal.and.thisPhoton%stellar) then
              write(*,*) "error both thermal and stellar set"
              stop
           endif

           if (thisPhoton%scattered) then
              write(*,*) "error not direct produced", thisPhoton%scattered
              stop
           endif


           observedLambda = thisPhoton%lambda
!           if (thisPhoton%contPhoton) then
!
!              meanr_cont = meanr_cont + modulus(thisPhoton%position)*thisPhoton%stokes%i
!              wtot_cont = wtot_cont + thisPhoton%stokes%i
!           else
!              meanr_line = meanr_line + modulus(thisPhoton%position)*thisPhoton%stokes%i
!              wtot_line = wtot_line + thisPhoton%stokes%i
!           endif


           iLambda = findIlambda(thisPhoton%lambda, grid%lamArray, nLambda, ok)

           ! now we fire the photon direct to the observer
           lineResAbs = .false.  ! F for the zero-th scattering.



           if (doRaman) then
              
              call integratePath(gridUsesAMR, VoigtProf, &
                         thisPhoton%lambda, lamLine, &
                         thisPhoton%velocity, &
                         thisPhoton%position, outVec, grid, &
                         lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau , nTau, thin_disc_on, opaqueCore, &
                         escProb, thisPhoton%contPhoton,  &
                         nLambda, contTau, hitCore, thinLine, lineResAbs, .false., &
                         .false., nUpper, nLower, 0., 0., 0., junk,&
                         intPathError, &
                         useInterp, grid%Rstar1, coolStarPosition, nSource, source)
              if (intPathError == -10) then 
                 tooFewSamples = tooFewSamples + 1  
                 goto 999
              endif
              if (intPathError == -20) then 
                 boundaryProbs = boundaryProbs + 1
                 goto 999 
              endif
              if (intPathError == -70) then 
                 negativeOpacity = negativeOpacity + 1 
                 goto 999 
              endif

              thisLam = real(thisPhoton%lambda + (thisPhoton%velocity .dot. viewVec) * 1031.928)
              j = findIlambda(thisLam, o6xArray, no6pts, ok)
              if (ok) then
                 obs_weight = real(oneOnFourPi * exp(-tauExt(nTau)))
!$OMP ATOMIC
                 o6yArray(j) = o6yArray(j) + obs_weight
                 thisVel = (thisLam-lamLine)/lamLine

                 call addPhotonToImage(viewVec, rotationAxis,o6Image(1), 1, &
                      thisPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre)
              endif
           endif

           if (.not. screened) then


              ! find optical depths to observer
              if (.not.fastIntegrate) then
                 call integratePath(gridUsesAMR, VoigtProf, &
                      thisPhoton%lambda, lamLine, &
                      thisPhoton%velocity, &
                      thisPhoton%position, outVec, grid, &
                      lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, &
                      escProb, thisPhoton%contPhoton,  &
                      nLambda, contTau, hitCore, thinLine, lineResAbs, .false.,&
                      .false., nUpper, nLower, 0., 0., 0., junk,&
                      intPathError, &
                      useInterp, grid%Rstar1, coolStarPosition, nSource, source, &
                      startOctal=sourceOctal, startSubcell=sourceSubcell)
              else
                 
                 call tauAlongPathFast(ilambdaPhoton, grid, thisPhoton%position, outvec, finalTau, &
                      startOctal = sourceOctal, startSubcell=sourceSubcell , nTau=nTau, xArray=lambda, tauArray=tauExt)
              endif


!              if (thisPhoton%contPhoton.and.thisPhoton%stellar) &
!                   write(*,*) "Optical depth to observer: ",tauExt(ntau),thisPhoton%lambda, grid%lamArray(iLambdaPhoton), ntau, &
!                   modulus(thisPhoton%position)*1.d10/rsol, modulus(thisPhoton%position+(dble(lambda(ntau))*outVec))

!              if (thisPhoton%contPhoton)then
!                 octVec = thisPhoton%position
!                 CALL findSubcellTD(octVec,grid%octreeRoot,thisOctal,subcell)
!                 write(*,*) "Optical depth to observer: ",tauExt(ntau),1.d0/thisOctal%biascont3d(subcell), &
!                      thisPhoton%stokes%i
!              endif

               if (intPathError == -10) then
                    tooFewSamples = tooFewSamples + 1  
                    goto 999
               endif
               if (intPathError == -20) then
                    boundaryProbs = boundaryProbs + 1
                    goto 999
               endif
               if (intPathError == -70) then 
                  negativeOpacity = negativeOpacity + 1
                  goto 999
               endif


              if (.not. gridUsesAMR) call getIndices(grid,thisPhoton%position,i1,i2,i3,t1,t2,t3)


              if (thisPhoton%resonanceLine) escProb = 1.

              if (thisPhoton%linePhoton) then
                    
!
!                 if (tauExt(nTau) < 0.1)  then
!!                    do j = 1,ntau
!!                       write(*,*) lambda(j),tauExt(j)
!!                    enddo
!!                    write(*,*) " "
!                    tauExt(Ntau) = 1.e20
!                    tauAbs(nTau) = 1.e20
!                    tauSca(nTau) = 1.e-20
!                 endif

!!                 if (tauExt(ntau) < 10.) then
!!                    write(*,*) sqrt(thisPhoton%position%x**2 + thisPhoton%position%y**2),thisPhoton%position%z,dopshift
!                    write(99,*) sqrt(thisPhoton%position%x**2 + thisPhoton%position%y**2),thisPhoton%position%z,dopshift
!                    if ((sqrt(thisPhoton%position%x**2 + thisPhoton%position%y**2)<grid%rCore).and. &
!                         (thisPhoton%position%z < 0.)) then
!                       write(*,*) "SCREW UP"
!                    endif
!                 endif
!                 if (tauExt(nTau) < 1.) escProb = 0.
                 vray = -(thisPhoton%velocity .dot. outVec)
                 vovercsqr = modulus(thisPhoton%velocity)**2
                 fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                 observedlambda = real(thisPhoton%lambda / fac)
                 tau_tmp = tauExt(nTau)
                 exp_minus_tau = EXP(-tau_tmp)
                 obs_weight = real(oneOnFourPi * exp_minus_tau * escProb)
              endif

              if (thisPhoton%contPhoton.and.(.not.contWindPhoton)) then
                 t = modulus(thisPhoton%position)
                 if (t /= 0.) then
                    rHat = thisPhoton%position / t
                 else
                    rHat = thisPhoton%direction
                 endif
                 t = rHat .dot. outVec

!                 if (t < 0.0) then 
                    ! The photon directs inward of the star...
                    ! this should be taken care by integratepth routine
                    ! by producing a large tau but just in case this was not taken care.
                    ! -- (RK)
!                    obs_weight = 1.0e-30
!                 else
                    if (lineEmission) then
                       obs_weight = real(abs(t)*exp(-tauExt(nTau)) / pi)
                    else
                       obs_weight = real(oneOnFourPi*abs(t)*exp(-tauExt(nTau)))
                    endif
!                 end if
                 observedlambda = thisPhoton%lambda
              endif

              if (.not.flatSpec.and.(.not.hitCore)) then
                 if (contWindPhoton) then
                    obs_weight = real(oneOnFourPi*exp(-tauExt(nTau)))
                 else 
                    if (.not.pointSource) then
                       obs_weight = real((outVec.dot.rHatinStar)*exp(-tauExt(nTau))/pi)
                    else
                       obs_weight = real(oneOnFourPi * exp(-tauExt(nTau)))
                    endif
                    if (obs_weight < 0.) obs_weight = 0.
                 endif
                 
                 iLambda = findIlambda(observedlambda, grid%lamArray,  nLambda, ok)
                 if (ok) then
!$OMP CRITICAL ( updateYarray )
!                    write(*,*) "addingPhoton",thisPhoton%stokes%i,obs_weight
                    yArray(iLambda) = yArray(iLambda) + thisPhoton%stokes * obs_weight

                    if (thisPhoton%stellar) then
                       if (thisPhoton%scattered) then
                          yArrayStellarScattered(iLambda) = yArrayStellarScattered(iLambda) + (thisPhoton%stokes * obs_weight)
                       else
                          yArrayStellarDirect(iLambda) = yArrayStellarDirect(iLambda) + (thisPhoton%stokes * obs_weight)
                       endif
                    endif

                    if (thisPhoton%thermal) then
                       if (thisPhoton%scattered) then
                          yArrayThermalScattered(iLambda) = yArrayThermalScattered(iLambda) + (thisPhoton%stokes * obs_weight)
                       else
                          yArrayThermalDirect(iLambda) = yArrayThermalDirect(iLambda) + (thisPhoton%stokes * obs_weight)
                       endif
                    endif
!$OMP END CRITICAL ( updateYarray )
                 endif
                 if (stokesImage) then
                    thisVel = 0. ! no velocity for dust continuum
                    call addPhotonToImage(viewVec, rotationAxis, obsImageSet, nImageLocal, thisPhoton,&
                         thisVel, obs_weight, filters, grid%octreeRoot%centre)

                 endif
                 if (doPVimage) then
                    do iSlit = 1, nSlit
                       call addPhotontoPVimage(pvImage(iSlit), thisPhoton, viewVec, rotationAxis, thisVel, &
                                 obs_weight, gridDistance)
                    enddo
                 endif
              endif
              
              if (flatSpec.and.(.not.hitCore)) then
                 if (thisPhoton%linePhoton) then
                    fac2 = 1.
!                    if (thisPHoton%resonanceline) fac2 = directionalWeight
                    iLambda = findIlambda(observedlambda, grid%lamArray, nLambda, ok)
                    if (ok) then
!$OMP CRITICAL ( updateYarray )
                       yArray(iLambda) = yArray(iLambda) + &
                            (thisPhoton%stokes * (fac2 * oneOnFourPi * escProb * exp(-tauExt(nTau))))
!$OMP END CRITICAL ( updateYarray )
                    endif

                    obs_weight = real(oneOnFourPi * escProb * exp(-tauExt(nTau))*fac2)

                    meanr0_line = real(meanr0_line + modulus(thisPhoton%position) * &
                         (thisPhoton%stokes%i * obs_weight))
                    wtot0_line = real(wtot0_line + (thisPhoton%stokes%i * obs_weight))
                    

                    thisVel = (observedLambda-lamLine)/lamLine
                    if (stokesImage) then
                       call addPhotonToImage(viewVec, rotationAxis, obsImageSet, nImageLocal, &
                            thisPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre)
                    endif
                    if (dopvImage) then
                       do iSlit = 1, nSlit
                          call addPhotontoPVimage(pvImage(iSlit), thisPhoton, viewVec,  rotationAxis,thisVel, &
                            obs_weight, gridDistance)
                       enddo
                    endif

                 else  ! contunuum photon
                    fac1 = 1.
                    if (.not.contWindPhoton) then
                       fac1 = real(2.*abs(thisPhoton%originalNormal.dot.outVec)/twoPi)
                    else
                       fac1 = oneOnfourPi
                    endif
                    
                    vray = -(thisPhoton%velocity .dot. outVec)
                    vovercsqr = modulus(thisPhoton%velocity)**2
                    fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                    observedlambda = real(thisPhoton%lambda / fac)
                    
                    i1 = 0
!
                    do iLambda = 1, nLambda

                       thisLam = (lamLine-observedlambda) + grid%lamArray(iLambda)


                       call hunt(grid%lamArray,nLambda,thisLam,i1)
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
                       obs_weight = (fac1 * exp(-(tauExt(ntau)+fac3)))*fac2

!$OMP CRITICAL ( updateYarray )
                       yArray(iLambda) = yArray(iLambda) + (thisPhoton%stokes * obs_weight)

                       if (thisPhoton%stellar) then
                          if (thisPhoton%scattered) then
                             yArrayStellarScattered(iLambda) = yArrayStellarScattered(iLambda) + (thisPhoton%stokes * obs_weight)
                          else
                             yArrayStellarDirect(iLambda) = yArrayStellarDirect(iLambda) + (thisPhoton%stokes * obs_weight)
                          endif
                       endif
                       
                       if (thisPhoton%thermal) then
                          if (thisPhoton%scattered) then
                             yArrayThermalScattered(iLambda) = yArrayThermalScattered(iLambda) + (thisPhoton%stokes * obs_weight)
                          else
                             yArrayThermalDirect(iLambda) = yArrayThermalDirect(iLambda) + (thisPhoton%stokes * obs_weight)
                          endif
                       endif
!$OMP END CRITICAL ( updateYarray )

                       meanr0_cont = real(meanr0_cont + modulus(thisPhoton%position) &
                            * thisPhoton%stokes%i*obs_weight)
                       wtot0_cont = real(wtot0_cont + thisPhoton%stokes%i*obs_weight)
                       thisVel = (observedLambda-lamLine)/lamLine
                       if (stokesImage) then
                          call addPhotonToImage(viewVec,  rotationAxis,obsImageSet, nImageLocal, &
                               thisPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre, &
                               grid%lamArray(iLambda))
                       endif
                       if (dopvImage) then
                          do iSlit = 1, nSlit
                             call addPhotontoPVimage(pvImage(islit), obsPhoton, viewVec,  rotationAxis, thisVel, &
                               obs_weight, gridDistance)
                          enddo
                       endif

                    enddo
                 endif
              endif

           endif



           if (.not.fastIntegrate) then
              call integratePath(gridUsesAMR, VoigtProf, &
                   thisPhoton%lambda, lamLine, &
                   thisPhoton%velocity, &
                   thisPhoton%position, &
                   thisPhoton%direction, grid, &
                   lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, &
                   escProb, thisPhoton%contPhoton, &
                   nLambda, contTau, hitCore, thinLine, lineResAbs, .false.,  &
                   .false., nUpper, nLower, 0., 0., 0., &
                   junk,intPathError, &
                   useInterp, grid%Rstar1, coolStarPosition, nSource, source, startOctal=sourceOctal, startSubcell=sourceSubcell)
           else
              call tauAlongPathFast(ilambdaPhoton, grid, thisPhoton%position, thisPhoton%direction, finaltau,&
                   startOctal = sourceOctal, startSubcell=sourceSubcell , nTau=nTau, xArray=lambda, tauArray=tauExt)
           endif


              if (intPathError == -10) then 
                 tooFewSamples = tooFewSamples + 1  
                 goto 999 
              endif
              if (intPathError == -20) then 
                 boundaryProbs = boundaryProbs + 1
                 goto 999 
              endif
              if (intPathError == -70) then 
                 negativeOpacity = negativeOpacity + 1
                 goto 999 
              endif

            if (forceFirstScat) then ! this option is off by default.
               ! Useful for optically thin case (to get correct scatterng flux).
               tau_bnd = tauExt(nTau) ! Optical depth to outer boundary.
               fac=1.0d0-exp(-tau_bnd)
               call randomNumberGenerator(getReal=r1)
               thisTau = real(min(-log(1.d0-r1*fac),tau_bnd)  )
               escaped = .false.!
               ! This is done so to force the first scattering if 
               ! MAXSCAT>0. 
            else
               call randomNumberGenerator(getReal=r1)
               fac = 1.
               thistau = -log(max(1.e-20,(1. - r1)))
               if (thistau .gt. tauExt(nTau)) then
                  escaped = .true.  
               else
                  escaped = .false.  
               endif
            end if



           if (thisPhoton%contPhoton) escProb = 1.

           if (thisPhoton%resonanceLine) escProb = 1.
              
           ! Correcting the weight here... 
           ! escProb -- allowing for local emission
           thisPhoton%stokes = thisPhoton%stokes * real(fac) * escProb

!           write(*,*) "after first scat ", thisPhoton%stokes%i
           currentScat = 0


!           escaped = .false.
           absorbed = .false.
           if (doRaman) hitCore = .false.


           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
           ! scattering loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
           lineResAbs = .true.
           do while ((.not.escaped) .and. (.not.absorbed) .and. &
                     (nScat < maxScat)   )
              ! Now we do not force the scattering exept for the first scattering.
              ! maxScat specifies the max number of scattering ALLOWED, but 
              ! not the number of scattering FORCED!
              

              nScat = nScat + 1

!              write(*,*) i,nscat,thisPhoton%stokes%i,thisPhoton%lambda, modulus(thisPhoton%position)*1.e10/drinner

              currentScat = currentScat  + 1

              ! find position of next interaction

              call locate(tauExt, nTau, thisTau, j)


              if (j < nTau) then
                 t = 0.
                 if ((tauExt(j+1) - tauExt(j)) /= 0.) then
                    t = (thisTau - tauExt(j)) / (tauExt(j+1) - tauExt(j))
                 endif
                 dlambda = real(lambda(j) + (lambda(j+1)-lambda(j))*t)
              else
                 dlambda = lambda(nTau)
              endif




! Here I have added a new algorithm to increase the S/N of scattered light images. We take the tausca and tauabs
! arrays and perform "peel offs" towards the observer at every index along the ray. This vastly increases
! the computational effort, but should hugely increase the S/N of the scattered light image

              if (doIntensivePeelOff) then
                 testPhoton = thisPhoton
                 redRegion = .false.
                 iLambda = findIlambda(testPhoton%lambda, grid%lamArray, nLambda, ok)

                 CALL findSubcellTD(testPhoton%position,grid%octreeRoot,tempOctal,tempsubcell)

                 do iStep = 2, ntau

                    testPhoton%position = thisPhoton%position + real(lambda(iStep),kind=oct)*thisPhoton%direction



                    if (tauExt(iStep) > 10.) exit
                    call findSubcellLocal(testPhoton%position, tempOctal, tempsubcell)


                    fac = exp(-tauExt(iStep))
                    dtau = tauExt(iStep) - tauExt(iStep-1)
                    fac= fac * (1.0d0-exp(-dtau)) * albedo

                    call scatterPhoton(grid, testPhoton, outVec, obsPhoton, mie, &
                         miePhase, nDustType, nLambda, grid%lamArray, nMuMie, ttau_disc_on, ttauri_disc, &
                         currentOctal=tempOctal, currentSubcell=tempsubcell)
                    

                    call tauAlongPath2(ilambda, grid, obsPhoton%position, obsPhoton%direction, thistaudble, tauMax=20.d0, &
                         startOctal=tempOctal, startSubcell=tempSubcell)

                    call amrGridValues(grid%octreeRoot, testPhoton%position, grid=grid, startOctal=tempOctal, &
                         actualSubcell=tempSubcell, iLambda=iLambda, &
                         kappaAbs = thisChi, kappaSca = thisSca)

                    albedo = thisSca / (thisSca + thisChi)

                    if (noScattering) albedo = 0.

                    obs_weight = real(fac*oneOnFourPi*exp(-thisTaudble))


!                    write(*,*) iStep, lambda(istep),obs_weight,tauExt(istep),tauExtObs(nTauObs)
                    if (stokesImage) then
                       thisVel = 0. ! no velocity for dust continuum emission
                       call addPhotonToImage(viewVec,  rotationAxis, obsImageSet, nImageLocal,  &
                            obsPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre)
                    endif
                    if (dopvImage) then
                       do iSlit = 1, nSlit
                          call addPhotontoPVimage(pvImage(iSlit), obsPhoton, viewVec,  rotationAxis,thisVel, &
                               obs_weight, gridDistance)
                       enddo
                    endif
                 enddo
              endif



                    
              ! New photon position 
              thisPhoton%position = thisPhoton%position + real(dlambda,kind=oct)*thisPhoton%direction
              ! adjusting the photon weights 
              if ((.not. mie) .and. (.not. thisPhoton%linePhoton)) then
                 if (j < nTau) then
                    contWeightArray(1:nLambda) = real(contWeightArray(1:nLambda) *  &
                         EXP(-(contTau(j,1:nLambda) + t*(contTau(j+1,1:nLambda)-contTau(j,1:nLambda))) ))
                 else
                    contWeightArray(1:nLambda) = contWeightArray(1:nLambda)*EXP(-(contTau(nTau,1:nLambda)))
                 end if
              end if


              if (flatspec) then
                 iLambda = 1
              else
                 iLambda = findIlambda(thisPhoton%lambda, grid%lamArray, nLambda, ok)
              endif

              if (grid%adaptive) then
                 positionOc = thisPhoton%position
                 call amrGridValues(grid%octreeRoot, positionOc, grid=grid, iLambda=iLambda, &
                      kappaAbs = thisChi, kappaSca = thisSca, foundOctal=currentOctal, foundSubcell=currentSubcell)
              else
                 call getIndices(grid, thisPhoton%position, i1, i2, i3, t1, t2, t3)
                 if (.not.grid%oneKappa) then
                    if (.not.flatspec) then
                       thisChi = interpGridKappaAbs(grid, i1, i2, i3, iLambda, t1, t2, t3) 
                       thisSca = interpGridKappaSca(grid, i1, i2, i3, iLambda, t1, t2, t3)
                    else
                       thisChi = interpGridKappaAbs(grid, i1, i2, i3, 1, t1, t2, t3)
                       thisSca = interpGridKappaSca(grid, i1, i2, i3, 1, t1, t2, t3)
                    endif
                 else
                    r = interpGridScalar2(grid%rho,grid%na1,grid%na2,grid%na3,i1,i2,i3,t1,t2,t3)
                    thisChi = grid%oneKappaAbs(1,iLambda) * r
                    thisSca = grid%oneKappaSca(1,iLambda) * r
                 endif
              end if
              
              if (contPhoton) then
                 if ((thisChi+thisSca) >= 0.) then
                    albedo = thisSca / (thisChi + thisSca)
                 else
                    albedo = 0.
                    write(*,*) "Error:: thisChi+thisSca < 0 in torusMain."
                    !                 stop
                 endif
              else
                 albedo = linePhotonAlbedo(j)
              endif
              if (noScattering) albedo = 0.

              if (grid%resonanceLine) albedo = 1.

              if (hitCore) then
                 if (grid%geometry /= "binary") then

                    if (grid%rCore /= 0.) then
                       r = real(modulus(thisPhoton%position)/grid%rCore)
                    else 
                       r = 10.
                    endif
                    if (r < 1.01) then
                       absorbed = .true.
                       albedo = 0.
                    endif
                 else
                    r1 = real(modulus(thisPhoton%position - grid%starpos1)/dble(grid%rStar1))
                    r2 = real(modulus(thisPhoton%position - grid%starpos2)/dble(grid%rStar2))
                    r = min (r1, r2)
                    if (r < 1.01) then
                       albedo = 0.5
                    endif
                 endif
              endif
              if (noScattering) albedo = 0.


              call randomNumberGenerator(getReal=r)
              if (r > albedo) then
                 absorbed = .true.
              endif


!              thisPhoton%stokes = thisPhoton%stokes * albedo  ! weight adjusted here!!!

              if (thisPhoton%stokes%i < reallySmall) then
                 absorbed = .true.
!                 write(*,*) "! Small photon weight",thisPhoton%stokes%i,thisPhoton%lambda,albedo,hitcore
              endif




              ! towards observer

              if (.not.absorbed) then

                 tempPhoton = thisPhoton
                 
                 call scatterPhoton(grid, tempPhoton, outVec, obsPhoton, mie, &
                       miePhase, nDustType, nLambda, grid%lamArray, nMuMie, ttau_disc_on, ttauri_disc, &
                       currentOctal, currentSubcell)

                 ! the o6 photon might get scattered towards the observer by a rayleigh scattering

                 if (doRaman) then
                    call integratePath(gridUsesAMR, VoigtProf, &
                            obsPhoton%lambda, lamLine, obsPhoton%velocity, &
                            obsPhoton%position, obsPhoton%direction, grid, &
                            lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau,  thin_disc_on, opaqueCore, &
                            escProb, obsPhoton%contPhoton, nLambda, contTau, hitCore, &
                            thinLine, lineResAbs, .false., .false., nUpper, nLower, 0., 0., 0.,&
                            junk, intPathError, &
                            useinterp, grid%Rstar1, coolStarPosition, nSource, source)


                   if (intPathError == -10) then 
                      tooFewSamples = tooFewSamples + 1  
                      goto 999 
                   endif
                   if (intPathError == -20) then 
                      boundaryProbs = boundaryProbs + 1
                      goto 999 
                   endif
                   if (intPathError == -70) then 
                      negativeOpacity = negativeOpacity + 1
                      goto 999 
                   endif

                    vray = -(obsPhoton%velocity .dot. outVec)
                    vovercsqr = modulus(thisPhoton%velocity)**2
                    fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                    observedlambda = real(obsPhoton%lambda / fac)
                    
                    obs_weight = real(oneOnFourPi*exp(-tauExt(nTau)) * 34./(6.6+34.))

                    j = findIlambda(observedLambda, o6xArray, no6pts, ok)
                    if (ok) then
!! !$OMP ATOMIC
                       o6yArray(j) = o6yArray(j) + obs_weight
                    endif

                 endif  ! raman




                 redRegion = .false.

                 if (doRaman) redRegion = .true.

                 if (.not.fastIntegrate) then
                    call integratePath(gridUsesAMR, VoigtProf, &
                         obsPhoton%lambda, lamLine, &
                         obsPhoton%velocity, &
                         obsPhoton%position, obsPhoton%direction, grid, &
                         lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau,  thin_disc_on, opaqueCore, &
                         escProb, obsPhoton%contPhoton, &
                         nLambda, contTau, hitCore, &
                         thinLine, lineResAbs, redRegion, &
                         .false., nUpper, nLower, 0., 0.,0.,junk, intPathError, &
                         useInterp, grid%Rstar1, coolStarPosition, nSource, source, &
                         startOctal=currentOctal, startSubcell=currentSubcell)                 
                 else
                    call tauAlongPathFast(ilambdaPhoton, grid,obsPhoton%position, obsPhoton%direction, finalTau,&
                         startOctal = currentOctal, startSubcell=currentSubcell , nTau=nTau, xArray=lambda, tauArray=tauExt)
                 endif

                 if (intPathError == -10) then 
                    tooFewSamples = tooFewSamples + 1  
                    goto 999 
                 endif
                 if (intPathError == -20) then 
                    boundaryProbs = boundaryProbs + 1
                    goto 999 
                 endif
                 if (intPathError == -70) then 
                    negativeOpacity = negativeOpacity + 1
                    goto 999 
                 endif



                 if (obsPhoton%linePhoton) then

                    
                    vray = -(obsPhoton%velocity .dot. outVec)
                    vovercsqr = modulus(obsPhoton%velocity)**2
                    fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                    if (.not.doRaman) then
                       observedlambda = real(obsPhoton%lambda / fac)
                    else
                       observedlambda = real(obsPhoton%redlambda / fac)
                    endif


                 else
                    observedLambda = obsPhoton%lambda
                 endif

                 thruStar = .false.


                 if (.not.flatSpec.and.(.not.hitCore).and.(.not.thruStar)) then

                    ramanWeight = 1.
                    if (doRaman) ramanWeight = 6.6 / (34.+6.6)
                    obs_weight = real(oneOnFourPi*exp(-tauExt(nTau)) * ramanWeight)


                    iLambda = findIlambda(observedlambda, grid%lamArray, nLambda, ok)
                      

                    if (ok) then
!$OMP CRITICAL ( updateYarray )
                       yArray(iLambda) = yArray(iLambda) + obsPhoton%stokes*obs_weight


                       if (obsPhoton%stellar) then
                          if (obsPhoton%scattered) then
                             yArrayStellarScattered(iLambda) = yArrayStellarScattered(iLambda) + (obsPhoton%stokes * obs_weight)
                          else
                             yArrayStellarDirect(iLambda) = yArrayStellarDirect(iLambda) + (obsPhoton%stokes * obs_weight)
                          endif
                       endif

                       if (obsPhoton%thermal) then
                          if (obsPhoton%scattered) then
                             yArrayThermalScattered(iLambda) = yArrayThermalScattered(iLambda) + (obsPhoton%stokes * obs_weight)
                             call randomNumberGenerator(getReal=r)
!                             if (obsPhoton%stokes%i*obs_weight > 1.e6) then
!                             write(*,'(i5,6f15.4)') myrankglobal,obsPhoton%stokes%i*obs_weight,obsPhoton%stokes%i, &
!                                  obs_weight,tauext(ntau), sqrt(obsPhoton%position%x**2+obsPhoton%position%y**2)/rinner,r
!                             endif
                          else
                             yArrayThermalDirect(iLambda) = yArrayThermalDirect(iLambda) + (obsPhoton%stokes * obs_weight)
                          endif
                       endif
!                       write(*,*) obsphoton%lambda,obsPhoton%stokes%i, obs_weight, nscat
!$OMP END CRITICAL ( updateYarray )
                    endif
                    if (doRaman) then
                       thisVel = (1./(1./observedLambda  + 1./1215.67))/1031.928-1.
                    else
                       thisVel = observedLambda
                    endif



                    if (stokesImage) then
                       thisVel = 0. ! no velocity for dust continuum emission
                       call addPhotonToImage(viewVec,  rotationAxis, obsImageSet, nImageLocal,  &
                            obsPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre)
                    endif
                    if (dopvImage) then
                       do iSlit = 1, nSlit
                          call addPhotontoPVimage(pvImage(iSlit), obsPhoton, viewVec,  rotationAxis,thisVel, &
                            obs_weight, gridDistance)
                       enddo
                    endif

                    if (doRaman) then
                       thisPhoton%stokes = thisPhoton%stokes*(1.-ramanWeight)
                    endif
                 endif



                 if (flatSpec.and.(.not.hitCore)) then



                    if (obsPhoton%linePhoton) then
                       iLambda = findIlambda(observedlambda, grid%lamArray, nLambda, ok)
                       obs_weight = real(oneOnFourPi*exp(-tauExt(nTau)))


                       if (obsPhoton%resonanceLine) obs_weight = obs_weight * directionalWeight

                       if (ok) then
!$OMP CRITICAL ( updateYarray )
                          yArray(iLambda) = yArray(iLambda) + obsPhoton%stokes*obs_weight
!$OMP END CRITICAL ( updateYarray )
                       endif

                       thisVel = observedLambda
                       thisVel = (observedLambda-lamLine)/lamLine
                       if (stokesImage) then
                          call addPhotonToImage(viewVec,  rotationAxis, obsImageSet, nImageLocal, &
                               obsPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre)
                       endif
                       if (dopvImage) then
                          do iSlit = 1 , nSlit
                             call addPhotontoPVimage(pvImage(iSlit), obsPhoton, viewVec,  rotationAxis, thisVel, &
                               obs_weight, gridDistance)
                          enddo
                       endif


                    else


                       vray = -(obsPhoton%velocity .dot. outVec)
                       vovercsqr = modulus(obsPhoton%velocity)**2
                       fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                       observedlambda = real(obsPhoton%lambda / fac)

                       do iLambda = 1,nLambda

                          thisLam = (lamLine-observedlambda) + grid%lamArray(iLambda)

                          i1 = findILambda(thisLam, grid%lamArray, nLambda, ok)
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


                          obs_weight = real(exp(-(tauExt(nTau)+fac3)) * oneOnFourpi & 
                          * fac2 * directionalWeight * contWeightArray(i1))

!$OMP CRITICAL ( updateYarray )
                          if (ok) then
                             yArray(iLambda) = yArray(iLambda) + &
                                  (obsPhoton%stokes*obs_weight)

                             if (obsPhoton%stellar) then
                                if (obsPhoton%scattered) then
                                   yArrayStellarScattered(iLambda) = yArrayStellarScattered(iLambda) + &
                                        (obsPhoton%stokes * obs_weight)
                                else
                                   yArrayStellarDirect(iLambda) = yArrayStellarDirect(iLambda) + &
                                        (obsPhoton%stokes * obs_weight)
                                endif
                             endif

                             if (obsPhoton%thermal) then
                                if (obsPhoton%scattered) then
                                   yArrayThermalScattered(iLambda) = yArrayThermalScattered(iLambda) + &
                                        (obsPhoton%stokes * obs_weight)
                                else
                                   yArrayThermalDirect(iLambda) = yArrayThermalDirect(iLambda) + &
                                        (obsPhoton%stokes * obs_weight)
                                endif
                             endif
                          else 
                             yarray(ilambda) = yArray(ilambda) + &
                                  obsPhoton%stokes*(oneOnFourpi* directionalweight * exp(-tauExt(ntau)))
                             if (obsPhoton%stellar) then
                                if (obsPhoton%scattered) then
                                   yArrayStellarScattered(iLambda) = yArrayStellarScattered(iLambda) + &
                                        (obsPhoton%stokes * (oneOnFourpi* directionalweight * exp(-tauExt(ntau))))
                                else
                                   yArrayStellarDirect(iLambda) = yArrayStellarDirect(iLambda) + &
                                        (obsPhoton%stokes * (oneOnFourpi* directionalweight * exp(-tauExt(ntau))))
                                endif
                             endif
                             
                             if (obsPhoton%thermal) then
                                if (obsPhoton%scattered) then
                                   yArrayThermalScattered(iLambda) = yArrayThermalScattered(iLambda) + &
                                        (obsPhoton%stokes * (oneOnFourpi* directionalweight * exp(-tauExt(ntau))))
                                else
                                   yArrayThermalDirect(iLambda) = yArrayThermalDirect(iLambda) + &
                                        (obsPhoton%stokes * (oneOnFourpi* directionalweight * exp(-tauExt(ntau))))
                                endif
                             endif
                          endif
!$OMP END CRITICAL ( updateYarray )


                          thisVel = (observedlambda-lamLine)/lamLine
                          if (stokesImage) then
                             call addPhotonToImage(viewVec,  rotationAxis, obsImageSet, nImageLocal, &
                                  obsPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre, &
                                  grid%lamArray(iLambda))
                          endif
                          if (dopvImage) then
                             do iSlit = 1, nSlit
                                call addPhotontoPVimage(pvImage(iSlit), obsPhoton, viewVec,  rotationAxis, thisVel, &
                                  obs_weight, gridDistance)
                             enddo
                          endif


                       enddo
                    endif
                 endif
                 
                 call scatterPhoton(grid,thisPhoton, zeroVector, outPhoton, mie, &
                       miePhase, nDustType, nLambda, grid%lamArray, nMuMie, ttau_disc_on, ttauri_disc, &
                       currentOctal, currentSubcell)
                 thisPhoton = outPhoton

                 if (.not.fastIntegrate) then
                 call integratePath(gridUsesAMR, VoigtProf, &
                            thisPhoton%lambda, lamLine, &
                            thisPhoton%velocity, thisPhoton%position, &
                            thisPhoton%direction, grid, lambda, tauExt, tauAbs, &
                            tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, &
                            escProb, thisPhoton%contPhoton, &
                            nLambda, contTau, hitCore, thinLine, lineResAbs, .false., &
                            .false., nUpper, nLower, 0.,&
                            0., 0., junk, intPathError, &
                            useInterp, grid%Rstar1, coolStarPosition, nSource, source, &
                            startOctal=currentOctal, startSubcell=currentSubcell)
                 else
                    call tauAlongPathFast(ilambdaPhoton, grid, thisPhoton%position, thisPhoton%direction, finalTau,&
                         startOctal = currentOctal, startSubcell=currentSubcell , nTau=nTau, xArray=lambda, tauArray=tauExt)
                 endif


                 if (intPathError == -10) then 
                    tooFewSamples = tooFewSamples + 1  
                    goto 999 
                 endif
                 if (intPathError == -20) then 
                    boundaryProbs = boundaryProbs + 1
                    goto 999 
                 endif
                 if (intPathError == -70) then 
                    negativeOpacity = negativeOpacity + 1
                    goto 999 
                 endif


                  call randomNumberGenerator(getReal=r1)
                  thistau = -log(max(1.e-20,(1. - r1)))
                  if (thistau .gt. tauExt(nTau)) then
                     escaped = .true.  
                  else
                     escaped = .false.  
                  endif
!                  write(*,*) myrankGlobal," Second scatter tauext, thisTau, escaped ",tauExt(ntau), thisTau, escaped


!                  if (maxScat == 1) then
!                     call randomNumberGenerator(getDouble=r1)
!                     thistau = -log(max(1.e-10,(1. - r1)))
!                     if (thistau .gt. tauExt(nTau)) then
!                        ! one scattering forced already, 
!                        ! so let it escape the second time if it must.
!                        escaped = .true.  
!                     else
!                        nscat = 0
!                     endif
!                  else  
!                     ! This should be never reached if MAXSCAT<2.
!                     tau_bnd = tauExt(nTau) ! Optical depth to outer boundary.
!                     if(abs(tau_bnd) < 0.01)then  ! appoximation
!                        fac=tau_bnd*(  1.0d0-0.5d0*tau_bnd*( 1.0d0-tau_bnd/3.0d0*  &
!                             (1.0d0-0.25d0*tau_bnd) )  )
!                     else
!                        fac=1.0d0-exp(-tau_bnd)
!                     end if
!                     call randomNumberGenerator(getDouble=r1)
!                     thisTau = -log(1.d0-r1*fac)
!                     ! (RK) commented out the next line.
!                     ! if (thisTau > tauExt(nTau)) print *, 'Fixed thisTau in scattering loop.'
!                     thisTau = min(dble(thisTau),tau_bnd)
!                     thisPhoton%stokes = thisPhoton%stokes * real(fac)
!                  endif


               endif  ! (.not. absorbed)

               
            enddo
            
            if ((nscat >= maxScat).and.(maxScat > 1)) then
               write(*,*) "Photon aborted after ",maxScat, " scatterings"
            endif
           nTot = nTot + nScat

999  continue  ! escape route for a bad photon

#ifdef _OPENMP
           deallocate(lambda)
           deallocate(tauSca)
           deallocate(tauExt)
           deallocate(tauAbs)
           deallocate(linePhotonalbedo)
           deallocate(contTau)
           deallocate(contWeightArray)
#endif
           
!           call plotspec(grid%lamArray, yArray, nLambda)

        enddo innerPhotonLoop



!        write(*,*) "Fraction of photons from envelope: ", real(nFromEnv)/real(nInnerLoop)

!$OMP END DO
!$OMP END PARALLEL

#ifdef MPI
!     write (*,'(A,I3,A,I3,A,I3,A)') 'Process ',myRankGlobal, &
!                      ' waiting to sync spectra... (',iOuterLoop,'/',nOuterLoop,')' 
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

   ! we have to syncronize the 'yArray' after each inner photon loop to
   !   get the statistics right 
     allocate(tempDoubleArray(SIZE(yArray)))
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArray%i,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArray%i = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArray%q,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArray%q = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArray%u,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArray%u = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArray%v,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArray%v = tempDoubleArray 
     deallocate(tempDoubleArray)

     allocate(tempDoubleArray(SIZE(yArrayStellarDirect)))
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayStellarDirect%i,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayStellarDirect%i = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayStellarDirect%q,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayStellarDirect%q = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayStellarDirect%u,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayStellarDirect%u = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayStellarDirect%v,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayStellarDirect%v = tempDoubleArray 
     deallocate(tempDoubleArray)

     allocate(tempDoubleArray(SIZE(yArrayStellarScattered)))
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayStellarScattered%i,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayStellarScattered%i = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayStellarScattered%q,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayStellarScattered%q = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayStellarScattered%u,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayStellarScattered%u = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayStellarScattered%v,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayStellarScattered%v = tempDoubleArray 
     deallocate(tempDoubleArray)

     allocate(tempDoubleArray(SIZE(yArrayThermalDirect)))
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayThermalDirect%i,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayThermalDirect%i = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayThermalDirect%q,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayThermaldirect%q = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayThermalDirect%u,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayThermalDirect%u = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayThermalDirect%v,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayThermalDirect%v = tempDoubleArray 
     deallocate(tempDoubleArray)

     allocate(tempDoubleArray(SIZE(yArrayThermalScattered)))
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayThermalScattered%i,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayThermalScattered%i = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayThermalScattered%q,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayThermalScattered%q = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayThermalScattered%u,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayThermalScattered%u = tempDoubleArray 
     tempDoubleArray = 0.0_db
     call MPI_REDUCE(yArrayThermalScattered%v,tempDoubleArray,SIZE(yArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     yArrayThermalScattered%v = tempDoubleArray 
     deallocate(tempDoubleArray)



#endif
     call torus_mpi_barrier !('finished syncing output. Waiting to continue...') 
     write(message,'(i10,a)') int(real(iOuterLoop)/real(nOuterLoop)*real(nPhotons)+0.5)," photons done"
     call writeInfo(message, TRIVIAL)


        call torus_mpi_barrier

      end subroutine do_one_outer_photon_loop


  subroutine define_rotation_axis

    type(VECTOR) :: tempVec
    type(VECTOR) :: normToRotation

    rotationAxis = zHat

    !
    ! Redefining the rotation axis here.
    !
    select case (grid%geometry)
    case ("ttauri")
       rotationAxis = grid%diskNormal
       if (writeoutput) write(*,*) "rotation axis",rotationAxis
    case("luc_cir3d","cmfgen") 
       rotationAxis = VECTOR(0.,0.,1.)
       if (writeoutput) write(*,*) "rotation axis",rotationAxis
    case("romanova" )
       rotationAxis = grid%diskNormal
       if (writeoutput) write(*,*) "rotation axis",rotationAxis
    case ("donati")
       rotationAxis = rotateX(rotationAxis, -dble(dipoleOffset))
    case("jets")
       rotationAxis = VECTOR(1., 0., 0.)
    end select

    if (distortionType == "spiral") then
       rotationAxis = rotateX(rotationAxis, dble(dipoleOffset))
    endif

    normToRotation = rotationAxis .cross. yHat
    call normalize(normToRotation)

! Second argument was inclination parameter in version 1. Removed in version 2 as 
! I don't think it is required. DMA, June 2010. 
    tempVec  = arbitraryRotate(rotationAxis, 0.0_db, normToRotation)
    originalViewVec =  (-1.d0)*tempVec

  end subroutine define_rotation_axis

end subroutine do_phaseloop

!-----------------------------------------------------------------------------------------------------------------------

subroutine choose_view ( geometry, nPhase, distortionType, doRaman, &
       rotateview, rotateDirection, tiltView)

  implicit none

! Intent in arguments
  character(len=*), intent(in) :: geometry
  integer, intent(in)          :: nPhase
  character(len=*), intent(in) :: distortionType
  logical, intent(in)          :: doRaman

! Intent out arguments
  logical, intent(out) :: rotateView
  real, intent(out)    :: rotateDirection
  logical, intent(out) :: tiltView

! Set default values of intent out arguments
  rotateView      = .false.
  rotateDirection = -1.0
  tiltView        = .false. 

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


  if (geometry == "warpeddisc") rotateView = .true.

  if ((geometry == "ttauri").or.(geometry == "magstream")) rotateView = .true.

  if ((geometry == "toruslogo")) then
     rotateView = .true.
     rotateDirection = 1.
  endif
  if ((geometry == "luc_cir3d")) rotateView = .true.

  if ((geometry == "cmfgen")) rotateView = .true.

  if ((geometry == "romanova")) rotateView = .true.

  if ((geometry(1:4) == "jets")) then
     rotateView = .true.
     tiltView = .true.
  end if

  if (doRaman) then
        rotateView = .true.
        rotateDirection = 1.
  endif

  if (geometry == "rolf") rotateView = .true.

  if (geometry == "disk") rotateView = .true.

  if (geometry == "starburst") rotateView = .true.

  if (geometry == "wr104") then
     rotateView = .true.
     rotateDirection = 1.
  endif

  if (geometry(1:7) == "cluster") then
     rotateView = .true.
!     tiltView = .true.
  end if

end subroutine choose_view

end module phaseloop_mod

