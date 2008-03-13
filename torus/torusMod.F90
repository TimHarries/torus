!
!            T   O   R   U  S
!
! torus is a general purpose monte-carlo radiative transfer code used to
! produce polarized spectra of hot and cool stellar winds
!

! written by tjh
! Extracted from torusMain by D. Acreman (March 2008)

  module torus_mod

    implicit none 

    public :: torus

    private

    contains

 subroutine torus(b_idim,b_npart,b_nactive,b_xyzmh,b_rho,b_iphase, &
                  b_nptmass,b_listpm,b_udist,b_umass, &
                  b_utime,b_time,b_gaspartmass, b_num_gas, b_temp)


  use constants_mod          ! physical constants
  use kind_mod               ! variable type KIND parameters
  use vector_mod             ! vector math
  use photon_mod             ! photon manipulation
  use gridtype_mod           ! type definition for the 3-d grid
  use grid_mod               ! opacity grid routines
  use phasematrix_mod        ! phase matrices and stokes vectors
  use math_mod               ! misc maths subroutines
  use blob_mod               ! clumps initialization and movement
  use distortion_mod         ! for distorting opacity grid
  use image_mod              ! stokes image
  use utils_mod
  use stateq_mod
  use inputs_mod
  use TTauri_mod
  use luc_cir3d_class
  use cmfgen_class
  use romanova_class
  use unix_mod
  use path_integral
  use puls_mod
  use input_variables         ! variables filled by inputs subroutine
  use lucy_mod
  use amr_mod
  use jets_mod
  use dust_mod
  use source_mod
  use spectrum_mod
  use wr104_mod
  use sph_data_class
  use cluster_class
  use timing
  use isochrone_class
  use filter_set_class
  use cluster_utils
  use surface_mod
  use disc_hydro_mod
  use disc_class
  use discwind_class
  use jet_class
  use gas_opacity_mod
  use diffusion_mod
  use messages_mod
  use photoion_mod
  use photoionAMR_mod
  use formal_solutions
  use starburst_mod
  use molecular_mod
  use modelatom_mod
  use cmf_mod
  use hydrodynamics_mod
  use mpi_global_mod
  use parallel_mod

  implicit none
#ifdef MPI
   include 'mpif.h'  
#endif

  integer :: nOuterLoop = 10
  integer :: nInnerLoop
  integer :: iLambdaPhoton
  integer :: iOuterLoop
  integer :: nScat
  integer :: cpuTime, startTime
  integer :: nSource
  type(SOURCETYPE), allocatable :: source(:)
  type(SOURCETYPE) a_star
  real(double) :: lCore
  real ::  weightDust=1.0, weightPhoto=1.0
  real :: inclination, cos_inc, cos_inc_first, cos_inc_last, d_cos_inc
  integer :: iInclination  
  
  character(len=80) :: thisdevice

  real(double) :: energyPerPhoton
  real(double) :: objectDistance

  ! variables for the grid

  type(GRIDTYPE) :: grid

  integer :: iPhase
  real(double) :: totLineEmission
  real(double) :: totContinuumEmission
  real(double) :: totCoreContinuumEmission
  real(double) :: cOrecontinuumflux
  real(double) :: totCoreContinuumEmission1
  real(double) :: totCoreContinuumEmission2
  real(double) :: totWindContinuumEmission
  real :: probLinePhoton 
  real :: weightContPhoton, weightLinePhoton
  real :: chanceLine, chanceContinuum
!  real :: sTot

  integer :: isize
  integer, allocatable :: iseed(:)

  ! optical depth variables
  !
  ! for monte calro with high optical depth models
  !  integer, parameter :: maxTau = 2000000
  ! for monte calro with optically thin model.  
  integer, parameter :: maxTau = 20000
  ! for direct integration
  !  integer, parameter :: maxTau = 8000
  !
  ! for Chris
  !  integer, parameter :: maxTau = 600000
  integer :: nTau
  real, allocatable :: contTau(:,:)

  real :: scaleFac
  real :: dopShift  ! velocity shift in terms of thermal line width
  real, allocatable :: tauExt(:)
  real, allocatable :: tauAbs(:)
  real, allocatable :: tauSca(:)
  real, allocatable :: linePhotonAlbedo(:)
  real, allocatable :: lambda(:)

  real  :: sigmaExt0, sigmaAbs0, sigmaSca0  ! cross section at the line centre
  real :: dlambda, thisTau

  ! variables to do with dust

  integer :: itestlam, ismoothlam
  integer, parameter :: nXmie = 20, nMuMie = 1000
  type(PHASEMATRIX), allocatable :: miePhase(:,:, :)

  ! torus images

  type(IMAGETYPE) :: o6image(1)

  type(IMAGETYPE), allocatable :: obsImageSet(:)
  integer           :: nImage  ! number of images in obsImageSet
  type(filter_set)  :: filters ! a set of filters used for imaging
  character(LEN=30) :: name_filter
  real(double) :: bandwidth   ! Band width of a filter[A]
  real(double) :: lambda_eff  ! Effective wavelength of a filter[A]

  type(PVIMAGETYPE), allocatable :: pvimage(:)
  real :: imageSize
  integer :: iSlit
  type(VECTOR) :: slitPosition

  ! intrinsic profile variables

  integer, parameter :: maxIntPro = 1000
  integer :: nIntPro, iSpline
  real :: lamIntPro(maxIntPro), intPro(maxIntPro)
  real, allocatable :: splineArray(:)

  real, allocatable :: statArray(:)
  real, allocatable :: sourceSpectrum(:)
  real, allocatable :: sourceSpectrum2(:)
!  logical :: resonanceLine


  ! variables for clumped wind models
  
  integer, parameter :: maxBlobs = 10000
  integer :: nCurrent
    type(BLOBTYPE), allocatable :: blobs(:)
  real, parameter :: blobTime = 1000.
  real :: timeEnd = 24.*60.*60.
  real :: timeStart = 0.
  real :: dTime

  ! filenames

  character(len=80) :: filename, specFile, obsfluxfile
  character(len=80) :: originalOutFile
  character(len=80) :: phasePopFilename

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
  type(OCTALVECTOR) :: rHat, rVec, rHatinStar

! Used in commented out call to findSubcellTD
!  type(octal), pointer :: thisOctal
!  integer :: subcell

  ! output arrays

  integer :: iLambda
  type(STOKESVECTOR), allocatable :: yArray(:)
  type(STOKESVECTOR), allocatable :: errorArray(:,:)
  real, allocatable :: xArray(:)
  real, allocatable :: contWeightArray(:)

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

  real :: meanDustParticleMass
!  real :: getMeanMass2  ! the call to this function is currently commented out


  real :: foreground = 0., background = 0.  ! plotting intensities
  real(oct) :: t1, t2, t3

  logical :: velocitySpace

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

  real(double) :: thisChi, thisSca, albedo
  integer :: currentScat
  logical :: redRegion
  real :: r
  logical :: contWindPhoton
  real :: ramanWeight
  real :: thisVel
  real :: fac1, fac2, fac3
  real(double) :: vRay, vOverCsqr
  real :: directionalWeight
  real :: obs_weight, tau_tmp, exp_minus_tau
  real :: r1, r2
  real(oct) :: t
  integer :: i, j
  integer :: i1,i2,i3
  integer :: nTot  
  real :: observedLambda
  real :: thisLam
  real :: nu
!  real :: fac
  real(double) :: fac
  real(double) :: tau_bnd
  real :: escProb
  logical :: contPhoton
  integer :: nContPhotons
  real :: phi
  real :: rStar
  character(len=80) :: tempChar
  logical :: lineResAbs    ! T if you want to include absorption
  !                        !   of line by line resonance zones.


  !
  real(double) :: totalMass
  real(double) :: T_ave   ! average temperature of cluster
  real(double) :: T_mass  ! mass weighted temperature
!  real(double) :: tau_max, tau_min, tau_ave

  real(double) :: Laccretion
  real :: Taccretion, fAccretion, sAccretion
  real :: theta1, theta2, chanceHotRing
  type(SURFACETYPE) :: starSurface

  ! Spot stuff
  
  real :: chanceSpot                     ! chance of spot
  logical :: spotPhoton                  ! photon from spot?

  real :: chanceDust = 0.
  real(double) :: totDustContinuumEmission, totEnvelopeEmission

  ! binary parameters

  integer :: nVec
  type(VECTOR), allocatable :: distortionVec(:)

  ! adaptive grid stuff

  type(OCTALVECTOR) :: amrGridCentre ! central coordinates of grid
  type(OCTALVECTOR) :: octVec
  real :: ang
  real(double) :: kabs, eta
  integer :: nt
  integer           :: nOctals       ! number of octals in grid
  integer           :: nVoxels       ! number of unique voxels in grid
                                     !   (i.e. the number of childless subcells)
  logical :: gridConverged           ! true when adaptive grid structure has 
                                     !   been finalised
  integer :: intPathError            ! error code from integratePathAMR
  type(OCTALVECTOR) :: positionOc    ! photon position position
!  type(octal),pointer :: octalLocation! octal located by search routine 
!  integer           :: subcellLocation! subcell located by search routine
  integer           :: tooFewSamples ! number of errors from integratePathAMR
  integer           :: boundaryProbs ! number of errors from integratePathAMR
  integer           :: negativeOpacity ! number of errors from integratePathAMR
  real(double)           :: Ne         ! for testing
  real(double),dimension(statEqMAxLevels) :: levelPops  ! for testing
  integer                         :: level      ! for testing
  character(len=80) :: newContFluxFile ! modified flux file (i.e. with accretion)
  real :: infallParticleMass         ! for T Tauri infall models
  logical :: alreadyDoneInfall = .false. ! whether we have already done an infall calculation
  type(alpha_disc)  :: ttauri_disc       ! parameters for ttauri disc
  type(discwind)    :: ttauri_discwind   ! parameters for ttauri disc wind
  type(jet)         :: ttauri_jet        ! parameters for ttauri jets

!  ! For luc_cir3d geometry case
!  type(luc_cir3d) :: cir3d ! parameters and data for luc_cir3d geometry

  ! For romanova geometry case
  type(romanova) :: romData ! parameters and data for romanova geometry

  !
  ! SPH data of Matthew
  type(sph_data) :: sphData

  ! Used for multiple sources (when geometry=cluster)
  type(cluster)   :: young_cluster

  ! Used in "plot_AMR_planes" and "plot_AMR_values"
  integer :: nmarker           ! number of markers
  real, allocatable    :: xmarker(:)  ! position of x
  real, allocatable    :: ymarker(:)  ! position of y
  real, allocatable    :: zmarker(:)  ! position of z
  real    :: width_3rd_dim         ! Use this to restrict the markers to be plotted..  
  real val_3rd_dim

  ! Name of the file to output various message from torus
  character(LEN=7), parameter :: messageFile = "MESSAGE"
  character(len=80) :: message
   real :: h !, totFrac
!  integer :: nFrac   ! In commented out call to sublimateDust

! molecular line stuff
  type(MOLECULETYPE) :: co
  type(DATACUBE) :: cube

!  real(double) :: temp(50)

  real(double), allocatable :: flux(:)

! Variables used when linking to sph code
  type(isochrone)       :: isochrone_data
  integer, intent(in)   :: b_idim,b_npart,b_nactive,b_nptmass
  integer*1, intent(in) :: b_iphase(b_idim)
  integer, intent(in)   :: b_listpm(b_nptmass)
  real*8, intent(in)    :: b_xyzmh(5,b_idim)
  real*4, intent(in)    :: b_rho(b_idim)
  real*8, intent(in)    :: b_udist, b_umass, b_utime, b_time, b_gaspartmass
  integer, intent(in)   :: b_num_gas           ! Number of gas particles
  real*8, intent(inout) :: b_temp(b_num_gas)   ! Temperature of gas particles
  integer :: ngaspart
  integer, save :: num_calls = 0

  type(modelatom), allocatable :: thisAtom(:)

  type(STREAMTYPE) :: thisStream(2000), bigStream
  integer :: nStreams
  
  integer :: nFromEnv
  logical :: photonFromEnvelope

  integer :: nRBBTrans
  integer :: indexRBBTrans(1000), indexAtom(1000)

  integer :: iInner_beg, iInner_end ! beginning and end of the innerPhotonLoop index.

  real(double) :: tempArray(10)
#ifdef MPI
  ! For MPI implementations =====================================================
  integer ::   ierr           ! error flag
  integer ::   tempInt        !
  real, dimension(:), allocatable :: tempRealArray
  real, dimension(:), allocatable :: tempRealArray2
  real(double), dimension(:), allocatable :: tempDoubleArray
  real(double), dimension(:), allocatable :: tempDoubleArray2
  integer, dimension(:), allocatable :: photonBelongsRank
  integer, parameter :: tag = 0
  logical :: rankComplete
  
! Begin executable statements --------------------------------------------------

  ! FOR MPI IMPLEMENTATION=======================================================

  ! Set up amrCOMMUNICATOR and global mpi groups
  call setupAMRCOMMUNICATOR

  call unixGetHostname(tempChar, tempInt) 
  print *, 'Process ', myRankGlobal,' running on host ',TRIM(ADJUSTL(tempChar))

  !===============================================================================

#endif

  writeoutput    = .true.
  outputwarnings = .true.
  outputinfo     = .true.
  doTuning       = .true.
  myRankIsZero   = .true.
#ifdef MPI
  if (myRankGlobal/=1) writeoutput  = .false.
  if (myRankGlobal/=1) doTuning     = .false.
  if (myRankGlobal/=0) myRankIsZero = .false.
#endif
  
  ! For time statistics
  if (doTuning) call tune(6, "Torus Main") ! start a stopwatch  

  ! set up a random seed
  
  call random_seed

  call random_seed(size=iSize)
  allocate(iSeed(1:iSize))
  call random_seed(get=iSeed)

  call torus_mpi_barrier
#ifdef MPI
  call MPI_BCAST(iSeed, iSize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call random_seed(put=iseed)
! write(*,*) myRankGlobal,"random seed",iseed(1:isize)
#endif
  deallocate(iSeed)


  ! initialize

  useNdf = .true.
  sed = .false.
  jansky = .false.
  SIsed = .false.
  velocitySpace = .false.
  movie = .false.
  thinLine = .false.
  flatspec = .false.
  greyContinuum = .false.
  secondSource = .false.
  firstplot = .true.
  doRaman = .false.
  enhance = .false.
  lineResAbs = .false.  

  contFluxFile = "none"
  intProFilename = "none"

  inputKappaSca = 0.
  inputKappaAbs = 0.

  zeroVec = VECTOR(0.,0.,0.)

  nMarker = 0
  tooFewSamples = 0 
  boundaryProbs = 0 
  negativeOpacity = 0 
  lucyRadiativeEq = .false. ! this has to be initialized here
  num_calls = num_calls + 1

  ! hardwired stuff
  do i = 1, no6pts
     o6xArray(i) = o6start + (o6end-o6start)*real(i-1)/real(no6pts-1)
     o6yarray(i) = 1.e-10
  enddo


  ! get the model parameters

  call inputs() ! variables are passed using the input_variables module
  if (.not.inputOK) goto 666

!  call test_profiles()  ! Testing Lorentz profile with Voigt profile


  if(geometry == "planetgap") then
     ! takes the gapWidth and calculates what mPlanet should be
     if (planetgap) then
        call calcPlanetMass
        write (*,*) "mPlanet set to ", mPlanet
     end if
  end if

  if (molecular) then
     if (geometry .eq. "molebench") then
        call readbenchmarkMolecule(co, "hco_benchmark.mol")
     elseif((geometry .eq. "h2obench1") .or. (geometry .eq. "h2obench2")) then
        call readMolecule(co, "fakewater.mol")
     elseif(geometry .eq. "agbstar") then
        call readMolecule(co, "p-h2o.dat")
     elseif(geometry .eq. "molefil") then
        call readMolecule(co, "c18o.mol")
     else
        call readMolecule(co, "co.mol")
     endif
  endif

  if (cmf) then
     allocate(thisAtom(1:nAtom))
     do i = 1, nAtom
        call readAtom(thisAtom(i),atomFilename(i))
        call stripAtomLevels(thisAtom(i), 5)
     enddo
    call createRBBarrays(nAtom, thisAtom, nRBBtrans, indexAtom, indexRBBTrans)
  endif



  objectDistance = griddistance * pctocm


  grid%photoionization = photoionization
  ! time elipsed since the beginning of the model
  grid%timeNow = phaseTime * REAL(nStartPhase-1)
  

  grid%resonanceLine = resonanceLine


  amrGridCentre = OCTALVECTOR(amrGridCentreX, amrGridCentreY, amrGridCentreZ)

  
  ! For plotting routines used later in this program
  if (plane_for_plot == "x-y") then
     val_3rd_dim = amrGridCentre%z
  else if (plane_for_plot == "y-z") then
     val_3rd_dim = amrGridCentre%x
  else if (plane_for_plot == "z-x") then
     val_3rd_dim = amrGridCentre%y
  else if (plane_for_plot == "x-z") then
     val_3rd_dim = amrGridCentre%y
  else
     val_3rd_dim = 0.0d0        
  end if

  !

  if (geometry.eq."raman") then
     hotSourcePosition = 0.75*secondSourcePosition
     coolStarPosition = (-0.25)*secondSourcePosition
     secondSourcePosition = hotSourcePosition
  endif

  if (doRaman) screened = .true.

  if (dopvimage) then
     allocate(pvimage(1:nSlit))
  endif

!  if (mie .or. (geometry == "ttauri" .and. ttau_disc_on)) then
!     meanDustParticleMass = getMeanMass2(amin(1), amax(1), a0(1), qDist(1), pdist(1), graintype(1))
!     if (meanDustParticleMass <=0.0) then 
!        write(message,*) "Error: meanDustParticleMass <=0 in torusMain."
!        call writeFatal(message)
!        stop
!     end if
!  endif

  probLinePhoton = 1. - probContPhoton

  if (allocated(inclinations)) then
     inclination = inclinations(1)
  else
     inclination = firstInclination
  end if
  if (inclination .lt. 1.e-6) then
     call writeWarning("inclination less than 1.e-6, now set to 1.e-6.")
     call writeWarning("Did you specify an inclination of 0?")
  end if
  inclination = max(inclination, 1.e-4)

  scale = scale * rSol

  rStar = 100.*rSol

  ! Choose whether to rotate the view
  call choose_view( geometry,   nPhase,          distortionType, doRaman, & ! Intent in
                    rotateview, rotateDirection, tiltView)                  ! Intent out

  ! switches for line emission

  if (lineEmission) then
     flatspec = .true.
     greyContinuum = .true.
     velocitySpace = .true.
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
  viewVec = originalViewVec

  !=====================================================================
  ! INIIALIZATIONS BEFORE CALLING INITAMR ROUTINE SHOULD BE HERE
  !=====================================================================

  !
  ! Special case
  !

  if (geometry == "cluster") then

     ! The total number of gas particles is the total number of active particles less the number of point masses.
     ngaspart = b_nactive-b_nptmass
     call init_sph_data2(sphData, b_udist, b_umass, b_utime, ngaspart, b_time, b_nptmass, &
          b_gaspartmass, b_npart, b_idim, b_iphase, b_xyzmh, b_rho, b_temp)
! Communicate particle data. Non-mpi case has a stubbed routine. 
     call gather_sph_data(sphData)

        ! Writing basic info of this data
     if (myRankIsZero) call info(sphData, "info_sph.dat")

     ! reading in the isochrone data needed to build an cluster object.
     call new(isochrone_data, "dam98_0225")   
     call read_isochrone_data(isochrone_data)
     
     ! making a cluster object
     call new(young_cluster, sphData, dble(amrGridSize), disc_on)
     call build_cluster(young_cluster, sphData, dble(lamstart), dble(lamend), isochrone_data)
    
     ! Wrting the stellar catalog readble for a human
     if (myRankIsZero) call write_catalog(young_cluster, sphData)

      ! Finding the inclinations of discs seen from +z directions...
!     if (myRankIsZero) call find_inclinations(sphData, 0.0d0, 0.0d0, 1.0d0, "inclinations_z.dat")

  end if

  if (geometry == "wr104") then
     if (.not.(readPops.or.readlucy)) then
        call readWR104Particles("harries_wr104.txt", sphData, objectDistance)
        call info(sphData,"*")
     endif

  elseif (geometry == "ttauri".or.(geometry == "magstream")) then
     onekappa=.false.
     if (ttau_disc_on)  &
          call new(ttauri_disc, dble(TTauriDiskRin*TTauriRstar/1.0e10), 1.5d3*100.d0, &
          dble(TTauriMstar/100.0), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, dble(TTauriMstar/mSol))
     if (ttau_discwind_on) &
          call new(ttauri_discwind, DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
          DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, &
          dble(TTauriMstar/mSol),  dble(TTauriDiskHeight/1.0e10) )
     if (ttau_jet_on) &
!          call new(ttauri_jet,  dble(TTauriRouter/1.0e10), dble(amrgridsize), &
          call new(ttauri_jet,  JET_Rmin, dble(amrgridsize), &
          JET_theta_j, dble(TTauriMstar/mSol), JET_Mdot, JET_a_param, JET_b_param,  &
          JET_Vbase, JET_Vinf, JET_beta, JET_gamma, limitscalar2, JET_T)

  elseif (geometry == "luc_cir3d") then
     onekappa=.false.
     !                    Rmin [10^10cm]         
     call new(CIR_Rstar*Rsol*1.0e-10, CIR_Twind, CIR_Mdot_scale,  &
     "zeus_rgrid.dat", "zeus_rho_vel.dat")     

  elseif (geometry == "cmfgen") then
     onekappa=.false.
     call read_cmfgen_data("OPACITY_DATA")
     call put_cmfgen_Rmin(CMFGEN_Rmin)  ! in [10^10cm]
     call put_cmfgen_Rmax(dble(amrGridSize)/2.0d0)  ! in [10^10cm]

  elseif (geometry == "romanova") then
     onekappa=.false.
     ! create the romanova data object and readin the data at the same time.
     ! After this call, the obejct will be ready to use. 
     call new(romData, ROM_Rs*DBLE(Rsol)*1.0d-10, 70.0d0, &
          ROM_Mass*DBLE(Msol), ROM_isoT, ROM_T_flow, ROM_r_ref, &
          ROM_rho_ref, ROM_T_ref, ROM_v_ref, ROM_tilt, ROM_period, ROM_datafile)

     ! this is a temporary solution... It should be called something 
     ! else .
     TTauriRouter = REAL(get_dble_parameter(romData, "Rmax")*1.0d10)           

  end if
  !==========================================================================
  !==========================================================================




  ! allocate the grid - this might crash out through memory problems

  if (gridUsesAMR) then
     ! any AMR allocation stuff goes here
     call initAMRGrid(greyContinuum, &
                        newContFluxFile,flatspec,grid,ok,theta1,theta2)
        if (.not.ok) goto 666

  else
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
  end if

  ! a few dynamic arrays

!  if (mie) sed = .true.

  grid%splitOverMpi = .false.

#ifdef MPI
 if (hydrodynamics) grid%splitOverMpi = .true.
 if (photoionization) grid%splitOverMpi = .true.
#endif

  grid%resonanceLine = resonanceLine
  grid%lambda2 = lamline

  ! Note: the first index should be either lambda or mu
  !       in order to speedup the array operations!!!  (RK) 
  allocate(miePhase(1:nDustType,1:nLambda,1:nMumie)) 
  allocate(statArray(1:nLambda))
  allocate(sourceSpectrum(1:nLambda))
  allocate(sourceSpectrum2(1:nLambda))

  grid%doRaman = doRaman

  statArray = 0.
  sourceSpectrum = 1.
  sourceSpectrum2 = 1.


  ! allocate the output arrays

  if (mie) then
     nOuterLoop = nLambda
  endif

  allocate(xArray(1:nLambda))
  allocate(yArray(1:nLambda))
  allocate(errorArray(1:nOuterLoop,1:nLambda))


  if (photoionization) then
     call addIons(grid%ion, grid%nion)
  endif

! Set up the wavelength arrays with either linear or log spaced values.
  call set_up_lambda_array

  !
  ! Setting the number of opacity (kappa) arrays in grid.
  !
  if (flatspec) then
     grid%nopacity = 1
  else
     grid%nopacity = nLambda
  end if


  errorArray(1:nOuterLoop,1:nLambda) = STOKESVECTOR(0.,0.,0.,0.)

  if (doRaman) then
     yarray(1:nLambda)%i = 1.e-20
     errorArray(1:nOuterLoop,1:nLambda)%i = 1.e-20
  endif


  call  createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)

  if (noScattering) then
     if (writeoutput) write(*,*) "! WARNING: Scattering opacity turned off in model"
     grid%oneKappaSca(1:nDustType,1:nLambda) = TINY(grid%oneKappaSca)
  endif



!  tempArray(1) = 1.d0
!  call testMiePhase(xarray(1), xArray, nLambda, miePhase, nDustType, nMuMie, temparray)
  if (writeoutput) then
     open(76, file="phasematrix.dat",status="unknown",form="formatted")
     ilambda = findIlambda(10000.0, xArray, nLambda, ok)
     do i = 1, nMuMie
	ang =  pi*real(i-1)/real(nMuMie-1)
        write(76,*) 180.*ang/pi, miephase(1, ilambda, i)%element(1,1)
     enddo
     close(76)
     tempArray(1) = 0.
     do i = 1, nMumie
        tempArray(1) = tempArray(1) + miePhase(1, ilambda, i)%element(1,1)/real(nMuMie)
     enddo
     write(*,*) "phase normalization ",tempArray(1)


  endif
!  stop

  if (includeGasOpacity) then

     ! This routine may cause trouble with MPI routine! 
     ! Propably multiple nodes are trying to write to a same file at a same time.

     call readTsujiPPTable()
     call readTsujiKPTable()
     call createAllMolecularTables(20, 1.e-6, 20000., grid%nLambda, xArray)

  end if

  ! if the grid uses an adaptive mesh, create it
     

!  if (gridUsesAMR) then
!     if ((mie) .or. (geometry == "ttauri" .and. ttau_disc_on)) then 
!        call setKappaTest(grid, scale, aMin, aMax, a0, qDist, pDist, grainType, &
!             ngrain, X_grain, grainname, lambdaTau)
!     end if
!  end if


  if (mie) then
     call locate(grid%lamArray, nLambda,lambdaTau,itestlam)
     call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)
     grid%itestlam = iTestLam
     write(message,*) "Test wavelength index: ",itestlam,ismoothlam
     call writeInfo(message, TRIVIAL)
  end if

  ! chris
  if (hydroWarp) then
     print *, "Creating and reading hydroGrid..."
     allocate(grid%hydroGrid)
     call initAMRGrid(greyContinuum, &
                        newContFluxFile,flatspec,grid%hydroGrid,ok,theta1,theta2)
     call readAMRgrid("hydrogrid.dat",readFileFormatted,grid%hydroGrid)
     print *, "hydroGrid read..."
     call hydroWarpFitSplines(grid%hydroGrid)
     print *, "splines fitted..."
     print *, "hydroGrid setup done."
  end if


  !===============================================================
  if (gridUsesAMR) then  !========================================
  !===============================================================


! below is really yucky code - please change me

     if (molecular .and. readmol .and. (.not. lucyRadiativeEq)) then
        if(.not. openlucy) then
           call readAMRgrid(molfilenamein,.false.,grid)
        else
           call readAMRgrid(lucyFileNamein,.false.,grid)
        endif
        goto 667
     endif

     if (cmf.and.readLucy) then
        call readAMRgrid("atom_tmp.grid",.false.,grid)
     else if (photoionization.and.readlucy) then
        continue
     else
        ! Set up the AMR grid
        call amr_grid_setup

        ! Write the information on the grid to file using a routine in grid_mod.f90
        if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")

        ! Plotting the various values stored in the AMR grid.
        call do_amr_plots
     end if

  !=================================================================
  else ! grid is not adaptive ======================================     
  !=================================================================


     print *, "ERROR: adaptive grid required when using torusMod"
     goto 666

  !=================================================================
  endif ! (gridusesAMR)
  !=================================================================


  outVec = (-1.)* originalViewVec

  ! set up the sources
  call set_up_sources

     if (geometry == "wr104") then
!        call IntegratePathAMR(lambdatau,  lamLine, VECTOR(1.,1.,1.), zeroVec, &
!             VECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
!             tauSca, maxTau, nTau, opaqueCore, escProb, .false. , &
!             lamStart, lamEnd, nLambda, contTau, hitCore, thinLine, .false., &
!             .false., nUpper, nLower, 0., 0., 0., junk,&
!             sampleFreq,intPathError)
!        if (intPathError < 0) then
!           write(*,*) '   Error encountered in cross-sections!!! (error = ',intPathError,')'
!        end if
        totalMass = 0.d0
        call findTotalMass(grid%octreeRoot, totalMass)
        scaleFac = massEnvelope / totalMass
        if (writeoutput) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
        call scaleDensityAMR(grid%octreeRoot, scaleFac)

     else if (geometry == "cluster") then
!        call find_average_temperature(grid, T_ave, T_mass, totalMass)
!        if (writeoutput) write(*,*) " "
!        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Total dust mass in cluster  : ", TotalMass, " [g]"
!        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Ave. temperature of cluster : ", T_ave,    " [K]"
!        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Mass weighted ave. temperature of cluster : ",T_mass,  " [K]"
!        if (writeoutput) write(*,*) " "
!        ! computing the max, min and average tau at 2 micron (of all
!        ! active cells)
!        ilambda = findIlambda(20000.0, xArray, nLambda, ok)
!        call find_max_min_ave_tau(grid, ilambda, tau_max, tau_min,  tau_ave)
!        if (writeoutput) write(*,*) " "
!        if (writeoutput) write(*,'(a,1pe12.4, a11)') "Max tau(1 micron) = ", tau_max, " [contiuum]"
!        if (writeoutput) write(*,'(a,1pe12.4, a11)') "Min tau(1 micron) = ", tau_min, " [contiuum]"
!        if (writeoutput) write(*,'(a,1pe12.4, a11)') "Ave tau(1 micron) = ", tau_ave, " [contiuum]"
!        if (writeoutput) write(*,*) " "

     else if (geometry == "ttauri") then
        call find_average_temperature(grid, T_ave, T_mass, totalMass)
        if (writeoutput) write(*,*) " "
        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Total mass in computational domain : ", TotalMass, " [g]"
        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Ave. temperature in computational domain : ", T_ave,    " [K]"
        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Mass weighted ave. temperature in computational domain : ",T_mass,  " [K]"
        if (writeoutput) write(*,*) " "
     end if


     if (associated(grid%octreeRoot)) then
        totalMass =0.d0
        call findTotalMass(grid%octreeRoot, totalMass)
        write(message,*) "Mass of envelope: ",totalMass/mSol, " solar masses"
        call writeInfo(message, IMPORTANT)
     endif
        

     if (mie) then
        if (geometry == "shakara" .and. geometry .eq. 'iras04158') then

!        sigma0 = totalMass / (twoPi*(rOuter*1.e10-rInner*1.e10)*1.*real(autocm)) ! defined at 1AU

!           sigma0 = totalMass * real(betaDisc-alphaDisc+2) / & 
!                (twoPi*((rOuter*1.e10)**(betaDisc-alphaDisc+2) - (rInner*1.e10)**(betaDisc-alphaDisc+2))*1.*auToCm)

!           sigma0 = rho0 * (rinner*1.e10/(1.*autocm))**alphaDisc * &
!           (height*1.e10) * (autocm / (100.d0*autocm))**betaDisc * sqrt(twopi)

           sigma0 = rho0 * (height*1.e10) * ((rinner*1.e10) / (100.d0*autocm))**betaDisc * sqrt(twopi)

           if (writeoutput) write(*,*) "Sigma0: ",sigma0
        endif

        if (geometry == "ppdisk") then
           sigma0 = totalMass / ((amrgridsize**2 - twoPi*0.4)*1.e10*1.*autocm) ! defined at 1AU
           if (writeoutput) write(*,*) "Sigma0: ",sigma0
        endif

        if (geometry == "planetgap") then
           rho0  = mDisc *(betaDisc-alphaDisc+2.) / ( twoPi**1.5 * height * (rCore*1.e10) &
                * (rCore*1.e10)**(alphaDisc-betaDisc) * &
                (((rOuter*1.e10)**(betaDisc-alphaDisc+2.)-(rInner*1.e10)**(betaDisc-alphaDisc+2.))) )
           h = height * rCore * 1.e10
           sigma0 = rho0 * sqrt(twoPi) * h

           if (writeoutput) write(*,*) "Sigma0: ",sigma0
        endif

        if (geometry == "warpeddisc") then
           rho0  = mDisc *(betadisc-alphadisc+2.) / ( twoPi**1.5 * (height*1.e10)  &
                   * (rOuter*1.d10)**(alphadisc-betadisc) * ( &
                   ((router*1.d10)**(betadisc-alphadisc+2.)-(rInner*1.d10)**(betadisc-alphadisc+2.))) )
           sigma0 = rho0 * sqrt(twoPi) * height * 1.e10
        endif

     endif

     call torus_mpi_barrier

667 continue
     call random_seed

     if (cmf) then
        if (movie) call plotAMRthreeDMovie(grid, source, nsource)
        if (.not.readlucy) call atomLoop(grid, nAtom, thisAtom, nsource, source)
!        call atomLoop(grid, nAtom, thisAtom, nsource, source)

           if (writeoutput) then
              call plot_AMR_values(grid, "n=3", plane_for_plot, val_3rd_dim,  &
                   "etaline.ps/vcps", .true., .false., & 
                   nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim)
              !
              
              fac = 120.*degtorad
              do i = 1, 1
                 ang = real(i-1)/40. * twoPi
                 octVec = OCTALVECTOR(sin(fac)*cos(ang), sin(fac)*sin(ang), cos(fac))
                 write(thisdevice,'(a,i3.3,a)') "col",i,".ps/vcps"
!                 call columnDensityPlot(grid, source, nsource, octVec, thisDevice)
              enddo
           endif

        do i = 1, 50
           ang = twoPi * real(i-1)/51. 
!           ang = 45.*degtorad
           t = 120.d0*degtorad
           call calculateAtomSpectrum(grid, thisAtom, nAtom, iTransAtom, iTransLine, &
                OCTALVECTOR(sin(t)*cos(ang), sin(t)*sin(ang), cos(t)), 100.d0*pctocm/1.d10, &
                source, nsource,i)
        enddo
        stop
     endif

#ifdef MPI
     if (photoIonization) then
        grid%splitOverMPI = .true.
        grid%photoionization = .true.
           call radiationHydro(grid, source, nSource, nLambda, xArray, readlucy, writelucy, &
              lucyfilenameout, lucyfilenamein)
           call torus_mpi_barrier
     endif

     if (hydrodynamics) then

        grid%splitOverMPI = .true.
	

        if (myRankGlobal /= 0) then
           if (grid%octreeRoot%twoD) then
              call doHydrodynamics2d(grid)
           else if (grid%octreeRoot%oneD) then
              call doHydrodynamics1d(grid)
           else if (grid%octreeRoot%threeD) then
              call doHydrodynamics3d(grid)
           endif
        endif
     endif
#endif


  if (lucyRadiativeEq) then
     if (doTuning) call tune(6, "LUCY Radiative Equilbrium")  ! start a stopwatch
 
     if (.not.grid%adaptive) then
        call lucyRadiativeEquilibrium(grid, miePhase, nDustType, nMuMie, nLambda, xArray, dble(teff), nLucy)
     else
        if (readLucy .and. .not. redoLucy) then
           continue
        else

           if (solveVerticalHydro) then
              call verticalHydrostatic(grid, mCore, sigma0, rInner, miePhase, nDustType, nMuMie, nLambda, xArray, &
                   source, nSource, nLucy, massEnvelope, tThresh, .false., mDisc)
           else
              call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, & 
                   nLambda, xArray, source, nSource, nLucy, massEnvelope, tthresh, &
                   lucy_undersampled, .false., IterLucy, plot_i=num_calls)
           endif

        endif

        if (myRankIsZero .and. writeLucy) call writeAMRgrid(lucyFilenameOut,writeFileFormatted,grid)
     endif

     if (doTuning) call tune(6, "LUCY Radiative Equilbrium")  ! stop a stopwatch

! If this is an SPH run then finish here ---------------------------------------
     call update_sph_temperature (b_idim, b_npart, b_iphase, b_xyzmh, sphData, grid, b_temp)
     goto 666

     if (grid%geometry(1:7) == "cluster") then

        if (myRankIsZero) then
           ! Finding apperant magnitudes and colors of the stars in cluster
           write (*,*) " "
           write (*,*) "Computing the magnitudes and colors of stars ..."
           ! -- note grid distance here is in [pc]
           call analyze_cluster(young_cluster,s2o(outVec),dble(gridDistance),grid)
        end if

        call torus_mpi_barrier

        ! restricting the sed and images calculations to the star 
        ! with ID number = idx_restrict_star. 
        ! By default idx_restrict_star =0 which means all stars the cluster
        ! are included in the calculations.  Specifty a value in your
        ! parameter file to restrict to one star.  (10 AU cylinder
        ! zone is used to restrict the effective computational domain.
        call restrict(grid%octreeroot, idx_restrict_star, nsource, &
                      source,  s2o(outVec), 7.5d4) ! the last value is 50 AU (in 10^10cm)
	!                      source,  s2o(outVec), 1.5d4) ! the last value is 10 AU (in 10^10cm)

        call reassign_10K_temperature(grid%octreeroot)

        ! delete the cluster object since it won't be used any more.
        call kill_all(young_cluster)
     end if
     


     ! Plotting the slices of planes
     if (myRankIsZero .and. plot_maps) then
       call plot_AMR_planes(grid, "temperature", plane_for_plot, 3, "temperature", &
            .true., .false., nmarker, xmarker, ymarker, zmarker, show_value_3rd_dim)
       call plot_AMR_planes(grid, "etaCont", plane_for_plot, 3, "etaCont", .true., .false., &
            nmarker, xmarker, ymarker, zmarker, show_value_3rd_dim)
    end if

     call torus_mpi_barrier

  endif

  if (molecular) then
     if (writemol) call molecularLoop(grid, co)
     if (readmol) call calculateMoleculeSpectrum(grid, co)
!        call createDataCube(cube, grid, OCTALVECTOR(0.d0, 1.d0, 0.d0), co, 1)

     if (myRankIsZero) call plotDataCube(cube, 'cube.ps/vcps')
     stop
  endif



!  if (grid%geometry == "shakara") then
!     call defineDiffusionZone(grid, .false., .false.)
!  endif

 
! Tidy up and finish the run 

666 continue

if (doTuning) call tune(6, "Torus Main") ! stop a stopwatch  

call writeInfo("TORUS exiting", FORINFO)

call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
call freeGrid(grid)

deallocate(miePhase) 
deallocate(statArray)
deallocate(sourceSpectrum)
deallocate(sourceSpectrum2)
deallocate(xArray)
deallocate(yArray)
deallocate(errorArray)

call torus_mpi_barrier

CONTAINS
!-----------------------------------------------------------------------------------------------------------------------

  subroutine set_up_lambda_array

    real :: deltaLambda
    real :: loglamStart, logLamEnd

    if (lamLinear) then
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

    else

       logLamStart = log10(lamStart)
       logLamEnd   = log10(lamEnd)
       
       do i = 1, nLambda
          xArray(i) = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
          xArray(i) = 10.**xArray(i)
          yArray(i) = STOKESVECTOR(0.,0.,0.,0.)
       enddo

    endif


       if (lamFile) then
          call writeInfo("Reading wavelength points from file.", TRIVIAL)
          open(77, file=lamfilename, status="old", form="formatted")
          nLambda = 0
333       continue
          nLambda = nLambda + 1
          read(77,*,end=334) xArray(nLambda)
          if (writeoutput) write(*,*) nlambda,xArray(nlambda)
          goto 333
334       continue
          nLambda = nLambda - 1
          close(77)
       endif

    !
    ! Copying the wavelength array to the grid
    do i = 1, nLambda
       grid%lamArray(i) = xArray(i)
    enddo
    grid%nLambda = nLambda

  end subroutine set_up_lambda_array

!-----------------------------------------------------------------------------------------------------------------------

  subroutine amr_grid_setup

    real(double) :: mass_scale, mass_accretion_old, mass_accretion_new

    if (doTuning) call tune(6, "AMR grid construction.")  ! start a stopwatch

    if (readPops .or. readPhasePops .or. readLucy) then 

       if (readLucy) call readAMRgrid(lucyFilenameIn,readFileFormatted,grid)

       if (readPhasePops) then ! need to get the right file for the phase
          write(tempChar,'(i3.3)') nStartPhase
          phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
          call readAMRgrid(phasePopFilename,readFileFormatted,grid)
       else ! just read the normal pops file
          if (readpops) call readAMRgrid(popFilename,readFileFormatted,grid)
       end if

       if (forceLineChange) then
          write (message,'(A,I1,A,I1,A)') 'Recalculating for n =',nUpper,' /',nLower,'levels...'
          call writeInfo(message, IMPORTANT)
          call generateOpacitiesAMR(grid, nLower, nUpper)
          do i = 1, nLambda
             grid%lamArray(i) = xArray(i)
          enddo
          print *, '...level change done.'
       end if

       ! for some geometries, it is valid to alter the 'dipoleOffset parameter
       !   for a grid we've already created. 
       if (geometry == "ttauri".or.(geometry == "magstream")) then 
          ! if (ABS(grid%dipoleOffset/dipoleoffset-1.0) > 0.01) then 
          write(message,'(a,f5.2,a)') 'Using new dipole offset value (',&
               dipoleOffset*radToDeg,'deg)'
          call writeInfo(message, IMPORTANT)
          grid%dipoleOffset = dipoleOffset                         
          grid%diskNormal = VECTOR(0.,0.,1.)
          grid%diskNormal = rotateX(grid%diskNormal,grid%dipoleOffSet)
!          end if
       elseif (geometry == "romanova") then
!           write(*,'(a,f5.2,a)') 'Using new dipole offset value (',&
!                dipoleOffset*radToDeg,'deg)'
!           grid%dipoleOffset = dipoleOffset                         
          grid%diskNormal = VECTOR(0.,0.,1.)
!           grid%diskNormal = rotateY(grid%diskNormal,-grid%dipoleOffSet)
!           grid%diskNormal = VECTOR(0.,0.,1.)
       end if

        ! In case, a user is using a new nlambda (number of wavelength bins)
        ! we have to update the wavelength array in the grid.
        ! This may not work for dust calculation...
        ! Copying the wavelength array to the grid
!        if (grid%nlambda /= nlambda) then
!           if (mie .or. (grid%geometry == "ttauri" .and. ttau_disc_on)) then              
!              write(*,*) "Error:: The number of the wavelength bins in the data file "
!              write(*,*) "        does not macth with that specifed in your parameter file."
!              write(*,*) "nlambda(old) = ", grid%nlambda
!              write(*,*) "nlambda(new) = ", nlambda
!              write(*,*) " "
!              write(*,*) "Make the new nlambda same as the old one, otherwise you cannot use"
!              write(*,*) "this data file. This should be changed in future...."
!!              stop
!           else 
!              write(*,*) "Warning:: The number of the wavelength bins in the data file "
!              write(*,*) "      does not macth with that specifed in your parameter file."
!              write(*,*) "nlambda(old) = ", grid%nlambda
!              write(*,*) "nlambda(new) = ", nlambda
!              write(*,*) "==> We recompute the wavelength array with new nlambda value, and continue."
              ! update the wavelength array.
              ! It over writes whatever the wavelength you had in your file!
       deallocate(grid%lamArray)
       allocate(grid%lamArray(nlambda))
       do i = 1, nLambda
          grid%lamArray(i) = xArray(i)
       enddo
       grid%nLambda = nLambda
!           end if
!        end if

        !
        ! If the grid read from file contains ttauri disc or jet, you can turn it
        ! off by setting ttau_trun_off_disc = .true. and/or ttaur_turn_off_jet
        ! in your parameter file.
       if (grid%geometry=="ttauri" .and. ttau_turn_off_disc) then
          write(*,*) " "
          write(*,*) "Turning off the alpha disc read in from file, but keeping"
          write(*,*) "the magnetorsphere alive."
          if (.not.ttau_disc_on)   then
              ! we need to create the disc parameter object 
             call new(ttauri_disc, dble(TTauriDiskRin*TTauriRstar/1.0e10), 1.5d3*100.d0, &
                  dble(TTauriMstar/100.0), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, dble(TTauriMstar/mSol))
          end if
          call turn_off_disc(grid%octreeroot, grid, ttauri_disc)
       end if

       if (grid%geometry=="ttauri" .and. ttau_turn_off_jet) then
          write(*,*) " "
          write(*,*) "Turning off the jet read in from file, but keeping"
          write(*,*) "the magnetorsphere alive."
          if (.not.ttau_jet_on) then
              ! the parameter file has to created here.
             call new(ttauri_jet,  dble(TTauriRouter/1.0e10), dble(amrgridsize), &
                  JET_theta_j, dble(TTauriMstar/mSol), JET_Mdot, JET_a_param, JET_b_param,  &
                  JET_Vbase, JET_Vinf, JET_beta, JET_gamma, limitscalar, JET_T)
          end if
          call turn_off_jet(grid%octreeroot, grid, ttauri_jet)
       end if

       if (grid%geometry=="ttauri" .and. ttau_turn_off_acc) then
          write(*,*) " "
          write(*,*) "Turning off the magnetosphere read in from file."
          write(*,*) " "
           ! using the routine in amr_mod.f90
          call turn_off_magnetosphere(grid%octreeroot, grid, dble(TTauriRouter/1.0e10))
       end if

       if (myRankIsZero) then   !-----------------------------------------------
          if (writePhasePops) then
             write(tempChar,'(i3.3)') nStartPhase
             phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
             call writeAMRgrid(phasePopFilename,writeFileFormatted,grid)
          else if (writePops) then
             call writeAMRgrid(popFilename,writeFileFormatted,grid)
          end if
       end if  ! (myRankIsZero) --------------------------------------------------------

       call torus_mpi_barrier

    else  ! not reading a population file

       amrGridCentre = octalVector(amrGridCentreX,amrGridCentreY,amrGridCentreZ)
       call writeInfo("Starting initial set up of adaptive grid...", TRIVIAL)
       
       select case (geometry)
       case("cluster")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, sphData, young_cluster, nDustType)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, sphData, young_cluster)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
!           call fill_in_empty_octals(young_cluster,grid%octreeRoot,sphData)
          call estimateRhoOfEmpty(grid, grid%octreeRoot, sphData)	
           !Removing the cells within 10^14 cm from the stars.
          call remove_too_close_cells(young_cluster,grid%octreeRoot,1.0d4)

       case("wr104")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, sphData, young_cluster, nDustType)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,sphData)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
          call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
          do
             gridConverged = .true.
             call myScaleSmooth(smoothFactor, grid%octreeRoot, grid, &
                  gridConverged,  inheritProps = .false., interpProps = .false.)
             if (gridConverged) exit
          end do
          call writeInfo("...grid smoothing complete", TRIVIAL)
          
       case("starburst")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, sphData, young_cluster, nDustType)
          gridconverged = .false.
          do while(.not.gridconverged) 
             call splitGridFractal(grid%octreeRoot, real(100.*mHydrogen), 0.1, grid, gridconverged)
          enddo
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

       case("magstream")
          nSource = 1
          allocate(source(1:1))
          source(1)%luminosity = grid%lCore
          source(1)%radius = ttaurirStar/1.d10
          source(1)%teff = 4000.
          source(1)%position = VECTOR(0.,0.,0.)
          call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(100.), dble(2.e8), 200)
          call normalizedSpectrum(source(1)%spectrum)
          call buildSphere(grid%starPos1, grid%rCore, source(1)%surface, 400, contFluxFile)
          nu =1.d15
!           call createMagStreamSurface(source(1)%surface, grid, nu, coreContinuumFlux, fAccretion)
!           call testSurface(source(1)%surface)

          call readStreams(thisStream,nStreams,"stream.dat")

          call allocateStream(bigStream, nStreams*200)
          do i = 1, nStreams
             do j = 1, thisStream(i)%nsamples
                bigStream%nSamples = bigStream%nSamples + 1
                bigStream%rho(bigStream%nSamples) = thisStream(i)%rho(j)
                bigStream%temperature(bigStream%nSamples) = thisStream(i)%temperature(j)
                bigStream%position(bigStream%nSamples) = thisStream(i)%position(j)
                bigStream%velocity(bigStream%nSamples) = thisStream(i)%velocity(j)
             enddo
          enddo

          write(*,*) "stream read . now refining"
          write(*,*) "splitting on stream"
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, sphData, young_cluster, nDustType, &
               stream=bigStream)

          if (doTuning) call tune(6, "Magstream grid construction.") ! start a stopwatch


          write(*,*) "Bigstream has ",bigStream%nsamples, " samples" 
          call splitGridOnStream3(grid%octreeRoot,  grid, bigstream) 
          call freeStream(bigStream)

          nOctals = 0
          nVoxels = 0
          call countVoxels(grid%octreeRoot,nOctals,nVoxels)
          write(message,*) "Adaptive grid currently contains: ",nOctals," octals"
          call writeInfo(message, TRIVIAL)
          write(message,*) "                                : ",nVoxels," unique voxels"
          call writeInfo(message, TRIVIAL)
           
          if (writeoutput) then
             fac = 110.*degtorad
             do i = 1, 1
                ang = real(i-1)/40. * twoPi
                octVec = OCTALVECTOR(sin(fac)*cos(ang), sin(fac)*sin(ang), cos(fac))
                write(thisdevice,'(a,i3.3,a)') "col",i,".ps/vcps"
                call columnDensityPlot(grid, source, nsource, octVec, thisDevice)
             enddo
          endif


          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
          if (doTuning) call tune(6, "Magstream grid construction.") ! stop a stopwatch
!           call plotAMRthreeDMovie(grid, source, nsource)

       case DEFAULT
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, &
               sphData, young_cluster, nDustType, romData=romData)
          call writeInfo("First octal initialized.", TRIVIAL)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,romData=romData)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

           ! This section is getting rather long. Maybe this should be done in 
           ! wrapper subroutine in amr_mod.f90.
          if (geometry=="ttauri") then 
              ! Finding the total mass in the accretion flow
             mass_accretion_old = 0.0d0
             call TTauri_accretion_mass(grid%octreeRoot, grid, mass_accretion_old)
             write(message,*) "Total mass in accretion flow is ",  mass_accretion_old, "[g]"
             call writeInfo(message,FORINFO)

             if (ttau_fuzzy_edge) then
                call writeInfo("Fuzzy edge for TTauri accreation stream will be used.",FORINFO)
                call ttauri_fuzzy_edge(grid%octreeRoot)
                
                 ! We now have to correct the density of the accretion flow 
                 ! due to the fuzziness introduced above.
                mass_accretion_new = 0.0d0
                call TTauri_accretion_mass(grid%octreeRoot, grid, mass_accretion_new)
                write(message,*) "Total mass in accretion flow (after fuzzy edge) is ",  mass_accretion_new, "[g]"
                call writeInfo(message, FORINFO)
                mass_scale = mass_accretion_old/mass_accretion_new
                write(message,*) "Scaling the density by ", mass_scale , "[g]"
                call writeInfo(message, FORINFO)

                call TTauri_accretion_scale_density(grid%octreeRoot, grid, mass_scale) 
                 ! Check to see if it has scaled correctly.
                mass_accretion_new = 0.0d0
                call TTauri_accretion_mass(grid%octreeRoot, grid, mass_accretion_new)                 
                if (writeoutput) write(*,*) "Total mass in accretion flow after rescaling: ", mass_accretion_new, "[g]"
                if (writeoutput) write(*,*) " "  
             end if

             if (ttau_discwind_on) then
                if (writeoutput) write(*,*) " "
                if (writeoutput) write(*,*) "Adding TTauri disc wind to magnetosphere model... "
                call add_discwind(grid%octreeRoot, grid, ttauri_discwind, limitscalar2)
                if (writeoutput) write(*,*) " "
             elseif (ttau_jet_on) then
                if (writeoutput) write(*,*) " "
                if (writeoutput) write(*,*) "Adding TTauri jets to magnetosphere model... "
                call add_jet(grid%octreeRoot, grid, ttauri_jet)
                call finish_grid_jet(grid%octreeroot, ttauri_jet)
                if (writeoutput) write(*,*) " "
             end if


             if (grid%geometry=="ttauri" .and. ttau_turn_off_disc) then
                write(*,*) " "
                write(*,*) "Turning off the alpha disc read in from file, but keeping"
                write(*,*) "the magnetorsphere alive."
                if (.not.ttau_disc_on)   then
                   ! we need to create the disc parameter object 
                    call new(ttauri_disc, dble(TTauriDiskRin*TTauriRstar/1.0e10), 1.5d3*100.d0, &
                         dble(TTauriMstar/100.0), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, dble(TTauriMstar/mSol))
                 end if
                 call turn_off_disc(grid%octreeroot, grid, ttauri_disc)
              end if
              
              if (grid%geometry=="ttauri" .and. ttau_turn_off_jet) then
                 write(*,*) " "
                 write(*,*) "Turning off the jet read in from file, but keeping"
                 write(*,*) "the magnetorsphere alive."
                 if (.not.ttau_jet_on) then
                    ! the parameter file has to created here.
                    call new(ttauri_jet,  dble(TTauriRouter/1.0e10), dble(amrgridsize), &
                         JET_theta_j, dble(TTauriMstar/mSol), JET_Mdot, JET_a_param, JET_b_param,  &
                         JET_Vbase, JET_Vinf, JET_beta, JET_gamma, limitscalar, JET_T)
                 end if
                 call turn_off_jet(grid%octreeroot, grid, ttauri_jet)
              end if


           end if  ! gemoetry == "ttaruri"

           ! 
           if (doSmoothGrid) then
	      call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
              do
                 gridConverged = .true.
                 ! The following is Tim's replacement for soomthAMRgrid.
                 call myScaleSmooth(smoothFactor, grid%octreeRoot, grid, &
                      gridConverged,  inheritProps = .false., &
		      interpProps = .false.,  &
                      sphData=sphData, stellar_cluster=young_cluster, romData=romData)
                 if (gridConverged) exit
              end do
              call writeInfo("...grid smoothing complete", TRIVIAL)
           endif

           ! Smooth the grid with respect to optical depth, if requested

           if (doSmoothGridTau.and.mie) then
              call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
              do j = iSmoothLam, iSmoothlam !nLambda, 5
                 do
                    gridConverged = .true.
                    call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, inheritProps = .false., interpProps = .false.)
                    if (gridConverged) exit
                 end do
              enddo
              call writeInfo("...grid smoothing complete", TRIVIAL)

              ! The tau smoothing may result in large differences in the size
              ! of neighbouring octals, so we smooth the grid again.
              
              if (doSmoothgrid) then
                 call writeInfo("Smoothing adaptive grid structure (again)...", TRIVIAL)
                 do
                    gridConverged = .true.
                 ! The following is Tim's replacement for soomthAMRgrid.
                    call myScaleSmooth(smoothFactor, grid%octreeRoot, grid, &
                         gridConverged,  inheritProps = .false., interpProps = .false., &
                         sphData=sphData, stellar_cluster=young_cluster, romData=romData)
                    if (gridConverged) exit
                 end do
                 call writeInfo("...grid smoothing complete", TRIVIAL)
              endif
           end if
        end select
   
!        call writeInfo("Unrefining optically thin cells...", TRIVIAL)
!        gridconverged = .false.
!        do while(.not.gridconverged)
!           gridconverged = .true.
!           call unrefineThinCells(grid%octreeRoot, grid, ismoothlam, gridconverged)
!        end do
!        call writeInfo("Done.", TRIVIAL)

!        call writeInfo("Unrefining optically thick cells...", TRIVIAL)
!        gridconverged = .false.
!        do while(.not.gridconverged)
!           gridconverged = .true.
!           call unrefineThickCells(grid%octreeRoot, grid, ismoothlam, gridconverged)
!        end do
!        call writeInfo("Done.", TRIVIAL)


        nOctals = 0
        nVoxels = 0
        call countVoxels(grid%octreeRoot,nOctals,nVoxels)
        write(message,*) "Adaptive grid contains: ",nOctals," octals"
        call writeInfo(message, TRIVIAL)
        write(message,*) "                      : ",nVoxels," unique voxels"
        call writeInfo(message, TRIVIAL)
        grid%nOctals = nOctals

        call writeInfo("Calling routines to finalize the grid variables...",TRIVIAL)
        gridConverged = .false.
     
        do
           call finishGrid(grid%octreeRoot,grid,gridConverged,romData=romData)
           if (gridConverged) exit
        end do        

        call writeInfo("...final adaptive grid configuration complete",TRIVIAL)

!        call grid_info(grid, "*")
        if (myRankIsZero) call grid_info(grid, "info_grid.dat")

        if (writeoutput) call writeAMRgrid("after_creation.grid",.false.,grid)

        if ((geometry == "shakara").and.(nDustType>1)) then
           call fillDustShakara(grid, grid%octreeRoot)
        endif

        if ((geometry == "whitney").and.(nDustType ==4)) then
           call	fillDustWhitney(grid, grid%octreeRoot)	
        endif

        if (((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry=="warpeddisc")).and.(nDustType>1)) then
           call fillDustUniform(grid, grid%octreeRoot)
        endif
        
!        if (variableDustSublimation) then
!           call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac)
!        endif


        
!     ! plotting column density (new routine in grid_mod.f90)
!     if (grid%geometry == "luc_cir3d") then 
!        call plot_column_density(grid, plane_for_plot,  "column_density.ps/vcps", &
!            nmarker, xmarker, ymarker, zmarker, &
!            width_3rd_dim, show_value_3rd_dim,val_3rd_dim)
!     end if


        if (geometry == "ttauri" .or. geometry == "luc_cir3d" .or. &
             geometry == "cmfgen" .or. geometry == "romanova".or.geometry == "magstream") then
           nu = cSpeed / (lamLine * angstromtocm)
           call contread(contFluxFile, nu, coreContinuumFlux)
           call buildSphere(grid%starPos1, grid%rCore, starSurface, 400, contFluxFile)
           if (geometry == "ttauri") then
              call createTTauriSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion) 
           elseif (geometry == "magstream") then
!              call createMagStreamSurface(source(1)%surface, grid, nu, coreContinuumFlux, fAccretion)
!              call testSurface(source(1)%surface)

           elseif (geometry == "romanova") then
              call createTTauriSurface2(starSurface, grid, romData, nu, coreContinuumFlux,fAccretion) 
           else
              call createSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion) 
           end if
!           call testSurface(starSurface)
        endif

        if (lineEmission.and.(.not.cmf)) then
           !  calculate the statistical equilibrium (and hence the emissivities 
           !  and the opacities) for all of the subcells in an
           !  adaptive octal grid.
           !  Using a routine in stateq_mod module.
           if (writeoutput) write(*,*) "Calling statistical equilibrium routines..."
           if (doTuning) call tune(6, "amrStateq") ! start a stopwatch  
           if (grid%geometry=="ttauri" .or. grid%geometry(1:8)=="windtest" .or. &
                grid%geometry(1:9)=="luc_cir3d" .or. geometry == "romanova") then
              call amrStateq(grid, newContFluxFile,lte, nLower, nUpper, &
                             starSurface, recalcPrevious=.false.)

              call torus_mpi_barrier('waiting for other amrStatEq calls to return...')

              ! the 'lte' setting is intended so that a model with infall
              ! enhancement is first set up with LTE values. Only after we have 
              ! done the first time-dependent change so we calculate non-LTE 
              ! values.

              
              ! RK added the followings ------------------------------------------
              !  adding alpha disc to the grid
              if (ttau_disc_on) then
                 if (writeoutput) write(*,*) "Adding an accretion disc ... "
                 call add_alpha_disc(grid%octreeroot, grid, ttauri_disc)
                 ! finding the photon-dust interaction x-secions
                 call mieCrossSection(sigmaExt0, sigmaAbs0, sigmaSca0,  &
                      aMin(1), aMax(1), a0(1), qDist(1), pDist(1), grainType(1), &
                      ngrain, X_grain, grainname, lamLine)
                 if (writeoutput) write(*,*) " "
                 if (writeoutput) write(*,*) "Photon-dust cross section at lambda = ", lamLine, " [A]"
                 if (writeoutput) write(*,*) "    sigma(tot) = ", sigmaExt0, " [cm^2]"
                 if (writeoutput) write(*,*) "    sigma(abs) = ", sigmaAbs0, " [cm^2]"
                 if (writeoutput) write(*,*) "    sigma(sca) = ", sigmaSca0, " [cm^2]"
                 if (writeoutput) write(*,*) " "
                 call finish_grid(grid%octreeroot, grid, ttauri_disc, 1.0, &
                      sigmaAbs0, sigmaSca0, meanDustParticleMass)
                 !
                 !
!                 allocate(lambda(1:maxTau),tauExt(1:maxTau),tauAbs(1:maxTau),tauSca(1:maxTau),&
!                          contTau(1:maxTau,1:nLambda), linePhotonalbedo(1:maxtau))
!                 call integratePath(gridUsesAMR, VoigtProf, &
!                      lambdatau,  lamLine, OCTALVECTOR(1.0d-20,1.0d-20,1.0d-20), OCTALVECTOR(150.0,0,-200.), &
!                      OCTALVECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
!                      tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .false. , &
!                      lamStart, lamEnd, nLambda, contTau, hitCore, thinLine, lineResAbs, .false., &
!                      .false., nUpper, nLower, 0., 0., 0., junk,&
!                      sampleFreq,intPathError, useInterp, grid%Rstar1, coolStarPosition)              
!                 write(*,*) " ----  Vertical optical depth at r = 0.1 AU is ", tauExt(ntau)
!                 write(*,*) " ----                           at lambda =  ", lambdatau
!                 write(*,*) "Rescaling scattering and absorption coefficients ... "
!                 call finish_grid(grid%octreeroot, grid, ttauri_disc, tauExt(ntau)/100.0)
!                 call integratePath(gridUsesAMR, VoigtProf, &
!                      lambdatau,  lamLine, OCTALVECTOR(1.,1.,1.), OCTALVECTOR(75.0,0,-50.), &
!                      OCTALVECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
!                      tauSca, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .false. , &
!                      lamStart, lamEnd, nLambda, contTau, hitCore, thinLine, lineResAbs, .false., &
!                      .false., nUpper, nLower, 0., 0., 0., junk,&
!                      sampleFreq,intPathError, useInterp, grid%Rstar1, coolStarPosition)              
!                 write(*,*) " ---- New Optical depth  at r = 0.05 AU    is ", tauExt(ntau)
!                 write(*,*) " ----                           at lambda =  ", lambdatau
!                 deallocate(lambda,tauExt,tauAbs,tauSca,contTau,linePhotonalbedo)
!                 write(*,*) ".. Finished rescaling scattering and absorption coefficients ... "
                 if (writeoutput) write(*,*) ".. Finished adding an acreation disc ... "
              end if
              !-------------------------------------------------------------------

           elseif (geometry == "cmfgen") then
              ! simply map CMFGEN opacity data to the AMR grid
              call map_cmfgen_opacities(grid)

           else

              call amrStateq(grid, newContFluxFile, lte, nLower, nUpper, &
                      starSurface, recalcPrevious=.false., ion_name=ion_name, ion_frac=ion_frac)
!              if (ttau_disc_on) then
!                 ! amrStateq will have messed up the disc, so we reset those cells
!                 call finish_grid(grid%octreeroot, grid, ttauri_disc, 1.0, sigmaAbs0, sigmaSca0)
!              end if

              call torus_mpi_barrier('waiting for other amrStatEq calls to return...')

           end if ! (grid%geometry=="ttaur" .or. ...)

           if (doTuning) call tune(6, "amrStateq") ! stop a stopwatch  
           if (writeoutput) write(*,*) "... statistical equilibrium routines complete"



        end if ! (lineEmission)

        !
        ! cleaning up unused memory here ....
        if ((geometry(1:7) == "cluster").or.(geometry(1:5)=="wr104")) then
           ! using the routine in sph_data_class
           call kill(sphData)
           ! using the routine in amr_mod.f90

          if (myRankIsZero) call delete_particle_lists(grid%octreeRoot)

        elseif (geometry == "luc_cir3d") then
           call deallocate_zeus_data()

! *** This part is commented out for now since romData is needed 
! *** for creating the surface elements later..... 
! *** This can be changed later (RK)
!
!        elseif (geometry == "romanova") then
!           call kill(romData)
!
        end if


     end if ! (readPops .or. readPhasePops)


  if (myRankIsZero)  then
     ! writing pop files after statEq routines
     if (writePhasePops) then
        write(tempChar,'(i3.3)') nStartPhase
        phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
        call writeAMRgrid(phasePopFilename,writeFileFormatted,grid)
     end if
     if (writePops) then
        call writeAMRgrid(popFilename,writeFileFormatted,grid)
     end if

  end if
  call torus_mpi_barrier ! sync here
  
  if (doTuning) call tune(6, "AMR grid construction.") ! stop a stopwatch

!     if (geometry == "ttauri".and.VoigtProf) then
!        write(*,*) "Setting biases for stark broadening..."
!        call setBiasChil(grid%octreeroot)
!        write(*,*) "done..."
!     endif

  if (geometry(1:6) == "ttauri" .and. myRankIsZero) then
     call writeHartmannValues(grid,'hartmann_logNH')
     call writeHartmannValues(grid,'hartmann_logNe')
     call writeHartmannValues(grid,'hartmann_temperature')
     call writeHartmannValues(grid,'hartmann_velPol')
     call writeHartmannValues(grid,'hartmann_velAz')
     call writeHartmannValues(grid,'hartmann_line')
     !call writeHartmannValues(grid,'hartmann_Nelectron')
     call writeHartmannValues(grid,'hartmann_Nlevel2')
     call writeHartmannValues(grid,'hartmann_NH')
     !call writeHartmannValues(grid,'hartmann_departCoeff')
     call writeHartmannValues(grid,'hartmann_N')
  end if

end subroutine amr_grid_setup

!-----------------------------------------------------------------------------------------------------------------------

! Plot values of variables on AMR grid. Extracted from torusMain by D. Acreman 

subroutine do_amr_plots

  !     if (grid%geometry == "jets"  .or. &
!          grid%geometry == "ttauri"  .or.  grid%geometry == "testamr" ) then
!        call draw_cells_on_density(grid, plane_for_plot, "cells_on_density.ps/vcps")
!        !     call draw_cells_on_density(grid, plane_for_plot, device)
!     end if
     

  ! Do some preparation for the arrays used in plot_AMR_* which will be used later
  if (grid%geometry(1:7) == "cluster") then
     nmarker = get_nstar(young_cluster)
     allocate(xMarker(1:nMarker))
     allocate(yMarker(1:nMarker))
     allocate(zMarker(1:nMarker))
     do i = 1, nmarker
        a_star = get_a_star(young_cluster, i)
        xmarker(i)= a_star%position%x
        ymarker(i)= a_star%position%y
        zmarker(i)= a_star%position%z
     end do
  else
     nmarker = 0
     ALLOCATE(xmarker(nmarker), ymarker(nmarker), zmarker(nmarker))
  end if
  width_3rd_dim = amrGridSize

  if ( plot_maps .and. myRankIsZero ) then

  ! Plot desired AMR grid value here... This is more generalized
  ! version of fancyAmrPlot.
  !
  ! See grid_mod.f90 for details.
  !
  ! Plotting some grid values
     call plot_AMR_values(grid, "rho", plane_for_plot, val_3rd_dim, &
          "rho_grid",.true., .true., nmarker, xmarker, ymarker, zmarker, &
          width_3rd_dim, show_value_3rd_dim, suffix="default", index=num_calls, &
          fixValMin=sph_rho_min, fixValMax=sph_rho_max, useFixedRange=.true. )
     call plot_AMR_values(grid, "rho", plane_for_plot, val_3rd_dim, &
          "rho_zoom",.true., .true., nmarker, xmarker, ymarker, zmarker, &
          width_3rd_dim, show_value_3rd_dim, boxfac=zoomFactor, suffix="default", index=num_calls, &
          fixValMin=sph_rho_min, fixValMax=sph_rho_max, useFixedRange=.true. )
     if ((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry=="warpeddisc")) then
        call plot_AMR_values(grid, "rho", plane_for_plot, val_3rd_dim, &
             "rho_ultrazoom",.true., .true., nmarker, xmarker, ymarker, zmarker, &
             width_3rd_dim, show_value_3rd_dim, boxfac=0.005, suffix="default", index=num_calls)
     end if
     call plot_AMR_planes(grid, "rho", plane_for_plot, 3, "rho", .true., .false., &
          nmarker, xmarker, ymarker, zmarker, show_value_3rd_dim)
     if (grid%lineEmission) then
        call plot_AMR_values(grid, "Vx", plane_for_plot, val_3rd_dim,  &
             "Vx", .false., .false.,  &
             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim, suffix="default")
        call plot_AMR_values(grid, "Vy", plane_for_plot, val_3rd_dim, &
             "Vy", .false., .false., &
             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim, suffix="default")
        call plot_AMR_values(grid, "Vz", plane_for_plot, val_3rd_dim, &
             "Vz", .false., .false., &
             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim, suffix="default")
     end if
  !        call plot_AMR_values(grid, "etaCont", plane_for_plot, val_3rd_dim,  &
  !             "etacont.ps/vcps", .true., .false.,  &
  !             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim)
  !        call plot_AMR_values(grid, "etaCont", plane_for_plot, val_3rd_dim,  &
  !             "etacont_zoom.ps/vcps", .true., .false.,  &
  !             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, &
  !             show_value_3rd_dim, boxfac=0.0004)
     call plot_AMR_values(grid, "temperature", plane_for_plot, val_3rd_dim, &
          "temperature", .true., .false., nmarker, xmarker, ymarker, zmarker, &
          width_3rd_dim, show_value_3rd_dim, suffix="default", index=num_calls, &
          fixValMin=sph_tem_min, fixValMax=sph_tem_max, useFixedRange=.false. )
  !        call plot_AMR_values(grid, "temperature", "x-y", 0., &
  !             "temperature2.ps/vcps", .true., .false., &
  !             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim)
  !        if (lineEmission) then
  !           call plot_AMR_values(grid, "etaLine", plane_for_plot, val_3rd_dim,  &
  !                "etaline.ps/vcps", .true., .false., & 
  !                nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim)
  !           call plot_AMR_values(grid, "chiLine", plane_for_plot, val_3rd_dim,  &
  !                "chiline.ps/vcps", .true., .false., &
  !                nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim)
  !           call plot_AMR_values(grid, "dV_dR", plane_for_plot, val_3rd_dim,  &
  !                "dv_dr.ps/vcps", .true., .false., &
  !                nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim)
  !        end if
  
     if (mie) then
        do i = 1, nDustType
           write(message,'(a,i1)') "dusttype",i
           write(filename,'(a,i1)') "dusttype",i
           call plot_AMR_values(grid, message, "x-z", 0., &
                trim(filename),.false., .false., nmarker, xmarker, ymarker, zmarker, &
                width_3rd_dim, show_value_3rd_dim, boxfac=zoomfactor, suffix="default", index=num_calls)
        enddo
     endif

  !     ! plotting column density (new routine in grid_mod.f90)
  !     if (grid%geometry(1:7) == "cluster") then 
  !       call plot_column_density(grid, plane_for_plot,  "column_density.ps/vcps", &
  !             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim,val_3rd_dim)
  !     end if

  endif

end subroutine do_amr_plots

!-----------------------------------------------------------------------------------------------------------------------
subroutine set_up_sources

  integer      :: nstar
  real(double) :: d1, d2, massRatio
  real         :: tmp

  ! set up the sources
  nSource = 0

  select case(geometry)
    case("testamr","benchmark")
       nSource = 1
       allocate(source(1:1))
       source(1)%luminosity = grid%lCore
       source(1)%radius = grid%rCore
       source(1)%teff = teff
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(lamstart), dble(lamEnd), nLambda)
       call normalizedSpectrum(source(1)%spectrum)

    case("magstream")
       nSource = 1
       if (.not.allocated(source)) then
          allocate(source(1:1))
          source(1)%luminosity = grid%lCore
          source(1)%radius = ttaurirStar/1.d10
          source(1)%teff = 4000.
          source(1)%position = VECTOR(0.,0.,0.)
          call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(100.), dble(2.e8), 200)
          call normalizedSpectrum(source(1)%spectrum)
          call buildSphere(grid%starPos1, grid%rCore, source(1)%surface, 400, contFluxFile)
       endif
       nu =1.d15
       call genericAccretionSurface(source(1)%surface, grid, nu, coreContinuumFlux, fAccretion)
       call testSurface(source(1)%surface)

    case("melvin")
       nSource = 1
       teff = 30000.
       allocate(source(1:1))
       source(1)%luminosity = fourPi * (10.*rsol)*(10.*rsol) * stefanBoltz * teff**4
       source(1)%radius = 10.*rSol / 1.e10
       source(1)%teff = teff
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(lamStart), dble(lamEnd),nLambda)
       call normalizedSpectrum(source(1)%spectrum)
       rstar = source(1)%radius

    case("whitney")
       nSource = 1
       allocate(source(1:1))
       source(1)%luminosity = fourPi * rStellar**2 * stefanBoltz * teff**4
       source(1)%radius = rStellar / 1.e10
       source(1)%teff = teff
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(lamStart), dble(lamEnd),nLambda)
       call normalizedSpectrum(source(1)%spectrum)
       rstar = source(1)%radius

    case("toruslogo")
       nSource = 1
       allocate(source(1:1))
       source(1)%luminosity = lsol
       source(1)%radius = rSol / 1.e10
       source(1)%teff = 6000.
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(6000),  dble(lamStart), dble(lamEnd),nLambda)
       call normalizedSpectrum(source(1)%spectrum)
       rstar = source(1)%radius


    case("starburst")
       call random_seed(size=iSize)
       allocate(iSeed(1:iSize))
       iseed = 2343245
       call random_seed(put=iSeed)
       deallocate(iSeed)

       allocate(source(1:10000))
       call createSources(nSource, source, "instantaneous", 1.d6, 1.d3, 1.d0)
       call random_seed


    case("wr104")
       nSource = 2

       allocate(source(1:nSource)) 
       source(1)%teff = 30000.  ! o star
       source(1)%radius = 10. * rSol / 1.e10
       source(1)%position = (VECTOR(1.,0.,0.)*real(autocm))/1.e10
       source(1)%luminosity = fourPi * stefanBoltz * (10.*rSol)**2.0 * (source(1)%teff)**4
       call readSpectrum(source(1)%spectrum, "ostar.flx", ok)
       call normalizedSpectrum(source(1)%spectrum)

       source(2)%teff = 40000.             ! wr star 
       source(2)%radius = 20. * rSol / 1.e10
       source(2)%position = (VECTOR(-1.,0.,0.)*real(autocm))/1.e10
       source(2)%luminosity = 0.5 * source(1)%luminosity
       call readSpectrum(source(2)%spectrum, "wr.flx", ok)
       call normalizedSpectrum(source(2)%spectrum)

    case("lexington","fractal")
       nSource = 1
       allocate(source(1:nSource)) 
       source(1)%teff = 40000.  ! o star
       source(1)%radius = 18.67 * rSol / 1.e10
       source(1)%position = VECTOR(0.,0.,0.)
       source(1)%luminosity = fourPi * stefanBoltz * (source(1)%radius*1.e10)**2.0 * (source(1)%teff)**4
       if (writeoutput) write(*,*) "Lexington source: ",source(1)%luminosity/1.e37
       fac = 1.e8*cspeed/5.d16
       call fillSpectrumBB(source(1)%spectrum, dble(source(1)%teff), fac, 1000.d4,1000)
       call normalizedSpectrum(source(1)%spectrum)

    case("symbiotic")
       nSource = 2
       allocate(source(1:nSource)) 
       source(1)%teff = 2500.  
       source(1)%radius = 100. * rSol / 1.e10
       source(1)%position = VECTOR(-250.*rSol/1.e10,0.,0.)
       source(1)%luminosity = fourPi * stefanBoltz * (source(1)%radius*1.e10)**2.0 * (source(1)%teff)**4
       call fillSpectrumBB(source(1)%spectrum, dble(source(1)%teff), 10.d0, 1000.d4,1000)
       call normalizedSpectrum(source(1)%spectrum)

       source(2)%teff = 150000.  
       source(2)%radius = 5000.e5 / 1.e10 ! 5000 km
       source(2)%position = VECTOR(250.*rSol/1.e10,0.,0.)
       source(2)%luminosity = fourPi * stefanBoltz * (source(2)%radius*1.e10)**2.0 * (source(2)%teff)**4
       call fillSpectrumBB(source(2)%spectrum, dble(source(2)%teff), 10.d0, 1000.d4,1000)
       call normalizedSpectrum(source(2)%spectrum)
       write(*,*) "!!!!!!! light ratio",source(1)%luminosity/source(2)%luminosity



    case ("cluster")
       ! Extract some info from cluster object.
       nstar = get_nstar(young_cluster)  ! number of stars in the cluster
       nSource =  n_stars_in_octal(young_cluster, grid%octreeRoot)
       
       ! copy the star over in the array.
       ! This is ugly. Maybe lucyRadiativeEquilibriumAMR should be changed to take
       ! an cluster_class object as an input variable in future.
       ALLOCATE(source(nSource))
       
       ! Restricting the source to be within the root cell (in case the root cell is 
       ! is smaller than the sph model space!
       j = 0 
       do i = 1, nstar
          a_star = get_a_star(young_cluster, i)
          ! using a function in source_mod
          if ( source_within_octal(a_star, grid%octreeRoot) ) then
             j = j+1
             source(j) = a_star
          end if
       end do

    case("shakara","clumpydisc","wrshell","warpeddisc","iras04158")
       nSource = 1
       allocate(source(1:1))
       source(1)%radius = grid%rCore
       source(1)%teff = teff   
       source(1)%position = VECTOR(0.,0.,0.)
       tmp = source(1)%radius * 1.e10  ! [cm]
       source(1)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(1)%teff)**4
       if (contFluxfile .eq. "blackbody") then
          call fillSpectrumBB(source(1)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd),nLambda, lamArray=xArray)
       else
          call buildSphere(o2s(source(1)%position), real(source(1)%radius), source(1)%surface, 400, contFluxFile)
          call readSpectrum(source(1)%spectrum, contfluxfile, ok)
       endif
       call normalizedSpectrum(source(1)%spectrum)

    case("gammavel")
       nSource = 2
       allocate(source(1:2))


       massRatio = mass1/mass2

       d1 = binarySep * (1./(massRatio+1.))
       d2 = binarySep - d1

! source 1 is the WR star
       source(1)%radius = rstar1
       source(1)%teff = teff1
       source(1)%position = VECTOR(0.,0.,-d1)
       tmp = source(1)%radius * 1.e10  ! [cm]
       source(1)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(1)%teff)**4

       call buildSphere(o2s(source(1)%position), real(source(1)%radius), source(1)%surface, 400, contFluxFile1)


       if (contFluxfile1 .eq. "blackbody") then
          call fillSpectrumBB(source(1)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd), nLambda, lamArray=xArray)
       else
          call readSpectrum(source(1)%spectrum, contfluxfile1, ok)
       endif
       call normalizedSpectrum(source(1)%spectrum)

! source 2 is the WR star
       source(2)%radius = rstar2
       source(2)%teff = teff2
       source(2)%position = VECTOR(0.,0.,d2)
       tmp = source(2)%radius * 1.e10  ! [cm]
       source(2)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(2)%teff)**4

       call buildSphere(o2s(source(2)%position), real(source(2)%radius), source(2)%surface, 400, contFluxFile2)

       if (contFluxfile2 .eq. "blackbody") then
          call fillSpectrumBB(source(2)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd), nLambda)
       else
          call readSpectrum(source(2)%spectrum, contfluxfile2, ok)
       endif
       call normalizedSpectrum(source(2)%spectrum)

    ! chris (26/05/04)
    case ("ppdisk")
       nSource = 1
       allocate(source(1:nSource))
!       source(1)%teff = 5780.
       source(1)%teff = Teff
!       source(1)%luminosity = 3.83d+32 
!       source(1)%luminosity = 10.**((4.75-Mbol)/2.5) * lSol
! Baraffe, et al. use 4.64 for solar bolometric magnitude
       source(1)%luminosity = 10.**((4.64-Mbol)/2.5) * lSol
!       source(1)%radius = 6.96d+0
       source(1)%radius = sqrt(source(1)%luminosity/(4.*pi*stefanBoltz*Teff**4))/1.d10
       source(1)%position = VECTOR(0.,0.,0.)

       call fillSpectrumBB(source(1)%spectrum, dble(teff), dble(lamStart), dble(lamEnd), nLambda)
       call normalizedSpectrum(source(1)%spectrum)

    case ("planetgap")
       nSource = 1
       allocate(source(1:nSource))
       source(1)%teff = Teff
       source(1)%luminosity = fourPi*(rcore*1.e10)**2 * stefanBoltz*teff**4
       source(1)%radius = rCore
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(teff), dble(lamStart), dble(lamEnd), nLambda)
       call normalizedSpectrum(source(1)%spectrum)


    case default
       ! Allocating the source array with size =0 to avoid, non-allocated array passed problem
       ! in subroutine initPhoton...
       if (.not. allocated(source)) ALLOCATE(source(0))
       nSource = 0  ! This must be zero!
       
    end select

  if (nSource > 0) then
     call randomSource(source, nSource, i, xArray, nLambda, initialize=.true.)  
     call writeInfo("Sources set up.",TRIVIAL)
  endif


end subroutine set_up_sources

!-----------------------------------------------------------------------------------------------------------------------

end subroutine torus

end module torus_mod

!-----------------------------------------------------------------------------------------------------------------------
subroutine choose_view ( geometry, nPhase, distortionType, doRaman, &
       rotateview, rotateDirection, tiltView)

  use input_variables, only: forceRotate, forceNoRotate

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

  if (forceRotate)   rotateView = .true.
  if (forceNoRotate) rotateView = .false.

end subroutine choose_view

!-------------------------------------------------------------------------------
! Name:    update_sph_temperature
! Purpose: Update the temperatures of the SPH particles using the torus temperature field
!          Each process works on its own subset of the total number of particles.
! Author:  D. Acreman, November 2007

  subroutine update_sph_temperature (b_idim, b_npart, b_iphase, b_xyzmh, sphData, grid, b_temp)

    USE vector_mod, only:     octalVector
    USE amr_mod, only:        amrGridValues
    USE gridtype_mod, only:   gridType
    USE sph_data_class, only: sph_data, get_udist
    USE messages_mod

    implicit none

#ifdef MPI
! MPI specific variables
    include 'mpif.h'
    real    :: mpi_max_deltaT, mpi_sum_deltaT
    integer :: ierr, mpi_iiigas
#endif

! Arguments 
    integer, intent(in)   :: b_idim, b_npart
    integer*1, intent(in) :: b_iphase(b_idim)
    real*8, intent(in)    :: b_xyzmh(5,b_idim)
    real*8, intent(inout) :: b_temp(b_idim)
    type(sph_data), intent(in) :: sphData
    type(GRIDTYPE), intent(in) :: grid

! Local variables
    character(len=80) :: message
    integer      :: iiigas, i
    real*8       :: xgas, ygas, zgas
    real(double) :: sphDistFac
    real         :: tgas
    real         :: deltaT, sum_deltaT, mean_deltaT, max_deltaT
    type(OCTALVECTOR) :: octVec

! Begin executable statements

! 1. Calculate conversion from sph distance units to torus distance units
       sphDistfac  = get_udist(sphData) ! [cm]
       sphDistfac = sphDistfac / 1.e10  ! to torus units

! 2. Update particle temperatures and calculate some statistics 
       iiigas = 0
       sum_deltaT = 0.0
       max_deltaT = 0.0
       do i=1, b_npart
          if (b_iphase(i) == 0) then
             iiigas = iiigas + 1 
! 2.1 Determine position of gas particle on the amr grid and extract the temperature value
             xgas = b_xyzmh(1,i) * sphDistFac
             ygas = b_xyzmh(2,i) * sphDistFac
             zgas = b_xyzmh(3,i) * sphDistFac
             octVec = OCTALVECTOR(xgas,ygas,zgas)
             call amrGridValues(grid%octreeRoot, octVec, temperature=tgas, grid=grid)
! 2.2 Calculate statistics of temperature change
             deltaT     = tgas - b_temp(iiigas)
             sum_deltaT = sum_deltaT + deltaT
             max_deltaT = MAX(max_deltaT, deltaT)
! 2.3 Update the gas particle temperature to pass back to sph code
             b_temp(iiigas) = tgas
          endif
       enddo

! 3. Calculate mean temperature change and perform MPI communication if required
!    MPI processes which did not do the SPH step will not have any particles to update
#ifdef MPI
! Calculate global values
       call MPI_ALLREDUCE( iiigas,     mpi_iiigas,     1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
       call MPI_ALLREDUCE( max_deltaT, mpi_max_deltaT, 1, MPI_REAL,    MPI_MAX, MPI_COMM_WORLD, ierr )
       call MPI_ALLREDUCE( sum_deltaT, mpi_sum_deltaT, 1, MPI_REAL,    MPI_SUM, MPI_COMM_WORLD, ierr )
! Update values to be output
       iiigas      = mpi_iiigas
       max_deltaT  = mpi_max_deltaT
       mean_deltaT = mpi_sum_deltaT / real(mpi_iiigas)
#else
       mean_deltaT = sum_deltaT / real(iiigas)
#endif

! 4. Write out the temperature change statistics
       write(message, *) "Number of particles= ", iiigas
       call writeInfo(message, FORINFO)
       write(message, *) "Maximum temperature change= ", max_deltaT
       call writeInfo(message, FORINFO)
       write(message, *) "Mean temperature change=  ", mean_deltaT
       call writeInfo(message, FORINFO)

  end subroutine update_sph_temperature

!-------------------------------------------------------------------------------  

!!! vim:set filetype=fortran :                                !!!  
!!! otherwise vim won't recognize a file with the suffix .raw !!!
