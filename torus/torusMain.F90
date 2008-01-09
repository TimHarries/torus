!
!            T   O   R   U  S
!
! torus is a general purpose monte-carlo radiative transfer code used to
! produce polarized spectra of hot and cool stellar winds
!

! written by tjh
! broken by chris

! v1.0 on 16/09/99

! raman scattering stuff added 3/3/2000

! OMP parallelization calls added 1/7/2001

! TJH: 19/7/02  Jorick's spotty star stuff added
! NHS: 02/9/02  adaptive mesh code merged
! RK : 12/02/03 Parallelized major loops in lucyStateEquilibriumAMR, stateqAMR, torusMain.
! RK : 07/07/05 Added the observed flux solver in oberver's frame. 

#ifdef SPH
 subroutine torus(b_idim,b_npart,b_nactive,b_xyzmh,b_rho,b_iphase, &
                  b_nptmass,b_listpm,b_udist,b_umass, &
                  b_utime,b_time,b_gaspartmass,b_temp)
#else
program torus
#endif


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

  real, allocatable :: kappaAbs(:), kappaSca(:), kappaExt(:)
  real  :: sigmaExt0, sigmaAbs0, sigmaSca0  ! cross section at the line centre
  real :: dlambda, thisTau

  ! variables to do with dust

  integer :: itestlam, ismoothlam
  integer, parameter :: nXmie = 20, nMuMie = 20
  type(PHASEMATRIX), pointer :: miePhase(:,:, :)

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

  integer,parameter :: IterLucy=5   ! Minimum number of iterations of LUCYRADIATIVEEQUILIBRIUMAMR
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
#ifdef SPH
  type(isochrone)       :: isochrone_data
  logical, parameter    :: ll_sphFromFile = .false.   ! Read in sph data from a file?
  integer, intent(in)   :: b_idim,b_npart,b_nactive,b_nptmass
  integer*1, intent(in) :: b_iphase(b_idim)
  integer, intent(in)   :: b_listpm(b_nptmass)
  real*8, intent(in)    :: b_xyzmh(5,b_idim)
  real*4, intent(in)    :: b_rho(b_idim)
  real*8, intent(in)    :: b_udist, b_umass, b_utime, b_time, b_gaspartmass
  real*8, intent(inout) :: b_temp(b_idim)
  integer :: ngaspart
  logical, parameter :: ll_sph = .true.
#else
  logical, parameter :: ll_sph = .false.
#endif
  integer, save :: num_calls = 0

  type(modelatom), allocatable :: thisAtom(:)

  type(STREAMTYPE) :: thisStream(2000), bigStream
  integer :: nStreams
  
  integer :: nFromEnv
  logical :: photonFromEnvelope

  integer :: nRBBTrans
  integer :: indexRBBTrans(1000), indexAtom(1000)

  integer :: iInner_beg, iInner_end ! beginning and end of the innerPhotonLoop index.

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
  !  initialize the system for running MPI
  if (.not. ll_sph) call MPI_INIT(ierr) 

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
  if (ll_sph) num_calls = num_calls + 1

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

  if (molecular .and. geometry == "molebench") then
     call readbenchmarkMolecule(co, "hco_benchmark.dat")
!     call readMolecule(co, "hco_plus.mol")
  endif

  if (geometry == "molefil") then

     call readMolecule(co, "c18o.dat")
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
#ifdef SPH
  if (geometry == "cluster") then

     ! HSC
     if  (ll_sphFromFile) then

        ! read in the sph data from a file
        call read_sph_data(sphData, "sph.dat")
        call read_stellar_disc_data(sphData, "stellar_disc.dat")

     else

     ! The total number of gas particles is the total number of active particles less the number of point masses.
        ngaspart = b_nactive-b_nptmass
        call init_sph_data2(sphData, b_udist, b_umass, b_utime, ngaspart, b_time, b_nptmass, &
             b_gaspartmass, b_npart, b_idim, b_iphase, b_xyzmh, b_rho, b_temp)
! Communicate particle data. Non-mpi case has a stubbed routine. 
        call gather_sph_data(sphData)

     end if

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
     if (myRankIsZero .and. ll_sphFromFile) call find_inclinations(sphData, 0.0d0, 0.0d0, 1.0d0, "inclinations_z.dat")

  end if
#endif

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
  allocate(kappaExt(1:nLambda))
  allocate(kappaAbs(1:nLambda))
  allocate(kappaSca(1:nLambda))
  allocate(statArray(1:nLambda))
  allocate(sourceSpectrum(1:nLambda))
  allocate(sourceSpectrum2(1:nLambda))

  grid%doRaman = doRaman

  statArray = 0.
  sourceSpectrum = 1.
  sourceSpectrum2 = 1.


  if (plezModelOn) then
     call readTioCrossSection(xArray, nLambda, kappaExt, kappaAbs, kappaSca)
  endif


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

     if (molecular.and.readmol) then
        call readAMRgrid(molfilenamein,.false.,grid)
        goto 667
     endif

     if (cmf.and.readLucy) then
        call readAMRgrid("atom_tmp.grid",.false.,grid)
     else if (photoionization.and.readlucy) then
        continue
     else
        ! Set up the AMR grid
        call amr_grid_setup

        if (geometry(1:7) .eq. "testamr") call testamr
        
        if (grid%geometry(1:8) == 'windtest') call windtest

        ! Write the information on the grid to file using a routine in grid_mod.f90
        if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")

        ! Plotting the various values stored in the AMR grid.
        if ( plot_maps .and. myRankIsZero ) call do_amr_plots
     end if

  !=================================================================
  else ! grid is not adaptive ======================================     
  !=================================================================

     ! fill up the grid with the appropriate opacities
     call fill_opacities_noamr( ok )
     if ( .not. ok ) goto 666

     ! Distort grid if required
     call distort_noamr

     if ( enhance ) then
        nVec = 100
        nPhi = 360
        allocate(distortionVec(nVec))
        timeStart = 0.
        timeEnd = 12.*3600.
        dTime = (timeEnd - timeStart)/real(nPhase)
        call initInfallEnhancement(distortionVec, nVec, nPhi, infallParticleMass)
        do i = 1, nStartPhase-1
           call infallEnhancment(grid, distortionVec, nVec, nPhi, dTime, &
                .false., infallParticleMass, alreadyDoneInfall)
        enddo
     endif

     ! Write the information on the grid to file using a routine in grid_mod.f90
     if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")

  !=================================================================
  endif ! (gridusesAMR)
  !=================================================================

  ! The source spectrum is normally a black body
!  if (.not.grid%lineEmission) then
!     stot = 0.
!     do i = 1, nLambda
!        sourceSpectrum(i) = bLambda(dble(xArray(i)), dble(teff))
!        stot = stot + sourceSpectrum(i)
!     enddo
!     sourceSpectrum  = sourceSpectrum / stot
!  endif

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
        if (geometry == "shakara") then

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

     if (molecular) then
        if (.not.readlucy) call  molecularLoop(grid, co)
        call calculateMoleculeSpectrum(grid, co)
        call createDataCube(cube, grid, OCTALVECTOR(0.d0, 1.d0, 0.d0), co, 1)
        if (myRankIsZero) call plotDataCube(cube, 'cube.ps/vcps')
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


#ifdef SPH
! If this is an SPH run then finish here ---------------------------------------
     call update_sph_temperature (b_idim, b_npart, b_iphase, b_xyzmh, sphData, grid, b_temp)
     goto 666
#endif

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


!  if (grid%geometry == "shakara") then
!     call defineDiffusionZone(grid, .false., .false.)
!  endif

  
  !
  ! Setting the emission bias.
  !
  if (useBias .and. .not. formalsol) then
     write(*,*) "Setting emission bias..."
     select case(grid%geometry)
        case("testamr")
           call setBiasAMR(grid%octreeRoot, grid)
        case("cluster")
           ! Computing the emission bias based on the optical depth at
           ! 2 microns. 
           call assign_emission_bias(grid%octreeroot, grid, 20000.0, xArray, nLambda)
        case("ttauri")
           call set_bias_ttauri(grid%octreeRoot, grid, lamline, outVec)
        case("cmfgen")
           call set_bias_cmfgen(grid%octreeRoot, grid, lamline)
        case("shakara","ppdisk","planetgap")
           if (forcedWavelength) then
              call locate(grid%lamArray, nLambda, usePhotonWavelength, itestlam)
              call set_bias_shakara(grid%octreeRoot, grid, ilam=itestlam, ross=.false.)
           else
              call locate(grid%lamArray, nLambda, usePhotonWavelength, itestlam)
              call set_bias_shakara(grid%octreeRoot, grid, ilam=1, ross=.true.)
!              call setBiasOnTau(grid)
           endif
!        case("ppdisk")
!           call locate(grid%lamArray, nLambda, usePhotonWavelength, itestlam)
!           call setDiscPhotosphereBias(grid, itestlam)
        case("whitney")
           call set_bias_whitney(grid%octreeRoot, grid)
        case("warpeddisc")
           call set_bias_rosseland(grid%octreeRoot, grid)
!	   call set_bias_radius(grid%octreeRoot, 3)
        case DEFAULT
           continue
     end select
     write(*,*) "Done."
  end if

  if (myRankIsZero .and. geometry == "shakara") call polardump(grid)
  
  ! initialize the blobs if required

    if (nBlobs > 0) then
  
       allocate(blobs(1:maxBlobs))
       if (freshBlobs) then
          do i = 1 , maxBlobs
             blobs(i)%inUse = .false.
          enddo
          if (writeoutput) write(*,'(a)') "Running blobs for five days..."
   
          ! now we run a few of days worth of blobs
  
          t1 = 0.
          t2 = 240. * 60. * 60.
          dTime = (t2 - t1) / 100.
          do i = 1, 100
             call addNewBlobs(grid, maxBlobs, blobs, blobTime, dTime, &
                              nCurrent, blobContrast)
             call moveBlobs(maxBlobs, blobs, 0., dTime, grid)
             if (writeoutput) write(*,*) i,nCurrent
          enddo
  
          j = 0 
          do i = 1, maxBlobs
             if (blobs(i)%inUse) j = j + 1
          enddo
  
          if (writeoutput) write(*,'(a,i3)') "Number of blobs after 5 days: ",j
  
  
          dTime = 0.
          if (nPhase /= 1) then
             dTime = (timeEnd - timeStart) / real(nPhase-1)
          endif
  
          open(50,file="files.lis",status="unknown",form="formatted")
  
          if (writeoutput) write(*,'(a)') "Running blobs and writing configuration files"
          do i = 1, nPhase
  
             grid%timeNow = timeStart + (timeEnd-timeStart)*real(i-1)/real(nPhase-1)	    
             write(specFile,'(a,i3.3,a)') trim(outfile),i,".dat"
             write(50,*) specFile(1:30), grid%timeNow
  
             call addNewBlobs(grid, maxBlobs, blobs, blobTime, dTime, nCurrent, blobContrast)
             call moveBlobs(maxBlobs, blobs, 0., dTime, grid)
             write(filename,"(a,i3.3,a)") "run",i,".blob"
             call writeBlobs(filename, maxBlobs, blobs)
  
          enddo
  
          close(50)
  
       endif
    endif

  if (myRankIsZero .and. plotVelocity) then
     write(*,*) ' Plotting velocity vectors...'
     call plotVelocityVectors(grid, "1.gif/GIF")
  endif


  if (nPhase /= 1) then
     dTime = (timeEnd - timeStart) / real(nPhase-1)
  endif

  if (doRaman) then
     ramanSourceVelocity = (ramVel/cSpeed)*VECTOR(0.,-1., 0.)
     if (writeoutput) write(*,'(a,f7.1,a)') "Raman source speed: ",modulus(ramanSourceVelocity)*cSpeed/1.e5, " km/s"
  endif


  !
  ! Redefining the rotation axis here.
  !
  if (geometry == "ttauri") then
     rotationAxis = grid%diskNormal
     if (writeoutput) write(*,*) "rotation axis",rotationAxis
  endif

  if (geometry == "luc_cir3d" .or. geometry == "cmfgen") then
     rotationAxis = VECTOR(0.,0.,1.)
     if (writeoutput) write(*,*) "rotation axis",rotationAxis
  endif

  if (geometry == "romanova" ) then
     rotationAxis = grid%diskNormal
     if (writeoutput) write(*,*) "rotation axis",rotationAxis
  endif

  if (geometry == "donati") then
     rotationAxis = rotateX(rotationAxis, -dipoleOffset)
  endif

  if (geometry == "jets") then
     rotationAxis = VECTOR(1., 0., 0.)
  end if 

  if (distortionType == "spiral") then
     rotationAxis = rotateX(rotationAxis, dipoleOffset)
  endif




  normToRotation = rotationAxis .cross. yAxis
  call normalize(normToRotation)

  tempVec  = arbitraryRotate(rotationAxis, inclination, normToRotation)
  originalViewVec =  (-1.)*tempVec


  if (modulus(thermalElectronVelocity(10000.)) == 0.) then
     if (writeoutput) write(*,'(a)') "THERMAL ELECTRON BROADENING IS OFF!!!!!!!!!!!!!!!!!!!!!!!!!!"
  endif


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


  
  if (geometry(1:6) == "ttauri" .or. geometry(1:9) == "luc_cir3d" .or. &
       geometry == "cmfgen" .or. geometry == "romanova") then
    call emptySurface(starSurface)
  end if

  if (forceRotate)   rotateView = .true.
  if (forceNoRotate) rotateView = .false.
  

! From here we do multiple runs if required !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  tiltView = .false.

  phaseLoop: do iPhase = nStartPhase, nEndPhase
          
     grid%timeNow = phaseTime * real(iPhase-1)

     call unixTimes(cpuTime, startTime)

     viewVec = originalViewVec
     outVec = (-1.)*viewVec
     thisVec = viewVec


     ! we rotate the view by an appropriate amount


     if (rotateView.and.(nPhase /= 1)) then
        if (writeoutput) write(*,'(a)') "Rotating view..."
        phi = -rotateDirection * twoPi * real(iPhase-1)/real(nPhase-1)
        viewVec =  arbitraryRotate(thisVec, phi, rotationAxis)
        outVec = (-1.) * viewVec
        if (writeoutput) write(*,'(a,f5.2,a,f5.2,a,f5.2,a)') "View vector: (",viewVec%x,",", &
             viewVec%y,",",viewVec%z,")"
     endif

     if (phaseOffset /= 0.) then
        if (writeoutput) write(*,'(a,f5.2)') "Rotating by phase offset: ",phaseOffset
        phi = twoPi * phaseOffset
        viewVec =  arbitraryRotate(viewVec, phi, rotationAxis)
        outVec = (-1.) * viewVec
     endif


     if (tiltView.and.(nPhase /= 1)) then
        phi = -twoPi * real(iPhase-1)/real(nPhase-1)
        viewVec =  arbitraryRotate(thisVec, phi, xAxis)
        outVec = (-1.) * viewVec
     endif


     ! refill the grids 
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!              
     if (.not.plezModelOn .and. .not. gridUsesAMR) then
        select case(geometry)
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
                     curtains, dipoleOffset, nLower, nUpper, theta1, theta2)
                     
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
           if (writeoutput) write(*,*) "! Unrecognised grid geometry: ",trim(geometry)
           goto 666
        end select
     endif


     if (geometry(1:6) == "ttauri" .or. geometry(1:9) == "luc_cir3d" .or. &
          geometry =="cmfgen" .or. geometry == "romanova") then      
        ! Nu must be set again here since it is not assigned when the population/grid 
        ! file is read from a file!  (RK changed here.)
        nu = cSpeed / (lamLine * angstromtocm)
        call contread(contFluxFile, nu, coreContinuumFlux)
        call buildSphere(grid%starPos1, grid%rCore, starSurface, 400, contFluxFile)
        if (geometry == "ttauri") then
           call createTTauriSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion) 
        elseif (geometry == "romanova") then
           call createTTauriSurface2(starSurface, grid, romData, nu, coreContinuumFlux,fAccretion) 
        else
           call createSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion)            
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
             call amrStateq(grid, newContFluxFile, lte, nLower, nUpper,  &
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

     

     if (geometry .eq. "rolf") then
        if (gridUsesAMR) then
           !      is this correct vvvvvvvvvvvvvvvvvvvvvv ?
           secondSourcePosition = grid%octreeRoot%centre
        else
           secondSourcePosition = VECTOR(grid%xAxis(grid%nx/2), &
             grid%yAxis(grid%ny/2), &
             grid%zAxis(grid%nz/2))
        end if
     endif

     call  createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)

     ! if we are producing a movie then plot this phase


     allocate(lambda(1:maxTau))
     allocate(tauExt(1:maxTau))
     allocate(tauAbs(1:maxTau))
     allocate(tauSca(1:maxTau))
     allocate(linePhotonalbedo(1:maxTau))
     allocate(contTau(1:maxTau,1:nLambda))
     allocate(contWeightArray(1:nLambda))


     if (myRankIsZero .and. movie) then
        if (device(2:3) == "xs") then
           call plotGrid(grid, viewVec,  opaqueCore, &
                device,contTau, foreground, background, coolStarPosition, firstPlot)
        endif
!
!! COMMENTED OUT FOR DEBUG
!!        write(filename,"(a,i3.3,a)") "frame",iPhase,".gif/gif"
        write(filename,"(a,i3.3,a)") "frame",iPhase,".ps/vcps"
        call plotGrid(grid, viewVec,  opaqueCore, &
             filename, contTau, foreground, background, coolStarPosition, firstPlot)
     endif

     if (geometry == "wr104") then
        write(*,*) "stripping dust!!!!!!!!!!!"
        call stripDustAway(grid%octreeRoot, 1.d-20, 4.d4)
     endif

     call plot_AMR_values(grid, "etaCont", plane_for_plot, val_3rd_dim,  &
          "etacont_zoom.ps/vcps", .true., .false.,  &
          nmarker, xmarker, ymarker, zmarker, width_3rd_dim, &
          show_value_3rd_dim, boxfac=zoomfactor)


     if (writeoutput) write(*,*) " "
     if (writeoutput) write(*,'(a)') "Some basic model parameters"
     if (writeoutput) write(*,'(a)') "---------------------------"
     if (writeoutput) write(*,*) " "

     ! THE FOLLOWING STATEMENT IS DANGEROUS. ACCORDING TO INPUT_MOD.F90 
     ! NOT ALL THE GEOMETRY HAS RCORE VALUES, AND SAME UNITS!
     ! THIS SHOULD BE DONE IN INITAMRGRID ROUTINE AS SOME GEOMETRY HAS
     ! DONE SO.
     if (geometry /= "cmfgen") then	
        grid%rStar1 = rcore/1.e10
     end if
     !
     ! Performs various optical depth tests here...
     !

     if ( (.not. hydrodynamics) .and. (.not. formalsol) .and. myRankIsZero ) then
        call test_optical_depth(gridUsesAMR, VoigtProf, &
             amrGridCentre, sphericityTest,  &
             outVec, lambdatau,  lambdatau, grid, thin_disc_on, opaqueCore, lamStart, lamEnd,  &
             thinLine, lineResAbs, nUpper, nLower, sampleFreq, useinterp, grid%Rstar1, coolStarPosition, maxTau, nSource, source)
      end if

      call torus_mpi_barrier

     weightLinePhoton = 0.
     weightContPhoton = 1.
     weightPhoto = 1.

     if (mie .and. (.not. useDust)) then

        call calcContinuumEmissivityLucy(grid, grid%octreeRoot , nlambda, xArray)

        call computeProbDist(grid, totLineEmission, &
             totDustContinuumEmission,lamline, .false.)
        totDustContinuumEmission = totdustContinuumEmission 
        lcore = grid%lCore
        if (nSource > 0) then
           write(*,*) "!!!!! summing source over",xarray(1),xarray(nlambda)
           lCore = sumSourceLuminosity(source, nsource, xArray(1), xArray(nLambda))
        endif

        totEnvelopeEmission = totDustContinuumEmission
        chanceDust = totDustContinuumEmission/(totDustContinuumEmission+lCore/1.e30)
        if (writeoutput) write(*,*) "totdustemission",totdustcontinuumemission
        if (writeoutput) write(*,'(a,f7.2)') "Chance of continuum emission from dust: ",chanceDust

        weightDust = chanceDust / probDust
        weightPhoto = (1. - chanceDust) / (1. - probDust)

        if (writeoutput) write(*,*) "WeightDust",weightDust
        if (writeoutput) write(*,*) "WeightPhoto",weightPhoto
        if (writeoutput) write(*,*) "core + envelope luminosity",lCore+totEnvelopeEmission*1.d30
        energyPerPhoton =  ((lCore + totEnvelopeEmission*1.d30) / dble(nPhotons))/1.d20
        if (writeoutput) write(*,*) "Energy per photon: ", energyPerPhoton

     endif

     



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
        if (writeoutput) write(*,*) "Total Raman Line Emission: ",totLineEmission
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


!        write(*,'(a)') "mu axis"
!        do i = 1, grid%nMu
!           write(*,*) i,acos(grid%muAxis(i))*radtodeg, &
!                             grid%muProbDistLine(grid%nr/2,i)
!        enddo
!        write(*,'(a)') "Phi axis"
!        do i = 1, grid%nPhi
!           write(*,*) i,grid%phiAxis(i)*radtodeg, grid%phiProbDistLine(grid%nr/2,grid%nmu/2,i)
!        enddo



        if (geometry == "donati") totWindContinuumEmission = 0.

        if (geometry == "binary") totWindContinuumEmission = 0.

        if (grid%resonanceLine) totWindContinuumEmission = 0.

        !if (geometry == "ttauri" .or. geometry == "windtest") then
        if (geometry == "windtest") then
           totWindContinuumEmission = 0.
           if (writeoutput) write(*,'(a)') "! Wind continuum emission switched off."
        endif

        if (lineOff) then
           totLineEmission = 0.
        endif

        nu = cSpeed / (lamLine * angstromtocm)
        if (writeoutput) write(*,'(a,e12.3)') "Line emission: ",totLineEmission
        select case(geometry)
           case("binary")
              call contread(contFluxFile, nu, totCoreContinuumEmission1)
              call contread(contFluxFile2, nu, totCoreContinuumEmission2)
           case("puls")
              totCoreContinuumEmission = pi * blackBody(0.77 * tEff, lamLine)
           case DEFAULT
              call contread(contFluxFile, nu, coreContinuumFlux)
              totCoreContinuumEmission = coreContinuumFlux
        end select


        ! reads the intrinsic core absorption profile for core-halo models
        ! Note: The profile should be normaised (to contiuum) and the x axis shoudl be in Angstrome.

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
        if (writeoutput) write(*,'(a,e12.3)') "Wind cont emission: ",totWindContinuumEmission


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

           if (writeoutput) write(*,*) "Binary luminosity ratio (p/s): ",totCoreContinuumEmission1/totCoreContinuumEmission2
        endif
           
        if (geometry == "ttauri") then
           write (*,'(a,e12.3)') 'T Tauri star: continuum emission: ',totCoreContinuumEmission
           fAccretion = fAccretion * (nuStart-nuEnd)
           write (*,'(a,e12.3)') 'T Tauri accretion: continuum emission: ',fAccretion
           write (*,'(a,f12.3)') 'Accretion continuum / stellar continuum: ',fAccretion/totCoreContinuumEmission
           totCoreContinuumEmission = totCoreContinuumEmission + fAccretion
           write (*,'(a,e12.3)') 'T Tauri total core continuum emission: ',totCoreContinuumEmission
           ! chanceHotRing = fAccretion/totCoreContinuumEmission ! no longer used
        endif

        if (geometry == "luc_cir3d" .or. geometry == "cmfgen" .or. geometry == "romanova") then
           write (*,'(a,e12.3)') 'Star: continuum emission: ',totCoreContinuumEmission
           fAccretion = fAccretion * (nuStart-nuEnd)
           write (*,'(a,e12.3)') 'Accretion: continuum emission: ',fAccretion
           write (*,'(a,f12.3)') 'Accretion continuum / stellar continuum: ',fAccretion/totCoreContinuumEmission
           totCoreContinuumEmission = totCoreContinuumEmission + fAccretion
        endif


        
        chanceSpot = 0.
        if ((geometry == "disk").and.(nSpot > 0)) then
           chanceSpot = fSpot * blackBody(tSpot, 6562.8) / &
                ((1.-fSpot)*blackBody(tEff, 6562.8) + &
                (fSpot * blackBody(tSpot, 6562.8)))
           if (writeoutput) write(*,'(a,f5.3)') "Spot chance at 6563A: ",chanceSpot
        endif




        if (writeoutput) write(*,'(a,e12.3)') "Total Core Continuum Emission: ",totCorecontinuumEmission


        totContinuumEmission = totCoreContinuumEmission + totWindContinuumEmission

        if (writeoutput) write(*,'(a,e12.3)') "Total Continuum emission: ", totContinuumEmission


        if ((totContinuumEmission + totLineEmission) /= 0.) then
           chanceLine = totLineEmission/(totContinuumEmission + totLineEmission)
           chanceContinuum = totContinuumEmission / &
                (totContinuumEmission + totLineEmission)
           grid%chanceWindOverTotalContinuum = totWindContinuumEmission &
                / max(1.d-30,totContinuumEmission)
        else
           chanceLine =0.
           chanceContinuum = 1.
           grid%chanceWindOverTotalContinuum = 0.
        endif



        if (writeoutput) write(*,'(a,e12.3)') "Chance continuum emission in wind: ", grid%chanceWindOverTotalContinuum

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


        energyPerPhoton =  (totLineEmission + totContinuumEmission) / dble(nPhotons)

        write(*,*) "Energy per photon: ", energyPerPhoton


        write(*,*) "chance line",chanceline


        ! compute the probability distributions


     endif
     if (.not. formalsol) then
        deallocate(lambda)
        deallocate(tauSca)
        deallocate(tauExt)
        deallocate(tauAbs)
        deallocate(contTau)
        deallocate(linePhotonAlbedo)
        deallocate(contWeightArray)
     end if

     nContPhotons = nint(probContPhoton * real(nPhotons) / real(nOuterLoop))

     !
     ! here we may loop over different inclinations
     !
     
  incLoop: do iInclination = 1, nInclination, 1 
       ! NB This loop is not indented in the code!
     
     if (nInclination > 1) then  
        if (allocated(inclinations)) then
           inclination = max(inclinations(iInclination),1.e-4)
        else
           ! Now the inclinations are equally spaced in  cos(incination).
           ! Asuuming that firstInclination < lastInclination <= pi/2
           cos_inc_first = COS(firstInclination)
           cos_inc_last = COS(lastInclination)
           d_cos_inc = (cos_inc_first - cos_inc_last)/ REAL(nInclination-1)
           cos_inc = cos_inc_first - d_cos_inc * REAL(iInclination-1)
           inclination = max(ACOS(cos_inc),1.e-4)


!          inclination = firstInclination + REAL(iInclination-1) * &
!                          (lastInclination-firstInclination)/REAL(nInclination-1)
        end if


       write(*,*) " "
       write(*,*) "Inclination = ",inclination*radToDeg,' degrees'
       write(*,*) " "

       if (iPhase == nStartPhase .and. iInclination == 1) originalOutFile = outFile
         
       write(tempChar,'(i3.3)') NINT(inclination*radToDeg)
       outFile = trim(originalOutFile)//'_inc'//TRIM(tempChar)
       
       viewVec%x = 0.
       viewVec%y = -sin(inclination)
       viewVec%z = -cos(inclination)
       outVec = (-1.)*viewVec
       thisVec = viewVec

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

     end if

     write(*,*) "Viewing vector: ",viewVec

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
        nImage = get_nfilter(filters)

        ! Allocate the image array
        if (allocated(obsImageSet)) deallocate(obsImageSet)
        allocate(obsImageSet(nImage))
        

        ! Initializing the images ...

        if (setImageSize == 0.) then
           if (grid%adaptive) then
              if (amr2d) then
                 imageSize = 2.*grid%octreeRoot%subcellSize
              else
                 imageSize = grid%octreeRoot%subcellSize          
              endif
           else if (grid%cartesian) then
              imageSize = grid%xAxis(grid%nx)-grid%xAxis(1)
           else
              imagesize = 2.*grid%rAxis(grid%nr)
           endif
        else
           imageSize = setImageSize / 1.e10
        endif

        if (imageInArcsec) then
           imagesize = objectdistance * (imageSizeInArcsec / 3600.) * degtorad / 1.e10
        endif


!     if (molecular) then
!        if (writemol) call  molecularLoop(grid, co)
!        call calculateMoleculeSpectrum(grid, co)
!        call createDataCube(cube, grid, OCTALVECTOR(0.d0, 1.d0, 0.d0), co, 1)
!        if (myRankIsZero) call plotDataCube(cube, 'cube.ps/vcps')
!        stop
!     endif


        ! Change the image size accoring to the scale 
        imageSize = imageSize * imagescale



        ! Initializing the images ...
        do i = 1, nImage           
           if (grid%cartesian) then
              obsImageSet(i) = initImage(npix, npix, imageSize, imageSize, vmin, vmax)
              if (doRaman) then
                 o6image(1) = initImage(npix, npix, imageSize, imageSize, vmin, vmax)
              endif
           else if (grid%adaptive) then
              obsImageSet(i) = initImage(npix, npix, imageSize, imageSize, vmin, vmax)
           else   
              select case (geometry)
              case("disk")
                 obsImageSet(i) = initImage(npix, npix, 2.*grid%rAxis(grid%nr), 2.*grid%rAxis(grid%nr),vMin, vMax)
              case("flared")
                 obsImageSet(i) = initImage(npix, npix, 8.*grid%rAxis(1),8.*grid%rAxis(1), vMin, vMax)
              case DEFAULT
                 obsImageSet(i) = initImage(npix,npix,  min(5.*grid%rAxis(1),grid%rAxis(grid%nr)), &
                      min(5.*grid%rAxis(1),grid%rAxis(grid%nr)), vMin, vMax)
              end select
           endif
        end do
        
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

     if (writeoutput) then
        write(*,*) " "
        write(*,'(a)') "Run-time messages"
        write(*,'(a)') "-----------------"
        write(*,*) " "
     endif
     
     ntot = 0

     ! now we loop 10 times over a tenth of the photons - this will be used
     ! to help compute the errors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call unixTimes(cpuTime, startTime)



     !====================================================================================
     ! Perform formal integration
     !===================================================================================
     if ((geometry == "ttauri" .or. geometry == "romanova") &
          .and. formalsol) then
        if (doTuning) call tune(6, "One Formal Sol") ! start a stopwatch  
        write(*,*) " "
        write(*,*) "Started formal integration of flux...."
        if (.not. ALLOCATED(flux)) ALLOCATE(flux(grid%nlambda))
        if (nInclination > 1) then
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
             OCTALVECTOR(0.0d0, 0.0d0, 0.0d0), starSurface, &
             form_nr_wind, form_nr_acc, form_nr_core, form_nphi,  &
             s2o(outVec), objectDistance, grid, sampleFreq, opaqueCore,  &
             flux, grid%lamArray, grid%nlambda, obsfluxfile, &
             thin_disc_on, (ttau_disc_on .and. .not. ttau_turn_off_disc), &
             (ttau_jet_on .and. .not. ttau_turn_off_jet), ttau_discwind_on, npix, &
             (ttau_disc_on .and. .not. ttau_turn_off_disc), ttauri_disc, thinLine, lineOff, &
             do_pos_disp) 
        write(*,*) "...Finished  formal integration of flux...."
        write(*,*) " "
        if (doTuning) call tune(6, "One Formal Sol") ! stop a stopwatch
        ! jumps to the end of inclination loops  (this should be replaced with while loop later)
        goto 777  
     end if

     !  These should be zero-ed for each viewing angle!
     tooFewSamples = 0 
     boundaryProbs = 0 
     negativeOpacity = 0 

     if (doTuning) call tune(6, "All Photon Loops")  ! Start a stopwatch
     
     call random_seed


     call randomSource(source, nSource, i, xArray, nLambda, initialize=.true.)  

     if (mie) then
        nInnerLoop = nPhotons / nOuterLoop
     endif


     outerPhotonLoop: do iOuterLoop = 1, nOuterLoop

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(i, contPhoton, contWindPhoton, r, nScat) &
!$OMP PRIVATE(thisPhoton, directionalWeight) &
!$OMP PRIVATE( ilambda, sourceSpectrum, ok) &
!$OMP PRIVATE(hitCore, junk, thisLam, j, obs_weight, thisVel) &
!$OMP PRIVATE(i1, i2, i3, t1, t2, t3, vray, vovercsqr, fac, observedLambda) &
!$OMP PRIVATE(t, rHat, islit, fac1, fac2, fac3, obsPhoton, r1, r2, thisTau) &
!$OMP PRIVATE(escaped, currentScat, absorbed, dlambda, thisChi, thisSca) &
!$OMP PRIVATE(albedo, tempPhoton, redRegion, thrustar, ramanWeight) &
!$OMP PRIVATE(outPhoton,intPathError) &
!$OMP PRIVATE(nTau, escProb, spotPhoton) &
!$OMP PRIVATE(lambda, tauExt, tauSca, tauAbs, contTau, contWeightArray) &
!$OMP PRIVATE(rHatinStar, positionOc, linePhotonalbedo, dopShift, lineResAbs, tau_bnd) &
!$OMP SHARED(grid) &
!$OMP SHARED(meanr_Cont, wtot_cont,meanr_line,wtot_line, ntot) &
!$OMP SHARED(nContPhotons, nPhotons, lineEmission, lamLine, nLambda) &
!$OMP SHARED(weightLinePhoton, flatSpec, vRot, secondSource, secondSourcePosition) &
!$OMP SHARED(ramanSourceVelocity, vO6, doRaman) &
!$OMP SHARED(weightContPhoton, useBias, pencilBeam ,outVec)&
!$OMP SHARED(opaqueCore, lamStart, lamEnd, thinLine, rStar, coolStarPosition) &
!$OMP SHARED(viewVec, o6xArray,o6yArray, rotationAxis, o6image, screened) &
!$OMP SHARED(xArray, yArray, statArray, stokesImage, obsImageSet, doPvimage) &
!$OMP SHARED(nSlit, pvimage, gridDistance, meanr0_line, wtot0_line) &
!$OMP SHARED(sourceSpectrum2, meanr0_cont,wtot0_cont, maxScat, mie) &
!$OMP SHARED(miePhase, zeroVec, theta1, theta2, chanceHotRing) &
!$OMP SHARED(nSpot, chanceSpot, thetaSpot, phiSpot, fSpot, chanceDust) &
!$OMP SHARED(narrowBandImage, vMin, vMax, gridUsesAMR) &
!$OMP SHARED(sampleFreq, useInterp, photLine,tooFewSamples,boundaryProbs) &
!$OMP SHARED(probDust, WeightDust, WeightPhoto, source, nsource) &
!$OMP SHARED(energyPerPhoton, filters, nUpper, nLower, nImage) &
!$OMP SHARED(negativeOpacity, iInner_beg, iInner_end) &
!$OMP SHARED(curtains, starSurface, VoigtProf, nDustType, ttauri_disc, ttau_disc_on) &
!$OMP SHARED(forcedWavelength, usePhotonWavelength, thin_disc_on, forceFirstScat)


        if (mie) then

           iLambdaPhoton = iOuterLoop

           call calcContinuumEmissivityLucyMono(grid, grid%octreeRoot , nlambda, xArray, iLambdaPhoton)

           call setBiasOnTau(grid, iLambdaPhoton)

           call computeProbDist(grid, totLineEmission, &
                totDustContinuumEmission,lamline, .false.)
           totDustContinuumEmission = totdustContinuumEmission 
           lcore = grid%lCore
           if (nSource > 0) then              
              lCore = sumSourceLuminosityMonochromatic(source, nsource, dble(xArray(iLambdaPhoton)))
           endif
           
           totEnvelopeEmission = totDustContinuumEmission
           chanceDust = totDustContinuumEmission/(totDustContinuumEmission+lCore/1.e30)
           if (writeoutput) write(*,*) "totdustemission",totdustcontinuumemission
           if (writeoutput) write(*,*) "totcontemission",lcore/1.e30
           if (writeoutput) write(*,'(a,f9.7)') "Chance of continuum emission from dust: ",chanceDust
           
!           weightDust = chanceDust / probDust
!           weightPhoto = (1. - chanceDust) / (1. - probDust)

           probDust = chanceDust
           weightDust = 1.
           weightPhoto = 1.

           energyPerPhoton =  (totDustContinuumEmission*1.d30 + lCore)/1.d20/dble(nInnerLoop)
           if (writeoutput) write(*,*) "WeightDust",weightDust
           if (writeoutput) write(*,*) "WeightPhoto",weightPhoto
           if (writeoutput) write(*,*) "core + envelope luminosity",lCore+totEnvelopeEmission*1.d30
           if (writeoutput) write(*,*) "Energy per photon: ", energyPerPhoton

        endif

        if (doTuning) call tune(6, "One Outer Photon Loop") ! Start a stop watch

        ! default inner loop indices
        iInner_beg = 1
        iInner_end = nInnerLoop

#ifdef MPI
  !====================================================================================
  ! Splitting the innerPhoton loop for multiple processors.
  if (myRankGlobal == 0) then
     print *, ' '
     print *, 'innerPhotonLoop computed by ', nThreadsGlobal-1, ' processors.'
     print *, ' '
  endif
  ! No need to use some processors if there are more processors
  ! than the number of photons....
  if (myRankGlobal > nPhotons/nOuterLoop - 1)  goto 666
    
  if (myRankGlobal == 0) then
     ! we will use an array to store the rank of the process
     !   which will calculate each photon
     allocate(photonBelongsRank(nPhotons/nOuterLoop))
    
     call mpiBlockHandout(nThreadsGlobal,photonBelongsRank,blockDivFactor=40,tag=tag,&
                          setDebug=.false.)
     deallocate(photonBelongsRank) ! we don't really need this here. 
  end if
  !====================================================================================

    
    
  if (myRankGlobal /= 0) then
    mpiBlockLoop: do  
      call mpiGetBlock(myRankGlobal,iInner_beg, iInner_end,rankComplete,tag,setDebug=.false.)  
      if (rankComplete) exit mpiBlockLoop  
    
#endif

!$OMP DO SCHEDULE(DYNAMIC)

        nFromEnv = 0
        innerPhotonLoop: do i = iInner_beg, iInner_end

#ifdef MPI
 !  if (MOD(i,nThreadsGlobal) /= myRankGlobal) cycle innerPhotonLoop
#endif
           ! The following six arrays must be allocated and deallocated for each 
           ! innerPhotonLoop to make the program work with OpenMP! (RK)
           allocate(lambda(1:maxTau))
           allocate(tauExt(1:maxTau))
           allocate(tauAbs(1:maxTau))
           allocate(tauSca(1:maxTau))
           allocate(linePhotonalbedo(1:maxTau))
           allocate(contTau(1:maxTau,1:nLambda)) 
           allocate(contWeightArray(1:nLambda))
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
              call random_number(r)
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
                 call initPhoton(thisPhoton, grid, nLambda, xArray, sourceSpectrum, &
                      lineEmission, lamLine, weightLinePhoton, &
                      weightContPhoton, contPhoton, flatspec, vRot, pencilBeam, &
                      secondSource, secondSourcePosition, &
                      ramanSourceVelocity, vo6, contWindPhoton, directionalWeight, useBias, &
                      theta1, theta2, chanceHotRing,  &
                      nSpot, chanceSpot, thetaSpot, phiSpot, fSpot, spotPhoton,  probDust, weightDust, weightPhoto,&
                      narrowBandImage, vMin, vMax, source, nSource, rHatinStar, energyPerPhoton, filters, mie,&
                      curtains, starSurface, forcedWavelength, usePhotonWavelength, iLambdaPhoton,VoigtProf, &
                      outVec, photonfromEnvelope, dopShift=dopShift)
                 if (thisPhoton%resonanceLine) then
                    r1 = real(i)/real(nPhotons/nOuterLoop)
                    thisPhoton%lambda = xArray(1) + r1*(xArray(nLambda)-xArray(1))
                 endif

           end select
           if (photonFromEnvelope) then
              nFromEnv = nFromEnv + 1
           endif


!           write(*,*) "after init",thisPhoton%stokes%i,contwindphoton
           observedLambda = thisPhoton%lambda
!           if (thisPhoton%contPhoton) then
!
!              meanr_cont = meanr_cont + modulus(thisPhoton%position)*thisPhoton%stokes%i
!              wtot_cont = wtot_cont + thisPhoton%stokes%i
!           else
!              meanr_line = meanr_line + modulus(thisPhoton%position)*thisPhoton%stokes%i
!              wtot_line = wtot_line + thisPhoton%stokes%i
!           endif



           iLambda = findIlambda(thisPhoton%lambda, xArray, nLambda, ok)

           ! now we fire the photon direct to the observer
           lineResAbs = .false.  ! F for the zero-th scattering.
           if (doRaman) then
              
              call integratePath(gridUsesAMR, VoigtProf, &
                         thisPhoton%lambda, lamLine, &
                         s2o(thisPhoton%velocity), &
                         thisPhoton%position, s2o(outVec), grid, &
                         lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau , nTau, thin_disc_on, opaqueCore, &
                         escProb, thisPhoton%contPhoton, lamStart, lamEnd, &
                         nLambda, contTau, hitCore, thinLine, lineResAbs, .false., &
                         .false., nUpper, nLower, 0., 0., 0., junk,&
                         sampleFreq,intPathError, &
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

              thisLam = thisPhoton%lambda + (thisPhoton%velocity .dot. viewVec) * 1031.928
              j = findIlambda(thisLam, o6xArray, no6pts, ok)
              if (ok) then
                 obs_weight = oneOnFourPi * exp(-tauExt(nTau))
                 o6yArray(j) = o6yArray(j) + obs_weight
                 thisVel = (thisLam-lamLine)/lamLine

                 call addPhotonToImage(viewVec, rotationAxis,o6Image(1), 1, &
                      thisPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre)
              endif
           endif

           if (.not. screened) then


              ! find optical depths to observer
              call integratePath(gridUsesAMR, VoigtProf, &
                         thisPhoton%lambda, lamLine, &
                         s2o(thisPhoton%velocity), &
                         thisPhoton%position, s2o(outVec), grid, &
                         lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, &
                         escProb, thisPhoton%contPhoton, lamStart, lamEnd, &
                         nLambda, contTau, hitCore, thinLine, lineResAbs, .false.,&
                         .false., nUpper, nLower, 0., 0., 0., junk,&
                         sampleFreq,intPathError, &
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
                 observedlambda = thisPhoton%lambda / fac
                 tau_tmp = tauExt(nTau)
                 exp_minus_tau = EXP(-tau_tmp)
                 obs_weight = oneOnFourPi * exp_minus_tau * escProb
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
                       obs_weight = abs(t)*exp(-tauExt(nTau)) / pi
                    else
                       obs_weight = oneOnFourPi*abs(t)*exp(-tauExt(nTau))
                    endif
!                 end if
                 observedlambda = thisPhoton%lambda
              endif

              if (.not.flatSpec.and.(.not.hitCore)) then
                 if (contWindPhoton) then
                    obs_weight = oneOnFourPi*exp(-tauExt(nTau))
                 else 
                    obs_weight = (outVec.dot.rHatinStar)*exp(-tauExt(nTau))/pi
                    if (obs_weight < 0.) obs_weight = 0.
                 endif
                 
                 iLambda = findIlambda(observedlambda, xArray,  nLambda, ok)
                 if (ok) then

!                    write(*,*) "addingPhoton",thisPhoton%stokes%i,obs_weight
                    yArray(iLambda) = yArray(iLambda) + &
                         (thisPhoton%stokes * obs_weight)
                    statArray(iLambda) = statArray(iLambda) + 1.
                 endif
                 if (stokesImage) then
                    thisVel = 0. ! no velocity for dust continuum
                    call addPhotonToImage(viewVec, rotationAxis, obsImageSet, nImage, thisPhoton,&
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
                    iLambda = findIlambda(observedlambda, xArray, nLambda, ok)
                    if (ok) then
                       yArray(iLambda) = yArray(iLambda) + &
                            (thisPhoton%stokes * (fac2 * oneOnFourPi * escProb * exp(-tauExt(nTau))))
                       statArray(iLambda) = statArray(iLambda) + 1.
                    endif

                    obs_weight = oneOnFourPi * escProb * exp(-tauExt(nTau))*fac2

                    meanr0_line = meanr0_line + modulus(thisPhoton%position) * &
                         (thisPhoton%stokes%i * obs_weight)
                    wtot0_line = wtot0_line + (thisPhoton%stokes%i * obs_weight)
                    

                    thisVel = (observedLambda-lamLine)/lamLine
                    if (stokesImage) then
                       call addPhotonToImage(viewVec, rotationAxis, obsImageSet, nImage, &
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
                       fac1 = 2.*abs(thisPhoton%originalNormal.dot.outVec)/twoPi
                    else
                       fac1 = oneOnfourPi
                    endif
                    
                    vray = -(thisPhoton%velocity .dot. outVec)
                    vovercsqr = modulus(thisPhoton%velocity)**2
                    fac = (1.d0-0.5d0*vovercsqr*(1.d0-0.5d0*vovercsqr))/(1.d0+vray)
                    observedlambda = thisPhoton%lambda / fac
                    
                    i1 = 0
!
                    do iLambda = 1, nLambda

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
                       obs_weight = (fac1 * exp(-(tauExt(ntau)+fac3)))*fac2
                       yArray(iLambda) = yArray(iLambda) + &
                                           (thisPhoton%stokes * obs_weight)
                       statArray(iLambda) = statArray(iLambda) + 1.

                       meanr0_cont = meanr0_cont + modulus(thisPhoton%position) &
                            * thisPhoton%stokes%i*obs_weight
                       wtot0_cont = wtot0_cont + thisPhoton%stokes%i*obs_weight
                       thisVel = (observedLambda-lamLine)/lamLine
                       if (stokesImage) then
                          call addPhotonToImage(viewVec,  rotationAxis,obsImageSet, nImage, &
                               thisPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre, &
                               xArray(iLambda))
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

           call integratePath(gridUsesAMR, VoigtProf, &
              thisPhoton%lambda, lamLine, &
              s2o(thisPhoton%velocity), &
              thisPhoton%position, &
              thisPhoton%direction, grid, &
              lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, &
              escProb, thisPhoton%contPhoton, lamStart, lamEnd, &
              nLambda, contTau, hitCore, thinLine, lineResAbs, .false.,  &
              .false., nUpper, nLower, 0., 0., 0., &
              junk,sampleFreq,intPathError, &
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

            if (forceFirstScat) then ! this option is off by default.
               ! Useful for optically thin case (to get correct scatterng flux).
               tau_bnd = tauExt(nTau) ! Optical depth to outer boundary.
               fac=1.0d0-exp(-tau_bnd)
               call random_number(r1)
               thisTau = min(-log(1.d0-r1*fac),tau_bnd)  
               escaped = .false.!
               ! This is done so to force the first scattering if 
               ! MAXSCAT>0. 
            else
               call random_number(r1)
               fac = 1.
               thistau = -log(max(1.e-10,(1. - r1)))
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

           currentScat = 0


!           escaped = .false.
           absorbed = .false.
           if (doRaman) hitCore = .false.


           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
           ! scattering loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
           lineResAbs = .true.
           do while (.not.escaped .and. .not.absorbed .and. &
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
                 dlambda = lambda(j) + (lambda(j+1)-lambda(j))*t
              else
                 dlambda = lambda(nTau)
              endif
                    
              ! New photon position 
              thisPhoton%position = thisPhoton%position + real(dlambda,kind=oct)*thisPhoton%direction
              ! adjusting the photon weights 
              if (.not. thisPhoton%linePhoton) then
                 if (j < nTau) then
                    contWeightArray(1:nLambda) = contWeightArray(1:nLambda) *  &
                         EXP(-(contTau(j,1:nLambda) + t*(contTau(j+1,1:nLambda)-contTau(j,1:nLambda))) )
                 else
                    contWeightArray(1:nLambda) = contWeightArray(1:nLambda)*EXP(-(contTau(nTau,1:nLambda)))
                 end if
              end if


              if (flatspec) then
                 iLambda = 1
              else
                 iLambda = findIlambda(thisPhoton%lambda, xArray, nLambda, ok)
              endif

              if (grid%adaptive) then
                 positionOc = thisPhoton%position
                 call amrGridValues(grid%octreeRoot, positionOc, grid=grid, iLambda=iLambda, &
                      kappaAbs = thisChi, kappaSca = thisSca)
              else
                 call getIndices(grid, thisPhoton%position, i1, i2, i3, t1, t2, t3)
                 if (.not.grid%oneKappa) then
                    if (.not.flatspec) then
                       thisChi = interpGridKappaAbs(grid, i1, i2, i3, iLambda, real(t1), real(t2), real(t3))
                       thisSca = interpGridKappaSca(grid, i1, i2, i3, iLambda, real(t1), real(t2), real(t3))
                    else
                       thisChi = interpGridKappaAbs(grid, i1, i2, i3, 1,  real(t1), real(t2), real(t3))
                       thisSca = interpGridKappaSca(grid, i1, i2, i3, 1,  real(t1), real(t2), real(t3))
                    endif
                 else
                    r = interpGridScalar2(grid%rho,grid%na1,grid%na2,grid%na3,i1,i2,i3,real(t1),real(t2),real(t3))
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

              if (grid%resonanceLine) albedo = 1.

              if (hitCore) then
                 if (grid%geometry /= "binary") then
                    r = modulus(thisPhoton%position)/grid%rCore
                    if (r < 1.01) then
                       absorbed = .true.
                       albedo = 0.
                    endif
                 else
                    r1 = modulus(o2s(thisPhoton%position) - grid%starpos1)/grid%rStar1
                    r2 = modulus(o2s(thisPhoton%position) - grid%starpos2)/grid%rStar2
                    r = min (r1, r2)
                    if (r < 1.01) then
                       albedo = 0.5
                    endif
                 endif
              endif
                 

              thisPhoton%stokes = thisPhoton%stokes * albedo  ! weight adjusted here!!!

              if (thisPhoton%stokes%i < reallySmall) then
                 absorbed = .true.
!                 write(*,*) "! Small photon weight",thisPhoton%stokes%i,thisPhoton%lambda,albedo,hitcore
              endif


              ! towards observer

              if (.not.absorbed) then

                 tempPhoton = thisPhoton
                 
                 call scatterPhoton(grid, tempPhoton, outVec, obsPhoton, mie, &
                       miePhase, nDustType, nLambda, nMuMie, ttau_disc_on, ttauri_disc)

                 ! the o6 photon might get scattered towards the observer by a rayleigh scattering

                 if (doRaman) then
                    call integratePath(gridUsesAMR, VoigtProf, &
                            obsPhoton%lambda, lamLine, s2o(obsPhoton%velocity), &
                            obsPhoton%position, obsPhoton%direction, grid, &
                            lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau,  thin_disc_on, opaqueCore, &
                            escProb, obsPhoton%contPhoton, lamStart, lamEnd, nLambda, contTau, hitCore, &
                            thinLine, lineResAbs, .false., .false., nUpper, nLower, 0., 0., 0.,&
                            junk, sampleFreq,intPathError, &
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
                    observedlambda = obsPhoton%lambda / fac
                    
                    obs_weight = oneOnFourPi*exp(-tauExt(nTau)) * 34./(6.6+34.)

                    j = findIlambda(observedLambda, o6xArray, no6pts, ok)
                    if (ok) then
                       o6yArray(j) = o6yArray(j) + obs_weight
                    endif

                 endif  ! raman




                 redRegion = .false.

                 if (doRaman) redRegion = .true.

                 call integratePath(gridUsesAMR, VoigtProf, &
                      obsPhoton%lambda, lamLine, &
                      s2o(obsPhoton%velocity), &
                      obsPhoton%position, obsPhoton%direction, grid, &
                      lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau,  thin_disc_on, opaqueCore, &
                      escProb, obsPhoton%contPhoton, lamStart, lamEnd, &
                      nLambda, contTau, hitCore, &
                      thinLine, lineResAbs, redRegion, &
                      .false., nUpper, nLower, 0., 0.,0.,junk,sampleFreq,intPathError, &
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
                    obs_weight = oneOnFourPi*exp(-tauExt(nTau)) * ramanWeight


                    iLambda = findIlambda(observedlambda, xArray, nLambda, ok)
                      

                    if (ok) then
                       yArray(iLambda) = yArray(iLambda) + obsPhoton%stokes*obs_weight

!                       write(*,*) obsphoton%lambda,obsPhoton%stokes%i, obs_weight, nscat
                       statArray(iLambda) = statArray(iLambda) + 1.
                    endif
                    if (doRaman) then
                       thisVel = (1./(1./observedLambda  + 1./1215.67))/1031.928-1.
                    else
                       thisVel = observedLambda
                    endif



                    if (stokesImage) then
                       thisVel = 0. ! no velocity for dust continuum emission
                       call addPhotonToImage(viewVec,  rotationAxis, obsImageSet, nImage,  &
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
                       iLambda = findIlambda(observedlambda, xArray, nLambda, ok)
                       obs_weight = oneOnFourPi*exp(-tauExt(nTau))


                       if (obsPhoton%resonanceLine) obs_weight = obs_weight * directionalWeight

                       if (ok) then
                          yArray(iLambda) = yArray(iLambda) + obsPhoton%stokes*obs_weight
                          statArray(iLambda) = statArray(iLambda) + 1.
                       endif

                       thisVel = observedLambda
                       thisVel = (observedLambda-lamLine)/lamLine
                       if (stokesImage) then
                          call addPhotonToImage(viewVec,  rotationAxis, obsImageSet, nImage, &
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
                       observedlambda = obsPhoton%lambda / fac

                       do iLambda = 1,nLambda

                          thisLam = (lamLine-observedlambda) + xArray(iLambda)

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


                          obs_weight = exp(-(tauExt(nTau)+fac3)) * oneOnFourpi & 
                          * fac2 * directionalWeight * contWeightArray(i1)

                          
                          if (ok) then
                             yArray(iLambda) = yArray(iLambda) + &
                                  (obsPhoton%stokes*obs_weight)
                          else 
                             yarray(ilambda) = yArray(ilambda) + &
                                  obsPhoton%stokes*(oneOnFourpi* directionalweight * exp(-tauExt(ntau)))
                          endif


                          thisVel = (observedlambda-lamLine)/lamLine
                          if (stokesImage) then
                             call addPhotonToImage(viewVec,  rotationAxis, obsImageSet, nImage, &
                                  obsPhoton, thisVel, obs_weight, filters, grid%octreeRoot%centre, &
                                  xArray(iLambda))
                          endif
                          if (dopvImage) then
                             do iSlit = 1, nSlit
                                call addPhotontoPVimage(pvImage(iSlit), obsPhoton, viewVec,  rotationAxis, thisVel, &
                                  obs_weight, gridDistance)
                             enddo
                          endif


                          statArray(iLambda) = statArray(iLambda) + 1.
                       enddo
                    endif
                 endif

                 call scatterPhoton(grid,thisPhoton, zeroVec, outPhoton, mie, &
                       miePhase, nDustType, nLambda, nMuMie, ttau_disc_on, ttauri_disc)
                 thisPhoton = outPhoton

                 call integratePath(gridUsesAMR, VoigtProf, &
                            thisPhoton%lambda, lamLine, &
                            s2o(thisPhoton%velocity), thisPhoton%position, &
                            thisPhoton%direction, grid, lambda, tauExt, tauAbs, &
                            tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, &
                            escProb, thisPhoton%contPhoton, lamStart, lamEnd, &
                            nLambda, contTau, hitCore, thinLine, lineResAbs, .false., &
                            .false., nUpper, nLower, 0.,&
                            0., 0., junk,sampleFreq,intPathError, &
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


                  call random_number(r1)
                  thistau = -log(max(1.e-10,(1. - r1)))
                  if (thistau .gt. tauExt(nTau)) then
                     escaped = .true.  
                  else
                     escaped = .false.  
                  endif



!                  if (maxScat == 1) then
!                     call random_number(r1)
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
!                     call random_number(r1)
!                     thisTau = -log(1.d0-r1*fac)
!                     ! (RK) commented out the next line.
!                     ! if (thisTau > tauExt(nTau)) print *, 'Fixed thisTau in scattering loop.'
!                     thisTau = min(dble(thisTau),tau_bnd)
!                     thisPhoton%stokes = thisPhoton%stokes * real(fac)
!                  endif


               endif  ! (.not. absorbed)


            enddo
            
            if ((nscat >= maxScat).and.(maxScat .ne. 0)) then
               write(*,*) "Photon aborted after ",maxScat, " scatterings"
            endif
           nTot = nTot + nScat

999  continue  ! escape route for a bad photon

           deallocate(lambda)
           deallocate(tauSca)
           deallocate(tauExt)
           deallocate(tauAbs)
           deallocate(linePhotonalbedo)
           deallocate(contTau)
           deallocate(contWeightArray)
           
!           call plotspec(xArray, yArray, nLambda)

        enddo innerPhotonLoop

!        write(*,*) "Fraction of photons from envelope: ", real(nFromEnv)/real(nInnerLoop)

!$OMP END DO
!$OMP END PARALLEL


#ifdef MPI
 if (.not.blockHandout) exit mpiblockloop        
    end do mpiBlockLoop  
  end if ! (myRankGlobal /= 0)
     write (*,'(A,I3,A,I3,A,I3,A)') 'Process ',myRankGlobal, &
                      ' waiting to sync spectra... (',iOuterLoop,'/',nOuterLoop,')' 
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
#endif
     call torus_mpi_barrier('finished syncing output. Waiting to continue...') 
     if(myRankIsZero) write(*,'(i8,a,f7.3)') iOuterLoop*nPhotons/nOuterLoop," photons done"

        errorArray(iOuterLoop,1:nLambda) = yArray(1:nLambda)

        call torus_mpi_barrier

        if (doTuning) call tune(6, "One Outer Photon Loop") ! Stop a stop watch        

        yArray(1:nLambda) = STOKESVECTOR(0.,0.,0.,0.)



     end do outerPhotonLoop ! outer photon loop

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
   do i = 1, nImage
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

     write(*,*) " "
     write(*,'(a)') "Model summary"
     write(*,'(a)') "-------------"
     write(*,*) " "

     call systemInfo(startTime,nPhotons)
     
     if (grid%adaptive) then
        print *, tooFewSamples, ' rays had 2 or less samples.'
        print *, BoundaryProbs, ' rays had numerical problems with octal boundaries.'
        print *, negativeOpacity, ' rays had problems with negative opacity values.'
     end if
     write(*,*) " "
     write(*,*) "Average # of scattering per photon:: ", real(nTot)/real(nPhotons)
     write(*,*) " " 

 end if 

!     if (.not.grid%cartesian.and.(grid%rCore /= 0.)) then
!        if (wtot_line /= 0.) write(*,*) "Mean radius of line formation",meanr_line/wtot_line/grid%rCore
!        if (wtot0_line /= 0.) write(*,*) "Mean radius of line zero",meanr0_line/wtot0_line/grid%rCore
!        if (wtot_cont /= 0.) write(*,*) "Mean radius of cont formation",meanr_cont/wtot_cont/grid%rCore
!        if (wtot0_cont /=0.) write(*,*) "Mean radius of cont zero",meanr0_cont/wtot0_cont/grid%rCore
!     endif

 if (myRankIsZero) then 
     if (nPhase == 1) then

        call writeSpectrum(outFile,  nLambda, xArray, yArray,  errorArray, nOuterLoop, &
     .false., useNdf, sed, objectDistance, jansky, SIsed, .false., lamLine)

        if (velocitySpace) then
           specFile = trim(outfile)//"_v"
           call writeSpectrum(specFile,  nLambda, xArray, yArray,  errorArray, nOuterLoop, &
             .true., useNdf, sed, objectDistance, jansky, SIsed, velocitySpace, lamLine)
        endif

     else
        write(tempChar,'(i3.3)') iPhase
        specFile = trim(outfile)//trim(tempChar)

        call writeSpectrum(specFile,  nLambda, xArray, yArray,  errorArray, nOuterLoop, &
      .false., useNdf, sed, objectDistance, jansky, SIsed, velocitySpace, lamLine)

        if (velocitySpace) then
           tempChar = trim(specFile)//"_v"
           call writeSpectrum(tempChar,  nLambda, xArray, yArray,  errorArray, nOuterLoop, &
             .true., useNdf, sed, objectDistance, jansky, SIsed, velocitySpace, lamLine)
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

     close(22)

     if (stokesimage) then
        do i = 1, nImage
           name_filter = get_filter_name(filters, i)
           bandwidth = 0.5*FWHM_filters(filters, i)  ! 1/2 of FWHM  [A]
           lambda_eff = lambda_eff_filters(filters, i) ! Effective wavelength of filter in [A]   
           write(specFile,'(a,a,a,i3.3)') trim(outfile),"_"//trim(name_filter),"_image",iPhase
           call writeImage(obsImageSet(i), specfile, objectDistance, imageInArcsec, lambda_eff, bandwidth)
           write(specFile,'(a,a,a,i3.3,a)') trim(outfile),"_"//trim(name_filter),"_image",iPhase,".fits"
           call writeFitsImage(obsImageSet(i), trim(specfile))
        end do
        if (doRaman) then
           write(specFile,'(a,a,i3.3)') trim(outfile),"_o6image",iPhase
           !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           ! Check the the effectieve wavelength and the 
           ! bandwith of O6 image here later and replace 5000 and 1.0d0 below!
           !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           call writeImage(o6Image(1), specfile, objectDistance, imageInArcsec, 5000d0, 1.0d0)
        endif
        if (get_filter_set_name(filters) == "pn") then
           write(specFile,'(a,a,i3.3,a)') trim(outfile),"_image",iPhase,".ppm"
           call writeFalseColourPPM(trim(specfile), obsImageSet, 3)
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
           ! doing this only for the first image in obsImageSet
           call plotSlitOnImage(obsImageSet(1), PVimage(iSlit), plotfile, gridDistance)
        enddo
     endif
  end if ! (myRankIsZero)

     if (stokesImage) then
        do i = 1, nImage
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
 if (myRankGlobal /= 0 .and. .not.noPhaseUpdate) call freeGrid(grid)
#endif
  call torus_mpi_barrier('waiting inside end of phase loop...')
  enddo phaseLoop

! Tidy up and finish the run 

666 continue

if (doTuning) call tune(6, "Torus Main") ! stop a stopwatch  

call writeInfo("TORUS exiting", FORINFO)
 
call torus_mpi_barrier
#ifdef MPI
if (.not. ll_sph) call MPI_FINALIZE(ierr)
#endif

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

!     if (photoionization) then
!        xArray(1) = lamStart
!        xArray(2) = lamEnd
!        nCurrent = 2
!        call refineLambdaArray(xArray, nCurrent, grid)
!        nt = nLambda - nCurrent
!        do i = 1, nt
!           fac = logLamStart + real(i)/real(nt+1)*(logLamEnd - logLamStart)
!           fac = 10.**fac
!           nCurrent=nCurrent + 1
!           xArray(nCurrent) = fac
!           call sort(nCurrent, xArray)
!        enddo
!!        do i = 1, nlambda
!!           write(*,*) xArray(i)
!!        enddo
!     endif

       if (mie) then
          if ((lambdaTau > Xarray(1)).and.(lambdaTau < xArray(nLambda))) then
             call locate(xArray, nLambda, lambdaTau, i)
             t1 = (lambdaTau - xArray(i))/(xArray(i+1)-xArray(i))
             if (t1 > 0.5) then
                write(message,*) "Replacing ",xArray(i+1), " wavelength step with ",lambdaTau
                call writeInfo(message, TRIVIAL)
                xArray(i+1) = lambdaTau
             else
                write(message,*) "Replacing ",xArray(i), " wavelength step with ",lambdaTau
                call writeInfo(message, TRIVIAL)
                xArray(i) = lambdaTau
             endif
          endif
       endif

    endif

    !
    ! Copying the wavelength array to the grid
    do i = 1, nLambda
       grid%lamArray(i) = xArray(i)
    enddo
    grid%nLambda = nLambda

  end subroutine set_up_lambda_array

!-----------------------------------------------------------------------------------------------------------------------

  subroutine windtest

    integer, parameter :: nrGrid = 1000
    real :: rGrid(nrGrid), drGrid(nrgrid)
    real,dimension(statEqMAxLevels) :: meanDepart ! for testing
    real :: treal

    do i = 1, 1000
       r = log10(grid%rInner) + (log10(grid%rOuter)-log10(grid%rInner))* real(i-1)/999.
       r = 10.**r
       r1 = log10(grid%rInner) + (log10(grid%rOuter)-log10(grid%rInner))* real(i)/999.
       r1 = 10.**r1
       rGrid(i) = r
       drGrid(i) = r1 - r
    enddo
    
    open(21,file="rDepart.dat",status="unknown",form="formatted")

    do i = 1, 1000
       meanDepart = 0.
       nt = 0
       outVec = VECTOR(1., 0., 0.)
       rVec = rGrid(i) * outVec 
       r = rgrid(i)
!       call integratePathAMR(10.e4,  lamLine, VECTOR(1.,1.,1.), &
!            rVec, outVec, grid, lambda, &
!            tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
!            .false., lamStart, lamEnd, nLambda, contTau, &
!            hitCore, thinLine,.false.,  &
!            .false., nUpper, nLower, 0., 0., 0., junk,sampleFreq,intPathError)

       do j = 1, 100
          ang = twoPi * real(j-1)/100.
          octVec = OCTALVECTOR(r*cos(ang), r*sin(ang),0.)
          call amrGridValues(grid%octreeRoot, octVec, temperature=tReal, &
              ilambda=1, N=levelPops, Ne=Ne)
          t1 = dble(treal)
          do level = 1, statEqMaxLevels, 1
             meanDepart(level) = meanDepart(level) + real(levelPops(level))/boltzSaha(level,Ne,dble(treal))
          end do
          nt = nt + 1
       enddo
       if (nt > 0) then
          if (writeoutput) write(21,*) r/grid%rInner,meanDepart/real(nt)
       endif
    Enddo
    close(21)

  end subroutine windtest

!-----------------------------------------------------------------------------------------------------------------------

  subroutine testamr

    real :: meant, meaneta
    integer, parameter :: nrGrid = 1000
    real :: rGrid(nrGrid), drGrid(nrgrid)
    real :: treal

    do i = 1, 1000
       meant = 0.
       meaneta =0.
       nt = 0
       r = log10(grid%rInner) + (log10(grid%rOuter)-log10(grid%rInner))* real(i-1)/999.
       r = 10.**r
       r1 = log10(grid%rInner) + (log10(grid%rOuter)-log10(grid%rInner))* real(i)/999.
       r1 = 10.**r1
       rGrid(i) = r
       drGrid(i) = r1 - r
    enddo

    open(21,file="r.dat",status="unknown",form="formatted")
    call locate(xarray, nLambda, 10.e4, ilambda)

    do i = 1, 1000
       meant = 0.
       meaneta =0.
       nt = 0
       outVec = VECTOR(1., 0., 0.)
       rVec = rGrid(i) * outVec 
       r = rgrid(i)
!       call integratePathAMR(10.e4,  lamLine, VECTOR(1.,1.,1.), &
!            rVec, outVec, grid, lambda, &
!            tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
!            .false., lamStart, lamEnd, nLambda, contTau, &
!            hitCore, thinLine,.false.,  &
!            .false., nUpper, nLower, 0., 0., 0., junk,sampleFreq,intPathError)



       do j = 1, 100
          ang = twoPi * real(j-1)/100.
          octVec = OCTALVECTOR(r*cos(ang), r*sin(ang),0.)
          call amrGridValues(grid%octreeRoot, octVec, temperature=treal, &
              ilambda=1, grid=grid ,etacont = eta, kappaAbs=kabs)
          if (treal > 1.) then
             t1 = dble(treal)
             meant = meant + t1
             meaneta = meaneta + eta
             nt = nt + 1
          endif
       enddo
       if (nt > 0) then
          if (writeoutput) write(21,*) r/grid%rInner,meant/real(nt), meaneta * r**2 * drgrid(i)/ real(nt)
       endif
    Enddo
    close(21)

  end subroutine testamr

!-----------------------------------------------------------------------------------------------------------------------

  subroutine fill_opacities_noamr( ok )

    logical, intent(out) :: ok

         ! fill up the grid with the appropriate opacities
    ok = .true.
    if (.not.plezModelOn) then
       select case(geometry)
       case("torus")
          call fillGridTorus(grid, rho, rTorus, rOuter)
       case("sphere")
          call fillGridSpheriod(grid, rho, radius, kFac)
       case("flared")
          call fillGridFlaredDisk(grid, meanDustParticleMass)
       case("ellipse")
          call fillGridEllipse(Grid,rho, rMin, rMaj, rinner, teff)
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
           call fillGridPuls(grid, mDot, rcore, tEff, v0, vterm, beta, xfac, blobs, maxBlobs, .false., vContrast)
   !        call initGridStateq(grid, contFluxFile, contFluxFile2, popFilename, &
   !                            readPops, writePops, lte, nLower, nUpper)
        case("wind")
!           call fillGridWind(grid, mDot, rStar, tEff, v0, vterm, beta, &
!           lte, contFluxFile, writePops, readPops, popFilename, nLower, nUpper)

           call initgridstateq(grid, contfluxFile, " ", popFileName, &
                readPops, writePops, lte, nLower, nUpper)
           grid%etaCont = 1.e-30
           grid%kappaAbs = 1.e-30
        
        case("resonance")
           call fillGridResonance(grid, rCore, mDot, vTerm, beta, temp1)

        case DEFAULT
           if (writeoutput) write(*,*) "! Unrecognised grid geometry: ",trim(geometry)
           ok = .false.
           return
        end select
     endif ! (.not.plezModelOn)

     if (fillTio) then
        call fillGridTio(grid, scale)
     endif

     if (fillRayleighOpacity) then
        call fillGridRayleigh(grid,scale)
     endif 

     if (fillThomson) then
        call fillGridThomson(grid)
     endif


  end subroutine fill_opacities_noamr

!-----------------------------------------------------------------------------------------------------------------------

  subroutine distort_noamr

  ! the distortion types

     select case(distortionType)
     case("spiral")
        call distortGridSpiral(grid, vRot, nSpiral)
     case("rotation")
        call distortRotation(grid, vRot)
        if (writeoutput) write(*,'(a,f5.1,a)') "Grid distorted by a rotational velocity of ", &
                               vRot/1.e5," km/s"
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

  end subroutine distort_noamr

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


        if (writePhasePops) then
          write(tempChar,'(i3.3)') nStartPhase
          phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
          call writeAMRgrid(phasePopFilename,writeFileFormatted,grid)
        else if (writePops) then
          call writeAMRgrid(popFilename,writeFileFormatted,grid)
        end if


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


        if (noScattering) then
           if (writeoutput) write(*,*) "! WARNING: Scattering opacity turned off in model"
           grid%oneKappaSca(1:nDustType,1:nLambda) = TINY(grid%oneKappaSca)
        endif
        
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

  ! Plot desired AMR grid value here... This is more generalized
  ! version of fancyAmrPlot.
  !
  ! See grid_mod.f90 for details.
  !
  ! Plotting some grid values
  call plot_AMR_values(grid, "rho", plane_for_plot, val_3rd_dim, &
       "rho_grid",.true., .true., nmarker, xmarker, ymarker, zmarker, &
       width_3rd_dim, show_value_3rd_dim, suffix="default", index=num_calls, &
       fixValMin=sph_rho_min, fixValMax=sph_rho_max, useFixedRange=ll_sph )
  call plot_AMR_values(grid, "rho", plane_for_plot, val_3rd_dim, &
       "rho_zoom",.true., .true., nmarker, xmarker, ymarker, zmarker, &
       width_3rd_dim, show_value_3rd_dim, boxfac=zoomFactor, suffix="default", index=num_calls, &
       fixValMin=sph_rho_min, fixValMax=sph_rho_max, useFixedRange=ll_sph )
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
       width_3rd_dim, show_value_3rd_dim, suffix="default", index=num_calls)
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

    case("shakara","clumpydisc","wrshell","warpeddisc")
       nSource = 1
       allocate(source(1:1))
       source(1)%radius = grid%rCore
       source(1)%teff = teff   
       source(1)%position = VECTOR(0.,0.,0.)
       tmp = source(1)%radius * 1.e10  ! [cm]
       source(1)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(1)%teff)**4
       if (contFluxfile .eq. "blackbody") then
          call fillSpectrumBB(source(1)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd),nLambda)
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
               dble(lamStart), dble(lamEnd),nLambda)
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
               dble(lamStart), dble(lamEnd),nLambda)
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

#ifdef SPH
end subroutine torus
#else
end program torus
#endif

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
