  !   a large number of variables are passed to the torus main
  !   program by the inputs subroutine. they are all collected
  !   here so that they only need to be defined once.  nhs
  !
  !   These now are inlcuded in inputs_mod and can be declared as 
  !   'protected'. This is a F2003 feature which means that only 
  !   the module containing the variable can change it - other modules
  !   have read-only access (i.e. it is like intent(in)). D. Acreman 

!----------------------------------
! Physical ingredients of the model
!----------------------------------

  logical :: molecularPhysics
  logical :: photoionPhysics
  logical :: atomicPhysics
  logical :: nbodyPhysics
  logical :: dustPhysics

!--------------------------
! Type of model calculation
!--------------------------

  logical :: radiationHydrodynamics 
  logical :: radiativeEquilibrium
  logical :: statisticalEquilibrium
  logical :: photoionEquilibrium
  logical :: xrayCalc
  logical :: pdrCalc
  logical :: hydrodynamics
  logical :: timeDependentRT
  logical :: donBodyOnly

!pdr param
  integer :: hlevel !level of healpix refinement
  real(double) :: v_turb
!-------------------------------
! Type of output to be generated
!-------------------------------

  logical :: calcDataCube
  logical :: calcPhotometry
  logical :: calcImage
  logical :: calcMovie
  Logical :: calcSpectrum
  logical :: calcBenchmark

  logical :: dowriteRadialFile
  character(len=200) :: radialFilename

!-----------------------------------
! Write/read a grid for a warm start
!-----------------------------------
 
  logical :: readGrid, writeGrid  ! Do we read/write the AMR grid from/to a file
  logical :: doSetupAMRgrid       ! sometimes (rarely) we don't want a grid at all
  logical :: multimodels          ! perform calculation on multiple models
  integer :: nModelStart, nModelEnd ! start and end numbers for multiple models
  integer :: iModel               ! current model number
  logical :: justDump             !Dump a vtk file for the read in grid and exit
  logical :: rhofromtable
  character(len=80) :: rhofile
  integer :: nrholines
  logical :: densitySpectrum      !Dump a density spectrum for the read in grid and exit
  real(double) :: normfac         ! density spectrum normalization factor
  logical :: dumpBisbas           !Dump the grid for use in the Bisbas+ 3D-PDR code
  logical :: singleMegaPhoto             !Do an extensive photoionization calculation
  character(len=80) :: gridInputFilename, gridOutputFilename ! File names for reading/writing AMR grid

!----------------------------------
! Physical units of the calculation
!----------------------------------

  real(double) :: massUnit, timeUnit, lengthUnit !Code units

!-----------------
! Photoionisation 
!-----------------

  logical :: photoionization   !Perform a photoionization calculation
  logical :: hOnly             !Hydrogen only model (no Helium)
  logical :: massiveStars      !only use stars with M>20Msol in photo calc
  logical :: simpleMu          !simpified mean particle mass calc
  logical :: caseB             !use a case b recombination coefficient
  logical :: startFromNeutral  !Start photoionization loop from neutral
  logical :: usemetals         !Include species heavier than Helium
  logical :: usexraymetals     !Include x-ray metals
  logical :: useionparam       !use ionization parameter(T) in x-ray calc
  logical :: xrayonly
  logical :: checkForPhoto     !Check whether or not a photoionization loop is necessary
  logical :: monochromatic     !Use a monochromatic radiation field
  real(double) :: inputEV     !energy of monochromatic photons
  logical :: quickThermal      !Use a simplified thermal balance calculation
  logical :: mergeBoundSinks   ! merge gravitationally bound sinks
  logical :: addSinkParticles  ! add new sink particles
  logical :: stellarWinds      ! include stellar winds in hydro calc
  logical :: dumpregularVTUS   !dump vtu after every photo step
  ! Parameters  specific to domain decomposed photoionisation 

  !Stack optimization
  logical :: optimizeStack     !Perform run-time photon stack size optimization
  integer :: dStack            !Stack increment for run-time photon stack size optimization 
  integer :: stackLimit        !Maximum stack size for run-time photon stack size optimization
  integer, allocatable :: stackLimitArray(:)   !Maximum stack size for run-time photon stack size optimization
  logical :: customStacks      !use custom stack sizes

  integer :: bufferCap         !Number of photon stacks to accomodate in buffer

  logical :: binPhotons        !Dump a spectrum of propagated photon packets
  logical :: noDiffuseField    !Use the on the spot approximation
  logical :: dustOnly          !Consider dust physics only
  integer(bigInt) :: inputNMonte       !Number of photon packets to propagate
  integer :: maxPhotoIonIter   !Maximimum interation number

  logical :: periodicX, periodicY, periodicZ ! Periodic photon boundary conditions
  logical :: zBoundaryReflecting ! zboundary reflects photons for 2d cylindrical hydro
  logical :: doPhotoRefine     !Refine the AMR grid between iterations of the photo loop

!---------------
! Hydrodynamics
!---------------

  integer, protected :: nHydroThreadsInput !Number of hydrothreads for domain decomposition
  logical :: rhieChow                 !Use Rhie-Chow interpolation
  logical :: doSelfGrav               !Do self gravity calculation
  logical :: simpleGrav               !Do self gravity calculation
  logical :: doGasGravity             ! Include gas gravity in calculation
  logical :: dirichlet                !Use dirichlet boundary conditions - otherwise periodic used by default
  logical :: severeDamping            ! Turn on damping
  logical :: includePressureTerms     ! include pressure source terms
  logical :: dumpRadial               ! write a text radial cut each dump
  logical :: radiationPressure        ! use radiation pressure terms
  real :: cflNumber                   !Courant-Friedrichs-Lewy constant
  real(double) :: rhoFloor            !min density in grid
  real(double) :: etaViscosity        !Artificial viscosity parameter
  logical :: useTensorViscosity       ! Use tensor form for artificial viscosity
  logical :: cylindricalHydro         ! perform the hydrodynamics in cylindrical coordinates
  logical :: sphericalHydro           ! perform the hydrodynamics in spherical coordinates
  real(double) :: tStart, tEnd, tDump !Start, end and dump interval times
  real(double) :: rhoThreshold        ! threshold density for sink creation
  real(double) :: hydroSpeedLimit     ! fudge to limit hydrodynamic advection speed
  logical :: hydrovelocityConv        !Convert input velocity vector into simulation velocities 
  logical :: doRefine, doUnrefine     !Allow the AMR grid to refine/unrefine
  logical :: useViscosity             !Use artificial viscosity
  logical :: fluxinterp               !Interpolate fluxes at fine to coarse cell advections
  logical :: fixedRhoBound            !Use fixed density boundary conditions
  real(double) :: rho_const           !Density of fixed density boundary conditions
  character(len=20) :: limiterType    !Flux limiter type
  integer :: idump                    !Hydrodynamics time step dump number
  integer :: vtuToGrid
  real(double) :: gridDistanceScale   !Scale of grid
  integer :: CD_version               !Which version of contact discontinuity test to run? (1,2,3 or 4, StarBench)

  logical :: readTurb !read in a turbulent velocity grid
  character(len=20) :: turbvelfilex    !file for tubrulent velocty field (x)
  character(len=20) :: turbvelfiley    !file for tubrulent velocty field (y)
  character(len=20) :: turbvelfilez    !file for tubrulent velocty field (z)
  integer :: nTurbLines                !No. of lines in turbvel files

!inflow condition parameteres
  real(double) :: inflowPressure
  real(double) :: inflowRho
  real(double) :: inflowRhoE
  real(double) :: inflowEnergy
  real(double) :: inflowMomentum
  real(double) :: inflowTemp
  real(double) :: inflowSpeed

!gradient direction of inflow conditions
  logical :: xslope
  logical :: yslope
  logical :: zslope

!gradient value of inflow conditions
  real(double) :: momgrad
  real(double) :: egrad
  real(double) :: rhoegrad
  real(double) :: pgrad
  real(double) :: rhograd


  !Boundary conditions (strings)
  character(len=20) :: xminusboundString, yminusboundString, zminusboundString
  character(len=20) :: xplusboundString, yplusboundString, zplusboundString

  !Boundary conditions (integer codes)
  integer :: xminusbound, yminusbound, zminusbound 
  integer :: xplusbound, yplusbound, zplusbound 

!---------------------------------------
! Atomic physics (co-moving frame, CMF)
!---------------------------------------

  logical :: cmf, sobolev
  logical :: opticallyThickContinuum     
  logical :: lineOff               
  logical :: statEq2d       ! whether statEq can be run in 2-D
  integer :: nAtom      
  character(len=20) :: atomFilename(10)
  real(double) :: xAbundance, yAbundance 
  integer :: iTransLine, iTransAtom ! Not set !!!
  logical :: statEq1stOctant ! T if do stateqAMR routine for the octals in first octant
                             ! then later values are mapped to other octants. Not Set !!!

!----------------------------------
! Molecular physics (molecular_mod)
!----------------------------------

  logical :: molecular
  logical, protected :: constantAbundance
  logical, protected :: removeHotMolecular
  real, protected :: molAbundance
  character(len=80), protected :: molFilename
  character(len=80), protected :: moleculefile
  logical, protected :: restart
  logical, protected :: molRestartTest
  logical, protected :: addnewmoldata
  logical, protected :: dongstep ! controls Ng Acceleration
  logical, protected :: modelwashydro
  logical, protected :: forceIniRay
  logical, protected :: zeroghosts
  logical, protected :: renewinputrays
  logical, protected :: plotlevels   
  logical, protected :: writetempfits
  logical, protected :: doCOchemistry
  logical :: isinlte ! assume grid is in LTE
  logical, protected :: setupMolecularLteOnly ! Set up LTE level populations then exit molecular loop
  logical :: lowmemory ! if memory is tight
  logical, protected :: noturb ! Subsonic turbulent velocity
  real, protected :: tolerance ! maximum acceptable fractional change for J levels between iterations in molecular_mod
  integer, protected :: initnray ! number of rays to use in fixed ray case (stage 1)
  real(double), protected :: rotateViewAboutX, rotateViewAboutY, rotateViewAboutZ
  integer, protected :: setmaxlevel ! Maximum level to be converged
  real, protected    :: setmaxleveltemp ! Temperature of maximum level to be converged
  logical :: lineimage
  logical :: maxrhocalc
  logical, protected :: quasi ! use quasirandom numbers
  logical, protected :: rgbCube ! reverse velocity axis 
  logical, protected :: densitysubsample ! do you want to subsample density in your images?
  logical, protected :: wanttau ! Output tau cube? (also used in angularImage_mod)
  real(double), protected :: centrevecx, centrevecy, centrevecz ! cube centre (also in AngularImage_mod)
  real, protected :: beamSize
  logical, protected :: getdepartcoeffs, gettau
  logical, protected :: outputconvergence, dotune
  real, protected :: x_D, x_0 ! CO drop model
  logical ::  realdust ! molecular_mod includes continuum emission
  integer, protected :: nsubpixels ! No. of sub-pixels for ray trace
  integer :: itrans
  integer, protected :: observerpos ! position of observer in molecular_mod image generation


!-----------------------------
! Drabek model parameters    |
!-----------------------------
  integer :: density_Code
  real(double) :: peakRho
  real(double) :: meanT
  real(double) :: nCol
  real(double) :: n2max
  real(double) :: tKinetic

!-----------------------------

!----------------------------  
! Lucy radiative equilibrium
!----------------------------

  logical :: lucyRadiativeEq
  logical :: solveVerticalHydro
  integer :: nHydro
  integer :: nLucy
  integer :: iterLucy
  logical :: narrowBandImage
  real    :: lucy_undersampled  
  logical :: convergeOnUndersampled
  integer :: minCrossings
  logical :: forceLucyConv
  logical :: multiLucyFiles
  logical :: polarizationImages

  ! For diffusion approximation
  real :: tauDiff, tauForce
  logical :: resetDiffusion
  real :: eDensTol

!-----------------------------
! Wavelength array parameters
!-----------------------------

  integer :: nLambda
  logical :: oneKappa
  logical :: lamfile              ! Read in wavelengths from a file? 
  character(len=80):: lamfilename ! File to use if lamfile=T

!------------------- 
! Source parameters
!-------------------

  logical :: readSources
  logical :: moveSources
  character(len=80) :: sourceFilename
  logical :: sourceHistory
  character(len=80) :: sourceHistoryFilename
  integer :: inputNsource
  real(double) :: sourceTeff(10), sourceMass(10), sourceRadius(10), sourceProb(10), sourceMdot(10)
  real(double) :: sourceLum(10)
  real(double) :: accretionRadius
  logical :: stellarSource(10)
  character(len=10) :: diffuseType(10)
  type(VECTOR) :: sourcePos(10), sourceVel(10)
  character(len=80) :: inputContFluxFile(10)
  character(len=80) :: contFluxFile

  ! variables for a second source of radiation
  logical :: secondSource                    ! second photon source?
  type(VECTOR) :: secondSourcePosition       ! the position of it
  real :: binarySep
  real :: momRatio
  real(double) :: inputEps ! input gravity softening length

  ! Variance Reduction (spectrum_mod)
  logical :: biasToLyman
  real(double) :: biasMagnitude
  logical :: screened
  logical :: postsublimate

!-----------------
! Dust parameters
!-----------------

  integer, parameter :: maxDustTypes = 10
  character(len=80) :: grainType(maxDustTypes) ! sil_ow, sil_oc, sil_dl, amc_hn, sic_pg, gr1_dl, gr2_dl  
  real :: grainFrac(maxDustTypes)
  real :: grainDensity(maxDustTypes)
  logical :: mie
  logical :: storescattered
  real :: scatteredLightWavelength
  logical :: readMiePhase       ! used in dust_mod, not set !!!
  logical :: writeMiePhase      ! used in dust_mod, not set !!!
  logical :: useOldMiePhaseCalc ! used in dust_mod, not set !!!
  logical :: includeGasOpacity
  logical :: isotropicScattering
  real :: dusttogas
  logical :: dustfile
  character(len=80) :: dustfilename(10) ! used in dust_mod, not set !!!
  logical :: variableDustSublimation
  real(double) :: tSub ! variable dust sublimation temperature factor
  logical :: dustSettling
  integer :: nDustType
  logical :: readDustFromFile, writeDustToFile
  logical :: useDust
  ! abundances of different types of dust grains. These will be used when 
  ! the graintype assigned is "mixed."
  integer :: ngrain 
  real :: X_grain(10)    ! abundaunce 
  character(LEN=30) :: grainname(10)
  !
  ! size distribution of dust grain is now assumed to have 
  ! the following form:
  !  
  !   n(a) =  const * a^-q * e^((-a/a0)^p)
  real :: aMin(maxDustTypes)    !  The maximun size in microns. 
  real :: aMax(maxDustTypes)    !  The minimum size in microns.
  real :: qDist(maxDustTypes)   !  q exponet in the equation above.
  real :: a0(maxDustTypes)      !  scale length in the equation above.
  real :: pDist(maxDustTypes)   !  p exponent in the equation above.
  real(double) :: dustHeight(maxDustTypes)   !  p exponent in the equation above.
  real(double) :: dustBeta(maxDustTypes)   !  p exponent in the equation above.

!------------------------------
! Inclinations for SEDs/images 
! Now held in sed_mod
!------------------------------

  real :: thisInclination  ! Inclination when atomicPhysics=T (calculateAtomSpectrum and compute_obs_line_flux)
  logical :: freefreeSed   !include free-free emission in SED
  logical :: recombinationSed   !include recombination line emission in SED
  logical :: forbiddenSed   !include forbidden line emission in SED
  logical :: resolveSilicateFeature ! add points to SED to full resolve silicate feature
  logical :: dumpCut    !Dump pixel values in a cut across the image into a file
  character(len=30) :: cutType    !direction of cut (horizontal or vertical)
  integer :: sliceIndex !pixel through which the cut runs (x or y depends on cutType)

!--------------------------------
! Image and data cube parameters 
!--------------------------------

  logical :: monteCarloRT
  real(double) :: minVel, maxVel ! for molecular_mod and angularImage_mod
  real :: imageside    ! data cube size for molecular_mod and angularImage_mod
  integer :: ncubes ! number of data cubes
  integer :: nv ! number of velocity channels
  real :: vMinSpec, vMaxSpec ! For atomicDataCube option
  real :: gridDistance ! distance of observer for images, SEDs etc. 
  character(LEN=30) :: filter_set_name  ! name of filter set used for images (phaseloop_mod)

!----------------------------------------------------------
! Geometry specific parameters (this is a large section !)
!-----------------------------------------------------------
  real(double) :: slabwidth ! width of slab in bubble calcs
  real(double) :: lXOverLbol ! ??? Sources for Ttauri geometry
  real(double) :: maxCellMass ! ??? Refinement control for Ttauri geometry
  real(double) :: mStarburst, clusterRadius ! ??? for starburst geometry
  logical :: hosokawaTracks ! use Hosokawa tracks for protostellar properties
  logical :: usePacketSplitting ! use photon packet splitting method
  integer :: inputNSmallPackets
  character(len=80) :: intextFilename, outtextFilename ! ??? Fogel geometry
  logical :: doSpiral ! For a Shakara Sunyaev Disc
  type(VECTOR) :: sphereVelocity, spherePosition ! unisphere geometry
  real(double) :: sphereMass, sphereRadius       ! unisphere geometry
  real(double) :: omega
  integer :: nSphereSurface ! number of points on spherical surface
  real :: rpower ! radial density power  r^-rpower. For proto geometry

  ! whitney stuff
  real :: erInner, erOuter, drInner, drOuter, cavAngle, cavDens
  real(double) :: rCavity
  real :: rStellar, mdotEnv, mEnv

  ! variables for clumped wind models
  integer :: nBlobs
  logical :: freshBlobs
  real :: blobContrast

  ! Filenames for intrinsic core absorption profile for core-halo models
  character(len=80) :: intProFilename
  character(len=80) :: intProFilename2

  ! model parameters
  !
  real :: height
  real :: heightSplitFac
  real :: sigma0
  real :: rMin, rMaj
  real :: shellFrac
  real :: Teff, Teff1, Teff2
  real :: twind
  real :: rstar1, rstar2, mstar1, mstar2
  real :: mdot1, mdot2
  real :: Mbol
  real :: radius, kfac, xfac
  real :: rCore, rInner
  real :: rTorus, rOuter, rSublimation
  real :: rho, rho0
  real(double) :: extMass
  real :: scale, rscale
  real :: mCore, diskTemp, mDisc
  real :: epsilonDisc
  real :: alphaDisc, betaDisc, alphaDiscTemp
  real :: scaleDensity
  real :: vRot
  real :: mdot
  real :: beta
  real :: vterm
  real :: v0
  logical :: useHartmannTemp ! use T Tauri accretion stream temperatures
                             !   from Hartmann paper
  real :: maxHartTemp        ! maximum temperature of hartmann distribution                             
  logical :: isoTherm        ! use isothermal T Tauri accretion stream 

!limb darkening:
  real(double) :: sourcelimbaB, sourcelimbbB, sourcelimbaV, sourcelimbbV

  ! pp disk stuff
  real :: rSmooth, rHeight, sigmaPower
  real :: flaringPower, gapViscAlpha
  real :: rGap, mPlanet, gapWidth
  real :: rGapInner, rGapOuter, rhoGap, deltaCav
  logical :: planetGap

  ! single dust blob parameters (WR137 type model)
  !
  real :: dustBlobDistance, phiDustBlob
  real :: xDustBlobSize, yDustBlobSize, zDustBlobSize

  ! Spot stuff
  !
  integer :: nSpot                       ! number of spots
  real :: fSpot                          ! factional area coverage of spots
  real :: tSpot                          ! spot temperatures
  real :: thetaSpot, phiSpot             ! spot coords

  ! binary parameters
  !
  real :: period
  real :: shockWidth, shockFac
  real :: vNought1, vNought2
  real :: beta1, beta2
  real :: vTerm1 , vTerm2
  real :: radius1, radius2
  logical :: readPops, writePops
  logical :: readPhasePops, writePhasePops
  real :: tThresh
  logical :: readLucy
  real :: mass1, mass2, massRatio, streamFac
  real :: deflectionAngle
  logical :: lte
  logical :: LyContThick
  logical :: curtains, enhance
  real :: dipoleOffset

  ! T Tauri parameters ----------------------------------------------------
  character(len=80) :: MdotType ! variable accretion rate model in use
  real :: MdotParameter1, MdotParameter2, MdotParameter3, MdotParameter4
  real :: MdotParameter5, MdotParameter6
  real :: TTauriRinner, TTauriRouter ! disc sizes (in R_star units)
  real :: Thotspot
  logical :: ttauriwarp, ttauriwind, ttauridisc, ttauriMagnetosphere
  real :: TTauriRstar ! stellar radius (in R_sol units)
  real(double) :: holeRadius ! radius of inner hole to geometrically thin, optically thick disc
  real :: TTauriMstar ! stellar mass   (in M_sol units)
  real :: TTauriDiskHeight ! (in R_star units)
  real :: TTauriDiskRin    ! (in R_star units)
  real(double) :: hOverR
  real :: ThinDiskRin      ! (in R_star units)
  real :: curtainsPhi1s ! accretion curtains from (s)tart... 
  real(double) :: phiRefine, dphiRefine, minPhiResolution
  real :: curtainsPhi1e ! ... to (e)nd angle
  real :: curtainsPhi2s ! (all in degrees)
  real :: curtainsPhi2e ! must be: phi1s<phi1e<phi2s<phi2e
  ! The following two are used for constantcurtain geometry (RK)
  integer :: curtain_number ! nuumber of curtains
  real    :: curtain_width  ! Width of curtain in degrees.

  ! suboption for ttauri geometry
  logical :: ttau_acc_on        ! T to include magnetosphere
  logical :: ttau_disc_on       ! T to include disc
  logical :: ttau_discwind_on   ! T to include disc wind.
  logical :: ttau_jet_on        ! T to include jets.
  logical :: ttau_fuzzy_edge    ! T to use fuzzy edge for accretion flow.

  logical :: formalsol          ! T to perform formal solution
  integer :: form_nphi          ! # of angular poins for formal integration
  integer :: form_nr_core       ! # of radial poins for formal integration (core)
  integer :: form_nr_acc        ! # of radial poins for formal integration (accretion)
  integer :: form_nr_wind       ! # of radial poins for formal integration (wind)
  logical :: do_pos_disp        ! if T perform position displace ment calculation

  !--------------------------------------------------------------------

  ! Use this parameter to turn off the alpha disc, jets and magnetosphere accretion
  logical :: ttau_turn_off_disc  
  logical :: ttau_turn_off_jet
  logical :: ttau_turn_off_acc  ! magenetoshere

  logical :: pointSource, pointSourceArray(50)
  real(double) :: biasPhiDirection, biasPhiInterval, biasPhiProb

  !------ The disc wind parameters follows here -----------------------
  real(double) :: DW_d           ![10^10cm] displacement of souce point from the center of star
  real(double) :: DW_Rmin        ! the inner most radius of the disc [10^10cm]
  real(double) :: DW_Rmax        ! outer limit of the disc [10^10cm]
  real(double) :: DW_theta
  !
  ! Temperature : T(R)= Tmax*(R/Rmin)^gamma where R is the distance from the center
  !               along the disc along the disc
  real(double) :: DW_Tmax        ! [K] Temperature at the inner edge of the disc
  real(double) :: DW_Temperature
  real(double) :: DW_gamma       ! exponet in the temperature power low: 
  !
  ! mass loss rate per unit area (from the disc)
  !                                       Mdot*R^delta
  !     mdot_per_area =   ---------------------------------------------
  !                         4Pi*(Rmax ^(delta+2) - Rmin^(delta+2))
  ! where delta = 4*alpha*gamma
  !
  real(double) :: DW_Mdot        ! [Msun/yr] total mass-loss rate 
  real(double) :: DW_alpha       ! [-] exponent in the mass-loss rate per unit area
  !
  ! modefied  beta-velocity low
  !                                                  Rs
  !  V(r) = Cs(R) + ( f*Vesc(R) - Cs(R) ) * ( 1 - ------- )^beta
  !                                                s - Rs
  ! 
  !  Cs -- speed of sound
  !  f  -- scaling of the asymptotic terminal velocity
  !  Vesc(R) -- escape velocity from R.
  !  s -- distance from the disc along a stream line. Note this l in the paper.
  !  Rs -- constant effective accerelation length
  real(double) :: DW_beta  ! [-]
  real(double) :: DW_Rs    ! [10^10 cm]  usually 50 times of Rmin
  real(double) :: DW_f     ! [-]  usually 2.0
  !
  ! temperature of the disc wind 
  !  -- set to be isothermal for now.
  real(double) :: DW_Twind     ! [Kelvin] Isothermal temperature of the wind
  !
  real(double) :: DW_Mstar ! [M_sun]  mass of the central object
  !-----------------------------------------------------------------------------


  !----- For T Tauri Jets -----------------------------------------------------
  real(double) :: JET_Rmin    !  [10^10cm]  The minimum raidus of jet
  real(double) :: JET_theta_j !  [radian]  jet opening angle
  !
  real(double) :: JET_Mdot    ! [Msun/yr] mass loss rate in the jets
  real(double) :: JET_a_param ! [-] a parameter in density function 
  real(double) :: JET_b_param ! [-] a parameter in density function
  real(double) :: JET_Vbase   ! [km/s] Base velocity of jets
  real(double) :: JET_Vinf    ! [km/s] Terminal velocity of jets
  real(double) :: JET_beta    ! [-] a parameter in velocity function
  real(double) :: JET_gamma   ! [-] a parameter in velocity function
  !
  real(double) :: JET_T       ! [K]  Isothermal temperature of jets
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------


  ! For luc_cir3d geometry ------------------------------------------------------
  real(double)    :: CIR_Rstar       ! radius of central star  [R_sun]
  real(double)    :: CIR_Mass        ! [M_sun]  Mass of the star
  real(double)    :: CIR_Twind       ! [K]  Isothemal temperature of the stellar wind
  real(double)    :: CIR_Mdot_scale  ! [-]  Scaling factor for CIR density and Mdot
  !---------------------------------------------------------------------------------


  !---------------------------------------------------------------------------------
  logical         :: uniformStar
  real            :: bigomega, eddingtonGamma, alphaCAK
  real(double)    :: CMFGEN_Rmin       ! radius of central star  [10^10cm]
  real(double)    :: CMFGEN_Rmax        ! max radius
  !---------------------------------------------------------------------------------
  

  ! For "romanova" geometry --------------------------------------------------------
  real(double)      :: ROM_Rs          ! radius of central star  [Rsun]. 
  real(double)      :: ROM_Mass        ! [Msun]  Mass of the star
  logical           :: ROM_isoT        ! if T isothermal othewise use data
  real(double)      :: ROM_T_flow      ! [K]  Isothemal temperature of the flow
  ! Reference values for dimensionless units used in Romanova's model
  real(double)      :: ROM_r_ref       ! [cm] Reference length value
  real(double)      :: ROM_rho_ref     ! [g/cm^3]  Reference density value
  real(double)      :: ROM_T_ref       ! [T]       Reference temperature value
  real(double)      :: ROM_v_ref       ! [cm/s]    Reference speed 
  ! The tile angle of the magnetic axis 
  real(double)      :: ROM_tilt        ! [degrees]  will be changed to [rad].
  real(double)      :: ROM_Period      ! [day]  will be changed to [sec].
  !
  character(LEN=60) :: ROM_datafile    ! Name of the romanova's data file
  !---------------------------------------------------------------------------------

  !  Input parameters for bipolar jets geometry.
  !
  !  Bipolar flow. Flow into cones with openning angle theta_o around
  !  z-axis as the symetry axis.
  !
  !  The velocity is assumed to follow the classical
  !  beta velocity law with an offset:
  !
  !  In this formulation:
  !
  !                    Mdot
  !  rho(r) = ------------------------------
  !            4Pi V(r) r^2 (1-Cos(theta_o))
  !
  !  and
  !                         r*
  !   V(r) =  Vinf * ( 1 - --- )^beta + Vo
  !                         r
  !
  !   where Vinf is the terminal speed and
  !         Vo is the offeset speed.
  ! theta_o  is "half" opening angle.
  !
  ! For temperature structure of the bipolar is assumed as:
  !                Rmin
  !    T = Tcore*(------)^e6
  !                  R
  !  
  !===============================================================  
  real :: Rmin_bp    ! radius of central star  [10^10cm]
  real :: Rmax_bp    ! cutoff radius in [10^10 cm]
  real :: Vo_bp      ! a small offset in V in  [km/s]
     
  ! For jets
  real :: Vinf_jets    ! terminal velocity in  [kms]
  real :: beta_jets    ! beta in the beta-belocity law [-]
  real :: Mdot_jets    ! mass loss rate in jets [M_solar/yr]
  real :: theta_o_jets ! half opening angle [degrees]     
  real :: Tcore_jets   ! temperature at the core at Rmin in [10^4 K]
  real :: e6_jets      ! an exponet of in tempereture eq. [-]
  
  ! For disk wind
  real :: Vinf_disk    ! terminal velocity in  [kms]
  real :: beta_disk    ! beta in the beta-belocity law [-]
  real :: Mdot_disk    ! mass loss rate in jets [M_solar/yr]
  real :: theta_o_disk ! half opening angle [degrees]     
  real :: Tcore_disk   ! temperature at the core at Rmin in [10^4 K]
  real :: e6_disk      ! an exponet of in tempereture eq. [-]
  
  ! For equatorial disk
  real :: Rdisk_min  ! The minimum radius of the disk. in 10^10 cm.
  real :: Rdisk_max  ! The minimum  of the disk in 10^10 cm.
  real :: h_disk     ! Thickness of the disk in 10^10 cm.
  real :: rho_scale  ! The density in the units of rho max for disk wind.     

  real :: tauExtra   ! foreground optical depth
  real :: tauExtra2  ! foreground optical depth

  real :: warpRadius, warpSigma, warpFracHeight, warpAngle
  logical :: hydroWarp

  ! "magstream" parameters ---------------------------------------------
  logical :: limitSpotTemp ! whether to limit hot spot temperatures
  real :: maxSpotTemp ! maximum hot spot temperature
  logical :: scaleFlowRho ! whether to rescale accretion flow densities
  real :: flowRhoScale ! rescaling factor for accretion flow densities
  character(len=80) :: magStreamFile ! filename for accretion stream data
  logical :: isothermStream ! is accretion stream isotherma;
  real :: isothermTemp ! temperature for accretion stream
  logical :: magStreamFileDegrees ! does input file use degrees? (otherwise radians)
  !--------------------------------------------------------------------

  ! Dimensionless cutoff radius for Bonnor-Ebert Sphere
  real(double) :: zetacutoff 

  ! empty geometry
  real(double) :: centralMass 


!----------------------
! adaptive mesh stuff 
!----------------------

  logical :: gridUsesAMR    ! true if grid is adaptive
  logical :: splitOverMPI   ! true if grid is domain decomposed 
  logical :: amr1d, amr2d, amr3d
  logical :: gridshuffle    !repeatedly reset the starting grid to achieve better refinement
  integer :: minDepthAMR, maxDepthAMR
  real(double) :: limitScalar  ! value for controlling grid subdivision 
  real(double) :: limitScalar2 ! value for controlling grid subdivision 
  real(double) :: vturbmultiplier ! value for controlling grid subdivision  
  real :: amrGridSize          ! length of each side of the (cubic) grid 
  real(double) :: smallestCellSize     ! size of smallest cell
  real(double) :: amrGridCentreX       ! x-coordinate of grid centre 
  real(double) :: amrGridCentreY       ! y-coordinate of grid centre 
  real(double) :: amrGridCentreZ       ! z-coordinate of grid centre 
  logical :: doSmoothGrid   ! whether to correct large differences in the size
                            !   of adjacent grid cells
  logical :: doSmoothGridTau! smooth according to chris's algorithm
  real :: smoothFactor      ! maximum ratio between adjacent cell sizes before
                            !   smoothing is applied 
  real :: sampleFreq        ! maximum number of samples made per subcell
  logical :: amr2dOnly      ! only use cells in 2D plane through grid
  logical :: cylindrical    ! 3d grid is cylindrical
  logical :: spherical      ! spherical 1d grid
  real(double) :: alphaViscosity
  logical :: captureShocks !shock capturing
  logical :: refineOnTemperature !refine grid using temperature gradient
  logical :: refineOnMass !refine grid using cell mass
  logical :: refineOnRhoe !refine grid using cell rhoe
  logical :: refineOnJeans !refine grid using cell mass vs local jeans mass
  logical :: refineOnSpeed !refine grid using cell speed
  real(double) :: massTol !cell mass tolerance
  logical :: refineOnIonization !refine grid using ionization gradient
  real(double) :: amrTolerance !maximum gradient before AMR grid refines
  real(double) :: amrUnrefineTolerance !minimum gradient before AMR grid unrefines
  real(double) :: amrTemperatureTol !maximum temperature gradient before AMR grid refines 
  real(double) :: amrSpeedTol !maximum speed grad before AMR grid refines
  real(double) :: amrIonFracTol !maximum ion frac grad before AMR grid refines
  real(double) :: amrRhoeTol !Maximum rhoe grad before AMR grid refines

  ! Grid smoothing based on optical depth
  real :: lambdaSmooth
  real :: tauSmoothMax, tauSmoothMin

!---------------------------------------------------
! Parameters for setting up a run from SPH particles
!---------------------------------------------------

  character(len=80), protected :: sphdataFilename
  character(len=80), protected :: inputFileFormat 
  real, protected    :: hcritPercentile
  real, protected    :: hmaxPercentile
  real, protected    :: sph_norm_limit
  integer, protected :: kerneltype
  logical, protected :: useHull
  logical, protected :: refineCentre  ! switch on extra grid refinement for SPH-Torus discs 
  logical, protected :: SphOnePerCell ! Split to one particle per cell for galactic plane survey
  logical :: doVelocitySplit ! Should grid be split based on velocity values of SPH particles? 
  logical, protected :: convertRhoToHI ! Convert density to HI
  integer, protected :: ih2frac        ! column of SPH file which contains H2 fraction
  logical, protected :: sphwithchem    ! SPH has chemistry data which needs to be read
  logical, protected :: discardSinks   ! Don't store sink particles
  logical :: guessNe                   !guess the electron number density based on temperature

!------------------
! Other parameters 
!------------------

! Parameters specific to phaseloop_mod
  real :: lambdatau  ! Wavelength at which testing optical depth computed [A]
  logical :: fastIntegrate 
  logical :: noScattering
  logical :: forceFirstScat
  integer :: nPhase
  integer :: nStartPhase, nEndPhase
  real    :: phaseTime ! time of each phase of simulation (seconds)
  type(VECTOR) :: bondiCentre

! Parameters which control Torus behaviour
  logical           :: debug
  logical           :: suppressWarnings
  logical           :: useBinaryXMLVTKfiles
  logical           :: noVtkGrid        ! Don't write out VTK files of the grid 
  logical           :: vtkIncludeGhosts ! include ghosts in VTK output
  logical           :: parallelVTUFiles
  character(len=80) :: absolutePath
  integer(bigInt)   :: maxMemoryAvailable
  logical           :: blockHandout ! Enable MPI block handout
  character(len=10) :: geometry
  integer(kind=bigInt) :: nPhotons ! number of photons to use (phaseloop, photoion, photoionAMR, timedep)
  logical :: radPressureTest ! perform on the spot absorption for radiation pressure tests
  logical :: UV_vector
! Other physical parameters
  real    :: vturb        ! Subsonic turbulent velocity
  real    :: TMinGlobal   ! globally applied minimum temperature

! Some other parameters   
  character(len=10) :: object
  integer :: nr, nphi
  logical :: lineEmission
  logical :: smoothInnerEdge
  logical :: opaqueCore
  logical :: thin_disc_on       ! T to include disc
  logical :: pencilBeam
  logical :: useBias
  logical :: thinLine
  logical :: doRaman ! raman scattering model 
  real    :: lamLine 
  real    :: massEnvelope

