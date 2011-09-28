  !   a large number of variables are passed to the torus main
  !   program by the inputs subroutine. they are all collected
  !   here so that they only need to be defined once.  nhs
  !   These now are inlcuded in inputs_mod. D. Acreman 

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
  logical :: hydrodynamics
  logical :: timeDependentRT

!-------------------------------
! Type of output to be generated
!-------------------------------

  logical :: calcDataCube
  logical :: calcPhotometry
  logical :: calcImage
  logical :: calcSpectrum
  logical :: calcBenchmark

!-----------------------------------
! Write/read a grid for a warm start
!-----------------------------------
 
  logical :: readGrid, writeGrid  ! Do we read/write the AMR grid from/to a file
  logical :: justDump             !Dump a vtk file for the read in grid and exit
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
  logical :: usemetals         !Include species heavier than Helium
  logical :: checkForPhoto     !Check whether or not a photoionization loop is necessary
  logical :: monochromatic     !Use a monochromatic radiation field
  logical :: quickThermal      !Use a simplified thermal balance calculation

! Abundances used in ion_mod
  real :: h_abund, he_abund, c_abund, n_abund, o_abund, ne_abund, s_abund  

  ! Parameters  specific to domain decomposed photoionisation 

  !Stack optimization
  logical :: optimizeStack     !Perform run-time photon stack size optimization
  integer :: dStack            !Stack increment for run-time photon stack size optimization 
  integer :: stackLimit        !Maximum stack size for run-time photon stack size optimization

  integer :: bufferCap         !Number of photon stacks to accomodate in buffer

  logical :: binPhotons        !Dump a spectrum of propagated photon packets
  logical :: noDiffuseField    !Use the on the spot approximation
  logical :: dustOnly          !Consider dust physics only
  integer :: inputNMonte       !Number of photon packets to propagate
  integer :: maxPhotoIonIter   !Maximimum interation number

  logical :: periodicX, periodicY, periodicZ ! Periodic photon boundary conditions
  logical :: doPhotoRefine     !Refine the AMR grid between iterations of the photo loop

!---------------
! Hydrodynamics
!---------------

  logical :: rhieChow                 !Use Rhie-Chow interpolation
  logical :: doSelfGrav               !Do self gravity calculation
  logical :: severeDamping            ! Turn on damping
  logical :: dumpRadial               ! write a text radial cut each dump
  real :: cflNumber                   !Courant-Friedrichs-Lewy constant
  real(double) :: etaViscosity        !Artificial viscosity parameter
  real(double) :: tStart, tEnd, tDump !Start, end and dump interval times
  logical :: hydrovelocityConv        !Convert input velocity vector into simulation velocities 
  logical :: doRefine, doUnrefine     !Allow the AMR grid to refine/unrefine
  logical :: useViscosity             !Use artificial viscosity
  logical :: fluxinterp               !Interpolate fluxes at fine to coarse cell advections
  logical :: fixedRhoBound            !Use fixed density boundary conditions
  real(double) :: rho_const           !Density of fixed density boundary conditions
  character(len=20) :: limiterType    !Flux limiter type
  integer :: idump                    !Hydrodynamics time step dump number
  real(double) :: gridDistanceScale   !Scale of grid

  real(double) :: inflowPressure
  real(double) :: inflowRho
  real(double) :: inflowRhoE
  real(double) :: inflowEnergy
  real(double) :: inflowMomentum
  real(double) :: inflowSpeed

  !Boundary conditions (strings)
  character(len=20) :: xminusboundString, yminusboundString, zminusboundString
  character(len=20) :: xplusboundString, yplusboundString, zplusboundString

  !Boundary conditions (integer codes)
  integer :: xminusbound, yminusbound, zminusbound 
  integer :: xplusbound, yplusbound, zplusbound 

  real :: x1, x2

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
  logical :: constantAbundance
  real :: molAbundance
  character(len=80) :: molFilename
  character(len=80) :: moleculefile
  logical :: restart
  logical :: addnewmoldata
  logical :: dongstep ! controls Ng Acceleration
  logical :: plotlevels   
  logical :: writetempfits
  logical :: doCOchemistry
  logical :: isinlte ! assume grid is in LTE
  logical :: lowmemory ! if memory is tight
  logical :: noturb ! Subsonic turbulent velocity
  real :: tolerance ! maximum acceptable fractional change for J levels between iterations in molecular_mod
  integer :: initnray ! number of rays to use in fixed ray case (stage 1)
  real(double) :: rotateViewAboutX, rotateViewAboutY, rotateViewAboutZ
  integer :: setmaxlevel ! Subsonic turbulent velocity
  logical :: lineimage
  logical :: maxrhocalc
  logical :: quasi ! use quasirandom numbers
  logical :: rgbCube ! reverse velocity axis 
  logical :: densitysubsample ! do you want to subsample density in your images?
  logical :: wanttau ! Output tau cube? (also used in angularImage_mod)
  real(double) :: centrevecx, centrevecy, centrevecz ! cube centre (also in AngularImage_mod)
  real :: beamSize
  logical :: getdepartcoeffs, gettau
  logical :: outputconvergence, dotune
  real :: x_D, x_0 ! CO drop model
  logical ::  realdust ! molecular_mod includes continuum emission
  integer :: nsubpixels ! No. of sub-pixels for ray trace
  integer :: itrans

!------------------------------------------
! Galactic plane survey (angularImage_mod) 
!------------------------------------------

  logical :: h21cm ! also used in molecular_mod
  logical :: internalView ! observer is internal to the galaxy
  logical :: obsVelFromGrid ! Is observer velocity taken from Torus grid?
  logical :: thermalLineWidth ! Use thermal line width? 
  real(double) :: intPosX, intPosY, intPosZ ! Position of internal observer
  real(double) :: intDeltaVx, intDeltaVy, intDeltaVz ! Additional velocity applied to internal observer
  real(double) :: galaxyInclination, galaxyPositionAngle 
  integer :: dssMinSample ! Minimum number of samples per cell when density subsample is used

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
  real :: lamStart, lamEnd
  logical :: oneKappa
  logical :: lamfile              ! Read in wavelengths from a file? 
  character(len=80):: lamfilename ! File to use if lamfile=T

!------------------- 
! Source parameters
!-------------------

  logical :: readSources
  character(len=80) :: sourceFilename
  integer :: inputNsource
  real(double) :: sourceTeff(10), sourceMass(10), sourceRadius(10), sourceProb(10)
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

!-----------------
! Dust parameters
!-----------------

  integer, parameter :: maxDustTypes = 10
  character(len=80) :: grainType(maxDustTypes) ! sil_ow, sil_oc, sil_dl, amc_hn, sic_pg, gr1_dl, gr2_dl  
  real :: grainFrac(maxDustTypes)
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

!----------------
! SED parameters
!----------------

  logical, private :: sed, jansky, SIsed
  real    :: SEDlamMin, SEDlamMax
  logical :: SEDwavLin
  integer :: SEDnumLam

!------------------------------
! Inclinations for SEDs/images 
!------------------------------

  integer, private :: nInclination  ! number of inclinations (phaseloop_mod)
  real, private :: firstInclination ! first inclination angle (phaseloop_mod)
  real, private :: lastInclination=80.0  ! last inclination angle (phaseloop_mod)
  real, allocatable, private :: inclinations(:)
  real :: thisInclination  ! Inclination when atomicPhysics=T (calculateAtomSpectrum and compute_obs_line_flux)

!--------------------------------
! Image and data cube parameters 
!--------------------------------

  real :: positionAngle(10), inclinationArray(10)
  logical :: monteCarloRT
  real(double) :: dataCubeVelocityOffset ! Velocity offset for data cube
  logical :: SplitCubes ! Split cube into +ve and -ve contributions? 
  integer :: FitsBitpix ! bitpix parameter for FITS representation of the data
  real(double) :: minVel, maxVel ! for molecular_mod and angularImage_mod
  real :: vMin, vMax             ! for createLucyImage
  real :: imageside    ! data cube size for molecular_mod and angularImage_mod
  real :: setImageSize ! image size for createLucyImage and phaseloop_mod
  integer :: nv ! number of velocity channels
  integer :: npixels ! number of pixels 
  integer :: nImage ! Number of images to generate
  logical :: stokesImage ! use in phaselooop_mod
  real :: vMinSpec, vMaxSpec ! For atomicDataCube option
  real :: gridDistance ! distance of observer for images, SEDs etc. 
  integer :: observerpos ! position of observer in molecular_mod image generation
  logical :: inclineX, inclineY, inclineZ ! which inclination to use for photoionisation images
  real :: singleInclination ! Inclination angle to use for for inclineX, inclineY, inclineZ
  character(LEN=30) :: filter_set_name  ! name of filter set used for images (phaseloop_mod)
  logical :: imageInArcsec  ! Are images in arcsec or physical units? (image_mod and phaseloop_mod)
  real :: lambdaImage(10)  ! Wavelengths for images. 

  ! file names 
  character(len=80) :: outFile           ! used in phaseloop_mod
  character(len=80) :: datacubeFilename


!----------------------------------------------------------
! Geometry specific parameters (this is a large section !)
!-----------------------------------------------------------

  real(double) :: lXOverLbol ! ??? Sources for Ttauri geometry
  real(double) :: maxCellMass ! ??? Refinement control for Ttauri geometry
  real(double) :: mStarburst, clusterRadius ! ??? for starburst geometry
  character(len=80) :: intextFilename, outtextFilename ! ??? Fogel geometry
  logical :: doSpiral ! For a Shakara Sunyaev Disc
  type(VECTOR) :: sphereVelocity, spherePosition ! unisphere geometry
  real(double) :: sphereMass, sphereRadius       ! unisphere geometry
  real :: rpower ! radial density power  r^-rpower. For proto geometry

  ! whitney stuff
  real :: erInner, erOuter, drInner, drOuter, cavAngle, cavDens
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
  real :: scale, rscale
  real :: mCore, diskTemp, mDisc
  real :: epsilonDisc
  real :: alphaDisc, betaDisc
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
  logical :: ttauriwarp, ttauriwind, ttauridisc
  real :: TTauriRstar ! stellar radius (in R_sol units)
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


  !
  ! Voigt profile prameters
  !
  real :: C_rad    ! Damping constant (radiation)     in [A]
  real :: C_vdw    ! Damping constant (van der Waals) in [A]
  real :: C_stark  ! Damping constant (Stark)         in [A]


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
  integer :: minDepthAMR, maxDepthAMR
  real(double) :: limitScalar  ! value for controlling grid subdivision 
  real(double) :: limitScalar2 ! value for controlling grid subdivision 
  real(double) :: vturbmultiplier ! value for controlling grid subdivision  
  real :: amrGridSize          ! length of each side of the (cubic) grid 
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
  logical :: cylindrical
  logical :: refineOnTemperature !refine grid using temperature gradient
  logical :: refineOnMass !refine grid using cell mass
  logical :: refineOnRhoe !refine grid using cell rhoe
  logical :: refineOnJeans !refine grid using cell mass vs local jeans mass
  logical :: refineOnSpeed !refine grid using cell speed
  real(double) :: massTol !cell mass tolerance
  logical :: refineOnIonization !refine grid using ionization gradient
  real(double) :: amrTolerance !maximum gradient before AMR grid refines

  ! Grid smoothing based on optical depth
  real :: lambdaSmooth
  real :: tauSmoothMax, tauSmoothMin

!---------------------------------------------------
! Parameters for setting up a run from SPH particles
!---------------------------------------------------

  character(len=80) :: sphdataFilename
  character(len=80) :: inputFileFormat 
  real    :: hcritPercentile
  real    :: hmaxPercentile
  real    :: sph_norm_limit
  integer :: kerneltype
  logical :: useHull
  logical :: refineCentre  ! switch on extra grid refinement for SPH-Torus discs 
  logical :: SphOnePerCell ! Split to one particle per cell for galactic plane survey
  logical :: doVelocitySplit ! Should grid be split based on velocity values of SPH particles? 
  logical :: convertRhoToHI ! Convert density to HI
  integer :: ih2frac        ! column of SPH file which contains H2 fraction


!------------------
! Other parameters 
!------------------

! Parameters specific to phaseloop_mod
  real :: lambdatau  ! Wavelength at which testing optical depth computed [A]
  logical :: fastIntegrate 
  logical :: noScattering
  integer :: nPhase
  integer :: nStartPhase, nEndPhase
  real    :: phaseTime ! time of each phase of simulation (seconds)

! Parameters which control Torus behaviour
  logical           :: debug
  logical           :: suppressWarnings
  logical           :: useBinaryXMLVTKfiles
  logical           :: parallelVTUFiles
  character(len=80) :: absolutePath
  integer(bigInt)   :: maxMemoryAvailable
  logical           :: blockHandout ! Enable MPI block handout
  character(len=10) :: geometry
  integer(kind=bigInt) :: nPhotons ! number of photons to use (phaseloop, photoion, photoionAMR, timedep)

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

