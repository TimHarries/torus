module input_variables
  ! a large number of variables are passed to the torus main
  !   program by the inputs subroutine. they are all collected
  !   here so that they only need to be defined once.  nhs

  use kind_mod
  use vector_mod, only: vector
  
  implicit none

  public

 ! globally applied minimum temperature
  real    :: TMinGlobal     

  ! variables for the grid

  logical :: debug
  logical :: suppressWarnings
  integer :: idump ! hydrodynamics time step dump number
  logical :: dumpInnerEdge
  logical :: doSelfGrav
  real(double) :: gridDistanceScale
  real :: cflNumber
  integer :: minDepthAMR, maxDepthAMR
  character(len=10) :: geometry
  character(len=10) :: object
  real :: bigomega, eddingtonGamma, alphaCAK
  logical :: uniformStar
  character(len=80) :: absolutePath
  integer :: nx,ny,nz
  integer :: nr, nmu, nphi
  logical :: photoionization
  logical :: quickThermal
  logical :: dustOnly
  logical :: molecular
  logical :: dumpdi
  logical :: h21cm 
  real :: beamSize
  real :: eccentricity
  real(double) :: rotateViewAboutX, rotateViewAboutY, rotateViewAboutZ
  logical :: internalView ! observer is internal to the galaxy
  real(double) :: intPosX, intPosY, intPosZ ! Position of internal observer
  real(double) :: intDeltaVx, intDeltaVy, intDeltaVz ! Additional velocity applied to internal observer
  logical :: cmf, sobolev
  integer :: iTransLine, iTransAtom
  integer :: nAtom
  character(len=20) :: atomFilename(10)
  integer :: nClumps
  character(len=80) :: distortionType
  integer :: nPhase
  integer :: nStartPhase, nEndPhase
  real    :: phaseTime ! time of each phase of simulation (seconds)
  real :: tolerance ! maximum acceptable fractional change for J levels between iterations in molecular_mod
  integer :: setmaxlevel ! Subsonic turbulent velocity
  real :: vturb ! Subsonic turbulent velocity
  logical :: noturb ! Subsonic turbulent velocity
  logical :: lineEmission
  logical :: smoothInnerEdge
  logical :: lineimage
  logical :: quasi ! use quasirandom numbers
  logical :: dongstep ! controls Ng Acceleration
  logical :: densitysubsample ! do you want to subsample density in your images?
  logical :: lowmemory ! if memory is tight
  logical :: maxrhocalc
  logical :: resonanceLine
  real :: probContPhoton
  ! If T, the first scattering is forced in torusMain. (Useful for optically thin case.)
  logical :: forceFirstScat  
  logical :: doIntensivePeelOff
  logical :: fastIntegrate
  real(double) :: etaViscosity
  logical :: opticallyThickContinuum, onTheSpot
  logical :: sed, jansky, SIsed, inArcsec, dustySED

  ! Data cube parameters
  real(double) :: galaxyInclination, galaxyPositionAngle 
  real(double) :: dataCubeVelocityOffset ! Velocity offset for data cube
  logical :: SplitCubes ! Split cube into +ve and -ve contributions? 

  ! variables to do with dust

  integer :: nInclination  ! number of inclinations
  real :: firstInclination ! first inclination angle
  real :: lastInclination  !
  real :: thisInclination
  real, allocatable :: inclinations(:)

  ! variables to do with dust
  
  integer, parameter :: maxDustTypes = 10
  character(len=80) :: grainType(maxDustTypes) ! sil_ow, sil_oc, sil_dl, amc_hn, sic_pg, gr1_dl, gr2_dl  
  real :: grainFrac(maxDustTypes)
  real :: grainSize
  logical :: mie
  logical :: storescattered
  real :: scatteredLightWavelength
  logical :: readMiePhase
  logical :: writeMiePhase
  logical :: useOldMiePhaseCalc
  logical :: includeGasOpacity
  logical :: isotropicScattering
  real :: dusttogas
  logical :: dustfile
  logical :: variableDustSublimation
  real :: probDust
  integer :: nDustType
  logical :: forcedWavelength
  real :: usePhotonWavelength
  logical :: useDust, realdust ! molecular_mod includes continuum emission
  logical :: doCOchemistry
  logical :: isinlte ! assume grid is in LTE
  real :: r0, rhoC ! Filamentry parameters
  real :: x_D

  real(double) :: minVel, maxVel
  real :: imageside ! molecular_mod image parameters
  integer :: itrans, npixels, nsubpixels !nv already exists
  real(double) :: centrevecx, centrevecy, centrevecz
  logical :: wanttau

  logical :: blockHandout


  real :: diffDepth  ! depth of diffusion zone - in units of rosseland optical depth
  !
  ! abundances of different types of dust grains. These will be used when 
  ! the graintype assigned is "mixed."
  integer :: ngrain 
  real :: X_grain(10)    ! abundaunce 
  character(LEN=30) :: grainname(10)

  real :: h_abund, he_abund, c_abund, n_abund, &
       o_abund, ne_abund, s_abund  



  ! torus images

  real :: slitPA, slitWidth, slitLength
  real :: imageSizeinArcsec
  real :: vfwhm, pfwhm, vSys
  integer :: nSlit, np, nv
  type(VECTOR) :: slitPosition1, slitPosition2
  logical :: stokesImage
  real :: setImageSize
  real :: imageScale
  real :: vMin, vMax
  real :: gridDistance
  integer :: observerpos

  integer :: npix    ! Number of pixels for polimages  
  character(LEN=30) :: filter_set_name  ! name of filter set used for images
  ! if T, the dimension of the images will be in arcsec otherwise in phyiscal unit of length.
  logical :: imageInArcsec  

  logical :: forceRotate
  logical :: forceNoRotate

  ! variables to do with dust
  !
  integer :: nSpiral
  ! size distribution of dust grain is now assumed to have 
  ! the following form:
  !  
  !   n(a) =  const * a^-q * e^((-a/a0)^p)
  real :: aMin(maxDustTypes)    !  The maximun size in microns. 
  real :: aMax(maxDustTypes)    !  The minimum size in microns.
  real :: qDist(maxDustTypes)   !  q exponet in the equation above.
  real :: a0(maxDustTypes)      !  scale length in the equation above.
  real :: pDist(maxDustTypes)   !  p exponet in the equation above.


  ! Flag to include accreation disc in sph dust model or not
  logical :: disc_on

  ! restrcting a calculation to a certain star in cluster geometry
  integer :: idx_restrict_star 

  ! whitney stuff

  real :: erInner, erOuter, drInner, drOuter, cavAngle, cavDens
  real :: rStellar, mdotEnv, mEnv

  ! variables for clumped wind models

  integer :: nBlobs
  logical :: freshBlobs
  real :: blobContrast


  ! filenames
  character(len=80) :: outFile
  character(len=80) :: intProFilename
  character(len=80) :: intProFilename2
  character(len=80) :: opacityDataFile
  character(len=80) :: dustfilename(10)
  logical :: useNdf 
  

  logical :: plot_maps
  logical :: plotlevels
  logical :: writetempfits

  ! output arrays

  integer :: nLambda, nLambdaInput
  real :: lamStart, lamEnd
  logical :: lamLinear
  logical :: oneKappa
  logical :: lamfile
  character(len=80):: lamfilename

  real :: lambdatau  ! Wavelength at which testing optical depth computed [A]
  real :: lambdaSmooth
  real :: tauSmoothMax, tauSmoothMin
  real :: tauRad
  real :: tauDiff, tauForce
  logical :: resetDiffusion
  real :: eDensTol

  ! variables for a second source of radiation

  logical :: secondSource                    ! second photon source?
  type(VECTOR) :: secondSourcePosition       ! the position of it
  real :: binarySep
  real :: momRatio


  real :: rpower ! radial density power  r^-rpower

  ! model flags
  
  logical :: fillTio
  logical :: opaqueCore
  logical :: thin_disc_on       ! T to include disc
  logical :: pencilBeam
  logical :: plezModelOn
  logical :: useBias
  logical :: fillThomson
  logical :: screened
  logical :: thinLine
  logical :: VoigtProf


  ! model parameters
  
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
  real    :: isoThermTempLogical    ! If isoTherm is true, use this temperature in [K]

  ! pp disk stuff
  real :: rSmooth, rHeight, sigmaPower
  real :: flaringPower, gapViscAlpha
  real :: rGap, mPlanet, gapWidth
  logical :: planetGap


  ! single dust blob parameters (WR137 type model)

  real :: dustBlobDistance, phiDustBlob
  real :: xDustBlobSize, yDustBlobSize, zDustBlobSize

  real :: massEnvelope

  ! raman scattering model parameters

  logical :: doRaman
  real :: vo6
  real :: o6width
  character(len=20) :: ramanDist          ! raman distortion type
  real :: ramVel


  ! core emission line parameters

  logical :: coreEmissionLine
  real :: velWidthCoreEmission


  real :: vContrast

  ! misc

  character(len=80) :: misc

  character(len=80) :: contFluxFile, contFluxFile1, contFluxFile2
  real :: phaseOffset
  logical :: lineOff, starOff
  real :: logMassLossRate
  real :: lamLine
  character(LEN=20) :: ion_name  ! Name of the ion (see opacity_lte_mod.f90 for the list of a valid name.)
  real ::  ion_frac              ! n_ion/n_specie
  integer :: nLower, nUpper

  

 ! Spot stuff
  
  integer :: nSpot                       ! number of spots
  real :: fSpot                          ! factional area coverage of spots
  real :: tSpot                          ! spot temperatures
  real :: thetaSpot, phiSpot             ! spot coords
  logical :: photLine                    ! photospheric line production


! SPH parameters
  character(len=80) :: sphdataFilename
  real :: hcritPercentile
  real :: hmaxPercentile
  real :: sph_norm_limit
  integer :: kerneltype

  logical :: refineCentre  ! switch on extra grid refinement for SPH-Torus discs 

  ! binary parameters

  real :: period
  real :: shockWidth, shockFac
  real :: vNought1, vNought2
  real :: beta1, beta2
  real :: vTerm1 , vTerm2
  real :: radius1, radius2
  logical :: readPops, writePops
  logical :: readPhasePops, writePhasePops
  logical :: readLucy, writeLucy, openLucy
  logical :: readMol, writeMol ! molecular_mod equivalents of read/writelucy
  real :: tThresh
  logical :: redolucy
  character(len=80) :: lucyFilenameIn
  character(len=80) :: lucyFilenameOut
  character(len=80) :: molFilename, molFilenameIn, molFilenameOut
  character(len=80) :: moleculefile
  logical :: restart
  logical :: addnewmoldata
  real :: mass1, mass2, massRatio, streamFac
  real :: temp1, temp2
  character(len=80) :: popFilename
  real :: deflectionAngle
  logical :: lte
  logical :: LyContThick
  logical :: curtains, enhance
  real :: dipoleOffset
  logical :: amr1d, amr2d, amr3d

  ! T Tauri parameters ----------------------------------------------------
  character(len=80) :: MdotType ! variable accretion rate model in use
  real :: MdotParameter1, MdotParameter2, MdotParameter3, MdotParameter4
  real :: MdotParameter5, MdotParameter6
  real :: TTauriRinner, TTauriRouter ! disc sizes (in R_star units)
  logical :: ttauriwarp, ttauriwind, ttauridisc
  real :: TTauriRstar ! stellar radius (in R_sol units)
  real :: TTauriMstar ! stellar mass   (in M_sol units)
  real :: TTauriDiskHeight ! (in R_star units)
  real :: TTauriDiskRin    ! (in R_star units)
  real(double) :: hOverR
  real :: ThinDiskRin      ! (in R_star units)
  real :: curtainsPhi1s ! accretion curtains from (s)tart... 
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

  logical :: pointSource

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


  ! adaptive mesh stuff ---------------------------------------------------------
  logical :: gridUsesAMR    ! true if grid is adaptive
  real(double) :: limitScalar  ! value for controlling grid subdivision 
  real(double) :: limitScalar2 ! value for controlling grid subdivision 
  real(double) :: vturbmultiplier ! value for controlling grid subdivision  
  logical      :: doVelocitySplit ! Should grid be split based on velocity values of SPH particles? 
  real :: amrGridSize          ! length of each side of the (cubic) grid 
  real(double) :: amrGridCentreX       ! x-coordinate of grid centre 
  real(double) :: amrGridCentreY       ! y-coordinate of grid centre 
  real(double) :: amrGridCentreZ       ! z-coordinate of grid centre 
  logical :: doSmoothGrid   ! whether to correct large differences in the size
                            !   of adjacent grid cells
  logical :: doSmoothGridTau! smooth according to chris's algorithm
  real :: smoothFactor      ! maximum ratio between adjacent cell sizes before
                            !   smoothing is applied 
  logical :: suppressLucySmooth ! Suppress grid smoothing in lucy_mod? 
  real :: sampleFreq        ! maximum number of samples made per subcell
  logical :: readFileFormatted  ! whether 'grid' input  file is formatted
  logical :: writeFileFormatted ! whether 'grid' output file is formatte
  logical :: statEq2d       ! whether statEq can be run in 2-D
  logical :: noPhaseUpdate  ! disable updating AMR grid at each phase
  logical :: amr2dOnly      ! only use cells in 2D plane through grid
  logical :: cylindrical
  logical :: forceLineChange ! recalculate opacities for new transition 
  logical :: statEq1stOctant ! T if do stateqAMR routine for the octals in first octant
                             ! then later values are mapped to other octants.

  integer(kind=bigInt) :: nPhotons
  integer :: maxScat
  logical :: noScattering
  

  ! sphericity test

  logical :: sphericityTest
 

  ! some variables had different names in inputs_mod and torusMain. To avoid conflict 
  !   and to retain backwards compatibility, we will adopt the name used in torusMain.
  !   The names are changed in inputs_mod, but keep the old labels so that the  
  !   input files don't have to be changed. 
  
  logical :: fillRayleighOpacity       ! previously: 'fillRayleigh'
  character(len=10) :: gridCoords      ! previously: 'gridtype'
  logical :: doPVimage                 ! previously: 'pvimage'
  real :: relIntCoreEmissionLine       ! previously: 'relIntCoreEmission'
  real :: velWidthCoreEmissionLine     ! previously: 'velWidthCoreEmission'

  
  ! lucy radiative equ

  logical :: lucyRadiativeEq
  logical :: solveVerticalHydro
  integer :: nHydro
  integer :: nLucy
  integer :: iterLucy
  logical :: narrowBandImage
  logical :: useInterp
  real    :: lucy_undersampled  
  integer :: minCrossings
  logical :: forceLucyConv
  logical :: multiLucyFiles

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

  real :: molAbundance

!hydro stuff
  logical :: hydrodynamics
  integer :: xminusbound, yminusbound, zminusbound
  integer :: xplusbound, yplusbound, zplusbound
  real :: x1, x2

  logical :: timeDependentRT
  
end module input_variables


