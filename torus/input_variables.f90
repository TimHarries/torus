module input_variables
  ! a large number of variables are passed to the torus main
  !   program by the inputs subroutine. they are all collected
  !   here so that they only need to be defined once.  nhs

  use utils_mod

  use constants_mod
  use vector_mod
  use unix_mod
  
  implicit none

  public

  ! variables for the grid

  character(len=10) :: geometry
  integer :: nx,ny,nz
  integer :: nr, nmu, nphi

  integer :: nClumps
  character(len=80) :: distortionType
  integer :: nPhase
  integer :: nStartPhase, nEndPhase

  logical :: lineEmission
  logical :: resonanceLine
  real :: probContPhoton

  logical :: sed, jansky, SIsed, inArcsec, dustySED

  ! variables to do with dust

  character(len=80) :: grainType ! sil_ow, sil_oc, sil_dl, amc_hn, sic_pg, gr1_dl, gr2_dl
  real :: grainSize
  logical :: mie
  real :: dusttogas
  logical :: dustfile
  real :: probDust
  integer :: nDustType
  logical :: forcedWavelength
  real :: usePhotonWavelength


  !
  ! abundances of different types of dust grains. These will be used when 
  ! the graintype assigned is "mixed."
  integer, parameter :: ngrain  = 7 ! number of grain types implemented.
  real :: X_grain(ngrain)    ! abundaunce 
  character(LEN=30) :: grainname(ngrain)
  data grainname /"sil_ow", "sil_oc", "sil_dl", "amc_hn", &
       &          "sic_pg", "gr1_dl", "gr2_dl"/
  ! Note ::
  !   X_grain(1) => sil_ow 
  !   X_grain(2) => sil_oc
  !   X_grain(3) => sil_dl
  !   X_grain(4) => amc_hn
  !   X_grain(5) => sic_pg
  !   X_grain(6) => gr1_dl
  !   X_grain(7) => gr2_dl
  



  ! torus images

  real :: slitPA, slitWidth, slitLength
  real :: imageSizeinArcsec
  real :: vfwhm, pfwhm, vSys
  integer :: nSlit, np, nv
  type(VECTOR) :: slitPosition1, slitPosition2
  logical :: stokesImage
  real :: setImageSize
  real :: vMin, vMax
  real :: gridDistance

  integer :: npix    ! Number of pixels for polimages  
  character(LEN=30) :: filter_set_name  ! name of filter set used for images

  ! variables to do with dust
  !
  integer :: nSpiral
  ! size distribution of dust grain is now assumed to have 
  ! the following form:
  !  
  !   n(a) =  const * a^-q * e^((-a/a0)^p)
  real :: aMin    !  The maximun size in microns. 
  real :: aMax    !  The minimum size in microns.
  real :: qDist   !  q exponet in the equation above.
  real :: a0      !  scale length in the equation above.
  real :: pDist   !  p exponet in the equation above.


  ! Flag to include accreation disc in sph dust model or not
  logical :: disc_on

  ! restrcting a calculation to a certain star in cluster geometry
  integer :: idx_restrict_star 

  ! variables for clumped wind models

  integer :: nBlobs
  logical :: freshBlobs
  real :: blobContrast


  ! filenames
  character(len=80) :: outFile
  character(len=80) :: intProFilename
  character(len=80) :: intProFilename2
  character(len=80) :: device, opacityDataFile
  character(len=80) :: dustfilename(10)
  logical :: useNdf 
  
  ! plane to be plotted by plot_AMR_planes and plot_AMR_values
  ! Choose one of the following: "x-y", "y-z", "z-x" or "x-z"
  logical :: plot_maps 
  character(len=3)  :: plane_for_plot  
  logical :: show_value_3rd_dim    ! If T, the value of the third dimension will be shown.


  ! output arrays

  integer :: nLambda
  real :: lamStart, lamEnd
  logical :: lamLinear
  logical :: oneKappa
  logical :: lamfile
  character(len=80):: lamfilename

  real :: lambdatau  ! Wavelength at which testing optical depth computed [A]
  real :: lambdaSmooth
  real :: tauSmoothMax, tauSmoothMin
  real :: tauRad

  ! variables for a second source of radiation

  logical :: secondSource                    ! second photon source?
  type(VECTOR) :: secondSourcePosition       ! the position of it
  real :: binarySep
  real :: momRatio


  real :: rpower ! radial density power  r^-rpower

  ! model flags
  
  logical :: fillTio
  logical :: opaqueCore
  logical :: pencilBeam
  logical :: plezModelOn
  logical :: useBias
  logical :: fillThomson
  logical :: movie
  logical :: plotVelocity
  logical :: screened
  logical :: thinLine
  logical :: inputOK


  ! model parameters
  
  real :: height
  real :: sigma0
  real :: rMin, rMaj
  real :: shellFrac
  real :: Teff
  real :: radius, kfac, xfac
  real :: inclination
  real :: contrast
  real :: rCore, rInner
  real :: rTorus, rOuter
  real :: rho, rho0
  real :: scale
  real :: mCore, diskTemp, mDisc
  real :: scaleDensity
  real :: vRot
  real :: mdot
  real :: beta
  real :: vterm
  real :: v0


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

  character(len=80) :: contFluxFile, contFluxFile2
  real :: phaseOffset
  logical :: lineOff
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


  ! binary parameters

  real :: period
  real :: shockWidth, shockFac
  real :: vNought1, vNought2
  real :: beta1, beta2
  real :: vTerm1 , vTerm2
  real :: radius1, radius2
  real :: mdot1, mdot2
  logical :: readPops, writePops
  logical :: readLucy, writeLucy
  real :: tThresh
  character(len=80) :: lucyFilename
  logical :: redolucy
  character(len=80) :: lucyFilenameIn
  character(len=80) :: lucyFilenameOut
  real :: mass1, mass2
  real :: temp1, temp2
  character(len=80) :: popFilename
  real :: deflectionAngle
  logical :: lte
  logical :: curtains, enhance
  real :: dipoleOffset
  logical :: twoD

  ! adaptive mesh stuff 
  logical :: gridUsesAMR    ! true if grid is adaptive
  logical :: amr2d          ! a two-d AMR grid only
  real(kind=doubleKind) :: limitScalar  ! value for controlling grid subdivision 
  real(kind=doubleKind) :: limitScalar2 ! value for controlling grid subdivision 
  real :: amrGridSize          ! length of each side of the (cubic) grid 
  real(kind=doublekind) :: amrGridCentreX       ! x-coordinate of grid centre 
  real(kind=doublekind) :: amrGridCentreY       ! y-coordinate of grid centre 
  real(kind=doublekind) :: amrGridCentreZ       ! z-coordinate of grid centre 
  logical :: doSmoothGrid   ! whether to correct large differences in the size
                            !   of adjacent grid cells
  logical :: doSmoothGridTau! smooth according to chris's algorithm
  real :: smoothFactor      ! maximum ratio between adjacent cell sizes before
                            !   smoothing is applied 
  real :: sampleFreq        ! maximum number of samples made per subcell
  logical :: readFileFormatted  ! whether 'grid' input  file is formatted
  logical :: writeFileFormatted ! whether 'grid' output file is formatted


  integer :: nPhotons
  integer :: maxScat
  

  ! sphericity test

  logical :: sphericityTest
 

  ! some variables had different names in inputs_mod and torusMain. To avoid conflict 
  !   and to retain backwards compatibility, we will adopt the name used in torusMain.
  !   The names are changed in inputs_mod, but keep the old labels so that the  
  !   input files don't have to be changed. 
  
  real :: inputKappaSca, inputKappaAbs ! previously: 'kappaSca' and 'kappaAbs' 
  logical :: fillRayleighOpacity       ! previously: 'fillRayleigh'
  character(len=10) :: gridCoords      ! previously: 'gridtype'
  logical :: doPVimage                 ! previously: 'pvimage'
  real :: relIntCoreEmissionLine       ! previously: 'relIntCoreEmission'
  real :: velWidthCoreEmissionLine     ! previously: 'velWidthCoreEmission'

  
  ! lucy radiative equ

  logical :: lucyRadiativeEq
  logical :: solveVerticalHydro
  integer :: nLucy
  logical :: narrowBandImage
  logical :: useInterp
  real    :: lucy_undersampled  


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
  
end module input_variables


