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
  
  character(len=80) :: distortionType
  integer :: nPhase
  integer :: nStartPhase, nEndPhase

  logical :: lineEmission
  real :: probContPhoton


  ! variables to do with dust
  
  character(len=80) :: grainType
  real :: grainSize
  logical :: mie


  ! torus images

  real :: slitPA, slitWidth, slitLength
  real :: vfwhm, pfwhm, vSys
  integer :: nSlit, np, nv
  type(VECTOR) :: slitPosition1, slitPosition2
  logical :: stokesImage
  real :: vMin, vMax
  real :: gridDistance
  
  
  ! variables to do with dust
  
  integer :: nSpiral
  real :: aMin, aMax, qDist


  ! variables for clumped wind models

  integer :: nBlobs
  logical :: freshBlobs
  real :: blobContrast


  ! filenames
  character(len=80) :: outFile
  character(len=80) :: intProFilename
  character(len=80) :: intProFilename2
  character(len=80) :: device, opacityDataFile
  logical :: useNdf 


  ! output arrays

  integer :: nLambda
  real :: lamStart, lamEnd
  logical :: lamLinear
 
 
  ! variables for a second source of radiation

  logical :: secondSource                    ! second photon source?
  type(VECTOR) :: secondSourcePosition       ! the position of it
  real :: lumRatio                           ! lumonsity ratio
  real :: binarySep
  real :: momRatio


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
  real :: rMin, rMaj
  real :: shellFrac
  real :: Teff
  real :: radius, kfac, xfac
  real :: inclination
  real :: contrast
  real :: rCore, rInner
  real :: rTorus, rOuter
  real :: rho
  real :: scale
  real :: mCore, diskTemp
  real :: scaleDensity
  real :: vRot
  real :: mdot
  real :: beta
  real :: vterm
  real :: v0


  ! single dust blob parameters (WR137 type model)

  real :: dustBlobDistance, phiDustBlob
  real :: xDustBlobSize, yDustBlobSize, zDustBlobSize


  ! raman scattering model parameters

  logical :: doRaman
  real :: vo6
  real :: o6width
  character(len=20) :: ramanDist          ! raman distortion type
  real :: ramVel


  ! core emission line parameters

  logical :: coreEmissionLine
  real :: velWidthCoreEmission


  ! misc

  character(len=80) :: misc

  character(len=80) :: contFluxFile, contFluxFile2
  real :: phaseOffset
  logical :: lineOff
  real :: logMassLossRate
  real :: lamLine
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
  real :: mass1, mass2
  real :: temp1, temp2
  character(len=80) :: popFilename
  real :: deflectionAngle
  logical :: lte
  logical :: curtains, enhance
  real :: dipoleOffset

  ! adaptive mesh stuff 
  logical :: gridUsesAMR    ! true if grid is adaptive
  real :: limitScalar       ! value for controlling grid subdivision 
  real :: amrGridSize          ! length of each side of the (cubic) grid 
  real :: amrGridCentreX       ! x-coordinate of grid centre 
  real :: amrGridCentreY       ! y-coordinate of grid centre 
  real :: amrGridCentreZ       ! z-coordinate of grid centre 
  logical :: doSmoothGrid   ! whether to correct large differences in the size
                            !   of adjacent grid cells
  real :: smoothFactor      ! maximum ratio between adjacent cell sizes before
                            !   smoothing is applied 
  real :: sampleFreq        ! maximum number of samples made per subcell


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

  logical :: narrowBandImage
  logical :: useInterp
  
end module input_variables
