module gridtype_mod

! this module contains the description of the 3d grid of opacity arrays.
! the definition has been moved from the grid_mod module to resolve 
!   some dependancy problems. nhs

  use kind_mod
  use vector_mod                      ! vector math
  use octal_mod                       ! octal type for amr  
  use gaussian_mod

  implicit none

  integer, parameter :: statEqMaxLevels =   14   ! number of stateq levels
!  integer, parameter :: statEqMaxLevels =   6    ! number of stateq levels
                                                 ! (specified as a parameter for speed)

  public

  type GRIDTYPE
     integer :: nx                               ! size of cartesian grid in x
     integer :: ny                               ! size of cartesian grid in y
     integer :: nz                               ! size of cartesian grid in z
     integer :: nr                               ! size of polar grid in r
     integer :: nmu                              ! size of polar grid in mu
     integer :: nphi                             ! size of polar grid in phi
     integer :: na1, na2, na3                    ! generic 3d grid sizes
     integer :: nlambda                          ! number of wavelength points
     integer :: nopacity                         ! number of opacity array (kappaAbs and KappaSca)
     logical :: flatspec                         ! flat spectrum being computed
     logical :: adaptive                         ! are the cells adaptive, rather than fixed
     logical :: cartesian                        ! is the grid cartesian?
     logical :: polar
     logical :: isotropic                        ! are the axes evenly spaced?
     logical :: hitcore
     real :: diskRadius
     logical :: oneKappa                        
     real, pointer :: oneKappaAbs(:,:) => null()
     real, pointer :: oneKappaSca(:,:) => null()
     real :: kappaTest
     integer :: itestlam
     type(VECTOR) :: diskNormal
     real :: dipoleOffset
     character(len=20) :: geometry               ! type of geometry
     type(vector), pointer :: dVbydr(:,:,:)            ! velocity gradient

     real(double),pointer :: oneProbLine(:)
     real(double),pointer :: oneProbCont(:)
     integer :: nProb
     integer, pointer :: cellIndex(:,:)
     real, pointer :: rho(:,:,:) => null()       ! density grid
     real, pointer :: kappaAbs(:,:,:,:) => null()! cont absorption opacities
     real, pointer :: kappaSca(:,:,:,:) => null()! scattering opacities
     real, pointer :: kappaAbsRed(:,:,:,:) => null()! cont absorption opacities
     real, pointer :: kappaScaRed(:,:,:,:)=> null() ! scattering opacities
     real, pointer :: chiLine(:,:,:) => null()   ! line opacity
     real, pointer :: etaLine(:,:,:) => null()   ! line emissivity
     real, pointer :: etaCont(:,:,:) => null()   ! line emissivity
     real, pointer :: sigma(:,:,:)  => null()    ! radial velocity gradient
     type(VECTOR), pointer :: velocity(:,:,:) => null() ! velocity grid
     real, pointer :: rAxis(:) => null()                  ! r-axis

     real, pointer :: biasCont(:) => null()      ! radial bias for cont
     real, pointer :: biasLine(:) => null()      ! radial bias for line

     real, pointer :: muAxis(:) => null()                 ! mu-axis
     real, pointer :: phiAxis(:)                 ! phi-axis
     real          :: dMu, dPhi                  ! spacing of mu and phi axes
     real, pointer :: xAxis(:) => null()         ! x-axis
     real, pointer :: yAxis(:) => null()         ! y-axis
     real, pointer :: zAxis(:) => null()         ! z-axis
     real, pointer :: lamArray(:) => null()      ! the wavelength array
     real, pointer :: rProbDistLine(:) => null() ! emissivity probability
     real, pointer :: muProbDistLine(:,:) => null()   ! distributions
     real, pointer :: phiProbDistLine(:,:,:) => null()!
     real, pointer :: xProbDistLine(:) => null()      ! emissivity probability
     real, pointer :: yProbDistLine(:,:) => null()    ! distributions
     real, pointer :: zProbDistLine(:,:,:) => null()  !
     real, pointer :: rProbDistCont(:) => null()      ! emissivity probability
     real, pointer :: muProbDistCont(:,:) => null()   ! distributions
     real, pointer :: phiProbDistCont(:,:,:) => null()!
     real, pointer :: xProbDistCont(:) => null()      ! emissivity probability
     real, pointer :: yProbDistCont(:,:) => null()    ! distributions
     real, pointer :: zProbDistCont(:,:,:) => null()  !
     real, pointer :: temperature(:,:,:) => null()    ! grid cell temperatures
     real, pointer :: biasLine3D(:,:,:) => null()     ! grid bias distribution
     real, pointer :: biasCont3D(:,:,:) => null()     ! grid bias distribution
     real :: rCore                                   ! core radius
     real(double) :: lCore = -9.9e9         ! core luminosity
     real :: chanceWindOverTotalContinuum = -9.9e9   ! chance of continuum photon being produced
     ! in the wind rather than at the core

     logical :: lineEmission                        ! is there line emission?
     logical :: contEmission                        ! is there continuum emission?

     logical :: doRaman                             ! is this a raman-scattering model?
     logical :: resonanceLine

     real :: rStar1, rStar2
     real :: lumRatio = -9.9e9

     integer :: nDustType
     real :: densityScaleFac

     real :: tempSource = -9.9e9

     type(VECTOR) :: starPos1, starPos2

     real :: lambda2 = -9.9e9
     real :: rInner = -9.9e9
     real :: rOuter = -9.9e9


     real(double), pointer :: N(:,:,:,:) => null() ! stateq level pops
     real(double), pointer :: Ne(:,:,:) => null()  ! electron density
     real(double), pointer :: nTot(:,:,:) => null()! total density
     logical(kind=1), pointer :: inStar(:,:,:) => null()
     logical(kind=1), pointer :: inUse(:,:,:) => null()

     integer :: maxLevels                            ! number of stateq levels 
     
     ! adaptive mesh refinement stuff
     type(octal), POINTER :: octreeRoot => NULL()    ! root of adaptive mesh octree 
     integer :: maxDepth                             ! maximum depth of octree
     real(oct) :: halfSmallestSubcell 
       ! 'halfSmallestSubcell' is half the vertex length of the grid's smallest
       !    subcells. we store this because it is used frequently in calculations.
     integer :: nOctals                              ! total number of octals 
     real    :: smoothingFactor                      ! inter-cell maximum ratio before smooth
     integer :: ng
     type(GAUSSIAN), pointer :: gArray(:)

     logical :: statEq2d                             ! can do statEq in 2-D
     logical :: amr2dOnly                            ! only do statEq in 2-D
     !
     REAL    :: timeNow ! elapsed time since start of model   ! [seconds]
     

  end type GRIDTYPE

end module gridtype_mod
