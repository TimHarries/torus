module gridtype_mod

! this module contains the description of the 3d grid of opacity arrays.
! the definition has been moved from the grid_mod module to resolve 
!   some dependancy problems. nhs

  use kind_mod
  use vector_mod                      ! vector math
  use octal_mod                       ! octal type for amr  

  implicit none

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
     logical :: flatspec                         ! flat spectrum being computed
     logical :: adaptive                         ! are the cells adaptive, rather than fixed
     logical :: cartesian                        ! is the grid cartesian?
     logical :: isotropic                        ! are the axes evenly spaced?
     logical :: hitcore
     real :: diskRadius
     type(VECTOR) :: diskNormal
     real :: dipoleOffset
     character(len=20) :: geometry               ! type of geometry
     real, pointer :: rho(:,:,:)                 ! density grid
     real, pointer :: kappaAbs(:,:,:,:)          ! cont absorption opacities
     real, pointer :: kappaSca(:,:,:,:)          ! scattering opacities
     real, pointer :: kappaAbsRed(:,:,:,:)          ! cont absorption opacities
     real, pointer :: kappaScaRed(:,:,:,:)          ! scattering opacities
     real, pointer :: chiLine(:,:,:)             ! line opacity
     real, pointer :: etaLine(:,:,:)             ! line emissivity
     real, pointer :: etaCont(:,:,:)             ! line emissivity
     real, pointer :: sigma(:,:,:)               ! radial velocity gradient
     type(VECTOR), pointer :: velocity(:,:,:)    ! velocity grid
     real, pointer :: rAxis(:)                   ! r-axis

     real, pointer :: biasCont(:)                ! radial bias for cont
     real, pointer :: biasLine(:)                ! radial bias for line

     real, pointer :: muAxis(:)                  ! mu-axis
     real, pointer :: phiAxis(:)                 ! phi-axis
     real          :: dMu, dPhi                  ! spacing of mu and phi axes
     real, pointer :: xAxis(:)                   ! x-axis
     real, pointer :: yAxis(:)                   ! y-axis
     real, pointer :: zAxis(:)                   ! z-axis
     real, pointer :: lamArray(:)                ! the wavelength array
     real, pointer :: rProbDistLine(:)               ! emissivity probability
     real, pointer :: muProbDistLine(:,:)            ! distributions
     real, pointer :: phiProbDistLine(:,:,:)         !
     real, pointer :: xProbDistLine(:)               ! emissivity probability
     real, pointer :: yProbDistLine(:,:)             ! distributions
     real, pointer :: zProbDistLine(:,:,:)           !
     real, pointer :: rProbDistCont(:)               ! emissivity probability
     real, pointer :: muProbDistCont(:,:)            ! distributions
     real, pointer :: phiProbDistCont(:,:,:)         !
     real, pointer :: xProbDistCont(:)               ! emissivity probability
     real, pointer :: yProbDistCont(:,:)             ! distributions
     real, pointer :: zProbDistCont(:,:,:)           !
     real, pointer :: temperature(:,:,:)             ! grid cell temperatures
     real, pointer :: biasLine3D(:,:,:)              ! grid bias distribution
     real, pointer :: biasCont3D(:,:,:)              ! grid bias distribution
     real :: rCore                                   ! core radius
     real(kind=doublekind) :: lCore                  ! core luminosity
     real :: chanceWindOverTotalContinuum            ! chance of continuum photon being produced
     ! in the wind rather than at the core

     logical :: lineEmission                        ! is there line emission?
     logical :: contEmission                        ! is there continuum emission?

     logical :: doRaman                             ! is this a raman-scattering model?
     logical :: resonanceLine

     real :: rStar1, rStar2, lumRatio

     real :: tempSource

     type(VECTOR) :: starPos1, starPos2

     real :: lambda2

     real(kind=doubleKind), pointer :: N(:,:,:,:)           ! stateq level pops
     real(kind=doubleKind), pointer :: Ne(:,:,:)            ! electron density
     real(kind=doubleKind), pointer :: nTot(:,:,:)          ! total density
     logical(kind=1), pointer :: inStar(:,:,:)
     logical(kind=1), pointer :: inUse(:,:,:)

     ! adaptive mesh refinement stuff
     integer :: maxLevels                            ! number of stateq levels 
     type(octal), POINTER :: octreeRoot => NULL()    ! root of adaptive mesh octree 
     integer :: maxDepth                             ! maximum depth of octree
     real(kind=octalKind) :: halfSmallestSubcell 
       ! 'halfSmallestSubcell' is half the vertex length of the grid's smallest
       !    subcells. we store this because it is used frequently in calculations.

  end type GRIDTYPE

end module gridtype_mod
