
! this module contains the description of the 3d grid of opacity arrays
! and contains subroutines to initialize the grid with a variety of 
! geometries.

! written by tjh

! v1.0 on 13/08/99

! v1.1 on 08/01/00
!   added colliding winds and bipolar geometries

! v1.2 on 19/05/00
! added disk stuff in LTE



module grid_mod


  use kind_mod
  use constants_mod                   ! physical constants
  use vector_mod                      ! vector math
  use atom_mod                        ! LTE atomic physics
  use utils_mod

  implicit none

  public

  ! this is the grid type

  type GRIDTYPE
     integer :: nx                               ! size of cartesian grid in x
     integer :: ny                               ! size of cartesian grid in y
     integer :: nz                               ! size of cartesian grid in z
     integer :: nr                               ! size of polar grid in r
     integer :: nmu                              ! size of polar grid in mu
     integer :: nphi                             ! size of polar grid in phi
     integer :: na1, na2, na3                    ! generic 3d grid sizes
     integer :: nlambda                          ! number of wavelength points
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
     real(kind=doublekind) :: lCore                                   ! core luminosity
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

     real(kind=doubleKind), pointer :: N(:,:,:,:)             ! stateq level pops
     real(kind=doubleKind), pointer :: Ne(:,:,:)            ! electron density
     real(kind=doubleKind), pointer :: nTot(:,:,:)          ! total density
     logical(kind=1), pointer :: inStar(:,:,:)
     logical(kind=1), pointer :: inUse(:,:,:)


  end type GRIDTYPE



contains

  ! function to initialize a cartesian grid

  function initCartesianGrid(nx, ny, nz, nLambda, lamStart, lamEnd, &
       flatspec, ok)

    type(GRIDTYPE) :: initCartesianGrid   ! the grid
    integer :: nx, ny, nz                 ! x,y,z sizes
    integer :: nLambda                    ! no of wavelength points
    integer :: ierr                       ! allocation error status
    integer :: i,ilambda                  ! counters
    real    :: lamStart, lamEnd           ! start and end wavelengths
    logical :: ok                         ! function done ok?
    logical :: flatspec                   ! is the spectrum flat

    iLambda = nLambda

    initCartesianGrid%lineEmission = .false. ! the default

    ! if the spectrum is flat one only needs on wavelength point

    if (flatspec) ilambda = 1

    ! ok for now

    ok = .true.

    ! allocate the grid arrays, and error out of there isn't
    ! enough memory

    allocate(initCartesianGrid%rho(1:nx, 1:ny, 1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initCartesianGrid%kappaAbs(1:nx, 1:ny, 1:nz, 1:iLambda),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initCartesianGrid%kappaSca(1:nx, 1:ny, 1:nz, 1:iLambda),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initCartesianGrid%velocity(1:nx, 1:ny, 1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initCartesianGrid%temperature(1:nx, 1:ny, 1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initCartesianGrid%chiLine(1:nx, 1:ny, 1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initCartesianGrid%etaLine(1:nx, 1:ny, 1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    allocate(initCartesianGrid%etaCont(1:nx, 1:ny, 1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif


    allocate(initCartesianGrid%sigma(1:nx, 1:ny, 1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initCartesianGrid%biasLine3D(1:nx, 1:ny, 1:nz), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initCartesianGrid%biasCont3D(1:nx, 1:ny, 1:nz), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    initCartesianGrid%biasLine3D = 1.
    initCartesianGrid%biasCont3D = 1.


    ! allocate axes - should really check error status of these, too

    allocate(initCartesianGrid%xAxis(1:nx))
    allocate(initCartesianGrid%yAxis(1:ny))
    allocate(initCartesianGrid%zAxis(1:nz))
    allocate(initCartesianGrid%lamArray(1:nlambda))
    initCartesianGrid%nx = nx
    initCartesianGrid%ny = ny
    initCartesianGrid%nz = nz

    initCartesianGrid%na1 = nx
    initCartesianGrid%na2 = ny
    initCartesianGrid%na3 = nz

    

    initCartesianGrid%nLambda = nLambda

    allocate(initCartesianGrid%xProbDistLine(1:nx),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    allocate(initCartesianGrid%yProbDistLine(1:nx,1:ny),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    allocate(initCartesianGrid%zProbDistLine(1:nx,1:ny,1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    allocate(initCartesianGrid%xProbDistCont(1:nx),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    allocate(initCartesianGrid%yProbDistCont(1:nx,1:ny),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    allocate(initCartesianGrid%zProbDistCont(1:nx,1:ny,1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    allocate(initCartesianGrid%inStar(1:nx,1:ny,1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    initCartesiangrid%inStar = .false.

    allocate(initCartesianGrid%inUse(1:nx,1:ny,1:nz),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    initCartesiangrid%inUse = .true.


    ! initialize the arrays with zeros

    initCartesianGrid%rho = 1.e-30
    initCartesianGrid%kappaAbs = 1.e-30
    initCartesianGrid%kappaSca = 1.e-30
    initCartesianGrid%cartesian = .true.

    initCartesianGrid%resonanceLine = .false.

    ! allocate a few arrays to prevent passing uninitialized pointers

    allocate(initCartesianGrid%rProbDistLine(1:1))
    allocate(initCartesianGrid%muProbDistLine(1:1,1:1))
    allocate(initCartesianGrid%phiProbDistLine(1:1,1:1,1:1))
    allocate(initCartesianGrid%rProbDistCont(1:1))
    allocate(initCartesianGrid%muProbDistCont(1:1,1:1))
    allocate(initCartesianGrid%phiProbDistCont(1:1,1:1,1:1))


    allocate(initCartesianGrid%rAxis(1:1))
    allocate(initCartesianGrid%muAxis(1:1))
    allocate(initCartesianGrid%phiAxis(1:1))

666 continue
  end function initCartesianGrid

  ! analogous subroutine for a spherical polar grid

  function initPolarGrid(nr, nmu, nphi, nlambda, lamStart, lamEnd, &
       flatspec, ok)
    type(GRIDTYPE) :: initPolarGrid
    integer :: nr, nmu, nphi
    integer :: nLambda
    integer :: ierr
    integer :: i, ilambda
    real    :: lamStart, lamEnd
    logical :: ok, flatspec

    initPolarGrid%lineEmission = .false. ! the default

    iLambda = nLambda

    if (flatspec) iLambda = 1

    ok = .true.
    allocate(initPolarGrid%rho(1:nr, 1:nmu, 1:nphi),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%kappaAbs(1:nr, 1:nmu, 1:nphi, 1:iLambda),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif


    allocate(initPolarGrid%kappaSca(1:nr, 1:nmu, 1:nphi, 1:iLambda),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%velocity(1:nr, 1:nmu, 1:nphi),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%temperature(1:nr, 1:nmu, 1:nphi),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%chiLine(1:nr, 1:nmu, 1:nphi),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%etaLine(1:nr, 1:nmu, 1:nphi),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%etaCont(1:nr, 1:nmu, 1:nphi),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%sigma(1:nr, 1:nmu, 1:nphi),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif


    allocate(initPolarGrid%rProbDistLine(1:nr),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%muProbDistLine(1:nr, 1:nMu), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%phiProbDistLine(1:nr, 1:nMu, 1:nPhi), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%rProbDistCont(1:nr),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%muProbDistCont(1:nr, 1:nMu), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%phiProbDistCont(1:nr, 1:nMu, 1:nPhi), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%inUse(1:nr,1:nMu,1:nPhi),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif
    initPolargrid%inUse = .true.


    allocate(initpolargrid%inStar(1:nr, 1:nmu, 1:nphi), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    initPolarGrid%inStar = .false.



    allocate(initPolarGrid%rAxis(1:nr))
    allocate(initPolarGrid%muAxis(1:nmu))
    allocate(initPolarGrid%phiAxis(1:nphi))
    allocate(initPolarGrid%lamArray(1:nlambda))

    allocate(initPolarGrid%biasLine(1:nr))
    allocate(initPolarGrid%biasCont(1:nr))

    initPolarGrid%biasLine = 1.
    initPolarGrid%biasCont = 1.


    initPolarGrid%nr = nr
    initPolarGrid%nmu = nmu
    initPolarGrid%nphi = nphi

    initPolarGrid%na1 = nr
    initPolarGrid%na2 = nmu
    initPolarGrid%na3 = nphi


    initPolarGrid%nlambda = nlambda

    initPolarGrid%rho = 0.
    initPolarGrid%kappaSca = 0.
    initPolarGrid%kappaAbs = 0.
    initPolarGrid%cartesian = .false.

    initPolarGrid%resonanceLine = .false.

    ! allocate a few arrays to prevent passing uninitialized pointers

    allocate(initPolarGrid%xProbDistLine(1:1))
    allocate(initPolarGrid%yProbDistLine(1:1,1:1))
    allocate(initPolarGrid%zProbDistLine(1:1,1:1,1:1))
    allocate(initPolarGrid%xProbDistCont(1:1))
    allocate(initPolarGrid%yProbDistCont(1:1,1:1))
    allocate(initPolarGrid%zProbDistCont(1:1,1:1,1:1))

    allocate(initPolarGrid%xAxis(1:1))
    allocate(initPolarGrid%yAxis(1:1))
    allocate(initPolarGrid%zAxis(1:1))

666 continue
  end function initPolarGrid

  ! this subroutine sets up a shell grid of constant density

  subroutine fillGridShell(grid, radius, shellFrac, rho, kfac)

    implicit none
    type(GRIDTYPE) :: grid                  ! the opacity grid
    real :: rho                             ! the density
    real :: radius                          ! the radius
    real :: shellFrac                       ! size of shell in terms of radius
    real :: kfac                            ! latitudinal density parameter
    real :: rMin, rMax                      ! inner and outer radii
    integer :: i                            ! counter

    ! this is a shell geometry

    grid%geometry = "shell"

    ! set up outer an inner radii

    rMin = radius*(1.-shellFrac)
    rMax = radius

    ! logarithmic radial grid

    do i = 1 , grid%nr
       grid%rAxis(i) = log(rMin) + (log(rMax)-log(Rmin)) * &
            real(i-1)/real(grid%nr-1)
    enddo
    grid%rAxis = exp(grid%rAxis(1:grid%nr))

    ! evenly spaced mu and phi axes

    do i = 1, grid%nmu
       grid%muAxis(i) = 2.*real(i-1)/real(grid%nmu-1)-1.
    enddo

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nphi-1)
    enddo

    ! density grid has latitudinal dependence

    do i = 1, grid%nMu
       grid%rho(1:grid%nr, i, 1:grid%nPhi) = rho * (1.-kfac*grid%muAxis(i)**2)
    enddo

  end subroutine fillGridShell


  ! fills grid with a latitudinal density grid

  subroutine fillGridSpheriod(thisGrid, rho,  radius, kfac)

    type(GRIDTYPE) :: thisGrid               ! opacity grid
    integer, parameter :: nRad = 100         ! number of radial points
    type(VECTOR) :: rVec                     ! radius vector
    real :: theta, phi, r, w                 ! spherical polar coords
    real :: rho                              ! density and scale
    real,intent(in) :: radius
    real :: kfac                             ! latitudinal density factor
    integer, parameter :: nTheta = 101, nPhi = 200
    integer :: i, j, k                       ! counters
    integer :: i1, i2, i3

    ! set up the grid type

    thisgrid%geometry = "spheriod"

    ! only for cartesians

    if (.not.thisGrid%cartesian) then
       write(*,'(a)') "Spheriod structure for cartesian grids only"
       stop
    endif

    ! axes run from -1 to +1

    do i = 1, thisGrid%nx
       thisGrid%xAxis(i) = 2.*real(i-1)/real(thisGrid%nx-1) - 1.
    enddo


    do i = 1, thisGrid%ny
       thisGrid%yAxis(i) = 2.*real(i-1)/real(thisGrid%ny-1) - 1.
    enddo


    do i = 1, thisGrid%nz
       thisGrid%zAxis(i) = 2.*real(i-1)/real(thisGrid%nz-1) - 1.
    enddo


    ! now loop over a spherical polar coord sphere

    do i = 1, nTheta
       theta = pi*real(i-1)/real(nTheta-1)
       do j = 1, nPhi
          phi = (2.*real(j-1)/real(nPhi-1)-1.)*pi
          do k = 1, nRad
             r = radius * real(k-1)/real(nRad-1)
             rVec%x = r*sin(theta)*cos(phi)
             rVec%y = r*sin(theta)*sin(phi)
             rVec%z = r*cos(theta)
             w = cos(theta)
             call locate(thisGrid%xAxis,thisGrid%nx,rVec%x, i1)
             call locate(thisGrid%yAxis,thisGrid%ny,rVec%y, i2)
             call locate(thisGrid%zAxis,thisGrid%nz,rVec%z, i3)

             ! set up density grid

             thisGrid%rho(i1,i2,i3) = rho*(1.-kfac*w*w)

          enddo
       enddo
    enddo ! loop over spherical polar sphere

  end subroutine fillGridSpheriod


  ! sets up an ellipsoidal of constant density - suitable for dusty
  ! envelopes

  subroutine fillGridEllipse(grid, rho, rMin, rMaj)


    implicit none
    integer, parameter :: nRmaj = 100
    type(GRIDTYPE) :: grid
    real :: x,y,z,r
    real :: phi
    real :: rho
    real :: rMaj, rMin
    integer, parameter :: nPhi = 200
    integer :: i, j, k
    integer :: i1, i2, i3, i4

    grid%geometry = "ellipse"

    grid%rCore = rSol / 1.e10
    grid%lCore = lSol
    grid%inUse = .false.

    grid%temperature = 100.

    ! only for cartesians

    if (.not.grid%cartesian) then
       write(*,'(a)') "Spheriod structure for cartesian grids only"
       stop
    endif

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1) - 1.
    enddo


    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1) - 1.
    enddo


    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1) - 1.
    enddo

    grid%xAxis = grid%xAxis * rMaj * 1.1
    grid%yAxis = grid%yAxis * rMaj * 1.1
    grid%zAxis = grid%zAxis * rMaj * 1.1


    ! use spherical polars to compute ellipsoid

    do i = 1, nRmaj

       r = rMaj*real(i-1)/real(nRmaj-1)

       z = sqrt(rMin*rMin*(max(0.,1.-(r/rMaj)**2)))

       ! find height of ellipse above x axis

       ! volume of revolution...

       do j = 1, nPhi

          phi = 2.*pi*real(j-1)/real(nPhi-1)
          x = r * cos(phi)
          y = r * sin(phi)

          call locate(grid%xAxis,grid%nx,x,i1)
          call locate(grid%yAxis,grid%ny,y,i2)
          call locate(grid%zAxis,grid%nz,z,i3)
          call locate(grid%zAxis,grid%nz,-z,i4)
          do k = min(i3,i4), max(i3,i4)
             Grid%rho(i1,i2,k) = rho    
             grid%inUse(i1,i2,k) = .true.
          enddo
       enddo
    enddo

  end subroutine fillGridEllipse


  subroutine fillGridDisk(grid, rhoNought, rCore, rInner, rOuter, height, &
       mCore, diskTemp)


    implicit none
    type(GRIDTYPE) :: grid
    integer, parameter :: nRad = 100
    real :: rOuter, rInner, rCore
    type(VECTOR) :: rVec, vVec
    real :: phi, r
    real :: rho, rhoNought
    real :: diskTemp
    real :: height, theta
    real :: vel, mCore
    type(VECTOR) :: spinAxis
    integer :: i, j, k, nMu
    real :: sinTheta

    grid%geometry = "disk"
    grid%lineEmission = .true.
    grid%rCore = rCore


    spinAxis = VECTOR(0.,0.,1.)

    grid%rho = 0.
    grid%etaLine = 1.e-30
    grid%etaCont = 1.e-30
    grid%chiLine = 1.e-30
    grid%kappaSca = 1.e-30
    grid%kappaAbs = 1.e-30
    grid%velocity = VECTOR(1.e-30,1.e-30,1.e-30)
    grid%temperature = diskTemp

    nMu = grid%nMu/2

    if (nMu == nint(real(grid%nMu)/2.)) then
       write(*,'(a)') "! nmu should be an odd number"
       stop
    endif

    theta = asin(3.*height/rOuter)

    write(*,'(a,f7.3,a)') "Using a finer latitude grid within ",theta*radtodeg," deg of the equator"

    call fillmuAxisDisk(grid%muAxis, grid%nMu, theta)

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi * real(i-1)/real(grid%nphi-1)
    enddo


    

    do i = 1, grid%nr
       grid%rAxis(i) = real(i - 1)/real(grid%nr - 1)
    enddo
    grid%rAxis = (10.**grid%rAxis - 1.)/9.
    grid%rAxis = grid%rAxis * (rOuter - rInner) + rInner

    call writeAxes(grid)

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             sinTheta = sqrt(1.-grid%muAxis(j)**2)
             rVec = VECTOR(grid%rAxis(i)*cos(grid%phiAxis(k))*sinTheta, &
                           grid%rAxis(i)*sin(grid%phiAxis(k))*sinTheta, &
                           grid%rAxis(i)*grid%muAxis(j))
             r = grid%rAxis(i)
             rho = rhoNought * (rInner/r)**2
             if (sqrt(rVec%x**2 + rVec%y**2) > rInner) then
                grid%rho(i,j,k) = rho * exp(-abs(rVec%z/height))
                rVec = rVec / r
                vVec = rVec  .cross. spinAxis
                if (modulus(vVec) /= 0.) then
                   call normalize(vVec)
                   vel = sqrt(bigG * mCore / r) / cSpeed
                   vVec = vel  * vVec 
                   grid%velocity(i,j,k) = vVec
                endif
             endif
          enddo
       enddo
    enddo

    
    do i = 1, grid%nr
       do j = 1, grid%nmu
          do k = 1, grid%nphi
             if (grid%rho(i,j,k) /= 0.) then
                grid%kappaSca(i,j,k,1) = max(1.e-30,grid%rho(i,j,k) * sigmaE)
             endif
          enddo
       enddo
    enddo

  end subroutine fillGridDisk


  subroutine fillGridDisk2(grid, rhoNought, rCore, rInner, rOuter, height, &
       mCore, diskTemp)


    implicit none
    type(GRIDTYPE) :: grid
    integer, parameter :: nRad = 1000
    real :: rOuter, rInner, rCore
    type(VECTOR) :: rVec, vVec
    real :: phi, r
    real :: rho, rhoNought
    real :: diskTemp
    real :: height, theta
    real :: vel, mCore
    type(VECTOR) :: spinAxis
    integer :: i, j, k, nMu
    real :: sinTheta

    grid%geometry = "disk"
    grid%lineEmission = .true.
    grid%rCore = rCore


    spinAxis = VECTOR(0.,0.,1.)

    grid%rho = 0.
    grid%etaLine = 1.e-30
    grid%etaCont = 1.e-30
    grid%chiLine = 1.e-30
    grid%kappaSca = 1.e-30
    grid%kappaAbs = 1.e-30
    grid%velocity = VECTOR(1.e-30,1.e-30,1.e-30)
    grid%temperature = diskTemp

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1)-1.
    enddo
    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1)-1.
    enddo
    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1)-1.
    enddo

    grid%xAxis(1:grid%nx) = grid%xAxis(1:grid%nx) * rOuter
    grid%yAxis(1:grid%ny) = grid%yAxis(1:grid%ny) * rOuter
    grid%zAxis(1:grid%nz) = grid%zAxis(1:grid%nz) * 5.*height
    do i = 1, grid%nx
       write(*,*) i ,grid%xAxis(i), grid%yaxis(i), grid%zaxis(i)
    enddo

    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))
             r = modulus(rVec)
             if ((r > rInner).and. (r < rOuter)) then
                rho = rhoNought * (rInner/r)**2
                grid%rho(i,j,k) = rho * exp(-abs(grid%zAxis(k)/height))
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))
                rVec = rVec / r
                vVec = rVec  .cross. spinAxis
                if (modulus(vVec) /= 0.) then
                   call normalize(vVec)
                   vel = sqrt(bigG * mCore / r) / cSpeed
                   vVec = vel  * vVec 
                   grid%velocity(i,j,k) = vVec
                endif
             endif
          enddo
       enddo
    enddo

    
    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             if (grid%rho(i,j,k) /= 0.) then
                grid%kappaSca(i,j,k,1) = max(1.e-30,grid%rho(i,j,k) * sigmaE)
             endif
          enddo
       enddo
    enddo

  end subroutine fillGridDisk2



  subroutine fillGridTorus(Grid, rho,  rTorus, rOuter)


    implicit none
    type(GRIDTYPE) :: Grid
    integer, parameter :: nRad = 50
    real :: phi
    real :: rTorus
    real :: rOuter
    real :: theta
    type(VECTOR) :: r0Vec, rUp, rDown, aVec, bVec, rHat
    type(VECTOR) :: torusAxis, normVec
    integer, parameter :: nAng = 360
    integer :: i, j, k
    integer :: i1, i2, i3, i4
    real :: rho

    grid%geometry = "torus"

    grid%rCore = rSol / 1.e10
    grid%lCore = lSol

    grid%rho = 1.e-20

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1) - 1.
    enddo


    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1) - 1.
    enddo


    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1) - 1.
    enddo

    torusAxis = VECTOR(0.,0.,1.)
    grid%inUse = .false.

    grid%xAxis = grid%xAxis * (rTorus + rOuter) * 1.1
    grid%yAxis = grid%yAxis * (rTorus + rOuter) * 1.1
    grid%zAxis = grid%zAxis * (rTorus + rOuter) * 1.1

    do i = 1, nAng

       theta  = real(i-1)/real(nAng-1)*twoPi
       r0Vec%x = rTorus * cos(theta)
       r0Vec%y = rTorus * sin(theta)
       r0Vec%z = 0.
       normVec = r0Vec .cross. torusAxis
       call normalize(normVec)

       rHat = r0Vec
       call normalize(rHat)

       do j = 1, nAng
	  phi  = real(j-1)/real(nAng-1)*Pi
          aVec = router * (cos(phi) * rHat)
	  bVec = router * (sin(phi) * torusAxis)
	  rUp = r0Vec + (aVec + bVec)
 	  rDown = r0Vec + (aVec - bVec)
	  call locate(grid%xAxis,grid%nx,rUp%x, i1)
	  call locate(grid%yAxis,grid%ny,rUp%y, i2)
	  call locate(grid%zAxis,grid%nz,rUp%z, i3)
	  call locate(grid%zAxis,grid%nz,rDown%z, i4)

          do k = min(i3,i4),max(i3,i4)
             Grid%rho(i1,i2,k) = rho
             grid%inUse(i1,i2,k) = .true.
          enddo
       enddo
    enddo

    grid%xAxis = grid%xAxis / 1.e10
    grid%yAxis = grid%yAxis / 1.e10
    grid%zAxis = grid%zAxis / 1.e10


  end subroutine fillGridTorus


  subroutine fillGridCollide(grid, rho, momRatio, binarySep, dust, meanDustParticleMass, logMassLossRate)

    type(GRIDTYPE) :: grid
    real :: rho, r
    real :: momRatio
    real :: binarySep
    real :: mDot
    real :: ang
    real :: logMassLossRate
    real, parameter :: v = 2000.e5
    real :: stagPoint, openingAngle, rotAng
    real :: minAng, maxAng
    integer :: i, j, k
    integer :: i1, i2, i3
    type(VECTOR) :: rVec
    real :: fac
    logical :: dust
    real :: meanDustParticleMass
    real(kind=doubleKind) :: vElement
    real(kind=doubleKind) :: totalDustMass
    logical, allocatable :: done(:,:,:)

    grid%geometry(1:7) = "collide"

    fac = 2.
    if (dust) fac=50.
    mdot = 10.**logMassLossRate * mSol / (365.25 * 24. * 60. * 60.)

    do i = 1, grid%nx
       grid%xAxis(i) = (fac*real(i-1)/real(grid%nx-1) - fac/2.)*binarySep
    enddo


    do i = 1, grid%ny
       grid%yAxis(i) = (fac*real(i-1)/real(grid%ny-1) - fac/2.)*binarySep
    enddo


    do i = 1, grid%nz
       grid%zAxis(i) = (fac*real(i-1)/real(grid%nz-1) - fac/2.)*binarySep
    enddo

    vElement = dble((grid%xAxis(2)-grid%xAxis(1)))**3

    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             r = sqrt(grid%xAxis(i)**2 + grid%yAxis(j)**2 + grid%zAxis(k)**2)
             rho = mdot / (4. * pi * r**2 * v)
             if (.not.dust) then
                grid%rho(i,j,k) = rho/mHydrogen
             endif
          enddo
       enddo
    enddo

    allocate(done(1:grid%nx,1:grid%ny,1:grid%nz))
    done = .false.

    ! two relationships from Eichler & Usov, 1993, ApJ, 402, 271

    stagPoint = binarySep * sqrt(momRatio) / (1. + sqrt(momRatio))

    openingAngle = 2.1*(1 - (momRatio**(2./5.))/4.)*momRatio**(1./3.)


    !    ymin = stagPoint*pi/2. * tan(openingAngle*0.5)
    !    ymax = stagPoint*pi/2. * tan(openingAngle*1.5)
    !    do k = 1, 10
    !       ang = openingAngle + openingAngle * (real(k-1)/9.-0.5)
    !       y = ymin + (ymax-ymin)*real(k-1)/9.
    !       rs = 2.*y/pi
    !       do j = 1 , 200
    !          rotAng = twoPi * real(j-1)/199.
    !          do i = 1, 100
    !             chi  = (pi/2.)  * real(i-1)/99.
    !             rCone = rs*chi/sin(chi)
    !             rVec%x = -rCone * cos(chi)
    !             rVec%y = rCone * sin(chi)
    !             rVec%z = 0.
    !             rVec = rotateX(rVec, rotAng)
    !             rVec%x = rVec%x + binarySep
    !             call locate(grid%xAxis,Grid%nx, rVec%x, i1)
    !             call locate(grid%yAxis,Grid%ny, rVec%y, i2)
    !             call locate(grid%zAxis,Grid%nz, rVec%z, i3)
    !             if (.not.done(i1,i2,i3)) then
    !                grid%rho(i1,i2,i3) = grid%rho(i1,i2,i3) * 4.
    !                done(i1,i2,i3) = .true.
    !             endif
    !          enddo
    !       enddo
    !    enddo

    minAng = openingAngle * 0.5
    maxAng = openingAngle * 1.5
    if (dust) then
       minAng = 0.
       maxAng = openingAngle
    endif

    totalDustMass = 0.d0


    write(*,*) "stagpoint/binarysep",stagpoint/binarysep
    do j = 1 , 200
       rotAng = twoPi * real(j-1)/199.
       do i = 1, 100
          do k = 1, 100
             ang = minAng + (maxAng-minAng)*(real(k-1)/99.)
             rVec%x = grid%xAxis(grid%nx) * real(i-1)/99.
             rVec%y = rVec%x * tan(ang)
             rVec%z = 0.
             rVec = rotateX(rVec, rotAng)
             rVec%x = rVec%x + binarySep  - stagPoint
             call locate(grid%xAxis,Grid%nx, rVec%x, i1)
             call locate(grid%yAxis,Grid%ny, rVec%y, i2)
             call locate(grid%zAxis,Grid%nz, rVec%z, i3)
             r = modulus(rVec)
             rho = mdot / (4. * pi * r**2 * v)
             if (.not.done(i1,i2,i3)) then
                if (dust) then

                   if ((r > 2.*binarySep).and.(r < 100.*binarySep)) then
                      rho = rho / (4.*mHydrogen)  ! helium number density
                      rho = rho * 0.2             ! carbon number density
                      rho = rho * (12.*mHydrogen) ! carbon mass density
                      rho = 0.1 * rho / meanDustParticleMass
                      totalDustMass = totalDustMass + &
                           dble(rho)*vElement*dble(meanDustParticleMass)
                   else
                      rho = 0.
                   endif
                endif
                if (r /= 0.) then
                   grid%rho(i1,i2,i3) = rho
                   done(i1,i2,i3) = .true.
                endif
             endif
          enddo
       enddo
    enddo
    if (dust) then
       write(*,'(a,1pe8.2)') "Total dust mass (solar masses): ",totalDustMass/mSol
    endif
    deallocate(done)

  end subroutine fillGridCollide

  subroutine fillGridBipolar(Grid, rho,  openingAng)


    implicit none
    type(GRIDTYPE) :: Grid
    real :: rho, openingAng, ang
    real :: rInner, rOuter
    real :: zMin, zMax, zDash
    real :: xDash, yDash, rDash
    integer :: i, j, k, i1, i2, i3
    integer, parameter :: nCircle = 100, nRad = 20, nZ = 20

    grid%geometry = "bipolar"

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1) - 1.
    enddo


    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1) - 1.
    enddo


    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1) - 1.
    enddo

    ! bipolar lobes

    write(*,'(a,f4.1)') "Filling bipolar with opening angle of ",openingAng
    do i = 1, grid%nz
       rInner = abs(sqrt(abs(grid%zAxis(i)))*sin(openingAng*degToRad))
       rOuter = rInner * 0.1
       do j=1,nCircle
          ang = twoPi * real(j-1)/real(nCircle-1)
          do k = 1, nRad
             rDash = rInner + (rOuter - rInner) * (real(k-1)/real(nRad-1))
             xDash = rDash * cos(ang)
             yDash = rDash * sin(ang)
             call locate(grid%xAxis,Grid%nx, xDash, i1)
             call locate(grid%yAxis,Grid%ny, yDash, i2)
             grid%rho(i1,i2,i) = rho
          enddo
       enddo
    enddo

    ! obsuring disk

    zMin = -0.1
    zMax = 0.1
    rInner = 0.2
    rOuter = 0.5

    do k = 1, nZ
       zDash = zMin + (zMax - zMin)*real(k-1)/real(nZ-1)
       call locate(grid%zAxis, grid%nz, zDash, i3)
       do i = 1, nCircle
          ang = twoPi * real(i-1)/real(nCircle-1)
          do j = 1, nRad
             rDash = rInner + (rOuter - rInner) * (real(j-1)/real(nRad-1))
             xDash = rDash * cos(ang)
             yDash = rDash * sin(ang)
             call locate(grid%xAxis,Grid%nx, xDash, i1)
             call locate(grid%yAxis,Grid%ny, yDash, i2)
             grid%rho(i1,i2,i3) = rho*1.e10
          enddo
       enddo
    enddo

  end subroutine fillGridBipolar


  subroutine fillGridDustBlob(grid, distance, rho, phi, xSize, ySize, zSize)

    type(GRIDTYPE) :: grid
    real :: distance, rho
    real :: phi
    real :: xSize, ySize, zSize
    real :: xCen, yCen, zCen
    real :: xDash, yDash,zDash
    integer :: i, j, k
    integer :: i1, i2, i3

    grid%geometry = "dustblob"
    if (.not.Grid%cartesian) then
       write(*,*) "! Need to use cartesian grid for dust blob geometry"
       stop
    endif

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1) - 1.
    enddo


    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1) - 1.
    enddo


    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1) - 1.
    enddo


    xCen = distance * cos(phi)
    yCen = distance * sin(phi)
    zCen = 0.

    write(*,*) "distance",distance
    do i = 1, 20
       do j = 1, 20
          do k= 1, 20
             xDash = xCen - xSize/2. + xSize * real(i-1)/19.
             yDash = yCen - ySize/2. + ySize * real(j-1)/19.
             zDash = zCen - zSize/2. + zSize * real(k-1)/19.
             call locate(grid%xAxis, grid%nx, xDash, i1)
             call locate(grid%yAxis, grid%ny, yDash, i2)
             call locate(grid%zAxis, grid%nz, zDash, i3)
             grid%rho(i1,i2,i3) = rho
          enddo
       enddo
    enddo

  end subroutine fillGridDustBlob





  subroutine fillGridStar(Grid, radius, mdot, vel, kfac, scale)

    real :: scale
    type(GRIDTYPE) :: Grid
    type(VECTOR) :: rHat
    real :: radius, mdot, rMin, rMax, vel, kfac, vr, u, v, w, t, fac
    real :: dV
    integer :: i, j, k

    grid%geometry = "star"


    if (Grid%cartesian) then
       write(*,*) "! Need to use polar grid for star geometry"
       stop
    endif

    rMin = 1.
    rMax = 10.


    do i = 1 , grid%nr
       grid%rAxis(i) = log(rMin) + (log(rMax)-log(Rmin))*real(i-1)/real(grid%nr-1)
    enddo

    grid%rAxis = exp(grid%rAxis(1:grid%nr))

    grid%rAxis = radius * grid%rAxis

    do i = 1, grid%nmu
       grid%muAxis(i) = 2.*real(i-1)/real(grid%nmu-1)-1.
    enddo

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nphi-1)
    enddo


    do i = 1, grid%nr
       do j = 1, grid%nmu
          do k = 1 , grid%nphi
             vr = 10.e5+(vel-10.e5)*(1.-grid%rAxis(1)/grid%rAxis(i))
             grid%rho(i,j,k) = log(mdot)-log(fourPi)-log(vr)- &
                  2.*log(grid%rAxis(i)*scale)-log(mHydrogen)+log(1.-kfac*grid%muAxis(j)**2)
             grid%rho(i,j,k) = exp(grid%rho(i,j,k))
             w = grid%muAxis(j)
             t = sqrt(1.-w*w)
             u = t * cos(grid%phiAxis(k))
             v = t * sin(grid%phiAxis(k))
             rHat = VECTOR(u,v,w)
             grid%velocity(i,j,k) = (vr / cspeed) * rHat

             if (i < grid%nr) then
                dV = grid%rAxis(i)**2 * (grid%rAxis(i+1)-grid%rAxis(i)) * &
                     sqrt(1.-grid%muAxis(j)**2)*(grid%muAxis(2)-grid%muAxis(1)) * &
                     (grid%phiAxis(2)-grid%phiAxis(1))
             else
                dV = grid%rAxis(i)**2 * (grid%rAxis(i)-grid%rAxis(i-1)) * &
                     sqrt(1.-grid%muAxis(j)**2)*(grid%muAxis(2)-grid%muAxis(1)) * &
                     (grid%phiAxis(2)-grid%phiAxis(1))
             endif
             dV = abs(dV)

             fac = sqrt(1.- (grid%rAxis(1)**2/(grid%rAxis(i)**2)))

             grid%etaLine(i,j,k) = ((grid%rho(i,j,k))**2)*fac

          enddo
       enddo
    end do



    grid%etaLine = grid%etaLine / SUM(grid%etaLine)



  end subroutine fillGridStar

  subroutine fillGridSpiral(Grid, radius, mdot, vel, scale)

    real :: scale
    type(GRIDTYPE) :: Grid
    real :: radius, mdot, rMin, rMax, vel, v,  r, x
    integer :: i, j, k
    integer :: nSpiral, iSpiral
    integer :: i1,i2
    real :: muStart,muEnd
    real :: thickness,  woundFac
    integer :: j1,j2
    real :: phi, theta, phaseOffset, spiralPhase
    type(VECTOR) :: rHat, perp, zAxis, startVec, endVec, thisVec

    grid%geometry = "spiral"

    zAxis=VECTOR(0.,0.,1.)

    if (Grid%cartesian) then
       write(*,*) "! Need to use polar grid for spiral geometry"
    endif

    rMin = radius
    rMax = 10.*radius

    do i = 1 , grid%nr
       grid%rAxis(i) = log(rMin) + (log(rMax)-log(Rmin))*real(i-1)/real(grid%nr-1)
    enddo

    grid%rAxis = exp(grid%rAxis(1:grid%nr))

    do i = 1, grid%nmu
       grid%muAxis(i) = 2.*real(i-1)/real(grid%nmu-1)-1.
    enddo

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nphi-1)
    enddo




    do i = 1, grid%nr
       do j = 1, grid%nmu
          do k = 1 , grid%nphi
             v = 10.e5+(vel-10.e5)*(1.- radius/grid%rAxis(i))**1
             grid%rho(i,j,k) = log(mdot)-log(fourPi)-log(v)-2.*log(grid%rAxis(i)*scale)-log(mHydrogen)
             grid%rho(i,j,k) = exp(grid%rho(i,j,k))
          enddo
       enddo
    end do

    muStart = 0.5
    muEnd = -0.5
    thickness = 0.2
    woundFac = 0.75
    phaseOffset = 0.25*twoPi
    nSpiral = 2.

    call locate(grid%muAxis,grid%nMu,muStart,j1)
    call locate(grid%muAxis,grid%nMu,muEnd,j2)

    do iSpiral = 1 , nSpiral
       spiralPhase = real(iSpiral-1)*twoPi/real(nSpiral)

       do i = 1, grid%nr
          r = grid%rAxis(i)
          phi = twoPi*woundFac*real(i-1)/real(grid%nr-1) + phaseOffSet + spiralPhase
          if (phi > twoPi) phi = phi - twoPi 

          rHat = VECTOR(cos(phi),sin(phi),0.)
          perp = rHat .cross. zAxis

          startVec = (r*rHat) - ((thickness/2.)*r)*perp
          endVec = (r*rHat) + ((thickness/2.)*r)*perp


          do i1 = 1,100
             x = real(i1-1)/99.
             thisVec = startVec + (x*(endVec-startVec))
             call getPolar(thisVec, r, theta, phi)
             call locate(grid%rAxis,grid%nr,r,i2)
             call locate(grid%phiAxis,grid%nPhi,phi,k)

             do j = min(j1,j2), max(j1,j2)
                v = 10.e5+(vel-10.e5)*(1.- radius/grid%rAxis(i2))**1
                grid%rho(i2,j,k) = log(mdot)-log(fourPi)-log(v)-2.*log(grid%rAxis(i2)*scale)-log(mHydrogen) + log(5.)
                grid%rho(i2,j,k) = exp(grid%rho(i2,j,k))
             enddo
          enddo
       enddo


    enddo
  end subroutine fillGridSpiral



  subroutine plotGrid(Grid, viewVec, opaqueCore, device, contTau, fg, bg, coolStarPosition, firstPlot)


    implicit none
    logical :: opaqueCore, firstPlot
    type(GRIDTYPE) :: grid
    real :: x, y
    real :: minTau, maxTau
    real :: junk
    real :: imageExtent
    logical :: hitcore
    real :: contTau(2000,*)
    real :: rMax
    type(VECTOR) :: viewVec, posVec
    type(VECTOR) :: xProj, zAxis
    type(VECTOR) :: yProj, intersection, coolStarPosition
    integer :: resFac = 2
    integer :: i, j, npix
    real, allocatable :: plane(:,:), axis(:)
    logical, allocatable :: inuse(:,:)
    real :: tr(6)
    real :: escprob
    real :: fg, bg, dx, dz
    integer :: i1, i2, i3
    real :: r, mu, phi
    real :: modz
    integer :: nTau
    character(len=*) :: device
    logical :: ok

    if (grid%cartesian) then

       allocate(plane(1:grid%nx,1:grid%nz))

       plane = grid%rho(1:grid%nx,grid%ny/2,1:grid%nz)
       
       dx = grid%xAxis(2) - grid%xAxis(1)
       dz = grid%zAxis(2) - grid%zAxis(1)
       tr(1) =  - dx + grid%xAxis(1)
       tr(2) = dx
       tr(3) = 0.
       tr(4) =  - dz + grid%zAxis(1)
       tr(5) = 0.
       tr(6) = dz
       
       plane = log10(plane)
       fg = MAXVAL(plane, mask=grid%inUse(1:grid%nx,grid%ny/2,1:grid%nz))
       bg = MINVAL(plane, mask=grid%inUse(1:grid%nx,grid%ny/2,1:grid%nz))
       call pgbegin(0,device,1,1)
       
       call pgvport(0.1,0.9,0.1,0.9)
       
       call pgwnad(grid%xAxis(1), grid%xAxis(grid%nx), grid%zAxis(1), grid%zAxis(grid%nz))
       
       call pgimag(plane, grid%nx, grid%nz, 1, grid%nx, 1, grid%nz, fg, bg, tr)
       
       call pgbox('bcnst',0,0,'bcnst',0,0)
       call pgend
       
    else

       npix = 4*grid%nr
       allocate(plane(1:npix,1:npix))
       allocate(inuse(1:npix,1:npix))
       allocate(axis(1:npix))
       plane = 1.e-30
       inuse = .false.
       do i = 1, npix
          axis(i) = grid%rAxis(grid%nr)*(2.*real(i-1)/real(npix-1)-1.)
       enddo

       do i = 1, npix
          do j = 1, nPix
             r = sqrt(axis(i)**2 + axis(j)**2)
             mu = axis(j) / r
             phi = atan2(0.,axis(i))
             if (phi < 0.) phi = phi + twoPi
             call locate(grid%rAxis,grid%nr, r, i1)
             call locate(grid%muAxis,grid%nmu, mu, i2)
             call locate(grid%phiAxis,grid%nphi, phi, i3)
             if (grid%inUse(i1,i2,i3)) then
                plane(i,j) = grid%rho(i1,i2,i3)
                inUse(i,j) = .true.
             endif
          enddo
       enddo

       
       dx = axis(2) - axis(1)
       tr(1) =  - dx + axis(1)
       tr(2) = dx
       tr(3) = 0.
       tr(4) =  - dx + axis(1)
       tr(5) = 0.
       tr(6) = dx
       
       plane = log10(plane)
       fg = MAXVAL(plane, mask=inuse)
       bg = MINVAL(plane, mask=inuse)
       write(*,*) "fg, bg",fg,bg
       call pgbegin(0,device,1,1)
       
       call pgvport(0.1,0.9,0.1,0.9)
       
       call pgwnad(axis(1), axis(npix), axis(1), axis(npix))
       
       call pgimag(plane, npix, npix, 1, npix, 1, npix, fg, bg, tr)
       
       call pgbox('bcnst',0,0,'bcnst',0,0)
       call pgend




       write(*,'(a)') "Finished plotting grid..."

    endif


  end subroutine plotGrid


  type(GRIDTYPE) function plezModel(filename, lamStart, lamEnd, nLambda, kfac)



    character(len=*) filename
    real :: junk1, junk2, junk3, junk4
    integer :: npts, i, j, k
    integer :: nLambda
    real :: lamStart, lamEnd
    integer,parameter :: nMu = 20, nPhi = 2
    real :: kfac
    logical :: ok
    real :: r, temp, rho, scale

    plezmodel%geometry = "plez"

    open(20, file=filename,status="old",form="formatted")
    npts = 1
    do
       read(20,*,end=10) junk1
       npts = npts + 1
    enddo
10  continue
    rewind(20)
    npts = npts - 1
    plezModel = initPolarGrid(npts, nmu, nphi, nLambda, lamStart, lamEnd, .false., ok)

    do i = npts, 1, -1
       read(20, *) junk1, junk2, r, temp, junk3, junk4, rho
       plezModel%rAxis(i) = r
       do j = 1, nMu 
          do k = 1, nPhi
             plezModel%rho(i,j,k) = 10.**rho
          enddo
       enddo
    enddo

    close(20)


    scale = plezModel%rAxis(1)
    plezModel%rAxis = plezModel%rAxis / scale

    do i = 1, plezModel%nmu
       plezModel%muAxis(i) = 2.*real(i-1)/real(plezModel%nmu-1)-1.
    enddo

    do i = 1, plezModel%nPhi
       plezModel%phiAxis(i) = twoPi*real(i-1)/real(plezModel%nphi-1)
    enddo

    do i = 1, nLambda
       plezModel%lamArray(i) = lamStart + &
            (lamEnd-lamStart)*real(i-1)/real(nLambda-1)
    enddo


    call fillGridTio(plezModel, kfac, scale)
    call fillGridRayleigh(plezModel, kfac, scale)  


  end function plezModel


  subroutine writeGrid(grid)

    type(GRIDTYPE) :: grid
    integer :: i, j, k

    do i = 1 , grid%nr
       do j = 1 , grid%nmu
          do k = 1 , grid%nphi

             write(*,*) i,j,k,grid%rAxis(i),grid%muAxis(j), grid%phiAxis(k), grid%rho(i,j,k)
          enddo
       enddo
    enddo

  end subroutine writeGrid


  subroutine fillGridRaman(grid, size, mdot, vinf, rcool, coolStarPosition, beta)

    type(GRIDTYPE) :: grid
    real :: size
    real :: mdot, vinf
    real :: rCool
    real :: beta, mu,fac
    integer :: i,j,k
    real :: rho, r, v
    integer :: iErr
    type(VECTOR) :: rVec,  rHat, coolStarPosition
    real :: sigmaAbsBlue, sigmaScaBlue, sigmaAbsRed, sigmaScaRed

    grid%geometry = "raman"
    grid%lineEmission = .false.

    grid%rCore = rCool
    sigmaScaBlue = (34.+6.6)*sigmaE
    sigmaAbsBlue = 1.e-5 * sigmaE

    sigmaScaRed = 1.e-5 * sigmaE
    sigmaAbsRed = 1.e-5 * sigmaE

    write(*,'(a)') "Filling raman grid..."
    allocate(grid%kappaAbsRed(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       stop
    endif
    allocate(grid%kappaScaRed(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       stop
    endif

    do i = 1, grid%nx
       grid%xAxis(i) = size*(2.*real(i-1)/real(grid%nx-1) - 1.)
    enddo


    do i = 1, grid%ny
       grid%yAxis(i) = size*(2.*real(i-1)/real(grid%ny-1) - 1.)
    enddo


    do i = 1, grid%nz
       grid%zAxis(i) = size*(2.*real(i-1)/real(grid%nz-1) - 1.)
    enddo


    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             rVec = VECTOR(grid%xAxis(i), grid%yAxis(j), grid%zAxis(k))-coolStarPosition
             rHat = rVec
             call normalize(rHat)
             r = modulus(rVec)
             mu = (grid%zAxis(k)-coolStarPosition%z) / r
             fac = 1.
             !             fac = abs(mu)+1.
             if (r > rCool) then
                v =  max(1.e5,vinf*(1.-rCool/r)**beta)
                grid%velocity(i,j,k) = (v/cSpeed) * rHat
                !                grid%velocity(i,j,k)%z = grid%velocity(i,j,k)%z * fac**2
                rho = mdot / (4.* pi * r**2 * v)
                rho = rho / mHydrogen                
                grid%kappaAbs(i,j,k,1:1) = sigmaAbsBlue * rho
                grid%kappaSca(i,j,k,1:1) = sigmaScaBlue * rho
                grid%kappaAbsRed(i,j,k,1:1) = sigmaAbsRed * rho
                grid%kappaScaRed(i,j,k,1:1) = sigmaScaRed * rho
             else
                !                rho = mdot / (4.* pi * rCool**2 * vinf)
                !                rho = rho / mHydrogen
                !                rho = rho  * exp (rCool/min(r,0.1*rCool)-1.)
                rho = 1.e12
                grid%kappaAbs(i,j,k,1:1) = sigmaAbsBlue * rho
                grid%kappaSca(i,j,k,1:1) = sigmaScaBlue * rho
                grid%kappaAbsRed(i,j,k,1:1) = sigmaAbsRed * rho
                grid%kappaScaRed(i,j,k,1:1) = sigmaScaRed * rho
                grid%velocity(i,j,k) = VECTOR(0.,0.,0.)
             endif
          enddo
       enddo
    enddo
    grid%kappaAbs(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1) = 1.e-30
    grid%kappaAbsRed(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1) = 1.e-30
    grid%kappaScaRed(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1) = 1.e-30
  end subroutine fillGridRaman



  subroutine fillGridStateq(grid, filename, xfac, scaleDensity)

    character(len=*) :: filename
    type(GRIDTYPE) :: grid
    integer :: nr, i, j, k
    real :: scaleDensity
    real :: sinTheta,tmp
    real, allocatable :: r(:), sigma(:), t(:), v(:), eta(:), chi(:)
    real, allocatable :: etal(:), chil(:), escat(:)
    real :: kfac, pfac, xfac,fac
    type(VECTOR) :: rHat

    scaleDensity = 1.e0

    grid%etaLine = 1.e-20
    grid%etaCont = 1.e-30
    grid%chiLine = 1.e-20
    grid%kappaSca = 1.e-20
    grid%kappaAbs = 1.e-20

    grid%geometry = "stateq"

    grid%lineEmission = .true.

    open(20,file=filename, status="old", form="formatted")
    read(20,*) nr

    allocate(r(1:nr))
    allocate(sigma(1:nr))
    allocate(t(1:nr))
    allocate(v(1:nr))
    allocate(eta(1:nr))
    allocate(etal(1:nr))
    allocate(chi(1:nr))
    allocate(chil(1:nr))
    allocate(escat(1:nr))

    pfac = 2.e0
    kfac = 1.0e0/(1.0e0-xFac/(1.0e0+pFac))

    if (grid%nr /= nr) then
       write(*,*) "! warning - wrong nr value"
    endif

    write(*,*) "Reading opacities from: ",trim(filename)
    do j = 1, nr
       i = nr - j + 1
       read(20,*) r(i), t(i), sigma(i), v(i), eta(i), chi(i), escat(i), &
            etal(i), chil(i)
    enddo
    close(20)


    do i = 1, grid%nmu
       grid%muAxis(i) = 2.e0*real(i-1)/real(grid%nmu-1)-1.
    enddo

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nphi-1)
    enddo

    grid%isotropic = .true.
    grid%dMu = grid%muAxis(2) - grid%muAxis(1)
    grid%dPhi = grid%phiAxis(2) - grid%phiAxis(1)



    grid%rAxis = r

!        do i = 1, nr
!           grid%rAxis(i*2-1) = r(i)
!        enddo
!        do i = 1, nr-1
!           grid%rAxis(i*2) = 0.5*(grid%rAxis(i*2-1)+grid%rAxis(i*2+1))
!        enddo


    do i = 1, grid%nr
       do j = 1, grid%nMu
          sinTheta = sqrt(1.e0 - grid%muAxis(j)**2)
          fac = kFac * (1.e0-xfac*abs(grid%muAxis(j))**pFac)
          do k = 1, grid%nPhi
             rHat = VECTOR(cos(grid%phiAxis(k))*sinTheta, &
                  sin(grid%phiAxis(k))*sinTheta, &
                  grid%muAxis(j))
             grid%velocity(i,j,k) = (1.e5*logInterp(v,nr,r,grid%rAxis(i))/cSpeed) * rHat
             grid%temperature(i,j,k) = logInterp(t, nr, r, grid%rAxis(i))*1.e4
             grid%chiLine(i,j,k) = (logInterp(chil, nr, r, grid%rAxis(i)))*fac**2
             grid%etaLine(i,j,k) = (logInterp(etal, nr,r, grid%rAxis(i)))*fac**2
             grid%sigma(i, j, k) = sigma(i)
             grid%kappaSca(i,j,k,1:grid%nLambda) = logInterp(escat, nr, r, grid%rAxis(i))*fac
             tmp = logInterp(chi, nr, r,grid%rAxis(i))-grid%kappaSca(i,j,k,1)
             grid%kappaAbs(i,j,k,1) = max(1.e-20,tmp)
          enddo
       enddo
    enddo

    write(*,*) "Finished reading opacities and filling grid"

    write(*,*) "Inner radius (Rsol): ",grid%rAxis(1)*1.e10/rSol
    write(*,*) "Outer radius (Rsol): ",grid%rAxis(grid%nr)*1.e10/rsol

    grid%rCore = grid%rAxis(1)

    deallocate(r)
    deallocate(sigma)
    deallocate(t)
    deallocate(v)
    deallocate(eta)
    deallocate(etal)
    deallocate(chi)
    deallocate(chil)
    deallocate(escat)



  end subroutine fillGridStateq



  subroutine plotVelocityVectors(grid, device)

    character(len=*) :: device
    type(GRIDTYPE) :: grid
    integer :: i,j,i1,i2,i3
    integer :: nx = 25, ny = 25
    real, allocatable :: xArray(:,:), yArray(:,:)
    real :: x,y,tr(6),r,phi
    real :: rMax,rMin,fac
    logical, allocatable :: usedArray(:,:)


    call pgbegin(0,device,1,1)


    allocate(xArray(1:nx,1:ny))
    allocate(yArray(1:nx,1:ny))
    allocate(usedArray(1:nx,1:ny))

    usedArray = .false.
    xArray = 0.
    yArray = 0.

    if (grid%cartesian) then
       tr(1) = grid%xAxis(1)
       tr(2) = (grid%xAxis(grid%nx)-grid%xAxis(1)) / real(nx)
       tr(3) = 0.
       tr(4) = grid%yAxis(1)
       tr(5) = 0.
       tr(6) = (grid%yAxis(grid%ny)-grid%yAxis(1)) / real(ny)
       do i = 1, nx
          do j = 1, ny
             x = tr(1) + tr(2)*i
             y = tr(1) + tr(2)*j
             call locate(grid%xAxis,grid%nx,x,i1)
             call locate(grid%zAxis,grid%nz,y,i3)
             xArray(i,j) = grid%velocity(i1,grid%ny/2,i3)%x*cSpeed
             yArray(i,j) = grid%velocity(i1,grid%ny/2,i3)%z*cSpeed
          enddo
       enddo


    else

       rMin = grid%rAxis(1)
       rMax = 2.*grid%rAxis(1)
       do i = 1, nx
          x = -rMax + 2.*rMax*real(i-1)/real(nx-1)
          do j = 1, ny
             y = -rMax + 2.*rMax*real(j-1)/real(ny-1)

             r = sqrt(x**2 + y**2)
             phi = atan2(y,x)
             if (phi < 0.) phi = phi + twoPi
             if (r > rMin) then
                call locate(grid%rAxis,grid%nr,r,i1)
                call locate(grid%phiAxis,grid%nphi,phi,i3)
                i2 = grid%nMu/2

                xArray(i,j) = grid%velocity(i1,i2,i3)%x
                yArray(i,j) = grid%velocity(i1,i2,i3)%y
                usedArray(i,j) = .true.
             endif

          enddo
       enddo

       tr(1) = -rMax
       tr(2) = 2.*rMax / real(nx)
       tr(3) = 0.
       tr(4) = -rMax
       tr(5) = 0.
       tr(6) = 2.*rMax / real(ny)
    endif



    fac = 0.
    call pgenv(tr(1),tr(1)+real(nx)*tr(2),tr(4),tr(4)+real(ny)*tr(6),1,0)
    call pgsch(0.5)
    call pgvect(xArray, yArray, nx, ny, 1, nx, 1, ny, fac, 0, tr, 0.)

    deallocate(xArray)
    deallocate(yArray)


    call pgend

  end subroutine plotVelocityVectors


  subroutine freeGrid(grid)

    type(GRIDTYPE) :: grid

    if (associated(grid%rho)) deallocate(grid%rho)
    if (associated(grid%kappaAbs)) deallocate(grid%kappaAbs)
    if (associated(grid%kappaSca)) deallocate(grid%kappaSca)
    if (associated(grid%chiLine)) deallocate(grid%chiLine)
    if (associated(grid%etaLine)) deallocate(grid%etaLine)
    if (associated(grid%sigma)) deallocate(grid%sigma)
    if (associated(grid%velocity)) deallocate(grid%velocity)
    if (associated(grid%rAxis)) deallocate(grid%rAxis)
    if (associated(grid%muAxis)) deallocate(grid%muAxis)
    if (associated(grid%phiAxis)) deallocate(grid%phiAxis)
    if (associated(grid%xAxis)) deallocate(grid%xAxis)
    if (associated(grid%yAxis)) deallocate(grid%yAxis)
    if (associated(grid%zAxis)) deallocate(grid%zAxis)
    if (associated(grid%lamArray)) deallocate(grid%lamArray)
    if (associated(grid%rProbDistLine)) deallocate(grid%rProbDistLine)
    if (associated(grid%muProbDistLine)) deallocate(grid%muProbDistLine)
    if (associated(grid%phiProbDistLine)) deallocate(grid%phiProbDistLine)
    if (associated(grid%rProbDistCont)) deallocate(grid%rProbDistCont)
    if (associated(grid%muProbDistCont)) deallocate(grid%muProbDistCont)
    if (associated(grid%phiProbDistCont)) deallocate(grid%phiProbDistCont)
    if (associated(grid%xProbDistLine)) deallocate(grid%xProbDistLine)
    if (associated(grid%yProbDistLine)) deallocate(grid%yProbDistLine)
    if (associated(grid%zProbDistLine)) deallocate(grid%zProbDistLine)
    if (associated(grid%xProbDistCont)) deallocate(grid%xProbDistCont)
    if (associated(grid%yProbDistCont)) deallocate(grid%yProbDistCont)
    if (associated(grid%zProbDistCont)) deallocate(grid%zProbDistCont)

  end subroutine freeGrid

  logical function insideGrid(grid, rVec)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec

    insideGrid = .false.

    if (grid%cartesian) then

       if ( (rVec%x > grid%xAxis(1)) .and. &
            (rVec%x < grid%xAxis(grid%nx)) .and. &
            (rVec%y > grid%yAxis(1)) .and. &
            (rVec%y < grid%yAxis(grid%ny)) .and. &
            (rVec%z > grid%zAxis(1)) .and. &
            (rVec%z < grid%zAxis(grid%nz)) ) then
          insideGrid = .true.
       endif

    else

       if (modulus(rVec) < grid%rAxis(grid%nr)) then
          insideGrid = .true.
       endif
    endif
  end function insideGrid



  subroutine fillGridBinary(grid, radius1, radius2, mass1, mass2, &
       temp1, temp2, vNought1, vNought2, vTerm1, vTerm2, beta1, beta2, period,&
       mdot1, mdot2, deflectionAngle, shockWidth, shockFac)


    type(GRIDTYPE) :: grid
    integer ::  i, j, k
    real :: rho
    real :: mdot1, mdot2

    real :: radius1, radius2
    real :: temp1, temp2
    real :: vNought1, vNought2
    real :: vTerm1, vTerm2
    real :: beta1, beta2


    real :: massRatio

    real :: middle


    type(VECTOR) :: rVecFromStar1, rVecFromStar2, position1, position2
    type(VECTOR) :: rHat1, rHat2, rVec, vVec, zAxis
    type(VECTOR) :: depthVec

    real :: shockWidth, shockFac

    real :: vel

    real :: r1, r2

    real :: period
    real :: gridSize

    real :: mass1, mass2

    real :: stagPoint, curlyR, momRatio
    real :: d1, d2

    integer :: i1,i2,i3

    real :: rotationAngle

    real :: dybydx,dx

    integer, parameter :: nSteps = 1000
    real :: xDist(nSteps), yDist(nSteps)
    type(VECTOR) :: direction(nSteps)

    real :: binarySep

    real :: t1, t2, t3

    real :: phi, mu, sinTheta

    real :: fac, xTemp

    real :: deflectionAngle 

    type(VECTOR) :: stagVec

    logical, allocatable :: shockCone(:,:,:)
    logical, allocatable :: leftOfCone(:,:,:)
    type(VECTOR), allocatable :: shockDirection(:,:,:)
    type(VECTOR) ::  thisDirection, thisVec

    real :: gridStep

    logical :: foundCone


    grid%inStar = .false.

    zAxis=VECTOR(0., 0., 1.)

    ! keplers third law - mass in solar masses, period in years
    ! gives binary separation in AU

    binarySep = (((mass1 + mass2)/mSol) * (period*secsToyears)**2)**0.333333

    binarySep = binarySep * auToCm / 1.e10 ! binary sep in 10^10cm

    shockWidth = shockWidth * binarySep

    massRatio = mass1/mass2

    d1 = binarySep * (1./(massRatio+1.))
    d2 = binarySep - d1


    position1 = VECTOR(-d1,0.,0.)
    position2 = VECTOR(d2,0.,0.)

    grid%starPos1 = position1
    grid%starPos2 = position2

    grid%velocity = VECTOR(0.,0.,0.)
    grid%etaLine = 1.e-30
    grid%etaCont = 0.
    grid%chiLine = 1.e-30
    grid%kappaSca = 1.e-30
    grid%kappaAbs = 1.e-30

    grid%geometry = "binary"

    grid%rCore = 0.

    grid%lineEmission = .true.

    if (.not.grid%cartesian) then
       write(*,'(a)') "! Binary geometry only in cartesian frame"
       stop
    endif

    gridSize = modulus(position1 - position2)  * 4.
    do i = 1, grid%nx
       grid%xAxis(i) = real(i-1)/real(grid%nx-1) * gridSize - gridSize / 2.
    enddo
    do i = 1, grid%ny
       grid%yAxis(i) = real(i-1)/real(grid%ny-1) * gridSize - gridSize / 2.
    enddo
    do i = 1, grid%nz
       grid%zAxis(i) = real(i-1)/real(grid%nz-1) * gridSize - gridSize / 2.
    enddo

    gridStep = grid%xAxis(2) - grid%xAxis(1)

    gridStep = 0.

    write(*,'(a)') "Filling grid with a binary..."

    grid%rStar1 = radius1 / 1.e10
    grid%rStar2 = radius2 / 1.e10

    grid%temperature = temp1

    momRatio = (mdot1 * vterm1) / (mdot2 * vterm2)

    write(*,'(a,f8.1)') "Wind momentum ratio: ", momRatio

    curlyR = sqrt(momRatio)         ! = d1/d2 (Equ 1.) 
    ! Stevens, Blondin & Pollock 1992
    ! ApJ 386 265

    stagPoint = binarySep * sqrt(momRatio) / (1. + sqrt(momRatio))

    stagVec = position1 + VECTOR(stagPoint,0.,0.)

    allocate(shockCone(1:grid%nx, 1:grid%ny, 1:grid%nz))
    allocate(leftOfCone(1:grid%nx, 1:grid%ny, 1:grid%nz))
    allocate(shockDirection(1:grid%nx, 1:grid%ny, 1:grid%nz))

    shockDirection = VECTOR(0.,0.,0.)
    shockCone(1:grid%nx,1:grid%ny,1:grid%nz) = .false.

    dx = (grid%xAxis(grid%nx)-stagVec%x)/real(nsteps)
    yDist(1) = dx*10.
    xDist(1) = stagVec%x 
    direction(1) = VECTOR(0.,1.,0.)
    do j = 2,nsteps
       xDist(j) = stagVec%x + real(j-1)*(grid%xAxis(grid%nx)-stagVec%x)/real(nSteps-1)
       d1 = sqrt((xDist(j-1)-position1%x)**2 + yDist(j-1)**2)
       d2 = sqrt((xDist(j-1)-position2%x)**2 + yDist(j-1)**2)

       dybydx = ((curlyR*d2**2 + d1**2)*yDist(j-1)) / &
            (curlyR *d2**2*(xDist(j-1)-position1%x) + d1**2*(xDist(j-1)-position2%x))

       yDist(j) = yDist(j-1) + dx * dyBydx

       direction(j) = VECTOR(dx,dybydx,0.)
       call normalize(direction(j))
    enddo

    do k = 1, 10
       depthVec = VECTOR(shockWidth*(real(k-1)/9.-0.5),0.,0.)
       do i = 1, 300
          rotationAngle = twoPi * real(i-1)/299.
          do j = 1, 10000
             xTemp = (xDist(2) - xDist(1))*real(j-1)/9999. + xDist(1)
             rVec = VECTOR(xTemp, (yDist(2)-yDist(1))*real(j-1)/9999.+yDist(1), 0.)
             thisVec = rotateX(rVec, rotationAngle)
             thisVec = thisVec - grid%starPos2
             thisVec = rotateZ(thisVec, -deflectionAngle)
             thisDirection = rotateX(direction(2), rotationAngle)
             thisVec = thisVec + grid%starPos2 + depthVec
             thisDirection = rotateZ(thisDirection, deflectionAngle)
             if (.not.outSideGrid(thisVec, grid)) then
                call getIndices(grid, thisVec, i1,i2,i3,t1,t2,t3)
                shockCone(i1,i2,i3) = .true.
                shockDirection(i1,i2,i3) = thisdirection
             endif
          enddo
       enddo
    enddo


    do k = 1, 10
       depthVec = VECTOR(shockWidth*(real(k-1)/9.-0.5),0.,0.)
       do i = 1, 300
          rotationAngle = twoPi * real(i-1)/299.
          do j = 2, nSteps
             rVec = VECTOR(xDist(j), yDist(j), 0.)
             thisVec = rotateX(rVec, rotationAngle)
             thisVec = thisVec - grid%starPos2
             thisVec = rotateZ(thisVec, -deflectionAngle)
             thisDirection = rotateX(direction(j), rotationAngle)
             thisVec = thisVec + grid%starPos2 + depthVec
             thisDirection = rotateZ(thisDirection, deflectionAngle)
             if (.not.outSideGrid(thisVec, grid)) then
                call getIndices(grid, thisVec, i1,i2,i3,t1,t2,t3)
                shockCone(i1,i2,i3) = .true.
                shockDirection(i1,i2,i3) = thisdirection
             endif
          enddo
       enddo
    enddo


    leftOfCone(1:grid%nx,1:grid%ny,1:grid%nz) = .false.
    do j = 1, grid%ny
       do k = 1, grid%nz
          foundCone = .false.
          i = 1
          do while((.not.shockCone(i,j,k)).and.(i<=grid%nx))
             leftOfCone(i,j,k) = .true.
             i = i + 1
          enddo
       enddo
    enddo



    write(*,'(a,f5.3)') "Primary radius/binary sep: ",grid%rStar1/binarySep
    write(*,'(a,f5.3)') "Secondary radius/binary sep: ",grid%rStar2/binarySep

    ! now we have to map the opacities onto the cartesian grid

    middle = 0.5*modulus(position2 + position1)

    do j = 1, grid%ny
       do k = 1, grid%nz
          do i = 1, grid%nx

             rVec = VECTOR(grid%xAxis(i), grid%yAxis(j), grid%zAxis(k))

             rVecFromStar1 =  rVec- &
                  position1

             rVecFromStar2 = rVec - &
                  position2

             ! need to fill in the opaque cores


             r1 = modulus(rVecFromStar1)
             r2 = modulus(rVecFromStar2)

             if (r1 < grid%rStar1) grid%inStar(i,j,k) = .true.
             if (r2 < grid%rStar2) grid%inStar(i,j,k) = .true.


             rHat1 = rVecFromStar1 / modulus(rVecFromStar1)
             rHat2 = rVecFromStar2 / modulus(rVecFromStar2)


             if (leftOfCone(i,j,k)) then
                if (r1 > grid%rStar1) then
                   vel = vNought1 + (vTerm1 - vNought1)*(1. - grid%rStar1/r1)**beta1
                   grid%velocity(i,j,k) = (vel / cSpeed) * rHat1
                   grid%temperature(i,j,k) = temp1
                   rho = mdot1 / (fourPi * (r1*1.e10)**2 * vel)
                   grid%rho(i,j,k) = rho
                endif
             else
                if (r2 > grid%rStar2) then
                   vel = vNought2 + (vTerm2 - vNought2)*(1. - grid%rStar2/r2)**beta2
                   grid%velocity(i,j,k) = (vel / cSpeed) * rHat2
                   grid%temperature(i,j,k) = temp2
                   rho = mdot2 / (fourPi * (r2*1.e10)**2 * vel)
                   grid%rho(i,j,k) = rho
                endif
             endif
          enddo
       enddo
    enddo

    do i = 1, 300
       phi = twoPi * real(i-1)/299.
       do j = 1, 200
          mu = -1. + 2.*real(j-1)/199.
          sinTheta = sqrt(1.-mu*mu)
          r1 = grid%rStar1
          rVec = VECTOR(r1*sinTheta*cos(phi),r1*sinTheta*sin(phi),r1*mu)+position1
          call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
          grid%velocity(i1,i2,i3) = VECTOR(0.,0.,0.)
          r1 = grid%rStar2
          rVec = VECTOR(r1*sinTheta*cos(phi),r1*sinTheta*sin(phi),r1*mu)+position2
          call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
          grid%velocity(i1,i2,i3) = VECTOR(0.,0.,0.)
       enddo
    enddo






    fac = shockFac

    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             rVec = VECTOR(grid%xAxis(i), grid%yAxis(j), 0.)
             vVec = rVec .cross. zAxis
             call normalize(vVec)
             vel = (twoPi * (modulus(rVec)*1.e10) / period)/cSpeed
             grid%velocity(i,j,k) = grid%velocity(i,j,k) + (vel*vVec)
             if (shockCone(i,j,k)) then
                grid%velocity(i,j,k) = (grid%velocity(i,j,k) .dot. shockDirection(i,j,k))*shockDirection(i,j,k)
                grid%rho(i,j,k) = grid%rho(i,j,k) * fac
             endif
          enddo
       enddo
    enddo



    deallocate(shockCone)
    deallocate(leftOfCone)
    deallocate(shockDirection)

    write(*,'(a)') "Grid filling complete..."


  end subroutine fillGridBinary

  logical function outsideGrid(posVec, grid)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: posVec

    outSideGrid = .false.

    if ((posVec%x > grid%xAxis(grid%nx)) .or. &
         (posVec%y > grid%yAxis(grid%ny)) .or. &
         (posVec%z > grid%zAxis(grid%nz)) .or. &
         (posVec%x < grid%xAxis(1)) .or. &
         (posVec%y < grid%yAxis(1)) .or. &
         (posVec%z < grid%zAxis(1))) outsideGrid = .true.

  end function outsideGrid


  subroutine getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec
    integer :: i1, i2, i3
    real :: t1, t2 ,t3
    real :: r, theta, phi, mu

    if (grid%cartesian) then

       !       if (rVec%x < grid%xAxis(1)) then
       !          write(*,*) rVec%x,grid%xAxis(1)
       !          stop
       !       endif
       !       if (rVec%y < grid%yAxis(1)) then
       !          write(*,*) rVec%y,grid%yAxis(1)
       !          stop
       !       endif
       !       if (rVec%z < grid%zAxis(1)) then
       !          write(*,*) rVec%z,grid%zAxis(1)
       !          stop
       !       endif
       !       if (rVec%x > grid%xAxis(grid%nx)) then
       !          write(*,*) rVec%x,grid%xAxis(grid%nx)
       !          stop
       !       endif
       !       if (rVec%y > grid%yAxis(grid%ny)) then
       !          write(*,*) rVec%y,grid%yAxis(grid%ny)
       !          stop
       !       endif
       !       if (rVec%z > grid%zAxis(grid%nz)) then
       !          write(*,*) rVec%z,grid%zAxis(grid%nz)
       !          stop
       !       endif


       call hunt(grid%xAxis, grid%nx, rvec%x, i1)
       call hunt(grid%yAxis, grid%ny, rvec%y, i2)
       call hunt(grid%zAxis, grid%nz, rvec%z, i3)
       if (i1 == grid%nx) i1 = i1 - 1
       if (i2 == grid%ny) i2 = i2 - 1
       if (i3 == grid%nz) i3 = i3 - 1

       t1 = (rVec%x-grid%xAxis(i1))/(grid%xAxis(i1+1)-grid%xAxis(i1))
       t2 = (rVec%y-grid%yAxis(i2))/(grid%yAxis(i2+1)-grid%yAxis(i2))
       t3 = (rVec%z-grid%zAxis(i3))/(grid%zAxis(i3+1)-grid%zAxis(i3))

!       if ((abs(t1) > 1.) .or. (abs(t2) > 1.) .or. (abs(t3) > 1.)) then
!          write(*,*) "bug in getindices"
!          write(*,*) "t1,t2,t3",t1, t2, t3
!          write(*,*) "rVec",rVec
!          write(*,*) "xaxis",grid%xAxis(1),grid%xAxis(grid%nx)
!          write(*,*) "yaxis",grid%yAxis(1),grid%yAxis(grid%ny)
!          write(*,*) "zaxis",grid%zAxis(1),grid%zAxis(grid%nz)
!       endif

    else

       call getPolar(rVec, r, theta, phi)
       mu = rVec%z/r
       call hunt(grid%rAxis, grid%nr, r, i1)
       if (i1 == 0) then
          i1 = 1
       endif
       if (i1 == grid%nr) i1 = grid%nr-1

       t1 = (r-grid%rAxis(i1))/(grid%rAxis(i1+1)-grid%rAxis(i1))


       call hunt(grid%muAxis, grid%nMu, mu, i2)
       if (i2 == grid%nMu) i2 = grid%nMu-1
       t2 = (mu-grid%muAxis(i2))/(grid%muAxis(i2+1)-grid%muAxis(i2))

       call hunt(grid%phiAxis, grid%nPhi, phi, i3)
       if (i3 == grid%nPhi) i3=i3-1
       t3 = (phi-grid%phiAxis(i3))/(grid%phiAxis(i3+1)-grid%phiAxis(i3))

    endif

  end subroutine getIndices

  subroutine writeGridPopulations(filename, grid, nLevels)

    type(GRIDTYPE) :: grid
    character(len=*) :: filename
    integer ::  nLevels, i, j, k, m

    write(*,'(a,a)') "Writing populations file to: ",trim(filename)
    open(20, file=filename, form="formatted", status="unknown")
    write(20,'(3i5)') grid%na1, grid%na2, grid%na3
    write(20,'(i5)') nLevels
    do i = 1, grid%na1
       do j = 1, grid%na2
          do k = 1, grid%na3
             write(20,'(e12.5)') grid%rho(i,j,k)
             do m = 1, nLevels
                write(20,'(e12.5)') grid%N(i,j,k,m)
             enddo
             write(20,'(e12.5)') grid%Ne(i,j,k)
             write(20,'(e12.5)') grid%temperature(i,j,k)
          enddo
       enddo
    enddo
    close(20)

  end subroutine writeGridPopulations

  subroutine readGridPopulations(filename, grid, nLevels)

    type(GRIDTYPE) :: grid
    character(len=*) :: filename
    integer :: nLevels
    integer :: na1, na2, na3, i, j, k, m

    open(20, file=filename, form="formatted", status="old")
    read(20,'(3i5)',err=100) na1, na2, na3
    read(20,'(i5)') nLevels
    if ( (na1 /= grid%na1).or.(na2 /= grid%na2).or.(na3 /= grid%na3) ) then
       write(*,'(a)') "! grid wrong size for populations file"
       stop
    endif
    write(*,'(a,a)') "Reading populations file from: ",trim(filename)
    do i = 1, grid%na1
       do j = 1, grid%na2
          do k = 1, grid%na3
             read(20,'(e12.5)') grid%rho(i,j,k)
             do m = 1, nLevels
                read(20,'(e12.5)') grid%N(i,j,k,m)
             enddo
             read(20,'(e12.5)') grid%Ne(i,j,k)
             read(20,'(e12.5)') grid%temperature(i,j,k)
          enddo
       enddo
    enddo
    close(20)
    goto 666

100 continue
    close(20)
    write(*,'(a)') "! Doesn't look like a text file, using unformatted input"
    open(20, file=filename, form="unformatted", status="old")
    read(20) na1, na2, na3
    read(20) nLevels
    if ( (na1 /= grid%na1).or.(na2 /= grid%na2).or.(na3 /= grid%na3) ) then
       write(*,'(a)') "! grid wrong size for populations file"
       stop
    endif
    write(*,'(a,a)') "Reading populations file from: ",trim(filename)
    do i = 1, grid%na1
       do j = 1, grid%na2
          do k = 1, grid%na3
             read(20) grid%rho(i,j,k)
             do m = 1, nLevels
                read(20) grid%N(i,j,k,m)
             enddo
             read(20) grid%Ne(i,j,k)
             read(20) grid%temperature(i,j,k)
          enddo
       enddo
    enddo
    close(20)


666 continue
  end subroutine readGridPopulations


  subroutine fillGridRolf(grid, mdot, vterm, o6width)

    type(GRIDTYPE) :: grid
    character(len=80) :: tmp
    integer :: i, j, k, ix, iy, iz, ierr
    real :: sigmaScaBlue, sigmaAbsBlue
    real :: sigmaScaRed, sigmaAbsRed
    real :: mdot, vterm, rho
    real :: t1, t2, t3, o6width
    integer, parameter :: maxr = 200
    
    integer :: kPrime
    real, allocatable :: ionPattern(:,:,:,:)
    real, allocatable :: rArray(:), yArray(:), yarray2(:)
    integer, allocatable :: ir(:,:)
    integer :: nTheta, nPhi, nr
    integer :: iTheta, iPhi, i1
    real :: fac, fac2, r, theta, phi
    type(VECTOR) :: starPos1, starPos2, rVec, rHat
    character(len=80) :: filename


    grid%tempSource = (o6width/2.35)**2*(16.*mHydrogen)/kErg

    write(*,*) "Source temperature: ",grid%tempSource


    grid%geometry = "rolf"
    grid%etaCont = 0.

    grid%inUse = .true.

    starPos1 = VECTOR(0.,2.09484334e13,0.)
    starPos2 = VECTOR(0.,-0.89779e13,0.)


    sigmaScaBlue = (34.+6.6)*sigmaE
    sigmaAbsBlue = 1.e-10 * sigmaScaBlue

    sigmaScaRed = 1.e-10 * sigmaE
    sigmaAbsRed = 1.e-10 * sigmaE


    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1)-1.
    enddo
    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1)-1.
    enddo
    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1)-1.
    enddo

    grid%xAxis = grid%xAxis * 1.e14
    grid%yAxis = grid%yAxis * 1.e14
    grid%zAxis = grid%zAxis * 1.e14


    write(*,'(a)') "Reading hydro grids..."
    open(20, file="density.dat", status="old", form="formatted")
10  continue
    read(20,'(a)') tmp
    if (tmp(1:5) /= "*DATA") goto 10
    read(20,*) ix, iy, iz
    if ((grid%nx /= ix).or.(grid%ny /= iy).or.(grid%nz /= iz)) then
       write(*,'(a)') "! grid is not correct size for hydro data"
       write(*,'(a,3i4)') "! it should be: ",ix,iy,iz
    endif
    read(20,*) (((grid%rho(i,j,k),i=1,ix),j=1,iy),k=1,iz)
    close(20)


    allocate(grid%kappaAbsRed(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       stop
    endif
    allocate(grid%kappaScaRed(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1),stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       stop
    endif

    open(20, file="velocity_field.dat", status="old", form="formatted")
20  continue
    read(20,'(a)') tmp
    if (tmp(1:5) /= "*DATA") goto 20
    read(20,*) ix, iy, iz
    if ((grid%nx /= ix).or.(grid%ny /= iy).or.(grid%nz /= iz)) then
       write(*,'(a)') "! grid is not correct size for hydro data"
       write(*,'(a,3i4)') "! it should be: ",ix,iy,iz
    endif
    read(20,*) (((grid%velocity(i,j,k)%x, grid%velocity(i,j,k)%y, &
         grid%velocity(i,j,k)%z, i = 1, ix), j = 1, iy), k = 1, iz)
    close(20)
    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             grid%velocity(i,j,k) = grid%velocity(i,j,k) / cSpeed
          enddo
       enddo
    enddo
    write(*,'(a)') "Velocity done..."


    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             if (grid%rho(i,j,k) > 0.) then
                grid%rho(i,j,k) = grid%rho(i,j,k)
                grid%kappaAbs(i,j,k,1) = grid%rho(i,j,k) * sigmaAbsBlue
                grid%kappaSca(i,j,k,1) = grid%rho(i,j,k) * sigmaScaBlue
                grid%kappaAbsRed(i,j,k,1) = grid%rho(i,j,k) * sigmaAbsRed
                grid%kappaScaRed(i,j,k,1) = grid%rho(i,j,k) * sigmaScaRed
                grid%etaLine(i,j,k) = (grid%rho(i,j,k)**2)*1.e-22
             else
                if (grid%yAxis(j) > 0.) then
                   rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))-starPos1
                   r = modulus(rVec)
                   call normalize(rVec)
                   grid%velocity(i,j,k) = (vterm / cSpeed) * rVec
                   if (r > 0.) then
                      rho = (mdot / (fourPi * r**2 * vterm))/mHydrogen
                   else
                      rho = 1.e13
                   endif
                   grid%rho(i,j,k) = rho
                   grid%kappaAbs(i,j,k,1) = grid%rho(i,j,k) * sigmaAbsBlue
                   grid%kappaSca(i,j,k,1) = grid%rho(i,j,k) * sigmaScaBlue
                   grid%kappaAbsRed(i,j,k,1) = grid%rho(i,j,k) * sigmaAbsRed
                   grid%kappaScaRed(i,j,k,1) = grid%rho(i,j,k) * sigmaScaRed
                   grid%etaLine(i,j,k) = (grid%rho(i,j,k)**2)*1.e-22
                   grid%velocity(i,j,k) = VECTOR(0.,0.,0.)
                else
                   grid%rho(i,j,k) = 1.e14
                   grid%kappaAbs(i,j,k,1) = grid%rho(i,j,k) * sigmaAbsBlue
                   grid%kappaSca(i,j,k,1) = grid%rho(i,j,k) * sigmaScaBlue
                   grid%kappaAbsRed(i,j,k,1) = grid%rho(i,j,k) * sigmaAbsRed
                   grid%kappaScaRed(i,j,k,1) = grid%rho(i,j,k) * sigmaScaRed
                   grid%etaLine(i,j,k) = (grid%rho(i,j,k)**2)*1.e-22
                   grid%inUse(i,j,k) = .false.
                endif
             endif
          enddo
       enddo
    enddo
    write(*,'(a)') "Density done..."


    write(*,'(a)') "Reading ionization structure ray files..."
    nTheta = 51
    nPhi = 99
    nr = maxr
    allocate(ionPattern(1:nTheta, 1:nPhi, 1:nr, 3))
    allocate(iR(1:nTheta, 1:nPhi))
    do i = 1, nTheta
       if (i < ntheta) then
          iphi = nphi
       else
          iphi = 1
       endif
       do j = 1, iPhi
          write(filename,'(a,i2.2,a,i2.2)') "ion/fort.ion.it=",i,".ip=",j
          open(20, file=filename, status="old",form="formatted")
          nr = 0 
50        continue
          read(20, '(a)',end=60) tmp
          if (tmp(2:4) == "r =") then
             nr = nr + 1
             if (nr > maxr) then
                write(*,'(a)') "! maxr is too small"
             endif
             read(tmp(5:),*) ionPattern(i,j,nr,1)
             read(20,'(a)') tmp
             read(tmp(5:),*) ionPattern(i,j,nr,2)
             do while(tmp(1:3) /= " O ")
                read(20,'(a)') tmp
             enddo
             read(tmp(65:),*) ionPattern(i,j,nr,3)       ! O VII
          endif
          goto 50
60        continue
          close(20)
          ir(i,j) = nr
!          write(*,'(a,a,a,i3)') "Read file: ",trim(filename), " nr=",nr
!          write(*,*) ionPattern(i,j,1,1),ionPattern(i,j,1,2), ionPattern(i,j,1,3)
       enddo
    enddo

    allocate(rArray(1:maxr))
    allocate(yArray(1:maxr))
    allocate(yArray2(1:maxr))


    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             rVec = VECTOR(grid%xAxis(i), grid%yAxis(j), grid%zAxis(k)) - starPos1
             call getPolar(rVec, r, theta, phi)
             if (theta <= pi/2.) then
                theta = pi/2. - theta
             else
                theta = theta - pi/2.
             endif


             iTheta = min(nint( (theta/(pi/2.))*real(nTheta) ) + 1, nTheta)
             if (iTheta /= nTheta) then
                iPhi = min(nint( (phi/twoPi) * real(nPhi)) + 1, nPhi)
             else
                iPhi = 1
             endif


             rArray(1:ir(itheta,iphi)) = ionPattern(itheta,iphi,1:ir(itheta,iphi),1)
             yArray(1:ir(itheta,iphi)) = ionPattern(itheta,iphi,1:ir(itheta,iphi),2)
             yArray2(1:ir(itheta,iphi)) = ionPattern(itheta,iphi,1:ir(itheta,iphi),3)

             if ( (r > rArray(1)) .and. (r < rArray(ir(itheta,iphi)))) then
                call locate(rArray, ir(itheta,iphi), r, i1)
                fac = logint(r, rArray(i1), rArray(i1+1), yArray(i1), yArray(i1+1))
                fac2 = logint(r, rArray(i1), rArray(i1+1), yArray2(i1), yArray2(i1+1))
             endif
             if (r > rArray(ir(itheta,iphi))) then
                fac = 1.
                fac2 = yArray2(ir(iTheta,iPhi))
             endif
             if (r < rArray(1)) then
                fac = 1.e-10
                fac2 = yArray2(1)
             endif
             grid%kappaSca(i,j,k,1) = max(1.e-30,grid%kappaSca(i,j,k,1) * fac)
             grid%kappaAbs(i,j,k,1) = max(1.e-30,grid%kappaAbs(i,j,k,1) * fac)
             grid%kappaScaRed(i,j,k,1) = max(1.e-30,grid%kappaScaRed(i,j,k,1) * fac)
             grid%kappaAbsRed(i,j,k,1) = max(1.e-30,grid%kappaAbsRed(i,j,k,1) * fac)
             grid%etaLine(i,j,k) = grid%etaLine(i,j,k)*fac2**2
             grid%rho(i,j,k) = grid%rho(i,j,k) * fac
          enddo
       enddo
    enddo
    deallocate(ionPattern)
    deallocate(ir)
    deallocate(rArray)
    deallocate(yArray)
    deallocate(yArray2)


    write(*,'(a)') "Ionization structure completed"
    write(*,'(a)') "Ensuring symmetry..."

    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz/2
             kPrime = grid%nz - k + 1
             grid%rho(i,j,k) = grid%rho(i,j,kPrime)
             grid%kappaSca(i,j,k,1) = grid%kappaSca(i,j,kPrime,1)
             grid%kappaScaRed(i,j,k,1) = grid%kappaScaRed(i,j,kPrime,1)
             grid%kappaAbs(i,j,k,1) = grid%kappaAbs(i,j,kPrime,1)
             grid%kappaAbsRed(i,j,k,1) = grid%kappaAbsRed(i,j,kPrime,1)
             grid%etaLine(i,j,k) = grid%etaLine(i,j,kPrime)
             grid%velocity(i,j,k) = grid%velocity(i,j,kPrime)
             grid%velocity(i,j,k)%z = -grid%velocity(i,j,k)%z
!             write(*,'(3i3,1p,6e12.5)') i,j,k,(cSpeed/1.e5)*grid%velocity(i,j,k),(cSpeed/1.e5)*grid%velocity(i,j,kPrime)
          enddo
       enddo
    enddo

    grid%rho = grid%rho * mHydrogen

    write(*,'(a,1pe12.5)') "Neutral hydrogen mass: ",integratedDensity(grid)

!    grid%etaLine = 1.e-20
!    call getIndices(grid, starPos1, i, j, k, t1, t2, t3)
!    write(*,*) "etaLine",i,j,k
!    grid%etaLine(i-1:i+1,j-1:j+1,k-1:k+1) = 1.e10
!    grid%velocity(i-1:i+1,j-1:j+1,k-1:j+1) = VECTOR(0.,0.,0.)

!    write(*,'(a)') "Fiddling velocity field."
!    do i = 1, grid%nx
!       do j = 1, grid%ny
!          do k = 1, grid%nz
!             rHat =VECTOR(grid%xAxis(i), grid%yAxis(j),grid%zAxis(k))
!             call normalize(rHat)
!             grid%velocity(i,j,k) = (100.e5/cSpeed) * rHat
!          enddo
!       enddo
!    enddo




    write(*,'(a)') "Done."


  end subroutine fillGridRolf

  subroutine fillGridWR137(grid, rCore, mDot, vTerm, beta, temp)

    type(GRIDTYPE) :: grid
    real :: rCore, mDot, vTerm, beta, temp
    integer :: i,j,k
    real :: sinTheta, vel, tot, v0, fac
    type(VECTOR) :: rHat
    integer :: i1, i2
    real :: lineTot, contTot

    grid%lineEmission = .true.

    do i = 1, grid%nr
       grid%rAxis(i) = log10(rCore) + 2.*real(i-1)/real(grid%nr)
       grid%rAxis(i) = 10.**grid%rAxis(i)
    enddo


    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nPhi-1)
    enddo
    do i = 1, grid%nMu
       grid%muAxis(i) = 2.*real(i-1)/real(grid%nMu-1)-1.
    enddo

    grid%rCore = rCore

    v0 = 0.01 * vTerm

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             sinTheta = sqrt(1.-grid%muAxis(j)**2)
             rHat = VECTOR(sinTheta*cos(grid%phiAxis(k)), sinTheta*sin(grid%phiAxis(k)), grid%muAxis(j))
             vel = v0 + (vterm-v0)*(1.-grid%rCore / grid%rAxis(i))**beta
             grid%rho(i,j,k) = (mdot / (fourPi * grid%rAxis(i)**2 * vel))/mHydrogen
             grid%velocity(i,j,k) = (vel/cSpeed) * rHat
             grid%kappaAbs(i,j,k,1) = (1.e-20*grid%rho(i,j,k))**2
             grid%kappaSca(i,j,k,1) = grid%rho(i,j,k) * sigmaE / 4.
          enddo
       enddo
    enddo
    grid%temperature = temp

    grid%rAxis(1:grid%nr) =  grid%rAxis(1:grid%nr) / 1.e10
    grid%rCore = grid%rCore / 1.e10
    grid%kappaSca(1:grid%nr,1:grid%nMu,1:grid%nPhi,1) = &
         grid%kappaSca(1:grid%nr,1:grid%nMu,1:grid%nPhi,1) * 1.e10


    tot = 0.
    do i = 1 , grid%nr-1
       tot = tot + (grid%rAxis(i+1)-grid%rAxis(i))*grid%kappaAbs(i,1,1,1)
    enddo
    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             grid%kappaabs(i,j,k,1) = grid%kappaAbs(i,j,k,1) * 20./tot
          enddo
       enddo
    enddo
    tot = 0.
    do i = 1 , grid%nr-1
       tot = tot + (grid%rAxis(i+1)-grid%rAxis(i))*grid%kappaAbs(i,1,1,1)
    enddo


    grid%etaCont(1:grid%nr,1:grid%nMu,1:grid%nPhi) = grid%kappaAbs(1:grid%nr,1:grid%nMu,1:grid%nPhi,1)


    call locate(grid%rAxis, grid%nr,  5.*grid%rCore, i1)
    call locate(grid%rAxis, grid%nr, 10.*grid%rCore, i2)

    grid%etaLine = 1.e-30
    grid%etaLine(i1:i2,1:grid%nMu,1:grid%nPhi) = 1.

    call integratedEmission(grid, lineTot, contTot)


    fac = contTot*1.e12/lineTot
    grid%etaLine(1:grid%nr, 1:grid%nmu, 1:grid%nPhi) = &
         grid%etaLine(1:grid%nr, 1:grid%nmu, 1:grid%nPhi) * fac 

    grid%chiLine  = 1.e-20

  end subroutine fillGridWR137


  subroutine integratedEmission(grid, line, cont)
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    real :: dV,dx,dy,dz
    real :: dr, dTheta, dPhi, sinTheta, mu
    real :: wtot, w1, w2
    real :: line, cont
    line = 0.
    cont = 0.

    if (grid%cartesian) then
       do i = 1, grid%nx-1
          do j = 1, grid%ny - 1
             do k = 1, grid%nz - 1
                dx = grid%xAxis(i+1)-grid%xAxis(i)
                dy = grid%yAxis(j+1)-grid%yAxis(j)
                dz = grid%zAxis(k+1)-grid%zAxis(k)
                cont = cont + grid%etaCont(i,j,k)*dx*dy*dz
                line = line + grid%etaLine(i,j,k)*dx*dy*dz
             enddo
          enddo
       enddo
    else
       do i = 1,grid%nr-1
          do j = 1, grid%nMu-1
             do k = 1, grid%nPhi-1
                dTheta = acos(grid%muAxis(j+1))-acos(grid%muAxis(j))
                dPhi = grid%phiAxis(k+1) - grid%phiAxis(k)
                mu = 0.5*(grid%muAxis(j+1)+grid%muAxis(j))
                sinTheta = sqrt(1. - mu**2)
                dr = grid%rAxis(i+1) - grid%rAxis(i)
                dV = abs(sinTheta * dr * dTheta * dPhi)
                w1 = grid%rAxis(i)**2
                w2 = grid%rAxis(i+1)**2
                wtot = w1 + w2
                w1 = w1 / wtot
                w2 = w2 / wtot
                cont = cont + (w1*grid%etaCont(i,j,k)+w2*grid%etaCont(i+1,j,k)) * dV
                line = line + (w1*grid%etaCont(i,j,k)+w2*grid%etaLine(i+1,j,k)) * dV
             enddo
          enddo
       enddo
    endif

  end subroutine integratedEmission

  subroutine fillGridPlanet(grid)
    type(GRIDTYPE) :: grid
    real :: rPlanet, scaleHeight
    integer :: i, j, k
    real :: albedo = 0.5
    real :: fac, tot

    rPlanet = 6000.e3
    scaleHeight = 10.e3

    grid%geometry = "planet"

    do i = 1, grid%nr
       grid%rAxis(i) = log10(rPlanet) + log10((rPlanet+20.*scaleHeight)/rPlanet)*real(i-1)/real(grid%nr-1)
       grid%rAxis(i) = 10.**grid%rAxis(i)
       write(*,*) i, grid%rAxis(i)/grid%rAxis(1)
    enddo

    grid%rCore = grid%rAxis(1)

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nPhi-1)
    enddo
    do i = 1, grid%nMu
       grid%muAxis(i) = 2.*real(i-1)/real(grid%nMu-1)-1.
    enddo

    do i = 1, grid%nr
       fac = exp(-(grid%rAxis(i) - grid%rAxis(1))/scaleHeight)
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             grid%velocity(i,j,k) = VECTOR(0.,0.,0.)
             grid%rho(i,j,k) = max(fac,1.e-20)
             grid%kappaSca(i,j,k,1) = grid%rho(i,j,k)
             grid%kappaAbs(i,j,k,1) = grid%kappaSca(i,j,k,1)*(1.-albedo)/albedo
          enddo
       enddo
    enddo
    tot = 0.
    do i = 1, grid%nr-1
       tot = tot + (grid%rAxis(i+1)-grid%rAxis(i))*(grid%kappaSca(i,1,1,1)+grid%kappaAbs(i,1,1,1))
    enddo
    fac = 10./tot
    grid%kappaSca(1:grid%nr,1:grid%nmu,1:grid%nphi,1) = grid%kappaSca(1:grid%nr,1:grid%nmu,1:grid%nphi,1) * fac
    grid%kappaAbs(1:grid%nr,1:grid%nmu,1:grid%nphi,1) = grid%kappaAbs(1:grid%nr,1:grid%nmu,1:grid%nphi,1) * fac
  end subroutine fillGridPlanet

  subroutine fillGridHourglass(grid)
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    real :: gridExtent = 1.e6
    integer :: nAng = 400, nr = 20
    real :: rInner, rOuter
    integer :: i1, i2, i3
    real :: t1, t2, t3
    real :: dz
    real :: phi1, phi2, phi3, x, y, z, x2y2
    real :: zp, rw, h, ang, r
    real :: dr
    type(VECTOR) :: sVec, norm, sVec1, sVec2, norm2
    logical :: radial, test

    radial = .true.
    test = .false.
    grid%geometry = "hourglass"

    rInner = 0.5 * gridExtent
    rOuter = 0.6 * gridExtent

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1)-1.
    enddo
    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1)-1.
    enddo
    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1)-1.
    enddo


    grid%xAxis = grid%xAxis * 2.9e7
    grid%yAxis = grid%yAxis * 2.9e7
    grid%zAxis = grid%zAxis * 2.9e7

    grid%kappaSca = 1.e-20
    grid%kappaAbs = 1.e-20
    grid%chiLine = 1.e-20
    grid%etaLine = 1.e-20
    grid%etaCont = 0.
    grid%velocity = VECTOR(0.,0.,0.)

!    do m = 1, nr
!       hourGlassRadius = rInner + (rOuter-rInner) * real(m-1)/real(nr-1)
!       do i = 1, nAng
!          rotAng = twoPi * real(i-1)/real(nAng-1)
!          do k = 1, grid%nz
!             z = grid%zAxis(k)
!             dz = abs(rOuter - abs(z))
!             if (abs(z) < rOuter) then
!                if (dz < hourGlassRadius) then
!                   y = sqrt(hourGlassRadius**2 - dz**2)
!                   rVec = VECTOR(0., y, z)
!                   if (z >= 0.) then
!                      sVec = rVec - VECTOR(0.,0.,rOuter)
!                      vVec = sVec .cross. VECTOR(-1.,0.,0.)
!                      call normalize(vVec)
!                   else
!                      sVec = rVec - VECTOR(0.,0.,-rOuter)
!                      vVec = sVec .cross. VECTOR(1.,0.,0.)
!                      call normalize(vVec)
!                   endif
!                endif
!             else
!                y = hourglassRadius
!                rVec = VECTOR(0., y, z)
!                if (z >= 0.) then
!                   vVEC = VECTOR(0.,0.,1.)
!                else
!                   vVec = VECTOR(0.,0.,-1.)
!                endif
!             endif
!
!
!             sVec = rotateZ(rVec, rotAng)
!             tmp = rotateZ(vVec, rotAng)
!             vVec = tmp
!             call getIndices(grid, sVec, i1, i2, i3, t1, t2, t3)
!             grid%etaLine(i1,i2,i3) = 1./max(modulus(sVec),grid%zAxis(2)-grid%zAxis(1))**2
!             vel = 90.e5/cSpeed
!             grid%velocity(i1,i2,i3) = vel * vVec
!          enddo
!       enddo
!    enddo
!
!! radial velocity field
!
!    if (radial) then
!       do i = 1, grid%nx
!          do j = 1, grid%ny
!             do k = 1, grid%nz
!                rVec = VECTOR(grid%xAxis(i), grid%yAxis(j), grid%zAxis(k))
!                vel = modulus(rVec)/sqrt(3.*grid%xAxis(grid%nx)**2)
!                call normalize(rVec)
!                grid%velocity(i,j,k) = (vel*90.e5/cSpeed) * rVec
!             enddo
!          enddo
!       enddo
!    endif

    rw = 2e6
    h = 2.9e7
    phi1 = 44. * degToRad
    phi2 = 62. * degToRad
    phi3 = 85. * degToRad
    dr = 1.6e6

    zp  = (h * (tan(phi3)/tan(phi2)) - 1.)/((tan(phi3)/tan(phi1))-1.)
    
    do k = 1, grid%nz
       z = grid%zAxis(k)
       if (abs(z) < zp) then
          x2y2 = sqrt(((zp / tan(phi1)**2) - (rw**2/zp))*abs(z) + rw**2)
          dz = (((zp / tan(phi1)**2) - (rw**2/zp))/2.)* & 
               (((zp / tan(phi1)**2) - (rw**2/zp))*abs(z) + rw**2)**(-0.5)
       else
          x2y2 = (abs(z) / tan(phi3)) - zp*(1./tan(phi3) - 1./tan(phi1))
          dz = 1./tan(phi3)
       endif
       if (z < 0.) dz = -dz
       norm = VECTOR(1.,0.,dz)
       norm = norm .cross. VECTOR(0.,1.,0.)
       call normalize(norm)
       
       do i = 1, nAng
          ang = twoPi*real(i-1)/real(nAng-1)
          x = cos(ang)*x2y2
          y = sin(ang)*x2y2
          norm2 = rotateZ(norm, ang)
          sVec = VECTOR(x,y,z)
          sVec1 = sVec - (dr/2.)*norm2
          sVec2 = sVec + (dr/2.)*norm2
          do j = 1, nr
             sVec = sVec1 + (real(j-1)/real(nr-1))*(sVec2-sVec1)
             call getIndices(grid, sVec, i1, i2, i3, t1, t2, t3)
             r = modulus(sVec)
             if ((r >= 2.e6).and.(r < 3.e6)) then
                grid%etaLine(i1,i2,i3) = (r / 3.e6)**(-4.)
             else if ((r >= 3.e6).and.(r < 1.3e7)) then
                grid%etaLine(i1,i2,i3) = (r / 1.3e7)**(-4.)
             else if ((r >= 1.3e7).and.(r < 2.e7)) then
                grid%etaLine(i1,i2,i3) = (r / 1.3e7)**(+3.)
             else if (r >= 2.e7) then
                grid%etaLine(i1,i2,i3) = (r / 2.e7)**(-3.)
             endif
             call normalize(sVec)
             grid%velocity(i1,i2,i3) = ((24.e5)*((r/1.3e7)**0.6)/cSpeed) * sVec
           enddo
       enddo
    enddo
    if (test) then
       grid%etaLine(1:grid%nx,1:45,1:grid%nz) = 1.e-20
       grid%etaLine(1:grid%nx,55:100,1:grid%nz) = 1.e-20
    endif



  end subroutine fillGridHourglass


  subroutine fillGridBetaCep(grid)


    implicit none
    type(GRIDTYPE) :: grid
    integer, parameter :: nRad = 100
    real :: rOuter, rInner, rCore
    type(VECTOR) :: rVec
    real :: phi, r
    real :: height
    integer, parameter :: nTheta = 101, nPhi = 200
    real :: x, y, z, mCore
    type(VECTOR) :: spinAxis
    integer :: i, j, k
    integer :: nH = 100
    integer :: nr = 100
    integer :: i1, m1, m2, m3
    real :: t1, t2, t3
    real :: heightArray(6), radiiArray(6)

    radiiArray(1) = 3.
    radiiArray(2) = 4.
    radiiArray(3) = 5.
    radiiArray(4) = 6.
    radiiArray(5) = 7.
    radiiArray(6) = 8.

    heightArray(1) = 0.1
    heightArray(2) = 0.06
    heightArray(3) = 0.025
    heightArray(4) = 0.012
    heightArray(5) = 0.007
    heightArray(6) = 0.004


    rCore = 7.*rSol
    mCore = 12.*mSol

    grid%geometry = "betacep"
    grid%lineEmission = .true.
    grid%rCore = rCore


    spinAxis = VECTOR(0.,0.,1.)

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1) - 1.
    enddo

    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1) - 1.
    enddo

    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1) - 1.
    enddo

    rInner = 3.*rCore
    rOuter = 8.*rCore

    grid%rho = 0.
    grid%etaLine = 1.e-20
    grid%etaCont = 1.e-20
    grid%chiLine = 1.e-20
    grid%kappaSca = 1.e-20
    grid%kappaAbs = 1.e-20
    grid%xAxis = grid%xAxis * rOuter*1.1
    grid%yAxis = grid%yAxis * rCore
    grid%zAxis = grid%zAxis * rOuter*1.1
    grid%velocity = VECTOR(1.e-30,1.e-30,1.e-30)


    grid%temperature = 0.

    do i = 1, nPhi

       phi = twoPi*real(i-1)/real(nPhi-1)

       do j = 1, nr
          r = rInner+(rOuter-rInner)*real(j-1)/real(nr-1)
          call locate(radiiArray, 6, r/rCore, i1)
          height = heightArray(i1) + &
                  (r/rCore - radiiArray(i1)) * &
         (heightArray(i1+1)-heightArray(i1))/(radiiArray(i1+1)-radiiArray(i1))
          do k = 1, nH
             y = height*(2.*real(k-1)/real(nh-1)-1.)*rCore/2.
             x = cos(phi)*r
             z = sin(phi)*r
             rVec = VECTOR(x,y,z)
             call getIndices(grid, rVec, m1, m2, m3, t1, t2, t3)
             grid%rho(m1,m2,m3) = 3.e12
          enddo
       enddo
    enddo



    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             if (grid%rho(i,j,k) /= 0.) then
                grid%kappaSca(i,j,k,1) = max(1.e-30,grid%rho(i,j,k) * sigmaE)
             endif
          enddo
       enddo
    enddo

  end subroutine fillGridBetaCep


  logical function oddBall(grid,i1,i2,i3)
    type(GRIDTYPE) :: grid
    integer :: i1, i2, i3
    integer :: j1, j2, j3
    integer :: k1, k2, k3
    integer :: i, j, k
    oddBall = .true.
    
    j1 = max(i1-1,1)
    j2 = max(i2-1,1)
    j3 = max(i3-1,1)

    k1 = min(i1+1,grid%na1)
    k2 = min(i2+1,grid%na2)
    k3 = min(i3+1,grid%na3)

    do i = j1, k1
       do j = j2, k2
          do k = j3, k3
             if (grid%inUse(i,j,k)) then
                oddball = .false.
                EXIT
             endif
          enddo
          if (.not.oddBall) EXIT
       enddo
       if (.not.oddBall) EXIT
    enddo
  end function oddBall



  subroutine fillGridDonati(grid, resonanceLine)
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    integer :: i1, i2, i3, nMu, j1, j2
    real :: zPrime, yPrime, r, r1, fac, ang
    logical :: resonanceLine
    type(VECTOR) :: rHat, rVec, vHat
    real :: surfaceVel, vel, theta
    integer :: nz, ny
    real :: t1, t2, t3, sinTheta
    real :: thickness
    real, allocatable :: z(:), y(:), rho(:,:), t(:,:), phi(:,:), v(:,:)
    character(len=80) :: junkline
    integer, parameter :: ndisk = 8
    real :: rStar
    real :: diskDensity
    real :: xDisk(ndisk) = (/1.13, 1.22, 1.33, 1.49, 1.70, 2.00, 2.50, 3.0/)
    real :: yDisk(ndisk) = (/0.021, 0.055, 0.052, 0.031, 0.018, 0.010, 0.0046, &
                                                                   0.0025 /)
    real :: rhoDisk(ndisk) = (/5., 7.6, 9.5, 10.1, 10.0, 9.4, 8.5, 7.6 /)
    


    grid%geometry = "donati"
    grid%lineEmission = .true.


    diskDensity = 1.5e12

    rStar = 8.2 * rSol

    surfaceVel = 250.*1e5

    write(*,'(a)') "Reading Donati maps..."

    write(*,'(a)') "Density..."
    open(20, file = "rho35_closed.map", status = "old", form="formatted")
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    allocate(z(1:nz))
    allocate(y(1:ny))
    allocate(rho(ny, nz))
    allocate(t(ny, nz))
    allocate(v(ny, nz))
    allocate(phi(ny, nz))
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((rho(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Velocity..."
    open(20, file = "vel35_closed.map", status = "old", form="formatted")
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((v(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Phi..."
    open(20, file = "phi35_closed.map", status = "old", form="formatted")
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((phi(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Temperature..."
    open(20, file = "temp35_closed.map", status = "old", form="formatted")
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((t(i,j),i=1,nz),j=1,ny)
    close(20)
    write(*,'(a)') "Done."

    z = 10.**z
    v = 10.**v
    rho = 10.**rho
    t = 10.**t


    do i = 1, grid%nr
       grid%rAxis(i) = 1. + 2.*real(i-1)/real(grid%nr-1)
    enddo
    
    theta = atan(0.5*yDisk(2)/xDisk(2))
    call fillmuAxisDisk(grid%muAxis, grid%nMu, theta)


    theta = atan(0.5*yDisk(8)/xDisk(8))
    grid%muAxis(grid%nMu/2+1-1) = cos(pi/2. + theta)
    grid%muAxis(grid%nMu/2+1+1) = cos(pi/2. - theta)

    theta = atan(0.5*yDisk(7)/xDisk(7))
    grid%muAxis(grid%nMu/2+1-2) = cos(pi/2. + theta)
    grid%muAxis(grid%nMu/2+1+2) = cos(pi/2. - theta)
    
    


    call writeAxes(grid)

    y(1) = 0.
    z(1) = 0.

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi * real(i-1)/real(grid%nphi-1)
    enddo

    grid%rho = 1.e-30
    grid%temperature = 1.e-30

    grid%inUse = .false.

    write(*,'(a)') "Re-mapping grids..."
    do i1 = 1, grid%nr
       do i2 = 1, grid%nMu
          do i3 = 1, grid%nPhi
             r = grid%rAxis(i1)
             if (r > 1.1) then
                grid%inUse(i1,i2,i3) = .true.
                zPrime = abs(r * grid%muAxis(i2))
                yPrime = r * sqrt(1.d0 - grid%muAxis(i2)**2)
                call hunt(y, ny, yPrime, j)
                call hunt(z, nz, zPrime, i)
                grid%rho(i1, i2, i3) = rho(i,j)
                grid%temperature(i1, i2, i3) = t(i,j)
                ang = grid%phiAxis(i3)
                rHat%z = cos(phi(i,j))
                if ((r*grid%muAxis(i2)) < 0.) rHat%z = -rHat%z
                fac = sqrt(1.d0 - rHat%z**2)
                rHat%x = fac * cos(ang)
                rHat%y = fac * sin(ang)
                grid%velocity(i1, i2, i3) = (v(i,j)/cSpeed) * rHat
             endif
          enddo
       enddo
    enddo
    where((grid%temperature < 1.).or.(grid%rho < 1.e-19))
       grid%rho = 1.e-30
       grid%inUse = .false.
    endwhere
    write(*,'(a)') "Done."

 
    write(*,'(a)') "Adding disk..."

 

    do i1 = 1, grid%nr
       r = grid%rAxis(i1)
       if (r >= xDisk(1)) then
          call locate(xDisk, ndisk, r, j)
          thickness = yDisk(j) + (r - xDisk(j))*(yDisk(j+1)-yDisk(j))/ &
               (xDisk(j+1)-xDisk(j))
          diskDensity = rhoDisk(j) + (r - xDisk(j))*(rhoDisk(j+1)-rhoDisk(j))/ &
               (xDisk(j+1)-xDisk(j))

          rVec = VECTOR(r * cos(grid%phiAxis(i2)), &
               r1 * sin(grid%phiAxis(i2)), &
               -thickness/2.)
          call getIndices(grid, rVec, i, j1, k, t1, t2, t3)
          rVec = VECTOR(r * cos(grid%phiAxis(i2)), &
               r1 * sin(grid%phiAxis(i2)), &
               thickness/2.)
          call getIndices(grid, rVec, i, j2, k, t1, t2, t3)
          write(*,*) i1,j1,j2

          do i2 = 1, grid%nPhi
             do j = j1+1, j2
                grid%inUse(i1,j,i2) = .true.
                grid%rho(i1,j,i2) = (diskDensity*1.e12) * mHydrogen 
                grid%temperature(i1,j,i2) = 40000.
                zPrime = r * grid%muAxis(j)
                if (zPrime > 0.) then
                   vHat = VECTOR(0.,0.,-1.)
                else
                   vHat = VECTOR(0.,0.,+1.)
                endif
                vel = abs(zPrime) * (2.e3 * 1.e5) ! 2000km/s per rstar
                grid%velocity(i1,j,i2) = (vel/cSpeed) * vHat
             enddo
          enddo
       endif
    enddo

    write(*,'(a)') "Done."


    where (grid%temperature > 1.e6)
       grid%temperature = 1.e6
    end where

    grid%rAxis = grid%rAxis*rStar/1.e10

    grid%rCore = rStar/1.e10

    write(*,'(a,1p,2e12.5)') "Maximum density: ",&
         MAXVAL(grid%rho, mask=grid%inUse)/mHydrogen,MAXVAL(grid%rho, mask=grid%inUse)
    write(*,'(a,1p,2e12.5)') "Minimum density: ", &
         MINVAL(grid%rho, mask=grid%inUse)/mHydrogen,MINVAL(grid%rho, mask=grid%inUse)
    write(*,'(a,1p,2e12.5)') "Maximum/min temp: ",&
         MAXVAL(grid%temperature, mask=grid%inUse),MINVAL(grid%temperature, mask=grid%inUse)

    grid%rstar1 = grid%rCore
    grid%rstar2 = 0.
    grid%starpos1 = VECTOR(0.,0.,0.)
    grid%starpos2 = VECTOR(1.e20,0.,0.)

!    nMu = grid%nMu/2+1
!    do i = 1, nMu-1
!       grid%rho(1:grid%nr,i,1:grid%nPhi) = &
!            grid%rho(1:grid%nr,grid%nMu-i+1,1:grid%nPhi) 
!    enddo


    do i = 1, grid%nMu
       write(*,*) acos(grid%muAxis(i))*radtodeg, grid%rho(grid%nr/2,i,1), &
            grid%velocity(grid%nr/2,i,1)%z*cSpeed/1.e5
    enddo


    grid%etaLine = 1.e-30
    grid%etaCont = 1.e-30
    grid%chiLine = 1.e-30
    grid%kappaSca = 1.e-30
    grid%kappaAbs = 1.e-30

    write(*,'(a)') "Done."


    if (resonanceLine) then

       grid%etaLine = 1.e-30
       grid%etaCont = 1.e-30
       grid%chiLine = 1.e-30
       grid%kappaSca = 1.e-30
       grid%kappaAbs = 1.e-30
       grid%resonanceLine = .true.
       grid%lambda2 = 1548.202

       grid%chiLine = grid%rho*1.e26
       write(*,'(a)') "This is a resonance line model."

    endif



  end subroutine fillGridDonati

  subroutine fillGridDonati2(grid, resonanceLine, misc)
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    integer :: i1, i2, i3, nMu, j1, j2
    real :: zPrime, yPrime, r, r1, fac, ang
    logical :: resonanceLine
    character(len=*) :: misc
    type(VECTOR) :: rHat, rVec, vHat
    real :: surfaceVel, vel, theta
    integer :: nz, ny
    real :: t1, t2, t3, sinTheta
    real :: thickness
    real, allocatable :: z(:), y(:), rho(:,:), t(:,:), phi(:,:), v(:,:)
    character(len=80) :: junkline
    integer, parameter :: ndisk = 8
    real :: rStar
    real :: diskDensity
    real :: xDisk(ndisk) = (/1.13, 1.22, 1.33, 1.49, 1.70, 2.00, 2.50, 3.0/)
    real :: yDisk(ndisk) = (/0.021, 0.055, 0.052, 0.031, 0.018, 0.010, 0.0046, &
                                                                   0.0025 /)
    real :: rhoDisk(ndisk) = (/5., 7.6, 9.5, 10.1, 10.0, 9.4, 8.5, 7.6 /)

    real :: yDisk2(nDisk) = (/0.004, 0.011, 0.025, 0.060, 0.090, 0.050, 0.023, 0.0125 /)
    real :: rhoDisk2(nDisk) = (/ 1.0, 1.5, 1.9, 2.0, 2.0, 1.9, 1.7, 1.5 /)


    grid%geometry = "donati"
    grid%lineEmission = .true.


    diskDensity = 1.5e12

    rStar = 8.2 * rSol

    surfaceVel = 250.*1e5

    write(*,'(a)') "Reading Donati maps..."

    write(*,'(a)') "Density..."
    select case (misc)
       case("norm")
          open(20, file = "rho35_closed.map", status = "old", form="formatted")
       case("low")
          open(20, file = "rho35_m-0.7_closed.map", status = "old", form="formatted")
       case DEFAULT
          write(*,'(a)') "! unrecognised misc"
    end select
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    allocate(z(1:nz))
    allocate(y(1:ny))
    allocate(rho(ny, nz))
    allocate(t(ny, nz))
    allocate(v(ny, nz))
    allocate(phi(ny, nz))
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((rho(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Velocity..."
    select case (misc)
       case("norm")
          open(20, file = "vel35_closed.map", status = "old", form="formatted")
       case("low")
          open(20, file = "vel35_m-0.7_closed.map", status = "old", form="formatted")
       case DEFAULT
          write(*,'(a)') "! unrecognised misc"
    end select
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((v(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Phi..."
    select case (misc)
       case("norm")
          open(20, file = "phi35_closed.map", status = "old", form="formatted")
       case("low")
          open(20, file = "phi35_m-0.7_closed.map", status = "old", form="formatted")
       case DEFAULT
          write(*,'(a)') "! unrecognised misc"
    end select

    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((phi(i,j),i=1,nz),j=1,ny)
    close(20)

    write(*,'(a)') "Temperature..."
    select case (misc)
       case("norm")
          open(20, file = "temp35_closed.map", status = "old", form="formatted")
       case("low")
          open(20, file = "temp35_m-0.7_closed.map", status = "old", form="formatted")
       case DEFAULT
          write(*,'(a)') "! unrecognised misc"
    end select
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,'(a)') junkline
    read(20,*) nz, ny
    read(20,*) (z(i),i=1,nz)
    read(20,*) (y(i),i=1,ny)
    read(20,*) ((t(i,j),i=1,nz),j=1,ny)
    close(20)
    write(*,'(a)') "Done."

    z = 10.**z
    v = 10.**v
    rho = 10.**rho
    t = 10.**t




    y(1) = 0.
    z(1) = 0.


    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1)-1.
    enddo
    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1)-1.
    enddo
    grid%xAxis = grid%xAxis * 3.
    grid%yAxis = grid%yAxis * 3.

    i1 = grid%nz/2 + 1



    grid%zAxis(i1) = 0.
    grid%zAxis(i1+1) = 0.001
    fac = 1.1
    do i = i1+2,grid%nz
       grid%zAxis(i) = grid%zAxis(i-1) * fac
    enddo


    do i = 1,i1-1
       grid%zAxis(i) = -grid%zAxis(grid%nz-i+1)
    enddo

    grid%zAxis = 3.*grid%zAxis/ abs(grid%zAxis(1))

    do i = 1, grid%nz
       write(*,*) i, grid%zAxis(i)
    enddo
       



    grid%rho = 1.e-30
    grid%temperature = 1.e-30

    grid%inUse = .false.

    write(*,'(a)') "Re-mapping grids..."
    do i1 = 1, grid%nx
       do i2 = 1, grid%ny
          do i3 = 1, grid%nz
             r = sqrt(grid%xAxis(i1)**2 + grid%yAxis(i2)**2 + grid%zAxis(i3)**2)
             if ((r > xDisk(1)).and.(r < 3.)) then
                grid%inUse(i1,i2,i3) = .true.
                r = sqrt(grid%xAxis(i1)**2 + grid%yAxis(i2)**2)
                call hunt(y, ny, r, j)
                call hunt(z, nz, abs(grid%zAxis(i3)), i)
                grid%rho(i1, i2, i3) = rho(i,j)
                grid%temperature(i1, i2, i3) = t(i,j)
                ang = atan2(grid%yAxis(i2),grid%xAxis(i1))
                rHat%z = cos(phi(i,j))
                if (grid%zAxis(i3) < 0.) rHat%z = -rHat%z
                fac = sqrt(1.d0 - rHat%z**2)
                rHat%x = fac * cos(ang)
                rHat%y = fac * sin(ang)
                grid%velocity(i1, i2, i3) = (v(i,j)/cSpeed) * rHat
             endif
          enddo
       enddo
    enddo

    where((grid%temperature < 1.).or.(grid%rho < 1.e-19))
       grid%rho = 1.e-30
       grid%inUse = .false.
    endwhere
    write(*,'(a)') "Done."

 
    write(*,'(a)') "Adding disk..."
    do i1 = 1, grid%nx
       do i2 = 1, grid%ny
          r = sqrt(grid%xAxis(i1)**2 + grid%yAxis(i2)**2)
          if ((r >= xDisk(1)).and.(r <= xDisk(nDisk))) then

             select case(misc)
                case("norm")
                   call locate(xDisk, ndisk, r, j)
                   thickness = yDisk(j) + (r - xDisk(j))*(yDisk(j+1)-yDisk(j))/ &
                        (xDisk(j+1)-xDisk(j))
                   diskDensity = rhoDisk(j) + (r - xDisk(j))*(rhoDisk(j+1)-rhoDisk(j))/ &
                        (xDisk(j+1)-xDisk(j))
                case("low")
                   call locate(xDisk, ndisk, r, j)
                   thickness = yDisk2(j) + (r - xDisk(j))*(yDisk2(j+1)-yDisk2(j))/ &
                        (xDisk(j+1)-xDisk(j))
                   diskDensity = rhoDisk2(j) + (r - xDisk(j))*(rhoDisk2(j+1)-rhoDisk2(j))/ &
                        (xDisk(j+1)-xDisk(j))
             end select

             rVec = VECTOR(grid%xAxis(i1),grid%yAxis(i2),-thickness/2.)
             call getIndices(grid, rVec, i, j, j1, t1, t2, t3)
             rVec = VECTOR(grid%xAxis(i1),grid%yAxis(i2), thickness/2.)
             call getIndices(grid, rVec, i, j, j2, t1, t2, t3)
             do j = j1, j2
                grid%inUse(i1,i2,j) = .true.
                grid%rho(i1,i2,j) = (diskDensity*1.e12) * mHydrogen 
                grid%temperature(i1,i2,j) = 40000.
                zPrime = grid%zAxis(j)
                if (zPrime > 0.) then
                   vHat = VECTOR(0.,0.,-1.)
                else
                   vHat = VECTOR(0.,0.,+1.)
                endif
                vel = abs(zPrime) * (2.e3 * 1.e5) ! 2000km/s per rstar
                grid%velocity(i1,i2,j) = (vel/cSpeed) * vHat
             enddo
          endif
       enddo
    enddo
    
    write(*,'(a)') "Done."


    where (grid%temperature > 1.e6)
       grid%temperature = 1.e6
    end where


    grid%xAxis = grid%xAxis * rStar / 1.e10
    grid%yAxis = grid%yAxis * rStar / 1.e10
    grid%zAxis = grid%zAxis * rStar / 1.e10



    grid%rCore = rStar/1.e10

    write(*,'(a,1p,2e12.5)') "Maximum density: ",&
         MAXVAL(grid%rho, mask=grid%inUse)/mHydrogen,MAXVAL(grid%rho, mask=grid%inUse)
    write(*,'(a,1p,2e12.5)') "Minimum density: ", &
         MINVAL(grid%rho, mask=grid%inUse)/mHydrogen,MINVAL(grid%rho, mask=grid%inUse)
    write(*,'(a,1p,2e12.5)') "Maximum/min temp: ",&
         MAXVAL(grid%temperature, mask=grid%inUse),MINVAL(grid%temperature, mask=grid%inUse)

    grid%rstar1 = grid%rCore
    grid%rstar2 = 0.
    grid%starpos1 = VECTOR(0.,0.,0.)
    grid%starpos2 = VECTOR(1.e20,0.,0.)


    grid%etaLine = 1.e-30
    grid%etaCont = 1.e-30
    grid%chiLine = 1.e-30
    grid%kappaSca = 1.e-30
    grid%kappaAbs = 1.e-30

    write(*,'(a)') "Done."


    if (resonanceLine) then

       grid%etaLine = 1.e-30
       grid%etaCont = 1.e-30
       grid%chiLine = 1.e-30
       grid%kappaSca = 1.e-30
       grid%kappaAbs = 1.e-30
       grid%resonanceLine = .true.
       grid%lambda2 = 1548.202

       grid%chiLine = grid%rho*1.e26
       write(*,'(a)') "This is a resonance line model."

    endif



  end subroutine fillGridDonati2



  subroutine fillGridWR104(grid, rho, teff)
    
    implicit none
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, tubeVec, offset

    integer :: i, j, k
    real :: t1, t2, t3
    integer :: i1, i2, i3
    real :: gridSize
    real :: tubeRadius
    real :: outflowSpeed = 1220. * 10e5
    real :: temperature
    real :: rho
    real :: thisDist, thisTime
    real :: theta, phi
    integer , parameter :: nAng = 360
    integer , parameter :: nRad = 20
    real :: r
    real :: alpha
    real :: teff
    real :: kfac

    write(*,'(a)') "Filling with WR104 model..."
    gridSize = 2000. * AUtoCM / 1.e10
    tubeRadius = 100. * AUtoCM 
    grid%inUse = .false.
    grid%temperature = 10.
    grid%etaCont = 0.
    temperature = 1000.
    do i = 1, grid%nx
       grid%xAxis(i) = 2.*real(i-1)/real(grid%nx-1) - 1.
    enddo

    do i = 1, grid%ny
       grid%yAxis(i) = 2.*real(i-1)/real(grid%ny-1) - 1.
    enddo

    do i = 1, grid%nz
       grid%zAxis(i) = 2.*real(i-1)/real(grid%nz-1) - 1.
    enddo

    grid%xAxis = grid%xAxis * gridSize
    grid%yAxis = grid%yAxis * gridSize
    grid%zAxis = grid%zAxis * gridSize

    offSet = VECTOR(2.*tubeRadius,0.,0.)
    do i = 1, nAng
       theta = real(i-1)/real(nAng-1) * twoPi * 2.
       thisTime = real(i-1)/real(nAng-1) * 220. * 24.* 60. * 60. * 2.
       thisDist = thisTime * outflowSpeed 
       kFac = 220. * 24.* 60. * 60. * 2. * outflowSpeed / fourPi
       do j = 1, nAng
          phi = real(j-1)/real(nAng-1) * twoPi
          do k = 1, nRad
             r = tubeRadius * real(k-1)/real(nRad-1)
             alpha = atan2(kFac,thisDist)
             rVec = VECTOR(thisDist, 0., 0.)
             rVec = rotateZ(rVec, theta)
             tubeVec = VECTOR(r * cos(phi), 0., r*sin(phi))
             tubeVec = rotateZ(tubeVec, theta)
             tubeVec = rotateZ(tubeVec, -alpha)
             rVec = rVec + tubeVec 
             rVec = rVec + offSet
             rVec = rVec / 1.e10
             
             if (.not.outsideGrid(rVec, grid)) then
                call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
                grid%rho(i1,i2,i3) = rho
                grid%inUse(i1,i2,i3) = .true.
                grid%temperature(i1,i2,i3) = temperature
             endif
          enddo
       enddo
    enddo
!    rVec=VECTOR(0., 0., 0.)
!    call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
!    grid%rho(i1,i2,i3) = 1.e-30
!    grid%temperature(i1,i2,i3) = 1.
!    grid%inUse(i1,i2,i3) = .false.


    grid%rCore = 20.*rSol / 1.e10
    grid%lCore = fourPi * stefanBoltz * (20.*rSol)**2 * teff**4

    write(*,'(a)') "done."

  end subroutine fillGridWR104



  subroutine fillGridWind(grid, mDot, rStar, tEff, v0, vterm, beta, &
       lte, contFluxFile, writePops, readPops, popFilename, nLower, nUpper)

    implicit none
    real :: mDot, rStar, tEff, v0, vTerm, beta
    integer :: i, j, k
    character(len=*) :: popFilename, contFluxFile
    logical :: readPops, writePops, lte
    integer :: nLower, nUpper
    real :: sinTheta
    logical :: ok
    type(VECTOR) :: rHat, rVec
    real :: vel
    integer, parameter :: maxLevels = 9
    real :: x, dx
    real(kind=doubleKind) :: phiT, ne1, ne2, ntot
    integer ::m
    real :: v, b2, b3
    type(GRIDTYPE) :: grid

    grid%geometry = "wind"

    do i = 1, grid%nmu
       grid%muAxis(i) = 2.e0*real(i-1)/real(grid%nmu-1)-1.
    enddo

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nphi-1)
    enddo

    x = 1.
    dx = 1.e-3
    do i = 1, grid%nr
       grid%rAxis(i) = x * rStar / 1.e10
       write(*,*) i, x
       x = x + dx
       dx = dx * 1.2
    enddo

    
    grid%rCore = grid%rAxis(1)
    grid%rStar1 = grid%rCore
    grid%rStar2 = 0.
    grid%starPos1 = VECTOR(0.,0.,0.)
    grid%starPos2 = VECTOR(1.e30,0.,0.)
    grid%temperature = 0.9 * tEff
    grid%lineEmission = .true.
    
    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             sinTheta = sqrt(1. - grid%muAxis(j)**2)
             rVec = VECTOR(grid%rAxis(i)*cos(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*sin(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*grid%muAxis(j))
             rHat = rVec
             call normalize(rHat)
             vel = v0 + (vTerm - v0) * (1. - grid%rAxis(1)/grid%rAxis(i))**beta
             grid%velocity(i,j,k) = (vel / cSpeed) * rHat
             grid%rho(i,j,k) = mDot / (fourPi * vel * grid%rAxis(i)**2 * 1.e20)
          enddo
       enddo
    enddo

!    grid%rho = grid%rho * 0.375
!    write(*,'(a)') "! Assuming N(H)/N(He) = 1.5"
    grid%etaCont = 1.e-30

  end subroutine fillGridWind



  
  subroutine writeAxes(grid)
    type(GRIDTYPE) :: grid
    integer :: i, j, k


    if (.not.grid%cartesian) then
       write(*,*) "Radial grid (solar radii)"
       write(*,*) "-------------------------"
       write(*,*) " "
       do i = 1, grid%nr
          write(*,'(i4,f8.1)') i, grid%rAxis(i)/rSol
       enddo
       write(*,*) " "


       write(*,*) "cos(theta) grid (deg)"
       write(*,*) "---------------------"
       write(*,*) " "
       do i = 1, grid%nmu
          write(*,'(i4,f9.3)') i, acos(grid%muAxis(i))*radToDeg
       enddo
       write(*,*) " "
    endif
  end subroutine writeAxes


  real function integratedDensity(grid)
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    real :: tot, dV, dTheta, dPhi, sinTheta, r, phi, dr
    real :: dx, dy ,dz
    type(VECTOR) :: rVec
    real :: fac, mu
    integer :: i1, i2, i3
    real :: t1, t2, t3


    if (grid%cartesian) then
       tot = 0.
       do i = 1, grid%nx-1
          do j = 1, grid%ny-1
             do k = 1, grid%nz-1
                dx = (grid%xAxis(i+1)-grid%xAxis(i))/1.e10 
                dy = (grid%yAxis(j+1)-grid%yAxis(j))/1.e10 
                dz = (grid%zAxis(k+1)-grid%zAxis(k))/1.e10 
                dv = dx * dy * dz
                if (grid%inUse(i,j,k)) tot = tot + grid%rho(i,j,k) * dV
             enddo
          enddo
       enddo
       tot = tot * 1.e30
    else
       tot = 0.
       do i = 2,grid%nr
          do j = 1, grid%nMu-1
             do k = 1, grid%nPhi-1
                dTheta = acos(grid%muAxis(j+1))-acos(grid%muAxis(j))
                mu = 0.5*(grid%muAxis(j+1) + grid%muAxis(j))
                dPhi = grid%phiAxis(k+1) - grid%phiAxis(k)
                sinTheta = sqrt(1.d0 - mu**2)
                r = 0.5*(grid%rAxis(i) + grid%rAxis(i-1))
                phi = 0.5*(grid%phiAxis(k+1) + grid%phiAxis(k))
                dr = grid%rAxis(i) - grid%rAxis(i-1)
                dV = abs(sinTheta * dr * dTheta * dPhi)*r**2
                rVec = VECTOR(r*sinTheta*cos(phi),r*sinTheta*sin(phi),r*mu)
                call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
                fac = grid%rho(i1,i2,i3) !interpGridScalar3(grid%rho,grid%nr,grid%nmu,grid%nphi,i1, i2, i3, t1, t2, t3)
                tot = tot +  fac * dV
             enddo
          enddo
       enddo
    endif
    integratedDensity = tot
  end function integratedDensity


  subroutine fillGridResonance(grid, rCore, mDot, vTerm, beta, temp)

    type(GRIDTYPE) :: grid
    real :: rCore, mDot, vTerm, beta, temp
    integer :: i,j,k
    real :: sinTheta, vel, tot, v0, fac
    real :: x, dx
    type(VECTOR) :: rHat

    grid%lineEmission = .true.

    grid%resonanceLine = .true.


    x = 1.
    dx = 1.e-4
    do i = 1, grid%nr
       grid%rAxis(i) = x * rCore
       x = x + dx
       dx = dx * 1.25
    enddo


    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nPhi-1)
    enddo
    do i = 1, grid%nMu
       grid%muAxis(i) = 2.*real(i-1)/real(grid%nMu-1)-1.
    enddo

    grid%rCore = rCore

    v0 = 0.001 * vTerm

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             sinTheta = sqrt(1.-grid%muAxis(j)**2)
             rHat = VECTOR(sinTheta*cos(grid%phiAxis(k)), sinTheta*sin(grid%phiAxis(k)), grid%muAxis(j))
             vel = v0 + (vterm-v0)*(1.-grid%rCore / grid%rAxis(i))**beta
             grid%rho(i,j,k) = (mdot / (fourPi * grid%rAxis(i)**2 * vel))/mHydrogen
             grid%velocity(i,j,k) = (vel/cSpeed) * rHat
             grid%chiLine(i,j,k) = grid%rho(i,j,k)*10.
             grid%kappaSca(i,j,k,1) = 1.e-30
          enddo
       enddo
    enddo
    grid%temperature = temp

    grid%rAxis(1:grid%nr) =  grid%rAxis(1:grid%nr) / 1.e10
    grid%rCore = grid%rCore / 1.e10
    grid%etaLine = 1.e-30
    grid%kappaAbs = 1.e-30
    grid%etaCont = 1.e-30


  end subroutine fillGridResonance

  subroutine fillmuAxisDisk(muAxis, nMu, openingAngle)
    integer :: nMu
    real :: openingAngle
    real :: muAxis(*)
    integer :: i, n, n1, n2
    real :: cosTheta1, cosTheta2, cosThetaStart, cosThetaEnd
    
    cosThetaStart = -1.
    cosThetaEnd = 1.
    
    cosTheta1 = cos(pi/2+openingAngle)
    cosTheta2 = cos(pi/2.-openingAngle)

    n1 = nMu / 4
    n2 = 3 * nMu / 4
    do i = 1, n1
       muAxis(i) = cosThetaStart + (cosTheta1 - cosThetaStart) * real(i-1)/real(n1)
    enddo
    do i = n1+1,n2
       muAxis(i) = cosTheta1 + (cosTheta2 - cosTheta1) * real(i - n1 -1)/real(n2 - n1)
    enddo
    do i = n2 + 1, nMu
       muAxis(i) = cosTheta2 + (cosThetaEnd - cosTheta2) * real(i-n2-1)/real(nMu-n2-1)
    enddo
  end subroutine fillmuAxisDisk
end module grid_mod


