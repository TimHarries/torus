! this module contains the description of the 3d grid of opacity arrays
! and contains subroutines to initialize the grid with a variety of 
! geometries.

! written by tjh

! v1.0 on 13/08/99

! v1.1 on 08/01/00
!   added colliding winds and bipolar geometries

! v1.2 on 19/05/00
! added disk stuff in LTE

! extended to allow adaptive mesh refinement grid. nhs


module grid_mod

  use gridtype_mod                    ! type definition for the 3-d grid
  use kind_mod
  use constants_mod                   ! physical constants
  use vector_mod                      ! vector math
  use atom_mod                        ! LTE atomic physics
  use utils_mod
  use octal_mod                       ! octal type for amr
  use amr_mod
  use density_mod                     ! to use generic density function
  use pixplot_module                  ! To use some plotting routines.
  use cluster_class
  use cmfgen_class
  use messages_mod
  use ion_mod

  implicit none

  public

  private :: writeReal1D, writeReal2D, writeDouble2D
  private :: readReal1D,  readReal2D,  readDouble2D


  interface getIndices
     module procedure getIndices_single
     module procedure getIndices_octal
  end interface

  interface outsideGrid
     module procedure outsideGrid_single
     module procedure outsideGrid_double
     module procedure outsideGrid_octal
  end interface
     

contains

  ! function to initialize a cartesian grid

  function initCartesianGrid(nx, ny, nz, nLambda, lamStart, lamEnd, &
       flatspec, ok)

    type(GRIDTYPE) :: initCartesianGrid   ! the grid
    integer :: nx, ny, nz                 ! x,y,z sizes
    integer :: nLambda                    ! no of wavelength points
    integer :: ierr                       ! allocation error status
    integer :: ilambda                   ! counters
    real    :: lamStart, lamEnd           ! start and end wavelengths
    logical :: ok                         ! function done ok?
    logical :: flatspec                   ! is the spectrum flat

    iLambda = nLambda

    initCartesianGrid%oneKappa = .false.
    initCartesianGrid%lineEmission = .false. ! the default

    ! if the spectrum is flat one only needs on wavelength point

    initCartesianGrid%flatspec = flatspec

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
    initCartesianGrid%polar = .false.
    initCartesianGrid%adaptive = .false.

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
    integer :: ilambda
    real    :: lamStart, lamEnd
    integer :: i, j, k
    logical :: ok, flatspec

    initPolarGrid%lineEmission = .false. ! the default

    iLambda = nLambda

    initPolarGrid%flatspec = flatspec
    initPolarGrid%oneKappa = .false.

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

    allocate(initPolarGrid%dvbydr(1:nr,1:nmu,1:nphi))


    allocate(initPolarGrid%rAxis(1:nr))
    allocate(initPolarGrid%muAxis(1:nmu))
    allocate(initPolarGrid%phiAxis(1:nphi))
    allocate(initPolarGrid%lamArray(1:nlambda))

    allocate(initPolarGrid%biasLine(1:nr))
    allocate(initPolarGrid%biasCont(1:nr))


    allocate(initPolarGrid%oneProbCont(1:initPolarGrid%nProb))
    allocate(initPolarGrid%oneProbLine(1:initPolarGrid%nProb))
    k = nr*nMu*nPhi
    allocate(initPolarGrid%cellIndex(1:k,1:3))
    
    initPolargrid%nProb = 0
    do i = 1, nr
       do j = 1, nMu
          do k = 1, nPhi
             initPolargrid%nProb = initPolargrid%nProb + 1
             initPolargrid%cellIndex(initPolargrid%nProb,1) = i
             initPolargrid%cellIndex(initPolargrid%nProb,2) = j
             initPolargrid%cellIndex(initPolargrid%nProb,3) = k
          enddo
       enddo
    enddo
    allocate(initPolarGrid%oneProbCont(1:initPolargrid%nProb))
    allocate(initPolarGrid%oneProbLine(1:initPolargrid%nProb))


    initPolarGrid%biasLine = 1.
    initPolarGrid%biasCont = 1.


    initPolarGrid%nr = nr
    initPolarGrid%nmu = nmu
    initPolarGrid%nphi = nphi

    initPolarGrid%na1 = nr
    initPolarGrid%na2 = nmu
    initPolarGrid%na3 = nphi

    allocate(initPolarGrid%biasLine3D(1:nr, 1:nmu, 1:nphi), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    allocate(initPolarGrid%biasCont3D(1:nr, 1:nmu, 1:nphi), stat=ierr)
    if (ierr /=0) then
       write(*,'(a)') "! Cannot allocate grid memory"
       ok = .false.
       goto 666
    endif

    initPolarGrid%biasLine3D = 1.
    initPolarGrid%biasCont3D = 1.
    


    initPolarGrid%nlambda = nlambda

    initPolarGrid%rho = 0.
    initPolarGrid%kappaSca = 0.
    initPolarGrid%kappaAbs = 0.
    initPolarGrid%cartesian = .false.
    initPolarGrid%polar = .true.
    initPolarGrid%adaptive = .false.

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
  ! function to initialize a cartesian grid


  subroutine initAMRgrid(greyContinuum,newContFile,flatspec,grid,ok,theta1,theta2)

    use input_variables

    implicit none

    ! grid%timeNow must be assigned before this routine is called!

    logical, intent(in) :: greyContinuum
    character(len=80), intent(out) :: newContFile
    logical, intent(in) :: flatspec        ! is the spectrum flat
    logical, intent(inout) :: ok           ! function done ok?
    type(GRIDTYPE), intent(out) :: grid                 ! the grid
    real, intent(out)   :: theta1, theta2
    
    integer :: i,ilambda                   ! counters
    real :: rStar
    
    ! ok for now
    ok = .true.

    iLambda = nLambda
    grid%oneKappa = oneKappa

    newContFile = " "
    
    if (oneKappa) then
       allocate(grid%oneKappaAbs(nDustType,1:nLambda))
       allocate(grid%oneKappaSca(nDustType,1:nLambda))
       grid%nTempRossArray = 100
       allocate(grid%kappaRossArray(nDustType,1:grid%nTempRossArray))
       allocate(grid%tempRossArray(1:grid%nTempRossArray))
    endif

    grid%lineEmission = lineEmission
    grid%maxLevels = statEqMaxLevels

    ! if the spectrum is flat one only needs on wavelength point

    grid%flatspec = flatspec
    grid%statEq2d = statEq2d
    if (flatspec) ilambda = 1

    grid%amr2dOnly = amr2dOnly
    grid%geometry = geometry

    select case (geometry)
    
    case("toruslogo")
       grid%geometry = "toruslogo"
       grid%rStar1 = rSol/1.e10

    case("whitney")
       grid%geometry = "whitney"
       grid%rCore = rStellar/1.e10
       grid%lCore = fourPi * rCore**2 * stefanBoltz * teff**4 * 1.e20
       grid%rInner = erInner/1.e10
       grid%rOuter = erOuter/1.e10
       grid%starPos1 = vector(1.e-6,1.e-6,1.e-6)


    case("benchmark")
       grid%geometry = "benchmark"
       grid%rCore = rCore
       grid%lCore = fourPi * rCore**2 * stefanBoltz * teff**4 * 1.e20
       grid%rInner = rInner
       grid%rOuter = rOuter

    case("molebench")
       grid%geometry = "molebench"

    case("shakara")
       grid%geometry = "shakara"
       grid%rCore = rCore
       grid%lCore = fourPi * rCore**2 * stefanBoltz * teff**4 * 1.e20
       grid%rInner = rInner
       grid%rOuter = rOuter

    case("warpeddisc")
       grid%geometry = "warpeddisc"
       grid%rCore = rCore
       grid%lCore = fourPi * rCore**2 * stefanBoltz * teff**4 * 1.e20
       grid%rInner = rInner
       grid%rOuter = rOuter

    case("clumpydisc")
       grid%geometry = "clumpydisc"
       grid%rCore = rCore
       grid%lCore = fourPi * rCore**2 * stefanBoltz * teff**4 * 1.e20
       grid%rInner = rInner
       grid%rOuter = rOuter

    case("aksco")
       grid%geometry = "aksco"
       grid%rCore = rCore
       grid%lCore = fourPi * rCore**2 * stefanBoltz * teff**4 * 1.e20
       grid%rInner = rInner
       grid%rOuter = rOuter

    case("symbiotic")
       grid%geometry = "symbiotic"
       grid%oneKappa = .true.
       oneKappa = .true.

    case("lexington","fractal")
       grid%rCore = 18.67 * rSol / 1.e10
       grid%rInner = rinner
       grid%rOuter = 2.e09
       grid%oneKappa = .true.
       oneKappa = .true.

    case ("ttauri") 
       call initTTauriAMR(grid,theta1,theta2)
    case ("windtest") 
       call initWindTestAMR(grid)

    case ("jets") 
       call initJetsAMR(grid)

    case ("luc_cir3d") 
       rStar  = CIR_Rstar*Rsol/1.0d10   ! in [10^10cm] 
       grid%rCore = rStar
       grid%rStar1 = rStar
       grid%starPos1 = vector(0.,0.,0.)

    case ("cmfgen") 
       rStar  = get_cmfgen_data_array_element("R", 1)   ! in [10^10cm] 
       grid%rCore = rStar
       grid%rStar1 = rStar
       grid%starPos1 = vector(0.,0.,0.)

    case ("romanova") 
       rStar  = ROM_Rs*Rsol/1.0d10   ! in [10^10cm] 
       grid%rCore = rStar
       grid%rStar1 = rStar
       grid%rStar2 = 0.
       grid%starPos1 = vector(0.,0.,0.)
       grid%diskRadius = ThinDiskRin  ! [10^10cm]
       grid%diskNormal = VECTOR(0.,0.,1.)

    case ("testamr")
       call initTestAMR(grid)

    case("starburst")
       grid%rCore = 18.67 * rSol / 1.e10
       grid%rInner = rinner
       grid%rOuter = 2.e09
       grid%oneKappa = .true.
       oneKappa = .true.
       grid%geometry = "starburst"
       grid%lineEmission = .false.


    case ("proto")
       call initProtoAMR(grid)

    case("cluster","wr104")
       call initClusterAMR(grid)
       
    case("spiralwind")
       grid%rCore = rCore / 1.e10

    case("ppdisk")
       grid%geometry = "ppdisk"
       grid%lineEmission = .false.
       grid%rInner = rInner
       grid%rOuter = rOuter
       grid%rCore = rCore

    case("planetgap")
       grid%geometry = "planetgap"
       grid%lineEmission = .false.
       grid%rInner = rInner
       grid%rOuter = rOuter
       grid%rCore = rCore

    case("wrshell")
       grid%geometry = "wrshell"
       grid%lineEmission = .false.
       grid%rInner = rInner
       grid%rOuter = rOuter
       grid%rCore = rCore


    case DEFAULT
       print *, '!!!WARNING: The ''',geometry,''' geometry may not yet have been implemented'
       print *, '            for use with an adaptive grid.'
    end select
    


    nullify(grid%rho)
    nullify(grid%kappaAbs)
    nullify(grid%kappaSca)
    nullify(grid%velocity)
    nullify(grid%temperature)
    nullify(grid%chiLine)
    nullify(grid%etaLine)
    nullify(grid%etaCont)
    nullify(grid%sigma)
    nullify(grid%biasLine3D)
    nullify(grid%biasCont3D)
    nullify(grid%xAxis)
    nullify(grid%yAxis)
    nullify(grid%zAxis)
    nullify(grid%xProbDistLine)
    nullify(grid%yProbDistLine)
    nullify(grid%zProbDistLine)
    nullify(grid%xProbDistCont)
    nullify(grid%yProbDistCont)
    nullify(grid%zProbDistCont)
    nullify(grid%inStar)
    nullify(grid%inUse)
    nullify(grid%rProbDistLine)
    nullify(grid%muProbDistLine)
    nullify(grid%phiProbDistLine)
    nullify(grid%rProbDistCont)
    nullify(grid%muProbDistCont)
    nullify(grid%phiProbDistCont)
    nullify(grid%rAxis)
    nullify(grid%muAxis)
    nullify(grid%phiAxis)


    allocate(grid%lamArray(1:nlambda))

    grid%nLambda = nLambda


    grid%adaptive = .true.
    grid%cartesian = .false.
    grid%polar = .false.

    grid%resonanceLine = .false.
    grid%smoothingFactor = smoothFactor

  end subroutine initAMRgrid


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

  subroutine fillGridEllipse(grid, rho, rMin, rMaj, rInner, teff)


    implicit none
    integer, parameter :: nRmaj = 200
    type(GRIDTYPE) :: grid
    real :: x,y,z,r,r1
    real :: phi
    real :: rho
    real :: teff
    real :: rMaj, rMin, rInner
    integer, parameter :: nPhi = 500
    integer :: i, j, k
    integer :: i1, i2, i3, i4

    grid%geometry = "ellipse"

    grid%rCore = rSol / 1.e10
    grid%lCore = fourPi * stefanBoltz * grid%rCore**2 * 1.e20 * teff**4
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

    grid%inUse = .false.

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
             r1 = sqrt(x**2 + y**2 + grid%zAxis(k)**2)
             if (r1 > rInner) then
                Grid%rho(i1,i2,k) = rho * (rInner / r1)**2
                grid%inUse(i1,i2,k) = .true.
             endif
          enddo
       enddo
    enddo

    open(22,file="temps.dat",status="old",form="unformatted",err=666)
    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             read(22) grid%temperature(i,j,k),grid%etaCont(i,j,k)
          enddo
       enddo
    enddo
    close(22)

666 continue

  end subroutine fillGridEllipse


  subroutine fillGridDisk(grid, rhoNought, rCore, rInner, rOuter, height, &
       mCore, diskTemp)


    implicit none
    type(GRIDTYPE) :: grid
    integer, parameter :: nRad = 100
    real :: rOuter, rInner, rCore
    type(VECTOR) :: rVec, vVec
    real :: r 
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
                grid%kappaSca(i,j,k,1) = max(1.e-30,grid%rho(i,j,k) * real(sigmaE))
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
    real :: r
    real :: rho, rhoNought
    real :: diskTemp
    real :: height
    real :: vel, mCore
    type(VECTOR) :: spinAxis
    integer :: i, j, k 

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
                grid%kappaSca(i,j,k,1) = max(1.e-30,grid%rho(i,j,k) * real(sigmaE))
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
    grid%temperature = 100.
    grid%lCore = fourPi * stefanBoltz * (rsol*rsol) * 4000.**4
    grid%rCore = 2.* rsol / 1.e10

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
    real(double) :: vElement
    real(double) :: totalDustMass
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
    real, dimension(:,:), intent(in) :: contTau
    type(VECTOR) :: viewVec
    type(VECTOR) :: coolStarPosition
    integer :: i, j, npix
    real, allocatable :: plane(:,:), axis(:)
!    real(double), allocatable :: plane_dble(:,:)
    logical, allocatable :: inuse(:,:)
    real :: tr(6)
    real :: fg, bg, dx, dz
    integer :: resolution
    real :: sampleFreq
    integer :: i1, i2, i3
    real :: r, mu, phi
    character(len=*) :: device
    real(oct) :: x, y, z
    real :: smallOffset
    type(octalVector) :: startPoint
   integer :: pgbegin

    if (grid%adaptive) then
       
       resolution = 250 
       allocate(plane(resolution,resolution))
!       allocate(plane_dble(resolution,resolution))

       ! we will try to slightly offset our sampling points w.r.t.
       !   the subcell walls to reduce numerical problems.
       y = grid%octreeRoot%centre%y - 0.995 * grid%octreeRoot%subcellSize 
       
       tr(1) = grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize
       tr(2) = 0.995 * grid%octreeRoot%subcellSize * 2.0 / real(resolution)
       tr(3) = 0.
       tr(4) = grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize
       tr(5) = 0.
       tr(6) = 0.995 * grid%octreeRoot%subcellSize * 2.0 / real(resolution)
       do i = 1, resolution
          do j = 1, resolution
             x = tr(1) + tr(2)*i
             z = tr(1) + tr(2)*j
             startPoint = octalVector(x,y,z)
             plane(i,j) = columnDensity(grid,startPoint,yHatOctal,1.0)
          enddo
       enddo
       bg = 1.e30
       fg = -1.e30 
       do i = 1, resolution
          do j = 1, resolution
             if (plane(i,j) > 0.) then
                plane(i,j) = log10(plane(i,j))
                if (plane(i,j) > fg) fg = plane(i,j)
                if (plane(i,j) < bg) bg = plane(i,j)
             else
                plane(i,j) = -30.
             endif
          enddo
       enddo
       write(*,*) "foreground",fg,"background",bg

       i = pgbegin(0,device,1,1)
       
       call pgvport(0.1,0.9,0.1,0.9)
       
       call pgwnad(REAL(grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize),&
                   REAL(grid%octreeRoot%centre%x + grid%octreeRoot%subcellSize),&
                   REAL(grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize),&
                   REAL(grid%octreeRoot%centre%z + grid%octreeRoot%subcellSize))
       
       call pgimag(plane, resolution, resolution, 1, resolution, 1, resolution, bg, fg, tr)
       
       call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
       call pgend

       deallocate(plane)

    else if (grid%cartesian) then

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
       i =  pgbegin(0,device,1,1)
       
       call pgvport(0.1,0.9,0.1,0.9)
       
       call pgwnad(grid%xAxis(1), grid%xAxis(grid%nx), grid%zAxis(1), grid%zAxis(grid%nz))
       
       call pgimag(plane, grid%nx, grid%nz, 1, grid%nx, 1, grid%nz, fg, bg, tr)
       
       call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
       call pgend
       deallocate(plane)
       
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
       i =  pgbegin(0,device,1,1)
       
       call pgvport(0.1,0.9,0.1,0.9)
       
       call pgwnad(axis(1), axis(npix), axis(1), axis(npix))
       
       call pgimag(plane, npix, npix, 1, npix, 1, npix, fg, bg, tr)
       
       call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
       call pgend



       deallocate(plane)
       deallocate(inuse)
       deallocate(axis)


    endif

    write(*,'(a)') "Finished plotting grid..."

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
    sigmaScaBlue = (34.+6.6)*sigmaE * 1.e10
    sigmaAbsBlue = 1.e-5 * sigmaE * 1.e10

    sigmaScaRed = 1.e-5 * sigmaE * 1.e10
    sigmaAbsRed = 1.e-5 * sigmaE * 1.e10

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
                rho = mdot / (4.* pi * (r*1.e10)**2 * v)
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
                grid%kappaAbs(i,j,k,1:1) = 100.
                grid%kappaSca(i,j,k,1:1) = 100.
                grid%kappaAbsRed(i,j,k,1:1) = 100.
                grid%kappaScaRed(i,j,k,1:1) = 100.
                grid%velocity(i,j,k) = VECTOR(0.,0.,0.)
             endif
          enddo
       enddo
    enddo
!    grid%kappaAbs(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1) = 1.e-30
!    grid%kappaAbsRed(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1) = 1.e-30
!    grid%kappaScaRed(1:grid%nx, 1:grid%ny, 1:grid%nz, 1:1) = 1.e-30
    grid%etaLine(1:grid%nx, 1:grid%ny, 1:grid%nz) = 1.e-30
    grid%chiLine(1:grid%nx, 1:grid%ny, 1:grid%nz) = 1.e-30
    grid%etaCont(1:grid%nx, 1:grid%ny, 1:grid%nz) = 1.e-30
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
             if (i < nr) then
                grid%velocity(i,j,k) = (1.e5*logInterp(v,nr,r,grid%rAxis(i))/cSpeed) * rHat
                grid%temperature(i,j,k) = logInterp(t, nr, r, grid%rAxis(i))*1.e4
                grid%chiLine(i,j,k) = (logInterp(chil, nr, r, grid%rAxis(i)))*fac**2
                grid%etaLine(i,j,k) = (logInterp(etal, nr,r, grid%rAxis(i)))*fac**2
                grid%etaCont(i,j,k) = (logInterp(eta, nr,r, grid%rAxis(i)))*fac**2
                grid%sigma(i, j, k) = sigma(i)
                grid%kappaSca(i,j,k,1) = logInterp(escat, nr, r, grid%rAxis(i))*fac
                tmp = logInterp(chi, nr, r,grid%rAxis(i))-grid%kappaSca(i,j,k,1)
                grid%kappaAbs(i,j,k,1) = max(1.e-20,tmp)
             else
                grid%velocity(i,j,k) = (1.e5*v(i)/cSpeed) * rHat
                grid%temperature(i,j,k) = t(i)*1.e4
                grid%chiLine(i,j,k) = chil(i)*fac**2
                grid%etaLine(i,j,k) = etal(i)*fac**2
                grid%sigma(i, j, k) = sigma(i)
                grid%kappaSca(i,j,k,1) = escat(i)*fac
                tmp = logInterp(chi, nr, r,grid%rAxis(i))-grid%kappaSca(i,j,k,1)
                grid%etaCont(i,j,k) = eta(i)*fac**2
                grid%kappaAbs(i,j,k,1) = max(1.e-20,tmp)
             endif
          enddo
       enddo
    enddo

!    grid%etaLine(1:grid%nr, 1:grid%nMu, 1:grid%nPhi) = grid%etaLine(1:grid%nr, 1:grid%nMu, 1:grid%nPhi) * fourPi

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
    real :: smallOffset
    type(octalVector) :: samplePoint
    type(vector) :: velocityVector
    integer :: pgbegin

    i = pgbegin(0,device,1,1)


    allocate(xArray(1:nx,1:ny))
    allocate(yArray(1:nx,1:ny))
    allocate(usedArray(1:nx,1:ny))

    usedArray = .false.
    xArray = 0.
    yArray = 0.

    if (grid%adaptive) then

       ! we will try to slightly offset our sampling points w.r.t.
       !   the subcell walls to reduce numerical problems.
       smallOffset = grid%halfSmallestSubcell * 0.1
         
       ! also need to make sure we do not try to sample outside the grid.  
       tr(1) = grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize
       tr(2) = 0.995 * grid%octreeRoot%subcellSize * 2.0 / real(nx)
       tr(3) = 0.
       tr(4) = grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize
       tr(5) = 0.
       tr(6) = 0.995 * grid%octreeRoot%subcellSize * 2.0 / real(ny)
       do i = 1, nx
          do j = 1, ny
             x = tr(1) + tr(2)*i
             y = tr(1) + tr(2)*j
             samplePoint = octalVector(x,grid%octreeRoot%centre%y+smallOffset,y)
             velocityVector = amrGridVelocity(grid%octreeRoot,samplePoint)
             xArray(i,j) = velocityVector%x*cSpeed
             yArray(i,j) = velocityVector%z*cSpeed
          enddo
       enddo

       

    elseif (grid%cartesian) then
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

    print *, 'Deallocating grid structure'
    
!    if (associated(grid%octreeRoot)) then
!      call deleteOctreeBranch(grid%octreeRoot,grid)   
!      deallocate(grid%octreeRoot)
!      nullify(grid%octreeRoot)
!    end if
      
    if (associated(grid%oneKappaAbs)) deallocate(grid%oneKappaAbs)
       nullify(grid%oneKappaAbs)
    if (associated(grid%oneKappaSca)) deallocate(grid%oneKappaSca)
       nullify(grid%oneKappaSca)
    if (associated(grid%rho)) deallocate(grid%rho)
       nullify(grid%rho)
    if (associated(grid%kappaAbs)) deallocate(grid%kappaAbs)
       nullify(grid%kappaAbs)
    if (associated(grid%kappaSca)) deallocate(grid%kappaSca)
       nullify(grid%kappaSca)
    if (associated(grid%kappaAbsRed)) deallocate(grid%kappaAbsRed)
       nullify(grid%kappaAbsRed)
    if (associated(grid%kappaScaRed)) deallocate(grid%kappaScaRed)
       nullify(grid%kappaScaRed)
    if (associated(grid%chiLine)) deallocate(grid%chiLine)
       nullify(grid%chiLine)
    if (associated(grid%etaLine)) deallocate(grid%etaLine)
       nullify(grid%etaLine)
    if (associated(grid%etaCont)) deallocate(grid%etaCont)
       nullify(grid%etaCont)
    if (associated(grid%sigma)) deallocate(grid%sigma)
       nullify(grid%sigma)
    if (associated(grid%velocity)) deallocate(grid%velocity)
       nullify(grid%velocity)
    if (associated(grid%rAxis)) deallocate(grid%rAxis)
       nullify(grid%rAxis)
    if (associated(grid%biasCont)) deallocate(grid%biasCont)
       nullify(grid%biasCont)
    if (associated(grid%biasLine)) deallocate(grid%biasLine)
       nullify(grid%biasLine)
    if (associated(grid%muAxis)) deallocate(grid%muAxis)
       nullify(grid%muAxis)
    if (associated(grid%phiAxis)) deallocate(grid%phiAxis)
       nullify(grid%phiAxis)
    if (associated(grid%xAxis)) deallocate(grid%xAxis)
       nullify(grid%xAxis)
    if (associated(grid%yAxis)) deallocate(grid%yAxis)
       nullify(grid%yAxis)
    if (associated(grid%zAxis)) deallocate(grid%zAxis)
       nullify(grid%zAxis)
    if (associated(grid%lamArray)) deallocate(grid%lamArray)
       nullify(grid%lamArray)
    if (associated(grid%rProbDistLine)) deallocate(grid%rProbDistLine)
       nullify(grid%rProbDistLine)
    if (associated(grid%muProbDistLine)) deallocate(grid%muProbDistLine)
       nullify(grid%muProbDistLine)
    if (associated(grid%phiProbDistLine)) deallocate(grid%phiProbDistLine)
       nullify(grid%phiProbDistLine)
    if (associated(grid%xProbDistLine)) deallocate(grid%xProbDistLine)
       nullify(grid%xProbDistLine)
    if (associated(grid%yProbDistLine)) deallocate(grid%yProbDistLine)
       nullify(grid%yProbDistLine)
    if (associated(grid%zProbDistLine)) deallocate(grid%zProbDistLine)
       nullify(grid%zProbDistLine)
    if (associated(grid%rProbDistCont)) deallocate(grid%rProbDistCont)
       nullify(grid%rProbDistCont)
    if (associated(grid%muProbDistCont)) deallocate(grid%muProbDistCont)
       nullify(grid%muProbDistCont)
    if (associated(grid%phiProbDistCont)) deallocate(grid%phiProbDistCont)
       nullify(grid%phiProbDistCont)
    if (associated(grid%xProbDistCont)) deallocate(grid%xProbDistCont)
       nullify(grid%xProbDistCont)
    if (associated(grid%yProbDistCont)) deallocate(grid%yProbDistCont)
       nullify(grid%yProbDistCont)
    if (associated(grid%zProbDistCont)) deallocate(grid%zProbDistCont)
       nullify(grid%zProbDistCont)
    if (associated(grid%temperature)) deallocate(grid%temperature)
       nullify(grid%temperature)
    if (associated(grid%biasLine3D)) deallocate(grid%biasLine3D)
       nullify(grid%biasLine3D)
    if (associated(grid%biasCont3D)) deallocate(grid%biasCont3D)
       nullify(grid%biasCont3D)
    if (associated(grid%N)) deallocate(grid%N)
       nullify(grid%N)
    if (associated(grid%Ne)) deallocate(grid%Ne)
       nullify(grid%Ne)
    if (associated(grid%nTot)) deallocate(grid%nTot)
       nullify(grid%nTot)
    if (associated(grid%inStar)) deallocate(grid%inStar)
       nullify(grid%inStar)
    if (associated(grid%inUse)) deallocate(grid%inUse)
       nullify(grid%inUse)

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

  logical pure function outsideGrid_single(posVec, grid)    
    type(GRIDTYPE), intent(in) :: grid
    type(VECTOR), intent(in)   :: posVec

    outSideGrid_single = .false.

    if ((posVec%x > grid%xAxis(grid%nx)) .or. &
         (posVec%y > grid%yAxis(grid%ny)) .or. &
         (posVec%z > grid%zAxis(grid%nz)) .or. &
         (posVec%x < grid%xAxis(1)) .or. &
         (posVec%y < grid%yAxis(1)) .or. &
         (posVec%z < grid%zAxis(1))) outsideGrid_single = .true.

  end function outsideGrid_single

  logical pure function outsideGrid_double(posVec, grid)
    type(GRIDTYPE), intent(in) :: grid
    type(DOUBLEVECTOR), intent(in)   :: posVec

    outSideGrid_double = .false.

    if ((posVec%x > grid%xAxis(grid%nx)) .or. &
         (posVec%y > grid%yAxis(grid%ny)) .or. &
         (posVec%z > grid%zAxis(grid%nz)) .or. &
         (posVec%x < grid%xAxis(1)) .or. &
         (posVec%y < grid%yAxis(1)) .or. &
         (posVec%z < grid%zAxis(1))) outsideGrid_double = .true.

  end function outsideGrid_double

  logical pure function outsideGrid_octal(posVec, grid)
    type(GRIDTYPE), intent(in) :: grid
    type(OCTALVECTOR), intent(in)   :: posVec

    outSideGrid_octal = .false.
    
    if ((posVec%x > grid%xAxis(grid%nx)) .or. &
         (posVec%y > grid%yAxis(grid%ny)) .or. &
         (posVec%z > grid%zAxis(grid%nz)) .or. &
         (posVec%x < grid%xAxis(1)) .or. &
         (posVec%y < grid%yAxis(1)) .or. &
         (posVec%z < grid%zAxis(1))) outsideGrid_octal = .true.
    
  end function outsideGrid_octal


  subroutine getIndices_single(grid, rVec, i1, i2, i3, t1, t2, t3)
    ! making this PURE may cause problems with XL Fortran
    type(GRIDTYPE), intent(in) :: grid
    type(VECTOR), intent(in)   :: rVec
    integer, intent(out)       :: i1, i2, i3
    real, intent(out)          :: t1, t2 ,t3
    real                       :: r, theta, phi, mu

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

  end subroutine getIndices_single



  subroutine getIndices_octal(grid, rVec, i1, i2, i3, t1, t2, t3)
    ! making this PURE may cause problems with XL Fortran
    type(GRIDTYPE), intent(in)          :: grid
    type(OCTALVECTOR), intent(in)       :: rVec
    integer, intent(out)                :: i1, i2, i3
    real(oct), intent(out)   :: t1, t2 ,t3
    real(oct)                :: r, theta, phi, mu

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


       call hunt(grid%xAxis, grid%nx, real(rvec%x), i1)
       call hunt(grid%yAxis, grid%ny, real(rvec%y), i2)
       call hunt(grid%zAxis, grid%nz, real(rvec%z), i3)
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
       call hunt(grid%rAxis, grid%nr, real(r), i1)
       if (i1 == 0) then
          i1 = 1
       endif
       if (i1 == grid%nr) i1 = grid%nr-1

       t1 = (r-grid%rAxis(i1))/(grid%rAxis(i1+1)-grid%rAxis(i1))


       call hunt(grid%muAxis, grid%nMu, real(mu), i2)
       if (i2 == grid%nMu) i2 = grid%nMu-1
       t2 = (mu-grid%muAxis(i2))/(grid%muAxis(i2+1)-grid%muAxis(i2))

       call hunt(grid%phiAxis, grid%nPhi, real(phi), i3)
       if (i3 == grid%nPhi) i3=i3-1
       t3 = (phi-grid%phiAxis(i3))/(grid%phiAxis(i3+1)-grid%phiAxis(i3))

    endif

  end subroutine getIndices_octal



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

  
  subroutine writeAMRgrid(filename,fileFormatted,grid)
    use input_variables, only : molecular, cmf
    ! writes out the 'grid' for an adaptive mesh geometry  

    implicit none
  
    character(len=*)           :: filename
    logical, intent(in)        :: fileFormatted
    type(GRIDTYPE), intent(in) :: grid
    
    integer, dimension(8) :: timeValues ! system date and time
    integer               :: error      ! error code

    if (fileFormatted) then 
       open(unit=20,iostat=error, file=filename, form="formatted", status="replace")
    else 
       open(unit=20,iostat=error, file=filename, form="unformatted", status="replace")
    end if        
    call writeInfo("Writing AMR grid file to: "//trim(filename),TRIVIAL)
    
    call date_and_time(values=timeValues)
    
    if (fileFormatted) then 
            
       ! write a time stamp to the file
       write(unit=20,fmt=*,iostat=error) timeValues(:)
       if (error /=0) then
         print *, 'Panic: write error in writeAMRgrid (formatted timeValues)' ; stop
       end if
           
       ! write the variables that are stored in the top-level 'grid' structure
       write(unit=20,fmt=*,iostat=error) grid%nLambda, grid%flatSpec, grid%adaptive,& 
               grid%cartesian, grid%isotropic, grid%hitCore, grid%diskRadius, &
               grid%diskNormal, grid%DipoleOffset, grid%geometry, grid%rCore, &
               grid%lCore, grid%chanceWindOverTotalContinuum,                 &
               grid%lineEmission, grid%contEmission, grid%doRaman,            &
               grid%resonanceLine, grid%rStar1, grid%rStar2, grid%lumRatio,   &
               grid%tempSource, grid%starPos1, grid%starPos2, grid%lambda2,   &
               grid%maxLevels, grid%maxDepth, grid%halfSmallestSubcell,       &
               grid%nOctals, grid%smoothingFactor, grid%oneKappa, grid%rInner,&
               grid%rOuter, grid%amr2dOnly, grid%photoionization
       if (error /=0) then
         print *, 'Panic: write error in writeAMRgrid (formatted variables)' ; stop
       end if
               
    else
            
       ! write a time stamp to the file
       write(unit=20,iostat=error) timeValues(:)
       if (error /=0) then
         print *, 'Panic: write error in writeAMRgrid (unformatted timeValues)' ; stop
       end if
    
       ! write the variables that are stored in the top-level 'grid' structure
       write(unit=20,iostat=error) grid%nLambda, grid%flatSpec, grid%adaptive,& 
               grid%cartesian, grid%isotropic, grid%hitCore, grid%diskRadius, &
               grid%diskNormal, grid%DipoleOffset, grid%geometry, grid%rCore, &
               grid%lCore, grid%chanceWindOverTotalContinuum,                 &
               grid%lineEmission, grid%contEmission, grid%doRaman,            &
               grid%resonanceLine, grid%rStar1, grid%rStar2, grid%lumRatio,   &
               grid%tempSource, grid%starPos1, grid%starPos2, grid%lambda2,   &
               grid%maxLevels, grid%maxDepth, grid%halfSmallestSubcell,       &
               grid%nOctals, grid%smoothingFactor, grid%oneKappa, grid%rInner,& 
               grid%rOuter, grid%amr2dOnly, grid%photoionization
       if (error /=0) then
         print *, 'Panic: write error in writeAMRgrid (unformatted variables)' ; stop
       end if
               
    end if 
    
    call writeReal1D(grid%lamarray,fileFormatted)
    call writeReal2D(grid%oneKappaAbs,fileFormatted)
    call writeReal2D(grid%oneKappaSca,fileFormatted)
    call writeClumps(fileFormatted)

    ! now we call the recursive subroutine to store the tree structure 
    if (associated(grid%octreeRoot)) then
       if (fileFormatted) then 
          write(unit=20,fmt=*) .true.
       else
          write(unit=20) .true.
       end if
       call writeOctreePrivate(grid%octreeRoot,fileFormatted, grid)
    else 
       if (fileFormatted) then 
          write(unit=20,fmt=*) .false.
       else
          write(unit=20) .false.
       end if
    end if

    endfile 20
    close(unit=20)
    
  contains
  
    recursive subroutine writeOctreePrivate(thisOctal,fileFormatted, grid)
       ! writes out an octal from the grid octree

       type(octal), intent(in), target :: thisOctal
       logical, intent(in)             :: fileFormatted
       type(gridtype) :: grid
       type(octal), pointer :: thisChild
       integer              :: iChild
       
       if (fileFormatted) then 
          write(iostat=error,fmt=*,unit=20) thisOctal%nDepth, thisOctal%nChildren,&
                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
                  thisOctal%subcellSize,thisOctal%threed, thisOctal%twoD,    &
                  thisOctal%maxChildren, thisOctal%dustType,  &
                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
                  thisOctal%dPhi, thisOctal%r
       else
          write(iostat=error,unit=20) thisOctal%nDepth, thisOctal%nChildren, &
                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
                  thisOctal%probDistCont, thisOctal%Ne, thisOctal%nh, thisOctal%nTot,      &
                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
                  thisOctal%subcellSize, thisOctal%threeD, thisOctal%twoD,   &
                  thisOctal%maxChildren, thisOctal%dustType, &
                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
                  thisOctal%dPhi, thisOctal%r
       end if 
       if (.not.grid%oneKappa) then
!          call writeReal2D(thisOctal%kappaAbs,fileFormatted)
!          call writeReal2D(thisOctal%kappaSca,fileFormatted)
          call writeDouble2D(thisOctal%kappaAbs,fileFormatted)
          call writeDouble2D(thisOctal%kappaSca,fileFormatted)
       endif
       if (grid%photoionization) then
          call writeDouble2D(thisOctal%ionFrac,fileFormatted)
       endif
       if (molecular) then
          call writeDouble2D(thisOctal%molecularLevel,fileFormatted)
          if (fileformatted) then
             write(unit=20,iostat=error,fmt=*) thisOctal%microturb, thisOctal%nh2
          else
             write(unit=20,iostat=error) thisOctal%microturb, thisOctal%nh2
          endif
       endif
       call writeDouble2D(thisOctal%N,fileFormatted)
       call writeReal2D(thisOctal%departCoeff,fileFormatted)
       call writeDouble2D(thisOctal%dustTypeFraction, fileFormatted)


       if (cmf) then
          call writeDouble3D(thisOctal%atomLevel,fileFormatted)
          if (fileformatted) then
             write(unit=20,iostat=error,fmt=*) thisOctal%microturb
          else
             write(unit=20,iostat=error) thisOctal%microturb
          endif
       endif
       
       if (thisOctal%nChildren > 0) then 
          do iChild = 1, thisOctal%nChildren, 1
             thisChild => thisOctal%child(iChild)
             call writeOctreePrivate(thisChild,fileFormatted,grid)
          end do
       end if

    end subroutine writeOctreePrivate
    
    
  end subroutine writeAMRgrid


  
  subroutine readAMRgrid(filename,fileFormatted,grid)
    ! reads in a previously saved 'grid' for an adaptive mesh geometry  

    use input_variables, only: geometry,dipoleOffset,amr2dOnly,statEq2d, molecular, cmf
    implicit none

    character(len=*)            :: filename
    logical, intent(in)         :: fileFormatted
    type(GRIDTYPE), intent(inout) :: grid
    
    integer, dimension(8) :: timeValues    ! system date and time
    integer               :: dummy         
    integer               :: error         ! status code
    logical               :: octreePresent ! true if grid has an octree
    integer :: nOctal
    character(len=80) :: absolutePath, inFile
    integer :: iJunk
    real, pointer :: junk(:), junk2(:,:), junk3(:,:)

    error = 0


  call unixGetEnv("TORUS_JOB_DIR",absolutePath)
  inFile = trim(absolutePath)//trim(filename)

    if (fileFormatted) then
       open(unit=20, iostat=error, file=inFile, form="formatted", status="old")
       if (error /=0) then
         print *, 'Panic: file open error in readAMRgrid, file:',trim(inFile) ; stop
       end if
       ! read the file's time stamp
       read(unit=20,fmt=*,iostat=error) timeValues 
       if (error /=0) then
         print *, 'Panic: read error in readAMRgrid (formatted timeValues)' ; stop
       end if
    else
       open(unit=20, iostat=error, file=inFile, form="unformatted", status="old")
       if (error /=0) then
         print *, 'Panic: file open error in readAMRgrid, file:',trim(inFile) ; stop
       end if
       ! read the file's time stamp
       read(unit=20,iostat=error) timeValues
       if (error /=0) then
         print *, 'Panic: read error in readAMRgrid (unformatted timeValues)' ; stop
       end if
    end if

    if (writeoutput) write(*,'(a,a)') "Reading populations file from: ",trim(filename)
    if (writeoutput) write(*,'(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') ' - data file written at: ', &
                          timeValues(1),'/',timeValues(2),'/',&
                          timeValues(3),'  ',timeValues(5),':',timeValues(6)
                          
    ! read the variables to be stored in the top-level 'grid' structure
    if (fileFormatted) then
       read(unit=20,fmt=*,iostat=error) ijunk, grid%flatSpec, grid%adaptive,& 
               grid%cartesian, grid%isotropic, grid%hitCore, grid%diskRadius, &
               grid%diskNormal, grid%DipoleOffset, grid%geometry, grid%rCore, &
               grid%lCore, grid%chanceWindOverTotalContinuum,                 &
               grid%lineEmission, grid%contEmission, grid%doRaman,            &
               grid%resonanceLine, grid%rStar1, grid%rStar2, grid%lumRatio,   &
               grid%tempSource, grid%starPos1, grid%starPos2, grid%lambda2,   &
               grid%maxLevels, grid%maxDepth, grid%halfSmallestSubcell,       &
               grid%nOctals, grid%smoothingFactor, grid%oneKappa, grid%rInner,&
               grid%rOuter, grid%amr2dOnly, grid%photoionization
    else
       read(unit=20,iostat=error) ijunk, grid%flatSpec, grid%adaptive, & 
               grid%cartesian, grid%isotropic, grid%hitCore, grid%diskRadius, &
               grid%diskNormal, grid%DipoleOffset, grid%geometry, grid%rCore, &
               grid%lCore, grid%chanceWindOverTotalContinuum,                 &
               grid%lineEmission, grid%contEmission, grid%doRaman,            &
               grid%resonanceLine, grid%rStar1, grid%rStar2, grid%lumRatio,   &
               grid%tempSource, grid%starPos1, grid%starPos2, grid%lambda2,   &
               grid%maxLevels, grid%maxDepth, grid%halfSmallestSubcell,       &
               grid%nOctals, grid%smoothingFactor, grid%oneKappa, grid%rInner,& 
               grid%rOuter, grid%amr2dOnly, grid%photoionization
    end if    

    if (error /=0) then
       print *, 'Panic: read error in readAMRgrid 1'
       stop
    end if
!    call readReal1D(grid%lamarray,fileFormatted)
!    call readReal2D(grid%oneKappaAbs,fileFormatted)
!    call readReal2D(grid%oneKappaSca,fileFormatted)


    call readReal1D(junk,fileFormatted)
    call readReal2D(junk2,fileFormatted)
    call readReal2D(junk3,fileFormatted)
    deallocate(junk, junk2, junk3)

    call readClumps(fileFormatted)
    

    ! now we call the recursive subroutine to read the tree structure 
    if (fileFormatted) then
       read(unit=20,fmt=*) octreePresent
    else
       read(unit=20) octreePresent
    end if
     
    if (octreePresent) then
       allocate(grid%octreeRoot)
       nOctal = 0
       call readOctreePrivate(grid%octreeRoot,null(),fileFormatted, nOctal, grid)
!       write(*,*) noctal,"octals read"
    end if

    ! check that we are at the end of the file
    if (fileFormatted) then
       read(unit=20,fmt=*, iostat=error) dummy
    else
       read(unit=20, iostat=error) dummy
    end if
    if (error == 0) then
       print *, 'Panic: read error (expected end of file) in readAMRgrid' ; stop
    end if
    
    close(unit=20)

    if (writeoutput) then
       print *, 'setting ''geometry'':',trim(geometry),' previously:',trim(grid%geometry)
       print *, 'setting ''dipoleOffset'':',dipoleOffset,' previously:',grid%dipoleOffset
       print *, 'setting ''amr2donly'':',amr2donly,' previously:',grid%amr2donly
       print *, 'setting ''statEq2d'':',' previously:', grid%statEq2d
    endif
    grid%geometry = trim(geometry)
    grid%dipoleOffset = dipoleOffset
    grid%amr2donly = amr2donly
    grid%statEq2d = statEq2d

  contains
   
    recursive subroutine readOctreePrivate(thisOctal,parent,fileFormatted, noctal, grid)
       ! read in an octal to the grid octree

       implicit none
       type(octal), pointer :: thisOctal
       type(octal), pointer :: parent
       type(gridtype) :: grid

       logical, intent(in)  :: fileFormatted
       integer :: nOctal

       type(octal), pointer :: thisChild
       integer              :: iChild

       nOctal = nOctal+1
       thisOctal%parent => parent
       
       if (fileFormatted) then
          read(unit=20,fmt=*,iostat=error) thisOctal%nDepth, thisOctal%nChildren,  &
                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
                  thisOctal%subcellSize, thisOctal%threeD, thisOctal%twoD,   &
                  thisOctal%maxChildren, thisOctal%dustType, &
                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
                  thisOctal%dPhi, thisOctal%r
       else
          read(unit=20,iostat=error) thisOctal%nDepth, thisOctal%nChildren,  &
                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
                  thisOctal%subcellSize, thisOctal%threeD,  thisOctal%twoD,  &
                  thisOctal%maxChildren, thisOctal%dustType, &
                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
                  thisOctal%dPhi, thisOctal%r

       end if


       if (.not.grid%oneKappa) then
!          call readReal2D(thisOctal%kappaAbs,fileFormatted)
!          call readReal2D(thisOctal%kappaSca,fileFormatted)
          call readDouble2D(thisOctal%kappaAbs,fileFormatted)
          call readDouble2D(thisOctal%kappaSca,fileFormatted)
       endif
       if (grid%photoionization) then
          call readDouble2D(thisOctal%ionFrac,fileFormatted)
       endif
       if (molecular) then
          call readDouble2D(thisOctal%molecularLevel,fileFormatted)
          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%microturb, thisOctal%nh2
          else
             read(unit=20,iostat=error) thisOctal%microturb, thisOctal%nh2
          endif

       endif
       call readDouble2D(thisOctal%N,fileFormatted)
       call readReal2D(thisOctal%departCoeff,fileFormatted)
       call readDouble2D(thisOctal%dustTypeFraction, fileFormatted)

       if (cmf) then
          call readDouble3D(thisOctal%atomLevel,fileFormatted)
          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%microturb
          else
             read(unit=20,iostat=error) thisOctal%microturb
          endif

       endif


       if (thisOctal%nChildren > 0) then 
          allocate(thisOctal%child(1:thisOctal%nChildren)) 
          do iChild = 1, thisOctal%nChildren, 1
             thisChild => thisOctal%child(iChild)
             call readOctreePrivate(thisChild,thisOctal,fileFormatted, nOctal, grid)
          end do
       end if

    end subroutine readOctreePrivate
 
  end subroutine readAMRgrid

  
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
    if (writeoutput) write(*,'(a,a)') "Reading populations file from: ",trim(filename)
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
    if (writeoutput) write(*,'(a,a)') "Reading populations file from: ",trim(filename)
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
    real :: o6width
    integer, parameter :: maxr = 200
    
    integer :: kPrime
    real, allocatable :: ionPattern(:,:,:,:)
    real, allocatable :: rArray(:), yArray(:), yarray2(:)
    integer, allocatable :: ir(:,:)
    integer :: nTheta, nPhi, nr
    integer :: iTheta, iPhi, i1
    real :: fac, fac2, r, theta, phi
    type(VECTOR) :: starPos1, starPos2, rVec
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
             grid%velocity(i,j,k) = grid%velocity(i,j,k) / real(cSpeed)
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
                grid%kappaSca(i,j,k,1) = max(1.e-30,grid%rho(i,j,k) * real(sigmaE))
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
    integer :: i1, i2, i3, j1, j2
    real :: zPrime, yPrime, r, fac, ang
    logical :: resonanceLine
    type(VECTOR) :: rHat, rVec, vHat
    real :: surfaceVel, vel, theta
    integer :: nz, ny
    real :: t1, t2, t3
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
               r * sin(grid%phiAxis(i2)), &
               -thickness/2.)
          call getIndices(grid, rVec, i, j1, k, t1, t2, t3)
          rVec = VECTOR(r * cos(grid%phiAxis(i2)), &
               r * sin(grid%phiAxis(i2)), &
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
    integer :: i, j
    integer :: i1, i2, i3, j1, j2
    real :: zPrime, r, fac, ang
    logical :: resonanceLine
    character(len=*) :: misc
    type(VECTOR) :: rHat, rVec, vHat
    real :: surfaceVel, vel
    integer :: nz, ny
    real :: t1, t2, t3
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





  subroutine fillGridWind(grid, mDot, rStar, tEff, v0, vterm, beta, &
       lte, contFluxFile, writePops, readPops, popFilename, nLower, nUpper, vrot)

    implicit none
    real :: mDot, rStar, tEff, v0, vTerm, beta, vrot
    integer :: i, j, k
    character(len=*) :: popFilename, contFluxFile
    logical :: readPops, writePops, lte
    integer :: nLower, nUpper
    real :: sinTheta
    type(VECTOR) :: rHat, rVec, vvec, spinAxis
    real :: vel, beta1, beta2, vext, v1, r
    real :: x, dx, dv, mcore
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
       x = x + dx
       dx = dx * 1.2
    enddo

    
    grid%rCore = grid%rAxis(1)
    grid%rStar1 = grid%rCore
    grid%rStar2 = 0.
    grid%starPos1 = VECTOR(0.,0.,0.)
    grid%starPos2 = VECTOR(1.e10,0.,0.)
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
             v1 = vterm* (3.*abs(grid%muAxis(j))+1.)/4.
             vel = v0 + (v1 - v0) * (1. - grid%rAxis(1)/grid%rAxis(i))**beta
             dv = grid%rAxis(1)*((vTerm-v0)*beta*(1. - grid%rAxis(1)/grid%rAxis(i))**(beta-1.))/grid%rAxis(i)**2

!             beta1 = 1.5
!             beta2 = 3.
!             vext = 500.*1.e5
!             vel = v0 + (vterm-vext-v0)*(1. -grid%rAxis(1)/grid%rAxis(i))**beta1 + &
!                  vext*(1. -grid%rAxis(1)/grid%rAxis(i))**beta2
!            dv = (grid%rAxis(1)*beta1*(vTerm-Vext-v0)*(1. -grid%rAxis(1)/grid%rAxis(i))**(beta1-1.))/grid%rAxis(i)**2 &
!                  + (grid%rAxis(1)*beta2*vext*(1.-grid%rAxis(1)/grid%rAxis(i))**(beta2-1.))/grid%rAxis(i)**2

             grid%velocity(i,j,k) = (vel / cSpeed) * rHat
             grid%dVbyDr(i,j,k) = dv * rHat
             grid%rho(i,j,k) = mDot / (fourPi * vel * grid%rAxis(i)**2 * 1.e20)
             grid%inStar(i,j,k) = .false.
             grid%inUse(i,j,k) = .true.
          enddo
       enddo
    enddo

    if (vrot /= 0.) then
       mCore = 10. * mSol
       spinAxis = VECTOR( 0., 0., 1.)
       do i = 1, grid%nr
          do j = 1, grid%nMu
             do k = 1, grid%nPhi
                sinTheta = sqrt(1.-grid%muAxis(j)**2)
                rVec = VECTOR(grid%rAxis(i)*cos(grid%phiAxis(k))*sinTheta, &
                     grid%rAxis(i)*sin(grid%phiAxis(k))*sinTheta, &
                     grid%rAxis(i)*grid%muAxis(j))
                r = grid%rAxis(i)
                rVec = rVec / r
                vVec = rVec  .cross. spinAxis
                if (modulus(vVec) /= 0.) then
                   call normalize(vVec)
                   vel = sqrt(bigG * mCore / (r*1.e10)) / cSpeed
                   vVec = vel  * vVec 
                   grid%velocity(i,j,k) = grid%velocity(i,j,k) +  vVec
                endif
             enddo
          enddo
       enddo
    endif



!    grid%rho = grid%rho * 0.375
!    write(*,'(a)') "! Assuming N(H)/N(He) = 1.5"
    grid%etaCont = 1.e-30

  end subroutine fillGridWind



  
  subroutine writeAxes(grid)
    type(GRIDTYPE) :: grid
    integer :: i


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
    real :: sinTheta, vel, v0
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
    real :: muAxis(:)
    integer :: i,  n1, n2
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

  subroutine writeReal1D(variable,fileFormatted)

     real,dimension(:),pointer :: variable
     logical, intent(in)       :: fileFormatted
     integer :: error

     if (fileFormatted) then
        if (associated(variable)) then
           write(unit=20,fmt=*,iostat=error) .true.
           if (error /=0) then
             print *, 'Panic: write error in writeReal1D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) SIZE(variable)
           if (error /=0) then
             print *, 'Panic: write error in writeReal1D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then
             print *, 'Panic: write error in writeReal1D' ; stop
           end if
        else
           write(unit=20,fmt=*,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeReal1D' ; stop
           end if
        end if
     else
        if (associated(variable)) then
           write(unit=20,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeReal1D' ; stop
           end if
           write(unit=20,iostat=error) SIZE(variable)
           if (error /=0) then 
             print *, 'Panic: write error in writeReal1D' ; stop
           end if
           write(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: write error in writeReal1D' ; stop
           end if
        else
           write(unit=20,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeReal1D' ; stop
           end if
        end if
     end if
             
  end subroutine writeReal1D
    
  subroutine writeReal2D(variable,fileFormatted)

     real,dimension(:,:),pointer :: variable
     logical, intent(in)         :: fileFormatted
     integer :: error

     if (fileFormatted) then
        if (associated(variable)) then
           write(unit=20,fmt=*,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeReal2D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) SIZE(variable,1),SIZE(variable,2)
           if (error /=0) then 
             print *, 'Panic: write error in writeReal2D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: write error in writeReal2D' ; stop
           end if
        else
           write(unit=20,fmt=*,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeReal2D' ; stop
           end if
        end if
     else
        if (associated(variable)) then
           write(unit=20,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeReal2D' ; stop
           end if
           write(unit=20,iostat=error) SIZE(variable,1),SIZE(variable,2)
           if (error /=0) then 
             print *, 'Panic: write error in writeReal2D' ; stop
           end if
           write(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: write error in writeReal2D' ; stop
           end if
        else
           write(unit=20,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeReal2D' ; stop
           end if
        end if
     end if 
    
  end subroutine writeReal2D

  subroutine writeDouble2D(variable,fileFormatted)

     real(double),dimension(:,:),pointer :: variable
     logical, intent(in)                          :: fileFormatted
     integer :: error

     if (fileFormatted) then
        if (associated(variable)) then
           write(unit=20,fmt=*,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble2D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) SIZE(variable,1),SIZE(variable,2)
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble2D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble2D' ; stop
           end if
        else
           write(unit=20,fmt=*,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble2D' ; stop
           end if
        end if
     else
        if (associated(variable)) then
           write(unit=20,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble2D' ; stop
           end if
           write(unit=20,iostat=error) SIZE(variable,1),SIZE(variable,2)
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble2D' ; stop
           end if
           write(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble2D' ; stop
           end if
        else
           write(unit=20,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble2D' ; stop
           end if
        end if
     end if
    
  end subroutine writeDouble2D

  subroutine writeDouble3D(variable,fileFormatted)

     real(double),dimension(:,:,:),pointer :: variable
     logical, intent(in)                          :: fileFormatted
     integer :: error

     if (fileFormatted) then
        if (associated(variable)) then
           write(unit=20,fmt=*,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble3D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) SIZE(variable,1),SIZE(variable,2), SIZE(variable,3)
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble3D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble3D' ; stop
           end if
        else
           write(unit=20,fmt=*,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble3D' ; stop
           end if
        end if
     else
        if (associated(variable)) then
           write(unit=20,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble3D' ; stop
           end if
           write(unit=20,iostat=error) SIZE(variable,1),SIZE(variable,2), SIZE(variable,3)
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble3D' ; stop
           end if
           write(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble3D' ; stop
           end if
        else
           write(unit=20,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble3D' ; stop
           end if
        end if
     end if
    
  end subroutine writeDouble3D

  subroutine writeClumps(fileFormatted)

     use clump_mod
     logical, intent(in)       :: fileFormatted
     integer :: error

     if (fileFormatted) then
        if (allocated(clumps)) then
           write(unit=20,fmt=*,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeClumps' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) SIZE(clumps)
           if (error /=0) then 
             print *, 'Panic: write error in writeClumps' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) clumps
           if (error /=0) then 
             print *, 'Panic: write error in writeClumps' ; stop
           end if
        else
           write(unit=20,fmt=*,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeClumps' ; stop
           end if
        end if
     else
        if (allocated(clumps)) then
           write(unit=20,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeClumps' ; stop
           end if
           write(unit=20,iostat=error) SIZE(clumps)
           if (error /=0) then 
             print *, 'Panic: write error in writeClumps' ; stop
           end if
           write(unit=20,iostat=error) clumps
           if (error /=0) then 
             print *, 'Panic: write error in writeClumps' ; stop
           end if
        else
           write(unit=20,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeClumps' ; stop
           end if
        end if
     end if
             
  end subroutine writeClumps
  
  subroutine readClumps(fileFormatted)
 
     use clump_mod
  
     logical, intent(in)       :: fileFormatted
     logical                   :: present 
     integer                   :: length
     integer                   :: error

     if (allocated(clumps)) then
       print *, 'Clearing existing ''clumps'' variable'
       deallocate(clumps)
     end if
     
     if (fileFormatted) then
        read(unit=20,fmt=*,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readClumps' ; stop
        end if
        if (present) then
           read(unit=20,fmt=*,iostat=error) length
           if (error /=0) then 
             print *, 'Panic: read error in readClumps' ; stop
           end if
           allocate(clumps(length),stat=error)
           if (error /=0) then 
             print *, 'Panic: in readClumps, can''t allocate clumps(',length,')' ; stop
           end if
           read(unit=20,fmt=*,iostat=error) clumps
           if (error /=0) then 
             print *, 'Panic: read error in readClumps' ; stop
           end if
        end if
     else
        read(unit=20,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readClumps' ; stop
        end if
        if (present) then
           read(unit=20,iostat=error) length
           if (error /=0) then 
             print *, 'Panic: read error in readClumps' ; stop
           end if
           allocate(clumps(length),stat=error)
           if (error /=0) then 
             print *, 'Panic: in readClumps, can''t allocate clumps(',length,')' ; stop
           end if
           read(unit=20,iostat=error) clumps
           if (error /=0) then 
             print *, 'Panic: read error in readClumps' ; stop
           end if
        end if
     end if
    
  end subroutine readClumps
 
  subroutine readReal1D(variable,fileFormatted)

     real,dimension(:),pointer :: variable
     logical, intent(in)       :: fileFormatted
     logical                   :: present 
     integer                   :: length
     integer                   :: error

     if (fileFormatted) then
        read(unit=20,fmt=*,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readReal1D' ; stop
        end if
        if (present) then
           read(unit=20,fmt=*,iostat=error) length
           if (error /=0) then 
             print *, 'Panic: read error in readReal1D' ; stop
           end if
           allocate(variable(length))
           read(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readReal1D' ; stop
           end if
        end if
     else
        read(unit=20,iostat=error) present
        if (present) then
           read(unit=20,iostat=error) length
           if (error /=0) then 
             print *, 'Panic: read error in readReal1D' ; stop
           end if
           allocate(variable(length))
           read(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readReal1D' ; stop
           end if
        end if
     end if
    
  end subroutine readReal1D
    
  subroutine readReal2D(variable,fileFormatted)

     real,dimension(:,:),pointer :: variable
     logical, intent(in)       :: fileFormatted
     logical                   :: present 
     integer                   :: length1
     integer                   :: length2
     integer                   :: error

     if (fileFormatted) then
        read(unit=20,fmt=*,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readReal2D' ; stop
        end if
        if (present) then
           read(unit=20,fmt=*,iostat=error) length1, length2
           if (error /=0) then 
             print *, 'Panic: read error in readReal2D' ; stop
           end if
           allocate(variable(length1,length2))
           read(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readReal2D' ; stop
           end if
        end if
     else
        read(unit=20,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readReal2D' ; stop
        end if
        if (present) then
           read(unit=20,iostat=error) length1, length2
           if (error /=0) then 
             print *, 'Panic: read error in readReal2D' ; stop
           end if
           allocate(variable(length1,length2))
           read(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readReal2D' ; stop
           end if
        end if
     endif
    
  end subroutine readReal2D

  subroutine readDouble2D(variable,fileFormatted)

     real(double),dimension(:,:),pointer :: variable
     logical, intent(in)       :: fileFormatted
     logical                   :: present 
     integer                   :: length1
     integer                   :: length2
     integer                   :: error
     
     if (fileFormatted) then
        read(unit=20,fmt=*,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readDouble2D' ; stop
        end if
        if (present) then
           read(unit=20,fmt=*,iostat=error) length1, length2
           if (error /=0) then 
             print *, 'Panic: read error in readDouble2D' ; stop
           end if
           allocate(variable(length1,length2))
           read(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readDouble2D' ; stop
           end if
        end if
     else
        read(unit=20,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readDouble2D' ; stop
        end if
        if (present) then
           read(unit=20,iostat=error) length1, length2
           if (error /=0) then 
             print *, 'Panic: read error in readDouble2D' ; stop
           end if
           allocate(variable(length1,length2))
           read(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readDouble2D' ; stop
           end if
        end if
     end if

  end subroutine readDouble2D

  subroutine readDouble3D(variable,fileFormatted)

     real(double),dimension(:,:,:),pointer :: variable
     logical, intent(in)       :: fileFormatted
     logical                   :: present 
     integer                   :: length1
     integer                   :: length2
     integer                   :: length3
     integer                   :: error
     
     if (fileFormatted) then
        read(unit=20,fmt=*,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readDouble3D' ; stop
        end if
        if (present) then
           read(unit=20,fmt=*,iostat=error) length1, length2, length3
           if (error /=0) then 
             print *, 'Panic: read error in readDouble3D' ; stop
           end if
           allocate(variable(length1,length2,length3))
           read(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readDouble3D' ; stop
           end if
        end if
     else
        read(unit=20,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readDouble2D' ; stop
        end if
        if (present) then
           read(unit=20,iostat=error) length1, length2, length3
           if (error /=0) then 
             print *, 'Panic: read error in readDouble3D' ; stop
           end if
           allocate(variable(length1,length2, length3))
           read(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readDouble3D' ; stop
           end if
        end if
     end if

  end subroutine readDouble3D

  subroutine initTestAmr(grid)
    use input_variables
    type(GRIDTYPE) :: grid

    grid%geometry = "testamr"

    grid%rCore = rCore / 1.e10
    grid%rStar1 = rCore / 1.e10
    grid%starpos1 = VECTOR(1.e-30,1.e-30,1.e-30)
    grid%lCore = fourPi * stefanBoltz * grid%rCore**2 * 1.e20 * teff**4
    grid%rInner = rInner / 1.e10
    grid%rOuter = rOuter / 1.e10
    grid%lineEmission = .false.
  end subroutine initTestAmr

  subroutine initProtoAmr(grid)
    use input_variables
    type(GRIDTYPE) :: grid

    grid%geometry = "proto"

    grid%rCore = rCore / 1.e10
    grid%lCore = fourPi * stefanBoltz * grid%rCore**2 * 1.e20 * teff**4
    grid%rInner = rInner / 1.e10
    grid%rOuter = rOuter / 1.e10
    grid%lineEmission = .false.
  end subroutine initProtoAmr



  subroutine fancyAMRplot(grid, device)
    type(GRIDTYPE) :: grid
    character(len=*) device
    integer :: i, pgbeg

    i =  pgbeg(0,device, 2, 2)
    call pgpage
    call dofancyAMRplot(grid, device, 1.)
    call pgpage
    call dofancyAMRplot(grid, device, 0.1)
    call pgpage
    call dofancyAMRplot(grid, device, 0.01)
    call pgpage
    call dofancyAMRplot(grid, device, 0.001)
    call pgend
  end subroutine fancyAMRplot


  subroutine dofancyAMRplot(grid, device, scaleFac)
    type(GRIDTYPE) :: grid
    character(len=*) device
    type(OCTAL), pointer ::  thisOctal
    real :: scaleFac
    integer :: oldSubcell
    real :: rhoMin, rhoMax
    integer :: ilo, ihi, i

    write(*,*) "Fancy plotting to: ",trim(device)

    rhoMin = 1.e30
    rhoMax = -1.e30
    thisOctal => grid%octreeRoot
    call minMaxDensity(thisOctal, rhoMin, rhoMax)

    write(*,*) "Density range",rhoMin,rhoMax

!    call pgbegin(0,device,1,1)

    call pgvport(0.1,0.9,0.1,0.9)

    call pgwnad(real(-grid%octreeRoot%subcellsize)*scaleFac, real(grid%octreeRoot%subcellsize)*scaleFac, &
         real(-grid%octreeRoot%subcellsize)*scaleFac, real(grid%octreeRoot%subcellsize)*scaleFac)
    oldSubcell =0


    call pgqcir(ilo, ihi)
    call pgscir(ilo, ihi)
    call palette(2)
    if (device == "/vps") then
       i = ilo
       ilo = ihi
       ihi = i 
    endif
    thisOctal => grid%octreeRoot
    call plotZplane(thisOctal, rhoMax, rhoMin, ilo, ihi)
    call pgsci(1)
    call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
    call pgend
  end subroutine dofancyAMRplot

  recursive subroutine plotZplane(thisOctal, rhoMax, rhoMin, ilo, ihi)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(octalvector) :: rvec
  real :: rho
  real :: rhoMin, rhoMax
  integer :: subcell, i, ilo, ihi, idx
  real :: t
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call plotZplane(child, rhoMax, rhoMin, ilo, ihi)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          if (((rVec%z + thisOctal%subcellSize/2.) >= 0.).and. &
              ((rVec%z - thisOctal%subcellSize/2.) <= 0.)) then
             
             rho = thisOctal%rho(subcell)
             if (rho == 0.) rho = rhoMin
             t = (log(rho)-log(rhoMin))/(log(rhoMax)-log(rhoMin))
             if (t < 0.) t = 0.
             idx = int(t * real(ihi - ilo) + real(ilo))
             call pgsci(idx)
             call pgrect(real(rVec%x-thisOctal%subcellSize/2.), real(rVec%x + thisOctal%subcellSize/2.), &
                  real(rVec%y-thisOctal%subcellSize/2.),real(rVec%y+thisOctal%subcellSize/2.))
             call pgsci(1)
             call pgmove(real(rVec%x-thisOctal%subcellSize/2.), &
                  real(rVec%y-thisOctal%subcellSize/2.))
             call pgdraw(real(rVec%x+thisOctal%subcellSize/2.), &
                  real(rVec%y-thisOctal%subcellSize/2.))
             call pgdraw(real(rVec%x+thisOctal%subcellSize/2.), &
                  real(rVec%y+thisOctal%subcellSize/2.))
             call pgdraw(real(rVec%x-thisOctal%subcellSize/2.), &
                  real(rVec%y+thisOctal%subcellSize/2.))
             call pgdraw(real(rVec%x-thisOctal%subcellSize/2.), &
                  real(rVec%y-thisOctal%subcellSize/2.))
          endif
       endif
    enddo
    
  end subroutine plotZplane

  
  
  ! It overlay the density structure and the grid structure using PGPLOT
  ! routines. 
  ! For plotting density, we use a generic "density" function which returns
  ! the density for a given postion and the type of density.
  !
  subroutine draw_cells_on_density(grid, plane, device)
    implicit none
    type(GRIDTYPE), intent(in) :: grid
    ! must be 'x-y', 'y-z' , 'z-x' or 'x-z'
    character(LEN=*), intent(in)  :: plane
    character(len=*), intent(in)  :: device
    type(OCTAL), pointer :: root

    !integer :: thisSubcell, oldSubcell
    real :: d
    
    ! For plotting density
    integer :: pgbeg
    integer :: n, ncol, nlev
    parameter (n=1500, ncol=32, nlev=10)
    integer :: i,j,ci1,ci2
    real :: f(n,n),fmin,fmax,tr(6)
    real(double)  :: xmap(n), ymap(n), zmap(n), x1, y1, z1
    real :: box_size
    real :: offset
    real(double)  :: tmp
    type(octalVector) :: pos,dir
    
    
    write(*,*) "draw_cells_on_density plotting to: ",trim(device)

    ! retriving the address to the root of the tree.
    root => grid%octreeRoot
    
    ! setup the density 
    box_size = REAL(grid%octreeRoot%subcellsize)*2.0
    d = box_size

    ! The transformation matrix TR is used to calculate the world
    ! coordinates of the center of the cell that represents each
    ! array element. The world coordinates of the center of the cell
    ! corresponding to array element f(I,J) are given by
    !        X = TR(1) + TR(2)*I + TR(3)*J
    !        Y = TR(4) + TR(5)*I + TR(6)*J
    
    TR = 0.0
    TR(1) = -d/2.0
    TR(4) = -d/2.0
    TR(2) = d/real(n)
    TR(6) = d/real(n)
    
    
    x1 =  - d/2.0d0
    y1 =  - d/2.0d0
    z1 =  - d/2.0d0
    
    do i = 1, n
       xmap(i) = x1 + dble(i-1)*d/dble(n-1) 
       ymap(i) = y1 + dble(i-1)*d/dble(n-1) 
       zmap(i) = z1 + dble(i-1)*d/dble(n-1) 
    end do
    
    fmin = 1.0+200.0
    fmax = -1.0-200.0
    
    if (plane(1:3) == 'x-y') then  
       do j = 1, n
          do i = 1, n


             pos = Vector(xmap(i),ymap(j),-d/2.0d0)
             dir = Vector(0.0d0, 1.0d0, 1.0d0)


             ! using the function in density_mod.f90
             tmp = ABS( density(pos, grid) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
             


             pos = Vector(xmap(i),ymap(j),0.0d0)
             ! using the function in density_mod.f90
             tmp = ABS( density(pos, grid) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
          end do
       end do
    elseif (plane(1:3) == 'y-z') then  
       do j = 1, n
          do i = 1, n
             pos = Vector(0.0d0, ymap(i), zmap(j))
             ! using the function in density_mod.f90
            tmp = ABS( density(pos, grid) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
          end do
       end do
    elseif (plane(1:3) == 'z-x') then  
       do j = 1, n
          do i = 1, n
             pos = Vector(xmap(j),0.0d0, zmap(i))
             ! using the function in density_mod.f90
             tmp = ABS( density(pos, grid) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
          end do
       end do
    elseif (plane(1:3) == 'x-z') then  
       do j = 1, n
          do i = 1, n
             pos = Vector(xmap(i),0.0d0, zmap(j))
             ! using the function in density_mod.f90
             tmp = ABS( density(pos, grid) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
          end do
       end do
    end if

    ! Just for safty
    if(fmin == fmax) then
       fmin = fmin-fmin/10.0
       fmax = fmax+fmax/10.0
    end if

    ! Open plot device and set up coordinate system. We will plot the
    !image within a unit square.
    !
    IF (PGBEG(0,device,1,1) .NE. 1) STOP
    CALL PGQCOL(CI1, CI2)
    IF (CI2.LT. 15+NCOL) THEN
       WRITE (*,*) 'This program requires a device with at least',&
            15+NCOL,' colors'
       STOP
    END IF
    CALL PGPAGE

    !    CALL PGSCR(0, 0.0, 0.3, 0.2)
    !    CALL PGSVP(0.05,0.95,0.05,0.95)

    call setvp

    if (plane(1:3) == 'x-y') then
       CALL PGWNAD(real(xmap(1)), real(xmap(n)), real(ymap(1)), real(ymap(n)))
    elseif (plane(1:3) == 'y-z') then
       CALL PGWNAD(real(ymap(1)), real(ymap(n)), real(zmap(1)), real(zmap(n)))
    elseif (plane(1:3) == 'z-x') then
       CALL PGWNAD(real(zmap(1)), real(zmap(n)), real(xmap(1)), real(xmap(n)))
    elseif (plane(1:3) == 'x-z') then
       CALL PGWNAD(real(xmap(1)), real(xmap(n)), real(zmap(1)), real(zmap(n)))
    end if
    
        
    ! Setting the color range and etc..e
    call palett(2, 1.0, 0.5)
    
    ! Plot the image
    CALL PGIMAG(f,N,N, 1, N, 1, N,  fmin, fmax, TR)


    !
    ! Annotation.
    !
    CALL PGSCI(1)
    CALL PGMTXT('t',1.0,0.0,0.0,' ')
    !CALL PGBOX('bcnts',0.0,0,'bcnts',0.0,0)
    CALL PGBOX('bcntsi',0.0,0,'bcntsiv',0.0,0)
    
    if (plane(1:3) == 'x-y') then
       CALL PGMTXT('B',3.0,1.0,1.0,'x')
       CALL PGMTXT('L',3.0,1.0,1.0,'y')
    elseif (plane(1:3) == 'y-z') then
       CALL PGMTXT('B',3.0,1.0,1.0,'y')
       CALL PGMTXT('L',3.0,1.0,1.0,'z')
    elseif (plane(1:3) == 'z-x') then 
       CALL PGMTXT('B',3.0,1.0,1.0,'z')
       CALL PGMTXT('L',3.0,1.0,1.0,'x')
    elseif (plane(1:3) == 'x-z') then 
       CALL PGMTXT('B',3.0,1.0,1.0,'x')
       CALL PGMTXT('L',3.0,1.0,1.0,'z')
    end if
    
    ! draw boxies
    call PGSFS(2)  ! we don't want to fill in a box

    ! this is a recursive routine
    call draw_rectangle(root, plane)

    ! routines in pixplot_module    
    
    ! draw wedge
    offset = 0
    CALL PGWEDG('BI', 4.0, 5.0, FMIN, FMAX, 'pixel value')
    CALL PGSCH(0.6)

!    CALL FIDDLE
    CALL PGASK(.FALSE.)  

    call PGEND
    
    ! deallocate(a_root)
  end subroutine draw_cells_on_density



  !
  ! Plots column density  (g/cm^3) on log scale on a specified plane
  ! (This is based on draw_cells_on_density routine above.)
  subroutine plot_column_density(grid, plane, device,  &
       nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim, value_3rd_dim) 
    implicit none
    type(GRIDTYPE), intent(in) :: grid
    ! must be 'x-y', 'y-z' , 'z-x' or 'x-z'
    character(LEN=*), intent(in)  :: plane
    character(len=*), intent(in)  :: device
    ! To put the markers to indicate the important points in the plots.
    integer, intent(in) :: nmarker           ! number of markers
    real, intent(in)    :: xmarker(nmarker)  ! position of x
    real, intent(in)    :: ymarker(nmarker)  ! position of y
    real, intent(in)    :: zmarker(nmarker)  ! position of z
    real, intent(in)    :: width_3rd_dim     ! Use this to restrict the markers to be plotted..  
    logical, intent(in) :: show_value_3rd_dim! If T, the value of the third dimension will be shown.
    real, intent(in)              :: value_3rd_dim 
    !
    type(OCTAL), pointer :: root
    real :: d, dd, off
    real :: xc, yc , zc, v3

    ! For plotting density
    integer :: pgbeg
    integer :: n, ncol, nlev
!    parameter (n=1000, ncol=32, nlev=10)
    parameter (n=100, ncol=32, nlev=10)
    integer :: i,j,ci1,ci2
    real :: f(n,n),fmin,fmax,tr(6)
    real(double)  :: xmap(n), ymap(n), zmap(n), x1, y1, z1
    real :: box_size
    real :: offset
    real(double)  :: tmp
    type(octalVector) :: pos, dir
    
    
    write(*,*) "plot_column_density plotting to: ",trim(device)

    ! retriving the address to the root of the tree.
    root => grid%octreeRoot
    
    ! setup the density 
    box_size = REAL(grid%octreeRoot%subcellsize)*2.0
    d = box_size
    off=0.01
    dd= d*(1.0-off)


    ! position of the root node
    xc = real(root%centre%x)
    yc = real(root%centre%y)
    zc = real(root%centre%z)


    ! The transformation matrix TR is used to calculate the world
    ! coordinates of the center of the cell that represents each
    ! array element. The world coordinates of the center of the cell
    ! corresponding to array element f(I,J) are given by
    !        X = TR(1) + TR(2)*I + TR(3)*J
    !        Y = TR(4) + TR(5)*I + TR(6)*J
    
    TR = 0.0
    TR(1) = -d/2.0
    TR(4) = -d/2.0
    TR(2) = d/real(n)
    TR(6) = d/real(n)
    
    
    x1 =  - dd/2.0d0+xc
    y1 =  - dd/2.0d0+yc
    z1 =  - dd/2.0d0+zc
    
    do i = 1, n
       xmap(i) = x1 + dble(i-1)*dd/dble(n-1) + off/2.0
       ymap(i) = y1 + dble(i-1)*dd/dble(n-1) + off/2.0
       zmap(i) = z1 + dble(i-1)*dd/dble(n-1) + off/2.0
    end do
    
    fmin = 1.0e20
    fmax = -1.0e20
    
    if (plane(1:3) == 'x-y') then  
       do j = 1, n
          do i = 1, n
             dir = Vector(0.0d0, 0.0d0, 1.0d0)
             pos = Vector(xmap(i),ymap(j),-dd/2.0d0+zc)
             ! using the function in amr_mod.f90
             tmp = ABS( columnDensity(grid, pos, dir, 1.0) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
!             f(i,j) = tmp
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
          end do
       end do
    elseif (plane(1:3) == 'y-z') then  
       do j = 1, n
          do i = 1, n
             dir = Vector(1.0d0, 0.0d0, 0.0d0)
             pos = Vector(-dd/2.0d0+xc, ymap(i), zmap(j))
             ! using the function in amr_mod.f90
             tmp = ABS( columnDensity(grid, pos, dir, 1.0) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
!             f(i,j) = tmp
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
          end do
       end do
    elseif (plane(1:3) == 'z-x') then  
       do j = 1, n
          do i = 1, n
             dir = Vector(0.0d0, 1.0d0, 0.0d0)
             pos = Vector(xmap(j),-dd/2.0d0+yc, zmap(i))
             ! using the function in amr_mod.f90
             tmp = ABS( columnDensity(grid, pos, dir, 1.0) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
!             f(i,j) = tmp
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
          end do
       end do
    elseif (plane(1:3) == 'x-z') then  
       do j = 1, n
          do i = 1, n
             dir = Vector(0.0d0, -1.0d0, 0.0d0)
             pos = Vector(xmap(i),dd/2.0d0+yc, zmap(j))
             ! using the function in amr_grid.f90
             tmp = ABS( columnDensity(grid, pos, dir, 1.0) )
             if (tmp<=0) tmp=1.0e-36
             f(i,j) = LOG10(tmp)
!             f(i,j) = tmp
             fmin = min(f(i,j),fmin)
             fmax = max(f(i,j),fmax)
          end do
       end do
    end if



    ! Just for safty
    if(fmin == fmax) then
       fmin = fmin-fmin/10.0
       fmax = fmax+fmax/10.0
    end if

    !
    ! FOR DEBUG
    !
!    fmin = -7.0


    ! Open plot device and set up coordinate system. We will plot the
    !image within a unit square.
    !
    IF (PGBEG(0,device,1,1) .NE. 1) STOP
    CALL PGQCOL(CI1, CI2)
    IF (CI2.LT. 15+NCOL) THEN
       WRITE (*,*) 'This program requires a device with at least',&
            15+NCOL,' colors'
       STOP
    END IF
    CALL PGPAGE

    !    CALL PGSCR(0, 0.0, 0.3, 0.2)
    !    CALL PGSVP(0.05,0.95,0.05,0.95)

    call setvp



    do i = 1, n
       xmap(i) = -d/2.0d0 + dble(i-1)*d/dble(n-1)
       ymap(i) = -d/2.0d0 + dble(i-1)*d/dble(n-1)
       zmap(i) = -d/2.0d0 + dble(i-1)*d/dble(n-1)
    end do


    if (plane(1:3) == 'x-y') then
       CALL PGWNAD(real(xmap(1)), real(xmap(n)), real(ymap(1)), real(ymap(n)))
    elseif (plane(1:3) == 'y-z') then
       CALL PGWNAD(real(ymap(1)), real(ymap(n)), real(zmap(1)), real(zmap(n)))
    elseif (plane(1:3) == 'z-x') then
       CALL PGWNAD(real(zmap(1)), real(zmap(n)), real(xmap(1)), real(xmap(n)))
    elseif (plane(1:3) == 'x-z') then
       CALL PGWNAD(real(xmap(1)), real(xmap(n)), real(zmap(1)), real(zmap(n)))
    end if
    
        
    ! Setting the color range and etc..e
    call palett(2, 1.0, 0.5)
    
    ! Plot the image
    CALL PGIMAG(f,N,N, 1, N, 1, N,  fmin, fmax, TR)


    !
    ! Annotation.
    !
    CALL PGSCI(1)

    CALL PGMTXT('t',1.0,0.0,0.0,' ')
    !CALL PGBOX('bcnts',0.0,0,'bcnts',0.0,0)
    CALL PGBOX('bcntsi',0.0,0,'bcntsiv',0.0,0)
    
    if (plane(1:3) == 'x-y') then
       CALL PGMTXT('B',3.0,1.0,1.0,'x [10\u10\dcm]')
       CALL PGMTXT('L',3.0,1.0,1.0,'y [10\u10\dcm]')
    elseif (plane(1:3) == 'y-z') then
       CALL PGMTXT('B',3.0,1.0,1.0,'y [10\u10\dcm]')
       CALL PGMTXT('L',3.0,1.0,1.0,'z [10\u10\dcm]')
    elseif (plane(1:3) == 'z-x') then 
       CALL PGMTXT('B',3.0,1.0,1.0,'z [10\u10\dcm]')
       CALL PGMTXT('L',3.0,1.0,1.0,'x [10\u10\dcm]')
    elseif (plane(1:3) == 'x-z') then 
       CALL PGMTXT('B',3.0,1.0,1.0,'x [10\u10\dcm]')
       CALL PGMTXT('L',3.0,1.0,1.0,'z [10\u10\dcm]')
    end if
    
    ! draw boxies
    call PGSFS(2)  ! we don't want to fill in a box

!    ! this is a recursive routine
!    call draw_rectangle(root, plane)


    !
    ! Draw the markers on the plot with the value of the third dimension.
    ! -- using a routine in this module.    
    call draw_markers(nmarker, xmarker, ymarker, zmarker, &
         plane, value_3rd_dim, width_3rd_dim, show_value_3rd_dim)




    ! routines in pixplot_module    
    
    ! draw wedge
    offset = 0
    CALL PGWEDG('BI', 4.0, 5.0, FMIN, FMAX+offset, "LOG10( [g/cm\u3\d] )")
    CALL PGSCH(0.6)

!    CALL FIDDLE
    CALL PGASK(.FALSE.)  

    call PGEND
    
    ! deallocate(a_root)
  end subroutine plot_column_density




  ! recursively draws a rectangle by calliing a pgplot rouitne
  ! must be only used in draw_cells_on_density routine in this module.
  recursive subroutine draw_rectangle(thisOctal, plane, val_3rd_dim)
    implicit none
    type(octal), intent(in)           :: thisOctal
    character(LEN=*), intent(in)      :: plane
    ! The value of the third dimension.
    ! For example, if plane = "x-y", the third dimension is the value of z.
    ! Then, when  val_3rd_dim = 0.0, this will plot the density (and the grid)
    ! on the z=0 plane, and so on...
    real, optional, intent(in)        :: val_3rd_dim 
    !
    type(octal), pointer  :: child 
    type(octalvector) :: rvec
    integer :: subcell, i
    !real :: t
    
    real  :: x,y,z,h,w, eps
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call draw_rectangle(child, plane)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          h = thisOctal%subcellSize  ! height
          w = thisOctal%subcellSize  ! width
          if (.not.thisOctal%cylindrical) then
             x = rVec%x; y = rVec%y; z = rVec%z
          else
             x  = sign(sqrt(rVec%x**2 + rvec%y**2), rvec%x)
             y  = sign(sqrt(rVec%x**2 + rvec%y**2), rvec%y)
             z  = rVec%z
          endif

          if (present(val_3rd_dim)) then
!            if ( plane(1:3) == "x-y" .and.  &
!               rVec%z > (val_3rd_dim-h) .and.  rVec%z < (val_3rd_dim+h)) then
!               call PGRECT(x-w/2.0, x+w/2.0, y-h/2.0, y+h/2.0)
!            elseif ( plane(1:3) == "y-z" .and. &
!               rVec%x > (val_3rd_dim-h) .and.  rVec%x < (val_3rd_dim+h)) then
!               call PGRECT(y-w/2.0, y+w/2.0, z-h/2.0, z+h/2.0)
!            elseif ( plane(1:3) == "z-x" .and.  &
!               rVec%y > (val_3rd_dim-h) .and.  rVec%y < (val_3rd_dim+h)) then
!               call PGRECT(z-w/2.0, z+w/2.0, x-h/2.0, x+h/2.0)
!            elseif ( plane(1:3) == "x-z" .and.  &
!               rVec%y > (val_3rd_dim-h) .and.  rVec%y < (val_3rd_dim+h)) then
!               call PGRECT(x-w/2.0, x+w/2.0, z-h/2.0, z+h/2.0)
!            end if
             
             eps = h/100.0
             if ((rVec%z > 0.d0).and. plane(1:3) == "x-y" &
                  .and. (ABS(rVec%z-val_3rd_dim)-h/2.0) < eps) then
                if (.not.thisOctal%cylindrical) then
                   call PGRECT(x-w/2.0, x+w/2.0, y-h/2.0, y+h/2.0)
                else
                   call draw_segment(thisOctal, subcell, .false.)
                endif
             elseif ( plane(1:3) == "y-z" .and. (ABS(rVec%x-val_3rd_dim)-h/2.0) < eps) then 
                call PGRECT(y-w/2.0, y+w/2.0, z-h/2.0, z+h/2.0)
             elseif ( plane(1:3) == "z-x" .and. (ABS(rVec%y-val_3rd_dim)-h/2.0) < eps) then 
                call PGRECT(z-w/2.0, z+w/2.0, x-h/2.0, x+h/2.0)
             elseif ( plane(1:3) == "x-z" .and. (ABS(rVec%y-val_3rd_dim)-h/2.0) < eps) then 
                call PGRECT(x-w/2.0, x+w/2.0, z-h/2.0, z+h/2.0)
             end if
          else
             if (plane(1:3) == 'x-y' .and.  ABS(rVec%z) < h/2.d0 ) then
                if (.not.thisOctal%cylindrical) then
                   call PGRECT(x-w/2.0, x+w/2.0, y-h/2.0, y+h/2.0)
                else
                   call draw_segment(thisOctal, subcell, .false.)
                endif
             elseif (plane(1:3) == 'y-z' .and. ABS(rVec%x) < h ) then
                call PGRECT(y-w/2.0, y+w/2.0, z-h/2.0, z+h/2.0)
             elseif (plane(1:3) == 'z-x' .and. ABS(rVec%y) < h ) then
                call PGRECT(z-w/2.0, z+w/2.0, x-h/2.0, x+h/2.0)
             elseif (plane(1:3) == 'x-z' .and. ABS(rVec%y) < h ) then
                call PGRECT(x-w/2.0, x+w/2.0, z-h/2.0, z+h/2.0)
             end if
          end if

          
       endif
    enddo
  end subroutine draw_rectangle


  subroutine draw_segment(thisOctal, subcell, filled)
    type(OCTAL) :: thisOctal
    integer :: subcell
    real :: phi1,phi2,r1,r2,ang
    type(OCTALVECTOR) :: rVec
    integer :: i
    integer, parameter :: nAng = 100
    integer ::npts
    logical :: filled
    real :: x(2*nAng+10), y(2*nAng+10)
    real :: x1, y1

    rVec = subcellCentre(thisOctal, subcell)
    r1 = sqrt(rVec%x**2+rVec%y**2)-thisOctal%subcellsize/2.d0
    r2 = sqrt(rVec%x**2+rVec%y**2)+thisOctal%subcellsize/2.d0


    phi1 = atan2(rvec%y,rvec%x)
    phi1 = phi1 - returndPhi(thisOctal)
    phi2 = atan2(rvec%y,rvec%x)
    phi2 = phi2 + returndPhi(thisOctal)

!    if (thisOctal%splitAzimuthally) then
!       phi1 = atan2(rvec%y,rvec%x)
!       if (phi1  < 0.d0) phi1 = phi1 + twoPi
!       phi1 = phi1 - thisOctal%dPhi/4.d0
!       phi2 = atan2(rvec%y,rvec%x)
!       if (phi2  < 0.d0) phi2 = phi2 + twoPi
!       phi2 = phi2 + thisOctal%dPhi/4.d0
!    else
!       phi1 = atan2(rvec%y,rvec%x) - thisOctal%dPhi/2.d0
!       phi2 = atan2(rvec%y,rvec%x) + thisOctal%dPhi/2.d0
!    endif
    npts = 0

    npts = npts + 1
    x(npts) = r1 * cos(phi1)
    y(npts) = r1 * sin(phi1)
    npts = npts + 1
    x(npts) = r2 * cos(phi1)
    y(npts) = r2 * sin(phi1)
    do i = 1, nang
       ang = phi1+(phi2-phi1)*real(i-1)/real(nAng-1)
       npts = npts + 1
       x(npts) = r2 * cos(ang)
       y(npts) = r2 * sin(ang)
    enddo
    npts = npts + 1
    x(npts) = r1 * cos(phi2)
    y(npts) = r1 * sin(phi2)
    do i = 1, nAng
       ang = phi2+(phi1-phi2)*real(i-1)/real(nAng-1)
       npts = npts + 1
       x(npts) = r1 * cos(ang)
       y(npts) = r1 * sin(ang)
    enddo
    if (filled) then
       call pgpoly(npts,x, y)
    else
       call pgline(npts, x, y)
    endif


!    x1  = r1 * cos(thisOctal%phi)
!    y1  = r1 * sin(thisOctal%phi)
!    call pgqci(i)
!    call pgmove(x1,y1)
!    call pgsci(2)
!    x1  = r2 * cos(thisOctal%phi)
!    y1  = r2 * sin(thisOctal%phi)
!    call pgdraw(x1,y1)
!    call pgsci(i)



  end subroutine draw_segment




  !
  ! This overlays the AMR grid values on the grid structure using PGPLOT
  ! routines. 
  !
  subroutine plot_AMR_values(grid, name, plane, value_3rd_dim,  device, logscale, withgrid, &
       nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim, boxfac, ilam, xStart, &
       xEnd, yStart, yEnd, fixValMin,fixValMax)
    implicit none
    type(gridtype), intent(in) :: grid
    character(len=*), intent(in)  :: name   ! "rho", "temperature", chiLine", "etaLine", 
    !                                       ! "etaCont", "Vx", "Vy" or "Vz"
    character(len=*), intent(in)  :: plane  ! must be 'x-y', 'y-z', 'z-x' or' x-z'
    character(len=10) :: thisPlane
    real(double) , optional :: fixValMin, fixValMax
    ! The value of the third dimension.
    ! For example, if plane = "x-y", the third dimension is the value of z.
    ! Then, when  valuel_3rd_dim = 0.0, this will plot the density (and the grid)
    ! on the z=0 plane, and so on...
    real, intent(in)              :: value_3rd_dim 
    character(len=*), intent(in)  :: device ! PGPLOT plotting device
    logical, intent(in)           :: logscale ! logscale if T, linear if not
    logical, intent(in)           :: withgrid ! plot
    real, intent(in),optional     :: boxfac
    real, intent(in), optional :: xStart, xEnd, yStart, yEnd
    integer, intent(in),optional     :: ilam
    ! To put the markers to indicate the important points in the plots.
    integer, intent(in),optional :: nmarker           ! number of markers
    real, intent(in),optional    :: xmarker(:)  ! position of x
    real, intent(in),optional    :: ymarker(:)  ! position of y
    real, intent(in),optional    :: zmarker(:)  ! position of z
    real, intent(in),optional    :: width_3rd_dim     ! Use this to restrict the markers to be plotted..  
    logical, intent(in),optional :: show_value_3rd_dim! If T, the value of the third dimension will be shown.

    !
    type(octal), pointer :: root

    !integer :: thisSubcell, oldSubcell
    real :: d
    real :: xc, yc , zc
    real(double) :: valuemax, valuemin
    
    ! For plotting density
    integer :: pgbeg
    integer :: n, ncol, nlev
    parameter (n=1500, ncol=32, nlev=10)
    integer :: ci1,ci2, ilo, ihi 
        
    real :: box_size, fac
    character(LEN=30) :: char_val
    character(LEN=30) :: label_wd
    integer, parameter :: luout = 26
    integer :: iPlane, nPlanes
    character(len=80) :: message
    character(LEN=50) :: filename_prof
    real   :: v3

    
    call writeInfo("plot_AMR_values plotting to: "//trim(device), TRIVIAL)

    ! retriving the address to the root of the tree.
    root => grid%octreeRoot
    
    !
    box_size = REAL(grid%octreeRoot%subcellsize)*2.0
    if (grid%octreeRoot%cylindrical) box_size = box_size * 2.0
    d = box_size
    if (PRESENT(boxfac)) d = d * boxfac

    ! position of the root node
    xc = real(root%centre%x)
    yc = real(root%centre%y)
    zc = real(root%centre%z)


!    if (name.eq."temperature") then
!       valueMax = 1500.
!       valueMin = 100.
!    endif
    
!    !=================================================
!    ! Just for hardwaring manupulations.
!    ! Comment the lines below after a use.
!    valueMax = 10.0**1.6
!    valueMin = 10.0**0.9
    
    ! Open plot device and set up coordinate system. We will plot the
    !image within a unit square.
    !

    if (grid%octreeRoot%twod) then
       nPlanes = 1
    else
       nPlanes = 3
    endif

    if (nPlanes == 1) then
       IF (PGBEG(0,trim(device),1,1) .NE. 1) STOP
    else
       IF (PGBEG(0,trim(device),3,1) .NE. 1) STOP
    endif

    do iPlane = 1, nPlanes

       if (nPlanes == 1) then
           thisPlane = plane
       else
          select case(iPlane)
             case(1)
                thisPlane = "x-y"
             case(2)
                thisPlane = "x-z"
             case(3)
                thisPlane = "y-z"
          end select
       endif

    ! Finding the plotting range
    valueMin = 1.d30
    valueMax = -1.d30
    call minMaxValue(grid%octreeRoot, name, plane, ValueMin, ValueMax, grid, ilam)
        
    if (PRESENT(fixValMax)) valueMax = fixValMax
    if (PRESENT(fixValMin)) valueMin = fixValMin


    if (logscale) then
       if (valueMax<=0) valueMax=tiny(valueMax)
       if (valueMin<=0) valueMin=tiny(valueMin)
       if ((log10(valueMax)-log10(valueMin)) > 20.d0) valueMin = 1.d-20 * valueMax
!       if (name.eq."rho") valuemin = valuemax * 1.e-10
       write(message,*) "Value Range: ",log10(valueMin), " -- ", log10(valueMax)
       call writeInfo(message,TRIVIAL)
    else
       write(message,*) "Value Range: ",valueMin, " -- ", valueMax
       call writeInfo(message,TRIVIAL)
    end if

    !
    ! For safty
    if (valueMin == valueMax) then
       if (valueMax /= 0.0) then
          valueMax = valueMax + valueMax/10.0
       else
          valueMax = valueMax + 1.0e6
       end if
    end if

 
       if (iplane == 1) then
          CALL PGQCOL(CI1, CI2)
          IF (CI2.LT. 15+NCOL) THEN
             WRITE (*,*) 'This program requires a device with at least',&
                  15+NCOL,' colors'
             STOP
          END IF
       endif


    CALL PGPAGE

    !    CALL PGSCR(0, 0.0, 0.3, 0.2)
    !    CALL PGSVP(0.05,0.95,0.05,0.95)

    call setvp

    ! Setting the plot boundaries.

    if (.not.grid%octreeRoot%cylindrical) then
       select case(thisplane)
       case("x-y")
          CALL PGWNAD(-d/2.0+xc, d/2.0+xc, -d/2.0+yc, d/2.0+yc)
          v3 = root%centre%z
       case("y-z")
          CALL PGWNAD(-d/2.0+yc, d/2.0+yc, -d/2.0+zc, d/2.0+zc)
          v3 = root%centre%x
       case("z-x")
          CALL PGWNAD(-d/2.0+zc, d/2.0+zc, -d/2.0+xc, d/2.0+xc)
          v3 = root%centre%y
       case("x-z")
          CALL PGWNAD(-d/2.0+xc, d/2.0+xc, -d/2.0+zc, d/2.0+zc)
          if (present(boxfac).and.grid%octreeRoot%twod) then
             CALL PGWNAD(0., d, -d/2.0, d/2.0)
          endif
          v3 = root%centre%y
          if (present(xStart).and.present(xEnd).and.present(yStart).and.present(yEnd)) then
             call pgwnad(xStart, xEnd, yStart, yEnd)
          endif
       end select
    else
       select case(thisplane)
       case("x-y")
          CALL PGWNAD(-d/2.0, d/2.0, -d/2.0, d/2.0)
          v3 = root%centre%z
       case("y-z")
          CALL PGWNAD(-d/2.0, d/2.0, -d/2.0, d/2.0)
          v3 = root%centre%x
       case("z-x")
          CALL PGWNAD(-d/2.0, d/2.0, -d/2.0, d/2.0)
          v3 = root%centre%y
       case("x-z")
          CALL PGWNAD(-d/2.0, d/2.0, -d/2.0, d/2.0)
          if (present(boxfac).and.grid%octreeRoot%twod) then
             CALL PGWNAD(0., d, -d/2.0, d/2.0)
          endif
          v3 = root%centre%y
          if (present(xStart).and.present(xEnd).and.present(yStart).and.present(yEnd)) then
             call pgwnad(xStart, xEnd, yStart, yEnd)
          endif
       end select
    endif





!    if (iPlane == 1) then
       call palett(2, 1.0, 0.5)
       call pgqcir(ilo, ihi)
!    endif


    !
    ! Calling a recursive function in this module
    call plot_values(root, name, thisplane, value_3rd_dim, logscale, valueMax, valueMin, &
         ilo, ihi, grid, ilam, withgrid)


! Comment out the following if you want the pix value as a function of distance from 
! the center of the gird. Useful for some cases.

!    ! writing the radial profile of the values on the specified plane.
!    filename_prof = 'profile_'//TRIM(ADJUSTL(name))//'_'//TRIM(ADJUSTL(plane))//'.dat' 
!    open(unit=luout, file = TRIM(ADJUSTL(filename_prof)), status='replace')
!!    write(luout, '(a, 2x, 1PE18.4)') '#  The 3rd dimension value = ', value_3rd_dim
!!    write(luout, '(a)') '#  format :  distance [10^10cm]  --  values [?]'  
!
!    call radial_profile(root, name, plane, v3, luout, root%centre, grid)
!    close(luout)
    !
    ! Annotation.
    !
    CALL PGSCI(1)
!    CALL PGMTXT('t',1.0,0.5,0.0, TRIM(name))
    !CALL PGBOX('bcnts',0.0,0,'bcnts',0.0,0)
    CALL PGBOX('bcntsi',0.0,0,'bcntsiv',0.0,0)

    write(char_val, '(1PE9.1)') value_3rd_dim
    
    if (thisplane(1:3) == 'x-y') then
       CALL PGMTXT('B',3.0,1.0,1.0,'X [10\u10\dcm]')
       CALL PGMTXT('L',3.0,1.0,1.0,'Y [10\u10\dcm]')
!       CALL PGMTXT('T',1.0,0.0,0.0, 'Z='//TRIM(char_val))
    elseif (thisplane(1:3) == 'y-z') then
       CALL PGMTXT('B',3.0,1.0,1.0,'Y [10\u10\dcm]')
       CALL PGMTXT('L',3.0,1.0,1.0,'Z [10\u10\dcm]')
!       CALL PGMTXT('T',1.0,0.0,0.0, 'X='//TRIM(char_val))
    elseif (thisplane(1:3) == 'z-x') then 
       CALL PGMTXT('B',3.0,1.0,1.0,'Z [10\u10\dcm]')
       CALL PGMTXT('L',3.0,1.0,1.0,'X [10\u10\dcm]')
!       CALL PGMTXT('T',1.0,0.0,0.0, 'Y='//TRIM(char_val))
    elseif (thisplane(1:3) == 'x-z') then 
       CALL PGMTXT('B',3.0,1.0,1.0,'X [10\u10\dcm]')
       CALL PGMTXT('L',3.0,1.0,1.0,'Z [10\u10\dcm]')
!       CALL PGMTXT('T',1.0,0.0,0.0, 'Y='//TRIM(char_val))
    end if

    
    !
    ! Draw grids..
    !
!    if (withgrid) then
!       ! draw boxies
!       call PGSFS(2)  ! we don't want to fill in a box this time
!       call PGSCI(1) ! changing the color index.
!       ! this is a recursive routine in this module
!       call draw_rectangle(root, thisplane, value_3rd_dim)
!       call PGSCI(1) ! change color index to default.
!       call PGSFS(1)
!    end if


    
    !
    ! Draw the markers on the plot with the value of the third dimension.
    ! -- using a routine in this module.    
    if (present(nmarker)) then
       call draw_markers(nmarker, xmarker, ymarker, zmarker, &
            thisplane, value_3rd_dim, width_3rd_dim, show_value_3rd_dim)
    endif


    !
    ! draw wedge
    !

    select case (name)
    case ("rho")
       if (logscale) then
         label_wd = "LOG10( [g/cm\u3\d] )"
       else
         label_wd = "[g/cm\u3\d]"
       end if
    case ("temperature")
       if (logscale) then
         label_wd = "LOG10( [Kelvin] )"
       else
         label_wd = "[Kelvin]"
       end if
    case default
       label_wd = "Pixel Value"
    end select
       

    if(logscale) then
       CALL PGWEDG('BI', 4.0, 5.0, real(log10(valueMin)), real(log10(valueMax)), TRIM(ADJUSTL(label_wd)) )
    else
       CALL PGWEDG('BI', 4.0, 5.0, real(valueMin), real(valueMax), TRIM(ADJUSTL(label_wd)))
    end if
    
!    CALL PGSCH(0.6)

!    CALL FIDDLE
    CALL PGASK(.FALSE.)  

 enddo
    call PGEND

    
    ! deallocate(a_root)
  end subroutine plot_AMR_values


  !
  ! This is modefied from plotZplane
  !
  recursive subroutine plot_values(thisOctal, name, plane, val_3rd_dim, &
       logscale, valueMax, valueMin, ilo, ihi, grid, ilam, withgrid)
    implicit none
    !
    type(octal), pointer   :: thisOctal
    character(len=*), intent(in)  :: name     ! "rho", "temperature", chiLine", "etaLine",  
    !                                         ! "etaCont", "Vx", "Vy" or "Vz", "dV_dR"
    character(len=*), intent(in)  :: plane    ! must be 'x-y', 'y-z', 'z-x', 'x-z'
    ! The value of the third dimension.
    ! For example, if plane = "x-y", the third dimension is the value of z.
    ! Then, when  val_3rd_dim = 0.0, this will plot the density (and the grid)
    ! on the z=0 plane, and so on...
    real, intent(in)              :: val_3rd_dim 
    logical, intent(in)           :: logscale ! logscale if T, linear if not
    real(double), intent(in) :: valueMin, valueMax
    integer, intent(in) :: ilo, ihi
    logical :: withgrid
    integer, intent(in),optional :: ilam
    !
    !
    type(octal), pointer  :: child 
    type(octalvector) :: rvec, rhat
    type(gridtype) :: grid
    real(double) :: value
    real(double) :: kabs,ksca
    TYPE(vector)   :: Velocity

    integer :: subcell, i, idx
    real :: t
    real :: xp, yp, xm, ym, zp, zm
    real(double) :: d, L, eps, phi, phiFromplane, dphi
    logical :: plot_this_subcell
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call plot_values(child, name, plane, val_3rd_dim, &
                     logscale, valueMax, valueMin, ilo, ihi, grid, ilam, withgrid)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          L = thisOctal%subcellSize
          d = L/2.0d0
!          eps = d/100.0d0
          eps = l/2.d0
!          eps = d/5.0d0
!          eps = 0.0d0
          if (.not.thisOctal%cylindrical) then
             xp = REAL(rVec%x + d)
             xm = REAL(rVec%x - d)
             yp = REAL(rVec%y + d)
             ym = REAL(rVec%y - d)
          else
             xp = sign(sqrt(rVec%x**2+rVec%y**2) + d,rVec%x)
             xm = sign(sqrt(rVec%x**2+rVec%y**2) - d,rVec%x)
             yp = sign(sqrt(rVec%x**2+rVec%y**2) + d,rVec%y)
             ym = sign(sqrt(rVec%x**2+rVec%y**2) - d,rVec%y)
          endif
          zp = REAL(rVec%z + d)
          zm = REAL(rVec%z - d)     
          
          plot_this_subcell = .false.
          
          if ( plane(1:3) == "x-y" .and. ((ABS(rVec%z-val_3rd_dim)-d) < eps).and.(rVec%z>0.d0)) then
!             write(*,*) rVec%z,val_3rd_dim,abs(rVec%z-val_3rd_dim)-d,eps
             plot_this_subcell = .true.      
          elseif ( plane(1:3) == "y-z" .and. (ABS(rVec%x-val_3rd_dim)-d) < eps) then 
             plot_this_subcell = .true.
          elseif ( plane(1:3) == "z-x" .and. (ABS(rVec%y-val_3rd_dim)-d) < eps) then 
             plot_this_subcell = .true.
          elseif ( plane(1:3) == "x-z" .and. (ABS(rVec%y-val_3rd_dim)-d) < eps) then 
             plot_this_subcell = .true.
          else
             plot_this_subcell = .false.
          end if

!          if (plane(1:3)=="x-y") plot_this_subcell = .true.

          if (thisOctal%cylindrical.and.(plane.ne."x-y")) then
             plot_this_subcell = .false.
             phi = atan2(rVec%y,rVec%x)
             if (phi < 0.) phi = phi + twoPi
             dphi = returndPhi(thisOctal)
!             if (thisOctal%splitAzimuthally) then
!                dphi = thisOctal%dPhi / 4.d0
!             else
!                dphi = thisOctal%dPhi / 2.d0
!             endif
             select case (plane)
                case("y-z")
                   phiFromPlane = min(abs(pi/2.d0-phi),abs(3.d0*pi/2.d0-phi))
                case("z-x","x-z")
                   phiFromPlane = min(phi,abs(phi-pi))
             end select
             if (phiFromPlane <= 2.*dphi) then
                plot_this_subcell = .true.
             endif
          endif
!          if ( plane(1:3) == "x-y" .and.  &
!               rVec%z > (val_3rd_dim-L) .and.  rVec%z < (val_3rd_dim+L)) then
!             plot_this_subcell = .true.      
!          elseif ( plane(1:3) == "y-z" .and. &
!               rVec%x > (val_3rd_dim-L) .and.  rVec%x < (val_3rd_dim+L)) then
!             plot_this_subcell = .true.
!          elseif ( plane(1:3) == "z-x" .and.  &
!               rVec%y > (val_3rd_dim-L) .and.  rVec%y < (val_3rd_dim+L)) then
!             plot_this_subcell = .true.
!          elseif ( plane(1:3) == "x-z" .and.  &
!             rVec%y > (val_3rd_dim-L) .and.  rVec%y < (val_3rd_dim+L)) then
!             plot_this_subcell = .true.
!          else
!             plot_this_subcell = .false.
!          end if

          
          if (plot_this_subcell) then
             select case (name)
             case("rho")
                value = thisOctal%rho(subcell)
                if (thisOctal%diffusionApprox(subcell)) value = 1.
             case("ionization")
                value = thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon))
             case("photocoeff")
                value = thisOctal%photoIonCoeff(subcell,1)
             case("temperature")
                value = thisOctal%temperature(subcell)
             case("chiLine")
                value = thisOctal%chiLine(subcell)
             case("etaLine")
                value = thisOctal%etaLine(subcell)
             case("etaCont")
                value = thisOctal%etaCont(subcell)
             case("dusttype")
                value = thisOctal%dustTypeFraction(subcell,1)
             case("dusttype1")
                value = thisOctal%dustTypeFraction(subcell,1)
             case("dusttype2")
                value = thisOctal%dustTypeFraction(subcell,2)
             case("dusttype3")
                value = thisOctal%dustTypeFraction(subcell,3)
             case("dusttype4")
                value = thisOctal%dustTypeFraction(subcell,4)
             case("crossings")
                value = thisOctal%nCrossings(subcell)
                if (thisOctal%diffusionApprox(subcell)) value = 1.e6
             case("direct")
                if (thisOctal%ncrossings(subcell) > 0) then
                   value = real(thisOctal%nDirectPhotons(subcell)) / real(thisOctal%nCrossings(subcell))
                else
                   value = 0.
                endif

             case("Vx")
!                value = thisOctal%velocity(subcell)%x * cSpeed/1.0d5 ![km/s]
                velocity = amrGridVelocity(thisOctal, rvec)
                value = velocity%x * cSpeed/1.0d5 ![km/s]
             case("Vy")
!                value = thisOctal%velocity(subcell)%y * cSpeed/1.0d5 ![km/s]
                velocity = amrGridVelocity(thisOctal, rvec)
                value = velocity%y * cSpeed/1.0d5 ![km/s]
             case("Vz")
!                value = thisOctal%velocity(subcell)%z * cSpeed/1.0d5 ![km/s]
                velocity = amrGridVelocity(thisOctal, rvec)
                value = velocity%z * cSpeed/1.0d5 ![km/s]
             case("dV_dR")
                ! The directional derivative the velocity field in radial direction.
                rhat = rvec
                call normalize(rhat)
                value = amrGridDirectionalDeriv(grid,rvec, o2s(rhat), thisOctal) * (cSpeed_dbl/1.0d5)  ![km/s]
             case("tau")
                call returnKappa(grid, thisOctal, subcell, rosselandKappa=kabs)
                value = thisOctal%subcellsize * kabs * thisOctal%rho(subcell) * 1.e10
                if (thisOctal%diffusionApprox(subcell)) value = 1.e10

             case default
                write(*,*) "Error:: unknown name passed to grid_mod::plot_values."
                stop
             end select

             ! for safety
             if (logscale .and. value <=0 ) value = 1.0e-30

             if (logscale) then            
                t = (log10(value)-log10(valueMin))/(log10(valueMax)-log10(valueMin))
             else ! linear
                t = (value-valueMin)/(valueMax-valueMin)
             end if
             
             if (t < 0.) t = 0.
             if (t > 1.) t = 1.
             idx = int(t * real(ihi - ilo) + real(ilo))

             if (.not.thisOctal%inFlow(subcell)) idx = ilo

             ! set color
             call pgsci(idx)
             

             select case (plane)
             case ("x-y")
                if (.not.thisOctal%cylindrical) then
                   call pgrect(xm, xp, ym, yp)
                else
                   call draw_segment(thisOctal, subcell, .true.)
                   if (withgrid) then
                      call pgqci(i)
                      call PGSFS(2)  ! we don't want to fill in a box this time
                      call PGSCI(1) ! changing the color index.
                      call draw_segment(thisOctal, subcell, .false.)
                      call PGSFS(1)
                      call pgsci(i)
                   endif
                endif
             case ("y-z")
                call pgrect(ym, yp, zm, zp)
                if (withgrid) then
                   call pgqci(i)
                   call PGSFS(2)  ! we don't want to fill in a box this time
                   call PGSCI(1) ! changing the color index.
                   call pgrect(ym, yp, zm, zp)
                   call PGSCI(1) ! change color index to default.
                   call PGSFS(1)
                   call pgsci(i)
                endif
                if ((thisOctal%cylindrical).and.(thisOctal%dphi > 1.001d0*pi)) then
                   ym = -ym
                   yp = -yp
                   call pgrect(ym, yp, zm, zp)
                   if (withgrid) then
                      call pgqci(i)
                      call PGSFS(2)  ! we don't want to fill in a box this time
                      call PGSCI(1) ! changing the color index.
                      call pgrect(ym, yp, zm, zp)
                      call PGSCI(1) ! change color index to default.
                      call PGSFS(1)
                      call pgsci(i)
                   endif
                endif

             case ("z-x")
                call pgrect(zm, zp, xm, xp)
                if (withgrid) then
                      call pgqci(i)
                   call PGSFS(2)  ! we don't want to fill in a box this time
                   call PGSCI(1) ! changing the color index.
                   call pgrect(zm, zp, xm, xp)
                   call PGSFS(1)
                   call pgsci(i)
                endif
                if ((thisOctal%cylindrical).and.(thisOctal%dphi > 1.001d0*pi)) then
                   xm = -xm
                   xp = -xp
                   call pgrect(zm, zp, xm, xp)
                   if (withgrid) then
                      call pgqci(i)
                      call PGSFS(2)  ! we don't want to fill in a box this time
                      call PGSCI(1) ! changing the color index.
                      call pgrect(zm, zp, xm, xp)
                      call PGSCI(i) ! change color index to default.
                      call PGSFS(1)
                   endif
                endif
             case ("x-z")
                call pgrect(xm, xp, zm, zp)
                if (withgrid) then
                   call pgqci(i)
                   call PGSFS(2)  ! we don't want to fill in a box this time
                   call PGSCI(1) ! changing the color index.
                   call pgrect(xm, xp, zm, zp)
                   call PGSCI(i) ! change color index to default.
                   call PGSFS(1)
                endif
                if ((thisOctal%cylindrical).and.(thisOctal%dphi > 1.001d0*pi)) then
                   xm = -xm
                   xp = -xp
                   call pgrect(xm, xp, zm, zp)
                   if (withgrid) then
                      call pgqci(i)
                      call PGSFS(2)  ! we don't want to fill in a box this time
                      call PGSCI(1) ! changing the color index.
                      call pgrect(xm, xp, zm, zp)
                      call PGSCI(i) ! change color index to default.
                      call PGSFS(1)
                   endif
                endif


             end select

          end if
       end if

    end do

  end subroutine plot_values


  
  recursive subroutine minMaxDensity(thisOctal, rhoMin, rhoMax)

  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: rhoMin, rhoMax
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call minMaxDensity(child, rhoMin, rhoMax)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell)) then
             if (thisOctal%rho(subcell) > rhoMax) rhoMax = thisOctal%rho(subcell)
             if (thisOctal%rho(subcell) < rhoMin) rhoMin = thisOctal%rho(subcell)
          endif
       endif
    enddo

  end subroutine minMaxDensity



  !
  !  Recursively finds the min and max values of "name"d value.
  !
  recursive subroutine minMaxValue(thisOctal, name, plane, valueMin, valueMax, grid, ilam)
    implicit none
    type(octal), pointer   :: thisOctal
    type(gridtype) :: grid
    character(LEN=*), intent(in) :: name ! See the options in plot_AMR_values
    character(len=*), intent(in)  :: plane    ! must be 'x-y', 'y-z' or 'z-x' 
    real(double), intent(inout) :: valueMin, valueMax
    !
    type(octal), pointer  :: child 
    integer :: subcell, i
    integer, intent(in), optional :: ilam
    real(double) :: value
    !
    real :: xp, yp, xm, ym, zp, zm
    real(double) :: ksca, kabs
    real(double) :: d
    logical :: use_this_subcell, update
    type(octalvector) :: rvec, rhat, velocity
    
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call minMaxValue(child, name, plane, valueMin, valueMax, grid, ilam)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          d = thisOctal%subcellSize/2.d0
          xp = REAL(rVec%x + d)
          xm = REAL(rVec%x - d)
          yp = REAL(rVec%y + d)
          ym = REAL(rVec%y - d)
          zp = REAL(rVec%z + d)
          zm = REAL(rVec%z - d)          

          use_this_subcell = .false.
          if ( plane(1:3) == "x-y" .and. ABS(rVec%z) < d*2.0d0) then
             use_this_subcell = .true.
          elseif ( plane(1:3) == "y-z" .and. ABS(rVec%x) < d*2.0d0) then
             use_this_subcell = .true.
          elseif ( plane(1:3) == "z-x" .and. ABS(rVec%y) < d*2.0d0) then
             use_this_subcell = .true.
          else
             use_this_subcell = .false.
          end if

          use_this_subcell = .true.

          if (use_this_subcell) then
             update =.true.
             select case (name)
             case("rho")
                value = thisOctal%rho(subcell)
                if (value < 1.e-25) update = .false.
             case("ionization")
                value = thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon))
             case("photocoeff")
                value = thisOctal%photoIonCoeff(subcell,1)
             case("temperature")
                value = thisOctal%temperature(subcell)
                if (value < 3.)  update = .false.
             case("dusttype")
                value = thisOctal%dustTypeFraction(subcell,1)
             case("dusttype1")
                value = thisOctal%dustTypeFraction(subcell,1)
             case("dusttype2")
                value = thisOctal%dustTypeFraction(subcell,2)
             case("dusttype3")
                value = thisOctal%dustTypeFraction(subcell,3)
             case("dusttype4")
                value = thisOctal%dustTypeFraction(subcell,4)
             case("chiLine")
                value = thisOctal%chiLine(subcell)
             case("etaLine")
                value = thisOctal%etaLine(subcell)
             case("etaCont")
                value = thisOctal%etaCont(subcell)
             case("crossings")
                value = thisOctal%ncrossings(subcell)
                if (thisOctal%diffusionApprox(subcell)) value = 1.e6
             case("direct")
                if (thisOctal%ncrossings(subcell) > 0) then
                   value = real(thisOctal%nDirectPhotons(subcell)) / real(thisOctal%nCrossings(subcell))
                else
                   value = 0.
                endif
             case("Vx")
                value = thisOctal%velocity(subcell)%x * cSpeed/1.0d5 ![km/s]
             case("Vy")
                value = thisOctal%velocity(subcell)%y * cSpeed/1.0d5 ![km/s]
             case("Vz")
                value = thisOctal%velocity(subcell)%z * cSpeed/1.0d5 ![km/s]
             case("dV_dR")
                ! The directional derivative the velocity field in radial direction.
                rhat = rvec
                call normalize(rhat)
                value = amrGridDirectionalDeriv(grid,rvec, o2s(rhat), thisOctal) * (cSpeed_dbl/1.0d5)  ![km/s]
             case("tau")
                call returnKappa(grid, thisOctal, subcell, rosselandKappa=kabs)
                value = thisOctal%subcellsize * kabs * thisOctal%rho(subcell) * 1.e10
                if (thisOctal%diffusionApprox(subcell)) update = .false.

             case default
                write(*,*) "Error:: unknow name passed to MinMaxValue.",trim(name)
                stop
             end select
             if (update) then
                valueMax = MAX(valueMax, value)
                valueMin = Min(valueMin, value)
             endif
          end if
       end if
    enddo

  end subroutine minMaxValue


  !
  ! Print out the basic infomation about the grid.
  !
  
  ! if filename is '*' then it prints on screen.
  subroutine grid_info(thisGrid, filename)
    implicit none
    type(gridtype), intent(in) :: thisGrid
    character(LEN=*), intent(in) :: filename
    integer :: UN
    real(double) :: tmp,fac
    integer :: nOctals,nVoxels
    
    if (filename(1:1) == '*') then
       UN = 6   ! prints on screen
    else
       UN = 69
       open(unit=UN, file = TRIM(filename), status = 'replace',form='formatted')
    end if


    if (thisGrid%adaptive) then
       nOctals=0; nVoxels=0
       call countVoxels(thisGrid%octreeRoot,nOctals,nVoxels)
    
       write(UN,'(a)') ' '
       write(UN,'(a)') '######################################################'
       write(UN,'(a)') 'Grid info :'
       write(UN,'(a)') ' '
       write(UN,*)     'geometry             = ', thisGrid%geometry
       write(UN,*)     'maxDepth             = ', thisGrid%maxDepth
       write(UN,*)     'halfSmallestSubcell  = ', thisGrid%halfSmallestSubcell, ' [10^10 cm]'
       write(UN,*)     'nOctals              = ', nOctals
       write(UN,*)     'nVoxels              = ', nVoxels
       write(UN,*)     'smoothingFactor      = ', thisGrid%smoothingFactor
       write(UN,*)     'grid center          =', thisGrid%octreeRoot%centre
       write(UN,*)     'Size of largest cell =', thisGrid%octreeRoot%subcellSize*2.0, ' [10^10 cm]'
       write(UN,'(a)') '#######################################################'
       write(UN,'(a)') ' '
       write(Un,'(a)') ' '
    
       fac = 1.d0/(2.d0**dble(thisGrid%maxDepth))
       if (fac < epsilon(1.d0)) then
          write(UN,'(a)') "**** WARNING: Grid cell depth is so great numerical problems may occur****"
       endif
       write(*,*) fac,epsilon(1.d0)
    else
       write(UN,'(a)') ' '
       write(UN,'(a)') '######################################################'
       write(UN,'(a)') 'Grid info :'
       write(UN,'(a)') ' '
       write(UN,*)    'geometry  = ', thisGrid%geometry
       write(UN,*)    'cartesian = ', thisGrid%cartesian
       write(UN,*)    'polar     = ', thisGrid%polar
       write(UN,*)    'nx        = ', thisGrid%nx
       write(UN,*)    'ny        = ', thisGrid%ny
       write(UN,*)    'nz        = ', thisGrid%nz
       write(UN,*)    'na1       = ', thisGrid%na1
       write(UN,*)    'na2       = ', thisGrid%na2
       write(UN,*)    'na3       = ', thisGrid%na3
       write(UN,*)    'dipoleOffset  = ',thisGrid%dipoleOffset
       write(UN,'(a)') '#######################################################'
       write(UN,'(a)') ' '
       write(Un,'(a)') ' '
              
    end if
    
    
    if (filename(1:1) /= '*')  close(UN)

  end subroutine grid_info
    


  !
  ! It overlay the AMR grid values for a give number of planes sequentially.
  ! Uses the plot_AMR_values routine in this module.
  subroutine plot_AMR_planes(grid, name, plane,  nplane, filename, logscale, withgrid, &
       nmarker, xmarker, ymarker, zmarker, show_value_3rd_dim)
    implicit none
    type(gridtype), intent(in) :: grid
    character(len=*), intent(in)  :: name     ! "rho", "temperature", chiLine", "etaLine", 
    !                                         ! "etaLine", or "etaCont"
    character(len=*), intent(in)  :: plane    ! must be 'x-y', 'y-z' or 'z-x'
    integer, intent(in)           :: nplane   ! the number of planes to be plotted.
    character(len=*), intent(in)  :: filename ! the head of filename
    logical, intent(in)           :: logscale ! logscale if T, linear if not
    logical, intent(in)           :: withgrid ! plot

    ! To put the markers to indicate the important points in the plots.
    integer, intent(in) :: nmarker           ! number of markers
    real, intent(in)    :: xmarker(nmarker)  ! position of x
    real, intent(in)    :: ymarker(nmarker)  ! position of y
    real, intent(in)    :: zmarker(nmarker)  ! position of z
    logical, intent(in) :: show_value_3rd_dim! If T, the value of the third dimension will be shown.    

    !
    integer :: i
    real    :: val_3rd_dim
    real    :: d   ! the size of the largest cell.
    character(LEN=50) :: device
    real :: d0  ! off set for the cell center.


    d = grid%octreeRoot%subcellSize*2.0

    select case(plane)
    case("x-y")
       d0=grid%octreeRoot%centre%z
    case("y-z")
       d0=grid%octreeRoot%centre%x
    case("z-x", "x-z")
       d0=grid%octreeRoot%centre%y
    case default
       write(*,*) "Error:: Unknown name for plane passed to plot_AMR_planes."
       stop
    end select

    do i = 1, nplane
       val_3rd_dim = d0 -d/2.0 + real(i-1)*d/real(nplane-1)

       ! Setting up the name for the output file....
       ! using a function in utils_mod.f90
       device = tail_num_to_char(filename, i)
       device = TRIM(device)//".ps/vcps"
       
       call plot_AMR_values(grid, name, plane, val_3rd_dim, &
            device, logscale, withgrid, &
            nmarker, xmarker, ymarker, zmarker, d, show_value_3rd_dim)
    end do

  end subroutine plot_AMR_planes



  ! 
  ! Given a set of position arrays (xmarker, ymarker, zmarker), this 
  !routine plots the points only when the value of the third dimension is 
  ! within in the width of the third dimension specified in an input 
  ! parameter (width_3rd_dim) from value_3rd_dim.
  !
  ! Note:This routine assumes that a PGPLOT device is already open.
  !
  subroutine draw_markers(nmarker, xmarker, ymarker, zmarker, &
       plane, value_3rd_dim, width_3rd_dim, show_value_3rd_dim) 
    implicit none
    integer, intent(in) :: nmarker  ! number of markers
    real, intent(in) :: xmarker(nmarker)
    real, intent(in) :: ymarker(nmarker)
    real, intent(in) :: zmarker(nmarker)
    character(LEN=*), intent(in) :: plane
    real, intent(in) :: value_3rd_dim
    real, intent(in) :: width_3rd_dim
    ! If T, write the values of the thrid dimension in text (next to the point).
    logical, intent(in) :: show_value_3rd_dim  
    !
    character(LEN=30) :: text
    integer :: i
    

    !
    ! Draw the markers on the plot with the value of the third dimension.
    !   

    ! Draw points first..
    if (nmarker /=0) then  ! plots markers
       if (plane(1:3) == 'x-y') then           
          do i = 1, nmarker
             if ( ABS(zmarker(i)-value_3rd_dim) < width_3rd_dim/2.0 ) then
                call PGPT1(xmarker(i), ymarker(i), 21) 
                if (show_value_3rd_dim) then
                   write(text, '(1PE13.5)') zmarker(i)
                   call PGTEXT(xmarker(i), ymarker(i), TRIM(text))
                end if
             end if
          end do
       elseif (plane(1:3) == 'y-z') then
          do i = 1, nmarker
             if ( ABS(xmarker(i)-value_3rd_dim) < width_3rd_dim/2.0 ) then
                call PGPT1(ymarker(i), zmarker(i), 21) 
                if (show_value_3rd_dim) then
                   write(text, '(1PE13.5)') xmarker(i)
                   call PGTEXT(ymarker(i), zmarker(i), TRIM(text))
                end if
             end if
          end do
       elseif (plane(1:3) == 'z-x') then 
          do i = 1, nmarker
             if ( ABS(ymarker(i)-value_3rd_dim) < width_3rd_dim/2.0 ) then                
                call PGPT1(zmarker(i), xmarker(i), 21) 
                if (show_value_3rd_dim) then
                   write(text, '(1PE13.5)') ymarker(i)
                   call PGTEXT(zmarker(i), xmarker(i), TRIM(text))
                end if
             end if
          end do
       elseif (plane(1:3) == 'x-z') then 
          do i = 1, nmarker
             if ( ABS(ymarker(i)-value_3rd_dim) < width_3rd_dim/2.0 ) then
                call PGPT1(xmarker(i), zmarker(i), 21)
                if (show_value_3rd_dim) then
                   write(text, '(1PE13.5)') ymarker(i)
                   call PGTEXT(xmarker(i), zmarker(i), TRIM(text))
                end if
             end if
          end do
       end if
    else
       ! do nothing
       continue
    end if

  end subroutine draw_markers



  !
  ! Routine to shift the poistion of stars in a cluster to be at the center of leaf node
  ! assuming (hoping) that no two stars are in the same cell. This is done to avoid 
  ! or minimize the asymmetry in the temperature map computed later. 
  ! Note: this must be done only after computing the octree tree without the disc.
  !       If you want to include the disc, you must recompute the grid with new postion...
  ! 

  subroutine move_stars_to_cell_center(a_cluster, grid)
    implicit none
    type(cluster), intent(inout) :: a_cluster
    type(gridtype), intent(in) :: grid

    integer :: nstar
    integer :: i 
    type(sourcetype) :: a_star
    type(octalvector) :: position, newposition
    type(OCTAL), pointer :: thisOctal

    
    nstar = get_nstar(a_cluster)

    write (*, *) " "
    write (*, *) "Shifting the star positions...."
    write (*, *) " "
    write (*, *) "-- i --- old position and new position "     
    do i = 1, nstar       
       ! Finding the position of the stars in the list,
       a_star = get_a_star(a_cluster,i)         
       position = OCTALVECTOR(a_star%position%x, a_star%position%y, a_star%position%z)

       ! Finding the octal which contains the position.
       call amrGridValues(grid%octreeRoot, position, foundOctal=thisOctal)

       newposition = thisOctal%centre
       
       write(*,*) i, position 
       write(*,*) i, newposition 
       
       ! shifting the position to new position (the cell center).
       a_star%position = newposition 

       ! reinserting the star in the cluster.
       call put_a_star(a_cluster, a_star, i)


    end do

    write(*,*) " "
    write(*,*) "Finished shifting the stars ..."
    write(*,*) " "
    
              
  end subroutine move_stars_to_cell_center



  !  
  !  This subrouine recursively writes out the distance from the center of the
  !  given plane,  with a give value of the 3rd dimension, and the correcsponding 
  !  value stored in a octal to a file specified by the logical unit number. 
  !  
  recursive subroutine radial_profile(thisOctal, name, plane, val_3rd_dim, luout, center, grid)
    implicit none
    !
    type(octal), pointer   :: thisOctal
    type(gridtype) :: grid
    character(len=*), intent(in)  :: name     ! "rho", "temperature", chiLine", "etaLine",  
    !                                         ! "etaCont", "Vx", "Vy" or "Vz"
    character(len=*), intent(in)  :: plane    ! must be 'x-y', 'y-z', 'z-x', 'x-z'
    ! The value of the third dimension.
    ! For example, if plane = "x-y", the third dimension is the value of z.
    ! Then, when  val_3rd_dim = 0.0, this will plot the density (and the grid)
    ! on the z=0 plane, and so on...
    real, intent(in)              :: val_3rd_dim 
    integer, intent(in) :: luout  ! file unit number for the output
    type(octalvector), intent(in) :: center   ! position of the center of the root node
    !
    !
    type(octal), pointer  :: child 
    type(octalvector) :: rvec
    real :: value

    integer :: subcell, i, idx
    real :: t
    real :: xp, yp, xm, ym, zp, zm
    real(double) :: kabs, ksca
    integer :: ilam
    real(double) :: d, L, eps, distance
    logical :: use_this_subcell

  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call radial_profile(child, name, plane, val_3rd_dim, &
                     luout, center, grid)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          L = thisOctal%subcellSize
          d = L/2.0d0
          eps = d/1000.0d0
          xp = REAL(rVec%x + d)
          xm = REAL(rVec%x - d)
          yp = REAL(rVec%y + d)
          ym = REAL(rVec%y - d)
          zp = REAL(rVec%z + d)
          zm = REAL(rVec%z - d)     
          
          use_this_subcell = .false.
          if ( plane(1:3) == "x-y" .and. (ABS(rVec%z-val_3rd_dim)-d) < eps) then
             use_this_subcell = .true.      
          elseif ( plane(1:3) == "y-z" .and. (ABS(rVec%x-val_3rd_dim)-d) < eps) then 
             use_this_subcell = .true.
          elseif ( plane(1:3) == "z-x" .and. (ABS(rVec%y-val_3rd_dim)-d) < eps) then 
             use_this_subcell = .true.
          elseif ( plane(1:3) == "x-z" .and. (ABS(rVec%y-val_3rd_dim)-d) < eps) then 
             use_this_subcell = .true.
          else
             use_this_subcell = .false.
          end if
          
          if ((plane(1:3)=="x-y").and.(thisOctal%twoD)) then
             use_this_subcell = .false.
             if  ( ((rVec%z-L/2.d0) < 0.d0).and.((rvec%z+L/2.d0) >= 0.)) use_this_subcell = .true.
          endif

          
          if (use_this_subcell) then
             select case (name)
             case("rho")
                value = thisOctal%rho(subcell)
             case("temperature")
                value = thisOctal%temperature(subcell)
             case("chiLine")
                value = thisOctal%chiLine(subcell)
             case("etaLine")
                value = thisOctal%etaLine(subcell)
             case("etaCont")
                value = thisOctal%etaCont(subcell)
             case("Vx")
                value = thisOctal%velocity(subcell)%x * cSpeed/1.0d5 ![km/s]
             case("Vy")
                value = thisOctal%velocity(subcell)%y * cSpeed/1.0d5 ![km/s]
             case("Vz")
                value = thisOctal%velocity(subcell)%z * cSpeed/1.0d5 ![km/s]
!             case("tau")
!                call returnKappa(grid, thisOctal, subcell, ilam, grid%lamArray(ilam), kappaSca=ksca, kappaAbs=kabs)
!                value = thisOctal%subcellsize * (kSca+kAbs)
             case default
                value = 666.
             end select

!             distance = modulus(rvec-center)   ! length of the vector
             distance = modulus(rVec)
             write(luout, '(2(2x, 1PE18.4))')  distance, value


          end if
       end if

    end do

  end subroutine radial_profile

  subroutine dullemondplot(grid, rinner, router, angMax, nSource, source)
    use input_variables, only : rcore
    type(GRIDTYPE) :: grid
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    real :: rinner, router, angMax
    real(double) :: r, ang
    real :: tr(6)
    integer :: nc
    real :: tc(100), rc(100)
    integer :: nx, ny
    real, allocatable :: rhoimage(:,:),tempimage(:,:)
    real :: dx, dy, fg, bg
    integer :: i, j, k
    integer :: ilo, ihi
    real(double) :: rtemp
    real :: ttemp
    real(double), allocatable  :: kabs(:), ksca(:), kros(:)   ! (size=maxTau)
    real, allocatable :: lambda(:), dlambda(:)
    real, allocatable :: tauSca(:), tauAbs(:), tauExt(:)
    integer :: ilambda
    integer :: maxTau = 10000
    integer :: nr = 1000
    integer :: error, nTau
    logical :: hitcore
    real :: rtau(1000)
    real, allocatable :: dv(:), temp(:)
    real,allocatable :: chiline(:)
    real(double), allocatable :: rho(:)
    real(double), allocatable :: Levelpop(:,:),ne(:)
    logical, allocatable :: inflow(:)
    type(VECTOR),allocatable :: vel(:)
    logical :: first
    real(double),allocatable :: dummy(:)

    type(OCTALVECTOR) :: octVec, Uhatoctal, avecoctal
    integer :: pgbegin

    allocate(lambda(1:maxTau),dlambda(1:maxTau),ksca(1:maxtau), &
         kabs(1:maxtau), kros(1:maxtau))
    allocate(tauAbs(1:maxTau), tauSca(1:maxTau), tauExt(1:maxTau))
    allocate(vel(1:maxTau), dv(1:maxTau),temp(1:maxtau), rho(1:maxtau))
    allocate(chiline(1:maxTau),ne(1:maxTau),levelPop(1:maxtau,grid%maxLevels))
    allocate(inflow(1:maxTau))
    allocate(dummy(1:maxTau))
    nx = 1000
    ny = 1000

    nc = 20
    do i = 1, nc
       tc(i) = real(i)*100.
    enddo
    do i = 1, nc
       rc(i) = -8.-0.5*real(i-1)
    enddo



    dx = log10(rOuter/rInner)/real(nx)
    dy = angMax/real(ny)
    tr(1) =  - dx + log10(rInner)
    tr(2) = dx
    tr(3) = 0.
    tr(4) =  - dy + 0.
    tr(5) = 0.
    tr(6) = dy


    allocate(rhoImage(1:nx,1:ny),tempImage(1:nx,1:ny))
    
    do i = 1, nx
       do j = 1, ny
          r = log10(rinner)+real(i-1)*log10(rOuter/rInner)/real(nx-1)
          r = 10.d0**r
          ang = + pi /2. - angMax * real(j-1)/real(ny-1) 
          octVec = OCTALVECTOR(r*sin(ang),0.d0,r*cos(ang))
          call amrGridValues(grid%octreeRoot,octVec,rho=rtemp, temperature=ttemp)
          rhoimage(i,j) = log10(max(1.d-30,rtemp))
          tempimage(i,j) = ttemp
       enddo
    enddo

    do i = 1, nr
       ang = real(i-1)/real(nr-1) * pi/2.
       uHatOctal = OCTALVECTOR(sin(ang), 0.d0, cos(ang))
       aVecOctal = dble(rCore) * uhatOctal
       ilambda = 32
    
       CALL startReturnSamples2 (aVecOctal,uHatOctal,grid,1.,nTau,       &
            maxTau,.false.,.true.,hitcore,.false.,iLambda,error,&
            lambda,nSource,source,kappaAbs=kAbs,kappaSca=kSca,velocity=vel,velocityderiv=dv, &
            temperature=temp, &
            chiLine=chiLine,    &
            levelPop=levelPop,rho=rho, &
            Ne=Ne, inflow=inflow, etaLine=dummy, etaCont=dummy)

       dlambda(1:nTau-1) = lambda(2:nTau) - lambda(1:nTau-1)  


       tauAbs(1:2) = 0.
       tauSca(1:2) = 0.

       do j = 2, nTau, 1
          tauSca(j) = tauSca(j-1) + dlambda(j-1)*0.5*(ksca(j-1)+ksca(j))
          tauAbs(j) = tauAbs(j-1) + dlambda(j-1)*0.5*(kabs(j-1)+kabs(j))
       enddo
       tauExt(1:nTau) = tauSca(1:nTau)+tauAbs(1:nTau)
       if (tauExt(nTau) > 1.) then
          call locate(tauExt, nTau, 1., k)
          rtau(i) = lambda(k) + (lambda(k+1)-lambda(k))*(1.-tauExt(k))/(tauExt(k+1)-tauExt(k))
       else
          rtau(i) = 0.
       endif
    enddo

    
    fg = -8.
    bg = -18.
    i = pgbegin(0,"dull.ps/cps",1,1)
    call pgvport(0.1,0.9,0.1,0.9)
    call pgqcir(ilo, ihi)
    call pgscir(ilo, ihi)
    call palette(1)

    call pgwindow(log10(rinner),log10(router), 0., angMax)
    call pgimag(rhoimage, nx, ny, 1, nx, 1, ny, fg, bg, tr)
    call pgsci(2)
    call pgcont(rhoimage, nx, ny, 1, nx, 1, ny, rc, nc, tr)
    call pgsci(3)
    call pgcont(tempimage, nx, ny, 1, nx, 1, ny, tc, nc, tr)
    call pgsci(4)
    first = .true.
    do i = 2, nr
       ang = pi/2. - real(i-1)/real(nr-1) * pi/2.
       if (rTau(i) /= 0.) then
          if (first) then
             call pgmove(real(log10(rtau(i))), real(ang))
             first = .false.
          else
             call pgdraw(real(log10(rtau(i))), real(ang))
          endif
       endif
    enddo
    call pgsci(1)
    call pgwindow(log10(rinner/1.5e3),log10(router/1.5e3), 0., angMax)
    call pgbox('bclnst',0.0,0,'bcnst',0.0,0)
    call pglab("R [AU]","\gp/2 -\gt", " ")
    call pgend
    
    deallocate(rhoimage,tempimage, ksca, kabs, lambda, dlambda)
    deallocate(dummy)
  end subroutine dullemondplot

  subroutine nattaplot(grid, nSource, source)
    use input_variables, only : rcore, rinner, router
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    type(GRIDTYPE) :: grid
    real :: angMax
    real(double) :: r, ang
    real :: tr(6)
    integer :: nc
    real :: tc(100), rc(100)
    integer :: nx, ny
    real, allocatable :: rhoimage(:,:),tempimage(:,:)
    real :: dx, dy, fg, bg
    integer :: i, j, k
    integer :: ilo, ihi
    real(double) :: rtemp
    real :: ttemp
    real(double), allocatable  :: kros(:),ksca(:),kabs(:)   ! (size=maxTau)
    real, allocatable :: lambda(:), dlambda(:)
    real, allocatable :: tauSca(:), tauAbs(:), tauExt(:)
    integer :: ilambda
    integer :: maxTau = 10000
    integer :: nr = 1000
    integer :: error, nTau
    logical :: hitcore
    real :: rtau(1000)
    real, allocatable :: dv(:), temp(:)
    real,allocatable :: chiline(:)
    real(double), allocatable :: rho(:)
    real(double), allocatable :: Levelpop(:,:),ne(:)
    logical, allocatable :: inflow(:)
    type(VECTOR),allocatable :: vel(:)
    logical :: first
    real :: x, z
    type(OCTALVECTOR) :: octVec, Uhatoctal, avecoctal
    integer :: pgbegin
    real(double), allocatable :: dummy(:)

    allocate(lambda(1:maxTau),dlambda(1:maxTau),kros(1:maxtau))
    allocate(kabs(1:maxTau),ksca(1:maxTau))
    allocate(tauAbs(1:maxTau), tauSca(1:maxTau), tauExt(1:maxTau))
    allocate(vel(1:maxTau), dv(1:maxTau),temp(1:maxtau), rho(1:maxtau))
    allocate(chiline(1:maxTau),ne(1:maxTau),levelPop(1:maxtau,grid%maxLevels))
    allocate(inflow(1:maxTau))
    allocate(dummy(maxtau))
    nx = 1000
    ny = 1000

    nc = 20
    do i = 1, nc
       tc(i) = real(i)*100.
    enddo
    do i = 1, nc
       rc(i) = -8.-0.5*real(i-1)
    enddo

    angmax = pi/2.

    dx = log10(rOuter/rInner)/real(nx)
    dy = angMax/real(ny)
    tr(1) =  - dx + log10(rInner)
    tr(2) = dx
    tr(3) = 0.
    tr(4) =  - dy + 0.
    tr(5) = 0.
    tr(6) = dy

    open(21, file="natta.dat", status="unknown", form="formatted")

    do i = 1, nr
       ang = real(i-1)/real(nr-1) * pi/2.
       uHatOctal = OCTALVECTOR(sin(ang), 0.d0, cos(ang))
       aVecOctal = dble(rCore) * uhatOctal
       ilambda = 32
    
       CALL startReturnSamples2 (aVecOctal,uHatOctal,grid,1.,nTau,       &
            maxTau,.false.,.true.,hitcore,.false.,iLambda,error,&
            lambda, nSource, source, kappaAbs=kabs, kappaSca=ksca, &
            velocity=vel,velocityderiv=dv, &
            temperature=temp, &
            chiLine=chiLine,    &
            levelPop=levelPop,rho=rho, &
            Ne=Ne, inflow=inflow, etaCont=dummy, etaLine=dummy)

       dlambda(1:nTau-1) = lambda(2:nTau) - lambda(1:nTau-1)  


       tauExt(1:2) = 0.

       do j = 2, nTau, 1
          tauExt(j) = tauExt(j-1) + dlambda(j-1)*0.5*(kabs(j-1)+kabs(j))
       enddo
       if (tauExt(nTau) > 1.) then
          call locate(tauExt, nTau, 1., k)
          rtau(i) = lambda(k) + (lambda(k+1)-lambda(k))*(1.-tauExt(k))/(tauExt(k+1)-tauExt(k))
       else
          rtau(i) = 0.
       endif
    enddo
    
    i = pgbegin(0,"natta.ps/cps",1,1)
    call pgvport(0.1,0.9,0.1,0.9)

    call pgwnad(-0.05*real(autocm)/1.e10+rinner, rinner+0.2*real(autocm)/1.e10, 0., 0.25*real(autocm)/1.e10)
    first = .true.
    do i = 2, nr
       ang = pi/2. - real(i-1)/real(nr-1) * pi/2.
       x = rtau(i) * cos(ang)
       z = rtau(i) * sin(ang)
       write(21,*) 1.e10*x/autocm, 1.e10*z/autocm
       if (rTau(i) /= 0.) then
          if (first) then
             call pgmove(x, z)
             first = .false.
          else
             call pgdraw(x,z)
          endif
       endif
    enddo
    call pgsci(1)
    call pgwnad(0., 0.2, 0., 0.2)
    call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
    call pglab("R [AU]","Z [AU]", " ")
    call pgend
    close(21)
    
    deallocate(kros, kabs, ksca, lambda, dlambda)
    deallocate(dummy)
  end subroutine nattaplot

    
end module grid_mod


