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

  use kind_mod
  use constants_mod                   ! physical constants
  use messages_mod
  use vector_mod                      ! vector math
  use amr_mod, only: deleteOctreeBranch, deleteOctal, deallocateOctalDynamicAttributes, &
       octalOnThread
  use octal_mod, only: OCTAL
  use gridtype_mod, only: GRIDTYPE    ! type definition for the 3-d grid
  use mpi_global_mod, only: myRankGlobal, nThreadsGlobal
  use utils_mod, only: locate

  implicit none

  public

#ifdef HYDRO
  private :: writeReal1D, writeReal2D, writeDouble2D
  private :: readReal1D,  readReal2D,  readDouble2D
#endif

  interface getIndices
     module procedure getIndices_single
     module procedure getIndices_octal
  end interface
     

contains



  subroutine initAMRgrid(grid)

    use inputs_mod
    use cmfgen_class, only: get_cmfgen_data_array_element
    use jets_mod, only: initJetsAMR
    use amr_mod, only: initTTauriAMR, initWindTestAMR
    use gridtype_mod, only: statEqMaxLevels

    implicit none

    ! grid%timeNow must be assigned before this routine is called!

    type(GRIDTYPE), intent(inout) :: grid                 ! the grid
    real :: rStar
    
    grid%oneKappa = oneKappa
    
    if (oneKappa) then
       allocate(grid%oneKappaAbs(nDustType,1:nLambda))
       allocate(grid%oneKappaSca(nDustType,1:nLambda))
       grid%oneKappaAbs = 0.d0; grid%oneKappaSca = 0.d0
       grid%nTempRossArray = 1000
       allocate(grid%kappaRossArray(nDustType,1:grid%nTempRossArray))
       allocate(grid%tempRossArray(1:grid%nTempRossArray))
    endif

    grid%lineEmission = lineEmission
    grid%maxLevels = statEqMaxLevels

    if ( dustPhysics ) then 
       grid%flatspec = .false.
    else
       grid%flatspec = .true.
    end if

    grid%statEq2d = statEq2d

    grid%amr2dOnly = amr2dOnly
    grid%geometry = geometry

    grid%idump = 0
    grid%currentTime = 0.d0


    select case (geometry)
    
    case("toruslogo")
       grid%geometry = "toruslogo"
       grid%rStar1 = real(rSol/1.e10)

    case("whitney")
       grid%geometry = "whitney"
       grid%rCore = rStellar/1.e10
       grid%lCore = fourPi * rCore**2 * stefanBoltz * teff**4 * 1.e20
       grid%rInner = real(erInner,si)/1.e10
       grid%rOuter = real(erOuter,si)/1.e10
       grid%starPos1 = vector(1.e-6,1.e-6,1.e-6)

    case("benchmark")
       grid%geometry = "benchmark"
       grid%rCore = rCore
       grid%lCore = fourPi * rCore**2 * stefanBoltz * sourceTeff(1)**4 * 1.e20
       grid%rInner = rInner
       grid%rOuter = rOuter

    case("circumbin")
       grid%geometry = "circumbin"
       grid%rInner = rInner
       grid%rOuter = rOuter
       grid%rCore = real(rSol/1.e10)

    case("shakara")
       grid%geometry = "shakara"
       grid%rCore = rCore
       grid%lCore = fourPi * rCore**2 * stefanBoltz * teff**4 * 1.e20
       grid%rInner = rInner
       grid%rOuter = rOuter

   case("iras04158")
       grid%geometry = "iras04158"
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

    case("lexington","fractal","runaway")
       grid%rCore = real(18.67 * rSol / 1.e10)
       grid%rInner = rinner
       grid%rOuter = 2.e09
       grid%oneKappa = .true.
       oneKappa = .true.

!    case ("ttauri","magstream") 
!       call initTTauriAMR(grid)
    case ("windtest") 
       call initWindTestAMR(grid)

    case ("jets") 
       call initJetsAMR(grid)

    case ("luc_cir3d") 
       rStar  = real(CIR_Rstar*Rsol/1.0d10)   ! in [10^10cm] 
       grid%rCore = rStar
       grid%rStar1 = rStar
       grid%starPos1 = vector(0.,0.,0.)

    case ("cmfgen") 
       rStar  = real(get_cmfgen_data_array_element("R", 1))   ! in [10^10cm] 
       grid%rCore = rStar
       grid%rStar1 = rStar
       grid%starPos1 = vector(0.,0.,0.)

    case ("romanova") 
       rStar  = real(ROM_Rs*Rsol/1.0d10)   ! in [10^10cm] 
       grid%rCore = rStar
       grid%rStar1 = rStar
       grid%rStar2 = 0.
       grid%starPos1 = vector(0.,0.,0.)
       grid%diskRadius = ThinDiskRin  ! [10^10cm]
       grid%diskNormal = VECTOR(0.,0.,1.)

    case ("testamr")
       call initTestAMR(grid)

    case("starburst")
       grid%rCore = real(18.67 * rSol / 1.e10)
       grid%rInner = rinner
       grid%rOuter = 2.e09
       grid%oneKappa = .true.
       oneKappa = .true.
       grid%geometry = "starburst"
       grid%lineEmission = .false.

    case ("proto")
       call initProtoAMR(grid)

    case("sphfile","cluster","wr104","molcluster", "theGalaxy")
       grid%lineEmission = .false.
       grid%rCore        = rCore
       
    case("spiralwind","wind")
       grid%rCore = rCore 

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
       ! In the defalt case nothing special needs to be done
       grid%geometry = geometry

    end select
    


    nullify(grid%rho)
    nullify(grid%kappaAbs)
    nullify(grid%kappaSca)
    nullify(grid%velocity)
    nullify(grid%temperature)
    nullify(grid%chiLine)
    nullify(grid%etaLine)
    nullify(grid%etaCont)
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
       grid%phiAxis(i) = real(twoPi*real(i-1)/real(grid%nphi-1))
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
       theta = real(pi*real(i-1)/real(nTheta-1))
       do j = 1, nPhi
          phi = real((2.*real(j-1)/real(nPhi-1)-1.)*pi)
          do k = 1, nRad
             r = radius * real(k-1)/real(nRad-1)
             rVec%x = r*sin(theta)*cos(phi)
             rVec%y = r*sin(theta)*sin(phi)
             rVec%z = r*cos(theta)
             w = cos(theta)
             call locate(thisGrid%xAxis,thisGrid%nx,real(rVec%x), i1)
             call locate(thisGrid%yAxis,thisGrid%ny,real(rVec%y), i2)
             call locate(thisGrid%zAxis,thisGrid%nz,real(rVec%z), i3)

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

    grid%rCore = real(rSol / 1.e10)
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

          phi = real(2.*pi*real(j-1)/real(nPhi-1))
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
    mdot = real(10.**logMassLossRate * mSol / (365.25 * 24. * 60. * 60.))

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
             rho = real(mdot / (4. * pi * r**2 * v))
             if (.not.dust) then
                grid%rho(i,j,k) = real(rho/mHydrogen)
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
       rotAng = real(twoPi * real(j-1)/199.)
       do i = 1, 100
          do k = 1, 100
             ang = minAng + (maxAng-minAng)*(real(k-1)/99.)
             rVec%x = grid%xAxis(grid%nx) * real(i-1)/99.
             rVec%y = rVec%x * tan(ang)
             rVec%z = 0.
             rVec = rotateX(rVec, dble(rotAng))
             rVec%x = rVec%x + binarySep  - stagPoint
             call locate(grid%xAxis,Grid%nx, real(rVec%x), i1)
             call locate(grid%yAxis,Grid%ny, real(rVec%y), i2)
             call locate(grid%zAxis,Grid%nz, real(rVec%z), i3)
             r = real(modulus(rVec))
             rho = real(mdot / (4. * pi * r**2 * v))
             if (.not.done(i1,i2,i3)) then
                if (dust) then

                   if ((r > 2.*binarySep).and.(r < 100.*binarySep)) then
                      rho = real(rho / (4.*mHydrogen))  ! helium number density)
                      rho = rho * 0.2             ! carbon number density
                      rho = real(rho * (12.*mHydrogen)) ! carbon mass density
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
       write(*,'(a,1pe9.2)') "Total dust mass (solar masses): ",totalDustMass/mSol
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
       rInner = real(abs(sqrt(abs(grid%zAxis(i)))*sin(openingAng*degToRad)))
       rOuter = rInner * 0.1
       do j=1,nCircle
          ang = real(twoPi * real(j-1)/real(nCircle-1))
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
          ang = real(twoPi * real(i-1)/real(nCircle-1))
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
       grid%phiAxis(i) = real(twoPi*real(i-1)/real(grid%nphi-1))
    enddo


    do i = 1, grid%nr
       do j = 1, grid%nmu
          do k = 1 , grid%nphi
             vr = 10.e5+(vel-10.e5)*(1.-grid%rAxis(1)/grid%rAxis(i))
             grid%rho(i,j,k) = real(log(mdot)-log(fourPi)-log(vr)- &
                  2.*log(grid%rAxis(i)*scale)-log(mHydrogen)+log(1.-kfac*grid%muAxis(j)**2))
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
    real :: radius, mdot, rMin, rMax, vel, v,  x
    real(double) :: r
    integer :: i, j, k
    integer :: nSpiral, iSpiral
    integer :: i1,i2
    real :: muStart,muEnd
    real(double) :: thickness,  woundFac
    integer :: j1,j2
    real(double) :: phi, theta, phaseOffset, spiralPhase
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
       grid%phiAxis(i) = real(twoPi*real(i-1)/real(grid%nphi-1))
    enddo




    do i = 1, grid%nr
       do j = 1, grid%nmu
          do k = 1 , grid%nphi
             v = 10.e5+(vel-10.e5)*(1.- radius/grid%rAxis(i))**1
             grid%rho(i,j,k) = real(log(mdot)-log(fourPi)-log(v)-2.*log(grid%rAxis(i)*scale)-log(mHydrogen))
             grid%rho(i,j,k) = exp(grid%rho(i,j,k))
          enddo
       enddo
    end do

    muStart = 0.5
    muEnd = -0.5
    thickness = 0.2
    woundFac = 0.75
    phaseOffset = 0.25*twoPi
    nSpiral = 2

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

          startVec = (r*rHat) - ((thickness/2.d0)*r)*perp
          endVec = (r*rHat) + ((thickness/2.d0)*r)*perp


          do i1 = 1,100
             x = real(i1-1)/99.
             thisVec = startVec + (dble(x)*(endVec-startVec))
             call getPolar(thisVec, r, theta, phi)
             call locate(grid%rAxis,grid%nr,real(r),i2)
             call locate(grid%phiAxis,grid%nPhi,real(phi),k)

             do j = min(j1,j2), max(j1,j2)
                v = 10.e5+(vel-10.e5)*(1.- radius/grid%rAxis(i2))**1
                grid%rho(i2,j,k) = real(log(mdot)-log(fourPi)-log(v)-2.*log(grid%rAxis(i2)*scale)-log(mHydrogen) + log(5.))
                grid%rho(i2,j,k) = exp(grid%rho(i2,j,k))
             enddo
          enddo
       enddo


    enddo
  end subroutine fillGridSpiral




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
    sigmaScaBlue = real((34.+6.6)*sigmaE * 1.e10)
    sigmaAbsBlue = real(1.e-5 * sigmaE * 1.e10)

    sigmaScaRed = real(1.e-5 * sigmaE * 1.e10)
    sigmaAbsRed = real(1.e-5 * sigmaE * 1.e10)

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
             r = real(modulus(rVec))
             mu = real((grid%zAxis(k)-coolStarPosition%z) / r)
             fac = 1.
             !             fac = abs(mu)+1.
             if (r > rCool) then
                v =  max(1.e5,vinf*(1.-rCool/r)**beta)
                grid%velocity(i,j,k) = (v/cSpeed) * rHat
                !                grid%velocity(i,j,k)%z = grid%velocity(i,j,k)%z * fac**2
                rho = real(mdot / (4.* pi * (r*1.e10)**2 * v))
                rho = real(rho / mHydrogen                )
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
    use utils_mod, only: logInterp

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
       grid%phiAxis(i) = real(twoPi*real(i-1)/real(grid%nphi-1))
    enddo

    grid%isotropic = .true.
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
                grid%kappaSca(i,j,k,1) = logInterp(escat, nr, r, grid%rAxis(i))*fac
                tmp = logInterp(chi, nr, r,grid%rAxis(i))-grid%kappaSca(i,j,k,1)
                grid%kappaAbs(i,j,k,1) = max(1.e-20,tmp)
             else
                grid%velocity(i,j,k) = (1.e5*v(i)/cSpeed) * rHat
                grid%temperature(i,j,k) = t(i)*1.e4
                grid%chiLine(i,j,k) = chil(i)*fac**2
                grid%etaLine(i,j,k) = etal(i)*fac**2
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




  subroutine freeGrid(grid)

    type(GRIDTYPE) :: grid

    call writeInfo("Deallocating grid structure", TRIVIAL)
    
    if (associated(grid%octreeRoot)) then
      call deleteOctal(grid%octreeRoot, deleteChildren=.true.,adjustparent=.false.)
      call deallocateOctalDynamicAttributes(grid%octreeRoot)
      deallocate(grid%octreeRoot)
      nullify(grid%octreeRoot)
    end if
      
    if (associated(grid%oneKappaAbs)) deallocate(grid%oneKappaAbs)
       nullify(grid%oneKappaAbs)
    if (associated(grid%oneKappaSca)) deallocate(grid%oneKappaSca)
       nullify(grid%oneKappaSca)
    if (associated(grid%kappaRossArray)) deallocate(grid%kappaRossArray)
       nullify(grid%kappaRossArray)
    if (associated(grid%tempRossArray)) deallocate(grid%tempRossArray)
       nullify(grid%tempRossArray)
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

    call writeinfo("done.",TRIVIAL)
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



  logical pure function outsideGrid(posVec, grid)    
    type(GRIDTYPE), intent(in) :: grid
    type(VECTOR), intent(in)   :: posVec

    outSideGrid = .false.

    if ((posVec%x > grid%xAxis(grid%nx)) .or. &
         (posVec%y > grid%yAxis(grid%ny)) .or. &
         (posVec%z > grid%zAxis(grid%nz)) .or. &
         (posVec%x < grid%xAxis(1)) .or. &
         (posVec%y < grid%yAxis(1)) .or. &
         (posVec%z < grid%zAxis(1))) outsideGrid = .true.

  end function outsideGrid


  subroutine getIndices_single(grid, rVec, i1, i2, i3, t1, t2, t3)
    use utils_mod, only: hunt
    ! making this PURE may cause problems with XL Fortran
    type(GRIDTYPE), intent(in) :: grid
    type(VECTOR), intent(in)   :: rVec
    integer, intent(out)       :: i1, i2, i3
    real, intent(out)          :: t1, t2 ,t3
    real(double)                       :: r, theta, phi, mu

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

       t1 = real((rVec%x-grid%xAxis(i1))/(grid%xAxis(i1+1)-grid%xAxis(i1)))
       t2 = real((rVec%y-grid%yAxis(i2))/(grid%yAxis(i2+1)-grid%yAxis(i2)))
       t3 = real((rVec%z-grid%zAxis(i3))/(grid%zAxis(i3+1)-grid%zAxis(i3)))

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

       t1 = real((r-grid%rAxis(i1))/(grid%rAxis(i1+1)-grid%rAxis(i1)))


       call hunt(grid%muAxis, grid%nMu, real(mu), i2)
       if (i2 == grid%nMu) i2 = grid%nMu-1
       t2 = real((mu-grid%muAxis(i2))/(grid%muAxis(i2+1)-grid%muAxis(i2)))

       call hunt(grid%phiAxis, grid%nPhi, real(phi), i3)
       if (i3 == grid%nPhi) i3=i3-1
       t3 = real( (phi-grid%phiAxis(i3))/(grid%phiAxis(i3+1)-grid%phiAxis(i3)))

    endif

  end subroutine getIndices_single



  subroutine getIndices_octal(grid, rVec, i1, i2, i3, t1, t2, t3)
    use utils_mod, only: hunt
    ! making this PURE may cause problems with XL Fortran
    type(GRIDTYPE), intent(in)          :: grid
    type(VECTOR), intent(in)       :: rVec
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

#ifdef HYDRO
  subroutine writeAMRgridMPI(filename,fileFormatted,grid,mpiFlag)
    use inputs_mod, only : molecular, cmf, hydrodynamics
    ! writes out the 'grid' for an adaptive mesh geometry  

    implicit none
  
    character(len=*)           :: filename
    logical, intent(in)        :: fileFormatted
    type(GRIDTYPE), intent(in) :: grid
    
    integer, dimension(8) :: timeValues ! system date and time
    integer               :: error      ! error code
    logical, optional :: mpiFlag
    logical :: doMpiFlag
    logical :: append
    character(len=20) :: fileStatus, positionStatus

    doMpiFlag = .false.
    if (PRESENT(mpiFlag)) doMpiFlag = mpiFlag

    append = .false.
    fileStatus = "replace"
    positionStatus="rewind"

    if (dompiFlag .and. (myRankGlobal /= 1)) then
       append = .true.
       fileStatus = "old"
       positionStatus = "append"
    endif


    if (fileFormatted) then 
       open(unit=20,iostat=error, file=filename, form="formatted", status=fileStatus, position=positionStatus)
    else 
       open(unit=20,iostat=error, file=filename, form="unformatted", status=fileStatus, position=positionStatus)
    end if        
    call writeInfo("Writing AMR grid file to: "//trim(filename),TRIVIAL)
    
    call date_and_time(values=timeValues)
    

    if (.not.append) then

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
               grid%rOuter, grid%amr2dOnly, grid%photoionization, grid%iDump, grid%currentTime
          if (error /=0) then
             print *, 'Panic: write error in writeAMRgrid (formatted variables)' ; stop
          end if
          
          !       if (cmf) then
          !          write(unit=20,fmt=*) grid%nFreqArray, grid%freqArray(1:2000)
          !       endif
          
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
               grid%rOuter, grid%amr2dOnly, grid%photoionization, grid%iDump, grid%currentTime
          if (error /=0) then
             print *, 'Panic: write error in writeAMRgrid (unformatted variables)' ; stop
          end if
          
          !       if (cmf) then
          !          write(unit=20) grid%nFreqArray, grid%freqArray(1:2000)
          !       endif
          
          
       end if
       
       
       call writeReal1D(grid%lamarray,fileFormatted)
       call writeReal2D(grid%oneKappaAbs,fileFormatted)
       call writeReal2D(grid%oneKappaSca,fileFormatted)
       call writeClumps(fileFormatted)

    endif

!    if (dompiFlag.and.(myRankGlobal==0)) then
!       goto 666
!    endif
       
    ! now we call the recursive subroutine to store the tree structure 
    if (associated(grid%octreeRoot)) then
       if (dompiFlag.and.(myRankGlobal ==1)) then
          if (fileFormatted) then 
             write(unit=20,fmt=*) .true.
          else
             write(unit=20) .true.
          end if
       endif
       call writeOctreePrivate(grid%octreeRoot,fileFormatted, grid, doMpiFlag)
    else 
       if (fileFormatted) then 
          write(unit=20,fmt=*) .false.
       else
          write(unit=20) .false.
       end if
    end if
    
!666 continue
    endfile 20
    close(unit=20)
    
  contains
  
    recursive subroutine writeOctreePrivate(thisOctal,fileFormatted, grid, doMpiFlag)
       ! writes out an octal from the grid octree

       type(octal),  pointer :: thisOctal
       logical, intent(in)             :: fileFormatted
       type(gridtype) :: grid
       type(octal), pointer :: thisChild
       integer              :: iChild
       logical :: doMpiFlag
       logical :: writeOctal
       integer :: tempNChildren
       integer :: tempIndexChild(8), i
       logical :: tempHasChild(8)

       writeOctal = .true.



!       if (doMpiFlag.and.(thisOctal%nDepth >=2 ).and.octalOnThread(thisOctal, 1, myRankGlobal)) then
!          writeOctal=.true.
!       else
!          writeOctal=.false.
!       endif
       

       tempNChildren = thisOctal%nChildren
       tempIndexChild = thisOctal%indexChild
       tempHasChild = thisOctal%hasChild

       if (nThreadsGlobal == 9) then
          if (doMpiFlag.and.(thisOctal%nDepth < 2) .and. (myRankGlobal==1)) then
             writeOctal = .true.
             tempNChildren = 8
             do i = 1, 8
                tempIndexChild(i) = i
             enddo
             tempHasChild = .true.
          endif
          if (doMpiFlag.and.(thisOctal%nDepth==1).and.(myRankGlobal /=1)) writeOctal = .false.
       endif

       if (nThreadsGlobal == 65) then
          if (doMpiFlag.and.(thisOctal%nDepth < 3) .and. (myRankGlobal==1)) then
             writeOctal = .true.
             tempNChildren = 8
             do i = 1, 8
                tempIndexChild(i) = i
             enddo
             tempHasChild = .true.
          endif
          if (doMpiFlag.and.(thisOctal%nDepth < 3).and.(myRankGlobal /=1)) writeOctal = .false.

          if (doMpiFlag.and.(thisOctal%nDepth>=3).and. .not.(octalOnThread(thisOctal,1, myRankGlobal))) writeOctal = .false.

       endif
       
!       if (writeOctal) then
!          write(*,*) myRankGlobal, " writing depth ", thisOctal%nDepth, " writeoctal ", writeoctal, &
!               octalOnThread(thisOctal, 1, myRankGlobal), " label ",thisOctal%label, " mpithread ", &
!               thisOctal%mpithread
!       endif
       if (writeOctal) then

          if (fileFormatted) then 
             write(iostat=error,fmt=*,unit=20) thisOctal%nDepth, tempnChildren,&
                  tempindexChild, temphasChild, thisOctal%centre,& 
                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
                  thisOctal%subcellSize,thisOctal%threed, thisOctal%twoD,    &
                  thisOctal%maxChildren, thisOctal%dustType,  &
                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
          else
             write(iostat=error,unit=20) thisOctal%nDepth, tempnChildren, &
                  tempindexChild, temphasChild, thisOctal%centre,& 
                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
                  thisOctal%probDistCont, thisOctal%Ne, thisOctal%nh, thisOctal%nTot,      &
                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
                  thisOctal%subcellSize, thisOctal%threeD, thisOctal%twoD,   &
                  thisOctal%maxChildren, thisOctal%dustType, &
                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
          end if
          if (.not.grid%oneKappa) then
             !          call writeReal2D(thisOctal%kappaAbs,fileFormatted)
             !          call writeReal2D(thisOctal%kappaSca,fileFormatted)
             call writeDouble2D(thisOctal%kappaAbs,fileFormatted)
             call writeDouble2D(thisOctal%kappaSca,fileFormatted)
          endif
          if (grid%photoionization) then
             call writeDouble2D(thisOctal%ionFrac,fileFormatted)
             call writeDouble2D(thisOctal%photoIonCoeff,fileFormatted)
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
             call writeDouble2D(thisOctal%atomAbundance, fileFormatted)
             call writeDouble3D(thisOctal%atomLevel,fileFormatted)
             call writeDouble2D(thisOctal%jnuCont,fileFormatted)
             call writeDouble2D(thisOctal%jnuLine,fileFormatted)
             if (fileformatted) then
                write(unit=20,iostat=error,fmt=*) thisOctal%microturb
             else
                write(unit=20,iostat=error) thisOctal%microturb
             endif
          endif

          if (hydrodynamics) then
             if (fileformatted) then
                write(unit=20,iostat=error,fmt=*) thisOctal%boundaryCondition
             else
                write(unit=20,iostat=error) thisOctal%boundaryCondition
             endif
          endif
          
          if (grid%splitOverMpi) then
             if (fileformatted) then
                write(unit=20,iostat=error,fmt=*) thisOctal%mpiThread
             else
                write(unit=20,iostat=error) thisOctal%mpiThread
             endif
          endif
          
       endif

       if (thisOctal%nChildren > 0) then 
          do iChild = 1, thisOctal%nChildren, 1
             thisChild => thisOctal%child(iChild)
             call writeOctreePrivate(thisChild,fileFormatted,grid,dompiflag)
          end do
       end if
     end subroutine writeOctreePrivate
    
    
  end subroutine writeAMRgridMPI


  subroutine readAMRgridMPI(filename,fileFormatted,grid, mpiflag)
    ! reads in a previously saved 'grid' for an adaptive mesh geometry  

    use unix_mod, only: unixGetenv
    use inputs_mod, only: geometry,dipoleOffset,amr2dOnly,statEq2d, molecular, cmf, hydrodynamics
    implicit none

    character(len=*)            :: filename
    logical, intent(in)         :: fileFormatted
    type(GRIDTYPE), intent(inout) :: grid
    
    integer, dimension(8) :: timeValues    ! system date and time
    integer               :: error         ! status code
    logical               :: octreePresent ! true if grid has an octree
    integer :: nOctal
    character(len=80) :: absolutePath, inFile
    integer :: iJunk
    logical, optional :: mpiFlag
    logical :: doMpiFlag
    real, pointer :: junk(:) => null(), junk2(:,:) => null(), junk3(:,:)=>null()
    absolutePath = " "
    error = 0

    dompiFlag = .false.
    if (present(mpiFlag)) doMpiFlag = mpiFlag


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
               grid%rOuter, grid%amr2dOnly, grid%photoionization, grid%iDump, grid%currentTime
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
               grid%rOuter, grid%amr2dOnly, grid%photoionization, grid%iDump, grid%currentTime
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
    if (associated(junk)) deallocate(junk)
    if (associated(junk2)) deallocate(junk2)
    if (associated(junk3)) deallocate(junk3)

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
       call readOctreePrivate(grid%octreeRoot,null(),fileFormatted, nOctal, grid, dompiFlag)
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
    write(*,*) myRankGlobal, " read in successfully"

  contains
   
    recursive subroutine readOctreePrivate(thisOctal,parent,fileFormatted, noctal, grid, doMpiFlag)
       ! read in an octal to the grid octree

       implicit none
       type(octal), pointer :: thisOctal
       type(octal), pointer :: parent
       type(gridtype) :: grid

       logical, intent(in)  :: fileFormatted
       integer :: nOctal
       logical :: doMpiFlag
       type(octal), pointer :: thisChild, tempOctal
       integer              :: iChild
       logical :: loopAround !, onThread
       logical :: deleteBranch, found
       integer :: i


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
                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
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
                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell

       end if


!       write(*,*) myrankGlobal, " reading in octal number ", nOctal, " depth ",thisOctal%nDepth, " label ", &
!            thisOctal%label, " split ", grid%splitOvermpi,dompiflag, " nchild ",thisOctal%nChildren
!       write(*,*) "flags ",grid%onekappa,grid%photoionization,molecular,cmf
       thisOctal%oneD = .not.(thisOctal%twoD.or.thisOctal%threeD)


       if (.not.grid%oneKappa) then
!          call readReal2D(thisOctal%kappaAbs,fileFormatted)
!          call readReal2D(thisOctal%kappaSca,fileFormatted)
          call readDouble2D(thisOctal%kappaAbs,fileFormatted)
          call readDouble2D(thisOctal%kappaSca,fileFormatted)
       endif
       if (grid%photoionization) then
          call readDouble2D(thisOctal%ionFrac,fileFormatted)
          call readDouble2D(thisOctal%photoIonCoeff,fileFormatted)
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
          call readDouble2D(thisOctal%atomAbundance, fileFormatted)
          call readDouble3D(thisOctal%atomLevel,fileFormatted)
          call readDouble2D(thisOctal%jnuCont,fileFormatted)
          call readDouble2D(thisOctal%jnuLine,fileFormatted)

          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%microturb
          else
             read(unit=20,iostat=error) thisOctal%microturb
          endif

       endif

       if (hydrodynamics) then
          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%boundaryCondition
          else
             read(unit=20,iostat=error) thisOctal%boundaryCondition
          endif
       endif
          

       if (grid%splitOverMpi) then
          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%mpiThread
          else
             read(unit=20,iostat=error) thisOctal%mpiThread

          endif
       endif

!       write(*,*) myrankGlobal, " rank has mpithread ", thisOctal%mpiThread
       if (myRankGlobal == 0) then
          thisOctal%nChildren = 0
          thisOctal%hasChild = .false.
          thisOctal%indexChild = -999
          goto 666
       endif

       loopAround = .true.

       do while (loopAround)

          if (nThreadsGlobal == 9) then
             if (thisOctal%nDepth == 1) then
                thisOctal%nChildren = 1
                thisOctal%hasChild = .false.
                thisOctal%hasChild(myRankGlobal) = .true.
                thisOctal%indexChild(1) = myRankGlobal
             endif
          endif

          if (nThreadsGlobal == 65) then

             if (thisOctal%nDepth == 2) then
                thisOctal%nChildren = 1
                allocate(thisOctal%child(1:thisOctal%nChildren)) 
                do i = 1, thisOctal%nChildren
                   thisChild => thisOctal%child(i)
                   call readOctreePrivate(thisChild,thisOctal,fileFormatted, nOctal, grid, dompiflag)
                enddo
                goto 666
             endif

             if (thisOctal%nDepth == 3) then

                write(*,*) "rank ", myRankGlobal, " depth ",thisOctal%nDepth, " thread ", thisOctal%mpiThread
                found = .false.
                if (octalOnThread(thisOctal, 1, myRankGlobal)) then
                   found = .true.
                   allocate(thisOctal%child(1:thisOctal%nChildren)) 
                   do iChild = 1, myRankGlobal
                      do i = 1, thisOctal%nChildren
                         thisChild => thisOctal%child(i)
                         call readOctreePrivate(thisChild,thisOctal,fileFormatted, nOctal, grid, dompiflag)
                      enddo
                      if (iChild /= myRankGlobal) then
                         tempOctal => thisChild
                         call deleteOctreeBranch(tempOctal,onlyChildren = .false., &
                              adjustParent = .false. , grid = grid, adjustGridInfo = .false.)
                      endif
                   enddo
                else
                   thisOctal%nChildren = 0
                endif
                if (.not.found) thisOctal%nChildren = 0
                goto 666
             endif
          endif



          loopAround = .false.
          
          if (thisOctal%nChildren > 0) then 
             allocate(thisOctal%child(1:thisOctal%nChildren)) 
             do iChild = 1, thisOctal%nChildren, 1
                thisChild => thisOctal%child(iChild)
!                write(*,*) myRankGlobal, " calling readoctreeprivate. depth ", thisOctal%ndepth, " ichild ", ichild, &
!                     " thisoctal%nchildren ",thisOctal%nChildren

                call readOctreePrivate(thisChild,thisOctal,fileFormatted, nOctal, grid, dompiflag)             

!                write(*,*) myRankGlobal, " has called readoctreeprivate. depth ", thisOctal%ndepth
                deleteBranch = .false.

                if (nThreadsGlobal == 9) then
                   if ((thisOctal%nDepth == 1).and.(thisChild%mpiThread(1) /= myRankGlobal)) deleteBranch = .true.
                endif


             
                if (deleteBranch) then

                   loopAround = .true.

                   tempOctal => thisChild
                   write(*,*) myrankglobal, " split over mpi ", grid%splitOverMpi
                   write(*,*) myRankGlobal, " attempting to delete branch. mpithread ",thisChild%mpithread
                   write(*,*) myrankGlobal, " child depth is ", thisChild%nDepth
                   write(*,*) myrankglobal, " this octal has nchildren set to ", thisOctal%nChildren 
                   call deleteOctreeBranch(tempOctal,onlyChildren = .false., &
                        adjustParent = .false. , grid = grid, adjustGridInfo = .false.)
                   
                   CALL deleteOctal(tempOctal, deleteChildren=.false.,     &
                        adjustParent=.false., grid=grid, & 
                        adjustGridInfo=.false.)
                   deallocate(thisOctal%child)
                   write(*,*) myRankGlobal, " deleted branch ", thisOctal%nChildren
                   write(*,*) myRankGlobal, " thisoctal ndepth ", thisOctal%nDepth
                   
                endif
                
             end do

          endif
       enddo
666    continue
    end subroutine readOctreePrivate



  end subroutine readAMRgridMPI
#endif
  
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




  subroutine fillGridWR137(grid, rCore, mDot, vTerm, beta, temp)

    type(GRIDTYPE) :: grid
    real :: rCore, mDot, vTerm, beta, temp
    integer :: i,j,k
    real :: sinTheta, vel, tot, v0, fac
    type(VECTOR) :: rHat
    integer :: i1, i2
    real :: lineTot, contTot
    contTot = 0.
    lineTot = 0.
    grid%lineEmission = .true.

    do i = 1, grid%nr
       grid%rAxis(i) = log10(rCore) + 2.*real(i-1)/real(grid%nr)
       grid%rAxis(i) = 10.**grid%rAxis(i)
    enddo


    do i = 1, grid%nPhi
       grid%phiAxis(i) = real(twoPi*real(i-1)/real(grid%nPhi-1))
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
             grid%rho(i,j,k) = real((mdot / (fourPi * grid%rAxis(i)**2 * vel))/mHydrogen)
             grid%velocity(i,j,k) = (vel/cSpeed) * rHat
             grid%kappaAbs(i,j,k,1) = (1.e-20*grid%rho(i,j,k))**2
             grid%kappaSca(i,j,k,1) = real(grid%rho(i,j,k) * sigmaE / 4.)
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
                sinTheta = real(sqrt(1.d0 - mu**2))
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




#ifdef HYDRO
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

  subroutine writeDouble1D(variable,fileFormatted)

     real(double),dimension(:),pointer :: variable
     logical, intent(in)       :: fileFormatted
     integer :: error

     if (fileFormatted) then
        if (associated(variable)) then
           write(unit=20,fmt=*,iostat=error) .true.
           if (error /=0) then
             print *, 'Panic: write error in writeDouble1D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) SIZE(variable)
           if (error /=0) then
             print *, 'Panic: write error in writeDouble1D' ; stop
           end if
           write(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then
             print *, 'Panic: write error in writeDouble1D' ; stop
           end if
        else
           write(unit=20,fmt=*,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble1D' ; stop
           end if
        end if
     else
        if (associated(variable)) then
           write(unit=20,iostat=error) .true.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble1D' ; stop
           end if
           write(unit=20,iostat=error) SIZE(variable)
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble1D' ; stop
           end if
           write(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble1D' ; stop
           end if
        else
           write(unit=20,iostat=error) .false.
           if (error /=0) then 
             print *, 'Panic: write error in writeDouble1D' ; stop
          end if
        end if
     end if
             
  end subroutine writeDouble1D
    
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

  subroutine readDouble1D(variable,fileFormatted)

     real(double),dimension(:),pointer :: variable
     logical, intent(in)       :: fileFormatted
     logical                   :: present 
     integer                   :: length
     integer                   :: error

     if (fileFormatted) then
        read(unit=20,fmt=*,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readDouble1D' ; stop
        end if
        if (present) then
           read(unit=20,fmt=*,iostat=error) length
           if (error /=0) then 
             print *, 'Panic: read error in readDouble1D' ; stop
           end if
           allocate(variable(length))
           read(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readDouble1D' ; stop
           end if
        end if
     else
        read(unit=20,iostat=error) present
        if (present) then
           read(unit=20,iostat=error) length
           if (error /=0) then 
             print *, 'Panic: read error in readDouble1D' ; stop
           end if
           allocate(variable(length))
           read(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readDouble1D' ; stop
           end if
        end if
     end if
    
  end subroutine readDouble1D
    
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
          print *, 'Panic: read error in readDouble2D',myrankglobal ; stop
        end if
        if (present) then
           read(unit=20,fmt=*,iostat=error) length1, length2
           if (error /=0) then 
             print *, 'Panic: read error in readDouble2D',myrankglobal ; stop
           end if
           allocate(variable(length1,length2))
           read(unit=20,fmt=*,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readDouble2D',myrankglobal ; stop
           end if
        end if
     else
        read(unit=20,iostat=error) present
        if (error /=0) then 
          print *, 'Panic: read error in readDouble2D',myrankglobal ; stop
        end if
        if (present) then
           read(unit=20,iostat=error) length1, length2
           if (error /=0) then 
             print *, 'Panic: read error in readDouble2D',myrankglobal ; stop
           end if
           allocate(variable(length1,length2))
           read(unit=20,iostat=error) variable
           if (error /=0) then 
             print *, 'Panic: read error in readDouble2D',myrankglobal ; stop
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
#endif

  subroutine initTestAmr(grid)
    use inputs_mod
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
    use inputs_mod
    type(GRIDTYPE) :: grid

    grid%geometry = "proto"

    grid%rCore = rCore / 1.e10
    grid%lCore = fourPi * stefanBoltz * grid%rCore**2 * 1.e20 * teff**4
    grid%rInner = rInner / 1.e10
    grid%rOuter = rOuter / 1.e10
    grid%lineEmission = .false.
  end subroutine initProtoAmr




  !
  ! Print out the basic infomation about the grid.
  !
  
  ! if filename is '*' then it prints on screen.
  subroutine grid_info(thisGrid, filename)
    use amr_mod, only:  countVoxels
    implicit none
    type(gridtype), intent(in) :: thisGrid
    character(LEN=*), intent(in) :: filename
    integer :: UN
    real(double) :: fac
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

    
end module grid_mod

