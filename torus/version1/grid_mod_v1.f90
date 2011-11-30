  subroutine fillGridDisk(grid, rhoNought, rCore, rInner, rOuter, height, &
       mCore, diskTemp)


    implicit none
    type(GRIDTYPE) :: grid
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
                rVec = rVec / dble(r)
                vVec = rVec  .cross. spinAxis
                if (modulus(vVec) /= 0.) then
                   call normalize(vVec)
                   vel = sqrt(bigG * mCore / r) / cSpeed
                   vVec = dble(vel)  * vVec 
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


  subroutine fillGridTorus(Grid, rho,  rTorus, rOuter)


    implicit none
    type(GRIDTYPE) :: Grid
    real(double) :: phi
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
          aVec = dble(router) * (cos(phi) * rHat)
          bVec = dble(router) * (sin(phi) * torusAxis)
          rUp = r0Vec + (aVec + bVec)
          rDown = r0Vec + (aVec - bVec)
          call locate(grid%xAxis,grid%nx,real(rUp%x), i1)
          call locate(grid%yAxis,grid%ny,real(rUp%y), i2)
          call locate(grid%zAxis,grid%nz,real(rUp%z), i3)
          call locate(grid%zAxis,grid%nz,real(rDown%z), i4)

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
             thisVec = rotateX(rVec, dble(rotationAngle))
             thisVec = thisVec - grid%starPos2
             thisVec = rotateZ(thisVec, dble(-deflectionAngle))
             thisDirection = rotateX(direction(2), dble(rotationAngle))
             thisVec = thisVec + grid%starPos2 + depthVec
             thisDirection = rotateZ(thisDirection, dble(deflectionAngle))
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
             thisVec = rotateX(rVec, dble(rotationAngle))
             thisVec = thisVec - grid%starPos2
             thisVec = rotateZ(thisVec, dble(-deflectionAngle))
             thisDirection = rotateX(direction(j), dble(rotationAngle))
             thisVec = thisVec + grid%starPos2 + depthVec
             thisDirection = rotateZ(thisDirection, dble(deflectionAngle))
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
             grid%velocity(i,j,k) = grid%velocity(i,j,k) + (dble(vel)*vVec)
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

  subroutine fillGridRolf(grid, mdot, vterm, o6width)
    use utils_mod, only: logInt
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
    real(double) :: fac, fac2, r, theta, phi
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
                call locate(rArray, ir(itheta,iphi), real(r), i1)
                fac = logint(real(r), rArray(i1), rArray(i1+1), yArray(i1), yArray(i1+1))
                fac2 = logint(real(r), rArray(i1), rArray(i1+1), yArray2(i1), yArray2(i1+1))
             endif
             if (r > rArray(ir(itheta,iphi))) then
                fac = 1.
                fac2 = yArray2(ir(iTheta,iPhi))
             endif
             if (r < rArray(1)) then
                fac = 1.e-10
                fac2 = yArray2(1)
             endif
             grid%kappaSca(i,j,k,1) = max(1.e-30_db,grid%kappaSca(i,j,k,1) * fac)
             grid%kappaAbs(i,j,k,1) = max(1.e-30_db,grid%kappaAbs(i,j,k,1) * fac)
             grid%kappaScaRed(i,j,k,1) = max(1.e-30_db,grid%kappaScaRed(i,j,k,1) * fac)
             grid%kappaAbsRed(i,j,k,1) = max(1.e-30_db,grid%kappaAbsRed(i,j,k,1) * fac)
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
!!$
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
          norm2 = rotateZ(norm, dble(ang))
          sVec = VECTOR(x,y,z)
          sVec1 = sVec - (dr/2.d0)*norm2
          sVec2 = sVec + (dr/2.d0)*norm2
          do j = 1, nr
             sVec = sVec1 + (dble(j-1)/dble(nr-1))*(sVec2-sVec1)
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
    real :: rOuter, rInner, rCore
    type(VECTOR) :: rVec
    real :: phi, r
    real :: height
    integer, parameter :: nPhi = 200
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


  subroutine fillGridDonati(grid, resonanceLine)
    use utils_mod, only: hunt
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

  subroutine fillGridDisk2(grid, rhoNought, rCore, rInner, rOuter, height, &
       mCore, diskTemp)


    implicit none
    type(GRIDTYPE) :: grid
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
                rVec = rVec / dble(r)
                vVec = rVec  .cross. spinAxis
                if (modulus(vVec) /= 0.) then
                   call normalize(vVec)
                   vel = sqrt(bigG * mCore / r) / cSpeed
                   vVec = dble(vel)  * vVec 
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

  subroutine writeAMRgridOld(filename,fileFormatted,grid)
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
      use octal_mod, only: writeattributepointer, writeattributestatic
       ! writes out an octal from the grid octree

       type(octal), intent(in), target :: thisOctal
       logical, intent(in)             :: fileFormatted
       type(gridtype) :: grid
       type(octal), pointer :: thisChild
       integer              :: iChild

       call writeAttributeStatic(20, thisOctal%nDepth, fileFormatted)
       call writeAttributeStatic(20, thisOctal%nChildren, fileFormatted)
       call writeAttributeStatic(20, thisOctal%indexChild, fileFormatted)
       call writeAttributeStatic(20, thisOctal%hasChild, fileFormatted)
       call writeAttributeStatic(20, thisOctal%centre, fileFormatted)
       call writeAttributeStatic(20, thisOctal%rho, fileFormatted)
       call writeAttributeStatic(20, thisOctal%temperature, fileFormatted)
       call writeAttributeStatic(20, thisOctal%label, fileFormatted)
       call writeAttributeStatic(20, thisOctal%subcellSize, fileFormatted)
       call writeAttributeStatic(20, thisOctal%threeD, fileFormatted)
       call writeAttributeStatic(20, thisOctal%twoD, fileFormatted)
       call writeAttributeStatic(20, thisOctal%maxChildren, fileFormatted)
       call writeAttributeStatic(20, thisOctal%cylindrical, fileFormatted)
       call writeAttributeStatic(20, thisOctal%splitAzimuthally, fileFormatted)
       call writeAttributeStatic(20, thisOctal%phi, fileFormatted)
       call writeAttributeStatic(20, thisOctal%dphi, fileFormatted)
       call writeAttributeStatic(20, thisOctal%r, fileFormatted)
       call writeAttributeStatic(20, thisOctal%parentSubcell, fileFormatted)
       call writeAttributeStatic(20, thisOctal%inStar, fileFormatted)
       call writeAttributeStatic(20, thisOctal%inFlow, fileFormatted)
       call writeAttributeStatic(20, thisOctal%velocity, fileFormatted)
       call writeAttributeStatic(20, thisOctal%cornervelocity, fileFormatted)

       call writeAttributePointer(20, thisOctal%chiLine, fileFormatted)
       call writeAttributePointer(20, thisOctal%etaLine, fileFormatted)
       call writeAttributePointer(20, thisOctal%etaCont, fileFormatted)
       call writeAttributePointer(20, thisOctal%biasLine3D, fileFormatted)
       call writeAttributePointer(20, thisOctal%biasCont3D, fileFormatted)
       call writeAttributePointer(20, thisOctal%probDistLine, fileFormatted)
       call writeAttributePointer(20, thisOctal%probDistCont, fileFormatted)
       call writeAttributePointer(20, thisOctal%ne, fileFormatted)
       call writeAttributePointer(20, thisOctal%nH, fileFormatted)
       call writeAttributePointer(20, thisOctal%nTot, fileFormatted)
       call writeAttributePointer(20, thisOctal%dustType, fileFormatted)


       
!       if (fileFormatted) then 
!          write(iostat=error,fmt=*,unit=20) thisOctal%nDepth, thisOctal%nChildren,&
!                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
!                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
!                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
!                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
!                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
!                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
!                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
!                  thisOctal%subcellSize,thisOctal%threed, thisOctal%twoD,    &
!                  thisOctal%maxChildren, thisOctal%dustType,  &
!                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
!                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
!       else
!          write(iostat=error,unit=20) thisOctal%nDepth, thisOctal%nChildren, &
!                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
!                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
!                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
!                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
!                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
!                  thisOctal%probDistCont, thisOctal%Ne, thisOctal%nh, thisOctal%nTot,      &
!                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
!                  thisOctal%subcellSize, thisOctal%threeD, thisOctal%twoD,   &
!                  thisOctal%maxChildren, thisOctal%dustType, &
!                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
!                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
!       end if 
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
          call writeDouble2D(thisOctal%jnu,fileFormatted)
          call writeDouble2D(thisOctal%bnu,fileFormatted)
          call writeReal1D(thisOctal%molAbundance,fileFormatted)
          call writeDouble1D(thisOctal%nh2,fileFormatted)
          call writeDouble1D(thisOctal%microturb,fileFormatted)
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

       if (grid%splitOverMpi) then
          if (fileformatted) then
             write(unit=20,iostat=error,fmt=*) thisOctal%mpiThread
          else
             write(unit=20,iostat=error) thisOctal%mpiThread
          endif
       endif
       
       if (thisOctal%nChildren > 0) then 
          do iChild = 1, thisOctal%nChildren, 1
             thisChild => thisOctal%child(iChild)
             call writeOctreePrivate(thisChild,fileFormatted,grid)
          end do
       end if

    end subroutine writeOctreePrivate
    
    
  end subroutine writeAMRgridOld

  subroutine readAMRgridOld(filename,fileFormatted,grid)
    ! reads in a previously saved 'grid' for an adaptive mesh geometry  

    use input_variables, only: geometry,dipoleOffset,amr2dOnly,statEq2d, molecular, cmf, hydrodynamics
    use unix_mod, only: unixGetenv
    implicit none

    character(len=*)            :: filename
    character(len=80)            :: message
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
    absolutePath = " "
    error = 0
    junk => null()
    junk2 => null()
    junk3 => null()

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

    write(message,'(a,a)') "Reading populations file from: ",trim(filename)
    call writeInfo(message,TRIVIAL)
    write(message,'(a,i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') ' - data file written at: ', &
                          timeValues(1),'/',timeValues(2),'/',&
                          timeValues(3),'  ',timeValues(5),':',timeValues(6)
    call writeInfo(message,TRIVIAL)
                          
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

!    deallocate(junk, junk2, junk3)
   
    Call readClumps(fileFormatted)
    
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
    
    close(unit=20)    grid%geometry = trim(geometry)
    grid%dipoleOffset = dipoleOffset
    grid%amr2donly = amr2donly
    grid%statEq2d = statEq2d

  contains
   
    recursive subroutine readOctreePrivate(thisOctal,parent,fileFormatted, noctal, grid)
       ! read in an octal to the grid octree
      use octal_mod, only: readattributepointer, readattributestatic

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


       call readAttributeStatic(20, thisOctal%nDepth, fileFormatted)
       call readAttributeStatic(20, thisOctal%nChildren, fileFormatted)
       call readAttributeStatic(20, thisOctal%indexChild, fileFormatted)
       call readAttributeStatic(20, thisOctal%hasChild, fileFormatted)
       call readAttributeStatic(20, thisOctal%centre, fileFormatted)
       call readAttributeStatic(20, thisOctal%rho, fileFormatted)
       call readAttributeStatic(20, thisOctal%temperature, fileFormatted)
       call readAttributeStatic(20, thisOctal%label, fileFormatted)
       call readAttributeStatic(20, thisOctal%subcellSize, fileFormatted)
       call readAttributeStatic(20, thisOctal%threeD, fileFormatted)
       call readAttributeStatic(20, thisOctal%twoD, fileFormatted)
       call readAttributeStatic(20, thisOctal%maxChildren, fileFormatted)
       call readAttributeStatic(20, thisOctal%cylindrical, fileFormatted)
       call readAttributeStatic(20, thisOctal%splitAzimuthally, fileFormatted)
       call readAttributeStatic(20, thisOctal%phi, fileFormatted)
       call readAttributeStatic(20, thisOctal%dphi, fileFormatted)
       call readAttributeStatic(20, thisOctal%r, fileFormatted)
       call readAttributeStatic(20, thisOctal%parentSubcell, fileFormatted)
       call readAttributeStatic(20, thisOctal%inStar, fileFormatted)
       call readAttributeStatic(20, thisOctal%inFlow, fileFormatted)
       call readAttributeStatic(20, thisOctal%velocity, fileFormatted)
       call readAttributeStatic(20, thisOctal%cornervelocity, fileFormatted)

       call readAttributePointer(20, thisOctal%chiLine, fileFormatted)
       call readAttributePointer(20, thisOctal%etaLine, fileFormatted)
       call readAttributePointer(20, thisOctal%etaCont, fileFormatted)
       call readAttributePointer(20, thisOctal%biasLine3D, fileFormatted)
       call readAttributePointer(20, thisOctal%biasCont3D, fileFormatted)
       call readAttributePointer(20, thisOctal%probDistLine, fileFormatted)
       call readAttributePointer(20, thisOctal%probDistCont, fileFormatted)
       call readAttributePointer(20, thisOctal%ne, fileFormatted)
       call readAttributePointer(20, thisOctal%nH, fileFormatted)
       call readAttributePointer(20, thisOctal%nTot, fileFormatted)
       call readAttributePointer(20, thisOctal%dustType, fileFormatted)

       
!       if (fileFormatted) then
!          read(unit=20,fmt=*,iostat=error) thisOctal%nDepth, thisOctal%nChildren,  &
!                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
!                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
!                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
!                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
!                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
!                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
!                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
!                  thisOctal%subcellSize, thisOctal%threeD, thisOctal%twoD,   &
!                  thisOctal%maxChildren, thisOctal%dustType, &
!                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
!                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
!       else
!          read(unit=20,iostat=error) thisOctal%nDepth, thisOctal%nChildren,  &
!                  thisOctal%indexChild, thisOctal%hasChild, thisOctal%centre,& 
!                  thisOctal%rho, thisOctal%velocity,thisOctal%cornerVelocity,&
!                  thisOctal%temperature, thisOctal%chiLine,thisOctal%etaLine,&
!                  thisOctal%etaCont, thisOctal%biasLine3D,                   &
!                  thisOctal%biasCont3D, thisOctal%probDistLine,              &
!                  thisOctal%probDistCont, thisOctal%Ne,  thisOctal%nh,thisOctal%nTot,      &
!                  thisOctal%inStar, thisOctal%inFlow, thisOctal%label,       &
!                  thisOctal%subcellSize, thisOctal%threeD,  thisOctal%twoD,  &
!                  thisOctal%maxChildren, thisOctal%dustType, &
!                  thisOctal%cylindrical, thisOctal%splitAzimuthally, thisOctal%phi, &
!                  thisOctal%dPhi, thisOctal%r, thisOctal%parentSubcell
!
!       end if
!       
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
          call readDouble2D(thisOctal%jnu,fileFormatted)
          call readDouble2D(thisOctal%bnu,fileFormatted)
          call readReal1D(thisOctal%molAbundance,fileFormatted)

          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%nh2
             read(unit=20,iostat=error,fmt=*) thisOctal%microturb
          else
             read(unit=20,iostat=error) thisOctal%nh2
             read(unit=20,iostat=error) thisOctal%microturb
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
             write(unit=20,iostat=error,fmt=*) thisOctal%boundaryCondition
          else
             write(unit=20,iostat=error) thisOctal%boundaryCondition
          endif
       endif


       if (grid%splitOverMpi) then
          if (fileformatted) then
             read(unit=20,iostat=error,fmt=*) thisOctal%mpiThread
          else
             read(unit=20,iostat=error) thisOctal%mpiThread

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
 
  end subroutine readAMRgridOld

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



  subroutine fillGridDonati2(grid, resonanceLine, misc)
    use utils_mod, only: hunt
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


  !  
  !  This subrouine recursively writes out the distance from the center of the
  !  given plane,  with a give value of the 3rd dimension, and the correcsponding 
  !  value stored in a octal to a file specified by the logical unit number. 
  !  
  recursive subroutine radial_profile(thisOctal, name, plane, val_3rd_dim, luout, center, grid)
    use octal_mod, only: subcellCentre
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
    type(VECTOR), intent(in) :: center   ! position of the center of the root node
    !
    !
    type(octal), pointer  :: child 
    type(VECTOR) :: rvec
    real :: value

    integer :: subcell, i
    real :: xp, yp, xm, ym, zp, zm
!    real(double) :: kabs, ksca
!    integer :: ilam
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

