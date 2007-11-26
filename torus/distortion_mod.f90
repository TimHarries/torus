! 
! Module containing subroutines that distort the grid in some way
! Examples include adding rotational velocity fields and latitudinal
! density structures.
!

! written by tjh

! v1.0  on 13/08/99

module distortion_mod

  use gridtype_mod          ! opacity grid
  use grid_mod              ! opacity grid routines
  use constants_mod         ! physical constants
  use vector_mod            ! vector math
  use blob_mod              ! blob module

  implicit none

  public

contains

  ! this is a test distortion, which is a big blob placed at
  ! two stellar radii on the positive x-axis

  subroutine distortGridTest(grid)
  
    type(GRIDTYPE) :: grid                   ! the opacity grid
    integer, parameter :: maxBlobs = 1       ! just one blob
    type(BLOBTYPE) :: blob(maxBlobs)
    real :: r                                ! radius
    integer :: i                             ! counter

    ! inform the user that they've chosen the test

    write(*,'(a)') "Distorting grid with test blob..."

    ! two stellar radii

    r = grid%rAxis(1)*2.
    call locate(grid%rAxis, grid%nr, r, i)

    ! set up the blob

    blob(1)%position = VECTOR(r, 0., 0.)
    blob(1)%velocity = grid%velocity(i,grid%nmu/2,1)
    blob(1)%contrast = 100                            ! high contrast
    blob(1)%radius = grid%rAxis(1)/2.                 ! big
    blob(1)%inUse = .true.

    ! distort the grid

    call distortGridWithBlobs(grid, maxBlobs, blob)  

  end subroutine distortGridTest


  ! This distorts the grid by one (or many) spiral arms

  subroutine distortGridSpiral(Grid, Vrot, nSpiral)

    type(GRIDTYPE) :: Grid                        ! the opacity grid
    real :: vRot                                  ! rotational speed
    integer :: i, j, k                            ! counters
    integer :: nSpiral                            ! no of spiral arms
    integer :: iSpiral                            
    real, allocatable :: facGrid(:,:,:)           ! density contrast grid
    real :: x, w
    integer :: i1,i2
    real :: muStart,muEnd                    ! latitudinal extent of spiral
    real :: thickness                        ! thickness of arm
    real :: r1
    integer ::  j1,j2
    real :: phi, theta, r
    real :: phaseOffset, spiralPhase         ! phase offsets of arms
    type(VECTOR) :: rHat, perp, zAxis, startVec, endVec  ! vectors
    type(VECTOR) :: thisVec, posvec, vVec
    real :: rMin, rMax, r0, phi0, phi1
    real :: tot
    real :: x0, x1, y0, y1, dx, dy
    real :: vPhi

    write(*,*) "Distorting grid with spiral..."

    ! the z-axis (the direction of the angular momentum vector)

    zAxis = VECTOR(0.,0.,1.)

    ! only works for polar grids

    if (Grid%cartesian) then
       write(*,*) "! Need to use polar grid for spiral geometry"
       stop
    endif

    ! 60 < theta < 120 for the spiral

    muStart = 0.5
    muEnd = -0.5
    
    ! thickness of 0.1 stellar radii

    thickness = 0.1

    ! no phase offset

    phaseOffset = 0.

    ! radial limits - no need to go beyond 30 stellar radii

    rMin = grid%rAxis(1)
    rMax = 30. * grid%rAxis(1)

    ! allocate the density contrast grid

    allocate(facGrid(1:grid%nr, 1:grid%nMu, 1:grid%nPhi))

    ! should be mostly ones

    facGrid = 1.

    ! find the grid indices corresponding to the the latitudinal extent

    call locate(grid%muAxis,grid%nMu,muStart,j1)
    call locate(grid%muAxis,grid%nMu,muEnd,j2)

    ! loop over grid

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi

             ! the azimuthal speed
             
             vPhi = (vRot/cSpeed/(grid%rAxis(i)/rMin)) * &
                                     sqrt(1.-grid%muAxis(j)**2)

             ! current velocity - purely radial

             vVec = grid%velocity(i,j,k)

             rHat = vVec
             call normalize(rHat)
             
             ! radial unit vector crossed with the z-axis will give
             ! the azimuthal direction

             perp = rHat .cross. zAxis

             ! keep radial velocity but add a phi component

             grid%velocity(i,j,k) = vVec + vPhi * perp

          enddo
       enddo
    enddo  ! loop over grid


    ! do this for each spiral arm

    do iSpiral = 1, nSpiral

       ! work out the starting phase of this arm

       spiralPhase = twoPi * real(iSpiral-1)/real(nSpiral)

       ! loop over all radii

       do i = 1 , grid%nr-1


          r0 = grid%rAxis(i)
          r1 = grid%rAxis(i+1)

          tot = 0.
          do j = 1, 10000
             r = rMin + (r0 - rMin) * real(j-1)/9999.
             call locate(grid%rAxis, grid%nr, r, i1)
             w = modulus(grid%velocity(i1,1,1))/modulus(grid%velocity(grid%nr,1,1))
             tot = tot + ((r0-rMin)/(10000.*rMin)) * (1. - (rMin/r)**2) / w
          enddo
          phi0 = tot * vRot/modulus(grid%velocity(grid%nr,1,1))/cSpeed
          phi0 = phi0 + phaseOffSet + spiralPhase

          x0 = r0 * cos(phi0)
          y0 = r0 * sin(phi0)

          tot = 0.
          do j = 1, 10000
             r = rMin + (r1 - rMin) * real(j-1)/9999.
             call locate(grid%rAxis, grid%nr, r, i1)
             w = modulus(grid%velocity(i1,1,1))/modulus(grid%velocity(grid%nr,1,1))
             tot = tot + ((r1-rMin)/(10000.*rMin)) * (1. - (rMin/r)**2) / w
          enddo
          phi1 = -tot * vRot/modulus(grid%velocity(grid%nr,1,1))/cSpeed
          phi1 = phi1 + phaseOffSet + spiralPhase

          x1 = r1 * cos(phi1)
          y1 = r1 * sin(phi1)

          dx = x1 - x0
          dy = y1 - x0

          rHat = VECTOR(dx,dy,0.)
          call normalize(rHat)
          perp = rHat .cross. zAxis

          posVec = r0*VECTOR(cos(phi0),sin(phi0),0.)


          startVec = posVec - ((thickness/2.)*r0)*perp
          endVec = posVec + ((thickness/2.)*r0)*perp

          do i1 = 1,100
             x = real(i1-1)/99.
             thisVec = startVec + (x*(endVec-startVec))
             call getPolar(thisVec, r, theta, phi)
             call locate(grid%rAxis,grid%nr,r,i2)
             call locate(grid%phiAxis,grid%nPhi,phi,k)
             do j = min(j1,j2), max(j1,j2)
                if (facGrid(i2,j,k) == 1.) then
                   facGrid(i2,j,k) = 20.
                endif
             enddo
          enddo
       enddo

    enddo

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             grid%etaLine(i,j,k) = grid%etaLine(i,j,k) * facGrid(i,j,k)**2
             grid%chiLine(i,j,k) = grid%chiLine(i,j,k) * facGrid(i,j,k)**2
             grid%kappaSca(i,j,k,1) = grid%kappaSca(i,j,k,1) * facGrid(i,j,k)
          enddo
       enddo
    enddo



    deallocate(facGrid) ! free grid
    write(*,*) "done."

  end subroutine distortGridSpiral


  ! this subroutine distorts the grid by a rotational velocity field

  subroutine distortRotation(Grid, Vrot)

    type(GRIDTYPE) :: Grid                      ! the opacity grid
    real :: vRot                                ! the rotational speed
    integer :: i, j, k                          ! counters
    real :: rMin                                ! minimum radius
    type(VECTOR) :: rHat, zAxis,  vVec, phiHat  ! vectors
    real :: vPhi                                ! azimuthal speed

    ! the rotation axis
    
    zAxis=VECTOR(0.,0.,1.)

    ! only works for polar grids

    if (Grid%cartesian) then
       write(*,*) "! Need to use polar grid for rotational geometry"
       stop
    endif


    rMin = grid%rAxis(1)

    ! inform the user

    write(*,'(a)') "Adding rotational velocity fields..."


    ! loop over grid

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi

             ! azimuthal speed

             vPhi = (vRot/cSpeed/(grid%rAxis(i)/rMin)) * &
                                       sqrt(1.-grid%muAxis(j)**2)


             ! initial field is purely radial

             vVec = grid%velocity(i,j,k)
             rHat = vVec
             call normalize(rHat)

             ! find the perpendicular, which points in the phiHat direction
             
             phiHat = rHat .cross. zAxis
             if (modulus(phiHat) /= 0.) then
                call normalize(phiHat)
             endif

             grid%velocity(i,j,k) = vVec + vPhi * phiHat

          enddo
       enddo
    enddo   ! loop over grid

    ! tell user we've finished

    write(*,'(a)') "Done."

  end subroutine distortRotation

  subroutine distortWRdisk(grid)
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    real :: height,z,fac


    height = 2.*grid%rAxis(1)

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi

             z = grid%rAxis(i)*grid%muAxis(j)
             fac = max(exp(-abs(z/height)),1.e-10)
             grid%etaLine(i,j,k) = grid%etaLine(i,j,k) * fac
             grid%etaCont(i,j,k) = grid%etaCont(i,j,k) * fac
             grid%kappaSca(i,j,k,1) = grid%kappaSca(i,j,k,1) * fac

          enddo
       enddo
    enddo
  end subroutine distortWRdisk
             


  ! this subroutine knocks out a couple of the radial shells

  subroutine distortRaman(Grid)

    type(GRIDTYPE) :: Grid                      ! the opacity grid
    integer :: i, j, k, i1
    real :: vRaman


    vRaman = 500.*1.e5/cSpeed



    if (Grid%cartesian) then
       write(*,*) "! Need to use polar grid for this"
       stop
    endif


!    call locate(grid%rAxis,grid%nr, vRaman, i)

    i = 14

    write(*,'(a,i2)') "Doing Raman distortion...",i


    ! loop over grid

    do i1= i-1,i+1
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             grid%chiLine(i1,j,k) = grid%chiLine(i1,j,k)*1.e-4
             grid%etaLine(i1,j,k) = grid%etaLine(i1,j,k)*1.e-4
          enddo
       enddo
    enddo   ! loop over grid

    ! tell user we've finished

    write(*,'(a)') "Done."

  end subroutine distortRaman

  subroutine distortStrom(grid, hotSourcePosition, concave, flatten, zScale, coolStarPosition, &
       ramanDist)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: hotSourcePosition
    logical :: concave, convex, flatten
    integer :: i,j,k
    real :: r, rStrom, rMin, rMax
    type(VECTOR) :: rVec, rHat, coolDirection, zAxis, vHat
    type(VECTOR) :: dr, oldrVec, coolStarPosition
    real :: phi, cosPhi, x,y,z,zScale,fac,v,rDisk,mHot
    logical :: bipolar, spiral, circle, plane, reverse, disk
    real :: xCen, xDash, yDash, mu
    real :: ang, openingAng, cosAng, rho, rHole
    logical, allocatable :: done(:,:,:)
    integer :: i1,i2,i3
    character(len=*) :: ramanDist

    mHot = 0.6 * mSol
    rDisk = modulus(hotSourcePosition-coolStarPosition)/2.

    zAxis = VECTOR(0.,0.,1.)

    coolDirection = VECTOR(-1.,0.,0.)

    concave = .false.
    convex = .false.
    flatten = .false.
    bipolar = .false.
    spiral = .false.
    circle = .false.
    plane = .false.
    reverse = .false.
    disk = .false.

    zScale = grid%rCore

    if (flatten) then
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                fac = abs(grid%zAxis(k)/zScale)
                fac = exp(-fac)
                grid%kappaSca(i,j,k,1) = grid%kappaSca(i,j,k,1) * fac
             enddo
          enddo
       enddo
    endif

 

    select case (ramanDist)

    case("convex")

       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))-hotSourcePosition
                rHat = rVec
                call normalize(rHat)
                r = modulus(rVec)
                cosPhi = rHat .dot. coolDirection
                phi = acos(cosPhi)
                rStrom = modulus(hotSourcePosition-coolStarPosition)/ (1. + 2.*cos(2.*phi))
                if (((rVec%x > 0.).or.((rStrom > 0.) .and.(r < rStrom)).or.(rStrom < 0.))) then
                   grid%kappaAbs(i,j,k,1) = 1.e-30
                   grid%kappaSca(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                endif
             enddo
          enddo
       enddo

    case("concave")

       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))-hotSourcePosition
                rHat = rVec
                call normalize(rHat)
                r = modulus(rVec)
                cosPhi = rHat .dot. coolDirection
                phi = acos(cosPhi)
                fac = 2./3.
                rStrom = modulus(hotSourcePosition-coolStarPosition)/ (fac*(1. + 2.*cos(phi)))
                if (((rStrom < 0.) .or.(r < rStrom))) then
                   grid%kappaAbs(i,j,k,1) = 1.e-30
                   grid%kappaSca(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                endif
             enddo
          enddo
       enddo

    case("hole")

       rHole = 0.5*modulus(coolStarPosition - hotSourcePosition)
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))-hotSourcePosition
                r = modulus(rVec)
                if (r < rHole) then
                   grid%kappaAbs(i,j,k,1) = 1.e-30
                   grid%kappaSca(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                endif
             enddo
          enddo
       enddo





    case ("spiral")
       allocate(done(1:grid%nx,1:grid%ny,1:grid%nz))
       done = .false.
       oldrVec = VECTOR(0.,0.,0.)
       do j = 1, 100
          mu = 2.*real(j-1)/99.-1.
          do i = 1, 1000       
             ang = twoPi*real(i)/999.
             r = 2.*ang/twoPi*modulus(hotSourcePosition-coolStarPosition)
             fac = sqrt(1.-mu**2)
             rVec = VECTOR(fac*r*cos(ang),-fac*r*sin(ang),r*mu)
             dr = rVec - oldrVec 
             rVec = rVec + coolStarPosition
             call normalize(dr)
             if (insideGrid(grid, rVec)) then
                call hunt(grid%xAxis,grid%nx,rvec%x, i1)
                call hunt(grid%yAxis,grid%ny,rvec%y, i2)
                call hunt(grid%zAxis,grid%nz,rvec%z, i3)
                if (.not.done(i1,i2,i3)) then
                   rho = 1.e11
                   grid%kappaSca(i1,i2,i3,1) = grid%kappaSca(i1,i2,i3,1) * 100.
                   grid%velocity(i1,i2,i3) = (50.e5/cSpeed) * dr
                   done(i1,i2,i3) = .true.
                endif
                oldrVec = rVec
             endif
          enddo
       enddo
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                             if (.not.done(i,j,k)) grid%kappaSca(i,j,k,1) = grid%kappaSca(i,j,k,1) * 1.e-30
             enddo
          enddo
       enddo
       deallocate(done)

    case ("circle")
       allocate(done(1:grid%nx,1:grid%ny,1:grid%nz))
       done = .false.
       oldrVec = VECTOR(0.,0.,0.)
       rMin = 0.1*modulus(hotSourcePosition - coolStarPosition)
       rMax = 0.5*modulus(hotSourcePosition - coolStarPosition)
       do k = 1,100
          do i = 1, 1000       
             ang = twoPi*real(i)/999.
             r = rMin + (rMax-rMin)*real(k-1)/99.
             rVec = VECTOR(r*cos(ang),r*sin(ang),0.)
             dr = rVec .cross. zAxis
             rVec = hotSourcePosition - rVec
             call normalize(dr)
             call hunt(grid%xAxis,grid%nx,rvec%x, i1)
             call hunt(grid%yAxis,grid%ny,rvec%y, i2)
             call hunt(grid%zAxis,grid%nz,rvec%z, i3)
             if (.not.done(i1,i2,i3)) then
                rho = 1.e10
                grid%kappaSca(i1,i2,i3,1) = rho * sigmaE * (34. + 6.6)
                grid%velocity(i1,i2,i3) = (50.e5/cSpeed) * dr
                done(i1,i2,i3) = .true.
             endif
             oldrVec = rVec
          enddo
       enddo
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                if (.not.done(i,j,k)) then
                   grid%kappaSca(i,j,k,1) = grid%kappaSca(i,j,k,1) * 1.e-10
                   grid%velocity(i,j,k) = VECTOR(0.,0.,0.)
                endif
             enddo
          enddo
       enddo
       
       
       deallocate(done)


    case ("disk")
       allocate(done(1:grid%nx,1:grid%ny,1:grid%nz))
       done = .false.
       do i = 1, 1000       
          ang = twoPi*real(i-1)/999.
          do j = 1, 1000
             rVec = (real(j-1)/999.)*VECTOR(rDisk*cos(ang),rDisk*sin(ang),0.)
             r = modulus(rVec)
             vHat = rVec .cross. zAxis
             if (modulus(vHat) /= 0.) call normalize(vHat)
             rVec = hotSourcePosition + rVec
             call normalize(dr)
             call hunt(grid%xAxis,grid%nx,rvec%x, i1)
             call hunt(grid%yAxis,grid%ny,rvec%y, i2)
             call hunt(grid%zAxis,grid%nz,rvec%z, i3)
             if (.not.done(i1,i2,i3)) then
                rho = 1.e10
                grid%kappaSca(i1,i2,i3,1) = rho * sigmaE * (34. + 6.6)
                if (r /= 0.) then
                   v = sqrt(bigG * mHot / r)/cSpeed
                else
                   v = 0.
                endif
                grid%velocity(i1,i2,i3) = v * vHat
                done(i1,i2,i3) = .true.
             endif
          enddo
       enddo
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                if (.not.done(i,j,k)) then
                   grid%kappaSca(i,j,k,1) = grid%kappaSca(i,j,k,1) * 1.e-10
                   grid%velocity(i,j,k) = VECTOR(0.,0.,0.)
                endif
             enddo
          enddo
       enddo
       deallocate(done)

    case ("plane")
       xCen = 0.5*(hotSourcePosition%x + coolStarPosition%x)
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                x = grid%xAxis(i)
                y = grid%yAxis(j)
                z = grid%zAxis(k)
                yDash = y
                xDash = xCen !0.1*(y-Xcen)
                if (x > xDash) then
                   grid%kappaAbs(i,j,k,1) = 1.e-30
                   grid%kappaSca(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                endif
             enddo
          enddo
       enddo

    case("reverse")
       do i = 1, grid%nz
          do j = 1, grid%ny
             do k = 1, grid%nz
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j), grid%zAxis(k)) - coolStarPosition
                r = modulus(rVec)
                rHat = rVec
                call normalize(rHat)
                if (r > grid%rCore) then
                   grid%velocity(i,j,k) = (50.e5*sqrt(grid%rCore/r)/cSpeed) * rHat
                endif
             enddo
          enddo
       enddo

    case ("bipolar")
       openingAng = 20.*degToRad
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz                
                rho = 1.e9
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))-hotSourcePosition
                rHat = rVec
                call normalize(rHat)                
                cosAng = abs(rHat .dot. zAxis)
                ang = acos(cosAng)

                if (ang < openingAng) then
                   grid%velocity(i,j,k) = &
                        (100.e5/cSpeed)*rHat

!                   fac = 1. ! (abs(grid%zAxis(k)/grid%zAxis(grid%nz)))**2
!                   grid%kappaAbs(i,j,k,1) = rho * fac * 0.
!                   grid%kappaSca(i,j,k,1) = rho * fac * (34.+6.6)*sigmaE
!                   grid%kappaAbsRed(i,j,k,1) = 0.
!                   grid%kappaAbsRed(i,j,k,1) = 0.
                endif

             enddo
          enddo
       enddo
       case DEFAULT
          write(*,'(a,a)') "! Unrecognised distortion type (ramandist): ",trim(ramanDist)

    end select

    if (bipolar) then
      openingAng = 20.*degToRad
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz                
                rho = 1.e9
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))-hotSourcePosition
                rHat = rVec
                call normalize(rHat)                
                cosAng = abs(rHat .dot. zAxis)
                ang = acos(cosAng)

                if (ang < openingAng) then
                   grid%velocity(i,j,k) = &
                        (40.e5/cSpeed)*rHat

                   grid%kappaAbs(i,j,k,1) = rho * fac * 1.e-30
                   grid%kappaSca(i,j,k,1) = rho * fac * (34.+6.6)*sigmaE
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                   grid%kappaAbsRed(i,j,k,1) = 1.e-30
                endif

             enddo
          enddo
       enddo
    endif

  end subroutine distortStrom

  subroutine distortWindCollision(grid, momRatio, binarySep)
    
    type(GRIDTYPE) :: grid
    real :: momratio
    real :: binarySep, binarySep10
    real :: stagPoint
    real :: openingAngle
    real :: rotationAngle
    real :: xDist, yDist, y
    type(VECTOR) :: rVec
    real :: t1, t2, t3
    integer :: i1, i2, i3, i, j, k, iMin

    binarySep10 = binarySep/1.e10

    stagPoint = binarySep10 * sqrt(momRatio) / (1. + sqrt(momRatio))

    openingAngle = 2.1*(1 - (momRatio**(2./5.))/4.)*momRatio**(1./3.)



    if (grid%cartesian) then
       write(*,'(a)') "! wind collision distortion only works on polar grid"
       stop
    endif
    write(*,'(a)') "Removing a wind-wind collision volume..."


    write(*,'(a,f5.1)') "Stagnation point is at wind radius (rCore): ",(binarySep10-stagPoint)/grid%rCore
    write(*,'(a,f5.1)') "Opening angle (degrees): ",openingAngle*radtodeg

    rVec = VECTOR(binarySep10-stagPoint,0.,0.)
    call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
    iMin = i1

    do i = 1, 100
       rotationAngle = twoPi * real(i-1)/99.
       do j = iMin,grid%nr
          xDist = max(grid%rAxis(j),binarySep10-stagPoint)
          yDist = (xDist -(binarySep10-stagPoint)) * tan(openingAngle)
          do k = 1, 100
             y = yDist * real(k-1)/99.
             rVec = VECTOR(xDist, y, 0.)
             rVec = rotateX(rVec, rotationAngle)
             call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
!             write(*,*) i1, i2, i3
             grid%etaLine(i1,i2,i3) = 1.e-30
             grid%etaCont(i1,i2,i3) = 1.e-30
             grid%chiLine(i1,i2,i3) = 1.e-30
             grid%kappaAbs(i1,i2,i3,1) = 1.e-30
             grid%kappaSca(i1,i2,i3,1) = 1.e-30
          enddo
       enddo
    enddo
  end subroutine distortWindCollision

end module distortion_mod
