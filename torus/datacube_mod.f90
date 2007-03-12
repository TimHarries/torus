module datacube_mod

  use kind_mod
  use vector_mod
  use messages_mod

  implicit none

  type DATACUBE
     character(len=80) :: label
     character(len=10) :: vUnit ! units for velocity
     character(len=10) :: xUnit ! units for space
     character(len=10) :: IntensityUnit ! units for intensity

     integer, pointer :: nsubpixels(:,:,:) ! contains resolution information 
     integer, pointer :: converged(:,:,:)  ! contains convergence information (should take 1 or 0)
     real(double), pointer :: weight(:,:)     ! Weighting for integration (used to find spectra)
     integer :: nx 
     integer :: ny
     integer :: nv

     real(double) :: obsdistance ! observation distance for use with instrument functions
     real(double), pointer :: xAxis(:)
     real(double), pointer :: yAxis(:)
     real(double), pointer :: vAxis(:)
     real(double), pointer :: intensity(:,:,:)
  end type DATACUBE

contains

! Initialises cube - sets intensity for cube to 0 
  subroutine initCube(thisCube, nx, ny, nv)
    type(DATACUBE) :: thisCube
    integer :: nx, ny, nv

    thisCube%nx = nx
    thisCube%ny = ny
    thisCube%nv = nv
    allocate(thisCube%xAxis(1:nx))
    allocate(thisCube%yAxis(1:ny))
    allocate(thisCube%vAxis(1:nv))
    allocate(thisCube%intensity(1:nx,1:ny,1:nv))
    allocate(thisCube%nsubpixels(1:nx,1:ny,1:nv))
    allocate(thisCube%converged(1:nx,1:ny,1:nv))
    allocate(thisCube%weight(1:nx,1:ny))

    thisCube%intensity = 0.d0
    thisCube%nsubpixels = 0.d0
    thisCube%converged = 0
    thisCube%weight = 1

  end subroutine initCube

! Set spatial axes for datacube - Equally spaced (linearly) between min and max
  subroutine addSpatialAxes(cube, xMin, xMax, yMin, yMax)
    type(DATACUBE) :: cube
    real(double) :: xMin, xMax, yMax, yMin
    integer :: i

    do i = 1, cube%nx
       cube%xAxis(i) = xmin + (xmax-xmin)*dble(i-1)/dble(cube%nx)
    enddo
    Do i = 1, cube%ny
       cube%yAxis(i) = ymin + (ymax-ymin)*dble(i-1)/dble(cube%ny)
    enddo
  end subroutine addSpatialAxes

! Set velocity axis for datacube - Equally spaced (linearly) between min and max
  subroutine addVelocityAxis(cube, vMin, vMax)
    type(DATACUBE) :: cube
    real(double) :: vMin, vMax
    integer :: i

    do i = 1, cube%nv
       cube%vAxis(i) = vmin + (vmax-vmin)*dble(i-1)/dble(cube%nv)
    enddo
  end subroutine addVelocityAxis

  subroutine plotDataCube(cube, device, withspec)
    type(DATACUBE) :: cube
    character(len=*) :: device
    logical, optional :: withSpec
    integer :: i, j, k
    integer :: pgbegin
    real, allocatable :: image(:,:)
    integer :: nx, ny
    real :: iMin, iMax
    real :: tr(6)
    real :: dx, dy
    real(double), allocatable :: spec(:)
    real :: vxs, vxe, vys, vye
    real :: x1, x2, y1, y2
    real(double) :: sMax, Smin
    real :: range
    integer :: nstep 
    logical :: doSpec

    if (present(withSpec)) then
       doSpec = withSpec
    else
       doSpec = .true.
    endif

    nx = cube%nx
    ny = cube%ny

    dx = (cube%Xaxis(nx) - cube%xAxis(1))/real(nx-1)
    dy = (cube%yAxis(ny) - cube%yAxis(1))/real(ny-1)

    tr(1) = cube%xAxis(1)-dx
    tr(2) = dx
    tr(3) = 0.
    tr(4) = cube%yAxis(1)-dy
    tr(5) = 0.
    tr(6) = dy

    allocate(image(1:nx, 1:ny))

    do i = 1, nx
       do j = 1, ny
          if ( sum(cube%intensity(i,j,1:cube%nv)) .gt. 0) then
             image(i,j) = log10(sum(cube%intensity(i,j,1:cube%nv)))
          else
             image(i,j) = -30.
          endif
       enddo
    enddo

! Useful for visualising this 'flattened' cube

    open(42, file="image.dat",status="unknown",form="formatted")
    do i = 1, nx
          write(42, *) image(:,i)
    enddo
    close(42)

    iMin = MINVAL(image(1:nx,1:ny))
    iMax = MAXVAL(image(1:nx,1:ny))
    write(*,*) "min/max",imin,imax
    i =  pgbegin(0,device,1,1)
    write(*,*) "opening ",trim(device),i
    call pgvport(0.1, 0.9, 0.1, 0.9)
    call pgwnad(real(cube%xAxis(1))-dx/2., real(cube%xAxis(nx))+dx/2., &
         real(cube%yAxis(1))-dy/2., real(cube%yAxis(ny))+dy/2.)

    call palette(3)
 
    call pgimag(image, nx, ny, 1, nx, 1, ny, imin, imax, tr)
       
    call pgbox('bcnst',0.0,0,'bcnst',0.0,0)

    call pgqvp(0, x1, x2, y1, y2)

    allocate(spec(1:cube%nv))

    if (doSpec) then

       nStep = nx / 5
       
       smin = 1.e30
       smax = -1.e30
       do i = 1, nx-1, nstep
          do j = 1, ny-1, nstep
             call getSpectrum(cube, i, i+nstep-1, j, j+nstep-1, spec)
             sMax = MAX(sMax,MAXVAL(spec))
             sMin = MIN(sMin,MINVAL(spec))
          enddo
       enddo

       range = sMax - sMin
       sMin = sMin - 0.2*range
       sMax = sMax + 0.2*range
       write(*,*) "min/max",smin,smax
       
       do i = 1, nx-1, nstep
          do j = 1, ny-1, nstep
             call getSpectrum(cube, i, i+nstep-1, j, j+nstep-1, spec)
             vxs = x1 + (x2-x1)*real(i-1)/real(nx)
             vxe = x1 + (x2-x1)*real(i+nstep-1)/real(nx)
             vys = y1 + (y2-y1)*real(j-1)/real(ny)
             vye = y1 + (y2-y1)*real(j+nstep-1)/real(ny)
             call pgsci(3)
             call pgvport(vxs, vxe, vys, vye)
             !          call pgbox('bc',0.0,0,'bc',0.0,0)
             call pgwindow(real(cube%vAxis(1))-0.1, &
                  real(cube%vAxis(cube%nv))+0.1, &
                  real(smin), real(smax))
             
             !write(*,*) spec(1:cube%nv)
             call pgline(cube%nv, real(cube%vAxis), &
                  real(spec))
             
             call pgsci(1)
          enddo
       enddo
    end if
    call pgend
    call getSpectrum(cube, 1, cube%nx, 1, cube%ny, spec)
    
    write(*,*) "SPECTRUM"
    open(43, file="imagespectrum.dat",status="unknown",form="formatted")
    
    do k = 1, cube%nv
       write(*,*) cube%vAxis(k), spec(k)
       write(43,*) cube%vAxis(k), spec(k)
    enddo
    close(43)

  end subroutine plotDataCube
  subroutine getSpectrum(cube, ix1, ix2, iy1, iy2, spec)
    type(DATACUBE) :: cube
    integer :: ix1, ix2, iy1, iy2
    real(double) ::  spec(:)
    integer :: i,j, n
    n = 0
    spec = 0.d0
    do i = ix1, ix2
       do j = iy1, iy2
          n = n + 1
          spec(1:cube%nv) = spec(1:cube%nv) + cube%intensity(i,j,1:cube%nv)
       enddo
    enddo
!    spec = spec / dble(n)
  end subroutine getSpectrum

  subroutine getWeightedSpectrum(cube, ix1, ix2, iy1, iy2, spec)
    type(DATACUBE) :: cube
    integer :: ix1, ix2, iy1, iy2
    real(double) ::  spec(:) , tot
    integer :: i,j
    tot = 0.d0
    spec = 0.d0
    do i = ix1, ix2
       do j = iy1, iy2
          tot =  tot + cube%weight(i,j)
          spec(1:cube%nv) = spec(1:cube%nv) + cube%intensity(i,j,1:cube%nv) * cube%weight(i,j)
       enddo
    enddo
!    spec = spec / dble(tot)
  end subroutine getWeightedSpectrum

  subroutine convolveCube(cube, beamSize)
    type(DATACUBE) :: cube
    real(double) :: beamSize ! beamsize in arcsec
    real(double), allocatable :: newArray(:,:)
    real(double) :: r, rinArcSec, weight,fac 
    integer :: ix, iy, iv, i, j
    integer :: i1, j1
    real(double) :: sigma, dx, dy, tot, flux, background
    real(double) :: deltaX, deltaY

    sigma = beamsize /2.35d0

    dx = 3600.d0*((cube%xAxis(2) - cube%xAxis(1))/(cube%obsDistance/1.d10))*180.d0/pi
    dy = 3600.d0*((cube%yAxis(2) - cube%yAxis(1))/(cube%obsDistance/1.d10))*180.d0/pi

    allocate(newArray(1:cube%nx, 1:cube%ny))
    call writeInfo("Convolving data cube with beam size", TRIVIAL)
    write(*,*) "cube%obsdist",cube%obsDistance/pctocm

    do iv = 1, cube%nv

       background = cube%intensity(cube%nx,cube%ny,iv)
       newArray = 0.d0

       do ix = 1, cube%nx
          do iy = 1, cube%ny
             tot = 0.d0
             do i = ix - cube%nx, ix + cube%nx
                do j = iy - cube%ny, iy + cube%ny
                   if ((i > 0).and.(i<=cube%nx).and.(j > 0).and.(j <= cube%ny)) then
                      flux = cube%intensity(i,j,iv)
                   else
                      flux = background
                   endif

                   deltaX = dble(ix-i)*dx
                   deltaY = dble(iy-j)*dy
                   rInArcSec = sqrt(deltaX**2 + deltaY**2)
             
                   fac = (1.d0/(twoPi*sigma**2))*exp(-0.5d0*(rInArcSec**2/sigma**2))*dx*dy

                   newArray(ix,iy) = newArray(ix, iy) + flux*fac
                   tot = tot + fac
                enddo
             enddo
!             write(*,*) "Weight",tot
          enddo
!          write(*,*) newArray,tot
       enddo

       cube%intensity(1:cube%nx, 1:cube%ny, iv) = newArray(1:cube%nx, 1:cube%ny)!/dble(cube%nx*cube%ny)
    enddo
    deallocate(newArray)
    call writeInfo("Done.",TRIVIAL)
  end subroutine convolveCube
    
subroutine TranslateCubeIntensity(cube,constant)

  type(DATACUBE) :: cube
  real(double) :: constant
  
  cube%intensity = cube%intensity + constant
    
end subroutine TranslateCubeIntensity

end module datacube_mod

