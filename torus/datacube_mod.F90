module datacube_mod

  use kind_mod
  use vector_mod
  use messages_mod

  implicit none

  type TELESCOPE
     
     character(len=10) :: label
     real :: diameter
     real :: beamsize ! in arcsecs
     
  end type TELESCOPE

  type DATACUBE
     TYPE(TELESCOPE) :: telescope

     character(len=80) :: label
     character(len=10) :: vUnit ! units for velocity
     character(len=10) :: xUnit ! units for space
     character(len=10) :: IntensityUnit ! units for intensity
     character(len=10) :: FluxUnit ! units for flux

     integer, pointer :: nsubpixels(:,:,:) ! contains resolution information 
     integer, pointer :: converged(:,:,:)  ! contains convergence information (should take 1 or 0)
     real(double),pointer :: weight(:,:)     ! Weighting for integration (used to find spectra)
     integer :: nx 
     integer :: ny
     integer :: nv

     real(double) :: obsdistance ! observation distance for use with instrument functions
     real(double), pointer :: xAxis(:)
     real(double), pointer :: yAxis(:)
     real(double), pointer :: vAxis(:)
     real(double), pointer :: intensity(:,:,:)
     real(double), pointer :: flux(:,:,:)
  end type DATACUBE

contains

  subroutine writeDataCube(thisCube, filename)

    implicit none
    
    type(DATACUBE), intent(in) :: thisCube

    character(len=*) :: filename

#ifdef USECFITSIO
    integer :: status,unit,blocksize,bitpix,naxis
    integer, dimension(5) :: naxes
    integer :: group,fpixel,nelements
    logical :: simple, extend
    
    status=0
    
    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou ( unit, status )
    
    !  Create the new empty FITS file.
    blocksize=1
    call ftinit(unit,trim(filename),blocksize,status)

    !  Initialize parameters about the FITS image
    simple=.true.
    extend=.true.
    group=1
    fpixel=1

    ! 1st HDU : intensity
    bitpix=-64
    naxis=3
    naxes(1)=thisCube%nx
    naxes(2)=thisCube%ny
    naxes(3)=thisCube%nv
    nelements=naxes(1)*naxes(2)*naxes(3)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    
    !  Write the array to the FITS file.
    call ftpprd(unit,group,fpixel,nelements,thisCube%flux,status)


    !  Write keywords to the header.
    call ftpkyj(unit,'LABEL',1,thisCube%label,status) 
    call ftpkyj(unit,'VUNIT',1,thisCube%vUnit,status)
    call ftpkyj(unit,'XUNIT',1,thisCube%xUnit,status)
    call ftpkyj(unit,'IUNIT',1,thisCube%IntensityUnit,status)
    call ftpkyd(unit,'DISTANCE',thisCube%obsdistance,-3,'observation distance',status)

    ! 2nd HDU : converged
    call FTCRHD(unit, status)
    bitpix=32
    naxis=3
    naxes(1)=thisCube%nx
    naxes(2)=thisCube%ny
    naxes(3)=thisCube%nv
    nelements=naxes(1)*naxes(2)*naxes(3)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    !  Write the array to the FITS file.
    call ftpprj(unit,group,fpixel,nelements,thisCube%converged,status)

    ! 3rd HDU : weight
    call FTCRHD(unit, status)
    bitpix=-64
    naxis=2
    naxes(1)=thisCube%nx
    naxes(2)=thisCube%ny
    nelements=naxes(1)*naxes(2)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    !  Write the array to the FITS file.
    call ftpprd(unit,group,fpixel,nelements,thisCube%weight,status)

    ! 4th HDU : xAxis
    call FTCRHD(unit, status)
    bitpix=-64
    naxis=1
    naxes(1)=thisCube%nx     
    nelements=naxes(1)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    !  Write the array to the FITS file.
    call ftpprd(unit,group,fpixel,nelements,thisCube%xAxis,status)

    ! 5th HDU : yAxis
    call FTCRHD(unit, status)
    bitpix=-64
    naxis=1
    naxes(1)=thisCube%ny    
    nelements=naxes(1)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    !  Write the array to the FITS file.
    call ftpprd(unit,group,fpixel,nelements,thisCube%yAxis,status)

    ! 6th HDU : vAxis
    call FTCRHD(unit, status)
    bitpix=-64 
    naxis=1
    naxes(1)=thisCube%nv        
    nelements=naxes(1)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    !  Write the array to the FITS file.
    call ftpprd(unit,group,fpixel,nelements,thisCube%vAxis,status)

    ! 7th HDU : nsubpixels
    call FTCRHD(unit, status)
    bitpix=32
    naxis=3
    naxes(1)=thisCube%nx
    naxes(2)=thisCube%ny
    naxes(3)=thisCube%nv
    nelements=naxes(1)*naxes(2)*naxes(3)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    !  Write the array to the FITS file.
    call ftpprj(unit,group,fpixel,nelements,thisCube%nsubpixels,status)


    !  Close the file and free the unit number.
    call ftclos(unit, status)
    call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
    if (status > 0) then
       call print_error(status)
    end if
    
#endif
    return
    
  end subroutine writeDataCube

  !**********************************************************************

  subroutine readDataCube(thisCube)

    implicit none
    
#ifdef USECFITSIO
    type(DATACUBE), intent(out) :: thisCube
    integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels, hdutype, hdu, junk
    character(len=512) :: filename
    character(len=80) :: comment
    real :: nullval
    integer, dimension(4) :: naxes
    logical :: anynull

    status=0
    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou(unit,status)

    filename="datacube.fits.gz"
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)

    group=1
    firstpix=1
    nullval=-999
     

    ! 1st HDU : intensity
    hdu=1
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%intensity,anynull,status)


    !  Read keywords from the header.
    call FTGKYJ(unit,"LABEL", junk,comment,status)
    thisCube%label = comment
    call FTGKYJ(unit,"VUNIT", junk,comment,status)
    thisCube%vUnit = comment
    call FTGKYJ(unit,"XUNIT", junk,comment,status)
     thisCube%xUnit = comment
    call FTGKYJ(unit,"IUNIT", junk,comment,status)
    thisCube%IntensityUnit = comment
    call FTGKYD(unit,"DISTANCE", thisCube%obsdistance,comment,status)

    ! 2nd HDU : converged
    hdu=2
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)
    nbuffer=npixels
    ! read_image
    call ftgpvj(unit,group,firstpix,nbuffer,nullval,thisCube%converged,anynull,status)

    ! 3rd HDU : weight
    hdu=3
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
    if (nfound /= 2) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.' 
       stop
    endif
    npixels=naxes(1)*naxes(2)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%weight,anynull,status)

    ! 4th HDU : xAxis
    hdu=4
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'     
       stop
    endif
    npixels=naxes(1)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%xAxis,anynull,status)
    
    ! 5th HDU : yAxis
    hdu=5
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'   
       stop
    endif
    npixels=naxes(1)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%yAxis,anynull,status)
    
    ! 6th HDU : vAxis
    hdu=6
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.' 
       stop
    endif
    npixels=naxes(1)
    nbuffer=npixels
    ! read_image
    call ftgpvd(unit,group,firstpix,nbuffer,nullval,thisCube%vAxis,anynull,status)
    
    ! 7th HDU : intensity
    hdu=7
    call ftmahd(unit,hdu,hdutype,status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU', hdu,'. Exiting.'
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)
    nbuffer=npixels
    ! read_image
    call ftgpvj(unit,group,firstpix,nbuffer,nullval,thisCube%nsubpixels,anynull,status)   

#else
    type(DATACUBE) :: thisCube

#endif
    return

  end subroutine readDataCube

  !**********************************************************************


  subroutine print_error(status)
    ! PRINT_ERROR prints out the FITSIO error messages to the user.
    
    integer status
#ifdef USECFITSIO
    character ( len = 30 ) errtext
    character ( len = 80 ) errmessage

    !  Check if status is OK (no error); if so, simply return.
    if (status <= 0) then
       return
    end if

    !  Get the text string which describes the error
    call ftgerr(status,errtext)
    print *,'FITSIO Error Status =',status,': ',errtext

    !  Read and print out all the error messages on the FITSIO stack
    call ftgmsg(errmessage)
    do while (errmessage .ne. ' ')
       print *,errmessage
       call ftgmsg(errmessage)
    end do
#endif    

    return
  end subroutine print_error

  !***********************************************************

! Initialises cube - sets intensity for cube to 0 
  subroutine initCube(thisCube, nx, ny, nv, mytelescope)
    type(DATACUBE) :: thisCube
    type(TELESCOPE), optional :: mytelescope
    integer :: nx, ny, nv

    if(present(mytelescope)) then
       
       thisCube%telescope = mytelescope

    else

       thisCube%telescope%label=" "
       thisCube%telescope%diameter = 1.
       thisCube%telescope%beamsize = 1d10 !arcsecs - set large to give unconvolved image

    endif

    thisCube%label=" "
    thisCube%vUnit=" "
    thisCube%xUnit=" "
    thisCube%IntensityUnit=" "
    thisCube%FluxUnit=" "

    thisCube%nx = nx
    thisCube%ny = ny
    thisCube%nv = nv
    allocate(thisCube%xAxis(1:nx))
    allocate(thisCube%yAxis(1:ny))
    allocate(thisCube%vAxis(1:nv))
    allocate(thisCube%intensity(1:nx,1:ny,1:nv))
    allocate(thisCube%flux(1:nx,1:ny,1:nv))
    allocate(thisCube%nsubpixels(1:nx,1:ny,1:nv))
    allocate(thisCube%converged(1:nx,1:ny,1:nv))
    allocate(thisCube%weight(1:nx,1:ny))

    thisCube%intensity = 0.d0
    thisCube%flux = 0.d0
    thisCube%nsubpixels = 0.d0
    thisCube%converged = 0
    thisCube%weight = 1.d0

  end subroutine initCube

! Set spatial axes for datacube - Equally spaced (linearly) between min and max
  subroutine addSpatialAxes(cube, xMin, xMax, yMin, yMax)
    type(DATACUBE) :: cube
    real(double) :: xMin, xMax, yMax, yMin, dx, dy
    integer :: i

    dx = (xMax - xMin)/dble(cube%nx)
    dy = (yMax - yMin)/dble(cube%ny)
    do i = 1, cube%nx
       cube%xAxis(i) = xmin + dx/2.d0 + dble(i-1)*dx 

    enddo

    do i = 1, cube%ny
       cube%yAxis(i) = ymin + dy/2.d0 + dble(i-1)*dx 
    enddo

  end subroutine addSpatialAxes

! Set velocity axis for datacube - Equally spaced (linearly) between min and max
  subroutine addVelocityAxis(cube, vMin, vMax)
    type(DATACUBE) :: cube
    real(double) :: vMin, vMax, dv
    integer :: i
    
    dv = (vMax - vMin) / dble(cube%nv)

    do i = 1, cube%nv
       cube%vAxis(i) = vmin + dv/2.d0 + dv * dble(i-1)
    enddo
  end subroutine addVelocityAxis

  subroutine plotDataCube(cube, device, withspec, gotPixels, plotflux, twoPanels)
    type(DATACUBE) :: cube
    character(len=*) :: device
    character(len=40) :: message
    logical, optional :: withSpec, twoPanels, gotPixels, plotFlux
    logical :: doTwoPanels, doplotFlux
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
    logical :: doSpec !,PixelCount


    write(*,*) "Starting plotdatacube ",trim(device)
    if (present(withSpec)) then
       doSpec = withSpec
    else
       doSpec = .true.
    endif

   if (present(plotflux)) then
       doplotflux = plotflux
    else
       doplotflux = .false.
    endif

    if (present(twopanels)) then
       dotwopanels = twopanels
    else
       doTwoPanels = .false.
    endif

    nx = cube%nx
    ny = cube%ny

    dx = cube%Xaxis(2) - cube%xAxis(1)
    dy = cube%yAxis(2) - cube%yAxis(1)

    tr(1) = cube%xAxis(1)-dx
    tr(2) = dx
    tr(3) = 0.
    tr(4) = cube%yAxis(1)-dy
    tr(5) = 0.
    tr(6) = dy

    allocate(image(1:nx, 1:ny))

    if(doplotflux) then
       write(*,*) "Plotting Flux?", doplotflux
       do i = 1, nx
          do j = 1, ny
             if (sum(cube%flux(i,j,1:cube%nv)) .gt. 0.) then
                image(i,j) = log10(sum(cube%flux(i,j,1:cube%nv)))
             else
                image(i,j) = -30.
             endif
          enddo
       enddo
    else
       write(*,*) "Plotting Intensity!"
       do i = 1, nx
          do j = 1, ny
             if (sum(cube%intensity(i,j,1:cube%nv)) .gt. 0) then
                image(i,j) = log10(sum(cube%intensity(i,j,1:cube%nv)))
             else
                image(i,j) = -30.
             endif
          enddo
       enddo
    endif

    allocate(spec(1:cube%nv))

! Useful for visualising this 'flattened' cube

    open(42, file="image.dat",status="unknown",form="formatted")
    do i = 1, nx
       write(42, *) image(:,i)
    enddo
    close(42)

!    if (present(gotPixels)) then
!       PixelCount = .true.
!    endif

!    if(PixelCount) then
!       doSpec = .false.
!       do i = 1, nx
!          do j = 1, ny
!             image(i,j) = sum(cube%nsubpixels(i,j,1:cube%nv))
!          enddo
!       enddo
!    endif
    
!    open(41, file="pixeldensity.dat",status="unknown",form="formatted")
!    do i = 1, nx
!       write(41, *) image(:,i)
!    enddo
!    close(41)


    iMin = MINVAL(image(1:nx,1:ny))
    iMax = MAXVAL(image(1:nx,1:ny))

    imin = imax-5.
    write(*,*) "min/max",imin,imax
    i =  pgbegin(0,device,1,1)
    write(*,*) "opening ",trim(device),i

    write(message,*) 'Telescope: ',cube%telescope%label
    call pgmtxt('T',2.0,0.5,0.5,message)

    write(message,*) 'Diameter: ',cube%telescope%diameter/100.,' metres'
    call pgmtxt('T',1.0,0.5,0.5,message)

    write(message,*) 'Beamsize: ',cube%telescope%beamsize,' arcseconds'
    call pgmtxt('T',0.0,0.5,0.5,message)

    write(message,*) 'Max: 10**',iMax
    call pgmtxt('B',0.0,0.0,0.0,message)

    write(message,*) 'Min: 10**',iMin
    call pgmtxt('B',1.0,0.0,0.0,message)

    call pgvport(0.1, 0.9, 0.1, 0.9)

    if (doTwoPanels) then
       call pgvport(0.1, 0.95, 0.1, 0.95)
       call getSpectrum(cube, 1, cube%nx, 1, cube%ny, spec)
       spec = spec / spec(1)
       call pgwindow(real(cube%vAxis(1)), real(cube%vAxis(cube%nv)), &
            0., 10.)
       call pgsci(3)
       call pgline(cube%nv, real(cube%vAxis), real(spec))
       call pgsci(1)
       call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
       call pgvport(0.65, 0.9, 0.65, 0.9)

    endif

    call pgqvp(0, x1, x2, y1, y2)

!    call pgvport(0.1, 0.9, 0.1, 0.9)

    call pgwnad(real(cube%xAxis(1))-dx/2., real(cube%xAxis(nx))+dx/2., &
         real(cube%yAxis(1))-dy/2., real(cube%yAxis(ny))+dy/2.)

    call palette(3)
 
    call pgimag(image, nx, ny, 1, nx, 1, ny, imin, imax, tr)
       
!    call pgbox('bcnst',10.0,10,'bcnst',10.0,10)

    call pgqvp(0, x1, x2, y1, y2)

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
    write(message,*) "Telescope Name:", cube%telescope%label
    call writeinfo(message,FORINFO)
    call pgend
    call getSpectrum(cube, 1, cube%nx, 1, cube%ny, spec)
    
    write(*,*) "SPECTRUM"
    open(43, file="imagespectrum.dat",status="unknown",form="formatted")
    
    do k = 1, cube%nv
!       write(*,*) cube%vAxis(k), spec(k)
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
          spec(1:cube%nv) = spec(1:cube%nv) + cube%flux(i,j,1:cube%nv)
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
    real(double) :: rrinArcSec, fac 
    integer :: ix, iy, iv, i, j
    real(double) :: sigma, sigma2, dx, dy, tot, flux, background
    real(double) :: deltaX, deltaY

    sigma = beamsize/2.35d0 ! changed from 2.35 (FWHM) as a result of reading IRAM-30m paper
    sigma2 = sigma**2
    dx = 3600.d0*((cube%xAxis(2) - cube%xAxis(1))/(cube%obsDistance/1.d10))*180.d0/pi
    dy = 3600.d0*((cube%yAxis(2) - cube%yAxis(1))/(cube%obsDistance/1.d10))*180.d0/pi

    allocate(newArray(1:cube%nx, 1:cube%ny))
    call writeInfo("Convolving data cube with beam size", TRIVIAL)
    write(*,*) "cube%obsdist",cube%obsDistance/pctocm, "dx",dx,"dy",dy

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
                   rrInArcSec = deltaX**2 + deltaY**2
             
                   fac = (1.d0/(twoPi*sigma2))*exp(-0.5d0*(rrInArcSec/sigma2))*dx*dy

                   newArray(ix,iy) = newArray(ix, iy) + flux*fac
                   tot = tot + fac
                enddo
             enddo
          enddo
          !write(*,*) "Weight",tot
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

subroutine freeDataCube(thiscube)
  type(DATACUBE) :: thiscube


    if (associated(thisCube%xAxis)) deallocate(thiscube%xAxis)
    if (associated(thisCube%yAxis)) deallocate(thiscube%yAxis)
    if (associated(thisCube%vAxis)) deallocate(thiscube%vAxis)

    if (associated(thisCube%intensity)) deallocate(thiscube%intensity)
    if (associated(thisCube%flux)) deallocate(thiscube%flux)
    if (associated(thisCube%nsubpixels)) deallocate(thiscube%nSubpixels)
    if (associated(thisCube%converged)) deallocate(thiscube%converged)
    if (associated(thisCube%weight)) deallocate(thiscube%weight)

  end subroutine freeDataCube

! This subroutine was raising a warning about intensityUnit being used but not set. As this 
! subroutine is never called I have commented it out (Dave Acreman, 16/11/07).
!
!!$subroutine ConvertUnits(cube)
!!$  type(DATACUBE) :: cube
!!$  character(len=10) :: intensityUnit ! System 1 - cgs | System 2 - SI | System 3 - Janskys
!!$
!!$  ! - convert to brightness temperature - c**2/2*nu**2*k
!!$
!!$  select case(cube%intensityUnit)
!!$
!!$  case('erg/cm2/Hz')
!!$       if(intensityUnit .eq. 'W/m2') cube%intensity = cube%intensity * 1d-3
!!$       if(intensityUnit .eq. 'Janskys') cube%intensity = cube%intensity * 1d-23
!!$!       if(intensityUnit .eq. 'Tbright') cube%intensity = cube%intensity *
!!$!       if(intensityUnit .eq. 'Tantenna') cube%intensity = cube%intensity *
!!$       
!!$       cube%intensityUnit = 'erg/cm2/Hz'
!!$
!!$    case('W/m2')
!!$       if(intensityUnit .eq. 'erg/cm2/Hz') cube%intensity = cube%intensity * 1d3
!!$       if(intensityUnit .eq. 'Janskys') cube%intensity = cube%intensity * 1d-26
!!$!       if(intensityUnit .eq. 'Tbright') cube%intensity = cube%intensity * 
!!$!       if(intensityUnit .eq. 'Tantenna') cube%intensity = cube%intensity *
!!$       
!!$       cube%intensityUnit = 'W/cm2'
!!$
!!$    case('Janskys')
!!$       if(intensityUnit .eq. 'erg/cm2/Hz') cube%intensity = cube%intensity * 1d23
!!$       if(intensityUnit .eq. 'W/m2') cube%intensity = cube%intensity * 1d26
!!$!       if(intensityUnit .eq. 'Tbright') cube%intensity = cube%intensity * 
!!$!       if(intensityUnit .eq. 'Tantenna') cube%intensity = cube%intensity * 
!!$       
!!$       cube%intensityUnit = 'Janskys'
!!$
!!$!    case('Tbright')
!!$!       if(intensityUnit .eq. 'erg/cm2/Hz') cube%intensity = cube%intensity * 
!!$!       if(intensityUnit .eq. 'W/m2') cube%intensity = cube%intensity * 
!!$!       if(intensityUnit .eq. 'Janskys') cube%intensity = cube%intensity * 
!!$!       if(intensityUnit .eq. 'Tantenna') cube%intensity = cube%intensity * 
!!$       
!!$!       cube%intensityUnit = 'Tbright'
!!$
!!$!    case('Tantenna')
!!$!       if(intensityUnit .eq. 'erg/cm2/Hz') cube%intensity = cube%intensity * 
!!$!       if(intensityUnit .eq. 'Janskys') cube%intensity = cube%intensity * 
!!$!       if(intensityUnit .eq. 'Tbright') cube%intensity = cube%intensity * 
!!$!       if(intensityUnit .eq. 'W/m2') cube%intensity = cube%intensity * 
!!$       
!!$!       cube%intensityUnit = 'Tantenna'
!!$
!!$       case default
!!$          cube = cube
!!$        
!!$end select
!!$
!!$end subroutine ConvertUnits  

end module datacube_mod

