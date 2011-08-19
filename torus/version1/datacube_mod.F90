  subroutine convolveCube(cube, beamSize)
    use constants_mod, only: pcTocm, pi, twoPi
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

       background = cube%intensity(cube%nx,cube%ny,iv) ! This needs to be in flux units not Intensity
       newArray = 0.d0

       do ix = 1, cube%nx
          do iy = 1, cube%ny
             tot = 0.d0
             do i = ix - cube%nx, ix + cube%nx
                do j = iy - cube%ny, iy + cube%ny
                   if ((i > 0).and.(i<=cube%nx).and.(j > 0).and.(j <= cube%ny)) then
                      flux = cube%intensity(i,j,iv) ! This needs to be in flux units not Intensity
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
       ! For flux calculations this needs to be in flux units not Intensity
    enddo
    deallocate(newArray)
    call writeInfo("Done.",TRIVIAL)
  end subroutine convolveCube

#ifdef USECFITSIO

  subroutine readDataCube(thisCube)

    implicit none
    
    type(DATACUBE), intent(out) :: thisCube
    integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels, hdutype, hdu, junk
    character(len=512) :: filename
    character(len=80) :: comment
    real :: nullval
    integer, dimension(4) :: naxes
    logical :: anynull
    character(len=80) :: keyword

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

    call ftgcrd(unit,keyword,status)

    write(*,*) "Keyword", keyword

    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of datacube.fits.gz file HDU',hdu,'. Exiting.'
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)
    nbuffer=npixels
    ! read_image
    call ftgpve(unit,group,firstpix,nbuffer,nullval,thisCube%intensity,anynull,status)


    !  Read keywords from the header.
    call FTGKYJ(unit,"LABEL", junk,comment,status)
    thisCube%label = comment(1:10)
    call FTGKYJ(unit,"VUNIT", junk,comment,status)
    thisCube%vUnit = comment(1:10)
    call FTGKYJ(unit,"XUNIT", junk,comment,status)
     thisCube%xUnit = comment(1:10)
    call FTGKYJ(unit,"IUNIT", junk,comment,status)
    thisCube%IntensityUnit = comment(1:10)
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


    return

  end subroutine readDataCube
#endif

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


