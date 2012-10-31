! Check the image produced by phaseloop_mod for the disc benchmark
! D. Acreman, October 2012

program check_disc_image

  implicit none

  character(len=*), parameter :: image_file="test_22um.fits"
  logical :: found_file

! FITS file parameters and variables
  integer :: unit, blocksize, status, npixels, nfound
  integer, parameter :: readwrite=0 ! Open file read only 
  integer, parameter :: group=0
  logical :: anynull
  character(len=30 ) :: errtext
  character(len=80 ) :: errmessage

  integer :: npix(2)
  real, allocatable :: image(:,:)
  real :: image_sum

  status = 0

! Firstly make sure that the file exists and fail if it doesn't 
  inquire (file=image_file, exist=found_file)
  if ( found_file ) then
     write(*,*) "Found image file "//image_file
  else
     write(*,*) "Did not find image file "//image_file
     write(*,*) "TORUS: Test failed"
     STOP
  end if

  ! Get a free LUN then open the file
  call ftgiou(unit, status)
  call ftopen(unit, image_file,readwrite,blocksize,status)

! Find the axis sizes
  call ftgknj(unit,"NAXIS",1,2,npix,nfound,status)
  write(*,*) "Axis sizes: ", npix

! Read in the pixel values
  allocate (image(npix(1),npix(2)))
  npixels = npix(1) * npix(2)
  call ftgpve(unit,group,1,npixels,1e-33,image,anynull,status)

! Report basic stats
  image_sum = SUM(image)
  write(*,*) "Sum of pixel values= ", image_sum

! Close the file and free the LUN
  call ftclos(unit, status)
  call ftfiou(unit, status)

! Were there any FITS errors? 
  if (status /= 0) then
     write(*,*) "Exit status of FITS calls is non-zero"

      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext
      !
      !  Read and print out all the error messages on the FITSIO stack
      !
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
         print *,errmessage
         call ftgmsg(errmessage)
      end do

     write(*,*) "TORUS: Test failed"
     STOP
  endif

! We've reached the end so everything is OK
  write(*,*) "TORUS: Test successful"

  deallocate (image)

end program check_disc_image
