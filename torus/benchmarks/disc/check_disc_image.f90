! Check the images produced by phaseloop_mod for the disc benchmark
! D. Acreman, October 2012

program check_disc_image

  implicit none

  integer, parameter :: nimages=2
  character(len=*), parameter :: image_file(nimages)=(/"test_22um.fits ","test2_22um.fits"/)
  logical :: found_file
  real, parameter :: expectedVal(nimages)  = (/5.60e6, 2.14e6/) ! Expected values for sum of pixels
  real, parameter :: tolerance = 0.05 ! Fractional tolerance 
  real :: diff

! FITS file parameters and variables
  integer :: unit, blocksize, status, npixels, nfound
  integer, parameter :: readwrite=0 ! Open file read only 
  integer, parameter :: group=0
  logical :: anynull
  character(len=30 ) :: errtext
  character(len=80 ) :: errmessage

  integer :: npix(2), i
  real, allocatable :: image(:,:)
  real :: image_sum

  status = 0

images:  do i=1,nimages

     ! Firstly make sure that the file exists and fail if it doesn't 
     inquire (file=trim(image_file(i)), exist=found_file)
     if ( found_file ) then
        write(*,*) "Found image file "//image_file(i)
     else
        write(*,*) "Did not find image file "//image_file(i)
        write(*,*) "TORUS: Test failed"
        STOP
     end if

     ! Get a free LUN then open the file
     call ftgiou(unit, status)
     call ftopen(unit, trim(image_file(i)),readwrite,blocksize,status)

     ! Find the axis sizes
     call ftgknj(unit,"NAXIS",1,2,npix,nfound,status)
     write(*,*) "Axis sizes: ", npix

     ! Read in the pixel values
     allocate (image(npix(1),npix(2)))
     npixels = npix(1) * npix(2)
     call ftgpve(unit,group,1,npixels,1e-33,image,anynull,status)

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

     ! Report basic stats
     image_sum = SUM(image) 
     write(*,*) "Sum of pixel values= ", image_sum

     diff = abs(image_sum - expectedVal(i)) / expectedVal(i)
     write(*,'(a,f6.1,a)') "Difference from expected value= ", diff*100.0, " %"
     write(*,'(a,f6.1,a)') "Tolerance= ", tolerance*100.0, " %"
     if (diff > tolerance) then
        write(*,*) "TORUS: Test failed"
        STOP
     endif

     deallocate (image)

  end do images

! We've reached the end so everything is OK
  write(*,*) "TORUS: Test successful"

end program check_disc_image
