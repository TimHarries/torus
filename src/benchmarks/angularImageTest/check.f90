! Checking code for angular image test 
! Build with nagfor -lcfitsio -L/Users/acreman/cfitsio_nagfor check.f90 
! D. Acreman, July 2013

program check

  implicit none

  integer, parameter :: db = selected_real_kind(15, 307)

! Constants
  real(db), parameter :: pi = 3.1415926535897932_db
  real(db), parameter :: mHydrogen = 1.6733d-24      ! g
  real(db), parameter :: kpcToCm = 3.0856776d21

! Properties of the model galaxy and test case
  real(db), parameter :: r_gal = 10.0 !(kpc)
  real(db), parameter :: y_obs = 2.2e22/kpctocm ! (kpc)
  real(db), parameter :: rho   = 1.0e-23

! Dimensions (in pixels) of column density image
  integer :: npix(2)
! Column density from FITS file
  real, allocatable :: nCol_torus(:,:)

! This is a not a monte-carlo test so set the tolerance just above the error seen 
! when this test was first developed. If it gets any worse then something has gone wrong.
! This is for AMR depth 6:
  real(db), parameter :: tolerance=2.4e21
! This is for AMR depth 7:
!  real(db), parameter :: tolerance=1.1e21

  logical :: fail
  fail=.false.

  write(*,*) "Checking angular image results"

! Check that all the fits files have been generated successfully
  call checkFile("nCol_angImgTest.fits")
  call checkFile("intensity_angImgTest.fits")
  call checkFile("intensity_neg_angImgTest.fits")
  call checkFile("intensity_pos_angImgTest.fits")
  call checkFile("tau_angImgTest.fits")

! Check column densities
  call readNcol
  call checkNcol

  if (fail) then 
     write(*,*) "TORUS: Test failed"
  else
     write(*,*) "TORUS: Test successful"
  endif

contains

! This subroutine simply checks that thisFile exists
  subroutine checkFile(thisFile)

    character(len=*), intent(in) :: thisFile
    logical :: foundFile

    inquire(file=thisFile, exist=foundFile)
    if (foundFile) then 
       write(*,*) "Found "//thisFile
    else
       write(*,*) "Did not find "//thisFile
       fail=.true.
    endif

  end subroutine checkFile

! Read in column density from the Torus FITS file
  subroutine readNcol

! FITS file parameters and variables
    integer :: unit, blocksize, status, npixels, nfound
    integer, parameter :: readwrite=0 ! Open file read only 
    integer, parameter :: group=0
    logical :: anynull
    character(len=30 ) :: errtext
    character(len=80 ) :: errmessage

    character(len=*), parameter :: filename = "nCol_angImgTest.fits"

    write(*,*) "Reading in column density from "//filename

    ! Get a free LUN then open the file
    call ftgiou(unit, status)
    call ftopen(unit, filename,readwrite,blocksize,status)

    ! Find the axis sizes
    call ftgknj(unit,"NAXIS",1,2,npix,nfound,status)
    write(*,'(a,2(2x,i3))') "Axis sizes: ", npix

    ! Read in the pixel values
    allocate (nCol_torus(npix(1),npix(2)))
    npixels = npix(1) * npix(2)
    call ftgpve(unit,group,1,npixels,1e-33,nCol_torus,anynull,status)

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

       fail=.true.

    endif

  end subroutine readNcol

! Check the Torus column density against an analytical expression
  subroutine checkNcol

    integer :: i         ! loop index
    
    real(db) :: theta    ! radians 
    real(db) :: theta_deg! degrees
    real(db) :: a,b,c    ! polynomial co-efficients
    real(db) :: x_outer  ! x co-ordindate of the outer edge of the galaxy
    real(db) :: y_outer  ! y co-ordindate of the outer edge of the galaxy
    real(db) :: path_len ! path length to edge of galaxy
    real(db) :: nCol     ! Column density

    real(db) :: deviation, mean_deviation, nCol_torus_av
    integer :: npoints

    integer :: yPix ! Select this row to compare with
    integer, parameter :: pixPerDegree = 4
    integer :: xPix_min
    integer :: xPix_max

    yPix = npix(2) / 2
    write(*,'(a,i3)') "Comparing with row y= ", yPix

    deviation = 0.0
    npoints = 0

    open(unit=70, status="replace", file="angImgnCol.dat")

    do i = 0,179
       theta_deg = real(i)
       theta = theta_deg * pi/180.0
       a = 1.0 + (tan(theta))**2
       b = 2.0 * tan(theta) * y_obs
       c = y_obs**2 - r_gal**2

! Take the appropriate quadratic root depending on which quadrant we are in
       if ( (theta_deg) < 90.0) then 
          x_outer = ( -1.0*b + (sqrt(b**2 - 4.0*a*c))) / (2.0*a)
       else
          x_outer = ( -1.0*b - (sqrt(b**2 - 4.0*a*c))) / (2.0*a)
       endif

       y_outer = sqrt(r_gal**2 - x_outer**2)

! Path length to the edge of the galaxy. The observer x is always zero.
       path_len = sqrt( (y_outer - y_obs)**2 + x_outer**2 )

! The test case has a constant density 
       nCol = path_len * kpctocm * (rho/mHydrogen)

       xPix_min = theta_deg*pixPerDegree + 1
       xPix_max = theta_deg*pixPerDegree + pixPerDegree

       nCol_torus_av = SUM(nCol_torus(xPix_min:xPix_max,yPix)) / real(pixPerDegree) 

       write(70,*) theta*180.0/pi, nCol, nCol_torus_av
       deviation = deviation + abs( nCol - nCol_torus_av)
       npoints = npoints +1 

    end do

    mean_deviation = deviation/real(npoints)
    write(*,'(a,es10.4)') "Mean deviation= ", mean_deviation
    write(*,'(a,es10.4)') "Tolerance=      ", tolerance
    if (mean_deviation > tolerance) then 
       fail=.true.
    endif

    close(70)

  end subroutine checkNcol

end program check

