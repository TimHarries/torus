program comparespec
  implicit none

  integer :: n1, n2
  character(len=*), parameter :: torus_file = "speca.dat"
  character(len=*), parameter :: ref_file   = "specb.dat"
  double precision, allocatable :: x1(:), x2(:), y1(:), y2(:)
  double precision :: thisVal
  double precision :: loginterp, fac, mean
  integer :: i, nzero, nGoodPoints, status
  integer :: file_line_count
  character(len=100) :: header

! Read Torus file discarding 1 header line
  n1 = file_line_count(torus_file, 1)
  if ( n1 < 1 ) call file_read_error(torus_file) 
  write(*,*) "Reading", n1, " lines from "//torus_file
  allocate (x1(n1), y1(n1))
  open(20,file=torus_file, form="formatted", status="old",iostat=status)
  if ( status /= 0 ) call file_read_error(torus_file)
  read(20,*,iostat=status) header
  if ( status /= 0 ) call file_read_error(torus_file)
  do i = 1, n1
     read(20,*,iostat=status) x1(i), y1(i)
     if ( status /= 0 ) call file_read_error(torus_file)
  enddo
  close(20)

! Read reference SED
  n2 = file_line_count(ref_file, 0)
  if ( n2 < 1 ) call file_read_error(ref_file)
  write(*,*) "Reading", n2, " lines from "//ref_file
  allocate (x2(n2), y2(n2))
  open(20,file=ref_file, form="formatted", status="old", iostat=status)
  if ( status /= 0 ) call file_read_error(ref_file)
  do i = 1, n2
     read(20,*,iostat=status) x2(i), y2(i)
     if ( status /= 0 ) call file_read_error(ref_file)
  enddo
  close(20)

! Check Torus SED against reference SED
  nZero = 0
  mean = 0.
  do i = 1, n1
     if (y1(i) /= 0.) then
        thisval = logInterp(y2, n2, x2, x1(i))
        fac = abs(thisVal-y1(i))/y1(i)
	mean = mean  + fac
     else
        nZero = nZero + 1
     endif
  enddo

  nGoodPoints = n1-nZero
  mean        = mean / dble(nGoodPoints)

! Write out the results
  write(*,*) "Number of comparision points = ", nGoodPoints
  write(*,'(a,f10.2)') "Minimum wavelength for comparison = ", MINVAL(x1(:))
  write(*,'(a,f10.2)') "Maximum wavelength for comparison = ", MAXVAL(x1(:))

  if (nZero == n1) then
     write(*,*) "Spectrum full of zeros..."
     write(*,*) "TORUS: Test failed"
  elseif (mean > 0.1) then
     write(*,*) "Mean relative difference > 10%", mean*100.
     write(*,*) "TORUS: Test failed"
  else
     write(*,*) "TORUS: Test successful"
     write(*,*) "Mean relative difference= ", mean*100, "%"
  end if

  deallocate (x1, x2, y1, y2)

  contains 

    subroutine file_read_error(filename)

      character(len=*), intent(in) :: filename

      write(*,*) "Error reading from "//trim(filename)
      write(*,*) "TORUS: Test failed"
      STOP

    end subroutine file_read_error

end program comparespec

!-------------------------------------------------------------------------------

  PURE SUBROUTINE LOCATE(XX,N,X,J)
    double precision, intent(in)    :: XX(*)
    integer,intent(in)  :: n
    double precision,intent(in)     :: x
    integer,intent(out) :: j
    integer :: jl, ju,jm
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GE.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF

      ! Will force to be between 1 and the array size -1.
      if(x <= xx(1))then
        j=1
      else if(x>=xx(n))then
        j=n-1
      else
        j=jl
      end if

      return

    END SUBROUTINE LOCATE

!-------------------------------------------------------------------------------

  double precision function logInterp(y, ny, x, xi)
    double precision, intent(in)    :: y(*), x(*), xi
    integer, intent(in) :: ny
    integer, save       :: i
    double precision    :: t

    call locate(x, ny, xi, i)
    
    t = (xi - x(i))/(x(i+1)-x(i))

    logInterp = 10.e0**(log10(y(i)) + t * (log10(y(i+1))-log10(y(i))))

  end function logInterp

!-------------------------------------------------------------------------------

! Count the number of real pairs which can be read in from the specified file

  integer function file_line_count(filename, nheader)

    character(len=*), intent(in) :: filename
    integer, intent(in) :: nheader
    real :: dummy1, dummy2
    integer :: status, i
    character(len=100) :: header

    open(unit=30, status="old", file=filename)

! Ignore any header lines
    do i=1, nheader
       read(30,*) header
    end do

    file_line_count = 0
    do
       read(30,*,iostat=status) dummy1, dummy2
       if ( status /= 0 ) exit
       file_line_count = file_line_count + 1 
    end do

    close(30)

  end function file_line_count

!-------------------------------------------------------------------------------
