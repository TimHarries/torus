program comparespec
  implicit none

! Number of lines in Torus SED
  integer, parameter :: n1 = 200
! Number of lines in reference SED
  integer, parameter :: n2 = 61

  double precision :: x1(n1), x2(n2), y1(n1), y2(n2)
  double precision :: thisVal
  double precision :: loginterp, fac, mean
  integer :: i, nzero
  character(len=100) :: header

  open(20,file="speca.dat", form="formatted", status="old",err=666)
  read(20,*) header
  do i = 1, n1
     read(20,*) x1(i), y1(i)
  enddo
  close(20)

  open(20,file="specb.dat", form="formatted", status="old", err=666)
  do i = 1, n2
     read(20,*) x2(i), y2(i)
  enddo
  close(20)


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

  if (nZero == n1) then
     write(*,*) "Spectrum full of zeros..."
     goto 666
  endif

  mean = mean / dble(n1-nZero)
  if (mean > 0.1) then
     write(*,*) "Mean relative difference > 10%",mean*100.
     goto 666
  endif

  write(*,*) "TORUS: Test successful"
  write(*,*) "Mean relative difference= ", mean*100, "%"
  goto 999


666 continue
  write(*,*) "TORUS: Test failed"

999 continue
end program comparespec

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

  double precision function logInterp(y, ny, x, xi)
    double precision, intent(in)    :: y(*), x(*), xi
    integer, intent(in) :: ny
    integer, save       :: i
    double precision    :: t

    call locate(x, ny, xi, i)
    
    t = (xi - x(i))/(x(i+1)-x(i))

    logInterp = 10.e0**(log10(y(i)) + t * (log10(y(i+1))-log10(y(i))))

  end function logInterp
