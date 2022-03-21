! written by tjh


module utils_mod

  use kind_mod
  use vector_mod          ! vector maths
  use constants_mod, only: angstromToCm, cSpeed, kErg, hCgs, ergtoEv, mHydrogen, twoPi, pi, mElectron
  use messages_mod
  use octal_mod

  implicit none

 public

  interface stripSimilarValues
     module procedure stripSimilarValuesDouble
     module procedure stripSimilarValuesSingle
  end interface

  interface locate
     module procedure locate_single
!     module procedure locate_octal
     module procedure locate_double
  end interface

  interface locate_f90
     module procedure locate_single_f90
!     module procedure locate_octal
     module procedure locate_double_f90
  end interface

  interface hunt
     module procedure hunt_single
     module procedure hunt_octal
  end interface

  interface sort
     module procedure sortinteger
     module procedure sortsingle
     module procedure sortdouble
  end interface

  interface logint
     module procedure logint_single
     module procedure logint_dble
  end interface

  interface spline
     module procedure spline_single
     module procedure spline_double
  end interface

  interface splint
     module procedure splint_single
     module procedure splint_double
  end interface

  interface indexx
     module procedure indexx_int
     module procedure indexx_single
     module procedure indexx_double
  end interface

  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE

  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  INTEGER(I4B), PARAMETER, private :: NPAR_ARTH=16,NPAR2_ARTH=8

contains

  function ilog2(val ) result( res )
    integer, intent(in) :: val
    integer             :: res
    integer             :: tmp

    res = -1
    ! Negativ values not allowed
    if ( val < 1 ) return

    tmp = val
    do while (tmp > 0)
      res = res + 1
      tmp = shiftr( tmp, 1 )
    enddo
  end function ilog2

  integer function median(iArray)
    integer :: n, iArray(:)
    integer, allocatable :: iT(:)
    n = size(iArray)
    allocate(it(1:n))
    call sort(n, it)
    median = it(n/2)
  end function median

  real(double) function getGridValue(i, r1, r2, nr, logarithmic)
    integer :: i, nr
    real(double) :: r1, r2
    logical :: logarithmic

    if (.not.logarithmic) then
       if (nr > 1) then
          getGridValue = r1 + (r2-r1) * dble(i-1)/dble(nr-1)
       else
          getGridValue = r1
       endif
    else
       if (nr > 1) then
          getGridValue = log10(r1) + log10(r2/r1) * dble(i-1)/dble(nr-1)
          getGridValue = 10.d0**getGridValue
       else
          getGridValue = r1
       endif
    endif
  end function getGridValue

  ! solve a quadratic equation


  subroutine reverse(x)
    real(double) :: x(:)
    real(double), allocatable :: t(:)
    integer :: n, i

    n = SIZE(x)
    allocate(t(1:n))
    do i = 1 , n
       t(i) = x(n-i+1)
    enddo
    x = t
    deallocate(t)
  end subroutine reverse

  subroutine findMultiFilename(rootFilename, iNumber, actualFilename)
    character(len=*) :: rootFilename, actualFilename
    character(len=80) :: beginString, endString
    character(len=20) :: fortranFormat
    integer :: iNumber
    integer :: i1, i2, nAsterix

    if (index(rootfilename, "*") == 0) then
       actualFilename = rootFilename
       goto 666
    endif

    i1 = index(rootfilename,"*")
    i2 = index(rootfilename, "*", back=.true.)
    beginString = rootFilename(1:i1-1)
    endString = trim(rootFilename(i2+1:))
    nAsterix = i2 - i1 + 1
    write(fortranFormat, '(a,i1,a,i1,a)') "(a,i",nAsterix,".",nAsterix,",a)"
    write(actualFilename, fortranFormat) trim(beginString), iNumber, trim(endString)
666 continue
  end subroutine findMultiFilename

  subroutine solveQuad(a, b, c, x1, x2,ok)
    implicit none
    real,intent(in)     :: a, b, c
    real,intent(out)    :: x1, x2
    real                :: y
    logical,intent(out) :: ok

    ok = .true.

    y = b*b-4.*a*c

    if (y >= 0.) then
       x1 = (-b + sqrt(abs(y)))/(2.*a)
       x2 = (-b - sqrt(abs(y)))/(2.*a)
    else
       write(*,*) "! quad solver failed",y,a,b,c
       ok = .false.
       x1 = -b/(2.*a)
       x2 = -b/(2.*a)
    endif

  end subroutine solveQuad

  subroutine solveQuadDble(a, b, c, x1, x2,ok)
    implicit none
    real(double),intent(in) :: a, b, c
    real(double),intent(out):: x1, x2
    real(double)            :: y
    logical,intent(out)              :: ok

    ok = .true.

    ! special case:
    if (a == 0.d0) then
       if (b/=0.d0) then
          x1 = -c/b; x2 = -c/b
       else
          write(*,*) "! quad solver failed (1):   (a,b,c) = ", a, b, c
          ok = .false.
          x1 = 0.0d0; x2 =0.0d0
       end if
    end if


    y = b*b-4.d0*a*c

    if (y >= 0.d0) then
       x1 = (-b + sqrt(y))/(2.d0*a)
       x2 = (-b - sqrt(y))/(2.d0*a)
    else
       write(*,*) "! quad solver failed (2).:  (y,a,b,c) = ", y,a,b,c
       ok = .false.
       x1 = -b/(2.d0*a)
       x2 = -b/(2.d0*a)
       y = sqrt(y)
       stop
    endif

  end subroutine solveQuadDble

  ! return a blackbody function [B_nu]

  function toPerAngstrom(flux, freq) result(out)
    real(double) :: flux, freq, out, lambda
    lambda = cSpeed / freq
    out = flux * cSpeed / lambda**2 ! cm-1
    out = out * angstromToCm ! per A
  end function toPerAngstrom

  function returnMagnitude(flux, band) result(out)
    real(double) :: flux, out, offset
    character(len=*) :: band
    select case (band)
       case("B")
          offset = 6.4d-9
       case("V")
          offset = 3.75d-9
       case("R")
          offset = 1.8e-9
       case("J")
          offset = 3.133e-10 ! 2mass from gemini website
       case("H")
          offset = 1.111e-10
       case("K")
          offset = 4.288e-11
       case DEFAULT
          call writeWarning("Band "//trim(band)//" not recognised in returnMagntiude")
          offset = 0.d0
    end select
    out = -2.5d0 * log10(flux) + 2.5d0*log10(offset)
  end function returnMagnitude

  function returnFlux(magnitude, band) result(out)
    real(double) :: magnitude, out, offset
    character(len=*) :: band
    select case (band)
       case("B")
          offset = 6.4d-9
       case("V")
          offset = 3.75d-9
       case("R")
          offset = 1.8e-9
       case("J")
          offset = 3.133e-10 ! 2mass from gemini website
       case("H")
          offset = 1.111e-10
       case("K")
          offset = 4.288e-11
       case DEFAULT
          call writeWarning("Band "//trim(band)//" not recognised in returnMagntiude")
          offset = 0.d0
    end select
    out = offset * 10.d0**(-0.4d0*magnitude)
  end function returnFlux

  real elemental function blackBody(temperature, wavelength)

    real,intent(in) :: temperature
    real,intent(in) :: wavelength
    real            :: nu, fac

    nu = real(cSpeed / (wavelength*angstromToCm))
    fac= real(log(2.)+log(hCgs) + 3.*log(nu) - 2.*log(cSpeed))
    if (hCgs*nu/(kerg*temperature) < 33.) then
       blackBody = real(exp(fac)/(exp(hCgs*nu/(kErg*temperature))-1.))
    else
       blackBody = 0.
    endif
  end function blackBody

  ! return random  poissonion deviate

  REAL FUNCTION POIDEV(XM)
    REAL XM, EM, G, OLDM,T,SQ,ALXM, R1
    REAL Y
    REAL PI
    PARAMETER (PI=3.141592654)
    DATA OLDM /-1./
    save
    !$OMP THREADPRIVATE(oldm)
    IF (XM.LT.12.)THEN
       IF (XM.NE.OLDM) THEN
          OLDM=XM
          G=EXP(-XM)
       ENDIF
       EM=-1
       T=1.
2      EM=EM+1.
       call randomNumberGenerator(getReal=r1)
       T=T*r1
       IF (T.GT.G) GO TO 2
    ELSE
       IF (XM.NE.OLDM) THEN
          OLDM=XM
          SQ=SQRT(2.*XM)
          ALXM=ALOG(XM)
          G=XM*ALXM-GAMMLN(XM+1.)
       ENDIF
1      continue
       call randomNumberGenerator(getReal=r1)
       Y=TAN(PI*r1)
       EM=SQ*Y+XM
       IF (EM.LT.0.) GO TO 1
       EM=INT(EM)
       T=0.9*(1.+Y**2)*EXP(EM*ALXM-GAMMLN(EM+1.)-G)
       call randomNumberGenerator(getReal=r1)
       IF (r1.GT.T) GO TO 1
    ENDIF
    POIDEV=EM
  END FUNCTION POIDEV


  REAL FUNCTION GAMMLN(XX)
    REAL XX
    INTEGER  J
    real(double):: COF(6),STP,HALF,ONE,FPF,X,TMP,SER
    DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0, &
         -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
    DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
    X=XX-ONE
    TMP=X+FPF
    TMP=(X+HALF)*LOG(TMP)-TMP
    SER=ONE
    DO  J=1,6
       X=X+ONE
       SER=SER+COF(J)/X
    enddo
    GAMMLN=real(TMP+LOG(STP*SER))
    RETURN
  END FUNCTION GAMMLN

  ! logarithmic interpolation

  REAL FUNCTION LOGINT_single(X,X1,X2,Y1,Y2)
    IMPLICIT NONE
    REAL X,X1,X2,Y1,Y2,ANS
    REAL(double) ::  LX,LX1,LX2,LY1,LY2,GR
    if ( (x.le.0.0).or.(x1.le.0.0).or.(x2.le.0.0)) then
       write(*,*) 'f.up in logint',x,x1,x2,y1,y2
       do;enddo
    endif
    LX=LOG(X)
    LX1=LOG(X1)
    LX2=LOG(X2)
    LY1=LOG(MAX(Y1,1.E-20))
    LY2=LOG(MAX(Y2,1.E-20))
    if (lx1.eq.lx2) then
       write(*,*) 'Error:: Bad x in logint.'
       write(*,*) 'lx1 =', lx1
       write(*,*) 'lx2 =', lx2
       stop
    endif
    GR=(LY2-LY1)/(LX2-LX1)
    ANS=real(LY1+GR*(LX-LX1))
    LOGINT_single=EXP(ANS)
  END FUNCTION LOGINT_SINGLE

  REAL FUNCTION LOGINT_dble(X,X1,X2,Y1,Y2)
    IMPLICIT NONE
    REAL(double) ::  X,X1,X2,Y1,Y2,ANS
    REAL(double) ::  LX,LX1,LX2,LY1,LY2,GR
    if ( (x.le.0.0d0).or.(x1.le.0.0d0).or.(x2.le.0.0d0)) then
       write(*,*) 'f.up in logint_dble',x,x1,x2,y1,y2
              do;enddo

    endif

    LX=LOG(X)
    LX1=LOG(X1)
    LX2=LOG(X2)
    LY1=LOG(MAX(Y1,1.d-30))
    LY2=LOG(MAX(Y2,1.d-30))
    if (lx1.eq.lx2) then
       write(*,*) 'Error:: Bad x in logint_dble.'
       write(*,*) 'lx1 =', lx1
       write(*,*) 'lx2 =', lx2
       stop
    endif
    GR=(LY2-LY1)/(LX2-LX1)
    ANS=LY1+GR*(LX-LX1)
    LOGINT_dble=real(EXP(ANS))
  END FUNCTION LOGINT_DBLE

  SUBROUTINE dswap(a,b)
    REAL(DP), INTENT(INOUT) :: a,b
    REAL(DP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE dswap

  SUBROUTINE dmasked_swap(a,b,mask)
    REAL(DP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(DP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE dmasked_swap


  ! sort an array

  SUBROUTINE dquicksort(arr)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
    INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
    REAL(DP) :: a
    INTEGER(I4B) :: n,k,i,j,jstack,l,r
    INTEGER(I4B), DIMENSION(NSTACK) :: istack
    n=size(arr)
    jstack=0
    l=1
    r=n
    do
       if (r-l < NN) then
          do j=l+1,r
             a=arr(j)
             do i=j-1,l,-1
                if (arr(i) <= a) exit
                arr(i+1)=arr(i)
             end do
             arr(i+1)=a
          end do
          if (jstack == 0) RETURN
          r=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2
       else
          k=(l+r)/2
          call dswap(arr(k),arr(l+1))
          call dmasked_swap(arr(l),arr(r),arr(l)>arr(r))
          call dmasked_swap(arr(l+1),arr(r),arr(l+1)>arr(r))
          call dmasked_swap(arr(l),arr(l+1),arr(l)>arr(l+1))
          i=l+1
          j=r
          a=arr(l+1)
          do
             do
                i=i+1
                if (arr(i) >= a) exit
             end do
             do
                j=j-1
                if (arr(j) <= a) exit
             end do
             if (j < i) exit
             call dswap(arr(i),arr(j))
          end do
          arr(l+1)=arr(j)
          arr(j)=a
          jstack=jstack+2
          if (jstack > NSTACK) call writeFatal('sort: NSTACK too small')
          if (r-i+1 >= j-l) then
             istack(jstack)=r
             istack(jstack-1)=i
             r=j-1
          else
             istack(jstack)=j-1
             istack(jstack-1)=l
             l=i
          end if
       end if
    end do
  END SUBROUTINE dquicksort

  SUBROUTINE SORTsingle(N,RA)
    INTEGER N, L, IR, I, J
    REAL RRA
    REAL RA(:)
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
    ELSE
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).LT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.LT.RA(J))THEN
          RA(I)=RA(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    GO TO 10
  END subroutine sortsingle

  SUBROUTINE SORTinteger(N,RA)
    INTEGER N, L, IR, I, J
    integer RRA
    integer RA(:)
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
    ELSE
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).LT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.LT.RA(J))THEN
          RA(I)=RA(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    GO TO 10
  END subroutine sortinteger

  SUBROUTINE SORTdouble(N,RA)
    INTEGER N, L, IR, I, J
    real(double) RRA
    real(double) RA(:)
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
    ELSE
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).LT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.LT.RA(J))THEN
          RA(I)=RA(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    GO TO 10
  END subroutine sortdouble

  SUBROUTINE SORTdouble2(N,RA,rb)
    INTEGER N, L, IR, I, J
    real(double) RRA, rrb
    real(double) RA(N),rb(n)
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
       rrb = rb(L)
    ELSE
       RRA=RA(IR)
       rrb=rb(ir)
       RA(IR)=RA(1)
       rb(ir) = rb(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          rb(1) = rrb
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).LT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.LT.RA(J))THEN
          RA(I)=RA(J)
          rb(i)=rb(j)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    rb(i) = rrb
    GO TO 10
  END subroutine sortdouble2

  SUBROUTINE SORTdouble2index(RA, rb)
    INTEGER :: N, L, IR, I, J
    real(double) :: RRA
    integer :: rrb,rb(:)
    real(double) :: RA(:)

    N = size(RA)

    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
       rrb = rb(L)
    ELSE
       RRA=RA(IR)
       rrb=rb(ir)
       RA(IR)=RA(1)
       rb(ir) = rb(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          rb(1) = rrb
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).LT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.LT.RA(J))THEN
          RA(I)=RA(J)
          rb(i)=rb(j)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    rb(i) = rrb
    GO TO 10
  END subroutine sortdouble2index



  !
  !
  ! Taken from Numerical recipes
  ! See the comments in the book.
  SUBROUTINE polint(xa,ya,n,x,y,dy)
    INTEGER n,NMAX
    REAL dy,x,y,xa(:),ya(:)
    PARAMETER (NMAX=10)
    INTEGER i,m,ns
    REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=abs(x-xa(1))
    do i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          ns=i
          dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
    end do
    y=ya(ns)
    ns=ns-1
    do  m=1,n-1
       do  i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             print *, 'Error:: failure in polint.'
             stop
          end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       end do
       if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    end do
    return
  end SUBROUTINE polint


  SUBROUTINE polint2(xa,ya,n,x,y,dy)
    INTEGER n,NMAX
    REAL :: x,xa(:)
    real(double) :: y, ya(:), dy
    PARAMETER (NMAX=10)
    INTEGER i,m,ns
    REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=abs(x-xa(1))
    do i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          ns=i
          dif=dift
       endif
       c(i)=real(ya(i))
       d(i)=real(ya(i))
    end do
    y=ya(ns)
    ns=ns-1
    do  m=1,n-1
       do  i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             print *, 'Error:: failure in polint2.'
             stop
          end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       end do
       if (2*ns.lt.n-m)then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       endif
       y=y+dy
    end do
    return
  end SUBROUTINE polint2

  ! spline interpolation

  SUBROUTINE SPLINE_SINGLE(X,Y,N,YP1,YPN,Y2)
    INTEGER NMAX,N, I, K
    REAL SIG, P, QN, UN
    PARAMETER (NMAX=10000)
    REAL X(:),Y(:),Y2(:),U(NMAX),YP1,YPN
    IF (YP1.GT..99E30) THEN
       Y2(1)=0.
       U(1)=0.
    ELSE
       Y2(1)=-0.5
       U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
    ENDIF
    DO I=2,N-1
       SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
       P=SIG*Y2(I-1)+2.
       Y2(I)=(SIG-1.)/P
       U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
            /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
    ENDDO
    IF (YPN.GT..99E30) THEN
       QN=0.
       UN=0.
    ELSE
       QN=0.5
       UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
    ENDIF
    Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
    DO  K=N-1,1,-1
       Y2(K)=Y2(K)*Y2(K+1)+U(K)
    ENDDO
  END SUBROUTINE SPLINE_SINGLE


  SUBROUTINE SPLINT_SINGLE(XA,YA,Y2A,N,X,Y)
    INTEGER N, KLO, KHI, K
    REAL H, A, B
    REAL XA(:),YA(:),Y2A(:),X,Y
    KLO=1
    KHI=N
1   IF (KHI-KLO.GT.1) THEN
       K=(KHI+KLO)/2
       IF(XA(K).GT.X)THEN
          KHI=K
       ELSE
          KLO=K
       ENDIF
       GOTO 1
    ENDIF
    H=XA(KHI)-XA(KLO)
    IF (H.EQ.0.) then
       write(*,*) 'Bad XA input in splint_single.'
       stop
    endif
    A=(XA(KHI)-X)/H
    B=(X-XA(KLO))/H
    Y=A*YA(KLO)+B*YA(KHI)+ &
         ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
  END SUBROUTINE SPLINT_SINGLE

  ! From Numerical Recipes in F77
  ! Modified by chris
  subroutine spline_double(x,y,n,yp1,ypn,y2)
    integer :: n
    real(double) :: yp1, ypn, x(:), y(:), y2(:)
    integer, parameter :: nmax=5000
    integer :: i, k
    real(double) :: p, qn, sig, un, u(nmax)

    if (yp1 .gt. .99d30) then
        y2(1) = 0.d0
        u(1) = 0.d0
    else
        y2(1) = -0.5d0
        u(1) = (3.d0/(x(2)-x(1))) * ((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
!    write(*,*) "y2(1)",y2(1)
!    write(*,*) "u(1)",u(1)

    do i=2, n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2.d0
        y2(i) = (sig-1.d0)/p
        u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
                / (x(i+1)-x(i-1)) - sig*u(i-1))/p

!        write(*,*) "LOOP", "so...", sig,p,y2(i),u(i)
!        write(*,*) "LOOP", "so...", 6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)),(x(i+1)-x(i-1)) - sig*u(i-1))/p

    Enddo

    if (ypn .gt. .99d30) then
        qn = 0.d0
        un = 0.d0
    else
        qn = 0.5d0
        un = (3.d0/(x(n)-x(n-1))) * (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif

    y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.d0)
!    write(*,*) "y2(n)",y2(n),un

    do k=n-1, 1, -1
        y2(k) = y2(k)*y2(k+1)+u(k)
!        write(*,*) u(k), "u(K)"
!        write(*,*) y2(k)
    enddo

    return
  end subroutine spline_double

  ! From Numerical Recipes in F77
  ! Modified by chris
  subroutine splint_double(xa,ya,y2a,n,x,y)
    integer :: n
    real(double) :: x, y, xa(:), y2a(:), ya(:)
    integer :: k, khi, klo
    real(double) :: a, b, h

    klo = 1
    khi = n

1   if (khi-klo .gt. 1) then
        k = (khi+klo)/2
        if (xa(k) .gt. x) then
            khi = k
        else
            klo = k
        endif
    goto 1
    endif

    h = xa(khi) - xa(klo)
    if (h .eq. 0.d0) call writeFatal("Bad xa input in splint_double.")
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi)) &
        * (h**2) / 6.d0

    return
  end subroutine splint_double


  ! locate in a grid via bisection but starting at jlo

  pure SUBROUTINE HUNT_single(XX,N,X,JLO)
    REAL, INTENT(IN) :: XX(:)
    INTEGER, INTENT(IN) :: N
    REAL, INTENT(IN)    :: X
    INTEGER, INTENT(INOUT) :: JLO
    INTEGER :: JHI, INC, JM
    LOGICAL ASCND
    IF (X .LE. XX(1)) THEN
       JLO = 1
       RETURN
    ENDIF
    ASCND=XX(N).GT.XX(1)
    IF(JLO.LE.0.OR.JLO.GT.N)THEN
       JLO=0
       JHI=N+1
       GO TO 3
    ENDIF
    INC=1
    IF(X.GE.XX(JLO).EQV.ASCND)THEN
1      JHI=JLO+INC
       IF(JHI.GT.N)THEN
          JHI=N+1
       ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
       ENDIF
    ELSE
       JHI=JLO
2      JLO=JHI-INC
       IF(JLO.LT.1)THEN
          JLO=0
       ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
       ENDIF
    ENDIF
3   IF(JHI-JLO.EQ.1)RETURN
    JM=(JHI+JLO)/2
    IF(X.GT.XX(JM).EQV.ASCND)THEN
       JLO=JM
    ELSE
       JHI=JM
    ENDIF
    GO TO 3
  END SUBROUTINE HUNT_single


  pure SUBROUTINE HUNT_octal(XX,N,X,JLO)
    real(oct), INTENT(IN) :: XX(:)
    INTEGER, INTENT(IN) :: N
    real(oct), INTENT(IN)    :: X
    INTEGER, INTENT(INOUT) :: JLO
    INTEGER :: JHI, INC, JM
    LOGICAL ASCND
    IF (X .LE. XX(1)) THEN
       JLO = 1
       RETURN
    ENDIF
    ASCND=XX(N).GT.XX(1)
    IF(JLO.LE.0.OR.JLO.GT.N)THEN
       JLO=0
       JHI=N+1
       GO TO 3
    ENDIF
    INC=1
    IF(X.GE.XX(JLO).EQV.ASCND)THEN
1      JHI=JLO+INC
       IF(JHI.GT.N)THEN
          JHI=N+1
       ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
       ENDIF
    ELSE
       JHI=JLO
2      JLO=JHI-INC
       IF(JLO.LT.1)THEN
          JLO=0
       ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
       ENDIF
    ENDIF
3   IF(JHI-JLO.EQ.1)RETURN
    JM=(JHI+JLO)/2
    IF(X.GT.XX(JM).EQV.ASCND)THEN
       JLO=JM
    ELSE
       JHI=JM
    ENDIF
    GO TO 3
  end SUBROUTINE HUNT_octal


  real function logInterp(y, ny, x, xi)
    real, intent(in)    :: y(:), x(:), xi
    integer, intent(in) :: ny
    integer, save       :: i
    real                :: t
    !$OMP THREADPRIVATE (i)
    call locate(x, ny, xi, i)

    t = (xi - x(i))/(x(i+1)-x(i))

    logInterp = 10.e0**(log10(y(i)) + t * (log10(y(i+1))-log10(y(i))))

  end function logInterp

  real(double) function logInterp_dble(y, n, x, xi)
    integer, intent(in) :: n
    real(double), intent(in)    :: y(:), x(:), xi
    integer, save       :: i
    real(double)                :: t
    !$OMP THREADPRIVATE (i)

    call locate(x, n, xi, i)

    t = (xi - x(i))/(x(i+1)-x(i))

    logInterp_dble = 10.d0**(log10(y(i)) + t * (log10(y(i+1))-log10(y(i))))

  end function logInterp_dble

  PURE SUBROUTINE LOCATE_single(XX,N,X,J)
    real, intent(in)    :: XX(:)
    integer,intent(in)  :: n
    real,intent(in)     :: x
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

    END SUBROUTINE LOCATE_single


    PURE SUBROUTINE LOCATE_double(XX,N,X,J)
    real(double), intent(in) :: XX(:)
    integer, intent(in)              :: n
    real(double), intent(in) :: x
    integer,intent(out)              :: j
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

    END SUBROUTINE LOCATE_DOUBLE

    PURE SUBROUTINE LOCATE_octal(XX,N,X,J)
    real(oct), intent(in) :: XX(:)
    integer, intent(in)              :: n
    real(oct), intent(in) :: x
    integer,intent(out)              :: j
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

    END SUBROUTINE LOCATE_OCTAL

    pure subroutine locate_single_f90(xx,x,locate)
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: xx
      REAL, INTENT(IN) :: x
      INTEGER,intent(OUT) :: locate
      INTEGER :: n,jl,jm,ju
      LOGICAL :: ascnd

      n=size(xx)
      ascnd = (xx(n) >= xx(1))
      jl=0
      ju=n+1
      do
         if (ju-jl <= 1) exit
         jm=(ju+jl)/2
         if (ascnd .eqv. (x >= xx(jm))) then
            jl=jm
         else
            ju=jm
         end if
      end do
      if (x == xx(1)) then
         locate=1
      else if (x == xx(n)) then
         locate=n-1
      else
         locate=jl
      end if
    END subroutine locate_single_f90

    pure subroutine locate_double_f90(xx,x,locate)
      IMPLICIT NONE
      REAL(double), DIMENSION(:), INTENT(IN) :: xx
      REAL(double), INTENT(IN) :: x
      INTEGER, intent(OUT) :: locate
      INTEGER :: n,jl,jm,ju
      LOGICAL :: ascnd
      n=size(xx)
      ascnd = (xx(n) >= xx(1))
      jl=0
      ju=n+1
      do
         if (ju-jl <= 1) exit
         jm=(ju+jl)/2
         if (ascnd .eqv. (x >= xx(jm))) then
            jl=jm
         else
            ju=jm
         end if
      end do
      if (x == xx(1)) then
         locate=1
      else if (x == xx(n)) then
         locate=n-1
      else
         locate=jl
      end if
    END subroutine locate_double_f90

    real function gasdev()
      implicit none
      integer :: iset
      real :: r1, r2, gset
      real :: v1, v2, r, fac
      data iset/0/
      save iset,gset

!$OMP THREADPRIVATE (iset, gset)

      if (iset.eq.0) then
1        continue
         call randomNumberGenerator(getReal=r1)
         call randomNumberGenerator(getReal=r2)
         v1=2.*r1-1.
         v2=2.*r2-1.
         r=v1**2+v2**2
         if(r.ge.1.)go to 1
         fac=sqrt(-2.*log(r)/r)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
    end function gasdev

    pure real function gauss(sigma, dx)
      real, intent(in) :: sigma, dx
      real :: fac
      fac = real(1./(sigma*sqrt(twoPi)))
      gauss = fac * exp(-(dx**2)/(2.*sigma**2))
    end function gauss

    pure real(double) function unNormalizedGauss(sigma, dx)
      real(double), intent(in) :: sigma, dx
      unNormalizedGauss = exp(-(dx**2)/(2.*sigma**2))
    end function unNormalizedGauss


    type(VECTOR) function maxwellianVelocity(mass, temperature)


! Calculates a random maxwellian velocity using the
! rejection method.

      real(double), intent(in) :: mass
      real, intent(in) ::  temperature
      real(double) :: x, y, t, u, vel
      logical :: ok

      ok = .false.
      do while(.not.ok)

         call randomNumberGenerator(getDouble=x)
         x = 3.d0 * x
         call randomNumberGenerator(getDouble=y)

         t = 4.d0/sqrt(pi) * x**2 * exp(-x**2)

         if (y < t) then
            ok = .true.
            u = x
         endif

      end do

      vel = u * sqrt((2*kErg*temperature)/mass)

      maxwellianVelocity = vel * randomUnitVector()
    end function maxwellianVelocity

    real function fvmaxwellian(v, mass, temperature)

! calculated value of f(v) for a maxwellian speed distribution

      real, intent(in) :: v, mass, temperature
      real :: fac, fac2

      fac = real((4./sqrt(pi))*(mass/(2.*kErg*temperature))**1.5)
      fac2 = real(v**2 * exp(-mass*v**2 / (2.*Kerg*temperature)))
      fvmaxwellian = fac * fac2
    end function fvmaxwellian


    ! Computes a random Lorentzian frequency using rejection method
    ! -- The output in [Hz]
    !
    !  The form of the profile is assumed to be
    !                   (Gamma/4Pi)^2
    !  phi(nu) = ------------------------------
    !              (nu-nu0)^2 +  (Gamma/4Pi)^2
    !
    real(double) function random_Lorentzian_frequency(nu0, Gamma)
      real(double), intent(in) :: nu0   !  line center freq.  [Hz]
      real(double), intent(in) :: Gamma !  Damping Constant  [Hz]
      !
      real(double) :: x, y, t, a, a2
      logical :: ok
      a = Gamma/4.0d0/Pi
      a2 = a*a

      ok = .false.

      if (a<=0.0d0) then  ! special case. No shift in line frequency
         x=0.0d0
         ok = .true.
      end if

      do while(.not.ok)
         ! trial value of (nu-nu0)
         call randomNumberGenerator(getDouble=x)
         ! random number between -6a to 6a
         x =6.0d0*a*(-1.0d0 + 2.0d0*x)
         ! trial value of phi(nu)
         call randomNumberGenerator(getDouble=y)

         t = a2/(x*x +a2)  ! This should be always less than 1.

         if (y < t) then
            ok = .true.
         endif

      end do

      ! x = nu - nu0, so
      random_Lorentzian_frequency = x + nu0  ! [Hz]

    end function random_Lorentzian_frequency


    integer function findIlambda(lambda, xArray, nLambda, ok)
      implicit none
      integer :: nlambda, i
      real :: lambda
      real :: xArray(:)
      logical, intent(out) :: ok

      ok = .true.

      if (lambda < (xArray(1))) then
         findiLambda = 1
         ok = .false.
         goto 666
      endif
      if (lambda > (xArray(nLambda))) then
         findiLambda = nLambda
         ok = .false.
         goto 666
      endif

      call locate(xArray, nLambda, lambda, i)
      if (i /= nLambda) then
         if (lambda > 0.5*(xArray(i)+xArray(i+1))) then
            i= i+1
         endif
      endif
      findilambda = min(i,nLambda)
666   continue
    end function findilambda

    integer function findIlambdaDouble(lambda, xArray, nLambda, ok)
      implicit none
      integer :: nlambda, i
      real(double) :: lambda
      real(double) :: xArray(:)
      logical, intent(out) :: ok

      ok = .true.

      if (lambda < (xArray(1))) then
         findiLambdaDouble = 1
         ok = .false.
         goto 666
      endif
      if (lambda > (xArray(nLambda))) then
         findiLambdaDouble = nLambda
         ok = .false.
         goto 666
      endif

      call locate(xArray, nLambda, lambda, i)
      if (i /= nLambda) then
         if (lambda > 0.5*(xArray(i)+xArray(i+1))) then
            i= i+1
         endif
      endif
      findilambdaDouble = min(i,nLambda)
666   continue
    end function findilambdaDouble

    real function interpLinearSingle(xArray, yArray, n, x)
      real :: xarray(:), yArray(:)
      real :: x, t
      integer :: n, i
      logical :: ok
      ok = .false.

      i = findIlambda(x, xArray, n, ok)
      if (i == 1) then
         t = (x - xArray(i))/(xArray(i+1)-xArray(i))
         interpLinearSingle = yArray(i) + t * (yArray(i+1)-yArray(i))
      else if (i == n) then
         t = (x - xArray(i-1))/(xArray(i)-xArray(i-1))
         interpLinearSingle = yArray(i-1) + t * (yArray(i)-yArray(i-1))
      else
         if (x < xArray(i)) then
            t = (x - xArray(i-1))/(xArray(i)-xArray(i-1))
            Interplinearsingle = yArray(i-1) + t * (yArray(i)-yArray(i-1))
         else
            t = (x - xArray(i))/(xArray(i+1)-xArray(i))
            interpLinearSingle = yArray(i) + t * (yArray(i+1)-yArray(i))
         endif
      endif
    end function interpLinearSingle

    real(double) function interpLinearDouble(xArray, yArray, n, x)
      real(double) :: xarray(:), yArray(:)
      real(double) :: x, t
      integer :: n, i
      logical :: ok
      ok = .false.
      i = findIlambda(real(x), real(xArray), n, ok)
      if (i == 1) then
         t = (x - xArray(i))/(xArray(i+1)-xArray(i))
         interpLinearDouble = yArray(i) + t * (yArray(i+1)-yArray(i))
      else if (i == n) then
         t = (x - xArray(i-1))/(xArray(i)-xArray(i-1))
         interpLinearDouble = yArray(i-1) + t * (yArray(i)-yArray(i-1))
      else
         if (x < xArray(i)) then
            t = (x - xArray(i-1))/(xArray(i)-xArray(i-1))
            interpLinearDouble = yArray(i-1) + t * (yArray(i)-yArray(i-1))
         else
            t = (x - xArray(i))/(xArray(i+1)-xArray(i))
            interpLinearDouble = yArray(i) + t * (yArray(i+1)-yArray(i))
         endif
      endif
    end function interpLinearDouble

    real(double) function interpLogLinearDouble(xArray, yArray, n, x)
      real(double) :: xarray(:), yArray(:)
      real(double) :: x, t
      integer :: n, i
      logical :: ok
      ok = .false.
      i = findIlambda(real(x), real(xArray), n, ok)
      if (i == 1) then
         t = (x - xArray(i))/(xArray(i+1)-xArray(i))
         interpLogLinearDouble = 10.d0**(log10(yArray(i)) + t * (log10(yArray(i+1))-log10(yArray(i))))
      else if (i == n) then
         t = (x - xArray(i-1))/(xArray(i)-xArray(i-1))
         interpLogLinearDouble = 10.d0**(log10(yArray(i-1)) + t * (log10(yArray(i))-log10(yArray(i-1))))
      else
         if (x < xArray(i)) then
            t = (x - xArray(i-1))/(xArray(i)-xArray(i-1))
            interpLogLinearDouble = 10.d0**(log10(yArray(i-1)) + t * (log10(yArray(i))-log10(yArray(i-1))))
         else
            t = (x - xArray(i))/(xArray(i+1)-xArray(i))
            interpLogLinearDouble = 10.d0**(log10(yArray(i)) + t * (log10(yArray(i+1))-log10(yArray(i))))
         endif
      endif
    end function interpLogLinearDouble



    function convertToJanskies(flux, wavelength) result (newflux)
      real(double) :: flux   ! in erg/s/cm^2/A
      real(double) :: wavelength ! in A
      real(double) :: newFlux, nu

      nu  = cspeed / (wavelength * angstromtocm)
      newFlux = flux * (cSpeed * 1.e8)/ nu**2 ! from /A to /Hz

      newFlux = newFlux * 1.e23 ! to janskies
    end function convertToJanskies



    !
    ! -----Computes UNNORMALIZED voigtn function: (normalization = sqrt(pi))
    !
    FUNCTION VOIGTN(AA,VV)
      IMPLICIT NONE
      ! Modefied: 05-Nov-2004 :: Now treats a = 0 as a special case.   (R. Kurosawa)
      ! Imported from JDH's CMFGEN code: 27-Oct-2004    (R. Kurosawa)
      ! Cleaned : 20-Nov-2000 (JDH)
      !
      REAL(double) :: VOIGTN,AA,VV
      !
      REAL(double) ::  H(25)
      !
      REAL(double) :: V,V2,V4,V6
      REAL(double) ::  A,A2,A4,A6
      REAL(double) ::  Z,Z2
      REAL(double) ::  X,W
      REAL(double) ::  C1,C2
      DATA C1,C2/1.128379167095512D0  ,5.64189583547756D-1/
      SAVE C1,C2
      !$OMP THREADPRIVATE(c1,c2)
      !
      INTEGER(bigint) :: I, J
      !
      !REAL*8 DAWSON
      !EXTERNAL DAWSON
      !
      V=ABS(VV)
      A=AA
      V2=V*V
      A2=A*A
      Z=A2+V2
     IF(A .EQ. 0.0D0  ) THEN !!uncommented out - tjgw201

        ! ---- Normal Doppler brodening.
        IF (V2 < 1.D-04) THEN
           VOIGTN = 1.0D0 + V2*(1.0D0 + 0.5D0*V2*(1.0D0 + V2/6.0D0))
        ELSE
           VOIGTN = EXP(-V2)
        END IF
        RETURN
     END IF
      IF(A .LE. 0.5D0  )GOTO 20
      IF(Z .LT. 10.D0  )GOTO 50
      !-----ASYMPTOTIC EXPANSION FOR LARGE MODULUS
10    Z2=Z*Z
      V4=V2*V2
      V6=V4*V2
      A4=A2*A2
      A6=A4*A2
      VOIGTN=C2*A* (1.D0  +  ((1.875D0  *(7.D0  *V6-35.D0  *A2*V4+21.D0*A4*V2-A6) &
           /Z2+ 0.75D0  *(5.D0  *V4-10.D0  *A2*V2+A4))/Z2+1.5D0  *V2-0.5D0*A2)/Z2)/Z
      RETURN
      !-----HARRIS EXPANSION
20    IF(V .GT. 5.D0  ) GOTO 10
      W=DAWSON(V)
      H(1)= EXP(-V2)
      H(2)=-C1*(1.0D0  -2.0D0  *V*W)
      H(3)=(1.0D0  -2.0D0  *V2)*H(1)
      H(4)=-C1*(2.D0  *(1.D0  -V2)/3.D0  -2.D0  *V*W*(1.D0  -2.D0  *V2/3.D0  ))
      !-----HIGHER TERMS BY RECURSION
      DO I=5,11
         X=I-1
         H(I)= (2.0D0 *(2.0D0*X -3.0D0 -2.0D0*V2) *H(I-2) -4.0D0*H(I-4))/(X*(X-1.0D0  ))
      END DO
      VOIGTN=H(11)
      DO I=1,10
         J=11-I
         VOIGTN=H(J)+A*VOIGTN
      END DO
      RETURN
      !-----GRONWALL EXPANSION
50    X=1.0D0  /(1.0D0  +3.275911D-1*A)
      H(1)=((((1.061405429D0  *X-1.453152027D0  )*X+1.421413741D0  )*X &
           -2.84496736D-1)*X+2.54829592D-1)*X
      DO I=2,25
         X=I-1
         H(I)=2.0D0  *A*(C2-A*H(I-1))/(2.0D0  *X-1.0D0  )
      END DO
      VOIGTN=0.0D0
      DO I=1,24
         J=26-I
         X=J-1
         VOIGTN=(VOIGTN+H(J))*V2/X
      END DO
      VOIGTN= EXP(-V2)*(VOIGTN+H(1))
      RETURN
    END FUNCTION VOIGTN



    !
    !-----DAWSON*S INTEGRAL USING ANL ALGORITHM, MATH COMP,1970,171
    !
    FUNCTION DAWSON(XX)
      IMPLICIT NONE
      ! Imported from JDH's CMFGEN code: 27-Oct-2004    (R. Kurosawa)
      !
      REAL(double) ::  DAWSON,XX
      REAL(double) ::  X,U
      REAL(double) ::  UP,DOWN
      !
      X=XX
      U=X*X
      IF(X .LT. 5.0D0  ) GOTO 10
      !-----X GREATER THAN 5
      DAWSON= ((5.0000000167450D-1+7.4999919056701D-1/      &
           (-2.5001711668562D0  +U-2.4878765880441D0  /     &
           (-4.6731202214124D0  +U-4.1254406560831D0  /     &
           (-1.1195216423662D1+U))))/U+1.0D0  )/(2.0D0  *X)
      RETURN
      !-----X ON (3.5D0,5.0D0)
10    IF(X .LT. 3.5D0  ) GOTO 20
      DAWSON=(5.00001538408193D-1 +2.49811162845499D-1/     &
           (-1.53672069271915D0  +U-6.53419359860764D-1/    &
           (-1.77068693717670D1 +U+2.04866410976332D2/      &
           (7.49584016278357D0   +U-2.298758419286D0  /     &
           (4.02187490205698D1+U+2.53388006963558D3/        &
           (-5.9391591850032E1+U))))))/X
      RETURN
      !-----X ON (2.5,3.5)
20    IF(X .LT. 2.5D0  ) GOTO 30
      DAWSON=(5.0140106611704D-1+1.8897553014354D-1/      &
           (-7.4499050579364D0  +U+7.0204980729194D1/     &
           (7.5077816490106D0  +U+4.1821806337830D1/      &
           (-2.6629001073842D1 +U+3.7343084728334D1/      &
           (3.0984087863402D1+U+1.2599323546764D3/        &
           (-4.0847391212716D1+U))))))/X
      RETURN
      !-----X LESS THAN 2.5
   30 UP=(((((U*2.0846835103886D-2 -8.5410681195954D-1)*U  &
            +5.4616122556699D1)*U-4.3501160207595D2)*U     &
            +9.6696398191665D3)*U-2.9179464300780D4)*U+2.3156975201341D5
      DOWN=((((( U+2.9391995612556D1)*U +4.668490654511D2  &
           )*U+4.7447098440662D3)*U+3.1384620138163D4)*U   &
           +1.2520037031851D5)*U+2.3156975201425D5
      DAWSON=X*(UP/DOWN)
      RETURN
    END FUNCTION DAWSON



    subroutine voigt ( n, x, y, k )

      !     to calculate the faddeeva function with relative error less than 10^(-4).

      ! arguments
      integer n                                                         ! in   number of points
      real    x(0:n-1)                                                  ! in   input x array
      real    y                                                         ! in   input y value >=0.0
      real    k(0:n-1)                                                  ! out  real (voigt) array

      ! constants
      real        rrtpi                                                 ! 1/sqrt(pi)
      parameter ( rrtpi = 0.56418958 )
      real        y0,       y0py0,         y0q                          ! for cpf12 algorithm
      parameter ( y0 = 1.5, y0py0 = y0+y0, y0q = y0*y0  )
      real  c(0:5), s(0:5), t(0:5)
      save  c,      s,      t
      !     save preserves values of c, s and t (static) arrays between procedure calls
      data c / 1.0117281,     -0.75197147,        0.012557727, &
               0.010022008,   -0.00024206814,     0.00000050084806 /
      data s / 1.393237,       0.23115241,       -0.15535147, &
               0.0062183662,   0.000091908299,   -0.00000062752596 /
      data t / 0.31424038,     0.94778839,        1.5976826,&
               2.2795071,      3.0206370,         3.8897249 /

      ! local variables
      integer i, j                                                      ! loop variables
      integer rg1, rg2, rg3                                             ! y polynomial flags
      real abx, xq, yq, yrrtpi                                          ! |x|, x^2, y^2, y/sqrt(pi)
      real xlim0, xlim1, xlim2, xlim3, xlim4                            ! |x| on region boundaries
      real a0, d0, d2, e0, e2, e4, h0, h2, h4, h6                       ! w4 temporary variables
      real p0, p2, p4, p6, p8, z0, z2, z4, z6, z8
      real xp(0:5), xm(0:5), yp(0:5), ym(0:5)                           ! cpf12 temporary values
      real mq(0:5), pq(0:5), mf(0:5), pf(0:5)
      real d, yf, ypy0, ypy0q

      !$OMP THREADPRIVATE (c, s, t)
      !**** start of executable code *****************************************

      yq  = y*y                                                         ! y^2
      yrrtpi = y*rrtpi                                                  ! y/sqrt(pi)

      if ( y .ge. 70.55 ) then                                          ! all points
         do i = 0, n-1                                                   ! in region 0
            xq   = x(i)*x(i)
            k(i) = yrrtpi / (xq + yq)
         enddo
         return
      endif

      rg1 = 1                                                           ! set flags
      rg2 = 1
      rg3 = 1

      xlim0 = sqrt ( 15100.0 + y*(40.0 - y*3.6) )                       ! y<70.55
      if ( y .ge. 8.425 ) then
         xlim1 = 0.0
      else
         xlim1 = sqrt ( 164.0 - y*(4.3 + y*1.8) )
      endif
      xlim2 = 6.8 - y
      xlim3 = 2.4*y
      xlim4 = 18.1*y + 1.65
      if ( y .le. 0.000001 ) then                                       ! when y<10^-6
         xlim1 = xlim0                                                    ! avoid w4 algorithm
         xlim2 = xlim0
      endif

      do i = 0, n-1                                                     ! loop over all points
         abx = abs ( x(i) )                                               ! |x|
         xq  = abx*abx                                                    ! x^2
         if     ( abx .ge. xlim0 ) then                                   ! region 0 algorithm
            k(i) = yrrtpi / (xq + yq)

         elseif ( abx .ge. xlim1 ) then                                   ! humlicek w4 region 1
            if ( rg1 .ne. 0 ) then                                          ! first point in region 1
               rg1 = 0
               a0 = yq + 0.5                                                  ! region 1 y-dependents
               d0 = a0*a0
               d2 = yq + yq - 1.0
            endif
            d = rrtpi / (d0 + xq*(d2 + xq))
            k(i) = d*y   *(a0 + xq)

         elseif ( abx .gt. xlim2 ) then                                   ! humlicek w4 region 2
            if ( rg2 .ne. 0 ) then                                          ! first point in region 2
               rg2 = 0
               h0 =  0.5625 + yq*(4.5 + yq*(10.5 + yq*(6.0 + yq)))            ! region 2 y-dependents
               h2 = -4.5    + yq*(9.0 + yq*( 6.0 + yq* 4.0))
               h4 = 10.5    - yq*(6.0 - yq*  6.0)
               h6 = -6.0    + yq* 4.0
               e0 =  1.875  + yq*(8.25 + yq*(5.5 + yq))
               e2 =  5.25   + yq*(1.0  + yq* 3.0)
               e4 =  0.75*h6
            endif
            d = rrtpi / (h0 + xq*(h2 + xq*(h4 + xq*(h6 + xq))))
            k(i) = d*y   *(e0 + xq*(e2 + xq*(e4 + xq)))

         elseif ( abx .lt. xlim3 ) then                                   ! humlicek w4 region 3
            if ( rg3 .ne. 0 ) then                                          ! first point in region 3
               rg3 = 0
               z0 = 272.1014     + y*(1280.829 + y*(2802.870 + y*(3764.966&    ! region 3 y-dependents
                    + y*(3447.629 + y*(2256.981 + y*(1074.409 + y*(369.1989&
                    + y*(88.26741 + y*(13.39880 + y)))))))))
               z2 = 211.678      + y*(902.3066 + y*(1758.336 + y*(2037.310&
                    + y*(1549.675 + y*(793.4273 + y*(266.2987&
                    + y*(53.59518 + y*5.0)))))))
               z4 = 78.86585     + y*(308.1852 + y*(497.3014 + y*(479.2576&
                    + y*(269.2916 + y*(80.39278 + y*10.0)))))
               z6 = 22.03523     + y*(55.02933 + y*(92.75679 + y*(53.59518&
                    + y*10.0)))
               z8 = 1.496460     + y*(13.39880 + y*5.0)
               p0 = 153.5168     + y*(549.3954 + y*(919.4955 + y*(946.8970&
                    + y*(662.8097 + y*(328.2151 + y*(115.3772 + y*(27.93941&
                    + y*(4.264678 + y*0.3183291))))))))
               p2 = -34.16955    + y*(-1.322256+ y*(124.5975 + y*(189.7730&
                    + y*(139.4665 + y*(56.81652 + y*(12.79458&
                    + y*1.2733163))))))
               p4 = 2.584042     + y*(10.46332 + y*(24.01655 + y*(29.81482&
                    + y*(12.79568 + y*1.9099744))))
               p6 = -0.07272979  + y*(0.9377051+ y*(4.266322 + y*1.273316))
               p8 = 0.0005480304 + y*0.3183291
            endif
            d = 1.7724538 / (z0 + xq*(z2 + xq*(z4 + xq*(z6 + xq*(z8+xq)))))
            k(i) = d*(p0 + xq*(p2 + xq*(p4 + xq*(p6 + xq*p8))))

         else                                                             ! humlicek cpf12 algorithm
            ypy0 = y + y0
            ypy0q = ypy0*ypy0
            k(i) = 0.0
            do j = 0, 5
               d = x(i) - t(j)
               mq(j) = d*d
               mf(j) = 1.0 / (mq(j) + ypy0q)
               xm(j) = mf(j)*d
               ym(j) = mf(j)*ypy0
               d = x(i) + t(j)
               pq(j) = d*d
               pf(j) = 1.0 / (pq(j) + ypy0q)
               xp(j) = pf(j)*d
               yp(j) = pf(j)*ypy0
            enddo

            if ( abx .le. xlim4 ) then                                      ! humlicek cpf12 region i
               do j = 0, 5
                  k(i) = k(i) + c(j)*(ym(j)+yp(j)) - s(j)*(xm(j)-xp(j))
               enddo

            else                                                            ! humlicek cpf12 region ii
               yf   = y + y0py0
               do j = 0, 5
                  k(i) = k(i) &
                       + (c(j)*(mq(j)*mf(j)-y0*ym(j)) + s(j)*yf*xm(j)) / (mq(j)+y0q) &
                       + (c(j)*(pq(j)*pf(j)-y0*yp(j)) - s(j)*yf*xp(j)) / (pq(j)+y0q)
               enddo
               k(i) = y*k(i) + exp ( -xq )
            endif
         endif
      enddo
    end subroutine voigt




    subroutine resampleRay_tau(L, nTau, tau, dtau_max, maxtau, newL, newNTau,&
         inFlow, newInFlow)
      integer, intent(in) :: nTau, maxtau
      real, intent(in) :: L(:)
      real, intent(in) :: tau(:)
      ! maxmum line optical depth increment allowed
      real, intent(in) :: dtau_max
      integer, intent(inout) :: newNtau
      real, intent(inout)  :: newL(:)
      logical, intent(inout)  :: InFlow(:)
      logical, intent(inout)  :: newInFlow(:)
      !
      real :: dtau
      integer :: nAdd
      integer :: i, j
!      real, parameter :: dvel  = 10.e5/cSpeed
      real ::  dL

      ! Now we devide the line seqment if the optical depth
      ! is too large.

      ! Initial optical depth in each line sequemnt
      newNTau = 0
      do i = 1, nTau-1
         dtau = tau(i+1) - tau(i)
         if (dtau > dtau_max) then
!            if (dtau < 20.0) then
!            if (dtau < 2000.0) then
               nAdd = nint(dtau/dtau_max)  ! this should be at least 1
!            else
!               nAdd = 10  ! to avoid a large number of points added
!            end if
            dL = (L(i+1)-L(i))/real(nAdd+1)
            do j = 1, nAdd+1
               newNtau = newNtau + 1
               if (newntau > maxtau) then
                  print *, "Error: newntau > maxtau in resampleRay_tau. (1)"
                  print *, "       maxtau = ", maxtau
                  print *, "       newtau = ", newNtau
                  stop
               end if
               newL(newNTau) = real(j-1)*dL + L(i)
               newInFlow(newNTau) = inFlow(i)
            enddo
         else
            newNtau = newNtau + 1
            if (newntau > maxtau) then
               print *, "Error: newntau > maxtau in resampleRay_tau. (2)"
               stop
            end if
            newL(newNTau) = L(i)
            newInFlow(newNTau) = inFlow(i)
         endif
      enddo

      newNtau = newNtau + 1
      newL(newNTau) = L(nTau)
      newInFlow(newNTau) = inFlow(nTau)
    end subroutine resampleRay_tau

    subroutine linearResample(xArray, yArray, nX, nx_max, newXarray, newNx)
      real(single) :: xArray(:), yArray(:)
      integer :: nx, newNx
      integer, intent(in) :: nx_max
      real(single) :: newXarray(:)
      real(single), save, allocatable :: newYarray(:) ! newYarray(newNx) ! automatic array
      integer :: i, j
      logical, save :: first_time = .true.
      !$OMP THREADPRIVATE (first_time, newYarray)
      if (first_time) then
         first_time=.false.
         ALLOCATE(newYarray(nx_max))
      end if

      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         if (xArray(j+1) /= xArray(j)) then
            newYarray(i) = yArray(j) + (yArray(j+1)-yArray(j))*(newXarray(i)-xArray(j))/(xArray(j+1)-xArray(j))
         else
            newYarray(i) = yArray(j)
         endif
      enddo
      yArray(1:newNx) = newYArray(1:newNx)
!      deallocate(newYarray)
    end subroutine linearResample





    subroutine linearResample_dble(xArray, yArray, nX, nx_max, newXarray, newNx)
      real(single) :: xArray(:)
      real(double) :: yArray(:)
      integer :: nx, newNx
      integer, intent(in) :: nx_max
      real(single) :: newXarray(:)
      real(single), save, allocatable  :: newYarray(:) ! newYarray(newNx)  ! automatic array
      integer :: i, j
      logical, save :: first_time = .true.
      !$OMP THREADPRIVATE (first_time, newYarray)

      if (first_time) then
         first_time=.false.
         ALLOCATE(newYarray(Nx_max))
      end if


      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         if (xArray(j+1) /= xArray(j)) then
            newYarray(i) = real(yArray(j) + (yArray(j+1)-yArray(j))*(newXarray(i)-xArray(j))/(xArray(j+1)-xArray(j)))
         else
            newYarray(i) = real(yArray(j))
         endif
      enddo
      yArray(1:newNx) = newYArray(1:newNx)
    end subroutine linearResample_dble

  SUBROUTINE trapzd(func,a,b,s,n)
    IMPLICIT NONE
    REAL, INTENT(IN) :: a,b
    REAL, INTENT(INOUT) :: s
    INTEGER, INTENT(IN) :: n

    INTERFACE
      FUNCTION func(x)
        REAL, DIMENSION(:), INTENT(IN) :: x
        REAL, DIMENSION(size(x)) :: func
      END FUNCTION func
    END INTERFACE

    ! This routine computes the nth stage of refinement of an extended trapezoidal rule.
    ! func is input as the name of the function to be integrated between limits a and b,
    ! also input. When called with n=1, the routine returns as s the crudest estimate of
    ! b a f(x)dx. Subsequent calls with n=2,3,... (in that sequential order) will improve
    ! the accuracy of s by adding 2n-2 additional interior points. s should not be modified
    ! between sequential calls.
    REAL :: del,fsum
    INTEGER :: it
    if (n == 1) then
      s=0.5*(b-a)*sum(func( (/ a,b /) ))
    else
      it=2**(n-2)
      del=(b-a)/it ! This is the spacing of the points to be added.
      fsum=sum(func(arth(a+0.5*del,del,it)))
      s=0.5*(s+del*fsum) ! This replaces s by its re ned value.
    end if
  END SUBROUTINE trapzd


  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r
  !BL
  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d
  !BL
  FUNCTION arth_i(first,increment,n)
    INTEGER(I4B), INTENT(IN) :: first,increment,n
    INTEGER(I4B), DIMENSION(n) :: arth_i
    INTEGER(I4B) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

  FUNCTION qsimp(func,a,b)

   IMPLICIT NONE
   REAL, INTENT(IN) :: a,b
   REAL :: qsimp
   INTERFACE
     FUNCTION func(x)
       REAL, DIMENSION(:), INTENT(IN) :: x
       REAL, DIMENSION(size(x)) :: func
     END FUNCTION func
   END INTERFACE
   INTEGER, PARAMETER :: JMAX=20
   REAL, PARAMETER :: EPS=2.0e-4
   ! Returns the integral of the function func from a to b. The parameter EPS should be
   ! set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1
   ! is the maximum allowed number of steps. Integration is performed by Simpson's
   ! rule.
   INTEGER :: j
   REAL :: os,ost,st
   ost=0.0
   os= 0.0
   do j=1,JMAX
     call trapzd(func,a,b,st,j)
     qsimp=(4.0*st-ost)/3.0  ! Compare equation (4.2.4).
     if (j > 5) then !Avoid spurious early convergence.
       if (abs(qsimp-os) < EPS*abs(os) .or. (qsimp == 0.0 .and. os == 0.0)) RETURN
     end if
     os=qsimp
     ost=st
   end do
   print *, 'Too many steps in qsimp'
   stop

  END FUNCTION qsimp

  subroutine stripSimilarValuesSingle(x, nx, xtol)
    integer :: nx
    real :: x(:)
    real :: xtol
    integer :: i, newNx



    call sort(nx, x)

    newnx = 1
    do i = 2, nx
       if (abs(x(newnx)-x(i)) > xTol) then
          newnx = newnx  + 1
          x(newnx) = x(i)
       endif
    enddo

    nx = newnx
  
  end subroutine stripSimilarValuesSingle

  subroutine stripSimilarValuesDouble(x, nx, xtol)
    integer :: nx
    real(double) :: x(:)
    real(double) :: xtol
    real(double), allocatable :: xtemp(:)
    integer :: i, newNx


    call sort(nx, x)

    newnx = 1
  
    do i = 2, nx
       if (abs(x(newnx)-x(i)) > xTol) then
          newnx = newnx  + 1
          x(newnx) = x(i)
       endif
    enddo

  
    nx = newnx
   
  end subroutine stripSimilarValuesDouble


  subroutine convertByte4(iByte, ival)
    integer(kind=single) :: ibyte(4)
    integer(kind=single) :: ubyte(4)
    integer  :: ival

    ubyte = ibyte
    where (ubyte < 0) ubyte = ibyte + 256


    ival = ubyte(1) + 256*ubyte(2) + 65536 * ubyte(3) + 16777200 * ubyte(4)

  end subroutine convertByte4

  subroutine convertByte2(iByte, ival)
    integer(kind=single) :: ibyte(2)
    integer(kind=single) :: ubyte(2)
    integer(kind=double)  :: ival

    ubyte = ibyte
    where (ubyte < 0) ubyte = ibyte + 256
    ival = ubyte(1) + 256*ubyte(2)

  end subroutine convertByte2

  subroutine wavenumbertoEv(array, narray)
    integer :: nArray
    real :: array(:), wavelength
    integer :: i
    real :: freq,energy
    do i = 1, nArray
       if (array(i) /= 0.) then
          wavelength = 1./array(i)
          freq = real(cSpeed / wavelength)
          energy = real(freq * hcgs)
          array(i) = real(energy * ergtoEv)
       else
          array(i) = 1.e-20
       endif
    enddo
  end subroutine wavenumbertoEv


subroutine createRoman(n, roman)
  integer :: n
  character(len=*) :: roman

  select case(n)
     case(1)
        roman = "I"
     case(2)
        roman = "II"
     case(3)
        roman = "III"
     case(4)
        roman = "IV"
     case(5)
        roman = "V"
     case(6)
        roman = "VI"
     case(7)
        roman = "VII"
     case(8)
        roman = "VIII"
     case(9)
        roman = "IX"
     case(10)
        roman = "X"
     case(11)
        roman = "XI"
     case(12)
        roman = "XII"
     case(13)
        roman = "XIII"
     case(14)
        roman = "XIV"
     case(15)
        roman = "XV"
     case DEFAULT
        write(*,*) "Can't create roman numeral for: ",n
        roman = "???"
  end select
end subroutine createRoman

subroutine returnElement(z, element)
  character(len=2) :: element
  integer :: z
  select case(z)
     case(1)
        element = "H"
     case(2)
        element = "He"
     case(3)
        element = "Li"
     case(4)
        element = "Be"
     case(5)
        element = "B"
     case(6)
        element = "C"
     case(7)
        element = "N"
     case(8)
        element = "O"
     case(9)
        element = "F"
     case(10)
        element = "Ne"
     case(11)
        element = "Na"
     case(12)
        element = "Mg"
     case(13)
        element = "Al"
     case(14)
        element = "Si"
     case(16)
        element = "S"
     case(18)
        element = "Ar"
     case(20)
        element = "Ca"
     case(26)
        element = "Fe"
     case DEFAULT
        write(*,*) "returnElement :: Unknown element: ",z
        element = "??"
  end select
end subroutine returnElement


SUBROUTINE GAUSSJ(A,N,NP,B,M,MP, ok)
  implicit none
  integer, parameter :: nmax=100
  integer :: n, m, np, mp, i, j, k, l, icol, irow
  real(double) :: dum, pivinv
  integer :: ll
  real(double) A(NP,NP),B(NP,MP), big
  integer :: IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
  logical :: ok
  logical, save :: firsttime = .true.
  !$OMP THREADPRIVATE (firstTime)

  ok = .true.
  DO J=1,N
     IPIV(J)=0
  enddo
  DO I=1,N
     BIG=0.d0
     DO  J=1,N
        IF(IPIV(J).NE.1)THEN
           DO  K=1,N
              IF (IPIV(K).EQ.0) THEN
                 IF (ABS(A(J,K)).GE.BIG)THEN
                    BIG=ABS(A(J,K))
                    IROW=J
                    ICOL=K
                 ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                 write(*,*) '! Singular matrix: pivot'
                 ok = .false.
                 goto 666
              ENDIF
           enddo
        ENDIF
     enddo
     IPIV(ICOL)=IPIV(ICOL)+1
     IF (IROW.NE.ICOL) THEN
        DO L=1,N
           DUM=A(IROW,L)
           A(IROW,L)=A(ICOL,L)
           A(ICOL,L)=DUM
        enddo
        DO  L=1,M
           DUM=B(IROW,L)
           B(IROW,L)=B(ICOL,L)
           B(ICOL,L)=DUM
        enddo
     ENDIF
     INDXR(I)=IROW
     INDXC(I)=ICOL
     IF (A(ICOL,ICOL).EQ.0.d0) then
        if (firsttime) then
           write(*,*) 'Singular matrix.',icol
           stop
           firsttime = .false.
        endif
        ok = .false.
        goto 666
     endif
     PIVINV=1.D0/A(ICOL,ICOL)
     A(ICOL,ICOL)=1.d0
     DO L=1,N
        A(ICOL,L)=A(ICOL,L)*PIVINV
     enddo
     DO L=1,M
        B(ICOL,L)=B(ICOL,L)*PIVINV
     enddo
     DO LL=1,N
        IF(LL.NE.ICOL)THEN
           DUM=A(LL,ICOL)
           A(LL,ICOL)=0.d0
           DO L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
           enddo
           DO L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
           enddo
        ENDIF
     enddo
  enddo
  DO  L=N,1,-1
     IF(INDXR(L).NE.INDXC(L))THEN
        DO K=1,N
           DUM=A(K,INDXR(L))
           A(K,INDXR(L))=A(K,INDXC(L))
           A(K,INDXC(L))=DUM
        enddo
     ENDIF
  enddo
666 continue
END SUBROUTINE GAUSSJ

  ! this procedure performs the solution of linear equations
  subroutine luSlv(a, b, ok)
    implicit none

    real(double), intent(inout) :: a(:,:)
    real(double), intent(inout) :: b(:)
    logical, intent(out) :: ok
    integer :: i, n
    logical, save :: firstTime = .true.
    !$OMP THREADPRIVATE(firstTime)


! lured can change diagonal terms in a so check for
! zeros in the subroutine iteslf
    call lured(a,ok)
    if (.not. ok) then
       if (firstTime) then
          write(*,*) "LU solver: lured not ok"
          firstTime = .false.
       endif
       write(*,*) a
       return
    endif

! reslv only modifies b so check for zeros here
    n = size(a,1)
    do i=1,n
       if (a(i,i)==0.0) then
          ok=.false.
          write(*,*) "LU solver: zeros after lured"
          return
       endif
    end do

    call reslv(a,b)

  end subroutine luSlv

  subroutine lured(a,ok)
    implicit none

    real(double), intent(inout)  :: a(:,:)
    logical, intent(out) :: ok

    ! local variables
    integer          :: i, j, k, n                    ! counters

    real(double) :: factor                     ! general calculation factor

    ok = .true.
    n = size(a,1)

    if (n == 1) return

i_loop: do i = 1, n-1
       do k = i+1, n
          if (a(i,i)==0.d0) then
             ok=.false.
             exit i_loop
          endif
          factor = a(k,i)/a(i,i)
          do j = i+1, n
             a(k, j) = a(k, j) - a(i, j) * factor
          end do
       end do
    end do i_loop

  end subroutine lured

  subroutine reslv(a,b)
    implicit none

    real(double), intent(inout) :: a(:,:)
    real(double), intent(inout) :: b(:)

    ! local variables
    integer    :: i, j, k, l              ! counters
    integer    :: n

    n = size(b)

    if (n == 1) then
       b(n) = b(n) / a(n,n)
       return
    end if

    do i = 1, n-1
       do j = i+1, n
          b(j) = b(j) - b(i)*a(j, i)/ a(i, i)
       end do
    end do
    b(n) = b(n) / a(n,n)
    do i = 1, n-1
       k = n-i
       l = k+1
       do j = l, n
          b(k) = b(k) - b(j)*a(k, j)
       end do
       b(k) = b(k) / a(k,k)
    end do
  end subroutine reslv

  subroutine insertBin(xArray, nx, x, dx)
    real :: x, xarray(:)
    integer :: nx
    real :: dx !, xplus1, xplus2
    integer :: i, j

    if ((x > xArray(1)).and.(x < xArray(nx))) then
!       call locate(xArray, nx, x, i)
!       xplus2 = min(xArray(i+1)-0.01,x+dx/2.)
!       xplus1 = max(xArray(i)+0.01,x-dx/2.)
       xArray(nx+1) = x+dx/2.
       xArray(nx+2) = x-dx/2.
       xArray(nx+3) = x
       nx = nx + 3
       call sort(nx, xArray)
       i = 1
       do while (i < nx)
          if (xArray(i) == xArray(i+1)) then
             do j = i, nx-1
                xArray(j) = xArray(j+1)
             enddo
             nx = nx - 1
          else
             i = i + 1
          endif
       enddo
    endif
  end subroutine insertBin


  PURE SUBROUTINE INDEXX_int(N,ARRIN,INDX)
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(IN)     :: ARRIN(:)
    INTEGER, INTENT(OUT) :: INDX(:)
    INTEGER              :: J, L, IR, I, INDXT
    REAL                 :: Q
    DO  J=1,N
       INDX(J)=J
    ENDDO
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       INDXT=INDX(L)
       Q=ARRIN(INDXT)
    ELSE
       INDXT=INDX(IR)
       Q=ARRIN(INDXT)
       INDX(IR)=INDX(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          INDX(1)=INDXT
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
       ENDIF
       IF(Q.LT.ARRIN(INDX(J)))THEN
          INDX(I)=INDX(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    INDX(I)=INDXT
    GO TO 10
  END subroutine indexx_int

  PURE SUBROUTINE INDEXX_single(N,ARRIN,INDX)
    INTEGER, INTENT(IN)  :: N
    REAL, INTENT(IN)     :: ARRIN(:)
    INTEGER, INTENT(OUT) :: INDX(:)
    INTEGER              :: J, L, IR, I, INDXT
    REAL                 :: Q
    DO  J=1,N
       INDX(J)=J
    ENDDO
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       INDXT=INDX(L)
       Q=ARRIN(INDXT)
    ELSE
       INDXT=INDX(IR)
       Q=ARRIN(INDXT)
       INDX(IR)=INDX(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          INDX(1)=INDXT
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
       ENDIF
       IF(Q.LT.ARRIN(INDX(J)))THEN
          INDX(I)=INDX(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    INDX(I)=INDXT
    GO TO 10
  END subroutine indexx_single

  PURE SUBROUTINE INDEXX_double(N,ARRIN,INDX)
    INTEGER, INTENT(IN)  :: N
    REAL(double), INTENT(IN)     :: ARRIN(:)
    INTEGER, INTENT(OUT) :: INDX(:)
    INTEGER              :: J, L, IR, I, INDXT
    REAL(double)                 :: Q
    DO  J=1,N
       INDX(J)=J
    ENDDO
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       INDXT=INDX(L)
       Q=ARRIN(INDXT)
    ELSE
       INDXT=INDX(IR)
       Q=ARRIN(INDXT)
       INDX(IR)=INDX(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          INDX(1)=INDXT
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
       ENDIF
       IF(Q.LT.ARRIN(INDX(J)))THEN
          INDX(I)=INDX(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    INDX(I)=INDXT
    GO TO 10
  END subroutine indexx_double


  subroutine  convertToFnu(nuarray, fArray, n)
    real(double) :: nuArray(:), fArray(:), nu
    real, allocatable :: nuTemp(:), fTemp(:)
    integer :: i, n


    allocate(nuTemp(1:n), fTemp(1:n))
    if (nuArray(1) < 1.d5) then
       call writeWarning("Surface flux file appears to be in wavelength space: converting")
       do i = 1, n
          nu = cSpeed / (1.d-8 * nuArray(i)) ! hz
          nuTemp(i) = real(nu)
          fTemp(i) = real(fArray(i)*1.d8) ! erg/s/cm^2/angs -> erg/s/cm^2/cm
          fTemp(i) = real(fTemp(i) * cspeed / nu**2) ! erg/s/cm^2/Hz
       end do
       do i = 1, n
          nuArray(i) = nuTemp(n-i+1)
          fArray(i) = fTemp(n-i+1)
       enddo
       deallocate(fTemp, nuTemp)
!       do i = 1, n
!          write(45,*) nuarray(i), fArray(i)
!       enddo
!       stop
    endif
  end subroutine convertToFnu

  subroutine regular_tri_quadint(t1,t2,t3,weights)

    real(double) :: t1, t2, t3
    real(double) :: xu,xv,xw,yu,yv,yw,zu,zv,zw
    real(double) :: xbasis1, xbasis2, ybasis1, ybasis2, zbasis1, zbasis2
    real(double) :: xyweights(9)
    real(double) :: weights(27)
! OpenMP has problems with saved variables so set reuse=.false.


    xbasis1 = t1 - 1.d0
    xbasis2 = 2.d0 * t1 - 1.d0
    ybasis1 = t2 - 1.d0
    ybasis2 = 2.d0 * t2 - 1.d0
    zbasis1 = t3 - 1.d0
    zbasis2 = 2.d0 * t3 - 1.d0

    xu = xbasis1 * xbasis2
    xv = -4.d0 * xbasis1 * t1
    xw = xbasis2 * t1

    yu = ybasis1 * ybasis2
    yv = -4.d0 * ybasis1 * t2
    yw = ybasis2 * t2

    zu = zbasis1 * zbasis2
    zv = -4.d0 * zbasis1 * t3
    zw = zbasis2 * t3

    ! base levels

    xyweights(:) = (/ xu*yu, xv*yu, xw*yu, xu*yv, xv*yv, xw*yv, xu*yw, xv*yw, xw*yw /)

    weights(1:9) = zu * xyweights
    weights(10:18) = zv * xyweights
    weights(19:27) = zw * xyweights


  end subroutine regular_tri_quadint

   function general_quadint(xi,yi,x) result(y)
     implicit none

     real(double), intent(in) :: xi(:),yi(:),x
     real(double) :: y
     real(double) :: xxi(3),a(3), det

     xxi = xi * xi

     det = xxi(1)*(xi(2)-xi(3))+&
           xxi(2)*(xi(3)-xi(1))+&
           xxi(3)*(xi(1)-xi(2))

     a(1) = yi(1)*(xi(2)-xi(3))+&
            yi(2)*(xi(3)-xi(1))+&
            yi(3)*(xi(1)-xi(2))

     a(2) = yi(1)*(xxi(3)-xxi(2))+&
            yi(2)*(xxi(1)-xxi(3))+&
            yi(3)*(xxi(2)-xxi(1))

     a(3) = yi(1)*(xxi(2)*xi(3)-xxi(3)*xi(2))+&
            yi(2)*(xxi(3)*xi(1)-xxi(1)*xi(3))+&
            yi(3)*(xxi(1)*xi(2)-xxi(2)*xi(1))

     a = a / det

     y = ((a(1) * x) + a(2)) * x + a(3)

   end function general_quadint

      SUBROUTINE LINFIT(X,Y,SIGMAY,NPTS,MODE,A,SIGMAA,B,SIGMAB,R)
        implicit none
        integer :: npts, mode, i
        real(double) ::  X(NPTS),Y(NPTS),SIGMAY(NPTS)
        real(double) :: sum, sumx, sumy, sumx2,sumy2, sumxy, xi, yi, weight
        real(double) :: c, varnce, delta
        real(double), intent(out) :: a, b, sigmaa, sigmab,r
!
!     MAKES LEAST-SQUARES FIT TO DATA WITH A STRAIGHT LINE  Y = A + B*X
!
!     X      - ARRAY OF DATA POINTS FOR INDEPENDENT VARIABLE
!     Y      - ARRAY OF DATA POINTS FOR DEPENDENT VARIABLE
!     SIGMAY - ARRAY OF STANDARD DEVIATIONS FOR Y DATA POINTS
!     NPTS   - NUMBER OF PAIRS OF DATA POINTS
!     MODE   - DETERMINES METHOD OF WEIGHTING LEAST-SQUARES FIT
!              +1 (INSTRUMENTAL) WEIGHT(I) = 1/SIGMAY(I)**2
!              0  (NO WEIGHTING) WEIGHT(I) = 1
!              -1 (STATISTICAL)  WEIGHT(I) = 1/Y(I)
!     A      - Y INTERCEPT OF FITTED STRAIGHT LINE
!     SIGMAA - STANDARD DEVIATION OF A
!     B      - SLOPE OF FITTED STRAIGHT LINE
!     SIGMAB - STANDARD DEVIATION OF B
!     R      - LINEAR CORRELATION COEFFICIENT
!
!     ACCUMULATE WEIGHTED SUMS
!
      SUM = 0.
      SUMX = 0.
      SUMY = 0.
      SUMX2 = 0.
      SUMY2 = 0.
      SUMXY = 0.
      VARNCE = 1.
      DO I = 1,NPTS
        XI = X(I)
        YI = Y(I)
        IF (MODE .LT. 0) THEN
          IF (YI .LT. 0) THEN
            WEIGHT = 1/(-YI)
          ELSEIF (YI .EQ. 0) THEN
            WEIGHT = 1
          ELSE
            WEIGHT = 1/YI
          ENDIF
        ELSEIF (MODE .EQ. 0) THEN
          WEIGHT = 1
        ELSE
          WEIGHT = 1/SIGMAY(I)**2
        ENDIF
        SUM = SUM + WEIGHT
        SUMX = SUMX + WEIGHT*XI
        SUMY = SUMY + WEIGHT*YI
        SUMX2 = SUMX2 + WEIGHT*XI*XI
        SUMY2 = SUMY2 + WEIGHT*YI*YI
        SUMXY = SUMXY + WEIGHT*XI*YI
      ENDDO
!
!     CALCULATE COEFFICIENTS AND STANDARD DEVIATIONS
!
      DELTA = SUM*SUMX2 - SUMX*SUMX
      A = (SUMX2*SUMY-SUMX*SUMXY)/DELTA
      B = (SUMXY*SUM-SUMX*SUMY)/DELTA
      IF (MODE .NE. 0) THEN
        VARNCE = 1.
      ELSE
        IF (NPTS .GT. 2) THEN
           C = NPTS - 2
           VARNCE = (SUMY2+A*A*SUM+B*B*SUMX2 &
                -2*(A*SUMY+B*SUMXY-A*B*SUMX))/C
        ENDIF
      ENDIF
      SIGMAA = SQRT(ABS(VARNCE*SUMX2/DELTA))
      SIGMAB = SQRT(ABS(VARNCE*SUM/DELTA))
      R = (SUM*SUMXY-SUMX*SUMY)/SQRT(DELTA*(SUM*SUMY2-SUMY*SUMY))
    END SUBROUTINE LINFIT

!!! Vector sequence accleration routine by Ng J. Chem. Phys. (61) 7, 1974 and OAB JQSRT (35) 431, 1986
!!! coded by DAR Feb 09
!!! final version? April 2010 - got rid of all of the allocating and deallocating.
!!! Added capability to refine on just part of the vector
subroutine ngStep(out, qorig, rorig, sorig, torig, weight, doubleweight, length)

  implicit none

!#ifdef USEIEEEISNAN
!  use ieee_arithmetic, isnan=> ieee_is_nan
!#endif

    integer :: length

    real(double) :: q(length), r(length), s(length), t(length), oneoverT(length)
    real(double), optional :: weight(:)
    real(double) :: diff1(length), diff2(length)
    real(double) :: diff01(length), diff02(length), diff(length)
    real(double) :: a,b,d,e,f ! matrix coefficients - matrix is symmetric (at least)
    real(double) :: c1, c2
    real(double) :: out(:)
    real(double) :: qorig(:), rorig(:), sorig(:), torig(:)
    real(double) :: det
    real(double) :: OABweight(length)

    integer :: vectorsize
    logical, optional :: doubleweight
    logical :: dodoubleweight

    vectorsize = size(torig)
    out = torig
! Catch any nasty 0s
    if(any(torig(1:length) .eq. 0.d0)) then
       out = torig
       goto 666
    endif

    if(present(doubleweight)) then
       dodoubleweight = doubleweight
    else
       dodoubleweight = .false.
    endif

    if(present(weight)) then
       OABweight(1:length) = weight(1:length)
    else
       OABweight(1:length) = 1.d0
    endif

    q = qorig(1:length)
    r = rorig(1:length)
    s = sorig(1:length)
    t = torig(1:length)

    OneOverT(1:length) = 1.d0 / t(1:length)

    q(1:length) = q(1:length) * OneOverT(1:length)
    r(1:length) = r(1:length) * OneOverT(1:length)
    s(1:length) = s(1:length) * OneOverT(1:length)
    t(1:length) = 1.d0

      if(dodoubleweight) then
         q = q * 0.125d0
         r = r * 0.25d0
         s = s * 0.5d0
         t = t
      endif

      diff  = (t - s)
      diff1 = (s - r)
      diff2 = (r - q)

      diff01 = (diff - diff1)
      diff02 = (diff - diff2)

      a = sum(diff01**2 * OABweight)
      b = sum(diff01 * diff02 * OABweight)
      d = sum(diff02**2 * OABweight)
      e = sum(diff * diff01 * OABweight)
      f = sum(diff * diff02 * OABweight)

      det = a*d - b*b

      if(det .eq. 0.d0) then
         out = torig
         goto 666
      endif

      c1 = (e*d-b*f) / det
      c2 = (a*f-e*b) / det

      out(1:length) = (1d0 - c1 - c2) * torig(1:length) + c1 * sorig(1:length) + c2 * rorig(1:length)

! Uncomment this code if you find yourself trapping NaNs
!      if(any(isnan(out))) then
!         outarray(1:vectorsize) = torig(1:vectorsize)
!         out => outarray(1:vectorsize)
!      endif
666 continue
    end subroutine ngStep

  !********************************************************************************
  !** FICHE F.35.  THE VORONOI CONSTRUCTION IN 2D AND 3D.                        **
  !** This FORTRAN code is intended to illustrate points made in the text.       **
  !** To our knowledge it works correctly.  However it is the responsibility of  **
  !** the user to test it, if it is to be used in a research application.        **
  !********************************************************************************

  !    *******************************************************************
  !    ** TWO SEPARATE PARTS: TWO AND THREE DIMENSIONAL VERSIONS.       **
  !    *******************************************************************



  !    *******************************************************************
  !    ** FICHE F.35  -  PART A                                         **
  !    ** THE VORONOI CONSTRUCTION IN 2D.                               **
  !    *******************************************************************



  SUBROUTINE VORON2(N, RX, RY, BOX, AREA, success)
    implicit none

    !    *******************************************************************
    !    ** CONSTRUCTION OF THE VORONOI POLYGON IN 2D.                    **
    !    **                                                               **
    !    ** THIS PROGRAM TAKES IN A CONFIGURATION IN A SQUARE BOX WITH    **
    !    ** CONVENTIONAL PERIODIC BOUNDARY CONDITIONS AND FOR EACH ATOM   **
    !    ** OBTAINS THE SURROUNDING VORONOI POLYGON, DEFINED AS THAT      **
    !    ** REGION OF SPACE CLOSER TO THE CHOSEN ATOM THAN TO ANY OTHER.  **
    !    ** NEIGHBOURING POLYGONS DEFINE NEIGHBOURING ATOMS.              **
    !    ** THE PROGRAM IS SLOW BUT ESSENTIALLY FOOLPROOF.                **
    !    ** WE USE THE MINIMUM IMAGE CONVENTION AND SET A CUTOFF BEYOND   **
    !    ** WHICH ATOMS ARE ASSUMED NOT TO BE NEIGHBOURS: BOTH OF THESE   **
    !    ** MEASURES ARE DANGEROUS FOR SMALL AND/OR RANDOM SYSTEMS.       **
    !    ** WE DELIBERATELY DO NOT USE PREVIOUSLY-FOUND NEIGHBOURS IN     **
    !    ** CONSTRUCTING NEIGHBOUR LISTS, SO THAT AN INDEPENDENT CHECK    **
    !    ** MAY BE MADE AT THE END.                                       **
    !    ** HERE WE SIMPLY PRINT OUT THE GEOMETRICAL INFORMATION AT THE   **
    !    ** END.  THE OUTPUT IS QUITE LENGTHY.  IN PRACTICE, IT WOULD     **
    !    ** PROBABLY BE ANALYZED DIRECTLY WITHOUT PRINTING OUT.           **
    !    ** NB: BEWARE DEGENERATE CONFIGURATIONS, I.E. ONES IN WHICH MORE **
    !    ** THAN THREE VORONOI DOMAINS SHARE A VERTEX. THE SQUARE LATTICE **
    !    ** IS AN EXAMPLE.                                                **
    !    **                                                               **
    !    ** PRINCIPAL VARIABLES:                                          **
    !    **                                                               **
    !    ** INTEGER N                        NUMBER OF ATOMS              **
    !    ** REAL    RX(N),RY(N)              POSITIONS                    **
    !    ** REAL    PX(MAXCAN),PY(MAXCAN)    CANDIDATE RELATIVE POSITIONS **
    !    ** REAL    PS(MAXCAN)               SQUARED RELATIVE DISTANCES   **
    !    ** INTEGER NVER                     NUMBER OF VERTICES FOUND     **
    !    ** INTEGER NEDGE                    NUMBER OF EDGES FOUND        **
    !    ** INTEGER VERTS(MAXCAN)            VERTICES FOR EACH CANDIDATE  **
    !    **                                  = 0 IF NOT A NEIGHBOUR       **
    !    **                                  = 2 ( 1 EDGE ) IF NEIGHBOUR  **
    !    ** REAL    RXVER(MAXVER)            VERTEX RELATIVE X-COORD      **
    !    ** REAL    RYVER(MAXVER)            VERTEX RELATIVE Y-COORD      **
    !    ** INTEGER IVER(MAXVER)             ATOMIC INDICES TAGGING       **
    !    ** INTEGER JVER(MAXVER)             .. EACH VERTEX OF POLYGON    **
    !    **                                                               **
    !    ** ROUTINES REFERENCED:                                          **
    !    **                                                               **
    !    ** SUBROUTINE READCN ( CNFILE, N, BOX )                          **
    !    **    READS IN CONFIGURATION, NUMBER OF ATOMS, BOX SIZE          **
    !    ** SUBROUTINE SORT ( MAXCAN, PX, PY, PS, TAG, NCAN )             **
    !    **    SORTS NEIGHBOUR DETAILS INTO ASCENDING DISTANCE ORDER      **
    !    ** SUBROUTINE WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE,          **
    !    **     PX, PY, PS, VERTS, RXVER, RYVER, IVER, JVER )             **
    !    **    CARRIES OUT THE VORONOI CONSTRUCTION                       **
    !    *******************************************************************

    INTEGER     MAXN, MAXCAN, MAXVER
    PARAMETER ( MAXN = 10000, MAXCAN = 10000, MAXVER = 500 )
    logical :: success
    REAL(double)        RX(:), RY(:), AREA(:)

    REAL(double)        PX(MAXCAN), PY(MAXCAN), PS(MAXCAN)
    INTEGER     TAG(MAXCAN), VERTS(MAXCAN)

    REAL(double)        RXVER(MAXVER), RYVER(MAXVER)
    INTEGER     IVER(MAXVER), JVER(MAXVER)
    INTEGER     NNAB(MAXN), INAB, JNAB

! NABLST has been made allocatable to avoid segfaults with OpenMP
    integer, allocatable :: NABLST(:,:)

    INTEGER     NCAN, NVER, NCOORD, NEDGE
    INTEGER     I, J, CAN, VER, N
    REAL(double)        BOX, BOXINV, RCUT, RCUTSQ, COORD
    REAL(double)        RXJ, RYJ, RXIJ, RYIJ, RIJSQ
    real(double)        xc(maxn), yc(maxn), xp(maxn), yp(maxn), a,totArea
!    CHARACTER   CNFILE*30
    LOGICAL     OK

    !    *******************************************************************

!    WRITE(*,'(1H1,'' **** PROGRAM VORON2 ****                  '')')
!    WRITE(*,'(//1X,''VORONOI CONSTRUCTION IN 2D                '')')

    !    ** BASIC PARAMETERS **

!    WRITE(*,'('' ENTER CONFIGURATION FILENAME                  '')')
!    READ (*,'(A)') CNFILE
!    WRITE(*,'('' CONFIGURATION FILENAME '',A)') CNFILE

    !    ** READCN MUST READ IN INITIAL CONFIGURATION  **

!    CALL READCN ( CNFILE, N, BOX )
    RCUT = 1.d10
!    WRITE(*,'(1X,I5,''-ATOM CONFIGURATION'')') N
!    WRITE(*,'('' BOX LENGTH = '',F10.5)') BOX
!    WRITE(*,'('' ENTER NEIGHBOUR CUTOFF IN SAME UNITS '')')
!    READ (*,*) RCUT
!    WRITE(*,'('' NEIGHBOUR CUTOFF = '',F10.5)') RCUT

    RCUTSQ = RCUT ** 2
    BOXINV = 1.0 / BOX

    !    ** ZERO ACCUMULATORS **


    NNAB = 0

    allocate ( NABLST(MAXVER,MAXN) )
    NABLST = 0

    !    *******************************************************************
    !    ** MAIN LOOP STARTS                                              **
    !    *******************************************************************

    totArea = 0.

!    if (writeoutput) open(33, file="plot.gnu",status="unknown",form="formatted")

    if (N > MAXN) call writeFatal("VORON2: N > MAXN")

    DO J = 1, N

       IF ( MOD ( J, 2 ) .EQ. 0 ) THEN

!          WRITE(*,'(///1X,''RESULTS FOR ATOM '',I5)') J

       ELSE

!          WRITE(*,'(1H1,''RESULTS FOR ATOM '',I5)') J

       ENDIF

       RXJ = RX(J)
       RYJ = RY(J)
       CAN = 0

       !       ** SELECT CANDIDATES **

       DO I = 1, N

          IF ( I .NE. J ) THEN

             RXIJ = RX(I) - RXJ
             RYIJ = RY(I) - RYJ
             !                 RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
             !                 RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
             RIJSQ  = RXIJ ** 2 + RYIJ ** 2

             IF ( RIJSQ .LT. RCUTSQ ) THEN

                CAN = CAN + 1

                IF ( CAN .GT. MAXCAN ) THEN

                   WRITE(*,'('' TOO MANY CANDIDATES '')')
                   write(*,*) can,maxcan
                   STOP

                ENDIF

                PX(CAN)  = RXIJ
                PY(CAN)  = RYIJ
                PS(CAN)  = RIJSQ
                TAG(CAN) = I

             ENDIF

          ENDIF

       enddo

          CAN = CAN + 1
          PX(CAN) = -2.*RX(j)
          PY(CAN) = 0.d0
          PS(CAN) = PX(CAN)**2 + PY(CAN)**2
          TAG(CAN) = J

          CAN = CAN + 1
          PX(CAN) = 0.d0
          PY(CAN) = -2.*RY(j)
          PS(CAN) = PX(CAN)**2 + PY(CAN)**2
          TAG(CAN) = J

          CAN = CAN + 1
          PX(CAN) = 2.d0*(BOX - RX(J))
          PY(CAN) = 0.d0
          PS(CAN) = PX(CAN)**2 + PY(CAN)**2
          TAG(CAN) = J

          CAN = CAN + 1
          PX(CAN) = 0.d0
          PY(CAN) = 2.d0*(BOX - RY(J))
          PS(CAN) = PX(CAN)**2 + PY(CAN)**2
          TAG(CAN) = J



          !       ** CANDIDATES HAVE BEEN SELECTED **

          NCAN = CAN

          !       ** SORT INTO INCREASING DISTANCE ORDER **
          !       ** THIS SHOULD IMPROVE EFFICIENCY      **

          CALL SORT_VORON ( MAXCAN, PX, PY, PS, TAG, NCAN )

          !       ** PERFORM VORONOI CONSTRUCTION **

          CALL WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE, &
                           PX, PY, PS, VERTS, &
                           RXVER, RYVER, IVER, JVER , SUCCESS)

          if (.not.success) then
!             do i = 1, n
!                write(*,*) "rx ",rx(i),ry(i)
!             enddo
             goto 666
          endif

          !       ** WRITE OUT RESULTS **

!          WRITE(*,'(/1X,''NUMBER OF NEIGHBOURS '',I5)') NEDGE
!          WRITE(*,'(/1X,''NEIGHBOUR LIST '')')
!          WRITE(*,10001)



          DO 800 CAN = 1, NCAN

             IF ( VERTS(CAN) .NE. 0 ) THEN

                PS(CAN) = SQRT ( PS(CAN) )
!                WRITE(*,'(1X,I5,3X,I5,3X,2F12.5,3X,F12.5)') &
!                              TAG(CAN), VERTS(CAN), PX(CAN), PY(CAN), PS(CAN)
                NNAB(J) = NNAB(J) + 1
                NABLST(NNAB(J),J) = TAG(CAN)

             ENDIF

800          CONTINUE

!             WRITE(*,'(/1X,''NUMBER OF VERTICES   '',I5)') NVER
!             WRITE(*,'(/1X,''VERTEX LIST '')')
!             WRITE(*,10002)

             DO VER = 1, NVER

!                WRITE(*,'(1X,2I5,3X,2F12.5)') &
!                           TAG(IVER(VER)), TAG(JVER(VER)), &
!                           RXVER(VER), RYVER(VER)

                xp(ver) = rx(j)+rxver(ver)
                yp(ver) = ry(j)+ryver(ver)

             enddo

             call  sortConvex(nver, xp, yp, xc, yc)
!             if (writeoutput) then
!                write(33,*) "plot ""-"" using 1:2 w l"
!                do ver = 1, nver
!                   write(33,*) xc(ver), yc(ver)
!                enddo
!                write(33,*) xc(1), yc(1)
!                write(33,*) "end"
!                write(33,*) " "
!             endif
             call areaPolygon(nver,xc,yc,a)
             totArea = totArea + a
             area(j) = a

          enddo
!          if (writeoutput) then
!             close(33)
!             stop
!          endif

!          write(*,*) "total area ", totArea

          !    *******************************************************************
          !    ** MAIN LOOP ENDS                                                **
          !    *******************************************************************

!          WRITE(*,'(1H1,''FINAL SUMMARY'')')
!          WRITE(*,10003)

          NCOORD = 0

          DO J = 1, N

             NCOORD = NCOORD + NNAB(J)

!             WRITE(*,'(1X,I5,3X,I5,3X,30I3)') J, NNAB(J), &
!                    ( NABLST(INAB,J), INAB = 1, NNAB(J) )

             !       ** CHECK THAT IF I IS A NEIGHBOUR OF J **
             !       ** THEN J IS ALSO A NEIGHBOUR OF I     **

             DO INAB = 1, NNAB(J)

                I = NABLST(INAB,J)

                OK = .FALSE.
                JNAB = 1

1200            IF ( ( .NOT. OK ) .AND. ( JNAB .LE. NNAB(I) ) ) THEN

                   OK = ( J .EQ. NABLST(JNAB,I) )
                   JNAB = JNAB + 1
                   GOTO 1200

                ENDIF

                IF ( .NOT. OK ) THEN

!                   WRITE(*,'(1X,I3,'' IS NOT A NEIGHBOUR OF '',I3)') J, I

                ENDIF

             enddo

          enddo

          COORD = dble ( NCOORD ) / dble ( N )

!          WRITE(*,'(/1X,'' AVERAGE COORDINATION NUMBER = '',F10.5)') COORD


!10001     FORMAT(/1X,'ATOM ',3X,'EDGE ', &
!                   /1X,'INDEX',3X,'VERTS',3X, &
!                   '      RELATIVE POSITION   ',3X,'  DISTANCE  ')
!10002     FORMAT(/1X,'   INDICES         RELATIVE POSITION ')
!10003     FORMAT(/1X,'INDEX    NABS    ... NEIGHBOUR INDICES ... ')

666 continue

          deallocate(NABLST)

       end subroutine voron2

!!$       SUBROUTINE READCN ( CNFILE, N, BOX )
!!$
!!$         COMMON / BLOCK1 / RX, RY
!!$
!!$         !    *******************************************************************
!!$         !    ** SUBROUTINE TO READ IN INITIAL CONFIGURATION                   **
!!$         !    *******************************************************************
!!$
!!$         INTEGER     MAXN
!!$         PARAMETER ( MAXN = 1080 )
!!$
!!$         REAL        RX(MAXN), RY(MAXN), BOX
!!$         INTEGER     N
!!$
!!$         CHARACTER   CNFILE*(*)
!!$
!!$         INTEGER     CNUNIT, I
!!$         PARAMETER ( CNUNIT = 10 )
!!$
!!$         !    *******************************************************************
!!$
!!$         OPEN ( UNIT = CNUNIT, FILE = CNFILE, &
!!$                  STATUS = 'OLD', FORM = 'UNFORMATTED' )
!!$
!!$         READ ( CNUNIT ) N, BOX
!!$         IF ( N .GT. MAXN ) STOP ' N TOO LARGE '
!!$         READ ( CNUNIT ) ( RX(I), I = 1, N ), ( RY(I), I = 1, N )
!!$
!!$         CLOSE ( UNIT = CNUNIT )
!!$
!!$         RETURN
!!$       END SUBROUTINE READCN


       SUBROUTINE WORK ( MAXCAN, MAXV, NN, NV, NE, RX, RY, RS, VERTS, &
            VX, VY, IV, JV , SUCCESS)

         !    *******************************************************************
         !    ** ROUTINE TO PERFORM VORONOI ANALYSIS                           **
         !    **                                                               **
         !    ** WE WORK INITIALLY ON DOUBLE THE CORRECT SCALE,                **
         !    ** I.E. THE EDGES OF THE POLYGON GO THROUGH THE POINTS.          **
         !    *******************************************************************

         INTEGER     MAXCAN, NN, MAXV, NV, NE
         INTEGER     VERTS(MAXCAN)
         REAL(double)        RX(MAXCAN), RY(MAXCAN), RS(MAXCAN)
         REAL(double)        VX(MAXV), VY(MAXV)
         INTEGER     IV(MAXV), JV(MAXV)

         LOGICAL     OK
         INTEGER     I, J, L, NN1, N, V
         REAL(double)        AI, BI, CI, AJ, BJ, CJ, DET, DETINV
         REAL(double)        VXIJ, VYIJ
         REAL(double)        TOL
         logical :: success
         PARAMETER ( TOL = 1.d-10 )

         !    *******************************************************************

         !    ** IF THERE ARE LESS THAN 3 POINTS GIVEN **
         !    ** WE CANNOT CONSTRUCT A POLYGON         **

         success = .true.

         IF ( NN .LT. 3 ) THEN

            WRITE(*,'('' LESS THAN 3 POINTS GIVEN TO WORK '',I5)') NN
            Stop

         ENDIF

         NN1 = NN - 1
         V = 0

         !    ** WE AIM TO EXAMINE EACH POSSIBLE VERTEX  **
         !    ** DEFINED BY THE INTERSECTION OF 2 EDGES  **
         !    ** EACH EDGE IS DEFINED BY RX,RY,RS.       **

         DO I = 1, NN1

            AI =  RX(I)
            BI =  RY(I)
            CI = -RS(I)

            DO  J = I + 1, NN

               AJ =  RX(J)
               BJ =  RY(J)
               CJ = -RS(J)

               DET = AI * BJ - AJ * BI

               IF ( ABS ( DET ) .GT. TOL ) THEN

                  !             ** THE EDGES INTERSECT **

                  DETINV = 1.d0 / DET

                  VXIJ = ( BI * CJ - BJ * CI ) * DETINV
                  VYIJ = ( AJ * CI - AI * CJ ) * DETINV

                  !             ** NOW WE TAKE SHOTS AT THE VERTEX **
                  !             ** USING THE REMAINING EDGES ..... **

                  OK = .TRUE.
                  L  = 1

100               IF ( OK .AND. ( L .LE. NN ) ) THEN

                     IF ( ( L .NE. I ) .AND. ( L .NE. J ) ) THEN

                        OK = ( RX(L) * VXIJ + RY(L) * VYIJ ) .LE. RS(L)

                     ENDIF

                     L = L + 1
                     GOTO 100

                  ENDIF

                  !             ** IF THE VERTEX MADE IT      **
                  !             ** ADD IT TO THE HALL OF FAME **
                  !             ** CONVERT TO CORRECT SCALE   **

                  IF ( OK ) THEN

                     V = V + 1
                     IF ( V .GT. MAXV ) THEN
                        WRITE(*,*) 'TOO MANY VERTICES'
                        STOP
                     ENDIF
                     IV(V)  = I
                     JV(V)  = J
                     VX(V) = 0.5 * VXIJ
                     VY(V) = 0.5 * VYIJ

                  ENDIF

               ENDIF

            enddo

         enddo

         !    ** THE SURVIVING VERTICES DEFINE THE VORONOI POLYGON **

         NV = V

         IF ( NV .LT. 3 ) THEN

            WRITE(*,'('' LESS THAN 3 VERTICES FOUND IN WORK '',I5)') NV
!            do  i = 1, nn
!               write(*,*) i, rx(i), ry(i)
!            enddo

            SUCCESS = .false.
            goto 666

         ENDIF

         !    ** IDENTIFY NEIGHBOURING POINTS **

         DO N = 1, NN

            VERTS(N) = 0

         enddo

         DO V = 1, NV

            VERTS(IV(V)) = VERTS(IV(V)) + 1
            VERTS(JV(V)) = VERTS(JV(V)) + 1

         enddo

         !    ** POINTS WITH NONZERO VERTS ARE NEIGHBOURS **
         !    ** IF NONZERO, VERTS SHOULD BE EQUAL TO 2   **

         !    ** CHECK RESULT AND COUNT EDGES **

         OK = .TRUE.
         NE = 0

         DO N = 1, NN

            IF ( VERTS(N) .GT. 0 ) THEN

               NE = NE + 1

               IF ( VERTS(N) .NE. 2 ) THEN

                  OK = .FALSE.

               ENDIF

            ENDIF

         enddo

         IF ( .NOT. OK ) THEN

!            WRITE (*,'('' **** VERTEX ERROR: DEGENERACY ? **** '')')

         ENDIF

         IF ( NE .NE. NV ) THEN

!            WRITE(*,'('' **** EDGE   ERROR: DEGENERACY ? ****  '')')

         ENDIF

666 continue
       end SUBROUTINE WORK


       SUBROUTINE SORT_VORON ( MAXCAN, RX, RY, RS, TAG, NN )

         !    *******************************************************************
         !    ** ROUTINE TO SORT NEIGHBOURS INTO INCREASING ORDER OF DISTANCE  **
         !    **                                                               **
         !    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL.    **
         !    *******************************************************************

         INTEGER MAXCAN, NN
         REAL(double)    RX(MAXCAN), RY(MAXCAN), RS(MAXCAN)
         INTEGER TAG(MAXCAN)

         LOGICAL CHANGE
         INTEGER I, ITOP, I1, TAGI
         REAL(double)    RXI, RYI, RSI

         !    *******************************************************************

         CHANGE = .TRUE.
         ITOP = NN - 1

1000     IF ( CHANGE .AND. ( ITOP .GE. 1 ) ) THEN

            CHANGE = .FALSE.

            DO  I = 1, ITOP

               I1 = I + 1

               IF ( RS(I) .GT. RS(I1) ) THEN

                  RXI = RX(I)
                  RYI = RY(I)
                  RSI = RS(I)
                  TAGI = TAG(I)

                  RX(I) = RX(I1)
                  RY(I) = RY(I1)
                  RS(I) = RS(I1)
                  TAG(I) = TAG(I1)

                  RX(I1) = RXI
                  RY(I1) = RYI
                  RS(I1) = RSI
                  TAG(I1) = TAGI

                  CHANGE = .TRUE.

               ENDIF

            enddo

            ITOP = ITOP - 1
            GOTO 1000

         ENDIF

       END SUBROUTINE SORT_VORON

       subroutine sortConvex(n, xp, yp, xc, yc)
         implicit none
         integer n
         real(double) xp(*), yp(*)
         real(double) xc(*), yc(*)
         real(double) cosTheta(100)
         real(double) xs, ys, ts
         logical swap
         integer i, ip
         ip = 1
         do i = 2, n
            if ((yp(i) .lt. yp(ip))) then
               ip = i
            endif
         enddo
         do i = 1, n
            cosTheta(i) = atan2((yp(i)-yp(ip)) ,(xp(i)-xp(ip)))
            xc(i) = xp(i)
            yc(i) = yp(i)
         enddo



10       continue
         swap = .false.
         do i = 2, n
            if (cosTheta(i-1) .gt. cosTheta(i)) then
               swap = .true.
               xs = xc(i)
               ys = yc(i)
               ts = cosTheta(i)
               xc(i) = xc(i-1)
               yc(i) = yc(i-1)
               cosTheta(i) = cosTheta(i-1)
               xc(i-1) = xs
               yc(i-1) = ys
               cosTheta(i-1) = ts
            endif
         enddo
         if (swap) goto 10
       end subroutine sortConvex

       subroutine areaPolygon(n,x,y,area)
         implicit none
         integer  n, i, nt
         real(double) x(*), y(*), xt(1000), yt(1000)
         real(double) area, sum1, sum2
         do i = 1, n
            xt(i) = x(i)
            yt(i) = y(i)
         enddo
         nt = n + 1
         xt(nt) = x(1)
         yt(nt) = y(1)

         sum1 = 0.
         sum2 = 0.
         do i = 1, n
            sum1 = sum1 + xt(i) * yt(i+1)
            sum2 = sum2 + yt(i) * xt(i+1)
         enddo
         area = 0.5 * (sum1 - sum2)
       end subroutine areaPolygon

  subroutine bondi(x, y, z)
    real(double) :: x, y, y1, y2, z
    real(double), parameter :: lambda = 0.25d0*exp(1.d0)**(3.d0/2.d0)

    if (x > 0.5) then
       y1 = 0.01d0
    else
       y1 = 1.1d0
    endif

    y2 = 1.d30
    do while (abs((y2-y1)/y1) > 1.d-6)
       y2 = y1
       if (x > 0.5) then
          y1 = exp(-(-log(lambda)+2.d0*log(x) + 1.d0/x - 0.5d0*y1**2))
       else
          y1 = sqrt(2.d0*(-log(lambda)+2.d0*log(x) + 1.d0/x + log(y1)))
       endif
    enddo
    y = y1
    z = lambda / (x**2 * y1)
  end subroutine bondi

  function hfunc(y, x, lambda)
    real(double) :: x, lambda, hfunc, y

    hfunc = 0.5d0*y**2 - log(y) + log(lambda) - (1.d0/x + 2.d0*log(x))
  end function hfunc


  real(double) function macLaurinPhi(x, y, z, a1, a2, a3) result(phi)
    real(double) :: x, y, z, a1, a2, a3
    real(double) :: lam1, lam2, lam_mid, y1, y2, ym, z1
    real(double) :: bigA1, bigA3, h
    real(double), parameter :: rho = 1.d0
    logical :: converged, inside
    real(double) :: e
    e = sqrt(1.d0-(a3/a1)**2)


    inside = .false.
    if (sqrt(x**2 + y**2) < a1) then
       z1 = sqrt(a3**2 * (1.d0 - (x**2 + y**2)/a1**2))
       if (abs(z) <= z1) then
          inside = .true.
       endif
    endif


    if (inside) then
          bigA1 = sqrt(1.d0-e**2)/e**3 * asin(e) - (1.d0-e**2)/e**2
          bigA3 = (2.d0/e**2) - (2.d0*sqrt(1.d0-e**2)/e**3) * asin(e)
          phi = -pi * bigG * rho * (2.d0*bigA1 * a1**2 - bigA1 * (x**2+y**2) + bigA3 * (a3**2 - z**2))
       else

          lam1 = 0.d0
          lam2 = 1.d25
          converged = .false.
          do while(.not.converged)
             lam_mid = 0.5d0*(lam1 + lam2)
             y1 = lambdaFunc(x, y, z, a1, a2, a3, lam1)
             y2 = lambdaFunc(x, y, z, a1, a2, a3, lam2)
             ym = lambdaFunc(x, y, z, a1, a2, a3, lam_mid)


             if (y1*ym < 0.d0) then
                lam1 = lam1
                lam2 = lam_mid
             else if (y2*ym < 0.d0) then
                lam1 = lam_mid
                lam2 = lam2
             else
                converged = .true.
                lam_mid = 0.5d0*(lam1+lam2)
             endif

             if (abs((lam1-lam2)/lam2) .le. 1.d-6) then
                converged = .true.
             endif
          enddo

          h = a1*e /(sqrt(a3**2 + lam_mid))

          phi = -(2.d0*pi*bigG*rho*a1*a3/e) * (atan(h) - (1.d0/(2.d0*a1**2*e**2)) * ( (x**2+y**2)*(atan(h) - (h/(1.d0+h**2))) + &
               2.d0 * z**2 * (h - atan(h)) ) )

       endif

     end function macLaurinPhi

  real(double) function lambdaFunc(x, y, z, a1, a2, a3, lambda)
    real(double) :: x, y, z, a1, a2, a3, lambda

    lambdaFunc = x**2/(a1**2 + lambda) + y**2 / (a2**2 + lambda) + z**2 / (a3**2 + lambda) - 1.d0
  end function lambdaFunc

  subroutine myInterp(p,nPoints, x, y, z, a, x0, y0, z0, a0)
    integer :: nPoints
    real(kind=double) :: x(:), y(:),z(:),a(:)
    integer :: p
    real(kind=double) :: x0, y0, z0, a0
    real(kind=double), allocatable :: r(:), w(:)

    allocate(r(1:nPoints), w(1:npoints))

    r(1:nPoints) = sqrt((x(1:npoints)-x0)**2 + (y(1:npoints)-y0)**2 + (z(1:npoints)-z0)**2)
    w(1:nPoints) = (1.d0/r(1:nPoints))**p
    w(1:npoints) = w(1:nPoints)/SUM(w(1:nPoints))

    a0 = SUM(w(1:npoints)*a(1:nPoints))
    deallocate(r, w)
  end subroutine myInterp


!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+
!     File:        ccubsolv.f (Fortran 77 source for the ccubsolv routine)
!     Author:      Fredrik Jonsson <fj@phys.soton.ac.uk>
!     Date:        December 29, 2005
!     Last change: December 29, 2006
!     Description: The CCUBSOLV(A,Z) routine solves the cubic polynomial
!
!                        z^3+c[2]*z^2+c[1]*z+c[0]=0
!
!                  for general complex coefficients c[k]. The routine
!                  takes a vector A(1..3) of COMPLEX*8 floating point
!                  precision as input, containing the c[k] coefficients
!                  as conforming to the convention
!
!                        A(1)=c[0],  A(2)=c[1],  A(3)=c[2],
!
!                  and returns the three complex-valued roots for z in the
!                  vector Z(1..3), also being of COMPLEX*8 floating point
!                  precision.
!
!     Copyright (C) 2006, Fredrik Jonsson <fj@phys.soton.ac.uk>
!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+
  subroutine ccubsolv(a,z)
    complex(double) a(3),z(3)
    real(double) pi,phi,onethird,twothird,absrr,absqq,absrr23,absqq32
    complex(double) p,q,rr,qq,cc,w
    integer k

    pi=3.14159265358979323846d0
    onethird=1.0/3.0
    twothird=2.0/3.0
    p=a(2)-onethird*(a(3)**2)
    q=2.0*(onethird*a(3))**3-onethird*a(2)*a(3)+a(1)
    rr=-0.5*q
    qq=onethird*p

    !     *********************************************************************
    !     The following construction is made in order to avoid numeric overflow
    !     in evaluation of the discriminant, in a a manner similar to that used
    !     in standard routines for evaluation of sqrt(x^2+y^2), for example.
    !     *********************************************************************
    absrr=abs(rr)
    absqq=abs(qq)
    absrr23=absrr**twothird
    if ((absrr23).lt.(absqq)) then
       absqq32=absqq**1.5
       cc=rr+absqq32*sqrt((rr/absqq32)**2+(qq/absqq)**3)
    else
       cc=rr+absrr*sqrt((rr/absrr)**2+(qq/absrr23)**3)
    endif

    !     *********************************************************************
    !     Evaluate the solutions for w from the binomial equation w^3=c, and
    !     form the solutions z(k) from the three obtained roots.
    !     *********************************************************************
    cc=cc**onethird
    do k=0,2
       phi=k*twothird*pi
       w=cc*cmplx(cos(phi),sin(phi), kind=double)
       z(k+1)=w-onethird*(a(3)+p/w)
    enddo

  end subroutine ccubsolv

  complex(double) function Ylm(l, m, theta, phi)
    real(double) :: theta, phi
    integer :: l, m
    ylm = sqrt( (dble((2*l + 1)*factorial(l-m))/(fourPi * dble(factorial(l+m))) )) * &
         plgndr(l,m,cos(theta))*cmplx(cos(dble(m)*phi), sin(dble(m)*phi),kind=double)
  end function Ylm
  integer function factorial(n)
    integer :: i, n
    factorial = 1
    do i = 1, n
       factorial = factorial * i
    enddo
  end function factorial

  real(double) FUNCTION PLGNDR(L,M,X)
    integer :: l, m,i, ll
    real(double) :: x, pmm, fact, somx2, pmmp1,  pll
    IF (M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.) call writeFatal('bad arguments')
    PMM=1.
    IF(M.GT.0) THEN
       SOMX2=SQRT((1.-X)*(1.+X))
       FACT=1.
       DO I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.
       enddo
    ENDIF
    IF(L.EQ.M) THEN
       PLGNDR=PMM
    ELSE
       PMMP1=X*(2*M+1)*PMM
       IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
       ELSE
          DO LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
         enddo
         PLGNDR=PLL
      ENDIF
   ENDIF
 END FUNCTION PLGNDR

 ! returns a random value from a Gaussian PDF between +/- 5*sigma of the mean
 function randomValueGaussian(sigma, xMean) result (xOut)
   real(double) :: xMin, xMax, xOut, sigma, xMean
   integer, parameter :: nx = 100
   integer :: i
   real(double) :: xArray(nx), prob(nx), r, t

   xMin = xMean - 5.d0*sigma
   xMax = xMean + 5.d0*sigma

   do i = 1, nx
      xArray(i) = xMin + (dble(i-1)/dble(nx-1)) * (xMax-xMin)
   enddo

   ! gaussian pdf
   prob = 1.d0/(sigma*sqrt(twoPi)) * exp(-(xArray(:)-xMean)**2/(2.d0*sigma**2))
   do i = 2, nx
      prob(i) = prob(i) + prob(i-1)
   enddo
   prob(1:nx) = prob(1:nx) - prob(1)
   prob(1:nx) = prob(1:nx) / prob(nx)
   call randomNumberGenerator(getDouble=r)
   call locate(prob, nx, r, i)
   t = (r - prob(i))/(prob(i+1)-prob(i))
   xOut = xArray(i) + t * (xArray(i+1)-xArray(i))
 end function randomValueGaussian

end module utils_mod
