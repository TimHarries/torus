! written by tjh


module utils_mod

  use kind_mod
  use vector_mod          ! vector maths
  use constants_mod, only: angstromToCm, cSpeed, kErg, hCgs, ergtoEv, mHydrogen, twoPi, pi, mElectron
  use messages_mod
  use nrtype
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


contains

  subroutine testFileExists(thisFile)
    character(len=*) :: thisFile
    integer :: error
    open(53, file=thisFile, status="old",form="formatted", iostat=error)

    if (error /= 0) then
       call writeFatal("Error opening file: "//trim(thisFile))
       stop
    else
       close(53)
    endif
  end subroutine testFileExists

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
    endif

  end subroutine solveQuadDble

!  subroutine solveQuadDble(a, b, c, x1, x2,ok)
!    implicit none
!    real(double),intent(in) :: a, b, c
!    real(double),intent(out):: x1, x2
!    real(double)            :: y
!    logical,intent(out)              :: ok
    
!    real(double) :: OneOverTwoA, sqrty, MinusBOverTwoA

!    ok = .true.

!    ! special case:
!    if (a == 0) then
!       if (b/=0) then
!          x1 = -c/b; x2 = -c/b 
!       else
!          write(*,*) "! quad solver failed (1):   (a,b,c) = ", a, b, c
!          ok = .false.
!          x1 = 0.0d0; x2 =0.0d0 
!       end if
!    end if

!    OneOverTwoA = 1.d0 / (2.d0 * a)
!    MinusBOverTwoA = -b * OneOverTwoA

!    y = b*b-4.d0*a*c
!    sqrty = y**0.5d0
!    sqrty = sqrty * OneOverTwoA

!    if (y >= 0.) then
!       x1 = MinusBOverTwoA - sqrty !(-b + sqrt(y))/(2.d0*a)
!       x2 = x1 + (2.d0 * sqrty)
!    else
!       write(*,*) "! quad solver failed (2).:  (y,a,b,c) = ", y,a,b,c
!       ok = .false.
!       x1 = -b/(2.d0*a)
!       x2 = -b/(2.d0*a)
!    endif
!
!  end subroutine solveQuadDble

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
          offset = 6.7d-9
       case("V") 
          offset = 3.75d-9
       case("R")
          offset = 1.8e-9
       case DEFAULT
          call writeWarning("Band "//trim(band)//" not recognised in returnMagntiude")
          offset = 0.d0
    end select
    out = -2.5d0 * log10(flux) + 2.5d0*log10(offset)
  end function returnMagnitude

  real elemental function blackBody(temperature, wavelength)

    real,intent(in) :: temperature
    real,intent(in) :: wavelength
    real            :: nu, fac

    nu = cSpeed / (wavelength*angstromToCm)
    fac= log(2.)+log(hCgs) + 3.*log(nu) - 2.*log(cSpeed)
    if (hCgs*nu/(kerg*temperature) < 33.) then
       blackBody = exp(fac)/(exp(hCgs*nu/(kErg*temperature))-1.)
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
    GAMMLN=TMP+LOG(STP*SER)
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
    ANS=LY1+GR*(LX-LX1)
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
    LOGINT_dble=EXP(ANS)
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
    USE nrtype; USE nrutil, ONLY : nrerror
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
          if (jstack > NSTACK) call nrerror('sort: NSTACK too small')
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

  SUBROUTINE dquicksort2(arr,slave)
    USE nrtype; USE nrutil, ONLY : assert_eq
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr,slave
    INTEGER :: ndum
    INTEGER, DIMENSION(size(arr)) :: index
    ndum=assert_eq(size(arr),size(slave),'sort2')
    call indexx(size(arr), arr,index)
    arr=arr(index)
    slave=slave(index)
  END SUBROUTINE dquicksort2

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

  SUBROUTINE polint3(xa,ya,n,x,y,dy)
    INTEGER n,NMAX
    REAL :: x,xa(n), ya(n)
    real(double) :: y, dy
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
             print *, 'Error:: failure in polint3.'
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
  end SUBROUTINE polint3




  ! spline interpolation

  SUBROUTINE SPLINE_SINGLE(X,Y,N,YP1,YPN,Y2)
    INTEGER NMAX,N, I, K
    REAL SIG, P, QN, UN
    PARAMETER (NMAX=1000)
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
    integer, parameter :: nmax=500
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



  pure subroutine getDerivs(x, y, n, derivs)
    implicit none
    integer,intent(in) :: n
    real,intent(in)    :: x(:), y(:)           
    real,intent(out)   :: derivs(:)
    integer            :: i

    do i = 2, n-1
       derivs(i) = log(y(i-1)/y(i+1)) / log(x(i-1)/x(i+1)) * y(i)/x(i)
    enddo

    derivs(1) = log(y(1)/y(2)) / log(y(1)/y(2))
    if ( (x(1)-x(3))/(x(1)-x(2)) > 3.) derivs(2) = derivs(1)*y(2)/x(2)

    derivs(1) = derivs(1) * y(1) / x(1)

    derivs(n) = log(y(n-1)/y(n)) / log(x(n-1)/x(n))
    if ( (x(n-2)-x(n))/(x(n-1)-x(1)) > 3.) &
         derivs(n-1)  = derivs(n)*y(n-1)/x(n-1)
    derivs(n) = derivs(n)*y(n)/x(n)
  end subroutine getDerivs

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

    real(double) function GaussianComplementRand()

      real(double) :: x,y, ytest
      logical :: done

      done = .false.

      do while(.not. done)
         call randomNumberGenerator(getDouble=x)
         x = 6. * (x - 0.5d0) 
         y = 1. - 0.5 * unNormalizedGauss(1.d0,x)
         
         call randomNumberGenerator(getDouble=ytest)         
         if(ytest .lt. y) then
            GaussianComplementRand = x
            done = .true.
         else
            done = .false.
         endif
      end do

    end function GaussianComplementRand

    pure real function gauss(sigma, dx)
      real, intent(in) :: sigma, dx
      real :: fac
      fac = 1./(sigma*sqrt(twoPi))
      gauss = fac * exp(-(dx**2)/(2.*sigma**2))
    end function gauss

    pure real(double) function unNormalizedGauss(sigma, dx)
      real(double), intent(in) :: sigma, dx
      unNormalizedGauss = exp(-(dx**2)/(2.*sigma**2))
    end function unNormalizedGauss

    subroutine timeString(iTime,time)
      character(len=*) :: time
      integer :: iTime

      integer :: iHour, iMin, iSec, ir

      ir = iTime

      iHour = int(real(ir)/3600.)
      ir = ir - iHour*3600
      iMin = int(real(ir)/60.)
      ir = ir - iMin * 60
      iSec = ir

      write(time,'(i2,a,i2.2,a,i2.2)') iHour,":",iMin,":",iSec
    end subroutine timeString      


    subroutine replaceDots(cString, done)
      character(len=*) :: cString
      integer :: i
      logical :: done

      done = .false.
      i = index(cString, ".")
      do while (i /= 0)
         done = .true.
         cString(i:i) = "_"
         i = index(cString, ".")
      end do
    end subroutine replaceDots

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

      fac = (4./sqrt(pi))*(mass/(2.*kErg*temperature))**1.5
      fac2 = v**2 * exp(-mass*v**2 / (2.*Kerg*temperature))
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
    






    ! a functions to convert characters to integer 
    function char2int(a) RESULT(out)
      implicit none 
      integer :: out
      character, intent(in),dimension(:) :: a
      integer :: ierr
      
      read(a, *, IOSTAT = ierr) out
      if (ierr /= 0) then
         print *, 'Error : Non-numerical characters passed to [char2int] &
              & function in char_function_class module.'
         stop
      end if
      
    end function char2int
  
  
    ! a function to convert integer to strings
    function int2char(num) RESULT(out)
      implicit none
      character(LEN=50) :: out
      integer, intent(in) :: num
      
      write(out,*) num
    end function int2char
    
    

    ! A function to attach numbers at the end of strings
    ! The length of output character is 50.
    function tail_num_to_char(strings, number) RESULT(out)
      character(LEN=50) :: out
      character(LEN=*), intent(in) :: strings
      integer, intent(in) :: number
      character(LEN=50) :: char_number
      integer :: length, n
      
      ! convert number into strings .. using a function in this module
      char_number = int2char(number)
      char_number = ADJUSTL(char_number)
      length = LEN( TRIM(char_number) )
      
      ! Stick them togeter.
      n = LEN_TRIM(strings)
      out = strings(1:n)//char_number(1:length)
      
    end function tail_num_to_char


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
!      IF(A .EQ. 0.0D0  ) THEN
!         ! ---- Normal Doppler brodening.         
!         IF (V2 < 1.D-04) THEN
!            VOIGTN = 1.0D0 + V2*(1.0D0 + 0.5D0*V2*(1.0D0 + V2/6.0D0))
!         ELSE
!            VOIGTN = EXP(-V2)
!         END IF
!         RETURN
!      END IF
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

    real function bigLambda(N_HII, temperature, Ne)
      real :: N_HII
      real :: temperature
      real :: Ne

      real, parameter :: C_rad  = 1.
      real, parameter :: C_vdw = 0.
      real, parameter :: C_stark =1.

! see Muzerolle et al. 2001 ApJ 550 944

      bigLambda = C_rad + C_vdw*(N_HII / 1.e16)*(temperature/5000.)*0.3 + C_stark*(Ne/1.e12)**0.6666
    end function bigLambda

    !
    ! Damping contant in a Voigt Profile in [1/s]
    !
    real function bigGamma(N_HI, temperature, Ne, nu)      
      use input_variables,  only:  C_rad, C_vdw, C_stark
      real(double), intent(in) :: N_HI         ! [#/cm^3]  number density of HI
      real(double), intent(in) :: temperature  ! [Kelvins]
      real(double), intent(in) :: Ne           ! [#/cm^3]  nunmber density of electron
      real(double), intent(in) :: nu           ! [1/s]  line center frequency 

      !------------------------------------------------------------
      ! Quadratic Stark broadening (?): 
      ! -- Good for most lines especially hot stars (Gray's comment)
      !------------------------------------------------------------
      !! see Muzerolle et al. 2001 ApJ 550 944
      !bigGamma = &
      !     &    C_rad  &
      !     &  + C_vdw*(N_HI / 1.e16)*(temperature/5000.)**0.3 &
      !     &  + C_stark*(Ne/1.e12)**0.6666  ! [Angstrom]


      !----------------------------------------------------------
      ! Linear Stark broadening:
      ! -- Good for Hydrogen lines
      !---------------------------------------------------------
      ! See Luttermoser & Johnson, 1992, ApJ, 388, 579
      !  Note: smallGamma = bigGamma/(4*pi)
      !        e.g. For H-alpha: C_rad   = 8.16e-3 [A], 
      !                          C_vdw   = 5.53e-3 [A], 
      !                          C_stark = 1.47e-2 [A], 

      bigGamma = &
           &    C_rad  &
           &  + C_vdw*(N_HI / 1.e16)*(temperature/5000.)**0.3 &
           &  + C_stark*(Ne/1.e12)  ! [Angstrom]



      ! convert units 
      bigGamma = (bigGamma*1.e-8) * nu**2   / cSpeed  ! [1/s]
      !                  [cm]       * [1/s^2] / [cm/s]
    end function bigGamma


    
    subroutine resampleRay(lambda, nTau, projVel, maxtau, newLambda, newNTau, &
         inFlow, newInFlow)
      use  input_variables, only: lamstart, lamend, nlambda, lamline
      integer, intent(in) :: nTau, maxtau
      real, intent(in) :: lambda(nTau)
      real(double), intent(in) :: projVel(nTau)
      integer, intent(inout) :: newNtau
      real, intent(inout)  :: newLambda(maxtau)
      logical, intent(inout)  :: InFlow(nTau)
      logical, intent(inout)  :: newInFlow(maxtau)
      !
      real(double) :: dProjVel  ! automatic array
      integer :: nAdd 
      integer :: i, j
!      real, parameter :: dvel  = 10.e5/cSpeed
      real :: dvel, dlam

      ! -- using the values in input_variables module
      dvel = (lamend-lamstart)/lamline/real(nlambda-1)  ! should be in [c]
      dvel = dvel/2.0  ! to be safe  
      dvel = 1.e5/cspeed ! 1 km/s

      newNTau = 0
      do i = 2, nTau
         dProjVel = projVel(i) - projVel(i-1)
         ! if (i==2) we always add points
         if (abs(dProjVel) > dVel .or. i==2) then
            nAdd = nint(abs(dProjVel)/dVel)
            dlam = (lambda(i)-lambda(i-1))/real(nAdd+1)
            do j = 1, nAdd+1
               newNtau = newNtau + 1
               newLambda(newNTau) = real(j-1)*dlam + lambda(i-1)
               newInFlow(newNTau) = inFlow(i-1)
            enddo
         else
            newNtau = newNtau + 1
            newLambda(newNTau) = lambda(i-1)
            newInFlow(newNTau) = inFlow(i-1)
         endif
      enddo
      newNtau = newNtau + 1
      newLambda(newNTau) = lambda(nTau)
      newInFlow(newNTau) = inFlow(nTau)
    end subroutine resampleRay


    !
    ! Adding extra points near the resonace zone if any. 
    !
    subroutine resampleRay2(lambda, nTau, projVel, thisVel, maxTau, &
         newLambda, newNTau, inFlow, newInFlow)
      integer, intent(in) :: nTau
      integer, intent(in) :: maxtau
      real, intent(in) :: lambda(nTau)  ! ray sequments, but not the wavelenghth!
      real(double), intent(in) :: projVel(nTau)
      real, intent(in) :: thisVel       ! velocity at the photon emission 
      integer, intent(inout) :: newNtau
      real, intent(inout) :: newLambda(maxTau)
      logical, intent(inout)  :: InFlow(nTau)
      logical, intent(inout)  :: newInFlow(maxtau)
      integer,parameter :: nAdd = 50
      integer,parameter :: nclose =10
!      integer, parameter :: nAdd = 10
      integer :: i, j
      real :: dlam
      logical :: in_res_zone
      
      ! loop through the projected velocities to find the resonance zones
      newNTau = 0

      ! check for the resonance zone(s)
      do i = 2, nTau 
         if ( ( projVel(i-1) <= thisVel  .and. thisVel < projVel(i) ) &
              .or. &
              ( projVel(i-1) >= thisVel  .and. thisVel > projVel(i) )  )  then
            in_res_zone = .true.
         else
            in_res_zone = .false.
         end if

         if (in_res_zone .or. i < nclose ) then
            ! adding extra points (linearly interporating between points).
            dlam = (lambda(i)-lambda(i-1))/real(nAdd+1)
            do j = 1, nAdd+1
               newNtau = newNtau + 1
               newLambda(newNTau) = real(j-1)*dlam + lambda(i-1)
               newInFlow(newNTau) = inFlow(i-1)
            enddo
         else
            newNtau = newNtau + 1
            newLambda(newNTau) = lambda(i-1)
            newInFlow(newNTau) = inFlow(i-1)
         endif
      end do
      ! last points
      newNtau = newNtau + 1
      newLambda(newNTau) = lambda(nTau)
      newInFlow(newNTau) = inFlow(nTau)
    end subroutine resampleRay2

    !
    ! Add points where velocity is slowly changing.
    ! NEED TO ADD ALPHA_DISC CHECK SECTION LATER! (IMPORTNAT IF DISC IS PRESENT)
    subroutine resampleRay3(lambda, nTau, projVel, maxtau, newLambda, newNTau, &
         inFlow, newInFlow)
      integer, intent(in) :: nTau, maxtau
      real, intent(in) :: lambda(nTau)
      real(double), intent(in) :: projVel(nTau)
      integer, intent(inout) :: newNtau
      real, intent(inout)  :: newLambda(maxtau)
      logical, intent(inout)  :: InFlow(nTau)
      logical, intent(inout)  :: newInFlow(maxtau)
      !
      real(double) :: dProjVel  ! automatic array
      integer,parameter :: nAdd = 5
      integer :: i, j
      real :: dvel, dlam

      dvel = 10.e5/cspeed !10 km/s

      newNTau = 0
      do i = 2, nTau
         dProjVel = projVel(i) - projVel(i-1)
         if (abs(dProjVel) > dVel) then
            dlam = (lambda(i)-lambda(i-1))/real(nAdd+1)
            do j = 1, nAdd+1
               newNtau = newNtau + 1
               newLambda(newNTau) = real(j-1)*dlam + lambda(i-1)
               newInFlow(newNTau) = inFlow(i-1)
            enddo
         else
            newNtau = newNtau + 1
            newLambda(newNTau) = lambda(i-1)
            newInFlow(newNTau) = inFlow(i-1)
         endif
      enddo
      newNtau = newNtau + 1
      newLambda(newNTau) = lambda(nTau)
      newInFlow(newNTau) = inFlow(nTau)
    end subroutine resampleRay3

    
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



    !
    ! Inserting additional points in ray (lambda) array.
    ! 
    subroutine insert_points_in_ray(lambda, nTau, nadd, newLambda, newNTau)
      real, intent(in)     :: lambda(:)  ! ray segments
      integer, intent(in)  :: nTau       ! number of points along ray
      integer, intent(in)  :: nadd       ! number of points to be inserted in between points
      integer, intent(out) :: newNtau      ! new number of points along ray
      real, intent(out)    :: newLambda(:) ! new ray segments.
      integer :: i, j
      real :: dlam

      newNTau = 0
      do i = 2, nTau
         dlam= (lambda(i)-lambda(i-1))/real(nAdd+1)
         do j = 1, nAdd+1
            newNtau = newNtau + 1
            newLambda(newNTau) = lambda(i-1) +  real(j-1)*dlam
         end do
      end do
      newNtau = newNtau + 1
      newLambda(newNTau) = lambda(nTau)
    end subroutine insert_points_in_ray


               
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
            newYarray(i) = yArray(j) + (yArray(j+1)-yArray(j))*(newXarray(i)-xArray(j))/(xArray(j+1)-xArray(j))
         else
            newYarray(i) = yArray(j)
         endif
      enddo
      yArray(1:newNx) = newYArray(1:newNx)
    end subroutine linearResample_dble




    subroutine log_log_resample(xArray, yArray, nX, newXarray, newNx)
      implicit none
      real, intent(in)    :: xArray(:)
      real, intent(inout) :: yArray(:)
      integer, intent(in) :: nx, newNx
      real, intent(in)    :: newXarray(:)
      !
      real :: newYarray(newNx) ! automatic array
      integer :: i, j
      real:: log_y2, log_y1, log_x2, log_x1, T1

      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         if (j >= Nx) j = j -1
         log_x1 = LOG(xArray(j))
         log_x2 = LOG(xArray(j+1))
         if (log_x1 /= log_x2) then
            log_y1 = LOG(yArray(j))
            log_y2 = LOG(yArray(j+1))
            T1 = (log_y2-log_y1)/(log_x2-log_x1)
            newYarray(i) = log_y1 + T1*(LOG(newXarray(i))-log_x1)
            newYarray(i) = EXP(newYarray(i))
         else
            newYarray(i) = yArray(j)
         endif
      enddo
      yArray(1:newNx) = newYArray(1:newNx)
    end subroutine log_log_resample



    subroutine log_log_resample_dble(xArray, yArray, nX, newXarray, newNx)
      implicit none      
      real, intent(in)            :: xArray(:)
      real(double), intent(inout) :: yArray(:)
      integer, intent(in)         :: nx, newNx
      real, intent(in)            :: newXarray(:)
      !
      real(double) :: newYarray(newNx) ! automatic array
      integer      :: i, j
      real(double) :: log_y2, log_y1, log_x2, log_x1, T1

      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         if (j >= Nx) j = j -1
         log_x1 = LOG(xArray(j)); log_x2 = LOG(xArray(j+1))
         if (log_x1 /= log_x2) then
            log_y1 = LOG(yArray(j))
            log_y2 = LOG(yArray(j+1))
            
            T1 = (log_y2-log_y1)/(log_x2-log_x1)
            newYarray(i) = log_y1 + T1*(LOG(newXarray(i))-log_x1)
            newYarray(i) = EXP(newYarray(i))
         else
            newYarray(i) = yArray(j)
         endif
      enddo
      yArray(1:newNx) = newYArray(1:newNx)
    end subroutine log_log_resample_dble


    subroutine quadraticResample(xArray, yArray, nx, nx_max,newXarray, newNx)
      real :: xArray(:), yArray(:)
      integer :: nx, newNx
      integer, intent(in) :: nx_max
      real :: newXarray(:)
      integer :: i, j
      real :: y
      real :: dy  ! error estimate
      real, save, allocatable  :: newYarray(:) ! newYarray(newNx)  ! automatic array
      logical, save :: first_time = .true.
      !$OMP THREADPRIVATE (first_time, newYarray)
       y = 0.
       dy = 0.
      if (first_time) then
         first_time=.false.
         ALLOCATE(newYarray(nx_max))
      end if


    ! polynomial interpolation 
      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         if (j==1) then
            j=2
         elseif (j==nx) then
            j=nx-1
         end if
         call polint(xArray(j-1:j+1), yArray(j-1:j+1), 3, &
              newXarray(i), y, dy)
         newyarray(i) = y
      enddo

      yArray(1:newNx) = newYArray(1:newNx)
      
    end subroutine quadraticResample


    subroutine quadraticResample_dble(xArray, yArray, nx, nx_max,  newXarray, newNx)
      real, intent(in)            :: xArray(:)
      real(double), intent(inout) :: yArray(:)
      integer, intent(in)         :: nx, newNx
      integer, intent(in)         :: nx_max
      real, intent(in)            :: newXarray(:)
      !
      integer      :: i, j
      real(double) :: dy  ! error estimate
      real(double) :: y
      real(double), save, allocatable  :: newYarray(:) ! newYarray(newNx)  ! automatic array
      logical, save :: first_time = .true.
      !$OMP THREADPRIVATE (first_time, newyarray)

      dy = 0.; y = 0.
      if (first_time) then
         first_time=.false.
         ALLOCATE(newYarray(nx_max))
      end if


    ! polynomial interpolation 
      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         if (j<=1) then
            j=2
         elseif (j>=nx) then
            j=nx-1
         end if
         call polint2(xArray(j-1:j+1), yArray(j-1:j+1), 3, &
              newXarray(i), y, dy)
         newYarray(i) = y
      enddo

      yArray(1:newNx) = newYArray(1:newNx)
      
    end subroutine quadraticResample_dble


!already in nrutil.f90
 
!  FUNCTION arth(first,increment,n)
!    ! Array function returning an arithmetic progression. 
!    REAL, INTENT(IN) :: first,increment 
!    INTEGER, INTENT(IN) :: n 
!    REAL, DIMENSION(n) :: arth 
!    INTEGER :: k,k2 
!    REAL :: temp 
!    INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
!    if (n > 0) arth(1)=first 
!    if (n <= NPAR_ARTH) then 
!      do k=2,n
!        arth(k)=arth(k-1)+increment 
!      end do
!    else
!      do k=2,NPAR2_ARTH 
!        arth(k)=arth(k-1)+increment 
!      end do 
!      temp=increment*NPAR2_ARTH 
!      k=NPAR2_ARTH
!      do
!        if (k >= n) exit 
!        k2=k+k 
!        arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
!        temp=temp+temp 
!        k=k2
!      end do 
!    end if
!  END FUNCTION arth 
  
  SUBROUTINE trapzd(func,a,b,s,n) 
    use nrutil, only: arth
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
    real, allocatable :: xtemp(:)
    integer :: i, newNx

    allocate(xtemp(1:nx))
    
    call sort(nx, x)

    newnx = 1
    xtemp(newnx) = x(1)
    do i = 2, nx
       if (abs(xtemp(newnx)-x(i)) > xTol) then
          newnx = newnx  + 1
          xtemp(newnx) = x(i)
       endif
    enddo

    x(1:newnx) = xtemp(1:newnx)
    nx = newnx
    deallocate(xtemp)
  end subroutine stripSimilarValuesSingle

  subroutine stripSimilarValuesDouble(x, nx, xtol)
    integer :: nx
    real(double) :: x(:)
    real(double) :: xtol
    real(double), allocatable :: xtemp(:)
    integer :: i, newNx

    allocate(xtemp(1:nx))
    
    call sort(nx, x)

    newnx = 1
    xtemp(newnx) = x(1)
    do i = 2, nx
       if (abs(xtemp(newnx)-x(i)) > xTol) then
          newnx = newnx  + 1
          xtemp(newnx) = x(i)
       endif
    enddo

    x(1:newnx) = xtemp(1:newnx)
    nx = newnx
    deallocate(xtemp)
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
          freq = cSpeed / wavelength
          energy = freq * hcgs
          array(i) = energy * ergtoEv
       else
          array(i) = 1.e-20
       endif
    enddo
  end subroutine wavenumbertoEv

  !
  !
  ! A safer way to write message.
  subroutine write_message(filename, message)
    implicit none
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: message
    integer, parameter :: luout = 82
    logical, save :: first_time =  .true.
    !$OMP THREADPRIVATE (first_Time)
    !
    
    if (first_time) then
       open(unit=luout, file=filename, status="replace")
       first_time = .false.
    else
       open(unit=luout, file=filename, status="old", position="append")
    end if


    write(luout, "(a)") TRIM(ADJUSTL(message))

    close(luout)

  end subroutine write_message



  !
  !
  ! A safer way to write message with an integer value.
  subroutine write_message_int(filename, message, val)
    implicit none
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: message
    integer, intent(in) :: val
    integer, parameter :: luout = 82
    logical, save :: first_time =  .true.
    !$OMP THREADPRIVATE (first_Time)
    !
    
    if (first_time) then
       open(unit=luout, file=filename, status="replace")
       first_time = .false.
    else
       open(unit=luout, file=filename, status="old", position="append")
    end if


    write(luout, *) TRIM(ADJUSTL(message)), "==>", val

    close(luout)

  end subroutine write_message_int


  !
  !
  ! A safer way to write message with an real value.
  subroutine write_message_real(filename, message, val)
    implicit none
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: message
    real, intent(in) :: val
    integer, parameter :: luout = 82
    logical, save :: first_time =  .true.
    !$OMP THREADPRIVATE (first_Time)
    !
    
    if (first_time) then
       open(unit=luout, file=filename, status="replace")
       first_time = .false.
    else
       open(unit=luout, file=filename, status="old", position="append")
    end if


    write(luout, *) TRIM(ADJUSTL(message)), "==>", val

    close(luout)

  end subroutine write_message_real


  !
  !
  ! A safer way to write message with an real value.
  subroutine write_message_dble(filename, message, val)
    implicit none
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: message
    real(double), intent(in) :: val
    integer, parameter :: luout = 82
    logical, save :: first_time =  .true.
    !$OMP THREADPRIVATE (first_Time)

    !
    
    if (first_time) then
       open(unit=luout, file=filename, status="replace")
       first_time = .false.
    else
       open(unit=luout, file=filename, status="old", position="append")
    end if


    write(luout, *) TRIM(ADJUSTL(message)), "==>", val

    close(luout)

  end subroutine write_message_dble


  !
  !
  ! A safer way to write message with an real value.
  subroutine write_message_vec(filename, message, val)
    implicit none
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: message
    type(vector), intent(in) :: val
    integer, parameter :: luout = 82
    logical, save :: first_time =  .true.
    !$OMP THREADPRIVATE (first_Time)

    !
    
    if (first_time) then
       open(unit=luout, file=filename, status="replace")
       first_time = .false.
    else
       open(unit=luout, file=filename, status="old", position="append")
    end if


    write(luout, *) TRIM(ADJUSTL(message)), "==>", val

    close(luout)

  end subroutine write_message_vec


  !
  !
  ! A safer way to write message with an real value.
  subroutine write_message_logical(filename, message, val)
    implicit none
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: message
    logical, intent(in) :: val
    integer, parameter :: luout = 82
    logical, save :: first_time =  .true.
    !$OMP THREADPRIVATE (first_Time)

    !
    
    if (first_time) then
       open(unit=luout, file=filename, status="replace")
       first_time = .false.
    else
       open(unit=luout, file=filename, status="old", position="append")
    end if


    write(luout, *) TRIM(ADJUSTL(message)), "==>", val

    close(luout)

  end subroutine write_message_logical


  !
  !
  ! A safer way to write message with an real value.
  subroutine write_message_char(filename, message, val)
    implicit none
    character(LEN=*), intent(in) :: filename
    character(LEN=*), intent(in) :: message
    character(LEN=*), intent(in) :: val
    integer, parameter :: luout = 82
    logical, save :: first_time =  .true.
    !$OMP THREADPRIVATE (first_Time)

    !
    
    if (first_time) then
       open(unit=luout, file=filename, status="replace")
       first_time = .false.
    else
       open(unit=luout, file=filename, status="old", position="append")
    end if


    write(luout, *) TRIM(ADJUSTL(message)), "==>", TRIM(ADJUSTL(val))

    close(luout)

  end subroutine write_message_char

! Count the number of lines in a file. 
! Should return the same answer as wc -l i.e. includes blank lines in the total
! D. Acreman, June 2010
  integer function file_line_count(filename)

    character(len=*) :: filename
    character(len=1) :: dummy
    integer :: status 

    file_line_count = 0
    open(unit=30, status="old", file=filename)
    do
       read(30,'(a1)',iostat=status) dummy
       if ( status /= 0 ) exit
       file_line_count = file_line_count + 1 
    end do
    close(30)

  end function file_line_count

!
!     ******************************************************************
!
     SUBROUTINE POLFIT(X,Y,SIGMAY,NPTS,NTERMS,MODE,A,CHISQR)
!
!  Bevingtons polynomial least-sqaures fitting routine
!
!
!     DECLARE EVERYTHING TO BE REAL*8
!
      IMPLICIT none
      real(double) :: chisqr, xi, yi, weight, xterm, yterm, delta, free
      INTEGER :: NMAX, NTERMS, N, J, I, NPTS, MODE, K, L
      real(double) :: chisq
!
!     SET UP ARRAYS
!
      real(double) ::  X(:), Y(:), SIGMAY(:), A(:)
      real(double) ::  ARRAY(100,100), SUMX(10000), SUMY(10000)
!
!     ACCUMULATE WEIGHTED SUMS
!
      NMAX  =  2*NTERMS - 1
      DO N  =  1, NMAX
        SUMX(N)  =  0.0D0
      ENDDO
      DO J  =  1, NTERMS
        SUMY(J)  =  0.0D0
      ENDDO
      CHISQ  =  0.0D0
      DO I  =  1, NPTS
        XI  =  X(I)
        YI  =  Y(I)
        IF (MODE.LT.0) THEN
          IF (YI.LT.0) THEN
            WEIGHT  =  -1.0D0/YI
          ELSEIF (YI.EQ.0) THEN
            GOTO 50
          ELSE
            WEIGHT  =  1.0D0/YI
          ENDIF
        ELSEIF (MODE.EQ.0) THEN
          GOTO 50
        ELSE
          WEIGHT  =  1.0D0/SIGMAY(I)**2
        ENDIF
        GOTO 100
   50   CONTINUE
        WEIGHT  =  1.0D0
  100   CONTINUE
        XTERM  =  WEIGHT
        DO N  =  1, NMAX
          SUMX(N)  =  SUMX(N) + XTERM
          XTERM  =  XTERM*XI
        ENDDO
        YTERM  =  WEIGHT*YI
        DO N  =  1, NTERMS
          SUMY(N)  =  SUMY(N) + YTERM
          YTERM  =  YTERM*XI
        ENDDO
        CHISQ  =  CHISQ + WEIGHT*YI**2
      ENDDO
!
!     CONSTRUCT MATRICES AND CALCULATE COEFFICIENTS
!
      DO J  =  1, NTERMS
        DO K  =  1, NTERMS
          N  =  J + K - 1
          ARRAY(J,K)  =  SUMX(N)
        ENDDO
      ENDDO
      DELTA  =  DETERM(ARRAY,NTERMS)
      IF (DELTA.NE.0.d0) THEN
        DO L  =  1, NTERMS
          DO J  =  1, NTERMS
            DO K  =  1, NTERMS
              N  =  J + K - 1
              ARRAY(J,K)  =  SUMX(N)
            ENDDO
            ARRAY(J,L)  =  SUMY(J)
          ENDDO
          A(L)  =  DETERM(ARRAY,NTERMS)/DELTA
        ENDDO
!
!     CALCULATE CHI SQUARE
!
        DO J  =  1, NTERMS
          CHISQ  =  CHISQ - 2.0D0*A(J)*SUMY(J)
          DO K  =  1, NTERMS
            N  =  J + K - 1
            CHISQ  =  CHISQ + A(J)*A(K)*SUMX(N)
          ENDDO
        ENDDO
        FREE  =  dble(NPTS-NTERMS)
        CHISQR  =  CHISQ/FREE
      ELSE
        CHISQR  =  0.0D0
        DO J  =  1, NTERMS
          A(J)  =  0.0D0
        ENDDO
      ENDIF
    END SUBROUTINE POLFIT

!     ******************************************************************
!     ****  FUNCTION DETERM  *******************************************
!     ******************************************************************
!
! Bevington again
!
      real(double) FUNCTION DETERM(ARRAY,NORDER)
!
!     DECLARE EVERYTHING TO BE REAL*8
!
      IMPLICIT none
      real(double) :: ARRAY(100,100), save
      INTEGER :: I, J, K, K1, NORDER
      DETERM  =  1.0D0

      DO K  =  1, NORDER
!
!     INTERCHANGE COLUMNS IF DIAGONAL ELEMENT IS ZERO
!
        IF (ARRAY(K,K).NE.0) GOTO 100
        DO J  =  K, NORDER
          IF (ARRAY(K,J).NE.0) GOTO 50
        ENDDO
        DETERM  =  0.0D0
        GOTO 200
   50   CONTINUE
        DO I  =  K, NORDER
          SAVE  =  ARRAY(I,J)
          ARRAY(I,J)  =  ARRAY(I,K)
          ARRAY(I,K)  =  SAVE
        ENDDO
        DETERM  =  -DETERM
!
!     SUBTRACT ROW K FROM THE LOWER ROWS TO GET A DIAGONAL MATRIX
!
  100   CONTINUE
        DETERM  =  DETERM*ARRAY(K,K)
        IF (K.LT.NORDER) THEN
          K1  =  K + 1
          DO I  =  K1, NORDER
            DO J  =  K1, NORDER
              ARRAY(I,J)  =  ARRAY(I,J) - &
                            ARRAY(I,K)*ARRAY(K,J)/ARRAY(K,K)
            ENDDO

          ENDDO
        ENDIF
      ENDDO

  200 CONTINUE
    END FUNCTION DETERM

    real(double) function calcPoly(aPoly, nPoly, x) result (y)
      integer :: npoly
      real(double) :: aPoly(:), x
      integer :: i

      y = apoly(nPoly)
      do i = nPoly-1,1,-1
         y = y * x + aPoly(i)
      enddo


    end function calcPoly
      

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
  subroutine luSlv(a, b)
    implicit none
    
    real(double), intent(inout) :: a(:,:)
    real(double), intent(inout) :: b(:)    

    call lured(a)
    
    call reslv(a,b)

  end subroutine luSlv
  
  subroutine lured(a)
    implicit none
    
    real(double), intent(inout)  :: a(:,:)
    
    ! local variables
    integer          :: i, j, k, n                    ! counters
    
    real(double) :: factor                     ! general calculation factor

    n = size(a,1)

    if (n == 1) return
    
    do i = 1, n-1
       do k = i+1, n
          factor = a(k,i)/a(i,i)
          do j = i+1, n
             a(k, j) = a(k, j) - a(i, j) * factor
          end do
       end do
    end do
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
    

  !
  ! Compares the Lorentz profile and Voigt Profile
  !
  subroutine test_profiles()
    implicit none
    
    integer, parameter :: nbin = 100
    real(double) :: Lorentz(nbin), lambda(nbin), freq(nbin), Voigt_prof(nbin)
    integer :: i, j  
    integer :: nsample = 10000000
    real(double) :: nu0, nu_rand, lam0, lam_min, lam_max, Gamma
    real(double) :: lam_rand
    real(double) :: c=2.99792458e10  ! cm/s
    real(double) :: tmp, dlam, a, doppler, dnu, dv, T, N_HI, Ne
    
    lam_min = 6550.0d-8  ! cm
    lam0    = 6562.8d-8  ! cm
    lam_max = 6580.0d-8  ! cm
    
    nu0 = c/lam0
    
    ! sets up wavelength array
    tmp = (lam_max-lam_min)/dble(nbin-1)
    do i = 1, nbin
       lambda(i) = lam_min + tmp*dble(i-1)
       freq(i) = c/lambda(i)
    end do
    
     
    
    ! take random samples
!    Gamma = 3.48d11
    T = 6000.0; N_HI=1.e-19; Ne=1.e-9
    Gamma = bigGamma(N_HI, T, Ne, nu0)
    !Gamma = 0.0d0

    
    dlam = lambda(2) - lambda(1)
    
    Lorentz(1:nbin) = 0.0d0
    do i = 1, nsample
       
       nu_rand = random_Lorentzian_frequency(nu0, Gamma)
       
       lam_rand = c/nu_rand + dlam/2.0
       
       call locate(lambda, nbin, lam_rand, j)
       
       Lorentz(j) = Lorentz(j) + 1.0d0
       
    end do
    
    
    ! No computes the Voigt profile
    doppler = nu0/cSpeed * sqrt(2.*kErg*T/mHydrogen)
    a = Gamma/4.0d0/3.141593d0/doppler
    do i = 1, nbin
       dnu = nu0 - freq(i)
       dv = dnu/doppler
       voigt_prof(i) = voigtn(a, dv)
    end do
    
  
    ! Renomarlize the profiles so that the max is 1
    tmp = MAXVAL(Lorentz(1:nbin)) 
    if (tmp/=0) Lorentz(:) = Lorentz(:)/tmp

    tmp = MAXVAL(voigt_prof(1:nbin)) 
    if (tmp/=0) voigt_prof(:) = voigt_prof(:)/tmp
  
    
    ! writing the results in a file
    
    open(unit=66, file="lorentz_voigt.dat", status="replace")
    do i = 1, nbin
       write(66,*) lambda(i)*1.0d8, lorentz(i), voigt_prof(i)
       !             [a]              
    end do
    close(66)
    
  end subroutine test_profiles

  real(double) function median(y, n)
    real(double) :: y(:)
    integer :: n
    real(double), allocatable :: t(:)

    allocate(t(1:n))

    t(1:n) = y(1:n)
    
    call sort(n, t)
    if (mod(n, 2) == 0) then ! even case
       median = (t(n/2)+t(n/2+1))/2.
    else
       median = t(n/2+1)
    endif
    deallocate(t)
  end function median


      subroutine ludcmp(a,n,np,indx,d)
        integer :: nMax,n,np,indx(:)
        real(double) :: vtiny,d
      parameter (nmax=100)
      real(double) :: a(np,np),vv(nmax)
      integer :: i, j,k,imax
      real(double) :: aamax, sum, dum
      vtiny = tiny(vtiny)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) call writeFatal("ludcmp: singular matrix.")
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.d0)a(j,j)=vtiny
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.d0)a(n,n)=vtiny
      end subroutine

      subroutine lubksb(a,n,np,indx,b)
        integer :: np,n,ii,i,ll,j
        integer :: indx(n)
       real(double) ::  a(np,np),b(n)
       real(double) :: sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b(j)
13        continue
        endif
        b(i)=sum/a(i,i)
14    continue
     end subroutine

      subroutine svdcmp(a,m,n,mp,np,w,v)
        integer :: nMax, m, n, mp, np
      parameter (nmax=100)
      real(double) ::  a(mp,np),w(np),v(np,np),rv1(nmax)
      real(double) :: g, c, f, h, s, scale, anorm, x, y, z
      integer :: i, j, k, l, its, nm
      g=0.d0
      scale=0.d0
      anorm=0.d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.d0
        s=0.d0
        scale=0.d0
        if (i.le.m) then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if (scale.ne.0.d0) then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            if (i.ne.n) then
              do 15 j=l,n
                s=0.d0
                do 13 k=i,m
                  s=s+a(k,i)*a(k,j)
13              continue
                f=s/h
                do 14 k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
14              continue
15            continue
            endif
            do 16 k= i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.d0
        s=0.d0
        scale=0.d0
        if ((i.le.m).and.(i.ne.n)) then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if (scale.ne.0.d0) then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            if (i.ne.m) then
              do 23 j=l,m
                s=0.d0
                do 21 k=l,n
                  s=s+a(j,k)*a(i,k)
21              continue
                do 22 k=l,n
                  a(j,k)=a(j,k)+s*rv1(k)
22              continue
23            continue
            endif
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if (i.lt.n) then
          if (g.ne.0.d0) then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.d0
            v(j,i)=0.d0
31        continue
        endif
        v(i,i)=1.d0
        g=rv1(i)
        l=i
32    continue
      do 39 i=n,1,-1
        l=i+1
        g=w(i)
        if (i.lt.n) then
          do 33 j=l,n
            a(i,j)=0.d0
33        continue
        endif
        if (g.ne.0.d0) then
          g=1.d0/g
          if (i.ne.n) then
            do 36 j=l,n
              s=0.d0
              do 34 k=l,m
                s=s+a(k,i)*a(k,j)
34            continue
              f=(s/a(i,i))*g
              do 35 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
35            continue
36          continue
          endif
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.d0
38        continue
        endif
        a(i,i)=a(i,i)+1.d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if ((abs(rv1(l))+anorm).eq.anorm)  go to 2
            if ((abs(w(nm))+anorm).eq.anorm)  go to 1
41        continue
1         c=0.d0
          s=1.d0
          do 43 i=l,k
            f=s*rv1(i)
            if ((abs(f)+anorm).ne.anorm) then
              g=w(i)
              h=sqrt(f*f+g*g)
              w(i)=h
              h=1.d0/h
              c= (g*h)
              s=-(f*h)
              do 42 j=1,m
                y=a(j,nm)
                z=a(j,i)
                a(j,nm)=(y*c)+(z*s)
                a(j,i)=-(y*s)+(z*c)
42            continue
            endif
43        continue
2         z=w(k)
          if (l.eq.k) then
            if (z.lt.0.d0) then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            go to 3
          endif
          if (its.eq.30) call writeFatal("no convergence in 30 iterations")
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=sqrt(f*f+1.d0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.d0
          s=1.d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=sqrt(f*f+h*h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 nm=1,n
              x=v(nm,j)
              z=v(nm,i)
              v(nm,j)= (x*c)+(z*s)
              v(nm,i)=-(x*s)+(z*c)
45          continue
            z=sqrt(f*f+h*h)
            w(j)=z
            if (z.ne.0.d0) then
              z=1.d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 nm=1,m
              y=a(nm,j)
              z=a(nm,i)
              a(nm,j)= (y*c)+(z*s)
              a(nm,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.d0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
   end subroutine svdcmp

      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      INTEGER m,mp,n,np,NMAX
      REAL(double) :: b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=500)
      INTEGER i,j,jj
      REAL(double) s,tmp(NMAX)
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      end subroutine svbksb

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
    real :: nuArray(:), fArray(:), nu
    real, allocatable :: nuTemp(:), fTemp(:)
    integer :: i, n


    allocate(nuTemp(1:n), fTemp(1:n))
    if (nuArray(1) < 1.d5) then
       call writeWarning("Surface flux file appears to be in wavelength space: converting")
       do i = 1, n
          nu = cSpeed / (1.d-8 * nuArray(i)) ! hz
          nuTemp(i) = nu
          fTemp(i) = (fArray(i)*1.d8) ! erg/s/cm^2/angs -> erg/s/cm^2/cm
          fTemp(i) = fTemp(i) * cspeed / nu**2 ! erg/s/cm^2/Hz
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

  subroutine bonnorEbertRun(t, mu, rho0,  nr, r, rho)
    use constants_mod
    implicit none
    real(double) :: t, rho0
    integer :: nr
    real(double) :: r(:), rho(:)
    real(double), allocatable :: zeta(:), phi(:)
    real(double) ::  r0
    real(double) :: mu, soundSpeed, mass, zinner, eThermal, eGrav, dv
    integer :: i
    real(double), parameter :: zetaCrit = 5.d0
    real(double) :: dr, drhodr, d2rhodr2
    character(len=80) :: message

    zinner = 0.001d0
    allocate(zeta(1:nr), phi(1:nr))
    zeta = 0.d0
    phi = 0.d0



    do i = 1, nr
      zeta(i) = log10(zinner) + (log10(zetaCrit)-log10(zinner))*dble(i-1)/dble(nr-1)
   enddo
   zeta(1:nr) = 10.d0**zeta(1:nr)
   zeta(1) = 0.d0
   soundSpeed = sqrt((kErg*t) /(mu*mHydrogen))

   r0 = 0.89d0 * soundSpeed / sqrt(bigG * rho0)
!   r0 = 1.6d0*pctocm
   do i = 1, nr
      r(i) = zeta(i) * r0
   enddo
   
   drhodr  = 0.d0
   rho(1) = rho0
   d2rhodr2 = (-fourPi*bigG*rho(1)**2 * (mu*mHydrogen) / (kerg * t)) 
 
   do i = 2, nr
!   do i = 1, nr
      dr = r(i) - r(i-1)
      rho(i) = rho(i-1) + drhodr * dr
      d2rhodr2 = (-fourPi*bigG*rho(i)**2 * (mu*mHydrogen) / (kerg * t)) - (2.d0/r(i))*drhodr + (1.d0/rho(i))*drhodr**2
      drhodr = drhodr + d2rhodr2 * dr
   enddo
      

   mass = 0.d0
   eThermal = 0.d0
   eGrav = 0.d0
   do i = 2, nr
!   do i = 1, nr
      dv = fourPi*r(i)**2*(r(i)-r(i-1))
      mass = mass + rho(i)*dv
      eGrav = eGrav + bigG*dv*rho(i)*mass/r(i)
      eThermal = eThermal + (dv*rho(i)/(mu*mHydrogen))*kerg*t
   enddo
   
   write(message,'(a,f5.2)') "Outer radius of Bonnor-Ebert sphere is (in pc): ",r0/pctocm
    call writeInfo(message, TRIVIAL)
    !print *, "mass", mass
    !print *, "msol", msol
    !write(message,'(a,f5.2)') "Mass contained in Bonnor-Ebert sphere is: ",mass/msol
   call writeInfo(message, TRIVIAL)
   write(message,'(a,1pe12.3)') "Gravitational p.e.  contained in Bonnor-Ebert sphere is: ",eGrav
   call writeInfo(message, TRIVIAL)
   write(message,'(a,1pe12.3)') "Thermal energy contained in Bonnor-Ebert sphere is: ",eThermal
   call writeInfo(message, TRIVIAL)
   write(message,'(a,f6.3)') "Ratio of thermal enery/grav energy: ",eThermal/eGrav
   call writeInfo(message, TRIVIAL)

 end subroutine bonnorEbertRun
      
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
        C = NPTS - 2
        VARNCE = (SUMY2+A*A*SUM+B*B*SUMX2 &
                -2*(A*SUMY+B*SUMXY-A*B*SUMX))/C
      ENDIF
      SIGMAA = SQRT(VARNCE*SUMX2/DELTA)
      SIGMAB = SQRT(VARNCE*SUM/DELTA)
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

end module utils_mod

