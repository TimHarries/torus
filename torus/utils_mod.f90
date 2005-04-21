!
! written by tjh


module utils_mod

  use kind_mod
  use vector_mod          ! vector maths
  use constants_mod       ! physical constants
  use unix_mod
  
  implicit none

  public

  interface locate
     module procedure locate_single
     module procedure locate_octal
  end interface

  interface hunt
     module procedure hunt_single
     module procedure hunt_octal
  end interface

  interface sort
     module procedure sortsingle
     module procedure sortdouble
  end interface

  interface logint
     module procedure logint_single
     module procedure logint_dble
  end interface


contains

  ! solve a quadratic equation


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

    y = b*b-4.d0*a*c

    if (y >= 0.) then
       x1 = (-b + sqrt(y))/(2.d0*a)
       x2 = (-b - sqrt(y))/(2.d0*a)
    else
       write(*,*) "! quad solver failed",y,a,b,c
       ok = .false.
       x1 = -b/(2.d0*a)
       x2 = -b/(2.d0*a)
    endif

  end subroutine solveQuadDble

  ! return a blackbody function

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
    IF (XM.LT.12.)THEN
       IF (XM.NE.OLDM) THEN
          OLDM=XM
          G=EXP(-XM)
       ENDIF
       EM=-1
       T=1.
2      EM=EM+1.

       call random_number(r1)
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
       call random_number(r1)
       Y=TAN(PI*r1)
       EM=SQ*Y+XM
       IF (EM.LT.0.) GO TO 1
       EM=INT(EM)
       T=0.9*(1.+Y**2)*EXP(EM*ALXM-GAMMLN(EM+1.)-G)
       call random_number(r1)
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
       stop
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
       stop
    endif
    LX=LOG(X)
    LX1=LOG(X1)
    LX2=LOG(X2)
    LY1=LOG(MAX(Y1,1.d-90))
    LY2=LOG(MAX(Y2,1.d-90))
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


  ! sort an array

  SUBROUTINE SORTsingle(N,RA)
    INTEGER N, L, IR, I, J
    REAL RRA
    REAL RA(N)
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

  SUBROUTINE SORTdouble(N,RA)
    INTEGER N, L, IR, I, J
    real(double) RRA
    real(double) RA(N)
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

  ! sort an array by indexing

  PURE SUBROUTINE INDEXX(N,ARRIN,INDX)
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
  END subroutine indexx


  !
  !
  ! Taken from Numerical recipes
  ! See the comments in the book.
  SUBROUTINE polint(xa,ya,n,x,y,dy)
    INTEGER n,NMAX
    REAL dy,x,y,xa(n),ya(n)
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
    REAL :: x,xa(n)
    real(double) :: y, ya(n), dy
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

  SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
    INTEGER NMAX,N, I, K
    REAL SIG, P, QN, UN
    PARAMETER (NMAX=1000)
    REAL X(N),Y(N),Y2(N),U(NMAX),YP1,YPN
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
  END SUBROUTINE SPLINE


  SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
    INTEGER N, KLO, KHI, K
    REAL H, A, B
    REAL XA(N),YA(N),Y2A(N),X,Y
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
       write(*,*) 'Bad XA input.'
       stop
    endif
    A=(XA(KHI)-X)/H
    B=(X-XA(KLO))/H
    Y=A*YA(KLO)+B*YA(KHI)+ &
         ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
  END SUBROUTINE SPLINT


  ! locate in a grid via bisection but starting at jlo

  pure SUBROUTINE HUNT_single(XX,N,X,JLO)
    REAL, INTENT(IN) :: XX(*)
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
    real(oct), INTENT(IN) :: XX(*)
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

    call locate(x, ny, xi, i)
    
    t = (xi - x(i))/(x(i+1)-x(i))

    logInterp = 10.e0**(log10(y(i)) + t * (log10(y(i+1))-log10(y(i))))

  end function logInterp

  real(double) function logInterp_dble(y, n, x, xi)
    integer, intent(in) :: n
    real(double), intent(in)    :: y(n), x(n), xi
    integer, save       :: i
    real(double)                :: t

    call locate(x, n, xi, i)
    
    t = (xi - x(i))/(x(i+1)-x(i))

    logInterp_dble = 10.d0**(log10(y(i)) + t * (log10(y(i+1))-log10(y(i))))

  end function logInterp_dble

  PURE SUBROUTINE LOCATE_single(XX,N,X,J)
    real, intent(in)    :: XX(*)
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
  

    PURE SUBROUTINE LOCATE_octal(XX,N,X,J)
    real(oct), intent(in) :: XX(*)
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

    END SUBROUTINE LOCATE_octal

    real function gasdev()
      implicit none
      integer :: iset
      real :: r1, r2, gset
      real :: v1, v2, r, fac
      data iset/0/
      save iset,gset
      if (iset.eq.0) then
1        continue
         call random_number(r1)
         call random_number(r2)
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

    real function gauss(sigma, dx)
      real :: sigma, dx, fac
      fac = 1./(sigma*sqrt(twoPi))
      gauss = fac * exp(-(dx**2)/(2.*sigma**2))
    end function gauss


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
      

    subroutine systemInfo(startTime,nPhotons)
      integer :: nPhotons
      character(len=80) :: time, login, hostname
      integer :: cpuTime, endTime, startTime, i

      call unixTimes(cpuTime, endTime)
      call unixGethostname(hostname,i)
      call unixGetlogin(login,i)

      write(*,'(a)') " "
      write(*,'(a)') "System information"
      write(*,'(a)') "------------------"
      write(*,'(a)') " "
      write(*,'(a,a)') "Job run by: ",trim(login)
      write(*,'(a,a)') "Job run on: ",trim(hostname)
      write(*,'(a)') " "
      call timeString(endTime - startTime, time)
      write(*,'(a,a)') "Time elapsed: ",time
      call timeString(cpuTime, time)
      write(*,'(a,a)') "CPU time used: ",time
      write(*,'(a)') " "
      write(*,'(a,f7.1)') "Photons/second: ",real(nPhotons)/real(endTime-startTime)

    end subroutine systemInfo

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

      real, intent(in) :: mass, temperature
      real :: x, y, t, u, vel
      logical :: ok

      ok = .false.
      do while(.not.ok)

         call random_number(x)
         x = 3. * x
         call random_number(y)
         
         t = 4./sqrt(pi) * x**2 * exp(-x**2)
         
         if (y < t) then
            ok = .true.
            u = x
         endif

      end do

      vel = u * sqrt((2*kErg*temperature)/mass)

      maxwellianVelocity = vel * randomUnitVector()
    end function maxwellianVelocity



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
      real(double) :: x, y, t, u, vel, a, a2
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
         call random_number(x)
         ! random number between -6a to 6a
         x =6.0d0*a*(-1.0d0 + 2.0d0*x) 
         ! trial value of phi(nu)
         call random_number(y)

         t = a2/(x*x +a2)  ! This should be always less than 1.
         
         if (y < t) then
            ok = .true.
            u = x
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
      real :: xArray(nLambda)
      logical :: ok
      
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
      REAL*8 VOIGTN,AA,VV
      !
      REAL*8 H(25)
      DATA C1,C2/1.128379167095512D0  ,5.64189583547756D-1/
      SAVE C1,C2
      !
      REAL*8 V,V2,V4,V6
      REAL*8 A,A2,A4,A6
      REAL*8 Z,Z2                       
      REAL*8 X,W
      REAL*8 C1,C2
      !
      INTEGER*4 I,J
      !
      !REAL*8 DAWSON
      !EXTERNAL DAWSON
      !
      V=ABS(VV)
      A=AA
      V2=V*V
      A2=A*A
      Z=A2+V2
      IF(A .EQ. 0.0D0  ) THEN
         ! ---- Normal Doppler brodening.         
         IF (V2 < 1.D-10) THEN
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
      REAL*8 DAWSON,XX
      REAL*8 X,U
      REAL*8 UP,DOWN
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

      ! The following are set in input_variables module.
      !  
      ! real, parameter :: C_rad  =  0.   ! [A]  
      ! real, parameter :: C_vdw =   0.   ! [A]
      ! real, parameter :: C_stark = 0.   ! [A]

      ! see Muzerolle et al. 2001 ApJ 550 944

      bigGamma = &
           &    C_rad  &
           &  + C_vdw*(N_HI / 1.e16)*(temperature/5000.)**0.3 &
           &  + C_stark*(Ne/1.e12)**0.6666  ! [Angstrom]

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
      real :: dProjVel(nTau)  ! automatic array
      integer :: nAdd 
      integer :: i, j
!      real, parameter :: dvel  = 10.e5/cSpeed
      real :: dvel, dlam

      ! -- using the values in input_variables module
      dvel = (lamend-lamstart)/lamline/real(nlambda-1)  ! should be in [c]
      dvel = dvel/4.0  ! to be safe  
!      dvel = dvel/10.0  ! to be safe  
!      dProjVel(1)  = 0.
      dProjVel(1:nTau-1) = projVel(2:nTau) - projVel(1:nTau-1)
      dProjVel(nTau)  = 0.0

      newNTau = 0
      do i = 2, nTau
         ! if (i==2) we always add points
         if (abs(dProjVel(i)) > dVel .or. i==2) then
!            nAdd = nint(abs(ProjVel(i))/dVel)
            nAdd = nint(abs(dProjVel(i))/dVel)
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
      integer :: nClose = 10
      integer :: nAdd = 10
      integer :: i, j, n
      real :: dlam
      
      ! loop through the projected velocities to find the resonance zones
      newNTau = 0

      ! always add points between points near the current location of photon
      n = min(nClose, nTau)
      ! check for the resonance zone(s)
      do i = 2, nTau         
         if ( ( i<n ) &
              .or.    &            
            ( (projVel(i-1) <= thisVel) .and. (thisVel < projVel(i)) ) &
              .or. &
              ( (projVel(i-1) >= thisVel) .and. (thisVel > projVel(i)) )   ) then
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


    
    subroutine resampleRay_tau(L, nTau, tau, dtau_max, maxtau, newL, newNTau,&
         inFlow, newInFlow)
      integer, intent(in) :: nTau, maxtau
      real, intent(in) :: L(nTau)
      real, intent(in) :: tau(nTau)
      ! maxmum line optical depth increment allowed
      real, intent(in) :: dtau_max
      integer, intent(inout) :: newNtau
      real, intent(inout)  :: newL(maxtau)
      logical, intent(inout)  :: InFlow(nTau)
      logical, intent(inout)  :: newInFlow(maxtau)
      !
      real, allocatable, save :: dtau(:) ! dtau(maxtau)  ! automatic array
      integer :: nAdd 
      integer :: i, j
!      real, parameter :: dvel  = 10.e5/cSpeed
      real ::  dL
      logical, save  :: first_time = .true.      
      !
      ! allocates memotry for work array for the first time.
      if (first_time) then
         first_time=.false.
         ALLOCATE(dtau(maxtau))
      end if


      ! Initial optical depth in each line sequemnt
      dtau(1:nTau-1) = tau(2:nTau) - tau(1:nTau-1)

      newNTau = 0
      do i = 1, nTau-1
         if (dtau(i) > dtau_max) then
            nAdd = nint(dtau(i)/dtau_max)  ! this should be at least 1
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
      real :: xArray(:), yArray(:)
      integer :: nx, newNx
      integer, intent(in) :: nx_max
      real :: newXarray(:)
      real, save, allocatable :: newYarray(:) ! newYarray(newNx) ! automatic array
      integer :: i, j
      logical, save :: first_time = .true.

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
      real :: xArray(:)
      real(double) :: yArray(:)
      integer :: nx, newNx
      integer, intent(in) :: nx_max
      real :: newXarray(:)
      real, save, allocatable  :: newYarray(:) ! newYarray(newNx)  ! automatic array
      integer :: i, j
      logical, save :: first_time = .true.

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



    subroutine log_linear_resample(xArray, yArray, nX, newXarray, newNx)
      implicit none
      real, intent(in)    :: xArray(:)
      real, intent(inout) :: yArray(:)
      integer, intent(in) :: nx, newNx
      real, intent(in)    :: newXarray(:)
      !
      real :: newYarray(newNx) ! automatic array
      integer :: i, j
      real:: log_y2, log_y1, x2, x1, T1

      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         if (j >= Nx) j = j -1
         x1 = xArray(j); x2 = xArray(j+1)
         if (x1 /= x2) then
            log_y1 = LOG(yArray(j))
            log_y2 = LOG(yArray(j+1))
            T1 = (log_y2-log_y1)/(x2-x1)
            newYarray(i) = log_y1 + T1*(newXarray(i)-x1)
            newYarray(i) = EXP(newYarray(i))
         else
            newYarray(i) = yArray(j)
         endif
      enddo
      yArray(1:newNx) = newYArray(1:newNx)
    end subroutine log_linear_resample



    subroutine log_linear_resample_dble(xArray, yArray, nX, newXarray, newNx)
      implicit none      
      real, intent(in)            :: xArray(:)
      real(double), intent(inout) :: yArray(:)
      integer, intent(in)         :: nx, newNx
      real, intent(in)            :: newXarray(:)
      !
      real(double) :: newYarray(newNx) ! automatic array
      integer      :: i, j
      real(double) :: log_y2, log_y1, x2, x1, T1

      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         if (j >= Nx) j = j -1
         x1 = xArray(j); x2 = xArray(j+1)
         if (x2 /= x1) then
            log_y1 = LOG(yArray(j))
            log_y2 = LOG(yArray(j+1))
            T1 = (log_y2-log_y1)/(x2-x1)
            newYarray(i) = log_y1 + T1*(newXarray(i)-x1)
            newYarray(i) = EXP(newYarray(i))
         else
            newYarray(i) = yArray(j)
         endif
      enddo
      yArray(1:newNx) = newYArray(1:newNx)
    end subroutine log_linear_resample_dble


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



 
  FUNCTION arth(first,increment,n)
    ! Array function returning an arithmetic progression. 
    REAL, INTENT(IN) :: first,increment 
    INTEGER, INTENT(IN) :: n 
    REAL, DIMENSION(n) :: arth 
    INTEGER :: k,k2 
    REAL :: temp 
    INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
    if (n > 0) arth(1)=first 
    if (n <= NPAR_ARTH) then 
      do k=2,n
        arth(k)=arth(k-1)+increment 
      end do
    else
      do k=2,NPAR2_ARTH 
        arth(k)=arth(k-1)+increment 
      end do 
      temp=increment*NPAR2_ARTH 
      k=NPAR2_ARTH
      do
        if (k >= n) exit 
        k2=k+k 
        arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
        temp=temp+temp 
        k=k2
      end do 
    end if
  END FUNCTION arth 
  
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

  subroutine stripSimilarValues(x, nx, xtol)
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
       if (abs((xtemp(newnx)-x(i))/x(i)) > xTol) then
          newnx = newnx  + 1
          xtemp(newnx) = x(i)
       endif
    enddo

    x(1:newnx) = xtemp(1:newnx)
    nx = newnx
  end subroutine stripSimilarValues
       

  subroutine convertByte4(iByte, ival)
    integer(kind=1) :: ibyte(4), iswap
    integer(kind=2) :: ubyte(4)
    integer  :: ival

    ubyte = ibyte
    where (ubyte < 0) ubyte = ibyte + 256


    ival = ubyte(1) + 256*ubyte(2) + 65536 * ubyte(3) + 16777200 * ubyte(4)

  end subroutine convertByte4

  subroutine convertByte2(iByte, ival)
    integer(kind=1) :: ibyte(2), iswap
    integer(kind=2) :: ubyte(2)
    integer(kind=2)  :: ival

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


end module utils_mod




