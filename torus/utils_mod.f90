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
    real(kind=doubleKind),intent(in) :: a, b, c
    real(kind=doubleKind),intent(out):: x1, x2
    real(kind=doubleKind)            :: y
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

  real pure function blackBody(temperature, wavelength)

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
    REAL(kind=doubleKind):: COF(6),STP,HALF,ONE,FPF,X,TMP,SER
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

  REAL FUNCTION LOGINT(X,X1,X2,Y1,Y2)
    IMPLICIT NONE
    REAL X,X1,X2,Y1,Y2,ANS
    REAL LX,LX1,LX2,LY1,LY2,GR
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
       write(*,*) 'bad x in logint'
       stop
    endif
    GR=(LY2-LY1)/(LX2-LX1)
    ANS=LY1+GR*(LX-LX1)
    LOGINT=EXP(ANS)
  END FUNCTION LOGINT


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
    REAL(kind=doubleKind) RRA
    REAL(kind=doubleKind) RA(N)
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
    REAL, INTENT(IN)     :: ARRIN(*)
    INTEGER, INTENT(OUT) :: INDX(*)
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
    REAL(kind=octalKind), INTENT(IN) :: XX(*)
    INTEGER, INTENT(IN) :: N
    REAL(kind=octalKind), INTENT(IN)    :: X
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
    real,intent(in)    :: x(*), y(*)           
    real,intent(out)   :: derivs(*)
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
    real, intent(in)    :: y(*), x(*), xi
    integer, intent(in) :: ny
    integer, save       :: i
    real                :: t

    call locate(x, ny, xi, i)
    
    t = (xi - x(i))/(x(i+1)-x(i))

    logInterp = 10.e0**(log10(y(i)) + t * (log10(y(i+1))-log10(y(i))))

  end function logInterp

  PURE SUBROUTINE LOCATE_single(XX,N,X,J)
    real, intent(in)    :: XX(*)
    integer,intent(in)  :: n
    real,intent(in)     :: x
    integer,intent(out) :: j
    integer :: jl, ju,jm
      JL=1
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
      J=JL
    END SUBROUTINE LOCATE_single
  

    PURE SUBROUTINE LOCATE_octal(XX,N,X,J)
    real(kind=octalKind), intent(in) :: XX(*)
    integer, intent(in)              :: n
    real(kind=octalKind), intent(in) :: x
    integer,intent(out)              :: j
    integer :: jl, ju,jm
      JL=1
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
      J=JL
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


    ! a functions to convert characters to integer 
    function char2int(a) RESULT(out)
      implicit none 
      integer :: out
      character(LEN=*), intent(in) :: a
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

    real(kind=doubleKind) function interpLinearDouble(xArray, yArray, n, x)
      real(kind=doubleKind) :: xarray(:), yArray(:)
      real(kind=doubleKind) :: x, t
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

    real(kind=doubleKind) function interpLogLinearDouble(xArray, yArray, n, x)
      real(kind=doubleKind) :: xarray(:), yArray(:)
      real(kind=doubleKind) :: x, t
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
      real(kind=doubleKind) :: flux   ! in erg/s/cm^2/A
      real(kind=doubleKind) :: wavelength ! in A
      real(kind=doubleKind) :: newFlux, nu

      nu  = cspeed / (wavelength * angstromtocm)
      newFlux = flux * (cSpeed * 1.e8)/ nu**2 ! from /A to /Hz

      newFlux = newFlux * 1.e23 ! to janskies
    end function convertToJanskies

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


    subroutine resampleRay(lambda, nTau, projVel, newLambda, newNTau)
      real :: lambda(:)
      real :: projVel(:)
      integer :: nTau, newNtau
      real :: newLambda(:)
      real,allocatable :: dProjVel(:)
      integer :: nAdd = 5
      integer :: i, j
      allocate(dProjVel(1:nTau))
      dProjVel  = 0.
      dProjVel(2:nTau) = projVel(2:nTau) - projVel(1:nTau-1)

      newNTau = 0
      do i = 2, nTau
         if (abs(dProjVel(i)) > 50.) then
            do j = 1, nAdd
               newNtau = newNtau + 1
               newLambda(newNTau) = real(j-1)/real(nAdd) * real(lambda(i)-lambda(i-1)) + lambda(i-1)
            enddo
         else
            newNtau = newNtau + 1
            newLambda(newNTau) = lambda(i-1)
         endif
      enddo
      newNtau = newNtau + 1
      newLambda(newNTau) = lambda(nTau)
      deallocate(dProjVel)
    end subroutine resampleRay
               
    subroutine linearResample(xArray, yArray, nX, newXarray, newNx)
      real :: xArray(:), yArray(:)
      integer :: nx, newNx
      real :: newXarray(:)
      real, allocatable :: newYarray(:)
      integer :: i, j

      allocate(newYarray(1:newNx))
      do i = 1, newNx
         call hunt(xArray, nx, newXarray(i), j)
         newYarray(i) = yArray(j) + yArray(j+1)*(newXarray(i)-xArray(j))/(xArray(j+1)-xArray(j))
      enddo
      yArray(1:newNx) = newYArray(1:newNx)
      deallocate(newYarray)
    end subroutine linearResample



    function spiraldist(x,y,k,phase) result(d)
      integer :: i, nr = 100, imin
      real :: x,y,k,d,theta,ktheta
      real :: r,t1,t2,theta1,theta2,phase
      logical :: converged
      d = 1.e30
      theta1 = -phase
      theta2 = twoPi-phase
      converged = .false.
      do while(.not.converged)
         do i = 1, nr
            theta = theta1+(theta2-theta1)*real(i-1)/real(nr-1)
            ktheta = k * (theta+phase)
            r = (x-ktheta*cos(theta+phase))**2 + (y-ktheta*sin(theta+phase))**2
            if (r < d) then
               d=r
               imin = i
            endif
         enddo
         t1 = theta1+(theta2-theta1)*real(imin-2)/real(nr-1)
         t2 = theta1+(theta2-theta1)*real(imin)/real(nr-1)
         theta1 = max(t1,-phase)
         theta2 = t2
         if (abs(theta1-theta2)/twopi < 1.e-4) converged = .true.
      enddo
      d = sqrt(d)
    end function spiraldist


end module utils_mod




