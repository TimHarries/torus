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
     module procedure locate_double
  end interface

  interface hunt
     module procedure hunt_single
     module procedure hunt_double
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

    if (y >= -1.e-2) then
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
       stop
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

  SUBROUTINE SORT(N,RA)
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
  END subroutine sort

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

  pure SUBROUTINE HUNT_double(XX,N,X,JLO)
    REAL(kind=doubleKind), INTENT(IN) :: XX(*)
    INTEGER, INTENT(IN) :: N
    REAL(kind=doubleKind), INTENT(IN)    :: X
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
  end SUBROUTINE HUNT_double

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
  
    PURE SUBROUTINE LOCATE_double(XX,N,X,J)
    real(kind=doubleKind), intent(in) :: XX(*)
    integer, intent(in)               :: n
    real(kind=doubleKind), intent(in) :: x
    integer,intent(out)               :: j
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
    END SUBROUTINE LOCATE_double

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
      character*(*), intent(in) :: a
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
      character*(*), intent(in) :: strings
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

    

end module utils_mod




