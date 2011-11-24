  !
  ! Compares the Lorentz profile and Voigt Profile
  !
  subroutine test_profiles()
    use atom_mod, only: bigGamma
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


