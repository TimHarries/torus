module atom_mod
  use constants_mod, only: cSpeed, hConst, mElectron, eCharge, pi, kConst, hCgs, kErg
  use kind_mod, only: double
  implicit none
  public

contains

  subroutine getLTEopacities(rne, rtemp, retal, rchil, reta, rchi, resec, m, n)
    implicit none
    
    integer :: k
    integer :: i
    real, intent(in) :: rne, rtemp
    real, intent(inout) :: retal, rchil, reta, rchi, resec
    integer, intent(in), optional :: m,n
    
    real(double) ::  thresh
    real(double) ::  transe
    real(double) ::  nsaha(15),freq
    real(double) ::  temp
    real(double) ::  ne
    real(double) ::  chi
    real(double) ::  eta
    real(double) ::  esec
    real(double) ::  etal
    real(double) ::  chil
    

    
    real(double) :: kappa
    real(double) ::  ah(6,6)
    real(double) ::  bh(6,6)
    real(double) ::  gh(15)
    real(double) ::  eh(10)
    real(double) ::  fh(6,6)
    real(double) :: lambdah(6,6)

    DATA EH &
    /  0.000D0, 10.199D0, 12.088D0, 12.749D0, 13.055D0, 13.221D0, &
     13.321D0, 13.386D0, 13.431D0, 13.463D0 /

      DATA LAMBDAH &
     / 0.00000D-8 , 1215.67D-8 , 1025.72D-8 , 992.537D-8 , 949.743D-8 , 937.803D-8 , &
       1215.67D-8 , 0000000D-8 , 6562.80D-8 , 4861.32D-8 , 4340.36D-8 , 4101.73D-8 , &
       1025.72D-8 , 6562.80D-8 , 0000000D-8 , 18751.0D-8 , 12818.1D-8 , 10938.1D-8 , &
       992.537D-8 , 4861.32D-8 , 18751.0D-8 , 0.00000D-8 , 40512.0D-8 , 26252.0D-8 , &
       949.743D-8 , 4340.46D-8 , 12818.1D-8 , 40512.0D-8 , 0000000D-8 , 74578.0D-8 , &
       937.803D-8 , 4101.73D-8 , 10938.1D-8 , 26252.0D-8 , 74578.0D-8 , 0000000D-8 / 

      DATA AH &
     / 0.000D0 , 4.699D8 , 5.575D7 , 1.278D7 , 4.125D6 , 1.644D6 , &
       0.000D0 , 0.000D0 , 4.410D7 , 8.419D6 , 2.530D6 , 9.732D5 ,&
       0.000D0 , 0.000D0 , 0.000D0 , 8.986D6 , 2.201D6 , 7.783D5 ,&
       0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 2.699D6 , 7.711D5 ,&
       0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 1.025D6 ,&
       0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 /
      DATA BH &
     /0.000D+00,0.212D+10,0.151D+09,0.315D+08,0.889D+07,0.341D+07,&
      0.850D+10,0.000D+00,0.314D+11,0.243D+10,0.521D+09,0.169D+09,&
      0.136D+10,0.706D+11,0.000D+00,0.149D+12,0.117D+11,0.256D+10,&
      0.503D+09,0.974D+10,0.265D+12,0.000D+00,0.452D+12,0.351D+11,&
      0.222D+09,0.325D+10,0.324D+11,0.706D+12,0.000D+00,0.107D+13,&
      0.123D+09,0.152D+10,0.103D+11,0.790D+11,0.154D+13,0.000D+00 /

      DATA FH &
     / 0.000D00,-1.041D00,-0.879D-2,-0.181D-2,-0.558D-2,-0.217D-3, &
       4.162D-1, 0.000D00,-0.285D00,-0.298D-1,-0.715D-2,-0.245D-2, &
       7.910D-2, 6.407D-1, 0.000D00,-0.474D00,-0.542D-1,-0.140D-1, &
       2.899D-2, 1.193D-1, 8.421D-1, 0.000D00,-0.664D00,-0.797D-1, &
       1.394D-2, 4.467D-2, 1.506D-1, 1.038D00, 0.000D00,-0.855D00, &
       7.799D-3, 2.209D-2, 5.584D-2, 1.793D-1, 1.231D00, 0.000D00/

!    m = 2 
!    n = 3

    ne = dble(rne)
    temp = dble(rtemp)

   do k=1,15
       gh(k)=2.d0*dble(k*k)
    enddo

    do k=2,6
       do i=1,k-1
          lambdah(i,k)=lambdah(k,i)
          ah(k,i)=ah(k,i)
          freq=cSpeed/lambdah(i,k)
          bh(k,i)=((ah(k,i)*cSpeed*cSpeed)/(2.d0*hConst*freq**3))
          bh(i,k)=bh(k,i)*gh(k)/gh(i)
          fh(k,i)=-gh(i)*fh(i,k)/gh(k)
       enddo
    enddo



    freq=cSpeed/lambdah(m,n)
    transe=abs(eh(n)-eh(m))

    nsaha = 0.
    call boltz_saha(nsaha,ne,temp)

    !
    esec=6.65d-25*ne
    !
    ! calculate the line opacity and emissivity mihalas II eq 10.2
    !
    chil=((pi*eCharge**2)/(mElectron*cSpeed))*fh(m,n)*nsaha(m)
    chil = chil * (1.d0 - exp(-hConst*freq / (kConst * temp)))
    etal=chil*bnu(freq, temp)

    !
    ! continuous opacity.. bound-free and free-free processes (+es)
    !


    ! bound-free

    kappa = 0.d0

    do i = 1, 15
        thresh =( 13.598d0-eh(i))
        if (transe > thresh) then
           kappa = kappa + nsaha(i)*annudouble(i,freq)
        endif
    enddo

    kappa = kappa + ne*ne*alpkka(freq, temp)

    chi = kappa
    
    eta = kappa * bnu(freq, temp)


    retal = real(etal)
    rchil = real(chil)
    reta = real(eta)
    rchi = real(chi)
    resec = real(esec)

!    write(*,*) retal,rchil,reta,rchi,resec
  end subroutine getLTEopacities

  subroutine boltz_saha(nsaha,ne,te)
    !     
    !     this subroutine calculates the level populations from the saha-
    !     boltzmann equation given the temperature and the electron density.
    !     (see mihalas ii equation 5-14)
    !     
    implicit none
    integer maxsaha
    parameter (maxsaha=15)
    real(double) ::  nsaha(maxsaha)
    integer i                        ! loop counter
    real(double) ::  ne
    real(double) ::  ipot            ! hydrogen ionization potential
    real(double) ::  ci              ! cgs constant
    real(double) ::  te              ! wind temperature
    real(double), parameter :: kev =8.6171E-5
    real(double) :: eh(10), gh(10)
      DATA GH &
     / 2.D0,  8.D0,  18.D0,  32.D0,  50.D0, &
      72.D0, 98.D0, 128.D0, 162.D0, 200.D0 /

    DATA EH &
    /  0.000D0, 10.199D0, 12.088D0, 12.749D0, 13.055D0, 13.221D0, &
     13.321D0, 13.386D0, 13.431D0, 13.463D0 /

    !
    ! setup the variables...
    !
    ipot=13.598
    ci=2.07d-16
    !
    ! calculate the level populations (mihalas ii equ 5-14).
    !
    do i=1,15
       nsaha(i)=(ne**2)*gh(i)*ci* &
       (exp( (ipot-eh(i)) / (kev*te) )) /(te**1.5d0)
    enddo
    !
  end subroutine boltz_saha


  real(double)   function annudouble(n,nu)
    !
    ! this subroutine calculates the photoionization x-section
    ! for hydrogen from the n-th level for a given freq photon.
    !
    implicit none
    integer n                   ! the level
    real(double) ::  nu         ! the photon frequency
    real(double) ::  lam        ! the photon wavelength
    lam=cSpeed/nu
    lam=lam*1.d8
    annudouble=1.044d-26*giidble(n,1.d0,lam)*(lam**3)/dble(n)**5
  end function annudouble

  real(double)  function giidble(n,z,wl)
    !
    ! returns bound-free gaunt factors (nicked from idh)
    !
    implicit none
    real(double) ::  z,wl,coeff(6)
    real(double) ::  efree,sum,a,alam
    integer n,i

    data coeff /-0.338276d0, 0.808398d0, -0.59941d0, 0.104292d0, &
                -6.61998d-3, 1.15609d-4 /

    efree=(911.76d0/wl-1.d0/(n*n))/(z*z)
    if (efree.le.(1.d0+2.d0/n)) then
      giidble = giiadble(n,z,wl)
       goto 500
    elseif (efree.gt.10.d0) then
       sum=coeff(1)
       a=1.d0
       efree=log10(efree)
       do i=2,6
          a=a*efree
          sum=sum+coeff(i)*a
       enddo
       giidble=10.d0**sum
       goto 500
    else
       alam=911.76d0/(z*z*(1.d0+2.d0/n)+1.d0/(n*n))
       alam=giiadble(n,z,alam)
       sum=log10(1.d0+2.d0/n)
       efree=log10(efree)
       sum=(efree-sum)/(1.d0-sum)
       giidble=(1.d0-sum)*alam+0.93d0*sum+0.2d0*sum*(1.d0-sum)
       goto 500
    endif
500 continue
  end function giidble

  real(double)  function giiadble(n,z,wl)
    implicit none
    integer n
    real(double) ::  z,wl,a,u,term
    real(double) :: giib
    a=dble(n)
    u=a*a*911.76d0/(wl*z*z)-1.d0
    giib=1.d0+0.1728d0*(u-1.d0)/((n*(u+1.d0))**(2.d0/3.d0))
    term=0.0496d0*(1.d0+u*(u+1.333d0))/(n*(u+1.d0)**(4.d0/3.d0))
    if ((term/giib).le.0.25d0) then
       giiadble=giib-term
       goto 500
    else
       write(*,*) '!screw up in giia'
       giiadble=0.d0
       goto 500
    endif
500 continue
  end function giiadble

  real(double)  function alpkka(freq,t)
    !
    ! this function returns the free-free absorption coefficient for hydrogen
    !
    implicit none
    real(double) ::  freq,t,wav,gauntf
    wav=1.d8*cSpeed/freq
    gauntf=giiia(1.d0,dble(t),dble(wav))
    alpkka=gauntf*3.69d8/((freq**3)*sqrt(t))
  end function alpkka
  !
  ! free-free gaunt factor routine courtesy of idh
  !
  real(double) function giiia (z, t, wl)
    !
    !   ferland's fabulous functional fits
    !
    real(double) :: coeff(28), a(7)
    real(double) :: u, wl, t, ulog, gam2, z, frac,b,c,sum1,sum2,d
    integer i,k,m
    !
    data coeff &
     /1.102D0       ,-0.1085D0     ,0.09775D0     ,-0.01125D0     , &
      1.2D0         ,-0.24016667D0 ,0.07675D0     ,-0.01658333D0  , &
      1.26D0        ,-0.313166667D0,0.15075D0     ,0.00241667D0   , &
      1.29D0        ,-0.4518333D0  ,0.12925D0     ,0.00258333D0   , &
      1.27D0        ,-0.579D0      ,0.092D0       ,-0.003D0       , &
      1.16D0        ,-0.707333D0   ,0.112D0       ,0.0053333333D0 , &
      0.883D0       ,-0.76885D0    ,0.190175D0    ,0.022675D0     / 
    data a &
     /100.D0, 10.D0, 3.D0, 1.D0, 0.3D0, 0.1D0, 0.001D0/ 
    !
    u = 1.44D+8 / (wl*t)
    ulog = log10(u)
    gam2 = 1.58D+5 * z*z/t
    if (gam2.gt.a(7)) go to 10
    i = 7
    k = 7
    frac = 0.5
    go to 60
10  continue
    if (gam2.lt.a(1)) go to 20
    i = 1
    k = 1
    frac = 0.5
    go to 60
20  continue
    do  i = 2, 7
       if (gam2.gt.a(i)) go to 40
    end do
40  continue
    k = i - 1

    b = log10(a(k))
    c = log10(a(i))
    gam2 = log10(gam2)
    frac = abs ((gam2-b) / (b-c))
60  continue
    k = (k-1)*4
    sum1 = coeff(k+1)
    d = 1.0
    do m = 2, 4
       d = d*ulog
       sum1 = sum1 + coeff(k+m)*d
    enddo
    sum1 = sum1 * (1.0 - frac)
    i = (i-1)*4
    sum2 = coeff(i+1)
    d = 1.0
    do  m = 2, 4
       d = d*ulog
       sum2 = sum2 + coeff(i+m)*d
    enddo
    sum2 = sum2 * frac
    giiia = sum1 + sum2
  end function giiia

  real(double) function bNu(nu,T)
    
    real(double) :: fac1, fac2, fac3, nu, T

    fac1 = (2.d0*dble(hCgs)*nu**3)/dble(cSpeed)**2
    fac3 =  (dble(hCgs)*nu)/ (dble(kErg) * T) 
    if (fac3 > 100.d0) then
       fac2 = 0.d0
    else
       fac2 = 1.d0/(exp(fac3) - 1.d0)
    endif
    bNu = fac1 * fac2
  end function bNu

  real(double) function dbNubydT(nu,T)
    
    real(double) :: fac1, fac2, fac3, nu, T

    fac1 = (2.d0*(hcgs*nu**2)**2)/(cSpeed**2 * kErg  * T**2)
    fac3 =  (hCgs*nu)/ (kErg * T) 
    if (fac3 > 100.d0) then
       fac2 = 0.d0
    else
       fac2 = exp(fac3)/(exp(fac3) - 1.d0)**2
    endif
    dbNubydT = fac1 * fac2
!    if (myRankGlobal == 1) write(*,*) nu,T,fac1,fac2,fac3,dbnubydt
  end function dbNubyDt


!!$  real(double) function bLambda(lambda,T)
!!$    
!!$    real(double) :: fac1, fac2, fac3,  T, lambda
!!$
!!$    fac1 = (2.*hCgs*cSpeed**2)/(lambda *1.d-8)**5
!!$    fac3 =  (hCgs * cSpeed)/ (lambda * 1.d-8 * kErg * T) 
!!$    if (fac3 > 100.d0) then
!!$       fac2 = 0.d0
!!$    else
!!$       fac2 = 1.d0/(exp(fac3) - 1.d0)
!!$    endif
!!$    bLambda = fac1 * fac2
!!$  end function bLambda
  
  !
  ! Plancks function B_lambda(T)
  ! in [erg cm^-2 s^-2 cm^-1 sr^-1]
  real(double) function bLambda(lambda,T)
    implicit none
    real(double), intent(in)  :: T       ! temperature in Kelvinslambda
    real(double), intent(in)  :: lambda  ! wavelength in Angstrom
    real(double) :: x, y, lambda_cm

    lambda_cm = lambda *1.d-8  ! wavelength in cm
    x =  (hCgs * cSpeed)/ (lambda_cm * kErg * T) 
    y = (2.0d0*hCgs*cSpeed*cSpeed)/(lambda_cm)**5
    !               ^^^^^^^^^^^^^
    !      It is written this way because of a problem with g95 compiler (RK)
    !

    if (x > 100.d0) then      ! applying Wien's law
       bLambda = y * EXP(-x)
    elseif ( x < 0.01) then   ! applying Raylegh-Jeans law
       bLambda = 2.0d0*(cSpeed/lambda_cm**4) *kErg*T
       ! -- [erg cm^-2 s^-2 cm^-1 sr^-1]
    else
       blambda = y / (exp(x) - 1.d0)
       ! -- [erg cm^-2 s^-2 cm^-1 sr^-1]
    endif
    
  end function bLambda


  real(double) function dbLambdabydT(lambda,T)
    
    real(double) :: fac1, fac2, fac3,  T, lambda

    fac1 = (2.d0*dble(hCgs)**2*dble(cSpeed)**3)/((lambda *1.d-8)**6 * dble(Kerg) * T**2)
    fac3 =  (hCgs * cSpeed)/ (lambda * 1.d-8 * kErg * T) 
    fac2 = exp(fac3)/(exp(fac3) - 1.d0)**2
    dbLambdabydT = fac1 * fac2
  end function dbLambdabyDt

end module atom_mod
