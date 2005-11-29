module subs
  ! used for plot_discwind

  public:: &
       mdot_on_disc, &
       temperature_disc, &
       speed_of_sound, &
       escape_velocity, &
       Vp, &
       Vphi, &
       discwind_density
  


contains
  
  !
  !
  ! Mass loss rate per unit are on the disc surface [g/s/cm^2]
  function mdot_on_disc(w, Rmin, Rmax, alpha, gamma) RESULT(out)
    implicit none
    real*8 :: out
    !
    real*8, intent(in) :: w      ! [10^10cm]  cylindical diatance from the center
    real*8, intent(in) :: Rmin   ! [10^10 cm]
    real*8, intent(in) :: Rmax   ! [10^10 cm]
    real*8, intent(in) :: alpha  ! power law index
    real*8, intent(in) :: gamma  ! power law index
    !
    real*8 :: Mdot, R_max, R_min, R, delta, fac, p
    real*8, parameter ::  pi = 3.1415927d0
    
    delta =4.0d0*alpha*gamma
    p = delta + 2.0d0
    ! converting everthing in cgs
    R = w*1.d10  !  [cm]
    R_min = Rmin * 1.0d10 ! converting from [10^10cm] to [cm]
    R_max = Rmax * 1.0d10 ! converting from [10^10cm] to [cm]
    Mdot = Mdot*1.989d33/3.1536d07   ! Converting from [Msun/yr] to [g/s]

    fac = p / ( R_max**p - R_min**p )

    ! just checking
    if (fac < 0.0d0) then
       write(*,*) "Error:: fac <0.0 in [discwind_class::mdot_on_disc]."
       stop
    end if

    out = Mdot*(R**delta)*fac/(4.0d0*pi)   ! [g/s/cm^2]

  end function mdot_on_disc


  !
  ! disc temperature (at z=0) power law  [K]
  ! Temperature : T(R)= Tmax*(R/Rmin)^gamma where R is the distance from the center
  ! (private)
  function temperature_disc(R, Rmin, Tmax, gamma) RESULT(out)
    implicit none
    real*8 :: out
    real*8, intent(in) :: R      ! [10^10cm] of distance from z-axis
    real*8, intent(in) :: Rmin   ! [10^10cm] the inner radius of the disc
    real*8, intent(in) :: Tmax   ! [K] temperature of the disc at the inner edge of the disc.
    real*8, intent(in) :: gamma  ! index of the temperature power law.
    
    out = Tmax * (r/Rmin)**gamma  ! [K]

    if (out < 3.0d0) out = 3.0d0  ! setting the minimum temperature

  end function temperature_disc



  !
  ! Approximated speed of sound in [cm/s]
  ! (private)
  real*8 function speed_of_sound(r, Rmin, Tmax, gamma)  
    implicit none
    real*8, intent(in) :: R ! [10^10cm] of distance from z-axis
    real*8, intent(in) :: Rmin   ! [10^10cm] the inner radius of the disc
    real*8, intent(in) :: Tmax   ! [K] temperature of the disc at the inner edge of the disc.
    real*8, intent(in) :: gamma  ! index of the temperature power law.
        
    ! using the function in this module
    speed_of_sound = 10.0d0*SQRT(temperature_disc(R, Rmin, Tmax, gamma)/1.0d4) ! [km/s]
    speed_of_sound = speed_of_sound*1.d5  ! [cm/s]

  end function speed_of_sound



  !
  ! Vesc(R) = SQRT(2GM/R)
  ! (private) 
  real*8 function escape_velocity(r, Mstar) 
    implicit none
    real*8, intent(in)  ::  r ! [10^10cm] of distance from the center of the star
    real*8, intent(in)  :: Mstar ! mass of the star in solar mass

    real*8 :: Mass, Radius
    
    Radius = r*1.0d8  ! [m]
    Mass = Mstar*1.98892d30   ! [kg]

    escape_velocity = SQRT(2.0d0*6.673d-11*Mass/Radius) * 1.d2 ! [cm/s]
    

  end function escape_velocity



  real*8 function Vp(l, wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta, Mstar)
    implicit none
    !
    real*8, intent(in) :: l      ! [10^10cm]  distance from the disc along a stream line
    real*8, intent(in) :: wi     ! [10^10cm]  cylindical diatance from the center to the foot point
    real*8, intent(in) :: f      ! [-]  
    real*8, intent(in) :: Rmin   ! [10^10 cm]
    real*8, intent(in) :: Rmax   ! [10^10 cm]
    real*8, intent(in) :: Tmax   ! [K] temperature of the disc at the inner edge of the disc.
    real*8, intent(in) :: alpha  ! power law index
    real*8, intent(in) :: gamma  ! power law index
    real*8, intent(in) :: d      ! [10^10cm]  teh distance from the origin to a source.
    real*8, intent(in) :: beta   ! [-] wind acceralation parameter
    real*8, intent(in)  :: Mstar ! mass of the star in solar mass
    !
    real*8 :: r  ! [10^10cm]  distance from the origin.
    real*8 :: s  ! [10^10cm]  distance from the disc
    real*8 :: wi_sq   ! [10^20 cm]
    real*8 :: Cs     ! speed of sound in [cm/s]
    real*8 :: fac, tmp, Rs
    wi_sq = wi*wi
    ! 
    Cs = speed_of_sound(wi, Rmin, Tmax, gamma)  ! [cm/s]    
    Rs = 10.0d0*Rmin
    fac = 1.0d0 -  Rs/(l+Rs)  ! [-]
    fac = fac**beta


    tmp = Cs + (f*escape_velocity(wi, Mstar) - Cs)*fac  ! [cm/s]
    Vp = tmp
 
  end function Vp




  !
  ! Assuming that angular momentum around z-axis conserves...
  ! The output is in cm/s.
  real*8 function Vphi(l, wi, d, Mstar)
    implicit none
    !
    real*8, intent(in) :: l      ! [10^10cm]  distance from the disc along a stream line
    real*8, intent(in) :: wi     ! [10^10cm]  cylindical diatance from the center to the foot point
    real*8, intent(in) :: d      ! [10^10cm]  teh distance from the origin to a source.
    real*8, intent(in)  :: Mstar ! mass of the star in solar mass
    !
    real*8 :: Vk, w, lp, omega
    real*8, parameter :: G = 6.67259d-8 ! in cgs
    real*8 :: M_star


    ! Now compute the Keplerian orbital speed at this radius


    M_star = Mstar * 1.9891d33  ! converting from [Msun] to [g]

    if (wi/=0) then
       Vk = SQRT( G*M_star / (wi*1.0d10) )  ! [cm/s]
    else
       Vk = 0.0d0
    end if
    
    lp = SQRT(wi*wi + d*d)
    omega = (l/lp)*wi    
    w = wi + omega

    ! Now using the conservation of angular momentum around z...   
    if (w /= 0.0d0) then
       Vphi = Vk * (wi/w)   ! [cm/s]
    else
       Vphi = 0.0d0  ! [cm/s]
    end if
 
  end function Vphi



  ! Given a position (x, y and z) in 10^10 cm, this routine
  ! will return the density [g/cm^3] at the point.
  !
  function discwind_density(l, wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta, Mstar) RESULT(out)
    implicit none
    !
    real*8 :: out  ! in [g/cm^3]
    !
    real*8, intent(in) :: l      ! [10^10cm]  distance from the disc along a stream line
    real*8, intent(in) :: wi      ! [10^10cm]  cylindical diatance from the center to the foot point
    real*8, intent(in) :: f      ! [-]  
    real*8, intent(in) :: Rmin   ! [10^10 cm]
    real*8, intent(in) :: Rmax   ! [10^10 cm]
    real*8, intent(in) :: Tmax   ! [K] temperature of the disc at the inner edge of the disc.
    real*8, intent(in) :: alpha  ! power law index
    real*8, intent(in) :: gamma  ! power law index
    real*8, intent(in) :: d      ! [10^10cm]  teh distance from the origin to a source.
    real*8, intent(in) :: beta   ! [-] wind acceralation parameter
    real*8, intent(in)  :: Mstar ! mass of the star in solar mass
    !
    real*8, parameter :: pi = 3.14159265d0
    !
    real*8, parameter :: density_min = tiny(density_min) ! [g/cm^3]
    real*8, parameter :: density_max = 1.0d10            ! [g/cm^3]

    real*8 :: lp, omega, h, cos_delta, fac, Q

    lp = SQRT(wi*wi + d*d)
    h = (l/lp)*d
    omega = (l/lp)*wi
    cos_delta = d/lp
    Q = l + lp

    fac = d/(Q*cos_delta)

    out = mdot_on_disc(wi, Rmin, Rmax, alpha, gamma)*fac*fac &
         /Vp(l, wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta, Mstar)  &
         /cos_delta

    if (out > density_max) out=density_max

    
  end function discwind_density



end module subs
  
