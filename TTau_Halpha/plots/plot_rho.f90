program plot_rho

  implicit none

  real*8, parameter :: pi =3.141592654d0
  ! Parameters 
  real*8 :: m, theta_wind, Mdot_wind, V0, Vinf, R_mo, beta, r_star

  real*8 :: R_min, R_max
  integer, parameter :: nr = 200
  integer, parameter :: nt = 3 
  real*8 :: r(nr), rmin, rmax, log10_rmin, log10_rmax, log10_dr
  real*8 :: theta(nt)
  real*8 :: density(nr,nt), rho0, r_cm, gamma
  integer :: i, j
!  real*8, parameter :: G = 6.67259d-8 ! in cgs
!  real*8 :: M_star

  ! setup parameters
  R_star = 2.0d0  ! [R_sun]
  m = 2.0d0                           ! [-]
  theta_wind = 80.0d0 * (pi/180.0d0)  ! [rad]
  Mdot_wind  = 1.0d-8 * 1.98892d33 / 31536000.0d0  ! [grams/s]
!  V0 = 20.0       * 1.0d5   ! [cm/s]
!  V0 = 2.0       * 1.0d5   ! [cm/s]
  V0 = 10.0       * 1.0d5   ! [cm/s]
  Vinf = 200.0d0  * 1.0d5   ! [cm/s]
  R_mo = 3.0d0*R_star * 6.960d10  ! [cm]
  beta = 2.0d0
!  gamma = 0.05
!  M_star = 2.0 * 1.9891d33  !  [g]

  ! setting up r and theta
  rmin = R_mo
  rmax = 1000.0d0*R_mo
  log10_rmin = log10(rmin); log10_rmax = log10(rmax)
  log10_dr = (log10_rmax-log10_rmin)/dble(nr-1)
  do i = 1, nr 
     r(i) = log10_rmin + log10_dr*(i-1)
     r(i) = 10.0d0 ** r(i)
  end do

  theta(1) = 0.0d0
  theta(2) = pi/3.0
  theta(3) = theta_wind - theta_wind/100.0d0  

  do j = 1, nt
     do i = 1, nr
        density(i,j) = rho(r(i), theta(j), m, theta_wind, Mdot_wind, &
             V0, Vinf, R_mo, beta) ! [g/cm^2]
     end do
  end do

  ! normalize the density at r=R_mo and theta = 0
!  rho0 = rho(R_mo, 0.0d0, m, theta_wind, Mdot_wind, &
!       V0, Vinf, R_mo, beta)
  rho0 = rho(R_mo, 0.0d0, m, theta_wind, Mdot_wind, &
       V0, Vinf, R_mo, 1.0d0) ! normaizing with beta=1.0 density

  density(:,:) = density(:,:)/rho0  

  ! normalise r value with r_star
  r(:) = r(:) / R_mo
  
  
  ! now write results
  open(unit = 66, file="density_wind_pol.dat", status="replace")
  open(unit = 67, file="density_wind_mid.dat", status="replace")
  open(unit = 68, file="density_wind_vr.dat", status="replace")
  open(unit = 69, file="density_wind_dvr_dr.dat", status="replace")
  open(unit = 70, file="density_wind_pol_cgs.dat", status="replace")
  open(unit = 71, file="bipolar_vphi.dat", status="replace")
22  format(4(1PE15.3, 1x))
23  format(2(1PE15.3, 1x))
  do i = 1, nr
     write(66,22) r(i), density(i,1)   
     write(67,22) r(i), density(i,2)
     write(68,22) r(i), Vr(r(i)*R_mo, V0, Vinf, R_mo, beta)*1.0d-5  ! Vr in [km/s]
     r_cm = r(i)*R_mo  ! [cm]
     write(69,22) r(i), beta * (Vinf/r_cm) * (1.0d0/r(i))**2 * (1.0d0-(1.0d0/r(i)))**(beta-1.0d0)  ! [cm/s^2]
     write(70,22) r(i), density(i,1)*rho0  ! [g/cm^3]
  end do


  !
  ! Output Vr here.


contains
  
  real*8 function norm_fac(m, theta_wind)
    implicit none
    real*8, intent(in) :: m            ! [-]
    real*8, intent(in) :: theta_wind   ! [rad]
    !
    real*8, parameter :: pi =3.141592654d0
    real*8 :: a, b

    if (theta_wind ==0.0d0 .or. theta_wind == pi) then
       print *, "Error:: theta_wind ==0.0d0 .or. theta_wind == pi in norm_fac."
       stop
    else
       a = (1.0d0+m)
       b = COS(theta_wind)
       norm_fac = a/(1.0d0-b**a)
    end if

  end function norm_fac



  real*8 function Vr(r, V0, Vinf, R_mo, beta)
    implicit none
    real*8, intent(in) :: r      ! [cm]
    real*8, intent(in) :: V0           ! [cm/s]
    real*8, intent(in) :: Vinf         ! [cm/s]
    real*8, intent(in) :: R_mo         ! [R_mo]  Outer radius of the magnetosphere
    real*8, intent(in) :: beta         ! [-]

    Vr = V0 + Vinf*(1.0d0 - R_mo/r)**beta

  end function Vr

  real*8 function rho(r, theta, m, theta_wind, Mdot_wind, V0, Vinf, R_mo, beta)
    implicit none
    real*8, intent(in) :: r      ! [cm]
    real*8, intent(in)  :: theta  ! [rad]
    real*8, intent(in) :: m            ! [-]
    real*8, intent(in) :: theta_wind   ! [rad]
    real*8, intent(in) :: Mdot_wind    ! [grams/s]
    real*8, intent(in) :: V0           ! [cm/s]
    real*8, intent(in) :: Vinf         ! [cm/s]
    real*8, intent(in) :: R_mo         ! [R_mo]  Outer radius of the magnetosphere
    real*8, intent(in) :: beta         ! [-]
    !
    real*8, parameter :: pi =3.141592654d0

    rho = (norm_fac(m,theta_wind) * COS(theta)**m * Mdot_wind)   &
         /                                                       &
         (4.0d0*pi*r*r*Vr(r,V0,Vinf,R_mo,beta))
    
  end function rho



end program plot_rho
