program plot_discwind

  use subs

  implicit none

  real*8, parameter :: pi =3.141592654d0
  ! Parameters 
  real*8 :: Rmin, Rmax, Tmax, alpha, gamma, d, Mstar

  real*8 :: R_min, R_max
  integer, parameter :: nl = 200
  integer, parameter :: nbeta = 4
  real*8 :: beta(nbeta), rho0(nbeta)
  real*8 :: l(nl), lmin, lmax, log10_lmin, log10_lmax, log10_dl
  real*8 :: density(nl,nbeta)
  integer :: i, j
  real*8 :: wi, f, AU

  ! setting up the parameters====================
  Rmin = 43.0d0  ! [10^10 cm]
  Rmax = 1500d0  ! [10^10 cm]
  Tmax = 1600d0  ! [K]
  alpha = 0.33   ! [-]
  gamma = -1.15  ! [-]
  Mstar = 2.0d0  ! [Msun]
  d = 300.d0     ! [10^10cm]
  f = 2.0d0
  AU = 1.5e3     ! [10^10 cm]  1 AU
  !==============================================

  ! set up L values
!  lmin = 1.0d-3           ! [10^10 cm]
  lmin = 1.0d-6           ! [10^10 cm]
!  lmax = 36.0d0*1000.d0  ! [10^10 cm]
  lmax = 150000.d0  ! [10^10 cm]
  log10_lmin = LOG10(lmin); log10_lmax = LOG10(lmax)
  log10_dl = (log10_lmax-log10_lmin)/dble(nl-1)

  do i = 1, nl
     l(i) = log10_lmin + log10_dl*dble(i-1)
     l(i) = 10.0d0**l(i)
  end do

  ! set up the values of beta
  beta(1) = 0.2
  beta(2) = 0.5
  beta(3) = 1.0
  beta(4) = 2.0


  ! now compute the density
  wi = (Rmin + Rmax)/2.0d0  ! mid point on the disc
  
  do j = 1, nbeta
     do i = 1, nl
        density(i,j) = discwind_density(l(i), wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta(j), Mstar)
     end do
  end do


  ! normalizing density
  do j = 1, nbeta     
     rho0(j) = discwind_density(0.0d0, wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta(3), Mstar)
  end do


  ! now computes the velocity components:
  


  ! wrting the results to files
  open(unit = 55, file ="dw_rho.dat", status="replace")
  open(unit = 56, file ="dw_vp.dat", status="replace")
  open(unit = 57, file ="dw_vphi.dat", status="replace")
22  format(5(1PE15.3, 1x))
24  format(2(1PE15.3, 1x))
  wi = (Rmin + Rmax)/2.0d0  ! mid point on the disc
  do i = 1, nl
     write(55,22) l(i)/AU, density(i,1)/rho0(1), density(i,2)/rho0(2), &
          density(i,3)/rho0(3), density(i,4)/rho0(4)
     write(56,22) l(i)/AU, &
          Vp(l(i), wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta(1), Mstar)*1.0d-5, &
          Vp(l(i), wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta(2), Mstar)*1.0d-5, &
          Vp(l(i), wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta(3), Mstar)*1.0d-5, &
          Vp(l(i), wi, f, Rmin, Rmax, Tmax, alpha, gamma, d, beta(4), Mstar)*1.0d-5      ! [km/s]
     write(57,24) l(i)/AU,  Vphi(l(i), wi, d,  Mstar)*1.0d-5  ! [km/s]
  end do

  close(55)
  close(56)
  close(57)



end program plot_discwind
