! A test of SPH-Torus which uses particles to represent the benchmark disc 

program sphbench

use torus_mod, only: torus

  implicit none

! number of particle spaces
  integer, parameter :: npart=1e6

! loop and particle index
  integer :: ipart

! Define double precision 
  integer, parameter :: db = selected_real_kind(15,307)

! Cylindrical polar co-ordinates
  real(db) :: r, theta, z

! Cartesian co-ordinates used with z above
  real(db) :: x, y 

! Torus arguments
  integer, parameter :: b_idim  = npart + 1 ! Set maximum array sizes
  integer, parameter :: b_npart = npart + 1 ! Gas particles + one point mass
  integer, parameter :: b_nactive=-99       ! Not used
  integer, parameter :: b_nptmass = 1       ! Number of point masses
  integer, parameter :: b_num_gas = npart   ! Number of gas particles
  real(kind=8) :: b_xyzmh(5,b_idim)
  real(kind=4) :: b_rho(b_idim)
  integer(kind=1) :: b_iphase(b_idim) 
  real(kind=8), parameter :: b_udist=1.0 ! unit of distance
  real(kind=8), parameter :: b_umass=1.0 ! units of mass
  real(kind=8), parameter :: b_utime=1.0 ! Units of time
  real(kind=8), parameter :: year = 60.0*60.0*24.0*365.25 
  real(kind=8), parameter :: b_time=182.0e6_db * year ! Current time, used as age of star
  real(kind=8) :: b_temp(b_num_gas)
  real(kind=8), parameter :: temp_min=3.0
  real(kind=8), parameter :: total_gas_mass=0.011  ! Taken from 2D benchmark
  character(len=11), parameter :: file_tag = "sphbench   "

! Source parameters
  real(db), parameter :: source_x = 0.0
  real(db), parameter :: source_y = 0.0
  real(db), parameter :: source_z = 0.0

  real(db), parameter :: mSol = 1.9891e33_db
  real(db), parameter :: source_mass = 1.0 * mSol

! disc parameters
  real(db), parameter :: pi = 3.1415926535897932_db
  real(db), parameter :: minusPiByFour = ( pi / 4.0_db) * (-1.0_db)
  real(db), parameter :: auToCm = 1.495979d13
  real(db), parameter :: disc_r_outer = 1000.0 * auToCm
  real(db), parameter :: disc_r_inner = 1.0 * auToCm
  real(db), parameter :: r_d = 0.5 * disc_r_outer
  real(db), parameter :: z_d = 0.25 * r_d
  real(db), parameter :: rho_zero = 0.81614E-17_db
  real(db), parameter :: rho_bg   = 1.001e-30_db ! Background density

  real(db) :: f1, f2, hr, z_over_h

  real :: ran_num

! Begin executable statments -------------

! 1. Set up gas particles 

! Set up gas particle information
  part_loop:  do ipart=1, npart

     call random_number(ran_num)
     r     = disc_r_inner * exp ( log (disc_r_outer/disc_r_inner) * ran_num ) 
     call random_number(ran_num)
     theta = ran_num * 2.0_db * pi
     call random_number(ran_num)
     z     =  (ran_num-0.5) * 2.0_db * disc_r_outer 
           
     if ( r > disc_r_outer .or. r < disc_r_inner ) then
        b_rho(ipart) = rho_bg
     else

        hr = z_d * ( ( r / r_d ) ** 1.125 )
        z_over_h = z / hr
        f1 = ( r / r_d ) ** (-1) 
        f2 = exp ( minusPiByFour * ( z_over_h **2 ) )

        b_rho(ipart) = f1 * f2 * rho_zero

     endif

! Set positions of the gas particles
     x = r * sin(theta) 
     y = r * cos(theta) 
     b_xyzmh(1,ipart) = x
     b_xyzmh(2,ipart) = y
     b_xyzmh(3,ipart) = z

  end do part_loop

  where ( b_rho < rho_bg ) 
     b_rho = rho_bg
  end where


   b_xyzmh(4,1:npart) = 0.0
   b_xyzmh(5,1:npart) = 0.0      

! Initialise phase flag. Gas particles are denoted by zero.
   b_iphase(1:npart) = 0 

   b_temp(1:npart) = temp_min

! 2. Set up point mass properties

! Set properties of the point mass
   b_iphase(npart+1)  = 1
   b_xyzmh(1,npart+1) = 0.0
   b_xyzmh(2,npart+1) = 0.0
   b_xyzmh(3,npart+1) = 0.0
   b_xyzmh(4,npart+1) = source_mass

! 3. Call torus

! Call torus
  call torus(b_idim,  b_npart,       b_nactive, b_nptmass, b_num_gas, &
             b_xyzmh, b_rho,         b_iphase,                        &
             b_udist, b_umass,       b_utime,   b_time,    b_temp,    &
            temp_min, total_gas_mass, file_tag )


end program sphbench
