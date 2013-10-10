! A test of SPH-Torus which uses particles to represent the benchmark disc 

program sphbench

!use torus_mod, only: torus
use particle_pos_mod, only: particle_pos

  implicit none

! Total number of particles
  integer, parameter :: npart=1e5

! loop and particle index
  integer :: ipart

! Define double precision 
  integer, parameter :: db = selected_real_kind(15,307)

! Cylindrical polar co-ordinates
  real(db), allocatable :: r(:), z(:)
  real(db) :: theta

! Cartesian co-ordinates used with z above
  real(db) :: x, y 

! Torus arguments
  integer            :: b_idim              ! Maximum array size
  integer            :: b_npart
  integer            :: b_nptmass           ! Number of point masses
  integer            :: b_num_gas           ! Number of gas particles
  real(kind=8),allocatable     :: b_xyzmh(:,:)
  real(kind=4),allocatable     :: b_rho(:)
  integer(kind=1), allocatable :: b_iphase(:) 
  real(kind=8), parameter :: b_udist=1.0 ! unit of distance
  real(kind=8), parameter :: b_umass=1.0 ! units of mass
  real(kind=8), parameter :: b_utime=1.0 ! Units of time
  real(kind=8), parameter :: year = 60.0*60.0*24.0*365.25 
  real(kind=8), parameter :: b_time=182.0e6_db * year ! Current time, used as age of star
  real(kind=8), allocatable :: b_temp(:)
  real(kind=8), parameter :: total_disc_mass=0.011  ! Taken from 2D benchmark
  real(kind=8)            :: total_gas_mass
  character(len=11), parameter :: file_tag = "sphbench   "

! Source parameters
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

! Smoothing length factor. See Price and Bate, 2007
  real, parameter :: eta_smooth = 1.2

  real(db) :: f1, f2, hr, z_over_h

  real :: ran_num

! Begin executable statments -------------

! 1. Set up gas particles 

! Set up the random number generator  
  call random_seed

  total_gas_mass = total_disc_mass
  b_num_gas      = npart

  write(*,*) "SPH benchmark: generating", npart, "particles"
  b_idim = b_num_gas + 1

  allocate ( r(npart) )
  allocate ( z(npart) )

  call particle_pos( npart, r, z)

  r(:) = r(:) * auToCm
  z(:) = z(:) * auToCm

  allocate ( b_xyzmh(5,b_idim) )
  allocate ( b_rho(b_idim)     )
  allocate ( b_iphase(b_idim)  )
  allocate ( b_temp(b_num_gas) )

! Set up gas particle information
  part_loop:  do ipart=1, b_num_gas

     call random_number(ran_num)
     theta = ran_num * 2.0_db * pi
           
     if ( r(ipart) > disc_r_outer .or. r(ipart) < disc_r_inner ) then
        b_rho(ipart) = rho_bg
     else

        hr = z_d * ( ( r(ipart) / r_d ) ** 1.125 )
        z_over_h = z(ipart) / hr
        f1 = ( r(ipart) / r_d ) ** (-1) 
        f2 = exp ( minusPiByFour * ( z_over_h **2 ) )

        b_rho(ipart) = f1 * f2 * rho_zero
        if (b_rho(ipart) < rho_bg) b_rho(ipart) = rho_bg

     endif

! Set positions of the gas particles
     x = r(ipart) * sin(theta) 
     y = r(ipart) * cos(theta) 
     b_xyzmh(1,ipart) = x
     b_xyzmh(2,ipart) = y
     b_xyzmh(3,ipart) = z(ipart) 
! Set particle mass assuming equal mass for all particles 
     b_xyzmh(4,ipart) = msol * total_gas_mass / real(b_num_gas, kind=db)
! Smoothing length based on particle mass and density
     b_xyzmh(5,ipart) = eta_smooth * ( b_xyzmh(4,ipart) / b_rho(ipart) ) ** (1.0/3.0)

  end do part_loop

  deallocate ( r, z )

! Initialise phase flag. Gas particles are denoted by zero.
   b_iphase(1:b_num_gas) = 0 
   b_temp(1:b_num_gas) = 3.0

! 2. Set up point mass properties

! Set properties of the point mass
   b_npart                = b_num_gas + 1
   b_nptmass              = 1
   b_iphase(b_num_gas+1)  = 1
   b_xyzmh(1,b_num_gas+1) = 0.0
   b_xyzmh(2,b_num_gas+1) = 0.0
   b_xyzmh(3,b_num_gas+1) = 0.0
   b_xyzmh(4,b_num_gas+1) = source_mass

  open(unit=61, file="part_out.dat", status='replace')
! Gas particles
  do ipart=1, b_num_gas
     write(61,'(1x,10(e14.7,2x),i4)')  b_xyzmh(:,ipart), b_rho(ipart), b_temp(ipart), 0.0, 0.0, 0.0, 1
  end do
! Sink particle
   write(61,'(1x,10(e14.7,2x),i4)') b_xyzmh(:,b_num_gas+1), 1e-20, 1000.0, 0.0, 0.0, 0.0, 3
  close(61)

  deallocate ( b_xyzmh  )
  deallocate ( b_rho    )
  deallocate ( b_iphase )
  deallocate ( b_temp   ) 

end program sphbench

