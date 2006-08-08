module zeus_data_class

  use kind_mod
  use vector_mod
  use utils_mod, only : hunt
  ! -----------------------------------------------------------------
  ! Class definition for the data from grid based hydrodynamics code.
  ! e.g. Zeus, and so on. 
  !------------------------------------------------------------------
  ! 
  ! Modefined : May-11-2006 Added a few routines to handle Romanova's 
  !                         MHD simulation data. (R. Kurosawa) 
  !                         
  ! Created   : Nov-04-2004 (R. Kurosawa)
  ! 

  implicit none

  public:: &
       new,   &
       kill,  &
       find_num_of_grid_points, &
       put_r, put_theta, put_phi, &
       put_density, put_temperature, put_velocity, &
       get_r, get_theta, get_phi, &       
       get_nr, get_ntheta, get_nphi, &       
       get_density, get_temperature, get_velocity, &
       interp_density, &
       interp_temperature, &
       interp_velocity, &
       read_zeus_data, &
       read_romanova_data, &
       scale_coordinates, &
       scale_density, &
       scale_temperature, &
       scale_velocity, &
       max_density

  private::  &
       init_zeus_data, &
       kill_zeus_data

  ! For one time slice!
  type zeus_data
     private
     logical :: alive    ! T if the data arrays are allocated
     integer :: nr       ! number of radial points
     integer :: nphi     ! number of azimuth angles
     integer :: ntheta   ! number of polar angles
     real(double), pointer :: r(:)      ! radis array  [R*]
     real(double), pointer :: theta(:)  ! polar angle array   [rad]
     real(double), pointer :: phi(:)    ! azimuth angle array [rad]
     ! the components of the velocity are either cartisian coordinates
     !  (Vx, Vy, Vz) or spherical coordinates for now (Vr, Vtheta, Vphi)
     character(LEN=50)     :: coordinates_v    ! Either "cartisian" or "spherical"
     type(vector), pointer :: velocity(:,:,:)  ! dimension::(nr, ntheta, nphi) [km/s]
     !
     real(double), pointer :: density(:,:,:)   ! dimension::(nr, ntheta, nphi) [g/cm^3]
     ! This is for a use later. 
     real(double), pointer :: temperature(:,:,:)  ! dimension::(nr, ntheta, nphi)
     
  end type zeus_data


  !
  ! Interfaces 
  !
  interface new
     module procedure init_zeus_data
     module procedure init_zeus_data_specified
  end interface

  interface kill
     module procedure kill_zeus_data
  end interface
  
  
contains  


  !  
  ! Find the number of grid points of this data from 
  ! the data file(s), and the results will be stored in 
  ! this object.
  subroutine find_num_of_grid_points(this)
    type(zeus_data), intent(inout) :: this
    !
    ! In future, you should find out these numbers
    ! from the data files, but for now the numbers 
    ! are manually assigned here. (RK)
    this%nr = 204
    this%ntheta = 60
    this%nphi = 60
  end subroutine find_num_of_grid_points
    
  !========================================================
  ! CONSTRUCTORS
  !========================================================

  ! 
  ! Initializes an object with parameters (if possible).
  ! 
  subroutine init_zeus_data(this)
    type(zeus_data), intent(inout) :: this
    !
    ! finding the numer of grid points in the data
    call find_num_of_grid_points(this)
    
    ! allocating memory for arrays
    ALLOCATE(this%r(this%nr),  this%theta(this%ntheta),  this%phi(this%nphi))
    ALLOCATE(this%velocity(this%nr, this%ntheta, this%nphi))
    ALLOCATE(this%density(this%nr, this%ntheta, this%nphi))
    ALLOCATE(this%temperature(this%nr, this%ntheta, this%nphi))
    
    ! initialize the values for safty
    this%r(:)=0.0; this%theta(:)=0.0; this%phi(:)=0.0
    this%density(:,:,:) = 0.0
    this%temperature(:,:,:)=0.0

    ! I am not initilaizing velocity array....

    ! indicating that data arrays have been allocated.
    this%alive = .true.
    
  end subroutine init_zeus_data

  ! 
  ! Initializes an object with numer of grid points specified.
  ! 
  subroutine init_zeus_data_specified(this, nr, ntheta, nphi)
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: nr
    integer, intent(in) :: ntheta
    integer, intent(in) :: nphi
    
    !
    ! Saving the number of grid points to this object
    this%nr = nr
    this%ntheta = ntheta
    this%nphi = nphi
    
    ! allocating memory for arrays
    ALLOCATE(this%r(this%nr),  this%theta(this%ntheta),  this%phi(this%nphi))
    ALLOCATE(this%velocity(this%nr, this%ntheta, this%nphi))
    ALLOCATE(this%density(this%nr, this%ntheta, this%nphi))
    ALLOCATE(this%temperature(this%nr, this%ntheta, this%nphi))

    
    ! initialize the values for safty
    this%r(:)=0.0; this%theta(:)=0.0; this%phi(:)=0.0
    this%density(:,:,:) = 0.0; this%temperature(:,:,:) = 0.0

    ! I am not initilaizing velocity array....

    ! indicating that data arrays have been allocated.
    this%alive = .true.
    

  end subroutine init_zeus_data_specified


  !========================================================
  ! DESTRUCTORS
  !========================================================
  
  ! deallocates memory
  subroutine kill_zeus_data(this)
    type(zeus_data),  intent(inout) :: this
    !
    DEALLOCATE(this%r, this%theta, this%phi)
    DEALLOCATE(this%velocity)
    DEALLOCATE(this%density)
    DEALLOCATE(this%temperature)

    !
    NULLIFY(this%r, this%theta, this%phi)
    NULLIFY(this%velocity)
    NULLIFY(this%density)
    NULLIFY(this%temperature)

    ! I hope all the elements of the arrays are nullify.

    ! indicating that data arrays have been deallocated.
    this%alive = .false.

  end subroutine kill_zeus_data



  !========================================================
  ! TOOLS
  !=======================================================

  ! stores values in the components
  subroutine put_r(this, i, value) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i 
    real(double), intent(in) :: value
    this%r(i) = value
  end subroutine put_r

  !
  subroutine put_theta(this, i, value) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i 
    real(double), intent(in) :: value
    this%theta(i) = value
  end subroutine put_theta
  
  !
  subroutine put_phi(this, i, value) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i 
    real(double), intent(in) :: value
    this%phi(i) = value
  end subroutine put_phi

  !
  subroutine put_velocity(this, i, j, k, vel) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i, j, k
    type(vector), intent(in) :: vel
    this%velocity(i,j,k) = vel
  end subroutine put_velocity
  
  !
  subroutine put_density(this, i, j, k, value) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i, j, k
    real(double), intent(in) :: value
    this%density(i,j,k) = value
  end subroutine put_density

  !
  subroutine put_temperature(this, i, j, k, value) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i, j, k
    real(double), intent(in) :: value
    this%temperature(i,j,k) = value
  end subroutine put_temperature


  !
  ! access stored components
  !
  real(double) function get_r(this, i) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i 
    get_r = this%r(i)
  end function get_r

  !
  real(double) function get_theta(this, i) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i 
    get_theta = this%theta(i)
  end function get_theta
  
  !
  real(double) function get_phi(this, i) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i 
    get_phi = this%phi(i)
  end function get_phi

  !
  integer function get_nr(this) 
    type(zeus_data), intent(inout) :: this
    get_nr = this%nr
  end function get_nr

  !
  integer function get_ntheta(this) 
    type(zeus_data), intent(inout) :: this
    get_ntheta = this%ntheta
  end function get_ntheta

  !
  integer function get_nphi(this) 
    type(zeus_data), intent(inout) :: this
    get_nphi = this%nphi
  end function get_nphi


  !
  type(vector) function get_velocity(this, i, j, k) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i, j, k
    get_velocity = this%velocity(i,j,k)
  end function get_velocity
  
  !
  real(double) function get_density(this, i, j, k) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i, j, k
    get_density = this%density(i,j,k)
  end function get_density
  
  !
  real(double) function get_temperature(this, i, j, k) 
    type(zeus_data), intent(inout) :: this
    integer, intent(in) :: i, j, k
    get_temperature = this%temperature(i,j,k)
  end function get_temperature

  

  !
  !  Given a position in the model space (r, theta, phi), this 
  !  function returns a value of density by interpolarting the array values
  !  stored in this object.
  ! -- Use linear interpolations in phi and theta and ln(r).
  ! -- Output is in [g/cm^3]   I suppoese.
  real(double) function interp_density(this, r, theta, phi)
    type(zeus_data), intent(in) :: this
    real(double), intent(in) :: r      ! [R*] 
    real(double), intent(in) :: theta  ! [rad] 
    real(double), intent(in) :: phi    ! [rad] 
    !
    integer :: ir, itheta, iphi  ! index
    ! intermediate values
    real(double) :: d3(2,2,2), d2(2,2), d1(2), t
!    real(double) :: ln_r1, ln_r2, ln_r
    integer :: i, j, k

    ! Finding the indecies
    ! using a routine in utils_mod module
    call hunt(this%r, this%nr, r, ir)
    call hunt(this%theta, this%ntheta, theta, itheta)
    call hunt(this%phi, this%nphi, phi, iphi)
    
    ! minor adjustments
    if (ir == this%nr) ir = ir -1
    if (itheta == this%ntheta) itheta = itheta -1
    if (iphi == this%nphi) iphi = iphi -1
    if (ir<=0) ir=1
    if (itheta<=0) itheta=1
    if (iphi<=0) iphi=1

    ! Finding the eight values on the corners 
    do k = 1, 2
       do j = 1, 2
          do i = 1, 2
             d3(i,j,k) = this%density(ir+i-1, itheta+j-1, iphi+k-1)
          end do
       end do
    end do

    ! now linear interpolate values in phi 
    do i = 1, 2
       do j = 1, 2
          t = ( d3(i,j,2) - d3(i,j,1) )  /  (this%phi(iphi+1) - this%phi(iphi))
          d2(i,j) = t*(phi - this%phi(iphi)) + d3(i,j,1)
       end do
    end do

    ! now linear interpolate values in theta
    do i = 1, 2
       t = ( d2(i,2) - d2(i,1) )  /  (this%theta(itheta+1) - this%theta(itheta))
       d1(i) = t*(theta - this%theta(itheta)) + d2(i,1)
    end do


!    ! now do linear interpolate in ln(r)
!    ln_r = log(r); ln_r1 = log(this%r(ir)); ln_r2 = log(this%r(ir+1))
!    t = ( d1(2)-d1(1) ) / ( ln_r2 - ln_r1 )

    ! now do linear interpolate in r
    t = ( d1(2)-d1(1) ) / ( this%r(ir+1) - this%r(ir) )    
    interp_density = t*(r-this%r(ir)) + d1(1)   ! [g/cm^3]

    
  end function interp_density



  !
  !  Given a position in the model space (r, theta, phi), this 
  !  function returns a value of temperature by interpolarting the array values
  !  stored in this object.
  ! -- Use linear interpolations in phi and theta and ln(r).
  ! -- Output is in [g/cm^3]   I suppoese.
  real(double) function interp_temperature(this, r, theta, phi)
    type(zeus_data), intent(in) :: this
    real(double), intent(in) :: r      ! [R*] 
    real(double), intent(in) :: theta  ! [rad] 
    real(double), intent(in) :: phi    ! [rad] 
    !
    integer :: ir, itheta, iphi  ! index
    ! intermediate values
    real(double) :: d3(2,2,2), d2(2,2), d1(2), t
!    real(double) :: ln_r1, ln_r2, ln_r
    integer :: i, j, k

    ! Finding the indecies
    ! using a routine in utils_mod module
    call hunt(this%r, this%nr, r, ir)
    call hunt(this%theta, this%ntheta, theta, itheta)
    call hunt(this%phi, this%nphi, phi, iphi)
    
    ! minor adjustments
    if (ir == this%nr) ir = ir -1
    if (itheta == this%ntheta) itheta = itheta -1
    if (iphi == this%nphi) iphi = iphi -1
    if (ir<=0) ir=1
    if (itheta<=0) itheta=1
    if (iphi<=0) iphi=1

    ! Finding the eight values on the corners 
    do k = 1, 2
       do j = 1, 2
          do i = 1, 2
             d3(i,j,k) = this%temperature(ir+i-1, itheta+j-1, iphi+k-1)
          end do
       end do
    end do

    ! now linear interpolate values in phi 
    do i = 1, 2
       do j = 1, 2
          t = ( d3(i,j,2) - d3(i,j,1) )  /  (this%phi(iphi+1) - this%phi(iphi))
          d2(i,j) = t*(phi - this%phi(iphi)) + d3(i,j,1)
       end do
    end do

    ! now linear interpolate values in theta
    do i = 1, 2
       t = ( d2(i,2) - d2(i,1) )  /  (this%theta(itheta+1) - this%theta(itheta))
       d1(i) = t*(theta - this%theta(itheta)) + d2(i,1)
    end do


!    ! now do linear interpolate in ln(r)
!    ln_r = log(r); ln_r1 = log(this%r(ir)); ln_r2 = log(this%r(ir+1))
!    t = ( d1(2)-d1(1) ) / ( ln_r2 - ln_r1 )

    ! now do linear interpolate in r
    t = ( d1(2)-d1(1) ) / ( this%r(ir+1) - this%r(ir) )    
    interp_temperature = t*(r-this%r(ir)) + d1(1)   ! [g/cm^3]

    
  end function interp_temperature







  !--------------------------------------------------------------
  !  Given a position in the model space (r, theta, phi), this 
  !  function returns a velocity by interpolarting the array values
  !  stored in this object. Do real(double) calculation internally.
  ! -- Use linear interpolations in phi and theta and ln(r).
  ! -- The output velocity is in cartisian coordinates, but the 
  !    the velocity stored in this object can be either 
  !    in spherical or cartisian coordinates. 
  !  
  ! === THE RESULT WILL BE GIVEN IN CARTISIAN COORDINATES =====
  type(vector) function interp_velocity(this, r, theta, phi)
    type(zeus_data), intent(in) :: this
    real(double), intent(in) :: r      ! [R*] 
    real(double), intent(in) :: theta  ! [rad] 
    real(double), intent(in) :: phi    ! [rad] 
    ! 
    integer :: ir, itheta, iphi  ! index
    ! intermediate values
    type(doublevector) :: d3(2,2,2), d2(2,2), d1(2), t, Vip
!    real(double) :: ln_r1, ln_r2, ln_r
    integer :: i, j, k

    ! Finding the indecies
    ! using a routine in utils_mod module
    call hunt(this%r, this%nr, r, ir)
    call hunt(this%theta, this%ntheta, theta, itheta)
    call hunt(this%phi, this%nphi, phi, iphi)
    
    ! minor adjustments
    if (ir == this%nr) ir = ir -1
    if (itheta == this%ntheta) itheta = itheta -1
    if (iphi == this%nphi) iphi = iphi -1    
    if (ir<=0) ir=1
    if (itheta<=0) itheta=1
    if (iphi<=0) iphi=1


    ! Finding the eight values on the corners 
    do k = 1, 2
       do j = 1, 2
          do i = 1, 2
             d3(i,j,k) = this%velocity(ir+i-1, itheta+j-1, iphi+k-1)
          end do
       end do
    end do

    ! now linear interpolate values in phi 
    do i = 1, 2
       do j = 1, 2
          t = ( d3(i,j,2) - d3(i,j,1) )  /  (this%phi(iphi+1) - this%phi(iphi))
          d2(i,j) = t*(phi - this%phi(iphi)) + d3(i,j,1)
       end do
    end do

    ! now linear interpolate values in theta
    do i = 1, 2
          t = ( d2(i,2) - d2(i,1) )  /  (this%theta(itheta+1) - this%theta(itheta))
          d1(i) = t*(theta - this%theta(itheta)) + d2(i,1)
    end do

!    ! now do linear interpolate in ln(r)
!    ln_r = log(r); ln_r1 = log(this%r(ir)); ln_r2 = log(this%r(ir+1))
!    t = ( d1(2)-d1(1) ) / ( ln_r2 - ln_r1 )
    ! now do linear interpolate in lr
    t = ( d1(2)-d1(1) ) / ( this%r(ir+1) - this%r(ir) )
    Vip = t*(r-this%r(ir)) + d1(1)  ! [cm/s]  

    if (this%coordinates_v == "spherical") then
       ! Now converts to the cartisian coordinates
       interp_velocity%x = Vip%x*SIN(theta)*COS(phi)  &  ! R componet
            + Vip%y*COS(theta)*COS(phi)  &               ! Theta componet
            - Vip%z*SIN(phi)                             ! phi comonent
       
       interp_velocity%y = Vip%x*SIN(theta)*SIN(phi)  &  ! R componet
            + Vip%y*COS(theta)*SIN(phi)  &               ! Theta componet
            + Vip%z*COS(phi)                             ! phi comonent
       

       interp_velocity%z = Vip%x*COS(theta)  &  ! R componet
            - Vip%y*SIN(theta)                  ! Theta componet

    elseif (this%coordinates_v == "cartisian")  then
       ! no need to do coordinate transformation
       interp_velocity = Vip
    else
       print *, "Error: Unknown coordinates specified for the velocity."
       print *, "       this%coordinates_v = ", TRIM(this%coordinates_v)
       print *, "        ... [zeus_data_class::interp_velocity]."
       stop
    end if
    
  end function interp_velocity



  !
  ! Read in the data and store them in this data structure.
  ! Run init_zeus_data or new before calling this routine!
  ! The radius, theta and phi grids are also assigned to this
  ! obejct here.
  !
  subroutine read_zeus_data(this, rgrid_filename, rho_vel_filename)
    type(zeus_data), intent(inout) :: this
    character(LEN=*),intent(in) :: rgrid_filename   ! File name containig r grid.
    character(LEN=*),intent(in) :: rho_vel_filename ! File name for velocity and density
    !
    integer, parameter :: LU_R    = 55
    integer, parameter :: LU_RHOV = 56
    integer :: i, j , k
    real(double) :: dum_d, v_r, v_theta, v_phi
    real(double) :: pi_over_2, dtheta, dphi
    !
    !
    open(unit=LU_R, file=rgrid_filename, status="old")
    open(unit=LU_RHOV, file=rho_vel_filename, status="old")
    
    ! reading the radius grid data file
!    ! The first and last 2 are values for boudary zone
!    read(LU_R,*) dum_d
!    read(LU_R,*) dum_d
    do i = 1, this%nr
       read(LU_R,*) this%r(i)   ! [R_*]
    end do
    close(LU_R)

    ! now reading the velcoity and density
    do i = 1, this%nr
       do j = 1, this%ntheta
          do k = 1, this%nphi
             read(LU_RHOV, *) v_r, v_theta, v_phi, this%density(i,j,k)
             this%velocity(i,j,k) = VECTOR(v_r, v_theta, v_phi)  ! spherical components
          end do
       end do
    end do
    close(LU_RHOV)

    !
    ! Indicate the velocity is in spherical coordinates.
    this%coordinates_v = "spherical" 
    write(*,*) " "
    write(*,*) "Note: The velocity componets of the zeus data is in spherical coordinates."
    write(*,*) "             ---  [zeus_data_clas::read_zeus_data]."
    write(*,*) " "



    ! 
    ! Sets theta and phi grids  
    ! NB: In general, this should be read from a file also, 
    !     but for now we just assign them by hand here.
    pi_over_2 = ACOS(0.0)
    dtheta = pi_over_2/dble(this%ntheta)
    do j = 1, this%ntheta
       ! staring from dtheta but not from 0!
       this%theta(j) = dtheta/2.0d0 + dble(j-1)*dtheta
    end do
 
    dphi = pi_over_2/dble(this%nphi)
    do i = 1, this%nphi
       ! staring from dtheta but not from 0!
       this%phi(i) = dphi/2.0d0 + dble(i-1)*dphi
    end do
    

  end subroutine read_zeus_data


  !-----------------------------------------------------------------
  ! Read in the data from a file, allocate the array memory, and store the
  ! the number of grid points in this object. Also sets up the 
  ! coordinates arrays. 
  !-----------------------------------------------------------------
  ! Specific to the data of romanova. You should write a similar 
  ! routine to this if you want to read other types of data. 
  !
  subroutine read_romanova_data(this, filename)
    implicit none
    type(zeus_data), intent(inout) :: this
    character(LEN=*), intent(in)  :: filename
    !   
    integer, parameter  :: LUIN = 20 ! logical unit # of the data file
    integer :: nr, nphi, ntheta, ndata
    real(single), allocatable :: dummy_x(:,:,:)   ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_y(:,:,:)   ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_z(:,:,:)   ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_rho(:,:,:) ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_T(:,:,:)   ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_Vx(:,:,:)  ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_Vy(:,:,:)  ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_Vz(:,:,:)  ! size = (nr,nphi,ntheta)

    character(LEN=20) :: dum_a
    integer :: i,j,k
    real(double) :: v_x, v_y, v_z, pi, rmax, rmin, tmp

    pi = 2.0d0*ACOS(0.0d0)

    open(unit=LUIN, file=TRIM(filename), status="old")
    

    ! reading headers:
    READ(LUIN,"(a)") dum_a
    READ(LUIN,"(a)") dum_a
    READ(LUIN,"(a)") dum_a
    ! reading the number of grid points
    read(LUIN,*) nr
    read(LUIN,*) nphi
    read(LUIN,*) ntheta
    

    !initializing the data arrays.
    call new(this, nr, ntheta, nphi) 


    ! assiging the coordonates.
    ! Units for rmax and rmin are in romanova's 
    ! code units. See her paper 2004
    rmin = 0.35d0  !
!    rmax = 8.686006d0  ! full range
    rmax = 3.0d0   ! 
    tmp = (rmax-rmin)/dble(nr-1)
    do i = 1, nr
       this%r(i) = rmin + dble(i-1)*tmp
    end do
    !
    tmp = 2.0d0*pi/dble(nphi-1)
    do i = 1, nphi
       this%phi(i) = dble(i-1)*tmp
    end do
    !
    tmp = pi/dble(ntheta-1)
    do i = 1, ntheta
       this%theta(i) = dble(i-1)*tmp
    end do


    !
    ! reading the positions and velocity components
    !
    ! NOTE: The dimensions of Romanova's arrays are (nr, nphi, ntheta) while
    !       the arrays in this objects are in different order (nr, ntheta, nphi).
    !       So, we have to rearrange the Romanova's data when storing in this object.
    
    ALLOCATE(dummy_x(nr, nphi, ntheta))
    ALLOCATE(dummy_y(nr, nphi, ntheta))
    ALLOCATE(dummy_z(nr, nphi, ntheta))
    ALLOCATE(dummy_rho(nr, nphi, ntheta))
    ALLOCATE(dummy_T(nr, nphi, ntheta))
    ALLOCATE(dummy_Vx(nr, nphi, ntheta))
    ALLOCATE(dummy_Vy(nr, nphi, ntheta))
    ALLOCATE(dummy_Vz(nr, nphi, ntheta))

    write(*,*) ' '
    write(*,*) 'Reading Romanova''s 3D data....'
    write(*,*) ' '
    read(LUIN,*) dummy_x, dummy_y, dummy_z, dummy_rho, &
         dummy_T, dummy_Vx, dummy_Vy, dummy_Vz
    close(LUIN)


    ! save them this objects
    do k = 1, ntheta
       do j = 1, nphi
          do i = 1, nr
             ! Mapping
             this%density(i,k,j) = dummy_rho(i,j,k)
             this%temperature(i,k,j) = dummy_T(i,j,k)
             V_x = dble(dummy_Vx(i,j,k))
             V_y = dble(dummy_Vy(i,j,k))
             V_z = dble(dummy_Vz(i,j,k))
             ! velocity is in cartisian coordinates!
             this%velocity(i,k,j) = VECTOR(V_x, V_y, V_z) 
          end do
       end do
    end do


    !
    ! Indicate the velocity is in cartisian coordinates.
    this%coordinates_v = "cartisian" 
    write(*,*) " "
    write(*,*) "Note: The velocity componets of the romanova data is in CARTISIAN coordinates."
    write(*,*) "             ---  [zeus_data_clas::read_romanova_data]."
    write(*,*) " "



    ! cleaning up
    DEALLOCATE(dummy_x)
    DEALLOCATE(dummy_y)
    DEALLOCATE(dummy_z)
    DEALLOCATE(dummy_rho)
    DEALLOCATE(dummy_T)
    DEALLOCATE(dummy_Vx)
    DEALLOCATE(dummy_Vy)
    DEALLOCATE(dummy_Vz)

    write(*,*) ' '
    write(*,*) '... done.'
    write(*,*) ' '

  end subroutine read_romanova_data



  !
  ! routines to scale the data arrays
  ! with a constant. 
  
  ! Scaling the compponents of the coordinates
  subroutine scale_coordinates(this, scale) 
    implicit none
    type(zeus_data), intent(inout) :: this
    real(double), intent(in) :: scale 
    
    this%r(:) = this%r(:)*scale

  end subroutine scale_coordinates

  !
  ! scaling the density
  subroutine scale_density(this, scale)
    implicit none
    type(zeus_data), intent(inout) :: this
    real(double), intent(in) :: scale 
    
    this%density(:,:,:) = this%density(:,:,:)*scale
    
  end subroutine scale_density

  !
  ! scaling the temperature
  subroutine scale_temperature(this, scale)
    implicit none
    type(zeus_data), intent(inout) :: this
    real(double), intent(in) :: scale 
    
    this%temperature(:,:,:) = this%temperature(:,:,:)*scale
    
  end subroutine scale_temperature



  !
  !scaling velocity
  subroutine scale_velocity(this, scale) 
    implicit none
    type(zeus_data), intent(inout) :: this
    real(double), intent(in) :: scale 
    real(single) :: dummy
    integer :: i, j, k
    dummy = scale
    do k = 1, this%nphi
       do j = 1, this%ntheta
          do i = 1, this%nr
             this%velocity(i,j,k) = this%velocity(i,j,k)*dummy
          end do
       end do
    end do


  end subroutine scale_velocity



  ! 
  ! Given this object, this find the max density 
  !
  real(double) function max_density(this)
    implicit none
    type(zeus_data), intent(in) :: this
    
    max_density = MAXVAL(this%density)
  end function max_density


end module zeus_data_class
    
