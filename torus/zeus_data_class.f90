module zeus_data_class

  use kind_mod
  use vector_mod
  use utils_mod, only : hunt
  ! 
  ! Class definition for the data from grid based hydrodynamics code.
  ! e.g. Zeus, and so on. 
  !
  ! Created on Nov-04-2004 (R. Kurosawa)
  ! 

  implicit none

  public:: &
       new,   &
       kill,  &
       find_num_of_grid_points, &
       put_r, put_theta, put_phi, &
       put_density, put_velocity, &
       get_r, get_theta, get_phi, &       
       get_nr, get_ntheta, get_nphi, &       
       get_density, get_velocity, &
       interp_density, &
       interp_velocity

  private::  &
       init_zeus_data, &
       kill_zeus_data

  ! For one time slice!
  type zeus_data
     private
     integer :: nr       ! number of radial points
     integer :: nphi     ! number of azimuth angles
     integer :: ntheta   ! number of polar angles
     real(double), pointer :: r(:)      ! radis array  [R*]
     real(double), pointer :: theta(:)  ! polar angle array   [rad]
     real(double), pointer :: phi(:)    ! azimuth angle array [rad]
     ! the components of the velocity are in spherical coordinates for now (Vr, Vtheta, Vphi)
     type(vector), pointer :: velocity(:,:,:)  ! dimension::(nr, ntheta, nphi) [km/s]
     real(double), pointer :: density(:,:,:)   ! dimension::(nr, ntheta, nphi) [g/cm^3]
     ! This is for a use later. 
     !real(double), pointer :: temperature(:,:,:)   ! dimension::(nr, ntheta, nphi)
  end type zeus_data


  !
  ! Interfaces 
  !
  interface new
     module procedure init_zeus_data
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
    
    ! initialize the values for safty
    this%r(:)=0.0; this%theta(:)=0.0; this%phi(:)=0.0
    this%density(:,:,:) = 0.0
    ! I am not initilaizing velocity array....

  end subroutine init_zeus_data


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
    !
    NULLIFY(this%r, this%theta, this%phi)
    NULLIFY(this%velocity)
    NULLIFY(this%density)
    ! I hope all the elements of the arrays are nullify.
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
  !  function returns a velocity by interpolarting the array values
  !  stored in this object. Do real(double) calculation internally.
  ! -- Use linear interpolations in phi and theta and ln(r).
  ! -- The output is in [cm/s] .... I think.
  ! === THE RESULT WILL BE GIVEN IN CARTISIAN COORDINATES =====
  type(vector) function interp_velocity(this, r, theta, phi)
    type(zeus_data), intent(in) :: this
    real(double), intent(in) :: r      ! [R*] 
    real(double), intent(in) :: theta  ! [rad] 
    real(double), intent(in) :: phi    ! [rad] 
    ! 
    integer :: ir, itheta, iphi  ! index
    ! intermediate values
    type(doublevector) :: d3(2,2,2), d2(2,2), d1(2), t, Vsp
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
    Vsp = t*(r-this%r(ir)) + d1(1)  ! [cm/s]  in spherical coordinates

    ! Now converts to the cartisian coordinates
    interp_velocity%x = Vsp%x*SIN(theta)*COS(phi)  &  ! R componet
         + Vsp%y*COS(theta)*COS(phi)  &               ! Theta componet
         - Vsp%z*SIN(phi)                             ! phi comonent
         
    interp_velocity%y = Vsp%x*SIN(theta)*SIN(phi)  &  ! R componet
         + Vsp%y*COS(theta)*SIN(phi)  &               ! Theta componet
         + Vsp%z*COS(phi)                             ! phi comonent
         

    interp_velocity%z = Vsp%x*COS(theta)  &  ! R componet
         - Vsp%y*SIN(theta)                  ! Theta componet
         

    
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



end module zeus_data_class
    
