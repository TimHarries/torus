module irreg_grid_data_class

  use kind_mod
  ! ------------------------------------------------------------
  ! Base on sph_data_class. Simplified verion of sph_data_class.
  ! The data class which holds density, velocity at irregular or 
  ! scattered points in 3D. The position (specified in cartesian
  ! coordinates) are also stored.  
  ! ------------------------------------------------------------
  ! LOGS: 
  ! Created on Apr-04-2006 (Ryuichi Kurosawa)
  ! -----------------------------------------------------------


  public:: &
       new, &
       kill, &
       read_data_romanova, &
       get_ndata, &
       get_position, &
       put_position, &
       get_rho, &
       put_rho, &
       get_min_value, &
       get_max_value, &
       info, &
       scale_xyz, &
       scale_rho, &
       scale_velocity
  
  
  private:: &
       init_irreg_grid_data, &
       kill_irreg_grid_data
  


  ! At a given time (time)
  type irreg_grid_data
!     private  ! Believe me. It's better to be private!    
     logical      :: alive   ! T if array components are allocated
     integer      :: ndata                  ! Total number data points
     ! The following arrays should be allocated with ndata
     ! Positions of data points
     real(double), pointer, dimension(:) :: x,y,z
     ! Density 
     real(double), pointer, dimension(:) :: rho
     ! Velocity components
     real(double), pointer, dimension(:) :: Vx, Vy, Vz
  end type irreg_grid_data


  !
  interface new
     module procedure init_irreg_grid_data
  end interface
  !
  interface kill
     module procedure kill_irreg_grid_data
  end interface
  !
  
  
contains  


  ! 
  ! Initializes an object with parameters (if possible).
  ! 
  subroutine init_irreg_grid_data(this,  ndata)
    implicit none
    type(irreg_grid_data), intent(inout) :: this
    integer, intent(in) :: ndata     ! Number data points

    ! save these values in this object
    this%ndata = ndata

    ! allocate arrays
    ! -- 
    if (this%alive) then 
       print *, "Error :: the data is already allocated [irregular_data_class::new]."
       stop
    else
       ALLOCATE(this%x(ndata))
       ALLOCATE(this%y(ndata))
       ALLOCATE(this%z(ndata))
       ALLOCATE(this%rho(ndata))
       ALLOCATE(this%Vx(ndata))
       ALLOCATE(this%Vy(ndata))
       ALLOCATE(this%Vz(ndata))
    end if
    
  end subroutine init_irreg_grid_data

  
  !
  !
  ! Read in the data from a file, allocate the array memory, and store the
  ! the number of gas and stars in this object.
  ! Specific to the data of romanova. You should write a similar 
  ! routine to this if you want to read other types of data. 
  !
  subroutine read_data_romanova(this, filename)
    implicit none
    type(irreg_grid_data), intent(inout) :: this
    character(LEN=*), intent(in)  :: filename
    !   
    integer, parameter  :: LUIN = 20 ! logical unit # of the data file
    integer :: nr, nphi, ntheta, ndata
    real(single), allocatable :: dummy_x(:,:,:)   ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_y(:,:,:)   ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_z(:,:,:)   ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_rho(:,:,:) ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_Vx(:,:,:)  ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_Vy(:,:,:)  ! size = (nr,nphi,ntheta)
    real(single), allocatable :: dummy_Vz(:,:,:)  ! size = (nr,nphi,ntheta)

    character(LEN=20) :: dum_a
    integer :: i,j,k, m

    open(unit=LUIN, file=TRIM(filename), status="old")
    

    ! reading headers:
    READ(LUIN,"(a)") dum_a
    READ(LUIN,"(a)") dum_a
    READ(LUIN,"(a)") dum_a
    ! reading the number of grid points
    read(LUIN,*) nr
    read(LUIN,*) nphi
    read(LUIN,*) ntheta
    
    ndata = nr*nphi*ntheta  ! total number of data

    !initializing the data arrays.
    call new(this, ndata) 


    ! reading the positions and velocity components
    
    ALLOCATE(dummy_x(nr, nphi, ntheta))
    ALLOCATE(dummy_y(nr, nphi, ntheta))
    ALLOCATE(dummy_z(nr, nphi, ntheta))
    ALLOCATE(dummy_rho(nr, nphi, ntheta))
    ALLOCATE(dummy_Vx(nr, nphi, ntheta))
    ALLOCATE(dummy_Vy(nr, nphi, ntheta))
    ALLOCATE(dummy_Vz(nr, nphi, ntheta))

    write(*,*) ' '
    write(*,*) 'Reading Romanova''s IRREG_GRID data....'
    write(*,*) ' '
    read(LUIN,*) dummy_x, dummy_y, dummy_z, dummy_rho, &
         dummy_Vx, dummy_Vy, dummy_Vz
    close(LUIN)

    ! save them this objects
    m = 0
    do k = 1, ntheta
       do j = 1, nphi
          do i = 1, nr
             m = m+1 
             ! maping to 1 D array
             this%x(m) = dummy_x(i,j,k)
             this%y(m) = dummy_y(i,j,k)
             this%z(m) = dummy_z(i,j,k)
             this%rho(m) = dummy_rho(i,j,k)
             this%Vx(m) = dummy_Vx(i,j,k)
             this%Vy(m) = dummy_Vy(i,j,k)
             this%Vz(m) = dummy_Vz(i,j,k)
          end do
       end do
    end do

    ! cleaning up
    ALLOCATE(dummy_x(nr, nphi, ntheta))
    ALLOCATE(dummy_y(nr, nphi, ntheta))
    ALLOCATE(dummy_z(nr, nphi, ntheta))
    ALLOCATE(dummy_rho(nr, nphi, ntheta))
    ALLOCATE(dummy_Vx(nr, nphi, ntheta))
    ALLOCATE(dummy_Vy(nr, nphi, ntheta))
    ALLOCATE(dummy_Vz(nr, nphi, ntheta))

    write(*,*) ' '
    write(*,*) '... done.'
    write(*,*) ' '

    
  end subroutine read_data_romanova



  !
  !
  ! accessors
  !
  
  ! returns the number of gas particles
  function get_ndata(this) RESULT(out)
    implicit none
    integer ::out 
    type(irreg_grid_data), intent(in) :: this
    out = this%ndata
  end function get_ndata
    

    
    

  !
  ! Returns the position of the i_th data point
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_position(this, i, x, y, z)
    implicit none    
    type(irreg_grid_data), intent(in) :: this
    integer, intent(in) :: i
    real(double), intent(out) :: x, y, z
    
    x = this%x(i) 
    y = this%y(i)
    z = this%z(i)
    
  end subroutine get_position

  !
  ! saves the position of the i_th data point
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine put_position(this, i, x, y, z)
    implicit none    
    type(irreg_grid_data), intent(inout) :: this
    integer, intent(in) :: i
    real(double), intent(in) :: x, y, z
    
    this%x(i) = x
    this%y(i) = y
    this%z(i) = z
    
  end subroutine put_position


  ! Returns the density of i-th data point.
  
  function get_rho(this, i) RESULT(out)
    implicit none
    real(double) :: out 
    type(irreg_grid_data), intent(in) :: this
    integer, intent(in) :: i

    out  = this%rho(i)
    
  end function get_rho


  ! Assigns the density to i-th data point
  ! i-th particle.  
  subroutine put_rho(this, i, value)
    implicit none
    type(irreg_grid_data), intent(inout) :: this
    integer, intent(in) :: i
    real(double), intent(in) :: value

    this%rho(i) = value
    
  end subroutine put_rho
   

  
  
  
  

  !
  !  destructor
  ! 
  !  Deallocates the array memories

  subroutine kill_irreg_grid_data(this)
    implicit none
    type(irreg_grid_data), intent(inout) :: this
 
    if (this%alive) then 
       ! cleaning up
       DEALLOCATE(this%x, this%y, this%z)
       DEALLOCATE(this%rho)
       DEALLOCATE(this%Vx, this%Vy, this%Vz)
       ! for safty
       NULLIFY(this%x)
       NULLIFY(this%y)
       NULLIFY(this%z)
       NULLIFY(this%rho)
       NULLIFY(this%Vx)
       NULLIFY(this%Vy)
       NULLIFY(this%Vz)

       ! now the data is dead
       this%alive = .false.
    end if

  end subroutine kill_irreg_grid_data
    



  !
  ! retuns the minimum value of double array
  ! in the coomponets of this object
  !
  function get_min_value(this, name) RESULT(out)
    implicit none
    real(double) :: out
    type(irreg_grid_data), intent(in) :: this
    character(LEN=*), intent(in)  :: name 

    select case (name)
    case ("x")
       out = MINVAL(this%x)
    case ("y")
       out = MINVAL(this%y)
    case ("z")
       out = MINVAL(this%z)
    case ("rho")
       out = MINVAL(this%rho)
    case ("Vx")
       out = MINVAL(this%Vx)
    case ("Vy")
       out = MINVAL(this%Vy)
    case ("Vz")
       out = MINVAL(this%Vz)
    case DEFAULT
       print *, "Error: Unknown name passed to irregular_data_class::get_min_value."
       stop
    end select
  end function get_min_value


  ! retuns the maximum value of double array
  ! in the coomponets of this object
  !
  function get_max_value(this, name) RESULT(out)
    implicit none
    real(double) :: out
    type(irreg_grid_data), intent(in) :: this
    character(LEN=*), intent(in)  :: name 

    select case (name)
    case ("x")
       out = MAXVAL(this%x)
    case ("y")
       out = MAXVAL(this%y)
    case ("z")
       out = MAXVAL(this%z)
    case ("rho")
       out = MAXVAL(this%rho)
    case ("Vx")
       out = MAXVAL(this%Vx)
    case ("Vy")
       out = MAXVAL(this%Vy)
    case ("Vz")
       out = MAXVAL(this%Vz)
    case DEFAULT
       print *, "Error: Unknown name passed to irregular_data_class::get_max_value."
       stop
    end select
  end function get_max_value

  

  !
  ! Prints basic infomation
  !
  ! if filename is '*' then it prints on screen.
  subroutine info(this, filename)
    implicit none
    type(irreg_grid_data), intent(in) :: this
    character(LEN=*), intent(in) :: filename
    integer :: UN
    
    if (filename(1:1) == '*') then
       UN = 6   ! prints on screen
    else
       UN = 69
       open(unit=UN, file = TRIM(filename), status = 'replace')
    end if

    
    write(UN,'(a)') ' '
    write(UN,'(a)') '######################################################'
    write(UN,'(a)') 'IRREG_GRID data info :'
    write(UN,'(a)') ' '    
    write(UN,*)     'Number of data points      : ', get_ndata(this), ' [-]'
    write(UN,'(a)') ' '
    write(UN,'(a)') '#######################################################'
    write(UN,'(a)') ' '
    write(Un,'(a)') ' '
    
    if (filename(1:1) /= '*')  close(UN)

  end subroutine info


  !
  ! Given a "scale", this routine multiplies the scaling values 
  ! to every component of x, y, and z arrays in this class. 
  !
  subroutine scale_xyz(this, scale) 
    implicit none
    type(irreg_grid_data), intent(inout) :: this
    real(double), intent(in) :: scale 
    ! 
    this%x(:) = this%x(:)*scale
    this%y(:) = this%y(:)*scale
    this%z(:) = this%z(:)*scale
  end subroutine scale_xyz

  !
  ! Given a "scale", this routine multiplies the scaling values 
  ! to every component of rho array in this class. 
  !
  subroutine scale_rho(this, scale) 
    implicit none
    type(irreg_grid_data), intent(inout) :: this
    real(double), intent(in) :: scale 
    ! 
    this%rho(:) = this%rho(:)*scale
  end subroutine scale_rho

  !
  ! Given a "scale", this routine multiplies the scaling values 
  ! to every component of Vx, Vy, and Vz arrays in this class. 
  !
  subroutine scale_velocity(this, scale) 
    implicit none
    type(irreg_grid_data), intent(inout) :: this
    real(double), intent(in) :: scale 
    ! 
    this%Vx(:) = this%Vx(:)*scale
    this%Vy(:) = this%Vy(:)*scale
    this%Vz(:) = this%Vz(:)*scale
  end subroutine scale_velocity

    
end module irreg_grid_data_class



    

    
