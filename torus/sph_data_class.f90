module sph_data_class
  ! 
  ! Class definition for Mathew's SPH data.
  ! 
  ! LOGS: 
  ! Created on Jan-10-2003 (R. Kurosawa)
  ! 


  public:: &
       read_sph_data, &
       get_udist, &
       get_umass, &
       get_utime, &
       get_npart, &
       get_nptmass, &
       get_time, &
       get_gaspartmass, &
       get_position_gas_partilce, &
       put_position_gas_partilce, &
       get_rhon, &
       put_rhon, &
       get_position_pt_mass, &
       get_pt_mass, &
       get_rhon_min, &
       get_rhon_max, &
       max_distance, &
       info


  ! At a given time (time)
  type sph_data
     private  ! Believe me. It's better to be private!    
     double precision :: udist, umass, utime    ! Units of distance, mass, time in cgs
     !                                          ! (umass is M_sol, udist=0.1 pc)
     integer          :: npart                  ! Number of gas particles
     double precision :: time                   ! Time of sph data dump (in units of utime)
     integer          :: nptmass                ! Number of stars/brown dwarfs
     double precision :: gaspartmass            ! Mass of each gas particle
     ! Positions of gas particles
     double precision, pointer, dimension(:) :: xn,yn,zn
     ! Density of the gas particles
     double precision, pointer, dimension(:) :: rhon
     ! Positions of stars
     double precision, pointer, dimension(:) :: x,y,z
     !
     double precision, pointer, dimension(:) :: ptmass ! Masses of stars
     
  end type sph_data

  
  
contains  


  ! 
  ! Initializes an object with parameters (if possible).
  ! 
  subroutine init_sph_data(this, udist, umass, utime, npart, time, nptmass, gaspartmass)
    implicit none
    type(sph_data), intent(inout) :: this
    double precision, intent(in)  :: udist, umass, utime    ! Units of distance, mass, time in cgs
    !                                                       ! (umass is M_sol, udist=0.1 pc)
    integer, intent(in)           :: npart                  ! Number of gas particles
    double precision, intent(in)  :: time                   ! Time of sph data dump (in units of utime)
    integer, intent(in)           :: nptmass                ! Number of stars/brown dwarfs
    double precision, intent(in)  :: gaspartmass            ! Mass of each gas particle

    ! save these values in this object
    this%udist = udist
    this%umass = umass
    this%utime = utime
    this%npart = npart
    this%time = time
    this%nptmass = nptmass
    this%gaspartmass =gaspartmass


    ! allocate arrays
    ! -- for gas particles
    ALLOCATE(this%xn(npart))
    ALLOCATE(this%yn(npart))
    ALLOCATE(this%zn(npart))
    ALLOCATE(this%rhon(npart))


    ! -- for star positions
    ALLOCATE(this%x(nptmass))
    ALLOCATE(this%y(nptmass))
    ALLOCATE(this%z(nptmass))
    
    ! -- for mass of stars
    ALLOCATE(this%ptmass(nptmass))
    
    
  end subroutine init_sph_data

  
  !
  !
  ! Read in the data from a file, allocate the array memory, and store the
  ! the number of gas and stars in this object.
  !
  subroutine read_sph_data(this, filename)
    implicit none
    type(sph_data), intent(inout) :: this
    character(LEN=*), intent(in)  :: filename
    !   
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
    double precision :: udist, umass, utime,  time,  gaspartmass
    integer :: npart,  nptmass
    double precision, allocatable :: dummy(:)     

!    ! for debug
!    double precision :: d, r
!    integer :: i
    
    open(unit=LUIN, file=TRIM(filename), form='unformatted')
    

    ! reading in the first line
    READ(LUIN) udist, umass, utime, npart, time, nptmass, gaspartmass


    ! initilaizing the sph_data object (allocating arrays, saving parameters and so on....)
    call init_sph_data(this, udist, umass, utime, npart, time, nptmass, gaspartmass)


    ! reading the positions  of gas particles and stars,
    ALLOCATE(dummy(npart))

    write(*,*) ' '
    write(*,*) 'Reading Matthew''s SPH data....'
    write(*,*) ' '
    READ(LUIN) this%xn
    READ(LUIN) this%yn
    READ(LUIN) this%zn

    READ(LUIN) dummy   ! Vx
    READ(LUIN) dummy   ! Vy
    READ(LUIN) dummy   ! Vz

    READ(LUIN) this%rhon
    
    READ(LUIN) this%x
    READ(LUIN) this%y
    READ(LUIN) this%z

    READ(LUIN) this%ptmass

    DEALLOCATE(dummy)

    !
    ! JUST for testting...
    ! Should be deleted here later.
    !this%npart=1000

!    !
!    ! For debug...
!    !
!    ! reassign the position of gas particles randomly
!    d=max_distance(this)*4.0d0
!    do i=1,npart
!       call random_number(r)
!       this%xn(i) = -d/2.0d0 + d*r
!       call random_number(r)
!       this%yn(i) = -d/2.0d0 + d*r
!       call random_number(r)
!       this%zn(i) = -d/2.0d0 + d*r
!    end do
    
  end subroutine read_sph_data
    

  !
  !
  ! accessors
  !
  
  ! returns program units of distance in cm 
  function get_udist(this) RESULT(out)
    implicit none
    double precision :: out
    type(sph_data), intent(in) :: this
    out = this%udist
  end function get_udist

  ! returns program units of mass in g
  function get_umass(this) RESULT(out)
    implicit none
    double precision :: out
    type(sph_data), intent(in) :: this
    out = this%umass
  end function get_umass

  ! returns program units of time in s
  function get_utime(this) RESULT(out)
    implicit none
    double precision :: out
    type(sph_data), intent(in) :: this
    out = this%utime
  end function get_utime
    

  ! returns the number of gas particles
  function get_npart(this) RESULT(out)
    implicit none
    integer ::out 
    type(sph_data), intent(in) :: this
    out = this%npart
  end function get_npart
    

  ! returns the number of point masses
  function get_nptmass(this) RESULT(out)
    implicit none
    integer :: out
    type(sph_data), intent(in) :: this
    out = this%nptmass
  end function get_nptmass
    

  ! returns the time of dump time in the unit of [utime]
  function get_time(this) RESULT(out)
    implicit none
    double precision :: out
    type(sph_data), intent(in) :: this
    out = this%time
  end function get_time
    

  ! returns the mass of each gass particle in [umass]
  function get_gaspartmass(this) RESULT(out)
    implicit none
    double precision :: out
    type(sph_data), intent(in) :: this
    out = this%gaspartmass
  end function get_gaspartmass
    

  !
  ! Returns the position of the i_th gas particle. 
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_position_gas_particle(this, i, x, y, z)
    implicit none    
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i
    double precision, intent(out) :: x, y, z
    
    x = this%xn(i) 
    y = this%yn(i)
    z = this%zn(i)
    
  end subroutine get_position_gas_particle

  !
  ! saves the position of the i_th gas particle in this object. 
  !
  ! "name" must be one of the following: "x", "y", "z"
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine put_position_gas_particle(this, i, name, value)
    implicit none    
    type(sph_data), intent(inout) :: this
    integer, intent(in) :: i
    character(LEN=*), intent(in) :: name
    double precision, intent(in) :: value
    
    select case(name)
    case ("x", "X")
       this%xn(i) = value
    case ("y", "Y") 
       this%yn(i) = value
    case ("z", "Z") 
       this%zn(i) = value
    case default
       write(*,*) "Error: Unknown name passed to sph_data_class::put_position_gas_particle."
       stop
    end select
    
  end subroutine put_position_gas_particle


  ! Returns the density of gas particle at the postion of
  ! i-th particle.
  
  function get_rhon(this, i) RESULT(out)
    implicit none
    double precision :: out 
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = this%rhon(i)
    
  end function get_rhon

  ! Assigns the density of gas particle at the postion of
  ! i-th particle.
  
  subroutine put_rhon(this, i, value)
    implicit none
    type(sph_data), intent(inout) :: this
    integer, intent(in) :: i
    double precision, intent(in) :: value

    this%rhon(i) = value
    
  end subroutine put_rhon
   

  
  !
  ! Rerurns the postions of the i-th point mass.
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_position_pt_mass(this, i, x, y, z)
    implicit none    
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i
    double precision, intent(out) :: x, y, z
    
    x = this%x(i) 
    y = this%y(i)
    z = this%z(i)
    
  end subroutine get_position_pt_mass
 

  !
  ! Returns the mass of the point mass (star) in
  ! the i-th position of data.
  !
  function get_pt_mass(this, i) RESULT(out)
    implicit none
    double precision  :: out 
    type(sph_data), intent(in) :: this
    integer, intent(in):: i

    out = this%ptmass(i)  !in program units [umass]. See above.

  end function get_pt_mass
  
  
  

  !
  !  destructor
  ! 
  !  Deallocates the array memories

  subroutine kill(this)
    implicit none
    type(sph_data), intent(inout) :: this
    
    DEALLOCATE(this%xn, this%yn, this%zn)
    DEALLOCATE(this%rhon)
    DEALLOCATE(this%x, this%y, this%z)
    DEALLOCATE(this%ptmass)

    NULLIFY(this%xn, this%yn, this%zn)
    NULLIFY(this%x, this%y, this%z)
    NULLIFY(this%ptmass)

  end subroutine kill
    




  !
  ! find the maximum distance between the pt masses
  !

  function max_distance(this) RESULT(out)
    implicit none
    double precision :: out
    type(sph_data), intent(in) :: this
    !
    double precision :: d_max, d, x, y, z
    integer :: i, j, n
    
    d_max= -1.0
    n=get_nptmass(this) ! function in this moudle
    
    do i = 1, n-2
       do j = i+1, n
          call get_position_pt_mass(this, j, x, y, z)
          d = (x*x+y*y+z*z)  ! omit SQRT here cus it costs too much.
          d_max = MAX(d_max, d)
       end do
    end do 
    d_max = SQRT(d_max)
    out = d_max
  end function max_distance


  !
  ! retuns the minimum value of rhon
  !
  function get_rhon_min(this) RESULT(out)
    implicit none
    double precision :: out
    type(sph_data), intent(in) :: this

    out = MINVAL(this%rhon)
    
  end function get_rhon_min

  
  !
  ! retuns the maximum value of rhon
  !
  function get_rhon_max(this) RESULT(out)
    implicit none
    double precision :: out
    type(sph_data), intent(in) :: this

    out = MAXVAL(this%rhon)
    
  end function get_rhon_max


  !
  ! Prints basic infomation
  !
  ! if filename is '*' then it prints on screen.
  subroutine info(this, filename)
    implicit none
    type(sph_data), intent(in) :: this
    character(LEN=*), intent(in) :: filename
    integer :: UN
    double precision :: tmp
    
    if (filename(1:1) == '*') then
       UN = 6   ! prints on screen
    else
       UN = 69
       open(unit=UN, file = TRIM(filename), status = 'replace')
    end if

    ! time of the data dump
    tmp = get_time(this)*get_utime(this)/(60.0d0*60.0d0*24.0d0*365.0d0*1.0d6)
    
    write(UN,'(a)') ' '
    write(UN,'(a)') '######################################################'
    write(UN,'(a)') 'SPH data info :'
    write(UN,'(a)') ' '    
    write(UN,*)     'Units of length         : ', get_udist(this), ' [cm]'
    write(UN,*)     'Units of mass           : ', get_umass(this), ' [g]'
    write(UN,*)     'Units of time           : ', get_utime(this),  ' [s]' 
    write(UN,'(a)') ' '    
    write(UN,*)     'Number of stars         : ',  get_nptmass(this)
    write(UN,*)     'Number of gas particles : ',  get_npart(this)   
    write(UN,*)     'Time of data dump       : ',  tmp, ' [Myr]'    
    write(UN,*)     'Mass (gas particle)     : ',  get_gaspartmass(this)*get_umass(this), ' [g]'    
    write(UN,'(a)') '#######################################################'
    write(UN,'(a)') ' '
    write(Un,'(a)') ' '
    
    if (filename(1:1) /= '*')  close(UN)

  end subroutine info
    
end module sph_data_class



    

    
