module sph_data_class

  use kind_mod
  use vector_mod
  use messages_mod
  ! 
  ! Class definition for Mathew's SPH data.
  ! 
  ! LOGS: 
  ! Created on Jan-10-2003 (R. Kurosawa)
  ! 

  implicit none

  public:: &
       kill, &
       read_sph_data, &
       get_udist, &
       get_umass, &
       get_utime, &
       get_npart, &
       get_nptmass, &
       get_time, &
       get_position_gas_particle, &
       put_position_gas_particle, &
       get_rhon, &
       put_rhon, &
       get_position_pt_mass, &
       get_pt_mass, &
       get_rhon_min, &
       get_rhon_max, &
       get_spins, &
       max_distance, &
       info, &
       get_stellar_disc_parameters, &
       stellar_disc_exists, &
       find_inclinations, &
       isAlive
  
    real(double), parameter :: sph_rho_min=1.0e-25_db ! Lower limit for rho plots from sphtorus
    real(double), parameter :: sph_rho_max=1.0e-12_db ! Upper limit for rho plots from sphtorus
    real(double), parameter :: sph_tem_min=10.0_db ! Lower limit for temperature plots from sphtorus
    real(double), parameter :: sph_tem_max=3.0e2_db ! Upper limit for temperature plots from sphtorus

  private:: &
       kill_sph_data
  


  ! At a given time (time)
  type sph_data
!     private  ! Believe me. It's better to be private!    
     logical      :: inUse=.false.          ! Flag to indicate if this object is in use.
     real(double) :: udist, umass, utime    ! Units of distance, mass, time in cgs
     !                                          ! (umass is M_sol, udist=0.1 pc)
     integer          :: npart                  ! Total number of gas particles (field+disc)
     real(double) :: time                   ! Time of sph data dump (in units of utime)
     integer          :: nptmass                ! Number of stars/brown dwarfs
     real(double), pointer, dimension(:) :: gasmass            ! Mass of each gas particle ! DAR changed to allow variable mass
     ! Positions of gas particles
     real(double), pointer, dimension(:) :: xn,yn,zn
     real(double), pointer, dimension(:) :: vxn,vyn,vzn
     ! Density of the gas particles
     real(double), pointer, dimension(:) :: rhon
     ! Temperature of the gas particles
     real(double), pointer, dimension(:) :: temperature
     ! Positions of stars
     real(double), pointer, dimension(:) :: x,y,z
     real(double), pointer, dimension(:) :: vx,vy,vz
     !
     real(double), pointer, dimension(:) :: ptmass ! Masses of stars
     ! 
     !
     ! Some extra data for the stellar disk.
     ! Note: the units used here are diffrent from the ones used above!     
     logical :: have_stellar_disc            ! T if the following data are assigned
     real(double), pointer, dimension(:) :: discrad ! in [10^10cm]
     real(double), pointer, dimension(:) :: discmass   ! in [g]
     real(double), pointer, dimension(:) :: spinx, spiny, spinz
          
  end type sph_data


  !
  !
  interface kill
     module procedure kill_sph_data
  end interface
  
  
contains  


  ! 
  ! Initializes an object with parameters (if possible).
  ! 
  subroutine init_sph_data(sphdata, udist, umass, utime, npart,  time, nptmass)
    implicit none
    type(sph_data), intent(inout) :: sphdata
    real(double), intent(in)  :: udist, umass, utime    ! Units of distance, mass, time in cgs
    !                                                       ! (umass is M_sol, udist=0.1 pc)
    integer, intent(in)           :: npart                  ! Number of gas particles (field+disc)
    real(double), intent(in)  :: time                   ! Time of sph data dump (in units of utime)
    integer, intent(in)           :: nptmass                ! Number of stars/brown dwarfs

    ! Indicate that this object is in use
    sphdata%inUse = .true.

    ! save these values in this object
    sphdata%udist = udist
    sphdata%umass = umass
    sphdata%utime = utime
    sphdata%npart = npart
    sphdata%time = time
    sphdata%nptmass = nptmass


    ! allocate arrays
    ! -- for gas particles
    ALLOCATE(sphdata%xn(npart))
    ALLOCATE(sphdata%yn(npart))
    ALLOCATE(sphdata%zn(npart))


    ALLOCATE(sphdata%vxn(npart))
    ALLOCATE(sphdata%vyn(npart))
    ALLOCATE(sphdata%vzn(npart))

    ALLOCATE(sphdata%rhon(npart))
    ALLOCATE(sphdata%temperature(npart))
    ALLOCATE(sphdata%gasmass(npart))

    ! -- for star positions
    ALLOCATE(sphdata%x(nptmass))
    ALLOCATE(sphdata%y(nptmass))
    ALLOCATE(sphdata%z(nptmass))

    ALLOCATE(sphdata%vx(nptmass))
    ALLOCATE(sphdata%vy(nptmass))
    ALLOCATE(sphdata%vz(nptmass))
    
    ! -- for mass of stars
    ALLOCATE(sphdata%ptmass(nptmass))

    !
    sphdata%have_stellar_disc = .false. 

  end subroutine init_sph_data

  ! 
  ! Initializes an object with parameters when torus is called as a subroutine from sphNG.
  ! 
  subroutine init_sph_data2(this, udist, umass, utime, npart,  time, nptmass, &
        b_npart, b_idim, b_iphase, b_xyzmh, b_rho, b_temp)
    implicit none

! Arguments --------------------------------------------------------------------
    type(sph_data), intent(inout) :: this
    real(double), intent(in)  :: udist, umass, utime    ! Units of distance, mass, time in cgs
    !                                                   ! (umass is M_sol, udist=0.1 pc)
    integer, intent(in)           :: npart              ! Number of gas particles (field+disc)
    real(double), intent(in)  :: time                   ! Time of sph data dump (in units of utime)
    integer, intent(in)           :: nptmass            ! Number of stars/brown dwarfs

    integer, intent(in)   :: b_npart, b_idim
    integer*1, intent(in) :: b_iphase(b_idim)
    real*8, intent(in)    :: b_xyzmh(5,b_idim)
    real*4, intent(in)    :: b_rho(b_idim)
    real*8, intent(in)    :: b_temp(b_idim)
! Local variables --------------------------------------------------------------
    integer :: iii, iiipart, iiigas

! Begin executable statements --------------------------------------------------

    ! Indicate that this object is in use
    this%inUse = .true.

    ! save these values in this object
    this%udist = udist
    this%umass = umass
    this%utime = utime
    this%npart = npart
    this%time = time
    this%nptmass = nptmass
    this%nptmass     = nptmass

    ! allocate arrays
    ! -- for gas particles
    ALLOCATE(this%xn(npart))
    ALLOCATE(this%yn(npart))
    ALLOCATE(this%zn(npart))
    ALLOCATE(this%rhon(npart))
    ALLOCATE(this%temperature(npart))


    ! -- for star positions
    ALLOCATE(this%x(nptmass))
    ALLOCATE(this%y(nptmass))
    ALLOCATE(this%z(nptmass))
    
    ! -- for mass of stars
    ALLOCATE(this%ptmass(nptmass))

    !
    this%have_stellar_disc = .false. 
    
    iiipart=0
    iiigas=0

    ! Set up density and position of gas and point mass particles
    do iii=1,b_npart
       if (b_iphase(iii) == 0) then
          iiigas=iiigas+1
          if (iiigas > npart) then
             write (*,*) 'Error: iiigas>npart',iiigas, npart
             cycle
          endif
          this%rhon(iiigas)        = b_rho(iii)
          this%xn(iiigas)          = b_xyzmh(1,iii)
          this%yn(iiigas)          = b_xyzmh(2,iii)
          this%zn(iiigas)          = b_xyzmh(3,iii)
! b_temp contains gas particle data only
          this%temperature(iiigas) = b_temp(iiigas)
       elseif (b_iphase(iii) > 0) then
          iiipart=iiipart+1
          if (iiipart > nptmass) then
             write (*,*) 'Error: iiipart>nptmass',iiipart,nptmass
             cycle
          endif
          this%x(iiipart)      = b_xyzmh(1,iii)
          this%y(iiipart)      = b_xyzmh(2,iii)
          this%z(iiipart)      = b_xyzmh(3,iii)
          this%ptmass(iiipart) = b_xyzmh(4,iii)
       endif
    end do

    ! Check the number of gas particles and point masses are as expected
    IF ( iiigas /= npart ) THEN
       print *,'Error: expecting to find', npart, 'gas particles but found', iiigas
       do; end do
    ENDIF
    IF ( iiipart /= nptmass ) THEN
       print *,'Error: expecting to find', nptmass, 'point masses but found', iiipart
       do; end do 
    ENDIF

  end subroutine init_sph_data2

  
  !
  !
  ! Read in the data from a file, allocate the array memory, and store the
  ! the number of gas and stars in this object.
  !
  ! Note: This routine only works with the 'old' (unknown how old) dump file
  ! format. Use new_read_sph_data to read in ASCII created by SPLASH.

  subroutine read_sph_data(this, filename)
    implicit none
    type(sph_data), intent(inout) :: this
    character(LEN=*), intent(in)  :: filename
    !   
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
!    real(double) :: udist, umass, utime,  time,  gaspartmass, discpartmass
!    integer*4 :: npart,  nsph, nptmass
    real(double) :: udist, umass, utime
    real(double) :: gaspartmass, time
    integer :: npart, nptmass, n1, n2
    real(double), allocatable :: dummy(:)     


    open(unit=LUIN, file=TRIM(filename), form='unformatted')
  

    ! reading in the first line
    READ(LUIN) udist, umass, utime, npart, n1, n2, time, nptmass, gaspartmass

    ! initilaizing the sph_data object (allocating arrays, saving parameters and so on....)
    call init_sph_data(this, udist, umass, utime, npart, time, nptmass)


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

  end subroutine read_sph_data

  subroutine new_read_sph_data(sphdata, filename)
    implicit none
    type(sph_data), intent(inout) :: sphdata
    character(LEN=*), intent(in)  :: filename
    !   
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
!    real(double) :: udist, umass, utime,  time,  gaspartmass, discpartmass
!    integer*4 :: npart,  nsph, nptmass
    real(double) :: udist, umass, utime,  time
    real(double) :: xn, yn, zn, vx, vy, vz, gaspartmass, rhon, masscounter, u
    integer :: itype, ipart, icount, iptmass, igas, idead
    integer :: npart, nptmass, n1, n2
    real(double) junk
    character(LEN=1)  :: junkchar
    character(LEN=100) :: message

    open(unit=LUIN, file=TRIM(filename), form="formatted")

    read(LUIN,*) 
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*) junkchar, time, utime
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*) junkchar, npart, n1, nptmass, n2
    npart = npart + n1 ! no. of particles we're interested in
    read(LUIN,*)
    read(LUIN,*) junkchar, udist, junk, junk, umass
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*)

    call init_sph_data(sphdata, udist, umass, utime, npart, time, nptmass)

    npart = npart + n2 + nptmass ! npart now equal to no. lines - 12 = sum of particles dead or alive

    write(message,*) "Reading SPH data from ASCII...."
    call writeinfo(message, TRIVIAL)

    iptmass = 0
    icount = 12
    igas = 0
    idead = 0
    masscounter = 0.d0

    do ipart=1, npart

       read(LUIN,*) xn, yn, zn, gaspartmass, junk, rhon, vx, vy, vz, u, junk, junk, junk, itype
       icount = icount + 1

       if(mod(icount,10000) .eq. 1) then
          write(message,*) "Reading", icount,"th particle"
          call writeinfo(message, TRIVIAL)
       endif

       if(itype .eq. 1) then ! .or. itype .eq. 4) then
          igas = igas + 1

          sphdata%xn(igas) = xn
          sphdata%yn(igas) = yn
          sphdata%zn(igas) = zn

          sphdata%gasmass(igas) = gaspartmass
          sphdata%rhon(igas) = rhon

          sphdata%vxn(igas) = vx
          sphdata%vyn(igas) = vy
          sphdata%vzn(igas) = vz

!         sphdata%temperature = 2. * 2.46 * (u * 1d-7) / (3. * 8.314472) ! 8.31 is gas constant
          sphdata%temperature = 1.9725e-8 * u
          
          masscounter = masscounter + gaspartmass
          
       Elseif(itype .eq. 3) then

          iptmass = iptmass+1

          sphdata%x(iptmass) = xn
          sphdata%y(iptmass) = yn
          sphdata%z(iptmass) = zn

          sphdata%vx(iptmass) = vx
          sphdata%vy(iptmass) = vy
          sphdata%vz(iptmass) = vz

          sphdata%ptmass(iptmass) = gaspartmass 

          masscounter = masscounter + gaspartmass

          write(message,*) "Sink Particle number", iptmass," - mass", gaspartmass, " Msol"
          call writeinfo(message, TRIVIAL)
       else

          idead = idead + 1

       endif

    enddo

    write(message,*) "Read ",icount, " lines"
    call writeinfo(message, TRIVIAL)

    write(message,*)  iptmass,"are sink particles and ",igas,"are gas particles and", idead, "are dead"
    call writeinfo(message, TRIVIAL)

    write(message,*) "Total Mass in all particles, ", masscounter, "Msol"
    call writeinfo(message, TRIVIAL)

  end subroutine new_read_sph_data

  !
  !
  ! Read in the data from a file, allocate the array memory, and store the
  ! the number of gas and stars in this object.
  !
  subroutine read_stellar_disc_data(this, filename)
    implicit none
    type(sph_data), intent(inout) :: this
    character(LEN=*), intent(in)  :: filename
    ! 
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
    integer :: nstar
    integer :: i, dum_i
    character(LEN=1) :: dum_a
    real(double), parameter :: M_sun = 1.989e33 ! [grams]

    open(unit=LUIN, file=TRIM(filename), status='old')
    

    nstar = this%nptmass

    ! reading in the header
    do i = 1, 6
       READ(LUIN, *) dum_a
    end do


    ! allocating arrays
    ! --- for stellar discs
    ALLOCATE(this%discrad(nstar))
    ALLOCATE(this%discmass(nstar))
    ALLOCATE(this%spinx(nstar))
    ALLOCATE(this%spiny(nstar))
    ALLOCATE(this%spinz(nstar))


    write(*,*) ' '
    write(*,*) 'Reading Matthew''s stellar disc data....'
    write(*,*) ' '

    
    do i=1, nstar 
       read(luin, *) dum_i, this%discrad(i), this%discmass(i), &
            this%spinx(i), this%spiny(i), this%spinz(i)
       ! convert the mass into grams
       this%discmass(i) = this%discmass(i)*M_sun  ![g]
    end do


    this%have_stellar_disc = .true.    


  end subroutine read_stellar_disc_data
  

    

  !
  !
  ! accessors
  !
  
  ! returns program units of distance in cm 
  function get_udist(this) RESULT(out)
    implicit none
    real(double) :: out
    type(sph_data), intent(in) :: this
    out = this%udist
  end function get_udist

  ! returns program units of mass in g
  function get_umass(this) RESULT(out)
    implicit none
    real(double) :: out
    type(sph_data), intent(in) :: this
    out = this%umass
  end function get_umass

  ! returns program units of time in s
  function get_utime(this) RESULT(out)
    implicit none
    real(double) :: out
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
    real(double) :: out
    type(sph_data), intent(in) :: this
    out = this%time
  end function get_time
    

  !
  ! Returns the position of the i_th gas particle. 
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_position_gas_particle(this, i, x, y, z)
    implicit none    
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i
    real(double), intent(out) :: x, y, z
    
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
    real(double), intent(in) :: value
    
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
    real(double) :: out 
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = this%rhon(i)
    
  end function get_rhon

  function get_mass(this, i) RESULT(out)
    implicit none
    real(double) :: out 
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = this%gasmass(i)
    
  end function get_mass

  ! Returns the temperature of gas particle at the postion of
  ! i-th particle.
  
  function get_temp(this, i) RESULT(out)
    implicit none
    real(double) :: out 
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = this%temperature(i)
    
  end function get_temp

  function get_vel(this, i) RESULT(out)
    implicit none
    type(VECTOR) :: out 
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = VECTOR(this%vxn(i),this%vyn(i),this%vzn(i))
    
  end function get_vel

  ! Assigns the density of gas particle at the postion of
  ! i-th particle.
  
  subroutine put_rhon(this, i, value)
    implicit none
    type(sph_data), intent(inout) :: this
    integer, intent(in) :: i
    real(double), intent(in) :: value

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
    real(double), intent(out) :: x, y, z
    
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
    real(double)  :: out 
    type(sph_data), intent(in) :: this
    integer, intent(in):: i

    out = this%ptmass(i)  !in program units [umass]. See above.

  end function get_pt_mass
  
  
  

  !
  !  destructor
  ! 
  !  Deallocates the array memories

  subroutine kill_sph_data(this)
    implicit none
    type(sph_data), intent(inout) :: this
    
    DEALLOCATE(this%xn, this%yn, this%zn)
    DEALLOCATE(this%rhon, this%temperature)
    DEALLOCATE(this%x, this%y, this%z)
    DEALLOCATE(this%ptmass)

    NULLIFY(this%xn, this%yn, this%zn)
    NULLIFY(this%rhon, this%temperature)
    NULLIFY(this%x, this%y, this%z)
    NULLIFY(this%ptmass)
    
    this%inUse = .false.

  end subroutine kill_sph_data
    




  !
  ! find the maximum distance between the pt masses
  !

  function max_distance(this) RESULT(out)
    implicit none
    real(double) :: out
    type(sph_data), intent(in) :: this
    !
    real(double) :: d_max, d, x, y, z
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
    real(double) :: out
    type(sph_data), intent(in) :: this

    out = MINVAL(this%rhon)
    
  end function get_rhon_min

  
  !
  ! retuns the maximum value of rhon
  !
  function get_rhon_max(this) RESULT(out)
    implicit none
    real(double) :: out
    type(sph_data), intent(in) :: this

    out = MAXVAL(this%rhon)
    
  end function get_rhon_max



  !
  ! retuns the spin direction of i_th star
  !
  subroutine get_spins(this, i, sx, sy, sz) 
    implicit none
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i 
    real(double), intent(out) :: sx, sy, sz 
    
    sx = this%spinx(i)
    sy = this%spiny(i)
    sz = this%spinz(i)
    
  end subroutine get_spins



  !
  ! Prints basic infomation
  !
  ! if filename is '*' then it prints on screen.
  subroutine info(this, filename)
    implicit none
    type(sph_data), intent(in) :: this
    character(LEN=*), intent(in) :: filename
    integer :: UN
    real(double) :: tmp
    
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
    write(UN,*)     'Units of length            : ', get_udist(this), ' [cm]'
    write(UN,*)     'Units of mass              : ', get_umass(this), ' [g]'
    write(UN,*)     'Units of time              : ', get_utime(this),  ' [s]' 
    write(UN,'(a)') ' '    
    write(UN,*)     '# of stars                 : ',  get_nptmass(this)
    write(UN,*)     '# of gas particles (total) : ',  get_npart(this)   
    write(UN,*)     'Time of data dump          : ',  tmp, ' [Myr]'    
    write(UN,'(a)') '#######################################################'
    write(UN,'(a)') ' '
    write(Un,'(a)') ' '
    
    if (filename(1:1) /= '*')  close(UN)

  end subroutine info



  !
  ! Returns the parameters for stellar disc of i-th star
  ! in the data.
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_stellar_disc_parameters(this, i, discrad, discmass, &
       spinx, spiny, spinz)
    implicit none    
    type(sph_data), intent(in) :: this
    integer, intent(in) :: i
    real(double), intent(out) :: discrad    ! in [10^10cm] 
    real(double), intent(out) :: discmass   ! in [grams] 
    real(double), intent(out) :: spinx, spiny, spinz
    logical, save :: first_time = .true.
    
    ! quick check
    if (first_time) then
       if (this%have_stellar_disc) then
          first_time = .false.
       else
          write(*,*) "Error:: You did not read in the stellar disc data! &
               &[get_stellar_disc_parameters]."
          stop
       end if
    end if
    
    discrad = this%discrad(i)    ! in [10^10cm] 
    discmass = this%discmass(i)  ! in [grams] 
    spinx = this%spinx(i)
    spiny = this%spiny(i)
    spinz = this%spinz(i)
        
  end subroutine get_stellar_disc_parameters


  !
  !
  !
  function stellar_disc_exists(this) RESULT(out)
    implicit none
    logical :: out
    type(sph_data), intent(in) :: this
    out = this%have_stellar_disc
  end function stellar_disc_exists




  !
  ! Compute the inclination angles  angle between the spin ]
  ! axis and the observer. The results will be written in 
  ! a file. 
  !
  subroutine find_inclinations(this, obs_x, obs_y, obs_z, outfilename)
    implicit none
    type(sph_data), intent(in) :: this
    real(double), intent(in) :: obs_x,  obs_y, obs_z  ! directions cosines of of observer.
    character(LEN=*), intent(in) :: outfilename 
    real(double) :: r1, r2, inc, dp, pi
    real(double) :: sx, sy, sz ! spins
    integer :: i

    open(unit=43, file=TRIM(ADJUSTL(outfilename)), status="replace")

    pi = 2.0d0*acos(0.0d0)

    write(43, '(a)') "#   star ID    ------   inclination [deg]"

    ! just in case the vectors are not normalized.
    r1 = SQRT(obs_x*obs_x + obs_y*obs_y + obs_z*obs_z)
    do i = 1, this%nptmass       
       call get_spins(this, i , sx, sy, sz)
       ! just in case the vectors are not normalized.
       r2 = SQRT(sx*sx+sy*sy+sz*sz)
    
       ! inclinations
       dp = obs_x*sx + obs_y*sy + obs_z*sz  ! dot product

       inc = ACOS(dp/r1/r2)*(180.d0/Pi)   ! degrees
       if (inc>180.d0) inc = inc-180.d0
       if (inc>90.d0) inc = 90.0 - (inc-90.d0)
       write(43, '(1x, i5, 2x, f6.1)') i, inc
       
    end do
    

    close(43)

  end subroutine find_inclinations
 

  !
  ! Interface function to get the status of SPH data object
  !
  pure function isAlive(this) RESULT(out)
    implicit none
    logical :: out
    type(sph_data), intent(in) :: this
    out = this%inUse
  end function isAlive
   
end module sph_data_class
    
