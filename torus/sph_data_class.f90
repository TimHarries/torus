module sph_data_class

  use kind_mod
  use vector_mod
  use messages_mod
  use utils_mod
  use timing
  use gridtype_mod
  use math_mod2
  use constants_mod, only : OneOnFourPi

!  use parallel_mod, only : torus_abort
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
       ClusterParameter, &
       isAlive

  ! At a given time (time)
  type sph_data
!     private  ! Believe me. It's better to be private!    
     logical      :: inUse=.false.          ! Flag to indicate if this object is in use.
     real(double) :: udist, umass, utime    ! Units of distance, mass, time in cgs
     real(double) :: codeVelocitytoTORUS    ! Conversion from SPH code velocity units to Torus units
     !                                          ! (umass is M_sol, udist=0.1 pc)
     integer          :: npart                  ! Total number of gas particles (field+disc)
     real(double) :: time                   ! Time of sph data dump (in units of utime)
     integer          :: nptmass                ! Number of stars/brown dwarfs
     real(double), pointer, dimension(:) :: gasmass            ! Mass of each gas particle ! DAR changed to allow variable mass
     real(double)                        :: totalgasmass       ! Total gas mass summed over all SPH particles. 
     ! Positions of gas particles
     real(double), pointer, dimension(:) :: xn,yn,zn
     real(double), pointer, dimension(:) :: vxn,vyn,vzn
     real(double), pointer, dimension(:) :: hn                 ! Smoothing length
     ! Density of the gas particles
     real(double), pointer, dimension(:) :: rhon
     ! Temperature of the gas particles
     real(double), pointer, dimension(:) :: temperature
     ! Smoothing lengths of the gas/sink particles
     real(double), pointer, dimension(:) :: h
     real(double), pointer, dimension(:) :: hpt
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
          
     ! Does the temperature of the SPH particles get used?
     logical :: useSphTem 

  end type sph_data

  real(double), allocatable :: PositionArray(:,:), OneOverHsquared(:), oneoverHcubed(:), &
                               RhoArray(:), TemArray(:), VelocityArray(:,:), MassArray(:), Harray(:)
  type(sph_data), save :: sphdata
  integer, save :: npart
  real(double), allocatable :: partArray(:), tempPosArray(:,:), q2array(:)
  integer, allocatable  :: indexArray(:)

  private:: &
       kill_sph_data

  !
  !
  interface kill
     module procedure kill_sph_data
  end interface
  
  
contains  


  ! 
  ! Initializes an object with parameters (if possible).
  ! 
  subroutine init_sph_data(udist, umass, utime,  time, nptmass)
    implicit none

    real(double), intent(in)  :: udist, umass, utime    ! Units of distance, mass, time in cgs
    !                                                       ! (umass is M_sol, udist=0.1 pc)
!    integer, intent(in)           :: npart                  ! Number of gas particles (field+disc)
    real(double), intent(in)  :: time                   ! Time of sph data dump (in units of utime)
    integer, intent(in)           :: nptmass                ! Number of stars/brown dwarfs

    ! Indicate that this object is in use
    sphdata%inUse = .true.

! Do not use temperature from SPH particles to initialise temperature grid
    sphData%useSphTem = .false. 

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
    ALLOCATE(sphdata%hn(npart))

    ! -- for star positions
    ALLOCATE(sphdata%x(nptmass))
    ALLOCATE(sphdata%y(nptmass))
    ALLOCATE(sphdata%z(nptmass))

    ALLOCATE(sphdata%vx(nptmass))
    ALLOCATE(sphdata%vy(nptmass))
    ALLOCATE(sphdata%vz(nptmass))
    ALLOCATE(sphdata%hpt(nptmass))
    
    ! -- for mass of stars
    ALLOCATE(sphdata%ptmass(nptmass))

    !
    sphdata%have_stellar_disc = .false. 

  end subroutine init_sph_data

  ! 
  ! Initializes an object with parameters when torus is called as a subroutine from sphNG.
  ! 
  subroutine init_sph_data2(udist, umass, utime, b_num_gas, time, nptmass, &
        b_npart, b_idim, b_iphase, b_xyzmh, b_rho, b_temp, b_totalgasmass)
    implicit none

! Arguments --------------------------------------------------------------------
    real(double), intent(in)  :: udist, umass, utime    ! Units of distance, mass, time in cgs
    !                                                   ! (umass is M_sol, udist=0.1 pc)
    integer, intent(in)           :: b_num_gas          ! Number of gas particles (field+disc)
    real(double), intent(in)  :: time                   ! Time of sph data dump (in units of utime)
    integer, intent(in)           :: nptmass            ! Number of stars/brown dwarfs

    integer, intent(in)   :: b_npart, b_idim
    integer*1, intent(in) :: b_iphase(b_idim)
    real*8, intent(in)    :: b_xyzmh(5,b_idim)
    real*4, intent(in)    :: b_rho(b_idim)
    real*8, intent(in)    :: b_temp(b_idim)
    real*8, intent(in)    :: b_totalgasmass
! Local variables --------------------------------------------------------------
    integer :: iii, iiipart, iiigas

! Begin executable statements --------------------------------------------------

    npart = b_num_gas

    ! Indicate that this object is in use
    sphData%inUse = .true.

! Use temperature from SPH particles to initialise temperature grid
    sphData%useSphTem = .true. 

    ! save these values in this object
    sphData%udist = udist
    sphData%umass = umass
    sphData%utime = utime
    sphData%npart = npart
    sphData%time = time
    sphData%nptmass = nptmass
    sphData%totalgasmass = b_totalgasmass


    ! allocate arrays
    ! -- for gas particles
    ALLOCATE(sphData%xn(npart))
    ALLOCATE(sphData%yn(npart))
    ALLOCATE(sphData%zn(npart))
    ALLOCATE(sphData%hn(npart))
    ALLOCATE(sphData%rhon(npart))
    ALLOCATE(sphData%temperature(npart))
    ALLOCATE(sphData%gasmass(npart))


    ! -- for star positions
    ALLOCATE(sphData%x(nptmass))
    ALLOCATE(sphData%y(nptmass))
    ALLOCATE(sphData%z(nptmass))
    
    ! -- for mass of stars
    ALLOCATE(sphData%ptmass(nptmass))

    !
    sphData%have_stellar_disc = .false. 
    
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
          sphData%rhon(iiigas)        = b_rho(iii)
          sphData%temperature(iiigas) = b_temp(iii)
          sphData%xn(iiigas)          = b_xyzmh(1,iii)
          sphData%yn(iiigas)          = b_xyzmh(2,iii)
          sphData%zn(iiigas)          = b_xyzmh(3,iii)
          sphData%gasmass(iiigas)     = b_xyzmh(4,iii)
          sphData%hn(iiigas)          = b_xyzmh(5,iii)

       elseif (b_iphase(iii) > 0) then
          iiipart=iiipart+1
          if (iiipart > nptmass) then
             write (*,*) 'Error: iiipart>nptmass',iiipart,nptmass
             cycle
          endif
          sphData%x(iiipart)      = b_xyzmh(1,iii)
          sphData%y(iiipart)      = b_xyzmh(2,iii)
          sphData%z(iiipart)      = b_xyzmh(3,iii)
          sphData%ptmass(iiipart) = b_xyzmh(4,iii)
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
    integer :: nptmass, n1, n2
    real(double), allocatable :: dummy(:)     


    open(unit=LUIN, file=TRIM(filename), form='unformatted')
  

    ! reading in the first line
    READ(LUIN) udist, umass, utime, npart, n1, n2, time, nptmass, gaspartmass

    ! initilaizing the sph_data object (allocating arrays, saving parameters and so on....)
    call init_sph_data(udist, umass, utime, time, nptmass)


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

  subroutine new_read_sph_data(filename)
    implicit none

    character(LEN=*), intent(in)  :: filename
    !   
    integer, parameter  :: LUIN = 10 ! logical unit # of the data file
!    real(double) :: udist, umass, utime,  time,  gaspartmass, discpartmass
!    integer*4 :: npart,  nsph, nptmass
    real(double) :: udist, umass, utime,  time
    real(double) :: xn, yn, zn, vx, vy, vz, gaspartmass, rhon, masscounter, u, h
!    real(double), allocatable :: xarray(:), harray(:)
!    INTEGER, allocatable :: ind(:)
    integer :: itype, ipart, icount, iptmass, igas, idead
    integer :: nptmass, n1, n2, nlines
    real(double) junk
    character(LEN=1)  :: junkchar
    character(LEN=150) :: message
!    real(double) :: hcrit

    open(unit=LUIN, file=TRIM(filename), form="formatted")

    read(LUIN,*) 
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*) junkchar, time, utime
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*) junkchar, npart, n1, nptmass, n2
    read(LUIN,*)
!    read(LUIN,*) junkchar, udist, junk, junk, umass, junk, urho, uvel, junk, junk
    read(LUIN,*) junkchar, udist, junk, junk, umass
    read(LUIN,*)
    read(LUIN,*)
    read(LUIN,*)

    write(message,*) "Allocating ", npart, " gas particles and ", nptmass, " sink particles"
    call writeinfo(message, TRIVIAL)
    call init_sph_data(udist, umass, utime, time, nptmass)
    ! velocity unit is derived from distance and time unit (converted to seconds from years)
    sphdata%codeVelocitytoTORUS = (udist / (utime * 31536000.)) / cspeed 

    nlines = npart + n2 + nptmass + n1 ! npart now equal to no. lines - 12 = sum of particles dead or alive

    write(message,*) "Reading SPH data from ASCII...."
    call writeinfo(message, TRIVIAL)

    iptmass = 0
    icount = 12
    igas = 0
    idead = 0
    masscounter = 0.d0

    do ipart=1, nlines

       read(LUIN,*) xn, yn, zn, gaspartmass, h, rhon, vx, vy, vz, u, junk, junk, junk, itype
       icount = icount + 1

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
          sphdata%temperature(igas) = 1.9725e-8 * u

          sphdata%hn(igas) = h

          sphdata%totalgasmass = sphdata%totalgasmass + gaspartmass
          
       Elseif(itype .eq. 3 ) then

          iptmass = iptmass + 1

          sphdata%x(iptmass) = xn
          sphdata%y(iptmass) = yn
          sphdata%z(iptmass) = zn

          sphdata%vx(iptmass) = vx
          sphdata%vy(iptmass) = vy
          sphdata%vz(iptmass) = vz

          sphdata%ptmass(iptmass) = gaspartmass 
          sphdata%hpt(iptmass) = h

          sphdata%totalgasmass = sphdata%totalgasmass + gaspartmass

          write(message,*) "Sink Particle number", iptmass," - mass", gaspartmass, " Msol - Index", iptmass + igas
          call writeinfo(message, TRIVIAL)
          write(98,*) iptmass, xn*udist*1e-10, yn*udist*1e-10, zn*udist*1e-10

       else

          idead = idead + 1

       endif

    enddo

    write(message,*) "Read ",icount, " lines"
    call writeinfo(message, TRIVIAL)

    write(message,*)  iptmass," are sink particles and ",igas," are gas particles and ", idead, " are dead"
    call writeinfo(message, TRIVIAL)

    write(message,*) "Total Mass in all particles, ", sphdata%totalgasmass, " Msol"
    call writeinfo(message, TRIVIAL)
    
    close(LUIN)
   
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
!  function get_udist(this) RESULT(out)
  function get_udist() RESULT(out)
    implicit none
    real(double) :: out
!    type(sph_data), intent(in) :: this
    out = sphdata%udist
  end function get_udist

  ! returns program units of mass in g
!  function get_umass(this) RESULT(out)
  function get_umass() RESULT(out)
    implicit none
    real(double) :: out
!    type(sph_data), intent(in) :: this
    out = sphdata%umass
  end function get_umass

  ! returns program units of time in s
!  function get_utime(this) RESULT(out)
  function get_utime() RESULT(out)
    implicit none
    real(double) :: out
!    type(sph_data), intent(in) :: this
    out = sphdata%utime
  end function get_utime
    

  ! returns the number of gas particles
  function get_npart() RESULT(out)
    implicit none
    integer ::out 
!    type(sph_data), intent(in) :: this
    out = sphdata%npart
  end function get_npart
    

  ! returns the number of point masses
!  function get_nptmass(this) RESULT(out)
  function get_nptmass() RESULT(out)
    implicit none
    integer :: out
!    type(sph_data), intent(in) :: this
    out = sphdata%nptmass
  end function get_nptmass
    

  ! returns the time of dump time in the unit of [utime]
 ! function get_time(this) RESULT(out)
 function get_time() RESULT(out)
    implicit none
    real(double) :: out
!    type(sph_data), intent(in) :: this
    out = sphdata%time
  end function get_time
    

  !
  ! Returns the position of the i_th gas particle. 
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_position_gas_particle(i, x, y, z)
    implicit none    
!    type(sph_data), intent(in) :: this
    integer, intent(in) :: i
    real(double), intent(out) :: x, y, z
    
    x = sphdata%xn(i) 
    y = sphdata%yn(i)
    z = sphdata%zn(i)
    
  end subroutine get_position_gas_particle

  !
  ! saves the position of the i_th gas particle in this object. 
  !
  ! "name" must be one of the following: "x", "y", "z"
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine put_position_gas_particle(i, name, value)
    implicit none    
!    type(sph_data), intent(inout) :: this
    integer, intent(in) :: i
    character(LEN=*), intent(in) :: name
    real(double), intent(in) :: value
    
    select case(name)
    case ("x", "X")
       sphdata%xn(i) = value
    case ("y", "Y") 
       sphdata%yn(i) = value
    case ("z", "Z") 
       sphdata%zn(i) = value
    case default
       write(*,*) "Error: Unknown name passed to sph_data_class::put_position_gas_particle."
       stop
    end select
    
  end subroutine put_position_gas_particle


  ! Returns the density of gas particle at the postion of
  ! i-th particle.
  
  function get_rhon(i) RESULT(out)
    implicit none
    real(double) :: out 
!    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = sphdata%rhon(i)
    
  end function get_rhon

  function get_mass(i) RESULT(out)
    implicit none
    real(double) :: out 
!    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = sphdata%gasmass(i)
    
  end function get_mass

  ! Returns the temperature of gas particle at the postion of
  ! i-th particle.
  
  function get_temp(i) RESULT(out)
    implicit none
    real(double) :: out 
!    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = sphdata%temperature(i)
    
  end function get_temp

  function get_vel(i) RESULT(out)
    implicit none
    type(VECTOR) :: out 
!    type(sph_data), intent(in) :: this
    integer, intent(in) :: i

    out  = VECTOR(sphdata%vxn(i),sphdata%vyn(i),sphdata%vzn(i))
    
  end function get_vel

  ! Assigns the density of gas particle at the postion of
  ! i-th particle.
  
  subroutine put_rhon(i, value)
    implicit none
!    type(sph_data), intent(inout) :: this
    integer, intent(in) :: i
    real(double), intent(in) :: value

    sphdata%rhon(i) = value
    
  end subroutine put_rhon
     
  !
  ! Rerurns the postions of the i-th point mass.
  !
  ! x, y, z are in [udist] ... See the type definition section.
  subroutine get_position_pt_mass(i, x, y, z)
    implicit none    
!    type(sph_data), intent(in) :: this
    integer, intent(in) :: i
    real(double), intent(out) :: x, y, z
    
    x = sphdata%x(i) 
    y = sphdata%y(i)
    z = sphdata%z(i)
    
  end subroutine get_position_pt_mass
 

  !
  ! Returns the mass of the point mass (star) in
  ! the i-th position of data.
  !
  function get_pt_mass(i) RESULT(out)
    implicit none
    real(double)  :: out 
!    type(sph_data), intent(in) :: this
    integer, intent(in):: i

    out = sphdata%ptmass(i)  !in program units [umass]. See above.

  end function get_pt_mass
  
  
  

  !
  !  destructor
  ! 
  !  Deallocates the array memories

  subroutine kill_sph_data()
    implicit none
!    type(sph_data), intent(inout) :: this
    
    DEALLOCATE(sphdata%xn, sphdata%yn, sphdata%zn)
    DEALLOCATE(sphdata%rhon, sphdata%temperature)
    DEALLOCATE(sphdata%x, sphdata%y, sphdata%z)
    DEALLOCATE(sphdata%ptmass)

    NULLIFY(sphdata%xn, sphdata%yn, sphdata%zn)
    NULLIFY(sphdata%rhon, sphdata%temperature)
    NULLIFY(sphdata%x, sphdata%y, sphdata%z)
    NULLIFY(sphdata%ptmass)

    if ( ASSOCIATED(sphdata%gasmass) ) then 
       DEALLOCATE(sphdata%gasmass)
       NULLIFY(sphdata%gasmass)
    end if

    if ( ASSOCIATED(sphdata%hn) ) then 
       DEALLOCATE(sphdata%hn)
       NULLIFY(sphdata%hn)
    end if

    sphdata%inUse = .false.

  end subroutine kill_sph_data
    




  !
  ! find the maximum distance between the pt masses
  !

  function max_distance() RESULT(out)
    implicit none
    real(double) :: out
!    type(sph_data), intent(in) :: this
    !
    real(double) :: d_max, d, x, y, z
    integer :: i, j, n
    
    d_max= -1.0
    n=get_nptmass() ! function in this moudle
    
    do i = 1, n-2
       do j = i+1, n
          call get_position_pt_mass(j, x, y, z)
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
  function get_rhon_min() RESULT(out)
    implicit none
    real(double) :: out
!    type(sph_data), intent(in) :: this

    out = MINVAL(sphdata%rhon)
    
  end function get_rhon_min

  
  !
  ! retuns the maximum value of rhon
  !
  function get_rhon_max() RESULT(out)
    implicit none
    real(double) :: out
!    type(sph_data), intent(in) :: this

    out = MAXVAL(sphdata%rhon)
    
  end function get_rhon_max



  !
  ! retuns the spin direction of i_th star
  !
  subroutine get_spins(i, sx, sy, sz) 
    implicit none
!    type(sph_data), intent(in) :: this
    integer, intent(in) :: i 
    real(double), intent(out) :: sx, sy, sz 
    
    sx = sphdata%spinx(i)
    sy = sphdata%spiny(i)
    sz = sphdata%spinz(i)
    
  end subroutine get_spins



  !
  ! Prints basic infomation
  !
  ! if filename is '*' then it prints on screen.
  subroutine info(filename)
    implicit none
!    type(sph_data), intent(in) :: this
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
    tmp = get_time()*get_utime()/(60.0d0*60.0d0*24.0d0*365.0d0*1.0d6)
    
    write(UN,'(a)') ' '
    write(UN,'(a)') '######################################################'
    write(UN,'(a)') 'SPH data info :'
    write(UN,'(a)') ' '    
    write(UN,*)     'Units of length            : ', get_udist(), ' [cm]'
    write(UN,*)     'Units of mass              : ', get_umass(), ' [g]'
    write(UN,*)     'Units of time              : ', get_utime(),  ' [s]' 
    write(UN,'(a)') ' '    
    write(UN,*)     '# of stars                 : ',  get_nptmass()
    write(UN,*)     '# of gas particles (total) : ',  get_npart()   
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
  subroutine get_stellar_disc_parameters(i, discrad, discmass, &
       spinx, spiny, spinz)
    implicit none    
!    type(sph_data), intent(in) :: this
    integer, intent(in) :: i
    real(double), intent(out) :: discrad    ! in [10^10cm] 
    real(double), intent(out) :: discmass   ! in [grams] 
    real(double), intent(out) :: spinx, spiny, spinz
    logical, save :: first_time = .true.
    
    ! quick check
    if (first_time) then
       if (sphdata%have_stellar_disc) then
          first_time = .false.
       else
          write(*,*) "Error:: You did not read in the stellar disc data! &
               &[get_stellar_disc_parameters]."
          stop
       end if
    end if
    
    discrad = sphdata%discrad(i)    ! in [10^10cm] 
    discmass = sphdata%discmass(i)  ! in [grams] 
    spinx = sphdata%spinx(i)
    spiny = sphdata%spiny(i)
    spinz = sphdata%spinz(i)
        
  end subroutine get_stellar_disc_parameters


  !
  !
  !
  function stellar_disc_exists() RESULT(out)
    implicit none
    logical :: out
!    type(sph_data), intent(in) :: this
    out = sphdata%have_stellar_disc
  end function stellar_disc_exists




  !
  ! Compute the inclination angles  angle between the spin ]
  ! axis and the observer. The results will be written in 
  ! a file. 
  !
  subroutine find_inclinations(obs_x, obs_y, obs_z, outfilename)
    implicit none
!    type(sph_data), intent(in) :: this
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
    do i = 1, sphdata%nptmass       
       call get_spins(i , sx, sy, sz)
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
  pure function isAlive() RESULT(out)
    implicit none
    logical :: out
!    type(sph_data), intent(in) :: this
    out = sphdata%inUse
  end function isAlive

  subroutine sortbyx(xarray,ind)
    real(double) :: xarray(:)
    integer :: ind(:)
    integer :: i
      
    do i=1,npart
       ind(i) = i
    enddo
    
    call sortdouble2index(xarray,ind)

  end subroutine sortbyx

  subroutine FindCriticalValue(array,critval,percentile,output)

    real(double) :: array(:), critval
    real(double), intent(in) :: percentile
    integer :: i
    logical, optional :: output 

    character(len=100) :: message

!    npart = size(array)

    if(output .and. writeoutput) then
    
       call writeinfo("Start sort", TRIVIAL)
       call tune(6, "Sort")  ! start a stopwatch
       call dquicksort(array)
       call tune(6, "Sort")  ! stop a stopwatch
       call writeinfo("Start sort", TRIVIAL)

       open(unit=50, file='sortedarray.dat', status="replace")
       open(unit=49, file='percentiles.dat', status="replace")
    
       do i=1,npart
          write(50,*) array(i)
       enddo

       do i=100,1,-1
          write(49,*) i,"%",array(nint((npart*i/100.)))
       enddo
    
       write(message,*) "Maximum Value",array(npart)
       call writeinfo(message, FORINFO)
       write(message,*) "Minimum Value",array(1)
       call writeinfo(message, FORINFO)
    else
       call dquicksort(array)
    endif

    critval = array(nint(percentile*npart))

  end subroutine FindCriticalValue

  TYPE(vector)  function Clusterparameter(point, grid, theparam, isdone, shouldreuse, d, RhoMin, RhoMax)
    USE input_variables, only: hcritPercentile, hmaxPercentile

    type(vector), intent(in) :: point
    type(GRIDTYPE), intent(in), optional :: grid
    type(vector) :: posvec

    integer :: i

    integer, optional :: theparam
    integer :: param

    logical, save :: firsttime = .true.
    logical, optional :: isdone
    logical :: done
    
    integer, allocatable, save :: ind(:)
    real(double), save :: hcrit, hmax, rcrit, OneOverHcrit, OneOverhMax, rmax
    real(double) :: codeVelocitytoTORUS, codeLengthtoTORUS, codeDensitytoTORUS, udist, umass, utime
    real(double) :: r
    real(double), optional :: d
    real(double) :: fac
    real(double), save :: sumWeight
    real(double) :: paramValue(4)
    integer,save :: nparticles
    
    character(len=100) :: message

    logical, save :: notfound
    logical, optional :: shouldreuse
    logical :: reuse

    real(double), optional:: RhoMin, RhoMax

    if(present(rhomin)) then
       rhomin = 1d30
       rhomax = -1d30
    endif

    if(present(shouldreuse)) then
       reuse = shouldreuse
    else
       reuse = .false.
    endif

    if(present(theparam)) then
       param = theparam
    else
       param = 1
    endif

    if(present(isdone)) then
       done = isdone
    else
       done = .false.
    endif

    if(firsttime) then
!       call new_read_sph_data(tempsphdata, "newsph.dat.ascii") ! read in sphdata

       udist = get_udist()
       utime = get_utime()
       umass = get_umass()
       codeLengthtoTORUS = udist * 1d-10
       codeVelocitytoTORUS = sphdata%codeVelocitytoTORUS
       codeDensitytoTORUS = umass / ((udist) ** 3)

!       npart = get_npart() ! total gas particles

       allocate(PositionArray(npart,3)) ! allocate memory
       allocate(q2Array(npart))
       allocate(RhoArray(npart))
       allocate(TemArray(npart))
       allocate(MassArray(npart))
       allocate(Harray(npart))
       allocate(ind(npart))
       allocate(OneOverHsquared(npart))
       allocate(OneOverHcubed(npart))
       allocate(tempPosArray(npart,3))

       PositionArray = 0.d0; MassArray = 0.d0; hArray = 0.d0; ind = 0; tempPosArray = 0.d0; q2array = 0.d0

       PositionArray(:,1) = sphdata%xn(:) * codeLengthtoTORUS! fill with x's to be sorted

       call sortbyx(PositionArray(:,1),ind(:)) ! sort the x's and recall their indices

       PositionArray(:,2) = sphdata%yn(ind(:)) * codeLengthtoTORUS ! y's go with their x's
       PositionArray(:,3) = sphdata%zn(ind(:)) * codeLengthtoTORUS ! z's go with their x's

       write(message, *) "Max abs(x)", max(abs(PositionArray(1,1)),abs(PositionArray(npart,1)))
       call writeinfo(message, TRIVIAL)
       write(message, *) "Max abs(y)", maxval(PositionArray(:,2))
       call writeinfo(message, TRIVIAL)
       write(message, *) "Max abs(z)", maxval(PositionArray(:,3))
       call writeinfo(message, TRIVIAL)

! Decide if we need to set velocities for this configuration
       if ( associated(sphData%vxn) ) then 

          allocate(VelocityArray(3,npart))
          VelocityArray(:,:) = 0.d0
          VelocityArray(1,:) = sphdata%vxn(ind(:)) * codeVelocitytoTORUS ! velocities

          if ( associated(sphData%vyn) ) then 
             VelocityArray(2,:) = sphdata%vyn(ind(:)) * codeVelocitytoTORUS 
          else
             STOP
!             call torus_abort("Error in function Clusterparameter: vxn is associated but vyn is not")
          end if

          if ( associated(sphData%vzn) ) then 
             VelocityArray(3,:) = sphdata%vzn(ind(:)) * codeVelocitytoTORUS
          else
             STOP
!             call torus_abort("Error in function Clusterparameter: vxn is associated but vzn is not")
          end if

       end if

       RhoArray(:) = sphdata%rhon(ind(:)) * codeDensitytoTORUS
       TemArray(:) = sphdata%temperature(ind(:))

       MassArray(:) = sphdata%gasmass(ind(:)) * umass! in case of unequal mass
       Harray(:) = sphdata%hn(ind(:)) ! fill h array

       call FindCriticalValue(harray, hcrit, real(hcritPercentile,kind=db), output = .false.) ! find hcrit as percentile of total h
       call FindCriticalValue(harray, hmax,  real(hmaxPercentile, kind=db), output = .false.) ! find hmax as percentile of total h
       Harray(:) = sphdata%hn(ind(:)) ! fill h array

       write(message, *) "Critical smoothing Length in code units", hcrit
       call writeinfo(message, TRIVIAL)
       write(message, *) "Maximum smoothing Length in code units", hmax
       call writeinfo(message, TRIVIAL)

       hcrit = hcrit * codeLengthtoTORUS
       OneOverHcrit = 1.d0 / hcrit

       hmax = hmax * codeLengthtoTORUS
       OneOverHmax = 1.d0 / hmax

       Harray(:) = Harray(:) * codeLengthtoTORUS  ! fill h array
       OneOverHsquared(:) = 1.d0 / (Harray(:)**2)
       OneOverHcubed(:) = 1.d0 / (Harray(:)**3)

       write(message,*) "Critical smoothing Length in 10^10cm", hcrit
       call writeinfo(message, TRIVIAL)
       write(message,*) "Maximum smoothing Length in 10^10cm", hmax
       call writeinfo(message, TRIVIAL)

       rcrit = 2.d0 * hcrit ! edge of smoothing sphere
       rmax = 2.d0 * hmax ! edge of smoothing sphere

!       call kill() ! don't need harray anymore
       allocate(partarray(npart), indexarray(npart))

       firsttime = .false.
    endif

    if(done) then

       deallocate(PositionArray, MassArray, harray, RhoArray, Temarray, ind, tempPosArray, q2Array)

       deallocate(OneOverHsquared, OneOverHcubed)
       deallocate(partarray, indexarray)

       if (allocated(VelocityArray)) deallocate (VelocityArray)
       firsttime = .true.
       return
    endif

    posVec = point
    r = rcrit ! CHECK HERE!!!

    if(present(d)) then
       r = 1.75d0 * d + rcrit ! 1.75 is like sqrt(3)!
!       r = rmax
    endif

    if(reuse) then
       call doweights(posvec, nparticles, sumweight, qpresent = .false.)
    else
       notfound = .false.
       call findNearestParticles(posvec, nparticles, r, expkernel = .true.)
       call doweights(posvec, nparticles, sumweight, qpresent = .true.)
    endif

    if(sumweight .le. 0.d0) notfound = .true.

    paramvalue(:) = 0.d0

    if(.not. notfound) then

       if(sumweight .gt. 0.3d0) then
          fac = 1.d0 / sumWeight
       else
          fac = 1.d0
       endif

       if(param .eq. 1) then
          do i = 1, nparticles
             paramValue(1) = paramValue(1) + partArray(i) * VelocityArray(1, indexArray(i)) ! Vx
             paramValue(2) = paramValue(2) + partArray(i) * VelocityArray(2, indexArray(i)) ! Vy
             paramValue(3) = paramValue(3) + partArray(i) * VelocityArray(3, indexArray(i)) ! Vz
          enddo

          Clusterparameter = VECTOR(paramValue(1) * fac, paramValue(2) * fac, paramValue(3) * fac) ! Velocity 

!          write(69,*) nparticles, sumweight

       elseif(param .eq. 2) then

          do i = 1, nparticles
             paramValue(3) = paramValue(3) + partArray(i) * TemArray(indexArray(i)) ! Temperature
             paramValue(4) = paramValue(4) + partArray(i) * RhoArray(indexArray(i)) ! rho
!             if(present(rhomin)) then ! have to have rhomax with rhomin !!!
!                RhoMin = min(Rhomin, RhoArray(indexArray(i))) ! rhomin
!                RhoMax = max(Rhomax, RhoArray(indexArray(i))) ! rhomax
!             endif
          enddo
          
          Clusterparameter = VECTOR(paramValue(4)*fac, paramValue(3)*fac, 0.d0)  ! density ! stays as vector for moment
       endif
    else
       if(param .eq. 1) then
          if(.not. reuse) then
             call findNearestParticles(posvec, nparticles, rmax, expkernel = .true.) ! redo but with exponential kernel
          endif

          call doweights(posvec, nparticles, sumweight, qpresent = .true.)
          
          if(sumweight .ne. 0.d0) then
             fac = 1.d0 / sumWeight
          else
             fac = 1.d0
          endif

          do i = 1, nparticles
             paramValue(1) = paramValue(1) + partArray(i) * VelocityArray(1, indexArray(i)) ! Vx
             paramValue(2) = paramValue(2) + partArray(i) * VelocityArray(2, indexArray(i)) ! Vy
             paramValue(3) = paramValue(3) + partArray(i) * VelocityArray(3, indexArray(i)) ! Vz
          enddo
!          write(69,*) "0", sumweight
          Clusterparameter = VECTOR(paramValue(1) * fac, paramValue(2) * fac, paramValue(3) * fac) ! Velocity 
!          Clusterparameter = VECTOR(-1.d20,1.d20,-1.d20) ! Velocity

       elseif(param .eq. 2) then
          if(.not. reuse) then
             call findNearestParticles(posvec, nparticles, rmax, expkernel = .true.) ! redo but with exponential kernel
          endif

          call doweights(posvec, nparticles, sumweight, qpresent = .true.)
  
          do i = 1, nparticles
             paramValue(3) = paramValue(3) + partArray(i) * TemArray(indexArray(i)) ! Temperature
             paramValue(4) = paramValue(4) + partArray(i) * RhoArray(indexArray(i)) ! rho
!             if(present(rhomin)) then ! have to have rhomax with rhomin !!!
!                RhoMin = min(Rhomin, RhoArray(indexArray(i))) ! rhomin
!                RhoMax = max(Rhomax, RhoArray(indexArray(i))) ! rhomax
!             endif
          enddo

          paramvalue(4) = max(1d-60, paramvalue(4))
          Clusterparameter = VECTOR(paramValue(4), paramValue(3), 0.d0) ! density ! stays as vector for moment
       endif
    endif

  end function Clusterparameter

  subroutine findnearestparticles(pos, partcount, r, expkernel)
    type(VECTOR) :: pos
    real(double) :: x,y,z
    integer :: i
    integer, save :: nupper, nlower
    integer, save :: closestXindex, testIndex
    integer, intent(out) :: partcount
    real(double) :: weightFac
    real(double) :: r2test, q2test, qtest, r
    real(double) :: ydiff, zdiff, rr
    logical, optional :: expkernel
    logical :: doexpkernel
    
    real(double), parameter :: OneOversqrtPiCubed = 1.d0 / 5.568328000d0
    real(double), parameter :: num = 125.d0 / 216.d0 ! (5/6)^3 ! (1/1.2^3)

    integer :: stepsize, sense
    logical :: test, prevtest, up

    if(present(expkernel)) then
       doexpkernel = expkernel
    else
       doexpkernel = .false.
    endif
  
    x = pos%x
    y = pos%y
    z = pos%z

    call locate_double_f90(PositionArray(:,1), x, closestXindex) ! find the nearest particle to your point
  
    nupper = 1
    nlower = -1

    rr = r**2
    
    closestxIndex = min(max(1, closestXindex), npart)

    if(closestxIndex .lt. npart) then
  
    up = .true.
    nupper = 1
    stepsize = 1
    sense = 1
    prevtest = .true.

    do while (stepsize .ge. 1)
       test = abs(PositionArray(min(npart,closestXindex + nupper),1) - x) .le. r
 
       if(test .and. (nupper .eq. npart - closestXindex)) exit
       if(closestXindex + nupper .ge. npart) nupper = npart - closestXindex

       if(.not. (test .eqv. prevtest)) then
          sense = -sense ! forwards or backwards
          up = .false.
       endif
       
       if(up) then
          stepsize = 2 * stepsize
       else
          stepsize = stepsize / 2
       endif
    
       if(stepsize .lt. 1) then
          if (test) then
             nupper = nupper
          else
             nupper = nupper - 1 ! always has to be lower
          endif
          exit
       endif

       nupper = min(nupper + sense * stepsize, npart - closestXindex)

       prevtest = test
    enddo

    else
       nupper = 0
    endif
! repeat for nlower
    if(closestXindex .gt. 1) then

    up = .true.
    nlower = -1
    stepsize = 1
    sense = -1
    prevtest = .true.

    do while (stepsize .ge. 1)
       test = abs(PositionArray(max(1,closestXindex + nlower),1) - x) .le. r
 
       if(test .and. (nlower .eq. 1 - closestXindex)) exit
       if(closestXindex + nlower .le. 1) nlower = 1 - closestXindex

       if(.not. (test .eqv. prevtest)) then
          sense = -sense ! forwards or backwards
          up = .false.
       endif
       
       if(up) then
          stepsize = 2 * stepsize
       else
          stepsize = stepsize / 2
       endif
    
       if(stepsize .lt. 1) then
          if (test) then
             nlower = nlower
          else
             nlower = nlower + 1 ! always has to be higher
          endif
          exit
       endif

       nlower = max(nlower + sense * stepsize, 1 - closestXindex)

       prevtest = test
    enddo

    else
       nlower = 0
    endif
       
    partcount = 0

    do i = nlower, nupper ! search over all those particles we just found
       
       testIndex = closestXindex + i
       
       zdiff = (PositionArray(testIndex,3) - z) ** 2

       if(zdiff .le. rr) then ! if it's near in y then
          ydiff = (PositionArray(testIndex,2) - y)**2

          if(ydiff .le. rr) then !only if it's near in z as well then work out contribution

             r2test = (PositionArray(testIndex,1) - x)**2 + & !The kernel will do the rest of the work for those outside the sphere
                  ydiff + zdiff

             q2test = r2test * OneOverHsquared(testIndex) ! dimensionless parameter that we're interested in

             if(q2test .lt. 4.d0) then
                partcount = partcount + 1
!                Weightfac = num * SmoothingKernel3d(qtest) ! normalised contribution from this particle
!                Weightfac = num * OneOversqrtPiCubed * exp(-q2test)
                indexArray(partcount) = testIndex
                tempPosArray(partcount,1:3) = PositionArray(partcount,1:3)
                q2array(partcount) = q2test
             else
                if(q2test .lt. 100.d0) then
                   if(doexpkernel) then
                      partcount = partcount + 1
                      indexArray(partcount) = testIndex
                      tempPosArray(partcount,1:3) = PositionArray(partcount,1:3)
                      q2array(partcount) = q2test
                   endif
                endif
             endif
          endif
       endif
    enddo

  end subroutine findnearestparticles

  subroutine doWeights(posvec, partcount, sumweight, qpresent)

    type(VECTOR) :: posvec
    integer :: partcount
    real(double) :: x,y,z
    logical, optional :: qpresent
    logical :: q
    real(double) :: r2test(3500000)
    real(double), parameter :: OneOversqrtPiCubed = 1.d0 / 5.568328000d0
    real(double), parameter :: num = 125.d0 / 216.d0 ! (5/6)^3 ! (1/1.2^3)
    real(double) :: sumweight

    if(present(qpresent)) then
       q = qpresent
    else
       q = .false.
    endif

    if(partcount .gt. 0) then       
    
       x = posvec%x
       y = posvec%y
       z = posvec%z
       
       if(.not. q) then
          r2test(1:partcount) = (tempPosArray(1:partcount,1) - x)**2 + & !The kernel will do the rest of the work for those outside the sphere
               (tempPosArray(1:partcount,2) - y)**2 + (tempPosArray(1:partcount,3) - z)**2
          q2array(1:partcount) = r2test(1:partcount) * OneOverHsquared(indexArray(1:partcount)) ! dimensionless parameter that we're interested in
       endif
       
       partarray(1:partcount) = num * OneOversqrtPiCubed * exp(-q2array(1:partcount))
       sumWeight = sum(partarray(1:partcount))
    else
       sumweight = 0.d0
    endif
    
  end subroutine doWeights

  real(double) pure function SmoothingKernel3d(q) result(Weight)
      
    real(double), intent(in) :: q
    real(double) :: qminus2, qminus1

    qminus2 = q - 2.d0
    qminus1 = q - 1.d0
   
    if(qminus2 .lt. 0.d0) then
       Weight = (-qminus2) ** 3
       if(qminus1 .lt. 0.d0) then
          Weight = Weight + 4.d0 * ((qminus1) ** 3)
       endif
       Weight = OneOnFourPi * Weight
    else
       Weight = 0.d0
    endif
    
  end function SmoothingKernel3d
  
end module sph_data_class
    
