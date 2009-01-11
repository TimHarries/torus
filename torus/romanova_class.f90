module romanova_class
  ! --------------------------------------------------------------------
  !  The data class to handle the data output from Romanova's tilted
  !  dipole magnetosphere MHD simulation.  
  !---------------------------------------------------------------------
  !  Created: 05-april-2006  (R.Kurosawa)
  !---------------------------------------------------------------------

  use kind_mod
  use constants_mod, only: pi, cSpeed_dbl, cSpeed
  use octal_mod
  use zeus_data_class
  use vector_mod

  IMPLICIT NONE
    
  public :: &
       new, &
       kill, &
       get_dble_parameter, &
       get_log_parameter, &
       romanova_data_alive, &
       romanova_density, &
       romanova_velocity, &
       calc_romanova_mass_velocity, &
       calc_romanova_temperature

  private :: &
       init_romanova, &
       get_coordinates_romanova, &
       get_spherical_coordinates, &
       VelocityCorners

  
  
  type romanova
     private
     real(double)    :: Rs          ! radius of central star  [10^10cm]
     real(double)    :: Rmax        ! cutoff radius in [10^10 cm]
     real(double)    :: Mass        ! [g]  Mass of the star
     ! If T, T_flow in all space. If F, the temperature from Romanova's
     ! data file will be used.
     logical         :: isothermal  
     real(double)    :: T_flow       ! [K]  Isothemal temperature of the flow
     ! Reference values for dimensionless units used in Romanova's model
     real(double)    :: r_ref       ! [cm] Reference length value
     real(double)    :: rho_ref     ! [g/cm^3]  Reference density value
     real(double)    :: T_ref       ! [K]       Reference temperature
     real(double)    :: v_ref       ! [cm/s]    Reference speed 
     ! The tile angle of the magnetic axis 
     real(double)    :: tilt        ! [rad.]
     real(double)    :: Period      ! [sec] rotational period
     !
     type(zeus_data) :: romanova_output ! Data structure for romanova's data
     ! flag to indicate if romanova_output memory for is allocated.
     logical         :: data_alive = .false.  
  end type romanova

  !
  !
  !
  interface new
     module procedure init_romanova
  end interface

  interface kill
     module procedure deallocate_romanova_data
  end interface


!  ! This should be used only in this module.
!  ! This should be assigned in input paramerter later.
!  real(double), parameter ::  rom_rho_min = 2.5d-13 ! [g/cm^3]

  
contains


  ! 
  ! CONSTRUCTOR
  !
  ! -- 1. initializes the romanova's data, 2. assign parameters
  subroutine init_romanova(this, Rs, Rmax, Mass, isothermal, T_flow, &
       r_ref, rho_ref, T_ref, v_ref, tilt, period, filename)
    type(romanova), intent(inout) :: this
    real(double), intent(in) :: Rs          ! radius of central star  [10^10cm
    real(double), intent(in) :: Rmax        ! max radius of geometry  [10^10cm]
    real(double), intent(in) :: Mass        ! [g]  Mass of the star
    logical,      intent(in) :: isothermal  ! If T, uses T_flow; oterwise uses data
    real(double), intent(in) :: T_flow      ! [K]  Isothermal temperature of flow
    !
    real(double), intent(in) :: r_ref       ! [cm] Reference length value
    real(double), intent(in) :: rho_ref     ! [g/cm^3]  Reference density value
    real(double), intent(in) :: T_ref       ! [K]       Reference temperature
    real(double), intent(in) :: v_ref       ! [cm/s]    Reference speed
    real(double), intent(in) :: tilt        ! [rad] Tilt angle of the magnetic axis
    real(double), intent(in) :: Period      ! [sec] rotational period
    !
    character(LEN=*), intent(in) :: filename   ! File name containig the data

    !
    ! Allocating memory for data and reading the data
    ! -- using a routine in irreg_grid_data_class
    call read_romanova_data(this%romanova_output, filename)
    this%data_alive = .true.  

    ! converting units of the data in this%romanova_data    
    ! --using routine in zeus_data_class
    call scale_coordinates(this%romanova_output, r_ref*1.0d-10)      !  [10^10 cm] after scaling
    call scale_density(this%romanova_output, rho_ref)                !  [g/cm^3]
    call scale_temperature(this%romanova_output, T_ref)              !  [K]
    call scale_velocity(this%romanova_output, v_ref)                 !  [cm/s] in speed of light
    !
    !Rmax = get_r(this%romanova_output, get_nr(this%romanova_output))

    !
    ! Save parameters in romanova_data objects
    this%Rs = Rs                  ! [10^10cm]
    this%Rmax = Rmax              ! [10^10cm] of the original data but not the model space of our model
    this%isothermal = isothermal  !
    this%T_flow = T_flow          ! [K]
    this%r_ref = r_ref            ! [cm]
    this%rho_ref = rho_ref        ! 
    this%T_ref = T_ref            ! 
    this%v_ref = v_ref
    this%tilt = tilt
    this%Period = Period

  end subroutine init_romanova



  !
  ! -- Deallocate memory for the romanova data 
  !    but keeping the rest of the data
  subroutine deallocate_romanova_data(this)
    type(romanova), intent(inout) :: this
    !
    ! using a routine in irreg_grid_data_class
    call kill(this%romanova_output)
    this%data_alive = .false.  
  end subroutine deallocate_romanova_data



  !
  ! access parameters
  !
  function get_dble_parameter(this,name)  RESULT(out)
    implicit none    
    real(double) :: out
    type(romanova), intent(in) :: this
    character(LEN=*) :: name
    
    select case(name)
    case ("Rs")
       out = this%Rs
    case ("Rmax")
       out = this%Rmax
    case ("Mass")
       out = this%Mass
    case ("T_flow")
       out = this%T_flow
    case ("r_ref")
       out = this%r_ref
    case ("rho_ref")
       out = this%rho_ref
    case ("T_ref")
       out = this%T_ref
    case ("tilt")
       out = this%tilt
    case ("Period", "period")
       out = this%Period
    case ("v_ref")
       out = this%v_ref
    case default
       print *, 'Unknown name passed romanova_clas::get_dble_parameter.'
       print *, 'name =', TRIM(name)
       print *, 'Exiting the program ...'
       stop
    end select

  end function get_dble_parameter


  function get_log_parameter(this,name)  RESULT(out)
    implicit none    
    logical :: out
    type(romanova), intent(in) :: this
    character(LEN=*) :: name
    
    select case(name)
    case ("isothermal")
       out = this%isothermal
    case default
       print *, 'Unknown name passed romanova_clas::get_logical_parameter.'
       print *, 'name =', TRIM(name)
       print *, 'Exiting the program ...'
       stop
    end select

  end function get_log_parameter



  ! 
  ! To check if the memory if allocated for romanova_output
  logical function  romanova_data_alive(this)
    implicit none
    type(romanova), intent(in) :: this
    !
    romanova_data_alive = this%data_alive
    !
  end function romanova_data_alive



  !===================================================================
  ! Given a position vector in cartisian coordinates (x,y,z) and the tilt angle (tilt, 
  ! rotation around y axis), this routine find the components the points in 
  ! the original unlitlted spherical coordinates (r, phi, theta) of 
  ! the data original data ( R in R*, theta is polar angle in rad, phi 
  ! azimuth angle in rad)
  subroutine get_coordinates_romanova(this, position, r, phi, theta)
    
    IMPLICIT NONE

    type(romanova), intent(in)    :: this
    type(VECTOR), intent(in)  :: position       ! [10^10cm]
    real*8, intent(out) :: r, phi, theta ! r in R*, theta and phi are in rad, 
    !
    real*8 :: x, y, z       ! [10^10cm]
    real*8 :: xp, yp, zp, t 
    real*8, parameter :: PI = 3.141592653589793d0

    ! rotate the vector back to the original coodinates
    x = position%x;  y = position%y;  z = position%z

    !
    t = -this%tilt 
    yp  =  y
    xp  =  x*cos(t) + z*sin(t)
    zp  = -x*sin(t) + z*cos(t)

        
    ! Changing coordinates to spherical caoordinates
    
    r=SQRT(xp*xp+yp*yp+zp*zp)  ! [10^10cm]

    theta = ACOS(zp/r)    ! should be between 0 and PI
    phi = ATAN2(yp,xp)     ! is between -PI/2 and PI/2

    ! this part should be changed accodingly. 
    if (phi < 0.0) phi = phi + 2.0d0*pi  ! now between 0 and 2*pi
    !
    if (theta<0.0) then
       write(*,*) "Error:: theta < 0.0 in romanova_class::get_coordinates_romanova."
       write(*,*) "        theta = ", theta
       stop
    end if
    if (phi<0.0) then
       write(*,*) "Error:: phi < 0.0 in Romanova_class::get_coordinates_romanova."
       write(*,*) "        phi = ", phi
       stop
    end if
       
  end subroutine get_coordinates_romanova




  !===================================================================
  ! Given a position vector in cartisian coordinates (x,y,z),
  ! this routine find the components the points in 
  ! the spherical coordinates (r, phi, theta) 
  subroutine get_spherical_coordinates(position, r, phi, theta)
    
    IMPLICIT NONE

    type(VECTOR), intent(in)  :: position       ! [10^10cm]
    real*8, intent(out) :: r, phi, theta ! r in R*, theta and phi are in rad, 
    !
    real*8 :: x, y, z       ! [10^10cm]
    real*8, parameter :: PI = 3.141592653589793d0

    ! rotate the vector back to the original coodinates
    x = position%x;  y = position%y;  z = position%z

        
    ! Changing coordinates to spherical caoordinates
    
    r=SQRT(x*x+y*y+z*z)  ! [10^10cm]

    theta = ACOS(z/r)    ! should be between 0 and PI
    phi = ATAN2(y,x)     ! is between -PI/2 and PI/2

    ! this part should be changed accodingly. 
    if (phi < 0.0) phi = phi + 2.0d0*pi  ! now between 0 and 2*pi
    !
    if (theta<0.0) then
       write(*,*) "Error:: theta < 0.0 in romanova_class::get_spherical_coordinates."
       write(*,*) "        theta = ", theta
       stop
    end if
    if (phi<0.0) then
       write(*,*) "Error:: phi < 0.0 in Romanova_class::get_spherical_coordinates."
       write(*,*) "        phi = ", phi
       stop
    end if
       
  end subroutine get_spherical_coordinates




    
  !
  !
  !===================================================================
  ! Given a point and grid, this function returns the density in [g] 
  ! at the point by interpolating the romanova's original data.
  !    
  FUNCTION romanova_density(this, point) RESULT(out)
    
    IMPLICIT NONE

    REAL(double)                  :: out     ! [g/cm^3]
    type(romanova), intent(in)    :: this
    TYPE(VECTOR), INTENT(IN) :: point

    !
    real(double), parameter ::  rom_rho_min = 1.0d-17 ! [g/cm^3]     

    ! Point in spherical coordinates used in romanova's data.
    real(double) :: r, theta, phi ! r in R*, theta and phi are in rad, 


    ! Finding the spherical coordinates used in Romanova's data 
    ! file.
    call get_coordinates_romanova(this, point, r, phi, theta)



    ! now phi and theta are in the range of the orignal data, so do interpolations.    
    ! using a function in zeus_data_class
    if (r<this%Rs) then ! inside of the star 
       out = rom_rho_min  ! just assign a small value and return
    else       
 !      ! for testing =================================
 !      theta = 0.0d0; phi=0.0d0
 !      !==============================================

       out = interp_density(this%romanova_output, r, theta, phi)  ! in [g/cm^3] 

    end if

    ! In original data, the density inside the star may be set to a very high density
    ! or in case something goes wrong, we set the density to be small....
!    if (out > 1.0e-5) then
!    out = rom_rho_min  
!    end if
    if (out < rom_rho_min) then
       out = rom_rho_min  
    end if
    

  END FUNCTION romanova_density




  !===================================================================
  ! Given a point and grid, this function returns
  ! the wind velocity (as a vector) in units of speed of light [c]. 
  ! 
  !====================================================================
  FUNCTION romanova_velocity(this, point) RESULT(out)
    
    IMPLICIT NONE

    type(vector)                  :: out   ! [c]
    type(romanova), intent(in)    :: this
    TYPE(VECTOR), INTENT(IN) :: point
    type(vector) :: v, rot
    real :: Vrot
    real, parameter :: pi=3.14159265

    ! Point in spherical coordinates
    real(double) :: r, theta, phi ! r in R*, theta and phi are in rad, 


    ! Finding the spherical coordinates used in Romanova's data 
    ! file.
    call get_coordinates_romanova(this, point, r, phi, theta)


    ! using a function in zeus_data_class
    ! The result is in cartisian coordinates.
    if (r<this%Rs) then ! inside of the star
       v = VECTOR(1.0, 1.0, 1.0)
    else
!       ! for testing =================================
!       theta = 0.0d0
!       !==============================================

       v = interp_velocity(this%romanova_output, r, theta, phi)  ! [cm/s]

       ! adding rotational velocity componets
       ! Since the original velocity's in romanova's data is in 
       ! the frame of rotating star. 
       Vrot = 2.0*PI*REAL(this%Rs)*1.0e10/REAL(this%Period)  ! [cm/s]
       call get_spherical_coordinates(point, r, phi, theta)
       rot = VECTOR(-Vrot*SIN(phi), Vrot*COS(phi), 0.0)

       ! vector addition
       v = v + rot           ! [cm/s]
    end if


    ! changing units to speed of light
    out = v/cSpeed                                   ! [c]

    
  END FUNCTION Romanova_velocity




  !
  !
  !
  !===============================================================
  !  Returns given point and grid, this function return 
  !  the wind temperature in [K].
  !
  !===============================================================  
  FUNCTION romanova_temperature(this, point)  RESULT(out)
    
    IMPLICIT NONE
    REAL                          :: out 
    type(romanova), intent(in)    :: this
    TYPE(VECTOR), INTENT(IN) :: point
    !
    real(double) :: r, theta, phi
    ! The following should be set as a componet of this obejct
    ! later. 
    real(double), parameter :: Rdisc = 100.0d0  ! [10^10cm] 
    real(double), parameter :: Tdisc = 3000.d0  ! [K] tempearture of disc
    real(double), parameter :: Tmin = 3500.d0   ! [K] Minimum temperature
    
    ! Finding the spherical coordinates used in Romanova's data 
    ! file.
    call get_coordinates_romanova(this, point, r, phi, theta)


    ! This 1.05 is important, and it should be there.
    ! if you change this number, you should change it also
    ! in calc_romanova_mass_velocity
    if (r<this%Rs*1.05d0 .or. r > this%Rmax) then 
       ! inside of the star or radius too large
       out = Tdisc  ! for now is just isothermal flow.
    elseif ( (this%Rs < r .and. r < this%Rmax*1.05d0) .and. r < Rdisc) then 
       ! in funnel flow
       if (this%isothermal) then
          out = this%T_flow  ! for now is just isothermal flow.
       else
          out = interp_temperature(this%romanova_output, r, theta, phi)  ! in [g/cm^3] 
       end if
    else 
       ! in the disc 
       out = Tdisc  ! for now is just isothermal flow.
    end if

    ! 
    if (out < Tmin) out = Tmin


  END FUNCTION Romanova_temperature


  
  !
  ! The interface routine to be used in generic routine in amr_mod.f90
  !
  
  SUBROUTINE calc_romanova_mass_velocity(this, thisOctal,subcell)
    !  Assigings the density and velocity values to thisOctal at a given
    !  subcell using the density and the velocity functions as descrubed
    !  below.

    !    use input_variables
        
    IMPLICIT NONE

    type(romanova), intent(in)    :: this
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    !    
    TYPE(VECTOR) :: point    
    REAL(oct) :: r 
    real(double), parameter :: rho_contrast = 0.015d0
!    real(double), parameter :: rho_contrast = 0.060d0
    real(double), save :: rhomin, rhomax
    logical, save :: first_time = .true.
    real(double) :: rho

    if (first_time) then
       ! find the max value of the density
       ! -- using the routine in zeus_data_class.f90
       rhomax = max_density(this%romanova_output)
       ! setting the minimum density 
       ! All the cells which has smaller density than this value 
       ! will be flagged as .not. inflow.
       rhomin = rho_contrast*rhomax
!       rhomin = rom_rho_min  !using the value assigned near the top of this module.

       first_time = .false.
    end if

    point = subcellCentre(thisOctal,subcell)
    
    r = modulus( point )   ! [10^10cm]
                      
    ! test if the point lies within the star
    IF ( r > this%Rs*1.05d0 .AND. r < this%Rmax) THEN 
       ! This 1.05 is important, and it should be there.
       ! if you change this number, you should change it also
       ! in romanova_temperature
       thisOctal%inFlow(subcell) = .TRUE.       
       thisOctal%inStar(subcell) = .FALSE.
    ELSEIF (r <= this%Rs ) THEN
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%inStar(subcell) = .TRUE.
    ELSE
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%inStar(subcell) = .FALSE.            
    END IF

    ! assiging density using a function in this module    
    rho = romanova_density(this,point)   ! [g/cm^3]
    if (rho > rhomin) then
       thisOctal%rho(subcell)  = rho
    else
       thisOctal%rho(subcell)  = rhomin
       thisOctal%inFlow(subcell) = .FALSE.
    end if


    ! assiging velocity using a fuction in this module
    thisOctal%velocity(subcell) = romanova_velocity(this, point)  ! [c]

    ! assigining corner velocities using a function in this module
    IF (subcell == 8) CALL VelocityCorners(this, thisOctal)
    
    

  END SUBROUTINE calc_romanova_mass_velocity







  !
  ! The interface routine to be used in generic routine in amr_mod.f90
  !  
  SUBROUTINE calc_romanova_temperature(this, thisOctal,subcell)
  
    IMPLICIT NONE

    type(romanova), intent(in)    :: this
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    
    TYPE(VECTOR) :: point 

    point = subcellCentre(thisOctal,subcell)

    if ( thisOctal%inFlow(subcell) ) THEN
       ! using the function in this module    
       thisOctal%temperature(subcell) = romanova_temperature(this, point)  ! [K]
       ! we will initialise the bias distribution
       !  --- I am not sure why this should be done here....
       thisOctal%biasLine3D(subcell) = 1.0d0
       thisOctal%biasCont3D(subcell) = 1.0d0
    else
       thisOctal%biasLine3D(subcell) = 1.0d-150
       thisOctal%biasCont3D(subcell) = 1.0d-150
    end if

  END SUBROUTINE calc_romanova_temperature



  SUBROUTINE VelocityCorners(this, thisOctal)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE

    type(romanova), intent(in)    :: this
    TYPE(octal), INTENT(INOUT) :: thisOctal

    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    

    ! we first store the values we use to assemble the position vectors
    
    x1 = thisOctal%centre%x - thisOctal%subcellSize
    x2 = thisOctal%centre%x
    x3 = thisOctal%centre%x + thisOctal%subcellSize

    y1 = thisOctal%centre%y - thisOctal%subcellSize
    y2 = thisOctal%centre%y
    y3 = thisOctal%centre%y + thisOctal%subcellSize
    
    z1 = thisOctal%centre%z - thisOctal%subcellSize
    z2 = thisOctal%centre%z
    z3 = thisOctal%centre%z + thisOctal%subcellSize

    ! now store the 'base level' values
    
    thisOctal%cornerVelocity(1) = Romanova_velocity(this, VECTOR(x1,y1,z1))
    thisOctal%cornerVelocity(2) = Romanova_velocity(this, VECTOR(x2,y1,z1))
    thisOctal%cornerVelocity(3) = Romanova_velocity(this, VECTOR(x3,y1,z1))
    thisOctal%cornerVelocity(4) = Romanova_velocity(this, VECTOR(x1,y2,z1))
    thisOctal%cornerVelocity(5) = Romanova_velocity(this, VECTOR(x2,y2,z1))
    thisOctal%cornerVelocity(6) = Romanova_velocity(this, VECTOR(x3,y2,z1))
    thisOctal%cornerVelocity(7) = Romanova_velocity(this, VECTOR(x1,y3,z1))
    thisOctal%cornerVelocity(8) = Romanova_velocity(this, VECTOR(x2,y3,z1))
    thisOctal%cornerVelocity(9) = Romanova_velocity(this, VECTOR(x3,y3,z1))

    ! middle level
  
    thisOctal%cornerVelocity(10) = Romanova_velocity(this, VECTOR(x1,y1,z2))
    thisOctal%cornerVelocity(11) = Romanova_velocity(this, VECTOR(x2,y1,z2))
    thisOctal%cornerVelocity(12) = Romanova_velocity(this, VECTOR(x3,y1,z2))
    thisOctal%cornerVelocity(13) = Romanova_velocity(this, VECTOR(x1,y2,z2))
    thisOctal%cornerVelocity(14) = Romanova_velocity(this, VECTOR(x2,y2,z2))
    thisOctal%cornerVelocity(15) = Romanova_velocity(this, VECTOR(x3,y2,z2))
    thisOctal%cornerVelocity(16) = Romanova_velocity(this, VECTOR(x1,y3,z2))
    thisOctal%cornerVelocity(17) = Romanova_velocity(this, VECTOR(x2,y3,z2))
    thisOctal%cornerVelocity(18) = Romanova_velocity(this, VECTOR(x3,y3,z2))

    ! top level
    
    thisOctal%cornerVelocity(19) = Romanova_velocity(this, VECTOR(x1,y1,z3))
    thisOctal%cornerVelocity(20) = Romanova_velocity(this, VECTOR(x2,y1,z3))
    thisOctal%cornerVelocity(21) = Romanova_velocity(this, VECTOR(x3,y1,z3))
    thisOctal%cornerVelocity(22) = Romanova_velocity(this, VECTOR(x1,y2,z3))
    thisOctal%cornerVelocity(23) = Romanova_velocity(this, VECTOR(x2,y2,z3))
    thisOctal%cornerVelocity(24) = Romanova_velocity(this, VECTOR(x3,y2,z3))
    thisOctal%cornerVelocity(25) = Romanova_velocity(this, VECTOR(x1,y3,z3))
    thisOctal%cornerVelocity(26) = Romanova_velocity(this, VECTOR(x2,y3,z3))
    thisOctal%cornerVelocity(27) = Romanova_velocity(this, VECTOR(x3,y3,z3))
    
  END SUBROUTINE VelocityCorners





  
end module romanova_class



