module luc_cir3d_class
  ! --------------------------------------------------------------------------------
  !  Created: nov-08-2004  (R.Kurosawa)
  !---------------------------------------------------------------------------------

  use kind_mod
  use constants_mod, only: pi
  use octal_mod
  use zeus_data_class
  use vector_mod

  IMPLICIT NONE
    
  public :: &
       new, &
       deallocate_zeus_data, &
       get_dble_param, &
       zeus_data_alive, &
       luc_cir3d_density, &
       luc_cir3d_velocity, &
       calc_cir3d_mass_velocity, &
       calc_cir3d_temperature


  private :: &
       init_luc_cir3d, &
       spherical_coordinates, &
       VelocityCorners

  
  
  type luc_cir3d
     private
     real(double)    :: Rs          ! radius of central star  [10^10cm]
     real(double)    :: Rmax        ! cutoff radius in [10^10 cm]
     real(double)    :: Mass        ! [g]  Mass of the star
     real(double)    :: Teff        ! [K]  Effective temperature of the star
     real(double)    :: Luminosity  ! [erg/s]
     real(double)    :: Twind       ! [K]  Isothemal temperature of the stellar wind
     real(double)    :: Mdot_scale  ! [-]  Denisty is scaled by this factor
     type(zeus_data) :: zeus_output ! Data structure for zeus data
     ! flag to indicate if zeus_output memory for is allocated.
     logical         :: data_alive= .false.  
  end type luc_cir3d

  !
  !
  !
  interface new
     module procedure init_luc_cir3d
  end interface


  !
  ! One instance to be used ONLY IN THIS MODULE.
  ! This is a work around since to use Neil's generic rouines (calcValuesAMR and  finishGrid  
  ! rouine. 
  type(luc_cir3d), save :: cir3d_data  


  
contains


  ! 
  ! CONSTRUCTOR
  !
  ! -- 1. initializes the zeus data, 2. assign parameters
!  subroutine init_luc_cir3d(this, Rmin, Twind, rgrid_filename, rho_vel_filename)
  subroutine init_luc_cir3d(Rs, Twind, Mdot_scale, rgrid_filename, rho_vel_filename)
!    type(luc_cir3d), intent(inout) :: this
    real(double), intent(in) :: Rs        ! radius of central star  [10^10cm]
    real(double), intent(in) :: Twind       ! [K]  Isothermal temperature of the stellar wind
    real(double), intent(in) :: Mdot_scale  ! [-]  Denisty is scaled by this factor
    character(LEN=*), intent(in) :: rgrid_filename   ! File name containig r grid.
    character(LEN=*), intent(in) :: rho_vel_filename ! File name for velocity and density    
    !
    integer :: nr
    real(double) :: Rmax
    !
    ! Allocateing memory for data 
    ! -- using a routine in zeus_data_class
    call new(cir3d_data%zeus_output)
    
    ! reading the data from file
    ! -- using a routine in zeus_data_class
    call read_zeus_data(cir3d_data%zeus_output, rgrid_filename, rho_vel_filename)
    cir3d_data%data_alive = .true.  

    nr = get_nr(cir3d_data%zeus_output)
    Rmax = get_r(cir3d_data%zeus_output, nr)  ! [Rs]
    !
    ! Save parameters in cir3d_data objects
    cir3d_data%Rs = Rs               ! [10^10cm]
    cir3d_data%Rmax = Rs*Rmax        ! [10^10cm]
    cir3d_data%Twind = Twind
    cir3d_data%Mdot_scale = Mdot_scale
  end subroutine init_luc_cir3d



  !
  ! -- Deallocate memory for the zeus data 
  !
  subroutine deallocate_zeus_data()
!    type(luc_cir3d), intent(inout) :: this
    !
    ! using a routine in zeus_data_class
    call kill( cir3d_data%zeus_output )
    cir3d_data%data_alive = .false.  
  end subroutine deallocate_zeus_data



  !
  ! access parameters
  !
  function get_dble_param(this,name)  RESULT(out)
    implicit none    
    real(double) :: out
    type(luc_cir3d), intent(in) :: this
    character(LEN=*) :: name
    
    select case(name)
    case ("Rs")
       out = this%Rs
    case ("Rmax")
       out = this%Rmax
    case ("Twind")
       out = this%Twind
    case ("Mdot_scale")
       out = this%Mdot_scale
    case default
       print *, 'Unknown name passed luc_cir3d_clas::get_dble_param.'
       print *, 'name =', TRIM(name)
       print *, 'Exiting the program ...'
       stop
    end select

  end function get_dble_param



  ! 
  ! To check if the memory if allocated for zeus_output
  logical function  zeus_data_alive()
    implicit none
!    type(luc_cir3d), intent(in) :: this
    !
    zeus_data_alive = cir3d_data%data_alive
    !
  end function zeus_data_alive




  !===================================================================
  ! Given a point in vector object (cartisian coordinates), the this routine 
  ! returns the corresponding sperrical componets : 
  ! ( R in R*, theta in rad, phi in rad)
  ! Note: phi and theta are forced to be in [0,pi/2].
  subroutine spherical_coordinates(point, r, theta, phi, phi_original)
    
    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    real(double), intent(out) :: r, theta, phi ! r in R*, theta and phi are in rad, 
    real(double), optional, intent(out) :: phi_original ! Original phi value in rad, 
    ! Point in cartisian coordinates
    real(double) :: x, y, z ! [10^10cm]
!    real(double) :: tmp

    if (.not.zeus_data_alive()) then
       write(*,*) "Error: cir3d_data is not alive (not allocated). " 
       write(*,*) "       ... [luc_cir3d_class::spherical_coordinates]." 
       stop 
    end if
    

    ! First find the componets in the spherical must be found
    x = point%x; y=point%y; z=point%z  ! [10^10cm]
    
    r=SQRT(x*x+y*y+z*z)  ! [10^10cm]

    theta = ACOS(z/r)    
    phi = ATAN2(y,x)

    ! converting units to R_* (R_min)
    r = r/(cir3d_data%Rs)

    ! now do range checks 
    do while (theta >=pi/2.0) 
       theta = pi - theta
    end do
    ! Note: If the periodicity around the z-axis changes from pi/2 to pi
    ! this part should be changed accodingly. 
    if (phi < 0.0) phi = phi + 2.0*pi
    if (PRESENT(phi_original)) phi_original = phi  
    do while (phi >=pi/2.0) 
       phi = phi - pi/2.0  
    end do
    !
    if (theta<0.0) then
       write(*,*) "Error:: theta < 0.0 in luc_cir3d_class::spherical_coordinates."
       write(*,*) "        theta = ", theta
       stop
    end if
    if (phi<0.0) then
       write(*,*) "Error:: phi < 0.0 in Luc_cir3d_class::spherical_coodrinates."
       write(*,*) "        phi = ", phi
       stop
    end if
       
  end subroutine spherical_coordinates


    
  !
  !
  !===================================================================
  ! Given a point and grid, this function returns the density in [g] 
  ! at the point by interpolating the original zeus data.
  !    
  FUNCTION luc_cir3d_density(point) RESULT(out)
    
    IMPLICIT NONE

    REAL(double)                  :: out
    TYPE(octalVector), INTENT(IN) :: point

    !
    real(double), parameter ::  rho_min = 1.0e-18 ! [g/cm^3]     

    ! Point in spherical coordinates
    real(double) :: r, theta, phi ! r in R*, theta and phi are in rad, 

    ! Changing coodninates to spherical coordinates
    call spherical_coordinates(point, r,    theta,     phi)
    !                                 [R*]   [rad]     [rad]

    ! now phi and theta are in the range of the orignal data, so do interpolations.    
    ! using a function in zeus_data_class
!    if (r<1.001) then ! inside of the star 
    if (r<0.996178765880217587) then ! inside of the star 
       out = rho_min  ! just assign a small value and return
    else       
 !      ! for testing =================================
 !      theta = 0.0d0; phi=0.0d0
 !      !==============================================

       out = interp_density(cir3d_data%zeus_output, r, theta, phi)  ! in [g/cm^3] 
       ! multiplying it by the scale factor
       out = out*cir3d_data%Mdot_scale
    end if

    ! In original data, the density inside the star may be set to a very high density
    ! or in case something goes wrong, we set the density to be small....
    if (out > 1.0e-5) out = rho_min  
    



  END FUNCTION luc_cir3d_density


  !===================================================================
  ! Given a point and grid, this function returns
  ! the wind velocity (as a vector) in units of speed of light [c]. 
  ! 
  !====================================================================
  FUNCTION luc_cir3d_velocity(point) RESULT(out)
    
    IMPLICIT NONE

    type(vector)                  :: out   ! [c]
    TYPE(octalVector), INTENT(IN) :: point
    type(vector) :: v

!    real(double) :: Vr, Vtheta, Vphi

    ! Point in spherical coordinates
    real(double) :: r, theta, phi ! r in R*, theta and phi are in rad, 
    real(double) :: phi_original, phi_rot

    ! Changing coodninates to spherical coordinates
    call spherical_coordinates(point, r,    theta,     phi,     phi_original)
    !                                 [R*]   [rad]     [rad]       [rad]

    ! using a function in zeus_data_class
    ! The result is in cartisian coordinates.
!    if (r<1.001d0) then ! inside of the star
    if (r<0.996178765880217587) then ! inside of the star
       v = VECTOR(1.0, 1.0, 1.0)
    else
!       ! for testing =================================
!       theta = 0.0d0
!       !==============================================

       v = interp_velocity(cir3d_data%zeus_output, r, theta, phi)  ! [cm/s]
    end if
    ! now rotate the vector zeus data is only for phi = 0 -- pi/2
    phi_rot = phi_original - phi
    ! using a routine in vector_mod
    ! (pass rotation around z axis thus minus sign)
    out = rotateZ(v, REAL(-phi_rot))


    

    ! changing units to speed of light
    out = out/cSpeed                                  ! [c]


    ! Vz flips sign if z is negative.
    if (point%z<0.0) out%z = -1.0*out%z  
       
    
  END FUNCTION Luc_cir3d_velocity




  !
  !
  !
  !===============================================================
  !  Returns given point and grid, this function return 
  !  the wind temperature in [K].
  !
  !===============================================================  
  FUNCTION luc_cir3d_temperature(point)  RESULT(out)
    
    IMPLICIT NONE
    REAL                          :: out 
    TYPE(octalVector), INTENT(IN) :: point
    !
    real(double) :: r, theta, phi
    
    ! Changing coodninates to spherical coordinates
    call spherical_coordinates(point, r,    theta,     phi)
    !                                 [R*]   [rad]     [rad]

    out = cir3d_data%Twind  ! for now is just isothermal wind!

  END FUNCTION Luc_cir3d_temperature


  
  !
  ! The interface routine to be used in generic routine in amr_mod.f90
  !
  
  SUBROUTINE calc_cir3d_mass_velocity(thisOctal,subcell)
    !  Assigings the density and velocity values to thisOctal at a given
    !  subcell using the density and the velocity functions as descrubed
    !  below.

    !    use input_variables
        
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    !    
    TYPE(octalVector) :: point
    
    REAL(oct) :: r    
        

    point = subcellCentre(thisOctal,subcell)
    r = modulus( point )   ! [10^10cm]
                      
    ! test if the point lies within the star
    IF ( r > cir3d_data%Rs*1.001 .AND. r < cir3d_data%Rmax) THEN
       thisOctal%inFlow(subcell) = .TRUE.       
       thisOctal%inStar(subcell) = .FALSE.
    ELSEIF (r <= cir3d_data%Rs ) THEN
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%inStar(subcell) = .TRUE.
    ELSE
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%inStar(subcell) = .FALSE.            
    END IF

    ! assiging density using a function in this module    
    thisOctal%rho(subcell) = luc_cir3d_density(point)   ! [g/cm^3]

    ! assiging velocity using a fuction in this module
    thisOctal%velocity(subcell) = luc_cir3d_velocity(point)  ! [c]

    ! assigining corner velocities using a function in this module
    IF (subcell == 8) CALL VelocityCorners(thisOctal)
    
    

  END SUBROUTINE calc_cir3d_mass_velocity


  !
  ! The interface routine to be used in generic routine in amr_mod.f90
  !  
  SUBROUTINE calc_cir3d_temperature(thisOctal,subcell)
  
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    
    TYPE(octalVector) :: point 

    point = subcellCentre(thisOctal,subcell)

    ! using the function in this module    
    thisOctal%temperature(subcell) = luc_cir3d_temperature(point)  ! [K]
    ! we will initialise the bias distribution
    !  --- I am not sure why this should be done here....
    thisOctal%biasLine3D(subcell) = 1.0
    
  END SUBROUTINE calc_cir3d_temperature



  SUBROUTINE VelocityCorners(thisOctal)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
  
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
    
    thisOctal%cornerVelocity(1) = Luc_cir3d_velocity(octalVector(x1,y1,z1))
    thisOctal%cornerVelocity(2) = Luc_cir3d_velocity(octalVector(x2,y1,z1))
    thisOctal%cornerVelocity(3) = Luc_cir3d_velocity(octalVector(x3,y1,z1))
    thisOctal%cornerVelocity(4) = Luc_cir3d_velocity(octalVector(x1,y2,z1))
    thisOctal%cornerVelocity(5) = Luc_cir3d_velocity(octalVector(x2,y2,z1))
    thisOctal%cornerVelocity(6) = Luc_cir3d_velocity(octalVector(x3,y2,z1))
    thisOctal%cornerVelocity(7) = Luc_cir3d_velocity(octalVector(x1,y3,z1))
    thisOctal%cornerVelocity(8) = Luc_cir3d_velocity(octalVector(x2,y3,z1))
    thisOctal%cornerVelocity(9) = Luc_cir3d_velocity(octalVector(x3,y3,z1))

    ! middle level
  
    thisOctal%cornerVelocity(10) = Luc_cir3d_velocity(octalVector(x1,y1,z2))
    thisOctal%cornerVelocity(11) = Luc_cir3d_velocity(octalVector(x2,y1,z2))
    thisOctal%cornerVelocity(12) = Luc_cir3d_velocity(octalVector(x3,y1,z2))
    thisOctal%cornerVelocity(13) = Luc_cir3d_velocity(octalVector(x1,y2,z2))
    thisOctal%cornerVelocity(14) = Luc_cir3d_velocity(octalVector(x2,y2,z2))
    thisOctal%cornerVelocity(15) = Luc_cir3d_velocity(octalVector(x3,y2,z2))
    thisOctal%cornerVelocity(16) = Luc_cir3d_velocity(octalVector(x1,y3,z2))
    thisOctal%cornerVelocity(17) = Luc_cir3d_velocity(octalVector(x2,y3,z2))
    thisOctal%cornerVelocity(18) = Luc_cir3d_velocity(octalVector(x3,y3,z2))

    ! top level
    
    thisOctal%cornerVelocity(19) = Luc_cir3d_velocity(octalVector(x1,y1,z3))
    thisOctal%cornerVelocity(20) = Luc_cir3d_velocity(octalVector(x2,y1,z3))
    thisOctal%cornerVelocity(21) = Luc_cir3d_velocity(octalVector(x3,y1,z3))
    thisOctal%cornerVelocity(22) = Luc_cir3d_velocity(octalVector(x1,y2,z3))
    thisOctal%cornerVelocity(23) = Luc_cir3d_velocity(octalVector(x2,y2,z3))
    thisOctal%cornerVelocity(24) = Luc_cir3d_velocity(octalVector(x3,y2,z3))
    thisOctal%cornerVelocity(25) = Luc_cir3d_velocity(octalVector(x1,y3,z3))
    thisOctal%cornerVelocity(26) = Luc_cir3d_velocity(octalVector(x2,y3,z3))
    thisOctal%cornerVelocity(27) = Luc_cir3d_velocity(octalVector(x3,y3,z3))
    
  END SUBROUTINE VelocityCorners

  
  
end module luc_cir3d_class



