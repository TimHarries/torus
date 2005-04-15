module cmfgen_class
  ! --------------------------------------------------------------------------------
  !  This module uses the OPACITY_DATA output file from CMFGEN of J. Hillier to
  !  to map the opacity and denisty (it does not matter, but assuming proportional 
  !  to the square root of the thermal opacity).
  !  AND MORE COMMENTS TO BE INSERTED LATER.
  !---------------------------------------------------------------------------------
  !  Created: mar-24-2005  (R.Kurosawa)
  !---------------------------------------------------------------------------------

  use kind_mod
  use constants_mod, only: pi, cspeed_dbl
  use octal_mod
  use vector_mod
  use utils_mod

  IMPLICIT NONE
    
  public :: &
       read_cmfgen_data, &
       deallocate_cmfgen_data, &
       cmfgen_density, &
       cmfgen_velocity, &
       cmfgen_mass_velocity, &
       calc_cmfgen_temperature, &
       get_cmfgen_nd,&
       get_cmfgen_data_array, &
       get_cmfgen_data_array_element, &
       put_cmfgen_Rmin, &
       get_cmfgen_Rmin, &
       put_cmfgen_Rmax, &
       get_cmfgen_Rmax
       


  private :: &
       VelocityCorners

  
  
  type cmfgen_data
     private
     real(double)      :: Rmin          ! Minimum radius  [10^10cm]
     real(double)      :: Rmax          ! Maximum radius  [10^10cm]
     ! header info in the OPACITY_DATA
     character(LEN=12) :: format_date ! Revised format date
     character(LEN=12) :: model_id    ! Model ID
     character(LEN=12) :: trans_id    ! Transition ID
     logical           :: diff_aprox  ! if T diffusiton approx use at inner boundary
     real(double)      :: lambda      ! wavelength of the line
     real(double)      :: amass       ! atomic mass
     real(double)      :: Ic          ! Schuster intensity
     integer           :: nd          ! Number of depth points
     ! The following arrays should be allocated with size = nd
     real(double), pointer  :: R(:)      ! radial distance  [10^cm]
     real(double), pointer  :: T(:)      ! temperature [K]
     real(double), pointer  :: sigma(:)  ! 
     real(double), pointer  :: V(:)      ! radial velocity [km/s]
     real(double), pointer  :: eta(:)    ! thermal emissivity  [ ]
     real(double), pointer  :: chi_th(:) ! thermal opacity     [ ]
     real(double), pointer  :: esec(:)   ! electron scattering opacity []
     real(double), pointer  :: etal(:)   ! line emissivity []
     real(double), pointer  :: chil(:)   ! line emissivity []     
  end type cmfgen_data


  !
  ! One instance to be used ONLY IN THIS MODULE.
  ! This is a work around since to use Neil's generic rouines (calcValuesAMR and  finishGrid  
  ! rouine.  (This acts like an common block! This should be stopped in the future!)
  type(cmfgen_data), save :: cmfgen_opacity


  
contains


  ! 
  ! CONSTRUCTOR
  !
  ! -- Intialize the object and read the data.
  subroutine read_cmfgen_data(filename)
    implicit none
    character(LEN=*), intent(in) :: filename ! File name for velocity and density    
    !
    character(LEN=12) :: dum_a
    character(LEN=120) :: dum_a_long
    integer :: i, nd
    integer :: luin =88
    real(double) :: freq

    open(unit=luin, file=filename, status="old")

    ! Reading header and writing it to the standard output for record.
    write(*,'(a)') " "
    write(*,'(a)') " "
    write(*,'(a)') "======= Reading OPACITY_DATA of CMFGEN =================================== "
    do i = 1, 9
       read(luin,'(a)') dum_a_long
       write(*,'(a)') dum_a_long
    end do
    write(*,'(a)') "============================================================================== "
    write(*,'(a)') " "
    write(*,'(a)') " "


    ! now reading the data and storing them
    rewind(luin)
    read(luin,*) cmfgen_opacity%format_date
    read(luin,*) cmfgen_opacity%model_id
    read(luin,*) cmfgen_opacity%trans_id
    read(luin,*) cmfgen_opacity%diff_aprox
    read(luin,*) freq
    read(luin,*) cmfgen_opacity%lambda
    read(luin,*) cmfgen_opacity%amass
    read(luin,*) cmfgen_opacity%Ic
    read(luin,*) cmfgen_opacity%nd
    read(luin,*) dum_a  ! skip header    
    ! now allocate memory for data arrays
    nd = cmfgen_opacity%nd
    ALLOCATE(cmfgen_opacity%r(nd), cmfgen_opacity%T(nd), cmfgen_opacity%sigma(nd), cmfgen_opacity%V(nd))
    ALLOCATE(cmfgen_opacity%eta(nd), cmfgen_opacity%chi_th(nd), cmfgen_opacity%esec(nd), &
         cmfgen_opacity%etal(nd), cmfgen_opacity%chil(nd))

    do i = 1, nd
       read(luin, *)              &
            cmfgen_opacity%r(nd+1-i),       &
            cmfgen_opacity%T(nd+1-i),       &
            cmfgen_opacity%sigma(nd+1-i),   &
            cmfgen_opacity%V(nd+1-i),       &
            cmfgen_opacity%eta(nd+1-i),     &
            cmfgen_opacity%chi_th(nd+1-i),  &
            cmfgen_opacity%esec(nd+1-i),    &
            cmfgen_opacity%etal(nd+1-i),    &
            cmfgen_opacity%chil(nd+1-i)
    end do

    ! changing units of T from 10^4 K to K
    cmfgen_opacity%T(1:nd) = 1.0d4*cmfgen_opacity%T(1:nd)

    ! to be consistent with the notation in stateq_mod::generateOpacitiesAMR
    cmfgen_opacity%eta(1:nd) = 1.0d10*cmfgen_opacity%eta(1:nd)
    cmfgen_opacity%etal(1:nd) = 1.0d10*cmfgen_opacity%etal(1:nd)

    close(luin)

  end subroutine read_cmfgen_data
  


  !
  ! -- Deallocate memory for the data arrays
  !
  subroutine deallocate_cmfgen_data()
    implicit none
    !
    if (ASSOCIATED(cmfgen_opacity%r)) DEALLOCATE(cmfgen_opacity%r); NULLIFY(cmfgen_opacity%r)
    if (ASSOCIATED(cmfgen_opacity%T)) DEALLOCATE(cmfgen_opacity%T); NULLIFY(cmfgen_opacity%T)
    if (ASSOCIATED(cmfgen_opacity%sigma)) DEALLOCATE(cmfgen_opacity%sigma); NULLIFY(cmfgen_opacity%sigma)
    if (ASSOCIATED(cmfgen_opacity%V)) DEALLOCATE(cmfgen_opacity%V); NULLIFY(cmfgen_opacity%V)
    if (ASSOCIATED(cmfgen_opacity%eta)) DEALLOCATE(cmfgen_opacity%eta); NULLIFY(cmfgen_opacity%eta)
    if (ASSOCIATED(cmfgen_opacity%chi_th)) DEALLOCATE(cmfgen_opacity%chi_th); NULLIFY(cmfgen_opacity%chi_th)
    if (ASSOCIATED(cmfgen_opacity%esec)) DEALLOCATE(cmfgen_opacity%esec); NULLIFY(cmfgen_opacity%esec)
    if (ASSOCIATED(cmfgen_opacity%etal)) DEALLOCATE(cmfgen_opacity%etal); NULLIFY(cmfgen_opacity%etal)
    if (ASSOCIATED(cmfgen_opacity%chil)) DEALLOCATE(cmfgen_opacity%chil); NULLIFY(cmfgen_opacity%chil)
    !
  end subroutine deallocate_cmfgen_data



  !
  ! ACCESSORS
  !
  !

  ! reruns the number depth
  function get_cmfgen_nd() RESULT(out)
    integer :: out
    out = cmfgen_opacity%nd
  end function get_cmfgen_nd


  ! rerturn named data array
  subroutine get_cmfgen_data_array(name, array)
    implicit none
    character(LEN=*), intent(in) :: name 
    real(double), intent(inout) :: array(:)
    integer :: nd
    nd = size(array)
    if (nd< cmfgen_opacity%nd) then
       print *, "Error:: nd < cmfgen_opacity%nd in [cmfgen_class::get_cmfgen_data_array]."
       print *, "                       nd = ", nd
       print *, "        cmfgen_opacity%nd = ", cmfgen_opacity%nd
       stop
    end if
    nd = cmfgen_opacity%nd
    !
    select case (name)
    case ("R")
       array(1:nd) = cmfgen_opacity%R(1:nd)
    case ("T")
       array(1:nd) = cmfgen_opacity%T(1:nd)
    case ("sigma")
       array(1:nd) = cmfgen_opacity%sigma(1:nd)
    case ("V")
       array(1:nd) = cmfgen_opacity%V(1:nd)
    case ("eta")
       array(1:nd) = cmfgen_opacity%eta(1:nd)
    case ("chi_th")
       array(1:nd) = cmfgen_opacity%chi_th(1:nd)
    case ("esec")
       array(1:nd) = cmfgen_opacity%esec(1:nd)
    case ("etal")
       array(1:nd) = cmfgen_opacity%etal(1:nd)
    case ("chil")
       array(1:nd) = cmfgen_opacity%chil(1:nd)
    case default
       print *, "Error:: Unknown name passed to [cmfgen_class::get_cmfgen_data_array]."
       print *, "       name = ", name
       stop
    end select

  end subroutine get_cmfgen_data_array

  !
  ! get i-th componet of a named array
  function get_cmfgen_data_array_element(name, i) RESULT(out)
    implicit none
    real(double) :: out
    character(LEN=*), intent(in) :: name 
    integer, intent(in) :: i 

    !
    select case (name)
    case ("R")
       out = cmfgen_opacity%R(i)
    case ("T")
       out = cmfgen_opacity%T(i)
    case ("sigma")
       out = cmfgen_opacity%sigma(i)
    case ("V")
       out = cmfgen_opacity%V(i)
    case ("chi_th")
       out = cmfgen_opacity%chi_th(i)
    case ("esec")
       out = cmfgen_opacity%esec(i)
    case ("etal")
       out = cmfgen_opacity%etal(i)
    case ("chil")
       out = cmfgen_opacity%chil(i)
    case default
       print *, "Error:: Unknown name passed to [cmfgen_class::get_cmfgen_data_array]."
       print *, "       name = ", name
       stop
    end select

  end function get_cmfgen_data_array_element



       
  !
  !
  !===================================================================
  ! Given a point and grid, this function returns the density in [g] 
  ! at the point by interpolating the original zeus data.
  !    
  FUNCTION cmfgen_density(point) RESULT(out)
    
    IMPLICIT NONE

    REAL(double)                  :: out
    TYPE(octalVector), INTENT(IN) :: point    
    !
    real(double), parameter ::  rho_min = 1.0e-18 ! [g/cm^3]     
    real(double), parameter ::  rho_max = 1.0e-2  ! [g/cm^3]     
    !
    real(double) :: ri

    ri = MODULUS(point)

    if (ri<cmfgen_opacity%Rmin .or. ri>cmfgen_opacity%Rmax) then 
       out = rho_min  ! just assign a small value and return
    else       
       ! using a routine in utils_mod
       out = loginterp_dble(cmfgen_opacity%esec, cmfgen_opacity%nd, cmfgen_opacity%R, ri)    
    end if

  END FUNCTION cmfgen_density



  !===================================================================
  ! Given a point and grid, this function returns
  ! the wind velocity (as a vector) in units of speed of light [c]. 
  ! 
  !====================================================================
  FUNCTION cmfgen_velocity(point) RESULT(out)
    
    IMPLICIT NONE

    type(vector)                  :: out   ! [c]
    TYPE(octalVector), INTENT(IN) :: point
    type(vector) :: v

    ! Point in spherical coordinates
    real(double) :: r, Vr

    r = modulus(point)
    !
    if (r<cmfgen_opacity%Rmin) then ! inside of the star
       v = VECTOR(1.0e-10, 1.0e-10, 1.0e-10)  ! [c]
    else
       Vr = loginterp_dble(cmfgen_opacity%V, cmfgen_opacity%nd, cmfgen_opacity%R, r)    
       v = (point/r) * (Vr*1.0e5/cspeed_dbl)   ! [c]
       !     \hat(r) x  Vr/c
    end if

    out = v
    
  END FUNCTION cmfgen_velocity


  !
  !
  !
  !===============================================================
  !  Returns given point and grid, this function return 
  !  the wind temperature in [K].
  !
  !===============================================================  
  FUNCTION cmfgen_temperature(point)  RESULT(out)
    
    IMPLICIT NONE
    REAL                          :: out 
    TYPE(octalVector), INTENT(IN) :: point
    !
    real(double) :: r
    
    r = modulus(point)
    if (r<cmfgen_opacity%Rmin) then ! inside of the star
       out = 1000.0d0
    else
       out = loginterp_dble(cmfgen_opacity%T, cmfgen_opacity%nd, cmfgen_opacity%R, r)        
    end if

  END FUNCTION cmfgen_temperature


  
  !
  ! The interface routine to be used in generic routine in amr_mod.f90
  !
  
  SUBROUTINE cmfgen_mass_velocity(thisOctal,subcell)
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
    IF ( r >= cmfgen_opacity%Rmin .AND. r <= cmfgen_opacity%Rmax) THEN
       thisOctal%inFlow(subcell) = .TRUE.       
       thisOctal%inStar(subcell) = .FALSE.
    ELSEIF (r < cmfgen_opacity%Rmin) THEN
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%inStar(subcell) = .TRUE.
    ELSEIF (r > cmfgen_opacity%Rmax) THEN
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%inStar(subcell) = .FALSE.
    ELSE
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%inStar(subcell) = .FALSE.            
    END IF

    ! assiging density using a function in this module    
    thisOctal%rho(subcell) = cmfgen_density(point)   ! [g/cm^3]

    ! assiging velocity using a fuction in this module
    thisOctal%velocity(subcell) = cmfgen_velocity(point)  ! [c]

    ! assigining corner velocities using a function in this module
  IF ((thisoctal%threed).and.(subcell == 8)) &
       CALL VelocityCorners(thisOctal)
  IF ((thisoctal%twod).and.(subcell == 4)) &
       CALL VelocityCorners(thisOctal)
    

  END SUBROUTINE cmfgen_mass_velocity


  !
  ! The interface routine to be used in generic routine in amr_mod.f90
  !  
  SUBROUTINE calc_cmfgen_temperature(thisOctal,subcell)
  
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    
    TYPE(octalVector) :: point 

    point = subcellCentre(thisOctal,subcell)

    ! using the function in this module    
    thisOctal%temperature(subcell) = cmfgen_temperature(point)  ! [K]
    ! we will initialise the bias distribution
    !  --- I am not sure why this should be done here....
    thisOctal%biasLine3D(subcell) = 1.0
    
  END SUBROUTINE calc_cmfgen_temperature



  SUBROUTINE VelocityCorners(thisOctal)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
  
    TYPE(octal), INTENT(INOUT) :: thisOctal

    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    


    if (thisOctal%threeD) then
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
       
       thisOctal%cornerVelocity(1) = cmfgen_velocity(octalVector(x1,y1,z1))
       thisOctal%cornerVelocity(2) = cmfgen_velocity(octalVector(x2,y1,z1))
       thisOctal%cornerVelocity(3) = cmfgen_velocity(octalVector(x3,y1,z1))
       thisOctal%cornerVelocity(4) = cmfgen_velocity(octalVector(x1,y2,z1))
       thisOctal%cornerVelocity(5) = cmfgen_velocity(octalVector(x2,y2,z1))
       thisOctal%cornerVelocity(6) = cmfgen_velocity(octalVector(x3,y2,z1))
       thisOctal%cornerVelocity(7) = cmfgen_velocity(octalVector(x1,y3,z1))
       thisOctal%cornerVelocity(8) = cmfgen_velocity(octalVector(x2,y3,z1))
       thisOctal%cornerVelocity(9) = cmfgen_velocity(octalVector(x3,y3,z1))
       
       ! middle level
       
       thisOctal%cornerVelocity(10) = cmfgen_velocity(octalVector(x1,y1,z2))
       thisOctal%cornerVelocity(11) = cmfgen_velocity(octalVector(x2,y1,z2))
       thisOctal%cornerVelocity(12) = cmfgen_velocity(octalVector(x3,y1,z2))
       thisOctal%cornerVelocity(13) = cmfgen_velocity(octalVector(x1,y2,z2))
       thisOctal%cornerVelocity(14) = cmfgen_velocity(octalVector(x2,y2,z2))
       thisOctal%cornerVelocity(15) = cmfgen_velocity(octalVector(x3,y2,z2))
       thisOctal%cornerVelocity(16) = cmfgen_velocity(octalVector(x1,y3,z2))
       thisOctal%cornerVelocity(17) = cmfgen_velocity(octalVector(x2,y3,z2))
       thisOctal%cornerVelocity(18) = cmfgen_velocity(octalVector(x3,y3,z2))
       
       ! top level
       
       thisOctal%cornerVelocity(19) = cmfgen_velocity(octalVector(x1,y1,z3))
       thisOctal%cornerVelocity(20) = cmfgen_velocity(octalVector(x2,y1,z3))
       thisOctal%cornerVelocity(21) = cmfgen_velocity(octalVector(x3,y1,z3))
       thisOctal%cornerVelocity(22) = cmfgen_velocity(octalVector(x1,y2,z3))
       thisOctal%cornerVelocity(23) = cmfgen_velocity(octalVector(x2,y2,z3))
       thisOctal%cornerVelocity(24) = cmfgen_velocity(octalVector(x3,y2,z3))
       thisOctal%cornerVelocity(25) = cmfgen_velocity(octalVector(x1,y3,z3))
       thisOctal%cornerVelocity(26) = cmfgen_velocity(octalVector(x2,y3,z3))
       thisOctal%cornerVelocity(27) = cmfgen_velocity(octalVector(x3,y3,z3))
    
    else  ! 2D case
       
       ! we first store the values we use to assemble the position vectors
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values
       
       thisOctal%cornerVelocity(1) = cmfgen_velocity(octalVector(x1,0.d0,z1))
       thisOctal%cornerVelocity(2) = cmfgen_velocity(octalVector(x2,0.d0,z1))
       thisOctal%cornerVelocity(3) = cmfgen_velocity(octalVector(x3,0.d0,z1))
       thisOctal%cornerVelocity(4) = cmfgen_velocity(octalVector(x1,0.d0,z2))
       thisOctal%cornerVelocity(5) = cmfgen_velocity(octalVector(x2,0.d0,z2))
       thisOctal%cornerVelocity(6) = cmfgen_velocity(octalVector(x3,0.d0,z2))
       thisOctal%cornerVelocity(7) = cmfgen_velocity(octalVector(x1,0.d0,z3))
       thisOctal%cornerVelocity(8) = cmfgen_velocity(octalVector(x2,0.d0,z3))
       thisOctal%cornerVelocity(9) = cmfgen_velocity(octalVector(x3,0.d0,z3))
    endif

  END SUBROUTINE VelocityCorners


  !
  ! Stores Rmin value (this does not have to be same as R(1) value.
  !
  subroutine put_cmfgen_Rmin(Rmin)
    implicit none
    real(double), intent(in) :: Rmin
    cmfgen_opacity%Rmin = Rmin
  end subroutine put_cmfgen_Rmin

  !
  !
  !
  function get_cmfgen_Rmin() RESULT(Rmin)
    implicit none
    real(double) :: Rmin
    Rmin = cmfgen_opacity%Rmin
  end function get_cmfgen_Rmin


  subroutine put_cmfgen_Rmax(Rmax)
    implicit none
    real(double), intent(in) :: Rmax
    cmfgen_opacity%Rmax = Rmax
  end subroutine put_cmfgen_Rmax

  !
  !
  !
  function get_cmfgen_Rmax() RESULT(Rmax)
    implicit none
    real(double) :: Rmax
    Rmax = cmfgen_opacity%Rmax
  end function get_cmfgen_Rmax


end module cmfgen_class



