module density_mod

  !
  ! Generic function to return a density for a given point
  ! and geometry type.
  !
  ! If you add a new geometry, please don't forget to add the name
  ! in the list in 
  !
  
  use vector_mod
  use gridtype_mod
  use jets_mod
  use wr104_mod
  
  public :: density, TTauriDensity
  ! the specific definition of density functions should really be private, 
  ! but for now they are public... Better yet, they should be in their own module.
  ! See jets_mod for example. 
  
  private :: print_geometry_list
  
contains

  !
  !  For a given point (as a vector) with each component in
  !  10^10 cm and gridtype object, it will return the density in g/cm^3
  function density(r_vec, grid) RESULT(out)
    implicit none
    double precision :: out
    type(octalVector), intent(in) :: r_vec
    type(gridtype), intent(in) :: grid

    select case (grid%geometry)
       
    case ("jets")
       ! using a routine in jets_mod.f90
       out = JetsDensity(r_vec, grid)*1000.d0 
       !                              ^^^^^^^
       !                          Converting kg/m^3 to g/cm^3
       ! For debugging, 
       out = out*1.0d20

    case ("ttauri")
       !  using a routine this module
       out = TTauriDensity(r_vec, grid) ! [g/cm^3]
       ! Double check the units

    case("testamr")
       out = testDensity(r_vec,grid)


       
    case default
       print *, 'Error:: Geometry option passed to [density_mod::density] '&
            &' (through grid) has not been implemented yet. Sorry.'
       print *, 'Your geometry = ', grid%geometry
       print *, 'Avilable gemetries are: '
       call print_geometry_list()
       print *, 'Terminating the program ...'
       stop 
    end select

  end function density
    


  !
  ! List the available geometry
  ! 
  subroutine print_geometry_list()    
    implicit none
    ! # of option available currently.
    integer, parameter :: n = 4
    character(LEN=30) :: name(n)
    integer :: i

    ! List here
    name(1) =  'jets'
    name(2) =  'wr104'
    name(3) =  'ttauri'
    name(4) =  'testamr'
    
    do i = 1, n
       write(*, *)  '   ', i, '. ', name(i)
    end do
       
  end subroutine print_geometry_list


  !
  !
  !
  PURE FUNCTION TTauriDensity(point,grid) RESULT(rho)
    ! calculates the density at a given point for a model of a T Tauri 
    !   star with magnetospheric accretion
    !   see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    use parameters_mod

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    REAL                          :: rho

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec

    REAL :: r, rM, theta, y

    starPosn = grid%starPos1
    pointVec = (point - starPosn) * 1.e10_oc
    r = modulus( pointVec ) 

    ! test if the point lies within the star
    IF ( r < TTauriRstar ) THEN
      rho = 1.e-25
      RETURN
    END IF
    
    ! test if the point lies too close to the disk
    IF ( ABS(pointVec%z) < 4.0e-2 * TTauriRstar) THEN
      ! we will fake the coordinates so that the density
      !  remains constant within this region 
      pointVec%z = SIGN( 4.0e-2 * TTauriRstar, REAL(pointVec%z) )
      r = modulus( pointVec ) 
    END IF
   
    theta = ACOS( pointVec%z  / r )
    IF (ABS(MODULO(theta,pi)) > 1.e-10 ) THEN 
      rM  = r / SIN(theta)**2
    ELSE
      rM = HUGE(rM)
    END IF
     
    ! test if the point lies outside the accretion stream
    IF  ((rM > TTauriRinner) .AND. (rM < TTauriRouter )) THEN
      y = SIN(theta)**2 
  
      rho = (TTauriMdot * TTauriRstar) / (4.0 * pi * &
              (TTauriRStar/TTauriRinner - TTauriRstar/TTauriRouter))
      rho = rho * r**(-5.0/2.0) / SQRT( 2.0 * bigG * TTauriMstar ) 
      rho = rho * SQRT( 4.0 - 3.0*y) / SQRT( 1.0 - y) 
    ELSE
      rho = 1.e-25
    END IF
    
  END FUNCTION TTauriDensity


  function testDensity(point, grid)
    use input_variables
    real :: testDensity
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real :: r
    r = modulus(point)
    testDensity = 1.e-30
    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       testDensity = rho * (grid%rInner / r)**2
    endif
  end function testDensity

end module density_mod 
