module analytical_velocity_mod

  use vector_mod
  implicit none
  

contains

  type(VECTOR) function analyticalVelocity(point, i)
    use discwind_class,only : globalDiscWind, discwind_velocity
    type(VECTOR) :: point
    integer :: i

    analyticalVelocity = VECTOR(0.d0, 0.d0, 0.d0)

    select case(i)
       case(0)
          call writeFatal("Analytical velocity called with type of zero")
       case(1) ! disc wind
          analyticalVelocity =  discwind_velocity(globalDiscWind, point)
       case DEFAULT
          call writeFatal("Error in logic in amrGridVelocity")
       end select
     end function analyticalVelocity


   end module analytical_velocity_mod
