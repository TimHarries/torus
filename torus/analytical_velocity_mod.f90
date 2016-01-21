module analytical_velocity_mod

  use vector_mod
  use magnetic_mod
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
       case(2) ! ttauri keplerian
          analyticalVelocity =   TTauriKeplerianVelocity(point)
       case(3) ! ttauri stellar wind
          analyticalVelocity = TTauriStellarWindVelocity(point)



       case DEFAULT
          call writeFatal("Error in logic in amrGridVelocity")
       end select
     end function analyticalVelocity




     type(VECTOR) function TTauriStellarWindVelocity(point)
       use inputs_mod, only : ttauriMstar, ttaurirStar, SW_vMin, &
            SW_vMax, SW_Rmin, SW_beta, SW_rmax
       type(VECTOR) :: point
       real(double) :: r, v, vEsc, vterm
       type(VECTOR) :: vHat

       r = modulus(point)
       if ((r > SW_rMin).and.(r < SW_rMax)) then
          vEsc = sqrt(2.d0 * bigG * ttauriMstar / ttaurirStar)
          vterm = SW_vMax * vEsc
          v = SW_vMin + (vTerm - SW_Vmin) * (1.d0 - SW_Rmin/r)**SW_beta
          vHat = point
          call normalize(vHat)
          TTauriStellarWindVelocity = (v/cSpeed) * vHat
       else
          TTauriStellarWindVelocity = VECTOR(0.d0, 0.d0, 0.d0)
       endif

     end function TTauriStellarWindVelocity



   end module analytical_velocity_mod
