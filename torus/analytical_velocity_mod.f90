module analytical_velocity_mod

  use vector_mod
  use constants_mod
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
            SW_vMax, SW_Rmin, SW_beta, SW_rmax, SW_protation
       type(VECTOR) :: point, rHat, phiHat
       real(double) :: r, v_r, vEsc, vterm, v_eq, v_phi, x, y, z, phi
       real(double) :: vx, vy, vz, theta
       type(VECTOR) :: vHat
       type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.0d0)

       r = modulus(point)
       x=point%x; y=point%y; z=point%z

       if ((r > SW_rMin).and.(r < SW_rMax)) then
          vEsc = sqrt(2.d0 * bigG * ttauriMstar / ttaurirStar)
          vterm = SW_vMax * vEsc
          v_r = SW_vMin + (vTerm - SW_Vmin) * (1.d0 - SW_Rmin/r)**SW_beta
          v_eq = twoPi * ttaurirstar / sw_protation ! [cm/s]
          phi=ATAN2(y,x)
          theta=ACOS(z/r)
          if (phi < 0.d0) phi=phi+twoPi
!          print *, "x", x, "y", y, "z", z, "r", r
          v_phi = v_eq * (SW_Rmin / r) * sin(theta)

          phiHat = point.cross.zAxis
          call normalize(phiHat)
!       Vx = (V_r*SIN(theta)*COS(phi) - V_phi*SIN(phi))/cSpeed
!       Vy = (+1.d0*V_r*SIN(theta)*SIN(phi) + V_phi*COS(phi))/cSpeed
!       Vz = (V_r*COS(theta))/cSpeed
!       if (z < 0.d0) Vz = -Vz
!       print *, "Vx", Vx
!       print *, "Vy", Vy
!       print *, "Vz", Vz
!          print *, "cspeed", cSpeed, "Vesc", vEsc, "vterm", vterm, "v_eq", v_eq
!          print *, "theta", theta
!          print *, "V_r", v_r, "V_r / c", v_r/cSpeed
!          print *, "v_phi", v_phi, "Vphi /c", v_phi/cSpeed

          
          rHat = point
          call normalize(rHat)
!          TTauriStellarWindVelocity = (v_r/cSpeed) * vHat
          TTauriStellarWindVelocity = (v_r/cspeed) * rHat + (v_phi/cspeed) * phiHat
       else
          TTauriStellarWindVelocity = VECTOR(0.d0, 0.d0, 0.d0)
       endif

     end function TTauriStellarWindVelocity

   end module analytical_velocity_mod
