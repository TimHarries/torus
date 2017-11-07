module analytical_velocity_mod

  use vector_mod
  use constants_mod
  use magnetic_mod
  use density_mod
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
       case(4) ! whitney
          analyticalVelocity = WhitneyVelocity(point)
       case(5) ! whitney
          analyticalVelocity = ulrichVelocity(point)



       case DEFAULT
          call writeFatal("Error in logic in amrGridVelocity")
       end select
     end function analyticalVelocity

     type(VECTOR) function ulrichVelocity(point)
       use utils_mod, only : ccubsolv
       use inputs_mod, only : erinner, erouter
       use inputs_mod, only : mcore
       type(VECTOR) :: point
       real(Double) :: vr, vt, vp, mu, r, mu_0
       real(double) :: costheta, sintheta, sinphi, cosphi, sintheta0
       real(double) :: vx, vy, vz, v
       complex(double) :: a(3),z(3)
       r = modulus(point)*1.d10

       ulrichVelocity = VECTOR(0.d0, 0.d0, 0.d0)
       mu = (point%z*1.e10) /r
       if ((r > erInner).and.(r < erOuter)) then

          a(1) = cmplx(-mu * (r/erinner),0.d0,kind=double)
          a(2) = cmplx((r/erinner-1.d0), 0.d0,kind=double)
          a(3) = 0.d0
          call ccubsolv(a,z)
          mu_0 = real(z(1))
          if ((imag(z(1)) /= 0.d0).or.(imag(z(2))==0.d0).or.(imag(z(3))==0.d0)) then
             write(*,*)"problem with cubic solver"
          endif
          sintheta0 = sqrt(1.d0 - mu_0**2)
          costheta = point%z*1.d10/r
          sintheta = sqrt(1.d0- costheta**2)
          cosphi = point%x/sqrt(point%x**2 + point%y**2)
          sinphi = point%y/sqrt(point%x**2 + point%y**2)
          
          v = sqrt(bigG * mcore/r)/cSpeed

          vr = -v * sqrt(1.d0 + costheta/mu_0)
          vt = v * (mu_0 - costheta) * sqrt(abs((mu_0+costheta)/(mu_0*costheta)))
          vp = v * (sintheta0/sintheta) * sqrt(1.d0-costheta/mu_0)

          vx = vr * cosphi * sintheta + vt * cosphi * costheta - vp * sinphi
          vy = vr * sinphi * sintheta + vt * sinphi * costheta + vp * cosphi
          vz = vr * costheta - vt  * sintheta

          ulrichVelocity = VECTOR(vx, vy, vz)



       endif

     end function ulrichVelocity


     type(VECTOR) function TTauriStellarWindVelocity(point)
       use inputs_mod, only : ttauriMstar, ttaurirStar, SW_vMin, &
            SW_vMax, SW_Rmin, SW_beta, SW_rmax, SW_protation, SW_Veq

       type(VECTOR), intent(in) :: point
       type(VECTOR) :: rHat, phiHat
       real(double) :: r, v_r, vEsc, vterm, v_eq, v_phi, x, y, z, phi
       real(double) :: theta
       type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.0d0)

       r = modulus(point)
       x=point%x; y=point%y; z=point%z

       if ((r > SW_rMin).and.(r < SW_rMax)) then
          vEsc = sqrt(2.d0 * bigG * ttauriMstar / ttaurirStar)
          vterm = SW_vMax * vEsc
          v_r = SW_vMin + (vTerm - SW_Vmin) * (1.d0 - SW_Rmin/r)**SW_beta
          if (SW_Protation > 0.d0) then
             v_eq = twoPi * ttaurirstar / sw_protation ! [cm/s]
          else
             v_eq = SW_Veq
          endif
          phi=ATAN2(y,x)
          theta=ACOS(z/r)
          if (phi < 0.d0) phi=phi+twoPi
!          print *, "x", x, "y", y, "z", z, "r", r
          v_phi = v_eq * (SW_Rmin / r) * sin(theta)

          rHat = point
          call normalize(rHat)

          phiHat = point.cross.zAxis
          if (modulus(phiHat) /= 0.d0) then
             call normalize(phiHat)
             TTauriStellarWindVelocity = (v_r/cspeed) * rHat + (v_phi/cspeed) * phiHat
          else
             TTauriStellarWindVelocity = (v_r/cspeed) * rHat
          endif

       else
          TTauriStellarWindVelocity = VECTOR(0.d0, 0.d0, 0.d0)
       endif

     end function TTauriStellarWindVelocity

   end module analytical_velocity_mod
