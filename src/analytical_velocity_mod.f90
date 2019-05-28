MODULE analytical_velocity_mod

  USE vector_mod
  USE constants_mod
  USE magnetic_mod
  USE density_mod
  IMPLICIT NONE


CONTAINS

  TYPE(VECTOR) FUNCTION analyticalVelocity(point, i)
    USE discwind_class,ONLY : globalDiscWind, discwind_velocity
    TYPE(VECTOR) :: point
    INTEGER :: i

    analyticalVelocity = VECTOR(0.d0, 0.d0, 0.d0)

    SELECT CASE(i)
    CASE(0)
       CALL writeFatal("Analytical velocity called with type of zero")
    CASE(1) ! disc wind
       analyticalVelocity =  discwind_velocity(globalDiscWind, point)
    CASE(2) ! ttauri keplerian
       analyticalVelocity =   TTauriKeplerianVelocity(point)
    CASE(3) ! ttauri stellar wind
       analyticalVelocity = TTauriStellarWindVelocity(point)
    CASE(4) ! whitney
       analyticalVelocity = WhitneyVelocity(point)
    CASE(5) ! whitney
       analyticalVelocity = ulrichVelocity(point)



    CASE DEFAULT
       CALL writeFatal("Error in logic in amrGridVelocity")
    END SELECT
  END FUNCTION analyticalVelocity

  TYPE(VECTOR) FUNCTION ulrichVelocity(point)
    USE utils_mod, ONLY : ccubsolv
    USE inputs_mod, ONLY : erinner, erouter
    USE inputs_mod, ONLY : mcore
    TYPE(VECTOR) :: point
    REAL(DOUBLE) :: vr, vt, vp, mu, r, mu_0
    REAL(DOUBLE) :: costheta, sintheta, sinphi, cosphi, sintheta0
    REAL(DOUBLE) :: vx, vy, vz, v
    COMPLEX(DOUBLE) :: a(3),z(3)
    r = modulus(point)*1.d10

    ulrichVelocity = VECTOR(0.d0, 0.d0, 0.d0)
    mu = (point%z*1.e10) /r
    IF ((r > erInner).AND.(r < erOuter)) THEN

       a(1) = CMPLX(-mu * (r/erinner),0.d0,kind=DOUBLE)
       a(2) = CMPLX((r/erinner-1.d0), 0.d0,kind=DOUBLE)
       a(3) = 0.d0
       CALL ccubsolv(a,z)
       mu_0 = REAL(z(1))
       IF ((imag(z(1)) /= 0.d0).OR.(imag(z(2))==0.d0).OR.(imag(z(3))==0.d0)) THEN
          WRITE(*,*)"problem with cubic solver"
       ENDIF
       sintheta0 = SQRT(1.d0 - mu_0**2)
       costheta = point%z*1.d10/r
       sintheta = SQRT(1.d0- costheta**2)
       cosphi = point%x/SQRT(point%x**2 + point%y**2)
       sinphi = point%y/SQRT(point%x**2 + point%y**2)

       v = SQRT(bigG * mcore/r)/cSpeed

       vr = -v * SQRT(1.d0 + costheta/mu_0)
       vt = v * (mu_0 - costheta) * SQRT(ABS((mu_0+costheta)/(mu_0*costheta)))
       vp = v * (sintheta0/sintheta) * SQRT(1.d0-costheta/mu_0)

       vx = vr * cosphi * sintheta + vt * cosphi * costheta - vp * sinphi
       vy = vr * sinphi * sintheta + vt * sinphi * costheta + vp * cosphi
       vz = vr * costheta - vt  * sintheta

       ulrichVelocity = VECTOR(vx, vy, vz)



    ENDIF

  END FUNCTION ulrichVelocity


  !      type(VECTOR) function TTauriStellarWindVelocity(point)
  !        use inputs_mod, only : ttauriMstar, ttaurirStar, SW_vMin, &
  !             SW_vMax, SW_Rmin, SW_beta, SW_rmax, SW_protation, SW_Veq
  !
  !        type(VECTOR), intent(in) :: point
  !        type(VECTOR) :: rHat, phiHat
  !        real(double) :: r, v_r, vEsc, vterm, v_eq, v_phi, x, y, z, phi
  !        real(double) :: theta
  !        type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.0d0)
  !
  !        r = modulus(point)
  !        x=point%x; y=point%y; z=point%z
  !
  !        if ((r > SW_rMin).and.(r < SW_rMax)) then
  !           vEsc = sqrt(2.d0 * bigG * ttauriMstar / ttaurirStar)
  !           vterm = SW_vMax * vEsc
  !           v_r = SW_vMin + (vTerm - SW_Vmin) * (1.d0 - SW_Rmin/r)**SW_beta
  !           if (SW_Protation > 0.d0) then
  !              v_eq = twoPi * ttaurirstar / sw_protation ! [cm/s]
  !           else
  !              v_eq = SW_Veq
  !           endif
  !           phi=ATAN2(y,x)
  !           theta=ACOS(z/r)
  !           if (phi < 0.d0) phi=phi+twoPi
  ! !          print *, "x", x, "y", y, "z", z, "r", r
  !           v_phi = v_eq * (SW_Rmin / r) * sin(theta)
  !
  !           rHat = point
  !           call normalize(rHat)
  !
  !           phiHat = point.cross.zAxis
  !           if (modulus(phiHat) /= 0.d0) then
  !              call normalize(phiHat)
  !              TTauriStellarWindVelocity = (v_r/cspeed) * rHat + (v_phi/cspeed) * phiHat
  !           else
  !              TTauriStellarWindVelocity = (v_r/cspeed) * rHat
  !           endif
  !
  !        else
  !           TTauriStellarWindVelocity = VECTOR(0.d0, 0.d0, 0.d0)
  !        endif
  !
  !      end function TTauriStellarWindVelocity




TYPE (VECTOR) FUNCTION TTauriStellarWindVelocity(point)
  USE inputs_mod, ONLY : ttauriRouter, dipoleOffset, &
       ttauriRstar, SW_openAngle, ttauriMstar, SW_vMin, &
       SW_vMax, SW_rMin, SW_beta, SW_rMax, SW_protation, SW_Veq, &
       SW_eqGap, SW_alfven

  TYPE(VECTOR), INTENT(IN) :: point
  TYPE(VECTOR) :: rVec, rVecDash, wind, phiHat

  REAL(DOUBLE) :: theta, beta
  REAL(DOUBLE) :: rBoundary, openAngleDash
  REAL(DOUBLE) :: r, rDash, thetaDash, phiDash
  REAL(DOUBLE) :: rMaxMax, sin2theta0dash
  REAL(DOUBLE) :: y, rAlfven, vAlfven
  REAL(DOUBLE) :: vMax, radV, Veq, vPhi

  rVec = point
  r = modulus(rVec)
  rVec = rVec * 1.d10
  wind = vector(0.d0, 0.d0, 0.d0)

  beta = dipoleOffset
  openAngleDash = SW_openAngle - beta
  rVecDash = rotateY(rVec, -beta)
  rDash = modulus(rVec)
  theta = ACOS(rVec%z/rDash)
  thetaDash = ACOS(rVecDash%z/rDash)
  IF (thetaDash == 0.d0) thetaDash = 1.d-20
  phiDash = ATAN2(rVecDash%y, rVecDash%x)
  IF (phiDash == 0.d0) phiDash = 1.d-20
  sin2theta0dash = (1.d0 + TAN(beta)**2.d0 * COS(phiDash)**2.d0)**(-1.d0)
  rMaxMax = SW_eqGap / sin2theta0dash

  rBoundary = rMaxMax * SIN(SW_openAngle)**2.d0
  rAlfven = SW_alfven * rBoundary * SIN(SW_openAngle) / 1.d10
  IF (rDash <= rBoundary) THEN
     y = SIN(thetaDash)**2.d0
     wind = vector(3.d0 * SQRT(y) * SQRT(1.d0-y) / SQRT(4.d0 - (3.d0*y)), &
          0.d0, (2.d0 - 3.d0 * y) / SQRT(4.d0 - 3.d0 * y))

     IF ((rVecDash%z/rDash) < 0.d0) wind%z = (-1.d0)*wind%z
     IF ((rVecDash%z/rDash) < 0.d0) wind = (-1.d0)*wind
     IF (rVecDash%z < 0.d0) wind = (-1.d0)*wind
     wind = rotateZ(wind, -phiDash)
     wind = rotateY(wind, beta)
  ELSE
     wind = point
  ENDIF
  CALL normalize(wind)
  !!max speed is multiple of escape velocity
  vMax = SW_vMax * SQRT(2.d0*bigG*ttauriMstar/ttauriRstar)
  !!radial velocity given by beta law
  radV = SW_vMin + (vMax - SW_vMin)*(1.d0 - SW_rMin/r)**SW_beta

  IF ((SW_Protation /= 0.d0).OR.(SW_Veq /= 0.d0)) THEN
     IF (SW_Protation > 0.d0) THEN
        Veq = twoPi * tTauriRstar / SW_Protation ![cm/s]
     ELSE
        Veq = SW_Veq
     END IF
     vAlfven = Veq * (rAlfven/SW_rMin)
     r = r*SIN(theta)
     IF (r <= rAlfven) THEN
       vPhi = Veq * (r/SW_rMin)
     ELSE
       vPhi = vAlfven * (rAlfven / r)
     END IF
     phiHat = point.cross.vector(0.d0,0.d0,1.d0)
     CALL normalize(phiHat)
     TTauriStellarWindVelocity = ((radV * wind) + (vPhi * phiHat)) / cspeed
  ELSE
     TTauriStellarWindVelocity = (radV / cSpeed) * wind
  END IF
END FUNCTION TTauriStellarWindVelocity


END MODULE analytical_velocity_mod
