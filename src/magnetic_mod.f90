! Magnetic field stuff added in order to compute zeeman profiles
! TJH in December 2003
! This is for the AMR grid style only

module magnetic_mod

  use gridtype_mod
  use vector_mod
  use kind_mod
  use constants_mod
  implicit none

  public

  interface inflowMahdavi
     module procedure inflowMahdaviSingle
     module procedure inflowMahdaviArray
  end interface

  interface inflowBlandfordPayne
     module procedure inflowBlandfordPayneSingle
     module procedure inflowBlandfordPayneArray
  end interface

contains


  logical function inFlowMahdaviArray(vArray)
    type(VECTOR) :: vArray(:)
    integer :: i

    inFlowMahdaviArray = .true.
    do i = 1, SIZE(vArray)
       inFlowMahdaviArray = inflowMahdaviSingle(vArray(i))
       if (.not.inFlowMahdaviArray) exit
    enddo
  end function inFlowMahdaviArray

  real(double) function pMahdavi(rVec)
!    use inputs_mod, only : ttauriRinner, ttauriRouter, dipoleOffset, ttauriRstar, &
!         TTauriDiskHeight
    use inputs_mod, only : dipoleOffset
    type(VECTOR) :: rVec
    real(double) :: r, theta, phi
    real(double) :: rDash, phiDash, beta
    real(double) :: thisY, thisX
    beta = dipoleOffset
    r = modulus(rVec)

    thisY = rVec%y
    thisX = rVec%x
    if ((thisY == 0.d0).and.(rVec%x > 0.d0)) thisY = 1.d-10
    if ((thisY == 0.d0).and.(rVec%x < 0.d0)) thisY = -1.d-5
    if ((thisX == 0.d0).and.(rVec%y > 0.d0)) thisX = -1.d-10
    phi = atan2(thisY,thisX)
    rDash = r
    theta = acos(rVec%z/r)
    phiDash = atan2((sin(phi)*sin(theta)),(cos(phi)*sin(theta)*cos(beta)-cos(theta)*sin(beta)))

    pMahdavi = 0.5d0*(1.d0+tan(beta)**2 * cos(phiDash)**2)**(-1.d0)

    if (cos(phiDash) > 0.d0) pMahdavi = 1.d0 - pMahdavi
    pMahdavi = 0.5d0

  end function pMahdavi


  !!Functon determines if there is stellar outflow from poles (ttauri) at the position given by rVec
  !!tjgw201 14/03/19 based on function: inFlowMahdaviSingle(rVec)
  logical function outFlow(rVec)
    use inputs_mod, only : ttauriRouter, dipoleOffset, ttauriRstar
    type(VECTOR) :: rVec
    real(double) :: r, theta, phi, beta, xi
    real(double) :: thisr, thisTheta
    real(double) :: thetaDash, phiDash, thisThetaFinal
    real(double) :: thetaFinal

    beta = dipoleOffset
    r = modulus(rVec)
    xi = ttauriRstar / ttauriRouter
    if (r == 0.d0) then
       outFlow = .false.
       goto 666
    endif

    thisTheta = acos(rVec%z/r)
    phi = atan2(rVec%y,rVec%x)+1.d-3

    !theta_0 = pi/2 therefore can remove all theta terms
    phiDash = atan2(sin(phi),(cos(phi)*cos(beta)))
    if (phiDash == 0.d0) phiDash = 1.d-20

    thetaDash = acos(cos(phi)*sin(beta))
    thetaFinal = asin(sqrt(xi*sin(thetaDash)**2.))

    thisThetaFinal = asin(max(-1.d0,min(1.d0,sin(phi)*sin(thisTheta)/sin(phidash))))
    if(thisThetaFinal .LE. thetaFinal) then
      outflow = .True.
    else
      outflow = .False.
    endif
    !
    ! if (beta /= 0.d0) then
    !    if (cos(theta) > 0.d0) then
    !       if (cos(phidash) < 0.d0) outFlow = .false.
    !    else
    !       if (cos(phidash) > 0.d0) outFlow = .false.
    !    endif
    ! endif

    666 continue
  end function outFlow


!!returns the half opening angle between accretion hot spots - tjgw201 15/03/19
  function capHalfAngle()
    use inputs_mod, only : ttauriRouter, dipoleOffset, ttauriRstar
    real(double) :: theta0dash, xi
    real(double) :: capHalfAngle

    xi = ttauriRstar / ttauriRouter

    theta0dash = acos(sin(dipoleOffset))
    capHalfAngle = asin(sqrt(xi*sin(theta0Dash)**2.))
    return
  end function capHalfAngle


  logical function inFlowMahdaviSingle(rVec)
    use inputs_mod, only : ttauriRinner, ttauriRouter, dipoleOffset, ttauriRstar, &
         TTauriDiskHeight
    type(VECTOR) :: rVec
    real(double) :: r, theta, phi
    real(double) :: rDash, thetaDash, phiDash, beta
    real(double) :: thisRmax, sin2theta0dash,rmaxmin, rmaxmax,thisY, thisX

    beta = dipoleOffset
    r = modulus(rVec)
    if (r == 0.d0) then
       inflowMahdaviSingle = .false.
       goto 666
    endif
    theta = acos(rVec%z/r)

    thisY = rVec%y
    thisX = rVec%x
!    if ((thisY == 0.d0).and.(rVec%x > 0.d0)) thisY = 1.d-10
!    if ((thisY == 0.d0).and.(rVec%x < 0.d0)) thisY = -1.d-5
!    if ((thisX == 0.d0).and.(rVec%y > 0.d0)) thisX = -1.d-10
    phi = atan2(thisY,thisX)+1.d-3

    rDash = r
    phiDash = atan2((sin(phi)*sin(theta)),(cos(phi)*sin(theta)*cos(beta)-cos(theta)*sin(beta)))
    if (phiDash == 0.d0) phiDash = 1.d-20
    thetaDash = asin(max(-1.d0,min(1.d0,sin(phi)*sin(theta)/sin(phiDash))))


    sin2theta0dash = (1.d0 + tan(beta)**2 * cos(phiDash)**2)**(-1.d0)

    if (sin(thetaDash) == 0.d0) then
       inFlowMahdaviSingle = .false.
       goto 666
    endif
    thisrMax = rDash / sin(thetaDash)**2

    rMaxMin = ttauriRinner / sin2theta0dash
    rMaxMax = ttauriRouter / sin2theta0dash

    inflowMahdaviSingle = .false.
    if ((thisRmax >= rMaxMin).and.(thisRmax <= rMaxMax)) inflowMahdaviSingle = .true.
    if (r < ttauriRstar) inflowMahdaviSingle = .false.
    if (abs(rVec%z) < TTauriDiskHeight) inflowMahdaviSingle = .false.

!    pMahdavi = 0.5d0*(1.d0+tan(beta)**2 * cos(phiDash)**2)**(-1.d0)
!    if (pMahdavi < 1.d0) inflowMahdaviSingle = .false.


    if (beta /= 0.d0) then
       if (cos(theta) > 0.d0) then
          if (cos(phidash) < 0.d0) inFlowMahdaviSingle = .false.
       else
          if (cos(phidash) > 0.d0) inFlowMahdaviSingle = .false.
       endif
    endif


666 continue
  end function inFlowMahdaviSingle



  type (VECTOR) function velocityMahdavi(point)
    use inputs_mod, only : dipoleOffset, ttauriRInner, ttauriRouter, ttauriMstar, &
         ttaurirstar
    type(VECTOR), intent(in) :: point
    type(VECTOR) :: rvec, vp, magneticAxis, rVecDash, vSolid
    real(double) :: r, rDash, phi, phiDash, theta,thetaDash,sin2theta0dash, beta
    real(double) :: deltaU, y, modVp, thisRmax, cosThetaDash, rTrunc, rMaxMin,rMaxMax
    real(double) :: velMagAtCorotation


    rVec = point*1.d10

    velocityMahdavi = VECTOR(0.d0, 0.d0, 0.d0)
    if (modulus(rVec) < ttaurirstar) goto 666
    if (.not.inFlowMahdavi(rVec)) goto 666

    beta = dipoleOffset
    r = modulus(rVec)
    theta = acos(rVec%z/r)
    phi = atan2(rVec%y, rVec%x)
!    if (phi < 0.d0) phi = phi + twoPi

    rVecDash = rotateY(rVec, -beta)
    rDash = modulus(rVecDash)
    cosThetaDash = rVecDash%z / rDash
    thetaDash = acos(rVecDash%z/rDash)
    phiDash = atan2(rVecDash%y, rVecDash%x)

    sin2theta0dash = (1.d0 + tan(beta)**2 * cos(phiDash)**2)**(-1.d0)

    thisrMax = rDash / sin(thetaDash)**2
    rMaxMin = ttauriRinner / sin2theta0dash
    rMaxMax = ttauriRouter / sin2theta0dash
    rTrunc = ttauriRinner + (ttauriRouter-ttauriRinner)*(thisrMax-rMaxmin)/(rMaxMax-rMaxMin)

    deltaU =  bigG * TTauriMstar * (1.d0/r - 1.d0/rTrunc)
    modvp = sqrt(2.d0*abs(deltaU))
    if (deltaU < 0.d0) modvp = -modvp
    y = SIN(thetaDash)**2

    vP = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
         0.0, &
         (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))
    vP = (-modVp/cSpeed) * vP
    IF (costhetaDash < 0.d0) vP%z = -vp%z
    IF (costhetaDash < 0.d0) vP = (-1.d0)*vp

    if (rVec%z < 0.d0) vP = (-1.d0)*vp
    magneticAxis = VECTOR(0.d0, sin(beta), cos(beta))
    vp = rotateZ(vp, -phiDash)

    vp = rotateY(vp, beta)

    velMagAtCorotation = sqrt(bigG*ttauriMstar/ttauriRouter)/cSpeed
    rVec = VECTOR(point%x,point%y,0.d0)
    vSolid = rVec .cross. VECTOR(0.d0, 0.d0, 1.d0)
    call normalize(vSolid)
    vSolid = (modulus(rVec)/ttauriRouter*velMagAtCorotation) * vSolid


    velocityMahdavi = vp + vSolid


666 continue

  end function velocityMahdavi


!Equation give in Hartmann et al. (1994) - tjgw201 10/04/19
  real (double) function densityMahdavi(point, mdot)
    use inputs_mod, only : dipoleOffset, ttauriRInner, ttauriRouter, ttauriMstar, &
         ttaurirstar
    type(VECTOR), intent(in) :: point
    type(VECTOR) :: rVec
    real(double) :: r, theta, mdot
    real(double) :: y, rho

    rVec = point * 1.d10
    r = modulus(rVec)
    theta = acos(rVec%z/r)
    y = sin(theta)**2.d0

    rho = mdot*sqrt(r**(-5.d0))*sqrt((4.d0 - 3.d0*y)/ (1.d0-y)) / sqrt(2.d0*bigG*ttauriMstar)
    densityMahdavi = rho / (fourpi*((1.d0/ttauriRInner) - (1.d0/ttauriRouter)))

  end function densityMahdavi


  type (VECTOR) function velocityAlphadisc(point)
    type(VECTOR), intent(in) :: point

    velocityAlphaDisc = point
    velocityAlphaDisc = ttauriKeplerianVelocity(point)
  end function velocityAlphaDisc

  subroutine safierFits(solution, chi, zeta0dash, h0, zeta, psi, eta)
    character(len=1) :: solution
    real(double) :: eta, chi, zeta, psi ! safier 1993 ApJ 408 115
    real(double) :: a(8), b(0:9)
    real(double) :: psi1, psi2, zeta1, zeta2, psi0
    real(double) :: zeta1dash, zeta2dash, zeta0dash, zetadash
    real(double) :: h0
    select case (solution)
       case("C")
          a(1) = 0.21; a(2) = 0.30; a(3) = 1.23; a(4) = 0.21
          a(5) = 1.27; a(6) = 0.92; a(7) = 0.04; a(8) = 0.28
          b(0) = 0.035; b(1) = 12.16; b(2) = 0.65; b(3) = 0.33
          b(4) = 1.0; b(5) = 0.40; b(6) = 0.50; b(7) = 0.0
          b(8) = 1.25; b(9) = 0.6
       case DEFAULT
          write(*,*) "Solution ",solution, " not found in safierFits"
          stop
    end select
    chi = max(1.d-20,chi)
    psi0 = 0.03
    zeta1 = (1.d0 + a(1)*(chi-h0) + a(2)*(chi-h0)**a(3))*exp(-a(4)*(chi-h0))
    psi1 = (b(0) + b(6) * (chi - h0) + b(7)*(chi - h0)**2)*exp(-b(8)*(chi-h0)**b(9))
    zeta2 = a(5) * (chi - h0)**a(6) * exp(-4.d0*a(7)*(chi-h0)**a(8))
    psi2 = b(1)*exp(-1.d0/max(1.d-20,(b(2)*(chi-h0)**b(3)))) * exp(-b(4)/max(1.d-20,((chi-h0)*b(5))))

    zeta1dash = -a(4)*exp(-a(4)*(chi-h0)) * (1.d0 + a(1)*(chi-h0)+a(2)*(chi-h0)**a(3)) + &
         (a(1)+a(2)*a(3)*(chi-h0)*(a(3)-1.d0))*exp(-a(4)*(chi-h0))

    zeta2dash = a(5)*a(6)*(chi-h0)**(a(6)-1.d0) * exp(-4.d0*a(7)*(chi-h0)**a(8)) + &
         a(5)*(chi-h0)*a(6)*(-4.d0*a(8)*a(7)*(chi-h0)**(a(8)-1.d0)*exp(-4.d0 * a(7) * (chi-h0)**a(8)))
    zetadash = zeta1dash + zeta2dash

    zeta = zeta1 + zeta2
    psi = psi1 + psi2

    eta = psi0 * (1.d0 - h0*zeta0dash)/(zeta * psi * (zeta - chi * zetaDash))
  end subroutine safierFits

  type (VECTOR) function TTauriKeplerianVelocity(point)
    use inputs_mod, only : sourceMass
    type(VECTOR), intent(in) :: point
    type(VECTOR) :: vVec, rVec
    real(double) :: v
    type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.d0)

    rVec = VECTOR(point%x, point%y, 0.d0)
    vVec = rVec .cross. zAxis
    call normalize(vVec)
    v = 0.d0
    if (modulus(rVec) /= 0.d0) then
       v = sqrt(bigG*sourceMass(1)/(modulus(rVec)*1.d10))/cSpeed
    endif
    TTauriKeplerianVelocity = v * vVec
  end function TTauriKeplerianVelocity

  logical function inflowBlandfordPayneArray(rVec)
    type(VECTOR) :: rVec(:)
    integer :: i

    inflowBlandfordPayneArray = .true.
    do i = 1, size(rVec)
       inflowBlandfordPayneArray = inflowBlandfordPayneSingle(rVec(i))
       if (.not.inflowBlandfordPayneArray) exit
    enddo
  end function inflowBlandfordPayneArray

  logical function inflowBlandfordPayneSingle(rVec)
    use inputs_mod, only : DW_Rmin, DW_Rmax, DW_theta
    type(VECTOR) :: rVec
    real(double) :: r0, r, z, zMax, zMin

    inFlowBlandfordPayneSingle = .false.
    r0 = sqrt(rVec%x**2 + rVec%y**2)
    if (r0 < DW_Rmin) goto 666
    r = modulus(rVec)
    z = abs(rvec%z)
    zMax = (r0-DW_Rmin) * tan(DW_theta)
    Zmin = 0.d0
    if (r0 > DW_Rmax) zMin = (r0-DW_Rmax)*tan(DW_theta)
    if ((z >= zMin).and.(z <= zMax)) then
       inFlowBlandfordPayneSingle = .true.
    endif
666 continue
  end function inflowBlandfordPayneSingle


  type (VECTOR) function velocityBlandfordPayne(point)
    use inputs_mod, only : DW_theta, DW_rMin, ttauriMstar, ttauriRstar
    type(VECTOR), intent(in) :: point
    type(VECTOR) :: rvec
    real(double) :: phi, r0, r, Vesc, vel, x



    velocityBlandfordPayne = VECTOR(0.d0, 0.d0, 0.d0)
    if (.not.inflowBlandfordPayne(point)) goto 666

    phi = atan2(point%y, point%x)
    rVec = VECTOR(cos(DW_theta), 0.d0, sin(DW_theta))
    rVec = rotateZ(rVec, phi)
    r0 = sqrt(point%x**2 + point%y**2)
    r = modulus(point)
    if (point%z < 0.d0) rVec%z = -rVec%z

    Vesc = sqrt(2.d0*bigG * ttauriMstar/ttauriRstar)

    x = r0/DW_rMin

    if ((r /= 0.d0).and.(x >= 0.d0)) then
       if (r0/r < 1.d0) then
          vel = vEsc * (1.d0/sqrt(x))*sqrt(1.d0 - r0/r) ! Kwan Edwards & Fischer
!          vel = max(20.d5, vel)
          velocityBlandfordPayne = (vel/cSpeed) * rVec
          velocityBlandfordPayne = rotateZ(velocityBlandfordPayne, -phi)
          velocityBlandfordPayne = velocityBlandfordPayne + ttauriKeplerianVelocity(point)
       endif
    endif

666 continue
  end function velocityBlandfordPayne
!!!finds the accretion area by randomyl sampling vectors in the source - tjgw201 26/02/19
  function accretingAreaMahdavi() result (accretingarea)
    use inputs_mod, only : ttauriRstar
    type(VECTOR) :: rVec
    real(double) :: accretingArea
    integer :: i ,j
    integer, parameter :: nLines = 1000000
    j = 0
    do i = 1, nLines
       rVec = (ttauriRstar * 1.01d0)*randomUnitVector()
       if (inFlowMahdavi(rVec)) j = j + 1
    enddo
    accretingArea = fourPi * ttauriRstar**2 * dble(j)/dble(nLines)
  end function accretingAreaMahdavi

   function rhoBlandfordPayne(rVec) result(rho)
     use inputs_mod, only : DW_rmax, DW_rmin, DW_mdot
     type(VECTOR) :: rVec
     real(double) :: rho, kconst,vel,mdot

     mdot = DW_mdot*mSol/(365.25d0*24.d0*3600.d0)
     kconst = 0.5d0*(mdot/pi)/((DW_rMax**2-DW_rmin**2)*1.d20)
     rho = 1.d-25
     if (inflowBlandFordPayne(rVec)) then
        vel = modulus(velocityBlandfordPayne(rVec))*cSpeed
        rho = kconst / vel
     endif

   end function rhoBlandfordPayne

   function rhoAlphaDisc(grid, rVec) result(rho)
     use inputs_mod, only : rinner, router, mdisc, alphaDisc, betaDisc, height
     type(GRIDTYPE) :: grid
     type(VECTOR) :: rVec
     real(double) :: rho, r, fac, r0
     real(double) :: h, rho0
    logical :: test

    test=grid%octreeRoot%threed

     r0 = 100.d0 * autocm/1.d10
     fac = betaDisc - alphaDisc + 2.d0
     rho0 = mDisc / (twoPi**1.5 * height * r0**(-betaDisc) * rInner**alphaDisc * &
          (rOuter**fac - rInner**fac))
     rho0 = rho0 / 1.d30
     r = sqrt(rvec%x**2 + rvec%y**2)
      rho = 1.d-30
!      write(*,*) "r ",r*1d10/rsol,rinner*1.d10/rsol,router*1.d10/rsol
      if ( (r > Rinner).and.(r < rOuter) ) then
         r = sqrt(rVec%x**2 + rVec%y**2)
         h = height * (r / (100.d0*autocm/1.d10))**betaDisc
         fac = -0.5d0 * (dble(rVec%z)/h)**2
         rho = dble(rho0) * (dble(rInner/r))**dble(alphaDisc) * exp(fac)
      endif


   end function rhoAlphaDisc


end module magnetic_mod
