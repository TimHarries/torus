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

contains

  logical function inFlowMahdavi(rVec) 
    use input_variables, only : ttauriRinner, ttauriRouter, dipoleOffset, ttauriRstar
    type(VECTOR) :: rVec
    real(double) :: r, theta, phi
    real(double) :: rDash, thetaDash, phiDash, beta
    real(double) :: thisRmax, sin2theta0dash,rmaxmin, rmaxmax
    beta = dipoleOffset
    r = modulus(rVec)    
    theta = acos(rVec%z/r)

    phi = atan2(rVec%y, rVec%x)
    if (phi < 0.d0) phi = phi + twoPi

    rDash = r
    phiDash = atan((sin(phi)*sin(theta))/(cos(phi)*sin(theta)*cos(beta)-cos(theta)*sin(beta)))
    if (phiDash /= 0.d0) then
       thetaDash = asin(max(-1.d0,min(1.d0,sin(phi)*sin(theta)/sin(phiDash))))
    else
       thetaDash = theta
    endif

    sin2theta0dash = (1.d0 + tan(beta)**2 * cos(phiDash)**2)**(-1.d0)

    if (sin(thetaDash) == 0.d0) then
       inFlowMahdavi = .false.
       goto 666
    endif
    thisrMax = rDash / sin(thetaDash)**2

    rMaxMin = ttauriRinner / sin2theta0dash
    rMaxMax = ttauriRouter / sin2theta0dash

    inflowMahdavi = .false.
    if ((thisRmax > rMaxMin).and.(thisRmax < rMaxMax)) inflowMahdavi = .true.
    if ((rVec%z < 0.d0).and.(rVec%x > 0.d0)) inflowMahdavi = .false.
    if ((rVec%z > 0.d0).and.(rVec%x < 0.d0)) inflowMahdavi = .false.
    if (r < ttauriRstar) inflowMahdavi = .false.
666 continue
  end function inFlowMahdavi

  type (VECTOR) function velocityMahdavi(point, grid)
    use input_variables, only : dipoleOffset, ttauriRInner, ttauriRouter, ttauriMstar, &
         ttaurirstar
    type(GRIDTYPE), intent(in) :: grid
    type(VECTOR), intent(in) :: point
    type(VECTOR) :: rvec, vp, magneticAxis, rVecDash
    real(double) :: r, rDash, phi, phiDash, theta,thetaDash,sin2theta0dash, beta
    real(double) :: deltaU, y, modVp, thisRmax, cosThetaDash, rTrunc, rMaxMin,rMaxMax

    rVec = point*1.d10

    velocityMahdavi = VECTOR(0.d0, 0.d0, 0.d0)
    if (modulus(rVec) < ttaurirstar) goto 666
    if (.not.inFlowMahdavi(rVec)) goto 666

    beta = dipoleOffset
    r = modulus(rVec)
    theta = acos(rVec%z/r)
    phi = atan2(rVec%y, rVec%x)
    if (phi < 0.d0) phi = phi + twoPi

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
    velocityMahdavi = vp

666 continue
  end function velocityMahdavi


  logical function inflowBlandfordPayne(rVec)
    use input_variables, only : DW_Rmin, DW_Rmax, DW_theta
    type(VECTOR) :: rVec
    real(double) :: r0, r, z, zMax, zMin
       
    inFlowBlandfordPayne = .false.
    r0 = sqrt(rVec%x**2 + rVec%y**2)
    if (r0 < DW_Rmin) goto 666
    r = modulus(rVec)
    z = abs(rvec%z)
    zMax = (r0-DW_Rmin) * tan(DW_theta)
    Zmin = 0.d0
    if (r0 > DW_Rmax) zMin = (r0-DW_Rmax)*tan(DW_theta)
    if ((z >= zMin).and.(z <= zMax)) then
       inFlowBlandfordPayne = .true.
    endif
666 continue
  end function inFlowBlandfordPayne
    

  type (VECTOR) function velocityBlandfordPayne(point, grid)
    use input_variables, only : DW_theta, DW_rMin, ttauriMstar, ttauriRstar
    type(GRIDTYPE), intent(in) :: grid
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
    
    vel = vEsc * (1.d0/sqrt(x))*sqrt(1.d0 - r0/r) ! Kwan Edwards & Fischer

    velocityBlandfordPayne = (vel/cSpeed) * rVec 

    velocityBlandfordPayne = rotateZ(velocityBlandfordPayne, -phi)


666 continue
  end function velocityBlandfordPayne

  function accretingAreaMahdavi(grid) result (accretingarea)
    use input_variables, only : ttauriRstar
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec
    real(double) :: accretingArea
    integer :: i ,j
    integer :: nLines
    nLines = 100000
    j = 0
    do i = 1, nLines
       rVec = (ttauriRstar + 1.d10 *grid%halfSmallestSubcell)*randomUnitVector() 
       if (inFlowMahdavi(rVec)) j = j + 1
    enddo
    accretingArea = fourPi * ttauriRstar**2 * dble(j)/dble(nLines)
  end function accretingAreaMahdavi

   function rhoBlandfordPayne(grid, rVec) result(rho)
     use input_variables, only : DW_rmax, DW_rmin, DW_mdot
     type(GRIDTYPE) :: grid
     type(VECTOR) :: rVec
     real(double) :: rho, kconst,vel,mdot

     mdot = DW_mdot*mSol/(365.25d0*24.d0*3600.d0)
     kconst = 0.5d0*(mdot/pi)/((DW_rMax**2-DW_rmin**2)*1.d20)
     rho = 1.d-25
     if (inflowBlandFordPayne(rVec)) then
        vel = modulus(velocityBlandfordPayne(rVec,grid))*cSpeed
        rho = kconst / vel
     endif

   end function rhoBlandfordPayne

end module magnetic_mod
