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
    thetaDash = asin(max(-1.d0,min(1.d0,sin(phi)*sin(theta)/sin(phiDash))))

    sin2theta0dash = (1.d0 + tan(beta)**2 * cos(phiDash)**2)**(-1.d0)
    
    thisrMax = rDash / sin(thetaDash)**2

    rMaxMin = ttauriRinner / sin2theta0dash
    rMaxMax = ttauriRouter / sin2theta0dash

    inflowMahdavi = .false.
    if ((thisRmax > rMaxMin).and.(thisRmax < rMaxMax)) inflowMahdavi = .true.
    if ((rVec%z < 0.d0).and.(rVec%x > 0.d0)) inflowMahdavi = .false.
    if ((rVec%z > 0.d0).and.(rVec%x < 0.d0)) inflowMahdavi = .false.
    if (r < ttauriRstar) inflowMahdavi = .false.
  end function inFlowMahdavi

  type (VECTOR) function velocityMahdavi(point, grid)
    use input_variables, only : dipoleOffset, mCore, ttauriRInner, ttauriRouter, ttauriMstar
    type(GRIDTYPE) :: grid
    type(VECTOR) :: point, rvec, vp, magneticAxis, rVecDash
    real(double) :: r, rDash, phi, phiDash, theta,thetaDash,sin2theta0dash, beta
    real(double) :: deltaU, v, y, modVp, thisRmax, cosThetaDash, rTrunc, rMaxMin,rMaxMax

    rVec = point*1.d10

    velocityMahdavi = VECTOR(0.d0, 0.d0, 0.d0)
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

  subroutine assignDensitiesMahdavi(grid, mdot)
    use amr_mod, only : findSubcellTD, findsubcellLocal, inOctal, cellVolume
    use input_variables, only : ttauriRstar, dipoleOffset
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    type(VECTOR) :: rVec, rVecDash, thisRvec
    real(double) :: mdot
    integer :: nLines, iline
    type(OCTAL), pointer :: thisOctal
    integer :: subcell, nr
    real(double) :: accretingArea, vstar, astar, beta
    real(double) :: thetaDash, phiDash, rDash, r, cosThetaDash
    real(double) :: thisrMax, sin2theta0dash, thisr,thisrho, thisPhi, thisTheta
    beta = dipoleOffset
    nLines = 1000000
    j = 0
    do i = 1, nLines
       rVec = ttauriRstar * (1.01d0)*randomUnitVector() 
       if (inFlowMahdavi(rVec)) j = j + 1
    enddo
    accretingArea = fourPi * ttauriRstar**2 * dble(j)/dble(nLines)
    write(*,*) "Fractional area of accretion (%): ",100.d0 *  dble(j)/dble(nLines)
    nLines = 10000
    aStar = accretingArea / dble(nLines)
    do iLine = 1, nLines
       rVec = ttauriRstar * (1.01d0)*randomUnitVector() 
       do while(.not.inFlowMahdavi(rVec))
          rVec = ttauriRstar * (1.01d0)*randomUnitVector() 
       enddo
       call findSubcellTD(rVec/1.d10, grid%octreeRoot, thisOctal, subcell)
       thisOctal%velocity(subcell) = velocityMahdavi(rVec/1.d10, grid)
       vStar = modulus(thisOctal%velocity(subcell)*cspeed)
       rVecDash = rotateY(rVec, -beta)
       rDash = modulus(rVecDash)
       cosThetaDash = rVecDash%z / rDash
       thetaDash = acos(rVecDash%z/rDash)
       phiDash = atan2(rVecDash%y, rVecDash%x)
       sin2theta0dash = (1.d0 + tan(beta)**2 * cos(phiDash)**2)**(-1.d0)
       thisrMax = rDash / sin(thetaDash)**2

       nr = 1000
       thisOctal => grid%octreeRoot
       do i = 1, nr
          thisr = rDash + (thisRmax-rDash) * dble(i-1)/dble(nr-1)
          thisRho = (ttauriRstar/thisR)**3 * mdot /(aStar * vStar)
          thisphi = phiDash
          thisTheta = asin(sqrt(thisr/thisRmax))
          if (cosThetaDash  < 0.d0) thisTheta = pi - thisTheta
          thisRVec = VECTOR(thisR*cos(thisPhi)*sin(thisTheta),thisR*sin(thisPhi)*sin(thisTheta), &
               thisR * cos(thisTheta))
          thisRvec = rotateY(thisRvec, beta)
          if (inflowMahdavi(thisRVec)) then
             if (inOctal(grid%octreeRoot, thisRvec/1.d10)) then
                call findSubcellLocal(thisRvec/1.d10, thisOctal, subcell)
                thisOctal%rho(Subcell) = thisRho*(thisOctal%subcellSize*1.d10*aStar) / &
                     (cellVolume(thisOctal,subcell)*1.d30)
             endif
          endif
       enddo
    end do
  end subroutine assignDensitiesMahdavi

          
       
    

end module magnetic_mod
