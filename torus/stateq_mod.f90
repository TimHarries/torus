!
! a module to do 3d statistical equilibrium according to Klein & Castor (1978)
! - lots of it based on my "stateq" PhD code
!
module stateq_mod

  use gridtype_mod
  use grid_mod
  use math_mod
  use vector_mod
  use constants_mod
  use path_integral
  use jets_mod
  use opacity_lte_mod

  implicit none
  public

!  public initGridStateq, expint, gii, giia, beta_mn, cijt

  real(kind=doubleKind) :: bEinstein(9, 9) 
  
  real(kind=doubleKind), parameter :: eTrans(15) =                      &
       (/  0.000e0, 10.199e0, 12.088e0, 12.749e0, 13.055e0, 13.221e0,   &
       13.321e0, 13.386e0, 13.431e0, 13.463e0, 13.486e0, 13.504e0,      &
       13.518e0, 13.529e0, 13.538e0 /)

  real(kind=doubleKind), parameter :: gDegen(15) =                      &
       (/ 2.0, 8.0, 18.0, 32.0, 50.0, 72.0, 98.0, 128.0, 162.0,&
       200.0, 242.0, 288.0, 338.0, 392.0, 450.0 /)

  real(kind=doubleKind), parameter :: aEinstein(9, 9) = reshape( source=&
       (/ 0.000e0, 4.699e8, 5.575e7, 1.278e7, 4.125e6, 1.644e6, &
       7.568e5, 3.869e5, 2.143e5,  &
       0.000e0, 0.000e0, 4.410e7, 8.419e6, 2.530e6, 9.732e5, &
       4.389e5, 2.215e5, 1.216e5, &
       0.000e0, 0.000e0, 0.000e0, 8.986e6, 2.201e6, 7.783e5, &
       3.358e5, 1.651e5, 8.905e4, &
       0.000e0, 0.000e0, 0.000e0, 0.000e0, 2.699e6, 7.711e5, &
       3.041e5, 1.424e5, 7.459e4, &
       0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, 1.025e6, &
       3.253e5, 1.388e5, 6.908e4, &
       0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, &
       4.561e5, 1.561e5, 7.065e4, &
       0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, &
       0.000e0, 2.272e5, 8.237e4, &
       0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, &
       0.000e0, 0.000e0, 1.233e5,  &
       0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, 0.000e0, &
       0.000e0, 0.000e0, 0.000e0 /), shape=(/9,9/))
  
  real(kind=doubleKind) :: fStrength(9,9) = reshape( source=& 
       (/ 0.000e00, 0.000e00, 0.000e00, 0.000e00, 0.000e00, &
       0.000e00, 0.000e00, 0.000e00, 0.000e00, &
       4.162e-1, 0.000e00, 0.000e00, 0.000e00, 0.000e00, &
       0.000e00, 0.000e00, 0.000e00, 0.000e00, &
       7.910e-2, 6.407e-1, 0.000e00, 0.000e00, 0.000e00, &
       0.000e00, 0.000e00, 0.000e00, 0.000e00, &
       2.899e-2, 1.193e-1, 8.421e-1, 0.000e00, 0.000e00, &
       0.000e00, 0.000e00, 0.000e00, 0.000e00, &
       1.394e-2, 4.467e-2, 1.506e-1, 1.038e00, 0.000e00, &
       0.000e00, 0.000e00, 0.000e00, 0.000e00, &
       7.799e-3, 2.209e-2, 5.584e-2, 1.793e-1, 1.231e00, &
       0.000e00, 0.000e00, 0.000e00, 0.000e00, &
       4.814e-3, 1.270e-2, 2.768e-2, 6.549e-2, 2.069e-1, &
       1.424e00, 0.000e00, 0.000e00, 0.000e00, &
       3.183e-3, 8.036e-3, 1.604e-2, 3.230e-2, 7.448e-2, &
       2.340e-1, 1.616e00, 0.000e00, 0.000e00, &
       2.216e-3, 5.429e-3, 1.023e-2, 1.870e-2, 3.645e-2, &
       8.315e-2, 2.609e-1, 1.807e00, 0.000e00 /), shape=(/9,9/))
       
  real(kind=doubleKind) :: lambdaTrans(9,9) = reshape( source=&
       (/000000.e-8,1215.67e-8,1025.72e-8,992.537e-8,949.743e-8, &
       937.803e-8,930.748e-8,926.226e-8,923.150e-8, &
       0.00000e-8,0000000e-8,6562.80e-8,4861.32e-8,4340.36e-8, &
       4101.73e-8,3970.07e-8,3889.05e-8,3835.38e-8, &
       0.00000e-8,0.00000e-8,0000000e-8,18751.0e-8,12818.1e-8, &
       10938.1e-8,10049.4e-8,9545.98e-8,9229.02e-8, &
       0.00000e-8,0.00000e-8,0.00000e-8,0000000e-8,40512.0e-8, &
       26252.0e-8,21655.0e-8,19445.6e-8,18174.1e-8, &
       0.00000e-8,0.00000e-8,0.00000e-8,0.00000e-8,0000000e-8, &
       74578.0e-8,46525.0e-8,37395.0e-8,32961.0e-8, &
       0.00000e-8,0.00000e-8,0.00000e-8,0.00000e-8,0.00000e-8, &
       0.00000e-8,123680.e-8,75005.0e-8,59066.0e-8, &
       0.00000e-8,0000000e-8,0000000e-8,0000000e-8,0000000e-8, &
       0.00000e-8,0000000e-8,190570.e-8,113060.e-8, &
       0.00000e-8,0000000e-8,0000000e-8,0000000e-8,0000000e-8, &
       0.00000e-8,0000000e-8,0000000e-8,277960.e-8, &
       0.00000e-8,0000000e-8,0000000e-8,0000000e-8,0000000e-8, &
       0.00000e-8,0000000e-8,0000000e-8,0.00000e-8 /), shape=(/9,9/))
  
contains

  real function beta_mn(m, n, rVec, i1, i2, i3, grid, thisOctal, thisSubcell)

    type(GRIDTYPE), intent(in)   :: grid
    type(VECTOR), intent(in)     :: rVec
    integer, intent(in)          :: i1, i2, i3
    integer, intent(in)          :: m,n
    type(OCTAL), pointer,optional:: thisOctal
    integer, intent(in),optional :: thisSubcell 
    
    type(VECTOR)      :: direction
    integer           :: i, j
    real              :: theta, phi
    integer,parameter :: nTheta = 10
    integer           :: nPhi
    real              :: dTheta, dPhi, dOmega
    real              :: totomega
    real              :: escprob,  tau_mn
    type(octalVector) :: rVecOctal
    type(OCTAL),pointer :: octalCopy
    
!integer :: labels(8)
!integer :: label1
!label1 = thisOctal%label(1)                
!labels = thisOctal%label
    dtheta = pi / real(ntheta-1)
    dphi = twopi / real(nphi-1)
    beta_mn = 0.
    totomega = 0.
    do i = 1, ntheta
       theta = pi*real(i-1)/real(ntheta-1)
       nphi = max(2,nint(real(ntheta)*sin(theta)))
       dphi = twopi / real(nphi-1)
       do j = 1, nphi-1
          phi = twopi * real(j-1)/real(nphi-1)
          direction = vector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          domega = sin(theta)*dtheta*dphi
          totomega = totomega + domega

!          call integratePath(real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                   grid%velocity(i1,i2,i3), &
!                   rVec, direction, grid, &
!                   lambdaArray, tauExt, tauAbs, tauSca, maxTau, nTau, .false., &
!                   escProb, .false., real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                   1, contTau, hitCore, .false.,.false., 0.,&
!                   VECTOR(0.,0.,0.), 1., .true., m, n, real(fstrength(m,n)), real(gDegen(m)), real(gDegen(n)), &
!                    tau_mn)

          tau_mn = (pi*echarge**2)/(melectron*cspeed)
          tau_mn = tau_mn * gDegen(m) * fStrength(m,n)
          tau_mn = tau_mn *  lambdaTrans(m,n) / cSpeed
          if (grid%adaptive) then
             if (grid%geometry(1:4) == "jets")  then
                tau_mn = tau_mn * abs((thisOctal%N(thisSubcell,m)/gDegen(m)) - &
                     (thisOctal%N(thisSubcell,n)/gDegen(n)))
                tau_mn = tau_mn/dV_dn_jets(rVec, direction)
             else
                tau_mn = tau_mn * abs((thisOctal%N(thisSubcell,m)/gDegen(m)) - &
                                      (thisOctal%N(thisSubcell,n)/gDegen(n))) ! eq 5.
                rVecOctal = rVec
                octalCopy => thisOctal
                tau_mn = tau_mn / (amrGridDirectionalDeriv(grid,rVecOctal,direction, &
                                                        startOctal=thisOctal) / 1.e10)
!if (label1 /= thisOctal%label(1)) print *, "Label change during amrGridDirectionalDeriv" 		
             end if
          else
             tau_mn = tau_mn * abs((grid%n(i1,i2,i3,m)/gDegen(m)) - &
                                   (grid%n(i1,i2,i3,n)/gDegen(n))) ! eq 5.
             tau_mn = tau_mn / (directionalderiv(grid,rvec,i1,i2,i3,direction)/1.e10)
          end if
          if (tau_mn < 0.) then
             tau_mn = 1.e-10 
          endif

        if (tau_mn < 0.1) then
          escProb = 1.0-tau_mn*0.5*(1.0 - tau_mn/3.0*(1. - tau_mn*0.25*(1.0 - 0.20*tau_mn)))
        else if (tau_mn < 15.) then
          escProb = (1.0-exp(-tau_mn))/tau_mn
        else
          escProb = 1./tau_mn
        end if
        beta_mn = beta_mn + escprob * domega

       enddo
    enddo
   beta_mn = beta_mn / totOmega 

!   if ((m == 1).and.(n == 2)) beta_mn = beta_mn * 100.
!    write(*,*) totOmega/fourPi
!   write(*,*) m,n,beta_mn

!if (SUM(labels-thisOctal%label) /= 0) print *, "            leaving beta_mn, labels = ", labels    
  end function beta_mn


  real function beta_cmn(m,n,rVec,i1,i2,i3,grid,nstar,thisOctal,thisSubcell)

    integer, intent(in)         :: m,n
    type(gridtype), intent(in)  :: grid
    type(vector), intent(in)    :: rVec
    integer, intent(in)         :: i1, i2, i3
    integer, parameter          :: ntheta=10
    integer, parameter          :: nphi=10
    integer, intent(in)         :: nstar
    type(OCTAL),pointer,optional:: thisOctal
    integer,intent(in),optional :: thisSubcell 
    
    real :: r, thetatostar, phiTostar
    real :: h, totomega, disttooccult
    integer :: i, j
    type(vector) :: tostar, occultposition, starposition
    type(vector) :: tooccult, direction
    real :: starradius, occultradius
    real :: disttostar
    real :: cosang, sinang, ang
    logical :: occulted, full
    real :: tau_cmn
    real :: dtheta, dphi
    real :: theta, phi, domega
    real :: escprob, dotprod
    integer, parameter :: maxTau = 4000
    integer :: nTau
    real, allocatable :: contTau(:,:)
    real :: tauExt(maxTau)
    real :: tauAbs(maxTau)
    real :: tauSca(maxTau)
    real :: lambdaArray(maxTau)
    logical :: hitCore
    real :: tauLocal
    type(octalVector) :: rVecOctal
    type(OCTAL),pointer :: octalCopy

!integer :: labels(8)
!integer :: label1
!labels = thisOctal%label
    full = .false.

    if (nstar == 1) then
       starradius = grid%rstar1
       occultradius = grid%rstar2
       starposition = grid%starpos1
       occultposition = grid%starpos2
    else
       starradius = grid%rstar2
       occultradius = grid%rstar1
       starposition = grid%starpos2
       occultposition = grid%starpos1
    endif

    tostar = starposition - rVec
    disttostar = modulus(tostar)
    tostar = tostar / disttostar

    tooccult = occultposition - rVec
    disttooccult = modulus(tooccult)
    tooccult = tooccult / disttooccult

    h  = sqrt(max(0.,(disttostar**2 - starradius**2)))
    cosang = h / disttostar
    ang = acos(min(1.,max(-1.,cosAng)))

    call getPolar(toStar, r, thetaTostar, phiToStar)

    dtheta = pi / real(ntheta-1)
    dphi = twopi / real(nphi-1)
    beta_cmn = 0.
    totomega = 0.
    do i = 1, ntheta
       theta = thetaToStar + (2.*real(i-1)/real(nTheta-1)-1.)*ang
       if (theta > pi) theta = theta - pi
       if (theta < 0.) theta = theta + pi
       dTheta = 2.*ang/real(nTheta-1)
       do j = 1, nphi
          phi = phiToStar + (2.*real(j-1)/real(nPhi-1)-1.)*ang
          if (phi < 0.) phi = phi + twoPi
          if (phi > twoPi) phi = phi - twoPi
          dphi = 2.*ang/real(nPhi-1)
          occulted = .false.
          direction = vector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          domega = sin(theta)*dtheta*dphi
          dotprod = direction .dot. tostar
          escProb = 0.
          if ((dotprod > 0.) .and. (acos(min(1.,max(-1.,dotprod))) < ang)) then

             dotprod = direction .dot. tooccult
             if ((dotprod > 0.).and.(distToOccult < distToStar)) then
                sinang = sqrt(max(0.,(1.-dotprod*dotprod)))
                if ((sinang*disttooccult) < occultradius) then 
                   occulted = .true.
                endif
             endif


             if (.not.occulted) then

                if (full) then
!                   call integratePath(real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                        grid%velocity(i1,i2,i3), &
!                        rVec, direction, grid, &
!                        lambdaArray, tauExt, tauAbs, tauSca, maxTau, nTau, .true., &
!                        escProb, .false., real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                        1, contTau, hitCore, .false.,.false., 0.,&
!                        VECTOR(0.,0.,0.), 1., .true., m, n, real(fstrength(m,n)), real(gDegen(m)), real(gDegen(n)), &
!                        tauLocal)
!                   tau_cmn = tauLocal + tauExt(nTau-1)
                else


                   tau_cmn = (pi*echarge**2)/(melectron*cspeed)
                   tau_cmn = tau_cmn * gDegen(m) * fStrength(m,n)
                   if (grid%adaptive) then
                      tau_cmn = tau_cmn * abs((thisOctal%N(thisSubcell,m)/gDegen(m)) - &
                                              (thisOctal%N(thisSubcell,n)/gDegen(n))) ! eq 5.
                   else 
                      tau_cmn = tau_cmn * abs((grid%n(i1,i2,i3,m)/gDegen(m)) - &
                                              (grid%n(i1,i2,i3,n)/gDegen(n))) ! eq 5.
                   endif
                   if (tau_cmn < 0.) then
                      tau_cmn = 0. !abs(tau_cmn)
                   endif
                   tau_cmn = tau_cmn  / (cSpeed / lambdaTrans(m,n))
                   if (grid%adaptive) then
                      if (grid%geometry(1:4) == "jets")  then
                         tau_cmn = tau_cmn/dV_dn_jets(rVec, direction)
                      else
                         rVecOctal = rVec
                         
                         octalCopy => thisOctal
!label1 = thisOctal%label(1)                
                         tau_cmn = tau_cmn / (amrGridDirectionalDeriv(grid,rVecOctal,direction, &
                                                                   startOctal=octalCopy) / 1.e10)
!if (label1 /= thisOctal%label(1)) print *, "Label change during amrGridDirectionalDeriv"                
                      end if
                   else
                      tau_cmn = tau_cmn / (directionalderiv(grid,rVec,i1,i2,i3,direction)/1.e10)
                   end if
                endif

                if (tau_cmn < 0.1) then
                   escProb = 1.d0-tau_cmn*0.5d0*(  1.0d0 - tau_cmn/3.0d0*( 1.0d0 - &
                        tau_cmn*0.25d0*(1.0d0 -0.20d0*tau_cmn) )  )
                else
                   escprob = (1.d0 - exp(-tau_cmn))/tau_cmn
                endif
                beta_cmn = beta_cmn + escprob * domega
                totOmega = totOmega + dOmega
             endif
          endif
       enddo
    enddo
    beta_cmn =   beta_cmn / fourpi
!    write(*,*) "calc",fourPi*(pi*starradius**2)/(fourPi * disttostar**2)
!    write(*,*) "found",totOmega
!if (SUM(labels-thisOctal%label) /= 0) print *, "            leaving beta_cmn, labels = ", labels    
  end function beta_cmn

  subroutine beta_cmn_sub(m,n,i1,i2,i3,grid,nstar)

    integer :: m,n
    type(gridtype) :: grid
    type(vector) :: rvec
    integer :: i1, i2, i3
    integer :: i, j, nphi=5, ntheta=5
    real :: r, thetatostar, phiTostar,betacmn
    real :: h, totomega, disttooccult
    integer nstar
    type(vector) :: tostar, occultposition, starposition
    type(vector) :: tooccult, direction
    real :: starradius, occultradius
    real :: disttostar
    real :: cosang, sinang, ang
    logical :: occulted, full
    real :: tau_cmn
    real :: dtheta, dphi
    real :: theta, phi, domega
    real :: escprob, dotprod
    integer, parameter :: maxTau = 4000
    integer :: nTau
    real, allocatable :: contTau(:,:)
    real :: tauExt(maxTau)
    real :: tauAbs(maxTau)
    real :: tauSca(maxTau)
    real :: lambdaArray(maxTau)
    logical :: hitCore
    real :: tauLocal

    full = .false.

    rVec = VECTOR(grid%xAxis(i1),grid%yAxis(i2),grid%zAxis(i3))

    if (nstar == 1) then
       starradius = grid%rstar1
       occultradius = grid%rstar2
       starposition = grid%starpos1
       occultposition = grid%starpos2
    else
       starradius = grid%rstar2
       occultradius = grid%rstar1
       starposition = grid%starpos2
       occultposition = grid%starpos1
    endif

    tostar = starposition - rvec
    disttostar = modulus(tostar)
    tostar = tostar / disttostar

    tooccult = occultposition - rvec
    disttooccult = modulus(tooccult)
    tooccult = tooccult / disttooccult

    h  = sqrt(max(0.,(disttostar**2 - starradius**2)))
    cosang = h / disttostar
    ang = acos(min(1.,max(-1.,cosAng)))

    call getPolar(toStar, r, thetaTostar, phiToStar)

    dtheta = pi / real(ntheta-1)
    dphi = twopi / real(nphi-1)
    betacmn = 0.
    totomega = 0.
    do i = 1, ntheta
       theta = thetaToStar + (2.*real(i-1)/real(nTheta-1)-1.)*ang
       if (theta > pi) theta = theta - pi
       if (theta < 0.) theta = theta + pi
       dTheta = 2.*ang/real(nTheta-1)
       do j = 1, nphi
          phi = phiToStar + (2.*real(j-1)/real(nPhi-1)-1.)*ang
          if (phi < 0.) phi = phi + twoPi
          if (phi > twoPi) phi = phi + twoPi
          dphi = 2.*ang/real(nPhi-1)
          occulted = .false.
          direction = vector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          domega = sin(theta)*dtheta*dphi
          dotprod = direction .dot. tostar
          escProb = 0.
          if ((dotprod > 0.) .and. (acos(min(1.,max(-1.,dotprod))) < ang)) then

             dotprod = direction .dot. tooccult
             if ((dotprod > 0.).and.(distToOccult < distToStar)) then
                sinang = sqrt(max(0.,(1.-dotprod*dotprod)))
                if ((sinang*disttooccult) < occultradius) then 
                   occulted = .true.
                endif
             endif


             if (.not.occulted) then

                if (full) then
!                   call integratePath(real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                        grid%velocity(i1,i2,i3), &
!                        rVec, direction, grid, &
!                        lambdaArray, tauExt, tauAbs, tauSca, maxTau, nTau, .true., &
!                        escProb, .false., real(lambdaTrans(m,n))*1.e8, real(lambdaTrans(m,n))*1.e8, &
!                        1, contTau, hitCore, .false.,.false., 0.,&
!                        VECTOR(0.,0.,0.), 1., .true., m, n, real(fstrength(m,n)), real(gDegen(m)), real(gDegen(n)), &
!                        tauLocal)
!                   tau_cmn = tauLocal + tauExt(nTau-1)
                else


                   tau_cmn = (pi*echarge**2)/(melectron*cspeed)
                   tau_cmn = tau_cmn * gDegen(m) * fStrength(m,n)
                   tau_cmn = tau_cmn * abs((grid%n(i1,i2,i3,m)/gDegen(m)) - (grid%n(i1,i2,i3,n)/gDegen(n))) ! eq 5.
                   if (tau_cmn < 0.) then
                      tau_cmn = 0. !abs(tau_cmn)
                   endif
                   tau_cmn = tau_cmn  / (cSpeed / lambdaTrans(m,n))
                   tau_cmn = tau_cmn / (directionalderiv(grid,rvec,i1,i2,i3,direction)/1.e10)
                endif

                if (tau_cmn < 0.1) then
                   escProb = 1.-tau_cmn*0.5*(  1.0 - tau_cmn/3.0*( 1.0 - &
                        tau_cmn*0.25*(1.0 -0.20*tau_cmn) )  )
                else
                   escprob = (1. - exp(-tau_cmn))/tau_cmn
                endif
                betacmn = betacmn + escprob * domega
                totOmega = totOmega + dOmega
             endif
          endif
       enddo
    enddo
    betacmn =   betacmn / fourpi
!    write(*,*) "calc",fourPi*(pi*starradius**2)/(fourPi * disttostar**2)
!    write(*,*) "found",totOmega
  end subroutine beta_cmn_sub



  ! this subroutine fills the grid with nlte level populations
  ! the grid should have temperatures/densities at all points

  subroutine initgridstateq(grid, contfile1, contfile2, popFileName, &
       readPops, writePops, lte, nLower, nUpper)
    implicit none
    type(gridtype) :: grid
    integer, parameter :: maxLevels = statEqMaxLevels
    integer :: nLower, nUpper
    logical :: lte, debugInfo
    real :: sinTheta
    logical :: readPops, writePops
    character(len=*) :: popFileName
    real :: visFrac1, visFrac2
    logical :: isBinary
    integer :: i, j, k, m
    real(kind=doubleKind) :: freq
    real :: hnu1(2000), hnu2(2000)
    real :: nuarray1(2000), nuarray2(2000)
    real :: percentDone
    integer :: nnu1, nnu2
    character(len=*) :: contfile1, contfile2
    real(kind=doublekind), parameter :: tolx = 1.d-3
    real(kind=doublekind), parameter :: tolf = 1.d-3
    real(kind=doublekind) :: x(15), nTot
    real :: t1, t2, t3
    integer :: i1, i2, i3
    type(vector) :: rvec, rHat, thisVec
    real :: departCoeff(maxLevels), ang, r
    real, allocatable :: departCoeffall(:)
    logical :: oneD, twoD, threeD, firstTime

    integer :: nIter, iIter
    real :: oldLevels(maxLevels+1)
    real(kind=doubleKind), allocatable :: xall(:)
    real (kind=doubleKind):: ne1, ne2
    real :: temp, crate
    real(kind=doubleKind) :: phiT
    logical :: ok
    logical :: aboutX, aboutZ
    integer :: ierr


    debugInfo = .true.
    oneD = .false.
    twoD = .false.
    threeD = .true.
    aboutX = .false.
    aboutZ = .true.


    isBinary = .false.
    if (grid%geometry == "binary") then
       isBinary = .true.
       threeD = .false.
       twoD = .true.
       aboutX = .true.
       aboutZ = .false.
    endif

    if (grid%geometry == "wind") then
       threeD = .false.
       oneD = .true.
       aboutZ = .false.
    endif
    
    if (grid%geometry == "donati") then
       twoD = .true.
       threeD = .false.
       aboutZ = .true.
       aboutX = .false.
    endif

    if (grid%geometry == "ttauri") then
       twoD = .true.
       threeD = .false.
       aboutZ = .true.
       aboutX = .false.
    endif

    if (grid%geometry == "puls") then
       twoD = .true.
       threeD = .false.
       aboutZ = .true.
       aboutX = .false.
    endif

    if (lte) then
       nIter = 1
    else
       nIter = 1
    endif

    do k = 2, maxlevels
       do i=1, k-1
          lambdaTrans(i, k) = lambdaTrans(k, i)
          ! aEinstein(k, i) = aEinstein(k, i) commented out by nhs
          freq = cspeed / lambdaTrans(i, k)
          bEinstein(k, i) = ((dble(aEinstein(k,i))*dble(cspeed)**2) / (2.d0*dble(hcgs)*dble(freq)**3))
          bEinstein(i, k) = bEinstein(k,i) * gDegen(k) / gDegen(i)
          fStrength(k, i) = -gDegen(i) * fStrength(i,k) / gDegen(k)
       enddo
    enddo





    allocate(grid%n(1:grid%na1, 1:grid%na2, 1:grid%na3,1:maxlevels))
    allocate(grid%ne(1:grid%na1, 1:grid%na2, 1:grid%na3))
    allocate(grid%nTot(1:grid%na1, 1:grid%na2, 1:grid%na3))

    grid%n = 1.e-20
    grid%ne = 1.

    if (readPops) then
       call readGridPopulations(popFilename, grid, maxLevels)
    else
       do i = 1, grid%na1
          do j = 1, grid%na2 
             do k = 1, grid%na3
                if ((.not.grid%instar(i,j,k)).and.grid%inUse(i,j,k)) then
                   nTot = dble(grid%rho(i,j,k) / mhydrogen)


                   if (grid%temperature(i,j,k) < 1.e5) then
                      phiT = 1.5d0*log10(2.d0*pi*mElectron*kErg*grid%temperature(i,j,k)) - &
                           3.d0*log10(hCgs) + log10(exp(-13.6d0/(kEV*grid%temperature(i,j,k))))
                     
                      phiT = 10.d0**phiT
                      call solveQuadDble(1.d0, phiT, &
                           -1.d0*phiT*dble(nTot), ne1, ne2, ok)
                      grid%ne(i,j,k) = min(max(ne1,ne2),nTot)   
                      grid%ne(i,j,k) = max(grid%ne(i,j,k),1.d0)
                   else
                      grid%ne(i,j,k) = dble(ntot)
                   endif

                   grid%nTot(i,j,k) = nTot 
                   do m = 1, maxlevels
                      grid%n(i,j,k,m) = boltzsaha(m, grid%ne(i,j,k), &
                           dble(grid%temperature(i,j,k)))
                   enddo
                endif
             enddo
          enddo
       enddo




       open(20,file=contfile1,status="old",form="formatted")
       nnu1 = 1
10     continue
       read(20,*,end=20) nuarray1(nnu1), hnu1(nnu1)
       nnu1 = nnu1 + 1
       goto 10
20     continue
       nnu1 = nnu1  - 1
9      close(20)
       hnu1(1:nnu1) = hnu1(1:nnu1) / fourPi   ! Converts from flux to flux momemnt

       if (grid%geometry == "binary") then
          open(20,file=contfile2,status="old",form="formatted")
          nnu2 = 1
30        continue
          read(20,*,end=40) nuarray2(nnu2), hnu2(nnu2)
          nnu2 = nnu2 + 1
          goto 30
40        continue
          nnu2 = nnu2  - 1
          close(20)
          hnu2(1:nnu2) = hnu2(1:nnu2) / fourPi   ! Converts from flux to flux momemnt
       endif

       if (threeD) then

          percentDone = 0.
          departCoeff = 1.

          !$OMP PARALLEL DO &
          !$OMP DEFAULT(NONE) &
          !$OMP PRIVATE(i1, i2, i3) &
          !$OMP SHARED(grid) &
          !$OMP SHARED(lte, tolx, tolf, hnu1, hnu2, nuArray1, nnu1, nuarray2, nnu2) &
          !$OMP SHARED(isBinary, debugInfo) &
          !$OMP PRIVATE(i, visFrac1, visFrac2, rVec, nTot, percentDone) &
          !$OMP PRIVATE(departCoeffAll, xAll) 
          do i1 = 1, grid%na1
             allocate(departCoeffAll(1:maxlevels))
             departCoeffAll = 1.d0
             do i2 = 1, grid%na2
                do i3 = 1, grid%na3
                   allocate(xall(1:15))
                   if (grid%inUse(i1,i2,i3)) then
                      rvec = vector(grid%xaxis(i1), grid%yaxis(i2), grid%zaxis(i3))
                      if (.not.grid%inStar(i1,i2,i3)) then
                         do i = 1, maxlevels
                            xall(i) = grid%n(i1,i2,i3,i) * dble(departCoeffall(i))
                         enddo
                         xall(maxlevels+1) = grid%ne(i1,i2,i3)
                         if (.not.lte) then

                            if (grid%geometry == "binary") then
                               call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                                    grid%starPos2, grid%rStar2, visFrac1)
                               call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                                    grid%starPos1, grid%rStar1, visFrac2)
                            else
                               visFrac1 = 1.
                               visFrac2 = 0.
                            endif


                            call mnewt(grid, 20,xall,maxlevels+1,tolx,tolf, hnu1, &
                                 nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
                            if ((.not.lte).and.debugInfo) write(*,'(a,3i3,f8.1)') "grid cell",i1,i2,i3,grid%temperature(i1,i2,i3)
                            write(*,*) visFrac1, visFrac2
                            grid%ne(i1,i2,i3) = xall(maxlevels+1)
                            do i = 1 , maxLevels
                               departCoeffall(i) = real(xall(i))/boltzSaha(i, grid%Ne(i1,i2,i3), dble(grid%temperature(i1,i2,i3)))
                               grid%n(i1,i2,i3,i) = xall(i)
                               if ((.not.lte).and.debugInfo) then
                                  write(*,'(i3,1p,e12.3,e12.3)') &
                                       i,departCoeffall(i), xall(i)
                               endif

                            enddo
                            if ((.not.lte).and.debugInfo) write(*,'(a,f6.2)') "log Ne: ",log10(grid%ne(i1,i2,i3))
                            nTot = dble(grid%rho(i1,i2,i3) / mhydrogen)
                            if ((.not.lte).and.debugInfo) write(*,'(a,1p,e12.5)') "Ne / Ntot: ",grid%ne(i1,i2,i3)/nTot
                         endif
                      endif
                   endif
                   deallocate(xall)
                enddo
             enddo
             percentDone = 100.*(real(i1)/real(grid%nx))
             if (.not.lte) write(*,'(a,f5.1)') "Percentage complete: ",percentDone
             deallocate(departCoeffAll)

          enddo
          !$OMP END PARALLEL DO

       endif

       ! if we are only computing in 2D then compute the populations in the z=0 plane
       ! and remap the populations
       ! onto the full 3D grid via rotation about the x-axis

       if (twoD) then


          if (aboutX) then
             do iIter = 1, nIter
                do i1 = 1, grid%nx
                   departCoeff = 1.
                   do i2 = 1, grid%ny
                      rvec = vector(grid%xaxis(i1), grid%yaxis(i2), 0.)
                      call hunt(grid%zAxis, grid%nz, 0., i3)

                      if (.not.grid%inStar(i1,i2,i3)) then
                         if (modulus(rVec-grid%starpos1) /= 0.) then
                            if (.not.lte) then
                               departCoeff(1) = 1./(0.5*(1.-sqrt(max(0.,(1.-grid%rStar1**2/modulus(rVec-grid%starpos1)**2)))))
                            endif
                         endif
                         if (iIter == 1) then
                            do i = 1, maxlevels
                               x(i) = grid%n(i1,i2,i3,i) * dble(departCoeff(i))
                            enddo
                            x(maxlevels+1) = grid%ne(i1,i2,i3)
                            oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                            oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                         else
                            oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                            oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                            do i = 1, maxlevels
                               x(i) = grid%n(i1,i2,i3,i)
                            enddo
                            x(maxlevels+1) = grid%ne(i1,i2,i3)
                         endif


                         if (.not.lte) then
                            if (grid%geometry == "binary") then
                               call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                                    grid%starPos2, grid%rStar2, visFrac1)
                               call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                                    grid%starPos1, grid%rStar1, visFrac2)
                            else
                               visFrac1 = 1.
                               visFrac2 = 0.
                            endif

                            call mnewt(grid, 20,x,maxlevels+1,tolx,tolf, hnu1, &
                                 nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
                         endif
                         if (.not.lte) write(*,*) "Grid: ",i1,i2,i3
                         do i = 1 , maxLevels
                            departCoeff(i) = real(x(i))/boltzSaha(i, grid%Ne(i1,i2,i3), dble(grid%temperature(i1,i2,i3)))

                            if (.not.lte) write(*,'(i3,1p,e12.3,e12.3)') i,departCoeff(i), x(i)/oldLevels(i)
                            grid%n(i1,i2,i3,i) = x(i)
                         enddo
                         grid%ne(i1,i2,i3) = x(maxLevels+1)
                         write(*,'(a,1p,e12.3,e12.3)') "Ne ",grid%ne(i1,i2,i3),grid%ne(i1,i2,i3)/oldLevels(maxLevels+1)
                      endif
                   enddo
                   percentDone = 100.*(real(i1)/real(grid%nx))
                   if (.not.lte) write(*,'(a,f5.1)') "Percentage complete: ",percentDone
                enddo
             enddo

             ! now perform the remapping

             do i = 1, grid%nx
                do j = 1, grid%ny
                   do k = 1, grid%nz
                      if (k /= i3) then
                         rVec = vector(grid%xaxis(i), grid%yaxis(j), grid%zaxis(k))
                         ang = atan2(rVec%z, rVec%y)
                         thisVec = rotateX(rVec, ang)
                         call getIndices(grid, thisVec, i1, i2, i3, t1, t2, t3)
                         if (.not.grid%inStar(i1,i2,i3)) then
                            grid%ne(i,j,k) = grid%ne(i1,i2,i3)
                            do m = 1, maxLevels
                               grid%n(i,j,k,m) = grid%n(i1,i2,i3,m)
                            enddo
                         endif
                      endif
                   enddo
                enddo
             enddo

          endif

          if (aboutZ) then
             if (.not.grid%cartesian) then


                WRITE(*,*) "about Z in polar"
                do iIter = 1, nIter
                   !$OMP PARALLEL DO &
                   !$OMP DEFAULT(NONE) &
                   !$OMP PRIVATE(i1, i2, i3) &
                   !$OMP SHARED(grid) &
                   !$OMP SHARED(lte, tolx, tolf, hnu1, hnu2, nuArray1, nnu1, nuarray2, nnu2) &
                   !$OMP SHARED(isBinary, debugInfo, iiter) &
                   !$OMP PRIVATE(i, visFrac1, visFrac2, rVec, nTot, percentDone) &
                   !$OMP PRIVATE(departCoeffAll, xAll, oldLevels, sinTheta, ang) 

                   do i1 = 1, grid%nr
                      allocate(departCoeffAll(1:maxlevels))
                      allocate(xall(1:maxLevels+1))
                      do i2 = 1, grid%nmu
                         departCoeffAll = 1.d0
                         ang = 0.
                         sinTheta = sqrt(1.d0 - grid%muAxis(i2)**2)
                         rvec = vector(grid%raxis(i1)*sinTheta*cos(ang), &
                              grid%raxis(i1)*sinTheta*sin(ang), &
                              grid%rAxis(i1)*grid%muAxis(i2))
                         call hunt(grid%phiAxis, grid%nphi, ang, i3)



                         if (.not.grid%inStar(i1,i2,i3).and.grid%inUse(i1,i2,i3)) then
                            if (iIter == 1) then
                               do i = 1, maxlevels
                                  xall(i) = grid%n(i1,i2,i3,i) * dble(departCoeffAll(i))
                               enddo
                               xall(maxlevels+1) = grid%ne(i1,i2,i3)
                               oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                               oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                            else
                               oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                               oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                               do i = 1, maxlevels
                                  xall(i) = grid%n(i1,i2,i3,i)
                               enddo
                               xall(maxlevels+1) = grid%ne(i1,i2,i3)
                            endif


                            if (.not.lte) then
                               if (grid%geometry == "binary") then
                                  call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                                       grid%starPos2, grid%rStar2, visFrac1)
                                  call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                                       grid%starPos1, grid%rStar1, visFrac2)
                               else
                                  visFrac1 = 1.
                                  visFrac2 = 0.
                               endif


                               call mnewt(grid, 20,xall,maxlevels+1,tolx,tolf, hnu1, &
                                    nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
                            endif
                            if ((.not.lte).and.debuginfo) write(*,*) "Grid: ",i1,i2,i3,grid%temperature(i1,i2,i3)
                            do i = 1 , maxLevels
                               departCoeffall(i) = real(xall(i))/boltzSaha(i, xall(maxlevels+1), dble(grid%temperature(i1,i2,i3)))

                               if ((.not.lte).and.debugInfo) then
                                  write(*,'(i3,1p,e12.3,e12.3, e12.3)') &
                                       i,departCoeffall(i), xall(i), abs(xall(i)-oldLevels(i))/xall(i)
                               endif
                               grid%n(i1,i2,i3,i) = xall(i)
                            enddo
                            grid%ne(i1,i2,i3) = xall(maxLevels+1)
                            if ((.not.lte).and.debugInfo) then
                               write(*,'(a,1p,e12.3,e12.3)') "Ne ",grid%ne(i1,i2,i3),grid%ne(i1,i2,i3)/oldLevels(maxLevels+1)
                            endif
                         if (minval(xall(1:maxlevels)) < 0.) then
                            write(*,'(a)') "Negative population..."
                            stop
                         endif
                      endif
                   enddo
                   percentDone = 100.*(real(i1)/real(grid%nr))
                   if ((.not.lte).and.debugInfo) write(*,'(a,f5.1)') "Percentage complete: ",percentDone
                   deallocate(departCoeffAll)
                   deallocate(xall)

                enddo
!$OMP END PARALLEL DO

             enddo

          else

             do iIter = 1, nIter
                do i1 = 1, grid%nx
                   departCoeff = 1.
                   do i3 = 1, grid%nz
                      rvec = vector(grid%xaxis(i1), 0., grid%zaxis(i3))
                      call hunt(grid%yAxis, grid%ny, 0., i2)

                      if (.not.grid%inStar(i1,i2,i3).and.grid%inUse(i1,i2,i3)) then
                         departCoeff = 1. !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                         if (iIter == 1) then
                            do i = 1, maxlevels
                               x(i) = grid%n(i1,i2,i3,i) * dble(departCoeff(i))
                            enddo
                            x(maxlevels+1) = grid%ne(i1,i2,i3)
                            oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                            oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                         else
                            oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                            oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                            do i = 1, maxlevels
                               x(i) = grid%n(i1,i2,i3,i)
                            enddo
                            x(maxlevels+1) = grid%ne(i1,i2,i3)
                         endif


                         if (.not.lte) then
                            if (grid%geometry == "binary") then
                               call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                                    grid%starPos2, grid%rStar2, visFrac1)
                               call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                                    grid%starPos1, grid%rStar1, visFrac2)
                            else
                               visFrac1 = 1.
                               visFrac2 = 0.
                            endif

                            call mnewt(grid, 20,x,maxlevels+1,tolx,tolf, hnu1, &
                                 nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
                         endif
                         if (.not.lte) write(*,*) "Grid: ",i1,i2,i3,grid%temperature(i1,i2,i3)
                         do i = 1 , maxLevels
                            departCoeff(i) = real(x(i))/boltzSaha(i, x(maxlevels+1), dble(grid%temperature(i1,i2,i3)))

                            if (.not.lte) then
                               write(*,'(i3,1p,e12.3,e12.3)') &
                                    i,departCoeff(i), x(i)
                            endif
                            grid%n(i1,i2,i3,i) = x(i)
                         enddo
                         grid%ne(i1,i2,i3) = x(maxLevels+1)
                         if (.not.lte) then
                            write(*,'(a,1p,e12.3,e12.3)') "Ne ",grid%ne(i1,i2,i3),grid%ne(i1,i2,i3)/oldLevels(maxLevels+1)
                         endif
                      endif
                   enddo
                   percentDone = 100.*(real(i1)/real(grid%nx))
                   if (.not.lte) write(*,'(a,f5.1)') "Percentage complete: ",percentDone
                enddo
             enddo
          endif
             
             ! now perform the remapping

          if (.not.lte) then
          if (grid%cartesian) then
             do i = 1, grid%nx
                do j = 1, grid%ny
                   do k = 1, grid%nz
                      if (j /= i2) then
                         rVec = vector(grid%xaxis(i), grid%yaxis(j), grid%zaxis(k))
                         ang = atan2(rVec%y, rVec%x)
                         thisVec = rotateZ(rVec, ang)
                         call getIndices(grid, thisVec, i1, i2, i3, t1, t2, t3)
                         if (.not.grid%inStar(i1,i2,i3)) then
                            grid%ne(i,j,k) = grid%ne(i1,i2,i3)
                            do m = 1, maxLevels
                               grid%n(i,j,k,m) = grid%n(i1,i2,i3,m)
                            enddo
                         endif
                      endif
                   enddo
                enddo
             enddo
             else

                write(*,'(a)') "Remapping about Z axis..."
                do k = 2, grid%nphi
                   grid%ne(1:grid%nr,1:grid%nMu,k) = &
                        grid%ne(1:grid%nr,1:grid%nMu,1)
                   grid%n(1:grid%nr,1:grid%nmu,k,1:maxLevels) = &
                        grid%n(1:grid%nr,1:grid%nMu,1,1:maxLevels)
                enddo
             endif
       endif
    endif
    
 endif



 if (oneD) then
    
    nIter = 1.
    
    departCoeff = 1.
    do iIter = 1, nIter
       do i1 = 1,grid%na1
          i2 = 1
          i3 = 1
          sinTheta = sqrt(1.d0 - grid%muAxis(i2)**2)
          rvec = vector(grid%raxis(i1)*sinTheta*cos(ang), &
               grid%raxis(i1)*sinTheta*sin(ang), &
               grid%rAxis(i1)*grid%muAxis(i2))
          
          where (departCoeff < 0.) 
             departCoeff = 1.
          end where
          departCoeff(1) = 1./(0.5*(1.-sqrt(max(0.,(1.-grid%rStar1**2/modulus(rVec-grid%starpos1)**2)))))
          
          
          if (.not.grid%inStar(i1,i2,i3).and.grid%inUse(i1,i2,i3)) then
             
             if (iIter == 1) then
                do i = 1, maxlevels
                   x(i) = grid%n(i1,i2,i3,i) * dble(departCoeff(i))
                enddo
                x(maxlevels+1) = grid%ne(i1,i2,i3)
                oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
             else
                oldLevels(1:maxLevels) = grid%n(i1,i2,i3,1:maxLevels)
                oldLevels(maxLevels+1) = grid%ne(i1,i2,i3)
                do i = 1, maxlevels
                   x(i) = grid%n(i1,i2,i3,i)
                enddo
                x(maxlevels+1) = grid%ne(i1,i2,i3)
             endif
             
             
             if (.not.lte) then
                if (grid%geometry == "binary") then
                   call occultTest(grid, i1, i2, i3, grid%starPos1, grid%rStar1, &
                        grid%starPos2, grid%rStar2, visFrac1)
                   call occultTest(grid, i1, i2, i3, grid%starPos2, grid%rStar2, &
                        grid%starPos1, grid%rStar1, visFrac2)
                else
                   visFrac1 = 1.
                   visFrac2 = 0.
                endif
                
                call mnewt(grid, 20,x,maxlevels+1,tolx,tolf, hnu1, &
                     nuarray1, nnu1, hnu2, nuarray2, nnu2, rvec, i1, i2, i3, visFrac1, visFrac2, isBinary)
             endif
             if (.not.lte) write(*,*) "Grid: ",i1,i2,i3,grid%temperature(i1,i2,i3)
             do i = 1 , maxLevels
                departCoeff(i) = real(x(i))/boltzSaha(i, x(maxlevels+1), dble(grid%temperature(i1,i2,i3)))
                
                if (.not.lte) then
                   write(*,'(i3,1p,e12.3,e12.3,e12.3,e12.3)') &
                        i,departCoeff(i), x(i),log10(departCoeff(i)),log10(x(2)/x(1))
                endif
                grid%n(i1,i2,i3,i) = x(i)
             enddo
             grid%ne(i1,i2,i3) = x(maxLevels+1)
             if (.not.lte) then
                write(*,'(a,1p,e12.3,e12.3)') "Ne ",grid%ne(i1,i2,i3),grid%ne(i1,i2,i3)/grid%nTot(i1,i2,i3)
             endif

             endif
             enddo
          enddo
             ! now perform the remapping
             rHat = VECTOR(1., 1., 1.)
             call normalize(rHat)

             do i = 1, grid%na1
                do j = 1, grid%na2
                   do k = 1, grid%na3
                      i1 = i
                      i2 = 1
                      i3 = 1
                      if (.not.grid%inStar(i1,i2,i3)) then
                         grid%ne(i,j,k) = grid%ne(i1,i2,i3)
                         do m = 1, maxLevels
                            grid%n(i,j,k,m) = grid%n(i1,i2,i3,m)
                         enddo
                      endif
                   enddo
                enddo
             enddo

          endif







          if (writePops) then
             call writeGridPopulations(popFilename, grid, maxLevels)
          endif
    endif

    write(*,'(a,f8.1)') "Generating opacities for ",lambdaTrans(nLower, nUpper)*1.e8


!    write(*,*) "!!!!!!!!!!!!!!! depart coeff of level 3 is 2"
!    grid%n(:,:,:,3) = grid%n(:,:,:,3) * 100.
!    grid%n(:,:,:,5) = grid%n(:,:,:,5) * 30.
    call generateOpacities(grid, nLower, nUpper)

    deallocate(grid%n)
    deallocate(grid%ne)
    deallocate(grid%nTot)



  end subroutine initGridStateq



  real pure function BoltzSaha(m, Ne, t)
  
    integer, intent(in)              :: m
    real(kind=doubleKind), intent(in):: Ne, t
    
    real(kind=doubleKind), parameter :: ci = 2.07d-16
    real(kind=doubleKind), parameter :: iPot = 13.598d0

    BoltzSaha = Ne**2 * gDegen(m) * ci * exp( (iPot-eTrans(m))/(kev * t) ) / t**1.5
    
  end function BoltzSaha


  real pure function boltzmann(m, N0, t)
  
    integer, intent(in)              :: m
    real(kind=doubleKind), intent(in):: N0, t
    real                             :: z0

    z0 = SUM(gDegen(1:15)*exp(-eTrans(1:15)/(kev*t)))
    boltzmann = (gDegen(m)/z0)*n0*exp(-eTrans(m)/(kev*t))

  end function boltzmann


  real(kind=doubleKind) pure function cikt(i,t)
    !
    ! this function calculates the collisional ionization rate (see k&c and
    ! mihalas 1967 (apj 149 169) ).
    ! 
    integer, intent(in)              :: i        ! the level
    real(kind=doubleKind),intent(in) :: t        ! the temperature
    
    real(kind=doubleKind) :: t1                  
    real(kind=doubleKind) :: gammait             ! see k&c
    real(kind=doubleKind) :: lgt
    real(kind=doubleKind) :: chi                 ! potential from i to k
    real(kind=doubleKind) :: cint(5,10)
    ! making cint a PARAMETER may cause problems with XL Fortran
    cint = reshape(source=                                                                   &
         (/-0.4350000e0_db, 0.3000000e00_db,0.0000000e0_db, 0.0000000e0_db,0.00000000e0_db,  &
            1.9987261e1_db,-5.8906298e-5_db,0.0000000e0_db,-2.8185937e4_db,5.44441600e7_db,  &
            1.3935312e3_db,-1.6805859e02_db,0.0000000e0_db,-2.5390000e3_db,0.00000000e0_db,  &
            2.0684609e3_db,-3.3415820e02_db,0.0000000e0_db, 0.0000000e0_db,-7.6440625e3_db,  &
            3.2174844e3_db,-5.5882422e02_db,0.0000000e0_db, 0.0000000e0_db,-6.8632500e3_db,  &
            5.7591250e3_db,-1.5163125e03_db,8.1750000e1_db, 0.0000000e0_db,0.00000000e0_db,  &
            1.4614750e4_db,-4.8283750e03_db,3.9335938e2_db, 0.0000000e0_db,0.00000000e0_db,  &
            2.8279250e4_db,-1.0172750e04_db,9.1967968e2_db, 0.0000000e0_db,0.00000000e0_db,  &
            4.6799250e4_db,-1.7627500e04_db,1.6742031e3_db, 0.0000000e0_db,0.00000000e0_db,  &
           -7.4073000e4_db, 6.8599375e03_db,0.0000000e0_db, 2.0169800e5_db,0.00000000e0_db/),&
                                                                             shape=(/5,10/))

    t1 = min(t,1.5_db)
    lgt=log10(t1)
    if (i .ne. 2) then 
       gammait=cint(1,i)+cint(2,i)*lgt+cint(3,i)*(lgt**2)+ &
            (cint(4,i)/lgt)+(cint(5,i)/(lgt**2))
    else
       gammait=cint(1,i)+cint(2,i)*t1+(cint(4,i)/t1)+(cint(5,i)/t1**2)
    endif
    chi=13.598e0-eTrans(i)
    cikt=(5.465e-11)*sqrt(t1)*exp(-chi/(kev*t1))*gammait

  end function cikt


  real(kind=doubleKind) pure function cijt(i,j,t)
    !
    ! this function calculates the excitation rate according to crandell et al
    ! (apj 191 789).
    !
    integer, intent(in)              :: i,j  ! lower/upper level
    real(kind=doubleKind),intent(in) :: t    ! temperature
    
    integer :: r                             ! counter
    real(kind=doubleKind) :: t1              ! temperature
    real(kind=doubleKind) :: bsum(9)         ! summation
    real(kind=doubleKind) :: chi             ! excitation energy
    real(kind=doubleKind) :: x
    integer               :: i1,j1           ! lower/upper level
    real(kind=doubleKind) :: avals(7) 
    ! making avals a PARAMETER may cause problems with XL Fortran
    avals = (/ 2.579997e-10_db, -1.629166e-10_db, 7.713069e-11_db, -2.668768e-11_db, &
               6.642513e-12_db, -9.422885e-13_db, 0.e0_db                            /)

    t1 = min(t,1.e5_db)
    i1 = min(i,j)
    j1 = max(i,j)

    x=log10(t1)-4.e0
    chi=abs(eTrans(j)-eTrans(i))
    chi=chi/(kev*t1)
    if ((i1.eq.1).and.(j1.eq.2)) then
       bsum(9)=0.0
       bsum(8)=0.0
       do r=7,1,-1
          bsum(r)=2.0*x*bsum(r+1)-bsum(r+2)+avals(r)
       enddo
       cijt=4.e0*sqrt(t1)*exp(-chi)*0.5e0*(bsum(1)-bsum(3))
    else
       cijt=abs(5.465e-11*sqrt(t1)*4.e0_db*fStrength(i1,j1)/ &
            (( (1.e0/real(i1)**2) - (1.0/real(j1)**2))**2))*chi* &
            (expint(1,chi)+0.148d0*chi*expint(5,chi))
    endif

    if (j < i) then
       cijt = cijt / (exp(-chi) * 2.)
    endif


  end function cijt


  real(kind=doubleKind) pure function annu(n,nu)
    !
    ! this subroutine calculates the photoionization x-section
    ! for hydrogen from the n-th level for a given freq photon.
    !
    integer, intent(in)              :: n     ! the level
    real(kind=doubleKind), intent(in):: nu    ! the photon frequency
    real(kind=doubleKind)            :: lam,e ! the photon wavelength
    
    lam=cSpeed/nu
    lam=lam * 1.e8_db

    e = hCgs * nu * ergToEv

    if (e > (13.598 - eTrans(n))) then
       annu=1.044e-26_db*gii(n,1.e0_db,lam)*(lam**3)/dble(n)**5
    else
       annu = 1.e-30_db
    endif

  end function annu


  real pure function expint(n,x1)
    ! 
    ! exponential integrals. numerical approximations given by gray.
    !
    integer, intent(in)              :: n
    real(kind=doubleKind),intent(in) :: x1
    real(kind=doubleKind)            :: x,ep(5)
    integer                          :: i
    real(kind=doubleKind)            :: a1(4)
    real(kind=doubleKind)            :: b1(4)
    ! making a1 and b1 PARAMETERs may cause problems with XL Fortran
    a1=(/ 0.2677737343e0_db,8.6347608925e0_db, 18.059016973e0_db, 8.5733287401e0_db /)
    b1=(/ 3.9584969228e0_db,21.0996530827e0_db,25.6329561486e0_db,9.5733223454e0_db /)

    x=abs(x1)
    if (x.le.1.0) then
       ep(1)=-log(x)-0.57721566e0_db+0.99999193e0_db*x-0.24991055e0_db*x**2+ &
            0.05519968e0_db*x**3-0.00976004*x**4+0.00107857*x**5
    else
       ep(1)=(x**4+a1(4)*x**3+a1(3)*x**2+a1(2)*x+a1(1))
       ep(1)=ep(1)/(x**4+b1(4)*x**3+b1(3)*x**2+b1(2)*x+b1(1))
       ep(1)=ep(1)/(x*exp(x))
    endif
    if (n.gt.1) then
       do i=1,n-1
          ep(i+1)=(exp(-x)-x*ep(i))/real(i)
       enddo
    endif
    expint=ep(n)
  end function expint



  real(kind=doubleKind) pure function gii(n,z,wl)
    !
    ! returns bound-free gaunt factors (nicked from idh)
    !
    integer, intent(in)   :: n
    real(kind=doubleKind), intent(in) ::  z,wl
    real(kind=doubleKind) ::  coeff(6)
    real(kind=doubleKind) ::  efree,sum,ag,alam
    integer               :: i
    ! making coeff  a PARAMETER may cause problems with XL Fortran
    coeff = (/-0.338276d0, 0.808398d0, -0.59941d0, 0.104292d0, &
              -6.61998d-3, 1.15609d-4 /)


    efree=(911.76e0/wl-1.e0/real(n*n))/(z*z)
    if (efree.le.(1.e0+2.e0/n)) then
       gii = giia(n,z,wl)
       goto 500
    elseif (efree.gt.10.e0) then
       sum=coeff(1)
       ag=1.e0
       efree=log10(efree)
       do i=2,6
          ag=ag*efree
          sum=sum+coeff(i)*ag
       enddo
       gii=10.e0_db**sum
       goto 500
    else
       alam=911.76e0/(z*z*(1.e0+2.e0/n)+1.e0/(n*n))
       alam=giia(n,z,alam)
       sum=log10(1.e0+2.e0/n)
       efree=log10(efree)
       sum=(efree-sum)/(1.e0_db-sum)
       gii=(1.e0-sum)*alam+0.93e0*sum+0.2e0*sum*(1.e0-sum)
       goto 500
    endif
500 continue
  end function gii
  

  real pure function giia(n,z,wl)
    integer, intent(in)               :: n
    real(kind=doubleKind), intent(in) :: z, wl
    real(kind=doubleKind)             :: ag, u, gii, term
    ag=real(n)
    u=ag*ag*911.76e0_db/(wl*z*z)-1.e0_db
    gii=1.e0_db+0.1728e0_db*(u-1.e0_db)/((n*(u+1.e0_db))**(2.e0_db/3.e0_db))
    term=0.0496e0_db*(1.e0_db+u*(u+1.333e0_db))/(n*(u+1.e0_db)**(4.e0_db/3.e0_db))
    if ((term/gii).le.0.25e0_db) then
       giia=gii-term
       goto 500
    else
       giia=1.0
       goto 500
    endif
500 continue
  end function giia


  real(kind=doubleKind) function equation8(n, nPop, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, &
       rVec, i1, i2, i3, grid, visFrac1, visFrac2, binary, thisOctal, thisSubcell)
    type(GRIDTYPE), intent(in)  :: grid
    integer, intent(in)         :: nPop, nNu1, nNu2
    logical, intent(in)         :: binary
    real, intent(in)            :: visFrac1, visFrac2
    real, intent(in), dimension(:) :: Hnu1, nuArray1
    real, intent(in), dimension(:) :: Hnu2, nuArray2
    integer, intent(in)         :: i1, i2, i3
    type(VECTOR), intent(in)    :: rVec
    type(OCTAL),pointer,optional:: thisOctal
    integer,intent(in),optional :: thisSubcell
    
    integer :: m,n
    real(kind=doubleKind) :: nstar
    real(kind=doubleKind) :: freq
    integer, parameter :: debug = 0
    integer :: i
    real(kind=doubleKind) :: Inu
    real(kind=doubleKind) :: tot
    real(kind=doubleKind) :: fac1, fac2

!integer :: labels(8)
!labels = thisOctal%label
    if (n == debug) write(*,*) " "

    ! we have one IF statement to decide whether we use the block containin the
    !   AMR version of the code. 

    if (grid%adaptive) then
   
       tot = 0.e0_db
       do m = 1, n - 1
          freq = cSpeed/lambdaTrans(m,n)
          if ((freq >= nuArray1(1)).and.(freq <= nuArray1(nnu1))) then
             call hunt(nuArray1, nNu1, real(freq), i)
             Inu = 4.*logint(real(freq), nuArray1(i),nuArray1(i+1),hnu1(i), hnu1(i+1))
          else
             inu = 0.
          endif
   
   
          tot = tot  + (thisOctal%N(thisSubcell,m)*bEinstein(m,n) &
               - thisOctal%N(thisSubcell,n)*bEinstein(n,m)) * &
               beta_cmn(m, n, rVec, i1, i2, i3, grid, 1, thisOctal, thisSubcell) * Inu
   
          if (binary) then
             if ((freq >= nuArray2(1)).and.(freq <= nuArray2(nnu2))) then
                call hunt(nuArray2, nNu2, real(freq), i)
                Inu = 4.*logint(real(freq), nuArray2(i),nuArray2(i+1),hnu2(i), hnu2(i+1))
             else
                inu = 0.
             endif
             tot = tot  + (thisOctal%N(thisSubcell,m)*bEinstein(m,n) - &
                  thisOctal%N(thisSubcell,n)*bEinstein(n,m)) * &
                  beta_cmn(m, n, rVec, i1, i2, i3, grid, 2, thisOctal, thisSubcell) * Inu
          endif
   
          tot = tot - thisOctal%N(thisSubcell,n)*aEinstein(n,m)*&
                   beta_mn(m, n, rVec, i1, i2, i3, grid, thisOctal, thisSubcell)
          tot = tot + thisOctal%Ne(thisSubcell) * &
               (thisOctal%N(thisSubcell,m)*cijt(m,n,real(thisOctal%temperature(thisSubcell),kind=db)) &
               - thisOctal%N(thisSubcell,n) * cijt(n,m,real(thisOctal%temperature(thisSubcell),kind=db)))
       enddo
   
       if (n == debug) then
          write(*,*) "m < n",tot
       endif
   
   
       do m = n+1, nPop
          tot = tot + thisOctal%N(thisSubcell,m)*aEinstein(m,n)*&
                   beta_mn(n,m,rVec, i1, i2, i3,grid, thisOctal, thisSubcell)
          if (n == debug) then
             write(*,*) "spont",m,tot
          endif
   
          freq = cSpeed/lambdaTrans(m,n)
          if ((freq >= nuArray1(1)).and.(freq <= nuArray1(nnu1))) then
             call hunt(nuArray1, nNu1, real(freq), i)
             Inu = 4.*logint(real(freq), nuArray1(i),nuArray1(i+1),hnu1(i), hnu1(i+1))
          else
             inu = 0.
          endif
   
   
          tot = tot  - (thisOctal%N(thisSubcell,n)*bEinstein(n,m) - &
               thisOctal%N(thisSubcell,m)*bEinstein(m,n)) * &
               beta_cmn(n, m, rVec, i1, i2, i3, grid, 1, thisOctal, thisSubcell) * Inu
   
          if (n == debug) then
             write(*,*) m,"star 1",tot
          endif
   
   
          if (binary) then
             if ((freq >= nuArray2(1)).and.(freq <= nuArray2(nnu2))) then
                call hunt(nuArray2, nNu2, real(freq), i)
                Inu = 4.*logint(real(freq), nuArray2(i),nuArray2(i+1),hnu2(i), hnu2(i+1))
             else
                inu = 0.
             endif
             tot = tot  -(thisOctal%N(thisSubcell,n)*bEinstein(n,m) - &
                  thisOctal%N(thisSubcell,m)*bEinstein(m,n)) * &
                  beta_cmn(n, m, rVec, i1, i2, i3, grid, 2, thisOctal, thisSubcell) * Inu
          endif
   
          if (n == debug) then
             write(*,*) m,"star 2",tot
          endif
   
   
          tot = tot + thisOctal%Ne(thisSubcell) * (thisOctal%N(thisSubcell,m) &
           *cijt(m,n,real(thisOctal%temperature(thisSubcell),kind=db)) - thisOctal%N(thisSubcell,n)&
           *cijt(n,m,real(thisOctal%temperature(thisSubcell),kind=db)))
          if (n == debug) then
             write(*,*) m,"collisional",tot
          endif
   
   
       enddo
   
   
       if ( n == debug) then
          write(*,*) "m > n",tot
       endif

       NStar = boltzSaha(n, thisOctal%Ne(thisSubcell),dble(thisOctal%temperature(thisSubcell)))
   
       fac1 = integral1(n,hnu1, nuArray1, nNu1, rVec, i1, i2, i3, grid, 1, thisOctal, thisSubcell)*visFrac1
       if (binary) then
          fac2 = integral1(n,hnu2, nuArray2, nNu2, rVec, i1, i2, i3, grid, 2, thisOctal, thisSubcell)*visFrac2
       else
          fac2 = 0.
       endif
       tot = tot + Nstar * (fac1 + fac2 + thisOctal%Ne(thisSubcell)* &
          cikt(n,real(thisOctal%temperature(thisSubcell),kind=db)))
   
       if (n == debug) then
          write(*,*) "Recombination: ",tot
          write(*,*) nstar, fac1, fac2, thisOctal%Ne(thisSubcell)* &
             cikt(n,real(thisOctal%temperature(thisSubcell),kind=db))
       endif
       
       
   
       fac1 = integral2(n,hnu1, nuArray1, nNu1, rVec, grid, 1)*visFrac1
       if (binary) then
          fac2 = integral2(n,hnu2, nuArray2, nNu2, rVec, grid, 2)*visFrac2
       else
          fac2 = 0.
       endif
       tot = tot - thisOctal%N(thisSubcell,n)*(fac1 + fac2 + thisOctal%Ne(thisSubcell) * &
          cikt(n, real(thisOctal%temperature(thisSubcell),kind=db)))
   
       if (n == debug) then
          write(*,*) "Ionization: ", &
               thisOctal%N(thisSubcell,n)*(fac1 + fac2 + thisOctal%Ne(thisSubcell) * &
                  cikt(n, real(thisOctal%temperature(thisSubcell),kind=db)))
          write(*,*) thisOctal%N(thisSubcell,n),fac1,fac2,thisOctal%Ne(thisSubcell) * &
             cikt(n, real(thisOctal%temperature(thisSubcell),kind=db))
       endif 
       
       
    else ! grid is not adaptive
       
       
       tot = 0.e0_db
       do m = 1, n - 1
          freq = cSpeed/lambdaTrans(m,n)
          if ((freq >= nuArray1(1)).and.(freq <= nuArray1(nnu1))) then
             call hunt(nuArray1, nNu1, real(freq), i)
             Inu = 4.*logint(real(freq), nuArray1(i),nuArray1(i+1),hnu1(i), hnu1(i+1))
          else
             inu = 0.
          endif
   
   
          tot = tot  + (grid%N(i1,i2,i3,m)*bEinstein(m,n) &
               - grid%N(i1,i2,i3,n)*bEinstein(n,m)) * &
               beta_cmn(m, n, rVec, i1, i2, i3, grid, 1) * Inu
   
          if (binary) then
             if ((freq >= nuArray2(1)).and.(freq <= nuArray2(nnu2))) then
                call hunt(nuArray2, nNu2, real(freq), i)
                Inu = 4.*logint(real(freq), nuArray2(i),nuArray2(i+1),hnu2(i), hnu2(i+1))
             else
                inu = 0.
             endif
             tot = tot  + (grid%N(i1,i2,i3,m)*bEinstein(m,n) - &
                  grid%N(i1,i2,i3,n)*bEinstein(n,m)) * &
                  beta_cmn(m, n, rVec, i1, i2, i3, grid, 2) * Inu
          endif
   
          tot = tot - grid%N(i1,i2,i3,n)*aEinstein(n,m)*beta_mn(m, n, rVec, i1, i2, i3, grid)
          tot = tot + grid%Ne(i1,i2,i3) * &
               (grid%N(i1,i2,i3,m)*cijt(m,n,dble(grid%temperature(i1,i2,i3))) - grid%N(i1,i2,i3,n) * &
                cijt(n,m,dble(grid%temperature(i1,i2,i3))))
       enddo
   
       if (n == debug) then
          write(*,*) "m < n",tot
       endif
   
   
       do m = n+1, nPop
          tot = tot + grid%N(i1,i2,i3,m)*aEinstein(m,n)*beta_mn(n,m,rVec, i1, i2, i3,grid)
          if (n == debug) then
             write(*,*) "spont",m,tot
          endif
   
          freq = cSpeed/lambdaTrans(m,n)
          if ((freq >= nuArray1(1)).and.(freq <= nuArray1(nnu1))) then
             call hunt(nuArray1, nNu1, real(freq), i)
             Inu = 4.*logint(real(freq), nuArray1(i),nuArray1(i+1),hnu1(i), hnu1(i+1))
          else
             inu = 0.
          endif
   
   
          tot = tot  - (grid%N(i1,i2,i3,n)*bEinstein(n,m) - &
               grid%N(i1,i2,i3,m)*bEinstein(m,n)) * &
               beta_cmn(n, m, rVec, i1, i2, i3, grid, 1) * Inu
   
          if (n == debug) then
             write(*,*) m,"star 1",tot
          endif
   
   
          if (binary) then
             if ((freq >= nuArray2(1)).and.(freq <= nuArray2(nnu2))) then
                call hunt(nuArray2, nNu2, real(freq), i)
                Inu = 4.*logint(real(freq), nuArray2(i),nuArray2(i+1),hnu2(i), hnu2(i+1))
             else
                inu = 0.
             endif
             tot = tot  -(grid%N(i1,i2,i3,n)*bEinstein(n,m) - &
                  grid%N(i1,i2,i3,m)*bEinstein(m,n)) * &
                  beta_cmn(n, m, rVec, i1, i2, i3, grid, 2) * Inu
          endif
   
          if (n == debug) then
             write(*,*) m,"star 2",tot
          endif
   
   
          tot = tot + grid%Ne(i1,i2,i3) * (grid%N(i1,i2,i3,m) &
           *cijt(m,n,dble(grid%temperature(i1,i2,i3))) - grid%N(i1,i2,i3,n)*cijt(n,m,dble(grid%temperature(i1,i2,i3))))
          if (n == debug) then
             write(*,*) m,"collisional",tot
          endif
   
   
       enddo
   
   
       if ( n == debug) then
          write(*,*) "m > n",tot
       endif
   
       NStar = boltzSaha(n, grid%Ne(i1,i2,i3),dble(grid%temperature(i1,i2,i3)))
   
       fac1 = integral1(n,hnu1, nuArray1, nNu1, rVec, i1, i2, i3, grid, 1, thisOctal, thisSubcell)*visFrac1
       if (binary) then
          fac2 = integral1(n,hnu2, nuArray2, nNu2, rVec, i1, i2, i3, grid, 2, thisOctal, thisSubcell)*visFrac2
       else
          fac2 = 0.
       endif
       tot = tot + Nstar * (fac1 + fac2 + grid%Ne(i1,i2,i3)*cikt(n,dble(grid%temperature(i1,i2,i3))))
   
       if (n == debug) then
          write(*,*) "Recombination: ",tot
          write(*,*) nstar, fac1, fac2, grid%Ne(i1,i2,i3)*cikt(n,dble(grid%temperature(i1,i2,i3)))
       endif
       
       
   
       fac1 = integral2(n,hnu1, nuArray1, nNu1, rVec, grid, 1)*visFrac1
       if (binary) then
          fac2 = integral2(n,hnu2, nuArray2, nNu2, rVec, grid, 2)*visFrac2
       else
          fac2 = 0.
       endif
       tot = tot - grid%N(i1,i2,i3,n)*(fac1 + fac2 + grid%Ne(i1,i2,i3) * cikt(n, dble(grid%temperature(i1,i2,i3))))
   
       if (n == debug) then
          write(*,*) "Ionization: ", &
               grid%N(i1,i2,i3,n)*(fac1 + fac2 + grid%Ne(i1,i2,i3) * cikt(n, dble(grid%temperature(i1,i2,i3))))
          write(*,*) grid%N(i1,i2,i3,n),fac1,fac2,grid%Ne(i1,i2,i3) * cikt(n, dble(grid%temperature(i1,i2,i3)))
       endif 
       
    end if ! (grid%adaptive)

    if (n == debug) write(*,*) "Final total:",tot
    equation8 = tot

!if (SUM(labels-thisOctal%label) /= 0) print *, "          leaving equation8, labels = ", labels    

  end function equation8
 

  real (kind=doubleKind) pure function equation14(nPop, grid, i1, i2, i3, thisOctal, thisSubcell)

    type(GRIDTYPE), intent(in)   :: grid
    integer, intent(in)          :: i1, i2, i3, nPop
    type(OCTAL),pointer,optional :: thisOctal
    integer,intent(in),optional  :: thisSubcell
    
    real(kind=doubleKind) :: tot
    integer :: i
    integer, parameter :: maxLevels = statEqMaxLevels

    if (grid%adaptive) then

      tot = 0.e0_db
       do i = 1, nPop
          tot = tot + dble(thisOctal%N(thisSubcell,i))
       enddo
       
       do i = nPop+1, maxLevels+3
          tot = tot + boltzSaha(i, dble(thisOctal%Ne(thisSubcell)), dble(thisOctal%temperature(thisSubcell)))
       enddo
   !    write(*,*) "Neutral",tot
       tot = tot + thisOctal%Ne(thisSubcell)
   !    write(*,*) "Ne",grid%ne(i1,i2,i3)
       tot = tot - dble(thisOctal%rho(thisSubcell) / mHydrogen)
   !    write(*,*) "Ntot",grid%rho(i1,i2,i3) / mHydrogen
       equation14 = tot
   !    write(*,*) "tot",tot

    else ! grid not adaptive

      tot = 0.e0_db
       do i = 1, nPop
          tot = tot + dble(grid%n(i1,i2,i3,i))
       enddo
       
       do i = nPop+1, maxLevels+3
          tot = tot + boltzSaha(i, grid%ne(i1,i2,i3), dble(grid%temperature(i1,i2,i3)))
       enddo
   !    write(*,*) "Neutral",tot
       tot = tot + grid%ne(i1,i2,i3)
   !    write(*,*) "Ne",grid%ne(i1,i2,i3)
       tot = tot - dble(grid%rho(i1,i2,i3) / mHydrogen)
   !    write(*,*) "Ntot",grid%rho(i1,i2,i3) / mHydrogen
       equation14 = tot
   !    write(*,*) "tot",tot

    end if ! (grid%adaptive)

  end function equation14



  real function integral1(n,hnu, nuArray, nNu, rVec, i1, i2, i3, grid, nStar, thisOctal, thisSubcell)
    real, intent(in), dimension(:) :: hnu, nuArray
    integer, intent(in)      :: nNu, nStar, n
    type(VECTOR), intent(in) :: rVec
    type(GRIDTYPE),intent(in):: grid
    integer,intent(in)       :: i1, i2, i3
    type(OCTAL),pointer,optional:: thisOctal
    integer,intent(in),optional :: thisSubcell
    
    real    :: r
    integer :: i, imin
    real    :: fac1, fac2
    real    :: w, x1, freq, tot, jnu

    if (nStar == 1) then
       r = modulus(rVec-grid%starPos1)
       x1 = sqrt(max(0.,(1. - grid%rStar1**2 / r**2)))
       w = 0.5*(1. - x1)
    else
       r = modulus(rVec-grid%starPos2)
       x1 = sqrt(max(0.,1. - grid%rStar2**2 / r**2))
       w = 0.5*(1. - x1)
    endif

    freq = ((13.598-eTrans(n))*1.602192e-12)/hcgs
    if ((freq < nuArray(1)) .or.(freq > nuArray(nNu))) then
       write(*,*) "error in integral1",nNu,freq,nuArray(1),nuArray(nNu)
do ; enddo       
       jnu = 1.e-28
       iMin = 1
    else
       call hunt(nuArray, nNu, freq, iMin)
       jnu = 4.*w*logint(freq,nuArray(iMin), nuArray(iMin+1), hnu(imin), hnu(imin+1))
    endif
    tot = 0.

    if (grid%adaptive) then 
            
       fac1  = (fourPi/(hCgs*freq))*annu(n,dble(freq))* &
 ((2.*real(dble(hCgs)*dble(freq)**3) / cSpeed**2) + jnu)*exp(-hCgs*freq/(kErg*thisOctal%temperature(thisSubcell)))
       jnu = 4.*w*hnu(imin+1)
       fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))* &
            ((2.*real(dble(hcgs)*dble(nuArray(imin+1))**3) / cSpeed**2) + jnu) &
            *exp(-hcgs*nuArray(imin+1)/(kErg*thisOctal%temperature(thisSubcell)))
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
       do i = iMin+1, nNu-1
          Jnu = 4.*w*hnu(i)
          fac1  = (fourPi/(hCgs*nuArray(i)))*annu(n,dble(nuArray(i)))* &
 ((2.*real(dble(hCgs)*dble(nuArray(i))**3) / cSpeed**2) + jnu)*exp(-hCgs*nuArray(i)/(kErg*thisOctal%temperature(thisSubcell)))
          Jnu = 4.*w*hnu(i+1)
          fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))* &
 ((2.*real(dble(hcgs)*dble(nuArray(i+1))**3) / cSpeed**2) + jnu)*exp(-hcgs*nuArray(i+1)/(kErg*thisOctal%temperature(thisSubcell)))
          tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
       enddo

    else ! not adaptive
            
       fac1  = (fourPi/(hCgs*freq))*annu(n,dble(freq))* &
 ((2.*real(dble(hCgs)*dble(freq)**3) / cSpeed**2) + jnu)*exp(-hCgs*freq/(kErg*grid%temperature(i1,i2,i3)))
       jnu = 4.*w*hnu(imin+1)
       fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))* &
            ((2.*real(dble(hcgs)*dble(nuArray(imin+1))**3) / cSpeed**2) + jnu) &
            *exp(-hcgs*nuArray(imin+1)/(kErg*grid%temperature(i1,i2,i3)))
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
       do i = iMin+1, nNu-1
          Jnu = 4.*w*hnu(i)
          fac1  = (fourPi/(hCgs*nuArray(i)))*annu(n,dble(nuArray(i)))* &
 ((2.*real(dble(hCgs)*dble(nuArray(i))**3) / cSpeed**2) + jnu)*exp(-hCgs*nuArray(i)/(kErg*grid%temperature(i1,i2,i3)))
          Jnu = 4.*w*hnu(i+1)
          fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))* &
 ((2.*real(dble(hcgs)*dble(nuArray(i+1))**3) / cSpeed**2) + jnu)*exp(-hcgs*nuArray(i+1)/(kErg*grid%temperature(i1,i2,i3)))
          tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
       enddo

    end if ! (grid%adaptive)
    
    integral1 = tot
    
  end function integral1

  
  real function integral2(n,hnu, nuArray, nNu, rVec, grid, nStar)
    real, dimension(:) :: hnu, nuArray
    integer :: nNu, nStar, n
    type(VECTOR) :: rVec
    type(GRIDTYPE) :: grid
    real :: r
    integer :: i, iMin
    real :: fac1, fac2
    real :: w, x1, freq, tot, jnu

    if (nStar == 1) then
       r = modulus(rVec-grid%starPos1)
       x1 = sqrt(max(0.,(1. - grid%rStar1**2 / r**2)))
       w = 0.5*(1. - x1)
    else
       r = modulus(rVec-grid%starPos2)
       x1 = sqrt(max(0.,1. - grid%rStar2**2 / r**2))
       w = 0.5*(1. - x1)
    endif

    tot = 0.
    freq = ((13.598-eTrans(n))*1.602192e-12)/hcgs
    if ((freq < nuArray(1)).or.(freq > nuArray(nNu))) then
       write(*,*) "error in integral2",freq,nuArray(1),nuArray(nNu)
       iMin = 1
       jnu = 1.e-28
    else
       call hunt(nuArray, nNu, freq, iMin)
       jnu = 4.*w*logint(freq,nuArray(imin), nuArray(iMin+1), hnu(imin), hnu(imin+1))
    endif
    fac1  = (fourPi/(hcgs*freq))*annu(n,dble(freq))*jnu
    jnu = 4.*w*hnu(imin+1)
    fac2  = (fourPi/(hcgs*nuArray(imin+1)))*annu(n,dble(nuArray(imin+1)))*jnu
    tot = tot + 0.5*(fac1 + fac2)*(nuArray(imin+1)-freq)
    do i = iMin+1, nNu-1
       Jnu = 4.*w*hnu(i)
       fac1  = (fourPi/(hcgs*nuArray(i)))*annu(n,dble(nuArray(i)))*jnu
       Jnu = 4.*w*hnu(i+1)
       fac2  = (fourPi/(hcgs*nuArray(i+1)))*annu(n,dble(nuArray(i+1)))*jnu
       tot = tot + 0.5*(fac1 + fac2)*(nuArray(i+1)-nuArray(i))
    enddo
    integral2 = tot
  end function integral2


  subroutine setupMatrices(x, alpha, beta, np, rVec, i1, i2, i3, grid, &
       hnu1, nuArray1, nnu1, hnu2, nuarray2, nnu2, visFrac1, visFrac2, &
       isBinary, thisOctal, thisSubcell)
    integer, intent(in)                :: np
    type(GRIDTYPE), intent(inout)      :: grid
    logical, intent(in)                :: isBinary
    real, intent(in)                   :: visFrac1, visFrac2
    real(kind=doubleKind),intent(in),dimension(:) :: x
    real(kind=doubleKind), intent(out) :: alpha(:,:)
    real(kind=doubleKind), intent(out) :: beta(:)
    real, intent(in), dimension(:)     :: hnu1, hnu2
    real, intent(in), dimension(:)     :: nuarray1, nuarray2
    integer, intent(in)                :: nNu1, nNu2
    type(VECTOR), intent(in)           :: rVec
    integer, intent(in)                :: i1, i2, i3
    type(octal), pointer, optional     :: thisOctal 
    integer, intent(in),optional       :: thisSubcell 
    
    real(kind=doubleKind)              :: tmp, incr
    real(kind=doubleKind), parameter   :: fac = 1.e-3_db
    integer                            :: i, j
    integer, parameter                 :: maxLevels = statEqMaxLevels
    
!integer :: labels(8)
!labels = thisOctal%label
    if (grid%adaptive) then 

       thisOctal%N(thisSubcell,1:maxLevels) = x(1:maxLevels)
       thisOctal%Ne(thisSubcell) = x(maxLevels+1)
   
       do i = 1, maxLevels
          beta(i) = -equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, &
               rVec, i1, i2, i3, grid, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell)
       enddo
       beta(maxLevels+1) = -equation14(maxLevels, grid, i1, i2, i3, thisOctal, thisSubcell)
   
   !    write(*,*) "beta",real(beta(1:7))
       do i = 1, maxLevels
          do j = 1, maxLevels
             tmp = thisOctal%N(thisSubcell,j)
             thisOctal%N(thisSubcell,j) = thisOctal%N(thisSubcell,j) * (1.e0_db + fac)
             incr = equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, rVec, &
                  i1, i2, i3, grid, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell)
             alpha(i,j) = (beta(i) + incr) / (tmp*fac)
             if (alpha(i,j) == 0.) then
   !             write(*,*) i,j,beta(i),incr,tmp,grid%n(i1,i2,i3,j),grid%temperature(i1,i2,i3)
             endif
             thisOctal%N(thisSubcell,j) = tmp
          enddo
       enddo
   
   !    do j = 1, maxLevels
   !       tmp = grid%N(i1,i2,i3,j)
   !       grid%N(i1,i2,i3,j) = grid%N(i1,i2,i3,j) * (1.d0 + fac)
   !       incr = equation14(maxLevels, grid, i1, i2, i3)
   !       alpha(maxLevels+1,j) = 1.d0 !(beta(maxLevels+1) + incr) / (tmp*fac) ! 1.d0
   !       grid%N(i1,i2,i3,j) = tmp
   !    enddo
   
       alpha(maxLevels+1,1:maxLevels) = 1.e0_db
   
       tmp = thisOctal%Ne(thisSubcell)
       thisOctal%Ne(thisSubcell) = thisOctal%Ne(thisSubcell) * (1.e0_db+fac)
       do i = 1, maxLevels
          incr =  equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, rVec, &
               i1, i2, i3, grid, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell)
   !       alpha(i, maxLevels+1) = (beta(maxLevels+1)+incr)/(tmp*fac)
          alpha(i, maxLevels+1) = (beta(i)+incr)/(tmp*fac)
       enddo
   
       incr = equation14(maxLevels, grid, i1, i2, i3, thisOctal, thisSubcell)
       incr = incr + beta(maxLevels+1)
       alpha(maxLevels+1,maxLevels+1) = incr/(tmp*fac)
       thisOctal%Ne(thisSubcell) = tmp
   
   !    do i = 1, maxLevels+1
   !       do j = 1, maxLevels+1
   !          write(*,*) i,j,real(alpha(i,j))
   !       enddo
   !    enddo
   
    else ! grid not adaptive
           
       do i = 1, maxLevels
          grid%N(i1,i2,i3,i) = x(i)
       enddo
       grid%Ne(i1,i2,i3) = x(maxLevels+1)
   
       do i = 1, maxLevels
          beta(i) = -equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, &
               rVec, i1, i2, i3, grid, visFrac1, visFrac2, isBinary)
       enddo
       beta(maxLevels+1) = -equation14(maxLevels, grid, i1, i2, i3)
   
   !    write(*,*) "beta",real(beta(1:7))
       do i = 1, maxLevels
          do j = 1, maxLevels
             tmp = grid%N(i1,i2,i3,j)
             grid%N(i1,i2,i3,j) = grid%N(i1,i2,i3,j) * (1.e0_db + fac)
             incr = equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, rVec, &
                  i1, i2, i3, grid, visFrac1, visFrac2, isBinary)
             alpha(i,j) = (beta(i) + incr) / (tmp*fac)
             if (alpha(i,j) == 0.) then
   !             write(*,*) i,j,beta(i),incr,tmp,grid%n(i1,i2,i3,j),grid%temperature(i1,i2,i3)
             endif
             grid%N(i1,i2,i3,j) = tmp
          enddo
       enddo
   
   !    do j = 1, maxLevels
   !       tmp = grid%N(i1,i2,i3,j)
   !       grid%N(i1,i2,i3,j) = grid%N(i1,i2,i3,j) * (1.d0 + fac)
   !       incr = equation14(maxLevels, grid, i1, i2, i3)
   !       alpha(maxLevels+1,j) = 1.d0 !(beta(maxLevels+1) + incr) / (tmp*fac) ! 1.d0
   !       grid%N(i1,i2,i3,j) = tmp
   !    enddo
   
       alpha(maxLevels+1,1:maxLevels) = 1.e0_db
   
       tmp = grid%Ne(i1,i2,i3)
       grid%Ne(i1,i2,i3) = grid%Ne(i1,i2,i3) * (1.e0_db+fac)
       do i = 1, maxLevels
          incr =  equation8(i, maxLevels, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, rVec, &
               i1, i2, i3, grid, visFrac1, visFrac2, isBinary)
   !       alpha(i, maxLevels+1) = (beta(maxLevels+1)+incr)/(tmp*fac)
          alpha(i, maxLevels+1) = (beta(i)+incr)/(tmp*fac)
       enddo
   
       incr = equation14(maxLevels, grid, i1, i2, i3)
       incr = incr + beta(maxLevels+1)
       alpha(maxLevels+1,maxLevels+1) = incr/(tmp*fac)
       grid%Ne(i1,i2,i3) = tmp
   
   !    do i = 1, maxLevels+1
   !       do j = 1, maxLevels+1
   !          write(*,*) i,j,real(alpha(i,j))
   !       enddo
   !    enddo

    end if ! (grid%adaptive)
  
!if (SUM(labels-thisOctal%label) /= 0) print *, "        leaving setupMatrices, labels = ", thisOctal%label    

  end subroutine setupMatrices

  subroutine mnewt(grid, ntrial,x,n,tolx,tolf, Hnu1, nuArray1, nNu1, Hnu2, nuArray2, nNu2, &
                   rVec, i1, i2, i3, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell)
    implicit none

    integer, parameter             :: np=15
    type(GRIDTYPE), intent(inout)  :: grid
    integer, parameter             :: maxLevels = statEqMaxLevels
    integer, intent(in)            :: i1, i2, i3
    logical, intent(in)            :: isBinary
    real, intent(in)               :: visFrac1, visFrac2
    type(VECTOR), intent(in)       :: rVec
    real, intent(in), dimension(:) :: nuarray1, nuarray2
    real, intent(in), dimension(:) :: hnu1, hnu2
    integer, intent(in)            :: nNu1, nNu2
    real(kind=doubleKind),intent(inout),dimension(:) :: x
    integer, intent(in)            :: nTrial
    integer, intent(in)            :: n
    integer                        :: i, k
    integer                        :: indx(np)
    real(kind=doubleKind)          :: tolx, tolf
    real(kind=doubleKind)          :: errf, d, errx
    real(kind=doubleKind)          :: alpha(np,np),beta(np)
    type(octal), pointer, optional :: thisOctal 
    integer, intent(in),optional   :: thisSubcell 
      
!integer :: labels(8)
!labels = thisOctal%label
      do k = 1, ntrial
        call setupMatrices(x,alpha,beta,np,rVec, i1, i2, i3, grid,Hnu1, nuArray1, nNu1, Hnu2, &
             nuArray2, nNu2, visFrac1, visFrac2, isBinary, thisOctal, thisSubcell)
        errf=0.
        do i = 1, n
          errf=errf+abs(beta(i))
        end do 
!          write(*,*) "errf,tolf",errf,tolf
        if(errf.le.tolf)return

        call ludcmp(alpha,n,np,indx,d) 

        call lubksb(alpha,n,np,indx,beta)
        errx=0.
        do i = 1, n
          errx=errx+abs(beta(i))
          x(i)=x(i)+beta(i)
        end do 
!          write(*,*) "errx,tolx",errx,tolx
        if (errx.le.tolx) return
      end do
!if (SUM(labels-thisOctal%label) /= 0)  print *, "      leaving mnewt, labels = ", thisOctal%label    

  end subroutine 

      subroutine ludcmp(a,n,np,indx,d)
      real(kind=doubleKind), parameter :: tiny=1.0e-20_db
      integer, parameter :: nmax=200
      integer :: n, np
      real(kind=doubleKind) :: d
      real(kind=doubleKind) ::  a(np,np),vv(nmax)
      integer :: indx(n)
      integer :: i
      real(kind=doubleKind) :: aamax
      real(kind=doubleKind) :: sum, dum
      integer :: j, k
      integer :: imax
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
           write(*,*) 'singular matrix.'
           stop
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        if (j.gt.1) then
          do 14 i=1,j-1
            sum=a(i,j)
            if (i.gt.1)then
              do 13 k=1,i-1
                sum=sum-a(i,k)*a(k,j)
13            continue
              a(i,j)=sum
            endif
14        continue
        endif
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          if (j.gt.1)then
            do 15 k=1,j-1
              sum=sum-a(i,k)*a(k,j)
15          continue
            a(i,j)=sum
          endif
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(j.ne.n)then
          if(a(j,j).eq.0.)a(j,j)=tiny
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      if(a(n,n).eq.0.)a(n,n)=tiny
      end subroutine

      subroutine lubksb(a,n,np,indx,b2)
        integer :: np,n
        real(kind=doubleKind) ::  a(np,np)
        integer :: indx(n)
        real(kind=doubleKind) :: b2(n)
        integer :: ii, i, ll, j
        real(kind=doubleKind) :: sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b2(ll)
        b2(ll)=b2(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b2(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b2(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b2(i)
        if(i.lt.n)then
          do 13 j=i+1,n
            sum=sum-a(i,j)*b2(j)
13        continue
        endif
        b2(i)=sum/a(i,i)
14    continue
      end subroutine
   

  real(kind=doubleKind) pure function alpkk(freq,t)
     !
     ! this function returns the free-free absorption coefficient for hydrogen
     !
     real(kind=doubleKind), intent(in) :: freq,t
     real(kind=doubleKind)             :: wav,gauntf
      
     wav=1.e8_db*cSpeed/freq
     gauntf=giii(1.e0_db,t,wav)
     alpkk=gauntf*real(3.6d8/((dble(freq)**3)*sqrt(dble(t))))
     
  end function alpkk


  real(kind=doubleKind) pure function giii (z, t, wl)
     !
     !   ferland's fabulous functional fits
     !
     
     real(kind=doubleKind), intent(in) :: wl, t, z
     real(kind=doubleKind) :: c, u, ulog, gam2
     integer               :: i,j,k, m
     real(kind=doubleKind) :: b2
     real(kind=doubleKind) :: frac, sum1, sum2, d
     ! making coeff and a2 PARAMETERs may cause problems with XL Fortran
     real(kind=doubleKind) :: coeff(28) 
     real(kind=doubleKind) :: a2(7)
     coeff =                                                           &
        (/1.102d0       ,-0.1085d0     ,0.09775d0     ,-0.01125d0     ,&
          1.2d0         ,-0.24016667d0 ,0.07675d0     ,-0.01658333d0  ,&
          1.26d0        ,-0.313166667d0,0.15075d0     ,0.00241667d0   ,&
          1.29d0        ,-0.4518333d0  ,0.12925d0     ,0.00258333d0   ,&
          1.27d0        ,-0.579d0      ,0.092d0       ,-0.003d0       ,&
          1.16d0        ,-0.707333d0   ,0.112d0       ,0.0053333333d0 ,&
          0.883d0       ,-0.76885d0    ,0.190175d0    ,0.022675d0     /)
     a2 = (/100.d0, 10.d0, 3.d0, 1.d0, 0.3d0, 0.1d0, 0.001d0/)
       
       u = 1.44e+8 / (wl*t)
       ulog = log10(u)
       gam2 = 1.58e+5 * z*z/t
       if (gam2.gt.a2(7)) go to 10
         i = 7
         j = 7
         k = 7
         frac = 0.5
         go to 60
   10  continue
       if (gam2.lt.a2(1)) go to 20
         i = 1
         j = 1
         k = 1
         frac = 0.5
         go to 60
   20  continue
       do 30 i = 2, 7
           if (gam2.gt.a2(i)) go to 40
   30  continue
   40  continue
       k = i - 1
   50  continue
       b2 = log10(a2(k))
       c = log10(a2(i))
       gam2 = log10(gam2)
       frac = abs ((gam2-b2) / (b2-c))
   60  continue
       k = (k-1)*4
       sum1 = coeff(k+1)
       d = 1.0
       do 70 m = 2, 4
           d = d*ulog
           sum1 = sum1 + coeff(k+m)*d
   70  continue
       sum1 = sum1 * (1.0 - frac)
       i = (i-1)*4
       sum2 = coeff(i+1)
       d = 1.0
       do 80 m = 2, 4
           d = d*ulog
           sum2 = sum2 + coeff(i+m)*d
   80  continue
       sum2 = sum2 * frac

       giii = sum1 + sum2
  end function giii


    subroutine occultTest(grid, i1, i2, i3, starPos, starRadius, &
         occultPos, occultRadius, visFrac)

      type(GRIDTYPE) :: grid
      integer :: i1, i2, i3
      real :: visFrac
      type(VECTOR) :: toStar, toOccult, direction, rVec, starPos, occultPos
      real :: distToStar, distToOccult, h
      real :: starRadius, occultRadius
      real :: thetaToStar, phiToStar
      logical :: occulted
      real :: dotProd, sinAng
      real :: r, occultedFrac
      integer :: nTheta = 10, nPhi = 10
      real :: cosAng, ang, theta, phi, dtheta, dphi, domega
      integer :: i, j

      rVec = VECTOR(grid%xAxis(i1), grid%yAxis(i2), grid%zAxis(i3))
      toStar = starPos - rVec
      toOccult = occultPos - rVec
      distToStar = modulus(toStar)
      distToOccult = modulus(toOccult)
      toStar = toStar / distToStar
      toOccult = toOccult / distToOccult
      
      if (distToOccult > distToStar) then
         visFrac = 1.
         goto 666
      endif

      
      visFrac = 0.
      occultedFrac = 0.

      h  = sqrt(max(0.,distToStar**2 - starRadius**2))
      cosang = h / distToStar
      ang = acos(min(1.,max(-1.,cosAng)))
      call getPolar(toStar, r, thetaTostar, phiToStar)
      dtheta = pi / real(ntheta-1)
      dphi = twopi / real(nphi-1)
      do i = 1, ntheta
       theta = thetaToStar + (2.*real(i-1)/real(nTheta-1)-1.)*ang
       if (theta > pi) theta = theta - pi
       if (theta < 0.) theta = theta + pi
       dTheta = 2.*ang/real(nTheta-1)
       do j = 1, nphi
          phi = phiToStar + (2.*real(j-1)/real(nPhi-1)-1.)*ang
          if (phi < 0.) phi = phi + twoPi
          if (phi > twoPi) phi = phi + twoPi
          dphi = 2.*ang/real(nPhi-1)

          occulted = .false.
          direction = vector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
          domega = sin(theta)*dtheta*dphi
          dotprod = direction .dot. tostar
          if ((dotprod > 0.) .and. (acos(min(1.,max(-1.,dotprod))) < ang)) then

             visFrac = visFrac + dOmega
             dotprod = direction .dot. tooccult
             if (dotprod > 0.) then
                sinang = sqrt(max(0.,(1.-dotprod*dotprod)))
                if ((sinang*disttooccult) < occultradius) then 
                   occultedFrac = occultedFrac + dOmega
                endif
             endif
          endif
       enddo
    enddo
    visFrac = (visFrac-occultedFrac)/visFrac
666 continue
  end subroutine occultTest

  
  subroutine generateOpacities(grid, m, n)

    type(GRIDTYPE) :: grid
    integer, parameter :: maxLevels = statEqMaxLevels
    integer :: m, n
    integer :: i1, i2, i3
    real :: chil, fac
    real :: transe, thresh
    integer :: i,j,k
    real :: chi, eta
    real(kind=doubleKind) :: freq

    write(*,'(a)') "Generating opacities..."

    do k = 2, maxlevels
       do i=1, k-1
          lambdaTrans(i, k) = lambdaTrans(k, i)
          ! aEinstein(k, i) = aEinstein(k, i) commented out by nhs
          freq = cspeed / lambdaTrans(i, k)
          bEinstein(k, i) = ((dble(aEinstein(k,i))*dble(cspeed)**2) / (2.d0*dble(hcgs)*dble(freq)**3))
          bEinstein(i, k) = bEinstein(k,i) * gDegen(k) / gDegen(i)
          fStrength(k, i) = -gDegen(i) * fStrength(i,k) / gDegen(k)
       enddo
    enddo

    transe = abs(eTrans(n)-eTrans(m))
    freq = cspeed/lambdaTrans(m,n)

    grid%etaLine = 1.e-20
    grid%etaCont = 1.e-20
    grid%kappaAbs = 1.e-20
    grid%kappaSca = 1.e-20
    grid%chiLine = 1.e-20

    do i1 = 1, grid%na1
       do i2 = 1, grid%na2
          do i3 = 1, grid%na3
             if (.not.grid%inStar(i1,i2,i3).and.grid%inUse(i1,i2,i3)) then
                grid%kappasca(i1,i2,i3,1) = grid%ne(i1,i2,i3) * sigmaE * 1.e10
                !
                ! calculate the line opacity and emissivity
                !
                chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength(m,n)
                chil = chil* (grid%n(i1,i2,i3,m)-( ( gDegen(m) / gDegen(n)) * grid%n(i1,i2,i3,n)))
                grid%chiLine(i1,i2,i3) = 1.e10*chil
                
                if (grid%n(i1,i2,i3,n) == 0.d0) then
                   write(*,*) i1,i2,i3,n
                   write(*,*) grid%n(i1,i2,i3,1:maxLevels)
                   write(*,*) grid%Ne(i1,i2,i3)
                   write(*,*) grid%temperature(i1,i2,i3)
                   stop
                endif
                fac=( ((grid%n(i1,i2,i3,m)*gDegen(n))/(grid%n(i1,i2,i3,n)*gDegen(m)))-1.d0)
                grid%etaLine(i1,i2,i3)=1.e10*chil*real((2.d0*dble(hcgs)*dble(freq)**3)/(dble(cSpeed**2)))/fac

!                write(*,*) i1,i2,i3,grid%etaline(i1,i2,i3)/grid%chiline(i1,i2,i3)
                !
                ! continuous opacity.. bound-free and free-free processes (+es)
                !
                chi=0.d0
                do j=1,maxLevels
                   thresh=(13.598d0-eTrans(j))
                   if (transe.ge.thresh) then
                      chi=chi+(grid%n(i1,i2,i3,j)- &
                           boltzsaha(j, grid%ne(i1,i2,i3), dble(grid%temperature(i1,i2,i3)))* &
                           exp((-hcgs*freq)/(kerg*grid%temperature(i1,i2,i3))))* annu(j,dble(freq))
                   endif
                enddo
                
                chi=chi+real(grid%ne(i1,i2,i3))**2*alpkk(freq,dble(grid%temperature(i1,i2,i3)))*&
                     (1.d0-exp((-hcgs*freq)/(kerg*grid%temperature(i1,i2,i3))))
!                chi=chi+grid%ne(i1,i2,i3)*sigmaE
                
                grid%kappaabs(i1,i2,i3,1) = chi * 1.e10
                !
                ! continuous emissivity...bf and ff
                ! 
                eta=0.d0
                do j=1,15
                   thresh=(13.598d0-eTrans(j))
                   if (transe.ge.thresh) then
                      eta=eta+boltzsaha(j, grid%ne(i1,i2,i3),dble(grid%temperature(i1,i2,i3))) &
                           *annu(j,freq)*exp(-(hcgs*freq)/(kerg*grid%temperature(i1,i2,i3)))
                   endif
                enddo
                
                eta=eta + (grid%ne(i1,i2,i3)**2) * alpkk(freq,dble(grid%temperature(i1,i2,i3)))* &
                     exp(-(hcgs*freq)/(kerg*grid%temperature(i1,i2,i3)))
                
                eta=eta*real((2.0*dble(hcgs)*dble(freq)**3)/(dble(cspeed)**2))
                
                grid%etacont(i1,i2,i3) = eta*1.e10

                if (grid%chiLine(i1,i2,i3) < 0.) then
                   write(*,*) i1, i2, i3, chil
                   grid%chiLine(i1,i2,i3) = 1.e-20
                   grid%etaLine(i1,i2,i3) = 1.e-20
                   grid%etaCont(i1,i2,i3) = 1.e-20
                   grid%kappaAbs(i1,i2,i3,1) = 1.e-20
                   grid%kappaSca(i1,i2,i3,1) = 1.e-20
                endif

             endif
          enddo
       enddo
    enddo
    write(*,'(a)') "Done."


    where (grid%kappaabs < 0.) grid%kappaabs = 1.e-20
  end subroutine generateOpacities

  subroutine amrStateq(grid, contfile, lte, nLower, nUpper,  ion_name, ion_frac)
    ! calculate the statistical equilibrium for all of the subcells in an
    !   adaptive octal grid.

    implicit none
    
    type(GRIDTYPE),intent(inout):: grid      
    character(len=*),intent(in) :: contfile      ! filename for continuum flux
    logical,intent(in)          :: lte           ! true if lte conditions
    integer,intent(in)          :: nLower, nUpper! level populations
    ! Name of the ion (see opacity_lte_mod.f90 for the list of a valid name.)
    character(LEN=*),intent(in), optional :: ion_name
    real,intent(in), optional             :: ion_frac      ! n_ion/n_specie
    
    integer, parameter          :: maxLevels = statEqMaxLevels ! number of levels to compute
    integer                     :: iOctal        ! loop counter
    integer                     :: iSubcell      ! loop counter
    integer                     :: nOctal        ! number of octals in grid
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    real(kind=doublekind),parameter :: tolx = 1.e-3  ! tolerance
    real(kind=doublekind),parameter :: tolf = 1.e-3  ! tolerance
    integer                     :: returnVal     ! status value
    real                        :: chil
    real                        :: chi
    real                        :: fac
    real                        :: eta
    real                        :: etal
    real                        :: chi_es
    real                        :: thresh
    real(kind=doubleKind),allocatable,dimension(:) :: xall
    real(kind=doublekind)       :: nTot
    real(kind=doubleKind)       :: phiT
    real                        :: visFrac1      ! visible fraction of star 1
    real                        :: visFrac2      ! visible fraction of star 2
    logical                     :: isBinary      ! true for binary system
    integer                     :: nNu1, nNu2
    real                        :: hNu1(2000)
    real                        :: hNu2(2000)
    real                        :: nuarray1(2000)
    real                        :: nuarray2(2000)
    integer                     :: i, j, m       ! loop counters
    real                        :: transe
    real(kind=doubleKind)       :: freq
    type(vector)                :: rVec
    real (kind=doubleKind)      :: Ne1, Ne2
    logical                     :: ok
    type(octal), pointer        :: thisOctal => null()
    real(kind=doubleKind)       :: tmp0, tmp1, tmp2   ! used for debug
    real                        :: departCoeff(maxLevels)
    logical                     :: debugInfo = .true.

    ! Initialize the data arrays (lambdaTrans, bEinstein, fStrength) defined at the top of this module.
    call map_transition_arrays(maxLevels)
    
    transe = abs(eTrans(nUpper)-eTrans(nLower))
    freq = cspeed/lambdaTrans(nLower,nUpper)
    
    open(20,file=contfile,status="old",form="formatted")
    nnu1 = 1
    do
       read(20,*,iostat=returnVal) nuarray1(nNu1), hnu1(nNu1)
       if (returnVal /= 0) exit
       nNu1 = nNu1 + 1
    end do
    nNu1 = nNu1  - 1
    close(20)
    hNu1(1:nNu1) = hnu1(1:nNu1) / fourPi   ! Converts from flux to flux momemnt
print *, contfile, 'nnu1 = ',nnu1,'  nuarray1(nnu1) = ',nuarray1(nnu1)

    allocate(octalArray(grid%nOctals))
    
    ! get an array of octals comprising the entire tree
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)

    ! now loop over all octals
    print *, "   calculating LTE values..." 
    do iOctal = 1, nOctal
       thisOctal => octalArray(iOctal)%content
       do iSubcell = 1, 8
          if (octalArray(iOctal)%inUse(iSubcell)) then
             nTot = dble(thisOctal%rho(iSubcell) / mHydrogen)
             
             phiT = 1.5e0_db*log10(2.e0_db*pi*mElectron*kErg*thisOctal%temperature(iSubcell)) - &
                      3.e0_db*log10(hCgs) + &
                      log10(exp(-13.6e0_db/(kEV*real(thisOctal%temperature(iSubcell),kind=db))))
             phiT = 10.e0_db**phiT             
             
             call solveQuadDble(1.e0_db, phiT, -1.e0_db*phiT*real(nTot,kind=doublekind), Ne1, Ne2, ok)
             thisOctal%Ne(iSubcell) = min(max(Ne1,Ne2),nTot)
             thisOctal%nTot(iSubcell) = nTot 
             do m = 1, maxlevels
                thisOctal%n(iSubcell,m) =               &
                   boltzsaha(m,thisOctal%ne(iSubcell),real(thisOctal%temperature(iSubcell),kind=db))
             enddo
          endif
       enddo

    enddo
    print *, "   ...LTE calculations complete." 


    if (.not.lte) then
       print *, "   calculating non-LTE values..." 
            
       allocate(xall(1:maxLevels+1))
       visFrac1 = 1.
       visFrac2 = 0.
       isBinary = .false.
       
       do iOctal = 1, nOctal
       if (debugInfo) print *, 'Octal #',iOctal
          thisOctal => octalArray(iOctal)%content
          do iSubcell = 1, 8, 1
             if (octalArray(iOctal)%inUse(iSubcell) .and. &
                (thisOctal%nTot(iSubcell) > 1.0)) then

                xAll(1:maxLevels) = thisOctal%N(iSubcell,1:maxLevels)
                xAll(maxLevels+1) = thisOctal%Ne(iSubcell)
                rVec = subcellCentre(thisOctal,iSubcell)

                call mnewt(grid, 20, xAll, maxlevels+1, tolx, tolf, hNu1, nuArray1, nNu1, &
                           hNu2, nuArray2, nNu2, rVec, 1, 1, 1, visFrac1, visFrac2,&
                           isBinary, thisOctal, iSubcell)
                           
                if (debugInfo) then
                   write (*,'(a12,i1,a15,f5.0,a11,e10.1)') '   subcell #', iSubcell, '  temperature: ', & 
                                      thisOctal%temperature(iSubcell),'  density: ',thisOctal%rho(iSubcell) 
                   do i = 1 , maxLevels
                      departCoeff(i) = real(xall(i))/boltzSaha(i, thisOctal%Ne(iSubcell),          &
                                                         real(thisOctal%temperature(iSubcell),kind=db))
                      write(*,'(a5,i3,1p,e12.3,e12.3)') '     ',i,departCoeff(i),xall(i)
                   enddo
                end if   
                
             endif
          enddo
       enddo
       
       deallocate(xall)
       print *, "   non-LTE calculations complete..." 
       
    endif

    write(*,'(a,f8.1)') "Generating opacities for ",lambdaTrans(nLower, nUpper)*1.e8

    do iOctal = 1, nOctal
       do iSubcell = 1, 8
          if (octalArray(iOctal)%inUse(iSubcell)) then
                  
             octalArray(iOctal)%content%kappaSca(iSubcell,1) = &
                octalArray(iOctal)%content%Ne(iSubcell) * sigmae * 1.e10
             

             ! calculate the line opacity and emissivity
             
             chil=( (pi*eCharge**2) / (mElectron*cSpeed) ) * fStrength(nLower,nUpper)
             chil = chil * (octalArray(iOctal)%content%N(iSubcell,nLower) - &
                       ((gDegen(nLower) / gDegen(nUpper)) * octalArray(iOctal)%content%N(iSubcell,nUpper)))
             octalArray(iOctal)%content%chiLine(iSubcell) = 1.e10 * chil
             
             if (octalArray(iOctal)%content%n(iSubcell,nUpper) == 0.e0_db) then
                write(*,*) 'In amrStatEq, octalArray(iOctal)%content%n(iSubcell,nUpper) == 0.d0'
                write(*,*) nUpper
                write(*,*) octalArray(iOctal)%content%N(iSubcell,1:maxLevels)
                write(*,*) octalArray(iOctal)%content%Ne(iSubcell)
                write(*,*) octalArray(iOctal)%content%temperature(iSubcell)
                stop
             endif
             
             fac=(((octalArray(iOctal)%content%N(iSubcell,nLower)* gDegen(nUpper)) / &
                 (octalArray(iOctal)%content%N(iSubcell,nUpper)*gDegen(nLower)))-1.e0_db)
             octalArray(iOctal)%content%etaLine(iSubcell) = &
                  1.e10*chil*real((2.e0_db*dble(hcgs)*dble(freq)**3)/(dble(cSpeed**2)))/fac
             
             
             ! continuous opacity.. bound-free and free-free processes (+es)
             
             chi=0.e0_db
             do j=1,maxLevels
                thresh=(13.598e0_db-eTrans(j))
                if (transe.ge.thresh) then
                   chi=chi+(octalArray(iOctal)%content%N(iSubcell,j)- &
                        boltzsaha(j, real(octalArray(iOctal)%content%Ne(iSubcell),kind=db), &
                                  real(octalArray(iOctal)%content%temperature(iSubcell),kind=db))* &
                        exp((-hcgs*freq)/(kerg*octalArray(iOctal)%content%temperature(iSubcell)))) * &
                        annu(j,real(freq,kind=db))
                endif
             enddo
             
             chi=chi+real(octalArray(iOctal)%content%Ne(iSubcell))**2 * &
                    alpkk(freq,real(octalArray(iOctal)%content%temperature(iSubcell),kind=db))*&
                    (1.e0_db-exp((-hcgs*freq)/(kerg*octalArray(iOctal)%content%temperature(iSubcell))))
             !                chi=chi+grid%ne(i1,i2,i3)*sigmaE
             
             octalArray(iOctal)%content%kappaAbs(iSubcell,1) = chi * 1.e10

             
             ! continuous emissivity...bf and ff
              
             eta=0.e0_db
             do j=1,15
                thresh=(13.598e0_db-eTrans(j))
                if (transe.ge.thresh) then
                   eta=eta+boltzsaha(j, real(octalArray(iOctal)%content%Ne(iSubcell),kind=db),&
                                     real(octalArray(iOctal)%content%temperature(iSubcell),kind=db)) * &
                        annu(j,freq)*exp(-(hcgs*freq)/(kerg*octalArray(iOctal)%content%temperature(iSubcell)))
                endif
             enddo
             
             eta=eta + (octalArray(iOctal)%content%Ne(iSubcell)**2) * &
                        alpkk(freq,real(octalArray(iOctal)%content%temperature(iSubcell),kind=db))* &
                        exp(-(hcgs*freq)/(kerg*octalArray(iOctal)%content%temperature(iSubcell)))
             
             eta=eta*real((2.0*dble(hcgs)*dble(freq)**3)/(dble(cspeed)**2))
             
             octalArray(iOctal)%content%etaCont(iSubcell) = eta*1.e10
             
             if (octalArray(iOctal)%content%chiLine(iSubcell) < 0.) then
                write(*,*) iOCtal, chil
                octalArray(iOctal)%content%chiLine(iSubcell)    = 1.e-20
                octalArray(iOctal)%content%etaLine(iSubcell)    = 1.e-20
                octalArray(iOctal)%content%etaCont(iSubcell)    = 1.e-20
                octalArray(iOctal)%content%kappaAbs(iSubcell,1) = 1.e-20
                octalArray(iOctal)%content%kappaSca(iSubcell,1) = 1.e-20
             endif

          endif
       enddo
    enddo

    deallocate(octalArray)

  end subroutine amrStateq


  ! Initialize the data arrays (lambdaTrans, bEinstein, fStrength) defined at the top of this module.
  subroutine map_transition_arrays(maxLevels)
    implicit none
    integer, intent(in)   :: maxLevels      
    integer               :: i, k
    real(kind=doubleKind) :: freq    

        do k = 2, maxlevels
           do i=1, k-1
              lambdaTrans(i, k) = lambdaTrans(k, i)
              ! aEinstein(k, i) = aEinstein(k, i) commented out by nhs
              freq = cspeed / lambdaTrans(i, k)
              bEinstein(k, i) = ((dble(aEinstein(k,i))*dble(cspeed)**2)  &
                   &   / (2.e0_db*dble(hcgs)*dble(freq)**3))
              bEinstein(i, k) = bEinstein(k,i) * gDegen(k) / gDegen(i)
              fStrength(k, i) = -gDegen(i) * fStrength(i,k) / gDegen(k)
           end do
        end do

      end subroutine map_transition_arrays

end module stateq_mod






