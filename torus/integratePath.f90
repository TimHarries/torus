module path_integral

  use vector_mod             ! vector maths 
  use grid_mod               ! opacity grid
  use math_mod               ! misc maths routines
  use photon_mod             ! the photon module
  use constants_mod          ! physical constants
  implicit none
!
! This subroutine performs an integration through the grid in order
! to compute the scattering and absorption optical depths.
!
  contains

subroutine integratePath(wavelength,  lambda0, vVec, aVec, uHat, Grid,  lambda, &
     tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton, &
     lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, redRegion, rStar, &
     coolStarPosition, iRev, usePops, mLevel, nLevel, fStrength, gM, gN, localTau, interp)


  type(GRIDTYPE) :: grid                            ! the opacity grid

  real :: localTau
  logical :: usePops
  integer :: mLevel, nLevel
  real :: fStrength, gM, GN
  real :: iRev
  real :: rStar
  logical :: interp
  type(VECTOR) :: aVec                              ! starting position vector
  type(VECTOR) :: uHat                              ! direction
  type(VECTOR) :: rVec                              ! position vector
  type(VECTOR) :: rHat                              ! unit radius vector
  type(VECTOR) :: vVec                              ! velocity vector
  type(VECTOR) :: fVec
  type(VECTOR) :: coolStarPosition                  ! pos of cool giant star
  type(VECTOR) :: rFromStar

  logical :: redRegion                              ! use raman scattering red region opacities

  real :: sinTheta, dr, dtheta, dphi

  logical :: thinLine                               ! ignore line absorption of continuum

  real :: lambda0                                   ! rest wavelength of line

  real :: escProb                                   ! the escape probability

  logical :: hitDisk
  real :: distToDisk

  real :: t1, t2, t3                                ! interpolation factors

  real, parameter :: precision = 1.e-4              ! relative precision

  integer :: maxTau

  real :: wavelength                                ! the wavelength and frequency
  real :: nu

  logical :: contPhoton                             ! is this a continuum photon?

  real ::  mu, cosTheta                    ! polar coords
  real :: r,  phi  
  real :: theta

  integer :: iStart                                 ! starting index

  real :: tauSob, tauSob1, tauSob2                  ! Sobolev optical depths

  real :: thisVel
  integer :: ilambda                                ! wavelength index




  integer, parameter :: maxLambda = 200              ! max size of tau arrays
  integer :: nLambda
  real :: lambda(1:maxTau),dlambda                         ! distance arrays
  real :: tauExt(1:maxTau), tauAbs(1:maxTau), tauSca(1:maxTau)           ! optical depths
  real :: tauCont(maxTau, maxLambda)


  real :: fudgeFactor = 1.00001                     ! fudge factor
  real :: projVel(maxTau)
  integer :: posIndex(maxTau,3)                     ! position indices in grid
  real :: posOffset(maxtau,3)                       ! position interpolation offsets
  real :: kabs(maxtau), ksca(maxtau), dlam(maxtau)  ! opacities
  type(VECTOR) :: rVecArray(maxTau)                 ! array of position vectors
  real :: x,x1,x2                                 ! multipliers

  logical :: endLoop                                ! is the loop finished?

  integer :: nTau                                   ! size of optical depth arrays
  real :: chil                                      ! line opacity
  real :: dx, dy, dz                                ! cartesian increments
  real :: lamStart, lamEnd
  real :: r1, r2
  
  real :: lambda2, thisVel2

  integer :: i, i1, i2, i3,  j                   ! counters
  logical :: opaqueCore                             ! is the core opaque
  logical :: hitcore                                ! has the photon hit the core
  logical :: ok                                     ! is everything ok?
  real    :: resFac                                 ! resolution of the integration
  logical :: outwards
  real :: kscaInterp, kabsInterp
  real :: minDist
  type(VECTOR) :: tVec

  ! initialize variables

  hitcore = .false.
  rVec = aVec
  escProb = 1.

  ! locate this wavelength in the grid of wavelengths

  if (grid%flatspec.or.(grid%doRaman)) then
     iLambda = 1
  else
     call hunt(grid%lamArray, grid%nLambda, wavelength, iLambda)
  endif


  ! There are two types of integration - one for cartesian grids (n.b. this needs updating) and one for polars

  if (grid%cartesian.and.(.not.grid%lineEmission)) then


     ! size of the increments

     dx = grid%xAxis(2)-grid%xAxis(1)
     dy = grid%yAxis(2)-grid%yAxis(1)
     dz = grid%zAxis(2)-grid%zAxis(1)

     ! size of the distance increment - this won't work will for non-uniform grids

     dLambda = min(dx,dy,dz) 

     lambda(1:1000) = 0.
     tauExt(1:1000) = 0.
     tauAbs(1:1000) = 0.
     tauSca(1:1000) = 0.

     nTau = 1

     ! first optical depths are all zero of course

     lambda(nTau) = 0.
     tauExt(nTau) = 0.
     tauAbs(nTau) = 0.
     tauSca(nTau) = 0.


     ! add on a small distance

     rVec = rVec + dlambda * uHat

     if ((modulus(rVec) /= 0.).and.(modulus(rVec-coolStarPosition) < rStar)) then
        rHat = rVec
        call normalize(rHat)
        outwards = .false.
        if ((rHat .dot. uHat) > 0.) outwards = .true.
     endif


     ! while inside the grid

     do  while (.not.outsideGrid(rVec, grid))

        call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)


        ! some grid points are empty

        nTau = nTau + 1

        ! if the overspecified array is still too small...

        if (nTau > maxTau) then
           write(*,*) "! tau array full"
           write(*,*) "nTau",nTau
           write(*,*) "rVec",rVec
           write(*,*) grid%xAxis(grid%nx),grid%yAxis(grid%ny),grid%zAxis(grid%nz)
           write(*,*) "dlambda",dlambda
           write(*,*) "red region",redRegion
           nTau = maxTau
           EXIT
        endif

        ! add on this small optical depth




        if (.not.redRegion) then


           if (interp) then
              kabsInterp = interpGridKappaAbs(grid,i1,i2,i3,iLambda,t1,t2,t3) 
              kscaInterp = interpGridKappaSca(grid,i1,i2,i3,iLambda,t1,t2,t3) 
           else
              kabsInterp = grid%kappaAbs(i1,i2,i3,iLambda)
              kscaInterp = grid%kappaSca(i1,i2,i3,iLambda)
           endif
           tauExt(nTau) = tauExt(nTau-1) + &
                (kabsInterp + kscaInterp)* dlambda
           tauSca(nTau) = tauSca(nTau-1) + kscaInterp * dlambda
           tauAbs(nTau) = tauAbs(nTau-1) + kabsInterp * dlambda
           lambda(nTau) = lambda(nTau-1) + dlambda
        else
           kabsInterp = interpGridkappaAbsRed(grid,i1,i2,i3,iLambda,t1,t2,t3) 
           kscaInterp = interpGridkappaScaRed(grid,i1,i2,i3,iLambda,t1,t2,t3) 

           tauSca(nTau) = tauSca(nTau-1) + kscaInterp * dlambda
           tauAbs(nTau) = tauAbs(nTau-1) + kabsInterp * dlambda
           tauExt(nTau) = tauAbs(ntau) + tauSca(ntau)
           lambda(nTau) = lambda(nTau-1) + dlambda
        endif
        rVec = rVec + dlambda *  uHat ! add on small distance

     enddo ! until out of the grid


  else

!!!!!!!! picky piecemeal integration - slow but accurate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     ! two steps for each radial step size - good enough. Higher resFac means
     ! a slower code but greater accuracy.

     resFac = 2.
     Ntau = 1



     ! set opacities and optical depths to zero

     ksca(nTau) = 0.
     kabs(nTau) = 0.
     dlam(nTau) = 0.

     lambda(nTau) = 0.
     tauExt(nTau) = 0.
     tauAbs(nTau) = 0.
     tauSca(nTau) = 0.

     ! some geometries produce photons at the centre of the grid

     if (modulus(rVec) == 0.)  then
        if (grid%rCore == 0.) then
           rVec = 1.0001 * grid%rAxis(1)*uHat
        else
           rVec = 1.0001 * grid%rCore*uHat
        endif
     endif


!     if (.not.grid%cartesian) then
!        r = modulus(rVec)
!        if (r < grid%rAxis(1)) then
!           rHat = (-1.)*(rVec / r)
!           cosTheta  = uHat .dot. rHat
!           call solveQuad(1.,-2.*r*cosTheta,r*r-grid%rAxis(1)*grid%rAxis(1),x1,x2,ok)
!           dlambda = max(x1,x2) * fudgeFactor
!           rVec = rVec +  dlambda * uHat   
!           lambda(nTau) = dlambda
!        endif
!     endif



     ! find the initial position of the position vector


     call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)


     ! find the project velocity

     thisVel = (lambda0 - wavelength)/lambda0
     thisVel = thisVel  + (uHat .dot. vVec)  

     if (grid%resonanceLine) then
        lambda2 = grid%lambda2
        thisVel2 = (lambda2 - wavelength)/lambda2
        thisVel2 = thisVel2  + (uHat .dot. vVec) 
     endif


     escProb = 1.

     ! opacities at first grid point


     if (interp) then
        kabs(nTau) = interpGridKappaAbs(grid,i1,i2,i3,iLambda,t1,t2,t3) 
        ksca(nTau) = interpGridKappaSca(grid,i1,i2,i3,iLambda,t1,t2,t3) 
     else
        kabs(nTau) = grid%kappaAbs(i1,i2,i3,iLambda)
        ksca(nTau) = grid%kappaSca(i1,i2,i3,iLambda)
     endif

     ! first position held in indices


     posIndex(nTau,1) = i1
     posIndex(nTau,2) = i2
     posIndex(nTau,3) = i3

     posOffset(nTau,1) = t1
     posOffset(nTau,2) = t2
     posOffset(nTau,3) = t3

     ! first position held as vector

     rVecArray(nTau) = rVec



     ! compute the escape probability if this is a line photon

     if (.not.contPhoton) then


        nu = cSpeed / (lambda0*angstromtocm)

        if (.not. usePops) then
           chil = interpGridChil(grid,i1,i2,i3,t1,t2,t3)
        else
           chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
           chil = 1.e10*chil* (grid%n(i1,i2,i3,mLevel)-( ( gM / gN) * grid%n(i1,i2,i3,nLevel)))
        endif


        tauSob = chil  / nu
        tauSob = tauSob / directionalDeriv(grid,rVec,i1,i2,i3,uHat)

        if (tauSob < 0.1) then
          escProb = 1.0d0-tauSob*0.5d0*(1.0d0 - tauSob/3.0d0*(1. - tauSob*0.25d0*(1.0d0 - 0.20d0*tauSob)))
        else if (tauSob < 15.) then
          escProb = (1.0d0-exp(-tauSob))/tauSob
        else
          escProb = 1.d0/tauSob
        end if
        localTau = tauSob
     endif


     ! first projected velocity 

     if (contPhoton) then
        projVel(1) = interpGridVelocity(grid,i1,i2,i3,t1,t2,t3) .dot. uHat
     else
        projVel(1) = thisVel
     endif


     ! first length step

     if (grid%cartesian) then
        if (i1 < grid%nx) then
           dx = grid%xAxis(i1+1)-grid%xAxis(i1)
        else
           dx = grid%xAxis(i1)-grid%xAxis(i1-1)
        endif
        if (i2 < grid%ny) then
           dy = grid%yAxis(i2+1)-grid%yAxis(i2)
        else
           dy = grid%yAxis(i2)-grid%yAxis(i2-1)
        endif
        if (i3 < grid%nz) then
           dz = grid%zAxis(i3+1)-grid%zAxis(i3)
        else
           dz = grid%zAxis(i3)-grid%zAxis(i3-1)
        endif

        fVec = VECTOR(dx, dy, dz)
        dlambda = abs(fvec%x*uHat%x) + abs(fVec%y*uHat%y) + abs(fVec%z*uHat%z)

     else
	r1 = grid%rAxis(i1)
	if (i1 < grid%nr) then
	   r2 = grid%rAxis(i1+1)
	else
	   r2 = grid%rAxis(i1-1)
	endif
	dr = abs(r2 - r1)
	sinTheta = sqrt(1.-grid%muAxis(i2)**2)
	if (i3 < grid%nphi) then
	   dx = r1 * sinTheta * cos(grid%phiAxis(i3)) - r2 * sinTheta * cos(grid%phiAxis(i3+1))
	   dy = r1 * sinTheta * sin(grid%phiAxis(i3)) - r2 * sinTheta * sin(grid%phiAxis(i3+1))
	else
	   dx = r1 * sinTheta * cos(grid%phiAxis(i3)) - r2 * sinTheta * cos(grid%phiAxis(i3-1))
	   dy = r1 * sinTheta * sin(grid%phiAxis(i3)) - r2 * sinTheta * sin(grid%phiAxis(i3-1))
	endif
	
	if (i2 < grid%nMu) then
	   dz = r1 * grid%muAxis(i2) - r2 * grid%muAxis(i2+1)
	else
	   dz = r1 * grid%muAxis(i2) - r2 * grid%muAxis(i2-1)
	endif
	fVec = VECTOR(dx, dy, dz)
	dlambda  = min(modulus(fVec),dr)
!	dlambda = abs(fvec%x*uHat%x) + abs(fVec%y*uHat%y) + abs(fVec%z*uHat%z)

     endif



     dlambda = dlambda / resFac


     endLoop = .false.

     hitDisk = .false.


     if (grid%geometry == "ttauri") then
        tVec = intersectionLinePlane(rVec, uHat, grid%diskNormal, 0., ok)
        if (ok) then
           minDist = modulus(tVec)
           if (minDist > grid%diskRadius) then
              hitCore = .true.
              hitDisk = .true.
              distToDisk = modulus(tVec-rVec)
           endif
        endif
     endif

     ! this is the big loop




     do while(.not.endLoop)



        ! add on the small distance step to the position

        rVec = rVec + dlambda * uHat
        call getPolar(rVec, r, theta, phi)
        mu = rVec%z/r

        ! if the photon is inside the inner boundary then it has hit the core

        ! if the core is opaque then the loop stops here, otherwise the
        ! photon cross the core and appears on the other side

        if (grid%geometry /= "binary") then
           if ((r < grid%rCore) .and. opaqueCore) then
              nTau = nTau + 1
              ksca(nTau) = 0.
              kabs(nTau) = 1.e20
              lambda(nTau) = lambda(nTau-1)
              hitCore = .true.
              endLoop = .true.
              posIndex(nTau,1:3) = posIndex(nTau-1,1:3)
              posOffset(nTau,1:3) = posOffset(nTau-1,1:3)
              rvecarray(ntau) = rVecArray(ntau-1)
              dlam(nTau) = dLam(nTau-1)
           else if ((r < grid%rCore).and.(.not.opaqueCore)) then
              rHat = (-1.)*(rVec / r)
              cosTheta  = uHat .dot. rHat
              call solveQuad(1.,-2.*r*cosTheta,r*r-grid%rCore*grid%rCore,x1,x2,ok)
              dlambda = max(x1,x2) * fudgeFactor
              rVec = rVec +  dlambda * uHat   
              r = modulus(rVec)
              mu = rVec%z/r
              call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)
           endif
        else
           rFromStar = rVec - grid%starPos1
           r = modulus(rFromStar)
           if ((r < grid%rStar1) .and. opaqueCore) then
              nTau = nTau + 1
              ksca(nTau) = 0.
              kabs(nTau) = 1.e20
              lambda(nTau) = lambda(nTau-1)
              hitCore = .true.
              endLoop = .true.
              posIndex(nTau,1:3) = posIndex(nTau-1,1:3)
              posOffset(nTau,1:3) = posOffset(nTau-1,1:3)
              rvecarray(ntau) = rVecArray(ntau-1)
              dlam(nTau) = dLam(nTau-1)
           endif

           rFromStar = rVec - grid%starPos2
           r = modulus(rFromStar)
           if ((r < grid%rStar2) .and. opaqueCore) then
              nTau = nTau + 1
              ksca(nTau) = 0.
              kabs(nTau) = 1.e20
              lambda(nTau) = lambda(nTau-1)
              hitCore = .true.
              endLoop = .true.
              posIndex(nTau,1:3) = posIndex(nTau-1,1:3)
              posOffset(nTau,1:3) = posOffset(nTau-1,1:3)

              dlam(nTau) = dLam(nTau-1)
           endif
        endif


        ! if the photon is outside the grid then the loop's finished

        if (grid%cartesian) then
           if (.not.( ((rVec%x >= grid%xAxis(1)) .and. (rVec%x <= grid%xAxis(grid%nx))) .and. &
                ((rVec%y >= grid%yAxis(1)) .and. (rVec%y <= grid%yAxis(grid%ny))) .and. &
                ((rVec%z >= grid%zAxis(1)) .and. (rVec%z <= grid%zAxis(grid%nz))) )) then
              endLoop = .true.
           endif



           if (grid%geometry == "ttauri") then
              if (hitDisk.and.(lambda(nTau)>distToDisk)) then
                 nTau = nTau + 1
                 ksca(nTau) = 1.e20
                 kabs(nTau) = 1.e20
                 lambda(nTau) = lambda(nTau-1)
                 endLoop = .true.
                 posIndex(nTau,1:3) = posIndex(nTau-1,1:3)
                 posOffset(nTau,1:3) = posOffset(nTau-1,1:3)
                 dlam(nTau) = dLam(nTau-1)
              endif
           endif

        else
           if ((r > grid%rAxis(grid%nr))) then
              endLoop = .true.
              if (nTau > 1) then
                 ksca(nTau) = ksca(nTau-1)
                 kabs(nTau) = kabs(nTau-1)
                 lambda(nTau) = lambda(nTau-1)
                 posIndex(nTau,1:3) = posIndex(nTau-1,1:3)
                 posOffset(nTau,1:3) = posOffset(nTau-1,1:3)
                 rvecarray(ntau) = rVecArray(ntau-1)
                 dlam(nTau) = dLam(nTau-1)
              endif
           endif
        endif


        ! if we are still in the grid

        if (.not.endLoop) then

           ! find the position in the grid

           call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)

           ! next optical depth point

           ! opacities are interpolated in three dimensions

           nTau = nTau + 1

           if (nTau > maxTau) then
              write(*,'(a)') "! maxTau exceeded"
              ntau = ntau - 1
              endloop = .true.
           endif

           if (interp) then
              kabs(nTau) = interpGridKappaAbs(grid,i1,i2,i3,iLambda,t1,t2,t3) 
              ksca(nTau) = interpGridKappaSca(grid,i1,i2,i3,iLambda,t1,t2,t3) 
           else
              kabs(nTau) = grid%kappaAbs(i1,i2,i3,iLambda)
              ksca(nTau) = grid%kappaSca(i1,i2,i3,iLambda)
           endif


           ! store the position indices and projected velocities etc

           posIndex(nTau,1) = i1
           posIndex(nTau,2) = i2
           posIndex(nTau,3) = i3
           posOffset(nTau,1) = t1
           posOffset(nTau,2) = t2
           posOffset(nTau,3) = t3

           rVecArray(nTau) = rVec

           lambda(nTau) = lambda(nTau-1) + dlambda
           dlam(nTau-1) = dLambda
           projVel(nTau) = interpGridVelocity(grid,i1,i2,i3,t1,t2,t3) .dot. uHat


           ! get the new length increment
           if (grid%cartesian) then
              if (i1 < grid%nx) then
                 dx = grid%xAxis(i1+1)-grid%xAxis(i1)
              else
                 dx = grid%xAxis(i1)-grid%xAxis(i1-1)
              endif
              if (i2 < grid%ny) then
                 dy = grid%yAxis(i2+1)-grid%yAxis(i2)
              else
                 dy = grid%yAxis(i2)-grid%yAxis(i2-1)
              endif
              if (i3 < grid%nz) then
                 dz = grid%zAxis(i3+1)-grid%zAxis(i3)
              else
                 dz = grid%zAxis(i3)-grid%zAxis(i3-1)
              endif

              fVec = VECTOR(dx, dy, dz)
              dlambda = abs(fvec%x*uHat%x) + abs(fVec%y*uHat%y) + abs(fVec%z*uHat%z)

!              dlambda = abs(fvec .dot. uHat)
!              dlambda = min(dx,dy,dz)

           else
	      r1 = grid%rAxis(i1)
	      if (i1 < grid%nr) then
		 r2 = grid%rAxis(i1+1)
	      else
		 r2 = grid%rAxis(i1-1)
	      endif
	      dr = abs (r2 - r1)
	      sinTheta = sqrt(1.-grid%muAxis(i2)**2)
	      if (i3 < grid%nphi) then
		 dx = r1 * sinTheta * cos(grid%phiAxis(i3)) - r2 * sinTheta * cos(grid%phiAxis(i3+1))
		 dy = r1 * sinTheta * sin(grid%phiAxis(i3)) - r2 * sinTheta * sin(grid%phiAxis(i3+1))
	      else
		 dx = r1 * sinTheta * cos(grid%phiAxis(i3)) - r2 * sinTheta * cos(grid%phiAxis(i3-1))
		 dy = r1 * sinTheta * sin(grid%phiAxis(i3)) - r2 * sinTheta * sin(grid%phiAxis(i3-1))
	      endif
	      
	      if (i2 < grid%nMu) then
		 dz = r1 * grid%muAxis(i2) - r2 * grid%muAxis(i2+1)
	      else
		 dz = r1 * grid%muAxis(i2) - r2 * grid%muAxis(i2-1)
	      endif
	      fVec = VECTOR(dx, dy, dz)
	      dlambda = min(modulus(fVec),dr)
!	      dlambda = abs(fvec%x*uHat%x) + abs(fVec%y*uHat%y) + abs(fVec%z*uHat%z)
           endif

           dlambda = dlambda / resFac


        endif

     enddo



     ! now compute the optical depths from the opacities

     do i = 2, nTau
        tauSca(i) = tauSca(i-1) + dlam(i-1) * (0.5*(ksca(i)+ksca(i-1)))
        tauAbs(i) = tauAbs(i-1) + dlam(i-1) * (0.5*(kabs(i)+kabs(i-1)))
        tauExt(i) = tauSca(i) + tauAbs(i)
     enddo



     iStart = 1

     if (grid%lineEmission) then
        if (nTau > 0) then

           ! line photons are treated as individuals

           if (.not.contPhoton) then

              ! loop through the projected velocities to find the resonance zones

              do i = 2, nTau-1

                 if ( ((projVel(i) < thisVel) .and. (thisVel <= projVel(i+1))) .or. &
                      ((projVel(i+1) < thisVel) .and. (thisVel <= projVel(i))) ) then

                    ! compute sobolev optical depth at this point and the next one

                    i1 = posIndex(i,1)
                    i2 = posIndex(i,2)
                    i3 = posIndex(i,3)

                    t1 = posOffset(i,1)
                    t2 = posOffset(i,2)
                    t3 = posOffset(i,3)

                    nu = cSpeed / (lambda0*angstromtocm)

                    if (.not.usePops) then
                       tauSob = interpGridChil(grid,i1,i2,i3,t1,t2,t3) / nu
                    else
                       chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
                       chil = 1.e10*chil* (grid%n(i1,i2,i3,mLevel)-( ( gM / gN) * grid%n(i1,i2,i3,nLevel)))
                       tauSob = chil / nu
                    endif

                    tauSob1 = tauSob / directionalDeriv(grid,rVecArray(i),i1,i2,i3,uHat)



                    i1 = posIndex(i+1,1)
                    i2 = posIndex(i+1,2)
                    i3 = posIndex(i+1,3)

                    t1 = posOffset(i+1,1)
                    t2 = posOffset(i+1,2)
                    t3 = posOffset(i+1,3)

                    nu = cSpeed / (lambda0*angstromtocm)


                    if (.not.usePops) then
                       tauSob = interpGridChil(grid,i1,i2,i3,t1,t2,t3) / nu
                    else
                       chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
                       chil = 1.e10*chil* (grid%n(i1,i2,i3,mLevel)-( ( gM / gN) * grid%n(i1,i2,i3,nLevel)))
                       tauSob = chil / nu
                    endif

                    tauSob2 = tauSob / directionalDeriv(grid,rVecArray(i+1),i1,i2,i3,uHat)

                    ! interpolate in velocity space to get sobolev depth

                    x = 1.
                    if ((projVel(i+1)-projVel(i)) /= 0.) then
                       x = (thisVel-projVel(i))/(projVel(i+1)-projVel(i))
                    endif
                    tauSob = tauSob1  + x * (tauSob2 - tauSob1)

                    ! somethings gone wrong - usually numerical problems 

                    if (tauSob < 0.) then
                       write(*,*) "tau sob",tausob
                       write(*,*) "this Vel",thisVel*cSpeed/1.e5
                       write(*,*) "modulus vvec",modulus(vVec)*cSpeed/1.e5
                       write(*,*) "wavelength",wavelength,lambda0
                       write(*,*) "proj",projVel(i)*cSpeed/1.e5,projVel(i+1)*cSpeed/1.e5
                       write(*,*) "m,n",mLevel, nLevel
                       write(*,*) "g",gM,gN
                       write(*,*) "f",fStrength
                       write(*,*) "popm",grid%n(i1,i2,i3,mLevel)
                       write(*,*) "popn",grid%n(i1,i2,i3,nLevel)
                    endif

                    ! add depths

                    tauExt(i:nTau) = tauExt(i:nTau) + tauSob
                    tauAbs(i:nTau) = tauAbs(i:nTau) + tauSob


                 endif

                 if (grid%resonanceLine) then


                    if ( ((projVel(i) < thisVel2) .and. (thisVel2 <= projVel(i+1))) .or. &
                         ((projVel(i+1) < thisVel2) .and. (thisVel2 <= projVel(i))) ) then
                       
                       ! compute sobolev optical depth at this point and the next one
                       
                       i1 = posIndex(i,1)
                       i2 = posIndex(i,2)
                       i3 = posIndex(i,3)
                       
                       t1 = posOffset(i,1)
                       t2 = posOffset(i,2)
                       t3 = posOffset(i,3)

                       nu = cSpeed / (lambda2*angstromtocm)

                       if (.not.usePops) then
                          tauSob = interpGridChil(grid,i1,i2,i3,t1,t2,t3) / nu
                       else
                          chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
                          chil = 1.e10*chil* (grid%n(i1,i2,i3,mLevel)-( ( gM / gN) * grid%n(i1,i2,i3,nLevel)))
                          tauSob = chil / nu
                       endif
                       
                       tauSob1 = tauSob / directionalDeriv(grid,rVecArray(i),i1,i2,i3,uHat)
                       

                       
                       i1 = posIndex(i+1,1)
                       i2 = posIndex(i+1,2)
                       i3 = posIndex(i+1,3)
                       
                       t1 = posOffset(i+1,1)
                       t2 = posOffset(i+1,2)
                       t3 = posOffset(i+1,3)

                       nu = cSpeed / (lambda2*angstromtocm)


                       if (.not.usePops) then
                          tauSob = interpGridChil(grid,i1,i2,i3,t1,t2,t3) / nu
                       else
                          chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
                          chil = 1.e10*chil* (grid%n(i1,i2,i3,mLevel)-( ( gM / gN) * grid%n(i1,i2,i3,nLevel)))
                          tauSob = chil / nu
                       endif
                       
                       tauSob2 = tauSob / directionalDeriv(grid,rVecArray(i+1),i1,i2,i3,uHat)
                       
                       ! interpolate in velocity space to get sobolev depth
                       
                       x = 1.
                       if ((projVel(i+1)-projVel(i)) /= 0.) then
                          x = (thisVel2-projVel(i))/(projVel(i+1)-projVel(i))
                       endif
                       tauSob = tauSob1  + x * (tauSob2 - tauSob1)

                       ! somethings gone wrong - usually numerical problems 

                       if (tauSob < 0.) then
                          write(*,*) "tau sob",tausob
                          write(*,*) "this Vel",thisVel2*cSpeed/1.e5
                          write(*,*) "modulus vvec",modulus(vVec)*cSpeed/1.e5
                          write(*,*) "wavelength",wavelength,lambda2
                          write(*,*) "proj",projVel(i)*cSpeed/1.e5,projVel(i+1)*cSpeed/1.e5
                          write(*,*) "m,n",mLevel, nLevel
                          write(*,*) "g",gM,gN
                          write(*,*) "f",fStrength
                          write(*,*) "popm",grid%n(i1,i2,i3,mLevel)
                          write(*,*) "popn",grid%n(i1,i2,i3,nLevel)
                       endif

                       ! add depths
                       
                       tauExt(i:nTau) = tauExt(i:nTau) + tauSob
                       tauAbs(i:nTau) = tauAbs(i:nTau) + tauSob


                    endif
                 endif

              enddo


           else

              ! for continuum photons we treat an array of optical depths corresponding
              ! to different velocities/frequencies

              ! initialize array

              tauCont(1:nTau,1:nLambda) = 0.


              ! we can ignore this if the line is very thin

              if (.not.thinLine) then

                 ! loop over all wavelength bins


                 do j = 1, nLambda

                    ! compute projected velocity of this bin

                    thisVel = real(j-1)/real(nLambda-1)
                    thisVel = lamStart + thisVel*(lamEnd-lamStart)
                    thisVel = (thisVel-lambda0)/lambda0
                    thisVel = -thisVel  + (uHat .dot. vVec) - (wavelength - lambda0)/lambda0

                    ! now find resonances

                    do i = 1, nTau - 1

                       if ( ((projVel(i) < thisVel) .and. (thisVel <= projVel(i+1))) .or. &
                            ((projVel(i+1) < thisVel) .and. (thisVel <= projVel(i))) ) then

                          ! this bit is the same as for the line photon - linear interpolation
                          ! in velocity space to get sobolev optical depth

                          i1 = posIndex(i,1)
                          i2 = posIndex(i,2)
                          i3 = posIndex(i,3)

                          t1 = posOffset(i,1)
                          t2 = posOffset(i,2)
                          t3 = posOffset(i,3)

                          nu = cSpeed / (lambda0*angstromtocm)
                          if (.not. usePops) then
                             tauSob = interpGridChil(grid,i1,i2,i3,t1,t2,t3) / nu
                          else
                             chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
                             chil = 1.e10*chil* (grid%n(i1,i2,i3,mLevel)-( ( gM / gN) * grid%n(i1,i2,i3,nLevel)))
                             tauSob = chil / nu
                          endif

                          tauSob1 = tauSob / directionalDeriv(grid,rVecArray(i),i1,i2,i3,uHat)

                          i1 = posIndex(i+1,1)
                          i2 = posIndex(i+1,2)
                          i3 = posIndex(i+1,3)

                          t1 = posOffset(i+1,1)
                          t2 = posOffset(i+1,2)
                          t3 = posOffset(i+1,3)

                          nu = cSpeed / (lambda0*angstromtocm)
                          if (.not.usePops) then
                             tauSob = interpGridChil(grid,i1,i2,i3,t1,t2,t3) / nu
                          else
                             chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
                             chil = 1.e10*chil* (grid%n(i1,i2,i3,mLevel)-( ( gM / gN) * grid%n(i1,i2,i3,nLevel)))
                             tauSob = chil / nu
                          endif

                          tauSob2 = tauSob / directionalDeriv(grid,rVecArray(i+1),i1,i2,i3,uHat)

                          x = 1.
                          if ((projVel(i+1)-projVel(i)) /= 0.) then
                             x = (thisVel-projVel(i))/(projVel(i+1)-projVel(i))
                          endif
                          tauSob = tauSob1  + x * (tauSob2 - tauSob1)


                          ! might have gone wrong due to rounding errors

                          if (tauSob < 0.) then
                             write(*,*) "tau sob",tausob,tausob1,tausob2
                             write(*,*) "this Vel",thisVel*cSpeed/1.e5
                             write(*,*) "modulus vvec",modulus(vVec)*cSpeed/1.e5
                             write(*,*) "wavelength",wavelength,lambda0
                             write(*,*) "proj",projVel(i)*cSpeed/1.e5
                             write(*,*) modulus(rVecArray(i+1)-grid%starpos1)/grid%rStar1, &
modulus(rVecArray(i)-grid%starpos1)/grid%rStar1
                             write(*,*) modulus(rVecArray(i+1)-grid%starpos2)/grid%rStar2, &
modulus(rVecArray(i)-grid%starpos2)/grid%rStar2
                             write(*,*) "inuse",grid%inuse(i1,i2,i3)
                             write(*,*) posindex(i,1:3),posoffset(i,1:3)
                             write(*,*) posindex(i+1,1:3),posoffset(i+1,1:3)
                          endif

                          ! add optical depths


                          tauCont(i:nTau,j) = tauCont(i:nTau,j) + tauSob

                       endif ! resonance zone

                    enddo ! over optical depth array

                 enddo   ! over wavelength array

              endif  ! thin line?

           endif ! continuum photon?

        endif

     endif ! stateq grid?


  endif  ! cartesian grid?




666 continue

end subroutine integratePath

   
    
   
end module path_integral
