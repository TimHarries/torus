module path_integral

  use vector_mod             ! vector maths 
  use gridtype_mod           ! opacity grid
  use grid_mod               ! opacity grid routines
  use math_mod               ! misc maths routines
  use photon_mod             ! the photon module
  use constants_mod          ! physical constants
  use amr_mod
  
  implicit none
!
! This subroutine performs an integration through the grid in order
! to compute the scattering and absorption optical depths.
!
  contains

subroutine integratePath(wavelength,  lambda0, vVec, aVec, uHat, Grid,  lambda, &
     tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton, &
     lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, redRegion, rStar, &
     coolStarPosition, usePops, mLevel, nLevel, fStrength, gM, gN, localTau, interp)


  type(GRIDTYPE) :: grid                            ! the opacity grid

  real :: localTau
  logical :: usePops
  integer :: mLevel, nLevel
  real :: fStrength, gM, GN
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

  real :: sinTheta, dr

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
  real :: lambda(:),dlambda                         ! distance arrays
  real :: tauExt(:), tauAbs(:), tauSca(:)           ! optical depths
  real, dimension(:,:) :: tauCont


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

  logical :: firstTime = .true.

  if (firstTime) then
     print *, 'integratePath may need some bug fixes for out-by-one errors'
     print *, '  changes were made by nhs for the ttauri geometry, but these may'
     print *, '  not be correct. Other geometries were not checked. The changes were '
     print *, '  made around lines 553,554,492,906.'
     firstTime = .false.
  end if

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
              if (.not.grid%oneKappa) then
                 kabsInterp = interpGridKappaAbs(grid,i1,i2,i3,iLambda,t1,t2,t3) 
                 kscaInterp = interpGridKappaSca(grid,i1,i2,i3,iLambda,t1,t2,t3) 
              else
                 r = interpGridScalar2(grid%rho,grid%na1,grid%na2,grid%na3,i1,i2,i3,t1,t2,t3)
                 kabsInterp = grid%oneKappaAbs(1,iLambda) * r
                 kscaInterp = grid%oneKappaSca(1,iLambda) * r
              endif
           else
              if (.not.grid%oneKappa) then
                 kabsInterp = grid%kappaAbs(i1,i2,i3,iLambda)
                 kscaInterp = grid%kappaSca(i1,i2,i3,iLambda)
              else
                 kabsInterp = grid%oneKappaAbs(1,iLambda)*grid%rho(i1,i2,i3)
                 kscaInterp = grid%oneKappaSca(1,iLambda)*grid%rho(i1,i2,i3)
              endif
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
        if (.not.grid%oneKappa) then
           kabs(nTau) = interpGridKappaAbs(grid,i1,i2,i3,iLambda,t1,t2,t3) 
           ksca(nTau) = interpGridKappaSca(grid,i1,i2,i3,iLambda,t1,t2,t3) 
        else
           r = interpGridScalar2(grid%rho,grid%na1,grid%na2,grid%na3,i1,i2,i3,t1,t2,t3)
           kabs(nTau) = grid%oneKappaAbs(1,iLambda) * r
           ksca(nTau) = grid%oneKappaSca(1,iLambda) * r
        endif
     else
        if (.not.grid%oneKappa) then
           kabs(nTau) = grid%kappaAbs(i1,i2,i3,iLambda)
           ksca(nTau) = grid%kappaSca(i1,i2,i3,iLambda)
        else
           kabs(nTau) = grid%oneKappaAbs(1,iLambda) * grid%rho(i1,i2,i3)
           ksca(nTau) = grid%oneKappaSca(1,iLambda) * grid%rho(i1,i2,i3)
        endif
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
        tauSob = tauSob / directionalDeriv(grid,s2o(rVec),i1,i2,i3,uHat)

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
        dlambda  = abs(dr)!!!!!!!!!!!!!!!1min(modulus(fVec),dr)
!       dlambda = abs(fvec%x*uHat%x) + abs(fVec%y*uHat%y) + abs(fVec%z*uHat%z)


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
              if (ntau > 2) then
                 dlam(nTau-1) = dLam(nTau-2)
              end if
              projVel(nTau) = projVel(nTau-1)
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
                 if (nTau > 2) then 
                    dlam(nTau-1) = dLam(nTau-2)
                 end if
                 projVel(nTau) = projVel(nTau-1)
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
              if (.not.grid%oneKappa) then
                 kabs(nTau) = interpGridKappaAbs(grid,i1,i2,i3,iLambda,t1,t2,t3) 
                 ksca(nTau) = interpGridKappaSca(grid,i1,i2,i3,iLambda,t1,t2,t3) 
              else
                 r = interpGridScalar2(grid%rho,grid%na1,grid%na2,grid%na3,i1,i2,i3,t1,t2,t3)
                 kabs(nTau) = grid%oneKappaAbs(1,iLambda) * r
                 ksca(nTau) = grid%oneKappaSca(1,iLambda) * r
              endif
           else
              if (.not.grid%oneKappa) then
                 kabs(nTau) = grid%kappaAbs(i1,i2,i3,iLambda)
                 ksca(nTau) = grid%kappaSca(i1,i2,i3,iLambda)
              else
                 kabs(nTau) = grid%oneKappaAbs(1,iLambda) * grid%rho(i1,i2,i3)
                 ksca(nTau) = grid%oneKappaSca(1,iLambda) * grid%rho(i1,i2,i3)
              endif
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
              dlambda = abs(dr)!!!!!!!!!!!!!!!!!min(modulus(fVec),dr)
!             dlambda = abs(fvec%x*uHat%x) + abs(fVec%y*uHat%y) + abs(fVec%z*uHat%z)
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

                    tauSob1 = tauSob / directionalDeriv(grid,s2o(rVecArray(i)),i1,i2,i3,uHat)



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

                    tauSob2 = tauSob / directionalDeriv(grid,s2o(rVecArray(i+1)),i1,i2,i3,uHat)

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
                       
                       tauSob1 = tauSob / directionalDeriv(grid,s2o(rVecArray(i)),i1,i2,i3,uHat)
                       

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
                       
                       tauSob2 = tauSob / directionalDeriv(grid,s2o(rVecArray(i+1)),i1,i2,i3,uHat)
                       
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


           else ! contPhoton

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

                    do i = 2, nTau - 1

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

                          tauSob1 = tauSob / directionalDeriv(grid,s2o(rVecArray(i)),i1,i2,i3,uHat)

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

                          tauSob2 = tauSob / directionalDeriv(grid,s2o(rVecArray(i+1)),i1,i2,i3,uHat)

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


subroutine integratePathAMR(wavelength,  lambda0, vVec, aVec, uHat, Grid, &
     lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb,&
     contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, &
     redRegion, usePops, mLevel, nLevel, &
     fStrength, gM, gN, localTau,sampleFreq,error)

  ! should we add the 'interp' argument and implement it?

  ! error codes are:  -10: too few samples made (the code should be fixed to 
  !                          stop this being a problem)
  !                   -20: the photon passed too close to a cell boundary and 
  !                          returnSamples decided to abort it 
 

  implicit none

  real, intent(in)          :: wavelength             ! the wavelength 
  real, intent(in)          :: lambda0                ! rest wavelength of line
  type(OCTALVECTOR), intent(in)  :: vVec                   ! velocity vector
  type(OCTALVECTOR), intent(in)  :: aVec                   ! starting position vector
  type(OCTALVECTOR), intent(in)  :: uHat                   ! direction
  type(GRIDTYPE), intent(in):: grid                   ! the opacity grid
  integer, intent(in)       :: maxTau 
  real, dimension(:), intent(out) :: lambda           ! path distance array
  real, dimension(:), intent(out) :: tauExt           ! optical depth
  real, dimension(:), intent(out) :: tauAbs           ! optical depth
  real, dimension(:), intent(out) :: tauSca           ! optical depth
  integer, intent(out)      :: nTau                   ! size of optical depth arrays
  logical, intent(in)       :: opaqueCore             ! is the core opaque
  real, intent(out)         :: escProb                ! the escape probability
  real, intent(in)          :: lamStart, lamEnd
  integer, intent(in)       :: nLambda
  real,intent(out),dimension(:,:) :: tauCont
  logical, intent(out)      :: hitcore                ! has the photon hit the core
  logical, intent(in)       :: thinLine               ! ignore line absorption of continuum
  logical, intent(in)       :: redRegion              ! use raman scattering red region opacities
  logical, intent(in)       :: usePops
  integer, intent(in)       :: mLevel, nLevel
  real, intent(in)          :: fStrength, gM, gN
  real, intent(out)         :: localTau
  real, intent(in)          :: sampleFreq             ! max. samples per grid cell
  integer, intent(out)      :: error                  ! error code returned


  real                      :: rho(1:maxTau)          ! density
  real                      :: temperature(1:maxTau)  ! temperature
  type(octalVector)         :: aVecOctal              ! octalVector version of 'aVec'
  type(octalVector)         :: uHatOctal              ! octalVector version of 'uHat'
  type(OCTALVECTOR)         :: rVec                   ! position vector
  logical                   :: contPhoton             ! is this a continuum photon?
  real, dimension(1:maxTau) :: chiLine                ! line opacities
  real, dimension(1:maxTau,1:grid%maxLevels) :: levelPop  ! stateq level pops
  real, dimension(1:maxTau) :: tauSob, tauSob2        ! Sobolev optical depths
  real, dimension(1:maxTau) :: escProbArray           ! escape probabilities
  type(vector), dimension(1:maxTau) :: velocity       ! 
  real, dimension(1:maxTau) :: velocityDeriv          ! directional derivative
  
  real :: nu, nu2                                        ! frequencies

  integer :: iStart                                 ! starting index

  real :: tauSobScalar,tauSobScalar1,tauSobScalar2  ! Sobolev optical depths

  real :: thisVel
  integer :: ilambda                                ! wavelength index


  integer, parameter :: maxLambda = 2000             ! max size of tau arrays
  real :: dlambda(maxTau)                           ! distance increment array

  real :: projVel(maxTau)
  real :: kabs(maxtau), ksca(maxtau)                ! opacities
  real :: x                                         ! multipliers

  real :: chil                                      ! line opacity
  
  real :: lambda2, thisVel2

  integer :: i, j                   ! counters

  ! initialize variables

  hitcore = .false.
  rVec = aVec
  escProb = 1.
  error = 0

  ! locate this wavelength in the grid of wavelengths

  if (grid%flatspec.or.(grid%doRaman)) then
     iLambda = 1
  else
     call locate(grid%lamArray, grid%nLambda, wavelength, iLambda)
  endif

  if (.not. grid%adaptive) then
     print *, 'integratePathAMR called on a non-adaptive grid!'
     stop
  end if

  nTau = 0
  
  ! convert some variables to the appropriate kind for startReturnSamples
  aVecOctal = aVec
  uHatOctal = uHat
     
  if (.not. grid%lineEmission) then

     ! sample the grid along the path of the ray
     CALL startReturnSamples (aVecOctal,uHatOctal,grid,sampleFreq,nTau,       &
                              maxTau,opaqueCore,hitcore,usePops,iLambda,error,&
                              lambda,kappaAbs=kAbs,kappaSca=kSca,velocity=velocity,&
                              velocityDeriv=velocityDeriv,chiLine=chiLine,    &
                              levelPop=levelPop,rho=rho, temperature=temperature)

     if (.not.redRegion) then
        ! first optical depths are all zero of course

        tauExt(1) = 0.
        tauAbs(1) = 0.
        tauSca(1) = 0.
                     
        
        do i = 2, nTau, 1 
           dlambda(i) = lambda(i)-lambda(i-1)
           tauSca(i) = tauSca(i-1) + dlambda(i)*0.5*(ksca(i-1)+ksca(i))
           tauAbs(i) = tauAbs(i-1) + dlambda(i)*0.5*(kabs(i-1)+kabs(i))
           if ((ksca(i) < 0.).or.(kabs(i)<0.)) then
              write(*,*) "negative opacity"
              do;enddo
           endif
        enddo

        tauExt(1:nTau) = tauSca(1:nTau)+tauAbs(1:nTau)

        else
           print *, "No code for handling redRegion"
           stop
 
        endif

  else ! grid%lineEmission

     ! some geometries produce photons at the origin

     if (modulus(aVecOctal) < 0.01 * grid%halfSmallestSubcell) then
        if (grid%rCore < 1.e-10) then
           aVecOctal = (1.0001_oc * grid%halfSmallestSubcell) * uHatOctal
        else 
           aVecOctal = (1.0001_oc * grid%rCore) * uHatOctal
        end if
     end if

     ! find the project velocity

     thisVel = (lambda0 - wavelength)/lambda0
     thisVel = thisVel  + (uHat .dot. vVec)  

     if (grid%resonanceLine) then
        lambda2 = grid%lambda2
        thisVel2 = (lambda2 - wavelength)/lambda2
        thisVel2 = thisVel2  + (uHat .dot. vVec) 
     endif


     escProb = 1.

     CALL startReturnSamples (aVecOctal,uHatOctal,grid,sampleFreq,nTau,       &
                              maxTau,opaqueCore,hitcore,usePops,iLambda,error,&
                              lambda,kappaAbs=kAbs,kappaSca=kSca,velocity=velocity,&
                              velocityDeriv=velocityDeriv,chiLine=chiLine,    &
                              levelPop=levelPop,rho=rho,temperature=temperature)

     if (nTau <= 2) then
        !print *, 'The code does not yet include routines for handling ',&
        !         'the case where nTau <= 2 (in integratePathAMR)'
        error = -10
        return
     end if


     ! compute the escape probability if this is a line photon

     if (.not.contPhoton) then

        if (usePops) then
           chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
           chiLine(1:nTau) = 1.e10 * chil * (levelPop(1:nTau,mLevel)&
                          -( ( gM / gN) * (levelPop(1:nTau,nLevel))))
        endif

        nu  = cSpeed / (lambda0*angstromtocm)
        tauSob(1:nTau) = chiLine(1:nTau)  / (nu  * velocityDeriv(1:nTau))
     
        if (grid%resonanceLine) then
           nu2 = cSpeed / (lambda2*angstromtocm)
           tauSob2(1:nTau) = chiLine(1:nTau)  / (nu2  * velocityDeriv(1:nTau))
        end if


        where (tauSob(1:nTau) < 0.1)
           escProbArray(1:nTau) = 1.0d0-tauSob(1:nTau)*0.5d0*(1.0d0 - tauSob(1:nTau)/3.0d0* &
                       (1. - tauSob(1:nTau)*0.25d0*(1.0d0 - 0.20d0*tauSob(1:nTau))))
        elsewhere (tauSob(1:nTau) < 15.) 
           escProbArray(1:nTau)  = (1.0d0-exp(-tauSob(1:nTau)))/tauSob(1:nTau)
        elsewhere  
           escProbArray(1:nTau) = 1.d0/tauSob(1:nTau)
        end where
        localTau = tauSob(1)

     endif


     ! projected velocities 

     if (contPhoton) then
        forall (i = 1:nTau)
           projVel(i) = velocity(i) .dot. uHat
        end forall
     else
        projVel(1:nTau) = thisVel
     endif


     dlambda(1:nTau-1) = lambda(2:nTau) - lambda(1:nTau-1)

     ! now compute the optical depths from the opacities

     tauExt(1:2) = 0.
     tauAbs(1:2) = 0.
     tauSca(1:2) = 0.

     do i = 2, nTau, 1
        tauSca(i) = tauSca(i-1) + dlambda(i-1)*0.5*(ksca(i-1)+ksca(i))
        tauAbs(i) = tauAbs(i-1) + dlambda(i-1)*0.5*(kabs(i-1)+kabs(i))
           if ((ksca(i) < 0.).or.(kabs(i)<0.)) then
              write(*,*) "negative opacity"
              do;enddo
           endif

     enddo

     tauExt(1:nTau) = tauSca(1:nTau) + tauAbs(1:nTau)

     iStart = 1

     if (grid%lineEmission) then
        if (nTau > 0) then

           ! line photons are treated as individuals

           if (.not.contPhoton) then

              ! loop through the projected velocities to find the resonance zones

              do i = 2, nTau-1, 1

                 if ( ((projVel(i) < thisVel) .and. (thisVel <= projVel(i+1))) .or. &
                      ((projVel(i+1) < thisVel) .and. (thisVel <= projVel(i))) ) then

                    ! compute sobolev optical depth at this point and the next one

                    tauSobScalar1 = tauSob(i)

                    tauSobScalar2 = tauSob(i+1)

                    ! interpolate in velocity space to get sobolev depth

                    x = 1.
                    if ((projVel(i+1)-projVel(i)) /= 0.) then
                       x = (thisVel-projVel(i))/(projVel(i+1)-projVel(i))
                    endif
                    tauSobScalar = tauSobScalar1  + x * (tauSobScalar2 - tauSobScalar1)




                    if (tauSobScalar < 0. .and.  abs(tauSobScalar) < 1.0e-3) then
                       ! somethings gone wrong - usually numerical problems 
                       ! Correction forced (R.K. Oct-15-2002)
                       tauSobScalar = 1.0e-10
                    elseif (tauSobScalar < 0.) then
                       write(*,*) "tau sob",tausobScalar
                       write(*,*) "this Vel",thisVel*cSpeed/1.e5
                       write(*,*) "modulus vvec",modulus(vVec)*cSpeed/1.e5
                       write(*,*) "wavelength",wavelength,lambda0
                       write(*,*) "proj",projVel(i)*cSpeed/1.e5,projVel(i+1)*cSpeed/1.e5
                       write(*,*) "m,n",mLevel, nLevel
                       write(*,*) "g",gM,gN
                       write(*,*) "f",fStrength
                       write(*,*) "popm",levelPop(i,mLevel)
                       write(*,*) "popn",levelPop(i,nLevel)
                    endif

                    ! add depths

                    tauExt(i:nTau) = tauExt(i:nTau) + tauSobScalar
                    tauAbs(i:nTau) = tauAbs(i:nTau) + tauSobScalar


                 endif

                 if (grid%resonanceLine) then


                    if ( ((projVel(i) < thisVel2) .and. (thisVel2 <= projVel(i+1))) .or. &
                         ((projVel(i+1) < thisVel2) .and. (thisVel2 <= projVel(i))) ) then
                       
                       ! compute sobolev optical depth at this point and the next one
                       
                       tauSobScalar1 = tauSob2(i)

                       tauSobScalar2 = tauSob2(i+1)
                                             ! interpolate in velocity space to get sobolev depth
                       
                       x = 1.
                       if ((projVel(i+1)-projVel(i)) /= 0.) then
                          x = (thisVel2-projVel(i))/(projVel(i+1)-projVel(i))
                       endif
                       tauSobScalar = tauSobScalar1  + x * (tauSobScalar2 - tauSobScalar1)

                       ! somethings gone wrong - usually numerical problems 
                       if (tauSobScalar < 0. .and.  abs(tauSobScalar) < 1.0e-3) then
                          ! somethings gone wrong - usually numerical problems 
                          ! Correction forced (R.K. Oct-15-2002)
                          tauSobScalar = 1.0e-10
                       elseif (tauSobScalar < 0.) then
                          write(*,*) "tau sob",tausobScalar
                          write(*,*) "this Vel",thisVel*cSpeed/1.e5
                          write(*,*) "modulus vvec",modulus(vVec)*cSpeed/1.e5
                          write(*,*) "wavelength",wavelength,lambda0
                          write(*,*) "proj",projVel(i)*cSpeed/1.e5,projVel(i+1)*cSpeed/1.e5
                          write(*,*) "m,n",mLevel, nLevel
                          write(*,*) "g",gM,gN
                          write(*,*) "f",fStrength
                          write(*,*) "popm",levelPop(i,mLevel)
                          write(*,*) "popn",levelPop(i,nLevel)
                       endif

                       ! add depths
                       
                       tauExt(i:nTau) = tauExt(i:nTau) + tauSobScalar
                       tauAbs(i:nTau) = tauAbs(i:nTau) + tauSobScalar


                    endif
                 endif ! (grid%resonanceLine)

              enddo


           else ! contPhoton

              ! for continuum photons we treat an array of optical depths corresponding
              ! to different velocities/frequencies

              ! initialize array

              tauCont(1:nTau,:) = 0.


              ! we can ignore this if the line is very thin

              if (.not.thinLine) then

                 if (usePops) then
                    chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
                    chiLine(1:nTau) = 1.e10 * chil * (levelPop(1:nTau,mLevel)&
                                   -( ( gM / gN) * (levelPop(1:nTau,nLevel))))
                 endif
                    
                 nu  = cSpeed / (lambda0*angstromtocm)
                 tauSob(1:nTau) = chiLine(1:nTau)  / (nu  * velocityDeriv(1:nTau))
                 
                 ! loop over all wavelength bins

                 do j = 1, nLambda, 1

                    ! compute projected velocity of this bin

                    thisVel = real(j-1)/real(nLambda-1)
                    thisVel = lamStart + thisVel*(lamEnd-lamStart)
                    thisVel = (thisVel-lambda0)/lambda0
                    thisVel = -thisVel  + (uHat .dot. vVec) - (wavelength - lambda0)/lambda0

                    
                    ! loop through the projected velocities to find the resonance zones
 
                    do i = 2, nTau-1, 1

                       if ( ((projVel(i) < thisVel) .and. (thisVel <= projVel(i+1))) .or. &
                            ((projVel(i+1) < thisVel) .and. (thisVel <= projVel(i))) ) then

                          ! this bit is the same as for the line photon - linear interpolation
                          !   in velocity space to get sobolev optical depth
                          
                          ! compute sobolev optical depth at this point and the next one
   
                          tauSobScalar1 = tauSob(i)
   
                          tauSobScalar2 = tauSob(i+1)
   
   
                          ! interpolate in velocity space to get sobolev depth
   
                          x = 1.
                          if ((projVel(i+1)-projVel(i)) /= 0.) then
                             x = (thisVel-projVel(i))/(projVel(i+1)-projVel(i))
                          endif
                          tauSobScalar = tauSobScalar1  + x * (tauSobScalar2 - tauSobScalar1)

                          if (tauSobScalar < 0. .and.  abs(tauSobScalar) < 1.0e-3) then
                             ! somethings gone wrong - usually numerical problems 
                             ! Correction forced (R.K. Oct-15-2002)
                             tauSobScalar = 1.0e-10
                          elseif (tauSobScalar < 0.) then
                             write(*,*) "tau sob",tausobScalar
                             write(*,*) "this Vel",thisVel*cSpeed/1.e5
                             write(*,*) "modulus vvec",modulus(vVec)*cSpeed/1.e5
                             write(*,*) "wavelength",wavelength,lambda0
                             write(*,*) "proj",projVel(i)*cSpeed/1.e5,projVel(i+1)*cSpeed/1.e5
                             write(*,*) "m,n",mLevel, nLevel
                             write(*,*) "g",gM,gN
                             write(*,*) "f",fStrength
                             write(*,*) "popm",levelPop(i,mLevel)
                             write(*,*) "popn",levelPop(i,nLevel)
                          endif

                          ! add optical depths

                          tauCont(i:nTau,j) = tauCont(i:nTau,j) + tauSobScalar

                       endif ! resonance zone

                    enddo ! over optical depth array

                 enddo   ! over wavelength array

              endif  ! thin line?

           endif ! continuum photon?

        endif ! nTau > 0

     endif ! grid%lineEmission


  endif ! .not. grid%lineEmission


end subroutine integratePathAMR

subroutine integratePathAMR2D(wavelength,  lambda0, vVec, aVec, uHat, Grid, &
     lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb,&
     contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, &
     redRegion, usePops, mLevel, nLevel, &
     fStrength, gM, gN, localTau,sampleFreq,error)

  ! should we add the 'interp' argument and implement it?

  ! error codes are:  -10: too few samples made (the code should be fixed to 
  !                          stop this being a problem)
  !                   -20: the photon passed too close to a cell boundary and 
  !                          returnSamples decided to abort it 
 

  implicit none

  logical :: escaped
  type(OCTALVECTOR) :: octVec, octvec2d
  integer :: subcell
  real(kind=octalKind) :: tval

  real, intent(in)          :: wavelength             ! the wavelength 
  real, intent(in)          :: lambda0                ! rest wavelength of line
  type(OCTALVECTOR), intent(in)  :: vVec                   ! velocity vector
  type(OCTALVECTOR), intent(in)  :: aVec                   ! starting position vector
  type(OCTALVECTOR), intent(in)  :: uHat                   ! direction
  type(GRIDTYPE), intent(in):: grid                   ! the opacity grid
  integer, intent(in)       :: maxTau
  real,intent(out),dimension(:) :: lambda             ! path distance array
  real,intent(out),dimension(:) :: tauExt             ! optical depth
  real,intent(out),dimension(:) :: tauAbs             ! optical depth
  real,intent(out),dimension(:) :: tauSca             ! optical depth
  integer, intent(out)      :: nTau                   ! size of optical depth arrays
  real, intent(out)         :: escProb                ! the escape probability
  real, intent(in)          :: lamStart, lamEnd
  integer, intent(in)       :: nLambda
  real,intent(out),dimension(:,:) :: tauCont
  logical, intent(out)      :: hitcore                ! has the photon hit the core
  logical, intent(in)       :: thinLine               ! ignore line absorption of continuum
  logical, intent(in)       :: redRegion              ! use raman scattering red region opacities
  logical, intent(in)       :: usePops
  integer, intent(in)       :: mLevel, nLevel
  real, intent(in)          :: fStrength, gM, gN
  real, intent(out)         :: localTau
  real, intent(in)          :: sampleFreq             ! max. samples per grid cell
  integer, intent(out)      :: error                  ! error code returned


  real                      :: rho(1:maxTau)          ! density
!  type(octalVector)         :: aVecOctal              ! octalVector version of 'aVec'
!  type(octalVector)         :: uHatOctal              ! octalVector version of 'uHat'
  type(octalvector)         :: rVec                   ! position vector
  logical                   :: contPhoton             ! is this a continuum photon?
  real, dimension(1:maxTau,1:grid%maxLevels) :: levelPop  ! stateq level pops
  type(vector), dimension(1:maxTau) :: velocity       ! 
  real, dimension(1:maxTau) :: velocityDeriv          ! directional derivative
  
  real :: nu, nu2                                        ! frequencies


  integer :: ilambda                                ! wavelength index


  integer, parameter :: maxLambda = 2000             ! max size of tau arrays
  real :: dlambda(maxTau)                           ! distance increment array

  real :: projVel(maxTau)
  real :: kabs(maxtau), ksca(maxtau)                ! opacities
  real :: x                                         ! multipliers

  real :: chil                                      ! line opacity

   type(OCTAL), pointer :: thisOctal
   type(OCTAL),pointer :: oldOctal

  real :: kappaScaReal, kappaAbsReal
  
  real :: lambda2, thisVel2

  integer :: i                   ! counters

  logical :: opaqueCore

  ! initialize variables

  hitcore = .false.
  rVec = aVec
  escProb = 1.
  error = 0

  ! locate this wavelength in the grid of wavelengths

  if (grid%flatspec.or.(grid%doRaman)) then
     iLambda = 1
  else
     call locate(grid%lamArray, grid%nLambda, wavelength, iLambda)
  endif

  if (.not. grid%adaptive) then
     print *, 'integratePathAMR called on a non-adaptive grid!'
     stop
  end if

  nTau = 1
  tauSca(nTau) = 0.
  tauAbs(nTau) = 0.
  lambda(nTau) = 0.
  escaped = .false.
  !call amrGridValues(grid%octreeRoot, octVec, iLambda=iLambda, &
  !     foundOctal=thisOctal, foundSubcell=subcell, &
  !     kappaAbs=kappaAbsReal,kappaSca=kappaScaReal, grid=grid)
  oldOctal => grid%octreeRoot
    do while (.not.escaped) 
       octVec = rVec

       if (.not.inOctal(grid%octreeRoot, octVec)) escaped = .true.

       if (.not.escaped) then
          call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal,iLambda=iLambda, &
               foundOctal=thisOctal, foundSubcell=subcell, &
               kappaAbs=kappaAbsReal,kappaSca=kappaScaReal, grid=grid)


          call intersectCubeAMR2D(grid, octVec, uHat, tval)

          nTau = nTau + 1
          if (nTau > maxTau) then
             write(*,*) "maxTau exceeded in integratepathamr2d"
             nTau = maxTau
             escaped = .true.
          endif
          if (thisOctal%inflow(subcell)) then
             tauSca(nTau) = tauSca(nTau-1) + tval * kappaScaReal
             tauAbs(nTau) = tauAbs(nTau-1) + tval * kappaAbsReal
          else
             tauSca(nTau) = tauSca(nTau-1) + 1.e-20
             tauAbs(nTau) = tauAbs(nTau-1) + 1.e-20
          endif
          lambda(nTau) = lambda(nTau-1) + tval
          rVec = rVec + tval * uHat
          oldOctal => thisOctal
       endif
    end do
    tauExt(1:nTau) = tauSca(1:nTau) + tauAbs(1:nTau)
  end subroutine integratePathAMR2D


  subroutine intersectCubeAMR(grid, posVec, direction, tval)
   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTALVECTOR), intent(in) :: posVec
   type(OCTALVECTOR), intent(in) :: direction
   real(kind=octalkind), intent(out) :: tval
   !
   type(OCTALVECTOR) :: norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(OCTALVECTOR) :: subcen, point
   integer :: subcell
   
   real(kind=octalKind) :: t(6),denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)


   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.

   norm(1) = OCTALVECTOR(1., 0., 0.)
   norm(2) = OCTALVECTOR(0., 1., 0.)
   norm(3) = OCTALVECTOR(0., 0., 1.)
   norm(4) = OCTALVECTOR(-1., 0., 0.)
   norm(5) = OCTALVECTOR(0., -1., 0.)
   norm(6) = OCTALVECTOR(0., 0., -1.)

   p3(1) = OCTALVECTOR(subcen%x+thisOctal%subcellsize/2., subcen%y, subcen%z)
   p3(2) = OCTALVECTOR(subcen%x, subcen%y+thisOctal%subcellsize/2. ,subcen%z)
   p3(3) = OCTALVECTOR(subcen%x,subcen%y,subcen%z+thisOctal%subcellsize/2.)
   p3(4) = OCTALVECTOR(subcen%x-thisOctal%subcellsize/2., subcen%y,  subcen%z)
   p3(5) = OCTALVECTOR(subcen%x,subcen%y-thisOctal%subcellsize/2., subcen%z)
   p3(6) = OCTALVECTOR(subcen%x,subcen%y,subcen%z-thisOctal%subcellsize/2.)

   thisOk = .true.

   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.
      endif
      if (t(i) < 0.) thisOk(i) = .false.
!      if (denom > 0.) thisOK(i) = .false.
   enddo




  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo

  if (j == 0) ok = .false.
   
  if (.not.ok) then
     write(*,*) "Error: j=0 (no intersection???) in intersectCubeAMR. "
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     do;enddo
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble(thisOctal%subCellSize/1000.))


  if (tval == 0.) then
     write(*,*) posVec
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > sqrt(3.)*thisOctal%subcellsize) then
!     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif

end subroutine intersectCubeAMR



!
! Performs various optical depth tests.
!
! Check the optical depth in diffrent direction from the center of the model
! space.  The results will be written to a standard output and files.
!
subroutine test_optical_depth(gridUsesAMR, amrGridCentre, sphericityTest,  &
     dir_obs, wavelength,  lambda0, grid, opaqueCore, lamStart, lamEnd,   &
     ThinLine, nUpper, nLower, sampleFreq, useinterp, Rstar, coolStarPosition)

  implicit none
  
  logical, intent(in)           :: gridUsesAMR    !
  type(OCTALVECTOR), intent(in) :: amrGridCentre  ! central coordinates of grid
  logical, intent(in)           :: sphericityTest !

  type(vector), intent(in)  :: dir_obs            ! direction to the observer (unit vector)
  real, intent(in)          :: wavelength         ! the wavelength 
  real, intent(in)          :: lambda0            ! rest wavelength of line
  type(GRIDTYPE), intent(in):: grid               ! the opacity grid
  !

  logical, intent(in)       :: opaqueCore         ! is the core opaque
  real, intent(in)          :: lamStart, lamEnd
  !
  logical, intent(in)       :: thinLine           ! ignore line absorption of continuum
  integer, intent(in)       :: nUpper, nLower
  real, intent(in)          :: sampleFreq         ! max. samples per grid cell
  !
  logical, intent(in)       :: useinterp
  real, intent(in)          :: Rstar
  type(VECTOR), intent(in)  :: coolStarPosition

  integer, parameter :: maxTau  = 1000
  integer, parameter :: nLambda = 1000
  ! local variables
  real :: lambda(maxTau)    ! path distance array
  real :: tauExt(maxTau)    ! optical depth
  real :: tauAbs(maxTau)    ! optical depth
  real :: tauSca(maxTau)    ! optical depth
  real :: tauCont(maxTau,nLambda)
  !

  integer  :: nTau        ! size of optical depth arrays
  real     :: escProb     ! the escape probability
  logical  :: hitcore     ! has the photon hit the core
  real     :: localTau
  integer  :: error       ! error code returned
  integer  :: i 
  type(vector) :: tempVec
  type(vector) :: zerovec
  real :: junk
  type(OCTALVECTOR) :: octVec, position
  logical :: contPhoton = .true.

  integer :: ntest = 301
  real :: theta
  real(kind=octalkind) :: x1, x2
  integer, parameter :: UNOUT1 = 23
  integer, parameter :: UNOUT2 = 24
  
  real(kind=octalkind) :: R
  



  
  write(*,'(a,f7.1,a)') "Cross-sections at ",wavelength, " angstroms"
  write(*,'(a)') "------------------------------------------"
  write(*,*) " "

  zeroVec = VECTOR(0.,0.,0.)
  

  if (gridUsesAMR) then
     if (grid%octreeRoot%threed) then
        call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.),  &
             amrGridCentre,  &
             OCTALVECTOR(1.,0.,0.), grid, lambda, tauExt, tauAbs, &
             tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,Error)
     else
        call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.),  &
             amrGridCentre,  &
             OCTALVECTOR(1.,0.,0.), grid, lambda, tauExt, tauAbs, &
             tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,Error)
     endif
     if (error < 0) then
        write(*,*) '   Error encountered in cross-sections!!! (error = ',error,')'
     end if
  else
     call integratePath(wavelength,  lambda0, VECTOR(1.,1.,1.), zeroVec, &
          VECTOR(1.,0.,0.), grid, lambda, tauExt, tauAbs, &
          tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
          lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., rStar, &
          coolStarPosition,.false., nUpper, nLower, 0., 0., 0., junk, useInterp)
  end if
  
  if (nTau > 2) then 
     write(*,'(a,1pe10.3)') "Optical depth in x-axis from centre: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth in x-axis from centre: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth in x-axis from centre: ",tauSca(ntau)
     write(*,'(a,i4)') "Number of samples: ",nTau
     write(*,*) " "
  end if

  
  if (gridUsesAMR) then
     if (grid%octreeRoot%threed) then
        call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.), &
             amrGridCentre, &
             OCTALVECTOR(0.,1.,0.), grid, lambda, tauExt, tauAbs, &
             tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,error)
     else
        call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.), &
             amrGridCentre, &
             OCTALVECTOR(0.,1.,0.), grid, lambda, tauExt, tauAbs, &
             tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,error)
     endif
        if (error < 0) then
        write(*,*) '   Error encountered in cross-sections!!! (error = ',error,')'
        end if
     else
        call integratePath(wavelength, lambda0, VECTOR(1.,1.,1.), zeroVec, &
             VECTOR(0.,1.,0.), grid, lambda, tauExt, tauAbs, &
             tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., rStar, &
             coolStarPosition,.false., nUpper, nLower, 0., 0., 0., junk, useInterp)
     end if
     
     if (nTau > 2) then 
        write(*,'(a,1pe10.3)') "Optical depth in y-axis from centre: ",tauExt(ntau)
        write(*,'(a,1pe10.3)') "Absorption depth in y-axis from centre: ",tauAbs(ntau)
        write(*,'(a,1pe10.3)') "Scattering depth in y-axis from centre: ",tauSca(ntau)
        write(*,'(a,i4)') "Number of samples: ",nTau
        write(*,*) " "
     end if

     if (gridUsesAMR) then
        if (grid%octreeRoot%threed) then
           call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.), &
                amrGridCentre, &
                OCTALVECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
                tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
                lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
                .false., nUpper, nLower, 0., 0., 0., junk,&
                sampleFreq,error)
        else
           call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.), &
                amrGridCentre, &
                OCTALVECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
                tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
                lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
                .false., nUpper, nLower, 0., 0., 0., junk,&
                sampleFreq,error)
        endif

        if (error < 0) then
           write(*,*) '   Error encountered in cross-sections!!! (error = ',error,')'
        end if
     else
        call integratePath(wavelength,  lambda0, VECTOR(1.,1.,1.), zeroVec, &
             VECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
             tauSca, maxTau, nTau, opaqueCore, escProb, contPhoton , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., rStar, &
             coolStarPosition,.false., nUpper, nLower, 0., 0., 0., junk, useInterp)
     end if
     
     if (nTau > 2) then 
        write(*,'(a,1pe10.3)') "Optical depth in z-axis from centre: ",tauExt(ntau)
        write(*,'(a,1pe10.3)') "Absorption depth in z-axis from centre: ",tauAbs(ntau)
        write(*,'(a,1pe10.3)') "Scattering depth in z-axis from centre: ",tauSca(ntau)
        write(*,'(a,i4)') "Number of samples: ",nTau
        write(*,*) " "
     end if
     
     
     if (gridUsesAMR) then
        if (grid%octreeRoot%threed) then
           call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.), &
                amrGridCentre, s2o(dir_obs), grid, lambda, &
                tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
                contPhoton, lamStart, lamEnd, nLambda, tauCont, &
                hitCore, thinLine,.false.,  &
                .false., nUpper, nLower, 0., 0., 0., junk,sampleFreq,error)
        else
           call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.), &
                amrGridCentre, s2o(dir_obs), grid, lambda, &
                tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
                contPhoton, lamStart, lamEnd, nLambda, tauCont, &
                hitCore, thinLine,.false.,  &
                .false., nUpper, nLower, 0., 0., 0., junk,sampleFreq,error)
        endif
        if (error < 0) then
           write(*,*) '   Error encountered in test towards observer!!! (error = ',error,')'
        end if
!          open(20,file="tau.dat", status="unknown",form="formatted")
!          do i = 1, nTau
!             write(20,*) lambda(i)/grid%rInner, tauExt(i)
!          enddo
!          close(20)

     else
        call integratePath(wavelength,  lambda0, VECTOR(1.,1.,1.), &
             zeroVec, dir_obs, grid, lambda, &
             tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
             contPhoton, lamStart, lamEnd, nLambda, tauCont, &
             hitCore, thinLine,.false., rStar, coolStarPosition, &
             .false., nUpper, nLower, 0., 0., 0., junk, useInterp)
     end if
     
     write(*,'(a,1pe10.3)') "Optical depth to observer: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth to observer: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth to observer: ",tauSca(ntau)
     write(*,'(a,i4)') "Number of samples: ",nTau
     write(*,*) " "


     if (sphericityTest) then
        write(*,'(a)') "Sphericity test - 100 random directions:"
        do i = 1, 100
           tempVec = randomUnitVector()
           if (gridUsesAMR) then
              call IntegratePathAMR(lambda0,  lambda0, OCTALVECTOR(1.,1.,1.), &
                   amrGridCentre, s2o(tempVec), grid, lambda, &
                   tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
                   contPhoton, lamStart, lamEnd, nLambda, tauCont, &
                   hitCore, thinLine,.false., &
                   .false., nUpper, nLower, 0., 0., 0., junk, sampleFreq,error)
              if (error < 0) then
                 write(*,*) 'Error handling not implemented for ''sphericityTest''!!! (error = ',error,')'
              end if
           else
              call integratePath(lambda0,  lambda0, VECTOR(1.,1.,1.), &
                   zeroVec, tempVec, grid, lambda, &
                   tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
                   contPhoton, lamStart, lamEnd, nLambda, tauCont, &
                   hitCore, thinLine,.false., rStar, coolStarPosition, &
                   .false., nUpper, nLower, 0., 0., 0., junk,useInterp)
           end if
           write(*,'(a,1pe10.3,1pe10.3,1pe10.3)') "Optical depths: ",tauExt(nTau),tauAbs(nTau),tauSca(nTau)
        enddo
     endif


     !
     ! Compute the optical depth on x-y and z-x planes
     !  
     
     R = 1.0d-3   ! converting it to octal...
     if (gridUsesAMR) then


        ! 
        ! X-Y plane 
        !
        open(unit=UNOUT1, file = 'tau_xy_plane.dat', status ='replace')
        open(unit=UNOUT2, file = 'taul_xy_plane.dat', status ='replace')
        
        write(UNOUT1, '(a)') '#    phi [deg]  -- tau(total) -- tau(abs) -- tau(scat) -- ntau'
        write(UNOUT2, '(a)') '#    phi [deg]  -- taul(total) -- taul(abs) -- taul(scat) -- ntau'

        ! direction of the beam...
        do i = 1, ntest
           theta = 2.0d0*Pi*real(i-1)/real(ntest-1) + 0.01
           x1 = cos(theta); x2 = sin(theta)
           octVec = OCTALVECTOR(x1, x2, 0.0)
!           call Normalize(octVec)  ! just in case ..

           ! position of emission
           position = (octVec*R) + amrGridCentre
           
           ! continuum
           if (grid%octreeRoot%threed) then
              call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.),  &
                   position, octVec, grid, lambda, tauExt, tauAbs, &
                   tauSca, maxTau, nTau, opaqueCore, escProb, .true. , &
                   lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
                   .false., nUpper, nLower, 0., 0., 0., junk,&
                   sampleFreq,Error)
           else
              call IntegratePathAMR2D(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.),  &
                   position, octVec, grid, lambda, tauExt, tauAbs, &
                   tauSca, maxTau, nTau, opaqueCore, escProb, .true. , &
                   lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
                   .false., nUpper, nLower, 0., 0., 0., junk,&
                   sampleFreq,Error)
           endif
           
           write(UNOUT1, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau

           ! line 
           if (grid%octreeroot%threed) then
              call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.),  &
                   position, octVec, grid, lambda, tauExt, tauAbs, &
                   tauSca, maxTau, nTau, opaqueCore, escProb, .false. , &
                   lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
                   .false., nUpper, nLower, 0., 0., 0., junk,&
                   sampleFreq,Error)
           
              write(UNOUT2, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau
           endif

        end do

        close(UNOUT1)
        close(UNOUT2)


        ! 
        ! Z-X plane 
        !
        open(unit=UNOUT1, file = 'tau_zx_plane.dat', status ='replace')
        open(unit=UNOUT2, file = 'taul_zx_plane.dat', status ='replace')
        
        write(UNOUT1, '(a)') '#    phi [deg]  -- tau(total) -- tau(abs) -- tau(scat) -- ntau'        
        write(UNOUT2, '(a)') '#    phi [deg]  -- taul(total) -- taul(abs) -- taul(scat) -- ntau'

        ! direction of the beam...
        do i = 1, ntest
           theta = 2.0d0*Pi*real(i-1)/real(ntest-1) + 0.01
           x1 = cos(theta); x2 = sin(theta)
           octVec = OCTALVECTOR(x2, 0.0, x1)
 !          call Normalize(octVec)  ! just in case ..

           ! position of emission
           position = (octVec*R)  + amrGridCentre

           ! continuum
           if (grid%octreeRoot%threed) then
              call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.),  &
                   position, octVec, grid, lambda, tauExt, tauAbs, &
                   tauSca, maxTau, nTau, opaqueCore, escProb, .true. , &
                   lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
                   .false., nUpper, nLower, 0., 0., 0., junk,&
                   sampleFreq,Error)
           else
              call IntegratePathAMR2D(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.),  &
                   position, octVec, grid, lambda, tauExt, tauAbs, &
                   tauSca, maxTau, nTau, opaqueCore, escProb, .true. , &
                   lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
                   .false., nUpper, nLower, 0., 0., 0., junk,&
                   sampleFreq,Error)
           endif

           
           write(UNOUT1, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau

           if (grid%octreeroot%threed) then
              ! line 
              call IntegratePathAMR(wavelength,  lambda0, OCTALVECTOR(1.,1.,1.),  &
                   position, octVec, grid, lambda, tauExt, tauAbs, &
                   tauSca, maxTau, nTau, opaqueCore, escProb, .false. , &
                   lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, .false., &
                   .false., nUpper, nLower, 0., 0., 0., junk,&
                   sampleFreq,Error)
           
              write(UNOUT2, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau
           endif

        end do


        close(UNOUT1)
        close(UNOUT2)
        

       
     end if

     
     
   end subroutine test_optical_depth




subroutine integratePathAMRStark(wavelength,  lambda0, vVec, aVec, uHat, Grid, &
     lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb,&
     contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, &
     redRegion, usePops, mLevel, nLevel, &
     fStrength, gM, gN, localTau,sampleFreq,error)

  ! should we add the 'interp' argument and implement it?

  ! error codes are:  -10: too few samples made (the code should be fixed to 
  !                          stop this being a problem)
  !                   -20: the photon passed too close to a cell boundary and 
  !                          returnSamples decided to abort it 

  ! this routine does not employ the sobolev approximation but includes
  ! pressure broadening according to the algorithm of Muzerolle et al 2001, ApJ, 550, 944

  use constants_mod

  implicit none

  real, intent(in)          :: wavelength             ! the wavelength 
  real, intent(in)          :: lambda0                ! rest wavelength of line
  type(VECTOR), intent(in)  :: vVec                   ! velocity vector
  type(VECTOR), intent(in)  :: aVec                   ! starting position vector
  type(VECTOR), intent(in)  :: uHat                   ! direction
  type(GRIDTYPE), intent(in):: grid                   ! the opacity grid
  integer, intent(in)       :: maxTau 
  real, dimension(:), intent(out) :: lambda           ! path distance array
  real, dimension(:), intent(out) :: tauExt           ! optical depth
  real, dimension(:), intent(out) :: tauAbs           ! optical depth
  real, dimension(:), intent(out) :: tauSca           ! optical depth
  integer, intent(out)      :: nTau                   ! size of optical depth arrays
  logical, intent(in)       :: opaqueCore             ! is the core opaque
  real, intent(out)         :: escProb                ! the escape probability
  real, intent(in)          :: lamStart, lamEnd
  integer, intent(in)       :: nLambda
  real, intent(out)         :: tauCont(maxTau, nLambda)
  logical, intent(out)      :: hitcore                ! has the photon hit the core
  logical, intent(in)       :: thinLine               ! ignore line absorption of continuum
  logical, intent(in)       :: redRegion              ! use raman scattering red region opacities
  logical, intent(in)       :: usePops
  integer, intent(in)       :: mLevel, nLevel
  real, intent(in)          :: fStrength, gM, gN
  real, intent(out)         :: localTau
  real, intent(in)          :: sampleFreq             ! max. samples per grid cell
  integer, intent(out)      :: error                  ! error code returned


  real                      :: rho(1:maxTau)          ! density
  type(octalVector)         :: aVecOctal              ! octalVector version of 'aVec'
  type(octalVector)         :: uHatOctal              ! octalVector version of 'uHat'
  type(VECTOR)              :: rVec                   ! position vector
  logical                   :: contPhoton             ! is this a continuum photon?
  real, dimension(1:maxTau) :: chiLine                ! line opacities
  real, dimension(1:maxTau,1:grid%maxLevels) :: levelPop  ! stateq level pops
  real, dimension(1:maxTau) :: tauSob, tauSob2        ! Sobolev optical depths
  real, dimension(1:maxTau) :: escProbArray           ! escape probabilities
  type(vector), dimension(1:maxTau) :: velocity       ! 
  real, dimension(1:maxTau) :: velocityDeriv          ! directional derivative

  real :: nu, nu2                                        ! frequencies

  integer :: iStart                                 ! starting index

  real :: tauSobScalar,tauSobScalar1,tauSobScalar2  ! Sobolev optical depths

  real :: thisVel
  integer :: ilambda                                ! wavelength index


  integer, parameter :: maxLambda = 2000             ! max size of tau arrays
  real :: dlambda(maxTau)                           ! distance increment array

  real :: projVel(maxTau)
  real :: kabs(maxtau), ksca(maxtau)                ! opacities
  real :: temperature(maxtau)
  real :: newLambda(maxtau)
  real :: Ne(maxtau)
  real :: N_HI(maxtau)
  real :: x                                         ! multipliers

  real :: chil                                      ! line opacity

  real :: lambda2, thisVel2

  integer :: i, j                   ! counters


  real :: deltaNu
  real :: meanMoleMass
  real :: a, y
  real :: Hay(0:0), dv(0:0)

  integer :: newNtau

  ! initialize variables

  hitcore = .false.
  rVec = aVec
  escProb = 1.
  error = 0

  ! locate this wavelength in the grid of wavelengths

  if (grid%flatspec.or.(grid%doRaman)) then
     iLambda = 1
  else
     call locate(grid%lamArray, grid%nLambda, wavelength, iLambda)
  endif

  if (.not. grid%adaptive) then
     print *, 'integratePathAMR called on a non-adaptive grid!'
     stop
  end if

  nTau = 0

  ! convert some variables to the appropriate kind for startReturnSamples
  aVecOctal = aVec
  uHatOctal = uHat


  ! some geometries produce photons at the origin

  if (modulus(aVecOctal) < 0.01 * grid%halfSmallestSubcell) then
     if (grid%rCore < 1.e-10) then
        aVecOctal = (1.0001_oc * grid%halfSmallestSubcell) * uHatOctal
     else 
        aVecOctal = (1.0001_oc * grid%rCore) * uHatOctal
     end if
  end if

  ! find the project velocity

  thisVel = (lambda0 - wavelength)/lambda0
  thisVel = thisVel  + (uHat .dot. vVec)  

  escProb = 1.

  CALL startReturnSamples (aVecOctal,uHatOctal,grid,sampleFreq,nTau,       &
       maxTau,opaqueCore,hitcore,usePops,iLambda,error,&
       lambda,kappaAbs=kAbs,kappaSca=kSca,velocity=velocity,&
       velocityDeriv=velocityDeriv,chiLine=chiLine,    &
       levelPop=levelPop,rho=rho)

  if (nTau <= 2) then
     !print *, 'The code does not yet include routines for handling ',&
     !         'the case where nTau <= 2 (in integratePathAMR)'
     error = -10
     return
  end if


  ! projected velocities 

  forall (i = 1:nTau)
     projVel(i) = velocity(i) .dot. uHat
  end forall

  dlambda(1:nTau-1) = lambda(2:nTau) - lambda(1:nTau-1)

  ! now compute the optical depths from the opacities

  tauExt(1:2) = 0.
  tauAbs(1:2) = 0.
  tauSca(1:2) = 0.

  do i = 2, nTau, 1
     tauSca(i) = tauSca(i-1) + dlambda(i-1)*0.5*(ksca(i-1)+ksca(i))
     tauAbs(i) = tauAbs(i-1) + dlambda(i-1)*0.5*(kabs(i-1)+kabs(i))
     if ((ksca(i) < 0.).or.(kabs(i)<0.)) then
        write(*,*) "negative opacity"
        do
        enddo
     endif
        
  enddo


  iStart = 1

  if (usePops) then
     chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
     chiLine(1:nTau) = 1.e10 * chil * (levelPop(1:nTau,mLevel)&
          -( ( gM / gN) * (levelPop(1:nTau,nLevel))))
  endif


! Now we need to resample 

  call resampleRay(lambda, nTau, projVel, newLambda, newNTau)

! Now interpolate on to newly sampled ray

  call linearResample(lambda, projVel, nTAu, newLambda, newNtau)
  call linearResample(lambda, tauExt, nTAu, newLambda, newNtau)
  call linearResample(lambda, tauSca, nTAu, newLambda, newNtau)
  call linearResample(lambda, chiline, nTAu, newLambda, newNtau)
  call linearResample(lambda, temperature, nTAu, newLambda, newNtau)

! finish off

  nTau = newNtau
  lambda(1:nTau) = newLambda(1:nTau)

  dlambda(1:nTau-1) = lambda(2:nTau) - lambda(1:nTau-1)

                    
! line photons

  do i = 2, nTau
     deltaNu = (cSpeed/lambda0*angstromToCm)/cSpeed * sqrt(2.*kErg*temperature(i)/meanMoleMass) !eq 7
     a = bigLambda(N_HI(i), temperature(i), Ne(i)) / (fourPi * deltaNu)
     dv(0) = abs(projVel(i)-thisVel)
     call voigt(1, dv, a, Hay)
     chil = chiLine(i) / (sqrt(pi)*deltaNu) * Hay(0) ! equation 5
     tauAbs(i) = tauAbs(i) + chil * dlambda(i)
  enddo

  tauExt(1:nTau) = tauSca(1:nTau) + tauAbs(1:nTau)


! continuum photons


  do j = 1, nLambda, 1
     
     ! compute projected velocity of this bin
     
     thisVel = real(j-1)/real(nLambda-1)
     thisVel = lamStart + thisVel*(lamEnd-lamStart)
     thisVel = (thisVel-lambda0)/lambda0
     thisVel = -thisVel  + (uHat .dot. vVec) - (wavelength - lambda0)/lambda0
     

     do i = 2, nTau
        deltaNu = (cSpeed/lambda0*angstromToCm)/cSpeed * sqrt(2.*kErg*temperature(i)/meanMoleMass)
        a = bigLambda(N_HI(i), temperature(i), Ne(i)) / (fourPi * deltaNu)
        dv(0) = abs(projVel(i)-thisVel) / sqrt(2.*kErg*temperature(i)/meanMoleMass)
        call voigt(1, dv, a, Hay)
        chil = chiLine(i) / (sqrt(pi)*deltaNu) * Hay(0)
        tauCont(i:nTau,j) = tauCont(i:nTau,j) + chil * dlambda(i)
     enddo
                    

  enddo   ! over wavelength array



end subroutine integratePathAMRStark



  subroutine intersectCubeAMR2D(grid, posVec, direction, tval)

! this is to find a cell intersection for a 2D AMR grid
! which is essentially a 2D-grid that projects into
! a 3D cylindrical coordinate system


   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTALVECTOR), intent(in) :: posVec
   type(OCTALVECTOR), intent(in) :: direction
   real(kind=octalkind), intent(out) :: tval
   !
   type(OCTALVECTOR) :: norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(OCTALVECTOR) :: subcen, point
   real :: dx, dz ! directions in cylindrical coords
   integer :: subcell
   real(kind=doubleKind) :: compX, compZ,currentX,currentZ
   real(kind=doubleKind) :: distToZBoundary, distToXboundary
   real(kind=octalKind) :: r1,r2,d,cosmu,x1,x2,distTor1,distTor2, theta, mu
   integer :: i,j
   logical :: ok, thisOk(6)
   type(OCTALVECTOR) :: xHat, zHAt

   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.


   r1 = subcen%x - thisOctal%subcellSize/2.
   r2 = subcen%x + thisOctal%subcellSize/2.
   d = sqrt(point%x**2+point%y**2)
   xHat = VECTOR(point%x, point%y,0.d0)
   call normalize(xHat)

   cosmu =((-1.d0)*xHat).dot.direction
   call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
   if (.not.ok) then
      write(*,*) "Quad solver failed in intersectcubeamr2d"
      do;enddo
   endif
   distTor2 = max(x1,x2)

   theta = asin(max(-1.d0,min(1.d0,r1 / d)))
   cosmu = xHat.dot.direction
   mu = acos(max(-1.d0,min(1.d0,cosmu)))
   distTor1 = 1.e30
   if (mu  < theta ) then
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in intersectcubeamr2d"
         do;enddo
      endif
      distTor1 = max(x1,x2)
   endif
         
   distToXboundary = min(distTor1, distTor2)


   zHat = VECTOR(0.d0, 0.d0, 1.d0)
   compZ = zHat.dot.direction
   currentZ = point%z

   if (compZ /= 0.d0 ) then
      if (compZ > 0.d0) then
         distToZboundary = (subcen%z + thisOctal%subcellsize/2. - currentZ ) / compZ
      else
         distToZboundary = abs((subcen%z - thisOctal%subcellsize/2. - currentZ ) / compZ)
      endif
   else
      disttoZboundary = 1.e30
   endif

   tVal = min(distToZboundary, distToXboundary)
   tval = max(tval * 1.0001d0,dble(thisOctal%subCellSize/10000.d0))
   if (tVal > 1.e29) then
      write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "x,z",currentX,currentZ
   endif
   if (tval < 0.) then
      write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "x,z",currentX,currentZ
   endif

  end subroutine intersectCubeAMR2D

end module path_integral
