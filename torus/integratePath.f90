module path_integral

  use vector_mod             ! vector maths 
  use gridtype_mod           ! opacity grid
  use grid_mod               ! opacity grid routines
  use math_mod               ! misc maths routines
  use photon_mod             ! the photon module
  use constants_mod          ! physical constants
  use amr_mod
  
  implicit none

  public::  &
       integratePath, &
       test_optical_depth, &
       intersectCubeAMR
       
  private:: &
       integratePathCaresian, &
       integratePathAMR, &
       integratePathVoigtProf
       
!
! This subroutine performs an integration through the grid in order
! to compute the scattering and absorption optical depths.
!
  contains


    ! 
    ! This is the interface routine for three different 
    ! integratePath subroutines defined in this module.
    !  Use this to compunicate outside of this module!
    !
    subroutine integratePath(gridUsesAMR, VoigtProf, &
         wavelength,  lambda0, vVec, aVec, uHat, Grid, &
         lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb,&
         contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, &
         lineResAbs, redRegion, usePops, mLevel, nLevel, &
         fStrength, gM, gN, localTau,sampleFreq,error, interp, &
         rStar, coolStarPosition, nSource, source)
      
      implicit none
      integer :: nSource
      type(SOURCETYPE) :: source(:)
      logical, intent(in) :: gridUsesAMR                  ! T if AMR grid 
      logical, intent(in) :: VoigtProf              ! T if use voigt profile
      real, intent(in)          :: wavelength             ! the wavelength 
      real, intent(in)          :: lambda0                ! rest wavelength of line
      type(OCTALVECTOR), intent(in)  :: vVec              ! velocity vector
      type(OCTALVECTOR), intent(in)  :: aVec              ! starting position vector
      type(OCTALVECTOR), intent(in)  :: uHat              ! direction
      type(GRIDTYPE), intent(in)     :: grid              ! the opacity grid
      logical, intent(in)            :: contPhoton        ! is this a continuum photon?
      integer, intent(in)            :: maxTau 
      real, dimension(:), intent(out) :: lambda           ! path distance array
      real, dimension(:), intent(out) :: tauExt           ! optical depth
      real, dimension(:), intent(out) :: tauAbs           ! optical depth
      real, dimension(:), intent(out) :: tauSca           ! optical depth
      real, intent(inout)       :: linePhotonAlbedo(1:maxTau) ! the line photon albedo along the ray
      integer, intent(out)      :: nTau                   ! size of optical depth arrays
      logical, intent(in)       :: thin_disc_on      ! T to include thin disc
      logical, intent(in)       :: opaqueCore             ! is the core opaque
      real, intent(out)         :: escProb                ! the escape probability
      real, intent(in)          :: lamStart, lamEnd
      integer, intent(in)       :: nLambda
      real,intent(inout),dimension(maxTau, nLambda) :: tauCont
      logical, intent(out)      :: hitcore                ! has the photon hit the core
      logical, intent(in)       :: thinLine               ! ignore line absorption of continuum
      logical, intent(in)       :: lineResAbs             ! T if you want to include absorption
      !                                                   !   of line by line resonance zones.
      logical, intent(in)       :: redRegion              ! use raman scattering red region opacities
      logical, intent(in)       :: usePops
      integer, intent(in)       :: mLevel, nLevel
      real, intent(in)          :: fStrength, gM, gN
      real, intent(out)         :: localTau
      real, intent(in)          :: sampleFreq             ! max. samples per grid cell
      integer, intent(out)      :: error                  ! error code returned
      logical, intent(in)       :: interp                 ! ????
      type(VECTOR), optional :: coolStarPosition                  ! pos of cool giant star
      real, optional :: rStar


      if (gridUsesAMR) then
         if (VoigtProf) then
            call integratePathVoigtProf(wavelength,  lambda0, o2s(vVec), o2s(aVec), o2s(uHat), Grid, &
                 lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb,&
                 contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, &
                 redRegion, usePops, mLevel, nLevel, &
                 fStrength, gM, gN, sampleFreq,error)
         else
            call integratePathAMR(wavelength,  lambda0, vVec, aVec, uHat, Grid, &
                 lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb,&
                 contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, &
                 lineResAbs, redRegion, usePops, mLevel, nLevel, &
                 fStrength, gM, gN, localTau,sampleFreq,error, nSource, source)
         end if

      else  ! grid uses a Cartisian grid
         if (.not. PRESENT(rstar) ) then
            write(*,*) "Error: rstar must be present. [integratePath]"
            stop
         end if
         if (.not. PRESENT(coolStarPosition) ) then
            write(*,*) "Error: coolStarPosition must be present. [integratePath]"
            stop
         end if
         call integratePathCaresian(wavelength,  lambda0, o2s(vVec), o2s(aVec), o2s(uHat), Grid,  lambda, &
              tauExt, tauAbs, tauSca,  linePhotonAlbedo, maxTau, nTau, opaqueCore, escProb, contPhoton, &
              lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, redRegion, rStar, &
              coolStarPosition, usePops, mLevel, nLevel, fStrength, gM, gN, localTau, interp)
      end if

    end subroutine integratePath


!
!
!
!
         
subroutine integratePathCaresian(wavelength,  lambda0, vVec, aVec, uHat, Grid,  lambda, &
     tauExt, tauAbs, tauSca, linePhotonalbedo, maxTau, nTau, opaqueCore, escProb, contPhoton, &
     lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, redRegion, rStar, &
     coolStarPosition, usePops, mLevel, nLevel, fStrength, gM, gN, localTau, interp)


  type(GRIDTYPE) :: grid                            ! the opacity grid

  real :: localTau
  logical :: usePops
  integer :: mLevel, nLevel
  real :: fStrength, gM, GN
  real :: rStar
  logical :: interp

  type(VECTOR), intent(in)  :: vVec                   ! velocity vector
  type(VECTOR), intent(in)  :: aVec                   ! starting position vector
  type(VECTOR), intent(in)  :: uHat                   ! direction
  type(VECTOR) :: rVec                              ! position vector
  type(VECTOR) :: rHat                              ! unit radius vector
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

  real:: tauSob, tauSob1, tauSob2                  ! Sobolev optical depths

  real :: thisVel
  integer :: ilambda                                ! wavelength index




!  integer, parameter :: maxLambda = 200              ! max size of tau arrays
  integer :: nLambda
  real :: lambda(:),dlambda                         ! distance arrays
  real :: tauExt(:), tauAbs(:), tauSca(:)           ! optical depth
  real, intent(inout)       :: linePhotonAlbedo(1:maxTau) ! the line photon albedo along the ray

  real, dimension(maxTau, nLambda) :: tauCont


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
           linePhotonAlbedo(nTau) = kScaInterp/ (kScaInterp + kAbsInterp)
        else
           kabsInterp = interpGridkappaAbsRed(grid,i1,i2,i3,iLambda,t1,t2,t3) 
           kscaInterp = interpGridkappaScaRed(grid,i1,i2,i3,iLambda,t1,t2,t3) 

           tauSca(nTau) = tauSca(nTau-1) + kscaInterp * dlambda
           tauAbs(nTau) = tauAbs(nTau-1) + kabsInterp * dlambda
           tauExt(nTau) = tauAbs(ntau) + tauSca(ntau)
           linePhotonAlbedo(nTau) = kScaInterp/ (kScaInterp + kAbsInterp)
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

        if (tauSob < 0.01) then
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
              dlambda = min(modulus(fVec),dr)
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

end subroutine integratePathCaresian


subroutine integratePathAMR(wavelength,  lambda0, vVec, aVec, uHat, Grid, &
     lambda, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb,&
     contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs,&
     redRegion, usePops, mLevel, nLevel, &
     fStrength, gM, gN, localTau,sampleFreq,error, nSource, source)

  ! should we add the 'interp' argument and implement it?

  ! error codes are:  -10: too few samples made (the code should be fixed to 
  !                          stop this being a problem)
  !                   -20: the photon passed too close to a cell boundary and 
  !                          returnSamples decided to abort it 
 
  use constants_mod

  implicit none

  integer :: nSource
  type(SOURCETYPE) :: source(:)
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
  real, intent(inout)       :: linePhotonAlbedo(1:maxTau) ! the line photon albedo along the ray
  integer, intent(out)      :: nTau                   ! size of optical depth arrays
  logical, intent(in)       :: thin_disc_on           ! T to include thin disc
  logical, intent(in)       :: opaqueCore             ! is the core opaque
  real, intent(out)         :: escProb                ! the escape probability
  real, intent(in)          :: lamStart, lamEnd
  integer, intent(in)       :: nLambda
  real,intent(inout),dimension(maxtau, nLambda) :: tauCont
  logical, intent(out)      :: hitcore                ! has the photon hit the core
  logical, intent(in)       :: thinLine               ! ignore line absorption of continuum
  logical, intent(in)       :: lineResAbs       ! T if you want to include absorption
  !                                                   !   of line by line resonance zones.
  logical, intent(in)       :: redRegion              ! use raman scattering red region opacities
  logical, intent(in)       :: usePops
  integer, intent(in)       :: mLevel, nLevel
  real, intent(in)          :: fStrength, gM, gN
  real, intent(out)         :: localTau
  real, intent(in)          :: sampleFreq             ! max. samples per grid cell
  integer, intent(out)      :: error                  ! error code returned
  !
  ! WORK ARRAYS
  logical, save :: first_time = .true.
  real(double), allocatable, save :: rho(:)         ! density (size=maxTau)
  real, allocatable, save :: temperature(:) ! temperature (size=maxTau)
!  real(double), allocatable, save :: Ne(:)  ! electron density (size=maxTau)
!  real(double), allocatable, save :: etaLine(:)  ! line emissivity(size=maxTau)
!  real(double), allocatable, save :: etaCont(:)  ! continuum emissivity (size=maxTau)
  logical, allocatable, save :: inFlow(:)   ! inFlow flags (size=maxTau)   
  real(double), allocatable, save :: levelPop(:,:)  ! (size=maxTau x grid%maxLevels)
  type(vector), allocatable, save :: velocity(:) ! size=maxTau       ! 
  real, allocatable, save  :: chiLine(:)         ! line optical depth (size=maxTau)
  real, allocatable, save  :: tauSob(:), tauSob2(:) ! (size=maxTau) Sobolev optical depths
  real, allocatable, save  :: velocityDeriv(:)   ! directional derivative (size=maxTau)
  real, allocatable, save  :: dlambda(:)         ! distance increment array (size=maxTau)
  real(double), allocatable, save :: projVel(:)  ! (size=maxTau)
  real(double), allocatable, save  :: kabs(:), ksca(:)  ! (size=maxTau)
  real(double), allocatable, save :: dummy(:)  ! Work array (size=maxTau)
  !
  !
  type(octalVector)         :: aVecOctal              ! octalVector version of 'aVec'
  type(octalVector)         :: uHatOctal              ! octalVector version of 'uHat'
  type(OCTALVECTOR)         :: rVec                   ! position vector
  logical                   :: contPhoton             ! is this a continuum photon?

  real :: nu, nu2                                        ! frequencies

  integer :: iStart                                 ! starting index

  real :: tauSobScalar,tauSobScalar1,tauSobScalar2  ! Sobolev optical depths
  real(double) :: tauDouble

  real :: thisVel
  integer :: ilambda                                ! wavelength index


  real :: x                                         ! multipliers

  real :: chil                                      ! line opacity
  
  real :: lambda2, thisVel2

  integer :: i, j                   ! counters

  ! initialize variables

  hitcore = .false.
  rVec = aVec
  escProb = 1.
  error = 0


    ! Allocates memory for work arrays for the first time
    ! This should be faster than using automatic arrays which allocates
    ! and deallocates memory every single time.  These array are SAVED
    ! so they should be available next time as well.
    if (first_time) then
       first_time = .false.
       ALLOCATE(rho(maxTau))         
       ALLOCATE(temperature(maxTau)) 
       ALLOCATE(inFlow(maxTau))
       ALLOCATE(levelPop(maxTau,grid%maxLevels))  
       ALLOCATE(velocity(maxTau)) 
       ALLOCATE(chiLine(maxTau))
       ALLOCATE(tauSob(maxTau), tauSob2(maxTau)) 
       ALLOCATE(kabs(maxTau), ksca(maxTau))
       ALLOCATE(velocityDeriv(maxTau)) 
       ALLOCATE(dLambda(maxTau))
       ALLOCATE(projVel(maxTau))
       ALLOCATE(dummy(maxTau))
    end if




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
     CALL startReturnSamples2 (aVecOctal,uHatOctal,grid,sampleFreq,nTau,       &
                              maxTau,thin_disc_on,opaqueCore,hitcore,usePops,iLambda,error,&
                              lambda, nSource, source, kappaAbs=kAbs,kappaSca=kSca, velocity=velocity,&
                              velocityDeriv=velocityDeriv,chiLine=chiLine,    &
                              levelPop=levelPop,rho=rho, &
                              temperature=temperature, Ne=dummy, inflow=inflow, etaCont=dummy, etaLine=dummy)

     dlambda(1:nTau-1) = lambda(2:nTau) - lambda(1:nTau-1)        


     if (.not.redRegion) then
        ! first optical depths are all zero of course

        tauExt(1) = 0.
        tauAbs(1) = 0.
        tauSca(1) = 0.

        do i = 2, nTau, 1 
! RK CHANGED THE FOLLOWINGS
           dlambda(i) = lambda(i)-lambda(i-1)
           tauSca(i) = tauSca(i-1) + dlambda(i)*0.5*(ksca(i-1)+ksca(i))
           tauAbs(i) = tauAbs(i-1) + dlambda(i)*0.5*(kabs(i-1)+kabs(i))
!           tauSca(i) = tauSca(i-1) + dlambda(i-1)*ksca(i-1)
!           tauAbs(i) = tauAbs(i-1) + dlambda(i-1)*kabs(i-1)

           if ((ksca(i) < 0.).or.(kabs(i)<0.)) then
              write(*,*) "negative opacity"
              error = -70
              return
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

     ! find the projected velocity

     thisVel = (lambda0 - wavelength)/lambda0  
     thisVel = thisVel  + (uHat .dot. vVec)    ! (+ve when toward observer!)

     if (grid%resonanceLine) then
        lambda2 = grid%lambda2
        thisVel2 = (lambda2 - wavelength)/lambda2
        thisVel2 = thisVel2  + (uHat .dot. vVec) 
     endif


     escProb = 1.

     CALL startReturnSamples2 (aVecOctal,uHatOctal,grid,sampleFreq,nTau,       &
                              maxTau,thin_disc_on,opaqueCore,hitcore,usePops,iLambda,error,&
                              lambda,nSource, source, kappaAbs=kAbs,kappaSca=kSca,velocity=velocity,&
                              velocityDeriv=velocityDeriv,chiLine=chiLine,    &
                              levelPop=levelPop,rho=rho, &
                              temperature=temperature, Ne=dummy, inflow=inflow, etaCont=dummy, etaLine=dummy)


!     linePhotonAlbedo(1:nTau) = kSca(1:nTau) / (kSca(1:nTau)+kAbs(1:nTau)+chiLine(1:nTau))
     linePhotonAlbedo(1:nTau) = kSca(1:nTau) / (kSca(1:nTau)+kAbs(1:nTau))



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
           chil = 1.e10 * chil * (levelPop(1,mLevel)&
                          -( ( gM / gN) * (levelPop(1,nLevel))))
        else
          chil = chiLine(1)
        endif

        nu  = cSpeed / (lambda0*angstromtocm)
        tauDouble = chil  / (nu  * velocityDeriv(1))
        
        if (tauDouble < 0.01_db) then
          escProb = 1.0_db-tauDouble*0.5_db* &
            (1.0_db - tauDouble/3.0_db*(1.0_db - tauDouble*0.25_db*(1.0_db - 0.2_db*tauDouble)))
        else if (tauDouble < 15.0_db) then
          escProb = (1.0_db-exp(-tauDouble))/tauDouble
        else
          escProb = 1.0_db/tauDouble
        end if
        localTau = real(tauDouble)

     endif


     ! projected velocities 

!     if (contPhoton) then
!        forall (i = 1:nTau)
!           projVel(i) = velocity(i) .dot. uHat
!        end forall
!     else
!        projVel(1:nTau) = thisVel
!     endif

     forall (i = 1:nTau)
        projVel(i) = velocity(i) .dot. uHat  ! (+ve when moving toward!)
     end forall


     dlambda(1:nTau-1) = lambda(2:nTau) - lambda(1:nTau-1)  

     ! now compute the optical depths from the opacities

     tauAbs(1:2) = 0.
     tauSca(1:2) = 0.

     do i = 2, nTau, 1
        tauSca(i) = tauSca(i-1) + dlambda(i-1)*0.5*(ksca(i-1)+ksca(i))
        tauAbs(i) = tauAbs(i-1) + dlambda(i-1)*0.5*(kabs(i-1)+kabs(i))
!        tauSca(i) = tauSca(i-1) + dlambda(i-1)*ksca(i-1)
!        tauAbs(i) = tauAbs(i-1) + dlambda(i-1)*kabs(i-1)
           if ((ksca(i) < 0.).or.(kabs(i)<0.)) then
              write(*,*) "negative opacity"
              error = -70
              return
           endif

     enddo

     tauExt(1:nTau) = tauSca(1:nTau) + tauAbs(1:nTau)



     iStart = 1

     if (grid%lineEmission) then
        if (nTau > 0) then

           ! line photons are treated as individuals

           if (.not.contPhoton  .and. lineResAbs) then

              ! Setting up the Tau_sov 
              if (usePops) then
                 chil=((pi*eCharge**2)/(mElectron*cSpeed))*fStrength
                 chiLine(1:nTau) = 1.e10 * chil * (levelPop(1:nTau,mLevel)&
                      -( ( gM / gN) * (levelPop(1:nTau,nLevel))))
              endif              
              nu  = cSpeed / (lambda0*angstromtocm)
              tauSob(1:nTau) = chiLine(1:nTau)  / (nu  * velocityDeriv(1:nTau))
              
              ! for resonance line case
              if (grid%resonanceLine) then
                 nu2  = cSpeed / (lambda2*angstromtocm)
                 tauSob2(1:nTau) = chiLine(1:nTau)  / (nu2  * velocityDeriv(1:nTau))
              end if

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
                          write(*,*) "Error:: integratepathamr "
                          write(*,*) "   tauSobScalar < 0. .and.  abs(tauSobScalar) < 1.0e-3"
                          write(*,*) "tausobScalar",tausobScalar
                          write(*,*) "thisVel",thisVel*cSpeed/1.e5
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
!                    thisVel = real(j-1)/real(nLambda-1)
!                    thisVel = lamStart + thisVel*(lamEnd-lamStart)
!                    thisVel = (thisVel-lambda0)/lambda0
!                    thisVel = -thisVel  + (uHat .dot. vVec) - (wavelength - lambda0)/lambda0

                    thisVel = real(j-1)/real(nLambda-1)
                    thisVel = grid%lamArray(j)
                    thisVel = (lambda0-thisVel)/lambda0  + (uHat .dot. vVec)  &
                         - (wavelength - lambda0)/lambda0 ! (+ve if moving toward)
                    
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

! subroutine integratePathAMR2(wavelength,  lambda0, vVec, aVec, uHat, Grid, &
!      lambda, tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb,&
!      contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, &
!      redRegion, usePops, mLevel, nLevel, &
!      fStrength, gM, gN, sampleFreq,error)

!   ! should we add the 'interp' argument and implement it?

!   ! error codes are:  -10: too few samples made (the code should be fixed to 
!   !                          stop this being a problem)
!   !                   -20: the photon passed too close to a cell boundary and 
!   !                          returnSamples decided to abort it 
 
!   use constants_mod

!   implicit none

!   logical :: escaped
!   type(OCTALVECTOR) :: octVec
!   integer :: subcell
!   real(oct) :: tval

!   real, intent(in)          :: wavelength             ! the wavelength 
!   real, intent(in)          :: lambda0                ! rest wavelength of line
!   type(VECTOR), intent(in)  :: vVec                   ! velocity vector
!   type(VECTOR), intent(in)  :: aVec                   ! starting position vector
!   type(VECTOR), intent(in)  :: uHat                   ! direction
!   type(GRIDTYPE), intent(in):: grid                   ! the opacity grid
!   integer, intent(in)       :: maxTau
!   real,intent(out),dimension(:) :: lambda             ! path distance array
!   real,intent(out),dimension(:) :: tauExt             ! optical depth
!   real,intent(out),dimension(:) :: tauAbs             ! optical depth
!   real,intent(out),dimension(:) :: tauSca             ! optical depth
!   integer, intent(out)      :: nTau                   ! size of optical depth arrays
!   logical, intent(in)       :: opaqueCore             ! is the core opaque
!   real, intent(out)         :: escProb                ! the escape probability
!   real, intent(in)          :: lamStart, lamEnd
!   integer, intent(in)       :: nLambda
!   real,intent(inout),dimension(maxTau, nLambda) :: tauCont
!   logical, intent(out)      :: hitcore                ! has the photon hit the core
!   logical, intent(in)       :: thinLine               ! ignore line absorption of continuum
!   logical, intent(in)       :: redRegion              ! use raman scattering red region opacities
!   logical, intent(in)       :: usePops
!   integer, intent(in)       :: mLevel, nLevel
!   real, intent(in)          :: fStrength, gM, gN
! !  real, intent(out)         :: localTau
!   real, intent(in)          :: sampleFreq             ! max. samples per grid cell
!   integer, intent(out)      :: error                  ! error code returned


!   real                      :: rho(1:maxTau)          ! density
!   type(octalVector)         :: aVecOctal              ! octalVector version of 'aVec'
!   type(octalVector)         :: uHatOctal              ! octalVector version of 'uHat'
!   type(VECTOR)              :: rVec                   ! position vector
!   logical                   :: contPhoton             ! is this a continuum photon?
!   real, dimension(1:maxTau) :: chiLine                ! line opacities
!   real, dimension(1:maxTau,1:grid%maxLevels) :: levelPop  ! stateq level pops
!   real, dimension(1:maxTau) :: tauSob, tauSob2        ! Sobolev optical depths
! !  real, dimension(1:maxTau) :: escProbArray           ! escape probabilities
!   type(vector), dimension(1:maxTau) :: velocity       ! 
!   real, dimension(1:maxTau) :: velocityDeriv          ! directional derivative
  
!   real :: nu, nu2                                        ! frequencies

!   integer :: iStart                                 ! starting index

!   real :: tauSobScalar,tauSobScalar1,tauSobScalar2  ! Sobolev optical depths

!   real :: thisVel
!   integer :: ilambda                                ! wavelength index


!   integer, parameter :: maxLambda = 2000             ! max size of tau arrays
!   real :: dlambda(maxTau)                           ! distance increment array

!   real :: projVel(maxTau)
!   real :: kabs(maxtau), ksca(maxtau)                ! opacities
!   real :: x                                         ! multipliers

!   real :: chil                                      ! line opacity

!    type(OCTAL), pointer :: thisOctal
!    type(OCTAL),pointer :: oldOctal

!   real :: kappaScaReal, kappaAbsReal
  
!   real :: lambda2, thisVel2

!   integer :: i, j                   ! counters

!   ! initialize variables

!   hitcore = .false.
!   rVec = aVec
!   escProb = 1.
!   error = 0

!   ! locate this wavelength in the grid of wavelengths

!   if (grid%flatspec.or.(grid%doRaman)) then
!      iLambda = 1
!   else
!      call locate(grid%lamArray, grid%nLambda, wavelength, iLambda)
!   endif

!   if (.not. grid%adaptive) then
!      print *, 'integratePathAMR called on a non-adaptive grid!'
!      stop
!   end if

!   nTau = 1
!   tauSca(nTau) = 0.
!   tauAbs(nTau) = 0.
!   lambda(nTau) = 0.
!   escaped = .false.
!   !call amrGridValues(grid%octreeRoot, octVec, iLambda=iLambda, &
!   !     foundOctal=thisOctal, foundSubcell=subcell, &
!   !     kappaAbs=kappaAbsReal,kappaSca=kappaScaReal, grid=grid)
!   oldOctal => grid%octreeRoot
!     do while (.not.escaped) 
!        octVec = rVec
!        if (.not.inOctal(grid%octreeRoot, octVec)) then
!           escaped = .true.
!        endif
!        if (.not.escaped) then
!           call intersectCubeAMR(grid, s2o(rVec), s2o(uHat), tVal)
!           call amrGridValues(grid%octreeRoot, octVec, startOctal=oldOctal,iLambda=iLambda, &
!                foundOctal=thisOctal, foundSubcell=subcell, &
!                kappaAbs=kappaAbsReal,kappaSca=kappaScaReal, grid=grid)
!           nTau = nTau + 1
!           tauSca(nTau) = tauSca(nTau-1) + tval * kappaScaReal
!           tauAbs(nTau) = tauAbs(nTau-1) + tval * kappaAbsReal
!           lambda(nTau) = lambda(nTau-1) + tval
!           rVec = rVec + tval * uHat
!           oldOctal => thisOctal
!        endif
!     end do
!     tauExt(1:nTau) = tauSca(1:nTau) + tauAbs(1:nTau)
!   end subroutine integratePathAMR2


  !
  !
  !
  subroutine integratePathVoigtProf(wavelength,  lambda0, vVec, aVec, uHat, Grid, &
       L, tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb,&
       contPhoton, lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, &
       redRegion, usePops, mLevel, nLevel, &
       fStrength, gM, gN, sampleFreq,error)
    
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
    real, dimension(:), intent(out) :: L           ! path distance array
    real, dimension(:), intent(out) :: tauExt           ! total optical depth
    real, dimension(:), intent(out) :: tauAbs           ! thermal optical depth
    real, dimension(:), intent(out) :: tauSca           ! e-scattering optical depth
    integer, intent(out)      :: nTau                   ! size of optical depth arrays
    logical, intent(in)       :: thin_disc_on           ! T to include thin disc
    logical, intent(in)       :: opaqueCore             ! is the core opaque
    real, intent(out)         :: escProb                ! the escape probability
    real, intent(in)          :: lamStart, lamEnd
    integer, intent(in)       :: nLambda
    real, intent(inout)       :: tauCont(maxTau, nLambda)
    real, intent(inout)       :: linePhotonAlbedo(1:maxTau) ! the line photon albedo along the ray
    logical, intent(out)      :: hitcore                ! has the photon hit the core
    logical, intent(in)       :: thinLine               ! ignore line absorption of continuum
    logical, intent(in)       :: redRegion              ! use raman scattering red region opacities
    logical, intent(in)       :: usePops
    integer, intent(in)       :: mLevel, nLevel
    real, intent(in)          :: fStrength, gM, gN
    real, intent(in)          :: sampleFreq             ! max. samples per grid cell
    integer, intent(out)      :: error                  ! error code returned

    !    
    type(octalVector)         :: aVecOctal              ! octalVector version of 'aVec'
    type(octalVector)         :: uHatOctal              ! octalVector version of 'uHat'
    type(VECTOR)              :: rVec                   ! position vector
    logical                   :: contPhoton             ! is this a continuum photon?    
    real(double) :: nu, nu2, nu_p, nu0, nu0_p           ! frequencies   
    real :: thisVel
    integer :: ilambda                                ! wavelength index    
    real :: kappa_total
    real(double) :: chil                                       ! line opacity
    integer :: i, j                   ! counters
    
    
    real(double) :: deltaNu, DopplerWidth
    real :: meanMoleMass = mHydrogen
    real(double) :: a
    real(double) :: Hay, dv, Vn1, Vn2, Vrel
    real(double) :: T_mid, Ne_mid, N_HI_mid, chiline_mid, projVel_mid
    integer :: newNtau
    real :: T, lam, tmp
    real :: sqrt_pi
    real(double) :: dtau, taul, tau_tot
    real :: dtau_max
    !
    ! WORK ARRAYS
    logical, save :: first_time = .true.
    real(double), allocatable, save :: rho(:)         ! density (size=maxTau)
    real, allocatable, save :: temperature(:) ! temperature (size=maxTau)
    real(double), allocatable, save :: Ne(:)  ! electron density (size=maxTau)
    logical, allocatable, save :: inFlow(:)   ! inFlow flags (size=maxTau)   
    real(double), allocatable, save :: levelPop(:,:)  ! (size=maxTau x grid%maxLevels)
    type(vector), allocatable, save :: velocity(:) ! size=maxTau       ! 
    real, allocatable, save :: chiLine(:)          ! line optical depth (size=maxTau)
    real, allocatable, save  :: tauAbsLine(:)      ! (size=maxTau)
    real, allocatable, save  :: velocityDeriv(:)   ! directional derivative (size=maxTau)
    real, allocatable, save :: dL(:)     ! distance increment array (size=maxTau)
    real(double), allocatable, save :: projVel(:)   ! (size=maxTau)
    real(double), allocatable, save :: ksca(:)              ! (size=maxTau)
    real(double), allocatable, save :: kabs(:)              ! (size=maxTau)
    real, allocatable, save :: newL(:)              ! (size=maxTau)
    logical, allocatable, save :: newInFlow(:)      ! size=maxtau
    real, allocatable, save :: N_HI(:)              ! size = maxTau
    real(double), allocatable, save :: dummy(:)     ! work array (size=maxTau)

    !
    ! initialize variables
    sqrt_pi = SQRT(pi)
    hitcore = .false.
    rVec = aVec
    escProb = 1.
    error = 0
   

    if (usePops) then
       print *, "Error:: usePop=F option has not been implemented in integratePathVoigtProf."
       stop
    end if



    ! Allocates memory for work arrays for the first time
    ! This should be faster than using automatic arrays which allocates
    ! and deallocates memory every single time.  These array are SAVED
    ! so they should be available next time as well.
    if (first_time) then
       first_time = .false.
       ALLOCATE(rho(maxTau))         
       ALLOCATE(temperature(maxTau)) 
       ALLOCATE(Ne(maxTau))  
       ALLOCATE(inFlow(maxTau))
       ALLOCATE(levelPop(maxTau,grid%maxLevels))  
       ALLOCATE(velocity(maxTau)) 
       ALLOCATE(chiLine(maxTau)) 
       ALLOCATE(tauAbsLine(maxTau)) 
       ALLOCATE(velocityDeriv(maxTau)) 
       ALLOCATE(dL(maxTau))
       ALLOCATE(projVel(maxTau))
       ALLOCATE(kabs(maxTau))       
       ALLOCATE(ksca(maxTau))
       ALLOCATE(newL(maxTau))   
       ALLOCATE(newInFlow(maxTau))
       ALLOCATE(N_HI(maxTau))
       ALLOCATE(dummy(maxTau))
    end if



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


    ! find the projected velocity [c]
    Vn1 = uHat .dot. vVec         ! projected velocity of the local gas at emission location
    thisVel = Vn1 + (lambda0-wavelength)/lambda0
!    thisVel = (wavelength-lambda0)/lambda0


    escProb = 1.
    
    CALL startReturnSamples2 (aVecOctal,uHatOctal,grid,sampleFreq,nTau,       &
         maxTau,thin_disc_on,opaqueCore,hitcore,usePops,iLambda,error,&
         L,kappaAbs=kAbs,kappaSca=kSca, velocity=velocity,&
         velocityDeriv=velocityDeriv,chiLine=chiLine,    &
         levelPop=levelPop,rho=rho, &
         temperature=temperature, Ne=Ne, inFlow=inFlow, etaCont=dummy, etaLine=dummy)



    !---------------------------------------------------------------
    L(1) = 1.0e-25
    
    ! Some elements of L array could be dupulicated, so we removed them since
    ! they could potentially cause problems in interpolation routines used later.
    newntau = 0
    do i = 1, ntau-1       
       if (L(i) /= L(i+1)) then
          newNtau = newNtau + 1
          L(newNtau) = L(i)
          kabs(newNtau) = kabs(i)
          ksca(newNtau) = ksca(i)
          velocity(newNtau) = velocity(i)
          chiline(newNtau) = chiline(i)
          rho(newNtau) = rho(i)
          Ne(newNtau) = Ne(i)
          temperature(newNtau) = temperature(i)
          inFlow(newNtau) = inFlow(i)
       end if
    end do
    ! and .. the last point
    if (L(ntau-1) /= L(ntau) .and. ntau >=2) then
       newNtau = newNtau + 1
       L(newNtau) = L(ntau)
       kabs(newNtau) = kabs(ntau)
       ksca(newNtau) = ksca(ntau)
       velocity(newNtau) = velocity(ntau)
       chiline(newNtau) = chiline(ntau)
       rho(newNtau) = rho(ntau)
       Ne(newNtau) = Ne(ntau)
       temperature(newNtau) = temperature(ntau)
       inFlow(newNtau) = inFlow(ntau)
    end if
    !
    ntau= newNtau    

!    if (nTau <= 2) then
    if (nTau <= 3) then
       !       print *, 'The code does not yet include routines for handling ',&
       !            'the case where nTau <= 2 (in integratePathAMR)'
       error = -10
       return
    end if



    !----------------------------------------------------------------------------
    Vn2 = velocity(1) .dot. uHat  ! projected velocity of the local gas at this location
    
    ! projected velocities 
    forall (i = 1:nTau)
       ! (+ve when moving toward in the direction of photon.)
       projVel(i) = dble(velocity(i) .dot. uHat)  
    end forall
    
    dL(1:nTau-1) = L(2:nTau) - L(1:nTau-1)
    ! Number density of HI N_H = N_HI + Np, but Np=Ne
    N_HI(1:ntau) = rho(1:ntau)/meanMoleMass - Ne(1:nTau)


!---------------------------------------------------------------
!---------------------------------------------------------------
    !
    ! Finding the initial optical depth scales...
    ! and resample points along the ray.
    nu0 = cSpeed / (lambda0*angstromtocm)    ! line center frequency
    nu = cSpeed / (wavelength*angstromtocm)  ! freq of this photon
    nu_p = nu  ! freq in the rest frame of local gas
!    nu_p = nu/(1.0d0-Vn1)  ! the photon frequency in observer's frame
!    nu_p = nu*(1.0d0+Vn1)  ! binominal expansion of above


    if (.not.contPhoton) then
       tauAbs(1:2) = 1.0e-25
       tauSca(1:2) = 1.0e-25
       tauExt(1:2) = 1.0e-25
       tauAbsLine(1:2) = 1.0e-25
       linePhotonAlbedo(1:2) = 1.0e-25
       do i = 2, nTau
!       do i = 3, nTau
          if (inflow(i-1)) then
          T_mid = 0.5d0*(temperature(i-1)+temperature(i))
          Ne_mid = 0.5d0*(Ne(i-1)+Ne(i))
          N_HI_mid = 0.5d0*(N_HI(i-1)+N_HI(i))
          chiline_mid = 0.5d0*(chiline(i-1)+chiline(i))
          projVel_mid = 0.5d0*(projVel(i-1)+projVel(i))
          !----------------------------------------------
!          T_mid = temperature(i-1)
!          Ne_mid = Ne(i-1)
!          N_HI_mid = N_HI(i-1)
!          chiline_mid = chiline(i-1)
!          projVel_mid = projVel(i-1)
          !----------------------------------------------
          T_mid = MAX(T_mid, 10000.0d0) ! [K]  To avoid a tiny Doppler width
          
          ! relative velocity wrt the location of photon (CMF)
          Vrel = projVel_mid -  Vn1

!          ! relative velocity wrt the observer
!          Vrel = projVel_mid 
          
          ! The line centre of absorption profile shifted by Doppler.
          nu0_p = nu0/(1.0d0-Vrel)  ! [Hz] 
!          nu0_p = nu0*(1.0d0+Vrel)  ! [Hz]  binomial expansion
                    
          DopplerWidth = nu0_p/cSpeed * sqrt(2.*kErg*T_mid/meanMoleMass) !eq 7  [Hz]
          
          a = bigGamma(N_HI_mid, T_mid, Ne_mid, nu0_p) / (fourPi * DopplerWidth) ! [-]
          deltaNu = nu_p - nu0_p     !  [Hz]
!          deltaNu = nu - nu0    !  [Hz]
          dv = deltaNu/DopplerWidth  ! [-]
          Hay = voigtn(a,dv)
          chil = chiLine_mid / (sqrt_pi*DopplerWidth) * Hay ! equation 5
!          ! transform this back to observer's frame value
          chil = chil/(1.0d0+projVel_mid)
          dtau = chil*dL(i-1)
          ! line albedo added here - tjh
          kappa_total = (kSca(i-1) + kAbs(i-1) + chil)
          !             kappa_total = (kSca(i) + kAbs(i))  ! don't include line opacity for albedo
          if (kappa_total >0.0) then
             linePhotonAlbedo(i) = kSca(i-1) / kappa_total
          else
             linePhotonAlbedo(i) = 0.0
          end if

          tauAbsLine(i) = tauAbsLine(i-1) +  abs(dtau)
          tauSca(i) = tauSca(i-1) + dL(i-1)*(0.5*(ksca(i)+ksca(i-1)))
          tauAbs(i) = tauAbs(i-1) + dL(i-1)*(0.5*(kabs(i)+kabs(i-1)))
          !------------------------------------------------------------
!          tauSca(i) = tauSca(i-1) + dL(i-1)*ksca(i-1)
!          tauAbs(i) = tauAbs(i-1) + dL(i-1)*kabs(i-1)

          tauExt(i) = tauSca(i) + tauAbs(i) + tauAbsLine(i)

          else
             tauAbsLine(i) = tauAbsLine(i-1) 
             tauSca(i) = tauSca(i-1) 
             tauAbs(i) = tauAbs(i-1) 
             tauExt(i) = tauExt(i-1)
             linePhotonAlbedo(i) = linePhotonAlbedo(i-1)
          end if ! (inflow(i-1)
       enddo
       taul =  tauAbsLine(nTau)
       tau_tot =  tauExt(nTau)
       escProb = 1.0
    else  ! continuum photon
       escProb = 1.0
       tauAbs(1:2) = 1.0e-25
       tauSca(1:2) = 1.0e-25
       tauExt(1:2) = 1.0e-25
       do i = 2, nTau
          if (inflow(i-1)) then
             tauSca(i) = tauSca(i-1) + dL(i-1)*(0.5*(ksca(i)+ksca(i-1)))
             tauAbs(i) = tauAbs(i-1) + dL(i-1)*(0.5*(kabs(i)+kabs(i-1)))
             !------------------------------------------------------------
!             tauSca(i) = tauSca(i-1) + dL(i-1)*ksca(i-1)
!             tauAbs(i) = tauAbs(i-1) + dL(i-1)*kabs(i-1)

             tauExt(i) = tauSca(i) + tauAbs(i)
          else
             tauAbsLine(i) = tauAbsLine(i-1) 
             tauSca(i) = tauSca(i-1) 
             tauAbs(i) = tauAbs(i-1) 
             tauExt(i) = tauExt(i-1)
          end if ! (inflow(i-1) 
       end do
       tauAbsLine(1:ntau)=1.0e-25
    endif 



    ! Now resample rays using tauExt values
!    dtau_max = 10.0
    dtau_max = 0.05
    call resampleRay_tau(L, nTau, tauExt, dtau_max, maxTau, newL, newNTau, &
         inflow, newInFlow)
!    call resampleRay_tau(L, nTau, tauAbsLine, dtau_max, maxTau, newL, newNTau, &
!         inflow, newInFlow)


    ! Now interpolate on to newly sampled ray    
!    call quadraticResample_dble(L, projVel, nTAu, maxTau, newL, newNtau)
!    call quadraticResample(L, kSca, nTAu, maxTau, newL, newNtau)
!    call quadraticResample(L, kAbs, nTAu, maxTau, newL, newNtau)
!    call quadraticResample(L, chiline, nTAu, maxTau, newL, newNtau)
!    call quadraticResample(L, temperature, nTAu, maxTau, newL, newNtau)
!    call quadraticResample(L, N_HI, nTAu, maxTau, newL, newNtau)
!    call quadraticResample_dble(L, Ne, nTAu, maxTau, newL, newNtau)

    call linearResample_dble(L, projVel, nTAu, maxTau, newL, newNtau)
!    call linearResample(L, kSca, nTAu, maxTau, newL, newNtau)
!    call linearResample(L, kAbs, nTAu, maxTau, newL, newNtau)
    call linearResample_dble(L, kSca, nTAu, maxTau, newL, newNtau)
    call linearResample_dble(L, kAbs, nTAu, maxTau, newL, newNtau)
    call linearResample(L, chiline, nTAu, maxTau, newL, newNtau)
    call linearResample(L, temperature, nTAu, maxTau, newL, newNtau)
    call linearResample(L, N_HI, nTAu, maxTau, newL, newNtau)
    call linearResample_dble(L, Ne, nTAu, maxTau, newL, newNtau)


    ! updating the ray
    nTau = newNtau
    L(1:nTau) = newL(1:nTau)
    L(1) = 1.0e-25
    dL(1:nTau-1) = L(2:nTau) - L(1:nTau-1)
    ! updating the inFlow flags.
    inFlow(1:nTau) = newInFlow(1:nTau)  

    Ne(1:ntau) = ABS(Ne(1:ntau))  ! just for safty

!     !---------------------------------------------------------------
!     ! additional resampling for the place where the velocity
!     ! is changing slowly.
!     call resampleRay3(L, nTau, projVel, maxTau, newL, newNTau, inflow, newInFlow)

!     ! Now interpolate on to newly sampled ray    
!     call linearResample_dble(L, projVel, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, kSca, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, kAbs, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, chiline, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, temperature, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, N_HI, nTAu, maxTau, newL, newNtau)
!     call linearResample_dble(L, Ne, nTAu, maxTau, newL, newNtau)


!     ! updating the ray
!     nTau = newNtau
!     L(1:nTau) = newL(1:nTau)
!     L(1) = 1.0e-25
!     dL(1:nTau-1) = L(2:nTau) - L(1:nTau-1)
!     ! updating the inFlow flags.
!     inFlow(1:nTau) = newInFlow(1:nTau)  

!     Ne(1:ntau) = ABS(Ne(1:ntau))  ! just for safty


!---------------------------------------------------------------
!---------------------------------------------------------------
  
!     ! Now resample rays using velocity gradiants
! !    call resampleRay(L, nTau, projVel, maxTau, newL, newNTau, inflow, newInFlow)

!     thisVel = (lambda0 - wavelength)/lambda0  
!     thisVel = thisVel  + (uHat .dot. vVec)    ! (+ve when toward observer!)
!     call resampleRay2(L, nTau, projVel, thisVel, maxTau, newL, newNTau, inflow, newInFlow)


!     ! Now interpolate on to newly sampled ray    
!     call linearResample_dble(L, projVel, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, kSca, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, kAbs, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, chiline, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, temperature, nTAu, maxTau, newL, newNtau)
!     call linearResample(L, N_HI, nTAu, maxTau, newL, newNtau)
!     call linearResample_dble(L, Ne, nTAu, maxTau, newL, newNtau)


!     ! updating the ray
!     nTau = newNtau
!     L(1:nTau) = newL(1:nTau)
!     L(1) = 1.0e-25
!     dL(1:nTau-1) = L(2:nTau) - L(1:nTau-1)
!    ! updating the inFlow flags.
!    inFlow(1:nTau) = newInFlow(1:nTau)  

!     Ne(1:ntau) = ABS(Ne(1:ntau))  ! just for safty

!---------------------------------------------------------------
!---------------------------------------------------------------




    !
    !
    ! Now we compute optical depth again with improved samping
    !
    !    


    ! line photons
    nu0 = cSpeed / (lambda0*angstromtocm)    ! line center frequency
    nu = cSpeed / (wavelength*angstromtocm)  ! freq of this photon
    nu_p = nu  ! freq in the rest frame of local gas
!    nu_p = nu/(1.0d0-Vn1)
!    nu_p = nu*(1.0d0+Vn1)  ! binominal expansion of above    
 
    if (.not.contPhoton) then
       tauAbs(1:2) = 1.0e-25
       tauSca(1:2) = 1.0e-25
       tauExt(1:2) = 1.0e-25
       tauAbsLine(1:2) = 1.0e-25
       linePhotonAlbedo(1:2) = 1.0e-25
       do i = 2, nTau
          if (inflow(i-1)) then
             ! Evaluating the values in the mid point
             T_mid = 0.5d0*(temperature(i-1)+temperature(i))
             Ne_mid = 0.5d0*(Ne(i-1)+Ne(i))
             N_HI_mid = 0.5d0*(N_HI(i-1)+N_HI(i))
             chiline_mid = 0.5d0*(chiline(i-1)+chiline(i))
             projVel_mid = 0.5d0*(projVel(i-1)+projVel(i))
             !----------------------------------------------
!             T_mid = temperature(i-1)
!             Ne_mid = Ne(i-1)
!             N_HI_mid = N_HI(i-1)
!             chiline_mid = chiline(i-1)
!             projVel_mid = projVel(i-1)
             !----------------------------------------------

             T_mid = MAX(T_mid, 1000.0d0) ! [K]  To avoid a tiny Doppler width
             
             ! relative velocity wrt the location of photon (CMF)
             Vrel = projVel_mid -  Vn1

!             ! relative velocity wrt the observer
!             Vrel = projVel_mid
          
             ! The line centre of absorption profile shifted by Doppler.
             nu0_p = nu0/(1.0d0-Vrel)  ! [Hz] 
                    
             DopplerWidth = nu0_p/cSpeed * sqrt(2.*kErg*T_mid/meanMoleMass) !eq 7  [Hz]
          
             a = bigGamma(N_HI_mid, T_mid, Ne_mid, nu0_p) / (fourPi * DopplerWidth) ! [-]
             deltaNu = nu_p - nu0_p     !  [Hz]
!             deltaNu = nu - nu0    !  [Hz]
             dv = deltaNu/DopplerWidth  ! [-]
             Hay = voigtn(a,dv)
             chil = chiLine_mid / (sqrt_pi*DopplerWidth) * Hay ! equation 5
!             ! transform this back to observer's frame value
             chil = chil/(1.0d0+projVel_mid)
             dtau = chil*dL(i-1)
             ! line albedo added here - tjh
             kappa_total = (kSca(i-1) + kAbs(i-1) + chil)
             !             kappa_total = (kSca(i) + kAbs(i))  
             if (kappa_total >0.0) then
                linePhotonAlbedo(i) = kSca(i-1) / kappa_total
             else
                linePhotonAlbedo(i) = 0.0
             end if
          
             tauAbsLine(i) = tauAbsLine(i-1) +  abs(dtau)
             tauSca(i) = tauSca(i-1) + dL(i-1)*(0.5*(ksca(i)+ksca(i-1)))
             tauAbs(i) = tauAbs(i-1) + dL(i-1)*(0.5*(kabs(i)+kabs(i-1)))
             !------------------------------------------------------------
!             tauSca(i) = tauSca(i-1) + dL(i-1)*ksca(i-1)
!             tauAbs(i) = tauAbs(i-1) + dL(i-1)*kabs(i-1)

             tauExt(i) = tauSca(i) + tauAbs(i) + tauAbsLine(i)

          else
             tauAbsLine(i) = tauAbsLine(i-1) 
             tauSca(i) = tauSca(i-1) 
             tauAbs(i) = tauAbs(i-1) 
             tauExt(i) = tauExt(i-1)
             linePhotonAlbedo(i) = linePhotonAlbedo(i-1)
          end if ! inflow(i-1)) then

       enddo
       taul =  tauAbsLine(nTau)
       tau_tot =  tauExt(nTau)
       escProb = 1.0
    else  ! continuum photon
       escProb = 1.0
       tauAbs(1:2) = 1.0e-25
       tauSca(1:2) = 1.0e-25
       tauExt(1:2) = 1.0e-25
       do i = 2, nTau
          if ( inflow(i-1) ) then
             tauSca(i) = tauSca(i-1) + dL(i-1)*(0.5*(ksca(i)+ksca(i-1)))
             tauAbs(i) = tauAbs(i-1) + dL(i-1)*(0.5*(kabs(i)+kabs(i-1)))
             !------------------------------------------------------------
!             tauSca(i) = tauSca(i-1) + dL(i-1)*ksca(i-1)
!             tauAbs(i) = tauAbs(i-1) + dL(i-1)*kabs(i-1)

             tauExt(i) = tauSca(i) + tauAbs(i)
          else
             tauAbsLine(i) = tauAbsLine(i-1) 
             tauSca(i) = tauSca(i-1) 
             tauAbs(i) = tauAbs(i-1) 
             tauExt(i) = tauExt(i-1)
          end if ! inflow(i-1)
       end do
    endif 

    

    ! continuum photons
    ! If the line is optically thick, we consider the absorption of 
    ! of the contiuum photons by line....
    nu0 = cSpeed / (lambda0*angstromtocm)    ! line center frequency
    if (contPhoton .and. .not.thinLine) then
       tauCont(1:2,1:nlambda)=0.d0
       do j = 1, nLambda
          ! compute projected velocity of this bin  
!          thisVel = (lambda0-grid%lamArray(j))/lambda0  + (uHat .dot. vVec)  &
!               - (wavelength - lambda0)/lambda0 ! (+ve if moving toward)
!          thisVel = (lambda0-grid%lamArray(j))/lambda0   &
!               - (wavelength - lambda0)/lambda0 ! (+ve if moving toward)

          lam = (wavelength-lambda0) + grid%lamArray(j) 
          nu = cSpeed / (lam*angstromtocm)   ! freq of this photon
          nu_p = nu
!          nu_p = nu/(1.0d0-Vn1)             ! freq in the rest frame of local gas
!          nu_p = nu/(1.0d0-(Vn1-Vn2))       ! the photon frequency in observer's frame
!          nu_p = nu0/(1.0d0-thisVel)         ! freq in the rest frame of local gas
!          nu_p = nu0*(1.0d0+thisVel)         ! binominal expansion of above.

          do i = 2, nTau
!          do i = 3, nTau
             if (inflow(i-1)) then
                ! Evaluating the values in the mid point
                T_mid = 0.5d0*(temperature(i-1)+temperature(i))
                Ne_mid = 0.5d0*(Ne(i-1)+Ne(i))
                N_HI_mid = 0.5d0*(N_HI(i-1)+N_HI(i))
                chiline_mid = 0.5d0*(chiline(i-1)+chiline(i))
                projVel_mid = 0.5d0*(projVel(i-1)+projVel(i))
                !----------------------------------------------
!                T_mid = temperature(i-1)
!                Ne_mid = Ne(i-1)
!                N_HI_mid = N_HI(i-1)
!                chiline_mid = chiline(i-1)
!                projVel_mid = projVel(i-1)
             
                T_mid = MAX(T_mid, 1000.0d0) ! [K]  To avoid a tiny Doppler width                
                ! relative velocity wrt the location of photon (CMF)
                Vrel = projVel_mid -  Vn1

!                ! relative velocity wrt the observer
!                Vrel = projVel_mid
             
                ! The line centre of absorption profile shifted by Doppler.
                nu0_p = nu0/(1.0d0-Vrel)  

             
                DopplerWidth = nu0_p/cSpeed * sqrt(2.*kErg*T_mid/meanMoleMass) !eq 7  [Hz]
             
                a = bigGamma(N_HI_mid, T_mid, Ne_mid, nu0_p) / (fourPi * DopplerWidth) ! [-]
                deltaNu = nu_p - nu0_p     !  [Hz]
!                deltaNu = nu - nu0     !  [Hz]
                dv = deltaNu/DopplerWidth  ! [-]
                Hay = voigtn(a,dv)
                chil = chiLine_mid / (sqrt_pi*DopplerWidth) * Hay ! equation 5
                ! transform this back to observer's frame value
                chil = chil/(1.0d0+projVel_mid)
                dtau = chil*dL(i-1)
                tauCont(i,j) = tauCont(i-1,j) + abs(dtau)
             else
                tauCont(i,j) = tauCont(i-1,j)
             endif ! (inflow(i-1))
          enddo
       enddo   ! over wavelength array
    else  ! for optically thin line option
       tauCont(1:nTau,1:nLambda) = 0.
    end if



  end subroutine integratePathVoigtProf
  


  !
  !
  !
  subroutine intersectCubeAMR(grid, posVec, direction, tval)
   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTALVECTOR), intent(in) :: posVec
   type(OCTALVECTOR), intent(in) :: direction
   real(oct), intent(out) :: tval
   !
   type(OCTALVECTOR) :: norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(OCTALVECTOR) :: subcen, point
   integer :: subcell
   
   real(oct) :: t(6),denom(6)
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
     stop
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
subroutine test_optical_depth(gridUsesAMR, VoigtProf, &
     amrGridCentre, sphericityTest,  &
     dir_obs, wavelength,  lambda0, grid, thin_disc_on, opaqueCore, lamStart, lamEnd,   &
     ThinLine, lineResAbs, nUpper, nLower, sampleFreq, useinterp, Rstar, coolStarPosition, maxTau, nSource, source)

  implicit none
  
  integer :: nsource
  type(SOURCETYPE) :: source(:)
  logical, intent(in)           :: gridUsesAMR     !
  logical, intent(in)           :: VoigtProf       ! T to use Voigt Profile
  type(OCTALVECTOR), intent(in) :: amrGridCentre  ! central coordinates of grid
  logical, intent(in)           :: sphericityTest !

  type(vector), intent(in)  :: dir_obs            ! direction to the observer (unit vector)
  real, intent(in)          :: wavelength         ! the wavelength 
  real, intent(in)          :: lambda0            ! rest wavelength of line
  type(GRIDTYPE), intent(in):: grid               ! the opacity grid
  !
  logical, intent(in)       :: thin_disc_on      ! T to include thin disc
  logical, intent(in)       :: opaqueCore         ! is the core opaque
  real, intent(in)          :: lamStart, lamEnd
  !
  logical, intent(in)       :: thinLine           ! ignore line absorption of continuum
  logical, intent(in)       :: lineResAbs    ! T if you want to include absorption
  !                                               !   of line by line resonance zones.
  integer, intent(in)       :: nUpper, nLower
  real, intent(in)          :: sampleFreq         ! max. samples per grid cell
  !
  logical, intent(in)       :: useinterp
  real, intent(in)          :: Rstar
  type(VECTOR), intent(in)  :: coolStarPosition
  integer, intent(in) :: maxTau
  
  integer :: nLambda
  ! local variables
  real, allocatable :: lambda(:)    ! path distance array  (SIZE=maxTau)
  real, allocatable :: tauExt(:)    ! optical depth  (SIZE=maxTau)
  real, allocatable :: tauAbs(:)    ! optical depth  (SIZE=maxTau)
  real, allocatable :: tauSca(:)    ! optical depth  (SIZE=maxTau)
  real, allocatable :: tauCont(:,:) !tauCont(maxTau,nLambda)
  real, allocatable  :: linePhotonAlbedo(:) ! the line photon albedo along the ray


  integer  :: nTau        ! size of optical depth arrays
  real     :: escProb     ! the escape probability
  logical  :: hitcore     ! has the photon hit the core
  real     :: localTau
  integer  :: error       ! error code returned
  integer  :: i 
  type(vector) :: tempVec
  type(octalvector) :: zerovec
  real :: junk
  type(OCTALVECTOR) :: octVec, position
  logical :: contPhoton = .true.

  integer :: ntest = 301
  real :: theta
  real(oct) :: x1, x2
  integer, parameter :: UNOUT1 = 23
  integer, parameter :: UNOUT2 = 24
  
  real(oct) :: R


  allocate(lambda(maxTau), tauExt(maxTau), tauAbs(maxTau), tauSca(maxTau), linePhotonAlbedo(maxtau))
  
  nlambda = grid%nlambda
  allocate(tauCont(maxTau,nLambda))

  ! chuck out some useful information to the user
    
  write(*,'(a,f7.1,a)') "Cross-sections at ",wavelength, " angstroms"
  write(*,'(a)') "------------------------------------------"
  write(*,*) " "

  zeroVec = OCTALVECTOR(0.,0.,0.)
  
  !
  ! The distance from the center from which a test ray starts.
  if (grid%geometry == "cluster") then
     R = 1.0d-3
  else
!     R = 1.001*rStar
     R = 0.01d0*rStar
  end if
  if (gridUsesAMR) then

  
  !
  ! test along x-axis 
  !
  octVec = OCTALVECTOR(1.d0, 0.d0, 0.d0)
!  position = (octVec*R)
  R = max(R, 1.e-3*rSol/1.e10)
  position = (octVec*R) + grid%starPos1
  position%x = max(position%x, grid%halfSmallestSubcell)
  position%y = max(position%y, grid%halfSmallestSubcell)
  position%z = max(position%z, grid%halfSmallestSubcell)
  call integratePath(gridUsesAMR, VoigtProf, &
       wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  position, &
       octVec, grid, lambda, tauExt, tauAbs, &
       tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, contPhoton , &
       lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
       .false., nUpper, nLower, 0., 0., 0., junk,&
       sampleFreq, error, useInterp, rStar, coolStarPosition, nSource, source)
  if (error < 0) then
     write(*,*) '   Error encountered in cross-sections!!! (error = ',error,')'
  end if
  
  if (nTau > 2) then 
     write(*,'(a,1pe10.3)') "Optical depth in x-axis from centre: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth in x-axis from centre: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth in x-axis from centre: ",tauSca(ntau)
     write(*,'(a,i4)') "Number of samples: ",nTau
     write(*,*) " "
  end if

  !
  ! test along y-axis 
  !
  octVec = OCTALVECTOR(0.d0, 1.d0, 0.d0)
!  position = (octVec*R) 
  R = max(R, 1.e-3*rSol/1.e10)
  position = (octVec*R) + grid%starPos1
  position%x = max(position%x, grid%halfSmallestSubcell)
  position%y = max(position%y, grid%halfSmallestSubcell)
  position%z = max(position%z, grid%halfSmallestSubcell)
  call integratePath(gridUsesAMR, VoigtProf, &
       wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  position, &
       octVec, grid, lambda, tauExt, tauAbs, &
       tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, contPhoton , &
       lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
       .false., nUpper, nLower, 0., 0., 0., junk,&
       sampleFreq,error, useInterp, rStar, coolStarPosition, nSource, source)
  if (error < 0) then
     write(*,*) '   Error encountered in cross-sections!!! (error = ',error,')'
  end if
  
  if (nTau > 2) then 
     write(*,'(a,1pe10.3)') "Optical depth in y-axis from centre: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth in y-axis from centre: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth in y-axis from centre: ",tauSca(ntau)
     write(*,'(a,i4)') "Number of samples: ",nTau
     write(*,*) " "
  end if

  !
  ! test along z-axis 
  !
  octVec = OCTALVECTOR(1.d-3, 1.2d-4, 1.d0)
!  position = (octVec*R)
  R = max(R, 1.e-3*rSol/1.e10)
  position = (octVec*R) + grid%starPos1
  position%x = max(position%x, grid%halfSmallestSubcell)
  position%y = max(position%y, grid%halfSmallestSubcell)
  position%z = max(position%z, grid%halfSmallestSubcell)
  call integratePath(gridUsesAMR, VoigtProf, &
       wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  position, &
       octVec, grid, lambda, tauExt, tauAbs, &
       tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, contPhoton , &
       lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
       .false., nUpper, nLower, 0., 0., 0., junk,&
       sampleFreq,error, useInterp, rStar, coolStarPosition, nSource, source)
  if (error < 0) then
     write(*,*) '   Error encountered in cross-sections!!! (error = ',error,')'
  end if
  
  if (nTau > 2) then 
     write(*,'(a,1pe10.3)') "Optical depth in z-axis from centre: ",tauExt(ntau)
     write(*,'(a,1pe10.3)') "Absorption depth in z-axis from centre: ",tauAbs(ntau)
     write(*,'(a,1pe10.3)') "Scattering depth in z-axis from centre: ",tauSca(ntau)
     write(*,'(a,i4)') "Number of samples: ",nTau
     write(*,*) " "
  end if

  !
  ! test towards the observer
  !
  octVec = s2o(dir_obs)
!  position = (octVec*R) + amrGridCentre
  R = max(R, 1.e-3*rSol/1.e10)
  position = (octVec*R) + grid%starPos1
  position%x = max(position%x, grid%halfSmallestSubcell)
  position%y = max(position%y, grid%halfSmallestSubcell)
  position%z = max(position%z, grid%halfSmallestSubcell)
  call integratePath(gridUsesAMR, VoigtProf, &
       wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5), &
       position, s2o(dir_obs), grid, lambda, &
       tauExt, tauAbs, tauSca, linePhotonalbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, &
       contPhoton, lamStart, lamEnd, nLambda, tauCont, &
       hitCore, thinLine,lineResAbs, .false.,  &
       .false., nUpper, nLower, 0., 0., 0., junk,sampleFreq,error, &
       useinterp, rStar, coolStarPosition, nSource, source)

  if (error < 0) then
     write(*,*) '   Error encountered in test towards observer!!! (error = ',error,')'
  end if     
     
  write(*,'(a,1pe10.3)') "Optical depth to observer: ",tauExt(ntau)
  write(*,'(a,1pe10.3)') "Absorption depth to observer: ",tauAbs(ntau)
  write(*,'(a,1pe10.3)') "Scattering depth to observer: ",tauSca(ntau)
  write(*,'(a,i4)') "Number of samples: ",nTau
  write(*,*) " "


  !
  !
  !
  if (sphericityTest) then
     write(*,'(a)') "Sphericity test - 100 random directions:"
     do i = 1, 100
        tempVec = randomUnitVector()
        call IntegratePath(gridUsesAMR, VoigtProf, &
             lambda0,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5), &
             amrGridCentre, s2o(tempVec), grid, lambda, &
             tauExt, tauAbs, tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, &
             contPhoton, lamStart, lamEnd, nLambda, tauCont, &
             hitCore, thinLine, lineResAbs, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk, sampleFreq,error, &
             useinterp, rStar, coolStarPosition, nSource, source)
        if (error < 0) then
           write(*,*) 'Error handling not implemented for ''sphericityTest''!!! (error = ',error,')'
        end if
        write(*,'(a,1pe10.3,1pe10.3,1pe10.3)') "Optical depths: ",tauExt(nTau),tauAbs(nTau),tauSca(nTau)
     enddo
  endif


  !
  ! Compute the optical depth on x-y and z-x planes
  !  
!  goto 999
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
!        position = (octVec*R) + amrGridCentre
        position = (octVec*R) + grid%starPos1
        position%x = max(position%x, grid%halfSmallestSubcell)
        position%y = max(position%y, grid%halfSmallestSubcell)
        position%z = max(position%z, grid%halfSmallestSubcell)
        ! continuum
        call IntegratePath(gridUsesAMR, VoigtProf, &
             wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  &
             position, octVec, grid, lambda, tauExt, tauAbs, &
             tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .true. , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,Error,&
             useinterp, rStar, coolStarPosition, nSource, source)

        write(UNOUT1, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau
        
        ! line 
        call IntegratePath(gridUsesAMR, VoigtProf, &
             wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  &
             position, octVec, grid, lambda, tauExt, tauAbs, &
             tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .false. , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,Error, &
             useinterp, rStar, coolStarPosition, nSource, source)

        write(UNOUT2, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau
        
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
!        position = (octVec*R)  + amrGridCentre
        position = (octVec*R) + grid%starPos1
        position%x = max(position%x, grid%halfSmallestSubcell)
        position%y = max(position%y, grid%halfSmallestSubcell)
        position%z = max(position%z, grid%halfSmallestSubcell)
        ! continuum
        call IntegratePath(gridUsesAMR, VoigtProf, &
             wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  &
             position, octVec, grid, lambda, tauExt, tauAbs, &
             tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .true. , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,Error,&
             useinterp, rStar, coolStarPosition, nSource, source)
        
        write(UNOUT1, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau

        ! line 
        call IntegratePath(gridUsesAMR, VoigtProf, &
             wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  &
             position, octVec, grid, lambda, tauExt, tauAbs, &
             tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .false. , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,Error, &
             useinterp, rStar, coolStarPosition, nSource, source)
        
        write(UNOUT2, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau
        
     end do


     close(UNOUT1)
     close(UNOUT2)

!=== FOR DEBUG ONLY REMOVE THIS SECTION LATER ===========================================
     ! 
     ! Z-X plane 
     !
     open(unit=UNOUT1, file = 'tau_zx_plane_test.dat', status ='replace')
     open(unit=UNOUT2, file = 'taul_zx_plane_test.dat', status ='replace')
        
     write(UNOUT1, '(a)') '#    phi [deg]  -- tau(total) -- tau(abs) -- tau(scat) -- ntau'        
     write(UNOUT2, '(a)') '#    phi [deg]  -- taul(total) -- taul(abs) -- taul(scat) -- ntau'
     
     ! direction of the beam...
     do i = 1, ntest
        theta = 2.0d0*Pi*real(i-1)/real(ntest-1) + 0.01
        x1 = cos(theta); x2 = sin(theta)
        octVec = OCTALVECTOR(x2, 0.0, x1)
        !          call Normalize(octVec)  ! just in case ..
        
        ! position of emission
        position =  OCTALVECTOR(9.5d0, 0.0d0, 11.0d0) ! in a middle of stream
!        position =  OCTALVECTOR(34.0d0, 0.0d0, 6.0d0) ! in a middle of stream
!        position = (octVec*R) + grid%starPos1
        position%x = max(position%x, grid%halfSmallestSubcell)
        position%y = max(position%y, grid%halfSmallestSubcell)
        position%z = max(position%z, grid%halfSmallestSubcell)

        ! continuum
        call IntegratePath(gridUsesAMR, VoigtProf, &
             wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  &
             position, octVec, grid, lambda, tauExt, tauAbs, &
             tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .true. , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,Error,&
             useinterp, rStar, coolStarPosition, nSource, source)
        
        write(UNOUT1, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau

        ! line 
        call IntegratePath(gridUsesAMR, VoigtProf, &
             wavelength,  lambda0, OCTALVECTOR(1.0e-5,1.0e-5,1.0e-5),  &
             position, octVec, grid, lambda, tauExt, tauAbs, &
             tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .false. , &
             lamStart, lamEnd, nLambda, tauCont, hitCore, thinLine, lineResAbs, .false., &
             .false., nUpper, nLower, 0., 0., 0., junk,&
             sampleFreq,Error, &
             useinterp, rStar, coolStarPosition, nSource, source)
        
        write(UNOUT2, *)    theta*180.0/Pi,  tauExt(ntau), tauAbs(ntau), tauSca(ntau), ntau
        
     end do


     close(UNOUT1)
     close(UNOUT2)

!=== FOR DEBUG ONLY REMOVE THIS SECTION LATER ===========================================
       
  end if ! gridUsesAMR
999 continue

     
end subroutine test_optical_depth





end module path_integral
