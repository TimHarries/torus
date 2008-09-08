!
! this module defines and describes the photon type
! and includes subroutines for photon initialization
! and scattering
!

! written by tjh

! v1.0 on 16/08/99
! adaptive mesh stuff added 2002/08/06. nhs

module photon_mod

  use photoion_mod ! test
  use vector_mod            ! the vector maths
  use constants_mod         ! physical constants
  use phasematrix_mod       ! phase matrix manipulation
  use gridtype_mod          ! type definition for the 3-d grid
  use grid_mod              ! opacity grid routines
  use math_mod              ! maths utils
  use amr_mod               ! adaptive grid routines
  use utils_mod
  use phasematrix_mod
  use disc_class
  use source_mod
  use filter_set_class
  use surface_mod
  use kind_mod

  implicit none

  public

  ! the description of the photon type


  type PHOTON
     type(STOKESVECTOR) :: stokes         ! the stokes intensities
     real :: weight                       ! photon weight
     real :: lambda                       ! wavelength
     type(VECTOR) :: normal               ! scattering normal
     type(VECTOR) :: oldNormal            ! the last scattering normal
     type(OCTALVECTOR) :: position             ! the photon position
     type(OCTALVECTOR) :: direction            ! the photon direction
     type(VECTOR) :: velocity             ! the photon 'velocity'
     type(VECTOR) :: originalNormal
     logical :: contPhoton                ! continuum photon? or
     logical :: linePhoton                ! line photon?
     logical :: redPhoton                 ! raman photon in red region
     real :: redLambda                    ! raman wavelength in red
     logical :: fromStar1
     logical :: fromStar2
     logical :: resonanceLine             ! resonance line photon?
     logical :: thermal ! thermal emission?
     logical :: scattered ! has photon been scattered?
     logical :: stellar ! stellar emission
  end type PHOTON



contains

  ! this subroutine rotates the stokes vectors of a photon - this
  ! subroutine should probably moved from here to just work
  ! on a stokes vector (see phasematrix_mod)

  subroutine rotate(thisPhoton, cosx, sinx)

    type(PHOTON), intent(inout) :: thisPhoton          ! the photon
    real, intent(in) :: cosx, sinx       ! x should be twice theta
    real :: qdash, udash                 ! the rotated stokes vector

    ! perform the rotation

! COMMENTED OUT by RK for DEBUGGING PURPOSE
!    if ( (sinx**2 + cosx**2) < 0.999) then
!       write(*,*) "! non-normalized rotation attempt",sinx**2+cosx**2
!       stop
!    endif

    qdash = cosx * thisPhoton%stokes%q + sinx * thisPhoton%stokes%u
    udash =-sinx * thisPhoton%stokes%q + cosx * thisPhoton%stokes%u

    thisPhoton%stokes%q = qdash
    thisPhoton%stokes%u = udash

  end subroutine rotate


  ! this subroutine performs a scattering

  subroutine scatterPhoton(grid, thisPhoton, givenVec, outPhoton, mie, &
        miePhase, nDustType, nLambda, lamArray, nMuMie, ttau_disc_on, alpha_disc_param)

    
    type(GRIDTYPE) :: grid                       ! the opacity grid
    type(PHOTON) :: thisPhoton, outPhoton        ! current/output photon
    type(VECTOR) :: incoming, outgoing           ! directions
    type(VECTOR) :: obsNormal, refNormal         ! scattering normals
    integer :: nDustType
    type(VECTOR) :: sVec, givenVec      ! vectors
    type(VECTOR) :: zAxis                        ! the z-axis
    type(PHASEMATRIX) :: rayleighPhase           ! rayleigh phase matrix
    logical, intent(in) :: mie                   ! is this a mie scattering?
    integer :: i, j, k                              ! counters
    integer :: nLambda                           ! size of wavelength array
    real :: lamArray(:)
    integer :: nMumie                            ! number of mu angles for mie
    type(PHASEMATRIX), intent(in) :: miePhase(nDustType, nLambda, nMumie) ! mie phase matrices   
    type(PHASEMATRIX),save,allocatable :: miePhaseTemp(:)
    ! if the system has accretion disc around the obeject
    logical, intent(in) :: ttau_disc_on          
    ! to find if scattering occurs in the accretion disc
    type(alpha_disc), intent(in)  :: alpha_disc_param
    type(octal), pointer :: thisOctal
    integer :: subcell
    real :: costheta                             ! cos scattering angle
    real :: ang                                  ! scattering angle
    real :: r1, r2                               ! radii
    integer :: i1, i2, i3                        ! position indices
    real(oct) :: t1,t2,t3             ! interp factors
    real :: u,v,w,t                              ! direction components
    real :: sinGamma, cosGamma, sin2Gamma, cos2Gamma
    logical :: randomDirection                   ! is this a random direction
    real :: vRay, vOverCsqr
    real :: fac
    real :: nuRest
    type(octalVector) :: pointOctalVec
    type(octal), pointer :: octalLocation
    integer :: subcellLocation
    logical :: mie_scattering
    real :: weight
    logical, save :: firstTime = .true.
    real, allocatable, save :: cosArray(:)

    if (firstTime)  then
       allocate(cosArray(1:nMuMie))
       do i = 1, nMumie
          cosArray(i) = -1.d0 + 2.d0*dble(i-1)/dble(nMuMie)
       enddo
       allocate(miePhaseTemp(1:nDustType))
       firstTime = .false.
    endif
!    real :: dx

!    dx = abs(grid%xAxis(2) - grid%xAxis(1))

    ! initialize variables

    outPhoton = thisPhoton 
    pointOctalVec = outPhoton%position
    incoming = outPhoton%direction

    zAxis = VECTOR(0.,0.,1.)

    outgoing = givenVec

    randomDirection = .false.


    !
    ! Quick check to see if the location of scattering is
    ! in accretion disc. 
    if (grid%geometry == "ttauri" .and. ttau_disc_on) then
       if (in_alpha_disc(alpha_disc_param, pointOctalVec))  then
          mie_scattering = .true.
       else
          mie_scattering = .false.
       end if
    else
       ! use the value from the parameter passed to this routine
       mie_scattering = mie
    end if
       
    call amrgridvalues(grid%octreeRoot, pointOctalVec, foundOctal=thisOctal, foundSubcell=subcell)


    ! if the outgoing vector is the zero vector then this flags that
    ! we are going to scattering into a random direction

    if (modulus(outgoing) == 0.) then

       if (.not.mie) then
          randomDirection = .true.
          call random_number(r1)
          
          ! determine cos phi
          
          w = 2.d0*r1 - 1.d0
          t = sqrt(1.d0-w*w)
          call random_number(r2)
          ang = Pi * (2.*r2-1.d0)
          u = t*cos(ang)
          v = t*sin(ang)
          outgoing%x = u
          outgoing%y = v
          outgoing%z = w
       else

          outgoing = newDirectionMie(incoming, thisPhoton%lambda, lamArray, nLambda, &
               miePhase, nDustType, nMuMie, thisOctal%dustTypeFraction(subcell,:), weight)
          outPhoton%stokes = outPhoton%stokes * (1./weight)

!          outgoing = randomUnitVector()
       endif
    endif


    ! set up the scattering normals (see Hillier 1991)

    obsNormal = incoming .cross. outgoing
    call normalize(obsNormal)

    ! cos theta is the scattering angle

    costheta = incoming .dot.  outgoing

    sVec = obsNormal .cross. incoming
    call normalize(sVec)
    cosGamma = obsNormal .dot. outPhoton%normal
    sinGamma = incoming .dot. (obsNormal .cross. outPhoton%normal)
    cos2Gamma = 2.*cosGamma*cosGamma - 1.d0
    sin2Gamma = 2.*sinGamma*cosGamma

    ! rotate to new normal

    call rotate(outPhoton,cos2Gamma, sin2Gamma)



    if (.not.thisPhoton%resonanceLine) then
       if (.not.mie_scattering) then
          
          ! set up the rayleigh phase matrix and apply it

          rayleighPhase = fillRayleigh(costheta)
          outPhoton%stokes = apply(rayleighPhase, outPhoton%stokes)
       else
          
          ! do the same for a mie scattering
          
          call locate(grid%lamArray, nLambda, outPhoton%lambda, i)
          call locate(cosArray, nMuMie, costheta, j)
          if (j < nMuMie) then
             fac = (cosTheta - cosArray(j))/(cosArray(j+1)-cosArray(j))
          else
             j = nMumie
             fac =1.
          endif
          

          do k = 1, nDustType
             miePhaseTemp(k) = miePhase(k, i, j) + fac * &
                  (miePhase(k, i, j+1) - miePhase(k, i, j))
          enddo

          outPhoton%stokes = applyMean(miePhaseTemp(1:nDustType), &
               thisOctal%dustTypeFraction(subcell,1:nDustType), nDustType, outPhoton%stokes)

          
       endif
    else
       outPhoton%stokes%q = 0.
       outPhoton%stokes%u = 0.
       outPhoton%stokes%v = 0.
    endif



    if (.not.randomDirection) then

       ! since this is a scattering towards to observer we have
       ! to rotate the polarization so that is measured w.r.t. to 
       ! reference direction

       refNormal = zAxis .cross. ((-.1)*outgoing)
       call normalize(refNormal)
       cosGamma = refNormal .dot. obsNormal
       sinGamma = outgoing .dot. (refNormal .cross. obsNormal)
       cos2Gamma = 2.*cosGamma*cosGamma - 1.d0
       sin2Gamma = 2.*sinGamma*cosGamma
       call rotate(outPhoton, cos2Gamma, sin2Gamma)
       outPhoton%normal = refNormal .cross. obsNormal
    else

       ! otherwise we need to set up the scattering normal

       outPhoton%normal = sVec .cross. incoming
       call normalize(outPhoton%normal)
       cosGamma = outPhoton%normal .dot. outPhoton%oldNormal
       sinGamma = outgoing .dot. (outPhoton%normal .cross. outPhoton%oldNormal)
       cos2Gamma = 2.*cosGamma*cosGamma - 1.d0
       sin2Gamma = 2.*sinGamma*cosGamma
       outPhoton%oldNormal = outPhoton%normal
    endif

    outPhoton%direction = outgoing


    ! now find the position of the scatter in the grid and adjust the velocity
    
    if (.not.grid%doRaman) then
       
       if (grid%adaptive) then
          pointOctalVec = outPhoton%position
          if (.not.grid%resonanceLine) then
                  
             outPhoton%velocity = amrGridVelocity(grid%octreeRoot,pointOctalVec, &
                         foundOctal=octalLocation,foundSubcell=subcellLocation) 


             if (.not.mie_scattering) then
                outPhoton%velocity = outPhoton%velocity + thermalElectronVelocity( &
                     amrGridTemperature(grid%octreeRoot,pointOctalVec,&
                     startOctal=octalLocation,actualSubcell=subcellLocation))

             endif
             
          else
             outPhoton%velocity = amrGridVelocity(grid%octreeRoot,pointOctalVec)
          end if
       else
          call getIndices(grid, outPhoton%position, i1, i2, i3, t1, t2, t3)

          if (.not.grid%resonanceLine) then
             outPhoton%velocity = interpGridVelocity(grid,i1,i2,i3,real(t1),real(t2),real(t3)) + &
                  thermalElectronVelocity(grid%temperature(i1,i2,i3))
          else
             outPhoton%velocity = interpGridVelocity(grid,i1,i2,i3,real(t1),real(t2),real(t3))
          endif
       endif

       if (.not.mie_scattering) then
          vray = (outPhoton%velocity-thisPhoton%velocity) .dot. incoming
          vovercsqr = (outPhoton%velocity-thisPhoton%velocity) .dot. &
               (outPhoton%velocity-thisPhoton%velocity)
          fac = (1.d0 - 0.5d0*vovercsqr*(1.d0-0.25d0*vovercsqr))/(1.d0 + vray)
          outPhoton%lambda = outPhoton%lambda  / fac
       endif


    else

       outPhoton%redPhoton = .true.
       
       if (grid%adaptive) then
          pointOctalVec = outPhoton%position
          outPhoton%Velocity = amrGridVelocity(grid%octreeRoot,pointOctalVec)
       else
          call getIndices(grid, outPhoton%position, i1, i2, i3, t1, t2, t3)
          outPhoton%velocity = interpGridVelocity(grid,i1,i2,i3,real(t1),real(t2),real(t3))
       endif

       vray = (outPhoton%velocity-thisPhoton%velocity) .dot. incoming
       vovercsqr = (outPhoton%velocity-thisPhoton%velocity) .dot. &
            (outPhoton%velocity-thisPhoton%velocity)
       fac = (1.d0 - 0.5*vovercsqr*(1.d0-0.25*vovercsqr))/(1.d0 + vray)

       outPhoton%lambda = outPhoton%lambda / fac

       nuRest = cSpeed / (outPhoton%lambda*angstromtoCm)
       nuRest = nuRest - 2.466067749E15
       vray = outPhoton%velocity .dot. outgoing
       nuRest = nuRest * (1.d0 + vRay)
       outPhoton%redlambda = 1.e8 * cSpeed / nuRest   ! in Angs


    endif 

       
  end subroutine scatterPhoton

  ! this subroutine initializes a photon



  subroutine initPhoton(thisPhoton, Grid, nLambda, lambda, sourceSpectrum, &
       lineEmission, lamLine,  &
       weightLinePhoton, weightContPhoton, contPhoton, flatspec, vRot, &
       pencilBeam, secondSource, secondSourcePosition, &
       ramanSourceVelocity, vo6, contWindPhoton, directionalWeight,useBias, &
       theta1,theta2, chanceHotRing, &
       nSpot, chanceSpot, thetaSpot, phiSpot, fSpot, spotPhoton, probDust, weightDust, weightPhoto, &
       narrowBandImage, narrowBandMin, narrowBandMax, source, nSource, rHatInStar, energyPerPhoton, &
       filterSet, mie, curtains, starSurface, forcedWavelength, usePhotonWavelength, iLambdaPhoton,&
       VoigtProf, dirObs,  photonFromEnvelope, dopShift)
    use input_variables, only : photoionization

    implicit none

    integer :: nSource, thisSource
    type(SOURCETYPE) :: source(:)
    type(SURFACETYPE) :: starSurface
    type(PHOTON) :: thisPhoton                 ! the photon
    type(GRIDTYPE) :: grid                     ! the opacity grid
    integer :: nLambda, iLambda                ! wavelength indices
    logical :: useBias
    real :: weightContPhoton, weightLinePhoton ! the photon weights
    integer :: iLambdaPhoton
    real, optional :: dopShift
    real :: vTherm
    logical :: contWindPhoton                  ! is this continuum photon produced in the wind
    logical :: ok
    logical :: narrowBandImage
    real :: narrowBandMin, narrowBandMax  ! parameters for a narrow band image
    real(double) :: energyPerPhoton
    real :: vo6
    real :: x,y,z
    real(oct) :: xOctal, yOctal, zOctal
    type(octalVector) :: octalCentre, octVec
    type(OCTAL), pointer :: thisOctal
    type(filter_set) :: filterSet
    real :: directionalWeight
    real :: lambda(:), sourceSpectrum(:)       ! wavelength array/spectrum
    real :: dlam(1000)
    real :: r1,r2,r3                              ! radii
    real(oct)  :: r1_oct,r2_oct,r3_oct ! radii
    real(double) :: randomDouble         ! a real(double) random number
    real :: u, v, w, t                            ! direction components
    real :: r, mu, phi                         ! spherical polar coords
    real :: sinTheta      
    real(oct) :: t1, t2, t3                         ! multipliers
    integer, parameter :: nv = 100
!    real :: dv(nv)
!    real :: vArray(nv), pv(nv), vbias(nv)
    real :: vPhi                               ! azimuthal velocities
    real(double) :: kabs
    real :: ang                                ! angle
    logical :: contPhoton                      ! is this a continuum photon?
    logical :: lineEmission                    ! is there line emission?
    real :: lamLine                            ! wavelength of the line
    logical :: flatspec                        ! is the spectrum flat
    real :: vRot                               ! rotational velocity
    integer :: i1, i2, i3                      ! position indices
    type(VECTOR) :: rHat, perp                 ! radial unit vector
    type(OCTALVECTOR) :: rHatInStar
    type(VECTOR) :: rotatedVec                 ! rotated vector
    logical :: pencilBeam                      ! beamed radiation
    type(VECTOR), parameter :: zAxis = VECTOR(0.0,0.0,1.0) ! the z axis
    logical :: secondSource                    ! second photon source?
    type(VECTOR) :: secondSourcePosition       ! the position of it
    type(VECTOR) :: ramanSourceVelocity        ! what it says
    type(VECTOR) :: rVel
    type(octalVector) :: octalPoint            ! 
    type(octal), pointer :: sourceOctal        ! randomly selected octal
!    type(octal), pointer :: foundOctal       
    integer :: subcell
    real :: probDust, weightDust, weightPhoto


    logical :: oneProb = .false.

    ! Spot stuff
  
    integer :: nSpot                       ! number of spots
    real :: fSpot                          ! factional area coverage of spots
    real :: thetaSpot, phiSpot             ! spot coords
    logical :: spotPhoton
  
    real :: fac

    type(VECTOR) :: rSpot, tVec

    real :: chanceSpot
    real :: maxTheta
    real :: cosThisTheta
    real :: rotAngle
    real :: tot
    real :: weight

    real(double) :: biasWeight

    logical :: forcedWavelength
    real :: usePhotonWavelength
    logical, intent(in) ::  VoigtProf
    type(vector), intent(in) :: dirObs  ! direction to an observer

!    real :: tempr

    real :: theta1, theta2                     ! defines hot ring of accretion for TTaus
    real :: chanceHotRing                      ! chance of core photon in ttauri accretion ring
    logical, intent(in) :: curtains  ! T Tauri accretion curtains
!    real :: curtain1size             ! angular size of first accretion curtain
    real :: thisTheta, thisPhi
    real :: x1, x2
    real :: tempXProbDist(2000)
    real :: tempYProbDist(2000)
    real :: tempZProbDist(2000)
    real :: temprProbDistLine(2000)
    real :: tempmuProbDistLine(2000)
    real :: tempphiProbDistLine(2000)

    logical :: photonFromEnvelope, mie
    real(double) :: tempSpectrum(2000), prob(5000), weightarray(5000)
    real(double) :: rd, bias(5000),totDouble, lambias(5000), dlambias(5000)
    integer :: i, nbias

!    type(octalVector) :: positionOctal     ! octalVector type version of thisPhoton%position

    type(octalVector) :: octalvec_tmp
    type(Vector) :: vec_tmp
    ! For Voigt Profile 
    real :: temperature,  N_HI
    real(double) :: rho
    real(double) :: nu_shuffled, lambda_shuffled, nu, Gamma, Ne
    type(Vector) ::  velocity
!    type(octal), pointer :: thisOctal

    ! set up the weights and the stokes intensities (zero at emission)

    photonFromEnvelope = .false.

    thisPhoton%resonanceLine = .false.

    thisPhoton%velocity = VECTOR(0.,0.,0.)

    if (grid%resonanceLine) thisPhoton%resonanceLine = .true.

    vo6 = 0.
    ramanSourceVelocity = VECTOR(0.,0.,0.)

    directionalWeight = 1.
    biasWeight = 1.
    thisPhoton%fromstar1 = .false.
    thisPhoton%fromstar2 = .false.

    thisPhoton%weight = 1.d0
    thisPhoton%stokes%i = 1.d0 * energyPerPhoton

    thisPhoton%scattered = .false.
    thisPhoton%thermal = .false.
    thisPhoton%stellar = .false.

    thisPhoton%stokes%q = 0.
    thisPhoton%stokes%u = 0.
    thisPhoton%stokes%v = 0.

    ! if we are doing this by sources then find out which source we are using
    
    if (nSource > 0) then
       call randomSource(source, nSource, thisSource)
    endif

    do i = 2, nLambda-1
       dlam(i) = 0.5*((lambda(i+1)+lambda(i))-(lambda(i)+lambda(i-1)))
    enddo
    dlam(1) = lambda(2)-lambda(1)
    dlam(nLambda) = lambda(nlambda)-lambda(nLambda-1)


    ! is this a raman photon
    if (grid%doRaman) thisphoton%redPhoton = .false.

    ! if there is line emission then set up the appropriate flags

    if (lineEmission.or.grid%doRaman) then

       if (.not.contPhoton) then
          thisPhoton%linePhoton = .true.
          thisPhoton%contPhoton = .false.
       else
          thisPhoton%linePhoton = .false.
          thisPhoton%contPhoton = .true.
       endif

    else

       ! of course if there is no line emission then all photons
       ! are continuum ones

       thisPhoton%linePhoton = .false.
       thisPhoton%contPhoton = .true.

    endif

    contwindphoton = .false.

    if (grid%lineEmission) then
       call random_number(r)
       if (r < grid%chanceWindOverTotalContinuum) then
          contWindPhoton = .true.
       endif
    endif


    if (mie) then
       call random_number(r)
       if (r < probDust) then
          photonFromEnvelope = .true.
          contWindPhoton = .true.
          thisPhoton%stokes = thisPhoton%stokes * weightDust
       else
          thisPhoton%stokes = thisPhoton%stokes * weightPhoto
       endif
    endif



    ! if it is a continuum photon initialize its position

    if (thisPhoton%contPhoton) then


       if ((grid%cartesian.or.grid%adaptive) .and.     &
                         (.not.grid%lineEmission)) then

          thisPhoton%position = grid%rCore * randomUnitVector()
          
          if (nSource > 0) then
             call getPhotonPositionDirection(source(thisSource), thisPhoton%position, thisPhoton%direction, rHatInStar)
             thisPhoton%stellar = .true.
          endif

                  
          if (photonFromEnvelope) then
             thisPhoton%direction = randomUnitVector()
             thisPhoton%thermal = .true.
             if (grid%adaptive) then 

                do  ! dummy loop, in case we pick a position inside a star

                  call random_number(randomDouble)
                  ! we search through the tree to find the subcell that contains the
                  !   probability value 'randomDouble'
                  sourceOctal => grid%octreeRoot
                  call locateContProbAMR(randomDouble,sourceOctal,subcell)
                  if (.not.sourceOctal%inFlow(subcell)) then
                    write(*,'(a)') "! Photon in cell that's not in flow. Screw-up in locatecontProbAmr"
                    stop
                  endif


                  thisPhoton%position = randomPositionInCell(sourceOctal, subcell)

!                  octalCentre = subcellCentre(sourceOctal,subcell)
!
!                  !!! we will just choose a random point within the subcell.
!                  !!! this *should* be done in a better way.
!
!                  call random_number(r1)
!                  r1 = r1 - 0.5  ! shift value mean value to zero
!                  r1 = r1 * 0.9999 ! to avoid any numerical accuracy problems
!                  xOctal = r1 * sourceOctal%subcellSize + octalCentre%x
!
!                  call random_number(r2)
!                  r2 = r2 - 0.5                                  
!                  r2 = r2 * 0.9999                                          
!                  yOctal = r2 * sourceOctal%subcellSize + octalCentre%y
!
!                
!                call random_number(r3)
!                r3 = r3 - 0.5                                  
!                r3 = r3 * 0.9999                                          
!                zOctal = r3 * sourceOctal%subcellSize + octalCentre%z
!                
!                thisPhoton%position = vector(xOctal,yOctal,zOctal)
!
!                if (sourceOctal%twod) then
!                   call random_number(ang)
!                   ang = ang * twoPi
!                   thisPhoton%position = rotateZ(thisPhoton%position, dble(ang))
!                endif

                if (grid%geometry(1:7) == "ttauri"    .or.  &
                    grid%geometry(1:9) == "luc_cir3d" .or.  &
                    grid%geometry(1:6) == "cmfgen"    .or.  &
                    grid%geometry(1:8) == "romanova")  then
                   ! need to check the position is not inside the star
                   if ((modulus(thisPhoton%position-(s2o(grid%starPos1)))) > grid%rStar1) then
                      exit
                   else
                      continue ! pick another one
                   endif
                else if (inOctal(grid%octreeRoot, thisPhoton%position)) then
                   exit
                end if

                end do 

                !!! need to call an interpolation routine, rather than
                !!!   use subcell central value
                if (useBias) then
                   biasWeight = biasWeight * 1.0_db / sourceOctal%biasCont3D(subcell)
                else
                   biasWeight = 1.0_db
                end if

             else ! grid is cartesian
                ok = .false.
                do while(.not.ok)
                   call random_number(r1)
                   tempXprobdist(1:grid%nx) = grid%xProbDistCont(1:grid%nx)
                   call locate(grid%xProbDistCont, grid%nx, r1, i1)
                   t1 = (r1 - grid%xProbDistCont(i1))/(grid%xProbDistCont(i1+1)-grid%xProbDistCont(i1))
                   x = grid%xAxis(i1) + t1 * (grid%xAxis(i1+1)-grid%xAxis(i1))


                   tempYProbDist(1:grid%ny) = grid%yProbDistCont(i1,1:grid%ny) + t1 * &
                        (grid%yProbDistCont(i1+1,1:grid%ny) - grid%yProbDistCont(i1,1:grid%ny))
                   tempYProbDist(1:grid%ny) = tempYprobDist(1:grid%ny) / tempYprobDist(grid%ny)

                   call random_number(r2)
                   call locate(tempYProbDist, grid%ny, r2, i2)
                   t2 = (r2 - tempYProbDist(i2)) / &
                        (tempYProbDist(i2+1) - tempYProbDist(i2))
                   y = grid%yAxis(i2) + t2 * (grid%yAxis(i2+1)-grid%yAxis(i2))

                   tempZProbDist(1:grid%nz) = &
                        (1.d0-t1)*(1.d0-t2) * grid%zProbDistCont(i1  , i2  , 1:grid%nz) +&
                        (   t1)*(1.d0-t2) * grid%zProbDistCont(i1+1, i2  , 1:grid%nz) +&
                        (1.d0-t1)*(   t2) * grid%zProbDistCont(i1  , i2+1, 1:grid%nz) +&
                        (   t1)*(   t2) * grid%zProbDistCont(i1+1, i2+1, 1:grid%nz)
                   tempZProbDist(1:grid%nz) = tempZprobDist(1:grid%nz) / tempZprobDist(grid%nz)


                   call random_number(r3)
                   call locate(tempZProbDist, grid%nz, r3, i3)
                   t3 = (r3 - tempZProbDist(i3)) / &
                        (tempZProbDist(i3+1) - tempZProbDist(i3))
                   z = grid%zAxis(i3) + t3 * (grid%zAxis(i3+1)-grid%zAxis(i3))


                   thisPhoton%position  = VECTOR(x,y,z)

                   if (outSideGrid(thisPhoton%position,grid).or.(.not.grid%inUse(i1,i2,i3))) then
                      !                      write(*,*) "Mistake in initPhoton",thisPhoton%position
                      !                      write(*,*) i1,i2,i3,t1,t2,t3
                      ok = .false.
                   else
                      ok = .true.
                   endif


                enddo
                biasWeight =  biasWeight * &
                  (1.d0/ interpGridScalar2(grid%biasCont3d,grid%nx,grid%ny,grid%nz,i1,i2,i3,real(t1),real(t2),real(t3)))

             end if ! cartesian or adaptive
          
          end if ! (photonFromEnvelope)
          
       else ! grid has polar coordinates  or grid%lineEmission==.true.

          ! if there is line emission, then the continuum photons can be produced either
          ! at the core, or in the wind itself.

          if (contWindPhoton) then
             if (grid%cartesian) then
                ok = .false.
                do while(.not.ok)
                   call random_number(r1)
                   tempXprobdist(1:grid%nx) = grid%xProbDistCont(1:grid%nx)
                   call locate(grid%xProbDistCont, grid%nx, r1, i1)
                   t1 = (r1 - grid%xProbDistCont(i1))/(grid%xProbDistCont(i1+1)-grid%xProbDistCont(i1))
                   x = grid%xAxis(i1) + t1 * (grid%xAxis(i1+1)-grid%xAxis(i1))


                   tempYProbDist(1:grid%ny) = grid%yProbDistCont(i1,1:grid%ny) + t1 * &
                        (grid%yProbDistCont(i1+1,1:grid%ny) - grid%yProbDistCont(i1,1:grid%ny))
                   tempYProbDist(1:grid%ny) = tempYprobDist(1:grid%ny) / tempYprobDist(grid%ny)

                   call random_number(r2)
                   call locate(tempYProbDist, grid%ny, r2, i2)
                   t2 = (r2 - tempYProbDist(i2)) / &
                        (tempYProbDist(i2+1) - tempYProbDist(i2))
                   y = grid%yAxis(i2) + t2 * (grid%yAxis(i2+1)-grid%yAxis(i2))

                   tempZProbDist(1:grid%nz) = &
                        (1.d0-t1)*(1.d0-t2) * grid%zProbDistCont(i1  , i2  , 1:grid%nz) +&
                        (   t1)*(1.d0-t2) * grid%zProbDistCont(i1+1, i2  , 1:grid%nz) +&
                        (1.d0-t1)*(   t2) * grid%zProbDistCont(i1  , i2+1, 1:grid%nz) +&
                        (   t1)*(   t2) * grid%zProbDistCont(i1+1, i2+1, 1:grid%nz)
                   tempZProbDist(1:grid%nz) = tempZprobDist(1:grid%nz) / tempZprobDist(grid%nz)


                   call random_number(r3)
                   call locate(tempZProbDist, grid%nz, r3, i3)
                   t3 = (r3 - tempZProbDist(i3)) / &
                        (tempZProbDist(i3+1) - tempZProbDist(i3))
                   z = grid%zAxis(i3) + t3 * (grid%zAxis(i3+1)-grid%zAxis(i3))


                   thisPhoton%position  = VECTOR(x,y,z)

                   if (outSideGrid(thisPhoton%position,grid).or.(.not.grid%inUse(i1,i2,i3))) then
                      !                      write(*,*) "Mistake in initPhoton",thisPhoton%position
                      !                      write(*,*) i1,i2,i3,t1,t2,t3
                      ok = .false.
                   else
                      ok = .true.
                   endif


                enddo
                biasWeight =  biasWeight * &
                     (1.d0/ interpGridScalar2(grid%biasCont3d,grid%nx,grid%ny,grid%nz,i1,i2,i3,real(t1),real(t2),real(t3)))
                     
             elseif (grid%adaptive) then
                   
                do ! dummy loop, in case we pick a position inside a star
                  call random_number(randomDouble)

                  ! we search through the tree to find the subcell that contains the
                  !   probability value 'randomDouble'
                  sourceOctal => grid%octreeRoot
                  call locateContProbAMR(randomDouble,sourceOctal,subcell)


!                  octalCentre = subcellCentre(sourceOctal,subcell)
!
!                  !!! we will just choose a random point within the subcell.
!                  !!! this *should* be done in a better way.
!
!                  call random_number(r1)
!                  r1 = r1 - 0.5  ! shift value mean value to zero
!                  r1 = r1 * 0.995 ! to avoid any numerical accuracy problems
!                  xOctal = r1 * sourceOctal%subcellSize + octalCentre%x
!
!                  call random_number(r2)
!                  r2 = r2 - 0.5                                  
!                  r2 = r2 * 0.995                                           
!                  yOctal = r2 * sourceOctal%subcellSize + octalCentre%y
!
!                  call random_number(r3)
!                  r3 = r3 - 0.5                                  
!                  r3 = r3 * 0.995                                           
!                  zOctal = r3 * sourceOctal%subcellSize + octalCentre%z
!
!                octalCentre = subcellCentre(sourceOctal,subcell)
!                
!                !!! we will just choose a random point within the subcell.
!                !!! this *should* be done in a better way.
!                
!                call random_number(r1)
!                r1 = r1 - 0.5  ! shift value mean value to zero
!                r1 = r1 * 0.995 ! to avoid any numerical accuracy problems
!                xOctal = r1 * sourceOctal%subcellSize + octalCentre%x
!                
!                if (sourceOctal%threed) then
!                   call random_number(r2)
!                   r2 = r2 - 0.5                                  
!                   r2 = r2 * 0.995                                           
!                   yOctal = r2 * sourceOctal%subcellSize + octalCentre%y
!                else
!                   yOctal = 0.
!                endif
!                
!                call random_number(r3)
!                r3 = r3 - 0.5                                  
!                r3 = r3 * 0.995                                           
!                zOctal = r3 * sourceOctal%subcellSize + octalCentre%z
!                
!                thisPhoton%position = vector(xOctal,yOctal,zOctal)
!                                
!                if (sourceOctal%twod) then
!                   call random_number(ang)
!                   ang = ang * twoPi
!                   thisPhoton%position = rotateZ(thisPhoton%position, dble(ang))
!                endif


                thisPhoton%position = randomPositionInCell(sourceOctal, subcell)



                  if (grid%geometry(1:7) == "ttauri"    .or. &
                      grid%geometry(1:9) == "luc_cir3d" .or. &
                      grid%geometry(1:6) == "cmfgen"    .or. &
                      grid%geometry(1:8) == "romanova"          ) then
                    ! need to check the position is not inside the star
                    if ((modulus(thisPhoton%position-s2o(grid%starPos1))) > grid%rStar1) exit
                  else if (.not. inOctal(grid%octreeRoot, thisPhoton%position)) then
                    exit
                  end if

                end do 

                !!! need to call an interpolation routine, rather than
                !!!   use subcell central value
                if (useBias) then
                   biasWeight = biasWeight * 1.0_db / sourceOctal%biasCont3D(subcell)
                else
                   biasWeight = 1.0_db
                end if
               
             else ! grid is polar

980             continue
                call random_number(r1)
                call locate(grid%rProbDistCont, grid%nr, r1, i1)
                t1 = (r1-grid%rProbDistCont(i1))/(grid%rProbDistCont(i1+1)-grid%rProbDistCont(i1))
                r = grid%rAxis(i1) + t1 * (grid%rAxis(i1+1)-grid%rAxis(i1))

                call random_number(r2)
                call locate(grid%muProbDistCont(i1,1:grid%nmu), grid%nmu, r2, i2)
                t2 = (r2-grid%muProbDistCont(i1,i2))/(grid%muProbDistCont(i1,i2+1)-grid%muProbDistCont(i1,i2))
                mu = grid%muAxis(i2) + t2 * (grid%muAxis(i2+1)-grid%muAxis(i2))

                call random_number(r3)
                call locate(grid%phiProbDistCont(i1,i2,1:grid%nphi), grid%nphi, r3, i3)
                t3 = (r3-grid%phiProbDistCont(i1,i2,i3))/(grid%phiProbDistCont(i1,i2,i3+1)-grid%phiProbDistCont(i1,i2,i3))
                phi = grid%phiAxis(i3) + t3 * (grid%phiAxis(i3+1)-grid%phiAxis(i3))

                if (.not.grid%inUse(i1,i2,i3)) then
                   goto 980
                endif


                x = r * sqrt(1.d0-mu*mu) * cos(phi)
                y = r * sqrt(1.d0-mu*mu) * sin(phi)
                z = r * mu
                if (useBias.and.(.not.grid%cartesian)) then
                   biasWeight = interpGridScalar2(grid%biasCont3d,grid%nr,grid%nmu,grid%nphi,i1,i2,i3,real(t1),real(t2),real(t3))
!                   biasWeight = grid%biasCont3d(i1,i2,i3)
                   biasWeight = 1.d0/biasWeight
                endif

                thisPhoton%position = VECTOR(x, y, z)

             endif



          else ! not (contWindPhoton)


             r = grid%rCore*1.0001d0
             select case(grid%geometry)

                case("disk")
                   if (nSpot > 0) then
                      rSpot = VECTOR(cos(phiSpot)*sin(thetaSpot),sin(phiSpot)*sin(thetaSpot),cos(thetaSpot))
                      if (nSpot == 1) then
                         maxTheta = Pi * fSpot   
                      else
                         maxTheta = Pi * fSpot * 0.5
                      endif

                      spotPhoton = .false.

                      call random_number(r1)
                      if (r1 < chanceSpot) then
                         spotPhoton = .true.
                         call random_number(r1)
                         cosThisTheta = r1 * (1. - cos(maxTheta)) + cos(maxTheta)    
                         thisTheta = acos(cosThisTheta)
                         call random_number(r1)
                         thisPhi = twoPi * r1
                         rHat = VECTOR(cos(thisPhi)*sin(thisTheta),sin(thisPhi)*sin(thisTheta),cos(thisTheta))
                         
                         tVec = zAxis .cross. rSpot  
                         call normalize(tVec)
                         rotAngle = zAxis .dot. rSpot
                         rotAngle = acos(rotAngle)
                         rHat = arbitraryrotate(rHat,rotAngle,tVec)
                         if (nSpot == 2) then
                            call random_number(r1)
                            if (r1 < 0.5) then
                               rHat = (-1.) * rHat
                               rSpot = (-1.) * rSpot
                            endif
                         endif
                         thisPhoton%position = r * rHat
                      else
                         tVec = randomUnitVector()
                         ang = tVec .dot. rSpot
                         ang = acos(ang)
                         do while (ang < maxTheta)
                            tVec = randomUnitVector()
                            ang = tVec .dot. rSpot
                            ang = acos(ang)
                         enddo
                         thisPhoton%position = r*tVec
                      endif
                   else
                      thisPhoton%position = (r*randomUnitVector())
                   endif

                case("binary")
                   call random_number(r)
                   if (r < grid%lumRatio) then
                      thisPhoton%position = grid%starPos1 + (1.01*(grid%rStar1) * randomUnitVector())
                      thisPhoton%fromStar1 = .true.
                      !call getIndices(grid,thisPhoton%position,i1,i2,i3,t1,t2,t3)
                   else
                      thisPhoton%position = grid%starPos2 + (1.01*(grid%rStar2) * randomUnitVector())
                      thisPhoton%fromStar2 = .true.
                      !call getIndices(grid,thisPhoton%position,i1,i2,i3,t1,t2,t3)
                   endif

                case("ttauri", "luc_cir3d",  "cmfgen", "romanova")

!                   call random_number(r1)
!                   if (r1 < chanceHotRing) then
!                      call random_number(r1)
!                      thisTheta = r1*(theta2 - theta1)+theta1
!                      call random_number(r1)
!                      if (r1 < 0.5) thisTheta = pi - thisTheta
!                      call random_number(r1)
!                      if (curtains) then 
!                        ! if we have accretion curtains, we have to restrict the
!                        ! emission region.
!                        curtain1size = curtainsPhi1e - curtainsPhi1s
!                        thisPhi =  r1 * (curtain1size+(curtainsPhi2e-curtainsPhi2s))    &
!                                     + curtainsPhi1s
!                        if (thisPhi > curtainsPhi1e) then
!                          thisPhi = thisPhi + (curtainsPhi2s-curtainsPhi1e)
!                        end if
!                      else
!                        thisPhi = twoPi * r1
!                      end if
!                      rHat = VECTOR(cos(thisPhi)*sin(thisTheta),sin(thisPhi)*sin(thisTheta),cos(thisTheta))
!                      thisPhoton%position = r * rHat
!                   else
!                      thisPhoton%position = (r*randomUnitVector())
!                   endif


                   call getPhotoVec(starSurface, thisPhoton%position, thisPhoton%direction)


                case("jets")
                   ! emission from the surface
!                   r = get_jets_parameter("Rmin")
                   r = grid%rStar1
                   rHat = vector(x,y,z)
                   thisPhoton%position = r*randomUnitVector()

                case DEFAULT
                   thisPhoton%position = (r*randomUnitVector())
             end select


             thisPhoton%originalNormal = thisPhoton%position
             call normalize(thisPhoton%originalNormal)

          endif ! (contWindPhoton)
          
       endif ! grid cartesian/adaptive/polar

       ! get wavelength
       
       if (narrowBandImage) then
          call createWeightArrays(filterSet, Lambias, dlambias, nbias, bias, weightArray, dble(lambda(1)), dble(lambda(nLambda)))
          prob(1) = 0.
          do i = 2, nbias
             prob(i) = prob(i-1) + bias(i)
          enddo
          prob(1:nbias) = prob(1:nbias) / prob(nbias)
          call random_number(rd)
          call locate(prob, nbias, rd, i)
          t = (rd-prob(i))/(prob(i+1)-prob(i))
!          do i = 1, nbias
!             write(99,*) lambias(i),prob(i)
!          enddo
!          stop
          thisPhoton%lambda = lambias(i) + t * (lambias(i+1)-lambias(i))
          thisPhoton%stokes = thisPhoton%stokes * real(weightArray(i) +t*(weightArray(i+1)-weightArray(i)))
          call locate(lambda, nLambda, thisPhoton%lambda, ilambda)
       else
          if (mie) then
             tot = SUM(dlam(1:nLambda))
             call random_number(r1)
             if (forcedWavelength) then
                call locate(lambda, nlambda, usePhotonWavelength, ilambda)
             else
                iLambda = int(r1 * real(nLambda)) + 1
             endif

             iLambda = iLambdaPhoton

             weight =  dlam(iLambda) !tjh 21/3/07
 

             thisPhoton%lambda = lambda(ilambda)
             thisPhoton%stokes = thisPhoton%stokes * weight ! * real(nLambda) tjh 21/3/07
             call random_number(r)
             if (iLambda == 1) then
                x1 = lambda(1)
                x2 = 0.5*(lambda(1)+lambda(2))
             else if (iLambda == nLambda) then
                x1 = 0.5*(lambda(nLambda)+lambda(nLambda-1))
                x2 = lambda(nLambda)
             else
                x1 = 0.5*(lambda(ilambda-1)+lambda(ilambda))
                x2 = 0.5*(lambda(ilambda+1)+lambda(ilambda))
             endif
!             call random_number(r)
!          thisPhoton%lambda = x1+r*(x2-x1)
          elseif (grid%lineEmission .and. thisPhoton%contPhoton) then
             !pick line center wavelength
             thisPhoton%lambda = lamLine
          else
             call random_number(r1)
             iLambda = int(r1 * real(nLambda)) + 1
             thisPhoton%lambda = lambda(ilambda)
          endif
       endif

       
       if (.not.flatspec) then
          if (photonFromEnvelope) then

             totDouble = 0.d0
             
             if (grid%adaptive) then

                if (.not.photoionization) then
!                   positionOctal = thisPhoton%position
!                   call amrGridvalues(grid%octreeRoot,positionOctal,&
!                        foundOctal=foundOctal,foundSubcell=subcell, temperature=tempr, kappaAbs=kabs, grid=grid, iLambda=ilambda)
!                   do i = 1, nLambda
!                      call amrGridvalues(grid%octreeRoot,positionOctal,&
!                           foundOctal=foundOctal,foundSubcell=subcell, temperature=tempr, kappaAbs=kabs, grid=grid, iLambda=i)
!                      tempSpectrum(i)= blambda(dble(lambda(i)), dble(tempr)) * dble(kabs) !/ dble(lambda(i))
!                      totDouble = totDouble + tempSpectrum(i) * dlam(i)
!                   enddo
!                   if (totDouble == 0.d0) then
!                      do i = 1, nLambda
!                         call amrGridvalues(grid%octreeRoot,positionOctal,&
!                              foundOctal=foundOctal,foundSubcell=subcell, temperature=tempr, kappaAbs=kabs, grid=grid, iLambda=i)
!                         tempSpectrum(i)= blambda(dble(lambda(i)), dble(tempr)) * dble(kabs) !/ dble(lambda(i))
!                         totDouble = totDouble + tempSpectrum(i) * dlam(i)
!                         write(*,*) i,lambda(i),dlam(i),tempr,kabs,blambda(dble(lambda(i)), dble(tempr)), &
!                              foundOctal%dustTypeFraction(subcell, 1),foundOctal%etaCont(subcell)
!                      enddo
!                      totDouble = 10.
!                   endif
!
!
!                   tempSpectrum(1:nLambda) = tempSpectrum(1:nLambda) / totDouble
!
!                   if (totDouble<=0.0) totDouble = 1.0e-28   ! for safty
!
!TJH 21/3/07 - comment this out 
!                   thisPhoton%stokes = thisPhoton%stokes * real(tempSpectrum(iLambda))

                else

                   octVec = thisPhoton%position
                   call amrgridvalues(grid%octreeRoot, octVec, foundOctal=thisOctal, foundSubcell=subcell)
                   call  getWavelengthBiasPhotoion(grid, thisOctal, subcell, lambda, dlam, nLambda, ilambda, fac, .true.)
                   thisPhoton%lambda = lambda(ilambda)
                   thisPhoton%stokes = thisPhoton%stokes * fac

                  
                endif


!                thisPhoton%stokes = thisPhoton%stokes * interplinearSingle(real(lambda), &
!                     real(tempSpectrum), nLambda, thisPhoton%lambda)

             else 
                do i = 1, nLambda
                   if (.not.grid%oneKappa) then
                      kabs = grid%kappaAbs(i1,i2,i3,iLambda)
                   else
                      kabs = grid%oneKappaAbs(1,iLambda) * grid%rho(i1,i2,i3)
                   endif
                   tempSpectrum(i) =  blambda(dble(lambda(i)),dble(grid%temperature(i1,i2,i3))) &
                         *  dble(kabs) 
                   totDouble = totDouble + tempSpectrum(i)
                enddo
                tempSpectrum(1:nLambda)  = tempSpectrum(1:nLambda) / totDouble

                thisPhoton%stokes = thisPhoton%stokes * interplinearSingle(real(lambda), &
                     real(tempSpectrum), nLambda, thisPhoton%lambda)
             end if
          else

             if (nSource == 0) then
                thisPhoton%stokes = thisPhoton%stokes * &
                     (sourceSpectrum(iLambda) * weightContPhoton)
             else
!TJH 21/3/07 - comment this out 
!               thisPhoton%stokes = thisPhoton%stokes * &
 !                    real(returnNormValue2(source(thissource)%spectrum, &
 !                    dble(thisPhoton%lambda), dble(lambda(1)), dble(lambda(nlambda))))
             endif
          endif
          
       else
          thisPhoton%stokes = thisPhoton%stokes * weightContPhoton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!          thisPhoton%lambda = lamLine
       endif


       ! get direction

       if (.not.pencilBeam) then
          r3_oct = modulus(thisPhoton%position)
          if (r3_oct /= 0.) then
             if (.not.grid%geometry == "binary") then
                rHat = thisPhoton%position / r3_oct
             else
                if (thisPhoton%fromStar1) then
                   octalvec_tmp = grid%starPos1
                   rHat = (thisPhoton%position-octalvec_tmp) 
                   call normalize(rHat)
                   thisPhoton%originalNormal = rHat
                else
                   octalvec_tmp = grid%starPos2
                   rHat = (thisPhoton%position-octalvec_tmp) 
                   call normalize(rHat)
                   thisPhoton%originalNormal = rHat
                endif
             endif
          endif

! the line below introduced a bug for photospheric photons
! whose direction had been set earlier via 
! a call to getPhotonPositionDirection

!          thisPhoton%direction = randomUnitVector()

       else
          call random_number(r1)
          w = 1.d0 - r1*sin(5.*degToRad)
          call random_number(r2)
          if (r2 < 0.5) w = -w
          t = sqrt(1.d0-w*w)
          call random_number(r2)
          ang = pi*(2.*r2-1.d0)
          u = t*cos(ang)
          v = t*sin(ang)
          thisPhoton%direction%x = u
          thisPhoton%direction%y = v
          thisPhoton%direction%z = w
          vec_tmp = thisPhoton%direction
          rotatedVec = rotateY(vec_tmp, 30.*real(degToRad))
          thisPhoton%direction = rotatedVec

       endif



       if (.not.contWindPhoton) then
          if (.not.grid%geometry=="rolf") then
             t = (thisPhoton%direction .dot. rHat)

             ! must be outwards from the photosphere if this is a core continuum photon

             if ((t < 0.).and.(r3_oct /= 0.).and.(grid%lineEmission).and.(.not.contWindPhoton)) then
                thisPhoton%direction = (-1.d0)*thisPhoton%direction
             endif
             t = (thisPhoton%direction .dot. rHat)
             directionalWeight = abs(2.*t)
          endif




       endif

       ! this does the photon velocity for a rotational velocity field
       ! not the local grid velocity, since the radial component must
       ! be zero

       if (grid%geometry == "binary") then
          if (grid%adaptive) then 
             octalPoint = thisPhoton%position  
             thisPhoton%velocity = amrGridVelocity(grid%octreeRoot,octalPoint)
          else  
            call getIndices(grid, thisPhoton%position, i1, i2, i3, t1, t2, t3)
            thisPhoton%velocity = grid%velocity(i1,i2,i3)
          endif
       endif


       if (.not.grid%cartesian .and. .not. grid%adaptive) then
          if (grid%rCore /= 0.) then
             vPhi = (vRot/cSpeed)*sqrt(max(0.d0,1.d0-(thisPhoton%position%z/grid%rCore)**2))
             perp = thisPhoton%position .cross. zAxis
             call normalize(perp)
             thisPhoton%velocity = vPhi*perp
          endif
       endif

       if (grid%geometry == "disk") then
          if (.not.contWindPhoton) thisPhoton%velocity = VECTOR(0.,0.,0.)
       endif


    endif ! (thisPhoton%contPhoton)


    ! now for a line photon

    if (thisPhoton%linePhoton) then

       ! for the polar grid we use the probability distributions to get
       ! r, theta and phi for the emission

       

       if (grid%cartesian) then

          ok = .false.

          if (secondSource) then
             thisPhoton%position = secondSourcePosition
             ok = .true.
             call getIndices(grid,thisPhoton%position,i1,i2,i3,t1,t2,t3)   
          endif


          do while(.not.ok)
             call random_number(r1)
             tempXprobdist(1:grid%nx) = grid%xProbDistLine(1:grid%nx)
             call locate(grid%xProbDistLine, grid%nx, r1, i1)
             t1 = (r1 - grid%xProbDistLine(i1))/(grid%xProbDistLine(i1+1)-grid%xProbDistLine(i1))
             x = grid%xAxis(i1) + t1 * (grid%xAxis(i1+1)-grid%xAxis(i1))

             if (grid%geometry == "rolf") then
                x = x + (grid%xAxis(2) - grid%xAxis(1))/2.
                call locate(grid%xAxis, grid%nx, x, i1)
                t1 = (x - grid%xAxis(i1))/(grid%xAxis(i1+1)-grid%xAxis(i1))
             endif


             tempYProbDist(1:grid%ny) = grid%yProbDistLine(i1,1:grid%ny) + t1 * &
                  (grid%yProbDistLine(i1+1,1:grid%ny) - grid%yProbDistLine(i1,1:grid%ny))
             tempYProbDist(1:grid%ny) = tempYprobDist(1:grid%ny) / tempYprobDist(grid%ny)
             call random_number(r2)
             call locate(tempYProbDist, grid%ny, r2, i2)
             t2 = (r2 - tempYProbDist(i2)) / &
                  (tempYProbDist(i2+1) - tempYProbDist(i2))
             y = grid%yAxis(i2) + t2 * (grid%yAxis(i2+1)-grid%yAxis(i2))

             if (grid%geometry == "rolf") then
                y = y + (grid%yAxis(2) - grid%yAxis(1))/2.
                call locate(grid%yAxis, grid%ny, y, i2)
                t2 = (y - grid%yAxis(i2))/(grid%yAxis(i2+1)-grid%yAxis(i2))
             endif


             tempZProbDist(1:grid%nz) = &
                  (1.d0-t1)*(1.d0-t2) * grid%zProbDistLine(i1  , i2  , 1:grid%nz) +&
                  (   t1)*(1.d0-t2) * grid%zProbDistLine(i1+1, i2  , 1:grid%nz) +&
                  (1.d0-t1)*(   t2) * grid%zProbDistLine(i1  , i2+1, 1:grid%nz) +&
                  (   t1)*(   t2) * grid%zProbDistLine(i1+1, i2+1, 1:grid%nz)
             tempZProbDist(1:grid%nz) = tempZprobDist(1:grid%nz) / tempZprobDist(grid%nz)

             call random_number(r3)
             call locate(tempZProbDist, grid%nz, r3, i3)
             t3 = (r3 - tempZProbDist(i3)) / &
                  (tempZProbDist(i3+1) - tempZProbDist(i3))
             z = grid%zAxis(i3) + t3 * (grid%zAxis(i3+1)-grid%zAxis(i3))


             if (grid%geometry == "rolf") then
                z = z + (grid%zAxis(2) - grid%zAxis(1))/2.
                call locate(grid%zAxis, grid%nz, z, i3)
                t3 = (z - grid%zAxis(i3))/(grid%zAxis(i3+1)-grid%zAxis(i3))
             endif


             thisPhoton%position%x = x
             thisPhoton%position%y = y
             thisPhoton%position%z = z


             if (outSideGrid(thisPhoton%position,grid).or.(.not.grid%inUse(i1,i2,i3))) then
                !                write(*,*) "Mistake in initPhoton",thisPhoton%position
                !                write(*,*) i1,i2,i3,t1,t2,t3
                ok = .false.
             else
                ok = .true.
             endif

          enddo
          biasWeight =  biasWeight * &
               (1.d0/ interpGridScalar2(grid%biasLine3D,grid%nx,grid%ny,grid%nz,i1,i2,i3,real(t1),real(t2),real(t3)))

       elseif (grid%adaptive) then
    
          ok = .false.
          
          if (secondSource) then
             thisPhoton%position = secondSourcePosition
             ok = .true.
          endif

          if (grid%geometry == "rolf") then
             write(*,*) "Panic: in initPhoton, there is no code for handling ",&
                        "grid%geometry==rolf with an adaptive grid"
             stop
          endif
          ok = .false.
          if (.not. ok) then

           do  ! dummy loop, in case we pick a position inside a star ===========

             ! we search through the tree to find the subcell that contains the
             !   probability value 'randomDouble'
             call random_number(randomDouble)
             sourceOctal => grid%octreeRoot
             call locateLineProbAMR(randomDouble,sourceOctal,subcell)
             if (.not.sourceOctal%inFlow(subcell)) then
                write(*,'(a)') "Photon in cell that's not in flow. Screw-up in locateLineProbAmr!"
                stop
             endif


             octalCentre = subcellCentre(sourceOctal,subcell)

            !!! need to call an interpolation routine, rather than
            !!!   use subcell central value
             if (useBias) then
                biasWeight = biasWeight * 1.0_db / sourceOctal%biasLine3D(subcell)
             else
                biasWeight = 1.0_db
             end if
             !!! we will just choose a random point within the subcell.
             !!! this *should* be done in a better way.
             call random_number(r1_oct)
             r1_oct = r1_oct - 0.5d0  ! shift value mean value to zero
             r1_oct = r1_oct * 0.995d0 ! to avoid any numerical accuracy problems
             xOctal = r1_oct * sourceOctal%subcellSize + octalCentre%x
               
             if (sourceOctal%threed) then
                call random_number(r2_oct)
                r2_oct = r2_oct - 0.5d0                                  
                r2_oct = r2_oct * 0.995d0                                           
                yOctal = r2_oct * sourceOctal%subcellSize + octalCentre%y
             else
                yOctal = 0.0d0
             endif

                
             call random_number(r3_oct)
             r3_oct = r3_oct - 0.5d0                                  
             r3_oct = r3_oct * 0.995d0                                           
             zOctal = r3_oct * sourceOctal%subcellSize + octalCentre%z
             

!!just for testing ... 
!             thisPhoton%position = octalCentre
             thisPhoton%position = vector(xOctal,yOctal,zOctal)

             if (grid%geometry(1:7) == "ttauri"      .or.  &
                  grid%geometry(1:9) == "luc_cir3d"  .or.  &
                  grid%geometry(1:6) == "cmfgen"     .or.  &
                  grid%geometry(1:8) == "romanova"          ) then
                ! need to check the position is not inside the star
                if ((modulus(thisPhoton%position-(s2o(grid%starPos1)))) > grid%rStar1) exit
             else
                ! exit
                ! pick another one
                continue
             end if

            end do
                

             if (sourceOctal%twod) then
                call random_number(ang)
                ang = ang * twoPi
                thisPhoton%position = rotateZ(thisPhoton%position, dble(ang))
             endif

          end if ! (.not. OK)

       else ! grid has polar coordinates


          if (thisPhoton%resonanceLine) then
             r = grid%rCore*1.00001d0
             rHat = randomUnitVector()
             thisPhoton%position = r * rHat
             thisPhoton%originalNormal = thisPhoton%position
             thisPhoton%direction = randomUnitVector()
             t = (thisPhoton%direction .dot. rHat)

             ! must be outwards from the photosphere if this is a core continuum photon

             if (t < 0.) then
                thisPhoton%direction = (-1.0d0)*thisPhoton%direction
             endif
             t = (thisPhoton%direction .dot. rHat)
             directionalWeight = abs(2.*t)
             call normalize(thisPhoton%originalNormal)
             thisPhoton%velocity = VECTOR(1.e-30,1.e-30,1.e-30)
             !             call random_number(r1)
             !             thisPhoton%lambda = lambda(1)+r1*(lambda(nLambda)-lambda(1))

             ! The resonance line wavelength is given in torus main according to the photon loop
             ! in order to evenly distribute photons across the input wavelength range

          else


             if (oneProb) then
                
                call getLinePosition(grid, r, mu, phi, i1, i2, i3)
                
             else

990             continue

             call random_number(r1)
             temprProbDistLine(1:grid%nr) = grid%rProbDistLine(1:grid%nr)
             call locate(temprProbDistLine, grid%nr, r1, i1)
             t1 = (r1-temprProbDistLine(i1))/(temprProbDistLine(i1+1)-temprProbDistLine(i1))
             r = grid%rAxis(i1) + t1 * (grid%rAxis(i1+1)-grid%rAxis(i1))


             call random_number(r1)
             tempMuProbDistLine(1:grid%nMu) = grid%muProbDistLine(i1,1:grid%nMu) + t1 * &
                  (grid%muProbDistLine(i1+1,1:grid%nMu) - grid%muProbDistLine(i1,1:grid%nMu))
             call locate(tempmuProbDistLine, grid%nMu, r1, i2)
             t2 = (r1-tempmuProbDistLine(i2))/(tempmuProbDistLine(i2+1)-tempmuProbDistLine(i2))
             mu = grid%muAxis(i2) + t2 * (grid%muAxis(i2+1)-grid%muAxis(i2))


             call random_number(r1)
             tempphiProbDistLine(1:grid%nphi) = &
                  (1.d0-t1)*(1.d0-t2) * grid%phiProbDistLine(i1  , i2  , 1:grid%nphi) +&
                  (   t1)*(1.d0-t2) * grid%phiProbDistLine(i1+1, i2  , 1:grid%nphi) +&
                  (1.d0-t1)*(   t2) * grid%phiProbDistLine(i1  , i2+1, 1:grid%nphi) +&
                  (     t1)*(   t2) * grid%phiProbDistLine(i1+1, i2+1, 1:grid%nphi)
             tempphiProbDistLine(1:grid%nphi) = tempphiprobDistLine(1:grid%nphi) / tempphiprobDistLine(grid%nPhi)
             call locate(tempphiProbDistLine(1:grid%nPhi), grid%nPhi, r1, i3)
             t3 = (r1-tempphiProbDistLine(i3)) / & 
                  (tempphiProbDistLine(i3+1)-tempphiProbDistLine(i3))
             phi = grid%phiAxis(i3) + t3 * (grid%phiAxis(i3+1)-grid%phiAxis(i3))

                if (.not.grid%inUse(i1,i2,i3)) then
                   goto 990
                endif
             ! set up the position

          endif

             sinTheta = sqrt(1.d0-mu**2)
             thisPhoton%position%x = r*cos(phi)*sinTheta
             thisPhoton%position%y = r*sin(phi)*sinTheta
             thisPhoton%position%z = r*mu
                if (useBias.and.(.not.grid%cartesian)) then
                   biasWeight = interpgridscalar2(grid%biasline3d, grid%nr,grid%nmu,grid%nphi,i1,i2,i3,real(t1),real(t2),real(t3))
                   biasWeight = 1.d0/biasWeight
                endif
             endif

       endif

       ! interpolate to get the velocity

       if (.not.thisPhoton%resonanceLine) then
               
          if (grid%adaptive) then
             octalPoint = thisPhoton%position  
             
             !  Assiging the velocity at the emission location
!             if (VoigtProf) then
                CALL amrGridValues(grid%octreeRoot,octalPoint,  &
                     velocity=velocity, temperature=temperature, &
                     Ne=Ne, rho=rho)
                ! Setting the velocity at the emission location + offset 
                !  velocity (rVel) by the thermal motion of gas.
                rVel = thermalHydrogenVelocity(temperature)  ! [C] (vector)
                vTherm = sqrt(2.* kErg * temperature/mHydrogen)/cSpeed
!                if (thisPhoton%linePhoton) then
!
!                   do i = 1, nv
!                      vArray(i) = real(i-1)/real(nv-1) * 5. * vTherm
!                      if (vArray(i) /= 0.) then
!                         vbias(i) = 1./fvmaxwellian(vArray(i)*cSpeed, mHydrogen, temperature)
!                      else
!                         vbias(i) = 1.
!                      endif
!                   enddo
!
!                   fac = 0.
!                   pv(1) = 0.
!                   do i  = 2, nv
!                      dv(i) = vArray(i)-vArray(i-1)
!                      pv(i) = fvmaxwellian(vArray(i)*cSpeed, mHydrogen, temperature)
!                      fac = fac + pv(i) * dv(i)
!                      pv(i) = pv(i-1) + pv(i) * vbias(i) * dv(i)
!                   enddo
!                   vBias(1:nv) = vBias(1:nv) * fac / pv(nv)
!                   pv(1:nv) = pv(1:nv) / pv(nv)
!
!                   call random_number(r)
!                   call locate(pv, nv, r, i)
!                   
!
!                   rVel = vArray(i) *  randomUnitVector()
!                   thisPhoton%stokes = thisPhoton%stokes  * (1./vbias(i))
!                endif
                dopShift = modulus(rVel)/vTherm - 1.
                thisPhoton%velocity = velocity + rVel
                
                ! We shuffle the emission frequency acoording to the shape 
                ! of the Voigt profile.
                N_HI = MAX((rho/mHydrogen - Ne), 1.d-25) ! number density of HI.
                nu = cSpeed_dbl/dble(lamline*angstromtocm)  ! [Hz]
                if (VoigtProf) then
                   Gamma = bigGamma(dble(N_HI), dble(temperature), dble(Ne), nu)
                else
                   Gamma = 0.0d0
                end if
                nu_shuffled = random_Lorentzian_frequency(nu, Gamma) ! [Hz]                
                lambda_shuffled = (cSpeed_dbl/nu_shuffled)*1.e8  ! [A]
!                lambda_shuffled = lamline  ! [A]
                ! ==> This will be assigend to the photon later.
!             else
!                thisPhoton%velocity = amrGridVelocity(grid%octreeRoot,octalPoint)!,&
!                !startOctal=sourceOctal,actualSubcell=subcell)
!                rVel = ((gasdev()*20.e5)/cSpeed)*randomUnitVector()
!                thisPhoton%velocity = thisPhoton%velocity + rVel
!             end if
             
          
          else
             thisPhoton%velocity = &
                  ((1.d0-t1)*(1.d0-t2)*(1.d0-t3))*grid%velocity(i1,i2,i3) + &
                  ((1.d0-t1)*(1.d0-t2)*t3     )*grid%velocity(i1,i2,i3+1) + &
                  ((1.d0-t1)*t2     *t3     )*grid%velocity(i1,i2+1,i3+1) + &
                  (t1     *(1.d0-t2)*(1.d0-t3))*grid%velocity(i1+1,i2,i3) + &
                  ((1.d0-t1)*t2     *(1.d0-t3))*grid%velocity(i1,i2+1,i3) + &
                  (t1     *(1.d0-t2)* t3    )*grid%velocity(i1+1,i2,  i3+1) + &
                  (t1     *t2     *(1.d0-t3))*grid%velocity(i1+1,i2+1,i3  ) + &
                  (t1     *t2     *t3     )*grid%velocity(i1+1,i2+1,i3+1)
          end if
       endif

       if (grid%geometry == "donati") then
          if (grid%adaptive) then
             thisPhoton%velocity = thisPhoton%velocity + &
                                   thermalHydrogenVelocity(amrGridTemperature(grid%octreeRoot,octalPoint))
          else
             thisPhoton%velocity = thisPhoton%velocity + thermalHydrogenVelocity(grid%temperature(i1,i2,i3))
          end if
       endif


       if (secondSource) thisPhoton%velocity = ramanSourceVelocity

       ! Assigining the photon wavelength 
       if (.not.thisPhoton%resonanceLine) then
!          if (StarkBroadening) then 
             ilambda = int(real(nLambda)* &
                  (lambda_shuffled-lambda(1))/(lambda(nLambda)-lambda(1)))+1
             thisPhoton%lambda = lambda_shuffled
!          else
!             ilambda = int(real(nLambda)* &
!                  (lamLine-lambda(1))/(lambda(nLambda)-lambda(1)))+1             
!             thisPhoton%lambda = lamLine
!          end if
       endif

       thisPhoton%stokes = thisPhoton%stokes * weightLinePhoton


       ! random direction
       thisPhoton%direction = randomUnitVector()

       if (grid%doRaman) then
          thisPhoton%velocity = &
               maxwellianVelocity(16.d0*mHydrogen, grid%tempSource)/real(cSpeed)
       endif


    endif ! (thisPhoton%linePhoton) 



    ! set up the two normals

    thisPhoton%normal = VECTOR(0., 0., 0.)
    do while(modulus(thisPhoton%normal) == 0.)
       thisPhoton%normal = randomUnitVector() .cross. thisPhoton%direction
    enddo
    call normalize(thisPhoton%normal)


    thisPhoton%oldnormal = VECTOR(0., 0., 0.)
    do while(modulus(thisPhoton%oldnormal) == 0.)
       thisPhoton%oldnormal = randomUnitVector() .cross. thisPhoton%direction
    enddo
    call normalize(thisPhoton%oldnormal)


    thisPhoton%stokes = thisPhoton%stokes * biasWeight


  end subroutine initPhoton





  ! this subroutine finds the grid array indices of a particular
  ! photon - OUTDATED


  subroutine findPosition(thisPhoton, grid, i1,i2,i3)

    type(PHOTON) :: thisPhoton
    type(GRIDTYPE) :: grid
    real :: r, mu, phi
    integer :: i1,i2,i3
    if (grid%cartesian) then
       
       ! cartesian case is easy

       call locate(grid%xAxis, grid%nx, real(thisPhoton%position%x), i1)
       call locate(grid%yAxis, grid%ny, real(thisPhoton%position%y), i2)
       call locate(grid%zAxis, grid%nz, real(thisPhoton%position%z), i3)

    else

       ! polar case assumes evenly space mu and phi axes

       r = modulus(thisPhoton%position)
       mu = thisPhoton%position%z / r
       if ((thisPhoton%position%x == 0.) .and. (thisPhoton%position%y == 0.)) then
          phi = 0.
       else
          phi = atan2(thisPhoton%position%y, thisPhoton%position%x)
       endif
       if (phi < 0.) phi = phi + twoPi
       call locate(grid%rAxis, grid%nr, r, i1)
       i2 = int(real(grid%nMu)*0.5*(mu+1.d0))+1
       i3 = int(real(grid%nphi)*phi/twoPi)+1
    endif

  end subroutine findPosition


  subroutine getDiskVelocity(vel)
    type(VECTOR) :: vel, zAxis, rHat, vHat
    real :: mHot, ang, radius, r, v, rInner, vMax

    zAxis = VECTOR(0.,0.,1.)
    mHot = 0.6 * mSol
    vMax = 30.e5

    rInner = (bigG *  mHot)/(vMax**2)
    call random_number(r)
    radius = rInner + rSol * sqrt(r)

    v = sqrt(bigG * mHot/radius)


    call random_number(ang)
    ang = 2.*pi*ang

    rHat = VECTOR(cos(ang),sin(ang),0.)
    vHat = rHat .cross. zAxis
    call normalize(vHat)

    call random_number(r)
    if (r < 0.5) then
       vHat = VECTOR(1.,0.,0.)
    else
       vHat = VECTOR(-1.,0.,0.)
    endif

    v = 30.e5

    vel = (v/cSpeed) * vHat
!    vel = VECTOR(0.,0.,0.)

  end subroutine getDiskVelocity

  subroutine initPlanetPhoton(thisPhoton, grid, lamLine)
    type(GRIDTYPE) :: grid
    type(PHOTON) :: thisPhoton
    real ::  r
    type(VECTOR) :: direction, position, toPlanet
    real :: cosTheta
    real :: lamLine, x, y, d, x1, x2
    logical :: ok
    

    r = 2.*grid%rAxis(grid%nr)
    do while ( r > grid%rAxis(grid%nr))
       call random_number(x)
       x = (2.*x-1.) * grid%rAxis(grid%nr)
       call random_number(y)
       y = (2.*y-1.) * grid%rAxis(grid%nr)
       r = sqrt(x**2 + y**2)
    enddo
    position = VECTOR(grid%rAxis(grid%nr), x, y)
    direction = VECTOR(-1.,0.,0.)
    d = modulus(position)
    toPlanet = (-1./d)*position
    cosTheta = toPlanet .dot. direction
    call solveQuad(1.,-2.*d*cosTheta,d*d-grid%rAxis(grid%nr)*grid%rAxis(grid%nr),x1,x2,ok)
    x = min(x1,x2)
    position = position + x*direction

    thisPhoton%position = position
    thisPhoton%direction = direction

    thisPhoton%weight = 1.
    thisPhoton%stokes%i = 1.

    thisPhoton%stokes%q = 0.
    thisPhoton%stokes%u = 0.
    thisPhoton%stokes%v = 0.

    thisPhoton%velocity = VECTOR(0.,0.,0.)

    thisPhoton%contPhoton = .true.
    thisPhoton%lambda = lamLine

    ! set up the two normals

    thisPhoton%normal = randomUnitVector() .cross. thisPhoton%direction
    call normalize(thisPhoton%normal)



    thisPhoton%oldnormal = randomUnitVector() .cross. thisPhoton%direction
    call normalize(thisPhoton%oldnormal)


    thisPhoton%linePhoton = .false.
    thisPhoton%contPhoton = .true.

  end subroutine initPlanetPhoton

end module photon_mod

   
