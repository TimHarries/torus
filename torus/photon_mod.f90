!
! this module defines and describes the photon type
! and includes subroutines for photon initialization
! and scattering
!

! written by tjh

! v1.0 on 16/08/99

module photon_mod

  use vector_mod            ! the vector maths
  use constants_mod         ! physical constants
  use phasematrix_mod       ! phase matrix manipulation
  use grid_mod              ! the opacity grid
  use math_mod              ! maths utils
  use utils_mod
  implicit none

  public

  ! the description of the photon type


  type PHOTON
     type(STOKESVECTOR) :: stokes         ! the stokes intensities
     real :: weight                       ! photon weight
     real :: lambda                       ! wavelength
     type(VECTOR) :: normal               ! scattering normal
     type(VECTOR) :: oldNormal            ! the last scattering normal
     type(VECTOR) :: position             ! the photon position
     type(VECTOR) :: direction            ! the photon direction
     type(VECTOR) :: velocity             ! the photon 'velocity'
     type(VECTOR) :: originalNormal
     logical :: contPhoton                ! continuum photon? or
     logical :: linePhoton                ! line photon?
     logical :: redPhoton                 ! raman photon in red region
     real :: redLambda                    ! raman wavelength in red
     logical :: fromStar1
     logical :: fromStar2
     logical :: resonanceLine             ! resonance line photon?
  end type PHOTON

contains

  ! this subroutine rotates the stokes vectors of a photon - this
  ! subroutine should probably moved from here to just work
  ! on a stokes vector (see phasematrix_mod)

  subroutine rotate(thisPhoton, cosx, sinx)

    type (PHOTON) :: thisPhoton          ! the photon
    real :: qdash, udash                 ! the rotated stokes vector
    real :: cosx, sinx                   ! x should be twice theta

    ! perform the rotation

    if ( (sinx**2 + cosx**2) < 0.999) then
       write(*,*) "! non-normalized rotation attempt",sinx**2+cosx**2
       stop
    endif

    qdash = cosx * thisPhoton%stokes%q + sinx * thisPhoton%stokes%u
    udash =-sinx * thisPhoton%stokes%q + cosx * thisPhoton%stokes%u

    thisPhoton%stokes%q = qdash
    thisPhoton%stokes%u = udash

  end subroutine rotate


  ! this subroutine performs a scattering

  subroutine scatterPhoton(grid, thisPhoton, givenVec, outPhoton, mie, &
        miePhase, nLambda, nMuMie, lamStart, lamEnd)

    type(GRIDTYPE) :: grid                       ! the opacity grid
    type(PHOTON) :: thisPhoton, outPhoton        ! current/output photon
    type(VECTOR) :: incoming, outgoing           ! directions
    type(VECTOR) :: obsNormal, refNormal         ! scattering normals
    type(VECTOR) :: sVec, givenVec      ! vectors
    type(VECTOR) :: zAxis                        ! the z-axis
    type(PHASEMATRIX) :: rayleighPhase           ! rayleigh phase matrix
    logical :: mie                               ! is this a mie scattering?
    integer :: i, j                              ! counters
    integer :: nLambda                           ! size of wavelength array
    integer :: nMumie                            ! number of mu angles for mie
    type(PHASEMATRIX) :: miePhase(nLambda, nMumie) ! mie phase matrices
    real :: lamStart, lamEnd                     ! start and end wavelengths
    real :: costheta                             ! cos scattering angle
    real :: ang                                  ! scattering angle
    real :: r1, r2, r                            ! radii
    integer :: i1, i2, i3                        ! position indices
    real :: t1,t2,t3                             ! interp factors
    real :: u,v,w,t                              ! direction components
    real :: sinGamma, cosGamma, sin2Gamma, cos2Gamma
    logical :: randomDirection                   ! is this a random direction
    real :: vRay, vOverCsqr
    real :: fac
    real :: nuRest
    real :: theta, phi, mu
!    real :: dx

!    dx = abs(grid%xAxis(2) - grid%xAxis(1))

    ! initialize variables

    outPhoton = thisPhoton                    
    incoming = outPhoton%direction

    zAxis = VECTOR(0.,0.,1.)

    outgoing = givenVec

    randomDirection = .false.

    ! if the outgoing vector is the zero vector then this flags that
    ! we are going to scattering into a random direction

    if (modulus(outgoing) == 0.) then
       randomDirection = .true.
       call random_number(r1)
       w = 2.d0*r1 - 1.d0
       t = sqrt(1.d0-w*w)
       call random_number(r2)
       ang = Pi * (2.*r2-1.d0)
       u = t*cos(ang)
       v = t*sin(ang)
       outgoing%x = u
       outgoing%y = v
       outgoing%z = w
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
       if (.not.mie) then
          
          ! set up the rayleigh phase matrix and apply it

          rayleighPhase = fillRayleigh(costheta)
          outPhoton%stokes = apply(rayleighPhase, outPhoton%stokes)
       else
          
          ! do the same for a mie scattering
          
          i = int(real(nLambda)*(thisPhoton%lambda-lamStart)/(lamEnd-lamStart))+1
          j = int(0.5*(costheta+1.d0)*real(nmumie))+1
          outPhoton%stokes = apply(miePhase(i,j), outPhoton%stokes)
          
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
       call findPosition(outPhoton, grid, i1, i2, i3)
       if (grid%cartesian) then
          t1 = (outPhoton%position%x - grid%xAxis(i1))/(grid%xAxis(i1+1)-grid%xAxis(i1))
          t2 = (outPhoton%position%y - grid%yAxis(i2))/(grid%yAxis(i2+1)-grid%yAxis(i2))
          t3 = (outPhoton%position%z - grid%zAxis(i3))/(grid%zAxis(i3+1)-grid%zAxis(i3))
       else
          call getPolar(outPhoton%position, r, theta, phi)
          mu = cos(theta)
          t1 = (r - grid%rAxis(i1))/(grid%rAxis(i1+1)-grid%rAxis(i1))
          t2 = (mu - grid%muAxis(i2))/(grid%muAxis(i2+1)-grid%muAxis(i2))
          t3 = (phi - grid%phiAxis(i3))/(grid%phiAxis(i3+1)-grid%phiAxis(i3))
       endif

       if (.not.grid%doRaman) then
          if (.not.grid%resonanceLine) then
             outPhoton%velocity = interpGridVelocity(grid,i1,i2,i3,t1,t2,t3) + &
                  thermalElectronVelocity(grid%temperature(i1,i2,i3))
          else
             outPhoton%velocity = interpGridVelocity(grid,i1,i2,i3,t1,t2,t3)
          endif
       else
          outPhoton%velocity = interpGridVelocity(grid,i1,i2,i3,t1,t2,t3) !+ &
!               thermalHydrogenVelocity(grid%temperature(i1,i2,i3))
       endif

       vray = (outPhoton%velocity-thisPhoton%velocity) .dot. incoming
       vovercsqr = (outPhoton%velocity-thisPhoton%velocity) .dot. &
            (outPhoton%velocity-thisPhoton%velocity)
       fac = (1.d0 - 0.5d0*vovercsqr*(1.d0-0.25d0*vovercsqr))/(1.d0 + vray)
       outPhoton%lambda = outPhoton%lambda  / fac

    else

       outPhoton%redPhoton = .true.
       call findPosition(outPhoton, grid, i1, i2, i3)
       t1 = (outPhoton%position%x - grid%xAxis(i1))/(grid%xAxis(i1+1)-grid%xAxis(i1))
       t2 = (outPhoton%position%y - grid%yAxis(i2))/(grid%yAxis(i2+1)-grid%yAxis(i2))
       t3 = (outPhoton%position%z - grid%zAxis(i3))/(grid%zAxis(i3+1)-grid%zAxis(i3))
       outPhoton%velocity = interpGridVelocity(grid,i1,i2,i3,t1,t2,t3)

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
       pencilBeam, secondSource, secondSourcePosition, lumRatio, &
       ramanSourceVelocity, vo6, contWindPhoton, directionalWeight,useBias,theta1,theta2, chanceHotRing)

    type(PHOTON) :: thisPhoton                 ! the photon
    type(GRIDTYPE) :: grid                     ! the opacity grid
    integer :: nLambda, iLambda                ! wavelength indices
    logical :: useBias
    real :: weightContPhoton, weightLinePhoton ! the photon weights
    logical :: contWindPhoton                  ! is this continuum photon produced in the wind
    logical :: ok
    real :: lumRatio
    real :: vo6
    real :: x,y,z
    real :: directionalWeight
    real :: lambda(:), sourceSpectrum(*)       ! wavelength array/spectrum
    real :: r1,r2,r3                           ! radii
    real :: u, v, w, t                         ! direction components
    real :: r, mu, phi                         ! spherical polar coords
    real :: sinTheta      
    real :: t1, t2, t3                         ! multipliers
    real :: vPhi                               ! azimuthal velocities
    real :: ang                                ! angle
    logical :: contPhoton                      ! is this a continuum photon?
    logical :: lineEmission                    ! is there line emission?
    real :: lamLine                            ! wavelength of the line
    logical :: flatspec                        ! is the spectrum flat
    real :: vRot                               ! rotational velocity
    integer :: i1, i2, i3                      ! position indices
    integer :: j, k
    type(VECTOR) :: rHat, perp                 ! radial unit vector
    type(VECTOR) :: rotatedVec                 ! rotated vector
    logical :: pencilBeam                      ! beamed radiation
    type(VECTOR), parameter :: zAxis = VECTOR(0.,0.,1.) ! the z axis
    logical :: secondSource                    ! second photon source?
    type(VECTOR) :: secondSourcePosition       ! the position of it
    type(VECTOR) :: ramanSourceVelocity        ! what it says
    real :: biasWeight
    integer :: i

    real :: theta1, theta2                     ! defines hot ring of accretion for TTaus
    real :: chanceHotRing                      ! chance of core photon in ttauri accretion ring
    real :: thisTheta, thisPhi
    
    real :: tempXProbDist(2000)
    real :: tempYProbDist(2000)
    real :: tempZProbDist(2000)
    real :: temprProbDistLine(2000)
    real :: tempmuProbDistLine(2000)
    real :: tempphiProbDistLine(2000)

    ! set up the weights and the stokes intensities (zero at emission)

    thisPhoton%resonanceLine = .false.

    if (grid%resonanceLine) thisPhoton%resonanceLine = .true.

    vo6 = 0.
    ramanSourceVelocity = VECTOR(0.,0.,0.)

    directionalWeight = 1.
    biasWeight = 1.
    thisPhoton%fromstar1 = .false.
    thisPhoton%fromstar2 = .false.

    thisPhoton%weight = 1.d0
    thisPhoton%stokes%i = 1.d0

    thisPhoton%stokes%q = 0.
    thisPhoton%stokes%u = 0.
    thisPhoton%stokes%v = 0.



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


    ! if it is a continuum photon initialize its position

    if (thisPhoton%contPhoton) then




       if (grid%cartesian.and.(.not.grid%lineEmission)) then

          ! cartesian photons are produced at the grid centre, or at
          ! the second source position

          thisPhoton%position = VECTOR(0.,0.,0.)


       else

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
                        (1.d0/ interpGridScalar2(grid%biasCont3d,grid%nx,grid%ny,grid%nz,i1,i2,i3,t1,t2,t3))

             else


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

                x = r * sqrt(1.d0-mu*mu) * cos(phi)
                y = r * sqrt(1.d0-mu*mu) * sin(phi)
                z = r * mu
                if (useBias.and.(.not.grid%cartesian)) then
                   biasWeight = grid%biasCont(i1) + t1 * (grid%biasCont(i1+1)-grid%biasCont(i1))
                   biasWeight = 1.d0/biasWeight
                endif

             endif
             thisPhoton%position = VECTOR(x, y, z)


          else
             if (grid%geometry /= "binary") then
                r = grid%rCore*1.00001d0
		if (grid%geometry /= "ttauri") then
		   thisPhoton%position = (r*randomUnitVector())
		else
		   call random_number(r1)
		   if (r1 < chanceHotRing) then
		      call random_number(r1)
		      thisTheta = r1*(theta2 - theta1)+theta1
		      call random_number(r1)
		      if (r1 < 0.5) thisTheta = pi - thisTheta
		      call random_number(r1)
		      thisPhi = twoPi * r1
		      rHat = VECTOR(cos(thisPhi)*sin(thisTheta),sin(thisPhi)*sin(thisTheta),cos(thisTheta))
		      thisPhoton%position = r * rHat
		   else
		      thisPhoton%position = (r*randomUnitVector())
		   endif
		endif

                if (r /= 0.) then		   
                   thisPhoton%originalNormal = thisPhoton%position
                else
                   thisPhoton%originalNormal = randomUnitVector()
                endif
                call normalize(thisPhoton%originalNormal)
                thisPhoton%velocity = VECTOR(1.e-30,1.e-30,1.e-30)
                

             else
                call random_number(r)
                if (r < grid%lumRatio) then
                   thisPhoton%position = grid%starPos1 + (1.01*(grid%rStar1) * randomUnitVector())
                   thisPhoton%fromStar1 = .true.
                   call getIndices(grid,thisPhoton%position,i1,i2,i3,t1,t2,t3)
                else
                   thisPhoton%position = grid%starPos2 + (1.01*(grid%rStar2) * randomUnitVector())
                   thisPhoton%fromStar2 = .true.
                   call getIndices(grid,thisPhoton%position,i1,i2,i3,t1,t2,t3)
                endif
             endif

          endif

       endif

       ! get wavelength

       call random_number(r1)
       iLambda = int(r1*real(nLambda))+1
       thisPhoton%lambda = lambda(iLambda)
       if (.not.flatspec) then
          thisPhoton%stokes = thisPhoton%stokes * &
               (sourceSpectrum(iLambda) * weightContPhoton)
       else
          thisPhoton%stokes = thisPhoton%stokes * weightContPhoton
          thisPhoton%lambda = lamLine
       endif


       ! get direction

       if (.not.pencilBeam) then
          r3 = modulus(thisPhoton%position)
          if (r3 /= 0.) then
             if (.not.grid%geometry == "binary") then
                rHat = thisPhoton%position / r3
             else
                if (thisPhoton%fromStar1) then
                   rHat = (thisPhoton%position-grid%starPos1) 
                   call normalize(rHat)
                   thisPhoton%originalNormal = rHat
                else
                   rHat = (thisPhoton%position-grid%starPos2) 
                   call normalize(rHat)
                   thisPhoton%originalNormal = rHat
                endif
             endif
          endif
          thisPhoton%direction = randomUnitVector()
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
          rotatedVec = rotateY(thisPhoton%direction, 30.*degToRad)
          thisPhoton%direction = rotatedVec
       endif



       if (.not.contWindPhoton) then
          if (.not.grid%geometry=="rolf") then
             t = (thisPhoton%direction .dot. rHat)

             ! must be outwards from the photosphere if this is a core continuum photon

             if ((t < 0.).and.(r3 /= 0.).and.(grid%lineEmission).and.(.not.contWindPhoton)) then
                thisPhoton%direction = (-1.)*thisPhoton%direction
             endif
             t = (thisPhoton%direction .dot. rHat)
             directionalWeight = abs(2.*t)
          endif



       endif

       ! this does the photon velocity for a rotational velocity field
       ! not the local grid velocity, since the radial component must
       ! be zero

       if (grid%geometry == "binary") then
          call getIndices(grid, thisPhoton%position, i1, i2, i3, t1, t2, t3)
          thisPhoton%velocity = grid%velocity(i1,i2,i3)
       endif


       if (.not.grid%cartesian) then
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


    endif


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
               (1.d0/ interpGridScalar2(grid%biasLine3D,grid%nx,grid%ny,grid%nz,i1,i2,i3,t1,t2,t3))

       else


          if (thisPhoton%resonanceLine) then
             r = grid%rCore*1.00001d0
             rHat = randomUnitVector()
             thisPhoton%position = r * rHat
             thisPhoton%originalNormal = thisPhoton%position
             thisPhoton%direction = randomUnitVector()
             t = (thisPhoton%direction .dot. rHat)

             ! must be outwards from the photosphere if this is a core continuum photon

             if (t < 0.) then
                thisPhoton%direction = (-1.)*thisPhoton%direction
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



             call random_number(r1)
             temprProbDistLine(1:grid%nr) = grid%rProbDistLine(1:grid%nr)
             call locate(temprProbDistLine, grid%nr, r1, i1)
             t1 = (r1-temprProbDistLine(i1))/(temprProbDistLine(i1+1)-temprProbDistLine(i1))
             r = grid%rAxis(i1) + t1 * (grid%rAxis(i1+1)-grid%rAxis(i1))


             call random_number(r1)
             tempMuProbDistLine = grid%muProbDistLine(i1,1:grid%nMu) + t1 * &
                  (grid%muProbDistLine(i1+1,1:grid%nMu) - grid%muProbDistLine(i1,1:grid%nMu))
             call locate(tempmuProbDistLine, grid%nMu, r1, i2)
             t2 = (r1-tempmuProbDistLine(i2))/(tempmuProbDistLine(i2+1)-tempmuProbDistLine(i2))
             mu = grid%muAxis(i2) + t2 * (grid%muAxis(i2+1)-grid%muAxis(i2))


             call random_number(r1)
             tempphiProbDistLine(1:grid%nphi) = &
                  (1.d0-t1)*(1.d0-t2) * grid%phiProbDistLine(i1  , i2  , 1:grid%nphi) +&
                  (   t1)*(1.d0-t2) * grid%phiProbDistLine(i1+1, i2  , 1:grid%nphi) +&
                  (1.d0-t1)*(   t2) * grid%phiProbDistLine(i1  , i2+1, 1:grid%nz) +&
                  (     t1)*(   t2) * grid%phiProbDistLine(i1+1, i2+1, 1:grid%nphi)
             tempphiProbDistLine(1:grid%nphi) = tempphiprobDistLine(1:grid%nphi) / tempphiprobDistLine(grid%nPhi)
             call locate(tempphiProbDistLine(1:grid%nPhi), grid%nPhi, r1, i3)
             t3 = (r1-tempphiProbDistLine(i3)) / & 
                  (tempphiProbDistLine(i3+1)-tempphiProbDistLine(i3))
             phi = grid%phiAxis(i3) + t3 * (grid%phiAxis(i3+1)-grid%phiAxis(i3))
             
             ! set up the position

             sinTheta = sqrt(1.d0-mu**2)
             thisPhoton%position%x = r*cos(phi)*sinTheta
             thisPhoton%position%y = r*sin(phi)*sinTheta
             thisPhoton%position%z = r*mu
             if (useBias.and.(.not.grid%cartesian)) then
                biasWeight = log10(grid%biasLine(i1)) + t1 * (log10(grid%biasLine(i1+1))-log10(grid%biasLine(i1)))
                biasWeight = 10.**biasWeight
                biasWeight = 1.d0/biasWeight
             endif
          endif

       endif

       ! interpolate to get the velocity

       if (.not.thisPhoton%resonanceLine) then
          thisPhoton%velocity = &
               ((1.d0-t1)*(1.d0-t2)*(1.d0-t3))*grid%velocity(i1,i2,i3) + &
               ((1.d0-t1)*(1.d0-t2)*t3     )*grid%velocity(i1,i2,i3+1) + &
               ((1.d0-t1)*t2     *t3     )*grid%velocity(i1,i2+1,i3+1) + &
               (t1     *(1.d0-t2)*(1.d0-t3))*grid%velocity(i1+1,i2,i3) + &
               ((1.d0-t1)*t2     *(1.d0-t3))*grid%velocity(i1,i2+1,i3) + &
               (t1     *(1.d0-t2)* t3    )*grid%velocity(i1+1,i2,  i3+1) + &
               (t1     *t2     *(1.d0-t3))*grid%velocity(i1+1,i2+1,i3  ) + &
               (t1     *t2     *t3     )*grid%velocity(i1+1,i2+1,i3+1)
       endif

       if (grid%geometry == "donati") then
          thisPhoton%velocity = thisPhoton%velocity + thermalHydrogenVelocity(grid%temperature(i1,i2,i3))
       endif


       if (secondSource) thisPhoton%velocity = ramanSourceVelocity


       if (.not.thisPhoton%resonanceLine) then
          ilambda = int(real(nLambda)* &
               (lamLine-lambda(1))/(lambda(nLambda)-lambda(1)))+1
          thisPhoton%lambda = lamLine
       endif

       thisPhoton%stokes = thisPhoton%stokes * weightLinePhoton



    endif


    if (thisPhoton%linePhoton) then

       ! random direction

       thisPhoton%direction = randomUnitVector()

       if (grid%doRaman) then
          thisPhoton%velocity = &
               maxwellianVelocity(16.*mHydrogen, grid%tempSource)/cSpeed
       endif


    endif


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

       call locate(grid%xAxis, grid%nx, thisPhoton%position%x, i1)
       call locate(grid%yAxis, grid%ny, thisPhoton%position%y, i2)
       call locate(grid%zAxis, grid%nz, thisPhoton%position%z, i3)

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

   
