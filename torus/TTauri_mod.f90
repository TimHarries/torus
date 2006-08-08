! this module contains stuff about filling the grid with T Tauri style
! models.

! v1.0 8 March 2001   @Exeter!

module ttauri_mod

  use kind_mod
  use constants_mod                   ! physical constants
  use vector_mod                      ! vector math
  use atom_mod                        ! LTE atomic physics
  use utils_mod
  use gridtype_mod
  use grid_mod
  use stateq_mod
  use math_mod
  use octal_mod
  use parameters_mod
  use amr_mod

  implicit none
  
  public :: fillGridMagneticAccretion, &
       &    infallEnhancment, &
       &    initInfallEnhancement, &
       &    fillGridFlaredDisk


contains
  
  subroutine fillGridMagneticAccretion(grid,contfile1, popFileName, &
         readPops, writePops, lte, lamLine, Laccretion, Taccretion, &
         sAccretion, curtains, dipoleOffset, nLower, nUpper, theta1, theta2)

    type(GRIDTYPE) :: grid
    character(len=*) :: contFile1, popFileName
    character(len=80) :: contfile2, newContfile
    logical :: readPops, writePops, lte
    logical :: curtains
    real :: dipoleOffset
    type(VECTOR) :: rVec
    integer, parameter :: nPhi = 360
    integer, parameter :: nTheta = 201
    integer, parameter :: nR = 100
    integer :: i, j, k
    integer :: i1, i2, i3
    integer :: nLower, nUpper
    real :: ang
    real :: t1,t2,t3, lamLine, nu
    real :: rInner, rOuter
    real :: theta1, theta2
    real :: thetaStart, thetaEnd
    real :: mDot, rStar, mStar
    real :: r, theta, phi
    real :: rM, y, rho
    real :: nuArray(3000),fnu(3000), fac
    integer :: nNu
    real :: tot
    type(VECTOR) :: vP, posVec
    real :: modVp
    real(double) :: Laccretion
    real :: Taccretion, sAccretion
    real :: minRho, maxRho
    !type(VECTOR) :: rHat, vHat
    real :: mLoss, vTerm, v0 ! vel, dist, rMag
    integer, parameter :: nCalvet = 22
    !integer :: iCalvet
    !real :: tCalvet(nCalvet) = (/10000.,10000.,10000.,10000.,10000.,10000.,10000., &
    !                   9600., 9000., 8000.,7300.,6600., 6000., 5400., &
    !                   5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000./)
    !real :: rCalvet(nCalvet) = (/0.190, 0.270, 0.310, 0.360, 0.400, 0.450, &
    !                             0.569, 0.742, 0.938, 1.12, 1.31, 1.91, &
    !                             2.50, 3.10, 3.50, 3.90, 4.88, 5.89, 6.89, &
    !                             8.09, 9.04, 10.48 /)


!    if (.not.grid%cartesian) then
!       write(*,'(a)') "! this geometry is only suitable for cartesians"
!       stop
!    endif

    grid%geometry = "ttauri"
    mDot = 1.e-7 * mSol / (365.25 * 24. * 3600.)

    mLoss = 1.e-8 * mSol / (365.25 * 24. * 3600.)
    vTerm = 200. * 1.e5
    v0 = 10. * 1e5
    rStar = 2.*rSol
    rInner = 2.2 * rStar
    rOuter = 3. * rStar
    mStar = 0.8 * mSol
    grid%rCore = rStar /1.e10

    grid%dipoleOffset = dipoleOffset
    grid%etaCont = 1.e-20
    grid%etaLine = 1.e-20
    grid%chiLine = 1.e-30
    grid%kappaSca = 1.e-30
    grid%kappaAbs = 1.e-30

    grid%diskRadius = rInner/1.e10

    grid%diskNormal = VECTOR(0.,0.,1.)
    grid%diskNormal = rotateX(grid%diskNormal, dipoleOffSet)

    grid%inUse = .false.

    theta1 = asin(sqrt(rStar/rOuter))
    theta2 = asin(sqrt(rStar/rInner))



    thetaStart = theta1
    thetaEnd = pi - theta1


    grid%rho = 1.e-25
    grid%velocity = VECTOR(1.e-20,1.e-20,1.e-20)

    grid%starPos1 = VECTOR(0.,0.,0.)
    grid%rStar1 = rStar/1.e10
    grid%lineEmission = .true.

    grid%biasLine3D = 1.
    grid%biasCont3D = 1.

    minRho = 1.e20
    maxRho =-1.e20

    if (grid%cartesian) then
       do i = 1, grid%nx
          grid%xAxis(i) = 2.*(real(i-1)/real(grid%nx-1)) - 1.
       enddo
       do i = 1, grid%ny
          grid%yAxis(i) = 2.*(real(i-1)/real(grid%ny-1)) - 1.
       enddo
       do i = 1, grid%nz
          grid%zAxis(i) = 2.*(real(i-1)/real(grid%nz-1)) - 1.
       enddo
       grid%xAxis(1:grid%nx) = grid%xAxis(1:grid%nx) * rOuter * 1.05 / 1.e10
       grid%yAxis(1:grid%ny) = grid%yAxis(1:grid%ny) * rOuter * 1.05 / 1.e10
       grid%zAxis(1:grid%nz) = grid%zAxis(1:grid%nz) * rOuter * 1.05 / 1.e10
    else
       do i = 1, grid%nr
          grid%rAxis(i) = grid%rCore + ((rOuter/1.e10)-grid%rCore)*real(i-1)/real(grid%nr-1)
       enddo
       do i = 1, grid%nMu
          grid%muAxis(i) = 2.*(real(i-1)/real(grid%nMu-1)) - 1.
       enddo
       do i = 1, grid%nPhi
          grid%phiAxis(i) = twoPi * real(i-1)/real(grid%nPhi-1)
       enddo
    endif

    grid%inStar = .false.

    if (grid%cartesian) then
       grid%instar = .false.
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                if (modulus(VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))) < grid%rCore) grid%inStar(i,j,k) = .true.
             enddo
          enddo
       enddo
    endif

    write(*,'(a)') "Filling grid with magnetic accretion model..."

    do i = 1, nPhi
       phi = twoPi * real(i-1)/real(nPhi-1)
       do j = 1, nR
          fac = real(j-1)/real(nR-1)
          rM = rInner + (rOuter - rInner)*fac
          

          do k = 1, nTheta
             theta = thetaStart + (thetaEnd-thetaStart)*real(k-1)/real(nTheta-1)

             r = rM * sin(theta)**2

             y = r/rM

             vP = VECTOR(3.*sqrt(y)*sqrt(1.-y)/sqrt(4.-3.*y),0., &
                  (2.-3.*y)/sqrt(4.-3.*y))
             modVp = sqrt((2.*bigG*mStar / rStar)*(rStar/r - rStar/rM))

             vP = (-1.*(modVp/cSpeed)) * vP

             rho = (mdot * rStar)/(4.*pi*(rStar/rInner - rStar/rOuter))
             rho = rho * r**(-5./2.) / sqrt(2.*bigG * mStar)
             rho = rho * sqrt(4.- 3.*y) / sqrt(1.-min(y,0.99))

             posVec = VECTOR(r * sin(theta), 0., r*cos(theta)) / 1.e10
             posVec = rotateZ(posVec, phi)
             vP = rotateZ(vP, phi)
             if (theta > pi/2.) vP%z = -vP%z

             call getIndices(grid, posVec, i1, i2, i3, t1, t2 ,t3)
             grid%velocity(i1,i2,i3) = vP
             grid%rho(i1,i2,i3) = rho
             minRho = min(minRho, rho/REAL(mHydrogen))
             maxRho = max(maxRho, rho/REAL(mHydrogen))
             grid%inUse(i1,i2,i3) = .true.
          enddo
       enddo
    enddo

    grid%temperature = 6500.

    do i = 1, grid%na1
       do j = 1, grid%na2
          do k = 1, grid%na3
             if (grid%inUse(i,j,k).and.(.not.grid%inStar(i,j,k))) then
                grid%temperature(i,j,k) = max(5000.d0,7500.-2500.*(grid%rho(i,j,k)/mHydrogen - minRho)/(maxRho - minRho))
             endif
          enddo
       enddo
    enddo


    write(*,*) "rho",minRho,maxRho

    Laccretion = (bigG*mStar*mDot/rStar)*(1.-(2.*rStar/(rOuter+rInner)))
    Taccretion = Laccretion / ((fourPi * rStar**2)*stefanBoltz* &
         abs(cos(theta1)-cos(theta2)))

    sAccretion = (fourPi * rStar**2)*abs(cos(theta1)-cos(theta2))/1.e20
    Taccretion = Taccretion**0.25

    write(*,*) "accretion lum/temp",Laccretion/Lsol, Taccretion

    open(20,file=contfile1,status="old",form="formatted")
    nnu = 1
10  continue
    read(20,*,end=20) nuarray(nnu), fnu(nnu)
    nnu = nnu + 1
    goto 10
20  continue
    nnu = nnu  - 1
    close(20)
    tot = 0.
    do i = 1, nnu-1
       tot = tot + 0.5*(nuArray(i+1)-nuArray(i))*(fnu(i+1)+fnu(i))
    enddo
    write(*,*) (fourPi*rStar**2)*tot/lSol," solar luminosities"   

    write(*,*) (fourPi*rStar**2)*tot/(fourPi*rStar**2*stefanBoltz*4000.**4)

    ! add the accretion luminosity spectrum to the stellar spectrum,
    ! write it out and pass it to the stateq routine.
    
    open(22,file="star_plus_acc.dat",form="formatted",status="unknown")
    do i = 1, nNu
       fNu(i) = fNu(i) + blackbody(tAccretion, 1.e8*REAL(cSpeed)/ nuArray(i))*(1.e20*sAccretion/(fourPi*rStar**2))
       write(22,*) nuArray(i), fNu(i)
    enddo
    close(22)
    newContfile="star_plus_acc.dat"



!    write(*,*) "adding wind..."
!    do i = 1, grid%nx
!       do j = 1, grid%ny
!          do k = 1, grid%nz
!             rVec = VECTOR(grid%xAxis(i), grid%yAxis(j), grid%zAxis(k))
!             rHat = rVec
!             call normalize(rHat)
!             r = modulus(rVec)
!             theta = acos(rVec%z/r)
!             if ((theta < thetaStart).or.(theta > thetaEnd)) then
!                rMag = 1.05 * rStar/1.e10
!             else
!                rMag = 1.05 * (rOuter*sin(theta)**2)/1.e10
!             endif
!             if (r > rMag) then
!                dist = r/(rStar/1.e10)
!                call hunt(rCalvet, nCalvet, dist, iCalvet)
!                
!                grid%temperature(i,j,k)= tCalvet(iCalvet) + &
!                      ((dist-rCalvet(iCalvet)) / &
!                     (rCalvet(iCalvet+1)-rCalvet(iCalvet))) * &
!                     (tCalvet(iCalvet+1)-tCalvet(iCalvet))
!
!
!                vel = v0 + (vTerm-v0)*(1.-(rStar/1.e10)/r)
!                grid%rho(i,j,k) = mLoss / (fourPi*(r*1.e10)**2*vel)
!                grid%velocity(i,j,k) = (vel/cSpeed) * rHat
!                grid%inUse(i,j,k) = .true.
!             endif
!          enddo
!       enddo
!    enddo
!    write(*,*) "done."




    call initgridstateq(grid, newContfile, contfile2, popFileName, &
         readPops, writePops, lte, nLower, nUpper)


    grid%etaLine = grid%etaLine * 0.2
    
    nu = cSpeed/(lamLine * angstromtocm )


    if (curtains) then
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                if (grid%chiLine(i,j,k) == 0.) grid%chiLine(i,j,k) = 1.e-30
                ang = atan2(grid%yAxis(j),grid%xAxis(i))
                if (ang < 0.) ang = ang + twoPi
                ang = ang * radToDeg
                if (((ang > 30.).and.(ang < 150.)).or.((ang > 210.).and.(ang < 330.))) then
                   grid%etaLine(i,j,k) = 1.e-30
                   grid%etaCont(i,j,k) = 1.e-30
                   grid%chiLine(i,j,k) = 1.e-30
                   grid%kappaAbs(i,j,k,1) = 1.e-30
                   grid%kappaSca(i,j,k,1) = 1.e-30
                   grid%inUse(i,j,k) = .false.
                endif
             enddo
          enddo
       enddo
    endif

    write(*,'(a)') "Computing 3D bias distribution..."
    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             if (grid%inUse(i,j,k).and. (.not.grid%inStar(i,j,k))) then
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))
                grid%biasLine3D(i,j,k) = 1. !max(1.e-5, beta_mn(nLower, nUpper, rVec, i, j, k, grid))
                grid%biasCont3D(i,j,k) = 1. !max(1.e-2, beta_mn(nLower, nUpper, rVec, i, j, k, grid))
             endif
          enddo
       enddo
    enddo
    write(*,'(a)') "done."


  end subroutine fillGridMagneticAccretion
            
  subroutine infallEnhancment(grid, distortionVec, nVec, nPhi, timeStep, doDistortion,&
                              particleMass, alreadyDoneInfall)

    type(GRIDTYPE), intent(inout) :: grid
    integer, intent(in) :: nVec, nPhi
    logical, intent(in) :: doDistortion
    type(VECTOR), intent(inout) :: distortionVec(nVec)
    logical, allocatable :: done(:,:,:)
    real, intent(in) :: timeStep
    real, intent(in) :: particleMass
    logical, intent(inout) :: alreadyDoneInfall
    real :: dTime
    real, parameter :: etaFac = 9
    real, parameter :: chiFac = 1
    type(VECTOR) :: thisVel, thisVec

    integer, parameter :: nTimes = 1000
    real :: phi
    integer :: i, j
    integer :: i1, i2, i3
    real :: t1, t2, t3

    if (grid%adaptive) then 
       call infallEnhancmentAMR(grid, distortionVec, nVec, timeStep, doDistortion, &
                               particleMass,alreadyDoneInfall)    
    else
      
      write(*,*) "Time stepping vectors..."
      dTime = timeStep/real(nTimes)
      do j = 1, nVec
        do i = 1, nTimes
          thisVec = distortionVec(j)/1.e10
          call getIndices(grid, thisVec, i1, i2, i3, t1, t2, t3)
          thisVel = interpGridVelocity(grid,i1, i2, i3, t1, t2, t3)
          thisVel = cSpeed * thisVel
          distortionVec(j) = distortionVec(j) + (dTime * thisVel)
        enddo
      enddo
      write(*,*) "done."

      allocate(done(1:grid%nx, 1:grid%ny, 1:grid%nz))
      done = .false.


      if (doDistortion) then
        write(*,*) "Distorting grid..."    
        do i = 1, nVec
          do j = 1, nPhi
             phi = twoPi * real(j-1)/real(nPhi-1)
             thisVec = rotateZ(distortionVec(i), phi)
             thisVec = thisVec / 1.e10
             call getIndices(grid, thisVec, i1, i2, i3, t1, t2, t3)
             if ((.not.done(i1,i2,i3)).and.(.not.grid%inStar(i1,i2,i3)).and. &
                  grid%inUse(i1,i2,i3)) then
                grid%etaLine(i1,i2,i3) = grid%etaLine(i1,i2,i3) * etaFac
                grid%chiLine(i1,i2,i3) = grid%chiLine(i1,i2,i3) * chiFac
                done(i1,i2,i3) = .true.
             endif
             thisVec%z = -thisVec%z
             call getIndices(grid, thisVec, i1, i2, i3, t1, t2, t3)
             if ((.not.done(i1,i2,i3)).and.(.not.grid%inStar(i1,i2,i3)).and. &
                  grid%inUse(i1,i2,i3)) then
                grid%etaLine(i1,i2,i3) = grid%etaLine(i1,i2,i3) * etaFac
                grid%chiLine(i1,i2,i3) = grid%chiLine(i1,i2,i3) * chiFac
                done(i1,i2,i3) = .true.
             endif
          enddo
        enddo
        write(*,*) "done."
      endif
      deallocate(done)

    end if

  end subroutine infallEnhancment
       

  subroutine initInfallEnhancement(distortionVec, nVec, nPhi, particleMass)
    integer, intent(in) :: nVec, nPhi
    type(VECTOR), intent(inout) :: distortionVec(nVec)
    real, intent(inout) :: particleMass
    real ::  rStar, rinner, rOuter, fac, theta, rM
    integer :: j

    rStar = 2.*rSol
    rInner = 2.2 * rStar
    rOuter = 3. * rStar
    theta = 85.*degtoRad

    write(*,*) "initializing vectors..."
    do j = 1, nVec
       fac = real(j-1)/real(nVec-1)
       rM = rInner + (rOuter - rInner)*fac
       distortionVec(j) = VECTOR(rM * sin(theta), 0., rM * cos(theta))
       write(*,*) distortionVec(j)
    enddo
    write(*,*) "done."

    !particleMass =  sqrt(2.*TTauriRinner/((bigG * TTauriMstar) / TTauriRinner**2)) * &
    !                 TTauriMdot * 1.0 !0.05
    !particleMass = particleMass / real(nVec * nPhi)
    !print *, 'Particle Mass = ', particleMass

    particleMass = 1.0 ! assiging a dummy value  for safty (RK added this)

  end subroutine initInfallEnhancement
       

  subroutine fillGridTTauriWind(grid,contfile1, popFileName, &
         readPops, writePops, lte, nLower, nUpper)

    type(GRIDTYPE) :: grid
    character(len=*) :: contFile1, popFileName
    character(len=80) :: contfile2
    logical :: readPops, writePops, lte
    integer :: i, j, k
    real :: rStar, vTerm, v0, v, r, mdot
    type(VECTOR) :: rVec, rHat
    integer :: nLower, nUpper
    grid%geometry = "ttwind"


    grid%lineEmission = .true.

    


    rStar = sqrt((20.*lSol)/(fourPi * stefanBoltz * 4500.**4))

    write(*,*) "rStar: ",rStar/rSol

    grid%rCore = rStar / 1.e10

    v0 = 20. * 1.e5
    vTerm = 320. * 1.e5
    

    mdot = 1.e-7 * mSol / (365.25 * 24. * 3600.)

    do i = 1, grid%nx
       grid%xAxis(i) = 2.*(real(i-1)/real(grid%nx-1)) - 1.
    enddo
    do i = 1, grid%ny
       grid%yAxis(i) = 2.*(real(i-1)/real(grid%ny-1)) - 1.
    enddo
    do i = 1, grid%nz
       grid%zAxis(i) = 2.*(real(i-1)/real(grid%nz-1)) - 1.
    enddo

    grid%xAxis = grid%xAxis * 10.*grid%rCore
    grid%yAxis = grid%yAxis * 10.*grid%rCore
    grid%zAxis = grid%zAxis * 10.*grid%rCore

    grid%inStar = .false.
    grid%inUse = .false.

    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             rVec = VECTOR(grid%xAxis(i), grid%yAxis(j), grid%zAxis(k))
             r = modulus(rVec)
             if (r > grid%rCore) then
                rHat = rVec
                call normalize(rHat)
                
                v = v0 + (vTerm-v0) * (1. - (grid%rCore/r)**2)
                
                grid%velocity(i, j, k) = (v / cSpeed) * rHat
                grid%rho(i, j, k) = mDot / (fourPi * (r*1.e10)**2 * v)
                grid%inUse(i,j,k) = .true.
             else
                grid%inStar(i,j,k) = .true.
             endif
          enddo
       enddo
    enddo
    
    grid%temperature = 8000.

    call initgridstateq(grid, contfile1, contfile2, popFileName, &
         readPops, writePops, lte, nLower, nUpper)


    write(*,'(a)') "Computing 3D bias distribution..."
    do i = 1, grid%nx
       do j = 1, grid%ny
          do k = 1, grid%nz
             if (grid%inUse(i,j,k)) then
                rVec = VECTOR(grid%xAxis(i),grid%yAxis(j),grid%zAxis(k))
                grid%biasLine3D(i,j,k) = max(1.e-6, sqrt(beta_mn(nLower, nUpper, rVec, i, j, k, grid)))
                grid%biasCont3D(i,j,k) = max(1.e-6, sqrt(beta_mn(nLower, nUpper, rVec, i, j, k, grid)))
             endif
          enddo
       enddo
    enddo
    write(*,'(a)') "done."




  end subroutine fillGridTTauriWind


  subroutine fillGridFlaredDisk(grid,meanParticleMass)
    type(GRIDTYPE) :: grid
    real :: diskMass, bigH, rho, bigR
    real :: rOuter, rHole
    integer :: i, j, k, nMu, n1, n2
    real :: rhoNought, sinTheta
    type(VECTOR) :: rVec
    real :: scaleFac, thisDiskMass
    real :: meanParticleMass

    grid%rho = 1.e-30

    diskMass = 0.01 * 0.01 * mSol
    rOuter = 100. * auTocm
    rHole = 10. * rSol
    grid%rCore = rSol

    nMu = grid%nMu/2

    if (nMu == nint(real(grid%nMu)/2.)) then
       write(*,'(a)') "! nmu should be an odd number"
       stop
    endif

    do i = 1, nMu
       grid%muAxis(i) = -2. *real(i-1)/real(nMu-1) 
    enddo
    grid%muAxis(1:nMu) = 10.**grid%muAxis(1:nMu) 
    grid%muAxis(nMu+1) = 0.
    do i = nMu+2, grid%nMu
       grid%muAxis(i) = -grid%muAxis(grid%nMu-i+1)
    enddo
    grid%muAxis = -grid%muAxis

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi * real(i-1)/real(grid%nphi-1)
    enddo

    do i = 1, grid%nr/2
       grid%rAxis(i) = real(i - 1)/real(grid%nr/2 - 1)
    enddo
    grid%rAxis(1:(grid%nr/2)) = (10.**grid%rAxis(1:(grid%nr/2)) - 1.)/9.
    grid%rAxis(1:(grid%nr/2)) = grid%rAxis(1:(grid%nr/2)) * (4.*rHole) + rHole
    n1 = grid%nr/2+1
    n2 = grid%nr
    do i = n1, n2
       grid%rAxis(i) = real(i - n1)/real(n2 - n1)
    enddo
    grid%rAxis(n1:n2) = (10.**grid%rAxis(n1:n2) - 1.)/9.
    grid%rAxis(n1:n2) = grid%rAxis(n1:n2) * (rOuter - 5.*rHole) + 6.*rHole


    call writeAxes(grid)

    rhoNought = 1.e-10
    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             sinTheta = sqrt(1.-grid%muAxis(j)**2)
             rVec = VECTOR(grid%rAxis(i)*cos(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*sin(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*grid%muAxis(j))
             bigR = sqrt(rVec%x**2 + rVec%y**2)
             if (bigR > rHole) then
                bigH = 0.05 * bigR
                rho = rhoNought * (bigR/grid%rAxis(1))**(-3./2.)
                rho = rho * exp ( -(rVec%z**2/(2.*bigH**2)))
                grid%rho(i,j,k) = rho
             endif
          enddo
       enddo
    enddo
    
    thisDiskMass = integratedDensity(grid)

    scaleFac = diskMass/thisDiskMass

    grid%rho= grid%rho * scalefac

    write(*,*) "Disk mass (solar) :",integratedDensity(grid)/mSol
    grid%rho = grid%rho / meanParticleMass

  end subroutine fillGridFlaredDisk

  

end module ttauri_mod


