module timedep_mod

! time dependent RT module

  use kind_mod
  use constants_mod
  use vector_mod
  use utils_mod

  implicit none

  private
  
  public :: timeDependentRT

contains

  subroutine timeDependentRT()
    real(double) :: currentTime, endTime
    integer, parameter :: nx = 100
    integer, parameter :: nProb = nx + 1
    real(double) :: xCen(nx)
    real(double) :: etaCont(nx)
    real(double) :: aDot(nx)
    real(double) :: rho(nx)
    real(double) :: kappa(nx)
    real(double) :: temperature(nx)
    real(double) :: distanceGrid(nx)
    real(double) :: prob(nProb)
    real(double) :: xArray(nProb)
    real(double) :: uDens(nx)
    real(double) :: uDensAnalytical(nx)
    real(double) :: tauBox
    real(double) :: luminosity
    real(double) :: photonTime
    logical :: absorbed, finished
    real(double) :: distToBoundary, timeToBoundary, tauToBoundary, r
    real(double) :: distanceToEvent, deltaT
    real(double) :: tau, dx, k
    real(double) :: xSize
    real(double) :: photonSpeed
    real(double) :: xPhoton
    integer, parameter :: maxStack = 10000000

    integer :: currentNStack, oldnStack
    real(double) :: currenttimeStack(maxStack)
    type(VECTOR) :: currentpositionStack(maxStack)
    type(VECTOR) :: currentdirectionStack(maxStack)
    real(double) :: currentepsStack(maxStack)

    real(double) :: oldtimeStack(maxStack)
    type(VECTOR) :: oldpositionStack(maxStack)
    type(VECTOR) :: olddirectionStack(maxStack)
    real(double) :: oldepsStack(maxStack)

    real(double) :: epsOverDeltaT
    logical :: outOFtime
    type(VECTOR) :: rVec, uHat
    integer :: i, iMonte, nMonte, iPos,iter
    character(len=30) :: outfile

    nMonte = 10000000
    tauBox = 100.d0
    deltaT = 1.d0

    xSize = 1.d0
    photonSpeed = 1.d0
    dx = xSize / dble(nx)
    xCen(1) = dx/2.d0
    do i = 2, nx
       xCen(i) = xCen(i-1) + dx
    enddo
    do i = 1, nx
       !       uDens(i) = cos((pi/2.d0)*xCen(i)/xSize)
       uDens(i) = max(0.d0,cos(2.d0*xCen(i)))
       temperature(i) = (uDens(i)/aRad)**0.25d0
    enddo
    do i = 1, nProb
       xArray(i) = xSize * dble(i-1)/dble(nProb-1)
    enddo


    rho = 1.d0
    kappa = tauBox / (rho*xSize)
    write(*,*) "total tau ",SUM(kappa*rho*dx)

    uDensAnalytical = uDens
    k = cSpeed/(kappa(1)*rho(1))

    currentTime = 0.d0
    endTime = 1.d0

    iter = 0
    currentnStack = 0
    oldnStack = 0
    do while (.true.)
       iter = iter + 1

       deltaT = 0.3d0*(dx**2/(2.d0*k))

       distanceGrid = 0.d0

       call solveNewUdens(xCen, uDensAnalytical, nx, k, deltaT)

       etaCont = fourPi * stefanBoltz * temperature**4 * kappa * rho

       call calculateProb(prob, xArray, nProb, xCen, etaCont, nx)
       luminosity =SUM(etaCont*dx)
       epsOverDeltaT = luminosity / dble(nMonte)
       write(*,*) "energy monte ",SUM(uDens*dx)
       write(*,*) "energy analy ",SUM(uDensAnalytical*dx)
       write(*,*) "light distance ", deltaT*cSpeed

       do iMonte = 1, nMonte + oldnStack

          if (oldnStack > 0)  then
             call getPhotonFromStack(oldnStack, rVec, uHat, photonTime, epsOverDeltaT, &
                  oldpositionStack, olddirectionStack, oldtimeStack, oldepsStack)
          else
             call getX(xArray, prob, nProb, xPhoton)
             rVec%x = xPhoton
             rVec%y = 0.d0
             rVec%z = 0.d0
             photonTime = 0.d0
             uHat = randomUnitVector()
          endif

          absorbed = .false.
          finished = .false.
          outOfTime = .false.
          do while (.not.finished)

             if (absorbed) then
                uHat = randomUnitVector()
                absorbed = .false.
             endif

             call findArrayIndex(xCen, dx, nx, rVec%x, iPos)
             if (uHat%x > 0.d0) then
                distToBoundary = ((xCen(iPos)+dx/2.d0) - rVec%x)/uHat%x
             else
                distToBoundary = abs( ((xCen(iPos)-dx/2.d0) - rVec%x )/uHat%x )
             endif

             timeToBoundary = distToBoundary / cSpeed

             tauToBoundary = distToBoundary * kappa(iPos) * rho(iPos)

             !             write(*,*) "tautoboundary ", tautoboundary, iPos, rVec%x, xCen(iPos)-dx/2.d0, xCen(iPos)+dx/2.d0
             call random_number(r)
             tau = -log(1.d0 - r)


             if (tau > tauToBoundary) then
                distanceToEvent = distToBoundary
                absorbed = .false.
             else
                distanceToEvent = distToBoundary * tau/tauToBoundary
                absorbed = .true.
                finished = .true.
             endif

!             if (absorbed) then
!                if ((photonTime + distanceToEvent/cSpeed) > deltaT) then
!                   absorbed = .false.
!                   distanceToEvent = (deltaT - photonTime) * cSpeed
!                   outOfTime = .true.
!                   finished = .true.
!                endif
!             else
!                if ((photonTime + distanceToEvent/cSpeed) > deltaT) then
!                   distanceToEvent = (deltaT - photonTime) * cSpeed
!                   outOfTime = .true.
!                   finished = .true.
!                endif
!             endif

             photonTime = photonTime + distanceToEvent / cSpeed

             distanceGrid(iPos) = distanceGrid(iPos) + &
                  distanceToEvent * kappa(iPos) * rho(iPos) * epsOverDeltaT

             rVec = rVec + (distanceToEvent + 1.d-6) * uHat
             if (rVec%x < xCen(1)-dx/2.d0) then
                rVec%x =  (xCen(1)-dx/2.d0)+1.d-6
                uHat%x = -1.d0 * uHat%x
             endif
             if (rVec%x > xCen(nx)+dx/2.d0) then
                rVec%x = (xCen(nx)+dx/2.d0)-1.d-6
                uHat%x = -1.d0 * uHat%x
             endif
          end do
          if (outOfTime) then
             call addPhotonToStack(currentnStack, rVec, uHat, photonTime, &
                  epsOverDeltaT, currentpositionStack, currentdirectionStack, &
                  currenttimeStack, currentepsStack)
          endif
       end do
       call calculateADot(distanceGrid, aDot, xCen, dx, nx)

       do i = 1, nx
          uDens(i) = uDens(i) + (aDot(i) - etaCont(i)) * deltaT
          uDens(i) = max (0.d0, uDens(i))
          temperature(i) = (uDens(i)/aRad)**0.25d0
       enddo

       write(outfile,'(a,i6.6,a)') "udens",iter,".dat"
       write(*,*) iter
       open(21, file=outfile, status="unknown", form="formatted")
       do i = 1, nx
          write(21, *) xCen(i), uDens(i), uDensAnalytical(i)
       enddo
       close(21)

       oldnStack = currentnStack
       oldPositionStack(1:currentnStack) = currentPositionStack(1:currentnStack)
       oldDirectionStack(1:currentnStack) = currentDirectionStack(1:currentnStack)
       oldTimeStack(1:currentnStack) = currentTimeStack(1:currentnStack)
       oldEpsStack(1:currentnStack) = currentEpsStack(1:currentnStack)
    enddo
  end subroutine timeDependentRT

  subroutine addPhotonToStack(nStack, rVec, uVec, photonTime, &
       epsOverDeltaT, positionStack, directionStack, timeStack, epsStack)
    integer :: nStack
    real(double) :: photonTime, epsOverDeltaT
    type(VECTOR) :: rVec, uVec, positionStack(:), directionStack(:)
    real(double) :: timeStack(:), epsStack(:)
    
    nStack = nStack + 1
!    write(*,*) "adding photon to stack ", nStack
    if (nStack > SIZE(positionStack)) then
       write(*,*) "nStack too large for array ",nStack
       stop
    endif
    positionStack(nStack) = rVec
    directionStack(nStack) = uVec
    timeStack(nStack) = photonTime
    epsStack(nStack) = epsOverDeltaT
  end subroutine addPhotonToStack

  subroutine getPhotonFromStack(nStack, rVec, uHat, photonTime, epsOverDeltaT, &
                  positionStack, directionStack, timeStack, epsStack)
    integer :: nStack
    type(VECTOR) :: rVec
    type(VECTOR) :: uHat
    real(double) :: photonTime
    real(double) :: epsOverDeltaT
    type(VECTOR) :: positionStack(:)
    type(VECTOR) :: directionStack(:)
    real(double) :: timeStack(:)
    real(double) :: epsStack(:)

!    write(*,*) "getting photon from stack ",nStack
    rVec = positionStack(nStack)
    uHat = directionStack(nStack)
    photonTime = timeStack(nStack)
    epsOverDeltaT = epsStack(nStack)
    nStack = nStack - 1
  end subroutine getPhotonFromStack


  subroutine calculateADot(distanceGrid, aDot, xCen, dx, nx)
    real(double) :: distanceGrid(:), aDot(:), xCen(:), dx
    integer :: nx, i


    do i = 1, nx
       aDot(i) = (1.d0 / dx) * distancegrid(i)
    enddo
  end subroutine calculateADot

  subroutine calculateProb(prob, xArray, nProb, xCen, etaCont, nx)
    real(double) :: prob(:), xCen(:), etaCont(:), xArray(:)
    integer :: nx, i, nProb

    prob(1) = 0.d0
    do i = 2, nProb
       prob(i) = prob(i-1) + etaCont(i-1) * (xCen(2)-xCen(1))
    enddo
    prob(1:nProb) = prob(1:nProb) / prob(nProb)
  end subroutine calculateProb
       
  subroutine getX(xArray, prob, nProb, x)
    real(double) :: xArray(:), prob(:), x
    integer :: nProb
    real(double) :: r
    integer :: i

    call random_number(r)
    call locate(prob, nprob, r, i)
    x = xArray(i) + ((r - prob(i))/(prob(i+1)-prob(i)))*(xArray(i+1) - xArray(i))
  end subroutine getX

  subroutine findArrayIndex(xCen, dx, nx, x, iPos)
    real(double) :: xCen(:), dx, x
    integer :: iPos, nx
    if (x < xCen(1)-dx/2.d0) then
       write(*,*) "x value outside array ",x
       stop
    endif
    if (x > xCen(nx)+dx/2.d0) then
       write(*,*) "x value outside array ",x
       stop
    endif
    if (x < xCen(1)) then
       iPos = 1
    else if (x > xCen(nx)) then
       iPos = nx
    else
       iPos = nint((x-dx/2.d0)/dx)+1
    endif
  end subroutine findArrayIndex

  subroutine solveNewUdens(xCen, u, nx, k, deltaT)
    real(double) :: xCen(:), u(:), k, deltaT,dx
    real(double) :: uxx, uPrime(1000)
    integer :: nx, i
    dx = xCen(2)-xCen(1)
    

!    uxx = (u(3) - 2.d0* u(2) + u(1)) / dx**2

!    uprime(1) = u(1) + k * uxx * deltaT

!    uxx = (u(nx) - 2.d0* u(nx-1) + u(nx-2)) / dx**2

!    uprime(nx) = u(nx) + k * uxx * deltaT

    do i = 2, nx-1
       uxx = (u(i+1) - 2.d0* u(i) + u(i-1)) / dx**2
       uprime(i) = u(i) + k * uxx * deltaT
    enddo
    uprime(1) = uprime(2)
    uprime(nx) = uprime(nx-1)

    u = uprime
  end subroutine solveNewUdens

end module timedep_mod
