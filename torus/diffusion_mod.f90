! this is a module to apply the diffusion approximation to
! the very optically thick central regions of circumstellar stellar discs.


module diffusion_mod

use constants_mod
use gridtype_mod
use amr_mod
use vector_mod

implicit none


contains




  subroutine solveDiffusion(grid, zArray, xPos, temperature, rho,  diffApprox, nz)
    implicit none
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: nz
    real :: zArray(*)
    real :: xPos
    type(OCTALVECTOR) :: octVec
    real :: temperature(*)
    real :: rho(*)
    real :: flux(1000)
    real :: jd(1000), kappa_ross(1000), kappap(1000)
    real :: tau(1000)
    logical :: diffApprox(*)
    real, allocatable :: dTdZ(:), newTemperature(:), frac(:)
    real :: kappa
    integer :: i, nIter
    logical :: converged
    real(double) :: tReal(1000),zReal(1000),rhoReal(1000)
    integer :: nReal
    real :: bigJ
    real :: dtdz1, dfdz(1000)
    integer :: istart
    real :: safetemp
    real :: dJdz(1000)
    real :: t0, t1, fac

    allocate(dTdZ(1:nz))
    allocate(newTemperature(1:nz))
    allocate(frac(1:nz))

! now take the first nreal points in the temperature/z run
! which correspond to temperatures that have be calculated
! via the lucy algorithm

    write(*,*) "starting diffusion calculation..."

    nReal = 0
    do i = 1, nz
       !          if (.not.diffusionApprox(i)) then
       if (temperature(i) > 100.) then
          nReal = nReal + 1
          tReal(nReal) = temperature(i)
          zReal(nReal) = zArray(i)
          rhoReal(nReal) = rho(i)
       else
          iStart = i ! first element of diffusion array
          exit
       endif
    enddo
    
 
    dTdz1 = getTempGradient(nReal, tReal, zReal)
!    write(*,*) (treal(nreal)-treal(nreal-1))/(zReal(nreal)-zreal(nreal-1)) /1.e10
!    dtdz1 = (treal(nreal)-treal(nreal-1))/(zReal(nreal)-zreal(nreal-1)) /1.e10 
    
    converged = .false.

    octVec = VECTOR(xPos, 0., zReal(nReal))
    call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
         foundSubcell=subcell, rosselandKappa=kappa, grid=grid, ilambda = 32)

    flux(iStart-1) = -16.* stefanBoltz * treal(nreal)**3 * dtdz1 / (3.* kappa * rhoreal(nreal))
    write(*,*) flux(istart-1)

!    flux(istart-1) = -stefanBoltz*treal(nreal)**4 
    jd(iStart-1) = flux(iStart-1) / twoPi


    t0 = 100.
    t1 = treal(nreal)
    fac = log(t1/t0)/zreal(nreal)
    do i = iStart,nz
       temperature(i) = t0 * exp(fac*zArray(i))
    enddo


    nIter = 0
    do while(.not.converged)
       

       do i = 1, nz+1
          newTemperature(i) = temperature(i)
       enddo

       
       ! first calculate the opacities

       do i = iStart, nz+1
          octVec = VECTOR(xPos, 0., zArray(i))
          call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
               foundSubcell=subcell, rosselandKappa=kappa_ross(i), kappap=kappap(i), grid=grid, &
               atthistemperature=newtemperature(i))
          write(*,*) i,zArray(i),kappa_ross(i),kappap(i)
       enddo
       

       flux(nz+1) = 0.
       jd(nz+1) = -0.9*(stefanBoltz*newtemperature(nz+1)**4 / pi)
       djdz(nz+1) = -3. * kappa_ross(nz+1) * rho(nz+1) * flux(nz+1) / fourPi
       dfdz(nz+1) = fourPi * kappap(nz+1) * rho(nz+1) * ( (stefanBoltz*newtemperature(nz+1)**4 /pi)-jd(nz+1))
       do i = nz, iStart,-1
          flux(i) = flux(i+1) + dfdz(i+1) * (zArray(i)-zArray(i-1)) * 1.e10
          jd(i) = jd(i+1) + djdz(i+1)  * (zArray(i)-zArray(i-1)) * 1.e10
          dfdz(i) = fourPi * kappap(i) * rho(i) * ( (stefanBoltz*newtemperature(i)**4 /pi)-jd(i))
          djdz(i) = 3. * kappa_ross(i) * rho(i) * flux(i) / fourPi
       enddo

!       do i = iStart, nz
!          dfdz(i) = fourPi * kappap(i) * rho(i) * ( (stefanBoltz*newtemperature(i)**4 /pi)-jd(i))
!          flux(i) = flux(i-1) + dfdz(i) * (zArray(i)-zArray(i-1)) * 1.e10
!          djdz(i) = -3. * kappa_ross(i) * rho(i) * flux(i) / fourPi
!          jd(i) = jd(i-1) + djdz(i)  * (zArray(i)-zArray(i-1)) * 1.e10
!       enddo



       temperature(1:nz) = newTemperature(1:nz)

       do i = nz+1,istart,-1
          write(*,'(i3,1p,6e12.3)') i, zArray(i), temperature(i),flux(i), dfdz(i),jd(i), djdz(i)
       enddo


       nIter = nIter + 1000
       if (nIter > 10) converged = .true.
       stop
    enddo

    deallocate(newTemperature, frac, dTdZ)

  end subroutine solveDiffusion

  subroutine throughoutMidplaneDiff(grid)

    type(GRIDTYPE) :: grid
    real :: zAxis(10000), rho(10000), temperature(10000), subcellsize(10000)
    integer :: nz
    real :: xpos, ypos, radius, drho, smallestSubcell
    logical :: diffApprox(10000)
    logical :: converged 
    real :: xAxis(100000)
    integer :: nx, i
    real :: flux

    nx = 0
    call getxValuesdiff(grid%octreeRoot,nx,xAxis)
    call stripSimilarValues(xAxis,nx,real(1.e-5*grid%halfSmallestSubcell))
    xAxis(1:nx) = xAxis(1:nx) + 1.e-5*grid%halfSmallestSubcell

    xPos = grid%halfSmallestSubcell - grid%octreeRoot%subcellSize
    smallestSubcell = 2. * grid%halfSmallestSubcell
    
    xPos = grid%halfSmallestSubcell
    yPos = 0.

    do i = 1, nx
       radius = xAxis(i)
       xPos = xAxis(i)
       call getTemperatureDensityRunDiff(grid, zAxis, subcellsize, rho, temperature, diffApprox, xPos, yPos, nz, -1.)

       if (nz > 1) then
          call solveDiffusion(grid, zAxis, xPos, temperature, rho, diffapprox, nz)
!          

          call putTemperatureRunDiff(grid, zAxis, temperature, nz, xPos, yPos, -1.)
       endif

!       call getTemperatureDensityRunDiff(grid, zAxis, subcellsize, rho, temperature, diffApprox, xPos, yPos, nz, +1.)
!
!       if (nz > 1) then
!          call solveDiffusion(grid, zAxis, xPos,  temperature, rho, diffApprox, nz, flux)
!          
!          call putTemperatureRunDiff(grid, zAxis, temperature, nz, xPos, yPos, +1.)
!       endif

    enddo

    call copyChilineTotemperature(grid%octreeRoot)

  end subroutine throughoutMidplaneDiff


  subroutine getTemperatureDensityRunDiff(grid, zAxis, subcellsize, rho, temperature, diffApprox, xPos, yPos, nz, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer :: nz
    real :: rho(*), temperature(*), zAxis(*), subcellsize(*)
    logical :: diffApprox(*)
    real :: xPos, yPos
    integer :: subcell
    real :: rhotemp, temptemp
    real :: direction
    type(OCTALVECTOR) :: currentPos, temp
    real :: halfSmallestSubcell
    logical :: someDiffusion
    integer :: i

    nz = 0
    halfSmallestSubcell = grid%halfSmallestSubcell

    currentPos = OCTALVECTOR(xPos, yPos, -1.*direction*grid%ocTreeRoot%subcellsize)

    do while((-1.*direction*currentPos%z) > 0. )
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp, temperature=temptemp)
          temp = subCellCentre(thisOctal, subcell)
          if (thisOctal%inflow(subcell)) then
             nz = nz + 1
             temperature(nz) = temptemp
             diffApprox(nz) = thisOctal%diffusionApprox(subcell)
             rho(nz) = rhotemp
             zAxis(nz) = temp%z
             subcellsize(nz) = thisOctal%subcellsize
          endif
          currentPos = OCTALVECTOR(xPos, yPos, temp%z+0.5*direction*thisOctal%subcellsize+direction*halfSmallestSubcell)
    end do


    someDiffusion = .false.
    do i = 1, nz
       if (diffApprox(i)) someDiffusion = .true.
    enddo
    if (.not.someDiffusion) nz = 0
    

  end subroutine getTemperatureDensityRunDiff

  subroutine putTemperatureRunDiff(grid, zAxis, temperature, nz, xPos, yPos, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer :: nz
    real :: xPos, yPos
    real :: zAxis(:), temperature(:)
    integer :: i, subcell
    real :: direction
    type(octalvector) :: currentPos
    
    do i = 1, nz
       currentPos = OCTALVECTOR(xPos, yPos, zAxis(i))
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell)
       if (thisOctal%diffusionApprox(subcell)) then
          thisOctal%chiline(subcell) = temperature(i) ! use chiline as temporary temperature storage
       endif
    enddo
  end subroutine putTemperatureRunDiff

  recursive subroutine getxValuesdiff(thisOctal, nx, xAxis)

    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(octalvector) :: rVec
    integer :: nx, subcell, i
    real :: xAxis(:)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getxValuesdiff(child, nx, xAxis)
                exit
             end if
          end do
       else

          rVec = subcellCentre(thisOctal, subcell)
          nx = nx + 1
          xAxis(nx) = rVec%x
       end if
    end do

  end subroutine getxValuesdiff

  subroutine defineDiffusionZone(grid)

    type(GRIDTYPE) :: grid
    real :: zAxis(10000), rho(10000), temperature(10000), subcellsize(10000)
    integer :: nz
    real :: xpos, ypos, radius, drho, smallestSubcell
    logical :: converged 
    real :: xAxis(1000000)
    integer :: nx, i
    real :: flux
    real :: rosselandOpticalDepth(10000)
    integer :: j
    type(OCTALVECTOR) :: octVec
    real :: kappa
    type(OCTAL), pointer :: thisOctal, boundaryOctal
    integer :: subcell
    real :: xStart
    integer :: iBoundary, boundarySubcell
    real :: diffDepth = 10
    real(double) :: tot


    call getxAxisRun(grid, xAxis, subcellSize, 0., 0., nx, +1.)


    write(*,*) "Defining diffusion zone..."

    call zeroDiffusionProb(grid%octreeRoot)

    rosselandOpticalDepth(1) = 0.
    do j = 2, nx
       octVec = OCTALVECTOR(xAxis(j), 0., 0.)
       call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
            foundSubcell=subcell, rosselandKappa=kappa, grid=grid)
       rosselandOpticalDepth(j) = rosselandOpticalDepth(j-1) + kappa * thisOctal%rho(subcell) *  subcellsize(j)*1.e10
       if (rosselandOpticalDepth(j) > diffDepth) then
          xStart = xAxis(j)
          exit
       endif
    enddo

    

    nx = 0
    call getxValuesdiff(grid%octreeRoot,nx,xAxis)
    call stripSimilarValues(xAxis,nx,real(1.e-5*grid%halfSmallestSubcell))
    xAxis(1:nx) = xAxis(1:nx) + 1.e-5*grid%halfSmallestSubcell

    call locate(xAxis, nx, xStart, j)

    iBoundary = j-1

    do i = j, nx
       xPos = xAxis(i)

       call getzAxisRun(grid, zAxis, subcellSize, xPos, yPos, nz, -1.)
       

       rosselandOpticalDepth(1) = 0.
       do j = 2, nz
          octVec = OCTALVECTOR(xpos, 0., zAxis(j))
          call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
               foundSubcell=subcell, rosselandkappa=kappa, grid=grid)
          rosselandOpticalDepth(j) = rosselandOpticalDepth(j-1) + kappa * thisOctal%rho(subcell) * subcellsize(j)*1.e10
          if (rosselandOpticalDepth(j) > diffDepth) then
             thisOctal%diffusionApprox(subcell) = .true.

             ! put in the lefthand boundary of the diffusion zone if necessary

             if (i == (iBoundary + 1)) then
                octVec = OCTALVECTOR(xAxis(iBoundary), 0., zAxis(j))
                call amrGridValues(grid%octreeRoot, octVec, foundOctal=boundaryOctal, &
                     foundSubcell=boundarySubcell)
                call putDiffusionProb(grid, boundaryOctal, boundarySubcell)
             endif
          else
             thisOctal%diffusionApprox(subcell) = .false.
          endif
       enddo

       ! now put in boundary at top of diffusion zone

       call locate(rosselandOpticalDepth, nz, diffDepth, j)
       octVec = OCTALVECTOR(xpos, 0., zAxis(j-1))
       call amrGridValues(grid%octreeRoot, octVec, foundOctal=boundaryOctal, &
            foundSubcell=boundarySubcell)
       call putDiffusionProb(grid, boundaryOctal, boundarySubcell)
       


       call getzAxisRun(grid, zAxis, subcellSize, xPos, yPos, nz, +1.)

       rosselandOpticalDepth(1) = 0.
       do j = 2, nz
          octVec = OCTALVECTOR(xpos, 0., zAxis(j))
          call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
               foundSubcell=subcell, rosselandKappa=kappa, grid=grid, ilambda=32)
          rosselandOpticalDepth(j) = rosselandOpticalDepth(j-1) + kappa * thisOctal%rho(subcell) * subcellsize(j)*1.e10
          if (rosselandOpticalDepth(j) > diffDepth) then
             thisOctal%diffusionApprox(subcell) = .true.
          else
             thisOctal%diffusionApprox(subcell) = .false.
          endif
       enddo


       ! now put in boundary at bottom of diffusion zone

       call locate(rosselandOpticalDepth, nz, diffDepth, j)
       octVec = OCTALVECTOR(xpos, 0., zAxis(j-1))
       call amrGridValues(grid%octreeRoot, octVec, foundOctal=boundaryOctal, &
            foundSubcell=boundarySubcell)
       call putDiffusionProb(grid, boundaryOctal, boundarySubcell)

    enddo

    tot = 0.d0
    call sumDiffusionProb(grid%octreeRoot, tot)
    call normDiffusionProb(grid%octreeRoot, tot)


    write(*,*) "Done."
  end subroutine defineDiffusionZone



  subroutine getxAxisRun(grid, xAxis, subcellSize, zPos, yPos, nx, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer :: nx
    real :: xAxis(:), subcellsize(:)
    real :: zPos, yPos
    integer :: subcell
    real :: rhotemp, temptemp
    real :: direction
    type(OCTALVECTOR) :: currentPos, temp
    real :: halfSmallestSubcell

    nx = 0
    halfSmallestSubcell = grid%halfSmallestSubcell

    currentPos = OCTALVECTOR(0., yPos, zPos)

    do while((currentPos%x) < grid%octreeRoot%subcellSize )
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp, temperature=temptemp, grid=grid)
       nx = nx + 1
       temp = subCellCentre(thisOctal, subcell)
       xAxis(nx) = temp%x
       subcellsize(nx) = thisOctal%subcellsize
       currentPos = OCTALVECTOR(xAxis(nx)+0.5*direction*thisOctal%subcellsize+direction*halfSmallestSubcell, yPos, zPos)
    end do
    
 
  end subroutine getxAxisRun

  subroutine getzAxisRun(grid, zAxis, subcellSize, xPos, yPos, nz, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer :: nz
    real :: zAxis(:), subcellsize(:)
    real :: xPos, yPos
    integer :: subcell
    real :: rhotemp, temptemp
    real :: direction
    type(OCTALVECTOR) :: currentPos, temp
    real :: halfSmallestSubcell

    nz = 0
    halfSmallestSubcell = grid%halfSmallestSubcell

    currentPos = OCTALVECTOR(xpos, yPos, -1.*direction*grid%octreeRoot%subcellsize)

    do while((-1.*direction*currentPos%z) > 0.)
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp, temperature=temptemp, grid=grid)
       nz = nz + 1
       temp = subCellCentre(thisOctal, subcell)
       zAxis(nz) = temp%z
       subcellsize(nz) = thisOctal%subcellsize
       currentPos = OCTALVECTOR(xpos, yPos, zAxis(nz)+0.5*direction*thisOctal%subcellsize+direction*halfSmallestSubcell)
    end do
    
  end subroutine getzAxisRun


  real function getTempGradient(nz, t, z)
    integer :: nz
    real(double) :: t(*), z(*), chisq, apoly(4)
    integer :: nterms
    real(double) :: revT(1000), revZ(1000)
    integer :: i

    ! fit a polynomial

    do i =1, nz
       revZ(i) = z(nz-i+1)
       revT(i) = t(nz-i+1)
    enddo


    nTerms = 3

    call polfit(revz, revt, t, nz ,nterms , -1, apoly, chisq)

    getTempGradient = apoly(2) + 2.*revZ(1)*apoly(3)

    getTempGradient = getTempGradient / 1.e10 ! [K cm-1]

  end function getTempGradient


  recursive subroutine copychilinetoTemperature(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call copychilinetoTemperature(child)
                exit
             end if
          end do
       else
          if (thisOctal%diffusionApprox(subcell)) then
             thisOctal%temperature(subcell) = max(thisOctal%chiline(subcell),3.d0)
          endif
       endif
    enddo
  end subroutine copychilinetoTemperature

  recursive subroutine zeroDiffusionProb(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroDiffusionProb(child)
                exit
             end if
          end do
       else
          thisOctal%diffusionProb(subcell) = 0.d0
          thisOctal%diffusionApprox = .false.
       endif
    enddo
  end subroutine zeroDiffusionProb


  subroutine putDiffusionProb(grid, thisOctal, subcell)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    type(OCTALVECTOR) :: rVec
    integer :: subcell
    integer :: nFreq, iLam
    real(double) :: freq(1000), dnu(1000), thisLam, norm, r1, r2, v, kappaP
    real :: kabsArray(1000)
    integer :: i


    nFreq = grid%nLambda
    do i = 1, nFreq
       freq(nFreq-i+1) = cSpeed / (grid%lamArray(i)*1.e-8)
    enddo


    do i = 2, nFreq-1
       dnu(i) = 0.5*((freq(i+1)+freq(i))-(freq(i)+freq(i-1)))
    enddo
    dnu(1) = freq(2)-freq(1)
    dnu(nFreq) = freq(nFreq)-freq(nFreq-1)

    call amrGridValues(grid%octreeRoot, subcellCentre(thisOctal,subcell), startOctal=thisOctal, &
         actualSubcell=subcell, kappaAbsArray=kAbsArray, grid=grid)


    do i = 1, nFreq
       thisLam = (cSpeed / freq(i)) * 1.e8
       call hunt(grid%lamArray, grid%nLambda, real(thisLam), iLam)
       if ((iLam >=1) .and. (iLam <= grid%nLambda)) then
          kappaP = kappaP + dble(kabsArray(ilam)) * &
               dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell)))) * dble(dnu(i)) 
          norm = norm + dble(bnu(dble(freq(i)),dble(thisOctal%temperature(subcell)))) * dble(dnu(i)) 
       endif
    enddo

    kappaP = kappaP / norm / 1.e10
    rVec = subcellCentre(thisOctal,subcell)
    r1 = rVec%x-thisOctal%subcellSize/2.d0
    r2 = rVec%x+thisOctal%subcellSize/2.d0
    v = dble(pi) * (r2**2 - r1**2) * thisOctal%subcellSize * 1.d30


    thisOctal%diffusionProb(subcell) = fourPi * kappaP * (stefanBoltz/pi) * &
         (thisOctal%temperature(subcell)**4) * V 

  end subroutine putDiffusionProb


  recursive subroutine sumDiffusionProb(thisOctal, tot)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  real(double) :: tot
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sumDiffusionProb(child, tot)
                exit
             end if
          end do
       else
          tot = tot + thisOctal%diffusionProb(subcell)
       endif
       thisOctal%diffusionProb(subcell) = tot
    enddo
  end subroutine sumDiffusionProb

  recursive subroutine normDiffusionProb(thisOctal, tot)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  real(double) :: tot
  
  if (thisOctal%nChildren > 0) then
            
     ! call this subroutine recursively on each of its children
     do subcell = 1, thisOctal%nChildren, 1 
        child => thisOctal%child(subcell)
        call normDiffusionProb(child, tot)
     end do
  end if

  thisOctal%diffusionProb  = thisOctal%diffusionProb /tot


  end subroutine normDiffusionProb



  RECURSIVE SUBROUTINE locateDiffusionProbAMR(probability,thisOctal,subcell) 
    ! finds the subcell that contains a given value of 'probability'.
    ! each subcell of the tree's octals has a value for line emission 
    !   probability which is an upper bound of the subcell's value in the cumulative probability
    !   distribution for the 

    IMPLICIT NONE

    real(double), INTENT(IN) :: probability
    TYPE(octal), POINTER              :: thisOctal
    INTEGER, INTENT(OUT)              :: subcell

    INTEGER              :: i, j 
   
   
    ! we need to treat the first subcell as a special case
    
    IF (probability < thisOctal%diffusionProb(1)) THEN

      IF (thisOctal%hasChild(1)) THEN
      
        ! find the child
        DO j = 1, thisOctal%nChildren, 1
          IF (thisOctal%indexChild(j) == 1) THEN
            thisOctal => thisOctal%child(j)
            CALL locateDiffusionProbAMR(probability,thisOctal,subcell)
            RETURN
          END IF
        END DO
        
      ELSE 
        subcell = 1
        RETURN
        
      END IF
    END IF   
  
   
    DO i = 2, thisOctal%maxChildren, 1
      IF (probability > thisOctal%diffusionProb(i-1) .AND. &
          probability < thisOctal%diffusionProb(i)) THEN
      
        IF (thisOctal%hasChild(i)) THEN
          
          ! find the child
          DO j = 1, thisOctal%nChildren, 1
            IF (thisOctal%indexChild(j) == i) THEN
              thisOctal => thisOctal%child(j)
              CALL locateDiffusionProbAMR(probability,thisOctal,subcell)
              RETURN
            END IF
          END DO
          
        ELSE 
          subcell = i
          RETURN
        
        END IF
      END IF
    END DO

      
  END SUBROUTINE locateDiffusionProbAMR



end module diffusion_mod






