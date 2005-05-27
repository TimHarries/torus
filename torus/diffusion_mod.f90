! this is a module to apply the diffusion approximation to
! the very optically thick central regions of circumstellar stellar discs.


module diffusion_mod

use constants_mod
use gridtype_mod
use amr_mod
use vector_mod

implicit none


contains




  subroutine solveDiffusion(grid, zArray, xPos, temperature, rho,  diffApprox, nz, ok, debugoutput)
    implicit none
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: nz
    real :: zArray(*)
    real :: xPos
    logical :: debugoutput
    type(OCTALVECTOR) :: octVec
    real :: temperature(*)
    real :: rho(*)
    real :: flux(1000)
    real :: jd(1000), kappa_ross(1000), kappap(1000)
    real :: tau(1000)
    logical :: diffApprox(*)
    real(double) :: dTdZ(1000), newTemperature(1000), frac(1000)
    real :: oldVal, thisVal
    real :: kappa
    integer :: i, nIter, j
    logical :: converged
    real(double) :: tReal(1000),zReal(1000),rhoReal(1000)
    integer :: nReal
    real :: bigJ
    real :: dtdz1, dfdz(1000)
    real(double) :: lambdad(1000),gammad(1000),dfraddz(1000)
    integer :: istart
    real :: safetemp
    real(double) :: dJdz(1000)
    real :: nextJd, nextfd, nextfdiff
    real :: t0, t1, fac
    real(double) :: gammairr(1000),frad(1000)
    real(double) :: fd(1000), dfddz(1000), firr(1000)
    real(double) :: jirr(1000)
    real :: influx
    logical :: ok
    
    ok = .true.


! now take the first nreal points in the temperature/z run
! which correspond to temperatures that have be calculated
! via the lucy algorithm


    nReal = 0
    do i = 1, nz
       if (.not.diffApprox(i)) then
          nReal = nReal + 1
          tReal(nReal) = temperature(i)
          zReal(nReal) = zArray(i)
          rhoReal(nReal) = rho(i)
       else
          iStart = i ! first element of diffusion array
          exit
       endif
    enddo
    
 
    converged = .false.

!    do i = 1, iStart-1
       octVec = VECTOR(xPos, 0., zArray(istart-1))
       call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
            foundSubcell=subcell, rosselandKappa=kappa, grid=grid)
       if (thisOctal%undersampled(subcell)) then
          ok = .false.
          goto 666
       endif
!    enddo

    

    zarray(nz+1) = 0.
    t1 = treal(nreal)
    t0 = 0.99 * t1
    fac = log(t1/t0)/abs(zreal(nreal))
    do i = iStart,nz+1
       temperature(i) = t0 * exp(fac*abs(zArray(i)))
       dtdz(i) = t0 * fac *  exp(fac*abs(zArray(i))) / 1.e10
    enddo

    rho(nz+1) = rho(nz)

    
    octVec = VECTOR(xPos, 0., zArray(istart))
    call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
            foundSubcell=subcell,  grid=grid)

    if (thisOctal%nDiffusion(subcell)  == 0.) then
       temperature(istart:nz) = temperature(istart-1)
       goto 666
    endif

    inFlux = thisOctal%incidentflux(subcell)


    octVec = VECTOR(xPos, 0., zArray(istart-1))
    call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
         foundSubcell=subcell, rosselandKappa=kappa, grid=grid)

    firr(istart-1) = -inFlux

    nIter = 0
    do while(.not.converged)

       do i = 1, nz+1
          newTemperature(i) = temperature(i)
       enddo

       
       ! first calculate the opacities

       do i = iStart-1, nz+1
          octVec = VECTOR(xPos, 0., zArray(i))
          call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, &
               foundSubcell=subcell, rosselandKappa=kappa_ross(i), kappap=kappap(i), grid=grid, &
               atthistemperature=real(newtemperature(i)))
       enddo

       tau(istart-1) = 0.
       do i = istart, nz
          tau(i) = tau(i-1) + kappa_ross(i) *  rho(i) * abs(zArray(i)-zArray(i-1)) * 1.e10
       enddo
       do i = istart-1, nz
          firr(i) = -inFlux * exp(-tau(i))
          jirr(i) = influx * exp(-tau(i)) / fourPi
       enddo

       
       do i = iStart-1, nz
          gammairr(i) = fourPi * kappaP(i) * rho(i) * jirr(i)
       enddo
       

       fd(istart-1) = -firr(istart-1)
       frad(iStart-1) = firr(iStart-1) + fd(iStart-1)
       jd(istart-1) = (stefanBoltz * temperature(istart-1)**4 / pi) - jirr(istart-1)
       djdz(istart-1) = (3./fourPi) * kappa_ross(istart-1) * rho(istart-1) * fd(istart-1)
       lambdad(istart-1) = fourPi * kappap(istart-1)*rho(istart-1)*(stefanBoltz*newtemperature(istart-1)**4 /pi)
       gammad(istart-1) = fourPi * kappap(istart-1)*rho(istart-1)*jd(istart-1)
       dfddz(iStart-1) = lambdad(istart-1) - gammad(istart-1) -gammairr(istart-1)
       zArray(nz+1) = 0.
!       if (debugoutput) write(*,'(i3,8a12)') 0,"z","t","dfddz","fd","firr","jd","djdz","jdiff"
        i = istart-1
!        if (debugoutput) write(*,'(i3,1p,8e12.3)') &
!             i, zArray(i), newtemperature(i), dfddz(i), fd(i), firr(i), jd(i), djdz(i), jirr(i)
       do i = iStart, nz

          
          fd(i) = -firr(i)
          jd(i) = jd(i-1) + djdz(i-1)  * (zArray(i)-zArray(i-1)) * 1.e10             

          newTemperature(i) = (pi*(jd(i) + jirr(i))/stefanBoltz)
          if (newtemperature(i) > 0.) then
             newtemperature(i) = newTemperature(i)**0.25
             if (newTemperature(i) < 50.) then
                temperature(istart:nz) = temperature(istart-1)
                goto 666
             endif
          else
             temperature(istart:nz) = temperature(istart-1)
             goto 666
          endif

          lambdad(i) = fourPi * kappap(i)*rho(i)*(stefanBoltz*newtemperature(i)**4 /pi)
          gammad(i) = fourPi * kappap(i)*rho(i)*jd(i)
          dfddz(i) = lambdad(i) - gammad(i) - gammairr(i)
          djdz(i) = (3./fourPi) * kappa_ross(i) * rho(i) * fd(i)
!          if (debugoutput) &
!               write(*,'(i3,1p,8e12.3)') i, zArray(i), newtemperature(i), dfddz(i), fd(i), firr(i), jd(i), djdz(i), jirr(i)
       enddo

555 continue
       temperature(1:nz) = newTemperature(1:nz)

       do i = nz - 1, istart-1
          dtdz(i) = (temperature(i+1)-temperature(i))/(zArray(i+1)-Zarray(i)) / 1.e10
       enddo


       nIter = nIter + 1
       if (nIter > 3) then
          converged = .true.
!          stop
       endif
    enddo

666 continue

  end subroutine solveDiffusion

  subroutine throughoutMidplaneDiff(grid, epsoverdt, debugoutput)

    type(GRIDTYPE) :: grid
    real :: zAxis(10000), rho(10000), temperature(10000), subcellsize(10000)
    integer :: nz
    real(oct) :: epsOverdt
    real :: xpos, ypos, radius, drho, smallestSubcell
    logical :: diffApprox(10000)
    logical :: converged, debugoutput
    real :: xAxis(100000)
    integer :: nx, i, j
    real :: flux
    logical :: ok


    call calcIncidentFlux(grid, epsoverdt)

    call copytemperaturetoChiline(grid%octreeRoot)


    nx = 0
    call getxValuesdiff(grid%octreeRoot,nx,xAxis)
    call stripSimilarValues(xAxis,nx,real(1.e-5*grid%halfSmallestSubcell))
    xAxis(1:nx) = xAxis(1:nx) + 1.e-5*grid%halfSmallestSubcell

    xPos = grid%halfSmallestSubcell - grid%octreeRoot%subcellSize
    smallestSubcell = 2. * grid%halfSmallestSubcell
    
    xPos = grid%halfSmallestSubcell
    yPos = 0.

!    call locate(xaxis, nx, 14.*rsol/1.e10, j)
!    write(*,*) "Solving from x=2*rinner only!!!!!!!!!!!!!!!!!!"
    j = 1

    do i = j, nx
       radius = xAxis(i)
       xPos = xAxis(i)

       call getTemperatureDensityRunDiff(grid, zAxis, subcellsize, rho, temperature, diffApprox, xPos, yPos, nz, -1.)

       if (nz > 1) then
          call solveDiffusion(grid, zAxis, xPos, temperature, rho, diffapprox, nz, ok, debugoutput)

          if (ok) then
             call putTemperatureRunDiff(grid, zAxis, temperature, nz, xPos, yPos, -1.)
          endif
       endif

       call getTemperatureDensityRunDiff(grid, zAxis, subcellsize, rho, temperature, diffApprox, xPos, yPos, nz, +1.)

       if (nz > 1) then
          call solveDiffusion(grid, zAxis, xPos, temperature, rho, diffapprox, nz, ok, debugoutput)

          if (ok) then
             call putTemperatureRunDiff(grid, zAxis, temperature, nz, xPos, yPos, +1.)
          endif
       endif


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
             zAxis(nz) = abs(temp%z)
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
       currentPos = OCTALVECTOR(xPos, yPos, -1.*direction*zAxis(i))
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
       rosselandOpticalDepth(j) = rosselandOpticalDepth(j-1) + &
            kappa * thisOctal%rho(subcell) *  subcellsize(j)*1.e10
       if (rosselandOpticalDepth(j) > diffDepth) then
          xStart = xAxis(j)
          exit
       endif
    enddo
    

    nx = 0
    call getxValuesdiff(grid%octreeRoot,nx,xAxis)
    call stripSimilarValues(xAxis,nx,real(1.e-5*grid%halfSmallestSubcell))
    xAxis(1:nx) = xAxis(1:nx) + 1.e-5*grid%halfSmallestSubcell


    call locate(xAxis, nx, xStart, iBoundary)


    do i = iBoundary, nx
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

             if (i == iBoundary) then
                octVec = OCTALVECTOR(xAxis(iBoundary), 0., zAxis(j))
                call amrGridValues(grid%octreeRoot, octVec, foundOctal=boundaryOctal, &
                     foundSubcell=boundarySubcell)
                boundaryOctal%leftHandDiffusionBoundary(boundarySubcell) = .true.
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
          rosselandOpticalDepth(j) = rosselandOpticalDepth(j-1) + &
               kappa * thisOctal%rho(subcell) * subcellsize(j)*1.e10
          if (rosselandOpticalDepth(j) > diffDepth) then
             thisOctal%diffusionApprox(subcell) = .true.

             if (i == iBoundary) then
                octVec = OCTALVECTOR(xAxis(iBoundary), 0., zAxis(j))
                call amrGridValues(grid%octreeRoot, octVec, foundOctal=boundaryOctal, &
                     foundSubcell=boundarySubcell)
                boundaryOctal%leftHandDiffusionBoundary(boundarySubcell) = .true.
             endif

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


  subroutine diffPhotonPosition(grid, thisOctal, subcell, rVec)
    type(GRIDTYPE) :: grid
    integer :: subcell, tempsubcell
    type(OCTAL), pointer :: thisOctal, tempOctal, oldOctal
    type(OCTALVECTOR) :: rVec, cVec
    real(double) :: phi
    logical :: found
     
    rVec = subcellCentre(thisOctal, subcell)
    phi = atan2(rvec%y, rvec%z)
    cVec = rVec
    oldOctal => thisOctal
    tempOctal => thisOctal
    found = .false.

    if (thisOctal%leftHandDiffusionBoundary(subcell)) then
       rVec = rVec - xHatOctal * (0.5d0*tempOctal%subcellSize+grid%halfsmallestsubcell*0.01d0)
    else
       do while (.not.found)
          
          if (cVec%z > 0.d0) then
             rVec = rVec + zHatOctal * (0.5d0*tempOctal%subcellSize+grid%halfsmallestsubcell*0.01d0)
          else
             rVec = rVec - zHatOctal * (0.5d0*tempOctal%subcellSize+grid%halfsmallestsubcell*0.01d0)
          endif
          call amrGridValues(grid%octreeRoot, rVec, startOctal=oldOctal, &
               foundOctal=tempOctal, foundSubcell=tempsubcell)
          if (.not.tempOctal%diffusionApprox(tempSubcell)) then
             found = .true.
          else
             oldOctal => tempOctal
             rVec = subcellCentre(tempOctal, tempSubcell)
          endif
          
       end do
    endif
    rVec = rotateZ(rVec, phi)
  end subroutine diffPhotonPosition
       






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
    integer :: tnz
    real(double) :: t(*), z(*), chisq, apoly(4)
    integer :: nterms
    real(double) :: revT(1000), revZ(1000)
    integer :: i

    real :: xp(1000), yp(1000)
    ! fit a polynomial

    do i =1, nz
       revZ(i) = z(nz-i+1)
       revT(i) = t(nz-i+1)
    enddo

    tnz = min(4, nz)

    nTerms = 2

    call polfit(revz(2:tnz), revt(2:tnz), revt(2:tnz), tnz-1 ,nterms , -1, apoly, chisq)

!    getTempGradient = apoly(2) + 2.*revZ(1)*apoly(3)

    getTempGradient = apoly(2) * revZ(1)

    getTempGradient = getTempGradient / 1.e10 ! [K cm-1]

!    xp(1:nz) = revz(1:nz)
!    yp(1:nz) = revT(1:nz)
!    call pgbegin(0,"/xs",1,1)
!    call pgask(.true.)
!    call pgenv(xp(1),xp(nz), 6.,yp(nz)*1.1,0.,0.)
!    call pgbin(nz, xp, yp,.true.)
!    do i = 1, nz
!       yp(i) = apoly(1) + apoly(2)*xp(i)
!    enddo
!    call pgsci(2)
!    call pgline(tnz,xp,yp)
!    read(*,*) i
!    call pgend


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

  recursive subroutine copyTemperaturetochiline(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call copytemperaturetochiline(child)
                exit
             end if
          end do
       else
          if (thisOctal%diffusionApprox(subcell)) then
             thisOctal%chiline(subcell) = thisOctal%temperature(subcell)
          endif
       endif
    enddo
  end subroutine copyTemperaturetochiline

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


    kappap = 0.
    norm = 0.
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

  subroutine calcIncidentFlux(grid, epsoverdt)
    type(GRIDTYPE) :: grid
    type(OCTALVECTOR) :: octVec, cVec
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(oct) :: epsoverdt
    integer :: nx
    real :: xAxis(100000), zAxis(10000), rho(10000), temperature(10000)
    real :: subcellSize(10000)
    logical :: diffApprox(100000)
    real :: xpos, ypos
    integer :: i, j, nz
    real :: area
    integer :: iBoundary
    nx = 0
    call getxValuesdiff(grid%octreeRoot,nx,xAxis)
    call stripSimilarValues(xAxis,nx,real(1.e-5*grid%halfSmallestSubcell))
    xAxis(1:nx) = xAxis(1:nx) + 1.e-5*grid%halfSmallestSubcell

    iBoundary = 0
    do i = 1, nx
       xpos = xAxis(i)
       ypos = 0.
       call getTemperatureDensityRunDiff(grid, zAxis, subcellsize, rho, temperature, diffApprox, xPos, yPos, nz, -1.)
       do j = 1, nz
          if (diffApprox(j)) then
             iBoundary = j 
             exit
          endif
       end do
       if (iBoundary /= 0) then
          octVec=OCTALVECTOR(xpos, ypos, zAxis(iBoundary))
          call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, foundSubcell=Subcell, grid=grid)
          cVec = subcellCentre(thisOctal, subcell)
          area = pi*((cVec%x+thisOctal%subcellSize)**2-(cVec%x-thisOctal%subcellSize)**2) * 1.e20
          thisOctal%incidentflux(subcell) = thisOctal%nDiffusion(subcell)*epsOverdT / area 
       endif
       
       iBoundary = 0
       call getTemperatureDensityRunDiff(grid, zAxis, subcellsize, rho, temperature, diffApprox, xPos, yPos, nz, +1.)
       do j = 1, nz
          if (diffApprox(j)) then
             iBoundary = j 
             exit
          endif
       end do
       if (iBoundary /= 0) then
          octVec=OCTALVECTOR(xpos, ypos, -1.*zAxis(iBoundary))
          call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, foundSubcell=Subcell, grid=grid)
          cVec = subcellCentre(thisOctal, subcell)
          area = pi*((cVec%x+thisOctal%subcellSize/2.)**2-(cVec%x-thisOctal%subcellSize/2.)**2) * 1.e20
          thisOctal%incidentflux(subcell) = thisOctal%nDiffusion(subcell)*epsOverdT / area 
       endif
    enddo
  end subroutine calcIncidentFlux

  real function returnFluxDiff(thisTemp, ival, jval, fd, jd, gammadiff, &
       kappa_ross, kappap, rho, zAxis, temperature, nz, &
       nextJd, nextFd, nextfDiff)
       real :: thisTemp
       integer :: ival, jval
       real(double) :: fd(*), gammadiff(*)
       real :: jd(*), kappa_ross(*), kappap(*)
       real :: rho(*), zAxis(*), temperature(*)
       integer :: nz
       real :: djdz, dfddz, dtdz
       real :: nextjd, nextfd, nextfdiff


       dfddz = fourPi * kappap(ival) * rho(ival) * (stefanBoltz*temperature(ival)**4/pi - jd(ival))
       nextFd = fd(ival) + dfddz * (zAxis(jval)-zAxis(ival)) * 1.e10
       dtdz = (thisTemp-temperature(ival))/(zAxis(jval)-zAxis(ival)) / 1.e10
       nextFdiff = (-16.*stefanBoltz*thisTemp**3)/(3.*kappa_ross(jval)*rho(jval)) * dtdz
       djdz = -3./fourPi * kappa_ross(ival) * rho(ival) * fd(ival)
       nextJd = jd(ival) +  djdz * (zAxis(jval)-zAxis(ival)) * 1.e10

       returnFluxdiff = nextFdiff + nextFd
       write(*,*) thisTemp,temperature(ival),nextfd,nextfDiff,returnFluxdiff
  end function returnFluxDiff
end module diffusion_mod






