! this module solves the vertical density structure of a disc
! by repeated calls to the lucy algorithm and iteratively solving the
! equation of hydrostatic equilibrium

! written by tjh

! v1.0 on 16/3/04

module disc_hydro_mod

  use kind_mod
  use constants_mod
  use vector_mod
  use messages_mod
  use input_variables
  use gridtype_mod, only: GRIDTYPE
  use octal_mod, only: OCTAL, subcellCentre
  use lucy_mod, only: lucyRadiativeEquilibriumAMR

  implicit none

  public

contains
  
  subroutine solveHydro(temperature, zAxis, subcellsize, rho, nz, &
       radius, mStar, sigma0,  converged, drho)

    use density_mod, only: fractgap, fractgap2 
    use input_variables, only: alphaDisc, betaDisc, rinner, router, geometry, planetgap, rho0
    integer, intent(in) :: nz ! number of vertical grid points
    real(double),intent(in) :: subcellSize(1:nz)  ! size of this subcell
    real(single),intent(in) :: temperature(1:nz)  ! temperature
    real(double),intent(inout) :: rho(1:nz)         ! density
    real(double),intent(in) :: zAxis(1:nz)       ! z grid
    real :: mStar                     ! stellar mass, disk mass and radius
    real :: radius, radiusAU ! radius at disc midplane
    real :: sigma0                                   ! surface density at rinner
    integer :: i
    real ::  smoothscalelength
    logical :: converged
    real :: mu ! mean molecular weight
    real, allocatable :: oldRho(:)
    real :: dz
    real :: fac, fac2
    real(double) :: scalefac
    real(double), allocatable :: lnrho(:), doubleRho(:)
    real :: dtdz ! temperature gradient
    real(double) :: sigma
    real :: drho

    mu = 2.3 ! mean mole mass assuming pure H_2/He
    radiusAU = (radius * 1.e10) / autocm


    allocate(oldRho(1:nz), lnrho(1:nz),doubleRho(1:nz))
    oldRho(1:nz) = rho(1:nz)

    converged = .false.

    lnrho(1) = 0.
    do i = 2, nz
       dz = (zAxis(i)-zAxis(i-1))
       dTdz = (temperature(i)-temperature(i-1))/dz
       fac = (mu * mHydrogen * bigG * mStar * zAxis(i)/ (kErg * radius**3))/1.e30
       fac2 = -1.d0*(dTdz + fac)/temperature(i)
       lnrho(i) = lnrho(i-1) + fac2 * dz
    enddo

    ! now integrate over the column to get the surface density
    ! 

    doublerho(1:nz) = exp(lnrho(1:nz))

    sigma = 0.
    do i = 1,nz
       sigma = sigma + doublerho(i) * dble(subcellSize(i)) * 1.d10
    enddo

    ! rescale the density to retain the 1/r surface density fall-off

    scalefac = 0.5 * (sigma0*(autocm/(radius*1.e10))) / sigma ! factor of 0.5 cos only integrating over half disc vertically
    
    fac = 1.

    if ((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry == "warpeddisc")) then
       if (radius < rInner) then
          fac = ((rInner - radius)/(0.01*rInner))**2
          fac = exp(-fac)
       endif
    endif


    scalefac = 0.5 * (sigma0 * fac * (radius/rinner)**(betaDisc-alphaDisc)) / sigma


    if (geometry == "ppdisk") then
       smoothScaleLength = height * rHeight * (rSmooth/rHeight)**flaringPower
       scalefac = fractgap(dble(radiusAU)) / (1. + 81.d0**(dble(rSmooth-radiusAU)/smoothScaleLength)) &
            * (rHeight/radiusAU)**sigmaPower
       scalefac = scalefac * rho0 * Msol / autocm**2

!       scalefac in g/cm^3
       write (137,*) radiusAU, scalefac

       scalefac = scalefac * 0.5 / sigma
    endif


    if (geometry == "planetgap") then

       fac =  1.d0-min(dble(radius - rInner)/(0.02d0*rinner),1.d0)
       fac = exp(-fac*10.d0)
       scalefac = sigma0 * (radius/rCore)**(betaDisc-alphaDisc) * fac 
       if (planetGap) then
          scalefac = scalefac * fractGap2(dble(radiusAU))
       endif

!       scalefac in g/cm^3
       write (137,*) radiusAU, scalefac

       scalefac = scalefac * 0.5 / sigma

    endif

    if (geometry == "warpeddisc") then
       fac =  1.d0-min(dble(radius - rInner)/(0.05d0*rinner),1.d0)
       fac = exp(-fac*10.d0)
       scalefac = sigma0 * (radius/rOuter)**(betaDisc-alphaDisc) * fac

       write (137,*) radiusAU, scalefac

       scalefac = scalefac * 0.5 / sigma
    end if


!    scalefac = fractgap(dble(radiusAU)) / (sqrt(radiusAU) * (1 + exp(-100.*(radiusAU-rSmooth))) * sigma)
!    scalefac = scalefac * 0.5 * rho0 * Msol / autocm**2

!    write(*,*) radius, sigma0, sigma, &
!	0.5 * (sigma0 * (radius*1.e10/autoCm)**(betaDisc-alphaDisc)), scalefac

    rho(1:nz) = doublerho(1:nz) * scalefac

    ! now check for convergence. maximum fractional change in
    ! density should be less than 1%.

    do i = 1, nz
       if (oldRho(i) /= 0.) then
          drho = max(drho, abs(real(rho(i))-oldrho(i)))
       endif
    enddo

    deallocate(oldRho, doublerho)
  end subroutine solveHydro


  subroutine getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, xPos, yPos, nz, direction)
    use amr_mod, only: amrGridValues
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer, intent(out) :: nz
    real(double) :: rho(:)
    real :: temperature(:)
    real(double) :: zAxis(:), subcellsize(:)
    real :: xPos, yPos
    integer :: subcell
    real(double) :: rhotemp
    real :: temptemp
    real :: direction
    type(VECTOR) :: currentPos, temp
    real :: halfSmallestSubcell

    nz = 0
    halfSmallestSubcell = grid%halfSmallestSubcell

    currentPos = VECTOR(xPos, yPos, direction*halfSmallestSubcell)

    do while(abs(currentPos%z) < grid%ocTreeRoot%subcellsize)
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp, temperature=temptemp)
       thisOctal%chiLine(subcell) = 1.e-30
!       if (thisOctal%inFlow(subcell)) then
          nz = nz + 1
          temperature(nz) = temptemp
          rho(nz) = rhotemp
          temp = subCellCentre(thisOctal, subcell)
          zAxis(nz) = temp%z
          subcellsize(nz) = thisOctal%subcellsize
!       endif
          currentPos = VECTOR(xPos, yPos, zAxis(nz)+0.5*direction*thisOctal%subcellsize+direction*halfSmallestSubcell)
!       else
!          currentPos = VECTOR(xPos, yPos, grid%octreeRoot%subcellsize+halfSmallestSubcell)
!       endif
    end do
    zAxis(1:nz) = abs(zAxis(1:nz)) * 1.d10  ! convert to cm
  end subroutine getTemperatureDensityRun

  subroutine putDensityRun(grid, zAxis, rho, nz, xPos, yPos, direction)
    use amr_mod, only: amrGridValues
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer :: nz
    real :: xPos, yPos
    real(double) :: zAxis(:)
    real(double) :: rho(:)
    integer :: i, subcell
    real :: direction
    type(VECTOR) :: currentPos
    
    do i = 1, nz
       currentPos = VECTOR(xPos, yPos, direction*zAxis(i)/1.d10)
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell)
! original code moved to realPutDensity
!          if (thisOctal%rho(subcell) > 1.e-30) then
!             thisOctal%rho(subcell) = rho(i)
!             thisOctal%inFlow(subcell) = .true. 
!          else
!             thisOctal%inFlow(subcell) = .false.
!          endif
          thisOctal%chiLine(subcell) = rho(i)
    enddo
  end subroutine putDensityRun

  recursive subroutine realPutDensity(grid, thisOctal)

    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
!    real(double) :: deltaRho
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call realPutDensity(grid, child)
                exit
             end if
          end do
       else
          if (thisOctal%chiLine(subcell) > 1.e-30) then
!             deltaRho = thisOctal%chiLine(subcell) - thisOctal%rho(subcell)
!             thisOctal%rho(subcell) = thisOctal%rho(subcell) + 0.5d0 *deltaRho
             thisOctal%rho(subcell) = thisOctal%chiLine(subcell)
!             thisOctal%inFlow(subcell) = .true. 
          else
!             thisOctal%rho(subcell) = thisOctal%chiLine(subcell)
             thisOctal%rho(subcell) = 1.e-30
!             thisOctal%inFlow(subcell) = .false.
          endif
       end if
    end do

  end subroutine realPutDensity

  subroutine throughoutMidpane(grid, mStar, sigma0,  drho)
    use amr_mod, only: getxValues
    use utils_mod, only: stripSimilarValues
    type(GRIDTYPE) :: grid
    integer, parameter :: maxvals = 100000
    real(double) :: zAxis(maxVals)
    real(double) :: rho(maxVals)
    real :: temperature(maxVals)
    real(double) :: subcellsize(maxvals)
    real :: mStar, sigma0
    integer :: nz
    real :: xpos, ypos, radius, drho, smallestSubcell
    logical :: converged 
    real(double) :: xAxis(maxVals*10)
    integer :: nx, i
    converged = .false.; nz = 0; temperature = 0.; zAxis = 0.; subcellSize = 0.

    nx = 0
    call getxValues(grid%octreeRoot,nx,xAxis)
    call stripSimilarValues(xAxis,nx,1.d-5*grid%halfSmallestSubcell)
    xAxis(1:nx) = xAxis(1:nx) + 1.d-5*grid%halfSmallestSubcell

    xPos = grid%halfSmallestSubcell - grid%octreeRoot%subcellSize
    smallestSubcell = 2. * grid%halfSmallestSubcell
    



  if (grid%octreeRoot%twod) then
    xPos = grid%halfSmallestSubcell
    yPos = 0.
!    do while (xPos < (2.* grid%octreeRoot%subcellSize))
!       radius = xPos
!

    do i = 1, nx
       radius = xAxis(i)
       xPos = xAxis(i)
       if ((radius > grid%rInner).and.(radius < grid%rOuter)) then
          call getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, xPos, yPos, nz, +1.)

          if (nz > 1) then
             call solveHydro(temperature, zAxis, subcellsize, rho, nz, radius, mStar, &
                  sigma0, converged, drho)                

             call putDensityRun(grid, zAxis, rho, nz, xPos, yPos, +1.)
          endif

          call getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, xPos, yPos, nz, -1.)

          if (nz > 1) then
             call solveHydro(temperature, zAxis, subcellsize, rho, nz, radius, mStar, sigma0,  converged, drho)
             
             call putDensityRun(grid, zAxis, rho, nz, xPos, yPos, -1.)
          endif
       endif
!       xPos = xPos + smallestSubcell
!    end do
    enddo
  else if (grid%octreeRoot%threed) then
    do while (xPos < grid%octreeRoot%subcellSize)
       yPos = grid%halfSmallestSubcell - grid%octreeRoot%subcellSize
       do while (yPos < grid%octreeRoot%subcellSize)
          radius = sqrt(xPos**2 + yPos**2)
          if ((radius > grid%rInner).and.(radius < grid%rOuter)) then
             call getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, xPos, yPos, nz, +1.)

             if (nz > 1) then
                call solveHydro(temperature, zAxis, subcellsize, rho, &
                     nz, radius, mStar, sigma0, converged, drho)                

                call putDensityRun(grid, zAxis, rho, nz, xPos, yPos, +1.)
             endif

             call getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, xPos, yPos, nz, -1.)

             if (nz > 1) then
                call solveHydro(temperature, zAxis, subcellsize, rho, &
                     nz, radius, mStar, sigma0, converged, drho)
             
                call putDensityRun(grid, zAxis, rho, nz, xPos, yPos, -1.)
             endif
          endif
          yPos = yPos + smallestSubcell
       enddo
       xPos = xPos + smallestSubcell
    enddo
   end if
  end subroutine throughoutMidpane

  subroutine verticalHydrostatic(grid, mStar, sigma0,  miePhase, nDustType, nMuMie, nLambda, lamArray, &
       source, nSource, nLucy, massEnvelope)
    use input_variables, only : variableDustsublimation, rGap
    use messages_mod, only: myRankIsZero
    use parallel_mod, only: torus_mpi_barrier
    use source_mod, only: SOURCETYPE
    use phasematrix_mod, only: PHASEMATRIX
    use lucy_mod, only: refineDiscGrid, getSublimationRadius, putTau, unrefineBack
    use vtk_mod, only: writeVtkFile
    use amr_mod, only: myScaleSmooth, myTauSmooth, findTotalMass
    use utils_mod, only: locate

    type(GRIDTYPE) :: grid
    real :: mStar, sigma0
    integer :: nDustType
    integer :: iSmoothLam
    logical :: gridConverged
    logical :: converged
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nLucy
    integer :: nLambda, nMuMie
    type(PHASEMATRIX):: miePhase(:,:,:)
    real  :: lamArray(:)
    real :: massEnvelope
    real :: drho
    integer :: nIter, j
    real(double) totalMass
    real(double) :: rSub, betaEstimate, heightEstimate
    real :: temp
    real :: rGapCM
    integer :: nUnrefine
    real :: lamSmoothArray(5)
    character(len=80) :: message, plotfile
    
    lamSmoothArray = (/5500., 1.e4, 2.e4, 5.e4, 10.e4/)

    sigma0 = rho0 * (height*1.e10) * ((rinner*1.e10) / (100.d0*autocm))**betaDisc * sqrt(twopi)


    converged = .false.
    nIter = 1
    totalMass = 0.
    call findTotalMass(grid%octreeRoot, totalMass)


    if(myRankIsZero) &
         write(*,*) "Total disc mass: ",totalMass/msol," solar masses"



    do while(.not.converged)

       if(myRankIsZero) then
          write(*,*) "Starting iteration number",nIter
          write(*,*) "Calling the lucy algorithm to get temperature..."
       endif

       totalMass = 0.
       call findTotalMass(grid%octreeRoot, totalMass)

       call torus_mpi_barrier()

       call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, & 
            nLambda, lamArray, source, nSource, nLucy, massEnvelope, lucy_undersampled)


       if(myRankIsZero) &
            write(*,*) "Zeroing new density values..."
       call zeroChiline(grid%octreeRoot)       

       if(myRankIsZero) then
          if ((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry=="warpeddisc")) &
               open(unit=137, file='surface_density.dat', status='replace', form='formatted')
       endif

       if(myRankIsZero) &
            write(*,*) "Solving the vertical hydrostatic equilibrium..."
       drho = 0.
       call throughoutMidpane(grid, mStar, sigma0,  drho)

       if(myRankIsZero) then
          if ((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry=="warpeddisc")) close(137)
       endif

       if(myRankIsZero) &
            write(*,*) "Updating cell densities..."
       call realPutDensity(grid, grid%octreeRoot)

       if(myRankIsZero) &
            write(*,*) "Maximum absolute change in density: ",drho



       totalMass = 0.
       call findTotalMass(grid%octreeRoot, totalMass)
       if (myRankisZero) write(*,*) "Total disc mass: ",totalMass/msol," solar masses"

       call getBetaValue(grid, betaEstimate, heightEstimate)
       call getSublimationRadius(grid, rSub)
       write(message, '(a, f7.3,a )') "After hydro adjustment: Dust Sublimation radius is: ",(1.d10*rSub/rSol), " solar radii"
       call writeInfo(message, FORINFO)
       write(message, '(a, f7.3,a )') "After hydro adjustment: Dust Sublimation radius is: ",(rSub/rCore), " core radii"
       call writeInfo(message, FORINFO)

       if (writeoutput) write(*,*) "Refining on new estimates of disc parameters"
       gridconverged = .false.
       do while(.not.gridconverged)
          gridconverged = .true.
          nUnrefine = 0
          call refineDiscGrid(grid%octreeRoot, grid, betaEstimate, heightEstimate, rSub, gridconverged, &
               inheritprops = .false., interpProps = .true.)
       end do
       if (writeoutput) then
          write(*,*) "done."
       endif


       

       if (writeoutput) write(*,*) "Unrefining back to current gridding..."
       gridconverged = .false.
       do while(.not.gridconverged)
          gridconverged = .true.
          nUnrefine = 0
          call unrefineBack(grid%octreeRoot, grid, betaEstimate, heightEstimate, rSub, nUnrefine, gridconverged)
          if (writeoutput) write(*,*) "Unrefined ",nUnrefine, " cells on this pass"
       end do

       if (writeoutput) then
          write(*,*) "done."
       endif

       if ((.not.variableDustSublimation).and.grid%octreeRoot%twoD) then
          call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)
          
          call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
          do j = iSmoothLam, nLambda, 2
             write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
             call writeInfo(message, TRIVIAL)
             do
                gridConverged = .true.
                call putTau(grid, grid%lamArray(j))
                call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                     inheritProps = .false., interpProps = .true., photosphereSplit = .true.)
                
                if (gridConverged) exit
             end do
          enddo
       endif

       call writeInfo("Smoothing adaptive grid structure (again)...", TRIVIAL)
       do
          gridConverged = .true.
          call myScaleSmooth(smoothFactor,  grid, &
               gridConverged,  inheritProps = .false., interpProps = .true.)
          if (gridConverged) exit
       end do
       call writeInfo("...grid smoothing complete", TRIVIAL)

       call writeVtkFile(grid, "aftersmooth.vtk", &
            valueTypeString=(/"rho        ", "temperature", "tau        ", "crossings  ", "etacont    " , &
            "dust1      ", "deltaT     ", "etaline    "/))
       

       if(myRankIsZero) then

          ! chris
          if ((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry=="warpeddisc")) then
             if (geometry.ne."warpeddisc") then
                rGapCM = rGap * autocm / 1.d10
             end if

          end if
       endif

       ! chris (26/05/04)
       ! Smooth the grid with respect to optical depth, if requested


       if (.not.variableDustSublimation) then

          call locate(grid%lamArray, nLambda, lambdaSmooth,ismoothlam)
          call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
          do j = iSmoothLam, nLambda
             write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
             call writeInfo(message, TRIVIAL)
             do
                call putTau(grid, grid%lamArray(j))
                gridConverged = .true.
                call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                     inheritProps = .false., interpProps = .true., photospheresplit = .true.)

                if (gridConverged) exit
             end do
          enddo
          call writeInfo("...grid smoothing complete", TRIVIAL)

          call writeInfo("Smoothing adaptive grid structure (again)...", TRIVIAL)
          do
             gridConverged = .true.
             call myScaleSmooth(smoothFactor, grid, &
                  gridConverged,  inheritProps = .false., interpProps = .true.)
             if (gridConverged) exit
          end do
          call writeInfo("...grid smoothing complete", TRIVIAL)

          if(myRankIsZero) &
               write(*,*) "...grid smoothing complete"

          if(myRankIsZero) &
               write(*,*) "Total disc mass: ",totalMass/msol," solar masses"

       end if
       
       write(plotfile,'(a,i3.3,a)') "hydrostep_",nIter,".vtk"
       call writeVtkFile(grid, plotfile, &
            valueTypeString=(/"rho        ", "temperature", "tau        ", "crossings  ", "etacont    " , &
            "dust1      ", "deltaT     ", "etaline    "/))


    nIter = nIter + 1

    if (nIter > nHydro) then
       if(myRankIsZero) &
            write(*,*) "Maximum number of iterations exceeded. Aborting."
       converged = .true.
    else
       temp = 20.
       call setTemperature(grid%octreeRoot, temp)
    endif

 enddo


 ! final call to sort out temperature

 temp = 20.
 call setTemperature(grid%octreeRoot, temp)
 call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, & 
      nLambda, lamArray, source, nSource, nLucy, massEnvelope, lucy_undersampled,  finalpass = .true.)




 ! chris
 if(myRankIsZero) then
    if ((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry=="warpeddisc")) then

    end if
 endif

end subroutine verticalHydrostatic


  recursive subroutine setTemperature(thisOctal, temperature)

    use input_variables, only : rinner, router
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real :: temperature
    type(VECTOR) :: rvec
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setTemperature(child, temperature)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          if ((rVec%x > rinner).and.(rVec%x< rOuter)) then
             thisOctal%temperature(subcell) = temperature
          endif
       end if
    end do

  end subroutine setTemperature

  recursive subroutine zeroChiline(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroChiline(child)
                exit
             end if
          end do
       else
          thisOctal%chiline(subcell) = tiny(thisOctal%chiline(subcell))
       endif
    enddo
  end subroutine zerochiline


  subroutine getBetaValue(grid, beta, heightat100AU)
    use input_variables, only : rinner, router
    use utils_mod, only: linfit 

    type(GRIDTYPE) :: grid
    real(double), intent(out) :: beta, heightat100AU
    real(double) :: height(100), r(100)
    real(double) :: a, sigmaa, b, sigmab, rCoeff, rhoMid, rhoScale
    real(double), allocatable :: zAxis(:), rho(:), subcellsize(:)
    real, allocatable :: temperature(:)
    integer :: nz
    integer :: i, j
    integer, parameter :: m =  100000

    allocate(zAxis(m), rho(m), temperature(m), subcellsize(m))
    do i = 1, 100
       r(i) = log10(rInner*1.1) + (log10(rOuter*0.99) - log10(Rinner*1.1))*dble(i-1)/99.d0
       r(i) = 10.d0**r(i)
       call getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, real(r(i)), 0., nz, 1.)
       
       rhoMid = rho(1)
       rhoScale = rhoMid * exp(-0.5d0)
       do j = 1, nz-1
          if ((rhoScale < rho(j)).and.(rhoScale >= rho(j+1))) then
             height(i) = zAxis(j+1) + (zAxis(j)-zAxis(j+1)) * (rhoScale - rho(j+1))/(rho(j) - rho(j+1))
             exit
          endif
       enddo
    enddo
    r = log10(r)
    height = log10(height)-10.d0 ! factor 10^10

    
    call LINFIT(r,height,height,100, 0, A, SIGMAA, B, SIGMAB, Rcoeff)
    heightAt100AU = (10.d0**a) * (100.d0*autocm/1.d10)**b
    beta = b

    if (writeoutput) then
       write(*,*) "Disc beta value: ",beta
       write(*,*) "Height at 100 AU: ",1.d10*heightAt100AU/autocm
    endif
  end subroutine getBetaValue

  

end module disc_hydro_mod

