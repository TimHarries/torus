
! hydrodynamics module added by TJH on 18th June 2007

module hydrodynamics_mod

  use kind_mod
  use constants_mod
  use amr_mod

  implicit none


contains

  recursive subroutine fluxLimiter(thisOctal, limiterType)
    character(len=*) :: limiterType
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: a, b, dq
    character(len=80) :: message
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fluxLimiter(child, limiterType)
                exit
             end if
          end do
       else
          dq = thisOctal%q_i(subcell) - thisOctal%q_i_minus_1(subcell)
          if (abs(dq) > 0.d0) then
             if (thisOctal%u_interface(subcell) .ge. 0.d0) then
                thisOctal%rLimit(subcell) = (thisOctal%q_i_minus_1(subcell) - thisOctal%q_i_minus_2(subcell)) / dq
             else
                thisOctal%rLimit(subcell) = (thisOctal%q_i_plus_1(subcell) - thisOctal%q_i(subcell)) / dq
             endif
          else
             thisOctal%rLimit(subcell) = 0.d0
          endif
          select case (limiterType)
             case("superbee")
               a = min(1.d0, 2.d0*thisOctal%rLimit(subcell))
               b = min(2.d0, thisOctal%rLimit(subcell))
               thisOctal%phiLimit(subcell) = max(0.d0, a, b)
             case DEFAULT
                write(message,'(a,a)') "Flux limiter not recognised: ", trim(limiterType)
                call writeFatal(message)
                stop
           end select
!           write(*,*) "philimit",thisOctal%phiLimit(subcell)
           thisOctal%philimit(subcell) = 1.d0
       endif
    enddo
  end subroutine fluxLimiter

  recursive subroutine constructFlux(thisOctal, dt)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call constructFlux(child, dt)
                exit
             end if
          end do
       else
          if (thisOctal%u_interface(subcell).ge.0.d0) then
             thisOctal%flux_i(subcell) = thisOctal%u_interface(subcell) * thisOctal%q_i_minus_1(subcell)
          else
             thisOctal%flux_i(subcell) = thisOctal%u_interface(subcell) * thisOctal%q_i(subcell)
          endif
          thisOctal%flux_i(subcell) = thisOctal%flux_i(subcell) + &
               0.5d0 * abs(thisOctal%u_interface(subcell)) * &
               (1.d0 - abs(thisOctal%u_interface(subcell) * dt / &
               (thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell)))) * &
               thisOctal%phiLimit(subcell) * (thisOctal%q_i(subcell) - thisOctal%q_i_minus_1(subcell))
          write(*,*) "construct flux",thisOctal%flux_i(subcell)
       endif
    enddo
  end subroutine constructFlux

  recursive subroutine updateCellQ(thisOctal, dt)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call updateCellQ(child, dt)
                exit
             end if
          end do
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%q_i(subcell) = thisOctal%q_i(subcell) - dt * &
                  (thisOctal%flux_i_plus_1(subcell) - thisOctal%flux_i(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i(subcell))
             write(*,*) "in updatecellq", thisOctal%flux_i_plus_1(subcell) - thisOCtal%flux_i(subcell)
          endif
       endif
    enddo
  end subroutine updateCellQ

  recursive subroutine setupX(thisOctal, grid, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupX(child, grid, direction)
                exit
             end if
          end do
       else
          thisOctal%x_i(subcell) = subcellCentre(thisOctal, subcell) .dot. direction
       endif
    enddo
  end subroutine setupX

  recursive subroutine setupQX(thisOctal, grid, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupQX(child, grid, direction)
                exit
             end if
          end do
       else
          thisOctal%x_i_minus_1(subcell) = 0.d0
          thisOctal%x_i_plus_1(subcell) = 0.d0
          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%x_i_plus_1(subcell) = subcellCentre(neighbourOctal, neighbourSubcell) .dot. direction
             thisOctal%q_i_plus_1(subcell) = neighbourOctal%q_i(neighbourSubcell)
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%x_i_minus_1(subcell) = subcellCentre(neighbourOctal, neighbourSubcell) .dot. direction
             thisOctal%q_i_minus_1(subcell) = neighbourOctal%q_i(neighbourSubcell)
             
             locator = subcellCentre(neighbourOctal, neighboursubcell) - &
                  direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%q_i_minus_2(subcell) = neighbourOctal%q_i(neighbourSubcell)
          endif
       endif
    enddo
  end subroutine setupQX

  recursive subroutine setupUi(thisOctal, grid, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupUi(child, grid, direction)
                exit
             end if
          end do
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             rho_i_minus_1 = neighbourOctal%rho(neighbourSubcell)
             rhou_i_minus_1 = neighbourOctal%rhou(neighbourSubcell)
             thisOctal%u_interface(subcell) = 0.5d0 * &
                  (thisOctal%rhou(subcell)/thisOctal%rho(subcell) + rhou_i_minus_1/rho_i_minus_1)
             write(*,*) "setupi",thisOctal%u_interface(subcell)
          endif
       endif
    enddo
  end subroutine setupUI


  recursive subroutine setupFlux(thisOctal, grid, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupFlux(child, grid, direction)
                exit
             end if
          end do
       else

          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%x_i(subcell) = subcellCentre(thisOctal, subcell) .dot. direction
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%flux_i_plus_1(subcell) = neighbourOctal%flux_i(neighbourSubcell)
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%flux_i_minus_1(subcell) = neighbourOctal%flux_i(neighbourSubcell)
          endif
       endif
    enddo
  end subroutine setupFlux



   recursive subroutine setupPressure(thisOctal, gamma, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: gamma
    type(OCTALVECTOR) :: direction, locator
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupPressure(child, gamma, direction)
                exit
             end if
          end do
       else
          
          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%pressure_i_plus_1(subcell) = neighbourOctal%pressure_i(neighbourSubcell)
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%pressure_i_minus_1(subcell) = neighbourOctal%pressure_i(neighbourSubcell)
          endif

       endif
    enddo
  end subroutine setupPressure

   recursive subroutine computePressure(thisOctal, gamma, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: gamma
    type(OCTALVECTOR) :: direction, locator
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computePressure(child, gamma, direction)
                exit
             end if
          end do
       else
          
          thisOctal%pressure_i(subcell) = (gamma - 1.d0) * thisOctal%rho(subcell) * thisOctal%energy(subcell)

       endif
    enddo
  end subroutine computePressure

  recursive subroutine pressureForce(thisOctal, dt)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call pressureForce(child, dt)
                exit
             end if
          end do
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))
             if (isnan(thisOctal%rhou(subcell))) then
                write(*,*) "bug"
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureForce



  recursive subroutine copyRhoToQ(thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call copyRhoToQ(child)
                exit
             end if
          end do
       else
  
          thisOctal%q_i(subcell) = thisOctal%rho(subcell)
        
       endif
    enddo
  end subroutine copyRhoToQ

  recursive subroutine copyRhoUToQ(thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call copyRhoUToQ(child)
                exit
             end if
          end do
       else
  
          thisOctal%q_i(subcell) = thisOctal%rhoU(subcell)
        
       endif
    enddo
  end subroutine copyRhoUToQ

  recursive subroutine copyQtoRho(thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call copyQtoRho(child)
                exit
             end if
          end do
       else
  
          thisOctal%rho(subcell) = thisOctal%q_i(subcell)
        
       endif
    enddo
  end subroutine copyQtoRho

  recursive subroutine copyQtoRhoU(thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call copyQtoRhoU(child)
                exit
             end if
          end do
       else
  
          thisOctal%rhou(subcell) = thisOctal%q_i(subcell)
        
       endif
    enddo
  end subroutine copyQtoRhoU


  subroutine advectRho(grid, direction, dt)

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRho(grid%octreeRoot)

  end subroutine advectRho

  subroutine advectRhoU(grid, direction, dt)

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoUToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRhoU(grid%octreeRoot)

  end subroutine advectRhoU

  subroutine hydroStep(grid, gamma, dt, direction, boundaryCondition)
    type(GRIDTYPE) :: grid
    real(double) :: gamma, dt
    type(OCTALVECTOR) :: direction
    character(len=*) :: boundaryCondition

!    call imposeBoundaryConditions(grid, boundarycondition)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt)
    call advectRhoU(grid, direction, dt)
    call computePressure(grid%octreeRoot, gamma, direction)
    call setupPressure(grid%octreeRoot, gamma, direction)
    call  pressureForce(grid%octreeRoot, dt)
!    call imposeBoundaryConditions(grid, boundarycondition)
  end subroutine hydroStep


  recursive subroutine computeCourantTime(thisOctal, tc, gamma)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: tc, dx, cs, gamma
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computeCourantTime(child, tc, gamma)
                exit
             end if
          end do
       else
  
          if (.not.thisOctal%ghostCell(subcell)) then
             cs = sqrt(gamma*(gamma-1.d0)*thisOctal%energy(subcell))
             dx = abs(thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell))
             tc = min(tc, dx / (cs + abs(thisOctal%rhou(subcell)/thisOctal%rho(subcell)) ) )
          endif
        
       endif
    enddo
  end subroutine computeCourantTime

  subroutine doHydrodynamics(grid)
    type(gridtype) :: grid
    real(double) :: dt, tc, cfl, gamma, mu
    integer :: i, pgbegin
    type(OCTALVECTOR) :: direction

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    gamma = 5.d0 / 3.d0
    cfl = 0.5d0
    mu = 2.d0
    i = pgbegin(0,"/xs",1,1)
    call pgenv(0., 100., 0., 1.5, 0, 0)
!    call pgenv(0., 100., -1., 1., 0, 0)

    call setGhostCells(grid, direction)
    call calculateRhoU(grid%octreeRoot, direction)
    call calculateEnergy(grid%octreeRoot, gamma, mu)
    do while(.true.)
       tc = 1.d30
       call computeCourantTime(grid%octreeRoot, tc, gamma)
       dt = tc * cfl
       write(*,*) "dt",dt
       call hydroStep(grid, gamma, dt, direction, "periodic")
       call plotHydroResults(grid)
    enddo
    call pgend
  end subroutine doHydrodynamics

  recursive subroutine calculateEnergy(thisOctal, gamma, mu)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: gamma, mu
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateEnergy(child, gamma, mu)
                exit
             end if
          end do
       else
 
!          thisOctal%energy(subcell) = ((thisOctal%rho(subcell)/(mu*mhydrogen))*kerg*thisOctal%temperature(subcell)) / &
!               ((gamma-1.d0)*thisOctal%rho(subcell))
          thisOctal%energy(subcell) = 1.d0
        
       endif
    enddo
  end subroutine calculateEnergy

  recursive subroutine calculateRhoU(thisOctal, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(OCTALVECTOR) :: direction
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateRhoU(child, direction)
                exit
             end if
          end do
       else 
          thisOctal%rhoU(subcell) = thisOctal%rho(subcell) * &
               (cSpeed * (thisOctal%velocity(subcell).dot.direction))
       endif
    enddo
  end subroutine calculateRhoU

  subroutine plotHydroResults(grid)

    type(GRIDTYPE) :: grid
    real(double) :: x(1000), rho(1000), v(1000), rhou(1000)
    real :: xr(1000), yr(1000)
    integer :: n,i,pgbegin
    
    n = 0
    call getArray(grid%octreeRoot, x, rho, rhou, v, n)
    xr(1:n) = real(x(1:n))
    yr(1:n) = real(rho(1:n))
!    yr(1:n) = real(rhou(1:n))



!    write(*,*) "xr",xr(1:n)
!    write(*,*) "yr",yr(1:n)
    call pgline(n, xr, yr)
  end subroutine plotHydroResults

  recursive subroutine getArray(thisOctal, x, rho, rhou, v, n)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: x(:), rho(:), v(:), rhou(:)
    integer :: n
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getArray(child,  x, rho, rhou, v, n)
                exit
             end if
          end do
       else
          n = n + 1
          rho(n) = thisOctal%rho(subcell)
          rhou(n) = thisOctal%rhou(subcell)
          v(n) = thisOctal%rhou(subcell)/thisOctal%rho(subcell)
          x(n) = thisOctal%x_i(subcell)
      endif
    enddo
  end subroutine getArray


  subroutine setGhostCells(grid, direction)
    use input_variables, only : x1, x2
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, bOctal
    integer :: subcell, bSubcell
    type(OCTALVECTOR) :: locator, direction
    real(double) :: dx

    locator = OCTALVECTOR(dble(x1)+ 0.1d0*grid%halfSmallestSubcell, 0.d0, 0.d0)
    thisOctal => grid%octreeRoot
    dx = thisOctal%subcellSize
    call findSubcellLocal(locator, thisOctal, subcell)
    thisOctal%ghostCell(subcell) = .true.

    locator = OCTALVECTOR(dble(x1)+dx+0.1d0*grid%halfSmallestSubcell, 0.d0, 0.d0)
    thisOctal => grid%octreeRoot
    call findSubcellLocal(locator, thisOctal, subcell)
    thisOctal%ghostCell(subcell) = .true.

    locator = OCTALVECTOR(dble(x2) - 0.1d0*grid%halfSmallestSubcell, 0.d0, 0.d0)
    thisOctal => grid%octreeRoot
    dx = thisOctal%subcellSize
    call findSubcellLocal(locator, thisOctal, subcell)
    thisOctal%ghostCell(subcell) = .true.

    locator = OCTALVECTOR(dble(x2) - dx - 0.1d0*grid%halfSmallestSubcell, 0.d0, 0.d0)
    thisOctal => grid%octreeRoot
    dx = thisOctal%subcellSize
    call findSubcellLocal(locator, thisOctal, subcell)
    thisOctal%ghostCell(subcell) = .true.
    
  end subroutine setGhostCells

end module hydrodynamics_mod
