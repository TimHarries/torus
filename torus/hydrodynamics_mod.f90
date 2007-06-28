
! hydrodynamics module added by TJH on 18th June 2007

module hydrodynamics_mod

  use kind_mod
  use constants_mod
  use amr_mod
  use grid_mod

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
             thisOctal%rLimit(subcell) = 1.d0
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
!             thisOctal%q_i(subcell) = thisOctal%q_i(subcell) - dt * &
!                  (thisOctal%flux_i_plus_1(subcell) - thisOctal%flux_i(subcell)) / &
!                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i(subcell))
             thisOctal%q_i(subcell) = thisOctal%q_i(subcell) - dt * &
                  (thisOctal%flux_i_plus_1(subcell) - thisOctal%flux_i(subcell)) / &
                  (thisOctal%subcellSize)

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
    real(double) :: rhou, rhov, rho, q
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

             if (thisOctal%nDepth >= neighbourOctal%nDepth) then ! same level
                thisOctal%q_i_plus_1(subcell) = neighbourOctal%q_i(neighbourSubcell)
             else
                call averageValue(thisOctal, subcell, neighbourOctal, q, rhou, rhov, rho)
                thisOctal%q_i_plus_1(subcell) = q
             endif
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%x_i_minus_1(subcell) = subcellCentre(neighbourOctal, neighbourSubcell) .dot. direction


             if (thisOctal%nDepth >= neighbourOctal%nDepth) then ! same level
                thisOctal%q_i_minus_1(subcell) = neighbourOctal%q_i(neighbourSubcell)
             else
                call averageValue(thisOctal, subcell, neighbourOctal, q, rhou, rhov, rho)
                thisOctal%q_i_minus_1(subcell) = q
             endif

             
             locator = subcellCentre(neighbourOctal, neighboursubcell) - &
                  direction * (neighbourOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             if (thisOctal%nDepth >= neighbourOctal%nDepth) then ! same level
                thisOctal%q_i_minus_2(subcell) = neighbourOctal%q_i(neighbourSubcell)
             else
                call averageValue(thisOctal, subcell, neighbourOctal, q, rhou, rhov, rho)
                thisOctal%q_i_minus_2(subcell) = q
             endif


          endif
       endif
    enddo
  end subroutine setupQX

  recursive subroutine setupUi(thisOctal, grid, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, q
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

             if (thisOctal%nDepth >= neighbourOctal%nDepth) then ! same level
                rho_i_minus_1 = neighbourOctal%rho(neighbourSubcell)
                rhou_i_minus_1 = neighbourOctal%rhou(neighbourSubcell)
                thisOctal%u_interface(subcell) = 0.5d0 * &
                     (thisOctal%rhou(subcell)/thisOctal%rho(subcell) + rhou_i_minus_1/rho_i_minus_1)
             else
                call averageValue(thisOctal, subcell, neighbourOctal, q, rhou, rhov, rho)
                thisOctal%q_i_minus_2(subcell) = q
                rho_i_minus_1 = rho
                rhou_i_minus_1 = rhou
                thisOctal%u_interface(subcell) = 0.5d0 * &
                     (thisOctal%rhou(subcell)/thisOctal%rho(subcell) + rhou_i_minus_1/rho_i_minus_1)
             endif

          endif
       endif
    enddo
  end subroutine setupUI

  recursive subroutine setupVi(thisOctal, grid, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, q
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupVi(child, grid, direction)
                exit
             end if
          end do
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             if (thisOctal%nDepth >= neighbourOctal%nDepth) then ! same level
                rho_i_minus_1 = neighbourOctal%rho(neighbourSubcell)
                rhou_i_minus_1 = neighbourOctal%rhov(neighbourSubcell)
                thisOctal%u_interface(subcell) = 0.5d0 * &
                     (thisOctal%rhov(subcell)/thisOctal%rho(subcell) + rhou_i_minus_1/rho_i_minus_1)
             else
                call averageValue(thisOctal, subcell, neighbourOctal, q, rhou, rhov, rho)
                thisOctal%q_i_minus_2(subcell) = q
                rho_i_minus_1 = rho
                rhou_i_minus_1 = rhov
                thisOctal%u_interface(subcell) = 0.5d0 * &
                     (thisOctal%rhov(subcell)/thisOctal%rho(subcell) + rhou_i_minus_1/rho_i_minus_1)
             endif
          endif
       endif
    enddo
  end subroutine setupVI


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


  recursive subroutine setupPressure(thisOctal, grid, direction)
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
                call setupPressure(child, grid, direction)
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

  recursive subroutine setupUpm(thisOctal, grid, direction)
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
                call setupUpm(child, grid, direction)
                exit
             end if
          end do
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%u_i_plus_1(subcell) = neighbourOctal%rhou(neighbourSubcell)/neighbourOctal%rho(neighbourSubcell)
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%u_i_minus_1(subcell) = neighbourOctal%rhou(neighbourSubcell)/neighbourOctal%rho(neighbourSubcell)

          endif
       endif
    enddo
  end subroutine setupUpm

  recursive subroutine setupVpm(thisOctal, grid, direction)
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
                call setupVpm(child, grid, direction)
                exit
             end if
          end do
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%u_i_plus_1(subcell) = neighbourOctal%rhov(neighbourSubcell)/neighbourOctal%rho(neighbourSubcell)
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             thisOctal%u_i_minus_1(subcell) = neighbourOctal%rhov(neighbourSubcell)/neighbourOctal%rho(neighbourSubcell)

          endif
       endif
    enddo
  end subroutine setupVpm

   recursive subroutine computePressureU(thisOctal, gamma, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: gamma, eThermal, eKinetic, eTot, u, bigGamma,eta
    type(OCTALVECTOR) :: direction, locator

    eta = 3.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computePressureU(child, gamma, direction)
                exit
             end if
          end do
       else
          
          u = thisOctal%rhou(subcell) / thisOctal%rho(subcell)
          eKinetic = u**2 / 2.d0
          eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
          eThermal = eTot - eKinetic
          thisOctal%pressure_i(subcell) = (gamma - 1.d0) * thisOctal%rho(subcell) * eThermal
!          thisOctal%pressure_i(subcell) = (gamma - 1.d0) * thisOctal%rho(subcell) * thisOctal%energy(subcell)

!          write(*,*) thisOctal%pressure_i(subcell), thisOctal%rho(subcell),ethermal, thisOCtal%energy(subcell)
          bigGamma = 0.d0
          if (.not.thisOctal%ghostCell(subcell)) then
             if (thisOctal%u_i_plus_1(subcell) .le. thisOctal%u_i_minus_1(subcell)) then
                bigGamma = 0.25d0 * eta**2 * (thisOctal%u_i_plus_1(subcell) - thisOctal%u_i_minus_1(subcell))**2 &
                     * thisOctal%rho(subcell)
             endif
          endif
          thisOctal%pressure_i(subcell) = thisOctal%pressure_i(subcell) + bigGamma
       endif
    enddo
  end subroutine computePressureU

   recursive subroutine computePressureV(thisOctal, gamma, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: gamma, eThermal, eKinetic, eTot, u, bigGamma,eta
    type(OCTALVECTOR) :: direction, locator

    eta = 3.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computePressureV(child, gamma, direction)
                exit
             end if
          end do
       else
          
          u = thisOctal%rhov(subcell) / thisOctal%rho(subcell)
          eKinetic = u**2 / 2.d0
          eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
          eThermal = eTot - eKinetic
          thisOctal%pressure_i(subcell) = (gamma - 1.d0) * thisOctal%rho(subcell) * eThermal
!          thisOctal%pressure_i(subcell) = (gamma - 1.d0) * thisOctal%rho(subcell) * thisOctal%energy(subcell)
!          write(*,*) thisOctal%pressure_i(subcell), thisOctal%rho(subcell),ethermal, thisOCtal%energy(subcell)
          bigGamma = 0.d0
          if (.not.thisOctal%ghostCell(subcell)) then
             if (thisOctal%u_i_plus_1(subcell) .le. thisOctal%u_i_minus_1(subcell)) then
                bigGamma = 0.25d0 * eta**2 * (thisOctal%u_i_plus_1(subcell) - thisOctal%u_i_minus_1(subcell))**2 &
                     * thisOctal%rho(subcell)
             endif
          endif
          thisOctal%pressure_i(subcell) = thisOctal%pressure_i(subcell) + bigGamma
       endif
    enddo
  end subroutine computePressureV

  recursive subroutine pressureForceU(thisOctal, dt)
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
                call pressureForceU(child, dt)
                exit
             end if
          end do
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))

             thisOctal%rhoE(subcell) = thisOctal%rhoE(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
                  thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))

             

             if (isnan(thisOctal%rhou(subcell))) then
                write(*,*) "bug"
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureForceU

  recursive subroutine pressureForceV(thisOctal, dt)
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
                call pressureForceV(child, dt)
                exit
             end if
          end do
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%rhov(subcell) = thisOctal%rhov(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))

             thisOctal%rhoE(subcell) = thisOctal%rhoE(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
                  thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))

             

             if (isnan(thisOctal%rhov(subcell))) then
                write(*,*) "bug"
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureForceV



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

  recursive subroutine copyRhoEToQ(thisOctal)
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
                call copyRhoEToQ(child)
                exit
             end if
          end do
       else
  
          thisOctal%q_i(subcell) = thisOctal%rhoE(subcell)
        
       endif
    enddo
  end subroutine copyRhoEToQ

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

  recursive subroutine copyRhoVToQ(thisOctal)
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
                call copyRhoVToQ(child)
                exit
             end if
          end do
       else
  
          thisOctal%q_i(subcell) = thisOctal%rhoV(subcell)
        
       endif
    enddo
  end subroutine copyRhoVToQ

  recursive subroutine copyQtoRhoV(thisOctal)
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
                call copyQtoRhoV(child)
                exit
             end if
          end do
       else
  
          thisOctal%rhoV(subcell) = thisOctal%q_i(subcell)
        
       endif
    enddo
  end subroutine copyQtoRhoV

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

  recursive subroutine copyQtoRhoE(thisOctal)
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
                call copyQtoRhoE(child)
                exit
             end if
          end do
       else
  
          thisOctal%rhoE(subcell) = thisOctal%q_i(subcell)
        
       endif
    enddo
  end subroutine copyQtoRhoE

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

  subroutine advectRhoE(grid, direction, dt)

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoEToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRhoE(grid%octreeRoot)

  end subroutine advectRhoE

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

  subroutine advectRhoV(grid, direction, dt)

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoVToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRhoV(grid%octreeRoot)

  end subroutine advectRhoV


  subroutine hydroStep(grid, gamma, dt, direction, boundaryCondition)
    type(GRIDTYPE) :: grid
    real(double) :: gamma, dt
    type(OCTALVECTOR) :: direction
    character(len=*) :: boundaryCondition

    call imposeBoundary(grid%octreeRoot, boundarycondition)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt)
    call advectRhoU(grid, direction, dt)
    call advectRhoE(grid, direction, dt)
    call setupUpm(grid%octreeRoot, grid, direction)
    call computePressureU(grid%octreeRoot, gamma, direction)
    call setupPressure(grid%octreeRoot, grid, direction)
    call pressureForceU(grid%octreeRoot, dt)
    call imposeBoundary(grid%octreeRoot, boundarycondition)
  end subroutine hydroStep

  subroutine hydroStep2d(grid, gamma, dt, boundaryCondition)
    type(GRIDTYPE) :: grid
    real(double) :: gamma, dt, subdt
    type(OCTALVECTOR) :: direction
    character(len=*) :: boundaryCondition
    integer :: i


    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)

    call imposeBoundary(grid%octreeRoot, boundarycondition)

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt/2.d0)
    call advectRhoU(grid, direction, dt/2.d0)
    call advectRhoV(grid, direction, dt/2.d0)
    call advectRhoE(grid, direction, dt/2.d0)
    call setupUpm(grid%octreeRoot, grid, direction)
    call computePressureU(grid%octreeRoot, gamma, direction)
    call setupPressure(grid%octreeRoot, grid, direction)
    call pressureForceU(grid%octreeRoot, dt/2.d0)

    direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
    call setupVi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt)
    call advectRhoV(grid, direction, dt)
    call advectRhoU(grid, direction, dt)
    call advectRhoE(grid, direction, dt)
    call setupVi(grid%octreeRoot, grid, direction)
    call setupVpm(grid%octreeRoot, grid, direction)
    call computePressureV(grid%octreeRoot, gamma, direction)
    call setupPressure(grid%octreeRoot, grid, direction)
    call pressureForceV(grid%octreeRoot, dt)

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt/2.d0)
    call advectRhoU(grid, direction, dt/2.d0)
    call advectRhoV(grid, direction, dt/2.d0)
    call advectRhoE(grid, direction, dt/2.d0)
    call setupUpm(grid%octreeRoot, grid, direction)
    call computePressureU(grid%octreeRoot, gamma, direction)
    call setupPressure(grid%octreeRoot, grid, direction)
    call pressureForceU(grid%octreeRoot, dt/2.d0)

    call imposeBoundary(grid%octreeRoot, boundarycondition)

 
  end subroutine hydroStep2d


  recursive subroutine computeCourantTime(thisOctal, tc, gamma)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: tc, dx, cs, gamma, eThermal, eTot, eKinetic, u
  
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

             u = max(abs(thisOctal%rhou(subcell) / thisOctal%rho(subcell)), &
                  abs(thisOctal%rhov(subcell) / thisOctal%rho(subcell)))
             eKinetic = u**2 / 2.d0
             eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
             eThermal = eTot - eKinetic
             cs = sqrt(gamma*(gamma-1.d0)*eThermal)
             dx = thisOctal%subcellSize
             tc = min(tc, dx / (cs + abs(thisOctal%rhou(subcell)/thisOctal%rho(subcell)) ) )
          endif
        
       endif
    enddo
  end subroutine computeCourantTime

  subroutine doHydrodynamics1d(grid)
    type(gridtype) :: grid
    real(double) :: dt, tc, cfl, gamma, mu
    real(double) :: currentTime
    integer :: i, pgbegin
    type(OCTALVECTOR) :: direction

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    cfl = 0.5d0
    mu = 2.d0
    i = pgbegin(0,"/xs",1,1)
!    call pgenv(0., 100., -1., 1., 0, 0)
    call pgenv(0., 1., 0., 1., 0, 0)


    call setupGhostCells(grid%octreeRoot, grid, "mirror")
    call calculateRhoU(grid%octreeRoot, direction)
    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    currentTime = 0.d0
!    do while(currentTime < 0.2d0)
    do while(.true.)
       tc = 1.d30
       call computeCourantTime(grid%octreeRoot, tc, gamma)
       dt = tc * cfl
       call hydroStep(grid, gamma, dt, direction, "mirror")
       currentTime = currentTime + dt
       write(*,*) "current time ",currentTime
       call plotHydroResults(grid)
    enddo
    call pgend
  end subroutine doHydrodynamics1d

  subroutine doHydrodynamics2d(grid)
    type(gridtype) :: grid
    real(double) :: dt, tc, cfl, gamma, mu
    real(double) :: currentTime
    integer :: i, pgbegin, it
    character(len=20) :: plotfile
    real(double) :: tDump, nextDumpTime
    type(OCTALVECTOR) :: direction
    logical :: gridConverged

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    cfl = 0.3d0
    mu = 2.d0
    i = pgbegin(0,"/xs",1,1)
    call pgenv(0., 1., 0., 1., 0, 0)



    call setupGhostCells(grid%octreeRoot, grid, "mirror")
    call refineGridGeneric(grid%octreeRoot, grid, "test")
    call writeInfo("Refine done", TRIVIAL)
    do
       gridConverged = .true.
       call  evenUpGrid(grid%octreeRoot, grid,  gridconverged)
       if (gridConverged) exit
    end do
    call writeInfo("...grid smoothing complete", TRIVIAL)

       call plot_AMR_values(grid, "rho", "x-z", 0., &
            "/xs",.false., .true., fixvalmin=0.d0, fixvalmax=1.d0)

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoV(grid%octreeRoot, direction)

    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    currentTime = 0.d0
    it = 0
    nextDumpTime = 0.d0
    tDump = 0.005d0


!    do while(currentTime < 0.2d0)
    do while(.true.)
       call plot_AMR_values(grid, "rho", "x-z", 0., &
            "/xs",.false., .false., fixvalmin=0.d0, fixvalmax=1.d0)


       tc = 1.d30
       call computeCourantTime(grid%octreeRoot, tc, gamma)
       dt = tc * cfl
       call hydroStep2d(grid, gamma, dt, "mirror")

       call refineGridGeneric(grid%octreeRoot, grid, "test")
       call writeInfo("Refine done", TRIVIAL)
       do
          gridConverged = .true.
          call  evenUpGrid(grid%octreeRoot, grid,  gridconverged)
          if (gridConverged) exit
       end do
       call writeInfo("...grid smoothing complete", TRIVIAL)
       call setupGhostCells(grid%octreeRoot, grid, "mirror")
       call plot_AMR_values(grid, "rho", "x-z", 0., &
            "/xs",.false., .false., fixvalmin=0.d0, fixvalmax=1.d0)

       currentTime = currentTime + dt
       write(*,*) "current time ",currentTime,dt
!       call plotHydroResults(grid)
       if (currentTime .gt. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
          write(plotfile,'(a,i4.4,a)') "image",it,".png/png"
          call plot_AMR_values(grid, "rho", "x-z", 0., &
               plotfile,.false., .false., fixvalmin=0.d0, fixvalmax=1.d0)
       endif
    enddo
  end subroutine doHydrodynamics2d

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
          thisOctal%energy(subcell) = thisOctal%pressure_i(subcell) / ((gamma-1.d0)*thisOctal%rho(subcell))
!          thisOctal%energy(subcell) = 1.d0!thisOctal%energy(subcell) + 0.5d0 * (modulus(thisOctal%velocity(subcell))*cspeed)**2
        
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

  recursive subroutine calculateRhoV(thisOctal, direction)
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
                call calculateRhoV(child, direction)
                exit
             end if
          end do
       else 
          thisOctal%rhoV(subcell) = thisOctal%rho(subcell) * &
               (cSpeed * (thisOctal%velocity(subcell).dot.direction))
       endif
    enddo
  end subroutine calculateRhoV

  recursive subroutine calculateRhoE(thisOctal, direction)
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
                call calculateRhoE(child, direction)
                exit
             end if
          end do
       else 
          thisOctal%rhoE(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       endif
    enddo
  end subroutine calculateRhoE

  subroutine plotHydroResults(grid)

    type(GRIDTYPE) :: grid
    real(double) :: x(100000), rho(100000), v(100000), rhou(100000)
    real,save :: xr(100000), yr(100000)
    real,save :: xd(100000), yd(100000)
    logical, save :: firstTime = .true.
    integer,save :: n
    integer :: i
    
    if (firstTime) then
       open(20,file="sod.dat", form="formatted",status="old")
       do i = 1 , 19
          read(20,*) xd(i), yd(i)
       enddo
       close(20)
    endif


    if (.not.firstTime) then
       call pgsci(0)
       call pgpoint(n, xr, yr, 30)
    endif
    firstTime = .false.

    n = 0
    call getArray(grid%octreeRoot, x, rho, rhou, v, n)
    xr(1:n) = real(x(1:n))
    yr(1:n) = real(rho(1:n))
!    yr(1:n) = real(rhou(1:n))



!    write(*,*) "xr",xr(1:n)
!    write(*,*) "yr",yr(1:n)
    call pgsci(2)
    call pgline(19,xd,yd)

    call pgsci(1)
    call pgpoint(n, xr, yr, 30)
  end subroutine plotHydroResults

  recursive subroutine getArray(thisOctal, x, rho, rhou, v, n)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: x(:), rho(:), v(:), rhou(:)
    integer :: n
    type(OCTALVECTOR) :: rVec
  
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
          rVec = subcellCentre(thisOctal, subcell)
          x(n) = rVec%x
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

    locator = (dble(x1) + 0.1d0*grid%halfSmallestSubcell) * direction
    thisOctal => grid%octreeRoot
    call findSubcellLocal(locator, thisOctal, subcell)
    dx = thisOctal%subcellSize
    thisOctal%ghostCell(subcell) = .true.
    thisOctal%boundaryPartner(subcell) = subcellCentre(thisOctal, subcell) + (3.d0*dx)*direction
    thisOctal%boundaryCondition(subcell) = "mirror"

    locator = (dble(x1)+dx+0.1d0*grid%halfSmallestSubcell) * direction
    thisOctal => grid%octreeRoot
    call findSubcellLocal(locator, thisOctal, subcell)
    thisOctal%ghostCell(subcell) = .true.
    thisOctal%boundaryPartner(subcell) = subcellCentre(thisOctal, subcell) + (1.d0*dx)*direction
    thisOctal%boundaryCondition(subcell) = "mirror"

    locator = (dble(x2) - 0.1d0*grid%halfSmallestSubcell) * direction
    thisOctal => grid%octreeRoot
    call findSubcellLocal(locator, thisOctal, subcell)
    thisOctal%ghostCell(subcell) = .true.
    thisOctal%boundaryPartner(subcell) = subcellCentre(thisOctal, subcell) - (3.d0*dx)*direction
    thisOctal%boundaryCondition(subcell) = "mirror"

    locator = (dble(x2) - dx - 0.1d0*grid%halfSmallestSubcell) * direction 
    thisOctal => grid%octreeRoot
    call findSubcellLocal(locator, thisOctal, subcell)
    thisOctal%ghostCell(subcell) = .true.
    thisOctal%boundaryPartner(subcell) = subcellCentre(thisOctal, subcell) - (1.d0*dx)*direction
    thisOctal%boundaryCondition(subcell) = "mirror"
    
  end subroutine setGhostCells


  recursive subroutine imposeBoundary(thisOctal, boundaryCondition)
    type(octal), pointer   :: thisOctal, bOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, bSubcell
    character(len=*) :: boundaryCondition
    type(OCTALVECTOR) :: locator, dir
    real(double) :: gamma, machNumber, Pr, rhor

    machNumber = 2.d0
    gamma = 7.d0/5.d0
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call imposeBoundary(child,  boundaryCondition)
                exit
             end if
          end do
       else
          if (thisOctal%ghostCell(subcell)) then
             select case(thisOctal%boundaryCondition(subcell))
                case("mirror")
                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocal(locator, bOctal, bSubcell)
!                   if (bOctal%ghostCell(bSubcell)) then
!                      write(*,*)subcellCentre(bOctal, bSubcell)
!                   endif

                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)
                   thisOctal%rho(subcell) = bOctal%rho(bSubcell)
                   thisOctal%rhoE(subcell) = bOctal%rhoE(bSubcell)
! NB confusion regarding 2d being x,z rather than x,y
                   if (abs(dir%x) > 0.9d0) then
                      thisOctal%rhou(subcell) = -bOctal%rhou(bSubcell)
                      thisOctal%rhov(subcell) = bOctal%rhov(bSubcell)
                   endif
                   if (abs(dir%z) > 0.9d0) then
                      thisOctal%rhou(subcell) = bOctal%rhou(bSubcell)
                      thisOctal%rhov(subcell) = -bOctal%rhov(bSubcell)
                   endif
                   if ((abs(dir%x) > 0.2d0).and.abs(dir%z) > 0.2d0) then
                      thisOctal%rhou(subcell) = -bOctal%rhou(bSubcell)
                      thisOctal%rhov(subcell) = -bOctal%rhov(bSubcell)
                   endif
                     
                case("shock")
                   rhor = 1.d0
                   Pr = 0.1d0
                   thisOctal%rho(subcell) = rhor * (gamma+1.d0)*MachNumber**2 / ((gamma-1.d0)*machNumber**2 + 2.d0)
                   thisOctal%pressure_i(subcell) = Pr * (2.d0 * gamma * machNumber**2 - (gamma-1.d0))/(gamma + 1.d0)
                   thisOctal%rhou(subcell) = thisOctal%rho(subcell)*machNumber*sqrt(gamma*0.1d0 / 1.d0) * &
                         (thisOctal%rho(subcell)-rhor)/thisOctal%rho(subcell)

                   thisOctal%energy(subcell) = 0.5d0 * (thisOctal%rhou(subcell) / thisOctal%rho(subcell))**2
                   thisOctal%energy(subcell) = thisOctal%energy(subcell) + &
                        thisOctal%pressure_i(subcell)/((gamma-1.d0)*thisOctal%rho(subcell))
                   thisOctal%rhoe(subcell) =  thisOctal%rho(subcell)*thisOctal%energy(subcell)

                case DEFAULT
                   write(*,*) "Unrecognised boundary condition"
             end select
          endif
      endif
    enddo
  end subroutine imposeBoundary

  recursive subroutine setupGhostCells(thisOctal, grid, boundaryCondition)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    character(len=*) :: boundaryCondition
    type(OCTALVECTOR) :: locator, rVec
    integer :: nProbes, iProbe
    type(OCTALVECTOR) :: probe(6), direction
    real(double) :: dx
    integer :: nProbeOutside
    logical :: corner

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupGhostCells(child,  grid, boundaryCondition)
                exit
             end if
          end do
       else
          if (thisOctal%oned) then
             nProbes = 2
             probe(1) = OCTALVECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
          endif
          if (thisOctal%twod) then
             nProbes = 4
             probe(1) = OCTALVECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
             probe(3) = OCTALVECTOR(0.d0, 0.d0, 1.d0)
             probe(4) = OCTALVECTOR(0.d0, 0.d0, -1.d0)
          endif
          if (thisOctal%threeD) then
             nProbes = 6
             probe(1) = OCTALVECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
             probe(3) = OCTALVECTOR(0.d0, 1.d0, 0.d0)
             probe(4) = OCTALVECTOR(0.d0, -1.d0, 0.d0)
             probe(5) = OCTALVECTOR(0.d0, 0.d0, 1.d0)
             probe(6) = OCTALVECTOR(0.d0, 0.d0, -1.d0)
          endif
          rVec = subcellCentre(thisOctal, subcell)
          nProbeOutside = 0
          do iProbe = 1, nProbes
             locator = rVec + &
                  (thisOctal%subcellsize/2.d0 + 0.1d0*grid%halfSmallestSubcell)*probe(iProbe)
             if (.not.inOctal(grid%octreeRoot, locator).or. &
                  (locator%x < (grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize))) &
                  nProbeOutside = nProbeOutside + 1
          enddo
          corner = .false.
          if (thisOctal%twoD.and.(nProbeOutside > 1)) corner = .true.
          if (thisOctal%threeD.and.(nProbeOutside > 2)) corner = .true.
          if (corner) write(*,*) "corner found"
          if (.not.corner) then
             do iProbe = 1, nProbes
                locator = rVec + &
                     (thisOctal%subcellsize/2.d0 + 0.1d0*grid%halfSmallestSubcell)*probe(iProbe)
                ! this is fudged here because AMR mod assumes that the grid is one-d spherical....
                if (.not.inOctal(grid%octreeRoot, locator).or. &
                     (locator%x < (grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize))) then  ! this is a boundary
                   ! setup the the two ghosts near the boundary
                   
                   thisOctal%ghostCell(subcell) = .true.
                   locator = rVec - &
                        (thisOctal%subcellsize/2.d0 + 0.1d0*grid%halfSmallestSubcell)*probe(iProbe)
                   neighbourOctal => thisOctal
                   call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                   neighbourOctal%ghostCell(neighbourSubcell) = .true.
                   
                   
                   ! now a case to determine the boundary cell relations
                   select case (boundaryCondition)
                   case("mirror")
                      dx = thisOctal%subcellSize
                      thisOctal%ghostCell(subcell) = .true.
                      thisOctal%boundaryPartner(subcell) = subcellCentre(thisOctal, subcell) - (3.d0*dx)*probe(iProbe)
                      thisOctal%boundaryCondition(subcell) = "mirror"
                      
                      neighbourOctal%boundaryPartner(neighboursubcell) = &
                           subcellCentre(neighbourOctal, neighbourSubcell) - (1.d0*dx)*probe(iProbe)
                      neighbourOctal%boundaryCondition(neighboursubcell) = "mirror"
                      
                   case DEFAULT
                      write(*,*) "unknown boundary condition: ",trim(boundaryCondition)
                   end select
                endif
             enddo
          else ! this is a corner
             direction = subcellCentre(thisOctal, subcell) - grid%octreeRoot%centre
             call normalize(direction)

                   
             thisOctal%ghostCell(subcell) = .true.
             locator = rVec - &
                  (sqrt(2.d0)*thisOctal%subcellsize/2.d0 + 0.1d0*grid%halfSmallestSubcell)*direction
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             neighbourOctal%ghostCell(neighbourSubcell) = .true.
                   
                   
             ! now a case to determine the boundary cell relations
             select case (boundaryCondition)
             case("mirror")
                dx = sqrt(2.d0)*thisOctal%subcellSize
                thisOctal%ghostCell(subcell) = .true.
                thisOctal%boundaryPartner(subcell) = subcellCentre(thisOctal, subcell) - (3.d0*dx)*direction
                thisOctal%boundaryCondition(subcell) = "mirror"
                
                neighbourOctal%boundaryPartner(neighboursubcell) = &
                     subcellCentre(neighbourOctal, neighbourSubcell) - (1.d0*dx)*direction
                neighbourOctal%boundaryCondition(neighboursubcell) = "mirror"
                      
             case DEFAULT
                write(*,*) "unknown boundary condition: ",trim(boundaryCondition)
             end select


          endif
      endif
    enddo
  end subroutine setupGhostCells


  subroutine averageValue(thisOctal, subcell, neighbourOctal, q, rhou, rhov, rho)
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell
    integer :: nSubcell(4)
    real(double) :: q, rho, rhov, rhou
    type(OCTALVECTOR) :: direction

    direction = neighbourOctal%centre - subcellCentre(thisOctal, subcell)
    call normalize(direction)
    write(*,*) direction

    if (thisOctal%oneD) then
       if (direction%x > 0.9d0) then
          nSubcell(1) = 1
       else
          nSubcell(2) = 1
       endif
       q = neighbourOctal%q_i(nSubcell(1))
       rho = neighbourOctal%rho(nSubcell(1))
       rhou = neighbourOctal%rhou(nSubcell(1))
       rhov = neighbourOctal%rhov(nSubcell(1))
     else if (thisOctal%twoD) then
       if (direction%x > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 3
       else if (direction%x < -0.9d0) then
          nSubcell(1) = 2
          nSubcell(2) = 4
       else if (direction%y > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 2
       else
          nSubcell(1) = 3
          nSubcell(2) = 4
       endif
       q = 0.5d0*(neighbourOctal%q_i(nSubcell(1)) + neighbourOctal%q_i(nSubcell(2)))
       rho = 0.5d0*(neighbourOctal%rho(nSubcell(1)) + neighbourOctal%rho(nSubcell(2)))
       rhou = 0.5d0*(neighbourOctal%rhou(nSubcell(1)) + neighbourOctal%rhou(nSubcell(2)))
       rhov = 0.5d0*(neighbourOctal%rhov(nSubcell(1)) + neighbourOctal%rhov(nSubcell(2)))
    endif
  end subroutine averageValue


  recursive subroutine refineGridGeneric(thisOctal, grid, criterion)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    type(OCTALVECTOR) :: rVec, corner
    integer :: i, j, subcell, iStream
    logical :: split
    integer :: n
    character(len=*) :: criterion

    subcell = 1
    do while (subcell <= thisOctal%maxChildren)
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          children : do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call refineGridGeneric(child, grid, criterion)
                exit children
             end if
          end do children
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             call  splitCondition(thisOctal, grid, subcell, criterion, split)
             
             if (split.and.(thisOctal%nDepth < 8)) then
                call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                     inherit=.false., interp=.false.)
                subcell = subcell - 1
             endif
          end if
       endif
       subcell = subcell + 1
    enddo

  end subroutine refineGridGeneric

  
  subroutine splitCondition(thisOctal, grid, subcell, criterion, split)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    logical :: split
    character(len=*) :: criterion
    integer :: nProbes
    type(OCTALVECTOR) :: locator, probe(6)
    integer :: i
    real(double) :: grad, maxGradient
    
    Split = .false.
    if (thisOctal%oned) then
       nProbes = 2
       probe(1) = OCTALVECTOR(1.d0, 0.d0, 0.d0)
       probe(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
    endif
    if (thisOctal%twod) then
       nProbes = 4
       probe(1) = OCTALVECTOR(1.d0, 0.d0, 0.d0)
       probe(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
       probe(3) = OCTALVECTOR(0.d0, 0.d0, 1.d0)
       probe(4) = OCTALVECTOR(0.d0, 0.d0, -1.d0)
    endif
    if (thisOctal%threeD) then
       nProbes = 6
       probe(1) = OCTALVECTOR(1.d0, 0.d0, 0.d0)
       probe(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
       probe(3) = OCTALVECTOR(0.d0, 1.d0, 0.d0)
       probe(4) = OCTALVECTOR(0.d0, -1.d0, 0.d0)
       probe(5) = OCTALVECTOR(0.d0, 0.d0, 1.d0)
       probe(6) = OCTALVECTOR(0.d0, 0.d0, -1.d0)
    endif

    maxGradient = 1.d-30
    do i = 1, nProbes
       locator = subcellCentre(thisOctal, subcell) + &
            (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell) * probe(i)
       if (inOctal(grid%octreeRoot, locator)) then
          neighbourOctal => thisOctal
          call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
          grad = abs((thisOctal%rho(subcell) - neighbourOctal%rho(neighbourSubcell)) / thisOctal%rho(subcell))
          maxGradient = max(grad, maxGradient)
       endif
    end do
    if (maxGradient > 0.5d0) split = .true.
  end subroutine splitCondition



  recursive subroutine evenUpGrid(thisOctal, grid,  converged)

    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    !
    integer :: subcell, i, ilambda
    logical :: converged, converged_tmp
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    converged = .true.
    converged_tmp=.true.

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call evenUpGrid(child, grid,  converged)
                if (.not.converged) converged_tmp = converged
                exit
             end if
          end do
          if (.not.converged_tmp) converged=converged_tmp
       else

          r = thisOctal%subcellSize/2.d0 + 0.1d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell)
          if (thisOctal%threed) then
             nDir = 6
             dirVec(1) = OCTALVECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = OCTALVECTOR( 0.d0,+1.d0,  0.d0)
             dirVec(3) = OCTALVECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(4) = OCTALVECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(5) = OCTALVECTOR( 0.d0,-1.d0,  0.d0)
             dirVec(6) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = OCTALVECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = OCTALVECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = OCTALVECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = OCTALVECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
          endif

          if (.not.thisOctal%ghostcell(subcell)) then
             do j = 1, nDir
                octVec = centre + r * dirvec(j)
                if (inOctal(grid%octreeRoot, octVec)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)
                   if ((neighbourOctal%nDepth-thisOctal%nDepth) > 1) then
                      call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                           inherit=.false., interp=.false.)
                      converged = .false.
                      exit
                   endif
                   if ((thisOctal%nDepth-neighbourOctal%nDepth) > 1) then
                      if (.not.neighbourOctal%ghostCell(neighbourSubcell)) then
                         call addNewChild(neighbourOctal,neighboursubcell,grid,adjustGridInfo=.TRUE., &
                              inherit=.false., interp=.false.)
                         converged = .false.
                         exit
                      endif
                   endif
                endif
             enddo
          endif
          if (.not.converged) exit
       endif
    end do

  end subroutine evenUpGrid


 
end module hydrodynamics_mod
