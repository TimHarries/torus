! hydrodynamics module added by TJH on 18th June 2007

module hydrodynamics_mod
#ifdef MPI

  use kind_mod
  use constants_mod
  use amr_mod
  use grid_mod
  use source_mod
  use timing
  use mpi_amr_mod

  implicit none


contains


  recursive subroutine fluxLimiter(thisOctal, limiterType)
    include 'mpif.h'
    integer :: myRank, ierr
    character(len=*) :: limiterType
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: a, b, dq
    character(len=80) :: message
 
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

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
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
  
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

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (thisOctal%u_interface(subcell).ge.0.d0) then
             thisOctal%flux_i(subcell) = thisOctal%u_interface(subcell) * thisOctal%q_i_minus_1(subcell)
          else
             thisOctal%flux_i(subcell) = thisOctal%u_interface(subcell) * thisOctal%q_i(subcell)
          endif

          if (thisOctal%x_i(subcell) == thisOctal%x_i_minus_1(subcell)) then
             write(*,*) "problem with the x_i values"
             stop
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
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
  
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

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%q_i(subcell) = thisOctal%q_i(subcell) - dt * &
                  (thisOctal%flux_i_plus_1(subcell) - thisOctal%flux_i(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i(subcell))
             write(*,*) "fluxes ", thisOctal%flux_i_plus_1(subcell),thisOctal%flux_i(subcell)
             write(*,*) "rho ",thisOctal%rho(subcell) 
!             thisOctal%q_i(subcell) = thisOctal%q_i(subcell) - dt * &
!                  (thisOctal%flux_i_plus_1(subcell) - thisOctal%flux_i(subcell)) / &
!                  (thisOctal%subcellSize*1.d10)

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
          thisOctal%x_i(subcell) = (subcellCentre(thisOctal, subcell) .dot. direction) * 1.d10
       endif
    enddo
  end subroutine setupX

  recursive subroutine setupQX(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator, reverseDirection
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

!          if (associated(thisOctal%mpiBoundaryStorage)) then
!             if (myrank == 1) then
!                write(*,*) "x_i", &
!                  thisOctal%x_i(subcell),thisOctal%mpiBoundaryStorage(subcell, 1:6, 7)
!                write(*,*) "subcell x_i ",thisOctal%x_i(subcell)
!             endif
!          endif

          thisOctal%x_i_minus_1(subcell) = 0.d0
          thisOctal%x_i_plus_1(subcell) = 0.d0
          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)

             thisOctal%x_i_plus_1(subcell) = x
             thisOctal%q_i_plus_1(subcell) = q

             reverseDirection = (-1.d0) * direction
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%x_i_minus_1(subcell) = x
             thisOctal%q_i_minus_1(subcell) = q
             thisOctal%q_i_minus_2(subcell) = qnext


!             write(*,*) "q: ", thisOctal%q_i_minus_2(subcell), thisOctal%q_i_minus_1(subcell), thisOctal%q_i(subcell), &
!                  thisOctal%q_i_plus_1(subcell)
             if (thisOctal%x_i_plus_1(subcell) == thisOctal%x_i_minus_1(subcell)) then
                write(*,*) myrank," error in setting up x_i values"
                write(*,*) thisOctal%x_i_plus_1(subcell),thisOctal%x_i_minus_1(subcell)
             endif

          endif
       endif
    enddo
  end subroutine setupQX

  recursive subroutine setupUi(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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

          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)

!             if (rho < 1.d-5) then
!                write(*,*) "rho",rho, thisOctal%nDepth,neighbourOctal%nDepth
!                if (associated(thisOctal%mpiBoundaryStorage)) then
!                   write(*,*) "depth",thisOctal%mpiBoundaryStorage(subcell,1:6,9)
!                   write(*,*) "rho",thisOctal%mpiBoundaryStorage(subcell,1:6,2)
!                endif
!             endif
             rho_i_minus_1 = rho
             rhou_i_minus_1 = rhou
             thisOctal%u_interface(subcell) = 0.5d0 * &
                  (thisOctal%rhou(subcell)/thisOctal%rho(subcell) + rhou_i_minus_1/rho_i_minus_1)

          endif
       endif
    enddo
  end subroutine setupUI

  recursive subroutine setupVi(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)

                rho_i_minus_1 = rho
                rhou_i_minus_1 = rhov
                thisOctal%u_interface(subcell) = 0.5d0 * &
                     (thisOctal%rhov(subcell)/thisOctal%rho(subcell) + rhou_i_minus_1/rho_i_minus_1)

          endif
       endif
    enddo
  end subroutine setupVI

  recursive subroutine setupWi(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupWi(child, grid, direction)
                exit
             end if
          end do
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)

                rho_i_minus_1 = rho
                rhou_i_minus_1 = rhow
                thisOctal%u_interface(subcell) = 0.5d0 * &
                     (thisOctal%rhow(subcell)/thisOctal%rho(subcell) + rhou_i_minus_1/rho_i_minus_1)

          endif
       endif
    enddo
  end subroutine setupWI



  recursive subroutine setupFlux(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%x_i(subcell) = (subcellCentre(thisOctal, subcell) .dot. direction) * 1.d10
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%flux_i_plus_1(subcell) = flux

             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%flux_i_minus_1(subcell) = flux
          endif
       endif
    enddo
  end subroutine setupFlux


  recursive subroutine setupPressure(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux
    Type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle


          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%pressure_i_plus_1(subcell) = pressure

             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%pressure_i_minus_1(subcell) = pressure

          endif
       endif
    enddo
  end subroutine setupPressure

  recursive subroutine setupUpm(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux
  

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%u_i_plus_1(subcell) = rhou/rho
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%u_i_minus_1(subcell) = rhou/rho

          endif
       endif
    enddo
  end subroutine setupUpm

  recursive subroutine setupVpm(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux
  

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal =>thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%u_i_plus_1(subcell) = rhov/rho
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%u_i_minus_1(subcell) = rhov/rho

          endif
       endif
    enddo
  end subroutine setupVpm

  recursive subroutine setupWpm(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux
  

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupWpm(child, grid, direction)
                exit
             end if
          end do
       else
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%u_i_plus_1(subcell) = rhow/rho
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux)
             thisOctal%u_i_minus_1(subcell) = rhow/rho

          endif
       endif
    enddo
  end subroutine setupWpm


   recursive subroutine computePressureU(thisOctal, gamma, direction, iEquationOfState)
     include 'mpif.h'
     integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: gamma, eThermal, eKinetic, eTot, u, bigGamma,eta, u2
    type(OCTALVECTOR) :: direction, locator
    integer :: iEquationOfState

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    eta = 3.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computePressureU(child, gamma, direction, iEquationOfState)
                exit
             end if
          end do
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          
          select case(iEquationOfState)
             case(0)
                u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
                     thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
                eKinetic = u2 / 2.d0
                eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
                eThermal = eTot - eKinetic
             case(1)
                eThermal = thisOctal%rhoe(subcell) / thisOctal%rho(subcell)
             case DEFAULT
                write(*,*) "Equation of state not recognised: ", iEquationOfState
          end select
                

          thisOctal%pressure_i(subcell) = (gamma - 1.d0) * thisOctal%rho(subcell) * eThermal
          bigGamma = 0.d0
          if (.not.thisOctal%ghostCell(subcell)) then
             if (thisOctal%u_i_plus_1(subcell) .le. thisOctal%u_i_minus_1(subcell)) then
                bigGamma = 0.25d0 * eta**2 * (thisOctal%u_i_plus_1(subcell) - thisOctal%u_i_minus_1(subcell))**2 &
                     * thisOctal%rho(subcell)
             else
                bigGamma = 0.d0
             endif
          endif
          thisOctal%pressure_i(subcell) = thisOctal%pressure_i(subcell) + bigGamma
          if (isnan(thisOctal%pressure_i(subcell))) then
             write(*,*) "pressureU has nan"
             write(*,*) thisOctal%rhou(subcell),thisOctal%rhov(subcell), thisOctal%rhow(subcell)
             write(*,*) thisOctal%rho(subcell)
             stop
          endif

       endif
    enddo
  end subroutine computePressureU

   recursive subroutine computePressureV(thisOctal, gamma, direction,  iEquationOfState)
     include 'mpif.h'
     integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: gamma, eThermal, eKinetic, eTot, u, bigGamma,eta,u2
    type(OCTALVECTOR) :: direction, locator
    integer ::  iEquationOfState

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    eta = 3.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computePressureV(child, gamma, direction, iEquationofState)
                exit
             end if
          end do
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          select case(iEquationOfState)
             case(0)
                u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
                     thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
                eKinetic = u2 / 2.d0
                eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
                eThermal = eTot - eKinetic
             case(1)
                eThermal = thisOctal%rhoe(subcell) / thisOctal%rho(subcell)
             case DEFAULT
                write(*,*) "Equation of state not recognised: ", iEquationOfState
          end select
          
!          u = thisOctal%rhov(subcell) / thisOctal%rho(subcell)
!          eKinetic = u**2 / 2.d0

          thisOctal%pressure_i(subcell) = (gamma - 1.d0) * thisOctal%rho(subcell) * eThermal
          bigGamma = 0.d0
          if (.not.thisOctal%ghostCell(subcell)) then
             if (thisOctal%u_i_plus_1(subcell) .le. thisOctal%u_i_minus_1(subcell)) then
                bigGamma = 0.25d0 * eta**2 * (thisOctal%u_i_plus_1(subcell) - thisOctal%u_i_minus_1(subcell))**2 &
                     * thisOctal%rho(subcell)
             else
                bigGamma = 0.d0
             endif
          endif
          thisOctal%pressure_i(subcell) = thisOctal%pressure_i(subcell) + bigGamma

          if (thisOctal%pressure_i(subcell) < 0.d0) then
             write(*,*) "pressure ", thisOctal%pressure_i(subcell)
             write(*,*) "rho ", thisOctal%rho(subcell)
             write(*,*) "biggamma ", biggamma
             write(*,*) "ethermal ", ethermal
             write(*,*) "EOS ", iEquationofState
          endif
          if (isnan(thisOctal%pressure_i(subcell))) then
             write(*,*) "pressureV has nan"
             write(*,*) thisOctal%rhou(subcell),thisOctal%rhov(subcell), thisOctal%rhow(subcell)
             write(*,*) thisOctal%rho(subcell)
             stop
          endif
       endif
    enddo
  end subroutine computePressureV

   recursive subroutine computePressureW(thisOctal, gamma, direction,  iEquationOfState)
     include 'mpif.h'
     integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: gamma, eThermal, eKinetic, eTot, u, bigGamma,eta,u2
    type(OCTALVECTOR) :: direction, locator
    integer :: iEquationOfState


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    eta = 3.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computePressureW(child, gamma, direction, iEquationOfState)
                exit
             end if
          end do
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          select case(iEquationOfState)
             case(0)
                u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
                     thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
                eKinetic = u2 / 2.d0
                eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
                eThermal = eTot - eKinetic
             case(1)
                eThermal = thisOctal%rhoe(subcell) / thisOctal%rho(subcell)
             case DEFAULT
                write(*,*) "Equation of state not recognised: ", iEquationOfState
          end select
          
!          u = thisOctal%rhow(subcell) / thisOctal%rho(subcell)
!          eKinetic = u**2 / 2.d0

          thisOctal%pressure_i(subcell) = (gamma - 1.d0) * thisOctal%rho(subcell) * eThermal
          bigGamma = 0.d0
          if (.not.thisOctal%ghostCell(subcell)) then
             if (thisOctal%u_i_plus_1(subcell) .le. thisOctal%u_i_minus_1(subcell)) then
                bigGamma = 0.25d0 * eta**2 * (thisOctal%u_i_plus_1(subcell) - thisOctal%u_i_minus_1(subcell))**2 &
                     * thisOctal%rho(subcell)
             else
                bigGamma = 0.d0
             endif
          endif
          thisOctal%pressure_i(subcell) = thisOctal%pressure_i(subcell) + bigGamma
          if (isnan(thisOctal%pressure_i(subcell))) then
             write(*,*) "pressureW has nan"
             write(*,*) thisOctal%rhou(subcell),thisOctal%rhov(subcell), thisOctal%rhow(subcell)
             write(*,*) thisOctal%rho(subcell)
             stop
          endif
       endif
    enddo
  end subroutine computePressureW


  recursive subroutine pressureForceU(thisOctal, dt, iEquationOfState)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, tmp
    integer :: iEquationOfState


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call pressureForceU(child, dt, iEquationOfState)
                exit
             end if
          end do
       else
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then


             if (thisOctal%x_i_plus_1(subcell) == thisOctal%x_i_minus_1(subcell)) then
                write(*,*) myrank," error in setting up x_i values"
                write(*,*) thisOctal%x_i_plus_1(subcell),thisOctal%x_i_minus_1(subcell), thisOctal%x_i(subcell)
                write(*,*) thisOctal%nDepth
                write(*,*) "centre ",subcellCentre(thisOctal,subcell)
                write(*,*) "thread ",thisOctal%mpiThread(subcell)
             endif

             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))


             if (thisOctal%rhou(subcell)/thisOctal%rho(subcell) > 1.d10) then
                write(*,*) "u ",thisOctal%rhou(subcell)/thisOctal%rho(subcell)/1.d5
                write(*,*) "rho ", thisOctal%rho(subcell)
                write(*,*) "pressure i + 1 " ,thisOctal%pressure_i_plus_1(subcell)
                write(*,*) "pressure i - 1 " ,thisOctal%pressure_i_minus_1(subcell)
                write(*,*) "x i + 1 ", thisOctal%x_i_plus_1(subcell)
                write(*,*) "x i - 1 ", thisOctal%x_i_minus_1(subcell)
             endif

                
!             Write(*,*) thisOctal%rhou(subcell), thisOctal%pressure_i_plus_1(subcell), &
!                  thisOctal%pressure_i_minus_1(subcell), thisOctal%x_i_plus_1(subcell), &
!                  thisOctal%x_i_minus_1(subcell)

             if (iEquationOfState == 0) then
                thisOctal%rhoE(subcell) = thisOctal%rhoE(subcell) - dt * &
                     (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
                     thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / &
                     (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))
             endif


             if (isnan(thisOctal%rhou(subcell))) then
                write(*,*) "bug",thisOctal%rhou(subcell), &
                     thisOctal%pressure_i_plus_1(subcell),thisOctal%pressure_i_minus_1(subcell)
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureForceU

  recursive subroutine pressureForceV(thisOctal, dt, iEquationOfState)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt
    integer :: iEquationOfState

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call pressureForceV(child, dt, iEquationOfState)
                exit
             end if
          end do
       else
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then


             if (thisOctal%x_i_plus_1(subcell) == thisOctal%x_i_minus_1(subcell)) then
                write(*,*) myrank," error in setting up x_i values"
                write(*,*) thisOctal%x_i_plus_1(subcell),thisOctal%x_i_minus_1(subcell), thisOctal%x_i(subcell)
                write(*,*) thisOctal%nDepth
                write(*,*) "centre ",subcellCentre(thisOctal,subcell)
                write(*,*) "thread ",thisOctal%mpiThread(subcell)
             endif

             thisOctal%rhov(subcell) = thisOctal%rhov(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))

!             write(*,*) thisOctal%rhou(subcell), thisOctal%pressure_i_plus_1(subcell), &
!                  thisOctal%pressure_i_minus_1(subcell), thisOctal%x_i_plus_1(subcell), &
!                  thisOctal%x_i_minus_1(subcell)

             if (iEquationOfState == 0) then
                thisOctal%rhoE(subcell) = thisOctal%rhoE(subcell) - dt * &
                     (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
                     thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / &
                     (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))
             endif

             if (isnan(thisOctal%rhov(subcell))) then
                write(*,*) "bug",thisOctal%pressure_i_plus_1(subcell),thisOctal%pressure_i_minus_1(subcell)
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureForceV

  recursive subroutine pressureForceW(thisOctal, dt, iEquationOfState)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt
    integer :: iEquationOfState

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call pressureForceW(child, dt, iEquationOfState)
                exit
             end if
          end do
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle


          if (.not.thisOctal%ghostCell(subcell)) then


             if (thisOctal%x_i_plus_1(subcell) == thisOctal%x_i_minus_1(subcell)) then
                write(*,*) myrank," error in setting up x_i values"
                write(*,*) thisOctal%x_i_plus_1(subcell),thisOctal%x_i_minus_1(subcell), thisOctal%x_i(subcell)
                write(*,*) thisOctal%nDepth
                write(*,*) "centre ",subcellCentre(thisOctal,subcell)
                write(*,*) "thread ",thisOctal%mpiThread(subcell)
             endif

             thisOctal%rhoW(subcell) = thisOctal%rhoW(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))

!             write(*,*) thisOctal%rhou(subcell), thisOctal%pressure_i_plus_1(subcell), &
!                  thisOctal%pressure_i_minus_1(subcell), thisOctal%x_i_plus_1(subcell), &
!                  thisOctal%x_i_minus_1(subcell)

             if (iEquationOfState == 0) then
                thisOctal%rhoE(subcell) = thisOctal%rhoE(subcell) - dt * &
                     (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
                     thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / &
                     (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))
             endif


             if (isnan(thisOctal%rhoW(subcell))) then
                write(*,*) "bug",thisOctal%pressure_i_plus_1(subcell),thisOctal%pressure_i_minus_1(subcell)
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureForceW




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

  recursive subroutine copyRhoWToQ(thisOctal)
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
                call copyRhoWToQ(child)
                exit
             end if
          end do
       else
  
          thisOctal%q_i(subcell) = thisOctal%rhoW(subcell)
        
       endif
    enddo
  end subroutine copyRhoWToQ

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

  recursive subroutine copyQtoRhoW(thisOctal)
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
                call copyQtoRhoW(child)
                exit
             end if
          end do
       else
  
          thisOctal%rhoW(subcell) = thisOctal%q_i(subcell)
        
       endif
    enddo
  end subroutine copyQtoRhoW

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
          if (thisOCtal%rho(subcell) < 0.d0) then
             write(*,*) "rho warning ", thisOctal%rho(subcell), subcellCentre(thisOctal, subcell)
          endif
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


  subroutine advectRho(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup
    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupQX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRho(grid%octreeRoot)

  end subroutine advectRho

  subroutine advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoEToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupQX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRhoE(grid%octreeRoot)

  end subroutine advectRhoE

  subroutine advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoUToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupQX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRhoU(grid%octreeRoot)

  end subroutine advectRhoU

  subroutine advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoVToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupQX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRhoV(grid%octreeRoot)

  end subroutine advectRhoV

  subroutine advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(OCTALVECTOR) :: direction

    call copyRhoWToQ(grid%octreeRoot)
    call setupX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupQX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)
    call copyQtoRhoW(grid%octreeRoot)

  end subroutine advectRhoW


  subroutine hydroStep(grid, gamma, dt, direction)
    type(GRIDTYPE) :: grid
    real(double) :: gamma, dt
    type(OCTALVECTOR) :: direction

!    call imposeBoundary(grid%octreeRoot, boundarycondition)
!    call setupUi(grid%octreeRoot, grid, direction)
!    call advectRho(grid, direction, dt)
!    call advectRhoU(grid, direction, dt)
!    call advectRhoE(grid, direction, dt)
!    call setupUpm(grid%octreeRoot, grid, direction)
!    call computePressureU(grid%octreeRoot, gamma, direction)
!    call setupPressure(grid%octreeRoot, grid, direction)
!    call pressureForceU(grid%octreeRoot, dt)
!    call imposeBoundary(grid%octreeRoot, boundarycondition)

  end subroutine hydroStep


  subroutine hydroStep3d(iEquationOfState, grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup
    integer :: iEquationOfState

    real(double) :: gamma, dt, subdt
    type(OCTALVECTOR) :: direction
    integer :: i


    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)

    call imposeBoundary(grid%octreeRoot)
    call transferTempStorage(grid%octreeRoot)


    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoV(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUpm(grid%octreeRoot, grid, direction)
    call computePressureU(grid%octreeRoot, gamma, direction, iEquationOfState)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceU(grid%octreeRoot, dt/2.d0, iEquationOfState)


    direction = OCTALVECTOR(0.d0, 1.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupVi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupVpm(grid%octreeRoot, grid, direction)
    call computePressureV(grid%octreeRoot, gamma, direction, iEquationOfState)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceV(grid%octreeRoot, dt, iEquationOfState)

    direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupWi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupWpm(grid%octreeRoot, grid, direction)
    call computePressureW(grid%octreeRoot, gamma, direction, iEquationOfState)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceW(grid%octreeRoot, dt, iEquationOfState)

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoV(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUpm(grid%octreeRoot, grid, direction)
    call computePressureU(grid%octreeRoot, gamma, direction, iEquationOfState)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceU(grid%octreeRoot, dt/2.d0, iEquationOfState)

    call imposeBoundary(grid%octreeRoot)
    call transferTempStorage(grid%octreeRoot)

 
  end subroutine hydroStep3d

  subroutine hydroStep2d(grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    real(double) :: gamma, dt, subdt
    type(OCTALVECTOR) :: direction
    integer :: i


    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)

!    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)

    call plotGridMPI(grid, "rhoafter.png/png", "x-z", "rho", plotgrid=.true.)

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUpm(grid%octreeRoot, grid, direction)
    call computePressureU(grid%octreeRoot, gamma, direction, 0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceU(grid%octreeRoot, dt/2.d0, 0)

    call plotGridMPI(grid, "rhoafterx.png/png", "x-z", "rho", plotgrid=.true.)


    direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupWi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupWpm(grid%octreeRoot, grid, direction)
    call computePressureW(grid%octreeRoot, gamma, direction, 0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceW(grid%octreeRoot, dt,0 )

    call plotGridMPI(grid, "rhoaftery.png/png", "x-z", "rho", plotgrid=.true.)

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUpm(grid%octreeRoot, grid, direction)
    call computePressureU(grid%octreeRoot, gamma, direction, 0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceU(grid%octreeRoot, dt/2.d0, 0)

    call plotGridMPI(grid, "rhoafterx2.png/png", "x-z", "rho", plotgrid=.true.)

    call periodBoundary(grid)
!    call imposeBoundary(grid%octreeRoot)
    call transferTempStorage(grid%octreeRoot)

    call plotGridMPI(grid, "rhofinal.png/png", "x-z", "rho", plotgrid=.true.)
 
  end subroutine hydroStep2d


  recursive subroutine computeCourantTime(thisOctal, tc, gamma, iEquationOfState)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: tc, dx, cs, gamma, eThermal, eTot, eKinetic, u, speed
    integer :: iEquationOfState
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computeCourantTime(child, tc, gamma, iEquationOfState)
                exit
             end if
          end do
       else
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then

             cs = soundSpeed(thisOctal, subcell, gamma, iEquationOfState)
!             if (myrank==1) write(*,*) "cs ", cs/1.d5, "km/s"
             dx = thisOctal%subcellSize * 1.d10
             speed = sqrt((thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 &
                  + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2)
                tc = min(tc, dx / (cs + speed) )

                if ((dx/(cs+speed)) < 1.d6) then
                   write(*,*) "dx ",dx, "speed ", speed/1.d5, "cs ",cs/1.d5
                   write(*,*) "x y z ",thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                        thisOctal%rhov(subcell)/thisOctal%rho(subcell), thisOctal%rhow(subcell)/thisOctal%rho(subcell) 
                endif

          endif
 
       endif
    enddo
  end subroutine computeCourantTime

  function soundSpeed(thisOctal, subcell, gamma, iEquationOfState) result (cs)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: cs, gamma
    real(double) :: u2, eKinetic, eTot, eThermal
    integer :: iEquationOfSate
    logical, save :: firstTime = .true.
    integer :: iEquationOfState

    select case(iEquationOfState)
       case(0) 
          if (thisOctal%threed) then
             u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
          else
             u2 = (thisOctal%rhou(subcell)**2 +  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
          endif
          eKinetic = u2 / 2.d0
          eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
          eThermal = eTot - eKinetic
       case(1)
          eThermal = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
    end select

    if (eThermal < 0.d0) then
       if (firstTime) then
          write(*,*) "eThermal problem: ",ethermal
          write(*,*) eTot, ekinetic, u2
          write(*,*) "rhoe, rho", thisOctal%rhoe(subcell), thisOctal%rho(subcell)
          write(*,*) thisOctal%rhou(subcell)/thisOctal%rho(subcell), thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
               thisOctal%rhow(subcell)/thisOctal%rho(subcell)
          write(*,*) "cen ",subcellCentre(thisOctal, subcell), thisOctal%ghostCell(subcell)
          write(*,*) "temp ",thisOctal%temperature(subcell)
          firstTime = .false.
       endif
       eThermal = tiny(eThermal)
    endif
    cs = sqrt(gamma*(gamma-1.d0)*eThermal)

  end function soundSpeed

  subroutine doHydrodynamics1d(grid)
    type(gridtype) :: grid
    real(double) :: dt,  tc, cfl, gamma, mu
    real(double) :: currentTime
    integer :: i, pgbegin
    type(OCTALVECTOR) :: direction
    logical :: gridConverged
    integer :: iEquationOfState = 0

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    cfl = 0.3d0
    mu = 2.d0
    i = pgbegin(0,"/xs",1,1)
!    call pgenv(0., 100., -1., 1., 0, 0)
    call pgenv(0., 1., 0., 1., 0, 0)

    call setupEdges(grid%octreeRoot, grid)
    call unsetGhosts(grid%octreeRoot)
    call setupGhostCells(grid%octreeRoot, grid)


    call zeroRefinedLastTime(grid%octreeRoot)
    call refineGridGeneric(grid%octreeRoot, grid, "test", .false.)
    call writeInfo("Refine done", TRIVIAL)
    do
       gridConverged = .true.
       call evenUpGrid(grid%octreeRoot, grid,  gridconverged, inherit=.false.)
       call setupEdges(grid%octreeRoot, grid)
       call unsetGhosts(grid%octreeRoot)
       call setupGhostCells(grid%octreeRoot, grid)
       if (gridConverged) exit
    end do
    call writeInfo("...grid smoothing complete", TRIVIAL)

    call calculateRhoU(grid%octreeRoot, direction)
    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    currentTime = 0.d0
    do while(currentTime < 0.2d0)
!    do while(.true.)
       tc = 1.d30
       call computeCourantTime(grid%octreeRoot, tc, gamma, iEquationOfState)
       dt = tc * cfl
       call hydroStep(grid, gamma, dt, direction)

!       write(*,*) "unrefining"
!       do
!          gridConverged = .true.
!          call unrefineCells(grid%octreeRoot, grid, gridconverged, gamma)
!          if (gridConverged) exit
!       end do
!       write(*,*) "done"

       call zeroRefinedLastTime(grid%octreeRoot)
       call refineGridGeneric(grid%octreeRoot, grid, "test", .true.)
       call writeInfo("Refine done", TRIVIAL)
       do
          gridConverged = .true.
          call  evenUpGrid(grid%octreeRoot, grid,  gridconverged, inherit=.true.)
          call setupEdges(grid%octreeRoot, grid)
          call unsetGhosts(grid%octreeRoot)
          call setupGhostCells(grid%octreeRoot, grid)
          if (gridConverged) exit
       end do
       call writeInfo("...grid smoothing complete", TRIVIAL)
       currentTime = currentTime + dt
       write(*,*) "current time ",currentTime
       call plotHydroResults(grid)
    enddo
    call pgend
  end subroutine doHydrodynamics1d


  subroutine doHydrodynamics3d(grid)
    include 'mpif.h'
    type(gridtype) :: grid
    real(double) :: dt, tc(64), temptc(64),cfl, gamma, mu
    real(double) :: currentTime
    integer :: i, pgbegin, it, iUnrefine
    integer :: myRank, ierr
    character(len=20) :: plotfile
    real(double) :: tDump, nextDumpTime, ang
    type(OCTALVECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: nSource = 0
    type(SOURCETYPE) :: source(1)
    integer :: nDependent, dependentThread(100)
    integer :: thread1(200), thread2(200), nBound(200), nPairs
    logical :: globalConverged(64), tConverged(64)
    integer :: nHydroThreads
    integer :: nGroup, group(200)
    integer :: iEquationOfState = 0

    nHydroThreads = nThreadsGlobal - 1

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    cfl = 0.3d0

    cfl = 0.01d0

    mu = 2.d0

    viewVec = OCTALVECTOR(-1.d0,0.d0,0.d0)
!    viewVec = rotateZ(viewVec, 20.d0*degtorad)
    viewVec = rotateY(viewVec, 25.d0*degtorad)
    

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    if (myRank == 1) write(*,*) "CFL set to ", cfl


    call writeInfo("Plotting grid", TRIVIAL)    
    call plotGridMPI(grid, "mpi.ps/vcps", "x-z", "mpi", plotgrid=.true.)



    call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)

    do i = 1, nPairs
       if (myrankglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
    enddo


    call writeInfo("Calling exchange across boundary", TRIVIAL)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call writeInfo("Done", TRIVIAL)

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)

    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)



    call writeInfo("Refining individual subgrids", TRIVIAL)
    if (.not.grid%splitOverMpi) then
       do
          gridConverged = .true.
          call setupEdges(grid%octreeRoot, grid)
          call refineEdges(grid%octreeRoot, grid,  gridconverged, inherit=.false.)
          call unsetGhosts(grid%octreeRoot)
          call setupGhostCells(grid%octreeRoot, grid)
          if (gridConverged) exit
       end do
    else
       call evenUpGridMPI(grid, .false.)
    endif

    call writeInfo("Refining grid", TRIVIAL)
    do
       gridConverged = .true.
       call refineGridGeneric2(grid%octreeRoot, grid,  gamma, gridconverged, iEquationOfState, inherit=.false.)
       if (gridConverged) exit
    end do
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)

    call writeInfo("Refining grid part 2", TRIVIAL)    
    do
       globalConverged(myRank) = .true.
       call writeInfo("Refining grid", TRIVIAL)    
       call refineGridGeneric2(grid%octreeRoot, grid,  gamma, globalConverged(myRank), iEquationOfState, inherit=.false.)
       call writeInfo("Exchanging boundaries", TRIVIAL)    
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreads, MPI_LOGICAL, MPI_LOR,amrCOMMUNICATOR, ierr)
       if (ALL(tConverged(1:nHydroThreads))) exit
    end do


    call writeInfo("Evening up grid", TRIVIAL)    
    call evenUpGridMPI(grid,.false.)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)




    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)
    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)

    currentTime = 0.d0
    it = 0
    nextDumpTime = 0.d0
    tDump = 0.005d0

    iUnrefine = 0
!    call writeInfo("Plotting grid", TRIVIAL)    
!    call plotGridMPI(grid, "mpi.ps/vcps", "x-z", "rhoe", 0., 1.)

    iUnrefine = 0

    call writeInfo("Plotting col density", TRIVIAL)    
    call columnDensityPlotAMR(grid, viewVec, "test.png/png", iminfix = 0., imaxfix = 0.2)

!      call plotGridMPI(grid, "/xs", "x-z", "rho", 0., 1.)

!    do i = 1, 20
!       ang = twoPi * dble(i-1)/20.d0
!       viewVec = OCTALVECTOR(cos(ang), sin(ang), 0.d0)
!       write(plotfile,'(a,i4.4,a)') "image",i,".gif/gif"
!       call columnDensityPlotAMR(grid, viewVec, plotfile, iminfix = 0., imaxfix = 1.)
!    enddo
!    stop

!    do while(currentTime < 0.2d0)
    do while(.true.)
       tc = 0.d0
       tc(myrank) = 1.d30
       call computeCourantTime(grid%octreeRoot, tc(myRank), gamma, iEquationOfState)
       call MPI_ALLREDUCE(tc, tempTc, 8, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
!       write(*,*) "tc", tc(1:8)
!       write(*,*) "temp tc",temptc(1:8)
       tc = tempTc
       dt = MINVAL(tc(1:8)) * cfl

       if (myrank == 1) write(*,*) "courantTime", dt
       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call writeInfo("calling hydro step",TRIVIAL)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       call hydroStep3d(0, grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup)
       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       call writeInfo("Refining grid", TRIVIAL)
       do
          gridConverged = .true.
          call refineGridGeneric2(grid%octreeRoot, grid,  gamma, gridconverged, iEquationOfState, inherit=.true.)
          if (gridConverged) exit
       end do
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       
       call writeInfo("Refining grid part 2", TRIVIAL)    
       do
          globalConverged(myRank) = .true.
          call writeInfo("Refining grid", TRIVIAL)    
          call refineGridGeneric2(grid%octreeRoot, grid,  gamma, globalConverged(myRank), iEquationOfState, inherit=.true.)
          call writeInfo("Exchanging boundaries", TRIVIAL)    
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          call MPI_ALLREDUCE(globalConverged, tConverged, 8, MPI_LOGICAL, MPI_LOR,amrCOMMUNICATOR, ierr)
          if (ALL(tConverged(1:8))) exit
       end do
       
       iUnrefine = iUnrefine + 1
       if (iUnrefine == 5) then
          if (myrank == 1)call tune(6, "Unrefine grid")
          call unrefineCells(grid%octreeRoot, grid, gamma, iEquationOfState)
          if (myrank == 1)call tune(6, "Unrefine grid")
          iUnrefine = 0
       endif

       call evenUpGridMPI(grid, .true.)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


!     call plot_AMR_values(grid, "rho", "x-y", 0., &
 !           "/xs",.false., .true., fixvalmin=0.d0, fixvalmax=1.d0, quiet=.true.)


!       call plotGridMPI(grid, "/xs", "x-z", "rho", 0., 1.)

       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt
       if (currentTime .gt. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
!          write(plotfile,'(a,i4.4,a)') "image",it,".png/png"
!          call columnDensityPlotAMR(grid, viewVec, plotfile, resetRangeFlag=.false.)
          write(plotfile,'(a,i4.4,a)') "rho",it,".png/png"
          call plotGridMPI(grid, plotfile, "x-z", "rho", 0., 1.)
       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)


    enddo
  end subroutine doHydrodynamics3d

  subroutine doHydrodynamics2d(grid)
    include 'mpif.h'
    type(gridtype) :: grid
    real(double) :: dt, tc(64), temptc(64),cfl, gamma, mu
    real(double) :: currentTime
    integer :: i, pgbegin, it, iUnrefine
    integer :: myRank, ierr
    character(len=20) :: plotfile
    real(double) :: tDump, nextDumpTime, ang
    type(OCTALVECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: nSource = 0
    type(SOURCETYPE) :: source(1)
    integer :: nDependent, dependentThread(100)
    integer :: thread1(100), thread2(100), nBound(100), nPairs
    integer :: nGroup, group(100)

    logical :: globalConverged(64), tConverged(64)
    integer :: nHydroThreads 
    integer :: iEquationOfState = 0
    

    nHydroThreads = nThreadsGlobal - 1


    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    cfl = 0.3d0

    mu = 2.d0

    viewVec = OCTALVECTOR(-1.d0,0.d0,0.d0)
!    viewVec = rotateZ(viewVec, 20.d0*degtorad)
    viewVec = rotateY(viewVec, 25.d0*degtorad)
    

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    if (myRank == 1) write(*,*) "CFL set to ", cfl


!    call writeInfo("Plotting grid", TRIVIAL)    
!    write(plotfile,'(a,i2.2,a)') "test",myRank,".png/png"
!    call plot_AMR_values(grid, "rho", "x-z", 0., &
!           plotfile,.false., .true.)

    call plotGridMPI(grid, "mpi.png/png", "x-z", "mpi", plotgrid=.true.)
    call plotGridMPI(grid, "chi.png/png", "x-z", "chi", plotgrid=.true.)



    call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)

    do i = 1, nPairs
       if (myrankglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
    enddo


    call writeInfo("Calling exchange across boundary", TRIVIAL)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call writeInfo("Done", TRIVIAL)

    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)

    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)



    call writeInfo("Refining individual subgrids", TRIVIAL)
    if (.not.grid%splitOverMpi) then
       do
          gridConverged = .true.
          call setupEdges(grid%octreeRoot, grid)
          call refineEdges(grid%octreeRoot, grid,  gridconverged, inherit=.false.)
          call unsetGhosts(grid%octreeRoot)
          call setupGhostCells(grid%octreeRoot, grid, flag=.true.)
          if (gridConverged) exit
       end do
    else
       call evenUpGridMPI(grid,.false.)
    endif

    

    call writeInfo("Refining grid", TRIVIAL)
    do
       gridConverged = .true.
       call refineGridGeneric2(grid%octreeRoot, grid,  gamma, gridconverged, iEquationOfState, inherit=.false.)
       if (gridConverged) exit
    end do
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)

    call writeInfo("Refining grid part 2", TRIVIAL)    
    do
       globalConverged(myRank) = .true.
       call writeInfo("Refining grid", TRIVIAL)    
       call refineGridGeneric2(grid%octreeRoot, grid,  gamma, globalConverged(myRank), iEquationOfState, inherit=.false.)
       call writeInfo("Exchanging boundaries", TRIVIAL)    
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       call MPI_ALLREDUCE(globalConverged, tConverged, 4, MPI_LOGICAL, MPI_LOR,amrCOMMUNICATOR, ierr)
       if (ALL(tConverged(1:4))) exit
    end do


    call writeInfo("Evening up grid", TRIVIAL)    
    call evenUpGridMPI(grid, .false.)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


    call plotGridMPI(grid, "rhostart.png/png", "x-z", "rho", plotgrid=.true.)


    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)
    direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)

    currentTime = 0.d0
    it = 0
    nextDumpTime = 0.d0
    tDump = 0.005d0

    iUnrefine = 0
!    call writeInfo("Plotting grid", TRIVIAL)    
!    call plotGridMPI(grid, "mpi.ps/vcps", "x-z", "rhoe", 0., 1.)

    iUnrefine = 0

!    call writeInfo("Plotting col density", TRIVIAL)    
!    call columnDensityPlotAMR(grid, viewVec, "test.png/png", iminfix = 0., imaxfix = 0.2)

!      call plotGridMPI(grid, "/xs", "x-z", "rho", 0., 1.)

!    do i = 1, 20
!       ang = twoPi * dble(i-1)/20.d0
!       viewVec = OCTALVECTOR(cos(ang), sin(ang), 0.d0)
!       write(plotfile,'(a,i4.4,a)') "image",i,".gif/gif"
!       call columnDensityPlotAMR(grid, viewVec, plotfile, iminfix = 0., imaxfix = 1.)
!    enddo
!    stop

!    do while(currentTime < 0.2d0)

    do while(.true.)
       tc = 0.d0
       tc(myrank) = 1.d30
       call computeCourantTime(grid%octreeRoot, tc(myRank), gamma, iEquationOfState)
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
!       write(*,*) "tc", tc(1:8)
!       write(*,*) "temp tc",temptc(1:8)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads)) * cfl

       if (myrank == 1) write(*,*) "courantTime", dt
       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call writeInfo("calling hydro step",TRIVIAL)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


       call hydroStep2d(grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup)



       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       call writeInfo("Refining grid", TRIVIAL)
       do
          gridConverged = .true.
          call refineGridGeneric2(grid%octreeRoot, grid,  gamma, gridconverged, iEquationOfState, inherit=.true.)
          if (gridConverged) exit
       end do
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       
       call writeInfo("Refining grid part 2", TRIVIAL)    
       do
          globalConverged(myRank) = .true.
          call writeInfo("Refining grid", TRIVIAL)    
          call refineGridGeneric2(grid%octreeRoot, grid,  gamma, globalConverged(myRank), iEquationOfState, inherit=.true.)
          call writeInfo("Exchanging boundaries", TRIVIAL)    
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreads, MPI_LOGICAL, MPI_LOR,amrCOMMUNICATOR, ierr)
          if (ALL(tConverged(1:nHydroThreads))) exit
       end do
       
       iUnrefine = iUnrefine + 1
       if (iUnrefine == 5) then
          call tune(6, "Unrefine grid")
          call unrefineCells(grid%octreeRoot, grid, gamma, iEquationOfState)
          call tune(6, "Unrefine grid")
          iUnrefine = 0
       endif

       call evenUpGridMPI(grid, .true.)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


!     call plot_AMR_values(grid, "rho", "x-y", 0., &
 !           "/xs",.false., .true., fixvalmin=0.d0, fixvalmax=1.d0, quiet=.true.)


!       call plotGridMPI(grid, "/xs", "x-z", "rho", 0., 1.)

       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt
       if (currentTime .gt. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
!          write(plotfile,'(a,i4.4,a)') "image",it,".png/png"
!          call columnDensityPlotAMR(grid, viewVec, plotfile, resetRangeFlag=.false.)
          write(plotfile,'(a,i4.4,a)') "rho",it,".png/png"
          call plotGridMPI(grid, plotfile, "x-z", "rho", 0., 2.,plotgrid=.false.)
!          call plotGridMPI(grid, "/xs", "x-z", "rhoe", plotgrid=.true.)
       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)


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


  recursive subroutine zeroRefinedLastTime(thisOctal)
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
                call zeroRefinedLastTime(child)
                exit
             end if
          end do
       else 
          thisOctal%refinedLastTime = .false.
          thisOctal%parent%refinedLastTime = .false.
       endif
    enddo
  end subroutine zeroRefinedLastTime


  recursive subroutine transferTempStorage(thisOctal)
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
                call transferTempStorage(child)
                exit
             end if
          end do
       else 
          if (thisOctal%ghostCell(subcell)) then
             if (.not.associated(thisOctal%tempStorage)) then
                call writeFatal("Tempstorage not allocated when it should have been")
                stop
             endif
             thisOctal%rho(subcell) = thisOctal%tempStorage(subcell,1)
             thisOctal%rhoe(subcell) = thisOctal%tempStorage(subcell,2)
             thisOctal%rhou(subcell) = thisOctal%tempStorage(subcell,3)
             thisOctal%rhov(subcell) = thisOctal%tempStorage(subcell,4)
             thisOctal%rhow(subcell) = thisOctal%tempStorage(subcell,5)
          endif
       endif
    enddo
    if (associated(thisOCtal%tempStorage)) Deallocate(thisOctal%tempStorage)
  end subroutine transferTempStorage

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

  recursive subroutine calculateRhoW(thisOctal, direction)
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
                call calculateRhoW(child, direction)
                exit
             end if
          end do
       else 
          thisOctal%rhoW(subcell) = thisOctal%rho(subcell) * &
               (cSpeed * (thisOctal%velocity(subcell).dot.direction))
       endif
    enddo
  end subroutine calculateRhoW


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
!          write(*,*) "rhoe", thisOctal%rhoe(subcell), thisOctal%rho(subcell), thisOctal%energy(subcell)
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


  recursive subroutine unsetGhosts(thisOctal)
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
                call unsetGhosts(child)
                exit
             end if
          end do
       else
          thisOctal%ghostCell(subcell) = .false.
          thisOctal%feederCell(subcell) = .false.
      endif
    enddo
  end subroutine unsetGhosts

  recursive subroutine feederCellCheck(thisOctal)
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
                call feederCellCheck(child)
                exit
             end if
          end do
       else
          if (thisOctal%feederCell(subcell).and.thisOctal%ghostCell(subcell)) then
             thisOctal%rho(subcell) = 0.8d0
          endif
      endif
    enddo
  end subroutine feederCellCheck

  recursive subroutine boundaryCondCheck(thisOctal)
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
                call boundaryCondCheck(child)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          write(*,*) "looking at ", thisOctal%boundaryCondition(subcell), " subcell ", subcell
          if (thisOctal%boundaryCondition(subcell) /= 1) then
             write(*,*) "boundary cond check failure: ",thisOctal%boundaryCondition(subcell), &
                  thisOctal%chiline(subcell), " subcell " , subcell
             write(*,*) " depth ",thisOctal%nDepth
             write(*,*) "rank ",myrankglobal, " mpithread ", thisOctal%mpiThread(subcell)
             write(*,*) "rho ", thisOctal%rho(subcell)
          endif
      endif
    enddo
  end subroutine boundaryCondCheck



  recursive subroutine imposeBoundary(thisOctal)
    type(octal), pointer   :: thisOctal, bOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, bSubcell
    logical :: firsttime
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
                call imposeBoundary(child)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%ghostCell(subcell)) then
             select case(thisOctal%boundaryCondition(subcell))
                case(1)
                   if (.not.associated(thisOctal%tempStorage)) allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:5))
                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocal(locator, bOctal, bSubcell)

!                   if (bOctal%ghostCell(bSubcell)) then
!                      write(*,*) "Error selecting boundary partner!!!"
!                   endif

                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)
                   Thisoctal%tempstorage(subcell,1) = bOctal%rho(bSubcell)
                   thisOctal%tempStorage(subcell,2) = bOctal%rhoE(bSubcell)

! NB confusion regarding 2d being x,z rather than x,y

                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                      endif
                      if ((abs(dir%x) > 0.5d0).and.abs(dir%z) > 0.5d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                      endif
                   else if (thisOctal%threed) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (abs(dir%y) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = -bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                      endif
                      if ((abs(dir%x) > 0.2d0).and.abs(dir%z) > 0.2d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = -bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                      endif
                   endif

                case(2)
                   if (.not.associated(thisOctal%tempStorage)) allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:5))
!                   locator = thisOctal%boundaryPartner(subcell)
!                   bOctal => thisOctal
!                   call findSubcellLocal(locator, bOctal, bSubcell)
!
!                   thisoctal%tempstorage(subcell,1) = bOctal%rho(bSubcell)
!                   thisOctal%tempStorage(subcell,2) = bOctal%rhoE(bSubcell)
!                   thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
!                   thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
!                   thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)

                   thisoctal%tempstorage(subcell,1) = thisOctal%rho(subcell)
                   thisOctal%tempStorage(subcell,2) = thisOctal%rhoE(subcell)
                   thisOctal%tempStorage(subcell,3) = thisOctal%rhou(subcell)
                   thisOctal%tempStorage(subcell,4) = thisOctal%rhov(subcell)
                   thisOctal%tempStorage(subcell,5) = thisOctal%rhow(subcell)



                case(3)
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
                   write(*,*) "Unrecognised boundary condition in impose boundary: ", thisOctal%boundaryCondition(subcell)
             end select
          endif
      endif
    enddo
  end subroutine imposeBoundary



  recursive subroutine setupGhostCells(thisOctal, grid, flag)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal, tempOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell, tempSubcell
    type(OCTALVECTOR) :: locator, rVec
    integer :: nProbes, iProbe
    type(OCTALVECTOR) :: probe(6), direction
    real(double) :: dx
    integer :: nProbeOutside
    logical :: corner
    logical, optional :: flag
    integer :: myRank, ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupGhostCells(child,  grid, flag)
                exit
             end if
          end do
       else
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle


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
                  (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)
             if (.not.inOctal(grid%octreeRoot, locator).or. &
                  (locator%x < (grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize))) &
                  nProbeOutside = nProbeOutside + 1
          enddo
          corner=.false.
          if (thisOctal%twoD.and.(nProbeOutside > 1)) corner = .true.
          if (thisOctal%threeD.and.(nProbeOutside > 2)) corner = .true.
          if (.not.corner) then
             do iProbe = 1, nProbes
                locator = rVec + &
                     (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)
                ! this is fudged here because AMR mod assumes that the grid is one-d spherical....
                if (.not.inOctal(grid%octreeRoot, locator).or. &
                     (locator%x < (grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize))) then  ! this is a boundary
                   ! setup the the two ghosts near the boundary
                   
                   thisOctal%ghostCell(subcell) = .true.
                   if (present(flag)) then
                      thisOCtal%rho(subcell) = 2.5d0
                   endif
                   locator = rVec - &
                        (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)
                   neighbourOctal => thisOctal
                   call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                   neighbourOctal%ghostCell(neighbourSubcell) = .true.

                   if (present(flag)) neighbourOctal%rho(neighbourSubcell) = 2.5d0
                   
                   ! now a case to determine the boundary cell relations
                   select case (thisOctal%boundaryCondition(subcell))
                   case(1)
                      dx = thisOctal%subcellSize
                      thisOctal%ghostCell(subcell) = .true.

                      call locatorToNeighbour(grid, thisOctal, subcell, (-1.d0)*probe(iProbe), 3, locator)
                      thisOctal%boundaryPartner(subcell) = locator

                      tempOctal => thisOctal
                      tempSubcell = 1
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      tempOctal%feederCell(tempsubcell) = .true.
                      
                      if (present(flag)) then
                         tempOctal%rho(tempSubcell) = 3.d0
                      endif


                      call locatorToNeighbour(grid, neighbourOctal, neighboursubcell, (-1.d0)*probe(iProbe), 1, locator)
                      neighbourOctal%boundaryPartner(neighboursubcell) = locator
                      tempOctal => thisOctal
                      tempSubcell = 1
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      tempOctal%feederCell(tempsubcell) = .true.
                      
                      if (present(flag)) then
                         tempOctal%rho(tempSubcell) = 3.d0
                      endif

                   case(2)
                      dx = thisOctal%subcellSize
                      thisOctal%ghostCell(subcell) = .true.

                      locator = subcellCentre(thisOctal, subcell) - &
                           (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)*probe(iProbe)
                      thisOctal%boundaryPartner(subcell) = locator

                      tempSubcell = 1
                      tempOctal => thisOctal
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      if (tempOctal%mpiThread(tempSubcell) == myRankGlobal) then
                         write(*,*) "locator problem ",myRankGlobal
                         stop
                      endif
                      locator = subcellCentre(neighbourOctal, neighboursubcell) - &
                           (grid%octreeRoot%subcellSize*2.d0-4.d0*neighbourOctal%subcellSize)*probe(iProbe)
                      neighbourOctal%boundaryPartner(neighboursubcell) = locator
                      tempSubcell = 1
                      tempOctal => thisOctal
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      if (tempOctal%mpiThread(tempSubcell) == myRankGlobal) then
                         write(*,*) "locator problem ",myRankGlobal
                         stop
                      endif



                   case DEFAULT
                      write(*,*) "unknown boundary condition in setupghostcells1: ",thisOCtal%boundaryCondition(subcell)
                      write(*,*) "rank ",myRankGlobal, " mpi thread ", thisOctal%mpithread(subcell)
                      write(*,*) "subcell ",subcell
                      write(*,*) "all ", thisOctal%boundaryCondition(subcell)
                   end select
                endif
             enddo
          else ! this is a corner
             direction = subcellCentre(thisOctal, subcell) - grid%octreeRoot%centre
             call normalize(direction)
             if (thisOctal%twoD) then
                direction = sqrt(2.d0)*direction
             else
                direction = sqrt(3.d0)*direction
             endif
                   
             thisOctal%ghostCell(subcell) = .true.
             if (present(flag))  thisOCtal%rho(subcell) = 1.5d0
             locator = rVec - &
                  (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*direction
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             neighbourOctal%ghostCell(neighbourSubcell) = .true.
             if (present(flag)) neighbourOctal%rho(neighbourSubcell) = 2.2d0
                   
                   
             ! now a case to determine the boundary cell relations
             select case (thisOctal%boundaryCondition(subcell))
             case(1)
                 dx = sqrt(2.d0)*thisOctal%subcellSize
                thisOctal%ghostCell(subcell) = .true.
                
                call locatorToNeighbour(grid, thisOctal, subcell, (-1.d0)*direction, 3, locator)
                thisOctal%boundaryPartner(subcell) = locator
                      tempOctal => thisOctal
                      tempSubcell = 1
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      tempOctal%feederCell(tempsubcell) = .true.
                      
                      if (present(flag)) then
                         tempOctal%rho(tempSubcell) = 0.8d0
                      endif
                
                call locatorToNeighbour(grid, neighbourOctal, neighboursubcell, (-1.d0)*direction, 1, locator)
                neighbourOctal%boundaryPartner(neighboursubcell) = locator
                      tempOctal => thisOctal
                      tempSubcell = 1
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      tempOctal%feederCell(tempsubcell) = .true.
                      
                      if (present(flag)) then
                         tempOctal%rho(tempSubcell) = 0.8d0
                      endif

                case(2)
                   dx = thisOctal%subcellSize
                   thisOctal%ghostCell(subcell) = .true.
                   
                   locator = subcellCentre(thisOctal, subcell) - &
                        (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)*direction
                   thisOctal%boundaryPartner(subcell) = locator
                      tempSubcell = 1
                      tempOctal => thisOctal
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      if (tempOctal%mpiThread(tempSubcell) == myRankGlobal) then
                         write(*,*) "locator problem ",myRankGlobal
                         stop
                      endif

                   
                   locator = subcellCentre(neighbourOctal, neighboursubcell) - &
                        (grid%octreeRoot%subcellSize*2.d0-4.d0*neighbourOctal%subcellSize)*direction
                   neighbourOctal%boundaryPartner(neighboursubcell) = locator

                      tempSubcell = 1
                      tempOctal => thisOctal
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      if (tempOctal%mpiThread(tempSubcell) == myRankGlobal) then
                         write(*,*) "locator problem ",myRankGlobal
                         stop
                      endif
                   
                case DEFAULT
                   write(*,*) "Unknown boundary condition in setupghostcells2: ",thisOctal%boundaryCondition(subcell)
             end select


          endif
      endif
    enddo
  end subroutine setupGhostCells

  recursive subroutine setupEdges(thisOctal, grid)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(OCTALVECTOR) :: locator, rVec
    integer :: nProbes, iProbe
    type(OCTALVECTOR) :: probe(6), direction
    real(double) :: dx
    integer :: nProbeOutside
    logical :: corner
    integer :: myRank, ierr


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupEdges(child,  grid)
                exit
             end if
          end do
       else
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          thisOctal%edgeCell(subcell) = .false.
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
                  (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)
             if (.not.inOctal(grid%octreeRoot, locator).or. &
                  (locator%x < (grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize))) &
                  nProbeOutside = nProbeOutside + 1
          enddo
          if (nProbeOutside >= 1) then
             thisOctal%edgeCell(subcell) = .true.
          endif
       endif
    enddo
  end subroutine setupEdges




  recursive subroutine refineGridGeneric(thisOctal, grid, criterion, inherit)
    use input_variables, only : maxDepthAMR
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    type(OCTALVECTOR) :: rVec, corner
    integer :: i, j, subcell, iStream
    logical :: inherit
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
                call refineGridGeneric(child, grid, criterion, inherit)
                exit children
             end if
          end do children
       else
          if (.not.thisOctal%ghostCell(subcell)) then
             call  splitCondition(thisOctal, grid, subcell, criterion, split)
             if (split.and.(thisOctal%nDepth < maxDepthAMR)) then
                call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                     inherit=inherit, interp=.false.)
!                subcell = subcell - 1
             endif
          endif
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
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux
    real(double) :: grad, maxGradient
    
    split = .false.
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
       probe(1) = OCTALVECTOR( 0.d0, 0.d0, +1.d0)
       probe(2) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
       probe(3) = OCTALVECTOR(-1.d0, 0.d0,  0.d0)
       probe(4) = OCTALVECTOR(+1.d0, 0.d0,  0.d0)
       probe(5) = OCTALVECTOR( 0.d0, 1.d0,  0.d0)
       probe(6) = OCTALVECTOR( 0.d0,-1.d0,  0.d0)
    endif

    maxGradient = 1.d-30
    do i = 1, nProbes
       locator = subcellCentre(thisOctal, subcell) + &
            (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * probe(i)
       if (inOctal(grid%octreeRoot, locator)) then
          neighbourOctal => thisOctal
          call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
          call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(i), q, rho, rhoe, &
               rhou, rhov, rhow, x, qnext, pressure, flux)
          grad = abs((thisOctal%rho(subcell)-rho) / &
               thisOctal%rho(subcell))
          maxGradient = max(grad, maxGradient)
          grad = abs((thisOctal%rhoe(subcell)-rhoe) / &
               thisOctal%rhoe(subcell))
          maxGradient = max(grad, maxGradient)
          if (.not.thisOctal%ghostCell(subcell)) then
             grad = abs((thisOctal%rhou(subcell)-rhou) / &
                  thisOctal%rhou(subcell))
             maxGradient = max(grad, maxGradient)
             grad = abs((thisOctal%rhov(subcell)-rhov) / &
                  thisOctal%rhov(subcell))
             maxGradient = max(grad, maxGradient)
             grad = abs((thisOctal%rhow(subcell)-rhow) / &
                  thisOctal%rhow(subcell))
             maxGradient = max(grad, maxGradient)
          endif
       endif
    end do
    if (maxGradient > 0.2d0) split = .true.
  end subroutine splitCondition



  recursive subroutine evenUpGrid(thisOctal, grid,  converged, inherit)

    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    !
    integer :: subcell, i, ilambda
    logical :: converged, converged_tmp
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat, loc
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    logical, optional :: inherit
    integer :: myRank, ierr
    converged = .true.
    converged_tmp=.true.


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call evenUpGrid(child, grid,  converged, inherit)
                if (.not.converged) return
                if (.not.converged) converged_tmp = converged
                exit
             end if
          end do
          if (.not.converged_tmp) converged=converged_tmp
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
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
             dirVec(1) = OCTALVECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
          endif

             do j = 1, nDir
                octVec = centre + r * dirvec(j)
                if (inOctal(grid%octreeRoot, octVec)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)

                   if (neighbourOctal%mpiThread(neighbourSubcell) == myRank) then

                      if ((neighbourOctal%nDepth-thisOctal%nDepth) > 1) then
                         call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                              inherit=inherit, interp=.false.)
                         converged = .false.
                         exit
                      endif
                      
                      if ((thisOctal%nDepth-neighbourOctal%nDepth) > 1) then
                         call addNewChild(neighbourOctal,neighboursubcell,grid,adjustGridInfo=.TRUE., &
                              inherit=inherit, interp=.false.)
                         converged = .false.
                         exit
                      endif
                      
                      if (thisOctal%edgeCell(subcell).and.(.not.neighbourOctal%edgeCell(neighboursubcell))  &
                           .and.(thisOctal%nDepth < neighbourOctal%nDepth)) then
                         call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                              inherit=inherit, interp=.false.)
                         converged = .false.
                         exit
                      endif
                   endif
                end if
             enddo
          if (.not.converged) exit
       endif
    end do

  end subroutine evenUpGrid

  recursive subroutine refineGridGeneric2(thisOctal, grid,  gamma, converged, iEquationOfState, inherit)
    use input_variables, only : maxDepthAMR
    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    !
    integer :: subcell, i, ilambda
    logical :: converged, converged_tmp
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat, loc, locator
    logical :: split
    integer :: neighbourSubcell, j, nDir
    real(double) :: r, grad, maxGradient
    logical, optional :: inherit
    real(double), parameter :: limit = 0.01d0
    real(double) :: gamma
    real(double) :: cs, rhocs
    integer :: myRank, ierr
    integer :: iEquationOfState
    converged = .true.
    converged_tmp=.true.


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call refineGridGeneric2(child, grid,  gamma, converged, iEquationOfState, inherit)
                if (.not.converged) return
                if (.not.converged) converged_tmp = converged
                exit
             end if
          end do
          if (.not.converged_tmp) converged=converged_tmp
       else


!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
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
             dirVec(1) = OCTALVECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
          endif


          do i = 1, nDir
             maxGradient = 1.d-30
             locator = subcellCentre(thisOctal, subcell) + &
                  (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
             if (inOctal(grid%octreeRoot, locator)) then
                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                
                grad = abs((thisOctal%rho(subcell)-neighbourOctal%rho(neighbourSubcell)) / &
                     thisOctal%rho(subcell))
                maxGradient = max(grad, maxGradient)

!                if (grad > limit) write(*,*) "rho",grad

                if (thisOctal%rhoe(subcell) /= 0.d0) then
                   grad = abs((thisOctal%rhoe(subcell)-neighbourOctal%rhoe(neighbourSubcell)) / &
                        thisOctal%rhoe(subcell))
                   maxGradient = max(grad, maxGradient)
                endif
!                if (grad > limit) write(*,*) "rhoe",grad, thisOctal%rhoe(subcell), neighbourOctal%rhoe(neighboursubcell),&
!                     thisOctal%mpiThread(subcell), neighbourOctal%mpithread(neighboursubcell)
                if (.not.thisOctal%ghostCell(subcell)) then
                   cs = soundSpeed(thisOctal, subcell, gamma, iEquationOfState)
                   rhocs = thisOctal%rho(subcell) * cs

                   if (rhocs /= 0.d0) then
                      grad = abs((thisOctal%rhou(subcell)-neighbourOctal%rhou(neighbourSubcell)) / &
                           rhocs)
                      maxGradient = max(grad, maxGradient)
!                      if (grad > limit) write(*,*) "rhou",grad
                      
                      grad = abs((thisOctal%rhov(subcell)-neighbourOctal%rhov(neighbourSubcell)) / &
                        rhocs)
                      maxGradient = max(grad, maxGradient)
!                      if (grad > limit) write(*,*) "rhov",grad
                      
                      grad = abs((thisOctal%rhow(subcell)-neighbourOctal%rhow(neighbourSubcell)) / &
                           rhocs)
                      maxGradient = max(grad, maxGradient)
!                      if (grad > limit) write(*,*) "rhow",grad
                   endif
                   if (maxGradient > limit) then
                      split = .true.
                   else
                      split = .false.
                   endif


                   if (split) then
                      if ((neighbourOctal%nDepth >= thisOctal%nDepth).and.(thisOctal%nDepth < maxDepthAMR)) then
                         call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                              inherit=inherit, interp=.false.)
                         converged = .false.
                         exit
                      endif

!                      if ((thisOctal%nDepth >= neighbourOctal%nDepth).and.(neighbourOctal%nDepth < maxDepth)) then
!                         call addNewChild(neighbourOctal,neighboursubcell,grid,adjustGridInfo=.TRUE., &
!                              inherit=inherit, interp=.false.)
!                         converged = .false.
!                         exit
!                      endif

                   endif
                endif
             endif
          end do
          if (.not.converged) exit
       endif
    end do

  end subroutine refineGridGeneric2

  recursive subroutine refineEdges(thisOctal, grid,  converged, inherit)

    use input_variables, only : maxDepthAMR
    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    integer :: subcell, i, ilambda
    logical :: converged, converged_tmp
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat, rvec
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    logical, optional :: inherit
    integer :: myRank, ierr
    converged = .true.
    converged_tmp=.true.

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call refineEdges(child, grid,  converged, inherit)
                if (.not.converged) return
                if (.not.converged) converged_tmp = converged
                exit
             end if
          end do
          if (.not.converged_tmp) converged=converged_tmp
       else
          
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle


          if ((thisOctal%edgeCell(subcell)) &
          .and.thisOctal%nDepth<maxDepthAMR) then
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=inherit, interp=.false.)
             converged = .false.
          endif
          if (.not.converged) exit
       endif
    end do

  end subroutine refineEdges

  recursive subroutine refineFeeders(thisOctal, grid,  converged, inherit)

    use input_variables, only : maxDepthAMR
    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    !
    integer :: subcell, i, ilambda
    logical :: converged, converged_tmp
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat, rvec
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    logical, optional :: inherit
    integer :: myRank, ierr
    converged = .true.
    converged_tmp=.true.

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren



       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call refineFeeders(child, grid,  converged, inherit)
                if (.not.converged) converged_tmp = converged
                exit
             end if
          end do
          if (.not.converged_tmp) converged=converged_tmp
       else
          
!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle


          if ((thisOctal%feederCell(subcell)) &
          .and.thisOctal%nDepth<maxDepthAMR) then
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=inherit, interp=.false.)
             converged = .false.
          endif
          if (.not.converged) exit
       endif
    end do

  end subroutine refineFeeders



  subroutine locatorToNeighbour(grid, thisOctal, subcell, direction, ncells, locator)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, tempOctal
    integer :: subcell, tempSubcell, i, nCells
    type(OCTALVECTOR) :: direction, locator, rVec
    tempOctal => thisOctal
    tempSubcell = subcell
    do i = 1, nCells
       rVec = subcellCentre(tempOctal, tempSubcell) + &
            (tempOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell) * direction 
       call findSubcellLocal(rVec, tempOctal, tempSubcell)
    enddo
    locator = rVec
  end subroutine locatorToNeighbour
       
  recursive subroutine unrefineCells(thisOctal, grid, gamma, iEquationOfState)
    use input_variables, only : minDepthAMR
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: ilambda
    real(double) :: gamma
    integer :: subcell, i, j
    logical :: unrefine, ok, converged
    integer :: nTau
    integer :: nVals, nc
    real(double) :: rhow(8), rhov(8), rhou(8), rho(8), rhoe(8), fac, limit
    real(double) :: cs(8)
    real(double) :: rhocs, rhomean, rhoemean
    logical :: refinedLastTime, ghostCell
    integer :: iEquationOfState
    limit  = 0.01d0

    unrefine = .true.
    refinedLastTime = .false.
    ghostCell = .false.
    nc = 0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unrefineCells(child, grid,  gamma, iEquationOfState)
                exit
             end if
          end do
       else
          nc = nc + 1
          rho(nc) = thisOctal%rho(subcell)
          rhoe(nc) = thisOctal%rhoe(subcell)
          rhou(nc) = thisOctal%rhou(subcell)
          rhov(nc) = thisOctal%rhov(subcell)
          rhow(nc) = thisOctal%rhow(subcell)

          cs(nc) = soundSpeed(thisOctal, subcell, gamma, iEquationOfState)
          if (thisOctal%ghostCell(subcell)) ghostCell=.true.
!          if (thisOctal%refinedLastTime(subcell)) refinedLastTime = .true.
       endif
    enddo


    unrefine = .false.
    if ((nc > 1).and.(.not.ghostCell)) then

       rhomean = SUM(rho(1:nc))/dble(nc)
       rhoemean = SUM(rhoe(1:nc))/dble(nc)
       rhocs = rhomean*(SUM(cs(1:nc))/dble(nc))
       unrefine = .true.
       
       fac = sigma(rho, nc) / rhoMean
       if (fac > limit) then
          unrefine = .false.
!          write(*,*) "rho",fac
       endif
      
       fac = sigma(rhoe, nc) / rhoeMean
       if (fac > limit) then
          unrefine = .false.
!          write(*,*) "rhoe",fac
       endif
       
       fac = sigma(rhou, nc) / rhocs
       if (fac > limit) then
          unrefine = .false.
!          write(*,*) "rhou",fac
       endif
       
       fac = sigma(rhov, nc) / rhocs
       if (fac > limit) then
          unrefine = .false.
!          write(*,*) "rhov",fac
       endif
       
       
       if (thisOctal%nDepth <= minDepthAMR) unrefine = .false.
    endif

    if ((thisOctal%nChildren == 0).and.unrefine) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
    endif
    
  end subroutine unrefineCells
 

  function sigma(x, n)
    real(double) :: x(:), sigma, mean
    integer :: n, i

!    mean = SUM(x(1:n))/dble(n)
!    do i = 1, n
!       sigma = sigma + (x(i) - mean)**2
!    enddo
!    sigma = sqrt(sigma / dble(n))
    sigma = maxval(abs(x(1:n)))-minval(abs(x(1:n)))
  end function sigma

  subroutine evenUpGridMPI(grid, inheritFlag)

    type(GRIDTYPE) :: grid
    logical :: gridConverged, inheritFlag
    integer :: myRank, nThreads
    include 'mpif.h'  
    integer :: ierr
    logical :: globalConverged, localChanged(64)
    type(OCTALVECTOR) :: locs(10000), eLocs(10000)
    integer :: nLocs(64), tempNlocs(10000)
    integer :: thread(10000), nLocsGlobal,i, depth(10000)
    integer :: nGroup, group(1000)

    real(double) :: temp(10000,4),tempsent(4)
    integer :: nTemp(1), nSent(1), eDepth(10000)
    integer :: iThread, nExternalLocs
    integer :: iter
    integer, parameter :: tag = 1
    logical :: globalChanged(64)
    integer :: status(MPI_STATUS_SIZE)
    character(len=20) :: plotfile
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)
    ! first even up the local part of the grid


    globalConverged = .false.
    iter = 0
    do while(.not.globalConverged)
       localChanged = .false.

       i = 1
       do
          gridConverged = .true.
          call setupEdges(grid%octreeRoot, grid)
          call refineEdges(grid%octreeRoot, grid,  gridconverged, inherit=inheritFlag)

!          write(plotfile,'(a,i4.4,a)') "grid",i,".png/png"
!          call plotGridMPI(grid, plotfile, "x-z", "rho", valueMinFlag = 1.5e-22, valueMaxFlag = 4.e-22, plotgrid=.true.)
!          i = i  + 1
          if (gridConverged) exit
       end do

       call unsetGhosts(grid%octreeRoot)
       call setupGhostCells(grid%octreeRoot, grid)

       do
          gridConverged = .true.
          call evenUpGrid(grid%octreeRoot, grid,  gridconverged, inherit=inheritFlag)
          call setupEdges(grid%octreeRoot, grid)
          call unsetGhosts(grid%octreeRoot)
          call setupGhostCells(grid%octreeRoot, grid)
          if (gridConverged) exit
       end do
       nLocs = 0
       call locatorsToExternalCells(grid%octreeRoot, grid, nLocs(myRank), locs, thread, depth)
       call MPI_ALLREDUCE(nLocs, tempnLocs, nThreadsglobal-1, MPI_INTEGER, MPI_SUM,amrCOMMUNICATOR, ierr)
       nLocsGlobal = SUM(tempnLocs)


       do iThread = 1, nThreads-1
          nTemp(1) = 0
          if (iThread /= myRank) then
             do i = 1 , nLocs(myRank)
                if (thread(i) == iThread) then
                   nTemp(1) = nTemp(1) + 1
                   temp(nTemp(1),1) = locs(i)%x
                   temp(nTemp(1),2) = locs(i)%y
                   temp(nTemp(1),3) = locs(i)%z
                   temp(nTemp(1),4) = dble(depth(i))
                endif
             enddo
             call mpi_send(nTemp, 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, ierr)
             if (nTemp(1) > 0) then
                do i = 1, nTemp(1)
                   call mpi_send(temp(i,:), 4, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
                enddo
             endif
          endif
       enddo
!       write(*,*) myrank, " entering receive loop"
       do iThread = 1, nThreads - 1
          if (iThread /= myRank) then
!                          write(*,*) "rank ",myRank," receiving message from ",ithread
             call mpi_recv(nSent, 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
!                          write(*,*) "received ",nSent
             if (nSent(1) > 0) then
                do i = 1, nSent(1)
                   call mpi_recv(tempsent, 4, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, status, ierr)
                   eLocs(i)%x = tempsent(1)
                   eLocs(i)%y = tempsent(2)
                   eLocs(i)%z = tempsent(3)
                   eDepth(i) = nint(tempsent(4))
                enddo
             endif
             nExternalLocs = nSent(1)
             
             call setAllUnchanged(grid%octreeRoot)
             
             do i = 1, nExternalLocs
                call splitAtLocator(grid, elocs(i), edepth(i), localChanged(myRank), inherit=inheritflag)
             enddo
!             if (localChanged(myrank)) then
                do
                   gridConverged = .true.
                   call evenUpGrid(grid%octreeRoot, grid,  gridconverged, inherit=inheritFlag)
                   call setupEdges(grid%octreeRoot, grid)
                   call unsetGhosts(grid%octreeRoot)
                   call setupGhostCells(grid%octreeRoot, grid)!, flag=.true.)
                   if (gridConverged) exit
                end do
!             endif
          endif
       enddo
       call MPI_ALLREDUCE(localChanged, globalChanged, nThreadsGlobal-1, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)

       if (ANY(globalChanged(1:nThreadsGlobal-1))) then
          globalConverged = .false.
       else
          globalConverged = .true.
       endif
       iter = iter + 1
       
    end do
  end subroutine evenUpGridMPI

  subroutine splitAtLocator(grid, locator, depth,  localchanged, inherit)
    use input_variables, only :  maxDepthAMR
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    type(OCTALVECTOR) :: locator
    logical, optional :: inherit
    integer :: depth
    integer :: subcell
    integer :: myRank, ierr
    logical :: localChanged 

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    call findSubcellTD(locator, grid%octreeRoot, thisOctal,subcell)
    
    if (myRank /= thisOctal%mpiThread(subcell)) then
       write(*,*) "There is a bug somewhere: ",myrank
       stop
    endif
!    WRITE(*,*) "ATTEMPTNG TO SPLIT",depth,thisOctal%nDepth, maxdepth
    if (((depth-thisOctal%nDepth) > 1).and.(thisOctal%nDepth < maxDepthAMR)) then
       call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
            inherit=inherit, interp=.false.)
       localChanged = .true.
!       write(*,*) "split at locator"
    endif

  end subroutine splitAtLocator

  recursive subroutine locatorsToExternalCells(thisOctal, grid, nLocs, loc, thread, depth)

    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    !
    integer :: subcell, i, ilambda
    logical :: converged, converged_tmp
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat, loc(:)
    integer :: thread(:), depth(:)
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    integer :: myRank, ierr
    integer :: nLocs
    converged = .true.
    converged_tmp=.true.


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call locatorsToExternalCells(child, grid, nLocs, loc, thread, depth)
                exit
             end if
          end do
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
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
             dirVec(1) = OCTALVECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
          endif

             do j = 1, nDir
                octVec = centre + r * dirvec(j)
                if (inOctal(grid%octreeRoot, octVec)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)

                   if (neighbourOctal%mpiThread(neighboursubcell) /= myRank) then
                      nLocs = nLocs + 1
                      loc(nLocs) = octVec
                      thread(nLocs) = neighbourOctal%mpiThread(neighbourSubcell)
                      depth(nLocs) = thisOctal%nDepth
                   endif

                end if
             enddo
       endif
    end do

  end subroutine locatorsToExternalCells



  recursive subroutine determineDependentThreads(thisOctal, grid, nDependent,  dependentthread)

    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    !
    integer :: subcell, i, ilambda
    logical :: converged, converged_tmp
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat
    integer :: dependentThread(:)
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    integer :: myRank, ierr
    integer :: nDependent
    logical :: alreadyInList
    integer :: iDep
    converged = .true.
    converged_tmp=.true.


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call determineDependentThreads(child, grid, nDependent, dependentThread)
                exit
             end if
          end do
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
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
             dirVec(1) = OCTALVECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
          endif

             do j = 1, nDir
                octVec = centre + r * dirvec(j)
                if (inOctal(grid%octreeRoot, octVec)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)

                   if (neighbourOctal%mpiThread(neighboursubcell) /= myRank) then
                      alreadyInList = .false.
                      do iDep = 1, nDependent
                         if (neighbourOctal%mpiThread(neighbourSubcell) == dependentThread(iDep)) then
                            alreadyinList = .true.
                         endif
                      enddo
                      if (.not.alreadyInList) then
                         nDependent = nDependent + 1
                         dependentThread(nDependent) = neighbourOctal%mpiThread(neighbourSubcell)
                      endif
                   endif

                end if
             enddo
       endif
    end do

  end subroutine determineDependentThreads





#endif
    
end module hydrodynamics_mod
