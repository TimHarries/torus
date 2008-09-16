! hydrodynamics module added by TJH on 18th June 2007

module hydrodynamics_mod
#ifdef MPI

  use input_variables
  use kind_mod
  use constants_mod
  use amr_mod
  use grid_mod
  use source_mod
  use timing
  use mpi_amr_mod
  use vtk_mod

  implicit none
  real(double) :: gridDistanceScale = 1.d10
!  real(double) :: gridDistanceScale = 1.d0

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
           if (thisOctal%ghostCell(subcell)) thisOctal%phiLimit(subcell) = 1.d0
!           thisOctal%PhiLimit(subcell) = 1.d0
      endif
    enddo
  end subroutine fluxLimiter


  recursive subroutine updateDensityTree(thisOctal)
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal, parentOctal, testOctal
    type(octal), pointer  :: child 
    integer :: i, n, m

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    if (thisOctal%nChildren > 0) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call updateDensityTree(child)
       end do
    else

          
       testOctal => thisOctal
          
       do while(testOctal%nDepth >= 2)
             
          parentOctal => testOctal%parent
          m = testOctal%parentSubcell
          if (testOctal%twoD) then
             n = 4
          else
             n = 8
          endif
          parentOctal%rho(m) = SUM(testOctal%rho(1:n))/dble(n)
          parentOctal%phi_i(m) = SUM(testOctal%phi_i(1:n))/dble(n)
          parentOctal%edgeCell(m) = ANY(testOctal%edgeCell(1:n))
          testOctal => parentOctal
       enddo
    endif
  end subroutine updateDensityTree

  recursive subroutine updatePhiTree(thisOctal, nDepth)
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, n, nDepth

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    if (thisOctal%nDepth == (nDepth+1)) then

       if (thisOctal%twoD) then
          n = 4
       else
          n = 8
       endif

       where(.not.thisOctal%edgeCell)
          thisOctal%phi_i = thisOctal%parent%phi_i(thisOctal%parentSubcell)
       end where
!       write(*,*) "phi ",thisOctal%phi_i(1:8)
    else

       do subcell = 1, thisOctal%maxChildren
          if (thisOctal%hasChild(subcell)) then
             ! find the child
             do i = 1, thisOctal%nChildren, 1
                if (thisOctal%indexChild(i) == subcell) then
                   child => thisOctal%child(i)
                   call updatephiTree(child, nDepth)
                   exit
                end if
             end do
          endif
       enddo
    endif

  end subroutine updatePhiTree


  recursive subroutine constructFlux(thisOctal, dt)
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx

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

          dx = (thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell))
          dx = thisOctal%subcellSize * gridDistanceScale

          if (.not.thisOctal%edgeCell(subcell)) then
             thisOctal%flux_i(subcell) = thisOctal%flux_i(subcell) + &
                  0.5d0 * abs(thisOctal%u_interface(subcell)) * &
                  (1.d0 - abs(thisOctal%u_interface(subcell) * dt / &
                   dx)) * &
                  thisOctal%phiLimit(subcell) * (thisOctal%q_i(subcell) - thisOctal%q_i_minus_1(subcell))
          endif
!          if (thisOctal%flux_i(subcell) > 1.d-10) write(*,*) "flux ",thisOctal%flux_i(subcell),thisOctal%u_interface(subcell), &
!               thisOctal%rho(subcell), thisOctal%rhou(subcell)
       endif
    enddo
  end subroutine constructFlux

  recursive subroutine updateCellQ(thisOctal, dt)
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx !, x_plus_half, x_minus_half

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
             dx = thisOctal%subcellSize * gridDistanceScale

!             x_minus_half = 0.5d0*(thisOctal%x_i(subcell)+thisOctal%x_i_minus_1(subcell))
!             x_plus_half = 0.5d0*(thisOctal%x_i(subcell)+thisOctal%x_i_plus_1(subcell))
!             dx = x_plus_half - x_minus_half

!             thisOctal%q_i(subcell) = thisOctal%q_i(subcell) - dt * &
!                  (thisOctal%flux_i_plus_1(subcell) - thisOctal%flux_i(subcell)) / &
!                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i(subcell))


!             dx = thisOctal%subcellSize * gridDistanceScale

             thisOctal%q_i(subcell) = thisOctal%q_i(subcell) - dt * &
                  (thisOctal%flux_i_plus_1(subcell) - thisOctal%flux_i(subcell)) / dx


          endif
       endif
    enddo
  end subroutine updateCellQ

  recursive subroutine synchronizeFluxes(thisOctal, dt, iDepth)
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iDepth
    real(double) :: dt !, dx

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call synchronizeFluxes(child, dt, iDepth)
                exit
             end if
          end do
       else

!          if (thisOctal%mpiThread(subcell) /= myRank) cycle
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (thisOctal%nDepth /= iDepth) cycle

          if (.not.thisOctal%ghostCell(subcell)) then


!             thisOctal%q_i(subcell) = thisOctal%q_i(subcell) - dt * &
!                  (thisOctal%flux_i_plus_1(subcell) - thisOctal%flux_i(subcell)) / dx


          endif
       endif
    enddo
  end subroutine synchronizeFluxes

  recursive subroutine setupX(thisOctal, grid, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: direction
  
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
          thisOctal%x_i(subcell) = (subcellCentre(thisOctal, subcell) .dot. direction) * gridDistanceScale
       endif
    enddo
  end subroutine setupX


  recursive subroutine dumpx(thisOctal, n)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, n
    type(VECTOR) :: locator
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call dumpX(child, n)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          if (myRankGlobal==1) then
             locator = subcellCentre(thisOctal,subcell)
             n = n + 1
             write(*,*) n , locator%x
          endif
       endif
    enddo
  end subroutine dumpx

  recursive subroutine setupQX(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi
    integer :: subcell, i, neighbourSubcell
    integer :: nd
    type(VECTOR) :: direction, locator, reverseDirection
  
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
          if (.not.thisOctal%edgeCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

             thisOctal%x_i_plus_1(subcell) = x
             thisOctal%q_i_plus_1(subcell) = q

             reverseDirection = (-1.d0) * direction
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
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

  recursive subroutine setupRhoPhi(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: direction, locator, reverseDirection
    integer :: nd
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupRhoPhi(child, grid, direction)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%edgeCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

             thisOctal%rho_i_plus_1(subcell) = rho
             thisOctal%phi_i_plus_1(subcell) = phi

             reverseDirection = (-1.d0) * direction
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisOctal%rho_i_minus_1(subcell) = rho
             thisOctal%phi_i_minus_1(subcell) = phi



          endif
       endif
    enddo
  end subroutine setupRhoPhi

  recursive subroutine setupUi(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, weight
    integer :: nd
  
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

          if (.not.thisOctal%edgeCell(subcell)) then !xxx changed from ghostcell
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

             rho_i_minus_1 = rho
             rhou_i_minus_1 = rhou

!             if (thisOctal%nDepth == nd) then
                weight = 0.5d0
!             else if (thisOctal%nDepth > nd) then ! fine to coarse
!                weight  = 0.666666666d0
!             else
!                weight  = 0.333333333d0 ! coarse to fine
!             endif


             thisOctal%u_interface(subcell) = &
                  weight*thisOctal%rhou(subcell)/thisOctal%rho(subcell) + (1.d0-weight)*rhou_i_minus_1/rho_i_minus_1

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
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, weight
    integer :: nd

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

          if (.not.thisOctal%edgeCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

                rho_i_minus_1 = rho
                rhou_i_minus_1 = rhov

!             if (thisOctal%nDepth == nd) then
                weight = 0.5d0
!             else if (thisOctal%nDepth > nd) then ! fine to coarse
!                weight  = 0.666666666d0
!             else
!                weight  = 0.333333333d0 ! coarse to fine
!             endif

                thisOctal%u_interface(subcell) = &
                     weight*thisOctal%rhov(subcell)/thisOctal%rho(subcell) + (1.d0-weight)*rhou_i_minus_1/rho_i_minus_1

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
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, weight
    integer :: nd
  
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

          if (.not.thisOctal%edgeCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

                rho_i_minus_1 = rho
                rhou_i_minus_1 = rhow

!             if (thisOctal%nDepth == nd) then
                weight = 0.5d0
!             else if (thisOctal%nDepth > nd) then ! fine to coarse
!                weight  = 0.666666666d0
!             else
!                weight  = 0.333333333d0 ! coarse to fine
!             endif


                thisOctal%u_interface(subcell) = &
                     weight*thisOctal%rhow(subcell)/thisOctal%rho(subcell) + (1.d0-weight)*rhou_i_minus_1/rho_i_minus_1

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
    type(VECTOR) :: direction, locator
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi
    integer :: nd
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

          if (.not.thisOctal%edgeCell(subcell)) then
!             thisOctal%x_i(subcell) = (subcellCentre(thisOctal, subcell) .dot. direction) * gridDistanceScale
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisOctal%flux_i_plus_1(subcell) = flux

             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisOctal%flux_i_minus_1(subcell) = flux
          endif
       endif
    enddo
  end subroutine setupFlux


  recursive subroutine setupPressure(thisOctal, grid, direction)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    Type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: direction, locator
    integer :: nd
  
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


          if (.not.thisOctal%edgeCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisOctal%pressure_i_plus_1(subcell) = pressure

             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
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
    type(VECTOR) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    integer :: nd
  

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

          if (.not.thisOctal%edgeCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisOctal%u_i_plus_1(subcell) = rhou/rho
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
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
    type(VECTOR) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    integer :: nd
  

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

          if (.not.thisOctal%edgeCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal =>thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisOctal%u_i_plus_1(subcell) = rhov/rho
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
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
    type(VECTOR) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    integer :: nd
  

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

          if (.not.thisOctal%edgeCell(subcell)) then
             locator = subcellCentre(thisOctal, subcell) + direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisOctal%u_i_plus_1(subcell) = rhow/rho
             
             locator = subcellCentre(thisOctal, subcell) - direction * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisOctal%u_i_minus_1(subcell) = rhow/rho

          endif
       endif
    enddo
  end subroutine setupWpm


   recursive subroutine computePressureU(thisOctal, gamma, direction, iEquationOfState)
     include 'mpif.h'
     integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: gamma, bigGamma,eta
    type(VECTOR) :: direction
    integer :: iEquationOfState
    logical :: useViscosity

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



          thisOctal%pressure_i(subcell) = getPressure(thisOctal, subcell, iEquationOfState, gamma)
          
          bigGamma = 0.d0
          if (.not.thisOctal%edgeCell(subcell)) then
             useViscosity = .false.
             if (thisOctal%u_i_plus_1(subcell) .le. thisOctal%u_i_minus_1(subcell)) useViscosity = .true.
!             if ((thisOctal%u_i_plus_1(subcell) < 0.d0).and.(thisOctal%u_i_plus_1(subcell) .ge. thisOctal%u_i_minus_1(subcell))) &
!                  useViscosity = .true.


             if (useViscosity) then
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
             write(*,*) "cen ",subcellCentre(thisOctal, subcell)
             stop
          endif

       endif
    enddo
  end subroutine computePressureU

   recursive subroutine computePressureV(thisOctal, gamma, direction,  iEquationOfState)
     include 'mpif.h'
     integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: gamma, bigGamma,eta
    type(VECTOR) :: direction
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

          thisOctal%pressure_i(subcell) = getPressure(thisOctal, subcell, iEquationOfState, gamma)

          bigGamma = 0.d0
          if (.not.thisOctal%edgeCell(subcell)) then
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
             write(*,*) "rhoe ", thisOctal%rhoe(subcell)
             write(*,*) "biggamma ", biggamma
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: gamma, bigGamma,eta
    type(VECTOR) :: direction
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

          thisOctal%pressure_i(subcell) = getPressure(thisOctal, subcell, iEquationOfState, gamma)

          bigGamma = 0.d0
          if (.not.thisOctal%edgeCell(subcell)) then
             if (thisOctal%u_i_plus_1(subcell) .le. thisOctal%u_i_minus_1(subcell)) then
                bigGamma = 0.25d0 * eta**2 * (thisOctal%u_i_plus_1(subcell) - thisOctal%u_i_minus_1(subcell))**2 &
                     * thisOctal%rho(subcell)
             else
                bigGamma = 0.d0
             endif
          endif
!          if (myrankGlobal == 1) write(*,*) biggamma/thisOctal%pressure_i(subcell), thisOctal%u_i_plus_1(subcell), &
!               thisOctal%u_i_minus_1(subcell), thisOctal%rho(subcell)
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, rhou, dx
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

             rhou = thisOctal%rhou(subcell)

             dx = (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))
             dx = thisOctal%subcellSize * gridDistanceScale


             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / dx

             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) - dt * & !gravity
                  thisOctal%rho(subcell) *(thisOctal%phi_i_plus_1(subcell) - &
                  thisOctal%phi_i_minus_1(subcell)) / dx

             if (isnan(thisOctal%rhou(subcell))) then
                write(*,*) "rhou ",thisOctal%rhou(subcell)
                write(*,*) "dt ",dt
                write(*,*) "phi ",thisOctal%phi_i_plus_1(subcell), thisOctal%phi_i_minus_1(subcell)
                write(*,*) "x ",thisOctal%x_i_plus_1(subcell), thisOctal%x_i_minus_1(subcell)
             endif

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

!             if (iEquationOfState == 0) then


             if (iEquationOfState /= 1) then 
                thisOctal%rhoE(subcell) = thisOctal%rhoE(subcell) - dt * &
                     (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
                     thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / dx
                     

                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * & !gravity
                  rhou  * (thisOctal%phi_i_plus_1(subcell) - thisOctal%phi_i_minus_1(subcell)) / dx
             endif

!             endif
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, rhou, dx
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

             rhou = thisOctal%rhov(subcell)

             dx = (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))
             dx = thisOctal%subcellSize * gridDistanceScale


             thisOctal%rhov(subcell) = thisOctal%rhov(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / dx

             thisOctal%rhov(subcell) = thisOctal%rhov(subcell) - dt * & !gravity
                  thisOctal%rho(subcell) *(thisOctal%phi_i_plus_1(subcell) - &
                  thisOctal%phi_i_minus_1(subcell)) / dx

             if (isnan(thisOctal%rhov(subcell))) then
                write(*,*) "rhov ",thisOctal%rhov(subcell)
                write(*,*) "dt ",dt
                write(*,*) "phi ",thisOctal%phi_i_plus_1(subcell), thisOctal%phi_i_minus_1(subcell)
                write(*,*) "x ",thisOctal%x_i_plus_1(subcell), thisOctal%x_i_minus_1(subcell)
             endif

             if (thisOctal%rhov(subcell)/thisOctal%rho(subcell) > 1.d10) then
                write(*,*) "u ",thisOctal%rhov(subcell)/thisOctal%rho(subcell)/1.d5
                write(*,*) "rho ", thisOctal%rho(subcell)
                write(*,*) "pressure i + 1 " ,thisOctal%pressure_i_plus_1(subcell)
                write(*,*) "pressure i - 1 " ,thisOctal%pressure_i_minus_1(subcell)
                write(*,*) "x i + 1 ", thisOctal%x_i_plus_1(subcell)
                write(*,*) "x i - 1 ", thisOctal%x_i_minus_1(subcell)
             endif

                
!             Write(*,*) thisOctal%rhou(subcell), thisOctal%pressure_i_plus_1(subcell), &
!                  thisOctal%pressure_i_minus_1(subcell), thisOctal%x_i_plus_1(subcell), &
!                  thisOctal%x_i_minus_1(subcell)

!             if (iEquationOfState == 0) then


             if (iEquationOfState /= 1) then 
                thisOctal%rhoE(subcell) = thisOctal%rhoE(subcell) - dt * &
                     (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
                     thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / dx
                     

                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * & !gravity
                  rhou  * (thisOctal%phi_i_plus_1(subcell) - thisOctal%phi_i_minus_1(subcell)) / dx

             endif
             if (isnan(thisOctal%rhov(subcell))) then
                write(*,*) "bug",thisOctal%rhov(subcell), &
                     thisOctal%pressure_i_plus_1(subcell),thisOctal%pressure_i_minus_1(subcell)
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureForceV


  recursive subroutine pressureForceW(thisOctal, dt, iEquationOfState)
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, rhow, dx
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

             rhow = thisOctal%rhow(subcell)


             dx = (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i_minus_1(subcell))
             dx = thisOctal%subcellSize * gridDistanceScale

             thisOctal%rhoW(subcell) = thisOctal%rhoW(subcell) - dt * &
                  (thisOctal%pressure_i_plus_1(subcell) - thisOctal%pressure_i_minus_1(subcell)) / dx

             thisOctal%rhow(subcell) = thisOctal%rhow(subcell) - dt * & !gravity
                  thisOctal%rho(subcell) *(thisOctal%phi_i_plus_1(subcell) - &
                  thisOctal%phi_i_minus_1(subcell)) / dx


!             write(*,*) thisOctal%rhou(subcell), thisOctal%pressure_i_plus_1(subcell), &
!                  thisOctal%pressure_i_minus_1(subcell), thisOctal%x_i_plus_1(subcell), &
!                  thisOctal%x_i_minus_1(subcell)

             if (iEquationOfState /= 1) then
                thisOctal%rhoE(subcell) = thisOctal%rhoE(subcell) - dt * &
                     (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
                     thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / dx

                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * & !gravity
                  rhow * (thisOctal%phi_i_plus_1(subcell) - thisOctal%phi_i_minus_1(subcell)) / dx
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
  
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%rhoV(subcell) = thisOctal%q_i(subcell)
          endif
        
       endif
    enddo
  end subroutine copyQtoRhoV

  recursive subroutine copyQtoRhoW(thisOctal)
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
  
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%rhoW(subcell) = thisOctal%q_i(subcell)
          endif
        
       endif
    enddo
  end subroutine copyQtoRhoW

  recursive subroutine copyQtoRho(thisOctal, direction)
    type(VECTOR) :: direction
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call copyQtoRho(child, direction)
                exit
             end if
          end do
       else
  
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%rho(subcell) = thisOctal%q_i(subcell)
          endif
          if (thisOCtal%rho(subcell) < 0.d0) then
             write(*,*) "rho warning ", thisOctal%rho(subcell), subcellCentre(thisOctal, subcell)
             write(*,*) "label ",thisOctal%label
             write(*,*) "direction ", direction
             write(*,*) "phi limit ", thisOctal%phiLimit(subcell)
             write(*,*) "flux ", thisOctal%flux_i_plus_1(subcell), thisOctal%flux_i(subcell), &
                  thisOctal%flux_i_minus_1(subcell)
             write(*,*) "u ", thisOctal%u_i_plus_1(subcell)/1.e5, thisOctal%u_interface(subcell)/1.e5, &
                  thisOctal%u_i_minus_1(subcell)/1.e5
             write(*,*) "x ", thisOctal%x_i_plus_1(subcell), thisOctal%x_i(subcell), thisOctal%x_i_minus_1(subcell)
             write(*,*) "dx ", thisOctal%x_i_plus_1(subcell) - thisOCtal%x_i(subcell), thisOctal%x_i(subcell) - &
                  thisOctal%x_i_minus_1(subcell)
             if (abs(thisOctal%x_i_plus_1(subcell) - thisOCtal%x_i(subcell)) > &
                 abs(thisOctal%x_i_minus_1(subcell) - thisOCtal%x_i(subcell))) then
                write(*,*) "fine to coarse"
             endif
             write(*,*) "rho ", thisOctal%rho_i_plus_1(subcell), thisOctal%rho(subcell), thisOctal%rho_i_minus_1(subcell)
             write(*,*) "uinterface ", thisOctal%u_interface(subcell)
             write(*,*) "q ", thisOctal%q_i_plus_1(subcell), thisOctal%q_i(subcell), thisOctal%q_i_minus_1(subcell)
             if (.not.associated(thisOctal%mpiBoundaryStorage)) &
                write(*,*) "Cell is NOT on a boundary"

          endif
       endif
    enddo
  end subroutine copyQtoRho

  recursive subroutine copyQtoRhoE(thisOctal)
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
  
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%rhoE(subcell) = thisOctal%q_i(subcell)
          endif
        
       endif
    enddo
  end subroutine copyQtoRhoE

  recursive subroutine copyQtoRhoU(thisOctal)
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
  
          if (.not.thisOctal%ghostCell(subcell)) then
             thisOctal%rhou(subcell) = thisOctal%q_i(subcell)
             if (isnan(thisOctal%rhou(subcell))) write(*,*) "nan in q to rhou"
          endif
        
       endif
    enddo
  end subroutine copyQtoRhoU


  subroutine advectRho(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup
    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(VECTOR) :: direction

    call copyRhoToQ(grid%octreeRoot)
    call advectQ(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    call copyQtoRho(grid%octreeRoot, direction)

  end subroutine advectRho

  subroutine advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(VECTOR) :: direction

    call copyRhoEToQ(grid%octreeRoot)
    call advectQ(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    call copyQtoRhoE(grid%octreeRoot)

  end subroutine advectRhoE

  subroutine advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(VECTOR) :: direction

    call copyRhoUToQ(grid%octreeRoot)
    call advectQ(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    call copyQtoRhoU(grid%octreeRoot)

  end subroutine advectRhoU

  subroutine advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(VECTOR) :: direction

    call copyRhoVToQ(grid%octreeRoot)
    call advectQ(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    call copyQtoRhoV(grid%octreeRoot)

  end subroutine advectRhoV

  subroutine advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(VECTOR) :: direction

    call copyRhoWToQ(grid%octreeRoot)
    call advectQ(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    call copyQtoRhoW(grid%octreeRoot)

  end subroutine advectRhoW

  subroutine advectQ(grid, direction, dt, nPairs, thread1, thread2, nbound, group, nGroup)
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(VECTOR) :: direction

    call setupX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupQX(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call fluxLimiter(grid%octreeRoot, "superbee")
    call constructFlux(grid%octreeRoot, dt)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupFlux(grid%octreeRoot, grid, direction)
    call updateCellQ(grid%octreeRoot, dt)

  end subroutine advectQ


  subroutine  hydroStep(grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    type(GRIDTYPE) :: grid
    real(double) :: gamma, dt
!    type(VECTOR) :: direction
    integer :: nPairs, thread1(:), thread2(:), nBound(:), group(:), nGroup

!    direction = VECTOR(1.d0, 0.d0, 0.d0)
!
!    call imposeBoundary(grid%octreeRoot)
!    call periodBoundary(grid)
!    call transferTempStorage(grid%octreeRoot)
!
!    direction = VECTOR(1.d0, 0.d0, 0.d0)
!    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!    call setupUi(grid%octreeRoot, grid, direction)
!    call setupUpm(grid%octreeRoot, grid, direction)
!    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
!    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
!    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
!    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!
!    call imposeBoundary(grid%octreeRoot)
!    call periodBoundary(grid)
!    call transferTempStorage(grid%octreeRoot)
!
!    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!    call setupUi(grid%octreeRoot, grid, direction)
!    call setupUpm(grid%octreeRoot, grid, direction)
!    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!    call computePressureU(grid%octreeRoot, gamma, direction, 0)
!    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!    call setupPressure(grid%octreeRoot, grid, direction)
!    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!    call pressureForceU(grid%octreeRoot, dt, 0, 1)
!
!
!    call imposeBoundary(grid%octreeRoot)
!    call periodBoundary(grid)
!    call transferTempStorage(grid%octreeRoot)



  end subroutine hydroStep


  recursive subroutine fullStep3d(iEquationOfState, grid, gamma, timeStep, nPairs, thread1, thread2, &
       nBound, group, nGroup, iDepth, minDepth, maxDepth)
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup
    integer :: iEquationOfState
    integer :: iDepth
    real(double) :: gamma, timeStep, nextTimeStep
    integer :: minDepth, maxDepth, nextDepth


    if (iDepth /= maxDepth) then
       nextTimeStep = timeStep / 2.d0
       nextDepth = iDepth + 1
       call fullStep3d(iEquationOfState, grid, gamma, nextTimeStep, nPairs, thread1, thread2, &
            nBound, group, nGroup, nextDepth, minDepth, maxDepth)
       call fullStep3d(iEquationOfState, grid, gamma, nextTimeStep, nPairs, thread1, thread2, &
            nBound, group, nGroup, nextDepth, minDepth, maxDepth)
    endif

!    call hydroStep3d(iEquationOfState, grid, gamma, timeStep, nPairs, thread1, thread2, nBound, group, &
!         nGroup, iDepth, minDepth, maxDepth)


  end subroutine fullStep3d

    


  subroutine hydroStep3d(iEquationOfState, grid, gamma, dt, nPairs, thread1, thread2, nBound, &
       group, nGroup,doSelfGrav)
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    logical, optional :: doSelfGrav
    logical :: selfGravity
    integer :: group(:), nGroup
    integer :: iEquationOfState
    real(double) :: gamma, dt
    type(VECTOR) :: direction
    integer :: minDepth, maxDepth
    direction = VECTOR(1.d0, 0.d0, 0.d0)

    selfGravity = .true.
    if (PRESENT(doSelfGrav)) selfgravity = doSelfGrav

    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)

    call minMaxDepth(grid%octreeRoot, minDepth, maxDepth)

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoV(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call computePressureU(grid%octreeRoot, gamma, direction, iEquationOfState)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceU(grid%octreeRoot, dt/2.d0, iEquationOfState)

!    call writeVtkFile(grid, "usweep1.vtk", "vtk.txt")


    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupVi(grid%octreeRoot, grid, direction)
    call setupVpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupVi(grid%octreeRoot, grid, direction)
    call setupVpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call computePressureV(grid%octreeRoot, gamma, direction, iEquationOfState)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceV(grid%octreeRoot, dt, iEquationOfState)

!    call writeVtkFile(grid, "vsweep.vtk", "vtk.txt")

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call computePressureW(grid%octreeRoot, gamma, direction, iEquationOfState)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceW(grid%octreeRoot, dt, iEquationOfState)

!    call writeVtkFile(grid, "wsweep.vtk", "vtk.txt")


    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)


    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoV(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call computePressureU(grid%octreeRoot, gamma, direction, iEquationOfState)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceU(grid%octreeRoot, dt/2.d0, iEquationOfState)

!    call writeVtkFile(grid, "usweep2.vtk", "vtk.txt")


    call periodBoundary(grid)
    call imposeBoundary(grid%octreeRoot)
    call transferTempStorage(grid%octreeRoot)


!    call writeVtkFile(grid, "beforeselfgrav.vtk", "vtk.txt")

    if (selfGravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!    call writeVtkFile(grid, "afterselfgrav.vtk", "vtk.txt")


  end subroutine hydroStep3d

  subroutine hydroStep2d(grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup, iEquationofState)
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup
    integer :: iEquationOfState

    real(double) :: gamma, dt
    type(VECTOR) :: direction


    direction = VECTOR(1.d0, 0.d0, 0.d0)

    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)



    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call computePressureU(grid%octreeRoot, gamma, direction, iEquationOfState)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceU(grid%octreeRoot, dt/2.d0, iEquationOfState)

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupWi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call computePressureW(grid%octreeRoot, gamma, direction, iEquationOfState)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceW(grid%octreeRoot, dt, iEquationOfState)


    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call computePressureU(grid%octreeRoot, gamma, direction, iEquationOfState)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call pressureForceU(grid%octreeRoot, dt/2.d0, iEquationOfState)


    call periodBoundary(grid)
    call imposeBoundary(grid%octreeRoot)
    call transferTempStorage(grid%octreeRoot)

 
  end subroutine hydroStep2d


  recursive subroutine computeCourantTime(thisOctal, tc, gamma, iEquationOfState)
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: tc, dx, cs, gamma, speed
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

          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then

             cs = soundSpeed(thisOctal, subcell, gamma, iEquationOfState)
!             if (myrank==1) write(*,*) "cs ", cs/1.d5, "km/s"
             dx = thisOctal%subcellSize * gridDistanceScale
!             speed = max(thisOctal%rhou(subcell)**2, thisOctal%rhov(subcell)**2, thisOctal%rhow(subcell)**2)
             speed = thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2
             speed = sqrt(speed)/thisOctal%rho(subcell)
                tc = min(tc, dx / (cs + speed) )


          endif
 
       endif
    enddo
  end subroutine computeCourantTime

  function soundSpeed(thisOctal, subcell, gamma, iEquationOfState) result (cs)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: cs, gamma
    real(double) :: u2, eKinetic, eTot, eThermal
    logical, save :: firstTime = .true.
    integer :: iEquationOfState
    real(double), parameter :: gamma2 = 1.4d0, rhoCrit = 1.d-14
    select case(iEquationOfState)
       case(0) ! adiabatic 
          if (thisOctal%threed) then
             u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
          else
             u2 = (thisOctal%rhou(subcell)**2 +  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
          endif
          eKinetic = u2 / 2.d0
          eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
          eThermal = eTot - eKinetic
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
       case(1) ! isothermal
          cs = sqrt(getPressure(thisOctal, subcell, iEquationOfState, gamma)/thisOctal%rho(subcell))
       case(2) ! barotropic
          if (thisOctal%rho(subcell) < rhoCrit) then
             cs = sqrt( getPressure(thisOctal, subcell, iEquationOfState, gamma)/thisOctal%rho(subcell))
          else
             cs = sqrt(gamma2 *  getPressure(thisOctal, subcell, iEquationOfState, gamma)/thisOctal%rho(subcell))
          endif
       case(3) ! polytropic
          cs = sqrt(gamma*thisOctal%pressure_i(subcell)/thisOctal%rho(subcell))
       case DEFAULT
          write(*,*) "Unknown equation of state passed to sound speed ", iEquationofState
          
    end select


  end function soundSpeed

  subroutine doHydrodynamics1d(grid)
    include 'mpif.h'
    type(gridtype) :: grid
    real(double) :: dt,  gamma, mu
    real(double) :: currentTime
    integer :: i
    type(VECTOR) :: direction
    logical :: gridConverged
    integer :: iEquationOfState = 0
    logical :: globalConverged(10)
    logical :: tConverged(10)
    integer :: nHydroThreads
    real(double) :: tc(10)
    integer :: it
    integer :: myRank, ierr
    integer :: nPairs, thread1(100), thread2(100), group(100), nBound(100), ngroup
    integer :: iUnrefine, nUnrefine
    real(double) :: nextDumpTime, tdump, temptc(10)
    character(len=80) :: plotfile

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    mu = 2.d0
    nHydroThreads = nThreadsGlobal - 1


    direction = VECTOR(1.d0, 0.d0, 0.d0)
    gamma = 5.d0 / 3.d0

    mu = 2.d0


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    if (myRank == 1) write(*,*) "CFL set to ", cflNumber

    call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)

    do i = 1, nPairs
       if (myrankglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
    enddo


    call writeInfo("Calling exchange across boundary", TRIVIAL)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call writeInfo("Done", TRIVIAL)

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 0.d0, 1.d0)
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
       call evenUpGridMPI(grid,.false.,.true.)
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
       if (ALL(tConverged(1:nHydrothreads))) exit
    end do


    call writeInfo("Evening up grid", TRIVIAL)    
    call evenUpGridMPI(grid, .false.,.true.)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)



    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)
    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)

    currentTime = 0.d0
    it = 0
    nextDumpTime = 0.d0
    tDump = 0.005d0

    iUnrefine = 0
    do while(.true.)
       tc = 0.d0
       tc(myrank) = 1.d30
       call computeCourantTime(grid%octreeRoot, tc(myRank), gamma, iEquationOfState)
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
!       write(*,*) "tc", tc(1:8)
!       write(*,*) "temp tc",temptc(1:8)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)

       if (myrank == 1) write(*,*) "courantTime", dt
       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call writeInfo("calling hydro step",TRIVIAL)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


       call hydroStep(grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup)

!       direction = VECTOR(1.d0, 0.d0, 0.d0)
!       call roeSolver(grid, direction, gamma, dt, nPairs, thread1, thread2, nBound, group, ngroup)

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
          nUnrefine = 0 
          call unrefineCells(grid%octreeRoot, grid, gamma, iEquationOfState, nUnrefine)
          write(*,*) "Unrefined ", nUnrefine, " cells"
          call tune(6, "Unrefine grid")
          iUnrefine = 0
       endif

       call evenUpGridMPI(grid, .true., .true.)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


!     call plot_AMR_values(grid, "rho", "x-y", 0., &
 !           "/xs",.false., .true., fixvalmin=0.d0, fixvalmax=1.d0, quiet=.true.)



       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt
       if (currentTime .gt. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
          grid%iDump = it
          grid%currentTime = currentTime

       endif


    enddo

  end subroutine doHydrodynamics1d


  subroutine doHydrodynamics3d(grid)
    use input_variables, only : zoomFactor
    include 'mpif.h'
    type(gridtype) :: grid
    real(double) :: dt, tc(64), temptc(64), gamma, mu
    real(double) :: currentTime, smallTime
    integer :: i, it, iUnrefine
    integer :: myRank, ierr
    character(len=20) :: plotfile
    real(double) :: tDump, nextDumpTime, tff !, ang
    type(VECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: nSource = 0
    integer :: nUnrefine
    integer :: thread1(200), thread2(200), nBound(200), nPairs
    logical :: globalConverged(64), tConverged(64)
    integer :: nHydroThreads
    integer :: nGroup, group(200)
    integer :: iEquationOfState = 2
    logical :: doRefine

    dorefine = .true.

    nHydroThreads = nThreadsGlobal - 1

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0

    mu = 2.d0

    viewVec = VECTOR(-1.d0,0.d0,0.d0)
!    viewVec = rotateZ(viewVec, 20.d0*degtorad)
    viewVec = rotateY(viewVec, 25.d0*degtorad)
    
    it = grid%iDump
    currentTime = grid%currentTime
    nextDumpTime = grid%currentTime
!    tDump = 0.01d0
    tff = 1.d0 / sqrt(bigG * (1d0*mSol/((4.d0/3.d0)*pi*7.d15**3)))
    tDump = 1.d-2 * tff
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    if (myRank == 1) write(*,*) "CFL set to ", cflNumber



    call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)

    do i = 1, nPairs
       if (myrankglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
    enddo


    call writeInfo("Calling exchange across boundary", TRIVIAL)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call writeInfo("Done", TRIVIAL)

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)

!    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)


    if (myrank == 1) call tune(6, "Initial refine")
    
!    call writeVtkFile(grid, "beforerefine.vtk")

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
       call evenUpGridMPI(grid, .false., dorefine)
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
    call evenUpGridMPI(grid,.false., dorefine)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)



    if (myrank == 1) call tune(6, "Initial refine")

!    call writeVtkFile(grid, "afterrefine.vtk")

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
!    call calculateEnergy(grid%octreeRoot, gamma, mu)
    call calculateRhoE(grid%octreeRoot, direction)
    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)




    if (myrankglobal==1) write(*,*) "Setting tdump to: ", tdump
    nextDumpTime = tdump + currentTime
    iUnrefine = 0
!    call writeInfo("Plotting grid", TRIVIAL)    

    iUnrefine = 0

    call writeInfo("Plotting col density", TRIVIAL)    

    if (myrank == 1) call tune(6, "Self-Gravity")
    if (myrank == 1) write(*,*) "Doing multigrid self gravity"
    call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.) 
    if (myrank == 1) write(*,*) "Done"
    if (myrank == 1) call tune(6, "Self-Gravity")

    do while(.true.)
       tc = 0.d0
       tc(myrank) = 1.d30
       call computeCourantTime(grid%octreeRoot, tc(myRank), gamma, iEquationOfState)
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
!       write(*,*) "tc", tc(1:8)
       if (myrank==1)write(*,*) "temp tc",temptc(1:nHydroThreads)
       dt = MINVAL(temptc(1:nHydroThreads)) * dble(cflNumber)

       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

!       smallTime = 0.3d0 * 1.d4*gridDistanceScale/2.e5
!       write(*,*) "myrank ",myrank, "smallest cell ",grid%halfSmallestSubcell
       if (myrank == 1) write(*,*) "courantTime", dt,cflnumber
!       if (myrank == 1) write(*,*) "smallTime", smallTime
       if (myrank == 1) call tune(6,"Hydrodynamics step")
!       if (dt > smallTime) dt = smallTime
       call writeInfo("calling hydro step",TRIVIAL)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


       call hydroStep3d(iEquationOfState, grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup, doSelfGrav=.true.)


       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)



       if (myrank == 1) call tune(6, "Loop refine")
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
          if (ALL(tConverged(1:nHydroThreads))) exit
       end do
       
       iUnrefine = iUnrefine + 1
       if (iUnrefine == 5) then
          if (myrank == 1)call tune(6, "Unrefine grid")
          nUnrefine = 0
          call unrefineCells(grid%octreeRoot, grid, gamma, iEquationOfState, nUnrefine)
          write(*,*) "Unrefined ", nUnrefine, " cells"
          if (myrank == 1)call tune(6, "Unrefine grid")
          iUnrefine = 0
       endif
!          
       call evenUpGridMPI(grid, .true., dorefine)

       call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)


       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       
       if (myrank == 1) call tune(6, "Loop refine")
!

       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt,nextDumpTime
!       stop

       if (currentTime .ge. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1


          write(plotfile,'(a,i4.4,a)') "rho",it,".vtk"
          call writeVtkFile(grid, plotfile, "vtk.txt")


          write(plotfile,'(a,i4.4,a)') "dump",it,".grid"
          grid%iDump = it
          grid%currentTime = currentTime
          call writeAmrGridMpiAll(plotfile,.false.,grid)

       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)


    enddo
  end subroutine doHydrodynamics3d

  subroutine doHydrodynamics2d(grid)
    include 'mpif.h'
    type(gridtype) :: grid
    real(double) :: dt, tc(64), temptc(64), gamma, mu
    real(double) :: currentTime
    integer :: i, it, iUnrefine
    integer :: myRank, ierr
    character(len=20) :: plotfile, titleString
    real(double) :: tDump, nextDumpTime, tff !, ang
    type(VECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: nSource = 0
    integer :: thread1(100), thread2(100), nBound(100), nPairs
    integer :: nGroup, group(100)

    logical :: globalConverged(64), tConverged(64)
    integer :: nHydroThreads 
    integer :: iEquationOfState = 1
    logical :: dorefine

    dorefine = .true.
    

    nHydroThreads = nThreadsGlobal - 1


    direction = VECTOR(1.d0, 0.d0, 0.d0)
    gamma = 5.d0 / 3.d0

    mu = 2.d0

    viewVec = VECTOR(-1.d0,0.d0,0.d0)
!    viewVec = rotateZ(viewVec, 20.d0*degtorad)
    viewVec = rotateY(viewVec, 25.d0*degtorad)
    

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    if (myRank == 1) write(*,*) "CFL set to ", cflNumber


!    call writeInfo("Plotting grid", TRIVIAL)    
!    write(plotfile,'(a,i2.2,a)') "test",myRank,".png/png"
!    call plot_AMR_values(grid, "rho", "x-z", 0., &
!           plotfile,.false., .true.)



    call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)

    do i = 1, nPairs
       if (myrankglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
    enddo




    call writeInfo("Calling exchange across boundary", TRIVIAL)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    call writeInfo("Done", TRIVIAL)

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)

!    call calculateEnergy(grid%octreeRoot, gamma, mu)
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
      call evenUpGridMPI(grid,.false., dorefine)
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
    call evenUpGridMPI(grid, .false.,dorefine)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)




    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)
!    call calculateEnergy(grid%octreeRoot, gamma, mu)


    call calculateRhoE(grid%octreeRoot, direction)
    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call calculateRhoU(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call calculateRhoV(grid%octreeRoot, direction)
    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call calculateRhoW(grid%octreeRoot, direction)

    currentTime = 0.d0
    it = 0
    nextDumpTime = 0.d0


    tc = 0.d0
    tc(myrank) = 1.d30
    call computeCourantTime(grid%octreeRoot, tc(myRank), gamma, iEquationOfState)
    call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
    tc = tempTc
    dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)

    tff = 1.d0 / sqrt(bigG * (0.1d0*mSol/((4.d0/3.d0)*pi*7.d15**3)))
    tDump = 1.d-2 * tff


    tdump = 1.d7
    write(*,*) "Setting tdump to: ", tdump

    iUnrefine = 0
!    call writeInfo("Plotting grid", TRIVIAL)    

    iUnrefine = 0



    do while(.true.)
       tc = 0.d0
       tc(myrank) = 1.d30
       call computeCourantTime(grid%octreeRoot, tc(myRank), gamma, iEquationOfState)
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
!       write(*,*) "tc", tc(1:8)
!       write(*,*) "temp tc",temptc(1:8)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)

       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

       if (myrank == 1) write(*,*) "courantTime", dt
       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call writeInfo("calling hydro step",TRIVIAL)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


       call hydroStep2d(grid, gamma, dt, nPairs, thread1, thread2, nBound, group, nGroup, iEquationOfState)



       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!       call writeInfo("Refining grid", TRIVIAL)
!       do
!          gridConverged = .true.
!          call refineGridGeneric2(grid%octreeRoot, grid,  gamma, gridconverged, iEquationOfState, inherit=.true.)
!          if (gridConverged) exit
!       end do
!       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!       
!       call writeInfo("Refining grid part 2", TRIVIAL)    
!       do
!          globalConverged(myRank) = .true.
!          call writeInfo("Refining grid", TRIVIAL)    
!          call refineGridGeneric2(grid%octreeRoot, grid,  gamma, globalConverged(myRank), iEquationOfState, inherit=.true.)
!          call writeInfo("Exchanging boundaries", TRIVIAL)    
!          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!          call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreads, MPI_LOGICAL, MPI_LOR,amrCOMMUNICATOR, ierr)
!          if (ALL(tConverged(1:nHydroThreads))) exit
!       end do
!       
!       iUnrefine = iUnrefine + 1
!       if (iUnrefine == 5) then
!          call tune(6, "Unrefine grid")
!          call unrefineCells(grid%octreeRoot, grid, gamma, iEquationOfState)
!          call tune(6, "Unrefine grid")
!          iUnrefine = 0
!       endif

       call evenUpGridMPI(grid, .true., dorefine)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


!     call plot_AMR_values(grid, "rho", "x-y", 0., &
 !           "/xs",.false., .true., fixvalmin=0.d0, fixvalmax=1.d0, quiet=.true.)



       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt
       if (currentTime .ge. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
          write(plotfile,'(a,i4.4,a)') "dump",it,".grid"
          grid%iDump = it
          grid%currentTime = currentTime
          call writeAmrGridMpiAll(plotfile,.false.,grid)

       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)


    enddo
  end subroutine doHydrodynamics2d

  recursive subroutine calculateEnergy(thisOctal, gamma, mu)
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
 
          thisOctal%energy(subcell) = thisOctal%pressure_i(subcell) / ((gamma-1.d0)*thisOctal%rho(subcell))

          thisOctal%energy(subcell) = thisOctal%energy(subcell) + 0.5d0 * (thisOctal%rhou(subcell)**2 + &
               thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
        
       endif
    enddo
  end subroutine calculateEnergy


  recursive subroutine zeroRefinedLastTime(thisOctal)
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
             thisOctal%energy(subcell) = thisOctal%tempStorage(subcell,6)
             thisOctal%pressure_i(subcell) = thisOctal%tempStorage(subcell,7)
          endif
       endif
    enddo
    if (associated(thisOCtal%tempStorage)) Deallocate(thisOctal%tempStorage)
  end subroutine transferTempStorage

  recursive subroutine calculateRhoV(thisOctal, direction)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: direction
  
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: direction
  
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: direction
  
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: direction
  
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

  recursive subroutine getArray(thisOctal, x, rho, rhou, v, n)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: x(:), rho(:), v(:), rhou(:)
    integer :: n
    type(VECTOR) :: rVec
  
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
    type(VECTOR) :: locator, dir
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

             if (.not.associated(thisOctal%tempStorage)) allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:7))

             select case(thisOctal%boundaryCondition(subcell))
                case(1)
                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocal(locator, bOctal, bSubcell)

                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)
                   Thisoctal%tempstorage(subcell,1) = bOctal%rho(bSubcell)
                   thisOctal%tempStorage(subcell,2) = bOctal%rhoE(bSubcell)

                   thisOctal%tempStorage(subcell,6) = bOctal%energy(bSubcell)
                   thisOctal%tempStorage(subcell,7) = bOctal%pressure_i(bSubcell)

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



                case(3) ! shock
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

                case(4) ! free outflow, no inflow

                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocal(locator, bOctal, bSubcell)

                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)
                   Thisoctal%tempstorage(subcell,1) = bOctal%rho(bSubcell)
                   thisOctal%tempStorage(subcell,2) = bOctal%rhoE(bSubcell)

                   thisOctal%tempStorage(subcell,6) = bOctal%energy(bSubcell)
                   thisOctal%tempStorage(subcell,7) = bOctal%pressure_i(bSubcell)

! NB confusion regarding 2d being x,z rather than x,y

                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (dir%x > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) =-abs(bOctal%rhou(bSubcell))
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (dir%x < -0.9d0) then
                         thisOctal%tempStorage(subcell,3) = abs(bOctal%rhou(bSubcell))
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (dir%z > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -abs(bOctal%rhow(bSubcell))
                      endif
                      if (dir%z < -0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = abs(bOctal%rhow(bSubcell))
                      endif
                   else if (thisOctal%threed) then
                      if (dir%x > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) =-abs(bOctal%rhou(bSubcell))
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (dir%x < -0.9d0) then
                         thisOctal%tempStorage(subcell,3) = abs(bOctal%rhou(bSubcell))
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (dir%y > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = -abs(bOctal%rhov(bSubcell))
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (dir%y < -0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = abs(bOctal%rhov(bSubcell))
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (dir%z > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -abs(bOctal%rhow(bSubcell))
                      endif
                      if (dir%z < -0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = abs(bOctal%rhow(bSubcell))
                      endif



                   endif

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
    type(VECTOR) :: locator, rVec
    integer :: nProbes, iProbe
    type(VECTOR) :: probe(6), direction
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
             probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          endif
          if (thisOctal%twod) then
             nProbes = 4
             probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
             probe(4) = VECTOR(0.d0, 0.d0, -1.d0)
          endif
          if (thisOctal%threeD) then
             nProbes = 6
             probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             probe(3) = VECTOR(0.d0, 1.d0, 0.d0)
             probe(4) = VECTOR(0.d0, -1.d0, 0.d0)
             probe(5) = VECTOR(0.d0, 0.d0, 1.d0)
             probe(6) = VECTOR(0.d0, 0.d0, -1.d0)
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

          corner = .false.

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

                   case(2) ! periodic
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


                   case(4)
                      dx = thisOctal%subcellSize
                      thisOctal%ghostCell(subcell) = .true.

                      call locatorToNeighbour(grid, thisOctal, subcell, (-1.d0)*probe(iProbe), 2, locator)
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
             case(1, 4)
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: locator, rVec
    integer :: nProbes, iProbe
    type(VECTOR) :: probe(6)
    integer :: nProbeOutside
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
             probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          endif
          if (thisOctal%twod) then
             nProbes = 4
             probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
             probe(4) = VECTOR(0.d0, 0.d0, -1.d0)
          endif
          if (thisOctal%threeD) then
             nProbes = 6
             probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             probe(3) = VECTOR(0.d0, 1.d0, 0.d0)
             probe(4) = VECTOR(0.d0, -1.d0, 0.d0)
             probe(5) = VECTOR(0.d0, 0.d0, 1.d0)
             probe(6) = VECTOR(0.d0, 0.d0, -1.d0)
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
    integer :: i, subcell
    logical :: inherit
    logical :: split
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
    type(VECTOR) :: locator, probe(6)
    integer :: i, nd
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi
    real(double) :: grad, maxGradient
    
    split = .false.
    if (thisOctal%oned) then
       nProbes = 2
       probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
       probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
    endif
    if (thisOctal%twod) then
       nProbes = 4
       probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
       probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
       probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
       probe(4) = VECTOR(0.d0, 0.d0, -1.d0)
    endif
    if (thisOctal%threeD) then
       nProbes = 6
       probe(1) = VECTOR( 0.d0, 0.d0, +1.d0)
       probe(2) = VECTOR( 0.d0, 0.d0, -1.d0)
       probe(3) = VECTOR(-1.d0, 0.d0,  0.d0)
       probe(4) = VECTOR(+1.d0, 0.d0,  0.d0)
       probe(5) = VECTOR( 0.d0, 1.d0,  0.d0)
       probe(6) = VECTOR( 0.d0,-1.d0,  0.d0)
    endif

    maxGradient = 1.d-30
    do i = 1, nProbes
       locator = subcellCentre(thisOctal, subcell) + &
            (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * probe(i)
       if (inOctal(grid%octreeRoot, locator)) then
          neighbourOctal => thisOctal
          call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
          call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(i), q, rho, rhoe, &
               rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
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



  recursive subroutine refineGridGeneric2(thisOctal, grid,  gamma, converged, iEquationOfState, inherit)
    use input_variables, only : maxDepthAMR
    include 'mpif.h'
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
    type(VECTOR) :: dirVec(6), centre, locator
    logical :: split
    integer :: neighbourSubcell, nDir
    real(double) :: r , grad, maxGradient
    logical, optional :: inherit
    real(double), parameter :: limit = 0.1d0
    real(double) :: gamma
    real(double) :: cs, rhocs
    integer :: myRank, ierr
    integer :: iEquationOfState
    logical :: refineOnMass, refineOnIonization
    converged = .true.
    converged_tmp=.true.


    refineOnMass = .true.
    refineOnIonization = .false.

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


          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell)
          if (thisOctal%threed) then
             nDir = 6
             dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
             dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
             dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          endif


!          do i = 1, nDir
!             maxGradient = 1.d-30
!             locator = subcellCentre(thisOctal, subcell) + &
!                  (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
!             if (inOctal(grid%octreeRoot, locator)) then
!                neighbourOctal => thisOctal
!                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
!                
!                grad = abs((thisOctal%rho(subcell)-neighbourOctal%rho(neighbourSubcell)) / &
!                     thisOctal%rho(subcell))
!                maxGradient = max(grad, maxGradient)
!
!                !                if (grad > limit) write(*,*) "rho",grad
!
!                !                if (thisOctal%rhoe(subcell) /= 0.d0) then
!                !                   grad = abs((thisOctal%rhoe(subcell)-neighbourOctal%rhoe(neighbourSubcell)) / &
!                !                        thisOctal%rhoe(subcell))
!                !                   maxGradient = max(grad, maxGradient)
!                !                endif
!                !                if (grad > limit) write(*,*) "rhoe",grad, thisOctal%rhoe(subcell), neighbourOctal%rhoe(neighboursubcell),&
!                !                     thisOctal%mpiThread(subcell), neighbourOctal%mpithread(neighboursubcell)
!                !                if (.not.thisOctal%ghostCell(subcell)) then
!                !                   cs = soundSpeed(thisOctal, subcell, gamma, iEquationOfState)
!                !                   rhocs = thisOctal%rho(subcell) * cs
!                !
!                !                   if (rhocs /= 0.d0) then
!                !                      grad = abs((thisOctal%rhou(subcell)-neighbourOctal%rhou(neighbourSubcell)) / &
!                !                           rhocs)
!                !                      maxGradient = max(grad, maxGradient)
!                !!                      if (grad > limit) write(*,*) "rhou",grad
!                !                      
!                !                      grad = abs((thisOctal%rhov(subcell)-neighbourOctal%rhov(neighbourSubcell)) / &
!                !                        rhocs)
!                !                      maxGradient = max(grad, maxGradient)
!                !!                      if (grad > limit) write(*,*) "rhov",grad
!                !                      
!                !                      grad = abs((thisOctal%rhow(subcell)-neighbourOctal%rhow(neighbourSubcell)) / &
!                !                           rhocs)
!                maxGradient = max(grad, maxGradient)
!                !                      if (grad > limit) write(*,*) "rhow",grad
!             endif
!             if (maxGradient > limit) then
!                split = .true.
!             else
!                split = .false.
!             endif
!
!
!             if (split) then
!                if ((neighbourOctal%nDepth >= thisOctal%nDepth).and.(thisOctal%nDepth < maxDepthAMR)) then
!                   call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
!                        inherit=inherit, interp=.false.)
!                   converged = .false.
!                   exit
!                endif
!
!                !                      if ((thisOctal%nDepth >= neighbourOctal%nDepth).and.(neighbourOctal%nDepth < maxDepth)) then
!                !                         call addNewChild(neighbourOctal,neighboursubcell,grid,adjustGridInfo=.TRUE., &
!                !                              inherit=inherit, interp=.false.)
!                !                         converged = .false.
!                !                         exit
!                !                      endif
!
!             endif
!          enddo
!
!
    if (refineOnMass) then
       if (((thisOctal%rho(subcell)*cellVolume(thisOctal, subcell)*1.d30) > 1.d-5*mSol) &
            .and.(thisOctal%nDepth < maxDepthAMR))  then
          call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
               inherit=inherit, interp=.false.)
          converged = .false.
       endif
    endif


    if (refineOnIonization) then

       do i = 1, nDir
          locator = subcellCentre(thisOctal, subcell) + &
               (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
          if (inOctal(grid%octreeRoot, locator)) then
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             if ((thisOctal%ionFrac(subcell,1) > 0.9d0).and.(neighbourOctal%ionFrac(neighbourSubcell,1) < 0.1d0) .and. &
                  (thisOctal%nDepth < maxDepthAMR) ) then
                call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                     inherit=inherit, interp=.false.)
                converged = .false.
                exit
             endif

             if ((thisOctal%ionFrac(subcell,1) < 0.1d0).and.(neighbourOctal%ionFrac(neighbourSubcell,1) > 0.9d0) .and. &
                  (thisOctal%nDepth < maxDepthAMR) ) then
                call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                     inherit=inherit, interp=.false.)
                converged = .false.
                exit
             endif

          endif
       enddo
    endif



    if (.not.converged) exit
 endif
end do

end subroutine refineGridGeneric2





  recursive subroutine refineEdges(thisOctal, grid,  converged, inherit)

    use input_variables, only : maxDepthAMR
    include 'mpif.h'
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i
    logical :: converged, converged_tmp
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
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
    type(VECTOR) :: direction, locator, rVec
    tempOctal => thisOctal
    tempSubcell = subcell
    do i = 1, nCells
       rVec = subcellCentre(tempOctal, tempSubcell) + &
            (tempOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell) * direction 
       call findSubcellLocal(rVec, tempOctal, tempSubcell)
    enddo
    locator = rVec
  end subroutine locatorToNeighbour
       
  recursive subroutine unrefineCells(thisOctal, grid, gamma, iEquationOfState, nUnrefine)
    use input_variables, only : minDepthAMR
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: gamma
    integer :: subcell, i
    logical :: unrefine
    integer :: nc, nUnrefine
    real(double) :: rhow(8), rhov(8), rhou(8), rho(8), rhoe(8) !, fac, limit
    real(double) :: cs(8), mass
!    real(double) :: rhocs, rhomean, rhoemean
    logical :: refinedLastTime, ghostCell
    integer :: iEquationOfState
!    limit  = 0.1d0

    unrefine = .true.
    refinedLastTime = .false.
    ghostCell = .false.
    nc = 0
    mass = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unrefineCells(child, grid,  gamma, iEquationOfState, nUnrefine)
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

          mass = mass +  thisOctal%rho(subcell)*cellVolume(thisOctal, subcell) * 1.d30

          cs(nc) = soundSpeed(thisOctal, subcell, gamma, iEquationOfState)
          if (thisOctal%ghostCell(subcell)) ghostCell=.true.
!          if (thisOctal%refinedLastTime(subcell)) refinedLastTime = .true.
       endif
    enddo


    unrefine = .false.

!    if ((nc > 1).and.(.not.ghostCell)) then
!
!       rhomean = SUM(rho(1:nc))/dble(nc)
!       rhoemean = SUM(rhoe(1:nc))/dble(nc)
!       rhocs = rhomean*(SUM(cs(1:nc))/dble(nc))
!       unrefine = .true.
!       
!       fac = sigma(rho, nc) / rhoMean
!       if (fac > limit) then
!          unrefine = .false.
!!          write(*,*) "rho",fac
!       endif
!      
!       fac = sigma(rhoe, nc) / rhoeMean
!       if (fac > limit) then
!          unrefine = .false.
!!          write(*,*) "rhoe",fac
!       endif
!       
!       fac = sigma(rhou, nc) / rhocs
!       if (fac > limit) then
!          unrefine = .false.
!!          write(*,*) "rhou",fac
!       endif
!       
!       fac = sigma(rhov, nc) / rhocs
!       if (fac > limit) then
!          unrefine = .false.
!!          write(*,*) "rhov",fac
!       endif
!
!       fac = sigma(rhow, nc) / rhocs
!       if (fac > limit) then
!          unrefine = .false.
!!          write(*,*) "rhov",fac
!       endif
!       

    if (mass  < 2.d-6*mSol) unrefine = .true.
       
    if (thisOctal%nDepth <= minDepthAMR) unrefine = .false.

    if ((thisOctal%nChildren == 0).and.unrefine) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
       nunrefine = nunrefine + 1
    endif
    
  end subroutine unrefineCells

  recursive subroutine unrefineCellsPhotoion(thisOctal, grid, gamma, iEquationOfState)
    use input_variables, only : minDepthAMR
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: gamma
    integer :: subcell, i
    logical :: unrefine
    integer :: nc
    real(double) :: ionfrac(8), limit
    real(double) :: mass
    logical :: refinedLastTime, ghostCell
    integer :: iEquationOfState
    limit  = 0.1d0

    unrefine = .true.
    refinedLastTime = .false.
    ghostCell = .false.
    nc = 0
    mass = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unrefineCellsPhotoIon(child, grid,  gamma, iEquationOfState)
                exit
             end if
          end do
       else
          nc = nc + 1
          ionfrac(nc) = thisOctal%ionFrac(subcell,1)
       endif
    enddo


    unrefine = .false.

    if (ALL(ionFrac(1:nc) < 0.01d0).or.ALL(ionFrac(1:nc)> 0.99d0)) unrefine=.true.

    if (thisOctal%nDepth <= minDepthAMR) unrefine = .false.

    if ((thisOctal%nChildren == 0).and.unrefine) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
    endif
    
  end subroutine unrefineCellsPhotoIon
 

  function sigma(x, n)
    real(double) :: x(:), sigma !, mean
    integer :: n !, i

!    mean = SUM(x(1:n))/dble(n)
!    do i = 1, n
!       sigma = sigma + (x(i) - mean)**2
!    enddo
!    sigma = sqrt(sigma / dble(n))
    sigma = maxval(abs(x(1:n)))-minval(abs(x(1:n)))
  end function sigma

  subroutine evenUpGridMPI(grid, inheritFlag, evenAcrossThreads)

    type(GRIDTYPE) :: grid
    logical :: gridConverged, inheritFlag
    integer :: myRank, nThreads
    include 'mpif.h'  
    integer :: ierr
    logical :: globalConverged, localChanged(64)
    type(VECTOR) :: locs(200000), eLocs(200000)
    integer :: nLocs(64), tempNlocs(20000)
    integer :: thread(100000), nLocsGlobal,i, depth(200000)

    real(double) :: temp(20000,4),tempsent(4)
    integer :: nTemp(1), nSent(1), eDepth(200000)
    integer :: iThread, nExternalLocs
    integer :: iter
    integer, parameter :: tag = 1
    logical :: globalChanged(64)
    integer :: status(MPI_STATUS_SIZE)
    logical :: evenAcrossThreads
!    character(len=20) :: plotfile
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
!          call refineEdges(grid%octreeRoot, grid,  gridconverged, inherit=inheritFlag)
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

       if (evenAcrossThreads) then
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
                      call mpi_send(temp(i,1:4), 4, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
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
       else
          globalConverged = .true.
       endif
    end do
  end subroutine evenUpGridMPI

  recursive subroutine evenUpGrid(thisOctal, grid,  converged, inherit)

    include 'mpif.h'
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
    type(VECTOR) :: dirVec(6), centre, octVec
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
             dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
             dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
             dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)  
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
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

  subroutine splitAtLocator(grid, locator, depth,  localchanged, inherit)
    use input_variables, only :  maxDepthAMR
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    type(VECTOR) :: locator
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
    type(VECTOR) :: dirVec(6), centre, octVec, loc(:)
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
             dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
             dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
             dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          endif

             do j = 1, nDir
                octVec = centre + r * dirvec(j)
                if (inOctal(grid%octreeRoot, octVec)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)

                   if (.not.octalOnThread(neighbourOctal, neighbourSubcell, myRank)) then
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
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
    type(VECTOR) :: dirVec(6), centre, octVec
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
             dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
             dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
             dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
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


  subroutine writeAMRgridMpiALL(filename,fileFormatted,grid)
    include 'mpif.h'
    character(len=*) :: filename
    type(GRIDTYPE) :: grid
    logical :: fileformatted
    integer :: iThread
    integer :: ierr

    do iThread = 1, nThreadsGlobal-1
       if (myRankGlobal == iThread) then
          call writeAMRgridMpi(filename,fileFormatted,grid,mpiFlag=.true.)
       endif
       call MPI_BARRIER(amrCommunicator, ierr)
    enddo
  end subroutine writeAMRgridMpiALL

  subroutine readAMRgridMpiALL(filename,fileFormatted,grid)
    include 'mpif.h'
    character(len=*) :: filename
    type(GRIDTYPE) :: grid
    logical :: fileformatted
    integer :: iThread
    integer :: ierr
    integer :: nOctals, nVoxels

    do iThread = 0, nThreadsGlobal-1
       if (myRankGlobal == iThread) then
          call readAMRgridMpi(filename,fileFormatted,grid,mpiFlag=.true.)
          call updateMaxDepth(grid, searchLimit = 100000)
          call setSmallestSubcell(grid)
          call countVoxels(grid%octreeRoot,nOctals,nVoxels)  
          grid%nOctals = nOctals
       endif
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    enddo
  end subroutine readAMRgridMpiALL



  recursive subroutine gSweep(thisOctal, grid, deltaT, fracChange)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6)
    integer :: n, ndir, nd
    real(double) :: x1, x2, g(6)
    real(double) :: deltaT, fracChange, gGrav, newPhi, frac !, d2phidx2(3), sumd2phidx2

    gGrav = bigG

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call gsweep(child, grid, deltaT, fracChange)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%twoD) then
             nDir = 4
             dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
             dir(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
             dir(4) = VECTOR(0.d0, 0.d0,-1.d0)
          else if (thisOctal%threed) then
             nDir = 6
             dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
             dir(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
             dir(4) = VECTOR(0.d0, 0.d0,-1.d0)
             dir(5) = VECTOR(0.d0, 1.d0, 0.d0)
             dir(6) = VECTOR(0.d0,-1.d0, 0.d0)
          endif



          if (.not.thisOctal%edgeCell(subcell)) then !xxx changed from ghostcell


             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + dir(n) * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, dir(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
                
                x1 = subcellCentre(thisOctal, subcell) .dot. dir(n)
                x2 = subcellCentre(neighbourOctal, neighboursubcell) .dot. dir(n)

!                g(n) =   (phi - thisOctal%phi_i(subcell))/(x2 - x1)
!                g(n) =   (phi - thisOctal%phi_i(subcell))/thisOctal%subcellSize !!!!!!!!!!!!!!!!!!!!!!!!!!

!                if (slow) then
                   g(n) = phi
!                endif

             enddo

             if (thisOctal%twoD) then
                newphi = 0.25d0*(SUM(g(1:4))) - fourpi*gGrav*thisOctal%rho(subcell)*deltaT
             else
                newphi = 0.16666666666667d0*(SUM(g(1:6))) - fourpi*gGrav*thisOctal%rho(subcell)*deltaT
             endif!
!             if (thisOctal%twoD) then
!                d2phidx2(1) = (g(1) - g(2)) / thisOctal%subcellSize
!                d2phidx2(2) = (g(3) - g(4)) / thisOctal%subcellSize
!                sumd2phidx2 = SUM(d2phidx2(1:2))
!             else
!                d2phidx2(1) = (g(1) - g(2)) / thisOctal%subcellSize
!                d2phidx2(2) = (g(3) - g(4)) / thisOctal%subcellSize
!                d2phidx2(3) = (g(5) - g(6)) / thisOctal%subcellSize
!                sumd2phidx2 = SUM(d2phidx2(1:3))
!             endif
!             newPhi = thisOctal%phi_i(subcell) + (deltaT * sumd2phidx2 - fourPi * gGrav * thisOctal%rho(subcell) * deltaT)
             if (thisOctal%phi_i(subcell) /= 0.d0) then
                frac = abs((newPhi - thisOctal%phi_i(subcell))/thisOctal%phi_i(subcell))
                fracChange = max(frac, fracChange)
             else
                fracChange = 1.d30
             endif
             thisOctal%phi_i(subcell) = newPhi
             
          endif
       endif
    enddo
  end subroutine gSweep

  recursive subroutine gSweep2(thisOctal, grid, deltaT, fracChange)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6), probe(6)
    integer :: n, ndir
    real(double) ::  g(6), dx
    real(double) :: deltaT, fracChange, gGrav, newPhi, frac, d2phidx2(3), sumd2phidx2
    integer :: nd

    gGrav = bigG

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call gsweep2(child, grid, deltaT, fracChange)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%twoD) then
             nDir = 4
             probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
             probe(4) = VECTOR(0.d0, 0.d0,-1.d0)
             dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
             dir(2) = VECTOR(1.d0, 0.d0, 0.d0)
             dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
             dir(4) = VECTOR(0.d0, 0.d0, 1.d0)
          else if (thisOctal%threed) then
             nDir = 6
             probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
             probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
             probe(4) = VECTOR(0.d0, 0.d0,-1.d0)
             probe(5) = VECTOR(0.d0, 1.d0, 0.d0)
             probe(6) = VECTOR(0.d0,-1.d0, 0.d0)
             dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
             dir(2) = VECTOR(1.d0, 0.d0, 0.d0)
             dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
             dir(4) = VECTOR(0.d0, 0.d0, 1.d0)
             dir(5) = VECTOR(0.d0, 1.d0, 0.d0)
             dir(6) = VECTOR(0.d0, 1.d0, 0.d0)
          endif



          if (.not.thisOctal%edgeCell(subcell)) then !xxx changed from ghostcell


             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + probe(n) * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
                
                x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)


                if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
                   x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                   x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                   dx = x2 - x1
                else
                   if (nd == thisOctal%nDepth) then ! coarse/coarse or fine/fine
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                      dx = x2 - x1
                      dx = sign(thisOctal%subcellSize, dx)
                   else if (nd > thisOctal%nDepth) then ! coarse cells with a fine boundary
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                      dx = x2 - x1
                      dx = sign(thisOctal%subcellSize, dx)
                   else
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n) ! fine cells
                      x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                      dx = x2 - x1
                      dx = sign(thisOctal%subcellSize, dx)
                   endif
                endif


                g(n) =   (phi - thisOctal%phi_i(subcell))/(dx*gridDistanceScale)

                
             enddo

             if (thisOctal%twoD) then
                d2phidx2(1) = (g(1) - g(2)) / (thisOctal%subcellSize*gridDistanceScale)
                d2phidx2(2) = (g(3) - g(4)) / (thisOctal%subcellSize*gridDistanceScale)
                sumd2phidx2 = SUM(d2phidx2(1:2))
             else
                d2phidx2(1) = (g(1) - g(2)) / (thisOctal%subcellSize*gridDistanceScale)
                d2phidx2(2) = (g(3) - g(4)) / (thisOctal%subcellSize*gridDistanceScale)
                d2phidx2(3) = (g(5) - g(6)) / (thisOctal%subcellSize*gridDistanceScale)
                sumd2phidx2 = SUM(d2phidx2(1:3))
             endif
             deltaT = (1.d0/6.d0)*(thisOctal%subcellSize*gridDistanceScale)**2

             newPhi = thisOctal%phi_i(subcell) + (deltaT * sumd2phidx2 - fourPi * gGrav * thisOctal%rho(subcell) * deltaT)

!             if (myrankglobal == 1) write(*,*) newphi,thisOctal%phi_i(subcell),deltaT, sumd2phidx2, thisOctal%rho(subcell)
             if (thisOctal%phi_i(subcell) /= 0.d0) then
                frac = abs((newPhi - thisOctal%phi_i(subcell))/thisOctal%phi_i(subcell))
                fracChange = max(frac, fracChange)
             else
                fracChange = 1.d30
             endif
             thisOctal%phi_i(subcell) = newPhi
             
          endif
       endif
    enddo
  end subroutine gSweep2

  recursive subroutine gSweepLevel(thisOctal, grid, deltaT, fracChange, nDepth)
    include 'mpif.h'
    integer :: myRank, ierr, nd
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: nDepth
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6)
    integer :: n, ndir
    real(double) :: x1, x2, g(6)
    real(double) :: deltaT, fracChange, gGrav, newPhi, frac !, d2phidx2(3), sumd2phidx2

    gGrav = bigG

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call gsweepLevel(child, grid, deltaT, fracChange, ndepth)
       end do
    else

       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%twoD) then
             nDir = 4
             dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
             dir(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
             dir(4) = VECTOR(0.d0, 0.d0,-1.d0)
          else if (thisOctal%threed) then
             nDir = 6
             dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
             dir(2) = VECTOR(-1.d0, 0.d0, 0.d0)
             dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
             dir(4) = VECTOR(0.d0, 0.d0,-1.d0)
             dir(5) = VECTOR(0.d0, 1.d0, 0.d0)
             dir(6) = VECTOR(0.d0,-1.d0, 0.d0)
          endif
          


          if (.not.thisOctal%edgeCell(subcell)) then !xxx changed from ghostcell


             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + dir(n) * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
                neighbourOctal => thisOctal
                call findSubcellLocalLevel(locator, neighbourOctal, neighbourSubcell, nDepth)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, dir(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
                
                x1 = subcellCentre(thisOctal, subcell) .dot. dir(n)
                x2 = subcellCentre(neighbourOctal, neighboursubcell) .dot. dir(n)

!                g(n) =   (phi - thisOctal%phi_i(subcell))/(x2 - x1)
!                g(n) =   (phi - thisOctal%phi_i(subcell))/thisOctal%subcellSize !!!!!!!!!!!!!!!!!!!!!!!!!!

!                if (slow) then
                   g(n) = phi
!                endif

                enddo

                if (thisOctal%twoD) then
                   newphi = 0.25d0*(SUM(g(1:4))) - fourpi*gGrav*thisOctal%rho(subcell)*deltaT
                else
                   newphi = 0.16666666666667d0*(SUM(g(1:6))) - fourpi*gGrav*thisOctal%rho(subcell)*deltaT
                endif!
                !             if (thisOctal%twoD) then
!                d2phidx2(1) = (g(1) - g(2)) / thisOctal%subcellSize
!                d2phidx2(2) = (g(3) - g(4)) / thisOctal%subcellSize
!                sumd2phidx2 = SUM(d2phidx2(1:2))
!             else
!                d2phidx2(1) = (g(1) - g(2)) / thisOctal%subcellSize
!                d2phidx2(2) = (g(3) - g(4)) / thisOctal%subcellSize
!                d2phidx2(3) = (g(5) - g(6)) / thisOctal%subcellSize
!                sumd2phidx2 = SUM(d2phidx2(1:3))
!             endif
!             newPhi = thisOctal%phi_i(subcell) + (deltaT * sumd2phidx2 - fourPi * gGrav * thisOctal%rho(subcell) * deltaT)
                if (thisOctal%phi_i(subcell) /= 0.d0) then
                   frac = abs((newPhi - thisOctal%phi_i(subcell))/thisOctal%phi_i(subcell))
                   fracChange = max(frac, fracChange)
                else
                   fracChange = 1.d30
                endif
             thisOctal%phi_i(subcell) = newPhi
             
          endif
       enddo
    endif
  end subroutine gSweepLevel
  
  subroutine selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid)
    use input_variables, only :  maxDepthAMR
    include 'mpif.h'
    type(gridtype) :: grid
    logical, optional :: multigrid
    integer, parameter :: maxThreads = 100
    integer :: iDepth
    integer :: nPairs, thread1(:), thread2(:), nBound(:), group(:), nGroup
    real(double) :: fracChange(maxthreads), tempFracChange(maxthreads), deltaT, dx
    integer :: nHydrothreads
    real(double), parameter :: tol = 1.d-5
    integer :: it, ierr
!    character(len=30) :: plotfile
    nHydroThreads = nThreadsGlobal - 1

!    if (myrankglobal == 1) call tune(6,"Complete self gravity")

    if (PRESENT(multigrid)) then

       call updateDensityTree(grid%octreeRoot)

       do iDepth = 3, maxDepthAmr-1

          dx = grid%octreeRoot%subcellSize/dble(2.d0**(iDepth-1))
          if (grid%octreeRoot%twoD) then
             deltaT =  (dx*gridDistanceScale)**2 / 4.d0
          else
             deltaT =  (dx*gridDistanceScale)**2 / 6.d0
          endif
          it = 0
          fracChange = 1.d30
          do while (ANY(fracChange(1:nHydrothreads) > tol))
             fracChange = 0.d0
             it = it + 1
             
             call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
             
             call gSweepLevel(grid%octreeRoot, grid, deltaT, fracChange(myRankGlobal), iDepth)
             call MPI_ALLREDUCE(fracChange, tempFracChange, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR, ierr)
             
             fracChange = tempFracChange
          enddo
          if (myRankGlobal == 1) write(*,*) "Gsweep of depth ", iDepth, " done in ", it, " iterations"
          call updatePhiTree(grid%octreeRoot, iDepth)
       enddo


    endif


    if (grid%octreeRoot%twoD) then
       deltaT =  (2.d0*grid%halfSmallestSubcell*gridDistanceScale)**2 / 4.d0
    else
       deltaT =  (2.d0*grid%halfSmallestSubcell*gridDistanceScale)**2 / 6.d0
    endif

    fracChange = 1.d30
    it = 0
    do while (ANY(fracChange(1:nHydrothreads) > tol))
       fracChange = 0.d0
       it = it + 1

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       
       call gSweep2(grid%octreeRoot, grid, deltaT, fracChange(myRankGlobal))
       call MPI_ALLREDUCE(fracChange, tempFracChange, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR, ierr)
       
       fracChange = tempFracChange
       !       write(plotfile,'(a,i4.4,a)') "grav",it,".png/png"
!e           if (myrankglobal == 1)   write(*,*) it,MAXVAL(fracChange(1:nHydroThreads))
    enddo
    if (myRankGlobal == 1) write(*,*) "Gravity solver completed after: ",it, " iterations"


!    if (myrankglobal == 1) call tune(6,"Complete self gravity")


  end subroutine selfGrav

  real(double) function getPressure(thisOctal, subcell, iEquationOfState, gamma)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: iEquationOfState
    real(double) :: eKinetic, eThermal, K, u2, gamma, eTot
    real(double), parameter :: gamma2 = 1.4d0, rhoCrit = 1.d-14

    select case(iEquationOfState)
       case(0) ! adiabatic
          u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
            thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
          eKinetic = u2 / 2.d0
          eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
          eThermal = eTot - eKinetic
          getPressure =  (gamma - 1.d0) * thisOctal%rho(subcell) * eThermal
       case(1) ! isothermal
          eThermal = thisOctal%rhoe(subcell) / thisOctal%rho(subcell)
          getPressure =  (gamma - 1.d0) * thisOctal%rho(subcell) * eThermal
       case(2) !  equation of state from Bonnell 1994
          if (thisOctal%rho(subcell) < rhoCrit) then
             K = 4.1317d8
             getPressure = K * thisOctal%rho(subcell)
!             write(*,*) "low density pressure called ",getPressure
          else
             K = (4.1317d8*rhoCrit)/(rhoCrit**gamma2)
             getPressure = K * thisOctal%rho(subcell)**gamma2
!             write(*,*) "high density pressure called ",getpressure
          endif

       case(3) ! polytropic
!          eThermal  = rk2 * thisOctal%rho(subcell)**(gamma-1.d0)
!          getPressure = (2.d0/3.d0) * eThermal * thisOctal%rho(subcell)
          

       case DEFAULT
          write(*,*) "Equation of state not recognised: ", iEquationOfState
          stop
    end select
                
  end function getPressure

  subroutine minMaxDepth(thisOctal, minDepth, maxDepth)
    type(octal), pointer :: thisOctal
    integer :: minDepth, maxDepth

    minDepth = 200
    maxDepth = 0
    call minMaxDepthPrivate(thisOctal, minDepth, maxDepth)

    contains

  RECURSIVE SUBROUTINE minMaxDepthPrivate(thisOctal, minDepth, maxDepth)
    !   to .FALSE.
  
    TYPE(octal), INTENT(INOUT) :: thisOctal 
    INTEGER :: iChild
    integer :: minDepth, maxDepth

    if ( (thisOctal%nChildren == 0) ) then
       minDepth = min(minDepth, thisOctal%nDepth)
       maxDepth = max(maxDepth, thisOctal%nDepth)
    endif

    DO iChild = 1, thisOctal%nChildren
      CALL minMaxDepthPrivate(thisOctal%child(iChild), minDepth, maxDepth)
    END DO

  end SUBROUTINE minMaxDepthPrivate
end subroutine minMaxDepth

  function phiInterp(rVec, thisOctal, subcell, direction) result(phi)
    real(double) :: phi
    type(VECTOR) :: rVec, direction
    type(OCTAL), pointer :: thisOctal
    integer :: subcell

    phi = 0.d0
  end function phiInterp
    

#endif
    
end module hydrodynamics_mod
