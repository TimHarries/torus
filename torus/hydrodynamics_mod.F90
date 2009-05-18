! hydrodynamics module added by tjh on 18th june 2007

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
  use gridio_mod
  use vtk_mod

  implicit none

  public
!  real(double) :: griddistancescale = 1.d10
!  real(double) :: griddistancescale = 1.d0

contains

  recursive subroutine fluxlimiter(thisoctal, limitertype)
    include 'mpif.h'
    integer :: myrank, ierr
    character(len=*) :: limitertype
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: a, b, dq
    character(len=80) :: message
 
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call fluxlimiter(child, limitertype)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle


!          if (.not.thisoctal%ghostcell(subcell)) then

          
             dq = thisoctal%q_i(subcell) - thisoctal%q_i_minus_1(subcell)
!             write(*,*) myrankglobal,dq,subcellcentre(thisoctal,subcell), &
!                  thisoctal%q_i(subcell), thisoctal%q_i_minus_1(subcell), thisoctal%u_interface(subcell)
             if (dq /= 0.d0) then
                if (thisoctal%u_interface(subcell) .ge. 0.d0) then
                   thisoctal%rlimit(subcell) = (thisoctal%q_i_minus_1(subcell) - thisoctal%q_i_minus_2(subcell)) / dq
                else
                   thisoctal%rlimit(subcell) = (thisoctal%q_i_plus_1(subcell) - thisoctal%q_i(subcell)) / dq
                endif
             else
                thisoctal%rlimit(subcell) = 1.d0
             endif
             select case (limitertype)
             case("superbee")
!                write(*,*) myrankglobal, "q_i ", thisoctal%q_i(subcell)
!                write(*,*) myrankglobal, "q_i-1", thisoctal%q_i_minus_1(subcell)
!                write(*,*) myrankglobal, "q_i-2", thisoctal%q_i_minus_2(subcell)
!                write(*,*) myrankglobal, "q_i+1", thisoctal%q_i_plus_1(subcell)
!                write(*,*) myrankglobal, "dq ",dq
!                write(*,*) myrankglobal, "r ", thisoctal%rlimit(subcell)

                a = min(1.d0, 2.d0*thisoctal%rlimit(subcell))
                b = min(2.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit(subcell) = max(0.d0, a, b)
                
             case default
                write(message,'(a,a)') "flux limiter not recognised: ", trim(limitertype)
                call writefatal(message)
                stop
             end select
             if (thisoctal%ghostcell(subcell)) thisoctal%philimit(subcell) = 1.d0
!                        thisoctal%philimit(subcell) = 1.d0
!          endif
      endif
    enddo
  end subroutine fluxlimiter


  recursive subroutine updatedensitytree(thisoctal)
    include 'mpif.h'
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal, parentoctal, testoctal
    type(octal), pointer  :: child 
    integer :: i, n, m

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)


    if (thisoctal%nchildren > 0) then
       do i = 1, thisoctal%nchildren, 1
          child => thisoctal%child(i)
          call updatedensitytree(child)
       end do
    else

          
       testoctal => thisoctal
          
       do while(testoctal%ndepth >= 2)
             
          parentoctal => testoctal%parent
          m = testoctal%parentsubcell
          if (testoctal%twod) then
             n = 4
          else
             n = 8
          endif
          parentoctal%rho(m) = sum(testoctal%rho(1:n))/dble(n)
          parentoctal%phi_i(m) = sum(testoctal%phi_i(1:n))/dble(n)
          parentoctal%edgecell(m) = any(testoctal%edgecell(1:n))
          testoctal => parentoctal
       enddo
    endif
  end subroutine updatedensitytree

  recursive subroutine updatephitree(thisoctal, ndepth)
    include 'mpif.h'
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, n, ndepth

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    if (thisoctal%ndepth == (ndepth+1)) then

       if (thisoctal%twod) then
          n = 4
       else
          n = 8
       endif

       where(.not.thisoctal%edgecell)
          thisoctal%phi_i = thisoctal%parent%phi_i(thisoctal%parentsubcell)
       end where
!       write(*,*) "phi ",thisoctal%phi_i(1:8)
    else

       do subcell = 1, thisoctal%maxchildren
          if (thisoctal%haschild(subcell)) then
             ! find the child
             do i = 1, thisoctal%nchildren, 1
                if (thisoctal%indexchild(i) == subcell) then
                   child => thisoctal%child(i)
                   call updatephitree(child, ndepth)
                   exit
                end if
             end do
          endif
       enddo
    endif

  end subroutine updatephitree


  recursive subroutine constructflux(thisoctal, dt)
    include 'mpif.h'
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call constructflux(child, dt)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             if (thisoctal%u_interface(subcell).ge.0.d0) then
                thisoctal%flux_i(subcell) = thisoctal%u_interface(subcell) * thisoctal%q_i_minus_1(subcell)
             else
                thisoctal%flux_i(subcell) = thisoctal%u_interface(subcell) * thisoctal%q_i(subcell)
             endif
             
             if (thisoctal%x_i(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) "problem with the x_i values"
                stop
             endif

             dx = (thisoctal%x_i(subcell) - thisoctal%x_i_minus_1(subcell))
             dx = thisoctal%subcellsize * griddistancescale

             thisoctal%flux_i(subcell) = thisoctal%flux_i(subcell) + &
                  0.5d0 * abs(thisoctal%u_interface(subcell)) * &
                  (1.d0 - abs(thisoctal%u_interface(subcell) * dt / &
                   dx)) * &
                  thisoctal%philimit(subcell) * (thisoctal%q_i(subcell) - thisoctal%q_i_minus_1(subcell))

          endif
!          if (thisoctal%flux_i(subcell) > 1.d-10) write(*,*) "flux ",thisoctal%flux_i(subcell),thisoctal%u_interface(subcell), &
!               thisoctal%rho(subcell), thisoctal%rhou(subcell)
       endif
    enddo
  end subroutine constructflux

  recursive subroutine updatecellq(thisoctal, dt)
    include 'mpif.h'
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call updatecellq(child, dt)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle


          if (.not.thisoctal%ghostcell(subcell)) then
             dx = thisoctal%subcellsize * griddistancescale

!             x_minus_half = 0.5d0*(thisoctal%x_i(subcell)+thisoctal%x_i_minus_1(subcell))
!             x_plus_half = 0.5d0*(thisoctal%x_i(subcell)+thisoctal%x_i_plus_1(subcell))
!             dx = x_plus_half - x_minus_half

!             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * &
!                  (thisoctal%flux_i_plus_1(subcell) - thisoctal%flux_i(subcell)) / &
!                  (thisoctal%x_i_plus_1(subcell) - thisoctal%x_i(subcell))


!             dx = thisoctal%subcellsize * griddistancescale

!             write(*,*) "old q ",thisoctal%q_i(subcell)
             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * &
                  (thisoctal%flux_i_plus_1(subcell) - thisoctal%flux_i(subcell)) / dx



!             write(*,*) "new q ", thisoctal%q_i(subcell),dt,thisoctal%flux_i_plus_1(subcell), &
!                  thisoctal%flux_i(subcell),dx
          endif
       endif
    enddo
  end subroutine updatecellq

  recursive subroutine synchronizefluxes(thisoctal, dt, idepth)
    include 'mpif.h'
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, idepth
    real(double) :: dt !, dx

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call synchronizefluxes(child, dt, idepth)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (thisoctal%ndepth /= idepth) cycle

          if (.not.thisoctal%ghostcell(subcell)) then


!             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * &
!                  (thisoctal%flux_i_plus_1(subcell) - thisoctal%flux_i(subcell)) / dx


          endif
       endif
    enddo
  end subroutine synchronizefluxes

  recursive subroutine setupx(thisoctal, grid, direction)
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(vector) :: direction
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupx(child, grid, direction)
                exit
             end if
          end do
       else
          thisoctal%x_i(subcell) = (subcellcentre(thisoctal, subcell) .dot. direction) * griddistancescale
       endif
    enddo
  end subroutine setupx


  recursive subroutine dumpx(thisoctal, n)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, n
    type(vector) :: locator
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call dumpx(child, n)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle
          if (myrankglobal==1) then
             locator = subcellcentre(thisoctal,subcell)
             n = n + 1
             write(*,*) n , locator%x
          endif
       endif
    enddo
  end subroutine dumpx

  recursive subroutine setupqx(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi
    integer :: subcell, i, neighboursubcell
    integer :: nd
    type(vector) :: direction, locator, reversedirection
  
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupqx(child, grid, direction)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

!          if (associated(thisoctal%mpiboundarystorage)) then
!             if (myrank == 1) then
!                write(*,*) "x_i", &
!                  thisoctal%x_i(subcell),thisoctal%mpiboundarystorage(subcell, 1:6, 7)
!                write(*,*) "subcell x_i ",thisoctal%x_i(subcell)
!             endif
!          endif

          thisoctal%x_i_minus_1(subcell) = 0.d0
          thisoctal%x_i_plus_1(subcell) = 0.d0
          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

             thisoctal%x_i_plus_1(subcell) = x
             thisoctal%q_i_plus_1(subcell) = q

             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%x_i_minus_1(subcell) = x
             thisoctal%q_i_minus_1(subcell) = q
             thisoctal%q_i_minus_2(subcell) = qnext


!             write(*,*) "q: ", thisoctal%q_i_minus_2(subcell), thisoctal%q_i_minus_1(subcell), thisoctal%q_i(subcell), &
!                  thisoctal%q_i_plus_1(subcell)
             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrank," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell)
             endif

          endif
       endif
    enddo
  end subroutine setupqx

  recursive subroutine setuprhophi(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator, reversedirection
    integer :: nd
  
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setuprhophi(child, grid, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

             thisoctal%rho_i_plus_1(subcell) = rho
             thisoctal%phi_i_plus_1(subcell) = phi

             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%rho_i_minus_1(subcell) = rho
             thisoctal%phi_i_minus_1(subcell) = phi



          endif
       endif
    enddo
  end subroutine setuprhophi

  recursive subroutine setupui(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, weight
    integer :: nd
  
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupui(child, grid, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle
!          if (thisoctal%mpithread(subcell) /= myrank) cycle

          if (.not.thisoctal%edgecell(subcell)) then !xxx changed from ghostcell
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)

             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

             rho_i_minus_1 = rho
             rhou_i_minus_1 = rhou

!             if (thisoctal%ndepth == nd) then
                weight = 0.5d0
!             else if (thisoctal%ndepth > nd) then ! fine to coarse
!                weight  = 0.666666666d0
!             else
!                weight  = 0.333333333d0 ! coarse to fine
!             endif


             thisoctal%u_interface(subcell) = &
                  weight*thisoctal%rhou(subcell)/thisoctal%rho(subcell) + (1.d0-weight)*rhou_i_minus_1/rho_i_minus_1

          endif
       endif
    enddo
  end subroutine setupui

  recursive subroutine setupvi(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, weight
    integer :: nd

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupvi(child, grid, direction)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)

             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

                rho_i_minus_1 = rho
                rhou_i_minus_1 = rhov

!             if (thisoctal%ndepth == nd) then
                weight = 0.5d0
!             else if (thisoctal%ndepth > nd) then ! fine to coarse
!                weight  = 0.666666666d0
!             else
!                weight  = 0.333333333d0 ! coarse to fine
!             endif

                thisoctal%u_interface(subcell) = &
                     weight*thisoctal%rhov(subcell)/thisoctal%rho(subcell) + (1.d0-weight)*rhou_i_minus_1/rho_i_minus_1

          endif
       endif
    enddo
  end subroutine setupvi

  recursive subroutine setupwi(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, weight
    integer :: nd
  
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupwi(child, grid, direction)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)

             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

                rho_i_minus_1 = rho
                rhou_i_minus_1 = rhow

!             if (thisoctal%ndepth == nd) then
                weight = 0.5d0
!             else if (thisoctal%ndepth > nd) then ! fine to coarse
!                weight  = 0.666666666d0
!             else
!                weight  = 0.333333333d0 ! coarse to fine
!             endif


                thisoctal%u_interface(subcell) = &
                     weight*thisoctal%rhow(subcell)/thisoctal%rho(subcell) + (1.d0-weight)*rhou_i_minus_1/rho_i_minus_1

          endif
       endif
    enddo
  end subroutine setupwi



  recursive subroutine setupflux(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi
    integer :: nd
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupflux(child, grid, direction)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
!             thisoctal%x_i(subcell) = (subcellcentre(thisoctal, subcell) .dot. direction) * griddistancescale
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

             if (nd >= thisoctal%ndepth) then ! this is a coarse-to-fine cell boundary
                thisoctal%flux_i_plus_1(subcell) = flux
             else
                ! now we need to do the fine-to-coarse flux

                if (thisoctal%u_i_plus_1(subcell) .ge. 0.d0) then ! flow is out of this cell into next
                   thisoctal%flux_i_plus_1(subcell) = thisoctal%q_i(subcell) * thisoctal%u_i_plus_1(subcell)
                else
                   thisoctal%flux_i_plus_1(subcell) = flux ! flow is from neighbour into this one
                endif
             endif




             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%flux_i_minus_1(subcell) = flux
          endif
       endif
    enddo
  end subroutine setupflux


  recursive subroutine setuppressure(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    integer :: nd
  
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setuppressure(child, grid, direction)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle


          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%pressure_i_plus_1(subcell) = pressure

             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%pressure_i_minus_1(subcell) = pressure

          endif
       endif
    enddo
  end subroutine setuppressure

  recursive subroutine setupupm(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    integer :: nd
  

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupupm(child, grid, direction)
                exit
             end if
          end do
       else
!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%u_i_plus_1(subcell) = rhou/rho
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%u_i_minus_1(subcell) = rhou/rho

          endif
       endif
    enddo
  end subroutine setupupm

  recursive subroutine setupvpm(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    integer :: nd
  

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupvpm(child, grid, direction)
                exit
             end if
          end do
       else
!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal =>thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%u_i_plus_1(subcell) = rhov/rho
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%u_i_minus_1(subcell) = rhov/rho

          endif
       endif
    enddo
  end subroutine setupvpm

  recursive subroutine setupwpm(thisoctal, grid, direction)
    include 'mpif.h'
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    integer :: nd
  

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupwpm(child, grid, direction)
                exit
             end if
          end do
       else
!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%u_i_plus_1(subcell) = rhow/rho
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
             thisoctal%u_i_minus_1(subcell) = rhow/rho

          endif
       endif
    enddo
  end subroutine setupwpm


   recursive subroutine computepressureu(thisoctal, direction)
     include 'mpif.h'
     integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: biggamma,eta
    type(vector) :: direction
    logical :: useviscosity

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    eta = 3.d0

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call computepressureu(child, direction)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle



          thisoctal%pressure_i(subcell) = getpressure(thisoctal, subcell)
          
          biggamma = 0.d0
          if (.not.thisoctal%edgecell(subcell)) then
             useviscosity = .false.
             if (thisoctal%u_i_plus_1(subcell) .le. thisoctal%u_i_minus_1(subcell)) useviscosity = .true.
!             if ((thisoctal%u_i_plus_1(subcell) < 0.d0).and.(thisoctal%u_i_plus_1(subcell) .ge. thisoctal%u_i_minus_1(subcell))) &
!                  useviscosity = .true.


             if (useviscosity) then
                biggamma = 0.25d0 * eta**2 * (thisoctal%u_i_plus_1(subcell) - thisoctal%u_i_minus_1(subcell))**2 &
                     * thisoctal%rho(subcell)
             else
                biggamma = 0.d0
             endif
          endif


          thisoctal%pressure_i(subcell) = thisoctal%pressure_i(subcell) + biggamma


          if (isnan(thisoctal%pressure_i(subcell))) then
             write(*,*) "pressureu has nan"
             write(*,*) thisoctal%rhou(subcell),thisoctal%rhov(subcell), thisoctal%rhow(subcell)
             write(*,*) thisoctal%rho(subcell)
             write(*,*) "cen ",subcellcentre(thisoctal, subcell)
             stop
          endif

       endif
    enddo
  end subroutine computepressureu

   recursive subroutine computepressurev(thisoctal, direction)
     include 'mpif.h'
     integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: biggamma,eta
    type(vector) :: direction

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    eta = 3.d0

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call computepressurev(child, direction)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          thisoctal%pressure_i(subcell) = getpressure(thisoctal, subcell)

          biggamma = 0.d0
          if (.not.thisoctal%edgecell(subcell)) then
             if (thisoctal%u_i_plus_1(subcell) .le. thisoctal%u_i_minus_1(subcell)) then
                biggamma = 0.25d0 * eta**2 * (thisoctal%u_i_plus_1(subcell) - thisoctal%u_i_minus_1(subcell))**2 &
                     * thisoctal%rho(subcell)
             else
                biggamma = 0.d0
             endif
          endif



          thisoctal%pressure_i(subcell) = thisoctal%pressure_i(subcell) + biggamma

          if (thisoctal%pressure_i(subcell) < 0.d0) then
             write(*,*) "pressure ", thisoctal%pressure_i(subcell)
             write(*,*) "rho ", thisoctal%rho(subcell)
             write(*,*) "rhoe ", thisoctal%rhoe(subcell)
             write(*,*) "rhoe ", thisoctal%energy(subcell)
             write(*,*) "biggamma ", biggamma
             write(*,*) "eos ", thisoctal%iequationofstate(subcell)
          endif
          if (isnan(thisoctal%pressure_i(subcell))) then
             write(*,*) "pressurev has nan"
             write(*,*) thisoctal%rhou(subcell),thisoctal%rhov(subcell), thisoctal%rhow(subcell)
             write(*,*) thisoctal%rho(subcell)
             stop
          endif
       endif
    enddo
  end subroutine computepressurev

   recursive subroutine computepressurew(thisoctal, direction)
     include 'mpif.h'
     integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: biggamma,eta
    type(vector) :: direction


    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    eta = 3.d0

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call computepressurew(child, direction)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          thisoctal%pressure_i(subcell) = getpressure(thisoctal, subcell)

          biggamma = 0.d0
          if (.not.thisoctal%edgecell(subcell)) then
             if (thisoctal%u_i_plus_1(subcell) .le. thisoctal%u_i_minus_1(subcell)) then
                biggamma = 0.25d0 * eta**2 * (thisoctal%u_i_plus_1(subcell) - thisoctal%u_i_minus_1(subcell))**2 &
                     * thisoctal%rho(subcell)
             else
                biggamma = 0.d0
             endif
          endif
!          if (myrankglobal == 1) write(*,*) biggamma/thisoctal%pressure_i(subcell), thisoctal%u_i_plus_1(subcell), &
!               thisoctal%u_i_minus_1(subcell), thisoctal%rho(subcell)


          thisoctal%pressure_i(subcell) = thisoctal%pressure_i(subcell) + biggamma

          if (isnan(thisoctal%pressure_i(subcell))) then
             write(*,*) "pressurew has nan"
             write(*,*) thisoctal%rhou(subcell),thisoctal%rhov(subcell), thisoctal%rhow(subcell)
             write(*,*) thisoctal%rho(subcell)
             stop
          endif
       endif
    enddo
  end subroutine computepressurew


  recursive subroutine pressureforceu(thisoctal, dt)
    include 'mpif.h'
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, rhou, dx


    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call pressureforceu(child, dt)
                exit
             end if
          end do
       else
!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle


          if (.not.thisoctal%ghostcell(subcell)) then


             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrank," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell), thisoctal%x_i(subcell)
                write(*,*) thisoctal%ndepth
                write(*,*) "centre ",subcellcentre(thisoctal,subcell)
                write(*,*) "thread ",thisoctal%mpithread(subcell)
             endif

             rhou = thisoctal%rhou(subcell)

             dx = (thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
             dx = thisoctal%subcellsize * griddistancescale


             thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - (dt/2.d0) * &!!!!!!!!!!!!!!!!!!!!!!!
                  (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx

             thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - (dt/2.d0) * & !gravity
                  thisoctal%rho(subcell) *(thisoctal%phi_i_plus_1(subcell) - &
                  thisoctal%phi_i_minus_1(subcell)) / dx

             if (isnan(thisoctal%rhou(subcell))) then
                write(*,*) "rhou ",thisoctal%rhou(subcell)
                write(*,*) "dt ",dt
                write(*,*) "phi ",thisoctal%phi_i_plus_1(subcell), thisoctal%phi_i_minus_1(subcell)
                write(*,*) "x ",thisoctal%x_i_plus_1(subcell), thisoctal%x_i_minus_1(subcell)
             endif

             if (thisoctal%rhou(subcell)/thisoctal%rho(subcell) > 1.d10) then
                write(*,*) "u ",thisoctal%rhou(subcell)/thisoctal%rho(subcell)/1.d5
                write(*,*) "rho ", thisoctal%rho(subcell)
                write(*,*) "rho i+1", thisoctal%rho_i_plus_1(subcell)
                write(*,*) "rho i-1", thisoctal%rho_i_minus_1(subcell)
                write(*,*) "pressure i + 1 " ,thisoctal%pressure_i_plus_1(subcell)
                write(*,*) "pressure i - 1 " ,thisoctal%pressure_i_minus_1(subcell)
                write(*,*) "x i + 1 ", thisoctal%x_i_plus_1(subcell)
                write(*,*) "x i - 1 ", thisoctal%x_i_minus_1(subcell)
             endif

                
!             write(*,*) thisoctal%rhou(subcell), thisoctal%pressure_i_plus_1(subcell), &
!                  thisoctal%pressure_i_minus_1(subcell), thisoctal%x_i_plus_1(subcell), &
!                  thisoctal%x_i_minus_1(subcell)

!             if (iequationofstate == 0) then


             if (thisoctal%iequationofstate(subcell) /= 1) then 
                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - (dt/2.d0) * &!!!!!!!!!!!
                     (thisoctal%pressure_i_plus_1(subcell) * thisoctal%u_i_plus_1(subcell) - &
                     thisoctal%pressure_i_minus_1(subcell) * thisoctal%u_i_minus_1(subcell)) / dx
                     

                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - (dt/2.d0) * & !gravity
                  rhou  * (thisoctal%phi_i_plus_1(subcell) - thisoctal%phi_i_minus_1(subcell)) / dx
             endif

!             endif
             if (isnan(thisoctal%rhou(subcell))) then
                write(*,*) "bug",thisoctal%rhou(subcell), &
                     thisoctal%pressure_i_plus_1(subcell),thisoctal%pressure_i_minus_1(subcell)
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureforceu

  recursive subroutine pressureforcev(thisoctal, dt)
    include 'mpif.h'
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, rhou, dx


    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call pressureforcev(child, dt)
                exit
             end if
          end do
       else
!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle


          if (.not.thisoctal%ghostcell(subcell)) then


             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrank," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell), thisoctal%x_i(subcell)
                write(*,*) thisoctal%ndepth
                write(*,*) "centre ",subcellcentre(thisoctal,subcell)
                write(*,*) "thread ",thisoctal%mpithread(subcell)
             endif

             rhou = thisoctal%rhov(subcell)

             dx = (thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
             dx = thisoctal%subcellsize * griddistancescale


             thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - (dt/2.d0) * &
                  (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx

             thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - (dt/2.d0) * & !gravity
                  thisoctal%rho(subcell) *(thisoctal%phi_i_plus_1(subcell) - &
                  thisoctal%phi_i_minus_1(subcell)) / dx

             if (isnan(thisoctal%rhov(subcell))) then
                write(*,*) "rhov ",thisoctal%rhov(subcell)
                write(*,*) "dt ",dt
                write(*,*) "phi ",thisoctal%phi_i_plus_1(subcell), thisoctal%phi_i_minus_1(subcell)
                write(*,*) "x ",thisoctal%x_i_plus_1(subcell), thisoctal%x_i_minus_1(subcell)
             endif

             if (thisoctal%rhov(subcell)/thisoctal%rho(subcell) > 1.d10) then
                write(*,*) "u ",thisoctal%rhov(subcell)/thisoctal%rho(subcell)/1.d5
                write(*,*) "rho ", thisoctal%rho(subcell)
                write(*,*) "rho i+1", thisoctal%rho_i_plus_1(subcell)
                write(*,*) "rho i-1", thisoctal%rho_i_minus_1(subcell)
                write(*,*) "pressure i + 1 " ,thisoctal%pressure_i_plus_1(subcell)
                write(*,*) "pressure i - 1 " ,thisoctal%pressure_i_minus_1(subcell)
                write(*,*) "x i + 1 ", thisoctal%x_i_plus_1(subcell)
                write(*,*) "x i - 1 ", thisoctal%x_i_minus_1(subcell)
             endif

                
!             write(*,*) thisoctal%rhou(subcell), thisoctal%pressure_i_plus_1(subcell), &
!                  thisoctal%pressure_i_minus_1(subcell), thisoctal%x_i_plus_1(subcell), &
!                  thisoctal%x_i_minus_1(subcell)

!             if (iequationofstate == 0) then


             if (thisoctal%iequationofstate(subcell) /= 1) then 
                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - (dt/2.d0) * &
                     (thisoctal%pressure_i_plus_1(subcell) * thisoctal%u_i_plus_1(subcell) - &
                     thisoctal%pressure_i_minus_1(subcell) * thisoctal%u_i_minus_1(subcell)) / dx
                     

                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - (dt/2.d0) * & !gravity
                  rhou  * (thisoctal%phi_i_plus_1(subcell) - thisoctal%phi_i_minus_1(subcell)) / dx

             endif
             if (isnan(thisoctal%rhov(subcell))) then
                write(*,*) "bug",thisoctal%rhov(subcell), &
                     thisoctal%pressure_i_plus_1(subcell),thisoctal%pressure_i_minus_1(subcell)
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureforcev


  recursive subroutine pressureforcew(thisoctal, dt)
    include 'mpif.h'
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, rhow, dx

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call pressureforcew(child, dt)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle


          if (.not.thisoctal%ghostcell(subcell)) then


             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrank," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell), thisoctal%x_i(subcell)
                write(*,*) thisoctal%ndepth
                write(*,*) "centre ",subcellcentre(thisoctal,subcell)
                write(*,*) "thread ",thisoctal%mpithread(subcell)
             endif

             rhow = thisoctal%rhow(subcell)


             dx = (thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
             dx = thisoctal%subcellsize * griddistancescale

             thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - (dt/2.d0) * &
                  (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx

             thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - (dt/2.d0) * & !gravity
                  thisoctal%rho(subcell) *(thisoctal%phi_i_plus_1(subcell) - &
                  thisoctal%phi_i_minus_1(subcell)) / dx


!             write(*,*) thisoctal%rhou(subcell), thisoctal%pressure_i_plus_1(subcell), &
!                  thisoctal%pressure_i_minus_1(subcell), thisoctal%x_i_plus_1(subcell), &
!                  thisoctal%x_i_minus_1(subcell)

             if (thisoctal%iequationofstate(subcell) /= 1) then
                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - (dt/2.d0) * &
                     (thisoctal%pressure_i_plus_1(subcell) * thisoctal%u_i_plus_1(subcell) - &
                     thisoctal%pressure_i_minus_1(subcell) * thisoctal%u_i_minus_1(subcell)) / dx

                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - (dt/2.d0) * & !gravity
                  rhow * (thisoctal%phi_i_plus_1(subcell) - thisoctal%phi_i_minus_1(subcell)) / dx
             endif


             if (isnan(thisoctal%rhow(subcell))) then
                write(*,*) "bug",thisoctal%pressure_i_plus_1(subcell),thisoctal%pressure_i_minus_1(subcell)
                do;enddo
                endif
             endif
       endif
    enddo
  end subroutine pressureforcew




  recursive subroutine copyrhotoq(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyrhotoq(child)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%rho(subcell)
        
       endif
    enddo
  end subroutine copyrhotoq

  recursive subroutine copyrhoetoq(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyrhoetoq(child)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%rhoe(subcell)
        
       endif
    enddo
  end subroutine copyrhoetoq

  recursive subroutine copyrhoutoq(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyrhoutoq(child)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%rhou(subcell)
        
       endif
    enddo
  end subroutine copyrhoutoq

  recursive subroutine copyrhovtoq(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyrhovtoq(child)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%rhov(subcell)
        
       endif
    enddo
  end subroutine copyrhovtoq

  recursive subroutine copyrhowtoq(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyrhowtoq(child)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%rhow(subcell)
        
       endif
    enddo
  end subroutine copyrhowtoq

  recursive subroutine copyqtorhov(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtorhov(child)
                exit
             end if
          end do
       else
  
          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rhov(subcell) = thisoctal%q_i(subcell)
          endif
        
       endif
    enddo
  end subroutine copyqtorhov

  recursive subroutine copyqtorhow(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtorhow(child)
                exit
             end if
          end do
       else
  
          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rhow(subcell) = thisoctal%q_i(subcell)
          endif
        
       endif
    enddo
  end subroutine copyqtorhow

  recursive subroutine copyqtorho(thisoctal, direction)
    type(vector) :: direction
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtorho(child, direction)
                exit
             end if
          end do
       else
  
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle

          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rho(subcell) = thisoctal%q_i(subcell)
          endif
          if (thisoctal%rho(subcell) < 0.d0) then
             write(*,*) "rho warning ", thisoctal%rho(subcell), subcellcentre(thisoctal, subcell)
             write(*,*) "label ",thisoctal%label
             write(*,*) "direction ", direction
             write(*,*) "phi limit ", thisoctal%philimit(subcell)
             write(*,*) "flux ", thisoctal%flux_i_plus_1(subcell), thisoctal%flux_i(subcell), &
                  thisoctal%flux_i_minus_1(subcell)
             write(*,*) "u ", thisoctal%u_i_plus_1(subcell)/1.e5, thisoctal%u_interface(subcell)/1.e5, &
                  thisoctal%u_i_minus_1(subcell)/1.e5
             write(*,*) "x ", thisoctal%x_i_plus_1(subcell), thisoctal%x_i(subcell), thisoctal%x_i_minus_1(subcell)
             write(*,*) "dx ", thisoctal%x_i_plus_1(subcell) - thisoctal%x_i(subcell), thisoctal%x_i(subcell) - &
                  thisoctal%x_i_minus_1(subcell)
             if (abs(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i(subcell)) > &
                 abs(thisoctal%x_i_minus_1(subcell) - thisoctal%x_i(subcell))) then
                write(*,*) "fine to coarse"
             endif
             write(*,*) "rho ", thisoctal%rho_i_plus_1(subcell), thisoctal%rho(subcell), thisoctal%rho_i_minus_1(subcell)
             write(*,*) "uinterface ", thisoctal%u_interface(subcell)
             write(*,*) "q ", thisoctal%q_i_plus_1(subcell), thisoctal%q_i(subcell), thisoctal%q_i_minus_1(subcell)
             if (.not.associated(thisoctal%mpiboundarystorage)) &
                write(*,*) "cell is not on a boundary"

          endif
       endif
    enddo
  end subroutine copyqtorho

  recursive subroutine copyqtorhoe(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtorhoe(child)
                exit
             end if
          end do
       else
  
          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rhoe(subcell) = thisoctal%q_i(subcell)
          endif
        
       endif
    enddo
  end subroutine copyqtorhoe

  recursive subroutine copyqtorhou(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtorhou(child)
                exit
             end if
          end do
       else
  
          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rhou(subcell) = thisoctal%q_i(subcell)
             if (isnan(thisoctal%rhou(subcell))) write(*,*) "nan in q to rhou"
          endif
        
       endif
    enddo
  end subroutine copyqtorhou


  subroutine advectrho(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhotoq(grid%octreeroot)
    call advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorho(grid%octreeroot, direction)

  end subroutine advectrho

  subroutine advectrhoe(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoetoq(grid%octreeroot)
    call advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhoe(grid%octreeroot)

  end subroutine advectrhoe

  subroutine advectrhou(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoutoq(grid%octreeroot)
    call advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhou(grid%octreeroot)

  end subroutine advectrhou

  subroutine advectrhov(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhovtoq(grid%octreeroot)
    call advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhov(grid%octreeroot)

  end subroutine advectrhov

  subroutine advectrhow(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhowtoq(grid%octreeroot)
    call advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhow(grid%octreeroot)

  end subroutine advectrhow

  subroutine advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call setupx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupqx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)

    call fluxlimiter(grid%octreeroot, "superbee")
    call constructflux(grid%octreeroot, dt)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupflux(grid%octreeroot, grid, direction)
    call updatecellq(grid%octreeroot, dt)

    
  end subroutine advectq


  subroutine  hydrostep1d(grid, dt, npairs, thread1, thread2, nbound, group, ngroup)
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    integer :: npairs, thread1(:), thread2(:), nbound(:), group(:), ngroup

    direction = vector(1.d0, 0.d0, 0.d0)

    call imposeboundary(grid%octreeroot)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)

    direction = vector(1.d0, 0.d0, 0.d0)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
 
    call setupui(grid%octreeroot, grid, direction)
    call setupupm(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
!
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call setupui(grid%octreeroot, grid, direction)
    call setupupm(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call computepressureu(grid%octreeroot, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call setuppressure(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call pressureforceu(grid%octreeroot, dt)
!
    call imposeboundary(grid%octreeroot)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)

  end subroutine hydrostep1d


  recursive subroutine fullstep3d(grid,timestep, npairs, thread1, thread2, &
       nbound, group, ngroup, idepth, mindepth, maxdepth)
    type(gridtype) :: grid
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: idepth
    real(double) :: timestep, nexttimestep
    integer :: mindepth, maxdepth, nextdepth


    if (idepth /= maxdepth) then
       nexttimestep = timestep / 2.d0
       nextdepth = idepth + 1
       call fullstep3d(grid,  nexttimestep, npairs, thread1, thread2, &
            nbound, group, ngroup, nextdepth, mindepth, maxdepth)
       call fullstep3d(grid, nexttimestep, npairs, thread1, thread2, &
            nbound, group, ngroup, nextdepth, mindepth, maxdepth)
    endif


  end subroutine fullstep3d

    


  subroutine hydrostep3d(grid, dt, nPairs, thread1, thread2, nBound, &
       group, nGroup,doSelfGrav)
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    logical, optional :: doSelfGrav
    logical :: selfGravity
    integer :: group(:), nGroup
    real(double) :: dt
    type(VECTOR) :: direction
    direction = VECTOR(1.d0, 0.d0, 0.d0)

    selfGravity = .true.
    if (PRESENT(doSelfGrav)) selfgravity = doSelfGrav

    if (myrankglobal == 1) call tune(6,"Boundary conditions")
    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)
    if (myrankglobal == 1) call tune(6,"Boundary conditions")


    if (myrankglobal == 1) call tune(6,"X-direction step")
    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoV(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call computePressureU(grid%octreeRoot, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call pressureForceU(grid%octreeRoot, dt/2.d0)
    if (myrankglobal == 1) call tune(6,"X-direction step")



    if (myrankglobal == 1) call tune(6,"Y-direction step")
    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call setupVi(grid%octreeRoot, grid, direction)
    call setupVpm(grid%octreeRoot, grid, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call setupVi(grid%octreeRoot, grid, direction)
    call setupVpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call computePressureV(grid%octreeRoot, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call pressureForceV(grid%octreeRoot, dt)
    if (myrankglobal == 1) call tune(6,"Y-direction step")


    if (myrankglobal == 1) call tune(6,"Z-direction step")
    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call computePressureW(grid%octreeRoot, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call pressureForceW(grid%octreeRoot, dt)
    if (myrankglobal == 1) call tune(6,"Z-direction step")



    if (myrankglobal == 1) call tune(6,"X-direction step")
    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoV(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call computePressureU(grid%octreeRoot, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call pressureForceU(grid%octreeRoot, dt/2.d0)
    if (myrankglobal == 1) call tune(6,"X-direction step")



    if (myrankglobal == 1) call tune(6,"Boundary conditions")
    call periodBoundary(grid)
    call imposeBoundary(grid%octreeRoot)
    call transferTempStorage(grid%octreeRoot)
    if (myrankglobal == 1) call tune(6,"Boundary conditions")



    if (myrankglobal == 1) call tune(6,"Self-gravity")
    if (selfGravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    if (myrankglobal == 1) call tune(6,"Self-gravity")


  end subroutine hydroStep3d

  subroutine hydroStep2d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup

    real(double) :: dt
    type(VECTOR) :: direction


    direction = VECTOR(1.d0, 0.d0, 0.d0)

    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)



    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call computePressureU(grid%octreeRoot, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call pressureForceU(grid%octreeRoot, dt/2.d0)

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call computePressureW(grid%octreeRoot, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call pressureForceW(grid%octreeRoot, dt)


    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, usethisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call computePressureU(grid%octreeRoot, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call pressureForceU(grid%octreeRoot, dt/2.d0)


    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)

 
  end subroutine hydroStep2d


  recursive subroutine computeCourantTime(thisOctal, tc)
    include 'mpif.h'
    integer :: myRank, ierr
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: tc, dx, cs, speed
  
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computeCourantTime(child, tc)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then

             cs = soundSpeed(thisOctal, subcell)
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

  function soundSpeed(thisOctal, subcell) result (cs)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: cs
    real(double) :: u2, eKinetic, eTot, eThermal
    logical, save :: firstTime = .true.
    real(double), parameter :: gamma2 = 1.4d0, rhoCrit = 1.d-14
    select case(thisOctal%iEquationOfState(subcell))
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
          cs = sqrt(thisOctal%gamma(subcell)*(thisOctal%gamma(subcell)-1.d0)*eThermal)
       case(1) ! isothermal
          cs = sqrt(getPressure(thisOctal, subcell)/thisOctal%rho(subcell))
       case(2) ! barotropic
          if (thisOctal%rho(subcell) < rhoCrit) then
             cs = sqrt( getPressure(thisOctal, subcell)/thisOctal%rho(subcell))
          else
             cs = sqrt(gamma2 *  getPressure(thisOctal, subcell)/thisOctal%rho(subcell))
          endif
       case(3) ! polytropic
          cs = sqrt(thisOctal%gamma(subcell)*thisOctal%pressure_i(subcell)/thisOctal%rho(subcell))
       case DEFAULT
          write(*,*) "Unknown equation of state passed to sound speed ", thisOctal%iEquationofState(subcell)
          
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
    logical :: globalConverged(10)
    logical :: tConverged(10)
    integer :: nHydroThreads
    real(double) :: tc(10)
    integer :: it
    integer :: myRank, ierr
    integer :: nPairs, thread1(100), thread2(100), group(100), nBound(100), ngroup
    integer :: iUnrefine, nUnrefine
    real(double) :: nextDumpTime, tdump, temptc(10), tend

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    mu = 2.d0
    nHydroThreads = nThreadsGlobal - 1


    direction = VECTOR(1.d0, 0.d0, 0.d0)

    mu = 2.d0


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    if (myRank == 1) write(*,*) "CFL set to ", cflNumber

    if (myrankGlobal /= 0) then
       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       do i = 1, nPairs
          if (myrankglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i), " group ", group(i)
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
          call evenUpGridMPI(grid,.false.,.true.)
       endif


       globalConverged = .false.
       call writeInfo("Refining grid part 2", TRIVIAL)    
       do
          globalConverged(myRank) = .true.
          call refineGridGeneric2(grid%octreeRoot, grid, globalConverged(myRank), inheritval=.false.)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreads, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
          if (ALL(tConverged(1:nHydroThreads))) exit
       end do

       call writeInfo("Evening up grid", TRIVIAL)    
       call evenUpGridMPI(grid, .false.,.true.)
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

    endif

    currentTime = 0.d0
    it = 0
    nextDumpTime = 0.d0
    tDump = 0.005d0
!    tdump = 0.2d0
    tend = 0.2d0
    iUnrefine = 0
    do while(currentTime <= tend)
       tc = 0.d0
       if (myrank /= 0) then
          tc(myrank) = 1.d30
          call computeCourantTime(grid%octreeRoot, tc(myRank))
       endif
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       !       write(*,*) "tc", tc(1:8)
       !       write(*,*) "temp tc",temptc(1:8)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)

       if (myrank == 1) write(*,*) "courantTime", dt
       if (myrank == 1) call tune(6,"Hydrodynamics step")
       call writeInfo("calling hydro step",TRIVIAL)

       if (myrankGlobal /= 0) then
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

          call hydroStep1d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)

          iUnrefine = iUnrefine + 1
          if (iUnrefine == 5) then
             if (myrankglobal == 1) call tune(6, "Unrefine grid")
             nUnrefine = 0 
!             call unrefineCells(grid%octreeRoot, grid, nUnrefine)
             !          write(*,*) "Unrefined ", nUnrefine, " cells"
             if (myrankglobal == 1) call tune(6, "Unrefine grid")
             iUnrefine = 0
          endif

          if (myrank == 1) call tune(6,"Hydrodynamics step")
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

          call writeInfo("Refining grid part 2", TRIVIAL)    
          globalConverged = .false.
          do
             globalConverged(myRank) = .true.
             call refineGridGeneric2(grid%octreeRoot, grid, globalConverged(myRank), inheritval=.false.)
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
             call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreads, MPI_LOGICAL, MPI_LOR,amrCOMMUNICATOR, ierr)
             if (ALL(tConverged(1:nHydroThreads))) exit
          end do


          call evenUpGridMPI(grid, .true., .true.)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       endif

       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt
       if (currentTime .gt. nextDumpTime) then
          call  dumpValuesAlongLine(grid, "sod.dat", VECTOR(0.d0,0.d0,0.0d0), VECTOR(1.d0, 0.d0, 0.0d0), 1000)
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
          grid%iDump = it
          grid%currentTime = currentTime
       endif


    enddo

  end subroutine doHydrodynamics1d


  subroutine doHydrodynamics3d(grid)
    include 'mpif.h'
    type(gridtype) :: grid
    real(double) :: dt, tc(64), temptc(64),  mu
    real(double) :: currentTime !, smallTime
    integer :: i, it, iUnrefine
    integer :: myRank, ierr
    character(len=80) :: plotfile
    real(double) :: tDump, nextDumpTime, tff,totalMass !, ang
    type(VECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: nSource = 0
    integer :: nUnrefine
    integer :: thread1(200), thread2(200), nBound(200), nPairs
    logical :: globalConverged(64), tConverged(64)
    integer :: nHydroThreads
    integer :: nGroup, group(200)
    logical :: doRefine
    logical :: doSelfGrav
    logical, save  :: firstStep = .true.

    doSelfGrav = .true.

    if (grid%geometry == "shakara") doSelfGrav = .false.

    dorefine = .true.

    nHydroThreads = nThreadsGlobal - 1

    direction = VECTOR(1.d0, 0.d0, 0.d0)



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

    tff = sqrt((3.d0*pi)/(32.d0*bigG*3.466d-21))
    tDump = 0.01d0 * tff



    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)



    if (myrankGlobal /= 0) then

       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       do i = 1, nPairs
          if (myrankglobal==1) &
            write(*,'(a,i3,i3,a,i3,a,i3,a,i3)') "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i), " group ",group(i)
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

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call writeInfo("Refining grid part 2", TRIVIAL)    
       globalConverged = .false.
       do
          globalConverged(myRank) = .true.
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call refineGridGeneric2(grid%octreeRoot, grid, globalConverged(myRank), inheritval=.true.)
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreads, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
          if (ALL(tConverged(1:nHydroThreads))) exit
       end do


       call writeInfo("Evening up grid", TRIVIAL)    
       call evenUpGridMPI(grid,.false., dorefine)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)



       if (myrank == 1) call tune(6, "Initial refine")


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


       if (doselfGrav) then
          if (myrank == 1) call tune(6, "Self-Gravity")
          if (myrank == 1) write(*,*) "Doing multigrid self gravity"
          call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.) 
          if (myrank == 1) write(*,*) "Done"
          if (myrank == 1) call tune(6, "Self-Gravity")
       endif

    endif

    tc = 0.d0
    if (myrank /= 0) then
    tc(myrank) = 1.d30
       call computeCourantTime(grid%octreeRoot, tc(myRank))
    endif
    call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)


    dt = MINVAL(temptc(1:nHydroThreads)) * dble(cflNumber)

    if (grid%geometry == "shakara") then
       tdump =  dt
    endif

    if (myRank == 1) write(*,*) "CFL set to ", cflNumber

    if (myrankglobal==1) write(*,*) "Setting tdump to: ", tdump
    nextDumpTime = tdump + currentTime
    iUnrefine = 0


    write(plotfile,'(a)') "start.vtk"
    call writeVtkFile(grid, plotfile, &
         valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          " /))

    do while(.true.)
       tc = 0.d0
       tc(myrank) = 1.d30
       if (myrank /= 0) call computeCourantTime(grid%octreeRoot, tc(myRank))
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       !       write(*,*) "tc", tc(1:8)
       !       if (myrank==1)write(*,*) "temp tc",temptc(1:nHydroThreads)

       if (firstStep) then 
          cflNumber = 1.d-2
          firstStep = .false.
       else
          cflNumber = 0.3d0
       endif

       dt = MINVAL(temptc(1:nHydroThreads)) * dble(cflNumber)

       if (myrank==1)write(*,*) "Current time is ",currentTime/tff, " free-fall times"
       call findMassoverAllThreads(grid, totalmass)
       if (myrank==1)write(*,*) "Current mass: ",totalmass/msol

       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

       if (myrankGlobal /= 0) then

          !       smallTime = 0.3d0 * 1.d4*gridDistanceScale/2.e5
          !       write(*,*) "myrank ",myrank, "smallest cell ",grid%halfSmallestSubcell
          if (myrank == 1) write(*,*) "courantTime", dt,cflnumber
          !       if (myrank == 1) write(*,*) "smallTime", smallTime
          if (myrank == 1) call tune(6,"Hydrodynamics step")
          !       if (dt > smallTime) dt = smallTime
          call writeInfo("calling hydro step",TRIVIAL)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


          call hydroStep3d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup, doSelfGrav=doSelfGrav)

          if (myrank == 1) call tune(6,"Hydrodynamics step")


          iUnrefine = iUnrefine + 1
          if (iUnrefine == 5) then
             if (myrank == 1)call tune(6, "Unrefine grid")
             nUnrefine = 0
             call unrefineCells(grid%octreeRoot, grid, nUnrefine)
             !          write(*,*) "Unrefined ", nUnrefine, " cells"
             if (myrank == 1)call tune(6, "Unrefine grid")
             iUnrefine = 0
          endif


          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call writeInfo("Refining grid part 2", TRIVIAL)    
          globalConverged = .false.
          do
             globalConverged(myRank) = .true.
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
             call refineGridGeneric2(grid%octreeRoot, grid, globalConverged(myRank), inheritval=.false.)
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
             call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreads, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
             if (ALL(tConverged(1:nHydroThreads))) exit
          end do

          !          
          call evenUpGridMPI(grid, .true., dorefine)

          !       if (doSelfGrav) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)


          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

          if (myrank == 1) call tune(6, "Loop refine")
          !


       endif

       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt,nextDumpTime
       if (myRank == 1) write(*,*) "percent to next dump ",100.*(nextDumpTime-currentTime)/tdump


       if (currentTime .ge. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1

          grid%iDump = it
          grid%currentTime = currentTime

          if (myrankGlobal /= 0) then
             write(plotfile,'(a,i4.4,a)') "radial",it,".dat"
             call  dumpValuesAlongLine(grid, plotfile, VECTOR(0.d0,0.d0,0.0d0), &
                  VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.0d0), 1000)
          endif

          write(plotfile,'(a,i4.4,a)') "dump",it,".grid"
          call writeAMRgrid(plotfile,.false. ,grid)

          write(plotfile,'(a,i4.4,a)') "output",it,".vtk"
          call writeVtkFile(grid, plotfile, &
               valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          " /))
          if (myrank==1) write(*,*) trim(plotfile), " written at ",currentTime/tff, " free-fall times"
       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)


    enddo
  end subroutine doHydrodynamics3d

  subroutine doHydrodynamics2d(grid)
    include 'mpif.h'
    type(gridtype) :: grid
    real(double) :: dt, tc(64), temptc(64), mu
    real(double) :: currentTime
    integer :: i, it, iUnrefine
    integer :: myRank, ierr
    character(len=80) :: plotfile
    real(double) :: tDump, nextDumpTime, tff !, ang
    type(VECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: nSource = 0
    integer :: thread1(100), thread2(100), nBound(100), nPairs
    integer :: nGroup, group(100)

    logical :: globalConverged(64), tConverged(64)
    integer :: nHydroThreads 
    logical :: dorefine
    integer :: nUnrefine, jt


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
       nHydroThreads = nThreadsGlobal - 1



    if (myrankGlobal /= 0) then

       dorefine = .true.



       direction = VECTOR(1.d0, 0.d0, 0.d0)
       mu = 2.d0

       viewVec = VECTOR(-1.d0,0.d0,0.d0)
       !    viewVec = rotateZ(viewVec, 20.d0*degtorad)
       viewVec = rotateY(viewVec, 25.d0*degtorad)


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


       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       !    call writeInfo("Refining grid", TRIVIAL)
       !    do
       !       gridConverged = .true.
       !       if (myrank == 1) call tune(6, "one refine")
       !       call refineGridGeneric2(grid%octreeRoot, grid, gridconverged, inherit=.false.)
       !       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       !       if (myrank == 1) call tune(6, "one refine")
       !       if (gridConverged) exit
       !    end do
       !    call MPI_BARRIER(amrCOMMUNICATOR, ierr)

       call writeInfo("Refining grid part 2", TRIVIAL)    
       globalConverged = .false.
       do
          globalConverged(myRank) = .true.
          call refineGridGeneric2(grid%octreeRoot, grid, globalConverged(myRank), inheritval=.false.)
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

    endif

       currentTime = 0.d0
       it = 0
       nextDumpTime = 0.d0
    call writeVTKfile(grid, "start.vtk")

    tc = 0.d0
    if (myrank /= 0) then
       tc(myrank) = 1.d30
       call computeCourantTime(grid%octreeRoot, tc(myRank))
    endif
    call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierr)
    tc = tempTc
    dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)

    tff = 1.d0 / sqrt(bigG * (0.1d0*mSol/((4.d0/3.d0)*pi*7.d15**3)))
    tDump = 1.d-2 * tff
    
    tdump = 5.d0 * dt

    write(*,*) "Setting tdump to: ", tdump

    iUnrefine = 0
       !    call writeInfo("Plotting grid", TRIVIAL)    


    jt = 0

    do while(.true.)

       jt = jt + 1
       tc = 0.d0
       if (myrank /= 0) then
          tc(myrank) = 1.d30
          call computeCourantTime(grid%octreeRoot, tc(myRank))
       endif
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierr)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)
       
       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

       if (myrankGlobal /= 0) then

          if (myrank == 1) write(*,*) "courantTime", dt
          if (myrank == 1) call tune(6,"Hydrodynamics step")
          call writeInfo("calling hydro step",TRIVIAL)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


          call hydroStep2d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)

          iUnrefine = iUnrefine + 1
          if (iUnrefine == 5) then
             if (myrankglobal == 1) call tune(6, "Unrefine grid")
             call unrefineCells(grid%octreeRoot, grid, nUnrefine)
             if (myrankglobal == 1) call tune(6, "Unrefine grid")
             iUnrefine = 0
          endif


          if (myrank == 1) call tune(6,"Hydrodynamics step")

          !       call writeInfo("Refining grid", TRIVIAL)
          !       do
          !          gridConverged = .true.
          !          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          !          call refineGridGeneric2(grid%octreeRoot, grid, gridconverged, inherit=.true.)
          !          if (gridConverged) exit
          !       end do
          !       call MPI_BARRIER(amrCOMMUNICATOR, ierr)


          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call writeInfo("Refining grid part 2", TRIVIAL)    
          globalConverged = .false.
          do
             globalConverged(myRank) = .true.
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
             call refineGridGeneric2(grid%octreeRoot, grid, globalConverged(myRank), inheritval=.true.)
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
             call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreads, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
             if (ALL(tConverged(1:nHydroThreads))) exit
          end do


          call evenUpGridMPI(grid, .true., dorefine)!, dumpfiles=jt)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)




       endif

       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt


       if (currentTime .ge. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
          grid%iDump = it
          grid%currentTime = currentTime
          !          write(plotfile,'(a,i4.4,a,i3.3,a)') "dump",it,"_rank_",myrankglobal,".grid"
          !          call writeAMRgrid(plotfile,.false. ,grid)
          write(plotfile,'(a,i4.4,a)') "output",it,".vtk"
          call writeVtkFile(grid, plotfile, &
               valueTypeString=(/"rho          ","velocity     ","rhoe         " ,"u_i          ","hydrovelocity" /))

       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)


    enddo
  end subroutine doHydrodynamics2d


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
             thisOctal%phi_i(subcell) = thisOctal%tempStorage(subcell,8)
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
    real(double) :: machNumber, Pr, rhor, gamma

    gamma = 5.d0/3.d0
    machNumber = 2.d0
  
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

             if (.not.associated(thisOctal%tempStorage)) then
                allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:8))
                thisOctal%tempStorage = 0.d0
             endif

             select case(thisOctal%boundaryCondition(subcell))
                case(1) ! reflecting
                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocal(locator, bOctal, bSubcell)

                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)
                   Thisoctal%tempstorage(subcell,1) = bOctal%rho(bSubcell)
                   thisOctal%tempStorage(subcell,2) = bOctal%rhoE(bSubcell)

                   thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                   thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                   thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)

                   thisOctal%tempStorage(subcell,6) = bOctal%energy(bSubcell)
                   thisOctal%tempStorage(subcell,7) = bOctal%pressure_i(bSubcell)
                   thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bSubcell)

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

                case(2) ! periodic



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
                   thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bSubcell)

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

!          corner = .false.

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
                      if (.not.inOctal(grid%octreeRoot, locator)) then
                         write(*,*) "locator problem for period boundary"
                         write(*,*) "locator ",locator
                         write(*,*) "probe ", probe(iprobe)
                         write(*,*) "neighbour octal centre ",subcellCentre(neighbourOctal, neighbourSubcell)
                         write(*,*) "rvec ",rvec
                         stop
                      endif
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
                     inherit=.false., interp=.false., amrHydroInterp = .true.)
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
    logical, intent(out) :: split
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



  recursive subroutine refineGridGeneric2(thisOctal, grid, converged, inheritval)
    use input_variables, only : maxDepthAMR, photoionization
    include 'mpif.h'
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    logical, optional :: inheritval
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
    type(VECTOR) :: dirVec(6), centre, locator
    logical :: split
    integer :: neighbourSubcell, nDir
    real(double) :: r, grad, maxGradient
    real(double), parameter :: limit = 0.1d0
!    real(double) :: cs, rhocs
    integer :: myRank, ierr
    logical :: refineOnMass, refineOnIonization, refineOnGradient
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi
    real(double) :: mach, speed, cs
    integer :: nd

    converged = .true.
    converged_tmp=.true.
    

    refineOnGradient = .not.photoionization
    refineOnMass = .false.
    refineOnIonization = photoionization


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call refineGridGeneric2(child, grid, converged, inheritval)
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


          if (refineOnGradient) then
          do i = 1, nDir
             maxGradient = 1.d-30
             locator = subcellCentre(thisOctal, subcell) + &
                  (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
             if (inOctal(grid%octreeRoot, locator)) then
                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                

                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, dirvec(i), q, rho, rhoe, &
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)

                cs = soundSpeed(neighbourOctal, neighbourSubcell)
                speed = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
                     thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)
                mach = speed/cs
                
                if (mach > 1.1d0) then
                   write(*,*) "Shock detected ", mach
                   split = .true.
                endif

                cs = soundSpeed(thisOctal, subcell)
                speed = sqrt(rhou**2 + rhov**2 + rhow**2) / rho
                mach = speed/cs
                
                if (mach > 1.1d0) then
                   write(*,*) "Shock detected ", mach
                   split = .true.
                endif
                   


                grad = abs((thisOctal%rho(subcell)-rho) / &
                     thisOctal%rho(subcell))
                maxGradient = max(grad, maxGradient)
                if (maxGradient > limit) then
                   split = .true.
                else
                   split = .false.
                endif
                
                
                if (split) then
                   if (thisOctal%nDepth < maxDepthAMR) then
                      call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                           inherit=.false., interp=.false., amrHydroInterp = .true.)
                      converged = .false.
                      exit
                   endif
                   
                   if (neighbourOctal%nDepth < maxDepthAMR) then
                      call addNewChild(neighbourOctal,neighboursubcell,grid,adjustGridInfo=.TRUE., &
                           inherit=.false., interp=.false., amrHydroInterp = .true.)
                      converged = .false.
                      exit
                   endif
                endif
             
          endif
       enddo
    endif
!

    if (converged.and.refineOnMass) then
       if (((thisOctal%rho(subcell)*cellVolume(thisOctal, subcell)*1.d30) > 1.d-5*mSol) &
            .and.(thisOctal%nDepth < maxDepthAMR))  then
          call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
               inherit=.false., interp=.false., amrHydroInterp = .true.)
          converged = .false.
          exit
       endif
    endif


    if (converged.and.refineOnIonization) then

       do i = 1, nDir
          locator = subcellCentre(thisOctal, subcell) + &
               (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
          if (inOctal(grid%octreeRoot, locator)) then
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

             if (octalOnThread(neighbourOctal, neighbourSubcell, myRank)) then
                if ((thisOctal%ionFrac(subcell,1) > 0.9d0).and.(neighbourOctal%ionFrac(neighbourSubcell,1) < 0.1d0) .and. &
                     (thisOctal%nDepth < maxDepthAMR) ) then
                   call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                        inherit=.false., interp=.false., amrHydroInterp = .true.)
                   converged = .false.
                   exit
                endif

                if ((thisOctal%ionFrac(subcell,1) < 0.4d0).and.(neighbourOctal%ionFrac(neighbourSubcell,1) > 0.6d0) .and. &
                     (neighbourOctal%nDepth < maxDepthAMR) ) then
                   call addNewChild(neighbourOctal,neighbourSubcell,grid,adjustGridInfo=.TRUE., &
                        inherit=.false., interp=.false., amrHydroInterp = .true.)
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
                  inherit=.false., interp=.false., amrHydroInterp = .true.)
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
                  inherit=.false., interp=.false., amrHydroInterp = .true.)
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
       
  recursive subroutine unrefineCells(thisOctal, grid, nUnrefine)
    use input_variables, only : minDepthAMR
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i
    logical :: unrefine
    integer :: nc
    integer, intent(inout) :: nUnrefine
    real(double) :: rhow(8), rhov(8), rhou(8), rho(8), rhoe(8) , fac, limit
    real(double) :: cs(8), mass
    real(double) :: rhocs, rhomean, rhoemean
    logical :: refinedLastTime, ghostCell

    limit  = 0.01d0

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
                call unrefineCells(child, grid, nUnrefine)
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

          cs(nc) = soundSpeed(thisOctal, subcell)
          if (thisOctal%ghostCell(subcell)) ghostCell=.true.
!          if (thisOctal%refinedLastTime(subcell)) refinedLastTime = .true.
       endif
    enddo


    unrefine = .false.

    if ((nc > 1)) then !.and.(.not.ghostCell)) then

       rhomean = SUM(rho(1:nc))/dble(nc)
       rhoemean = SUM(rhoe(1:nc))/dble(nc)
       rhocs = rhomean*(SUM(cs(1:nc))/dble(nc))
       unrefine = .true.
       
       fac = sigma(rho, nc) / rhoMean
       if (fac > limit) then
          unrefine = .false.
!          write(*,*) "rho",fac
       endif
      
       if (rhoeMean /= 0.d0) then
          fac = sigma(rhoe, nc) / rhoeMean
          if (fac > limit) then
             unrefine = .false.
             !          write(*,*) "rhoe",fac
          endif
       endif

       if (rhocs /= 0.d0) then
          fac = sigma(rhou, nc) / rhocs
          if (fac > limit) then
             unrefine = .false.
             !          write(*,*) "rhou",fac
          endif
       endif
!       
       if (rhocs /= 0.d0) then
          fac = sigma(rhov, nc) / rhocs
          if (fac > limit) then
             unrefine = .false.
             !          write(*,*) "rhov",fac
          endif
       endif
!
       if (rhocs /= 0.d0) then
          fac = sigma(rhow, nc) / rhocs
          if (fac > limit) then
             unrefine = .false.
             !          write(*,*) "rhov",fac
          endif
       endif
!       
    endif
!    if (mass  < 2.d-6*mSol) unrefine = .true.
       
    if (thisOctal%nDepth <= minDepthAMR) unrefine = .false.

    if ((thisOctal%nChildren == 0).and.unrefine) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
       nunrefine = nunrefine + 1
    endif
    
  end subroutine unrefineCells

  recursive subroutine unrefineCellsPhotoion(thisOctal, grid)
    use input_variables, only : minDepthAMR
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i
    logical :: unrefine
    integer :: nc
    real(double) :: ionfrac(8), limit
    real(double) :: mass
    logical :: refinedLastTime, ghostCell
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
                call unrefineCellsPhotoIon(child, grid)
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

  subroutine evenUpGridMPI(grid, inheritFlag, evenAcrossThreads, dumpfiles)

    type(GRIDTYPE) :: grid
    logical :: gridConverged, inheritFlag
    integer, optional :: dumpfiles
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
    character(len=30) :: vtkFilename

!    character(len=20) :: plotfile
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)
    ! first even up the local part of the grid


    globalConverged = .false.
    iter = 0

    if (PRESENT(dumpfiles)) then
       write(vtkfilename,'(a,i4.4,a,i4.4,a)') "start",dumpfiles,"_",iter,".vtk"
       call writeVtkFile(grid, vtkFilename)
    endif

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
          call setupGhostCells(grid%octreeRoot, grid)!, flag=.true.)
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

             endif
          enddo

          if (PRESENT(dumpfiles)) then
             write(vtkfilename,'(a,i4.4,a,i4.4,a)') "aftersplit",dumpfiles,"_",iter,".vtk"
             call writeVtkFile(grid, vtkFilename)
          endif

          do
             gridConverged = .true.
             call evenUpGrid(grid%octreeRoot, grid,  gridconverged, inherit=inheritFlag)
             if (.not.gridConverged) localChanged(myRank) = .true.
             call setupEdges(grid%octreeRoot, grid)
             call unsetGhosts(grid%octreeRoot)
             call setupGhostCells(grid%octreeRoot, grid) !, flag=.true.)
             if (gridConverged) exit
          end do

          if (PRESENT(dumpfiles)) then
             write(vtkfilename,'(a,i4.4,a,i4.4,a)') "aftereven",dumpfiles,"_",iter,".vtk"
             call writeVtkFile(grid, vtkFilename)
          endif

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
                              inherit=.false., interp=.false.,  amrHydroInterp = .true.)
                         converged = .false.
                         exit
                      endif
                      
                      if ((thisOctal%nDepth-neighbourOctal%nDepth) > 1) then
                         call addNewChild(neighbourOctal,neighboursubcell,grid,adjustGridInfo=.TRUE., &
                              inherit=.false., interp=.false., amrHydroInterp = .true.)
                         converged = .false.
                         exit
                      endif
                      
!                      if (thisOctal%edgeCell(subcell).and.(.not.neighbourOctal%edgeCell(neighboursubcell))  &
!                           .and.(thisOctal%nDepth < neighbourOctal%nDepth)) then
!                         call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
!                              inherit=inherit, interp=.false.)
!                         converged = .false.
!                         exit
!                      endif
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
            inherit=.false., interp=.false., amrHydroInterp = .true.)
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
    type(VECTOR) :: dirVec(6), centre, octVec
    type(VECTOR), intent(out) :: loc(:)
    integer, intent(out) :: thread(:), depth(:)
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    integer :: myRank, ierr
    integer, intent(out) :: nLocs
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

!             call imposeBoundary(grid%octreeRoot)
!             call periodBoundary(grid)
!             call transferTempStorage(grid%octreeRoot)
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

!       call imposeBoundary(grid%octreeRoot)
!       call periodBoundary(grid)
!       call transferTempStorage(grid%octreeRoot)

       call MPI_ALLREDUCE(fracChange, tempFracChange, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR, ierr)
       
       fracChange = tempFracChange
       !       write(plotfile,'(a,i4.4,a)') "grav",it,".png/png"
!e           if (myrankglobal == 1)   write(*,*) it,MAXVAL(fracChange(1:nHydroThreads))
    enddo
    if (myRankGlobal == 1) write(*,*) "Gravity solver completed after: ",it, " iterations"


!    if (myrankglobal == 1) call tune(6,"Complete self gravity")


  end subroutine selfGrav

  real(double) function getPressure(thisOctal, subcell)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: eKinetic, eThermal, K, u2, eTot
    real(double), parameter :: gamma2 = 1.4d0, rhoCrit = 1.d-14


    select case(thisOctal%iEquationOfState(subcell))
       case(0) ! adiabatic
          if (thisOctal%threed) then
             u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
          else
             u2 = (thisOctal%rhou(subcell)**2 +  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
          endif
          eKinetic = u2 / 2.d0
          eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
          eThermal = eTot - eKinetic


          getPressure =  (thisOctal%gamma(subcell) - 1.d0) * thisOctal%rho(subcell) * eThermal
!          write(*,*) "gamma, ethermal",thisOctal%gamma(subcell), ethermal
!          write(*,*) "rhou,rhow ", thisOCtal%rhou(subcell), thisOctal%rhow(subcell)
!          write(*,*) "etot, ekinetic, eThermal ",etot, ekinetic, ethermal
       case(1) ! isothermal
          eThermal = thisOctal%rhoe(subcell) / thisOctal%rho(subcell)
          getPressure =  (thisOctal%gamma(subcell) - 1.d0) * thisOctal%rho(subcell) * eThermal
!          getPressure = (thisOctal%rho(subcell)/(2.d0*mHydrogen))*kerg*thisOctal%temperature(subcell)

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
          getPressure = ((4.1317d9*1.d-20)/(1.d-20**2)) * thisOctal%rho(subcell)**2
          

       case DEFAULT
          write(*,*) "Equation of state not recognised: ", thisOctal%iEquationOfState(subcell)
          stop
       end select
!       if (myrankglobal == 1) write(*,*) "p, gamma, ios ", &
!       getPressure, thisOctal%gamma(subcell), thisOctal%iEquationOfState(subcell)         
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
    
  subroutine testCell(grid, rVec)
    type(VECTOR) :: rVec
    type(OCTAL), pointer :: thisOctal
    type(GRIDTYPE) :: grid
    integer :: subcell

    call findSubcellTD(rVec, grid%octreeRoot, thisOctal, subcell)
    
    if (octalonthread(thisOctal, subcell, myRankGlobal)) then
       write(*,*) "-----------------------------------------------"
       write(*,*) "Vector: ",rVec
       write(*,*) "Q_i-2: ", thisOctal%q_i_minus_2(subcell)
       write(*,*) "Q_i-1: ", thisOctal%q_i_minus_1(subcell)
       write(*,*) "Q_i: ", thisOctal%q_i(subcell)
       write(*,*) "Q_i+1: ", thisOctal%q_i_plus_1(subcell)
       write(*,*) "flux_i: ", thisOctal%flux_i(subcell)
       write(*,*) "flux_i+1: ", thisOctal%flux_i_plus_1(subcell)
       write(*,*) "u_i: ", thisOctal%u_interface(subcell)
       write(*,*) "limiter: ", thisOctal%phiLimit(subcell)
       write(*,*) "uvel: ", thisOctal%rhou(subcell)/thisOctal%rho(subcell)
       write(*,*) "wvel: ", thisOctal%rhow(subcell)/thisOctal%rho(subcell)
    endif
  end subroutine testCell

#endif
    
end module hydrodynamics_mod
