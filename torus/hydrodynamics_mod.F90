#ifdef HYDRO
! hydrodynamics module added by tjh on 18th june 2007


module hydrodynamics_mod
#ifdef MPI

  use inputs_mod
  use dimensionality_mod
  use kind_mod
  use constants_mod
  use amr_mod
  use grid_mod
  use source_mod
  use timing
  use mpi_amr_mod
  use gridio_mod
  use vtk_mod
  use nbody_mod
  implicit none

  public

  type(OCTALWRAPPER), allocatable :: globalChildlessOctalArray(:)
  integer :: nGlobalChildlessOctals


contains

  subroutine dohydrodynamics(grid)
    type(GRIDTYPE) :: grid

    if (grid%octreeRoot%twoD) then
       call doHydrodynamics2d(grid)
    else if (grid%octreeRoot%oneD) then
       call doHydrodynamics1d(grid)
    else if (grid%octreeRoot%threeD) then
       call doHydrodynamics3d(grid)
    endif
    
  end subroutine dohydrodynamics

  recursive subroutine fluxlimiter(thisoctal)
    use inputs_mod, only : limiterType
    use mpi
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: a, b, dq, dx
    character(len=80) :: message
 
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call fluxlimiter(child)
                exit
             end if
          end do
       else

!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

      !    if (.not.thisoctal%ghostcell(subcell)) then
             dq = thisoctal%q_i(subcell) - thisoctal%q_i_minus_1(subcell)
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
                a = min(1.d0, 2.d0*thisoctal%rlimit(subcell))
                b = min(2.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit(subcell) = max(0.d0, a, b)
                
             case("modified_superbee")
                if (dq /= 0.d0) then
                   if (thisoctal%u_interface(subcell) .ge. 0.d0) then
                      thisoctal%rlimit(subcell) = (thisoctal%q_i_minus_1(subcell) - thisoctal%q_i_minus_2(subcell)) / dq
                      dx = 0.5d0*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))

                      thisoctal%rlimit(subcell)=thisoctal%rlimit(subcell)*dx/(thisoctal%x_i_minus_1(subcell) - &
                           thisoctal%x_i_minus_2(subcell))

                   else
                      thisoctal%rlimit(subcell) = (thisoctal%q_i_plus_1(subcell) - thisoctal%q_i(subcell)) / dq
                      dx = 0.5d0*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
                
                      thisoctal%rlimit(subcell)=thisoctal%rlimit(subcell)*dx/(thisoctal%x_i_plus_1(subcell) - &
                           thisoctal%x_i(subcell))
                   endif
                else
                   thisoctal%rlimit(subcell) = 1.d0
                endif

                a = min(1.d0, 2.d0*thisoctal%rlimit(subcell))
                b = min(2.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit(subcell) = max(0.d0, a, b)

             case("minmod")
                
                a = min(1.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit(subcell) = max(0.d0, a)
                
             case("MC")
                a = min((1.d0+thisoctal%rlimit(subcell))/2.d0, 2.d0, (2.d0*thisoctal%rlimit(subcell)))
                thisoctal%philimit(subcell) = max(0.d0, a)

             case("fromm")
                thisoctal%philimit(subcell) = 0.5d0*(1.d0 + thisoctal%rlimit(subcell))
                

             case("vanleer")
                !write(*,*) "vanleer"
                thisoctal%philimit(subcell) = (thisoctal%rlimit(subcell) + &
                abs(thisoctal%rlimit(subcell))) / (1 + abs(thisoctal%rlimit(subcell)))
               
             case("donorcell")
                thisOctal%philimit(subcell) = 0.d0


             case default
                write(message,'(a,a)') "flux limiter not recognised: ", trim(limitertype)
                call writefatal(message)
                stop
             end select
!             if (thisoctal%ghostcell(subcell)) thisoctal%philimit(subcell) = 1.d0
!                        thisoctal%philimit(subcell) = 1.d0
         ! endif
          endif

       enddo
  end subroutine fluxlimiter


  recursive subroutine updatedensitytree(thisoctal)
    use mpi
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
          parentoctal%phi_gas(m) = sum(testoctal%phi_gas(1:n))/dble(n)
          parentoctal%edgecell(m) = any(testoctal%edgecell(1:n))
          testoctal => parentoctal
       enddo
    endif
  end subroutine updatedensitytree

  recursive subroutine updatephitree(thisoctal, ndepth)
    use mpi
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

       thisoctal%phi_gas = thisoctal%parent%phi_gas(thisoctal%parentsubcell)
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
    use mpi
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
             !print *, "thisOctal%flux(subcell) = ", thisOctal%flux_i(subcell)
        
             if (thisoctal%x_i(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) "problem with the x_i values"
                stop
             endif


             dx = thisOctal%subcellSize*griddistancescale

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
    use mpi
    integer :: myrank, ierr
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx, df

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

          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle

          if (.not.thisoctal%ghostcell(subcell)) then
          
!             dx = 0.5*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
!             dx = (thisoctal%x_i_plus_1(subcell) - thisoctal%x_i(subcell))
             dx = thisOctal%subcellSize*griddistancescale

             df = (thisoctal%flux_i_plus_1(subcell) - thisoctal%flux_i(subcell))

             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * (df/dx)


          endif
       endif
    enddo
  end subroutine updatecellq

  recursive subroutine synchronizefluxes(thisoctal, dt, idepth)
    use mpi
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


!THaw boundary partners should be directly opposite their ghosts
  recursive subroutine checkBoundaryPartners(thisOctal, grid)
    type(gridtype) :: grid
    type(octal), pointer :: thisOctal
    type(octal), pointer :: child
    integer :: subcell, i
    type(vector) :: rVec, bVec, direction

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child                                                                                                                                                                              
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call checkBoundaryPartners(child, grid)
                exit
             end if
          end do
       else
          if(thisOctal%ghostcell(subcell)) then
             rVec = subcellCentre(thisOctal, subcell)
             bVec = thisOctal%boundaryPartner(subcell)

             direction = bVec - rVec
             
!             if(abs(direction%x) > abs(direction%y)) then
!                direction%y = 0.d0
!                if(abs(direction%x) > abs(direction%z)) then !We are moving in the 悉 direction                                                                                                         
!                   direction%z = 0.d0
!                else                                         !Moving in the 您 direction                                                                                                                
!                   direction%x = 0.d0
!                end if
!             else
!                direction%x = 0.d0
!                if(abs(direction%y) > abs(direction%z)) then !We are moving in the 悠 direction                                                                                                         
!                   direction%z = 0.d0
!                else                                         !Moving in the 您 direction                                                                                                            $
!                   direction%y = 0.d0
!                end if
!             end if

             !Make it a unit vector                                                                                                                                                                     
             if(direction%x > 0.d0) direction%x = direction%x / direction%x
             if(direction%x < 0.d0) direction%x = -direction%x / direction%x
             if(direction%y > 0.d0) direction%y = direction%y / direction%y
             if(direction%y < 0.d0) direction%y = -direction%y / direction%y
             if(direction%z > 0.d0) direction%z = direction%z / direction%z
             if(direction%z < 0.d0) direction%z = -direction%z / direction%z

             
             if((abs(direction%x) + abs(direction%y) + abs(direction%z)) /= 1.d0) then
                print *, "boundary partner is at a diagonal!"
                print *, "direction ", direction
                print *, "rVec ", rVec
                print *, "bVec ", bVec
                stop
             end if

          end if
          
       end if
    end do

  end subroutine checkBoundaryPartners


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
          thisoctal%x_i(subcell) = returnCodeUnitLength(griddistancescale*(subcellcentre(thisoctal, subcell) .dot. direction))
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
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi, phigas
    real(double) :: xnext
    integer :: subcell, i, neighboursubcell
    integer :: nd
    type(vector) :: direction, locator, reversedirection


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
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          thisoctal%x_i_minus_1(subcell) = 0.d0
          thisoctal%x_i_plus_1(subcell) = 0.d0
          if (.not.thisoctal%edgecell(subcell)) then

             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)

             thisoctal%x_i_plus_1(subcell) = x
             thisoctal%q_i_plus_1(subcell) = q

             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%x_i_minus_1(subcell) = x
!             thisoctal%x_i_minus_2(subcell) = xnext
             thisoctal%q_i_minus_1(subcell) = q
             thisoctal%q_i_minus_2(subcell) = qnext
!
!             if (thisoctal%ndepth == nd) then
!                thisoctal%q_i_minus_1(subcell) = q
!
!             else if (thisoctal%ndepth > nd) then ! fine to coarse
!                thisoctal%q_i_minus_1(subcell) = (q + thisOctal%q_i!(subcell))*2.d0/3.d0
!                !
!             els!e!!
!
!                thisoctal%u_interface(subcell) = &
!                     (2.d0/3.d0)*(((thisoctal%rhou(subcell)/thisoctal%rho(subcell))/2.d0) + &
!                     (rhou_i_minus_1/rho_i_minus_1))
!             endif



!             write(*,*) "q: ", thisoctal%q_i_minus_2(subcell), thisoctal%q_i_minus_1(subcell), thisoctal%q_i(subcell), &
!                  thisoctal%q_i_plus_1(subcell)
             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrankGlobal," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell)
             endif

          endif
       endif
    enddo
  end subroutine setupqx

  recursive subroutine setuprhophi(thisoctal, grid, direction)
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi, phigas 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator, reversedirection
    integer :: nd
    real(double) :: xnext
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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)

             thisoctal%rho_i_plus_1(subcell) = rho
             thisoctal%phi_i_plus_1(subcell) = phi

             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%rho_i_minus_1(subcell) = rho
             thisoctal%phi_i_minus_1(subcell) = phi



          endif
       endif
    enddo
  end subroutine setuprhophi

  recursive subroutine setupui(thisoctal, grid, direction)
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas
    integer :: subcell, i, neighboursubcell, nd
    type(vector) :: direction, locator!, rotator
    real(double) :: rhou_i_minus_1, rho_i_minus_1
    real(double) :: xnext, weight

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    weight = 0.5d0
!    nOdd = 0
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

          !thaw - edgecell is correct.
          if (.not.thisoctal%edgecell(subcell)) then !xxx changed fromghostcell
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
            
             rho_i_minus_1 = rho
             rhou_i_minus_1 = rhou
          
             thisoctal%u_interface(subcell) = &
                  weight*thisoctal%rhou(subcell)/thisoctal%rho(subcell) + &
                  (1.d0-weight)*rhou_i_minus_1/rho_i_minus_1
     
          endif
       endif
    enddo

  end subroutine setupui

  recursive subroutine setupvi(thisoctal, grid, direction)
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, weight
    integer :: nd
    real(double) :: xnext
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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)

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
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, weight
    integer :: nd
    real(double) :: xnext
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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)

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
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi, phigas
    integer :: nd
    real(double) :: xnext, fac

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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)

   
             !thisoctal%flux_i_plus_1(subcell) = flux

             !For fine cells constructing a +1 flux from a coarser cell, it is worth checking the flux gradient    
             if(thisOctal%nDepth > nd .and. fluxinterp) then
                if(.not. thisOctal%oneD) then
                   neighbourOctal%flux_i(neighboursubcell) = flux
!                   probe = subcellCentre(neighbourOctal, neighbourSubcell)
                   call NormalFluxGradient(thisOctal, subcell, neighbourOctal, neighbourSubcell, grid, direction, fac)
!                   if(fac /= 0.d0) print *, "fac", fac
                   else
                   fac = 0.d0
                   end if
             else
                fac = 0.d0
             end if

             thisoctal%flux_i_plus_1(subcell) = flux + fac

             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%flux_i_minus_1(subcell) = flux

! new stuff added by tjh

!             if (thisoctal%ndepth >= nd) then ! this is a coarse-to-fine cell boundary
                thisoctal%flux_i_minus_1(subcell) = flux
!             else
!                ! now we need to do the fine-to-coarse flux
!                if (thisoctal%u_i_minus_1(subcell) .ge. 0.d0) then ! flow is out of this cell into next
!                   thisoctal%flux_i_minus_1(subcell) = thisoctal%q_i(subcell) * thisoctal%u_i_minus_1(subcell)
!                else
!                   thisoctal%flux_i_minus_1(subcell) = flux ! flow is from neighbour into this one
!                endif
!             endif

          endif
       endif
    enddo
  end subroutine setupflux


!THaw - Flux interpolation routine                                                                                                                                                                                                       
  subroutine normalFluxGradient(thisOctal, subcell, neighbourOctal, neighbourSubcell, grid, direction, fac)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer :: thisoctal
    type(octal), pointer :: neighbourOctal
    type(octal), pointer :: communityOctal
    integer :: subcell, neighbourSubcell, communitySubcell
    type(vector) :: direction
    type(vector) :: rVec
    type(vector), allocatable :: community(:)
    type(vector) :: locator
    type(vector) :: nVec
    real(double), intent(out) :: fac
    real(double), allocatable ::  xpos(:), f(:)
    integer :: myRank, ierr
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi, phigas
    integer :: nd, iTot, i
    real(double) :: xnext, m, dx, df

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    nVec = subcellCentre(neighbourOctal, neighbourSubcell)

    fac = 0.d0
    if(neighbourOctal%twoD) then
       iTot = 2
       allocate(community(iTot))
       allocate(xpos(iTot))
       allocate(f(iTot))
       !Get the direction of the neighbour's community cells                                                                                                                                                                             
       if(abs(direction%x) == 1.d0) then !Advecting in 悉-direction                                                                                                                                                                      
          community(1) = VECTOR(0.d0, 0.d0,1.d0)
          community(2) = VECTOR(0.d0, 0.d0,-1.d0)
       else if(abs(direction%z) == 1.d0) then !Advecting in 您 direction                                                                                                                                                                 
          community(1) = VECTOR(1.d0, 0.d0,0.d0)
          community(2) = VECTOR(-1.d0, 0.d0,0.d0)
       else
          print *, "Advecting in unknown direction!"
          stop
       end if

       do i = 1, iToT
          locator = subcellcentre(neighbouroctal, neighboursubcell) + community(i) * &
               (neighbouroctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)

          if(inOctal(grid%octreeRoot, locator)) then
             communityoctal => neighbouroctal
             call findsubcelllocal(locator, communityoctal, communitysubcell)
             call getneighbourvalues(grid, neighbouroctal, neighboursubcell, communityoctal, communitysubcell, direction, q, &
                  rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)

             rVec = subcellcentre(communityoctal, communitysubcell)
             f(i) = flux
             if(abs(direction%x) == 1.d0) then
                xpos(i) = rVec%z
             else if(abs(direction%z) == 1.d0) then
                xpos(i) = rVec%x
             else
                print *, "unrecognized direction!"
                stop
             end if

          else
             !At the edge just use the gradient over half the distance                                                                                                                                                                   
             f(i)  = neighbourOctal%flux_i(subcell)
             if(abs(direction%x) == 1.d0) then
                xpos(i) = nVec%z
             else if(abs(direction%z) == 1.d0) then
                xpos(i) = nVec%x
             else
                print *, "unrecognized direction!"
                stop
             end if


          end if
       end do

       !Flux gradient                                                                                                                                                                                                                    
       m = (f(1) - f(2))/(xpos(1) - xpos(2))

       dx = thisOctal%subcellSize*griddistancescale


       !Flux variation for coarse cell centre to fine cell centre, perpendicular to advection direction                                                                                                                                  
       df = abs(m*dx)/2.d0

       rVec = subcellCentre(thisOctal, subcell)!                                                                                                                                                                                         
       if(m > 0.d0) then !Positive gradient                                                                                                                                                                                              

          if(abs(direction%x) == 1.d0) then !悉 advection                                                                                                                                                                                
             if(rVec%z > nVec%z) then         !Upper cell                                                                                                                                                                                
                   fac = df
             else                             !lower cell                                                                                                                                                                                
                   fac = -df
             end if
          else                                !您 advection                                                                                                                                                                              

             if(rVec%x > nVec%x) then         !Right cell                                                                                                                                                                                
                fac = df
             else                             !Left cell                                                                                                                                                                                 
                   fac = -df
             end if
          end if
       else if (m < 0.d0) then   !Negative gradient                                                                                                                                                               
          if(abs(direction%x) == 1.d0) then !悉 advection                                                                                                                                                                                
             if(rVec%z > nVec%z) then         !Upper cell                                                                                                                                                                                
                fac = -df
             else                             !lower cell                                                                                                                                                                                
                fac = df
             end if
          else                                !您 advection                                                                                                                                                                              

             if(rVec%x > nVec%x) then         !Right cell                                                                                                                                                                                
                fac = -df
             else                             !Left cell                                                                                                                                                                                 
                fac = df
             end if

          end if

       else                      !No gradient present                                                                                                                                                                                    
          fac = 0.d0
       end if

    else !3D case                                                                                                                                                                                                                        

       !Only going to do the 2D properly for now                                                                                                                                                                                         

       fac = 0.d0

    end if

    deallocate(community)
    deallocate(xpos)
    deallocate(f)

  end subroutine normalFluxGradient

  recursive subroutine setuppressure(thisoctal, grid, direction)
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    integer :: nd
    real(double) :: xnext  

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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%pressure_i_plus_1(subcell) = pressure

             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%pressure_i_minus_1(subcell) = pressure

          endif
       endif
    enddo
  end subroutine setuppressure

  recursive subroutine setupupm(thisoctal, grid, direction)
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas
    integer :: nd
    real(double) :: xnext

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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%u_i_plus_1(subcell) = rhou/rho
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%u_i_minus_1(subcell) = rhou/rho

          endif
       endif
    enddo
  end subroutine setupupm

  recursive subroutine setupvpm(thisoctal, grid, direction)
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas
    integer :: nd
    real(double) :: xnext

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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%u_i_plus_1(subcell) = rhov/rho
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%u_i_minus_1(subcell) = rhov/rho

          endif
       endif
    enddo
  end subroutine setupvpm

  recursive subroutine setupwpm(thisoctal, grid, direction)
    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas
    integer :: nd
    real(double) :: xnext

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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%u_i_plus_1(subcell) = rhow/rho
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
             thisoctal%u_i_minus_1(subcell) = rhow/rho

          endif
       endif
    enddo
  end subroutine setupwpm


  !Thaw - Rhie-Chow interpolation 
  recursive subroutine rhiechowui(thisoctal, grid, direction, dt)
     !Subroutine to modify the cell interface velocities according to Rhie-Chow
     !interpolation

    use mpi
    integer :: myrank, ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child
    integer :: subcell, i
    type(vector) :: direction
    real(double) :: dt, dx
    
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call rhiechowui(child, grid, direction, dt)
                exit
             end if
          end do
       else

      
!          if (thisoctal%mpithread(subcell) /= myrank) cycle
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then
             !if (.not.thisoctal%edgecell(subcell)) then
             dx = returnCodeUnitLength(thisoctal%subcellsize*gridDistanceScale)
       !      dx = (thisoctal%x_i(subcell) - thisoctal%x_i_minus_1(subcell))
             thisoctal%u_interface(subcell) = thisoctal%u_interface(subcell) - dt * &
                  ((thisoctal%pressure_i_plus_1(subcell) + thisoctal%pressure_i(subcell))/2.d0 -  &
                  (thisoctal%pressure_i(subcell) + thisoctal%pressure_i_minus_1(subcell))/2.d0)/dx
             ! write(*,*), "ui", thisoctal%u_interface(subcell)
          endif
       endif
    enddo

  end subroutine rhiechowui

  recursive subroutine computepressureGeneral(grid, thisOctal, withViscosity)
     use inputs_mod, only : etaViscosity, useViscosity
     use mpi
     integer :: myrank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child
    integer :: subcell, i
    real(double) :: biggamma,eta
    logical :: withViscosity
!    logical :: useviscosity                                                                                                                                                                                             

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    eta = etaViscosity

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child                                                                                                                                                                                               
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call computepressureGeneral(grid, child, withViscosity)
                exit
             end if
          end do
       else
!          if (thisoctal%mpithread(subcell) /= myrank) cycle                                                                                                                                                             
          if (.not.octalonthread(thisoctal, subcell, myrank)) cycle


          thisoctal%pressure_i(subcell) = getpressure(thisoctal, subcell)
  !        if(thisOctal%pressure_i(subcell) == 0.d0) then                                                                                                                                                                
  !           print *, "ERROR, pressure is zero"                                                                                                                                                                         
  !           stop                                                                                                                                                                                                       
  !        end if                                                                                                                                                                                                        

          biggamma = 0.d0


          if(withViscosity) then
             if (.not.thisoctal%edgecell(subcell)) then
                useviscosity = .false.
                if (thisoctal%u_i_plus_1(subcell) .le. thisoctal%u_i_minus_1(subcell)) useviscosity = .true.
                if (useviscosity) then
                   !     if ((thisoctal%u_i_plus_1(subcell) < 0.d0).and.(thisoctal%u_i_plus_1(subcell) .ge. thisoctal%u_i_minus_1(subcell))) &                                                                              
                   !     useviscosity = .true.                                                                                                                                                                              
                   
                   biggamma = 0.25d0 * eta**2 * (thisoctal%u_i_plus_1(subcell) - thisoctal%u_i_minus_1(subcell))**2 &
                        * thisoctal%rho(subcell)
                else
                   biggamma = 0.d0
                endif
             endif
             
             thisoctal%pressure_i(subcell) = thisoctal%pressure_i(subcell) + biggamma
             
          end if

          if (isnan(thisoctal%pressure_i(subcell))) then
             write(*,*) "pressureu has nan"
             write(*,*) "velocity: ",thisoctal%rhou(subcell),thisoctal%rhov(subcell), thisoctal%rhow(subcell)
             write(*,*) "rho: ", thisoctal%rho(subcell)
             write(*,*) "cen ",subcellcentre(thisoctal, subcell)
             stop
          endif

       endif

    enddo

  end subroutine computepressureGeneral

  recursive subroutine pressureforceu(thisoctal, dt)
    use mpi
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

!             dx = (thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
             dx = thisoctal%subcellsize * griddistancescale


!             thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - (dt/2.d0) * &!!!!!!!!!!!!!!!!!!!!!!!
!                  (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx

  
             !Thaw - Rhie-Chow interpolation
!             !print *, "rhou", thisoctal%rhou(subcell)
             if(rhieChow) then
                thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - (dt/2.d0) * &
                     ((thisoctal%pressure_i_plus_1(subcell) + thisoctal%pressure_i(subcell))/2.d0 -  & 
                    (thisoctal%pressure_i(subcell) + thisoctal%pressure_i_minus_1(subcell))/2.d0)/dx
                !print *, "rhou", thisoctal%rhou(subcell)

!                write(*,*) "rhou pressure ", thisOctal%rhou(subcell),  (dt) * &
!                     ((thisoctal%pressure_i_plus_1(subcell) + thisoctal%pressure_i(subcell))/2.d0 -  & 
!                    (thisoctal%pressure_i(subcell) + thisoctal%pressure_i_minus_1(subcell))/2.d0)/dx, &
!                    thisOctal%pressure_i(subcell)
!

             else 
                thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - (dt/2.d0) * &!!!!!!!!!!!!!!!!!!!!!!!
                     (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx


             end if

             thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - (dt/2.d0) * & !gravity
                  thisoctal%rho(subcell) *(thisoctal%phi_i_plus_1(subcell) - &
                  thisoctal%phi_i_minus_1(subcell)) / dx

             !     write(*,*) "rhou grav ",thisOctal%rhou(subcell),(dt/2.d0) * & !gravity
!                  thisoctal%rho(subcell) *(thisoctal%phi_i_plus_1(subcell) - &
!                  thisoctal%phi_i_minus_1(subcell)) / dx 


             if (isnan(thisoctal%rhou(subcell))) then
                write(*,*) "rhou ",thisoctal%rhou(subcell)
                write(*,*) "dt ",dt
                write(*,*) "phi ",thisoctal%phi_i_plus_1(subcell), thisoctal%phi_i_minus_1(subcell)
                write(*,*) "x ",thisoctal%x_i_plus_1(subcell), thisoctal%x_i_minus_1(subcell)
             endif

!             if (thisoctal%rhou(subcell)/thisoctal%rho(subcell) > 1.d30) then
!                write(*,*) "u ",thisoctal%rhou(subcell)/thisoctal%rho(subcell)/1.d5
!                write(*,*) "rho ", thisoctal%rho(subcell)
!                write(*,*) "rho i+1", thisoctal%rho_i_plus_1(subcell)
!                write(*,*) "rho i-1", thisoctal%rho_i_minus_1(subcell)
!                write(*,*) "pressure i + 1 " ,thisoctal%pressure_i_plus_1(subcell)
!                write(*,*) "pressure i - 1 " ,thisoctal%pressure_i_minus_1(subcell)
!                write(*,*) "x i + 1 ", thisoctal%x_i_plus_1(subcell)
!                write(*,*) "x i - 1 ", thisoctal%x_i_minus_1(subcell)
!             endif

                
!             write(*,*) thisoctal%rhou(subcell), thisoctal%pressure_i_plus_1(subcell), &
!                  thisoctal%pressure_i_minus_1(subcell), thisoctal%x_i_plus_1(subcell), &
!                  thisoctal%x_i_minus_1(subcell)

!             if (iequationofstate == 0) then


             !Thaw - apply similar rhie-chow interpolation here?
             if (thisoctal%iequationofstate(subcell) /= 1) then             
                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - (dt/2.d0) * &!!!!!!!!!!!
                     (thisoctal%pressure_i_plus_1(subcell) * thisoctal%u_i_plus_1(subcell) - &
                     thisoctal%pressure_i_minus_1(subcell) * thisoctal%u_i_minus_1(subcell)) / dx
                     
                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - (dt/2.d0)* & !gravity
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
    use mpi
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

!             dx = 0.5d0*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
             dx = thisoctal%subcellsize * griddistancescale


         !    thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - (dt/2.d0) * &
         !         (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx

            if(rhieChow) then
               thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - (dt/2.d0) * &
                    ((thisoctal%pressure_i_plus_1(subcell) + thisoctal%pressure_i(subcell))/2.d0 -  &
                    (thisoctal%pressure_i(subcell) + thisoctal%pressure_i_minus_1(subcell))/2.d0)/dx
                    !print *, "rhou", thisoctal%rhou(subcell)

!                write(*,*) "rhov pressure ", thisOctal%rhov(subcell),  (dt) * &
!                     ((thisoctal%pressure_i_plus_1(subcell) + thisoctal%pressure_i(subcell))/2.d0 -  & 
!                    (thisoctal%pressure_i(subcell) + thisoctal%pressure_i_minus_1(subcell))/2.d0)/dx, &
!                    thisOctal%pressure_i(subcell)


            else
                    thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - (dt/2.d0) * &
                    (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx
            end if
             thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - (dt/2.d0) * & !gravity
                  thisoctal%rho(subcell) *(thisoctal%phi_i_plus_1(subcell) - &
                  thisoctal%phi_i_minus_1(subcell)) / dx

                    
             if (isnan(thisoctal%rhov(subcell))) then
                write(*,*) "rhov ",thisoctal%rhov(subcell)
                write(*,*) "dt ",dt
                write(*,*) "phi ",thisoctal%phi_i_plus_1(subcell), thisoctal%phi_i_minus_1(subcell)
                write(*,*) "x ",thisoctal%x_i_plus_1(subcell), thisoctal%x_i_minus_1(subcell)
             endif

!             if (thisoctal%rhov(subcell)/thisoctal%rho(subcell) > 1.d10) then
!                write(*,*) "u ",thisoctal%rhov(subcell)/thisoctal%rho(subcell)/1.d5
!                write(*,*) "rho ", thisoctal%rho(subcell)
!                write(*,*) "rho i+1", thisoctal%rho_i_plus_1(subcell)
!                write(*,*) "rho i-1", thisoctal%rho_i_minus_1(subcell)
!                write(*,*) "pressure i + 1 " ,thisoctal%pressure_i_plus_1(subcell)
!                write(*,*) "pressure i - 1 " ,thisoctal%pressure_i_minus_1(subcell)
!                write(*,*) "x i + 1 ", thisoctal%x_i_plus_1(subcell)
!                write(*,*) "x i - 1 ", thisoctal%x_i_minus_1(subcell)
!             endif

                
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
    use mpi
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


!             dx = 0.5d0 * (thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
             dx = thisoctal%subcellsize * griddistancescale

            ! thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - (dt/2.d0) * &
            !      (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx


            if(rhieChow) then
               thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - (dt/2.d0) * &
               ((thisoctal%pressure_i_plus_1(subcell) + thisoctal%pressure_i(subcell))/2.d0 -  &
               (thisoctal%pressure_i(subcell) + thisoctal%pressure_i_minus_1(subcell))/2.d0)/dx
               !print *, "rhou", thisoctal%rhou(subcell)

!                write(*,*) "rhow pressure ", thisOctal%rhow(subcell),  (dt) * &
!                     ((thisoctal%pressure_i_plus_1(subcell) + thisoctal%pressure_i(subcell))/2.d0 -  & 
!                    (thisoctal%pressure_i(subcell) + thisoctal%pressure_i_minus_1(subcell))/2.d0)/dx, &
!                    thisOctal%pressure_i(subcell)

            else
                thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - (dt/2.d0) * &
                  (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / dx

            end if

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

  recursive subroutine copyIonfractoq(thisoctal, iion)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iion
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyIonFractoq(child, iion)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%ionFrac(subcell, iion)
        
       endif
    enddo
  end subroutine copyIonfractoq

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
             write(*,*) "rhou, rhov, rhow ", thisOctal%rhou(subcell),thisOctal%rhov(subcell), thisOctal%rhow(subcell)
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

  recursive subroutine copyqtoIonfrac(thisoctal, direction, iion)
    type(vector) :: direction
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iion
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtoIonfrac(child, direction, iion)
                exit
             end if
          end do
       else
  
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle

          thisoctal%ionFrac(subcell,iion) = thisoctal%q_i(subcell)
          if(thisoctal%ionFrac(subcell,iion) > 1.d0) then
             thisoctal%ionFrac(subcell,iion) = 1.d0
          end if
       endif
    enddo
  end subroutine copyqtoIonfrac

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

  subroutine advectIonFrac(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    integer :: i

    do i = 1, grid%nion
       call copyIonfractoq(grid%octreeroot, i)
       call advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
       call copyqtoIonfrac(grid%octreeroot, direction, i)
    enddo

  end subroutine advectIonFrac

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
    call fluxlimiter(grid%octreeroot)
    call constructflux(grid%octreeroot, dt)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupflux(grid%octreeroot, grid, direction)
    call updatecellq(grid%octreeroot, dt)
    
  end subroutine advectq


recursive subroutine sumFluxes(thisOctal, dt, totalFlux)
  use mpi
  integer :: subcell, i
  type(octal), pointer :: thisOctal
  type(octal), pointer :: child
  real(double) :: totalFlux, dt
  integer :: myRank, ierr

  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

! totalFlux = 0.d0

  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call sumFluxes(child, dt, totalFlux)
              exit
           end if
        end do
     else
        if (.not.octalonthread(thisoctal, subcell, myrank)) cycle
        if (.not.thisoctal%edgecell(subcell)) then
        totalFlux = totalFlux + thisOctal%flux_i(subcell)
      !  print *, "thisflux = ", thisOctal%flux_i(subcell)
        end if
     end if
  end do
 ! print *, "total flux = ", totalFlux
end subroutine sumFluxes


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
    if (rhieChow) then
     call computepressureGeneral(grid, grid%octreeroot, .false.) !Thaw Rhie-Chow might just need setuppressureu
     call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
     call setuppressure(grid%octreeroot, grid, direction) !Thaw
     call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
     call rhiechowui(grid%octreeroot, grid, direction, dt)
     call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    end if

    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)


    call setupui(grid%octreeroot, grid, direction)
    call setupupm(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)

    call computepressureGeneral(grid, grid%octreeroot, .true.) 
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    if(rhieChow) then
       call rhiechowui(grid%octreeroot, grid, direction, dt)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    end if
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
    use inputs_mod, only : nBodyPhysics
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
   if (selfGravity) then
!       call periodBoundary(grid, justGrav = .true.)
!       call transferTempStorage(grid%octreeRoot, justGrav = .true.)
    endif
    if (myrankglobal == 1) call tune(6,"Boundary conditions")
    
    if (myrankglobal == 1) call tune(6,"X-direction step")
    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)

    if(rhieChow) then
       call computepressureGeneral(grid, grid%octreeroot, .false.)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
       call setuppressure(grid%octreeroot, grid, direction) !Thaw
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
       call rhiechowui(grid%octreeroot, grid, direction, dt/2.d0)
    end if 

    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoV(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    if(photoionPhysics .and. hydrodynamics) then
       call advectIonFrac(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    end if
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)

    call imposeboundary(grid%octreeroot)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)

    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!    call computePressureU(grid, grid%octreeRoot, direction)
    call computepressureGeneral(grid, grid%octreeroot, .true.)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    if(rhieChow) then
       call rhiechowui(grid%octreeroot, grid, direction, dt/2.d0)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    end if
    call pressureForceU(grid%octreeRoot, dt/2.d0)

    if (myrankglobal == 1) call tune(6,"X-direction step")


    if (myrankglobal == 1) call tune(6,"Y-direction step")
    direction = VECTOR(0.d0, 1.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call setupVi(grid%octreeRoot, grid, direction)
    call setupVpm(grid%octreeRoot, grid, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    if(rhieChow) then
       call computepressureGeneral(grid, grid%octreeroot, .false.)
!       call computepressurev(grid, grid%octreeroot, direction) !Thaw
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=5)
       call setuppressure(grid%octreeroot, grid, direction) !Thaw
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=5)
       call rhiechowui(grid%octreeroot, grid, direction, dt)
    end if
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=5)   
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    if(photoionPhysics .and. hydrodynamics) then
       call advectIonFrac(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    end if
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)

    call imposeboundary(grid%octreeroot)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=5)


    call setupVi(grid%octreeRoot, grid, direction)
    call setupVpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    !call computePressureV(grid, grid%octreeRoot, direction)
    call computepressureGeneral(grid, grid%octreeroot, .true.)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=5)
    if(rhieChow) then
       call rhiechowui(grid%octreeroot, grid, direction, dt)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=5)
    end if  
   call pressureForceV(grid%octreeRoot, dt)
    if (myrankglobal == 1) call tune(6,"Y-direction step")

    if (myrankglobal == 1) call tune(6,"Z-direction step")
    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    if(rhieChow) then
!       call computepressurew(grid, grid%octreeroot, direction) !Thaw
       call computepressureGeneral(grid, grid%octreeroot, .false.)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)
       call setuppressure(grid%octreeroot, grid, direction) !Thaw
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)
       call rhiechowui(grid%octreeroot, grid, direction, dt)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)   
    end if
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    if(photoionPhysics .and. hydrodynamics) then
       call advectIonFrac(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    end if
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)

    call imposeboundary(grid%octreeroot)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)


    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    !call computePressureW(grid, grid%octreeRoot, direction)
    call computepressureGeneral(grid, grid%octreeroot, .true.)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    if(rhieChow) then
       call rhiechowui(grid%octreeroot, grid, direction, dt)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)
    end if
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
    if(rhieChow) then
       call computepressureGeneral(grid, grid%octreeroot, .false.)
       !call computepressureu(grid, grid%octreeroot, direction) !Thaw 
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
       call setuppressure(grid%octreeroot, grid, direction) !Thaw
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
       call rhiechowui(grid%octreeroot, grid, direction, dt/2.d0)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    end if
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoV(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    if(photoionPhysics .and. hydrodynamics) then
      call advectIonFrac(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    end if
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)

    call imposeboundary(grid%octreeroot)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)

    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!    call computePressureU(grid, grid%octreeRoot, direction)
    call computepressureGeneral(grid, grid%octreeroot, .false.)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupPressure(grid%octreeRoot, grid, direction)
    if(rhieChow) then
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
       call rhiechowui(grid%octreeroot, grid, direction, dt/2.d0)
    end if
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
   
    call pressureForceU(grid%octreeRoot, dt/2.d0)
    if (myrankglobal == 1) call tune(6,"X-direction step")



    if (myrankglobal == 1) call tune(6,"Boundary conditions")
    call periodBoundary(grid)
    call imposeBoundary(grid%octreeRoot)
    call transferTempStorage(grid%octreeRoot)

    !THAW
    if (selfGravity) then
!       call periodBoundary(grid, justGrav = .true.)
!       call transferTempStorage(grid%octreeRoot, justGrav = .true.)
    endif


    if ((globalnSource > 0).and.(dt > 0.d0).and.nBodyPhysics) &
         call doAccretion(grid, globalsourceArray, globalnSource, dt)


    if ((globalnSource > 0).and.(dt > 0.d0).and.nBodyPhysics) &
         call updateSourcePositions(globalsourceArray, globalnSource, dt, grid)

    if (myrankglobal == 1) call tune(6,"Boundary conditions")

    if (selfGravity) then
       if (myrankglobal == 1) call tune(6,"Self-gravity")
       call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call zeroSourcepotential(grid%octreeRoot)
       if (globalnSource > 0) then
          call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, grid%halfSmallestSubcell)
       endif
       call sumGasStarGravity(grid%octreeRoot)
       if (myrankglobal == 1) call tune(6,"Self-gravity")
    endif


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
    if(rhieChow) then
     call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!     call computepressureu(grid, grid%octreeroot, direction) !Thaw
     call computepressureGeneral(grid, grid%octreeroot, .false.)
     call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    
     call setuppressure(grid%octreeroot, grid, direction) !Thaw
     call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
     call rhiechowui(grid%octreeroot, grid, direction, dt/2.d0)
    end if

    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)   
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
   
    !call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    !call imposeboundary(grid%octreeroot)
    !call periodboundary(grid)
    !call transfertempstorage(grid%octreeroot)

    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    !call computePressureU(grid, grid%octreeRoot, direction)
    call computepressureGeneral(grid, grid%octreeroot, .true.)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupPressure(grid%octreeRoot, grid, direction)
    if(rhieChow) then
     call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
     call rhiechowui(grid%octreeroot, grid, direction, dt)
    end if
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)

    call pressureForceU(grid%octreeRoot, dt/2.d0)

    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    if(rhieChow) then
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
!       call computepressurew(grid, grid%octreeroot, direction) !Thaw
       call computepressureGeneral(grid, grid%octreeroot, .false.)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
       call setuppressure(grid%octreeroot, grid, direction) !Thaw
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
       call rhiechowui(grid%octreeroot, grid, direction, dt)
    end if
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    
    !call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)
    !call imposeboundary(grid%octreeroot)
    !call periodboundary(grid)
    !call transfertempstorage(grid%octreeroot)

    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupWi(grid%octreeRoot, grid, direction)
    call setupWpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    !call computePressureW(grid, grid%octreeRoot, direction)
    call computepressureGeneral(grid, grid%octreeroot, .true.)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    call setupPressure(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=3)
    if(rhieChow) then
       call rhiechowui(grid%octreeroot, grid, direction, dt)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)
    end if
    call pressureForceW(grid%octreeRoot, dt)

    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    if(rhieChow) then
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!       call computepressureu(grid, grid%octreeroot, direction) !Thaw
       call computepressureGeneral(grid, grid%octreeroot, .false.)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
       call setuppressure(grid%octreeroot, grid, direction) !Thaw
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
       call rhiechowui(grid%octreeroot, grid, direction, dt/2.d0)
    end if
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call advectRho(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoW(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    
    !call imposeboundary(grid%octreeroot)
    !call periodboundary(grid)
    !call transfertempstorage(grid%octreeroot)
    !call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)


    call setupUi(grid%octreeRoot, grid, direction)
    call setupUpm(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!    call computePressureU(grid, grid%octreeRoot, direction)
    call computepressureGeneral(grid, grid%octreeroot, .false.)
    call setupRhoPhi(grid%octreeRoot, grid, direction)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call setupPressure(grid%octreeRoot, grid, direction)
    if(rhieChow) then
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
       call rhiechowui(grid%octreeroot, grid, direction, dt/2.d0)
    end if
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)

    call pressureForceU(grid%octreeRoot, dt/2.d0)

    call imposeBoundary(grid%octreeRoot)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)
 
  end subroutine hydroStep2d


  recursive subroutine computeCourantTime(grid, thisOctal, tc)
    use mpi
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
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
                call computeCourantTime(grid, child, tc)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then

             cs = soundSpeed(thisOctal, subcell)
!             if (myrank == 1) write(*,*) "cs ", returnPhysicalUnitSpeed(cs)/1.d5, " km/s ",cs, " code"
!             dx = returnCodeUnitLength(thisOctal%subcellSize*gridDistanceScale)
             dx = grid%halfSmallestsubcell *gridDistanceScale* 2.d0

!Use max velocity not average
             speed = max(thisOctal%rhou(subcell)**2, thisOctal%rhov(subcell)**2, thisOctal%rhow(subcell)**2)
!             speed = thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2
             speed = sqrt(speed)/thisOctal%rho(subcell)
!             if (myrank == 1) write(*,*) "speed ", returnPhysicalUnitSpeed(speed)/1.d5, " km/s ",speed, " code"
!             if (myrank == 1) write(*,*) "dx ", returnPhysicalUnitLength(dx), " cm ",dx," code"
    

            ! if(myRank == 1) then
            !    tc = min(tc, dx / (cs + speed))
            !    write(*,*) "tc = min(tc, dx / (cs + speed) )", tc!

             !   tc = returnPhysicalUnitLength(dx)/ (returnPhysicalUnitSpeed(cs) + returnPhysicalUnitSpeed(speed))


            !    write(*,*) "min(tc, returnPhysicalUnitLength(dx)/(returnPhysicalUnitSpeed(cs) +", &
            !         "returnPhysicalUnitSpeed(speed))", tc
            ! end if
             
             tc = min(tc, dx / (cs + speed) )
!             if (myrank == 1) write(*,*) "tc ",tc

          endif
 
       endif
    enddo
  end subroutine computeCourantTime

  recursive subroutine sumGasStarGravity(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sumgasStarGravity(child)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          thisOctal%phi_i(subcell) = thisOctal%phi_gas(subcell) + thisOctal%phi_stars(subcell)
       endif
    enddo
  end subroutine sumGasStarGravity

  function soundSpeed(thisOctal, subcell) result (cs)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: cs, rhoPhys
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
          !!print *, "getPressure(thisOctal, subcell)", getPressure(thisOctal, subcell)
          !print *, "thisOctal%rho(subcell)", thisOctal%rho(subcell)
          !print *, "cs = ", cs
          !stop
          !if(cs > 2200000.) then
          !   print *, "cs ", cs
          !   stop
          !end if
       case(2) ! barotropic
          rhoPhys = returnPhysicalUnitDensity(thisOctal%rho(subcell))
          if (rhoPhys < rhoCrit) then
             cs = sqrt(getPressure(thisOctal, subcell)/thisOctal%rho(subcell))
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
    use inputs_mod, only : tStart, tEnd, tDump, dorefine, amrTolerance
    use mpi
    type(gridtype) :: grid
    real(double) :: dt,  gamma, mu
    real(double) :: currentTime
    integer :: i
    type(VECTOR) :: direction
    integer :: nHydroThreads
    real(double) :: tc(10), totalMass, totalEnergy
    integer :: it
    integer :: myRank, ierr
    integer :: nPairs, thread1(100), thread2(100), group(100), nBound(100), ngroup
    integer :: iUnrefine, nUnrefine
    real(double) :: nextDumpTime, temptc(10)
    character(len=80) :: plotfile
    real(double) :: iniM, endM, iniE, endE

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

       
       call evenUpGridMPI(grid,.false.,dorefine)
       
       if(dorefine) then
!          call refineGridGeneric(grid, 1.d-2)
          call refineGridGeneric(grid, amrTolerance)
       call writeInfo("Evening up grid", TRIVIAL)    
       end if

       call evenUpGridMPI(grid, .false., dorefine)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call setupX(grid%octreeRoot, grid, direction)
       call setupQX(grid%octreeRoot, grid, direction)
       !    call calculateEnergy(grid%octreeRoot, gamma, mu)
!       call calculateRhoE(grid%octreeRoot, direction) 
       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call calculateRhoU(grid%octreeRoot, direction)
       direction = VECTOR(0.d0, 1.d0, 0.d0)
       call calculateRhoV(grid%octreeRoot, direction)
       direction = VECTOR(0.d0, 0.d0, 1.d0)
       call calculateRhoW(grid%octreeRoot, direction)

    endif

    currentTime = tStart
    it = 0
    nextDumpTime = 0.d0
    iUnrefine = 0

          call findEnergyOverAllThreads(grid, totalenergy)
          call findMassOverAllThreads(grid, totalmass)
          iniM = totalMass
          iniE = totalEnergy
    do while(currentTime <= tend)
               
       tc = 0.d0
       if (myrank /= 0) then
          tc(myrank) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRank))
       endif

       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)

!       if (myrank == 1) write(*,*) "courantTime", dt, it
       if (myrank == 1) call tune(6,"Hydrodynamics step")
!       call writeInfo("calling hydro step",TRIVIAL)

       if (myrankGlobal /= 0) then
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          
          call hydroStep1d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
         
          if(dounrefine) then
             iUnrefine = iUnrefine + 1
             
             if (iUnrefine == 100) then
                if (myrankglobal == 1) call tune(6, "Unrefine grid")
                call unrefineCells(grid%octreeRoot, grid, nUnrefine, 1.d-3)
                if (myrankglobal == 1) call tune(6, "Unrefine grid")
                iUnrefine = 0
             endif
          end if
          call evenUpGridMPI(grid, .true., dorefine)
          if (myrank == 1) call tune(6,"Hydrodynamics step")
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call zeroRefinedLastTime(grid%octreeRoot)
          if(dorefine) then
             call refineGridGeneric(grid, amrTolerance)
           end if
           call evenUpGridMPI(grid, .true., dorefine)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       else

       end if

!          if (myRank == 1) write(*,*) "current time ",currentTime,dt
       currentTime = currentTime + dt
       if (currentTime .gt. nextDumpTime) then
          !Thaw
          !  write(plotfile,'(a,i4.4,a)') "gaussian",it,".dat"
          !call  dumpValuesAlongLine(grid, plotfile, VECTOR(0.d0,0.d0,0.0d0), &
          !VECTOR(1.d0, 0.d0, 0.0d0), 1000)
          write(plotfile,'(a,i4.4,a)') "sod.dat"
          call  dumpValuesAlongLine(grid, plotfile, &
               VECTOR(0.d0,0.d0,0.0d0), VECTOR(1.d0, 0.d0, 0.0d0), 1000)
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
          grid%iDump = it
          grid%currentTime = currentTime
       endif
   
    enddo


    call findEnergyOverAllThreads(grid, totalenergy)
    call findMassOverAllThreads(grid, totalmass)
    
    endM = totalMass
    endE = totalEnergy
    
    if(myRank == 1) then
       print *, "dM (%) = ", (1.d0 - (endM/iniM))*100.d0
       print *, "dE (%) = ", (1.d0 - (endE/iniE))*100.d0
       print *, "endM", endM
       print *, "iniM", iniM
       print *, "endE", endE
       print *, "iniE", iniE
    
       open(20, file="conservation.dat", status="unknown")
       write(20,*) iniM, endM, iniE, endE
       close(20)
    end if

  end subroutine doHydrodynamics1d


  subroutine doHydrodynamics3d(grid)
    use vtk_mod, only : writeVtkFilenBody
    use inputs_mod, only : tdump, tend, doRefine, doUnrefine, amrTolerance, dumpRadial
    use mpi
    type(gridtype) :: grid
    real(double) :: dt, tc(64), temptc(64),  mu
    real(double) :: currentTime !, smallTime
    real(double) :: initialmass
    integer :: i, it, iUnrefine
    integer :: myRank, ierr
    character(len=80) :: plotfile
    real(double) :: nextDumpTime, tff,totalMass !, ang
    type(VECTOR) :: direction, viewVec
    integer :: nUnrefine
    integer :: thread1(200), thread2(200), nBound(200), nPairs
    integer :: nHydroThreads
    integer :: nGroup, group(200)
!    logical :: doRefine
    logical :: doSelfGrav
    logical, save  :: firstStep = .true.

    doSelfGrav = .true.

    if (writeoutput) then
       open(57, file="pos.dat", status="unknown", form="formatted")
       close(57)
    endif

    if (grid%geometry == "shakara") doSelfGrav = .false.
    if (grid%geometry == "rtaylor") doSelfGrav = .false.
    if (grid%geometry == "diagSod") doSelfGrav = .false.

!    dorefine = .true.

    nHydroThreads = nThreadsGlobal - 1

    direction = VECTOR(1.d0, 0.d0, 0.d0)

    mu = 2.d0

    viewVec = VECTOR(-1.d0,0.d0,0.d0)
    !    viewVec = rotateZ(viewVec, 20.d0*degtorad)
    viewVec = rotateY(viewVec, 25.d0*degtorad)



    tff = 1.d0 / sqrt(bigG * (1d0*mSol/((4.d0/3.d0)*pi*7.d15**3)))
    if (tdump == 0.d0) then
       !    tDump = 0.01d0
       !tDump = 1.d-2 * tff
       
!       tff = sqrt((3.d0*pi)/(32.d0*bigG*3.466d-21))
       !tDump = 0.01d0 * tff
       tDump = 1.d0 * tff
    endif

    it = grid%iDump
    currentTime = grid%currentTime
    nextDumpTime = grid%currentTime+tdump

    call setCodeUnit(time=timeUnit)
    call setCodeUnit(mass=massUnit)
    call setCodeUnit(length=lengthUnit)
    call writeCodeUnits()

    tdump = returnCodeUnitTime(tdump)
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
    endif

!    if (readgrid) then
!       if(myrankglobal /= 0) then
!          call zeroPhiGas(grid%octreeRoot)
!          call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.) 
!       endif
!    endif


    if (.not.readgrid) then

    if (myrankGlobal /= 0) then



       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call calculateRhoU(grid%octreeRoot, direction)
       direction = VECTOR(0.d0, 1.d0, 0.d0)
       call calculateRhoV(grid%octreeRoot, direction)
       direction = VECTOR(0.d0, 0.d0, 1.d0)
       call calculateRhoW(grid%octreeRoot, direction)

       call calculateRhoInCodeUnits(grid%octreeRoot)
       call calculateEnergyInCodeUnits(grid%octreeRoot)
       call calculatePhiInCodeUnits(grid%octreeRoot)


       call calculateRhoE(grid%octreeRoot, direction)
    end if

    if (myrankGlobal /= 0) then


       if(doRefine) then
           if (myrank == 1) call tune(6, "Initial refine")
           call refineGridGeneric(grid, amrTolerance)
           call writeInfo("Evening up grid", TRIVIAL)    
       end if
   end if
   
    if(myRankGlobal /= 0) then

       if(dorefine) then
          call evenUpGridMPI(grid,.false., dorefine)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          
          if (myrank == 1) call tune(6, "Initial refine")
       end if

    end if

    if(myRankGlobal /= 0) then
       
       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call setupX(grid%octreeRoot, grid, direction)
       call setupQX(grid%octreeRoot, grid, direction)

       if (doselfGrav) then

          call findMassOverAllThreads(grid, totalmass)
          if (writeoutput) write(*,*) "Total mass: ",totalMass/msol, " solar masses"

          if (myrank == 1) call tune(6, "Self-Gravity")
          if (myrank == 1) write(*,*) "Doing multigrid self gravity"
          call writeVtkFile(grid, "beforeselfgrav.vtk", &
               valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " /))
          call zeroPhiGas(grid%octreeRoot)
          call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.) 

          call zeroSourcepotential(grid%octreeRoot)
          if (globalnSource > 0) then
             call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, grid%halfSmallestSubcell)
          endif
          call sumGasStarGravity(grid%octreeRoot)

          call writeVtkFile(grid, "afterselfgrav.vtk", &
               valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", &
               "phigas       ","phi          " /))
          if (myrank == 1) write(*,*) "Done"
          if (myrank == 1) call tune(6, "Self-Gravity")
       endif

    endif
 endif

    tc = 0.d0
    if (myrank /= 0) then
       tc(myrank) = 1.d30
       call computeCourantTime(grid, grid%octreeRoot, tc(myRank))
    endif
    call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (firstStep) then
       firstStep = .false.
       dt = MINVAL(temptc(1:nHydroThreads)) * 1.d-5
    else
       dt = MINVAL(temptc(1:nHydroThreads)) * dble(cflNumber)
    endif
    if (writeoutput) write(*,*) "Courant time is ",dt

    if ((grid%geometry == "shakara").or.(grid%geometry=="rtaylor")) then
       tdump =  20.d0 * dt
    endif

    if (myRank == 1) write(*,*) "CFL set to ", cflNumber

    if (myrankglobal==1) write(*,*) "Setting tdump to: ", tdump
    if (it ==0) nextDumpTime = 0.
    iUnrefine = 0


!    write(plotfile,'(a)') "start.vtk"
!    call writeVtkFile(grid, plotfile, &
!         valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phi          " /))

    call findMassoverAllThreads(grid, initialMass)

       if (nBodyPhysics) then
          initialMass = initialMass + SUM(globalSourceArray(1:globalnSource)%mass)
       endif

    do while(currentTime < tend)
       if (myrank /= 0) then
          tc(myrank) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRank))
       endif
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       !       write(*,*) "tc", tc(1:8)
       !       if (myrank==1)write(*,*) "temp tc",temptc(1:nHydroThreads)

!       if (firstStep) then 
!          cflNumber = 1.d-2
!          firstStep = .false.
!       else
!          cflNumber = 0.3d0
!       endif

       dt = MINVAL(temptc(1:nHydroThreads)) * dble(cflNumber)
       if (writeoutput) write(*,*) "Courant time is ",dt
       if (myrank==1)write(*,*) "Current time is ",returnPhysicalUnitTime(currentTime)
       call findMassoverAllThreads(grid, totalmass)
       if (nBodyPhysics) then
          totalmass = totalMass + SUM(globalSourceArray(1:globalnSource)%mass)
       endif
       if (myrank==1)write(*,*) "Current mass: ",totalmass/initialmass,initialMass/msol

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

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


          call hydroStep3d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup, doSelfGrav=doSelfGrav)

          if (nbodyPhysics) call addSinks(grid, globalsourceArray, globalnSource)

          if (nbodyPhysics) call mergeSinks(grid, globalsourceArray, globalnSource)


          if (myrank == 1) call tune(6,"Hydrodynamics step")


          iUnrefine = iUnrefine + 1
          if(doUnRefine) then
             if (iUnrefine == 5) then
                if (myrank == 1)call tune(6, "Unrefine grid")
               ! nUnrefine = 0
                call unrefineCells(grid%octreeRoot, grid, nUnrefine, 5.d-3)
                                !          write(*,*) "Unrefined ", nUnrefine, " cells"
                if (myrank == 1)call tune(6, "Unrefine grid")
                iUnrefine = 0
             endif
          end if

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          if(doRefine) then
             call writeInfo("Refining grid part 2", TRIVIAL)    
             call refineGridGeneric(grid, amrTolerance)
          !          
             call writeInfo("Done the refine part", TRIVIAL)                 
             call evenUpGridMPI(grid, .true., dorefine)
             call writeInfo("Done the even up part", TRIVIAL)    

          !       if (doSelfGrav) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)

             if (myrank == 1) call tune(6, "Loop refine")
             
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
             
             if (myrank == 1) call tune(6, "Loop refine")
                                !
          end if

       endif

       currentTime = currentTime + dt
       if (myRank == 1) write(*,*) "current time ",currentTime,dt,nextDumpTime
       if (myRank == 1) write(*,*) "percent to next dump ",100.*(nextDumpTime-currentTime)/tdump


       if (currentTime .ge. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1

          grid%iDump = it
          grid%currentTime = currentTime

          if (dumpRadial) then
             write(plotfile,'(a,i4.4,a)') "radial",it,".dat"
             call  dumpValuesAlongLine(grid, plotfile, VECTOR(0.d0,0.d0,0.0d0), &
                  VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0),1000)
          endif
          write(plotfile,'(a,i4.4,a)') "dump",it,".grid"
          call writeAMRgrid(plotfile,.false. ,grid)

          if (writeoutput) then
             write(plotfile,'(a,i4.4,a)') "source",it,".dat"
             call writeSourceArray(plotfile)
          endif

          write(plotfile,'(a,i4.4,a)') "output",it,".vtk"
          call writeVtkFile(grid, plotfile, &
               valueTypeString=(/"rho          ",&
               "hydrovelocity", &
               "rhoe         ", &
               "u_i          ", &
               "phi          "/))

          write(plotfile,'(a,i4.4,a)') "nbody",it,".vtk"
          call writeVtkFilenBody(globalnSource, globalsourceArray, plotfile, grid)
          if (myrank==1) write(*,*) trim(plotfile), " written at ",currentTime/tff, " free-fall times"
       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)


    enddo
  end subroutine doHydrodynamics3d

  subroutine doHydrodynamics2d(grid)
    use inputs_mod, only : tEnd, tDump, doRefine, doUnrefine, amrTolerance
    use mpi
    type(gridtype) :: grid
    real(double) :: dt, tc(64), temptc(64), mu
    real(double) :: currentTime
    integer :: i, it, iUnrefine
    integer :: myRank, ierr
    character(len=80) :: plotfile
    real(double) :: nextDumpTime, tff!, ang
    real(double) :: totalEnergy, totalMass
    type(VECTOR) :: direction, viewVec
    integer :: thread1(100), thread2(100), nBound(100), nPairs
    integer :: nGroup, group(100)
    integer :: nHydroThreads 
    logical :: converged
    integer :: nUnrefine, jt

    nUnrefine = 0

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    nHydroThreads = nThreadsGlobal - 1

    it = grid%iDump
    currentTime = grid%currentTime
    nextDumpTime = 0.d0

    if (it /= 1) then
       call writeVTKfile(grid, "readin.vtk")
    endif
    if (myrankGlobal /= 0) then

       !famr
       call refineEdges(grid%octreeRoot, grid,  converged)

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


       if (it == 0) then
          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call calculateRhoU(grid%octreeRoot, direction)
          direction = VECTOR(0.d0, 1.d0, 0.d0)
          call calculateRhoV(grid%octreeRoot, direction)
          direction = VECTOR(0.d0, 0.d0, 1.d0)
          call calculateRhoW(grid%octreeRoot, direction)

          !    call calculateEnergy(grid%octreeRoot, gamma, mu)
          call calculateRhoE(grid%octreeRoot, direction)



          call evenUpGridMPI(grid,.false., dorefine)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          if(doRefine) then
             call refinegridGeneric(grid, amrTolerance)          
          end if
          call evenUpGridMPI(grid, .false.,dorefine)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call setupX(grid%octreeRoot, grid, direction)
          call setupQX(grid%octreeRoot, grid, direction)
          !    call calculateEnergy(grid%octreeRoot, gamma, mu)

          call writeInfo("Checking Boundary Partner Vectors", TRIVIAL)
          call checkBoundaryPartners(grid%octreeRoot, grid)
          call writeInfo("Initial Boundary Partner Check Passed", TRIVIAL)
       endif
    endif



    tc = 0.d0
    if (myrank /= 0) then
       tc(myrank) = 1.d30
       call computeCourantTime(grid, grid%octreeRoot, tc(myRank))
    endif
    call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierr)
    tc = tempTc
    dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)

    tff = 1.d0 / sqrt(bigG * (0.1d0*mSol/((4.d0/3.d0)*pi*7.d15**3)))

    if (writeoutput) write(*,*) "Setting tdump to: ", tdump

    if (it /= 0) then
       nextDumpTime = grid%currentTime + tdump
    endif


    iUnrefine = 0
    !    call writeInfo("Plotting grid", TRIVIAL)    


    jt = 0

    !Thaw - trace courant time history
    open (444, file="tcHistory.dat", status="unknown")

    do while(currentTime < tEnd)
       if (myrank == 1) write(*,*) "current time " ,currentTime
       jt = jt + 1
       tc = 0.d0
       if (myrank /= 0) then
          tc(myrank) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRank))
       endif
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierr)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)
       write(444, *) jt, MINVAL(tc(1:nHydroThreads)), dt

       !       if ((jt < 2000).and.(grid%geometry=="sedov")) then
       !          dt = MINVAL(tc(1:nHydroThreads)) * 1.d-4
       !       endif

       if ((currentTime + dt).gt.tEnd) then
          nextDumpTime = tEnd
          dt = nextDumpTime - currentTime
       endif

       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

       if (myrankGlobal /= 0) then
          if (myrank == 1) write(*,*) "courantTime", dt, it
          if (myrank == 1) call tune(6,"Hydrodynamics step")
          call writeInfo("calling hydro step",TRIVIAL)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

          call hydroStep2d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
       end if

       call findEnergyOverAllThreads(grid, totalenergy)
       if (writeoutput) write(*,*) "Total energy: ",totalEnergy
       call findMassOverAllThreads(grid, totalmass)
       if (writeoutput) write(*,*) "Total mass: ",totalMass

       if(doUnrefine) then
          iUnrefine = iUnrefine + 1
          if (iUnrefine == 20) then
             if (myrankglobal == 1) call tune(6, "Unrefine grid")
             call unrefineCells(grid%octreeRoot, grid, nUnrefine, 5.d-3)
             if (myrankglobal == 1) call tune(6, "Unrefine grid")
             iUnrefine = 0
          endif
       end if

       if (myrank == 1) call tune(6,"Hydrodynamics step")
!Thaw - temprorary extra evenup
       call evenUpGridMPI(grid, .true., dorefine) !, dumpfiles=jt)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       if(doRefine) then
          call refinegridGeneric(grid, amrTolerance)
       end if
       call evenUpGridMPI(grid, .true., dorefine) !, dumpfiles=jt)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       currentTime = currentTime + dt
       !       if (myRank == 1) write(*,*) "current time ",currentTime,dt


       !Perform another boundary partner check
       call checkBoundaryPartners(grid%octreeRoot, grid)


       if (currentTime .ge. nextDumpTime) then
          it = it + 1
          nextDumpTime = nextDumpTime + tDump
          grid%iDump = it
          grid%currentTime = currentTime
          call hydrovelocityConvert(grid%octreeRoot)
          write(plotfile,'(a,i4.4,a)') "dump",it,".grid"
          call writeAMRgrid(plotfile,.false. ,grid)
          write(plotfile,'(a,i4.4,a)') "output",it,".vtk"
          call writeVtkFile(grid, plotfile, &
               valueTypeString=(/"rho          ","velocity     ","rhoe         " , &
               "u_i          ", &
               "hydrovelocity", &
               "rhou         ", &
               "rhov         ", &
               "rhow         ", &
               "phi          "/))
          if (grid%geometry == "sedov") &
               call dumpValuesAlongLine(grid, "sedov.dat", VECTOR(0.5d0,0.d0,0.0d0), VECTOR(0.9d0, 0.d0, 0.0d0), 1000)

       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)


    enddo
    close(444)

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
       endif
    enddo
  end subroutine zeroRefinedLastTime


!  recursive subroutine printRefinedLastTime(thisOctal)
!    type(octal), pointer   :: thisOctal
!    type(octal), pointer  :: child
!    integer :: subcell, i!
!
!    do subcell = 1, thisOctal%maxChildren
!       if (thisOctal%hasChild(subcell)) then
!          ! find the child
!          do i = 1, thisOctal%nChildren, 1
!             if (thisOctal%indexChild(i) == subcell) then
!                child => thisOctal%child(i)
!                call zeroRefinedLastTime(child)
!                exit
!             end if
!          end do
!       else
!          if(thisOctal%refinedLastTime) then
!       endif
!    enddo
!  end subroutine printRefinedLastTime



  recursive subroutine zeroPhigas(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroPhiGas(child)
                exit
             end if
          end do
       else 
          thisOctal%phi_gas = 0.d0
       endif
    enddo
  end subroutine zeroPhigas


  recursive subroutine transferTempStorage(thisOctal, justGrav)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    logical, optional :: justGrav
    logical :: doJustGrav

    dojustgrav = .false.
    if (PRESENT(justGrav)) doJustGrav = justGrav
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call transferTempStorage(child, justGrav)
                exit
             end if
          end do
       else 
          if (thisOctal%ghostCell(subcell)) then
             if (.not.associated(thisOctal%tempStorage)) then
                call writeFatal("Tempstorage not allocated when it should have been")
                stop
             endif
             if (.not.doJustGrav) then
                thisOctal%rho(subcell) = thisOctal%tempStorage(subcell,1)
                thisOctal%rhoe(subcell) = thisOctal%tempStorage(subcell,2)
                thisOctal%rhou(subcell) = thisOctal%tempStorage(subcell,3)
                thisOctal%rhov(subcell) = thisOctal%tempStorage(subcell,4)
                thisOctal%rhow(subcell) = thisOctal%tempStorage(subcell,5)
                thisOctal%energy(subcell) = thisOctal%tempStorage(subcell,6)
                thisOctal%pressure_i(subcell) = thisOctal%tempStorage(subcell,7)
             else
                thisOctal%phi_i(subcell) = thisOctal%tempStorage(subcell,1)
             endif
          endif
       endif
    enddo
    if (associated(thisOCtal%tempStorage)) Deallocate(thisOctal%tempStorage)
  end subroutine transferTempStorage

  recursive subroutine transferTempStorageLevel(thisOctal, nDepth, justGrav)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: nDepth
    integer :: subcell, i
    logical, optional :: justGrav
    logical :: doJustGrav

    doJustGrav = .false.
    if (PRESENT(justGrav)) doJustGrav = justGrav

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call transferTempStorageLevel(child, nDepth, justGrav)
       end do
    else

       do subcell = 1, thisOctal%maxChildren
          if (thisOctal%ghostCell(subcell)) then
             if (.not.associated(thisOctal%tempStorage)) then
                call writeFatal("Tempstorage not allocated when it should have been")
                stop
             endif
             if (.not.doJustGrav) then
                thisOctal%rho(subcell) = thisOctal%tempStorage(subcell,1)
                thisOctal%rhoe(subcell) = thisOctal%tempStorage(subcell,2)
                thisOctal%rhou(subcell) = thisOctal%tempStorage(subcell,3)
                thisOctal%rhov(subcell) = thisOctal%tempStorage(subcell,4)
                thisOctal%rhow(subcell) = thisOctal%tempStorage(subcell,5)
                thisOctal%energy(subcell) = thisOctal%tempStorage(subcell,6)
                thisOctal%pressure_i(subcell) = thisOctal%tempStorage(subcell,7)
             else
                thisOctal%phi_i(subcell) = thisOctal%tempStorage(subcell,1)
             endif
          endif
       enddo
    endif
    if (associated(thisOCtal%tempStorage)) Deallocate(thisOctal%tempStorage)
  end subroutine transferTempStorageLevel

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

          thisOctal%rhoV(subcell) = returnCodeUnitRhoV(thisOctal%rhoV(subcell))

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

          thisOctal%rhoW(subcell) = returnCodeUnitRhoV(thisOctal%rhoW(subcell))

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

          thisOctal%rhoU(subcell) = returnCodeUnitRhoV(thisOctal%rhoU(subcell))


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

       endif
    enddo
  end subroutine calculateRhoE

  recursive subroutine calculateRhoInCodeUnits(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: i
  
    do i = 1, thisOctal%maxChildren
       thisOctal%rho(i) = returnCodeUnitDensity(thisOctal%rho(i))
    enddo


    IF ( thisOctal%nChildren > 0 ) THEN
       ! call this subroutine recursively on each of its children
       DO i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          CALL calculateRhoInCodeUnits(child)
       END DO
    END IF
  end subroutine calculateRhoInCodeUnits

  recursive subroutine calculatePhiInCodeUnits(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: i
  
    do i = 1, thisOctal%maxChildren
       thisOctal%Phi_i(i) = returnCodeUnitEnergy(thisOctal%phi_i(i))
    enddo


    IF ( thisOctal%nChildren > 0 ) THEN
       ! call this subroutine recursively on each of its children
       DO i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          CALL calculatePhiInCodeUnits(child)
       END DO
    END IF
  end subroutine calculatePhiInCodeUnits

  recursive subroutine calculateEnergyinCodeUnits(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateEnergyinCodeUnits(child)
                exit
             end if
          end do
       else 

          thisOctal%energy(subcell) = returnCodeUnitEnergy(thisOctal%energy(subcell))
       endif
    enddo
  end subroutine calculateEnergyInCodeUnits


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
    use inputs_mod, only : fixedRhoBound, rho_const
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

! NB confusion regarding 2d being x,z rather than x,y

                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (dir%x > 0.9d0) then 
                         thisOctal%tempStorage(subcell,3) = -abs(bOctal%rhou(bSubcell))
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

                case(5) ! constant density, not velocity reflection
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

! NB confusion regarding 2d being x,z rather than x,y

                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if ((abs(dir%x) > 0.5d0).and.abs(dir%z) > 0.5d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                   else if (thisOctal%threed) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (abs(dir%y) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if ((abs(dir%x) > 0.2d0).and.abs(dir%z) > 0.2d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                   endif

                case(6) ! inflow boundary condition

                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocal(locator, bOctal, bSubcell)

                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)
                   Thisoctal%tempstorage(subcell,1) = inflowRho
                   thisOctal%tempStorage(subcell,2) = inflowRhoE

                   thisOctal%tempStorage(subcell,6) = inflowEnergy
                   thisOctal%tempStorage(subcell,7) = inflowPressure

! NB confusion regarding 2d being x,z rather than x,y

                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = inflowMomentum
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = inflowMomentum
                      endif
                   else if (thisOctal%threed) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = inflowMomentum
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%y) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = inflowMomentum
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = inflowMomentum
                      endif
                   endif
                   

                case DEFAULT
                   write(*,*) "Unrecognised boundary condition in impose boundary: ", thisOctal%boundaryCondition(subcell)
             end select
             
             if(fixedRhoBound) then
                thisOctal%tempStorage(subcell,1) = rho_const
             end if

          endif
      endif
    enddo
  end subroutine imposeBoundary



  recursive subroutine setupGhostCells(thisOctal, grid, flag)
    use mpi
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
                   case(1, 5, 6)
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


                   ! now do the gravity locator
                   dx = thisOctal%subcellSize
                   thisOctal%ghostCell(subcell) = .true.
                   
                   locator = subcellCentre(thisOctal, subcell) - &
                        (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)*probe(iProbe)
                   thisOctal%GravboundaryPartner(subcell) = locator

                   tempSubcell = 1
                   tempOctal => thisOctal
                   call findSubcellLocal(locator, tempOctal, tempSubcell)
                   if (tempOctal%mpiThread(tempSubcell) == myRankGlobal) then
                      write(*,*) "locator problem ",myRankGlobal
                      stop
                   endif


                   locator = subcellCentre(neighbourOctal, neighboursubcell) - &
                        (grid%octreeRoot%subcellSize*2.d0-4.d0*neighbourOctal%subcellSize)*probe(iProbe)
                   neighbourOctal%GravboundaryPartner(neighboursubcell) = locator
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
             case(1, 4, 5, 6)
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
                   write(*,*) "Unknown boundary condition in setupghostcells2 A: ",thisOctal%boundaryCondition(subcell)
             end select

             ! gravity partner
             dx = thisOctal%subcellSize
             thisOctal%ghostCell(subcell) = .true.
             
             locator = subcellCentre(thisOctal, subcell) - &
                  (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)*direction
             thisOctal%GravboundaryPartner(subcell) = locator
             tempSubcell = 1
             tempOctal => thisOctal
             call findSubcellLocal(locator, tempOctal, tempSubcell)
             if (tempOctal%mpiThread(tempSubcell) == myRankGlobal) then
                write(*,*) "locator problem ",myRankGlobal
                stop
             endif
             
                   
             locator = subcellCentre(neighbourOctal, neighboursubcell) - &
                  (grid%octreeRoot%subcellSize*2.d0-4.d0*neighbourOctal%subcellSize)*direction
             neighbourOctal%GravboundaryPartner(neighboursubcell) = locator
             
             tempSubcell = 1
             tempOctal => thisOctal
             call findSubcellLocal(locator, tempOctal, tempSubcell)
             if (tempOctal%mpiThread(tempSubcell) == myRankGlobal) then
                write(*,*) "locator problem ",myRankGlobal
                stop
             endif

          endif
      endif
    enddo
  end subroutine setupGhostCells


  recursive subroutine setupEdges(thisOctal, grid)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: locator, rVec
    integer :: nProbes, iProbe
    type(VECTOR) :: probe(6)
    integer :: nProbeOutside
    integer :: myRank, ierr
    character(len=10) :: boundary

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    if (myrankGlobal == 0) goto 666
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
          if (.not.associated(thisOctal%corner)) then
             allocate(thisOctal%corner(1:thisOctal%maxChildren))
             thisOctal%corner = .false.
          endif
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
             if (.not.inOctal(grid%octreeRoot, locator)) then
                nProbeOutside = nProbeOutside + 1
                thisOctal%boundaryPartner(subcell) = (-1.d0)*probe(iProbe)


                if (thisOctal%oneD) then
                   if (iProbe == 1) then
                      boundary = "xplus"
                   else
                      boundary = "xminus"
                   endif
                endif
                if (thisOctal%twoD) then
                   select case (iProbe)
                      case(1) 
                         boundary = "xplus"
                      case(2)
                         boundary = "xminus"
                      case(3)
                         boundary = "zplus"
                      case(4)
                         boundary = "zminus"
                   end select
                endif
                if (thisOctal%threeD) then
                   select case (iProbe)
                      case(1) 
                         boundary = "xplus"
                      case(2)
                         boundary = "xminus"
                      case(3)
                         boundary = "yplus"
                      case(4)
                         boundary = "yminus"
                      case(5)
                         boundary  = "zplus"
                      case(6)
                         boundary = "zminus"
                   end select
                endif
             endif
          enddo
 
          if (nProbeOutside >= 1) then
             thisOctal%edgeCell(subcell) = .true.
             thisOctal%boundaryCondition(subcell) = getBoundary(boundary)

             if(thisOctal%twoD .and. nProbeOutside == 2) then
                thisOctal%corner(subcell) = .true.
             else if (thisOctal%threeD .and. nProbeOutside == 3) then
                thisOctal%corner(subcell) = .true.
             end if
                
          endif
       endif
    enddo
666 continue
  end subroutine setupEdges

  recursive subroutine setupEdgesLevel(thisOctal, grid, nDepth)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    type(VECTOR) :: locator, rVec
    integer :: nProbes, iProbe
    integer :: nDepth
    type(VECTOR) :: probe(6)
    integer :: nProbeOutside
    integer :: myRank, ierr
    character(len=10) :: boundary

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    if (myrankGlobal == 0) goto 666

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call  setupEdgesLevel(child, grid, nDepth)
       end do
       else
       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

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
             if (.not.inOctal(grid%octreeRoot, locator)) then
                nProbeOutside = nProbeOutside + 1
                thisOctal%boundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                thisOctal%gravboundaryPartner(subcell) = (-1.d0)*probe(iProbe)


                if (thisOctal%oneD) then
                   if (iProbe == 1) then
                      boundary = "xplus"
                   else
                      boundary = "xminus"
                   endif
                endif
                if (thisOctal%twoD) then
                   select case (iProbe)
                      case(1) 
                         boundary = "xplus"
                      case(2)
                         boundary = "xminus"
                      case(3)
                         boundary = "zplus"
                      case(4)
                         boundary = "zminus"
                   end select
                endif
                if (thisOctal%threeD) then
                   select case (iProbe)
                      case(1) 
                         boundary = "xplus"
                      case(2)
                         boundary = "xminus"
                      case(3)
                         boundary = "yplus"
                      case(4)
                         boundary = "yminus"
                      case(5)
                         boundary  = "zplus"
                      case(6)
                         boundary = "zminus"
                   end select
                endif
             endif
          enddo
          if (nProbeOutside >= 1) then
             thisOctal%edgeCell(subcell) = .true.
             thisOctal%boundaryCondition(subcell) = getBoundary(boundary)
          endif
       enddo
    endif
666 continue
  end subroutine setupEdgesLevel

  recursive subroutine setupGhosts(thisOctal, grid)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal, tempOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell, tempSubcell
    type(VECTOR) :: locator, rVec, tVec
    integer :: nProbes, iProbe
    type(VECTOR) :: probe(6), currentDirection
    integer :: nProbeOutside
    integer :: myRank, ierr


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    if (myrankGlobal ==0 ) goto 666
    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupGhosts(child,  grid)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          thisOctal%ghostCell(subcell) = .false.

          if (.not.thisOctal%edgeCell(subcell)) then
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
                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighboursubcell)
                if (neighbourOctal%edgeCell(neighbourSubcell)) then
                   tVec = subcellCentre(thisOctal, subcell) + &
                           (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)* &
                           ((-1.d0)*probe(iProbe))
                   if (inOctal(grid%octreeRoot, tVec)) then
                      thisOctal%ghostCell(subcell) = .true.
                      thisOctal%boundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                      thisOctal%gravboundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                      currentDirection = (-1.d0)*probe(iProbe)
                      thisOctal%boundaryCondition(subcell) = &
                           neighbourOctal%boundaryCondition(neighbourSubcell)
                      exit
                   endif
                endif
             enddo
          else
             thisOctal%ghostCell(subcell) = .true.

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

             do iProbe = 1, nProbes
                locator = rVec + &
                     (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)
                if (.not.inOctal(grid%octreeRoot, locator)) then
                   thisOctal%boundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                   thisOctal%gravboundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                   currentDirection = (-1.d0)*probe(iProbe)
                   exit
                endif
             enddo
             



          endif

          if (thisOctal%ghostCell(subcell)) then
             select case(thisOctal%boundaryCondition(subcell))
             case(1, 4, 5, 6)
                
                if (thisOctal%edgeCell(subcell)) then
                   call locatorToNeighbour(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 3, locator)
                   thisOctal%boundaryPartner(subcell) = locator
                   tempOctal => thisOctal
                   tempSubcell = 1
                   call findSubcellLocal(locator, tempOctal, tempSubcell)
                   tempOctal%feederCell(tempsubcell) = .true.
                endif
                if (.not.thisOctal%edgeCell(subcell)) then
                   call locatorToNeighbour(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 1, locator)
                   thisOctal%boundaryPartner(subcell) = locator
                   tempOctal => thisOctal
                   tempSubcell = 1
                   call findSubcellLocal(locator, tempOctal, tempSubcell)
                   tempOctal%feederCell(tempsubcell) = .true.
                endif


                case(2)
                   
                   if (thisOctal%edgeCell(subcell)) then
                      locator = subcellCentre(thisOctal, subcell) + &
                           (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize) &
                           * thisOctal%boundaryPartner(subcell)
                      thisOctal%boundaryPartner(subcell) = locator
                      if (.not.(inOctal(grid%octreeRoot,locator))) then
                         write(*,*) "BUG1: locator ",locator
                      endif

                      tempOctal => thisOctal
                      tempSubcell = 1
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      tempOctal%feederCell(tempsubcell) = .true.
                   endif

                   
                   if (.not.thisOctal%edgeCell(subcell)) then
                      tVec = thisOctal%boundaryPartner(Subcell)
                      locator = subcellCentre(thisOctal, subcell) + &
                           (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)* &
                           thisOctal%boundaryPartner(subcell)
                      thisOctal%boundaryPartner(subcell) = locator
                      if (.not.(inOctal(grid%octreeRoot,locator))) then
                         write(*,*) "BUG2: locator ",locator
                         write(*,*) "cell centre ", subcellCentre(thisOctal,subcell)
                         write(*,*) "direction ",tVec
                      endif
                      tempOctal => thisOctal
                      tempSubcell = 1
                      call findSubcellLocal(locator, tempOctal, tempSubcell)
                      tempOctal%feederCell(tempsubcell) = .true.
                   endif

                   
                case DEFAULT
                   write(*,*) "Unknown boundary condition in setupghostcells2 B: ",thisOctal%boundaryCondition(subcell)
             end select


             ! gravity boundary

             if (thisOctal%edgeCell(subcell)) then
                tVec = thisOctal%gravBoundaryPartner(subcell)
                locator = subcellCentre(thisOctal, subcell) + &
                     (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize) &
                     * currentDirection
                thisOctal%GravboundaryPartner(subcell) = locator

                if (.not.(inOctal(grid%octreeRoot,locator))) then
                   write(*,*) "ordinary BUG1: locator ",locator, tVec
                endif
                
                tempOctal => thisOctal
                tempSubcell = 1
                call findSubcellLocal(locator, tempOctal, tempSubcell)
                tempOctal%feederCell(tempsubcell) = .true.
             endif

                   
             if (.not.thisOctal%edgeCell(subcell)) then
                tVec = thisOctal%gravboundaryPartner(Subcell)
                locator = subcellCentre(thisOctal, subcell) + &
                     (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)* &
                     currentDirection
                thisOctal%GravboundaryPartner(subcell) = locator
                if (.not.(inOctal(grid%octreeRoot,locator))) then
                   write(*,*) "BUG2: locator ",locator
                   write(*,*) "cell centre ", subcellCentre(thisOctal,subcell)
                   write(*,*) "direction ",tVec
                endif
                tempOctal => thisOctal
                tempSubcell = 1
                call findSubcellLocal(locator, tempOctal, tempSubcell)
                tempOctal%feederCell(tempsubcell) = .true.
             endif


          endif
       endif
    enddo
666 continue
  end subroutine setupGhosts

  recursive subroutine setupGhostsLevel(thisOctal, grid, nDepth)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, neighbourOctal, tempOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell, tempSubcell
    type(VECTOR) :: locator, rVec, tVec, currentDirection
    integer :: nProbes, iProbe, nDepth
    type(VECTOR) :: probe(6)
    integer :: nProbeOutside
    integer :: myRank, ierr


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    if (myrankGlobal ==0 ) goto 666
    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call  setupGhostsLevel(child, grid, nDepth)
       end do
       else
       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (.not.thisOctal%edgeCell(subcell)) then
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
                neighbourOctal => thisOctal
                call findSubcellLocalLevel(locator, neighbourOctal, neighboursubcell, nDepth)
                if (neighbourOctal%edgeCell(neighbourSubcell)) then
                   tVec = subcellCentre(thisOctal, subcell) + &
                           (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)* &
                           ((-1.d0)*probe(iProbe))
                   if (inOctal(grid%octreeRoot, tVec)) then
                      thisOctal%ghostCell(subcell) = .true.
                      thisOctal%boundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                      thisOctal%gravboundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                      thisOctal%boundaryCondition(subcell) = &
                           neighbourOctal%boundaryCondition(neighbourSubcell)
                      exit
                   endif
                endif
             enddo
          else
             thisOctal%ghostCell(subcell) = .true.
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

             do iProbe = 1, nProbes
                locator = rVec + &
                     (thisOctal%subcellsize/2.d0 + 0.01d0*grid%halfSmallestSubcell)*probe(iProbe)
                if (.not.inOctal(grid%octreeRoot, locator)) then
                   thisOctal%boundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                   thisOctal%gravboundaryPartner(subcell) = (-1.d0)*probe(iProbe)
                   currentDirection = (-1.d0)*probe(iProbe)
                   exit
                endif
             enddo
          endif

          if (thisOctal%ghostCell(subcell)) then
             select case(thisOctal%boundaryCondition(subcell))
             case(1, 4, 5, 6)
                
                if (thisOctal%edgeCell(subcell)) then
                   call locatorToNeighbourLevel(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 3, locator, nDepth)
                   thisOctal%boundaryPartner(subcell) = locator
                   tempOctal => thisOctal
                   tempSubcell = 1
                   call findSubcellLocalLevel(locator, tempOctal, tempSubcell, nDepth)
                   tempOctal%feederCell(tempsubcell) = .true.
                endif
                if (.not.thisOctal%edgeCell(subcell)) then
                   call locatorToNeighbourLevel(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 1, locator, nDepth)
                   thisOctal%boundaryPartner(subcell) = locator
                   tempOctal => thisOctal
                   tempSubcell = 1
                   call findSubcellLocalLevel(locator, tempOctal, tempSubcell, nDepth)
                   tempOctal%feederCell(tempsubcell) = .true.
                endif


                case(2)
                   
                   if (thisOctal%edgeCell(subcell)) then
                      locator = subcellCentre(thisOctal, subcell) + &
                           (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize) &
                           * thisOctal%boundaryPartner(subcell)
                      thisOctal%boundaryPartner(subcell) = locator
                      if (.not.(inOctal(grid%octreeRoot,locator))) then
                         write(*,*) "BUG1: locator ",locator
                      endif

                      tempOctal => thisOctal
                      tempSubcell = 1
                      call findSubcellLocalLevel(locator, tempOctal, tempSubcell, nDepth)
                      tempOctal%feederCell(tempsubcell) = .true.
                   endif

                   
                   if (.not.thisOctal%edgeCell(subcell)) then
                      tVec = thisOctal%boundaryPartner(Subcell)
                      locator = subcellCentre(thisOctal, subcell) + &
                           (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)* &
                           thisOctal%boundaryPartner(subcell)
                      thisOctal%boundaryPartner(subcell) = locator
                      if (.not.(inOctal(grid%octreeRoot,locator))) then
                         write(*,*) "BUG2: locator ",locator
                         write(*,*) "cell centre ", subcellCentre(thisOctal,subcell)
                         write(*,*) "direction ",tVec
                      endif
                      tempOctal => thisOctal
                      tempSubcell = 1
                      call findSubcellLocalLevel(locator, tempOctal, tempSubcell, nDepth)
                      tempOctal%feederCell(tempsubcell) = .true.
                   endif

                   
                case DEFAULT
                   write(*,*) "Unknown boundary condition in setupghostcells2 C: ",thisOctal%boundaryCondition(subcell)
             end select


             if (thisOctal%edgeCell(subcell)) then
                locator = subcellCentre(thisOctal, subcell) + &
                     (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize) &
                     * thisOctal%gravboundaryPartner(subcell)
                thisOctal%gravboundaryPartner(subcell) = locator
                if (.not.(inOctal(grid%octreeRoot,locator))) then
                   write(*,*) "level BUG1: locator ",locator
                endif

                tempOctal => thisOctal
                tempSubcell = 1
                call findSubcellLocalLevel(locator, tempOctal, tempSubcell, nDepth)
                tempOctal%feederCell(tempsubcell) = .true.
             endif

                   
             if (.not.thisOctal%edgeCell(subcell)) then
                tVec = thisOctal%gravboundaryPartner(Subcell)
                locator = subcellCentre(thisOctal, subcell) + &
                     (grid%octreeRoot%subcellSize*2.d0-4.d0*thisOctal%subcellSize)* &
                     thisOctal%gravboundaryPartner(subcell)
                thisOctal%gravboundaryPartner(subcell) = locator
                if (.not.(inOctal(grid%octreeRoot,locator))) then
                   write(*,*) "level BUG2: locator ",locator
                   write(*,*) "cell centre ", subcellCentre(thisOctal,subcell)
                   write(*,*) "direction ",tVec
                endif
                tempOctal => thisOctal
                tempSubcell = 1
                call findSubcellLocalLevel(locator, tempOctal, tempSubcell, nDepth)
                tempOctal%feederCell(tempsubcell) = .true.
             endif


           endif
        enddo
     endif
666 continue
  end subroutine setupGhostsLevel

  integer function getBoundary(boundary) result(i)
    use inputs_mod, only : xplusbound, yplusbound, zplusbound, &
         xminusbound, yminusbound, zminusbound 
    character(len=*) :: boundary

    select case (boundary)
       case("xminus")
          i = xminusbound
       case("yminus")
          i = yminusbound
       case("zminus")
          i = zminusbound
       case("xplus")
          i = xplusbound
       case("yplus")
          i = yplusbound
       case("zplus")
          i = zplusbound
       case DEFAULT
          call writeFatal("Boundary not recognised in get boundary: "//trim(boundary))
          stop
     end select
   end function getBoundary

  subroutine refineGridGeneric(grid, tol)
    use inputs_mod, only : minDepthAMR, maxDepthAMR
    use mpi
    type(GRIDTYPE) :: grid
    integer :: iThread
    logical :: globalConverged(64), tConverged(64)
    integer :: ierr
    real(double) :: tol
    globalConverged = .false.
    if (myrankGlobal == 0) goto 666
    if (minDepthAMR == maxDepthAMR) goto 666 ! fixed grid

    do
       call setAllUnchanged(grid%octreeRoot)
       globalConverged(myRankGlobal) = .true.
       do iThread = 1, nThreadsGlobal-1
          if (myrankGlobal /= iThread) then 
             call hydroValuesServer(grid, iThread)
          else
             call refineGridGeneric2(grid%octreeRoot, grid, globalConverged(myRankGlobal), tol, inheritval=.false.)
             call shutdownServers()
          endif
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       enddo
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       call MPI_ALLREDUCE(globalConverged, tConverged, nThreadsGlobal-1, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
       if (ALL(tConverged(1:nthreadsGlobal-1))) exit
    enddo


    666 continue
  end subroutine refineGridGeneric


  recursive subroutine refineGridGeneric2(thisOctal, grid, converged, limit, inheritval)
    use inputs_mod, only : maxDepthAMR, photoionization, refineOnMass, refineOnTemperature
    use inputs_mod, only : refineonionization, massTol
    use mpi
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
    real(double) :: limit
    integer :: myRank, ierr
    logical :: refineOnGradient
    real(double) :: rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z
    integer :: nd

    converged = .true.
    converged_tmp=.true.
    
    refineOnGradient = .not.photoionization      

   if(refineonionization .and. .not. photoionization) then
      refineOnIonization = .false.
   end if


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call refineGridGeneric2(child, grid, converged, limit, inheritval)
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

                call getHydroValues(grid, locator, nd, rho, rhoe, rhou, rhov, rhow, energy, phi,x,y,z)
               
                split = .false.

                grad = abs((thisOctal%rho(subcell)-rho) / &
                     thisOctal%rho(subcell))
                maxGradient = max(grad, maxGradient)
                if (grad > limit) then
                   split = .true.
                endif

                if(thisOctal%corner(subcell) .and. thisOCtal%nDepth < maxDepthAMR) split = .true.


!                thisSpeed = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 &
!                     + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)
!                speed = sqrt(rhou**2 + rhov**2 + rhow**2)/rho
!                if (thisSpeed > 0.d0) then
!                   grad = abs(thisSpeed-speed) / thisSpeed
!                   maxGradient = max(grad, maxGradient)
!                   if (grad > limit) then
!                      split = .true.
!                   endif
!                endif
                
                if (split) then
 
                   if ((thisOctal%nDepth < maxDepthAMR).and.(thisOctal%nDepth <= nd)) then
                      call addNewChildWithInterp(thisOctal, subcell, grid)
!                      print *, "split A ", thisOctal%nDepth
                      exit
                   endif
                   

                   !Thaw - I am going to assume that this is handled in evenupgridMPI and rm it for now
!                   if ((neighbourOctal%nDepth < maxDepthAMR) .and. &
!                        octalOnThread(neighbourOctal, neighbourSubcell, myrankglobal).and. &
!                        (neighbourOctal%nDepth < thisOctal%nDepth)) then
!                      call addNewChildWithInterp(neighbourOctal, neighboursubcell, grid)
!                      print *, "split B ", neighbourOctal%nDepth, nd, thisOctal%nDepth
!                      converged = .false.!
!                      exit
!                   endif
   
                   if(thisOctal%corner(subcell) .and. thisOctal%nDepth < maxDepthAMR) then
!                     print *, "split C ", thisOctal%nDepth
                      call addNewChildWithInterp(neighbourOctal, neighboursubcell, grid)
                      converged = .false.
                      exit
                   end if
                
                endif
             endif
          enddo
       endif
!

    if (converged.and.refineOnMass) then
       if (((thisOctal%rho(subcell)*1.d30*thisOctal%subcellSize**3) > massTol) &
            .and.(thisOctal%nDepth < maxDepthAMR))  then
          if (thisOctal%oneD) then
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=.false., interp=.false., amrHydroInterp = .true.)
          else
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=.true., interp=.false.)
          endif
          converged = .false.
!        print *, "split D ", thisOctal%nDepth
!        write(*,*) masstol, (thisOctal%rho(subcell)*1.d30*thisOctal%subcellSize**3)
          exit
       endif
    endif

   if(converged .and. refineOnTemperature)then
       do i = 1, nDir
          maxGradient = 1.d-30
          locator = subcellCentre(thisOctal, subcell) + &
          (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
          if (inOctal(grid%octreeRoot, locator)) then

             neighbourOctal => thisOctal

             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

!Thaw - modded to use specifiable limit 

             grad = abs((thisOctal%temperature(subcell) - neighbourOctal%temperature(neighbourSubcell)) / & 
             (thisOctal%temperature(subcell)))

             maxGradient = max(grad, maxGradient)

             if (octalOnThread(neighbourOctal, neighbourSubcell, myRank)) then    
                if((maxGradient > limit) .and. (thisOctal%nDepth < maxDepthAMR)) then                                        
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   converged = .false.
                   exit
                endif
             endif
          end if
       enddo
      endif



    if (converged.and.refineOnIonization) then

       do i = 1, nDir
          
          maxGradient = 1.d-30

          locator = subcellCentre(thisOctal, subcell) + &
               (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
          if (inOctal(grid%octreeRoot, locator)) then
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

!Thaw - modded to use specifiable limit
             grad = abs((thisOctal%ionFrac(subcell, 1) - neighbourOctal%ionFrac(neighbourSubcell, 1)) / &
                  (thisOctal%ionFrac(subcell, 1)))

             maxGradient = max(grad, maxGradient)

             if (octalOnThread(neighbourOctal, neighbourSubcell, myRank)) then
!                if ((thisOctal%ionFrac(subcell,1) > 0.9d0).and.(neighbourOctal%ionFrac(neighbourSubcell,1) < 0.1d0) .and. &
                if((maxGradient > limit) .and. (thisOctal%nDepth < maxDepthAMR)) then
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   converged = .false.
                   exit
                endif
                
!                if ((thisOctal%ionFrac(subcell,1) < 0.4d0).and.(neighbourOctal%ionFrac(neighbourSubcell,1) > 0.6d0) .and. &
!                     (neighbourOctal%nDepth < maxDepthAMR) ) then
!                   call addNewChildWithInterp(neighbourOctal, neighboursubcell, grid)
!                   converged = .false.
!                   exit
!                endif
             endif

          endif
       enddo
    endif

    if (.not.converged) exit
 endif

end do

end subroutine refineGridGeneric2



  recursive subroutine refineEdges(thisOctal, grid,  converged, inherit)

    use inputs_mod, only : maxDepthAMR
    use mpi
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
             call addNewChildWithInterp(thisOctal, subcell, grid)
             converged = .false.
          endif
          if (.not.converged) exit
       endif
    end do

  end subroutine refineEdges

  recursive subroutine refineFeeders(thisOctal, grid,  converged, inherit)

    use inputs_mod, only : maxDepthAMR
    use mpi
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
             call addNewChildWithInterp(thisOctal, subcell, grid)
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

  subroutine locatorToNeighbourLevel(grid, thisOctal, subcell, direction, ncells, locator, nDepth)
    type(GRIDTYPE) :: grid
    integer :: nDepth
    type(OCTAL), pointer :: thisOctal, tempOctal
    integer :: subcell, tempSubcell, i, nCells
    type(VECTOR) :: direction, locator, rVec
    tempOctal => thisOctal
    tempSubcell = subcell
    do i = 1, nCells
       rVec = subcellCentre(tempOctal, tempSubcell) + &
            (tempOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell) * direction 
       call findSubcellLocalLevel(rVec, tempOctal, tempSubcell, nDepth)
    enddo
    locator = rVec
  end subroutine locatorToNeighbourLevel
       
  recursive subroutine unrefineCells(thisOctal, grid, nUnrefine, splitLimit)
    use inputs_mod, only : minDepthAMR
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    integer :: subcell, i, neighbourSubcell
    logical :: unrefine
    integer :: nc
    integer, intent(inout) :: nUnrefine
    real(double) :: rhow(8), rhov(8), rhou(8), rho(8), rhoe(8) , fac, limit
    real(double) :: cs(8), mass
    logical :: refinedLastTime, ghostCell, split
    real(double) :: r, maxGradient
    real(double) :: rho1, rhoe1, rhou1, rhov1, rhow1, grad, speed, thisSpeed, meanrho, meancs
    real(double) :: splitLimit, dv
    integer :: ndir
    logical :: debug, cornerCell
    type(VECTOR) :: dirvec(6), locator, centre

    debug = .false.
    limit  = 5.0d-3

    unrefine = .true.
    refinedLastTime = .false.
    ghostCell = .false.
    cornerCell = .false.
    nc = 0
    mass = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unrefineCells(child, grid, nUnrefine, splitLimit)
                exit
             end if
          end do
       else

          if (.not.octalonThread(thisOctal, subcell, myRankGlobal)) cycle
          if (thisOctal%refinedLastTime(subcell)) cycle
          nc = nc + 1
          rho(nc) = thisOctal%rho(subcell)
          rhoe(nc) = thisOctal%rhoe(subcell)
          rhou(nc) = thisOctal%rhou(subcell)
          rhov(nc) = thisOctal%rhov(subcell)
          rhow(nc) = thisOctal%rhow(subcell)
          cs(nc) = soundSpeed(thisOctal, subcell)

          if (thisOctal%threed) then
!             dv = cellVolume(thisOctal, Subcell) * 1.d30
             dv = thisOctal%subcellSize**3
          else if (thisOctal%twoD) then
             dv = thisOctal%subcellSize**2
          else if (thisOctal%oneD) then!
             dv = thisOctal%subcellSize
          endif

!          mass = mass +  thisOctal%rho(subcell)*cellVolume(thisOctal, subcell) * 1.d30
          mass = mass +  thisOctal%rho(subcell)*dv

          cs(nc) = soundSpeed(thisOctal, subcell)
          if (thisOctal%ghostCell(subcell)) ghostCell=.true.
          if (thisOctal%corner(subcell)) cornerCell=.true.
          
!          if (thisOctal%refinedLastTime(subcell)) refinedLastTime = .true.
       endif
    enddo


    unrefine = .false.


    if ((nc > 1).and..not.cornerCell) then

       unrefine = .true.
       meancs = SUM(cs(1:nc))/dble(nc)

       fac = MaxMinOverMean(rho,nc)
       if (fac > limit) then
          unrefine = .false.
       endif

       fac = maxMinOverMean(rhoe,nc)
       if (fac > limit) then
          unrefine = .false.
       endif

       fac = maxMinOverMean(rhou,nc)
       if ((fac > limit).and.(rhou(1)/meancs > 1.d-2)) then
          unrefine = .false.
       endif

       fac = maxMinOverMean(rhov,nc)
       if (fac > limit) then
          unrefine = .false.
       endif

       fac = maxMinOverMean(rhow,nc)
       if (fac > limit) then
          unrefine = .false.
       endif
      
       
    endif
       
    if (thisOctal%nDepth <= minDepthAMR) unrefine = .false.

    if (unrefine) then

       meanRho = sum(rho(1:nc))/dble(nc)
       thisSpeed = (sum(rhou(1:nc))/dble(nc))**2 + (sum(rhov(1:nc))/dble(nc))**2 + (sum(rhow(1:nc))/dble(nc))**2
       thisSpeed = sqrt(thisSpeed/meanrho**2)
       r = thisOctal%subcellSize + 0.01d0*grid%halfSmallestSubcell
       centre = thisOctal%centre
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

!       unrefine = .true.
       do i = 1, nDir
          maxGradient = 1.d-30
          locator = centre + r*dirVec(i)
          if (inOctal(grid%octreeRoot, locator)) then
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             
             if (octalOnThread(neighbourOctal, neighbourSubcell, myrankGlobal)) then

                rho1 = neighbourOctal%rho(neighbourSubcell)
                rhou1 = neighbourOctal%rhou(neighbourSubcell)
                rhov1 = neighbourOctal%rhov(neighbourSubcell)
                rhow1 = neighbourOctal%rhow(neighbourSubcell)
                rhoe1 = neighbourOctal%rhoe(neighbourSubcell)
                split = .false.
                
                grad = abs((meanrho-rho1) / meanRho)
                if (grad > splitLimit) then
                   split = .true.
                endif
                speed = sqrt(rhou1**2 + rhov1**2 + rhow1**2)/rho1
                if (thisSpeed > 0.d0) then
                   grad = abs(thisSpeed-speed) / thisSpeed
                   if ((grad > limit).and.(speed/meancs > 1d-2)) then
                      split = .true.
                   endif
                endif

                !THaw - added ".and. unrefine because previously it would only care about the last direction
                if (split) then
                   unrefine = .false.
                endif
             endif
          endif
       enddo
    endif


    if (unrefine.and.photoionization) then
       if ( (.not.all(thisOctal%ionFrac(1:thisOctal%maxChildren,2) < 0.01d0)).or. &
            (.not.all(thisOctal%ionFrac(1:thisOctal%maxChildren,2) < 0.99d0)) ) unrefine = .false.
    endif

    if ((thisOctal%nChildren == 0).and.unrefine) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
       nunrefine = nunrefine + 1
    endif
    
  end subroutine unrefineCells

  real(double) function maxminOverMean(vals, n)
    real(double) :: vals(:), mean
    integer :: n

    mean = SUM(vals(1:n))
    if (mean /= 0.d0) then
       maxminoverMean = abs((MAXVAL(vals(1:n)) - MINVAL(vals(1:n)))/mean)
    else
       maxminOverMean = 0.d0
    endif
  end function maxminOverMean
  
  recursive subroutine unrefineCellsPhotoion(thisOctal, grid)
    use inputs_mod, only : minDepthAMR
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
    use mpi  
    use inputs_mod, only : minDepthAMR, maxDepthAMR

    type(GRIDTYPE) :: grid
    logical :: inheritFlag
    integer, optional :: dumpfiles
    integer :: myRank, nThreads, thisThread
    integer :: ierr
    logical :: globalConverged(64), localChanged(64), tConverged(64)
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
    logical :: allThreadsConverged
    character(len=30) :: vtkFilename

!    character(len=20) :: plotfile
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)
    ! first even up the local part of the grid
    if (myrankGlobal == 0) goto 666

    globalConverged = .false.
    iter = 0

    if (PRESENT(dumpfiles)) then
       write(vtkfilename,'(a,i4.4,a,i4.4,a)') "start",dumpfiles,"_",iter,".vtk"
       call writeVtkFile(grid, vtkFilename)
    endif

    call unsetGhosts(grid%octreeRoot)
    call setupEdges(grid%octreeRoot, grid)
    call setupGhosts(grid%octreeRoot, grid)

    if (minDepthAMR == maxDepthAMR) goto 666 ! fixed grid

    allThreadsConverged = .false.
    do while(.not.allThreadsConverged)
       call setAllUnchanged(grid%octreeRoot)
       localChanged = .false.

       globalConverged = .false.
       do
          globalConverged(myRankGlobal) = .true.
          do iThread = 1, nThreadsGlobal-1
             if (myrankGlobal /= iThread) then
                call hydroValuesServer(grid, iThread)
             else
                call evenUpGrid(grid%octreeRoot, grid,  globalConverged(myrankGlobal), inherit=inheritFlag)
                call unsetGhosts(grid%octreeRoot)
                call setupEdges(grid%octreeRoot, grid)
                call setupGhosts(grid%octreeRoot, grid)!, flag=.true.)
                call shutdownServers()
             endif
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          enddo
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          call MPI_ALLREDUCE(globalConverged, tConverged, nThreadsGlobal-1, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
!          if (myrank == 1) write(*,*) "first evenup ",tConverged(1:nThreadsGlobal-1)
          if (ALL(tConverged(1:nthreadsGlobal-1))) exit
       enddo


       if (evenAcrossThreads) then
          nLocs = 0
          call locatorsToExternalCells(grid%octreeRoot, grid, nLocs(myRank), locs, thread, depth)
          call MPI_ALLREDUCE(nLocs, tempnLocs, nThreadsglobal-1, MPI_INTEGER, MPI_SUM,amrCOMMUNICATOR, ierr)
          nLocsGlobal = SUM(tempnLocs)

          do thisThread = 1, nThreads - 1
             if(thisThread == myRank) then
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

             else
                !                          write(*,*) "rank ",myRank," receiving message from ",ithread
                call mpi_recv(nSent, 1, MPI_INTEGER, thisThread, tag, MPI_COMM_WORLD, status, ierr)
                !                          write(*,*) "received ",nSent
                if (nSent(1) > 0) then
                   do i = 1, nSent(1)
                      call mpi_recv(tempsent, 4, MPI_DOUBLE_PRECISION, thisThread, tag, MPI_COMM_WORLD, status, ierr)
                      eLocs(i)%x = tempsent(1)
                      eLocs(i)%y = tempsent(2)
                      eLocs(i)%z = tempsent(3)
                      eDepth(i) = nint(tempsent(4))
                   enddo
                endif
                nExternalLocs = nSent(1)
             endif
          enddo

          do iThread = 1, nThreadsGlobal-1
             if (myrankGlobal /= iThread) then
                call hydroValuesServer(grid, iThread)
             else
                do i = 1, nExternalLocs
                   call splitAtLocator(grid, elocs(i), edepth(i), localChanged(myRank))
                enddo
                call shutdownServers()
             endif
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          enddo

          if (PRESENT(dumpfiles)) then
             write(vtkfilename,'(a,i4.4,a,i4.4,a)') "aftersplit",dumpfiles,"_",iter,".vtk"
             call writeVtkFile(grid, vtkFilename)
          endif

          globalConverged = .false.
          do
             globalConverged(myRankGlobal) = .true.
             do iThread = 1, nThreadsGlobal-1
                if (myrankGlobal /= iThread) then
                   call hydroValuesServer(grid, iThread)
                else
                   call evenUpGrid(grid%octreeRoot, grid,  globalConverged(myrankGlobal), inherit=inheritFlag)
                   call unsetGhosts(grid%octreeRoot)
                   call setupEdges(grid%octreeRoot, grid)
                   call setupGhosts(grid%octreeRoot, grid)!, flag=.true.)
                   call shutdownServers()
                endif
                call MPI_BARRIER(amrCOMMUNICATOR, ierr)
             enddo
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
             call MPI_ALLREDUCE(globalConverged, tConverged, nThreadsGlobal-1, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
!             if (myrank == 1) write(*,*) "second evenup ",tConverged(1:nThreadsGlobal-1)
             if (ALL(tConverged(1:nthreadsGlobal-1))) exit
          enddo
          if (PRESENT(dumpfiles)) then
             write(vtkfilename,'(a,i4.4,a,i4.4,a)') "aftereven",dumpfiles,"_",iter,".vtk"
             call writeVtkFile(grid, vtkFilename)
          endif

          call MPI_ALLREDUCE(localChanged, globalChanged, nThreadsGlobal-1, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)

!          if (myrank == 1) write(*,*) "globalChanged ",globalChanged(1:nThreadsGlobal-1)
          if (ANY(globalChanged(1:nThreadsGlobal-1))) then
             allThreadsConverged = .false.
          else
             allThreadsConverged = .true.
          endif
          iter = iter + 1
       else
          allThreadsConverged = .true.
       endif

    end do
666 continue
  end subroutine evenUpGridMPI

  recursive subroutine evenUpGrid(thisOctal, grid,  converged, inherit)
    use inputs_mod, only : maxDepthAMR
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    type(octal), pointer :: bOctal
    integer :: subcell, i, bSubcell
    logical :: converged, converged_tmp
    type(VECTOR) :: dirVec(6), centre, octVec, rVec, locator, bVec, direction, nVec
    type(vector) :: vecStore(3)
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    logical, optional :: inherit
    integer :: myRank, ierr
    real(double) :: rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z
    integer :: nd, nd2

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

          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell)

         
          !Thaw - force ghostcells to match their partner refinement
          if(thisOctal%ghostcell(subcell)) then
             locator = thisOctal%boundaryPartner(subcell)
             bOctal => thisOctal             
             rVec = subcellCentre(thisOctal, subcell)
             
             call findSubcellLocal(locator, bOctal, bSubcell)
             bVec = subcellCentre(bOctal, bSubcell)
             if(thisOctal%nDepth < bOctal%nDepth .and. thisOctal%nDepth < maxDepthAMR) then
                call addNewChildWithInterp(thisOctal, subcell, grid)
                rVec = subcellCentre(thisOctal, subcell)
                converged = .false.
                exit
             end if
             
             !THAW - need to provide a case for when boundary partner is not on thread (i.e. for periodics)
             if(octalOnThread(bOctal, bSubcell, myRank)) then
                if(thisOctal%nDepth > bOctal%nDepth .and. bOctal%nDepth < maxDepthAMR) then
                   call addNewChildWithInterp(bOctal, bsubcell, grid)
                   converged = .false.
                   exit
                end if
             end if
             

             !Periodic non-ghosts refinement (i.e. A and B, where |G|G|A|B|)
             if(thisOctal%boundaryCondition(subcell) == 2 .and. .not. thisOctal%edgecell(subcell)) then

                direction = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                
                !Make it a unit vector                                                                                                                                                                                                
                if(direction%x > 0.d0) direction%x = direction%x / direction%x
                if(direction%x < 0.d0) direction%x = -direction%x / direction%x
                if(direction%y > 0.d0) direction%y = direction%y / direction%y
                if(direction%y < 0.d0) direction%y = -direction%y / direction%y
                if(direction%z > 0.d0) direction%z = direction%z / direction%z
                if(direction%z < 0.d0) direction%z = -direction%z / direction%z
                
                locator = subcellCentre(thisOctal, subcell) + direction*(thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)

                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

                if(octalOnThread(neighbourOctal, neighbourSubcell, myRank)) then

                   if(neighbourOctal%nDepth > thisOctal%nDepth) then
                      call addNewChildWithInterp(thisOctal, subcell, grid)
                      converged = .false.
                      exit
                   else if (neighbourOctal%nDepth < thisOctal%nDepth) then
                      call addNewChildWithInterp(neighbourOctal, neighbourSubcell, grid)
                      converged = .false.
                      exit
                   end if
                end if
             end if
          end if


          !Ensure new ghosts plus partners are all equal refinement
          if(thisOctal%edgecell(subcell) .and. thisOctal%ghostcell(subcell)) then
             direction = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)


             if(abs(direction%x) > abs(direction%y)) then
                direction%y = 0.d0
                if(abs(direction%x) > abs(direction%z)) then !We are moving in the 悉 direction                                                                                                      
                   direction%z = 0.d0
                else                                         !Moving in the 您 direction                                                                                                             
                   direction%x = 0.d0
                   end if
             else
                direction%x = 0.d0
                if(abs(direction%y) > abs(direction%z)) then !We are moving in the 悠 direction                                                                                                      
                   direction%z = 0.d0
                   else                                         !Moving in the 您 direction                                                                                                             
                   direction%y = 0.d0
                   end if
             end if

             !Make it a unit vector
             if(direction%x > 0.d0) direction%x = direction%x / direction%x
             if(direction%x < 0.d0) direction%x = -direction%x / direction%x
             if(direction%y > 0.d0) direction%y = direction%y / direction%y
             if(direction%y < 0.d0) direction%y = -direction%y / direction%y
             if(direction%z > 0.d0) direction%z = direction%z / direction%z
             if(direction%z < 0.d0) direction%z = -direction%z / direction%z


             locator= subcellcentre(thisoctal, subcell) + direction * &
                  (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             
             neighbourOctal => thisOctal
             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
             
             if(neighbourOctal%nDepth > thisOctal%nDepth) then
                call addNewChildWithInterp(thisOctal, subcell, grid)
                converged = .false.
                exit
             end if

             !Thaw - special case for periodics                                                                                                                                                                                       
             if(thisOctal%boundaryCondition(subcell) == 2) then
                !direction will still be valid from above

                !Find the opposite edgecell
                locator = subcellCentre(thisOctal, subcell)+direction*((grid%octreeRoot%subcellSize*2.d0) - &
                     ((thisOctal%subcellSize/2.d0)+ grid%halfSmallestSubcell*0.01d0))

                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

                call getHydroValues(grid, locator, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z)
                if(nd > thisOctal%nDepth) then
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   converged = .false.
                   exit
                end if

             end if
                          
          end if
                       
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
                
                call getHydroValues(grid, octVec, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z)
                if(.not. thisOctal%ghostcell(subcell)) then
                   if(thisOctal%twoD) then
                      !                   if(.not. octalOnThread(neighbourOctal, neighbourSubcell, myrank)) then
                      rVec = subcellCentre(thisOctal, subcell)
                      bVec = subcellCentre(neighbourOctal, neighbourSubcell)
                      octVec = VECTOR(x, y, z) 
                      
                      !                      print *, "nd", nd, thisOctal%nDepth, myRank
                      if((nd - thisOctal%nDepth)==1) then
                         !THaw - need to ensure that both of the neighbours are of lower refinement(2D)
                         if(abs(dirVec(j)%z) == 1.d0) then
                            if(rVec%x > octVec%x) then                                  
!                               print *, "a1", octvec
                               octVec = octVec + dirVec(1)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
!                               print *, "a2", octvec
                            else
!                               print *, "b1", octvec
                               octVec = octVec + dirVec(2)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
!                               print *, "b2", octvec
                            end if
                            
                         else if (abs(dirVec(j)%x) == 1.d0) then
                            
                            if(rVec%z >octVec%z) then
!                               print *, "c1", octvec
                               octVec = octVec + dirVec(3)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
!                               print *, "c2", octvec
                            else
!                               print *, "d1", octvec
                               octVec = octVec + dirVec(4)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
!                               print *, "d2", octvec
                            end if
                            
                         else 
                            print *, "Unrecognized direction", dirVec(j)
                            stop
                         end if
                         
                         
                         call getHydroValues(grid, octVec, nd2, rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z)
                         
                         if(nd2 > nd) then
!                            print *, "nd1", nd
                            nd = nd2
!                            print *, "nd2", nd
                         end if
                      end if
                      !end if
                  
                   else if(thisOctal%threeD) then
                      !As with 2D, need to ensure that all cells in contact are checked across an mpi boundary
                      rVec = subcellCentre(thisOctal, subcell)
                      nVec = subcellCentre(neighbourOctal, neighbourSubcell)
                      octVec = VECTOR(x, y, z)
                      vecStore(1) = octVec
                      
!                      print *, "doing 3D evenup"

                      !Reminder of my labels
                      ! _______
                      !|i  |ii |   z        x        y
                      !|___|___|   ^        ^        ^
                      !|iii|iV |   |        |        |
                      !|___|___|   |____>y, |____>z, |____x
                      !                 

                      !If the result of initial check is inconclusive we need to check other cells of the neighbour
                      if((nd - thisOctal%nDepth)==1) then

                         !Work round the neighbouring cells clockwise and get their positions

                         !Note that this process may not actually reach each "quarter" first time
                         !If it encounters a cell of higher refinement it might get lost (but only within the neighbour cube)
                         !This is no problem since if it occurs it implies refinement is necessary 
                         
                         if(abs(dirVec(j)%x) == 1.d0) then
                            if(rVec%y > octVec%y) then
                               if(rVec%z > octVec%z) then
                                  !Have found (i)
                                  
                                  !goto (ii)
                                  vecStore(1) = vecStore(1) + dirVec(2)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iv)
                                  vecStore(2) = vecStore(1) + dirVec(6)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iii)
                                  vecStore(3) = vecStore(2) + dirVec(5)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)

                               else if(rVec%z < octVec%z) then
                                  !Have found (iii)

                                  !goto (i)                                                                                                                                                         
                                  vecStore(1) = vecStore(1) + dirVec(1)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(ii)                                                                                                                                                         
                                  vecStore(2) = vecStore(1) + dirVec(2)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iv)                                                                                         
                                  vecStore(3) = vecStore(2) + dirVec(6)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)

                               end if
                               
                            else if (rVec%y < octVec%y) then
                               if(rVec%z > octVec%z) then
                                  !have found ii
                                  !goto (iv)
                                  vecStore(1) = vecStore(1) + dirVec(6)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iii)
                                  vecStore(2) = vecStore(1) + dirVec(5)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(i)
                                  vecStore(3) = vecStore(2) + dirVec(1)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)

                               else if(rVec%z < octVec%z) then
                                  !have found iv
                                  !goto (iii)
                                  vecStore(1) = vecStore(1) + dirVec(5)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(i)
                                  !vecStore(2) = vecStore(1) + dirVec(1)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  vecStore(2) = vecStore(1) + dirVec(6)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(ii)
                                  vecStore(3) = vecStore(2) + dirVec(2)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)

                               end if
                            end if

                            !vecstore check                                                                                                                            
                            do i = 1, 3
                               locator = vector(rVec%x, vecStore(i)%y, vecStore(i)%z)
                               if(.not. inSubcell(thisOctal,subcell, locator)) then
                                  write(*,*) "Screw up in partner checks x"
                                  stop
                               end if
                            end do

                                                     
                         else if(abs(dirVec(j)%y) == 1.d0) then !x-->y, y-->z, z-->x 
                            if(rVec%z > octVec%z) then
                               if((rVec%x > octVec%x)) then
                               !Have found (i)       

                                  !goto (ii)
                                  vecStore(1) = vecStore(1) + dirVec(1)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iv)
                                  vecStore(2) = vecStore(1) + dirVec(4)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iii)
                                  vecStore(3) = vecStore(2) + dirVec(6)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                                                                                                                             
                               else if(rVec%x < octVec%x) then
                                  !Have found (iii)
                                  !goto (i)
                                  vecStore(1) = vecStore(1) + dirVec(3)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(ii)
                                  vecStore(2) = vecStore(1) + dirVec(1)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iv)
                                  vecStore(3) = vecStore(2) + dirVec(4)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                               end if

                            else if (rVec%z < octVec%z) then
                               if(rVec%x > octVec%x) then
                                  !have found ii
                                  !goto (iv)
                                  vecStore(1) = vecStore(1) + dirVec(4)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iii)
                                  vecStore(2) = vecStore(1) + dirVec(6)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(i)
                                  vecStore(3) = vecStore(2) + dirVec(3)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                               else if(rVec%x < octVec%x) then
                                  !have found iv
                                  !goto (iii)
                                  vecStore(1) = vecStore(1) + dirVec(6)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(i)
                                  !vecStore(2) = vecStore(1) + dirVec(3)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  vecStore(2) = vecStore(1) + dirVec(4)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(ii)
                                  vecStore(3) = vecStore(2) + dirVec(1)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                               end if
                            end if

                            !vecstore check                                                                                                                            
                            do i = 1, 3
                               locator = vector(vecStore(i)%x, rVec%y, vecStore(i)%z)
                               if(.not. inSubcell(thisOctal,subcell, locator)) then
                                  write(*,*) "Screw up in partner checks y"
                                  stop
                               end if
                            end do


                         else if(abs(dirVec(j)%z) == 1.d0) then !x-->y, y-->z, z-->x 
                            if(rVec%x > octVec%x) then
                               if(rVec%y > octVec%y) then
                                  !Have found (i)
                                  !goto (ii) 
                                  vecStore(1) = vecStore(1) + dirVec(3)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iv)
                                  vecStore(2) = vecStore(1) + dirVec(5)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iii)
                                  vecStore(3) = vecStore(2) + dirVec(4)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)

                               else if(rVec%y < octVec%y) then
                                  !Have found (iii)
                                  !goto (i)
                                  vecStore(1) = vecStore(1) + dirVec(2)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(ii)
                                  vecStore(2) = vecStore(1) + dirVec(3)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iv)
                                  vecStore(3) = vecStore(2) + dirVec(5)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                               end if

                            else if (rVec%x < octVec%x) then
                               if(rVec%y > octVec%y) then
                                  !have found ii
                                  !goto (iv)
                                  vecStore(1) = vecStore(1) + dirVec(5)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(iii)
                                  vecStore(2) = vecStore(1) + dirVec(4)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(i)
                                  vecStore(3) = vecStore(2) + dirVec(2)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                               else if(rVec%y < octVec%y) then
                                  !have found iv
                                  !goto (iii)
                                  vecStore(1) = vecStore(1) + dirVec(4)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(i)
                                  !vecStore(2) = vecStore(1) + dirVec(2)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  vecStore(2) = vecStore(1) + dirVec(5)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                                  !goto(ii)
                                  vecStore(3) = vecStore(2) + dirVec(3)*(thisOctal%subcellSize/4.d0+0.01d0*grid%halfsmallestsubcell)
                               end if
                            end if

                            !vecstore check
                            do i = 1, 3
                               locator = vector(vecStore(i)%x, vecStore(i)%y, rVec%z)
                               if(.not. inSubcell(thisOctal,subcell, locator)) then
                                  write(*,*) "Screw up in partner checks z"
                                  stop
                               end if
                            end do


                         else
                            print *, "Direction Error In EvenUpGrid (3D)"
                            stop
                                                      
                         end if

                         !Then use the cells to find the highest refinement neighbour
                         do i = 1, 3
                            call getHydroValues(grid, vecstore(i), nd2, rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z)
                            if(nd2 > nd) then
                               nd = nd2
                            end if
                         end do
                        
                      end if
                   end if
                end if
                
                !                   if (((neighbourOctal%nDepth-thisOctal%nDepth) > 1).and. (thisOCtal%ndepth < maxDepthAMR)) then
                if (((nd-thisOctal%nDepth) > 1).and. (thisOctal%nDepth < maxDepthAMR)) then
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   converged = .false.
                   exit
                endif
                
                if (octalOnThread(neighbourOctal, neighbourSubcell, myrank)) then         
                   !                      if (((thisOctal%nDepth-neighbourOctal%nDepth) > 1).and.(neighbourOctal%ndepth < maxDepthAMR)) then
                   if (((thisOctal%nDepth-neighbourOctal%nDepth) > 1).and.(neighbourOctal%ndepth < maxDepthAMR)) then
                      call addNewChildWithInterp(neighbourOctal, neighboursubcell, grid)
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
  
  subroutine splitAtLocator(grid, locator, depth,  localchanged)
    use inputs_mod, only :  maxDepthAMR
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    type(VECTOR) :: locator
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
       call addNewChildWithInterp(thisOctal,subcell,grid)
       localChanged = .true.
!       write(*,*) "split at locator"
    endif

  end subroutine splitAtLocator

  recursive subroutine locatorsToExternalCells(thisOctal, grid, nLocs, loc, thread, depth)

    use mpi
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

    use mpi
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
    use mpi
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
    use mpi
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
    use mpi
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6)
    integer :: n, ndir, nd
    real(double) :: x1, x2, g(6)
    real(double) :: deltaT, fracChange, gGrav, newPhi, frac !, d2phidx2(3), sumd2phidx2
    real(double) :: xnext
    gGrav = bigG  * lengthToCodeUnits**3 / (massToCodeUnits * timeToCodeUnits**2)

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



          if (.not.thisOctal%ghostCell(subcell)) then 


             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + dir(n) * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, dir(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
                
                x1 = subcellCentre(thisOctal, subcell) .dot. dir(n)
                x2 = subcellCentre(neighbourOctal, neighboursubcell) .dot. dir(n)

!                g(n) =   (phi - thisOctal%phi_i(subcell))/(x2 - x1)
!                g(n) =   (phi - thisOctal%phi_i(subcell))/thisOctal%subcellSize !!!!!!!!!!!!!!!!!!!!!!!!!!

!                if (slow) then
                   g(n) = phigas
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
             if (thisOctal%phi_gas(subcell) /= 0.d0) then
                frac = abs((newPhi - thisOctal%phi_gas(subcell))/thisOctal%phi_gas(subcell))
                fracChange = max(frac, fracChange)
             else
                fracChange = 1.d30
             endif
             thisOctal%phi_gas(subcell) = newPhi
             
          endif
       endif
    enddo
  end subroutine gSweep

  recursive subroutine gSweep2(thisOctal, grid, deltaT, fracChange)
    use mpi
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6), probe(6)
    integer :: n, ndir
    real(double) ::  g(6), dx
    real(double) :: deltaT, fracChange, gGrav, newPhi, newerPhi, frac, d2phidx2(3), sumd2phidx2
    integer :: nd
    real(double) :: xnext, oldphi
    real(double), parameter :: SOR = 1.2d0
    gGrav = bigG * lengthToCodeUnits**3 / (massToCodeUnits * timeToCodeUnits**2)

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



          if (.not.thisOctal%ghostCell(subcell)) then


             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + probe(n) * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
                neighbourOctal => thisOctal
                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
                
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


                g(n) =   (phigas - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))

                
             enddo

             if (thisOctal%twoD) then
                d2phidx2(1) = (g(1) - g(2)) / (returnCodeUnitLength(thisOctal%subcellSize*gridDistanceScale))
                d2phidx2(2) = (g(3) - g(4)) / (returnCodeUnitLength(thisOctal%subcellSize*gridDistanceScale))
                sumd2phidx2 = SUM(d2phidx2(1:2))
             else
                d2phidx2(1) = (g(1) - g(2)) / (returnCodeUnitLength(thisOctal%subcellSize*gridDistanceScale))
                d2phidx2(2) = (g(3) - g(4)) / (returnCodeUnitLength(thisOctal%subcellSize*gridDistanceScale))
                d2phidx2(3) = (g(5) - g(6)) / (returnCodeUnitLength(thisOctal%subcellSize*gridDistanceScale))
                sumd2phidx2 = SUM(d2phidx2(1:3))
             endif
             deltaT = (1.d0/6.d0)*(returnCodeUnitLength(thisOctal%subcellSize*gridDistanceScale))**2
!             deltaT = deltaT * timeToCodeUnits

             oldPhi = thisOctal%phi_gas(subcell)
             newerPhi = thisOctal%phi_gas(subcell) + (deltaT * sumd2phidx2 - fourPi * gGrav * thisOctal%rho(subcell) * deltaT) 

             newPhi = (1.d0-SOR)*oldPhi + SOR*newerPhi

!             if (myrankglobal == 1) write(*,*) newphi,thisOctal%phi_i(subcell),deltaT, sumd2phidx2, thisOctal%rho(subcell)
             if (thisOctal%phi_gas(subcell) /= 0.d0) then
                frac = abs((newPhi - thisOctal%phi_gas(subcell))/thisOctal%phi_gas(subcell))
                fracChange = max(frac, fracChange)
             else
                fracChange = 1.d30
             endif
             if (.not.associated(thisOctal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))
             thisOctal%chiline(subcell) = frac

             thisOctal%phi_gas(subcell) = newPhi
             
          endif
       endif
    enddo
  end subroutine gSweep2

  recursive subroutine gSweepLevel(thisOctal, grid, deltaT, fracChange, ghostFracChange, nDepth)
    use mpi
    integer :: myRank, ierr, nd
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: nDepth
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6)
    integer :: n, ndir
    real(double) :: x1, x2, g(6)
    real(double) :: deltaT, fracChange, gGrav, newPhi, frac, ghostFracChange !, d2phidx2(3), sumd2phidx2
    real(double) :: sorFactor, deltaPhi
    real(double) :: xnext
    sorFactor = 1.2d0

    gGrav = bigG * lengthToCodeUnits**3 / (massToCodeUnits * timeToCodeUnits**2)

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call gsweepLevel(child, grid, deltaT, fracChange, ghostFracChange, ndepth)
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
          


          if (.not.thisOctal%ghostCell(subcell)) then 


             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + dir(n) * (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell)
                neighbourOctal => thisOctal
                call findSubcellLocalLevel(locator, neighbourOctal, neighbourSubcell, nDepth)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, dir(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext)
                
                x1 = subcellCentre(thisOctal, subcell) .dot. dir(n)
                x2 = subcellCentre(neighbourOctal, neighboursubcell) .dot. dir(n)

!                g(n) =   (phi - thisOctal%phi_i(subcell))/(x2 - x1)
!                g(n) =   (phi - thisOctal%phi_i(subcell))/thisOctal%subcellSize !!!!!!!!!!!!!!!!!!!!!!!!!!

!                if (slow) then
                   g(n) = phigas
!                endif

                enddo

                if (thisOctal%twoD) then
                   newphi = 0.25d0*(SUM(g(1:4))) - fourpi*gGrav*thisOctal%rho(subcell)*deltaT
                else
                   newphi = 0.16666666666667d0*(SUM(g(1:6))) - fourpi*gGrav*thisOctal%rho(subcell)*deltaT
                endif

                if (thisOctal%ghostCell(subcell)) then
                   ghostFracChange = max(fracChange, ghostFracChange)
                endif
                deltaPhi = newPhi - thisOctal%phi_gas(subcell)
                newPhi = thisOctal%phi_gas(subcell) + sorFactor * deltaPhi

                if (thisOctal%phi_gas(subcell) /= 0.d0) then
                   frac = abs((newPhi - thisOctal%phi_gas(subcell))/thisOctal%phi_gas(subcell))
                   fracChange = max(frac, fracChange)
                else
                   fracChange = 1.d30
                endif
                thisOCtal%phi_gas(subcell) = newPhi
             
          endif
       enddo
    endif
  end subroutine gSweepLevel
  
  subroutine selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid)
    use inputs_mod, only :  maxDepthAMR
    use mpi
    type(gridtype) :: grid
    logical, optional :: multigrid
    integer, parameter :: maxThreads = 100
    integer :: iDepth
    integer :: nPairs, thread1(:), thread2(:), nBound(:), group(:), nGroup
    real(double) :: fracChange(maxthreads), ghostFracChange(maxthreads), tempFracChange(maxthreads), deltaT, dx
    integer :: nHydrothreads
    real(double), parameter :: tol = 1.d-4,  tol2 = 1.d-4
    real(double) :: thisFrac
    integer :: it, ierr, i, j
!    character(len=30) :: plotfile
    nHydroThreads = nThreadsGlobal - 1

!    if (myrankglobal == 1) call tune(6,"Complete self gravity")

!    call saveGhostCellPhi(grid%octreeRoot)

    if (PRESENT(multigrid)) then


       call updateDensityTree(grid%octreeRoot)

       do iDepth = 3, maxDepthAmr-1

          call unsetGhosts(grid%octreeRoot)
          call setupEdgesLevel(grid%octreeRoot, grid, iDepth)
          call setupGhostsLevel(grid%octreeRoot, grid, iDepth)

          if (myrankglobal == 1) call tune(6,"Dirichlet boundary conditions")
          call applyDirichlet(grid, iDepth)
          if (myrankglobal == 1) call tune(6,"Dirichlet boundary conditions")

          dx = returnCodeUnitLength(gridDistancescale*grid%octreeRoot%subcellSize/dble(2.d0**(iDepth-1)))
          if (grid%octreeRoot%twoD) then
             deltaT =  (dx)**2 / 4.d0
          else
             deltaT =  (dx)**2 / 6.d0
          endif
!          deltaT = deltaT * timeToCodeUnits
          it = 0
          fracChange = 1.d30
          do while (ANY(fracChange(1:nHydrothreads) > tol))
             fracChange = 0.d0
             ghostFracChange = 0.d0
             it = it + 1
              
!             if (myrankglobal == 1) call tune(6,"Boundary exchange")
             call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
!             if (myrankglobal == 1) call tune(6,"Boundary exchange")
              
!             if (myrankglobal == 1) call tune(6,"Gsweep level")
            call gSweepLevel(grid%octreeRoot, grid, deltaT, fracChange(myRankGlobal), ghostFracChange(myRankGlobal),iDepth)
!             if (myrankglobal == 1) call tune(6,"Gsweep level")

!             if (myrankGlobal == 1) write(*,*) "Multigrid iteration ",it, " maximum fractional change ", MAXVAL(fracChange(1:nHydroThreads))
!             if (myrankglobal == 1) call tune(6,"Periodic boundary")

            !Thaw - trying to remove expansion by removing periodic gravity
!             call periodBoundaryLevel(grid, iDepth, justGrav = .true.)
             !             call imposeboundary(grid%octreeroot)
!             call transferTempStorageLevel(grid%octreeRoot, iDepth, justGrav = .true.)
!             if (myrankglobal == 1) call tune(6,"Periodic boundary")

             call MPI_ALLREDUCE(fracChange, tempFracChange, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR, ierr)
             fracChange = tempFracChange

             do i = iDepth, maxDepthAMR-1
                call updatePhiTree(grid%octreeRoot, i)
             enddo


             if (myrankGlobal == 1) write(*,*) it,MAXVAL(fracChange(1:nHydroThreads))
          enddo
          if (myRankGlobal == 1) write(*,*) "Gsweep of depth ", iDepth, " done in ", it, " iterations"

          
          call updatePhiTree(grid%octreeRoot, iDepth)

       enddo


    endif

    call unsetGhosts(grid%octreeRoot)
    call setupEdges(grid%octreeRoot, grid)
    call setupGhosts(grid%octreeRoot, grid)

    if (myrankglobal == 1) call tune(6,"Dirichlet boundary conditions")
    call applyDirichlet(grid)
    if (myrankglobal == 1) call tune(6,"Dirichlet boundary conditions")

!    call reapplyGhostCellPhi(grid%octreeRoot)

    if (grid%octreeRoot%twoD) then
       deltaT =  (2.d0*returnCodeUnitLength(gridDistancescale*grid%halfSmallestSubcell))**2 / 4.d0
    else
       deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*grid%halfSmallestSubcell))**2 / 6.d0
    endif
!    deltaT = deltaT * timeToCodeUnits

    fracChange = 1.d30
    it =0 
    do while (ANY(fracChange(1:nHydrothreads) > tol2))
       fracChange = 0.d0
       it = it + 1
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!       thisFrac = 1.d30
!       do while (thisFrac > 1.d-6)
!          thisFrac = 0.d0
       do j = 1,10
          call gSweep2(grid%octreeRoot, grid, deltaT, thisFrac)
       enddo

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       fracChange = 0.d0
       call gSweep2(grid%octreeRoot, grid, deltaT, fracChange(myRankGlobal))
       call MPI_ALLREDUCE(fracChange, tempFracChange, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR, ierr)
       fracChange = tempFracChange


       !       write(plotfile,'(a,i4.4,a)') "grav",it,".png/png"
!           if (myrankglobal == 1)   write(*,*) it,MAXVAL(fracChange(1:nHydroThreads))

       if (myrankGlobal == 1) write(*,*) "Full grid iteration ",it, " maximum fractional change ", &
            MAXVAL(fracChange(1:nHydroThreads))


!             write(plotfile,'(a,i4.4,a)') "grav",it,".vtk"
!             call writeVtkFile(grid, plotfile, &
!                  valueTypeString=(/"chiline"/))

!       if (writeoutput) write(*,*) "frac change ",maxval(fracChange(1:nHydroThreads)),tol2
    enddo
    if (myRankGlobal == 1) write(*,*) "Gravity solver completed after: ",it, " iterations"


!    if (myrankglobal == 1) call tune(6,"Complete self gravity")


  end subroutine selfGrav

  real(double) function getPressure(thisOctal, subcell)
#ifdef PHOTOION
    use inputs_mod, only : photoionPhysics
    use ion_mod, only : nGlobalIon, globalIonArray, returnMu
#endif
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: eKinetic, eThermal, K, u2, eTot
    real(double), parameter :: gamma2 = 1.4d0, rhoCrit = 1.d-14
    real(double) :: mu, rhoPhys

    mu = 0.d0

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
             write(*,*) "eThermal problem: ",ethermal
             ethermal = TINY(ethermal)
          endif

          getPressure =  (thisOctal%gamma(subcell) - 1.d0) * thisOctal%rho(subcell) * eThermal
!          write(*,*) "gamma, ethermal",thisOctal%gamma(subcell), ethermal
!          write(*,*) "rhou,rhow ", thisOCtal%rhou(subcell), thisOctal%rhow(subcell)
!          write(*,*) "etot, ekinetic, eThermal ",etot, ekinetic, ethermal
       case(1) ! isothermal
#ifdef PHOTOION
          eThermal = thisOctal%rhoe(subcell) / thisOctal%rho(subcell)

          if (photoionPhysics) then
             mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
          else
             mu = 2.33d0
          endif

          getPressure =  thisOctal%rho(subcell)
          getPressure =  getpressure/((mu*mHydrogen))
          getPressure =  getpressure*kerg*thisOctal%temperature(subcell)

!          getPressure =  (thisOctal%rho(subcell)/(mu*mHydrogen))*kerg*thisOctal%temperature(subcell)

          !if(getPressure < 1.d-60) then
          !   print *, "LOW PRESSURE - halting"
         !    print *, "thisOctal%rho(subcell)", thisOctal%rho(subcell)
          !   print *, "mu ", mu
          !   print *, "thisOctal%temperature(subcell) ",thisOctal%temperature(subcell) 
          !   print *, "getPressure ", getPressure 
          !   print *, "mu*mHydrogen ", mu*mHydrogen!
          !   print *, "kerg*thisOctal%temperature(subcell) ", kerg*thisOctal%temperature(subcell)
          !    print *, "F", (thisOctal%rho(subcell)/mu*mHydrogen)*kerg*thisOctal%temperature(subcell)
          !   stop
          !end if
#else
          call writeFatal("hydrodynamics_mod: isothermal pressure needs photoionization")
          STOP
#endif

       case(2) !  equation of state from Bonnell 1994
          rhoPhys = returnPhysicalUnitDensity(thisOctal%rho(subcell))
          if (rhoPhys < rhoCrit) then
             K = 4.1317d8
             getPressure = K * rhoPhys
!             write(*,*) "low density pressure called ",getPressure
          else
             K = (4.1317d8*rhoCrit)/(rhoCrit**gamma2)
             getPressure = K * rhoPhys**gamma2
!             write(*,*) "high density pressure called ",getpressure
          endif
          getPressure = returnCodeUnitPressure(getPressure)
!          getpressure = getpressure * 1.d-5
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


  RECURSIVE SUBROUTINE getChildlessOctalArray(thisOctal,array,counter)
    ! returns an array of pointers to all of the subcells in the grid.
    ! NB because fortran cannot create arrays of pointers, the output
    !   array is actually of a derived type which *contains* the 
    !   pointer to an octal.
    ! counter should be set to 0 before this routine is called

    IMPLICIT NONE

    TYPE(octal), POINTER                            :: thisOctal
    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    INTEGER, INTENT(INOUT)                          :: counter 
  
    INTEGER              :: i
    TYPE(octal), POINTER :: child

    ! if this is the root of the tree, we initialize the counter
    IF (.NOT. ASSOCIATED(thisOctal%parent)) counter = 0

    if ((thisOctal%nChildren == 0).and.(octalOnThread(thisOctal,1,myRankGlobal))) then
       counter = counter + 1 
       array(counter)%content => thisOctal
       array(counter)%inUse = .TRUE. 
    endif
    
    IF ( thisOctal%nChildren > 0 ) THEN
      DO i = 1, thisOctal%nChildren, 1
        
        ! call this subroutine recursively on each of its children
        child => thisOctal%child(i)
        CALL getChildlessOctalArray(child,array,counter)
        
      END DO
    END IF

  END SUBROUTINE getChildlessOctalArray

  subroutine createChildlessOctalArray(grid)
    type(GRIDTYPE) :: grid
    integer :: nvoxels, noctals
    if (allocated(globalChildlessOctalArray)) deallocate(globalChildlessOctalArray)
    call countVoxels(grid%octreeRoot,nOctals,nVoxels)  
    allocate(globalChildlessOctalArray(1:nOctals))
    nGlobalChildlessOctals = 0
    call getChildlessOctalArray(grid%octreeRoot, globalChildlessOctalArray, nGlobalChildlessOctals)
  end subroutine createChildlessOctalArray


  real(double) function  legendre(n, x)
    integer :: n
    real(double) :: x

    select case(n)
       case(0)
          legendre = 1.d0
       case(1)
          legendre = x
       case(2)
          legendre = (3.d0*x**2 - 1.d0)/2.d0
       case(3)
          legendre = (5.d0*x**3 - 3.d0*x)/2.d0
       case(4)
          legendre = (35.d0*x**4 - 30.d0*x**2 +3)/8.d0
       case DEFAULT
          write(*,*) "legendre: n not found ", n
     end select
   end function legendre

   subroutine findCoM(grid, com)
     use mpi
     type(GRIDTYPE) :: grid
     type(VECTOR) :: com
     real(double) :: totalMass, temp(3), temp2(3)
     integer :: ierr

     call findMassOverAllThreads(grid, totalmass)
     
     com = VECTOR(0.d0, 0.d0, 0.d0)
     call sumCOM(grid%octreeRoot, com)

     temp(1) = com%x
     temp(2) = com%y
     temp(3) = com%z
     call MPI_ALLREDUCE(temp, temp2, 3, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
     com%x = temp2(1)
     com%y = temp2(2)
     com%z = temp2(3)
     com = com / totalmass
   end subroutine findCoM

   recursive subroutine sumCOM(thisOctal, com)
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com
     integer :: subcell, i
     do subcell = 1, thisoctal%maxchildren
        if (thisoctal%haschild(subcell)) then
           ! find the child
           do i = 1, thisoctal%nchildren, 1
              if (thisoctal%indexchild(i) == subcell) then
                 child => thisoctal%child(i)
                 call sumCom(child, com)
                 exit
              end if
           end do
        else
           com = com + subcellCentre(thisOctal,subcell) * &
                thisOctal%rho(subcell)*cellVolume(thisOctal,subcell)*1.d30
        endif
     enddo
   end subroutine sumCOM



   recursive subroutine multipoleExpansion(thisOctal, point, com, v)
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point, xVec,rVec
     real(double) :: v, r, x, cosTheta, dm
     integer :: subcell, i, ipole

     do subcell = 1, thisoctal%maxchildren
        if (thisoctal%haschild(subcell)) then
           ! find the child
           do i = 1, thisoctal%nchildren, 1
              if (thisoctal%indexchild(i) == subcell) then
                 child => thisoctal%child(i)
                 call multipoleExpansion(child, point, com, v)
                 exit
              end if
           end do
        else
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           if (.not.thisOctal%ghostCell(subcell)) then
              xVec = point - com
              rVec = subcellCentre(thisOctal, subcell) - com
              r = modulus(rVec)
              x = modulus(xVec)
              cosTheta = (xVec.dot.rVec) / (r*x)
              r = r * 1.d10
              x = x * 1.d10
              dm = thisOctal%rho(subcell) * cellVolume(thisOctal,subcell) * 1.d30
              do iPole = 0, 4
                 v = v + (-bigG/x)*(r/x)**ipole * legendre(ipole, cosTheta) * dm
              enddo
              
           endif
        endif
     enddo
   end subroutine multipoleExpansion

   recursive subroutine multipoleExpansionLevel(thisOctal, point, com, v, level)
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point, xVec,rVec
     integer :: level
     real(double) :: v, r, x, cosTheta, dm
     integer :: subcell, i, ipole

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < level)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call  multipoleExpansionLevel(child, point, com, v, level)
       end do
    else
       do subcell = 1, thisOctal%maxChildren
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           if (.not.thisOctal%ghostCell(subcell)) then
              xVec = point - com
              rVec = subcellCentre(thisOctal, subcell) - com
              r = modulus(rVec)
              x = modulus(xVec)
              cosTheta = (xVec.dot.rVec) / (r*x)
              dm = thisOctal%rho(subcell) * cellVolume(thisOctal,subcell) * 1.d30
              do iPole = 0, 4
                 v = v + (-bigG/(x*1.d10))*(r/x)**ipole * legendre(ipole, cosTheta) * dm
              enddo
              
           endif
        enddo
     endif
   end subroutine multipoleExpansionLevel

   subroutine applyDirichlet(grid,level)
     use mpi
     type(GRIDTYPE) :: grid
     type(VECTOR) :: com
     integer, optional :: level
     real(double) :: temp(3)
     integer :: iThread
     integer :: tag
     logical :: stillReceiving
     integer :: status(MPI_STATUS_SIZE)
     integer :: ierr, j
     real(double) :: v

     type(VECTOR) :: point
     tag = 94
     call findCoM(grid, com)
     call updateDensityTree(grid%octreeRoot)

     do iThread = 1, nHydroThreadsGlobal
        if (myRankGlobal == iThread) then
           if (present(level)) then
              call recursApplyDirichletLevel(grid, grid%octreeRoot, com, level)
           else
              call recursApplyDirichlet(grid, grid%octreeRoot, com)
           endif
           do j = 1, nHydroThreadsGlobal
              if (j /= iThread) then
                 temp(1) = 1.d30
                 call mpi_send(temp, 3, MPI_DOUBLE_PRECISION, j, tag, MPI_COMM_WORLD, ierr)
              endif
           enddo
        else
           stillReceiving = .true.
           do while (stillReceiving)
              call mpi_recv(temp, 3, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, status, ierr)
              if (temp(1) > 1.d29) then
                 stillReceiving = .false.
                 exit
              else
                 point%x = temp(1)
                 point%y = temp(2)
                 point%z = temp(3)
                 v = 0.d0
!                 if (present(level)) then
                    call multipoleExpansionLevel(grid%octreeRoot, point, com, v, level=4)
!                 else
!                    call multipoleExpansion(grid%octreeRoot, point, com, v)
!                 endif
                 call mpi_send(v, 1, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
              endif
           enddo
        endif
        call mpi_barrier(amrCommunicator, ierr)
     enddo
   end subroutine applyDirichlet


   recursive subroutine recursApplyDirichlet(grid, thisOctal, com)
     use mpi
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point
     real(double) :: v, vgrid
     real(double) :: temp(3)
     integer :: subcell, i
     integer :: ithread
     integer :: tag, ierr
     integer :: status(MPI_STATUS_SIZE)

     tag = 94
     
     do subcell = 1, thisoctal%maxchildren
        if (thisoctal%haschild(subcell)) then
           ! find the child
           do i = 1, thisoctal%nchildren, 1
              if (thisoctal%indexchild(i) == subcell) then
                 child => thisoctal%child(i)
                 call recursApplyDirichlet(grid, child, com)
                 exit
              end if
           end do
        else
           if (.not.OctalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           if (thisOctal%ghostCell(subcell)) then
              point = subcellCentre(thisOctal, subcell)
              v = 0.d0
              call multipoleExpansionLevel(grid%OctreeRoot, point, com, v, level=4)

              do iThread = 1, nHydroThreadsGlobal

                 temp(1) = point%x
                 temp(2) = point%y
                 temp(3) = point%z
                 if (iThread /= myRankGlobal) then
                    call mpi_send(temp, 3, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
                 endif
              enddo
              do iThread = 1, nHydroThreadsGlobal
                 if (iThread /= myRankGlobal) then
                    call mpi_recv(vgrid, 1, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, status, ierr)
                    v = v + vgrid
                 endif
              enddo

              thisOctal%phi_gas(subcell) = v
           endif
        endif
     enddo
   end subroutine recursApplyDirichlet

   recursive subroutine recursApplyDirichletLevel(grid, thisOctal, com, level)
     use mpi
     type(GRIDTYPE) :: grid
     integer :: level
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point
     real(double) :: v, vgrid
     real(double) :: temp(3)
     integer :: subcell, i
     integer :: ithread
     integer :: tag, ierr
     integer :: status(MPI_STATUS_SIZE)

     tag = 94

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < level)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call recursApplyDirichletLevel(grid, child, com, level) 
       end do
    else
       do subcell = 1, thisOctal%maxChildren

           if (thisOctal%ghostCell(subcell)) then
              point = subcellCentre(thisOctal, subcell)
              v = 0.d0
              call multipoleExpansionLevel(grid%OctreeRoot, point, com, v, level=3)

              do iThread = 1, nHydroThreadsGlobal

                 temp(1) = point%x
                 temp(2) = point%y
                 temp(3) = point%z
                 if (iThread /= myRankGlobal) then
                    call mpi_send(temp, 3, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
                 endif
              enddo
              do iThread = 1, nHydroThreadsGlobal
                 if (iThread /= myRankGlobal) then
                    call mpi_recv(vgrid, 1, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, status, ierr)
                    v = v + vgrid
                 endif
              enddo
              thisOctal%phi_gas(subcell) = v
           endif
        enddo
     endif
   end subroutine recursApplyDirichletLevel


   subroutine mergeSinks(grid, source, nSource)
     type(GRIDTYPE) :: grid
     type(SOURCETYPE) :: source(:)
     type(SOURCETYPE), allocatable :: newSource(:)
     type(VECTOR) :: newVelocity, newPosition
     real(double) :: newMass, newAge
     integer :: nSource
     integer :: i, j
     integer :: nGroup, newNGroup
     integer, allocatable :: group(:)
     real(double) :: rAcc, sep
     integer :: m, n
     integer :: newNSource
     logical :: converged

     rAcc = 4.d0 * 2.d0 * grid%halfSmallestSubcell

     allocate(group(1:nSource))


     nGroup = 0
     group = 0 
     converged = .false.
     do while (.not.converged)
        converged = .true.
        do i = 1, nSource
           do j = 1, nSource
              if (i /= j) then
                 sep = modulus(source(i)%position - source(j)%position)
                 if (sep < rAcc) then
                    if ((group(i) == 0).and.(group(j) == 0)) then ! new group
                       nGroup = nGroup + 1
                       group(i) = nGroup
                       group(j) = nGroup
                       converged = .false.
                    else if ((group(i) == 0).and.(group(j) /= 0)) then ! connects to existing group
                       group(i) = group(j)
                       converged = .false.
                    else if ((group(i) /= 0).and.(group(j) == 0)) then ! connects to existing group
                       group(j) = group(i)
                       converged = .false.
                    else if (group(i) /= group(j)) then! now connect two existing groups
                       m = MIN(group(i), group(j))
                       n = MAX(group(i), group(j))
                       where (group == n) group = m
                       converged = .false.
                    endif
                 endif
              endif
           enddo
        end do
     enddo
     newNGroup = 0
     do j = 1, nGroup
        if (ANY(group(1:nSource) == j)) then
           newNGroup = newNGroup + 1
           where (group == j) group = newNGroup
        endif
     enddo
     nGroup = newNgroup
     if (nGroup == 0) goto 666 ! no groups

     if (writeoutput) then
        m = 0 
        do i = 0, nGroup

           n = 0
           do j = 1, nSource
              if (group(j) == i) n = n + 1
           enddo
           write(*,*) "Group ", i, " has ",n," sources"
           m = m + n
        enddo
        write(*,*) "total number of sources ",m
     endif

     newNSource = 0
     allocate(newSource(1:100))
     do j = 1, nSource
        if (group(j) == 0) then
           newNSource = newNSource + 1
           newSource(newNSource)%mass = source(j)%mass
           newSource(newNSource)%velocity = source(j)%velocity
           newSource(newNSource)%position = source(j)%position
           newSource(newNSource)%radius = source(j)%radius
           newSource(newNSource)%age = source(j)%age
        endif
     enddo

     do i = 1, nGroup

        newMass = 0.d0
        newVelocity = VECTOR(0.d0, 0.d0, 0.d0)
        newPosition = VECTOR(0.d0, 0.d0, 0.d0)
        newAge = 0.d0
        n = 0

        do j = 1, nSource
           if (group(j) == i) then
              newMass = newMass + source(j)%mass
              newVelocity = newVelocity + source(j)%mass * source(j)%velocity
              newPosition = newPosition + source(j)%position * source(j)%mass
              newAge = newAge + source(j)%age
              n = n + 1
           endif
        enddo
        newNSource = newNSource + 1
        newSource(newNSource)%mass = newMass
        newSource(newnSource)%velocity = (1.d0/newMass) * newVelocity 
        newSource(newnSource)%position = (1.d0/newMass) * newPosition
        newSource(newnSource)%radius = rSol/1.d10
        newSource(newnSource)%age = newAge / dble(n)
     enddo

        

     nSource = newnSource
     source(1:nSource)%mass = newSource(1:nSource)%mass
     source(1:nSource)%velocity = newSource(1:nSource)%velocity
     source(1:nSource)%position = newSource(1:nSource)%position
     source(1:nSource)%radius = newSource(1:nSource)%radius
     source(1:nSource)%age = newSource(1:nSource)%age

     do j = 1, nSource
        call emptySurface(source(j)%surface)
        call buildSphereNbody(source(j)%position, grid%halfSmallestSubcell, source(j)%surface, 20)
     enddo

     if (writeoutput) then
        write(*,*) "Number of sources after merges: ",newnSource
        call writesourcelist(source, nSource)
        call writeVtkFilenBody(nsource, source, "nbody_temp.vtk", grid)
     endif
     deallocate(newSource)
666 continue



     deallocate(group)
   end subroutine mergeSinks
   subroutine addSinks(grid, source, nSource)
     use mpi
     type(GRIDTYPE) :: grid
     type(SOURCETYPE) :: source(:)
     real(double) :: temp(8),rhomax
     integer :: nSource
     integer :: iThread
     integer :: j, tag, ierr
     integer :: status(MPI_STATUS_SIZE)
     logical :: stillReceiving

     tag = 77

     do iThread = 1, nHydroThreadsGlobal
        if (myrankGlobal == iThread) then
           rhomax = 0.d0
           call recursaddSinks(grid%octreeRoot, grid, source, nSource, rhomax)
           write(*,*) "Maximum ratio of local density:Jeans density ",rhomax
           do j = 1, nHydroThreadsGlobal
              if (j /= myRankGlobal) then
                 temp(1) = 1.d30
                 call mpi_send(temp, 8, MPI_DOUBLE_PRECISION, j, tag, MPI_COMM_WORLD, ierr)
              endif
           enddo
        else
           stillReceiving = .true.
           do while (stillReceiving)
              call mpi_recv(temp, 8, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, status, ierr)
              if (temp(1) > 1.d29) then
                 stillReceiving = .false.
                 exit
              else
                 nSource = nSource + 1
                 source(nsource)%position%x = temp(1) 
                 source(nsource)%position%y = temp(2) 
                 source(nsource)%position%z = temp(3) 
                 source(nsource)%mass  = temp(4) 
                 source(nsource)%velocity%x = temp(5) 
                 source(nsource)%velocity%y = temp(6) 
                 source(nsource)%velocity%z = temp(7) 
                 source(nsource)%radius = temp(8) 
                 call buildSphereNbody(source(nsource)%position, grid%halfSmallestSubcell, source(nsource)%surface, 20)
              endif
           enddo
        endif
     enddo
   end subroutine addSinks
          

  recursive subroutine recursaddSinks(thisOctal, grid, source, nSource, rhomax)
    use mpi
    type(OCTAL), pointer :: thisOctal, child
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    type(VECTOR) :: vel
    real(double) :: rhoJeans, bigJ, cs, e, rhomax
    real(double) :: temp(8)
    integer :: nSource
    integer :: i, subcell, ierr, ithread, tag
    tag = 77
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call recursaddSinks(child, grid, source, nSource, rhomax)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          bigJ = 0.2d0

          cs = soundSpeed(thisOctal, subcell)
          rhoJeans = bigJ**2 * pi * cs**2 / (bigG * returnCodeUnitLength(thisOctal%subcellSize*1.d10)**2) ! krumholz eq 6


          rhoMax = max(thisOctal%rho(subcell) / rhoJeans, rhoMax)
          
          if (thisOctal%rho(subcell) > rhoJeans) then

             write(*,*) "Source created"
             write(*,*) "local density ",thisOctal%rho(subcell)
             write(*,*) "local jeans density ", rhoJeans
             nSource = nSource + 1
             source(nSource)%position = subcellCentre(thisOctal, subcell)
             source(nsource)%mass = (thisOctal%rho(subcell) - rhoJeans*bigJ)*thisOctal%subcellSize**3*1.d30
             source(nsource)%velocity%x = thisOctal%rhou(subcell)/thisOctal%rho(subcell)
             source(nsource)%velocity%y = thisOctal%rhov(subcell)/thisOctal%rho(subcell)
             source(nsource)%velocity%z = thisOctal%rhow(subcell)/thisOctal%rho(subcell)             
             source(nsource)%radius = rsol/1.d10
             call buildSphereNbody(source(nsource)%position, grid%halfSmallestSubcell, source(nsource)%surface, 20)

             vel = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
                  thisOctal%rhow(subcell)/thisOctal%rho(subcell))
             e = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
             thisOctal%rho(subcell) = bigJ * rhoJeans
             thisOctal%rhou(subcell) = bigJ * rhoJeans * vel%x
             thisOctal%rhov(subcell) = bigJ * rhoJeans * vel%x
             thisOctal%rhow(subcell) = bigJ * rhoJeans * vel%x
             thisOctal%rhoe(subcell) = bigJ * rhoJeans * e

             temp(1) = source(nsource)%position%x
             temp(2) = source(nsource)%position%y
             temp(3) = source(nsource)%position%z
             temp(4) = source(nsource)%mass 
             temp(5) = source(nsource)%velocity%x
             temp(6) = source(nsource)%velocity%y
             temp(7) = source(nsource)%velocity%z
             temp(8) = source(nsource)%radius
             
             do iThread = 1, nHydroThreadsGlobal
                if (iThread /= myRankGlobal) then
                   call mpi_send(temp, 8, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
                endif
             enddo
             
          endif

       endif
    enddo
  end subroutine recursaddSinks

  function BondiHoyleRadius(source, thisOctal, subcell) result (rBH)

    real(double) :: rBH
    type(SOURCETYPE) :: source
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: vInfty, cInfty
    
    if (.not.inSubcell(thisOCtal, subcell, source%position)) then
       write(*,*) "bug in bondi-hoyle radius"
    endif
    
    vInfty = sqrt((thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
         thisOctal%rhow(subcell)**2) / thisOctal%rho(subcell)**2)
    cInfty = soundSpeed(thisOctal, subcell)
    rBH = bigG * source%mass / (vinfty**2 + cInfty**2)

  end function BondiHoyleRadius

  function mDotBondiHoyle(rhobar, source, thisOctal, subcell) result (mdot)

    type(SOURCETYPE) :: source
    type(OCTAL), pointer :: thisOctal
    real(double) :: vInfty, cInfty, rBH, mdot, rhoBar, thisAlpha
    integer :: subcell
    real(double), parameter :: lambda = 1.120d0
    real(double) :: rhoInfty

    if (.not.inSubcell(thisOCtal, subcell, source%position)) then
       write(*,*) "bug in bondi-hoyle radius"
    endif

    vInfty = sqrt((thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
         thisOctal%rhow(subcell)**2) / thisOctal%rho(subcell)**2)
    cInfty = soundSpeed(thisOctal, subcell)
    rBH = BondiHoyleRadius(source, thisOctal, subcell)
    write(*,*) "r/rBH ",1.2d0*thisOctal%subcellSize*gridDistanceScale/rBH
    write(*,*) "delta x/rBH ",thisOctal%subcellSize*gridDistanceScale/(bigG*source%mass/soundSpeed(thisOctal,subcell)**2)
    write(*,*) "c, v ",cInfty, vInfty, vinfty/cinfty
    thisAlpha = alpha(1.2d0*thisOctal%subcellSize*gridDistanceScale/rBH)

    rhoInfty = rhoBar / thisalpha

    mDot = fourPi * rhoInfty * rBH**2 * (lambda**2 * Cinfty**2 + vInfty**2)**0.5d0

  end function mDotBondiHoyle
  
  
  function accretionKernelRadius(source, thisOctal, subcell) result (rK)

    use inputs_mod, only : gridDistanceScale
    real(double) :: rk
    type(SOURCETYPE) :: source
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: deltaX
    real(double) :: rBH
    real(double) :: rAcc

    deltaX = thisOctal%subcellSize * gridDistanceScale
    rBH = bondiHoyleRadius(source, thisOctal, subcell)
    rAcc = 4.d0 * deltaX
    if (rBH < deltaX/4.d0) then
       rK = deltaX/4.d0
    else if (rBH <= rAcc/2.d0) then
       rK = rBH
    else
       rK = rAcc/2.d0
    endif
  end function accretionKernelRadius



  recursive subroutine correctForRotationRecur(thisOctal, source, mdot, n, timestep)
    use mpi
    type(SOURCETYPE) :: source
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: mDot, timestep, massCell
    integer :: n, thisn
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call correctForRotationRecur(child, source, mDot, n, timestep)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

          if (.not.associated(thisoctal%etaline)) then
             allocate(thisOctal%etaLine(1:thisOctal%maxChildren))
             thisOctal%etaLine = 0.d0
          endif
          thisOctal%etaLine(subcell) = 0.d0
          if (.not.inSubcell(thisOctal, subcell, source%position)) then

             call correctMdot(source, thisOctal, subcell, mdot, thisn)
             n = max(thisN, n)

             massCell = cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell)
             if (thisOctal%etaline(subcell) * timestep > 0.25d0*massCell) then
                thisOctal%etaLine(subcell) = 0.25d0 * massCell / timestep
             endif

          endif

       endif
    enddo
  end subroutine correctForRotationRecur

  recursive subroutine sumAccretionRate(thisOctal, mDot)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: mDot
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call sumAccretionRate(child, mdot)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

          mDot = mDot + thisOctal%etaline(subcell)

       endif
    enddo
  end subroutine sumAccretionRate




    


  subroutine correctMdot(source, thisOctal, subcell, mdot, n)
    
    type(SOURCETYPE) :: source
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: i, j, k, n
    real(double) :: mDot
    type(VECTOR) :: cellCentre, pointVelocity, pointPosition
    real(double) :: x, y, z
    real(double) :: esp, jsp, deltaX
    
    cellCentre = subcellCentre(thisOctal, subcell)
    deltaX = thisOctal%subcellSize * gridDistanceScale
    n = 0 
    do i = 1, 8
       do j = 1, 8
          do k = 1, 8
             x = thisOctal%subcellSize/16.d0 + thisOctal%subcellSize * dble(i-1)/8.d0 - cellCentre%x - thisOctal%subcellSize/2.d0
             y = thisOctal%subcellSize/16.d0 + thisOctal%subcellSize * dble(j-1)/8.d0 - cellCentre%y - thisOctal%subcellSize/2.d0
             z = thisOctal%subcellSize/16.d0 + thisOctal%subcellSize * dble(k-1)/8.d0 - cellCentre%z - thisOctal%subcellSize/2.d0
             pointPosition = (VECTOR(x, y, z) - source%position)*gridDistanceScale
             pointVelocity = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell) - source%velocity%x, &
                                    thisOctal%rhov(subcell)/thisOctal%rho(subcell) - source%velocity%y, &
                                    thisOctal%rhow(subcell)/thisOctal%rho(subcell) - source%velocity%z)
             jsp = modulus(pointPosition.cross.pointVelocity)
             esp = 0.5d0*modulus(pointVelocity)**2 - bigG*source%mass/modulus(pointPosition)
             if (esp > 0.d0) then
                rMin = 1.d30
             else
                rMin = -((bigG*source%mass)/(2.d0 * esp)) * (1.d0 - sqrt(1.d0 + (2.d0 * jsp * esp)/(bigG*source%mass)**2))
             endif
             if (rMin > deltaX/4.d0) n = n + 1
          enddo
       enddo
    enddo
    if (.not.associated(thisOctal%etaline)) then
       allocate(thisOctal%etaline(1:thisOctal%maxChildren))
       thisOctal%etaLine(subcell) = 0.d0
    endif
    thisOctal%etaline(subcell) =  mdot * (1.d0 - dble(n)/8.d0**3) * thisOctal%chiline(subcell)
  end subroutine correctMdot

  subroutine correctMomentumOfGas(source, thisOctal, subcell, timestep, deltaMom)
    type(SOURCETYPE) :: source
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: cellMass, radialMomentum, timeStep
    type(VECTOR) :: cellCentre, rVec, totalMomentum, deltaMom
    real(double) :: u2, ekinetic, ethermal, eTot, deltaM

    cellCentre = subcellCentre(thisOctal, subcell)
    rVec = source%position - cellCentre
    call normalize(rVec)


    cellMass = cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell)
    deltaM = thisOctal%etaline(subcell) * timestep

    totalMomentum = (cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell)) * &
         VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell)-source%velocity%x, &
                thisOctal%rhov(subcell)/thisOctal%rho(subcell)-source%velocity%y, &
                thisOctal%rhow(subcell)/thisOctal%rho(subcell)-source%velocity%z)
    radialMomentum = totalMomentum .dot. rVec


    totalMomentum = totalMomentum - radialMomentum * rVec

    deltaMom = deltaMom + radialMomentum * rVec
    deltaMom = VECTOR(0.d0, 0.d0, 0.d0)

    radialMomentum = radialMomentum * (1.d0 - deltaM/cellMass)

    totalMomentum = totalMomentum + radialMomentum * rVec

 
    if (thisOctal%threed) then
       u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
    else
       u2 = (thisOctal%rhou(subcell)**2 +  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
    endif
    eKinetic = u2 / 2.d0
    eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
    eThermal = eTot - eKinetic

    cellMass = cellMass - deltaM
    thisOctal%rho(subcell) = cellMass / (cellVolume(thisOctal, subcell)*1.d30)

!    thisOctal%rhou(subcell) = totalMomentum%x / (cellVolume(thisOctal, subcell)*1.d30)
!    thisOctal%rhov(subcell) = totalMomentum%y / (cellVolume(thisOctal, subcell)*1.d30)
!    thisOctal%rhow(subcell) = totalMomentum%z / (cellVolume(thisOctal, subcell)*1.d30)

!    if (thisOctal%threed) then
!       u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
!    else
!       u2 = (thisOctal%rhou(subcell)**2 +  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
!    endif
!    eKinetic = u2 / 2.d0
!
!    Etot = eThermal + eKinetic

!    thisOctal%rhoe(subcell) = eTot * thisOctal%rho(subcell)


  end subroutine correctMomentumOfGas


  subroutine correctMomentumOfSink(source, deltaMomentum)
    type(SOURCETYPE) :: source
    type(VECTOR) :: deltaMomentum, sourceMomentum

    sourceMomentum = source%mass * source%velocity
    sourceMomentum = sourceMomentum - deltaMomentum
    source%velocity = (1.d0/source%mass) * sourceMomentum
  end subroutine correctMomentumOfSink


  subroutine doAccretion(grid, sourceArray, nSource, timeStep)
    use mpi
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: sourceArray(:)
    integer :: nSource
    real(double) :: timeStep
    integer :: iSource
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: n
    real(double) :: rAcc, rK, rhobar, temp, thisMdot, localMdot, totalMdot, massCell
    integer :: ierr, itemp
    real(double) :: v(3), tempa(3)
    type(VECTOR) :: deltaMom

    do iSource = 1, nSource
       
       thisOctal => grid%octreeRoot
       call findSubcellLocal(sourceArray(iSource)%position, thisOctal, subcell)

       if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
          rAcc = 4.d0 * thisOctal%subcellSize * gridDistanceScale
          rK = accretionKernelRadius(sourceArray(iSource), thisOctal, subcell)
       endif
       call MPI_BCAST(rK, 1, MPI_DOUBLE_PRECISION, thisOctal%mpiThread(subcell)-1, amrCOMMUNICATOR, ierr)
       call MPI_BCAST(rAcc, 1, MPI_DOUBLE_PRECISION, thisOctal%mpiThread(subcell)-1, amrCOMMUNICATOR, ierr)
          
       if (myrankGlobal == 1) then
          write(*,*) "rk, racc ", rk, racc
       endif
       call calculateWeights(grid%octreeRoot, sourceArray(isource), rK, rAcc)
       call normalizeWeights(grid)
       rhoBar = 0.d0
       call calculateRhobar(grid%octreeRoot, rhoBar)

       call MPI_ALLREDUCE(rhoBar, temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
       rhoBar = temp

       if (myrankGlobal == 1) then
          write(*,*) "rho bar ", rhobar
       endif

       if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
          thisMdot = mDotBondiHoyle(rhoBar, sourceArray(iSource), thisOctal, subcell)
       endif
       call MPI_BCAST(thisMdot, 1, MPI_DOUBLE_PRECISION, thisOctal%mpiThread(subcell)-1, amrCOMMUNICATOR, ierr)

       if (myrankGlobal == 1) then
          write(*,*) "mdot initial ", thisMdot/msol*365.25d0*24.d0*3600.d0
       endif

       n = 0
       call correctForRotationRecur(grid%octreeRoot, sourceArray(isource), thismdot, n, timestep)
       call MPI_ALLREDUCE(n, itemp, 1, MPI_INTEGER, MPI_MAX, amrCommunicator, ierr)
       n = itemp

       if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
          massCell = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell) * 1.d30
          if (.not.associated(thisOctal%etaline)) then
             allocate(thisOctal%etaline(1:thisOctal%maxChildren))
             thisOctal%etaline = 0.d0
          endif
          write(*,*) "sink is on thread ",myrankglobal
          write(*,*) "using n for correction: ",n
          thisOctal%etaLine(subcell) = thismdot * thisOctal%chiline(subcell) * (1.d0 - dble(n)/8.d0**3)
          if (thisOctal%etaline(subcell) * timestep > 0.25d0*massCell) then
             thisOctal%etaLine(subcell) = 0.25d0 * massCell / timestep
             write(*,*) "local cell limited"
          endif
          localMdot = thisOctal%etaLine(subcell)
       endif
       call MPI_BCAST(localMdot, 1, MPI_DOUBLE_PRECISION, thisOctal%mpiThread(subcell)-1, amrCOMMUNICATOR, ierr)

       totalMdot = 0.d0
       call sumAccretionRate(grid%octreeRoot, totalMdot)
       call MPI_ALLREDUCE(totalMdot, temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
       totalMdot = temp
       if (myrankGlobal == 1) then
          write(*,*) "mdot after corrections ", totalMdot/msol*365.25d0*24.d0*3600.d0
       endif
       sourceArray(isource)%mdot = totalMdot
       deltaMom = VECTOR(0.d0, 0.d0, 0.d0)
       call correctMomenta(grid%octreeroot, sourceArray(isource), timestep, deltaMom)

       v(1) = deltaMom%x
       v(2) = deltaMom%y
       v(3) = deltaMom%z
       call MPI_ALLREDUCE(v, tempa, 3, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
       deltaMom%x = tempa(1)
       deltaMom%y = tempa(2)
       deltaMom%z = tempa(3)

       call correctMomentumOfSink(sourceArray(isource), deltaMom)

       sourceArray(isource)%mass = sourceArray(isource)%mass + totalMdot * timestep
       if (myrankGlobal == 1) write(*,*) "source mass ", sourceArray(iSource)%mass/msol

    enddo

  end subroutine doAccretion

  recursive subroutine correctMomenta(thisOctal, source, timestep, deltaMom)
    use mpi
    type(SOURCETYPE) :: source
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: timestep
    type(VECTOR) :: deltaMom
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call correctMomenta(child, source, timestep, deltaMom)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle


          call correctMomentumOfGas(source, thisOctal, subcell, timestep, deltaMom)

       endif
    enddo
  end subroutine correctMomenta

  recursive subroutine calculateWeights(thisOctal, source, rK, rAcc)
    use mpi
    type(SOURCETYPE) :: source
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: rK, rAcc, r
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call calculateWeights(child, source, rK, rAcc)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

          if (.not.associated(thisOctal%chiline)) then
             allocate(thisOctal%chiLine(1:thisOctal%maxChildren))
             thisOctal%chiLine = 0.d0
          endif
          r = modulus(subcellCentre(thisOctal,subcell) - source%position)*gridDistanceScale
          if (r < rAcc) then
             thisOctal%chiLine(subcell) = exp(-r**2/rK**2)
          else
             thisOctal%chiLine(subcell) = 0.d0
          endif
       endif
    enddo
  end subroutine calculateWeights

  subroutine normalizeWeights(grid)
    use mpi
    type(GRIDTYPE) :: grid
    real(double) :: sum(1), temp(1)
    integer :: ierr

    sum = 0.d0
    call sumWeights(grid%octreeRoot, sum(1))

    call MPI_ALLREDUCE(sum, temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
    sum(1) = temp(1)

    call normWeights(grid%octreeRoot, sum(1))

  end subroutine normalizeWeights

  recursive subroutine calculateRhobar(thisOctal, rhobar)
    use mpi
    real(double) :: rhobar
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call calculateRhobar(child, rhobar)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          
          rhobar = rhoBar + thisOctal%chiline(subcell) * thisOctal%rho(subcell)

       endif
    enddo
  end subroutine calculateRhobar

  recursive subroutine sumWeights(thisOctal, sum)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) ::  sum
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call sumWeights(child, sum)
                exit
             end if
          end do
       else

          
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          sum = sum + thisOctal%chiline(subcell)
       endif
    enddo
  end subroutine sumWeights

  recursive subroutine normWeights(thisOctal, sum)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: sum
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call normWeights(child, sum)
                exit
             end if
          end do
       else

          
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          thisOctal%chiline(subcell) = thisOctal%chiline(subcell) / sum
       endif
    enddo
  end subroutine normWeights



  function gfunc(x, gamma) 
    real(double) :: x, gamma, gfunc

    gfunc = (x**(4.d0*(gamma-1.d0)/(gamma+1.d0)))/(gamma-1.d0) + x**(-(5.d0-3.d0*gamma)/(gamma+1.d0))
  end function gfunc

  function ffunc(u, gamma) 
    real(double) :: u, gamma, ffunc

    ffunc = (u**(4.d0/(gamma+1.d0)))*(0.5d0 + (1.d0/(gamma-1.d0)) * (1.d0/u**2))
  end function ffunc

  

#endif

end module hydrodynamics_mod

#endif
