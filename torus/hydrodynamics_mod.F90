#ifdef HYDRO
! Hydrodynamics module added by tjh on 18th june 2007


module hydrodynamics_mod
#ifdef MPI

  use inputs_mod
  use utils_mod
  use viscosity_mod
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

    call allocateHydrodynamicsAttributes(grid%octreeRoot)

    currentlyDoingHydroStep = .true.
    if (grid%octreeRoot%twoD) then
       if (.not.cylindricalHydro) then
          call doHydrodynamics2d(grid)
       else
          call doHydrodynamics2dCylindrical(grid)
       endif
    else if (grid%octreeRoot%oneD) then
       call doHydrodynamics1d(grid)
    else if (grid%octreeRoot%threeD) then
       call doHydrodynamics3d(grid)
    endif
    currentlyDoingHydroStep = .false.
    
  end subroutine dohydrodynamics



  subroutine checkMaclaurinBenchmark(grid)
    use mpi
    type(GRIDTYPE) :: grid
    real(double) :: maxDev, temp
    integer :: ierr
    
    if (myrankGlobal == 0) goto 666
    
    maxDev = 0.d0
    call checkDeviations(grid%octreeRoot, maxDev)
    call MPI_ALLREDUCE(maxDev, temp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, amrCommunicator, ierr)
    if (writeOutput) then
       write(*,*) "Maximum deviation of potential ",maxDev
       if (maxDev > 1.d-2) then
          write(*,*) "TORUS Maclaurin spheroid gravity test failed"
       else
          write(*,*) "TORUS Maclaurin spheroid gravity test successful"
       endif
    endif
    666 continue
    call torus_mpi_barrier
  end subroutine checkMaclaurinBenchmark

  recursive subroutine checkDeviations(thisoctal, maxDev)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: maxDev
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call checkDeviations(child, maxDev)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
          maxDev = max(maxDev, abs((thisOctal%phi_gas(subcell) - &

               thisOctal%biasLine3d(subcell))/thisOctal%biasline3d(subcell)))
       endif
    enddo
  end subroutine checkDeviations

  recursive subroutine setOnAxisRhoU(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    logical :: onAxis
    type(VECTOR) :: cellCentre

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setOnAxisRhoU(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
          cellCentre = subcellCentre(thisOctal, subcell)
          onAxis = (cellCentre%x - thisOctal%subcellSize/2.d0-0.1d0*smallestCellSize < 0.d0)
          if (onAxis.and.thisOctal%rhou(subcell) < 0.) then
             thisOctal%rhou(subcell) = 0.d0
          endif
       endif
    enddo
  end subroutine setOnAxisRhoU


 recursive subroutine perturbIfront(thisOctal, grid)
#ifdef PHOTOION
   use inputs_mod, only : photoionPhysics
   use ion_mod, only : nGlobalIon, globalIonArray, returnMu
#endif
   type(octal), pointer :: thisOctal
   type(octal), pointer :: child
   type(octal), pointer :: neighbourOCtal
   type(gridtype) :: grid   
   integer :: i, subcell, neighbourSubcell
   type(VECTOR) :: direction, rVec, locator
   real(double) :: getpressure, cs, A, mu
   
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child                 
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call perturbIfront(child, grid)
                exit
             end if
          end do
       else
          !hard wired for benchmark geometry
          rVec = subcellCentre(thisOctal, subcell)
          direction = VECTOR(1.d0, 0.d0, 0.d0)
          locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
          neighbouroctal => thisoctal
          if(inOctal(grid%octreeRoot, locator)) then
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             !know that in this hardwired benchmark I do not have to seek across mpi boundary
             !this is pretty rubbish but for a test it isn't worth overhauling loads of the rest of the code.
!             if(thisOctal%pressure_i(subcell) /= neighbourOctal%pressure_i(subcell)) then
             if (photoionPhysics) then
                mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
             else
                mu = 2.33d0
             endif
             getPressure =  thisOctal%rho(subcell)
             getPressure =  getpressure/((mu*mHydrogen))
             getPressure =  getpressure*kerg*thisOctal%temperature(subcell)
             cs = sqrt(getPressure/thisOctal%rho(subcell))
             A = 0.75d0*cs
             if((thisOctal%ionFrac(subcell, 1)- neighbourOctal%ionFrac(subcell, 1)) > 0.2 ) then
                thisOctal%rhou(subcell) = abs(thisOctal%rho(subcell)*(A*(sin(twopi*(rVec%z+(grid%octreeRoot%subcellSize/2.d0))* &
                     (grid%octreeroot%subcellsize/8.d0)))))
             end if
         end if
       end if
    end do

 end subroutine perturbIfront



  recursive subroutine fluxlimiter(thisoctal)
    use inputs_mod, only : limiterType
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: a, b, dq, dx
    character(len=80) :: message
 

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

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

      !    if (.not.thisoctal%ghostcell(subcell)) then
             dq = thisoctal%q_i(subcell) - thisoctal%q_i_minus_1(subcell)
! When dq is very small there could be a floating point overflow in thisoctal%rlimit so check size of dq
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

!superbee flux limiter is default
             case("superbee")
                a = min(1.d0, 2.d0*thisoctal%rlimit(subcell))
                b = min(2.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit(subcell) = max(0.d0, a, b)
                

!Use with caution, this is not in the literature, developed by THAW
!Tries to improve flux limiter by accounting for non-constant grid spacing
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
                abs(thisoctal%rlimit(subcell))) / (1.d0 + abs(thisoctal%rlimit(subcell)))
               
!i.e no flux limiter
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

  recursive subroutine fluxlimiter2(thisoctal)
    use inputs_mod, only : limiterType
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, j
    real(double) :: a, b, dq, dx
    character(len=80) :: message
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call fluxlimiter2(child)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle


          do j = 1, 2
             dq = thisoctal%q_i(subcell) - thisoctal%q_amr_i_minus_1(subcell,j)
! When dq is very small there could be a floating point overflow in thisoctal%rlimit so check size of dq
             if (dq /= 0.d0) then
                if (thisoctal%u_amr_interface(subcell,j) .ge. 0.d0) then
                   thisoctal%rlimit(subcell) = (thisoctal%q_amr_i_minus_1(subcell,j) - thisoctal%q_i_minus_2(subcell)) / dq
                else
                   thisoctal%rlimit(subcell) = (thisoctal%q_amr_i_plus_1(subcell,j) - thisoctal%q_i(subcell)) / dq
                endif
             else
                thisoctal%rlimit(subcell) = 1.d0
             endif


             select case (limitertype)

!superbee flux limiter is default
             case("superbee")
                a = min(1.d0, 2.d0*thisoctal%rlimit(subcell))
                b = min(2.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit_amr(subcell,j) = max(0.d0, a, b)
                

!Use with caution, this is not in the literature, developed by THAW
!Tries to improve flux limiter by accounting for non-constant grid spacing
             case("modified_superbee")
                if (dq /= 0.d0) then
                   if (thisoctal%u_interface(subcell) .ge. 0.d0) then
                      thisoctal%rlimit(subcell) = (thisoctal%q_amr_i_minus_1(subcell,j) - thisoctal%q_i_minus_2(subcell)) / dq
                      dx = 0.5d0*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))

                      thisoctal%rlimit(subcell)=thisoctal%rlimit(subcell)*dx/(thisoctal%x_i_minus_1(subcell) - &
                           thisoctal%x_i_minus_2(subcell))

                   else
                      thisoctal%rlimit(subcell) = (thisoctal%q_amr_i_plus_1(subcell,j) - thisoctal%q_i(subcell)) / dq
                      dx = 0.5d0*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
                
                      thisoctal%rlimit(subcell)=thisoctal%rlimit(subcell)*dx/(thisoctal%x_i_plus_1(subcell) - &
                           thisoctal%x_i(subcell))
                   endif
                else
                   thisoctal%rlimit(subcell) = 1.d0
                endif

                a = min(1.d0, 2.d0*thisoctal%rlimit(subcell))
                b = min(2.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit_amr(subcell,j) = max(0.d0, a, b)

             case("minmod")
                
                a = min(1.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit_amr(subcell,j) = max(0.d0, a)
                
             case("MC")
                a = min((1.d0+thisoctal%rlimit(subcell))/2.d0, 2.d0, (2.d0*thisoctal%rlimit(subcell)))
                thisoctal%philimit_amr(subcell,j) = max(0.d0, a)

             case("fromm")
                thisoctal%philimit_amr(subcell,j) = 0.5d0*(1.d0 + thisoctal%rlimit(subcell))
                

             case("vanleer")
                !write(*,*) "vanleer"
                thisoctal%philimit_amr(subcell,j) = (thisoctal%rlimit(subcell) + &
                abs(thisoctal%rlimit(subcell))) / (1.d0 + abs(thisoctal%rlimit(subcell)))
               
!i.e no flux limiter
             case("donorcell")
                thisOctal%philimit(subcell) = 0.d0


             case default
                write(message,'(a,a)') "flux limiter not recognised: ", trim(limitertype)
                call writefatal(message)
                stop
             end select
!             if (thisoctal%ghostcell(subcell)) thisoctal%philimit(subcell) = 1.d0
!                        thisoctal%philimit(subcell) = 1.d0
          enddo
       endif

    enddo
  end subroutine fluxlimiter2

  recursive subroutine fluxlimiter_amr(thisoctal)
    use inputs_mod, only : limiterType
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, j
    real(double) :: a, b, dq, dx
    character(len=80) :: message
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call fluxlimiter_amr(child)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle


          do j = 1, 4
             dq = thisoctal%q_i(subcell) - thisoctal%q_amr_i_minus_1(subcell,j)
! When dq is very small there could be a floating point overflow in thisoctal%rlimit so check size of dq
             if (dq /= 0.d0) then
                if (thisoctal%u_amr_interface(subcell,j) .ge. 0.d0) then
                   thisoctal%rlimit(subcell) = (thisoctal%q_amr_i_minus_1(subcell,j) - thisoctal%q_i_minus_2(subcell)) / dq
                else
                   thisoctal%rlimit(subcell) = (thisoctal%q_amr_i_plus_1(subcell,j) - thisoctal%q_i(subcell)) / dq
                endif
             else
                thisoctal%rlimit(subcell) = 1.d0
             endif


             select case (limitertype)

!superbee flux limiter is default
             case("superbee")
                a = min(1.d0, 2.d0*thisoctal%rlimit(subcell))
                b = min(2.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit_amr(subcell,j) = max(0.d0, a, b)
                

!Use with caution, this is not in the literature, developed by THAW
!Tries to improve flux limiter by accounting for non-constant grid spacing
             case("modified_superbee")
                if (dq /= 0.d0) then
                   if (thisoctal%u_interface(subcell) .ge. 0.d0) then
                      thisoctal%rlimit(subcell) = (thisoctal%q_amr_i_minus_1(subcell,j) - thisoctal%q_i_minus_2(subcell)) / dq
                      dx = 0.5d0*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))

                      thisoctal%rlimit(subcell)=thisoctal%rlimit(subcell)*dx/(thisoctal%x_i_minus_1(subcell) - &
                           thisoctal%x_i_minus_2(subcell))

                   else
                      thisoctal%rlimit(subcell) = (thisoctal%q_amr_i_plus_1(subcell,j) - thisoctal%q_i(subcell)) / dq
                      dx = 0.5d0*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
                
                      thisoctal%rlimit(subcell)=thisoctal%rlimit(subcell)*dx/(thisoctal%x_i_plus_1(subcell) - &
                           thisoctal%x_i(subcell))
                   endif
                else
                   thisoctal%rlimit(subcell) = 1.d0
                endif

                a = min(1.d0, 2.d0*thisoctal%rlimit(subcell))
                b = min(2.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit_amr(subcell,j) = max(0.d0, a, b)

             case("minmod")
                
                a = min(1.d0, thisoctal%rlimit(subcell))
                thisoctal%philimit_amr(subcell,j) = max(0.d0, a)
                
             case("MC")
                a = min((1.d0+thisoctal%rlimit(subcell))/2.d0, 2.d0, (2.d0*thisoctal%rlimit(subcell)))
                thisoctal%philimit_amr(subcell,j) = max(0.d0, a)

             case("fromm")
                thisoctal%philimit_amr(subcell,j) = 0.5d0*(1.d0 + thisoctal%rlimit(subcell))
                

             case("vanleer")
                !write(*,*) "vanleer"
                thisoctal%philimit_amr(subcell,j) = (thisoctal%rlimit(subcell) + &
                abs(thisoctal%rlimit(subcell))) / (1.d0 + abs(thisoctal%rlimit(subcell)))
               
!i.e no flux limiter
             case("donorcell")
                thisOctal%philimit(subcell) = 0.d0


             case default
                write(message,'(a,a)') "flux limiter not recognised: ", trim(limitertype)
                call writefatal(message)
                stop
             end select
!             if (thisoctal%ghostcell(subcell)) thisoctal%philimit(subcell) = 1.d0
!                        thisoctal%philimit(subcell) = 1.d0
          enddo
       endif

    enddo
  end subroutine fluxlimiter_amr

  recursive subroutine imposeAzimuthalVelocity(thisOctal)
    use inputs_mod, only : sourceMass
    integer :: i
    type(octal), pointer :: thisOctal, child
    integer :: subcell
    real(double) :: vkep, r, cs, n, eta, thisVel, x
    type(vector) :: rvec

   do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call imposeAzimuthalVelocity(child)
                exit
             end if
          end do
       else
          n = 1.125
          rVec = subcellCentre(thisOctal,subcell)
          r = real(modulus(rVec*1.d10))
          x = rvec%x*1.d10
          vkep = sqrt(bigG*sourceMass(1)/(x))
          cs = soundSpeed(thisOctal, subcell)
          eta = n *(cs**2/vkep**2)
          thisVel = vkep*(1.d0-eta*(abs(rVec%z*1.d10)/x**3)**2)
          thisOctal%rhov(subcell) = thisOctal%rho(subcell) * thisVel*x
!          thisOctal%rhov(subcell) = thisOctal%rho(subcell) * thisVel
       endif
    enddo
  end subroutine imposeAzimuthalVelocity

  recursive subroutine imposeFontVelocity(thisOctal)
    use inputs_mod, only : sourceMass, amrgridsize, maxdepthamr!, amrgridcentrez
    integer :: i
    type(octal), pointer :: thisOctal, child
    integer :: subcell
    real(double) :: r, thisVel, dx
    type(vector) :: rvec

   do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call imposefontVelocity(child)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          r = abs(rVec%x)*1.d10
          dx = amrgridsize/(2.**maxdepthamr)
!          upperbound = amrgridcentrez - 0.5d0*amrgridsize + dx
          if(abs(rvec%z) < dx) then
          thisVel = sqrt(bigG*sourceMass(1)/(r))
          thisOctal%rhov(subcell) = thisOctal%rho(subcell) * thisVel*r
          endif
!          thisOctal%rhov(subcell) = thisOctal%rho(subcell) * thisVel
       endif
    enddo
  end subroutine imposeFontVelocity


  recursive subroutine imposeKeplerianVelocity(thisOctal)
!    use inputs_mod, only : sourceMass, amrgridsize, maxdepthamr, amrgridcentrez
    integer :: i
    type(octal), pointer :: thisOctal, child
    integer :: subcell
    real(double) :: vkep,  r
    type(vector) :: rvec

   do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call imposeKeplerianVelocity(child)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          r = abs(rVec%x)*1.d10
  !        dx = amrgridsize/(2.**maxdepthamr)
          !          upperbound = amrgridcentrez - 0.5d0*amrgridsize + dx
!  !        if(abs(rvec%x) > 2.d0*dx) then
!   !       if(abs(rvec%x) > 3.d0*dx) then
          vkep = sqrt(bigG*sourceMass(1)/(r))
          thisOctal%rhov(subcell) = thisOctal%rho(subcell) * vkep*r
!          endif
             !          endif
          !          thisOctal%rhov(subcell) = thisOctal%rho(subcell) * thisVel


!          n = 1.125
 !         rVec = subcellCentre(thisOctal,subcell)
  !        r = real(modulus(rVec*1.d10))
   !       x = rvec%x*1.d10
    !      vkep = sqrt(bigG*sourceMass(1)/(x**3))
     !     cs = soundSpeed(thisOctal, subcell)
      !!    eta = n *(cs**2/vkep**2)
        !  thisVel = vkep*(1.d0-eta*(abs(rVec%z*1.d10)/x**3)**2)
         ! thisOctal%rhov(subcell) = thisOctal%rho(subcell) * thisVel*x



       endif
    enddo
  end subroutine imposeKeplerianVelocity

  recursive subroutine updatedensitytree(thisoctal)
    use mpi
    type(octal), pointer   :: thisoctal, parentoctal, testoctal
    type(octal), pointer  :: child 
    integer :: i, n, m, j
    real(double) :: mass(8), totMass, v(8)



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

          mass = 0.d0
          do j = 1, n
             v(j) = cellVolume(testOctal, j)
             mass(j) =  testoctal%rho(j) * v(j) * 1.d30
          enddo
          totMass = sum(mass(1:n))

          if (cellVolume(parentOctal,m) > 0.d0) then
             parentOctal%rho(m) = totmass / (cellVolume(parentOctal,m)*1.d30)
          endif


          parentoctal%phi_gas(m) = sum(testoctal%phi_gas(1:n))/dble(n)

                
          parentoctal%edgecell(m) = any(testoctal%edgecell(1:n))
          testoctal => parentoctal
       enddo
    endif
  end subroutine updatedensitytree

  recursive subroutine SumPhiUptree(thisoctal, nDepth)
    use mpi
    type(octal), pointer   :: thisoctal, parentoctal, testoctal
    type(octal), pointer  :: child 
    integer :: i, n, m, nDepth
    real(double) :: sumphi



    if (thisoctal%nchildren > 0) then
       do i = 1, thisoctal%nchildren, 1
          child => thisoctal%child(i)
          call sumPhiUptree(child, nDepth)
       end do
    else

          
       testoctal => thisoctal
       sumPhi = 0.d0
       do while(testoctal%ndepth > nDepth)
             
          parentoctal => testoctal%parent
          m = testoctal%parentsubcell
          if (testoctal%twod) then
             n = 4
          else
             n = 8
          endif

          sumPhi = sumPhi +  sum(testoctal%phi_gas(1:n))/dble(n)
          testoctal => parentoctal
          m = testOctal%parentSubcell
       enddo
       if (.not.testOctal%ghostCell(m)) then
          testOctal%phi_gas(m) = sumPhi
       endif
    endif
  end subroutine SumPhiUptree

  recursive subroutine zeroFlux(thisOctal)
    type(octal), pointer :: thisOctal, child
    integer :: i

    if (.not.associated(thisOctal%flux_i)) allocate(thisOctal%flux_i(1:thisOctal%maxChildren))
    thisOctal%flux_i = thisOctal%phi_gas

    do i = 1, thisOctal%nChildren
       child => thisOctal%child(i)
       call zeroFlux(child)
    enddo
  end subroutine zeroFlux
  

  recursive subroutine restrictResiduals(thisoctal, nDepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: i, nDepth

    if (thisOctal%nDepth == nDepth) then
       if (.not.associated(thisOctal%parent%chiline)) allocate(thisOctal%parent%chiline(1:thisOctal%parent%maxChildren))


       thisOctal%parent%chiLine(thisOctal%parentSubcell) = &
            SUM(thisOctal%chiline(1:thisOctal%maxChildren))/dble(thisOctal%maxChildren)
    else

       if (thisoctal%nchildren > 0) then
          do i = 1, thisoctal%nchildren, 1
             child => thisoctal%child(i)
             call restrictResiduals(child, nDepth)
          end do
       endif
    endif
  end subroutine restrictResiduals

  recursive subroutine addPhiFromBelow(thisoctal, nDepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: i, nDepth

    if (thisOctal%nDepth == (nDepth+1)) then


!       write(*,*) "parent ",thisOctal%parent%phi_gas(thisOctal%parentSubcell),thisOctal%parent%ndepth
!       write(*,*) "this ",thisOctal%phi_gas(1:thisOctal%maxChildren)
       thisOctal%parent%phi_gas(thisOctal%parentSubcell) = &
       thisOctal%parent%phi_gas(thisOctal%parentSubcell) + &
            SUM(thisOctal%phi_gas(1:thisOctal%maxChildren))/dble(thisOctal%maxChildren)
    else

       if (thisoctal%nchildren > 0) then
          do i = 1, thisoctal%nchildren, 1
             child => thisoctal%child(i)
             call addPhiFromBelow(child, nDepth)
          end do
       endif
    endif
  end subroutine addPhiFromBelow



  recursive subroutine setSourceToFourPiRhoG(thisoctal, nDepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: i, nDepth

    if (thisOctal%nDepth == nDepth) then
       if (.not.associated(thisOctal%source)) allocate(thisOctal%source(1:thisOctal%maxChildren))
       thisOctal%source(1:thisOctal%maxChildren) = thisOctal%rho(1:thisOctal%maxChildren) * bigG * fourPi
    else
       if (thisoctal%nchildren > 0) then
          do i = 1, thisoctal%nchildren, 1
             child => thisoctal%child(i)
             call setSourceToFourPiRhoG(child, nDepth)
          end do
       endif
    endif


  end subroutine setSourceToFourPiRhoG

  recursive subroutine setSourceToResiduals(thisoctal, nDepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: i, nDepth,subcell


    if (thisOctal%nDepth == nDepth) then
       if (.not.associated(thisOctal%source)) allocate(thisOctal%source(1:thisOctal%maxChildren))
       if (.not.associated(thisOctal%chiLine)) then
          allocate(thisOctal%chiLine(1:thisOctal%maxChildren))
          thisOctal%chiline = 0.d0
       endif
       do subcell = 1, thisOctal%maxChildren
          if (thisOctal%edgecell(subcell)) then
             thisOctal%source(subcell) = 0.d0
          else
             thisOctal%source(subcell) = thisOctal%chiLine(subcell)
          endif
       enddo
    else
       if (thisoctal%nchildren > 0) then
          do i = 1, thisoctal%nchildren, 1
             child => thisoctal%child(i)
             call setSourceToResiduals(child, nDepth)
          end do
       endif
    endif
  end subroutine setSourceToResiduals



  recursive subroutine prolongatePhiGas(thisoctal, ndepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: i, ndepth


    if (thisoctal%ndepth == ndepth) then

       thisOctal%phi_gas = thisOctal%phi_gas - SUM(thisOctal%phi_gas)/dble(thisOctal%maxChildren)


       do i = 1, thisOctal%maxChildren
          if (.not.thisOctal%edgeCell(i)) then
             thisOctal%phi_gas(i) = thisOctal%phi_gas(i) + thisOctal%parent%phi_gas(thisOctal%parentSubcell)
!             thisOctal%phi_gas(i) = thisOctal%parent%phi_gas(thisOctal%parentSubcell)
          endif
       enddo
    else

       do i = 1, thisOctal%nChildren
          child => thisOctal%child(i)
          call prolongatePhiGas(child, nDepth)
       enddo

    endif

  end subroutine prolongatePhiGas

  recursive subroutine addCorrections(thisoctal, ndepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, ndepth


    if (thisoctal%ndepth == ndepth) then

       do i = 1, thisOctal%maxChildren
          if (.not.thisOctal%edgeCell(i)) then
             thisOctal%phi_gas(i) = &
                  thisOctal%phi_gas(i) + thisOctal%parent%phi_gas(thisOctal%parentSubcell)
          endif
       enddo
    else

       do subcell = 1, thisoctal%maxchildren
          if (thisoctal%haschild(subcell)) then
             ! find the child
             do i = 1, thisoctal%nchildren, 1
                if (thisoctal%indexchild(i) == subcell) then
                   child => thisoctal%child(i)
                   call addCorrections(child, ndepth)
                   exit
                end if
             end do
          endif
       enddo
    endif

  end subroutine addCorrections
  


  recursive subroutine updatephitree(thisoctal, ndepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, n, ndepth


    if (thisoctal%ndepth == ndepth) then

       if (thisoctal%twod) then
          n = 4
       else
          n = 8
       endif

       thisoctal%phi_gas = thisoctal%parent%phi_gas(thisoctal%parentsubcell)
!       if (any(thisOctal%phi_gas > 0.d0)) write(*,*) myrankGlobal," phi ",thisoctal%phi_gas(1:4),thisOctal%ghostCell(1:4)
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

  recursive subroutine updatephitreeNew(thisoctal, ndepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, n, ndepth, j
    real(double) :: mass(8)


    if (thisoctal%ndepth == (ndepth+1)) then

       if (thisoctal%twod) then
          n = 4
       else
          n = 8
       endif
       do j = 1, n 
          mass(j) = thisOctal%rho(j) * cellVolume(thisOctal, j) * 1.d30
       enddo

!       thisOctal%phi_gas = thisOctal%phi_gas - (SUM(thisOctal%phi_gas(1:thisOctal%maxChildren)))/dble(thisOctal%maxChildren)

       thisOctal%phi_gas = thisOctal%phi_gas - &
            (SUM(thisOctal%phi_gas(1:thisOctal%maxChildren)*mass(1:thisOctal%maxChildren))/sum(mass(1:thisOctal%maxChildren)))

       thisoctal%phi_gas = thisoctal%parent%phi_gas(thisoctal%parentsubcell) + thisOctal%phi_gas
!       if (any(thisOctal%phi_gas > 0.d0)) write(*,*) myrankGlobal," phi ",thisoctal%phi_gas(1:4),thisOctal%ghostCell(1:4)
    else

       do subcell = 1, thisoctal%maxchildren
          if (thisoctal%haschild(subcell)) then
             ! find the child
             do i = 1, thisoctal%nchildren, 1
                if (thisoctal%indexchild(i) == subcell) then
                   child => thisoctal%child(i)
                   call updatephitreeNew(child, ndepth)
                   exit
                end if
             end do
          endif
       enddo
    endif

  end subroutine updatephitreeNew

  recursive subroutine zeroPhiGasLevel(thisoctal, ndepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, ndepth


    if (thisoctal%ndepth == ndepth) then

       if (.not.associated(thisOctal%adot)) then
          allocate(thisOctal%adot(1:thisOctal%maxChildren))
          thisOctal%adot = 0.d0
       endif

       do subcell = 1, thisOctal%maxChildren
          thisOctal%adot(subcell) = thisOctal%phi_gas(subcell)
          if (.not.thisOctal%edgeCell(subcell)) then
             thisoctal%phi_gas(subcell) = 0.d0
          endif
       enddo
    else

       do i = 1, thisoctal%nchildren, 1
          child => thisoctal%child(i)
          call zeroPhiGasLevel(child, ndepth)
       end do

    endif

  end subroutine zeroPhiGasLevel

  recursive subroutine zeroPhiGasLevelEverywhere(thisoctal, ndepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, ndepth


    if (thisoctal%ndepth == ndepth) then

       if (.not.associated(thisOctal%adot)) then
          allocate(thisOctal%adot(1:thisOctal%maxChildren))
          thisOctal%adot = 0.d0
       endif

       do subcell = 1, thisOctal%maxChildren
          thisOctal%adot(subcell) = thisOctal%phi_gas(subcell)
          thisoctal%phi_gas(subcell) = 0.d0
       enddo
    else

       do i = 1, thisoctal%nchildren, 1
          child => thisoctal%child(i)
          call zeroPhiGasLevelEverywhere(child, ndepth)
       end do

    endif

  end subroutine zeroPhiGasLevelEverywhere

  recursive subroutine restorePhiGas(thisoctal, ndepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: i, ndepth


    if (thisoctal%ndepth == ndepth) then

       thisOctal%phi_gas = thisOctal%adot
    else

       do i = 1, thisoctal%nchildren, 1
          child => thisoctal%child(i)
          call restorePhiGas(child, ndepth)
       end do
    endif

  end subroutine restorePhiGas

  recursive subroutine storePhiGas(thisoctal, ndepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: i, ndepth


    if (thisoctal%ndepth == ndepth) then

       thisOctal%adot  = thisOctal%phi_gas  
    else

       do i = 1, thisoctal%nchildren, 1
          child => thisoctal%child(i)
          call storePhiGas(child, ndepth)
       end do
    endif

  end subroutine storePhiGas


  subroutine setupAlphaViscosity(grid, alpha, HoverR)
!    use inputs_mod, only : amrGridSize, minDepthAMR, rhoThreshold
    use mpi
    type(GRIDTYPE) :: grid
    real(double) :: alpha, HOverR, r
    integer :: j
    integer, parameter :: tag = 54
    integer :: ierr
    real(double), allocatable :: rAxis(:), mr(:),thisRAxis(:)
    integer :: nr, i, nr1, thisNr
    integer :: status(MPI_STATUS_SIZE)

    allocate(rAxis(1:100000),mr(1:100000))
    allocate(thisrAxis(1:10000))

    thisnr = 0
    nr = 0
    call getRadialGrid(grid%octreeRoot, thisnr, thisrAxis)
    nr = thisNR
    rAxis(1:nr) = thisrAxis(1:nr)
    do i = 1, nHydroThreadsGlobal
       if (i == myrankGlobal) then
          do j = 1, nHydroThreadsGlobal
             if (j /= myrankGlobal) then
                call mpi_recv(nr1, 1, MPI_INTEGER, j, tag, localWorldCommunicator, status, ierr)
                call mpi_recv(rAxis(nr+1:nr+nr1), nr1, MPI_DOUBLE_PRECISION, j, tag, localWorldCommunicator, status, ierr)
                nr = nr + nr1
             endif
          enddo
       else
          call MPI_SEND(thisnr, 1, MPI_INTEGER, i, tag, localWorldCommunicator, ierr)
          call MPI_SEND(thisrAxis, thisnr, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, ierr)
       endif
    enddo

    rAxis(2:nr+1) = rAxis(1:nr)
    rAxis(1) = 0.d0
    nr = nr + 1
    call sort(nr, rAxis)

    if (myrankGlobal /= 1) then
       call findTotalMassWithinRServer(grid, 1)
    else
       do j = 1, nr
          call findMassOverAllThreadsWithinR(grid, mr(j), raxis(j))
       enddo
       do j = 2, nHydroThreadsGlobal
          r = 1.d30
          call MPI_SEND(r, 1, MPI_DOUBLE_PRECISION, j, tag, localWorldCommunicator, ierr)
       enddo
    endif
    call MPI_BCAST(mr, nr, MPI_DOUBLE_PRECISION, 0, amrCommunicator, ierr)
!    if (myrankGlobal == 1) then
!       do i = 1, nr
!          write(*,*) i,rAxis(i),mr(i)
!       enddo
!    endif
    call setupAlphaViscosityPrivate(grid, grid%octreeRoot, alpha, HoverR, rAxis, mr, nr)
    deallocate(mr, rAxis, thisRaxis)

  end subroutine setupAlphaViscosity

  recursive subroutine getRadialGrid(thisoctal, nr, rAxis)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, nr
    real(double) :: rAxis(:)
    type(VECTOR) :: rVec
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call getRadialGrid(child, nr, rAxis)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.thisOctal%ghostCell(subcell)) then
             rVec = subcellCentre(thisOctal, subcell)
             if ( ((rVec%z) > 0.d0).and. &
                  ((rVec%z-thisOctal%subcellSize/2.d0-0.1d0*smallestCellSize) < 0.d0) ) then
                nr = nr + 1
                rAxis(nr) = rVec%x
             endif
          endif
       endif
    enddo
  end subroutine getRadialGrid


  recursive subroutine setupAlphaViscosityPrivate(grid, thisoctal, alpha, HoverR, rAxis, mr, nr)
    use mpi
    type(octal), pointer   :: thisoctal, neighbourOctal
    type(octal), pointer  :: child 
    integer :: subcell, i, nr, neighbourSubcell
    type(GRIDTYPE) :: grid
    real(double) :: alpha, HoverR, r, mass, omegaK, rAxis(:), mr(:)
    real(double), parameter :: Qcrit = 5.d0
    real(double) :: toomreQ, cs, sigma, fac, thisRho
    type(VECTOR) :: rVec, rVec2
    

  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupAlphaViscosityPrivate(grid, child, alpha, HoverR, rAxis, mr, nr)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.associated(thisOctal%etaline)) then
             allocate(thisOctal%etaline(1:thisOctal%maxChildren))
          endif
          thisOCtal%etaline(subcell) = 0.d0
          if (.not.associated(thisOctal%etaCont)) then
             allocate(thisOctal%etaCont(1:thisOctal%maxChildren))
          endif
          thisOCtal%etaCont(subcell) = 0.d0
          if (.not.thisoctal%edgeCell(subcell)) then
             rVec = subcellCentre(thisOctal, subcell)
             r = sqrt(rVec%x**2 + rVec%y**2)
             if (r == 0.d0) then
                write(*,*) thisOctal%nDepth, rVec,thisOctal%ghostCell(subcell)
             endif
             call locate(rAxis, nr, r, i)
             mass = mr(i) + (mr(i+1)-mr(i))*(r - rAxis(i))/(rAxis(i+1)-rAxis(i))
             mass = mass + SUM(globalSourceArray(1:GlobalnSource)%mass)
             omegaK = sqrt(bigG * mass / (r*gridDistanceScale)**3)

             sigma = thisOctal%rho(subcell) * thisOctal%subcellSize * gridDistanceScale
             cs =soundSpeed(thisOctal, subcell)
             toomreQ = (cs * omegaK) / (pi * bigG * sigma)
             thisOctal%etaCont(subcell) = toomreQ
             fac = 1.d0
             if (.not.thisOctal%ghostCell(subcell)) then
                if (toomreQ < Qcrit) then
                   fac = min((Qcrit**2 / toomreQ**2),100.d0)
!                   write(*,*) "q ",toomreQ,fac
                endif
             endif

             thisRho = thisOctal%rho(subcell)
             rVec2 = rVec + VECTOR(thisOctal%subcellSize/2.d0+0.1d0*smallestCellsize, 0.d0, 0.d0)
             if ((modulus(rVec - globalSourceArray(1)%position)*1.d10 < globalSourceArray(1)%accretionRadius).and. &
                 (modulus(rVec2 - globalSourceArray(1)%position)*1.d10 > globalSourceArray(1)%accretionRadius)) then
                if ( ((rVec%z + thisOctal%subcellSize/2.d0)> 0.d0).and.((rVec%z - thisOctal%subcellSize/2.d0) < 0.d0)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(rVec2,neighbourOctal, neighbourSubcell)
                   if (neighbourOctal%rho(neighbourSubcell) > rhoThreshold) then
                      thisRho = neighbourOctal%rho(neighbourSubcell)
!                      write(*,*) "setting rho from ",thisOctal%rho(subcell), " to " ,thisRho
                   endif
                endif
             endif


             thisOctal%etaLine(subcell) = fac * alpha * omegaK * (r*gridDistanceScale)**2 * hOverR**2 * thisRho


             if (alpha < 0.d0) then
                thisOctal%etaLine(subcell) = -alpha
             endif

          endif
       endif
    enddo
  end subroutine setupAlphaViscosityPrivate


!Set up the i-1/2 flux for each cell on the grid
  recursive subroutine constructflux(thisoctal, dt)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx

  
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

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             if (thisoctal%u_interface(subcell).ge.0.d0) then
                thisoctal%flux_i(subcell) = thisoctal%u_interface(subcell) * thisoctal%q_i_minus_1(subcell)
             else
                thisoctal%flux_i(subcell) = thisoctal%u_interface(subcell) * thisoctal%q_i(subcell)
             endif
        
             if (thisoctal%x_i(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) "problem with the x_i values"
                stop
 !            else
!                print *, "x_i values were good"
             endif


             dx = thisOctal%subcellSize*griddistancescale !old
             dx = thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell)

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

!Set up the i-1/2 flux for each cell on the grid
  recursive subroutine constructfluxCylindrical(grid, thisoctal, dt, direction, writeDebug)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    type(VECTOR) :: direction, locator
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: dt, dx, speed
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, xnext, px, py, pz
    real(double) :: qViscosity(3,3), rm1, um1, pm1
    integer :: nd, nc
    logical :: debug, writeDebug
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call constructfluxCYlindrical(grid, child, dt, direction, writeDebug)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (thisOctal%edgeCell(subcell)) cycle
          locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*smallestCellsize)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
               rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

          ! need to construct flux correctly at this point
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


                

             dx = thisOctal%subcellSize*griddistancescale !old
             dx = thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell)
!
             thisoctal%flux_i(subcell) = thisoctal%flux_i(subcell) + &
                  0.5d0 * abs(thisoctal%u_interface(subcell)) * &
                  (1.d0 - abs(thisoctal%u_interface(subcell) * dt / &
                   dx)) * &
                  thisoctal%philimit(subcell) * (thisoctal%q_i(subcell) - thisoctal%q_i_minus_1(subcell))

          debug = .false.
          speed = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)/1.d5
          if ((speed > 2000.d0).and.(myHydroSetGlobal == 0)) then
             debug = .true.
          endif
          if (debug.and.writeDebug) then
             write(*,*) "thisOctal%phiLimit ", thisOctal%phiLimit(subcell)
             write(*,*) "q_i ",thisOctal%q_i(subcell)
             write(*,*) "q_i_minus_1 ",thisOctal%q_i_minus_1(subcell)
          endif

          endif
!          if (thisoctal%flux_i(subcell) > 1.d-10) write(*,*) "flux ",thisoctal%flux_i(subcell),thisoctal%u_interface(subcell), &
!               thisoctal%rho(subcell), thisoctal%rhou(subcell)
       endif
    enddo
  end subroutine constructfluxCylindrical

!Set up the i-1/2 flux for each cell on the grid
  recursive subroutine constructfluxCylindrical2(grid, thisoctal, dt, direction, writeDebug)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    type(VECTOR) :: direction, locator
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell, j
    real(double) :: dt, dx
    real(double) :: q(2), rho(2), rhoe(2), rhou(2), rhov(2), rhow(2), x, qnext, pressure(2), flux(2), phi, phigas, xnext, px, py, pz
    real(double) :: qViscosity(3,3), rm1, um1, pm1
    integer :: nd
    logical :: writeDebug
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call constructfluxCYlindrical2(grid, child, dt, direction, writeDebug)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (thisOctal%edgeCell(subcell)) cycle
          locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*smallestCellsize)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          call getneighbourvalues2(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
               rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

          ! need to construct flux correctly at this point
          if (.not.thisoctal%edgecell(subcell)) then


             do j = 1, 2
                if (thisoctal%u_amr_interface(subcell,j).ge.0.d0) then
                   thisoctal%flux_amr_i(subcell,j) = thisoctal%u_amr_interface(subcell,j) * thisoctal%q_amr_i_minus_1(subcell,j)
                else
                   thisoctal%flux_amr_i(subcell,j) = thisoctal%u_amr_interface(subcell,j) * thisoctal%q_i(subcell)
                endif
                
                if (thisoctal%x_i(subcell) == thisoctal%x_i_minus_1(subcell)) then
                   write(*,*) "problem with the x_i values"
                   stop
                endif


                

                dx = thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell)

                thisoctal%flux_amr_i(subcell,j) = thisoctal%flux_amr_i(subcell, j) + &
                     0.5d0 * abs(thisoctal%u_amr_interface(subcell,j)) * &
                     (1.d0 - abs(thisoctal%u_amr_interface(subcell,j) * dt / &
                     dx)) * &
                     thisoctal%philimit_amr(subcell ,j) * (thisoctal%q_i(subcell) - thisoctal%q_amr_i_minus_1(subcell,j))

             enddo
          endif
       endif
    enddo
  end subroutine constructfluxCylindrical2

!Set up the i-1/2 flux for each cell on the grid
  recursive subroutine constructflux_amr(grid, thisoctal, dt, direction, writeDebug)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    type(VECTOR) :: direction, locator
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell, j
    real(double) :: dt, dx
    real(double) :: q(4), rho(4), rhoe(4), rhou(4), rhov(4), rhow(4), x, qnext, pressure(4), flux(4), phi, phigas, xnext, px, py, pz
    real(double) :: qViscosity(3,3), rm1, um1, pm1
    integer :: nd
    logical :: writeDebug
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call constructflux_amr(grid, child, dt, direction, writeDebug)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (thisOctal%edgeCell(subcell)) cycle
          locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*smallestCellsize)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          call getneighbourvalues4(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
               rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

          ! need to construct flux correctly at this point
          if (.not.thisoctal%edgecell(subcell)) then


             do j = 1, 4
                if (thisoctal%u_amr_interface(subcell,j).ge.0.d0) then
                   thisoctal%flux_amr_i(subcell,j) = thisoctal%u_amr_interface(subcell,j) * thisoctal%q_amr_i_minus_1(subcell,j)
                else
                   thisoctal%flux_amr_i(subcell,j) = thisoctal%u_amr_interface(subcell,j) * thisoctal%q_i(subcell)
                endif
                
                if (thisoctal%x_i(subcell) == thisoctal%x_i_minus_1(subcell)) then
                   write(*,*) "problem with the x_i values"
                   stop
                endif


                

                dx = thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell)

                thisoctal%flux_amr_i(subcell,j) = thisoctal%flux_amr_i(subcell, j) + &
                     0.5d0 * abs(thisoctal%u_amr_interface(subcell,j)) * &
                     (1.d0 - abs(thisoctal%u_amr_interface(subcell,j) * dt / &
                     dx)) * &
                     thisoctal%philimit_amr(subcell ,j) * (thisoctal%q_i(subcell) - thisoctal%q_amr_i_minus_1(subcell,j))

             enddo
          endif
       endif
    enddo
  end subroutine constructflux_amr


!Set up the i-1/2 flux for each cell on the grid
  recursive subroutine constructfluxSpherical(grid, thisoctal, dt, direction, writeDebug)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    type(VECTOR) :: direction, locator
    type(octal), pointer  :: child 
    integer :: subcell, i, neighbourSubcell
    real(double) :: dt, dx, speed
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, xnext, px, py, pz
    real(double) :: qViscosity(3,3), rm1, um1, pm1
    integer :: nd, nc
    logical :: debug, writeDebug
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call constructfluxSpherical(grid, child, dt, direction, writeDebug)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (thisOctal%edgeCell(subcell)) cycle
          locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*smallestCellsize)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
               rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

          ! need to construct flux correctly at this point
          if (.not.thisoctal%edgecell(subcell)) then
             if (thisoctal%u_interface(subcell).ge.0.d0) then
                thisoctal%flux_i(subcell) = thisoctal%u_interface(subcell) * thisoctal%q_i_minus_1(subcell)
             else
                thisoctal%flux_i(subcell) = thisoctal%u_interface(subcell) * thisoctal%q_i(subcell)
             endif
        
             if (thisoctal%x_i(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) "problem with the x_i values ",myrankglobal,thisOctal%nDepth
                stop
             endif

             dx = thisOctal%subcellSize*griddistancescale !old
             dx = thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell)
!
             thisoctal%flux_i(subcell) = thisoctal%flux_i(subcell) + &
                  0.5d0 * abs(thisoctal%u_interface(subcell)) * &
                  (1.d0 - abs(thisoctal%u_interface(subcell) * dt / &
                   dx)) * &
                  thisoctal%philimit(subcell) * (thisoctal%q_i(subcell) - thisoctal%q_i_minus_1(subcell))

          debug = .false.
          speed = sqrt(thisOctal%rhou(subcell)**2)/thisOctal%rho(subcell)/1.d5
          if ((speed > 2000.d0).and.(myHydroSetGlobal == 0)) then
             debug = .true.
          endif
          if (debug.and.writeDebug) then
             write(*,*) "thisOctal%phiLimit ", thisOctal%phiLimit(subcell)
             write(*,*) "q_i ",thisOctal%q_i(subcell)
             write(*,*) "q_i_minus_1 ",thisOctal%q_i_minus_1(subcell)
          endif

          endif
!          if (thisoctal%flux_i(subcell) > 1.d-10) write(*,*) "flux ",thisoctal%flux_i(subcell),thisoctal%u_interface(subcell), &
!               thisoctal%rho(subcell), thisoctal%rhou(subcell)
       endif
    enddo
  end subroutine constructfluxSpherical

!Perform the actual advection
  recursive subroutine updatecellq(thisoctal, dt)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx, df

  
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

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if ((.not.thisoctal%ghostcell(subcell)).and.(.not.thisOctal%boundaryCell(subcell))) then
          
!             dx = 0.5*(thisoctal%x_i_plus_1(subcell) - thisoctal%x_i_minus_1(subcell))
             dx = thisOctal%subcellSize*griddistancescale !old
!             dx = (thisoctal%x_i_plus_1(subcell) - thisoctal%x_i(subcell)) !new

             df = (thisoctal%flux_i_plus_1(subcell) - thisoctal%flux_i(subcell))
             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * (df/dx)


!             dv = cellVolume(thisOctal, subcell)*1.d30
!             area_i_plus_1 = dx**2
!             area_i = dx**2
!             df = (thisoctal%flux_i_plus_1(subcell) * area_i_plus_1 - thisoctal%flux_i(subcell) * area_i)
!             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * (df/dv)


          endif
       endif
    enddo
  end subroutine updatecellq

!Perform the actual advection
  recursive subroutine updatecellqCylindrical(thisoctal, dt, direction)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    type(VECTOR) :: direction, rVec
    integer :: subcell, i
    logical :: radial
    real(double) :: dt, dv, df, area_i, area_i_plus_1

  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call updatecellqCylindrical(child, dt, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if ((.not.thisoctal%ghostcell(subcell))) then !.and.(.not.thisOctal%boundaryCell(subcell))) then

             radial = .true. 
             if (direction%x < 0.1d0) then
                radial = .false.
             endif
             rVec = subcellCentre(thisOctal,subcell)
             if (radial) then
                area_i = twoPi * ((thisOctal%x_i(subcell)-gridDistanceScale*thisOctal%subcellSize/2.d0) &
                     * thisOctal%subcellSize * gridDistanceScale)
                area_i_plus_1 = twoPi * ((thisOctal%x_i(subcell)+gridDistanceScale*thisOctal%subcellSize/2.d0) &
                     * thisOctal%subcellSize * gridDistanceScale)
             else
                area_i = pi * ((rVec%x*gridDistanceScale + thisOctal%subcellSize*gridDistanceScale/2.d0)**2 - &
                     (rVec%x*gridDistanceScale-thisOctal%subcellSize*gridDistanceScale/2.d0)**2)
                area_i_plus_1 = area_i
             endif

             dv = cellVolume(thisOctal, subcell)*1.d30

             df = (thisoctal%flux_i_plus_1(subcell) * area_i_plus_1 - thisoctal%flux_i(subcell) * area_i)

             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * (df/dv)

          endif
       endif
    enddo
  end subroutine updatecellqCylindrical

!Perform the actual advection
  recursive subroutine updatecellqCylindrical2(thisoctal, dt, direction)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    type(VECTOR) :: direction, rVec
    integer :: subcell, i, j
    logical :: radial, debug
    real(double) :: dt, dv, df, area_i(2), area_i_plus_1(2)

  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call updatecellqCylindrical2(child, dt, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if ((.not.thisoctal%ghostcell(subcell))) then


             debug = .false.
             if (inSubcell(thisOctal, subcell, VECTOR(20d3, 0.d0, 1.14e6))) debug=.true.
             if (inSubcell(thisOctal, subcell, VECTOR(10d3, 0.d0, 1.11e6))) debug=.true.
             radial = .true. 
             if (abs(direction%x) < 0.1d0) then
                radial = .false.
             endif
             rVec = subcellCentre(thisOctal,subcell)

             if (radial) then
                area_i(1:2) = twoPi * ((thisOctal%x_i(subcell)-gridDistanceScale*thisOctal%subcellSize/2.d0) &
                     * thisOctal%subcellSize * gridDistanceScale)/2.d0
!                if (thisOctal%x_i(subcell)-gridDistanceScale*thisOctal%subcellSize/2.d0 < 0.1d0*thisOctal%subcellSize*gridDistanceScale) then
!                   write(*,*) "area on axis"
!                   area_i(1:2) = twoPi * ((thisOctal%x_i(subcell)+gridDistanceScale*thisOctal%subcellSize/2.d0) &
!                        * thisOctal%subcellSize * gridDistanceScale)/2.d0
!                endif
                   
                area_i_plus_1(1:2) = twoPi * ((thisOctal%x_i(subcell)+gridDistanceScale*thisOctal%subcellSize/2.d0) &
                     * thisOctal%subcellSize * gridDistanceScale)/2.d0
             else
                area_i(1) = pi * ((rVec%x*gridDistanceScale)**2 - &
                     (rVec%x*gridDistanceScale-thisOctal%subcellSize*gridDistanceScale/2.d0)**2)
                area_i(2) = pi * ((rVec%x*gridDistanceScale+thisOctal%subcellSize*gridDistanceScale/2.d0)**2 - &
                     (rVec%x*gridDistanceScale)**2)

                area_i_plus_1(1:2) = area_i(1:2)
             endif

             dv = cellVolume(thisOctal, subcell)*1.d30
             do j = 1, 2
                df = (thisoctal%flux_amr_i_plus_1(subcell,j) * area_i_plus_1(j) - &
                     thisoctal%flux_amr_i(subcell,j) * area_i(j))
                thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * (df/dv)
!                if (debug.and..not.radial) then
!                   write(*,*) "j ", j, thisOctal%subcellSize
!                   write(*,*) "flux ", thisOctal%flux_amr_i_plus_1(subcell,j), thisOctal%flux_amr_i(subcell,j), area_i_plus_1(j), area_i(j)
!                   write(*,*) "u_i ",thisOctal%u_amr_interface(subcell,j)
!                Endif
             enddo

          endif
       endif
    enddo
  end subroutine updatecellqCylindrical2

!Perform the actual advection
  recursive subroutine updatecellq_amr(thisoctal, dt, direction)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    type(VECTOR) :: direction, rVec
    integer :: subcell, i, j
    real(double) :: dt, dv, df, area_i(4), area_i_plus_1(4)

  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call updatecellq_amr(child, dt, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if ((.not.thisoctal%ghostcell(subcell))) then


             debug = .false.
             rVec = subcellCentre(thisOctal,subcell)

!             if (inSubcell(thisOctal, subcell, VECTOR(-20d6, 1.d3, 7e6))) debug=.true.
!             if (inSubcell(thisOctal, subcell, VECTOR(-20d6, 1.d3, 8e6))) debug=.true.

             if (debug) then
                write(*,*) "vec ",rvec
                write(*,*) "flux_i        ", thisOctal%flux_amr_i(subcell,:)
                write(*,*) "flux_i_plus_1 ", thisOctal%flux_amr_i_plus_1(subcell,:)
                write(*,*) "u_amr_interface ", thisOctal%u_amr_interface(subcell,:) 
                write(*,*) "q_i   ", thisOctal%q_i(subcell)
                write(*,*) "q_i+1 ", thisOctal%q_amr_i_plus_1(subcell,:)
                write(*,*) "  "
             endif
             area_i(1:4) = ((thisOctal%subcellSize * gridDistanceScale)**2) /4.d0
             area_i_plus_1(1:4) = ((thisOctal%subcellSize * gridDistanceScale)**2)/4.d0
             dv = cellVolume(thisOctal, subcell)*1.d30
             do j = 1, 4
                df = (thisoctal%flux_amr_i_plus_1(subcell,j) * area_i_plus_1(j) - &
                     thisoctal%flux_amr_i(subcell,j) * area_i(j))
                thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * (df/dv)
             enddo

          endif
       endif
    enddo
  end subroutine updatecellq_amr


!Perform the actual advection
  recursive subroutine updatecellqSpherical(thisoctal, dt, direction)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    type(VECTOR) :: direction, rVec
    integer :: subcell, i
    real(double) :: dt, dv, df, area_i, area_i_plus_1

  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call updatecellqspherical(child, dt, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if ((.not.thisoctal%ghostcell(subcell))) then !.and.(.not.thisOctal%boundaryCell(subcell))) then



             rVec = subcellCentre(thisOctal,subcell)

             area_i = fourPi * (thisOctal%x_i(subcell)-gridDistanceScale*thisOctal%subcellSize/2.d0)**2

             area_i_plus_1 = fourPi * (thisOctal%x_i(subcell)+gridDistanceScale*thisOctal%subcellSize/2.d0)**2
             


             dv = cellVolume(thisOctal, subcell)*1.d30

             df = (thisoctal%flux_i_plus_1(subcell) * area_i_plus_1 - thisoctal%flux_i(subcell) * area_i)

             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * (df/dv)

          endif
       endif
    enddo
  end subroutine updatecellqSpherical
  
  recursive subroutine synchronizefluxes(thisoctal, dt, idepth)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, idepth
    real(double) :: dt !, dx

  
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
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (thisoctal%ndepth /= idepth) cycle

          if (.not.thisoctal%ghostcell(subcell)) then

!             thisoctal%q_i(subcell) = thisoctal%q_i(subcell) - dt * &
!                  (thisoctal%flux_i_plus_1(subcell) - thisoctal%flux_i(subcell)) / dx
          endif
       endif
    enddo
  end subroutine synchronizefluxes


!boundary partners should be directly opposite their ghosts, this routine checks that this is the case
  recursive subroutine checkBoundaryPartners(thisOctal, grid)
    type(gridtype) :: grid
    type(octal), pointer :: thisOctal, bOctal
    type(octal), pointer :: child
    integer :: subcell, i, bSubcell
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
             bOctal => thisOctal
             bSubcell = subcell
             call findsubcelllocal(bvec, boctal, bsubcell)
!             if (bOctal%ghostCell(bSubcell)) then
!                write(*,*) "Partner for ghost cell is in turn a ghostcell"
!                write(*,*) "cell is at ",rVec
!                write(*,*) "partner is at ",bvec
!             endif

             direction = bVec - rVec

             !Make it a unit vector                                                                                                                                                                     
             if(direction%x > 0.d0) direction%x = direction%x / direction%x
             if(direction%x < 0.d0) direction%x = -direction%x / direction%x
             if(direction%y > 0.d0) direction%y = direction%y / direction%y
             if(direction%y < 0.d0) direction%y = -direction%y / direction%y
             if(direction%z > 0.d0) direction%z = direction%z / direction%z
             if(direction%z < 0.d0) direction%z = -direction%z / direction%z
 
!Check that boundary partners are directly opposite one another
             if((abs(direction%x) + abs(direction%y) + abs(direction%z)) /= 1.d0) then
                print *, "boundary partner is at a diagonal!"
                print *, "direction ", direction
                print *, "rVec ", rVec,thisOctal%nDepth
                print *, "bVec ", bVec,bOctal%nDepth
!                stop
             end if

          end if
          
       end if
    end do

  end subroutine checkBoundaryPartners


!Set up cell centre values for given propagation direction
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
!          write(*,*) "x_i ",thisOctal%x_i(subcell),myrankglobal
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

!set up neighbour positions and values of the advecting quantity
  recursive subroutine setupqx(thisoctal, grid, direction)
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi, phigas
    real(double) :: xnext, px, py, pz, qViscosity(3,3), rm1, um1, pm1
    integer :: subcell, i, neighboursubcell
    integer :: nd, nc
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


          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          thisoctal%x_i_minus_1(subcell) = 0.d0
          thisoctal%x_i_plus_1(subcell) = 0.d0

          if (.not.thisoctal%edgecell(subcell)) then

             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

             thisoctal%x_i_plus_1(subcell) = x
             thisoctal%q_i_plus_1(subcell) = q

             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%x_i_minus_1(subcell) = x
!             thisoctal%x_i_minus_2(subcell) = xnext
             thisoctal%q_i_minus_1(subcell) = q
             thisoctal%q_i_minus_2(subcell) = qnext

             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrankGlobal," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell)
             endif

          endif
       endif
    enddo
  end subroutine setupqx

!set up neighbour positions and values of the advecting quantity
  recursive subroutine setupqx2(thisoctal, grid, direction)
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe(2), rhou(2), rhov(2), rhow(2), rho(2), q(2), x, qnext, pressure(2), flux(2), phi, phigas
    real(double) :: xnext, px, py, pz, qViscosity(3,3), rm1, um1, pm1
    integer :: subcell, i, neighboursubcell
    integer :: nd
    type(vector) :: direction, locator, reversedirection


    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupqx2(child, grid, direction)
                exit
             end if
          end do
       else


          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          thisoctal%x_i_minus_1(subcell) = 0.d0
          thisoctal%x_i_plus_1(subcell) = 0.d0

          if (.not.thisoctal%edgecell(subcell)) then

             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues2(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

             thisoctal%x_i_plus_1(subcell) = x
             thisoctal%q_amr_i_plus_1(subcell,1:2) = q(1:2)

             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues2(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%x_i_minus_1(subcell) = x
!             thisoctal%x_i_minus_2(subcell) = xnext
             thisoctal%q_amr_i_minus_1(subcell,1:2) = q(1:2)
             thisoctal%q_i_minus_2(subcell) = qnext

             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrankGlobal," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell)
             endif

          endif
       endif
    enddo
  end subroutine setupqx2


!set up neighbour positions and values of the advecting quantity
  recursive subroutine setupqx_amr(thisoctal, grid, direction)
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe(4), rhou(4), rhov(4), rhow(4), rho(4), q(4), x, qnext, pressure(4), flux(4), phi, phigas
    real(double) :: xnext, px, py, pz, qViscosity(3,3), rm1, um1, pm1
    integer :: subcell, i, neighboursubcell
    integer :: nd
    type(vector) :: direction, locator, reversedirection


    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupqx_amr(child, grid, direction)
                exit
             end if
          end do
       else


          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          thisoctal%x_i_minus_1(subcell) = 0.d0
          thisoctal%x_i_plus_1(subcell) = 0.d0

          if (.not.thisoctal%edgecell(subcell)) then

             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues4(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

             thisoctal%x_i_plus_1(subcell) = x
             thisoctal%q_amr_i_plus_1(subcell,1:4) = q(1:4)

             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues4(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%x_i_minus_1(subcell) = x
!             thisoctal%x_i_minus_2(subcell) = xnext
             thisoctal%q_amr_i_minus_1(subcell,1:4) = q(1:4)
             thisoctal%q_i_minus_2(subcell) = qnext

             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrankGlobal," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell)
             endif

          endif
       endif
    enddo
  end subroutine setupqx_amr


!set up neighbour densities and gravtiational potentials
  recursive subroutine setuprhophi(thisoctal, grid, direction)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi, phigas, qViscosity(3,3)
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator, reversedirection
    integer :: nd, nc
    real(double) :: xnext, px, py, pz, rm1, um1, pm1

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

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

             thisoctal%rho_i_plus_1(subcell) = rho
             thisoctal%phi_i_plus_1(subcell) = phi

             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%rho_i_minus_1(subcell) = rho
             thisoctal%phi_i_minus_1(subcell) = phi

          endif
       endif
    enddo
  end subroutine setuprhophi

!set up neighbour rhov 
  recursive subroutine setuprhorvplusminus1(thisoctal, grid)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rhoe, rhou, rhov, rhow, rho, q, x, qnext, pressure, flux, phi, phigas, qViscosity(3,3)
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator, reversedirection
    integer :: nd, nc
    real(double) :: xnext, px, py, pz, rm1, um1, pm1

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setuprhorvplusminus1(child, grid)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          direction = VECTOR(1.d0, 0.d0, 0.d0)

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

             thisoctal%rhorv_i_plus_1(subcell) = rhov
             reversedirection = (-1.d0) * direction
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%rhorv_i_minus_1(subcell) = rhov

          endif
       endif
    enddo
  end subroutine setuprhorvplusminus1

!set up neighbour densities and cell interface advecting velocities - x direction
  recursive subroutine setupui(thisoctal, grid, direction, dt)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas, qViscosity(3,3)
    integer :: subcell, i, neighboursubcell, nd, nc
    type(vector) :: direction, locator!, rotator
    real(double) :: rhou_i_minus_1, rho_i_minus_1, dt
    real(double) :: rhou_i_plus_1, x_interface_i_p_half
    real(double) :: xnext, weight, px, py, pz, x_interface
    real(double) :: pressure_i_plus_1, pressure_i, pressure_i_minus_1, pressure_i_minus_2
    real(double) :: rho_i, rho_i_minus_2, rho_i_plus_1
    real(double) :: rm1, um1, pm1
    real(double) :: x_i, x_i_plus_1, x_i_minus_1, x_i_minus_2
    real(double) :: dx
    real(double) :: u_i, u_i_minus_1
    real(double) :: dpdx_i, dpdx_i_minus_half, dpdx_i_minus_1, thisRhou

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupui(child, grid, direction, dt)
                exit
             end if
          end do
       else
          
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          !thaw - edgecell is correct.
          if (.not.thisoctal%edgecell(subcell)) then !xxx changed fromghostcell
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz,rm1,um1, pm1, qViscosity)

             rho_i_minus_1 = rho
             rhou_i_minus_1 = direction.dot.VECTOR(rhou,rhov,rhow)
             pressure_i_minus_1 = pressure             
             x_i_minus_1 = px

             rho_i = thisOctal%rho(subcell)
             pressure_i = thisOctal%pressure_i(subcell)
             x_i = thisOctal%x_i(subcell)

             rho_i_minus_2 = rm1
             pressure_i_minus_2 = pm1
             x_i_minus_2 = xnext

             
             x_interface = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistancescale/2.d0

             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz,rm1,um1, pm1, qViscosity)

             rho_i_plus_1 = rho
             rhou_i_plus_1 = direction.dot.VECTOR(rhou,rhov,rhow)
             pressure_i_plus_1 = pressure             
             x_i_plus_1 = px
             x_interface_i_p_half = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistancescale/2.d0



             if (thisOctal%x_i(subcell) == thisOctal%x_i_minus_1(subcell)) then
                write(*,*) "x_i bug ", thisOctal%x_i(subcell), thisOctal%x_i_minus_1(subcell)
                print *, "cen ", subcellCentre(thisOctal, subcell)
                if (.not.octalonthread(neighbouroctal, neighboursubcell, myrankGlobal)) then
                   print *, "AND IT WASNT ON THREAD"
                end if
!                print *, "EDGE", thisOctal%edgecell(subcell)
 !               print *, "GHOST", thisOctal%ghostcell(subcell)
             else
 !               print *, "good at u_i"
             endif
             weight = 1.d0 - (thisOctal%x_i(subcell) - x_interface) / (thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell))


             if ((rho_i_minus_1 == 0.d0).or.(thisOctal%rho(subcell) == 0.d0)) then
                write(*,*) "bug in setupui ",rho_i_minus_1, thisOctal%rho(subcell)
             endif
             !This is the velocity for the i'th cell at i-1/2
             thisRhou = direction.dot.VECTOR(thisOctal%rhou(subcell), thisOctal%rhov(subcell), thisOctal%rhow(subcell))
             thisoctal%u_interface(subcell) = &
                  weight*thisRhou/thisoctal%rho(subcell) + &
                  (1.d0-weight)*rhou_i_minus_1/rho_i_minus_1
             if (.not.associated(thisOctal%u_i)) allocate(thisOctal%u_i(1:thisOctal%maxChildren))
             thisOctal%u_i(subcell) = thisRhou / thisOctal%rho(subcell)


! now the cell at i+1/2

             weight = 1.d0 - (thisOctal%x_i_plus_1(subcell) - x_interface_i_p_half) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i(subcell))
             
             if (.not.associated(thisOctal%u_interface_i_plus_1)) allocate(thisOctal%u_interface_i_plus_1(1:thisOctal%maxChildren))
             thisoctal%u_interface_i_plus_1(subcell) = &
                  weight * rhou_i_plus_1/rho_i_plus_1+ &
                  (1.d0-weight)*thisRhou/thisoctal%rho(subcell)


             ! for info, see http://ses.library.usyd.edu.au/bitstream/2123/376/2/adt-NU20010730.12021505_chapter_4.pdf
             if(rhiechow .and. .not. thisOctal%ghostcell(subcell)) then

                locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*&
                     grid%halfsmallestsubcell)
                neighbouroctal => thisoctal
                call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
                call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, &
                     xnext, px, py, pz, rm1, um1, pm1, qViscosity)

                rho_i_plus_1 = rho
                pressure_i_plus_1 = pressure
                x_i_plus_1 = px

                !now we have i,i+1,i-1 and i-2 can apply rhie-chow
                !you basically subtract the local pressure gradient and add the average ambient one

                u_i = thisRhou / thisOctal%rho(subcell)
                u_i_minus_1 = rhou_i_minus_1 / rho_i_minus_1

                dx = thisOctal%subcellSize * gridDistanceScale

                dpdx_i_minus_half = (pressure_i - pressure_i_minus_1) / dx
                dpdx_i = (pressure_i_plus_1 - pressure_i_minus_1) / (2.d0*dx)
                dpdx_i_minus_1 = (pressure_i - pressure_i_minus_2) / (2.d0*dx)


                if(dpdx_i_minus_1 /= 0.d0 .and. dpdx_i_minus_half /= 0.d0 .and. &
                     dpdx_i /= 0.d0) then
!the approximation screws up if applied to cells that are symmetric on one side only
!i.e. in a model with zero gradient boundaries you get preferential flow in in - direction
                   thisOctal%u_interface(subcell) = 0.5d0*(u_i + u_i_minus_1) - &
                        (dt/(0.5d0*(rho_i_minus_1 + rho_i)))*dpdx_i_minus_half + &
                        0.5d0 * ((dt/rho_i)*dpdx_i + (dt/rho_i_minus_1)*dpdx_i_minus_1)
                end if

             endif

             
          end if
       endif
    enddo

  end subroutine setupui

!set up neighbour densities and cell interface advecting velocities - x direction
  recursive subroutine setupui_amr(thisoctal, grid, direction, dt)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho(2), rhou(2), rhov(2), rhow(2), q(2), qnext, x, rhoe(2), pressure(2), flux(2), phi, phigas, qViscosity(3,3)
    integer :: subcell, i, neighboursubcell, nd, j
    type(vector) :: direction, locator!, rotator
    real(double) :: rhou_i_minus_1(2), rho_i_minus_1(2), dt
    real(double) :: rhou_i_plus_1(2), x_interface_i_p_half
    real(double) :: xnext, weight, px, py, pz, x_interface
    real(double) :: pressure_i_plus_1(2), pressure_i(2), pressure_i_minus_1(2), pressure_i_minus_2(2)
    real(double) :: rho_i, rho_i_minus_2, rho_i_plus_1(2)
    real(double) :: rm1, um1, pm1
    real(double) :: x_i, x_i_plus_1, x_i_minus_1, x_i_minus_2
    real(double) :: thisRhou
    real(double) :: dx
    real(double) :: u_i(2), u_i_minus_1(2)
    real(double) :: dpdx_i(2), dpdx_i_minus_half(2), dpdx_i_minus_1(2)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupui_amr(child, grid, direction, dt)
                exit
             end if
          end do
       else
          
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          !thaw - edgecell is correct.
          if (.not.thisoctal%edgecell(subcell)) then !xxx changed fromghostcell
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues2(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd,  xnext, px, py, pz,rm1,um1, pm1, qViscosity)

             rho_i_minus_1(1:2) = rho(1:2)
             do j = 1, 2
                rhou_i_minus_1(j) = direction.dot.VECTOR(rhou(j),rhov(j),rhow(j))
             enddo
             pressure_i_minus_1(1:2) = pressure(1:2)
             x_i_minus_1 = px

             rho_i = thisOctal%rho(subcell)
             pressure_i = thisOctal%pressure_i(subcell)
             x_i = thisOctal%x_i(subcell)

             rho_i_minus_2 = rm1
             pressure_i_minus_2 = pm1
             x_i_minus_2 = xnext

             
             x_interface = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistancescale/2.d0

             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues2(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd,  xnext, px, py, pz,rm1,um1, pm1, qViscosity)

             rho_i_plus_1(1:2) = rho(1:2)
             do j = 1, 2
                rhou_i_plus_1(j) = direction.dot.VECTOR(rhou(j),rhov(j),rhow(j))
             enddo
             pressure_i_plus_1(1:2) = pressure(1:2)
             x_i_plus_1 = px
             x_interface_i_p_half = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistancescale/2.d0



             if (thisOctal%x_i(subcell) == thisOctal%x_i_minus_1(subcell)) then
                write(*,*) "x_i bug ", thisOctal%x_i(subcell), thisOctal%x_i_minus_1(subcell)
                print *, "cen ", subcellCentre(thisOctal, subcell)
                if (.not.octalonthread(neighbouroctal, neighboursubcell, myrankGlobal)) then
                   print *, "AND IT WASNT ON THREAD"
                end if
!                print *, "EDGE", thisOctal%edgecell(subcell)
 !               print *, "GHOST", thisOctal%ghostcell(subcell)
             else
 !               print *, "good at u_i"
             endif
             weight = 1.d0 - (thisOctal%x_i(subcell) - x_interface) / (thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell))

             if ((any(rho_i_minus_1(1:2) == 0.d0)).or.(thisOctal%rho(subcell) == 0.d0)) then
                write(*,*) "bug in setupui ",rho_i_minus_1, thisOctal%rho(subcell)
             endif
             !This is the velocity for the i'th cell at i-1/2
             thisRhou = direction.dot.VECTOR(thisOctal%rhou(subcell), thisOctal%rhov(subcell), thisOctal%rhow(subcell))
             thisoctal%u_amr_interface(subcell,1:2) = &
                  weight*thisRhou/thisoctal%rho(subcell) + &
                  (1.d0-weight)*rhou_i_minus_1(1:2)/rho_i_minus_1(1:2)
             if (.not.associated(thisOctal%u_i)) allocate(thisOctal%u_i(1:thisOctal%maxChildren))
             thisOctal%u_i(subcell) = thisRhou / thisOctal%rho(subcell)


! now the cell at i+1/2

             weight = 1.d0 - (thisOctal%x_i_plus_1(subcell) - x_interface_i_p_half) / (thisOctal%x_i_plus_1(subcell) - &
                  thisOctal%x_i(subcell))

             if (.not.associated(thisOctal%u_interface_i_plus_1)) allocate(thisOctal%u_interface_i_plus_1(1:thisOctal%maxChildren))
             thisoctal%u_amr_interface_i_plus_1(subcell,1:2) = &
                  weight * rhou_i_plus_1(1:2)/rho_i_plus_1(1:2)+ &
                  (1.d0-weight)*thisRhou/thisoctal%rho(subcell)

             ! for info, see http://ses.library.usyd.edu.au/bitstream/2123/376/2/adt-NU20010730.12021505_chapter_4.pdf
             if(rhiechow .and. .not. thisOctal%ghostcell(subcell)) then

                locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*&
                     grid%halfsmallestsubcell)
                neighbouroctal => thisoctal
                call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
                call getneighbourvalues2(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, &
                     xnext, px, py, pz, rm1, um1, pm1, qViscosity)


                rho_i_plus_1(1:2) = rho(1:2)
                do j = 1, 2
                   rhou_i_plus_1(j) = direction.dot.VECTOR(rhou(j),rhov(j),rhow(j))
                enddo
                pressure_i_plus_1(1:2) = pressure(1:2)
                x_i_plus_1 = px

!                !now we have i,i+1,i-1 and i-2 can apply rhie-chow
!                !you basically subtract the local pressure gradient and add the average ambient one
!
                u_i = thisRhou / thisOctal%rho(subcell)

                u_i_minus_1(1:2) = rhou_i_minus_1(1:2) / rho_i_minus_1(1:2)

                dx = thisOctal%subcellSize * gridDistanceScale

                dpdx_i_minus_half(1:2) = (pressure_i(1:2) - pressure_i_minus_1(1:2)) / dx
                dpdx_i(1:2) = (pressure_i_plus_1(1:2) - pressure_i_minus_1(1:2)) / (2.d0*dx)
                dpdx_i_minus_1(1:2) = (pressure_i(1:2) - pressure_i_minus_2(1:2)) / (2.d0*dx)


                if(dpdx_i_minus_1(1) /= 0.d0 .and. dpdx_i_minus_half(1) /= 0.d0 .and. &
                     dpdx_i(1) /= 0.d0.and.dpdx_i_minus_1(2) /= 0.d0 .and. dpdx_i_minus_half(2) /= 0.d0 .and. &
                     dpdx_i(2) /= 0.d0) then

!!the approximation screws up if applied to cells that are symmetric on one side only
!!i.e. in a model with zero gradient boundaries you get preferential flow in in - direction

                   thisOctal%u_amr_interface(subcell,1:2) = 0.5d0*(u_i(1:2) + u_i_minus_1(1:2)) - &
                        (dt/(0.5d0*(rho_i_minus_1(1:2) + rho_i)))*dpdx_i_minus_half(1:2) + &
                        0.5d0 * ((dt/rho_i)*dpdx_i(1:2) + (dt/rho_i_minus_1(1:2))*dpdx_i_minus_1(1:2))

                end if

          endif

             
          end if
       endif
    enddo

  end subroutine setupui_amr

!set up neighbour densities and cell interface advecting velocities - x direction
  recursive subroutine setupui_amr4(thisoctal, grid, direction, dt)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    real(double) :: rho(4), rhou(4), rhov(4), rhow(4), q(4), qnext, x, rhoe(4), pressure(4), flux(4), phi, phigas, qViscosity(3,3)
    integer :: subcell, i, neighboursubcell, nd, j
    type(vector) :: direction, locator!, rotator
    real(double) :: rhou_i_minus_1(4), rho_i_minus_1(4), dt
    real(double) :: rhou_i_plus_1(4), x_interface_i_p_half
    real(double) :: xnext, weight, px, py, pz, x_interface
    real(double) :: pressure_i_plus_1(4), pressure_i(4), pressure_i_minus_1(4), pressure_i_minus_2(4)
    real(double) :: rho_i, rho_i_minus_2, rho_i_plus_1(4)
    real(double) :: rm1, um1, pm1
    real(double) :: x_i, x_i_plus_1, x_i_minus_1, x_i_minus_2
    real(double) :: thisRhou
    real(double) :: dx
    real(double) :: u_i(4), u_i_minus_1(4)
    real(double) :: dpdx_i(4), dpdx_i_minus_half(4), dpdx_i_minus_1(4)

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupui_amr4(child, grid, direction, dt)
                exit
             end if
          end do
       else
          
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          !thaw - edgecell is correct.
          if (.not.thisoctal%edgecell(subcell)) then !xxx changed fromghostcell
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues4(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd,  xnext, px, py, pz,rm1,um1, pm1, qViscosity)

             rho_i_minus_1(1:4) = rho(1:4)
             do j = 1, 4
                rhou_i_minus_1(j) = direction.dot.VECTOR(rhou(j),rhov(j),rhow(j))
             enddo
             pressure_i_minus_1(1:4) = pressure(1:4)
             x_i_minus_1 = px

             rho_i = thisOctal%rho(subcell)
             pressure_i = thisOctal%pressure_i(subcell)
             x_i = thisOctal%x_i(subcell)

             rho_i_minus_2 = rm1
             pressure_i_minus_2 = pm1
             x_i_minus_2 = xnext

             
             x_interface = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistancescale/2.d0

             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues4(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz,rm1,um1, pm1, qViscosity)

             rho_i_plus_1(1:4) = rho(1:4)
             do j = 1, 4
                rhou_i_plus_1(j) = direction.dot.VECTOR(rhou(j),rhov(j),rhow(j))
             enddo
             pressure_i_plus_1(1:4) = pressure(1:4)
             x_i_plus_1 = px
             x_interface_i_p_half = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistancescale/2.d0



             if (thisOctal%x_i(subcell) == thisOctal%x_i_minus_1(subcell)) then
                write(*,*) "x_i bug ", thisOctal%x_i(subcell), thisOctal%x_i_minus_1(subcell)
                print *, "cen ", subcellCentre(thisOctal, subcell)
                if (.not.octalonthread(neighbouroctal, neighboursubcell, myrankGlobal)) then
                   print *, "AND IT WASNT ON THREAD"
                end if
!                print *, "EDGE", thisOctal%edgecell(subcell)
 !               print *, "GHOST", thisOctal%ghostcell(subcell)
             else
 !               print *, "good at u_i"
             endif
             weight = 1.d0 - (thisOctal%x_i(subcell) - x_interface) / (thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell))

             if ((any(rho_i_minus_1(1:4) == 0.d0)).or.(thisOctal%rho(subcell) == 0.d0)) then
                write(*,*) "bug in setupui ",rho_i_minus_1, thisOctal%rho(subcell)
             endif
             !This is the velocity for the i'th cell at i-1/2
             thisRhou = direction.dot.VECTOR(thisOctal%rhou(subcell), thisOctal%rhov(subcell), thisOctal%rhow(subcell))
             thisoctal%u_amr_interface(subcell,1:4) = &
                  weight*thisRhou/thisoctal%rho(subcell) + &
                  (1.d0-weight)*rhou_i_minus_1(1:4)/rho_i_minus_1(1:4)

             if (ANY(thisOctal%u_amr_interface(subcell,:) > 3d7)) then
                write(*,*) "u interface warning ",thisOctal%u_amr_interface(subcell,1:4)/1.d5
                write(*,*) "thisrhou ",thisRhou, " thisrho ",thisOctal%rho(subcell), " this u ", &
                     weight*thisRhou/thisoctal%rho(subcell)
                write(*,*) "rhou_i_minus_1 ",rhou_i_minus_1(1:4), " rho_i_minus_1 ",rho_i_minus_1(1:4), &
                     "u_i_minus_1 ",rhou_i_minus_1(1:4)/rho_i_minus_1(1:4)
             endif


             if (.not.associated(thisOctal%u_i)) allocate(thisOctal%u_i(1:thisOctal%maxChildren))
             thisOctal%u_i(subcell) = thisRhou / thisOctal%rho(subcell)


! now the cell at i+1/2

             weight = 1.d0 - (thisOctal%x_i_plus_1(subcell) - x_interface_i_p_half) / (thisOctal%x_i_plus_1(subcell) - &
                  thisOctal%x_i(subcell))

             if (.not.associated(thisOctal%u_interface_i_plus_1)) allocate(thisOctal%u_interface_i_plus_1(1:thisOctal%maxChildren))
             thisoctal%u_amr_interface_i_plus_1(subcell,1:4) = &
                  weight * rhou_i_plus_1(1:4)/rho_i_plus_1(1:4)+ &
                  (1.d0-weight)*thisRhou/thisoctal%rho(subcell)

             ! for info, see http://ses.library.usyd.edu.au/bitstream/2123/376/2/adt-NU20010730.12021505_chapter_4.pdf
             if(rhiechow .and. .not. thisOctal%ghostcell(subcell)) then

                locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*&
                     grid%halfsmallestsubcell)
                neighbouroctal => thisoctal
                call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
                call getneighbourvalues4(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, &
                     xnext, px, py, pz, rm1, um1, pm1, qViscosity)


                rho_i_plus_1(1:4) = rho(1:4)
                do j = 1, 4
                   rhou_i_plus_1(j) = direction.dot.VECTOR(rhou(j),rhov(j),rhow(j))
                enddo
                pressure_i_plus_1(1:4) = pressure(1:4)
                x_i_plus_1 = px

!                !now we have i,i+1,i-1 and i-2 can apply rhie-chow
!                !you basically subtract the local pressure gradient and add the average ambient one
!
                u_i = thisRhou / thisOctal%rho(subcell)

                u_i_minus_1(1:4) = rhou_i_minus_1(1:4) / rho_i_minus_1(1:4)

                dx = thisOctal%subcellSize * gridDistanceScale

                dpdx_i_minus_half(1:4) = (pressure_i(1:4) - pressure_i_minus_1(1:4)) / dx
                dpdx_i(1:4) = (pressure_i_plus_1(1:4) - pressure_i_minus_1(1:4)) / (2.d0*dx)
                dpdx_i_minus_1(1:4) = (pressure_i(1:4) - pressure_i_minus_2(1:4)) / (2.d0*dx)


                if(dpdx_i_minus_1(1) /= 0.d0 .and. dpdx_i_minus_half(1) /= 0.d0 .and. &
                     dpdx_i(1) /= 0.d0.and.dpdx_i_minus_1(2) /= 0.d0 .and. dpdx_i_minus_half(2) /= 0.d0 .and. &
                     dpdx_i(2) /= 0.d0) then

!!the approximation screws up if applied to cells that are symmetric on one side only
!!i.e. in a model with zero gradient boundaries you get preferential flow in in - direction

                   thisOctal%u_amr_interface(subcell,1:4) = 0.5d0*(u_i(1:4) + u_i_minus_1(1:4)) - &
                        (dt/(0.5d0*(rho_i_minus_1(1:4) + rho_i)))*dpdx_i_minus_half(1:4) + &
                        0.5d0 * ((dt/rho_i)*dpdx_i(1:4) + (dt/rho_i_minus_1(1:4))*dpdx_i_minus_1(1:4))

                end if

          endif

             
          end if
       endif
    enddo

  end subroutine setupui_amr4


!setup the flux at i+1/2 and if necessary modify flux at i-1/2
!note that thisOctal%flux_i is the flux at i-1/2
  recursive subroutine setupflux(thisoctal, grid, direction)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: nd,nc
    real(double) :: xnext, fac, px, py, pz, rm1, um1, pm1


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

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

   
             thisoctal%flux_i_plus_1(subcell) = flux

             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%flux_i_minus_1(subcell) = flux
             


!If this cell is finer than the cell from which material is streaming, then it is worth doing a coarse to fine flux interpolation

             if(thisOctal%nDepth > nd .and. fluxinterp .and. .not. thisOctal%ghostCell(subcell)) then
                if(.not. thisOctal%oneD) then
                   call NormalFluxGradient(thisOctal, subcell, neighbourOctal, neighbourSubcell, grid, direction, fac)
                   else
                      fac = 0.d0
                   end if
             else
                fac = 0.d0
             end if   

!Following the interpolation, modify the incoming flux by the appropriate factor
             thisOctal%flux_i(subcell) = thisOctal%flux_i(subcell) + fac

          endif
       endif
    enddo
  end subroutine setupflux

!setup the flux at i+1/2 and if necessary modify flux at i-1/2
!note that thisOctal%flux_i is the flux at i-1/2
  recursive subroutine setupfluxCylindrical(thisoctal, grid, direction)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: nd, nc
    real(double) :: xnext, fac, px, py, pz, rm1, um1, pm1


    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupfluxCylindrical(child, grid, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

   

             
             thisoctal%flux_i_plus_1(subcell) = flux


!             if (nd < thisOctal%nDepth) then ! neighbouring cell is coarser, use local u_interface for flux
!                if (thisoctal%u_interface_i_plus_1(subcell).ge.0.d0) then
!                   thisoctal%flux_i_plus_1(subcell) = thisoctal%u_interface_i_plus_1(subcell) * thisoctal%q_i(subcell)
!                else
!                   thisoctal%flux_i_plus_1(subcell) = thisoctal%u_interface_i_plus_1(subcell) * thisoctal%q_i_plus_1(subcell)
!                endif
!             endif


!             if (inSubcell(thisOctal, subcell, VECTOR(50.d3, 0.d0, 1.07d6)).and.(direction%z > 0.9d0)) then
!                write(*,*) "flux ", -flux, thisOctal%flux_i_plus_1(subcell), thisOctal%u_interface_i_plus_1(subcell)/1.d5
!                write(*,*) "neighbour u_i ",neighbourOctal%u_interface(neighboursubcell)/1.d5
!                write(*,*) "neighbour q_i ",neighbourOctal%q_i(neighboursubcell)
!                write(*,*) "neighbour flux ",neighbourOctal%q_i(neighboursubcell)*neighbourOctal%u_interface(neighboursubcell)
!                write(*,*) "q_i_p_1 ", thisOctal%q_i_plus_1(subcell)
!             endif



             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz,rm1, um1, pm1, qViscosity)
             thisoctal%flux_i_minus_1(subcell) = flux
             

!If this cell is finer than the cell from which material is streaming, then it is worth doing a coarse to fine flux interpolation

             if(thisOctal%nDepth > nd .and. fluxinterp .and. .not. thisOctal%ghostCell(subcell)) then
                if(.not. thisOctal%oneD) then
                   call NormalFluxGradient(thisOctal, subcell, neighbourOctal, neighbourSubcell, grid, direction, fac)
                   else
                      fac = 0.d0
                   end if
             else
                fac = 0.d0
             end if   

!Following the interpolation, modify the incoming flux by the appropriate factor
             thisOctal%flux_i(subcell) = thisOctal%flux_i(subcell) + fac

          endif
       endif
    enddo
  end subroutine setupfluxCylindrical

!setup the flux at i+1/2 and if necessary modify flux at i-1/2
!note that thisOctal%flux_i is the flux at i-1/2
  recursive subroutine setupfluxCylindrical2(thisoctal, grid, direction)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rho(2), rhoe(2), rhou(2), rhov(2), rhow(2), x, q(2), qnext, pressure(2), flux(2), phi, phigas,qViscosity(3,3)
    integer :: nd
    real(double) :: xnext, px, py, pz, rm1, um1, pm1


    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupfluxCylindrical2(child, grid, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues2(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd,  xnext, px, py, pz, rm1,um1, pm1, qViscosity)

             
             thisoctal%flux_amr_i_plus_1(subcell,1:2) = flux(1:2)

             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues2(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz,rm1, um1, pm1, qViscosity)

             thisoctal%flux_amr_i_minus_1(subcell,1:2) = flux(1:2)
             

          endif
       endif
    enddo
  end subroutine setupfluxCylindrical2

!setup the flux at i+1/2 and if necessary modify flux at i-1/2
!note that thisOctal%flux_i is the flux at i-1/2
  recursive subroutine setupflux_amr(thisoctal, grid, direction)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rho(4), rhoe(4), rhou(4), rhov(4), rhow(4), x, q(4), qnext, pressure(4), flux(4), phi, phigas,qViscosity(3,3)
    integer :: nd
    real(double) :: xnext, px, py, pz, rm1, um1, pm1


    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupflux_amr(child, grid, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues4(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

             
             thisoctal%flux_amr_i_plus_1(subcell,1:4) = flux(1:4)

             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues4(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz,rm1, um1, pm1, qViscosity)

             thisoctal%flux_amr_i_minus_1(subcell,1:4) = flux(1:4)
             

          endif
       endif
    enddo
  end subroutine setupflux_amr


!setup the flux at i+1/2 and if necessary modify flux at i-1/2
!note that thisOctal%flux_i is the flux at i-1/2
  recursive subroutine setupfluxSpherical(thisoctal, grid, direction)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: nd, nc
    real(double) :: xnext, fac, px, py, pz, rm1, um1, pm1


    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupfluxSpherical(child, grid, direction)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

   
             thisoctal%flux_i_plus_1(subcell) = flux

             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz,rm1, um1, pm1, qViscosity)
             thisoctal%flux_i_minus_1(subcell) = flux
             

!If this cell is finer than the cell from which material is streaming, then it is worth doing a coarse to fine flux interpolation

             if(thisOctal%nDepth > nd .and. fluxinterp .and. .not. thisOctal%ghostCell(subcell)) then
                if(.not. thisOctal%oneD) then
                   call NormalFluxGradient(thisOctal, subcell, neighbourOctal, neighbourSubcell, grid, direction, fac)
                   else
                      fac = 0.d0
                   end if
             else
                fac = 0.d0
             end if   

!Following the interpolation, modify the incoming flux by the appropriate factor
             thisOctal%flux_i(subcell) = thisOctal%flux_i(subcell) + fac

          endif
       endif
    enddo
  end subroutine setupfluxSpherical


!Flux interpolation routine                                                                                                                                                                                                       
!Determine the factor by which to modify the incoming flux in a cell
  subroutine normalFluxGradient(thisOctal, subcell, neighbourOctal, neighbourSubcell, grid, direction, fac)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer :: thisoctal
    type(octal), pointer :: neighbourOctal
    integer :: subcell, neighbourSubcell
    type(vector) :: direction
    type(vector) :: rVec
    type(vector), allocatable :: community(:), communitySubset(:)
    real(double), intent(out) :: fac
    real(double), allocatable ::  xpos(:), f(:)
    integer :: iTot
    real(double) :: m, dx, df
!    logical, intent(in) :: fineToCoarse
    logical :: upper, right
    real(double) :: fac_a, fac_b, m_a, m_b, df_a, df_b !3D parameters




    fac = 0.d0
    if(neighbourOctal%twoD) then

       iTot = 2
       allocate(community(iTot))
       allocate(communitySubset(2))
       allocate(xpos(iTot))
       allocate(f(iTot))

!The community cells are those perpendicular to the advection direction, they are not the neighbour cells.

       !Get the direction of the neighbour's community cells                                                                                                                                                                             
       if(abs(direction%x) == 1.d0) then !Advecting in ±x-direction                                                                                                                                                                      
          community(1) = VECTOR(0.d0, 0.d0,1.d0)
          community(2) = VECTOR(0.d0, 0.d0,-1.d0)
       else if(abs(direction%z) == 1.d0) then !Advecting in ±z direction                                                                                                                                                                 
          community(1) = VECTOR(1.d0, 0.d0,0.d0)
          community(2) = VECTOR(-1.d0, 0.d0,0.d0)
       else
          print *, "Advecting in unknown direction!"
          stop
       end if


!Get the perpendicular flux gradient
       call normalFluxGradientSingle(f, xpos, direction, community, upper, thisOctal, subcell, neighbourOctal, &
            neighbourSubcell, grid)
       
!Calculate the community flux gradient
       m = (f(1) - f(2))/(xpos(1) - xpos(2))

       dx = thisOctal%subcellSize*griddistancescale

       !Flux variation for coarse cell centre to fine cell centre, perpendicular to advection direction                                                                                                                                 
        df = abs(m*dx)/2.d0

       rVec = subcellCentre(thisOctal, subcell)!                                                                                                                                                                                         

!Give the modification factor the appropriate sign
       if(upper .and. m > 0.d0) then
          fac = df
       else if(upper .and. m < 0.d0) then
          fac = -df
       else if (.not. upper .and. m > 0.d0) then
          fac = -df
       else if(.not. upper .and. m < 0.d0) then
          fac = df
       else
          fac = 0.d0
       end if
    else !3D case                                                                                                                                                                                                                        

       fac = 0.d0

       iTot = 4
       allocate(community(iTot))
       allocate(communitySubset(2))
       allocate(xpos(2))
       allocate(f(2))
       !Get the direction of the neighbour's community cells
       if(abs(direction%x) == 1.d0) then !Advecting in ±x-direction
          community(1) = VECTOR(0.d0, 0.d0 ,1.d0)
          community(2) = VECTOR(0.d0, 0.d0, -1.d0)
          community(3) = VECTOR(0.d0, 1.d0, 0.d0)
          community(4) = VECTOR(0.d0, -1.d0, 0.d0)

       else if(abs(direction%z) == 1.d0) then !Advecting in ±z direction
          community(1) = VECTOR(1.d0, 0.d0, 0.d0)
          community(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          community(3) = VECTOR(0.d0, 1.d0, 0.d0)
          community(4) = VECTOR(0.d0, -1.d0, 0.d0)

       else if(abs(direction%y) == (1.d0)) then !advecting in the ±y direction
          community(1) = VECTOR(1.d0, 0.d0, 0.d0)
          community(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          community(3) = VECTOR(0.d0, 0.d0, 1.d0)
          community(4) = VECTOR(0.d0, 0.d0, -1.d0)

       else
          print *, "Advecting in unknown direction!"
          stop
       end if


!First get the community gradient in one direction
       communitySubset(1) = community(1)
       communitySubset(2) = community(2)

       call normalFluxGradientSingle(f, xpos, direction, community, upper, thisOctal, subcell, neighbourOctal, &
            neighbourSubcell, grid)
!Calculate the community gradient
       m_a = (f(1) - f(2))/(xpos(1) - xpos(2))


!Now repeat for the other direction
       communitySubset(1) = community(3)
       communitySubset(2) = community(4)

       right = .false.

       call normalFluxGradientSingle(f, xpos, direction, community, upper, thisOctal, subcell, neighbourOctal, &
            neighbourSubcell, grid)
       if(upper) right = .true.
       
       m_b = (f(1) - f(2))/(xpos(1) - xpos(2))

       dx = thisOctal%subcellSize*griddistancescale

       !Flux variation for coarse cell centre to fine cell centre, perpendicular to advection direction
        df_a = abs(m_a*dx)/2.d0
        df_b = abs(m_b*dx)/2.d0

       rVec = subcellCentre(thisOctal, subcell)

!Get the modification due to gradient in one direction
       if(upper .and. m_a > 0.d0) then
          fac_a = df_a
       else if(upper .and. m_a < 0.d0) then
          fac_a = -df_a
       else if (.not. upper .and. m_a > 0.d0) then
          fac_a = -df_a
       else if(.not. upper .and. m_a < 0.d0) then
          fac_a = df_a
       else
          fac_a = 0.d0
       end if

!Get the modification due to the other direction
       if(right .and. m_b > 0.d0) then
          fac_b = df_b
       else if(right .and. m_b < 0.d0) then
          fac_b = -df_b
       else if (.not. right .and. m_b > 0.d0) then
          fac_b = -df_b
       else if(.not. right .and. m_b < 0.d0) then
          fac_b = df_b
       else
          fac_b = 0.d0
       end if

!Total modifying factor is simply the sum of the factors for each perpendicular direction
       fac = fac_a + fac_b
       

    end if

    deallocate(community, communitySubset)
    deallocate(xpos)
    deallocate(f)

  end subroutine normalFluxGradient

!Find the flux gradient perpendicular to the advection direction
  subroutine   normalFluxGradientSingle(f, xpos, direction, community, upper, thisOctal, subcell, neighbourOctal, &
       neighbourSubcell, grid)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer :: thisoctal
    type(octal), pointer :: neighbourOctal
    type(octal), pointer :: communityOctal
    type(octal), pointer :: faceOctal
    type(octal), pointer :: mirrorOctal
    integer :: subcell, neighbourSubcell, communitySubcell, mirrorSubcell, faceSubcell
    type(vector) :: direction
    type(vector) :: rVec, cVec, mVec
    type(vector) :: community(:)
    type(vector) :: locator, probe
    type(vector) :: nVec
    real(double) ::  xpos(:), f(:)
    integer :: nCornerBound
    real(double) :: rho, rhoe, rhou, rhov, rhow, x, q, qnext, pressure, flux, phi, phigas,qViscosity(3,3)
    
    integer :: nd, i, ID(2), nc
    real(double) :: xnext, dx, px, py, pz, rm1, um1, pm1
    logical :: upper
    logical :: fineToCoarse = .false.

    ID = 0

    call getneighbourvalues(grid, thisOctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, &
         rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1 ,qViscosity)
         
    nVec = VECTOR(px, py, pz)

    
    do i = 1, 2
      
       locator = subcellCentre(thisOctal, subcell) + community(i) * &
            (thisOctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
       
       if(inOctal(grid%octreeRoot, locator)) then
          if(octalOnThread(neighbouroctal, neighboursubcell, myRankGlobal) .or. (.not. fineToCoarse)) then
             
             if(fineToCoarse) then
                communityoctal => neighbouroctal
                
                call findsubcelllocal(locator, communityoctal, communitysubcell)
                
                if(octalOnTHread(communityOctal, communitySubcell, myRankGlobal)) then
                   call getneighbourvalues(grid, neighbouroctal, neighboursubcell, communityoctal, communitysubcell, &
                        direction, q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, &
                        px, py, pz, rm1,um1, pm1, qViscosity)
                else
                   call getneighbourvalues(grid, neighbouroctal, neighboursubcell, communityoctal, communitysubcell, &
                        community(i), q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, &
                        xnext, px, py, pz,  rm1,um1, pm1, qViscosity)
                end if
                
                ID(i) = 1
                
             else
                !check we got the sufficiently far away cell
                communityoctal => thisoctal
                call findsubcelllocal(locator, communityoctal, communitysubcell)
                
                !Get first try values
                if(octalOnTHread(communityOctal, communitySubcell, myRankGlobal)) then
                   call getneighbourvalues(grid, thisoctal, subcell, communityoctal, communitySubcell, direction, q, &
                        rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, &
                        rm1,um1, pm1, qViscosity)
                else
                   call getneighbourvalues(grid, thisoctal, subcell, communityoctal, communitySubcell, community(i), q, &
                        rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, &
                        rm1,um1, pm1, qViscosity)
                end if
                
                rVec = VECTOR(px, py , pz)
                
                if(i == 1) then
                   upper = .true.
                else 
                   upper = .false.
                end if
                
                if(octalOnthread(communityoctal, communitysubcell, myRankGlobal) .and. (nd-thisOctal%nDepth == 0)) then
                   probe = subcellCentre(communityoctal, communitysubcell) - direction*(thisOctal%subcellSize/2.d0 &
                        + 0.01d0*grid%halfsmallestsubcell)
                   
                   mirrorOctal=>thisOctal
                   call findsubcelllocal(probe,mirrorOctal, mirrorSubcell)
                   
                   
                   if(octalOnTHread(mirrorOctal, mirrorSubcell, myRankGlobal)) then
                      Call getneighbourvalues(grid, communityOctal, communitySubcell, mirrorOctal, mirrorsubcell, direction, q, &
                           rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, &
                           nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                   else
                      Call getneighbourvalues(grid, communityOctal, communitySubcell, mirrorOctal, mirrorsubcell, community(i), &
                           q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, &
                           nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                   end if
                   
                   mVec = VECTOR(px, py, pz)
                   
                   if(mVec == nVec) then
                      !Found the neighbour again, need to move one more cell
                      
                      if(i == 1) then
                         upper = .false.
                      else 
                         upper = .true.
                      end if
                      dx = thisOctal%subcellSize*griddistancescale
                      
                      faceOctal=>communityOctal
                      
                      locator = rVec + community(i) * &
                           (thisOctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
                      
                      call findsubcelllocal(locator, faceOctal, faceSubcell)
                      
                      if(octalOnThread(faceOctal, faceSubcell, myRankGlobal)) then
                         call getneighbourvalues(grid, communityOctal, communitySubcell, faceOctal, facesubcell, direction, q, &
                              rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, &
                              nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                      else
                         call getneighbourvalues(grid, communityOctal, communitySubcell, faceOctal, facesubcell, community(i), &
                              q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, &
                              nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                      end if
                      
                      ID(i) = 2
                   else
                      if(i == 1) then
                         upper = .true.
                      else 
                         upper = .false.
                      end if
                      
                      call getneighbourvalues(grid, thisoctal, subcell, communityoctal, communitysubcell, direction, q, &
                           rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, &
                           rm1,um1, pm1, qViscosity)
                      ID(i) = 3
                   end if
                else
                   ID(i) = 1
                   end if
                end if
                rVec = VECTOR(px, py, pz)
                f(i) = flux


                if(abs(community(i)%z) == 1.d0) then
                   xpos(i) = rVec%z
                else if(abs(community(i)%x) == 1.d0) then
                   xpos(i) = rVec%x
                else if(abs(community(i)%y) == 1.d0) then
                   xpos(i) = rVec%y
                else
                   print *, "unrecognized direction!"
                   stop
                end if
                
             else
                
                !If the neighbour is not on the same thread then it will not hold the community values
                !Need to scan along the MPI boundary and probe into the next domain until the community cell is found
                !The community cell values can then be obtained from its the boundary partner on the native side of the domain. 

                if(.not. fineToCoarse) then
                   call torus_abort("Coarse to fine problem in flux interpolation")
                end if

                mirrorOctal => thisOctal
                !Move along the native side of the domain by 1 cell in the appropriate direction
                
                locator = subcellcentre(thisoctal, subcell) + community(i) * &
                (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)      

                call findsubcelllocal(locator, mirrorOctal, mirrorSubcell)


                if(octalOnThread(mirrorOctal,mirrorSubcell, myRankGlobal)) then
!Probe across the mpi boundary and check if the corresponding cell is the community subcell
                   locator = subcellcentre(mirrorOctal, mirrorSubcell) + direction * & 
                   (mirroroctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)  
                   
                   communityoctal => mirrorOctal
                   call findsubcelllocal(locator, communityOctal, communitySubcell)
      
                   call getneighbourvalues(grid, mirroroctal, mirrorsubcell, communityoctal, communitysubcell, direction, q, rho, & 
                        rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, &
                        pm1, qViscosity)
                   ID(i) = 2

                else
                   nCornerBound = getNCornerBoundFromDirection(direction, i)
                   if(.not. associated(thisOctal%mpiCornerStorage)) then
                      print *, direction
                      print *, community
                      print *, "centre ", subcellCentre(thisOctal, subcell)
                      print *, "nCornerBound{",i,")", nCornerBound
                      call torus_abort("Corner storage not allocated when it should have been")
                      
                   end if
                   ID(i) = 3
                   px = thisOctal%mpiCornerStorage(subcell, nCornerBound, 1)
                   py = thisOctal%mpiCornerStorage(subcell, nCornerBound, 2)
                   pz = thisOctal%mpiCornerStorage(subcell, nCornerBound, 3)
                   flux = thisOctal%mpiCornerStorage(subcell, nCornerBound, 4)
                end if

                cvec = VECTOR(px, py, pz)

                if(cVec%x /= nVec%x .or. cVec%z /= nVec%z) then
                   !Found the community subcell

!                   !Get the values ...
!                   call getneighbourvalues(grid, mirroroctal, mirrorsubcell, communityoctal, communitysubcell, direction, q, &
!                   rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz) 
                  
!                   rVec = subcellcentre(communityoctal, communitysubcell)
                   f(i) = flux

                   if(abs(community(i)%z) == 1.d0) then
                      xpos(i) = rVec%z
                   else if(abs(community(i)%x) == 1.d0) then
                      xpos(i) = rVec%x
                   else if(abs(community(i)%y) == 1.d0) then
                      xpos(i) = rVec%y
                   else
                      print *, "unrecognized direction!"
                      stop
                   end if
                   
                else
                   !Found the same cell
                   !Move along one more cell width

                   !THAW - pretty sure that it is this section that will be the origin of any woes

                   if(.not. fineToCoarse) then
                      call torus_abort("fine to coarse mistaken for coarse to fine")
                   end if

                   locator = subcellcentre(mirrorOctal, mirrorSubcell) + community(i) * &
                   (mirroroctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)  
                  
                   faceOctal=>mirrorOctal

                   call findsubcelllocal(locator, faceOctal, faceSubcell) 

!                   if(octalOnThread(faceOctal, faceSubcell, myRankGlobal)) then
                      call getneighbourvalues(grid, mirroroctal, mirrorsubcell, faceoctal, facesubcell, direction, q, rho, &
                           rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, &
                           qViscosity)
!                   else
!                      call getneighbourvalues(grid, mirroroctal, mirrorsubcell, faceoctal, facesubcell, community(i), q, rho, &
!                      rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz)
!                   End if

!Probe across the mpi boundary and check if the corresponding cell is the community subcell
                   locator = VECTOR(px, py, pz) + direction * & 
                      (mirroroctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)  

                   call findsubcelllocal(locator, communityOctal, communitySubcell)
                   call getneighbourvalues(grid, faceoctal, facesubcell, communityoctal, communitysubcell, direction, q, rho, &
                        rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, &
                        qViscosity)
                   
                   cvec = VECTOR(px, py, pz)
                
                   if(cVec%x /= nVec%x .or. cVec%z /= nVec%z) then
                    !Found the community subcell 
                    !Get the values ...  
!                      call getneighbourvalues(grid, mirroroctal, mirrorsubcell, communityoctal, communitysubcell, direction, q, &
!                      rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz)  
!                      
!                      rVec = subcellcentre(communityoctal, communitysubcell) 


                      ID(i) = 4
                      f(i) = flux


                      if(abs(community(i)%z) == 1.d0) then
                         xpos(i) = rVec%z
                      else if(abs(community(i)%x) == 1.d0) then
                         xpos(i) = rVec%x
                      else if(abs(community(i)%y) == 1.d0) then
                         xpos(i) = rVec%y
                      else
                         print *, "unrecognized direction!"
                         stop
                      end if


                   else
                      print *, "-----------------------"
                      print *, "cvec", cvec, myrankGlobal
                      print *, "nvec", nvec
                      print *, "subcellCentre(mirrorOctal, mirrorSubcell)", subcellCentre(mirrorOctal, mirrorSubcell)
                      print *, "-----------------------"
                      call torus_abort("Warning, Neighbour cells differ in refinement by more than one depth!")                     
                   end if
                end if
             end if
          end if
       end do
  end subroutine normalFluxGradientSingle


  logical function refineShock(thisOctal, subcell, grid)
    type(octal), pointer :: thisOCtal
    type(octal), pointer :: neighbourOctal
    type(gridtype) :: grid
    integer :: neighbourSubcell
    integer :: subcell
    type(vector) :: locator, dirVec(4), rVec
    integer :: nd, i, j, index
    real(double) :: rho, rhoe, rhou, rhov, rhow, energy, phi,x,y,z, pressure
    real(double), allocatable :: p(:), xArray(:)
    integer :: nDimensions
    real(double) :: m2, m4, Sp, fac

    if(thisOctal%threed) then
       nDimensions = 3
    else if (thisOctal%twoD) then
       nDimensions = 2
    else
       nDimensions = 1
    end if
    allocate(p(nDimensions*4))
    allocate(xArray(nDimensions*4))
    refineShock = .false.

    if(.not. thisOctal%ghostcell(subcell)) then
       refineShock = .false.
       p = 0.d0
       index = 1
       xArray = 0.d0
       p = 0.d0
       do i = 1, nDimensions
          if(i == 1) then
             dirVec(1) = VECTOR(-1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(1.d0, 0.d0, 0.d0)
             dirVec(3) = VECTOR(1.d0, 0.d0, 0.d0)
             dirVec(4) = VECTOR(1.d0, 0.d0, 0.d0)
          else if(i ==2) then
             dirVec(1) = VECTOR(0.d0, 0.d0, -1.d0)
             dirVec(2) = VECTOR(0.d0, 0.d0, 1.d0)
             dirVec(3) = VECTOR(0.d0, 0.d0, 1.d0)
             dirVec(4) = VECTOR(0.d0, 0.d0, 1.d0)
          else 
             dirVec(1) = VECTOR(0.d0, -1.d0, 0.d0)
             dirVec(2) = VECTOR(0.d0, 1.d0, 0.d0)
             dirVec(3) = VECTOR(0.d0, 1.d0, 0.d0)
             dirVec(4) = VECTOR(0.d0, 1.d0, 0.d0)
          end if
          
          do j = 1, 4             
             if(j==1) then
                locator = subcellCentre(thisOctal, subcell) + &
                     ((thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(j))
                if (inOctal(grid%octreeRoot, locator)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                   call getHydroValues(grid, locator, nd, rho, rhoe, rhou, rhov, rhow, energy, phi,x,y,z, pressure, .false.)
                   p(index) = pressure                   
                   if(i == 1) then
                      xArray(index) = x
                   else if(i == 2) then
                      xArray(index) = z
                   else
                      xArray(index) = y
                   end if
                else
                   refineShock = .false.
                   goto 123
!                   call torus_abort("point outside the grid")
                end if
             else if(j==2) then
                p(index)= thisOctal%pressure_i(subcell)
                rVec = subcellCentre(thisOctal, subcell)
                if(i == 1) then
                   xArray(index) = rVec%x
                else if(i == 2) then
                   xArray(index) = rVec%z
                else
                   xArray(index) = rVec%y
                end if
             else if (j==3) then
                locator = subcellCentre(thisOctal, subcell) + &
                     (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(j)
                if (inOctal(grid%octreeRoot, locator)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
!                   locator = subcellCentre(neighbourOctal, neighbourSubcell)
                   call getHydroValues(grid, locator, nd, rho, rhoe, rhou, rhov, rhow, energy, phi,x,y,z, pressure,.false.)   
                   p(index) = pressure
                   if(i == 1) then
                      xArray(index) = x
                   else if(i == 2) then
                      xArray(index) = z
                   else
                      xArray(index) = y
                   end if
                else
                   refineShock = .false.
                   goto 123
!                   call torus_abort("point outside the grid")                   
                end if
             else
                locator = subcellCentre(thisOctal, subcell) + &
                     (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(j)
                !need to advance a second cell
                if (inOctal(grid%octreeRoot, locator)) then
                   neighbourOctal => thisOctal
                   call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                   call getHydroValues(grid, locator, nd, rho, rhoe, rhou, rhov, rhow, energy, phi,x,y,z, pressure, .false.)
                   
                   if(nd > thisOctal%nDepth) then
                      fac = 0.5d0
                   else if(nd == thisOctal%nDepth) then
                      fac = 1.d0
                   else
                      fac = 2.d0
                   end if
                   locator = locator + ((fac*thisOctal%subcellSize)/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(j)
                   if (inOctal(grid%octreeRoot, locator)) then
                      neighbourOctal => thisOctal
                      call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                      call getHydroValues(grid, locator, nd, rho, rhoe, rhou, rhov, rhow, energy, phi,x,y,z, pressure, .false.)
                      p(index) = pressure                   
                      if(i == 1) then
                         xArray(index) = x
                      else if(i == 2) then
                         xArray(index) = z
                      else
                         xArray(index) = y
                      end if
                   else
                      refineShock = .false.
                      goto 123
!                      call torus_abort("point outside the grid")                   
                   end if
                end if
             end if
             index = index +1
          end do
       end do

    
    !Now we have pressures across 4 cells

!Identify shock by comparing pressure gradient over 2 and 4 cells
!see Fryxell et al 2000, ApJ, 131, 273 section 3.1.3
!Note that this is still pretty different to the above 

    do i = 1, ndimensions

       m2= (p(i*4 - 1) - p(i*4 -2))!/(xArray(i*4-1) - xArray(i*4 - 2))
       m4= (p(i*4) - p(i*4 - 3))!/(xArray(i*4) - xArray(i*4 - 3))


       if(m4 /= 0.d0) then
          Sp = m2/m4
       else
          Sp = 0.d0
       end if
       
       if(Sp > (4.d0/3.d0)) then
          refineShock = .true.
          exit
       else
          refineShock = .false.
       end if
    end do
 end if

 deallocate(p, xArray)

123 continue
  end function refineShock

!set up neighbour pressures
  recursive subroutine setuppressure(thisoctal, grid, direction)
    use mpi
    type(gridtype) :: grid
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, qViscosity(3,3)
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    integer :: nd, nc
    real(double) :: xnext, px, py, pz , rm1, um1, pm1


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
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle


          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%pressure_i_plus_1(subcell) = pressure

             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz,rm1,um1, pm1, qViscosity)
             thisoctal%pressure_i_minus_1(subcell) = pressure

          endif
       endif
    enddo
  end subroutine setuppressure

!set up i±1 neighbour central velocities (u_i)   - x direction
!note these are not advecting velocities (u_i-1/2) 
  recursive subroutine setupupm(thisoctal, grid, direction)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer   :: neighbouroctal
    type(octal), pointer  :: child 
    integer :: subcell, i, neighboursubcell
    type(vector) :: direction, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, qViscosity(3,3)
    integer :: nd, nc
    real(double) :: xnext, px, py, pz, rm1, um1, pm1


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
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisoctal%edgecell(subcell)) then
             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%u_i_plus_1(subcell) = (direction.dot.VECTOR(rhou, rhov, rhow))/rho
             
             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             thisoctal%u_i_minus_1(subcell) =  (direction.dot.VECTOR(rhou, rhov, rhow))/rho
             thisOctal%rho_i_minus_1(subcell) = rho
          endif
       endif
    enddo
  end subroutine setupupm


!  recursive subroutine setupupm(thisoctal, grid, direction)
!    use mpi
!    type(gridtype) :: grid
!    type(octal), pointer   :: thisoctal
!    type(octal), pointer   :: neighbouroctal
!    type(octal), pointer  :: child 
!    integer :: subcell, i, neighboursubcell
!    type(vector) :: direction, locator
!    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, q11,q22,q33
!    integer :: nd
!    real(double) :: xnext, px, py, pz
!
!
!    do subcell = 1, thisoctal%maxchildren
!       if (thisoctal%haschild(subcell)) then
!          ! find the child
!          do i = 1, thisoctal%nchildren, 1
!             if (thisoctal%indexchild(i) == subcell) then
!                child => thisoctal%child(i)
!                call setupupm(child, grid, direction)
!                exit
!             end if
!          end do
!       else
!!          if (thisoctal%mpithread(subcell) /= myrank) cycle
!          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
!
!          if (.not.thisoctal%edgecell(subcell)) then
!             locator = subcellcentre(thisoctal, subcell) + direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
!             neighbouroctal => thisoctal
!             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
!             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction, q, rho, rhoe, &
!                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, q11, q22, q33)
!             thisoctal%u_i_plus_1(subcell) = rhou/rho
!             
!             locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
!             neighbouroctal => thisoctal
!             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
!             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*direction, q, rho, rhoe, &
!                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, q11, q22, q33)
!             thisoctal%u_i_minus_1(subcell) = rhou/rho
!
!          endif
!       endif
!    enddo
!  end subroutine setupupm


!Calculate cell pressures, including the optional effects of artificial viscosity
  recursive subroutine computepressureGeneral(grid, thisOctal, withViscosity)
     use inputs_mod, only : etaViscosity, useViscosity, useTensorViscosity
     use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child
    integer :: subcell, i
    real(double) :: biggamma,eta
    logical :: withViscosity
!    logical :: useviscosity                                                                                                                                                                                             


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
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          if (.not.thisOctal%ghostCell(subcell)) then
             !Get pressure, independent of viscosity
             thisoctal%pressure_i(subcell) = getpressure(thisoctal, subcell)
             
             !          if(0 == 1) then
             if (.not.useTensorViscosity) then
                !Calculate viscosity contribution to cell pressure
                biggamma = 0.d0
                if(withViscosity) then
                   if (.not.thisoctal%edgecell(subcell)) then
                      useviscosity = .false.
                      if (thisoctal%u_i_plus_1(subcell) .le. thisoctal%u_i_minus_1(subcell)) useviscosity = .true.
                      !                   if (thisOctal%divV(subcell) < 0.d0) useViscosity = .true.
                      if (useviscosity) then
                         
                         biggamma = 0.25d0 * eta**2 * (thisoctal%u_i_plus_1(subcell) - thisoctal%u_i_minus_1(subcell))**2 &
                              * thisoctal%rho(subcell)
                         
                         !                      bigGamma = eta**2 * (thisOctal%subcellSize*gridDistancescale)**2 * &
                         !                           thisOctal%rho(subcell) * thisOctal%divV(subcell)**2
                      else
                         biggamma = 0.d0
                      endif
                   endif
                   !                if (bigGamma /= 0.d0) then
                   !                   write(*,*) thisOctal%pressure_i(subcell),bigGamma,biggamma2
                   !                endif
 !               else
 !                  biggamma = 0.d0
 !               end if

                   thisoctal%pressure_i(subcell) = thisoctal%pressure_i(subcell) + biggamma
                end if
             endif
          endif
          
          if (isnan(thisoctal%pressure_i(subcell))) then
             write(*,*) "pressureu has nan"
             write(*,*) "velocity: ",thisoctal%rhou(subcell),thisoctal%rhov(subcell), thisoctal%rhow(subcell)
             write(*,*) "rho: ", thisoctal%rho(subcell)
             write(*,*) "temperature ",thisOctal%temperature(subcell)
             write(*,*) "ios ",thisOctal%iEquationOfState(subcell)
             write(*,*) "cen ",subcellcentre(thisoctal, subcell)
             write(*,*) "rank ",myrankGlobal
             write(*,*) "mpithread ",thisOctal%mpiThread(subcell)
             write(*,*) "depth ",thisOctal%nDepth
             write(*,*) "ghost ", thisoctal%ghostcell(subcell)
             thisOctal%rho(subcell) = sqrt(-thisOctal%rho(subcell))
          endif
       endif
    enddo
  end subroutine computepressureGeneral

!Calculate the modification to cell velocity and energy due to the pressure gradient
  recursive subroutine pressureforce(thisoctal, dt, grid, direction)
    use mpi
    use inputs_mod, only  : radiationpressure
    type(octal), pointer   :: thisoctal
    type(gridtype) :: grid
    type(VECTOR) :: direction
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, rhou, dx, dv
    real(double) :: eps
    real(double) :: x_i_plus_half, x_i_minus_half, p_i_plus_half, p_i_minus_half, u_i_minus_half, u_i_plus_half
    real(double) :: phi_i_plus_half, phi_i_minus_half, fac1, fac2
    real(double) :: iniRhoe

    eps = smallestCellSize * 1.d10
    

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call pressureforce(child, dt, grid, direction)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then

             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrankGlobal," error in setting up x_i values in pressure force"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell), thisoctal%x_i(subcell)
                write(*,*) thisoctal%ndepth
                write(*,*) "centre ",subcellcentre(thisoctal,subcell)
                write(*,*) "thread ",thisoctal%mpithread(subcell)
             endif

             if (direction%x > 0.d0) then
                rhou = thisoctal%rhou(subcell)
             else if (direction%y > 0.d0) then
                rhou = thisOctal%rhov(subcell)
             else
                rhou = thisOctal%rhow(subcell)
             endif

             dx = thisoctal%subcellsize * griddistancescale

             dv = cellVolume(thisOctal, subcell) * 1.d30
  
!modify the cell velocity due to the pressure gradient
             iniRhoe = thisOctal%rhoe(subcell)

             x_i_plus_half = thisOctal%x_i(subcell) + thisOctal%subcellSize*gridDistanceScale/2.d0
             x_i_minus_half = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistanceScale/2.d0

             fac1 = (x_i_minus_half - thisOctal%x_i_minus_1(subcell)) / &
                  (thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell))
             fac2 = (x_i_plus_half - thisOctal%x_i(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i(subcell))

             p_i_minus_half = thisOctal%pressure_i_minus_1(subcell) + (thisOctal%pressure_i(subcell) - &
                  thisOctal%pressure_i_minus_1(subcell)) * fac1
             p_i_plus_half = thisOctal%pressure_i(subcell) + (thisOctal%pressure_i_plus_1(subcell) - &
                  thisOctal%pressure_i(subcell)) * fac2
             u_i_minus_half = thisOctal%u_i_minus_1(subcell) + (thisOctal%u_i(subcell) - &
                  thisOctal%u_i_minus_1(subcell)) * fac1
             u_i_plus_half = thisOctal%u_i(subcell) + (thisOctal%u_i_plus_1(subcell) - &
                  thisOctal%u_i(subcell)) * fac2 
             phi_i_minus_half = thisOctal%phi_i_minus_1(subcell) + (thisOctal%phi_i(subcell) - &
                  thisOctal%phi_i_minus_1(subcell)) * fac1
             phi_i_plus_half = thisOctal%phi_i(subcell) + (thisOctal%phi_i_plus_1(subcell) - &
                  thisOctal%phi_i(subcell)) * fac2

             dx = x_i_plus_half - x_i_minus_half

             if (direction%x > 0.d0) then
                thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - dt * &
                     (p_i_plus_half - p_i_minus_half) / dx
                if (useTensorViscosity) then
!                   thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - divQ(thisOctal, subcell, 1, grid)*dt
                endif
                thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - dt * & !gravity due to gas
                     thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx
                if (radiationPressure) then
                   thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + &
                        dt * thisOctal%kappaTimesFlux(subcell)%x/cspeed
                endif

             else if (direction%y > 0.d0) then
                thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - dt * &
                     (p_i_plus_half - p_i_minus_half) / dx
                if (useTensorViscosity) then
!                   thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - divQ(thisOctal, subcell, 1, grid)*dt
                endif
                thisoctal%rhov(subcell) = thisoctal%rhov(subcell) - dt * & !gravity due to gas
                     thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx
                if (radiationPressure) then
                   thisOctal%rhov(subcell) = thisOctal%rhov(subcell) + &
                        dt * thisOctal%kappaTimesFlux(subcell)%y/cspeed
                endif
             else
                thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - dt * &
                     (p_i_plus_half - p_i_minus_half) / dx
                if (useTensorViscosity) then
!                   thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - divQ(thisOctal, subcell, 1, grid)*dt
                endif
                thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - dt * & !gravity due to gas
                     thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx
                if (radiationPressure) then
                   thisOctal%rhow(subcell) = thisOctal%rhow(subcell) + &
                        dt * thisOctal%kappaTimesFlux(subcell)%z/cspeed
                endif

             endif


!Modify the cell rhoe due to pressure and gravitaitonal potential gradient
             if (thisoctal%iequationofstate(subcell) /= 1) then             

                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * &
                  (p_i_plus_half * u_i_plus_half - p_i_minus_half * u_i_minus_half) / dx


!                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * &
!                     (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
!                      thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / (2.d0*dx)


                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - dt* & !gravity
                  rhou  * (phi_i_plus_half - phi_i_minus_half) / dx


                if (useTensorviscosity.and.(direction%x > 0.d0)) then
                   thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt*thisOctal%chiline(subcell)
                endif


             endif
!
 !            if(grid%geometry == "SB_coolshk") then
!                print *, "ENFORCING COOLING "
  !              call enforceCooling(grid%octreeRoot, dt, grid, iniRhoe)
  !           end if

          endif
       endif
    enddo
  end subroutine pressureforce

!Calculate the modification to cell velocity and energy due to the pressure gradient
  recursive subroutine pressureforceCylindrical(thisoctal, dt, grid, direction)
    use mpi
    use inputs_mod, only : includePressureTerms, radiationpressure
    type(octal), pointer   :: thisoctal
    type(gridtype) :: grid
    type(VECTOR) :: direction, fVisc, rVec, gravForceFromSinks, cellCentre
    type(octal), pointer  :: child 
    integer :: subcell, i
    logical :: debug
    real(double) :: dt, rhou, dx, dv, speed
    real(double) :: eps
    real(double) :: x_i_plus_half, x_i_minus_half, p_i_plus_half, p_i_minus_half, u_i_minus_half, u_i_plus_half, r
    real(double) :: rhorv_i_plus_half, rhorv_i_minus_half
    real(double) :: phi_i_plus_half, phi_i_minus_half, fac1, fac2

    eps = smallestCellSize * 1.d10
    

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call pressureforceCylindrical(child, dt, grid, direction)
               exit
             end if
          end do
       else
          if (.not.associated(thisOctal%fViscosity)) then
             allocate(thisOctal%fViscosity(1:thisOctal%maxChildren))
             thisoctal%fViscosity = VECTOR(0.d0, 0.d0, 0.d0)
          endif

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then
             
             debug = .false.
             speed = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)/1.d5

             if (speed > 50.d0) debug = .true.
!             if (inSubcell(thisOctal, subcell, VECTOR(10.d3, 0.d0, 1.11d6))) then
!                debug = .true.
!             endif

             if (debug.and.(myHydroSetGlobal == 0)) then
                write(*,*) "speed before " ,speed, 1.d-5*thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                     1.d-5*thisOctal%rhow(subcell)/thisOctal%rho(subcell)
                write(*,*) "direction ",direction
                write(*,*) "u_interface ",thisOctal%u_amr_interface(subcell,1:2)/1.d5
                   write(*,*) "flux small i+1 ",thisOctal%flux_amr_i_plus_1(subcell,1:2)
                   write(*,*) "flux small i  ",thisOctal%flux_amr_i(subcell,1:2)
             endif

!             if (inSubcell(thisOctal, subcell, VECTOR(10.d3, 0.d0, 1.15d6))) then
!                if (myHydroSetGlobal == 0) then
!                   write(*,*) "big u_interface ",thisOctal%u_amr_interface(subcell,1:2)/1.d5
!                   write(*,*) "big flux ",thisOctal%flux_amr_i(subcell,1:2)
!                endif
!             endif
             

             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrankGlobal," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell), thisoctal%x_i(subcell)
                write(*,*) thisoctal%ndepth
                write(*,*) "centre ",subcellcentre(thisoctal,subcell)
                write(*,*) "thread ",thisoctal%mpithread(subcell)
             endif

             if (direction%x > 0.d0) then
                rhou = thisoctal%rhou(subcell)
             else if (direction%y > 0.d0) then
                rhou = thisOctal%rhov(subcell)
             else
                rhou = thisOctal%rhow(subcell)
             endif

             dx = thisoctal%subcellsize * griddistancescale

             dv = cellVolume(thisOctal, subcell) * 1.d30
  
!modify the cell velocity due to the pressure gradient


             rVec = subcellCentre(thisOctal, subcell)
             r = sqrt(rVec%x**2 + rVec%y**2) * gridDistanceScale

             x_i_plus_half = thisOctal%x_i(subcell) + thisOctal%subcellSize*gridDistanceScale/2.d0
             x_i_minus_half = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistanceScale/2.d0

             fac1 = (x_i_minus_half - thisOctal%x_i_minus_1(subcell)) / &
                  (thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell))
             fac2 = (x_i_plus_half - thisOctal%x_i(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i(subcell))

             p_i_minus_half = thisOctal%pressure_i_minus_1(subcell) + (thisOctal%pressure_i(subcell) - &
                  thisOctal%pressure_i_minus_1(subcell)) * fac1
             p_i_plus_half = thisOctal%pressure_i(subcell) + (thisOctal%pressure_i_plus_1(subcell) - &
                  thisOctal%pressure_i(subcell)) * fac2
             u_i_minus_half = thisOctal%u_i_minus_1(subcell) + (thisOctal%u_i(subcell) - &
                  thisOctal%u_i_minus_1(subcell)) * fac1
             u_i_plus_half = thisOctal%u_i(subcell) + (thisOctal%u_i_plus_1(subcell) - &
                  thisOctal%u_i(subcell)) * fac2 
             phi_i_minus_half = thisOctal%phi_i_minus_1(subcell) + (thisOctal%phi_i(subcell) - &
                  thisOctal%phi_i_minus_1(subcell)) * fac1
             phi_i_plus_half = thisOctal%phi_i(subcell) + (thisOctal%phi_i_plus_1(subcell) - &
                  thisOctal%phi_i(subcell)) * fac2

             rhorv_i_minus_half = thisOctal%rhorv_i_minus_1(subcell) + (thisOctal%rhov(subcell) - &
                  thisOctal%rhorv_i_minus_1(subcell)) * fac1
             rhorv_i_plus_half = thisOctal%rhov(subcell) + (thisOctal%rhorv_i_plus_1(subcell) - &
                  thisOctal%rhov(subcell)) * fac2

             fVisc =  newdivQ(thisOctal, subcell,  grid)

             call calculateForceFromSinks(thisOctal, subcell, globalsourceArray, globalnSource, &
                  2.d0 * smallestCellSize*gridDistanceScale, gravForceFromSinks)




             cellCentre = subcellCentre(thisOctal,subcell)
!             if (.not.((cellCentre%z-thisOctal%subcellSize/2.d0 <0.d0) &
!                  .and.(cellcentre%z+thisOctal%subcellSize/2.d0 > 0.d0))) then
!                fVisc = VECTOR(0.d0,0.d0,0.d0)
!             endif
             thisOctal%fViscosity(subcell) = fVisc * 1.d20

             cellCentre = subcellCentre(thisoctal, subcell)

!             if (modulus(fVisc) /= 0.d0) write(*,*) "fvisc ",fvisc

             dx = x_i_plus_half - x_i_minus_half

             if (direction%x > 0.d0) then


                if (includePressureTerms) then
                   thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - dt * &
                        (p_i_plus_half - p_i_minus_half) / dx

 !                  if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!!!                      write(*,*) "u speed over 200 after pressure forces"
  !                 endif
                   if (debug) then
                      if (myHydroSetGlobal == 0) then
                         write(*,*) "change in speed from pressure in x ",(- dt * &
                              (p_i_plus_half - p_i_minus_half) / dx)/(thisOctal%rho(subcell)*1.d5)
                         write(*,*) "alternative pressure gradient ",(-dt * &
                              (thisOctal%pressure_i_plus_1(subcell)-thisOctal%pressure_i_minus_1(subcell))/(2.d0*dx)/ &
                              (thisOctal%rho(subcell)*1.d5))
                      endif
                      
                   endif
                endif

! alpha viscosity

!                if (inSubcell(thisOctal, subcell, VECTOR(60.d3, 0.d0, 0.d0))) then
!                   write(*,*) "direction ",direction
!                   write(*,*) "fvisc innermost", fVisc/thisOctal%rho(subcell)
!                endif
!                if (inSubcell(thisOctal, subcell, VECTOR(75.d3, 0.d0, 0.d0))) then
!                   write(*,*) "fvisc middle", fVisc/thisOctal%rho(subcell)
!                endif
!                if (inSubcell(thisOctal, subcell, VECTOR(90.d3, 0.d0, 0.d0))) then
!                   write(*,*) "fvisc rightmost", fVisc/thisOctal%rho(subcell)
!                endif

                thisoctal%rhou(subcell) = thisoctal%rhou(subcell) + dt * fVisc%x

!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after viscous forces"
!                   endif

                if (debug) then
                   if (myHydroSetGlobal == 0) write(*,*) "change in speed from viscosity in x ", &
                        dt * fVisc%x/(thisOctal%rho(subcell)*1.d5)
                endif


                thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - dt * & !gravity due to gas
                     thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx

!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after gas gravity forces"
!                   endif


                thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + dt * gravForceFromSinks%x ! grav due to sinks

!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after sink gravity forces"
!                   endif

                if (debug) then
                   if (myHydroSetGlobal == 0) write(*,*) "change in speed  gravity in x ", (dt*gravforcefromsinks%x + dt * & !gravity due to gas
                     thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx) / (thisOctal%rho(subcell)*1.d5)
                endif


!                write(*,*) "fvisc ",fVisc%x, " phi ",thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx

                ! now centrifugal term


                thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + dt * (thisOctal%rhov(subcell)**2) &
                     / (thisOctal%rho(subcell)*thisOctal%x_i(subcell)**3)!/dx!**2

                if (debug) then
                   if (myHydroSetGlobal == 0) write(*,*) "change in speed from centrifugal term ", &
                        (dt * (thisOctal%rhov(subcell)**2) &
                     / (thisOctal%rho(subcell)*thisOctal%x_i(subcell)**3))/(thisOctal%rho(subcell)*1.d5)
                endif

!                thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + dt * (thisOctal%rhov(subcell)**2) &
 !                    / (thisOctal%rho(subcell)*thisOctal%x_i(subcell)**3)


                !THAW - NULLIFYING TORQUE FOR DEBUGGING PURPOSES
!                if(0 /= 0) then
                   thisoctal%rhov(subcell) = thisoctal%rhov(subcell) + dt * fVisc%y * r ! torque
 !               end if

!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after centrifugal forces"
!                      write(*,*) "rho ",thisOctal%rho(subcell)
!                   endif




!                if (thisOctal%rho(subcell) > 1.d-14) then
!                   write(*,*) "rhov ",thisOctal%rhov(subcell)
!                   write(*,*) "fvisc ",fVisc
!                   write(*,*) "change (%)",100.d0*dt*fVisc%y*r/thisOctal%rhov(subcell)
!                endif


                if (debug) then
                   if (myHydroSetGlobal == 0) write(*,*) "change in mom from centrifugal term in x (change in km/s): ", &
                        + dt * (thisOctal%rhov(subcell)**2) &
                     / (thisOctal%rho(subcell)**2*thisOctal%x_i(subcell)**3)/1.d5
                endif


                if (radiationPressure) then
                   thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + &
                        dt * thisOctal%kappaTimesFlux(subcell)%x/cspeed

                   if (debug) then
                      if (myHydroSetGlobal == 0) write(*,*) "change in speed from rad pressure in x ", &
                           (dt * thisOctal%kappaTimesFlux(subcell)%x/cspeed)/(thisOctal%rho(subcell)*1.d5)

                   endif
!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after rad pressure forces"
!                   endif

                endif


!Modify the cell rhoe due to pressure and gravitaitonal potential gradient
             if (thisoctal%iequationofstate(subcell) /= 1) then             


                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * &
                  (p_i_plus_half * u_i_plus_half * x_i_plus_half - p_i_minus_half * u_i_minus_half * x_i_minus_half) / &
                  (0.5d0 * (x_i_plus_half+x_i_minus_half) * (x_i_plus_half - x_i_minus_half))


                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - dt* & !gravity
                  rhou  * (phi_i_plus_half * x_i_plus_half - phi_i_minus_half * x_i_minus_half) / &
                  (0.5d0 * (x_i_plus_half+x_i_minus_half) * (x_i_plus_half - x_i_minus_half))

                if (useTensorviscosity.and.(direction%x > 0.d0)) then
                   thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt*thisOctal%chiline(subcell)
                endif

             endif
             speed = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)/1.d5

             if ((debug).and.(myHydroSetGlobal == 0)) then
                write(*,*) "speed after " ,speed, 1.d-5*thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                     1.d-5*thisOctal%rhow(subcell)/thisOctal%rho(subcell)
             endif



             else

                if (includePressureTerms) then
                   thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - dt * &
                        (p_i_plus_half - p_i_minus_half) / dx

!                   if (abs(thisOctal%rhow(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "w speed over 200 after pressure forces"
!                   endif


                endif
                if (debug) then
                   if (myHydroSetGlobal == 0) write(*,*) "change in speed from pressure in z ",(- dt * &
                     (p_i_plus_half - p_i_minus_half) / dx)/(thisOctal%rho(subcell)*1.d5)
                endif

! alpha viscosity
                thisoctal%rhow(subcell) = thisoctal%rhow(subcell) + dt * fVisc%z

!                   if (abs(thisOctal%rhow(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "w speed over 200 after viscous forces"
!                   endif


!                write(*,*) "ratio vr/vtheta: ",thisOctal%rhou(subcell)/(thisOctal%rhov(subcell)/r), thisOctal%rho(subcell),r/1e14


                if (debug) then
                   if (myHydroSetGlobal == 0) write(*,*) "change in speed from viscosity in z ", &
                        dt * fVisc%z/(thisOctal%rho(subcell)*1.d5)
                endif


                thisoctal%rhow(subcell) = thisoctal%rhow(subcell) - dt * & !gravity due to gas
                     thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx

!                   if (abs(thisOctal%rhow(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "w speed over 200 after gas gravity forces"
!                   endif


                thisOctal%rhow(subcell) = thisOctal%rhow(subcell) + dt * gravForceFromSinks%z ! grav due to sinks


                if (debug) then
                   if (myHydroSetGlobal == 0) write(*,*) "change in speed from gravity in z ", (dt * gravforcefromsinks%z + dt * & !gravity due to gas
                     thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx)/(thisoctal%rho(subcell)*1.d5)
                endif

!                   if (abs(thisOctal%rhow(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "w speed over 200 after sink gravity forces"
!                   endif

                if (radiationPressure) then
                   thisOctal%rhow(subcell) = thisOctal%rhow(subcell) + &
                        dt * thisOctal%kappaTimesFlux(subcell)%z/cspeed

!                   if (abs(thisOctal%rhow(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "w speed over 200 after rad pressure forces"
!                   endif

                   if (debug) then
                      if (myHydroSetGlobal == 0) write(*,*) "change in speed from rad pressure in x ", &
                           (dt * thisOctal%kappaTimesFlux(subcell)%z/cspeed)/(thisOctal%rho(subcell)*1.d5)
                   endif

                endif

!Modify the cell rhoe due to pressure and gravitaitonal potential gradient
             if (thisoctal%iequationofstate(subcell) /= 1) then             

                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * &
                  (p_i_plus_half * u_i_plus_half - p_i_minus_half * u_i_minus_half) / dx


!                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * &
!                     (thisOctal%pressure_i_plus_1(subcell) * thisOctal%u_i_plus_1(subcell) - &
!                      thisOctal%pressure_i_minus_1(subcell) * thisOctal%u_i_minus_1(subcell)) / (2.d0*dx)


                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - dt* & !gravity
                  rhou  * (phi_i_plus_half - phi_i_minus_half) / dx


                if (useTensorviscosity.and.(direction%x > 0.d0)) then
                   thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt*thisOctal%chiline(subcell)
                endif

             endif
             speed = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)/1.d5
             if ((debug).and.(myHydroSetGlobal == 0)) then
                write(*,*) "speed after " ,speed, 1.d-5*thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                     1.d-5*thisOctal%rhow(subcell)/thisOctal%rho(subcell)
             endif
          endif
       endif


    endif
 enddo
  end subroutine pressureforceCylindrical


!Calculate the modification to cell velocity and energy due to the pressure gradient
  recursive subroutine pressureforceSpherical(thisoctal, dt, grid, direction)
    use mpi
    use inputs_mod, only : includePressureTerms, radiationpressure
    type(octal), pointer   :: thisoctal
    type(gridtype) :: grid
    type(VECTOR) :: direction, fVisc, rVec, gravForceFromSinks
    type(octal), pointer  :: child 
    integer :: subcell, i
    logical :: debug
    real(double) :: dt, rhou, dx, dv, speed
    real(double) :: eps
    real(double) :: x_i_plus_half, x_i_minus_half, p_i_plus_half, p_i_minus_half, u_i_minus_half, u_i_plus_half, r
    real(double) :: rhorv_i_plus_half, rhorv_i_minus_half
    real(double) :: phi_i_plus_half, phi_i_minus_half, fac1, fac2

    eps = smallestCellSize * 1.d10
    

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call pressureforceSpherical(child, dt, grid, direction)
               exit
             end if
          end do
       else
          if (.not.associated(thisOctal%fViscosity)) then
             allocate(thisOctal%fViscosity(1:thisOctal%maxChildren))
             thisoctal%fViscosity = VECTOR(0.d0, 0.d0, 0.d0)
          endif

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then

          debug = .false.
          speed = sqrt(thisOctal%rhou(subcell)**2)/thisOctal%rho(subcell)/1.d5
          if ((speed > 200.d0).and.(myHydroSetGlobal == 0)) then
             write(*,*) "speed before" ,speed
             debug = .true.
          endif


             if (thisoctal%x_i_plus_1(subcell) == thisoctal%x_i_minus_1(subcell)) then
                write(*,*) myrankGlobal," error in setting up x_i values"
                write(*,*) thisoctal%x_i_plus_1(subcell),thisoctal%x_i_minus_1(subcell), thisoctal%x_i(subcell)
                write(*,*) thisoctal%ndepth
                write(*,*) "centre ",subcellcentre(thisoctal,subcell)
                write(*,*) "thread ",thisoctal%mpithread(subcell)
             endif

             rhou = thisoctal%rhou(subcell)
             dx = thisoctal%subcellsize * griddistancescale

             dv = cellVolume(thisOctal, subcell) * 1.d30
  
!modify the cell velocity due to the pressure gradient


             rVec = subcellCentre(thisOctal, subcell)
             r = sqrt(rVec%x**2) * gridDistanceScale

             x_i_plus_half = thisOctal%x_i(subcell) + thisOctal%subcellSize*gridDistanceScale/2.d0
             x_i_minus_half = thisOctal%x_i(subcell) - thisOctal%subcellSize*gridDistanceScale/2.d0

             fac1 = (x_i_minus_half - thisOctal%x_i_minus_1(subcell)) / &
                  (thisOctal%x_i(subcell) - thisOctal%x_i_minus_1(subcell))
             fac2 = (x_i_plus_half - thisOctal%x_i(subcell)) / &
                  (thisOctal%x_i_plus_1(subcell) - thisOctal%x_i(subcell))

             p_i_minus_half = thisOctal%pressure_i_minus_1(subcell) + (thisOctal%pressure_i(subcell) - &
                  thisOctal%pressure_i_minus_1(subcell)) * fac1
             p_i_plus_half = thisOctal%pressure_i(subcell) + (thisOctal%pressure_i_plus_1(subcell) - &
                  thisOctal%pressure_i(subcell)) * fac2
             u_i_minus_half = thisOctal%u_i_minus_1(subcell) + (thisOctal%u_i(subcell) - &
                  thisOctal%u_i_minus_1(subcell)) * fac1
             u_i_plus_half = thisOctal%u_i(subcell) + (thisOctal%u_i_plus_1(subcell) - &
                  thisOctal%u_i(subcell)) * fac2 
             phi_i_minus_half = thisOctal%phi_i_minus_1(subcell) + (thisOctal%phi_i(subcell) - &
                  thisOctal%phi_i_minus_1(subcell)) * fac1
             phi_i_plus_half = thisOctal%phi_i(subcell) + (thisOctal%phi_i_plus_1(subcell) - &
                  thisOctal%phi_i(subcell)) * fac2

             rhorv_i_minus_half = thisOctal%rhorv_i_minus_1(subcell) + (thisOctal%rhov(subcell) - &
                  thisOctal%rhorv_i_minus_1(subcell)) * fac1
             rhorv_i_plus_half = thisOctal%rhov(subcell) + (thisOctal%rhorv_i_plus_1(subcell) - &
                  thisOctal%rhov(subcell)) * fac2


             call calculateForceFromSinks(thisOctal, subcell, globalsourceArray, globalnSource, &
                  2.d0 * smallestCellSize*gridDistanceScale, gravForceFromSinks)
             thisOctal%fViscosity(subcell) = fVisc * 1.d20
!             if (modulus(fVisc) /= 0.d0) write(*,*) "fvisc ",fvisc

             dx = x_i_plus_half - x_i_minus_half



                if (includePressureTerms) then
                   thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - dt * &
                        (p_i_plus_half - p_i_minus_half) / dx

!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after pressure forces"
!                   endif
                   if (debug) then
                      write(*,*) "change in mom from pressure in x ",- dt * &
                           (p_i_plus_half - p_i_minus_half) / dx
                      write(*,*) "change in km/s ",1.d-5*(-dt * &
                           (p_i_plus_half - p_i_minus_half) / dx)/thisOctal%rho(subcell)
                      write(*,*) "x ",thisOctal%x_i(subcell)
                      write(*,*) "p_i_p_half ",p_i_plus_half
                      write(*,*) "p_i_m_half ",p_i_minus_half
                   endif
                endif

! alpha viscosity

!                thisoctal%rhou(subcell) = thisoctal%rhou(subcell) + dt * fVisc%x

!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after viscous forces"
!                   endif

                if (debug) then
                   write(*,*) "change in mom from viscosity in x ", dt * fVisc%x
                endif


!                thisoctal%rhou(subcell) = thisoctal%rhou(subcell) - dt * & !gravity due to gas
 !                    thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx
!
!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after gas gravity forces"
!                   endif
!
!
               thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + dt * gravForceFromSinks%x ! grav due to sinks

!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after sink gravity forces"
!                   endif

                if (debug) then
                   write(*,*) "change in mom from gravity in x ", dt * & !gravity due to gas
                     thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx
                endif


!                write(*,*) "fvisc ",fVisc%x, " phi ",thisOctal%rho(subcell) * (phi_i_plus_half - phi_i_minus_half) / dx

                ! now centrifugal term


  !              thisoctal%rhov(subcell) = thisoctal%rhov(subcell) + dt * fVisc%y * r ! torque


   !             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + dt * (thisOctal%rhov(subcell)**2) &
   !                  / (thisOctal%rho(subcell)*thisOctal%x_i(subcell)**3)

!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after centrifugal forces"
 !                     write(*,*) "rho ",thisOctal%rho(subcell)
 !                  endif




!                if (thisOctal%rho(subcell) > 1.d-14) then
!                   write(*,*) "rhov ",thisOctal%rhov(subcell)
!                   write(*,*) "fvisc ",fVisc
!                   write(*,*) "change (%)",100.d0*dt*fVisc%y*r/thisOctal%rhov(subcell)
!                endif


                if (debug) then
                   write(*,*) "change in mom from centrifugal term in x (change in km/s): ", + dt * (thisOctal%rhov(subcell)**2) &
                     / (thisOctal%rho(subcell)**2*thisOctal%x_i(subcell)**3)/1.d5
                endif


                if (radiationPressure) then
                   thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + &
                        dt * thisOctal%kappaTimesFlux(subcell)%x/cspeed

                        
!                   if (thisOctal%radiationMomentum(subcell)%x /= 0.d0) write(*,*) "from flux,mom ",thisOctal%kappaTimesFlux(subcell)%x/cspeed,thisOctal%radiationMomentum(subcell)

                   if (debug) then
                      write(*,*) "change in mom from rad pressure in x ", dt * thisOctal%kappaTimesFlux(subcell)%x/cspeed
                      
                   endif
!                   if (abs(thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*1.d5)) > 201.d0) then
!                      write(*,*) "u speed over 200 after rad pressure forces"
 !                  endif

                endif


!Modify the cell rhoe due to pressure and gravitaitonal potential gradient
!             if (thisoctal%iequationofstate(subcell) /= 1) then             !!
!
!
!                thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt * &
!                  (p_i_plus_half * u_i_plus_half * x_i_plus_half - p_i_minus_half * u_i_minus_half * x_i_minus_half) / &
!                  (0.5d0 * (x_i_plus_half+x_i_minus_half) * (x_i_plus_half - x_i_minus_half))!!!!

!
!                thisoctal%rhoe(subcell) = thisoctal%rhoe(subcell) - dt* & !gravity
!                  rhou  * (phi_i_plus_half * x_i_plus_half - phi_i_minus_half * x_i_minus_half) / &
!                  (0.5d0 * (x_i_plus_half+x_i_minus_half) * (x_i_plus_half - x_i_minus_half))!!
!
!                if (useTensorviscosity.and.(direction%x > 0.d0)) then
!                   thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) - dt*thisOctal%chiline(subcell)
!                endif!!!
!
!             endif
             speed = sqrt(thisOctal%rhou(subcell)**2)/thisOctal%rho(subcell)/1.d5
             if ((debug).and.(myHydroSetGlobal == 0)) then
                write(*,*) "speed after " ,speed
             endif

          endif
       endif
    enddo
  end subroutine pressureforceSpherical


!Calculate the modification to cell velocity and energy due to the pressure gradient
  recursive subroutine enforceCooling(thisoctal, dt, grid)!, inirhoe)
    use mpi
    type(octal), pointer   :: thisoctal
    type(gridtype) :: grid
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) ::  dt !, du
    logical, save :: firstTime = .true.!, uKineticSave=.true.
!    real(double), parameter :: Teq = 1.d0
    real(double) ::  u_A!, u_kinetic, u_B
    real, parameter :: temperatureUnit = real(2.33d0*mHydrogen/(kerg))!*(3.d0/2.d0)*((5.d0/3.d0)-1.d0)))
    real :: oldT
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call enforcecooling(child, dt, grid)!, inirhoe)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                

             if(dt /= 0.d0) then
               
                u_A = thisOctal%rhoe(subcell) / thisOctal%rho(subcell)

                oldT = real((thisOctal%gamma(subcell) - 1.d0)*((3.d0/2.d0)))!*temperatureUnit
                
                thisOctal%temperature(subcell) = real((thisOctal%gamma(subcell) - 1.d0)*u_A)

                if(thisOctal%temperature(subcell) < 1.0) then
                   thisOctal%temperature(subcell) = 1.0!temperatureUnit
                end if

                thisOctal%energy(subcell) = u_A - ((dt*256.d0)* &
                     (dble(thisOctal%temperature(subcell)) - dble(oldT))) !- u_kinetic

                if(thisOctal%energy(subcell) < 3.d0/2.d0) then
                   thisOctal%energy(subcell) = 3.d0/2.d0
                end if

                thisOctal%temperature(subcell) = real(thisOctal%temperature(subcell) * temperatureUnit)

                thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
!
                if(thisOCtal%temperature(subcell) /= oldT .and. firstTime) then
                   print *, "COOLING STATUS"
                   print *, "du ", thisOctal%energy(subcell) - u_A
                   print *, "U ", thisOctal%energy(subcell)
                   print *, "U_A", u_A
                   print *, "Teq ", oldT
                   print *, "Temperature ", thisOctal%temperature(subcell)/temperatureUnit
                   firstTime = .false.
                end if
                if(thisOctal%energy(subcell) < 0.0d0) then
                   print *, "negative temperature"
                   print *, "du ", THISOCTal%energy(subcell) - u_A
                   print *, "oldT ", oldT
!                   print *, "Teq ", Teq
                   print *, "Temperature ", thisOctal%temperature(subcell) 
                   print *, (dt*256.d0)*(dble(thisOctal%temperature(subcell)) - dble(oldT))
                endif
             end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
             
          end if
       end if
    end do

  end subroutine enforceCooling


!Use to damp oscillations/flow
  recursive subroutine damp(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
!    real(double) :: speed, fac
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call damp(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
  
!         speed = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
!               thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
!          speed = sqrt(speed)
!          if (speed > 0.d0) then
!             fac = min(speed, 10d6)/speed
!          else
!             fac = 1.d0
!          endif

!          thisOctal%rhou(subcell) = thisOctal%rhou(subcell) * fac
!          thisOctal%rhov(subcell) = thisOctal%rhov(subcell) * fac
!          thisOctal%rhow(subcell) = thisOctal%rhow(subcell) * fac
          if (thisOctal%rho(subcell) < 1.d-25) then
             thisoctal%rhou(subcell) = 1.d-26
             thisoctal%rhov(subcell) = 1.d-26
             thisoctal%rhow(subcell) = 1.d-26
             thisoctal%rhoe(subcell) = 1.d-26
          endif
       endif
    enddo
  end subroutine damp

!Use to damp oscillations/flow
  recursive subroutine computeDivV(thisoctal, grid)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call computeDivV(child, grid)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
  
          if (.not.associated(thisOctal%divV)) then
             allocate(thisOctal%divV(1:thisOctal%maxChildren))
          endif
          thisOctal%divV(subcell) = divV(thisOctal, subcell, grid)
       endif
    enddo
  end subroutine computeDivV

!Use to damp oscillations/flow
  recursive subroutine limitSpeed(thisoctal)
    use inputs_mod, only : hydroSpeedLimit
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    real(double) :: speed,fac
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call limitSpeed(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
  

          if (.not.cylindricalHydro) then
             speed = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + &
                  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
             speed = sqrt(speed)
             if (speed > hydroSpeedLimit) then
                fac = hydroSpeedLimit/speed
             else
                fac = 1.d0
             endif
             
             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) * fac
             thisOctal%rhov(subcell) = thisOctal%rhov(subcell) * fac
             thisOctal%rhow(subcell) = thisOctal%rhow(subcell) * fac
          else
!             print *, "rhou ", thisOctal%rhou(subcell)
!             print *, "rhow ", thisOctal%rhow(subcell)
             speed = (thisOctal%rhou(subcell)**2 + &
                  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
             speed = sqrt(speed)
             if (speed > hydroSpeedLimit) then
                fac = hydroSpeedLimit/speed
             else
                fac = 1.d0
             endif
             
             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) * fac
             thisOctal%rhow(subcell) = thisOctal%rhow(subcell) * fac
          endif

       endif
    enddo
  end subroutine limitSpeed

!Use to damp oscillations/flow
  recursive subroutine cutVacuum(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call cutVacuum(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
  
          if (thisOctal%rho(subcell) < 1.d-27) then
             thisoctal%rhou(subcell) = 1.d-30
             thisoctal%rhov(subcell) = 1.d-30
             thisoctal%rhow(subcell) = 1.d-30
             thisoctal%rhoe(subcell) = 1.d-30
!             thisOctal%ghostCell(subcell) = .true.
!             thisOctal%boundaryCondition(subcell) = 7
          endif
       endif
    enddo
  end subroutine cutVacuum

  recursive subroutine allocateHydrodynamicsAttributes(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call allocateHydrodynamicsAttributes(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle


       call allocateAttribute(thisOctal%q_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_minus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_minus_2,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%u_interface,thisOctal%maxchildren)


       if (thisOctal%twoD) then
          call allocateAttribute(thisOctal%q_amr_i_minus_1,thisOctal%maxchildren,2)
          call allocateAttribute(thisOctal%q_amr_i_plus_1,thisOctal%maxchildren,2)
          call allocateAttribute(thisOctal%u_amr_interface,thisOctal%maxchildren,2)
          call allocateAttribute(thisOctal%u_amr_interface_i_plus_1,thisOctal%maxchildren,2)
          call allocateAttribute(thisOctal%flux_amr_i,thisOctal%maxchildren,2)
          call allocateAttribute(thisOctal%flux_amr_i_plus_1,thisOctal%maxchildren,2)
          call allocateAttribute(thisOctal%flux_amr_i_minus_1,thisOctal%maxchildren,2)
          call allocateAttribute(thisOctal%phiLimit_amr,thisOctal%maxchildren, 2)
       else
          call allocateAttribute(thisOctal%q_amr_i_minus_1,thisOctal%maxchildren,4)
          call allocateAttribute(thisOctal%q_amr_i_plus_1,thisOctal%maxchildren,4)
          call allocateAttribute(thisOctal%u_amr_interface,thisOctal%maxchildren,4)
          call allocateAttribute(thisOctal%u_amr_interface_i_plus_1,thisOctal%maxchildren,4)
          call allocateAttribute(thisOctal%flux_amr_i,thisOctal%maxchildren,4)
          call allocateAttribute(thisOctal%flux_amr_i_plus_1,thisOctal%maxchildren,4)
          call allocateAttribute(thisOctal%flux_amr_i_minus_1,thisOctal%maxchildren,4)
          call allocateAttribute(thisOctal%phiLimit_amr,thisOctal%maxchildren, 4)
       endif

       

       call allocateAttribute(thisOctal%fviscosity,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%x_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%x_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%x_i_minus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%x_i_minus_2,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%u_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%u_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%u_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%flux_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%flux_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%flux_i_minus_1,thisOctal%maxchildren)


       call allocateAttribute(thisOctal%phiLimit,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%qViscosity,thisOctal%maxchildren,3,3)

       call allocateAttribute(thisOctal%ghostCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%corner,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%feederCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%edgeCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%refinedLastTime,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%pressure_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%pressure_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%pressure_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rho_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%divV,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rhou,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhov,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhow,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rhoe,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhoeLastTime,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%energy,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%phi_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_gas,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_stars,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rho_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rho_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rhorv_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhorv_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%boundaryCondition,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%boundaryCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%boundaryPartner,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%gravboundaryPartner,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%rLimit,thisOctal%maxChildren)

       call allocateAttribute(thisOctal%iEquationOfState,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%gamma,thisOctal%maxChildren)


       call allocateAttribute(thisOctal%radiationMomentum,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%kappaTimesFlux,thisOctal%maxChildren)



  
       endif
    enddo
  end subroutine allocateHydrodynamicsAttributes


!copy cell rho to advecting quantity q
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
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
  
          thisoctal%q_i(subcell) = thisoctal%rho(subcell)
        
       endif
    enddo
  end subroutine copyrhotoq

!copy cell rho to advecting quantity q
  recursive subroutine copyTemptoq(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyTemptoq(child)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%temperature(subcell)
        
       endif
    enddo
  end subroutine copyTemptoq

!copy cell rho to advecting quantity q
  recursive subroutine copyrhorvtoq(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyrhorvtoq(child)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%rhov(subcell)
        
       endif
    enddo
  end subroutine copyrhorvtoq

!copy cell ionfrac to advecting quantity q
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

!copy cell ionfrac to advecting quantity q
  recursive subroutine copyDusttoq(thisoctal, dustType)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, dustType
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyDusttoq(child, dustType)
                exit
             end if
          end do
       else
  
          thisoctal%q_i(subcell) = thisoctal%dustTypeFraction(subcell, dustType)
        
       endif
    enddo
  end subroutine copyDusttoq

!copy cell rhoe to advecting quantity q
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

!copy cell rhou to advecting quantity q
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

!copy cell rhov to advecting quantity q
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

!copy cell rhow to advecting quantity q
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
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle

          thisoctal%q_i(subcell) = thisoctal%rhow(subcell)
        
!          if (inSubcell(thisOctal, subcell, VECTOR(50.d3, 0.d0, 1.07d6))) then
!             if (myHydroSetGlobal == 0) write(*,*) "before advect rhow pos ",thisOctal%rhow(subcell)/thisOctal%rho(subcell)/1.d5
!          endif
!          if (inSubcell(thisOctal, subcell, VECTOR(50.d3, 0.d0, -990d3))) then
!             if (myHydroSetGlobal == 0) write(*,*) "before advect rhow neg ",thisOctal%rhow(subcell)/thisOctal%rho(subcell)/1.d5
!          endif
       endif
    enddo
  end subroutine copyrhowtoq

!copy advecting quantity q back to cell rhov
  recursive subroutine copyqtoTemp(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtoTemp(child)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%temperature(subcell) = real(thisoctal%q_i(subcell))
          endif
          
       endif
    enddo
  end subroutine copyqtoTemp


!copy advecting quantity q back to cell rhov
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

!copy advecting quantity q back to cell rhov
  recursive subroutine copyqtorhorv(thisoctal)
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtorhorv(child)
                exit
             end if
          end do
       else
  
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rhov(subcell) = thisoctal%q_i(subcell)
          endif
        
       endif
    enddo
  end subroutine copyqtorhorv

!copy advecting quantity q back to cell rhow
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
            if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle

          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rhow(subcell) = thisoctal%q_i(subcell)
!             if (inSubcell(thisOctal, subcell, VECTOR(50.d3, 0.d0, 1.07d6))) then
!                if (myHydroSetGlobal == 0) write(*,*) "after advect rhow pos ",thisOctal%rhow(subcell)/thisOctal%rho(subcell)/1.d5
!             endif
!             if (inSubcell(thisOctal, subcell, VECTOR(50.d3, 0.d0, -990d3))) then
!                if (myHydroSetGlobal == 0) write(*,*) "after advect rhow neg ",thisOctal%rhow(subcell)/thisOctal%rho(subcell)/1.d5
!             endif

!             if (inSubcell(thisOctal, subcell, VECTOR(140.d3, 0.d0, 625d3))) then
!                if (myHydroSetGlobal == 0) write(*,*) "after advect rhow cen ",thisOctal%rhow(subcell)/thisOctal%rho(subcell)/1.d5
!             endif

          endif
        
       endif
    enddo
  end subroutine copyqtorhow

!copy advecting quantity q back to cell rho
  recursive subroutine copyqtorho(thisoctal, direction)
    use inputs_mod, only : rhoFloor
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
             thisoctal%rho(subcell) = max(thisoctal%q_i(subcell),rhofloor)
          endif

          if ((thisoctal%rho(subcell) < 0.d0)) then
             write(*,*) "rho warning ", thisoctal%rho(subcell), subcellcentre(thisoctal, subcell)
             write(*,*) "rhou, rhov, rhow ", thisOctal%rhou(subcell),thisOctal%rhov(subcell), thisOctal%rhow(subcell)
             write(*,*) "label ",thisoctal%label
             write(*,*) "direction ", direction
             write(*,*) "phi limit ", thisoctal%philimit(subcell)
             write(*,*) "flux ", thisoctal%flux_amr_i_plus_1(subcell,1:4), thisoctal%flux_amr_i(subcell,1:4), &
                  thisoctal%flux_amr_i_minus_1(subcell,1:4)
             write(*,*) "u ", thisoctal%u_amr_interface(subcell,1:4)/1.e5
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
             
             thisOctal%rho(subcell) = SUM(thisOctal%rho(1:thisOctal%maxChildren), &
                  MASK = thisOctal%rho(1:thisOctal%maxChildren) > 0.d0) / dble(thisOctal%maxChildren)
             thisOctal%rhou(subcell) = SUM(thisOctal%rhou(1:thisOctal%maxChildren), & 
                  MASK = thisOctal%rho(1:thisOctal%maxChildren) > 0.d0) / dble(thisOctal%maxChildren)
             thisOctal%rhov(subcell) = SUM(thisOctal%rhov(1:thisOctal%maxChildren), & 
                  MASK = thisOctal%rho(1:thisOctal%maxChildren) > 0.d0) / dble(thisOctal%maxChildren)
             thisOctal%rhow(subcell) = SUM(thisOctal%rhow(1:thisOctal%maxChildren), &
                  MASK = thisOctal%rho(1:thisOctal%maxChildren) > 0.d0) / dble(thisOctal%maxChildren)

          endif
       endif
    enddo
  end subroutine copyqtorho

!copy advecting quantity q back to cell rho
  recursive subroutine copyqtorhoCylindrical(thisoctal, direction)
    use inputs_mod, only : rhoFloor
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
                call copyqtorhoCylindrical(child, direction)
                exit
             end if
          end do
       else
  
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle

          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rho(subcell) = max(thisoctal%q_i(subcell),rhoFloor)

          endif
       endif
    enddo
  end subroutine copyqtorhoCylindrical

!copy advecting quantity q back to cell rho
  recursive subroutine copyqtorhoSpherical(thisoctal, direction)
!    use inputs_mod, only : rhoFloor
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
                call copyqtorhoSpherical(child, direction)
                exit
             end if
          end do
       else
  
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle
          
          if (.not.thisoctal%ghostcell(subcell)) then
             thisoctal%rho(subcell) = max(thisoctal%q_i(subcell),rhoFloor)
          endif
       endif
    enddo
  end subroutine copyqtorhoSpherical

!copy advecting quantity q back to cell ionfrac
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

!copy advecting quantity q back to cell ionfrac
  recursive subroutine copyqtoDust(thisoctal, direction, dustType)
    type(vector) :: direction
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, dustType
  
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call copyqtoDust(child, direction, dustType)
                exit
             end if
          end do
       else
  
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle

          thisoctal%dustTypeFraction(subcell,dustType) = thisoctal%q_i(subcell)
          if(thisoctal%dustTypeFraction(subcell,dustType) > 1.d0) then
             thisoctal%dustTypeFraction(subcell,dustType) = 1.d0
          end if
       endif
    enddo
  end subroutine copyqtoDust

!copy advecting quantity q back to cell rhoe
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

!copy advecting quantity q back to cell rhou
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
!             if (isnan(thisoctal%rhou(subcell))) write(*,*) "nan in q to rhou"
          endif
        
       endif
    enddo
  end subroutine copyqtorhou

!copy cell rho to q, advect q, copy q back to cell rho
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

!copy cell rho to q, advect q, copy q back to cell rho
  subroutine advectTemperature(grid, direction, dt, npairs, thread1, &
       thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyTemptoq(grid%octreeroot)
    call advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtoTemp(grid%octreeroot)

  end subroutine advectTemperature

!copy cell rho to q, advect q, copy q back to cell rho
  subroutine advectrhoCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhotoq(grid%octreeroot)
    call advectqCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhoCylindrical(grid%octreeroot, direction)

  end subroutine advectrhoCylindrical

!copy cell rho to q, advect q, copy q back to cell rho
  subroutine advectrhoCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhotoq(grid%octreeroot)
    call advectqCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhoCylindrical(grid%octreeroot, direction)

  end subroutine advectrhoCylindrical_amr

!copy cell rho to q, advect q, copy q back to cell rho
  subroutine advectrho_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhotoq(grid%octreeroot)
    call advectq_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorho(grid%octreeroot, direction)

  end subroutine advectrho_amr

!copy cell rho to q, advect q, copy q back to cell rho
  subroutine advectrhoSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhotoq(grid%octreeroot)
    call advectqSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhoSpherical(grid%octreeroot, direction)

  end subroutine advectrhoSpherical

!copy cell rho to q, advect q, copy q back to cell rho
  subroutine advectrhorvCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhorvtoq(grid%octreeroot)
    call advectqCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhorv(grid%octreeroot)

  end subroutine advectrhorvCylindrical

!copy cell rho to q, advect q, copy q back to cell rho
  subroutine advectrhorvCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhorvtoq(grid%octreeroot)
    call advectqCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhorv(grid%octreeroot)

  end subroutine advectrhorvCylindrical_amr


!copy cell rho to q, advect q, copy q back to cell rho
  subroutine advectrhorvSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhorvtoq(grid%octreeroot)
    call advectqSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhorv(grid%octreeroot)

  end subroutine advectrhorvSpherical

!copy cell ionfrac to q, advect q, copy q back to cell ionfrac
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

!copy cell ionfrac to q, advect q, copy q back to cell ionfrac
  subroutine advectDust(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    use inputs_mod, only : nDustType
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    integer :: i

    do i = 1, nDustType
       call copyDusttoq(grid%octreeroot, i)
       call advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
       call copyqtoDust(grid%octreeroot, direction, i)
    enddo

  end subroutine advectDust


!copy cell rhoe to q, advect q, copy q back to cell rhoe
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

!copy cell rhoe to q, advect q, copy q back to cell rhoe
  subroutine advectrhoe_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoetoq(grid%octreeroot)
    call advectq_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhoe(grid%octreeroot)

  end subroutine advectrhoe_amr


!copy cell rhoe to q, advect q, copy q back to cell rhoe
  subroutine advectrhoeCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoetoq(grid%octreeroot)
    call advectqCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhoe(grid%octreeroot)

  end subroutine advectrhoeCylindrical

!copy cell rhoe to q, advect q, copy q back to cell rhoe
  subroutine advectrhoeCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoetoq(grid%octreeroot)
    call advectqCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhoe(grid%octreeroot)

  end subroutine advectrhoeCylindrical_amr

!copy cell rhoe to q, advect q, copy q back to cell rhoe
  subroutine advectrhoeSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoetoq(grid%octreeroot)
    call advectqSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhoe(grid%octreeroot)

  end subroutine advectrhoeSpherical


!copy cell rhou to q, advect q, copy q back to cell rhou
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

!copy cell rhou to q, advect q, copy q back to cell rhou
  subroutine advectrhouCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoutoq(grid%octreeroot)
    call advectqCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound, dowritedebug=.true.)
    call copyqtorhou(grid%octreeroot)
!    call setOnAxisRhoU(grid%octreeRoot)

  end subroutine advectrhouCylindrical

!copy cell rhou to q, advect q, copy q back to cell rhou
  subroutine advectrhouCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoutoq(grid%octreeroot)
    call advectqCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound, &
         dowritedebug=.true.)
    call copyqtorhou(grid%octreeroot)

  end subroutine advectrhouCylindrical_amr

!copy cell rhou to q, advect q, copy q back to cell rhou
  subroutine advectrhou_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoutoq(grid%octreeroot)
    call advectq_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound, &
         dowritedebug=.true.)
    call copyqtorhou(grid%octreeroot)

  end subroutine advectrhou_amr

!copy cell rhou to q, advect q, copy q back to cell rhou
  subroutine advectrhouSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhoutoq(grid%octreeroot)
    call advectqSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound, dowritedebug=.true.)
    call copyqtorhou(grid%octreeroot)

  end subroutine advectrhouSpherical

!copy cell rhov to q, advect q, copy q back to cell rhov
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

!copy cell rhov to q, advect q, copy q back to cell rhov
  subroutine advectrhov_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhovtoq(grid%octreeroot)
    call advectq_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhov(grid%octreeroot)

  end subroutine advectrhov_amr


!copy cell rhow to q, advect q, copy q back to cell rhow
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

!copy cell rhow to q, advect q, copy q back to cell rhow
  subroutine advectrhow_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhowtoq(grid%octreeroot)
    call advectq_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhow(grid%octreeroot)

  end subroutine advectrhow_amr

!copy cell rhow to q, advect q, copy q back to cell rhow
  subroutine advectrhowCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhowtoq(grid%octreeroot)
    call advectqCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhow(grid%octreeroot)

  end subroutine advectrhowCylindrical

!copy cell rhow to q, advect q, copy q back to cell rhow
  subroutine advectrhowCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhowtoq(grid%octreeroot)
    call advectqCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhow(grid%octreeroot)

  end subroutine advectrhowCylindrical_amr

!copy cell rhow to q, advect q, copy q back to cell rhow
  subroutine advectrhowSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound

    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call copyrhowtoq(grid%octreeroot)
    call advectqSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    call copyqtorhow(grid%octreeroot)

  end subroutine advectrhowSpherical

!Perform the advection on q
  subroutine advectq(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction

    call setupx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
!    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call setupqx(grid%octreeroot, grid, direction)
!   call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call fluxlimiter(grid%octreeroot)
    call constructflux(grid%octreeroot, dt)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupflux(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call updatecellq(grid%octreeroot, dt)

  end subroutine advectq

!Perform the advection on q
  subroutine advectqCylindrical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound, &
       dowritedebug)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    logical,optional :: dowriteDebug
    logical :: writeDebug

    writeDebug = .false.
    if (present(dowritedebug)) writeDebug = doWRiteDebug

    call setupx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupqx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call fluxlimiter(grid%octreeroot)
    call constructfluxCylindrical(grid, grid%octreeroot, dt, direction, writeDebug)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupfluxCylindrical(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call updatecellqCylindrical(grid%octreeroot, dt, direction)

  end subroutine advectqCylindrical

!Perform the advection on q
  subroutine advectqCylindrical_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound, &
       dowritedebug)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    logical,optional :: dowriteDebug
    logical :: writeDebug

    writeDebug = .false.
    if (present(dowritedebug)) writeDebug = doWRiteDebug

    call setupx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupqx2(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call fluxlimiter2(grid%octreeroot)
    call constructfluxCylindrical2(grid, grid%octreeroot, dt, direction, writeDebug)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupfluxCylindrical2(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call updatecellqCylindrical2(grid%octreeroot, dt, direction)

  end subroutine advectqCylindrical_amr

!Perform the advection on q
  subroutine advectq_amr(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound, &
       dowritedebug)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    logical,optional :: dowriteDebug
    logical :: writeDebug

    writeDebug = .false.
    if (present(dowritedebug)) writeDebug = doWRiteDebug

    call setupx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupqx_amr(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call fluxlimiter_amr(grid%octreeroot)
    call constructflux_amr(grid, grid%octreeroot, dt, direction, writeDebug)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupflux_amr(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call updatecellq_amr(grid%octreeroot, dt, direction)

  end subroutine advectq_amr

!Perform the advection on q
  subroutine advectqSpherical(grid, direction, dt, npairs, thread1, thread2, nbound, group, ngroup, usethisbound, &
       dowritedebug)
    integer :: npairs, thread1(:), thread2(:), nbound(:)
    integer :: group(:), ngroup
    integer :: usethisbound
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    logical,optional :: dowriteDebug
    logical :: writeDebug

    writeDebug = .false.
    if (present(dowritedebug)) writeDebug = doWRiteDebug

    call setupx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupqx(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call fluxlimiter(grid%octreeroot)
    call constructfluxSpherical(grid, grid%octreeroot, dt, direction, writeDebug)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call setupfluxSpherical(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, usethisbound=usethisbound)
    call updatecellqSpherical(grid%octreeroot, dt, direction)

  end subroutine advectqSpherical

!Sum all fluxes on the grid
recursive subroutine sumFluxes(thisOctal, dt, totalFlux)
  use mpi
  integer :: subcell, i
  type(octal), pointer :: thisOctal
  type(octal), pointer :: child
  real(double) :: totalFlux, dt


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
        if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
        if (.not.thisoctal%edgecell(subcell)) then
        totalFlux = totalFlux + thisOctal%flux_i(subcell)
        end if
     end if
  end do

end subroutine sumFluxes

!Perform a single hydrodynamics step, in the x direction, for the 1D case. 
  subroutine  hydrostep1d(grid, dt, npairs, thread1, thread2, nbound, group, ngroup)
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    integer :: npairs, thread1(:), thread2(:), nbound(:), group(:), ngroup
    
    direction = vector(1.d0, 0.d0, 0.d0)

!Boundary conditions
    call imposeboundary(grid%octreeroot, grid)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)

    direction = vector(1.d0, 0.d0, 0.d0)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call computeDivV(grid%octreeRoot, grid)
 
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)

!Set up the grid values
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call computepressureGeneral(grid, grid%octreeroot, .false.) !Thaw Rhie-Chow might just need setuppressureu
    call setupupm(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setupui(grid%octreeroot, grid, direction, dt)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)

!Make any necessary rhie-chow interpolation modifications to the advecting velocity
!    if (rhieChow) then
!     call computepressureGeneral(grid, grid%octreeroot, .false.) !Thaw Rhie-Chow might just need setuppressureu
!     call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
!     call setuppressure(grid%octreeroot, grid, direction) !Thaw!
!     call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
!     call rhiechowui(grid%octreeroot, grid, direction, dt)
!     call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
!    end if

!Advect rho, rhou, rhoe
    call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!    call advectTemperature(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    !if running a radiation hydrodynamics calculation, advect the ion fraction
    if(photoionPhysics .and. hydrodynamics) then
!       call advectIonFrac(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!       call advectDust(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    end if
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    

!impose boundary conditions again
!THAW - If the model is isothermal the temperature at the boundaries needs to be updated prior to the 
!new pressure calculations otherwise you use an out of date temperature
    call imposeboundary(grid%octreeroot, grid)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)

    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call computepressureGeneral(grid, grid%octreeroot, .false.) !Thaw Rhie-Chow might just need setuppressureu
    call setupupm(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setupui(grid%octreeroot, grid, direction, dt)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)


!calculate and set up grid pressures
    call computepressureGeneral(grid, grid%octreeroot, .true.) 
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
 !   if(rhieChow) then
 !      call rhiechowui(grid%octreeroot, grid, direction, dt)
 !      call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
 !   end if

!update cell velocities and rhoe's due to pressure gradient
    call pressureforce(grid%octreeroot, dt, grid, direction)
    if(grid%geometry == "SB_coolshk") then
       call enforceCooling(grid%octreeRoot, dt, grid)!, iniRhoe)
    end if
!impose boundary conditions again
    call imposeboundary(grid%octreeroot, grid)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)



  end subroutine hydrostep1d


!Perform a single hydrodynamics step, in x, y and z directions, for the 3D case.     
  subroutine hydrostep3d(grid, timestep, nPairs, thread1, thread2, nBound, &
       group, nGroup,doSelfGrav, perturbPressure)
    use mpi
    use inputs_mod, only : nBodyPhysics, severeDamping, dirichlet, doGasGravity, useTensorViscosity, &
         moveSources, hydroSpeedLimit, nbodyTest
    use starburst_mod, only : updateSourceProperties
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    logical, optional :: doSelfGrav
    logical :: selfGravity
    integer :: group(:), nGroup
    real(double) :: dt, timestep
    type(VECTOR) :: direction
    integer :: iDir, thisBound, i
!    character(len=80) :: cfile
    logical, optional :: perturbPressure

    selfGravity = .true.
    if (PRESENT(doSelfGrav)) selfgravity = doSelfGrav




!boundary conditions
    if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")
    call imposeBoundary(grid%octreeRoot, grid)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)
    if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")

!    if (selfGravity) then
!       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
!       if (dogasgravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!       call zeroSourcepotential(grid%octreeRoot)
!       if (globalnSource > 0) then
!          call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
!       endif
!       call sumGasStarGravity(grid%octreeRoot)
!       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
!   endif

!self gravity (periodic)
    if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")
    if (selfGravity .and. .not. dirichlet) then
       call periodBoundary(grid, justGrav = .true.)
       call transferTempStorage(grid%octreeRoot, justGrav = .true.)
    endif
    if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")

    if (.not.nBodyTest) then
    do iDir = 1, 4
       select case (iDir)
          case(1)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
             if (myrankWorldglobal == 1) call tune(6,"X-direction step")
          case(2)
             direction = VECTOR(0.d0, 1.d0, 0.d0)
             dt = timeStep
             thisBound = 5
             if (myrankWorldglobal == 1) call tune(6,"Y-direction step")
          case(3)
             direction = VECTOR(0.d0, 0.d0, 1.d0)
             dt = timeStep 
             thisBound = 3
             if (myrankWorldglobal == 1) call tune(6,"Z-direction step")
          case(4)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
             if (myrankWorldglobal == 1) call tune(6,"X-direction step")
       end select

       !set up the grid values
       call setupX(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupQX(grid%octreeRoot, grid, direction)
       
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       if(present(perturbPressure)) then
          continue ! do nothing but avoid compiler warning
!          call PerturbPressureGrid(grid%octreeRoot)
       end if 
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       
       call setupUi(grid%octreeRoot, grid, direction, dt)
       call setupUpm(grid%octreeRoot, grid, direction)
       call setupRhoPhi(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)


       !Advect rho, velocities and rhoe
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call advectRho(grid, direction,  dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoV(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)

       !if running a radiation hydrodynamics calculation, advect the ion fraction
       if(photoionPhysics .and. hydrodynamics) then
!          call advectIonFrac(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
!          call advectDust(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       end if
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)

       !boundary conditions
       call imposeboundary(grid%octreeroot, grid)
       call periodboundary(grid)
       call transfertempstorage(grid%octreeroot)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       !set up the grid values
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)


       call setupUi(grid%octreeRoot, grid, direction, dt)
       call setupUpm(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .true.)

       call setupRhoPhi(grid%octreeRoot, grid, direction)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupPressure(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)

       if (useTensorViscosity) then
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call setupViscosity(grid%octreeRoot, grid)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       endif

       !modify rhou and rhoe due to pressure/graviational potential gradient
       call pressureForce(grid%octreeRoot, dt, grid, direction)

       select case(idir)
         case(1)
            if (myrankWorldglobal == 1) call tune(6,"X-direction step")
         case(2)
            if (myrankWorldglobal == 1) call tune(6,"Y-direction step")
         case(3)
            if (myrankWorldglobal == 1) call tune(6,"Z-direction step")
         case(4)
            if (myrankWorldglobal == 1) call tune(6,"X-direction step")
      end select

!      write(cfile,'(a,i1.1,a)') "afterhydro",idir,".vtk"
!      call writeVtkFile(grid, cfile, &
!           valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!           "hydrovelocity","sourceCont   ","pressure     ","radmom       ",     "radforce     ", &
!           "diff         ","dust1        ","u_i          ",  &
!           "phi          ","rhou         ","rhov         ","rhow         ","rhoe         ", &
!           "vphi         ","jnu          ","mu           ", &
!           "fvisc1       ","fvisc2       ","fvisc3       ","crossings    "/))


   enddo
   endif
   if (severeDamping) call damp(grid%octreeRoot)
   if (hydroSpeedLimit /= 0.) call limitSpeed(grid%octreeRoot)
   if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")
   call periodBoundary(grid)
   call imposeBoundary(grid%octreeRoot, grid)
   call transferTempStorage(grid%octreeRoot)
   if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")

   !periodic self gravity boundary conditions
   if (selfGravity .and. .not. dirichlet) then
      call periodBoundary(grid, justGrav = .true.)
      call transferTempStorage(grid%octreeRoot, justGrav = .true.)
   endif
 
   if (myrankWorldglobal == 1) call tune(6,"Accretion onto sources")
   if ((globalnSource > 0).and.(timestep > 0.d0).and.nBodyPhysics) then
      call domyAccretion(grid, globalsourceArray, globalnSource, timestep)
   endif
   if (myrankWorldglobal == 1) call tune(6,"Accretion onto sources")

   if (myrankWorldglobal == 1) call tune(6,"Updating source positions")

   globalSourceArray(1:globalnSource)%age = globalSourceArray(1:globalnSource)%age + timestep*secstoyears
   if (.not.hosokawaTracks) then
      do i = 1, globalnSource
         call updateSourceProperties(globalsourcearray(i))
      enddo
   endif
   if ((globalnSource > 0).and.nBodyPhysics.and.moveSources) then
      call updateSourcePositions(globalsourceArray, globalnSource, timestep, grid)
   else
      globalSourceArray(1:globalnSource)%velocity = VECTOR(0.d0,0.d0,0.d0)
   endif
   if (myrankWorldglobal == 1) call tune(6,"Updating source positions")
   

   if (selfGravity) then
      if (writeoutput) call writeInfo("Solving self gravity...")
      if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
      if (dogasGravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)!, multigrid=.true.)
      call zeroSourcepotential(grid%octreeRoot)
      if (globalnSource > 0) then
         call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
      endif
      call sumGasStarGravity(grid%octreeRoot)
      if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
      if (writeoutput) call writeInfo("Done.")
   endif


 end subroutine hydroStep3d

!Perform a single hydrodynamics step, in x, y and z directions, for the 3D case.     
  subroutine hydrostep3d_amr(grid, timestep, nPairs, thread1, thread2, nBound, &
       group, nGroup,doSelfGrav, perturbPressure)
    use mpi
    use inputs_mod, only : nBodyPhysics, severeDamping, dirichlet, doGasGravity, useTensorViscosity, &
         moveSources, hydroSpeedLimit, nbodyTest
    use inputs_mod, only : pressureSupport
    use starburst_mod, only : updateSourceProperties
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    logical, optional :: doSelfGrav
    logical :: selfGravity
    integer :: group(:), nGroup
    real(double) :: dt, timestep !, totalMass
    type(VECTOR) :: direction
    integer :: iDir, thisBound, i
!    character(len=80) :: plotfile
    logical, optional :: perturbPressure

    selfGravity = .true.
    if (PRESENT(doSelfGrav)) selfgravity = doSelfGrav




!boundary conditions
    if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")
    call imposeBoundary(grid%octreeRoot, grid)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)
    if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")

!    if (selfGravity) then
!       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
!       if (dogasgravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!       call zeroSourcepotential(grid%octreeRoot)
!       if (globalnSource > 0) then
!          call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
!       endif
!       call sumGasStarGravity(grid%octreeRoot)
!       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
!   endif

!self gravity (periodic)
    if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")
    if (selfGravity .and. .not. dirichlet) then
       call periodBoundary(grid, justGrav = .true.)
       call transferTempStorage(grid%octreeRoot, justGrav = .true.)
    endif
    if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")




    if (.not.nBodyTest) then
    do iDir = 1, 4

!   call findMassOverAllThreads(grid, totalMass)
!   if (writeoutput) write(*,*) "Total mass on grid before step ",idir,"  is ", totalmass/mSol

!   write(plotfile,'(a,i1.1,a)') "beforestep",idir,".dat"
!       call writeVtkFile(grid, plotfile, &
!            valueTypeString=(/"rho          ","hydrovelocity"/))

       select case (iDir)
          case(1)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
             if (myrankWorldglobal == 1) call tune(6,"X-direction step")
          case(2)
             direction = VECTOR(0.d0, 1.d0, 0.d0)
             dt = timeStep
             thisBound = 5
             if (myrankWorldglobal == 1) call tune(6,"Y-direction step")
          case(3)
             direction = VECTOR(0.d0, 0.d0, 1.d0)
             dt = timeStep 
             thisBound = 3
             if (myrankWorldglobal == 1) call tune(6,"Z-direction step")
          case(4)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
             if (myrankWorldglobal == 1) call tune(6,"X-direction step")
       end select

       !set up the grid values
       call setupX(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupQX(grid%octreeRoot, grid, direction)
       
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       if(present(perturbPressure)) then
          continue ! do nothing but avoid compiler warning
!          call PerturbPressureGrid(grid%octreeRoot)
       end if 
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       
       call setupUi_amr4(grid%octreeRoot, grid, direction, dt)
       call setupUpm(grid%octreeRoot, grid, direction)
       call setupRhoPhi(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)


       !Advect rho, velocities and rhoe
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call advectRho_amr(grid, direction,  dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoU_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoV_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoW_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoE_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)

       !if running a radiation hydrodynamics calculation, advect the ion fraction
       if(photoionPhysics .and. hydrodynamics) then
!          call advectIonFrac(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
!          call advectDust(grid, direction, dt/2.d0, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       end if
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)

       !boundary conditions
       call imposeboundary(grid%octreeroot, grid)
       call periodboundary(grid)
       call transfertempstorage(grid%octreeroot)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       !set up the grid values
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)


       call setupUi_amr4(grid%octreeRoot, grid, direction, dt)
       call setupUpm(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .true.)

       call setupRhoPhi(grid%octreeRoot, grid, direction)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupPressure(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)

       if (useTensorViscosity) then
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call setupViscosity(grid%octreeRoot, grid)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       endif

       !modify rhou and rhoe due to pressure/graviational potential gradient
       call pressureForce(grid%octreeRoot, dt, grid, direction)
       if (pressuresupport) call calculateTemperatureFromEnergy(grid%octreeRoot)

       select case(idir)
         case(1)
            if (myrankWorldglobal == 1) call tune(6,"X-direction step")
         case(2)
            if (myrankWorldglobal == 1) call tune(6,"Y-direction step")
         case(3)
            if (myrankWorldglobal == 1) call tune(6,"Z-direction step")
         case(4)
            if (myrankWorldglobal == 1) call tune(6,"X-direction step")
      end select

!      write(cfile,'(a,i1.1,a)') "afterhydro",idir,".vtk"
!      call writeVtkFile(grid, cfile, &
!           valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!           "hydrovelocity","sourceCont   ","pressure     ","radmom       ",     "radforce     ", &
!           "diff         ","dust1        ","u_i          ",  &
!           "phi          ","rhou         ","rhov         ","rhow         ","rhoe         ", &
!           "vphi         ","jnu          ","mu           ", &
!           "fvisc1       ","fvisc2       ","fvisc3       ","crossings    "/))


!   call findMassOverAllThreads(grid, totalMass)
!   if (writeoutput) write(*,*) "Total mass on grid afer step ",idir,"  is ", totalmass/mSol

   enddo
   endif
   if (severeDamping) call damp(grid%octreeRoot)
   if (hydroSpeedLimit /= 0.) call limitSpeed(grid%octreeRoot)
   if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")
   call periodBoundary(grid)
   call imposeBoundary(grid%octreeRoot, grid)
   call transferTempStorage(grid%octreeRoot)
   if (myrankWorldglobal == 1) call tune(6,"Boundary conditions")

   !periodic self gravity boundary conditions
   if (selfGravity .and. .not. dirichlet) then
      call periodBoundary(grid, justGrav = .true.)
      call transferTempStorage(grid%octreeRoot, justGrav = .true.)
   endif
 
   if (myrankWorldglobal == 1) call tune(6,"Accretion onto sources")

!   call findMassOverAllThreads(grid, totalMass)
!   if (writeoutput) write(*,*) "Total mass on grid before accretion  is ", totalmass/mSol
   if ((globalnSource > 0).and.(timestep > 0.d0).and.nBodyPhysics) then
      call domyAccretion(grid, globalsourceArray, globalnSource, timestep)
   endif


!   call findMassOverAllThreads(grid, totalMass)
!   if (writeoutput) write(*,*) "Total mass on grid after accretion  is ", totalmass/mSol

   if (myrankWorldglobal == 1) call tune(6,"Accretion onto sources")

   if (myrankWorldglobal == 1) call tune(6,"Updating source positions")

   globalSourceArray(1:globalnSource)%age = globalSourceArray(1:globalnSource)%age + timestep*secstoyears
   if (.not.hosokawaTracks) then
      do i = 1, globalnSource
         call updateSourceProperties(globalsourcearray(i))
      enddo
   endif
   if ((globalnSource > 0).and.nBodyPhysics.and.moveSources) then
      call updateSourcePositions(globalsourceArray, globalnSource, timestep, grid)
   else
      globalSourceArray(1:globalnSource)%velocity = VECTOR(0.d0,0.d0,0.d0)
   endif
   if (myrankWorldglobal == 1) call tune(6,"Updating source positions")
   
!   call findMassOverAllThreads(grid, totalMass)
!   if (writeoutput) write(*,*) "Total mass on grid after accretion  is ", totalmass/mSol


   if (selfGravity) then

      if (writeoutput) call writeInfo("Solving self gravity...")
      if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
      if (dogasGravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)!, multigrid=.true.)
      call zeroSourcepotential(grid%octreeRoot)
      if (globalnSource > 0) then
         call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
      endif
      call sumGasStarGravity(grid%octreeRoot)
      if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
      if (writeoutput) call writeInfo("Done.")
!      call writeVtkFile(grid, "grav.vtk", &
!           valueTypeString=(/"phigas ", "rho    ","chiline","adot   "/))

   endif


 end subroutine hydroStep3d_amr

 recursive subroutine perturbPressureGrid(thisOctal)
   use mpi
   type(octal), pointer :: thisOctal, child
   real(double) :: rand
   integer :: i, subcell
   logical, save :: firsttime=.true.

   if(firsttime) then
      print *, "perturbing pressuregrid"
      firsttime = .false.
   end if

   do subcell = 1, thisOctal%maxChildren
      if (thisOctal%hasChild(subcell)) then
         ! find the child
         do i = 1, thisOctal%nChildren, 1
            if (thisOctal%indexChild(i) == subcell) then
               child => thisOctal%child(i)
               call perturbPressureGrid(child)
               exit
            end if
         end do
      else
         if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
         call randomNumberGenerator(getDouble=rand)
         rand = (2.d0*(rand - 0.5d0))/10.d0
         thisOctal%pressure_i(subcell) = thisOctal%pressure_i(subcell)*(1.d0+rand)
         
      end if
   end do
 end subroutine perturbPressureGrid
   
!Perform a single hydrodynamics step, in x and z directions, for the 2D case.     
  subroutine hydroStep2d(grid, timeStep, nPairs, thread1, thread2, nBound, group, &
       nGroup, perturbPressure)
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup
    real(double) :: timeStep, dt
    type(VECTOR) :: direction
    integer :: idir, thisBound
    logical, optional :: perturbPressure



    !boundary conditions
    call imposeBoundary(grid%octreeRoot, grid)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)
    
    do iDir = 1, 3
!       print *, " DOING DIRECTION ", iDir
       select case (iDir)
          case(1)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
!             dt = timeStep
             thisBound = 2
          case(2)
             direction = VECTOR(0.d0, 0.d0, 1.d0)
             dt = timeStep
             thisBound = 3
          case(3)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
       end select
       
       call setupX(grid%octreeRoot, grid, direction)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupQX(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       !set up grid values
!       call computepressureGeneral(grid, grid%octreeroot, .true.) 
       call computepressureGeneral(grid, grid%octreeroot, .false.) 

       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
    !boundary conditions
    call imposeBoundary(grid%octreeRoot, grid)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)


       if(present(perturbPressure)) then
          continue ! do nothing but avoid compiler warning
!          call PerturbPressureGrid(grid%octreeRoot)
       end if 
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
!       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       call setupUi(grid%octreeRoot, grid, direction, dt)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setupUpm(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)   

       !advect rho, velocities and rhoe
       call advectRho(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoU(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoW(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoE(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       !if running a radiation hydrodynamics calculation, advect the ion fraction
       if(photoionPhysics .and. hydrodynamics) then
!          call advectIonFrac(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
!          call advectDust(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       end if
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)

       !calculate and set up pressures   
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
    !boundary conditions
    call imposeBoundary(grid%octreeRoot, grid)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       call setupUi(grid%octreeRoot, grid, direction, dt)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupUpm(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .true.)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)

    !boundary conditions
    call imposeBoundary(grid%octreeRoot, grid)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupRhoPhi(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupPressure(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       if (useTensorViscosity) then
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call setupViscosity(grid%octreeRoot, grid)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       endif

       !modify rhou and rhoe due to pressure/gravitational potential gradient
       call pressureForce(grid%octreeRoot, dt, grid, direction)
    enddo
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
    call imposeBoundary(grid%octreeRoot, grid)
    call periodBoundary(grid)
    call transferTempStorage(grid%octreeRoot)

 
  end subroutine hydroStep2d

!Perform a single hydrodynamics step, in r and z directions, for the 2D cylindrical case.     
  subroutine hydroStep2dCylindrical(grid, timeStep, nPairs, thread1, thread2, nBound, group, nGroup)
    use inputs_mod, only : doselfGrav, doGasGravity, alphaViscosity
    type(GRIDTYPE) :: grid
    logical :: selfGravity
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup
    real(double) :: timeStep, dt
    type(VECTOR) :: direction
    integer :: idir, thisBound
    selfGravity = doSelfGrav


    !boundary conditions
    call computepressureGeneral(grid, grid%octreeroot, .false.) 
    call imposeBoundary(grid%octreeRoot, grid)
    call transferTempStorage(grid%octreeRoot)


    if (selfGravity) then
       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
       if (dogasgravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call zeroSourcepotential(grid%octreeRoot)
       call sumGasStarGravity(grid%octreeRoot)
       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
    endif

    if (myrankWorldglobal == 1) call tune(6,"Alpha viscosity")
    call setupAlphaViscosity(grid, alphaViscosity, 0.1d0)
    if (myrankWorldglobal == 1) call tune(6,"Alpha viscosity")




    do iDir = 1, 3
       select case (iDir)
          case(1)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
          case(2)
             direction = VECTOR(0.d0, 0.d0, 1.d0)
             dt = timeStep
             thisBound = 3
          case(3)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
       end select
       

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call setupCylindricalViscosity(grid%octreeRoot, grid)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


       call setupX(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupQX(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       !set up grid values
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call setupUi(grid%octreeRoot, grid, direction, dt)
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       call setupUi(grid%octreeRoot, grid, direction, dt)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)   
       call setupUpm(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       !advect rho, velocities and rhoe
       call advectRhoCylindrical(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoUCylindrical(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoWCylindrical(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoRVCylindrical(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoECylindrical(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)


       !calculate and set up pressures   
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call setupUi(grid%octreeRoot, grid, direction, dt)
       call setuppressure(grid%octreeroot, grid, direction)
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       call setupUi(grid%octreeRoot, grid, direction, dt)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupUpm(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .true.)
       call setupRhoPhi(grid%octreeRoot, grid, direction)       
       call setuprhorvplusminus1(grid%octreeRoot, grid)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupPressure(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)

       !modify rhou and rhoe due to pressure/gravitational potential gradient
       if (hydroSpeedLimit /= 0.) call limitSpeed(grid%octreeRoot)
       call pressureForceCylindrical(grid%octreeRoot, dt, grid, direction)
       if (hydroSpeedLimit /= 0.) call limitSpeed(grid%octreeRoot)


       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call imposeBoundary(grid%octreeRoot, grid)
       call transferTempStorage(grid%octreeRoot)

    enddo

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call setupUi(grid%octreeRoot, grid, direction, dt)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    


   if ((globalnSource > 0).and.(dt > 0.d0).and.nBodyPhysics) then
      call domyAccretion(grid, globalsourceArray, globalnSource, timestep)
   endif

   globalSourceArray(1:globalnSource)%age = globalSourceArray(1:globalnSource)%age + timestep * secstoYears

   if ((globalnSource > 0).and.(dt > 0.d0).and.nBodyPhysics.and.moveSources) then
      if (doselfGrav) then
         if (Writeoutput) write(*,*) "Updating source position"
         call updateSourcePositions(globalsourceArray, globalnSource, timestep, grid)
         if (Writeoutput) write(*,*) "Done."
      else
         if (globalnSource == 1) then
            globalSourceArray(1)%position =  globalSourceArray(1)%position + &
                 (timestep * globalSourceArray(1)%velocity/1.d10)
         endif
      endif
   else
      globalSourceArray(1:globalnSource)%velocity = VECTOR(0.d0,0.d0,0.d0)
   endif
   

 
  end subroutine hydroStep2dCylindrical

!Perform a single hydrodynamics step, in r and z directions, for the 2D cylindrical case.     
  subroutine hydroStep2dCylindrical_amr(grid, timeStep, nPairs, thread1, thread2, nBound, group, nGroup)
    use inputs_mod, only : doselfGrav, doGasGravity, alphaViscosity, advectHydro
    type(GRIDTYPE) :: grid
    logical :: selfGravity
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup
    real(double) :: timeStep, dt
    type(VECTOR) :: direction
    integer :: idir, thisBound
    selfGravity = doSelfGrav


    !boundary conditions
    call computepressureGeneral(grid, grid%octreeroot, .false.) 
    call imposeBoundary(grid%octreeRoot, grid)
    call transferTempStorage(grid%octreeRoot)


    if (selfGravity) then
       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
       if (dogasgravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call zeroSourcepotential(grid%octreeRoot)
       call sumGasStarGravity(grid%octreeRoot)
       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
    endif

    dt = timestep
    if (advectHydro) then

    if (myrankWorldglobal == 1) call tune(6,"Alpha viscosity")
    call setupAlphaViscosity(grid, alphaViscosity, 0.1d0)
    if (myrankWorldglobal == 1) call tune(6,"Alpha viscosity")




    do iDir = 1, 3
       select case (iDir)
          case(1)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
          case(2)
             direction = VECTOR(0.d0, 0.d0, 1.d0)
             dt = timeStep
             thisBound = 3
          case(3)
             direction = VECTOR(1.d0, 0.d0, 0.d0)
             dt = timeStep / 2.d0
             thisBound = 2
       end select
       

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call setupCylindricalViscosity(grid%octreeRoot, grid)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


       call setupX(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupQX(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       !set up grid values
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call setupUi_amr(grid%octreeRoot, grid, direction, dt)
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       call setupUi_amr(grid%octreeRoot, grid, direction, dt)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)   
       call setupUpm(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       !advect rho, velocities and rhoe
       call advectRhoCylindrical_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoUCylindrical_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoWCylindrical_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoRVCylindrical_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call advectRhoECylindrical_amr(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)


       !calculate and set up pressures   
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call setupUi_amr(grid%octreeRoot, grid, direction, dt)
       call setuppressure(grid%octreeroot, grid, direction)
       call setupupm(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)
       call setuppressure(grid%octreeroot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=thisBound)

       call setupUi_amr(grid%octreeRoot, grid, direction, dt)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupUpm(grid%octreeRoot, grid, direction)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call computepressureGeneral(grid, grid%octreeroot, .true.)
       call setupRhoPhi(grid%octreeRoot, grid, direction)       
       call setuprhorvplusminus1(grid%octreeRoot, grid)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=thisBound)
       call setupPressure(grid%octreeRoot, grid, direction)
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)

       !modify rhou and rhoe due to pressure/gravitational potential gradient
       if (hydroSpeedLimit /= 0.) call limitSpeed(grid%octreeRoot)
       call pressureForceCylindrical(grid%octreeRoot, dt, grid, direction)
       if (hydroSpeedLimit /= 0.) call limitSpeed(grid%octreeRoot)


       call computepressureGeneral(grid, grid%octreeroot, .false.) 
       call imposeBoundary(grid%octreeRoot, grid)
       call transferTempStorage(grid%octreeRoot)

    enddo

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
    call setupUi(grid%octreeRoot, grid, direction, dt)
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    
    endif

   if ((globalnSource > 0).and.(timestep > 0.d0).and.nBodyPhysics) then
      call domyAccretion(grid, globalsourceArray, globalnSource, timestep)
   endif

   globalSourceArray(1:globalnSource)%age = globalSourceArray(1:globalnSource)%age + timestep * secstoYears

   if ((globalnSource > 0).and.(timestep > 0.d0).and.nBodyPhysics.and.moveSources) then
      if (doselfGrav) then
         if (writeoutput) write(*,*) "Updating source positions..."
         call updateSourcePositions(globalsourceArray, globalnSource, timestep, grid)
         if (writeoutput) write(*,*) "Done."
      else
         if (globalnSource == 1) then
            globalSourceArray(1)%position =  globalSourceArray(1)%position + &
                 (timestep * globalSourceArray(1)%velocity/1.d10)
         endif
      endif
   endif
   if (.not.moveSources) globalSourceArray(1:globalnSource)%velocity = VECTOR(0.d0,0.d0,0.d0)
   
 
  end subroutine hydroStep2dCylindrical_amr


!Perform a single hydrodynamics step, in the x direction, for the 1D case. 
  subroutine  hydrostep1dSpherical(grid, dt, npairs, thread1, thread2, nbound, group, ngroup)
    use inputs_mod, only : doSelfGrav
    type(gridtype) :: grid
    real(double) :: dt
    type(vector) :: direction
    integer :: npairs, thread1(:), thread2(:), nbound(:), group(:), ngroup


    direction = vector(1.d0, 0.d0, 0.d0)

!Boundary conditions
    call imposeboundary(grid%octreeroot, grid)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)

    if (doSelfGrav) then
       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
       if (dogasgravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call zeroSourcepotential(grid%octreeRoot)
       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
    endif


    direction = vector(1.d0, 0.d0, 0.d0)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
 
    call setupX(grid%octreeRoot, grid, direction)
    call setupQX(grid%octreeRoot, grid, direction)

!Set up the grid values
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call computepressureGeneral(grid, grid%octreeroot, .false.) !Thaw Rhie-Chow might just need setuppressureu
    call setupupm(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setupui(grid%octreeroot, grid, direction, dt)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)


!Advect rho, rhou, rhoe
    call advectRhoSpherical(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoUSpherical(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    call advectRhoESpherical(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!    call advectTemperature(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    !if running a radiation hydrodynamics calculation, advect the ion fraction
    if(photoionPhysics .and. hydrodynamics) then
!       call advectIonFrac(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
!       call advectDust(grid, direction, dt, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    end if
    call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound=2)
    
!impose boundary conditions again
!THAW - If the model is isothermal the temperature at the boundaries needs to be updated prior to the 
!new pressure calculations otherwise you use an out of date temperature
    call imposeboundary(grid%octreeroot, grid)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)

    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call computepressureGeneral(grid, grid%octreeroot, .false.) !Thaw Rhie-Chow might just need setuppressureu
    call setupupm(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setupui(grid%octreeroot, grid, direction, dt)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)


!calculate and set up grid pressures
    call computepressureGeneral(grid, grid%octreeroot, .true.) 
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)

!update cell velocities and rhoe's due to pressure gradient
    call pressureforceSpherical(grid%octreeroot, dt, grid, direction)
    if(grid%geometry == "SB_coolshk") then
       call enforceCooling(grid%octreeRoot, dt, grid)!, iniRhoe)
    end if
!impose boundary conditions again

    if (severeDamping) call damp(grid%octreeRoot)
    
    if (myrankWorldglobal == 1) call tune(6,"Accretion onto sources")
    if ((globalnSource > 0).and.(dt > 0.d0).and.nBodyPhysics) then
       call domyAccretion(grid, globalsourceArray, globalnSource, dt)
    endif
    if (myrankWorldglobal == 1) call tune(6,"Accretion onto sources")
    
    if (myrankWorldglobal == 1) call tune(6,"Updating source positions")
    if ((globalnSource > 0).and.(dt > 0.d0).and.nBodyPhysics.and.moveSources) then
       call updateSourcePositions(globalsourceArray, globalnSource, dt, grid)
    else
       globalSourceArray(1:globalnSource)%velocity = VECTOR(0.d0,0.d0,0.d0)
    endif
    if (myrankWorldglobal == 1) call tune(6,"Updating source positions")
    
    call imposeboundary(grid%octreeroot, grid)
    call periodboundary(grid)
    call transfertempstorage(grid%octreeroot)


  end subroutine hydrostep1dSpherical

  subroutine computeCourantTimeNbody(grid, nSource, source, tc)
    use inputs_mod, only : maxDepthAMR
    type(GRIDTYPE) :: grid
    integer :: nSource, i, j
    type(SOURCETYPE) :: source(:)
    real(double) :: tc, tp, r, minr, acc
    real(double) :: smallest

    smallest = 1.d10*grid%octreeRoot%subcellSize/dble(2**maxdepthamr)
    tp = 1.d30
    do i = 1,nSource
       tp = min(tp, 0.5d0*smallest / max(1.d-20,modulus(source(i)%velocity)))
    enddo
    tc = min(tc, tp)
    if (nSource > 1) then
       tp = 1.d30
       minr = 1.d30
       do i = 1, nSource
          do j = 1, nSource
             if (i /= j) then
                r = 1.d10*modulus(source(i)%position - source(j)%position)
                minr = min(minr, r)
             endif
          enddo
       enddo
       minr = min(smallest, minr)
       do i = 1, nSource
          acc = modulus(source(i)%force)/source(i)%mass
          acc = max(1.d-30, acc)
          tp = min(tp, sqrt(minr/ acc))
       enddo
       tc = min(tc, tp)
    endif
  end subroutine computeCourantTimeNbody

  recursive subroutine computeCourantTimeGasSource(grid, thisOctal, nsource, source, tc)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: subcell, i, isource
    real(double) :: tc
    type(VECTOR) :: rVec, rHat, acc
    real(double) :: dx, rMod, eps
  

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computeCourantTimeGasSource(grid, child, nsource, source, tc)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          if (.not.thisOctal%ghostCell(subcell)) then


             dx= smallestCellSize * gridDistanceScale
             eps = 0.5d0*smallestCellSize * gridDistanceScale
    
             acc = VECTOR(0.d0, 0.d0, 0.d0)
             do isource = 1, globalnSource
                rVec = 1.d10*(subcellCentre(thisOctal, subcell)-source(isource)%position)
                rMod = modulus(rVec)
                rHat = rVec
                call normalize(rHat)
                acc = acc + ((-1.d0)*(bigG * source(isource)%mass*rMod/(rMod**2 + eps**2)**1.5d0))*rHat
             enddo
             tc = min(tc, sqrt(dx/max(modulus(acc),1.d-30)))

          endif
 
       endif
    enddo
  end subroutine computeCourantTimeGasSource

  recursive subroutine courantTimeGasSourceCylindrical(grid, thisOctal, nsource, source, tc)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: subcell, i, isource
    real(double) :: tc
    type(VECTOR) :: acc, force
    real(double) :: eps, dx
  

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call courantTimeGasSourceCylindrical(grid, child, nsource, source, tc)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          if (.not.thisOctal%ghostCell(subcell)) then


             eps = 2.5d0*smallestCellSize * gridDistanceScale
             dx = thisOctal%subcellSize * gridDistanceScale

             acc = VECTOR(0.d0, 0.d0, 0.d0)
             do isource = 1, globalnSource
                call calculateForceFromSinks(thisOctal, subcell, source, nSource, eps, force)
                acc = acc + (1.d0/thisOctal%rho(subcell))*force
             enddo
             tc = min(tc, sqrt(dx/max(modulus(acc),1.d-30)))

          endif
 
       endif
    enddo
  end subroutine courantTimeGasSourceCylindrical





!calculate the courant time - the shortest time step that can be taken with no material being 
!advect more than one grid cell. 

  recursive subroutine computeCourantTime(grid, thisOctal, tc)
    use mpi
!    use inputs_mod, only : hydroSpeedLimit
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: tc, dx, cs, speed
  

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
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          if (.not.thisOctal%ghostCell(subcell)) then

             cs = soundSpeed(thisOctal, subcell)

             dx= smallestCellSize * gridDistanceScale

!Use max velocity not average
             if (.not.cylindricalHydro) then
                speed = thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2
             else
                speed = thisOctal%rhou(subcell)**2 + thisOctal%rhow(subcell)**2
             endif
!             if(speed > hydroSpeedLimit) speed = hydrospeedlimit
             
             if(spherical) then
                speed = thisOctal%rhou(subcell)**2
             end if
             
             speed = sqrt(speed)/thisOctal%rho(subcell)
             tc = min(tc, dx / max(1.d-30,(cs + speed)) )
             
!             print *, "COURANT GAS", tc
!             print *, "speed ", speed
!             print *, "dx", dx
!             print *, "cs ", cs
!             print *, "thisOctal%rhou(subcell)", thisOctal%rhou(subcell)
          endif
 
       endif
    enddo
  end subroutine computeCourantTime

  recursive subroutine computeCourantV(grid, thisOctal, vbulk, vsound)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dx, cs, vBulk, vSound, speed

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computeCourantV(grid, child, vbulk, vsound)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          if (.not.thisOctal%ghostCell(subcell)) then

             cs = soundSpeed(thisOctal, subcell)

             dx= smallestCellSize * gridDistanceScale

!Use max velocity not average
             if (.not.cylindricalHydro) then
                speed = thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2
             else
                speed = thisOctal%rhou(subcell)**2 + thisOctal%rhow(subcell)**2
             endif
             if(spherical) then
                speed = thisoctal%rhou(subcell)**2
             end if
             speed = sqrt(speed)/thisOctal%rho(subcell)
             if (speed > vBulk) vBulk = speed
             if (cs > vSound) vSound = cs

          endif
 
       endif
    enddo
  end subroutine computeCourantV

  subroutine pressureGradientTimeStep(grid, dt, npairs,thread1,thread2,nbound,group,ngroup)
    use inputs_mod, only : amr3d
    integer :: nPairs, thread1(:), thread2(:), group(:), nBound(:), ngroup
    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(VECTOR) :: direction
    dt = 1.d30
    
    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call computepressureGeneral(grid, grid%octreeroot, .true.) 
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call setuprhoPhi(grid%octreeroot, grid, direction)    
    call pressureTimeStep(grid%octreeRoot, dt)

    if (amr3d) then
       direction = VECTOR(0.d0, 1.d0, 0.d0)
       call computepressureGeneral(grid, grid%octreeroot, .true.) 
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=5)
       call setuppressure(grid%octreeroot, grid, direction)
       call setuprhoPhi(grid%octreeroot, grid, direction)
       call pressureTimeStep(grid%octreeRoot, dt)
    endif

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call computepressureGeneral(grid, grid%octreeroot, .true.) 
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)
    call setuppressure(grid%octreeroot, grid, direction)
    call setuprhoPhi(grid%octreeroot, grid, direction)
    call pressureTimeStep(grid%octreeRoot, dt)

  end subroutine pressureGradientTimeStep

  recursive subroutine pressureTimeStep(thisoctal, dt)
    use inputs_mod, only : gridDistanceScale, smallestCellSize, includePressureTerms !, cylindricalHydro
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx, acc, acc_adj
    



    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call pressureTimestep(child, dt)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then

             
             dx = thisoctal%subcellsize * griddistancescale

             if (includePressureTerms) then
                ! calculate acc using pressure from cells either side of cell i
                acc =  (1.d0/thisOctal%rho(subcell)) * &
                     (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i_minus_1(subcell)) / (dx)
                
!                acc = acc + &
!                     (thisOctal%phi_i_plus_1(subcell) - thisOctal%phi_i_minus_1(subcell)) / (dx)
                
                acc = max(1.d-30, abs(acc))

                ! calculate acc using pressure from cell i and cell adjacent to cell i
                ! i, i-1
                acc_adj =  (1.d0/thisOctal%rho(subcell)) * &
                     (thisoctal%pressure_i(subcell) - thisoctal%pressure_i_minus_1(subcell)) / (dx)
                
!                acc_adj = acc_adj + &
!                     (thisOctal%phi_i(subcell) - thisOctal%phi_i_minus_1(subcell)) / (dx)
            
                acc = max(acc, abs(acc_adj))
                
                ! i+1, i
                acc_adj =  (1.d0/thisOctal%rho(subcell)) * &
                     (thisoctal%pressure_i_plus_1(subcell) - thisoctal%pressure_i(subcell)) / (dx)
                
!                acc_adj = acc_adj + &
!                     (thisOctal%phi_i_plus_1(subcell) - thisOctal%phi_i(subcell)) / (dx)
            
                acc = max(acc, abs(acc_adj))

               
                
!                dt = min(dt, 0.5d0*sqrt(dx/acc))
                dt = min(dt, sqrt(2.d0*dx/acc))

                if (dt == 0.d0) then
                   write(*,*) "dt is zero from pressure ",0.5d0*sqrt(smallestCellSize*gridDistanceScale/acc)
                endif
             endif


!             if (cylindricalHydro) then
!                acc = (1.d0/thisOctal%rho(subcell)) * (thisOctal%rhov(subcell)**2) &
!                     / (thisOctal%rho(subcell)*thisOctal%x_i(subcell)**3)
!                dt = min(dt, 0.5d0*sqrt(smallestCellSize*gridDistanceScale/max(acc,1.d-10)))
!             endif
             if (dt == 0.d0) then 
                write(*,*) myrankGlobal, " dt is zero ",thisOctal%rho(subcell),thisOctal%rhov(subcell), &
                     thisoctal%pressure_i_plus_1(subcell), &
                     thisOctal%pressure_i_minus_1(subcell),acc
             endif
          endif
       endif
    enddo
  end subroutine pressureTimeStep

  subroutine computeGravityTimeStep(grid, dt, npairs,thread1,thread2,nbound,group,ngroup)
    use inputs_mod, only : amr3d
    integer :: nPairs, thread1(:), thread2(:), group(:), nBound(:), ngroup
    type(GRIDTYPE) :: grid
    real(double) :: dt
    type(VECTOR) :: direction
    dt = 1.d30
    
    direction = VECTOR(1.d0, 0.d0, 0.d0)
    call computepressureGeneral(grid, grid%octreeroot, .true.) 
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=2)
    call setuppressure(grid%octreeroot, grid, direction)
    call setuprhoPhi(grid%octreeroot, grid, direction)    
    call gravityTimeStep(grid%octreeRoot, dt)

    if (amr3d) then
       direction = VECTOR(0.d0, 1.d0, 0.d0)
       call computepressureGeneral(grid, grid%octreeroot, .true.) 
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=5)
       call setuppressure(grid%octreeroot, grid, direction)
       call setuprhoPhi(grid%octreeroot, grid, direction)
       call gravityTimeStep(grid%octreeRoot, dt)
    endif

    direction = VECTOR(0.d0, 0.d0, 1.d0)
    call computepressureGeneral(grid, grid%octreeroot, .true.) 
    call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup, useThisBound=3)
    call setuppressure(grid%octreeroot, grid, direction)
    call setuprhoPhi(grid%octreeroot, grid, direction)
    call gravityTimeStep(grid%octreeRoot, dt)

  end subroutine computeGravityTimeStep

  recursive subroutine gravityTimeStep(thisoctal, dt)
    use inputs_mod, only : gridDistanceScale
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx, acc, acc_adj
    



    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call gravityTimestep(child, dt)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then

             
             dx = thisoctal%subcellsize * griddistancescale

             ! calculate acc using phi from cells either side of cell i
             acc = (thisOctal%phi_i_plus_1(subcell) - thisOctal%phi_i_minus_1(subcell)) / (dx)
             acc = max(1.d-30, abs(acc))

             ! calculate acc using phi from cell i and cell adjacent to cell i
             ! i, i-1
             acc_adj = (thisOctal%phi_i(subcell) - thisOctal%phi_i_minus_1(subcell)) / (dx)
             acc = max(acc, abs(acc_adj))
             
             ! i+1, i
             acc_adj = (thisOctal%phi_i_plus_1(subcell) - thisOctal%phi_i(subcell)) / (dx)
             acc = max(acc, abs(acc_adj))
             
             dt = min(dt, sqrt(2.d0*dx/acc))

          endif
       endif
    enddo
  end subroutine gravityTimeStep

!sum gas and star contributions to total graviational potential
  recursive subroutine sumGasStarGravity(thisOctal)
    use inputs_mod, only : doGasGravity
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

          if (doGasGravity) then
             thisOctal%phi_i(subcell) = thisOctal%phi_gas(subcell) + thisOctal%phi_stars(subcell)
          else
             thisOctal%phi_i(subcell) = thisOctal%phi_stars(subcell)
          endif
       endif
    enddo
  end subroutine sumGasStarGravity

!calculate the sound speed
  function soundSpeed(thisOctal, subcell) result (cs)
    use mpi
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: cs, rhoPhys
    real(double), parameter :: gamma2 = 1.4d0, rhoCrit = 1.d-14

    if(octalOnThread(thisOctal, subcell,myRankGlobal)) then
    select case(thisOctal%iEquationOfState(subcell))
       case(0) ! adiabatic 
          cs = sqrt(thisOctal%gamma(subcell)*getPressure(thisOctal, subcell)/thisOctal%rho(subcell))
      
       case(1) ! isothermal


!          cs = sqrt((thisOctal%gamma(subcell)*getPressure(thisOctal, subcell))&
!               /thisOctal%rho(subcell))
          cs = sqrt(getPressure(thisOctal, subcell)/thisOctal%rho(subcell))
          


       case(2) ! barotropic
          rhoPhys = returnPhysicalUnitDensity(thisOctal%rho(subcell))
          if (rhoPhys < rhoCrit) then
             cs = sqrt(getPressure(thisOctal, subcell)/thisOctal%rho(subcell))
          else
             cs = sqrt(gamma2 *  getPressure(thisOctal, subcell)/thisOctal%rho(subcell))
          endif
       case(3) ! polytropic
          cs = sqrt(thisOctal%gamma(subcell)*thisOctal%pressure_i(subcell)/thisOctal%rho(subcell))
       case(4) ! federrath
          cs = 0.166d5
       case DEFAULT
          write(*,*) "Unknown equation of state passed to sound speed ", thisOctal%iEquationofState(subcell)
          
    end select
    end if

  end function soundSpeed

!The main routine for 1D hydrodynamics
  subroutine doHydrodynamics1d(grid)
    use inputs_mod, only : tStart, tEnd, tDump, dorefine, amrTolerance
    use mpi
    type(gridtype) :: grid
    real(double) :: dt!,  gamma!, mu
    real(double) :: currentTime
    integer :: i
    type(VECTOR) :: direction
    integer :: nHydroThreads
    real(double) :: tc(512), totalMass, totalEnergy
    integer :: it
    integer :: nPairs, thread1(5120), thread2(5120), group(5120), nBound(5120), ngroup
    integer :: iUnrefine, nUnrefine
    real(double) :: nextDumpTime, temptc(512)
    character(len=80) :: plotfile
    real(double) :: iniM, endM, iniE, endE
    integer :: evenUpArray(nHydroThreadsGlobal)
    integer :: ierr

    direction = VECTOR(1.d0, 0.d0, 0.d0)
!!!!! What was this doing here
!    gamma = 7.d0 / 5.d0
!    mu = 2.d0
    nHydroThreads = nHydroThreadsGlobal

    direction = VECTOR(1.d0, 0.d0, 0.d0)

!    mu = 2.d0



    if (myRankWorldGlobal == 1) write(*,*) "CFL set to ", cflNumber

    if (myrankGlobal /= 0) then
!determine which mpi threads are in contact with one another (i.e. share a domain boundary)
       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       do i = 1, nPairs
          if (myrankWorldglobal==1) &
               write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i), " group ", group(i)
       enddo

       call writeInfo("Setting up even up array", TRIVIAL)
       call setupEvenUpArray(grid, evenUpArray)
       call writeInfo("Done", TRIVIAL)

!do initial exchange across boundaries. The exchange gives subdomain boundary cells information about their foreign neighbours
       call writeInfo("Calling exchange across boundary", TRIVIAL)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call writeInfo("Done", TRIVIAL)


!set up any initially defined velocities/rhoe
       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call calculateRhoU(grid%octreeRoot, direction)
       direction = VECTOR(0.d0, 1.d0, 0.d0)
       call calculateRhoV(grid%octreeRoot, direction)
       direction = VECTOR(0.d0, 0.d0, 1.d0)
       call calculateRhoW(grid%octreeRoot, direction)
       call calculateRhoE(grid%octreeRoot, direction)
       call assignRhoeLast(grid%octreeRoot)

!ensure that all cells are within one level of refinement of one another       
       call evenUpGridMPI(grid,.false.,dorefine, evenUpArray)

!refine the grid if necessary       
       if(dorefine) then
          call setAllUnchanged(grid%octreeRoot)
          call refineGridGeneric(grid, amrTolerance, evenuparray)
       call writeInfo("Evening up grid", TRIVIAL)    
       end if
       call evenUpGridMPI(grid, .false., dorefine, evenuparray)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!set up grid values
       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call setupX(grid%octreeRoot, grid, direction)
       call setupQX(grid%octreeRoot, grid, direction)
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

!get total mass/energy on grid (for diagnostics)
    call findEnergyOverAllThreads(grid, totalenergy)
    call findMassOverAllThreads(grid, totalmass)
    iniM = totalMass
    iniE = totalEnergy

!loop until simulation end time

    write(plotfile,'(a,i4.4,a)') "start.dat"
    
    call  dumpValuesAlongLine(grid, plotfile, &
         VECTOR(0.d0,0.d0,0.0d0), VECTOR(1.d0, 0.d0, 0.0d0), 1000)
    
    do while(currentTime < tend)
               
       tc = 0.d0

!find the largest allowed time step on the grid, such that no cell advects material to cells other than
!their nearest neighbours
       if (myrankGlobal /= 0) then
          tc(myrankGlobal) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRankGlobal))
       endif

       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)

       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads)) * dble(cflNumber)

       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

       if ((currentTime + dt).gt.tEnd) then
          nextDumpTime = tEnd
          dt = nextDumpTime - currentTime
       endif
       

       if (myrankWorldGlobal == 1) call tune(6,"Hydrodynamics step")

       if (myrankGlobal /= 0) then
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          
!perform a single hydrodynamics step 
          if (.not.spherical) then
             call hydroStep1d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
          else
             call hydroStep1dSpherical(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
          endif
         
!mass/energy conservation diagnostics
          call findEnergyOverAllThreads(grid, totalenergy)
          call findMassOverAllThreads(grid, totalmass)
          if(myRankWorldGlobal == 1) then
             endM = totalMass
             endE = totalEnergy
             print *, "dM (%) = ", (1.d0 - (endM/iniM))*100.d0
             print *, "dE (%) = ", (1.d0 - (endE/iniE))*100.d0
             print *, "endM", endM
             print *, "iniM", iniM
             print *, "endE", endE
             print *, "iniE", iniE             
          end if

!unrefine the grid where necessary
          if(dounrefine) then
             iUnrefine = iUnrefine + 1             
             if (iUnrefine == 100) then
                if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
                call unrefineCells(grid%octreeRoot, grid, nUnrefine, amrUnrefinetolerance)
                if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
                iUnrefine = 0
             endif
          end if
 
!ensure all celss are within one level of refinement of one another
          call setAllUnchanged(grid%octreeRoot)
          call evenUpGridMPI(grid, .true., dorefine, evenuparray)
          if (myrankWorldGlobal == 1) call tune(6,"Hydrodynamics step")
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call zeroRefinedLastTime(grid%octreeRoot)

!refine the grid where necessary and ensure that all cells are within one level of refinement of one another
          if(dorefine) then
             call setAllUnchanged(grid%octreeRoot)
             call refineGridGeneric(grid, amrTolerance, evenuparray)
          end if
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)          
          call evenUpGridMPI(grid, .true., dorefine, evenuparray)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)          
       else
       end if

!update the simulation time
       currentTime = currentTime + dt

!dump simulation data if dump time is reached!
       if (currentTime .ge. nextDumpTime) then
!       if(0 == 0) then
          if(grid%geometry == "hydro1d") then
             write(plotfile,'(a,i4.4,a)') "sod_",it,".dat"
          else
             write(plotfile,'(a,i4.4,a)') "torus1Dhydro_",it,".dat"
          end if

          if (writeoutput) then
             write(plotfile,'(a,i4.4,a)') "source",grid%idump,".dat"
             globalSourceArray(:)%time = grid%currentTime
             call writeSourceArray(plotfile)
          endif

         
          if(grid%geometry == "SB_isoshck" .or. grid%geometry == "SB_coolshk") then
             call  dumpValuesAlongLine(grid, plotfile, &
                  VECTOR(1.d0,0.d0,0.0d0), VECTOR(2.d0, 0.d0, 0.0d0), 1024)
          else if(grid%geometry == "SB_CD_1Da" .or. grid%geometry == "SB_CD_1Db") then
             call  dumpValuesAlongLine(grid, plotfile, &
                  VECTOR(0.00390625d0,0.d0,0.0d0), VECTOR(1.5d0, 0.d0, 0.0d0), 1000)             
          else
             call  dumpValuesAlongLine(grid, plotfile, &
                  VECTOR(0.d0,0.d0,0.0d0), VECTOR(1.5d0, 0.d0, 0.0d0), 1000)
          end if
          call assignRhoeLast(grid%octreeRoot)
          nextDumpTime = nextDumpTime + tDump
          it = it + 1
          grid%iDump = it
          grid%currentTime = currentTime
       endif
   
    enddo

!final mass/energy conservation diagnostic
    call findEnergyOverAllThreads(grid, totalenergy)
    call findMassOverAllThreads(grid, totalmass)
    
    endM = totalMass
    endE = totalEnergy
    
    if(myRankWorldGlobal == 1) then
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

!The main routine for 3D hydrodynamics
  subroutine doHydrodynamics3d(grid)
    use vtk_mod, only : writeVtkFilenBody
    use inputs_mod, only : tdump, tend, doRefine, doUnrefine, amrTolerance, dumpRadial
    use inputs_mod, only : addSinkParticles, vtutogrid, doSelfGrav
    use mpi
    type(gridtype) :: grid
    real(double) :: dt, tc(512), temptc(512),  mu
    real(double) :: tSourceSource, tGasSource, tPressureGrad, tCourant
    real(double) :: currentTime, tempDouble !, smallTime
    real(double) :: initialmass
    integer :: i, j, it, iUnrefine, iRefine
    character(len=80) :: plotfile
    real(double) :: nextDumpTime, tff,totalMass !, ang
    type(VECTOR) :: direction, viewVec, angMom
    integer :: nUnrefine
    integer :: thread1(5120), thread2(5120), nBound(5120), nPairs
    integer :: nHydroThreads
    integer :: nGroup, group(5120)
!    logical :: doRefine
    integer :: evenUpArray(nHydroThreadsGlobal)
    logical :: refinedSomeCells
    logical, save  :: firstStep = .true.
    integer :: ierr

!    doSelfGrav = .true.

    if (writeoutput) then
       open(57, file="pos.dat", status="unknown", form="formatted")
       close(57)
    endif

!This should not be hardwired ! 
!    if (grid%geometry == "shakara") doSelfGrav = .false.
    if (grid%geometry == "rtaylor") doSelfGrav = .false.
    if (grid%geometry == "diagSod") doSelfGrav = .false.
    if (grid%geometry == "kelvin" ) doSelfGrav = .false.  

    nHydroThreads = nHydroThreadsGlobal

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    mu = 2.d0
!    if(grid%geometry == "SB_gasmix") then
!       mu = 1.d0
 !   end if

    viewVec = VECTOR(-1.d0,0.d0,0.d0)
    viewVec = rotateY(viewVec, 25.d0*degtorad)

    tff = 1.d0 / sqrt(bigG * (1d0*mSol/((4.d0/3.d0)*pi*7.d15**3)))

    if (tdump == 0.d0) then
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

!determine which mpi threads are in contact with one another (i.e. share a domain boundary)
    if (myrankGlobal /= 0) then
       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       do i = 1, nPairs
          if (myrankglobal==1) &
            write(*,'(a,i5,i4,a,i4,a,i3,a,i3)') "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i), " group ",group(i)
       enddo
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)             

       call writeInfo("Setting up even up array", TRIVIAL)
       call setupEvenUpArray(grid, evenUpArray)

!       if (it /= 1) then
!       endif

    endif

!calculate initial values
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

          call evenUpGridMPI(grid,.false., dorefine, evenUpArray)
          call writeInfo("Done", TRIVIAL)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)             

       end if
!       call writeVTKfile(grid, "start.vtk")
!       call writeVtkFile(grid, "start.vtk", &
!            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
!            ,"mpithread    " /))

       if (myrankGlobal /= 0) then

!refine the grid where necessary
          if(doRefine) then
             call writeInfo("Initial refine", TRIVIAL)    
             if (myrankWorldGlobal == 1) call tune(6, "Initial refine")
             call setAllUnchanged(grid%octreeRoot)
             call refineGridGeneric(grid, amrTolerance, evenuparray)
!             call writeVtkFile(grid, "afterinitrefine.vtk", &
!                  valueTypeString=(/"rho          ","velocity     ","rhoe         " , &
!                  "u_i          ", &
!                  "hydrovelocity", &
!                  "rhou         ", &
!                  "rhov         ", &
!                  "rhow         ", &
!                  "phi          ", &
!                  "pressure     ", &
!                  "q_i          "/))

             call writeInfo("Evening up grid", TRIVIAL)    
          end if
       end if
       
       if(myRankGlobal /= 0) then

!          if(doUnRefine) then
!             do i = 1, 3
!                if (myrank == 1)call tune(6, "Unrefine grid")
!                ! nUnrefine = 0
!                call unrefineCells(grid%octreeRoot, grid, nUnrefine, amrtolerance)
!                call evenUpGridMPI(grid, .true., dorefine, evenUpArray)
!                !          write(*,*) "Unrefined ", nUnrefine, " cells"
!                if (myrank == 1)call tune(6, "Unrefine grid")
!                iUnrefine = 0
!             end do
!          end if


!evening up the grid ensures that no two neighbouring cells differ by more than one level of refinement          
          if(dorefine) then
!             call writeVtkFile(grid, "beforeinitevenup.vtk", &
!                  valueTypeString=(/"rho          ","velocity     ","rhoe         " , &
!                  "u_i          ", &
!                  "hydrovelocity", &
!                  "rhou         ", &
!                  "rhov         ", &
!                  "rhow         ", &
!                  "phi          ", &
!                  "pressure     ", &
!                  "q_i          "/))
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)             
             call evenUpGridMPI(grid,.false., dorefine, evenUpArray)
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)             
             if (myrankWorldGlobal == 1) call tune(6, "Initial refine")
!             call writeVtkFile(grid, "afterinitevenup.vtk", &
!                  valueTypeString=(/"rho          ","velocity     ","rhoe         " , &
!                  "u_i          ", &
!                  "hydrovelocity", &
!                  "rhou         ", &
!                  "rhov         ", &
!                  "rhow         ", &
!                  "phi          ", &
!                  "pressure     ", &
!                  "q_i          "/))
          end if
       end if

!do initial exchange across boundaries. The exchange gives subdomain boundary cells information about their foreign neighbours
       call writeInfo("Calling exchange across boundary", TRIVIAL)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call writeInfo("Done", TRIVIAL)

       
       if(myRankGlobal /= 0) then
          call evenUpGridMPI(grid,.false., dorefine, evenUpArray)
          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call setupX(grid%octreeRoot, grid, direction)
          call setupQX(grid%octreeRoot, grid, direction)
          if (doselfGrav) then
             

             call findMassOverAllThreads(grid, totalmass)
             if (writeoutput) write(*,*) "Total mass: ",totalMass/msol, " solar masses"


!initial gravity 
             if (myrankWorldGlobal == 1) call tune(6, "Self-Gravity")
             if (myrankWorldGlobal == 1) write(*,*) "Doing multigrid self gravity"
!             call writeVtkFile(grid, "beforeselfgrav.vtk", &
!                  valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
!                  ,"phi          " /))
             
             call zeroPhiGas(grid%octreeRoot)
             if (doGasGravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.) 
             
!for periodic self-gravity
             if(.not. dirichlet) then
                call periodBoundary(grid, justGrav = .true.)
                call transferTempStorage(grid%octreeRoot, justGrav = .true.)
                !             if (myrankglobal == 1) call tune(6,"Periodic boundary")
             end if
                         
             call zeroSourcepotential(grid%octreeRoot)
             if (globalnSource > 0) then
                call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
             endif
             call sumGasStarGravity(grid%octreeRoot)

!             call writeVtkFile(grid, "afterselfgrav.vtk", &
!                  valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", &
!                  "phigas       ","phi          "/))
             if (myrankWorldGlobal == 1) write(*,*) "Done"
             if (myrankWorldGlobal == 1) call tune(6, "Self-Gravity")
          endif          
       endif
    endif
    

       if (doSelfGrav .and. grid%geometry == "gravtest") then

          write(plotfile,'(a,i4.4,a)') "radial",1,".dat"
          call  dumpValuesAlongLine(grid, plotfile, VECTOR(0.d0,0.d0,0.0d0), &
               VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0),1000)
          goto 666

       end if




    tc = 0.d0

!calculate largest timestep that each cell on the grid can take without advecting a quantity further than their
!nearest neighbour
    if (myrankGlobal /= 0) then
       tc(myrankGlobal) = 1.d30
       tCourant = 1.d30
       call computeCourantTime(grid, grid%octreeRoot, tCourant)
       tSourceSource = 1.d30
       if (nbodyPhysics) call computeCourantTimeNbody(grid, globalnSource, globalsourceArray, tSourceSource)
       tGasSource = 1.d30
       if (nbodyPhysics) call computeCourantTimeGasSource(grid, grid%octreeRoot, globalnsource, globalsourceArray, tGasSource)
       tPressureGrad = 1.d30
       call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
       call computeDivV(grid%octreeRoot, grid)
       call pressureGradientTimeStep(grid, tPressureGrad, npairs,thread1,thread2,nbound,group,ngroup)
       if (useTensorViscosity) then
          call viscousTimeScale(grid%octreeRoot, grid, dt)
       endif
       tc(myRankGlobal) = min(tCourant, tSourceSource, tGasSource, tPressureGrad)
       call MPI_ALLREDUCE(tCourant, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, AMRCommunicator, ierr)
       tCourant = tempDouble
       call MPI_ALLREDUCE(tSourceSource, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, AMRCommunicator, ierr)
       tSourceSource = tempDouble
       call MPI_ALLREDUCE(tGasSource, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, AMRCommunicator, ierr)
       tGasSource = tempDouble
       call MPI_ALLREDUCE(tPressureGrad, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, AMRCommunicator, ierr)
       tPressureGrad = tempDouble
       if (writeOutput) then
          write(*,'(a,1pe12.3)') "Normal courant time ", tCourant
          write(*,'(a,1pe12.3)') "Source-source courant time ", tSourceSource
          write(*,'(a,1pe12.3)') "Gas-source courant time ", tGasSource
          write(*,'(a,1pe12.3)') "Pressure grad courant time ", tPressureGrad
       endif
    endif



    call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)
    if (firstStep) then
       firstStep = .false.
    endif
    dt = MINVAL(temptc(1:nHydroThreads))

    dt = dt * dble(cflNumber)

    if (writeoutput) write(*,*) "Courant time is ",dt

    if (grid%geometry=="rtaylor") then
       tdump =  20.d0 * dt
    endif

    if (myRankWorldGlobal == 1) write(*,*) "CFL set to ", cflNumber

    if (myrankglobal==1) write(*,*) "Setting tdump to: ", tdump
    if (it ==0) nextDumpTime = 0.
    iUnrefine = 0
    irefine = 0

    call findMassoverAllThreads(grid, initialMass)

       if (nBodyPhysics) then
          initialMass = initialMass + SUM(globalSourceArray(1:globalnSource)%mass)
       endif

       if (writeoutput.and.(globalnSource > 0)) then
          open(45, file="positions.dat", status="unknown",form="formatted")
       endif

!loop until end of simulation time
    do while(currentTime <= tend)
    tc = 0.d0

!calculate largest timestep that each cell on the grid can take without advecting a quantity further than their
!nearest neighbour
    if (myrankGlobal /= 0) then
       tc(myrankGlobal) = 1.d30
       tCourant = 1.d30
       call computeCourantTime(grid, grid%octreeRoot, tCourant)
       tSourceSource = 1.d30
       if (nbodyPhysics) call computeCourantTimeNbody(grid, globalnSource, globalsourceArray, tSourceSource)
       tGasSource = 1.d30
       if (nbodyPhysics) call computeCourantTimeGasSource(grid, grid%octreeRoot, globalnsource, globalsourceArray, tGasSource)
       tPressureGrad = 1.d30
       call pressureGradientTimeStep(grid, tPressureGrad, npairs,thread1,thread2,nbound,group,ngroup)
       if (useTensorViscosity) then
          call viscousTimeScale(grid%octreeRoot, grid, dt)
       endif
       tc(myRankGlobal) = min(tCourant, tSourceSource, tGasSource, tPressureGrad)

       if (nbodyTest) tc(myrankGlobal) = min(tSourceSource, tPressureGrad)

       call MPI_ALLREDUCE(tCourant, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, AMRCommunicator, ierr)
       tCourant = tempDouble
       call MPI_ALLREDUCE(tSourceSource, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, AMRCommunicator, ierr)
       tSourceSource = tempDouble
       call MPI_ALLREDUCE(tGasSource, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, AMRCommunicator, ierr)
       tGasSource = tempDouble
       call MPI_ALLREDUCE(tPressureGrad, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, AMRCommunicator, ierr)
       tPressureGrad = tempDouble
       if (writeOutput) then
          write(*,'(a,1pe12.3)') "Normal courant time ", tCourant
          write(*,'(a,1pe12.3)') "Source-source courant time ", tSourceSource
          write(*,'(a,1pe12.3)') "Gas-source courant time ", tGasSource
          write(*,'(a,1pe12.3)') "Pressure grad courant time ", tPressureGrad
       endif
    endif



    call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)
    if (firstStep) then
       firstStep = .false.
    endif
    dt = MINVAL(temptc(1:nHydroThreads))

    dt = dt * dble(cflNumber)

    if (writeoutput) write(*,*) "Courant time is ",dt
       if (myrankGlobal==1)write(*,*) "Current time is ",returnPhysicalUnitTime(currentTime)
       call findMassoverAllThreads(grid, totalmass)
       if (nBodyPhysics) then
          totalmass = totalMass + SUM(globalSourceArray(1:globalnSource)%mass)
       endif
       if (myrankGlobal==1)write(*,*) "Current mass: ",totalmass/initialmass,initialMass/msol

       call findAngMomoverAllThreads(grid, angMom, VECTOR(0.d0, 0.d0, 0.d0))
       if (nBodyPhysics.and.(globalnSource > 0)) then
          do i = 1, globalnSource
             angMom = angMom + globalSourceArray(i)%angMomentum
          enddo
       endif
       if (myrankGlobal==1)write(*,*) "Current ang mom: ",modulus(angMom), angmom

       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

       if (myrankGlobal /= 0) then
          if (myrankWorldGlobal == 1) write(*,*) "courantTime", dt,cflnumber
          if (myrankWorldGlobal == 1) call tune(6,"Hydrodynamics step")


          
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!perform a hydrodynamics step in the x, y and z directions
          call hydroStep3d_amr(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup, doSelfGrav=doSelfGrav)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


!add/merge sink particles where necessary
          if (nbodyPhysics.and.addSinkParticles) call addSinks(grid, globalsourceArray, globalnSource)


          if (nbodyPhysics) call mergeSinks(grid, globalsourceArray, globalnSource)

          do j = 1, globalnSource
             call emptySurface(globalsourceArray(j)%surface)
             call buildSphereNbody(globalsourceArray(j)%position, globalsourceArray(j)%accretionRadius/1.d10, &
                  globalsourceArray(j)%surface, 20)
          enddo

          if (myrankWorldGlobal == 1) call tune(6,"Hydrodynamics step")

!unrefine the grid where necessary

          iUnrefine = iUnrefine + 1
          if(doUnRefine) then
             if (iUnrefine == 5) then
                if (myrankWorldGlobal == 1)call tune(6, "Unrefine grid")
               ! nUnrefine = 0
                call setAllUnchanged(grid%octreeRoot)
                call unrefineCells(grid%octreeRoot, grid, nUnrefine, amrtolerance)
                call evenUpGridMPI(grid, .true., dorefine, evenUpArray)
                                !          write(*,*) "Unrefined ", nUnrefine, " cells"
                if (myrankWorldGlobal == 1)call tune(6, "Unrefine grid")
                iUnrefine = 0
             endif
          end if

!refine the grid where necessary
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


          if(doRefine) then
!             call writeVtkFile(grid, "beforerefine.vtk", &
!                  valueTypeString=(/"rho          ","velocity     ","rhoe         " , &
!                  "u_i          ", &
!                  "hydrovelocity", &
!                  "rhou         ", &
!                  "rhov         ", &
!                  "rhow         ", &
!                  "phi          ", &
!                  "pressure     ", &
!                  "q_i          "/))

!             iRefine = iRefine + 1 
!             if (irefine == 3) then
             if (myrankWorldGlobal == 1)call tune(6, "Refine grid")
                call writeInfo("Refining grid part 2", TRIVIAL)    
                call setAllUnchanged(grid%octreeRoot)
                call refineGridGeneric(grid, amrTolerance, evenuparray, refinedSomeCells=refinedSomeCells)
                if (refinedSomeCells) then
                   call writeInfo("Done the refine part", TRIVIAL)                 
                   call evenUpGridMPI(grid, .true., dorefine, evenUpArray)
                   call writeInfo("Done the even up part", TRIVIAL)    
                endif
             if (myrankWorldGlobal == 1)call tune(6, "Refine grid")
!                iRefine = 0
!             endif

!             call writeVtkFile(grid, "afterrefine.vtk", &
!                  valueTypeString=(/"rho          ","velocity     ","rhoe         " , &
!                  "u_i          ", &
!                  "hydrovelocity", &
!                  "rhou         ", &
!                  "rhov         ", &
!                  "rhow         ", &
!                  "phi          ", &
!                  "phigas       ", &
!                  "pressure     ", &
!                  "q_i          "/))

             if (myrankWorldGlobal == 1) call tune(6, "Loop refine")
             
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
             
             if (myrankWorldGlobal == 1) call tune(6, "Loop refine")                  
          end if
       endif

       call sendSinksToZerothThread(globalnSource, globalsourceArray)

       if (doSelfGrav) then
          call zeroSourcepotential(grid%octreeRoot)
          if (globalnSource > 0) then
             call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
          endif
          call sumGasStarGravity(grid%octreeRoot)
       endif


!update the simulation time
       currentTime = currentTime + dt
       if (nbodyPhysics) globalSourceArray(1:globalnSource)%time = currentTime
       if (myRankWorldGlobal == 1) write(*,*) "current time ",currentTime,dt,nextDumpTime
       if (myRankWorldGlobal == 1) write(*,*) "percent to next dump ",100.*(nextDumpTime-currentTime)/tdump

!dump simulation data if simulation time is greater than next dump time
       if (currentTime .ge. nextDumpTime) then
          nextDumpTime = nextDumpTime + tDump
          it = it + 1

          grid%iDump = it
          grid%currentTime = currentTime

          if(mod(real(grid%iDump), real(vtutogrid)) == 0.0) then
             write(plotfile,'(a,i4.4,a)') "dump_",it,".grid"
             call writeAMRgrid(plotfile,.false. ,grid)
          end if

          if (dumpRadial) then
             write(plotfile,'(a,i4.4,a)') "radial",it,".dat"
             call  dumpValuesAlongLine(grid, plotfile, VECTOR(0.d0,0.d0,0.0d0), &
                  VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0),1000)
          endif


          if (writeoutput) then
             write(plotfile,'(a,i4.4,a)') "source",it,".dat"
             call writeSourceArray(plotfile)
          endif

          write(plotfile,'(a,i4.4,a)') "dump_",it,".vtk"
          call writeVtkFile(grid, plotfile, &
               valueTypeString=(/"rho          ",&
               "hydrovelocity", &
               "rhoe         ", &  
               "rhou         ", &  
               "rhov         ", &  
               "rhow         ", &  
               "u_i          ", &
               "mpithread    ", &
               "phi          "/))

          write(plotfile,'(a,i4.4,a)') "nbody",it,".vtk"
          call writeVtkFilenBody(globalnSource, globalsourceArray, plotfile)
          if (myrankGlobal==1) write(*,*) trim(plotfile), " written at ",currentTime/tff, " free-fall times"

          if (writeoutput.and.(globalnSource>0)) &
               write(45,'(i6,1p,20e15.5,0p)') it, currentTime, globalSourceArray(1:globalnSource)%position


       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)

       if (currentTime  == tEnd) exit
    enddo
    if (Writeoutput.and.(globalnSource >0)) close(45)
666 continue
  end subroutine doHydrodynamics3d

!The main routine for 2D hydrodynamics
  subroutine doHydrodynamics2d(grid)
    use inputs_mod, only : tEnd, tDump, doRefine, doUnrefine, amrTolerance, modelwashydro
    use inputs_mod, only : vtutogrid
    use mpi
    type(gridtype) :: grid
    real(double) :: dt, tc(512), temptc(512), mu
    real(double) :: currentTime
    integer :: i, it, iUnrefine
    character(len=80) :: plotfile, filename
    real(double) :: nextDumpTime, tff!, ang
    real(double) :: totalEnergy, totalMass, tempdouble, dt_pressure, dt_viscous
    type(VECTOR) :: direction, viewVec
    integer :: thread1(512), thread2(512), nBound(512), nPairs
    integer :: nGroup, group(512)
!    integer :: cornerthread1(100), cornerthread2(100), ncornerBound(100), ncornerPairs
!    integer :: nCornerGroup, cornergroup(100)
    integer :: nHydroThreads 
!    logical :: converged
    integer :: nUnrefine, jt, count
    integer :: evenUpArray(nHydroThreadsGlobal)
    integer :: ierr
    logical :: openFile

    nUnrefine = 0
    count = 0
    nHydroThreads = nHydroThreadsGlobal

    it = grid%iDump
    currentTime = grid%currentTime
    nextDumpTime = 0.d0

    if (it /= 1) then
       call writeVTKfile(grid, "readin.vtk")
    endif
    if (myrankGlobal /= 0) then

       !famr
!       call refineEdges(grid%octreeRoot, grid,  converged)

       direction = VECTOR(1.d0, 0.d0, 0.d0)
       mu = 2.d0
!       if(grid%geometry == "SB_gasmix") then
!          mu = 1.d0
 !      end if


       viewVec = VECTOR(-1.d0,0.d0,0.d0)
       !    viewVec = rotateZ(viewVec, 20.d0*degtorad)
       viewVec = rotateY(viewVec, 25.d0*degtorad)

!determine which mpi threads are in contact with one another (i.e. share a domain boundary)
       if (myRankWorldGlobal == 1) write(*,*) "CFL set to ", cflNumber
       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       do i = 1, nPairs
          if (myrankglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
       enddo




!do initial exchange across boundaries. The exchange gives subdomain boundary cells information about their foreign neighbours
       call writeInfo("Doing initial exchange", TRIVIAL)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call writeInfo("Done", TRIVIAL)


!       call returnCornerPairs(grid, nCornerPairs, cornerthread1, cornerthread2, nCornerBound, cornerGroup, nCornerGroup)
!       do i = 1, nCornerPairs
!          if (myrankglobal==1)write(*,*) "Corner pair ", i, cornerthread1(i), " -> ", cornerthread2(i), " bound ", nCornerbound(i)
!       end do
!

!       call writeInfo("Calling exchange across boundary corners", TRIVIAL)
!       call exchangeAcrossMPICorner(grid, nCornerPairs, cornerThread1, cornerThread2, nCornerBound, cornerGroup, nCornerGroup)
!       call writeInfo("Done", TRIVIAL)

       call writeVTKfile(grid, "start1.vtk")
       call writeVtkFile(grid, "start1.vtk", &
            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
            ,"mpithread    ", "pressure     " /))

!set up initial values
       if (it == 0) then
          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call calculateRhoU(grid%octreeRoot, direction)
!          direction = VECTOR(0.d0, 1.d0, 0.d0)
!          call calculateRhoV(grid%octreeRoot, direction)
          direction = VECTOR(0.d0, 0.d0, 1.d0)
          call calculateRhoW(grid%octreeRoot, direction)
          call calculateRhoE(grid%octreeRoot, direction)
          
!ensure that all cells are within one level of refinement of one another
!          call writeInfo("Evening up", TRIVIAL)
!          call evenUpGridMPI(grid,.true., dorefine, evenUpArray)
!          call writeInfo("Done", TRIVIAL)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


    !      if(doRefine) then
    !         call setAllUnchanged(grid%octreeRoot)
    !         call writeInfo("Refining", TRIVIAL)
    !         call refinegridGeneric(grid, amrTolerance, evenuparray)          
    !         call writeInfo("Done", TRIVIAL)
    !      end if    
          call writeInfo("Evening up", TRIVIAL)
          call evenUpGridMPI(grid, .true.,dorefine, evenUpArray)
          call writeInfo("Done", TRIVIAL)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call setupX(grid%octreeRoot, grid, direction)
          call setupQX(grid%octreeRoot, grid, direction)
!          call computepressureGeneral(grid, grid%octreeroot, .false.)
          
          call writeInfo("Checking Boundary Partner Vectors", TRIVIAL)
          call checkBoundaryPartners(grid%octreeRoot, grid)
          call writeInfo("Initial Boundary Partner Check Passed", TRIVIAL)
       endif
    endif

    if(modelWasHydro) then
       print *, "MODEL WAS HYDRO"
       call populateHydroVelWithCornerVel(grid%octreeRoot, grid)
       call writeVtkFile(grid, "post_corners.vtk", &
            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " /))
    end if

       call writeVTKfile(grid, "start2.vtk")
       call writeVtkFile(grid, "start2.vtk", &
            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
            ,"mpithread    ", "pressure     " /))

!calculate largest timestep that each cell on the grid can take without advecting a quantity further than their
!nearest neighbour
    tc = 0.d0
    dt_pressure = 1.d30
       tc = 0.d0
       dt_pressure = 1.d30
       dt_viscous = 1.d30
       if (myrankGlobal /= 0) then
          tc(myrankGlobal) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRankGlobal))
          call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
          call computeDivV(grid%octreeRoot, grid)
          call pressureGradientTimeStep(grid, dt_pressure, npairs,thread1,thread2,nbound,group,ngroup)
          call viscousTimescale(grid%octreeRoot, grid, dt_viscous)
       endif
       call MPI_ALLREDUCE(dt_pressure, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_pressure = tempDouble
       call MPI_ALLREDUCE(dt_viscous, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_viscous = tempDouble
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads))
       dt = MIN(dt_pressure, dt)
       dt = MIN(dt_viscous, dt)
       if (writeoutput) then
          write(*,*) "Courant time from v ",dt
          write(*,*) "Courant time from P ",dt_pressure
          write(*,*) "Courant time from visc ",dt_viscous
       endif
       dt = dt * dble(cflNumber)

    tff = 1.d0 / sqrt(bigG * (0.1d0*mSol/((4.d0/3.d0)*pi*7.d15**3)))

    if (writeoutput) write(*,*) "Setting tdump to: ", tdump

    if (it /= 0) then
       nextDumpTime = grid%currentTime + tdump
    endif

    iUnrefine = 0

    jt = 0

    !trace courant time history
    open (444, file="tcHistory.dat", status="unknown")

!loop until simulation end time
    do while(currentTime < tEnd)

       if (myrankWorldGlobal == 1) write(*,*) "current time " ,currentTime
       jt = jt + 1
       tc = 0.d0

       tc = 0.d0
       dt_pressure = 1.d30
       dt_viscous = 1.d30
       if (myrankGlobal /= 0) then
          tc(myrankGlobal) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRankGlobal))
          call pressureGradientTimeStep(grid, dt_pressure, npairs,thread1,thread2,nbound,group,ngroup)
          call viscousTimescale(grid%octreeRoot, grid, dt_viscous)
       endif
       call MPI_ALLREDUCE(dt_pressure, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_pressure = tempDouble
       call MPI_ALLREDUCE(dt_viscous, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_viscous = tempDouble
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads))
       dt = MIN(dt_pressure, dt)
       dt= MIN(dt_viscous, dt)
       dt = dt * dble(cflNumber)
       if (writeoutput) then
          write(*,*) "Courant time from v ",dt
          write(*,*) "Courant time from P ",dt_pressure
          write(*,*) "Courant time from visc ",dt_viscous
       endif

       write(444, *) jt, MINVAL(tc(1:nHydroThreads)), dt

       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

       if ((currentTime + dt).gt.tEnd) then
          nextDumpTime = tEnd
          dt = nextDumpTime - currentTime
       endif



       if (myrankGlobal /= 0) then
          if (myrankWorldGlobal == 1) write(*,*) "courantTime", dt, it
          if (myrankWorldGlobal == 1) call tune(6,"Hydrodynamics step")
          call writeInfo("calling hydro step",TRIVIAL)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!perform a hydrodynamics step in the x and z directions
          call hydroStep2d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)!, nCornerPairs, cornerThread1, cornerThread2, &
!             nCornerBound, cornerGroup, nCornerGroup)
       end if
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!get total mass/energy on the grid, for diagnostics
       call findEnergyOverAllThreads(grid, totalenergy)
       if (writeoutput) write(*,*) "Total energy: ",totalEnergy
       call findMassOverAllThreads(grid, totalmass)
       if (writeoutput) write(*,*) "Total mass: ",totalMass


!unrefine the grid where necessary
       if(doUnrefine) then
          iUnrefine = iUnrefine + 1
          if (iUnrefine == 5) then


!          call writeVtkFile(grid, "beforeunrefine.vtk", &
!               valueTypeString=(/"rho          ", &
!                                 "mpistore     "/))

             if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
             call unrefineCells(grid%octreeRoot, grid, nUnrefine, amrtolerance)
             if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
             iUnrefine = 0
          endif
       end if
    
       if (myrankWorldGlobal == 1) call tune(6,"Hydrodynamics step")

       call evenUpGridMPI(grid, .true., dorefine, evenUpArray) !, dumpfiles=jt)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!not all values survive the refine - therefore need to re-add them if wanting to dump a vtk file
       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call setupX(grid%octreeRoot, grid, direction)
       call setupQX(grid%octreeRoot, grid, direction)
       call computepressureGeneral(grid, grid%octreeroot, .false.)

!refine the grid where necessary
       if(doRefine) then
          call setAllUnchanged(grid%octreeRoot)
          call refinegridGeneric(grid, amrTolerance, evenuparray)
       end if

       call evenUpGridMPI(grid, .true., dorefine, evenUpArray) !, dumpfiles=jt)
       
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       
       currentTime = currentTime + dt

       if (nbodyPhysics) globalSourceArray(1:globalnSource)%time = currentTime

       !Perform another boundary partner check
       call checkBoundaryPartners(grid%octreeRoot, grid)

!dump simulation data if the simulation time exceeds the next dump time
       if (currentTime .ge. nextDumpTime) then
          it = it + 1
          nextDumpTime = nextDumpTime + tDump
          grid%iDump = it
          grid%currentTime = currentTime
          call hydrovelocityConvert(grid%octreeRoot)

          if(mod(real(grid%iDump), real(vtutogrid)) == 0.0) then
             write(plotfile,'(a,i4.4,a)') "dump_",it,".grid"
             call writeAMRgrid(plotfile,.false. ,grid)
          end if
          write(plotfile,'(a,i4.4,a)') "dump_",it,".vtk"
          call writeVtkFile(grid, plotfile, &
               valueTypeString=(/"rho          ","velocity     ","rhoe         " , &
               "u_i          ", &
               "hydrovelocity", &
               "chiline      ", &
               "rhou         ", &
               "rhov         ", &
               "vrot         ", &
               "rhow         ", &
               "phi          ", &
               "pressure     ", &
               "q11          ", &
               "q22          ", &
               "q_i          "/))
               
          if(grid%geometry == "SB_CD_2Da" .or. grid%geometry == "SB_CD_2Db" .or. &
               grid%geometry == "SB_offCentre") then

             write(filename,'(a,i4.4,a)') "dump_",it,".txt"
             openFile = .true.

             if(myRankGlobal == 0) then
                call writePosRhoPressureVelZERO(filename, grid)
             else
                call writePosRhoPressureVel(grid%octreeRoot)
                call killZero()
             end if
          end if

          if (grid%geometry == "sedov") &
               call dumpValuesAlongLine(grid, "sedov.dat", VECTOR(0.5d0,0.d0,0.0d0), VECTOR(0.9d0, 0.d0, 0.0d0), 1000)
       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)
    enddo
    close(444)
  end subroutine doHydrodynamics2d

!The main routine for 2D hydrodynamics in cylindrical coordinates
  subroutine doHydrodynamics2dCylindrical(grid)
    use inputs_mod, only : tEnd, tDump, doRefine, doUnrefine, amrTolerance, modelwashydro, doselfgrav
    use inputs_mod, only : vtutogrid, advectHydro
    use mpi
    type(gridtype) :: grid
    real(double) :: dt, tc(512), temptc(512), mu
    real(double) :: currentTime
    integer :: i, it, iUnrefine
    character(len=80) :: plotfile
    real(double) :: nextDumpTime, tff!, ang
    real(double) :: totalEnergy, totalMass, tempdouble, dt_pressure, dt_viscous, vBulk, vSound, dt_grav, dt_nbody
    type(VECTOR) :: initAngMom
    
    type(VECTOR) :: direction, viewVec, totalAngMom
    integer :: thread1(512), thread2(512), nBound(512), nPairs
    integer :: nGroup, group(512)
!    integer :: cornerthread1(100), cornerthread2(100), ncornerBound(100), ncornerPairs
!    integer :: nCornerGroup, cornergroup(100)
    integer :: nHydroThreads 
!    logical :: converged
    integer :: nUnrefine, jt, count
    integer :: evenUpArray(nHydroThreadsGlobal)
    integer :: ierr

    nUnrefine = 0
    count = 0
    nHydroThreads = nHydroThreadsGlobal

    it = grid%iDump
    currentTime = grid%currentTime
    nextDumpTime = 0.d0


    if (it /= 1) then
       call writeVTKfile(grid, "readin.vtk")
    endif
    if (myrankGlobal /= 0) then


       !famr
!       call refineEdges(grid%octreeRoot, grid,  converged)

       direction = VECTOR(1.d0, 0.d0, 0.d0)
       mu = 2.d0
!       if(grid%geometry == "SB_gasmix") then
!          mu = 1.d0
 !      end if


       viewVec = VECTOR(-1.d0,0.d0,0.d0)
       !    viewVec = rotateZ(viewVec, 20.d0*degtorad)
       viewVec = rotateY(viewVec, 25.d0*degtorad)

!determine which mpi threads are in contact with one another (i.e. share a domain boundary)
       if (myRankWorldGlobal == 1) write(*,*) "CFL set to ", cflNumber
       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       do i = 1, nPairs
          if (myrankglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
       enddo


       call writeInfo("Setting up even up array", TRIVIAL)
       call setupEvenUpArray(grid, evenUpArray)
       call writeInfo("Done", TRIVIAL)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call writeInfo("Evening up", TRIVIAL)
       call evenUpGridMPI(grid, .true.,dorefine, evenUpArray)
       call writeInfo("Done", TRIVIAL)



!do initial exchange across boundaries. The exchange gives subdomain boundary cells information about their foreign neighbours
       call writeInfo("Doing initial exchange", TRIVIAL)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call writeInfo("Done", TRIVIAL)


!       call returnCornerPairs(grid, nCornerPairs, cornerthread1, cornerthread2, nCornerBound, cornerGroup, nCornerGroup)
!       do i = 1, nCornerPairs
!          if (myrankglobal==1)write(*,*) "Corner pair ", i, cornerthread1(i), " -> ", cornerthread2(i), " bound ", nCornerbound(i)
!       end do
!

!       call writeInfo("Calling exchange across boundary corners", TRIVIAL)
!       call exchangeAcrossMPICorner(grid, nCornerPairs, cornerThread1, cornerThread2, nCornerBound, cornerGroup, nCornerGroup)
!       call writeInfo("Done", TRIVIAL)

       call writeVTKfile(grid, "start1.vtk")
       call writeVtkFile(grid, "start1.vtk", &
            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
            ,"mpithread    ", "pressure     " /))

!set up initial values
       if (it == 0) then
          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call calculateRhoU(grid%octreeRoot, direction)
          direction = VECTOR(0.d0, 1.d0, 0.d0)
!          call calculateRhoRV(grid%octreeRoot, direction)
          direction = VECTOR(0.d0, 0.d0, 1.d0)
          call calculateRhoW(grid%octreeRoot, direction)
          call calculateRhoE(grid%octreeRoot, direction)
          
          
!ensure that all cells are within one level of refinement of one another
!          call writeInfo("Evening up", TRIVIAL)
!          call evenUpGridMPI(grid,.true., dorefine, evenUpArray)
!          call writeInfo("Done", TRIVIAL)

          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


    !      if(doRefine) then
    !         call setAllUnchanged(grid%octreeRoot)
    !         call writeInfo("Refining", TRIVIAL)
    !         call refinegridGeneric(grid, amrTolerance, evenuparray)          
    !         call writeInfo("Done", TRIVIAL)
    !      end if    
          call writeInfo("Evening up", TRIVIAL)
          call evenUpGridMPI(grid, .true.,dorefine, evenUpArray)
          call writeInfo("Done", TRIVIAL)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call setupX(grid%octreeRoot, grid, direction)
          call setupQX(grid%octreeRoot, grid, direction)
!          call computepressureGeneral(grid, grid%octreeroot, .false.)
          
!          call unsetGhosts(grid%octreeRoot)
!          call createGhostCells(grid)
          call writeVtkFile(grid, "ghost.vtk", &
            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
            ,"mpithread    ", "pressure     ","ghosts       ","edges        " /))

          call writeInfo("Checking Boundary Partner Vectors", TRIVIAL)
          call checkBoundaryPartners(grid%octreeRoot, grid)
          call writeInfo("Initial Boundary Partner Check Passed", TRIVIAL)

       if (doselfgrav) then
          if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
          call zeroPhiGas(grid%octreeRoot)
          if (dogasgravity) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.)

          call zeroSourcepotential(grid%octreeRoot)
          call sumGasStarGravity(grid%octreeRoot)
          if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
       endif

       endif
    endif



    if(modelWasHydro) then
       print *, "MODEL WAS HYDRO"
       call populateHydroVelWithCornerVel(grid%octreeRoot, grid)
       call writeVtkFile(grid, "post_corners.vtk", &
            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " /))
    end if

       call writeVTKfile(grid, "start2.vtk")
       call writeVtkFile(grid, "start2.vtk", &
            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
            ,"mpithread    ", "pressure     ","ghosts       ","edges        " /))

       if (doSelfGrav .and. grid%geometry == "gravtest") then

          write(plotfile,'(a,i4.4,a)') "radial",1,".dat"
          call  dumpValuesAlongLine(grid, plotfile, VECTOR(0.d0,0.d0,0.0d0), &
               VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0),1000)
           write(plotfile,'(a,i2.2,a)') "gravfull.vtk"
           call writeVtkFile(grid, plotfile, &
                valueTypeString=(/"phigas ", "rho    ","chiline","adot   "/))


          goto 666

       end if


       if (myrankGlobal /=0) then
          if (myrankWorldglobal == 1) call tune(6,"Alpha viscosity")
          call setupAlphaViscosity(grid, alphaViscosity, 0.1d0)
          if (myrankWorldglobal == 1) call tune(6,"Alpha viscosity")
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call setupCylindricalViscosity(grid%octreeRoot, grid)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       endif


!calculate largest timestep that each cell on the grid can take without advecting a quantity further than their
!nearest neighbour
    tc = 0.d0
    dt_pressure = 1.d30
       tc = 0.d0
       dt_pressure = 1.d30
       dt_viscous = 1.d30
       dt_nbody = 1.d30
       if (myrankGlobal /= 0) then
          tc(myrankGlobal) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRankGlobal))
          call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
          if (includePressureTerms) call pressureGradientTimeStep(grid, dt_pressure, npairs,thread1,thread2,nbound,group,ngroup)
          call viscousTimescaleCylindrical(grid%octreeRoot, grid, dt_viscous)
          if (nBodyPhysics) call computeCourantTimeNbody(grid, GlobalnSource, globalsourceArray, dt_nbody)
       endif
       call MPI_ALLREDUCE(dt_pressure, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_pressure = tempDouble
       call MPI_ALLREDUCE(dt_nbody, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_nbody = tempDouble
       call MPI_ALLREDUCE(dt_viscous, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_viscous = tempDouble
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads))
       dt = MIN(dt_pressure, dt)
       dt = MIN(dt_viscous, dt)
       dt = MIN(dt_nbody, dt)
       if (writeoutput) then
          write(*,*) "Courant time from v ",dt
          write(*,*) "Courant time from P ",dt_pressure
          write(*,*) "Courant time from visc ",dt_viscous
       endif
       dt = dt * dble(cflNumber)

    tff = 1.d0 / sqrt(bigG * (0.1d0*mSol/((4.d0/3.d0)*pi*7.d15**3)))

    if (writeoutput) write(*,*) "Setting tdump to: ", tdump

    if (it /= 0) then
       nextDumpTime = grid%currentTime + tdump
    endif

    iUnrefine = 0

    jt = 0
    call findAngMomOverAllThreads(grid, initAngMom, VECTOR(0.d0, 0.d0, 0.d0))
    initAngMom%z = initAngMom%z + globalSourceArray(1)%angMomentum%z

    !trace courant time history
    open (444, file="tcHistory.dat", status="unknown")

!loop until simulation end time
    do while(currentTime < tEnd)


       if (myrankWorldGlobal == 1) write(*,*) "current time " ,currentTime
       jt = jt + 1
       tc = 0.d0

       dt_pressure = 1.d30
       dt_viscous = 1.d30
       dt_grav = 1.d30
       vBulk = 0.d0
       vSound = 0.d0
       if (myrankGlobal /= 0) then
          tc(myrankGlobal) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRankGlobal))
          call computeCourantV(grid, grid%octreeRoot, vBulk, vSound)
          call pressureGradientTimeStep(grid, dt_pressure, npairs,thread1,thread2,nbound,group,ngroup)
          call viscousTimescaleCylindrical(grid%octreeRoot, grid, dt_viscous)
          call courantTimeGasSourceCylindrical(grid, grid%octreeRoot, globalnsource, globalsourcearray, dt_grav)
       endif

       call MPI_ALLREDUCE(vBulk, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MAX, localWorldCommunicator, ierr)
       vBulk = tempDouble
       call MPI_ALLREDUCE(vSound, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MAX, localWorldCommunicator, ierr)
       vSound = tempDouble
       call MPI_ALLREDUCE(dt_pressure, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_pressure = tempDouble
       call MPI_ALLREDUCE(dt_viscous, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_viscous = tempDouble
       call MPI_ALLREDUCE(dt_grav, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       dt_grav = tempDouble

       call MPI_ALLREDUCE(tc, tempTc, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)
       tc = tempTc
       dt = MINVAL(tc(1:nHydroThreads))
       dt = MIN(dt_pressure, dt)
       dt = MIN(dt_viscous, dt)
       dt = MIN(dt_grav, dt)
       dt = dt * dble(cflNumber)
       if (writeoutput) then
          write(*,*) "Courant time from v ",MINVAL(tc(1:nHydroThreads))
          write(*,*) "maximum bulk speed ",Vbulk/1.d5
          write(*,*) "maximum sound speed ",Vsound/1.d5
          write(*,*) "Courant time from P ",dt_pressure
          write(*,*) "Courant time from visc ",dt_viscous
          write(*,*) "Courant time from grav source on gas ",dt_grav
       endif

       write(444, *) jt, MINVAL(tc(1:nHydroThreads)), dt

       if (.not.advectHydro) dt = tdump
       if ((currentTime + dt).gt.tEnd) then
          nextDumpTime = tEnd
          dt = nextDumpTime - currentTime
       endif

       if ((currentTime + dt) .gt. nextDumpTime) then
          dt = nextDumpTime - currentTime
       endif

       if (myrankGlobal /= 0) then
          if (myrankWorldGlobal == 1) write(*,*) "courantTime", dt, it
          if (myrankWorldGlobal == 1) call tune(6,"Hydrodynamics step")
          call writeInfo("calling hydro step",TRIVIAL)



          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


!perform a hydrodynamics step in the x and z directions
!          call hydroStep2dCylindrical(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)

!       call writeVtkFile(grid, "beforestep.vtk", &
!            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
!            ,"mpithread    ", "pressure     ","ghosts       ","edges        " /))

          call hydroStep2dCylindrical_amr(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)

!       call writeVtkFile(grid, "afterstep.vtk", &
!            valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
!            ,"mpithread    ", "pressure     ","ghosts       ","edges        " /))

!          call calculateTemperatureFromEnergy(grid%octreeRoot)
       end if


       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!get total mass/energy on the grid, for diagnostics
       call findEnergyOverAllThreads(grid, totalenergy)
       if (writeoutput) write(*,*) "Total energy: ",totalEnergy
       call findMassOverAllThreads(grid, totalmass)
       if (writeoutput) write(*,*) "Total mass: ",totalMass
       if (cylindricalHydro) then
          call findAngMomOverAllThreads(grid, totalAngMom, VECTOR(0.d0, 0.d0, 0.d0))
          if (writeoutput) write(*,*) "Total angular momentum: ",totalAngMom%z+globalSourceArray(1)%angMomentum%z, &
               100.d0*(totalAngMom%z+globalSourceArray(1)%angMomentum%z-initAngMom%z)/initAngMom%z
       endif


!unrefine the grid where necessary
       if(doUnrefine) then
          iUnrefine = iUnrefine + 1
          if (iUnrefine == 5) then


!          call writeVtkFile(grid, "beforeunrefine.vtk", &
!               valueTypeString=(/"rho          ", &
!                                 "mpistore     "/))

             if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
             call unrefineCells(grid%octreeRoot, grid, nUnrefine, amrtolerance)
             if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
             iUnrefine = 0
          endif
       end if
    
       if (myrankWorldGlobal == 1) call tune(6,"Hydrodynamics step")

       call evenUpGridMPI(grid, .true., dorefine, evenUpArray) !, dumpfiles=jt)

!       call unsetGhosts(grid%octreeRoot)
!       call createGhostCells(grid)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!not all values survive the refine - therefore need to re-add them if wanting to dump a vtk file
       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call setupX(grid%octreeRoot, grid, direction)
       call setupQX(grid%octreeRoot, grid, direction)
       call computepressureGeneral(grid, grid%octreeroot, .false.)

!refine the grid where necessary
       if(doRefine) then
          call setAllUnchanged(grid%octreeRoot)
          call refinegridGeneric(grid, amrTolerance, evenuparray)
       end if

       call evenUpGridMPI(grid, .true., dorefine, evenUpArray) !, dumpfiles=jt)
!       call unsetGhosts(grid%octreeRoot)
!       call createGhostCells(grid)

       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       
       currentTime = currentTime + dt
       if (nbodyPhysics) globalSourceArray(1:globalnSource)%time = currentTime

       if (writeoutput) write(*,'(a,f7.2)') "Percent left until dump ",100.d0*(nextDumpTime - currentTime)/tdump

       !Perform another boundary partner check
       call checkBoundaryPartners(grid%octreeRoot, grid)



!dump simulation data if the simulation time exceeds the next dump time
       if (currentTime .ge. nextDumpTime) then
          it = it + 1
          nextDumpTime = nextDumpTime + tDump
          grid%iDump = it
          grid%currentTime = currentTime
          call hydrovelocityConvert(grid%octreeRoot)
          if(mod(real(grid%iDump), real(vtutogrid)) == 0.0) then
             write(plotfile,'(a,i4.4,a)') "dump_",it,".grid"
             call writeAMRgrid(plotfile,.false. ,grid)
          end if
          write(plotfile,'(a,i4.4,a)') "dump_",it,".vtk"
          call writeVtkFile(grid, plotfile, &
               valueTypeString=(/"rho          ","velocity     ","rhoe         " , &
               "u_i          ", &
               "hydrovelocity", &
               "chiline      ", &
               "mass         ", &
               "rhou         ", &
               "rhov         ", &
               "vrot         ", &
               "rhow         ", &
               "phi          ", &
               "pressure     ", &
               "ghosts       ", &
               "edges        ", &
               "q11          ", &
               "q22          ", &
               "temperature  ", &
               "fvisc1       ", &
               "fvisc2       ", &
               "fvisc3       ", &
               "q_i          "/))

!          write(plotfile,'(a,i4.4,a)') "visc_",it,".vtk"
!          call writeVtkFile(grid, plotfile, &
!               valueTypeString=(/"q11          ", &
!                                 "q12          ", &
!                                 "q13          ", &
!                                 "q21          ", &
!                                 "q22          ", &
!                                 "q23          ", &
!                                 "q31          ", &
!                                 "q32          ", &
!                                 "q33          ", &
!                                 "fvisc1       ", &
!                                 "vphi         ", &
!                                 "etaline      ", &
!                                 "fvisc2       ", &
!                                 "fvisc3       ", &
!                                 "rho          ", &
!                                 "rhou         ", &
!                                 "rhov         ", &
!                                 "rhow         " &
!                                 /))

          if (writeoutput) then
             write(*,*) "rank ",myrankglobal, " writing source array"
             write(plotfile,'(a,i4.4,a)') "source",it,".dat"
             call writeSourceArray(plotfile)
          endif
               
          write(plotfile,'(a,i4.4,a)') "nbody",it,".vtk"
          call writeVtkFilenBody(globalnSource, globalsourceArray, plotfile)


          if (dumpRadial) then
             write(plotfile,'(a,i4.4,a)') "radial",it,".dat"
             call  dumpValuesAlongLine(grid, plotfile, VECTOR(0.d0,0.d0,0.0d0), &
                  VECTOR(grid%octreeRoot%subcellSize*1.999999d0, 0.d0, 0.d0),1000)
          endif

          if (grid%geometry == "sedov") &
               call dumpValuesAlongLine(grid, "sedov.dat", VECTOR(0.5d0,0.d0,0.0d0), VECTOR(0.9d0, 0.d0, 0.0d0), 1000)
       endif
       viewVec = rotateZ(viewVec, 1.d0*degtorad)
    enddo
    close(444)
666 continue
  end subroutine doHydrodynamics2dCylindrical


  recursive  subroutine assignRhoeLast(thisOctal)
    use mpi
    type(octal), pointer :: thisOctal, child
    integer :: subcell
    integer :: i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call assignRhoeLast(child)
                exit
             end if
          end do
       else 
          thisOctal%rhoeLastTime(subcell) = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
       end if
    end do
  end subroutine assignRhoeLast


  recursive  subroutine writePosRhoPressureVel(thisOctal)
    use mpi
    type(octal), pointer :: thisOctal, child
    integer :: subcell
    integer :: ierr, i
    integer, parameter :: nStorage = 8
    real(double) :: tempstorage(nStorage)
    type(vector) :: rVec
    integer :: tag = 40

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call writePosRhoPressureVel(child)
                exit
             end if
          end do
       else 
          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
             rVec = subcellCentre(thisOctal, subcell)
             tempStorage(1) = rVec%x
             tempStorage(2) = rVec%y
             tempStorage(3) = rVec%z
             tempStorage(4) = thisOctal%rho(subcell)
             tempStorage(5) = thisOctal%pressure_i(subcell)
             tempStorage(6) = thisOctal%rhou(subcell)
             tempStorage(7) = thisOctal%rhov(subcell)
             tempStorage(8) = thisOctal%rhow(subcell)

             call mpi_send(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, &
                  tag, localWorldCommunicator, ierr)
          end if      
       end if
    end do   
 
  end subroutine writePosRhoPressureVel

  subroutine writePosRhoPressureVelZERO(fileName, grid)
    use mpi
    type(gridtype) :: grid
    character(len=*) :: fileName
    integer :: ier, ierr
    integer, parameter :: nStorage = 8
    real(double) :: tempstorage(nStorage)
    integer :: tag = 40
    logical :: cycling
    integer :: status(MPI_STATUS_SIZE)

    open(unit=912, file=fileName, status="replace", iostat=ier)
    write (912, '(6(a4, 3x))') "x", "z", "rho", "p", "vx", "vz"

    cycling = .true.
    do while (cycling)
       call mpi_recv(tempStorage, nStorage, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, &
            localWorldCommunicator, status, ierr)           

       if(tempStorage(4) > 1.d29) then
          cycling = .false.
       end if
       
       if(cycling) then
          if(grid%octreeroot%twod) then
             write(912, '(6(e12.3, 3x))') tempStorage(1), tempStorage(3), tempStorage(4), &
                  tempStorage(5), tempStorage(6), tempStorage(8)
          else
             write(912, '(6(e12.3, 3x))') tempStorage(1), tempStorage(2), tempStorage(3), &
                  tempStorage(4), tempStorage(5), tempStorage(6), tempStorage(8)
          end if
       end if
    end do
    close(912)
  end subroutine writePosRhoPressureVelZERO

  subroutine killZero()
    use mpi    
    integer :: ierr
    integer, parameter :: nStorage = 8
    real(double) :: tempstorage(nStorage)
    integer :: tag = 40

    call MPI_BARRIER(amrCommunicator, ierr)
    if(myRankGlobal == 1) then
       tempStorage = 0.d0
       tempStorage(4) = 1.d30
       
       call mpi_send(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, &
            tag, localWorldCommunicator, ierr)       

    end if
       

  end subroutine
  

!clear the memory of what cells werer refined in the last refinement sweep 
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

!clear the memory of what cells werer refined in the last refinement sweep 
  recursive subroutine populateHydroVelWithCornerVel(thisOctal, grid, startOctal)
    use amr_mod, only : amrGridVelocity
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(octal), pointer, optional :: startOctal
    type(gridtype) :: grid
    integer :: subcell, i
    type(vector) :: velocity, position

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call populateHydroVelWithCornerVel(child, grid, startOctal)
                exit
             end if
          end do
       else 
          position = subcellCentre(thisOctal, subcell)
          velocity = amrGridVelocity(grid%octreeRoot, position, startOctal = startOctal, &
               actualSubcell = subcell, linearinterp = .false.)

          thisOctal%rhou(subcell) = velocity%x
          thisOctal%rhov(subcell) = velocity%y
          thisOctal%rhow(subcell) = velocity%z
       endif
    enddo
  end subroutine populateHydroVelWithCornerVel


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

!set the gas contribution to the gravitational potential across the grid to zero
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

!Update ghost cells with their boundary condition imposed values
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
                thisOctal%temperature(subcell) = real(thisOctal%tempStorage(subcell,9))
                if (thisOctal%temperature(subcell) == 0.0) then
                   write(*,*) "warning: temperature is zero in transfertempstorage"
                   write(*,*) "pos ",subcellCentre(thisOctal,subcell)
                   write(*,*) "boundary ",thisOctal%boundaryCondition(subcell)
                endif
!                thisOctal%phi_i(subcell) =  thisOctal%tempStorage(subcell,8)
             else
                thisOctal%phi_gas(subcell) = thisOctal%tempStorage(subcell,1)
             endif
          endif
       endif
    enddo
    if (associated(thisOCtal%tempStorage)) Deallocate(thisOctal%tempStorage)
  end subroutine transferTempStorage

!Update ghost cells with their boundary condition imposed values
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
!                thisOctal%phi_i(subcell) = thisOctal%tempStorage(subcell,1)
                thisOctal%phi_gas(subcell) = thisOctal%tempStorage(subcell,1)
             endif
          endif
       enddo
    endif
    if (associated(thisOCtal%tempStorage)) Deallocate(thisOctal%tempStorage)
  end subroutine transferTempStorageLevel

  recursive subroutine zeroCorrections(thisOctal, nDepth)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: nDepth
    integer :: subcell, i


    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call zeroCorrections(child, nDepth)
       end do
    else
       if (.not.associated(thisOctal%phi_gas_corr)) allocate(thisOctal%phi_gas_corr(1:thisOctal%maxChildren))
       do subcell = 1, thisOctal%maxChildren
          thisOctal%phi_gas_corr(subcell) = 0.d0
       enddo
    endif
  end subroutine zeroCorrections

  recursive subroutine swapCorrectionsforPhiGas(thisOctal, nDepth)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: nDepth
    integer :: subcell, i

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call swapCorrectionsForPhiGas(child, nDepth)
       end do
    else

       do subcell = 1, thisOctal%maxChildren
          if (.not.associated(thisOctal%etaLine)) allocate(thisOctal%etaline(1:thisOctal%maxChildren))
          thisOctal%etaLine(subcell) = thisOctal%phi_gas(subcell)
          thisOctal%phi_gas(subcell) = thisOctal%phi_gas_corr(subcell)
       enddo
    endif
  end subroutine swapCorrectionsforPhiGas

  recursive subroutine swapPhiGasforCorrections(thisOctal, nDepth)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: nDepth
    integer :: subcell, i

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call swapPhiGasForCorrections(child, nDepth)
       end do
    else

       do subcell = 1, thisOctal%maxChildren
          thisOctal%phi_gas_corr(subcell) = thisOctal%phi_gas(subcell)
          thisOctal%phi_gas(subcell) = thisOctal%etaLine(subcell)
       enddo
    endif
  end subroutine swapPhiGasforCorrections

  recursive subroutine applyCorrections(thisOctal, nDepth)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: nDepth
    integer :: subcell, i


    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call applyCorrections(child, nDepth)
       end do
    else

       do subcell = 1, thisOctal%maxChildren
          thisOctal%phi_gas(subcell) = thisOctal%phi_gas(subcell) - thisOctal%phi_gas_corr(subcell)
       enddo
    endif
  end subroutine applyCorrections

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


recursive subroutine returnVelocityVector(thisOctal, grid, position, velocity)
!For use in molecular_mod.

integer :: subcell
integer :: i
type(octal), pointer :: thisOctal, child
type(vector) :: position 
type(vector) :: velocity
type(gridtype) :: grid
real(double) :: rho


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call returnVelocityVector(child, grid, position, velocity)
                exit
             end if
          end do
       else 
          if(inSubcell(thisOctal, subcell, position)) then
             rho = thisOctal%rho(subcell)
             velocity = VECTOR(thisOctal%rhou(subcell)/rho, thisOctal%rhov(subcell)/rho, &
                  thisOctal%rhow(subcell)/rho)
          end if
       endif
    enddo

  end subroutine returnVelocityVector





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

!unclassify all ghost cells as ghosts
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
          thisOctal%boundaryPartner(subcell) = subcellCentre(thisOctal, subcell)
          thisOctal%boundaryCell(subcell) = .false.
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

!set up the values that will be transferred to ghost cells in transferTempStorage
  recursive subroutine imposeBoundary(thisOctal, grid, phi_gas)
    use inputs_mod, only : fixedRhoBound, rho_const
    type(octal), pointer   :: thisOctal, bOctal
    logical, optional :: phi_gas
    type(octal), pointer  :: child
    logical :: doPhiGas
    type(GRIDTYPE) :: grid
    integer :: subcell, i, bSubcell
    type(VECTOR) :: locator, dir, rVec
    real(double) :: machNumber, Pr, rhor, gamma,  posFactor!, distance

    gamma = 5.d0/3.d0
    machNumber = 2.d0
  
    doPhiGas = .false.
    if (PRESENT(phi_gas)) doPhiGas = phi_gas

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call imposeBoundary(child, grid, phi_gas=doPhiGas)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%ghostCell(subcell)) then

             if (.not.associated(thisOctal%tempStorage)) then
                allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:9))
                thisOctal%tempStorage = 0.d0
             endif
             if (doPhiGas) then
                thisOctal%tempStorage(subcell,1) = thisOctal%phi_gas(subcell)
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
                   thisOctal%tempStorage(subcell,9) = bOctal%temperature(bSubcell)

! NB confusion regarding 2d being x,z rather than x,y

                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (cylindricalHydro.and.doPhiGas) then
                         thisOctal%tempStorage(subcell,1) = thisOctal%phi_gas(subcell)
                        if (dir%x > 0.9d0) then
                           thisOctal%tempStorage(subcell,1) = bOctal%phi_gas(bsubcell)
                        endif
                      endif

                     if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
!                      if ((abs(dir%x) > 0.5d0).and.abs(dir%z) > 0.5d0) then
!                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
!                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
!                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
!                      endif
                   else if (thisOctal%threed) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%y) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = -bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if ((abs(dir%x) > 0.2d0).and.abs(dir%z) > 0.2d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = -bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
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
                   if (.not.doPhiGas) Thisoctal%tempstorage(subcell,1) = bOctal%rho(bSubcell)
                   thisOctal%tempStorage(subcell,2) = bOctal%rhoE(bSubcell)

                   thisOctal%tempStorage(subcell,6) = bOctal%energy(bSubcell)
                   thisOctal%tempStorage(subcell,7) = bOctal%pressure_i(bSubcell)
                   thisOctal%tempStorage(subcell,9) = bOctal%temperature(bSubcell)

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
                   thisOctal%tempStorage(subcell,9) = inflowTemp
                   
! NB confusion regarding 2d being x,z rather than x,y
                   
!                   print *, "enforcing in mom of ", inflowmomentum
!                   print *, "enforcing in speed of ", inflowspeed
!                   print *, "enforcing in rho of ", inflowrho
                   


                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = inflowmomentum
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

                case(7) !inflow that varies as a function of position along the axis
                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocal(locator, bOctal, bSubcell)
                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)

                   rVec = subcellCentre(thisOctal, subcell)
                   if(xSlope) then
                      posFactor = rVec%x-(grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize)
                   else if(ySlope) then
                      posFactor = rVec%y-(grid%octreeRoot%centre%y-grid%octreeRoot%subcellSize)
                   else if(zSlope) then
                      posFactor = rVec%z-(grid%octreeRoot%centre%z-grid%octreeRoot%subcellSize)
                   end if


                   Thisoctal%tempstorage(subcell,1) = inflowRho+rhograd*posFactor
                   thisOctal%tempStorage(subcell,2) = inflowRhoE+rhoEgrad*posFactor

                   thisOctal%tempStorage(subcell,6) = inflowEnergy+Egrad*posFactor
                   thisOctal%tempStorage(subcell,7) = inflowPressure+Pgrad*posFactor


                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = inflowMomentum+momgrad*posFactor
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = inflowMomentum+momgrad*posFactor
                      endif
                   else if (thisOctal%threed) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = inflowMomentum+momgrad*posFactor
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%y) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = inflowMomentum+momgrad*posFactor
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = inflowMomentum+momgrad*posFactor
                      endif
                   endif

                case (8) ! stellar wind - do nothing added in addstellarwind


                case(9) !gas mix
                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   
                  
                    call findSubcellLocal(locator, bOctal, bSubcell)

                   rVec = subcellcentre(thisoctal, subcell)

                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)

                   if(sqrt(dir%x**2) > 2.d0*grid%halfsmallestsubcell + 1.d-10) then
                      if(locator%x > rvec%x) then
                         locator%x = locator%x - (grid%halfsmallestsubcell + 1.d-10*grid%halfsmallestsubcell)

                      else if(locator%x < rVec%x) then
                         locator%x = locator%x + (grid%halfsmallestsubcell + 1.d-10*grid%halfsmallestsubcell)
                      end if


                      elseif(sqrt(dir%z**2) > 2.d0*grid%halfsmallestsubcell + 1.d-10) then
                      if(locator%z > rvec%z) then
                         locator%z = locator%z - (grid%halfsmallestsubcell + 1.d-10*grid%halfsmallestsubcell)

                      else if(locator%z < rVec%z) then
                         locator%z = locator%z + (grid%halfsmallestsubcell + 1.d-10*grid%halfsmallestsubcell)
                      end if
                   end if

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
                   thisOctal%tempStorage(subcell,9) = bOctal%temperature(bSubcell)

! NB confusion regarding 2d being x,z rather than x,y

                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (cylindricalHydro.and.doPhiGas) then
                         thisOctal%tempStorage(subcell,1) = thisOctal%phi_gas(subcell)
                        if (dir%x > 0.9d0) then
                           thisOctal%tempStorage(subcell,1) = bOctal%phi_gas(bsubcell)
                        endif
                      endif

                     if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
!                      if ((abs(dir%x) > 0.5d0).and.abs(dir%z) > 0.5d0) then
!                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
!                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
!                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
!                      endif
                   else if (thisOctal%threed) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%y) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if ((abs(dir%x) > 0.2d0).and.abs(dir%z) > 0.2d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
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

!set up the values that will be transferred to ghost cells in transferTempStorage
  recursive subroutine imposeBoundaryLevel(thisOctal, grid, nDepth, phi_gas)
    use inputs_mod, only : fixedRhoBound, rho_const
    type(octal), pointer   :: thisOctal, bOctal
    logical, optional :: phi_gas
    integer :: nDepth
    type(octal), pointer  :: child
    logical :: doPhiGas
    type(GRIDTYPE) :: grid
    integer :: subcell, i, bSubcell
    type(VECTOR) :: locator, dir, rVec
    real(double) :: machNumber, Pr, rhor, gamma,  posFactor

    gamma = 5.d0/3.d0
    machNumber = 2.d0
  
    doPhiGas = .false.
    if (PRESENT(phi_gas)) doPhiGas = phi_gas

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call  imposeBoundaryLevel(child, grid, nDepth, phi_gas)
       end do
       else
       do subcell = 1, thisOctal%maxChildren

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%ghostCell(subcell)) then

             if (.not.associated(thisOctal%tempStorage)) then
                allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:8))
                thisOctal%tempStorage = 0.d0
             endif
             if (doPhiGas) then
                thisOctal%tempStorage(subcell,1) = thisOctal%phi_gas(subcell)
             endif

             select case(thisOctal%boundaryCondition(subcell))
                case(1) ! reflecting
                   locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocalLevel(locator, bOctal, bSubcell, nDepth)

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
                      if (cylindricalHydro.and.doPhiGas) then
                         thisOctal%tempStorage(subcell,1) = thisOctal%phi_gas(subcell)
                        if (dir%x > 0.9d0) then
                           thisOctal%tempStorage(subcell,1) = bOctal%phi_gas(bsubcell)
                        endif
!                        if (dir%z > 0.9d0) then
!                           thisOctal%tempStorage(subcell,1) = bOctal%phi_gas(bsubcell)
!                        endif
                      endif



                     if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
!                      if ((abs(dir%x) > 0.5d0).and.abs(dir%z) > 0.5d0) then
!                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
!                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
!                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
!                      endif
                   else if (thisOctal%threed) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%y) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = -bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
                      endif
                      if ((abs(dir%x) > 0.2d0).and.abs(dir%z) > 0.2d0) then
                         thisOctal%tempStorage(subcell,3) = -bOctal%rhou(bSubcell)
                         thisOctal%tempStorage(subcell,4) = -bOctal%rhov(bSubcell)
                         thisOctal%tempStorage(subcell,5) = -bOctal%rhow(bSubcell)
                         thisOctal%tempStorage(subcell,8) = bOctal%phi_i(bsubcell)
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
                   call findSubcellLocalLevel(locator, bOctal, bSubcell, nDepth)

                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)
                   if (.not.doPhiGas) Thisoctal%tempstorage(subcell,1) = bOctal%rho(bSubcell)
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
                   call findSubcellLocalLEvel(locator, bOctal, bSubcell, nDepth)

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
                   call findSubcellLocalLevel(locator, bOctal, bSubcell, nDepth)

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

                case(7) !inflow that varies as a function of position along the axis
                                      locator = thisOctal%boundaryPartner(subcell)
                   bOctal => thisOctal
                   call findSubcellLocalLevel(locator, bOctal, bSubcell, nDepth)
                   dir = subcellCentre(bOctal, bSubcell) - subcellCentre(thisOctal, subcell)
                   call normalize(dir)

                   rVec = subcellCentre(thisOctal, subcell)
                   if(xSlope) then
                      posFactor = rVec%x-(grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize)
                   else if(ySlope) then
                      posFactor = rVec%y-(grid%octreeRoot%centre%y-grid%octreeRoot%subcellSize)
                   else if(zSlope) then
                      posFactor = rVec%z-(grid%octreeRoot%centre%z-grid%octreeRoot%subcellSize)
                   end if


                   Thisoctal%tempstorage(subcell,1) = inflowRho+rhograd*posFactor
                   thisOctal%tempStorage(subcell,2) = inflowRhoE+rhoEgrad*posFactor

                   thisOctal%tempStorage(subcell,6) = inflowEnergy+Egrad*posFactor
                   thisOctal%tempStorage(subcell,7) = inflowPressure+Pgrad*posFactor


                   if (thisOctal%twod.or.thisOctal%oneD) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = inflowMomentum+momgrad*posFactor
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = inflowMomentum+momgrad*posFactor
                      endif
                   else if (thisOctal%threed) then
                      if (abs(dir%x) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = inflowMomentum+momgrad*posFactor
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%y) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = inflowMomentum+momgrad*posFactor
                         thisOctal%tempStorage(subcell,5) = 0.d0
                      endif
                      if (abs(dir%z) > 0.9d0) then
                         thisOctal%tempStorage(subcell,3) = 0.d0
                         thisOctal%tempStorage(subcell,4) = 0.d0
                         thisOctal%tempStorage(subcell,5) = inflowMomentum+momgrad*posFactor
                      endif
                   endif

                case (8) ! stellar wind - do nothing added in addstellarwind
                   
                case DEFAULT
                   write(*,*) "Unrecognised boundary condition in impose boundary: ", thisOctal%boundaryCondition(subcell)
             end select
             
             if(fixedRhoBound) then
                thisOctal%tempStorage(subcell,1) = rho_const
             end if

          endif
      enddo
    endif
  end subroutine imposeBoundaryLevel



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
    logical, save :: firsttime=.true.
    logical :: tmp

    tmp = currentlyDoingHydroStep
    currentlyDoingHydroStep = .true.

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
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle


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
                   case(1, 5, 6, 9)
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
                      if(firsttime) then
                         write(*,*) "unknown boundary condition in setupghostcells1: ",thisOCtal%boundaryCondition(subcell)
                         write(*,*) "rank ",myRankGlobal, " mpi thread ", thisOctal%mpithread(subcell)
                         write(*,*) "subcell ",subcell
                         write(*,*) "all ", thisOctal%boundaryCondition(subcell)
                         firsttime = .false.
                      end if
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
             case(1, 4, 5, 6, 9)
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
    currentlyDoingHydroStep = tmp
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
    character(len=10) :: boundary
    logical :: tmp

    tmp = currentlyDoingHydroStep
    if (myrankGlobal == 0) goto 666

    currentlyDoingHydroStep = .true.

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
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
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

!             if(rVec%x < 1.d8 .and. iProbe == 2) then
!                print *, "EDGE LOC IS ", locator
!                print *, "rVEC is ", rVec
!             end if
!

!             if(grid%octreeroot%oned .and. locator%x < 0.d0) then
!                nProbeOutside =  1
!                   if (iProbe == 1) then
!                      boundary = "xplus"
!                   else
!                      boundary = "xminus"
!                   endif
!             end if

!             print *, "checking ", locator
!             print *, "pos should be ", locator%x-6031250.0
             if (.not.inOctal(grid%octreeRoot, locator)) then
                nProbeOutside = nProbeOutside + 1
                thisOctal%boundaryPartner(subcell) = (-1.d0)*probe(iProbe)


                if (thisOctal%oneD) then
                   if (iProbe == 1) then
                      boundary = "xplus"
!                      print *, "found xplus"
                   else
                      boundary = "xminus"
!                      print *, "found xminus"
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
!             print *, "found an edge ", subcellCentre(thisOctal, subcell)
!          else if(rVec%x < 1.d1) then
!             print *, "EDGE SHOULD HAVE BEEN FOUND"
                
          endif
       endif
    enddo
666 continue
    currentlyDoingHydroStep = tmp
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
    character(len=10) :: boundary
    logical :: tmp

    tmp = currentlyDoingHydroStep 
    if (myrankGlobal == 0) goto 666
    currentlyDoingHydroStep = .true.

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call  setupEdgesLevel(child, grid, nDepth)
       end do
       else
       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle


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
    currentlyDoingHydroStep = tmp
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
    logical, save :: firsttime=.true.
    logical :: tmp

    tmp = currentlyDoingHydroStep
    if (myrankGlobal ==0 ) goto 666
    currentlyDoingHydroStep = .true.

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
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          thisOctal%ghostCell(subcell) = .false.
          if (cylindricalHydro.or.spherical) then
             rVec = subcellCentre(thisOctal,subcell)
             if (rvec%x < 0.d0) then
                thisOctal%ghostCell(subcell) = .true.
                thisOctal%boundaryCondition(subcell) = 1
!                thisOctal%boundaryCondition(subcell) = 4

! now look at next cell to the right

                thisOctal%boundaryPartner(subcell) = VECTOR(1.d0, 0.d0, 0.d0)
                call locatorToNeighbour(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 1, locator)
                thisOctal%boundaryPartner(subcell) = locator
                tempOctal => thisOctal
                tempSubcell = 1
                call findSubcellLocal(locator, tempOctal, tempSubcell)
                if (tempOctal%centre%x > 0.d0) cycle

! ok we get to here the cell is further to the left, check whether it only has one cell between it and the real grid


                thisOctal%boundaryPartner(subcell) = VECTOR(1.d0, 0.d0, 0.d0)
                call locatorToNeighbour(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 3, locator)
                thisOctal%boundaryPartner(subcell) = locator
                tempOctal => thisOctal
                tempSubcell = 1
                call findSubcellLocal(locator, tempOctal, tempSubcell)
                if (tempOctal%centre%x > 0.d0) cycle

! ok now we now that the ghostcell is deeply embedded - just set the boundary partner to be the ghostcell to its right


                thisOctal%boundaryPartner(subcell) = VECTOR(1.d0, 0.d0, 0.d0)
                call locatorToNeighbour(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 1, locator)
                thisOctal%boundaryPartner(subcell) = locator
                cycle

             endif
          endif

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
!                      print *, "found a ghost ", subcellCentre(thisOctal, subcell)
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
!                   print *, "EDGE IS GHOST ", subcellCentre(thisOctal, subcell)
                   exit
                endif
             enddo
             



          endif

          if (thisOctal%ghostCell(subcell)) then
             select case(thisOctal%boundaryCondition(subcell))
             case(1, 4, 5, 6, 9)
                
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

                case(7)

                   
                case DEFAULT
                   if(firsttime) then
                      write(*,*) "Unknown boundary condition in setupghostcells2 B: ",thisOctal%boundaryCondition(subcell)
                      firsttime = .false.
                   end if
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
    currentlyDoingHydroStep = tmp
  end subroutine setupGhosts

  subroutine createGhostCells(grid)
    type(GRIDTYPE) :: grid
    integer :: nProbes, iProbe
    type(VECTOR) :: probe(4)

    currentlyDoingHydroStep = .true.
    nProbes = 4
    probe(1) = VECTOR(1.d0, 0.d0, 0.d0)
    probe(2) = VECTOR(-1.d0, 0.d0, 0.d0)
    probe(3) = VECTOR(0.d0, 0.d0, 1.d0)
    probe(4) = VECTOR(0.d0, 0.d0, -1.d0)
    do iProbe = nProbes, 1, -1
       write(*,*) "calling setupghostsnew ",iprobe
       call setupGhostsNew(grid%octreeRoot, grid, probe(iProbe), dble(iProbe))
    enddo
    currentlyDoingHydroStep = .false.
  end subroutine createGhostCells

  recursive subroutine setupGhostsNew(thisOctal, grid, probe, flag)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal, octal1, octal2, octal3, octal4
    type(octal), pointer  :: child 
    integer :: subcell
    real(double) :: flag
    type(VECTOR) :: locator, rVec
    integer :: i, subcell1, subcell2, subcell3, subcell4
    type(VECTOR) :: probe

    if (myrankGlobal ==0 ) goto 666
    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupGhostsNew(child,  grid, probe, flag)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%edgeCell(subcell)) then
             thisOctal%ghostCell(subcell) = .true.
             rVec = subcellCentre(thisOctal, subcell)
             locator = rVec + (thisOctal%subcellSize/2.d0+grid%halfSmallestSubcell * 0.1d0) * probe
             if (inOctal(grid%octreeRoot, locator)) then
                octal1 => thisOctal
                subcell1 = subcell


                octal2 => octal1
                call findSubcellLocal(locator, octal2, subcell2)
                rVec = subcellCentre(octal2, subcell2)
                locator = rVec + (octal2%subcellSize/2.d0+grid%halfSmallestSubcell * 0.1d0) * probe

                if (.not.octalOnThread(octal2, subcell2, myRankGlobal)) cycle
                if (octal2%edgeCell(subcell2)) cycle
                if (octal2%ghostCell(subcell2)) cycle

                octal3 => octal2
                call findSubcellLocal(locator, octal3, subcell3)
                rVec = subcellCentre(octal3, subcell3)
                locator = rVec + (octal3%subcellSize/2.d0+grid%halfSmallestSubcell * 0.1d0) * probe

                octal4 => octal3
                call findSubcellLocal(locator, octal4, subcell4)

                octal1%ghostCell(subcell1) = .true.

                if (.not.octal4%ghostCell(subcell4)) then
                   octal1%boundaryPartner(subcell1) = subcellCentre(octal4, subcell4)
                   octal4%feederCell(subcell4) = .true.
!                   octal4%rho(subcell4) = -flag
                endif

                octal2%ghostCell(subcell2) = .true.

                if (.not.octal3%ghostCell(subcell3)) then
                   octal2%boundaryPartner(subcell2) = subcellCentre(octal3, subcell3)
                   octal3%feederCell(subcell3) = .true.
!                   octal3%rho(subcell3) = -flag
                endif

!                octal1%rho(subcell1) = flag
!                octal2%rho(subcell2) = flag

!                write(*,*) "flags set to ",flag

             endif
          endif
       endif
    enddo
666 continue
  end subroutine setupGhostsNew

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
    logical :: trigger
    logical, save :: firsttime=.true.

    if (myrankGlobal ==0 ) goto 666
    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call  setupGhostsLevel(child, grid, nDepth)
       end do
       else
       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle


!          if (cylindricalHydro) then
!             rVec = subcellCentre(thisOctal,subcell)
!             if (rvec%x < 0.d0) then
!                thisOctal%ghostCell(subcell) = .true.
!                thisOctal%boundaryCondition(subcell) = 1
!
!! now look at next cell to the right
!
 !               thisOctal%boundaryPartner(subcell) = VECTOR(1.d0, 0.d0, 0.d0)
 !               call locatorToNeighbourLevel(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 1, locator, nDepth)
 !               thisOctal%boundaryPartner(subcell) = locator
 !               tempOctal => grid%octreeRoot
 !               tempSubcell = 1
 !               call findSubcellLocalLevel(locator, tempOctal, tempSubcell, nDepth)
 !               if (tempOctal%centre%x > 0.d0) cycle
!
!! ok we get to here the cell is further to the left, check whether it only has one cell between it and the real grid
!
!
!                thisOctal%boundaryPartner(subcell) = VECTOR(1.d0, 0.d0, 0.d0)
!                call locatorToNeighbourLevel(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 3, locator, nDepth)
!                thisOctal%boundaryPartner(subcell) = locator
!                tempOctal => grid%octreeRoot
!                tempSubcell = 1
!                call findSubcellLocalLevel(locator, tempOctal, tempSubcell, nDepth)
!                if (tempOctal%centre%x > 0.d0) cycle
!
!! ok now we now that the ghostcell is deeply embedded - just set the boundary partner to be the ghostcell to its right
!
!
!                thisOctal%boundaryPartner(subcell) = VECTOR(1.d0, 0.d0, 0.d0)
!                call locatorToNeighbourLevel(grid, thisOctal, subcell, thisOctal%boundaryPartner(subcell), 1, locator, nDepth)
!                thisOctal%boundaryPartner(subcell) = locator
!                cycle
!
!             endif
!          endif



          trigger=.false.
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
                      trigger = .true.
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
                   trigger = .true.
                   exit
                endif
             enddo
          endif

          if (thisOctal%ghostCell(subcell) .and. .not. thisOctal%corner(subcell) .and. trigger) then
             select case(thisOctal%boundaryCondition(subcell))
             case(1, 4, 5, 6, 9)
                
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
                   if(firsttime) then
                      write(*,*) "Unknown boundary condition in setupghostcells2 C: ",thisOctal%boundaryCondition(subcell)
                      firsttime=.false.
                   end if
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

   integer function getEndLoop(evenuparray) result(max)
     integer :: evenUpArray(nHydroThreadsGlobal)
     integer :: i

     max = 0
     do i = 1, nhydrothreadsglobal
        if(evenuparray(i) >  max) then
           max = evenuparray(i)
        end if
     end do

   end function getEndLoop

   integer function getNWorking(evenuparray, k) result(n)
     integer :: evenUpArray(nHydroThreadsGlobal)
     integer :: i, k

     n = 0
     do i = 1, nhydrothreadsglobal
        if(evenuparray(i) == k) then
           n = n + 1
        end if
     end do

   end function getNWorking


  subroutine refineGridGeneric(grid, tol, evenupArray, refinedSomeCells)
    use inputs_mod, only : minDepthAMR, maxDepthAMR
    use mpi
    type(GRIDTYPE) :: grid
    logical, optional :: refinedSomeCells
    logical :: globalConverged(512), tConverged(512)
    integer :: ierr, k, j
    real(double) :: tol
    real(double) :: index = 1.d0
    integer :: evenuparray(nHydroThreadsGlobal), itemp
    integer :: endloop, nworking
    integer, allocatable ::  safe(:, :)
!    character(len=80) :: plotfile



    globalConverged = .false.
    if (PRESENT(refinedSomeCells)) refinedSomeCells = .false.
    if (myrankGlobal == 0) goto 666
    if (minDepthAMR == maxDepthAMR) goto 666 ! fixed grid


    endloop = getEndLoop(evenuparray)

    allocate(safe(endloop, nHydroThreadsGlobal))

    safe = 0
    do k = 1, endloop
       do j = 1, nHydroTHreadsGlobal
          if(evenupArray(j) == k) then
             safe(k, j) = j 
          end if
       end do
    end do

    itemp = 0
    do
       globalConverged(myRankGlobal) = .true.
!       do iThread = 1, nHydroThreadsGlobal
       do k = 1, endloop
          if(evenuparray(myRankGlobal) /= k) then
!nworking is the number of threads that will be working at once in this stage of the even up
             nworking = getNWorking(evenUpArray, k)
             !         if (myrankGlobal /= iThread) then 
             !call hydroValuesServer(grid, iThread)
!             call hydroValuesServer(grid, nworking)
             call hydroValuesServer2(grid, nworking, k, endloop, safe)
          else
             call setAllUnchanged(grid%octreeRoot)
             call refineGridGeneric2(grid%octreeRoot, grid, globalConverged(myRankGlobal), tol, index, inheritval=.false.)
             call shutdownServers2(safe, k, endloop)
          endif
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       enddo
       index = index + 1.d0
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreadsGlobal, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
       itemp = itemp+1

       if (ALL(tConverged(1:nHydrothreadsGlobal))) then
          exit
       else
          if (PRESENT(refinedSomeCells)) refinedSomeCells = .true.
       endif
    enddo


    deallocate(safe)

    666 continue
  end subroutine refineGridGeneric


  recursive subroutine refineGridGeneric2(thisOctal, grid, converged, limit, index, inheritval)
    use inputs_mod, only : maxDepthAMR, photoionization, refineOnMass, refineOnTemperature, refineOnJeans
    use inputs_mod, only : refineonionization, massTol, refineonrhoe, amrtemperaturetol, amrrhoetol
    use inputs_mod, only : amrspeedtol, amrionfractol,  captureshocks, dorefine
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    logical, optional :: inheritval
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
    type(VECTOR) :: dirVec(6), centre, locator
    integer :: iSource
    logical :: split
    integer :: neighbourSubcell, nDir
    real(double) :: r, grad, maxGradient
    real(double) :: limit
    logical :: refineOnGradient
    real(double) :: rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z, pressure
    real(double) :: speed1, speed2, cs, bigJ, rhoJeans
    integer :: nd
    integer :: index1, index2, step
    real(double) :: index


    converged = .true.
    converged_tmp=.true.
    if (.not.doRefine) goto 666
    
    refineOnGradient = .not.photoionization      


    if(photoionization .and. hydrodynamics) then 
       refineongradient = .true.
    end if

    if (amrTolerance > 1.d29) refineOnGradient = .false.

   if(refineonionization .and. .not. photoionization) then
      refineOnIonization = .false.
   end if




    if(mod(index, 2.d0) /= 0.d0) then
       index1 = 1
       index2 = thisOctal%maxChildren
       step = 1
    else
       index2 = 1
       index1 = thisOctal%maxChildren
       step = -1
    end if

    do subcell = index1, index2, step

!    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call refineGridGeneric2(child, grid, converged, limit, index, inheritval)
                if (.not.converged) return
                if (.not.converged) converged_tmp = converged
                exit
             end if
          end do
          if (.not.converged_tmp) converged=converged_tmp
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%ghostCell(subcell)) cycle

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
                call getHydroValues(grid, locator, nd, rho, rhoe, rhou, rhov, rhow, energy, phi,x,y,z, pressure, .true.)

                split = .false.

                grad = abs((thisOctal%rho(subcell)-rho) / &
                     thisOctal%rho(subcell))
                maxGradient = max(grad, maxGradient)
                if (grad > limit) then
                   split = .true.
!                   write(*,*) "refined on density gradient ",grad,limit
                endif

                if(thisOctal%corner(subcell) .and. thisOCtal%nDepth < maxDepthAMR) split = .true.

                if(refineonspeed) then
                   if(thisOctal%threeD) then
                      speed1 = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 &
                      + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)
                      speed2 = sqrt(rhou**2 + rhov**2 + rhow**2)/rho
                   else if(thisOctal%twoD) then
                     speed1 = sqrt(thisOctal%rhou(subcell)**2  &
                      + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)
                      speed2 = sqrt(rhou**2 + rhow**2)/rho
                   else
                      speed1 = sqrt(thisOctal%rhou(subcell)**2)/thisOctal%rho(subcell)
                   end if

                   if (Speed1 > 1.d-10 .and. speed2 > 1.d-10) then
                      grad = abs(speed1-speed2) / speed1
                      maxGradient = max(grad, maxGradient)
                      if (grad > amrspeedtol) then
                         split = .true.
                      endif
                   endif
                end if

                if(captureshocks) then !.and. .not. thisOctal%changed(subcell)) then
                   if(refineShock(thisOctal, subcell, grid)) then
                      split = .true.
                   end if
                end if

                if (split) then

 
                   if ((thisOctal%nDepth < maxDepthAMR).and.(.not.thisOctal%changed(subcell))) then
                      call addNewChildWithInterp(thisOctal, subcell, grid)
                      converged = .false.
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
   
                   if(thisOctal%corner(subcell) .and. (thisOctal%nDepth < maxDepthAMR).and.(.not.thisOctal%changed(subcell))) then
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

    if (converged.and.(globalnSource>0)) then
       centre = subcellCentre(thisOctal, subcell)
       do iSource = 1, globalNSource
          
          r = modulus(centre - globalsourceArray(iSource)%position)/smallestCellSize
          if ((r < 5.d0) .and. (thisOctal%nDepth < maxDepthAMR)) then
             write(*,*) myrankGlobal, " splitting cell near sink"
             call addNewChildWithInterp(thisOctal, subcell, grid)
             converged = .false.
             exit
          endif
       enddo
       if (.not.converged) exit
    endif

    if (converged.and.(globalnSource>0)) then
       do iSource = 1, globalNSource
          if (inSubcell(thisOctal,subcell, globalSourceArray(isource)%position) &
               .and. (thisOctal%nDepth < maxDepthAMR)) then
             call addNewChildWithInterp(thisOctal, subcell, grid)
             converged = .false.
             exit
          endif
       enddo
       if (.not.converged) exit
    endif
             

    if (converged.and.refineOnMass) then
       if (((thisOctal%rho(subcell)*1.d30*thisOctal%subcellSize**3) > massTol) &
            .and.(thisOctal%nDepth < maxDepthAMR).and.(.not.thisOctal%changed(subcell)))  then
          call addNewChildWithInterp(thisOctal, subcell, grid)
          converged = .false.
!        print *, "split D ", thisOctal%nDepth
!        write(*,*) masstol, (thisOctal%rho(subcell)*1.d30*thisOctal%subcellSize**3)
          exit
       endif
    endif

    if (converged.and.refineOnJeans) then
       bigJ = 0.25d0
       cs = soundSpeed(thisOctal, subcell)
       rhoJeans = max(1.d-30,bigJ**2 * pi * cs**2 / (bigG * returnCodeUnitLength(thisOctal%subcellSize*1.d10)**2)) ! krumholz eq 6
       massTol = rhoJeans*1.d30*cellVolume(thisOctal,subcell)
       if (((thisOctal%rho(subcell)*1.d30*cellVolume(thisOctal,subcell)) > massTol) &
            .and.(thisOctal%nDepth < maxDepthAMR).and.(.not.thisOctal%changed(subcell)))  then
          write(*,*)  myrankGlobal," splitting on mass: ",thisOctal%rho(subcell)*1.d30*thisOCtal%subcellSize**3 / masstol
!          write(*,*) "mass tol ",masstol

!          write(*,*) myrankglobal, " calling interp after split on jeans"
          call addNewChildWithInterp(thisOctal, subcell, grid)
          converged = .false.
!        print *, "split E ", thisOctal%nDepth
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

             if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then    

                if((maxGradient > amrtemperaturetol) .and. (thisOctal%nDepth < maxDepthAMR)) then 
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   converged = .false.
                   exit
                endif
             endif
          end if
       enddo
      endif


    if(converged .and. refineOnRhoe) then 
       do i = 1, nDir
          maxGradient = 1.d-30
          locator = subcellCentre(thisOctal, subcell) + &
          (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
          if (inOctal(grid%octreeRoot, locator)) then

             neighbourOctal => thisOctal

             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)

!Thaw - modded to use specifiable limit 

             grad = abs((thisOctal%rhoe(subcell) - neighbourOctal%rhoe(neighbourSubcell)) / & 
             (thisOctal%rhoe(subcell)))

             maxGradient = max(grad, maxGradient)

             if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then    
                if((maxGradient > amrRhoeTol) .and. (thisOctal%nDepth < maxDepthAMR)) then                   
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   converged = .false.
                   exit
                endif
             endif
          end if
       enddo
    end if

!    if(converged .and. refineOnSpeed) then 
!       do i = 1, nDir
!          maxGradient = 1.d-30
!          locator = subcellCentre(thisOctal, subcell) + &
!          (thisOctal%subcellSize/2.d0+0.01d0*grid%halfSmallestSubcell) * dirVec(i)
!          if (inOctal(grid%octreeRoot, locator)) then
!
!             neighbourOctal => thisOctal
!
!             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
!
!!Thaw - modded to use specifiable limit 
!
!
!             speed1 = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)
!             speed2 = sqrt(neighbourOctal%rhou(neighboursubcell)**2 + neighbourOctal%rhov(neighboursubcell)**2 + & 
!                neighbourOctal%rhow(neighboursubcell)**2)
!
!             grad = abs((speed1 - speed2) / speed1)
!            
!             maxGradient = max(grad, maxGradient)
!
!             if (octalOnThread(neighbourOctal, neighbourSubcell, myRank)) then    
!                if((maxGradient > amrspeedtol) .and. (thisOctal%nDepth < maxDepthAMR)) then                                        
!                   call addNewChildWithInterp(thisOctal, subcell, grid)
!                   converged = .false.
!                   exit
!                endif
!             endif
!          end if
!       enddo
 !   end if

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

             if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
!                if ((thisOctal%ionFrac(subcell,1) > 0.9d0).and.(neighbourOctal%ionFrac(neighbourSubcell,1) < 0.1d0) .and. &
                if((maxGradient > amrIonFracTol) .and. (thisOctal%nDepth < maxDepthAMR)) then
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
666 continue
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
    converged = .true.
    converged_tmp=.true.


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
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle


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
    converged = .true.
    converged_tmp=.true.

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
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle


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
    use inputs_mod, only : smallestCellSize
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, tempOctal
    integer :: subcell, tempSubcell, i, nCells
    type(VECTOR) :: direction, locator, rVec
    tempOctal => thisOctal
    tempSubcell = subcell
    do i = 1, nCells
       rVec = subcellCentre(tempOctal, tempSubcell) + &
            (tempOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize) * direction 
       tempOctal => grid%octreeRoot
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
    type(octal), pointer  :: child !, neighbourOctal
    integer :: subcell, i !, neighbourSubcell
    logical :: unrefine
    integer :: nc
    integer, intent(inout) :: nUnrefine
    real(double) :: rhow(8), rhov(8), rhou(8), rho(8), rhoe(8) , fac, limit
    real(double) :: cs(8), mass
    logical :: refinedLastTime, ghostCell !, split
    real(double) :: r !, maxGradient
    real(double) ::  meancs !,rhou1, rhov1, rhow1, speed, thisSpeed, rho1, rhoe1,, meanrho, grad
    real(double) :: splitLimit, dv
    integer :: isource !, ndir
    logical :: debug, cornerCell
!    real(double) :: x, xnext, qnext, qViscosity(3,3) , q , flux, phi, phigas, pressure, px, py, pz
!    real(double) rhot, rhoet, rhout, rhovt, rhowt, rm1, um1, pm1
!    integer :: nd
    type(VECTOR) ::  centre !, dirvec(6), locator

    debug = .false.
!    limit  = 5.0d-3
    limit = splitlimit

    unrefine = .true.
    refinedLastTime = .false.
    ghostCell = .false.
    cornerCell = .false.
    nc = 0
    mass = 0.d0


    
    IF ( thisOctal%nChildren > 0 ) THEN
       ! call this subroutine recursively on each of its children
       i = 1
       DO while(i <= thisOctal%nChildren)
          child => thisOctal%child(i)
          CALL unrefineCells(child, grid, nUnrefine, splitLimit)
          i = i + 1
       END DO
       goto 666
    END IF

    if (.not.octalonThread(thisOctal, 1, myRankGlobal)) goto 666

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%refinedLastTime(subcell)) cycle

          nc = nc + 1
          rho(nc) = max(thisOctal%rho(subcell),1.d-30)
          rhoe(nc) = max(thisOctal%rhoe(subcell),1.d-30)
          rhou(nc) = max(thisOctal%rhou(subcell),1.d-30)
          rhov(nc) = max(thisOctal%rhov(subcell),1.d-30)
          rhow(nc) = max(thisOctal%rhow(subcell),1.d-30)
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


          if (globalnSource > 0) then
             centre = subcellCentre(thisOctal, subcell)
             do iSource = 1, globalNSource
                r = modulus(centre - globalsourceArray(iSource)%position)
                if (r*1.d10 < globalSourceArray(iSource)%accretionRadius) then
                   unrefine = .false.
                   goto 666
                endif
             enddo
          endif



       enddo


    unrefine = .false.

    if(cornercell) unrefine=.false.

    if ((nc > 1).and..not.cornerCell .and. .not. ghostcell) then

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
       if (fac > limit) then !.and.(rhou(1)/meancs > 1.d-2)) then
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





!    if (unrefine) then
!
!       if (thisOctal%threed) then
!          nDir = 6
!          dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
!          dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
!          dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
!          dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
!          dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
!          dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
!       else if (thisOctal%twod) then
!          nDir = 4
!          dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
!          dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
!          dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
!          dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
!       else
!          nDir = 2
!          dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
!          dirVecx(2) = VECTOR(-1.d0, 0.d0, 0.d0)
!       endif
!       
!       do i = 1, nDir
!          do subcell = 1, thisOctal%maxChildren
!             locator = subcellCentre(thisOctal, subcell) + (thisOctal%subcellSize/2.d0 + 0.01d0*smallestCellSize)*dirVec(i)
!             neighbourOctal => thisOctal
!             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
!             call getNeighbourValues(grid, thisOctal, Subcell, neighbourOctal, neighbourSubcell, dirVec(i), q, rhot, rhoet, &
!                  rhout, rhovt, rhowt, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz,  rm1,um1, pm1, qViscosity)
!             if (thisOctal%nDepth < nd) then
!                unrefine = .false.
!                goto 666
!             endif
!          enddo
!       enddo
!
!
!       meanRho = sum(rho(1:nc))/dble(nc)
!       thisSpeed = (sum(rhou(1:nc))/dble(nc))**2 + (sum(rhov(1:nc))/dble(nc))**2 + (sum(rhow(1:nc))/dble(nc))**2
!       thisSpeed = sqrt(thisSpeed/meanrho**2)
!       r = thisOctal%subcellSize + 0.01d0*grid%halfSmallestSubcell
!       centre = thisOctal%centre
!
!!       unrefine = .true.
!       do i = 1, nDir
!          maxGradient = 1.d-30
!          locator = centre + r*dirVec(i)
!          if (inOctal(grid%octreeRoot, locator)) then
!             neighbourOctal => thisOctal
!             call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
!             
!             if (octalOnThread(neighbourOctal, neighbourSubcell, myrankGlobal)) then
!
!                rho1 = neighbourOctal%rho(neighbourSubcell)
!                rhou1 = neighbourOctal%rhou(neighbourSubcell)
!                rhov1 = neighbourOctal%rhov(neighbourSubcell)
!                rhow1 = neighbourOctal%rhow(neighbourSubcell)
!                rhoe1 = neighbourOctal%rhoe(neighbourSubcell)
!                split = .false.
!                
!                grad = abs((meanrho-rho1) / meanRho)
!                if (grad > splitLimit) then
!                   split = .true.
!                endif
!                speed = sqrt(rhou1**2 + rhov1**2 + rhow1**2)/rho1
!                if (thisSpeed > 0.d0) then
!                   grad = abs(thisSpeed-speed) / thisSpeed
!                   if ((grad > limit).and.(speed/meancs > 1d-2)) then
!                      split = .true.
!                   endif
!                endif
!
!                !THaw - added ".and. unrefine because previously it would only care about the last direction
!                if (split) then
!                   unrefine = .false.
!                endif
!             endif
!          endif
!       enddo
!    endif
!
!
!    if (unrefine.and.photoionization) then
!       if ( (.not.all(thisOctal%ionFrac(1:thisOctal%maxChildren,2) < 0.01d0)).or. &
!            (.not.all(thisOctal%ionFrac(1:thisOctal%maxChildren,2) < 0.99d0)) ) unrefine = .false.
!    endif
!
!    if(cornercell) unrefine=.false.
!
    if ((thisOctal%nChildren == 0).and.unrefine) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
       nunrefine = nunrefine + 1
    endif
666 continue    
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


   subroutine evenUpGridMPI(grid, inheritFlag, evenAcrossThreads, evenUpArray, dumpfiles)
    use mpi  
    use inputs_mod, only : minDepthAMR, maxDepthAMR

    type(GRIDTYPE) :: grid
    logical :: inheritFlag
    integer, optional :: dumpfiles
    integer :: nThreads, thisThread
    integer :: ierr
    logical :: globalConverged(512), localChanged(512), tConverged(512)
    type(VECTOR) :: locs(200000), eLocs(200000)
    integer :: nLocs(512), tempNlocs(20000)
    integer :: thread(100000), nLocsGlobal,i, depth(200000)
    integer, intent(in) :: evenUpArray(:)
    integer :: nSent, nTemp
    real(double) :: temp1d(80000), tempsent1d(80000)
!    real(double) :: temp(20000,4),tempsent(4)
    integer :: eDepth(200000)
    integer :: iThread, nExternalLocs
    integer :: iter, k, endloop, j, nworking
    integer, parameter :: tag = 1
    logical :: globalChanged(512)
    integer :: status(MPI_STATUS_SIZE)
    logical :: evenAcrossThreads, converged
    logical :: allThreadsConverged
    character(len=30) :: vtkFilename
    real(double) :: index = 1.d0
    integer, allocatable ::  safe(:, :)

!    Character(len=20) :: plotfile
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

    if (checkEven(grid)) then
       call writeInfo("No evenup necessary")
       goto 666
    endif


    allThreadsConverged = .false.

    endloop = getEndLoop(evenuparray)

    allocate(safe(endloop, nHydroThreadsGlobal))
    safe = 0
    do k = 1, endloop
       do j = 1, nHydroTHreadsGlobal
          if(evenupArray(j) == k) then
             safe(k, j) = j 
          end if
       end do
    end do


    do while(.not.allThreadsConverged)
       localChanged = .false.

       globalConverged = .false.
       do 
          globalConverged(myRankGlobal) = .true.
!          do iThread = 1, nHydroThreadsGlobal
          do k = 1, endloop
!             write(*,*) "k loop ",k,endloop
!             if (myrankGlobal /= iThread) then
             if(evenUpArray(myRankGlobal) /= k) then
                nworking = getNWorking(evenUpArray, k)
!                call hydroValuesServer(grid, iThread)
!                call hydroValuesServer(grid, nworking)
                call hydroValuesServer2(grid, nworking, k, endloop, safe)
             else
                call evenUpGrid(grid%octreeRoot, grid,  globalConverged(myrankGlobal), index, inherit=inheritFlag)

                call unsetGhosts(grid%octreeRoot)
                call setupEdges(grid%octreeRoot, grid)
                call setupGhosts(grid%octreeRoot, grid)!, flag=.true.)
!                call shutdownServers()
                call shutdownServers2(safe, k, endloop)
             endif
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!          enddo
             iter = iter + 1
             if (PRESENT(dumpfiles)) then
                write(vtkfilename,'(a,i4.4,a,i4.4,a)') "start",dumpfiles,"_",iter,".vtk"
                call writeVtkFile(grid, vtkFilename)
             endif


          end do

          index = index + 1.d0
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreadsGlobal, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
!          if (myrank == 1) write(*,*) "first evenup ",tConverged(1:nHydroThreadsGlobal)
          if (ALL(tConverged(1:nHydrothreadsGlobal))) exit
       enddo

       if (evenAcrossThreads) then             
             nLocs = 0
             nExternalLocs = 0
             call locatorsToExternalCells(grid%octreeRoot, grid, nLocs(myRankGlobal), locs, thread, depth)
             call MPI_ALLREDUCE(nLocs, tempnLocs, nHydroThreadsglobal, MPI_INTEGER, MPI_SUM,amrCOMMUNICATOR, ierr)
             
             nLocsGlobal = SUM(tempnLocs)
             do thisThread = 1, nHydroThreadsGlobal
                if(thisThread == myRankGlobal) then
                   do iThread = 1, nHydroThreadsGlobal
                      nTemp = 0
                      if (iThread /= myRankGlobal) then
                         do i = 1 , nLocs(myRankGlobal)
                            if (thread(i) == iThread) then
                               nTemp = nTemp + 1

!tjh
!                               temp(nTemp,1) = locs(i)%x
!                               temp(nTemp,2) = locs(i)%y
!                               temp(nTemp,3) = locs(i)%z
!                               temp(nTemp,4) = dble(depth(i))

                               temp1d(ntemp*4-3) = locs(i)%x
                               temp1d(ntemp*4-2) = locs(i)%y
                               temp1d(ntemp*4-1) = locs(i)%z
                               temp1d(ntemp*4  ) = dble(depth(i))

                            endif
                         enddo
                         call mpi_send(nTemp, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
                         if (nTemp > 0) then
                            call mpi_send(temp1d, ntemp*4, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
!                            do i = 1, nTemp(1)
!                               call mpi_send(temp(i,1:4), 4, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
!                            enddo
                         endif
                         
                      endif
                   enddo
                else

                   call mpi_recv(nSent, 1, MPI_INTEGER, thisThread, tag, localWorldCommunicator, status, ierr)
                         ! write(*,*) myRank, "received ",nSent, "from ", thisThread
                  
                   if (nSent > 0) then
                         call mpi_recv(tempsent1d, nsent*4, MPI_DOUBLE_PRECISION, thisThread, tag, localWorldCommunicator, &
                              status, ierr)

                      do i = 1, nSent

                         eLocs(nExternalLocs + i)%x = tempSent1d(i*4-3)
                         eLocs(nExternalLocs + i)%y = tempSent1d(i*4-2)
                         eLocs(nExternalLocs + i)%z = tempSent1d(i*4-1)
                         eDepth(nExternalLocs + i) = nint(tempSent1d(i*4   ))

!                         call mpi_recv(tempsent, 4, MPI_DOUBLE_PRECISION, thisThread, tag, localWorldCommunicator, status, ierr)
!                         eLocs(nExternalLocs + i)%x = tempsent(1)
!                         eLocs(nExternalLocs + i)%y = tempsent(2)
!                         eLocs(nExternalLocs + i)%z = tempsent(3)
!                         eDepth(nExternalLocs + i) = nint(tempsent(4))
                         
                      enddo
                   endif
                   nExternalLocs = nExternalLocs + nSent


               endif
             enddo
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!             do iThread = 1, nThreadsGlobal-1
!                if (myrankGlobal /= iThread) then
             do k = 1, endloop
                if(evenUpArray(myRankGlobal) /= k) then
!                   call hydroValuesServer(grid, nworking)
                   nworking = getNWorking(evenUpArray, k)
                   call hydroValuesServer2(grid, nworking, k, endloop, safe)
                else
!                   if(myRank == 2) then
!                      print *, "RANK TIME ", myRank, nExternalLocs, iThread
!                   end if
                   do i = 1, nExternalLocs
                      call splitAtLocator(grid, elocs(i), edepth(i), localChanged(myRankGlobal))
                      converged=.true.
                      if(ANY(localChanged)) converged = .false.
                   enddo
!                 call shutdownServers()
             call shutdownServers2(safe, k, endloop)
              endif
               call MPI_BARRIER(amrCOMMUNICATOR, ierr)
            enddo
!           end do
          
          if (PRESENT(dumpfiles)) then
             write(vtkfilename,'(a,i4.4,a,i4.4,a)') "aftersplit",dumpfiles,"_",iter,".vtk"
             call writeVtkFile(grid, vtkFilename)
          endif

          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          globalConverged = .false.
          do
             globalConverged(myRankGlobal) = .true.
!             do iThread = 1, nThreadsGlobal-1
             do k = 1, endloop
!                if (myrankGlobal /= iThread) then
                if(evenUpArray(myRankGlobal) /= k) then
!                   call hydroValuesServer(grid, iThread)
!                   call hydroValuesServer(grid, nworking)
                   nworking = getNWorking(evenUpArray, k)
                   call hydroValuesServer2(grid, nworking, k, endloop, safe)
                else
                   call evenUpGrid(grid%octreeRoot, grid,  globalConverged(myrankGlobal), index, inherit=inheritFlag)

                   call unsetGhosts(grid%octreeRoot)
                   call setupEdges(grid%octreeRoot, grid)
                   call setupGhosts(grid%octreeRoot, grid)!, flag=.true.)
 !                  call shutdownServers()
                   call shutdownServers2(safe, k, endloop)
                endif
                call MPI_BARRIER(amrCOMMUNICATOR, ierr)
             enddo
             index = index + 1.d0
             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
            call MPI_ALLREDUCE(globalConverged, tConverged, nHydroThreadsGlobal, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)
!             if (myrank == 1) write(*,*) "second evenup ",tConverged(1:nThreadsGlobal-1)
             if (ALL(tConverged(1:nHydrothreadsGlobal))) exit

          iter = iter + 1
          if (PRESENT(dumpfiles)) then
             write(vtkfilename,'(a,i4.4,a,i4.4,a)') "aftereven",dumpfiles,"_",iter,".vtk"
             call writeVtkFile(grid, vtkFilename)
          endif

          enddo

          call MPI_ALLREDUCE(localChanged, globalChanged, nHydroThreadsGlobal, MPI_LOGICAL, MPI_LOR, amrCOMMUNICATOR, ierr)

!          if (myrank == 1) write(*,*) "globalChanged ",globalChanged(1:nThreadsGlobal-1)
          if (ANY(globalChanged(1:nHydroThreadsGlobal)) .or. .not. converged) then
             allThreadsConverged = .false.
          else
             allThreadsConverged = .true.
          endif
      

       else
          allThreadsConverged = .true.
       endif

    end do

    deallocate(safe)

666 continue
  end subroutine evenUpGridMPI

  recursive subroutine evenUpGrid(thisOctal, grid,  converged, index, inherit)
    use inputs_mod, only : maxDepthAMR
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    type(octal), pointer :: bOctal
    integer :: subcell, i, bSubcell
    logical :: converged, converged_tmp
    type(VECTOR) :: dirVec(14), centre, octVec, rVec, locator, bVec, direction!, nVec
!    type(vector) :: vecStore(3)
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    logical, optional :: inherit
    real(double) :: rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z, pressure
    integer :: nd!, nd2
    real(double) :: index
    integer :: index1, index2, step

    converged = .true.
    converged_tmp=.true.



    if(mod(index, 2.d0) /= 0.d0) then
       index1 = 1
       index2 = thisOctal%maxChildren
       step = 1
    else
       index2 = 1
       index1 = thisOctal%maxChildren
       step = -1
    end if

    do subcell = index1, index2, step

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call evenUpGrid(child, grid,  converged, index, inherit)
                if (.not.converged) return
                if (.not.converged) converged_tmp = converged
                exit
             end if
          end do
          if (.not.converged_tmp) converged=converged_tmp
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell)


!          if (cylindricalHydro) then
!             if (thisOctal%edgeCell(subcell).and.(.not.inOctal(grid%octreeRoot, centre - VECTOR(r, 0.d0, 0.d0)))) then
!                if (thisOctal%nDepth < minDepthAMR) then 
!                   call addNewChildWithInterp(thisOctal, subcell, grid)
!                   converged = .false.
!                   exit
!                endif
!             endif
!             if (thisOctal%ghostcell(subcell)) then
!                if (thisOctal%nDepth < maxDepthAMR) then 
!                   call addNewChildWithInterp(thisOctal, subcell, grid)
!                   converged = .false.
!                   exit
!                endif
!             endif
!          endif

          !Thaw - force ghostcells to match their partner refinement
          if(thisOctal%ghostcell(subcell)) then
             locator = thisOctal%boundaryPartner(subcell)
             bOctal => thisOctal             
             rVec = subcellCentre(thisOctal, subcell)

             
             if (cylindricalHydro) then



                call findSubcellLocal(locator, bOctal, bSubcell)
                bVec = subcellCentre(bOctal, bSubcell)
                if(thisOctal%nDepth < bOctal%nDepth .and. thisOctal%nDepth < maxDepthAMR .and. &
                     (.not.bOCtal%ghostCell(bSubcell))) then
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   rVec = subcellCentre(thisOctal, subcell)
                   converged = .false.
                   exit
                end if

             else

                call findSubcellLocal(locator, bOctal, bSubcell)
                bVec = subcellCentre(bOctal, bSubcell)
                if(thisOctal%nDepth < bOctal%nDepth .and. thisOctal%nDepth < maxDepthAMR) then
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   rVec = subcellCentre(thisOctal, subcell)
                   converged = .false.
                   exit
                end if

             endif
             !THAW - need to provide a case for when boundary partner is not on thread (i.e. for periodics)
             if(octalOnThread(bOctal, bSubcell, myRankGlobal)) then
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

                if(octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then

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
                if(abs(direction%x) > abs(direction%z)) then !We are moving in the ±x direction
                   direction%z = 0.d0
                else                                         !Moving in the ±z direction       
                   direction%x = 0.d0
                   end if
             else
                direction%x = 0.d0
                if(abs(direction%y) > abs(direction%z)) then !We are moving in the ±y direction
                   direction%z = 0.d0
                   else                                         !Moving in the ±z direction    
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

                call getHydroValues(grid, locator, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z, pressure, .true.)
                if(nd > thisOctal%nDepth) then
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   converged = .false.
                   exit
                end if

             end if
                          
          end if
                       
          if (thisOctal%threed) then
             nDir = 14
             dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
             dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
             dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
             !corners
             dirVec(7) = VECTOR(-1.d0, 0.d0, 1.d0)
             dirVec(8) = VECTOR(1.d0, 0.d0, -1.d0)
             dirVec(9) = VECTOR(1.d0, 0.d0, 1.d0)
             dirVec(10) = VECTOR(-1.d0, 0.d0, -1.d0)
             dirVec(11) = VECTOR(1.d0, 1.d0, 1.d0)
             dirVec(12) = VECTOR(-1.d0, -1.d0, -1.d0)
             dirVec(13) = VECTOR(1.d0, -1.d0, 1.d0)
             dirVec(14) = VECTOR(-1.d0, 1.d0, -1.d0)

          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)  
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
             !corners

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
                
                call getHydroValues(grid, octVec, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z, pressure, .false.)
                
                !                   if (((neighbourOctal%nDepth-thisOctal%nDepth) > 1).and. (thisOCtal%ndepth < maxDepthAMR)) then
                if (((nd-thisOctal%nDepth) > 1).and. (thisOctal%nDepth < maxDepthAMR)) then
                   call addNewChildWithInterp(thisOctal, subcell, grid)
                   converged = .false.
                   exit
                endif
                
                if (octalOnThread(neighbourOctal, neighbourSubcell, myrankGlobal)) then         
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
    logical :: localChanged 


    call findSubcellTD(locator, grid%octreeRoot, thisOctal,subcell)
    
    if (myRankGlobal /= thisOctal%mpiThread(subcell)) then
       write(*,*) "There is a bug somewhere: ",myrankGlobal
       stop
    endif
!    WRITE(*,*) "ATTEMPTNG TO SPLIT", depth, thisOctal%nDepth,  myRank
    if (((depth-thisOctal%nDepth) > 1).and.(thisOctal%nDepth < maxDepthAMR)) then
       call addNewChildWithInterp(thisOctal,subcell,grid)
       localChanged = .true.
!       write(*,*) "split at locator"
!    else
!       WRITE(*,*) "ATTEMPTNG TO SPLIT FAILED", depth, thisOctal%nDepth,  myRank, maxDepthAMR
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
    integer, intent(out) :: nLocs
    converged = .true.
    converged_tmp=.true.



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
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

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

                   if (.not.octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
                      nLocs = nLocs + 1
!                      print *, "found ", nLocs, "on other threads", myRank
!                      print *, octVec
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
    integer :: nDependent
    logical :: alreadyInList
    integer :: iDep
    converged = .true.
    converged_tmp=.true.



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
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

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

                   if (neighbourOctal%mpiThread(neighboursubcell) /= myRankGlobal) then
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

    do iThread = 1, nHydroThreadsGlobal
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

    do iThread = 0, nHydroThreadsGlobal
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
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6)
    integer :: n, ndir, nd, nc
    real(double) :: x1, x2, g(6)
    real(double) :: deltaT, fracChange, gGrav, newPhi, frac !, d2phidx2(3), sumd2phidx2
    real(double) :: xnext, px, py, pz, rm1, um1, pm1
    gGrav = bigG  * lengthToCodeUnits**3 / (massToCodeUnits * timeToCodeUnits**2)


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
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                
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

  recursive subroutine gSweep2(thisOctal, grid, deltaT, fracChange,iter)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: iter
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6), probe(6)
    integer :: n, ndir
    real(double) :: x1, x2
    real(double) ::  g(6), dx, dxArray(6), g2(6), phiInterface(6)
    real(double) :: deltaT, fracChange, gGrav, newPhi, newerPhi, frac, d2phidx2(3), sumd2phidx2
    integer :: nd, nc
    real(double) :: tauMin
    real(double), parameter :: maxM = 100000.d0
    real(double) :: xnext, oldphi, px, py, pz, rm1, um1, pm1, thisR
    real(double), parameter :: SOR = 1.2d0
    gGrav = bigG * lengthToCodeUnits**3 / (massToCodeUnits * timeToCodeUnits**2)


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call gsweep2(child, grid, deltaT, fracChange, iter)
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
          if (.not.associated(thisOctal%adot)) then
             allocate(thisOctal%adot(1:thisOctal%maxChildren))
             thisOctal%adot = 0.d0
          endif

          if (.not.thisOctal%ghostCell(subcell) .and. .not. thisOctal%edgecell(subcell)) then

             locator = subcellCentre(thisOctal, subcell)
             thisR = returnCodeUnitLength(locator%x*gridDistanceScale)
             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + probe(n) * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
                neighbourOctal => thisOctal

                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                
                x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                

                if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
                   x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                   x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                   dx = x2 - x1
                else
                   if (nd == thisOctal%nDepth) then ! coarse/coarse or fine/fine
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1

                   else if (nd > thisOctal%nDepth) then ! coarse cells with a fine boundary
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                   else
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n) ! fine cells
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                   endif
                endif

                dxArray(n) = dx
                if (thisOctal%threed) then
                   g(n) =   (phigas - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                else
                   select case(n)
                   case(1, 2)
                      g(n) = gridDistanceScale * (x2 * phigas - x1 * thisOctal%phi_gas(subcell)) / &
                           (returnCodeUnitLength(dx*gridDistanceScale))
                      if (n == 1) thisOctal%adot(subcell) = thisOctal%adot(subcell) - phigas
                   case(3,4)
                      g(n) = (phigas - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                   end select
                endif
             enddo

!get the gravitational potential values at the cell interface
             do n = 1, nDir
                dx = thisOctal%subcellSize/2.d0
                if (thisOctal%threed) then
                   phiInterface(n) = thisOctal%phi_gas(subcell) +  g(n)*(&
                        returnCodeUnitLength(dx*gridDistanceScale))
                else
                   select case(n)
                   case(1, 2)
                   phiInterface(n) = thisOctal%phi_gas(subcell)*thisR  +  g(n)*(&
                        returnCodeUnitLength(dx*gridDistanceScale))
                   case(3,4)
                      phiInterface(n) = thisOctal%phi_gas(subcell) +  g(n)*(&
                           returnCodeUnitLength(dx*gridDistanceScale))
                   end select
                endif



             end do

!calculate the new gradient
             do n = 1, nDir
                dx = thisOctal%subcellSize/2.d0
                if (thisOctal%threed) then
                   g2(n) = (phiInterface(n) - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                else
                   select case(n)
                   case(1, 2)
                   g2(n) = (phiInterface(n) - thisR*thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                   case(3, 4)
                   g2(n) = (phiInterface(n) - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                   end select
                endif
             end do
!            tauMin  = (0.5d0/dble(2**maxDepthAMR))**2/4.d0
!            tauMax  = (0.5D0/dble(2**minDepthAMR))**2/4.d0
!            mVal = dble(iter)
!            epsilon = (1.d0/maxM) * (log(tauMax) - log(tauMin))
!            epsilon =  exp(epsilon)

!            deltaT = tauMin*(epsilon**mVal) * (gridDistanceScale * amrGridSize)**2
!            deltaT = tauMin * (gridDistanceScale * amrGridSize)**2

            tauMin  = (0.5d0/dble(2**maxDepthAMR))**2/4.d0
            deltaT = tauMin * (gridDistanceScale * amrGridSize)**2

!            if (Writeoutput) write(*,*) "deltaT ",deltaT, taumin,taumax,iter,epsilon
            if (thisOctal%twoD) then
               d2phidx2(1) = (1.d0/(thisR))*(g2(1) - g2(2)) / (returnCodeUnitLength(dx*gridDistanceScale))
               d2phidx2(2) = (g2(3) - g2(4)) / (returnCodeUnitLength(dx*gridDistanceScale))
               sumd2phidx2 = SUM(d2phidx2(1:2))
             else
                d2phidx2(1) = (g2(1) - g2(2)) / (returnCodeUnitLength(dx*gridDistanceScale))
                d2phidx2(2) = (g2(3) - g2(4)) / (returnCodeUnitLength(dx*gridDistanceScale))
                d2phidx2(3) = (g2(5) - g2(6)) / (returnCodeUnitLength(dx*gridDistanceScale))
                sumd2phidx2 = SUM(d2phidx2(1:3))
             endif

!             if(thisOctal%threed) then
!                deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*grid%halfSmallestSubcell))**2 / 6.d0
!             else
!                deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*grid%halfSmallestSubcell))**2 / 4.d0
!             endif

             oldPhi = thisOctal%phi_gas(subcell)
             newerPhi = thisOctal%phi_gas(subcell) + (deltaT * sumd2phidx2 - fourPi * gGrav * thisOctal%rho(subcell) * deltaT) 

             newPhi = (1.d0-SOR)*oldPhi + SOR*newerPhi

!             if (myrankglobal == 1) write(*,*) newphi,thisOctal%phi_gas(subcell),deltaT, sumd2phidx2, thisOctal%rho(subcell)
             if (.not.associated(thisOctal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))
             if (thisOctal%phi_gas(subcell) /= 0.d0) then
                frac = abs((newPhi - thisOctal%phi_gas(subcell))/thisOctal%phi_gas(subcell))
!                frac = abs(sumd2phidx2 - fourPi* gGrav * thisOctal%rho(subcell))/(fourPi * gGrav * thisOctal%rho(subcell))
                thisOctal%chiline(subcell) = frac

                fracChange = max(frac, fracChange)
             else
                fracChange = 1.d30
!                frac = fracChange
             endif
!             thisOctal%chiline(subcell) = frac

             thisOctal%phi_gas(subcell) = newPhi
             
          endif
       endif
    enddo
  end subroutine gSweep2

!This is for cylindrical hydro calcs?
  recursive subroutine gSweep2level(thisOctal, grid, fracChange, fracChange2, nDepth)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: nDepth
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6), probe(6)
    integer :: n, ndir
    real(double) :: x1, x2
    real(double) ::  g(6), dx, dxArray(6), g2(6), phiInterface(6)
    real(double) :: deltaT, fracChange, gGrav, newPhi, newerPhi, frac, d2phidx2(3), sumd2phidx2, fracChange2
    integer :: nd, nc
    real(double) :: tauMin, dfdrbyr
    real(double), parameter :: maxM = 100000.d0
    real(double) :: xnext, oldphi, px, py, pz, rm1, um1, pm1, thisR
    real(double), parameter :: SOR = 1.2d0
    gGrav = bigG * lengthToCodeUnits**3 / (massToCodeUnits * timeToCodeUnits**2)

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

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call gsweep2Level(child, grid, fracChange, fracChange2, ndepth)
       end do
       
       
    else

       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (.not.associated(thisOctal%adot)) then
             allocate(thisOctal%adot(1:thisOctal%maxChildren))
             thisOctal%adot = 0.d0
          endif

          if (.not.thisOctal%ghostCell(subcell) .and. .not. thisOctal%edgecell(subcell)) then

             locator = subcellCentre(thisOctal, subcell)
             thisR = returnCodeUnitLength(locator%x*gridDistanceScale)
             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + probe(n) * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
                neighbourOctal => grid%octreeRoot

                call findSubcellLocalLevel(locator, neighbourOctal, neighbourSubcell,nDepth)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

                
                x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                if (.not.thisOctal%twod) write(*,*) "twod error"

                if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
                   x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                   x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                   dx = x2 - x1
                else
                   if (nd == thisOctal%nDepth) then ! coarse/coarse or fine/fine
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1

                   else if (nd > thisOctal%nDepth) then ! coarse cells with a fine boundary
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                   else
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n) ! fine cells
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                   endif
                endif

                dxArray(n) = dx
                g(n) =   (phigas - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
             enddo

!get the gravitational potential values at the cell interface
             do n = 1, nDir
                dx = sign(thisOctal%subcellSize/2.d0,dxArray(n))
                phiInterface(n) = thisOctal%phi_gas(subcell) +  g(n)*(&
                     returnCodeUnitLength(dx*gridDistanceScale))
             end do

             if (thisOctal%twoD) then
                dfdrbyr = (phiInterface(1) - phiInterface(2))/(thisOctal%subcellSize*gridDistanceScale)
                dfdrbyr = dfdrbyr / thisR
             endif

!calculate the new gradient
             do n = 1, nDir
                dx = sign(thisOctal%subcellSize/2.d0,dxArray(n))
                g2(n) = (phiInterface(n) - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
             end do

             tauMin  = (0.5d0/dble(2**nDepth))**2/4.d0
             deltaT = tauMin * (gridDistanceScale * amrGridSize)**2

             dx = thisOctal%subcellSize
             if (thisOctal%twoD) then
               d2phidx2(1) = (g2(1) - g2(2)) / (returnCodeUnitLength(dx*gridDistanceScale))
               d2phidx2(2) = (g2(3) - g2(4)) / (returnCodeUnitLength(dx*gridDistanceScale))
               sumd2phidx2 = SUM(d2phidx2(1:2)) + dfdrbyr
             else
                d2phidx2(1) = (g2(1) - g2(2)) / (returnCodeUnitLength(dx*gridDistanceScale))
                d2phidx2(2) = (g2(3) - g2(4)) / (returnCodeUnitLength(dx*gridDistanceScale))
                d2phidx2(3) = (g2(5) - g2(6)) / (returnCodeUnitLength(dx*gridDistanceScale))
                sumd2phidx2 = SUM(d2phidx2(1:3))
             endif


             oldPhi = thisOctal%phi_gas(subcell)
             newerPhi = thisOctal%phi_gas(subcell) + (deltaT * sumd2phidx2 - thisOctal%source(subcell) * deltaT) 

             newPhi = (1.d0-SOR)*oldPhi + SOR*newerPhi

             if (.not.associated(thisOctal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))
             if (.not.associated(thisOctal%adot)) allocate(thisOctal%adot(1:thisOctal%maxChildren))
             if (thisOctal%phi_gas(subcell) /= 0.d0) then
                frac = abs((oldPhi - newPhi)/oldPhi)
                thisOctal%chiline(subcell) = frac
                thisOctal%adot(subcell) = frac
                fracChange = max(frac, fracChange)
                frac = abs(sumd2phidx2 - thisOctal%source(subcell))/thisOctal%source(subcell)
                fracChange2 = max(frac, fracChange2)

             else
                fracChange = 1.d30
                fracChange2 = 1.d30
             endif

             thisOctal%phi_gas(subcell) = newPhi
             
          endif
       enddo
    endif
  end subroutine gSweep2level

  recursive subroutine gaussSeidelSweep(thisOctal, grid, fracChange,nDepth, onlyCellsWithChildren)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: nDepth
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6)
    integer :: n, ndir
    real(double) ::  g(6), dx
    real(double) :: fracChange, sumd2phidx2
    integer :: nd, nc
    logical :: onlyCellsWithChildren
    real(double), parameter :: maxM = 100000.d0
    real(double) :: xnext, px, py, pz, rm1, um1, pm1
    real(double), parameter :: SOR = 1.2d0

    nDir = 6
    dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
    dir(2) = VECTOR(-1.d0, 0.d0, 0.d0)
    dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
    dir(4) = VECTOR(0.d0, 0.d0,-1.d0)
    dir(5) = VECTOR(0.d0, 1.d0, 0.d0)
    dir(6) = VECTOR(0.d0,-1.d0, 0.d0)
       
    if (thisOctal%nDepth == nDepth) then


       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (onlyCellsWithChildren .and. .not.thisOctal%hasChild(subcell)) cycle


          if (.not.associated(thisOctal%adot)) then
             allocate(thisOctal%adot(1:thisOctal%maxChildren))
             thisOctal%adot = 0.d0
          endif
          if (.not.associated(thisOctal%chiLine)) then
             allocate(thisOctal%chiline(1:thisOctal%maxChildren))
             thisOctal%chiLine = 0.d0
          endif


          if (.not.thisOctal%ghostCell(subcell)) then

             do n = 1, nDir                
                locator = subcellCentre(thisOctal, subcell) + dir(n) * (thisOctal%subcellSize/2.d0+0.01d0*smallestCellsize)
                neighbourOctal => grid%octreeRoot
                call findSubcellLocalLevel(locator, neighbourOctal, neighbourSubcell, nDepth)
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, dir(n), q, rho, rhoe, &
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                g(n) = phigas
             enddo
             dx = thisOCtal%subcellSize
             

             sumd2phidx2 = (sum(g(1:6)) - 6.d0*thisOctal%phi_gas(subcell))/(thisOctal%subcellSize*gridDistanceScale)**2
             thisOctal%chiline(subcell) = thisOctal%source(subcell) - sumd2phidx2 

             thisOctal%phi_gas(subcell) = 0.16666666666667d0*(SUM(g(1:6)) - &
                  thisOctal%source(subcell)*(returnCodeUnitLength(dx*gridDistanceScale))**2)

             if (thisOctal%source(subcell)/=0.d0) fracChange = max(fracChange,abs(sumd2phidx2 - &
                  thisOctal%source(subcell))/abs(thisOctal%source(subcell)))

          endif
       enddo
    else
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call gaussSeidelSweep(child, grid, fracChange, ndepth, onlyCellsWithChildren)
       end do
    endif

  end subroutine gaussSeidelSweep

  recursive subroutine calculateResiduals(thisOctal, grid, nDepth)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: nDepth
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6)
    integer :: n, ndir
    real(double) ::  g(6), dx
    real(double) :: sumd2phidx2
    integer :: nd, nc
    real(double), parameter :: maxM = 100000.d0
    real(double) :: xnext, px, py, pz, rm1, um1, pm1
    real(double), parameter :: SOR = 1.2d0

    nDir = 6
    dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
    dir(2) = VECTOR(-1.d0, 0.d0, 0.d0)
    dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
    dir(4) = VECTOR(0.d0, 0.d0,-1.d0)
    dir(5) = VECTOR(0.d0, 1.d0, 0.d0)
    dir(6) = VECTOR(0.d0,-1.d0, 0.d0)


    if (thisOctal%nDepth == nDepth) then


       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          thisOctal%chiline(subcell) = 0.d0


          if (.not.thisOctal%ghostCell(subcell)) then! .and. .not. thisOctal%edgecell(subcell)) then

             do n = 1, nDir                
                locator = subcellCentre(thisOctal, subcell) + dir(n) * (thisOctal%subcellSize/2.d0+0.01d0*smallestCellsize)
                neighbourOctal => grid%octreeRoot
                call findSubcellLocalLevel(locator, neighbourOctal, neighbourSubcell, nDepth)
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, dir(n), q, rho, rhoe, &
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                g(n) = phigas
             enddo
             dx = thisOCtal%subcellSize
             

             sumd2phidx2 = (sum(g(1:6)) - 6.d0*thisOctal%phi_gas(subcell))/(thisOctal%subcellSize*gridDistanceScale)**2
             thisOctal%chiline(subcell) = thisOctal%source(subcell) - sumd2phidx2 


          endif
       enddo

    else
       if (thisOctal%nChildren > 0) then
          do i = 1, thisOctal%nChildren, 1
             child => thisOctal%child(i)
             call calculateResiduals(child, grid, ndepth)
          end do
       endif
    endif
  end subroutine calculateResiduals

  recursive subroutine gSweep2new(thisOctal, grid, fracChange, fracChange2,it,deltaT,doOnlyChanged)
    use inputs_mod, only : smallestCellSize
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    logical :: doOnlyChanged
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6), probe(6)
    integer :: n, ndir
    real(double) :: x1, x2
    real(double) ::  g(6), dx, dxArray(6), g2(6), phiInterface(6),ptemp(6), frac2
    real(double) :: deltaT, fracChange, newPhi, newerPhi, frac, d2phidx2(3), sumd2phidx2, fracChange2
    integer :: nd
    real(double) :: dfdrbyr
    real(double), parameter :: maxM = 100000.d0
    real(double) :: xnext, oldphi, px, py, pz, rm1, um1, pm1, thisR
    real(double), parameter :: SOR = 1.2d0
    integer :: it, nc


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


    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call gsweep2new(child, grid, fracChange, fracChange2, it, deltaT, doOnlyChanged)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if ((.not.thisOctal%changed(subcell)).and.doOnlyChanged) cycle

          if (.not.associated(thisOctal%adot)) then
             allocate(thisOctal%adot(1:thisOctal%maxChildren))
             thisOctal%adot = 0.d0
          endif

          if (.not.thisOctal%ghostCell(subcell)) then

             locator = subcellCentre(thisOctal, subcell)
             thisR = locator%x*gridDistanceScale
             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + probe(n) * (thisOctal%subcellSize/2.d0+0.1d0*smallestCellSize)
                neighbourOctal => grid%octreeRoot

                call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

                ptemp(n) = phigas
                x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                x2 = vector(px,py,pz).dot.dir(n)

                if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
                   x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                   x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                   dx = x2 - x1
                   dx = dx * gridDistanceScale
                else
                   if (nd == thisOctal%nDepth) then ! coarse/coarse or fine/fine
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   else if (nd > thisOctal%nDepth) then ! coarse cells with a fine boundary
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   else
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n) ! fine cells
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   endif
                endif

                dxArray(n) = dx
                g(n) =   (phigas - thisOctal%phi_gas(subcell))/dxArray(n)
             enddo

!get the gravitational potential values at the cell interface
             do n = 1, nDir
                dx = sign(thisOctal%subcellSize/2.d0,dxArray(n))
                   phiInterface(n) = thisOctal%phi_gas(subcell) +  g(n)*(dx*gridDistanceScale)
             end do


             if (thisOctal%twoD) then
                dfdrbyr = (phiInterface(1) - phiInterface(2))/(thisOctal%subcellSize*gridDistanceScale)
                dfdrbyr = dfdrbyr / thisR
             endif



!calculate the new gradient
             do n = 1, nDir
                dx = sign(thisOctal%subcellSize/2.d0,dxArray(n))
                g2(n) = (phiInterface(n) - thisOctal%phi_gas(subcell))/(dx*gridDistanceScale)
             end do


             dx = thisOctal%subcellSize*gridDistanceScale
             if (thisOctal%twoD) then
               d2phidx2(1) = (g2(1) - g2(2)) / dx 
               d2phidx2(2) = (g2(3) - g2(4)) / dx
               sumd2phidx2 = SUM(d2phidx2(1:2)) + dfdrbyr
             else
                d2phidx2(1) = (g(1) - g(2)) / dx
                d2phidx2(2) = (g(3) - g(4)) / dx
                d2phidx2(3) = (g(5) - g(6)) / dx
                sumd2phidx2 = SUM(d2phidx2(1:3))
             endif

             oldPhi = thisOctal%phi_gas(subcell)
             newerPhi = thisOctal%phi_gas(subcell) + deltaT*(sumd2phidx2 - fourPiTimesbigG * thisOctal%rho(subcell))


             newPhi = (1.d0-SOR)*oldPhi + SOR*newerPhi
             
             if (.not.associated(thisOctal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))

             if (.not.associated(thisOctal%adot)) allocate(thisOctal%adot(1:thisOctal%maxChildren))
             if (thisOctal%phi_gas(subcell) /= 0.d0) then
                frac = abs((oldPhi - newPhi)/oldPhi)
                frac2 = abs(sumd2phidx2 - fourPiTimesbigG * thisOctal%rho(subcell))/ &
                     (fourPiTimesbigG * thisOctal%rho(subcell))

                thisOctal%chiLine(subcell) = frac2
                thisOctal%adot(subcell) = frac
                fracChange = max(frac, fracChange)
                fracChange2 = max(frac2, fracChange2)
             else
                fracChange = 1.d30
             endif

             thisOctal%phi_gas(subcell) = newPhi
             
          endif
       endif
    enddo
  end subroutine gSweep2new

  recursive subroutine calculateResiduals2(thisOctal, grid, nDepth)
    use inputs_mod, only : smallestCellSize
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6), probe(6)
    integer :: n, ndir
    integer :: nDepth
    real(double) :: x1, x2
    real(double) ::  g(6), dx, dxArray(6), g2(6), phiInterface(6),ptemp(6)
    real(double) :: deltaT, d2phidx2(3), sumd2phidx2
    integer :: nd
    real(double) :: dfdrbyr
    real(double), parameter :: maxM = 100000.d0
    real(double) :: xnext, px, py, pz, rm1, um1, pm1, thisR
    real(double), parameter :: SOR = 1.2d0
    integer :: nc


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

    if (thisOctal%nDepth == nDepth) then

       do subcell = 1, thisOctal%maxChildren

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (.not.associated(thisOctal%adot)) then
             allocate(thisOctal%adot(1:thisOctal%maxChildren))
             thisOctal%adot = 0.d0
          endif

          if(thisOctal%threed) then
             deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*thisOctal%subcellSize))**2 / 6.d0
          else
             deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*thisOctal%subcellSize))**2 / 4.d0
          endif



          if (.not.thisOctal%ghostCell(subcell)) then

             locator = subcellCentre(thisOctal, subcell)
             thisR = locator%x*gridDistanceScale
             do n = 1, nDir

                locator = subcellCentre(thisOctal, subcell) + probe(n) * (thisOctal%subcellSize/2.d0+0.1d0*smallestCellSize)
                neighbourOctal => grid%octreeRoot

                call findSubcellLocalLevel(locator, neighbourOctal, neighbourSubcell,ndepth)

                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)

                ptemp(n) = phigas
                x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                x2 = vector(px,py,pz).dot.dir(n)

                if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
                   x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                   x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                   dx = x2 - x1
                   dx = dx * gridDistanceScale
                else
                   if (nd == thisOctal%nDepth) then ! coarse/coarse or fine/fine
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   else if (nd > thisOctal%nDepth) then ! coarse cells with a fine boundary
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   else
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n) ! fine cells
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   endif
                endif

                dxArray(n) = dx
                g(n) =   (phigas - thisOctal%phi_gas(subcell))/dxArray(n)
             enddo

             !get the gravitational potential values at the cell interface
             do n = 1, nDir
                dx = sign(thisOctal%subcellSize/2.d0,dxArray(n))
                phiInterface(n) = thisOctal%phi_gas(subcell) +  g(n)*(dx*gridDistanceScale)
             end do


             if (thisOctal%twoD) then
                dfdrbyr = (phiInterface(1) - phiInterface(2))/(thisOctal%subcellSize*gridDistanceScale)
                dfdrbyr = dfdrbyr / thisR
             endif



             !calculate the new gradient
             do n = 1, nDir
                dx = sign(thisOctal%subcellSize/2.d0,dxArray(n))
                g2(n) = (phiInterface(n) - thisOctal%phi_gas(subcell))/(dx*gridDistanceScale)
             end do


             dx = thisOctal%subcellSize*gridDistanceScale
             if (thisOctal%twoD) then
                d2phidx2(1) = (g2(1) - g2(2)) / dx
                d2phidx2(2) = (g2(3) - g2(4)) / dx
                sumd2phidx2 = SUM(d2phidx2(1:2)) + dfdrbyr
             else
                d2phidx2(1) = (g(1) - g(2)) / dx
                d2phidx2(2) = (g(3) - g(4)) / dx
                d2phidx2(3) = (g(5) - g(6)) / dx
                sumd2phidx2 = SUM(d2phidx2(1:3))
             endif

             if (.not.associated(thisOctal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))

             if (.not.associated(thisOctal%adot)) allocate(thisOctal%adot(1:thisOctal%maxChildren))
             thisOctal%chiLine(subcell) = thisOctal%source(subcell) - sumd2phidx2

          endif
       enddo
    else
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call calculateResiduals2(child, grid, ndepth)
       end do
    endif


  end subroutine calculateResiduals2

  recursive subroutine gaussSeidelSweep2(thisOctal, grid, fracChange, ndepth, onlyCellsWithChildren)
    use inputs_mod, only : smallestCellSize
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6), probe(6)
    integer :: n, ndir
    integer :: nDepth
    logical :: onlyCellsWithChildren
    real(double) :: x1, x2
    real(double) ::  g(6), dx, dxArray(6), g2(6), phiInterface(6),ptemp(6)
    real(double) :: deltaT, fracChange, newPhi, newerPhi, d2phidx2(3), sumd2phidx2
    integer :: nd
    real(double) :: dfdrbyr
    real(double), parameter :: maxM = 100000.d0
    real(double) :: xnext, oldphi, px, py, pz, rm1, um1, pm1, thisR
    real(double), parameter :: SOR = 1.2d0
    integer :: nc


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

    if (thisOctal%nDepth == nDepth) then

       do subcell = 1, thisOctal%maxChildren

          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (onlyCellsWithChildren .and. .not.thisOctal%hasChild(subcell)) cycle
          if (.not.associated(thisOctal%adot)) then
             allocate(thisOctal%adot(1:thisOctal%maxChildren))
             thisOctal%adot = 0.d0
          endif
          if(thisOctal%threed) then
             deltaT =  0.1d0*(2.d0*returnCodeUnitLength(gridDistanceScale*thisOctal%subcellSize))**2 / 6.d0
          else
             deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*thisOctal%subcellSize))**2 / 4.d0
          endif


          if (.not.thisOctal%edgeCell(subcell)) then

             locator = subcellCentre(thisOctal, subcell)
             thisR = locator%x*gridDistanceScale
             do n = 1, nDir

                locator = subcellCentre(thisOctal, subcell) + probe(n) * (thisOctal%subcellSize/2.d0+0.1d0*smallestCellSize)
                neighbourOctal => grid%octreeRoot

                call findSubcellLocalLevel(locator, neighbourOctal, neighbourSubcell, nDepth)

                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                if ((nc > 0).and.(nDepth >= minDepthAMR)) then
                   call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                   call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                        rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, &
                        qViscosity)
                endif



                ptemp(n) = phigas
                x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                x2 = vector(px,py,pz).dot.dir(n)

                if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
                   x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                   x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                   dx = x2 - x1
                   dx = dx * gridDistanceScale

                else


                   if (nd == thisOctal%nDepth) then ! coarse/coarse or fine/fine
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   else if (nd > thisOctal%nDepth) then ! coarse cells with a fine boundary
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   else
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n) ! fine cells
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                      dx = dx * gridDistanceScale

                   endif
                endif

                dxArray(n) = dx
                g(n) =   (phigas - thisOctal%phi_gas(subcell))/dxArray(n)
             enddo

             !get the gravitational potential values at the cell interface
             do n = 1, nDir
                dx = sign(thisOctal%subcellSize/2.d0,dxArray(n))
                phiInterface(n) = thisOctal%phi_gas(subcell) +  g(n)*(dx*gridDistanceScale)
             end do


             if (thisOctal%twoD) then
                dfdrbyr = (phiInterface(1) - phiInterface(2))/(thisOctal%subcellSize*gridDistanceScale)
                dfdrbyr = dfdrbyr / thisR
             endif



             !calculate the new gradient
             do n = 1, nDir
                dx = sign(thisOctal%subcellSize/2.d0,dxArray(n))
                g2(n) = (phiInterface(n) - thisOctal%phi_gas(subcell))/(dx*gridDistanceScale)
             end do


             dx = thisOctal%subcellSize*gridDistanceScale
             if (thisOctal%twoD) then
                d2phidx2(1) = (g2(1) - g2(2)) / dx
                d2phidx2(2) = (g2(3) - g2(4)) / dx
                sumd2phidx2 = SUM(d2phidx2(1:2)) + dfdrbyr
             else
                d2phidx2(1) = (g(1) - g(2)) / dx
                d2phidx2(2) = (g(3) - g(4)) / dx
                d2phidx2(3) = (g(5) - g(6)) / dx
                sumd2phidx2 = SUM(d2phidx2(1:3))
             endif

             oldPhi = thisOctal%phi_gas(subcell)
             newerPhi = thisOctal%phi_gas(subcell) + deltaT*(sumd2phidx2 - thisOctal%source(subcell))


             newPhi = (1.d0-SOR)*oldPhi + SOR*newerPhi

             if (.not.associated(thisOctal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))

             if (.not.associated(thisOctal%adot)) allocate(thisOctal%adot(1:thisOctal%maxChildren))

          if (.not.associated(thisOctal%etaLine)) then
             allocate(thisOctal%etaline(1:thisOctal%maxChildren))
             thisOctal%etaLine = 0.d0
          endif

             thisOctal%phi_gas(subcell) = newPhi

!             if (nDepth == 5) write(*,*) "old phi ",oldphi, " newer phi ",newerPhi, " source ",thisOctal%source(subcell),g(1:6)

             if (sumd2phidx2 /= 0.) then
                fracChange = max(fracChange,abs((thisOctal%source(subcell) - sumd2phidx2)/thisOctal%source(subcell)))
                thisOctal%etaline(subcell) = abs((thisOctal%source(subcell) - sumd2phidx2)/thisOctal%source(subcell))
             endif

          endif
       enddo
    else
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call gaussSeidelSweep2(child, grid, fracChange, ndepth, onlyCellsWithChildren)
       end do
    endif


  end subroutine gaussSeidelSweep2

  recursive subroutine gSweep2residuals(thisOctal, grid, deltaT, fracChange, nDepth)
    use mpi
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas,qViscosity(3,3)
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6), probe(6)
    integer :: nDepth
    integer :: n, ndir
    real(double) :: x1, x2, thisDx
    real(double) ::  g(6), dx, dxArray(6), g2(6), phiInterface(6)
    real(double) :: deltaT, fracChange, gGrav, newPhi, newerPhi, frac, d2phidx2(3), sumd2phidx2
    integer :: nd, nc
    real(double) :: xnext, oldphi, px, py, pz, rm1, um1, pm1, thisR
    real(double), parameter :: SOR = 1.2d0
    gGrav = bigG * lengthToCodeUnits**3 / (massToCodeUnits * timeToCodeUnits**2)

    thisDx = amrGridSize /dble(2**(nDepth+1))

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call gsweep2residuals(child, grid, deltaT, fracChange, ndepth)
       end do
    else

       do subcell = 1, thisOctal%maxChildren

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

             locator = subcellCentre(thisOctal, subcell)
             thisR = returnCodeUnitLength(locator%x*gridDistanceScale)
             do n = 1, nDir
                
                locator = subcellCentre(thisOctal, subcell) + probe(n) * (thisOctal%subcellSize/2.d0+0.1d0*grid%halfSmallestSubcell)
                neighbourOctal => thisOctal

                call findSubcellLocalLevel(locator, neighbourOctal, neighbourSubcell, nDepth)
                
                call getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, probe(n), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                
                x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                

                if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then
                   x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                   x2 = subcellCentre(neighbourOctal, neighbourSubcell).dot.dir(n)
                   dx = x2 - x1
                else
                   if (nd == thisOctal%nDepth) then ! coarse/coarse or fine/fine
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1

                   else if (nd > thisOctal%nDepth) then ! coarse cells with a fine boundary
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n)
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                   else
                      x1 = subcellCentre(thisOctal, subcell).dot.dir(n) ! fine cells
                      x2 = VECTOR(px, py, pz).dot.dir(n)
                      dx = x2 - x1
                   endif
                endif
                dxArray(n) = dx
                if (thisOctal%threed) then
                   g(n) =   (phigas - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                else
                   select case(n)
                   case(1, 2)
                      g(n) = gridDistanceScale * (x2 * phigas - x1 * thisOctal%phi_gas(subcell)) / &
                           (returnCodeUnitLength(dx*gridDistanceScale))
                   case(3,4)
                      g(n) = (phigas - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                   end select
                endif
             enddo

!get the gravitational potential values at the cell interface
             do n = 1, nDir
                dx = thisOctal%subcellSize/2.d0
                if (thisOctal%threed) then
                   phiInterface(n) = thisOctal%phi_gas(subcell) +  g(n)*(&
                        returnCodeUnitLength(dx*gridDistanceScale))
                else
                   select case(n)
                   case(1, 2)
                   phiInterface(n) = thisOctal%phi_gas(subcell)*thisR  +  g(n)*(&
                        returnCodeUnitLength(dx*gridDistanceScale))
                   case(3,4)
                      phiInterface(n) = thisOctal%phi_gas(subcell) +  g(n)*(&
                           returnCodeUnitLength(dx*gridDistanceScale))
                   end select
                endif



             end do

!calculate the new gradient
             do n = 1, nDir
                dx = thisOctal%subcellSize/2.d0
                if (thisOctal%threed) then
                   g2(n) = (phiInterface(n) - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                else
                   select case(n)
                   case(1, 2)
                   g2(n) = (phiInterface(n) - thisR*thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                   case(3, 4)
                   g2(n) = (phiInterface(n) - thisOctal%phi_gas(subcell))/(returnCodeUnitLength(dx*gridDistanceScale))
                   end select
                endif
             end do

             dx = thisOctal%subcellSize
            if (thisOctal%twoD) then
               d2phidx2(1) = (1.d0/(thisR))*(g2(1) - g2(2)) / (returnCodeUnitLength(dx*gridDistanceScale))
               d2phidx2(2) = (g2(3) - g2(4)) / (returnCodeUnitLength(dx*gridDistanceScale))
               sumd2phidx2 = SUM(d2phidx2(1:2))
             else
                d2phidx2(1) = (g2(1) - g2(2)) / (returnCodeUnitLength(dx*gridDistanceScale))
                d2phidx2(2) = (g2(3) - g2(4)) / (returnCodeUnitLength(dx*gridDistanceScale))
                d2phidx2(3) = (g2(5) - g2(6)) / (returnCodeUnitLength(dx*gridDistanceScale))
                sumd2phidx2 = SUM(d2phidx2(1:3))
             endif

             if(thisOctal%threed) then
                deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*thisDx))**2 / 6.d0
             else
                deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*thisDx))**2 / 4.d0
             endif

             oldPhi = thisOctal%phi_gas(subcell)
             newerPhi =  (deltaT * sumd2phidx2 - thisOctal%chiLine(subcell) * deltaT) 

             newPhi = (1.d0-SOR)*oldPhi + SOR*newerPhi
             if (.not.associated(thisOCtal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))

!             if (myrankglobal == 1) write(*,*) newphi,thisOctal%phi_gas(subcell),deltaT, sumd2phidx2, thisOctal%rho(subcell)
             if (thisOctal%phi_gas(subcell) /= 0.d0) then

                frac = abs((sumd2phidx2 - Thisoctal%chiline(subcell))/(thisOctal%chiLine(subcell)))

                fracChange = max(frac, fracChange)



             else
                fracChange = 1.d30
!                frac = fracChange
             endif
             if (.not.associated(thisOctal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))
!             thisOctal%chiline(subcell) = frac

             thisOctal%phi_gas(subcell) = newPhi
             
          endif
       enddo
    endif
  end subroutine gSweep2residuals




  recursive subroutine gSweepLevel(thisOctal, grid, deltaT, fracChange, fracChange2, ghostFracChange, nDepth)
    use mpi
    integer :: nd
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer   :: neighbourOctal
    type(octal), pointer  :: child 
    integer :: nDepth
    real(double) :: rho, rhou, rhov, rhow, q, qnext, x, rhoe, pressure, flux, phi, phigas
    integer :: subcell, i, neighbourSubcell
    type(VECTOR) :: locator, dir(6)
    integer :: n, ndir, nc
    real(double) :: g(6) 
    real(double) :: deltaT, fracChange, gGrav, newPhi, frac, ghostFracChange !, d2phidx2(3), sumd2phidx2
    real(double) :: sorFactor, deltaPhi, dx, sumd2phidx2, fracChange2
    real(double) :: xnext, px, py, pz,qViscosity(3,3), rm1, um1, pm1
    sorFactor = 1.2d0

    gGrav = bigG * lengthToCodeUnits**3 / (massToCodeUnits * timeToCodeUnits**2)


    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call gsweepLevel(child, grid, deltaT, fracChange, fracChange2, ghostFracChange, ndepth)
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
                     rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
                g(n) = phigas
             enddo
             dx = thisOCtal%subcellSize
             
             if (thisOctal%twoD) then
                newphi = 0.25d0*(SUM(g(1:4)) - fourpi*gGrav*thisOctal%rho(subcell))
                
             else
                newphi = 0.16666666666667d0*(SUM(g(1:6)) - &
                     fourpi*gGrav*thisOctal%rho(subcell)*(returnCodeUnitLength(dx*gridDistanceScale))**2)
             endif
             
             deltaPhi = newPhi - thisOctal%phi_gas(subcell)
             newPhi = thisOctal%phi_gas(subcell) + sorFactor * deltaPhi
             
             sumd2phidx2 = (sum(g(1:6)) - 6.d0*newPhi)/(thisOctal%subcellSize*gridDistanceScale)**2
             
             if (thisOctal%phi_gas(subcell) /= 0.d0) then
!                frac = abs((sumd2phidx2 - fourPi*gGrav*thisOctal%rho(subcell))/(fourPi*gGrav*thisOctal%rho(subcell)))
                frac = abs((newPhi - thisOctal%phi_gas(subcell))/thisOctal%phi_gas(subcell))
             if (.not.associated(thisOctal%chiline)) allocate(thisOctal%chiline(1:thisOctal%maxChildren))
                thisOctal%chiLine(subcell) = frac
                fracChange = max(frac, fracChange)
             else
                fracChange = 1.d30
             endif
             frac = abs((sumd2phidx2 - fourPi*gGrav*thisOctal%rho(subcell))/(fourPi*gGrav*thisOctal%rho(subcell)))
             fracChange2 = max(frac, fracChange2)
             thisOCtal%phi_gas(subcell) = newPhi             
          endif
       enddo
    endif
  end subroutine gSweepLevel
  
  subroutine multiGridVcycle(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: nPairs, thread1(:), thread2(:), nBound(:), group(:), nGroup
    integer :: iDepth
    real(double) :: fracChange, tempFracChange, deltaT, temp
    integer :: ierr, iter, bigIter, j

    character(len=80) :: plotfile

    call updateDensityTree(grid%octreeRoot)
    call updatePhiTree(grid%octreeRoot, 4)

    call unsetGhosts(grid%octreeRoot)
    do iDepth = 3, maxDepthAMR
       call setupEdgesLevel(grid%octreeRoot, grid, iDepth)
       call setupGhostsLevel(grid%octreeRoot, grid, iDepth)
    enddo

    if (myrankWorldglobal == 1) call tune(6,"Dirichlet boundary conditions")
    call applyDirichlet(grid)
    if (myrankWorldglobal == 1) call tune(6,"Dirichlet boundary conditions")


    fracChange = 1d3
    bigIter = 0
    do while((fracChange > 0.01d0).or.(bigIter < 20))
       call writeInfo("Starting V-cycle", TRIVIAL)

       ! a few gauss-seidel sweeps at maximum depth


       call setSourceTofourPiRhoG(grid%octreeRoot, maxDepthAMR)
       call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, maxDepthAMR)

       fracChange = 1.d3
       iter = 0 

       do while(fracChange > 1.d-20)
          iter = iter + 1
          fracChange = 0.d0
          call gaussSeidelSweep2(grid%octreeRoot, grid, fracChange, maxDepthAMR, onlyCellsWithChildren = .false.)
          call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, maxDepthAMR)
          call MPI_ALLREDUCE(fracChange, tempFracChange, 1, MPI_DOUBLE_PRECISION, MPI_MAX, amrCOMMUNICATOR, ierr)
          fracChange = tempFracChange
          if (writeoutput) write(*,*) "Depth ",maxDepthAMR, " iter ",iter, " frac ",fracchange
          if (iter == 5) exit
       enddo
       call calculateResiduals2(grid%octreeRoot, grid, maxDepthAMR)




       ! now the down part of the V-cycle

       call writeInfo("Starting down part of V-cycle", TRIVIAL)
       do iDepth = maxDepthAMR-1, 3,-1

          if (iDepth > (minDepthAMR-1)) then

! first we need to iterate over residuals at this level for cells that have children

             call restrictResiduals(grid%octreeRoot, iDepth+1) ! from the depth deeper to this depth
             call setSourceToResiduals(grid%octreeRoot, iDepth)
             call zeroPhiGasLevelEverywhere(grid%octreeRoot, iDepth)
             
             call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
             call writeInfo("Iterating residuals", TRIVIAL)
             iter = 0 
             fracChange = 1.e3
             do while(fracChange > 1.d-2)
                iter = iter + 1
                fracChange = 0.d0
                call gaussSeidelSweep2(grid%octreeRoot, grid, fracChange, iDepth, onlyCellsWithChildren=.true.)
                call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
                call MPI_ALLREDUCE(fracChange, tempFracChange, 1, MPI_DOUBLE_PRECISION, MPI_MAX, amrCOMMUNICATOR, ierr)
                fracChange = tempFracChange
                if (writeoutput) write(*,*) "Depth ",iDepth, " iter ",iter, " frac ",fracchange
                if (iter == 50)exit
             enddo
             
             call addPhiFromBelow(grid%octreeRoot, iDepth)
             call restorePhiGas(grid%octreeRoot, iDepth)
             ! now cells at idepth have phi_gas as true phi_gas and not errors on phigas
             ! now we set source as fourpi rho G
             
             call setSourceTofourPiRhoG(grid%octreeRoot, iDepth)
             call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
             call writeInfo("Iterating phis", TRIVIAL)
             
             ! now iterate over this for ALL cells at this level and calculate the residuals
             iter = 0 
             fracChange = 1.d3
             do while(fracChange > 1.d-2)
                iter = iter + 1
                fracChange = 0.d0
                call gaussSeidelSweep2(grid%octreeRoot, grid, fracChange, iDepth, onlyCellsWithChildren = .false.)
                call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
                call MPI_ALLREDUCE(fracChange, tempFracChange, 1, MPI_DOUBLE_PRECISION, MPI_MAX, amrCOMMUNICATOR, ierr)
                fracChange = tempFracChange
                if (writeoutput) write(*,*) "Depth ",iDepth, " iter ",iter, " frac ",fracchange
                if (iter == 50) exit
             enddo
             call calculateResiduals2(grid%octreeRoot, grid, iDepth)
             
             ! now we repeat
             
          else
             ! since all cells are at this level of refinement or above we can do something simpler


             call restrictResiduals(grid%octreeRoot, iDepth+1) ! from the depth deeper to this depth
             call setSourceToResiduals(grid%octreeRoot, iDepth)
             call zeroPhiGasLevelEverywhere(grid%octreeRoot, iDepth)
             
             call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
             call writeInfo("Iterating residuals", TRIVIAL)
             iter = 0 
             fracChange = 1.d3
             do while(fracChange > 1.d-2)
                iter = iter + 1
                fracChange = 0.d0
                call gaussSeidelSweep2(grid%octreeRoot, grid, fracChange, iDepth, onlyCellsWithChildren=.false.)
                call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
                call MPI_ALLREDUCE(fracChange, tempFracChange, 1, MPI_DOUBLE_PRECISION, MPI_MAX, amrCOMMUNICATOR, ierr)
                fracChange = tempFracChange
                if (writeoutput) write(*,*) "Depth ",iDepth, " iter ",iter, " frac ",fracchange
                if (iter == 50) exit
             enddo
             call calculateResiduals2(grid%octreeRoot, grid, iDepth)

          endif

       enddo

       ! now the up part of the V-cycle
       call writeInfo("Beginning up part of V-cycle", TRIVIAL)

       do iDepth = 4, minDepthAMR

          fracChange = 1.d3
          iter = 0 

          call addCorrections(grid%octreeRoot, idepth) ! from depth above

          call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)

          fracChange = 1.d3
          do while(fracChange > 1.d-20)
             iter = iter + 1
             fracChange = 0.d0
             call gaussSeidelSweep2(grid%octreeRoot, grid, fracChange, iDepth, onlyCellsWithChildren=.false.)
             call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
             call MPI_ALLREDUCE(fracChange, tempFracChange, 1, MPI_DOUBLE_PRECISION, MPI_MAX, amrCOMMUNICATOR, ierr)
             fracChange = tempFracChange
             if (writeoutput) write(*,*) "Depth ",iDepth, " iter ",iter, " frac ",fracchange
             if (iter == 50) exit
          enddo

       enddo
       bigiter = bigiter+1
    enddo

    do

       call addCorrections(grid%octreeRoot, minDepthAMR) ! from depth above

!       do iDepth = minDepthAMR, maxDepthAMR
!          iter = 0 
!          fracChange = 1.d3
!          do while(fracChange > 0.01d0)
!             iter = iter + 1
!             fracChange = 0.d0
!             call gaussSeidelSweep2(grid%octreeRoot, grid, fracChange, iDepth, onlyCellsWithChildren = .false.)
!             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!             call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
!             call MPI_ALLREDUCE(fracChange, tempFracChange, 1, MPI_DOUBLE_PRECISION, MPI_MAX, amrCOMMUNICATOR, ierr)
!             fracChange = tempFracChange
!             if (writeoutput) write(*,*) "Depth ",iDepth, " iter ",iter, " frac ",fracchange
!             if (iter == 2) exit
!          enddo
!          call prolongatePhiGas(grid%octreeRoot, iDepth+1)! from depth above
!       enddo

       do iDepth = minDepthAMR, maxDepthAMR
          call setSourceToFourPiRhoG(grid%octreeRoot,iDepth)
       enddo

       call unsetGhosts(grid%octreeRoot)
       call setupEdges(grid%octreeRoot, grid)
       call setupGhosts(grid%octreeRoot, grid)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


       deltaT =  0.1d0*(2.d0*returnCodeUnitLength(gridDistanceScale*smallestCellSize))**2 / 6.d0

       fracChange = 1.d3
       iter = 0 
       do while(fracChange > 0.01)
          iter = iter + 1
          fracChange = 0.d0
!          call  gaussSeidelSweep2(grid%octreeRoot, grid, fracChange, maxDepthAMR, onlyCellsWithChildren=.false.)

          call gSweep2new(grid%octreeRoot, grid, temp, fracChange, j, deltaT, doOnlyChanged=.false.)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          call MPI_ALLREDUCE(fracChange, tempFracChange, 1, MPI_DOUBLE_PRECISION, MPI_MAX, amrCOMMUNICATOR, ierr)
          fracChange = tempFracChange
          if (writeoutput) write(*,*) "iteration ",iter, " fracchange ", fracchange
       write(plotfile,'(a,i4.4,a)') "smallgrav",bigiter*100+iter,".vtk"
       call writeVtkFile(grid, plotfile, &
            valueTypeString=(/"phigas   ", "rho      ","chiline  "/))
          if (iter == 10) exit
       enddo
       bigIter = bigIter + 1
       write(plotfile,'(a,i4.4,a)') "grav",bigiter,".vtk"
       call writeVtkFile(grid, plotfile, &
            valueTypeString=(/"phigas   ", "rho      ","chiline  "/))
       call writeInfo("Done.", TRIVIAL)
    enddo

  end subroutine multiGridVcycle
  
  
  recursive subroutine simpleGravity(thisOctal)
    use inputs_mod, only : sourceMass, sourcePos
    type(OCTAL), pointer :: thisOctal, child
    type(VECTOR) :: rVec
    integer :: subcell, i
    real(double) :: R, thisMass
    
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call simpleGravity(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          Rvec = subcellcentre(thisOCtal, subcell)
          R= modulus(Rvec - sourcePos(1))*1.d10
          thisMass = 1.d0
!          thisMass = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell)!*1.d30
!          thisOctal%phi_gas(subcell) = - bigG*sourcemass(1)/R
!          thisOctal%phi_i(subcell) = - bigG*sourcemass(1)/R
          thisOctal%phi_gas(subcell) = - bigG*sourcemass(1)/R**2
          thisOctal%phi_i(subcell) = - bigG*sourcemass(1)/R**2
          thisOctal%phi_stars(subcell) = 0.d0
       end if
    end do

  end subroutine simpleGravity

  subroutine selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid, onlyChanged)
    use inputs_mod, only :  maxDepthAMR, dirichlet, simpleGrav
    use mpi
    type(gridtype) :: grid
    logical, optional :: multigrid, onlyChanged
    logical :: doOnlyChanged
    integer, parameter :: maxThreads = 512
    integer :: iDepth
    integer :: nPairs, thread1(:), thread2(:), nBound(:), group(:), nGroup
    real(double) :: fracChange2(maxthreads), fracChange(maxthreads), &
         ghostFracChange(maxthreads), tempFracChange(maxthreads), deltaT, dx,taumin
    real(double)  :: tol = 1.d-4,  tol2 = 1.d-5
    integer :: it, ierr, minLevel
!    character(len=80) :: plotfile


    if (amr2d) then
       tol = 1.d-8
       tol2 = 1.d-8
    endif

    if (amr3d) then
       tol = 1.d-6
       tol2 = 1.d-5
    endif

    doOnlyChanged = .false.
    if (PRESENT(onlyChanged)) doOnlyChanged = onlyChanged

!    call multiGridVcycle(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!    goto 666


    if(simpleGrav) then
       call simpleGravity(grid%octreeRoot)
       if (myRankWorldGlobal == 1) write(*,*) "Done simplified self gravity calculation!"
       goto 666
    endif
! endif




    if (amr1d) then
       call selfGrav1D(grid)
       goto 666
    endif


    call updateDensityTree(grid%octreeRoot)

    if (PRESENT(multigrid)) then

       minLevel = 4
       do iDepth = minLevel, minDepthAMR

          call unsetGhosts(grid%octreeRoot)
          call setupEdgesLevel(grid%octreeRoot, grid, iDepth)
          call setupGhostsLevel(grid%octreeRoot, grid, iDepth)

         call updatePhiTree(grid%octreeRoot, iDepth)
         call setSourceTofourPiRhoG(grid%octreeRoot, iDepth)


          if (cylindricalHydro) then
             call imposeBoundaryLevel(grid%octreeRoot, grid, iDepth, phi_gas=.true.)
             call transferTempStorageLevel(grid%octreeRoot, iDepth, justGrav = .true.)
          endif



          if(dirichlet) then
             if (myrankWorldglobal == 1) call tune(6,"Dirichlet boundary conditions")
             call applyDirichlet(grid, iDepth)
             if (myrankWorldglobal == 1) call tune(6,"Dirichlet boundary conditions")
          end if

          dx = returnCodeUnitLength(gridDistancescale*grid%octreeRoot%subcellSize/dble(2.d0**(iDepth-1)))
          if (grid%octreeRoot%twoD) then
             deltaT =  (dx)**2 / 4.d0
          else
             deltaT =  (dx)**2 / 6.d0
          endif


!          deltaT = deltaT * timeToCodeUnits
          it = 0
          fracChange = 1.d30
          fracChange2 = 1.d30
          do while (ANY(fracChange(1:nHydrothreadsGlobal) > tol)) 
             fracChange = 0.d0
             fracChange2 = 0.d0
             ghostFracChange = 0.d0
             it = it + 1
              
!             if (myrankglobal == 1) call tune(6,"Boundary exchange")
             call exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, iDepth)
!             if (myrankglobal == 1) call tune(6,"Boundary exchange")
              
!             if (myrankglobal == 1) call tune(6,"Gsweep level")
             if (cylindricalHydro) then
                call  gSweep2level(grid%octreeRoot, grid, fracChange(myrankGlobal), fracChange2(myrankGlobal), idepth)
             else
                call gSweepLevel(grid%octreeRoot, grid, deltaT, fracChange(myRankGlobal), fracChange2(myrankGlobal), &
                     ghostFracChange(myRankGlobal),iDepth)
             endif
!             if (myrankglobal == 1) call tune(6,"Gsweep level")

!             if (myrankglobal == 1) call tune(6,"Periodic boundary")


            if (cylindricalHydro) then
               call imposeBoundaryLevel(grid%octreeRoot, grid, iDepth, phi_gas=.true.)
               call transferTempStorageLevel(grid%octreeRoot, iDepth, justGrav = .true.)
            endif


             call MPI_ALLREDUCE(fracChange, tempFracChange, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR &
                  , ierr)
             fracChange = tempFracChange


             call MPI_ALLREDUCE(fracChange2, tempFracChange, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR &
                  , ierr)
             fracChange2 = tempFracChange
!             if (myrankGlobal == 1) write(*,*) "Multigrid iteration ",it, " maximum fractional change ", &
!                  MAXVAL(fracChange(1:nHydroThreadsGlobal)),tol


!             if (writeoutput) write(*,*) it,MAXVAL(fracChange2(1:nHydroThreadsGlobal)),tol
          enddo



          if (myRankWorldGlobal == 1) write(*,*) "Gsweep of depth ", iDepth, " done in ", it, " iterations"

       enddo
    endif
 
! Used to bail out after multiGridVcycle      
!555 continue
 
    call unsetGhosts(grid%octreeRoot)
    call setupEdges(grid%octreeRoot, grid)
    call setupGhosts(grid%octreeRoot, grid)


    if(dirichlet.and.(.not.doOnlyChanged)) then
       if (myrankWorldglobal == 1) call tune(6,"Dirichlet boundary conditions")
       call applyDirichlet(grid)
       if (myrankWorldglobal == 1) call tune(6,"Dirichlet boundary conditions")
    end if


    if (cylindricalHydro) then
       call imposeBoundary(grid%octreeRoot, grid, phi_gas=.true.)
       call transferTempStorage(grid%octreeRoot, justGrav = .true.)
    endif



    if (grid%octreeRoot%twoD) then
       deltaT =  (2.d0*returnCodeUnitLength(gridDistancescale*grid%halfSmallestSubcell))**2 / 4.d0
    else
       deltaT =  (2.d0*returnCodeUnitLength(gridDistanceScale*grid%halfSmallestSubcell))**2 / 6.d0
    endif
!    deltaT = deltaT * timeToCodeUnits

    fracChange = 1.d30
    fracChange2 = 1.d30
    it =0 

    do while (ANY(fracChange(1:nHydrothreadsGlobal) > tol2))
       fracChange = 0.d0


       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       fracChange = 0.d0
       fracChange2 = 0.d0


       tauMin  =  0.5*(1.d0/dble(2**maxDepthAMR))**2
       if (amr2d) then
          deltaT = tauMin * (gridDistanceScale * amrGridSize)**2 / 2.d0
       else
          deltaT = tauMin * (gridDistanceScale * amrGridSize)**2 / 4.d0
       endif


       call gSweep2new(grid%octreeRoot, grid, fracChange(myRankGlobal), fracChange2(myRankGlobal), it, deltaT, doOnlyChanged)
       it = it + 1

       if (cylindricalHydro) then
          call imposeBoundary(grid%octreeRoot, grid, phi_gas=.true.)
          call transferTempStorage(grid%octreeRoot, justGrav = .true.)
       endif

       call MPI_ALLREDUCE(fracChange, tempFracChange, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR, ierr)
       fracChange = tempFracChange
       call MPI_ALLREDUCE(fracChange2, tempFracChange, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR, ierr)
       fracChange2 = tempFracChange

       !       write(plotfile,'(a,i4.4,a)') "grav",it,".png/png"
!           if (myrankglobal == 1)   write(*,*) it,MAXVAL(fracChange(1:nHydroThreadsGlobal))

!       if (myrankWorldGlobal == 1) write(*,*) "Full grid iteration ",it, " maximum fractional change ", &
!            MAXVAL(fracChange(1:nHydroThreadsGlobal))

!       if (mod(it,10) == 0) then
!          write(plotfile,'(a,i4.4,a)') "grav",it,".vtk"
!          call writeVtkFile(grid, plotfile, &
!               valueTypeString=(/"phigas ", "rho    ","chiline","adot   "/))
!       endif

!       if (writeoutput) write(*,*) it," frac change ",maxval(fracChange(1:nHydroThreadsGlobal)),tol2,maxval(fracChange2(1:nHydroThreadsGlobal))
!       if (writeoutput) write(*,*) it," frac change ",maxval(fracChange2(1:nHydroThreadsGlobal)),tol2
       if (it > 100000) then
          if (Writeoutput) write(*,*) "Maximum number of iterations exceeded in gravity solver",it
          exit
       endif
    enddo
    if (myRankWorldGlobal == 1) write(*,*) "Gravity solver completed after: ",it, " iterations"


!          call writeVtkFile(grid, "grav.vtk", &
!               valueTypeString=(/"phigas ", "rho    ","chiline"/))

666 continue

!    if (myrankglobal == 1) call tune(6,"Complete self gravity")


  end subroutine selfGrav

  real(double) function getPressure(thisOctal, subcell)
#ifdef PHOTOION
    use inputs_mod, only : photoionPhysics, honly, simpleMu, radpressureTest, mu, ndensity
    use ion_mod, only : nGlobalIon, globalIonArray, returnMu, returnmusimple
#endif
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: eKinetic, eThermal, K, u2, eTot
    real(double), parameter :: gamma2 = 1.4d0, rhoCrit = 1.d-14
    real(double) :: rhoPhys, gamma, cs, rhoNorm, rhov
    type(VECTOR) :: rVec
    logical, save :: firstTime = .true.

    select case(thisOctal%iEquationOfState(subcell))
       case(0) ! adiabatic
          if (thisOctal%threed) then
             u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
          else
             if (.not.cylindricalhydro) then
                u2 = (thisOctal%rhou(subcell)**2 +  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
             else
                rVec = subcellCentre(thisOctal, subcell)
                rhov = thisOctal%rhov(subcell) / (rVec%x * gridDistanceScale)
                u2 = (thisOctal%rhou(subcell)**2 +  rhov**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
             endif

          endif
          eKinetic = u2 / 2.d0
          eTot = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
          eThermal = eTot - eKinetic
          if (eThermal < 0.d0) then
             if (firstTime) then
                rVec = subcellCentre(thisOctal, subcell)
                write(*,'(a18, 1pe9.2)') "eThermal problem: ",ethermal
                write(*,'(1pe9.2,1pe9.2,1pe9.2)') eTot, ekinetic, u2
                write(*,'(a8,1pe9.2)') "3/2nkT ", 1.5d0*(1.d0/(2.33d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
!                write(*,'(a8,1pe9.2)') "1/2mv2 ", 0.5d0 * (thisOctal%rho(subcell) * cellVolume(thisOctal, subcell)*1.d30) * &
!                     (modulus(thisOctal%velocity(subcell))*cs)**2
                write(*,'(a10,1pe9.2,1pe9.2)') "rhoe, rho", thisOctal%rhoe(subcell), thisOctal%rho(subcell)
                write(*,'(1pe9.2,1pe9.2,1pe9.2)') thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                     thisOctal%rhov(subcell)/thisOctal%rho(subcell), thisOctal%rhow(subcell)/thisOctal%rho(subcell)
                write(*,'(a4,1pe9.2,1pe9.2,1pe9.2,a2)') "cen ", rvec%x, rvec%y, rvec%z, thisOctal%ghostCell(subcell)
                write(*,*) "temp ",thisOctal%temperature(subcell)
                firstTime = .false.
             endif

             ethermal = TINY(ethermal)
          endif

          getPressure =  (thisOctal%gamma(subcell) - 1.d0) * thisOctal%rho(subcell) * eThermal
!          write(*,*) "gamma, ethermal",thisOctal%gamma(subcell), ethermal
!          write(*,*) "rhou,rhow ", thisOCtal%rhou(subcell), thisOctal%rhow(subcell)
!          write(*,*) "etot, ekinetic, eThermal ",etot, ekinetic, ethermal
       case(1) ! isothermal
          eThermal = thisOctal%rhoe(subcell) / thisOctal%rho(subcell)

          if (photoionPhysics) then
             mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
             if(honly .and. simpleMu) mu = returnMuSimple(thisOctal, subcell)
          endif
          
!          if (writeoutput) write(*,*) "mu ",mu
!          mu = 1.d0

          getPressure =  thisOctal%rho(subcell)
          getPressure =  getpressure/((mu*mHydrogen))

!          if (radpressuretest.and.(thisOctal%rho(subcell)>1.d-30)) then
!             write(*,*) "N ",getPressure
!             write(*,*) "T ",thisOctal%temperature(subcell)
!             write(*,*) "P ",getpressure*kerg*thisOctal%temperature(subcell)
!          endif
          getPressure =  getpressure*kerg*thisOctal%temperature(subcell)

          if (radPressureTest) getPressure =  min(getPressure,nDensity * kerg * thisOctal%temperature(subcell))


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
!          call writeFatal("hydrodynamics_mod: isothermal pressure needs photoionization")
!          STOP

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

       case(4) ! federrath et al.
          rhoNorm = dble(thisOctal%rho(subcell))/1.d-15
          cs = 0.166d5
          if (rhoNorm <= 0.25d0) then
             gamma = 1.d0
          else if ((rhoNorm > 0.25d0).and.(rhoNorm <= 5.d0)) then
             gamma = 1.1d0
          else
             gamma = 4.d0/3.d0
          endif
          getPressure = cs**2 * thisOctal%rho(subcell)**gamma
          

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
    integer :: n, i
    real(double) :: x
    real(double) :: p(0:100)

     select case(n)
       case(0)
          legendre = 1.d0
       case(1)
          legendre = x
       case DEFAULT
          p(0) = 1.d0
          p(1) = x
          do i = 1, n-1
             p(i+1) = ((2.d0*dble(i)+1.d0)*x*p(i) - dble(i)*p(i-1))/dble(i+1)
          enddo
          legendre = p(n)
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
     type(VECTOR) :: com, rVec
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
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           if (.not.thisOctal%ghostCell(subcell)) then
              if (thisOctal%threed) then 
                 com = com + subcellCentre(thisOctal, subcell) * &
                      thisOctal%rho(subcell)*cellVolume(thisOctal,subcell)*1.d30
              else
                 rVec = subcellCentre(thisOctal,subcell)
                 if (rVec%x > 0.d0) then
                    com%z = com%z + rVec%z * &
                         thisOctal%rho(subcell)*cellVolume(thisOctal,subcell)*1.d30
                 endif
              endif
           endif
        endif
     enddo
   end subroutine sumCOM


   subroutine selfGrav1d(grid)
     use mpi
     type(GRIDTYPE) :: grid
     real(double) :: r !, mass, radius, localMass
     integer :: iThread
     integer :: ierr, j
     real(double) :: mass
     integer, parameter :: tag = 54, tag1 = 55

     call findMassOverAllthreads(grid, mass)

     do iThread = 1, nHydroThreadsGlobal
        if (iThread == myRankGlobal) then
           call setupPhi1d(grid, grid%octreeRoot)
           r = 1.d30
           do j = 1, nHydroThreadsGlobal
              if (j /= myrankGlobal) then
                 call MPI_SEND(r, 1, MPI_DOUBLE_PRECISION, j, tag, localWorldCommunicator, ierr)
              endif
           enddo
        else
           call findTotalMassWithinRServer(grid, iThread)
        endif
     enddo
!     do iThread = 1, nHydroThreadsGlobal
!        if (iThread == myRankGlobal) then
!           call integratePhi1d(grid, grid%octreeRoot, mass)
!           r = 1.d30
!           do j = 1, nHydroThreadsGlobal
!              if (j /= myrankGlobal) then
!                 call MPI_SEND(r, 1, MPI_DOUBLE_PRECISION, j, tag1, localWorldCommunicator, ierr)
!              endif
!           enddo
!        else
!           call findTotalPhiOutsideRServer(grid, iThread)
!        endif
!     enddo



   end subroutine selfGrav1d



   recursive subroutine setupPhi1D(grid, thisOctal)
     use mpi
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal, child
     integer :: subcell, i
     type(VECTOR) :: rVec
     real(double) :: r, mass, localmass
     integer :: ithread
     integer :: ierr, j
     integer, parameter :: tag = 54
     integer :: status(MPI_STATUS_SIZE)

     do subcell = 1, thisoctal%maxchildren
        if (thisoctal%haschild(subcell)) then
           ! find the child
           do i = 1, thisoctal%nchildren, 1
              if (thisoctal%indexchild(i) == subcell) then
                 child => thisoctal%child(i)
                 call setupPhi1d(grid, child)
                 exit
              end if
           end do
        else
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           rVec = subcellCentre(thisOctal, subcell)
           if (.not.rVec%x < 0.d0) then
              r = rVec%x
              mass = 0.d0
              localMass = 0.d0
              
              do iThread = 1, nHydroThreadsGlobal
                 if (myRankGlobal == iThread) then
                    localMass = 0.d0
                    call findTotalMassWithinRMPIPrivate(grid%octreeRoot, r, localMass)
                 else
                    call MPI_SEND(r, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                 endif
              enddo
              mass = mass + localMass
              
              do j = 1, nHydroThreadsGlobal
                 if (j /= myrankGlobal) then
                    call MPI_RECV(localMass, 1, MPI_DOUBLE_PRECISION, j, tag, localWorldCommunicator, status, ierr)
                    mass = mass + localMass
                 endif
              enddo
!              thisOctal%phi_gas(subcell) = (bigG * mass / (r * 1.d10)**2)*(thisOctal%subcellSize*1.d10)
              thisOctal%phi_gas(subcell) = -(bigG * mass / (r * 1.d10))
           endif
        endif
     enddo
   end subroutine setupPhi1d

   recursive subroutine integratePhi1D(grid, thisOctal, gridmass)
     use mpi
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal, child
     integer :: subcell, i
     type(VECTOR) :: rVec
     real(double) :: r, phi, localphi, gridMass
     integer :: ithread
     integer :: ierr, j
     integer, parameter :: tag = 55
     integer :: status(MPI_STATUS_SIZE)

     do subcell = 1, thisoctal%maxchildren
        if (thisoctal%haschild(subcell)) then
           ! find the child
           do i = 1, thisoctal%nchildren, 1
              if (thisoctal%indexchild(i) == subcell) then
                 child => thisoctal%child(i)
                 call integratePhi1d(grid, child, gridMass)
                 exit
              end if
           end do
        else
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           rVec = subcellCentre(thisOctal, subcell)
           if (.not.rVec%x < 0.d0) then
              r = rVec%x
              phi = 0.d0
              localphi = 0.d0
              
              do iThread = 1, nHydroThreadsGlobal
                 if (myRankGlobal == iThread) then
                    localPhi = 0.d0
                    call findTotalPhiOutsideRMPIPrivate(grid%octreeRoot, r, localPhi)
                 else
                    call MPI_SEND(r, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                 endif
              enddo
              phi = phi  + localPhi
              
              do j = 1, nHydroThreadsGlobal
                 if (j /= myrankGlobal) then
                    call MPI_RECV(localphi, 1, MPI_DOUBLE_PRECISION, j, tag, localWorldCommunicator, status, ierr)
                    phi = phi + localPhi
                 endif
              enddo
              thisOctal%phi_gas(subcell) = phi - bigG * gridMass / ((grid%octreeRoot%centre%x+grid%octreeRoot%subcellSize)*1.d10)
              thisOctal%phi_i(subcell) = phi - bigG * gridMass / ((grid%octreeRoot%centre%x+grid%octreeRoot%subcellSize)*1.d10)
           endif
        endif
     enddo
   end subroutine integratePhi1D

     


   recursive subroutine createDensityinEtaline(grid, thisOctal, nDepth)
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal, child, testOctal
     integer :: nDepth
     type(VECTOR) :: rVec
     integer :: subcell, i, testSubcell
     do subcell = 1, thisoctal%maxchildren
        if (thisoctal%haschild(subcell)) then
           ! find the child
           do i = 1, thisoctal%nchildren, 1
              if (thisoctal%indexchild(i) == subcell) then
                 child => thisoctal%child(i)
                 call createDensityinEtaline(grid, child, ndepth)
                 exit
              end if
           end do
        else
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           rVec = subcellCentre(thisOctal, subcell)
           testOctal => grid%octreeRoot
           call findSubcellLocalLevel(rVec, testOctal, testSubcell, nDepth)
           if (.not.associated(thisOctal%etaLine)) allocate(thisOctal%etaLine(1:thisOctal%maxChildren))
           if (.not.associated(thisOctal%chiLine)) allocate(thisOctal%chiLine(1:thisOctal%maxChildren))
           if (.not.associated(thisOctal%adot)) allocate(thisOctal%adot(1:thisOctal%maxChildren))
           thisOctal%etaLine(subcell) = testOctal%phi_gas(testSubcell)
           thisOctal%chiLine(subcell) = testOctal%phi_gas_corr(testSubcell)
        endif
     enddo
   end subroutine createDensityinEtaline



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
              do iPole = 0, 2
                 v = v + (-bigG/x)*(r/x)**ipole * legendre(ipole, cosTheta) * dm
              enddo
              
           endif
        endif
     enddo
   end subroutine multipoleExpansion

   recursive subroutine multipoleExpansionCylindrical(thisOctal, com, Ml, ipole)
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, rVec, rHat
     real(double) :: r, cosTheta, dv, Ml
     integer :: subcell, i, ipole

     do subcell = 1, thisoctal%maxchildren
        if (thisoctal%haschild(subcell)) then
           ! find the child
           do i = 1, thisoctal%nchildren, 1
              if (thisoctal%indexchild(i) == subcell) then
                 child => thisoctal%child(i)
                 call multipoleExpansionCylindrical(child, com, Ml, ipole)
                 exit
              end if
           end do
        else
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
           rVec = subcellCentre(thisOctal, subcell)

           if (.not.thisOctal%ghostCell(subcell).and.(rVec%x > 0.d0)) then

!              rVec = subcellCentre(thisOctal, subcell) - com
              r = modulus(rVec)

              rHat = rVec
              call normalize(rHat)
              cosTheta = VECTOR(0.d0, 0.d0, 1.d0).dot.rHat
              dv = cellVolume(thisOctal, subcell) 
              Ml = Ml + thisOctal%rho(subcell) * r**iPole * legendre(ipole, cosTheta) * dv
              
           endif
        endif
     enddo
   end subroutine multipoleExpansionCylindrical

   recursive subroutine multipoleExpansionCylindricalLevel(thisOctal, com, Ml, ipole, level)
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, rVec, rHat
     real(double) :: r, cosTheta, dv, Ml
     integer :: subcell, i, ipole, level

     if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < level)) then
        do i = 1, thisOctal%nChildren, 1
           child => thisOctal%child(i)
           call  multipoleExpansionCylindricalLevel(child, com, Ml, ipole, level)
        end do
     else
        do subcell = 1, thisOctal%maxChildren
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
           
           rVec = subcellCentre(thisOctal, subcell)
           if (.not.thisOctal%ghostCell(subcell).and.(rVec%x > 0.d0)) then
              
!              rVec = subcellCentre(thisOctal, subcell) - com
              r = modulus(rVec)
              
              rHat = rVec
              call normalize(rHat)
              cosTheta = VECTOR(0.d0, 0.d0, 1.d0).dot.rHat
              dv = cellVolume(thisOctal, subcell)
              Ml = Ml + thisOctal%rho(subcell) * r**iPole * legendre(ipole, cosTheta) * dv
              
          endif
        enddo
     endif
   end subroutine multipoleExpansionCylindricalLevel

   recursive subroutine multipoleExpansionLevel(thisOctal, point, com, v, m, level)
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point, xVec,rVec
     integer :: level
     real(double) :: v, r, x, cosTheta, dm, m
     integer :: subcell, i, ipole

    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < level)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call  multipoleExpansionLevel(child, point, com, v, m, level)
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
              r = r * 1.d10
              x = x * 1.d10
              dm = thisOctal%rho(subcell) * cellVolume(thisOctal,subcell) * 1.d30
              m = m + dm
              do iPole = 0, 2
                 v = v + (-bigG/x)*(r/x)**ipole * legendre(ipole, cosTheta) * dm
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
     real(double) :: v,m
     type(VECTOR) :: point


     if (cylindricalHydro) then
        call applyDirichletCylindrical(grid, level)
     else
        if (PRESENT(level)) then
           call  dirichletQuick(grid, level=level)
        else
           call  dirichletQuick(grid)
        endif
     endif
     goto 666
     
     tag = 94
     call findCoM(grid, com)
!     if (writeoutput) write(*,*) "Centre of Mass found at ",com


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
                 call mpi_send(temp, 3, MPI_DOUBLE_PRECISION, j, tag, localWorldCommunicator, ierr)
              endif
           enddo
        else
           stillReceiving = .true.
           do while (stillReceiving)
              call mpi_recv(temp, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
              if (temp(1) > 1.d29) then
                 stillReceiving = .false.
                 exit
              else
                 point%x = temp(1)
                 point%y = temp(2)
                 point%z = temp(3)
                 v = 0.d0
                 call multipoleExpansionLevel(grid%octreeRoot, point, com, v, m, level=4)
                 call mpi_send(v, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
              endif
           enddo
        endif
        call mpi_barrier(amrCommunicator, ierr)
     enddo
666 continue  
   end subroutine applyDirichlet

   subroutine applyDirichletCylindrical(grid, level)
     use mpi
     type(GRIDTYPE) :: grid
     type(VECTOR) :: com
     real(double) :: temp(1)
     integer, optional :: level
     integer :: iThread
     integer :: tag
     logical :: stillReceiving
     integer :: status(MPI_STATUS_SIZE)
     integer :: ierr, j, ipole
     real(double) :: v


     tag = 94
     call findCoM(grid, com)

     call recursApplyDirichletCylindrical(grid, grid%octreeRoot, com, reset = .true.)
     call recursApplyDirichletLevel2d(grid, grid%octreeRoot, com, 0, reset = .true.)


     do iThread = 1, nHydroThreadsGlobal
        if (myRankGlobal == iThread) then
           if (present(level)) then
!              if (myrankGlobal == 1) write(*,*) "calling cylindrical dirichlet with level ",level
!              if (myrankGlobal == 1) write(*,*) "centre of mass ",com
              call recursApplyDirichletLevel2d(grid, grid%octreeRoot, com, level)
           else
              call recursApplyDirichletCylindrical(grid, grid%octreeRoot, com)
           endif
           do j = 1, nHydroThreadsGlobal
              if (j /= iThread) then
                 temp(1) = 1.d30
                 call mpi_send(temp, 1, MPI_DOUBLE_PRECISION, j, tag, localWorldCommunicator, ierr)
              endif
           enddo
        else
           stillReceiving = .true.
           do while (stillReceiving)
              call mpi_recv(temp, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
              if (temp(1) > 1.d29) then
                 stillReceiving = .false.
                 exit
              else
                 ipole = nint(temp(1))
                 v = 0.d0
                 if (present(level)) then
                    call multipoleExpansionCylindricalLevel(grid%octreeRoot, com, v, ipole, level=level)
                 else
                    call multipoleExpansionCylindrical(grid%octreeRoot, com, v, ipole)
                 endif
                 call mpi_send(v, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
              endif
           enddo
        endif
        call mpi_barrier(amrCommunicator, ierr)
     enddo
   end subroutine applyDirichletCylindrical


   recursive subroutine recursApplyDirichlet(grid, thisOctal, com)
     use mpi
     type(GRIDTYPE) :: grid
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point
     real(double) :: v, vgrid
     real(double) :: temp(3),m,r
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

              if (thisOctal%threed.or.(thisOctal%twoD.and.(point%x > 0.d0))) then
                 do iThread = 1, nHydroThreadsGlobal
                    
                    temp(1) = point%x
                    temp(2) = point%y
                    temp(3) = point%z
                    if (iThread /= myRankGlobal) then
                       call mpi_send(temp, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                    endif
                 enddo
                 v = 0.d0
                 m = 0.d0
                 call multipoleExpansionLevel(grid%OctreeRoot, point, com, v, m, level=4)
                 do iThread = 1, nHydroThreadsGlobal
                    if (iThread /= myRankGlobal) then
                       call mpi_recv(vgrid, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
                       v = v + vgrid
                    endif
                 enddo
                 
!              write(*,*) "old/new phi gas boundary ",thisOctal%phi_gas(subcell), v
                 r = modulus(point)*1.d10
                 thisOctal%phi_gas(subcell) = v
!                 write(*,*) "v ", v, " analytical ", -bigG*msol/r, v/(-bigG*mSol/r)
              endif
           endif
        endif
     enddo
   end subroutine recursApplyDirichlet


   recursive subroutine recursApplyDirichletCylindrical(grid, thisOctal, com, reset)
     use mpi
     type(GRIDTYPE) :: grid
     logical, optional :: reset
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point, uHat
     real(double) :: mgrid, phiB
     real(double) :: temp(1),  muB, rB, panal
     integer :: subcell, i
     integer, parameter :: npole = 6
     integer :: ithread
     integer :: tag, ierr, ipole
     integer :: status(MPI_STATUS_SIZE)
     logical, save :: firstTime = .true.
     real(double), save :: ml(0:nPole)

     tag = 94
     
     if (PRESENT(reset)) then
        if (reset) then
           firstTime = .true.
        endif
        goto 666
     endif


     if (firstTime) then
        do iPole = 0, npole
           
           Ml(ipole) = 0.d0
           call multipoleExpansionCylindrical(grid%OctreeRoot, com, Ml(ipole), ipole)
           
           do iThread = 1, nHydroThreadsGlobal
              
              temp(1) = dble(iPole)
              if (iThread /= myRankGlobal) then
                 call mpi_send(temp, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
              endif
           enddo
           do iThread = 1, nHydroThreadsGlobal
              if (iThread /= myRankGlobal) then
                 call mpi_recv(mgrid, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
                 Ml(ipole) = Ml(ipole) + mgrid
!                 if (writeoutput.and.(ipole==1)) write(*,*) "ml ",mgrid,ml(ipole)
              endif
           enddo
        enddo
        firstTime = .false.
!        if (writeoutput) then
!           write(*,*) "ml ",ml(0:nPole)
!        endif
     endif


     do subcell = 1, thisoctal%maxchildren
        if (thisoctal%haschild(subcell)) then
           ! find the child
           do i = 1, thisoctal%nchildren, 1
              if (thisoctal%indexchild(i) == subcell) then
                 child => thisoctal%child(i)
                 call recursApplyDirichletCylindrical(grid, child, com)
                 exit
              end if
           end do
        else
           if (.not.OctalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           if (thisOctal%ghostCell(subcell)) then
              point = subcellCentre(thisOctal, subcell)

              if (thisOctal%threed.or.(thisOctal%twoD.and.(point%x > 0.d0))) then

                 PhiB = 0.d0
                 do iPole = 0, npole
                    rB = modulus(point)
                    uHat = point
                    call normalize(uHat)
                    muB = VECTOR(0.d0, 0.d0, 1.d0).dot. uHat
                    phiB = phiB - bigG * legendre(iPole,muB)*(rB**(-(ipole+1))) * Ml(ipole)
                 enddo
                 phiB = phiB * 1.d20
                 thisOctal%phi_gas(subcell) = phiB
                 panal = -bigG * 1.d0 * mSol /(rb*1.d10) - bigG *msol / (1.d10*modulus(point - VECTOR(0.d0,0.d0,5.d8)))
!                 if (abs((panal-phiB)/panal) > 0.01d0) then
!                    write(*,*) abs((panal-phiB)/panal),  " phiB ",phiB,panal,thisOctal%nDepth, point
!                 endif
!                 write(*,*) "phiB complete",phiB, " test ",-bigG * 1.d0 * mSol / rB, phiB/(-bigG * 1.d0 * mSol / rB)

              endif
           endif
        endif
     enddo
666 continue
   end subroutine recursApplyDirichletCylindrical

   recursive subroutine recursApplyDirichletlevel2d(grid, thisOctal, com, level, reset)
     use mpi
     type(GRIDTYPE) :: grid
     logical, optional :: reset
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point, uHat
     integer, parameter :: npole = 6
     integer :: level
     real(double) :: mgrid, phiB
     real(double) :: temp(1),  muB, rB, panal
     integer :: subcell, i
     integer :: ithread
     integer :: tag, ierr, ipole
     integer :: status(MPI_STATUS_SIZE)
     logical, save :: firstTime = .true.
     real(double), save :: ml(0:nPole)

     tag = 94
     
     if (PRESENT(reset)) then
        if (reset) then
           firstTime = .true.
           goto 666
        endif
     endif

     if (firstTime) then
        do iPole = 0, nPole
           
           Ml(ipole) = 0.d0
           call multipoleExpansionCylindricalLevel(grid%OctreeRoot, com, Ml(ipole), ipole, level)
!           if (writeoutput.and.(ipole==1)) write(*,*) "ml ",ml(ipole),ml(ipole)

           
           do iThread = 1, nHydroThreadsGlobal
              
              temp(1) = dble(iPole)
              if (iThread /= myRankGlobal) then
                 call mpi_send(temp, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
              endif
           enddo
           do iThread = 1, nHydroThreadsGlobal
              if (iThread /= myRankGlobal) then
                 call mpi_recv(mgrid, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
                 Ml(ipole) = Ml(ipole) + mgrid
!                 if (writeoutput.and.(ipole==1)) write(*,*) "ml ",mgrid,ml(ipole)
              endif
           enddo
                                  
        enddo
!           if (writeoutput) then
!              write(*,*) "ml ",ml(0:nPole)
!           endif
        firstTime = .false.
     endif
                 



     if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < level)) then
        do i = 1, thisOctal%nChildren, 1
           child => thisOctal%child(i)
           call recursApplyDirichletLevel2d(grid, child, com, level) 
        end do
     else
       do subcell = 1, thisOctal%maxChildren
           if (.not.OctalOnThread(thisOctal, subcell, myrankGlobal)) cycle

           if (thisOctal%ghostCell(subcell)) then
              point = subcellCentre(thisOctal, subcell)

              if (thisOctal%threed.or.(thisOctal%twoD.and.(point%x > 0.d0))) then

                 PhiB = 0.d0
                 do iPole = 0, nPole
                    rB = modulus(point)
                    uHat = point
                    call normalize(uHat)
                    muB = VECTOR(0.d0, 0.d0, 1.d0).dot. uHat
                    phiB = phiB - bigG * legendre(iPole,muB)*(rB**(-(ipole+1))) * Ml(ipole)
!		    if (writeoutput) write(*,*) "phib ",ipole, phib
                 enddo
                 phiB = phiB * 1.d20
                 
                 thisOctal%phi_gas(subcell) = phiB
                 panal = -bigG * 1.d0 * mSol /(rb*1.d10) - bigG *msol / (1.d10*modulus(point - VECTOR(0.d0,0.d0,5.d8)))
!                 if (abs((panal-phiB)/panal) > 0.01d0) then
!                    write(*,*) abs((panal-phiB)/panal),  "phiB ",phiB, panal,thisOctal%nDepth,point
!                 endif

              endif
           endif
        enddo
     endif
666 continue
   end subroutine recursApplyDirichletLevel2d


   recursive subroutine recursApplyDirichletLevel(grid, thisOctal, com, level)
     use mpi
     type(GRIDTYPE) :: grid
     integer :: level
     type(OCTAL), pointer :: thisOctal, child
     type(VECTOR) :: com, point
     real(double) :: v, vgrid
     real(double) :: temp(3),r,m
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
              do iThread = 1, nHydroThreadsGlobal

                 temp(1) = point%x
                 temp(2) = point%y
                 temp(3) = point%z
                 if (iThread /= myRankGlobal) then
                    call mpi_send(temp, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                 endif
              enddo
              v = 0.d0
              m = 0.d0
              call multipoleExpansionLevel(grid%OctreeRoot, point, com, v, m, level=4)
              do iThread = 1, nHydroThreadsGlobal
                 if (iThread /= myRankGlobal) then
                    call mpi_recv(vgrid, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
                    v = v + vgrid
                 endif
              enddo
              thisOctal%phi_gas(subcell) = v
              r = modulus(subcellCentre(thisOctal,subcell))*gridDistanceScale
!              write(*,*) "v ", v, " analytical level ", -bigG*msol/r, v/(-bigG*mSol/r)
           endif
        enddo
     endif
   end subroutine recursApplyDirichletLevel


   subroutine dirichletQuick(grid, level)
     use mpi
     type(GRIDTYPE) :: grid
     integer, optional :: level
     integer :: npoints
     type(VECTOR) :: com
     integer, allocatable :: nPointsByThread(:)
     type(VECTOR), allocatable :: points(:)
     real(double),allocatable :: v(:), temp(:)
     real(double) :: m
     integer :: ithread
     integer, parameter :: tag = 45
     integer :: status(MPI_STATUS_SIZE), ierr, np, i
     integer :: maxPoints, startPoint

     maxPoints = 8*(2**minDepthAmr)**2
     allocate(points(1:maxPoints), v(1:maxPoints), nPointsByThread(1:nHydroThreadsGlobal))

     call findCoM(grid, com)

     if (myrankGlobal == 1) then
        nPoints = 0
        if (.not.PRESENT(level)) then
           call getEdgePointsRecur(grid%octreeRoot, nPoints, points)
        else
           call getEdgePointsRecurLevel(grid%octreeRoot, nPoints, points, level)
        endif
        nPointsByThread(1) = nPoints
        do iThread = 2, nHydroThreadsGlobal
           call mpi_recv(np, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, status, ierr)
           call mpi_recv(points((nPoints+1):(nPoints+np))%x, np, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, &
                status, ierr)
           call mpi_recv(points((nPoints+1):(nPoints+np))%y, np, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, &
                status, ierr)
           call mpi_recv(points((nPoints+1):(nPoints+np))%z, np, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, &
                status, ierr)
           nPointsByThread(ithread) = np
           nPoints = nPoints + np
        enddo
     else
        nPoints = 0
        if (.not.PRESENT(level)) then
           call getEdgePointsRecur(grid%octreeRoot, nPoints, points)
        else
           call getEdgePointsRecurLevel(grid%octreeRoot, nPoints, points, level)
        endif
        call mpi_send(nPoints, 1, MPI_INTEGER, 1, tag, localWorldCommunicator, ierr)
        call mpi_send(points(1:nPoints)%x, nPoints, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, ierr)
        call mpi_send(points(1:nPoints)%y, nPoints, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, ierr)
        call mpi_send(points(1:nPoints)%z, nPoints, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, ierr)
     endif
     call mpi_barrier(amrCommunicator,ierr)
     call MPI_BCAST(npoints, 1, MPI_INTEGER, 0, amrCommunicator, ierr)
     call MPI_BCAST(npointsbythread, nHydroThreadsGlobal, MPI_INTEGER, 0, amrCommunicator, ierr)
     call MPI_BCAST(points%x, npoints, MPI_DOUBLE_PRECISION, 0, amrCommunicator, ierr)
     call MPI_BCAST(points%y, npoints, MPI_DOUBLE_PRECISION, 0, amrCommunicator, ierr)
     call MPI_BCAST(points%z, npoints, MPI_DOUBLE_PRECISION, 0, amrCommunicator, ierr)
     v = 0.d0
     do i = 1, nPoints
        call multipoleExpansionLevel(grid%octreeRoot, points(i), com, v(i), m, level=4)
     enddo
     allocate(temp(1:nPoints))
     call MPI_ALLREDUCE(v(1:nPoints), temp, npoints, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
     v(1:npoints) = temp(1:npoints)
     deallocate(temp)
     if (myrankGlobal == 1) then
        startPoint = 1
     else
        startPoint = SUM(nPointsByThread(1:myrankGlobal-1))+1
     endif
     if (.not.present(LEVEL)) then
         call putDirichletPotential(grid%OctreeRoot, v, startPoint)
     else
         call putDirichletPotentialLevel(grid%OctreeRoot, v, startPoint, level)
      endif

      deallocate(points, v, nPointsbythread)
   end subroutine dirichletQuick

   recursive subroutine putDirichletPotential(thisOctal, v, counter)
     type(OCTAL), pointer :: thisOctal, child
     integer :: counter, subcell, i
     real(double) :: v(:), r
     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call putDirichletPotential(child, v, counter)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          if (thisOctal%edgeCell(subcell)) then
             thisOctal%phi_gas(subcell) = v(counter)
             r = modulus(subcellCentre(thisOctal, subcell))
             counter = counter + 1
          endif
       endif
    enddo
  end subroutine putDirichletPotential

   recursive subroutine putDirichletPotentialLevel(thisOctal, v, counter, nDepth)
     type(OCTAL), pointer :: thisOctal, child
     integer :: counter, subcell, i, nDepth
     real(double) :: v(:), r
     if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
        do i = 1, thisOctal%nChildren, 1
           child => thisOctal%child(i)
           call putDirichletPotentialLevel(child, v, counter, nDepth)
        end do
     else

        do subcell = 1, thisOctal%maxChildren
           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
           if (thisOctal%edgeCell(subcell)) then
              thisOctal%phi_gas(subcell) = v(counter)
              r = modulus(subcellCentre(thisOctal, subcell))
              counter = counter + 1
           endif
        enddo
     endif
   end subroutine putDirichletPotentialLevel

  recursive subroutine getEdgePointsRecur(thisOctal, nPoints, points)
     type(OCTAL), pointer :: thisOctal, child
     integer :: nPoints, subcell, i
     type(VECTOR) :: points(:)
     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getEdgePointsRecur(child, nPoints, points)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          if (thisOctal%edgeCell(subcell)) then
             nPoints = npoints + 1
             points(npoints) = subcellCentre(thisOctal, subcell)
          endif
       endif
    enddo
  end subroutine getEdgePointsRecur

  recursive subroutine getEdgePointsRecurLevel(thisOctal, nPoints, points, nDepth)
     type(OCTAL), pointer :: thisOctal, child
     integer :: nPoints, subcell, i, nDepth
     type(VECTOR) :: points(:)

     if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
        do i = 1, thisOctal%nChildren, 1
           child => thisOctal%child(i)
           call  getEdgePointsRecurLevel(child, nPoints, points, nDepth)
        end do
     else

        do subcell = 1, thisOctal%maxChildren

           if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
           if (thisOctal%edgeCell(subcell)) then
              nPoints = npoints + 1
              points(npoints) = subcellCentre(thisOctal, subcell)
           endif
        enddo
     endif
   end subroutine getEdgePointsRecurLevel
  
             


   subroutine mergeSinksFF(grid, source, nSource)
     type(GRIDTYPE) :: grid
     type(SOURCETYPE) :: source(:)
     type(SOURCETYPE), allocatable :: newSource(:)
     type(VECTOR) :: newVelocity, newPosition, newAngMom
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
        newAngMom = VECTOR(0.d0, 0.d0, 0.d0)
        n = 0

        do j = 1, nSource
           if (group(j) == i) then
              newMass = newMass + source(j)%mass
              newVelocity = newVelocity + source(j)%mass * source(j)%velocity
              newPosition = newPosition + source(j)%position * source(j)%mass
              newAge = newAge + source(j)%age
              newAngMom = newAngMom + source(j)%angMomentum
              n = n + 1
           endif
        enddo
        newNSource = newNSource + 1
        newSource(newNSource)%mass = newMass
        newSource(newNSource)%angMomentum = newAngMom
        newSource(newnSource)%velocity = (1.d0/newMass) * newVelocity 
        newSource(newnSource)%position = (1.d0/newMass) * newPosition
        newSource(newnSource)%radius = rSol/1.d10
        newSource(newnSource)%age = newAge / dble(n)
        newSource(newnSource)%stellar = .true.
        newSource(newnSource)%diffuse = .false.
        newSource(newnSource)%accretionRadius = source(1)%accretionRadius
     enddo

        

     nSource = newnSource
     source(1:nSource)%mass = newSource(1:nSource)%mass
     source(1:nSource)%velocity = newSource(1:nSource)%velocity
     source(1:nSource)%position = newSource(1:nSource)%position
     source(1:nSource)%radius = newSource(1:nSource)%radius
     source(1:nSource)%age = newSource(1:nSource)%age
     source(1:nSource)%stellar = newSource(1:nSource)%stellar
     source(1:nSource)%diffuse = newSource(1:nSource)%diffuse

     do j = 1, nSource
        call emptySurface(source(j)%surface)
        call buildSphereNbody(source(j)%position, grid%halfSmallestSubcell, source(j)%surface, 20)
     enddo

     if (writeoutput) then
        write(*,*) "Number of sources after merges: ",newnSource
        call writesourcelist(source, nSource)
        call writeVtkFilenBody(nsource, source, "nbody_temp.vtk")
     endif
666 continue



     deallocate(newSource)
     deallocate(group)
   end subroutine mergeSinksFF

   subroutine mergeSinks(grid, source, nSource)
     use inputs_mod, only : mergeBoundSinks
     type(GRIDTYPE) :: grid
     type(SOURCETYPE) :: source(:)
     type(SOURCETYPE), allocatable :: newSource(:)
     type(VECTOR) :: newVelocity, newPosition, newAngMom
     real(double) :: newMass, newAge
     integer :: nSource
     integer :: i, j
     real(double) :: sep, eGrav, eKinetic, dotProd
     integer :: n
     integer :: newNSource
     logical :: converged
     logical :: merged(1000), bound

     if (.not.mergeBoundSinks) goto 666

     converged = .false.
     do while (.not.converged)
        converged = .true.
        newNSource = 0
        merged = .false.
        allocate(newSource(1:1000))

        do i = 1, nSource
           do j = 1, nSource
              if (i /= j) then
                 sep = modulus(source(i)%position - source(j)%position)
                 eGrav = -bigG * source(i)%mass * source(j)%mass / (sep*1.d10)
                 eKinetic = 0.5d0*source(i)%mass*modulus(source(i)%velocity)**2 + &
                      0.5d0 * source(j)%mass*modulus(source(j)%velocity)**2
                 dotProd = source(i)%velocity .dot. source(j)%velocity
                 bound = .false.
!                 if (writeoutput) write(*,*) "attempting to merge ",i,j,eGrav+eKinetic,dotprod,bound
                 if (((eGrav+eKinetic) < 0.d0).and.(dotProd < 0.d0)) bound = .true.

                 if ((sep < source(i)%accretionRadius).and.bound.and.(.not.(merged(i).or.merged(j)))) then
                    if (writeoutput) write(*,*) i,j," being merged ",eGrav,eKinetic
                    newMass = source(i)%mass + source(j)%mass
                    newVelocity = source(i)%mass * source(i)%velocity + source(j)%mass * source(j)%velocity
                    newPosition = source(i)%mass * source(i)%position + source(j)%mass * source(j)%position
                    newAge = source(i)%age + source(j)%age
                    newAngMom = source(i)%angMomentum + source(j)%angMomentum
                    n = 2
                    newNSource = newNSource + 1
                    newSource(newNSource)%mass = newMass
                    newSource(newNSource)%angMomentum = newAngMom
                    newSource(newnSource)%velocity = (1.d0/newMass) * newVelocity 
                    newSource(newnSource)%position = (1.d0/newMass) * newPosition
                    newSource(newnSource)%radius = rSol/1.d10
                    newSource(newnSource)%age = newAge / dble(n)
                    newSource(newnSource)%stellar = .true.
                    newSource(newnSource)%diffuse = .false.
                    newSource(newnSource)%accretionRadius = source(i)%accretionRadius
                    merged(i) = .true.
                    merged(j) = .true.
                    converged = .false.
                 endif
              endif
           enddo
        enddo
        do i = 1, nSource
           if (.not.merged(i)) then
              newNSource = newNSource + 1
              newSource(newNSource)%mass = source(i)%mass
              newSource(newNSource)%velocity = source(i)%velocity
              newSource(newNSource)%position = source(i)%position
              newSource(newNSource)%radius = source(i)%radius
              newSource(newNSource)%age = source(i)%age
              newSource(newnSource)%stellar = .true.
              newSource(newnSource)%diffuse = .false.
              newSource(newnSource)%accretionRadius = source(i)%accretionRadius
           endif
        enddo
        
        nSource = newnSource
        source(1:nSource)%mass = newSource(1:nSource)%mass
        source(1:nSource)%velocity = newSource(1:nSource)%velocity
        source(1:nSource)%position = newSource(1:nSource)%position
        source(1:nSource)%radius = newSource(1:nSource)%radius
        source(1:nSource)%age = newSource(1:nSource)%age
        source(1:nSource)%stellar = newSource(1:nSource)%stellar
        source(1:nSource)%diffuse = newSource(1:nSource)%diffuse
        source(1:nSource)%angMomentum = newSource(1:nSource)%angMomentum
        
        do j = 1, nSource
           call emptySurface(source(j)%surface)
           call buildSphereNbody(source(j)%position, grid%halfSmallestSubcell, source(j)%surface, 20)
        enddo
        
        if (writeoutput) then
           write(*,*) "Number of sources after merges: ",nSource
        endif
        deallocate(newSource)
     end do
     
666  continue



   end subroutine mergeSinks

   subroutine addSinks(grid, source, nSource)
     use mpi
     type(GRIDTYPE) :: grid
     type(SOURCETYPE) :: source(:)
     real(double) :: rhomax
     integer :: nSource
     integer :: iThread

     do iThread = 1, nHydroThreadsGlobal
        if (myrankGlobal == iThread) then
           rhomax = 0.d0
           call recursaddSinks(grid%octreeRoot, grid, source, nSource, rhomax)
           call shutdownRadiusServer()
        else
           call getPointsInAccretionRadiusServer(grid, nSource, source)
        endif
     enddo
   end subroutine addSinks
          

  recursive subroutine recursaddSinks(thisOctal, grid, source, nSource, rhomax)
    use mpi
    use inputs_mod, only : rhoThreshold, smallestCellSize, accretionRadius
    type(OCTAL), pointer :: thisOctal, child
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    type(VECTOR) :: centre
    real(double) :: bigJ, e, rhomax
    real(double), allocatable :: eGrav(:), eKinetic(:), eThermal(:)
    real(double) :: temp(13), racc
    integer :: nSource
    integer :: i, subcell, ierr, ithread, tag
    logical :: createSink
    type(VECTOR) :: thisVel
    integer, parameter :: maxPoints = 100
    integer :: npoints
    type(VECTOR) :: position(maxPoints), vel(maxPoints),rVec, vcom
    real(double) :: mass(maxPoints), phi(maxPoints), cs(maxPoints)
    real(double) :: cellMass, r, divV
    integer :: iPoint

    tag = 65

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


          if (thisOctal%ghostCell(subcell)) cycle

          centre = subcellCEntre(thisOctal, subcell)
          vel = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
               thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
               thisOctal%rhow(subcell)/thisOctal%rho(subcell))
          bigJ = 0.25d0
          cs = soundSpeed(thisOctal, subcell)
          
          createSink = .true.

          if (thisOctal%rho(subcell) < rhoThreshold) then
             createSink = .false.
             cycle
          endif


          call getPointsInAccretionRadius(thisOctal, subcell, 2.5d0*smallestCellSize, grid, npoints, position, vel, mass, phi, cs)

!          if (createSink) write(*,*) "Source creating passed jeans test ",thisOctal%rho(subcell)/rhoThreshold
          rhomax = max(rhomax, thisOctal%rho(subcell)/rhoThreshold)

          cellMass = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell) * 1.d30
          allocate(eKinetic(1:nPoints), eGrav(1:nPoints), eThermal(1:nPoints))

          vcom = VECTOR(0.d0,0.d0,0.d0)
          do iPoint = 1, nPoints
             vcom = vcom + (mass(iPoint)*vel(iPoint))
          enddo
          vcom = vcom / SUM(mass(1:nPoints))

          do iPoint = 1, nPoints
             eKinetic(iPoint) = 0.5d0 * mass(iPoint)* modulus(vel(ipoint)-vcom)**2
          enddo
          eGrav(1:nPoints) = mass(1:npoints) * phi(1:nPoints)
          eThermal(1:nPoints) = 0.5d0 * mass(1:nPoints) * cs(1:nPoints)**2

          if (SUM(eKinetic) + SUM(eGrav) + SUM(eThermal) > 0.d0) then
             createSink = .false.
             write(*,*) "bound test failed ",SUM(eKinetic)+SUM(eGrav)+SUM(eThermal),SUM(eKinetic), SUM(eGrav), SUM(eThermal)
             if (SUM(eKinetic)+SUM(eGrav)+SUM(eThermal) > 1.d60) then
                do iPoint = 1, nPoints
                   write(*,*) iPoint, " pos ",position(iPoint), " vel ",vel(iPoint), " mass ",mass(ipoint), &
                        " phi ",phi(iPoint), " cs ",cs(iPoint)
                enddo
             endif
          endif

!          endif
          if (abs(SUM(eGrav))< (2.d0*SUM(eThermal))) then
             createSink = .false.
!             write(*,*) "virial test failed ",abs(SUM(eGrav)),2.d0*SUM(ethermal)
          endif

          deallocate(eGrav, eKinetic, eThermal) 
          if (.not.createSink) cycle

!          if (createSink) write(*,*) "Source creating passed virial test ",abs(eGrav),2.d0*eThermal
          do i = 1, nSource
             racc = source(i)%accretionRadius
             if (modulus(centre - source(i)%position)*1.d10 < racc) then
                createsink = .false.
                exit
             endif
          enddo
          if (.not.createSink) cycle
!          if (createSink) write(*,*) "Source creating passed accretion radius test "

          if (createSink) then
             if (MINVAL(phi(2:nPoints)) < phi(1)) then
                createSink = .false.
!                write(*,*) "not phi minumum ",phi(1),minval(phi(1:nPoints))
                cycle
             endif
          endif



          if (createSink) then
             divV = 0.d0
             do iPoint = 2, nPoints
                rVec = position(iPoint) - position(1)
                r = modulus(rVec)
                call normalize(rVec)
                if (r < smallestCellSize*1.01d0) then
                   divV = divV + (rVec.dot.(vel(iPoint)-vel(1)))
                   if ((rVec.dot.(vel(iPoint)-vel(1))) > 0.d0) then
                      createSink = .false.
                      write(*,*) "not converging ",rVec.dot.((vel(ipoint)-vel(1))),rVec
                   endif
                endif
             enddo
             if (.not.createSink) then
                write(*,*) "Create source failed on converging flow test"
             endif
             if (divV < 0.d0) then
                write(*,*) "Source passed divV < 0 test so creating away!"
                createSink = .true.
             endif
          endif

          if (createSink) then
             nSource = nSource + 1
             write(*,*) "Source number ",nsource," created"
             source(nSource)%position = subcellCentre(thisOctal, subcell)
             source(nsource)%velocity%x = thisOctal%rhou(subcell)/thisOctal%rho(subcell)
             source(nsource)%velocity%y = thisOctal%rhov(subcell)/thisOctal%rho(subcell)
             source(nsource)%velocity%z = thisOctal%rhow(subcell)/thisOctal%rho(subcell)             
             source(nsource)%mass = (thisOctal%rho(subcell) - rhoThreshold)*thisOctal%subcellSize**3*1.d30
             source(nsource)%radius = rsol/1.d10
             source(nSource)%accretionRadius = accretionRadius * smallestCellSize * 1.d10
             source(nSource)%age = 0.d0
             source(nSource)%angMomentum = VECTOR(0.d0, 0.d0, 0.d0)
             call buildSphereNbody(source(nsource)%position, grid%halfSmallestSubcell, source(nsource)%surface, 20)

             thisvel = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
                  thisOctal%rhow(subcell)/thisOctal%rho(subcell))
             e = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)
             thisOctal%rho(subcell) =  rhoThreshold
             thisOctal%rhou(subcell) = rhoThreshold * thisvel%x
             thisOctal%rhov(subcell) = rhoThreshold * thisvel%y
             thisOctal%rhow(subcell) = rhoThreshold * thisvel%z
             thisOctal%rhoe(subcell) = rhoThreshold * e

             temp(1) = 1.d0
             temp(2) = source(nsource)%position%x
             temp(3) = source(nsource)%position%y
             temp(4) = source(nsource)%position%z
             temp(5) = source(nsource)%mass 
             temp(6) = source(nsource)%velocity%x
             temp(7) = source(nsource)%velocity%y
             temp(8) = source(nsource)%velocity%z
             temp(9) = source(nsource)%radius
             temp(10) = source(nsource)%accretionRadius
             temp(11) = source(nsource)%angMomentum%x
             temp(12) = source(nsource)%angMomentum%y
             temp(13) = source(nsource)%angMomentum%z
             
             do iThread = 1, nHydroThreadsGlobal
                if (iThread /= myRankGlobal) then
                   call mpi_send(temp, 13, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                endif
             enddo
             
          endif

       endif
    enddo
  end subroutine recursaddSinks

  subroutine broadcastSinks
    use mpi
    integer :: ierr
     call MPI_BCAST(globalnSource, 1, MPI_INTEGER, 0, localWorldCommunicator, ierr)
     if (globalnSource > 0) then
        call MPI_BCAST(globalsourceArray(1:globalnSource)%position%x, globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%position%y, globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%position%z, globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%velocity%x, globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%velocity%y, globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%velocity%z, globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%mass     , globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%luminosity, globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%age       , globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%radius   , globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%mdot     , globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%accretionradius     , globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%angMomentum%x     , globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%angMomentum%y     , globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%angMomentum%z     , globalnSource, MPI_DOUBLE_PRECISION, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%stellar     , globalnSource, MPI_LOGICAL, 0, &
             localWorldCommunicator, ierr)
        call MPI_BCAST(globalsourceArray(1:globalnSource)%diffuse     , globalnSource, MPI_LOGICAL, 0, &
             localWorldCommunicator, ierr)
     endif
end subroutine broadcastSinks
  
  subroutine sendSinksToZerothThread(nSource, source)
    use mpi
    integer :: nSource, iSource, ierr
    type(SOURCETYPE) :: source(:)
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 32
    if (myrankGlobal == 0) then
 
     call mpi_recv(nSource, 1, MPI_INTEGER, 1, tag, localWorldCommunicator, status, ierr)
 
       do iSource = 1, nSource
          call mpi_recv(source(1:nSource)%position%x, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%position%y, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%position%z, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%mass,       nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%luminosity, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%age,        nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%velocity%x, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%velocity%y, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%velocity%z, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%radius,     nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%mdot,      nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%accretionRadius, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, &
               status, ierr)
          call mpi_recv(source(1:nSource)%angMomentum%x, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, &
               status, ierr)
          call mpi_recv(source(1:nSource)%angMomentum%y, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, &
               ierr)
          call mpi_recv(source(1:nSource)%angMomentum%z, nSource, MPI_DOUBLE_PRECISION, 1, tag, localWorldCommunicator, status, &
               ierr)
          call mpi_recv(source(1:nSource)%stellar, nSource, MPI_LOGICAL, 1, tag, localWorldCommunicator, status, ierr)
          call mpi_recv(source(1:nSource)%diffuse, nSource, MPI_LOGICAL, 1, tag, localWorldCommunicator, status, ierr)
       enddo
  
    endif

    if (myrankGlobal == 1) then

       call mpi_send(nSource, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, ierr)
  
       do iSource = 1, nSource
          call mpi_send(source(1:nSource)%position%x, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%position%y, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%position%z, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%mass,       nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%luminosity, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%age,        nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%velocity%x, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%velocity%y, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%velocity%z, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%radius,     nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%mdot,       nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%accretionRadius, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%angMomentum%x, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%angMomentum%y, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%angMomentum%z, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%stellar, nSource, MPI_LOGICAL, 0, tag, localWorldCommunicator, ierr)
          call mpi_send(source(1:nSource)%diffuse, nSource, MPI_LOGICAL, 0, tag, localWorldCommunicator, ierr)
       enddo
    endif
  end subroutine sendSinksToZerothThread

!  subroutine returnControlVolumeCells(grid, sendToThisThread, pos, radius)
!    use mpi
!    type(GRIDTYPE) :: grid
!    real(double) :: rho(100), cs(100), vol(100), phi(100), xpos(100), ypos(100), zpos(100), vx(100),vy(100),vz(100)
!    real(double) :: rhoFromThread(100), csFromThread(100), volFromThread(100), phiFromThread(100)
!    real(double) :: xposFromThread(100), yposFromThread(100), zposFromThread(100), vxFromThread(100),vyFromThread(100),vzFromThread(100)
!    integer :: ncells
!    integer :: sendToThisThread
!    type(VECTOR) :: pos
!    real(double) :: radius
!    integer :: nCellsFromThread
!    integer :: status(MPI_STATUS_SIZE)
!    integer :: i, ierr, j
!    integer, parameter :: tag = 55
!    ncells = 0
!
!    call recursiveGetControlVolumeCells(grid%octreeRoot, pos, radius, &
!         ncells, rho, cs, vol, phi, xpos, ypos, zpos, vx, vy, vz)
!    if (myrankGlobal == sendToThisThread) then
!       do i = 1, nHydroThreadsGlobal
!          if (myrankGlobal /= i) then
 !             call mpi_recv(nCellsFromThread, 1, MPI_INTEGER, i, tag, localWorldCommunicator, status, ierr)
!             if (nCellsFromThread > 0) then
!                call mpi_recv(rhoFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(csFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(volFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(phiFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(xposFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(yposFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(zposFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(vxFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(vyFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!                call mpi_recv(vzFromThread, ncellsFromThread, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, status, ierr)
!
!                do j = 1, nCellsFromThread
!                   nCells = nCells + 1
!                   rho(nCells) = rhoFromThread(j)
!                   cs(nCells) = rhoFromThread(j)
!                   vol(nCells) = rhoFromThread(j)
!                   rho(nCells) = rhoFromThread(j)
!                   rho(nCells) = rhoFromThread(j)
!                   rho(nCells) = rhoFromThread(j)
!                enddo
!             endif
!          endif
!       enddo
!    else
!       call mpi_send(ncells, 1, MPI_INTEGER, sendToThisThread, tag, localWorldCommunicator, ierr)
!       if (ncells > 0) then
!          call mpi_send(rho, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(cs, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(vol, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(phi, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(xpos, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(ypos, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(zpos, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(vx, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(vy, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!          call mpi_send(vz, ncells, MPI_DOUBLE_PRECISION, sendToThisThread, tag, localWorldCommunicator, ierr)
!       endif
!    endif
!  end subroutine returnControlVolumeCells
!
!  recursive subroutine recursiveGetControlVolumeCells(thisOctal, pos, radius, &
!       ncells, rho, cs, vol, phi, xpos, ypos, zpos, vx, vy, vz)
!    type(OCTAL), pointer :: thisOctal, child
!    type(VECTOR) :: pos, centre 
!    real(double) :: radius, r
!    integer :: ncells
!    real(double) :: rho(:), cs(:), vol(:), phi(:), xpos(:), ypos(:), zpos(:)
!    real(double) :: vx(:), vy(:), vz(:)
!    integer :: subcell, i
!
!    do subcell = 1, thisoctal%maxchildren
!       if (thisoctal%haschild(subcell)) then
!          ! find the child
!          do i = 1, thisoctal%nchildren, 1
!             if (thisoctal%indexchild(i) == subcell) then
!                child => thisoctal%child(i)
!                call recursiveGetControlVolumeCells(child,  pos, radius, &
!                     ncells, rho, cs, vol, phi, xpos, ypos, zpos, vx, vy, vz)
!
!                exit
!             end if
!          end do
!       else
!          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
!          centre = subcellCentre(thisOctal, subcell)
!          r = modulus(pos - centre)
!          if (r < radius) then
!             nCells = nCells + 1
!             rho(nCells) = thisOctal%rho(subcell)
!             cs(ncells) = soundSpeed(thisOctal, subcell)
!             vol(ncells) = cellVolume(thisOctal, subcell) * 1.d30
!             phi(ncells) = thisOctal%phi_i(subcell)
!             xPos(ncells) = centre%x
!             yPos(ncells) = centre%y
!             zPos(ncells) = centre%z
!             vx(ncells) = thisOctal%rhou(subcell)/thisOctal%rho(subcell)
!             vy(ncells) = thisOctal%rhov(subcell)/thisOctal%rho(subcell)
!             vz(ncells) = thisOctal%rhow(subcell)/thisOctal%rho(subcell)
!          endif
!         
!       endif
!    enddo
!  end subroutine recursiveGetControlVolumeCells
!
!
  function BondiHoyleRadius(source, thisOctal, subcell) result (rBH)

    real(double) :: rBH
    type(SOURCETYPE) :: source
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: vInfty, cInfty
    
    if (.not.inSubcell(thisOCtal, subcell, source%position)) then
       write(*,*) "bug in bondi-hoyle radius"
    endif
    
    vInfty = (thisOctal%rhou(subcell)/thisOctal%rho(subcell)-source%velocity%x)**2 + &
         (thisOctal%rhov(subcell)/thisOctal%rho(subcell)-source%velocity%y)**2 + &
         (thisOctal%rhow(subcell)/thisOctal%rho(subcell)-source%velocity%z)**2
    vInfty = sqrt(vInfty)
    cInfty = soundSpeed(thisOctal, subcell)
    rBH = bigG * source%mass / (vinfty**2 + cInfty**2)

  end function BondiHoyleRadius

  
  
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
    integer, parameter :: nsub = 8 
    
    cellCentre = subcellCentre(thisOctal, subcell)
    deltaX = thisOctal%subcellSize * gridDistanceScale
    n = 0 
    do i = 1, nsub
       do j = 1, nsub
          do k = 1, nsub
             x = cellcentre%x + (dble(i-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0
             y = cellcentre%y + (dble(j-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0
             z = cellcentre%z + (dble(k-1)+0.5d0) * thisOctal%subcellSize / dble(nSub) - thisOctal%subcellSize/2.d0

             pointPosition = (VECTOR(x, y, z) - source%position)*gridDistanceScale
             pointVelocity = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell) - source%velocity%x, &
                                    thisOctal%rhov(subcell)/thisOctal%rho(subcell) - source%velocity%y, &
                                    thisOctal%rhow(subcell)/thisOctal%rho(subcell) - source%velocity%z)
             jsp = modulus(pointPosition.cross.pointVelocity)
             esp = 0.5d0*modulus(pointVelocity)**2 - bigG*source%mass/modulus(pointPosition)
             if (esp > 0.d0) then
                rMin = 1.e30
             else
                rMin = real(-((bigG*source%mass)/(2.d0 * esp)) * &
                     (1.d0 - sqrt(1.d0 + (2.d0 * jsp * esp)/(bigG*source%mass)**2)))
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
    type(VECTOR) :: cellCentre, rVec, initialMomentum, finalMomentum, deltaMom
    real(double) :: u2, ekinetic, ethermal, eTot, deltaM

    cellCentre = subcellCentre(thisOctal, subcell)
    rVec = cellCentre - source%Position
    call normalize(rVec)


    cellMass = cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell)
    deltaM = thisOctal%etaline(subcell) * timestep

    initialMomentum = cellMass * &
         VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell)-source%velocity%x, &
                thisOctal%rhov(subcell)/thisOctal%rho(subcell)-source%velocity%y, &
                thisOctal%rhow(subcell)/thisOctal%rho(subcell)-source%velocity%z)
    radialMomentum = initialMomentum .dot. rVec


    finalMomentum = initialMomentum - radialMomentum * rVec

    radialMomentum = radialMomentum * (1.d0 - deltaM/cellMass)

    finalMomentum = finalMomentum + radialMomentum * rVec

    deltaMom = deltaMom + (finalMomentum - initialMomentum)
 
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

    thisOctal%rhou(subcell) = finalMomentum%x / (cellVolume(thisOctal, subcell)*1.d30)
    thisOctal%rhov(subcell) = finalMomentum%y / (cellVolume(thisOctal, subcell)*1.d30)
    thisOctal%rhow(subcell) = finalMomentum%z / (cellVolume(thisOctal, subcell)*1.d30)

    if (thisOctal%threed) then
       u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
    else
       u2 = (thisOctal%rhou(subcell)**2 +  thisOctal%rhow(subcell)**2)/thisOctal%rho(subcell)**2
    endif
    eKinetic = u2 / 2.d0



    Etot = eThermal + eKinetic

    thisOctal%rhoe(subcell) = eTot * thisOctal%rho(subcell)


  end subroutine correctMomentumOfGas


  subroutine correctMomentumOfSink(source, deltaMomentum)
    type(SOURCETYPE) :: source
    type(VECTOR) :: deltaMomentum, sourceMomentum

    sourceMomentum = source%mass * source%velocity
    sourceMomentum = sourceMomentum - deltaMomentum
    source%velocity = (1.d0/source%mass) * sourceMomentum
  end subroutine correctMomentumOfSink

  subroutine gatherSinks()
    use mpi
    integer :: iThread
    integer :: iSink, j
    integer, parameter :: tag = 34
    integer :: status(MPI_STATUS_SIZE)
    integer :: ierr

    iSink = 0 
    if (myrankGlobal == 0) then
       do iThread = 1, nHydroThreadsGlobal
          call MPI_RECV(j, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, status, ierr)
           iSink = iSink + j
        enddo
        write(*,*) "Total number of sinks over all threads ",iSink
     else
        call MPI_SEND(globalnSource, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, ierr)
     endif
    if (myrankGlobal == 0) then
       globalnSource = iSink
       call freeGlobalSourceArray()
       allocate(globalSourceArray(1:globalnSource))
       iSink = 1
       do iThread = 1, nHydroThreadsGlobal
          call MPI_RECV(j, 1, MPI_INTEGER, iThread, tag, &
               localWorldCommunicator, status, ierr)
          call MPI_RECV(globalSourceArray(iSink:iSink+j-1)%mass, j, MPI_DOUBLE_PRECISION, &
               iThread, tag, localWorldCommunicator, status, ierr)
          call MPI_RECV(globalSourceArray(iSink:iSink+j-1)%position%x, j, MPI_DOUBLE_PRECISION, &
               iThread, tag, localWorldCommunicator, status, ierr)
          call MPI_RECV(globalSourceArray(iSink:iSink+j-1)%position%y, j, MPI_DOUBLE_PRECISION, &
               iThread, tag, localWorldCommunicator, status, ierr)
          call MPI_RECV(globalSourceArray(iSink:iSink+j-1)%position%z, j, MPI_DOUBLE_PRECISION, &
               iThread, tag, localWorldCommunicator, status, ierr)
          call MPI_RECV(globalSourceArray(iSink:iSink+j-1)%velocity%x, j, MPI_DOUBLE_PRECISION, &
               iThread, tag, localWorldCommunicator, status, ierr)
          call MPI_RECV(globalSourceArray(iSink:iSink+j-1)%velocity%y, j, MPI_DOUBLE_PRECISION, &
               iThread, tag, localWorldCommunicator, status, ierr)
          call MPI_RECV(globalSourceArray(iSink:iSink+j-1)%velocity%z, j, MPI_DOUBLE_PRECISION, &
               iThread, tag, localWorldCommunicator, status, ierr)
          iSink = iSink + j
        enddo
     else
        call MPI_SEND(globalnSource, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, ierr)
        call MPI_SEND(globalSourceArray(1:globalnSource)%mass, globalnSource, &
             MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
        call MPI_SEND(globalSourceArray(1:globalnSource)%position%x, globalnSource, &
             MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
        call MPI_SEND(globalSourceArray(1:globalnSource)%position%y, globalnSource, &
             MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
        call MPI_SEND(globalSourceArray(1:globalnSource)%position%z, globalnSource, &
             MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
        call MPI_SEND(globalSourceArray(1:globalnSource)%velocity%x, globalnSource, &
             MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
        call MPI_SEND(globalSourceArray(1:globalnSource)%velocity%y, globalnSource, &
             MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
        call MPI_SEND(globalSourceArray(1:globalnSource)%velocity%z, globalnSource, &
             MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
     endif
     if (myrankGlobal /= 0) then
        call freeGlobalSourceArray()
     endif
     call MPI_BCAST(globalnSource, 1, MPI_INTEGER, 0, localWorldCommunicator, ierr)
     if (myrankGlobal /= 0) then
        allocate(globalSourceArray(1:globalnSource))
     endif
     call MPI_BCAST(globalSourceArray(1:globalnSource)%mass, globalnSource, &
          MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
     call MPI_BCAST(globalSourceArray(1:globalnSource)%position%x, globalnSource, &
          MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
     call MPI_BCAST(globalSourceArray(1:globalnSource)%position%y, globalnSource, &
          MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
     call MPI_BCAST(globalSourceArray(1:globalnSource)%position%z, globalnSource, &
          MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
     call MPI_BCAST(globalSourceArray(1:globalnSource)%velocity%x, globalnSource, &
          MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
     call MPI_BCAST(globalSourceArray(1:globalnSource)%velocity%y, globalnSource, &
          MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
     call MPI_BCAST(globalSourceArray(1:globalnSource)%velocity%z, globalnSource, &
          MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
        

   end subroutine gatherSinks

  subroutine doMyAccretion(grid, sourceArray, nSource, timeStep)
    use mpi
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: sourceArray(:)
    integer :: nSource
    real(double) :: timeStep
    integer :: iSource
    integer :: ierr
    real(double) :: tempa(3), suma(3)
    type(VECTOR) :: sourceMom
    real(double) :: temp
    real(double), allocatable :: accretedMass(:)
    type(VECTOR), allocatable :: accretedLinMomentum(:), accretedAngMomentum(:)

    if (writeoutput) write(*,*) "Labelling accreting cells"
    call labelAccretingCells(grid%octreeRoot, sourceArray, nSource)
    if (writeoutput) write(*,*) "Done."

    allocate(accretedMass(1:nSource), accretedLinMomentum(1:nSource), accretedAngMomentum(1:nSource))

    accretedMass = 0.d0
    accretedLinMomentum = VECTOR(0.d0, 0.d0, 0.d0)
    accretedAngMomentum = VECTOR(0.d0, 0.d0, 0.d0)
    

    if (writeoutput) write(*,*) "Performing gas accretion"
    call performGasAccretion(grid%octreeRoot, accretedMass, accretedLinMomentum, accretedAngMomentum, &
         sourceArray, nSource)
    if (writeoutput) write(*,*) "Done."

    do iSource = 1, nSource
!       if (writeoutput) write(*,*) "Accreting gas onto source ",isource

!       write(*,*) "mass accretion rate on rank ",myrankGlobal, " is ",accretedMass(isource)/timestep/msol/secstoyears

       call MPI_ALLREDUCE(accretedMass(iSource), temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
       accretedMass(iSource) = temp

       tempa(1) = accretedLinMomentum(iSource)%x
       tempa(2) = accretedLinMomentum(iSource)%y
       tempa(3) = accretedLinMomentum(iSource)%z
    
       call MPI_ALLREDUCE(tempa, suma, 3, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)

       accretedLinMomentum(iSource)%x = suma(1)
       accretedLinMomentum(iSource)%y = suma(2)
       accretedLinMomentum(iSource)%z = suma(3)

       tempa(1) = accretedAngMomentum(iSource)%x
       tempa(2) = accretedAngMomentum(iSource)%y
       tempa(3) = accretedAngMomentum(iSource)%z
    
       call MPI_ALLREDUCE(tempa, suma, 3, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)

       accretedAngMomentum(iSource)%x = suma(1)
       accretedAngMomentum(iSource)%y = suma(2)
       accretedAngMomentum(iSource)%z = suma(3)


!       if (writeoutput) then
!
!          write(*,*) "source ",isource
!          write(*,*) "before vel ",sourceArray(isource)%velocity/1.d5
!          write(*,*) "before mass ",sourceArray(isource)%mass/msol
!
!          write(*,*) "accreted mass ", accretedMass(isource)
!          write(*,*) "accreted lin mom ", accretedLinMomentum(isource)
!          write(*,*) "accreted ang mom ", accretedAngmomentum(isource)
!       endif

       sourceMom = sourceArray(iSource)%mass * sourceArray(iSource)%velocity
       sourceMom = sourceMom + accretedLinMomentum(iSource)

       sourceArray(iSource)%mass = sourceArray(iSource)%mass + accretedMass(iSource)

       sourceArray(isource)%velocity = sourceMom / sourceArray(iSource)%mass

       sourceArray(iSource)%angMomentum = sourceArray(iSource)%angMomentum + accretedAngMomentum(iSource)

       sourceArray(iSource)%mdot = accretedMass(isource)/timestep
       if (myrankWorldGlobal == 1) then
       if (accretedMass(iSource) > 0.d0) write(*,*) "Accretion rate for source ",isource, ": ", &
            (accretedMass(isource)/timestep)/msol * (365.25d0*24.d0*3600.d0)
!          write(*,*)  "position ",sourceArray(isource)%position
!          write(*,*) "velocity ",sourceArray(isource)%velocity/1.d5
!          write(*,*) "mass (solar) ",sourceArray(isource)%mass/msol
       endif
    enddo
    deallocate(accretedMass, accretedLinMomentum, accretedAngMomentum)

  end subroutine doMyAccretion
recursive subroutine checkSetsAreTheSame(thisOctal)
  use mpi
  use ion_mod, only : nGlobalIon, globalIonArray, returnMu

  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i, ierr
  real(double), allocatable :: temp(:), temp2(:)
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call checkSetsAreTheSame(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          allocate(temp(1:nHydroSetsGlobal),temp2(1:nHydroSetsGlobal))
          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%temperature(subcell))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid temperatures differ"

          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%phi_i(subcell))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid phi_i differ ",temp2

          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%ionfrac(subcell,1))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid ionfrac1 differ ",temp2

          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%ionfrac(subcell,2))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid ionfrac2 differ ",temp2



          temp = 0.d0
          if (photoionPhysics) then
             temp(myHydroSetGlobal+1) = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
          else
             temp(myHydroSetGlobal+1)= 2.33d0
          endif
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid mu differ ",temp2

          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%pressure_i(subcell))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid pressure_i differ ",temp2

          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%rho(subcell))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid densities differ"

          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%rhou(subcell))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid rhou differ ",temp2

          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%rhov(subcell))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid rhov differ ",temp2

          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%rhow(subcell))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid rhow differ ",temp2


          deallocate(temp, temp2)
             
       endif
    enddo
  end subroutine checkSetsAreTheSame



  recursive subroutine labelAccretingCells(thisOctal, source, nsource)
    use mpi
    type(SOURCETYPE) :: source(:)
    type(VECTOR) :: cellCentre
    integer :: nSource
    type(VECTOR) :: cellVelocity
    real(double) :: cellMass, eGrav, Ekin, eMin
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    real(double) :: r , rv
    integer :: subcell, i
    integer :: nSourcesInAccretionRadius
    integer, allocatable :: iSource(:)
    integer :: thisSource

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call labelAccretingCells(child, source, nsource)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

          if (.not.associated(thisOctal%chiline)) then
             allocate(thisOctal%chiLine(1:thisOctal%maxChildren))
             thisOctal%chiLine = 0.d0
          endif

          cellCentre = subcellCentre(thisOctal, subcell)
          cellMass = cellVolume(thisOctal, subcell) * thisOctal%rho(subcell) * 1.d30

          if (thisOctal%threeD) then
             cellVelocity = VECTOR(thisOctal%rhou(subcell) / thisOctal%rho(subcell), &
                  thisOctal%rhov(subcell) / thisOctal%rho(subcell), &
                  thisOctal%rhow(subcell) / thisOctal%rho(subcell))
          else
             r = cellCentre%x
             rv = thisOctal%rhov(subcell) / (r * thisOctal%rho(subcell))
             cellVelocity = VECTOR(thisOctal%rhou(subcell) / thisOctal%rho(subcell), &
                  rv, &
                  thisOctal%rhow(subcell) / thisOctal%rho(subcell))
          endif

          nSourcesInAccretionRadius = 0
          allocate(iSource(1:nSource))
          do i = 1, nSource
             if (modulus(source(i)%position - cellCentre)*1.d10 < source(i)%accretionRadius) then
                nSourcesInAccretionRadius = nSourcesInAccretionRadius + 1
                iSource(nSourcesInAccretionRadius) = i
             endif
          enddo
          thisOctal%chiLine(subcell) = 0.d0
          if (nSourcesInAccretionRadius >= 1) then

             if (nSourcesInAccretionRadius == 1) then
                thisSource = iSource(1)
             else
                eMin = 1.d30
                thisSource = 0
                do i = 1, nSourcesInAccretionRadius
                   r = modulus(source(iSource(i))%position-cellCentre)*1.d10
                   eGrav = -bigG*cellMass*source(isource(i))%mass / r
                   eKin = 0.5d0 * cellMass * modulus(cellVelocity)**2
                   if ((eGrav + eKin) < eMin) then
                      thisSource = isource(i)
                      eMin = eGrav + eKin
                   endif
                enddo
             endif
             
             thisOctal%chiLine(subcell) = thisSource
          endif
          deallocate(iSource)
       endif
    enddo
  end subroutine labelAccretingCells

  recursive subroutine performGasAccretion(thisOctal, accretedMass, accretedLinMomentum, accretedAngMomentum, &
       source, nSource)
    use mpi
    use inputs_mod, only : rhoThreshold, geometry, smallestCellSize
    type(SOURCETYPE) :: source(:)
    real(double) :: accretedMass(:), thisCellVolume
    type(VECTOR) :: accretedLinMomentum(:), accretedAngMomentum(:), cellCentre, cellVelocity
    real(double) :: eGrav, eThermal, eKinetic, cellMass, rhoLocal, localAccretedMass, r, rv, localAccretedAngMom
    type(VECTOR) :: gasMom, localAccretedMom, localAngMom
    integer :: nSource
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i, iSource
    logical :: onAxis, cellBound
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call performGasAccretion(child, accretedMass, accretedLinMomentum, accretedAngMomentum, source, nSource)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          if (thisOctal%ghostCell(subcell)) cycle

          iSource = nint(thisOctal%chiLine(subcell))
          if (iSource > 0) then
             rhoLocal = thisOctal%rho(subcell)
             cellCentre = subcellCentre(thisOctal, subcell)
             
             thisCellVolume = cellVolume(thisOctal, subcell) * 1.d30
             cellMass = thisOctal%rho(subcell) * thisCellVolume






             onAxis = .false.
             if (cylindricalHydro) then
                onAxis = .not.(cellCentre%x - thisOctal%subcellSize/2.d0+0.1d0*smallestCellSize < 0.d0)
             endif

             if (thisOctal%threeD.or.spherical) then
                cellVelocity = VECTOR(thisOctal%rhou(subcell) / thisOctal%rho(subcell), &
                     thisOctal%rhov(subcell) / thisOctal%rho(subcell), &
                     thisOctal%rhow(subcell) / thisOctal%rho(subcell))
             else if (thisOctal%twoD) then
                r = cellCentre%x*1.d10
                rv = thisOctal%rhov(subcell) / (r * thisOctal%rho(subcell))
                cellVelocity = VECTOR(thisOctal%rhou(subcell) / thisOctal%rho(subcell), &
                     rv, &
                     thisOctal%rhow(subcell) / thisOctal%rho(subcell))
                cellVelocity = VECTOR(thisOctal%rhou(subcell) / thisOctal%rho(subcell), &
                     0.d0, &
                     thisOctal%rhow(subcell) / thisOctal%rho(subcell))
             else
                r = cellCentre%x*1.d10
                cellVelocity = VECTOR(thisOctal%rhou(subcell) / thisOctal%rho(subcell), &
                     0.d0, 0.d0)
             endif


             if (.not.(cylindricalHydro.or.spherical)) then
                eGrav = cellMass * thisOctal%phi_i(subcell)
             else
                 call calculatePotentialFromSinks(thisOctal, subcell, source, nSource, eGrav)
                 eGrav = eGrav + cellMass * thisOctal%phi_i(subcell)
              endif
             eThermal = 0.5d0 * cellMass * soundSpeed(thisOctal, subcell)**2
             if (geometry == "bondi") ethermal = 0.d0
             ekinetic = 0.5d0 * cellMass * modulus(cellVelocity-source(isource)%velocity)**2

             cellBound = .not.((eKinetic + eThermal + eGrav > 0.d0).and.(rhoLocal > rhoThreshold))
             if (cylindricalHydro.or.spherical) CellBound = .true.
             if (.not.cellBound) then
                write(*,*) "Cell in accretion radius but not bound ",eKinetic+eGrav+eThermal
                write(*,*) "eGrav ",eGrav
                write(*,*) "ethermal ",eThermal
                write(*,*) "eKinetic ",eKinetic
             endif
             localAccretedMass = 0.d0
             localAccretedMom = VECTOR(0.d0, 0.d0, 0.d0)
             localAngMom = VECTOR(0.d0, 0.d0, 0.d0)

             if (cylindricalHydro) then 
               cellBound = .true.
                localAccretedAngMom = thiscellVolume * thisOctal%rhov(subcell) 
                accretedAngMomentum(isource) = accretedAngMomentum(isource) + VECTOR(0.d0, 0.d0, localAccretedAngMom)
!                localAccretedMom = thiscellVolume * VECTOR(thisOctal%rhou(subcell), 0.d0, thisOctal%rhow(subcell))
!                accretedLinMomentum(isource) = accretedLinMomentum(isource) + localAccretedMom
!                thisOctal%rhou(subcell) = 0.d0
                thisOctal%rhov(subcell) = 0.d0
!                thisOctal%rhow(subcell) = 0.d0
             endif

             if (cellBound) then
                if (rhoLocal > rhoThreshold) then
                   
                   localaccretedMass = (rhoLocal - rhoThreshold) * cellVolume(thisOctal,Subcell)*1.d30
!                   write(*,*) myrankGlobal, " localaccretedmass ",localaccretedmass
                   accretedMass(isource) = accretedMass(isource) + localAccretedMass

                   if (thisOctal%threed.or.spherical) then
                      gasMom = (cellVolume(thisOctal, subcell) * 1.d30)  * &
                           VECTOR(thisOctal%rhou(subcell), thisOctal%rhov(subcell), thisOctal%rhow(subcell))
                   else if (spherical) then
                      gasMom = (cellVolume(thisOctal, subcell) * 1.d30)  * &
                           VECTOR(thisOctal%rhou(subcell), 0.d0, 0.d0)
                   else
                      gasMom = (cellVolume(thisOctal, subcell) * 1.d30)  * &
                           VECTOR(thisOctal%rhou(subcell), 0.d0, thisOctal%rhow(subcell))
                   endif

                   if (thisOctal%threed.or.spherical) then
                      localAccretedMom = localaccretedMass * VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                           thisOctal%rhov(subcell)/thisOctal%rho(subcell), thisOctal%rhow(subcell)/thisOctal%rho(subcell))
                   else if (spherical) then
                      localAccretedMom = localaccretedMass * VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), 0.d0, 0.d0)
                   else
                      if (onAxis) then
                         localAccretedMom = VECTOR(thisOctal%rhou(subcell), &
                              0.d0, localAccretedMass * thisOctal%rhow(subcell)/thisOctal%rho(subcell))
                      else
                         localAccretedMom = localaccretedMass * VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                              0.d0, thisOctal%rhow(subcell)/thisOctal%rho(subcell))
                      endif
                   endif

                   if (thisOctal%threed.or.spherical) then
                      localAngMom =  ((cellCentre - source(isource)%position)*1.d10) .cross. localAccretedMom
                      localAccretedAngMom = modulus(localAngMom)
                   endif

                   accretedLinMomentum(isource) = accretedLinMomentum(isource) + localAccretedMom
                   if (thisOctal%threed) then
                         accretedAngMomentum(isource) = accretedAngMomentum(isource) + VECTOR(0.d0, 0.d0, localAccretedAngMom)
                   endif
                   if (thisOctal%threed.or.spherical) then
                      gasMom = gasMom - localaccretedMom
                      thisOctal%rho(subcell) = rhoThreshold
                      cellMass = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell) * 1.d30
                      cellVelocity = gasMom / cellMass
                      thisOctal%rhou(subcell) =  thisOctal%rho(subcell) * cellVelocity%x
                      thisOctal%rhov(subcell) =  thisOctal%rho(subcell) * cellVelocity%y
                      thisOctal%rhow(subcell) =  thisOctal%rho(subcell) * cellVelocity%z
                   else if (spherical) then
                      gasMom = gasMom - localaccretedMom
                      thisOctal%rho(subcell) = rhoThreshold
                      cellMass = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell) * 1.d30
                      cellVelocity = gasMom / cellMass
                      thisOctal%rhou(subcell) =  thisOctal%rho(subcell) * cellVelocity%x
                   else
                      gasMom = gasMom - localaccretedMom
                      thisOctal%rho(subcell) = rhoThreshold
                      cellMass = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell) * 1.d30
                      cellVelocity = gasMom / cellMass
                      thisOctal%rhou(subcell) =  thisOctal%rho(subcell) * cellVelocity%x
                      thisOctal%rhow(subcell) =  thisOctal%rho(subcell) * cellVelocity%z
                      if (onAxis) thisOctal%rhou(subcell) = 0.
                      thisOctal%rhov(subcell) = 0.
                      if (thisOctal%rhov(subcell) < 0.d0) then
                         write(*,*) "warning negative rho v. vel = ",(1.d-5)*thisOctal%rhov(subcell)/thisOctal%rho(subcell)
                         write(*,*) "rho ",thisOctal%rho(subcell)
                      endif
                   endif
                   thisOctal%rhoe(subcell) = thisOctal%rhoe(subcell) * (rhoLocal-rhoThreshold)/rhoLocal
                endif
             endif
          endif
       endif
    enddo
  end subroutine performGasAccretion


  recursive subroutine correctMomenta(thisOctal, source, timestep,deltaMom)
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

  subroutine addStellarWind(grid, nSource, sourceArray, dt)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: nSource, i, ierr
    real(double) :: dt, totalVolume, temp
    type(SOURCETYPE) :: sourceArray(:)

    do i = 1, nSource
          totalVolume = 0.d0
          call findStellarWindVolume(grid%octreeRoot, sourceArray(i), totalVolume)
          call MPI_ALLREDUCE(totalVolume, temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
          totalVolume = temp
          call addStellarWindRecur(grid%octreeRoot, sourceArray(i), dt, totalVolume)
    enddo
  end subroutine addStellarWind

  subroutine addSupernovae(grid, nSource, sourceArray, dt, sneAdded)
    use mpi
    use starburst_mod, only : checkSourceSupernova, removeSource
    type(GRIDTYPE) :: grid
    integer :: nSource, i, ierr
    real(double) :: dt, totalVolume, temp
    logical :: sneAdded
    type(SOURCETYPE) :: sourceArray(:)
    integer :: nSupernova, supernovaIndex(10)
    real(double) :: ejectaMass(10), ke(10)

    sneAdded = .false.
    call checkSourceSupernova(nSource, sourceArray, nSupernova, supernovaIndex, ejectaMass, ke)

    if (nSupernova > 0) then
       sneAdded = .true.
       do i = 1, nSupernova
          totalVolume = 0.d0
          call findStellarWindVolume(grid%octreeRoot, sourceArray(supernovaIndex(i)), totalVolume)
          call MPI_ALLREDUCE(totalVolume, temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
          totalVolume = temp
          call addSupernovaRecur(grid%octreeRoot, sourceArray(supernovaIndex(i)), dt, totalVolume, ejectaMass(i),ke(i))
       enddo

!   need to sort supernovaindex into descending order at this point       

       do i = 1, nSupernova
          call removeSource(sourceArray, nSource, supernovaIndex(i))
       enddo


    endif

  end subroutine addSupernovae



  recursive subroutine addStellarWindRecur(thisOctal, source, dt, totalVolume)
    use mpi
    use inputs_mod, only : smallestCellSize
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    type(SOURCETYPE) :: source
    integer :: subcell, i
    real(double) :: r, v, vterm, maxSpeed, mdot, dt, totalVolume
    real(double) :: totalMassThisInterval, thisRho, thisMom
    type(VECTOR) :: cellCentre, rVec
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call addStellarWindRecur(child, source, dt, totalVolume)
                exit
             end if
          end do
       else

          
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          maxspeed = 0.2d5
          thisOctal%boundaryCell(subcell) = .false.
          cellCentre = subcellCentre(thisOctal, subcell)
          rVec = cellCentre - source%position
          r = modulus(rVec)
          rVec = rVec / r

          mdot = source%mDotWind

!           write(*,*) "adding stellar wind at rate ",(mdot/msol)*365.25d0*24.d0*3600.d0
          vterm = 2000.d0 * 1.d5

          totalMassThisInterval = mDot * dt

          if (r < smallestCellSize*5.d0) then
             thisOctal%boundaryCell(subcell) = .true.
             v = cellVolume(thisOctal, subcell) * 1.d30 
             thisRho = totalMassThisInterval / totalVolume
             thisMom = thisRho * vterm
             thisOctal%rho(subcell) = thisOctal%rho(subcell) + thisRho

             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + rVec%x * thisMom
             thisOctal%rhov(subcell) = thisOctal%rhov(subcell) + rVec%y * thisMom
             thisOctal%rhow(subcell) = thisOctal%rhow(subcell) + rVec%z * thisMom
          endif

       endif
    enddo
  end subroutine addStellarWindRecur

  recursive subroutine addSuperNovaRecur(thisOctal, source, dt, totalVolume, ejectaMass, ke)
    use mpi
    use inputs_mod, only : smallestCellSize
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    type(SOURCETYPE) :: source
    integer :: subcell, i
    real(double) :: r, v, vterm, maxSpeed, mdot, dt, totalVolume
    real(double) :: totalMassThisInterval, thisRho, thisMom, ke, ejectaMass
    type(VECTOR) :: cellCentre, rVec
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call addSuperNovaRecur(child, source, dt, totalVolume, ejectaMass, ke)
                exit
             end if
          end do
       else

          
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          maxspeed = 0.2d5
          thisOctal%boundaryCell(subcell) = .false.
          cellCentre = subcellCentre(thisOctal, subcell)
          rVec = cellCentre - source%position
          r = modulus(rVec)
          rVec = rVec / r

          mdot = 1.d-6 * mSol / (365.25d0 * 24.d0 * 3600.d0)
          vterm = 0.2d0 * 1.d5

          totalMassThisInterval = ejectaMass

          vterm = sqrt(2.d0 * ke / totalMassThisInterval)

          if (r < smallestCellSize*5.d0) then
             thisOctal%boundaryCell(subcell) = .true.
             v = cellVolume(thisOctal, subcell) * 1.d30
             thisRho = totalMassThisInterval / totalVolume
             thisMom = thisRho * vterm
             thisOctal%rho(subcell) = thisOctal%rho(subcell) + thisRho

             thisOctal%rhou(subcell) = thisOctal%rhou(subcell) + rVec%x * thisMom
             thisOctal%rhov(subcell) = thisOctal%rhov(subcell) + rVec%y * thisMom
             thisOctal%rhow(subcell) = thisOctal%rhow(subcell) + rVec%z * thisMom
          endif

       endif
    enddo
  end subroutine addSuperNovaRecur

  recursive subroutine findStellarWindVolume(thisOctal, source, totalVolume)
    use mpi
    use inputs_mod, only : smallestCellSize
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    type(SOURCETYPE) :: source
    integer :: subcell, i
    real(double) :: r, totalVolume
    type(VECTOR) :: cellCentre, rVec
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call findStellarWindVolume(child, source, totalVolume)
                exit
             end if
          end do
       else

          
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          cellCentre = subcellCentre(thisOctal, subcell)
          rVec = cellCentre - source%position
          r = modulus(rVec)
          rVec = rVec / r
          if (r < smallestCellSize*5.d0) then
             totalVolume = totalVolume + cellVolume(thisOctal, subcell)*1.d30
          endif

       endif
    enddo
  end subroutine findStellarWindVolume


  subroutine getPointsInAccretionRadius(thisOctal, subcell, radius, grid, npoints, position, vel, mass, phi, cs)
    use mpi
    real(double) :: radius
    type(GRIDTYPE) :: grid
    integer :: nPoints
    type(VECTOR) :: position(:), vel(:), rVec
    real(double) :: mass(:), phi(:), cs(:), r
    type(OCTAL), pointer :: thisOctal
    real(double) :: temp(13)
    integer :: iThread
    integer :: subcell
    integer :: ierr
    integer, parameter :: tag = 65
    integer :: nOther
    integer :: status(MPI_STATUS_SIZE)

    npoints = 1
    position(1) = subcellCentre(thisOctal, subcell)
    if (.not.cylindricalHydro) then
       vel(1) = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
            thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
            thisOctal%rhow(subcell)/thisOctal%rho(subcell))
    else
       rVec = subcellCentre(thisOctal, subcell)
       r = sqrt(rVec%x**2 + rVec%y**2) * gridDistanceScale
       vel(1) = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
            thisOctal%rhov(subcell)/(thisOctal%rho(subcell)*r), &
            thisOctal%rhow(subcell)/thisOctal%rho(subcell))
    endif
    mass(1) = cellVolume(thisOctal, subcell)*1.d30*thisOctal%rho(subcell)
    phi(1) = thisOctal%phi_i(subcell)
    cs(1) = soundSpeed(thisOctal, subcell)
    call getPointsInAccretionRadiusLocal(grid%octreeRoot, position(1), radius, grid, npoints, position, vel, mass, phi, cs)
    
    temp(1) = 2.d0
    temp(2) = position(1)%x
    temp(3) = position(1)%y
    temp(4) = position(1)%z
    temp(5) = radius
    temp(6) = thisOctal%mpiThread(subcell)
    do iThread = 1, nHydroThreadsGlobal
       if (iThread /= myRankGlobal) then
          call MPI_SEND(temp, 13, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
          call MPI_RECV(nOther, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, status, ierr)
          if (nOther > 0) then
             call MPI_RECV(position(nPoints+1:nPoints+nOther)%x, nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)
             call MPI_RECV(position(nPoints+1:nPoints+nOther)%y, nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)
             call MPI_RECV(position(nPoints+1:nPoints+nOther)%z, nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)
             
             call MPI_RECV(vel(nPoints+1:nPoints+nOther)%x, nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)
             call MPI_RECV(vel(nPoints+1:nPoints+nOther)%y, nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)
             call MPI_RECV(vel(nPoints+1:nPoints+nOther)%z, nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)
             
             call MPI_RECV(mass(nPoints+1:nPoints+nOther), nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)
             
             call MPI_RECV(phi(nPoints+1:nPoints+nOther), nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)
             
             call MPI_RECV(cs(nPoints+1:nPoints+nOther), nOther, MPI_DOUBLE_PRECISION, &
                  iThread, tag, localWorldCommunicator, status, ierr)

             nPoints = nPoints + nOther

          endif
       endif
    enddo

  end subroutine getPointsInAccretionRadius
  
             
  recursive subroutine getPointsinAccretionRadiusLocal(thisOctal, thisPoint, radius, grid, nPoints, position, vel, mass, phi, cs)
    type(VECTOR) :: thisPoint
    real(double) :: radius
    type(GRIDTYPE) :: grid
    integer :: nPoints
    type(VECTOR) :: position(:), vel(:)
    real(double) :: mass(:), phi(:), cs(:), r
    type(OCTAL), pointer :: thisOctal, child
    type(VECTOR) :: cen
    integer :: subcell, i
    
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call getPointsinAccretionRadiusLocal(child, thisPoint, radius, grid, nPoints, position, vel, mass, phi, cs)
                exit
             end if
          end do
       else

          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (thisOctal%ghostCell(subcell)) cycle

          if (inSubcell(thisOctal, subcell, thisPoint)) cycle

          cen = subcellCentre(thisOctal, subcell)
          if (modulus(cen - thisPoint) < radius) then
             nPoints = npoints + 1
             position(nPoints) = subcellCentre(thisOctal, subcell)
             if (.not.cylindricalHydro) then
                vel(nPoints) = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                     thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
                     thisOctal%rhow(subcell)/thisOctal%rho(subcell))
             else
                r = sqrt(cen%x**2+cen%y**2)*gridDistanceScale
                vel(nPoints) = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                     thisOctal%rhov(subcell)/(thisOctal%rho(subcell)*r), &
                     thisOctal%rhow(subcell)/thisOctal%rho(subcell))
             endif
             mass(nPoints) = cellVolume(thisOctal, subcell)*1.d30*thisOctal%rho(subcell)
             phi(nPoints) = thisOctal%phi_i(subcell)
             cs(nPoints) = soundSpeed(thisOctal, subcell)
          endif

       endif
    enddo
  end subroutine getPointsinAccretionRadiusLocal

  subroutine shutdownRadiusServer()
    use mpi
    real(double) :: loc(13)
    integer :: iThread
    integer :: ierr
    integer, parameter :: tag = 65

    do iThread = 1, nHydroThreadsGlobal
       if (myRankGlobal /= iThread) then
          loc(1) = 3.d0
          call MPI_SEND(loc, 13, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine shutdownRadiusServer
          
    

  subroutine getPointsInAccretionRadiusServer(grid, nSource, source)
    use mpi
    integer :: nSource
    type(sourcetype) :: source(:)
    real(double) :: temp(13)
    type(GRIDTYPE) :: grid
    integer :: nPoints
    integer :: ierr
    integer, parameter :: tag = 65
    integer :: status(MPI_STATUS_SIZE)
    logical :: stillServing
    type(VECTOR) :: thisPoint
    integer :: iThread
    integer, parameter :: maxPoints = 100
    type(VECTOR) :: position(maxPoints), vel(maxPoints)
    real(double) :: mass(maxPoints), phi(maxPoints), cs(maxPoints)
    integer :: isignal
    real(double) :: radius

    stillServing = .true.

    do while(stillServing) 
       call MPI_RECV(temp, 13, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, localWorldCommunicator, status, ierr)

       iSignal = int(temp(1))
       select case(isignal)
          case(1)  ! added new source
             nSource = nSource + 1
             source(nsource)%position%x = temp(2) 
             source(nsource)%position%y = temp(3) 
             source(nsource)%position%z = temp(4) 
             source(nsource)%mass  = temp(5) 
             source(nsource)%velocity%x = temp(6) 
             source(nsource)%velocity%y = temp(7) 
             source(nsource)%velocity%z = temp(8) 
             source(nsource)%radius = temp(9)   
             source(nSource)%stellar = .true.
             source(nSource)%diffuse = .false.
             source(nSource)%accretionRadius = temp(10)
             source(nSource)%angMomentum%x = temp(11)
             source(nSource)%angMomentum%y = temp(12)
             source(nSource)%angMomentum%z = temp(13)
             call fillSpectrumBB(source(nsource)%spectrum, 1000.d0, &
                  100.d0, 1.d7, 1000)
             call buildSphereNbody(source(nsource)%position, grid%halfSmallestSubcell, source(nsource)%surface, 20)
          case(2) ! want points in radius

             thisPoint%x = temp(2)
             thisPoint%y = temp(3)
             thisPoint%z = temp(4)
             radius = temp(5)
             iThread = int(temp(6))
             nPoints = 0
             call getPointsinAccretionRadiusLocal(grid%octreeRoot, thisPoint, radius, grid, nPoints, position, vel, mass, phi, cs)
             call MPI_SEND(nPoints, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
             if (nPoints > 0) then
                call MPI_SEND(position(1:nPoints)%x, nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                call MPI_SEND(position(1:nPoints)%y, nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                call MPI_SEND(position(1:nPoints)%z, nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                
                call MPI_SEND(vel(1:nPoints)%x, nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                call MPI_SEND(vel(1:nPoints)%y, nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                call MPI_SEND(vel(1:nPoints)%z, nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                
                call MPI_SEND(mass(1:nPoints), nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                
                call MPI_SEND(phi(1:nPoints), nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                
                call MPI_SEND(cs(1:nPoints), nPoints, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
             endif
          case(3) ! end server
             stillServing = .false.
          case DEFAULT
             write(*,*) "Unknown signal received in getpointsinaccretionradiusserver ",isignal
       end select
    end do
  end subroutine getPointsInAccretionRadiusServer









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


  recursive subroutine imposeTstructure(thisOctal)
    use inputs_mod, only : smallestCellSize
    real(double) :: rMod
    TYPE(OCTAL),pointer :: thisOctal
    TYPE(OCTAL),pointer :: child
    integer :: i, subcell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call imposeTstructure(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

          rMod = modulus(subcellCentre(thisOctal, subcell))
          thisOctal%temperature(subcell) = max(20.0,5000*real(sqrt(smallestCellSize/rMod)))



       end if
    end do
  end subroutine imposeTstructure


  recursive subroutine calculateTemperatureFromEnergy(thisOctal)
    use inputs_mod, only : photoionPhysics, honly, simplemu
    use ion_mod, only : nGlobalIon, globalIonArray, returnMu, returnMusimple
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: eThermal
    real(double) :: mu


    do subcell = 1, thisOctal%maxChildren
       if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateTemperatureFromEnergy(child)
                exit
             end if
          end do
       else
            if (photoionPhysics) then
               mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
               if(hOnly .and. simpleMu) then
                  mu = returnMuSimple(thisOctal, subcell)
               end if
            else
               mu = 2.33d0
            endif

          eThermal = thisOctal%rhoe(subcell)/thisOctal%rho(subcell)

          thisOctal%temperature(subcell) = real(mu * mHydrogen * eThermal / (1.5d0 * kerg))
          
          if (thisOctal%temperature(subcell) < 0.d0) write(*,*) "temp less than 0", ethermal, thisOctal%rhoe(subcell), mu, &
               subcellCentre(thisOctal,subcell)
       endif
    enddo
  end subroutine calculateTemperatureFromEnergy

  recursive subroutine setEquationOfState(thisOctal, eos)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    integer :: eos
 

    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setEquationOfState(child, eos)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle

          thisOctal%iEquationOfState(subcell) = eos
          if (myrankglobal == 1) write(*,*) subcell, " eos ", thisOctal%iEquationOfState(subcell)
       endif
    enddo
  end subroutine setEquationOfState


#endif

end module hydrodynamics_mod

#endif
