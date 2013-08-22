! this is a module to apply the diffusion approximation to
! the very optically thick central regions of circumstellar stellar discs.


module diffusion_mod

use constants_mod
use vector_mod
use messages_mod
use gridtype_mod, only: GRIDTYPE
use vtk_mod, only: writeVTKfile
use octal_mod, only: OCTAL, OCTALWRAPPER, subcellCentre, returndPhi
use amr_mod, only: returnKappa, tauAlongPath, inOctal, amrGridValues, &	
     countVoxels, getOctalArray, octalOnThread
implicit none


contains



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
          thisOctal%oldtemperature(subcell) = thisOctal%temperature(subcell)
          if (thisOctal%diffusionApprox(subcell).and.(.not.thisOctal%fixedTemperature(subcell))) then
             thisOctal%temperature(subcell) = real(max(thisOctal%chiline(subcell),3.d0))
          endif
       endif
    enddo
  end subroutine copychilinetoTemperature

  recursive subroutine countDiffusionCells(thisOctal, ndiff)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i, ndiff
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call countDiffusionCells(child, ndiff)
                exit
             end if
          end do
       else
          if (thisOctal%diffusionApprox(subcell).and.(.not.thisOctal%fixedTemperature(subcell))) then
             ndiff = ndiff + 1
          endif
       endif
    enddo
  end subroutine countDiffusionCells

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






  recursive subroutine setDiffusionCoeff(grid, thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(GRIDTYPE) :: grid
    real(double) :: kros
    integer :: subcell, i
    kros = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setDiffusionCoeff(grid, child)
                exit
             end if
          end do
       else
          thisOctal%temperature(subcell) = real((thisOctal%eDens(subcell)/arad)**0.25d0)
          call returnKappa(grid, thisOctal, subcell, rosselandKappa=kros)
          if (.not.associated(thisOctal%kappaRoss)) allocate(thisOctal%kappaRoss(1:thisOctal%maxChildren))
          if (.not.associated(thisOctal%diffusionCoeff)) allocate(thisOctal%diffusionCoeff(1:thisOctal%maxChildren))
          thisOctal%kappaRoss(subcell) = kRos
          thisOctal%diffusionCoeff(subcell) =  cSpeed / max(1.d-20,(kRos * thisOctal%rho(subcell)))
          thisOctal%oldeDens(subcell) = thisOctal%eDens(subcell)
       endif
    enddo
  end subroutine setDiffusionCoeff

  recursive subroutine seteDens(grid, thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(GRIDTYPE) :: grid
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call seteDens(grid, child)
                exit
             end if
          end do
       else
          if (.not.associated(thisOctal%eDens)) allocate(thisOctal%eDens(1:thisOctal%maxChildren))
          if (.not.associated(thisOctal%oldeDens)) allocate(thisOctal%oldeDens(1:thisOctal%maxChildren))
          thisOctal%eDens(subcell) = aRad * thisOctal%temperature(subcell)**4
       endif
    enddo
  end subroutine seteDens

  recursive subroutine checkConvergence(thisOctal, tol, dtmax, converged)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real :: deltaT, tol, dtmax
    integer :: subcell, i
    logical :: converged
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call checkConvergence(child, tol, dtMax, converged)
                exit
             end if
          end do
       else
          deltaT = real(abs(thisOctal%eDens(subcell)-thisOctal%oldeDens(subcell)) &
                / thisOctal%oldEdens(subcell))
          dtMax = max(dtMax, deltaT)
          if (deltaT > tol) converged = .false.
          
       endif
    enddo
  end subroutine checkConvergence

  subroutine gaussSeidelSweep(grid,  tol, demax, converged)
#ifdef MPI
    use inputs_mod, only : blockhandout
    use parallel_mod, only: mpiBlockHandout, mpiGetBlock
    use mpi_global_mod, only : myRankGlobal, nThreadsGlobal
    use mpi
#endif
    type(octal), pointer   :: thisOctal, neighbourOctal, startOctal
    type(GRIDTYPE) :: grid
    real(double) ::  deMax
    real :: tol
    logical :: converged
    integer :: subcell, neighbourSubcell
    real(double) :: eDens(-1:1,-1:1,-1:1)
    real(double) :: dCoeff(-1:1,-1:1,-1:1)
    real(double) :: dCoeffhalf(-1:1,-1:1,-1:1)
    real(double) :: eNplus1
    real(double) :: r, gradE, bigR, lambda
    real(double) :: DeltaT, DeltaX, deltaE
    logical :: firstTime
    type(VECTOR) :: octVec, rVec
    integer                     :: nOctal        ! number of octals in grid
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: iOctal
    integer :: iOctal_beg, iOctal_end
    real(double) :: phi, DeltaPhi, DeltaZ, deltaR
    real(double), parameter :: underCorrect = 0.8d0
    real(double), save :: overCorrect = 1.d0
#ifdef MPI
    integer :: my_rank, np !, isubcell
    real(double) :: globalDeMax
    integer, dimension(:), allocatable :: octalsBelongRank
    logical :: rankComplete
    integer :: ierr
    integer :: tag = 0
    real(double), allocatable :: eArray(:), tArray(:)
    integer :: nEdens, nVoxels
#endif

    firstTime = .true.
    deMax = -1.d30

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif

#ifdef MPI
    ! FOR MPI IMPLEMENTATION=======================================================
    np      = nThreadsGlobal
    my_rank = myRankGlobal

    ! we will use an array to store the rank of the process
    !   which will calculate each octal's variables
    allocate(octalsBelongRank(size(octalArray)))
    
    if (my_rank == 0) then
  !     print *, ' '
  !     print *, 'Gauss-Seidel sweep  computed by ', np-1, ' processors.'
  !     print *, ' '
       call mpiBlockHandout(np,octalsBelongRank,blockDivFactor=1,tag=tag,&
                            setDebug=.false.)
    
    endif
    ! ============================================================================
#endif
    
    ! default loop indices
    ioctal_beg = 1
    ioctal_end = nOctal

#ifdef MPI
 if (my_rank /= 0) then
  blockLoop: do     
 call mpiGetBlock(my_rank,iOctal_beg,iOctal_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit blockLoop 
#endif


!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(iOctal, eDens, dCoeff,subcell, thisOctal, octVec, startOctal, r) &
!$OMP PRIVATE(neighbourOctal, neighbourSubcell, enplus1, firsttime, deltax, deltat) &
!$OMP PRIVATE(deltaE, grade, dcoeffhalf, lambda, bigr) &
!$OMP PRIVATE(phi, DeltaPhi, rVec, DeltaR, DeltaZ) &
!$OMP SHARED(overCorrect, octalArray) & 
!$OMP SHARED(grid, tol, demax, ioctal_beg, ioctal_end) 
!$OMP DO SCHEDULE(static)
    do iOctal =  iOctal_beg, iOctal_end

       thisOctal => octalArray(iOctal)%content

       do subcell = 1, thisOctal%maxChildren

          if (.not.thisOctal%hasChild(subcell)) then

             eDens = 0.d0
             dCoeff = 0.d0

             if (thisOctal%diffusionApprox(subcell).and.(.not.thisOctal%fixedtemperature(subcell))) then

                if (thisOctal%twoD) then
                   r = thisOctal%subcellSize/2.d0 + grid%halfSmallestSubcell * 0.1d0

                   eDens(0, 0, 0) = thisOctal%eDens(subcell)
                   dCoeff(0, 0, 0) = thisOctal%diffusionCoeff(subcell)

                   ! positive x

                   octVec = subcellCentre(thisOctal, subcell) + VECTOR(r, 0.d0, 0.d0)
                   if (inOctal(grid%octreeRoot, octVec)) then
                      startOctal => thisOctal
                      call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                           foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                      eDens(1, 0, 0) = neighbourOctal%eDens(neighboursubcell)
                      dCoeff(1, 0, 0) = neighbourOctal%diffusionCoeff(neighboursubcell)
                   else
                      eDens(1, 0, 0) = eDens(0, 0, 0)
                      dCoeff(1, 0, 0) = dCoeff(0, 0, 0)
                   endif
                   ! negative x

                   octVec = subcellCentre(thisOctal, subcell) - VECTOR(r, 0.d0, 0.d0)
                   if (inOctal(grid%octreeRoot, octVec)) then
                      startOctal => thisOctal
                      call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                           foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                      eDens(-1, 0, 0) = neighbourOctal%eDens(neighbourSubcell)
                      dCoeff(-1, 0, 0) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                   else
                      eDens(-1, 0, 0) = eDens(0, 0, 0)
                      dCoeff(-1, 0, 0) = dCoeff(0, 0, 0)
                   endif

                   ! positive z

                   octVec = subcellCentre(thisOctal, subcell) + VECTOR(0.d0, 0.d0, r)
                   if (inOctal(grid%octreeRoot, octVec)) then
                      startOctal => thisOctal
                      call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                           foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                      eDens(0, 0, 1) = neighbourOctal%eDens(neighbourSubcell)
                      dCoeff(0, 0, 1) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                   else
                      eDens(0, 0, 1) = eDens(0, 0, 0)
                      dCoeff(0, 0, 1) = dCoeff(0, 0, 0)
                   endif

                   ! negative z

                   octVec = subcellCentre(thisOctal, subcell) - VECTOR(0.d0, 0.d0, r)
                   if (inOctal(grid%octreeRoot, octVec)) then
                      startOctal => thisOctal
                      call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                           foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                      eDens(0,  0, -1) = neighbourOctal%eDens(neighbourSubcell)
                      dCoeff(0, 0, -1) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                   else
                      eDens(0, 0, -1) = eDens(0, 0, 0)
                      dCoeff(0, 0, -1) = dCoeff(0, 0, 0)
                   endif

                else


                   if (.not.thisOctal%cylindrical) then
                      r = thisOctal%subcellSize/2.d0 + grid%halfSmallestSubcell * 0.1d0

                      eDens(0, 0, 0) = thisOctal%eDens(subcell)
                      dCoeff(0, 0, 0) = thisOctal%diffusionCoeff(subcell)

                      ! positive x
                      
                      octVec = subcellCentre(thisOctal, subcell) + VECTOR(r, 0.d0, 0.d0)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(1, 0, 0) = neighbourOctal%eDens(neighboursubcell)
                         dCoeff(1, 0, 0) = neighbourOctal%diffusionCoeff(neighboursubcell)
                      else
                         eDens(1, 0, 0) = eDens(0, 0, 0)
                         dCoeff(1, 0, 0) = dCoeff(0, 0, 0)
                      endif

                   ! negative x

                      octVec = subcellCentre(thisOctal, subcell) - VECTOR(r, 0.d0, 0.d0)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(-1, 0, 0) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(-1, 0, 0) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(-1, 0, 0) = eDens(0, 0, 0)
                         dCoeff(-1, 0, 0) = dCoeff(0, 0, 0)
                      endif
                      ! positive y

                      octVec = subcellCentre(thisOctal, subcell) + VECTOR(0.d0, r, 0.d0)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(0, 1, 0) = neighbourOctal%eDens(neighboursubcell)
                         dCoeff(0, 1, 0) = neighbourOctal%diffusionCoeff(neighboursubcell)
                      else
                         eDens(0, 1, 0) = eDens(0, 0, 0)
                         dCoeff(0, 1, 0) = dCoeff(0, 0, 0)
                      endif
                      ! negative y
                      
                      octVec = subcellCentre(thisOctal, subcell) - VECTOR(0.d0, r, 0.d0)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(0, -1, 0) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(0, -1, 0) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(0, -1, 0) = eDens(0, 0, 0)
                         dCoeff(0, -1, 0) = dCoeff(0, 0, 0)
                      endif

                      ! positive z

                      octVec = subcellCentre(thisOctal, subcell) + VECTOR(0.d0, 0.d0, r)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(0, 0, 1) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(0, 0, 1) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(0, 0, 1) = eDens(0, 0, 0)
                         dCoeff(0, 0, 1) = dCoeff(0, 0, 0)
                      endif

                      ! negative z

                      octVec = subcellCentre(thisOctal, subcell) - VECTOR(0.d0, 0.d0, r)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(0,  0, -1) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(0, 0, -1) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(0, 0, -1) = eDens(0, 0, 0)
                         dCoeff(0, 0, -1) = dCoeff(0, 0, 0)
                      endif

                   else


                      eDens(0, 0, 0) = thisOctal%eDens(subcell)
                      dCoeff(0, 0, 0) = thisOctal%diffusionCoeff(subcell)

                      r = thisOctal%subcellSize/2.d0 + grid%halfSmallestSubcell * 0.1d0


                      ! positive r
                      
                      octVec = subcellCentre(thisOctal, subcell) + &
                               VECTOR(r, grid%halfSmallestSubcell, grid%halfSmallestSubcell)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(1, 0, 0) = neighbourOctal%eDens(neighboursubcell)
                         dCoeff(1, 0, 0) = neighbourOctal%diffusionCoeff(neighboursubcell)
                      else
                         eDens(1, 0, 0) = eDens(0, 0, 0)
                         dCoeff(1, 0, 0) = dCoeff(0, 0, 0)
                      endif

                      ! negative r

                      octVec = subcellCentre(thisOctal, subcell) - &
                               VECTOR(r, grid%halfSmallestSubcell, grid%halfSmallestSubcell)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(-1, 0, 0) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(-1, 0, 0) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(-1, 0, 0) = eDens(0, 0, 0)
                         dCoeff(-1, 0, 0) = dCoeff(0, 0, 0)
                      endif


                      ! do the phi changes
                      
                      octVec = subcellCentre(thisOctal, subcell)
                      phi = atan2(octVec%y,octVec%x)
                      if (phi < 0.d0) phi = phi + twopi

                      DeltaPhi = returndPhi(thisOctal) * 1.01d0
!                      if (thisOctal%splitAzimuthally) then
!                         DeltaPhi = thisOctal%dPhi/4.d0 * 1.01d0
!                      else
!                         DeltaPhi = thisOctal%dPhi/2.d0 * 1.01d0
!                      endif

                      ! negative phi

                      rVec = octVec + VECTOR(grid%halfSmallestSubcell, 0.d0, grid%halfSmallestSubcell)
                      rVec = rotateZ(rVec,+DeltaPhi) ! remember that rotateZ rotates in the "wrong" sense
                      if (inOctal(grid%octreeRoot, rVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, rVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(0, -1, 0) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(0, -1, 0) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(0, -1, 0) = eDens(0, 0, 0)
                         dCoeff(0, -1, 0) = dCoeff(0, 0, 0)
                      endif

                      ! positive phi

                      rVec = octVec + VECTOR(grid%halfSmallestSubcell, 0.d0, grid%halfSmallestSubcell)
                      rVec =  rotateZ(rVec,-DeltaPhi)  ! remember that rotateZ rotates in the "wrong" sense
                      if (inOctal(grid%octreeRoot, rVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, rVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(0, 1, 0) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(0, 1, 0) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(0, 1, 0) = eDens(0, 0, 0)
                         dCoeff(0, 1, 0) = dCoeff(0, 0, 0)
                      endif


                      ! positive z

                      octVec = subcellCentre(thisOctal, subcell) + &
                               VECTOR(grid%halfSmallestSubcell, grid%halfSmallestSubcell, r)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(0, 0, 1) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(0, 0, 1) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(0, 0, 1) = eDens(0, 0, 0)
                         dCoeff(0, 0, 1) = dCoeff(0, 0, 0)
                      endif

                      ! negative z

                      octVec = subcellCentre(thisOctal, subcell) - &
                               VECTOR(grid%halfSmallestSubcell, grid%halfSmallestSubcell, r)
                      if (inOctal(grid%octreeRoot, octVec)) then
                         startOctal => thisOctal
                         call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                              foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)
                         eDens(0,  0, -1) = neighbourOctal%eDens(neighbourSubcell)
                         dCoeff(0, 0, -1) = neighbourOctal%diffusionCoeff(neighbourSubcell)
                      else
                         eDens(0, 0, -1) = eDens(0, 0, 0)
                         dCoeff(0, 0, -1) = dCoeff(0, 0, 0)
                      endif

                      DeltaR = 1.d10 * (thisOctal%subcellSize/2.d0 + grid%halfSmallestSubcell * 0.1d0)
                      DeltaZ = 1.d10 * (thisOctal%subcellSize/2.d0 + grid%halfSmallestSubcell * 0.1d0)

                      rVec = subcellCentre(thisOctal, subcell)
                      r = 1.d10 * sqrt(rVec%x**2 + rVec%y**2)
                   endif

                endif

                if (.not.thisOctal%cylindrical) then
                   gradE = sqrt((eDens(1,0,0)-eDens(-1,0,0))**2 + &
                        (eDens(0,1,0)-eDens(0,-1,0))**2 + &
                        (eDens(0,0,1)-eDens(0,0,-1))**2)
                   gradE = abs(gradE) / (2.d0*thisOctal%subcellsize*1.e10)
                else
                   gradE = ((eDens(1,0,0)-eDens(-1,0,0))/(2.d0*DeltaR))**2 + &
                        ((1.d0/r)*(eDens(0,1,0)-eDens(0,-1,0)/(2.d0*DeltaPhi)))**2 + &
                        ((eDens(0,0,1)-eDens(0,0,-1))/(2.d0*DeltaZ))**2
                   gradE = sqrt(gradE)
                endif

                bigR = gradE / (eDens(0,0,0) * (thisOctal%kappaRoss(subcell) * thisOctal%rho(subcell)))
                lambda = (2.d0 + bigR) / (6.d0 + 3.d0*bigR + bigR**2)

                lambda = 0.33333d0
                dCoeffHalf = 0.d0

                dCoeffHalf(1, 0, 0) = 0.5d0 * (dCoeff(1, 0, 0) + dCoeff(0, 0, 0))
                dCoeffHalf(-1, 0, 0) = 0.5d0 * (dCoeff(0, 0, 0) + dCoeff(-1, 0, 0))

                dCoeffHalf(0, 1, 0) = 0.5d0 * (dCoeff(0, 1, 0) + dCoeff(0, 0, 0))
                dCoeffHalf(0, -1, 0) = 0.5d0 * (dCoeff(0, 0, 0) + dCoeff(0, -1, 0))

                dCoeffHalf(0, 0, 1) = 0.5d0 * (dCoeff(0, 0, 1) + dCoeff(0, 0, 0))
                dCoeffHalf(0, 0, -1) = 0.5d0 * (dCoeff(0, 0, 0) + dCoeff(0, 0, -1))

                dCoeffHalf(-1:1,-1:1,-1:1) = dCoeffHalf(-1:1,-1:1,-1:1) * lambda

                if (.not.thisOctal%cylindrical) then
                   DeltaX = r * 1.d10
                   DeltaX = thisOctal%subcellsize * 1.d10
                   DeltaT = DeltaX**2 / (2.d0 * maxval(dcoeffHalf(-1:1,-1:1,-1:1)))
                   DeltaT = DeltaT * 0.1
                else
!                   deltaX = min(DeltaR, deltaZ, r * DeltaPhi)

                   deltaX = min(DeltaR, deltaZ) ! TJH Changed 21/1/11
                   DeltaT = DeltaX**2 / (2.d0 * maxval(dcoeffHalf(-1:1,-1:1,-1:1)))
                   DeltaT = DeltaT * 0.1
                endif

                if (thisOctal%twoD) then
                   enplus1 = eDens(0,0,0) + (DeltaT/DeltaX**2) * OverCorrect * (dCoeffHalf(1,0,0)*(eDens(1,0,0)-eDens(0,0,0)) &
                        - dCoeffHalf(-1,0,0)*(eDens(0,0,0)-eDens(-1,0,0)) &
                        + dCoeffHalf(0,0,1)*(eDens(0,0,1)-eDens(0,0,0)) &
                        - dCoeffHalf(0,0,-1)*(eDens(0,0,0)-eDens(0,0,-1)))
                else
                   if (.not.thisOctal%cylindrical) then
                      enplus1 = eDens(0,0,0) + (DeltaT/DeltaX**2) * &
                           (  dCoeffHalf(1,0,0)*(eDens(1,0,0)-eDens(0,0,0)) &
                           - dCoeffHalf(-1,0,0)*(eDens(0,0,0)-eDens(-1,0,0)) &
                           + dCoeffHalf(0,1,0)*(eDens(0,1,0)-eDens(0,0,0)) &
                           - dCoeffHalf(0,-1,0)*(eDens(0,0,0)-eDens(0,-1,0)) &
                           + dCoeffHalf(0,0,1)*(eDens(0,0,1)-eDens(0,0,0)) &
                           - dCoeffHalf(0,0,-1)*(eDens(0,0,0)-eDens(0,0,-1)))
                   else                      
!                      enplus1 = eDens(0,0,0) + deltaT * ( &
!                           (dCoeffHalf(1,0,0)*(eDens(1,0,0)-eDens(0,0,0))- &
!                           dCoeffHalf(-1,0,0)*(eDens(0,0,0)-eDens(-1,0,0)))/DeltaR**2 + &
!                           (1.d0/r)*(dCoeffHalf(0,0,0) * (eDens(1,0,0)-eDens(-1,0,0)) / (2.d0 * DeltaR)) + &
!                           (1.d0/r**2)*(dCoeffHalf(0,1,0)*(eDens(0,1,0)-eDens(0,0,0)) - &
!                           dCoeffHalf(0,-1,0)*(eDens(0,0,0)-eDens(0,-1,0)))/DeltaPhi**2 + &
!                           (dCoeffHalf(0,0,1)*(eDens(0,0,1)-eDens(0,0,0)) - &
!                           dCoeffHalf(0,0,-1)*(eDens(0,0,0)-eDens(0,0,-1)))/deltaZ**2 )


 ! TJH Changed 21/1/11
                   enplus1 = eDens(0,0,0) + (DeltaT/DeltaR**2) * OverCorrect * (dCoeffHalf(1,0,0)*(eDens(1,0,0)-eDens(0,0,0)) &
                        - dCoeffHalf(-1,0,0)*(eDens(0,0,0)-eDens(-1,0,0)) &
                        + dCoeffHalf(0,0,1)*(eDens(0,0,1)-eDens(0,0,0)) &
                        - dCoeffHalf(0,0,-1)*(eDens(0,0,0)-eDens(0,0,-1)))

                   endif
                endif
                overcorrect = min(4.d0,overcorrect + 0.1d0)

                if (enPlus1 < 0.d0) then
                   if (firstTime) write(*,*) "Warning: negative energy density."
                   overcorrect = 1.d0
                   firstTime = .false.
                   deltaE = (enPlus1-thisOctal%oldeDens(subcell))
                   enPlus1 = enPlus1 + undercorrect * deltaE ! undercorrect here
                   if (enPlus1 < 0.d0) then
                      enPlus1 = arad*(10.d0**4)
                   endif
                endif
                !$OMP CRITICAL
                thisOctal%eDens(subcell) = enPlus1
                thisOctal%temperature(subcell) = real(sqrt(sqrt(enPlus1 * OneOveraRad)))
                !$OMP END CRITICAL
                deltaE = abs(enPlus1-thisOctal%oldeDens(subcell)) &
                     / thisOctal%oldEdens(subcell)
                deMax = max(deMax, deltaE)
             endif
          endif
       enddo
    enddo
!$OMP END DO
!$OMP BARRIER
!$OMP END PARALLEL 

#ifdef MPI
 if (.not.blockHandout) exit blockloop
 end do blockLoop       
 end if ! (my_rank /= 0)


    call MPI_ALLREDUCE(deMax,globalDeMax,1,MPI_DOUBLE_PRECISION,&
         MPI_MAX,MPI_COMM_WORLD,ierr)
    deMax = globalDeMax

 !    print *,'Process ',my_rank,' waiting to update values in Gauss-Seidel sweep...' 
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     ! have to send out the 'octalsBelongRank' array
     call MPI_BCAST(octalsBelongRank,SIZE(octalsBelongRank),  &
                    MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     allocate(eArray(1:nVoxels))
     allocate(tArray(1:nVoxels))
     eArray = 0.d0
     call packEdens(octalArray, nEdens, eArray,octalsBelongRank)
     call MPI_ALLREDUCE(eArray,tArray,nEdens,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
     eArray = tArray
     call unpackEdens(octalArray, nEdens, eArray)
     deallocate(eArray, tArray)
 !      !
 !      ! Update the edens values of grid computed by all processors.
 !      !
 !      do iOctal = 1, SIZE(octalArray)
 !    !     print *,'Process ',my_rank,' starting octal ',iOctal 
 !
 !         thisOctal => octalArray(iOctal)%content
 !         
 !         do iSubcell = 1, thisOctal%maxChildren
 !            if (octalArray(iOctal)%inUse(iSubcell).and. &
 !                    octalArray(iOctal)%content%diffusionApprox(iSubcell)) then
 !               
 !               call MPI_BCAST(thisOctal%eDens(iSubcell), 1, MPI_DOUBLE_PRECISION,&
 !                    octalsBelongRank(iOctal), MPI_COMM_WORLD, ierr)
 !            end if
 !         
 !         end do
 !      end do
          
          
 !    print *,'Process ',my_rank,' finished updating values in Gauss-Seidel sweep...' 
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
#endif
!    if (writeoutput) write(*,*) "Demax ",demax
    if (deMax > tol) converged = .false.
    deallocate(octalArray)
#ifdef MPI
    deallocate(octalsBelongRank)
#endif
end subroutine gaussSeidelSweep

  subroutine solveArbitraryDiffusionZones(grid)
    use inputs_mod, only : eDensTol !, tauforce
    use messages_mod, only : myRankIsZero

    type(GRIDTYPE) :: grid
    character(len=80) :: message
    logical :: gridConverged
    real(double) :: deMax
    integer :: niter, i
    integer, parameter :: maxIter = 10000
    
!    logical, save :: firstTime = .true.

    call seteDens(grid, grid%octreeRoot)
    call setDiffusionCoeff(grid, grid%octreeRoot)



!    if (firstTime) then
!       call defineDiffusionOnRosseland(grid, grid%octreeRoot, 0.5)
!       firstTime = .false.
!    else
!       call defineDiffusionOnRosseland(grid, grid%octreeRoot, tauforce)


    call setDiffOnTau(grid)
!    if (writeoutput) call writeVtkFile(grid, "bias.vtk", valueTypeString=(/"chiline    ","temperature"/))

!    endif
!    call defineDiffusiononRho(grid%octreeRoot, 1.d-10)
!       call defineDiffusionOnUndersampled(grid%octreeRoot)
!    call resetDiffusionTemp(grid%octreeRoot, 100.)

    i  = 0
    call countDiffusionCells(grid%octreeRoot, i)
    write(message,*) "Number of cells in diffusion zone: ",i
    call writeInfo(message, TRIVIAL)

    gridconverged = .false.
    nIter = 0
     do while (.not.gridconverged)
        nIter = nIter + 1
        gridconverged = .true.
        deMax = -1.e30
        call gaussSeidelSweep(grid, edenstol, demax, gridConverged)
!        call copyEdens(grid%octreeRoot)
        call setDiffusionCoeff(grid, grid%octreeRoot)
!        write(message,*) nIter," Maximum relative change in eDens:",deMax
!	call writeInfo(message, TRIVIAL)
        if (nIter < 3) gridConverged = .false.
        if (nIter > maxIter) then
           if (myRankIsZero) write(*,*) "No solution found after ",maxIter," iterations"
           gridConverged = .true.
        endif
     enddo
     write(message,*) "Gauss-seidel took ",niter, " iterations"
     call writeInfo(message,TRIVIAL)
   end subroutine solveArbitraryDiffusionZones


  recursive subroutine defineDiffusionOnUndersampled(thisOctal, nDiff)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell
    integer, optional :: nDiff
    integer :: i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call defineDiffusionOnUndersampled(child, ndiff)
                exit
             end if
          end do
       else
          if (thisOctal%undersampled(subcell)) then
             thisOctal%diffusionApprox(subcell) = .true.
             if (PRESENT(nDiff)) then
                nDiff = nDiff + 1
             endif
          endif
       end if
    end do

  end subroutine defineDiffusionOnUndersampled

  recursive subroutine defineDiffusionOnRho(thisOctal, rhoLim)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: rhoLim
    integer :: subcell
    integer :: i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call defineDiffusionOnRho(child, rhoLim)
                exit
             end if
          end do
       else
          if (thisOctal%rho(subcell) > rhoLim) then
             thisOctal%diffusionApprox(subcell) = .true.
          endif
       end if
    end do

  end subroutine defineDiffusionOnRho

  recursive subroutine defineDiffusionOnRosseland(grid, thisOctal, taudiff, nDiff)
    use inputs_mod, only :  resetDiffusion
    real :: tauDiff
    type(GRIDTYPE) :: grid
    integer, optional :: nDiff
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell
    integer :: i
    real(double) :: kRos, tau
    kros = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call defineDiffusionOnRosseland(grid, child, taudiff, nDiff)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          call returnKappa(grid, thisOctal, subcell, rosselandKappa=kros)
          tau = kros * thisOctal%rho(subcell) * thisOctal%subcellSize *1.e10
          if (.not.associated(thisOctal%diffusionApprox)) then
             allocate(thisOCtal%diffusionApprox(1:thisOctal%maxChildren))
             thisOctal%diffusionApprox = .false.
          endif
          if (tau > taudiff) then
             thisOctal%diffusionApprox(subcell) = .true.
             if (PRESENT(ndiff))  nDiff = nDiff + 1
          else
             if (resetDiffusion) thisOctal%diffusionApprox(subcell) = .false.
          endif
       end if
    end do

  end subroutine defineDiffusionOnRosseland

  recursive subroutine defineDiffusionOnKappap(grid, thisOctal, taudiff, nDiff)
    use inputs_mod, only :  resetDiffusion
    real :: tauDiff
    type(GRIDTYPE) :: grid
    integer, optional :: nDiff
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell
    integer :: i
    real(double) :: tau
    real :: kappap
    kappap = 0.d0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call defineDiffusionOnKappap(grid, child, taudiff, nDiff)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          call returnKappa(grid, thisOctal, subcell, kappap=kappap)
          tau = dble(kappap) * thisOctal%subcellSize *1.e10
          if (.not.associated(thisOctal%diffusionApprox)) then
             allocate(thisOCtal%diffusionApprox(1:thisOctal%maxChildren))
             thisOctal%diffusionApprox = .false.
          endif
          if (tau > taudiff) then
             thisOctal%diffusionApprox(subcell) = .true.
             if (PRESENT(ndiff))  nDiff = nDiff + 1
          else
             if (resetDiffusion) thisOctal%diffusionApprox(subcell) = .false.
          endif
       end if
    end do

  end subroutine defineDiffusionOnKappap

  recursive subroutine resetDiffusionTemp(thisOctal, temp)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell
    integer :: i
    real :: temp
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call resetDiffusionTemp(child, temp)
                exit
             end if
          end do
       else
          if (thisOctal%diffusionApprox(subcell)) then
             thisOctal%temperature(subcell) = temp
             thisOctal%eDens(subcell) = arad * temp**4
          endif
       end if
    end do

  end subroutine resetDiffusionTemp

  recursive subroutine copyEdens(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell
    integer :: i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call copyEdens(child)
                exit
             end if
          end do
       else
          if (thisOctal%diffusionApprox(subcell)) then
             thisOctal%eDens(subcell) = thisOctal%chiLine(subcell)
          endif
       end if
    end do

  end subroutine copyEdens


  subroutine randomWalk(grid, startOctal, startSubcell,  endOctal, endSubcell, temp, ok)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: startOctal
    type(OCTAL), pointer :: endOctal
    type(OCTAL), pointer :: walkOctal, sOctal
    type(VECTOR) :: rVec
    type(VECTOR) :: xAxis, yAxis, zAxis, rHat
    logical, intent(out) :: ok

    real, intent(out) :: temp
    integer :: startSubcell, walkSubcell, nWalk
    integer, intent(out) :: endSubcell
    real :: r
    integer, parameter :: maxWalk = 10000

    ok = .true.

    xAxis = VECTOR(1.d0, 0.d0, 0.d0)
    yAxis = VECTOR(0.d0, 1.d0, 0.d0)
    zAxis = VECTOR(0.d0, 0.d0, 1.d0)
    walkOctal => startOctal
    walkSubcell = startSubcell
    rVec = subcellCentre(walkOctal, walkSubcell)
    nWalk  = 0

    call randomNumberGenerator(getReal=r)
    do while(walkOctal%diffusionApprox(walkSubcell))
       nWalk = nWalk + 1

       if (startOctal%twoD) then
          call randomNumberGenerator(getReal=r)
          rVec = subcellCentre(walkOctal, walkSubcell)
          
          if (r < 0.25) then
             rVec = rVec + (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*xAxis
          else if ((r >= 0.25).and.(r < 0.5)) then
             rVec = rVec - (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*xAxis
          else if ((r >= 0.5).and.(r < 0.75)) then
             rVec = rVec + (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*zAxis
          else if (r >= 0.75) then
             rVec = rVec - (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*zAxis
          endif
          if (nWalk > 1000) then
!             write(*,*) nWalk,r,rVec,inOctal(grid%octreeRoot,rVec)
             call randomNumberGenerator(getReal=r)             
          endif
          sOctal => walkOctal
          temp = walkOctal%temperature(walkSubcell)
          if (.not.inOctal(grid%octreeRoot, rVec)) then
             rVec = subcellCentre( walkOctal, walkSubcell)
             call randomNumberGenerator(getReal=r)
          endif
          call amrgridvalues(grid%octreeRoot, rVec, startOctal=sOctal, foundOctal=walkOctal, foundSubcell=walkSubcell)
       else

          if (.not.walkOctal%cylindrical) then
             call randomNumberGenerator(getReal=r)
             rVec = subcellCentre(walkOctal, walkSubcell)
             
             if (r < 0.166666) then
                rVec = rVec + (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*xAxis
             else if ((r >= 0.166666).and.(r < 0.333333)) then
                rVec = rVec - (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*xAxis
             else if ((r >= 0.333333).and.(r < 0.5)) then
                rVec = rVec + (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*yAxis
             else if ((r >= 0.5).and.(r < 0.666666)) then
                rVec = rVec - (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*yAxis
             else if ((r >= 0.6666666).and.(r < 0.833333)) then
                rVec = rVec + (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*zAxis
             else if (r >= 0.8333333) then
                rVec = rVec - (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*zAxis
             endif
             sOctal => walkOctal
             temp = walkOctal%temperature(walkSubcell)
             if (.not.inOctal(grid%octreeRoot, rVec)) then
                rVec = subcellCentre( walkOctal, walkSubcell)
                call randomNumberGenerator(getReal=r)
             endif
             call amrgridvalues(grid%octreeRoot, rVec, startOctal=sOctal, foundOctal=walkOctal, foundSubcell=walkSubcell)
          else
             call randomNumberGenerator(getReal=r)
             rVec = subcellCentre(walkOctal, walkSubcell)
             rHat = VECTOR(rVec%x,rVec%y, 0.d0)
             call normalize(rHat)
             if (r < 0.25) then
                rVec = rVec + (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*rHat
             else if ((r >= 0.25).and.(r < 0.5)) then
                rVec = rVec - (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*rHat
             else if ((r >= 0.5).and.(r < 0.75)) then
                rVec = rVec + (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*zAxis
             else if (r >= 0.75) then
                rVec = rVec - (walkOctal%subcellSize/2.d0+grid%halfSmallestsubcell)*zAxis
             endif
             sOctal => walkOctal
             temp = walkOctal%temperature(walkSubcell)
             if (.not.inOctal(grid%octreeRoot, rVec)) then
                rVec = subcellCentre( walkOctal, walkSubcell)
                call randomNumberGenerator(getReal=r)
             endif
             call amrgridvalues(grid%octreeRoot, rVec, startOctal=sOctal, foundOctal=walkOctal, foundSubcell=walkSubcell)
          endif
       endif
       if (nWalk > maxWalk) then
          call writeWarning("Photon exceeded random walk limit")
          ok = .false.
          exit
       endif
    enddo
    endOctal => walkOctal
    endSubcell = walkSubcell

  end subroutine randomWalk

  recursive subroutine unsetOnDirect(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell
    integer :: i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unsetOnDirect(child)
                exit
             end if
          end do
       else
          if (thisOctal%nDirectPhotons(subcell) > 0) thisOctal%diffusionApprox(subcell) = .false.
       end if
    end do

  end subroutine unsetOnDirect

  recursive subroutine unsetDiffusion(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell
    integer :: i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unsetDiffusion(child)
                exit
             end if
          end do
       else
          thisOctal%diffusionApprox(subcell) = .false.
       end if
    end do

  end subroutine unsetDiffusion

#ifdef MPI
      subroutine packEdens(octalArray, nEdens, eArray, octalsBelongRank)
        use mpi_global_mod, only : myRankGlobal

        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
        integer :: nEdens
        real(double) :: eArray(:)
        integer :: iOctal, iSubcell, my_rank
        type(OCTAL), pointer :: thisOctal

       !
       ! Update the edens values of grid computed by all processors.
       !
       my_rank = myRankGlobal
       nEdens = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
             if (octalArray(iOctal)%inUse(iSubcell).and. &
                     octalArray(iOctal)%content%diffusionApprox(iSubcell)) then
                
                 nEdens = nEdens + 1
                 if (octalsBelongRank(iOctal) == my_rank) then
                   eArray(nEdens) = octalArray(iOctal)%content%eDens(iSubcell)
                 else 
                   eArray(nEdens) = 0.d0
                 endif
             end if
          
          end do
       end do
     end subroutine packEdens

      subroutine unpackEdens(octalArray, nEdens, eArray)
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: nEdens
        real(double) :: eArray(:)
        integer :: iOctal, iSubcell
        type(OCTAL), pointer :: thisOctal

       !
       ! Update the edens values of grid computed by all processors.
       !
       nEdens = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
             if (octalArray(iOctal)%inUse(iSubcell).and. &
                     octalArray(iOctal)%content%diffusionApprox(iSubcell)) then
                
                 nEdens = nEdens + 1
                 octalArray(iOctal)%content%eDens(iSubcell) = eArray(nEdens)
             end if
          
          end do
       end do
     end subroutine unpackEdens
#endif

subroutine setDiffOnTau(grid)
    use inputs_mod, only : tauForce, cylindrical, rGapInner, rGapOuter
#ifdef MPI
    use inputs_mod, only : blockHandout
    use parallel_mod, only: mpiBlockHandout, mpiGetBlock
    use mpi
#endif
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: tau, thisTau
    type(VECTOR) :: rVec, direction
    integer :: i
    integer :: subcell
    integer :: iLambda
    integer                     :: nOctal        ! number of octals in grid
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: iOctal
    integer :: iOctal_beg, iOctal_end
    real(double) :: kappaAbs, r 
    type(VECTOR) :: arrayVec(6)
     integer :: nDir
#ifdef MPI
! Only declared in MPI case
     integer, dimension(:), allocatable :: octalsBelongRank
     logical :: rankComplete
     integer :: tag = 0
     logical, allocatable:: diffArray(:), tArray(:)
     real(double), allocatable :: tauArray(:), trArray(:)
     integer :: nVoxels, ierr
     integer :: nDiff
     integer :: np
     integer :: my_rank
     np = 1
     my_rank = 1
#endif


     iLambda = 0
     kappaAbs = 0.d0
     thisTau = 0.d0
    
    ndir = 4
     arrayVec(1) = VECTOR(1.d0, 1.d-10, 1.d-10)
     arrayVec(2) = VECTOR(-1.d0, 1.d-10, 1.d-10)
     arrayVec(3) = VECTOR(1.d-10, 1.d-10, 1.d0)
     arrayVec(4) = VECTOR(1.d-10, 1.d-10,-1.d0)

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif

#ifdef MPI
    ! FOR MPI IMPLEMENTATION=======================================================
    !  Get my process rank # 
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
    ! Find the total # of precessor being used in this run
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
    
    ! we will use an array to store the rank of the process
    !   which will calculate each octal's variables
    allocate(octalsBelongRank(size(octalArray)))
    
    if (my_rank == 0) then
!       print *, ' '
!       print *, 'Diffusion on Tau computed by ', np-1, ' processors.'
!       print *, ' '
       call mpiBlockHandout(np,octalsBelongRank,blockDivFactor=1,tag=tag,&
                            setDebug=.false.)
    
    endif
    ! ============================================================================
#endif
    
    ! default loop indices
    ioctal_beg = 1
    ioctal_end = nOctal

#ifdef MPI
 if (my_rank /= 0) then
  blockLoop: do     
 call mpiGetBlock(my_rank,iOctal_beg,iOctal_end,rankComplete,tag,setDebug=.false.)
   if (rankComplete) exit blockLoop 
#endif
    do iOctal =  iOctal_beg, iOctal_end
       thisOctal => octalArray(iOctal)%content

       do subcell = 1, thisOctal%maxChildren

          if (.not.thisOctal%hasChild(subcell)) then

             rVec = subcellCentre(thisOctal, subcell)

             r = modulus(rVec)
             if (grid%geometry == "shakara") then 
                if ((r > rGapInner).and.(r < rGapOuter)) then
                   thisOctal%diffusionApprox(subcell) = .false.
                   cycle
                endif
             end if

             if (thisOctal%threed) then
                rVec = rVec + 0.01d0*grid%halfSmallestSubcell * VECTOR(0.1d0,0.1d0,0.1d0)
             endif
	
             call returnKappa(grid, thisOctal, subcell,  rosselandKappa = kappaAbs)
             tau = thisOctal%subcellSize*kappaAbs*thisOctal%rho(subcell)*1.d10

             if (tau < 0.01d0) then
                thisOctal%diffusionApprox(subcell) = .false.
             else
                tau = 1.d30
                
                if (cylindrical) then
                   ndir = 4
                   arrayVec(1) = VECTOR(0.d0, 0.d0, 1.d0)
                   arrayVec(2) = VECTOR(0.d0, 0.d0,-1.d0)
                   arrayVec(3) = VECTOR(rVec%x, rVec%y,0.d0)
                   call normalize(arrayVec(3))
                   arrayVec(4) = (-1.d0)*arrayVec(3)
                   arrayVec(3) = rotateZ(arrayVec(3), 0.1d0*degtorad)
                   arrayVec(4) = rotateZ(arrayVec(4), 0.1d0*degtorad)
!                   arrayVec(5) = arrayVec(3).cross.arrayVec(1)
!                   call normalize(arrayVec(5))
!                   arrayVec(6) = (-1.d0)*arrayVec(5)
                endif


                

                do i = 1, nDir
                   direction = arrayVec(i)
                   call tauAlongPath(ilambda, grid, rVec, direction, thistau, 100.d0, ross=.true., stopAtGap=.true.)
                   tau = min(tau, thisTau)
                enddo
                thisOctal%chiLine(subcell) = tau
                if (tau > tauForce) then
                   thisOctal%diffusionApprox(subcell) = .true.
                else
                   thisOctal%diffusionApprox(subcell) = .false.
                endif


             endif

          endif

       enddo
    enddo

#ifdef MPI
 if (.not.blockHandout) exit blockloop
 end do blockLoop       
 end if ! (my_rank /= 0)


     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     ! have to send out the 'octalsBelongRank' array
     call MPI_BCAST(octalsBelongRank,SIZE(octalsBelongRank),  &
                    MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

     call countVoxels(grid%octreeRoot,nOctal,nVoxels)
     allocate(diffArray(1:nVoxels))
     allocate(tauArray(1:nVoxels))
     allocate(tArray(1:nVoxels))
     allocate(trArray(1:nVoxels))
     diffArray = .false.
     trArray = 0.d0
     tauArray = 0.d0
     call packDiff(octalArray, nDiff, diffArray, tauArray, octalsBelongRank)
     call MPI_ALLREDUCE(diffArray,tArray,nDiff,MPI_LOGICAL,&
         MPI_LOR,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(tauArray,trArray,nDiff,MPI_DOUBLE_PRECISION,&
         MPI_SUM,MPI_COMM_WORLD,ierr)
     diffArray = tArray
     tauArray = trArray
     call unpackDiff(octalArray, nDiff, diffArray, tauArray)
     deallocate(diffArray, tArray, trArray, tauArray)
#endif

    deallocate(octalArray)

  end subroutine setDiffOnTau

#ifdef MPI
      subroutine packDiff(octalArray, nDiff, diffArray, tauArray, octalsBelongRank)
        use mpi
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: octalsBelongRank(:)
        integer :: nDiff
        real(double) :: tauArray(:)
        logical :: diffArray(:)
        integer :: iOctal, iSubcell, my_rank, ierr
        type(OCTAL), pointer :: thisOctal

       call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
       !
       ! Update the bias values of grid computed by all processors.
       !
       nDiff = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
                
             if (.not.thisOctal%hasChild(iSubcell)) then
                 nDiff = nDiff + 1
                 if (octalsBelongRank(iOctal) == my_rank) then
                   diffArray(nDiff) = octalArray(iOctal)%content%diffusionApprox(iSubcell)
                   tauArray(nDiff) = octalArray(iOctal)%content%chiLine(isubcell)
                 else 
                   diffArray(nDiff) = .false.
                 endif
              endif
          
          end do
       end do
     end subroutine packDiff

      subroutine unpackDiff(octalArray, nDiff, diffArray, tauArray)
        type(OCTALWRAPPER) :: octalArray(:)
        integer :: nDiff
        logical :: diffArray(:)
        real(double) :: tauArray(:)
        integer :: iOctal, iSubcell
        type(OCTAL), pointer :: thisOctal

       !
       ! Update the bias values of grid computed by all processors.
       !
       nDiff = 0
       do iOctal = 1, SIZE(octalArray)

          thisOctal => octalArray(iOctal)%content
          
          do iSubcell = 1, thisOctal%maxChildren
                
             if (.not.thisOctal%hasChild(iSubcell)) then
                 nDiff = nDiff + 1
                 octalArray(iOctal)%content%diffusionApprox(iSubcell) = diffArray(nDiff)
                 octalArray(iOctal)%content%chiLine(iSubcell) = tauArray(nDiff)
             endif
          end do
       end do
     end subroutine unpackDiff
#endif

end module diffusion_mod




