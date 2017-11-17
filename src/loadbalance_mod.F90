#ifdef MPI

module loadbalance_mod

  use mpi_global_mod
  use gridio_mod
  use gridtype_mod
  use octal_mod, only: OCTAL
  use source_mod
  implicit none

  integer,pointer :: nLoadBalanceList(:) => null(), loadBalanceList(:,:) => null()
  integer,pointer :: listCounter(:) => null()
contains

  subroutine assignLoadBalancingThreads(grid, frac, method) 
    type(GRIDTYPE) :: grid
    real(double) :: frac(:)
    character(len=*) :: method
    integer :: i, iThread, j, nSolo
    logical, save :: firstTime=.true.
    logical :: writeList

    ! allocate and initilialise arrays
    if (associated(nLoadBalanceList)) then
       deallocate(nloadBalanceList)
       nullify(nloadBalanceList)
    endif
    if (associated(LoadBalanceList)) then
       deallocate(loadBalanceList)
       nullify(loadBalanceList)
    endif
    allocate(nLoadBalanceList(1:nHydroThreadsGlobal))
    nLoadBalanceList = 1

    if (associated(listCounter)) then
       deallocate(listCounter)
       nullify(listCounter)
    endif
    allocate(listCounter(1:nHydroThreadsGlobal))
    listCounter = 1

    allocate(loadBalanceList(1:nHydroThreadsGlobal,1:(nLoadBalancingThreadsGlobal+1)))
    do i = 1, nHydroThreadsGlobal
       loadBalanceList(i,1) = i
    enddo

    ! calculate how many LB threads should be assigned to each hydro thread based on frac array 
    call normaliseLoadBalanceThreads(nHydroThreadsGlobal,  nLoadBalancingThreadsGlobal, nLoadBalanceList, frac)

    ! don't want any lone LB threads 
    if (method == "crossings" .and. any(nLoadBalanceList > 2)) then
       nSolo = 0
       do i = 1, nHydroThreadsGlobal
          if (nLoadBalanceList(i) == 2) then
             nSolo = nSolo + 1
             nLoadBalanceList(i) = 1
          endif
       enddo
       i = 1
       do while (nSolo > 0)
          if (nLoadBalanceList(i) > 2) then
             nLoadBalanceList(i) = nLoadBalanceList(i) + 1
             nSolo = nSolo - 1
          endif
          if (i == nHydroThreadsGlobal) then
             i = 1
          else
             i = i+1
          endif
       enddo
    endif


    ! assign LB rank numbers to hydro threads 
    iThread = nHydroThreadsGlobal+1
    do i = 1, nHydroThreadsGlobal
       if (nLoadBalanceList(i) > 1) then
          do j  = 2, nLoadBalanceList(i)
             if (iThread > nThreadsGlobal-1) then
                write(*,*) "Error assigning load balancing threads ", iThread, myrankglobal 
                stop
             endif
             loadBalanceList(i,j) = iThread
             iThread = iThread + 1
          enddo
       endif
    enddo

    if (writeoutput) then
       write(*,*) "Load balancing thread list (balanced by "//method//")"

       writeList = .true.
       ! if balancing by cells, only need to write out the list the first time
       ! todo - add in condition that mindepth==maxdepth
       if (method == "cells") then
          if (firstTime) then
             firstTime = .false.
          else
             writeList = .false.
          endif
       endif

       if (writeList) then
          do i = 1, nHydroThreadsGlobal
             if (nLoadbalanceList(i) >= 1) &
                  write(*,'(20i4)') i, nLoadBalanceList(i), loadBalanceList(i,1:nLoadBalanceList(i))
          enddo
       else
          write(*,*) "(...)"
       endif
    end if

    ! mpi and grid copying
    call createLoadBalanceCommunicator
    call createLoadThreadDomainCopies(grid)

  end subroutine assignLoadBalancingThreads 

  subroutine setLoadBalancingThreadsBySources(grid)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: i
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer, allocatable :: numberOfSourcesonThread(:), itemp(:)
    integer :: ierr
    real(double), allocatable :: frac(:)

    allocate(numberOfSourcesOnThread(1:nHydroThreadsGlobal))
    numberOfSourcesOnThread = 0 


    if ((.not.loadBalancingThreadGlobal).and.(myrankGlobal /= 0)) then
       do i = 1, globalnSource
          thisOctal => grid%octreeRoot
          call findSubcellTD(globalSourcearray(i)%position, grid%octreeRoot,thisOctal, subcell)
          if (octalOnThread(thisOctal, subcell, myrankGlobal)) then
             numberOfSourcesOnThread(myrankGlobal) = &
                  numberOfSourcesOnThread(thisOctal%mpiThread(subcell)) + 1
          endif
       enddo
    endif
    if ((.not.loadBalancingThreadGlobal).and.(myrankGlobal /= 0)) then
       allocate(itemp(1:nHydroThreadsGlobal))
       call MPI_ALLREDUCE(numberOfSourcesOnThread, itemp, nHydroThreadsGlobal, MPI_INTEGER, MPI_SUM, amrCommunicator, ierr)
       numberOfSourcesOnThread = itemp
       deallocate(itemp)
    endif
    call MPI_BCAST(numberOfSourcesOnThread, nHydroThreadsGlobal, MPI_INTEGER, 1, localWorldCommunicator, ierr)

    allocate(frac(1:nHydroThreadsGlobal))
    frac = dble(numberOfSourcesOnThread(1:nHydroThreadsGlobal))/dble(globalnsource)

    call assignLoadBalancingThreads(grid, frac, "sources")

  end subroutine setLoadBalancingThreadsBySources

  subroutine setLoadBalancingThreadsByLum(grid)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: i
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double), allocatable :: luminosityOnThread(:), itemp(:)
    integer :: ierr
    real(double), allocatable :: frac(:)

    allocate(luminosityOnThread(1:nHydroThreadsGlobal))
    luminosityOnThread = 0.d0


    if ((.not.loadBalancingThreadGlobal).and.(myrankGlobal /= 0)) then
       do i = 1, globalnSource
          thisOctal => grid%octreeRoot
          call findSubcellTD(globalSourcearray(i)%position, grid%octreeRoot,thisOctal, subcell)
          if (octalOnThread(thisOctal, subcell, myrankGlobal)) then
             luminosityOnThread(myrankGlobal) = luminosityOnThread(myrankGlobal) + globalSourceArray(i)%luminosity
          endif
       enddo
       allocate(itemp(1:nHydroThreadsGlobal))
       call MPI_ALLREDUCE(luminosityOnThread, itemp, nHydroThreadsGlobal, & 
               MPI_DOUBLE_PRECISION, MPI_SUM, amrCommunicator, ierr)
       luminosityOnThread = itemp
       deallocate(itemp)
    endif
    call MPI_BCAST(luminosityOnThread, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, 1, localWorldCommunicator, ierr)

    allocate(frac(1:nHydroThreadsGlobal))
    frac = luminosityOnThread(1:nHydroThreadsGlobal)/sum(luminosityOnThread(1:nHydroThreadsGlobal))
    if (writeoutput) write(*,*) "LB total lum: ", sum(luminosityOnThread(1:nHydroThreadsGlobal))/msol

    call assignLoadBalancingThreads(grid, frac, "luminosity")

  end subroutine setLoadBalancingThreadsByLum

  subroutine checkLoadBalance(grid)
    type(GRIDTYPE) :: grid
    integer :: ithread
    integer(bigint) :: j 
    do ithread = 1, nHydroThreadsGlobal+nLoadBalancingThreadsGlobal
       if (iThread == myrankGlobal) then
          j = 0
          call sumCrossings(grid%octreeRoot,j)
          write(*,*) "Sum crossings on rank ",myrankglobal,": ",j
       endif
    enddo
  end subroutine checkLoadBalance
  
  subroutine setLoadBalancingThreadsByCrossings(grid)
    use mpi
    type(GRIDTYPE) :: grid
    integer(bigint), allocatable :: itemp(:)
    integer(bigint), allocatable :: numberOfCrossingsOnThread(:)
    real(double), allocatable :: frac(:)
    integer :: ierr


    allocate(numberOfCrossingsOnThread(1:nHydroThreadsGlobal))
    numberOfCrossingsOnThread = 0
    if ((.not.loadBalancingThreadGlobal).and.(myrankGlobal /=0)) then
       call sumCrossings(grid%octreeRoot, numberOfCrossingsOnThread(myRankGlobal))
       allocate(itemp(1:nHydroThreadsGlobal))
       call MPI_ALLREDUCE(numberOfCrossingsOnThread, itemp, nHydroThreadsGlobal, MPI_INTEGER8, MPI_SUM, amrCommunicator, ierr)
       numberOfCrossingsOnThread = itemp
       deallocate(itemp)
    endif

    call MPI_BCAST(numberOfCrossingsOnThread, nHydroThreadsGlobal, MPI_INTEGER8, 1, localWorldCommunicator, ierr)

    if (writeoutput) write(*,*) "Total number of crossings: ", sum(numberOfCrossingsOnThread)

    ! overflow check
    if (any(numberOfCrossingsOnThread < 0)) then 
       if (writeoutput) write(*,*) "WARNING: negative crossings on thread ", minloc(numberOfCrossingsOnThread), &
          minval(numberOfCrossingsOnThread)
    endif 
    
    ! if there are no crossings at all (e.g. grid hasn't had a photoion loop before)
    if (sum(numberOfCrossingsOnThread) <= 0) then  
       if (writeoutput) write(*,*) "Setting load balancing threads by cells instead."
       call setLoadBalancingThreadsByCells(grid)
       goto 666
    endif

    allocate(frac(1:nHydroThreadsGlobal))
    frac = dble(numberOfCrossingsOnThread)/dble(SUM(numberOfCrossingsOnThread))

    call assignLoadBalancingThreads(grid, frac, "crossings")

    deallocate(frac)
666 continue
    deallocate(numberOfCrossingsOnThread)
  end subroutine setLoadBalancingThreadsByCrossings

  subroutine setLoadBalancingThreadsByScatters(grid)
    use mpi
    type(GRIDTYPE) :: grid
    integer(bigint), allocatable :: itemp(:)
    integer(bigint), allocatable :: numberOfScattersOnThread(:)
    real(double), allocatable :: frac(:)
    integer :: ierr


    allocate(numberOfScattersOnThread(1:nHydroThreadsGlobal))
    numberOfScattersOnThread = 0
    if ((.not.loadBalancingThreadGlobal).and.(myrankGlobal /=0)) then
       call sumScatters(grid%octreeRoot, numberOfScattersOnThread(myRankGlobal))
       allocate(itemp(1:nHydroThreadsGlobal))
       call MPI_ALLREDUCE(numberOfScattersOnThread, itemp, nHydroThreadsGlobal, MPI_INTEGER8, MPI_SUM, amrCommunicator, ierr)
       numberOfScattersOnThread = itemp
       deallocate(itemp)
    endif

    call MPI_BCAST(numberOfScattersOnThread, nHydroThreadsGlobal, MPI_INTEGER8, 1, localWorldCommunicator, ierr)

    if (writeoutput) write(*,*) "Total number of scattering events: ", sum(numberOfScattersOnThread)

    ! overflow check
    if (any(numberOfScattersOnThread < 0)) then 
       if (writeoutput) write(*,*) "WARNING: negative scatters on thread ", minloc(numberOfScattersOnThread), &
          minval(numberOfScattersOnThread)
    endif 
    
    ! if there are no scatters at all (e.g. grid hasn't had a photoion loop before)
    if (sum(numberOfScattersOnThread) <= 0) then  
       if (writeoutput) write(*,*) "Setting load balancing threads by crossings instead."
       call setLoadBalancingThreadsByCrossings(grid)
       goto 777
    endif

    allocate(frac(1:nHydroThreadsGlobal))
    frac = dble(numberOfScattersOnThread)/dble(SUM(numberOfScattersOnThread))

    call assignLoadBalancingThreads(grid, frac, "scatters")

    deallocate(frac)
777 continue
    deallocate(numberOfScattersOnThread)
  end subroutine setLoadBalancingThreadsByScatters

  subroutine setLoadBalancingThreadsByCrossScat(grid)
    use mpi
    type(GRIDTYPE) :: grid
    integer(bigint), allocatable :: itemp(:)
    integer(bigint), allocatable :: numberOfScattersOnThread(:), numberOfCrossingsOnThread(:)
    real(double), allocatable :: frac(:), totalOnThread(:)
    integer :: ierr


    allocate(numberOfScattersOnThread(1:nHydroThreadsGlobal))
    allocate(numberOfCrossingsOnThread(1:nHydroThreadsGlobal))
    numberOfScattersOnThread = 0
    numberOfCrossingsOnThread = 0
    if ((.not.loadBalancingThreadGlobal).and.(myrankGlobal /=0)) then
       call sumScatters(grid%octreeRoot, numberOfScattersOnThread(myRankGlobal))
       call sumCrossings(grid%octreeRoot, numberOfCrossingsOnThread(myRankGlobal))
       allocate(itemp(1:nHydroThreadsGlobal))
       call MPI_ALLREDUCE(numberOfScattersOnThread, itemp, nHydroThreadsGlobal, MPI_INTEGER8, MPI_SUM, amrCommunicator, ierr)
       numberOfScattersOnThread = itemp
       call MPI_ALLREDUCE(numberOfCrossingsOnThread, itemp, nHydroThreadsGlobal, MPI_INTEGER8, MPI_SUM, amrCommunicator, ierr)
       numberOfCrossingsOnThread = itemp
       deallocate(itemp)
    endif

    call MPI_BCAST(numberOfScattersOnThread, nHydroThreadsGlobal, MPI_INTEGER8, 1, localWorldCommunicator, ierr)
    call MPI_BCAST(numberOfCrossingsOnThread, nHydroThreadsGlobal, MPI_INTEGER8, 1, localWorldCommunicator, ierr)

    allocate(totalOnThread(1:nHydroThreadsGlobal))
    totalOnThread = dble(numberOfScattersOnThread) + dble(numberOfCrossingsOnThread)

    if (writeoutput) write(*,*) "Total number of crossings plus scatters: ", sum(totalOnThread)

    ! overflow check
    if (any(numberOfScattersOnThread < 0)) then 
       if (writeoutput) write(*,*) "WARNING: negative scatters on thread ", minloc(numberOfScattersOnThread), &
          minval(numberOfScattersOnThread)
    endif 
    if (any(numberOfCrossingsOnThread < 0)) then 
       if (writeoutput) write(*,*) "WARNING: negative crossings on thread ", minloc(numberOfCrossingsOnThread), &
          minval(numberOfCrossingsOnThread)
    endif 
    
    ! if there are no scatters at all (e.g. grid hasn't had a photoion loop before)
    if (sum(numberOfScattersOnThread) <= 0) then  
       if (sum(numberOfCrossingsOnThread) <= 0) then  
          if (writeoutput) write(*,*) "Setting load balancing threads by cells instead."
          call setLoadBalancingThreadsByCells(grid)
       else
          if (writeoutput) write(*,*) "Setting load balancing threads by crossings instead."
          call setLoadBalancingThreadsByCrossings(grid)
       endif
       goto 777
    endif


    allocate(frac(1:nHydroThreadsGlobal))
    frac = totalOnThread/SUM(totalOnThread)

    call assignLoadBalancingThreads(grid, frac, "crossplusscat")

    deallocate(frac)
777 continue
    deallocate(numberOfScattersOnThread)
    deallocate(numberOfCrossingsOnThread)
  end subroutine setLoadBalancingThreadsByCrossScat

  subroutine setLoadBalancingThreadsByCells(grid)
    use mpi
    type(GRIDTYPE) :: grid
    integer, allocatable :: itemp(:)
    integer, allocatable :: numberOfCellsOnThread(:)
    real(double), allocatable :: frac(:)
    integer :: ierr

    allocate(numberOfCellsOnThread(1:nHydroThreadsGlobal))
    allocate(frac(1:nHydroThreadsGlobal))


    numberOfCellsOnThread = 0
    if ((.not.loadBalancingThreadGlobal).and.(myrankGlobal /=0)) then
       call sumCellsOnThread(grid%octreeRoot, numberOfCellsOnThread(myRankGlobal))
       allocate(itemp(1:nHydroThreadsGlobal))
       call MPI_ALLREDUCE(numberOfCellsOnThread, itemp, nHydroThreadsGlobal, MPI_INTEGER, MPI_SUM, amrCommunicator, ierr)
       numberOfCellsOnThread = itemp
       deallocate(itemp)
    endif


    call MPI_BCAST(numberOfCellsOnThread, nHydroThreadsGlobal, MPI_INTEGER, 1, localWorldCommunicator, ierr)
    frac = dble(numberOfCellsOnThread)/dble(SUM(numberOfCellsOnThread))

    call assignLoadBalancingThreads(grid, frac, "cells")

    deallocate(numberOfCellsOnThread)
    deallocate(frac)

  end subroutine setLoadBalancingThreadsByCells
       

  subroutine createLoadThreadDomainCopies(grid)
    use mpi
    use timing
    type(GRIDTYPE) :: grid
    integer :: i, ierr

    do i = 1, nHydroThreadsGlobal
       if (nLoadBalanceList(i) > 1) then
          if (any(loadBalanceList(i,1:nLoadBalanceList(i)) == myrankGlobal)) then
             call broadcastBranch(grid, loadBalanceList(i,1), loadBalanceCommunicator(i))
          endif
       endif
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  end subroutine createLoadThreadDomainCopies

  subroutine copyDomain(fromThread, toThread, grid)
    integer :: fromThread, toThread
    type(GRIDTYPE) :: grid

    if (myrankGlobal == toThread) then
       write(*,*) myrankGlobal, " getting domain from ", fromThread
       call deleteOctreeBranch(thisOctal=grid%octreeroot,                      &
               onlyChildren=.FALSE.,                                           &
               adjustParent=.FALSE.)

       call getBranchOverMpi(grid%OctreeRoot, null(), fromThread)
       copyOfThread = fromThread
       write(*,*) myrankGlobal, " done"
    endif
    if (myrankGlobal == fromThread) then
       write(*,*) myrankGlobal, " sending domain to ",toThread
       call sendBranchOverMpi(grid%octreeRoot, toThread)
       write(*,*) myrankGlobal, "done"
    endif
  end subroutine copyDomain

  recursive subroutine sendBranchOverMpi(thisoctal, iThread)
    use mpi
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: i, iThread

    call sendOctalViaMPI(thisOctal, iThread)

    IF ( thisOctal%nChildren > 0 ) THEN
       ! call this subroutine recursively on each of its children
       DO i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          CALL sendBranchOverMpi(child, iThread)
       END DO
    END IF


  end subroutine sendBranchOverMpi

  function loadBalancedThreadNumber(iThread) result(out)
    integer :: iThread, out

    if (.not.associated(nLoadBalanceList)) then
       out = iThread
       goto 666
    endif
    if (nLoadBalanceList(iThread) == 1) then
       out = ithread
    else
       out = loadBalanceList(iThread, listCounter(ithread))
       listCounter(iThread) = listCounter(ithread) + 1
       if (listCounter(iThread) > nLoadBalanceList(iThread)) then
          listCounter(iThread) = 1
       endif
    endif
666 continue
  end function loadBalancedThreadNumber


  subroutine createLoadBalanceCommunicator

    use mpi
    integer :: iThread, j
    integer, allocatable :: ranks(:)
    integer :: worldGroup, ierr, loadBalanceThreadGroup

    if (allocated(loadBalanceCommunicator)) then
       do iThread = 1, nHydroThreadsGlobal
          if (loadBalanceCommunicator(iThread) /= MPI_COMM_NULL) then
             call MPI_COMM_FREE(loadBalanceCommunicator(iThread), ierr)
          endif
       enddo
       deallocate(loadBalanceCommunicator)
    endif

    allocate(loadBalanceCommunicator(1:nHydroThreadsGlobal))
    
    call MPI_COMM_GROUP(MPI_COMM_WORLD, worldGroup, ierr)



    do iThread = 1, nHydroThreadsGlobal
       allocate(ranks(1:nLoadBalanceList(iThread)))
       do j = 1, nLoadBalanceList(iThread)
          ranks(j) = loadBalanceList(iThread,j)
       enddo
       call MPI_GROUP_INCL(worldGroup, nLoadBalanceList(iThread), ranks, loadBalanceThreadGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, loadBalanceThreadGroup, loadBalanceCommunicator(iThread), ierr)
       call MPI_GROUP_FREE(loadBalanceThreadGroup, ierr)
       deallocate(ranks)
    enddo
    call MPI_GROUP_FREE(worldGroup, ierr)


  end subroutine createLoadBalanceCommunicator

recursive subroutine sumDensity(thisOctal, density)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  real(double) :: density
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call sumDensity(child, density)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

        density = density + thisOctal%rho(subcell)
     end if
  end do
end subroutine sumDensity

recursive subroutine sumCrossings(thisOctal, n)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  integer(bigint) ::  n
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call sumCrossings(child, n)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

        if (thisOctal%nCrossings(subcell) < 0) then
           write(*,*) myrankglobal, " %ncrossings(subcell) negative ", thisOctal%nCrossings(subcell)
           write(*,*) myrankglobal, " negative location ", subcellCentre(thisOctal, subcell)
        endif
        n = n + thisOctal%nCrossings(subcell)
     end if
  end do
end subroutine sumCrossings

recursive subroutine sumScatters(thisOctal, n)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  integer(bigint) ::  n
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call sumScatters(child, n)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

        if (thisOctal%nScatters(subcell) < 0) then
           write(*,*) myrankglobal, " %nScatters(subcell) negative ", thisOctal%nScatters(subcell)
           write(*,*) myrankglobal, " negative location ", subcellCentre(thisOctal, subcell)
        endif
        n = n + thisOctal%nScatters(subcell)
     end if
  end do
end subroutine sumScatters

recursive subroutine sumCellsOnThread(thisOctal, n)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell, n
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call sumCellsOnThread(child, n)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
        n = n + 1
     end if
  end do
end subroutine sumCellsOnThread

subroutine normaliseLoadBalanceThreads(nHydroThreadsGlobal, nLoadBalancingThreadsGlobal, nLoadBalanceList, frac)
  integer :: nHydroThreadsGlobal, nLoadBalancingThreadsGlobal
  real(double) :: frac(:)
  integer :: nLoadBalanceList(:)
  integer :: iter
  logical, allocatable :: thisMask(:)

  nLoadBalanceList(1:nHydroThreadsGlobal) = 1 + &
       nint(dble(nLoadBalancingThreadsGlobal) * frac(1:nHydroThreadsGlobal))

  iter = 0
  do while ((SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal) /= nLoadBalancingThreadsGlobal)

     if ((SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal) > nLoadBalancingThreadsGlobal) then
        frac = frac / 1.0001d0
     else
        frac = frac * 1.0001d0
     endif
     nLoadBalanceList(1:nHydroThreadsGlobal) = 1 + &
          nint(dble(nLoadBalancingThreadsGlobal) * frac(1:nHydroThreadsGlobal))
     iter = iter + 1
     if (iter > 100000) exit
  end do


  allocate(thisMask(1:nHydroThreadsGlobal))
  thisMask = .true.
  do while ((SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal) /= nLoadBalancingThreadsGlobal)
     
     if ((SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal) > nLoadBalancingThreadsGlobal) then
        nLoadBalanceList(MAXLOC(nLoadBalanceList)) = nLoadBalanceList(MAXLOC(nLoadBalanceList)) - 1
     else
        nLoadBalanceList(MAXLOC(nLoadBalanceList,mask=thisMask)) = nLoadBalanceList(MAXLOC(nLoadBalanceList,mask=thisMask)) + 1
        thismask(MAXLOC(nLoadBalanceList,mask=thisMask)) = .false.
     endif
  end do
  deallocate(thisMask)

end subroutine normaliseLoadBalanceThreads
       
end module loadbalance_mod


#else

! Dummy module for non-MPI builds

module loadbalance_mod

end module loadbalance_mod

#endif
