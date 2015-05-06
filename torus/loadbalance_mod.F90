#ifdef MPI

module loadbalance_mod

  use mpi_global_mod
  use gridio_mod
  use gridtype_mod
  use octal_mod
  use source_mod
  implicit none

  integer,pointer :: nLoadBalanceList(:) => null(), loadBalanceList(:,:) => null()
  integer,pointer :: listCounter(:) => null()
contains


  subroutine setLoadBalancingThreadsBySources(grid)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: i, iThread, j
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer, allocatable :: numberOfSourcesonThread(:), itemp(:)
    integer :: ierr

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

    nLoadBalanceList(1:nHydroThreadsGlobal) = nLoadBalanceList(1:nHydroThreadsGlobal) + &
         int(dble(nLoadBalancingThreadsGlobal)*dble(numberOfSourcesOnThread(1:nHydroThreadsGlobal))/dble(globalnSource))


    do while ((SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal) /= nLoadBalancingThreadsGlobal)

       if ((SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal) > nLoadBalancingThreadsGlobal) then
          nLoadBalanceList(MAXLOC(nLoadBalanceList)) = nLoadBalanceList(MAXLOC(nLoadBalanceList)) - 1
       else
          nLoadBalanceList(MAXLOC(nLoadBalanceList)) = nLoadBalanceList(MAXLOC(nLoadBalanceList)) + 1
       endif
    end do

    iThread = nHydroThreadsGlobal+1
    do i = 1, nHydroThreadsGlobal
       if (nLoadBalanceList(i) > 1) then
          do j  = 2, nLoadBalanceList(i)
             if (iThread > nThreadsGlobal-1) then
                write(*,*) "Error assigning load balancing threads ",ithread
             endif
             loadBalanceList(i,j) = iThread
             iThread = iThread + 1
          enddo
       endif
    enddo

    if (writeoutput) then
       write(*,*) "Load balancing thread list"
       do i = 1, nHydroThreadsGlobal
          if (nLoadbalanceList(i) > 1) &
               write(*,'(20i4)') i, nLoadBalanceList(i), loadBalanceList(i,1:nLoadBalanceList(i))
       enddo
    end if

    call createLoadBalanceCommunicator
    call createLoadBalancingThreadDomainCopies(grid)

  end subroutine setLoadBalancingThreadsBySources

  subroutine checkLoadBalance(grid)
    type(GRIDTYPE) :: grid
    integer :: ithread,j 
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
    use utils_mod, only : median
    type(GRIDTYPE) :: grid
    integer :: i, iThread, j
    integer :: subcell
    integer, allocatable :: itemp(:)
    integer, allocatable :: numberOfCrossingsOnThread(:)
    real(double), allocatable :: frac(:)
    integer :: ierr, iter

    allocate(numberOfCrossingsOnThread(1:nHydroThreadsGlobal))
    allocate(frac(1:nHydroThreadsGlobal))


    numberOfCrossingsOnThread = 0
    if ((.not.loadBalancingThreadGlobal).and.(myrankGlobal /=0)) then
       call sumCrossings(grid%octreeRoot, numberOfCrossingsOnThread(myRankGlobal))
       allocate(itemp(1:nHydroThreadsGlobal))
       call MPI_ALLREDUCE(numberOfCrossingsOnThread, itemp, nHydroThreadsGlobal, MPI_INTEGER, MPI_SUM, amrCommunicator, ierr)
       numberOfCrossingsOnThread = itemp
       deallocate(itemp)
    endif


    call MPI_BCAST(numberOfCrossingsOnThread, nHydroThreadsGlobal, MPI_INTEGER, 1, localWorldCommunicator, ierr)
    frac = dble(numberOfCrossingsOnThread)/dble(SUM(numberOfCrossingsOnThread))

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
    nLoadBalanceList(1:nHydroThreadsGlobal) =  1 + &
         nint(dble(nLoadBalancingThreadsGlobal) * frac(1:nHydroThreadsGlobal))


    iter = 0
    do while ((SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal) /= nLoadBalancingThreadsGlobal)


       if ((SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal) > nLoadBalancingThreadsGlobal) then
          frac = frac / 1.001d0
       else
          frac = frac * 1.001d0
       endif
       nLoadBalanceList(1:nHydroThreadsGlobal) = 1 + &
            nint(dble(nLoadBalancingThreadsGlobal) * frac(1:nHydroThreadsGlobal))
       write(*,*) iter, (SUM(nLoadBalanceList(1:nHydroThreadsGlobal)) - nHydroThreadsGlobal), nLoadBalancingThreadsGlobal
       iter = iter + 1
       if (iter > 10000) exit
    end do


    iThread = nHydroThreadsGlobal+1
    do i = 1, nHydroThreadsGlobal
       if (nLoadBalanceList(i) > 1) then
          do j  = 2, nLoadBalanceList(i)
             if (iThread > nThreadsGlobal-1) then
                write(*,*) "Error assigning load balancing threads"
             endif
             loadBalanceList(i,j) = iThread
             iThread = iThread + 1
          enddo
       endif
    enddo

    if (writeoutput) then
       write(*,*) "Load balancing thread list"
       do i = 1, nHydroThreadsGlobal
          if (nLoadbalanceList(i) > 1) &
               write(*,'(20i4)') i, nLoadBalanceList(i), loadBalanceList(i,1:nLoadBalanceList(i))
       enddo
    end if

    call createLoadBalanceCommunicator
    call createLoadBalancingThreadDomainCopies(grid)

  end subroutine setLoadBalancingThreadsByCrossings
       
  subroutine createLoadBalancingThreadDomainCopies(grid)
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

  end subroutine createLoadBalancingThreadDomainCopies

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

    if (nLoadBalanceList(iThread) == 1) then
       out = ithread
    else
       out = loadBalanceList(iThread, listCounter(ithread))
       listCounter(iThread) = listCounter(ithread) + 1
       if (listCounter(iThread) > nLoadBalanceList(iThread)) then
          listCounter(iThread) = 1
       endif
    endif
  end function loadBalancedThreadNumber


  subroutine createLoadBalanceCommunicator

    use mpi
    integer :: iThread, j
    integer, allocatable :: ranks(:)
    integer :: worldGroup, ierr, loadBalanceThreadGroup

    if (allocated(loadBalanceCommunicator)) deallocate(loadBalanceCommunicator)

    allocate(loadBalanceCommunicator(1:nHydroThreadsGlobal))
    
    call MPI_COMM_GROUP(MPI_COMM_WORLD, worldGroup, ierr)


    do iThread = 1, nHydroThreadsGlobal
       allocate(ranks(1:nLoadBalanceList(iThread)))
       do j = 1, nLoadBalanceList(iThread)
          ranks(j) = loadBalanceList(iThread,j)
       enddo
       call MPI_GROUP_INCL(worldGroup, nLoadBalanceList(iThread), ranks, loadBalanceThreadGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, loadBalanceThreadGroup, loadBalanceCommunicator(iThread), ierr)
       deallocate(ranks)
    enddo


  end subroutine createLoadBalanceCommunicator

recursive subroutine sumCrossings(thisOctal, n)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell, n
  
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
        n = n + thisOctal%nCrossings(subcell)
     end if
  end do
end subroutine sumCrossings
       

end module loadbalance_mod

#else

! Dummy module for non-MPI builds

module loadbalance_mod

end module loadbalance_mod

#endif
