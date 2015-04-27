#ifdef MPI

module loadbalance_mod

  use mpi_global_mod
  use gridio_mod
  use gridtype_mod
  use octal_mod
  use source_mod


  integer,pointer :: nLoadBalanceList(:) => null(), loadBalanceList(:,:) => null()
  integer,pointer :: listCounter(:) => null()
contains


  subroutine setLoadBalancingThreadsBySources(grid)
    type(GRIDTYPE) :: grid
    integer :: i, iThread, j
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer, allocatable :: numberOfSourcesonThread(:)

    allocate(numberOfSourcesOnThread(1:nHydroThreadsGlobal))
    numberOfSourcesOnThread = 0 


    do i = 1, globalnSource
       thisOctal => grid%octreeRoot
       call findSubcellTD(globalSourcearray(i)%position, grid%octreeRoot,thisOctal, subcell)
       numberOfSourcesOnThread(thisOctal%mpiThread(subcell)) = &
            numberOfSourcesOnThread(thisOctal%mpiThread(subcell)) + 1
    enddo
    

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

    write(*,*) "nLoadBalanceList ",nLoadBalanceList
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
          write(*,'(10i4)') i, nLoadBalanceList(i), loadBalanceList(i,1:nLoadBalanceList(i))
       enddo
    end if

    call createLoadBalancingThreadDomainCopies(grid)

  end subroutine setLoadBalancingThreadsBySources
       
  subroutine createLoadBalancingThreadDomainCopies(grid)
    type(GRIDTYPE) :: grid
    integer :: i, j

    do i = 1, nHydroThreadsGlobal
       if (nLoadBalanceList(i) > 1) then
          do j = 2, nLoadBalanceList(i)
             call copyDomain(i, loadBalanceList(i,j), grid)
          enddo
       endif
    enddo
  end subroutine createLoadBalancingThreadDomainCopies

  subroutine copyDomain(fromThread, toThread, grid)
    integer :: fromThread, toThread
    type(GRIDTYPE) :: grid

    if (myrankGlobal == toThread) then
       write(*,*) myrankGlobal, " getting domain from ", fromThread
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


end module loadbalance_mod

#else

! Dummy module for non-MPI builds

module loadbalance_mod

end module loadbalance_mod

#endif
