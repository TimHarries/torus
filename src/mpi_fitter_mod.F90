module mpi_fitter_mod

#ifdef MPI

  
  use kind_mod
  use mpi_global_mod
  implicit none
  
  contains
  
  subroutine setupFitterCommunicators()
    use mpi
    integer :: i, ierr
    integer :: nThreadsPerFitter
    integer :: worldGroup, localWorldGroup
    integer, allocatable :: ranks(:)
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreadsGlobal, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRankWorldGlobal, ierr)

    call MPI_COMM_GROUP(MPI_COMM_WORLD, worldGroup, ierr)

    nThreadsPerFitter = 3

    myFitterSetGlobal = (myRankGlobal-1)/nThreadsPerFitter

    

    
    if (myrankGlobal /= 0) then
       allocate(ranks(1:nThreadsPerFitter))
       do i = 1, nThreadsPerFitter
          ranks(i) = myFitterSetGlobal * nThreadsPerFitter + i 
       enddo
       write(*,*) myrankglobal, " set ",myFitterSetGlobal," ranks ", ranks
       call MPI_GROUP_INCL(worldGroup, nThreadsPerFitter, ranks, localWorldGroup, ierr)
          call MPI_COMM_CREATE(MPI_COMM_WORLD, localWorldGroup, localWorldCOMMUNICATOR, ierr)
       deallocate(ranks)
    endif
       
    if (myrankGlobal == 0) then
       allocate(ranks(1:1))
       ranks(1) = 0
       call MPI_GROUP_INCL(worldGroup, 1, ranks, localWorldGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, localWorldGroup, localWorldCOMMUNICATOR, ierr)
       deallocate(ranks)
    endif

  end subroutine setupFitterCommunicators


  
#endif
    
end module mpi_fitter_mod
