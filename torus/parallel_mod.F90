module parallel_mod

! subprograms used for parallel computation

  use kind_mod
  use mpi_global_mod
  use messages_mod
  use utils_mod
  implicit none

contains  

#ifdef MPI
  subroutine mpiBlockHandout(nProc,unitBelongsRank,&
                             blockDivFactor,tag,&
                             minBlockSize,maxBlockSize,setDebug)
    ! this routine will divide a task into many small blocks, and
    !   hand them out to processes that request work.
  
    use input_variables, only: blockhandout
    implicit none
    include 'mpif.h'  
  
    integer, intent(in) :: nProc    ! the number of processes
    integer, intent(in) :: blockDivFactor ! used to set block size
    integer, dimension(:), intent(inout) :: unitBelongsRank
    integer, intent(in) :: tag ! MPI tag (probably zero)
    integer,intent(in),optional :: minBlockSize
    integer,intent(in),optional :: maxBlockSize
    logical,intent(in),optional :: setDebug ! print debug messages
  
    integer :: blockSize
    integer :: iErr     ! error flag
    integer :: nUnits
    integer :: startUnit, endUnit
    integer :: iRank
    logical, dimension(:), allocatable :: rankComplete
    logical :: debugInfo
    integer, dimension(MPI_STATUS_SIZE) :: mpiStatus    
    !character*(MPI_MAX_ERROR_STRING) :: errorString
        !    ^ outdated, but required to work with MPI wrapper 
    !integer :: errorLength
    integer :: count


    if (present(setDebug)) then
      debugInfo = setDebug
    else
      debugInfo = .false.
    end if
    
    if (present(minBlockSize)) then
      if (minBlockSize > nUnits) then
        print *, 'Warning: in ''mpiBlockHandout'', minBlockSize > nUnits'
      end if
    end if
    
    unitBelongsRank = -999
    nUnits = SIZE(unitBelongsRank)
    startUnit = 1
    
    blockSize = nUnits / (nProc * blockDivFactor)

    ! force blockSize to be within sensible limits
    blockSize = MAX(1,blockSize)
    if (present(minBlockSize)) blockSize = MAX(minBlockSize,blockSize) 
    if (present(maxBlockSize)) blockSize = MIN(maxBlockSize,blockSize)
  !  print *, 'Work unit blocksize = ', blockSize

    ! th added this

    if (.not.BlockHandout) then
          blockSize = nUnits / (nProc-1)
          iRank = 0
    endif
    endUnit = startUnit + (blockSize - 1)
       
    allocate(rankComplete(nProc-1)) ! we don't need a status value for rank 0
    rankComplete = .false.
       
    do 
      if (debugInfo) print *, 'At start of DO loop: ', startUnit, endUnit
      
      ! wait until we get sent the rank of a process needing a new block

    ! th added if statement

    if (blockHandout) then
      if (debugInfo) print *, 'Waiting for a block request...'
      call MPI_RECV(iRank,1,MPI_INTEGER,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,mpiStatus,iErr)
    else
      iErr = 0
      iRank = iRank + 1
      if (iRank .eq. (nProc-1)) endUnit=nUnits+1 ! last thread finishes rest
      if (iRank .eq. nProc) exit
    endif


      if (debugInfo .and. iErr /= 0) then
        print *, 'MPI return value = ',iErr
        print *, 'source = ', mpiStatus(MPI_SOURCE)
        print *, 'tag = ', mpiStatus(MPI_TAG)
        print *, 'error = ', mpiStatus(MPI_ERROR)
        !call MPI_ERROR_STRING(mpiStatus(MPI_ERROR),errorString,errorLength,iErr)
        !print *, 'errorString = ', errorString(1:errorLength), iErr
        call MPI_GET_COUNT(mpiStatus,MPI_INTEGER,count,iErr)
        print *, 'count = ',count, iErr
      end if
      
      if (debugInfo) print *, 'Received request from rank ', iRank   
      
      ! check whether we have ran out of blocks to give out
      if (startUnit > nUnits) then
        if (debugInfo) print *,'Have ran out of units', iRank
        
        call MPI_SEND(-999,1,MPI_INTEGER,iRank,tag,MPI_COMM_WORLD,iErr)
        if (debugInfo .and. iErr /= 0) print *, 'MPI return value = ',iErr
        call MPI_SEND(-999,1,MPI_INTEGER,iRank,tag,MPI_COMM_WORLD,iErr)
        if (debugInfo .and. iErr /= 0) print *, 'MPI return value = ',iErr
        rankComplete(iRank) = .true.
        
      else if (endUnit > nUnits) then
        ! this is the last block we have left
        if (debugInfo) print *,' (endUnit > nUnits) ', iRank, startUnit, nUnits
        
        endUnit = nUnits
        
        ! sanity check
        if (endUnit > nUnits .or. startUnit < 1) then
           print *, 'Trying to send an invalid work unit!' 
           print *, startUnit, endUnit, nUnits     
           stop
        end if
        
        if (ANY(unitBelongsRank(startUnit:endUnit) /= -999)) then
          print *, 'Trying to give out an already allocated block!'
          stop
        end if
        
        unitBelongsRank(startUnit:nUnits) = iRank
        call MPI_SEND(startUnit,1,MPI_INTEGER,iRank,tag,MPI_COMM_WORLD,iErr)
        if (debugInfo .and. iErr /= 0) print *, 'MPI return value = ',iErr
        call MPI_SEND(endUnit,1,MPI_INTEGER,iRank,tag,MPI_COMM_WORLD,iErr)
        if (debugInfo .and. iErr /= 0) print *, 'MPI return value = ',iErr
        
      else 
        ! we have (at least) one complete block still to give out
        if (debugInfo) print *,'Normal allocation:', iRank, startUnit, endUnit
        
        ! sanity check
        if (endUnit > nUnits .or. startUnit < 1) then
           print *, 'Trying to send an invalid work unit!' 
           print *, startUnit, endUnit, nUnits     
           stop
        end if
        
        if (ANY(unitBelongsRank(startUnit:endUnit) /= -999)) then
          print *, 'Trying to give out an already allocated block!'
          stop
        end if
        
        unitBelongsRank(startUnit:endUnit) = iRank
        call MPI_SEND(startUnit,1,MPI_INTEGER,iRank,tag,MPI_COMM_WORLD,iErr)
        if (debugInfo .and. iErr /= 0) print *, 'MPI return value = ',iErr
        call MPI_SEND(endUnit,1,MPI_INTEGER,iRank,tag,MPI_COMM_WORLD,iErr)
        if (debugInfo .and. iErr /= 0) print *, 'MPI return value = ',iErr
      end if  
    
      ! check if we're all done
      if (ALL(rankComplete)) exit  
    
      ! increment the counter variables
      startUnit = startUnit + blockSize 
      endUnit = endUnit + blockSize 
    
    end do

    deallocate(rankComplete)
  
  end subroutine mpiBlockHandout

  subroutine mpiGetBlock(myRank,startUnit,endUnit,allDone,tag,setDebug)
    ! this requests a work unit from the root process (which runs mpiBlockHandout)
  
    use input_variables, only: blockhandout
    implicit none
    include 'mpif.h'  
  
    integer, intent(in) :: myRank
    integer, intent(out) :: startUnit, endUnit
    logical, intent(out) :: allDone
    integer, intent(in) :: tag ! MPI tag (probably zero)
    logical, intent(in), optional :: setDebug
  
    integer :: iErr     ! error flag
    logical :: debugInfo
    integer, dimension(MPI_STATUS_SIZE) :: mpiStatus    
    !character, dimension(MPI_MAX_ERROR_STRING) :: errorString
    !character*(MPI_MAX_ERROR_STRING) :: errorString
    !integer :: errorLength
    integer :: count


    if (present(setDebug)) then
      debugInfo = setDebug
    else
      debugInfo = .false.
    end if
    
    ! send 'my_rank' to indicate we need another block 
    if (debugInfo) print *, 'Rank ', myRank, ' requesting a block'

    ! th added if statement

    iErr = 0
    if (blockHandout) then
    call MPI_SEND(myRank,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,ierr)
    if (debugInfo .and. iErr /= 0) print *, 'MPI return value = ',iErr
    if (debugInfo) print *, 'Rank ', myRank, ' request sent'
    endif
    
    ! wait until we receive the block boundaries
    call MPI_RECV(startUnit,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,mpiStatus,ierr)

      if (debugInfo .and. iErr /= 0) then
        print *, 'MPI return value = ',iErr
        print *, 'source = ', mpiStatus(MPI_SOURCE)
        print *, 'tag = ', mpiStatus(MPI_TAG)
        print *, 'error = ', mpiStatus(MPI_ERROR)
        !call MPI_ERROR_STRING(mpiStatus(MPI_ERROR),errorString,errorLength,iErr)
        !print *, 'errorString = ', errorString(1:errorLength), iErr
        call MPI_GET_COUNT(mpiStatus,MPI_INTEGER,count,iErr)
        print *, 'count = ',count, iErr
      end if

    call MPI_RECV(endUnit,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,mpiStatus,ierr)

      if (debugInfo .and. iErr /= 0) then
        print *, 'MPI return value = ',iErr
        print *, 'source = ', mpiStatus(MPI_SOURCE)
        print *, 'tag = ', mpiStatus(MPI_TAG)
        print *, 'error = ', mpiStatus(MPI_ERROR)
        !call MPI_ERROR_STRING(mpiStatus(MPI_ERROR),errorString,errorLength,iErr)
        !print *, 'errorString = ', errorString(1:errorLength), iErr
        call MPI_GET_COUNT(mpiStatus,MPI_INTEGER,count,iErr)
        print *, 'count = ',count, iErr
      end if
    
    if (debugInfo) print *, 'Rank ', myRank, ' received ',startUnit, endUnit
    
    if (startUnit == -999) then
      allDone = .true. 
    else
      allDone = .false. 
    end if

  end subroutine mpiGetBlock


!-------------------------------------------------------------------------------

  subroutine torus_mpi_barrier(message)

    USE mpi_global_mod, ONLY: myRankGlobal

    implicit none
    include 'mpif.h'

    character(len=*), optional, intent(in) :: message 
    integer :: ierr

    if ( present(message) ) then
       print *, "Process ", myRankGlobal, ": ", message
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

  end subroutine torus_mpi_barrier


! End of MPI routines ----------------------------------------------------------

#else

!-------------------------------------------------------------------------------

  subroutine torus_mpi_barrier(message)

    implicit none

    character(len=*), optional, intent(in) :: message 
    character(len=1) :: test
    if (PRESENT(message)) then
       test(1:1) = message(1:1)
    endif
  end subroutine torus_mpi_barrier

#endif

!-------------------------------------------------------------------------------

! Abort Torus. If the message argument is present this is used as the 
! error message, otherwise a generic abort message is given. 
!
! D. Acreman, September 2008

  subroutine torus_abort(message)

    implicit none

#ifdef MPI
    include 'mpif.h'
    integer :: ierr
#endif

    character(len=*), optional, intent(in) :: message

    if ( present(message) ) then
       write(*,*) message
    else
       write(*,*) "TORUS aborting"
    end if

#ifdef MPI
    call mpi_abort(MPI_COMM_WORLD, ierr)
#else
    STOP
#endif

  end subroutine torus_abort

end module parallel_mod


!!! vim:set filetype=fortran :                                !!!
!!! otherwise vim won't recognize a file with the suffix .raw !!!

!End of file -------------------------------------------------------------------
