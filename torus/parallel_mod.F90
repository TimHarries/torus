module parallel_mod

! subprograms used for parallel computation

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
! Name:    gather_sph_data
! Purpose: Perform the mpi communication required to ensure all processes have
!          access to the full list of sph particle data
! Author:  D. Acreman October 2007

  subroutine gather_sph_data

    USE kind_mod
    USE sph_data_class, only: sphData

    implicit none
    include 'mpif.h'

    integer :: my_rank, n_proc, ierr, i 
    integer :: npart_all, nptmass_all
    real(double)    :: totalgasmass_all
    integer, allocatable :: npart_arr(:),   nptmass_arr(:)
    integer, allocatable :: npart_displ(:), nptmass_displ(:)
    real(double), allocatable :: xn_tmp(:), yn_tmp(:), zn_tmp(:), rhon_tmp(:), temperature_tmp(:)
    real(double), allocatable :: gasmass_tmp(:), hn_tmp(:)
    real(double), allocatable :: x_tmp(:),  y_tmp(:),  z_tmp(:),  ptmass_tmp(:)
    character(len=3) :: char_my_rank
    logical, parameter :: ll_testwrite = .false.

! Begin executable statements -------

! 0. Preliminaries
  !  Get my process rank  
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
  ! Find the total number of processes being used in this run
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr)

! 1.1 Get total number of gas particles
  call MPI_ALLREDUCE(sphData%npart, npart_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

! 1.2 Get number of gas particles for each process
  ALLOCATE( npart_arr(n_proc) )
  CALL MPI_GATHER(sphData%npart, 1, MPI_INTEGER, npart_arr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

! 1.3 Get total gas mass
  call MPI_ALLREDUCE(sphData%totalgasmass, totalgasmass_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr) 
  sphData%totalgasmass = totalgasmass_all

! 2.1 Get total number of point masses
  call MPI_ALLREDUCE(sphData%nptmass, nptmass_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

! 2.2 Get total number of point masses for each proccess
  ALLOCATE( nptmass_arr(n_proc) )
  CALL MPI_GATHER(sphData%nptmass, 1, MPI_INTEGER, nptmass_arr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

! 3. Communicate arrays xn, yn, zn, rhon, x, y, z, ptmass and keep in tmp storage

! 3.1 Gas particles
  ALLOCATE ( npart_displ(n_proc) )
  npart_displ(1)=0
  DO i=2, n_proc
     npart_displ(i) = npart_displ(i-1) + npart_arr(i-1)
  END DO

  ALLOCATE ( xn_tmp(npart_all)          )
  ALLOCATE ( yn_tmp(npart_all)          )
  ALLOCATE ( zn_tmp(npart_all)          )
  ALLOCATE ( rhon_tmp(npart_all)        )
  ALLOCATE ( temperature_tmp(npart_all) )
  ALLOCATE ( gasmass_tmp(npart_all)     )
  ALLOCATE ( hn_tmp(npart_all)          )

! Gather the particle data on process 0 then broadcast to all. For some reason MPI_ALLGATHERV 
! doesn't work  for doing this. 
  CALL MPI_GATHERV(sphData%xn, sphData%npart, MPI_DOUBLE_PRECISION, xn_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(xn_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%yn, sphData%npart, MPI_DOUBLE_PRECISION, yn_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(yn_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%zn, sphData%npart, MPI_DOUBLE_PRECISION, zn_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(zn_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%rhon, sphData%npart, MPI_DOUBLE_PRECISION, rhon_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(rhon_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%temperature, sphData%npart, MPI_DOUBLE_PRECISION, temperature_tmp(:), npart_arr(:), &
       npart_displ(:), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(temperature_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%gasmass, sphData%npart, MPI_DOUBLE_PRECISION, gasmass_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(gasmass_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%hn, sphData%npart, MPI_DOUBLE_PRECISION, hn_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(hn_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

! 3.2 Point masses
  IF ( nptmass_all > 0 ) THEN

     ALLOCATE ( nptmass_displ(n_proc) )
     nptmass_displ(1)=0
     DO i=2, n_proc
        nptmass_displ(i) = nptmass_displ(i-1) + nptmass_arr(i-1)
     END DO

     ALLOCATE ( x_tmp(npart_all)      )
     ALLOCATE ( y_tmp(npart_all)      )
     ALLOCATE ( z_tmp(npart_all)      )
     ALLOCATE ( ptmass_tmp(npart_all) )

     CALL MPI_GATHERV(sphData%x, sphData%nptmass, MPI_DOUBLE_PRECISION, x_tmp(:), nptmass_arr(:), nptmass_displ(:), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST(x_tmp(:), nptmass_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     
     CALL MPI_GATHERV(sphData%y, sphData%nptmass, MPI_DOUBLE_PRECISION, y_tmp(:), nptmass_arr(:), nptmass_displ(:), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST(y_tmp(:), nptmass_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

     CALL MPI_GATHERV(sphData%z, sphData%nptmass, MPI_DOUBLE_PRECISION, z_tmp(:), nptmass_arr(:), nptmass_displ(:), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST(z_tmp(:), nptmass_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

     CALL MPI_GATHERV(sphData%ptmass, sphData%nptmass, MPI_DOUBLE_PRECISION, ptmass_tmp(:), nptmass_arr(:), nptmass_displ(:), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST(ptmass_tmp(:), nptmass_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  END IF

! Write out particle list for tests
  IF ( ll_testwrite ) THEN
     write(char_my_rank, '(i3)') my_rank
     open (unit=60, status='replace', file='mpi_test_tmp_'//TRIM(ADJUSTL(char_my_rank))//'.txt')
     do i=1, npart_all
        write(60,*) xn_tmp(i), yn_tmp(i), zn_tmp(i), rhon_tmp(i), temperature_tmp(i)
     end do
     IF ( nptmass_all > 0 ) THEN
        DO i=1, nptmass_all
           write(60, *) x_tmp(i), y_tmp(i), z_tmp(i), ptmass_tmp(i)
        END DO 
     END IF
     close(60)
  END IF

! 4. Deallocate and reAllocate arrays xn, yn, zn, rhon, x, y, z, ptmass

! 4.1 Gas particles
  DEALLOCATE ( sphData%xn          )
  DEALLOCATE ( sphData%yn          )
  DEALLOCATE ( sphData%zn          )
  DEALLOCATE ( sphData%rhon        )
  DEALLOCATE ( sphData%temperature )
  DEALLOCATE ( sphData%gasmass     )
  DEALLOCATE ( sphData%hn          ) 

  ALLOCATE   ( sphData%xn(npart_all)          )
  ALLOCATE   ( sphData%yn(npart_all)          )
  ALLOCATE   ( sphData%zn(npart_all)          )
  ALLOCATE   ( sphData%rhon(npart_all)        )
  ALLOCATE   ( sphData%temperature(npart_all) )
  ALLOCATE   ( sphData%gasmass(npart_all)     )
  ALLOCATE   ( sphData%hn(npart_all)          ) 

! 4.2 Point masses
  IF  ( nptmass_all > 0 ) THEN
     DEALLOCATE ( sphData%x   )
     DEALLOCATE ( sphData%y   )
     DEALLOCATE ( sphData%z   )
     DEALLOCATE ( sphData%ptmass )
     ALLOCATE   ( sphData%x(nptmass_all)      )
     ALLOCATE   ( sphData%y(nptmass_all)      )
     ALLOCATE   ( sphData%z(nptmass_all)      )
     ALLOCATE   ( sphData%ptmass(nptmass_all) )
  END IF

! 5. Populate new arrays with values from the tmp storage

! 5.1 Gas particles
  sphData%xn(:)          = xn_tmp(:)
  sphData%yn(:)          = yn_tmp(:)
  sphData%zn(:)          = zn_tmp(:)
  sphData%rhon(:)        = rhon_tmp(:)
  sphData%temperature(:) = temperature_tmp(:)
  sphData%npart          = npart_all
  sphData%gasmass(:)     = gasmass_tmp(:)
  sphData%hn(:)          = hn_tmp(:)

! 5.2 Point masses
  IF  ( nptmass_all > 0 ) THEN
     sphData%x(:)      = x_tmp(:)
     sphData%y(:)      = y_tmp(:)
     sphData%z(:)      = z_tmp(:)
     sphData%ptmass(:) = ptmass_tmp(:)
     sphData%nptmass   = nptmass_all
  END IF

! Write out particle list for tests
  If ( ll_testwrite ) THEN
     open (unit=60, status='replace', file='mpi_test_'//TRIM(ADJUSTL(char_my_rank))//'.txt')
     do i=1, sphData%npart
        write(60,*) sphData%xn(i), sphData%yn(i), sphData%zn(i), sphData%rhon(i), sphData%temperature(i)
     end do
     IF ( nptmass_all > 0 ) THEN
        DO i=1, sphData%nptmass
           write(60, *) sphData%x(i), sphData%y(i), sphData%z(i), sphData%ptmass(i)
        END DO 
     END IF
     close(60)
     open (unit=60, status='replace', file='smoothing_length_'//TRIM(ADJUSTL(char_my_rank))//'.txt')
     do i=1, sphData%npart
        write(60,*) ( sphData%gasmass(i) / sphData%rhon(i) ) ** (1.0/3.0),  sphData%hn(i) 
     end do
     close(60)
  END IF

! 6. Deallocate temporary storage
  DEALLOCATE ( rhon_tmp        )
  DEALLOCATE ( temperature_tmp )
  DEALLOCATE ( zn_tmp          )
  DEALLOCATE ( yn_tmp          )
  DEALLOCATE ( xn_tmp          )
  DEALLOCATE ( hn_tmp          )
  DEALLOCATE ( gasmass_tmp     ) 
  IF  ( nptmass_all > 0 ) THEN
     DEALLOCATE ( x_tmp      )
     DEALLOCATE ( y_tmp      )
     DEALLOCATE ( z_tmp      )
     DEALLOCATE ( ptmass_tmp )
  END IF
  DEALLOCATE ( npart_displ )
  DEALLOCATE ( nptmass_arr )
  DEALLOCATE ( npart_arr   )

! 7. Communicate data to processes which did not run the sph step
  CALL MPI_BCAST(sphData%udist,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  CALL MPI_BCAST(sphData%utime,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  CALL MPI_BCAST(sphData%umass,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  CALL MPI_BCAST(sphData%time,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 

  end subroutine gather_sph_data

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

! Dummy routines  for use in non-mpi case
  subroutine gather_sph_data

    USE sph_data_class
    implicit none

    return
  
  end subroutine gather_sph_data

!-------------------------------------------------------------------------------

  subroutine torus_mpi_barrier(message)

    implicit none

    character(len=*), optional, intent(in) :: message 

    return

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
