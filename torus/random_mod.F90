module random_mod

use kind_mod
use messages_mod
use mpi_global_mod
implicit none

public 
private :: seedFromClockTime, RAN3_double, localSort
contains

  subroutine randomNumberGenerator(reset, putIseed, getIseed, syncIseed, randomSeed, getReal, getDouble, &
       getRealArray,getDoubleArray, getRealArray2d,getDoubleArray2d)
    integer(bigInt),intent(in), optional :: putIseed
    integer(bigInt),intent(out), optional :: getIseed
    logical, optional, intent(in) :: reset
    logical, intent(in), optional :: syncIseed
    logical, intent(in), optional :: randomSeed
    real, intent(out), optional :: getReal
    real, intent(out), optional :: getRealArray(:)
    real, intent(out), optional :: getRealArray2d(:,:)
    real(double), intent(out), optional :: getDouble
    real(double), intent(out), optional :: getDoubleArray(:)
    real(double), intent(out), optional :: getDoubleArray2d(:,:)
    integer(bigInt), save :: iSeed = 0 
    integer(bigInt) :: iseed_syncval
    real :: r
    integer :: i, j

    !$OMP THREADPRIVATE (iSeed)

    if (PRESENT(reset)) then
       if (reset) then
          !$OMP PARALLEL
          r = ran3_double(iseed, reset=.true.)
          !$OMP END PARALLEL
       endif
    endif

    if (PRESENT(putiSeed)) then
       iSeed = putIseed
       r = ran3_double(iseed, reset=.true.)
    endif

    if (PRESENT(getiSeed)) then
       getIseed = iSeed
    endif

    if (PRESENT(syncIseed)) then
       if (syncIseed) then
#ifdef MPI
          call sync_random_seed(iSeed)
#endif

!$OMP PARALLEL DEFAULT (NONE) &
!$OMP SHARED (iseed_syncval) PRIVATE(r)
!$OMP MASTER
          iseed_syncval = iseed
!$OMP END MASTER
!$OMP BARRIER
          iseed = iseed_syncval
          r = ran3_double(iseed, reset=.true.)
!$OMP END PARALLEL 
       endif
    endif

    if (PRESENT(getReal)) then
       getReal = real(ran3_double(iseed))
    endif

    if (PRESENT(randomSeed)) then
!$OMP PARALLEL
       call seedFromClockTime(iSeed)
       r = ran3_double(iseed, reset=.true.)
!$OMP END PARALLEL
    endif

    if (PRESENT(getDouble)) then
       getDouble = ran3_double(iseed)
    endif

    if (PRESENT(getRealArray)) then
       do i = 1, SIZE(getRealArray)
          getRealArray(i) = real(ran3_double(iseed))
       enddo
    endif

    if (PRESENT(getRealArray2d)) then
       do i = 1, SIZE(getRealArray2d,1)
          do j = 1, SIZE(getRealArray2d,2)
             getRealArray2d(i,j) = real(ran3_double(iseed))
          enddo
       enddo
    endif

    if (PRESENT(getDoubleArray2d)) then
       do i = 1, SIZE(getDoubleArray2d,1)
          do j = 1, SIZE(getDoubleArray2d,2)
             getDoubleArray2d(i,j) = ran3_double(iseed)
          enddo
       enddo
    endif


    if (PRESENT(getDoubleArray2d)) then
       do i = 1, SIZE(getDoubleArray)
          getDoubleArray(i) = ran3_double(iseed)
       enddo
    endif

  end subroutine randomNumberGenerator

#ifdef MPI
! Syncronise random seed across all MPI processes. 
  subroutine sync_random_seed(iSeed)

    use mpi
    implicit none

    integer(bigInt) :: iseed
    integer :: ierr

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    call MPI_BCAST(iSeed, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)

  end subroutine sync_random_seed
#endif

! Run tests to check that random seeds are the same/different across
! MPI processes and OpenMP threads. 
  subroutine run_random_test_suite()
    implicit none

    if (torusSerial) return 

    call randomNumberGenerator(randomSeed=.true.)
    call test_random_hybrid()
    call randomNumberGenerator(synciseed=.true.)
    call test_same_hybrid()
    call randomNumberGenerator(reset=.true.)
    call test_same_hybrid()

  end subroutine run_random_test_suite

  subroutine test_random_hybrid()
#ifdef MPI
    use mpi
    use mpi_global_mod, only : nThreadsGlobal
#endif
    use mpi_global_mod, only : myRankGlobal
    integer :: nOmpThreads, iOmpThread, nTot
    integer, allocatable :: itest(:)
    logical :: different
    integer :: iThread, i
    real(double) :: r
#ifdef _OPENMP
    integer :: omp_get_num_threads, omp_get_thread_num
#endif
#ifdef MPI
    integer :: ierr
    integer, allocatable :: itemp(:)
#endif

    nOmpThreads = 1
    iOmpThread = 1

#ifdef _OPENMP
    !$OMP PARALLEL DEFAULT (NONE) &
    !$OMP SHARED(nOmpThreads)
    !$OMP MASTER
    nOmpThreads = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL
#endif

#ifdef MPI
    nTot = nThreadsGlobal * nOmpThreads
#else
    nTot = nOmpThreads
#endif
    allocate(itest(1:nTot))
    itest = 0

    !$OMP PARALLEL DEFAULT (NONE) &
    !$OMP PRIVATE(iOmpThread, ntot, ithread, r) &
    !$OMP SHARED(nOmpThreads, itest, nThreadsGlobal, myrankGlobal)


#ifdef _OPENMP
    iOmpThread = omp_get_thread_num() + 1
#endif
    
 

    iThread = myRankGlobal*nOmpThreads + iOmpThread
    call randomNumberGenerator(getDouble=r)
    itest(iThread) = nint(r * 100000000.d0)
    !$OMP END PARALLEL


#ifdef MPI
    allocate(itemp(SIZE(itest)))
    itemp = 0
    call MPI_REDUCE(itest,itemp,SIZE(itest),MPI_INTEGER,&
         MPI_SUM,0,MPI_COMM_WORLD,ierr)
    itest = itemp
    deallocate(itemp)
#endif
    
    if (myrankGlobal == 0) then
    different = .true.
    call localsort(size(itest),itest)
    different = .true.
    do i = 1, size(itest)-1
       if (itest(i) == itest(i+1)) different=.false.
    enddo
       if (.not.different) then
          write(*,*) "! Threads do not have independent random sequences"
       endif
          write(*,*) "Random sequence test:  Success - threads have independent random sequences"
    endif
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
#endif
    deallocate(itest)
  end subroutine test_random_hybrid
    

  subroutine test_same_hybrid()
#ifdef MPI
    use mpi
    use mpi_global_mod, only : nThreadsGlobal
#endif
    use mpi_global_mod, only : myRankGlobal
    integer :: nOmpThreads, iOmpThread, nTot
    integer, allocatable :: itest(:)
    logical :: different
    integer :: iThread, i
    real(double) :: r
#ifdef _OPENMP
    integer :: omp_get_num_threads, omp_get_thread_num
#endif
#ifdef MPI
    integer :: ierr
    integer, allocatable :: itemp(:)
#endif

    nOmpThreads = 1
    iOmpThread = 1

#ifdef _OPENMP
    !$OMP PARALLEL DEFAULT (NONE) &
    !$OMP SHARED (nOmpThreads)
    !$OMP MASTER
    nOmpThreads = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP END PARALLEL
#endif

#ifdef MPI
    nTot = nThreadsGlobal * nOmpThreads
#else
    nTot = nOmpThreads
#endif
    allocate(itest(1:nTot))
    itest = 0

    !$OMP PARALLEL DEFAULT (NONE) &
    !$OMP PRIVATE(iOmpThread, ntot, ithread, r) &
    !$OMP SHARED(nOmpThreads, itest, nThreadsGlobal, myrankGlobal)


#ifdef _OPENMP
    iOmpThread = omp_get_thread_num() + 1
#endif
    
 

    iThread = myRankGlobal*nOmpThreads + iOmpThread
    call randomNumberGenerator(getDouble=r)
    itest(iThread) = nint(r * 100000000.d0)
    !$OMP END PARALLEL


#ifdef MPI
    allocate(itemp(SIZE(itest)))
    itemp = 0
    call MPI_REDUCE(itest,itemp,SIZE(itest),MPI_INTEGER,&
         MPI_SUM,0,MPI_COMM_WORLD,ierr)
    itest = itemp
    deallocate(itemp)
#endif
    
    if (myrankGlobal == 0) then
    different = .false.
    call localsort(size(itest),itest)
    different = .false.
    do i = 1, size(itest)-1
       if (itest(i) /= itest(i+1)) different=.true.
    enddo
       if (different) then
          write(*,*) "! Threads do not have same random sequences"
          write(*,*) itest
       else
          write(*,*) "Synced seed test: Success - threads have same random number sequences"
       endif
    endif
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
#endif
    deallocate(itest)
  end subroutine test_same_hybrid
    
#ifdef MPI
! Test whether all threads are producing independent random numbers
  subroutine test_random_across_threads(debug)
    use mpi
    use mpi_global_mod, only : nThreadsGlobal, myRankGlobal
    implicit none
    real(double) :: r
    integer :: i, ierr
    integer, allocatable :: itest(:), itemp(:)
    logical :: different
    logical, optional :: debug
    call randomNumberGenerator(getDouble=r)
    i = nint(r * 100000000.d0)
    allocate(itest(1:nThreadsGlobal))
    itest = 0
    itest(myRankGlobal+1) = i
    allocate(itemp(SIZE(itest)))
     itemp = 0
     call MPI_REDUCE(itest,itemp,SIZE(itest),MPI_INTEGER,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     itest = itemp
     deallocate(itemp)
     if (myrankGlobal ==0) then
        if (present(debug)) then
           if (debug) then
              do i = 1, nThreadsGlobal
                 write(*,*) "Thread ",i-1," random integer ", itest(i)
              enddo
           endif
        endif
        different = .true.
        call localsort(size(itest),itest)
        different = .true.
        do i = 1, size(itest)-1
           if (itest(i) == itest(i+1)) different=.false.
        enddo
        if (.not.different) then
           call writeFatal("Thread do not have independent random sequences")
           stop
        endif
     endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     deallocate(itest)
   end subroutine test_random_across_threads
#endif

#ifdef _OPENMP
  subroutine testRandomOMP()
    integer :: omp_get_num_threads, omp_get_thread_num
    integer :: np
    integer :: nTest = 100
    real :: r
    integer :: i, j 
    integer, allocatable :: array(:,:)
    logical :: allTheSame

    allocate(array(1:nTest, 8))
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(r, i,j) &
    !$OMP SHARED(array, np, ntest)
    np =  omp_get_num_threads()
    j = omp_get_thread_num()+1
    do i = 1, nTest
       call randomNumberGenerator(getReal=r)
       array(i,j) = int(r*10000000.d0)
    enddo
    !$OMP END PARALLEL
    !$OMP BARRIER
    allTheSame = .true.
    do i = 1, ntest
       do j = 2, np
          if (array(i,1) /= array(i,j)) then
             allTheSame = .false.
          endif
       enddo
    enddo
    deallocate(array)
    if (allTheSame) then
       call writeFatal("OMP random sequences are the same")
    endif
  end subroutine testRandomOMP
#endif

  SUBROUTINE seedFromClockTime(ibig)
    use mpi_global_mod, only: myRankGlobal

    integer(bigint) :: ibig
    integer :: iValues(8)
    integer :: j, nt
#ifdef _OPENMP
     integer ::  omp_get_thread_num
     integer ::  omp_get_num_threads
#endif

    j = 0
    nt = 1
    CALL DATE_AND_TIME(values=iValues)
    ibig = ivalues(5) * 60 * 60 * 1000 +&
         ivalues(6) * 60 * 1000 + &
         ivalues(7) * 1000 +&
         ivalues(8)
#ifdef _OPENMP
     j = omp_get_thread_num()+1
     nt = omp_get_num_threads()
#endif
     
    ibig = (j+(myRankGlobal+1)*nt)*ibig
  END SUBROUTINE seedFromClockTime


  real(double) FUNCTION RAN3_double(IDUM, reset) result (ran3)
    implicit none
    real(double), parameter ::  MBIG=1000000000000000.d0,MSEED=161803398331234.d0,MZ=0.d0
    real(double), parameter :: fac = 1.d-15
    real(double),save  :: ma(55)
    integer(bigint),save :: iff = 0
    integer(bigint) :: idum, ii, i, k
    real(double) :: mj, mk
    integer(bigint), save :: inext, inextp
    logical, optional :: reset
    logical :: doreset
!$OMP THREADPRIVATE (ma, iff, inext, inextp)

    doreset = .false.
      if (present(reset)) doreset=reset 
     IF(doreset.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-ABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO I=1,54
          II=MOD(21*I,55_bigint)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
          enddo
        DO  K=1,4
          DO  I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55_bigint))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
         enddo
         enddo
        INEXT=0
        INEXTP=31
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      if ((ran3 < tiny(ran3)).or.(ran3 > 0.99999999999999d0)) then
         write(*,*) "ran3 bug ",ran3
         ran3 = min(0.99999999999999d0,max(ran3, tiny(ran3)))
      endif
    END FUNCTION RAN3_double

  SUBROUTINE localsort(N,RA)
    INTEGER N, L, IR, I, J
    integer RRA
    integer RA(N)
    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
    ELSE
       RRA=RA(IR)
       RA(IR)=RA(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).LT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.LT.RA(J))THEN
          RA(I)=RA(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    GO TO 10
  END subroutine localsort

end module random_mod
