module mpi_amr_mod
#ifdef MPI

  use kind_mod
  use amr_mod
  use mpi_global_mod
  implicit none

  real(double), allocatable :: buffer(:)
  integer :: nBuffer, maxBuffer

  interface packAttributeStatic
     module procedure packAttributeIntegerSingle
     module procedure packAttributeDoubleSingle
     module procedure packAttributeLogicalSingle
     module procedure packAttributeIntegerArray
     module procedure packAttributeLogicalArray
     module procedure packAttributeDoubleArray
     module procedure packAttributeRealArray
     module procedure packAttributeVectorSingle
  end interface


  interface packAttributePointer
     module procedure packAttributeDoubleArray1Dpointer
     module procedure packAttributeRealArray1Dpointer
     module procedure packAttributeVectorArray1Dpointer
     module procedure packAttributeIntegerArray1Dpointer
     module procedure packAttributeLogicalArray1Dpointer
     module procedure packAttributeDoubleArray2Dpointer
     module procedure packAttributeDoubleArray3Dpointer
  end interface

  interface unpackAttributePointer
     module procedure unpackAttributeDoubleArray1Dpointer
     module procedure unpackAttributeRealArray1Dpointer
     module procedure unpackAttributeVectorArray1Dpointer
     module procedure unpackAttributeLogicalArray1Dpointer
     module procedure unpackAttributeIntegerArray1Dpointer
     module procedure unpackAttributeDoubleArray2Dpointer
     module procedure unpackAttributeDoubleArray3Dpointer
  end interface

  interface unpackAttributeStatic
     module procedure unpackAttributeIntegerSingle
     module procedure unpackAttributeDoubleSingle
     module procedure unpackAttributeLogicalSingle
     module procedure unpackAttributeIntegerArray
     module procedure unpackAttributeLogicalArray
     module procedure unpackAttributeDoubleArray
     module procedure unpackAttributeRealArray
     module procedure unpackAttributeVectorSingle
  end interface


  private buffer, nbuffer, maxBuffer

contains

  subroutine setupAMRCOMMUNICATOR
    use mpi
    use inputs_mod, only : nHydroThreadsinput, splitOverMPI, amr1d, amr2d, amr3d, loadBalancing
    integer :: ierr, i, j, iset, k
    integer, allocatable :: ranks(:)
    integer :: worldGroup, amrGroup, localWorldGroup
    integer :: amrParallelGroup, zeroThreadGroup
    integer, parameter :: tag = 22
    character(len=40) :: proc, thisProc
    logical :: optimized
    integer :: status(MPI_STATUS_SIZE)
    character(len=80) :: message


    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreadsGlobal, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRankWorldGlobal, ierr)
    nHydroSetsGlobal = 1
    myHydroSetGlobal = 0
    nHydroThreadsGlobal = nHydroThreadsinput
    if (splitOverMPI) then

       if (loadBalancing) then
          if (nHydroThreadsInput == 0) then
             call torus_abort("nHydrothreads must be specified for load balancing")
          else
             nHydroThreadsGlobal  = nHydroThreadsInput
             nLoadBalancingThreadsGlobal = nThreadsGlobal - nHydroThreadsGlobal - 1
             nHydroSetsGlobal = 1
             myHydroSetGlobal = 0
             loadBalancingThreadGlobal = .false.
             if (myRankWorldGlobal > nHydroThreadsGlobal) loadBalancingThreadGlobal = .true.
          endif
       else

          if (nHydroThreadsInput == 0) then

             if (amr3d.and.(mod(nThreadsGlobal, 513) == 0)) then
                nHydroThreadsGlobal = 512
             else if (amr3d.and.(mod(nThreadsGlobal, 65) == 0)) then
                nHydroThreadsGlobal = 64
             else if (amr3d.and.(mod(nThreadsGlobal, 9) == 0)) then
                nHydroThreadsGlobal = 8
             else if (amr2d.and.(mod(nThreadsGlobal, 65) == 0)) then
                nHydroThreadsGlobal = 64
             else if (amr2d.and.(mod(nThreadsGlobal, 17) == 0)) then
                nHydroThreadsGlobal = 16
             else if (amr2d.and.(mod(nThreadsGlobal, 5) == 0)) then
                nHydroThreadsGlobal = 4
             else if (amr1d.and.(mod(nThreadsGlobal, 3) == 0)) then
                nHydroThreadsGlobal = 2
             else
                write(*,*) "Number of MPI threads is ",nThreadsGlobal
                write(*,*) "Can't figure out automatically what domain decomp is required"
                stop
             endif
          endif

          if (mod(nThreadsGlobal, (nHydroThreadsGlobal+1)) /= 0) then
             write(*,*) "Number of MPI threads is ",nThreadsGlobal
             write(*,*) "Number of threads per decomposed domain is ", nHydroThreadsGlobal+1
             write(*,*) "Can't distribute work properly"
             stop
          endif

          nHydroSetsGlobal = nThreadsGlobal / (nHydroThreadsGlobal+1)


          myHydroSetGlobal = myRankWorldGlobal / (nHydroThreadsGlobal+1)
       endif
       if (myrankWorldGlobal == 0) then
          write(*,*) " "
          write(*,*) "Parallel info:"
          write(*,*) "nThreadsGlobal: ",nThreadsGlobal
          write(*,*) "nHydrothreads: ",nHydroThreadsGlobal
          write(*,*) "nHydroSetsGlobal: ",nHydroSetsGlobal
          write(*,*) "nLoadBalancingThreadsGlobal: ",nLoadBalancingThreadsGlobal
          write(*,*) " "
       endif

       call MPI_COMM_GROUP(MPI_COMM_WORLD, worldGroup, ierr)

       allocate(ranks(1:(nHydroThreadsGlobal+1)))
       do i = 1, nHydroThreadsGlobal+1
          ranks(i) = myHydroSetGlobal * (nHydroThreadsGlobal) + i - 1
       enddo
       call MPI_GROUP_INCL(worldGroup, nHydroThreadsGlobal+1, ranks, localWorldGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, localWorldGroup, zeroplusAmrCOMMUNICATOR, ierr)
       deallocate(ranks)

       allocate(ranks(1:(nHydroThreadsGlobal+1+nLoadBalancingThreadsGlobal)))
       do i = 1, nHydroThreadsGlobal+1+nLoadBalancingThreadsGlobal
          ranks(i) = myHydroSetGlobal * (nHydroThreadsGlobal+1+nLoadBalancingThreadsGlobal) + i - 1
       enddo
       call MPI_GROUP_INCL(worldGroup, nHydroThreadsGlobal+1+nLoadBalancingThreadsGlobal, ranks, localWorldGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, localWorldGroup, localWorldCOMMUNICATOR, ierr)
       deallocate(ranks)




       allocate(ranks(1:1+nLoadBalancingThreadsGlobal))
       ranks(1) = 0
       do i = 1, nLoadBalancingThreadsGlobal
          ranks(1+i) = nHydroThreadsGlobal+i
       enddo
       call MPI_GROUP_EXCL(localWorldGroup, 1+nLoadBalancingThreadsGlobal, ranks, amrGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, amrGroup, amrCOMMUNICATOR, ierr)
       deallocate(ranks)

       allocate(ranks(1:1))
       ranks(1) = 0
       call MPI_GROUP_EXCL(localWorldGroup, 1, ranks, amrGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, amrGroup, allDomainsCOMMUNICATOR, ierr)
       deallocate(ranks)



       if (loadBalancingThreadGlobal) then
          call MPI_COMM_RANK(MPI_COMM_WORLD, myRankGlobal, ierr)
       else
          call MPI_COMM_RANK(localWorldCommunicator, myRankGlobal, ierr)
          copyOfThread = myrankGlobal
       endif
       allocate(amrParallelCommunicator(1:nHydroThreadsGlobal))
       allocate(ranks(1:nHydroSetsGlobal))
       do i = 1, nHydroThreadsGlobal
          do j = 1, nHydroSetsGlobal
             ranks(j) = (j-1) * (nHydroThreadsGlobal+1) + i
          enddo
          call MPI_GROUP_INCL(worldGroup, nHydroSetsGlobal, ranks, amrParallelGroup, ierr)
          call MPI_COMM_CREATE(MPI_COMM_WORLD, amrParallelGroup, amrParallelCommunicator(i), ierr)
       enddo
       deallocate(ranks)

       allocate(ranks(1:nHydroSetsGlobal))
       do j = 1, nHydroSetsGlobal
          ranks(j) = (j-1) * (nHydroThreadsGlobal+1) 
       enddo
       call MPI_GROUP_INCL(worldGroup, nHydroSetsGlobal, ranks, zeroThreadGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, zeroThreadGroup, zeroThreadCommunicator, ierr)
       deallocate(ranks)


       call mpi_get_processor_name(proc, j, ierr)
       do i = 0, nThreadsGlobal-1
          call mpi_barrier(MPI_COMM_WORLD, ierr)
          if (myRankWorldGlobal == i) then
             write(*,*) "Thread ",i, " is running on ",trim(proc)
          endif
          call mpi_barrier(MPI_COMM_WORLD, ierr)
       enddo

       if ((nHydroThreadsGlobal == 8).and.(.not.loadBalancingThreadGlobal)) then
          optimized = .true.
          do iSet = 1, nHydroSetsGlobal

             if (myHydroSetGlobal == (iSet-1)) then
                if (myRankGlobal == 0) then
                   call mpi_get_processor_name(proc, j, ierr)
                   do i = 1, nHydroThreadsGlobal
                      call mpi_recv(k, 1, MPI_INTEGER, i, tag, localWorldcommunicator, status, ierr)
                      call mpi_recv(thisProc, k, MPI_CHARACTER, i, tag, localWorldcommunicator, status, ierr)
                      if (thisProc(1:k) /= proc(1:j)) then
                         optimized = .false.
                      endif
                   enddo
                else
                   call mpi_get_processor_name(proc, j, ierr)
                   call mpi_send(j, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, ierr)
                   call mpi_send(proc, j, MPI_CHARACTER, 0, tag, localWorldCommunicator, ierr)
                endif
             endif
          enddo


          if (.not.optimized) then
             write(message,'(a)') "!!! Eight-way domain decomposition not optimized."
             if (myrankGlobal ==0) write(*,*) trim(message)
          endif
       endif


       !       if ((nHydroThreadsGlobal == 8).and.(nHydroSetsGlobal > nHydroThreadsGlobal)) then
       !          allocate(ranks(1:nhydroThreadsGlobal**2))
       !          n = 1
       !          ranks(1) = 0
       !          do i = 1, nHydroThreadsGlobal
       !             do j = 1, nHydroThreadsGlobal
       !                n = n + 1
       !                ranks(n) = (i-1) * (nHydroThreadsGlobal+1) + j
       !             enddo
       !          enddo
       !          call MPI_GROUP_INCL(worldGroup, nHydroThreadsGlobal**2+1, ranks, hydroGroup, ierr)
       !          call MPI_COMM_CREATE(MPI_COMM_WORLD, hydroGroup, hydroCommunicator, ierr)
       !          deallocate(ranks)
       !       endif


       !       call testMPIspeeds()
    else
       myrankGlobal = myrankWorldGlobal
    endif
  end subroutine setupAMRCOMMUNICATOR

  subroutine testMPIspeeds()
    use timing, only : tune
    use mpi
    integer, parameter :: ndata = 100000
    real(double) :: rdata(ndata)
    real(double) :: sdata(ndata)
    integer :: i, j, iThread, iset
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 54

    do i = 1, nData
       rData(i) = dble(i)
    enddo
    do iset = 1, nHydroSetsGlobal

       call mpi_barrier(MPI_COMM_WORLD, ierr)
       if (myrankWorldGlobal == 0) call tune(6, "Local send/receive test")  ! start a stopwatch
       if (myHydroSetGlobal == (iset-1)) then
          call mpi_barrier(localWorldCommunicator, ierr)
          do j = 1, 1000
             if (myrankGlobal == 0) then
                do iThread = 1, nHydroThreadsGlobal
                   call MPI_SEND(rData, nData, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
                enddo
             else
                call MPI_RECV(sData, nData, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
             endif
             call mpi_barrier(localWorldCommunicator, ierr)
          enddo
       endif
       call mpi_barrier(MPI_COMM_WORLD, ierr)
       if (myrankWorldGlobal == 0) call tune(6, "Local send/receive test")  ! stop a stopwatch
    enddo
  end subroutine testMPIspeeds

! Free communicator created by setupAMRCOMMUNICATOR
  subroutine freeAMRCOMMUNICATOR
    use mpi
    implicit none

    integer :: ierr

    if ( myRankGlobal /= 0 ) then 
       if (amrCommunicator /= MPI_COMM_NULL) then
          call MPI_COMM_FREE(amrCOMMUNICATOR, ierr)
       endif
    end if

  end subroutine freeAMRCOMMUNICATOR

  subroutine findMassOverAllThreads(grid, mass)
    use mpi
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: mass
    real(double), allocatable :: massOnThreads(:), temp(:), volumeOnThreads(:)
    integer :: ierr
!    real(double) :: totalVolume

    if (loadBalancingThreadGlobal) goto 666

    allocate(massOnThreads(1:nThreadsGlobal), temp(1:nThreadsGlobal))
    allocate(volumeOnThreads(1:nThreadsGlobal))
    massOnThreads = 0.d0
!    volumeOnThreads = 0.d0
    temp = 0.d0
    if (.not.grid%splitOverMpi) then
       call writeWarning("findMassOverAllThreads: grid not split over MPI")
       mass = 0.d0
       goto 666
    endif

    if (myRankGlobal /= 0) then
       call findtotalMassMPI(grid%octreeRoot, massOnThreads(myRankGlobal))
       call MPI_ALLREDUCE(massOnThreads, temp, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
       mass = SUM(temp(1:nThreadsGlobal))
       !call MPI_ALLREDUCE(volumeOnThreads, temp, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
       !totalVolume = sum(temp(1:nThreadsGlobal))
       !print *, "totalVolume = ", totalVolume
    end if
666 continue
  end subroutine findMassOverAllThreads

    
  recursive subroutine findTotalMassMPI(thisOctal, totalMass, minRho, maxRho)
    use inputs_mod, only : hydrodynamics, cylindricalHydro, spherical
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: rVec
  real(double) :: totalMass
  real(double),optional :: minRho, maxRho
  real(double) :: dv!, totalVolume
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findtotalMassMPI(child, totalMass, minRho, maxRho)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                dv = cellVolume(thisOctal, subcell)*1.d30
!                print *, "dv", dv, thisOctal%subcellSize**3

                if (hydrodynamics) then
                   if (thisOctal%twoD) then
                      if (cylindricalHydro) then
                         dv = cellVolume(thisOctal, subcell) * 1.d30
                         rVec = subcellCentre(thisOctal,subcell)
                         if (rVec%x < 0.d0) dv = 0.d0
                         if (thisOctal%ghostCell(subcell)) dv = 0.d0
                      else
                         dv = thisOctal%subcellSize**2
                      endif
                   else if (thisOctal%oned) then
                      if (spherical) then
                         dv = cellVolume(thisOctal, subcell) * 1.d30
                         rVec = subcellCentre(thisOctal,subcell)
                         if (rVec%x < 0.d0) dv = 0.d0
                         if (thisOctal%ghostCell(subcell)) dv = 0.d0
                      else
                         dv = thisOctal%subcellSize**2
                      endif
                   endif
                endif
                !             totalVolume = totalVolume + dv
                totalMass = totalMass + thisOctal%rho(subcell) * dv
                if (PRESENT(minRho)) then
                   if (thisOctal%rho(subcell) > 1.d-20) then
                      minRho = min(thisOctal%rho(subcell), minRho)
                   endif
                endif
                if (PRESENT(maxRho)) maxRho = max(dble(thisOctal%rho(subcell)), maxRho)
             endif
          endif
       end if
    enddo
  end subroutine findTotalMassMPI

  subroutine findkeOverAllThreads(grid, ke)
    use mpi
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: ke
    real(double), allocatable :: keOnThreads(:), temp(:)
    integer :: ierr
!    real(double) :: totalVolume

    allocate(keOnThreads(1:nThreadsGlobal), temp(1:nThreadsGlobal))
    keOnThreads = 0.d0
    temp = 0.d0
    if (.not.grid%splitOverMpi) then
       call writeWarning("findKEOverAllThreads: grid not split over MPI")
       ke = 0.d0
       goto 666
    endif

    if (myRankGlobal /= 0) then
       call findtotalKEMPI(grid%octreeRoot, keOnThreads(myRankGlobal))
       call MPI_ALLREDUCE(keOnThreads, temp, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
       ke = SUM(temp(1:nThreadsGlobal))
    end if
666 continue
  end subroutine findKEOverAllThreads

    
  recursive subroutine findTotalKEMPI(thisOctal, totalKE)
    use inputs_mod, only : hydrodynamics, cylindricalHydro, spherical
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: rVec
  real(double) :: totalKE
  real(double) :: dv!, totalVolume
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findtotalKEMPI(child, totalKE)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                dv = cellVolume(thisOctal, subcell)*1.d30
!                print *, "dv", dv, thisOctal%subcellSize**3

                if (hydrodynamics) then
                   if (thisOctal%twoD) then
                      if (cylindricalHydro) then
                         dv = cellVolume(thisOctal, subcell) * 1.d30
                         rVec = subcellCentre(thisOctal,subcell)
                         if (rVec%x < 0.d0) dv = 0.d0
                         if (thisOctal%ghostCell(subcell)) dv = 0.d0
                      else
                         dv = thisOctal%subcellSize**2
                      endif
                   else if (thisOctal%oned) then
                      if (spherical) then
                         dv = cellVolume(thisOctal, subcell) * 1.d30
                         rVec = subcellCentre(thisOctal,subcell)
                         if (rVec%x < 0.d0) dv = 0.d0
                         if (thisOctal%ghostCell(subcell)) dv = 0.d0
                      else
                         dv = thisOctal%subcellSize**2
                      endif
                   endif
                endif
                totalKE = totalKE  + 0.5d0 * (thisOctal%rho(subcell) * dv) * &
                     (modulus(thisOctal%velocity(subcell))*cspeed)**2
             endif
          endif
       end if
    enddo
  end subroutine findTotalKEMPI


  subroutine findMassOverAllThreadsWithinR(grid, mass, radius)
    use mpi
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: mass
    real(double) :: radius, localMass
    integer, parameter :: tag = 54
    integer :: ierr, iThread
    integer :: status(MPI_STATUS_SIZE)

    mass = 0.d0

    do iThread = 1, nHydroThreadsGlobal
       if (myRankGlobal == iThread) then
          localMass = 0.d0
          call findTotalMassWithinRMPIPrivate(grid%octreeRoot, radius, localMass)
       else
          call MPI_SEND(radius, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
!          call MPI_RECV(localMass, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
!       endif
!       mass = mass + localMass
    enddo

    mass = mass + localMass

    do iThread = 1, nHydroThreadsGlobal
       if (myRankGlobal /= iThread) then
          call MPI_RECV(localMass, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
          mass = mass + localMass
       endif
    enddo

  end subroutine findMassOverAllThreadsWithinR
    

  subroutine findTotalMassWithinRServer(grid, receiveThread)
    use mpi
    type(GRIDTYPE) :: grid
    real(double) :: radius, totalMass
    integer :: ierr
    integer :: receiveThread
    integer, parameter :: tag = 54
    logical :: stillServing
    integer :: status(MPI_STATUS_SIZE)

    stillServing = .true.

    call findTotalMassWithinRMPI(grid, radius, totalMass, reset=.true.)

    do while(stillServing)

       call MPI_RECV(radius, 1, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, status, ierr)
       if (radius > 1.d29) then
          stillServing = .false.
       else
          totalmass = 0.d0
          call findTotalMassWithinRMPI(grid, radius, totalMass)
          call MPI_SEND(totalMass, 1, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine findTotalMassWithinRServer

  subroutine findTotalPhiOutsideRServer(grid, receiveThread)
    use mpi
    type(GRIDTYPE) :: grid
    real(double) :: radius, totalPhi
    integer :: ierr
    integer :: receiveThread
    integer, parameter :: tag = 55
    logical :: stillServing
    integer :: status(MPI_STATUS_SIZE)

    stillServing = .true.

    call findTotalPhiOutsideRMPI(grid, radius, totalPhi, reset=.true.)

    do while(stillServing)

       call MPI_RECV(radius, 1, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, status, ierr)
       if (radius > 1.d29) then
          stillServing = .false.
       else
          totalphi = 0.d0
          call findTotalPhiOutsideRMPI(grid, radius, totalPhi)
          call MPI_SEND(totalPhi, 1, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine findTotalPhiOutsideRServer


  subroutine findTotalMassWithinRMPI(grid, radius, totalMass, reset)
    type(GRIDTYPE) :: grid
    real(double) :: radius
    real(double) :: totalMass
    logical, optional :: reset
    logical :: doReset
    logical, save :: firstTime = .true.
    real(double), save :: savedRadius, savedMass

    doReset = .false.
    if (PRESENT(reset)) then
       doReset = reset
    endif
    if (doReset) firstTime = .true.
    
    if (firstTime) then
       savedRadius = 0.d0
       call findMaxR(grid%octreeRoot, savedRadius)
       savedMass = 0.d0
       call findTotalMassWithinRMPIPrivate(grid%octreeRoot, savedRadius, savedMass)
       firstTime = .false.
    endif

    if (radius > savedRadius) then
       totalMass = savedMass
    else
       totalMass = 0.d0
       call findTotalMassWithinRMPIPrivate(grid%octreeRoot, radius, totalMass)
    endif
  end subroutine findTotalMassWithinRMPI

  subroutine findTotalPhiOutsideRMPI(grid, radius, totalPhi, reset)
    type(GRIDTYPE) :: grid
    real(double) :: radius
    real(double) :: totalPhi
    logical, optional :: reset
    logical :: doReset
    logical, save :: firstTime = .true.
    real(double), save :: savedRadius, savedPhi

    doReset = .false.
    if (PRESENT(reset)) then
       doReset = reset
    endif
    if (doReset) firstTime = .true.
    
    if (firstTime) then
       savedRadius = 1.d30
       call findMinR(grid%octreeRoot, savedRadius)
       savedPhi = 0.d0
       call findTotalPhiOutsideRMPIPrivate(grid%octreeRoot, savedRadius, savedPhi)
       firstTime = .false.
    endif

    if (radius < savedRadius) then
       totalPhi = savedPhi
    else
       totalPhi = 0.d0
       call findTotalPhiOutsideRMPIPrivate(grid%octreeRoot, radius, totalPhi)
    endif
  end subroutine findTotalPhiOutsideRMPI


  recursive subroutine findTotalMassWithinRMPIPrivate(thisOctal, radius, totalMass)
    use inputs_mod, only : hydrodynamics, cylindricalHydro
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: rVec
  real(double) :: totalMass, radius
  real(double) :: dv
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findtotalMassWithinRMPIPrivate(child, radius,totalMass)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                rVec = subcellCentre(thisOctal, subcell)
                if (modulus(rVec) < radius) then
                   dv = cellVolume(thisOctal, subcell)*1.d30
                   
                   if (hydrodynamics) then
                      if (thisOctal%twoD) then
                         if (cylindricalHydro) then
                            dv = cellVolume(thisOctal, subcell) * 1.d30
                            if (thisOctal%ghostCell(subcell)) dv = 0.d0
                         else
                            dv = thisOctal%subcellSize**2
                         endif
                      else if (thisOctal%oned) then
                         dv = cellVolume(thisOctal, subcell) * 1.d30
                         if (thisOctal%ghostCell(subcell)) dv = 0.d0
                      endif
                   endif
                   totalMass = totalMass + thisOctal%rho(subcell) * dv
                endif
             endif
          endif
       end if
    enddo
  end subroutine findTotalMassWithinRMPIPrivate

  recursive subroutine findTotalPhiOutsideRMPIPrivate(thisOctal, radius, totalPhi)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: rVec
  real(double) :: totalPhi, radius
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findtotalPhiOutsideRMPIPrivate(child, radius,totalPhi)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                rVec = subcellCentre(thisOctal, subcell)
                if (modulus(rVec) > radius) then
                   totalPhi = totalPhi - thisOctal%phi_gas(subcell) 
                endif
             endif
          endif
       end if
    enddo
  end subroutine findTotalPhiOutsideRMPIPrivate

  recursive subroutine findMaxR(thisOctal, radius)
!   use inputs_mod, only : hydrodynamics, cylindricalHydro
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: rVec
  real(double) :: radius
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findMaxR(child, radius)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if(.not. thisoctal%ghostcell(subcell)) then
             rVec = subcellCentre(thisOctal, subcell)
             if (modulus(rVec) > radius) then
                radius = modulus(rVec)
             endif
          endif
       end if
    enddo
  end subroutine findMaxR

  recursive subroutine findMinR(thisOctal, radius)
!   use inputs_mod, only : hydrodynamics, cylindricalHydro
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: rVec
  real(double) :: radius
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findMinR(child, radius)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if(.not. thisoctal%ghostcell(subcell)) then
             rVec = subcellCentre(thisOctal, subcell)
             if (modulus(rVec) < radius) then
                radius = modulus(rVec)
             endif
          endif
       end if
    enddo
  end subroutine findMinR

  subroutine findAngMomOverAllThreads(grid, angMom, centre)
    use mpi
    type(GRIDTYPE) :: grid
    type(VECTOR), intent(out) :: angMom
    type(VECTOR) :: centre
    real(double), allocatable :: temp(:), temp2(:)
    type(VECTOR), allocatable :: angMomOnThreads(:)
    integer :: ierr

    allocate(angMomOnThreads(1:nThreadsGlobal), temp(1:nThreadsGlobal), temp2(1:nThreadsGlobal))
    angMomOnThreads = VECTOR(0.d0,0.d0,0.d0)
    temp = 0.d0
    if (.not.grid%splitOverMpi) then
       call writeWarning("findAngMomOverAllThreads: grid not split over MPI")
       angMom = VECTOR(0.d0,0.d0,0.d0)
       goto 666
    endif

    if (myRankGlobal /= 0) then
       call findAngMomMPI(grid%octreeRoot, angMomOnThreads(myRankGlobal), centre)
       temp2 = angMomOnThreads%x
       call MPI_ALLREDUCE(temp2, temp, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
       angMom%x = SUM(temp(1:nThreadsGlobal))
       temp2 = angMomOnThreads%y
       call MPI_ALLREDUCE(temp2, temp, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
       angMom %y= SUM(temp(1:nThreadsGlobal))
       temp2 = angMomOnThreads%z
       call MPI_ALLREDUCE(temp2, temp, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
       angMom%z = SUM(temp(1:nThreadsGlobal))
    end if
    deallocate(angMomOnThreads, temp, temp2)
666 continue
  end subroutine findAngMomOverAllThreads
    
  recursive subroutine findAngMomMPI(thisOctal, angMom, centre)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: centre, cellCentre, rVec, angMom, vel
  real(double) :: dv, dm
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findangMomMPI(child, angMom, centre)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                if (thisOctal%threed) then
                   cellCentre = subcellCentre(thisOctal, subcell)
                   vel = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                        thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
                        thisOctal%rhow(subcell)/thisOctal%rho(subcell))
                   dv = cellVolume(thisOctal, subcell)*1.d30
                   dm = thisOctal%rho(subcell) * dv
                   rVec = 1.d10*(centre-cellCentre)
                   angMom = angMom + (rVec.cross.(dm*vel))
                else
                   angMom = angMom + VECTOR(0.d0, 0.d0, thisOctal%rhov(subcell) * cellVolume(thisOctal, subcell)*1.d30)
                endif
             endif
          endif
       end if
    enddo
  end subroutine findAngMomMPI

  subroutine findEnergyOverAllThreads(grid, energy)
    use mpi
!    use inputs_mod, only : cylindricalHydro
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: energy
    real(double), allocatable :: energyOnThreads(:), temp(:)
    integer :: ierr

    allocate(energyOnThreads(1:nThreadsGlobal), temp(1:nThreadsGlobal))
    energyOnThreads = 0.d0
    temp = 0.d0
    if (.not.grid%splitOverMpi) then
       call writeWarning("findEnergyOverAllThreads: grid not split over MPI")
       energy = 0.d0
       goto 666
    endif

    if (myRankGlobal /= 0) then
       call findtotalEnergyMPI(grid%octreeRoot, energyOnThreads(myRankGlobal))
       call MPI_ALLREDUCE(energyOnThreads, temp, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM,amrCOMMUNICATOR, ierr)
       energy = SUM(temp(1:nThreadsGlobal))
    end if
666 continue
  end subroutine findEnergyOverAllThreads
    
  recursive subroutine findTotalEnergyMPI(thisOctal, totalEnergy)
    use inputs_mod, only : cylindricalHydro
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: totalEnergy
  real(double) :: dv
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findtotalEnergyMPI(child, totalEnergy)
                exit
             end if
          end do
       else
          if(.not. thisoctal%ghostcell(subcell)) then
             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                if (thisOctal%threed) then
!                   dv = cellVolume(thisOctal, subcell) * 1.d30
                   dv = thisOctal%subcellSize**3
                else if (thisOctal%twoD) then
                   if (cylindricalhydro) then
                      dv = cellVolume(thisOctal, subcell) * 1.d30
                   else
                      dv = thisOctal%subcellSize**2
                   endif
                else if (thisOctal%oneD) then
                   dv = thisOctal%subcellSize
                endif
                totalEnergy = totalEnergy + thisOctal%rhoe(subcell) * dv
             endif
          endif
       end if
    enddo
  end subroutine findTotalEnergyMPI


! Causes warnings about precision loss with gfortran 4.6.x so these routines are commented 
! out as they are not used. 
! D. Acreman November 2011  

!!$  subroutine getSquares(grid, plane, valueName, nSquares, corners, value, speed, ang)
!!$    use mpi
!!$    type(GRIDTYPE) :: grid
!!$    character(len=*) :: plane, valueName
!!$    real, pointer :: corners(:,:), value(:), speed(:), ang(:)
!!$    integer :: myRank, ierr, nThreads, iThread
!!$    integer :: nSquares, n, tag=97,i
!!$    integer :: status(MPI_STATUS_SIZE)
!!$    integer :: j
!!$
!!$    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!!$    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
!!$    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)
!!$
!!$    if (myRank == 1) then
!!$       nSquares = 0
!!$       call recursGetSquares(grid%octreeRoot, grid, plane, valueName, nSquares, corners, value, speed, ang)
!!$       do iThread = 2, nThreads - 1
!!$          call MPI_RECV(n, 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
!!$          do i = 1, 4
!!$             call MPI_RECV(corners(nSquares+1:nSquares+n, i), n, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
!!$          enddo
!!$          call MPI_RECV(value(nSquares+1:nSquares+n), n, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
!!$          call MPI_RECV(speed(nSquares+1:nSquares+n), n, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
!!$          call MPI_RECV(ang(nSquares+1:nSquares+n), n, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
!!$          nSquares = nSquares+n
!!$       end do
!!$    else
!!$       nSquares = 0
!!$       call recursGetSquares(grid%octreeRoot, grid, plane, valueName, nSquares, corners, value, speed, ang)
!!$       call MPI_SEND(nSquares, 1, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, ierr)
!!$       do j = 1, 4
!!$          call MPI_SEND(corners(1:nSquares,j), nSquares, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
!!$       enddo
!!$       call MPI_SEND(value(1:nSquares), nSquares, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
!!$       call MPI_SEND(speed(1:nSquares), nSquares, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
!!$       call MPI_SEND(ang(1:nSquares), nSquares, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
!!$    endif
!!$    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!!$  end subroutine getSquares
!!$
!!$
!!$  recursive subroutine recursGetSquares(thisOctal, grid, plane, valueName, nSquares, corners, value, speed, ang)
!!$
!!$    use mpi
!!$    type(gridtype) :: grid
!!$    type(octal), pointer   :: thisOctal
!!$    type(octal), pointer  :: child
!!$    !
!!$    integer :: subcell, i
!!$    character(len=*) :: valueName, plane
!!$    integer :: nSquares
!!$    real, pointer :: corners(:,:)
!!$    real, pointer :: value(:)
!!$    real, pointer :: speed(:)
!!$    real, pointer :: ang(:)
!!$    real(double) ::tmp
!!$    real(double) :: eps, x
!!$    integer :: myRank, ierr
!!$    type(VECTOR) :: rVec
!!$
!!$    eps = 0.01d0 * grid%halfSmallestSubcell
!!$
!!$    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
!!$
!!$    do subcell = 1, thisOctal%maxChildren
!!$
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call recursGetSquares(child, grid, plane, valueName, nSquares, corners, value, speed, ang)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$
!!$!          if ((grid%splitOverMpi).and.(thisOctal%mpiThread(subcell) /= myRank)) then
!!$!             cycle
!!$!          endif
!!$          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle
!!$
!!$          rVec = subcellCentre(thisOctal, subcell)
!!$          select case(valuename)
!!$             case("rho")
!!$                tmp = thisOctal%rho(subcell)
!!$             case("cs")
!!$                tmp = sqrt(thisOctal%pressure_i(subcell)/thisOctal%rho(subcell))
!!$             case("mass")
!!$                tmp = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell) * 1.d30 / msol
!!$             case("pressure")
!!$                tmp = thisOctal%pressure_i(subcell)
!!$             case("u")
!!$                tmp = thisOctal%rhou(subcell)/thisOctal%rho(subcell)/1.d5
!!$             case("v")
!!$                tmp = thisOctal%rhov(subcell)/thisOctal%rho(subcell)/1.d5
!!$             case("w")
!!$                tmp = thisOctal%rhow(subcell)/thisOctal%rho(subcell)/1.d5
!!$             case("rhou")
!!$                tmp = thisOctal%rhou(subcell)
!!$             case("rhov")
!!$                tmp = thisOctal%rhov(subcell)
!!$             case("rhow")
!!$                tmp = thisOctal%rhow(subcell)
!!$             case("rhoe")
!!$                tmp = thisOctal%rhoe(subcell)
!!$             case("ionization")
!!$                tmp = max(1.d-20,thisOctal%ionFrac(subcell,1))
!!$             case("coeff")
!!$                tmp = max(1.d-20,thisOctal%photoionCoeff(subcell,1))
!!$             case("crossings")
!!$                tmp = thisOctal%nCrossings(subcell)
!!$             case("temperature")
!!$                tmp = thisOctal%temperature(subcell)
!!$             case("mpi")
!!$                tmp = real(thisOctal%mpiThread(subcell))
!!$             case("phi")
!!$                tmp = real(thisOctal%phi_i(subcell))
!!$             case("chi")
!!$                tmp = real(thisOctal%chiline(subcell))
!!$             case("u_i")
!!$                tmp = real(thisOctal%u_interface(subcell))
!!$             case("q_i")
!!$                tmp = real(thisOctal%q_i(subcell))
!!$             case("q_i_minus_1")
!!$                tmp = real(thisOctal%q_i_minus_1(subcell))
!!$             case("flux_i")
!!$                tmp = real(thisOctal%flux_i(subcell))
!!$             case("philim")
!!$                tmp = real(thisOctal%phiLimit(subcell))
!!$             case DEFAULT
!!$           end select
!!$           if (thisOctal%oneD) then
!!$              nSquares = nSquares + 1
!!$              corners(nsquares, 1) = rVec%x
!!$              value(nSquares) = tmp
!!$           else
!!$
!!$              select case(plane)
!!$              case("x-z")
!!$                 x = min(abs(rVec%y + thisOctal%subcellSize/2.d0), abs(rVec%y - thisOctal%subcellSize/2.d0))
!!$                 if ((x < eps).or.(thisOctal%twoD)) then
!!$                    nSquares = nSquares + 1
!!$                    corners(nSquares, 1) = rVec%x - thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 2) = rVec%x + thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 3) = rVec%z - thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 4) = rVec%z + thisOctal%subcellSize/2.d0
!!$                    value(nSquares) = tmp
!!$                    speed(nSquares) = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhow(subcell)**2) / thisOctal%rho(subcell)
!!$                    ang(nSquares) = atan2(thisOctal%rhow(subcell), thisOctal%rhou(subcell))
!!$                 endif
!!$              case("y-z")
!!$                 x = min(abs(rVec%x + thisOctal%subcellSize/2.d0), abs(rVec%x - thisOctal%subcellSize/2.d0))
!!$                 if (x < eps) then
!!$                    nSquares = nSquares + 1
!!$                    corners(nSquares, 1) = rVec%y - thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 2) = rVec%y + thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 3) = rVec%z - thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 4) = rVec%z + thisOctal%subcellSize/2.d0
!!$                    value(nSquares) = tmp
!!$                    speed(nSquares) = sqrt(thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2) / thisOctal%rho(subcell)
!!$                    ang(nSquares) = atan2(thisOctal%rhow(subcell), thisOctal%rhov(subcell))
!!$                 endif
!!$              case("x-y")
!!$                 x = min(abs(rVec%z + thisOctal%subcellSize/2.d0), abs(rVec%z - thisOctal%subcellSize/2.d0))
!!$                 if (x < eps) then
!!$                    nSquares = nSquares + 1
!!$                    corners(nSquares, 1) = rVec%x - thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 2) = rVec%x + thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 3) = rVec%y - thisOctal%subcellSize/2.d0
!!$                    corners(nSquares, 4) = rVec%y + thisOctal%subcellSize/2.d0
!!$                    value(nSquares) = tmp
!!$                    speed(nSquares) = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2) / thisOctal%rho(subcell)
!!$                    ang(nSquares) = atan2(thisOctal%rhov(subcell), thisOctal%rhou(subcell))
!!$                    if (thisOctal%ghostCell(subcell)) speed(nSquares) = 0.
!!$                 endif
!!$              case DEFAULT
!!$              end select
!!$           endif
!!$
!!$       endif
!!$    enddo
!!$  end subroutine recursGetSquares


  subroutine receiveAcrossMpiBoundary(grid, boundaryType, receiveThread, sendThread)
    use mpi
    use inputs_mod, only :  smallestCellSize

    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal, tOctal
    type(octal), pointer  :: neighbourOctal
    character(len=*) :: boundaryType
    integer :: receiveThread, sendThread, tsubcell
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: nOctals
    integer, parameter :: nStorage = 58
    type(VECTOR) :: octVec, direction, rVec, pVec,locator
    integer :: nBound
    integer :: iOctal
    integer :: subcell, neighbourSubcell
    integer :: tag = 77
    integer :: status(MPI_STATUS_SIZE)
    integer :: nDepth
    integer :: ierr
    integer :: i, j
    type(VECTOR), parameter :: xhat = VECTOR(1.d0, 0.d0, 0.d0), yHat = VECTOR(0.d0, 1.d0, 0.d0), zHat = VECTOR(0.d0, 0.d0, 1.d0)
    real(double) :: q , rho, rhoe, rhou, rhov, rhow, pressure, phi, flux, phigas
    real(double) :: temperature
    real(double) :: rm1, rum1, pm1, qViscosity(3,3)
    real(double), allocatable :: locStack(:), temp1d(:), tempStorage(:,:)
    integer, allocatable :: nDepthStorage(:)
    
    integer :: nLoc

    select case(boundaryType)
    case("left")
       direction = VECTOR(-1.d0, 0.d0, 0.d0)
       nBound = 1
    case("right")
       direction = VECTOR(1.d0, 0.d0, 0.d0)
       nBound = 2
    case("top")
       direction = VECTOR(0.d0, 0.d0, 1.d0)
       nBound = 3
    case("bottom")
       direction = VECTOR(0.d0, 0.d0, -1.d0)
       nBound = 4
    case("front")
       direction = VECTOR(0.d0, 1.d0, 0.d0)
       nBound = 5
    case("back")
       direction = VECTOR(0.d0, -1.d0, 0.d0)
       nBound = 6
    case DEFAULT
       write(*,*) "boundary type not recognised ",boundaryType
       stop
    end select
!    write(*,*) myrank, "boundary number is ",nbound
    

    if (myRankGlobal == receiveThread) then
       allocate(octalArray(grid%nOctals))
       nOctals = 0
       call getOctalArray(grid%octreeRoot,octalArray, nOctals)
       if (nOctals /= grid%nOctals) then
          write(*,*) "Screw up in get octal array", nOctals,grid%nOctals
          stop
       endif
!       write(*,*) myrank," generated ",nOctals, " array of octals"
       nLoc = 0
       allocate(locStack(1:nOctals*8*3))
       allocate(nDepthStorage(1:nOctals*8))
       do iOctal =  1, nOctals
          
          thisOctal => octalArray(iOctal)%content
          
          do subcell = 1, thisOctal%maxChildren             
             if (.not.thisOctal%hasChild(subcell)) then
                
                if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
                
                octVec = subcellCentre(thisOctal, subcell) + &
                     (thisOctal%subcellSize/2.d0+0.01d0 * smallestCellSize) * direction
                
                if (.not.inOctal(grid%octreeRoot, octVec)) then
                   write(*,*) "Grid doesn't have a ", boundaryType, " surface in this volume"
                   write(*,*) "centre",subcellCentre(thisOctal, subcell)
                   write(*,*) "depth",thisOctal%nDepth, thisOctal%haschild(1:8)
                   write(*,*) octVec
                   stop
                endif
                
                neighbourOctal => thisOctal
                call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)

                if (octalOnThread(neighbourOctal, neighbourSubcell, receiveThread)) cycle

                if (.not.octalOnThread(neighbourOctal, neighbourSubcell, sendThread)) then
                   write(*,*) "Error 1 Neighbour on ",boundaryType, " of ", myrankglobal, &
                        "  is not on thread ", sendThread, " but ", &
                   neighbourOctal%mpiThread(neighboursubcell), " depth ",neighbourOctal%ndepth
                   stop
                endif

                nLoc = nLoc + 1
                locStack(nLoc*3-2) = octVec%x
                locStack(nLoc*3-1) = octVec%y
                locStack(nLoc*3  ) = octVec%z
                nDepthStorage(nLoc) = thisOctal%nDepth


!                write(*,*) myRank, " has identified a boundary cell ", loc(1:3)



!                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

!                call MPI_SEND(thisOctal%nDepth, 1, MPI_INTEGER, sendThread, tag, localWorldCommunicator, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

!               call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
!                write(*,*) myrank, " received temp storage"
                if (associated(thisOctal%mpiBoundaryStorage).and.(size(thisOctal%mpiBoundaryStorage,3)/=nStorage)) then
                   deallocate(thisOctal%mpiBoundaryStorage)
                   allocate(thisOctal%mpiBoundaryStorage(1:thisOctal%maxChildren, 6, nStorage))
                   thisOctal%mpiBoundaryStorage = 0.d0
                endif

                if (.not.associated(thisOctal%mpiBoundaryStorage)) then
                   allocate(thisOctal%mpiBoundaryStorage(1:thisOctal%maxChildren, 6, nStorage))
                   thisOctal%mpiBoundaryStorage = 0.d0
                endif
!                write(*,*) myrank, " successfully stored"

             end if
          enddo
       enddo

       call MPI_SEND(nLoc, 1, MPI_INTEGER, sendThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(locStack, 3*nLoc, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(nDepthStorage, nLoc, MPI_INTEGER, sendThread, tag, localWorldCommunicator, ierr)

       allocate(tempStorage(nLoc, nStorage))
       allocate(temp1d(nLoc*nStorage))
       call MPI_RECV(temp1d, nLoc*nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
       tempStorage = RESHAPE(temp1d, shape=SHAPE(tempStorage))
       deallocate(temp1d)

       nLoc = 0

       do iOctal =  1, nOctals
          
          thisOctal => octalArray(iOctal)%content
          
          do subcell = 1, thisOctal%maxChildren             
             if (.not.thisOctal%hasChild(subcell)) then
                
                if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
                
                octVec = subcellCentre(thisOctal, subcell) + &
                     (thisOctal%subcellSize/2.d0+0.01d0 * smallestCellSize) * direction
                
                if (.not.inOctal(grid%octreeRoot, octVec)) then
                   write(*,*) "Grid doesn't have a ", boundaryType, " surface in this volume"
                   write(*,*) "centre",subcellCentre(thisOctal, subcell)
                   write(*,*) "depth",thisOctal%nDepth, thisOctal%haschild(1:8)
                   write(*,*) octVec
                   stop
                endif
                
                neighbourOctal => thisOctal
                call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)

                if (octalOnThread(neighbourOctal, neighbourSubcell, receiveThread)) cycle

                if (.not.octalOnThread(neighbourOctal, neighbourSubcell, sendThread)) then
                   write(*,*) "Error 1 Neighbour on ",boundaryType, " of ", myrankglobal, &
                        "  is not on thread ", sendThread, " but ", &
                   neighbourOctal%mpiThread(neighboursubcell), " depth ",neighbourOctal%ndepth
                   stop
                endif

                nLoc = nLoc + 1
                thisOctal%mpiBoundaryStorage(subcell, nBound, 1:nStorage) = tempStorage(nLoc,1:nStorage)

             end if
          enddo
       enddo

       deallocate(tempStorage)
       deallocate(octalArray)


       ! now send a finish signal to the sendThread
    else
       if (myRankGlobal /= sendThread) then
          write(*,*) "subroutine called within thread ", myRankGlobal, " but expecing to be ", sendthread
          stop
       endif

       call MPI_RECV(nLoc, 1, MPI_INTEGER, receiveThread, tag, localWorldCommunicator, status, ierr)
       allocate(locStack(1:nLoc*3))
       allocate(nDepthStorage(1:nLoc))
       allocate(tempStorage(1:nLoc,1:nStorage))
       call MPI_RECV(locStack, 3*nLoc, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, status, ierr)
       call MPI_RECV(nDepthStorage, nLoc, MPI_INTEGER, receiveThread, tag, localWorldCommunicator, status, ierr)

       do i = 1, nLoc
          octVec = VECTOR(locStack(i*3-2), locStack(i*3-1), locStack(i*3))
          nDepth = nDepthStorage(i)


          call findSubcellTD(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell)

          pVec = subcellCentre(neighbourOctal, neighbourSubcell)

          tempStorage(i,15) = pVec%x
          tempStorage(i,16) = pVec%y
          tempStorage(i,17) = pVec%z
          
          if (neighbourOctal%mpiThread(neighboursubcell) /= sendthread) then
             write(*,*) "trying to send on ",boundaryType, " but is not on thread ", sendThread
             write(*,*) "rather is on ", neighbourOctal%mpiThread(neighboursubcell)
             stop
          endif
          
          if (neighbourOctal%nDepth <= nDepth) then
             if ((nDepth - neighbourOctal%nDepth) > 1) then
                write(*,*) "Octal depth differs by more than 1 across boundary!!!"
                write(*,*) "ndepth ",nDepth, " neighbour%nDepth ",neighbourOctal%nDepth
                write(*,*) "myrank ",myrankGlobal
                write(*,*) "sendThread ",sendThread, " receivethread ",receivethread
                write(*,*) "at ", pVec%x, pVec%y, pVec%z
                stop
                
             endif
             tempStorage(i,1) = neighbourOctal%q_i(neighbourSubcell)
             tempStorage(i,2) = neighbourOctal%rho(neighbourSubcell)
             tempStorage(i,3) = neighbourOctal%rhoe(neighbourSubcell)
             tempStorage(i,4) = neighbourOctal%rhou(neighbourSubcell)
             tempStorage(i,5) = neighbourOctal%rhov(neighbourSubcell)
             tempStorage(i,6) = neighbourOctal%rhow(neighbourSubcell)
             tempStorage(i,7) = neighbourOctal%x_i(neighbourSubcell)
             rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
                  direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*smallestCellSize)
             tOctal => neighbourOctal
             tSubcell = neighbourSubcell
             call findSubcellLocal(rVec, tOctal, tSubcell)
             tempStorage(i,8) = tOctal%q_i(tsubcell)
             
             tempStorage(i,9) = dble(neighbourOctal%nDepth)
             tempStorage(i,10) = neighbourOctal%pressure_i(neighbourSubcell)
             tempStorage(i,11) = neighbourOctal%flux_i(neighbourSubcell)
             
             tempStorage(i,12) = neighbourOctal%phi_i(neighbourSubcell)
             
             tempStorage(i,13) = neighbourOctal%phi_gas(neighbourSubcell)
             tempStorage(i,14) = neighbourOctal%x_i_minus_1(neighbourSubcell)
             
             
             
             tempStorage(i,21) = neighbourOctal%rho_i_minus_1(neighbourSubcell)
             tempStorage(i,22) = neighbourOctal%u_i_minus_1(neighbourSubcell)
             tempStorage(i,23) = neighbourOctal%pressure_i_minus_1(neighbourSubcell)
             
             
             tempStorage(i,24) = neighbourOctal%qViscosity(neighbourSubcell,1,1)
             tempStorage(i,25) = neighbourOctal%qViscosity(neighbourSubcell,1,2)
             tempStorage(i,26) = neighbourOctal%qViscosity(neighbourSubcell,1,3)
             tempStorage(i,27) = neighbourOctal%qViscosity(neighbourSubcell,2,1)
             tempStorage(i,28) = neighbourOctal%qViscosity(neighbourSubcell,2,2)
             tempStorage(i,29) = neighbourOctal%qViscosity(neighbourSubcell,2,3)
             tempStorage(i,30) = neighbourOctal%qViscosity(neighbourSubcell,3,1)
             tempStorage(i,31) = neighbourOctal%qViscosity(neighbourSubcell,3,2)
             tempStorage(i,32) = neighbourOctal%qViscosity(neighbourSubcell,3,3)
             
             tempStorage(i,33) = neighbourOctal%temperature(neighbourSubcell)
             tempStorage(i,34) = neighbourOctal%flux_amr_i(neighbourSubcell,1)
             tempStorage(i,35) = neighbourOctal%flux_amr_i(neighbourSubcell,2)
             tempStorage(i,36) = neighbourOctal%flux_amr_i(neighbourSubcell,3)
             tempStorage(i,37) = neighbourOctal%flux_amr_i(neighbourSubcell,4)
             tempStorage(i,58) = dble(neighbourOctal%nChildren)

          else ! need to average (neighbour octal depth > this Octal depth)
             call averageValue(direction, neighbourOctal,  neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, &
                  flux, phi, phigas, rm1, rum1, pm1, qViscosity, temperature)
             tempStorage(i,1) = q
             tempStorage(i,2) = rho
             tempStorage(i,3) = rhoe
             tempStorage(i,4) = rhou
             tempStorage(i,5) = rhov
             tempStorage(i,6) = rhow
             tempStorage(i,7) = neighbourOctal%x_i(neighbourSubcell)
             rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
                  direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*smallestCellSize)
             tOctal => neighbourOctal
             tSubcell = neighbourSubcell
             call findSubcellLocal(rVec, tOctal, tSubcell)
             tempStorage(i,8) = tOctal%q_i(tsubcell)
             
             tempStorage(i,9) = dble(neighbourOctal%nDepth)
             tempStorage(i,10) = pressure
             tempStorage(i,11) = flux
             tempStorage(i,12) = phi
             tempStorage(i,13) = phigas
             tempstorage(i,14) = neighbourOctal%x_i_minus_1(neighbourSubcell)
             tempStorage(i,21) = rm1
             tempStorage(i,22) = rum1
             tempStorage(i,23) = pm1
             
             
             tempStorage(i,24) = qViscosity(1,1)
             tempStorage(i,25) = qViscosity(1,2)
             tempStorage(i,26) = qViscosity(1,3)
             tempStorage(i,27) = qViscosity(2,1)
             tempStorage(i,28) = qViscosity(2,2)
             tempStorage(i,29) = qViscosity(2,3)
             tempStorage(i,30) = qViscosity(3,1)
             tempStorage(i,31) = qViscosity(3,2)
             tempStorage(i,32) = qViscosity(3,3)
             
             tempStorage(i,33) = temperature
             tempStorage(i,34) = neighbourOctal%flux_amr_i(neighbourSubcell,1)
             tempStorage(i,35) = neighbourOctal%flux_amr_i(neighbourSubcell,2)
             tempStorage(i,36) = neighbourOctal%flux_amr_i(neighbourSubcell,3)
             tempStorage(i,37) = neighbourOctal%flux_amr_i(neighbourSubcell,4)
             
             tempStorage(i,58) = dble(neighbourOctal%nChildren)


          do j = 1, 4

! +ve x face                                                                                                                                                                                  
! -z -y 2  (1)                                                                                                                                                                                
! -z +y 4  (2)                                                                                                                                                                                
! +z -y 6  (3)                                                                                                                                                                                
! +z +y 8  (4)                                                                                                                                                                                
             if (abs(direction%x) > 0.9d0) then
                select case(j)
                   case(1)
                      locator = octVec &
                           - 0.5d0 * neighbourOctal%subcellSize * zHat &
                           - 0.5d0 * neighbourOctal%subcellSize * yHat
                   case(2)
                      locator = octVec &
                           - (0.5d0 * neighbourOctal%subcellSize) * zHat &
                           + (0.5d0 * neighbourOctal%subcellSize) * yHat
                   case(3)
                      locator = octVec &
                           + (0.5d0 * neighbourOctal%subcellSize) * zHat &
                           - (0.5d0 * neighbourOctal%subcellSize) * yHat
                   case(4)
                      locator = octVec &
                           + (0.5d0 * neighbourOctal%subcellSize) * zHat &
                           + (0.5d0 * neighbourOctal%subcellSize) * yHat
                end select
             endif
! -ve y face                                                                                                                                                                                  
! -x -z 1  (1)                                                                                                                                                                                
! -x +z 5  (2)                                                                                                                                                                                
! +x -z 2  (3)                                                                                                                                                                                
! +x +z 6  (4)                                                                                                                                                                                
             if (abs(direction%y) > 0.9d0) then
                select case(j)
                   case(1)
                      locator = octVec &
                           - (0.5d0 * neighbourOctal%subcellSize) * xHat &
                           - (0.5d0 * neighbourOctal%subcellSize) * zHat
                   case(2)
                      locator = octVec &
                           - (0.5d0 * neighbourOctal%subcellSize) * xHat &
                           + (0.5d0 * neighbourOctal%subcellSize) * zHat
                   case(3)
                      locator = octVec &
                           + (0.5d0 * neighbourOctal%subcellSize) * xHat &
                           - (0.5d0 * neighbourOctal%subcellSize) * zHat
                   case(4)
                      locator = octVec &
                           + (0.5d0 * neighbourOctal%subcellSize) * xHat &
                           + (0.5d0 * neighbourOctal%subcellSize) * zHat
                end select
             endif
! -ve z face                                                                                                                                                                                  
! -x -y 1  (1)                                                                                                                                                                                
! -x +y 3  (2)                                                                                                                                                                                
! +x -y 2  (3)                                                                                                                                                                                
! +x +y 4  (4)                                                                                                                                                                                
             if (abs(direction%z) > 0.9d0) then
                select case(j)
                   case(1)
                      locator = octVec &
                           - (0.5d0 * neighbourOctal%subcellSize) * xHat &
                           - (0.5d0 * neighbourOctal%subcellSize) * yHat
                   case(2)
                      locator = octVec &
                           - (0.5d0 * neighbourOctal%subcellSize) * xHat &
                           + (0.5d0 * neighbourOctal%subcellSize) * yHat
                   case(3)
                      locator = octVec &
                           + (0.5d0 * neighbourOctal%subcellSize) * xHat &
                           - (0.5d0 * neighbourOctal%subcellSize) * yHat
                   case(4)
                      locator = octVec &
                           + (0.5d0 * neighbourOctal%subcellSize) * xHat &
                           + (0.5d0 * neighbourOctal%subcellSize) * yHat
                end select
             endif




             tOctal => neighbourOctal
             call findSubcellLocal(locator, tOctal, tSubcell)

             tempStorage(i,42+j-1) = tOctal%q_i(tSubcell)
             tempStorage(i,38+j-1) = tOctal%rho(tSubcell)
             tempStorage(i,46+j-1) = tOctal%rhou(tSubcell)
             tempStorage(i,50+j-1) = tOctal%rhov(tSubcell)
             tempStorage(i,54+j-1) = tOctal%rhow(tSubcell)
             tempStorage(i,34+j-1) = tOctal%flux_amr_i(tSubcell,j)

          enddo
       endif
       

    enddo
       
    allocate(temp1d(nStorage*nLoc))
    temp1d = RESHAPE(tempStorage, SHAPE=SHAPE(temp1d))

    call MPI_SEND(tempStorage, nStorage*nLoc, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, ierr)
    deallocate(temp1d, tempStorage)
    endif
  end subroutine receiveAcrossMpiBoundary

  subroutine receiveAcrossMpiCorner(grid, boundaryType, receiveThread, sendThread)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: neighbourOctal
    character(len=*) :: boundaryType
    integer :: receiveThread, sendThread
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: nOctals
    integer, parameter :: nStorage = 4
    real(double) :: loc(3), tempStorage(nStorage)
    type(VECTOR) :: octVec, direction,  pVec
    integer :: nBound
    integer :: iOctal
    integer :: subcell, neighbourSubcell
    integer :: tag = 77
    integer :: status(MPI_STATUS_SIZE)
    logical :: sendLoop
    integer :: nDepth
    integer :: ierr
!    real(double) :: q , rho, rhoe, rhov, rhow, pressure, phi, flux, phigas
!
    select case(boundaryType)
    case("leftupper")
       direction = VECTOR(-1.d0, 0.d0, 1.d0)
       nBound = 1
    case("rightlower")
       direction = VECTOR(1.d0, 0.d0, -1.d0)
       nBound = 2
    case("rightupper")
       direction = VECTOR(1.d0, 0.d0, 1.d0)
       nBound = 3
    case("leftlower")
       direction = VECTOR(-1.d0, 0.d0, -1.d0)
       nBound = 4
    case("topupper")
       direction = VECTOR(1.d0, 1.d0, 1.d0)
       nBound = 5
    case("bottomlower")
       direction = VECTOR(-1.d0, -1.d0, -1.d0)
       nBound = 6
    case("toplower")
       direction = VECTOR(1.d0, -1.d0, 1.d0)
       nBound = 7
    case("bottomupper")
       direction = VECTOR(-1.d0, 1.d0, -1.d0)
       nBound = 8
    case DEFAULT
       write(*,*) "boundary type not recognised ",boundaryType
       stop
    end select
!    write(*,*) myrank, "boundary number is ",nbound
    
    if (myRankGlobal == receiveThread) then
       allocate(octalArray(grid%nOctals))
       nOctals = 0
       call getOctalArray(grid%octreeRoot,octalArray, nOctals)
       if (nOctals /= grid%nOctals) then
          write(*,*) "Screw up in get octal array", nOctals,grid%nOctals
          stop
       endif
!       write(*,*) myrank," generated ",nOctals, " array of octals"
       do iOctal =  1, nOctals
          
          thisOctal => octalArray(iOctal)%content
          
          do subcell = 1, thisOctal%maxChildren             
             if (.not.thisOctal%hasChild(subcell)) then
                
!                if (thisOctal%mpiThread(subcell) /= myRank) cycle
                if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
                
                octVec = subcellCentre(thisOctal, subcell) + &
                     ((sqrt(2.d0)*(thisOctal%subcellSize/2.d0))+0.01d0 * grid%halfSmallestSubcell) * direction
                
                if (.not.inOctal(grid%octreeRoot, octVec)) cycle !then
!                   write(*,*) "Grid doesn't have a ", boundaryType, " surface in this volume"
!                   write(*,*) "centre",subcellCentre(thisOctal, subcell)
!                   write(*,*) "depth",thisOctal%nDepth, thisOctal%haschild(1:8)
!                   write(*,*) "direction ", direction
!                   write(*,*) octVec
!                   stop
!                endif
                
                neighbourOctal => thisOctal
                call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)

                if (octalOnThread(neighbourOctal, neighbourSubcell, receiveThread)) cycle

                if(neighbourOctal%mpiThread(subcell) /= sendThread) cycle

                if (.not.octalOnThread(neighbourOctal, neighbourSubcell, sendThread)) then
                   write(*,*) "Error 1 Neighbour on ",boundaryType, " of ", myrankglobal, &
                        "  is not on thread ", sendThread, " but ", &
                   neighbourOctal%mpiThread(neighboursubcell), " depth ",neighbourOctal%ndepth
                   stop
                endif

                loc(1) = octVec%x
                loc(2) = octVec%y
                loc(3) = octVec%z

                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_SEND(thisOctal%nDepth, 1, MPI_INTEGER, sendThread, tag, localWorldCommunicator, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
!                write(*,*) myrank, " received temp storage"
                if (.not.associated(thisOctal%mpiCornerStorage)) then
                   allocate(thisOctal%mpiCornerStorage(1:thisOctal%maxChildren, 8, nStorage))
                   thisOctal%mpiCornerStorage = 0.d0
                endif
                thisOctal%mpiCornerStorage(subcell, nBound, 1:nStorage) = tempStorage(1:nStorage)
!                write(*,*) myrank, " successfully stored"

             end if
          enddo
       enddo
       deallocate(octalArray)
       loc(1) = HUGE(loc(1))
       loc(2) = 0.d0
       loc(3) = 0.d0
!       write(*,*) myrank, " sending a huge value to ", sendThread
       call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)

       ! now send a finish signal to the sendThread
    else
       if (myRankGlobal /= sendThread) then
          write(*,*) "subroutine called within thread ", myRankGlobal, " but expecing to be ", sendthread
          stop
       endif
       sendLoop = .true.
       do while (sendLoop)
          ! receive a locator
!          write(*,*) myrank, " waiting for a locator"
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, status, ierr)
!          write(*,*) myrank, " received a locator ", loc(1:3)
          if (loc(1) > 1.d20) then
             sendLoop = .false.
!             write(*,*) myRank, " found the signal to end the send loop"
          else

             octVec = VECTOR(loc(1), loc(2), loc(3))

             call MPI_RECV(nDepth, 1, MPI_INTEGER, receiveThread, tag, localWorldCommunicator, status, ierr)

             call findSubcellTD(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell)

             pVec = subcellCentre(neighbourOctal, neighbourSubcell)

             tempStorage(1) = pVec%x
             tempStorage(2) = pVec%y
             tempStorage(3) = pVec%z
             tempStorage(4) = neighbourOctal%flux_i(neighbourSubcell)

             if (neighbourOctal%mpiThread(neighboursubcell) /= sendthread) then
                write(*,*) "trying to send on ",boundaryType, " but is not on thread ", sendThread
                stop
             endif
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, ierr)
            
          endif
       end do
      endif
  end subroutine receiveAcrossMpiCorner
  
  subroutine exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: iPair, nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup, iGroup
    integer, optional :: useThisBound
    integer :: rBound
    logical :: doExchange
    character(len=10) :: boundaryType(6) = (/"left  ","right ", "top   ", "bottom", "front ", "back  "/)
    integer :: ierr
    if (myrankGlobal == 0) goto 666
    CALL MPI_BARRIER(amrCOMMUNICATOR, ierr)

!    print *, "nBound", nBound, myRank
!    print *, "nGroup", nGroup, myRank

    do iGroup = 1, nGroup

       do iPair = 1, nPairs

          if (group(iPair) == iGroup) then
             doExchange = .true.
             if (present(useThisBound)) then
                doExchange = .false.
                if (nBound(iPair) == useThisBound) doExchange = .true.
             endif
             if (doExchange) then
                !       write(*,*) "myrank ", iPair, thread1(iPair), thread2(iPair)
                if ((myRankGlobal == thread1(iPair)).or.(myRankGlobal == thread2(iPair))) then
                   !          write(*,*) myrank, " calling receive across boundary",iPair,nbound(ipair),boundaryType(nBound(iPair))
                   call receiveAcrossMpiBoundary(grid, boundaryType(nBound(iPair)), thread1(iPair), thread2(iPair))
                   !          write(*,*) myrank, " finished receive across boundary"
                   if      (nBound(iPair) == 1) then
                      rBound = 2
                   else if (nBound(iPair) == 2) then
                      rBound = 1
                   else if (nBound(iPair) == 3) then
                      rBound = 4
                   else if (nBound(iPair) == 4) then
                      rBound = 3
                   else if (nBound(iPair) == 5) then
                      rBound = 6
                   else if (nBound(iPair) == 6) then
                      rBound = 5
                   endif
                   call receiveAcrossMpiBoundary(grid, boundaryType(rBound), thread2(iPair), thread1(iPair))
                   !          if (myRank == 1) then
                   !             write(*,*) "Exchange done for direction ",nbound(ipair), " and ",rBound
                   !             write(*,*) "Between threads ", thread1(ipair), " and ", thread2(ipair)
                   !          endif
                endif
             endif
          endif
       enddo
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    enddo
666 continue
  end subroutine exchangeAcrossMPIboundary

  subroutine exchangeAcrossMPIcorner(grid, nCornerPairs, cornerThread1, cornerThread2, nCornerBound, cornerGroup, nCornerGroup, &
      useThisBound)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: iPair, nCornerPairs, cornerThread1(:), cornerThread2(:), nCornerBound(:)
    integer :: cornerGroup(:), nCornerGroup, iGroup
    integer, optional :: useThisBound
    integer :: rBound!, cOne, cTwo
    logical :: doExchange
    integer :: ierr
    character(len=12) :: cornerType(8) = (/"leftupper   ","rightlower  ", "rightupper  ", "leftlower   ", &
      "topupper    ", "bottomlower ", "toplower    ", "bottomupper " /)

    if (myrankGlobal == 0) goto 666
    CALL MPI_BARRIER(amrCOMMUNICATOR, ierr)

!      print *, "nCornerBound", nCornerBound, myRank
!      print *, "nCornerGroup", nCornerGroup, myRank

    do iGroup = 1, nCornerGroup
       do iPair = 1, nCornerPairs
          if (cornerGroup(iPair) == iGroup) then
             doExchange = .true.
             if (present(useThisBound)) then
                doExchange = .false.
                if (nCornerBound(iPair) == useThisBound) doExchange = .true.
             endif
             if (doExchange) then
                if ((myRankGlobal == cornerThread1(iPair)).or.(myRankGlobal == cornerThread2(iPair))) then

                   call receiveAcrossMpiCorner(grid, cornerType(nCornerBound(iPair)), &
                      cornerThread1(iPair), cornerThread2(iPair))
                   if      (nCornerBound(iPair) == 1) then
                      rBound = 2
                   else if (nCornerBound(iPair) == 2) then
                      rBound = 1
                   else if (nCornerBound(iPair) == 3) then
                      rBound = 4
                   else if (nCornerBound(iPair) == 4) then
                      rBound = 3
                   else if (nCornerBound(iPair) == 5) then
                      rBound = 6
                   else if (nCornerBound(iPair) == 6) then
                      rBound = 5
                   else if (nCornerBound(iPair) == 7) then
                      rBound = 8
                   else if (nCornerBound(iPair) == 8) then
                      rBound = 7
                   endif
                   call receiveAcrossMpiCorner(grid, cornerType(rBound), &
                      cornerThread2(iPair), cornerThread1(iPair))
!                             if (myRank == 1) then
!                                write(*,*) "Exchange done for direction ",nCornerbound(ipair), " and ",rBound
!                                write(*,*) "Between threads ", cornerThread1(ipair), " and ", cornerThread2(ipair)
!                             endif
                endif
             endif
          endif
       enddo
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    enddo
666 continue
  end subroutine exchangeAcrossMPIcorner


  subroutine receiveAcrossMpiBoundaryLevel(grid, boundaryType, receiveThread, sendThread, nDepth)

    use mpi
    !    use inputs_mod, only : useTensorViscosity
    integer :: ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal, tOctal
    type(octal), pointer  :: neighbourOctal
    character(len=*) :: boundaryType
    integer :: receiveThread, sendThread, tsubcell
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: nOctals
    integer, parameter :: nStorage = 58
    integer :: nLoc
    real(double),allocatable :: tempStorage(:,:),temp1d(:)
    real(double), allocatable :: locStack(:)
    integer :: i
    real(double) :: qViscosity(3,3)
    type(VECTOR) :: octVec, direction, rVec, pVec
    integer :: nBound
    integer :: iOctal
    integer :: nDepth !, thisnDepth
    integer :: subcell, neighbourSubcell
    integer :: tag = 77
    integer :: status(MPI_STATUS_SIZE)

    select case(boundaryType)
    case("left")
       direction = VECTOR(-1.d0, 0.d0, 0.d0)
       nBound = 1
    case("right")
       direction = VECTOR(1.d0, 0.d0, 0.d0)
       nBound = 2
    case("top")
       direction = VECTOR(0.d0, 0.d0, 1.d0)
       nBound = 3
    case("bottom")
       direction = VECTOR(0.d0, 0.d0, -1.d0)
       nBound = 4
    case("front")
       direction = VECTOR(0.d0, 1.d0, 0.d0)
       nBound = 5
    case("back")
       direction = VECTOR(0.d0, -1.d0, 0.d0)
       nBound = 6
    case DEFAULT
       write(*,*) "boundary type not recognised ",boundaryType
       stop
    end select
    !    write(*,*) myrank, "boundary number is ",nbound


    if (myRankGlobal == receiveThread) then
       allocate(octalArray(grid%nOctals))
       nOctals = 0
       call getOctalArrayLevel(grid%octreeRoot,octalArray, nOctals, nDepth)
       !       write(*,*) myrank," generated ",nOctals, " array of octals"

       allocate(locStack(1:nOctals*8*3))
       nLoc = 0

       do iOctal =  1, nOctals

          thisOctal => octalArray(iOctal)%content
          !          write(*,*) myrank, " doing octal of ", thisOctal%ndepth, " depth "

          do subcell = 1, thisOctal%maxChildren


             if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

             octVec = subcellCentre(thisOctal, subcell) + &
                  (thisOctal%subcellSize/2.d0+0.01d0 * grid%halfSmallestSubcell) * direction

             if (.not.inOctal(grid%octreeRoot, octVec)) then
                write(*,*) "Grid doesn't have a ", boundaryType, " surface in this volume"
                write(*,*) "centre",subcellCentre(thisOctal, subcell)
                write(*,*) "depth",thisOctal%nDepth, thisOctal%haschild(1:8)
                write(*,*) octVec
                stop
             endif

             neighbourOctal => grid%octreeRoot
             call findSubcellLocalLevel(octVec, neighbourOctal, neighbourSubcell, nDepth)

             if (octalOnThread(neighbourOctal, neighbourSubcell, receiveThread)) cycle

             if (.not.octalOnThread(neighbourOctal, neighbourSubcell, sendThread)) then
                write(*,*) "Error 2 Neighbour on ",boundaryType, " of ", myrankglobal, &
                     "  is not on thread ", sendThread, " but ", &
                     neighbourOctal%mpiThread(neighboursubcell), " depth ",neighbourOctal%ndepth
                stop
             endif

             nLoc = nLoc  +  1
             locStack(nLoc*3-2) = octVec%x
             locStack(nLoc*3-1) = octVec%y
             locStack(nLoc*3  ) = octVec%z

             if (.not.associated(thisOctal%mpiBoundaryStorage)) then
                allocate(thisOctal%mpiBoundaryStorage(1:thisOctal%maxChildren, 6, nStorage))
                thisOctal%mpiBoundaryStorage = 0.d0
             endif

          enddo
       enddo

       call MPI_SEND(nLoc, 1, MPI_INTEGER, sendThread, tag, localWorldCommunicator, ierr)
       call MPI_SEND(locStack, 3*nLoc, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)

       allocate(tempStorage(nLoc, nStorage))
       allocate(temp1d(nLoc*nStorage))
       call MPI_RECV(temp1d, nLoc*nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
       tempStorage = RESHAPE(temp1d, shape=SHAPE(tempStorage))
       deallocate(temp1d)

       nLoc = 0

       do iOctal =  1, nOctals

          thisOctal => octalArray(iOctal)%content
          !          write(*,*) myrank, " doing octal of ", thisOctal%ndepth, " depth "

          do subcell = 1, thisOctal%maxChildren


             if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

             octVec = subcellCentre(thisOctal, subcell) + &
                  (thisOctal%subcellSize/2.d0+0.01d0 * grid%halfSmallestSubcell) * direction

             if (.not.inOctal(grid%octreeRoot, octVec)) then
                write(*,*) "Grid doesn't have a ", boundaryType, " surface in this volume"
                write(*,*) "centre",subcellCentre(thisOctal, subcell)
                write(*,*) "depth",thisOctal%nDepth, thisOctal%haschild(1:8)
                write(*,*) octVec
                stop
             endif

             neighbourOctal => grid%octreeRoot
             call findSubcellLocalLevel(octVec, neighbourOctal, neighbourSubcell, nDepth)

             if (octalOnThread(neighbourOctal, neighbourSubcell, receiveThread)) cycle

             if (.not.octalOnThread(neighbourOctal, neighbourSubcell, sendThread)) then
                write(*,*) "Error 2 Neighbour on ",boundaryType, " of ", myrankglobal, &
                     "  is not on thread ", sendThread, " but ", &
                     neighbourOctal%mpiThread(neighboursubcell), " depth ",neighbourOctal%ndepth
                stop
             endif

             if (.not.associated(thisOctal%mpiBoundaryStorage)) then
                allocate(thisOctal%mpiBoundaryStorage(1:thisOctal%maxChildren, 6, nStorage))
                thisOctal%mpiBoundaryStorage = 0.d0
             endif
             nLoc = nLoc + 1
             thisOctal%mpiBoundaryStorage(subcell, nBound, 1:nStorage) = tempStorage(nLoc,1:nStorage)
          enddo
       enddo

    else


       call MPI_RECV(nLoc, 1, MPI_INTEGER, receiveThread, tag, localWorldCommunicator, status, ierr)
       allocate(locStack(1:nLoc*3))
       allocate(tempStorage(1:nLoc,1:nStorage))
       call MPI_RECV(locStack, 3*nLoc, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, status, ierr)

       do i = 1, nLoc
          octVec = VECTOR(locStack(i*3-2), locStack(i*3-1), locStack(i*3))


          call findSubcellTDLevel(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell, nDepth)

          pVec = subcellCentre(neighbourOctal, neighbourSubcell)

          if (neighbourOctal%mpiThread(neighboursubcell) /= sendthread) then
             write(*,*) "trying to send on ",boundaryType, " but is not on thread ", sendThread
             stop
          endif

          !             if (neighbourOctal%nDepth <= thisnDepth) then
          Tempstorage(i,1) = neighbourOctal%q_i(neighbourSubcell)
          tempStorage(i,2) = neighbourOctal%rho(neighbourSubcell)
          tempStorage(i,3) = neighbourOctal%rhoe(neighbourSubcell)
          tempStorage(i,4) = neighbourOctal%rhou(neighbourSubcell)
          tempStorage(i,5) = neighbourOctal%rhov(neighbourSubcell)
          tempStorage(i,6) = neighbourOctal%rhow(neighbourSubcell)
          tempStorage(i,7) = neighbourOctal%x_i(neighbourSubcell)
          rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
               direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
          tOctal => grid%octreeRoot
          tSubcell = neighbourSubcell
          call findSubcellLocalLevel(rVec, tOctal, tSubcell, nDepth)
          tempStorage(i,8) = tOctal%q_i(tsubcell)

          tempStorage(i,9) = dble(neighbourOctal%nDepth)
          tempStorage(i,10) = neighbourOctal%pressure_i(neighbourSubcell)
          tempStorage(i,11) = neighbourOctal%flux_i(neighbourSubcell)


          tempStorage(i,12) = neighbourOctal%phi_i(neighbourSubcell)

          tempStorage(i,13) = neighbourOctal%phi_gas(neighbourSubcell)

          tempStorage(i,15) = pVec%x
          tempStorage(i,16) = pVec%y
          tempStorage(i,17) = pVec%z

          tempStorage(i,21) = neighbourOctal%rho_i_minus_1(neighbourSubcell)
          tempStorage(i,22) = neighbourOctal%u_i_minus_1(neighbourSubcell)
          tempStorage(i,23) = neighbourOctal%pressure_i_minus_1(neighbourSubcell)

          tempStorage(i,24) = qViscosity(1,1)
          tempStorage(i,25) = qViscosity(1,2)
          tempStorage(i,26) = qViscosity(1,3)
          tempStorage(i,27) = qViscosity(2,1)
          tempStorage(i,28) = qViscosity(2,2)
          tempStorage(i,29) = qViscosity(2,3)
          tempStorage(i,30) = qViscosity(3,1)
          tempStorage(i,31) = qViscosity(3,2)
          tempStorage(i,32) = qViscosity(3,3)

          tempStorage(i,58) = dble(neighbourOctal%nChildren)

       enddo

       allocate(temp1d(nStorage*nLoc))
       temp1d = RESHAPE(tempStorage, SHAPE=SHAPE(temp1d))

       call MPI_SEND(tempStorage, nStorage*nLoc, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, ierr)
       deallocate(temp1d, tempStorage)
    endif
  end subroutine receiveAcrossMpiBoundaryLevel
  
  subroutine exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, nDepth)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: iPair, nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup, iGroup
    integer :: rBound
    integer :: nDepth
    character(len=10) :: boundaryType(6) = (/"left  ","right ", "top   ", "bottom", "front ", "back  "/)
    integer :: ierr
    
    CALL MPI_BARRIER(amrCOMMUNICATOR, ierr)

    do iGroup = 1, nGroup

       do iPair = 1, nPairs

          if (group(iPair) == iGroup) then
             !       write(*,*) "myrank ", iPair, thread1(iPair), thread2(iPair)
             if ((myRankGlobal == thread1(iPair)).or.(myRankGlobal == thread2(iPair))) then
                !          write(*,*) myrank, " calling receive across boundary",iPair,nbound(ipair),boundaryType(nBound(iPair))
                call receiveAcrossMpiBoundaryLevel(grid, boundaryType(nBound(iPair)), thread1(iPair), thread2(iPair), nDepth)
                !          write(*,*) myrank, " finished receive across boundary"
                if      (nBound(iPair) == 1) then
                   rBound = 2
                else if (nBound(iPair) == 2) then
                   rBound = 1
                else if (nBound(iPair) == 3) then
                   rBound = 4
                else if (nBound(iPair) == 4) then
                   rBound = 3
                else if (nBound(iPair) == 5) then
                   rBound = 6
                else if (nBound(iPair) == 6) then
                   rBound = 5
                endif
                call receiveAcrossMpiBoundaryLevel(grid, boundaryType(rBound), thread2(iPair), thread1(iPair), nDepth)
                !          if (myRank == 1) then
                !             write(*,*) "Exchange done for direction ",nbound(ipair), " and ",rBound
                !             write(*,*) "Between threads ", thread1(ipair), " and ", thread2(ipair)
                !          endif
             endif
          endif
       enddo
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    enddo

  end subroutine exchangeAcrossMPIboundaryLevel

!Calculate the number of neighbouring cells on different MPI threads                                                       
 recursive subroutine determineCornerPairs(thisOctal, grid, nCornerPairs, cornerThread1, cornerThread2, nCornerBound, &    
       CornerGroup, nCornerGroup, iThread)
!      use mpi                                                                                                             
      use vector_mod, only: modulus
      type(GRIDTYPE) :: grid
      integer, intent(out) :: ncornerPairs, cornerthread1(:), cornerthread2(:)                                             
      integer, intent(out) :: ncornerBound(:), ncornerGroup, cornergroup(:)                                                
      type(octal), pointer :: thisOctal
      type(octal), pointer :: neighbourOctal
      type(octal), pointer :: child
      integer :: subcell, neighbourSubcell
      integer :: iThread
      integer :: numMPINeighbours, nDir
      integer :: i, k, i1, i2, h, j
      type(VECTOR) :: dirVec(6), outsider(6)
      type(VECTOR) :: locator, v1, v2, cVec
      integer :: iCornerPair
      logical :: alreadyInList
      logical, save :: firstTime=.true.

!    call MPI_COMM_RANK(localWorldCommunicator, myRank, ierr)                                                                      

      if(firstTime) then
         nCornerPairs = 0
         firsttime = .false.
      end if

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
! find the child                                                                                                           
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call determineCornerPairs(child, grid, nCornerPairs, cornerThread1, cornerThread2, nCornerBound, &         
                   cornerGroup, nCornerGroup, iThread)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, iThread)) cycle

          if(.not. thisOctal%ghostcell(subcell)) then

             if (thisOctal%threed) then
                nDir = 6
                dirVec(1) = VECTOR( -1.d0, 0.d0,  0.d0)
                dirVec(2) = VECTOR( +1.d0, 0.d0,  0.d0)
                dirVec(3) = VECTOR( 0.d0,  0.d0,  +1.d0)
                dirVec(4) = VECTOR( 0.d0,  0.d0,  -1.d0)
                dirVec(5) = VECTOR( 0.d0, +1.d0,  0.d0)
                dirVec(6) = VECTOR( 0.d0, -1.d0, 0.d0)
             else if (thisOctal%twod) then
                nDir = 4
                dirVec(1) = VECTOR( -1.d0, 0.d0, 0.d0)
                dirVec(2) = VECTOR( +1.d0,0.d0, 0.d0)
                dirVec(3) = VECTOR( 0.d0, 0.d0,  +1.d0)
                dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
             else
                nDir = 2
                dirVec(1) = VECTOR( -1.d0, 0.d0, 0.d0)
                dirVec(2) = VECTOR( +1.d0, 0.d0, 0.d0)
             endif

             j = 1
             outsider = VECTOR(0.d0, 0.d0, 0.d0)
             v1 = VECTOR(0.d0, 0.d0, 0.d0)
             v2 = VECTOR(0.d0, 0.d0, 0.d0)
             cVec = VECTOR(0.d0, 0.d0, 0.d0)
             numMPIneighbours = 0
             do k = 1, nDir
                locator = subcellCentre(thisOctal, subcell) + &
                 dirVec(k)*(thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
                neighbourOctal => thisOctal
                if(inOctal(grid%octreeRoot, locator)) then
                   call findSubcellLocal(locator, neighbourOctal, neighbourSubcell) 
                   if(.not. octalOnThread(neighbourOctal, neighbourSubcell, iThread)) then
                      numMPIneighbours = numMPIneighbours + 1
                      outsider(j)= dirVec(k)
                      j = j + 1
                      if(numMPIneighbours > 1) then !lies on a corner
                         !find the last two vectors                  
                         do h = 6, 1, -1
                            if(modulus(outsider(h)) /= 0.d0) then
                               if(modulus(v1) == 0.d0) then
                                  v1 = outsider(h)
                               else
                                  v2 = outsider(h)
                                  exit
                               end if
                            end if
                         end do
!the new corner lies between the last two vectors
                         cVec = v1 + v2
                         locator = subcellCentre(thisOctal, subcell) + &
                            cVec*((sqrt(2.d0)*(thisoctal%subcellsize/2.d0))+0.01d0*grid%halfsmallestsubcell) 

                           if(inOctal(grid%octreeRoot, locator)) then
                               call findSubcellLocal(locator, neighbourOctal, neighbourSubcell)
                               if(.not. octalOnThread(neighbourOctal, neighbourSubcell, iThread)) then
                                  i1 = iThread
                                  i2 = neighbourOctal%mpiThread(neighbourSubcell)
                                  alreadyInList = .false.
                                  do iCornerPair = 1, nCornerPairs
                                     if(((i1 == cornerThread1(iCornerPair)).and.(i2 == cornerThread2(iCornerPair)).or. &    
                                        ((i2 == cornerThread1(iCornerPair)).and.(i1 == cornerThread2(iCornerPair))))) then  
                                       alreadyInList = .true.
                                     end if
                                  end do

                                  if(.not. alreadyInList) then
                                     nCornerPairs = nCornerPairs + 1
                                     cornerThread1(nCornerPairs) = i1
                                     cornerThread2(nCornerPairs) = i2
                                     nCornerBound(nCornerPairs) = getNCornerBoundFromDirection(cVec, 0)
                                  end if
                               else
                                  call torus_abort("Expected cell on a different thread not found")
                               end if 
                            end if   
                         end if  
                      end if   
                   end if      
                end do         
             end if            
          end if              
       end do             
   end subroutine determineCornerPairs

   subroutine returnCornerPairs(grid, nCornerPairs, cornerThread1, cornerThread2, nCornerBound, cornerGroup, nCornerGroup) 
      use mpi                    
      use utils_mod, only: indexx
      integer, allocatable :: indx(:), itmp(:) 
      integer :: list(1000), nList
      real, allocatable :: sort(:)
      type(GRIDTYPE) :: grid      
      integer :: i, iThread
      integer, intent(out) :: ncornerPairs, cornerthread1(:), cornerthread2(:) 
      integer, intent(out) :: ncornerBound(:), ncornerGroup, cornergroup(:)    
 
         nCornerPairs = 0                                 
         do iThread = 1, nHydroThreadsGlobal                       
            call determineCornerPairs(grid%octreeRoot, grid, nCornerPairs, cornerthread1, cornerthread2, nCornerBound, &        
            cornerGroup, nCornerGroup, iThread) 
         enddo                            
 
      if(nCornerPairs > 1) then          
                                        
         allocate(indx(1:nCornerPairs), sort(1:nCornerPairs), itmp(1:nCornerPairs))
         do i = 1, nCornerPairs 
            sort(i) = real(cornerThread1(i))*100. + real(cornerThread2(i))  
         enddo                      
         call indexx(nCornerPairs, sort, indx) 
 
         do i = 1, nCornerPairs 
            itmp(i) = cornerThread1(indx(i)) 
         enddo                        
         cornerthread1(1:nCornerPairs) = itmp(1:nCornerPairs) 
         do i = 1, nCornerPairs                 
            itmp(i) = cornerThread2(indx(i))    
         enddo                                  
         cornerthread2(1:nCornerPairs) = itmp(1:nCornerPairs) 
         do i = 1, nCornerPairs                          
            itmp(i) = nCornerBound(indx(i))              
         enddo                                          
         nCornerBound(1:nCornerPairs) = itmp(1:nCornerPairs) 
         deallocate(indx, sort, itmp)                      
      endif          

      nCornerGroup = 1
      cornerGroup = 0 
      nList = 0      
      do while(any(cornerGroup(1:nCornerPairs)==0))
         do i = 1, nCornerPairs                    
            if (cornerGroup(i) == 0) then          
               if (.not.inList(cornerthread1(i), list, nList).and.& 
               (.not.inList(cornerthread2(i), list, nList))) then   

               cornerGroup(i) = nCornerGroup                        
               list(nList+1) = cornerthread1(i)                     
               list(nList+2) = cornerthread2(i)                     
               nList = nList + 2                                   
            endif
         endif                                          
      enddo                                                 
      nCornerGroup = nCornerGroup + 1                             
      nList = 0                                                   
      enddo                                                        
 !      print *, "done"                                           
   end subroutine returnCornerPairs          

  recursive subroutine determineBoundaryPairs(thisOctal, grid, nPairs,  thread1, thread2, nBound, iThread)

    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    !
    integer :: subcell, i, iThread
    type(VECTOR) :: dirVec(6), centre, octVec
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    integer :: thread1(:), thread2(:), nBound(:), nPairs, iPair
    integer :: i1, i2
    logical :: alreadyInList


    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call determineBoundaryPairs(child, grid, nPairs,  thread1, thread2, nBound, iThread)
                exit
             end if
          end do
       else

!          if ((grid%splitOverMpi).and.(thisOctal%mpiThread(subcell) /= ithread)) then
!             cycle
!          endif
          if (.not.octalOnThread(thisOctal, subcell, iThread)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell)
          if (thisOctal%threed) then
             nDir = 6
             dirVec(1) = VECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(2) = VECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
             dirVec(5) = VECTOR( 0.d0, 1.d0,  0.d0)
             dirVec(6) = VECTOR( 0.d0,-1.d0,  0.d0)
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR(-1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(+1.d0, 0.d0, 0.d0)
             dirVec(3) = VECTOR(0.d0, 0.d0, +1.d0)
             dirVec(4) = VECTOR(0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR(-1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(+1.d0, 0.d0, 0.d0)
          endif
          do j = 1, nDir
             octVec = centre + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then
                neighbourOctal => thisOctal 
               call findSubcellTD(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell)
!                if (neighbourOctal%mpiThread(neighboursubcell) /= iThread) then
                   if (.not.octalOnThread(neighbourOctal, neighbourSubcell, iThread)) then
                   i1 = ithread
                   i2 = neighbourOctal%mpiThread(neighboursubcell)
                   alreadyInList = .false.
                   do iPair = 1, nPairs
                      if (((i1 == thread1(iPair)).and.(i2 == thread2(iPair)).or. &
                           ((i2 == thread1(iPair)).and.(i1 == thread2(iPair))))) then
                         alreadyinList = .true.
                      endif
                   enddo
                   if (.not.alreadyInList) then
                      nPairs = nPairs + 1
                      thread1(nPairs) = i1
                      thread2(nPairs) = i2
                      nBound(nPairs) = j
                   endif
                endif
                
             end if
          enddo
       endif
    end do

  end subroutine determineBoundaryPairs

  subroutine autoSetupEvenupArray(grid, evenUpArray)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: evenUpArray(:)
    integer, allocatable :: nNeighbours(:), neighbourList(:,:), tempInt(:)
    integer :: ierr, iThread, iPass, iThread1, iThread2, j
    logical :: checkNoSharedNeighbours
    allocate(nNeighbours(1:nHydroThreadsGlobal), neighbourList(1:nHydroThreadsGlobal,1:26))
    
    nNeighbours = 0
    neighbourList = 0
    call determineNeighbourList(grid%octreeRoot, grid, nNeighbours(myrankGlobal), neighbourList(myrankGlobal,:))
    allocate(tempInt(1:nHydrothreadsGlobal))
    call MPI_ALLREDUCE(nNeighbours, tempInt, nHydrothreadsGlobal, MPI_INTEGER, MPI_SUM, amrCommunicator, ierr)
    nNeighbours = tempInt
    deallocate(tempInt)

    allocate(tempInt(26))
    do ithread = 1, nHydrothreadsGlobal
       call MPI_ALLREDUCE(neighbourList(iThread,:), tempInt, 26, MPI_INTEGER, MPI_SUM, amrCommunicator, ierr)
       neighbourList(iThread,1:26) = tempInt(1:26)
    enddo
    deallocate(tempInt)

    if (writeoutput) then
       do ithread = 1, nHydrothreadsGlobal
          write(*,'(26i4)') ithread, nNeighbours(ithread), neighbourList(ithread,1:nNeighbours(ithread))
       enddo
    endif

    iPass = 1
    evenUpArray = 0


    do iThread1 = 1, nHydrothreadsGlobal
       if (evenUpArray(iThread1) == 0) then
          evenUpArray(iThread1) = iPass
       else
          cycle
       endif
       do iThread2 = iThread1, nHydrothreadsGlobal
          if (evenUpArray(iThread2) == 0) then

             checkNoSharedNeighbours = .true.
             do j = 1, nHydroThreadsGlobal
                if (evenupArray(j) == iPass) then
                   checkNoSharedNeighbours = checkNosharedNeighbours.and.&
                        noSharedNeighbours(j, iThread2, nHydroThreadsGlobal, nNeighbours, neighbourList)
                endif
             enddo
             if (checkNoSharedNeighbours) then
                evenUpArray(iThread2) = iPass
             endif
          endif
       enddo
       iPass = iPass + 1
    enddo
    if (Writeoutput) write(*,*) " evenupArray auto ",evenUpArray(1:nHydroThreadsGlobal)
  end subroutine autoSetupEvenupArray

  logical function noSharedNeighbours(i1, i2, nThreads, nNeighbours, neighbourList)
    integer, intent(in) :: nThreads
    integer :: i1, i2, nNeighbours(nThreads), neighbourList(nThreads,26), i

    noSharedNeighbours = .true.

    do i = 1, nNeighbours(i1)
       if (ANY(neighbourList(i2,1:nNeighbours(i2)) == neighbourList(i1,i))) then
          noSharedNeighbours = .false.
          exit
       endif
    enddo
    if (ANY(neighbourList(i1,1:nNeighbours(i1)) == i2).or. &
         ANY(neighbourList(i2,1:nNeighbours(i2)) == i1)) noSharedNeighbours = .false.

  end function noSharedNeighbours




  recursive subroutine determineNeighbourList(thisOctal, grid, nNeighbours, neighbourList)

    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    !
    integer :: subcell, i
    type(VECTOR) :: dirVec(26), centre, octVec
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    integer :: nNeighbours, neighbourList(:)


    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call determineNeighbourList(child, grid, nNeighbours, neighbourList)
                exit
             end if
          end do
       else

          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle

          r = thisOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell)
          if (thisOctal%threed) then
             nDir = 26
             dirVec(1) = VECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(2) = VECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
             dirVec(5) = VECTOR( 0.d0, 1.d0,  0.d0)
             dirVec(6) = VECTOR( 0.d0,-1.d0,  0.d0)

             dirVec(7)   = VECTOR( RootThree,  RootThree, -RootThree)
             dirVec(8)   = VECTOR( RootThree, -RootThree, -RootThree)
             dirVec(9)   = VECTOR( RootThree,  RootThree,  RootThree)
             dirVec(10)  = VECTOR( RootThree, -RootThree,  RootThree)
             dirVec(11)  = VECTOR(-RootThree,  RootThree, -RootThree)
             dirVec(12)  = VECTOR(-RootThree, -RootThree, -RootThree)
             dirVec(13)  = VECTOR(-RootThree,  RootThree,  RootThree)
             dirVec(14)  = VECTOR(-RootThree, -RootThree,  RootThree)

             dirVec(15) = VECTOR( RootTwo,  0.d0, RootTwo)
             dirVec(16) = VECTOR(-RootTwo,  0.d0,-RootTwo)
             dirVec(17) = VECTOR( RootTwo,  0.d0,-RootTwo)
             dirVec(18) = VECTOR(-RootTwo,  0.d0, RootTwo)

             dirVec(19) = VECTOR(  0.d0, RootTwo, RootTwo)
             dirVec(20) = VECTOR(  0.d0,-RootTwo,-RootTwo)
             dirVec(21) = VECTOR(  0.d0, RootTwo,-RootTwo)
             dirVec(22) = VECTOR(  0.d0,-RootTwo, RootTwo)

             dirVec(23) = VECTOR( RootTwo,  RootTwo,  0.d0)
             dirVec(24) = VECTOR(-RootTwo, -RootTwo,  0.d0)
             dirVec(25) = VECTOR( RootTwo, -RootTwo,  0.d0)
             dirVec(26) = VECTOR(-RootTwo,  RootTwo,  0.d0)

          else if (thisOctal%twod) then

             nDir = 8
             dirVec(1) = VECTOR(-1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(+1.d0, 0.d0, 0.d0)
             dirVec(3) = VECTOR(0.d0, 0.d0, +1.d0)
             dirVec(4) = VECTOR(0.d0, 0.d0, -1.d0)

             dirVec(5) = VECTOR( rootTwo, 0.d0, rootTwo)
             dirVec(6) = VECTOR(-rootTwo, 0.d0, rootTwo)
             dirVec(7) = VECTOR( rootTwo, 0.d0,-rootTwo)
             dirVec(8) = VECTOR(-rootTwo, 0.d0,-rootTwo)

          else
             nDir = 2
             dirVec(1) = VECTOR(-1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(+1.d0, 0.d0, 0.d0)
          endif
          do j = 1, nDir
             octVec = centre + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then
                neighbourOctal => thisOctal 
                call findSubcellTD(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell)
                if (.not.octalOnThread(neighbourOctal, neighbourSubcell, myrankGlobal)) then
                   if (nNeighbours == 0) then
                      nNeighbours = nNeighbours + 1
                      neighbourList(nNeighbours) = neighbourOctal%mpiThread(neighbourSubcell)
                   else if (.not.(ANY(neighbourList(1:nNeighbours) == neighbourOctal%mpiThread(neighbourSubcell)))) then
                      nNeighbours = nNeighbours + 1
                      neighbourList(nNeighbours) = neighbourOctal%mpiThread(neighbourSubcell)
                   endif
                endif
             end if
          enddo
       endif
    end do

  end subroutine determineNeighbourList

  subroutine returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    use utils_mod, only: indexx
    use mpi
    type(GRIDTYPE) :: grid
    integer, intent(out) :: nPairs, thread1(:), thread2(:), nBound(:), nGroup, group(:)
    integer :: iThread
    integer :: i
    integer, allocatable :: indx(:), itmp(:)
    real, allocatable :: sort(:)
    integer :: list(10000), nList


    nPairs = 0
    do iThread = 1, nHydroThreadsGlobal
       call determineBoundaryPairs(grid%octreeRoot, grid, nPairs,  thread1, thread2, nBound, iThread)
    enddo

    if (nPairs > 1) then
       allocate(indx(1:nPairs), sort(1:nPairs), itmp(1:nPairs))
       do i = 1, nPairs
          sort(i) = real(thread1(i))*512. + real(thread2(i))
       enddo
       call indexx(nPairs, sort, indx)
       
       
       do i = 1, nPairs
          itmp(i) = thread1(indx(i))
       enddo
       thread1(1:nPairs) = itmp(1:nPairs)
       do i = 1, nPairs
          itmp(i) = thread2(indx(i))
       enddo
       thread2(1:nPairs) = itmp(1:nPairs)
       do i = 1, nPairs
          itmp(i) = nBound(indx(i))
       enddo
       nBound(1:nPairs) = itmp(1:nPairs)
       deallocate(indx, sort, itmp)
    endif
    

    nGroup = 1
    group = 0
    nList = 0
    do while(any(group(1:nPairs)==0))
       do i = 1, nPairs
          if (group(i) == 0) then
             if (.not.inList(thread1(i), list, nList).and.(.not.inList(thread2(i), list, nList))) then
                group(i) = nGroup
                list(nList+1) = thread1(i)
                list(nList+2) = thread2(i)
                nList = nList + 2
             endif
          endif
       enddo
       nGroup = nGroup + 1
       nList = 0
    enddo
  end subroutine returnBoundaryPairs

  logical function inList(num, list, nList)
    integer :: num
    integer :: list(:)
    integer :: nList
    integer :: i 
    inList = .false.
    do i = 1, nList
       if (num == list(i)) inList = .true.
    enddo
  end function inList

  subroutine getNeighbourValues(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
       rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xplus, px, py, pz, rm1, rum1, pm1, qViscosity)
    use mpi
!    use inputs_mod, only : useTensorViscosity

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal, tOctal
    type(VECTOR) :: direction, rVec, pVec
    integer :: subcell, neighbourSubcell, tSubcell
    integer :: nBound, nDepth
    integer, intent(out) :: nd, nc

    real(double), intent(out) :: q, rho, rhoe, rhou, rhov, rhow, qnext, x, pressure, flux, phi, phigas
    real(double), intent(out) :: xplus, px, py, pz, qViscosity(3,3), rm1, rum1, pm1
    real(double) :: temperature

    nbound = getNboundFromDirection(direction)

    if ((thisOctal%twoD).and.((nBound == 5).or. (nBound == 6))) then
       write(*,*) "Bonndary error for twod: ",nbound
       x = -2.d0
       x = sqrt(x)
    endif
    if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then ! .or. &
!         thisOctal%mpiThread(subcell) == neighbourOctal%mpiThread(neighboursubcell)) then

       nd = neighbourOctal%nDepth
       nc = neighbourOctal%nChildren

       x = neighbourOctal%x_i(neighbourSubcell)
       pVEc = subcellCentre(neighbourOctal, neighbourSubcell)

       px = pVec%x
       py = pVec%y
       pz = pVec%z

       if (thisOctal%nDepth == neighbourOctal%nDepth) then ! same level

          q   = neighbourOctal%q_i(neighbourSubcell)
          rho = neighbourOctal%rho(neighbourSubcell)
          rhoe = neighbourOctal%rhoe(neighbourSubcell)
          rhou = neighbourOctal%rhou(neighbourSubcell)
          rhov = neighbourOctal%rhov(neighbourSubcell)
          rhow = neighbourOctal%rhow(neighbourSubcell)
          pressure = neighbourOctal%pressure_i(neighbourSubcell)
          flux = neighbourOctal%flux_i(neighbourSubcell)
          phi = neighbourOctal%phi_i(neighbourSubcell)
          phigas = neighbourOctal%phi_gas(neighbourSubcell)
          xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
          qViscosity = neighbourOctal%qViscosity(neighbourSubcell,:,:)
          rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
          rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
          pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)

       else if (thisOctal%nDepth > neighbourOctal%nDepth) then ! fine cells set to coarse cell fluxes (should be interpolated here!!!)
          q   = neighbourOctal%q_i(neighbourSubcell)
          rho = neighbourOctal%rho(neighbourSubcell)
          rhoe = neighbourOctal%rhoe(neighbourSubcell)
          rhou = neighbourOctal%rhou(neighbourSubcell)
          rhov = neighbourOctal%rhov(neighbourSubcell)
          rhow = neighbourOctal%rhow(neighbourSubcell)
          pressure = neighbourOctal%pressure_i(neighbourSubcell)
          phi = neighbourOctal%phi_i(neighbourSubcell)
          phigas = neighbourOctal%phi_gas(neighbourSubcell)
          xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
          qViscosity = neighbourOctal%qViscosity(neighbourSubcell,:,:)

          rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
          rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
          pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)

          if (thisOctal%oneD) then
             flux = neighbourOctal%flux_i(neighbourSubcell)
          else if (thisOctal%twoD) then
             flux = neighbourOctal%flux_i(neighbourSubcell)!/2.d0
!             flux = neighbourOctal%flux_i(neighbourSubcell)/2.d0
          else
             flux = neighbourOctal%flux_i(neighbourSubcell)!/4.d0
          endif
       else
          call averageValue(direction, neighbourOctal,  neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, &
               flux, phi, phigas, rm1, rum1, pm1, qViscosity, temperature) ! fine to coarse
          xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
          
       endif

       
       rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
            direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
       if (inOctal(grid%octreeRoot, rVec)) then
          tOctal => neighbourOctal
          tSubcell = neighbourSubcell
          call findSubcellLocal(rVec, tOctal, tSubcell)
          
          if (tOctal%mpiThread(tSubcell) == myRankGlobal) then
             qnext = tOctal%q_i(tSubcell)
          else
             if (associated(neighbourOctal%mpiBoundaryStorage)) then
                qNext = neighbourOctal%mpiBoundaryStorage(neighbourSubcell, nBound, 1)
             else
                qNext = 0.d0
             endif
          endif
       else
          qNext = 0.d0
       endif



    else


       if (.not.associated(thisOctal%mpiBoundaryStorage)) then
          write(*,*) "boundary storage not allocated when it should be!", myrankGlobal, &
               neighbourOctal%mpiThread(neighboursubcell), &
               thisOctal%mpiThread(subcell)
          write(*,*) "direction",  direction,nBound
          write(*,*) "depth ",thisOctal%nDepth
          write(*,*) "nChildren ",thisOctal%nChildren
          write(*,*) "rho ",thisOctal%rho(subcell)
          write(*,*) "this centre",subcellCentre(thisOctal, subcell)
          write(*,*) "neig centre",subcellCentre(neighbourOctal, neighboursubcell)
          x = -2.d0
          x = sqrt(x)
       endif

       x = thisOctal%mpiBoundaryStorage(subcell, nBound, 7)
       xplus = thisOctal%mpiBoundaryStorage(subcell, nBound, 14)

       nd =  nint(thisOctal%mpiBoundaryStorage(subcell, nBound,9))

       nDepth = nint(thisOctal%mpiBoundaryStorage(subcell, nBound,9))

       q   = thisOctal%mpiBoundaryStorage(subcell, nBound, 1)
       rho = thisOctal%mpiBoundaryStorage(subcell, nBound, 2)
       rhoe = thisOctal%mpiBoundaryStorage(subcell, nBound, 3)
       rhou = thisOctal%mpiBoundaryStorage(subcell, nBound, 4)
       rhov = thisOctal%mpiBoundaryStorage(subcell, nBound, 5)
       rhow = thisOctal%mpiBoundaryStorage(subcell, nBound, 6)
       qnext = thisOctal%mpiBoundaryStorage(subcell, nBound, 8)
       pressure = thisOctal%mpiBoundaryStorage(subcell, nBound, 10)
       flux =  thisOctal%mpiBoundaryStorage(subcell, nBound, 11)
       phi = thisOctal%mpiBoundaryStorage(subcell, nBound, 12)
       phigas = thisOctal%mpiBoundaryStorage(subcell, nBound, 13)
       px = thisOctal%mpiBoundaryStorage(subcell, nBound, 15)
       py = thisOctal%mpiBoundaryStorage(subcell, nBound, 16)
       pz = thisOctal%mpiBoundaryStorage(subcell, nBound, 17)
       rm1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 21)
       rum1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 22)
       pm1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 23)


       qViscosity(1,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 24)
       qViscosity(1,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 25)
       qViscosity(1,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 26)
       qViscosity(2,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 27)
       qViscosity(2,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 28)
       qViscosity(2,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 29)
       qViscosity(3,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 30)
       qViscosity(3,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 31)
       qViscosity(3,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 32)

       nc = nint(thisOctal%mpiBoundaryStorage(subcell, nBound, 58))

    endif
  end subroutine getNeighbourValues


  subroutine getNeighbourValues2(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
       rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xplus, px, py, pz, rm1, rum1, pm1, qViscosity)
    use mpi
    use inputs_mod, only : smallestCellSize

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal, tOctal
    type(VECTOR) :: direction, rVec, pVec, locator
    integer :: subcell, neighbourSubcell, tSubcell
    integer :: nBound, nDepth, j
    integer, intent(out) :: nd

    real(double), intent(out) :: q(2), rho(2), rhoe(2), rhou(2), rhov(2), rhow(2), qnext, x, pressure(2), flux(2), phi, phigas
    real(double), intent(out) :: xplus, px, py, pz, qViscosity(3,3), rm1, rum1, pm1

    nbound = getNboundFromDirection(direction)

    if ((thisOctal%twoD).and.((nBound == 5).or. (nBound == 6))) then
       write(*,*) "Bonndary error for twod: ",nbound
       x = -2.d0
       x = sqrt(x)
    endif
    if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then

       nd = neighbourOctal%nDepth

       x = neighbourOctal%x_i(neighbourSubcell)
       pVEc = subcellCentre(neighbourOctal, neighbourSubcell)

       px = pVec%x
       py = pVec%y
       pz = pVec%z

       if (thisOctal%nDepth == neighbourOctal%nDepth) then ! same level

          q(1:2)   = neighbourOctal%q_i(neighbourSubcell)
          rho(1:2) = neighbourOctal%rho(neighbourSubcell)
          rhoe(1:2) = neighbourOctal%rhoe(neighbourSubcell)
          rhou(1:2) = neighbourOctal%rhou(neighbourSubcell)
          rhov(1:2) = neighbourOctal%rhov(neighbourSubcell)
          rhow(1:2) = neighbourOctal%rhow(neighbourSubcell)
          pressure(1:2) = neighbourOctal%pressure_i(neighbourSubcell)
          flux(1:2) = neighbourOctal%flux_amr_i(neighbourSubcell,1:2)
          phi = neighbourOctal%phi_i(neighbourSubcell)
          phigas = neighbourOctal%phi_gas(neighbourSubcell)
          xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
          qViscosity = neighbourOctal%qViscosity(neighbourSubcell,:,:)

          rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
          rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
          pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)

       else if (thisOctal%nDepth > neighbourOctal%nDepth) then ! fine cells set to coarse cell fluxes (should be interpolated here!!!)
          q(1:2)   = neighbourOctal%q_i(neighbourSubcell)
          rho(1:2) = neighbourOctal%rho(neighbourSubcell)
          rhoe(1:2) = neighbourOctal%rhoe(neighbourSubcell)
          rhou(1:2) = neighbourOctal%rhou(neighbourSubcell)
          rhov(1:2) = neighbourOctal%rhov(neighbourSubcell)
          rhow(1:2) = neighbourOctal%rhow(neighbourSubcell)
          pressure(1:2) = neighbourOctal%pressure_i(neighbourSubcell)
          phi = neighbourOctal%phi_i(neighbourSubcell)
          phigas = neighbourOctal%phi_gas(neighbourSubcell)
          xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
          qViscosity = neighbourOctal%qViscosity(neighbourSubcell,:,:)

          rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
          rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
          pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)

          if (abs(direction%x) > 0.9d0) then ! radial
             if ((subcell == 1).or.(subcell == 2)) then
                flux(1:2) = neighbourOctal%flux_amr_i(neighbourSubcell,1)
             else
                flux(1:2) = neighbourOctal%flux_amr_i(neighbourSubcell,2)
             endif
          else ! vertical
             if ((subcell == 1).or.(subcell == 3)) then
                flux(1:2) = neighbourOctal%flux_amr_i(neighbourSubcell,1)
             else
                flux(1:2) = neighbourOctal%flux_amr_i(neighbourSubcell,2)
             endif
          endif

       else

          do j = 1, 2

             rVec = subcellCentre(thisOctal, subcell)
             if (abs(direction%x) > 0.9d0) then
                
                if (j == 1) then
                   locator = VECTOR(rVec%x + (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction%x, 0.d0, &
                        rVec%z - 0.25d0 * thisOctal%subcellSize)
                else
                   locator = VECTOR(rVec%x + (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction%x, 0.d0, &
                        rVec%z + 0.25d0 * thisOctal%subcellSize)
                endif
             else ! z direction
                if (j == 1) then
                   locator = VECTOR(rVec%x - 0.25d0 * thisOctal%subcellSize, 0.d0, &
                        rVec%z + (thisOctal%subcellSize/2.d0 + smallestCellSize * 0.1d0)*direction%z)
                else
                   locator = VECTOR(rVec%x + 0.25d0 * thisOctal%subcellSize, 0.d0, &
                        rVec%z + (thisOctal%subcellSize/2.d0 + smallestCellSize * 0.1d0)*direction%z)
                endif
             endif
             tOctal => thisOctal
             call findSubcellLocal(locator, tOctal, tSubcell)

             q(j)   = tOctal%q_i(tSubcell)
             rho(j) = tOctal%rho(tSubcell)
             rhoe(j) = tOctal%rhoe(tSubcell)
             rhou(j) = tOctal%rhou(tSubcell)
             rhov(j) = tOctal%rhov(tSubcell)
             rhow(j) = tOctal%rhow(tSubcell)
             pressure(j) = tOctal%pressure_i(tSubcell)
             phi = neighbourOctal%phi_i(neighbourSubcell)
             phigas = neighbourOctal%phi_gas(neighbourSubcell)
             xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
             qViscosity = neighbourOctal%qViscosity(neighbourSubcell,:,:)

             rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
             rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
             pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)

             flux(j) = tOctal%flux_amr_i(tSubcell,j)
          enddo
       endif

       
       rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
            direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
       if (inOctal(grid%octreeRoot, rVec)) then
          tOctal => neighbourOctal
          tSubcell = neighbourSubcell
          call findSubcellLocal(rVec, tOctal, tSubcell)
          
          if (tOctal%mpiThread(tSubcell) == myRankGlobal) then
             qnext = tOctal%q_i(tSubcell)
          else
             if (associated(neighbourOctal%mpiBoundaryStorage)) then
                qNext = neighbourOctal%mpiBoundaryStorage(neighbourSubcell, nBound, 1)
             else
                qNext = 0.d0
             endif
          endif
       else
          qNext = 0.d0
       endif



    else


       if (.not.associated(thisOctal%mpiBoundaryStorage)) then
          write(*,*) "boundary storage not allocated when it should be!", myrankGlobal, &
               neighbourOctal%mpiThread(neighboursubcell), &
               thisOctal%mpiThread(subcell)
          write(*,*) "direction",  direction,nBound
          write(*,*) "depth ",thisOctal%nDepth
          write(*,*) "nChildren ",thisOctal%nChildren
          write(*,*) "rho ",thisOctal%rho(subcell)
          write(*,*) "this centre",subcellCentre(thisOctal, subcell)
          write(*,*) "neig centre",subcellCentre(neighbourOctal, neighboursubcell)
          x = -2.d0
          x = sqrt(x)
       endif

       x = thisOctal%mpiBoundaryStorage(subcell, nBound, 7)
       xplus = thisOctal%mpiBoundaryStorage(subcell, nBound, 14)

       nd =  nint(thisOctal%mpiBoundaryStorage(subcell, nBound,9))

       nDepth = nint(thisOctal%mpiBoundaryStorage(subcell, nBound,9))

       q   = thisOctal%mpiBoundaryStorage(subcell, nBound, 1)
       rho = thisOctal%mpiBoundaryStorage(subcell, nBound, 2)
       rhoe = thisOctal%mpiBoundaryStorage(subcell, nBound, 3)
       rhou = thisOctal%mpiBoundaryStorage(subcell, nBound, 4)
       rhov = thisOctal%mpiBoundaryStorage(subcell, nBound, 5)
       rhow = thisOctal%mpiBoundaryStorage(subcell, nBound, 6)
       qnext = thisOctal%mpiBoundaryStorage(subcell, nBound, 8)
       pressure = thisOctal%mpiBoundaryStorage(subcell, nBound, 10)
       flux =  thisOctal%mpiBoundaryStorage(subcell, nBound, 11)
       phi = thisOctal%mpiBoundaryStorage(subcell, nBound, 12)
       phigas = thisOctal%mpiBoundaryStorage(subcell, nBound, 13)
       px = thisOctal%mpiBoundaryStorage(subcell, nBound, 15)
       py = thisOctal%mpiBoundaryStorage(subcell, nBound, 16)
       pz = thisOctal%mpiBoundaryStorage(subcell, nBound, 17)
       rm1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 21)
       rum1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 22)
       pm1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 23)


       qViscosity(1,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 24)
       qViscosity(1,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 25)
       qViscosity(1,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 26)
       qViscosity(2,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 27)
       qViscosity(2,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 28)
       qViscosity(2,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 29)
       qViscosity(3,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 30)
       qViscosity(3,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 31)
       qViscosity(3,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 32)

       flux(1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 34)
       flux(2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 35)

    endif
  end subroutine getNeighbourValues2

! +ve x face                                                                                                                                                                                  
! -z -y 2  (1)                                                                                                                                                                                
! -z +y 4  (2)                                                                                                                                                                                
! +z -y 6  (3)                                                                                                                                                                                
! +z +y 8  (4)                                                                                                                                                                                

! -ve x face                                                                                                                                                                                  
! -z -y 1  (1)                                                                                                                                                                                
! -z +y 3  (2)                                                                                                                                                                                
! +z -y 5  (3)                                                                                                                                                                                
! +z +y 7  (4)                                                                                                                                                                                

! -ve z face                                                                                                                                                                                  
! -x -y 1  (1)                                                                                                                                                                                
! -x +y 3  (2)                                                                                                                                                                                
! +x -y 2  (3)                                                                                                                                                                                
! +x +y 4  (4)                                                                                                                                                                                

! +ve z face                                                                                                                                                                                  
! -x -y 5  (1)                                                                                                                                                                                
! -x +y 7  (2)                                                                                                                                                                                
! +x -y 6  (3)                                                                                                                                                                                
! +x +y 8  (4)                                                                                                                                                                                

! -ve y face                                                                                                                                                                                  
! -x -z 1  (1)                                                                                                                                                                                
! -x +z 5  (2)                                                                                                                                                                                
! +x -z 2  (3)                                                                                                                                                                                
! +x +z 6  (4)                                                                                                                                                                                

! +ve y face                                                                                                                                                                                  
! -x -z 3  (1)                                                                                                                                                                                
! -x +z 7  (2)                                                                                                                                                                                
! +x -z 4  (3)                                                                                                                                                                                
! +x +z 8  (4)                  
!       y                  z                                                                                                                                                                  
!       |                 /                                                                                                                                                                   
!       |      __________/______                                                                                                                                                              
!       |     /        /       /|                                                                                                                                                             
!       |    /   7    /   8   / |                                                                                                                                                             
!       |   /________/_______/  |                                                                                                                                                             
!       |  /        /       /| 8|   Diagram showing the convention used here for                                                                                                              
!       | /   3    /   4   / |  |   numbering the subcells of each octal.                                                                                                                     
!       |/________/_______/  |  |                                                                                                                                                             
!       |        |        | 4| /|                                                                                                                                                             
!       |        |        |  |/ |                                                                                                                                                             
!       |    3   |   4    |  /  |                                                                                                                                                             
!       |        |        | /| 6|                                                                                                                                                             
!       |________|________|/ |  /                                                                                                                                                             
!       |        |        |  | /                                                                                                                                                              
!       |        |        | 2|/                                                                                                                                                               
!       |    1   |   2    |  /                                                                                                                                                                
!       |        |        | /                                                                                                                                                                 
!       |________|________|/________\ x                                                                                                                                                       
!                                   /                                                                                                                                                         

  subroutine getNeighbourValues4(grid, thisOctal, subcell, neighbourOctal, neighbourSubcell, direction, q, rho, rhoe, &
       rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xplus, px, py, pz, rm1, rum1, pm1, qViscosity)
    use mpi
    use inputs_mod, only : smallestCellSize

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal, tOctal
    type(VECTOR) :: direction, rVec, pVec, locator
    type(VECTOR), parameter :: xhat = VECTOR(1.d0, 0.d0, 0.d0), yHat = VECTOR(0.d0, 1.d0, 0.d0), zHat = VECTOR(0.d0, 0.d0, 1.d0)
    integer :: subcell, neighbourSubcell, tSubcell
    integer :: nBound, nDepth, j
    integer, intent(out) :: nd

    real(double), intent(out) :: q(4), rho(4), rhoe(4), rhou(4), rhov(4), rhow(4), qnext, x, pressure(4), flux(4), phi, phigas
    real(double), intent(out) :: xplus, px, py, pz, qViscosity(3,3), rm1, rum1, pm1

    nbound = getNboundFromDirection(direction)

    if ((thisOctal%twoD).and.((nBound == 5).or. (nBound == 6))) then
       write(*,*) "Bonndary error for twod: ",nbound
       x = -2.d0
       x = sqrt(x)
    endif
    if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then

       nd = neighbourOctal%nDepth

       x = neighbourOctal%x_i(neighbourSubcell)
       pVEc = subcellCentre(neighbourOctal, neighbourSubcell)

       px = pVec%x
       py = pVec%y
       pz = pVec%z

       if (thisOctal%nDepth == neighbourOctal%nDepth) then ! same level

          q(1:4)   = neighbourOctal%q_i(neighbourSubcell)
          rho(1:4) = neighbourOctal%rho(neighbourSubcell)
          rhoe(1:4) = neighbourOctal%rhoe(neighbourSubcell)
          rhou(1:4) = neighbourOctal%rhou(neighbourSubcell)
          rhov(1:4) = neighbourOctal%rhov(neighbourSubcell)
          rhow(1:4) = neighbourOctal%rhow(neighbourSubcell)
          pressure(1:4) = neighbourOctal%pressure_i(neighbourSubcell)
          flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell,1:4)
          phi = neighbourOctal%phi_i(neighbourSubcell)
          phigas = neighbourOctal%phi_gas(neighbourSubcell)
          xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
          qViscosity = neighbourOctal%qViscosity(neighbourSubcell,:,:)

          rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
          rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
          pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)

       else if (thisOctal%nDepth > neighbourOctal%nDepth) then ! fine cells set to coarse cell fluxes (should be interpolated here!!!)
          q(1:4)   = neighbourOctal%q_i(neighbourSubcell)
          rho(1:4) = neighbourOctal%rho(neighbourSubcell)
          rhoe(1:4) = neighbourOctal%rhoe(neighbourSubcell)
          rhou(1:4) = neighbourOctal%rhou(neighbourSubcell)
          rhov(1:4) = neighbourOctal%rhov(neighbourSubcell)
          rhow(1:4) = neighbourOctal%rhow(neighbourSubcell)
          pressure(1:4) = neighbourOctal%pressure_i(neighbourSubcell)
          phi = neighbourOctal%phi_i(neighbourSubcell)
          phigas = neighbourOctal%phi_gas(neighbourSubcell)
          xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
          qViscosity = neighbourOctal%qViscosity(neighbourSubcell,:,:)

          rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
          rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
          pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)


!       y                  z                                                                                                                                                                  
!       |                 /                                                                                                                                                                   
!       |      __________/______                                                                                                                                                              
!       |     /        /       /|                                                                                                                                                             
!       |    /   7    /   8   / |                                                                                                                                                             
!       |   /________/_______/  |                                                                                                                                                             
!       |  /        /       /| 8|   Diagram showing the convention used here for                                                                                                              
!       | /   3    /   4   / |  |   numbering the subcells of each octal.                                                                                                                     
!       |/________/_______/  |  |                                                                                                                                                             
!       |        |        | 4| /|                                                                                                                                                             
!       |        |        |  |/ |                                                                                                                                                             
!       |    3   |   4    |  /  |                                                                                                                                                             
!       |        |        | /| 6|                                                                                                                                                             
!       |________|________|/ |  /                                                                                                                                                             
!       |        |        |  | /                                                                                                                                                              
!       |        |        | 2|/                                                                                                                                                               
!       |    1   |   2    |  /                                                                                                                                                                
!       |        |        | /                                                                                                                                                                 
!       |________|________|/________\ x                                                                                                                                                       
!                                   /                                                                                                                                                         


! +ve x face                                                                                                                                                                                  
! -z -y 2  (1)                                                                                                                                                                                
! -z +y 4  (2)                                                                                                                                                                                
! +z -y 6  (3)                                                                                                                                                                                
! +z +y 8  (4)                                                                                                                                                                                



          if (direction%x > 0.9d0) then
             select case(subcell)
                case(2) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 1)
                case(4) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 2)
                case(6) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 3)
                case(8) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 4)
                case DEFAULT
                   write(*,*) "error in logic x +ve ",subcell
             end select
          endif

! -ve x face                                                                                                                                                                                  
! -z -y 1  (1)                                                                                                                                                                                
! -z +y 3  (2)                                                                                                                                                                                
! +z -y 5  (3)                                                                                                                                                                                
! +z +y 7  (4)                                                                                                                                                                                

          if (direction%x < -0.9d0) then
             select case(subcell)
                case(1) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 1)
                case(3) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 2)
                case(5) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 3)
                case(7) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 4)
                case DEFAULT
                   write(*,*) "error in logic x -ve ",subcell
             end select
          endif

! -ve z face                                                                                                                                                                                  
! -x -y 1  (1)                                                                                                                                                                                
! -x +y 3  (2)                                                                                                                                                                                
! +x -y 2  (3)                                                                                                                                                                                
! +x +y 4  (4)                                                                                                                                                                                

          if (direction%z <  -0.9d0) then
             select case(subcell)
                case(1) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 1)
                case(3) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 2)
                case(2) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 3)
                case(4) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 4)
                case DEFAULT
                   write(*,*) "error in logic z -ve ",subcell
             end select
          endif

! +ve z face                                                                                                                                                                                  
! -x -y 5  (1)                                                                                                                                                                                
! -x +y 7  (2)                                                                                                                                                                                
! +x -y 6  (3)                                                                                                                                                                                
! +x +y 8  (4)                                                                                                                                                                                

          if (direction%z > 0.9d0) then
             select case(subcell)
                case(5) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 1)
                case(7) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 2)
                case(6) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 3)
                case(8) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 4)
                case DEFAULT
                   write(*,*) "error in logic z +ve ",subcell
             end select
          endif

! -ve y face                                                                                                                                                                                  
! -x -z 1  (1)                                                                                                                                                                                
! -x +z 5  (2)                                                                                                                                                                                
! +x -z 2  (3)                                                                                                                                                                                
! +x +z 6  (4)                                                                                                                                                                                


          if (direction%y < -0.9d0) then
             select case(subcell)
                case(1) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 1)
                case(5) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 2)
                case(2) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 3)
                case(6) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 4)
                case DEFAULT
                   write(*,*) "error in logic y -ve ",subcell
             end select
          endif

! +ve y face                                                                                                                                                                                  
! -x -z 3  (1)                                                                                                                                                                                
! -x +z 7  (2)                                                                                                                                                                                
! +x -z 4  (3)                                                                                                                                                                                
! +x +z 8  (4)                  

          if (direction%y > 0.9d0) then
             select case(subcell)
                case(3) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 1)
                case(7) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 2)
                case(4) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 3)
                case(8) 
                   flux(1:4) = neighbourOctal%flux_amr_i(neighbourSubcell, 4)
                case DEFAULT
                   write(*,*) "error in logic y +ve ",subcell
             end select
          endif


       else

          do j = 1, 4

             rVec = subcellCentre(thisOctal, subcell)

! +ve x face                                                                                                                                                                                  
! -z -y 2  (1)                                                                                                                                                                                
! -z +y 4  (2)                                                                                                                                                                                
! +z -y 6  (3)                                                                                                                                                                                
! +z +y 8  (4)                                                                                                                                                                                
             if (abs(direction%x) > 0.9d0) then
                select case(j)
                   case(1)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           - 0.25d0 * thisOctal%subcellSize * zHat &
                           - 0.25d0 * thisOctal%subcellSize * yHat
                   case(2)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           - (0.25d0 * thisOctal%subcellSize) * zHat &
                           + (0.25d0 * thisOctal%subcellSize) * yHat
                   case(3)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           + (0.25d0 * thisOctal%subcellSize) * zHat &
                           - (0.25d0 * thisOctal%subcellSize) * yHat
                   case(4)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           + (0.25d0 * thisOctal%subcellSize) * zHat &
                           + (0.25d0 * thisOctal%subcellSize) * yHat
                end select
                locator = locator + rvec
             endif
! -ve y face                                                                                                                                                                                  
! -x -z 1  (1)                                                                                                                                                                                
! -x +z 5  (2)                                                                                                                                                                                
! +x -z 2  (3)                                                                                                                                                                                
! +x +z 6  (4)                                                                                                                                                                                
             if (abs(direction%y) > 0.9d0) then
                select case(j)
                   case(1)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           - (0.25d0 * thisOctal%subcellSize) * xHat &
                           - (0.25d0 * thisOctal%subcellSize) * zHat
                   case(2)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           - (0.25d0 * thisOctal%subcellSize) * xHat &
                           + (0.25d0 * thisOctal%subcellSize) * zHat
                   case(3)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           + (0.25d0 * thisOctal%subcellSize) * xHat &
                           - (0.25d0 * thisOctal%subcellSize) * zHat
                   case(4)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           + (0.25d0 * thisOctal%subcellSize) * xHat &
                           + (0.25d0 * thisOctal%subcellSize) * zHat
                end select
                locator = locator + rvec
             endif
! -ve z face                                                                                                                                                                                  
! -x -y 1  (1)                                                                                                                                                                                
! -x +y 3  (2)                                                                                                                                                                                
! +x -y 2  (3)                                                                                                                                                                                
! +x +y 4  (4)                                                                                                                                                                                
             if (abs(direction%z) > 0.9d0) then
                select case(j)
                   case(1)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           - (0.25d0 * thisOctal%subcellSize) * xHat &
                           - (0.25d0 * thisOctal%subcellSize) * yHat
                   case(2)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           - (0.25d0 * thisOctal%subcellSize) * xHat &
                           + (0.25d0 * thisOctal%subcellSize) * yHat
                   case(3)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           + (0.25d0 * thisOctal%subcellSize) * xHat &
                           - (0.25d0 * thisOctal%subcellSize) * yHat
                   case(4)
                      locator = (thisOctal%subcellSize/2.d0 + smallestCellSize*0.1)*direction &
                           + (0.25d0 * thisOctal%subcellSize) * xHat &
                           + (0.25d0 * thisOctal%subcellSize) * yHat
                end select
                locator = locator + rvec
             endif




             tOctal => thisOctal
             call findSubcellLocal(locator, tOctal, tSubcell)

             q(j)   = tOctal%q_i(tSubcell)
             rho(j) = tOctal%rho(tSubcell)
             rhoe(j) = tOctal%rhoe(tSubcell)
             rhou(j) = tOctal%rhou(tSubcell)
             rhov(j) = tOctal%rhov(tSubcell)
             rhow(j) = tOctal%rhow(tSubcell)
             pressure(j) = tOctal%pressure_i(tSubcell)
             phi = neighbourOctal%phi_i(neighbourSubcell)
             phigas = neighbourOctal%phi_gas(neighbourSubcell)
             xplus = neighbourOctal%x_i_minus_1(neighbourSubcell)
             qViscosity = neighbourOctal%qViscosity(neighbourSubcell,:,:)

             rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
             rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
             pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)

             flux(j) = tOctal%flux_amr_i(tSubcell,j)
          enddo
       endif

       
       rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
            direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
       if (inOctal(grid%octreeRoot, rVec)) then
          tOctal => neighbourOctal
          tSubcell = neighbourSubcell
          call findSubcellLocal(rVec, tOctal, tSubcell)
          
          if (tOctal%mpiThread(tSubcell) == myRankGlobal) then
             qnext = tOctal%q_i(tSubcell)
          else
             if (associated(neighbourOctal%mpiBoundaryStorage)) then
                qNext = neighbourOctal%mpiBoundaryStorage(neighbourSubcell, nBound, 1)
             else
                qNext = 0.d0
             endif
          endif
       else
          qNext = 0.d0
       endif



    else


       if (.not.associated(thisOctal%mpiBoundaryStorage)) then
          write(*,*) "boundary storage not allocated when it should be!", myrankGlobal, &
               neighbourOctal%mpiThread(neighboursubcell), &
               thisOctal%mpiThread(subcell)
          write(*,*) "direction",  direction,nBound
          write(*,*) "depth ",thisOctal%nDepth
          write(*,*) "nChildren ",thisOctal%nChildren
          write(*,*) "rho ",thisOctal%rho(subcell)
          write(*,*) "this centre",subcellCentre(thisOctal, subcell)
          write(*,*) "neig centre",subcellCentre(neighbourOctal, neighboursubcell)
          x = -2.d0
          x = sqrt(x)
       endif

       x = thisOctal%mpiBoundaryStorage(subcell, nBound, 7)
       xplus = thisOctal%mpiBoundaryStorage(subcell, nBound, 14)

       nd =  nint(thisOctal%mpiBoundaryStorage(subcell, nBound,9))

       nDepth = nint(thisOctal%mpiBoundaryStorage(subcell, nBound,9))

       if (nDepth <= thisOctal%nDepth) then

          q   = thisOctal%mpiBoundaryStorage(subcell, nBound, 1)
          rho = thisOctal%mpiBoundaryStorage(subcell, nBound, 2)
          rhoe = thisOctal%mpiBoundaryStorage(subcell, nBound, 3)
          rhou = thisOctal%mpiBoundaryStorage(subcell, nBound, 4)
          rhov = thisOctal%mpiBoundaryStorage(subcell, nBound, 5)
          rhow = thisOctal%mpiBoundaryStorage(subcell, nBound, 6)
          qnext = thisOctal%mpiBoundaryStorage(subcell, nBound, 8)
          pressure = thisOctal%mpiBoundaryStorage(subcell, nBound, 10)
          flux =  thisOctal%mpiBoundaryStorage(subcell, nBound, 11)
          phi = thisOctal%mpiBoundaryStorage(subcell, nBound, 12)
          phigas = thisOctal%mpiBoundaryStorage(subcell, nBound, 13)
          px = thisOctal%mpiBoundaryStorage(subcell, nBound, 15)
          py = thisOctal%mpiBoundaryStorage(subcell, nBound, 16)
          pz = thisOctal%mpiBoundaryStorage(subcell, nBound, 17)
          rm1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 21)
          rum1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 22)
          pm1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 23)


          qViscosity(1,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 24)
          qViscosity(1,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 25)
          qViscosity(1,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 26)
          qViscosity(2,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 27)
          qViscosity(2,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 28)
          qViscosity(2,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 29)
          qViscosity(3,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 30)
          qViscosity(3,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 31)
          qViscosity(3,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 32)




!       y                  z                                                                                                                                                                  
!       |                 /                                                                                                                                                                   
!       |      __________/______                                                                                                                                                              
!       |     /        /       /|                                                                                                                                                             
!       |    /   7    /   8   / |                                                                                                                                                             
!       |   /________/_______/  |                                                                                                                                                             
!       |  /        /       /| 8|   Diagram showing the convention used here for                                                                                                              
!       | /   3    /   4   / |  |   numbering the subcells of each octal.                                                                                                                     
!       |/________/_______/  |  |                                                                                                                                                             
!       |        |        | 4| /|                                                                                                                                                             
!       |        |        |  |/ |                                                                                                                                                             
!       |    3   |   4    |  /  |                                                                                                                                                             
!       |        |        | /| 6|                                                                                                                                                             
!       |________|________|/ |  /                                                                                                                                                             
!       |        |        |  | /                                                                                                                                                              
!       |        |        | 2|/                                                                                                                                                               
!       |    1   |   2    |  /                                                                                                                                                                
!       |        |        | /                                                                                                                                                                 
!       |________|________|/________\ x                                                                                                                                                       
!                                   /                                                                                                                                                         


! +ve x face                                                                                                                                                                                  
! -z -y 2  (1)                                                                                                                                                                                
! -z +y 4  (2)                                                                                                                                                                                
! +z -y 6  (3)                                                                                                                                                                                
! +z +y 8  (4)                                                                                                                                                                                



          if (direction%x > 0.9d0) then
             select case(subcell)
                case(2) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 34)
                case(4) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 35)
                case(6) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 36)
                case(8) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 37)
                case DEFAULT
                   write(*,*) "error in logic x +ve ",subcell
             end select
          endif

! -ve x face                                                                                                                                                                                  
! -z -y 1  (1)                                                                                                                                                                                
! -z +y 3  (2)                                                                                                                                                                                
! +z -y 5  (3)                                                                                                                                                                                
! +z +y 7  (4)                                                                                                                                                                                

          if (direction%x < -0.9d0) then
             select case(subcell)
                case(1) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 34)
                case(3) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 35)
                case(5) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 36)
                case(7) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 37)
                case DEFAULT
                   write(*,*) "error in logic x -ve ",subcell
             end select
          endif

! -ve z face                                                                                                                                                                                  
! -x -y 1  (1)                                                                                                                                                                                
! -x +y 3  (2)                                                                                                                                                                                
! +x -y 2  (3)                                                                                                                                                                                
! +x +y 4  (4)                                                                                                                                                                                

          if (direction%z <  -0.9d0) then
             select case(subcell)
                case(1) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 34)
                case(3) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 35)
                case(2) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 36)
                case(4) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 37)
                case DEFAULT
                   write(*,*) "error in logic z -ve ",subcell
             end select
          endif

! +ve z face                                                                                                                                                                                  
! -x -y 5  (1)                                                                                                                                                                                
! -x +y 7  (2)                                                                                                                                                                                
! +x -y 6  (3)                                                                                                                                                                                
! +x +y 8  (4)                                                                                                                                                                                

          if (direction%z > 0.9d0) then
             select case(subcell)
                case(5) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 34)
                case(7) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 35)
                case(6) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 36)
                case(8) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 37)
                case DEFAULT
                   write(*,*) "error in logic z +ve ",subcell
             end select
          endif

! -ve y face                                                                                                                                                                                  
! -x -z 1  (1)                                                                                                                                                                                
! -x +z 5  (2)                                                                                                                                                                                
! +x -z 2  (3)                                                                                                                                                                                
! +x +z 6  (4)                                                                                                                                                                                


          if (direction%y < -0.9d0) then
             select case(subcell)
                case(1) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 34)
                case(5) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 35)
                case(2) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 36)
                case(6) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 37)
                case DEFAULT
                   write(*,*) "error in logic y -ve ",subcell
             end select
          endif

! +ve y face                                                                                                                                                                                  
! -x -z 3  (1)                                                                                                                                                                                
! -x +z 7  (2)                                                                                                                                                                                
! +x -z 4  (3)                                                                                                                                                                                
! +x +z 8  (4)                  

          if (direction%y > 0.9d0) then
             select case(subcell)
                case(3) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 34)
                case(7) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 35)
                case(4) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 36)
                case(8) 
                   flux(1:4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 37)
                case DEFAULT
                   write(*,*) "error in logic y +ve ",subcell
             end select
          endif






       else


          rhoe = thisOctal%mpiBoundaryStorage(subcell, nBound, 3)
          qnext = thisOctal%mpiBoundaryStorage(subcell, nBound, 8)
          phi = thisOctal%mpiBoundaryStorage(subcell, nBound, 12)
          phigas = thisOctal%mpiBoundaryStorage(subcell, nBound, 13)
          px = thisOctal%mpiBoundaryStorage(subcell, nBound, 15)
          py = thisOctal%mpiBoundaryStorage(subcell, nBound, 16)
          pz = thisOctal%mpiBoundaryStorage(subcell, nBound, 17)
          rm1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 21)
          rum1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 22)
          pm1 = thisOctal%mpiBoundaryStorage(subcell, nBound, 23)


          qViscosity(1,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 24)
          qViscosity(1,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 25)
          qViscosity(1,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 26)
          qViscosity(2,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 27)
          qViscosity(2,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 28)
          qViscosity(2,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 29)
          qViscosity(3,1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 30)
          qViscosity(3,2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 31)
          qViscosity(3,3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 32)



          flux(1) = thisOctal%mpiBoundaryStorage(subcell, nBound, 34)
          flux(2) = thisOctal%mpiBoundaryStorage(subcell, nBound, 35)
          flux(3) = thisOctal%mpiBoundaryStorage(subcell, nBound, 36)
          flux(4) = thisOctal%mpiBoundaryStorage(subcell, nBound, 37)

          rho(1) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 38)
          rho(2) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 39)
          rho(3) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 40)
          rho(4) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 41)

          q(1) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 42)
          q(2) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 43)
          q(3) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 44)
          q(4) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 45)

          rhou(1) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 46)
          rhou(2) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 47)
          rhou(3) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 48)
          rhou(4) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 49)

          rhov(1) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 50)
          rhov(2) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 51)
          rhov(3) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 52)
          rhov(4) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 53)

          rhow(1) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 54)
          rhow(2) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 55)
          rhow(3) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 56)
          rhow(4) =  thisOctal%mpiBoundaryStorage(subcell, nBound, 57)

       endif

    endif
  end subroutine getNeighbourValues4


  subroutine averageValue(direction, neighbourOctal, neighbourSubcell, q, rhou, rhov, rhow, rho, &
       rhoe, pressure, flux, phi, phigas, rm1, rum1, pm1, qViscosity, temperature)
    use inputs_mod, only :  cylindricalHydro

    type(OCTAL), pointer ::  neighbourOctal
    integer :: nSubcell(4), neighbourSubcell
    real(double), intent(out) :: q, rho, rhov, rhou, rhow, pressure, flux, rhoe, phi, phigas, qViscosity(3,3)
    real(double) :: fac, rm1, rum1, pm1, temperature,fac1,fac2
    type(VECTOR) :: direction,rvec

!    direction = neighbourOctal%centre - subcellCentre(thisOctal, subcell)
!    call normalize(direction)

    if (neighbourOctal%oneD) then
       if (direction%x > 0.9d0) then
          nSubcell(1) = 1
       else
          nSubcell(1) = 2
       endif
       q = neighbourOctal%q_i(neighbourSubcell)
       rho = neighbourOctal%rho(neighbourSubcell)
       rhoe = neighbourOctal%rhoe(neighbourSubcell)
       rhou = neighbourOctal%rhou(neighbourSubcell)
       rhov = neighbourOctal%rhov(neighbourSubcell)
       rhow = neighbourOctal%rhow(neighbourSubcell)
       pressure = neighbourOctal%pressure_i(neighbourSubcell)
       temperature = neighbourOctal%temperature(neighbourSubcell)
       flux = neighbourOctal%flux_i(neighbourSubcell)
       phi = neighbourOctal%phi_i(neighbourSubcell)
       phigas = neighbourOctal%phi_gas(neighbourSubcell)
       qViscosity = neighbourOctal%qViscosity(neighbourSubcell, :, :)
       rm1 = neighbourOctal%rho_i_minus_1(neighbourSubcell)
       rum1 = neighbourOctal%u_i_minus_1(neighbourSubcell)
       pm1 = neighbourOctal%pressure_i_minus_1(neighbourSubcell)
       
     else if (neighbourOctal%twoD) then
        fac1 = 0.5d0
        fac2 = 0.5d0
       if (direction%x > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 3
       else if (direction%x < -0.9d0) then
          nSubcell(1) = 2
          nSubcell(2) = 4
       else if (direction%z > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 2
          if (cylindricalHydro) then
             rVec = subcellCentre(neighbourOctal,1)             
             fac1 = pi*((rVec%x+neighbourOctal%subcellSize/2.d0)**2 - (rVec%x-neighbourOctal%subcellsize/2.d0)**2)
             rVec = subcellCentre(neighbourOctal,2)             
             fac2 = pi*((rVec%x+neighbourOctal%subcellSize/2.d0)**2 - (rVec%x-neighbourOctal%subcellsize/2.d0)**2)
             fac = fac1+fac2
             fac1 = fac1/fac
             fac2 = fac2/fac
          endif
       else
          nSubcell(1) = 3
          nSubcell(2) = 4
          if (cylindricalHydro) then
             rVec = subcellCentre(neighbourOctal,3)             
             fac1 = pi*((rVec%x+neighbourOctal%subcellSize/2.d0)**2 - (rVec%x-neighbourOctal%subcellsize/2.d0)**2)
             rVec = subcellCentre(neighbourOctal,4)             
             fac2 = pi*((rvec%x+neighbourOctal%subcellSize/2.d0)**2 - (rVec%x-neighbourOctal%subcellsize/2.d0)**2)
             fac = fac1+fac2
             fac1 = fac1/fac
             fac2 = fac2/fac
          endif
       endif

!       if (inSubcell(neighbourOctal, 3, VECTOR(50.d3, 0.d0, 1.07d6))) then
!          if (myHydroSetGlobal ==0) write(*,*) "fac1,fac2" , fac1,fac2
!       endif

       q = fac1*neighbourOctal%q_i(nSubcell(1)) + fac2*neighbourOctal%q_i(nSubcell(2))
       rho = fac1*neighbourOctal%rho(nSubcell(1)) + fac2*neighbourOctal%rho(nSubcell(2))
       rhoe = fac1*neighbourOctal%rhoe(nSubcell(1)) + fac2*neighbourOctal%rhoe(nSubcell(2))
       rhou = fac1*neighbourOctal%rhou(nSubcell(1)) + fac2*neighbourOctal%rhou(nSubcell(2))
       rhov = fac1*neighbourOctal%rhov(nSubcell(1)) + fac2*neighbourOctal%rhov(nSubcell(2))
       rhow = fac1*neighbourOctal%rhow(nSubcell(1)) + fac2*neighbourOctal%rhow(nSubcell(2))
       pressure = fac1*neighbourOctal%pressure_i(nSubcell(1)) + fac2*neighbourOctal%pressure_i(nSubcell(2))
       temperature = fac1*neighbourOctal%temperature(nSubcell(1)) + fac2*neighbourOctal%temperature(nSubcell(2))
       flux = fac1*neighbourOctal%flux_i(nSubcell(1)) + fac2*neighbourOctal%flux_i(nSubcell(2))
       phi = fac1*neighbourOctal%phi_i(nSubcell(1)) + fac2*neighbourOctal%phi_i(nSubcell(2))
       phigas = fac1*neighbourOctal%phi_gas(nSubcell(1)) + fac2*neighbourOctal%phi_gas(nSubcell(2))
       rm1 = fac1*neighbourOctal%rho_i_minus_1(nSubcell(1)) + fac2*neighbourOctal%rho_i_minus_1(nSubcell(2))
       rum1 = fac1*neighbourOctal%u_i_minus_1(nSubcell(1)) + fac2*neighbourOctal%u_i_minus_1(nSubcell(2))
       pm1 = fac1*neighbourOctal%pressure_i_minus_1(nSubcell(1)) + fac2*neighbourOctal%pressure_i_minus_1(nSubcell(2))
       qViscosity = fac1*neighbourOctal%qViscosity(nSubcell(1),:,:) + fac2*neighbourOctal%qViscosity(nSubcell(2),:,:)

    else if (neighbourOctal%threed) then
       if (direction%x > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 3
          nSubcell(3) = 5
          nSubcell(4) = 7
       else if (direction%x < -0.9d0) then
          nSubcell(1) = 2
          nSubcell(2) = 4
          nSubcell(3) = 6
          nSubcell(4) = 8
       else if (direction%y > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 2
          nSubcell(3) = 5
          nSubcell(4) = 6
       else if (direction%y < -0.9d0) then
          nSubcell(1) = 3
          nSubcell(2) = 4
          nSubcell(3) = 7
          nSubcell(4) = 8
       else if (direction%z > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 2
          nSubcell(3) = 3
          nSubcell(4) = 4
       elseif (direction%z < -0.9d0) then
          nSubcell(1) = 5
          nSubcell(2) = 6
          nSubcell(3) = 7
          nSubcell(4) = 8
       endif
       fac = 0.25d0
       q = fac*(neighbourOctal%q_i(nSubcell(1)) + neighbourOctal%q_i(nSubcell(2)) + & 
            neighbourOctal%q_i(nSubcell(3)) + neighbourOctal%q_i(nSubcell(4)))

       rho = fac*(neighbourOctal%rho(nSubcell(1)) + neighbourOctal%rho(nSubcell(2)) + & 
            neighbourOctal%rho(nSubcell(3)) + neighbourOctal%rho(nSubcell(4)))

       rhoe = fac*(neighbourOctal%rhoe(nSubcell(1)) + neighbourOctal%rhoe(nSubcell(2)) + & 
            neighbourOctal%rhoe(nSubcell(3)) + neighbourOctal%rhoe(nSubcell(4)))

       rhou = fac*(neighbourOctal%rhou(nSubcell(1)) + neighbourOctal%rhou(nSubcell(2)) + & 
            neighbourOctal%rhou(nSubcell(3)) + neighbourOctal%rhou(nSubcell(4)))

       rhov = fac*(neighbourOctal%rhov(nSubcell(1)) + neighbourOctal%rhov(nSubcell(2)) + & 
            neighbourOctal%rhov(nSubcell(3)) + neighbourOctal%rhov(nSubcell(4)))

       rhow = fac*(neighbourOctal%rhow(nSubcell(1)) + neighbourOctal%rhow(nSubcell(2)) + & 
            neighbourOctal%rhow(nSubcell(3)) + neighbourOctal%rhow(nSubcell(4)))

       temperature = fac*(neighbourOctal%temperature(nSubcell(1)) + neighbourOctal%temperature(nSubcell(2)) + & 
            neighbourOctal%temperature(nSubcell(3)) + neighbourOctal%temperature(nSubcell(4)))

       pressure = fac*(neighbourOctal%pressure_i(nSubcell(1)) + neighbourOctal%pressure_i(nSubcell(2)) + & 
            neighbourOctal%pressure_i(nSubcell(3)) + neighbourOctal%pressure_i(nSubcell(4)))

       flux = fac * (neighbourOctal%flux_i(nSubcell(1)) + neighbourOctal%flux_i(nSubcell(2)) + & 
            neighbourOctal%flux_i(nSubcell(3)) + neighbourOctal%flux_i(nSubcell(4)))


       phi = fac*(neighbourOctal%phi_i(nSubcell(1)) + neighbourOctal%phi_i(nSubcell(2)) + & 
            neighbourOctal%phi_i(nSubcell(3)) + neighbourOctal%phi_i(nSubcell(4)))

       phigas = fac*(neighbourOctal%phi_gas(nSubcell(1)) + neighbourOctal%phi_gas(nSubcell(2)) + & 
            neighbourOctal%phi_gas(nSubcell(3)) + neighbourOctal%phi_gas(nSubcell(4)))

       rm1 = fac*(neighbourOctal%rho_i_minus_1(nSubcell(1)) + neighbourOctal%rho_i_minus_1(nSubcell(2)) + & 
            neighbourOctal%rho_i_minus_1(nSubcell(3)) + neighbourOctal%rho_i_minus_1(nSubcell(4)))

       rum1 = fac*(neighbourOctal%u_i_minus_1(nSubcell(1)) + neighbourOctal%u_i_minus_1(nSubcell(2)) + & 
            neighbourOctal%u_i_minus_1(nSubcell(3)) + neighbourOctal%u_i_minus_1(nSubcell(4)))

       pm1 = fac*(neighbourOctal%pressure_i_minus_1(nSubcell(1)) + neighbourOctal%pressure_i_minus_1(nSubcell(2)) + & 
            neighbourOctal%pressure_i_minus_1(nSubcell(3)) + neighbourOctal%pressure_i_minus_1(nSubcell(4)))

       qViscosity = fac*(neighbourOctal%qViscosity(nSubcell(1),:,:) + neighbourOctal%qViscosity(nSubcell(2),:,:) + & 
            neighbourOctal%qViscosity(nSubcell(3),:,:) + neighbourOctal%qViscosity(nSubcell(4),:,:))

    endif
    
  end subroutine averageValue



  function getNBoundFromDirection(direction) result (nBound)
    type(VECTOR) :: direction
    integer :: nBound

    if (direction%x > 0.9d0) then
       nBound = 2
    else if (direction%x < -0.9d0) then
       nBound = 1
    else if (direction%z < -0.9d0) then
       nBound = 4
    else if (direction%z > 0.9d0) then
       nBound = 3
    else if (direction%y > 0.9d0) then
       nBound = 5
    else if (direction%y < -0.9d0) then
       nBound = 6
    endif
  end function getNBoundFromDirection

  function getNCornerBoundFromDirection(direction, index) result (nCornerBound)
      type(VECTOR) :: direction
      integer :: index, nCornerBound

    if(index == 1) then
       if (direction%x > 0.9d0) then
          nCornerBound = 3
       else if (direction%x < -0.9d0) then
          nCornerBound = 1        
       else if (direction%z < -0.9d0) then
          nCornerBound = 2
       else if (direction%z > 0.9d0) then
          nCornerBound = 3
       else if (direction%y > 0.9d0) then
          nCornerBound = 5
       else if (direction%y < -0.9d0) then
          nCornerBound = 8
       endif
   else if (index == 2) then
       if (direction%x > 0.9d0) then
          nCornerBound = 2
       else if (direction%x < -0.9d0) then
          nCornerBound = 4
       else if (direction%z < -0.9d0) then
          nCornerBound = 4
       else if (direction%z > 0.9d0) then
          nCornerBound = 1
       else if (direction%y > 0.9d0) then
          nCornerBound = 7
       else if (direction%y < -0.9d0) then
          nCornerBound = 6
       endif
   else if(index == 0) then
      !Special case for when the diagonal is given
      if (direction%x < -0.9d0 .and. direction%z > 0.9d0) then
          nCornerBound = 1
      else if (direction%x > 0.9d0 .and. direction%z < -0.9d0) then
          nCornerBound = 2
      else if (direction%x > 0.9d0 .and. direction%z > 0.9d0) then
          nCornerBound = 3
      else if (direction%x > -0.9d0 .and. direction%z < -0.9d0) then
          nCornerBound = 4
      else if (direction%y > 0.9d0 .and. direction%z > 0.9d0) then
          nCornerBound = 5
      else if (direction%y < -0.9d0 .and. direction%z < -0.9d0) then
          nCornerBound = 6
      else if (direction%y < -0.9d0 .and. direction%z > 0.9d0) then
          nCornerBound = 7
      else if (direction%y > 0.9d0 .and. direction%z < -0.9d0) then
          nCornerBound = 8
      end if
   end if

  end function

  subroutine columnAlongPathAMR(grid, rVec, direction, sigma)
    use mpi
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, direction, currentPosition
    real(double) :: sigma, distToNextCell
    type(OCTAL), pointer :: thisOctal, sOctal
    real(double) :: fudgeFac = 1.d-3
    integer :: subcell
    real(double) ::  totDist

    sigma = 0.d0
    currentPosition = rVec
    totDist = 0.d0

    CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)


    do while (inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition,thisOctal,subcell)

       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
  
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       totDist = totDist + distToNextCell
       if (myrankGlobal == thisOctal%mpiThread(subcell)) then
          sigma = sigma + distToNextCell*thisOctal%rho(subcell)
       endif

    end do
  end subroutine columnAlongPathAMR

  subroutine countSubcellsMPI(grid, nSubcells, nSubcellArray, includeGhosts)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: nSubcells
    integer, optional :: nSubcellArray(:)
    integer :: ierr, iThread, nVoxels, n
    integer :: iBuffer(1)
    integer, allocatable :: iArray(:)
    integer :: myRank, nThreads
    integer :: tag = 81
    integer :: status(MPI_STATUS_SIZE)
    logical, optional :: includeGhosts
    logical :: doIncludeGhosts
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    call MPI_COMM_RANK(amrCOMMUNICATOR, myRank, ierr)
    call MPI_COMM_SIZE(amrCOMMUNICATOR, nThreads, ierr)

    allocate(iArray(1:nThreads))

    doIncludeGhosts = .true.
    if (present(IncludeGhosts)) doIncludeGhosts = includeGhosts

    nSubcells = 0
    if (myrank == 0) then
       call countVoxelsMPI(grid%octreeRoot,nVoxels, doIncludeGhosts)
       nSubcells = nSubcells + nVoxels
       iArray(1) = nSubcells
       do iThread = 1, nThreads - 1
          call MPI_RECV(n, 1, MPI_INTEGER, iThread, tag, amrCOMMUNICATOR, status, ierr)
          iArray(iThread+1) = n
          nSubcells = nSubcells + n
       end do
    else
       call countVoxelsMPI(grid%octreeRoot,nVoxels, doIncludeGhosts)
       call MPI_SEND(nVoxels, 1, MPI_INTEGER, 0, tag, amrCOMMUNICATOR, ierr)
    endif
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    iBuffer(1) = nSubcells
   call MPI_BCAST(iBuffer, 1, MPI_INTEGER, 0, amrCOMMUNICATOR, ierr)
    nSubcells = iBuffer(1)

    call MPI_BCAST(iArray, nThreads, MPI_INTEGER,0, amrCOMMUNICATOR, ierr)

    if (present(nSubcellArray)) nSubcellArray = iArray
    deallocate(iArray)

  end subroutine countSubcellsMPI

  SUBROUTINE countVoxelsMPI(thisOctal, nVoxels, includeGhosts)  
    ! count the number of octals in the current section of the grid.
    ! also counts the number of unique volume elements (voxels) i.e.
    !   those subcells that are not subdivided

    IMPLICIT NONE

    TYPE(OCTAL), POINTER  :: thisOctal 
    INTEGER,INTENT(INOUT) :: nVoxels   ! number of childless subcells
    logical :: includeGhosts
    nVoxels = 0
    CALL countVoxelsPrivate(thisOctal, nVoxels, includeGhosts)
    
    CONTAINS
    
      RECURSIVE SUBROUTINE countVoxelsPrivate(thisOctal, nVoxels, includeGhosts)
        use mpi
        use inputs_mod, only : hydrodynamics, cylindricalHydro
        type(VECTOR) :: rVec
        integer :: nVoxels
        integer :: subcell
        TYPE(OCTAL), POINTER  :: thisOctal 
        TYPE(OCTAL), POINTER  :: child
        INTEGER :: i
        logical :: includeGhosts



      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call countVoxelsPrivate(child, nVoxels, includeGhosts)
                  exit
               end if
            end do
         else

            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
            if (hydrodynamics) then
               rVec = subcellCentre(thisOctal, subcell)
               if (.not.includeGhosts) then
                  if (thisOctal%ghostCell(subcell)) then
                     if (cylindricalHydro) then
                        if (rVec%x < 0.d0) cycle
                     else
                        cycle
                     endif
                  endif
               endif
            endif
            nVoxels = nVoxels + 1
         endif
      enddo


      end SUBROUTINE countVoxelsPrivate

  END SUBROUTINE countVoxelsMPI

  subroutine periodBoundary(grid, justGrav)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: iThread
    integer :: ierr, i
    real(double) :: loc(3)
    integer :: tag = 78
    logical, optional :: justGrav

    logical :: doJustGrav

    doJustGrav = .false.
    if (PRESENT(justGrav)) doJustGrav = justGrav

    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    do iThread = 1, nHydroThreadsGlobal
       if (iThread /= myRankGlobal) then
          call periodBoundaryReceiveRequests(grid, iThread, doJustGrav)
       else
          call recursivePeriodSend(grid%octreeRoot, doJustGrav)
          loc(1) = 1.d30
          do i = 1, nHydroThreadsGlobal
             if (i /= iThread) then
                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, ierr)
             endif
          enddo
       endif
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    enddo
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
  end subroutine periodBoundary

  subroutine periodBoundaryLevel(grid, nDepth, justGrav)
    use mpi
    type(GRIDTYPE) :: grid
    integer :: iThread
    integer :: ierr, i, nDepth
    real(double) :: loc(3)
    integer :: tag = 78
    logical, optional :: justGrav

    logical :: doJustGrav

    doJustGrav = .false.
    if (PRESENT(justGrav)) doJustGrav = justGrav

    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    do iThread = 1, nhydroThreadsGlobal
       if (iThread /= myRankGlobal) then
!          write(*,*) myRankGlobal, " calling boundaryreceiverequests"
          call periodBoundaryReceiveRequestsLevel(grid, iThread, nDepth, doJustGrav)
       else
!          write(*,*) "now doing ", myRankGlobal
          call recursivePeriodSendLevel(grid%octreeRoot, nDepth, doJustGrav)
          loc(1) = 1.d30
          do i = 1, nHydroThreadsGlobal
             if (i /= iThread) then
!                write(*,*) myRankGlobal, " sending terminate to ", i
                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, i, tag, localWorldCommunicator, ierr)
             endif
          enddo
       endif
!       write(*,*) myrankGlobal, " waiting at barrier"
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!       write(*,*) myrankGlobal, " dropped through barrier"
    enddo
!    write(*,*) myRankGlobal, " HAS REACHED THE BARRIER"
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
  end subroutine periodBoundaryLevel

  recursive subroutine recursivePeriodSend(thisOctal, doJustGrav)

    use mpi
    type(octal), pointer   :: thisOctal, tOctal, child
    integer :: tSubcell
    real(double) :: loc(3), tempStorage(9)
    integer :: subcell, i
    logical :: doJustGrav
    integer :: tag1 = 78, tag2 = 79
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)


    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call recursivePeriodSend(child, doJustGrav)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if ( (thisOctal%ghostCell(subcell).and.thisOctal%boundaryCondition(subcell)==2) .or. &
               (thisOctal%ghostCell(subcell).and.doJustGrav) ) then
             if (.not.doJustGrav) then
                loc(1) = thisOctal%boundaryPartner(subcell)%x
                loc(2) = thisOctal%boundaryPartner(subcell)%y
                loc(3) = thisOctal%boundaryPartner(subcell)%z
             else
                loc(1) = thisOctal%gravboundaryPartner(subcell)%x
                loc(2) = thisOctal%gravboundaryPartner(subcell)%y
                loc(3) = thisOctal%gravboundaryPartner(subcell)%z
             endif

             tOctal => thisOctal
             tSubcell = 1
             if (.not.doJustGrav) then
                call findSubcellLocal(thisOctal%boundaryPartner(subcell), tOctal,tsubcell)
             else
                call findSubcellLocal(thisOctal%gravboundaryPartner(subcell), tOctal,tsubcell)
             endif
            ! write(*,*) "boundary partner ", thisOctal%boundaryPartner(subcell)
            ! write(*,*) myrankGlobal, " sending locator to ", tOctal%mpiThread(tsubcell)
             call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag1, localWorldCommunicator, ierr)
            ! write(*,*) myRankGlobal, " awaiting recv from ", tOctal%mpiThread(tsubcell)
             call MPI_RECV(tempStorage, 9, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), &
                  tag2, localWorldCommunicator, status, ierr)
            ! write(*,*) myrankglobal, " received from ",tOctal%mpiThread(tSubcell)
             if (.not.associated(thisOctal%tempStorage)) then
                if (.not.doJustGrav) then
                   allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:9))
                   thisOctal%tempStorage = 0.d0
                else
                   allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:1))
                   thisOctal%tempStorage = 0.d0
                endif
             endif
             thisOctal%tempStorage(subcell,1:SIZE(thisOctal%tempStorage,2)) = &
                  tempStorage(1:SIZE(thisOctal%tempStorage,2))
          endif
       endif
    enddo
  end subroutine recursivePeriodSend

  recursive subroutine recursivePeriodSendLevel(thisOctal, nDepth, doJustGrav)

    use mpi
    type(octal), pointer   :: thisOctal, tOctal, child
    integer :: tSubcell
    integer :: nDepth
    real(double) :: loc(3), tempStorage(9)
    logical :: doJustGrav
    integer :: subcell, i
    integer :: tag1 = 78, tag2 = 79
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)


    if ((thisOctal%nChildren > 0).and.(thisOctal%nDepth < nDepth)) then
       do i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          call recursivePeriodSendLevel(child, nDepth, doJustGrav)
       end do
    else

       do subcell = 1, thisOctal%maxChildren
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if ( (thisOctal%ghostCell(subcell).and.thisOctal%boundaryCondition(subcell)==2) .or. &
               (thisOctal%ghostCell(subcell).and.doJustGrav) ) then
             if (.not.doJustGrav) then
                loc(1) = thisOctal%boundaryPartner(subcell)%x
                loc(2) = thisOctal%boundaryPartner(subcell)%y
                loc(3) = thisOctal%boundaryPartner(subcell)%z
             else
                loc(1) = thisOctal%gravboundaryPartner(subcell)%x
                loc(2) = thisOctal%gravboundaryPartner(subcell)%y
                loc(3) = thisOctal%gravboundaryPartner(subcell)%z
             endif

             tOctal => thisOctal
             tSubcell = 1
             if (.not.dojustGrav) then
                call findSubcellLocalLevel(thisOctal%boundaryPartner(subcell), tOctal,tsubcell, nDepth)
             else
                call findSubcellLocalLevel(thisOctal%gravboundaryPartner(subcell), tOctal,tsubcell, nDepth)
             endif
             if (tOctal%mpiThread(tSubcell) == myRankGLobal) then
                write(*,*) "bug in recursiveperiodsendlevel ",thisOctal%gravBoundaryPartner(subcell)
             endif
!             write(*,*) "boundary partner ", thisOctal%boundaryPartner(subcell)
!             write(*,*) myrankGlobal, " sending locator to ", tOctal%mpiThread(tsubcell)
             call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag1, localWorldCommunicator, ierr)
!             write(*,*) myRankGlobal, " awaiting recv from ", tOctal%mpiThread(tsubcell)
             call MPI_RECV(tempStorage, 7, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag2, &
                  localWorldCommunicator, status, ierr)
!             write(*,*) myrankglobal, " received from ",tOctal%mpiThread(tSubcell)
             if (.not.associated(thisOctal%tempStorage)) then
                if (.not.doJustGrav) then
                   allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:9))
                   thisOctal%tempStorage = 0.d0
                else
                   allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:1))
                   thisOctal%tempStorage = 0.d0
                endif
             endif
             thisOctal%tempStorage(subcell,1:SIZE(thisOctal%tempStorage,2)) = &
                  tempStorage(1:SIZE(thisOctal%tempStorage,2))
          endif
       enddo
    endif
  end subroutine recursivePeriodSendLevel

  subroutine periodBoundaryReceiveRequests(grid, receiveThread, doJustGrav)
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    logical :: sendLoop
    real(double) :: loc(3), tempStorage(9)
    integer :: ierr, receiveThread
    integer :: status(MPI_STATUS_SIZE)
    integer :: subcell
    integer :: tag1 = 78, tag2 = 79
    logical :: doJustGrav
    type(VECTOR) :: octVec
    sendLoop = .true.
    !write(*,*) myrankGlobal, " waiting for a locator"
    do while (sendLoop)
       ! receive a locator
       
       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, receiveThread, tag1, localWorldCommunicator, status, ierr)
       !write(*,*) myrankglobal, " received a locator from ", receiveThread 
       if (loc(1) > 1.d20) then
          sendLoop = .false.
        !  write(*,*) myRankGlobal, " found the signal to end the send loop from ", receivethread
       else
          octVec = VECTOR(loc(1), loc(2), loc(3))
          thisOctal => grid%octreeRoot
          subcell = 1
         ! write(*,*) myrankglobal," calling subscell local"
          call findSubcellLocal(octVec, thisOctal, subcell)
         ! write(*,*) myrankglobal," subscell local done succesfully"
          
          if (.not.doJustGrav) then
             tempstorage(1) = thisOctal%rho(Subcell)
             tempStorage(2) = thisOctal%rhoE(Subcell)
             tempStorage(3) = thisOctal%rhou(Subcell)
             tempStorage(4) = thisOctal%rhov(Subcell)
             tempStorage(5) = thisOctal%rhow(Subcell)
             tempStorage(6) = thisOctal%energy(Subcell)
             tempStorage(7) = thisOctal%pressure_i(Subcell)
             tempStorage(9) = thisOctal%temperature(subcell)
             call MPI_SEND(tempStorage, 9, MPI_DOUBLE_PRECISION, receiveThread, tag2, localWorldCommunicator, ierr)
          else
!             tempStorage(1) = thisOctal%phi_i(Subcell)
             tempStorage(1) = thisOctal%phi_gas(Subcell)
             call MPI_SEND(tempStorage, 9, MPI_DOUBLE_PRECISION, receiveThread, tag2, localWorldCommunicator, ierr)
          endif
       endif
    enddo
    !write(*,*) myrankGlobal, " leaving receive requests ", sendLoop
  end subroutine periodBoundaryReceiveRequests


  subroutine periodBoundaryReceiveRequestsLevel(grid, receiveThread, nDepth, doJustGrav)
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    logical :: sendLoop
    integer :: nDepth
    real(double) :: loc(3), tempStorage(9)
    integer :: ierr, receiveThread
    integer :: status(MPI_STATUS_SIZE)
    integer :: subcell
    integer :: tag1 = 78, tag2 = 79
    logical :: doJustGrav
    type(VECTOR) :: octVec
    sendLoop = .true.
!    write(*,*) myrankGlobal, " waiting for a locator"
    do while (sendLoop)
       ! receive a locator
       
       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, receiveThread, tag1, localWorldCommunicator, status, ierr)
!       write(*,*) myrankglobal, " received a locator from ", receiveThread 
       if (loc(1) > 1.d20) then
          sendLoop = .false.
!          write(*,*) myRankGlobal, " found the signal to end the send loop from ", receivethread
       else
          octVec = VECTOR(loc(1), loc(2), loc(3))
          thisOctal => grid%octreeRoot
          subcell = 1
!          write(*,*) myrankglobal," calling subscell local"
          call findSubcellLocalLevel(octVec, thisOctal, subcell, nDepth)
!          write(*,*) myrankglobal," subscell local done succesfully"
          
          if (.not.doJustGrav) then
             tempstorage(1) = thisOctal%rho(Subcell)
             tempStorage(2) = thisOctal%rhoE(Subcell)
             tempStorage(3) = thisOctal%rhou(Subcell)
             tempStorage(4) = thisOctal%rhov(Subcell)
             tempStorage(5) = thisOctal%rhow(Subcell)
             tempStorage(6) = thisOctal%energy(Subcell)
             tempStorage(7) = thisOctal%pressure_i(Subcell)
             call MPI_SEND(tempStorage, 7, MPI_DOUBLE_PRECISION, receiveThread, tag2, localWorldCommunicator, ierr)
          else
!             tempStorage(1) = thisOctal%phi_i(Subcell)
             tempStorage(1) = thisOctal%phi_gas(Subcell)
             call MPI_SEND(tempStorage, 7, MPI_DOUBLE_PRECISION, receiveThread, tag2, localWorldCommunicator, ierr)
          endif
       endif
    enddo
!    write(*,*) myrankGlobal, " leaving receive requests ", sendLoop
  end subroutine periodBoundaryReceiveRequestsLevel


!THaw-to track evolution of I front with time
subroutine dumpStromgrenRadius(grid, thisFile, startPoint, endPoint, nPoints)
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    integer :: nPoints
    type(VECTOR) :: startPoint, endPoint, position, direction, cen
    real(double) :: loc(3), hi
    character(len=*) :: thisFile
    integer :: ierr
    integer, parameter :: nStorage = 5
    real(double) :: tempSTorage(nStorage), tval
    integer, parameter :: tag = 30
    integer :: status(MPI_STATUS_SIZE)
    logical :: stillLooping, done
    integer :: sendThread
    integer :: i
    logical, save :: firstTime =.true.
    real(double) :: oldR, newR, oldHi, r
    i = npoints

    thisOctal => grid%octreeRoot
    position = startPoint
    direction = endPoint - startPoint
    call normalize(direction)

    if (myHydroSetGlobal /= 0) goto 666

    if (myrankWorldGlobal == 0) then
       if(firstTime) then
          !Overwrite any existing file
          open(20, file=thisFile, form="formatted", status="unknown")
          write(20,*) " "
          close(20)
          firstTime = .false.
       end if
       !Append new positions to a new file
       open(20, file=thisFile, form="formatted", status="unknown", position="append")
       done = .false.
       !print *, "Starting dump of Stromgren data"

       oldR = 0.d0
       newR = 0.d0
       oldHi = 0.d0
       do while(inOctal(grid%octreeRoot, position))
          call findSubcellLocal(position, thisOctal, subcell)
          sendThread = thisOctal%mpiThread(subcell)
          loc(1) = position%x
          loc(2) = position%y
          loc(3) = position%z
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
          call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)

          cen%x = tempStorage(1)
          cen%y = tempStorage(2)
          cen%z = tempStorage(3)
          hi = tempStorage(4)
          tVal = tempStorage(5)
	  !print *, "hi ", hi
	  !print *, "cen: ", cen
          newR = modulus(cen-startPoint)
          if(hi > 0.5 .and. .not. done) then
              r = oldR + (newR-oldR) * (0.5 - oldHi)/(hi-oldhi)
	     !print *, "Edge found"
!             write(20,'(5e14.5)') grid%currentTime*secsToYears/1.d6, r*1.d10/(1000.d0*pctocm)
!              if(grid%octreeroot%oned .and. grid%geometry == "bonnor") then
!                 r = r + (3.2951d9 - grid%octreeroot%subcellsize)
!              end if
              write(20,'(5e14.5)') grid%currentTime*secsToYears, r*1.d10/pctocm                 
             done = .true.
             goto 555
          end if
          oldR = newR
          oldHi= hi
          position = cen
          position = position + (tVal+1.d-3*grid%halfSmallestSubcell)*direction
	  !print *, "POSITION ", position
	  !print *, "direction ", direction
	  !print *, "tVal ", tVal

       enddo
555    continue
       !Send escape trigger to other threads
       do sendThread = 1, nHydroThreadsGlobal
          loc(1) = 1.d30
          loc(2) = 1.d30
          loc(3) = 1.d30
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
       enddo
       close(20)
       goto 666


    else
       stillLooping = .true.
       !print *, "Starting to loop rank", myRankGlobal
       do while(stillLooping)
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
          position%x = loc(1)
          position%y = loc(2)
          position%z = loc(3)
          if (position%x > 1.d29) then
             stillLooping = .false.
          else
             call findSubcellLocal(position, thisOctal, subcell)
             sOctal => thisOctal
             cen = subcellCentre(thisOctal, subcell)
             call distanceToCellBoundary(grid, cen, direction, tVal, sOctal)
             tempStorage(1) = cen%x
             tempStorage(2) = cen%y
             tempStorage(3) = cen%z
             tempStorage(4) = thisOctal%ionFrac(subcell,1)
             tempStorage(5) = tVal
             !tempStorage(6) = thisOctal%rhou(subcell)             
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
	     
          endif
       enddo
    endif

 666   continue
 !print *, "Stromgren dump completed"
end subroutine dumpStromgrenRadius

subroutine write1dlist(grid, thisFile)
  use mpi
  type(GRIDTYPE) :: grid
  character(len=*) :: thisFile
  integer :: iThread, ierr

  if (myHydroSetGlobal == 0) then

     do iThread = 1, nHydrothreadsGlobal

        if (iThread == myRankGlobal) then
           if (iThread == 1) then
              open(55, file=thisfile, form="formatted", status="unknown")
           else
              open(55, file=thisfile, form="formatted", status="old",position="append")
           endif
           call write1dListToFile(grid, grid%octreeRoot, ithread)
           close(55)
        endif

        call MPI_BARRIER(localWorldCommunicator, ierr)

     enddo

  endif
end subroutine write1dlist

recursive subroutine write1dlisttoFile(grid, thisOctal, ithread)
!  use inputs_mod, only : sphereRadius
  use source_mod, only : globalSourceArray
  type(GRIDTYPE) :: grid
  real(double) :: radPress,v, kappaAbs, kappaSca, kappaExt,tau,fac,area,dtau
  real(double) :: kappaAbsDust, kappaScaDust
  type(OCTAL), pointer :: thisOctal, child
  type(VECTOR) :: rVec
  integer :: ithread, i, subcell

  if (thisOctal%nChildren > 0) then
     do i = 1, thisOctal%nChildren, 1
        child => thisOctal%child(i)
        call write1dlisttofile(grid, child, ithread)
     end do
  else
     
     do subcell = 1, thisOctal%maxChildren
        
        if (octalOnThread(thisOctal, subcell,ithread)) then
           rVec = subcellCentre(thisOctal, subcell)
           v = cellVolume(thisOctal,subcell) * 1.d30
           call returnKappa(grid, thisOctal, subcell, ilambda=1,&
                kappaScaDust=kappaScaDust, kappaAbsDust=kappaAbsDust, &
                kappaSca=kappaSca, kappaAbs=kappaAbs)
           kappaExt = (kappaScaDust + kappaAbsDust)
           tau = 0.d0 !1000d0*(rVec%x/sphereRadius)
           fac = exp(-tau) 
           dtau = kappaExt * thisOctal%subcellSize
           area = fourPi * rVec%x**2 * 1.d20
           radpress = fac*globalSourceArray(1)%luminosity * (kappaExt/1.d10) / (cSpeed * fourPi * rVec%x**2 * 1.d20)
           if (rVec%x >= 0.d0) write(55,'(1p,7e12.3)') rVec%x, &
                thisOctal%kappaTimesFlux(subcell)%x/cSpeed, &
                thisOctal%radiationMomentum(subcell)%x, radpress, thisOctal%rho(subcell), &
                tau,kappaExt
        endif
     enddo
  endif
end subroutine write1dlisttoFile


#ifdef PDR
!Dumps the grid to a file in a format that can be used in 3D-PDR (Bisbas et al. 2012)
recursive subroutine SendGridBisbas(thisOctal, grid)
  use mpi
  type(OCTAL), pointer :: thisOctal, child
  type(VECTOR) :: rVec
  type(GRIDTYPE) :: grid
  integer :: i, subcell,  ierr, nStorage
  integer, parameter :: tag = 10
  real(double) :: tempstorage(26)
  nStorage = 26

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  SendGridBisbas(child, grid)
                exit
             end if
          end do
       else 
!          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then             
          if(octalonthread(thisoctal, subcell, myrankglobal)) then
             rVec= subcellCentre(thisOctal, subcell)               
             tempStorage(1) = rvec%x
             tempStorage(2) = rvec%y
             tempStorage(3) = rvec%z
             tempStorage(4) = thisOctal%AV(subcell, 7)
             tempStorage(5) = thisOctaL%tLast(subcell)
             tempStorage(6) = thisOctal%dust_t(subcell)
             tempstorage(7) = thisOctal%abundance(subcell, 31) !H2
             tempstorage(8) = thisOctal%abundance(subcell, 32) !H
             tempStorage(9) = thisOctal%abundance(subcell, 11) !C+ 
             tempStorage(10) = thisOctal%abundance(subcell, 25) !C
             tempStorage(11) = thisOctal%abundance(subcell, 28) !CO
             tempStorage(12) = thisOctal%ciiline(subcell, 2, 1) !cii(2-->1) emissivity
             tempStorage(13) = thisOctal%ciline(subcell, 2, 1) !ci(2-->1) emissivity
             tempStorage(14) = thisOctal%ciline(subcell, 3, 1) !ci(3-->1) emissivity
             tempStorage(15) = thisOctal%ciline(subcell, 3, 2) !ci(3-->2) emissivity
             tempStorage(16) = thisOctal%oiline(subcell, 2, 1) !oi(2-->1) emissivity
             tempStorage(17) = thisOctal%oiline(subcell, 3, 1) !oi(3-->1) emissivity
             tempStorage(18) = thisOctal%oiline(subcell, 3, 2) !oi(3-->2) emissivity
             tempStorage(19) = thisOctal%c12oline(subcell, 2, 1) !co(2-->1) emissivity
             tempStorage(20) = thisOctal%c12oline(subcell, 3, 2) !co(3-->2) emissivity
             tempStorage(21) = thisOctal%c12oline(subcell, 4, 3) !co(4-->3) emissivity
             tempStorage(22) = thisOctal%c12oline(subcell, 5, 4) !co(5-->4) emissivity
             tempStorage(23) = thisOctal%c12oline(subcell, 6, 5) !co(6-->5) emissivity
             tempStorage(24) = thisOctal%c12oline(subcell, 7, 6) !co(7-->6) emissivity
             tempStorage(25) = thisOctal%c12oline(subcell, 8, 7) !co(8-->7) emissivity
             tempStorage(26) = thisOctal%UV(Subcell)
!             tempStorage(4) = thisOctal%AV(subcell, 1)
!             tempStorage(5) = thisOctaL%tLast(subcell)
!             tempStorage(6) = thisOctal%dust_t(subcell)
!             tempStorage(7) = thisOctal%UV(subcell)             
!!             tempStorage(8) = thisOctal%heatingRate(subcell, 12)
 !            tempStorage(9) = sum(thisOctal%coolingRate(subcell,:))
 !            tempStorage(10) = 1.d0 !dummy for tval
 !            tempStorage(11) = thisOctal%abundance(subcell, 11) !C+ 
 !            tempStorage(12) = thisOctal%abundance(subcell, 25) !C
 !            tempStorage(13) = thisOctal%abundance(subcell, 28) !CO
 !            tempStorage(14) = thisOctal%coolingRate(subcell, 1)
 !            tempStorage(15) = thisOctal%coolingRate(subcell, 2)
 !            tempStorage(16) = thisOctal%coolingRate(subcell, 3)
 !            tempStorage(17) = thisOctal%coolingRate(subcell, 4)
 !            tempStorage(18) = thisOctal%cii_pop(subcell, 1)
 !            tempStorage(19) = thisOctal%cii_pop(subcell, 2)
 !            tempstorage(20) = thisOctal%columnRho(subcell)
 !            tempstorage(21) = thisOctal%thisColRho(subcell, 1, 31)
 !            tempstorage(22) = thisOctal%thisColRho(subcell, 1, 32)
 !            tempstorage(23) = sum(thisOctal%thisColRho(subcell, 1, :))
 !            tempstorage(24) = thisOctal%abundance(subcell, 31) !
 !            tempstorage(25) = thisOctal%abundance(subcell, 32) !
             nstorage = 26
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          end if
       endif
    end do

  end subroutine SendGridBisbas

  subroutine terminateBisbas()
  use mpi
  integer :: ierr, nStorage, i
  integer, parameter :: tag = 10
  real(double) :: tempstorage(26)
  nStorage = 26

  do i = 1, nhydrothreadsglobal
     if(i == myrankglobal) then
        tempStorage = 0.d0
        tempStorage(5) = 2.d30
        call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
     end if
  end do

end subroutine terminateBisbas

  subroutine writeGridToBisbas()!(grid)
    use mpi
    integer, parameter :: nstorage=26
!    type(GRIDTYPE) :: grid
    integer :: numStillSending
    type(VECTOR) :: position
    real(double) :: tempstorage(nstorage)
    real(double) :: c, cplus, c12o!, cii_1, cii_2, H2Col, HCol
!    real(double) :: UVtemp, heating, cooling, column, totcol
    real(double) :: AVtemp, UV, tlast, dustt, uvtemp
    integer :: status, ierr
    integer, parameter :: tag = 10
    logical :: stillRecving
    real(double) :: H, H2, CII21, CI21, CI31, CI32, OI21, OI31, OI32, CO21, CO32, CO43, &
         CO54, CO65, CO76, CO87, tgas, tdust
    open(123, file="PDR_result.txt", form="formatted", status="unknown")
    stillRecving = .true.
    numStillSending = nHydroThreadsGlobal
    print *, "num still sending: ", numstillsending
    do while (stillRecving)
       
       call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, localWorldCommunicator, status, ierr)
       if(tempStorage(5) > 1.d30) then          
          numStillSending = numstillSending - 1
          print *, "num still sending: ", numstillsending
          if(numStillSEnding == 0) then
             stillRecving = .false.
             print *, "TIME TO STOP"
          end if
       else

       position%x = tempStorage(1)
       position%y = tempStorage(2)
       position%z = tempStorage(3)
       AVtemp = tempStorage(4)
       tlast = tempStorage(5)
       dustT = tempStorage(6)
       UVtemp = tempStorage(26)
       H = tempStorage(7)
       H2 = tempStorage(8)
       cplus = tempstorage(9)
       c = tempstorage(10)
       c12o = tempstorage(11)
       CII21 = tempstorage(12)
       CI21 = tempstorage(13)
       CI31 = tempstorage(14)
       CI32 = tempstorage(15)
       OI21 = tempstorage(16)
       OI31 = tempstorage(17)
       OI32 = tempstorage(18)
       CO21 = tempstorage(19)
       CO32 = tempstorage(20)
       CO43 = tempstorage(21)
       CO54 = tempstorage(22)
       CO65 = tempstorage(23)
       CO76 = tempstorage(24)
       CO87 = tempstorage(25)
       UV = tempstorage(26)
       
       write(123, '(1p,40e14.5)') ((position%x*1.d10/pctocm)), &
            ((position%y*1.d10/pctocm)), &
            ((position%z*1.d10/pctocm)), &
            AVtemp, UV, tGas, tDust, H, H2, Cplus, C, C12O, CII21, CI21, CI31, &
            CI32, OI21, OI31, OI32, CO21, CO32, CO43, CO54, CO65, CO76, CO87
            !AVtemp, tlast, dustT, UVtemp, heating, cooling, c, cplus, c12o, &
            !cool1, cool2, cool3, cool4, cii_1, cii_2, column, H2Col, HCol, totCol, H2pop, Hpop
!
!          write(20,'(1p,40e14.5)') cen%x

    end if
 end do
  end subroutine writeGridToBisbas
#endif

subroutine writeRadialFile(rootFilename, grid)
  use inputs_mod, only : iModel
  use utils_mod, only : findMultiFilename
  character(len=*) :: rootFilename
  character(len=80) :: thisFile
  type(GRIDTYPE) :: grid
  type(VECTOR) :: startPoint, endPoint
  integer :: nPoints

  call findMultiFilename(rootFilename, iModel, thisFile)
  startPoint = VECTOR(0.d0, 0.d0, 0.d0)
  endPoint = VECTOR(grid%octreeRoot%subcellSize, grid%octreeRoot%subcellSize, grid%octreeRoot%subcellsize)
  call  dumpValuesAlongLine(grid, thisFile, startPoint, endPoint, nPoints)
end subroutine writeRadialFile

  subroutine dumpValuesAlongLine(grid, thisFile, startPoint, endPoint, nPoints)
    use mpi
    use inputs_mod, only : inputgfac, dustPhysics, hydrodynamics
    use source_mod, only : globalSourceArray
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    integer :: nPoints
    type(VECTOR) :: startPoint, endPoint, position, direction, cen, rVec
    real(double) :: loc(3), rho, rhou , rhoe, p, phi_stars, phi_gas
    real(double) :: temperature, r
    character(len=*) :: thisFile
    integer :: ierr
    integer, parameter :: nStorage = 16
    real(double) :: tempSTorage(nStorage), tval, kappaTimesFlux, rpress, radmom
    real(double) :: kappaSca, kappaAbs, kappaAbsDust, kappaScaDust
    integer, parameter :: tag = 30
    integer :: status(MPI_STATUS_SIZE)
    logical :: stillLooping
    integer :: sendThread
    integer :: i
    logical :: ghostCell

    i = npoints

    thisOctal => grid%octreeRoot
    position = startPoint
    direction = endPoint - startPoint
    call normalize(direction)
    if (myHydroSetGlobal /= 0) goto 666

    if (myrankWorldGlobal == 0) then

       open(20, file=thisFile, form="formatted", status="unknown")
       do while(inOctal(grid%octreeRoot, position))
          call findSubcellLocal(position, thisOctal, subcell)
          sendThread = thisOctal%mpiThread(subcell)
          loc(1) = position%x
          loc(2) = position%y
          loc(3) = position%z
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
          call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)

          cen%x = tempStorage(1)
          cen%y = tempStorage(2)
          cen%z = tempStorage(3)
          rho = tempStorage(4)
          tval = tempStorage(5)
          rhou = tempStorage(6)
          rhoe = tempStorage(7)
          p = tempStorage(8)
          phi_stars = tempStorage(9)
          phi_gas = tempStorage(10)
          temperature = tempStorage(11)
          kappaTimesflux = tempStorage(12)
          radmom = tempStorage(13)
          kappaAbs = tempStorage(14)
          kappaSca = tempStorage(15)
          ghostCell = .false.
          if (hydrodynamics) then
             if (tempStorage(16) == 1) then
                ghostCell = .true.
             endif
          endif

          if (.not.ghostCell) then
          if(grid%geometry == "SB_CD_1Da" .or. grid%geometry == "SB_CD_1Db") then

             if(cen%x > 0.0078125d0 .and. cen%x < (1.d0+0.0078125d0)) then                
                write(20,'(1p,7e14.5)') modulus(cen), rho, p,  rhou/rho
             end if
          else if (grid%geometry == "SB_coolshk") then
             write(20,'(1p,7e14.5)') modulus(cen), rho, rhou/rho, p, temperature/(2.33d0*mHydrogen/kerg)
          else
             rpress = globalSourceArray(1)%luminosity * ((kappaAbs+(1.d0-inputgFac)*kappaSca)/1.d10)/ &
                  (cSpeed * fourPi * modulus(cen)**2 * 1.d20)
             write(20,'(1p,11e12.4)') modulus(cen), rho, rhou/rho, rhoe,p, phi_stars, phi_gas, kappaTimesFlux, radmom, rpress, &
                  temperature
!             write(20,'(1p,7e14.5)') modulus(cen), rho, rhou/rho, rhoe,p, temperature
          end if
          endif
          position = cen
          position = position + (tVal+1.d-3*grid%halfSmallestSubcell)*direction

       enddo
       !Send escape trigger to other threads
       do sendThread = 1, nHydroThreadsGlobal
          loc(1) = 1.d30
          loc(2) = 1.d30
          loc(3) = 1.d30
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
       enddo
       close(20)
       goto 666


    else
       stillLooping = .true.
       do while(stillLooping)
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
          position%x = loc(1)
          position%y = loc(2)
          position%z = loc(3)
          if (position%x > 1.d29) then
             stillLooping = .false.
          else
             call findSubcellLocal(position, thisOctal, subcell)
             sOctal => thisOctal
             cen = subcellCentre(thisOctal, subcell)
             call distanceToCellBoundary(grid, cen, direction, tVal, sOctal)
             tempStorage = 0.d0
             tempStorage(1) = cen%x
             tempStorage(2) = cen%y
             tempStorage(3) = cen%z
             tempStorage(4) = thisOctal%rho(subcell)
             tempStorage(5) = tVal
             tempStorage(6) = sqrt(thisOctal%rhou(subcell)**2+thisOctal%rhov(subcell)**2+thisOctal%rhow(subcell)**2)
             tempStorage(7) = thisOctal%rhoe(subcell)             
             tempStorage(8) = thisOctal%pressure_i(subcell)             
             tempStorage(9) = thisOctal%phi_stars(subcell)             
             tempStorage(10) = thisOctal%phi_gas(subcell)             
             tempStorage(11) = thisOctal%temperature(subcell)
             rVec = subcellCentre(thisOctal,subcell)
             r = modulus(rVec) * 1.d10
             tempStorage(12) = modulus(thisOctal%kappaTimesFlux(subcell))/cSpeed
             tempStorage(13) = modulus(thisOctal%radiationMomentum(subcell))

             if (dustPhysics) then
                call returnKappa(grid, thisOctal, subcell, ilambda=1,&
                     kappaScaDust=kappaScaDust, kappaAbsDust=kappaAbsDust, &
                     kappaSca=kappaSca, kappaAbs=kappaAbs)
                
                tempStorage(14) = kappaAbsDust
                tempStorage(15) = kappaScaDust
             endif
             if (hydrodynamics) then
                if (thisOctal%ghostCell(subcell)) then
                   tempStorage(16) = 1
                else
                   tempStorage(16) = 0
                endif
             endif
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          endif
       enddo
    endif

666 continue
  end subroutine dumpValuesAlongLine


#ifdef PDR
  subroutine dumpValuesAlongLinePDR(grid, thisFile, startPoint, endPoint, nPoints)
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    integer :: nPoints
    type(VECTOR) :: startPoint, endPoint, position, direction, cen
    real(double) :: loc(3)!, rho, rhou , rhoe, p, phi_stars, phi_gas
!    real(double) :: temperature
    character(len=*) :: thisFile
    integer :: ierr
    integer, parameter :: nStorage = 25
    real(double) :: c, cplus, c12o, cool1, cool2, cool3, cool4, cii_1, cii_2, H2Col, HCol
    real(double) :: tempSTorage(nStorage), UVtemp, heating, cooling, column, totcol
    real(double) :: tlast, dustT, AVtemp, tval, Hpop, H2pop
    integer, parameter :: tag = 30
    integer :: status(MPI_STATUS_SIZE)
    logical :: stillLooping
    integer :: sendThread
    integer :: i

    i = npoints

    thisOctal => grid%octreeRoot
    position = startPoint
    direction = endPoint - startPoint
    call normalize(direction)
    if (myHydroSetGlobal /= 0) goto 555

    if (myrankWorldGlobal == 0) then

       open(20, file=thisFile, form="formatted", status="replace")
!       open(21, file="hc_components.dat", form="formatted", status="replace")
       do while(inOctal(grid%octreeRoot, position))
          call findSubcellLocal(position, thisOctal, subcell)
          sendThread = thisOctal%mpiThread(subcell)
          loc(1) = position%x
          loc(2) = position%y
          loc(3) = position%z
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
          call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)

          cen%x = tempStorage(1)
          cen%y = tempStorage(2)
          cen%z = tempStorage(3)
          AVtemp = tempStorage(4)
          tlast = tempStorage(5)
          dustT = tempStorage(6)
          UVtemp = tempStorage(7)
          heating = tempStorage(8)
          cooling = tempStorage(9)
          tval =  tempStorage(10)
          cplus = tempstorage(11)
          c = tempstorage(12)
          c12o = tempstorage(13)
          cool1 = tempstorage(14)
          cool2 = tempstorage(15)
          cool3 = tempstorage(16)
          cool4 = tempstorage(17)
          cii_1 = tempstorage(18)
          cii_2 = tempstorage(19)
          column = tempstorage(20)
          H2Col = tempStorage(21)
          HCol = tempstorage(22)
          totCol = tempstorage(23)
          H2pop = tempStorage(24)
          Hpop = tempstorage(25)
          write(20,'(1p,40e14.5)') cen%x, AVtemp, tlast, dustT, UVtemp, heating, cooling, c, cplus, c12o, &
               cool1, cool2, cool3, cool4, cii_1, cii_2, column, H2Col, HCol, totCol, H2pop, Hpop

          position = cen
          position = position + (tVal+1.d-3*grid%halfSmallestSubcell)*direction

       enddo
       !Send escape trigger to other threads
       do sendThread = 1, nHydroThreadsGlobal
          loc(1) = 1.d30
          loc(2) = 1.d30
          loc(3) = 1.d30
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
       enddo
       close(20)
!       goto 555

    else
       stillLooping = .true.
       do while(stillLooping)
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
          position%x = loc(1)
          position%y = loc(2)
          position%z = loc(3)
          if (position%x > 1.d29) then
             stillLooping = .false.
          else
             call findSubcellLocal(position, thisOctal, subcell)
             sOctal => thisOctal
             cen = subcellCentre(thisOctal, subcell)
             call distanceToCellBoundary(grid, cen, direction, tVal, sOctal)
             tempStorage(1) = cen%x
             tempStorage(2) = cen%y
             tempStorage(3) = cen%z
             tempStorage(4) = thisOctal%AV(subcell, 1)
             tempStorage(5) = thisOctaL%tLast(subcell)
             tempStorage(6) = thisOctal%dust_t(subcell)
             tempStorage(7) = thisOctal%UV(subcell)             
             tempStorage(8) = thisOctal%heatingRate(subcell, 12)
             tempStorage(9) = sum(thisOctal%coolingRate(subcell,:))
             tempStorage(10) = tval
             tempStorage(11) = thisOctal%abundance(subcell, 11) !C+ 
             tempStorage(12) = thisOctal%abundance(subcell, 25) !C
             tempStorage(13) = thisOctal%abundance(subcell, 28) !CO
             tempStorage(14) = thisOctal%coolingRate(subcell, 1)
             tempStorage(15) = thisOctal%coolingRate(subcell, 2)
             tempStorage(16) = thisOctal%coolingRate(subcell, 3)
             tempStorage(17) = thisOctal%coolingRate(subcell, 4)
             tempStorage(18) = thisOctal%cii_pop(subcell, 1)
             tempStorage(19) = thisOctal%cii_pop(subcell, 2)
             tempstorage(20) = thisOctal%columnRho(subcell)
             tempstorage(21) = thisOctal%thisColRho(subcell, 1, 31)
             tempstorage(22) = thisOctal%thisColRho(subcell, 1, 32)
             tempstorage(23) = sum(thisOctal%thisColRho(subcell, 1, :))
             tempstorage(24) = thisOctal%abundance(subcell, 31) !
             tempstorage(25) = thisOctal%abundance(subcell, 32) !
!             tempStorage(14) = thisOctal%heatingRate(subcell, 1)


!             tempStorage(9) = thisOctal%ci_abund(subcell)             
!             tempStorage(10) = thisOctal%oi_abund(subcell)             
!             tempStorage(11) = thisOctal%c12o_abund(subcell)
             
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
             
          Endif
       enddo
    endif

555 continue
  end subroutine dumpValuesAlongLinePDR


  subroutine dumpFinalResultsAloneLine(grid, thisFile, startPoint, endPoint, nPoints)
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    integer :: nPoints
    type(VECTOR) :: startPoint, endPoint, position, direction, cen
    real(double) :: loc(3)!, rho, rhou , rhoe, p, phi_stars, phi_gas
!    real(double) :: temperature
    character(len=*) :: thisFile
    integer :: ierr
    integer, parameter :: nStorage = 26
    real(double) :: c, cplus, c12o!, cool1, cool2, cool3, cool4, cii_1, cii_2, H2Col, HCol
    real(double) :: tempSTorage(nStorage)!, heating, cooling, column, totcol
    real(double) :: AVtemp, tval!, Hpop, H2pop
    integer, parameter :: tag = 30
    integer :: status(MPI_STATUS_SIZE)
    logical :: stillLooping
    integer :: sendThread
    integer :: i
    real(double) :: H, H2, CII21, CI21, CI31, CI32, OI21, OI31, OI32, CO21, CO32, CO43, &
         CO54, CO65, CO76, CO87, tgas, tdust

    i = npoints

    thisOctal => grid%octreeRoot
    position = startPoint
    direction = endPoint - startPoint
    call normalize(direction)
    if (myHydroSetGlobal /= 0) goto 555

    if (myrankWorldGlobal == 0) then

       open(20, file=thisFile, form="formatted", status="replace")
!       open(21, file="hc_components.dat", form="formatted", status="replace")
       do while(inOctal(grid%octreeRoot, position))
          call findSubcellLocal(position, thisOctal, subcell)
          sendThread = thisOctal%mpiThread(subcell)
          loc(1) = position%x
          loc(2) = position%y
          loc(3) = position%z
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
          call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)

          cen%x = tempStorage(1)
          cen%y = tempStorage(2)
          cen%z = tempStorage(3)
          AVtemp = tempStorage(4)
          tGas = tempStorage(5)
          tDust = tempStorage(6)
          H = tempStorage(7)
          H2 = tempStorage(8)
          cplus = tempstorage(9)
          c = tempstorage(10)
          c12o = tempstorage(11)
          CII21 = tempstorage(12)
          CI21 = tempstorage(13)
          CI31 = tempstorage(14)
          CI32 = tempstorage(15)
          OI21 = tempstorage(16)
          OI31 = tempstorage(17)
          OI32 = tempstorage(18)
          CO21 = tempstorage(19)
          CO32 = tempstorage(20)
          CO43 = tempstorage(21)
          CO54 = tempstorage(22)
          CO65 = tempstorage(23)
          CO76 = tempstorage(24)
          CO87 = tempstorage(25)
          tval = tempstorage(26)
          write(20,'(1p,40e14.5)') cen%x, AVtemp, tGas, tDust, H, H2, Cplus, C, C12O, CII21, CI21, CI31, &
               CI32, OI21, OI31, OI32, CO21, CO32, CO43, CO54, CO65, CO76, CO87

          position = cen
          position = position + (tVal+1.d-3*grid%halfSmallestSubcell)*direction

       enddo
       !Send escape trigger to other threads
       do sendThread = 1, nHydroThreadsGlobal
          loc(1) = 1.d30
          loc(2) = 1.d30
          loc(3) = 1.d30
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
       enddo
       close(20)
!       goto 555

    else
       stillLooping = .true.
       do while(stillLooping)
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
          position%x = loc(1)
          position%y = loc(2)
          position%z = loc(3)
          if (position%x > 1.d29) then
             stillLooping = .false.
          else
             call findSubcellLocal(position, thisOctal, subcell)
             sOctal => thisOctal
             cen = subcellCentre(thisOctal, subcell)
             call distanceToCellBoundary(grid, cen, direction, tVal, sOctal)
             tempStorage(1) = cen%x
             tempStorage(2) = cen%y
             tempStorage(3) = cen%z
             tempStorage(4) = thisOctal%AV(subcell, 1)
             tempStorage(5) = thisOctaL%tLast(subcell)
             tempStorage(6) = thisOctal%dust_t(subcell)
             tempstorage(7) = thisOctal%abundance(subcell, 31) !H2
             tempstorage(8) = thisOctal%abundance(subcell, 32) !H
             tempStorage(9) = thisOctal%abundance(subcell, 11) !C+ 
             tempStorage(10) = thisOctal%abundance(subcell, 25) !C
             tempStorage(11) = thisOctal%abundance(subcell, 28) !CO
             tempStorage(12) = thisOctal%ciiline(subcell, 2, 1) !cii(2-->1) emissivity
             tempStorage(13) = thisOctal%ciline(subcell, 2, 1) !ci(2-->1) emissivity
             tempStorage(14) = thisOctal%ciline(subcell, 3, 1) !ci(3-->1) emissivity
             tempStorage(15) = thisOctal%ciline(subcell, 3, 2) !ci(3-->2) emissivity
             tempStorage(16) = thisOctal%oiline(subcell, 2, 1) !oi(2-->1) emissivity
             tempStorage(17) = thisOctal%oiline(subcell, 3, 1) !oi(3-->1) emissivity
             tempStorage(18) = thisOctal%oiline(subcell, 3, 2) !oi(3-->2) emissivity
             tempStorage(19) = thisOctal%c12oline(subcell, 2, 1) !co(2-->1) emissivity
             tempStorage(20) = thisOctal%c12oline(subcell, 3, 2) !co(3-->2) emissivity
             tempStorage(21) = thisOctal%c12oline(subcell, 4, 3) !co(4-->3) emissivity
             tempStorage(22) = thisOctal%c12oline(subcell, 5, 4) !co(5-->4) emissivity
             tempStorage(23) = thisOctal%c12oline(subcell, 6, 5) !co(6-->5) emissivity
             tempStorage(24) = thisOctal%c12oline(subcell, 7, 6) !co(7-->6) emissivity
             tempStorage(25) = thisOctal%c12oline(subcell, 8, 7) !co(8-->7) emissivity
             tempstorage(26) = tval
             
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
             
          Endif
       enddo
    endif

555 continue
  end subroutine dumpFinalResultsAloneLine

  subroutine dumpHeatingCooling(grid, thisFile, startPoint, endPoint, nPoints, iter)
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell, iter
    integer :: nPoints
    type(VECTOR) :: startPoint, endPoint, position, direction, cen
    real(double) :: loc(3)!, rho, rhou , rhoe, p, phi_stars, phi_gas
!    real(double) :: temperature
    character(len=*) :: thisFile
    integer :: ierr
    integer, parameter :: nStorage = 20
!    real(double) :: c, cplus, c12o
    real(double) :: tempSTorage(nStorage),  heating, cooling
    real(double) :: tlast,  AVtemp, tval, tlow, thigh, heatingRate(1:12)
    integer, parameter :: tag = 30
    integer :: status(MPI_STATUS_SIZE)
    logical :: stillLooping
    integer :: sendThread
    integer :: i
    logical, save :: hcfile = .true.
    i = npoints

    thisOctal => grid%octreeRoot
    position = startPoint
    direction = endPoint - startPoint
    call normalize(direction)
    if (myHydroSetGlobal /= 0) goto 555

    if (myrankWorldGlobal == 0) then
       if(hcfile) then
          open(20, file=thisFile, form="formatted", status="replace")
          hcfile = .false.
       else
          open(20, file=thisFile, form="formatted", status="old", position="append")
       endif
       !       do while(inOctal(grid%octreeRoot, position))
       call findSubcellLocal(position, thisOctal, subcell)
       sendThread = thisOctal%mpiThread(subcell)
       loc(1) = position%x
       loc(2) = position%y
       loc(3) = position%z
       call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
       call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
       heatingRate = 0.d0
       cen%x = tempStorage(1)
       cen%y = tempStorage(2)
       cen%z = tempStorage(3)
       AVtemp = tempStorage(4)
       tlast = tempStorage(5)
       thigh = tempStorage(6)
       tlow = tempStorage(7)
       heating = tempStorage(8)
       cooling = tempStorage(20)
!       tval =  tempStorage(10)
!       c = tempstorage(11)
!       cplus = tempstorage(12)
!!       c12o = tempstorage(13)
       heatingRate(1:12) = tempstorage(8:19)
       write(20,'(1p,30e14.5)') dble(iter), heatingRate(12), cooling, tLast, thigh, tlow, heatingRate(:)
       !          write(20,'(1p,10e14.5)') cen%x, AVtemp, tlast, dustT, UVtemp, heating, cooling, c, cplus, c12o
       position = cen
!       position = position + (tVal+1.d-3*grid%halfSmallestSubcell)*direction

!    enddo
    !Send escape trigger to other threads
    do sendThread = 1, nHydroThreadsGlobal
       loc(1) = 1.d30
       loc(2) = 1.d30
       loc(3) = 1.d30
       call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
    enddo
    close(20)
    !       goto 555
    
 else
       stillLooping = .true.
       do while(stillLooping)
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
          position%x = loc(1)
          position%y = loc(2)
          position%z = loc(3)
          if (position%x > 1.d29) then
             stillLooping = .false.
          else
             call findSubcellLocal(position, thisOctal, subcell)
             sOctal => thisOctal
             cen = subcellCentre(thisOctal, subcell)
             call distanceToCellBoundary(grid, cen, direction, tVal, sOctal)
             tempStorage(1) = cen%x
             tempStorage(2) = cen%y
             tempStorage(3) = cen%z
             tempStorage(4) = thisOctal%AV(subcell, 1)
             tempStorage(5) = thisOctaL%tLast(subcell)
             tempStorage(6) = thisOctal%thigh(subcell)
             tempStorage(7) = thisOctal%tlow(subcell)             
             tempStorage(8) = thisOctal%heatingRate(subcell, 1)
             tempStorage(9) = thisOctal%heatingRate(subcell, 2)
             tempStorage(10) = thisOctal%heatingRate(subcell, 3)
             tempStorage(11) = thisOctal%heatingRate(subcell, 4)
             tempStorage(12) = thisOctal%heatingRate(subcell, 5)
             tempStorage(13) = thisOctal%heatingRate(subcell, 6)
             tempStorage(14) = thisOctal%heatingRate(subcell, 7)
             tempStorage(15) = thisOctal%heatingRate(subcell, 8)
             tempStorage(16) = thisOctal%heatingRate(subcell, 9)
             tempStorage(17) = thisOctal%heatingRate(subcell, 10)
             tempStorage(18) = thisOctal%heatingRate(subcell, 11)
             tempStorage(19) = thisOctal%heatingRate(subcell, 12)
!            tempStorage(8) = thisOctal%heatingRate(subcell, 1)
             tempStorage(20) = sum(thisOctal%coolingRate(subcell,:))
!             tempStorage(10) = tval
!             tempStorage(11) = thisOctal%abundance(subcell, 11) !C+ 
!             tempStorage(12) = thisOctal%abundance(subcell, 25) !C
!             tempStorage(13) = thisOctal%abundance(subcell, 28) !CO

!             tempStorage(9) = thisOctal%ci_abund(subcell)             
!             tempStorage(10) = thisOctal%oi_abund(subcell)             
!             tempStorage(11) = thisOctal%c12o_abund(subcell)
             
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
             
          endif
       enddo
    endif

555 continue
  end subroutine dumpHeatingCooling


#endif

  subroutine dumpDensitySpectrumZero(filename, filenameB, time)
    use mpi 
    use inputs_mod, only : normFac
    implicit none
    real(double) :: rho
    character(len=*) :: filename, filenameB
    integer :: ierr
    integer, parameter :: nBins = 14
    real(double) :: tempSTorage(3), mass, logrho, totalmass
    integer, parameter :: tag = 50
    integer :: status(MPI_STATUS_SIZE)
    logical :: stillRecving!, found
    integer :: i    
    real(double) :: dv, binIDs(nBins+1), massSpec(nBins), dRho
    real(double) :: ionizedmass, neutralmass, time
    logical, save :: firsttime = .true.

    if (myHydroSetGlobal /= 0) goto 333

    if (myrankWorldGlobal == 0) then
       if(firstTime) then
          !Overwrite any existing file
          open(23, file=filenameB, form="formatted", status="unknown")
          write(23,*) " "
          close(23)
          firstTime = .false.
       end if
    end if
   ! print *, "ENTERED"
    if(myrankglobal == 0) then
       open(22, file=filename, form="formatted", status="replace")
!       open(23, file=filenameB, form="formatted", status="replace")
       open(23, file=filenameB, form="formatted", status="unknown", position="append")
       binIDs(1) = -25.
       dRho = 0.5
       do i = 2, nBins+1
          binIDs(i) = binIDs(i-1) + dRho
       end do
  !     print *, "BIN IDS SET UP"
       massSpec = 0.d0
       ionizedmass = 0.d0
       neutralmass = 0.d0
       stillRecving=.true.
       do while (stillRecving)
!          print *, "WAITING TO RECV"
          call MPI_RECV(tempStorage, 3, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
               tag, localWorldCommunicator, status, ierr)       
 !         print *, "RECVD ONE"
          rho = tempStorage(1)
          dv = tempStorage(2)
          mass = rho*dv
          if (rho > 1.d20) then
             stillRecving = .false.
          end if
          if(stillRecving) then
             logrho = log10(rho)
             
!             Found = .false.
             

 !            do while (.not. found)
             do i=1, nBins
                if(logrho > binIDs(i) .and. logrho < binIDs(i+1)) then
                   massSpec(i) = massSpec(i) + mass                   
                   if(tempstorage(3) < 0.5d0) then
                      ionizedmass = ionizedmass + mass
                   else
                      neutralmass = neutralmass + mass
                   end if
                   !                      found = .true.
                end if
             end do
             !              found = .true.
             !          end do
          end if
       end do
       totalmass = sum(massSpec)
       print *, "totalmass is ", totalmass
       if(normfac /= 1.d0) then 
          massspec = massspec/normfac
       else
          massspec = massspec/totalmass
       end if

       print *, "DONE, preparing to dump data"
       print *, binIDs
       print *, massSpec
       do i = 1, nBins
          write(22,'(1p,7e14.2)') binIDs(i), massSpec(i)
          write(22,'(1p,7e14.2)') binIDs(i+1), massSpec(i)
       end do
       write(23,'(1p,7e14.3)') time, ionizedmass, neutralmass
       close(23)
       close(22)
   end if
   
   333 continue

 end subroutine dumpDensitySpectrumZero


 recursive subroutine dumpDensitySpectrumOther(thisOctal)
    use mpi 
!    use inputs_mod, only : griddistancescale
    implicit none
    type(OCTAL), pointer :: thisOctal, childOctal
    integer :: subcell!, childSubcell
    integer :: ierr
    integer, parameter :: nBins = 10
    real(double) :: tempSTorage(3)!, dv
    integer, parameter :: tag = 50
    integer :: i    
!    real(double) :: binIDs(nBins+1), massSpec(nBins)

    if (myHydroSetGlobal /= 0) goto 444

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                childOctal => thisOctal%child(i)
                call  dumpDensitySpectrumOther(childOctal)
                exit
             end if
          end do
       else 
          if(.not. thisoctal%ghostcell(subcell)) then
             tempStorage(1) = thisOctal%rho(subcell)
!             if(thisOctal%twoD) then
!                tempStorage(2) = (thisOctal%subcellSize**2)*griddistancescale**2
!             else.
             tempstorage(2) = cellVolume(thisOctal, subcell)*1.d30
             tempstorage(3) = thisOctal%ionfrac(subcell, 1)
!             end if
!             print *, "SENDING"
             call MPI_SEND(tempStorage, 3, MPI_DOUBLE_PRECISION, 0, &
                  tag, localWorldCommunicator, ierr) 
!             print *, "SENT"
          end if
       end if
    end do

444 continue

!    print *, "DONE"
 end subroutine dumpDensitySpectrumOther

 subroutine killDensitySpectrumDumper()
   use mpi
   implicit none
   integer, parameter :: tag = 50
   integer :: ierr
   real(double) :: tempSTorage(3)

   if (myHydroSetGlobal /= 0) goto 555
   tempStorage(1) = 1.d25
   tempStorage(2) = 0.d0
   tempStorage(3) = 0.d0
   call MPI_SEND(tempStorage, 3, MPI_DOUBLE_PRECISION, 0, &
        tag, localWorldCommunicator, ierr)   
555 continue
 end subroutine killDensitySpectrumDumper

  subroutine tauRadius(grid, rVec, uHat, tauWanted, tauRad)
    use mpi
    use inputs_mod, only : smallestCellsize
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    type(VECTOR) :: uHat, position, rVec
    real(double) :: loc(4), tauRad, tauWanted, tau, flag
    integer, parameter :: nStorage = 5
    real(double) :: tempSTorage(nStorage), tval
    integer, parameter :: tag = 30
    integer :: status(MPI_STATUS_SIZE), ierr
    real(double) :: kappap
    logical :: stillLooping
    integer :: sendThread
    logical :: hitGrid

    if (loadBalancingThreadGlobal) goto 666

    thisOctal => grid%octreeRoot
    position = rVec + ((-10.d0*grid%octreeRoot%subcellSize) * uHat)
    tval = distanceToGridFromOutside(grid, position, uHat, hitGrid)
    if (hitGrid) then
       position = position + (tval+0.01d0*smallestcellsize)*uHat
    else
       stop
    endif
    tau = 0.d0
    if (myrankGlobal == 0) then

       do while(inOctal(grid%octreeRoot, position))
          call findSubcellLocal(position, thisOctal, subcell)
          sendThread = thisOctal%mpiThread(subcell)
          loc(1) = position%x
          loc(2) = position%y
          loc(3) = position%z
          loc(4) = tau
          call MPI_SEND(loc, 4, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
          call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
          flag = tempStorage(1)
          if (tempStorage(1) > 1.d29) then
             tauRad = tempStorage(2)
             exit
          else if (tempStorage(1) < -1.d29) then
             tauRad = 0.d0
             exit
          else
             tau = tau + tempStorage(2)
             position = VECTOR(tempStorage(3), tempStorage(4), tempStorage(5))
          endif
         
       enddo
       !Send escape trigger to other threads
       do sendThread = 1, nHydroThreadsGlobal
          loc(1) = 1.d30
          loc(2) = 1.d30
          loc(3) = 1.d30
          call MPI_SEND(loc, 5, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
       enddo
       close(20)
       goto 666


    else
       stillLooping = .true.
       do while(stillLooping)
          call MPI_RECV(loc, 5, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
          position%x = loc(1)
          position%y = loc(2)
          position%z = loc(3)
          if (position%x > 1.d29) then
             stillLooping = .false.
          else

             tau = loc(4)

             do
                call findSubcellLocal(position, thisOctal, subcell)
                if (.not.OctalOnThread(thisOctal, subcell, myRankGlobal)) then
                   tempStorage(1) = 0.d0
                   tempStorage(2) = tau
                   tempStorage(3) = position%x
                   tempStorage(4) = position%y
                   tempStorage(5) = position%z
                   exit
                endif


                sOctal => thisOctal
                call distanceToCellBoundary(grid, position, uHat, tVal, sOctal)
                call returnKappa(grid, thisOctal, subcell, kappap = kappap)
                tau = tau + dble(kappap) * tval * 1.d10 
                position = position + (tVal+0.01d0*thisOctal%subcellSize)*uHat
                if (tau > tauWanted) then
                   tauRad = modulus(position-rVec)
                   tempStorage(1) = 1.d30
                   tempStorage(2) = tauRad 
                   exit
                endif

                if (.not.inOctal(grid%octreeRoot, position)) then
                   tempStorage(1) = -1.d30
                   exit
                endif
             enddo
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          endif
       enddo
    endif

666 continue
  end subroutine tauRadius

  subroutine tauAlongPathMPI(grid, rVec, uHat, tauAbs, tauSca)
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    type(VECTOR) :: uHat, position, rVec
    real(double) :: loc(5), tauAbs, tauSca, flag
    integer, parameter :: nStorage = 6
    real(double) :: tempSTorage(nStorage), tval
    integer, parameter :: tag = 30
    integer :: status(MPI_STATUS_SIZE), ierr
    real(double) :: kappaAbs, kappaSca
    logical :: stillLooping
    integer :: sendThread

    if (loadBalancingThreadGlobal) goto 666
    thisOctal => grid%octreeRoot
    position = rVec
    tauAbs = 0.d0; tauSca = 0.d0;
    if (myrankGlobal == 0) then

       do while(inOctal(grid%octreeRoot, position))
          call findSubcellLocal(position, thisOctal, subcell)
          sendThread = thisOctal%mpiThread(subcell)
          loc(1) = position%x
          loc(2) = position%y
          loc(3) = position%z
          loc(4) = tauAbs
          loc(5) = tauSca
          call MPI_SEND(loc, 5, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
          call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
          flag = tempStorage(1)
          if (tempStorage(1) > 1.d29) then
             tauAbs = tempStorage(2)
             tauSca = tempStorage(3)
             exit
          else if (tempStorage(1) < -1.d29) then
             tauAbs = tempStorage(2)
             tauSca = tempStorage(3)
             exit
          else
             tauAbs = tauAbs + tempStorage(2)
             tauSca = tauSca + tempStorage(3)
             position = VECTOR(tempStorage(4), tempStorage(5), tempStorage(6))
          endif
         
       enddo
       !Send escape trigger to other threads
       do sendThread = 1, nHydroThreadsGlobal
          loc(1) = 1.d30
          loc(2) = 1.d30
          loc(3) = 1.d30
          call MPI_SEND(loc, 5, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
       enddo
       goto 666


    else
       stillLooping = .true.
       do while(stillLooping)
          call MPI_RECV(loc, 5, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
          position%x = loc(1)
          position%y = loc(2)
          position%z = loc(3)
          if (position%x > 1.d29) then
             stillLooping = .false.
          else

             tauAbs = loc(4)
             tauSca = loc(5)

             do
                call findSubcellLocal(position, thisOctal, subcell)
                if (.not.OctalOnThread(thisOctal, subcell, myRankGlobal)) then
                   tempStorage(1) = 0.d0
                   tempStorage(2) = tauAbs
                   tempStorage(3) = tauSca
                   tempStorage(4) = position%x
                   tempStorage(5) = position%y
                   tempStorage(6) = position%z
                   exit
                endif

                sOctal => thisOctal
                call distanceToCellBoundary(grid, position, uHat, tVal, sOctal)
                call returnKappa(grid, thisOctal, subcell, ilambda=1,&
                     kappaSca=kappaSca, kappaAbs=kappaAbs)
                tauAbs = tauAbs + dble(kappaAbs) * tval
                tauSca = tauSca + dble(kappaSca) * tval
                position = position + (tVal+0.01d0*thisOctal%subcellSize)*uHat

                if (.not.inOctal(grid%octreeRoot, position)) then
                   tempStorage(1) = -1.d30
                   tempStorage(2) = tauAbs
                   tempStorage(3) = tauSca
                   tempStorage(4) = position%x
                   tempStorage(5) = position%y
                   tempStorage(6) = position%z
                   exit
                endif
             enddo
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          endif
       enddo
    endif

666 continue
  end subroutine tauAlongPathMPI
    
  subroutine grid_info_mpi(thisGrid, filename)
    use amr_mod, only:  countVoxels
    use mpi
    type(gridtype), intent(in) :: thisGrid
    character(LEN=*), intent(in) :: filename
    integer :: UN
    real(double) :: fac
    integer :: nOctals,nVoxels
    integer :: tempInt(1)
    real(double) :: tempDouble(1)
    real(double) :: halfSmallestSubcell
    integer :: maxDepth
    integer :: ierr
    
    if (myrankWorldGlobal == 1) then
       if (filename(1:1) == '*') then
          UN = 6   ! prints on screen
       else
          UN = 69
          open(unit=UN, file = TRIM(filename), status = 'unknown',form='formatted')
       end if
    endif

    nOctals=0; nVoxels=0
    if (myRankGlobal /= 0) then
       call countSubcellsMPI(thisgrid, nVoxels)
       call MPI_REDUCE(thisgrid%maxDepth, tempInt, 1, MPI_INTEGER, MPI_MAX, 1, amrCommunicator, ierr)
       maxDepth = tempInt(1)
       call MPI_REDUCE(thisgrid%halfSmallestSubcell, tempDouble,1,MPI_DOUBLE_PRECISION, MPI_MIN, 1, amrCommunicator, ierr)
       halfSmallestSubcell = tempDouble(1)
    endif
    if (myRankWorldGlobal == 1) then
       write(UN,'(a)') ' '
       write(UN,'(a)') '######################################################'
       write(UN,'(a)') 'Grid info :'
       write(UN,'(a)') ' '
       write(UN,*)     'geometry             = ', thisGrid%geometry
       write(UN,*)     'maxDepth             = ', maxDepth
       write(UN,*)     'halfSmallestSubcell  = ', halfSmallestSubcell, ' [10^10 cm]'
       write(UN,*)     'nVoxels              = ', nVoxels
       write(UN,*)     'smoothingFactor      = ', thisGrid%smoothingFactor
       write(UN,*)     'grid center          =', thisGrid%octreeRoot%centre
       write(UN,*)     'Size of largest cell =', thisGrid%octreeRoot%subcellSize*2.0, ' [10^10 cm]'
       write(UN,'(a)') '#######################################################'
       write(UN,'(a)') ' '
       write(Un,'(a)') ' '
       
       fac = 1.d0/(2.d0**dble(maxDepth))
       if (fac < epsilon(1.d0)) then
          write(UN,'(a)') "**** WARNING: Grid cell depth is so great numerical problems may occur****"
       endif
    endif
    if (myrankWorldGlobal == 1) then
       if (filename(1:1) /= '*')  close(UN)
    endif
    
  end subroutine grid_info_mpi
  
#ifdef HYDRO
  subroutine addNewChildWithInterp(parent, iChild, grid, constantGravity)
    use inputs_mod, only : maxDepthAMR!, minDepthAmr
    use octal_mod, only: subcellRadius
    use utils_mod
    use mpi
    use qshep3d_mod
    use qshep2d_mod
    type(OCTAL), pointer :: parent, thisOctal, topOctal
    integer :: iChild
    logical, optional :: constantGravity
    type(GRIDTYPE) :: grid
    integer :: nChildren
    integer :: newChildIndex
    integer :: i
    integer :: iSubcell, parentSubcell, topOctalSubcell
    type(VECTOR) :: centre !, rVec
    real(double) :: massFactor, newMom, oldMom
    real(double) :: x1, x2, y1, y2, z1, z2, dv !, u, x, y, z
    real(double) :: oldMass, newMass, factor
    real(double) :: oldEnergy, newEnergy
!    integer :: npoints
!    integer :: ier
!    integer :: nr, nw, nq
!    integer, allocatable :: lnext(:), lcell(:,:,:)
!    integer, allocatable :: lcell2d(:,:)
!    real(double) :: xyzmin(3), xyzdel(3)
!    real(double) :: rMax
!    real(double), allocatable :: rsq(:), a(:,:)
!    integer, parameter :: maxpts = 10000
!    real(double) :: xPoint(maxpts)
!    real(double) :: yPoint(maxpts)
!    real(double) :: zPoint(maxpts)
!    real(double) :: rhoPoint(maxpts)
!    real(double) :: rhoePoint(maxpts)
!    real(double) :: uPoint(maxpts)
!    real(double) :: vPoint(maxpts)
!    real(double) :: wPoint(maxpts)
!    real(double) :: rhouPoint(maxpts)
!    real(double) :: rhovPoint(maxpts)
!    real(double) :: rhowPoint(maxpts)
!    real(double) :: phiPoint(maxpts)
!    real(double) :: energyPoint(maxpts)
!    real(double) :: pressurePoint(maxpts)
!    real(double) :: xmin, zmin, dx, dz
!    real(double) :: thisRho
!    integer :: counter
!    character(len=80) :: message
!    real(double) :: radius
    logical, save :: firstTime = .true.
    logical :: debug !, successful , doLogspace, triedLogSpace

    debug = .false.


    if(parent%edgecell(iChild)) then
       call addnewchild(parent,iChild, grid, adjustGridInfo=.true.)
    else


!    if (inSubcell(parent, ichild, testVec)) debug = .true.
    if (parent%ndepth == maxDepthAMR) then
       if (firstTime) then
          call writeWarning("Cell depth capped in addNewChildWithInterp")
          firstTime = .false.
       endif
       goto 666
    endif

    ! store the number of children that already exist
    nChildren = parent%nChildren

    ! safety checks of child array
    IF ( ASSOCIATED(parent%child) ) THEN
      IF ( ( nChildren == 0 ) .OR.                  &
           ( nChildren /= SIZE(parent%child) ) ) THEN
        PRINT *, 'Panic: in addNewChild, %child array wrong size'
        PRINT *, 'nChildren:',nChildren,' SIZE %child:', SIZE(parent%child)
        STOP
      END IF
    END IF
    IF ( (.NOT. ASSOCIATED(parent%child)) .AND. (nChildren > 0) ) THEN
      PRINT *, 'Panic: in addNewChild, %child array wrong size'
      PRINT *, 'nChildren:',nChildren,' ASSOCIATED %child:', ASSOCIATED(parent%child)
      STOP
    END IF

    ! check that new child does not already exist
    IF ( parent%hasChild(iChild) .EQV. .TRUE. ) THEN
      PRINT *, 'Panic: in addNewChild, attempted to add a child ',&
               '       that already exists'
      STOP
    ENDIF

    CALL growChildArray(parent, nNewChildren=1, grid=grid )

    ! update the bookkeeping
    nChildren = nChildren + 1
    newChildIndex = nChildren

    ! update the parent octal
    parent%nChildren = nChildren
    parent%hasChild(iChild) = .TRUE.
    parent%indexChild(newChildIndex) = iChild

    ! allocate any variables that need to be  
    IF (.NOT.grid%oneKappa) THEN
       ! The kappa arrays should be allocated with grid%nopacity instead of grid%nlambda
       ! because for line calculation, there is only one kappa needed.
       ! (but grid%nlambda is not 1). If you allocate the arrays with grid%nlambda,
       ! it will be a huge waste of RAM. ---  (RK) 
       ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nopacity))
       ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nopacity))
       ! ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nlambda))
       ! ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nlambda))
       parent%child(newChildIndex)%kappaAbs = 1.e-30
       parent%child(newChildIndex)%kappaSca = 1.e-30
    ENDIF
    NULLIFY(parent%child(newChildIndex)%child)


    parent%child(newChildIndex)%nDepth = parent%nDepth + 1

! setup mpiThread values

    if ( ((parent%twoD)  .and.((nHydroThreadsGlobal) == 4)) .or. &
         ((parent%threed).and.((nHydroThreadsGlobal) == 8)).or. &
         ((parent%oneD)  .and.((nHydroThreadsGlobal) == 2)) ) then
       parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
    else

       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          if (parent%oneD) then
             do i = 1, 2
                parent%child(newChildIndex)%mpiThread(i) = 2 * (parent%mpiThread(iChild) - 1) + i
             enddo
          else if (parent%twoD) then
             do i = 1, 4
                parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
             enddo
          else if (parent%Threed) then
             do i = 1, 8
                parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
             enddo
          endif
       endif
    endif

    if ((parent%threed).and.(nHydroThreadsGlobal) == 64) then
       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 8
             parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    if ((parent%threed).and.(nHydroThreadsGlobal) == 512) then
       if (parent%child(newChildIndex)%nDepth > 3) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 8
             parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    if ((parent%twod).and.(nHydroThreadsGlobal) == 64) then
       if (parent%child(newChildIndex)%nDepth > 3) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 4
             parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    ! set up the new child's variables
    parent%child(newChildIndex)%threeD = parent%threeD
    parent%child(newChildIndex)%twoD = parent%twoD
    parent%child(newChildIndex)%oneD = parent%oneD
    parent%child(newChildIndex)%maxChildren = parent%maxChildren
    parent%child(newChildIndex)%cylindrical = parent%cylindrical

       


    parent%child(newChildIndex)%inFlow = parent%inFlow
    parent%child(newChildIndex)%parent => parent
    parent%child(newChildIndex)%parentSubcell = iChild
    parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
    parent%child(newChildIndex)%hasChild = .false.
    parent%child(newChildIndex)%nChildren = 0
    parent%child(newChildIndex)%indexChild = -999 ! values are undefined
    parent%child(newChildIndex)%centre = subcellCentre(parent,iChild)
    if (parent%cylindrical) then
       parent%child(newChildIndex)%r = subcellRadius(parent,iChild)
    endif

    parent%child(newChildIndex)%xMin = parent%child(newChildIndex)%centre%x - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMin = parent%child(newChildIndex)%centre%y - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMin = parent%child(newChildIndex)%centre%z - parent%child(newChildIndex)%subcellSize

    parent%child(newChildIndex)%xMax = parent%child(newChildIndex)%centre%x + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMax = parent%child(newChildIndex)%centre%y + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMax = parent%child(newChildIndex)%centre%z + parent%child(newChildIndex)%subcellSize

    parentSubcell = iChild

    thisOctal => parent%child(newChildIndex)
    call allocateOctalAttributes(grid, thisOctal)

    thisOctal%boundaryCondition = parent%boundaryCondition(parentSubcell)
    thisOctal%temperature = parent%temperature(parentSubcell)
    thisOctal%iEquationOfState = parent%iEquationofState(parentSubcell)
    thisOctal%gamma = parent%gamma(parentSubcell)
    thisOctal%label = 666
    thisOctal%changed = .true.


    if (associated(parent%ionFrac)) then
       do iSubcell = 1, thisOctal%maxChildren
          thisOctal%ionFrac(isubcell,:) = parent%ionFrac(parentSubcell,:)
       enddo
    endif

    if (associated(parent%ne)) then
       do iSubcell = 1, thisOctal%maxChildren
          thisOctal%ne(isubcell) = parent%ne(parentSubcell)
       enddo
    endif



    topOctal => thisOctal%parent
    topOctalSubcell = thisOctal%parentsubcell
    do while(topOctal%changed(topOctalSubcell))
       topOctalSubcell = topOctal%parentSubcell
       topOctal => topOctal%parent
    enddo

    x1 = centre%x - topOctal%subcellSize/2.d0
    x2 = centre%x + topOctal%subcellSize/2.d0
    y1 = centre%y - topOctal%subcellSize/2.d0
    y2 = centre%y + topOctal%subcellSize/2.d0
    z1 = centre%z - topOctal%subcellSize/2.d0
    z2 = centre%z + topOctal%subcellSize/2.d0

    do iSubcell = 1, thisOctal%maxChildren

       if (thisOctal%threed) then
!          rVec = subcellcentre(thisOctal, iSubcell)
!          x = rVec%x
!          y = rVec%y
!          z = rVec%z


!          call myInterp(2,nPoints, xPoint, yPoint, zPoint, energyPoint, x, y, z, thisOctal%energy(iSubcell))
!          call myInterp(2,nPoints, xPoint, yPoint, zPoint, pressurePoint, x, y, z, thisOctal%pressure_i(iSubcell))
!          call myInterp(2,nPoints, xPoint, yPoint, zPoint, phiPoint, x, y, z, thisOctal%phi_gas(iSubcell))


         call interpTrilinear(grid, thisOctal, isubcell, thisOctal%rho(isubcell), thisOctal%rhoe(iSubcell), &
               thisOctal%rhou(iSubcell), thisOctal%rhov(iSubcell), thisOctal%rhow(iSubcell), thisOctal%energy(iSubcell), &
               thisOctal%phi_gas(iSubcell), thisOctal%pressure_i(iSubcell))


!
!          do while (nPoints < 33)
!             call getPointsInRadius(rVec, radius, grid, npoints, rhoPoint, rhoePoint, &
!                  rhouPoint, rhovPoint, rhowPoint, energyPoint, pressurePoint, phiPoint, xPoint, yPoint, zPoint)
!             radius = radius * 2.d0
!          enddo
!
!
!
!          uPoint(1:nPoints) = rhouPoint(1:nPoints)/rhoPoint(1:nPoints)
!          vPoint(1:nPoints) = rhovPoint(1:nPoints)/rhoPoint(1:nPoints)
!          wPoint(1:nPoints) = rhowPoint(1:nPoints)/rhoPoint(1:nPoints)
!
!          nq = 17 !min(40, nPoints - 1)
!          nw = 32 !min(40, nPoints - 1)
!          nr = max(1,nint((dble(nPoints)/3.d0)**0.333d0))
!          allocate(lCell(1:nr,1:nr,1:nr))
!          allocate(lnext(1:nPoints))
!          allocate(rsq(1:nPoints))
!          allocate(a(9,1:nPoints))
!          triedLogSpace = .false.
!          doLogSpace = .false.
!          successful = .false.
!          do while(.not.successful)
!
!             if (doLogSpace) then
!                rhoPoint(1:nPoints) = log10(rhoPoint(1:nPoints))
!             endif
!
!             call qshep3 (nPoints, xPoint, yPoint, zPoint, rhoPoint, nq, nw, nr, lcell, lnext, xyzmin, &
!                  xyzdel, rmax, rsq, a, ier ) 
!             if (ier /= 0) then
!                write(message,*) "Qshep3 returned an error ",ier
!                call writeWarning(message)
!                write(*,*) " npoints ",npoints
!                do i = 1, nPoints
!                   write(*,*) xPoint(i), ypoint(i),zpoint(i)
!                enddo
!             endif
!
!             thisRho = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, rhoPoint, nr, lcell, lnext, &
!                  xyzmin, xyzdel, rmax, rsq, a)
!             if (doLogSpace) then
!                thisRho = 10.d0**thisRho
!                triedLogSpace = .true.
!             endif
!             successful = .true.
!
!             if (thisRho < 0.d0) then
!                successful = .false.
!                write(*,*) "Negative density encountered in qs3val: ",thisRho
!                if (.not.triedLogspace) then
!                   write(*,*) "Trying log space interpolation"
!                   write(77,'(3e12.4)') x,y,z
!                   dologSpace = .true.
!                   do i = 1, nPoints
!                      write(*,*) xPoint(i), ypoint(i),zpoint(i),rhoPoint(i)
!                      write(77,'(3e12.4)') xPoint(i), ypoint(i),zpoint(i)
!                   enddo
!                   stop
!                else
!                   write(*,*) "Can't find sensible interpolation."
!                   write(77,'(3e12.4)') x,y,z
!                   do i = 1, nPoints
!                      write(*,*) xPoint(i), ypoint(i),zpoint(i)
!                      write(77,'(3e12.4)') xPoint(i), ypoint(i),zpoint(i)
!                   enddo
!                   stop
!                endif
!             endif
!          enddo
!
!          
!          thisOctal%rho(isubcell) = thisRho
!
!          call qshep3 (nPoints, xPoint, yPoint, zPoint, rhoePoint, nq, nw, nr, lcell, lnext, xyzmin, &
!               xyzdel, rmax, rsq, a, ier )
!          if (ier /= 0) then
!             write(message,*) "Qshep3 returned an error ",ier
!             call writeWarning(message)
!          endif
!
!          if (dologSpace) rhoePoint(1:nPoints) = log10(rhoePoint(1:nPoints))
!          thisOctal%rhoe(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, rhoePoint, nr, lcell, lnext, &
!               xyzmin, xyzdel, rmax, rsq, a)
!          if (dologSpace) thisOctal%rhoe(isubcell) = 10.d0**thisOctal%rhoe(isubcell)
!
!          call qshep3 (nPoints, xPoint, yPoint, zPoint, uPoint, nq, nw, nr, lcell, lnext, xyzmin, &
!               xyzdel, rmax, rsq, a, ier )
!          if (ier /= 0) then
!             write(message,*) "Qshep3 returned an error ",ier
!             call writeWarning(message)
!          endif
!
!          thisOctal%rhou(iSubcell) = thisRho * qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, uPoint, nr, lcell, lnext, &
!               xyzmin, xyzdel, rmax, rsq, a)
!
!          call qshep3 (nPoints, xPoint, yPoint, zPoint, vPoint, nq, nw, nr, lcell, lnext, xyzmin, &
!               xyzdel, rmax, rsq, a, ier )
!          if (ier /= 0) then
!             write(message,*) "Qshep3 returned an error ",ier
!             call writeWarning(message)
!          endif
!
!          thisOctal%rhov(iSubcell) = thisRho * qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, vPoint, nr, lcell, lnext, &
!               xyzmin, xyzdel, rmax, rsq, a)
!
!          call qshep3 (nPoints, xPoint, yPoint, zPoint, wPoint, nq, nw, nr, lcell, lnext, xyzmin, &
!               xyzdel, rmax, rsq, a, ier )
!          if (ier /= 0) then
!             write(message,*) "Qshep3 returned an error ",ier
!             call writeWarning(message)
!          endif
!
!          thisOctal%rhow(iSubcell) = thisRho * qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, wPoint, nr, lcell, lnext, &
!               xyzmin, xyzdel, rmax, rsq, a)
!
!          if (dologSpace) energyPoint(1:nPoints) = log10(energyPoint(1:nPoints))
!          call qshep3 (nPoints, xPoint, yPoint, zPoint, energyPoint, nq, nw, nr, lcell, lnext, xyzmin, &
!               xyzdel, rmax, rsq, a, ier )
!          if (ier /= 0) then
!             write(message,*) "Qshep3 returned an error ",ier
!             call writeWarning(message)
!          endif
!
!          thisOctal%energy(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, energyPoint, nr, lcell, lnext, &
!               xyzmin, xyzdel, rmax, rsq, a)
!          if (dologSpace) thisOctal%energy(isubcell) = 10.d0**thisOctal%energy(isubcell)
!
!
!          call qshep3 (nPoints, xPoint, yPoint, zPoint, phiPoint, nq, nw, nr, lcell, lnext, xyzmin, &
!               xyzdel, rmax, rsq, a, ier )
!          if (ier /= 0) then
!             write(message,*) "Qshep3 returned an error ",ier
!             call writeWarning(message)
!          endif
!
!          thisOctal%phi_gas(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, phiPoint, nr, lcell, lnext, &
!               xyzmin, xyzdel, rmax, rsq, a)
!
!          if (doLogSpace) pressurePoint(1:nPoints) = log10(pressurepoint(1:nPoints))
!          call qshep3 (nPoints, xPoint, yPoint, zPoint, pressurePoint, nq, nw, nr, lcell, lnext, xyzmin, &
!               xyzdel, rmax, rsq, a, ier )
!          if (ier /= 0) then
!             write(message,*) "Qshep3 returned an error ",ier
!             call writeWarning(message)
!          endif
!
!          thisOctal%pressure_i(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, pressurePoint, nr, lcell, lnext, &
!               xyzmin, xyzdel, rmax, rsq, a)
!          if (doLogSpace) thisOctal%pressure_i(iSubcell) = 10.d0**thisOctal%pressure_i(isubcell)
!
!          deallocate(lCell, lnext, rsq, a)
!          
       endif
       if (thisOctal%twod) then

         call interpBilinear(grid, thisOctal, isubcell, thisOctal%rho(isubcell), thisOctal%rhoe(iSubcell), &
               thisOctal%rhou(iSubcell), thisOctal%rhov(iSubcell), thisOctal%rhow(iSubcell), thisOctal%energy(iSubcell), &
               thisOctal%phi_gas(iSubcell), thisOctal%pressure_i(iSubcell))

!          rVec = subcellcentre(thisOctal, iSubcell)
!          x = rVec%x
!          y = 0.d0
!          z = rVec%z
!!          nPoints = 0
!          rVec%y = 0.d0
!
!!          radius = 4.d0*grid%octreeRoot%subcellSize / &
!!                                2.0_oc**REAL(minDepthAmr,kind=oct)
!
!          radius = thisOctal%subcellSize*4.d0
!
!          do while (nPoints < 12)
! !            call returnAddNewChildPointArrays(thisOctal, topOctal, topOctalSubcell, rhoPoint, rhoePoint, rhouPoint, &
! !                 rhovPoint, rhowPoint, phiPoint, energyPoint, pressurePoint)
!             call getPointsInRadius(rVec, radius, grid, npoints, rhoPoint, rhoePoint, &
!                  rhouPoint, rhovPoint, rhowPoint, energyPoint, pressurePoint, phiPoint, xPoint, yPoint, zPoint)
!             radius = radius * 2.d0
!          end do
!
!          ypoint = 0.d0
!          y = 0.d0
!!          rhovPoint = 0.d0
!                    
!          nq = min(40, nPoints - 1)
!          nw = min(40, nPoints - 1)
!          nr = max(1,nint((dble(nPoints)/3.d0)**0.5d0))
!          allocate(lCell2d(1:nr,1:nr))
!          allocate(lnext(1:nPoints))
!          allocate(rsq(1:nPoints))
!          allocate(a(5,1:nPoints))
!
!          call qshep2 (nPoints, xPoint, zPoint, rhoPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
!               dx, dz, rmax, rsq, a, ier )
!          if (ier /= 0) call writeWarning("Qshep2 returned an error for rho")
!          thisOctal%rho(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhoPoint, nr, lcell2d, lnext, &
!               xmin, zmin, dx, dz, rmax, rsq, a)
!
!
!          call qshep2 (nPoints, xPoint, zPoint, rhoePoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
!               dx, dz, rmax, rsq, a, ier )
!
!          if (ier /= 0) call writeWarning("Qshep2 returned an error for rhoe")
!          thisOctal%rhoe(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhoePoint, nr, lcell2d, lnext, &
!               xmin, zmin, dx, dz, rmax, rsq, a)
!
!          call qshep2 (nPoints, xPoint, zPoint, rhouPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
!               dx, dz, rmax, rsq, a, ier )
!
!          if (ier /= 0) call writeWarning("Qshep2 returned an error for rhou")
!          thisOctal%rhou(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhouPoint, nr, lcell2d, lnext, &
!               xmin, zmin, dx, dz, rmax, rsq, a)
!
!
!          call qshep2 (nPoints, xPoint, zPoint, rhovPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
!               dx, dz, rmax, rsq, a, ier )
!
!          if (ier /= 0) call writeWarning("Qshep2 returned an error for rhou")
!          thisOctal%rhov(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhovPoint, nr, lcell2d, lnext, &
!               xmin, zmin, dx, dz, rmax, rsq, a)
!
!
!          call qshep2 (nPoints, xPoint, zPoint, rhowPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
!               dx, dz, rmax, rsq, a, ier )
!
!          if (ier /= 0) call writeWarning("Qshep2 returned an error for rhow")
!          thisOctal%rhow(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhowPoint, nr, lcell2d, lnext, &
!               xmin, zmin, dx, dz, rmax, rsq, a)
!
!          call qshep2 (nPoints, xPoint, zPoint, energyPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
!               dx, dz, rmax, rsq, a, ier )
!
!          if (ier /= 0) call writeWarning("Qshep2 returned an error for energy")
!          thisOctal%energy(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, energyPoint, nr, lcell2d, lnext, &
!               xmin, zmin, dx, dz, rmax, rsq, a)
!
!          call qshep2 (nPoints, xPoint, zPoint, phiPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
!               dx, dz, rmax, rsq, a, ier )
!
!          if (ier /= 0) call writeWarning("Qshep2 returned an error for phi")
!          thisOctal%phi_gas(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, phiPoint, nr, lcell2d, lnext, &
!               xmin, zmin, dx, dz, rmax, rsq, a)
!
!
!          call qshep2 (nPoints, xPoint, zPoint, pressurePoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
!               dx, dz, rmax, rsq, a, ier )
!
!          if (ier /= 0) call writeWarning("Qshep2 returned an error for pressure")
!          thisOctal%pressure_i(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, pressurePoint, nr, lcell2d, lnext, &
!               xmin, zmin, dx, dz, rmax, rsq, a)
!
!
!          deallocate(lCell2d, lnext, rsq, a)

      endif

       if (thisOctal%oned) then
         call interplinear(grid, thisOctal, isubcell, thisOctal%rho(isubcell), thisOctal%rhoe(iSubcell), &
               thisOctal%rhou(iSubcell), thisOctal%rhov(iSubcell), thisOctal%rhow(iSubcell), thisOctal%energy(iSubcell), &
               thisOctal%phi_gas(iSubcell), thisOctal%pressure_i(iSubcell))
       Endif
       
    Enddo
    
    ! conservation normalizations

    ! mass
       if (thisOctal%threed) then
!          dv = cellVolume(thisOctal, parentSubcell) * 1.d30
          dv = parent%subcellSize**3
       else if (thisOctal%twoD) then
          dv = parent%subcellSize**2
       else if (thisOctal%oneD) then
          dv = parent%subcellSize
       endif
       
       oldMass = parent%rho(parentSubcell) * dv !cellVolume(parent, parentSubcell)
       !    oldMass = parent%rho(parentSubcell) * cellVolume(parent, parentSubcell)
       newMass = 0.d0
       
       do iSubcell = 1, thisOctal%maxChildren
          !THAW - redoing mass calculation
          if (thisOctal%threed) then
!             dv = cellVolume(thisOctal, iSubcell) * 1.d30
             dv = thisOctal%subcellSize**3
          else if (thisOctal%twoD) then
             dv = thisOctal%subcellSize**2
          else if (thisOctal%oneD) then
             dv = thisOctal%subcellSize
          endif
          !       print *, "dv", dv, cellVolume(thisOctal, iSubcell)
          
          !       newMass = newMass + thisOctal%rho(isubcell) * cellVolume(thisOctal, iSubcell)
          newMass = newMass + thisOctal%rho(isubcell) * dv
       enddo
       
!       write(*,*) "ndepth ",thisOctal%nDepth
!       write(*,*) "rho parent ",parent%rho(parentsubcell)
!       write(*,*) "old mass, new mass ",oldmass, newmass
       massfactor = oldMass / newMass
       thisOctal%rho(1:thisOctal%maxChildren) = thisOctal%rho(1:thisOctal%maxChildren) * massfactor
       
    if ( associated (thisOctal%nh) ) thisOctal%nh(1:thisOctal%maxChildren) = thisOctal%rho(1:thisOctal%maxChildren)/mHydrogen

    ! energy

    if (thisOctal%threed) then
!       dv = cellVolume(thisOctal, parentSubcell) * 1.d30
       dv = parent%subcellSize**3
    else if (thisOctal%twoD) then
       dv = parent%subcellSize**2
    else if (thisOctal%oneD) then
       dv = parent%subcellSize
    endif

    oldEnergy = parent%rhoe(parentSubcell) * dv !cellVolume(parent, parentSubcell)
!    oldEnergy = parent%rhoe(parentSubcell) * cellVolume(parent, parentSubcell)

    newEnergy = 0.d0

    do iSubcell = 1, thisOctal%maxChildren
!       newEnergy = newEnergy + thisOctal%rhoe(isubcell) * cellVolume(thisOctal, iSubcell)
       !THAW - redoing energy calculation
       if (thisOctal%threed) then
!          dv = cellVolume(thisOctal, iSubcell) * 1.d30
          dv = thisOctal%subcellSize**3
       else if (thisOctal%twoD) then
          dv = thisOctal%subcellSize**2
       else if (thisOctal%oneD) then!
          dv = thisOctal%subcellSize
       endif
        newEnergy = newEnergy + thisOctal%rhoe(isubcell) * dv
    enddo
    factor = 1.d0
    if (newEnergy /= 0.d0) then
       factor = oldEnergy/newEnergy
    endif
    thisOctal%rhoe(1:thisOctal%maxChildren) = thisOctal%rhoe(1:thisOctal%maxChildren) * factor

!
! momentum (u)
!
    if ( all(thisOctal%rhou(1:thisOctal%maxChildren) < 0.d0) .or. &
         all(thisOctal%rhou(1:thisOctal%maxChildren) > 0.d0) ) then
       oldMom = parent%rhou(parentSubcell) * cellVolume(parent,parentSubcell) 
       newMom = SUM(thisOctal%rhou(1:thisOctal%maxChildren)*cellVolume(thisOctal,1))
       factor = oldMom / newMom
       thisOctal%rhou(1:thisOctal%maxChildren) = thisOctal%rhou(1:thisOctal%maxChildren) * factor
    endif

!
! momentum (v)
!
    if ( all(thisOctal%rhov(1:thisOctal%maxChildren) < 0.d0) .or. &
         all(thisOctal%rhov(1:thisOctal%maxChildren) > 0.d0) ) then
       oldMom = parent%rhov(parentSubcell) * cellVolume(parent,parentSubcell) 
       newMom = SUM(thisOctal%rhov(1:thisOctal%maxChildren)*cellVolume(thisOctal,1))
       factor = oldMom / newMom
       thisOctal%rhov(1:thisOctal%maxChildren) = thisOctal%rhov(1:thisOctal%maxChildren) * factor
    endif

!
! momentum (w)
!
    if ( all(thisOctal%rhow(1:thisOctal%maxChildren) < 0.d0) .or. &
         all(thisOctal%rhow(1:thisOctal%maxChildren) > 0.d0) ) then
       oldMom = parent%rhow(parentSubcell) * cellVolume(parent,parentSubcell) 
       newMom = SUM(thisOctal%rhow(1:thisOctal%maxChildren)*cellVolume(thisOctal,1))
       factor = oldMom / newMom
       thisOctal%rhow(1:thisOctal%maxChildren) = thisOctal%rhow(1:thisOctal%maxChildren) * factor
    endif


    if (PRESENT(constantGravity)) then
       if (constantGravity) then
          do i = 1, thisOctal%maxChildren
             thisOctal%phi_i(i) =  getPhiValue(thisOctal, i, grid%geometry)
          enddo
       endif
    endif
    
    grid%nOctals = grid%nOctals + 1

    ! check for a new maximum depth 
    IF (parent%child(newChildIndex)%nDepth > grid%maxDepth) THEN
       grid%maxDepth = parent%child(newChildIndex)%nDepth
       CALL setSmallestSubcell(grid)
    END IF

 end if

666 continue
  end subroutine addNewChildWithInterp
#endif



  subroutine interpTrilinear(grid, thisOctal, subcell, rho, rhoe, rhou, rhov, rhow, energy, phi, pressure)
    use inputs_mod, only : maxDepthAMR!, minDepthAmr
    use mpi
    type(OCTAL), pointer :: thisOctal
    type(GRIDTYPE) :: grid
    integer :: iCorner, iDir, nCorner, nDir
    integer :: nd, j, subcell
    type(VECTOR) :: dir(8), corner(8), position, rVec, centre
    real(double) :: rho, rhoe, rhou, rhov, rhow, e, phi, pressure,r, energy
    real(double) :: u,v,w
    real(double) :: rhoCorner(8)
    real(double) :: rhoeCorner(8)
    real(double) :: rhouCorner(8)
    real(double) :: rhovCorner(8)
    real(double) :: rhowCorner(8)
    real(double) :: eCorner(8)
    real(double) :: phiCorner(8)
    real(double) :: pressureCorner(8)
    real(double) :: weight, totalWeight
    real(double) :: xh, yh, zh, smallDist
    logical :: debug

    debug=.false.

    smallDist = 0.01d0*grid%octreeRoot%subcellSize / &
         2.0_oc**REAL(maxDepthAmr,kind=oct)

    centre = subcellCentre(thisOctal%parent, thisOctal%parentSubcell)

    
    nDir = 8
    r = 0.1d0*smallDist
    dir(1) = VECTOR(-r, -r, -r)
    dir(2) = VECTOR(+r, -r, -r)
    dir(3) = VECTOR(-r, -r, +r)
    dir(4) = VECTOR(+r, -r, +r)
    dir(5) = VECTOR(-r, +r, -r)
    dir(6) = VECTOR(+r, +r, -r)
    dir(7) = VECTOR(-r, +r, +r)
    dir(8) = VECTOR(+r, +r, +r)
       
    nCorner = 8
    r = thisOctal%subcellSize
    corner(1) = centre + VECTOR(-r, -r, -r)
    corner(2) = centre + VECTOR(+r, -r, -r)
    corner(3) = centre + VECTOR(-r, -r, +r)
    corner(4) = centre + VECTOR(+r, -r, +r)
    corner(5) = centre + VECTOR(-r, +r, -r)
    corner(6) = centre + VECTOR(+r, +r, -r)
    corner(7) = centre + VECTOR(-r, +r, +r)
    corner(8) = centre + VECTOR(+r, +r, +r)


    rhoCorner = 0.d0
    rhoeCorner = 0.d0
    eCorner = 0.d0
    rhouCorner = 0.d0
    rhovCorner = 0.d0
    rhowCorner = 0.d0
    phicorner  = 0.d0
    pressureCorner = 0.d0

    do iCorner = 1, nCorner

       totalWeight = 0.d0
       j = 0

       do iDir = 1, nDir

          position = corner(iCorner) + dir(iDir)

          if (inOctal(grid%octreeRoot, position)) then
             call getHydroValues(grid, position, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, xh, yh, zh, pressure, .true.)
             weight = 1.d0/modulus(position-VECTOR(xh,yh,zh))
             totalWeight = totalWeight + weight
             rhoCorner(iCorner) = rhoCorner(iCorner) + weight * rho
             rhoeCorner(iCorner) = rhoeCorner(iCorner) + weight * rhoe
             rhouCorner(iCorner) = rhouCorner(iCorner) + weight * rhou
             rhovCorner(iCorner) = rhovCorner(iCorner) + weight * rhov
             rhowCorner(iCorner) = rhowCorner(iCorner) + weight * rhow
             eCorner(iCorner) = eCorner(iCorner) + weight * energy
             phiCorner(iCorner) = phiCorner(iCorner) + weight * phi
             pressureCorner(iCorner) = pressureCorner(iCorner) + weight * pressure
          endif

       enddo
       rhoCorner(iCorner) = rhoCorner(iCorner) / totalWeight
       rhoeCorner(iCorner) = rhoeCorner(iCorner) / totalWeight
       rhouCorner(iCorner) = rhouCorner(iCorner) / totalWeight
       rhovCorner(iCorner) = rhovCorner(iCorner) / totalWeight
       rhowCorner(iCorner) = rhowCorner(iCorner) / totalWeight
       eCorner(iCorner) = eCorner(iCorner) / totalWeight
       phiCorner(iCorner) = phiCorner(iCorner) / totalWeight
       pressureCorner(iCorner) = pressureCorner(iCorner) / totalWeight
    enddo
    
    rVec = subcellCentre(thisOctal, subcell)
    u = (rVec%x - corner(1)%x)/(corner(2)%x - corner(1)%x)
    v = (rVec%y - corner(4)%y)/(corner(5)%y - corner(4)%y)
    w = (rVec%z - corner(6)%z)/(corner(7)%z - corner(6)%z)
 
    rho = rhoCorner(1) * (1.d0 - u)*(1.d0 - v)*(1.d0 - w) + &
          rhoCorner(2) * (       u)*(1.d0 - v)*(1.d0 - w) + &
          rhoCorner(3) * (1.d0 - u)*(1.d0 - v)*(       w) + &
          rhoCorner(4) * (       u)*(1.d0 - v)*(       w) + &
          rhoCorner(5) * (1.d0 - u)*(       v)*(1.d0 - w) + &
          rhoCorner(6) * (       u)*(       v)*(1.d0 - w) + &
          rhoCorner(7) * (1.d0 - u)*(       v)*(       w) + &
          rhoCorner(8) * (       u)*(       v)*(       w) 

    rhoe = rhoeCorner(1) * (1.d0 - u)*(1.d0 - v)*(1.d0 - w) + &
           rhoeCorner(2) * (       u)*(1.d0 - v)*(1.d0 - w) + &
           rhoeCorner(3) * (1.d0 - u)*(1.d0 - v)*(       w) + &
           rhoeCorner(4) * (       u)*(1.d0 - v)*(       w) + &
           rhoeCorner(5) * (1.d0 - u)*(       v)*(1.d0 - w) + &
           rhoeCorner(6) * (       u)*(       v)*(1.d0 - w) + &
           rhoeCorner(7) * (1.d0 - u)*(       v)*(       w) + &
           rhoeCorner(8) * (       u)*(       v)*(       w) 

    rhou = rhouCorner(1) * (1.d0 - u)*(1.d0 - v)*(1.d0 - w) + &
           rhouCorner(2) * (       u)*(1.d0 - v)*(1.d0 - w) + &
           rhouCorner(3) * (1.d0 - u)*(1.d0 - v)*(       w) + &
           rhouCorner(4) * (       u)*(1.d0 - v)*(       w) + &
           rhouCorner(5) * (1.d0 - u)*(       v)*(1.d0 - w) + &
           rhouCorner(6) * (       u)*(       v)*(1.d0 - w) + &
           rhouCorner(7) * (1.d0 - u)*(       v)*(       w) + &
           rhouCorner(8) * (       u)*(       v)*(       w) 

    rhov = rhovCorner(1) * (1.d0 - u)*(1.d0 - v)*(1.d0 - w) + &
           rhovCorner(2) * (       u)*(1.d0 - v)*(1.d0 - w) + &
           rhovCorner(3) * (1.d0 - u)*(1.d0 - v)*(       w) + &
           rhovCorner(4) * (       u)*(1.d0 - v)*(       w) + &
           rhovCorner(5) * (1.d0 - u)*(       v)*(1.d0 - w) + &
           rhovCorner(6) * (       u)*(       v)*(1.d0 - w) + &
           rhovCorner(7) * (1.d0 - u)*(       v)*(       w) + &
           rhovCorner(8) * (       u)*(       v)*(       w) 

    rhow = rhowCorner(1) * (1.d0 - u)*(1.d0 - v)*(1.d0 - w) + &
           rhowCorner(2) * (       u)*(1.d0 - v)*(1.d0 - w) + &
           rhowCorner(3) * (1.d0 - u)*(1.d0 - v)*(       w) + &
           rhowCorner(4) * (       u)*(1.d0 - v)*(       w) + &
           rhowCorner(5) * (1.d0 - u)*(       v)*(1.d0 - w) + &
           rhowCorner(6) * (       u)*(       v)*(1.d0 - w) + &
           rhowCorner(7) * (1.d0 - u)*(       v)*(       w) + &
           rhowCorner(8) * (       u)*(       v)*(       w) 

    e    = eCorner(1) * (1.d0 - u)*(1.d0 - v)*(1.d0 - w) + &
           eCorner(2) * (       u)*(1.d0 - v)*(1.d0 - w) + &
           eCorner(3) * (1.d0 - u)*(1.d0 - v)*(       w) + &
           eCorner(4) * (       u)*(1.d0 - v)*(       w) + &
           eCorner(5) * (1.d0 - u)*(       v)*(1.d0 - w) + &
           eCorner(6) * (       u)*(       v)*(1.d0 - w) + &
           eCorner(7) * (1.d0 - u)*(       v)*(       w) + &
           eCorner(8) * (       u)*(       v)*(       w) 

    phi  = phiCorner(1) * (1.d0 - u)*(1.d0 - v)*(1.d0 - w) + &
           phiCorner(2) * (       u)*(1.d0 - v)*(1.d0 - w) + &
           phiCorner(3) * (1.d0 - u)*(1.d0 - v)*(       w) + &
           phiCorner(4) * (       u)*(1.d0 - v)*(       w) + &
           phiCorner(5) * (1.d0 - u)*(       v)*(1.d0 - w) + &
           phiCorner(6) * (       u)*(       v)*(1.d0 - w) + &
           phiCorner(7) * (1.d0 - u)*(       v)*(       w) + &
           phiCorner(8) * (       u)*(       v)*(       w) 

    pressure  = pressureCorner(1) * (1.d0 - u)*(1.d0 - v)*(1.d0 - w) + &
           pressureCorner(2) * (       u)*(1.d0 - v)*(1.d0 - w) + &
           pressureCorner(3) * (1.d0 - u)*(1.d0 - v)*(       w) + &
           pressureCorner(4) * (       u)*(1.d0 - v)*(       w) + &
           pressureCorner(5) * (1.d0 - u)*(       v)*(1.d0 - w) + &
           pressureCorner(6) * (       u)*(       v)*(1.d0 - w) + &
           pressureCorner(7) * (1.d0 - u)*(       v)*(       w) + &
           pressureCorner(8) * (       u)*(       v)*(       w) 



  end subroutine interpTrilinear

  subroutine interpBilinear(grid, thisOctal, subcell, rho, rhoe, rhou, rhov, rhow, energy, phi, pressure)
    use inputs_mod, only : maxDepthAMR!, minDepthAmr
    use mpi
    type(OCTAL), pointer :: thisOctal
    type(GRIDTYPE) :: grid
    integer :: iCorner, iDir, nCorner, nDir
    integer :: nd, j, subcell
    type(VECTOR) :: dir(8), corner(8), position, rVec, centre
    real(double) :: rho, rhoe, rhou, rhov, rhow, e, phi, pressure,r, energy
    real(double) :: u,w
    real(double) :: rhoCorner(8)
    real(double) :: rhoeCorner(8)
    real(double) :: rhouCorner(8)
    real(double) :: rhovCorner(8)
    real(double) :: rhowCorner(8)
    real(double) :: eCorner(8)
    real(double) :: phiCorner(8)
    real(double) :: pressureCorner(8)
    real(double) :: weight, totalWeight
    real(double) :: xh, yh, zh, smallDist
    logical :: debug

    debug=.false.

    smallDist = 0.01d0*grid%octreeRoot%subcellSize / &
         2.0_oc**REAL(maxDepthAmr,kind=oct)

    centre = subcellCentre(thisOctal%parent, thisOctal%parentSubcell)

    
    nDir = 4
    r = 0.1d0*smallDist
    dir(1) = VECTOR(-r, 0.d0,  -r)
    dir(2) = VECTOR(+r, 0.d0, -r)
    dir(3) = VECTOR(-r, 0.d0, +r)
    dir(4) = VECTOR(+r, 0.d0, +r)
       
    nCorner = 4
    r = thisOctal%subcellSize
    corner(1) = centre + VECTOR(-r, 0.d0, -r)
    corner(2) = centre + VECTOR(+r, 0.d0, -r)
    corner(3) = centre + VECTOR(-r, 0.d0, +r)
    corner(4) = centre + VECTOR(+r, 0.d0, +r)


    rhoCorner = 0.d0
    rhoeCorner = 0.d0
    eCorner = 0.d0
    rhouCorner = 0.d0
    rhovCorner = 0.d0
    rhowCorner = 0.d0
    phicorner  = 0.d0
    pressureCorner = 0.d0

    do iCorner = 1, nCorner

       totalWeight = 0.d0
       j = 0

       do iDir = 1, nDir

          position = corner(iCorner) + dir(iDir)

          if (inOctal(grid%octreeRoot, position)) then
             call getHydroValues(grid, position, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, xh, yh, zh, pressure, .true.)
             weight = 1.d0/modulus(position-VECTOR(xh,yh,zh))
             totalWeight = totalWeight + weight
             rhoCorner(iCorner) = rhoCorner(iCorner) + weight * rho
             rhoeCorner(iCorner) = rhoeCorner(iCorner) + weight * rhoe
             rhouCorner(iCorner) = rhouCorner(iCorner) + weight * rhou
             rhovCorner(iCorner) = rhovCorner(iCorner) + weight * rhov
             rhowCorner(iCorner) = rhowCorner(iCorner) + weight * rhow
             eCorner(iCorner) = eCorner(iCorner) + weight * energy
             phiCorner(iCorner) = phiCorner(iCorner) + weight * phi
             pressureCorner(iCorner) = pressureCorner(iCorner) + weight * pressure
          endif

       enddo
       rhoCorner(iCorner) = rhoCorner(iCorner) / totalWeight
       rhoeCorner(iCorner) = rhoeCorner(iCorner) / totalWeight
       rhouCorner(iCorner) = rhouCorner(iCorner) / totalWeight
       rhovCorner(iCorner) = rhovCorner(iCorner) / totalWeight
       rhowCorner(iCorner) = rhowCorner(iCorner) / totalWeight
       eCorner(iCorner) = eCorner(iCorner) / totalWeight
       phiCorner(iCorner) = phiCorner(iCorner) / totalWeight
       pressureCorner(iCorner) = pressureCorner(iCorner) / totalWeight
    enddo
    
    rVec = subcellCentre(thisOctal, subcell)
    u = (rVec%x - corner(1)%x)/(corner(2)%x - corner(1)%x)
    w = (rVec%z - corner(4)%z)/(corner(5)%z - corner(4)%z)

    rho = rhoCorner(1) * (1.d0 - u)*(1.d0 - w) + &
          rhoCorner(2) * (       u)*(1.d0 - w) + &
          rhoCorner(3) * (1.d0 - u)*(       w) + &
          rhoCorner(4) * (       u)*(       w) 

    rhoe = rhoeCorner(1) * (1.d0 - u)*(1.d0 - w) + &
           rhoeCorner(2) * (       u)*(1.d0 - w) + &
           rhoeCorner(3) * (1.d0 - u)*(       w) + &
           rhoeCorner(4) * (       u)*(       w)

    rhou = rhouCorner(1) * (1.d0 - u)*(1.d0 - w) + &
           rhouCorner(2) * (       u)*(1.d0 - w) + &
           rhouCorner(3) * (1.d0 - u)*(       w) + &
           rhouCorner(4) * (       u)*(       w) 

    rhov = rhovCorner(1) * (1.d0 - u)*(1.d0 - w) + &
           rhovCorner(2) * (       u)*(1.d0 - w) + &
           rhovCorner(3) * (1.d0 - u)*(       w) + &
           rhovCorner(4) * (       u)*(       w) 

    rhow = rhowCorner(1) * (1.d0 - u)*(1.d0 - w) + &
           rhowCorner(2) * (       u)*(1.d0 - w) + &
           rhowCorner(3) * (1.d0 - u)*(       w) + &
           rhowCorner(4) * (       u)*(       w) 

    e    = eCorner(1) * (1.d0 - u)*(1.d0 - w) + &
           eCorner(2) * (       u)*(1.d0 - w) + &
           eCorner(3) * (1.d0 - u)*(       w) + &
           eCorner(4) * (       u)*(       w) 

    phi  = phiCorner(1) * (1.d0 - u)*(1.d0 - w) + &
           phiCorner(2) * (       u)*(1.d0 - w) + &
           phiCorner(3) * (1.d0 - u)*(       w) + &
           phiCorner(4) * (       u)*(       w) 

    pressure  = pressureCorner(1) * (1.d0 - u)*(1.d0 - w) + &
           pressureCorner(2) * (       u)*(1.d0 - w) + &
           pressureCorner(3) * (1.d0 - u)*(       w) + &
           pressureCorner(4) * (       u)*(       w) 

  end subroutine interpBilinear

  subroutine interplinear(grid, thisOctal, subcell, rho, rhoe, rhou, rhov, rhow, energy, phi, pressure)
    use inputs_mod, only : maxDepthAMR!, minDepthAmr
    use mpi
    type(OCTAL), pointer :: thisOctal
    type(GRIDTYPE) :: grid
    integer :: iCorner, iDir, nCorner, nDir
    integer :: nd, j, subcell
    type(VECTOR) :: dir(8), corner(8), position, rVec, centre
    real(double) :: rho, rhoe, rhou, rhov, rhow, e, phi, pressure,r, energy
    real(double) :: u
    real(double) :: rhoCorner(8)
    real(double) :: rhoeCorner(8)
    real(double) :: rhouCorner(8)
    real(double) :: rhovCorner(8)
    real(double) :: rhowCorner(8)
    real(double) :: eCorner(8)
    real(double) :: phiCorner(8)
    real(double) :: pressureCorner(8)
    real(double) :: weight, totalWeight
    real(double) :: xh, yh, zh, smallDist
    logical :: debug

    debug=.false.

    smallDist = 0.01d0*grid%octreeRoot%subcellSize / &
         2.0_oc**REAL(maxDepthAmr,kind=oct)

    centre = subcellCentre(thisOctal%parent, thisOctal%parentSubcell)
    
    nDir = 2
    r = 0.1d0*smallDist
    dir(1) = VECTOR(-r, 0.d0, 0.d0)
    dir(2) = VECTOR(+r, 0.d0, 0.d0)
       
    nCorner = 2
    r = thisOctal%subcellSize
    corner(1) = centre + VECTOR(-r, 0.d0, 0.d0)
    corner(2) = centre + VECTOR(+r, 0.d0, 0.d0)


    rhoCorner = 0.d0
    rhoeCorner = 0.d0
    eCorner = 0.d0
    rhouCorner = 0.d0
    rhovCorner = 0.d0
    rhowCorner = 0.d0
    phicorner  = 0.d0
    pressureCorner = 0.d0

    do iCorner = 1, nCorner

       totalWeight = 0.d0
       j = 0

       do iDir = 1, nDir

          position = corner(iCorner) + dir(iDir)

          if (inOctal(grid%octreeRoot, position)) then
             call getHydroValues(grid, position, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, xh, yh, zh, pressure, .true.)
             weight = 1.d0/abs(position%x - xh)
             totalWeight = totalWeight + weight
             rhoCorner(iCorner) = rhoCorner(iCorner) + weight * rho
             rhoeCorner(iCorner) = rhoeCorner(iCorner) + weight * rhoe
             rhouCorner(iCorner) = rhouCorner(iCorner) + weight * rhou
             rhovCorner(iCorner) = rhovCorner(iCorner) + weight * rhov
             rhowCorner(iCorner) = rhowCorner(iCorner) + weight * rhow
             eCorner(iCorner) = eCorner(iCorner) + weight * energy
             phiCorner(iCorner) = phiCorner(iCorner) + weight * phi
             pressureCorner(iCorner) = pressureCorner(iCorner) + weight * pressure
          endif

       enddo
       rhoCorner(iCorner) = rhoCorner(iCorner) / totalWeight
       rhoeCorner(iCorner) = rhoeCorner(iCorner) / totalWeight
       rhouCorner(iCorner) = rhouCorner(iCorner) / totalWeight
       rhovCorner(iCorner) = rhovCorner(iCorner) / totalWeight
       rhowCorner(iCorner) = rhowCorner(iCorner) / totalWeight
       eCorner(iCorner) = eCorner(iCorner) / totalWeight
       phiCorner(iCorner) = phiCorner(iCorner) / totalWeight
       pressureCorner(iCorner) = pressureCorner(iCorner) / totalWeight
    enddo
    
    rVec = subcellCentre(thisOctal, subcell)
    u = (rVec%x - corner(1)%x)/(corner(2)%x - corner(1)%x)
 
    rho = rhoCorner(1) * (1.d0 - u) + &
          rhoCorner(2) * (       u)

    rhoe = rhoeCorner(1) * (1.d0 - u) + &
           rhoeCorner(2) * (       u)

    rhou = rhouCorner(1) * (1.d0 - u) + &
           rhouCorner(2) * (       u)


    e    = eCorner(1) * (1.d0 - u) + &
           eCorner(2) * (       u)

    phi  = phiCorner(1) * (1.d0 - u) + &
           phiCorner(2) * (       u)

    pressure  = pressureCorner(1) * (1.d0 - u) + &
           pressureCorner(2)

  end subroutine interplinear

  subroutine returnAddNewChildPointArrays(thisOctal, topOctal, topOctalSubcell,  grid, rhoPoint, rhoePoint, &
       rhouPoint, rhovPoint, rhowPoint, phiPoint, energyPoint, pressurePoint, npoints)
    use inputs_mod, only : maxDepthAMR!, minDepthAmr
    use mpi
    type(OCTAL), pointer :: thisOctal, topOctal
    type(GRIDTYPE) :: grid
    integer :: iCorner, iDir, nCorner, nDir
    integer :: nd, topOctalSubcell, j
    type(VECTOR) :: dir(8), corner(8), position, rVec, centre
    real(double) :: weight, totalWeight
    real(double) :: rho, rhoe, rhou, rhov, rhow, r, energy, phi, pressure
    real(double) :: xh, yh, zh
    real(double) :: smallDist
    integer :: npoints
    integer, parameter :: maxpts = 10000
    real(double) :: xPoint(maxpts)
    real(double) :: yPoint(maxpts)
    real(double) :: zPoint(maxpts)
    real(double) :: rhoPoint(maxpts)
    real(double) :: rhoePoint(maxpts)
    real(double) :: rhouPoint(maxpts)
    real(double) :: rhovPoint(maxpts)
    real(double) :: rhowPoint(maxpts)
    real(double) :: phiPoint(maxpts)
    real(double) :: energyPoint(maxpts)
    real(double) :: pressurePoint(maxpts)
    logical :: addLocal


    if (thisOctal%threed) then

       smallDist = 0.01d0*grid%octreeRoot%subcellSize / &
                                2.0_oc**REAL(maxDepthAmr,kind=oct)

       centre = subcellCentre(topOctal, topOctalSubcell)

       nDir = 8
       r = 0.1d0*smallDist
       dir(1) = VECTOR(-r, -r, -r)
       dir(2) = VECTOR(+r, -r, -r)
       dir(3) = VECTOR(-r, -r, +r)
       dir(4) = VECTOR(+r, -r, +r)
       dir(5) = VECTOR(-r, +r, -r)
       dir(6) = VECTOR(+r, +r, -r)
       dir(7) = VECTOR(-r, +r, +r)
       dir(8) = VECTOR(+r, +r, +r)
       
       r = topOctal%subcellSize/2.d0

       nCorner = 8
       corner(1) = centre + VECTOR(-r, -r, -r)
       corner(2) = centre + VECTOR(+r, -r, -r)
       corner(3) = centre + VECTOR(-r, -r, +r)
       corner(4) = centre + VECTOR(+r, -r, +r)
       corner(5) = centre + VECTOR(-r, +r, -r)
       corner(6) = centre + VECTOR(+r, +r, -r)
       corner(7) = centre + VECTOR(-r, +r, +r)
       corner(8) = centre + VECTOR(+r, +r, +r)
    endif

    if (thisOctal%twod) then
       nDir = 4
       r = 0.1d0*grid%halfSmallestSubcell
       dir(1) = VECTOR(-r, 0.d0, -r)
       dir(2) = VECTOR(+r, 0.d0, -r)
       dir(3) = VECTOR(+r, 0.d0, +r)
       dir(4) = VECTOR(-r, 0.d0, +r)
       centre = subcellCentre(topOctal, topOctalSubcell)


       nCorner = 4
       r = topOctal%subcellSize/2.d0
       corner(1) = thisOctal%centre + VECTOR(-r, 0.d0, -r)
       corner(2) = thisOctal%centre + VECTOR(+r, 0.d0, -r)
       corner(3) = thisOctal%centre + VECTOR(-r, 0.d0, +r)
       corner(4) = thisOctal%centre + VECTOR(+r, 0.d0, +r)
    endif

    if (thisOctal%oneD) then
       nDir = 2
       r = 0.1d0*grid%halfSmallestSubcell
       dir(1) = VECTOR(-r, 0.d0, 0.d0)
       dir(2) = VECTOR(+r, 0.d0, 0.d0)
       centre = subcellCentre(topOctal, topOctalSubcell)
       
       nCorner = 2
       r = topOctal%subcellSize/2.d0
       corner(1) = centre + VECTOR(-r, 0.d0, 0.d0)
       corner(2) = centre + VECTOR(+r, 0.d0, 0.d0)
    endif

    nPoints = 0
    addLocal = .true.
    do iCorner = 1, nCorner
       totalWeight = 0.d0
       j = 0
       do iDir = 1, nDir

          if (iDir == iCorner) cycle
          position = corner(iCorner) + dir(iDir)

          if (inOctal(grid%octreeRoot, position).and.(.not.inSubcell(topOctal, topOctalSubcell, position))) then
             call getHydroValues(grid, position, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, xh, yh, zh, pressure, .true.)

             weight = 1.d0
!            weight = abs(parent%ndepth - nd)+1.d0
             totalWeight = totalWeight + weight

             nPoints = nPoints + 1
             xPoint(nPoints) = xh
             yPoint(nPoints) = yh
             zPoint(nPoints) = zh
             rhoPoint(nPoints) = rho
             rhoePoint(nPoints) = rhoe
             rhouPoint(nPoints) = rhou
             rhovPoint(nPoints) = rhov
             rhowPoint(nPoints) = rhow
             energyPoint(nPoints) = energy
             phiPoint(nPoints) = phi
             pressurePoint(nPoints) = pressure
          else
             weight = 1.d0
             totalWeight = totalWeight + weight

             rVec = subcellCentre(topOctal, topOctalSubcell)

             if (addlocal) then
                nPoints = nPoints + 1
                xPoint(nPoints) = rVec%x
                yPoint(nPoints) = rVec%y
                zPoint(nPoints) = rVec%z

                rhoPoint(nPoints) = topOctal%rho(topOctalSubcell)
                rhoePoint(nPoints) = topOctal%rhoe(topOctalSubcell)
                rhouPoint(nPoints) = topOctal%rhou(topOctalSubcell)
                rhovPoint(nPoints) = topOctal%rhov(topOctalSubcell)
                rhowPoint(nPoints) = topOctal%rhow(topOctalSubcell)
                energyPoint(nPoints) = topOctal%energy(topOctalSubcell)
                phiPoint(nPoints) = topOctal%phi_i(topOctalSubcell)
                pressurePoint(nPoints) = topOctal%pressure_i(topOctalSubcell)
                addlocal = .false.
             endif
             j = j + 1
          endif
       enddo
    enddo
  end subroutine returnAddNewChildPointArrays


  subroutine getPointsInRadius(position, radius, grid, npoints, rho, rhoe, rhou, rhov, rhow, energy, pressure, phi, x, y, z)
    use mpi
    type(GRIDTYPE) :: grid
    real(double) :: rho(:), rhoe(:), rhou(:), rhov(:), rhow(:), energy(:), pressure(:), phi(:), x(:), y(:), z(:)
    type(VECTOR) :: position
    real(double) :: radius
    integer :: nPoints
    type(OCTAL), pointer :: thisOctal
    integer, parameter :: nStorage = 12
    real(double) :: storageArray(nStorage)
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr, ithread
    integer :: counter, nvals
    real(double) :: loc(7)
    integer :: evenUpArray(nHydroThreadsGlobal)
    integer :: i

    nPoints = 0
    thisOctal => grid%octreeRoot
    rho = 0.d0
    rhoe = 0.d0
    rhou = 0.d0
    rhov = 0.d0
    rhow = 0.d0
    energy = 0.d0
    pressure = 0.d0
    phi = 0.d0
    x = 0.d0
    y = 0.d0
    z = 0.d0
    call getPointsInRadiusLocal(position, radius, thisOctal, npoints, rho, rhoe, rhou, rhov, rhow, energy, pressure, phi, x, y, z)
!    if (abs(position%x-0.5d0)<0.05) then
!       write(*,*) myrankGlobal, " got ",nPoints, " locally"
!    endif


!    do i = 1, nPoints
!       do counter = 1, nPoints
!          if( i /= counter) then
!             if(x(i) == x(counter) .and. y(i) == y(counter) .and. z(i) == z(counter)) then
!                print *, "x", x(i), x(counter)
!                print *, "y", y(i), y(counter)
!                print *, "z", z(i), z(counter)
!                call torus_abort("DUPLICATE ENTRY A")
!             end if
!          end if
!       end do
!    end do

   
    call setupEvenUpArray(grid, evenUpArray)

    do i = 1, nHydroThreadsGlobal
       if(evenupArray(i) /= evenupArray(myRankGlobal)) then
          !Found a serving thread
          iThread = i
          exit
       end if
   end do

    loc(1) = position%x
    if(grid%octreeroot%twoD) then
       loc(2) = 0.d0
    else
       loc(2) = position%y
    end if
    loc(3) = position%z
    loc(4) = dble(myRankGlobal)
    loc(5) = 1.d0
    loc(6) = radius
    !7 is the identifier for whether or not to use topOctal. In this routine this will always be the case
    loc(7) = 1.d0
    call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
    call MPI_RECV(nvals, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, status, ierr) 
!    if (abs(position%x-0.5d0)<0.05) then
!       write(*,*) myrankGlobal, " got ",nVals, " from ",ithread
!    endif


    if (nVals > 0) then
       do counter = 1, nvals
          call MPI_RECV(storageArray, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
          nPoints = nPoints + 1
          rho(npoints) = storageArray(2)
          rhoe(npoints) = storageArray(3)
          rhou(npoints) = storageArray(4)
          rhov(npoints) = storageArray(5)
          rhow(npoints) = storageArray(6)
          energy(npoints) = storageArray(7)
          phi(npoints) = storageArray(8)
          pressure(npoints) = storageArray(12)
          x(npoints) = storageArray(9)
          y(npoints) = storageArray(10)
          z(npoints) = storageArray(11)
       end do
    endif


   
  end subroutine getPointsInRadius


  recursive subroutine getPointsInRadiusLocal(position, radius, thisOctal, npoints, rho, rhoe, rhou, rhov, rhow, &
       energy, pressure, phi, x, y, z)
    real(double) :: rho(:), rhoe(:), rhou(:), rhov(:), rhow(:), energy(:), pressure(:), phi(:), x(:), y(:), z(:)
    type(VECTOR) :: position
    real(double) :: radius, r
    integer :: nPoints
    type(OCTAL), pointer :: thisOctal, child, topOctal
    integer :: subcell, i, topOctalSubcell, k
    type(VECTOR) :: cen
    logical :: changed, check

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  getPointsInRadiusLocal(position, radius, child, npoints, rho, rhoe, rhou, rhov, rhow, energy, &
                     pressure, phi, x, y, z)
                exit
             end if
          end do
       else 
          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then


             changed = .false.
             topOctal => thisOctal
             topOctalSubcell = subcell
             do while(topOctal%changed(topOctalSubcell))
                topOctalSubcell = topOctal%parentSubcell
                topOctal => topOctal%parent
                changed = .true.
             enddo

             cen = subcellCentre(topOctal, topOctalsubcell)
             r = modulus(cen - position)
             if (r < radius) then
                check = .true.
                do k = 1, nPoints
                   if(cen%x == x(k) .and. cen%y == y(k) .and. cen%z == z(k)) then
                      check = .false.
                   end if
                end do
                if(check) then
                   nPoints = nPoints + 1
                   rho(nPoints) = topOctal%rho(topOctalsubcell)
                   rhoe(nPoints) = topOctal%rhoe(topOctalsubcell)
                   rhou(nPoints) = topOctal%rhou(topOctalsubcell)
                   rhov(nPoints) = topOctal%rhov(topOctalsubcell)
                   rhow(nPoints) = topOctal%rhow(topOctalsubcell)
                   energy(nPoints) = topOctal%energy(topOctalsubcell)
                   pressure(nPoints) = topOctal%pressure_i(topOctalsubcell)
                   phi(nPoints) = topOctal%phi_gas(topOctalsubcell)
                   x(nPoints) = cen%x
                   y(nPoints) = cen%y
                   z(nPoints) = cen%z
                end if
             endif
             if (changed) exit
          endif
       endif
    enddo
  end subroutine getPointsInRadiusLocal


  subroutine shutdownServers()
    use mpi
    integer :: iThread
    real(double) :: loc(3)
    integer, parameter :: tag = 50
    integer :: ierr

    do iThread = 1, nHydroThreadsGlobal
       if (iThread /= myrankGlobal) then
          loc = 1.d30
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine shutdownServers

  subroutine shutdownServersXRAY_PDR()
    use mpi
    integer :: iThread
    real(double) :: loc(7)
    integer, parameter :: tag = 50
    integer :: ierr

    do iThread = 1, nHydroThreadsGlobal
       if (iThread /= myrankGlobal) then
          loc = 1.d30
          loc(4) = 1.d0
          call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine shutdownServersXRAY_PDR

  subroutine shutdownServersXRAY()
    use mpi
    integer :: iThread
    real(double) :: loc(8)
    integer, parameter :: tag = 50
    integer :: ierr

    do iThread = 1, nHydroThreadsGlobal
       if (iThread /= myrankGlobal) then
          loc = 1.d30
          loc(4) = 1.d0
          call MPI_SEND(loc, 8, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine shutdownServersXRAY

  subroutine shutdownServers2(check, k, endloop)
    use mpi
    integer :: iThread
    real(double) :: loc(7)
    integer, parameter :: tag = 50
    integer :: ierr, endloop, k
    integer :: check(endloop, nHydroThreadsGlobal)
   
    do iThread = 1, nHydroThreadsGlobal
       if (iThread /= myrankGlobal .and. .not. ANY(iThread == check(k,1:nHydroThreadsGlobal))) then
          loc = 1.d30
          loc(4) = dble(myRankGlobal)
          loc(5) = 0.d0
          loc(6) = 0.d0
          loc(7) = 1.d0
          call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine shutdownServers2


!Set up an array of identifiers to speed up the even up grid process 
  subroutine setupEvenUpArray(grid, evenUpArray)
    type(gridtype) :: grid
    integer :: evenUpArray(nHydroThreadsGlobal)

!an evenuparray contains an element for every domain on the grid.
!the integer for each element specifies in which wave of the even
!up process that thread does its work.
!This allows threads that do not rely on each other to do their
!evening up at the same time
!e.g. (thread1, thread2, thread3) ---> (1, 2, 1)
!would mean that thread1 and thread3 do their evening up at the same time
!in wave 1 and thread 2 does its evening up separately in wave 2.

    if(grid%octreeRoot%twoD) then
       if((nHydrothreadsGlobal) == 4) then
          evenUpArray = (/1, 2, 3, 4/)
       else if((nHydrothreadsGlobal) == 16) then
!          evenUpArray = (/1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4/)
!new version so as to not screw up shephard's method
!          evenUpArray = (/1, 2, 3, 4, 5, 1, 6, 3, 7, 8, 1, 2, 9, 7, 5, 1/)
          evenUpArray = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)
!          evenUpArray = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)
       else if((nHydrothreadsGlobal) == 64) then
          !new version so as to not screw up shephard's method
          evenUpArray = (/1, 2, 3, 4, 5, 6,  7, 8, 9, 10, 11, 12, 13, 14, 15, 16, & !pane 1
               17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, & !pane 2           
               33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, & !pane 3           
               49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64/) !pane 4      
       end if
    else if(grid%octreeRoot%threeD) then
       if((nHydrothreadsGlobal) == 8) then
          evenUpArray = (/1, 2, 3, 4, 5, 6, 7, 8/)
       else if((nHydrothreadsGlobal) == 64) then
          !new version so as to not screw up shephard's method
          evenUpArray = (/1, 2, 3, 4, 10, 11, 12, 13, 5, 1, 6, 3, 14, 10, 15, 12, & !pane 1
               7, 8, 1, 2, 16, 17, 10, 11, 9, 7, 5, 1, 18, 16, 14, 10, & !pane 2           
               1, 2, 3, 4, 10, 11, 12, 13, 5, 1, 6, 3, 14, 10, 15, 12, & !pane 3           
               7, 8, 1, 2, 16, 17, 10, 11, 9, 7, 5, 1, 18, 16, 14, 10/) !pane 4      
               
       else if(nHydroTHreadsGlobal==512) then
          call autoSetupEvenupArray(grid, evenUpArray)

       else
          write(*,*) "nHydroThreadsGlobal ",nHydroThreadsGlobal
          call torus_abort("unknown no. of hydro threads in setupEvenUpArray")
       end if
    else
       evenUpArray = (/1, 2/)
    end if

  end subroutine setupEvenUpArray


  subroutine getHydroValues(grid, position, nd, rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z, pressure, useTop)
!       radially, searchRadius)
    use mpi
    type(GRIDTYPE) :: grid
    integer, intent(out) :: nd
    real(double), intent(out) :: rho, rhoe, rhou, rhov, rhow, energy, phi, x, y, z, pressure
    type(VECTOR) :: position, rVec
    real(double) :: loc(7)
    type(OCTAL), pointer :: thisOctal, topOctal
    integer :: iThread, topOctalSubcell
    integer, parameter :: nStorage = 12
    real(double) :: tempStorage(nStorage)
    integer :: subcell
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr
    logical :: useTop

    thisOctal => grid%octreeRoot
    call findSubcellLocal(position, thisOctal, subcell)
   
    if (octalOnThread(thisOctal, subcell, myrankGlobal)) then       

       if(useTop) then
          topOctal => thisOctal
          topOctalSubcell = subcell
          do while(topOctal%changed(topOctalSubcell))
             topOctalSubcell = topOctal%parentSubcell
             topOctal => topOctal%parent
          enddo
          rVec = subcellCentre(topOctal, topOctalsubcell)
          rho = topOctal%rho(topOctalsubcell)
          rhoe = topOctal%rhoe(topOctalsubcell)
          rhou = topOctal%rhou(topOctalsubcell)
          rhov = topOctal%rhov(topOctalsubcell)
          rhow = topOctal%rhow(topOctalsubcell)
          nd = thisOctal%nDepth ! this is needed for refinement
          energy = topOctal%energy(topOctalsubcell)
          phi = topOctal%phi_gas(topOctalsubcell)
          x = rVec%x
          y = rVec%y
          z = rVec%z
          pressure = topOctal%pressure_i(topOctalsubcell)
       else
          rVec = subcellCentre(thisOctal, subcell)
          rho = thisOctal%rho(subcell)
          rhoe = thisOctal%rhoe(subcell)
          rhou = thisOctal%rhou(subcell)
          rhov = thisOctal%rhov(subcell)
          rhow = thisOctal%rhow(subcell)
          nd = thisOctal%nDepth ! this is needed for refinement
          energy = thisOctal%energy(subcell)
          phi = thisOctal%phi_gas(subcell)
          x = rVec%x
          y = rVec%y
          z = rVec%z
          pressure = thisOctal%pressure_i(subcell)
       end if

    else

       iThread = thisOctal%mpiThread(subcell)
       loc(1) = position%x
       loc(2) = position%y
       loc(3) = position%z
       loc(4) = dble(myRankGlobal)
       loc(5) = 0.d0
       loc(6) = 0.d0
       if(useTop) then
          loc(7) = 1.d0
       else
          loc(7) = 0.d0
       end if
!       print *, "RANK ", myRankGlobal, "SENDING TO ", iThread, " usetop ",usetop,octalOnThread(thisOctal,subcell, myrankGlobal), thisOctal%nDepth
       call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
!       print *, "RANK ", myRankGlobal, "SENT"
       call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
  !     print *, "RANK ", myRankGlobal, "RECVING"
       nd = nint(tempStorage(1))
       rho = tempStorage(2)
       rhoe = tempStorage(3)
       rhou = tempStorage(4)
       rhov = tempStorage(5)
       rhow = tempStorage(6)
       energy = tempStorage(7)
       phi = tempStorage(8)
       x = tempstorage(9)
       y = tempstorage(10)
       z = tempstorage(11)
       pressure = tempstorage(12)
!       print *, "position ", position
!       print *, "rVec ", x, y, z
!       print *, "rho ", rho

    endif
  end subroutine getHydroValues

#ifdef PDR
  subroutine getRayTracingValuesPDR_TWO(grid, inputposition, direction, CII_POP, CI_POP, OI_POP, C12O_POP, tVal,&
       temperature)
!       HplusFrac)

    use mpi
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: CII_POP(5), CI_POP(5), OI_POP(5), C12O_POP(41), temperature
    type(VECTOR) :: position, rVec, direction
    type(VECTOR), intent(in) :: inputposition
    real(double) :: loc(7)
    type(OCTAL), pointer :: thisOctal
    integer :: iThread
    integer, parameter :: nStorage = 58
!    integer, parameter :: nStorage = 6
    real(double) :: tempStorage(nStorage), tval
    integer :: subcell
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr
!    logical :: useTop
    integer :: counter
    thisOctal => grid%octreeRoot
    position = inputposition

    call findSubcellLocal(position, thisOctal, subcell)
   
    if (octalOnThread(thisOctal, subcell, myrankGlobal)) then       
!58
       rVec = subcellCentre(thisOctal, subcell)
!       Hplusfrac = thisOctal%ionfrac(subcell, 2)
       temperature = thisOctal%temperature(subcell)
       CII_POP = thisOctal%CII_POP(subcell, :)!5
       CI_POP = thisOctal%CI_POP(subcell, :)!5
       OI_POP = thisOctal%OI_POP(subcell, :)!5
       C12O_POP = thisOctal%C12O_POP(subcell, :)!41
       call distanceToCellBoundary(grid, position, direction, tVal, thisOctal)
!       print *, "CII_POP ONTHREAD ", CII_POP

    else

       iThread = thisOctal%mpiThread(subcell)
!       print *, "myrankglobal ", myrankglobal, "to ", ithread
       loc(1) = position%x
       loc(2) = position%y
       loc(3) = position%z
       loc(4) = dble(myRankGlobal)
!       loc(4) = myRankGlobal
       loc(5) = direction%x
       loc(6) = direction%y
       loc(7) = direction%z

       call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
 !      print *, "sent "
       call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, &
            localWorldCommunicator, status, ierr)
  !     print *, "recvd"
       tval = tempstorage(1)
!       Hplusfrac = tempstorage(2)
       temperature = tempstorage(2)
       do counter = 1, 5
          CII_POP(counter) = tempstorage(counter + 2)
          CI_POP(counter) = tempstorage(counter + 7)
          OI_POP(counter) = tempstorage(counter + 12)
       enddo
       do counter = 1, 41
          C12O_POP(counter) = tempstorage(counter + 17)
       enddo


!       print *, "CII_POP NOTONTHREAD ", CII_POP
    endif
  end subroutine getRayTracingValuesPDR_TWO


  subroutine rayTracingServerPDR_TWO(grid)
    use mpi
    type(GRIDTYPE) :: grid
    logical :: stillServing
    real(double) :: loc(7)
    type(VECTOR) :: position, rVec, direction
    type(OCTAL), pointer :: thisOctal
!    type(OCTAL), pointer :: topOctal
    integer :: subcell!, nworking!, topOctalSubcell
    integer :: iThread!, servingArray!, workingTHreads(nworking)
    integer, parameter :: nStorage = 58
!    integer, parameter :: nStorage = 6
    real(double) :: tempStorage(nStorage), tval!, tmpthread
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr
    integer :: counter

    stillServing = .true.
!    servingArray = 0
!    workingTHreads = 0
    thisOctal => grid%octreeroot
    do while (stillServing)

!       do iThread = 1, nThreadsGlobal-1
!       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
       
       call MPI_RECV(loc, 7, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, localWorldCommunicator, status, ierr)
!       print *, "got data"
       position%x = loc(1)
       position%y = loc(2)
       position%z = loc(3)
       iThread = int(loc(4))

       direction%x = loc(5)
       direction%y = loc(6)
       direction%z = loc(7)
!       ithread = int(tmpThread)
       if (position%x > 1.d29) then          
          !          do j=1, nworking
          !             if(workingThreads(j) == 0) then
          !                workingThreads(j) = 1
          !                exit
          !             end if
          !         end do
          !         if(SUM(workingTHreads) == nworking) then
          !            workingThreads = 0
          stillServing= .false.
!       end if
       else

          call findSubcellLocal(position, thisOctal, subcell)
!          topOctal => thisOctal
!          topOctalSubcell = subcell
!          Do while(topOctal%changed(topOctalSubcell))
!             topOctalSubcell = topOctal%parentSubcell
!             topOctal => topOctal%parent
!          enddo
          call distanceToCellBoundary(grid, position, direction, tVal, thisOctal)
          rVec = subcellCentre(thisOctal, subcell)
          tempstorage(1) = tval
!          tempstorage(2) = thisOctal%ionfrac(subcell, 2)
          tempstorage(2) = thisOctal%temperature(subcell)
          do counter = 1, 5
             tempstorage(counter + 2) = thisOctal%CII_POP(subcell, counter)
             tempstorage(counter + 7) = thisOctal%CI_POP(subcell, counter)
             tempstorage(counter + 12) = thisOctal%oI_POP(subcell, counter)
          enddo
          do counter = 1, 41
             tempstorage(counter+17) = thisOctal%C12O_POP(subcell,counter)
          enddo      


          call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, &
               localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine rayTracingServerPDR_TWO



  subroutine getRayTracingValuesPDR(grid, inputposition, direction, rho, uvx, uvy, uvz, &
       uvx2, uvy2, uvz2, temperature, tval, abundanceArray)
!       radially, searchRadius)
    use mpi
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: rho, uvx, uvy, uvz, temperature, uvx2, uvy2, uvz2
    type(VECTOR) :: position, rVec, direction
    type(VECTOR), intent(in) :: inputposition
    real(double) :: loc(7)
    type(OCTAL), pointer :: thisOctal
    integer :: iThread
    integer, parameter :: nStorage = 42
!    integer, parameter :: nStorage = 6
    real(double) :: tempStorage(nStorage), tval, abundanceArray(1:33)
    integer :: subcell
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr
!    logical :: useTop
    logical :: notOnThread

    notonthread = .false.

    position = inputPosition

    if(inOctal(grid%octreeRoot, position)) then
       thisOctal => grid%octreeRoot
!       print *, "finding cell at ", position
       call findSubcellLocal(position, thisOctal, subcell)
    else
       print *, "OUTSIDE GRID ERROR"
       stop
       notOnThread = .true.
    endif

    if(notOnThread) then
       iThread = thisOctal%mpiThread(subcell)
       !       print *, "myrankglobal ", myrankglobal, "to ", ithread
       loc(1) = position%x
       loc(2) = position%y
       loc(3) = position%z
       loc(4) = dble(myRankGlobal)
       !       loc(4) = myRankGlobal
       loc(5) = direction%x
       loc(6) = direction%y
       loc(7) = direction%z
       
       call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       !      print *, "sent "
       call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, &
            localWorldCommunicator, status, ierr)
       !     print *, "recvd"
       rho = tempStorage(1)
       uvx = tempstorage(2)
       uvy = tempstorage(3)
       uvz = tempstorage(4)
       tval = tempstorage(5)
!       Hplusfrac = tempstorage(6)
       temperature = tempstorage(6)
       abundancearray(1:33) = tempstorage(7:39)
       uvx2 = tempstorage(40)
       uvy2 = tempstorage(41)
       uvz2 = tempstorage(42)


       if(tval == 0.d0) then
          print *, "tval zero A"
          stop
       endif

    elseif (octalOnThread(thisOctal, subcell, myrankGlobal)) then       

       rVec = subcellCentre(thisOctal, subcell)
       rho = thisOctal%rho(subcell)
!       uvx = thisOctal%uvvector(subcell)%x
!       uvy = thisOctal%uvvector(subcell)%y
!       uvz = thisOctal%uvvector(subcell)%z
       uvx = thisOctal%uvvectorPlus(subcell)%x
       uvy = thisOctal%uvvectorPlus(subcell)%y
       uvz = thisOctal%uvvectorPlus(subcell)%z
       uvx2 = thisOctal%uvvectorMinus(subcell)%x
       uvy2 = thisOctal%uvvectorMinus(subcell)%y
       uvz2 = thisOctal%uvvectorMinus(subcell)%z
       temperature = thisOctal%temperature(subcell)
       abundanceArray(:) = thisOctal%abundance(subcell, :)
       call distanceToCellBoundary(grid, position, direction, tVal, thisOctal)

       if(tval == 0.d0) then
          print *, "tval zero A"
          stop
       endif

    else
!       call torus_abort("error in pdr server")
       iThread = thisOctal%mpiThread(subcell)
!       print *, "myrankglobal ", myrankglobal, "to ", ithread
       loc(1) = position%x
       loc(2) = position%y
       loc(3) = position%z
       loc(4) = dble(myRankGlobal)
!       loc(4) = myRankGlobal
       loc(5) = direction%x
       loc(6) = direction%y
       loc(7) = direction%z

       call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
 !      print *, "sent "
       call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, &
            localWorldCommunicator, status, ierr)
  !     print *, "recvd"
       rho = tempStorage(1)
       uvx = tempstorage(2)
       uvy = tempstorage(3)
       uvz = tempstorage(4)
       tval = tempstorage(5)
       temperature = tempstorage(6)
       abundancearray(1:33) = tempstorage(7:39)
       uvx2 = tempstorage(40)
       uvy2 = tempstorage(41)
       uvz2 = tempstorage(42)

       if(tval == 0.d0) then
          print *, "tval zero A"
          stop
       endif

    endif
  end subroutine getRayTracingValuesPDR


  subroutine rayTracingServerPDR(grid)
    use mpi
    type(GRIDTYPE) :: grid
    logical :: stillServing
    real(double) :: loc(7)
    type(VECTOR) :: position, rVec, direction
    type(OCTAL), pointer :: thisOctal
!    type(OCTAL), pointer :: topOctal
    integer :: subcell!, nworking!, topOctalSubcell
    integer :: iThread!, servingArray!, workingTHreads(nworking)
    integer, parameter :: nStorage = 42
!    integer, parameter :: nStorage = 6
    real(double) :: tempStorage(nStorage), tval!, tmpthread
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr

    stillServing = .true.
!    servingArray = 0
!    workingTHreads = 0
    thisOctal => grid%octreeroot
    do while (stillServing)

!       do iThread = 1, nThreadsGlobal-1
!       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
       
       call MPI_RECV(loc, 7, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, localWorldCommunicator, status, ierr)
!       print *, "got data"
       position%x = loc(1)
       position%y = loc(2)
       position%z = loc(3)
       iThread = int(loc(4))

!       tmpThread = loc(4)
!!a temporary measure becau!se debug mode doesnt like int(dble) above...
!       if(tmpthread == 1.d!0) then
!          iThread = 1
!       elseif(tmpthread == 2.d0) then
!          iThread = 2
!       elseif(tmpthread == 3.d0) then
!          iThread = 3
!       elseif(tmpthread == 4.d0) then
!          iThread = 4
!       elseif(tmpthread == 5.d0) then
!          iThread = 5
!       elseif(tmpthread == 6.d0) then
!          iThread = 6
!       elseif(tmpthread == 7.d0) then
!          iThread = 7
!       elseif(tmpthread == 8.d0) then
!          iThread = 8
!       elseif(tmpthread == 9.d0) then
!          iThread = 9
!       else
!          print *, "tmpThread", tmpThread
!          call torus_abort("error in pdr server in debug mode")
!       endif

!       tmpThread = loc(4)
       direction%x = loc(5)
       direction%y = loc(6)
       direction%z = loc(7)
!       ithread = int(tmpThread)
       if (position%x > 1.d29) then          
          !          do j=1, nworking
          !             if(workingThreads(j) == 0) then
          !                workingThreads(j) = 1
          !                exit
          !             end if
          !         end do
          !         if(SUM(workingTHreads) == nworking) then
          !            workingThreads = 0
          stillServing= .false.
!       end if
       else

          call findSubcellLocal(position, thisOctal, subcell)
!          topOctal => thisOctal
!          topOctalSubcell = subcell
!          Do while(topOctal%changed(topOctalSubcell))
!             topOctalSubcell = topOctal%parentSubcell
!             topOctal => topOctal%parent
!          enddo
          call distanceToCellBoundary(grid, position, direction, tVal, thisOctal)
          rVec = subcellCentre(thisOctal, subcell)
          tempStorage(1) = thisOctal%rho(subcell)
!          tempStorage(2) = thisOctal%UVvector(subcell)%x
!          tempStorage(3) = thisOctal%UVvector(subcell)%y
!          tempStorage(4) = thisOctal%UVvector(subcell)%z
          tempStorage(2) = thisOctal%UVvectorPlus(subcell)%x
          tempStorage(3) = thisOctal%UVvectorPlus(subcell)%y
          tempStorage(4) = thisOctal%UVvectorPlus(subcell)%z
          tempStorage(5) = tval
          tempStorage(6) = thisOctal%temperature(subcell)
          tempStorage(7:39) = thisOctal%abundance(subcell, 1:33)
          tempStorage(40) = thisOctal%UVvectorMinus(subcell)%x
          tempStorage(41) = thisOctal%UVvectorMinus(subcell)%y
          tempStorage(42) = thisOctal%UVvectorMinus(subcell)%z
          call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, &
               localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine rayTracingServerPDR
#endif

  subroutine getRayTracingValuesXRAY(grid, position, direction, rho, tval, getTval)
!       radially, searchRadius)
    use mpi
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: rho
    type(VECTOR) :: position, rVec, direction
    real(double) :: loc(8)
    type(OCTAL), pointer :: thisOctal
    integer :: iThread
    integer, parameter :: nStorage = 2
    real(double) :: tempStorage(nStorage), tval
    integer :: subcell
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr
    logical, intent(in) :: getTval
!    logical :: useTop

    thisOctal => grid%octreeRoot
    call findSubcellLocal(position, thisOctal, subcell)
   
    if (octalOnThread(thisOctal, subcell, myrankGlobal)) then       

       rVec = subcellCentre(thisOctal, subcell)
       rho = thisOctal%rho(subcell)
       if(getTval) then
          call distanceToCellBoundary(grid, position, direction, tVal, thisOctal)
       else
          tVal = 1.d0
       endif
    else

       iThread = thisOctal%mpiThread(subcell)
       loc(1) = position%x
       loc(2) = position%y
       loc(3) = position%z
       loc(4) = dble(myRankGlobal)
       loc(5) = direction%x
       loc(6) = direction%y
       loc(7) = direction%z
       if(getTval) then
          loc(8) = 1.d0
       else
          loc(8) = 0.d0
       end if

       call MPI_SEND(loc, 8, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)

       rho = tempStorage(1)
       tval = tempstorage(2)
       if(.not.getTval) then
          tval = 1.d0
       endif

    endif
  end subroutine getRayTracingValuesXRAY


  subroutine rayTracingServerXRAY(grid)
    use mpi
    type(GRIDTYPE) :: grid
    logical :: stillServing
    real(double) :: loc(8)
    type(VECTOR) :: position, rVec, direction
    type(OCTAL), pointer :: thisOctal
!    type(OCTAL), pointer :: topOctal
    integer :: subcell!, nworking!, topOctalSubcell
    integer :: iThread!, servingArray, workingTHreads(nworking)
    integer, parameter :: nStorage = 2
    real(double) :: tempStorage(nStorage), tval
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr

    stillServing = .true.
 !   servingArray = 0
 !   workingTHreads = 0
    thisOctal => grid%octreeroot
    do while (stillServing)

!       do iThread = 1, nThreadsGlobal-1
!       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
       
       call MPI_RECV(loc, 8, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, localWorldCommunicator, status, ierr)
       position%x = loc(1)
       position%y = loc(2)
       position%z = loc(3)
       iThread = int(loc(4))
       direction%x = loc(5)
       direction%y = loc(6)
       direction%z = loc(7)
       if (position%x > 1.d29) then          
!          do j=1, nworking
!             if(workingThreads(j) == 0) then
!                workingThreads(j) = 1
!                exit
!             end if
!          end do
!          if(SUM(workingTHreads) == nworking) then
!             workingThreads = 0
             stillServing= .false.
!          end if
       else

          call findSubcellLocal(position, thisOctal, subcell)
!          topOctal => thisOctal
!          topOctalSubcell = subcell
!          do while(topOctal%changed(topOctalSubcell))
!             topOctalSubcell = topOctal%parentSubcell
!             topOctal => topOctal%parent
!          enddo
          if(loc(8) == 1.d0) then
             call distanceToCellBoundary(grid, position, direction, tVal, thisOctal)
             tempStorage(2) = tval
          else
             tempstorage(2) = 1.d0 
          endif
          rVec = subcellCentre(thisOctal, subcell)
          tempStorage(1) = thisOctal%rho(subcell)


          call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine rayTracingServerXRAY


  subroutine hydroValuesServer(grid, nworking)
    use mpi
    type(GRIDTYPE) :: grid
    logical :: stillServing
    real(double) :: loc(4)
    type(VECTOR) :: position, rVec
    type(OCTAL), pointer :: thisOctal
    type(OCTAL), pointer :: topOctal
    integer :: subcell, nworking, topOctalSubcell
    integer :: iThread, servingArray, workingTHreads(nworking)
    integer, parameter :: nStorage = 12
    real(double) :: tempStorage(nStorage)
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr, j

    stillServing = .true.
    servingArray = 0
    workingTHreads = 0
   thisOctal => grid%octreeroot
    do while (stillServing)

!       do iThread = 1, nThreadsGlobal-1
!       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
       
       call MPI_RECV(loc, 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, localWorldCommunicator, status, ierr)
       position%x = loc(1)
       position%y = loc(2)
       position%z = loc(3)
       iThread = int(loc(4))
       if (position%x > 1.d29) then          
          do j=1, nworking
             if(workingThreads(j) == 0) then
                workingThreads(j) = 1
                exit
             end if
          end do
          if(SUM(workingTHreads) == nworking) then
             workingThreads = 0
             stillServing= .false.
          end if
       else

          call findSubcellLocal(position, thisOctal, subcell)
          topOctal => thisOctal
          topOctalSubcell = subcell
          do while(topOctal%changed(topOctalSubcell))
             topOctalSubcell = topOctal%parentSubcell
             topOctal => topOctal%parent
          enddo

          rVec = subcellCentre(topOctal, topOctalsubcell)
          tempStorage(1) = thisOctal%nDepth
          tempStorage(2) = topOctal%rho(topOctalsubcell)
          tempStorage(3) = topOctal%rhoe(topOctalsubcell)             
          tempStorage(4) = topOctal%rhou(topOctalsubcell)             
          tempStorage(5) = topOctal%rhov(topOctalsubcell)             
          tempStorage(6) = topOctal%rhow(topOctalsubcell)        
          tempStorage(7) = topOctal%energy(topOctalsubcell)
          tempStorage(8) = topOctal%phi_gas(topOctalsubcell)
          tempStorage(9) = rVec%x
          tempStorage(10) = rVec%y
          tempStorage(11) = rVec%z
          tempstorage(12) = topOctal%pressure_i(topOctalsubcell)

          call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine hydroValuesServer

  subroutine hydroValuesServer2(grid, nworking, k, endloop, check)
    use mpi
    type(GRIDTYPE) :: grid
    logical :: stillServing
    real(double) :: loc(7)
    type(VECTOR) :: position, rVec
    type(OCTAL), pointer :: thisOctal
    type(OCTAL), pointer :: topOctal
    integer :: subcell, nworking, topOctalSubcell
    integer :: iThread, servingArray, workingTHreads(nworking)
    integer, parameter :: nStorage = 12
    real(double) :: tempStorage(nStorage)
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr, j, k, m, counter
    integer :: functionality
    integer :: endloop
    integer, parameter :: maxStorage = 1000
    real(double) :: storageArray(maxStorage,12), searchRadius
    real(double) :: tempStorageArray(maxStorage,12) 
    real(double) :: tempdoublearray(12*maxStorage)
    integer :: check(endloop, nHydroThreadsGlobal)
    integer :: nVals
    integer :: cellRef
!    integer :: foundID
    logical :: useTop
    real(double) :: topID

    stillServing = .true.
    servingArray = 0
    workingTHreads = 0
    thisOctal => grid%octreeroot

    do while (stillServing)
       storageArray = 0.d0
       tempStorageArray = 0.d0
       tempDoubleArray = 0.d0
       call MPI_RECV(loc, 7, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, localWorldCommunicator, status, ierr)
       position%x = loc(1)
       position%y = loc(2)
       position%z = loc(3)
       iThread = int(loc(4))
       functionality = int(loc(5))
       searchRadius = loc(6)
       topID = loc(7)

       if(loc(7) > 0.d0) then
          useTop = .true.
       else
          useTop = .false.
       end if


       if (position%x > 1.d29) then          
          do j=1, nworking
             if(workingThreads(j) == 0) then
                workingThreads(j) = 1
                exit
             end if
          end do
          if(SUM(workingTHreads) == nworking) then
             workingThreads = 0
             stillServing= .false.
          end if
          
       else if (functionality == 1) then
          !- this thread is going to host a gathering of values in the search radius
          loc(4) = myRankGlobal
          loc(5) = 2.d0
          do m = 1, nHydroThreadsGlobal
             if (m /= myrankGlobal .and. .not. ANY(m == check(k,1:nHydroThreadsGlobal))) then
                call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, m, tag, localWorldCommunicator, ierr)
             end if
          end do
          cellRef=0
          !Get the main serving thread's values within the request radius
          call getAllInRadius(grid%octreeRoot, position, searchRadius, storageArray, cellRef, useTop)
!          write(*,*) myrankGlobal, " called getallinradius and got ",cellref

          !send the order to the other serving threads to return values
          nvals = cellRef
          do m = 1, nHydroThreadsGlobal
             if (m /= myrankGlobal .and. .not. ANY(m == check(k,1:nHydroThreadsGlobal))) then 
                !Will recv a 1D array and translate it into 2D later
!                call MPI_RECV(tempStorageArray 12000, MPI_DOUBLE_PRECISION, m, tag, localWorldCommunicator, status, ierr)

                call MPI_RECV(cellRef, 1, MPI_INTEGER, m, tag, &
                     localWorldCommunicator, status, ierr)
                call MPI_RECV(tempDoubleArray, SIZE(tempDoubleArray), MPI_DOUBLE_PRECISION, m, tag, &
                     localWorldCommunicator, status, ierr)

                !Now convert it back to 2D
                tempStorageArray = reshape(tempDoubleArray, shape(tempStorageArray))


                storageArray((nVals+1):(nVals+cellRef),:) = tempStorageArray(1:cellRef,:)
                nVals = nVals + cellRef

                tempStorageArray = 0.d0
                tempDoubleArray = 0.d0
             end if
          end do

          !we should now have an array of values to send
          !first send the number of cells worth of data to expect
          call MPI_SEND(nvals, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr) 
          do counter = 1, nvals             
             call MPI_SEND(storageArray(counter,:), nStorage, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)             
          end do
          storageArray = 0.d0
       else if(functionality == 2)then
          !- responding to a request for values within the search radius
          cellRef = 0
          storageArray = 0.d0
          call getAllInRadius(grid%octreeRoot, position, searchRadius, storageArray, cellRef, useTop)
!          write(*,*) myrankGlobal, " called get all in radius ",cellRef
 
          !Don't want to send a 2D array
          tempDoubleArray = reshape(storageArray,(/SIZE(tempDoubleArray)/))
          call MPI_SEND(cellRef, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
          call MPI_SEND(tempDoubleArray, 12000, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
!          call MPI_SEND(storageArray, 1200, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)

       else
          call findSubcellLocal(position, thisOctal, subcell)
          if(useTop) then
             topOctal => thisOctal
             topOctalSubcell = subcell
             do while(topOctal%changed(topOctalSubcell))
                topOctalSubcell = topOctal%parentSubcell
                topOctal => topOctal%parent
             enddo
             rVec = subcellCentre(topOctal, topOctalsubcell)
             tempStorage(1) = topOctal%nDepth
             tempStorage(2) = topOctal%rho(topOctalsubcell)
             tempStorage(3) = topOctal%rhoe(topOctalsubcell)             
             tempStorage(4) = topOctal%rhou(topOctalsubcell)             
             tempStorage(5) = topOctal%rhov(topOctalsubcell)             
             tempStorage(6) = topOctal%rhow(topOctalsubcell)        
             tempStorage(7) = topOctal%energy(topOctalsubcell)
             tempStorage(8) = topOctal%phi_gas(topOctalsubcell)
             tempStorage(9) = rVec%x
             tempStorage(10) = rVec%y
             tempStorage(11) = rVec%z
             tempstorage(12) = topOctal%pressure_i(topOctalsubcell)
          else
             rVec = subcellCentre(thisOctal, subcell)
             tempStorage(1) = thisOctal%nDepth
             tempStorage(2) = thisOctal%rho(subcell)
             tempStorage(3) = thisOctal%rhoe(subcell)             
             tempStorage(4) = thisOctal%rhou(subcell)             
             tempStorage(5) = thisOctal%rhov(subcell)             
             tempStorage(6) = thisOctal%rhow(subcell)        
             tempStorage(7) = thisOctal%energy(subcell)
             tempStorage(8) = thisOctal%phi_gas(subcell)
             tempStorage(9) = rVec%x
             tempStorage(10) = rVec%y
             tempStorage(11) = rVec%z
             tempstorage(12) = thisOctal%pressure_i(subcell)
          end if
             


          call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
       endif
    enddo
  end subroutine hydroValuesServer2


  recursive  subroutine getAllInRadius(thisOctal, position, searchRadius, storageArray, cellRef, useTop)

    type(octal), pointer :: thisOctal, child, topOctal
    integer :: subcell, topOctalSubcell
    integer :: i
    type(VECTOR) :: position, rVec
    real(double) :: searchRadius
    integer, parameter :: maxStorage = 1000
    real(double) :: storageArray(maxStorage,12)
    integer :: cellRef, k
    logical :: changed, useTop, check

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getAllInRadius(child, position, searchRadius, storageArray, cellRef, useTop)
                exit
             end if
          end do
       else 
          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then

             changed = .false.
             topOctal => thisOctal
             topOctalSubcell = subcell
             do while(topOctal%changed(topOctalSubcell))
                topOctalSubcell = topOctal%parentSubcell
                topOctal => topOctal%parent
                changed = .true.
             enddo

             if(useTop) then
                rVec = subcellCentre(topOctal, topOctalsubcell)
             else
                rVec = subcellCentre(thisOCtal, subcell)
             end if

             if(modulus(rVec - position) < searchRadius) then
!                print *, "position ", position
!                print *, "rVec ", rVec
!                print *, "radius ", searchRadius
!                print *, "modulus(rVec - position)", modulus(rVec - position)
                if(useTop) then
                   check = .true.
                   do k = 1, cellRef
                      if(rVec%x == storageArray(k, 9) .and. rVec%y == storageArray(k, 10) &
                           .and. rVec%z == storageArray(k, 11)) then
                         check = .false.
                      end if
                   end do

                   if(check) then
                      cellRef = cellRef + 1
                      storageArray(cellRef, 1)  = topOctal%nDepth
                      storageArray(cellRef, 2)  = topOctal%rho(topOctalsubcell)
                      storageArray(cellRef, 3)  = topOctal%rhoe(topOctalsubcell)             
                      storageArray(cellRef, 4)  = topOctal%rhou(topOctalsubcell)             
                      storageArray(cellRef, 5)  = topOctal%rhov(topOctalsubcell)             
                      storageArray(cellRef, 6)  = topOctal%rhow(topOctalsubcell)        
                      storageArray(cellRef, 7)  = topOctal%energy(topOctalsubcell)
                      storageArray(cellRef, 8)  = topOctal%phi_gas(topOctalsubcell)
                      storageArray(cellRef, 9)  = rVec%x
                      storageArray(cellRef, 10)  = rVec%y
                      storageArray(cellRef, 11)  = rVec%z
                      storageArray(cellRef, 12)  = topOctal%pressure_i(topOctalsubcell)
                   end if
                else
                   cellRef = cellRef + 1
                   storageArray(cellRef, 1)  = thisOctal%nDepth
                   storageArray(cellRef, 2)  = thisOctal%rho(subcell)
                   storageArray(cellRef, 3)  = thisOctal%rhoe(subcell)             
                   storageArray(cellRef, 4)  = thisOctal%rhou(subcell)             
                   storageArray(cellRef, 5)  = thisOctal%rhov(subcell)             
                   storageArray(cellRef, 6)  = thisOctal%rhow(subcell)        
                   storageArray(cellRef, 7)  = thisOctal%energy(subcell)
                   storageArray(cellRef, 8)  = thisOctal%phi_gas(subcell)
                   storageArray(cellRef, 9)  = rVec%x
                   storageArray(cellRef, 10)  = rVec%y
                   storageArray(cellRef, 11)  = rVec%z
                   storageArray(cellRef, 12)  = thisOctal%pressure_i(subcell)
                end if
             end if
             if (changed) exit
          end if
       endif
    end do
end subroutine getAllInRadius
              
!  subroutine fillVelocityCornersFromHydro(grid)
!    use mpi
!    type(GRIDTYPE) :: grid
!    integer :: iThread
!    integer :: ierr
!
!       do iThread = 1, nThreadsGlobal-1
!          if (myrankGlobal /= iThread) then 
!             call hydroValuesServer(grid, iThread)
!          else
!             call  fillVelocityCornersFromCentres(grid%octreeRoot, grid)
!             call shutdownServers()
!          endif
!          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!       enddo
!
!     end subroutine fillVelocityCornersFromHydro
!
!  recursive subroutine fillVelocityCornersFromCentres(thisOctal, grid)
!    type(GRIDTYPE) :: grid
!    type(OCTAL), pointer :: thisOctal, child
!    integer :: i
!
!    if (thisOctal%nChildren > 0) then
!       do i = 1, thisOctal%nChildren
!          child => thisOctal%child(i)
!          call fillVelocityCornersFromCentres(child, grid)
!       enddo
!    else
!       call fillVelocityCorners(thisOctal, velocityCornerFixedGrid)
!    endif
!  end subroutine fillVelocityCornersFromCentres
!
!  function velocityCornerFixedGrid(rVec) result(meanVel)
!    type(OCTAL), pointer :: thisOctal, currentOctal
!    integer :: i, currentSubcell
!    type(VECTOR), intent(in) :: rVec
!    type(VECTOR) :: vel(8), cVec, pVec(8), meanVel
!    integer :: nd
!    real(double) :: rho, rhoe, rhou, rhov, rhow, energy, phi
!
!    pVec(1) = VECTOR(+0.1,+0.1,+0.1)
!    pVec(2) = VECTOR(-0.1,+0.1,+0.1)
!    pVec(3) = VECTOR(+0.1,-0.1,+0.1)
!    pVec(4) = VECTOR(-0.1,-0.1,+0.1)
!    pVec(5) = VECTOR(+0.1,+0.1,-0.1)
!    pVec(6) = VECTOR(-0.1,+0.1,-0.1)
!    pVec(7) = VECTOR(+0.1,-0.1,-0.1)
!    pVec(8) = VECTOR(-0.1,-0.1,-0.1)
!    currentOctal => grid%octreeRoot
!    call findSubcellLocal(rVec, currentOctal, currentSubcell)
!    thisOctal => currentOctal
!    meanVel = VECTOR(0.d0, 0.d0, 0.d0)
!    do i = 1, 8
!
!       cVec = rVec + grid%halfSmallestSubcell*pvec(i) 
!       if (inOctal(grid%octreeRoot, cVec)) then
!          call getHydroValues(grid, cVec, nd, rho, rhoe, rhou, rhov, rhow, energy, phi)
!          vel(i) = VECTOR(rhou/rho, rhov/rho, rhow/rho)/cSpeed
!       else
!          vel(i) = currentOctal%velocity(currentSubcell)
!       endif
!       meanVel = meanVel + vel(i)
!    enddo
!    meanVel  = meanVel / 8.d0
!  end function velocityCornerFixedGrid
!
#ifdef SPH
    subroutine distributeSphDataOverMPI()
      use sph_data_class, only : sphData, npart, init_sph_data
      use mpi
      integer :: ierr


      call MPI_BCAST(sphData%inUse, 1, MPI_LOGICAL, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%uDist, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%uMass, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%uTime, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%uVel , 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%uTemp, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%totalGasMass, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%codeVelocitytoTorus, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%codeEnergytoTemperature, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%npart, 1, MPI_INTEGER, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%time, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%nptmass, 1, MPI_INTEGER, 0, localWorldCommunicator, ierr)
      if (myrankGlobal /= 0) then
         npart = sphData%nPart
         call init_sph_data(sphdata%udist, sphdata%umass, sphdata%utime, sphdata%time, &
              sphdata%nptmass, sphdata%uvel, sphdata%utemp)
      endif

      call MPI_BCAST(sphData%gasMass, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

      call MPI_BCAST(sphData%xn, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%yn, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%zn, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

      call MPI_BCAST(sphData%vxn, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%vyn, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%vzn, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

      call MPI_BCAST(sphData%hn, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%rhon, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

!      call MPI_BCAST(sphData%rhoh2, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

      call MPI_BCAST(sphData%temperature, npart, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

      call MPI_BCAST(sphData%hpt, sphdata%nptmass, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

      call MPI_BCAST(sphData%x, sphdata%nptmass, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%y, sphdata%nptmass, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%z, sphdata%nptmass, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

      call MPI_BCAST(sphData%vx, sphdata%nptmass, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%vy, sphdata%nptmass, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
      call MPI_BCAST(sphData%vz, sphdata%nptmass, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

      call MPI_BCAST(sphData%ptmass, sphdata%nptmass, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
         
    end subroutine distributeSphDataOverMPI
#endif

!Check that appropriate number of threads have been used
    subroutine checkThreadNumber(grid)
      integer :: nDimensions
      type(gridtype) :: grid
      
      if(grid%octreeRoot%oneD) then
         nDimensions = 1
      else if(grid%octreeRoot%twoD) then
         nDimensions = 2
      else if(grid%octreeRoot%threeD) then
         nDimensions = 3
      else
         if(myRankWorldGlobal == 0) then
            write(*,*) "Can't comprehend more than 3, or zero, dimensions!"
         end if
         stop
      end if
      
      if((2**(nDimensions)) /= nHydroThreadsGlobal .and. &
         (4**(nDimensions)) /= nHydroThreadsGlobal .and. &
         (8**(nDimensions)) /= nHydroThreadsGlobal) then
         if(myRankWorldGlobal == 0) then
            write(*,*) "An incorrect number of threads has been used:"
            write(*,*) "For this model try: ", (2**(nDimensions) + 1), &
                 "or ", (4**(nDimensions) + 1)," or ", (8**(nDimensions) + 1), "threads"
         end if
         stop
      else
         if(myRankWorldGlobal == 0) then
            write(*,*) "Thread Check Complete"
         end if
      end if
    end subroutine checkThreadNumber


!When reading in a grid that uses a different number of mpi threads to that of the current run
!redistribute thread ID's to the appropriate cells
recursive subroutine distributeMPIthreadLabels(thisOctal)
!  use mpi_global_mod, only : myRankGlobal
  use mpi
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child
  integer :: subcell, i

!Go through grid and label MPI threads
!  write(*,*) "Rank ", myRankGlobal, "being allocated grid space"
  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        do i = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(i) == subcell) then
              call labelSingleSubcellMPI(thisOctal, subcell, i)
!              write(*,'(a,i3,a,i3,a,i3,a,i3,a,i3)'), &
!                   "rank ",myrankglobal," depth ",thisOctal%nDepth , " subcell ", subcell, &
!                   " thisOctal%mpiThread(subcell) ", thisOctal%mpiThread(subcell)
              child => thisoctal%child(i)
              call distributeMPIthreadLabels(child)
              exit
           end if
        end do
     end if
  end do

end subroutine distributeMPIthreadLabels


!Label an individual subcell with its MPI thread
subroutine labelSingleSubcellMPI(parent, iChild, newChildIndex)
  use mpi_global_mod, only: nHydroThreadsGlobal
  type(octal), pointer   :: parent
  integer ::  i, iChild, newChildIndex

    if ( ((parent%twoD)  .and.((nHydroThreadsGlobal) == 4)) .or. &
         ((parent%threed).and.((nHydroThreadsGlobal) == 8)).or. &
         ((parent%oneD)  .and.((nHydroThreadsGlobal) == 2)) ) then
       parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)

    else
       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          if (parent%oneD) then
             do i = 1, 2
                parent%child(newChildIndex)%mpiThread(i) = 2 * (parent%mpiThread(iChild) - 1) + i
             enddo
          else if (parent%twoD) then
             do i = 1, 4
                parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
             enddo
          else if (parent%Threed) then
             do i = 1, 8
                parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
             enddo
          endif
       endif
    endif
    if ((parent%twod).and.(nHydroThreadsGlobal) == 16) then
       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 4
             parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    if ((parent%threed).and.(nHydroThreadsGlobal) == 64) then
       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 8
             parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    if ((parent%threed).and.(nHydroThreadsGlobal) == 512) then
       if (parent%child(newChildIndex)%nDepth > 3) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 8
             parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    if ((parent%twoD).and.(nHydroThreadsGlobal) == 64) then
       if (parent%child(newChildIndex)%nDepth > 3) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 4
             parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

end subroutine labelSingleSubcellMPI


function shepardsMethod(xi, yi, zi, fi, n, x, y, z) result(out)
  
  real(double) :: xi(:), yi(:), zi(:), fi(:)
  real(double) :: x, y, z, out
  integer :: n, i
  real(double), allocatable :: w(:), h(:)
  real(double) :: R


  allocate(w(1:n), h(1:n))
      
  do i = 1, n
         h(i) = max(1.d-3,sqrt( (x-xi(i))**2 + (y-yi(i))**2 + (z-zi(i))**2))
      enddo
      R = maxval(h(1:n))
      do i = 1, n
         w(i) = ((R-h(i))/(R*h(i))**3) ! franke & nielson 1980
      enddo
      w = w / SUM(w)
      out = 0.d0
      do i = 1, n
         out = out + w(i) * fi(i)
      enddo
      deallocate(w, h)
      
    end function shepardsMethod


    
      


  recursive subroutine packBranch(thisOctal)
    type(OCTAL), pointer :: thisOctal, child
    integer :: i

    call packOctal(thisOctal)


    IF ( thisOctal%nChildren > 0 ) THEN
      DO i = 1, thisOctal%nChildren, 1
         ! call this subroutine recursively on each of its children
        child => thisOctal%child(i)
        CALL packBranch(child)
      END DO
    END IF
  end subroutine packBranch

  recursive subroutine unpackBranch(thisOctal, parent)
    type(OCTAL), pointer :: thisOctal, child, parent
    integer :: i

    call unpackOctal(thisOctal)



    if (thisOctal%nChildren > 0) allocate(thisOctal%child(1:thisOctal%nChildren))
    thisOctal%parent => parent

    IF ( thisOctal%nChildren > 0 ) THEN
      DO i = 1, thisOctal%nChildren, 1
         ! call this subroutine recursively on each of its children
        child => thisOctal%child(i)
        CALL unpackBranch(child, thisOctal)
      END DO
    END IF
  end subroutine unpackBranch
  

  subroutine packOctal(thisOctal)
    type(OCTAL), pointer :: thisOctal

    call packAttributeStatic(thisOctal%nDepth)
    call packAttributeStatic(thisOctal%maxChildren)
    call packAttributeStatic(thisOctal%mpiThread)
    call packAttributeStatic(thisOctal%nChildren)
    call packAttributeStatic(thisOctal%oneD)
    call packAttributeStatic(thisOctal%twoD)
    call packAttributeStatic(thisOctal%threeD)
    call packAttributeStatic(thisOctal%cylindrical)
    call packAttributeStatic(thisOctal%splitAzimuthally)
    call packAttributeStatic(thisOctal%hasChild)
    call packAttributeStatic(thisOctal%indexChild)
    call packAttributeStatic(thisOctal%parentSubcell)
    call packAttributeStatic(thisOctal%centre)
    call packAttributeStatic(thisOctal%xMax)
    call packAttributeStatic(thisOctal%xMin)
    call packAttributeStatic(thisOctal%yMax)
    call packAttributeStatic(thisOctal%yMin)
    call packAttributeStatic(thisOctal%zMax)
    call packAttributeStatic(thisOctal%zMin)
    call packAttributeStatic(thisOctal%subcellSize)
    call packAttributeStatic(thisOctal%rho)
    call packAttributeStatic(thisOctal%temperature)
    call packAttributeStatic(thisOctal%inflow)


    call packAttributePointer(thisOctal%dustTypeFraction)
    call packAttributePointer(thisOctal%diffusionApprox)
    call packAttributePointer(thisOctal%changed)
    call packAttributePointer(thisOctal%nDiffusion)
    call packAttributePointer(thisOctal%underSampled)
    call packAttributePointer(thisOctal%distanceGrid)
    call packAttributePointer(thisOctal%nCrossings)
    call packAttributePointer(thisOctal%nTot)
    call packAttributePointer(thisOctal%chiLine)
    call packAttributePointer(thisOctal%biasLine3D)
    call packAttributePointer(thisOctal%etaCont)
    call packAttributePointer(thisOctal%nh)
    call packAttributePointer(thisOctal%ne)
    call packAttributePointer(thisOctal%HHeating)
    call packAttributePointer(thisOctal%HeHeating)
    call packAttributePointer(thisOctal%radiationMomentum)
    call packAttributePointer(thisOctal%ionFrac)
    call packAttributePointer(thisOctal%photoionCoeff)
    call packAttributePointer(thisOctal%sourceContribution)
    call packAttributePointer(thisOctal%diffuseContribution)
    call packAttributePointer(thisOctal%normSourceContribution)
    call packAttributePointer(thisOctal%kappaTimesFlux)
    call packAttributePointer(thisOctal%uvvector)
    call packAttributePointer(thisOctal%oldfrac)

  end subroutine packOctal

  subroutine unpackOctal(thisOctal)
    type(OCTAL), pointer :: thisOctal

    call unpackAttributeStatic(thisOctal%nDepth)
    call unpackAttributeStatic(thisOctal%maxChildren)
    call unpackAttributeStatic(thisOctal%mpiThread)
    call unpackAttributeStatic(thisOctal%nChildren)
    call unpackAttributeStatic(thisOctal%oneD)
    call unpackAttributeStatic(thisOctal%twoD)
    call unpackAttributeStatic(thisOctal%threeD)
    call unpackAttributeStatic(thisOctal%cylindrical)
    call unpackAttributeStatic(thisOctal%splitAzimuthally)
    call unpackAttributeStatic(thisOctal%hasChild)
    call unpackAttributeStatic(thisOctal%indexChild)
    call unpackAttributeStatic(thisOctal%parentSubcell)
    call unpackAttributeStatic(thisOctal%centre)
    call unpackAttributeStatic(thisOctal%xMax)
    call unpackAttributeStatic(thisOctal%xMin)
    call unpackAttributeStatic(thisOctal%yMax)
    call unpackAttributeStatic(thisOctal%yMin)
    call unpackAttributeStatic(thisOctal%zMax)
    call unpackAttributeStatic(thisOctal%zMin)
    call unpackAttributeStatic(thisOctal%subcellSize)
    call unpackAttributeStatic(thisOctal%rho)
    call unpackAttributeStatic(thisOctal%temperature)
    call unpackAttributeStatic(thisOctal%inflow)

    call unpackAttributePointer(thisOctal%dustTypeFraction)
    call unpackAttributePointer(thisOctal%diffusionApprox)
    call unpackAttributePointer(thisOctal%changed)
    call unpackAttributePointer(thisOctal%nDiffusion)
    call unpackAttributePointer(thisOctal%underSampled)
    call unpackAttributePointer(thisOctal%distanceGrid)
    call unpackAttributePointer(thisOctal%nCrossings)
    call unpackAttributePointer(thisOctal%nTot)
    call unpackAttributePointer(thisOctal%chiLine)
    call unpackAttributePointer(thisOctal%biasLine3D)
    call unpackAttributePointer(thisOctal%etaCont)
    call unpackAttributePointer(thisOctal%nh)
    call unpackAttributePointer(thisOctal%ne)
    call unpackAttributePointer(thisOctal%HHeating)
    call unpackAttributePointer(thisOctal%HeHeating)
    call unpackAttributePointer(thisOctal%radiationMomentum)
    call unpackAttributePointer(thisOctal%ionFrac)
    call unpackAttributePointer(thisOctal%photoionCoeff)
    call unpackAttributePointer(thisOctal%sourceContribution)
    call unpackAttributePointer(thisOctal%diffuseContribution)
    call unpackAttributePointer(thisOctal%normSourceContribution)
    call unpackAttributePointer(thisOctal%kappaTimesFlux)
    call unpackAttributePointer(thisOctal%uvvector)
    call unpackAttributePointer(thisOctal%oldfrac)


  end subroutine unPackOctal


  subroutine broadcastBranch(grid, fromThread, communicator)
    use timing
    use mpi
    integer :: fromThread, nOctals
    type(GRIDTYPE) :: grid
    integer :: communicator
    integer :: ierr, i

    if (allocated(buffer)) deallocate(buffer)
    maxBuffer = 0
    if (myrankGlobal == fromThread) then
       call countVoxels(grid%octreeRoot, nOctals, maxBuffer)
       maxBuffer = maxBuffer * 1000
    endif
    call MPI_BCAST(maxBuffer, 1, MPI_INTEGER, 0, communicator, ierr)
    allocate(buffer(1:maxBuffer))
    nBuffer = 1

    if (myrankGlobal == fromThread) then
!       CALL checkAMRgrid(grid,checkNoctals=.FALSE.)
       call packBranch(grid%octreeRoot)
    endif

! NB fromThread is the zeroth rank of the communicator

    call MPI_BCAST(nBuffer, 1, MPI_INTEGER, 0, communicator, ierr)
    call MPI_BCAST(buffer(1:nBuffer), nBuffer, MPI_DOUBLE_PRECISION, 0, communicator, ierr)

    if (myrankGlobal /= fromThread) then

       call deleteOctreeBranch(thisOctal=grid%octreeroot,                      &
               onlyChildren=.FALSE.,                                           &
               adjustParent=.FALSE.)

       nBuffer = 1

       call unpackbranch(grid%octreeRoot, null())

       CALL updateMaxDepth(grid)
       CALL setSmallestSubcell(grid)
       call countVoxels(grid%octreeRoot, grid%nOctals, i)

!       CALL checkAMRgrid(grid,checkNoctals=.FALSE.)
       copyOfThread = fromThread
    endif
    call MPI_BARRIER(communicator,ierr)

  end subroutine broadcastBranch

  subroutine packAttributeDoubleArray2Dpointer(array)
    real(double) :: array(:,:)
    integer :: m, n, dim(1)
    m = SIZE(array,1)
    n = SIZE(array,2)
    buffer(nBuffer) = dble(m)
    nBuffer = nBuffer + 1
    buffer(nBuffer) = dble(n)
    nBuffer = nBuffer + 1
    dim(1) = n*m
    buffer(nBuffer:(nBuffer+m*n-1)) = reshape(array, dim)
    nBuffer = nBuffer + m*n
  end subroutine packAttributeDoubleArray2Dpointer

  subroutine packAttributeDoubleArray3Dpointer(array)
    real(double) :: array(:,:,:)
    integer :: m, n, o, dim(1)
    m = SIZE(array,1)
    n = SIZE(array,2)
    o = SIZE(array,3)
    buffer(nBuffer) = dble(m)
    nBuffer = nBuffer + 1
    buffer(nBuffer) = dble(n)
    nBuffer = nBuffer + 1
    buffer(nBuffer) = dble(o)
    nBuffer = nBuffer + 1
    dim(1) = n*m*o
    buffer(nBuffer:(nBuffer+m*n*o-1)) = reshape(array, dim)
    nBuffer = nBuffer + m*n*o
  end subroutine packAttributeDoubleArray3Dpointer

  subroutine packAttributeDoubleArray1Dpointer(array)
    real(double) :: array(:)
    integer :: m
    m = SIZE(array)
    buffer(nBuffer) = dble(m)
    nBuffer = nBuffer + 1
    buffer(nBuffer:(nBuffer+m-1)) = array
    nBuffer = nBuffer + m
  end subroutine packAttributeDoubleArray1Dpointer

  subroutine packAttributeRealArray1Dpointer(array)
    real :: array(:)
    integer :: m
    m = SIZE(array)
    buffer(nBuffer) = dble(m)
    nBuffer = nBuffer + 1
    buffer(nBuffer:(nBuffer+m-1)) = real(array)
    nBuffer = nBuffer + m
  end subroutine packAttributeRealArray1Dpointer

  subroutine packAttributeLogicalArray1Dpointer(array)
    logical :: array(:)
    integer :: m
    m = SIZE(array)
    buffer(nBuffer) = dble(m)
    nBuffer = nBuffer + 1
    buffer(nBuffer:(nBuffer+m-1)) = 0.d0
    where(array(1:m))
       buffer(nBuffer:(nBuffer+m-1)) = 1.d0
    end where

    nBuffer = nBuffer + m
  end subroutine packAttributeLogicalArray1Dpointer



  subroutine packAttributeVectorArray1Dpointer(array)
    type(VECTOR) :: array(:)
    integer :: m
    m = SIZE(array)
    buffer(nBuffer) = dble(m)
    nBuffer = nBuffer + 1

    buffer(nBuffer:(nBuffer+m-1)) = array(1:m)%x
    nBuffer = nBuffer + m
    buffer(nBuffer:(nBuffer+m-1)) = array(1:m)%y
    nBuffer = nBuffer + m
    buffer(nBuffer:(nBuffer+m-1)) = array(1:m)%z
    nBuffer = nBuffer + m
  end subroutine packAttributeVectorArray1Dpointer

  subroutine packAttributeIntegerArray1Dpointer(array)
    integer :: array(:)
    integer :: m
    m = SIZE(array)
    buffer(nBuffer) = dble(m)
    nBuffer = nBuffer + 1
    buffer(nBuffer:(nBuffer+m-1)) = dble(array)
    nBuffer = nBuffer + m
  end subroutine packAttributeIntegerArray1Dpointer

  subroutine unpackAttributeDoubleArray2Dpointer(array)
    real(double), pointer :: array(:,:)
    integer :: m, n, dim(2)
    m = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    n = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    allocate(array(1:m,1:n))
    dim(1) = m
    dim(2) = n
    array = reshape(buffer(nBuffer:(nBuffer+m*n-1)),dim)
    nBuffer = nBuffer + m*n
  end subroutine unpackAttributeDoubleArray2Dpointer

  subroutine unpackAttributeDoubleArray3Dpointer(array)
    real(double), pointer :: array(:,:, :)
    integer :: m, n, o, dim(3)
    m = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    n = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    o = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    allocate(array(1:m,1:n,1:o))
    dim(1) = m
    dim(2) = n
    dim(3) = o
    array = reshape(buffer(nBuffer:(nBuffer+m*n*o-1)),dim)
    nBuffer = nBuffer + m*n*o
  end subroutine unpackAttributeDoubleArray3Dpointer

  subroutine unpackAttributeDoubleArray1Dpointer(array)
    real(double), pointer :: array(:)
    integer :: m
    m = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    allocate(array(1:m))
    array = buffer(nBuffer:(nBuffer+m-1))
    nBuffer = nBuffer + m
  end subroutine unpackAttributeDoubleArray1Dpointer

  subroutine unpackAttributeRealArray1Dpointer(array)
    real, pointer :: array(:)
    integer :: m
    m = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    allocate(array(1:m))
    array = real(buffer(nBuffer:(nBuffer+m-1)))
    nBuffer = nBuffer + m
  end subroutine unpackAttributeRealArray1Dpointer

  subroutine unpackAttributeLogicalArray1Dpointer(array)
    logical, pointer :: array(:)
    integer :: m
    m = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    allocate(array(1:m))
    array = .false.
    where (buffer(nBuffer:(nBuffer+m-1)) > 0.d0)
       array = .true.
    end where
    nBuffer = nBuffer + m
  end subroutine unpackAttributeLogicalArray1Dpointer

  subroutine unpackAttributeVectorArray1Dpointer(array)
    type(VECTOR), pointer :: array(:)
    integer :: m
    m = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    allocate(array(1:m))
    array(1:m)%x = buffer(nBuffer:(nBuffer+m-1))
    nBuffer = nBuffer + m
    array(1:m)%y = buffer(nBuffer:(nBuffer+m-1))
    nBuffer = nBuffer + m
    array(1:m)%z = buffer(nBuffer:(nBuffer+m-1))
    nBuffer = nBuffer + m
  end subroutine unpackAttributeVectorArray1Dpointer

  subroutine unpackAttributeIntegerArray1Dpointer(array)
    integer, pointer :: array(:)
    integer :: m
    m = nint(buffer(nBuffer))
    nBuffer = nBuffer + 1
    allocate(array(1:m))
    array = nint(buffer(nBuffer:(nBuffer+m-1)))
    nBuffer = nBuffer + m
  end subroutine unpackAttributeIntegerArray1Dpointer




  subroutine packAttributeIntegerSingle(i)
    integer :: i
    buffer(nBuffer) = dble(i)
    nBuffer = nBuffer + 1
  end subroutine packAttributeIntegerSingle

  subroutine packAttributeDoubleSingle(i)
    real(double) :: i
    buffer(nBuffer) = i
    nBuffer = nBuffer + 1
  end subroutine packAttributeDoubleSingle

  subroutine packAttributeVectorSingle(v)
    type(VECTOR) :: v
    buffer(nBuffer) = v%x
    nBuffer = nBuffer + 1
    buffer(nBuffer) = v%y
    nBuffer = nBuffer + 1
    buffer(nBuffer) = v%z
    nBuffer = nBuffer + 1
  end subroutine packAttributeVectorSingle

  subroutine unpackAttributeVectorSingle(v)
    type(VECTOR) :: v
    v%x = buffer(nBuffer)
    nBuffer = nBuffer + 1
    v%y = buffer(nBuffer)
    nBuffer = nBuffer + 1
    v%z = buffer(nBuffer)
    nBuffer = nBuffer + 1
  end subroutine unpackAttributeVectorSingle

  subroutine packAttributeLogicalSingle(i)
    logical :: i
    if (i) then
       buffer(nBuffer) = 1.d0
    else
       buffer(nBuffer) = 0.d0
    endif
    nBuffer = nBuffer + 1
  end subroutine packAttributeLogicalSingle

  subroutine packAttributeIntegerArray(array)
    integer :: array(:)
    integer :: n
    n = size(array)
    buffer(nBuffer:(nBuffer+n-1)) = dble(array(1:n))
    nBuffer = nBuffer + n
  end subroutine packAttributeIntegerArray

  subroutine packAttributeLogicalArray(array)
    logical :: array(:)
    integer :: i(8)
    integer :: n
    n = size(array)
    i = 0 
    where(array)
       i = 1
    end where
    buffer(nBuffer:(nBuffer+n-1)) = dble(i(1:n))
    nBuffer = nBuffer + n
  end subroutine packAttributeLogicalArray

  subroutine packAttributeDoubleArray(array)
    real(double) :: array(:)
    integer :: n
    n = size(array)
    buffer(nBuffer:(nBuffer+n-1)) = array(1:n)
    nBuffer = nBuffer + n
  end subroutine packAttributeDoubleArray

  subroutine packAttributeRealArray(array)
    real :: array(:)
    integer :: n
    n = size(array)
    buffer(nBuffer:(nBuffer+n-1)) = dble(array(1:n))
    nBuffer = nBuffer + n
  end subroutine packAttributeRealArray


  subroutine unpackAttributeIntegerSingle(i)
    integer :: i
    i = nint(buffer(nbuffer))
    nBuffer = nBuffer + 1
  end subroutine unpackAttributeIntegerSingle


  subroutine unpackAttributeDoubleSingle(i)
    real(double) :: i
    i = buffer(nBuffer)
    nBuffer = nBuffer + 1
  end subroutine unpackAttributeDoubleSingle


  subroutine unpackAttributeLogicalSingle(i)
    logical:: i
    if (buffer(nBuffer) == 0.d0) then
       i = .false.
    else
       i = .true.
    endif
    nBuffer = nBuffer + 1
  end subroutine unpackAttributeLogicalSingle

  subroutine unpackAttributeIntegerArray(array)
    integer :: array(:)
    integer :: n
    n = size(array)
    array(1:n) = nint(buffer(nBuffer:(nBuffer+n-1)))
    nBuffer = nBuffer + n
  end subroutine unpackAttributeIntegerArray

  subroutine unpackAttributeDoubleArray(array)
    real(double) :: array(:)
    integer :: n
    n = size(array)
    array(1:n) = buffer(nBuffer:(nBuffer+n-1))
    nBuffer = nBuffer + n
  end subroutine unpackAttributeDoubleArray

  subroutine unpackAttributeRealArray(array)
    real :: array(:)
    integer :: n
    n = size(array)
    array(1:n) = real(buffer(nBuffer:(nBuffer+n-1)))
    nBuffer = nBuffer + n
  end subroutine unpackAttributeRealArray

  subroutine unpackAttributeLogicalArray(array)
    logical :: array(:)
    integer :: n
    n = size(array)
    array = .false.
    where(buffer(nBuffer:(nBuffer+n-1)) > 0.d0)
       array = .true.
    end where
    nBuffer = nBuffer + n
  end subroutine unpackAttributeLogicalArray


  logical function checkEven(grid)
    use mpi
    type(GRIDTYPE) :: grid
    logical, allocatable :: check(:), tc(:)
    integer :: ierr
    allocate(check(1:nHydroThreadsGlobal), tc(1:nHydroThreadsGlobal))

    check = .true.
    call recurCheckEven(grid, grid%octreeRoot, check(myRankGlobal))
    call MPI_ALLREDUCE(check, tc, nHydrothreadsGlobal, MPI_LOGICAL, MPI_LAND, amrCommunicator, ierr)
    check = tc
    checkEven = .false.
    if (ALL(check(1:nHydroThreadsGlobal))) checkEven = .true.
    deallocate(check, tc)
  end function checkEven


recursive subroutine recurCheckEven(grid, thisOctal, check)
  use inputs_mod, only : smallestCellSize
  type(GRIDTYPE) :: grid
  type(OCTAL), pointer :: thisOctal, neighbourOctal, child
  integer :: i, subcell, neighbourSubcell
  logical :: check
  type(VECTOR) :: locator,direction(6)
  integer :: nDir, nd, idir, nc
  real(double) :: q, rho, rhoe, rhov, rhou,rhow,x,qnext,pressure,&
       flux,phi,phigas,xnext,px,py,pz,rm1,um1,pm1,qviscosity(3,3)
  
  if (thisOctal%threed) then
     ndir = 6
  else if (thisOctal%twod) then
     ndir = 4
  else 
     ndir = 2
  endif
  direction(1) = VECTOR(+1.d0, 0.d0,  0.d0)
  direction(2) = VECTOR(-1.d0, 0.d0,  0.d0)
  direction(3) = VECTOR( 0.d0, 0.d0, -1.d0)
  direction(4) = VECTOR( 0.d0, 0.d0, +1.d0)
  direction(5) = VECTOR( 0.d0,+1.d0,  0.d0)
  direction(6) = VECTOR( 0.d0,-1.d0,  0.d0)


  if (thisOctal%nChildren > 0) then
     do i = 1, thisOctal%nChildren, 1
        child => thisOctal%child(i)
        call recurCheckEven(grid, child, check)
     end do
  else
     
     do subcell = 1, thisOctal%maxChildren
        
        if (.not.octalOnThread(thisOctal, subcell,myrankGlobal)) cycle

        do idir = 1, nDir
          locator = subcellcentre(thisoctal, subcell) + direction(idir) * (thisoctal%subcellsize/2.d0+0.01d0*smallestCellsize)
          if (inOctal(grid%octreeRoot, locator)) then
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, direction(idir), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, nc, xnext, px, py, pz, rm1,um1, pm1, qViscosity)
             if (abs(thisOctal%nDepth-nd) > 1) check = .false.
          endif

        enddo
     enddo
  endif
end subroutine recurCheckEven


#endif

  end module mpi_amr_mod
