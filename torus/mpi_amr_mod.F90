module mpi_amr_mod
#ifdef MPI

  use kind_mod
  use amr_mod
  use mpi_global_mod
  use ion_mod

  implicit none

contains

  subroutine setupAMRCOMMUNICATOR
    include 'mpif.h'
    integer :: ierr, ranks(1)
    integer :: worldGroup, amrGroup

    call MPI_COMM_GROUP(MPI_COMM_WORLD, worldGroup, ierr)
    ranks(1) = 0
    call MPI_GROUP_EXCL(worldGroup, 1, ranks, amrGroup, ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD, amrGroup, amrCOMMUNICATOR, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreadsGlobal, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRankGlobal, ierr)

    if (nThreadsGlobal == 8) nHydroThreadsGlobal = 8
    if (nThreadsGlobal == 9) nHydroThreadsGlobal = 8

    if (nThreadsGlobal == 64) nHydroThreadsGlobal = 64
    if (nThreadsGlobal == 65) nHydroThreadsGlobal = 64

  end subroutine setupAMRCOMMUNICATOR

  subroutine findMassOverAllThreads(grid, mass)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: mass
    real(double), allocatable :: massOnThreads(:), temp(:)
    integer :: ierr

    allocate(massOnThreads(1:nThreadsGlobal), temp(1:nThreadsGlobal))
    massOnThreads = 0.d0
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
    end if
666 continue
  end subroutine findMassOverAllThreads
    
  recursive subroutine findTotalMassMPI(thisOctal, totalMass, minRho, maxRho)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: totalMass
  real(double),optional :: minRho, maxRho
  real(double) :: dv
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

          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
             dv = cellVolume(thisOctal, subcell)
             totalMass = totalMass + (1.d30)*thisOctal%rho(subcell) * dv
             if (PRESENT(minRho)) then
                if (thisOctal%rho(subcell) > 1.d-20) then
                   minRho = min(thisOctal%rho(subcell), minRho)
                endif
             endif
             if (PRESENT(maxRho)) maxRho = max(dble(thisOctal%rho(subcell)), maxRho)
          endif
       endif
    enddo
  end subroutine findTotalMassMPI


  

  subroutine getSquares(grid, plane, valueName, nSquares, corners, value, speed, ang)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    character(len=*) :: plane, valueName
    real, pointer :: corners(:,:), value(:), speed(:), ang(:)
    integer :: myRank, ierr, nThreads, iThread
    integer :: nSquares, n, tag=97,i
    integer :: status(MPI_STATUS_SIZE)
    integer :: j

    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)

    if (myRank == 1) then
       nSquares = 0
       call recursGetSquares(grid%octreeRoot, grid, plane, valueName, nSquares, corners, value, speed, ang)
       do iThread = 2, nThreads - 1
          call MPI_RECV(n, 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
          do i = 1, 4
             call MPI_RECV(corners(nSquares+1:nSquares+n, i), n, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
          enddo
          call MPI_RECV(value(nSquares+1:nSquares+n), n, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
          call MPI_RECV(speed(nSquares+1:nSquares+n), n, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
          call MPI_RECV(ang(nSquares+1:nSquares+n), n, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
          nSquares = nSquares+n
       end do
    else
       nSquares = 0
       call recursGetSquares(grid%octreeRoot, grid, plane, valueName, nSquares, corners, value, speed, ang)
       call MPI_SEND(nSquares, 1, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, ierr)
       do j = 1, 4
          call MPI_SEND(corners(1:nSquares,j), nSquares, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
       enddo
       call MPI_SEND(value(1:nSquares), nSquares, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
       call MPI_SEND(speed(1:nSquares), nSquares, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
       call MPI_SEND(ang(1:nSquares), nSquares, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
    endif
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
  end subroutine getSquares


  recursive subroutine recursGetSquares(thisOctal, grid, plane, valueName, nSquares, corners, value, speed, ang)

    include 'mpif.h'
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    !
    integer :: subcell, i
    character(len=*) :: valueName, plane
    integer :: nSquares
    real, pointer :: corners(:,:)
    real, pointer :: value(:)
    real, pointer :: speed(:)
    real, pointer :: ang(:)
    real ::tmp
    real(double) :: eps, x
    integer :: myRank, ierr
    type(VECTOR) :: rVec

    eps = 0.01d0 * grid%halfSmallestSubcell

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call recursGetSquares(child, grid, plane, valueName, nSquares, corners, value, speed, ang)
                exit
             end if
          end do
       else

!          if ((grid%splitOverMpi).and.(thisOctal%mpiThread(subcell) /= myRank)) then
!             cycle
!          endif
          if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle

          rVec = subcellCentre(thisOctal, subcell)
          select case(valuename)
             case("rho")
                tmp = thisOctal%rho(subcell)
             case("cs")
                tmp = sqrt(thisOctal%pressure_i(subcell)/thisOctal%rho(subcell))
             case("mass")
                tmp = thisOctal%rho(subcell) * cellVolume(thisOctal, subcell) * 1.d30 / msol
             case("pressure")
                tmp = thisOctal%pressure_i(subcell)
             case("u")
                tmp = thisOctal%rhou(subcell)/thisOctal%rho(subcell)/1.d5
             case("v")
                tmp = thisOctal%rhov(subcell)/thisOctal%rho(subcell)/1.d5
             case("w")
                tmp = thisOctal%rhow(subcell)/thisOctal%rho(subcell)/1.d5
             case("rhou")
                tmp = thisOctal%rhou(subcell)
             case("rhov")
                tmp = thisOctal%rhov(subcell)
             case("rhow")
                tmp = thisOctal%rhow(subcell)
             case("rhoe")
                tmp = thisOctal%rhoe(subcell)
             case("ionization")
                tmp = max(1.d-20,thisOctal%ionFrac(subcell,1))
             case("coeff")
                tmp = max(1.d-20,thisOctal%photoionCoeff(subcell,1))
             case("crossings")
                tmp = thisOctal%nCrossings(subcell)
             case("temperature")
                tmp = thisOctal%temperature(subcell)
             case("mpi")
                tmp = real(thisOctal%mpiThread(subcell))
             case("phi")
                tmp = real(thisOctal%phi_i(subcell))
             case("chi")
                tmp = real(thisOctal%chiline(subcell))
             case("u_i")
                tmp = real(thisOctal%u_interface(subcell))
             case("q_i")
                tmp = real(thisOctal%q_i(subcell))
             case("q_i_minus_1")
                tmp = real(thisOctal%q_i_minus_1(subcell))
             case("flux_i")
                tmp = real(thisOctal%flux_i(subcell))
             case("philim")
                tmp = real(thisOctal%phiLimit(subcell))
             case DEFAULT
           end select
           if (thisOctal%oneD) then
              nSquares = nSquares + 1
              corners(nsquares, 1) = rVec%x
              value(nSquares) = tmp
           else

              select case(plane)
              case("x-z")
                 x = min(abs(rVec%y + thisOctal%subcellSize/2.d0), abs(rVec%y - thisOctal%subcellSize/2.d0))
                 if ((x < eps).or.(thisOctal%twoD)) then
                    nSquares = nSquares + 1
                    corners(nSquares, 1) = rVec%x - thisOctal%subcellSize/2.d0
                    corners(nSquares, 2) = rVec%x + thisOctal%subcellSize/2.d0
                    corners(nSquares, 3) = rVec%z - thisOctal%subcellSize/2.d0
                    corners(nSquares, 4) = rVec%z + thisOctal%subcellSize/2.d0
                    value(nSquares) = tmp
                    speed(nSquares) = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhow(subcell)**2) / thisOctal%rho(subcell)
                    ang(nSquares) = atan2(thisOctal%rhow(subcell), thisOctal%rhou(subcell))
                 endif
              case("y-z")
                 x = min(abs(rVec%x + thisOctal%subcellSize/2.d0), abs(rVec%x - thisOctal%subcellSize/2.d0))
                 if (x < eps) then
                    nSquares = nSquares + 1
                    corners(nSquares, 1) = rVec%y - thisOctal%subcellSize/2.d0
                    corners(nSquares, 2) = rVec%y + thisOctal%subcellSize/2.d0
                    corners(nSquares, 3) = rVec%z - thisOctal%subcellSize/2.d0
                    corners(nSquares, 4) = rVec%z + thisOctal%subcellSize/2.d0
                    value(nSquares) = tmp
                    speed(nSquares) = sqrt(thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2) / thisOctal%rho(subcell)
                    ang(nSquares) = atan2(thisOctal%rhow(subcell), thisOctal%rhov(subcell))
                 endif
              case("x-y")
                 x = min(abs(rVec%z + thisOctal%subcellSize/2.d0), abs(rVec%z - thisOctal%subcellSize/2.d0))
                 if (x < eps) then
                    nSquares = nSquares + 1
                    corners(nSquares, 1) = rVec%x - thisOctal%subcellSize/2.d0
                    corners(nSquares, 2) = rVec%x + thisOctal%subcellSize/2.d0
                    corners(nSquares, 3) = rVec%y - thisOctal%subcellSize/2.d0
                    corners(nSquares, 4) = rVec%y + thisOctal%subcellSize/2.d0
                    value(nSquares) = tmp
                    speed(nSquares) = sqrt(thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2) / thisOctal%rho(subcell)
                    ang(nSquares) = atan2(thisOctal%rhov(subcell), thisOctal%rhou(subcell))
	              if (thisOctal%ghostCell(subcell)) speed(nSquares) = 0.
                 endif
              case DEFAULT
              end select
           endif

       endif
    enddo
  end subroutine recursGetSquares



  subroutine receiveAcrossMpiBoundary(grid, boundaryType, receiveThread, sendThread)

    include 'mpif.h'
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal, tOctal
    type(octal), pointer  :: neighbourOctal
    character(len=*) :: boundaryType
    integer :: receiveThread, sendThread, tsubcell
    integer :: myRank, ierr
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: nOctals
    integer, parameter :: nStorage = 12
    real(double) :: loc(3), tempStorage(nStorage)
    type(VECTOR) :: octVec, direction, rVec
    integer :: nBound
    integer :: iOctal
    integer :: subcell, neighbourSubcell
    integer :: tag = 77
    integer :: status(MPI_STATUS_SIZE)
    logical :: sendLoop
    integer :: nDepth
    real(double) :: q , rho, rhoe, rhou, rhov, rhow, pressure, phi, flux

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
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
    

    if (myRank == receiveThread) then
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
                if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle
                
                octVec = subcellCentre(thisOctal, subcell) + &
                     (thisOctal%subcellSize/2.d0+0.01d0 * grid%halfSmallestSubcell) * direction
                
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
                   write(*,*) "Neighbour on ",boundaryType, " of ", myrankglobal,"  is not on thread ", sendThread, " but ", &
                   neighbourOctal%mpiThread(neighboursubcell)
                   stop
                endif

                loc(1) = octVec%x
                loc(2) = octVec%y
                loc(3) = octVec%z
!                write(*,*) myRank, " has identified a boundary cell ", loc(1:3)

                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_SEND(thisOctal%nDepth, 1, MPI_INTEGER, sendThread, tag, MPI_COMM_WORLD, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, status, ierr)
!                write(*,*) myrank, " received temp storage"
                if (.not.associated(thisOctal%mpiBoundaryStorage)) then
                   allocate(thisOctal%mpiBoundaryStorage(1:thisOctal%maxChildren, 6, nStorage))
                   thisOctal%mpiBoundaryStorage = 0.d0
                endif
                thisOctal%mpiBoundaryStorage(subcell, nBound, 1:nStorage) = tempStorage(1:nStorage)
!                write(*,*) myrank, " successfully stored"

             end if
          enddo
       enddo
       deallocate(octalArray)
       loc(1) = HUGE(loc(1))
       loc(2) = 0.d0
       loc(3) = 0.d0
!       write(*,*) myrank, " sending a huge value to ", sendThread
       call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, ierr)

       ! now send a finish signal to the sendThread
    else
       if (myRank /= sendThread) then
          write(*,*) "subroutine called within thread ", myRank, " but expecing to be ", sendthread
          stop
       endif
       sendLoop = .true.
       do while (sendLoop)
          ! receive a locator
!          write(*,*) myrank, " waiting for a locator"
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, receiveThread, tag, MPI_COMM_WORLD, status, ierr)
!          write(*,*) myrank, " received a locator ", loc(1:3)
          if (loc(1) > 1.d20) then
             sendLoop = .false.
!             write(*,*) myRank, " found the signal to end the send loop"
          else

             octVec = VECTOR(loc(1), loc(2), loc(3))

             call MPI_RECV(nDepth, 1, MPI_INTEGER, receiveThread, tag, MPI_COMM_WORLD, status, ierr)

             call findSubcellTD(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell)

             if (neighbourOctal%mpiThread(neighboursubcell) /= sendthread) then
                write(*,*) "trying to send on ",boundaryType, " but is not on thread ", sendThread
                stop
             endif

             if (neighbourOctal%nDepth <= nDepth) then
!                if ((nDepth - neighbourOctal%nDepth) > 1) then
!                   write(*,*) "Octal depth differs by more than 1 across boundary!!!"
!                   write(*,*) "ndepth ",nDepth, " neighbour%nDepth ",neighbourOctal%nDepth
!                   write(*,*) "myrank ",myrank
!                   write(*,*) "sendThread ",sendThread, " receivethread ",receivethread
!                   stop
!                endif
                tempStorage(1) = neighbourOctal%q_i(neighbourSubcell)
                tempStorage(2) = neighbourOctal%rho(neighbourSubcell)
                tempStorage(3) = neighbourOctal%rhoe(neighbourSubcell)
                tempStorage(4) = neighbourOctal%rhou(neighbourSubcell)
                tempStorage(5) = neighbourOctal%rhov(neighbourSubcell)
                tempStorage(6) = neighbourOctal%rhow(neighbourSubcell)
                tempStorage(7) = neighbourOctal%x_i(neighbourSubcell)
                rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
                     direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
                tOctal => neighbourOctal
                tSubcell = neighbourSubcell
                call findSubcellLocal(rVec, tOctal, tSubcell)
                tempStorage(8) = tOctal%q_i(tsubcell)
                
                tempStorage(9) = dble(neighbourOctal%nDepth)
                tempStorage(10) = neighbourOctal%pressure_i(neighbourSubcell)
                tempStorage(11) = neighbourOctal%flux_i(neighbourSubcell)

                tempStorage(12) = neighbourOctal%phi_i(neighbourSubcell)


!                write(*,*) myrank," set up tempstorage with ", &
!                     tempstorage(1:nStorage),neighbourOctal%nDepth, neighbourSubcell,neighbourOctal%ghostCell(neighbourSubcell), &
!                     neighbourOctal%edgeCell(neighbourSubcell)

             else ! need to average
                call averageValue(direction, neighbourOctal,  neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, flux, phi)
                tempStorage(1) = q
                tempStorage(2) = rho
                tempStorage(3) = rhoe
                tempStorage(4) = rhou
                tempStorage(5) = rhov
                tempStorage(6) = rhow
                tempStorage(7) = neighbourOctal%x_i(neighbourSubcell)
                rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
                     direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
                tOctal => neighbourOctal
                tSubcell = neighbourSubcell
                call findSubcellLocal(rVec, tOctal, tSubcell)
                tempStorage(8) = tOctal%q_i(tsubcell)
                
                tempStorage(9) = dble(neighbourOctal%nDepth)
                tempStorage(10) = pressure
                tempStorage(11) = flux
                tempStorage(12) = phi
             endif
!                          write(*,*) myRank, " sending temp storage ", tempStorage(1:nStorage)
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, receiveThread, tag, MPI_COMM_WORLD, ierr)
!                          write(*,*) myRank, " temp storage sent"

             


          endif
       end do
    endif
  end subroutine receiveAcrossMpiBoundary
  
  subroutine exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup, useThisBound)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    integer :: iPair, nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup, iGroup
    integer, optional :: useThisBound
    integer :: rBound
    integer :: myRank, ierr
    logical :: doExchange
    character(len=10) :: boundaryType(6) = (/"left  ","right ", "top   ", "bottom", "front ", "back  "/)

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    CALL MPI_BARRIER(amrCOMMUNICATOR, ierr)

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
                if ((myRank == thread1(iPair)).or.(myRank == thread2(iPair))) then
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

  end subroutine exchangeAcrossMPIboundary

  subroutine receiveAcrossMpiBoundaryLevel(grid, boundaryType, receiveThread, sendThread, nDepth)

    include 'mpif.h'
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal, tOctal
    type(octal), pointer  :: neighbourOctal
    character(len=*) :: boundaryType
    integer :: receiveThread, sendThread, tsubcell
    integer :: myRank, ierr
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: nOctals
    integer, parameter :: nStorage = 12
    real(double) :: loc(3), tempStorage(nStorage)
    type(VECTOR) :: octVec, direction, rVec
    integer :: nBound
    integer :: iOctal
    integer :: nDepth, thisnDepth
    integer :: subcell, neighbourSubcell
    integer :: tag = 77
    integer :: status(MPI_STATUS_SIZE)
    logical :: sendLoop

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
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
    

    if (myRank == receiveThread) then
       allocate(octalArray(grid%nOctals))
       nOctals = 0
       call getOctalArrayLevel(grid%octreeRoot,octalArray, nOctals, nDepth)
!       write(*,*) myrank," generated ",nOctals, " array of octals"
       do iOctal =  1, nOctals
          
          thisOctal => octalArray(iOctal)%content
!          write(*,*) myrank, " doing octal of ", thisOctal%ndepth, " depth "
          
          do subcell = 1, thisOctal%maxChildren
             
                
!                if (thisOctal%mpiThread(subcell) /= myRank) cycle
                if (.not.octalOnThread(thisOctal, subcell, myRank)) cycle
                
                octVec = subcellCentre(thisOctal, subcell) + &
                     (thisOctal%subcellSize/2.d0+0.01d0 * grid%halfSmallestSubcell) * direction
                
                if (.not.inOctal(grid%octreeRoot, octVec)) then
                   write(*,*) "Grid doesn't have a ", boundaryType, " surface in this volume"
                   write(*,*) "centre",subcellCentre(thisOctal, subcell)
                   write(*,*) "depth",thisOctal%nDepth, thisOctal%haschild(1:8)
                   write(*,*) octVec
                   stop
                endif
                
                neighbourOctal => thisOctal
                call findSubcellLocalLevel(octVec, neighbourOctal, neighbourSubcell, nDepth)

                if (octalOnThread(neighbourOctal, neighbourSubcell, receiveThread)) cycle

                if (.not.octalOnThread(neighbourOctal, neighbourSubcell, sendThread)) then
                   write(*,*) "Neighbour on ",boundaryType, " of ", myrankglobal,"  is not on thread ", sendThread, " but ", &
                   neighbourOctal%mpiThread(neighboursubcell)
                   stop
                endif

                loc(1) = octVec%x
                loc(2) = octVec%y
                loc(3) = octVec%z
!                write(*,*) myRank, " has identified a boundary cell ", loc(1:3)

                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_SEND(thisOctal%nDepth, 1, MPI_INTEGER, sendThread, tag, MPI_COMM_WORLD, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, status, ierr)
!                write(*,*) myrank, " received temp storage"
                if (.not.associated(thisOctal%mpiBoundaryStorage)) then
                   allocate(thisOctal%mpiBoundaryStorage(1:thisOctal%maxChildren, 6, nStorage))
                   thisOctal%mpiBoundaryStorage = 0.d0
                endif
                thisOctal%mpiBoundaryStorage(subcell, nBound, 1:nStorage) = tempStorage(1:nStorage)
!                write(*,*) myrank, " successfully stored"

             enddo
       enddo
       deallocate(octalArray)
       loc(1) = HUGE(loc(1))
       loc(2) = 0.d0
       loc(3) = 0.d0
!       write(*,*) myrank, " sending a huge value to ", sendThread
       call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, ierr)

       ! now send a finish signal to the sendThread
    else
       if (myRank /= sendThread) then
          write(*,*) "subroutine called within thread ", myRank, " but expecing to be ", sendthread
          stop
       endif
       sendLoop = .true.
       do while (sendLoop)
          ! receive a locator
!          write(*,*) myrank, " waiting for a locator"
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, receiveThread, tag, MPI_COMM_WORLD, status, ierr)
!          write(*,*) myrank, " received a locator ", loc(1:3)
          if (loc(1) > 1.d20) then
             sendLoop = .false.
!             write(*,*) myRank, " found the signal to end the send loop"
          else

             octVec = VECTOR(loc(1), loc(2), loc(3))

             call MPI_RECV(thisnDepth, 1, MPI_INTEGER, receiveThread, tag, MPI_COMM_WORLD, status, ierr)

             call findSubcellTDLevel(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell, nDepth)

             if (neighbourOctal%mpiThread(neighboursubcell) /= sendthread) then
                write(*,*) "trying to send on ",boundaryType, " but is not on thread ", sendThread
                stop
             endif

             if (neighbourOctal%nDepth <= thisnDepth) then
                tempStorage(1) = neighbourOctal%q_i(neighbourSubcell)
                tempStorage(2) = neighbourOctal%rho(neighbourSubcell)
                tempStorage(3) = neighbourOctal%rhoe(neighbourSubcell)
                tempStorage(4) = neighbourOctal%rhou(neighbourSubcell)
                tempStorage(5) = neighbourOctal%rhov(neighbourSubcell)
                tempStorage(6) = neighbourOctal%rhow(neighbourSubcell)
                tempStorage(7) = neighbourOctal%x_i(neighbourSubcell)
                rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
                     direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
                tOctal => neighbourOctal
                tSubcell = neighbourSubcell
                call findSubcellLocalLevel(rVec, tOctal, tSubcell, nDepth)
                tempStorage(8) = tOctal%q_i(tsubcell)
                
                tempStorage(9) = dble(neighbourOctal%nDepth)
                tempStorage(10) = neighbourOctal%pressure_i(neighbourSubcell)
                tempStorage(11) = neighbourOctal%flux_i(neighbourSubcell)
 

                tempStorage(12) = neighbourOctal%phi_i(neighbourSubcell)

!                          write(*,*) myRank, " sending temp storage"
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, receiveThread, tag, MPI_COMM_WORLD, ierr)
!                          write(*,*) myRank, " temp storage sent"

             else ! need to average

                write(*,*) "Error on receiveAcrossMpiBoundaryLevel"
                stop
             


             endif
          endif
       enddo
    endif
  end subroutine receiveAcrossMpiBoundaryLevel
  
  subroutine exchangeAcrossMPIboundaryLevel(grid, nPairs, thread1, thread2, nBound, group, nGroup, nDepth)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    integer :: iPair, nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup, iGroup
    integer :: rBound
    integer :: myRank, ierr, nDepth
    character(len=10) :: boundaryType(6) = (/"left  ","right ", "top   ", "bottom", "front ", "back  "/)

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    CALL MPI_BARRIER(amrCOMMUNICATOR, ierr)

    do iGroup = 1, nGroup

       do iPair = 1, nPairs

          if (group(iPair) == iGroup) then
             !       write(*,*) "myrank ", iPair, thread1(iPair), thread2(iPair)
             if ((myRank == thread1(iPair)).or.(myRank == thread2(iPair))) then
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

  recursive subroutine determineBoundaryPairs(thisOctal, grid, nPairs,  thread1, thread2, nBound, iThread)

    include 'mpif.h'
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal
    !
    integer :: subcell, i, iThread
    type(VECTOR) :: dirVec(6), centre, octVec
    integer :: neighbourSubcell, j, nDir
    real(double) :: r
    integer :: myRank, ierr
    integer :: thread1(:), thread2(:), nBound(:), nPairs, iPair
    integer :: i1, i2
    logical :: alreadyInList

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

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
               call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)
                
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

  subroutine returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    integer, intent(out) :: nPairs, thread1(:), thread2(:), nBound(:), nGroup, group(:)
    integer :: nThreads, iThread
    integer :: myRank, ierr, i
    integer, allocatable :: indx(:), itmp(:)
    real, allocatable :: sort(:)
    integer :: list(1000), nList


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)

    nPairs = 0
    do iThread = 1, nThreads-1
       call determineBoundaryPairs(grid%octreeRoot, grid, nPairs,  thread1, thread2, nBound, iThread)
    enddo

    if (nPairs > 1) then
       allocate(indx(1:nPairs), sort(1:nPairs), itmp(1:nPairs))
       do i = 1, nPairs
          sort(i) = real(thread1(i))*100. + real(thread2(i))
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
       rhou, rhov, rhow, x, qnext, pressure, flux, phi, nd)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal, tOctal
    type(VECTOR) :: direction, rVec
    integer :: subcell, neighbourSubcell, tSubcell
    integer :: nBound, nDepth
    integer, intent(out) :: nd

    real(double), intent(out) :: q, rho, rhoe, rhou, rhov, rhow, qnext, x, pressure, flux, phi
    integer :: myRank, ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    nbound = getNboundFromDirection(direction)

    if ((thisOctal%twoD).and.((nBound == 5).or. (nBound == 6))) then
       write(*,*) "Bonndary error for twod: ",nbound
    endif
    if (octalOnThread(neighbourOctal, neighbourSubcell, myRank)) then

       nd = neighbourOctal%nDepth

       x = neighbourOctal%x_i(neighbourSubcell)

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

       else if (thisOctal%nDepth > neighbourOctal%nDepth) then ! fine cells set to coarse cell fluxes (should be interpolated here!!!)
          q   = neighbourOctal%q_i(neighbourSubcell)
          rho = neighbourOctal%rho(neighbourSubcell)
          rhoe = neighbourOctal%rhoe(neighbourSubcell)
          rhou = neighbourOctal%rhou(neighbourSubcell)
          rhov = neighbourOctal%rhov(neighbourSubcell)
          rhow = neighbourOctal%rhow(neighbourSubcell)
          pressure = neighbourOctal%pressure_i(neighbourSubcell)
          phi = neighbourOctal%phi_i(neighbourSubcell)
          if (thisOctal%oneD) then
             flux = neighbourOctal%flux_i(neighbourSubcell)
          else if (thisOctal%twoD) then
             flux = neighbourOctal%flux_i(neighbourSubcell)!/2.d0
          else
             flux = neighbourOctal%flux_i(neighbourSubcell)!/4.d0
          endif
       else
          call averageValue(direction, neighbourOctal,  neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, flux, phi) ! fine to coarse
       endif

       
       rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
            direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
       if (inOctal(grid%octreeRoot, rVec)) then
          tOctal => neighbourOctal
          tSubcell = neighbourSubcell
          call findSubcellLocal(rVec, tOctal, tSubcell)
          
          if (tOctal%mpiThread(tSubcell) == myRank) then
             qnext = tOctal%q_i(tSubcell)
          else
             qNext = neighbourOctal%mpiBoundaryStorage(neighbourSubcell, nBound, 1)
          endif
       else
          qNext = 0.d0
       endif



    else


       if (.not.associated(thisOctal%mpiBoundaryStorage)) then
          write(*,*) "boundary storage not allocated when it should be!", myrank, neighbourOctal%mpiThread(neighboursubcell), &
               thisOctal%mpiThread(subcell)
          write(*,*) "direction",  direction,nBound
          write(*,*) "this centre",subcellCentre(thisOctal, subcell)
          write(*,*) "neig centre",subcellCentre(neighbourOctal, neighboursubcell)
          stop
       endif

       x = thisOctal%mpiBoundaryStorage(subcell, nBound, 7)

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

    endif
  end subroutine getNeighbourValues



  subroutine averageValue(direction, neighbourOctal, neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, flux, phi)
    type(OCTAL), pointer ::  neighbourOctal
    integer :: nSubcell(4), neighbourSubcell
    real(double), intent(out) :: q, rho, rhov, rhou, rhow, pressure, flux, rhoe, phi
    real(double) :: fac
    type(VECTOR) :: direction

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
       flux = neighbourOctal%flux_i(neighbourSubcell)
       phi = neighbourOctal%phi_i(neighbourSubcell)

     else if (neighbourOctal%twoD) then
       if (direction%x > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 3
       else if (direction%x < -0.9d0) then
          nSubcell(1) = 2
          nSubcell(2) = 4
       else if (direction%z > 0.9d0) then
          nSubcell(1) = 1
          nSubcell(2) = 2
       else
          nSubcell(1) = 3
          nSubcell(2) = 4
       endif
       q = 0.5d0*(neighbourOctal%q_i(nSubcell(1)) + neighbourOctal%q_i(nSubcell(2)))
       rho = 0.5d0*(neighbourOctal%rho(nSubcell(1)) + neighbourOctal%rho(nSubcell(2)))
       rhoe = 0.5d0*(neighbourOctal%rhoe(nSubcell(1)) + neighbourOctal%rhoe(nSubcell(2)))
       rhou = 0.5d0*(neighbourOctal%rhou(nSubcell(1)) + neighbourOctal%rhou(nSubcell(2)))
       rhov = 0.5d0*(neighbourOctal%rhov(nSubcell(1)) + neighbourOctal%rhov(nSubcell(2)))
       rhow = 0.5d0*(neighbourOctal%rhow(nSubcell(1)) + neighbourOctal%rhow(nSubcell(2)))
       pressure = 0.5d0*(neighbourOctal%pressure_i(nSubcell(1)) + neighbourOctal%pressure_i(nSubcell(2)))
       flux = 0.5d0*(neighbourOctal%flux_i(nSubcell(1)) + neighbourOctal%flux_i(nSubcell(2)))
       phi = 0.5d0*(neighbourOctal%phi_i(nSubcell(1)) + neighbourOctal%phi_i(nSubcell(2)))
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

       pressure = fac*(neighbourOctal%pressure_i(nSubcell(1)) + neighbourOctal%pressure_i(nSubcell(2)) + & 
            neighbourOctal%pressure_i(nSubcell(3)) + neighbourOctal%pressure_i(nSubcell(4)))

       flux = fac * (neighbourOctal%flux_i(nSubcell(1)) + neighbourOctal%flux_i(nSubcell(2)) + & 
            neighbourOctal%flux_i(nSubcell(3)) + neighbourOctal%flux_i(nSubcell(4)))


       phi = fac*(neighbourOctal%phi_i(nSubcell(1)) + neighbourOctal%phi_i(nSubcell(2)) + & 
            neighbourOctal%phi_i(nSubcell(3)) + neighbourOctal%phi_i(nSubcell(4)))

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

  subroutine columnAlongPathAMR(grid, rVec, direction, sigma)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, direction, currentPosition
    real(double) :: sigma, distToNextCell
    type(OCTAL), pointer :: thisOctal, sOctal
    real(double) :: fudgeFac = 1.d-3
    integer :: subcell
    real(double) ::  totDist

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
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
       if (myrank == thisOctal%mpiThread(subcell)) then
          sigma = sigma + distToNextCell*thisOctal%rho(subcell)
       endif

    end do
  end subroutine columnAlongPathAMR

  subroutine countSubcellsMPI(grid, nSubcells, nSubcellArray)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    integer :: nSubcells
    integer, optional :: nSubcellArray(:)
    integer :: ierr, myRank, iThread, nThreads, nVoxels, n
    integer :: iBuffer(1)
    integer, allocatable :: iArray(:)
    integer :: tag = 81
    integer :: status(MPI_STATUS_SIZE)


    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)
    allocate(iArray(1:nThreads-1))

    nSubcells = 0
    if (myrank == 1) then
       call countVoxelsMPI(grid%octreeRoot,nVoxels)
       nSubcells = nSubcells + nVoxels
       iArray(1) = nSubcells
       do iThread = 2, nThreads - 1
          call MPI_RECV(n, 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
          iArray(iThread) = n
          nSubcells = nSubcells + n
       end do
    else
       call countVoxelsMPI(grid%octreeRoot,nVoxels)
       call MPI_SEND(nVoxels, 1, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, ierr)
    endif
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    iBuffer(1) = nSubcells
    call MPI_BCAST(iBuffer, 1, MPI_INTEGER, 1, MPI_COMM_WORLD, ierr)
    nSubcells = iBuffer(1)

    call MPI_BCAST(iArray, nThreads-1, MPI_INTEGER, 1, MPI_COMM_WORLD, ierr)

    if (present(nSubcellArray)) nSubcellArray = iArray
    deallocate(iArray)

  end subroutine countSubcellsMPI

  SUBROUTINE countVoxelsMPI(thisOctal, nVoxels)  
    ! count the number of octals in the current section of the grid.
    ! also counts the number of unique volume elements (voxels) i.e.
    !   those subcells that are not subdivided

    IMPLICIT NONE

    TYPE(OCTAL), POINTER  :: thisOctal 
    INTEGER,INTENT(INOUT) :: nVoxels   ! number of childless subcells
    
    nVoxels = 0
    CALL countVoxelsPrivate(thisOctal, nVoxels)
    
    CONTAINS
    
      RECURSIVE SUBROUTINE countVoxelsPrivate(thisOctal, nVoxels)
        include 'mpif.h'
        integer :: nVoxels
        integer :: subcell
        TYPE(OCTAL), POINTER  :: thisOctal 
        TYPE(OCTAL), POINTER  :: child
        INTEGER :: i



      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call countVoxelsPrivate(child, nVoxels)
                  exit
               end if
            end do
         else

            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

            nVoxels = nVoxels + 1
         endif
      enddo


      end SUBROUTINE countVoxelsPrivate

  END SUBROUTINE countVoxelsMPI

  function octalOnThread(thisOctal, subcell, myRank) result(check)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: myRank
    logical :: check
    integer :: nHydroThreads
    integer :: nFirstLevel
    nHydroThreads = nThreadsGlobal - 1

    check = .true.

    if (thisOctal%twoD) then
       if (nHydroThreads == 4) then
          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif
       if (nHydroThreads == 16) then
          nFirstLevel = (myRank-1) / 4 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else
             if (thisOctal%mpiThread(subcell) == myRank) then
                check = .true.
             else
                check = .false.
             endif
          endif
       endif
    endif

    if (thisOctal%threeD) then
       if (nHydroThreads == 8) then
          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif
       if (nHydroThreads == 64) then
          nFirstLevel = (myRank-1) / 8 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else
             if (thisOctal%mpiThread(subcell) == myRank) then
                check = .true.
             else
                check = .false.
             endif
          endif
!          write(*,*) "thread ", myrankGlobal, " depth ",thisOctal%ndepth, " mpithread ", thisOctal%mpiThread(subcell), check
       endif
    endif

    if (thisOctal%oned) then
       if (nHydroThreads == 2) then
          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif
       if (nHydroThreads == 4) then
          nFirstLevel = (myRank-1) / 2 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else
             if (thisOctal%mpiThread(subcell) == myRank) then
                check = .true.
             else
                check = .false.
             endif
          endif
       endif
    endif
  end function octalOnThread

  subroutine periodBoundary(grid)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    integer :: iThread
    integer :: ierr, i
    real(double) :: loc(3)
    integer :: tag = 78

    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    do iThread = 1, nThreadsGlobal - 1
       if (iThread /= myRankGlobal) then
!          write(*,*) myRankGlobal, " calling boundaryreceiverequests"
          call periodBoundaryReceiveRequests(grid, iThread)
       else
!          write(*,*) "now doing ", myRankGlobal
          call recursivePeriodSend(grid%octreeRoot)
          loc(1) = 1.d30
          do i = 1, nThreadsGlobal-1
             if (i /= iThread) then
!                write(*,*) myRankGlobal, " sending terminate to ", i
                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, ierr)
             endif
          enddo
       endif
!       write(*,*) myrankGlobal, " waiting at barrier"
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!       write(*,*) myrankGlobal, " dropped through barrier"
    enddo
!    write(*,*) myRankGlobal, " HAS REACHED THE BARRIER"
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
  end subroutine periodBoundary

  recursive subroutine recursivePeriodSend(thisOctal)

    include 'mpif.h'
    type(octal), pointer   :: thisOctal, tOctal, child
    integer :: tSubcell
    real(double) :: loc(3), tempStorage(8)
    integer :: subcell, i
    integer :: tag1 = 78, tag2 = 79
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)


    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call recursivePeriodSend(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle

          if (thisOctal%ghostCell(subcell).and.thisOctal%boundaryCondition(subcell)==2) then
             loc(1) = thisOctal%boundaryPartner(subcell)%x
             loc(2) = thisOctal%boundaryPartner(subcell)%y
             loc(3) = thisOctal%boundaryPartner(subcell)%z


             tSubcell = 1
             tOctal => thisOctal
             call findSubcellLocal(thisOctal%boundaryPartner(subcell), tOctal, tSubcell)
             if (tOctal%mpiThread(tSubcell) == myRankGlobal) then
                write(*,*) "locator problem2 ",myRankGlobal
                stop
             endif

             tOctal => thisOctal
             tSubcell = 1
             call findSubcellLocal(thisOctal%boundaryPartner(subcell), tOctal,tsubcell)
!             write(*,*) "boundary partner ", thisOctal%boundaryPartner(subcell)
!             write(*,*) myrankGlobal, " sending locator to ", tOctal%mpiThread(tsubcell)
             call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag1, MPI_COMM_WORLD, ierr)
!             write(*,*) myRankGlobal, " awaiting recv from ", tOctal%mpiThread(tsubcell)
             call MPI_RECV(tempStorage, 8, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag2, MPI_COMM_WORLD, status, ierr)
!             write(*,*) myrankglobal, " received from ",tOctal%mpiThread(tSubcell)
             if (.not.associated(thisOctal%tempStorage)) then
                allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:8))
                thisOctal%tempStorage = 0.d0
             endif
             thisOctal%tempStorage(subcell,1:8) = tempStorage(1:8)
          endif
       endif
    enddo
  end subroutine recursivePeriodSend

  subroutine periodBoundaryReceiveRequests(grid, receiveThread)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    logical :: sendLoop
    real(double) :: loc(3), tempStorage(7)
    integer :: ierr, receiveThread
    integer :: status(MPI_STATUS_SIZE)
    integer :: subcell
    integer :: tag1 = 78, tag2 = 79
    type(VECTOR) :: octVec
    sendLoop = .true.
!    write(*,*) myrankGlobal, " waiting for a locator"
    do while (sendLoop)
       ! receive a locator
       
       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, receiveThread, tag1, MPI_COMM_WORLD, status, ierr)
!       write(*,*) myrankglobal, " received a locator from ", receiveThread 
       if (loc(1) > 1.d20) then
          sendLoop = .false.
!          write(*,*) myRankGlobal, " found the signal to end the send loop from ", receivethread
       else
          octVec = VECTOR(loc(1), loc(2), loc(3))
          thisOctal => grid%octreeRoot
          subcell = 1
!          write(*,*) myrankglobal," calling subscell local"
          call findSubcellLocal(octVec, thisOctal, subcell)
!          write(*,*) myrankglobal," subscell local done succesfully"
          
          tempstorage(1) = thisOctal%rho(Subcell)
          tempStorage(2) = thisOctal%rhoE(Subcell)
          tempStorage(3) = thisOctal%rhou(Subcell)
          tempStorage(4) = thisOctal%rhov(Subcell)
          tempStorage(5) = thisOctal%rhow(Subcell)
          tempStorage(6) = thisOctal%energy(Subcell)
          tempStorage(7) = thisOctal%pressure_i(Subcell)
!          tempStorage(8) = thisOctal%phi_i(Subcell)
!          write(*,*) myRankGlobal, " sending tempstorage to ", receiveThread
          call MPI_SEND(tempStorage, 8, MPI_DOUBLE_PRECISION, receiveThread, tag2, MPI_COMM_WORLD, ierr)
       endif
    enddo
!    write(*,*) myrankGlobal, " leaving receive requests ", sendLoop
  end subroutine periodBoundaryReceiveRequests


  subroutine dumpValuesAlongLine(grid, thisFile, startPoint, endPoint, nPoints)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    integer :: nPoints
    type(VECTOR) :: startPoint, endPoint, position, direction, cen
    real(double) :: loc(3), rho
    character(len=*) :: thisFile
    integer :: ierr
    integer, parameter :: nStorage = 5
    real(double) :: tempSTorage(nStorage), tval
    integer, parameter :: tag = 50
    integer :: status(MPI_STATUS_SIZE)
    logical :: stillLooping
    integer :: sendThread

    thisOctal => grid%octreeRoot
    position = startPoint
    direction = endPoint - startPoint
    call normalize(direction)


    if (myrankGlobal == 0) then

       open(20, file=thisFile, form="formatted", status="unknown")
       do while(inOctal(grid%octreeRoot, position))
          call findSubcellLocal(position, thisOctal, subcell)
          sendThread = thisOctal%mpiThread(subcell)
          loc(1) = position%x
          loc(2) = position%y
          loc(3) = position%z
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, ierr)
          call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, status, ierr)

          cen%x = tempStorage(1)
          cen%y = tempStorage(2)
          cen%z = tempStorage(3)
          rho = tempStorage(4)
          tval = tempStorage(5)
          write(20,*) cen%x, rho
          position = cen
          position = position + (tVal+1.d-3*grid%halfSmallestSubcell)*direction
       enddo
       do sendThread = 1, nThreadsGlobal-1
          loc(1) = 1.d30
          loc(2) = 1.d30
          loc(3) = 1.d30
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, ierr)
       enddo
       close(20)
       goto 666


    else
       stillLooping = .true.
       do while(stillLooping)
          call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, status, ierr)
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
             tempStorage(4) = thisOctal%rho(subcell)             
             tempStorage(5) = tVal
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, ierr)
          endif
       enddo
    endif

666 continue
  end subroutine dumpValuesAlongLine
    
  subroutine grid_info_mpi(thisGrid, filename)
    use amr_mod, only:  countVoxels
    include 'mpif.h'
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
    
    if (myrankGlobal == 0) then
       if (filename(1:1) == '*') then
          UN = 6   ! prints on screen
       else
          UN = 69
          open(unit=UN, file = TRIM(filename), status = 'replace',form='formatted')
       end if
    endif



    nOctals=0; nVoxels=0
    call countVoxels(thisGrid%octreeRoot,nOctals,nVoxels)
    if (myrankGlobal == 0) then
       nOctals = 0
       nVoxels = 0
    endif
    write(*,*) myrankGlobal, " has ", noctals, " octals and ",nvoxels, " voxels"
    call MPI_REDUCE(nOctals, tempInt,1,MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    nOctals = tempInt(1)
    call MPI_REDUCE(nVoxels, tempInt,1,MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    nVoxels = tempInt(1)
    call MPI_REDUCE(thisgrid%maxDepth, tempInt,1,MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
    maxDepth = tempInt(1)
    call MPI_REDUCE(thisgrid%halfSmallestSubcell, tempDouble,1,MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    halfSmallestSubcell = tempDouble(1)
    
    
    if (myRankGlobal == 0) then
       write(UN,'(a)') ' '
       write(UN,'(a)') '######################################################'
       write(UN,'(a)') 'Grid info :'
       write(UN,'(a)') ' '
       write(UN,*)     'geometry             = ', thisGrid%geometry
       write(UN,*)     'maxDepth             = ', maxDepth
       write(UN,*)     'halfSmallestSubcell  = ', halfSmallestSubcell, ' [10^10 cm]'
       write(UN,*)     'nOctals              = ', nOctals
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
    if (myrankGlobal == 0) then
       if (filename(1:1) /= '*')  close(UN)
    endif
    
  end subroutine grid_info_mpi
  

    

#else

  use messages_mod
  use octal_mod

  implicit none

  contains

! Stubbed version to allow non-MPI code to compile
  function octalOnThread(thisOctal, subcell, myRank) result(check)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: myRank
    logical :: check

! Set return value of function to prevent a compiler warning. 
    check=.false.
    call writefatal("octalOnThread called in non-MPI code")
    STOP

  end function octalOnThread



#endif
end module mpi_amr_mod
