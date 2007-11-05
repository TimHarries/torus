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


  end subroutine setupAMRCOMMUNICATOR

  

  subroutine plotGridMPI(grid, device, plane, valueName, valueMinFlag, valueMaxFlag, logFlag, plotgrid)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    character(len=*) :: device, plane, valueName
    character(len=3) planes(3)
    integer :: nPlanes, iPlane
    real, allocatable :: corners(:,:), value(:)
    integer :: nSquares
    integer :: pgbegin, i, j, idx
    integer :: iLo, iHi
    real :: valueMin, valueMax
    real,optional :: valueMinFlag, valueMaxFlag
    real :: xStart, xEnd, yStart, yEnd, t
    integer, parameter :: maxSquares = 100000
    integer :: myRank, ierr

    logical, optional :: logFlag
    logical, optional :: plotgrid
    logical :: logScale, doplotGrid

    logScale = .false.
    if (present(logflag)) logscale = logFlag

    doplotgrid = .false.
    if (present(plotgrid)) doplotgrid = plotgrid
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    allocate(corners(maxSquares, 4))
    allocate(value(maxSquares))

    if (grid%octreeRoot%twoD) then
       nplanes = 1
       planes(1) = "x-z"
    else if (grid%octreeRoot%threeD) then
       nPlanes = 3
       planes(1) = "x-z"
       planes(2) = "x-y"
       planes(3) = "y-z"
    endif

    if (myRank == 1) then
       if (grid%octreeRoot%twoD) then
          i = pgbegin(0, device, 1, 1)
       else if (grid%octreeRoot%threeD) then
          i = pgbegin(0, device, 3, 1)
       endif
    endif
    do iplane = 1, nPlanes
       call getSquares(grid, planes(iplane), valueName, nSquares, corners, value)
       
       if (myRank == 1) then
          xStart = grid%octreeRoot%centre%x-grid%octreeRoot%subcellSize
          yStart = grid%octreeRoot%centre%y-grid%octreeRoot%subcellSize
          xEnd = grid%octreeRoot%centre%x+grid%octreeRoot%subcellSize
          yEnd = grid%octreeRoot%centre%y+grid%octreeRoot%subcellSize
          if (iPlane == 1) then
             valueMin = MINVAL(value(1:nSquares))
             valueMax = MAXVAL(value(1:nSquares))
          endif
          call pgpage
          if (present(valueMinFlag)) valueMin = valueMinFlag
          if (present(valueMaxFlag)) valueMax = valueMaxFlag
          call pgvport(0.1, 0.9, 0.1, 0.9)
          call pgwnad(xStart, xEnd, yStart, yEnd)
          
          call pgqcir(ilo, ihi)
          call pgscir(ilo, ihi)
          call palette(2)
          
          
          do i = 1, nSquares
             if (.not.logscale) then
                if (valueMax /= valueMin) then
                   t = (value(i)-valueMin)/(valueMax-valueMin)
                else
                   t = 0.
                endif
             else
                t = (log10(value(i))-log10(valuemin))/(log10(valuemax)-log10(valuemin))
             endif
             if (t < 0.) t = 0.
             if (t > 1.) t = 1.
             idx = int(t * real(ihi - ilo) + real(ilo))
             
             call pgsci(idx)
             
             call pgrect(corners(i, 1), corners(i, 2), corners(i, 3), corners(i, 4))
             
             if (doplotgrid) then
                call pgqci(j)
                call PGSFS(2)  ! we don't want to fill in a box this time
                call PGSCI(1) ! changing the color index.
                call pgrect(corners(i, 1), corners(i, 2), corners(i, 3), corners(i, 4))
                call PGSCI(1) ! change color index to default.
                call PGSFS(1)
                call pgsci(j)
             endif
             
          enddo
          call pgsci(1)
          call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
          if(logscale) then
             CALL PGWEDG('BI', 4.0, 5.0, real(log10(valueMin)), real(log10(valueMax)), TRIM(ADJUSTL(valueName)) )
          else
             CALL PGWEDG('BI', 4.0, 5.0, real(valueMin), real(valueMax), TRIM(ADJUSTL(valueName)))
          end if
       endif

    enddo
    call pgend
    deallocate(corners, value)
  end subroutine plotGridMPI

  subroutine getSquares(grid, plane, valueName, nSquares, corners, value)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    character(len=*) :: plane, valueName
    real :: corners(:,:), value(:)
    integer :: myRank, ierr, nThreads, iThread
    integer :: nSquares, n, tag=97,i
    integer :: status(MPI_STATUS_SIZE)

    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)

    if (myRank == 1) then
       nSquares = 0
       call recursGetSquares(grid%octreeRoot, grid, plane, valueName, nSquares, corners, value)
       do iThread = 2, nThreads - 1
          call MPI_RECV(n, 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
          do i = 1, n
             nSquares = nSquares + 1
             call MPI_RECV(corners(nSquares, 1:4), 4, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
             call MPI_RECV(value(nSquares), 1, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
          end do
       end do
    else
       nSquares = 0 
       call recursGetSquares(grid%octreeRoot, grid, plane, valueName, nSquares, corners, value)
       call MPI_SEND(nSquares, 1, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, ierr)
       do i = 1, nSquares      
          call MPI_SEND(corners(i,1:4), 4, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
          call MPI_SEND(value(i), 1, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
       enddo
    endif
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
  end subroutine getSquares


  recursive subroutine recursGetSquares(thisOctal, grid, plane, valueName, nSquares, corners, value)

    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    !
    integer :: subcell, i
    character(len=*) :: valueName, plane
    integer :: nSquares
    real :: corners(:,:)
    real :: value(:), tmp
    real(double) :: eps, x
    integer :: myRank, ierr
    type(OCTALVECTOR) :: rVec

    eps = 0.01d0 * grid%halfSmallestSubcell

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    do subcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call recursGetSquares(child, grid, plane, valueName, nSquares, corners, value)
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
             case("chi")
                tmp = real(thisOctal%chiline(subcell))
             case DEFAULT
           end select
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
                endif
              case DEFAULT
              end select

       endif
    enddo
  end subroutine recursGetSquares



  subroutine receiveAcrossMpiBoundary(grid, boundaryType, receiveThread, sendThread)

    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal, tOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    character(len=*) :: boundaryType
    integer :: receiveThread, sendThread, tsubcell
    integer :: myRank, ierr
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: nOctals
    real(double) :: loc(3), tempStorage(15)
    type(OCTALVECTOR) :: octVec, direction, centre, rVec
    integer :: nBound
    integer :: iOctal
    integer :: subcell, neighbourSubcell
    integer :: tag = 77
    integer :: status(MPI_STATUS_SIZE)
    logical :: sendLoop
    integer :: nDepth
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, flux, pressure

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    select case(boundaryType)
    case("top")
       direction = OCTALVECTOR(0.d0, 0.d0, 1.d0)
       nBound = 1
    case("bottom")
       direction = OCTALVECTOR(0.d0, 0.d0, -1.d0)
       nBound = 2
    case("left")
       direction = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
       nBound = 3
    case("right")
       direction = OCTALVECTOR(1.d0, 0.d0, 0.d0)
       nBound = 4
    case("front")
       direction = OCTALVECTOR(0.d0, 1.d0, 0.d0)
       nBound = 5
    case("back")
       direction = OCTALVECTOR(0.d0, -1.d0, 0.d0)
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

                call MPI_RECV(tempStorage, 15, MPI_DOUBLE_PRECISION, sendThread, tag, MPI_COMM_WORLD, status, ierr)
!                write(*,*) myrank, " received temp storage"
                
                if (.not.associated(thisOctal%mpiBoundaryStorage)) then
                   allocate(thisOctal%mpiBoundaryStorage(1:thisOctal%maxChildren, 6, 15))
                   thisOctal%mpiBoundaryStorage = 0.d0
                endif
                thisOctal%mpiBoundaryStorage(subcell, nBound, 1:15) = tempStorage(1:15)
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

             octVec = OCTALVECTOR(loc(1), loc(2), loc(3))

             call MPI_RECV(nDepth, 1, MPI_INTEGER, receiveThread, tag, MPI_COMM_WORLD, status, ierr)

             call findSubcellTD(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell)

             if (neighbourOctal%mpiThread(neighboursubcell) /= sendthread) then
                write(*,*) "trying to send on ",boundaryType, " but is not on thread ", sendThread
                stop
             endif

             if (neighbourOctal%nDepth <= nDepth) then
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
             else ! need to average
                call averageValue(direction, neighbourOctal,  neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, flux)
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
             endif
!                          write(*,*) myRank, " sending temp storage"
             call MPI_SEND(tempStorage, 15, MPI_DOUBLE_PRECISION, receiveThread, tag, MPI_COMM_WORLD, ierr)
!                          write(*,*) myRank, " temp storage sent"

          endif
       end do
    endif
  end subroutine receiveAcrossMpiBoundary
  
  subroutine exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    integer :: iPair, nPairs, thread1(:), thread2(:), nBound(:)
    integer :: group(:), nGroup, iGroup
    integer :: rBound
    integer :: myRank, ierr
    character(len=10) :: boundaryType(6) = (/"top   ","bottom","left  ","right ", "front ", "back  "/)

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    CALL MPI_BARRIER(amrCOMMUNICATOR, ierr)

    do iGroup = 1, nGroup

       do iPair = 1, nPairs

          if (group(iPair) == iGroup) then
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
       enddo
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    enddo

  end subroutine exchangeAcrossMPIboundary

  recursive subroutine determineBoundaryPairs(thisOctal, grid, nPairs,  thread1, thread2, nBound, iThread)

    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    !
    integer :: subcell, i, iThread
    logical :: converged, converged_tmp
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat
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
             dirVec(1) = OCTALVECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
             dirVec(3) = OCTALVECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(4) = OCTALVECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(5) = OCTALVECTOR( 0.d0, 1.d0,  0.d0)
             dirVec(6) = OCTALVECTOR( 0.d0,-1.d0,  0.d0)
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = OCTALVECTOR(0.d0, 0.d0, +1.d0)
             dirVec(2) = OCTALVECTOR(0.d0, 0.d0, -1.d0)
             dirVec(3) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
             dirVec(4) = OCTALVECTOR(+1.d0, 0.d0, 0.d0)
          else
             nDir = 2
             dirVec(1) = OCTALVECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = OCTALVECTOR(-1.d0, 0.d0, 0.d0)
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
    integer :: nPairs, thread1(:), thread2(:), nBound(:)
    integer :: nThreads, iThread
    integer :: myRank, ierr, i
    integer, allocatable :: indx(:), itmp(:)
    real, allocatable :: sort(:)
    integer :: list(1000), nList, nGroup, group(:), iPair, iStart
    logical :: groupFound


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)

    nPairs = 0
    do iThread = 1, nThreads-1
       call determineBoundaryPairs(grid%octreeRoot, grid, nPairs,  thread1, thread2, nBound, iThread)
    enddo
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
    

    nGroup = 1
    group = 0
    group(1) = nGroup
    iStart = 1
    do while(any(group(1:nPairs)==0))
       nList = 2
       list(1) = thread1(iStart)
       list(2) = thread2(iStart)
       groupFound = .false.
       do i = iStart+1, nPairs
          if (group(i) == 0) then
             if (.not.inList(thread1(i), list, nList).and.(.not.inList(thread2(i), list, nList))) then
                group(i) = nGroup
                list(nList+1) = thread1(i)
                list(nList+2) = thread2(i)
                nList = nList + 2
                groupFound = .true.
             endif
          endif
       enddo
       if (.not.groupFound) then
          group(iStart) = nGroup
       endif
       nGroup = nGroup + 1
       do i = 1, nPairs
          if (group(i) == 0) then
             iStart = i
             exit
          endif
       enddo
    enddo
!    if (myRankGlobal == 1) then
!       do iPair = 1, nPairs
!          write(*,*) "Pair: ",iPair, thread1(iPair), " -> ", thread2(iPair), group(iPair)
!       enddo
!    endif
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
       rhou, rhov, rhow, x, qnext, pressure, flux)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal, tOctal
    type(OCTALVECTOR) :: direction, rVec
    integer :: subcell, neighbourSubcell, tSubcell
    integer :: nSubcell(6)
    integer :: nBound, nDepth
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, qnext, x, pressure, flux
    integer :: myRank, ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    nbound = getNboundFromDirection(direction)

    if ((thisOctal%twoD).and.((nBound == 5).or. (nBound == 6))) then
       write(*,*) "Bonndary error for twod: ",nbound
    endif
    if (neighbourOctal%mpiThread(neighbourSubcell) == myRank) then

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

       else if (thisOctal%nDepth > neighbourOctal%nDepth) then ! fine cells set to coarse cell fluxes (should be interpolated here!!!)
          q   = neighbourOctal%q_i(neighbourSubcell)
          rho = neighbourOctal%rho(neighbourSubcell)
          rhoe = neighbourOctal%rhoe(neighbourSubcell)
          rhou = neighbourOctal%rhou(neighbourSubcell)
          rhov = neighbourOctal%rhov(neighbourSubcell)
          rhow = neighbourOctal%rhow(neighbourSubcell)
          pressure = neighbourOctal%pressure_i(neighbourSubcell)
          if (thisOctal%oneD) then
             flux = neighbourOctal%flux_i(neighbourSubcell)
          else if (thisOctal%twoD) then
             flux = neighbourOctal%flux_i(neighbourSubcell)!/2.d0
          else
             flux = neighbourOctal%flux_i(neighbourSubcell)!/4.d0
          endif
       else
          call averageValue(direction, neighbourOctal,  neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, flux) ! fine to coarse
       endif

       
       rVec = subcellCentre(neighbourOctal, neighbourSubcell) + &
            direction * (neighbourOctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
       tOctal => neighbourOctal
       tSubcell = neighbourSubcell
       call findSubcellLocal(rVec, tOctal, tSubcell)
       
       if (tOctal%mpiThread(tSubcell) == myRank) then
          qnext = tOctal%q_i(tSubcell)
       else
          qNext = neighbourOctal%mpiBoundaryStorage(neighbourSubcell, nBound, 1)
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


    endif
  end subroutine getNeighbourValues



  subroutine averageValue(direction, neighbourOctal, neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, flux)
    type(OCTAL), pointer ::  neighbourOctal
    integer :: nSubcell(4), neighbourSubcell
    real(double) :: q, rho, rhov, rhou, rhow, pressure, flux, rhoe
    type(OCTALVECTOR) :: direction

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
          nSubcell(1) = 1
          nSubcell(2) = 2
          nSubcell(3) = 5
          nSubcell(4) = 6
       endif
       q = 0.25d0*(neighbourOctal%q_i(nSubcell(1)) + neighbourOctal%q_i(nSubcell(2)) + & 
            neighbourOctal%q_i(nSubcell(3)) + neighbourOctal%q_i(nSubcell(4)))

       rho = 0.25d0*(neighbourOctal%rho(nSubcell(1)) + neighbourOctal%rho(nSubcell(2)) + & 
            neighbourOctal%rho(nSubcell(3)) + neighbourOctal%rho(nSubcell(4)))

       rhoe = 0.25d0*(neighbourOctal%rhoe(nSubcell(1)) + neighbourOctal%rhoe(nSubcell(2)) + & 
            neighbourOctal%rhoe(nSubcell(3)) + neighbourOctal%rhoe(nSubcell(4)))

       rhou = 0.25d0*(neighbourOctal%rhou(nSubcell(1)) + neighbourOctal%rhou(nSubcell(2)) + & 
            neighbourOctal%rhou(nSubcell(3)) + neighbourOctal%rhou(nSubcell(4)))

       rhov = 0.25d0*(neighbourOctal%rhov(nSubcell(1)) + neighbourOctal%rhov(nSubcell(2)) + & 
            neighbourOctal%rhov(nSubcell(3)) + neighbourOctal%rhov(nSubcell(4)))

       rhow = 0.25d0*(neighbourOctal%rhow(nSubcell(1)) + neighbourOctal%rhow(nSubcell(2)) + & 
            neighbourOctal%rhow(nSubcell(3)) + neighbourOctal%rhow(nSubcell(4)))

       pressure = 0.25d0*(neighbourOctal%pressure_i(nSubcell(1)) + neighbourOctal%pressure_i(nSubcell(2)) + & 
            neighbourOctal%pressure_i(nSubcell(3)) + neighbourOctal%pressure_i(nSubcell(4)))

       flux = 0.25d0*(neighbourOctal%flux_i(nSubcell(1)) + neighbourOctal%flux_i(nSubcell(2)) + & 
            neighbourOctal%flux_i(nSubcell(3)) + neighbourOctal%flux_i(nSubcell(4)))
    endif

  end subroutine averageValue

  function getNBoundFromDirection(direction) result (nBound)
    type(OCTALVECTOR) :: direction
    integer :: nBound
    if       (direction%z > 0.9d0) then
       nBound = 1
    else if (direction%z < -0.9d0) then
       nBound = 2
    else if (direction%x < -0.9d0) then
       nBound = 3
    else if (direction%x > 0.9d0) then
       nBound = 4
    else if (direction%y > 0.9d0) then
       nBound = 5
    else if (direction%y < -0.9d0) then
       nBound = 6
    endif
  end function getNBoundFromDirection

  subroutine columnDensityPlotAMR(grid, viewVec, device, resetRangeFlag, iminfix, imaxfix)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(OCTALVECTOR) :: viewVec
    type(OCTALVECTOR) :: xProj, yProj, startVec
    character(len=*) :: device
    real, allocatable :: image(:,:)
    real(double) :: t
    integer :: nx, ny
    real, optional :: iMinFix, iMaxFix
    real(double) :: sigma(64), tempSigma(64)
    real :: tr(6)
    integer :: i, j
    real :: dx
    real, allocatable :: xAxis(:), yAxis(:)
    logical, optional :: resetRangeFlag
    logical :: resetRange
    real :: imageSize
    real, save :: iMax, iMin
    integer::pgbegin
    integer :: nHydroThreads

    nHydroThreads = nThreadsGlobal - 1


    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    resetRange = .true.
    if (PRESENT(resetRangeFlag)) then
       resetRange = resetRangeFlag
    endif

    nx = 100
    ny = 100
    allocate(image(1:nx,1:ny))
    allocate(xAxis(1:nx))
    allocate(yAxis(1:ny))

    imageSize = sqrt(3.d0)*2.d0*grid%octreeRoot%subcellSize
    dx = imageSize / real(nx)
    do i = 1, nx
       xAxis(i) = -imageSize/2. + dx/2. + real(i-1)*dx
    enddo
    do i = 1, ny
       yAxis(i) = -imageSize/2. + dx/2. + real(i-1)*dx
    enddo
       
    xProj = viewVec .cross. OCTALVECTOR(0.d0, 0.d0, 1.d0)  
    call normalize(xProj)
    yProj = xProj .cross. viewVec
    call normalize(yProj)


    do i = 1, nx
       do j = 1, ny
          startVec = dble(xAxis(i))*xProj + dble(yAxis(j))*yProj
          startVec = startVec - grid%octreeRoot%subcellSize*4.d0 * viewVec
          
          t = distanceToGridFromOutside(grid, startVec, viewVec)
          startVec = startVec + (t + 1.d-1*grid%halfSmallestSubcell) *viewVec

          sigma = 0.d0
          if (.not.inOctal(grid%octreeRoot, startVec)) then
             sigma(myrank) = 1.e-10
          else
             call columnAlongPathAMR(grid, startVec, viewVec, sigma(myRank))
          endif
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       call MPI_ALLREDUCE(sigma, tempSigma, nHydroThreads, MPI_DOUBLE_PRECISION, MPI_SUM, amrCOMMUNICATOR, ierr)
       image(i,j) = max(1.e-10,real(SUM(tempsigma(1:nHydroThreads))))
       enddo
    enddo
    
    tr(1) = xAxis(1)-dx
    tr(2) = dx
    tr(3) = 0.
    tr(4) = yAxis(1)-dx
    tr(5) = 0.
    tr(6) = dx

    if (myRank == 1) then
       i= pgbegin(0,device,1,1)

       call pgvport(0.1, 0.9, 0.1, 0.9)
       call pgwnad(xAxis(1)-dx/2., xAxis(nx)-dx/2.,yAxis(1)-dx/2.,yAxis(nx)-dx/2.)
       
       call palette(2)
 
       if (resetRange) then
          iMin = minval(image)
          iMax = maxVal(image)
       endif

       if (present(iMinFix)) iMin = iMinFix
       if (present(iMinFix)) iMax = iMaxFix

       write(*,*) "Column image: ",imin,imax
       call pgimag(image, nx, ny, 1, nx, 1, ny, imin, imax, tr)
       call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
       call pgend
    endif
  end subroutine columnDensityPlotAMR

  subroutine columnAlongPathAMR(grid, rVec, direction, sigma)
    include 'mpif.h'
    integer :: myRank, ierr
    type(GRIDTYPE) :: grid
    type(OCTALVECTOR) :: rVec, direction, currentPosition
    integer :: iLambda
    real(double) :: sigma, distToNextCell
    type(OCTAL), pointer :: thisOctal, sOctal
    real(double) :: fudgeFac = 1.d-3
    integer :: subcell
    real(double) ::  totDist
    logical :: hitSource

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


  function octalOnThread(thisOctal, subcell, myRank) result(check)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer :: myRank
    logical :: check
    integer :: nHydroThreads
    integer :: nFirstLevel
    nHydroThreads = nThreadsGlobal - 1

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
          call periodBoundaryReceiveRequests(grid, iThread)
       else
          write(*,*) "now doing ", myRankGlobal
          call recursivePeriodSend(grid%octreeRoot)
          loc(1) = 1.d30
          do i = 1, nThreadsGlobal-1
             if (i /= iThread) then
                write(*,*) myRankGlobal, " sending terminate to ", i
                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, ierr)
             endif
          enddo
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       endif
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    enddo
    write(*,*) myRankGlobal, " HAS REACHED THE BARRIER"
    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
  end subroutine periodBoundary

  recursive subroutine recursivePeriodSend(thisOctal)

    include 'mpif.h'
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal, tOctal, child
    integer :: tSubcell
    type(octal), pointer  :: neighbourOctal, startOctal
    real(double) :: loc(3), tempStorage(5)
    !
    integer :: subcell, i
    integer :: myrank
    integer :: tag1 = 78, tag2 = 79
    integer :: ierr
    integer :: status


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

          if (thisOctal%ghostCell(subcell)) then
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
             write(*,*) myrankGlobal, " sending locator to ", tOctal%mpiThread(tsubcell)
             call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag1, MPI_COMM_WORLD, ierr)
             write(*,*) myRankGlobal, " awaiting recv from ", tOctal%mpiThread(tsubcell)
             call MPI_RECV(tempStorage, 5, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag2, MPI_COMM_WORLD, status, ierr)
             write(*,*) myrankglobal, " received from ",tOctal%mpiThread(tSubcell)
             if (.not.associated(thisOctal%tempStorage)) allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:5))
             thisOctal%tempStorage(subcell,1:5) = tempStorage(1:5)
          endif
       endif
    enddo
  end subroutine recursivePeriodSend

  subroutine periodBoundaryReceiveRequests(grid, receiveThread)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, tOCtal
    logical :: sendLoop
    real(double) :: loc(3), tempStorage(5)
    integer :: status, ierr, receiveThread
    integer :: subcell
    integer :: tag1 = 78, tag2 = 79
    type(OCTALVECTOR) :: octVec
    sendLoop = .true.
    do while (sendLoop)
       ! receive a locator
       
       write(*,*) myrankGlobal, " waiting for a locator"
       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, receiveThread, tag1, MPI_COMM_WORLD, status, ierr)
       write(*,*) myrankglobal, " received a locator from ", receiveThread 
       if (loc(1) > 1.d20) then
          sendLoop = .false.
          write(*,*) myRankGlobal, " found the signal to end the send loop from ", receivethread
       else
          octVec = OCTALVECTOR(loc(1), loc(2), loc(3))
          thisOctal => grid%octreeRoot
          subcell = 1
          write(*,*) myrankglobal," calling subscell local"
          call findSubcellLocal(octVec, thisOctal, subcell)
          write(*,*) myrankglobal," subscell local done succesfully"
          
          tempstorage(1) = thisOctal%rho(Subcell)
          tempStorage(2) = thisOctal%rhoE(Subcell)
          tempStorage(3) = thisOctal%rhou(Subcell)
          tempStorage(4) = thisOctal%rhov(Subcell)
          tempStorage(5) = thisOctal%rhow(Subcell)
          write(*,*) myRankGlobal, " sending tempstorage to ", receiveThread
          call MPI_SEND(tempStorage, 5, MPI_DOUBLE_PRECISION, receiveThread, tag2, MPI_COMM_WORLD, ierr)
       endif
    enddo
  end subroutine periodBoundaryReceiveRequests


#endif
end module mpi_amr_mod
