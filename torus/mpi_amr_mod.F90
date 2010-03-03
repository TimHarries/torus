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


    if (nThreadsGlobal == 2) nHydroThreadsGlobal = 2
    if (nThreadsGlobal == 3) nHydroThreadsGlobal = 2

    if (nThreadsGlobal == 4) nHydroThreadsGlobal = 4
    if (nThreadsGlobal == 5) nHydroThreadsGlobal = 4


    if (nThreadsGlobal == 8) nHydroThreadsGlobal = 8
    if (nThreadsGlobal == 9) nHydroThreadsGlobal = 8

    if (nThreadsGlobal == 64) nHydroThreadsGlobal = 64
    if (nThreadsGlobal == 65) nHydroThreadsGlobal = 64

  end subroutine setupAMRCOMMUNICATOR

! Free communicator created by setupAMRCOMMUNICATOR
  subroutine freeAMRCOMMUNICATOR
    implicit none

    include 'mpif.h'
    integer :: ierr

    if ( myRankGlobal /= 0 ) then 
       call MPI_COMM_FREE(amrCOMMUNICATOR, ierr)
    end if

  end subroutine freeAMRCOMMUNICATOR

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
    use input_variables, only : hydrodynamics
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
             if (hydrodynamics) then
                if (thisOctal%twoD) then
                   dv = thisOctal%subcellSize**2/1.d30
                else if (thisOctal%oned) then
                   dv = thisOctal%subcellSize/1.d30
                endif
             endif
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

  subroutine findEnergyOverAllThreads(grid, energy)
    include 'mpif.h'
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

          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
             if (thisOctal%threed) then
                dv = cellVolume(thisOctal, subcell) * 1.d30
             else if (thisOctal%twoD) then
                dv = thisOctal%subcellSize**2
             else if (thisOctal%oneD) then
                dv = thisOctal%subcellSize
             endif
             totalEnergy = totalEnergy + thisOctal%rhoe(subcell) * dv
          endif
       endif
    enddo
  end subroutine findTotalEnergyMPI


  

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
    if (myrankGlobal == 0) goto 666

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
666 continue
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

  subroutine periodBoundary(grid, justGrav)
    include 'mpif.h'
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
    do iThread = 1, nThreadsGlobal - 1
       if (iThread /= myRankGlobal) then
!          write(*,*) myRankGlobal, " calling boundaryreceiverequests"
          call periodBoundaryReceiveRequests(grid, iThread, doJustGrav)
       else
!          write(*,*) "now doing ", myRankGlobal
          call recursivePeriodSend(grid%octreeRoot, doJustGrav)
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

  subroutine periodBoundaryLevel(grid, nDepth, justGrav)
    include 'mpif.h'
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
    do iThread = 1, nThreadsGlobal - 1
       if (iThread /= myRankGlobal) then
!          write(*,*) myRankGlobal, " calling boundaryreceiverequests"
          call periodBoundaryReceiveRequestsLevel(grid, iThread, nDepth, doJustGrav)
       else
!          write(*,*) "now doing ", myRankGlobal
          call recursivePeriodSendLevel(grid%octreeRoot, nDepth, doJustGrav)
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
  end subroutine periodBoundaryLevel

  recursive subroutine recursivePeriodSend(thisOctal, doJustGrav)

    include 'mpif.h'
    type(octal), pointer   :: thisOctal, tOctal, child
    integer :: tSubcell
    real(double) :: loc(3), tempStorage(8)
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
!             write(*,*) "boundary partner ", thisOctal%boundaryPartner(subcell)
!             write(*,*) myrankGlobal, " sending locator to ", tOctal%mpiThread(tsubcell)
             call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag1, MPI_COMM_WORLD, ierr)
!             write(*,*) myRankGlobal, " awaiting recv from ", tOctal%mpiThread(tsubcell)
             call MPI_RECV(tempStorage, 8, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag2, MPI_COMM_WORLD, status, ierr)
!             write(*,*) myrankglobal, " received from ",tOctal%mpiThread(tSubcell)
             if (.not.associated(thisOctal%tempStorage)) then
                if (.not.doJustGrav) then
                   allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:8))
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

    include 'mpif.h'
    type(octal), pointer   :: thisOctal, tOctal, child
    integer :: tSubcell
    integer :: nDepth
    real(double) :: loc(3), tempStorage(8)
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
             call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag1, MPI_COMM_WORLD, ierr)
!             write(*,*) myRankGlobal, " awaiting recv from ", tOctal%mpiThread(tsubcell)
             call MPI_RECV(tempStorage, 8, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag2, MPI_COMM_WORLD, status, ierr)
!             write(*,*) myrankglobal, " received from ",tOctal%mpiThread(tSubcell)
             if (.not.associated(thisOctal%tempStorage)) then
                if (.not.doJustGrav) then
                   allocate(thisOctal%tempStorage(1:thisOctal%maxChildren,1:8))
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
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    logical :: sendLoop
    real(double) :: loc(3), tempStorage(8)
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
          
          if (.not.doJustGrav) then
             tempstorage(1) = thisOctal%rho(Subcell)
             tempStorage(2) = thisOctal%rhoE(Subcell)
             tempStorage(3) = thisOctal%rhou(Subcell)
             tempStorage(4) = thisOctal%rhov(Subcell)
             tempStorage(5) = thisOctal%rhow(Subcell)
             tempStorage(6) = thisOctal%energy(Subcell)
             tempStorage(7) = thisOctal%pressure_i(Subcell)
             tempStorage(8) = thisOctal%phi_i(Subcell)
             call MPI_SEND(tempStorage, 8, MPI_DOUBLE_PRECISION, receiveThread, tag2, MPI_COMM_WORLD, ierr)
          else
             tempStorage(1) = thisOctal%phi_i(Subcell)
             call MPI_SEND(tempStorage, 8, MPI_DOUBLE_PRECISION, receiveThread, tag2, MPI_COMM_WORLD, ierr)
          endif
       endif
    enddo
!    write(*,*) myrankGlobal, " leaving receive requests ", sendLoop
  end subroutine periodBoundaryReceiveRequests


  subroutine periodBoundaryReceiveRequestsLevel(grid, receiveThread, nDepth, doJustGrav)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    logical :: sendLoop
    integer :: nDepth
    real(double) :: loc(3), tempStorage(8)
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
             tempStorage(8) = thisOctal%phi_i(Subcell)
             call MPI_SEND(tempStorage, 8, MPI_DOUBLE_PRECISION, receiveThread, tag2, MPI_COMM_WORLD, ierr)
          else
             tempStorage(1) = thisOctal%phi_i(Subcell)
             call MPI_SEND(tempStorage, 8, MPI_DOUBLE_PRECISION, receiveThread, tag2, MPI_COMM_WORLD, ierr)
          endif
       endif
    enddo
!    write(*,*) myrankGlobal, " leaving receive requests ", sendLoop
  end subroutine periodBoundaryReceiveRequestsLevel



  subroutine dumpValuesAlongLine(grid, thisFile, startPoint, endPoint, nPoints)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    integer :: nPoints
    type(VECTOR) :: startPoint, endPoint, position, direction, cen
    real(double) :: loc(3), rho, rhou , rhoe, p
    character(len=*) :: thisFile
    integer :: ierr
    integer, parameter :: nStorage = 8
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
          rhou = tempStorage(6)
          rhoe = tempStorage(7)
          p = tempStorage(8)
          write(20,'(5e14.5)') modulus(cen-startPoint), rho, rhou/rho, rhoe,p
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
             tempStorage(6) = thisOctal%rhou(subcell)             
             tempStorage(7) = thisOctal%rhoe(subcell)             
             tempStorage(8) = thisOctal%pressure_i(subcell)             
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
    
    if (myrankGlobal == 1) then
       if (filename(1:1) == '*') then
          UN = 6   ! prints on screen
       else
          UN = 69
          open(unit=UN, file = TRIM(filename), status = 'replace',form='formatted')
       end if
    endif



    nOctals=0; nVoxels=0
    if (myrankglobal /= 0) call countSubcellsMPI(thisgrid, nVoxels)
    call MPI_REDUCE(thisgrid%maxDepth, tempInt,1,MPI_INTEGER, MPI_MAX, 1, MPI_COMM_WORLD, ierr)
    maxDepth = tempInt(1)
    call MPI_REDUCE(thisgrid%halfSmallestSubcell, tempDouble,1,MPI_DOUBLE_PRECISION, MPI_MIN, 1, MPI_COMM_WORLD, ierr)
    halfSmallestSubcell = tempDouble(1)
    
    
    if (myRankGlobal == 1) then
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
    if (myrankGlobal == 1) then
       if (filename(1:1) /= '*')  close(UN)
    endif
    
  end subroutine grid_info_mpi
  
  subroutine addNewChildWithInterp(parent, iChild, grid, constantGravity)
    use input_variables, only : maxDepthAMR
    use octal_mod, only: subcellRadius
    include 'mpif.h'
    type(OCTAL), pointer :: parent, thisOctal
    integer :: iChild
    logical, optional :: constantGravity
    type(GRIDTYPE) :: grid
    integer :: nChildren
    integer :: newChildIndex
    integer :: i, iCorner, iDir, nCorner, nDir
    integer :: nd, iSubcell, parentSubcell
    type(VECTOR) :: dir(8), corner(8), position, rVec, testvec
    real(double) :: rhoCorner(8)
    real(double) :: rhoeCorner(8)
    real(double) :: rhouCorner(8)
    real(double) :: rhovCorner(8)
    real(double) :: rhowCorner(8)
    real(double) :: eCorner(8)
    real(double) :: phiCorner(8)
    real(double) :: weight, totalWeight
    real(double) :: rho, rhoe, rhou, rhov, rhow, r, energy, phi
    real(double) :: x1, x2, y1, y2, z1, z2, u, v, w, x, y, z
    real(double) :: oldMass, newMass, factor
    real(double) :: oldEnergy, newEnergy!, oldMom, newMom
    logical, save :: firstTime = .true.
    logical :: debug

    debug = .false.

    testVec = VECTOR(0.290d0, 0.d0, 0.253d0)
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



    if ( ((parent%twoD)  .and.((nThreadsGlobal - 1) == 4)) .or. &
         ((parent%threed).and.((nThreadsGlobal - 1) == 8)).or. &
         ((parent%oneD)  .and.((nThreadsGlobal - 1) == 2)) ) then
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

    if ((parent%threed).and.(nThreadsGlobal - 1) == 64) then
       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 8
             parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
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
    thisOctal%gamma = parent%gamma(parentSubcell)
    thisOctal%iEquationOfState = parent%iEquationofState(parentSubcell)

    if (associated(parent%ionFrac)) then
       do iSubcell = 1, thisOctal%maxChildren
          thisOctal%ionFrac(isubcell,:) = parent%ionFrac(parentSubcell,:)
       enddo
    endif



    if (thisOctal%threed) then
       nDir = 8
       r = 0.1d0*grid%halfSmallestSubcell
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
       corner(1) = thisOctal%centre + VECTOR(-r, -r, -r)
       corner(2) = thisOctal%centre + VECTOR(+r, -r, -r)
       corner(3) = thisOctal%centre + VECTOR(-r, -r, +r)
       corner(4) = thisOctal%centre + VECTOR(+r, -r, +r)
       corner(5) = thisOctal%centre + VECTOR(-r, +r, -r)
       corner(6) = thisOctal%centre + VECTOR(+r, +r, -r)
       corner(7) = thisOctal%centre + VECTOR(-r, +r, +r)
       corner(8) = thisOctal%centre + VECTOR(+r, +r, +r)
    endif

    if (thisOctal%twod) then
       nDir = 4
       r = 0.1d0*grid%halfSmallestSubcell
       dir(1) = VECTOR(-r, 0.d0, -r)
       dir(2) = VECTOR(+r, 0.d0, -r)
       dir(3) = VECTOR(+r, 0.d0, +r)
       dir(4) = VECTOR(-r, 0.d0, +r)
       
       nCorner = 4
       r = thisOctal%subcellSize
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
       
       nCorner = 2
       r = thisOctal%subcellSize
       corner(1) = thisOctal%centre + VECTOR(-r, 0.d0, 0.d0)
       corner(2) = thisOCtal%centre + VECTOR(+r, 0.d0, 0.d0)
    endif

    rhoCorner = 0.d0
    rhoeCorner = 0.d0
    eCorner = 0.d0
    rhouCorner = 0.d0
    rhovCorner = 0.d0
    rhowCorner = 0.d0

    if (debug) then
       write(*,*) "addnewchild with interp debug"
    endif
    do iCorner = 1, nCorner
       totalWeight = 0.d0
       do iDir = 1, nDir
          position = corner(iCorner) + dir(iDir)
          if (inOctal(grid%octreeRoot, position).and.(.not.inSubcell(parent, parentSubcell, position))) then
             call getHydroValues(grid, position, nd, rho, rhoe, rhou, rhov, rhow, energy, phi)
             weight = abs(parent%ndepth - nd)+1.d0

             totalWeight = totalWeight + weight
             rhoCorner(iCorner) = rhoCorner(iCorner) + weight * rho
             rhoeCorner(iCorner) = rhoeCorner(iCorner) + weight * rhoe
             rhouCorner(iCorner) = rhouCorner(iCorner) + weight * rhou
             if (debug) write(*,*) "from get hydro values corner ", icorner, " rhou ", rhou
             rhovCorner(iCorner) = rhovCorner(iCorner) + weight * rhov
             rhowCorner(iCorner) = rhowCorner(iCorner) + weight * rhow
             eCorner(iCorner) = eCorner(iCorner) + weight * energy
             phiCorner(iCorner) = phiCorner(iCorner) + weight * phi
          else
             weight = 1.d0
             totalWeight = totalWeight + weight
             rhoCorner(iCorner) = rhoCorner(iCorner) + parent%rho(parentSubcell)
             rhoeCorner(iCorner) = rhoeCorner(iCorner) + parent%rhoe(parentSubcell)
             rhouCorner(iCorner) = rhouCorner(iCorner) + parent%rhou(parentSubcell)
             if (debug) write(*,*) "from parent values corner ", icorner, " rhou ", parent%rhou(parentsubcell)
             rhovCorner(iCorner) = rhovCorner(iCorner) + parent%rhov(parentSubcell)
             rhowCorner(iCorner) = rhowCorner(iCorner) + parent%rhow(parentSubcell)
             eCorner(iCorner) = eCorner(iCorner) + parent%energy(parentSubcell)
             phiCorner(iCorner) = phiCorner(iCorner) + parent%phi_i(parentSubcell)
          endif
       enddo
       rhoCorner(iCorner) = rhoCorner(iCorner) / totalWeight
       rhoeCorner(iCorner) = rhoeCorner(iCorner) / totalWeight
       rhouCorner(iCorner) = rhouCorner(iCorner) / totalWeight
       rhovCorner(iCorner) = rhovCorner(iCorner) / totalWeight
       rhowCorner(iCorner) = rhowCorner(iCorner) / totalWeight
       eCorner(iCorner) = eCorner(iCorner) / totalWeight
       phiCorner(iCorner) = phiCorner(iCorner) / totalWeight
    enddo
    x1 = thisOctal%xMin
    x2 = thisOctal%xMax
    y1 = thisOctal%yMin
    y2 = thisOctal%yMax
    z1 = thisOctal%zMin
    z2 = thisOctal%zMax

    thisOctal%changed = .true.

    do iSubcell = 1, thisOctal%maxChildren

       if (thisOctal%threed) then
          rVec = subcellcentre(thisOctal, iSubcell)
          x = rVec%x
          y = rVec%y
          z = rVec%z
          u = (x - x1)/(x2 - x1)
          v = (y - y1)/(y2 - y1)
          w = (z - z1)/(z2 - z1)
          thisOctal%rho(iSubcell) = (1.d0 - u) * (1.d0 - v) * (1.d0 - w) * rhoCorner(1) + &
                                    (       u) * (1.d0 - v) * (1.d0 - w) * rhoCorner(2) + &
                                    (1.d0 - u) * (1.d0 - v) * (       w) * rhoCorner(3) + &
                                    (       u) * (1.d0 - v) * (       w) * rhoCorner(4) + &
                                    (1.d0 - u) * (       v) * (1.d0 - w) * rhoCorner(5) + &
                                    (       u) * (       v) * (1.d0 - w) * rhoCorner(6) + &
                                    (1.d0 - u) * (       v) * (       w) * rhoCorner(7) + &
                                    (       u) * (       v) * (       w) * rhoCorner(8) 

          thisOctal%rhoe(iSubcell) = (1.d0 - u) * (1.d0 - v) * (1.d0 - w) * rhoeCorner(1) + &
                                     (       u) * (1.d0 - v) * (1.d0 - w) * rhoeCorner(2) + &
                                     (1.d0 - u) * (1.d0 - v) * (       w) * rhoeCorner(3) + &
                                     (       u) * (1.d0 - v) * (       w) * rhoeCorner(4) + &
                                     (1.d0 - u) * (       v) * (1.d0 - w) * rhoeCorner(5) + &
                                     (       u) * (       v) * (1.d0 - w) * rhoeCorner(6) + &
                                     (1.d0 - u) * (       v) * (       w) * rhoeCorner(7) + &
                                     (       u) * (       v) * (       w) * rhoeCorner(8) 

          thisOctal%rhou(iSubcell) = (1.d0 - u) * (1.d0 - v) * (1.d0 - w) * rhouCorner(1) + &
                                    (       u) * (1.d0 - v) * (1.d0 - w) * rhouCorner(2) + &
                                    (1.d0 - u) * (1.d0 - v) * (       w) * rhouCorner(3) + &
                                    (       u) * (1.d0 - v) * (       w) * rhouCorner(4) + &
                                    (1.d0 - u) * (       v) * (1.d0 - w) * rhouCorner(5) + &
                                    (       u) * (       v) * (1.d0 - w) * rhouCorner(6) + &
                                    (1.d0 - u) * (       v) * (       w) * rhouCorner(7) + &
                                    (       u) * (       v) * (       w) * rhouCorner(8) 

          thisOctal%rhov(iSubcell) = (1.d0 - u) * (1.d0 - v) * (1.d0 - w) * rhovCorner(1) + &
                                    (       u) * (1.d0 - v) * (1.d0 - w) * rhovCorner(2) + &
                                    (1.d0 - u) * (1.d0 - v) * (       w) * rhovCorner(3) + &
                                    (       u) * (1.d0 - v) * (       w) * rhovCorner(4) + &
                                    (1.d0 - u) * (       v) * (1.d0 - w) * rhovCorner(5) + &
                                    (       u) * (       v) * (1.d0 - w) * rhovCorner(6) + &
                                    (1.d0 - u) * (       v) * (       w) * rhovCorner(7) + &
                                    (       u) * (       v) * (       w) * rhovCorner(8) 

          thisOctal%rhow(iSubcell) = (1.d0 - u) * (1.d0 - v) * (1.d0 - w) * rhowCorner(1) + &
                                    (       u) * (1.d0 - v) * (1.d0 - w) * rhowCorner(2) + &
                                    (1.d0 - u) * (1.d0 - v) * (       w) * rhowCorner(3) + &
                                    (       u) * (1.d0 - v) * (       w) * rhowCorner(4) + &
                                    (1.d0 - u) * (       v) * (1.d0 - w) * rhowCorner(5) + &
                                    (       u) * (       v) * (1.d0 - w) * rhowCorner(6) + &
                                    (1.d0 - u) * (       v) * (       w) * rhowCorner(7) + &
                                    (       u) * (       v) * (       w) * rhowCorner(8) 

          thisOctal%energy(iSubcell) = (1.d0 - u) * (1.d0 - v) * (1.d0 - w) * eCorner(1) + &
                                    (       u) * (1.d0 - v) * (1.d0 - w) * eCorner(2) + &
                                    (1.d0 - u) * (1.d0 - v) * (       w) * eCorner(3) + &
                                    (       u) * (1.d0 - v) * (       w) * eCorner(4) + &
                                    (1.d0 - u) * (       v) * (1.d0 - w) * eCorner(5) + &
                                    (       u) * (       v) * (1.d0 - w) * eCorner(6) + &
                                    (1.d0 - u) * (       v) * (       w) * eCorner(7) + &
                                    (       u) * (       v) * (       w) * eCorner(8) 

          thisOctal%phi_i(iSubcell) = (1.d0 - u) * (1.d0 - v) * (1.d0 - w) * phiCorner(1) + &
                                    (       u) * (1.d0 - v) * (1.d0 - w) * phiCorner(2) + &
                                    (1.d0 - u) * (1.d0 - v) * (       w) * phiCorner(3) + &
                                    (       u) * (1.d0 - v) * (       w) * phiCorner(4) + &
                                    (1.d0 - u) * (       v) * (1.d0 - w) * phiCorner(5) + &
                                    (       u) * (       v) * (1.d0 - w) * phiCorner(6) + &
                                    (1.d0 - u) * (       v) * (       w) * phiCorner(7) + &
                                    (       u) * (       v) * (       w) * phiCorner(8) 

       endif
       if (thisOctal%twod) then
          rVec = subcellcentre(thisOctal, iSubcell)
          x = rVec%x
          z = rVec%z
          u = (x - x1)/(x2 - x1)
          v = (z - z1)/(z2 - z1)
          thisOctal%rho(iSubcell) = (1.d0 - u) * (1.d0 - v) * rhoCorner(1) + &
                                 (       u) * (1.d0 - v) * rhoCorner(2) + &
                                 (1.d0 - u) * (       v) * rhoCorner(3) + &
                                 (       u) * (       v) * rhoCorner(4)
          
          
          thisOctal%rhoe(iSubcell) = (1.d0 - u) * (1.d0 - v) * rhoeCorner(1) + &
               (       u) * (1.d0 - v) * rhoeCorner(2) + &
               (1.d0 - u) * (       v) * rhoeCorner(3) + &
               (       u) * (       v) * rhoeCorner(4)
          
          
          thisOctal%rhou(iSubcell) = (1.d0 - u) * (1.d0 - v) * rhouCorner(1) + &
               (       u) * (1.d0 - v) * rhouCorner(2) + &
               (1.d0 - u) * (       v) * rhouCorner(3) + &
               (       u) * (       v) * rhouCorner(4)
          if (debug) write(*,*) "interped rhou ",thisOctal%rhou(iSubcell), " u ", u, " v ",v, &
               " corners ",rhoucorner(1:4)
          
          thisOctal%rhov(iSubcell) = (1.d0 - u) * (1.d0 - v) * rhovCorner(1) + &
               (       u) * (1.d0 - v) * rhovCorner(2) + &
               (1.d0 - u) * (       v) * rhovCorner(3) + &
               (       u) * (       v) * rhovCorner(4)
          
          thisOctal%rhow(iSubcell) = (1.d0 - u) * (1.d0 - v) * rhowCorner(1) + &
               (       u) * (1.d0 - v) * rhowCorner(2) + &
               (1.d0 - u) * (       v) * rhowCorner(3) + &
               (       u) * (       v) * rhowCorner(4)
          
          thisOctal%energy(iSubcell) = (1.d0 - u) * (1.d0 - v) * eCorner(1) + &
               (       u) * (1.d0 - v) * eCorner(2) + &
               (1.d0 - u) * (       v) * eCorner(3) + &
               (       u) * (       v) * eCorner(4)


          thisOctal%phi_i(iSubcell) = (1.d0 - u) * (1.d0 - v) * phiCorner(1) + &
               (       u) * (1.d0 - v) * phiCorner(2) + &
               (1.d0 - u) * (       v) * phiCorner(3) + &
               (       u) * (       v) * phiCorner(4)
       endif

       if (thisOctal%oned) then
          rVec = subcellcentre(thisOctal, iSubcell)
          x = rVec%x
          u = (x - x1)/(x2 - x1)
          thisOctal%rho(iSubcell) = (1.d0 - u) * rhoCorner(1) + &
                                    (       u) * rhoCorner(2)

          thisOctal%rhoe(iSubcell) = (1.d0 - u) * rhoeCorner(1) + &
                                    (       u) * rhoeCorner(2)


          thisOctal%rhou(iSubcell) = (1.d0 - u) * rhouCorner(1) + &
                                    (       u) *  rhouCorner(2)

          thisOctal%rhov(iSubcell) = (1.d0 - u) * rhovCorner(1) + &
                                    (       u) *  rhovCorner(2)

          thisOctal%rhow(iSubcell) = (1.d0 - u) * rhowCorner(1) + &
                                    (       u) *  rhowCorner(2)

          thisOctal%energy(iSubcell) = (1.d0 - u) * eCorner(1) + &
                                    (       u) *  eCorner(2)
       endif

    enddo
    
    ! conservation normalizations

    ! mass
    oldMass = parent%rho(parentSubcell) * cellVolume(parent, parentSubcell)
    newMass = 0.d0
    do iSubcell = 1, thisOctal%maxChildren
       newMass = newMass + thisOctal%rho(isubcell) * cellVolume(thisOctal, iSubcell)
    enddo
    factor = oldMass / newMass
    thisOctal%rho(1:thisOctal%maxChildren) = thisOctal%rho(1:thisOctal%maxChildren) * factor
    if ( associated (thisOctal%nh) ) thisOctal%nh(1:thisOctal%maxChildren) = thisOctal%rho(1:thisOctal%maxChildren)/mHydrogen


    ! energy

    oldEnergy = parent%rhoe(parentSubcell) * cellVolume(parent, parentSubcell)
    newEnergy = 0.d0
    do iSubcell = 1, thisOctal%maxChildren
       newEnergy = newEnergy + thisOctal%rhoe(isubcell) * cellVolume(thisOctal, iSubcell)
    enddo
    factor = oldEnergy/newEnergy
    thisOctal%rhoe(1:thisOctal%maxChildren) = thisOctal%rhoe(1:thisOctal%maxChildren) * factor
!
!    ! momentum (u)
!
!    oldMom = parent%rhou(parentSubcell)
!    newMom = SUM(thisOctal%rhou(1:thisOctal%maxChildren))/dble(thisOctal%maxChildren)
!    if (newMom /= 0.d0) then
!       factor = oldMom / newMom
!       thisOctal%rhou(1:thisOctal%maxChildren) = thisOctal%rhou(1:thisOctal%maxChildren) * factor
!    endif
!    
!
!    ! momentum (v)
!
!    oldMom = parent%rhov(parentSubcell)
!    newMom = SUM(thisOctal%rhov(1:thisOctal%maxChildren))/dble(thisOctal%maxChildren)
!    if (newMom /= 0.d0) then
!       factor = oldMom / newMom
!       thisOctal%rhov(1:thisOctal%maxChildren) = thisOctal%rhov(1:thisOctal%maxChildren) * factor
!    endif
!
!    ! momentum (w)
!
!    oldMom = parent%rhow(parentSubcell)
!    newMom = SUM(thisOctal%rhow(1:thisOctal%maxChildren))/dble(thisOctal%maxChildren)
!    if (newMom /= 0.d0) then
!       factor = oldMom / newMom
!       thisOctal%rhow(1:thisOctal%maxChildren) = thisOctal%rhow(1:thisOctal%maxChildren) * factor
!    endif

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
666 continue
  end subroutine addNewChildWithInterp

  subroutine shutdownServers()
    include 'mpif.h'
    integer :: iThread
    real(double) :: loc(3)
    integer, parameter :: tag = 50
    integer :: ierr

    do iThread = 1, nThreadsGlobal-1
       if (iThread /= myrankGlobal) then
          loc = 1.d30
          call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
       endif
    enddo
  end subroutine shutdownServers

  subroutine getHydroValues(grid, position, nd, rho, rhoe, rhou, rhov, rhow, energy, phi)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    integer, intent(out) :: nd
    real(double), intent(out) :: rho, rhoe, rhou, rhov, rhow, energy, phi
    type(VECTOR) :: position
    real(double) :: loc(3)
    type(OCTAL), pointer :: thisOctal, parent
    integer :: iThread
    integer, parameter :: nStorage = 8
    real(double) :: tempStorage(nStorage)
    integer :: subcell
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr

    thisOctal => grid%octreeRoot
    call findSubcellLocal(position, thisOctal, subcell)
    
    if (octalOnThread(thisOctal, subcell, myrankGlobal)) then
       if (.not.thisOctal%changed(subcell)) then
          rho = thisOctal%rho(subcell)
          rhoe = thisOctal%rhoe(subcell)
          rhou = thisOctal%rhou(subcell)
          rhov = thisOctal%rhov(subcell)
          rhow = thisOctal%rhow(subcell)
          nd = thisOctal%nDepth
          energy = thisOctal%energy(subcell)
          phi = thisOctal%phi_i(subcell)
       else
          parent => thisOctal%parent
          rho =  parent%rho(thisOctal%parentsubcell)
          rhoe = parent%rhoe(thisOctal%parentsubcell)
          rhou = parent%rhou(thisOctal%parentsubcell)
          rhov = parent%rhov(thisOctal%parentsubcell)
          rhow = parent%rhow(thisOctal%parentsubcell)
          nd = parent%nDepth
          energy = parent%energy(thisOctal%parentSubcell)
          phi = parent%phi_i(thisOctal%parentSubcell)
       endif
    else
       iThread = thisOctal%mpiThread(subcell)
       loc(1) = position%x
       loc(2) = position%y
       loc(3) = position%z
       call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
       call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, status, ierr)
       nd = nint(tempStorage(1))
       rho = tempStorage(2)
       rhoe = tempStorage(3)
       rhou = tempStorage(4)
       rhov = tempStorage(5)
       rhow = tempStorage(6)
       energy = tempStorage(7)
       phi = tempStorage(8)
    endif
  end subroutine getHydroValues

  subroutine hydroValuesServer(grid, iThread)
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    logical :: stillServing
    real(double) :: loc(3)
    type(VECTOR) :: position
    type(OCTAL), pointer :: thisOctal, parent
    integer :: subcell
    integer :: iThread
    integer, parameter :: nStorage = 8
    real(double) :: tempStorage(nStorage)
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 50
    integer :: ierr

    stillServing = .true.

    thisOctal => grid%octreeroot
    do while (stillServing)
       call MPI_RECV(loc, 3, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, status, ierr)
       position%x = loc(1)
       position%y = loc(2)
       position%z = loc(3)
       if (position%x > 1.d29) then
          stillServing= .false.
       else
          call findSubcellLocal(position, thisOctal, subcell)
          if (.not.thisOctal%changed(subcell)) then

             tempStorage(1) = thisOctal%nDepth
             tempStorage(2) = thisOctal%rho(subcell)
             tempStorage(3) = thisOctal%rhoe(subcell)             
             tempStorage(4) = thisOctal%rhou(subcell)             
             tempStorage(5) = thisOctal%rhov(subcell)             
             tempStorage(6) = thisOctal%rhow(subcell)        
             tempStorage(7) = thisOctal%energy(subcell)
             tempStorage(8) = thisOctal%phi_i(subcell)
          else
             parent => thisOctal%parent
             tempStorage(1) = parent%nDepth
             tempStorage(2) = parent%rho(thisOctal%parentSubcell)
             tempStorage(3) = parent%rhoe(thisOctal%parentsubcell)             
             tempStorage(4) = parent%rhou(thisOctal%parentsubcell)             
             tempStorage(5) = parent%rhov(thisOctal%parentsubcell)             
             tempStorage(6) = parent%rhow(thisOctal%parentsubcell)        
             tempStorage(7) = parent%energy(thisOctal%parentsubcell)        
             tempStorage(8) = parent%phi_i(thisOctal%parentsubcell)        
          endif
          call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD, ierr)
       endif
    enddo

  end subroutine hydroValuesServer

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
