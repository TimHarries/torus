module mpi_amr_mod
#ifdef MPI

  use kind_mod
  use amr_mod
  use mpi_global_mod
  implicit none

contains

  subroutine setupAMRCOMMUNICATOR
    use mpi
    use inputs_mod, only : nHydroThreadsinput, splitOverMPI
    integer :: ierr, i, j
    integer, allocatable :: ranks(:)
    integer :: worldGroup, amrGroup, localWorldGroup
    integer :: amrParallelGroup, zeroThreadGroup

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreadsGlobal, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRankWorldGlobal, ierr)
    nHydroSetsGlobal = 1
    myHydroSetGlobal = 0
    nHydroThreadsGlobal = nHydroThreadsinput
    if (splitOverMPI) then

       if (nHydroThreadsInput == 0) then

          if (mod(nThreadsGlobal, 513) == 0) then
             nHydroThreadsGlobal = 512
          else if (mod(nThreadsGlobal, 65) == 0) then
             nHydroThreadsGlobal = 64
          else if (mod(nThreadsGlobal, 9) == 0) then
             nHydroThreadsGlobal = 8
          else if (mod(nThreadsGlobal, 17) == 0) then
             nHydroThreadsGlobal = 16
          else if (mod(nThreadsGlobal, 5) == 0) then
             nHydroThreadsGlobal = 4
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
       
       call MPI_COMM_GROUP(MPI_COMM_WORLD, worldGroup, ierr)

       allocate(ranks(1:(nHydroThreadsGlobal+1)))
       do i = 1, nHydroThreadsGlobal+1
          ranks(i) = myHydroSetGlobal * (nHydroThreadsGlobal+1) + i - 1
       enddo
       call MPI_GROUP_INCL(worldGroup, nHydroThreadsGlobal+1, ranks, localWorldGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, localWorldGroup, localWorldCOMMUNICATOR, ierr)
       deallocate(ranks)

       allocate(ranks(1:1))
       ranks(1) = 0
       call MPI_GROUP_EXCL(localWorldGroup, 1, ranks, amrGroup, ierr)
       call MPI_COMM_CREATE(MPI_COMM_WORLD, amrGroup, amrCOMMUNICATOR, ierr)
       deallocate(ranks)
       call MPI_COMM_RANK(localWorldCommunicator, myRankGlobal, ierr)

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

          
          

    else
       myrankGlobal = myrankWorldGlobal
    endif
  end subroutine setupAMRCOMMUNICATOR

! Free communicator created by setupAMRCOMMUNICATOR
  subroutine freeAMRCOMMUNICATOR
    use mpi
    implicit none

    integer :: ierr

    if ( myRankGlobal /= 0 ) then 
       call MPI_COMM_FREE(amrCOMMUNICATOR, ierr)
    end if

  end subroutine freeAMRCOMMUNICATOR

  subroutine findMassOverAllThreads(grid, mass)
    use mpi
    type(GRIDTYPE) :: grid
    real(double), intent(out) :: mass
    real(double), allocatable :: massOnThreads(:), temp(:), volumeOnThreads(:)
    integer :: ierr
!    real(double) :: totalVolume

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
    use inputs_mod, only : hydrodynamics
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
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
                      dv = thisOctal%subcellSize**2
                   else if (thisOctal%oned) then
                      dv = thisOctal%subcellSize                      

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
                cellCentre = subcellCentre(thisOctal, subcell)
                vel = VECTOR(thisOctal%rhou(subcell)/thisOctal%rho(subcell), &
                     thisOctal%rhov(subcell)/thisOctal%rho(subcell), &
                     thisOctal%rhow(subcell)/thisOctal%rho(subcell))
                dv = cellVolume(thisOctal, subcell)*1.d30
                dm = thisOctal%rho(subcell) * dv
                rVec = 1.d10*(centre-cellCentre)
                angMom = angMom + (rVec.cross.(dm*vel))
             endif
          endif
       end if
    enddo
  end subroutine findAngMomMPI

  subroutine findEnergyOverAllThreads(grid, energy)
    use mpi
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
          if(.not. thisoctal%ghostcell(subcell)) then
             if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
                if (thisOctal%threed) then
!                   dv = cellVolume(thisOctal, subcell) * 1.d30
                   dv = thisOctal%subcellSize**3
                else if (thisOctal%twoD) then
                   dv = thisOctal%subcellSize**2
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
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal, tOctal
    type(octal), pointer  :: neighbourOctal
    character(len=*) :: boundaryType
    integer :: receiveThread, sendThread, tsubcell
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: nOctals
    integer, parameter :: nStorage = 20
    real(double) :: loc(3), tempStorage(nStorage)
    type(VECTOR) :: octVec, direction, rVec, pVec
    integer :: nBound
    integer :: iOctal
    integer :: subcell, neighbourSubcell
    integer :: tag = 77
    integer :: status(MPI_STATUS_SIZE)
    logical :: sendLoop
    integer :: nDepth
    integer :: ierr
    real(double) :: q , rho, rhoe, rhou, rhov, rhow, pressure, phi, flux, phigas, q11, q22, q33

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
       do iOctal =  1, nOctals
          
          thisOctal => octalArray(iOctal)%content
          
          do subcell = 1, thisOctal%maxChildren             
             if (.not.thisOctal%hasChild(subcell)) then
                
!                if (thisOctal%mpiThread(subcell) /= myRank) cycle
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
                
                neighbourOctal => thisOctal
                call findSubcellLocal(octVec, neighbourOctal, neighbourSubcell)

                if (octalOnThread(neighbourOctal, neighbourSubcell, receiveThread)) cycle

                if (.not.octalOnThread(neighbourOctal, neighbourSubcell, sendThread)) then
                   write(*,*) "Error 1 Neighbour on ",boundaryType, " of ", myrankglobal, &
                        "  is not on thread ", sendThread, " but ", &
                   neighbourOctal%mpiThread(neighboursubcell), " depth ",neighbourOctal%ndepth
                   stop
                endif

                loc(1) = octVec%x
                loc(2) = octVec%y
                loc(3) = octVec%z
!                write(*,*) myRank, " has identified a boundary cell ", loc(1:3)

                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_SEND(thisOctal%nDepth, 1, MPI_INTEGER, sendThread, tag, localWorldCommunicator, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
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

             tempStorage(15) = pVec%x
             tempStorage(16) = pVec%y
             tempStorage(17) = pVec%z

             if (neighbourOctal%mpiThread(neighboursubcell) /= sendthread) then
                write(*,*) "trying to send on ",boundaryType, " but is not on thread ", sendThread
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

                tempStorage(13) = neighbourOctal%phi_gas(neighbourSubcell)
                tempStorage(14) = neighbourOctal%x_i_minus_1(neighbourSubcell)


                tempStorage(18) = neighbourOctal%qViscosity(neighbourSubcell,1,1)
                tempStorage(19) = neighbourOctal%qViscosity(neighbourSubcell,2,2)
                tempStorage(20) = neighbourOctal%qViscosity(neighbourSubcell,3,3)

!                write(*,*) myrank," set up tempstorage with ", &
!                     tempstorage(1:nStorage),neighbourOctal%nDepth, neighbourSubcell,neighbourOctal%ghostCell(neighbourSubcell), &
!                     neighbourOctal%edgeCell(neighbourSubcell)

             else ! need to average
                call averageValue(direction, neighbourOctal,  neighbourSubcell, q, rhou, rhov, rhow, rho, rhoe, pressure, &
                     flux, phi, phigas, q11, q22, q33)
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
                tempStorage(13) = phigas
                tempstorage(14) = neighbourOctal%x_i_minus_1(neighbourSubcell)
                tempstorage(18) = q11
                tempstorage(19) = q22
                tempstorage(20) = q33
             endif
!                          write(*,*) myRank, " sending temp storage ", tempStorage(1:nStorage)
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, ierr)
!                          write(*,*) myRank, " temp storage sent"
             
          endif
       end do
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
!"
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
    integer :: ierr
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal, tOctal
    type(octal), pointer  :: neighbourOctal
    character(len=*) :: boundaryType
    integer :: receiveThread, sendThread, tsubcell
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: nOctals
    integer, parameter :: nStorage = 20
    real(double) :: loc(3), tempStorage(nStorage)
    type(VECTOR) :: octVec, direction, rVec, pVec
    integer :: nBound
    integer :: iOctal
    integer :: nDepth, thisnDepth
    integer :: subcell, neighbourSubcell
    integer :: tag = 77
    integer :: status(MPI_STATUS_SIZE)
    logical :: sendLoop

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
       do iOctal =  1, nOctals
          
          thisOctal => octalArray(iOctal)%content
!          write(*,*) myrank, " doing octal of ", thisOctal%ndepth, " depth "
          
          do subcell = 1, thisOctal%maxChildren
             
                
!                if (thisOctal%mpiThread(subcell) /= myRank) cycle
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
                
                neighbourOctal => thisOctal
                call findSubcellLocalLevel(octVec, neighbourOctal, neighbourSubcell, nDepth)

                if (octalOnThread(neighbourOctal, neighbourSubcell, receiveThread)) cycle

                if (.not.octalOnThread(neighbourOctal, neighbourSubcell, sendThread)) then
                   write(*,*) "Error 2 Neighbour on ",boundaryType, " of ", myrankglobal, &
                        "  is not on thread ", sendThread, " but ", &
                   neighbourOctal%mpiThread(neighboursubcell), " depth ",neighbourOctal%ndepth
                   stop
                endif

                loc(1) = octVec%x
                loc(2) = octVec%y
                loc(3) = octVec%z
!                write(*,*) myRank, " has identified a boundary cell ", loc(1:3)

                call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_SEND(thisOctal%nDepth, 1, MPI_INTEGER, sendThread, tag, localWorldCommunicator, ierr)
!                write(*,*) myrank, " sent the locator to ", sendThread

                call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)
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

             call MPI_RECV(thisnDepth, 1, MPI_INTEGER, receiveThread, tag, localWorldCommunicator, status, ierr)

             call findSubcellTDLevel(octVec, grid%octreeRoot, neighbourOctal, neighbourSubcell, nDepth)

             pVec = subcellCentre(neighbourOctal, neighbourSubcell)

             if (neighbourOctal%mpiThread(neighboursubcell) /= sendthread) then
                write(*,*) "trying to send on ",boundaryType, " but is not on thread ", sendThread
                stop
             endif

             if (neighbourOctal%nDepth <= thisnDepth) then
                Tempstorage(1) = neighbourOctal%q_i(neighbourSubcell)
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

                tempStorage(13) = neighbourOctal%phi_gas(neighbourSubcell)

                tempStorage(15) = pVec%x
                tempStorage(16) = pVec%y
                tempStorage(17) = pVec%z

                tempStorage(18) = neighbourOctal%qViscosity(neighbourSubcell,1,1)
                tempStorage(19) = neighbourOctal%qViscosity(neighbourSubcell,2,2)
                tempStorage(20) = neighbourOctal%qViscosity(neighbourSubcell,3,3)

!                          write(*,*) myRank, " sending temp storage"
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, receiveThread, tag, localWorldCommunicator, ierr)
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
    use utils_mod, only: indexx
    use mpi
    type(GRIDTYPE) :: grid
    integer, intent(out) :: nPairs, thread1(:), thread2(:), nBound(:), nGroup, group(:)
    integer :: iThread
    integer :: i
    integer, allocatable :: indx(:), itmp(:)
    real, allocatable :: sort(:)
    integer :: list(1000), nList


    nPairs = 0
    do iThread = 1, nHydroThreadsGlobal
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
       rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xplus, px, py, pz, q11, q22, q33)
    use mpi
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal, tOctal
    type(VECTOR) :: direction, rVec, pVec
    integer :: subcell, neighbourSubcell, tSubcell
    integer :: nBound, nDepth
    integer, intent(out) :: nd

    real(double), intent(out) :: q, rho, rhoe, rhou, rhov, rhow, qnext, x, pressure, flux, phi, phigas
    real(double), intent(out) :: xplus, px, py, pz, q11, q22, q33


    nbound = getNboundFromDirection(direction)

    if ((thisOctal%twoD).and.((nBound == 5).or. (nBound == 6))) then
       write(*,*) "Bonndary error for twod: ",nbound
    endif
    if (octalOnThread(neighbourOctal, neighbourSubcell, myRankGlobal)) then

       nd = neighbourOctal%nDepth

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
          q11 = neighbourOctal%qViscosity(neighbourSubcell,1,1)
          q22 = neighbourOctal%qViscosity(neighbourSubcell,2,2)
          q33 = neighbourOctal%qViscosity(neighbourSubcell,3,3)

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
          q11 = neighbourOctal%qViscosity(neighbourSubcell,1,1)
          q22 = neighbourOctal%qViscosity(neighbourSubcell,2,2)
          q33 = neighbourOctal%qViscosity(neighbourSubcell,3,3)

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
               flux, phi, phigas, q11, q22, q33) ! fine to coarse
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
             qNext = neighbourOctal%mpiBoundaryStorage(neighbourSubcell, nBound, 1)
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
          write(*,*) "this centre",subcellCentre(thisOctal, subcell)
          write(*,*) "neig centre",subcellCentre(neighbourOctal, neighboursubcell)
          stop
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
       q11 = thisOctal%mpiBoundaryStorage(subcell, nBound, 18)
       q22 = thisOctal%mpiBoundaryStorage(subcell, nBound, 19)
       q33 = thisOctal%mpiBoundaryStorage(subcell, nBound, 20)

    endif
  end subroutine getNeighbourValues



  subroutine averageValue(direction, neighbourOctal, neighbourSubcell, q, rhou, rhov, rhow, rho, &
       rhoe, pressure, flux, phi, phigas, q11, q22, q33)
    type(OCTAL), pointer ::  neighbourOctal
    integer :: nSubcell(4), neighbourSubcell
    real(double), intent(out) :: q, rho, rhov, rhou, rhow, pressure, flux, rhoe, phi, phigas,q11,q22,q33
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
       phigas = neighbourOctal%phi_gas(neighbourSubcell)
       q11 = neighbourOctal%qViscosity(neighbourSubcell, 1, 1)
       q22 = neighbourOctal%qViscosity(neighbourSubcell, 2, 2)
       q33 = neighbourOctal%qViscosity(neighbourSubcell, 3, 3)
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
       phigas = 0.5d0*(neighbourOctal%phi_gas(nSubcell(1)) + neighbourOctal%phi_gas(nSubcell(2)))


       q11 = 0.5d0*(neighbourOctal%qViscosity(nSubcell(1),1,1) + neighbourOctal%qViscosity(nSubcell(2),1,1))
       q22 = 0.5d0*(neighbourOctal%qViscosity(nSubcell(1),2,2) + neighbourOctal%qViscosity(nSubcell(2),2,2))
       q33 = 0.5d0*(neighbourOctal%qViscosity(nSubcell(1),3,3) + neighbourOctal%qViscosity(nSubcell(2),3,3))

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

       phigas = fac*(neighbourOctal%phi_gas(nSubcell(1)) + neighbourOctal%phi_gas(nSubcell(2)) + & 
            neighbourOctal%phi_gas(nSubcell(3)) + neighbourOctal%phi_gas(nSubcell(4)))


       q11 = fac*(neighbourOctal%qViscosity(nSubcell(1),1,1) + neighbourOctal%qViscosity(nSubcell(2),1,1) + & 
            neighbourOctal%qViscosity(nSubcell(3),1,1) + neighbourOctal%qViscosity(nSubcell(4),1,1))
       q22 = fac*(neighbourOctal%qViscosity(nSubcell(1),2,2) + neighbourOctal%qViscosity(nSubcell(2),2,2) + & 
            neighbourOctal%qViscosity(nSubcell(3),2,2) + neighbourOctal%qViscosity(nSubcell(4),2,2))
       q33 = fac*(neighbourOctal%qViscosity(nSubcell(1),3,3) + neighbourOctal%qViscosity(nSubcell(2),3,3) + & 
            neighbourOctal%qViscosity(nSubcell(3),3,3) + neighbourOctal%qViscosity(nSubcell(4),3,3))

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

  subroutine countSubcellsMPI(grid, nSubcells, nSubcellArray)
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

    call MPI_BARRIER(amrCOMMUNICATOR, ierr)
    call MPI_COMM_RANK(amrCOMMUNICATOR, myRank, ierr)
    call MPI_COMM_SIZE(amrCOMMUNICATOR, nThreads, ierr)

    allocate(iArray(1:nThreads))

    nSubcells = 0
    if (myrank == 0) then
       call countVoxelsMPI(grid%octreeRoot,nVoxels)
       nSubcells = nSubcells + nVoxels
       iArray(1) = nSubcells
       do iThread = 1, nThreads - 1
          call MPI_RECV(n, 1, MPI_INTEGER, iThread, tag, amrCOMMUNICATOR, status, ierr)
          iArray(iThread+1) = n
          nSubcells = nSubcells + n
       end do
    else
       call countVoxelsMPI(grid%octreeRoot,nVoxels)
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
        use mpi
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
    integer :: nFirstLevel

    check = .true.

    if (thisOctal%twoD) then
       if (nHydroThreadsGlobal == 4) then
          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif
       if (nHydroThreadsGlobal == 16) then
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
       if (nHydroThreadsGlobal == 8) then

          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif
       if (nHydroThreadsGlobal == 64) then
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

       if (nHydroThreadsGlobal == 512) then
          nFirstLevel = (myRank-1) / 64 + 1
          if (thisOctal%nDepth == 1) then
             if (thisOctal%mpiThread(subcell) == nFirstLevel) then
                check = .true.
             else
                check = .false.
             endif
          else if (thisOctal%nDepth == 2) then
             nFirstLevel = (myRank-1) / 8  + 1
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
       if (nHydroThreadsGlobal == 2) then
          if (thisOctal%mpiThread(subcell) == myRank) then
             check = .true.
          else
             check = .false.
          endif
       endif
       if (nHydroThreadsGlobal == 4) then
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
            ! write(*,*) "boundary partner ", thisOctal%boundaryPartner(subcell)
            ! write(*,*) myrankGlobal, " sending locator to ", tOctal%mpiThread(tsubcell)
             call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag1, localWorldCommunicator, ierr)
            ! write(*,*) myRankGlobal, " awaiting recv from ", tOctal%mpiThread(tsubcell)
             call MPI_RECV(tempStorage, 8, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), &
                  tag2, localWorldCommunicator, status, ierr)
            ! write(*,*) myrankglobal, " received from ",tOctal%mpiThread(tSubcell)
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

    use mpi
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
             call MPI_SEND(loc, 3, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag1, localWorldCommunicator, ierr)
!             write(*,*) myRankGlobal, " awaiting recv from ", tOctal%mpiThread(tsubcell)
             call MPI_RECV(tempStorage, 7, MPI_DOUBLE_PRECISION, tOctal%mpiThread(tSubcell), tag2, &
                  localWorldCommunicator, status, ierr)
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
    use mpi
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
             call MPI_SEND(tempStorage, 7, MPI_DOUBLE_PRECISION, receiveThread, tag2, localWorldCommunicator, ierr)
          else
!             tempStorage(1) = thisOctal%phi_i(Subcell)
             tempStorage(1) = thisOctal%phi_gas(Subcell)
             call MPI_SEND(tempStorage, 7, MPI_DOUBLE_PRECISION, receiveThread, tag2, localWorldCommunicator, ierr)
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
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, soctal
    integer :: subcell
    integer :: nPoints
    type(VECTOR) :: startPoint, endPoint, position, direction, cen
    real(double) :: loc(3), rho, rhou , rhoe, p, phi_stars, phi_gas
    character(len=*) :: thisFile
    integer :: ierr
    integer, parameter :: nStorage = 10
    real(double) :: tempSTorage(nStorage), tval
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
          write(20,'(1p,7e14.5)') modulus(cen), rho, rhou/rho, rhoe,p, phi_stars, phi_gas
          position = cen
          position = position + (tVal+1.d-3*grid%halfSmallestSubcell)*direction
          !print *, "POSITION 2 ", position
          !print *, "direction 2 ", direction
          !print *, "tVal 2 ", tVal

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
             call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
          endif
       enddo
    endif

666 continue
  end subroutine dumpValuesAlongLine
    
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
          open(unit=UN, file = TRIM(filename), status = 'replace',form='formatted')
       end if
    endif

    nOctals=0; nVoxels=0
    if (myRankGlobal /= 0) call countSubcellsMPI(thisgrid, nVoxels)
!    call MPI_BARRIER(localWorldCommunicator, ierr)
!    print *, "myRankGlobal ", myRankGlobal, "is counting subcells "
!    call countSubcellsMPI(thisgrid, nVoxels)
    call MPI_BARRIER(localWorldCommunicator, ierr)
       call MPI_REDUCE(thisgrid%maxDepth, tempInt, 1, MPI_INTEGER, MPI_MAX, 1, localWorldCommunicator, ierr)
       maxDepth = tempInt(1)
       call MPI_REDUCE(thisgrid%halfSmallestSubcell, tempDouble,1,MPI_DOUBLE_PRECISION, MPI_MIN, 1, localWorldCommunicator, ierr)
       halfSmallestSubcell = tempDouble(1)
    
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
    use inputs_mod, only : maxDepthAMR, minDepthAmr
    use octal_mod, only: subcellRadius
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
    type(VECTOR) :: rVec, centre
    real(double) :: newMom, oldMom
    real(double) :: rhoCorner(8)
    real(double) :: rhoeCorner(8)
    real(double) :: rhouCorner(8)
    real(double) :: rhovCorner(8)
    real(double) :: rhowCorner(8)
    real(double) :: eCorner(8)
    real(double) :: phiCorner(8)
    real(double) :: pressureCorner(8)
    real(double) :: x1, x2, y1, y2, z1, z2, u, x, y, z, dv
    real(double) :: oldMass, newMass, factor, massFactor
    real(double) :: oldEnergy, newEnergy
    integer :: npoints
    integer :: ier
    integer :: nr, nw, nq
    integer, allocatable :: lnext(:), lcell(:,:,:)
    integer, allocatable :: lcell2d(:,:)
    real(double) :: xyzmin(3), xyzdel(3)
    real(double) :: rMax
    real(double), allocatable :: a(:,:), rsq(:)
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
    real(double) :: dx, dz, xmin, zmin
!    integer :: counter
    character(len=80) :: message
    real(double) :: radius
    logical, save :: firstTime = .true.
    logical :: debug

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

    topOctal => thisOctal%parent
    topOctalSubcell = thisOctal%parentsubcell
    do while(topOctal%changed(topOctalSubcell))
       topOctalSubcell = topOctal%parentSubcell
       topOctal => topOctal%parent
    enddo

    if(thisOctal%oneD) then
       call returnAddNewChildCornerArrays(thisOctal, topOctal, topOctalSubcell, grid, rhoCorner, rhoeCorner, rhouCorner, &
            rhovCorner, rhowCorner, eCorner, phiCorner, pressureCorner)
    else 
       call returnAddNewChildPointArrays(thisOctal, topOctal, topOctalSubcell, grid, rhoPoint, rhoePoint, rhouPoint, &
            rhovPoint, rhowPoint, phiPoint, energyPoint, pressurePoint, npoints)
    end if

    x1 = centre%x - topOctal%subcellSize/2.d0
    x2 = centre%x + topOctal%subcellSize/2.d0
    y1 = centre%y - topOctal%subcellSize/2.d0
    y2 = centre%y + topOctal%subcellSize/2.d0
    z1 = centre%z - topOctal%subcellSize/2.d0
    z2 = centre%z + topOctal%subcellSize/2.d0

    do iSubcell = 1, thisOctal%maxChildren

       if (thisOctal%threed) then
          rVec = subcellcentre(thisOctal, iSubcell)
          x = rVec%x
          y = rVec%y
          z = rVec%z

          radius = thisOctal%subcellSize*4.d0
          nPoints = 0

          do while (nPoints < 10)
             call getPointsInRadius(rVec, radius, grid, npoints, rhoPoint, rhoePoint, &
                  rhouPoint, rhovPoint, rhowPoint, energyPoint, pressurePoint, phiPoint, xPoint, yPoint, zPoint)
             radius = radius * 2.d0
          enddo

          nq = min(40, nPoints - 1)
          nw = min(40, nPoints - 1)
          nr = max(1,nint((dble(nPoints)/3.d0)**0.333d0))
          allocate(lCell(1:nr,1:nr,1:nr))
          allocate(lnext(1:nPoints))
          allocate(rsq(1:nPoints))
          allocate(a(9,1:nPoints))

          call qshep3 (nPoints, xPoint, yPoint, zPoint, rhoPoint, nq, nw, nr, lcell, lnext, xyzmin, &
               xyzdel, rmax, rsq, a, ier ) 
          if (ier /= 0) then
             write(message,*) "Qshep3 returned an error ",ier
             call writeWarning(message)
             write(*,*) " npoints ",npoints
             do i = 1, nPoints
                write(*,*) xPoint(i), ypoint(i),zpoint(i)
             enddo
          endif

          thisOctal%rho(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, rhoPoint, nr, lcell, lnext, &
               xyzmin, xyzdel, rmax, rsq, a)
          if (thisOctal%rho(iSubcell) < 0.d0) then
             write(*,*) "Negative density after interp ",thisOctal%rho(iSubcell)
             thisOctal%rho(iSubcell) = SUM(rhoPoint(1:nPoints))/dble(nPoints)
             thisOctal%rhoe(iSubcell) = SUM(rhoePoint(1:nPoints))/dble(nPoints)
             thisOctal%rhou(iSubcell) = SUM(rhouPoint(1:nPoints))/dble(nPoints)
             thisOctal%rhov(iSubcell) = SUM(rhovPoint(1:nPoints))/dble(nPoints)
             thisOctal%rhow(iSubcell) = SUM(rhowPoint(1:nPoints))/dble(nPoints)
             thisOctal%energy(iSubcell) = SUM(energyPoint(1:nPoints))/dble(nPoints)
             thisOctal%phi_gas(iSubcell) = SUM(phiPoint(1:nPoints))/dble(nPoints)
             thisOctal%pressure_i(iSubcell) = SUM(pressurePoint(1:nPoints))/dble(nPoints)

          else
          call qshep3 (nPoints, xPoint, yPoint, zPoint, rhoePoint, nq, nw, nr, lcell, lnext, xyzmin, &
               xyzdel, rmax, rsq, a, ier )
          if (ier /= 0) then
             write(message,*) "Qshep3 returned an error ",ier
             call writeWarning(message)
          endif

          thisOctal%rhoe(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, rhoePoint, nr, lcell, lnext, &
               xyzmin, xyzdel, rmax, rsq, a)

          call qshep3 (nPoints, xPoint, yPoint, zPoint, rhouPoint, nq, nw, nr, lcell, lnext, xyzmin, &
               xyzdel, rmax, rsq, a, ier )
          if (ier /= 0) then
             write(message,*) "Qshep3 returned an error ",ier
             call writeWarning(message)
          endif

          thisOctal%rhou(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, rhouPoint, nr, lcell, lnext, &
               xyzmin, xyzdel, rmax, rsq, a)

          call qshep3 (nPoints, xPoint, yPoint, zPoint, rhovPoint, nq, nw, nr, lcell, lnext, xyzmin, &
               xyzdel, rmax, rsq, a, ier )
          if (ier /= 0) then
             write(message,*) "Qshep3 returned an error ",ier
             call writeWarning(message)
          endif

          thisOctal%rhov(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, rhovPoint, nr, lcell, lnext, &
               xyzmin, xyzdel, rmax, rsq, a)

          call qshep3 (nPoints, xPoint, yPoint, zPoint, rhowPoint, nq, nw, nr, lcell, lnext, xyzmin, &
               xyzdel, rmax, rsq, a, ier )
          if (ier /= 0) then
             write(message,*) "Qshep3 returned an error ",ier
             call writeWarning(message)
          endif

          thisOctal%rhow(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, rhowPoint, nr, lcell, lnext, &
               xyzmin, xyzdel, rmax, rsq, a)

          call qshep3 (nPoints, xPoint, yPoint, zPoint, energyPoint, nq, nw, nr, lcell, lnext, xyzmin, &
               xyzdel, rmax, rsq, a, ier )
          if (ier /= 0) then
             write(message,*) "Qshep3 returned an error ",ier
             call writeWarning(message)
          endif

          thisOctal%energy(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, energyPoint, nr, lcell, lnext, &
               xyzmin, xyzdel, rmax, rsq, a)


          call qshep3 (nPoints, xPoint, yPoint, zPoint, phiPoint, nq, nw, nr, lcell, lnext, xyzmin, &
               xyzdel, rmax, rsq, a, ier )
          if (ier /= 0) then
             write(message,*) "Qshep3 returned an error ",ier
             call writeWarning(message)
          endif

          thisOctal%phi_gas(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, phiPoint, nr, lcell, lnext, &
               xyzmin, xyzdel, rmax, rsq, a)

          call qshep3 (nPoints, xPoint, yPoint, zPoint, pressurePoint, nq, nw, nr, lcell, lnext, xyzmin, &
               xyzdel, rmax, rsq, a, ier )
          if (ier /= 0) then
             write(message,*) "Qshep3 returned an error ",ier
             call writeWarning(message)
          endif

          thisOctal%pressure_i(iSubcell) = qs3val(x, y, z, nPoints, xPoint, yPoint, zPoint, pressurePoint, nr, lcell, lnext, &
               xyzmin, xyzdel, rmax, rsq, a)

       endif
          deallocate(lCell, lnext, rsq, a)

       endif
       if (thisOctal%twod) then
          rVec = subcellcentre(thisOctal, iSubcell)
          x = rVec%x
          y = 0.d0
          z = rVec%z
!          nPoints = 0
          rVec%y = 0.d0

!          radius = 4.d0*grid%octreeRoot%subcellSize / &
!                                2.0_oc**REAL(minDepthAmr,kind=oct)

          radius = thisOctal%subcellSize*4.d0

          do while (nPoints < 12)
 !            call returnAddNewChildPointArrays(thisOctal, topOctal, topOctalSubcell, rhoPoint, rhoePoint, rhouPoint, &
 !                 rhovPoint, rhowPoint, phiPoint, energyPoint, pressurePoint)
             call getPointsInRadius(rVec, radius, grid, npoints, rhoPoint, rhoePoint, &
                  rhouPoint, rhovPoint, rhowPoint, energyPoint, pressurePoint, phiPoint, xPoint, yPoint, zPoint)
             radius = radius * 2.d0
          end do

          ypoint = 0.d0
          y = 0.d0
          rhovPoint = 0.d0
                    
          nq = min(40, nPoints - 1)
          nw = min(40, nPoints - 1)
          nr = max(1,nint((dble(nPoints)/3.d0)**0.5d0))
          allocate(lCell2d(1:nr,1:nr))
          allocate(lnext(1:nPoints))
          allocate(rsq(1:nPoints))
          allocate(a(5,1:nPoints))

          call qshep2 (nPoints, xPoint, zPoint, rhoPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
               dx, dz, rmax, rsq, a, ier )
          if (ier /= 0) call writeWarning("Qshep2 returned an error for rho")
          thisOctal%rho(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhoPoint, nr, lcell2d, lnext, &
               xmin, zmin, dx, dz, rmax, rsq, a)


          call qshep2 (nPoints, xPoint, zPoint, rhoePoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
               dx, dz, rmax, rsq, a, ier )

          if (ier /= 0) call writeWarning("Qshep2 returned an error for rhoe")
          thisOctal%rhoe(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhoePoint, nr, lcell2d, lnext, &
               xmin, zmin, dx, dz, rmax, rsq, a)

          call qshep2 (nPoints, xPoint, zPoint, rhouPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
               dx, dz, rmax, rsq, a, ier )

          if (ier /= 0) call writeWarning("Qshep2 returned an error for rhou")
          thisOctal%rhou(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhouPoint, nr, lcell2d, lnext, &
               xmin, zmin, dx, dz, rmax, rsq, a)

          call qshep2 (nPoints, xPoint, zPoint, rhovPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
               dx, dz, rmax, rsq, a, ier )

          if (ier /= 0) call writeWarning("Qshep2 returned an error for rhov")
          thisOctal%rhov(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhovPoint, nr, lcell2d, lnext, &
               xmin, zmin, dx, dz, rmax, rsq, a)

          call qshep2 (nPoints, xPoint, zPoint, rhowPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
               dx, dz, rmax, rsq, a, ier )

          if (ier /= 0) call writeWarning("Qshep2 returned an error for rhow")
          thisOctal%rhow(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, rhowPoint, nr, lcell2d, lnext, &
               xmin, zmin, dx, dz, rmax, rsq, a)

          call qshep2 (nPoints, xPoint, zPoint, energyPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
               dx, dz, rmax, rsq, a, ier )

          if (ier /= 0) call writeWarning("Qshep2 returned an error for energy")
          thisOctal%energy(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, energyPoint, nr, lcell2d, lnext, &
               xmin, zmin, dx, dz, rmax, rsq, a)

          call qshep2 (nPoints, xPoint, zPoint, phiPoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
               dx, dz, rmax, rsq, a, ier )

          if (ier /= 0) call writeWarning("Qshep2 returned an error for phi")
          thisOctal%phi_gas(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, phiPoint, nr, lcell2d, lnext, &
               xmin, zmin, dx, dz, rmax, rsq, a)


          call qshep2 (nPoints, xPoint, zPoint, pressurePoint, nq, nw, nr, lcell2d, lnext, xmin, zmin, &
               dx, dz, rmax, rsq, a, ier )

          if (ier /= 0) call writeWarning("Qshep2 returned an error for pressure")
          thisOctal%pressure_i(iSubcell) = qs2val(x, z, nPoints, xPoint, zPoint, pressurePoint, nr, lcell2d, lnext, &
               xmin, zmin, dx, dz, rmax, rsq, a)


          deallocate(lCell2d, lnext, rsq, a)

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
!!    oldEnergy = parent%rhoe(parentSubcell) * cellVolume(parent, parentSubcell)

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
    factor = oldEnergy/newEnergy
    thisOctal%rhoe(1:thisOctal%maxChildren) = thisOctal%rhoe(1:thisOctal%maxChildren) * factor

!
!!!    ! momentum (u)

    oldMom = abs(parent%rhou(parentSubcell))
    newMom = SUM(abs(thisOctal%rhou(1:thisOctal%maxChildren)))/dble(thisOctal%maxChildren)
    if (abs(newMom) > TINY(newMom)) then
       factor = oldMom / newMom
       thisOctal%rhou(1:thisOctal%maxChildren) = thisOctal%rhou(1:thisOctal%maxChildren) * factor
    endif
!
!    ! momentum (v)
!
    oldMom = abs(parent%rhov(parentSubcell))
    newMom = SUM(abs(thisOctal%rhov(1:thisOctal%maxChildren)))/dble(thisOctal%maxChildren)
    if (abs(newMom) > TINY(newMOM)) then
       factor = oldMom / newMom
       thisOctal%rhov(1:thisOctal%maxChildren) = thisOctal%rhov(1:thisOctal%maxChildren) * factor
    endif
!
    ! momentum (w)

    oldMom = abs(parent%rhow(parentSubcell))
    newMom = SUM(abs(thisOctal%rhow(1:thisOctal%maxChildren)))/dble(thisOctal%maxChildren)
    if (abs(newMom) > TINY(newMom)) then
       factor = oldMom / newMom
       thisOctal%rhow(1:thisOctal%maxChildren) = thisOctal%rhow(1:thisOctal%maxChildren) * factor
    endif

!
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

  subroutine returnAddNewChildCornerArrays(thisOCtal, topOctal, topOctalSubcell, grid, rhoCorner, rhoeCorner, rhouCorner, &
       rhovCorner, rhowCorner, eCorner, phiCorner, pressureCorner)
    use inputs_mod, only : maxDepthAMR, minDepthAmr
    use mpi
    type(OCTAL), pointer :: thisOctal, topOctal
    type(GRIDTYPE) :: grid
    integer :: iCorner, iDir, nCorner, nDir
    integer :: nd, topOctalSubcell, j
    type(VECTOR) :: dir(8), corner(8), position, rVec, centre
    real(double) :: rhoCorner(8)
    real(double) :: rhoeCorner(8)
    real(double) :: rhouCorner(8)
    real(double) :: rhovCorner(8)
    real(double) :: rhowCorner(8)
    real(double) :: eCorner(8)
    real(double) :: phiCorner(8)
    real(double) :: pressureCorner(8)
    real(double) :: weight, totalWeight
    real(double) :: rho, rhoe, rhou, rhov, rhow, r, energy, phi, pressure
    real(double) :: xh, yh, zh, smallDist
    logical :: debug, addLocal

    debug=.false.

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
       
       nCorner = 8
       r = topOctal%subcellSize/2.d0

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

    rhoCorner = 0.d0
    rhoeCorner = 0.d0
    eCorner = 0.d0
    rhouCorner = 0.d0
    rhovCorner = 0.d0
    rhowCorner = 0.d0
    phicorner  = 0.d0
    pressureCorner = 0.d0

    if (debug) then
       write(*,*) "addnewchild with interp debug"
    endif

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
             rhoCorner(iCorner) = rhoCorner(iCorner) + weight * rho
             if (debug) then
                write(*,'(a,i4,i4,3f13.4,1pe13.5)') "outside ",icorner,idir,position,rho
             endif

             rhoeCorner(iCorner) = rhoeCorner(iCorner) + weight * rhoe
             rhouCorner(iCorner) = rhouCorner(iCorner) + weight * rhou
             rhovCorner(iCorner) = rhovCorner(iCorner) + weight * rhov
             rhowCorner(iCorner) = rhowCorner(iCorner) + weight * rhow
             eCorner(iCorner) = eCorner(iCorner) + weight * energy
             phiCorner(iCorner) = phiCorner(iCorner) + weight * phi
             pressureCorner(iCorner) = pressureCorner(iCorner) + weight * pressure
          else
             weight = 1.d0
             totalWeight = totalWeight + weight
             rhoCorner(iCorner) = rhoCorner(iCorner) + topOctal%rho(topOctalSubcell)
             if (debug) then
                write(*,'(a,i4,i4,3f13.4,1pe13.5)') "inside  ",icorner,idir,position,topOCtal%rho(topOctalSubcell)
             endif
             rhoeCorner(iCorner) = rhoeCorner(iCorner) + topOctal%rhoe(topOctalSubcell)
             rhouCorner(iCorner) = rhouCorner(iCorner) + topOctal%rhou(topOCtalSubcell)
             rhovCorner(iCorner) = rhovCorner(iCorner) + topOctal%rhov(topOctalSubcell)
             rhowCorner(iCorner) = rhowCorner(iCorner) + topOctal%rhow(topOctalSubcell)
             eCorner(iCorner) = eCorner(iCorner) + topOctal%energy(topOctalSubcell)
             phiCorner(iCorner) = phiCorner(iCorner) + topOctal%phi_i(topOctalSubcell)
             pressureCorner(iCorner) = pressureCorner(iCorner) + topOctal%pressure_i(topOctalSubcell)
             rVec = subcellCentre(topOctal, topOctalSubcell)
             j = j + 1
          endif
       enddo
       if (totalWeight > 0.d0) then
          rhoCorner(iCorner) = rhoCorner(iCorner) / totalWeight
          rhoeCorner(iCorner) = rhoeCorner(iCorner) / totalWeight
          rhouCorner(iCorner) = rhouCorner(iCorner) / totalWeight
          rhovCorner(iCorner) = rhovCorner(iCorner) / totalWeight
          rhowCorner(iCorner) = rhowCorner(iCorner) / totalWeight
          eCorner(iCorner) = eCorner(iCorner) / totalWeight
          phiCorner(iCorner) = phiCorner(iCorner) / totalWeight
          pressureCorner(iCorner) = pressureCorner(iCorner) / totalWeight
       else
          rhoCorner(iCorner) = topOctal%rho(topOctalSubcell)
          rhoeCorner(iCorner) = topOctal%rhoe(topOctalSubcell)
          rhouCorner(iCorner) = topOctal%rhou(topOCtalSubcell)
          rhovCorner(iCorner) = topOctal%rhov(topOctalSubcell)
          rhowCorner(iCorner) = topOctal%rhow(topOctalSubcell)
          eCorner(iCorner) = topOctal%energy(topOctalSubcell)
          phiCorner(iCorner) = topOctal%phi_i(topOctalSubcell)
          pressureCorner(iCorner) = topOctal%pressure_i(topOctalSubcell)
       endif
    enddo

  end subroutine returnAddNewChildCornerArrays

  subroutine returnAddNewChildPointArrays(thisOctal, topOctal, topOctalSubcell,  grid, rhoPoint, rhoePoint, &
       rhouPoint, rhovPoint, rhowPoint, phiPoint, energyPoint, pressurePoint, npoints)
    use inputs_mod, only : maxDepthAMR, minDepthAmr
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

!    do i = 1, nPoints
!       do counter = 1, nPoints
!          if( i /= counter) then
!             if(x(i) == x(counter) .and. y(i) == y(counter) .and. z(i) == z(counter)) then
!                print *, "x", x(i), x(counter)
!                print *, "y", y(i), y(counter)
!                print *, "z", z(i), z(counter)
!                call torus_abort("DUPLICATE ENTRY B")
!             end if
!          end if
!       end do
!    end do


   
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

    if(grid%octreeRoot%twoD) then
       if((nHydrothreadsGlobal) == 4) then
          evenUpArray = (/1, 2, 3, 4/)
       else if((nHydrothreadsGlobal) == 16) then
          evenUpArray = (/1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4/)
       end if
    else if(grid%octreeRoot%threeD) then
       if((nHydrothreadsGlobal) == 8) then
          evenUpArray = (/1, 2, 3, 4, 5, 6, 7, 8/)
       else if((nHydrothreadsGlobal) == 64) then
          evenUpArray = (/1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 1
               1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 2           
               1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 3           
               1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8/)!pane 4             

       else if(nHydroTHreadsGlobal==512) then
            evenUpArray = (/1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, &
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 1  
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, &
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 2           
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, &
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 3           
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, &
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 4   
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, &
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 5   
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, &
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 6   
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, &
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, & !pane 7   
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8, &
                 1, 3, 2, 4, 5, 7, 6, 8, 1, 3, 2, 4, 5, 7, 6, 8/) !pane 8    
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
!       print *, "RANK ", myRankGlobal, "SENDING TO ", iThread
       call MPI_SEND(loc, 7, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
 !      print *, "RANK ", myRankGlobal, "SENT"
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
    integer :: ierr, j, k, m, counter, counter2
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

       if(loc(7) == 1.d0) then
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
          cellRef=1
          !Get the main serving thread's values within the request radius
          call getAllInRadius(grid%octreeRoot, position, searchRadius, storageArray, cellRef, useTop)

          !send the order to the other serving threads to return values
          nvals = 0
          do m = 1, nHydroThreadsGlobal
             if (m /= myrankGlobal .and. .not. ANY(m == check(k,1:nHydroThreadsGlobal))) then 
                !Will recv a 1D array and translate it into 2D later
!                call MPI_RECV(tempStorageArray 12000, MPI_DOUBLE_PRECISION, m, tag, localWorldCommunicator, status, ierr)
                call MPI_RECV(tempDoubleArray, SIZE(tempDoubleArray), MPI_DOUBLE_PRECISION, m, tag, &
                     localWorldCommunicator, status, ierr)

                !Now convert it back to 2D
                tempStorageArray = reshape(tempDoubleArray, shape(tempStorageArray))

                do counter = 1, maxStorage
                   if(tempstorageArray(counter,1) > 1.d0) then
                      do counter2 = 1, maxStorage+1
                         if(counter2 == maxStorage+1) then
                            call torus_abort("Fetched too many variables")
                         end if
                         if(storageArray(counter2, 1) == 0.d0) then
                            nVals = counter2
                            storageArray(counter2,:) = tempStorageArray(counter,:)
                            exit

                         end if                         
                      end do
                    end if
                end do
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
          cellRef = 1
          storageArray = 0.d0
          call getAllInRadius(grid%octreeRoot, position, searchRadius, storageArray, cellRef, useTop)
 
          !Don't want to send a 2D array
          tempDoubleArray = reshape(storageArray,(/SIZE(tempDoubleArray)/))
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
                      cellRef = cellRef + 1
                   end if
                else
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
                   cellRef = cellRef + 1
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
!              print *, "thisOctal%mpiThread(subcell) ", thisOctal%mpiThread(subcell)
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

end subroutine labelSingleSubcellMPI


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
      integer :: i
      i = thisOctal%nchildren
      i = subcell
      i = myrank
      ! Set return value of function to prevent a compiler warning. 
      check = .true.
!      check=.false.
!      call writefatal("octalOnThread called in non-MPI code")
!      STOP

    end function octalOnThread

#endif
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

  end module mpi_amr_mod
