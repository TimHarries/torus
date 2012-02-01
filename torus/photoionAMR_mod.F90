#ifdef PHOTOION
!Photoionization module - started on October 4th 2005 by th 

module photoionAMR_mod

#ifdef MPI
use constants_mod
use messages_mod

use parallel_mod
use photoion_utils_mod
use gridio_mod
use source_mod
use timing
use image_mod
use grid_mod
use amr_mod
use diffusion_mod
use mpi_amr_mod
use mpi_global_mod
use unix_mod, only: unixGetenv
use photon_mod
use phasematrix_mod
use phfit_mod, only : phfit2
use vtk_mod

implicit none

private
public :: photoIonizationloopAMR, createImagesplitgrid, ionizeGrid, &
     neutralGrid, resizePhotoionCoeff, resetNH, hasPhotoionAllocations, allocatePhotoionAttributes
#ifdef HYDRO
public :: radiationHydro
#endif

type PHOTONPACKET
    type(VECTOR) :: rVec
    type(VECTOR) :: uHat
    real(double) :: Freq
    real(double) :: tPhot
    real(double) :: ppw
    integer :: destination
    logical :: sourcePhoton
    logical :: crossedPeriodic
end type PHOTONPACKET

  real :: heIRecombinationFit(32, 3, 3)
  real :: heIRecombinationLambda(32)
  real :: heIRecombinationNe(3)
  real :: heIIrecombinationLines(3:30, 2:16)
  type(GAMMATABLE) :: gammaTableArray(3) ! H, HeI, HeII

contains

#ifdef HYDRO
  subroutine radiationHydro(grid, source, nSource, nLambda, lamArray)
    use inputs_mod, only : iDump, doselfgrav, readGrid, maxPhotoIonIter, tdump, tend, justDump !, hOnly
    use inputs_mod, only : dirichlet, amrtolerance, nbodyPhysics, amrUnrefineTolerance, smallestCellSize, dounrefine
    use dust_mod, only : emptyDustCavity
    use hydrodynamics_mod, only: hydroStep3d, calculaterhou, calculaterhov, calculaterhow, &
         calculaterhoe, setupedges, unsetGhosts, setupghostcells, evenupgridmpi, refinegridgeneric, &
         setupx, setupqx, computecouranttime, unrefinecells, selfgrav, sumgasstargravity, transfertempstorage, &
         zerophigas, zerosourcepotential, applysourcepotential, addStellarWind, cutVacuum, setupEvenUpArray, &
         pressureGradientTimestep, mergeSinks, addSinks, ComputeCourantTimenBody, &
         perturbIfront
    use dimensionality_mod, only: setCodeUnit
    use inputs_mod, only: timeUnit, massUnit, lengthUnit, readLucy, checkForPhoto, severeDamping
    use inputs_mod, only: singleMegaPhoto
    use parallel_mod, only: torus_abort
    use mpi
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nLambda
    real :: lamArray(:)
    character(len=80) :: mpiFilename, datFilename
    real(double) :: dt, tc(513), temptc(513),cfl, gamma, mu
    integer :: iUnrefine
    integer :: ierr
    real(double) ::  nextDumpTime
    type(VECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: thread1(512), thread2(512), nBound(512), nPairs
    integer :: group(1000), nGroup
!    logical :: globalConverged(512), tConverged(512)
    integer :: nHydroThreads
    logical :: dumpThisTime
    real(double) :: deltaTforDump, timeOfNextDump, loopLimitTime
    integer :: iRefine, nUnrefine
    logical :: startFromNeutral
    logical :: photoLoop, photoLoopGlobal=.false.
    integer :: i, status, tag=30, sign
    integer :: stageCounter=1,  nPhase, nstep
    real(double) :: timeSinceLastRecomb=0.d0
    logical :: noPhoto=.false.
    integer :: evenUpArray(nThreadsGlobal-1)
    real :: iterTime(3)
    integer :: iterStack(3)
    integer :: optID
    nHydroThreads = nThreadsGlobal-1
    dumpThisTime = .false.



    if (nbodyPhysics) then
       if (writeoutput) then
          open(57, file="pos.dat", status="unknown", form="formatted")
          close(57)
       endif
    endif

    direction = VECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    cfl = 0.3d0
    
    mu = 1.d0

    sign = 1
    viewVec = VECTOR(-1.d0,0.d0,0.d0)
    viewVec = rotateZ(viewVec, 40.d0*degtorad)
    viewVec = rotateY(viewVec, 30.d0*degtorad)

    if (writeoutput) write(*,*) "Ionizing photons per second ",ionizingFlux(source(1))

    call setCodeUnit(time=timeUnit)
    call setCodeUnit(mass=massUnit)
    call setCodeUnit(length=lengthUnit)
    call writeCodeUnits()

!    tdump = returnCodeUnitTime(tdump)




!    deltaTforDump = 3.14d10 !1kyr
!    if (grid%geometry == "hii_test") deltaTforDump = 2.d10



    if(.not. readGrid) then
       grid%currentTime = 0.d0
       grid%iDump = 0
       deltaTforDump = 3.14d10 !1kyr

       deltaTforDump = tdump

       if (grid%geometry == "hii_test") deltaTforDump = 9.d10
       if(grid%geometry == "bonnor") deltaTforDump = (1.57d11)!/5.d0 !5kyr
       if(grid%geometry == "radcloud") deltaTforDump = (1.57d11)!/5.d0 !5kyr
       if(grid%geometry == "starburst") deltaTforDump = tdump
       if(grid%geometry == "molefil") deltaTforDump = tdump
       if(grid%geometry == "sphere") deltaTforDump = tdump
       if(grid%geometry == "turbulence") then
!turbulence phase
          deltaTforDump = 1.57d12 !50kyr
          tend = 3.14d13        !1Myr
          noPhoto = .true.
       end if
       
       nextDumpTime = deltaTforDump
       timeofNextDump = 0.d0       
     else
       deltaTforDump = 1.57d11 !5kyr
       if(grid%geometry == "turbulence") then
!reset clocks
          grid%currentTime = 0.d0
          deltaTforDump = 1.57d11
          tend = 40.d0*deltaTforDump
       endif
       deltaTforDump = tDump
       
!       grid%currentTime = 0.d0
!       deltaTforDump = 3.14d11
!       tend = 50.d0*deltaTforDump

       if (grid%geometry(1:6) == "sphere") then
          call emptyDustCavity(grid%octreeRoot, VECTOR(0.d0, 0.d0, 0.d0), 1400.d0*autocm/1.d10)
       endif

       nextDumpTime = grid%currentTime + deltaTforDump
       timeofNextDump = nextDumpTime
       if(justDump) then 
          write(mpiFilename,'(a, i4.4, a)') "quickDump.vtk"
          call writeVtkFile(grid, mpiFilename, &
               valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
               "hydrovelocity","sourceCont   ","pressure     ","radmom       ", "rhou         " &
               , "rhov         ", "rhow         "/))
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call torus_abort("vtk dump completed. Aborting...")
       end if
       
!turn everything on and run a long calculation
       if(singleMegaPhoto) then
          maxPhotoionIter = 200
       end if
    end if

!    call writeVtkFile(grid, "ini.vtk", &
!         valueTypeString=(/"rho        ","HI         " ,"temperature", "sourceCont ", "mpithread  " /))

    iunrefine = 0
    startFromNeutral = .false.

    if (readlucy) then
       write(mpiFilename,'(a, i4.4, a)') "dump_", iDump, ".grid"
       write(*,*) "Reading from: ",trim(mpiFilename)
       call readAmrGrid(mpiFilename,.false.,grid)
    endif

    if (myRankWorldGlobal == 1) write(*,*) "CFL set to ", cfl

    if (myrankGlobal /= 0) then

       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       call writeInfo("Setting up even up array", TRIVIAL)
       call setupEvenUpArray(grid, evenUpArray)
       call writeInfo("Done", TRIVIAL)
       call resetnh(grid%octreeRoot)


!       do i = 1, nPairs
!          if (myrankWorldglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
!       enddo

       if(grid%currentTime == 0.d0) then
          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call calculateRhoU(grid%octreeRoot, direction)
          direction = VECTOR(0.d0, 1.d0, 0.d0)
          call calculateRhoV(grid%octreeRoot, direction)
          direction = VECTOR(0.d0, 0.d0, 1.d0)
          call calculateRhoW(grid%octreeRoot, direction)
          
          call calculateEnergyFromTemperature(grid%octreeRoot)
          
          call calculateRhoE(grid%octreeRoot, direction)


          call writeInfo("Refining individual subgrids", TRIVIAL)
          call setAllUnchanged(grid%octreeRoot)
          if (.not.grid%splitOverMpi) then
             do
                gridConverged = .true.
                call setupEdges(grid%octreeRoot, grid)
                !             call refineEdges(grid%octreeRoot, grid,  gridconverged, inherit=.false.)
                call unsetGhosts(grid%octreeRoot)
                call setupGhostCells(grid%octreeRoot, grid)
                if (gridConverged) exit
             end do
          else
             call evenUpGridMPI(grid, .false., .true., evenuparray)
          endif
          call resetnh(grid%octreeRoot)
          
          call writeInfo("Calling exchange across boundary", TRIVIAL)
          if (myRankWorldGlobal == 1) call tune(6,"Exchange across boundary")
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          if (myRankWorldGlobal == 1) call tune(6,"Exchange across boundary")
          call writeInfo("Done", TRIVIAL)
          
!          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          
          call evenUpGridMPI(grid,.false.,.true., evenuparray)      
          call refineGridGeneric(grid, amrtolerance, evenuparray)
          call writeInfo("Evening up grid", TRIVIAL)    
          call evenUpGridMPI(grid, .false.,.true., evenuparray)
          call resetnh(grid%octreeRoot)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call setupX(grid%octreeRoot, grid, direction)
          call setupQX(grid%octreeRoot, grid, direction)
          !    call calculateEnergy(grid%octreeRoot, gam

       endif
    end if

    if(grid%currentTime == 0.d0 .and. .not. readGrid .or. singleMegaPhoto) then
       call ionizeGrid(grid%octreeRoot)
       if (grid%geometry(1:6) == "sphere") &
            call emptyDustCavity(grid%octreeRoot, VECTOR(0.d0, 0.d0, 0.d0), 1400.d0*autocm/1.d10)

       if(.not. noPhoto) then

          looplimitTime = deltaTForDump
          looplimittime = 1.d30
          iterTime = 1.e30
          do irefine = 1, 1
             if (irefine == 1) then
                call writeInfo("Calling photoionization loop",TRIVIAL)
                call setupNeighbourPointers(grid, grid%octreeRoot)
                call resetnh(grid%octreeRoot)
                call photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, maxPhotoionIter, loopLimitTime, &
                     looplimittime, .false.,iterTime,.true., evenuparray, optID, iterStack)
                call writeInfo("Done",TRIVIAL)
             else
                call writeInfo("Calling photoionization loop",TRIVIAL)
                call setupNeighbourPointers(grid, grid%octreeRoot)
                call photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, maxPhotoionIter, loopLimitTime, &
                     looplimittime, .false.,iterTime,.true., evenuparray, optID, iterStack)
                call writeInfo("Done",TRIVIAL)
             endif
             
             call writeInfo("Dumping post-photoionization data", TRIVIAL)
             call writeVtkFile(grid, "start.vtk", &
             valueTypeString=(/"rho        ","HI         " ,"temperature", "sourceCont " /))
          end do
       end if

       if(singleMegaPhoto) then
          call torus_abort("Finished single photo calculation, ending...")
       end if

!       if (myrank /= 0) call addStellarWind(grid%octreeRoot, globalsourcearray(1))
          
          if (myrankGlobal /= 0) then
             call calculateEnergyFromTemperature(grid%octreeRoot)
             call calculateRhoE(grid%octreeRoot, direction)
          endif
                    
          if (myrankGlobal/=0) then
             call writeInfo("Refining individual subgrids", TRIVIAL)
             call setAllUnchanged(grid%octreeRoot)
             if (.not.grid%splitOverMpi) then
                do
                   gridConverged = .true.
                   call setupEdges(grid%octreeRoot, grid)
!             call refineEdges(grid%octreeRoot, grid,  gridconverged, inherit=.false.)
                   call unsetGhosts(grid%octreeRoot)
                   call setupGhostCells(grid%octreeRoot, grid)
                   if (gridConverged) exit
                end do
             else
                call evenUpGridMPI(grid, .false., .true., evenuparray)
             endif
             
             call evenUpGridMPI(grid,.false.,.true., evenuparray)      
             call refineGridGeneric(grid, amrtolerance, evenuparray)
             call writeInfo("Evening up grid", TRIVIAL)    
             call evenUpGridMPI(grid, .false.,.true., evenuparray)
             if (myrankGlobal /= 0) call addStellarWind(grid, globalnSource, globalsourcearray)
             call resetnh(grid%octreeRoot)
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!                call writeVtkFile(grid, "stellarwind.vtk", &
!                valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
!                ,"phi          " /))


             !THAW - self gravity

             if (doselfGrav) then
                if (myrankWorldGlobal == 1) call tune(6, "Self-Gravity")
                if (myrankWorldGlobal == 1) write(*,*) "Doing multigrid self gravity"
!                call writeVtkFile(grid, "beforeselfgrav.vtk", &
!                valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", "phigas       " &
!                ,"phi          " /))
                
                call zeroPhiGas(grid%octreeRoot)
                call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.) 
                
                if(.not. dirichlet) then
                   call periodBoundary(grid, justGrav = .true.)
                   call transferTempStorage(grid%octreeRoot, justGrav = .true.)
                end if
                                
                call zeroSourcepotential(grid%octreeRoot)
                if (globalnSource > 0) then
                   call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, grid%halfSmallestSubcell)
                endif
                call sumGasStarGravity(grid%octreeRoot)
                
!                call writeVtkFile(grid, "afterselfgrav.vtk", &
!                valueTypeString=(/"rho          ","hydrovelocity","rhoe         " ,"u_i          ", &
!                "phigas       ","phi          "/))
                if (myrankWorldGlobal == 1) write(*,*) "Done"
                if (myrankWorldGlobal == 1) call tune(6, "Self-Gravity")
              
             endif
             
          endif
          !       enddo
       
       if (myrankGlobal /= 0) then
          call calculateEnergyFromTemperature(grid%octreeRoot)
          call calculateRhoE(grid%octreeRoot, direction)          
       endif
    end if
    
    if(grid%geometry == "planar") then
       print *, "PERTURBING I FRONT "
       call perturbIfront(grid%octreeRoot, grid)
    end if
      
      if(grid%geometry == "hii_test") then
         tEnd = 3.14d13         !1x10^6 year
      else if(grid%geometry == "bonnor" .or. grid%geometry=="radcloud") then
         tEnd = 200.d0*3.14d10 !200kyr 
      end if

    nPhase = 1

    nstep = 0
    
!    !Thaw - trace courant time history                                                                                                                              
!    open (444, file="tcHistory.dat", status="unknown")

    do while(grid%currentTime < tEnd)
       nstep = nstep + 1

       
       tc = 0.d0
       if (myrankGlobal /= 0) then
          tc(myrankGlobal) = 1.d30
          call computeCourantTime(grid, grid%octreeRoot, tc(myRankGlobal))
          if (nbodyPhysics) call computeCourantTimeNbody(grid, globalnSource, globalsourceArray, tc(myrankGlobal))
          call pressureGradientTimeStep(grid, dt)
          tc(myRankGlobal) = min(tc(myrankGlobal), dt)
       endif
       call MPI_ALLREDUCE(tc, tempTc, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)
       dt = MINVAL(temptc(1:nHydroThreadsGlobal))
       dt = dt * dble(cfl)


!       write(444, *) jt, MINVAL(tc(1:nHydroThreads)), dt
       
!       if (nstep < 3) then
!          dt = dt * 0.01d0
!       endif
  
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       if (myrankWorldGlobal == 1) then
!          write(*,*) myHydroSetGlobal,temptc(1:nHydroThreadsGlobal)
          write(*,*) "courantTime", dt
       endif

!       if (myrankGlobal /= 0) call checkSetsAreTheSame(grid%octreeRoot)

       dumpThisTime = .false.

       if ((grid%currentTime + dt) >= timeOfNextDump) then
          dt =  timeofNextDump - grid%currentTime
          dumpThisTime = .true.
       endif

!       write(mpiFilename,'(a, i4.4, a)') "preStep.vtk"
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!            "hydrovelocity","sourceCont   ","pressure     ", "rhou         "/))



       if (myrankWorldGlobal == 1) write(*,*) "dump ",dumpThisTime, " current ", &
            grid%currentTime, " deltaTfordump ",deltaTforDump, " dt ", dt

       if (myrankWorldGlobal == 1) write(*,*) "Time step", dt
       if (myRankWorldGlobal == 1) write(*,*) "percent to next dump ",100.*(timeofNextDump-grid%currentTime)/deltaTforDump

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       if(checkForPhoto) then
          photoLoopGlobal = .false.
          photoLoop = .false.
          if(myRankGlobal /= 0) then
            ! call checkForPhotoLoop(grid, grid%octreeRoot, photoLoop, dt)
             call advancedCheckForPhotoLoop(grid, grid%octreeRoot, photoLoop, dt, timeSinceLastRecomb)
             call MPI_SEND(photoLoop, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator,  ierr)
             call MPI_RECV(photoLoopGlobal, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator, status, ierr)
          else
             do i = 1, nThreadsGlobal-1
                photoLoop = .false.
                call MPI_RECV(photoLoop, 1, MPI_LOGICAL, i, tag, localWorldCommunicator, status, ierr)
              
                if (photoLoop) then
                   photoLoopGlobal = .true.
                else if(.not. photoLoopGlobal) then
                   photoLoopGlobal = .false.
                end if
             end do
            
             do i = 1, nThreadsGlobal -1
                call MPI_SEND(photoLoopGlobal, 1, MPI_LOGICAL, i, tag, localWorldCommunicator,  ierr)
             end do
          end if
       
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          
       else
          photoLoopGlobal = .true.
       end if

       if(photoLoopGlobal .and. .not. noPhoto) then
          call writeInfo("Calling photoionization loop",TRIVIAL)
!          call ionizeGrid(grid%octreeRoot)
          if(dt /= 0.d0) then
             loopLimitTime = grid%currentTime+dt
          else
             looplimittime = deltaTForDump
          end if
          looplimittime = 1.d30
          call setupNeighbourPointers(grid, grid%octreeRoot)
          call resetnh(grid%octreeRoot)
          if (grid%geometry == "sphere") &
               call emptyDustCavity(grid%octreeRoot, VECTOR(0.d0, 0.d0, 0.d0), 1400.d0*autocm/1.d10)

!          call photoIonizationloopAMR(grid, source, nSource, nLambda,lamArray, 1, loopLimitTime, loopLimitTime, .false., iterTime, &
!               .true., evenuparray, sign)
          call photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, 1, loopLimitTime, &
               looplimittime, .false.,iterTime,.true., evenuparray, optID, iterStack)
          call writeInfo("Done",TRIVIAL)
          timeSinceLastRecomb = 0.d0
       else
          timeSinceLastRecomb = timeSinceLastRecomb + dt
         call writeInfo("Skipping photoionoization loop",TRIVIAL)
!          if(quickThermal) then
!             nTimes = 3
!             do i = 1, nTimes
!                call quickThermalCalc(grid%octreeRoot)
!             enddo
!          end if
       end if

          if (myrankGlobal /= 0) then
             call calculateEnergyFromTemperature(grid%octreeRoot)
             call calculateRhoE(grid%octreeRoot, direction)
          endif


       if (myRankWorldGlobal == 1) call tune(6,"Hydrodynamics step")
       call writeInfo("calling hydro step",TRIVIAL)
       
       if (myrankGlobal /= 0) call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

!       write(mpiFilename,'(a, i4.4, a)') "preStep.vtk"
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!            "hydrovelocity","sourceCont   ","pressure     "/))

       if (severeDamping) call cutVacuum(grid%octreeRoot)

!       write(mpiFilename,'(a, i4.4, a)') "preStep.vtk"
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!            "hydrovelocity","sourceCont   ","pressure     ", &
!            "mpithread    "/))
       if (myrankGlobal /= 0) call hydroStep3d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup,doSelfGrav=doselfGrav)

          if (myRankGlobal /= 0) then
             
             call setAllUnchanged(grid%octreeRoot)
             call evenUpGridMPI(grid,.false.,.true., evenuparray)

             call refineGridGeneric(grid, amrtolerance, evenuparray)

             call writeInfo("Evening up grid", TRIVIAL)
             call evenUpGridMPI(grid, .false.,.true., evenuparray)
             if (myrankGlobal /= 0) call addStellarWind(grid, globalnSource, globalsourcearray)
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

          endif

          if (severeDamping) call cutVacuum(grid%octreeRoot)



!add/merge sink particles where necessary
       if (myrankGlobal /= 0) then
          if (nbodyPhysics) call addSinks(grid, globalsourceArray, globalnSource)       
          if (nbodyPhysics) call mergeSinks(grid, globalsourceArray, globalnSource)
          do i = 1, globalnSource
             call emptySurface(globalsourceArray(i)%surface)
             call buildSphereNbody(globalsourceArray(i)%position, 2.5d0*smallestCellSize, &
                  globalsourceArray(i)%surface, 20)
          enddo

       endif


!       write(mpiFilename,'(a, i4.4, a)') "postStep.vtk"
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!            "hydrovelocity","sourceCont   ","pressure     "/))

       if (myRankWorldGlobal == 1) call tune(6,"Hydrodynamics step")
       if (myrankGlobal /= 0) call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       if (myrankGlobal /= 0) call resetNh(grid%octreeRoot)

       if (myrankGlobal /= 0) then
          if (dounrefine) then
             iUnrefine = iUnrefine + 1
             if (iUnrefine == 3) then
                if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
                call unrefineCells(grid%octreeRoot, grid, nUnrefine,amrUnrefinetolerance)
                if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
                iUnrefine = 0
             endif
             call evenUpGridMPI(grid, .true., .true., evenuparray)
             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          endif
       endif

       grid%currentTime = grid%currentTime + dt
       
       if (myRankWorldGlobal == 1) write(*,*) "Current time: ",grid%currentTime

!Track the evolution of the ionization front with time
!          write(datFilename, '(a, i4.4, a)') "Ifront.dat"
!          call dumpStromgrenRadius(grid, datFileName, VECTOR(-3.3d8,  -3.3d8, 3.3d8), &
!               VECTOR(3.3d8, 3.3d8, 3.3d8), 1000)


!       write(mpiFilename,'(a, i4.4, a)') "dump_", nPhase,".vtk"
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!            "hydrovelocity","sourceCont   ","pressure     "/))
!       nPhase = nPhase + 1
       if (dumpThisTime) then

          !Thaw to match time dumps of other codes
          if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d5)) then
             deltaTForDump = 1.d5
          end if
          if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d6)) then
             deltaTForDump = 1.d6
          end if
          if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d7)) then
             deltaTForDump = 1.d7
          end if
          if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d8)) then
             deltaTForDump = 1.d8
          end if
          if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d9)) then
             deltaTForDump = 1.d9
          end if
          if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d10)) then
             deltaTForDump = 1.d10
          end if
          if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d11)) then
             deltaTForDump = 1.d11
          end if

          timeOfNextDump = timeOfNextDump + deltaTForDump
          grid%iDump = grid%iDump + 1

          write(mpiFilename,'(a, i4.4, a)') "dump_", grid%iDump,".grid"
          call writeAmrGrid(mpiFilename, .false., grid)
          write(mpiFilename,'(a, i4.4, a)') "dump_", grid%iDump,".vtk"
          call writeVtkFile(grid, mpiFilename, &
               valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
               "hydrovelocity","sourceCont   ","pressure     ","radmom       "/))

          write(mpiFilename,'(a,i4.4,a)') "nbody",grid%iDump,".vtk"
          call writeVtkFilenBody(globalnSource, globalsourceArray, mpiFilename)

          if (writeoutput) then
             write(mpiFilename,'(a,i4.4,a)') "source",grid%idump,".dat"
             call writeSourceArray(mpifilename)
          endif

          
!Track the evolution of the ionization front with time
       if(grid%geometry == "hii_test") then
          write(datFilename, '(a, i4.4, a)') "hii_test",grid%iDump,".dat"
          call dumpValuesAlongLine(grid, datFileName, VECTOR(-1.5d9,  -1.5d9, 1.5d9), &
               VECTOR(1.5d9, 1.5d9, 1.5d9), 1000)

          write(datFilename, '(a, i4.4, a)') "Ifront.dat"     
          call dumpStromgrenRadius(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
               VECTOR(0.0d0, 0.0d0, 1.5d9), 1000)
       end if
 
    endif
    stageCounter = stageCounter + 1
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    close(444)
!    write(*,*) "myRank", myRankGlobal, "finishing loop. Time:", grid%currentTime, "tend ", tend
 enddo
end subroutine radiationHydro
#endif

  subroutine photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, maxIter, tLimit, deltaTime, timeDep, iterTime, &
       monteCheck, evenuparray, optID, iterStack, sublimate)
    use inputs_mod, only : quickThermal, inputnMonte, noDiffuseField, minDepthAMR, maxDepthAMR, binPhotons,monochromatic, &
         readGrid, dustOnly, minCrossings, bufferCap, doPhotorefine, hydrodynamics, doRefine, amrtolerance, hOnly, &
         optimizeStack, stackLimit, dStack, singleMegaPhoto, stackLimitArray, customStacks, tMinGlobal
    use hydrodynamics_mod, only: refinegridgeneric, evenupgridmpi
    use mpi
    implicit none
    integer(bigint) :: nMonte, nThreadMonte
    logical, optional :: sublimate
    logical :: doSublimate, timeDep, monteCheck
    real(double) :: tLimit
    real(double) :: deltaTime
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
!    character(len=80) :: message
    integer :: nlambda
    real :: lamArray(:)
    integer :: nSource
    type(SOURCETYPE) :: source(:), thisSource
    integer :: iSource
    type(VECTOR) :: rVec, uHat, rHat
    real(double) :: lCore
    integer(bigint) :: iMonte
    integer :: subcell
    integer :: i, j
    logical :: escaped
    real(double) :: wavelength, thisFreq
    real :: thisLam
    type(VECTOR) :: octVec
    real(double) :: r
    integer :: ilam
    integer(bigint) :: nInf
    integer(bigint) :: nEscapedGlobal=0
    real(double) :: kappaScadb, kappaAbsdb
    real(double) :: epsOverDeltaT
    integer :: nIter, nPeriodic
    logical :: converged, thisThreadConverged, failed
    real(double) :: photonPacketWeight
    real(double) :: tPhoton
    real(double) :: albedo
    integer :: maxIter
    type(SAHAMILNETABLE),save :: hTable, heTable
    type(RECOMBTABLE),save :: Hrecombtable
    real(double) :: freq(1006), dfreq(1006), spectrum(1006)
!    real(double) :: freq(1000+(2*grid%nIon)), dfreq(1000+(2*grid%nIon)), spectrum(1000+(2*grid%nIon)) 
    real(double) :: nuStart, nuEnd
    real(double) :: r1, kappaAbsGas, kappaAbsDust, escat
    integer, parameter :: nTemp = 9
    real(double) :: v, dustHeating
    real :: kappaP
!    integer, parameter :: nFreq = 1000
    integer :: nFreq

    logical, save :: firsttime = .true.
    !$OMP THREADPRIVATE(firstTime)

    integer(bigint) :: iMonte_beg, iMonte_end, nSCat
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: np, iOctal, iOctal_beg, iOctal_end, nOctal
    integer :: iThread
    logical :: endLoop
    logical :: crossedMPIboundary
    integer :: tag = 41
    integer(bigint) :: nTotScat, nPhot
    integer :: newThread, iSignal
    logical, save :: firstTimeTables = .true.
    integer(bigint) :: nEscaped=0
    integer, parameter :: nTimes = 3
    logical :: photonsStillProcessing
    integer(bigint), allocatable :: nEscapedArray(:)
    integer :: status(MPI_STATUS_SIZE)

    real(double) :: deltaT, fluctuationCheck 
    logical :: anyUndersampled, undersampledTOT
    logical :: underSamFailed, escapeCheck
    logical :: sourcePhoton

    integer :: dprCounter

    integer :: maxStackLimit
    integer :: sendStackLimit

    integer :: evenuparray(nHydroThreadsGlobal)
    integer :: k
    real(double) :: nuThresh

    !optimization variables
!    integer, parameter :: stackLimit=200
    real, intent(inout) :: iterTime(3)
    integer, intent(inout) :: iterStack(3)
    integer :: ZerothstackLimit, OldStackLimit
    real :: startTime, endTime, newTime
    integer, intent(inout) :: optID
    integer :: optCounter, optCount
    integer :: dstackNaught
    logical, save :: optConverged=.false.
    logical, save :: iniIterTime=.false.
    logical, save :: iniCustomTime=.true.
    real(double) :: m1, m2
  !  real :: oldTime = 1.e10
  !  integer :: newStackLimit= 0, oldStackLimit= 0

    integer :: thisPacket, sendCounter
    integer, allocatable :: nSaved(:)
    integer :: stackSize, p
    integer :: mpi_vector, mpi_photon_stack
    logical :: sendAllPhotons = .false., donePanicking = .true.
    real(double) :: nIonizingPhotons
    logical :: crossedPeriodic = .false.

!    type(PHOTONPACKET) :: photonPacketStack(stackLimit*nThreadsGlobal)
    type(PHOTONPACKET), allocatable :: photonPacketStack(:)
!    type(PHOTONPACKET) :: toSendStack(stackLimit), currentStack(stackLimit)
    type(PHOTONPACKET), allocatable :: toSendStack(:), currentStack(:)

   !Custom MPI type variables
    integer(MPI_ADDRESS_KIND) :: displacement(8)
    integer :: count = 8
    integer :: blocklengths(8) = (/ 1, 1, 1, 1, 1, 1, 1, 1/)
    integer :: oldTypes(8) 
    integer :: iDisp

    !Buffer send variables 
    integer :: bufferSize  !Send buffer size in bytes
    character, allocatable :: buffer(:)
    integer :: request
!    logical :: checkRequest
!    logical :: ready = .true.

    !OMP variables
    logical :: finished = .false.
    logical :: voidThread

    !Energy cons variables
    real(double) :: totalPower = 0.d0
    real(double) :: lams(1000) = 0.d0
    real(double) :: countArray(1000) = 0.d0
    real(double) :: tempDouble

    real, allocatable :: tempCell(:,:), temp1(:), temp2(:)

    integer :: n_rmdr, mOctal

    character(len=80) :: mpiFilename
    integer :: ierr

    !AMR
!    integer :: iUnrefine, nUnrefine


    !!Thaw - optimize stack will be run prior to a big job to ensure that the most efficient stack size is used
    !start with stack size of 1
    !if(optimizeStack) then 
    !   stackLimit = 1
    !   zerothStackLimit = 1
    !end if
    
    
    if(stacklimit == 0) then
       stacklimit = 200
    end if

    if(.not. iniIterTime) then
       iterTime = 1.e30
       optID = 1
       iniIterTime = .true.
    end if

    if(customStacks) then
       optimizeStack = .false.
       maxStackLimit = maxval(stackLimitArray(1:nHydroThreadsGlobal))
       sendStackLimit = stackLimitArray(myRankGlobal+1)
    else
       maxStackLimit = stackLimit
       sendStackLimit = stackLimit
    end if

    zerothstacklimit = stacklimit
    nPeriodic = 0
    bufferSize = 0 
    dStackNaught = dstack

!    stackLimit = 0
!    iUnrefine = 0

    
    allocate(photonPacketStack(maxStackLimit*nHydroThreadsGlobal))
    allocate(toSendStack(maxStackLimit))
    allocate(currentStack(maxStackLimit))

    !Custom MPI data types for easier send/receiving
    !MPI datatype for out TYPE(VECTOR) variables
    call MPI_TYPE_CONTIGUOUS(3, MPI_DOUBLE_PRECISION, MPI_VECTOR, ierr)
    call MPI_TYPE_COMMIT(MPI_VECTOR, ierr)

    !MPI datatype for the photon_stack data type
    oldTypes = (/ MPI_VECTOR, MPI_VECTOR, MPI_DOUBLE_PRECISION, &
         MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_LOGICAL/)

    call MPI_GET_ADDRESS(toSendStack(1)%rVec, displacement(1), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%uHat, displacement(2), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%freq, displacement(3), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%tPhot, displacement(4), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%ppw, displacement(5), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%destination, displacement(6), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%sourcePhoton, displacement(7), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%crossedPeriodic, displacement(8), ierr)

    do iDisp = 8, 1, -1
       displacement(iDisp) = displacement(iDisp) - displacement(1)
    end do

    call MPI_TYPE_CREATE_STRUCT(count, blockLengths, displacement, oldTypes, MPI_PHOTON_STACK, ierr )
    call MPI_TYPE_COMMIT(MPI_PHOTON_STACK, ierr)

    doSublimate = .true.
    if (PRESENT(sublimate)) doSublimate = sublimate
       


!Allocate memory for send buffer - fixes deadlocks
    !First, determine the no. of bytes used in a photonpacketstack array
    call MPI_PACK_SIZE(maxStackLimit, MPI_PHOTON_STACK, localWorldCommunicator, bufferSize, ierr)

    !Add some extra bytes for safety
    bufferSize = bufferCap*(bufferSize + MPI_BSEND_OVERHEAD)

    if(bufferSize < 0) then
       write(*,*) "warning negative buffer size", bufferSize
       bufferSize = -bufferSize
    end if

    allocate(buffer(bufferSize))



    allocate(nEscapedArray(1:nHydroThreadsGlobal))
    allocate(nSaved(nHydroThreadsGlobal))
    nescapedArray = 0    
    dprCounter = 0





!    if(.not. optimizeStack) then
       photonPacketStack%Freq = 0.d0
       photonPacketStack%Destination = 0
       photonPacketStack%tPhot = 0.d0
       photonPacketStack%crossedPeriodic = .false.
!    end if


       if(monochromatic) then
          nfreq = 2
          nuStart = (13.60001*evtoerg/hcgs)
          nuEnd =  (13.60001*evtoerg/hcgs)
          do i = 1, nFreq
             freq(i) = log10(nuStart) + dble(i-1)/dble(nFreq-1) * (log10(nuEnd)-log10(nuStart))
             freq(i) = 10.d0**freq(i)
             lams(nFreq-i+1) = cspeed/freq(i)
          enddo
          do i = 2, nFreq-1
             dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
          enddo
          dfreq(1) = 2.d0*(freq(2)-freq(1))
          dfreq(nfreq) = 2.d0*(freq(nfreq)-freq(nfreq-1))

       else
          nfreq = 1000
          nuStart = cSpeed / (1000.e4 * 1.d-8)
          nuEnd =  cSpeed / (10. * 1.d-8) ! 2.d0*maxval(grid%ion(1:grid%nIon)%nuThresh)
          do i = 1, nFreq
             freq(i) = log10(nuStart) + dble(i-1)/dble(nFreq-1) * (log10(nuEnd)-log10(nuStart))
             freq(i) = 10.d0**freq(i)
             lams(nFreq-i+1) = cspeed/freq(i)
          enddo
!          do i = 2, nFreq-1
!             dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
!          enddo
!          dfreq(1) = 2.d0*(freq(2)-freq(1))
!          dfreq(nfreq) = 2.d0*(freq(nfreq)-freq(nfreq-1))


          !H I, He I ionization edge refinement
          if(.not. hOnly) then
             k = 1
             do i = 3, 4
                !             if(i /= 2) then
                nuThresh = (grid%ion(i)%ipot)*evtoerg/hCgs
                freq(1000+k) = nuthresh*0.99d0
                k = k + 1
                freq(1000+k) = nuthresh*1.01d0
                k = k + 1
                nfreq = nfreq + 2
                !                end if
             end do
          end if
          call sort(nFreq, Freq)
          
          do i = 2, nFreq-1
             dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
          enddo
          dfreq(1) = 2.d0*(freq(2)-freq(1))
          dfreq(nfreq) = 2.d0*(freq(nfreq)-freq(nfreq-1))

!          print *, "nFreq ", nfreq

       end if

    if (firstTimeTables) then

       do i = 1, grid%nIon
          call addxSectionArray(grid%ion(i), nfreq, freq)
          !print *, "grid%ion(i)xsec",grid%ion(i)%xsec
     enddo
       call createGammaTable(gammaTableArray(1), 'gammaHI.dat')

       call createGammaTable(gammaTableArray(2), 'gammaHeI.dat')
       call createGammaTable(gammaTableArray(3), 'gammaHeII.dat')

       call createSahaMilneTables(hTable, heTable)

       call createRecombTable(Hrecombtable, 'e1b.d')

       call readHeIRecombinationLinesFit()

       call readHeIIrecombination()

       firstTimeTables = .false.
    endif

    lCore = 0.d0
    do i = 1, nSource
       !If outside the grid then correct for the flux attenuation due to distance of star from grid
       if(source(i)%outsideGrid) then
          lCore = lCore + source(i)%luminosity * (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 / &
               (fourPi*source(i)%distance**2)
       else
          lCore = lCore + source(i)%luminosity
       end if
    enddo

    if (myrankWorldglobal == 1) write(*,'(a,1pe12.5)') "Total source luminosity (lsol): ",lCore/lSol

    if (writeoutput) then
       write(*,'(a,1pe12.1)') "Ionizing photons per second ", ionizingFlux(source(1))
    endif

    nIonizingPhotons = ionizingFlux(source(1))

    !nmonte selector: Only works for fixed grids at present
    if (inputnMonte == 0) then
       
       if(.not. monteCheck .or. readGrid) then
          !waymaker photoionization loop
          if(grid%octreeRoot%twoD) then
                nMonte = int(4.d0**(maxDepthAMR))
             else if(grid%octreeRoot%threeD) then
                !nMonte = (8.d0**(maxDepthAMR))
                !nMonte = 5242880/2.
!                nMonte = 150000000
                nMonte = 209715200
               ! nMonte = 409600.0
             else
                !nMonte = 2**(maxDepthAMR)
                nMonte = 100000
             end if
          
       else 
          if(minDepthAMR == maxDepthAMR) then
             if(grid%octreeRoot%twoD) then
                nMonte = int(10.d0 * (4.d0**(maxDepthAMR)))
             else if(grid%octreeRoot%threeD) then
!                nMonte = 1000.d0
                nMonte = int(100.d0 * (8.d0**(maxDepthAMR)))
!                nMonte = 1.d0 * (8.d0**(maxDepthAMR))
                ! nMonte = 1.d0 * (8.d0**(maxDepthAMR))
                !nMonte = 5242880/2.
             else
                nMonte = int(500.d0 * 2**(maxDepthAMR))
             end if
          else
             call writeInfo("Non uniform grid, setting arbitrary nMonte", TRIVIAL)
             nMonte = 209715200
             write(*,*) "nMonte = ", nMonte
          end if
          
       end if
       
    else
       nMonte = inputnMonte
    endif

    nIter = 0

    converged = .false.
    if (nSource > 1) &
         call randomSource(source, nSource, iSource, photonPacketWeight, lamArray, nLambda, initialize=.true.)
    do while(.not.converged)


       if(optimizeStack .and. nIter > 0) then
          bufferSize = 0
          allocate(photonPacketStack(stackLimit*nHydroThreadsGlobal))
          allocate(toSendStack(stackLimit))
          allocate(currentStack(stackLimit))
          call MPI_PACK_SIZE(stackLimit, MPI_PHOTON_STACK, localWorldCommunicator, bufferSize, ierr)
          
          !Add some extra bytes for safety
          bufferSize = bufferCap*(bufferSize + MPI_BSEND_OVERHEAD)

          if(bufferSize < 0) then
             write(*,*) "warning negative buffer size", bufferSize
             print *, "bufferCap: ", bufferCap
             print *, "bufferSize: ", bufferSIze
             print *, "MPI_BSEND_OVERHEAD", MPI_BSEND_OVERHEAD
             bufferSize = -bufferSize
          end if

          allocate(buffer(bufferSize))
       end if
       if(customstacks) then
          zerothstacklimit = stacklimitArray(1)
       else
          zerothstacklimit = stacklimit
       end if

       !setup the buffer                
!       print *, "buffer ", buffer
!       print *, "bufferSize ", bufferSize
       call MPI_BUFFER_ATTACH(buffer,bufferSize, ierr)

       photonPacketStack%Freq = 0.d0
       photonPacketStack%Destination = 0
       photonPacketStack%tPhot = 0.d0
       photonPacketStack%crossedPeriodic = .false.
 
      nIter = nIter + 1
      nInf=0

       call clearContributions(grid%octreeRoot)

       if(monochromatic) then
          if (source(1)%outsidegrid) then
             
             epsoverdeltat = (((nIonizingPhotons*((13.60001)*evtoerg))*&
                  (2.d0*grid%octreeRoot%subcellSize*1.d10)**2)/dble(nMonte))                         
          else
             epsoverdeltat = (nIonizingPhotons*((13.60001)*evtoerg))/ &
                  dble(nMonte)
          endif
       else
          epsoverdeltat = lcore/dble(nMonte)
       end if
       
       if (myrankGlobal /= 0) call zeroDistanceGrid(grid%octreeRoot)

       if (myrankWorldGlobal == 1) write(*,*) "Running photoionAMR loop with ",nmonte," photons. Iteration: ",niter, maxIter

       if (myrankWorldGlobal == 1) call tune(6, "One photoionization itr")  ! start a stopwatch

       !Thaw - stack optimization
       if(optimizeStack) then
          call wallTime(startTime)
       end if

       nThreadMonte = nMonte / nHydroSetsGlobal
       iMonte_beg = 1
       iMonte_end = nMonte
       iMonte_end = nThreadMonte

       totalPower = 0.d0
       nTotScat = 0
       nPhot = 0
       nSaved = 0
       photonPacketStack%destination = 0
       photonPacketStack%freq = 0.d0
       toSendStack%freq = 0.d0
       toSendStack%destination = 0
       toSendStack%crossedPeriodic = .false.

       countArray = 0.d0

       if(optimizeStack .and. myRankGlobal == 0) then
          write(*,*) "DOING OPTIMIZATION"
          write(*,*) "StackLimit ", stacklimit
          write(*,*) "dstack", dstack
       end if

       if(customStacks .and. iniCustomTime) then
          write(*,*) "working with a custom stack array"
          write(*,*) "stack limit array: ", stacklimitarray
          iniCustomTime=.false.
       end if
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          if (myRankGlobal == 0) then
             mainloop: do iMonte = iMonte_beg, iMonte_end
                call randomSource(source, nSource, iSource, photonPacketWeight)
                thisSource = source(iSource)
                call getPhotonPositionDirection(thisSource, rVec, uHat,rHat,grid)

                !re-weighting for corner sources, edges still need work
                if(source(iSource)%onCorner) then
                   if (grid%octreeRoot%threeD) then
                      photonPacketWeight = photonPacketWeight * 1.d0/8.d0
                   else if (grid%octreeRoot%twoD) then
                      photonPacketWeight = photonPacketWeight * 1.d0 /4.d0
                   end if
                !else if(source(iSource)%onEdge .and. .not. source(iSource)%onCorner) then
                !   photonPacketWeight = 1.d0/2.d0
                else
                   photonPacketWeight = photonPacketWeight * 1.d0
                end if

                tPhoton = 0.d0

                if(monochromatic) then
                   thisFreq = ((13.60001)*evtoerg)/hcgs
                else
                   call getWavelength(thisSource%spectrum, wavelength, photonPacketWeight)                
                   thisFreq = cSpeed/(wavelength / 1.e8)
                end if

                totalPower = totalPower + epsOverDeltaT*photonPacketWeight
                
                if(binPhotons) then
                   do i = 1, nFreq
                      if((wavelength/1.e8) < lams(i+1) .and. (wavelength/1.e8) > lams(i)) then
                         countArray(i) = countArray(i) + 1.d0*photonPacketWeight
                         exit
                      end if
                   end do
                end if

                call findSubcellTD(rVec, grid%octreeRoot,thisOctal, subcell)
                iThread = thisOctal%mpiThread(subcell)

                !Create a bundle of photon packets, only modify the first available array space
                !available slots have frequency == 0.d0
                do optCounter = 1, (SIZE(photonPacketStack))
                   if(photonPacketStack(optCounter)%freq == 0.d0) then
                      photonPacketStack(optCounter)%rVec = rVec
                      photonPacketStack(optCounter)%uHat = uHat
                      photonPacketStack(optCounter)%freq = thisFreq
                      photonPacketStack(optCounter)%tPhot = tPhoton
                      photonPacketStack(optCounter)%ppw = photonPacketWeight
                      photonPacketStack(optCounter)%destination = iThread
                      photonPacketStack(optCounter)%sourcePhoton = .true.
                      photonPacketStack(optCounter)%crossedPeriodic = .false.
                      exit
                   end if
                end do

                !Keep track of how many packets are destined for each thread
                nSaved(iThread) = nSaved(iThread) + 1

                !Once the bundle for a specific thread has reached a critical size, send it to the thread for propagation
                do optCounter = 1, nHydroThreadsGlobal
                   if(nSaved(optCounter) /= 0) then
                      if(nSaved(optCounter) == (zerothstackLimit) .or. &
                           (nThreadMonte - nInf) < (zerothstackLimit*nHydroThreadsGlobal)) then
                         thisPacket = 1
                         do sendCounter = 1, (maxStackLimit*nHydroThreadsGlobal)
                            if(photonPacketStack(sendCounter)%destination == optCounter &
                                 .and. photonPacketStack(sendCounter)%freq /= 0.d0) then           
                               toSendStack(thisPacket)%rVec = photonPacketStack(sendCounter)%rVec
                               toSendStack(thisPacket)%uHat = photonPacketStack(sendCounter)%uHat
                               toSendStack(thisPacket)%freq = photonPacketStack(sendCounter)%freq
                               toSendStack(thisPacket)%tPhot = photonPacketStack(sendCounter)%tPhot
                               toSendStack(thisPacket)%ppw = photonPacketStack(sendCounter)%ppw
                               toSendStack(thisPacket)%destination = photonPacketStack(sendCounter)%destination
                               toSendStack(thisPacket)%sourcePhoton = photonPacketStack(sendCounter)%sourcePhoton
                               toSendStack(thisPacket)%crossedPeriodic = photonPacketStack(sendCounter)%crossedPeriodic 
                               thisPacket = thisPacket + 1
                               nInf = nInf + 1
                               
                               !Reset the photon frequency so that this array entry can be overwritten
                               photonPacketStack(sendCounter)%freq = 0.d0
                               
                               !Reset the destination so that it doesnt get picked up in subsequent sweeps!       
                               photonPacketStack(sendCounter)%destination = 0
                            
                            end if
                         end do

                         call MPI_SEND(toSendStack, zerothstackLimit, MPI_PHOTON_STACK, &
                              OptCounter, tag, localWorldCommunicator,  ierr)

                         !reset the counter for this thread's bundle recieve 
                         nSaved(optCounter) = 0
                         toSendStack%destination = 0
                         toSendStack%freq = 0.d0
                      end if
                   end if
                end do

             end do mainloop

             if (myrankGlobal == 0) then
                call mpi_allreduce(totalpower, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_SUM, zeroThreadCommunicator, ierr)
                totalPower = tempDouble
             endif

             if (myrankWorldGlobal == 0) then
                write(*,*) "lcore = ", lcore
                write(*,*) "totalPower = ", totalPower
                write(*,*) "ratio = ", lcore/totalPower
                write(*,*) "nmonte ",nMonte
                write(*,*) "nthreadmonte ",nThreadMonte
                write(*,*) "epsoverdt ",epsOverDeltaT
                write(*,*) "lcore / nMonte ",lCore/dble(nMonte)
                write(*,*) "lcore / nThreadMonte ",lCore/dble(nThreadMonte)
                write(*,*) "Rank 0 sent all initial bundles"
                write(*,*) "Telling Ranks to pass stacks ASAP "
             endif

             if(binPhotons) then
                open(333, file="bins.dat", status="unknown")
                do i=1, nFreq
                   write(333,*) lams(i)*1.e8, countArray(i)/dble(nMonte)
                end do
                close(333)
             end if
                
                do iThread = 1, nHydroThreadsGlobal
                   toSendStack(1)%destination = 500
                   call MPI_SEND(toSendStack, maxStackLimit, MPI_PHOTON_STACK, iThread, tag, localWorldCommunicator,  ierr)
                   call MPI_RECV(donePanicking, 1, MPI_LOGICAL, iThread, tag, localWorldCommunicator, status, ierr)
                   
                end do


             photonsStillProcessing = .true.
             
             do while(photonsStillProcessing)                   
                toSendStack%freq = 0.d0
                do iThread = 1, nHydroThreadsGlobal
                   toSendStack(1)%destination = 999
                   !toSendStack(1)%freq = 100.d0
                   
                   call MPI_SEND(toSendStack, maxStackLimit, MPI_PHOTON_STACK, iThread, tag, localWorldCommunicator,  ierr)
                   
                   call MPI_RECV(nEscapedArray(iThread), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, status, ierr)
                   
                enddo
                
                nEscaped = SUM(nEscapedArray(1:nHydroThreadsGlobal))
             
                nEscapedGlobal = nEscaped
                if (nEscaped == nThreadMonte) then
                   photonsStillProcessing = .false.
                else if(nEscaped > nThreadMonte) then
                   write(*,*) "nEscaped greater than nMonte, Exiting..."
                   do iThread = 1, nHydroThreadsGlobal
                      write(*,*) "nEscapedArray(iThread) ", iThread, nEscapedArray(iThread)
                   end do
                   stop
                end if
!                print *, "nEscaped ", nEscaped
             end do

             if (myrankWorldGlobal == 1) write(*,*) "Finishing iteration..."


             do iThread = 1, nHydroThreadsGlobal
                !tosendstack(1)%freq = 100.d0
                toSendStack(1)%destination = 888
!                toSendStack(2)%destination = 999
                call MPI_SEND(toSendStack, maxStackLimit, MPI_PHOTON_STACK, iThread, tag, localWorldCommunicator,  ierr)
             enddo

          else
             endLoop = .false.
             nEscaped = 0
             photonPacketStack%freq = 0.d0
             currentStack%freq = 0.d0
             stackSize = 0
             nSaved = 0
             sendAllPhotons = .false.
             !needNewPhotonArray = .true.  
             do while(.not.endLoop) 
                crossedMPIboundary = .false.

                !Get a new photon stack
                iSignal = -1
                if(stackSize == 0) then
                     
                   call MPI_RECV(toSendStack, maxStackLimit, MPI_PHOTON_STACK, MPI_ANY_SOURCE, &
                        tag, localWorldCommunicator, status, ierr)
                      currentStack = toSendStack
                      

                      !Check to see how many photons in stack are not null
                      do p = 1, maxStackLimit
                         if(currentStack(p)%freq /= 0.d0) then                            
                            stackSize = stackSize + 1
                         end if
                      end do

                   !Check to see if a special action is required
                   escapeCheck = .false.
                   !Escape checking
                   if(currentStack(1)%destination == 999) then
                      iSignal = 1
                      escapeCheck = .true.
                      stackSize = 0
                      !End photoionization loop
                   else if (currentStack(1)%destination == 888) then
                      iSignal = 0                    
                      currentStack%destination = 0                      
                   else if(currentStack(1)%destination == 500 .and. .not. sendAllPhotons) then
                      sendAllPhotons = .true.
                      stackSize = 0
                      !Evacuate everything currently in need of sending
                      do optCounter = 1, nHydroThreadsGlobal
                         if(optCounter /= myRankGlobal .and. nSaved(optCounter) /= 0) then
                            if(nSaved(optCounter) == sendStackLimit .or. sendAllPhotons) then
                               thisPacket = 1
                               
                               toSendStack%freq = 0.d0
                               do sendCounter = 1, (maxStackLimit*nHydroThreadsGlobal)
                                  if(photonPacketStack(sendCounter)%destination /= 0 .and. &
                                       photonPacketStack(sendCounter)%destination == optCounter) then
                                     
                                     toSendStack(thisPacket)%rVec = photonPacketStack(sendCounter)%rVec
                                    toSendStack(thisPacket)%uHat = photonPacketStack(sendCounter)%uHat
                                    toSendStack(thisPacket)%freq = photonPacketStack(sendCounter)%freq
                                    toSendStack(thisPacket)%tPhot = photonPacketStack(sendCounter)%tPhot
                                    toSendStack(thisPacket)%ppw = photonPacketStack(sendCounter)%ppw
                                    toSendStack(thisPacket)%destination = photonPacketStack(sendCounter)%destination
                                    toSendStack(thisPacket)%sourcePhoton = photonPacketStack(sendCounter)%sourcePhoton
                                    toSendStack(thisPacket)%crossedPeriodic = photonPacketStack(sendCounter)%crossedPeriodic
                                    thisPacket = thisPacket + 1
                                    !nInf = nInf + 1
                                    
                                    !Reset the photon frequency so that this array entry can be overwritten
                                    photonPacketStack(sendCounter)%freq = 0.d0
                                    
                                     !Reset the destination so that it doesnt get picked up in subsequent sweeps!
                                    photonPacketStack(sendCounter)%destination = 0
                                    
                                 end if
                              end do
                              call MPI_SEND(toSendStack, maxStackLimit, &
                                   MPI_PHOTON_STACK, OptCounter, tag, localWorldCommunicator,  ierr)
!                              call MPI_BSEND(toSendStack, stackLimit, MPI_PHOTON_STACK, OptCounter, tag, localWorldCommunicator, &
!                                   request, ierr)
                              !reset the counter for this thread's bundle recieve
                              nSaved(optCounter) = 0
                              toSendStack%freq = 0.d0
                              toSendStack%destination = 0
                           end if
                        end if
                     end do
                     call MPI_SEND(donePanicking, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator, ierr)
                     
                     goto 777
                  end if
                  Currentstack%destination = 0
               else if (stackSize < 0) then
                  write(*,*) "StackSize less than zero... aborting.", myRankGlobal
                  stop
               end if
               
               if (iSignal == 0) then
                  endLoop = .true.
                  goto 777
               endif
               if (iSignal == 1) then
                  call MPI_SEND(nEscaped, 1, MPI_INTEGER, 0, tag, localWorldCommunicator,  ierr)
                  goto 777
               endif

               !$OMP PARALLEL DEFAULT(NONE) &
               !$OMP PRIVATE(p, rVec, uHat, thisFreq, tPhoton, photonPacketWeight, sourcePhoton) &
               !$OMP PRIVATE(escaped, nScat, optCounter, octVec, ierr, thisLam, kappaabsdb) &
               !$OMP PRIVATE(kappascadb, albedo, r, kappaabsdust, thisOctal, subcell, sendStackLimit) &
               !$OMP PRIVATE(crossedMPIboundary, newThread, thisPacket, kappaabsgas, escat, tempcell ) &
               !$OMP PRIVATE(r1, finished, voidThread, crossedPeriodic, nperiodic, request, myrankworldglobal) &
               !$OMP SHARED(photonPacketStack, myRankGlobal, currentStack, escapeCheck) &
               !$OMP SHARED(tag, noDiffuseField, grid, epsoverdeltat, iSignal, MPI_PHOTON_STACK) &
               !$OMP SHARED(nlambda, lamarray, tlimit, nHydroThreadsGlobal, sendAllPhotons,toSendStack) &
               !$OMP SHARED(nTotScat, gammaTableArray, freq) &
               !$OMP SHARED(dfreq, iLam, endLoop, nIter, spectrum) &
               !$OMP SHARED(nSaved, maxStackLimit) &
               !$OMP SHARED(stackSize, nFreq) &
               !$OMP SHARED(nPhot, nEscaped, stackLimit, localWorldCommunicator, nhydrosetsglobal)
               
               finished = .false.
               escaped = .false.
               crossedMPIboundary=.false.
               voidThread = .false.
               
               !Get the next available photon in the stack for processing
               !$OMP CRITICAL (get_photon_stack)
               if(stackSize /= 0 .and. .not. escapeCheck) then
                  do p = 1, maxStackLimit
                     if(currentStack(p)%freq /= 0.d0) then
                        rVec = currentStack(p)%rVec
                        uHat = currentStack(p)%uHat
                        thisFreq = currentStack(p)%Freq
                        tPhoton = currentStack(p)%tPhot
                        photonPacketWeight = currentStack(p)%ppw
                        sourcePhoton = currentStack(p)%sourcePhoton
                        crossedPeriodic = currentStack(p)%crossedPeriodic
                        currentStack(p)%freq = 0.d0
                        exit
                     end if
                  end do
                  stackSize = stackSize - 1
               else
                  voidThread = .true.
               end if
               !$OMP END CRITICAL (get_photon_stack)

                ! here a new photon has been received, and we add its momentum

                call findSubcellTD(rVec, grid%octreeRoot,thisOctal, subcell)
                thisOctal%radiationMomentum(subcell) = thisOctal%radiationMomentum(subcell) + uHat * (epsOverDeltaT/cSpeed)
!                if (myrankWorldGlobal == 1) write(*,*) "mom add new ",thisOctal%radiationMomentum(subcell)

                if (.not.endLoop) then
                   nScat = 0
                   escaped = .false.
                   do while(.not.escaped)

                      call toNextEventPhoto(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray, &
                           photonPacketWeight, epsOverDeltaT, nfreq, freq, tPhoton, tLimit, &
                           crossedMPIboundary, newThread, sourcePhoton, crossedPeriodic)

                      if (crossedMPIBoundary) then    
                         !Create a bundle of photon packets, only modify the first available array space
                         !$OMP CRITICAL (crossed_mpi_boundary)
                         do optCounter = 1, (SIZE(photonPacketStack))
                            if(photonPacketStack(optCounter)%freq == 0.d0) then
                               photonPacketStack(optCounter)%rVec = rVec
                               photonPacketStack(optCounter)%uHat = uHat
                               photonPacketStack(optCounter)%freq = thisFreq
                               photonPacketStack(optCounter)%tPhot = tPhoton
                               photonPacketStack(optCounter)%ppw = photonPacketWeight
                               photonPacketStack(optCounter)%destination = newThread
                               photonPacketStack(optCounter)%sourcePhoton = sourcePhoton
                               photonPacketStack(optCounter)%crossedPeriodic = crossedPeriodic
                               if(crossedPeriodic) then
                                  nPeriodic = nPeriodic + 1
                               end if
                               exit
                            end if
                         end do

                         !Keep track of how many packets are destined for each thread
                         nSaved(newThread) = nSaved(newThread) + 1
                         
                         !Once the bundle for a specific thread has reached a critical size, send it to the thread for propagation
                         do optCounter = 1, nHydroThreadsGlobal
                            if(optCounter /= myRankGlobal .and. nSaved(optCounter) /= 0) then
                               if(nSaved(optCounter) == (sendStackLimit) .or. sendAllPhotons) then
                                  thisPacket = 1
                                  toSendStack%freq = 0.d0
                                  do sendCounter = 1, (maxStackLimit*nHydroThreadsGlobal)
                                     if(photonPacketStack(sendCounter)%destination /= 0 .and. &
                                          photonPacketStack(sendCounter)%destination == optCounter) then 
                                        
                                        toSendStack(thisPacket)%rVec = photonPacketStack(sendCounter)%rVec
                                        toSendStack(thisPacket)%uHat = photonPacketStack(sendCounter)%uHat
                                        toSendStack(thisPacket)%freq = photonPacketStack(sendCounter)%freq
                                        toSendStack(thisPacket)%tPhot = photonPacketStack(sendCounter)%tPhot
                                        toSendStack(thisPacket)%ppw = photonPacketStack(sendCounter)%ppw
                                        toSendStack(thisPacket)%destination = photonPacketStack(sendCounter)%destination
                                        toSendStack(thisPacket)%sourcePhoton = photonPacketStack(sendCounter)%sourcePhoton
                                        toSendStack(thisPacket)%crossedPeriodic = photonPacketStack(sendCounter)%crossedPeriodic
                                        thisPacket = thisPacket + 1
                                        !nInf = nInf + 1
                                     
                                        !Reset the photon frequency so that this array entry can be overwritten
                                        photonPacketStack(sendCounter)%freq = 0.d0
                                     
                                        !Reset the destination so that it doesnt get picked up in subsequent sweeps!
                                        photonPacketStack(sendCounter)%destination = 0
                                        
                                     end if
                                  end do

                                  if(thisPacket /= 1 .and. nSaved(optCounter) /= 0) then
!                                     call MPI_SEND(toSendStack, stackLimit, MPI_PHOTON_STACK, OptCounter, tag, localWorldCommunicator, &
!                                          ierr)
                                     call MPI_BSEND(toSendStack, maxStackLimit, MPI_PHOTON_STACK, OptCounter, tag, &
                                          localWorldCommunicator, request,ierr)

                                     nSaved(optCounter) = 0
                                     toSendStack%freq = 0.d0
                                     toSendStack%destination = 0
                                  end if
                                  
                               end if
                            end if
                         end do
                         !$OMP END CRITICAL (crossed_mpi_boundary)
                         finished = .true.
                      endif


                      if(.not. finished) then
                         if (noDiffuseField) escaped = .true.

                         if (escaped) then
                            !$OMP CRITICAL (update_escaped)
                            nEscaped = nEscaped + 1
                            !$OMP END CRITICAL (update_escaped)
                         end if

                         if (.not. escaped) then
                            thisLam = real(cSpeed / thisFreq) * 1.e8
                            call locate(lamArray, nLambda, real(thisLam), iLam)
                            octVec = rVec 
                            call amrGridValues(grid%octreeRoot, octVec, foundOctal=thisOctal, foundsubcell=subcell,iLambda=iLam, &
                                 lambda=real(thisLam), kappaSca=kappaScadb, kappaAbs=kappaAbsdb, grid=grid)
                            
                            albedo = kappaScaDb / (kappaAbsdb + kappaScadb)
                            
                            
                            call randomNumberGenerator(getDouble=r)
                            if (r < albedo) then
                               uHat = randomUnitVector() ! isotropic scattering
                               !                            if (myrank == 1) write(*,*) "new uhat scattering ",uhat
                            else
                               
                               
                               spectrum = 1.d-30
                               
                               call returnKappa(grid, thisOctal, subcell, ilambda=ilam, &
                                    kappaAbsDust=kappaAbsDust, kappaAbsGas=kappaAbsGas, &
                                    kappaSca=kappaScadb, kappaAbs=kappaAbsdb, kappaScaGas=escat)
                                               
                               

               
                               if ((thisFreq*hcgs*ergtoev) > 13.6) then ! ionizing photon
                                  call randomNumberGenerator(getDouble=r1)
                                  
                                  if (r1 < (kappaAbsGas / max(1.d-30,(kappaAbsGas + kappaAbsDust)))) then  ! absorbed by gas rather than dust
                                     call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
                                     call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
                                     call addHydrogenRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
                                     !                        call addHeRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
                                     call addForbiddenLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
                                  else
                                     call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, & 
                                          nlambda, lamArray)
                                  endif
                                  !Subsequent progress is now a diffuse photon
                                  sourcePhoton = .false.
                               else ! non-ionizing photon must be absorbed by dust
                                  call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
                               endif
                               if (firsttime.and.(myrankWorldglobal==1)) then
                                  firsttime = .false.
                                  open(67,file="pdf.dat",status="unknown",form="formatted")
                                  do i = 1, nfreq
                                     write(67,*) freq(i), spectrum(i)
                                  enddo
                                  close(67)
!                                  open(68,file="freq.dat",status="unknown",form="formatted")
!                                  do i = 1, nfreq
!                                     write(68,*) i, freq(i)
!                                  enddo
!                                  close(68)
                               endif
                               thisFreq =  getPhotonFreq(nfreq, freq, spectrum)
                               uHat = randomUnitVector() ! isotropic emission
                                                       
                               nScat = nScat + 1
                               
                                                        !                            if ((thisFreq*hcgs*ergtoev) < 13.6) escaped = .true.
                               
                               
                            endif
                         endif
                         if (nScat > 100000) then
                            write(*,*) "Nscat exceeded 10000, forcing escape"
                            write(*,*) 1.e8*cspeed/thisFreq
                            write(*,*) albedo, kappaScaDb, kappaAbsdb,escat
                            
                            thisLam = real(cSpeed / thisFreq) * 1.e8
                            call locate(lamArray, nLambda, real(thisLam), iLam)
                            
                            call returnKappa(grid, thisOctal, subcell, ilambda=ilam, &
                                 Kappaabsdust=kappaAbsDust, kappaAbsGas=kappaAbsGas, &
                                 kappaSca=kappaScadb, kappaAbs=kappaAbsdb, kappaScaGas=escat)
                            write(*,*) "thislam",thislam,ilam, lamArray(ilam)
                            write(*,*) lamArray(1:nLambda)
                            write(*,*) "kappaAbsDust",kappaAbsDust
                            write(*,*) "kappaAbsGas",kappaAbsGas
                            write(*,*) "kappaSca",kappaScadb
                            write(*,*) "kappaAbs",kappaAbsdb
                            write(*,*) "onekappaabs",grid%oneKappaAbs(1,ilam)
                            write(*,*) "onekappasca",grid%oneKappasca(1,ilam)
                  
                            escaped = .true.
                            
                         Endif
                      end if
                   enddo
                   nTotScat = nTotScat + nScat
                   nPhot = nPhot + 1
                end if
                !$OMP BARRIER
                !$OMP END PARALLEL
777             continue
             enddo


!             print *, "Rank ", myRank, "has nPeriodic : ", nPeriodic
             nPeriodic = 0

          endif

       if (myrankWorldGlobal == 1) call tune(6, "One photoionization itr")  ! stop a stopwatch

       !Get the time for the iteration and see if it has improved with a new stack size
       if(optimizeStack .and. .not. optConverged .and. myRankGlobal /= 0) then  
             call MPI_BARRIER(amrCommunicator, ierr)
             call wallTime(endTime)

             !Make sure we are working with the same time
             if(myRankGlobal == 1) then
                do optCount = 2, nHydroThreadsGlobal
                   call MPI_SEND(endTime, 1, MPI_INTEGER, optCount, tag, localWorldCommunicator,  ierr)
                end do
             else 
                call MPI_RECV(endTime, 1, MPI_INTEGER, 1, tag, localWorldCommunicator, status, ierr)
             end if

             oldStackLimit = stackLimit

             !get the time per photon packet
             newTime = endTime - startTime
             newTime = newTime/real(nthreadmonte)

             if(optID == 1) then
                optID = 2 
                iterTime(1) = newTime                
                iterStack(1) = stackLimit
                stackLimit = oldStackLimit + dStack
                
             else if(optID == 2) then
                optID = 3
                iterTime(2) = newTime
                iterStack(2) = stackLimit
                stackLimit = iterStack(1) - dStack

             else if(optID == 3) then
                optID = 1
                iterTime(3) = newTime
                iterStack(3) = stackLimit

                !do stuff
                !get the gradients
                m1 = (dble(iterTime(2) - iterTime(1)))/(dble(iterStack(2) - iterStack(1)))
                m2 = (dble(iterTime(3) - iterTime(1)))/(dble(iterStack(3) - iterStack(1)))

                if(m1 > 0.d0) then
                   !bigger stack took longer
                   if(m2 > 0.d0) then
                      !smaller stack was quicker
                      stackLimit = iterStack(3)
                   else
                      !smaller stack took longer/the same
                      stackLimit = iterStack(1)
                   end if
                else if(m1<0.d0) then
                   !bigger stack was quicker
                   if(m2 < 0.d0) then
                      !smaller stack was slower
                      stackLimit = iterStack(2)
                   else 
                      !smaller stack was also quicker
                      if(abs(m1) > abs(m2)) then
                         !but big stack is better
                         stackLimit = iterStack(2)
                      else
                         !and smaller stack is the best
                         stackLimit = iterStack(3)
                      end if
                   end if
                else
                   !bigger stack was the same
                   if(m2 == 0.d0) then
                      !smaller stack was the same
                      stackLimit = iterStack(1)
                   else if(m2 > 0.d0) then
                      !smaller stack was faster
                      stackLimit = iterStack(3)
                   else
                      !smaller stack was slower
                      stackLimit = iterStack(1)
                   end if
                end if

                !reset for next probing
                iterTime = 0.0
                iterStack = 0
             else
                call torus_abort("optID error")
             end if

!             print *, "startTime ", startTime
!             print *, "endTime ", endTime
!             print *, "newTime ", newTime
!             print *, "iterTime ", iterTime
             !
          if(stacklimit <= 0) then
             stacklimit = 1
          end if
          if(myRankGlobal == 1) then
             call MPI_SEND(stackLimit, 1, MPI_INTEGER, 0, tag, localWorldCommunicator,  ierr)
          end if
       end if
       
       if(myRankGlobal == 0 .and. optimizeStack) then
          call MPI_RECV(stackLimit, 1, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, status, ierr)
       end if

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!       epsOverDeltaT = (lCore) / dble(nMonte)
       if(monochromatic) then
          if (source(1)%outsidegrid) then

             epsoverdeltat = (((nIonizingPhotons*((13.60001)*evtoerg))*&
                  (2.d0*grid%octreeRoot%subcellSize*1.d10)**2)/dble(nMonte))

          else
             epsoverdeltat = (nIonizingPhotons*((13.60001)*evtoerg))/ &
                  dble(nMonte)
          endif
       else
          epsoverdeltat = lcore/dble(nMonte)
!          print *, "epsAlarm"
       end if
 

       if (myrankGlobal /= 0) then
          call MPI_ALLREDUCE(nTotScat, i, 1, MPI_INTEGER, MPI_SUM, amrCommunicator, ierr)
          nTotScat = i
          if (writeoutput) write(*,*) "Iteration had ",real(nTotScat)/real(nThreadMonte), " scatters per photon"
       endif


    np = 1
    firstTime = .true.

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif

    call  identifyUndersampled(grid%octreeRoot)

! now we need to reduce the estimators across all the hydro sets.

    if (nHydroSetsGlobal > 1) then
       do iThread = 1, nHydroThreadsGlobal
          if (myRankGlobal == iThread) then
!             write(*,*) myHydroSetGlobal, myrankGlobal, " calling updategridMpi"
             call updateGridMPIphoto(grid, amrParallelCommunicator(iThread))
          endif
          call mpi_barrier(MPI_COMM_WORLD, ierr)
       enddo
    endif



    if (myRankWorldGlobal == 0) write(*,*) "Ninf ",ninf
    if (myrankWorldGlobal == 1) call tune(6, "Temperature/ion corrections")

    if (writeoutput) &
         write(*,*) "Calculating ionization and thermal equilibria"

    if (myrankWorldGlobal == 1) then
       i = 0
       j = 0 
       call testforZero(grid%octreeRoot, i, j)
       write(*,*) "photoion(subcell,1) is zero for ", 100.d0*real(i)/real(j), "% of subcells"
    endif


!    call writeVtkFile(grid, "beforethermalcalc.vtk", &
!         valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!         "hydrovelocity","sourceCont   ","pressure     ", &
!         "crossings    "/))!


    if (myRankGlobal /= 0) then

       ! default loop indices
       ioctal_beg = 1
       ioctal_end = nOctal

       if (nHydroSetsGlobal > 1) then
          np = nHydroSetsGlobal
          n_rmdr = MOD(nOctal, np)
          mOctal = nOctal/np
          
          if (myHydroSetGlobal .lt. n_rmdr ) then
             iOctal_beg = (mOctal+1)*myHydroSetGlobal + 1
             iOctal_end = iOctal_beg + mOctal
          else
             iOctal_beg = mOctal*myHydroSetGlobal + 1 + n_rmdr
             iOctal_end = iOctal_beg + mOctal - 1
          end if
          allocate(tempCell(1:nOctal,1:8))
          tempCell = 0.
       endif

       !$OMP PARALLEL DEFAULT(NONE) &
       !$OMP PRIVATE(iOctal, thisOctal, subcell, v, kappap, i) &
       !$OMP PRIVATE(dustHeating, tempcell) &
       !$OMP SHARED(iOctal_beg, iOctal_end, dustOnly, octalArray, grid, epsOverDeltaT) &
       !$OMP SHARED(timedep, quickThermal, deltaTime, tminGlobal, myrankGlobal, nhydrosetsglobal)

       !$OMP DO SCHEDULE(DYNAMIC,2)
       do iOctal =  iOctal_beg, iOctal_end
          
          thisOctal => octalArray(iOctal)%content
          
          if (dustOnly) then
             do subcell = 1, thisOctal%maxChildren
                if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
                if (.not.thisOctal%hasChild(subcell)) then
                   v = cellVolume(thisOctal, subcell)
                   call returnKappa(grid, thisOctal, subcell, kappap=kappap)
                   dustHeating = (epsOverDeltaT / (v * 1.d30))*thisOctal%distanceGrid(subcell) ! equation 14 of Lucy 1999
                   thisOctal%temperature(subcell) = max(tMinGlobal,real((pi/stefanBoltz) * dustHeating / (fourPi * kappaP))**0.25e0)
                endif
             enddo

          else             
             do i = 1, nTimes
                if(timeDep) then
                   call calculateIonizationBalanceTimeDep(grid,thisOctal, epsOverDeltaT, deltaTime/dble(nTimes))
                   if (quickThermal) then
                      call quickThermalCalc(thisOctal)
                   else
                      call calculateThermalBalanceTimeDep(grid, thisOctal, epsOverDeltaT, deltaTime)
                   end if
                else
                   call calculateIonizationBalance(grid,thisOctal, epsOverDeltaT)
                   if (quickThermal) then
                      call quickThermalCalc(thisOctal)
                   else
                      call calculateThermalBalance(grid, thisOctal, epsOverDeltaT)

                   endif
                end if
             enddo

             
          endif
          if (nHydroSetsGlobal > 1) tempCell(iOctal,:) = thisOctal%temperature(:)
       enddo
       !$OMP END DO
       !$OMP END PARALLEL

       if (nHydroSetsGlobal > 1) then
          allocate(temp1(1:(nOctal*8)), temp2(1:nOctal*8))
          temp1 = RESHAPE(tempCell, (/SIZE(temp1)/))
          call MPI_ALLREDUCE(temp1, temp2, SIZE(temp1), MPI_REAL, MPI_SUM, amrParallelCommunicator(myRankGlobal),ierr)
          tempCell = RESHAPE(temp2,SHAPE(tempCell))
          do iOctal = 1, nOctal
             octalArray(iOctal)%content%temperature(:) = tempCell(iOctal,:)
          enddo
          deallocate(temp1, temp2, tempCell)
       endif

    endif

    !deallocate(octalArray)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (writeoutput) &
         write(*,*) "Finished calculating ionization and thermal equilibria"
    if (myrankWorldGlobal == 1) call tune(6, "Temperature/ion corrections")



    if(grid%geometry == "lexington") then
      call dumpLexingtonMPI(grid, epsoverdeltat, niter)
   end if

    if(grid%geometry == "point") then
       write(mpiFilename,'(a, i4.4, a)') "point.grid"
       call writeAmrGrid(mpiFilename, .false., grid)
    end if


    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    thisThreadConverged = .false.    
    failed = .false.

     anyUndersampled = .false.
     if(grid%geometry == "hii_test") then
        minCrossings = 1000
     else if(grid%geometry == "lexington") then
        minCrossings = 50000
     else if(singleMegaPhoto) then
        minCrossings = 50000 
     else
!        minCrossings = 5000
!        minCrossings = 1000
     end if

   !Thaw - auto convergence testing I. Temperature, will shortly make into a subroutine
       if (myRankGlobal /= 0) then
          iOctal_beg = 1
          iOctal_end = nOctal
          do iOctal =  iOctal_beg, iOctal_end
             thisOctal => octalArray(iOctal)%content
             do subcell = 1, thisOctal%maxChildren
                if (.not.thisOctal%hasChild(subcell)) then

                      if (niter == 1) then
                         thisOctal%TLastIter(subcell) = thisOctal%temperature(subcell)
                         thisThreadConverged = .false.
                         
                      else
                         
                         deltaT = (thisOctal%temperature(subcell)-thisOctal%TLastIter(subcell)) / &
                              thisOctal%TLastIter(subcell)  
                         deltaT = abs(deltaT)                 
                                                  
                         if(deltaT > 5.0d-2) then
                            if (thisOctal%nCrossings(subcell) /= 0 .and. thisOctal%nCrossings(subcell) < minCrossings) then
                               anyUndersampled = .true.
                            endif
                         end if
                         
                         if(deltaT < 5.0d-2 .and. .not. failed) then
                            thisThreadConverged = .true.
                         else 
                            if(niter > 2) then
                               fluctuationCheck = abs((thisOctal%temperature(subcell)-thisOctal%TLastLastIter(subcell))/ &
                                    thisOctal%TLastLastIter(subcell))
                               
                               if(fluctuationCheck < 5.0d-2 .and. .not. failed) then!
                                  thisThreadConverged = .true.
                               else
                                  thisThreadConverged = .false.                             
!                                  if(deltaT /= 0.d0 .and. .not. failed) then
!                                     write(*,*) "deltaT = ", deltaT
!                                     write(*,*) "thisOctal%temperature(subcell) ", thisOctal%temperature(subcell)
!                                     write(*,*) "thisOctal%TLastIter(subcell) ", thisOctal%TLastIter(subcell)
!                                     write(*,*) "thisOctal%TLastLastIter(subcell) ", thisOctal%TLastLastIter(subcell)
!                                     write(*,*) "cell center ", subcellCentre(thisOctal,subcell)
!                                     write(*,*) "nCrossings ", thisOctal%nCrossings(subcell)
!                                  end if
                                  failed = .true.
                               end if
                            else
                               thisThreadConverged = .false.  
!                               if(deltaT /= 0.d0 .and. .not. failed) then
!                                  write(*,*) "deltaT = ", deltaT
!                                  write(*,*) "thisOctal%temperature(subcell) ", thisOctal%temperature(subcell)
!                                  write(*,*) "thisOctal%TLastIter(subcell) ", thisOctal%TLastIter(subcell)
!                                  write(*,*) "cell center ", subcellCentre(thisOctal,subcell)
!                                  write(*,*) "nCrossings ", thisOctal%nCrossings(subcell)
!                               end if
                               failed = .true.
                            end if
                         end if
                         
                         !Check for temperature oscillations
                         if(niter > 1) then
                            thisOctal%TLastLastIter(subcell) = thisOctal%TLastIter(subcell)
                         end if
                         thisOctal%TLastIter(subcell) = thisOctal%temperature(subcell)
                      end if
                   end if
                   
                end do
                
                !Send result to master rank 
             end do
             
             !Send converged information
             call MPI_SEND(thisThreadConverged , 1, MPI_LOGICAL, 0, tag, localWorldCommunicator, ierr)
             call MPI_SEND(anyUndersampled, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator, ierr)

          else
             failed = .false.
             underSamFailed = .false.
             underSampledTOT = .false.
             
!!!Rank 0, collate results and decide if converged 
             converged = .false.
             do iThread = 1 , nHydroThreadsGlobal
                call MPI_RECV(thisThreadConverged,1, MPI_LOGICAL, iThread, tag, localWorldCommunicator, status, ierr )
                call MPI_RECV(anyUndersampled, 1, MPI_LOGICAL, iThread, tag, localWorldCommunicator, status, ierr)
                if(thisThreadConverged .and. .not. failed) then
                   converged = .true.
                else
                   converged = .false.
                   failed = .true.
                end if
                if(anyUndersampled) then
                   undersampledTOT=.true.
                end if
             end do
             
             do iThread = 1, nHydroThreadsGlobal
                call MPI_SEND(converged, 1, MPI_LOGICAL, iThread, tag, localWorldCommunicator, ierr) 
                call MPI_SEND(underSampledTOT, 1, MPI_LOGICAL, iThread, tag, localWorldCommunicator, ierr)   
             end do
          end if
          
          if(myRankGlobal /= 0) then
             call MPI_RECV(converged, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator, status, ierr)
             call MPI_RECV(underSampledTOT, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator, status, ierr)
             
          end if
          
     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     
     deallocate(octalArray)    
     
     !     converged = .false.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (niter >= maxIter) converged = .true. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if(converged) then
        if(myRankWorldGlobal == 0) then
           write(*,*) "photoionization loop converged at iteration ", niter
        end if
     else if(underSampledTOT .and. monteCheck) then
        if(myRankWorldGlobal == 0) then
           write(*,*) "Undersampled cell, increasing nMonte"
        end if
        nMonte = nMonte *2
     end if
     

!     write(mpiFilename,'(a, i4.4, a)') "photo", nIter,".vtk"!
!     call writeVtkFile(grid, mpiFilename, &
!          valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!          "hydrovelocity","sourceCont   ","pressure     ", &
!          "crossings    ", &
!          "dust1        "/))

     if(singleMegaPhoto) then

        write(mpiFilename,'(a, i4.4, a)') "photo_", grid%iDump,".grid"
        call writeAmrGrid(mpiFilename, .false., grid)
        write(mpiFilename,'(a, i4.4, a)') "photo", nIter,".vtk"!

!     if(hydrodynamics) then
        call writeVtkFile(grid, mpiFilename, &
             valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
             "OI           ","HeI          ","HeII         ", "OI           ", "OII          ", &
             "OIII         ","dust1        ","dust2        " /))
     end if
!
!     if(hydrodynamics) then
!        call writeVtkFile(grid, mpiFilename, &
!             valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!             "hydrovelocity","sourceCont   ","pressure     ", &
!             "crossings    "/))!
!!
 !    else
 !       call writeVtkFile(grid, mpiFilename, &
 !            valueTypeString=(/"HI           " , "temperature  ", &
 !            "sourceCont   "/))
!    else
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"HI           " , "temperature  ", &
!            "sourceCont   "/))
 !    end if

     if(myRankGlobal /= 0) then
!        iUnrefine = iUnrefine + 1
!        if(doUnRefine) then
!           if (iUnrefine == 5) then
!              if (myrank == 1)call tune(6, "Unrefine grid")
!               nUnrefine = 0                                                                                                            
!              call unrefineCells(grid%octreeRoot, grid, nUnrefine, 5.d-3)
!              !          write(*,*) "Unrefined ", nUnrefine, " cells"                                                   
!              if (myrank == 1)call tune(6, "Unrefine grid")
!              iUnrefine = 0
!           endif
!        end if
!        
!!!        call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
        if(doPhotoRefine) then!
           dprCounter = dprCounter + 1
           if(dprCounter == 4) then
              call refineGridGeneric(grid, amrtolerance, evenuparray) 
              call evenUpGridMPI(grid, .true., dorefine, evenuparray)
              dprCounter = 0
           end if
!           call writeInfo("Done the even up part", TRIVIAL)
!           
           !       if (doSelfGrav) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)                                    
           
!           if (myrank == 1) call tune(6, "Loop refine")
 !          
  !         call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
   !        
    !       if (myrank == 1) call tune(6, "Loop refine")
!           !                                                                                                         
        end if
     end if
!     
    

     call torus_mpi_barrier
     call MPI_BUFFER_DETACH(buffer,bufferSize, ierr)
     if(optimizeStack) then
        deallocate(photonPacketStack)     
        deallocate(currentStack)     
        deallocate(toSendStack)     
        deallocate(buffer)
     endif
  enddo


 call MPI_BARRIER(MPI_COMM_WORLD, ierr)

 call MPI_TYPE_FREE(mpi_vector, ierr)
 call MPI_TYPE_FREE(mpi_photon_stack, ierr)

 deallocate(nSaved)
 deallocate(nEscapedArray)
 if(allocated(buffer)) then
    deallocate(buffer)
 end if
end subroutine photoIonizationloopAMR

SUBROUTINE toNextEventPhoto(grid, rVec, uHat,  escaped,  thisFreq, nLambda, lamArray, photonPacketWeight, epsOverDeltaT, &
     nfreq, freq, tPhoton, tLimit, crossedMPIboundary, newThread, sourcePhoton, crossedPeriodic)
  use inputs_mod, only : periodicX, periodicY, periodicZ
  use mpi

   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec,uHat, octVec,thisOctVec, tvec, oldRvec, upperBound, lowerBound, iniVec
   type(OCTAL), pointer :: thisOctal, tempOctal
   type(OCTAL),pointer :: oldOctal
   type(OCTAL),pointer :: endOctal
   real(double) :: tPhoton, tLimit, epsOverDeltaT
   real(double) :: photonMomentum
   integer, intent(out) :: newThread
   integer :: endSubcell
   real(double) :: photonPacketWeight
   integer :: nFreq
   real(double) :: freq(:)
   integer :: subcell, tempSubcell
   real(oct) :: tval, tau, r
   real :: lamArray(:)
   integer :: nLambda
   logical :: stillinGrid
   logical :: escaped
   real(double) :: kappaScaDb, kappaAbsDb
   real(oct) :: thisTau
   real(oct) :: thisFreq
   real(oct) :: thisLam
   integer :: iLam
   logical ::inFlow
   real :: diffusionZoneTemp
!   real :: lambda
!   integer :: ilambda
   logical, intent(out) :: crossedMPIboundary
   type(OCTAL), pointer :: nextOctal
   integer :: nextSubcell
   logical :: ok, outofTime, neverInteracted
   logical :: sourcePhoton, crossedPeriodic


    stillinGrid = .true.
    escaped = .false.
    
    crossedMPIboundary = .false.
    outOfTime = .false.

    neverInteracted = .true.

    Photonmomentum = epsOverDeltaT / cSpeed
    thisLam = (cSpeed / thisFreq) * 1.e8

    call locate(lamArray, nLambda, real(thisLam), iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam,thisLam
    endif

! select an initial random tau and find distance to next cell

    call findSubcellTD(rVec, grid%octreeRoot,thisOctal, subcell)

    call randomNumberGenerator(getDouble=r)
    tau = -log(1.0-r)

!    write(*,*) "calling distancetocellboundary with "
!    write(*,*) "rVec ",rVec
!    write(*,*) "uHat ",uHat
    call distanceToCellBoundary(grid, rVec, uHat, tval, thisOctal, subcell)
    
    octVec = rVec
    thisOctVec = rVec

    call locate(lamArray, nLambda, real(thisLam), iLam)

    call amrGridValues(grid%octreeRoot, octVec,  iLambda=iLam, lambda=real(thisLam), startOctal=thisOctal, &
         actualSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
         grid=grid, inFlow=inFlow)
    oldOctal => thisOctal

    if (inFlow) then
       thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
    else
       thisTau = 1.0e-28
    end if

! if tau > thisTau then the photon has traversed a cell with no interactions

    do while(stillinGrid .and. (tau > thisTau) .and. .not. outOfTime) 

! add on the distance to the next cell

       oldrVec = rVec
       rVec = rVec + (tVal+1.d-3*grid%halfSmallestSubcell) * uHat ! rvec is now at face of thisOctal, subcell

! at this point the photon packet takes momentum from the cell

       thisOctal%radiationMomentum(subcell) = thisOctal%radiationMomentum(subcell) - uHat * photonMomentum
!       if (myrankGlobal == 1) write(*,*) "mom sub ",thisOctal%radiationMomentum(subcell)
       tPhoton = tPhoton + (tVal * 1.d10) / cSpeed
       if (tPhoton > tLimit) then
          escaped = .true.
          outOfTime = .true.
       endif

! check whether the photon has escaped from the grid
       if (inOctal(grid%octreeRoot,rVec)) then
          nextOctal => thisOctal
          nextSubcell = subcell
          call findSubcellLocal(rVec, nextOctal, nextSubcell)
       else
          if(.not. crossedPeriodic .and. neverInteracted) then
             !Give photon packets a second chance if they haven't interacted
             stillingrid = .false.
             crossedPeriodic = .true.
             crossedMPIBoundary = .true.
             
             !Calculate new position vector
             upperBound%x = grid%octreeRoot%centre%x + grid%octreeRoot%subcellSize
             lowerBound%x = grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize
             upperBound%y = grid%octreeRoot%centre%y + grid%octreeRoot%subcellSize
             lowerBound%y = grid%octreeRoot%centre%y - grid%octreeRoot%subcellSize
             upperBound%z = grid%octreeRoot%centre%z + grid%octreeRoot%subcellSize
             lowerBound%z = grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize

             iniVec = rVec

             if(rVec%x < lowerBound%x .and. periodicX)   rVec%x = rVec%x + grid%octreeRoot%subcellSize*2.d0
             if(rVec%y < lowerBound%y .and. periodicY)   rVec%y = rVec%y + grid%octreeRoot%subcellSize*2.d0
             if(rVec%z < lowerBound%z .and. periodicZ)   rVec%z = rVec%z + grid%octreeRoot%subcellSize*2.d0
             if(rVec%x > upperBound%x .and. periodicX)   rVec%x = rVec%x - grid%octreeRoot%subcellSize*2.d0
             if(rVec%y > upperBound%y .and. PeriodicY)   rVec%y = rVec%y - grid%octreeRoot%subcellSize*2.d0
             if(rVec%z > upperBound%z .and. periodicZ)   rVec%z = rVec%z - grid%octreeRoot%subcellSize*2.d0

             if(rVec == iniVec) then
                !no transfer occurs
                crossedPeriodic = .false.
                crossedMPIBoundary = .false.
             end if

             !This needs to be better
             !In 1D the vector may have been altered so that the pp escapes in the y or z direction
             !(not sure if this is a bug in the code) 
             if(grid%octreeRoot%oneD) then
                if(rvec%x == iniVec%x) then
                   !Escape was due to y or z
                   crossedPeriodic = .false.
                   crossedMPIBoundary = .false.
                   rVec = iniVec
                end if
             elseif(grid%octreeRoot%twoD) then
                if(rvec%y /= iniVec%y) then
                   !escape was due to y
                   rVec = iniVec
                end if 
             end if

             if(crossedPeriodic .and. (inOctal(grid%octreeRoot,rVec))) then
                nextOctal => thisOctal
                nextSubcell = subcell
                call findSubcellLocal(rVec, nextOctal, nextSubcell)
                newThread = nextOctal%mpiThread(nextSubcell)
!                print *, "rVec1 ", iniVec
!                print *, "rVec2 ", rVec
             else
                crossedPeriodic = .false.
                crossedMPIBoundary = .false.
                rVec = iniVec
   
             end if



          else
             stillingrid = .false.
          end if
          escaped = .true.
       endif
          
       octVec = rVec
       
! check whether the photon has  moved across an MPI boundary
        
      if (stillInGrid) then

          if (.not.octalOnThread(nextOctal, nextsubcell, myRankGlobal)) then
             stillInGrid = .false.
             escaped = .true.
             crossedMPIboundary = .true.
             newThread = nextOctal%mpiThread(nextSubcell)
          endif
       endif

! update the distance grid

       if (.not.outOfTime) then
          call updateGrid(grid, thisOctal, subcell, thisFreq, tVal, photonPacketWeight, ilam, nfreq, freq, sourcePhoton)
       endif
         
       if (stillinGrid) then
          
          thisOctal => nextOctal
          subcell = nextSubcell

! here the photon has entered a new cell on the current thread and we can add the momentum

          thisOctal%radiationMomentum(subcell) = thisOctal%radiationMomentum(subcell) + uHat * photonMomentum
!          if (myrankGlobal == 1) write(*,*) "mom add ",thisOctal%radiationMomentum(subcell)

          if (thisOctal%diffusionApprox(subcell)) then
             call randomWalk(grid, thisOctal, subcell,  endOctal, endSubcell, diffusionZoneTemp, ok)
             rVec = subcellCentre(endOctal,endSubcell)
             octVec = rVec
             call amrGridValues(grid%octreeRoot, rVec, foundOctal=tempOctal, &
                  foundSubcell=tempsubcell)
             uHat = randomUnitVector()
          endif
       endif


! now if the photon is in the grid choose a new random tau

       if (stillingrid) then
          neverInteracted = .false.
          call randomNumberGenerator(getDouble=r)
          tau = -log(1.0-r)
          call locate(lamArray, nLambda, real(thisLam), iLam)
          call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,lambda=real(thisLam), startOctal=thisOctal, &
               actualSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
               grid=grid, inFlow=inFlow)
          oldOctal => thisOctal
          thisOctVec = octVec

! calculate the distance to the next cell

          call distanceToCellBoundary(grid, rVec, uHat, tval, thisOctal, subcell)

          octVec = rVec

! calculate the optical depth to the next cell boundary

          if (inFlow) then
             thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
          else
             thisTau = 1.0e-28
          end if

          if (tVal == 0.0d0) then
             escaped = .true.
             stillingrid = .false.
          endif
       endif ! still in grid

! if photon is still in grid and  tau > tau_to_the_next_cell then loop round again
! choosing a new tau
       
    enddo

! the photon may have escaped the grid...

       if (inOctal(grid%octreeRoot,rVec)) then
          nextOctal => thisOctal
          nextSubcell = subcell
          call findSubcellLocal(rVec, nextOctal, nextSubcell)
       endif

 ! if not the photon must interact in this cell
       

    if (.not.escaped) then

       octVec = rVec
!       if (.not.inOctal(grid%octreeRoot, octVec)) then
!          write(*,*) "Error:: Photon location is out of boundaries, but its status is not ESCAPED."
!          write(*,*) "        .... lucy_mod::toNextEventAMR]"
!          write(*,*) "octVec-centre = ",octVec-grid%octreeRoot%centre
!          write(*,*) "cell size = ",grid%octreeRoot%subcellsize
!          stop
!       endif
       if (dble(tau)/thisTau .gt. 1.d0) then
          write(*,*) "tau prob",tau,thisTau
       endif

       call amrGridValues(grid%octreeRoot, octVec, startOctal=thisOctal, actualSubcell = subcell, &
            iLambda=iLam,lambda=real(thisLam),&
            kappaAbs=kappaAbsdb,kappaSca=kappaScadb, grid=grid, inFlow=inFlow)
       thisOctVec = octVec

       if (thisOctal%diffusionApprox(subcell)) then
          write(*,*) "Photon in diff zone but not escaped"
       endif

       if (.not.inFlow) kappaAbsdb =0.0d0


! update the distance grid

       if (thisTau > 0.d0) then

!          lambda = cSpeed*1.e8/thisFreq
!          call locate(lamArray, nLambda, lambda, ilambda)
          if (.not.outOfTime) then 
             call updateGrid(grid, thisOctal, subcell, thisFreq, &
                  dble(tval)*dble(tau)/thisTau, photonPacketWeight, ilam, nfreq, freq, sourcePhoton)
          endif

          oldOctal => thisOctal
          
       endif

       if (tau > thisTau) then
          write(*,*) "tau > thistau", tau, thisTau
          stop
       endif

! move the requisite distance within the cell and return

       tVec = rVec
       rVec = rVec + (dble(tVal)*dble(tau)/thisTau) * uHat
       tPhoton = tPhoton +  (dble(tVal*1.d10)*dble(tau)/thisTau) / cSpeed
       if (tPhoton > tLimit) then
          escaped = .true.
          outOfTime = .true.
        
       endif

       if (.not.inOctal(grid%octreeRoot, rVec)) then  ! this is only needed due to floating point boundary issues
          escaped = .true.
          goto 666
       endif

       octVec = rVec
       if (thisOctal%diffusionApprox(subcell)) then
          call randomWalk(grid, thisOctal, subcell,  endOctal, endSubcell, diffusionZoneTemp, ok)
          rVec = subcellCentre(endOctal,endSubcell)
          octVec = rVec
          call amrGridValues(grid%octreeRoot, rVec, foundOctal=tempOctal, &
            foundSubcell=tempsubcell)
          uHat = randomUnitVector()
       end if
    endif

666 continue

 end subroutine toNextEventPhoto

! recursive subroutine perturbIfront(thisOctal, grid)
!   type(octal) :: thisOctal
!   type(octal) :: child
!   type(gridtype) :: grid
!   integer :: i, subcell!!
!
!    do subcell = 1, thisOctal%maxChildren
!       if (thisOctal%hasChild(subcell)) then
!          ! find the child
!          do i = 1, thisOctal%nChildren, 1
!             if (thisOctal%indexChild(i) == subcell) then
!                child => thisOctal%child(i)
!                call perturbIfront(child, grid)
!                exit
!             end if
!          end do
!       else
!          locator = subcellcentre(thisoctal, subcell) - direction * (thisoctal%subcellsize/2.d0+0.01d0*grid%halfsmallestsubcell)
!          neighbouroctal => thisoctal
!          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
!          call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, reversedirection, q, rho, rhoe, &
!               rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz)!
!!
 !      end if
 !   end do
!
! end subroutine perturbIfront

  !Clear the tracers for photon contributions between iterations 
  !Removes contamination from noisy, early iterations
  recursive subroutine clearContributions(thisOctal)
    TYPE(OCTAL),pointer :: thisOctal
    TYPE(OCTAL),pointer :: child
    integer :: i, subcell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call clearContributions(child)
                exit
             end if
          end do
       else
          thisOctal%normSourceContribution(subcell, :) = 0.d0
          thisOctal%sourceContribution(subcell, :) = 0.d0
          thisOctal%diffuseContribution(subcell, :) = 0.d0
       end if
    end do
  end subroutine clearContributions

recursive subroutine advancedCheckForPhotoLoop(grid, thisOctal, photoLoop, dt, timeSinceLastRecomb)
    use constants_mod, only : cspeed
    use mpi
    TYPE(GRIDTYPE) :: grid
    TYPE(OCTAL),pointer :: thisOctal
    TYPE(OCTAL),pointer :: child
    integer :: i, subcell
    logical :: photoLoop
    real(double) :: chargeExchangeRecombination
    real(double) :: dt, timeSinceLastRecomb
    real(double) :: recombRateArray(grid%nIon), bigRecombRate(grid%nIon), recombTime(grid%nIon)
    integer :: k 
    
    photoLoop = .false.

    recombRateArray = 0.d0
    bigRecombRate = 0.d0
    recombTime = 0.d0
    k = 0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
                  if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call advancedCheckForPhotoLoop(grid, child, photoLoop, dt, timeSinceLastRecomb)
                exit
             end if
          end do
       else
          
          !Thaw -at present only considering H at ionization front
          if(thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon)) > 0.1d0 .and. &
               thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon)) < 0.90d0) then
 
                call getChargeExchangeRecomb(grid%ion(2), thisOctal%temperature(subcell), &
                     thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
                     chargeExchangeRecombination)

                recombRateArray(1) = recombRate(grid%ion(1), thisOctal%temperature(subcell))*thisOctal%ne(subcell) &
                     + chargeExchangeRecombination

             !Convert to recomb timescales
                recombTime(1) = 1.d0 / recombRateArray(1)

             
             !If the recomb timescale is less than the dynamical step, then a photoion loop is necessary
                if(recombTime(1) < dt) then
                   photoLoop = .true.
                   exit
                else if( timeSinceLastRecomb /= 0.d0) then
                   if(recombTime(1) < (timeSinceLastRecomb+dt)) then
                      photoLoop = .true.
                      timeSinceLastRecomb = 0.d0
                      exit
                   end if
                else
                   photoLoop = .false.
                end if
             end if
       end if
       if(photoLoop) exit
    end do

end subroutine advancedCheckForPhotoLoop

recursive subroutine checkForPhotoLoop(grid, thisOctal, photoLoop, dt)
  use constants_mod, only : cspeed
  use inputs_mod, only: lambdaSmooth, gridDistanceScale
  use mpi
  TYPE(GRIDTYPE) :: grid
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  logical :: photoLoop
  real(double) :: thisTau, ksca, kabs
  real(double) :: velocity, dt
  integer, save :: iLambda

  photoLoop = .false.
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call checkForPhotoLoop(grid, child, photoLoop, dt)
              exit
           end if
        end do
     else
        call locate(grid%lamArray, grid%nLambda, lambdaSmooth, ilambda)         
        !For cells that are not essentially completely ionized
        call returnKappa(grid, thisOctal, subcell, ilambda=ilambda,&
               kappaSca=ksca, kappaAbs=kabs)
        thisTau =  thisOctal%subcellSize * (kabs + ksca)* gridDistanceScale
        
        velocity = ((thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/ &
               (thisOctal%rho(subcell)**2))**(0.5)
        
        if(thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon)) > 0.01d0 .and.  & 
             thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon)) < 0.95d0) then
             
           if(velocity*dt <= (thisTau) .or. thisOctal%subcellSize < thisTau) then
              photoLoop = .true.
              exit
             else
                photoLoop = .false.
             end if
          end if
       end if
    end do
    
  end subroutine checkForPhotoLoop

  RECURSIVE SUBROUTINE resizePhotoionCoeff(thisOctal, grid)

    type(GRIDTYPE) :: grid
    TYPE(OCTAL), POINTER  :: thisOctal 
    TYPE(OCTAL), POINTER  :: child
    INTEGER :: i, j, k, m
    real(double), allocatable :: temp(:,:)

    IF ( thisOctal%nChildren > 0 ) THEN
       ! call this subroutine recursively on each of its children
       DO i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          CALL resizePhotoionCoeff(child, grid)
       END DO
    endif
    j = SIZE(thisOctal%photoIonCoeff,1)
    k = SIZE(thisOctal%photoIonCoeff,2)
    allocate(temp(1:j,1:k))
    temp = thisOctal%photoIonCoeff
    deallocate(thisOctal%photoionCoeff)
    allocate(thisOctal%photoionCoeff(j,grid%nion))
    thisOctal%photoionCoeff = 0.d0
    do m = 1, j
       if (k > grid%nion) then
          thisOctal%photoionCoeff(m, 1:grid%nion) = temp(m,1:grid%nion)
       else
          thisOctal%photoionCoeff(m, 1:k) = temp(m,1:k)
       endif
    enddo
    deallocate(temp)
    !Allocate contributors
    allocate(thisOctal%sourceContribution(j,grid%nion))
    allocate(thisOctal%diffuseContribution(j,grid%nion))
    allocate(thisOctal%normSourceContribution(j,grid%nion))

    thisOctal%sourceContribution = 0.d0
    thisOctal%diffuseContribution = 0.d0
    thisOctal%normSourceContribution = 0.d0


    j = SIZE(thisOctal%ionFrac,1)
    k = SIZE(thisOctal%ionFrac,2)
    allocate(temp(1:j,1:k))
    temp = thisOctal%ionFrac
    deallocate(thisOctal%ionFrac)
    allocate(thisOctal%ionFrac(j,grid%nion))
    thisOctal%photoionCoeff = 0.d0
    do m = 1, j
       if (k > grid%nion) then
          thisOctal%ionFrac(m, 1:grid%nion) = temp(m,1:grid%nion)
       else
          thisOctal%ionFrac(m, 1:k) = temp(m,1:k)
       endif
    enddo
    deallocate(temp)
          
  END SUBROUTINE resizePhotoionCoeff

  function hasPhotoionAllocations(grid) result(found)
    type(GRIDTYPE) :: grid
    logical :: found

    found = .true.
    call checkPhotoIonAllocationsSub(grid%octreeRoot, found)
  end function hasPhotoionAllocations

  recursive subroutine checkPhotoionAllocationsSub(thisOctal, found)
    integer :: subcell, i
    type(OCTAL), pointer :: thisOctal, child
    logical :: found
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call checkPhotoionAllocationsSub(child, found)
                exit
             end if
          end do
       else
          if (.not.associated(thisOctal%photoionCoeff)) found = .false.
       endif
    enddo
  end subroutine CheckPhotoionAllocationsSub


  RECURSIVE subroutine allocatePhotoionAttributes(thisOctal, grid)
    use inputs_mod, only : nDustType
    type(GRIDTYPE) :: grid
    TYPE(OCTAL), POINTER  :: thisOctal 
    TYPE(OCTAL), POINTER  :: child
    INTEGER :: i

    IF ( thisOctal%nChildren > 0 ) THEN
       ! call this subroutine recursively on each of its children
       DO i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          CALL allocatePhotoIonAttributes(child, grid)
       END DO
    endif

      call allocateAttribute(thisOctal%oldFrac, thisOctal%maxChildren)
      thisOctal%oldFrac = 1.d-30
      call allocateAttribute(thisOctal%dustType, thisOctal%maxChildren)
      thisOctal%dustType = 1
      call allocateAttribute(thisOctal%dustTypeFraction, thisOctal%maxChildren, nDustType)
      thisOctal%dustTypeFraction = 0.d0
      thisOctal%dustTypeFraction(:,:) = 1.d0

      call allocateAttribute(thisOctal%diffusionApprox, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%changed, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%nDiffusion, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%eDens, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%diffusionCoeff, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%oldeDens, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%nDirectPhotons, thisOctal%maxChildren)
       
      call allocateAttribute(thisOctal%underSampled, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%oldTemperature, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%kappaRoss, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%distanceGrid, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%nCrossings, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%biasLine3D, thisOctal%maxChildren)

      call allocateAttribute(thisOctal%etaCont, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%nh, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%ne, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%nhi, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%nhei, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%nhii, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%biasCont3D, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)

      call allocateAttribute(thisOctal%HHeating, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%HeHeating, thisOctal%maxChildren)
      call allocateAttribute(thisOctal%radiationMomentum,thisOctal%maxChildren)

      call allocateAttribute(thisOctal%ionFrac, thisOctal%maxchildren, grid%nIon)
      call allocateAttribute(thisOctal%photoionCoeff, thisOctal%maxchildren, grid%nIon)
      call allocateAttribute(thisOctal%sourceContribution, thisOctal%maxchildren, grid%nIon)
      call allocateAttribute(thisOctal%diffuseContribution, thisOctal%maxchildren, grid%nIon)
      call allocateAttribute(thisOctal%normSourceContribution, thisOctal%maxchildren, grid%nIon)
      thisOctal%inflow = .true.

    end SUBROUTINE allocatePhotoionAttributes



  recursive subroutine zeroDistanceGrid(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroDistanceGrid(child)
                exit
             end if
          end do
       else
          thisOctal%distanceGrid(subcell) = 0.d0
          thisOctal%nCrossings(subcell) = 0
          thisOctal%undersampled(subcell) = .false.
          thisOCtal%nDiffusion(subcell) = 0.
          thisOctal%Hheating(subcell) = 0.d0
          thisOctal%Heheating(subcell) = 0.d0
          thisOctal%photoIonCoeff(subcell,:) = 0.d0
          thisOctal%distanceGrid(subcell) = 0.d0
          thisOctal%radiationMomentum(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
       endif
    enddo
  end subroutine zeroDistanceGrid

  recursive subroutine testIonFront(thisOctal, currentTime)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: rVec
  real(double) :: currentTime
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call testIonFront(child, currentTime)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          if (rVec%x < -400.d6) then
             thisOctal%ionFrac(subcell,1) = 0.d0
             thisOctal%ionFrac(subcell,2) = 1.d0
          else
             thisOctal%ionFrac(subcell,1) = 1.d0
             thisOctal%ionFrac(subcell,2) = 0.d0
          endif
       endif
    enddo
  end subroutine testIonFront

  recursive subroutine ionizeGrid(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call ionizeGrid(child)
                exit
             end if
          end do
       else
          thisOctal%ionFrac(subcell,1) = 1.e-10
          thisOctal%ionFrac(subcell,2) = 1.
          thisOctal%ne(subcell) = thisOctal%nh(subcell)
          if (SIZE(thisOctal%ionFrac,2)>2) then
             thisOctal%ionFrac(subcell,3) = 1.e-10
             thisOctal%ionFrac(subcell,4) = 1.       
          endif
       endif
    enddo
  end subroutine ionizeGrid

  recursive subroutine neutralGrid(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call neutralGrid(child)
                exit
             end if
          end do
       else
          thisOctal%ionFrac(subcell,:) = 1.e-10
          thisOctal%ionFrac(subcell,1) = 1.
          thisOctal%ne(subcell) = 1.d-10
          thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
          if (SIZE(thisOctal%ionFrac,2)>2) then
             thisOctal%ionFrac(subcell,3) = 1.
             thisOctal%ionFrac(subcell,4) = 1.e-10       
          endif
       endif
    enddo
  end subroutine neutralGrid

  recursive subroutine resetNh(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call resetNh(child)
                exit
             end if
          end do
       else
          thisOctal%nh(subcell) = thisOctal%rho(subcell)/ mHydrogen
       endif
    enddo
  end subroutine resetNh

  recursive subroutine checkfor2d(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call checkfor2d(child)
                exit
             end if
          end do
       else
          if (thisOctal%twod) write(*,*) "2d cell found!"
       endif
    enddo
  end subroutine checkfor2d

  recursive subroutine testforZero(thisOctal, n, m)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i, n, m
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call testforzero(child, n, m)
                exit
             end if
          end do
       else
          m = m + 1
          if (thisOctal%photoIonCoeff(subcell,1) == 0.d0) then
             n = n + 1
          endif

        
       endif
    enddo
  end subroutine testforZero

  recursive subroutine getHbetaluminosity(thisOctal, grid, luminosity)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(GRIDTYPE) :: grid
  real(double) :: luminosity, v, hbeta
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getHbetaLuminosity(child, grid, luminosity)
                exit
             end if
          end do
       else
          v = cellVolume(thisOctal, subcell)
          hbeta = (10.d0**(-0.870d0*log10(thisOctal%temperature(subcell))+3.57d0)) * &
               thisOctal%ne(subcell) * thisOctal%ionFrac(subcell, 2) * &
               thisOctal%nh(subcell)*grid%ion(1)%abundance*1.d-25
          luminosity = luminosity + hbeta * (v*1.d30)
       endif
    enddo
  end subroutine getHbetaluminosity

  subroutine calculateIonizationBalance(grid, thisOctal, epsOverDeltaT)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT
    integer :: subcell
    
    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then

          if (thisOctal%inflow(subcell)) then
             if (.not.thisOctal%undersampled(subcell)) then
                call solveIonizationBalance(grid, thisOctal, subcell, thisOctal%temperature(subcell), epsOverdeltaT)
               
             else
                thisOctal%ionFrac(subcell, 1) = 1.d0
                thisOctal%ionFrac(subcell, 2) = 1.d-30
                if(thisOctal%nCrossings(subcell) /= 0) then
                   !write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
                end if
             endif
          endif
       endif
    enddo
  end subroutine calculateIonizationBalance

  subroutine calculateIonizationBalancetimeDep(grid, thisOctal, epsOverDeltaT, DeltaT)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT, DeltaT
    integer :: subcell
    
    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then

          if (thisOctal%inflow(subcell)) then
             if (.not.thisOctal%undersampled(subcell)) then
                call solveIonizationBalanceTimeDep(grid, thisOctal, subcell, thisOctal%temperature(subcell), epsOverdeltaT, deltaT)
                
             else
                thisOctal%ionFrac(subcell, 1) = 1.d0
                thisOctal%ionFrac(subcell, 2) = 1.d-30
                if(thisOctal%nCrossings(subcell) /= 0) then
                 !write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
                end if
             endif
          endif
       endif
    enddo
  end subroutine calculateIonizationBalanceTimeDep

  subroutine quickThermalCalc(thisOctal)
    use mpi
    type(OCTAL), pointer :: thisOctal
    integer :: subcell


    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then
          if(octalOnThread(thisOCtal, subcell, myRankGlobal)) then
             thisOctal%temperature(subcell) = real(10.d0 + (10000.d0-10.d0) * thisOctal%ionFrac(subcell,2))
          end if
       endif
       
    enddo
  end subroutine quickThermalCalc

  subroutine calculateThermalBalance(grid, thisOctal, epsOverDeltaT)
    use mpi
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT
    real(double) :: totalHeating
    integer :: subcell
    logical :: converged, found
    real :: t1, t2, tm
!    real, parameter :: Tlow = 100.
    real, parameter :: Tlow = 3.
    real(double) :: y1, y2, ym, Hheating, Heheating, dustHeating
    real :: deltaT
    real :: underCorrection = 1.
    integer :: nIter
    

    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then

 !THAW - trying to root out mpi things
          if(octalOnThread(thisOctal, subcell, myRankGlobal)) then
          
             if (thisOctal%inflow(subcell)) then
                
!             write(*,*) thisOctal%nCrossings(subcell),thisOctal%undersampled(subcell)
                
                if (.not.thisOctal%undersampled(subcell)) then
                   
                   
                   call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
!                if (thisOctal%ionFrac(subcell,1) > 0.9d0) write(*,*) "total heating ",totalheating
                   if (totalHeating < 1.d-30) then
                      thisOctal%temperature(subcell) = tLow
                   else
                      nIter = 0
                      converged = .false.
                      
                      t1 = tlow
                      t2 = 30000.
!
                      
!                   if (thisOctal%dustTypeFraction(subcell,1) > 0.9) then
!                      t1 = tlow
!                      t2 = 2000.
!                   endif
!                   else
!                      t1 = 100.
!                      t2 = 30000.
!                      t1 = 5000.
!                      t1 = 20000.
!                   endif
!                   t1 = tLow
!                   t2 = 50000.
                      found = .true.
                      
                      if (found) then
                         y1 = (HHecooling(grid, thisOctal, subcell, t1) &
                         - totalHeating)
                         y2 = (HHecooling(grid, thisOctal, subcell, t2) &
                         - totalHeating)
                         if (y1*y2 > 0.d0) then
                            if (HHecooling(grid, thisOctal, subcell, t1) > totalHeating) then
                               tm = t1
                            else
                               tm  = t2
!                            write(*,*) "Cell set to hot. Cooling is ", &
!                                 HHecooling(grid, thisOctal, subcell, epsOverDeltat,t1), " heating is ", totalHeating
                            endif
                            converged = .true.
                         endif
                         if (totalHeating < 1.d-30) then
                            tm = t1
                            converged = .true.
                         endif
                         
! Find root of heating and cooling by bisection
                         
                         do while(.not.converged)
                            tm = 0.5*(t1+t2)
                            y1 = (HHecooling(grid, thisOctal, subcell, t1) &
                            - totalheating)
                            y2 = (HHecooling(grid, thisOctal, subcell, t2) &
                            - totalheating)
                            ym = (HHecooling(grid, thisOctal, subcell, tm) &
                            - totalheating)
                            
                            if (y1*ym < 0.d0) then
                               t1 = t1
                               t2 = tm
                            else if (y2*ym < 0.d0) then
                               t1 = tm
                               t2 = t2
                            else
                               converged = .true.
                               tm = 0.5*(t1+t2)
                               write(*,*) t1, t2, y1,y2,ym
                            endif
                            
                            if (abs((t1-t2)/t1) .le. 1.e-2) then
                               converged = .true.
                            endif
                            niter = niter + 1
!                         if (myrankGlobal == 1) write(*,*) niter,tm, abs(t1-t2)/t1
                         enddo
                      endif
                      deltaT = tm - thisOctal%temperature(subcell)
                      thisOctal%temperature(subcell) = &
                      max(thisOctal%temperature(subcell) + underCorrection * deltaT,tLow)
!                                write(*,*) thisOctal%temperature(subcell), niter
                   endif
                else
                !                write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
                endif
             endif
          end if
       endif
      enddo
  end subroutine calculateThermalBalance

  subroutine calculateEquilibriumTemperature(grid, thisOctal, subcell, epsOverDeltaT, Teq)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT, Teq
    real(double) :: totalHeating
    integer :: subcell
    logical :: converged, found
    real, parameter :: Tlow = 3.
    real(double) :: y1, y2, ym, Hheating, Heheating, dustHeating
    integer :: nIter
    real :: t1, t2, tm
                
    call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
    !                if (thisOctal%ionFrac(subcell,1) > 0.9d0) write(*,*) "total heating ",totalheating
    nIter = 0
    converged = .false.
    
    t1 = 100.
    t2 = 50000.
    !    t1 = 3.
    !    t2 = 50000.
    found = .true.
       
    if (found) then
       y1 = (HHecooling(grid, thisOctal, subcell, t1) &
            - totalHeating)
       y2 = (HHecooling(grid, thisOctal, subcell, t2) &
            - totalHeating)
       if (y1*y2 > 0.d0) then
          if (HHecooling(grid, thisOctal, subcell, t1) > totalHeating) then
             tm = t1
          else
             tm  = t2
             !                            write(*,*) "Cell set to hot. Cooling is ", &
             !                                 HHecooling(grid, thisOctal, subcell, epsOverDeltat,t1), " heating is ", totalHeating
          endif
          converged = .true.
       endif
       if (totalHeating < 1.d-30) then
          tm = t1
          converged = .true.
       endif
       
       ! Find root of heating and cooling by bisection
       
       do while(.not.converged)
          tm = 0.5*(t1+t2)
          y1 = (HHecooling(grid, thisOctal, subcell, t1) &
               - totalheating)
          y2 = (HHecooling(grid, thisOctal, subcell, t2) &
               - totalheating)
          ym = (HHecooling(grid, thisOctal, subcell, tm) &
               - totalheating)
          
          if (y1*ym < 0.d0) then
             t1 = t1
             t2 = tm
          else if (y2*ym < 0.d0) then
             t1 = tm
             t2 = t2
          else
             converged = .true.
             tm = 0.5*(t1+t2)
             write(*,*) t1, t2, y1,y2,ym
          endif
          
          if (abs((t1-t2)/t1) .le. 1.e-2) then
             converged = .true.
          endif
             niter = niter + 1
             !                         if (myrankGlobal == 1) write(*,*) niter,tm, abs(t1-t2)/t1
          enddo
       endif
       Teq = tm
     end subroutine calculateEquilibriumTemperature
     
     subroutine calculateThermalBalanceTimeDep(grid, thisOctal, epsOverDeltaT, deltaTime)
       type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT, deltaTime, thermalTime, Teq
    real(double) :: totalHeating, totalCooling, Hheating, heheating,dustheating, mu, energy, smallDeltaT
    integer :: subcell, nTimes, i
    real, parameter :: Tlow = 100.
    

    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then
          
          if (thisOctal%inflow(subcell)) then

             if (.not.thisOctal%undersampled(subcell)) then
                
                call calculateEquilibriumTemperature(grid, thisOctal, subcell, epsOverDeltaT, Teq)

                call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)

                mu = returnMu(thisOctal, subcell, grid%ion, grid%nion)
                energy = 1.5d0*(thisOctal%rho(subcell)/(mu*mhydrogen))*kerg*thisOctal%temperature(subcell)
                totalCooling = HHecooling(grid, thisOctal, subcell, thisOctal%temperature(subcell)) 

                thermalTime = (energy / abs(totalCooling))

                if (deltaTime > thermalTime) then
                   thisOctal%temperature(subcell) = real(Teq)
                else

                   nTimes = 10
                
                   smallDeltaT = deltaTime /dble(nTimes)
                   
                   do i = 1, nTimes
                      energy = energy + (totalHeating - totalcooling)*smallDeltaT
                      thisOctal%temperature(subcell) = real((2.d0/3.d0)*energy*(mu*mhydrogen)/(thisOctal%rho(subcell)*kerg))
                      totalCooling = HHecooling(grid, thisOctal, subcell, thisOctal%temperature(subcell)) 
                   enddo
                endif

             endif
          endif
       endif
    enddo
  end subroutine calculateThermalBalanceTimeDep
  
  function HHeCooling(grid, thisOctal, subcell, temperature, debug) result (coolingRate)
    use inputs_mod, only : dustOnly, hOnly
    type(OCTAL),pointer :: thisOctal
    integer :: subcell
    type(GRIDTYPE) :: grid
    real(double) :: nHii, nHeii, ne, nh
    real :: temperature
    logical, optional :: debug
    real(double) :: coolingRate, crate, dustCooling
    real(double) :: gff
    real :: rootTbetaH(31) = (/ 8.287e-11, 7.821e-11, 7.356e-11, 6.982e-11, 6.430e-11, 5.971e-11, 5.515e-11, 5.062e-11, 4.614e-11, &
                               4.170e-11, 3.734e-11, 3.306e-11, 2.888e-11, 2.484e-11, 2.098e-11, 1.736e-11, 1.402e-11, 1.103e-11, &
                               8.442e-12, 6.279e-12, 4.539e-12, 3.192e-12, 2.185e-12, 1.458e-12, 9.484e-13, 6.023e-13, 3.738e-13, &
                               2.268e-13, 1.348e-13, 7.859e-14, 4.499e-14 /)

    real :: rootTbetaHe(18) = (/ 8.347e-11, 7.889e-11, 7.430e-11, 6.971e-11, 6.512e-11, 6.056e-11, 5.603e-11, 5.154e-11, 4.710e-11,&
                               4.274e-11, 3.847e-11, 3.431e-11, 3.031e-11, 2.650e-11, 2.291e-11, 1.960e-11, 1.660e-11, 1.394e-11 /)

    real :: logT(31) = (/ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, &
                      5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0 /)
    real(double) :: fac, thisRootTbetaH, betaH, betaHe, thisRootTbetaHe
    real :: thisLogT
    integer :: n
    real(double) :: log10te,betarec, coolrec, betaff, coolff
    real :: ch12, ch13, ex12, ex13, th12, th13, coolcoll, te4, teused
    real(double) :: becool
    real :: kappap
    real, parameter                :: hcRyd = &    ! constant: h*c*Ryd (Ryd at inf used) [erg]
         & 2.1799153e-11


    coolingRate = 0.d0

    if (.not.dustOnly) then
                
       nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
       if(.not. hOnly) then
          nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
       else
          nHeii = 0.d0
       end if
       nh = thisOctal%nh(subcell)
       ne = thisOctal%ne(subcell)
       
       becool = 0.
       
       gff = 1.1d0 + 0.34d0*exp(-(5.5d0 - log10(temperature))**2 / 3.d0)  ! Kenny equation 23
       
       coolingRate = 1.42d-27 * (nHii+nHeii) * ne * gff * sqrt(temperature)  ! Kenny equation 22 (free-free)
       
       thisLogT = log10(temperature)
       log10te = thisLogT
       teused = temperature
       
       ! cooling of gas due to FF radiation from H+ 
       ! fits to Hummer, MNRAS 268(1994) 109, Table 1. or  least square fitting to m=4
       
       betaFF = 1.0108464E-11 + 9.7930778E-13*log10Te - &
            & 6.6433144E-13*log10Te*log10Te + 2.4793747E-13*log10Te*log10Te*log10Te -&
            & 2.3938215E-14*log10Te*log10Te*log10Te*log10Te
       
       coolFF = Nhii*Ne*betaFF*kerg*TeUsed/sqrt(TeUsed)
       
       becool = becool + coolfF
       
       ! cooling of gas due to recombination of H+
       ! fits to Hummer, MNRAS 268(1994) 109, Table 1.  
       ! least square fitting to m=4
       betaRec = 9.4255985E-11 -4.04794384E-12*log10Te &
            & -1.0055237E-11*log10Te*log10Te +  1.99266862E-12*log10Te*log10Te*log10Te&
            & -1.06681387E-13*log10Te*log10Te*log10Te*log10Te

       
       coolRec = Nhii*Ne*betaRec*kerg*TeUsed/sqrt(TeUsed)
       
       becool = becool+ coolrec
              
       ! collisional excitation of Hydrogen
       ! Mathis, Ly alpha, beta
       ch12 = 2.47e-8
       ch13 = 1.32e-8
       ex12 = -0.228
       ex13 = -0.460
       th12 = 118338.
       th13 = 140252.
       te4 = temperature / 1.e4
       coolColl = 0.
       
       if (TeUsed > 5000.) then 
          coolColl = real( &
               (ch12*exp(-th12/TeUsed)*Te4**ex12 + &
               & ch13*exp(-th13/TeUsed)*Te4**ex13) * &
               & hcRyd*(nh-nhii)*Ne )
       else
          coolColl = 0.
       end if
       
       coolingrate = coolingrate + coolcoll
       
       becool = becool +coolcoll
       
       if (PRESENT(debug)) then
          if (debug) then
             if (becool /= 0.) then
                write(*,*) coolff/becool, coolrec/becool, coolcoll/becool
             endif
          endif
       endif
       
       if (coolingRate < 0.) then
          write(*,*) "negative ff cooling",nhii,nheii,ne,gff,sqrt(temperature)
       endif

       if(coolingRate < 0.) then
          print *, "coolingRate ", coolingRate
       end if

       call locate(logT, 31, thisLogT, n)
       fac = (thisLogT - logT(n))/(logT(n+1)-logT(n))
       
       thisRootTbetaH = rootTbetaH(n) + (rootTbetaH(n+1)-rootTbetaH(n)) * fac

       if(thisRootTbetaH < 0.d0) then
          print *, " "
          print *, "n ", n
          print *, "rootTbetaH(n) ", rootTbetaH(n)
          print *, "rootTbetaH(n+1) ", rootTbetaH(n+1)
          print *, "fac ", fac
          print *, "thisLogT ", thisLogT
          print *, "logT(n) ", logT(n)
          print *, "logT(n+1) ", logT(n+1)
          print *, " "
          print *, "thisRootBetaH < 0.d0"
       end if

       if (thisLogT < logT(18)) then
          thisRootTbetaHe = rootTbetaHe(n) + (rootTbetaHe(n+1)-rootTbetaHe(n)) * fac
       else
          thisRootTbetaHe = 0.d0
       endif
       
       betaH = thisrootTbetaH / sqrt(temperature)
       betaHe = thisrootTbetaHe / sqrt(temperature)

       coolingRate = coolingRate +  ne * nhii * kerg * temperature * betaH
              
       if (ne * nhii * kerg * temperature * betaH < 0.) then
          write(*,*) "negative H cooling",ne,nhii,kerg,temperature,betah
       endif
       
       if(.not. hOnly) then
          coolingRate = coolingRate + ne * nheii * kerg * temperature * betaHe
       end  if

       if (ne * nheii * kerg * temperature * betaHe < 0.) then
          write(*,*) "negative He cooling",ne,nheii,kerg,temperature,betaHe
       endif
              
       call metalcoolingRate(grid%ion, grid%nIon, thisOctal, subcell, nh, ne, temperature, crate)
       if (crate < 0.) then
          write(*,*) "total negative metal cooling",crate
       endif
!           write(*,*) log(coolingRate),log(crate),coolingRate/(coolingrate+crate)
       coolingRate = coolingRate + crate
 
    endif

    call returnKappa(grid, thisOctal, subcell, kappap=kappap)
       
    dustCooling = fourPi * kappaP * (stefanBoltz/pi) * temperature**4

    coolingRate = coolingRate + dustCooling

  end function HHeCooling

  subroutine updateGrid(grid, thisOctal, subcell, thisFreq, distance, photonPacketWeight, ilambda, nfreq, freq, sourcePhoton)
    use inputs_mod,only : dustOnly
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: nFreq, iFreq
    real(double) :: freq(:)
    integer :: subcell
    real(double) :: thisFreq, distance, kappaAbs,kappaAbsDust
    integer :: ilambda
    real(double) :: photonPacketWeight
    integer :: i 
    real(double) :: fac, xSec
    logical :: sourcePhoton

    
    kappaAbs = 0.d0; kappaAbsDust = 0.d0
    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1

    fac = distance * photonPacketWeight / (hCgs * thisFreq)

    call locate(freq, nFreq, thisFreq, iFreq)


!    e = hCgs * nuHydrogen * ergtoEv

    if (.not.dustOnly) then
       do i = 1, grid%nIon
          xSec = returnxSec(grid%ion(i), thisFreq, iFreq=iFreq)          
          
          if (xSec > 0.d0 .and. thisFreq > nuHydrogen) then
             thisOctal%photoIonCoeff(subcell,i) = thisOctal%photoIonCoeff(subcell,i) &
                  + fac * xSec
             if(i == 1) then
                if(sourcePhoton) then
                   thisOctal%sourceContribution(subcell, i) = thisOctal%sourceContribution(subcell, i)  &
                        + fac * xSec
                else
                   thisOctal%diffuseContribution(subcell, i) = thisOctal%diffuseContribution(subcell, i)  &
                        + fac * xSec
                end if
             end if
             
             if(thisOctal%SourceContribution(subcell,i) == 0.d0) then
                thisOctal%normSourceContribution(subcell,i) = 0.d0
             else
                thisOctal%normSourceContribution(subcell, i) = thisOctal%sourceContribution(subcell, i) / &
                     (thisOctal%sourceContribution(subcell, i) + thisOctal%diffuseContribution(subcell, i))
             end if
             
             !          thisOctal%photoIonCoeff(subcell,i) = thisOctal%photoIonCoeff(subcell,i) &
             !               + distance * dble(xsec) / (dble(hCgs) * thisFreq) * photonPacketWeight
          endif

       ! neutral h heating

          if ((grid%ion(i)%z == 1).and.(grid%ion(i)%n == 1)) then
             thisoctal%hheating(subcell) = thisoctal%hheating(subcell) &
                  + distance * xsec / (thisfreq * hcgs) &
                  * ((hcgs * thisfreq) - (hcgs * grid%ion(i)%nuthresh)) * photonpacketweight
          endif

       ! neutral he heating

          if ((grid%ion(i)%z == 2).and.(grid%ion(i)%n == 2)) then
             thisoctal%heheating(subcell) = thisoctal%heheating(subcell) &
                  + distance * xsec / (thisfreq * hcgs) &
                  * ((hcgs * thisfreq) - (hcgs * grid%ion(i)%nuthresh)) * photonpacketweight
          endif

       enddo
    endif

    call returnkappa(grid, thisoctal, subcell, ilambda=ilambda, kappaabsdust=kappaabsdust, kappaabs=kappaabs)

    thisoctal%distancegrid(subcell) = thisoctal%distancegrid(subcell) &
         + dble(distance) * dble(kappaabsdust) * photonPacketWeight

  end subroutine updategrid

subroutine solveIonizationBalance(grid, thisOctal, subcell, temperature, epsOverdeltaT)
  implicit none
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  real(double) :: epsOverDeltaT, V
  real :: temperature
  integer :: subcell
  integer :: i, k
  integer :: nIonizationStages
  integer :: iStart, iEnd, iIon
  real(double) :: chargeExchangeIonization, chargeExchangeRecombination
  real(double), allocatable :: xplus1overx(:)

  v = cellVolume(thisOctal, subcell)
  k = 1
  
  do while(k <= grid%nIon)

     iStart = k
     iEnd = k+1
     do while(grid%ion(istart)%z == grid%ion(iEnd)%z)
        iEnd = iEnd + 1
     enddo
     iEnd = iEnd - 1
     nIonizationStages = iEnd - iStart + 1
          
     allocate(xplus1overx(1:nIonizationStages-1))
     do i = 1, nIonizationStages-1
        iIon = iStart+i-1
        call getChargeExchangeRecomb(grid%ion(iion+1), temperature, &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             chargeExchangeRecombination)
        
        call getChargeExchangeIon(grid%ion(iion), temperature, &
             thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
             chargeExchangeIonization)
           
        xplus1overx(i) = ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization) / &
             max(1.d-50,(recombRate(grid%ion(iIon),temperature) * thisOctal%ne(subcell) + chargeExchangeRecombination))
!        if ((myRankGlobal==1).and.(grid%ion(iion)%species(1:2) =="H ")) &
!             write(*,*) i,xplus1overx(i), thisOctal%photoioncoeff(subcell,iion), &
!             thisOctal%ne(subcell), epsoverdeltat,v,chargeExchangeIonization, chargeExchangeRecombination
     enddo
     thisOctal%ionFrac(subcell, iStart:iEnd) = 1.d0
     do i = 1, nIonizationStages - 1
        iIon = iStart+i-1
        thisOctal%ionFrac(subcell,iIon+1) = thisOctal%ionFrac(subcell,iIon) * xplus1overx(i)
!        if ((myRankGlobal==1).and.grid%ion(iion)%species(1:2) =="H ") write(*,*) i,thisOctal%ionFrac(subcell,iIon)
     enddo
     
     if (SUM(thisOctal%ionFrac(subcell,iStart:iEnd)) /= 0.d0) then
        thisOctal%ionFrac(subcell,iStart:iEnd) = &
             max(1.d-50,thisOctal%ionFrac(subcell,iStart:iEnd))/SUM(thisOctal%ionFrac(subcell,iStart:iEnd))
     else
        thisOctal%ionFrac(subcell,iStart:iEnd) = 1.d-50
     endif

        deallocate(xplus1overx)

        k = iEnd + 1
     end do

     thisOctal%ne(subcell) = returnNe(thisOctal, subcell, grid%ion, grid%nion)

   end subroutine solveIonizationBalance

subroutine solveIonizationBalanceTimeDep(grid, thisOctal, subcell, temperature, epsOverdeltaT, deltaT)
  implicit none
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  real(double) :: epsOverDeltaT, V, deltaT
  real :: temperature
  integer :: subcell
  integer :: i, k
  integer :: nIonizationStages
  integer :: iStart, iEnd, iIon
  real(double) :: chargeExchangeIonization, chargeExchangeRecombination
  real(double), allocatable :: xplus1overx(:)
  real(double) :: photoIonRate, recombinationRate, recombTime
  real(double), allocatable :: ionFrac(:)
  v = cellVolume(thisOctal, subcell)
  k = 1
  allocate(ionFrac(1:SIZE(thisOctal%ionFrac,2)))
  
  do while(k <= grid%nIon)

     iStart = k
     iEnd = k+1
     do while(grid%ion(istart)%z == grid%ion(iEnd)%z)
        iEnd = iEnd + 1
     enddo
     iEnd = iEnd - 1
     nIonizationStages = iEnd - iStart + 1
          
     allocate(xplus1overx(1:nIonizationStages-1))
     do i = 1, nIonizationStages-1
        iIon = iStart+i-1
        call getChargeExchangeRecomb(grid%ion(iion+1), temperature, &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             chargeExchangeRecombination)
        
        call getChargeExchangeIon(grid%ion(iion), temperature, &
             thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
             chargeExchangeIonization)
        
        xplus1overx(i) = ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization) / &
             max(1.d-50,(recombRate(grid%ion(iIon),temperature) * thisOctal%ne(subcell) + chargeExchangeRecombination))
!        if ((myRankGlobal==1).and.(grid%ion(iion)%species(1:2) =="H ")) &
!             write(*,*) i,xplus1overx(i), thisOctal%photoioncoeff(subcell,iion), &
!             thisOctal%ne(subcell), epsoverdeltat,v,chargeExchangeIonization, chargeExchangeRecombination
     enddo
     ionFrac(iStart:iEnd) = 1.d0
     do i = 1, nIonizationStages - 1
        iIon = iStart+i-1
        ionFrac(iIon+1) = ionFrac(iIon) * xplus1overx(i)
        !if ((myRankGlobal==1).and.grid%ion(iion)%species(1:2) =="H ") write(*,*) i,thisOctal%ionFrac(subcell,iIon)
     enddo
     
     if (SUM(ionFrac(iStart:iEnd)) /= 0.d0) then
        ionFrac(iStart:iEnd) = &
             max(1.d-50,ionFrac(iStart:iEnd))/SUM(ionFrac(iStart:iEnd))
     else
        ionFrac(iStart:iEnd) = 1.d-50
     endif

        deallocate(xplus1overx)

        k = iEnd + 1
     end do

  k = 1
  
  do while(k <= grid%nIon)

     iStart = k
     iEnd = k+1
     do while(grid%ion(istart)%z == grid%ion(iEnd)%z)
        iEnd = iEnd + 1
     enddo
     iEnd = iEnd - 1
     nIonizationStages = iEnd - iStart + 1
          
     do i = 1, nIonizationStages-1
        iIon = iStart+i-1
        call getChargeExchangeRecomb(grid%ion(iion+1), temperature, &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             chargeExchangeRecombination)
        
        call getChargeExchangeIon(grid%ion(iion), temperature, &
             thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
             chargeExchangeIonization)
        
        photoIonRate = (epsOverDeltaT / (v*1.d30))*thisOctal%photoIonCoeff(subcell,iIon) + chargeExchangeIonization
        

        recombinationRate = recombRate(grid%ion(iion), temperature) * thisOctal%ne(subcell) + chargeExchangeRecombination

        recombTime = 1.d0/recombinationRate
        thisOctal%ionFrac(subcell,iion) = thisOctal%ionFrac(subcell,iion) + deltaT * (recombinationRate - photoIonRate)

     if (SUM(thisOctal%ionFrac(subcell,iStart:iEnd)) /= 0.d0) then
        thisOctal%ionFrac(subcell,iStart:iEnd) = &
             max(1.d-50,thisOctal%ionFrac(subcell,iStart:iEnd))/SUM(thisOctal%ionFrac(subcell,iStart:iEnd))
     else
        thisOctal%ionFrac(subcell,iStart:iEnd) = 1.d-50
     endif
     
        thisOctal%ionFrac(subcell, iIon) = min(max(thisOctal%ionFrac(subcell,iIon),ionFrac(iIon)),1.d0)

        if(thisOctal%ionFrac(subcell, iIon) <= 0.d0) then
           thisOctal%ionFrac(subcell, iIon) = 1.d-50
        end if

     enddo

     thisOctal%ionFrac(subcell, iEnd) = 1.d0 - SUM(thisOctal%ionFrac(subcell,iStart:iEnd-1))
     thisOctal%ionFrac(subcell, iend) = min(max(thisOctal%ionFrac(subcell,iend),ionFrac(iend)),1.d0)

     if(thisOctal%ionFrac(subcell, iend) <= 0.d0) then
        thisOctal%ionFrac(subcell, iend) = 1.d-50
     end if

!THAW - getting rid of -ve cooling rates
     if(thisOctal%ionFrac(subcell, istart) <= 0.d0) then
        thisOctal%ionFrac(subcell, istart) = 1.d-50
     end if

        k = iEnd + 1
     end do

     thisOctal%ne(subcell) = returnNe(thisOctal, subcell, grid%ion, grid%nion)
     deallocate(ionFrac)

   end subroutine solveIonizationBalanceTimeDep


function recombRate(thisIon, temperature) result (rate)
  type(IONTYPE) :: thisIon
  real :: temperature
  real(double) :: rate
  real(double) :: a, b, t0, t1, c, d, f, t, z

! recombinations INTO this species

  select case(thisIon%z)

     case(1)
        select case(thisIon%n)
           case(1) ! H I
              a = 7.982e-11
              b = 0.7480
              t0 = 3.148e0
              t1 = 7.036e5
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate = 0.
        end select
     case(2)
        select case(thisIon%n)
           case(2) ! He I 
              a = 3.294e-11
              b = 0.6910
              t0 = 1.544e1
              t1 = 3.676e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case(1) ! He II
              a = 1.891e-10
              b = 0.7524
              t0 = 9.370e00
              t1 = 2.774e6
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(6)
        select case(thisIon%n)
           case(6) ! C I 
              a = 0.0108
              b = -0.1075
              c = 0.2810
              d = -0.0193
              f = -0.1127
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 1
              a = 5.068
              b = -0.6192
              c = -0.0815
              d = 1.2910
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! C II
              a = 1.8267
              b = 4.1012
              c = 4.8443
              d = 0.2261
              f = 0.5960
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 5.434
              b = -0.6116
              c = 0.0694
              d = 0.7866
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! C III
              a = 2.3196
              b = 10.7328
              c = 6.8830
              d = -0.1824
              f = 0.4101
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.742
              b = -0.6167
              c = 0.2960
              d = 0.6167
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(3) ! C IV
              a = 8.540e-11
              b = 0.5247
              t0 = 5.014e2
              t1 = 1.479e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case(2) ! C V
              a = 2.765e-10
              b = 0.6858
              t0 = 1.535e2
              t1 = 2.556e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(7)
        select case(thisIon%n)
           case(7) ! N I 
              a = 0.0
              b = 0.6310
              c = 0.1990
              d = -0.0197
              f = 0.4398
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 1
              a = 3.874
              b = -0.6487
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! N II
              a = 0.0320
              b = -0.6624
              c = 4.3191
              d = 0.0003
              f = 0.5946
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 4.974
              b = -0.6209
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! N III
              a = -0.8806
              b = 11.2406
              c = 30.7066
              d = -1.1721
              f = 0.6127
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.750
              b = -0.5942
              c = 0.8452
              d = 2.8450
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! N IV
              a = 0.4134
              b = -4.6319
              c = 25.9172
              d = -2.2290
              f = 0.2360
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 4
              a = 4.626
              b = -0.9521
              c = 0.4729
              d = -0.4477
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(3) ! N V
              a = 1.169e-10
              b = 0.5470
              t0 = 6.793e2
              t1 = 1.650e7
              rate = vernerFerland(a, dble(temperature), t0, t1, b)
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select
     case(8)
        select case(thisIon%n)
           case(8) ! O I 
              t = temperature / 10000.
              if (t < 2.) then
                 a = -0.0001
                 b = 0.0001
                 c = 0.0956
                 d = 0.0193
                 f = 0.4106
              else
                 a = 0.3715
                 b = -0.0293
                 c = -0.0597
                 d = 0.0678
                 f = 0.7993
              endif
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 1
              a = 3.201
              b = -0.6880
              c = -0.0174
              d = 1.7070
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(7) ! O II
              a = -0.0036
              b = 0.7519
              c = 1.5252
              d =-0.0838
              f = 0.2769
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 4.092
              b = -0.6413
              c = 0.
              d = 1.
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! O III
              a = 0.0
              b = 21.8790
              c = 16.2730
              d = -0.7020
              f = 1.1899
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 4.890
              b = -0.6213
              c = 0.0184
              d = 1.5550
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! O IV
              t = temperature / 10000.
              if (t < 2.) then
                 a = -0.3648
                 b =  7.2698
                 c =  17.2187
                 d =  9.8335
                 f = -0.0166
              else
                 a = -2.5053
                 b = 3.4903
                 c = 67.4128
                 d = -3.4450
                 f = 0.8501
              endif
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 4
              a = 14.665
              b = -0.5140
              c = 2.7300
              d = 0.2328
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(4) ! O V
              a = -2.8425
              b = 0.2283
              c = 40.4072
              d = -3.4956
              f = 1.7558
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 5
              a = 4.626
              b = -0.9521
              c = 0.4729
              d = -0.4477
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select

     case(10)
        select case(thisIon%n) ! from nussbaumer and storey 1987
           case(10) ! Ne I
              z = 1
              a = 7.317
              b = -0.5358
              c = 2.4130
              d = 0.4176
              rate = ppb1991(z, a, b, c, d, dble(temperature))
           case(9) ! Ne II
              a = 0.0129
              b =-0.1779
              c = 0.9353
              d =-0.0682
              f = 0.4516
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 2
              a = 11.80
              b = -0.5404
              c = 3.0300
              d = 0.2050
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(8) ! Ne III
              a = 3.6781
              b = 14.1481
              c = 17.1175
              d = -0.5017
              f = 0.2313
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 3
              a = 5.841
              b = -0.5921
              c = 0.4498
              d = 0.6395
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(7) ! Ne IV
              a =-0.0254
              b = 5.5365
              c = 17.0727
              d = -0.7225
              f = 0.1702
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 4
              a = 15.550
              b = -0.4825
              c = 3.2740
              d = 0.3030
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(6) ! Ne V
              a = -0.0141
              b = 33.8479
              c = 43.1608
              d =-1.6072
              f = 0.1942
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 5
              a = 7.538
              b = -0.5540
              c = 1.2960
              d = 0.3472
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case(5) ! Ne VI
              a = 19.9280
              b = 235.0536
              c = 152.5096
              d = 9.1413
              f = 0.1282
              t = temperature / 10000.
              rate = nussbaumerStorey1983(t,a,b,c,d,f)
              z = 6
              a = 5.239
              b = -0.5966
              c = 0.7135
              d = 0.4537
              rate = rate + ppb1991(z, a, b, c, d, dble(temperature))
           case DEFAULT
              write(*,*) "No recombination rate for ",thisIon%species
              rate  = 0.
        end select

     case(16)
        select case(thisIon%n) ! from nussbaumer and storey 1987
           case(16) ! S I 
              rate = 3.e-13
              rate = rate + svs1982(dble(temperature), 4.10D-13, 6.30D-1)
           case(15) ! S II
              rate = 3.e-30 ! page 1344, para 2, kenny's paper
              rate = rate + svs1982(dble(temperature), 1.80D-12, 6.86D-1)
           case(14) ! S  III
              rate = 1.5e-11
              rate = rate + svs1982(dble(temperature), 2.70D-12, 7.45D-1)
           case(13) ! S IV
              rate = 2.5e-11
              rate = rate + svs1982(dble(temperature), 5.70D-12, 7.55D-1)
        end select
     case DEFAULT
        write(*,*) "No recombination rate for ",thisIon%species
        rate  = 0.
  end select
end function recombRate

!Thaw - dumpLexington is incompatiable with MPI - will possibly move this subroutine to mpi_amr_mod
subroutine dumpLexingtonMPI(grid, epsoverdt, nIter)
  use mpi
  type(GRIDTYPE) :: grid
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  integer :: i, j
  real(double) :: r, theta, phi
  real(double) :: t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
  real(double) :: oirate, oiirate, oiiirate, oivrate
  real(double) :: v, epsoverdt
  real :: fac
  real(double) :: hHeating, heHeating, totalHeating, heating, nh, nhii, nheii, ne
  real(double) :: cooling, dustHeating
  real(double) :: netot
  character(len=80) :: datFilename!, mpiFilename
  integer :: nIter

  !dumpLexingtonMPI specific variables
  integer :: ierr
  integer, parameter :: nStorage = 26, tag=50
  real(double) :: tempStorage(nStorage), tval
  integer :: status(MPI_STATUS_SIZE)
  integer :: sendThread
  integer, parameter :: nPoints = 500
  type(VECTOR) :: position, startPoint, endPoint, direction, octVec
  logical :: stillLooping
  logical, parameter :: useNiter=.false. ! Tag filename with iteration number?
  
  if ( useNiter ) then 
     write(datFilename,'(a,i2.2,a)') "lexington",niter,".dat"
  else
     write(datFilename,'(a,i2.2,a)') "lexington.dat"
  end if

!  write(mpiFilename,'(a, i4.4, a)') "lexington",niter,".grid"
!  call writeAmrGrid(mpiFilename, .false., grid)


  startPoint = vector(0.d0, 0.d0, 0.d0)
  endPoint = vector(4.4d9, 0.d0, 0.d0)

  thisOctal => grid%octreeRoot
  position = startPoint
  direction = endPoint - startPoint
  call normalize(direction)

  !Collate results to write to file in rank 0
  if(myRankGlobal == 0) then
        open(20,file=datFileName,form="formatted",status="unknown")
        open(21,file="orates.dat",form="formatted",status="unknown")
        open(22,file="ne.dat",form="formatted",status="unknown")

        do i=1, 500
           r = (1.+7.d0*dble(i-1)/499.d0)*pctocm/1.e10
           position = vector(r, 0.d0, 0.d0)
           call findSubcellLocal(position, thisOctal, subcell)
           sendThread = thisOctal%mpiThread(subcell)
           
           t=0;hi=0; hei=0;oii=0;oiii=0;cii=0;ciii=0;civ=0;nii=0;niii=0;niv=0;nei=0;neii=0;neiii=0;neiv=0;ne=0.
           oirate = 0; oiirate = 0; oiiirate = 0; oivrate = 0
           heating = 0.d0; cooling = 0.d0
           

           call MPI_SEND(r, 1, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
           call MPI_RECV(tempStorage, nStorage, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, status, ierr)

           hi = tempStorage(1)
           hei = tempStorage(2)
           oii = tempStorage(3)
           oiii = tempStorage(4)
           cii = tempStorage(5)
           ciii = tempStorage(6)
           civ = tempStorage(7)
           nii = tempStorage(8)
           niii = tempStorage(9)
           niv = tempStorage(10)
           nei = tempStorage(11)
           neii = tempStorage(12)
           neiii = tempStorage(13)
           neiv = tempStorage(14)
           ne = tempStorage(15)
           oirate = tempStorage(16)
           oiirate = tempStorage(17)
           oiiirate = tempStorage(18)
           oivrate = tempStorage(19)
           heating = tempStorage(20)
           cooling = tempStorage(21)
           netot = tempStorage(22)
           tVal = tempStorage(23)
           t = tempStorage(24)

        write(21,'(f6.3,1p,6e12.3,0p)') r*1.e10/pctocm,heating,cooling,oirate,oiirate,oiiirate,oivrate

        write(20,'(f6.3,f9.1,  14f8.3)') &
             r*1.e10/pctocm,t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
        write(22,*) r*1.e10/pctocm,netot
  
     end do
     do sendThread = 1, nHydroThreadsGlobal
        r = 1.d30
        call MPI_SEND(r, 1, MPI_DOUBLE_PRECISION, sendThread, tag, localWorldCommunicator, ierr)
     enddo
     close(20)
     close(21)
     close(22)
     goto 555

  !Other ranks send data to 0 for collation
  else
       stillLooping = .true.
       do while(stillLooping)
          call MPI_RECV(r, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
          if (r > 1.d29) then
             stillLooping = .false.
          else
           t=0;hi=0; hei=0;oii=0;oiii=0;cii=0;ciii=0;civ=0;nii=0;niii=0;niv=0;nei=0;neii=0;neiii=0;neiv=0;ne=0.
           oirate = 0; oiirate = 0; oiiirate = 0; oivrate = 0; netot = 0; tval = 0
           heating = 0.d0; cooling = 0.d0


     do j = 1, 100
        call randomNumberGenerator(getDouble=theta)
        theta = theta * Pi
        call randomNumberGenerator(getDouble=phi)
        phi = phi * twoPi

        octVec = VECTOR(r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta))

        call amrgridvalues(grid%octreeRoot, octVec,  foundOctal=thisOctal, foundsubcell=subcell)

        nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
        nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
        nh = thisOctal%nh(subcell)
        ne = thisOctal%ne(subcell)

        v = cellVolume(thisOctal, subcell)

        HI = HI + thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon))
        HeI = HeI + thisOctal%ionfrac(subcell,returnIonNumber("He I", grid%ion, grid%nIon))
        OII = OII + thisOctal%ionfrac(subcell,returnIonNumber("O II", grid%ion, grid%nIon))
        OIII = OIII + thisOctal%ionfrac(subcell,returnIonNumber("O III", grid%ion, grid%nIon))
        CII = CII + thisOctal%ionfrac(subcell,returnIonNumber("C II", grid%ion, grid%nIon))
        CIII = CIII + thisOctal%ionfrac(subcell,returnIonNumber("C III", grid%ion, grid%nIon))
        CIV = CIV + thisOctal%ionfrac(subcell,returnIonNumber("C IV", grid%ion, grid%nIon))
        NII = NII + thisOctal%ionfrac(subcell,returnIonNumber("N II", grid%ion, grid%nIon))
        NIII = NIII + thisOctal%ionfrac(subcell,returnIonNumber("N III", grid%ion, grid%nIon))
        NIV = 0.d0 !atomic data needs updated
        NeI = 0.d0 !atomic data needs updated
        NeII = NeII + thisOctal%ionfrac(subcell,returnIonNumber("Ne II", grid%ion, grid%nIon))
        NeIII = NeIII + thisOctal%ionfrac(subcell,returnIonNumber("Ne III", grid%ion, grid%nIon))
        NeIV = 0.d0 !atomic data needs updated  

        netot = netot + thisOctal%ne(subcell)
        call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDT)
        heating = heating + totalHeating
        fac = 1.
        oirate = oirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O I", grid%ion,grid%nIon)))
        oiirate = oiirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O II", grid%ion, grid%nIon)))
        oiiirate = oiiirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O III", grid%ion, grid%nIon)))
        t  = t + thisOctal%temperature(subcell)
     enddo

     hi = hi / 100.; hei = hei/100.; oii = oii/100.; oiii = oiii/100.; cii=cii/100.
     ciii = ciii/100; civ=civ/100.; nii =nii/100.; niii=niii/100.; niv=niv/100.
     nei=nei/100.;neii=neii/100.; neiii=neiii/100.; neiv=neiv/100.;t=t/100.
     netot = netot / 100.

     oirate = oirate / 100.
     oiirate = oiirate / 100.
     oiiirate = oiiirate / 100.
     oivrate = oivrate / 100.
     heating = heating / 100.
     cooling = cooling / 100.

     if(hi < 1.e-10) then
        write(*,*) "H I frac is low! ! ! ! ! !"
        write(*,*) "hi ", hi
     else if (hei < 1.e-10) then
        write(*,*) "HE I frac is low! ! ! ! ! !"
        write(*,*) "hei", hei
        
     end if
     hi = log10(max(hi, 1d-10))
     hei = log10(max(hei, 1d-10))
     oii = log10(max(oii, 1d-10))
     oiii = log10(max(oiii, 1d-10))
     cii = log10(max(cii, 1d-10))
     ciii = log10(max(ciii, 1d-10))
     civ = log10(max(civ, 1d-10))
     nii = log10(max(nii, 1d-10))
     niii = log10(max(niii, 1d-10))
     niv= log10(max(niv, 1d-10))
     nei = log10(max(nei, 1d-10))
     neii = log10(max(neii, 1d-10))
     neiii = log10(max(neiii, 1d-10))
     neiv = log10(max(neiv, 1d-10))
     ne = log10(max(ne,1.d-10))

     !Store quantities to send to rank 0
     !tempStorage(1) = cen%x
     !tempStorage(2) = cen%y
     !tempStorage(3) = cen%z
     tempStorage(1) = hi
     tempStorage(2) = hei
     tempStorage(3) = oii
     tempStorage(4) = oiii
     tempStorage(5) = cii
     tempStorage(6) = ciii
     tempStorage(7) = civ
     tempStorage(8) = nii
     tempStorage(9) = niii
     tempStorage(10) = niv
     tempStorage(11) = nei
     tempStorage(12) = neii
     tempStorage(13) = neiii
     tempStorage(14) = neiv
     tempStorage(15) = ne
     tempStorage(16) = oirate
     tempStorage(17) = oiirate
     tempStorage(18) = oiiirate
     tempStorage(19) = oivrate
     tempStorage(20) = heating
     tempStorage(21) = cooling
     tempStorage(22) = netot
     tempStorage(23) = tVal
     tempStorage(24) = t

     call MPI_SEND(tempStorage, nStorage, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
  endif
enddo

end if
555 continue
end subroutine dumpLexingtonMPI


subroutine getForbiddenLineLuminosity(grid, species, wavelength, luminosity)
  type(GRIDTYPE) :: grid
  character(len=*) :: species
  real(double) :: wavelength
  real(double) :: fac
  real(double), intent(out) :: luminosity
  integer :: iIon, iTrans, i

  iTrans = 0
  iIon = returnIonNumber(species, grid%ion, grid%nIon) 
  do i = 1, grid%ion(iIon)%nTransitions
     fac = grid%ion(iIon)%transition(i)%lambda
     fac = fac - wavelength
     fac = abs(fac)
     fac = fac / wavelength
     if (fac  < 0.001d0) then
        iTrans = i
        exit
     endif
  enddo
  if (iTrans == 0) then
     write(*,*) "No transition found at ",wavelength, "Angstroms"
     stop
  endif
  luminosity = 0.d0
  call sumLineLuminosity(grid%octreeroot, luminosity, iIon, iTrans, grid)
end subroutine getForbiddenLineLuminosity

recursive subroutine sumLineLuminosity(thisOctal, luminosity, iIon, iTrans, grid)
  type(GRIDTYPE) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i, iIon, iTrans
  real(double) :: luminosity, v, rate
  real :: pops(10)
  type(VECTOR) :: rvec
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sumLineLuminosity(child, luminosity, iIon, iTrans, grid)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          v = cellVolume(thisOctal, subcell)

          call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell))
               
          rate =  pops(grid%ion(iion)%transition(iTrans)%j) * grid%ion(iion)%transition(itrans)%energy * &
               grid%ion(iion)%transition(itrans)%a/ergtoev
          rate = rate * grid%ion(iion)%abundance * thisOctal%nh(subcell) * thisOctal%ionFrac(subcell, iion)
          luminosity = luminosity + rate * v * 1.d30
          
       endif
    enddo
  end subroutine sumLineLuminosity

subroutine metalcoolingRate(ionArray, nIons, thisOctal, subcell, nh, ne, temperature, total, debug)
  type(IONTYPE) :: ionArray(*)
  integer :: nIons, subcell
  type(OCTAL) :: thisOctal
  real(double) :: ne, nh
  real :: temperature
  real(double) :: rate
  real(double), intent(out) :: total
  real :: pops(10)
  integer :: i, j
  logical, optional :: debug

  total = 0.d0
  do j = 5, nIons
     if (ionArray(j)%nTransitions > 0) then
        call solvePops(ionArray(j), pops, ne, temperature)
        rate = 0.d0
        do i = 1, ionArray(j)%nTransitions
           rate = rate + pops(ionArray(j)%transition(i)%j)*ionArray(j)%transition(i)%energy*ionArray(j)%transition(i)%a/ergtoev
        enddo
        rate = rate * ionArray(j)%abundance * nh * thisOctal%ionFrac(subcell, j)
        if (rate < 0.) then
           write(100,'(a,a,a,1p,e12.4,a,0p, f10.1,a,1pe12.4)') "Negative contribution from ", &
                trim(ionArray(j)%species),":",rate," at T = ",temperature, &
                " ion frac ",thisOctal%ionFrac(subcell,j)
        endif
        if (present(debug)) then
           if (debug) then
                 write(100,'(a,a,a,1p,e12.4,a,0p, f10.1,a,1pe12.4)') "Contribution from ", &
                      trim(ionArray(j)%species),":",rate," at T = ",temperature, &
                      " ion frac ",thisOctal%ionFrac(subcell,j)
           endif
        endif
     total = total + rate
     endif
  enddo
end subroutine metalcoolingRate
  
subroutine getSahaMilneFreq(table,temperature, thisFreq)
  type(SAHAMILNETABLE) :: table
  real(double) :: temperature, thisfreq, r, t, fac
  integer :: i, j

  t = max(5000.d0, min(20000.d0, temperature))
  call locate(table%temp, table%nTemp, t, i)
  call randomNumberGenerator(getDouble=r)
  call locate(table%Clyc(i,1:table%nfreq), table%nFreq, r, j)
  fac = (r - table%Clyc(i,j))/(table%Clyc(i,j+1)-table%cLyc(i,j))
  thisFreq = table%freq(j) + fac * (table%freq(j+1)-table%freq(j))
end subroutine getSahaMilneFreq

subroutine twoPhotonContinuum(thisFreq)

! based on table ii of drake, victor, dalgarno, 1969, PhyRev Vol 180, pg 25

  real(double) :: thisFreq
  real :: y(21) = (/ 0., 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, &
       0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500 /)
  real :: hei(21) = (/ 0., 7.77e0, 2.52e1, 4.35e1, 5.99e1, 7.42e1, 8.64e1, 9.69e1, 1.06e2, 1.13e2, 1.20e2, 1.25e2, &
       1.30e2, 1.34e2, 1.37e2, 1.40e2, 1.42e2, 1.43e2, 1.45e2, 1.45e2, 1.45e2 /)
  real :: freq = 3.86e15, fac, r
  real :: prob(21)
  integer :: i
  prob(1) = 0.
  do i = 2, 21
     prob(i) = prob(i-1) + (y(i)-y(i-1)) * hei(i)
  enddo
  prob(1:21) = prob(1:21)/prob(21)
  thisFreq = 0.
  do while((thisFreq*hcgs*ergtoev) < 13.6)
     call randomNumberGenerator(getReal=r)
     call locate(prob, 21, r, i)
     fac = y(i) + ((r - prob(i))/(prob(i+1)-prob(i)))*(y(i+1)-y(i))
     thisFreq = (1.-fac)*freq
  enddo
end subroutine twoPhotonContinuum

subroutine getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  integer :: subcell
  real(double) :: v, epsOverDeltaT
  real(double), intent(out) :: hHeating, heHeating, totalHeating, dustHeating

  dustHeating  = 0.d0
  v = cellVolume(thisOctal, subcell)
  Hheating= thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,1) * grid%ion(1)%abundance &
       * (epsOverDeltaT / (v * 1.d30))*thisOctal%Hheating(subcell) ! equation 21 of kenny's
  if (grid%nion > 2) then
     Heheating= thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,3) * grid%ion(3)%abundance &
          * (epsOverDeltaT / (v * 1.d30))*thisOctal%Heheating(subcell) ! equation 21 of kenny's
  endif
  dustHeating = (epsOverDeltaT / (v * 1.d30))*thisOctal%distanceGrid(subcell) ! equation 14 of Lucy 1999
  totalHeating = (Hheating + HeHeating + dustHeating)
  
end subroutine getHeating

real(double) function returnGamma(table, temp, freq) result(out)
  use utils_mod, only: hunt
  type(GAMMATABLE) :: table
  real(double) :: temp , freq
  integer,save :: i, j
  real(double) :: tfac, ffac, gamma1, gamma2, logT, logF

  logT = log10(min(temp,table%temp(table%nTemp)))
  logF = log10(min(freq, table%freq(table%nfreq)))

  call hunt(table%temp, table%nTemp, logT, i)
  call hunt(table%freq, table%nFreq, logF, j)

  if (i == table%nTemp) i = i - 1
  if (j == table%nFreq) j = j - 1


  if (logF >= table%freq(1).and.logF <= table%freq(table%nFreq)) then
     tfac  = (logT - table%temp(i))/(table%temp(i+1) - table%temp(i))
     ffac  = (logF - table%freq(j))/(table%freq(j+1) - table%freq(j))
     
     gamma1 = table%gamma(j,i) + tfac*(table%gamma(j, i+1) - table%gamma(j,i))
     gamma2 = table%gamma(j+1,i) + tfac*(table%gamma(j+1, i+1) - table%gamma(j+1,i))
     out = gamma1 + ffac*(gamma2 - gamma1)
     out = 10.d0**out
  else
     out = tiny(out)
  endif

end function returnGamma

subroutine addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
  use inputs_mod, only: hOnly
  type(GRIDTYPE) :: grid
  TYPE(OCTAL) :: thisOctal
  integer :: subcell
  integer :: nFreq
  integer :: i, k, iIon, n1, n2
  real(double) :: freq(:), spectrum(:), dfreq(:)
  real(double) :: jnu
!  real :: e
  real(double) :: hxsec
  real(double), parameter :: statisticalWeight(3) = (/ 2.d0, 0.5d0, 2.d0 /)


  ! do Saha-Milne continua for H, HeI and HeII

  do k = 1 , 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (k == 1) iIon = 1
     if (k == 2 .and. .not. hOnly) iIon = 3
     if (k == 3 .and. .not. hOnly) iIon = 4

     call locate(freq, nfreq, grid%ion(iIon)%nuThresh, n1)
     n2 = nFreq
     if (iIon == 3) then
        call locate(freq, nfreq, grid%ion(iIon+1)%nuThresh, n2)
     endif
     do i = n1, n2

!        e = freq(i) * hcgs* ergtoev

       hxSec = returnxSec(grid%ion(iIon), freq(i), iFreq=i)

!        call phfit2(grid%ion(iIon)%z, grid%ion(iIon)%n, grid%ion(iIon)%outerShell , e , hxsec)

        jnu = tiny(jnu)

        if (hxSec > 0.) then
           if (thisOctal%temperature(subcell) > 100.) then
              jnu = statisticalWeight(k) * ((hcgs*freq(i)**3)/(cSpeed**2)) * &
                   ((hcgs**2) /(twoPi*mElectron*Kerg*thisOctal%temperature(subcell)))**(1.5d0) * &
                   (hxsec/1.d10) *  &
                   exp(-hcgs*(freq(i)-grid%ion(iIon)%nuThresh)/(kerg*thisOctal%temperature(subcell)))
              jnu = jnu * thisOctal%ne(subcell) *(thisOctal%nh(subcell) * &
                   thisOctal%ionFrac(subcell,iIon+1) * grid%ion(iIon)%abundance)
           else
              jnu = tiny(jnu)
           endif
        endif

        spectrum(i) = spectrum(i) + jnu * dFreq(i) * fourPi 
     enddo
  enddo

end subroutine addLymanContinua


subroutine addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, table)

! Ferland 1980 PASP 92 596

  type(GAMMATABLE) :: table(3)
  integer :: nFreq
  real(double) :: spectrum(:), freq(:), dfreq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: i, k, iEnd, iIon
  real(double) :: fac

  do k = 1 , 3

     if (k == 1) iIon = 1
     if (k == 2) iIon = 3
     if (k == 3) iIon = 4

     call locate(freq, nFreq, grid%ion(iIon)%nuThresh, iEnd)
     do i = 2, iEnd
        fac = returnGamma(table(k), dble(thisOctal%temperature(subcell)) , freq(i))
        fac = fac*1.d-40 ! units of 10^-40 erg/s/cm/cm/cm
        fac = fac * thisOctal%ne(subcell) * thisOctal%nh(subcell) &
             * thisOctal%ionFrac(subcell,iIon+1) * grid%ion(iIon)%abundance
        spectrum(i) = spectrum(i) + fac*dfreq(i)
     enddo
  enddo

end subroutine addHigherContinua

subroutine addHydrogenRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
  integer :: nFreq
  real(double) :: spectrum(:), freq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  real :: emissivity(3:15,2:8)
  integer :: ilow, iup
  real(double) :: lymanalpha, hBeta, LineEmissivity, lineFreq, energy
  integer :: i

  real(double) :: lambdaTrans(20,20) = reshape( source=&
    (/ 000000.d-8,1215.67d-8,1025.72d-8,992.537d-8,949.743d-8,937.803d-8,930.748d-8,926.226d-8,923.150d-8,920.963d00,&
       919.352d-8,918.129d-8,917.181d-8,916.429d-8,915.824d-8,915.329d-8,914.919d-8,914.576d-8,914.286d-8,914.039d-8,&  
       0.00000d-8,0000000d-8,6562.80d-8,4861.32d-8,4340.36d-8,4101.73d-8,3970.07d-8,3889.05d-8,3835.38d-8,3797.90d-8,&
       3770.63d-8,3750.15d-8,3734.37d-8,3721.94d-8,3711.97d-8,3703.85d-8,3697.15d-8,3691.55d-8,3686.83d-8,3682.81d-8,&  
       0.00000d-8,0.00000d-8,0000000d-8,18751.0d-8,12818.1d-8,10938.1d-8,10049.4d-8,9545.98d-8,9229.02d-8,9014.91d-8,&
       8862.79d-8,8750.47d-8,8665.02d-8,8598.39d-8,8545.39d-8,8502.49d-8,8467.26d-8,8437.96d-8,8413.32d-8,8392.40d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,40512.0d-8,26252.0d-8,21655.0d-8,19445.6d-8,18174.1d-8,17362.1d-8,&
       16806.5d-8,16407.2d-8,16109.3d-8,15880.5d-8,15700.7d-8,15556.5d-8,15438.9d-8,15341.8d-8,15260.6d-8,15191.8d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,74578.0d-8,46525.0d-8,37395.0d-8,32961.0d-8,30384.0d-8,&
       28722.0d-8,27575.0d-8,26744.0d-8,26119.0d-8,25636.0d-8,25254.0d-8,24946.0d-8,24693.0d-8,24483.0d-8,24307.0d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,123680.d-8,75005.0d-8,59066.0d-8,51273.0d-8,&
       46712.0d-8,43753.0d-8,41697.0d-8,40198.0d-8,39065.0d-8,38184.0d-8,37484.0d-8,36916.0d-8,36449.0d-8,36060.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,190570.d-8,113060.d-8,87577.0d-8,&
       75061.0d-8,67701.0d-8,62902.0d-8,59552.0d-8,57099.0d-8,55237.0d-8,53783.0d-8,52622.5d-8,51679.0d-8,50899.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,0000000d-8,277960.d-8,162050.d-8,&
       123840.d-8,105010.d-8,93894.0d-8,86621.0d-8,81527.0d-8,77782.0d-8,74930.0d-8,72696.0d-8,70908.0d-8,69448.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,388590.d-8,&  
       223340.d-8,168760.d-8,141790.d-8,125840.d-8,115360.d-8,108010.d-8,102580.d-8,98443.0d-8,95191.0d-8,92579.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       525200.d-8,298310.d-8,223250.d-8,186100.d-8,164070.d-8,149580.d-8,139380.d-8,131840.d-8,126080.d-8,121530.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,690500.d-8,388320.d-8,288230.d-8,238620.d-8,209150.d-8,189730.d-8,176030.d-8,165900.d-8,158120.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,887300.d-8,494740.d-8,364610.d-8,300020.d-8,261610.d-8,236260.d-8,218360.d-8,205090.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,111800.d-7,619000.d-8,453290.d-8,371000.d-8,322000.d-8,289640.d-8,266740.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,138600.d-7,762300.d-8,555200.d-8,452220.d-8,390880.d-8,350300.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,169400.d-7,926100.d-8,671200.d-8,544400.d-8,468760.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,204400.d-7,111200.d-7,802300.d-8,648200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,243800.d-7,132100.d-7,949200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,288200.d-7,155400.d-7,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,337400.d-7,&  
       0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /), shape=(/20,20/))


  emissivity(15, 2:8) = (/ 0.156E-01, 0.541E-02, 0.266E-02, 0.154E-02, 0.974E-03, 0.657E-03, 0.465E-03 /)
  emissivity(14, 2:8) = (/ 0.192E-01, 0.666E-02, 0.328E-02, 0.190E-02, 0.120E-02, 0.811E-03, 0.573E-03 /)
  emissivity(13, 2:8) = (/ 0.240E-01, 0.832E-02, 0.411E-02, 0.237E-02, 0.151E-02, 0.102E-02, 0.719E-03 /)
  emissivity(12, 2:8) = (/ 0.305E-01, 0.106E-01, 0.524E-02, 0.303E-02, 0.193E-02, 0.130E-02, 0.917E-03 /)
  emissivity(11, 2:8) = (/ 0.397E-01, 0.138E-01, 0.683E-02, 0.396E-02, 0.252E-02, 0.170E-02, 0.119E-02 /)
  emissivity(10, 2:8) = (/ 0.530E-01, 0.184E-01, 0.914E-02, 0.531E-02, 0.338E-02, 0.228E-02, 0.158E-02 /)
  emissivity(9, 2:8) = (/ 0.731E-01, 0.254E-01, 0.127E-01, 0.737E-02, 0.469E-02, 0.312E-02, 0.204E-02 /)
  emissivity(8,2:8) = (/ 0.105E+00, 0.366E-01, 0.183E-01, 0.107E-01, 0.673E-02, 0.425E-02, 0.000E+00 /)
  emissivity(7,2:8) = (/ 0.159E+00, 0.555E-01, 0.278E-01, 0.162E-01, 0.976E-02, 0.000E+00, 0.000E+00 /)
  emissivity(6,2:8) = (/ 0.259E+00, 0.904E-01, 0.455E-01, 0.255E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
  emissivity(5,2:8) = (/ 0.468E+00, 0.163E+00, 0.802E-01, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
  emissivity(4,2:8) = (/ 0.100E+01, 0.339E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
  emissivity(3,2:8) = (/ 0.286E+01, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)


  Hbeta = 10.d0**(-0.870d0*log10(thisOctal%temperature(subcell)) + 3.57d0)
  Hbeta = Hbeta * thisOctal%ne(subcell) *  thisOctal%nh(subcell) * &
       thisOctal%ionFrac(subcell, 2) * grid%ion(1)%abundance

  ! calculate emission due to HI recombination lines [e-25 ergs/s/cm^3]
  do iup = 15, 3, -1
     do ilow = 2, min0(8, iup-1)
        lineEmissivity = emissivity(iup, ilow) * Hbeta * 1.e-25
        energy = (hydE0eV / dble(iLow**2))-(hydE0eV / dble(iUp**2))
        lineFreq = cSpeed/lambdaTrans(iup, ilow)
        if ((lineFreq > freq(1)).and.(lineFreq < freq(nfreq))) then
           call locate(freq, nFreq, lineFreq, i)
           i = i + 1
           spectrum(i) = spectrum(i) + lineEmissivity
        endif
     end do
  end do

  LymanAlpha = 10.d0**(-0.897*log10(thisOctal%temperature(subcell)) + 5.05d0) * &
       thisOctal%ne(subcell) *(thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * &
       grid%ion(1)%abundance) * 1.d-25
  lineFreq = cSpeed/1215.67D-8
  if ((lineFreq > freq(1)).and.(lineFreq < freq(nfreq))) then
     call locate(freq, nFreq, lineFreq, i)
     i = i + 1
     spectrum(i) = spectrum(i) + lymanAlpha
  endif
  

end subroutine addHydrogenRecombinationLines

real(double) function getPhotonFreq(nfreq, freq, spectrum) result(Photonfreq)
  integer :: nFreq
  real(double) :: freq(:), spectrum(:), fac, r
  real(double), allocatable :: tSpec(:)
  integer :: i

  allocate(tSpec(1:nFreq))
  
  tSpec(1:nFreq) = spectrum(1:nFreq)

  do i = 2, nFreq
     tSpec(i) = tSpec(i) + tSpec(i-1)
  enddo
  tSpec(1:nFreq) = tSpec(1:nFreq) - tSpec(1)
  if (tSpec(nFreq) > 0.d0) then
     tSpec(1:nFreq) = tSpec(1:nFreq) / tSpec(nFreq)
     call randomNumberGenerator(getDouble=r)
     call locate(tSpec, nFreq, r, i)
     fac = (r - tSpec(i)) / (tSpec(i+1)-tSpec(i))
     photonFreq = freq(i) + fac * (freq(i+1)-freq(i))
  else
     photonFreq = cSpeed / (1.d-8 * 100.e4)
  endif
  deallocate(tSpec)
end function getPhotonFreq


subroutine addForbiddenLines(nfreq, freq,  spectrum, thisOctal, subcell, grid)

  integer :: nFreq
  real(double) :: spectrum(:), freq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: iIon, iTrans, j
  real :: pops(10)
  real(double) :: rate, lineFreq


  do iIon = 3, grid%nIon
     do iTrans = 1, grid%ion(iIon)%nTransitions
        call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell))
        rate =  pops(grid%ion(iion)%transition(iTrans)%j) * grid%ion(iion)%transition(itrans)%energy * &
             grid%ion(iion)%transition(itrans)%a/ergtoev
        rate = rate * grid%ion(iion)%abundance * thisOctal%nh(subcell) * thisOctal%ionFrac(subcell, iion)
        lineFreq  =  (grid%ion(iion)%transition(itrans)%energy / ergtoEV) / hCgs
        if ((lineFreq > freq(1)).and.(lineFreq < freq(nfreq))) then
           call locate(freq, nFreq, lineFreq, j)
           j = j + 1
           spectrum(j) = spectrum(j) + rate
        endif
     enddo
  enddo

end subroutine addForbiddenLines


subroutine findForbiddenLine(lambda, grid, iIon, iTrans)
  integer :: iIon, iTrans
  type(GRIDTYPE) :: grid
  real(double) :: lambda, lineLambda
  logical :: foundLine
  integer :: i, j

  foundLine = .false.
  do i = 3, grid%nIon
     do j = 1, grid%ion(iIon)%nTransitions
        lineLambda  =  cspeed * hcgs / (grid%ion(i)%transition(j)%energy / ergtoEV) / (angstromToCm)
        if (abs((lineLambda -lambda)/lambda) < 0.01d0) then
           foundLine = .true.
           iIon = i
           iTrans = j
         endif
     enddo
  enddo
  if (.not.foundLine) then
     write(*,*) "Error finding forbidden line at: ",lambda
     stop
  endif
end subroutine findForbiddenLine

       
subroutine addHeRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)

  integer :: nFreq
  real(double) :: spectrum(:), freq(:)
  type(OCTAL) :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: i, j, k
  real :: fac, t, aj ,bj, cj
  real(double) :: lineFreq !, lambda
  real(double) :: emissivity
!  real :: heII4686
!  integer :: ilow, iup
  integer,parameter :: nHeIILyman = 4
!  real(double) :: heIILyman(4)
!  real(double) :: freqheIILyman(4) = (/ 3.839530, 3.749542, 3.555121, 2.99963 /)

  ! HeI lines 

  call locate(heIrecombinationNe, 3, real(log10(thisOctal%ne(subcell))), i)
  fac = real ( &
       (log10(thisOctal%ne(subcell)) - heIrecombinationNe(i))/(heIrecombinationNe(i+1)-heIrecombinationNe(i)) )

  do j = 1, 32
     aj = heIrecombinationFit(j,i,1) + fac*(heIrecombinationfit(j,i+1,1)-heIrecombinationfit(j,i,1))
     bj = heIrecombinationFit(j,i,2) + fac*(heIrecombinationfit(j,i+1,2)-heIrecombinationfit(j,i,2))
     cj = heIrecombinationFit(j,i,3) + fac*(heIrecombinationfit(j,i+1,3)-heIrecombinationfit(j,i,3))
     t = thisOctal%temperature(subcell)/1.e4
     emissivity = aj * (t**bj) * exp(cj / t) ! Benjamin et al. 1999 ApJ 514 307
     emissivity = emissivity * thisOctal%ne(subcell) * thisOctal%nh(subcell) * &
          thisOctal%ionFrac(subcell, 3) * grid%ion(3)%abundance

     lineFreq = cspeed / (heiRecombinationLambda(j)*1.e-8)
     call locate(freq, nFreq, lineFreq, k)
     k = k + 1
     spectrum(k) = spectrum(k) + emissivity
  enddo

end subroutine addHeRecombinationLines

subroutine addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
  use utils_mod, only: hunt

  integer :: nFreq
  real(double) :: spectrum(:), freq(:), dfreq(:)
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  type(GRIDTYPE) :: grid
  integer :: nLambda
  real :: lamArray(:)
  integer :: i, iLam
  real :: thisLam
  real(double), allocatable :: kAbsArray(:)
  integer, allocatable, save :: indexLam(:)
  logical, save :: firstTime = .true.

  if (firstTime) then
     firstTime = .false.
     allocate(indexLam(1:nFreq))
     indexLam = 0
     do i = 1, nFreq
        thisLam = real(cSpeed / freq(i)) * 1.e8
        if ((thisLam >= lamArray(1)).and.(thisLam <= lamArray(nlambda))) then
           call hunt(lamArray, nLambda, real(thisLam), iLam)
           indexLam(i) = iLam
        endif
     enddo
  endif



  allocate(kAbsArray(1:nlambda))

  call returnKappa(grid, thisOctal, subcell, kappaAbsArray=kAbsArray)

  do i = 1, nFreq
     if (indexLam(i) /= 0) then
        spectrum(i) = spectrum(i) + bnu(freq(i), dble(thisOctal%temperature(subcell))) * &
             kAbsArray(indexLam(i)) *1.d-10* dFreq(i) * fourPi
     endif
  enddo

!  do i = 1, nFreq
!     thisLam = real(cSpeed / freq(i)) * 1.e8
!     if ((thisLam >= lamArray(1)).and.(thisLam <= lamArray(nlambda))) then
!        call hunt(lamArray, nLambda, real(thisLam), iLam)
!        spectrum(i) = spectrum(i) + bnu(freq(i), dble(thisOctal%temperature(subcell))) * &
!             kAbsArray(iLam) *1.d-10* dFreq(i) * fourPi
!     endif
!  enddo

  deallocate(kAbsArray)

end subroutine addDustContinuum


subroutine readHeIRecombinationLinesFit()
  character(len=200) :: filename, datadirectory
  integer :: i

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"bss1999.dat"

  open(40,file=filename,status="old",form="formatted")
  do i = 1, 32
     read(40, *) HeIRecombinationLambda(i), heIRecombinationFit(i,1,1:3), &
           heIRecombinationFit(i,2,1:3),  heIRecombinationFit(i,3,1:3)
  enddo
  close(40)
!  HeiRecombinationFit(1:32, 1:3, 1:3) = log10(HeiRecombinationFit(1:32, 1:3, 1:3))
  heIrecombinationNe(1) = 2.
  heIrecombinationNe(2) = 4.
  heIrecombinationNe(3) = 6.
end subroutine readHeIRecombinationLinesFit


subroutine readHeIIrecombination()
  character(len=200) :: filename, datadirectory
  integer :: iup, ilow, i

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"r2b0100.dat"

  open(40,file=filename,status="old",form="formatted")
  
  ! read in HeII recombination lines [e-25 ergs*cm^3/s]                                                                      
  ! (Storey and Hummer MNRAS 272(1995)41)                                                                                    
  open(40, file=filename, status = "old",form="formatted")
  do iup = 30, 3, -1
     read(40,*) (HeIIRecombinationLines(iup, ilow), ilow = 2, min0(16, iup-1))
  end do
  close(40)
end subroutine readHeIIrecombination

recursive subroutine packvalues(thisOctal,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: nCrossings(:)
  real(double) :: photoIonCoeff(:,:)
  real(double) :: hHeating(:)
  real(double) :: heHeating(:)
  real(double) :: distanceGrid(:)
  type(VECTOR) :: radMomVec(:)
  
  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call packvalues(child, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
        nIndex = nIndex + 1
        nCrossings(nIndex) = real(thisOctal%nCrossings(subcell))
        photoIonCoeff(nIndex, :) = thisOctal%photoIonCoeff(subcell, :)
        hHeating(nIndex) = thisOctal%hHeating(subcell)
        heHeating(nIndex) = thisOctal%heHeating(subcell)
        distanceGrid(nIndex) = thisOctal%distanceGrid(subcell)
        radMomVec(nIndex) = thisOctal%radiationMomentum(subcell)
     endif
  enddo
end subroutine packvalues

recursive subroutine unpackvalues(thisOctal,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: nCrossings(:)
  real(double) :: photoIonCoeff(:,:)
  real(double) :: hHeating(:)
  real(double) :: heHeating(:)
  real(double) :: distanceGrid(:)
  type(VECTOR) :: radMomVec(:)

  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unpackvalues(child, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          nIndex = nIndex + 1
          thisOctal%nCrossings(subcell) = int(nCrossings(nIndex))
          thisOctal%photoIonCoeff(subcell, :) = photoIonCoeff(nIndex, :)
          thisOctal%hHeating(subcell) = hHeating(nIndex) 
          thisOctal%heHeating(subcell) = heHeating(nIndex) 
          thisOctal%distanceGrid(subcell) = distanceGrid(nIndex) 
          thisOctal%radiationMomentum(subcell) = radMomVec(nIndex)
       endif
    enddo
  end subroutine unpackvalues

recursive subroutine countVoxelsOnThread(thisOctal, nVoxels)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: nVoxels
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call countVoxelsOnThread(child, nVoxels)
                exit
             end if
          end do
       else
        if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
        nVoxels = nVoxels + 1
       endif
    enddo
  end subroutine countVoxelsOnThread

recursive subroutine checkSetsAreTheSame(thisOctal)
  use mpi
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i, ierr
  real(double), allocatable :: temp(:), temp2(:)
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call checkSetsAreTheSame(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          allocate(temp(1:nHydroSetsGlobal),temp2(1:nHydroSetsGlobal))
          temp = 0.d0
          temp(myHydroSetGlobal+1) = dble(thisOctal%temperature(subcell))
          call mpi_allreduce(temp, temp2, nHydroSetsGlobal, &
               MPI_DOUBLE_PRECISION, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
               ierr)
!          if (abs(temp2(1)-temp2(2)) > epsilon(temp2)) write(*,*) "WARNING grid temperatures differ"


          deallocate(temp, temp2)
             
       endif
    enddo
  end subroutine checkSetsAreTheSame

  recursive subroutine  identifyUndersampled(thisOctal)
    use inputs_mod, only : minCrossings
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
!    minCrossings = 100

  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call identifyUndersampled(child)
                exit
             end if
          end do
       else
          if (thisOctal%nCrossings(subcell) < minCrossings) then
             thisOctal%undersampled(subcell) = .true.
          else
             thisOctal%undersampled(subcell) = .false.
          endif
       endif
    enddo
  end subroutine identifyUndersampled

  recursive subroutine  calcContinuumEmissivity(grid, thisOctal, nlambda, lamArray)
    type(GRIDTYPE) :: grid
    integer :: nFreq
    real(double), allocatable :: freq(:), spectrum(:), dfreq(:)
    integer :: nLambda
    real :: lamArray(:)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i


    nFreq = nlambda

    allocate(freq(1:nFreq), spectrum(1:nFreq), dfreq(1:nFreq))

    do i = 1, nFreq
       freq(i) = cSpeed/(lamArray(nFreq-i+1)*1.e-8)
    enddo
    do i = 2, nFreq-1
       dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
    enddo
    dfreq(1) = (freq(2)-freq(1))
    dfreq(nfreq) = (freq(nfreq)-freq(nfreq-1))



  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calcContinuumEmissivity(grid, child, nlambda, lamArray)
                exit
             end if
          end do
       else
          spectrum = 1.d-30
          thisOctal%etaCont(subcell) = 1.d-40
          if (thisOctal%temperature(subcell) > 1.5d0) then
             call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
             call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
             call addHydrogenRecombinationLines(nfreq,  freq, spectrum, thisOctal, subcell, grid)
             !         call addHeRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
             call addForbiddenLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
             call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
             
             do i = 1, nFreq
                thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + spectrum(i)
             enddo
          endif

       endif
    enddo

    deallocate(freq, spectrum, dfreq)

  end subroutine calcContinuumEmissivity


  function getRandomWavelengthPhotoion(grid, thisOctal, subcell, lamArray, nLambda) result(thisLambda)

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real :: thisLambda
    integer :: nFreq
    real(double), allocatable :: freq(:), spectrum(:), tspec(:),lamspec(:), dfreq(:)
    real(double) :: nuStart, nuEnd, r, fac
    integer :: nLambda
    real :: lamArray(:)
    integer :: i

    nFreq = 1000

    allocate(freq(1:nFreq), spectrum(1:nFreq), lamSpec(1:nFreq), dFreq(1:nFreq))
    nuStart = cSpeed / (1000.d4 * 1.d-8)
    nuEnd = cSpeed / (10.d0 * 1.d-8)

    do i = 1, nFreq
       freq(i) = log10(nuStart) + dble(i-1)/dble(nFreq-1) * (log10(nuEnd)-log10(nuStart))
       freq(i) = 10.d0**freq(i)
    enddo
    do i = 2, nFreq-1
       dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
    enddo
    dfreq(1) = (freq(2)-freq(1))
    dfreq(nfreq) = (freq(nfreq)-freq(nfreq-1))


    spectrum = 1.d-30

    call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
    call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
    call addHydrogenRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
    !                        call addHeRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
    call addForbiddenLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
    call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
    

    do i = 1, nFreq
       lamSpec(i) = cspeed/freq(i)
       spectrum(i) = spectrum(i) * cspeed/(lamspec(i)**2)
    enddo
    lamSpec = lamSpec * 1.d8

    allocate(tSpec(1:nFreq))
  
    tSpec(1:nFreq) = spectrum(1:nFreq)

    do i = 2, nFreq
       tSpec(i) = tSpec(i) + tSpec(i-1)
    enddo
    tSpec(1:nFreq) = tSpec(1:nFreq) - tSpec(1)
    if (tSpec(nFreq) > 0.d0) then
       tSpec(1:nFreq) = tSpec(1:nFreq) / tSpec(nFreq)
       call randomNumberGenerator(getDouble=r)
       call locate(tSpec, nFreq, r, i)
       fac = (r - tSpec(i)) / (tSpec(i+1)-tSpec(i))
       thisLambda = real(lamspec(i) + fac * (lamspec(i+1)-lamspec(i)))
    else
       thisLambda = 1000.e4
  endif
    
  deallocate(freq, spectrum, lamSpec, dfreq)

  end function getRandomWavelengthPhotoion


  recursive subroutine calculateEnergyFromTemperature(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: eThermal
    real(double) :: mu

    do subcell = 1, thisOctal%maxChildren
       if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateEnergyFromTemperature(child)
                exit
             end if
          end do
       else
          mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
          if (mu == 0.d0) then
             write(*,*) "nh ",thisOctal%nh(subcell)
             write(*,*) "ionfrac ",thisOctal%ionFrac(subcell,:)
             write(*,*) "thisOctal%ndepth ",thisOctal%nDepth
             write(*,*) "thisOctal%rho ",thisOctal%rho(subcell)
          endif
          eThermal = 1.5d0 * thisOctal%temperature(subcell)*kerg/(mu*mHydrogen)
          thisOctal%energy(subcell) =  eThermal
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)

       endif
    enddo
  end subroutine calculateEnergyFromTemperature

  subroutine createImageSplitGrid(grid, nSource, source, imageNum)
    use inputs_mod, only: nPhotons, gridDistance
    use image_utils_mod
    use mpi
    integer, intent(in) :: imageNum
    real :: lambdaImage
    type(GRIDTYPE) :: grid
    character(len=80) :: imageFilename 
    character(len=80) :: outputImageType
    integer :: nSource
    type(SOURCETYPE) :: source(:), thisSource
    type(PHOTON) :: thisPhoton, observerPhoton
    type(OCTAL), pointer :: thisOctal
    real(double) :: totalFlux
    real(double), allocatable :: totalFluxArray(:), tempTotalFlux(:)
    logical :: directFromSource
    integer :: subcell
    integer(kind=bigInt) :: iPhoton
    integer :: iSource
    integer :: iThread
    type(VECTOR) ::  rHat, observerDirection
    logical :: endLoop, addToiMage
    integer :: newthread
    type(IMAGETYPE) :: thisImage
    logical :: escaped, absorbed, crossedBoundary, photonsStillProcessing, stillSCattering
    real(double) :: totalEmission
    integer :: iLambdaPhoton
    real(double) :: lCore, r
    real(double), allocatable :: threadProbArray(:)
    integer :: np(10)
    integer :: nDone
    integer :: tag = 41
    integer, allocatable :: nDoneArray(:)
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    integer :: isignal
    real(double) :: powerPerPhoton
    logical :: freefreeImage
    real(double) :: chanceSource, probsource
! Weight of all sources relative to envelope
    real(double) :: weightSource
! Weight of envelope relative to sources
    real(double) :: weightEnv    
! weight of this particular source (for multiple sources)
    real(double) :: thisSourceWeight 
#ifdef USECFITSIO
    real(double) :: theoretical
#endif
    integer :: nLams
    character(len=80) :: message

    allocate(nDoneArray(1:nHydroThreadsGlobal))
    allocate(totalFluxArray(1:nHydroThreadsGlobal))
    allocate(tempTotalFlux(1:nHydroThreadsGlobal))
    allocate(threadProbArray(1:nHydroThreadsGlobal))

    lambdaImage     = getImageWavelength(imageNum)
    outputImageType = getImageType(imageNum)
    imageFilename   = getImageFilename(imageNum)
    observerDirection = getImageViewVec(imageNum)

    call randomNumberGenerator(randomSeed=.true.)
    totalFluxArray = 0.d0

    call zeroEtaCont(grid%octreeRoot)
    
    call quickSublimate(grid%octreeRoot) ! do dust sublimation

    call torus_mpi_barrier ! Why is there a barrier here?

    thisImage = initImage(imageNum)

    call setupGridForImage(grid, outputimageType, lambdaImage, iLambdaPhoton, nsource, source, lcore)

    call computeProbDistAMRMpi(grid, totalEmission, threadProbArray)
    if (myrankglobal == 0) write(*,*) "prob array ", threadProbArray(1:nHydroThreadsGlobal)
    totalEmission = totalEmission * 1.d30

    ! Probability that a photon comes from a source rather than the envelope
    chanceSource = lCore / (lCore + totalEmission)
    ! Fraction of photons actually used to sample the sources
    probSource   = 0.1d0
    ! Weight the source and envelope photons accordingly
    weightSource = chanceSource / probSource
    weightEnv    = (1.d0 - chanceSource) / (1.d0 - probSource) 

    write(message,*) "Total continuum emission ",totalEmission
    call writeInfo (message, FORINFO)
    write(message,*) "Total source emission ",lCore
    call writeInfo (message, FORINFO)
    write(message,*) "Total  emission ",totalEmission+lCore
    call writeInfo (message, FORINFO)
    write(message,*) "Ratio of source to total luminosity: ", chanceSource
    call writeInfo (message, FORINFO)
    write(message,*) "Probability of photon from sources: ", probSource
    call writeInfo (message, FORINFO)
    write(message,*) "Weight of sources: ", weightSource
    call writeInfo (message, FORINFO)
    write(message,*) "Envelope weight: ", weightEnv
    call writeInfo (message, FORINFO)
    if(chanceSource == 0.d0) then
       write(message,*) "Zero probability from source" 
       call writeWarning(message)
       write(message,*) "lCore :", lcore, " totalEmission : ", totalEmission
       call writeInfo(message)
    endif

    powerPerPhoton = (lCore + totalEmission) / dble(nPhotons)
    write(message,*) "power per photon ",powerperphoton
    call writeInfo (message, FORINFO)

    absorbed = .false.
    escaped = .false.

    if (outputimageType == "freefree") then
       freefreeImage = .true.
       nLams = 1
    else
       freefreeImage = .false.
    end if

    Call randomSource(source, nSource, iSource, thisSourceWeight, nLambda=nLams, initialize=.true., &
         lambdaMono=LambdaImage)

    if (myRankGlobal == 0) then
       np = 0
       mainloop: do iPhoton = 1, nPhotons

          thisPhoton%weight = 1.d0
          thisPhoton%stokes = STOKESVECTOR(1.d0*powerPerPhoton, 0.d0, 0.d0, 0.d0)
          thisPhoton%iLam = iLambdaPhoton
          thisPhoton%lambda = grid%lamArray(iLambdaPhoton)
          thisPhoton%observerPhoton = .false.
        
          call randomNumberGenerator(getDouble=r)

          if (r < probSource) then
             call randomSource(source, nSource, iSource, thisSourceWeight)
             thisSource = source(iSource)
             call getPhotonPositionDirection(thisSource, thisPhoton%position, thisPhoton%direction, rHat,grid)
             call findSubcellTD(thisPhoton%position, grid%octreeRoot,thisOctal, subcell)
             iThread = thisOctal%mpiThread(subcell)
             call sendPhoton(thisPhoton, iThread, endloop = .false.) 
             directFromSource = .true.
             thisPhoton%weight = weightSource * thisSourceWeight
          else
             call randomNumberGenerator(getDouble=r)
             if (r < threadProbArray(1)) then
                iThread = 1
             else
                call locate(threadProbArray, SIZE(threadProbArray), r, iThread)
                iThread = iThread + 1
             endif
             np(iThread) = np (iThread) + 1
             thisPhoton%lambda = grid%lamArray(iLambdaPhoton)
             thisPhoton%direction = randomUnitVector()
             call sendPhoton(thisPhoton, iThread, endloop = .false., getPosition=.true.) 
             call receivePhoton(thisPhoton, iSignal)
             directFromSource = .false.
             thisPhoton%weight = weightEnv
          endif

          if ((directFromSource.and.(.not.thisSource%outsideGrid)).or.(.not.directFromSource)) then
             observerPhoton = thisPhoton
             observerPhoton%observerPhoton = .true.
             observerPhoton%tau = 0.d0
             observerPhoton%direction = observerDirection
             call sendPhoton(observerPhoton, iThread, endLoop = .false.)
          endif

       end do mainloop

       photonsStillProcessing = .true.
       do while (photonsStillProcessing)
          do iThread = 1, nHydroThreadsGlobal
             call sendPhoton(thisPhoton, iThread, endLoop = .false., report=.true.)
             call MPI_RECV(nDoneArray(iThread), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, status, ierr)
          enddo
          nDone = SUM(nDoneArray(1:nHydroThreadsGlobal))
          write(*,*) myrankglobal, " thinks ", nDone, " photons have completed "
          if (nDone == nPhotons) photonsStillProcessing = .false.
       enddo
       do iThread = 1, nHydroThreadsGlobal
          call sendPhoton(thisPhoton, iThread, endLoop = .true.)
       enddo

    else
       nDone = 0
       endLoop = .false.
       do while (.not.endLoop)
          
          call receivePhoton(thisPhoton, iSignal)
          

          if (iSignal == 0) then
             thisOctal => grid%octreeRoot
             call randomNumberGenerator(getDouble=r)
             call locateContProbAMR(r,thisOctal,subcell)
             thisPhoton%position = randomPositionInCell(thisOctal, subcell)
             if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) then
                write(*,*) myrankGlobaL, &
                     " Bug: octal not on thread from locatecontprobamr ",thisphoton%position
             endif

             call sendPhoton(thisPhoton, 0, endLoop = .false.)
          endif


          if (iSignal == 1) then
             endLoop = .true.
             goto 777
          end if

          if (iSignal == 2) then
             call MPI_SEND(nDone, 1, MPI_INTEGER, 0, tag, localWorldCommunicator,  ierr)
             goto 777
          end if


          stillScattering = .true.
          do while ((.not.endLoop).and.stillScattering)

             if (thisPhoton%observerPhoton) then
                newThread = -1
                if (.not.freeFreeImage) then
                   call propagateObserverPhoton(grid, thisPhoton, addToImage, newThread)
                else
                   addtoImage = .true.
                endif
                if (addToImage) then
                   call addPhotonToPhotoionImage(observerDirection, thisImage, thisPhoton, totalFluxArray(myRankGlobal))
                   goto 777
                else
                   call sendPhoton(thisPhoton, newThread, endLoop = .false.)
                   goto 777
                endif
             else
                if (.not.freeFreeImage) then
                   call moveToNextScattering(grid, thisPhoton, escaped, absorbed, crossedBoundary, newThread)
                else
                   escaped = .true.
                endif
                
                if (escaped.or.absorbed) then
                   stillScattering = .false.
                   nDone = nDone + 1
                   exit
                endif

                if (.not.crossedBoundary) then
                   call scatterPhotonLocal(thisPhoton)
                   observerPhoton = thisPhoton
                   observerPhoton%observerPhoton = .true.
                   observerPhoton%tau = 0.d0
                   observerPhoton%direction = observerDirection
                   newThread = -2
                   call propagateObserverPhoton(grid, observerPhoton, addToImage, newThread)
                   if (addToImage) then
                      call addPhotonToPhotoionImage(observerDirection, thisImage, observerPhoton, totalFluxArray(myRankGlobal))
                   else
                      call sendPhoton(observerPhoton, newThread, endLoop = .false.)
                   endif
                else
                   call sendPhoton(thisPhoton, newThread, endLoop = .false.)
                   goto 777
                endif
             endif
          end do
777       continue
       enddo
    endif
    call collateImages(thisImage)

    call MPI_ALLREDUCE(totalFluxArray, tempTotalFlux, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, &
         MPI_SUM, localWorldCommunicator, ierr)
     totalFluxArray = tempTotalFlux
     totalFlux = SUM(totalFluxArray(1:nHydroThreadsGlobal))

#ifdef USECFITSIO
    if (myrankGlobal == 0) then
       if(grid%geometry == "point") then
          open (123, file="image_flux.dat", status="unknown")
          theoretical = (lCore/(fourpi * (griddistance*pctocm)**2))
          write(123, *) theoretical
          write(*,*) "theoretical ", theoretical
          close(123)
          call writeFitsImage(thisimage, imageFilename, griddistance*pctocm, "intensity", lambdaImage, &
               pointTest=.true.)
       else if(grid%geometry == "imgTest") then
          call writeFitsImage(thisimage, imageFilename, griddistance*pctocm, "intensity", lambdaImage, &
               cylinderTest=.true.)
       else
          call writeFitsImage(thisimage, imageFilename, griddistance*pctocm, "intensity", lambdaImage)
       end if
    endif
#else
    call writeInfo("FITS not enabled, not writing "//trim(imageFilename),FORINFO)
#endif
    call freeImage(thisImage)
  end subroutine createImageSplitGrid


  subroutine sendPhoton(thisPhoton, iThread, endLoop, report, getPosition)
    use mpi
    type(PHOTON) :: thisPhoton
    integer :: iThread
    logical :: endLoop
    logical, optional :: report, getPosition
    integer, parameter  :: nTemp = 16
    real(double) :: temp(nTemp)
    integer :: ierr
    integer :: tag = 42

    temp(1) = thisPhoton%stokes%i
    temp(2) = thisPhoton%stokes%q
    temp(3) = thisPhoton%stokes%u
    temp(4) = thisPhoton%stokes%v

    temp(5) = thisPhoton%position%x
    temp(6) = thisPhoton%position%y
    temp(7) = thisPhoton%position%z

    temp(8) = thisPhoton%direction%x
    temp(9) = thisPhoton%direction%y
    temp(10) = thisPhoton%direction%z

    temp(11) = thisPhoton%lambda

    temp(12) = 0.d0
    if (endLoop) then
       temp(12) = 2.1d0
    endif
    if (PRESENT(report)) then
       if (report) temp(12) = 3.1d0
    endif
    if (PRESENT(getPosition)) then
       if (getPosition) temp(12) = 1.1d0
    endif


    temp(13) = thisPhoton%tau
    temp(14) = thisPhoton%iLam

    if (thisPhoton%observerPhoton) then
       temp(15) = 1.d0
    else
       temp(15) = 0.d0
    endif

    temp(16) = thisPhoton%weight

    if (iThread == myRankGlobal) then
       write(*,*) "sending to self bug ", ithread
       stop
    endif

    call MPI_SEND(temp, nTemp, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator,  ierr)
    
  end subroutine sendPhoton

  subroutine receivePhoton(thisPhoton, iSignal)
    use mpi
    type(PHOTON) :: thisPhoton
    integer, parameter :: nTemp = 16
    real(double) :: temp(nTemp)
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    integer :: tag = 42
    integer, intent(out) :: iSignal

    call MPI_RECV(temp, nTemp, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, localWorldCommunicator, status, ierr)
    iSignal = -1
    if (temp(12) > 1.d0) then
       iSignal = 0 ! newPhoton
    endif

    if (temp(12) > 2.d0) then
       iSignal = 1 ! abort loop
    endif

    if (temp(12) > 3.d0) then
       iSignal = 2 ! report ndone
    endif

    thisPhoton%stokes%i =     temp(1)  
    thisPhoton%stokes%q =     temp(2)  
    thisPhoton%stokes%u =     temp(3)  
    thisPhoton%stokes%v =     temp(4)  

    thisPhoton%position%x =     temp(5)  
    thisPhoton%position%y =     temp(6)  
    thisPhoton%position%z =     temp(7)  

    thisPhoton%direction%x =     temp(8)  
    thisPhoton%direction%y =     temp(9)  
    thisPhoton%direction%z =     temp(10) 

    thisPhoton%lambda =     real(temp(11)) 

    thisPhoton%tau =     real(temp(13)) 
    thisPhoton%iLam =     nint(temp(14) )

    if (temp(15) > 0.d0) then
       thisPhoton%observerPhoton = .true.
    else
       thisPhoton%observerPhoton = .false.
    endif

    thisPhoton%weight = temp(16)
    
  end subroutine receivePhoton

  subroutine scatterPhotonLocal(thisPhoton)
    type(PHOTON) :: thisPhoton

    thisPhoton%direction = randomUnitVector() !isotropic scattering

  end subroutine scatterPhotonLocal

  subroutine propagateObserverPhoton(grid, thisPhoton, addToImage, newThread)
    type(GRIDTYPE) :: grid
    type(PHOTON) :: thisPhoton
    logical,intent(out) :: addToImage
    logical :: endLoop
    integer,intent(out) :: newThread
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: tVal
    real(double) :: kappaAbsDust, kappaScaDust, kappaExt

    thisOctal => grid%octreeRoot
    call findSubcellTD(thisPhoton%position, grid%octreeRoot, thisOctal, subcell)
    if (myRankGlobal /= thisOctal%mpiThread(subcell)) then
       write(*,*) "PropagateObserverPhoton call with wrong thread ", myRankGlobal, thisOctal%mpiThread(subcell),newThread
       stop
    endif

    addToImage = .false.
    endLoop = .false.
    do while (.not.endLoop)
       if (.not.inSubcell(thisOctal, subcell, thisPhoton%position)) then
          write(*,*) myrankGlobal, " bug in propagate observer ", thisphoton%position
       endif
       call distanceToCellBoundary(grid, thisPhoton%position, thisPhoton%direction, tval, thisOctal, subcell)
       call returnKappa(grid, thisOctal, subcell, ilambda=thisPhoton%ilam, &
            kappaAbs=kappaAbsDust, kappaSca=kappaScaDust)
       kappaExt = kappaAbsDust + kappaScaDust

!       if(grid%geometry == "imgTest") then
!          write(*,*) "kappaExt",kappaExt
!       end if

       thisPhoton%tau = thisPhoton%tau + real(tval * kappaExt)
       thisPhoton%position = thisPhoton%position + (tVal + 1.d-3*grid%halfSmallestSubcell) * thisPhoton%direction
       if (.not.inOctal(grid%octreeRoot, thisPhoton%position)) then
          addToImage = .true.
          endLoop = .true.
       else
          call findSubcellLocal(thisPhoton%position, thisOctal, subcell)
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) then
             newThread = thisOctal%mpiThread(subcell)
             endLoop = .true.
          endif
       endif
    enddo
  end subroutine propagateObserverPhoton

  subroutine moveToNextScattering(grid, thisPhoton, escaped, absorbed, crossedBoundary, newThread)
    type(GRIDTYPE) :: grid
    type(PHOTON) :: thisphoton
    logical, intent(out) :: escaped, absorbed, crossedBoundary
    logical :: scattered
    integer, intent(out) :: newThread
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: tau, thisTau, tVal, r, albedo
    real(double) :: kappaAbsDust, kappaScaDust,  kappaExt
    logical :: endLoop
    
    thisOctal => grid%octreeRoot
    call findSubcellTD(thisPhoton%position, grid%octreeRoot, thisOctal, subcell)
    
    if (myRankGlobal /= thisOctal%mpiThread(subcell)) then
       write(*,*) "moveToNextScattering call with wrong thread ", myRankGlobal, thisOctal%mpiThread(subcell)
       stop
    endif

    crossedboundary = .false.
    absorbed = .false.
    escaped = .false.
    endLoop = .false.
    scattered = .false.
    do while (.not.endLoop)
       if (.not.inSubcell(thisOctal, subcell, thisPhoton%position)) then
          write(*,*) myrankGlobal, " bug in move to next scattering ", thisPhoton%position
       endif
       call distanceToCellBoundary(grid, thisPhoton%position, thisPhoton%direction, tval, thisOctal, subcell)
       call returnKappa(grid, thisOctal, subcell, ilambda=thisPhoton%ilam, &
            kappaAbs=kappaAbsDust, kappaSca=kappaScaDust)
       kappaExt = kappaAbsDust + kappaScaDust

!       if(grid%geometry == "imgTest") then
!          write(*,*) "kappaExt",kappaExt
!       end if


       tau = kappaExt * tVal

       call randomNumberGenerator(getDouble=r)
       thisTau = -log(1.d0-r)

       if (thisTau > tau) then ! photon crosses to boundary
          thisPhoton%position = thisPhoton%position + (tVal + 1.d-3*grid%halfSmallestSubcell) * thisPhoton%direction
          if (.not.inOctal(grid%octreeRoot, thisPhoton%position)) then
             escaped = .true.
             endLoop = .true.
          else
             call findSubcellLocal(thisPhoton%position, thisOctal, subcell)
             if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) then
                newThread = thisOctal%mpiThread(subcell)
                crossedBoundary = .true.
                endLoop = .true.
             endif
          endif
       else

          call findSubcellLocal(thisPhoton%position, thisOctal, subcell)
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) then
             write(*,*) myrankGlobal , " has an error as photon position is on wrong thread. ",thisPhoton%position
          endif

          thisPhoton%position = thisPhoton%position + ((thisTau/tau)*tVal) * thisPhoton%direction

          call findSubcellLocal(thisPhoton%position, thisOctal, subcell)
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) then
             write(*,*) myrankGlobal , " has an error as photon moved outside subcell. ",thistau,tau,tval
             write(*,*) "halfsmallest ",grid%halfSmallestSubcell
          endif
          endLoop = .true.
          albedo = kappaScaDust/kappaExt
!          write(*,*) "test tau  ",tau, " thistau ",thistau, "  albedo ",albedo, " ext ",kappaext, " sca ",kappaScaDust, &
!               "dust ",thisOctal%dustTypeFraction(subcell,1)
          call randomNumberGenerator(getDouble=r)
          if (r < albedo) then
             scattered = .true.
          else
             absorbed = .true.
          endif

       endif
    enddo
  end subroutine moveToNextScattering

   subroutine  computeProbDistAMRMpi(grid, totalEmission, threadProbArray)
     use mpi
     type(GRIDTYPE) :: grid
     real(double) :: totalEmission, totalProb, biasCorrection
     real(double) :: threadProbArray(:)
     real(double), allocatable :: totalEmissionArray(:), totalProbArray(:), tArray(:)
     integer :: ierr, i

     totalEmission = 0.d0

     allocate(totalEmissionArray(1:nHydroThreadsGlobal), totalProbArray(1:nHydroThreadsGlobal), &
          tarray(1:nHydroThreadsGlobal))
     totalEmissionArray = 0.d0
     totalProbArray = 0.d0

     if (myrankGlobal /= 0) then
        call computeProbDist2AMRMpi(grid%octreeRoot,totalEmissionArray(myRankGlobal), totalProbArray(myRankGlobal))
        write(*,*) myrankGlobal, " total emission  ", totalEmissionArray(myrankGlobal)
     endif


     tArray = 0.d0
     call MPI_ALLREDUCE(totalEmissionArray, tArray, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, &
          localWorldCommunicator, ierr)
     totalEmissionArray = tArray

     tArray = 0.d0
     call MPI_ALLREDUCE(totalProbArray, tArray, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, &
          localWorldCommunicator, ierr)
     totalProb = SUM(tArray(1:nHydroThreadsGlobal))
     
     threadProbArray = tArray(1:nHydroThreadsGlobal)
     do i = 2, SIZE(threadProbArray)
        threadProbArray(i) = threadProbArray(i) + threadProbArray(i-1)
     enddo
     if (threadProbArray(SIZE(threadPRobArray)) /= 0.d0) then
        threadProbArray = threadProbArray / threadProbArray(SIZE(threadProbArray))
     else
        threadProbArray = 1.d0
     endif


     totalEmission = SUM(totalEmissionArray(1:nHydroThreadsGlobal))

     if (totalProb /= 0.0) then
        biasCorrection = totalEmission / totalProb
     else
        biasCorrection = 1.0
     end if
     
     if (myrankGlobal /= 0) call computeProbDist3AMRMpi(grid%octreeRoot, biasCorrection, totalProbArray(myRankGlobal))

     deallocate(totalEmissionArray, totalProbArray, tArray)
   end subroutine computeProbDistAMRMpi


  recursive subroutine computeProbDist2AMRMpi(thisOctal, totalEmission, totalProb)

    implicit none

    type(octal), pointer                 :: thisOctal
    real(double), intent(inout) :: totalEmission, totalProb
    
    type(octal), pointer  :: child 
    real(double)          :: dV 
    integer               :: subcell
    integer               :: i
    
    do subcell = 1, thisOctal%maxChildren, 1

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call computeProbDist2AMRMpi(child, totalEmission, totalProb)
                exit
             end if
          end do            
       else
          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then
             dv = cellVolume(thisOctal, subcell)
             
             totalProb = totalProb + dV * &
                  thisOctal%etaCont(subcell) * thisOctal%biasCont3D(subcell)

             totalEmission = totalEmission + dV * thisOctal%etaCont(subcell)

!             if (thisOctal%etaCont(subcell) > 0.d0) then
!                write(*,*) "teff ", thisOctal%temperature(subcell), thisOctal%dustTypeFraction(subcell,1), &
!                     thisOctal%nDepth
!             endif
          endif

       end if

       if (.not.associated(thisOctal%probDistCont)) allocate(thisOctal%probDistCont(1:thisOctal%maxChildren))
       thisOctal%probDistCont(subcell) = totalProb
      
    end do

  end subroutine computeProbDist2AMRMpi


  recursive subroutine computeProbDist3AMRMpi(thisOctal, biasCorrection, totalProb) 

    implicit none

    type(octal), pointer              :: thisOctal
    real(double), intent(in) :: biasCorrection
    real(double), intent(in) :: totalProb
    integer :: subcell
    type(octal), pointer  :: child 
    integer :: nSubcell

    if (thisOctal%nChildren > 0) then
       ! call this subroutine recursively on each of its children
       do subcell = 1, thisOctal%nChildren, 1 
          child => thisOctal%child(subcell)
          call computeProbDist3AMRMpi(child, biasCorrection, totalProb)
       end do 
    end if

    nSubcell = thisOctal%maxChildren

    if (totalProb /= 0.0d0) then
       thisOctal%probDistCont(1:nSubcell) = thisOctal%probDistCont(1:nSubcell) / totalProb
    else
       thisOctal%probDistCont(1:nSubcell) = 1.0d0
    end if

    ! probDist[Line|Cont] are set for all subcells, regardless of whether they
    ! have children or not. bias[Line|Cont]3D, on the other hand, are only set
    ! for childless subcells, so we must only 'correct' valid values.
    do subcell = 1, thisOctal%maxChildren, 1
       if (.not. thisOctal%hasChild(subcell)) then
          thisOctal%biasCont3D(subcell) = thisOctal%biasCont3D(subcell) * biasCorrection
       end if
    end do 
  end subroutine computeProbDist3AMRMpi

#endif    


#ifdef MPI
  subroutine updateGridMPIphoto(grid, amrParComm)
    use gridtype_mod, only : gridtype
    use mpi
    implicit none
    integer :: amrParComm

    type(gridtype) :: grid
    integer :: nVoxels, i
    real, allocatable :: nCrossings(:)
    real, allocatable :: tempRealArray(:)
    real(double), allocatable :: hHeating(:), heHeating(:)
    real(double), allocatable :: photoIonCoeff(:,:)
    real(double), allocatable :: tempDoubleArray(:)
    real(double), allocatable :: distanceGrid(:)
    type(VECTOR), allocatable :: radMomVec(:)
    integer :: ierr, nIndex

    ! FOR MPI IMPLEMENTATION=======================================================

    nVoxels = 0
    call countVoxelsOnThread(grid%octreeRoot,nVoxels)
    allocate(nCrossings(1:nVoxels))
    allocate(hHeating(1:nVoxels))
    allocate(heHeating(1:nVoxels))
    allocate(distanceGrid(1:nVoxels))
    allocate(photoIonCoeff(1:nVoxels, 1:grid%nIon))
    allocate(radMomVec(1:nVoxels))

    nIndex = 0
    call packValues(grid%octreeRoot,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec)


    allocate(tempDoubleArray(nVoxels))
    allocate(tempRealArray(nVoxels))

    do i = 1, grid%nIon
      tempDoubleArray = 0.d0
      call MPI_ALLREDUCE(photoIonCoeff(1:nVoxels,i),tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
          MPI_SUM, amrParComm ,ierr)
       photoIonCoeff(1:nVoxels, i) = tempDoubleArray 
    enddo

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(hHeating,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    hHeating = tempDoubleArray

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(heHeating,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm,ierr)
    heHeating = tempDoubleArray

    tempRealArray = 0.0
    call MPI_ALLREDUCE(nCrossings,tempRealArray,nVoxels,MPI_REAL,&
         MPI_SUM, amrParComm,ierr)
    nCrossings = tempRealArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(distanceGrid,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    distanceGrid = tempDoubleArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(radMomVec(1:nVoxels)%x,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    radMomVec(1:nVoxels)%x = tempDoubleArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(radMomVec(1:nVoxels)%y,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    radMomVec(1:nVoxels)%y = tempDoubleArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(radMomVec(1:nVoxels)%z,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    radMomVec(1:nVoxels)%z = tempDoubleArray 





    deallocate(tempRealArray, tempDoubleArray)
     
    call MPI_BARRIER(amrParComm, ierr) 
    
    nIndex = 0
    call unpackValues(grid%octreeRoot, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec)

    deallocate(nCrossings, photoIonCoeff, hHeating, heHeating, distanceGrid, radMomVec)

  end subroutine updateGridMPIphoto
#endif


end module photoionAMR_mod



#endif
