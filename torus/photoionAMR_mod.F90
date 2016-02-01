#ifdef PHOTOION
!Photoionization module - started on October 4th 2005 by th 


module photoionAMR_mod

#ifdef MPI
use constants_mod
use loadbalance_mod
use messages_mod
use parallel_mod
use photoion_utils_mod
use gridio_mod
use viscosity_mod
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
#ifdef XRAY
use xray_mod
#endif
use vtk_mod
#ifdef PDR
use pdr_mod
use inputs_mod, only: pdrcalc, dumpBisbas
#endif
implicit none

private
public :: photoIonizationloopAMR, createImagesplitgrid, ionizeGrid, &
     neutralGrid, resizePhotoionCoeff, resetNH, hasPhotoionAllocations, allocatePhotoionAttributes, &
     computeProbDistAMRMpi, putStarsInGridAccordingToDensity

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
    logical :: bigPhotonPacket
    logical :: smallPhotonPacket
    logical :: lastPhoton
end type PHOTONPACKET

  real :: heIRecombinationFit(32, 3, 3)
  real :: heIRecombinationLambda(32)
  real :: heIRecombinationNe(3)
  real :: heIIrecombinationLines(3:30, 2:16)
  type(GAMMATABLE) :: gammaTableArray(3) ! H, HeI, HeII
  real(double) :: globalEpsOverDeltaT
  real(double) :: kappaPArray(10,2000)
contains

#ifdef HYDRO
  subroutine radiationHydro(grid, source, nSource, nLambda, lamArray, miePhase, nMuMie)
    use inputs_mod, only : iDump, doselfgrav, readGrid, maxPhotoIonIter, tdump, tend, justDump !, hOnly
    use inputs_mod, only : dirichlet, amrtolerance, nbodyPhysics, amrUnrefineTolerance, smallestCellSize, dounrefine
    use inputs_mod, only : addSinkParticles, cylindricalHydro, vtuToGrid, timedependentRT,dorefine, alphaViscosity
    use inputs_mod, only : UV_vector, spherical, forceminrho
    use starburst_mod
    use viscosity_mod, only : viscousTimescale
    use dust_mod, only : emptyDustCavity, sublimateDust
    use hydrodynamics_mod, only: hydroStep3d, hydrostep3d_amr, calculaterhou, calculaterhov, calculaterhow, &
         calculaterhoe, setupedges, unsetGhosts, setupghostcells, evenupgridmpi, refinegridgeneric, &
         setupx, setupqx, computecouranttime, unrefinecells, selfgrav, sumgasstargravity, transfertempstorage, &
         zerophigas, applysourcepotential, addStellarWind, cutVacuum, setupEvenUpArray, &
         pressureGradientTimestep, computeGravityTimestep, mergeSinks, addSinks, ComputeCourantTimenBody, addSupernovae,&
         perturbIfront, checkSetsAreTheSame, computeCourantTimeGasSource,  hydroStep2dCylindrical, hydroStep2dCylindrical_amr, &
         computeCourantV, writePosRhoPressureVel, writePosRhoPressureVelZERO, killZero, hydrostep2d, checkBoundaryPartners, &
         hydrostep1d, setupAlphaViscosity, sendSinksToZerothThread, computePressureGeneral, hydrostep1dspherical, &
         imposeazimuthalvelocity, forcegascourant, imposefontvelocity, imposekeplerianvelocity, allocatehydrodynamicsAttributes, &
         broadCastSinks, setEquationOfState
    use nbody_mod, only : zerosourcepotential

    use dimensionality_mod, only: setCodeUnit
    use inputs_mod, only: timeUnit, massUnit, lengthUnit, readLucy, checkForPhoto, severeDamping, radiationPressure
    use inputs_mod, only: singleMegaPhoto, stellarwinds, useTensorViscosity, hosokawaTracks, startFromNeutral
    use inputs_mod, only: densitySpectrum, cflNumber, useionparam, xrayonly, isothermal, supernovae, &
         mstarburst, burstTime, starburst, inputseed, doSelfGrav, redoGravOnRead, nHydroperPhoto
    use parallel_mod, only: torus_abort
    use mpi
    integer :: nMuMie
    type(PHASEMATRIX) :: miePhase(:,:,:)
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nLambda
    integer, parameter :: tagsne = 123
    real :: lamArray(:)
    character(len=80) :: mpiFilename, datFilename, mpifilenameB
    real(double) :: dt, cfl, gamma, mu
    integer :: iUnrefine
    integer :: ierr
    real(double) ::  nextDumpTime
    type(VECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: thread1(5120), thread2(5120), nBound(5120), nPairs
    integer :: group(2000), nGroup
!    logical :: globalConverged(512), tConverged(512)
    logical :: dumpThisTime
    real(double) :: deltaTforDump, timeOfNextDump, loopLimitTime
    integer :: iRefine, nUnrefine
    logical :: photoLoop, photoLoopGlobal=.false.
    integer :: i, status(MPI_STATUS_SIZE), tag=30, sign
    integer :: stageCounter=1,  nPhase, nstep, nPhotoIter
    real(double) :: timeSinceLastRecomb=0.d0
    real(double) :: radDt, pressureDt, gravityDt, sourcesourceDt, gasSourceDt, gasDt, tempDouble, viscDt
    real(double) :: vBulk, vSound, recombinationDt, ionizationDt, thermalDt, fractionOfAccretionLum
    logical :: noPhoto=.false., tmpCylindricalHydro, refinedSomeCells, sneAdded, lArray(1)
    integer :: evenUpArray(nHydroThreadsGlobal)
    real :: iterTime(3)
    integer :: iterStack(3), itemp
    integer :: optID
    logical, save :: firstWN=.true.
    integer :: niter, nHydroCounter
    real(double) :: epsoverdeltat, totalMass, tauSca, tauAbs


    nPhotoIter = 1
!    real :: gridToVtu_value
!    integer :: gridVtuCounter
!    gridVtuCounter = 0
    dumpThisTime = .false.
    i = nsource


    call torus_mpi_barrier


    if (nbodyPhysics) then
       if (writeoutput) then
          open(57, file="pos.dat", status="unknown", form="formatted")
          close(57)
       endif
    endif



    direction = VECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    cfl = cflNumber
    globalEpsOverDeltaT = 0.d0
    mu = 1.d0

    sign = 1
    viewVec = VECTOR(-1.d0,0.d0,0.d0)
    viewVec = rotateZ(viewVec, 40.d0*degtorad)
    viewVec = rotateY(viewVec, 30.d0*degtorad)

    if ((globalnsource > 0).and.writeoutput) &
         write(*,*) "Ionizing photons per second ",ionizingFlux(source(1))

    call setCodeUnit(time=timeUnit)
    call setCodeUnit(mass=massUnit)
    call setCodeUnit(length=lengthUnit)
    call writeCodeUnits()

!    tdump = returnCodeUnitTime(tdump)



!    deltaTforDump = 3.14d10 !1kyr
!    if (grid%geometry == "hii_test") deltaTforDump = 2.d10

    fractionOfAccretionLum = 0.d0

    if(.not. readGrid) then
       grid%currentTime = 0.d0
       grid%iDump = -1
       deltaTforDump = tdump

!       if (grid%geometry == "hii_test") deltaTforDump = 9.d10
       if (grid%geometry == "hii_test") deltaTforDump = 1.d9
       if(grid%geometry == "bonnor") deltaTforDump = (1.57d11)!/5.d0 !5kyr
!       if(grid%geometry == "radcloud") deltaTforDump = (1.57d11)!/5.d0 !5kyr
       if(grid%geometry == "starburst") deltaTforDump = tdump
       if(grid%geometry == "sphere") deltaTforDump = tdump
       if(grid%geometry == "SB_WNHII") deltaTforDump = tdump
       if(grid%geometry == "SB_offCentre") deltaTforDump = tdump
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

!       if (grid%geometry(1:6) == "sphere") then
!          call emptyDustCavity(grid%octreeRoot, VECTOR(0.d0, 0.d0, 0.d0), 1400.d0*autocm/1.d10)
!       endif

       ! send grid%currentTime to load balancing threads (and others)
       call MPI_BCAST(grid%currentTime, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
       nextDumpTime = grid%currentTime + deltaTforDump
       timeofNextDump = nextDumpTime



       if(justDump) then 


          if(grid%geometry == "lexington" .or. grid%geometry == "lexpdr") then
             niter = 0
             epsoverdeltat = 1.d0
             call dumpLexingtonMPI(grid, epsoverdeltat, niter)
             call torus_abort("lexington.dat dump completed. Aborting...")
          end if

          if(densitySpectrum) then
  !           print *, "DUMPING DENSITY SPECTRUM"
             if(myRankGlobal == 0) then
                write(mpiFilename,'(a, i4.4, a)') "rhoSpectrum.dat"
                write(mpiFilenameB,'(a, i4.4, a)') "io_neut_mass.dat"
!                print *, "RANK ", myrankglobal, " GOING IN"
                call dumpDensitySpectrumZero(mpiFilename, mpifilenameB, grid%currentTime)
             else
 !               print *, "RANK ", myrankglobal, " GOING IN"
                call dumpDensitySpectrumOther(grid%octreeRoot)
                print *, "RANK ", myrankglobal, "IS OUT"
                call MPI_BARRIER(amrCommunicator, ierr)
                if(myRankGlobal == 1) then
                   call killDensitySpectrumDumper()
                   print *, "killing zero"
                end if
             end if
    !         print *, "RANK ", myrankglobal, "IS AT FINAL GATE"
             call MPI_BARRIER(MPI_COMM_WORLD, ierr)
             call torus_abort("Density spectrum dump completed. Aborting...")
          else
             write(mpiFilename,'(a, i4.4, a)') "quickDump.vtk"
             call writeVtkFile(grid, mpiFilename, &
                  valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
                  "hydrovelocity","sourceCont   ","pressure     ","radmom       ", "rhou         " &
                  , "rhov         ", "rhow         "/))
             
             call MPI_BARRIER(MPI_COMM_WORLD, ierr)
             call torus_abort("vtk dump completed. Aborting...")
          end if

       end if
#ifdef PDR       
       if(dumpBisbas) then
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          if(myrankglobal /= 0) then
             print *, "rank ", myrankglobal, "dumping grid to bisbas"
             call sendGridBisbas(grid%octreeroot, grid)
             call MPI_BARRIER(amrCommunicator, ierr)
             call terminateBisbas()
          else
          print *, "rank ", myrankglobal, "On writing duties"
             call writeGridToBisbas()
!             call writeGridToBisbas(grid)
          end if
          print *, "rank ", myrankglobal, "is done dumping bisbas"
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call torus_abort("3D-PDR grid dump completed. Aborting...")
       end if
#endif
!turn everything on and run a long calculation
       if(singleMegaPhoto) then
          maxPhotoionIter = 200

       else if(UV_vector) then
          maxPhotoionIter = 12
       end if
    end if

!    call writeVtkFile(grid, "ini.vtk", &
!         valueTypeString=(/"rho        ","HI         " ,"temperature", "sourceCont ", "mpithread  " /))

    iunrefine = 0


    if (readlucy) then
       write(mpiFilename,'(a, i4.4, a)') "dump_", iDump, ".grid"
       write(*,*) "Reading from: ",trim(mpiFilename)
       call readAmrGrid(mpiFilename,.false.,grid)
    endif

!    call writeVtkFile(grid, "preghosts.vtk", &
!         valueTypeString=(/"ghosts     ", "mpithread  " /))



    if (.not.loadBalancingThreadGlobal) call allocateHydrodynamicsAttributes(grid%octreeRoot)

    call setupKappaParrays(grid)


!    call writeVtkFile(grid, "inighosts.vtk", &
!         valueTypeString=(/"ghosts     " /))

    if (myRankWorldGlobal == 1) write(*,*) "CFL set to ", cfl

    if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then

       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)
!       do i = 1, nPairs
!          if (myrankWorldglobal==1) &
!               write(*,'(a,i4,i4,a,i4,a,i4,a,i4)') "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i), " group ", group(i)
!       enddo
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

    if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
       call writeInfo("Setting up even up array", TRIVIAL)
       call setupEvenUpArray(grid, evenUpArray)
!       call autoSetupEvenupArray(grid, evenUpArray)
       call writeInfo("Done", TRIVIAL)
!       call evenUpGridMPI(grid, .false., .true., evenuparray)
       call evenUpGridMPI(grid, .true.,dorefine, evenUpArray)
    endif


       call resetnh(grid%octreeRoot)


!       do i = 1, nPairs
!          if (myrankWorldglobal==1)write(*,*) "pair ", i, thread1(i), " -> ", thread2(i), " bound ", nbound(i)
!       enddo


       if(grid%currentTime == 0.d0) then
          if(grid%geometry == "RHDDisc") then
             call imposeazimuthalvelocity(grid%octreeroot)
          elseif(grid%geometry == "fontdisc") then
             call imposefontvelocity(grid%octreeroot)
          elseif(grid%geometry == "simpledisc") then
             call imposekeplerianvelocity(grid%octreeroot)
          endif
          direction = VECTOR(1.d0, 0.d0, 0.d0)
          call calculateRhoU(grid%octreeRoot, direction)
          direction = VECTOR(0.d0, 1.d0, 0.d0)
!          if (.not.cylindricalHydro .or. .not. grid%octreeroot%twod) then
          if(grid%octreeRoot%threeD) then
             call calculateRhoV(grid%octreeRoot, direction)
          endif
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
          
          if(dorefine) then
             call evenUpGridMPI(grid,.false.,.true., evenuparray)      
             call setAllUnchanged(grid%octreeRoot)
             call refineGridGeneric(grid, amrtolerance, evenuparray)
             call writeInfo("Evening up grid", TRIVIAL)    
             call evenUpGridMPI(grid, .false.,.true., evenuparray)
          endif
          call resetnh(grid%octreeRoot)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
          direction = VECTOR(1.d0, 0.d0, 0.d0)
!          print *, "SETUP X TEST A "
          call setupX(grid%octreeRoot, grid, direction)
!          direction = VECTOR(0.d0, 0.d0, 1.d0)
!          print *, "SETUP X TEST B "
 !         call setupX(grid%octreeRoot, grid, direction)
          call setupX(grid%octreeRoot, grid, direction)
          call setupQX(grid%octreeRoot, grid, direction)
          !    call calculateEnergy(grid%octreeRoot, gam

          call writeInfo("Checking Boundary Partner Vectors", TRIVIAL)
          call checkBoundaryPartners(grid%octreeRoot, grid)
          call writeInfo("Initial Boundary Partner Check Passed", TRIVIAL)

       endif
    end if

    call tauAlongPathMPI(grid, VECTOR(0.d0, 0.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), tauAbs, tauSca)
    if (writeoutput) then
       write(*,*) "Tau abs across grid ",tauAbs
       write(*,*) "Tau sca across grid ",tauSca
    endif



    if(grid%currentTime == 0.d0 .and. .not. readGrid .or. singleMegaPhoto .or. UV_vector) then
!       if (startFromNeutral) then
!          call neutralGrid(grid%octreeRoot) 
!       else
!          call ionizeGrid(grid%octreeRoot)
!       endif
       if (.not.timeDependentRT .and. .not. uv_vector .and. .not. useionparam .and. .not.loadBalancingThreadGlobal) then
          if (startFromNeutral) then
             call neutralGrid(grid%octreeRoot)
          else
             call ionizeGrid(grid%octreeRoot)
          endif
       endif

!       if (grid%geometry(1:6) == "sphere") &
!            call emptyDustCavity(grid%octreeRoot, VECTOR(0.d0, 0.d0, 0.d0), 1400.d0*autocm/1.d10)



       if(.not. noPhoto) then

          if ((myrankglobal==1)) then
             !                   firsttime = .false.
             open(67,file="spectrum.dat",status="replace",form="formatted")
             !       write(67,*) "% ",thisOctal%temperature(subcell)
             do i = 1, Source(1)%spectrum%nlambda
                write(67,*) (cspeed/Source(1)%spectrum%lambda(i))*1.e8, Source(1)%spectrum%flux(i)
             enddo
             close(67)
             open(68,file="lamspectrum.dat",status="replace",form="formatted")
             !       write(67,*) "% ",thisOctal%temperature(subcell)
             do i = 1, Source(1)%spectrum%nlambda
                write(68,*) (Source(1)%spectrum%lambda(i)), Source(1)%spectrum%flux(i)
             enddo
             close(68)
          endif


          looplimitTime = deltaTForDump
          looplimittime = 1.d30
          iterTime = 1.e30
          do irefine = 1, 1
             if (irefine == 1) then
                call writeInfo("Calling photoionization loop",TRIVIAL)
!                call setupNeighbourPointers(grid, grid%octreeRoot)
                call resetnh(grid%octreeRoot)
                if (nbodyPhysics.and.hosokawaTracks) then
                   call  setSourceArrayProperties(globalsourceArray, globalnSource, fractionOfAccretionLum)
                endif
                tmpcylindricalhydro=cylindricalhydro
                cylindricalHydro = .false.
                if (.not.timedependentRT .and. .not. xrayonly .and. (.not.isothermal)) &
                     call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, lamArray, &
                     maxPhotoionIter, &
                     loopLimitTime, &
                     looplimittime, .false.,iterTime,.true., evenuparray, optID, iterStack, miePhase, nMuMie)

                if (isoThermal) then
                   call neutralGrid(grid%octreeRoot)
                   call setPhotoionIsothermal(grid%octreeRoot)
                endif
                cylindricalHydro = tmpCylindricalHydro
                call writeInfo("Done",TRIVIAL)
                if(useionparam) then
                   call writeInfo("Calling x-ray step with ionization parameter",TRIVIAL)
                   call simpleXRAY(grid, source(1))
                   call writeInfo("Done",TRIVIAL)
                end if
#ifdef PDR
                if(pdrcalc)then
                   call pdr_main(grid, globalsourcearray, globalnsource)
                end if
#endif
             else
                call writeInfo("Calling photoionization loop",TRIVIAL)
                call setupNeighbourPointers(grid, grid%octreeRoot)
                if (nbodyPhysics.and.hosokawaTracks) then
                   call  setSourceArrayProperties(globalsourceArray, globalnSource, fractionOfAccretionLum)
                endif
                tmpcylindricalhydro=cylindricalhydro
                cylindricalHydro = .false.
                if (.not.timeDependentRT .and. .not. xrayonly .and. (.not. isothermal)) &
                     call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, lamArray, &
                     maxPhotoionIter, loopLimitTime, &
                     looplimittime, timeDependentRT,iterTime,.false., evenuparray, optID, iterStack, miePhase, nMuMie)

                if (isoThermal) then
                   call neutralGrid(grid%octreeRoot)
                   call setPhotoionIsothermal(grid%octreeRoot)
                endif

                cylindricalHydro = tmpCylindricalHydro
                call writeInfo("Done",TRIVIAL)
             endif
             
!             call writeInfo("Dumping post-photoionization data", TRIVIAL)
!             call writeVtkFile(grid, "start.vtk", &
!             valueTypeString=(/"rho        ","HI         " ,"temperature", "ghosts     " /))
!             if(grid%octreeroot%oned) then
!                write(datFilename, '(a, i4.4, a)') "start.dat"
!                call dumpValuesAlongLine(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
!                     VECTOR(3.86d8, 0.d0, 0.d0), 1000)
!                
!             end if
          end do
       end if


       if(singleMegaPhoto .or. UV_vector) then
          write(mpiFilename,'(a, i4.4, a)') "singlemegaphot.grid"
          call writeAmrGrid(mpiFilename, .false., grid)
          call torus_abort("Finished single photo calculation, ending...")
       end if

!       write(mpiFilename,'(a, i4.4, a)') "start.grid"
!       call writeAmrGrid(mpiFilename, .false., grid)

       
       if (.not.LoadBalancingThreadGlobal) then
          call calculateEnergyFromTemperature(grid%octreeRoot)
          call calculateRhoE(grid%octreeRoot, direction)
       endif
       
                    
          if ((myrankGlobal/=0).and.(.not.loadBalancingThreadGlobal)) then
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
             
             if(dorefine) then
                call evenUpGridMPI(grid,.false.,.true., evenuparray)      
                call setAllUnchanged(grid%octreeRoot)
                call refineGridGeneric(grid, amrtolerance, evenuparray)
                call writeInfo("Evening up grid", TRIVIAL)    
                call evenUpGridMPI(grid, .false.,.true., evenuparray)
             endif
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
!                if ((globalnSource > 0).and.(.not.cylindricalHydro)) then
                if ((globalnSource > 0)) then
                   call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, &
                        smallestCellSize/2.d0)
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
       
       if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
          call calculateEnergyFromTemperature(grid%octreeRoot)
          call calculateRhoE(grid%octreeRoot, direction)          
       endif
    end if


    if ((myrankGlobal /= 0).and.cylindricalHydro.and.(.not.loadBalancingThreadGlobal)) then
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       if (myrankWorldglobal == 1) call tune(6,"Alpha viscosity")
       call setupAlphaViscosity(grid, alphaViscosity, 0.1d0)
       if (myrankWorldglobal == 1) call tune(6,"Alpha viscosity")
       call setupCylindricalViscosity(grid%octreeRoot, grid)
    endif

    
    if(grid%geometry == "SB_instblt") then
       print *, "PERTURBING I FRONT "
       call perturbIfront(grid%octreeRoot, grid)
    end if
      
      if(grid%geometry == "hii_test") then
         tEnd = 3.14d13         !1x10^6 year
      else if(grid%geometry == "bonnor") then! .or. grid%geometry=="radcloud") then
         tEnd = 200.d0*3.14d10 !200kyr 
      end if

    nPhase = 1

    nstep = 0
    
    
    if ((.not.loadBalancingThreadglobal).and.(myrankGlobal/=0).and.redoGravOnRead.and.readGrid) then
       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
       call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.)
       call zeroSourcepotential(grid%octreeRoot)
       if (globalnSource > 0) then
          call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
       endif
       call sumGasStarGravity(grid%octreeRoot)
       if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
    endif
    

!    !Thaw - trace courant time history                                                                                                                              
!    open (444, file="tcHistory.dat", status="unknown")
    nHydroCounter = 0

    do while(grid%currentTime < tEnd)
       call MPI_BCAST(grid%currentTime, 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
       nstep = nstep + 1

       
!       tc = 0.d0
!       if (myrankGlobal /= 0) then
!          tc(myrankGlobal) = 1.d30
!          call computeCourantTime(grid, grid%octreeRoot, tc(myRankGlobal))
!          if (nbodyPhysics) call computeCourantTimeNbody(grid, globalnSource, globalsourceArray, tc(myrankGlobal))
!          if (nbodyPhysics) call computeCourantTimeGasSource(grid, grid%octreeRoot, globalnsource, &
!               globalsourceArray, tc(myrankGlobal))
!          raddt = 1.d30
!          call radpressureTimeStep(grid%octreeRoot, raddt)
!          tc(myRankGlobal) = min(tc(myrankGlobal), raddt)
!          call pressureGradientTimeStep(grid, dt)
!          tc(myRankGlobal) = min(tc(myrankGlobal), dt)
!       endif
!       call MPI_ALLREDUCE(tc, tempTc, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, localWorldCommunicator, ierr)
!       dt = MINVAL(temptc(1:nHydroThreadsGlobal))
!       dt = dt * dble(cfl)

       gasDt = 1.d30
       sourcesourcedt = 1.d30
       gasSourcedt = 1.d30
       raddt = 1.d30
       pressuredt = 1.d30
       gravityDt = 1.d30
       recombinationDt = 1.d30
       ionizationDt = 1.d30
       thermalDt = 1.d30
       viscDt = 1.d30
       vBulk = 0.d0
       vSound = 0.d0
       if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
          call computeCourantTime(grid, grid%octreeRoot, gasDt)
          call computeCourantV(grid, grid%octreeRoot, vBulk, vSound)
          if (nbodyPhysics) call computeCourantTimeNbody(grid, globalnSource, globalsourceArray, sourceSourceDt)
          if (nbodyPhysics) call computeCourantTimeGasSource(grid, grid%octreeRoot, globalnsource, &
               globalsourceArray, gasSourceDt)
          if (radiationPressure) call radpressureTimeStep(grid%octreeRoot, raddt)
          call exchangeacrossmpiboundary(grid, npairs, thread1, thread2, nbound, group, ngroup)
          call pressureGradientTimeStep(grid, pressureDt, nPairs, thread1, thread2, nBound, group, nGroup)
          call computeGravityTimeStep(grid, gravityDt, nPairs, thread1, thread2, nBound, group, nGroup)
          if (timeDependentRT) then
             call getRecombinationTime(grid%octreeRoot, recombinationDt)
             call getIonizationTime(grid%octreeRoot, grid, ionizationDt, globalepsOverDeltaT)
             if(.not. cart2d) then
                call getThermalTime(grid%octreeRoot, grid, thermalDt, globalepsOverDeltaT)
             end if
          endif


          if (useTensorViscosity) then
             call viscousTimescale(grid%octreeRoot, grid, viscDt)
          endif
!          if (cylindricalHydro) then
!             call viscousTimescaleCylindrical(grid%octreeRoot, grid, viscDt)
!          endif


       endif

       call MPI_ALLREDUCE(vBulk, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MAX, localWorldCommunicator, ierr)
       vBulk = tempDouble
       call MPI_ALLREDUCE(vSound, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MAX, localWorldCommunicator, ierr)
       vSound = tempDouble


       call MPI_ALLREDUCE(recombinationDt, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       recombinationDt = tempDouble

       call MPI_ALLREDUCE(ionizationDt, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       ionizationDt = tempDouble

       call MPI_ALLREDUCE(thermalDt, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       thermalDt = tempDouble

       call MPI_ALLREDUCE(gasdt, tempdouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       gasDt = tempDouble
       call MPI_ALLREDUCE(sourceSourcedt, tempdouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       sourceSourcedt = tempDouble
       call MPI_ALLREDUCE(gasSourcedt, tempdouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       gasSourceDt = tempDouble
       call MPI_ALLREDUCE(raddt, tempdouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       raddt = tempDouble
       call MPI_ALLREDUCE(pressuredt, tempdouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       pressureDt = tempDouble
       call MPI_ALLREDUCE(gravityDt, tempdouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       gravityDt = tempDouble
       call MPI_ALLREDUCE(viscdt, tempdouble, 1, MPI_DOUBLE_PRECISION, MPI_MIN, localWorldCommunicator, ierr)
       viscDt = tempDouble
       dt = MIN(gasDt, sourceSourceDt, gasSourcedt, raddt, pressureDt, gravityDt, viscDt, &
            thermalDt)

!       dt = MIN(gasDt, sourceSourceDt, gasSourcedt)

       if (writeoutput) then
          write(*,"(a30,1p,e12.3)") "Courant Time: ", dt* dble(cfl)
          write(*,"(a30,1p,e12.3)") "Gas courant Time: ", gasdt
          write(*,"(a30,1p,e12.3)") "Recombination time: ", recombinationDt
          write(*,"(a30,1p,e12.3)") "Ionization time: ", ionizationDt
          write(*,"(a30,1p,e12.3)") "Thermal time: ", thermalDt
          write(*,"(a30,1p,e12.3)") "Source-source courant Time: ", sourceSourceDT
          write(*,"(a30,1p,e12.3)") "Gas-source courant Time: ", gasSourceDt
          write(*,"(a30,1p,e12.3)") "Radiation pressure courant time: ", raddt
          write(*,"(a30,1p,e12.3)") "Pressure gradient courant time: ", pressuredt
          write(*,"(a30,1p,e12.3)") "Grav force courant time: ", gravityDt 
          write(*,"(a30,1p,e12.3)") "Max bulk velocity: ", vBulk/1.d5
          write(*,"(a30,1p,e12.3)") "Max sound speed: ", vSound/1.d5
          write(*,"(a30,1p,e12.3)") "Viscous time: ", viscdt
       endif
       !THAW
       if(forceGasCourant) then
          print *, "forcing gas courant"
          dt = pressuredt
       endif
       dt = dt * dble(cfl)


!       call findMassOverAllThreads(grid, totalMass)
!       if (writeoutput) write(*,*) "Total mass on grid 1  is ", totalmass/mSol


!       write(444, *) jt, MINVAL(tc(1:nHydroThreads)), dt
       
!       if (nstep < 3) then
!          dt = dt * 0.01d0
!       endif
  
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       if (myrankWorldGlobal == 1) then
!          write(*,*) myHydroSetGlobal,temptc(1:nHydroThreadsGlobal)
          write(*,*) "courantTime", dt
       endif



       dumpThisTime = .false.

!       if (.not.loadbalancingthreadglobal) then 
          if ((grid%currentTime + dt) >= timeOfNextDump) then
             dt =  timeofNextDump - grid%currentTime
             dumpThisTime = .true.
          endif
!       endif
       

       fractionOfAccretionlum = 1.d0 !min(1.d0, grid%currentTime / 3.1d10)



       if (myrankWorldGlobal == 1) write(*,*) "dump ",dumpThisTime, " current ", &
            grid%currentTime, " deltaTfordump ",deltaTforDump, " dt ", dt

       if (myrankWorldGlobal == 1) write(*,*) "Time step", dt
       if (myRankWorldGlobal == 1) write(*,*) "percent to next dump ",100.*(timeofNextDump-grid%currentTime)/deltaTforDump, &
            " next dump ", (grid%iDump + 1)

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       if (.not.loadBalancingThreadGlobal) then
          if(checkForPhoto) then
             photoLoopGlobal = .false.
             photoLoop = .false.
             if(myRankGlobal /= 0) then
                ! call checkForPhotoLoop(grid, grid%octreeRoot, photoLoop, dt)
                call advancedCheckForPhotoLoop(grid, grid%octreeRoot, photoLoop, dt, timeSinceLastRecomb)
                call MPI_SEND(photoLoop, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator,  ierr)
                call MPI_RECV(photoLoopGlobal, 1, MPI_LOGICAL, 0, tag, localWorldCommunicator, status, ierr)
             else
                do i = 1, nHydroThreadsGlobal
                   photoLoop = .false.
                   call MPI_RECV(photoLoop, 1, MPI_LOGICAL, i, tag, localWorldCommunicator, status, ierr)
                   
                   if (photoLoop) then
                      photoLoopGlobal = .true.
                   else if(.not. photoLoopGlobal) then
                      photoLoopGlobal = .false.
                   end if
                end do
            
                do i = 1, nHydroThreadsGlobal
                   call MPI_SEND(photoLoopGlobal, 1, MPI_LOGICAL, i, tag, localWorldCommunicator,  ierr)
                end do
             end if
             
             
          else
             photoLoopGlobal = .true.
          end if
       else
          photoLoopGlobal = .true.
       endif
       if(photoLoopGlobal .and. .not. noPhoto) then
          call writeInfo("Calling photoionization loop A",TRIVIAL)
!          call ionizeGrid(grid%octreeRoot)
          if(dt /= 0.d0) then
             loopLimitTime = grid%currentTime+dt
          else
             looplimittime = deltaTForDump
          end if
          if (timeDependentRT) loopLimitTime = dt
!         write(*,*) myrankGlobal, " calling setup neighbour pointers"
!          call setupNeighbourPointers(grid, grid%octreeRoot)
!          write(*,*) myrankGlobal, " calling resetnh"
          call resetnh(grid%octreeRoot)
!          if (grid%geometry == "sphere") &
!               call emptyDustCavity(grid%octreeRoot, VECTOR(0.d0, 0.d0, 0.d0), 1400.d0*autocm/1.d10)

          if (nbodyPhysics.and.hosokawaTracks) then
             call  setSourceArrayProperties(globalsourceArray, globalnSource, fractionOfAccretionLum)
          endif
!          call photoIonizationloopAMR(grid, source, nSource, nLambda,lamArray, 1, loopLimitTime, loopLimitTime, .false., iterTime, &
!               .true., evenuparray, sign)
                tmpcylindricalhydro=cylindricalhydro
                cylindricalHydro = .false.
!                nPhotoIter = 10
!                nPhotoIter = int(10 - grid%idump)
!                nphotoIter = max(1, nPhotoIter)
                if(.not. xrayonly) then
                   if (mod(nHydroCounter,nHydroPerPhoto) == 0) then
                      call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, &
                           lamArray, nPhotoIter, loopLimitTime, &
                           looplimittime, timeDependentRT,iterTime,.true., evenuparray, optID, iterStack, miePhase, nMuMie) 
                   endif
                endif

                if (isoThermal) then
                   call neutralGrid(grid%octreeRoot)
                   call setPhotoionIsothermal(grid%octreeRoot)
                endif

                cylindricalHydro = tmpCylindricalHydro

         call writeInfo("Done",TRIVIAL)
         
         if(useionparam) then
            call writeInfo("Calling x-ray step with ionization parameter",TRIVIAL)
            call simpleXRAY(grid, source(1))
            call writeInfo("Done",TRIVIAL)
         end if
#ifdef PDR
         if(pdrcalc)then
            call pdr_main(grid, globalsourcearray, globalnsource)
         end if
#endif
!         if(useionparam) then
!            call writeInfo("Calling x-ray step with ionization parameter",TRIVIAL)
!            
!            call writeInfo("Done",TRIVIAL)
!         end if
         
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

          if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
             call calculateEnergyFromTemperature(grid%octreeRoot)
             call calculateRhoE(grid%octreeRoot, direction)
          endif


!          call writeVtkFile(grid, "afterloop.vtk", &
!               valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!               "hydrovelocity","sourceCont   ","pressure     ","radmom       ",     "radforce     ", &
!               "diff         ","dust         ","u_i          ",  &
!               "phi          ","rhou         ","rhov         ","rhow         ","rhoe         ", &
!               "vphi         ","jnu          ","mu           ", &
!               "fvisc1       ","fvisc2       ","fvisc3       ","crossings    "/))


       if (myRankWorldGlobal == 1) call tune(6,"Hydrodynamics step")
       nHydroCounter = nHydroCounter + 1
       currentlyDoingHydroStep = .true.


       call writeInfo("calling hydro step",TRIVIAL)
       
!       if (myrankGlobal /= 0) call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

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


       if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
          if (cylindricalHydro) then
             call writeInfo("Calling cylindrical coord hydro step", TRIVIAL)
!             call hydroStep2dCylindrical(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
             call hydroStep2dCylindrical_amr(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
          else
             if(grid%currentTime==0.d0 .and. grid%geometry == "SB_WNHII")then
                if(grid%octreeRoot%twoD) then
                   call writeInfo("calling hydro step 2D",TRIVIAL)
                   if(firstWN) then
                      call hydroStep2d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup, &
                           perturbPressure=.true.)
                   else
                      call hydroStep2d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup, &
                           perturbPressure=.false.)
                   end if
                elseif(grid%octreeRoot%threeD) then
                   call writeInfo("calling hydro step 2D",TRIVIAL)
                   if(firstWN) then
                      call hydroStep3d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup,doSelfGrav=doselfGrav &
                           , perturbPressure=.true.)
                   else
                      call hydroStep3d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup,doSelfGrav=doselfGrav &
                           , perturbPressure=.false.)
                   end if
                   
                end if
                
             else if (grid%octreeRoot%threeD) then

!                call writeVtkFile(grid, "beforehydro.vtk", &
!                     valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!                     "hydrovelocity","sourceCont   ","pressure     ","radmom       ",     "radforce     ", &
!                     "diff         ","dust         ","u_i          ",  &
!                     "phi          ","rhou         ","rhov         ","rhow         ","rhoe         ", &
!                     "vphi         ","jnu          ","mu           ", &
!                     "fvisc1       ","fvisc2       ","fvisc3       ","crossings    "/))


!                call hydroStep3d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup,doSelfGrav=doselfGrav &
!                     , perturbPressure=.false.)

                call hydroStep3d_amr(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup,doSelfGrav=doselfGrav &
                     , perturbPressure=.false.)



!       call findMassOverAllThreads(grid, totalMass)
!       if (writeoutput) write(*,*) "Total mass on grid 2  is ", totalmass/mSol

             else if (grid%octreeroot%twod) then
                call hydroStep2d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup, &
                     perturbPressure=.false.)
             else
                if(spherical) then
                   call hydroStep1dspherical(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
                else
                   call hydroStep1d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup)
                end if
                !                call torus_abort("1d radhydro not supported")
             end if
             firstWN = .false.
             !          else
             !             call hydroStep3d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup,doSelfGrav=doselfGrav)
          end if
          if(grid%geometry == "RHDDisc") then
             call imposeazimuthalvelocity(grid%octreeroot)
          elseif(grid%geometry == "fontdisc") then
             call imposefontvelocity(grid%octreeroot)
             if (forceMinRho) then
                call enforceminrhofont(grid%octreeroot) 
             endif
          elseif(grid%geometry == "simpledisc") then
             call imposekeplerianvelocity(grid%octreeroot)
             if (forceMinRho) then
                call enforceminrhodisc(grid%octreeroot) 
             endif
          endif
       endif
!    endif

    currentlyDoingHydroStep = .false.

       if (nHydroSetsGlobal > 1) call checkSetsAreTheSame(grid%octreeRoot)


          if ((myRankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
             if(dorefine) then
                call setAllUnchanged(grid%octreeRoot)
                !             if (myRankWorldGlobal == 1) call tune(6,"Even up grid")
                !             call evenUpGridMPI(grid,.false.,.true., evenuparray)
                !             if (myRankWorldGlobal == 1) call tune(6,"Even up grid")
                if (myRankWorldGlobal == 1) call tune(6,"Refine grid")
                call refineGridGeneric(grid, amrtolerance, evenuparray, refinedSomeCells=refinedSomeCells)
                if (myRankWorldGlobal == 1) call tune(6,"Refine grid")
             endif


             if (refinedSomeCells) then
                if (myRankWorldGlobal == 1) call tune(6,"Even up grid")
                call writeInfo("Evening up grid", TRIVIAL)
                call evenUpGridMPI(grid, .false.,.true., evenuparray)
                if (myRankWorldGlobal == 1) call tune(6,"Even up grid")

                if (doselfGrav) then
                   if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
                   call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, onlyChanged=.true.)
                   call zeroSourcepotential(grid%octreeRoot)
                   if (globalnSource > 0) then
                      call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
                   endif
                   call sumGasStarGravity(grid%octreeRoot)
                   if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
                endif

             endif
!             call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

          endif



          if (severeDamping) call cutVacuum(grid%octreeRoot)

          if (writeoutput) write(*,*) "starburst ",starburst, "current ",grid%currentTime, &
               " burst time ",burstTime, " nSource ",globalnSource
          if (starburst.and.(grid%currentTime > burstTime).and.(globalnSource == 0)) then
             if (writeoutput) call writeInfo("Starting setting up sources.")
             call randomNumberGenerator(putISeed = inputSeed)
             call randomNumberGenerator(syncIseed=.true.)
             if (associated(globalSourceArray)) then
                deallocate(globalSourceArray)
                globalSourceArray => null()
             endif
             allocate(globalsourcearray(1:1000))
             globalsourceArray(:)%outsideGrid = .false.
             globalnSource = 0
             call createSources(globalnSource, globalSourceArray, "instantaneous", 0.d0, mStarburst, 0.d0)
             if (.not.loadBalancingThreadGlobal) then
                 call putStarsInGridAccordingToDensity(grid, globalnSource, globalsourceArray)
                 call dumpSources(globalsourceArray, globalnSource)

                 call writeVtkFilenBody(globalnSource, globalsourceArray, "starburst_initial.vtk")

                 call writeVtkFile(grid, "rho_at_starburst.vtk", &
                   valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
                   "hydrovelocity","sourceCont   ","pressure     ","radmom       ",     "radforce     ", &
                   "diff         ","dust         ","u_i          ",  &
                   "phi          ","rhou         ","rhov         ","rhow         ","rhoe         ", &
                   "vphi         ","jnu          ","mu           ", &
                   "fvisc1       ","fvisc2       ","fvisc3       ","crossings    "/))


                 if (doselfgrav.and.(myrankGlobal/=0)) then
                    if (myrankWorldglobal == 1) call writeInfo("Solving self gravity...")
                    if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
                    call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.)
                    if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
                    if (myrankWorldglobal == 1) call writeInfo("Done.")
                 endif
                 call randomNumberGenerator(randomSeed = .true.)
             endif
             if (writeoutput) call writeInfo("Setting up sources done.")
             call sendSinksToZerothThread(globalnSource, globalsourceArray)
             call broadcastSinks
             call photoIonizationloopAMR(grid, globalsourceArray, globalnSource, nLambda, &
                  lamArray, 10, loopLimitTime, &
                  looplimittime, timeDependentRT,iterTime,.true., evenuparray, optID, iterStack, miePhase, nMuMie) 


          endif

          if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal).and.stellarwinds) &
               call addStellarWind(grid, globalnSource, globalsourcearray, dt)

          sneAdded = .false.
          if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal).and.supernovae) then
             call addSupernovae(grid, globalnSource, globalsourcearray, dt, sneAdded)
             if (doselfgrav.and.sneAdded) call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup, multigrid=.true.)
             lArray(1) = sneAdded
             if (myrankGlobal==1) call MPI_SEND(lArray, 1, MPI_LOGICAL, 0, tagsne, localWorldCommunicator,  ierr)
          endif
          if ((myrankGlobal == 0).and.(supernovae)) then
             call MPI_RECV(lArray, 1, MPI_LOGICAL, 1, &
               tagsne, localWorldCommunicator, status, ierr)
             sneAdded = lArray(1)
          endif


          if (supernovae.and.(.not.loadBalancingThreadGlobal)) then
             if (sneAdded) then
                dumpThisTime = .true.
                timeOfNextDump = grid%currentTime + dt
                deltaTfordump = 1.d9
             endif
          endif

!add/merge sink particles where necessary
       if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
          if (nbodyPhysics.and.addSinkParticles) call addSinks(grid, globalsourceArray, globalnSource)       

!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
!          call freeglobalsourceArray()
!          globalnSource = 1
!          allocate(globalSourceArray(1:1))
!          globalSourceArray(1)%position = VECTOR(2d7,2d7,2d7)
!          globalSourceArray(1)%mass = 20.d0*msol
!          globalSourceArray(1)%accretionRadius = 1.d10


!          if (myRankWorldGlobal == 1) call tune(6,"Merging sinks")
!          if (nbodyPhysics) call mergeSinks(grid, globalsourceArray, globalnSource)
!          if (myRankWorldGlobal == 1) call tune(6,"Merging sinks")
          do i = 1, globalnSource
             call emptySurface(globalsourceArray(i)%surface)
             call buildSphereNbody(globalsourceArray(i)%position, 2.5d0*smallestCellSize, &
                  globalsourceArray(i)%surface, 20)
          enddo

       endif
       call sendSinksToZerothThread(globalnSource, globalsourceArray)
       call broadcastSinks
       if (nbodyPhysics.and.hosokawaTracks) then
          call  setSourceArrayProperties(globalsourceArray, globalnSource, fractionOfAccretionLum)
       endif


       

!       write(mpiFilename,'(a, i4.4, a)') "postStep.vtk"
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!            "hydrovelocity","sourceCont   ","pressure     "/))

       if (myRankWorldGlobal == 1) call tune(6,"Hydrodynamics step")
       if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) & 
            call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       if (myrankGlobal /= 0) call resetNh(grid%octreeRoot)

       if ((myrankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
          if (dounrefine) then
             nUnrefine = 0
             iUnrefine = iUnrefine + 1
             if (iUnrefine == 5) then
                if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
                call unrefineCells(grid%octreeRoot, grid, nUnrefine,amrUnrefinetolerance)
                call MPI_ALLREDUCE(nUnrefine, itemp, 1, MPI_INTEGER, MPI_SUM, amrCommunicator, ierr)
                nUnrefine = itemp
                if (writeoutput) write(*,*) "Number of unrefined cells = ",nUnrefine
                if (myrankWorldglobal == 1) call tune(6, "Unrefine grid")
                iUnrefine = 0
             endif


             if (nUnrefine > 0) then
                call evenUpGridMPI(grid, .true., .true., evenuparray)
                call allocateHydrodynamicsAttributes(grid%octreeRoot)
                if (doselfGrav) then
                   if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
                   call selfGrav(grid, nPairs, thread1, thread2, nBound, group, nGroup)
                   call zeroSourcepotential(grid%octreeRoot)
                   if (globalnSource > 0) then
                      call applySourcePotential(grid%octreeRoot, globalsourcearray, globalnSource, smallestCellSize)
                   endif
                   call sumGasStarGravity(grid%octreeRoot)
                   if (myrankWorldglobal == 1) call tune(6,"Self-gravity")
                endif
             endif
          endif
       endif

       grid%currentTime = grid%currentTime + dt

       if (nbodyPhysics) globalSourceArray(1:globalnSource)%time = grid%currentTime

       if (.not.loadBalancingThreadGlobal) call findMassOverAllThreads(grid, totalMass)
       if (writeoutput) write(*,*) "Total mass on grid is ", totalmass/mSol

       
       if (myRankWorldGlobal == 1) write(*,*) "Current time: ",grid%currentTime

!Track the evolution of the ionization front with time
!          write(datFilename, '(a, i4.4, a)') "Ifront.dat"
!          call dumpStromgrenRadius(grid, datFileName, VECTOR(-3.3d8,  -3.3d8, 3.3d8), &
!               VECTOR(3.3d8, 3.3d8, 3.3d8), 1000)


!       write(mpiFilename,'(a, i4.4, a)') "dump_", nPhase,".vtk"
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!            "hydrovelocity","sourceCont   ", "pressure     "/))
!       nPhase = nPhase + 1
       if ((.not.loadbalancingthreadglobal).and.(gasdt < 1.d-1)) then
         dumpThisTime = .true.
       endif

       
       if (dumpThisTime) then
 
          timeOfNextDump = timeOfNextDump + deltaTForDump
          grid%iDump = grid%iDump + 1

       endif
       if (dumpThisTime.and.(.not.loadBalancingThreadGlobal)) then
          if(mod(real(grid%iDump), real(vtuToGrid)) == 0.0) then
             write(mpiFilename,'(a, i4.4, a)') "dump_", grid%iDump,".grid"
             call writeAmrGrid(mpiFilename, .false., grid)
          end if


          if(grid%geometry == "SB_WNHII") then
             call dumpWhalenNormanTest(grid)
          end if
!          if(grid%geometry == "SB_offCen") then
          if(0 == 1) then
             if(myRankGlobal == 0) then
                write(mpiFilename,'(a, i4.4, a)') "dump_", grid%iDump,".txt"
                call writePosRhoPressureVelZERO(mpiFilename, grid)
             else
                call writePosRhoPressureVel(grid%octreeRoot)
                call killZero()
             end if
             
             
          end if

          if(grid%geometry == "SB_instblt") then
             call dumpIfrontTest(grid)
          end if

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


          if (myrankglobal /= 0) call computepressureGeneral(grid, grid%octreeRoot, .false.)
          write(mpiFilename,'(a, i4.4, a)') "dump_", grid%iDump,".vtk"
          call writeVtkFile(grid, mpiFilename, &
               valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
               "hydrovelocity","sourceCont   ","pressure     ","radmom       ",     "radforce     ", &
               "diff         ","dust         ","u_i          ",  &
               "phi          ","rhou         ","rhov         ","rhow         ","rhoe         ", &
               "vphi         ","jnu          ","mu           ", &
               "fvisc1       ","fvisc2       ","fvisc3       ","crossings    "/))


          if(densitySpectrum) then
             !           print *, "DUMPING DENSITY SPECTRUM"
             if(myRankGlobal == 0) then
                write(mpiFilename,'(a, i4.4, a)') "rhoSpectrum_",grid%iDump,".dat"
                write(mpiFilenameB,'(a, i4.4, a)') "io_neu_mass.dat"
                !                print *, "RANK ", myrankglobal, " GOING IN"
                call dumpDensitySpectrumZero(mpiFilename, mpiFilenameB, grid%currentTime)
             else
                !               print *, "RANK ", myrankglobal, " GOING IN"
                call dumpDensitySpectrumOther(grid%octreeRoot)
                print *, "RANK ", myrankglobal, "IS OUT"
                call MPI_BARRIER(amrCommunicator, ierr)
                if(myRankGlobal == 1) then
                   call killDensitySpectrumDumper()
                   print *, "killing zero"
                end if
             end if
             !         print *, "RANK ", myrankglobal, "IS AT FINAL GATE"
             call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          end if


          write(mpiFilename,'(a,i4.4,a)') "nbody",grid%iDump,".vtk"
          call writeVtkFilenBody(globalnSource, globalsourceArray, mpiFilename)

          if (writeoutput) then
             write(mpiFilename,'(a,i4.4,a)') "source",grid%idump,".dat"
             globalSourceArray(:)%time = grid%currentTime
             call writeSourceArray(mpifilename)
          endif

          if (spherical) then
             write(mpiFilename,'(a,i4.4,a)') "cells",grid%idump,".dat"
             call write1dlist(grid, mpifilename)
          endif
!          if (.not.CylindricalHydro) then
             write(mpiFilename,'(a,i4.4,a)') "radial",grid%idump,".dat"
             call  dumpValuesAlongLine(grid, mpiFilename, VECTOR(1.d0,0.d0,0.0d0), &
                  VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0),1000)
!          endif
           
          if(grid%geometry == "bonnor" .and. grid%octreeRoot%oneD) then
             write(mpiFilename,'(a,i4.4,a)') "1DRHD_torus",grid%idump,".dat"
             call  dumpValuesAlongLine(grid, mpiFilename, & 
                  VECTOR(-grid%octreeRoot%subcellSize,0.d0,0.0d0), &
                  VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0),1000)

             write(datFilename, '(a, i4.4, a)') "Ifront.dat"     
             call dumpStromgrenRadius(grid, datFileName,  & 
                  VECTOR(-grid%octreeRoot%subcellSize,0.d0,0.0d0), &
                  VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0),1000)
          end if


          
!Track the evolution of the ionization front with time
       if(grid%geometry == "hii_test" .or. grid%geometry == "SB_Dtype") then
          write(datFilename, '(a, i4.4, a)') "hii_test",grid%iDump,".dat"
          call dumpValuesAlongLine(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
               VECTOR(4.d9, 0.d0, 0.d0), 1000)


!          write(datFilename, '(a, i4.4, a)') "Ifront.dat"     
!          call dumpStromgrenRadius(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
!               VECTOR(3.86d8, 0.d0, 0.d0), 1000)!
          write(datFilename, '(a, i4.4, a)') "Ifront.dat"     
          call dumpStromgrenRadius(grid, datFileName, VECTOR(0.d0,  0.d0, 0.d0), &
               VECTOR(4.d9, 0.d0, 0.d0), 1000)!
!          write(datFilename, '(a, i4.4, a)') "hii_test",grid%iDump,".dat"
!          call dumpValuesAlongLine(grid, datFileName, VECTOR(1.75d9,  0.d0, 1.75d9), &
!               VECTOR(3.5d9, 0.d0, 1.75d9), 1000)!
!
!!          write(datFilename, '(a, i4.4, a)') "Ifront.dat"     
!!          call dumpStromgrenRadius(grid, datFileName, VECTOR(1.75d9,  0.d0, 1.75d9), &!
!               VECTOR(3.5d9, 0.0d0, 1.75d9), 1000)!
!
!          write(datFilename, '(a, i4.4, a)') "Ifront.dat"     
!          call dumpStromgrenRadius(grid, datFileName, VECTOR(2d9,  0.d0, 2d9), &
!               VECTOR(4.d9, 0.0d0, 4.d9), 1000)
       end if
 
    endif
    stageCounter = stageCounter + 1
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    close(444)
!    write(*,*) "myRank", myRankGlobal, "finishing loop. Time:", grid%currentTime, "tend ", tend
    if (gasdt < 1.d-1) then
      call torus_stop("tc less than 0.1 s")
    endif

 enddo

end subroutine radiationHydro
#endif

  subroutine photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, maxIter, tLimit, deltaTime, timeDep, iterTime, &
       monteCheck, evenuparray, optID, iterStack,  miePhase, nMuMie, sublimate)
    use inputs_mod, only : quickThermal, inputnMonte, noDiffuseField, minDepthAMR, maxDepthAMR, binPhotons,monochromatic, &
         readGrid, dustOnly, bufferCap, doPhotorefine, doRefine, amrtolerance, hOnly, &
         optimizeStack, stackLimit, dStack, singleMegaPhoto, stackLimitArray, customStacks, tMinGlobal, variableDustSublimation, &
         radPressureTest, justdump, uv_vector, inputEV, xrayCalc, useionparam, dumpregularVTUs


    use inputs_mod, only : usePacketSplitting, inputNSmallPackets, amr2d, amr3d, forceminrho, nDustType, readgrid, &
         loadBalancing, inputSeed, tsub

    use hydrodynamics_mod, only: refinegridgeneric, evenupgridmpi, checkSetsAreTheSame
    use dust_mod, only : sublimateDust, stripDustAway
    use diffusion_mod, only : defineDiffusionOnKappap
    use mpi
    implicit none
    integer :: nMuMie
    type(PHASEMATRIX) :: miePhase(:,:,:)
    integer(bigint) :: nMonte, nThreadMonte, nTotalMonte
    logical, optional :: sublimate
    logical :: doSublimate, timeDep, monteCheck
    real(double) :: tLimit
    real(double) :: deltaTime
    type(GRIDTYPE) :: grid
    ! Commented out variables beforeSubcell and beforeOctal need to be declared OMP private if re-instated
    type(OCTAL), pointer :: thisOctal!, beforeOctal
!    integer :: beforeSubcell
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
    integer :: i
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
    integer :: nIter, nPeriodic, loopSource
    logical :: converged, thisThreadConverged, failed
    real(double) :: photonPacketWeight, spectrumweight
    real(double) :: tPhoton
    real(double) :: albedo
    integer :: maxIter
    type(SAHAMILNETABLE),save :: hTable, heTable
    type(RECOMBTABLE),save :: Hrecombtable
    real(double) :: freq(1006), dfreq(1006), spectrum(1006)
!    real(double) :: freq(1000+(2*grid%nIon)), dfreq(1000+(2*grid%nIon)), spectrum(1000+(2*grid%nIon)) 
    real(double) :: nuStart, nuEnd
    real(double) :: kappaAbsGas, kappaAbsDust, escat 
    real(double) :: r1 ! add this to the private openmp section if reinstated
    integer, parameter :: nTemp = 9
    real(double) :: v, dustHeating
    real(double) :: kappaP
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
    integer, parameter :: tag = 41, finaltag = 55
    integer(bigint) :: nTotScat, nPhot
    integer :: newThread, iSignal
    logical, save :: firstTimeTables = .true.
    integer(bigint) :: nEscaped=0
    integer, parameter :: nTimes = 3
    logical :: photonsStillProcessing
    integer(bigint), allocatable :: nEscapedArray(:)
    integer :: status(MPI_STATUS_SIZE)

    real(double) :: maxDeltaT
    logical :: anyUndersampled
    logical :: escapeCheck
    logical :: sourcePhoton

    integer :: dprCounter
    integer :: recCounter
    integer :: threadCounter
    
    integer :: maxStackLimit
    integer :: sendStackLimit

    integer :: evenuparray(nHydroThreadsGlobal)
    integer :: k
    real(double) :: nuThresh

    !optimization variables
!    integer, parameter :: stackLimit=200
    real, intent(inout) :: iterTime(3)
    integer, intent(inout) :: iterStack(3)
    integer :: ZerothstackLimit
    real :: startTime, endTime, newTime
    integer, intent(inout) :: optID
    integer :: optCounter, optCount
    integer :: dstackNaught
    logical, save :: optConverged=.false.
    logical, save :: iniIterTime=.false.
    logical, save :: iniCustomTime=.true.
    real(double) :: m1, m2
    real :: tauMax
  !  real :: oldTime = 1.e10
  !  integer :: newStackLimit= 0, oldStackLimit= 0
    real(double) :: oldT(8), oldF(8), fracT(8), fracF(8)
    integer :: thisPacket, sendCounter
    integer, allocatable :: nSaved(:)
    integer :: stackSize, p
    integer :: mpi_vector, mpi_photon_stack
    logical :: sendAllPhotons = .false., donePanicking = .true.
    real(double) :: nIonizingPhotons, thisSourceFlux
    logical :: crossedPeriodic = .false.

!    type(PHOTONPACKET) :: photonPacketStack(stackLimit*nThreadsGlobal)
    type(PHOTONPACKET), allocatable :: photonPacketStack(:)
!    type(PHOTONPACKET) :: toSendStack(stackLimit), currentStack(stackLimit)
    type(PHOTONPACKET), allocatable :: toSendStack(:), currentStack(:)

   !Custom MPI type variables
    integer(MPI_ADDRESS_KIND) :: displacement(11)
    integer :: count = 11
    integer :: blocklengths(11) = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
    integer :: oldTypes(11)
    Integer :: iDisp

    !Buffer send variables 
    integer :: bufferSize  !Send buffer size in bytes
    character, allocatable :: buffer(:)
!    logical :: ready = .true.
    real(double) :: fac, luminosity1, luminosity2, luminosity3
    real(double) :: zTemp, thisLum



    !OMP variables
    logical :: finished = .false.
    logical :: voidThread

    logical :: doingSmallPackets
    logical :: startNewSmallPacket
    integer :: nSmallPackets
    logical :: bigPhotonPacket, smallPhotonPacket, lastPhoton
    integer :: iSmallPhotonPacket
    real(double) :: smallPhotonPacketWeight, smallPacketFreq, bigPhotonPacketWeight
    type(VECTOR) :: smallPacketOrigin, uHatDash, zHat
    type(VECTOR) :: uHatBefore, uHatAfter
    

    logical :: firstConverged=.true.


    !Energy cons variables
    real(double) :: totalPower = 0.d0
    real(double) :: lams(1000) = 0.d0
    real(double) :: countArray(1000) = 0.d0
    real(double) :: tempDouble

    real, allocatable :: tempCell(:,:), temp1(:), temp2(:), tempIon(:,:,:)
    real(double), allocatable :: temp1d(:)

    integer :: n_rmdr, mOctal, nNotEscaped
    character(len=80) :: mpiFilename, message
    integer :: ierr
    integer :: iter, nFrac, nNotEndLoop, nToNextEventPhoto
    real :: totFrac
    integer :: maxChild
    integer :: nScatBigPacket, j, nScatSmallPacket
    logical :: sourceInThickCell, tempLogical
    logical :: undersampled, flushBuffer, containsLastPacket, movedCells
    logical, save :: splitThisTime = .false.
    logical, save :: firstWarning = .true.
    logical, save :: firstLoadBalancing = .true.
    real(double) :: maxDiffRadius(1:100)
    real(double) :: maxDiffRadius1(1:100), maxDiffRadius2(1:100)
    real(double) :: maxDiffRadius3(1:100), tauWanted, photonMomentum
    type(VECTOR) ::  vec_tmp, uNew
    integer :: receivedStackSize, nToSend
    integer :: nDomainThreads, localRank, m, nBundles
    real :: FinishTime, WaitingTime, globalStartTime, globalTime
    real(double), allocatable :: efficiencyArray(:)
    real(double) :: medianEfficiency
    !xray stuff
    type(AUGER) :: augerArray(5, 5, 10)
    integer, save :: oldStackLimit = 0
    logical, save :: firstCall = .true.
    !AMR
!    integer :: iUnrefine, nUnrefine


    !!Thaw - optimize stack will be run prior to a big job to ensure that the most efficient stack size is used
    !start with stack size of 1
    !if(optimizeStack) then 
    !   stackLimit = 1
    !   zerothStackLimit = 1
    !end if

    if (firstCall) then
       firstCall = .false.
       oldStackLimit = stackLimit
    endif

    
    if (globalnSource == 0) goto 666

    allocate(efficiencyArray(1:nThreadsGlobal-1))

    if (readGrid) splitThisTime = .true.

    if(justdump .and. grid%geometry == "lexington") then
       niter = 0
       epsoverdeltaT = 0.5d0
      call dumpLexingtonMPI(grid, epsoverdeltat, niter)
   end if

   nDomainThreads = nHydroThreadsGlobal + nLoadBalancingThreadsGlobal
    


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
       maxStackLimit = maxval(stackLimitArray(1:nDomainThreads))
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

    if (forceMinRho) then
       call enforceminrho(grid%octreeroot) 
    endif

    
    allocate(photonPacketStack(maxStackLimit*nDomainThreads))
    allocate(toSendStack(maxStackLimit))
    allocate(currentStack(maxStackLimit))

    !Custom MPI data types for easier send/receiving
    !MPI datatype for out TYPE(VECTOR) variables
    call MPI_TYPE_CONTIGUOUS(3, MPI_DOUBLE_PRECISION, MPI_VECTOR, ierr)
    call MPI_TYPE_COMMIT(MPI_VECTOR, ierr)

    !MPI datatype for the photon_stack data type
    oldTypes = (/ MPI_VECTOR, MPI_VECTOR, MPI_DOUBLE_PRECISION, &
         MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_LOGICAL, &
         MPI_LOGICAL, MPI_LOGICAL, MPI_LOGICAL/)

    call MPI_GET_ADDRESS(toSendStack(1)%rVec, displacement(1), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%uHat, displacement(2), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%freq, displacement(3), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%tPhot, displacement(4), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%ppw, displacement(5), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%destination, displacement(6), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%sourcePhoton, displacement(7), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%crossedPeriodic, displacement(8), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%bigPhotonPacket, displacement(9), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%smallPhotonPacket, displacement(10), ierr)
    call MPI_GET_ADDRESS(toSendStack(1)%lastPhoton, displacement(11), ierr)

    do iDisp = 11, 1, -1
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



    allocate(nEscapedArray(1:nDomainThreads))
    allocate(nSaved(nDomainThreads))
    nescapedArray = 0    
    dprCounter = 0

    
!    if(uv_vector .or. singlemegaphoto .or. dumpRegularVTUS) then
!       write(mpiFilename,'(a, i4.4, a)') "pre", nIter,".vtk"!
!       call writeVtkFile(grid, mpiFilename, &
!            valueTypeString=(/"rho          ", "HI           " , "temperature  ", "uvvec        "/))
!    end if
    



!    if(.not. optimizeStack) then
       photonPacketStack%Freq = 0.d0
       photonPacketStack%Destination = 0
       photonPacketStack%tPhot = 0.d0
       photonPacketStack%crossedPeriodic = .false.
!    end if


       if(monochromatic) then
          nfreq = 2
          if(grid%geometry == "SB_WNHII") then
             nuStart = (18.60001*evtoerg/hcgs)
             nuEnd =  (18.60001*evtoerg/hcgs)
          elseif(uv_vector) then
             nustart = inputEV*evtoerg/hcgs
             nuend = (inputEV+1.d-10)*evtoerg/hcgs
          else
             nuStart = (13.60001*evtoerg/hcgs)
             nuEnd =  (13.60001*evtoerg/hcgs)
          end if
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

       if(xraycalc .and. .not. useionparam) then
          call setUpAugerData(augerArray)
       end if

       firstTimeTables = .false.
    endif

    lCore = 0.d0
    do i = 1, nSource
       !If outside the grid then correct for the flux attenuation due to distance of star from grid
       if(source(i)%outsideGrid) then
          if(cart2d) then             
             if(grid%octreeroot%oned) then
                lCore = lCore + source(i)%luminosity * (1.d20* &
                     (2.d0*grid%halfsmallestsubcell)**2) / (fourPi*source(i)%distance**2)
             else
                lCore = lCore + source(i)%luminosity * (2.d0*grid%octreeRoot%subcellSize*1.d20* &
                     2.d0*grid%halfsmallestsubcell) / (fourPi*source(i)%distance**2)
             end if
          else
             lCore = lCore + source(i)%luminosity * (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 / &
                  (fourPi*source(i)%distance**2)
          end if
       else
          lCore = lCore + source(i)%luminosity
       end if
    enddo

    if (myrankWorldglobal == 1) write(*,'(a,1pe12.5)') "Total source luminosity (lsol): ",lCore/lSol

    nionizingphotons = 0
    if(globalnsource > 1) then
       do loopsource = 1, globalnsource
          thisSourceFlux = ionizingFlux(source(loopSource))
          nIonizingPhotons = nionizingphotons + thisSourceFlux
          if (writeoutput .and. thisSourceFlux > 1.d40) then
            write(*,'(a, i3, 1pe12.1)') "Ionizing photons per sec for source ", loopsource, thisSourceFlux 
          endif 
       enddo
    else
       nIonizingPhotons = ionizingFlux(source(1))
    endif
    
    if (writeoutput) then
!       write(*,'(a,1pe12.1)') "Ionizing photons per second ", ionizingFlux(source(1))
       write(*,'(a,1pe12.1)') "Total ionizing photons per second ", nIonizingPhotons
    endif

    if (inputNMonte < 0) then
       if (writeoutput) write(*,*) "Skipping photoion loop as input nMonte < 0"
       goto 666
    endif

    !nmonte selector: Only works for fixed grids at present
    if (inputnMonte == 0) then
       
       if(.not. monteCheck .or. readGrid) then
          !waymaker photoionization loop
          if(grid%octreeRoot%twoD) then
                nTotalMonte = int(4.d0**(maxDepthAMR))*100
             else if(grid%octreeRoot%threeD) then
                !nMonte = (8.d0**(maxDepthAMR))
                !nMonte = 5242880/2.
!                nMonte = 150000000
                nTotalMonte = 209715200
               ! nMonte = 409600.0
             else
                !nMonte = 2**(maxDepthAMR)
                nTotalMonte = 100000
                if(usepacketsplitting) then
                   nTotalMonte = nTotalMonte * 100
                end if
             end if
          
       else 
          if(minDepthAMR == maxDepthAMR) then
             if(grid%octreeRoot%twoD) then
                if(hOnly) then
                   nTotalMonte = int(100.d0 * (4.d0**(maxDepthAMR)))
                else
                   nTotalMonte = int(2000.d0 * (4.d0**(maxDepthAMR)))
                end if
             else if(grid%octreeRoot%threeD) then
!                nMonte = 1000.d0
                nTotalMonte = int(100.d0 * (8.d0**(maxDepthAMR)))
!                nMonte = 1.d0 * (8.d0**(maxDepthAMR))
                ! nMonte = 1.d0 * (8.d0**(maxDepthAMR))
                !nMonte = 5242880/2.
             else
                nTotalMonte = int(500.d0 * 2**(maxDepthAMR))
                if(usepacketsplitting) then
                   nTotalMonte = nTotalMonte * 100
                end if
             end if
          else
             call writeInfo("Non uniform grid, setting arbitrary nMonte", TRIVIAL)
             nTotalMonte = 209715200
             write(*,*) "nMonte = ", nTotalMonte
          end if
          
       end if
       
    else
       nTotalMonte = inputnMonte
    endif

    nIter = 0

    converged = .false.

    if (nSource > 1) then
       call randomSource(source, nSource, iSource, photonPacketWeight, lamArray, nLambda, initialize=.true.)
    else
       iSource = 1
    endif

    if (maxiter > 1) then
       if (variableDustSublimation) then
          if (Writeoutput) write(*,*) "Stripping dust away."
          call stripDustAway(grid%octreeRoot, 1.d-20, 1.d30)
          if (Writeoutput) write(*,*) "Done."
       endif
    endif


    do while(.not.converged)


       if(optimizeStack .and. nIter > 0) then
          bufferSize = 0
          allocate(photonPacketStack(stackLimit*nDomainThreads))
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


! Set up the dust fraction
       if (variableDustSublimation) then


          if (nIter == 2) tauMax = 0.1
          if (nIter == 3) tauMax = 1.e0
          if (niter == 4) tauMax = 10.e0
          if (nIter == 5) tauMax = 100.e0
          if (nIter == 6) tauMax = 1000.e0
          if (nIter == 7) tauMax = 1.e30
          nFrac = 0
          totFrac = 0.d0
          if (nIter >= 3) then
             write(message,*) myrankWorldGlobal," Setting dust type max optical depth to ", tauMax
             call writeInfo(message, TRIVIAL)
             call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax, subTemp=1500.d0, minLevel=1.d-10)
          endif
!          call writeVTKfile(grid,"dust.vtk",valueTypeString=(/"dust"/))
          if (maxIter == 1) call  sublimateDust(grid, grid%octreeRoot, totFrac, nFrac, tauMax=1.e30, subTemp=tsub(1))
       end if

       sourceInThickCell = .false.
       maxDiffRadius = 0.d0
       nSmallPackets = 0

       call MPI_BARRIER(localWorldCommunicator,ierr)
       if (.not. cart2d) then
          maxDiffRadius3  = 1.d30
          tauWanted = 1.d0
          do isource = 1, globalnSource
             call tauRadius(grid, globalSourceArray(iSource)%position, VECTOR(-1.d0, 0.d0, 0.d0), tauWanted, &
                  maxDiffRadius1(iSource))
             call MPI_BARRIER(localWorldCommunicator,ierr)
             call MPI_BCAST(maxDiffRadius1(iSource), 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)

             if (amr2d.or.amr3d) then
                call tauRadius(grid, globalSourceArray(iSource)%position,VECTOR(0.d0, 0.d0, -1.d0), tauWanted, &
                     maxDiffRadius2(iSource))
                call MPI_BARRIER(localWorldCommunicator,ierr)
                call MPI_BCAST(maxDiffRadius2(iSource), 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
             endif
             if (amr3D) then
                call tauRadius(grid, globalSourceArray(iSource)%position,VECTOR(0.d0, -1.d0, 0.d0), tauWanted, &
                     maxDiffRadius3(iSource))
                call MPI_BARRIER(localWorldCommunicator,ierr)
                call MPI_BCAST(maxDiffRadius3(iSource), 1, MPI_DOUBLE_PRECISION, 0, localWorldCommunicator, ierr)
             endif
             call MPI_BARRIER(localWorldCommunicator, ierr)
          enddo
          if (splitThisTime) sourceInThickCell = .true.
          
          do isource = 1, globalnSource
             if (amr3d) then
!                if (writeoutput) write(*,*) "r1, r2, r3 ",maxDiffRadius1(isource), maxDiffRadius2(iSource), maxDiffRadius3(iSource)
                maxDiffRadius(isource) = (maxDiffRadius1(isource) +  maxDiffRadius2(iSource) + maxDiffRadius3(iSource))/3.d0
             else if (amr2d) then
!                if (writeoutput) write(*,*) "r1, r2 ",maxDiffRadius1(isource), maxDiffRadius2(iSource)
                maxDiffRadius(isource) = (maxDiffRadius1(isource) +  maxDiffRadius2(iSource))/2.d0
             else
                if (writeoutput) write(*,*) "r1 ",maxDiffRadius1(isource)
               maxDiffRadius(isource) = maxDiffRadius1(isource)
             endif
             if (writeoutput.and.(maxDiffRadius(iSource) /= 0.d0)) &
                 write(*,*) myrankGlobal," Max diffusion radius from tauRadius ",maxDiffRadius(iSource)
          enddo

          
          if (myrankGlobal /= 0) then
             call unsetDiffusion(grid%octreeRoot)
             do i = 1, globalnSource
                call setDiffusionZoneOnRadius(grid%octreeRoot, globalSourceArray(i)%position, maxDiffRadius(i))
                if (maxDiffRadius(i)*1.d10 > globalSourceArray(i)%accretionRadius) sourceInThickCell = .true.
             enddo
          endif
          

          call mpi_allreduce(sourceInThickCell, tempLogical, 1, MPI_LOGICAL, MPI_LOR, localWorldCommunicator, ierr)
          sourceInThickCell = tempLogical
          if (sourceInThickCell.and.usePacketSplitting) then
             nSmallPackets = inputNSmallPackets
             smallPhotonPacketWeight = 1.d0/(dble(nSmallPackets))
          endif
       end if
       if (radPressureTest) nSmallPackets = 0
    if (writeoutput) then
       write(*,*) myrankGlobal, " Setting nSmallPackets to ",nSmallPackets
    endif



      nInf=0
      nMonte = nTotalMonte / max(nSmallPackets,1)
      if (.not.loadBalancingThreadGlobal)  call clearContributions(grid%octreeRoot)

       if(monochromatic) then
          if (source(1)%outsidegrid) then
             if(cart2d) then
                if(grid%octreeroot%oned) then
                   epsoverdeltat = (((nIonizingPhotons*((13.60001)*evtoerg))*&
                        (1.d20*(2.d0*grid%halfsmallestsubcell)**2))/dble(nMonte))
                else
                   epsoverdeltat = (((nIonizingPhotons*((13.60001)*evtoerg))*&
                        (2.d0*grid%octreeRoot%subcellSize*1.d20*2.d0* &
                        grid%halfsmallestsubcell))/dble(nMonte))
                end if
             else
                epsoverdeltat = (((nIonizingPhotons*((inputEV)*evtoerg))*&
                     (2.d0*grid%octreeRoot%subcellSize*1.d10)**2)/dble(nMonte)) 
             end if
          else
             epsoverdeltat = (nIonizingPhotons*((inputEV)*evtoerg))/ &
                  dble(nMonte)
          endif
       else
          epsoverdeltat = lcore/dble(nMonte)
       end if
       
       globalEpsOverDeltaT = epsOverDeltaT




       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if (myrankWorldGlobal == 1) call tune(6, "Setting up load balance")  ! start a stopwatch

       if (firstLoadBalancing.and.(.not.readGrid)) then
          call setLoadBalancingThreadsByCells(grid)
          firstLoadBalancing = .false.
       else
          call  setLoadBalancingThreadsByCrossings(grid)
       endif

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if (myrankWorldGlobal == 1) call tune(6, "Setting up load balance")  ! start a stopwatch

       if (myrankGlobal /= 0) call zeroDistanceGrid(grid%octreeRoot)

       if (myrankGlobal /= 0) call setKappaP(grid%octreeRoot, grid)

       if (myrankWorldGlobal == 1) write(*,*) "Running photoionAMR loop with ",nmonte," photons. Iteration: ",niter, maxIter

       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       if (myrankWorldGlobal == 1) call tune(6, "One photoionization itr")  ! start a stopwatch

       !Thaw - stack optimization
       if(optimizeStack) then
          call wallTime(startTime)
       end if
       nThreadMonte = nMonte / nHydroSetsGlobal
       iMonte_beg = 1
       iMonte_end = nThreadMonte

       Photonmomentum = epsOverDeltaT / cSpeed

       totalPower = 0.d0
       nScatBigPacket = 0
       nScatSmallPacket = 0
       nTotScat = 0
       nPhot = 0
       nSaved = 0
       nBundles = 0
       photonPacketStack%destination = 0
       photonPacketStack%freq = 0.d0
       toSendStack%freq = 0.d0
       toSendStack%destination = 0
       toSendStack%crossedPeriodic = .false.

       countArray = 0.d0

       call zeroNFreq(grid%octreeRoot)

!       call writeVtkFile(grid, "beforeloop.vtk", &
!            valueTypeString=(/"rho          ","dust        ", "HI           " , "temperature  ", &
!            "hydrovelocity","sourceCont   ","pressure     ", &
!            "crossings    "/))!


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

       if (inputSeed /=0 ) then
          call randomNumberGenerator(putISeed = inputSeed) ! TJH added 21/9/15
          call randomNumberGenerator(syncIseed=.true.)
       endif


       if (myRankGlobal == 0) then
             if (myrankWorldGlobal == 0) call tune(6, "All photons sent from rank 0")  ! stop a stopwatch
             mainloop: do iMonte = iMonte_beg, iMonte_end
!                   if ((myHydroSetGlobal == 0).and.&
!                        (mod(iMonte,max(int(1,kind=bigint),(imonte_end-imonte_beg+1)/int(10,kind=bigint),&
!                        kind=bigint) == 0))) &
!                        write(*,*) "imonte ",imonte
!                   if ((myHydroSetGlobal == 1).and.&
!                        (mod(iMonte,(imonte_end-imonte_beg+1)/10) == 0)) write(*,*) "imonte1 ",imonte
!                if (mod(iMonte,(imonte_end-imonte_beg+1)/10) == 0) then
!                   write(*,*) myHydroSetGlobal, " imonte ",imonte
!                endif
                   if (iMonte == nThreadMonte) then
                      lastPhoton = .true.
!                      write(*,*) myrankWorldGlobal, " doing last photon"
                   else
                      lastPhoton = .false.
                   endif
                   
                   call randomSource(source, nSource, iSource, photonPacketWeight, lamArray, nLambda, initialize=.true.)
                thisSource = source(iSource)
                call getPhotonPositionDirection(thisSource, rVec, uHat,rHat,grid)
                if(cart2d .and. .not. source(1)%outsidegrid) then
                   call Pseudo3DUnitVector(uHat, photonPacketWeight,grid%halfsmallestsubcell,&
                        2.d0*grid%octreeRoot%subcellSize)
                end if
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
                spectrumweight = 1.d0

                
                if(monochromatic) then
                   thisFreq = ((inputEV)*evtoerg)/hcgs
                else
                   call getWavelength(thisSource%spectrum, wavelength, spectrumWeight)    
                   photonPacketWeight = photonPacketWeight*spectrumweight
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


! TJH this is where the load balancing thread stuff needs to go for rank 0 thread
! sending photons from stars

                iThread = loadBalancedThreadNumber(thisOctal%mpiThread(subcell))

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
                      photonPacketStack(optCounter)%bigPhotonPacket = .true.
                      photonPacketStack(optCounter)%smallPhotonPacket = .false.
                      photonPacketStack(optCounter)%lastPhoton = .false.
                      exit
                   end if
                end do

                !Keep track of how many packets are destined for each thread
                nSaved(iThread) = nSaved(iThread) + 1

                !Once the bundle for a specific thread has reached a critical size, send it to the thread for propagation
                do optCounter = 1, nDomainThreads
                   if(nSaved(optCounter) /= 0) then
                      if(nSaved(optCounter) == (zerothstackLimit) .or. &
                           (nThreadMonte - nInf) < (zerothstackLimit*nDomainThreads)) then
                         thisPacket = 1
                         do sendCounter = 1, (maxStackLimit*nDomainThreads)
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
                               toSendStack(thisPacket)%bigPhotonPacket = photonPacketStack(sendCounter)%bigPhotonPacket
                               toSendStack(thisPacket)%smallPhotonPacket = photonPacketStack(sendCounter)%smallPhotonPacket
                               toSendStack(thisPacket)%lastPhoton = photonPacketStack(sendCounter)%lastPhoton
                               if (toSendStack(thisPacket)%lastPhoton) then
!                                  write(*,*) myrankWorldGlobal, " sending last packet ",thisPacket, imonte
                               endif
                               thisPacket = thisPacket + 1
                               nInf = nInf + 1
                               
                               !Reset the photon frequency so that this array entry can be overwritten
                               photonPacketStack(sendCounter)%freq = 0.d0
                               
                               !Reset the destination so that it doesnt get picked up in subsequent sweeps!       
                               photonPacketStack(sendCounter)%destination = 0
                            
                            end if
                         end do
                         toSendStack(thisPacket - 1)%lastPhoton = .true.
! stuck here
                         nBundles = nBundles + 1
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
                write(*,*) "Rank 0 sent all initial bundles"
             endif

             if(binPhotons) then
                open(333, file="bins.dat", status="unknown")
                do i=1, nFreq
                   write(333,*) lams(i)*1.e8, countArray(i)/dble(nMonte)
                end do
                close(333)
             end if


             do i = 1, max(nSmallPackets,1)*nBundles
!               write(*,*) myHydroSetGlobal," looping from 1 to ",max(nsmallpackets,1), i
                call MPI_RECV(j, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     finaltag, localWorldCommunicator, status, ierr)
!                write(*,*) myrankWorldglobal, " receiving from ",j
             enddo

             do iThread = 1, nDomainThreads
                toSendStack(1)%destination = 500
                call MPI_SEND(toSendStack, maxStackLimit, MPI_PHOTON_STACK, iThread, tag, localWorldCommunicator,  ierr)
                call MPI_RECV(donePanicking, 1, MPI_LOGICAL, iThread, tag, localWorldCommunicator, status, ierr)  
             end do



             if (myrankWorldGlobal == 0) call tune(6, "All photons sent from rank 0")  ! stop a stopwatch
             if (myrankWorldGlobal==0)  write(*,*) "Telling Ranks to pass stacks ASAP "
                
             photonsStillProcessing = .true.
             
             i = 0

             do while(photonsStillProcessing)                   

!                do iThread = 1, nHydroThreadsGlobal
!                   toSendStack(1)%destination = 600
!                   call MPI_SEND(toSendStack, maxStackLimit, MPI_PHOTON_STACK, iThread, tag, localWorldCommunicator,  ierr)
!                   write(*,*) myrankWorldGlobal, " sent buffer flush signal to  ",ithread
!                end do


                i = i + 1
                toSendStack%freq = 0.d0
                do iThread = 1, nDomainThreads
                   toSendStack(1)%destination = 999
                   !toSendStack(1)%freq = 100.d0
                   call MPI_SEND(toSendStack, maxStackLimit, MPI_PHOTON_STACK, iThread, tag, localWorldCommunicator,  ierr)   
                   call MPI_RECV(nEscapedArray(iThread), 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, status, ierr)
!                   write(*,*) "thread ",ithread," escaped ",nEscapedArray(iThread)
                enddo
                
                nEscaped = SUM(nEscapedArray(1:nDomainThreads))
                nEscapedGlobal = nEscaped
!                if (myRankGlobal == 0) write(*,*) myhydrosetglobal," nEscaped ",nEscaped
                if (nEscaped == nThreadMonte*max(nSmallPackets,1)) then
                   photonsStillProcessing = .false.
                else if(nEscaped > nThreadMonte*max(1,nSmallPackets)) then
                   write(*,*) "nEscaped greater than nMonte, Exiting...",nEscaped
                   do iThread = 1, nDomainThreads
                      write(*,*) "nEscapedArray(iThread) ", iThread, nEscapedArray(iThread)
                   end do
                   photonsStillProcessing = .false.
                end if
!                print *, "nEscaped ", nEscapedGlobal
                if (i > 100000) then
                   write(*,*) "!!! Problem with photon conservation ",nEscaped 
                   exit
                endif
             end do

             if (myrankWorldGlobal == 0) write(*,*) "Finishing iteration..."


             do iThread = 1, nDomainThreads
                !tosendstack(1)%freq = 100.d0
                toSendStack(1)%destination = 888
!                toSendStack(2)%destination = 999
                call MPI_SEND(toSendStack, maxStackLimit, MPI_PHOTON_STACK, iThread, tag, localWorldCommunicator,  ierr)
             enddo



          else

! TJH this is where all the non-zero threads end up initially. Include the load balancing threads
! here


             endLoop = .false.
             nEscaped = 0
             photonPacketStack%freq = 0.d0
             currentStack%freq = 0.d0
             stackSize = 0
             nSaved = 0
             nNotEndLoop = 0
             waitingTime = 0.d0
             call wallTime(globalStartTime)
             sendAllPhotons = .false.
             !needNewPhotonArray = .true.  
             do while(.not.endLoop) 
                nNotEndLoop = nNotEndLoop + 1
                
                crossedMPIboundary = .false.

                !Get a new photon stack
                iSignal = -1
                if(stackSize == 0) then
                   flushBuffer = .true.
                   call walltime(startTime)
                   call MPI_RECV(toSendStack, maxStackLimit, MPI_PHOTON_STACK, MPI_ANY_SOURCE, &
                        tag, localWorldCommunicator, status, ierr)
                   call walltime(finishTime)
                   waitingTime = waitingTime + (finishTime-startTime)

                      currentStack = toSendStack
                      
!                      write(*,*) myrankGlobal, " received stack ",currentStack(1)%destination

                      !Check to see how many photons in stack are not null
                      containsLastPacket = .false.
                      do p = 1, maxStackLimit
                         if(currentStack(p)%freq /= 0.d0) then                            
                            stackSize = stackSize + 1
                            if (currentStack(p)%lastPhoton) containsLastPacket = .true.
                         end if
                      end do
                      receivedStackSize = stackSize
                      nToSend = 0

!                      write(*,*) myrankWorldGlobal, " received a new stack of size ",stacksize,currentStack(1)%destination
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
                   else if (currentStack(1)%destination == 600) then
                      stackSize = 0
                      !Evacuate everything currently in need of sending
                      do optCounter = 1, nDomainThreads
                         if(optCounter /= myRankGlobal .and. nSaved(optCounter) /= 0) then
                            thisPacket = 1
                            
                            toSendStack%freq = 0.d0
                            do sendCounter = 1, (maxStackLimit*nDomainThreads)
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
                                  toSendStack(thisPacket)%bigPhotonPacket = photonPacketStack(sendCounter)%bigPhotonPacket
                                  toSendStack(thisPacket)%smallPhotonPacket = photonPacketStack(sendCounter)%smallPhotonPacket
                                  toSendStack(thisPacket)%lastPhoton = photonPacketStack(sendCounter)%lastPhoton
                                  thisPacket = thisPacket + 1
                                  !nInf = nInf + 1
                                     
                                  !Reset the photon frequency so that this array entry can be overwritten
                                  photonPacketStack(sendCounter)%freq = 0.d0
                                  
                                  !Reset the destination so that it doesnt get picked up in subsequent sweeps!
                                  photonPacketStack(sendCounter)%destination = 0
                                  
                               end if
                            end do
!TJH this was originally mpi_send not mpi_bsend ! 11/9/2012
!                              call MPI_SEND(toSendStack, maxStackLimit, &
!                                   MPI_PHOTON_STACK, OptCounter, tag, localWorldCommunicator,  ierr)
!                               write(*,*) myrankWorldGlobal," flushing  stack to ",optCounter, " size ",thisPacket-1
                               call MPI_BSEND(toSendStack, maxstackLimit, MPI_PHOTON_STACK, OptCounter, tag, &
                                    localWorldCommunicator, ierr)
                              !reset the counter for this thread's bundle recieve
                               nSaved(optCounter) = 0
                               toSendStack%freq = 0.d0
                               toSendStack%destination = 0
                            end if
                         end do
                     
                         goto 777
                         Currentstack%destination = 0

                   else if(currentStack(1)%destination == 500 .and. .not. sendAllPhotons) then
                      sendAllPhotons = .true.
                      stackSize = 0
                      !Evacuate everything currently in need of sending
                      do optCounter = 1, nDomainThreads
                         if(optCounter /= myRankGlobal .and. nSaved(optCounter) /= 0) then
                            if(nSaved(optCounter) == sendStackLimit .or. sendAllPhotons) then
                               thisPacket = 1
                               
                               toSendStack%freq = 0.d0
                               do sendCounter = 1, (maxStackLimit*nDomainThreads)
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
                                    toSendStack(thisPacket)%bigPhotonPacket = photonPacketStack(sendCounter)%bigPhotonPacket
                                    toSendStack(thisPacket)%smallPhotonPacket = photonPacketStack(sendCounter)%smallPhotonPacket
                                    toSendStack(thisPacket)%lastPhoton = photonPacketStack(sendCounter)%lastPhoton
                                    thisPacket = thisPacket + 1
                                    !nInf = nInf + 1
                                    
                                    !Reset the photon frequency so that this array entry can be overwritten
                                    photonPacketStack(sendCounter)%freq = 0.d0
                                    
                                     !Reset the destination so that it doesnt get picked up in subsequent sweeps!
                                    photonPacketStack(sendCounter)%destination = 0
                                    
                                 end if
                              end do
!TJH this was originally mpi_send not mpi_bsend ! 11/9/2012
!                              call MPI_SEND(toSendStack, maxStackLimit, &
!                                   MPI_PHOTON_STACK, OptCounter, tag, localWorldCommunicator,  ierr)
                              call MPI_BSEND(toSendStack, maxstackLimit, MPI_PHOTON_STACK, OptCounter, tag, &
                                   localWorldCommunicator, &
                                   ierr)
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
               !$OMP PRIVATE(p, rVec, uHat, thisFreq, tPhoton, photonPacketWeight, sourcePhoton, r1) &
               !$OMP PRIVATE(escaped, nScat, optCounter, octVec, ierr, thisLam, kappaabsdb) &
               !$OMP PRIVATE(doingsmallpackets,startnewsmallpacket,ismallphotonpacket,bigphotonpacket) &
               !$OMP PRIVATE(smallphotonpacket,smallpacketorigin,smallpacketfreq,smallphotonpacketweight,kappap) &
               !$OMP PRIVATE(kappascadb, albedo, r, kappaabsdust, thisOctal, subcell) &
               !$OMP PRIVATE(crossedMPIboundary, newThread, thisPacket, kappaabsgas, escat, tempcell, lastPhoton) &
               !$OMP PRIVATE(finished, voidThread, crossedPeriodic, nperiodic,  myrankworldglobal) &
               !$OMP PRIVATE(bigPhotonPacketWeight, iLam, flushbuffer,ntosend) &
               !$OMP PRIVATE(uHatBefore, vec_tmp, unew, uhatafter, uHatDash, rHat, zHat, movedCells) & 
               !$OMP SHARED(photonPacketStack, myRankGlobal, currentStack, escapeCheck, cart2d) &
               !$OMP SHARED(noDiffuseField, grid, epsoverdeltat, iSignal, MPI_PHOTON_STACK) &
               !$OMP SHARED(nlambda, lamarray, tlimit, nHydroThreadsGlobal, sendAllPhotons,toSendStack) &
               !$OMP SHARED(nTotScat, nScatbigPacket, nScatSmallPacket, gammaTableArray, freq, nsmallpackets) &
               !$OMP SHARED(dfreq, endLoop, nIter, spectrum, sendStackLimit) &
               !$OMP SHARED(nSaved, maxStackLimit, source, uv_vector, containslastpacket) &
               !$OMP SHARED(stackSize, nFreq, radPressureTest, augerArray, firstwarning) &
               !$OMP SHARED(nPhot, nEscaped, stackLimit, localWorldCommunicator, nhydrosetsglobal, nToNextEventPhoto, nNotEscaped) &
               !$OMP SHARED(miephase, nmumie, ndusttype, photonmomentum, loadBalancing, nDomainThreads)
               
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
                        bigPhotonPacket = currentStack(p)%bigPhotonPacket
                        smallPhotonPacket = currentStack(p)%smallPhotonPacket
                        lastPhoton = currentStack(p)%lastPhoton
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

!               write(*,*) myrankGlobal, " received photon packet at ", rVec, &
!                     octalOnThread(thisOctal,subcell,myrankGlobal),thisOctal%rho(subcell), &
!                     thisOctal%nDepth, thisOctal%subcellsize,thisOctal%ionfrac(subcell,:)

                doingSmallPackets = .false.
                startNewSmallPacket = .false.
                nToNextEventPhoto = 0

                if (.not.endLoop) then
                   nScat = 0
                   escaped = .false.
                   nNotEscaped = 0
                   do while(.not.escaped)

                      if ((doingSmallPackets).and.(startNewSmallPacket) .and. .not. cart2d) then
                         startNewSmallPacket = .false.
                         iSmallPhotonPacket = iSmallPhotonPacket + 1
                         if (iSmallPhotonPacket > nSmallPackets) then
                            doingSmallPackets = .false.
!                            write(*,*) myrankWorldGlobal," Finished doing small photon packets"
                            goto 555
                         else
                            rVec = smallPacketOrigin
                            thisFreq = smallPacketFreq
                            smallPhotonPacketWeight=smallPhotonPacketWeight
                            photonPacketWeight = smallPhotonPacketWeight * bigPhotonPacketWeight
                            Uhat = randomUnitVector()
!                            if(cart2d) then
!                            if(cart2d .and. .not. source(1)%outsidegrid) then
!                               call Pseudo3DUnitVector(uHat, photonPacketWeight,grid%halfsmallestsubcell,&
!                                    2.d0*grid%octreeRoot%subcellSize)
!                            end if
!                            if(grid%twoD) then
!                               thisOctal%
!                            end if
                            smallPhotonPacket = .true.
                            bigPhotonPacket = .false.
                            call findSubcellTD(rVec, grid%octreeRoot,thisOctal, subcell)
                            kappaP = thisOctal%kappaP(subcell)
!                            call returnKappa(grid, thisOctal, subcell, kappap=kappap)
!                            write(*,*) myrankWorldGlobal, " creating small photon packet ",ismallPhotonPacket, &
!                                 thisOctal%subcellSize*kappap*1.d10, (cspeed/thisFreq)*1.d8
!                            write(*,*) "escaped ",escaped
                         endif
                      endif
                         

                      nToNextEventPhoto = nToNextEventPhoto + 1
                      if (nToNextEventPhoto > 80000) then
!                         write(*,*) "too next event photo called ",nToNextEventPhoto
 !                        write(*,*) "rank ",myrankWorldGlobal
 !                        write(*,*) "rVec ",rVec
 !                        write(*,*) "uHat ",uhat
 !                        write(*,*) "thisOctal%xmin ",thisOctal%xmin
 !                        write(*,*) "thisOctal%xmax ",thisOctal%xmax
                      endif
                      
!                      call findsubcellTd(rVec, grid%octreeRoot, beforeOctal, beforeSubcell)

                      call toNextEventPhoto(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray, &
                           photonPacketWeight, epsOverDeltaT, nfreq, freq, dFreq, tPhoton, tLimit, &
                           crossedMPIboundary, newThread, sourcePhoton, crossedPeriodic, miePhase, nMuMie, movedCells)

                      if (nToNextEventPhoto == 1) movedCells = .true. ! a new photon has effectively moved cells
                      

                      if (loadBalancing.and.crossedMPIboundary) newThread = loadBalancedThreadNumber(newThread)

!                      if (escaped.and.bigPHotonPacket) then
!                         write(*,*) myrankGlobal, " big photon packet escaped ",rVec
!                      endif

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
                               photonPacketStack(optCounter)%destination = newthread !loadBalancedThreadNumber(newThread)
                               photonPacketStack(optCounter)%sourcePhoton = sourcePhoton
                               photonPacketStack(optCounter)%crossedPeriodic = crossedPeriodic
                               photonPacketStack(optCounter)%bigPhotonPacket = bigPhotonPacket
                               photonPacketStack(optCounter)%smallPhotonPacket = smallPhotonPacket
                               photonPacketStack(optCounter)%lastPhoton = lastPhoton
                               if(crossedPeriodic) then
                                  nPeriodic = nPeriodic + 1
                               end if
                               exit
                            end if
                         end do

                         !Keep track of how many packets are destined for each thread
                         nSaved(newThread) = nSaved(newThread) + 1
                         nToSend = nToSend + 1
                         flushBuffer = containsLastPacket
                         if (flushBuffer) then
!                            write(*,*) myrankWorldGlobal, " flushing buffer with ",nToSend
                         endif

                         !Once the bundle for a specific thread has reached a critical size, send it to the thread for propagation
                         do optCounter = 1, nDomainThreads
                            if(optCounter /= myRankGlobal .and. nSaved(optCounter) /= 0) then
                               if(nSaved(optCounter) == (sendStackLimit) .or. sendAllPhotons.or.flushBuffer) then
                                  thisPacket = 1
                                  toSendStack%freq = 0.d0
                                  do sendCounter = 1, (maxStackLimit*nDomainThreads)
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
                                        toSendStack(thisPacket)%bigPhotonPacket = photonPacketStack(sendCounter)%bigPhotonPacket
                                        toSendStack(thisPacket)%smallPhotonPacket = photonPacketStack(sendCounter)%smallPhotonPacket
                                        toSendStack(thisPacket)%lastPhoton = photonPacketStack(sendCounter)%lastPhoton
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
                                          localWorldCommunicator, ierr)

                                     nSaved(optCounter) = 0
                                     toSendStack%freq = 0.d0
                                     toSendStack%destination = 0
                                  end if
                                  
                               end if
                            end if
                         end do
                         !$OMP END CRITICAL (crossed_mpi_boundary)
                         finished = .true.
!                         if (doingsmallpackets) &
!                              write(*,*) myrankGlobal, " set finished true ", ismallPhotonPacket,doingsmallPackets
                      endif

!                      if ((.not.crossedMPIboundary).and.finished.and.doingSmallPackets &
!                           .and.escaped) write(*,*) "ismall ",ismallphotonpacket

                      

                      if(.not. finished) then
                         if (noDiffuseField) escaped = .true.
                         
                         UhatBefore = uHat

                         if (escaped) then
                            !$OMP CRITICAL (update_escaped)
                            if (smallPhotonPacket) nEscaped = nEscaped + 1
                            if (bigPhotonPacket) nEscaped = nEscaped + max(nSmallPackets,1)
                            if (lastPhoton.and.smallPhotonPacket) then
!                               write(*,*) myrankWorldGlobal, " last small photon escaped ",rVec,inOctal(grid%octreeRoot, rVec)
                               call MPI_BSEND(myrankGlobal, 1, MPI_INTEGER, 0, finaltag, localWorldCommunicator,  ierr)
                            endif
                            if (lastPhoton.and.bigPhotonPacket) then
!                               write(*,*) myrankWorldGlobal, " last big photon escaped ",rVec,inOctal(grid%octreeRoot,rvec)
                               do i = 1, max(1,nSmallPackets)
                                  call MPI_BSEND(myrankGlobal, 1, MPI_INTEGER, 0, finaltag, localWorldCommunicator,  ierr)
                               enddo
                            endif
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

!                               uHat = randomUnitVector() ! isotropic scattering

                               vec_tmp = uHat
                               uNew = newDirectionMie(grid, thisOctal, subcell, vec_tmp, real(thisLam), lamArray, nLambda, &
                                    miePhase, nDustType, nMuMie)
!                               write(*,*) "scattering angle ",radtodeg*acos(uNew.dot.uHat)
                               uHat = uNew


!                               if(cart2d) then
                               if(cart2d) then !.and. .not. source(1)%outsidegrid) then
                                  photonpacketweight = 1.d0
                                  call Pseudo3DUnitVector(uHat, photonPacketWeight,grid%halfsmallestsubcell,&
                                       2.d0*grid%octreeRoot%subcellSize)
                               end if
!                               if (myrankglobal == 1) write(*,*) "new uhat scattering ",uhat, ntotscat
                            else
                               
                               if (thisOctal%nFreq(subcell) == 0) then ! only need to recalculate emissivity spectrum
                                  thisOctal%nFreq(subcell) = nFreq
                                  if (.not.associated(thisOctal%spectrum)) then
                                     allocate(thisOctal%spectrum(1:thisOctal%maxChildren, 1:nfreq))
                                  endif
                               spectrum = 1.d-30
                               
                               call returnKappa(grid, thisOctal, subcell, ilambda=ilam, &
                                    kappaAbsDust=kappaAbsDust, kappaAbsGas=kappaAbsGas, &
                                    kappaSca=kappaScadb, kappaAbs=kappaAbsdb, kappaScaGas=escat)
                                               
                                              
                               if ((thisFreq*hcgs*ergtoev) > 13.6) then ! ionizing photon
!                                  if (thisOctal%temperature(subcell) == 0.0) then
!                                     write(*,*) "temperature of cell is zero"
!                                     write(*,*) "cell is ghost " ,thisOctal%ghostCell(subcell)
!                                     write(*,*) "edge? ", thisOctal%edgeCell(subcell)
!                                  endif
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
                                  call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, &
                                       subcell, grid, nlambda, lamArray)
                               endif


!                               if (firsttime.and.(myrankWorldglobal==1)) then
!                                  firsttime = .false.
!                                  open(67,file="pdf.dat",status="unknown",form="formatted")
!                                  write(67,*) "% ",thisOctal%temperature(subcell)
!                                  do i = 1, nfreq
!                                     write(67,*) (cspeed/freq(i))*1.e8, spectrum(i), &
!                                          bnu(freq(i),dble(thisOctal%temperature(subcell)))
!                                  enddo
!                                  close(67)
!                                  open(68,file="freq.dat",status="unknown",form="formatted")
!                                  do i = 1, nfreq
!                                     write(68,*) i, freq(i)
!                                  enddo
!                                  close(68)
!                               endif

                               thisOctal%spectrum(subcell,1:nfreq) = spectrum(1:nFreq)
                            endif

                               thisFreq =  getPhotonFreq(nfreq, freq, thisOctal%spectrum(subcell,1:nFreq))

                               uHat = randomUnitVector() ! isotropic emission
!                               if(cart2d) then
                               if(cart2d) then !.and. .not. source(1)%outsidegrid) then
                                  photonpacketweight = 1.d0
                                  call Pseudo3DUnitVector(uHat, photonPacketWeight,grid%halfsmallestsubcell,&
                                       2.d0*grid%octreeRoot%subcellSize)
                               end if
                               if (radpressuretest) thisFreq = freq(1)

                            endif

                            UhatAfter = uHat
                            if (escaped) UhatAfter = VECTOR(0.d0, 0.d0, 0.d0)
                                                       

                            uHatDash = uHatBefore-uHatAfter
                            if (thisOctal%twoD .and. .not. cart2d) then
                               rHat = VECTOR(rVec%x, rVec%y, 0.d0)
                               call normalize(rHat)
                               zHat = VECTOR(0.d0, 0.d0, 1.d0)
                               uHatDash = VECTOR(rHat.dot.(uHatBefore-uHatAfter), 0.d0, zHat.dot.(uHatBefore-uHatAfter))
                            endif


                            if (thisOctal%oneD) then
                               rHat = VECTOR(rVec%x, rVec%y, rVec%z)
                               call normalize(rHat)
                               uHatDash = VECTOR(rHat.dot.(uHatBefore - uHatAfter), 0.d0, 0.d0)
                            endif



                            thisOctal%radiationMomentum(subcell) = thisOctal%radiationMomentum(subcell) + uHatDash * &
                                 photonMomentum * photonPacketWeight

                               nScat = nScat + 1

                               if (bigPhotonPacket) nScatBigPacket = nScatBigPacket + 1
                               if (smallPhotonPacket) nScatSmallPacket = nScatSmallPacket + 1
                               nTotScat = nTotScat + 1
                               thisOctal%chiLine(subcell) = thisOctal%chiline(subcell) + 1.d0
                               

! ok here we can split the current photon into multiple smaller photon packets if necessary

                               if ((bigPhotonPacket).and.(.not.thisOctal%diffusionApprox(subcell)).and. &
                                    (nSmallPackets > 0)) then
!                                  write(*,*) myrankWorldGlobal, " splitting bigphoton into small photons"
                                  bigPhotonPacket = .false.
                                  smallPhotonPacket = .true.
                                  iSmallPhotonPacket = 0
                                  doingSmallPackets = .true.
                                  bigphotonPacketWeight = photonPacketWeight
                                  smallPacketOrigin = rVec
                                  smallPacketFreq = thisFreq
                                  startNewSmallPacket = .true.
                                  thisLam = real(cSpeed / thisFreq) * 1.e8
                                  call locate(lamArray, nLambda, real(thisLam), iLam)
                                  
                                  call returnKappa(grid, thisOctal, subcell, ilambda=ilam, &
                                 kappaSca=kappaScadb, kappaAbs=kappaAbsdb)

!                                  write(*,*) "small packet wavelength ", thisLam
!                                  write(*,*) "optical depth of cell ",(kappaAbsDb+kappaScadb)*thisOctal%subcellSize
                               endif
                               
                               
                            endif
                         if (nScat > 1000000) then
                            write(*,*) "Nscat exceeded 1000000, forcing escape"
                            write(*,*) 1.e8*cspeed/thisFreq
                            write(*,*) albedo, kappaScaDb, kappaAbsdb,escat
                            
                            thisLam = real(cSpeed / thisFreq) * 1.e8
                            call locate(lamArray, nLambda, real(thisLam), iLam)
                            
                            call returnKappa(grid, thisOctal, subcell, ilambda=ilam, &
                                 Kappaabsdust=kappaAbsDust, kappaAbsGas=kappaAbsGas, &
                                 kappaSca=kappaScadb, kappaAbs=kappaAbsdb, kappaScaGas=escat)
!                            write(*,*) "thislam",thislam,ilam, lamArray(ilam)
!!                            write(*,*) lamArray(1:nLambda)
 !                           write(*,*) "kappaAbsDust",kappaAbsDust
 !                           write(*,*) "kappaAbsGas",kappaAbsGas
 !                           write(*,*) "kappaSca",kappaScadb
 !                           write(*,*) "kappaAbs",kappaAbsdb
 !                           write(*,*) "onekappaabs",grid%oneKappaAbs(1,ilam)
 !                           write(*,*) "onekappasca",grid%oneKappasca(1,ilam)
                  
                            escaped = .true.
                            
                         Endif
                      end if

                      if (escaped.and.doingSmallPackets) then
                         escaped = .false.
                         finished = .false.
                         startNewSmallPacket = .true.
!                         write(*,*) myrankWorldGlobal, " small photon packet escaped, getting new one"
                      endif
                      nNotEscaped = nNotEscaped + 1
!                      if (nNotEscaped > 10000) then
!                         write(*,*) "been around the not escaped loop too many times"
!                         write(*,*) "finished ",finished
!                         write(*,*) "rVec ",rVec
!                         write(*,*) "nScat ",nScat
!                         escaped =  .true.
!                      endif
                   enddo
!                   write(*,*) "photon escaped after ",nscat," scatterings"
                   nPhot = nPhot + 1
                end if
555             continue
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
          call MPI_RECV(stackLimit, 1, MPI_INTEGER, 1, tag, localWorldCommunicator, status, ierr)
       end if

       call MPI_BARRIER(localWorldCommunicator, ierr)

!       epsOverDeltaT = (lCore) / dble(nMonte)
       if(monochromatic) then
          if (source(1)%outsidegrid) then

             if(cart2d) then
                if(grid%octreeroot%oned) then
                   epsoverdeltat = (((nIonizingPhotons*((13.60001)*evtoerg))*&
                        (1.d20*2.d0*grid%halfsmallestsubcell)**2)/dble(nMonte))
                else
                   epsoverdeltat = (((nIonizingPhotons*((13.60001)*evtoerg))*&
                        (2.d0*grid%octreeRoot%subcellSize*1.d20*2.d0* &
                        grid%halfsmallestsubcell))/dble(nMonte))                   
                end if
             else
                epsoverdeltat = (((nIonizingPhotons*((inputEV)*evtoerg))*&
                     (2.d0*grid%octreeRoot%subcellSize*1.d10)**2)/dble(nMonte))
             end if
          else
             epsoverdeltat = (nIonizingPhotons*((inputEV)*evtoerg))/ &
                  dble(nMonte)
          endif
       else
          epsoverdeltat = lcore/dble(nMonte)
!          print *, "epsAlarm"
       end if

 

       efficiencyArray = 0.d0
       if (myrankGlobal /= 0) then
          call wallTime(globalTime)
          globalTime = globalTime - globalStartTime

          efficiencyArray(myrankGlobal) = 100.d0*(1.d0 - dble(waitingTime/GlobalTime))
          allocate(temp1d(1:nThreadsGlobal-1))
          call MPI_ALLREDUCE(efficiencyArray, temp1d, nThreadsGlobal-1, &
               MPI_DOUBLE_PRECISION, MPI_SUM, allDomainsCommunicator, ierr)
          efficiencyArray = temp1d
          deallocate(temp1d)
          call sort(nThreadsGlobal-1, efficiencyArray)
          if (myrankGlobal == 1) then
             do i = 1, nThreadsGlobal-1
                write(*,*) i, " efficiency ", efficiencyArray(i), " %"
             enddo
          endif
          if (writeoutput) write(*,*) "Median efficiency is ",efficiencyArray((nThreadsGlobal-1)/2), " %"
          medianEfficiency = efficiencyArray((nThreadsGlobal-1)/2)
          oldStackLimit = stackLimit

!          if (medianEfficiency > oldMedianEfficiency) then
!             stackLimit = stackLimit + 1
!          else
!             stackLimit = max(1, sackLimit-1)
!          endif
          

          


          call MPI_ALLREDUCE(nTotScat, i, 1, MPI_INTEGER, MPI_SUM, allDomainsCommunicator, ierr)
          nTotScat = i
          if (writeoutput) write(*,*) "Iteration had ",&
               real(nTotScat)/real(nThreadMonte*max(nSmallPackets,1)), " scatters per photon packet"
          call MPI_ALLREDUCE(nScatSmallPacket, i, 1, MPI_INTEGER, MPI_SUM, allDomainsCommunicator, ierr)
          nScatSmallPacket = i
          if (nSmallPackets > 0) then
             if (writeoutput) write(*,*) "Iteration had ",&
                  real(nScatSmallPacket)/real(nThreadMonte*nSmallPackets), " scatters per small photon packet"
          endif

          
          if (nSmallPackets > 0) then
             call MPI_ALLREDUCE(nScatBigPacket, i, 1, MPI_INTEGER, MPI_SUM, allDomainsCommunicator, ierr)
             nScatBigPacket = i
             if (writeoutput) write(*,*) "Iteration had ",&
                  real(nScatBigPacket)/real(nThreadMonte), " scatters per big photon packet"
             if (real(nScatBigPacket)/real(nThreadMonte) > 20.) then
                splitThisTime = .true.
             else
                splitThisTime = .false.
             endif
          endif

          if (nSmallPackets == 0) then
             if (real(nTotScat)/real(nThreadMonte) > 20.) then
                splitThisTime = .true.
             else
                splitThisTime = .false.
             endif
          endif


       endif


    np = 1
    firstTime = .true.



! now we need to reduce the estimators across all the hydro sets.
    call mpi_barrier(MPI_COMM_WORLD, ierr)


    if (myrankWorldGlobal == 1) call tune(6, "Update MPI grid")  ! start a stopwatch


    if (loadBalancing) then
       do iThread = 1, nHydroThreadsGlobal
          if (nLoadBalanceList(iThread) > 1) then
             if (ANY(loadBalanceList(iThread,1:nLoadBalanceList(iThread)) == myRankGlobal)) then
                call updateGridMPIPhoto(grid, loadBalanceCommunicator(iThread))
             endif
          endif
       enddo
    else
       if (nHydroSetsGlobal > 1) then
          do iThread = 1, nHydroThreadsGlobal
             if (myRankGlobal == iThread) then
                call updateGridMPIphoto(grid, amrParallelCommunicator(iThread))
             endif
             call mpi_barrier(MPI_COMM_WORLD, ierr)
          enddo
       endif
    endif




    call mpi_barrier(MPI_COMM_WORLD, ierr)
    if (myrankWorldGlobal == 1) call tune(6, "Update MPI grid")  ! stop a stopwatch

    call  identifyUndersampled(grid%octreeRoot)



    if (myRankWorldGlobal == 0) write(*,*) "Ninf ",ninf
    if (myrankWorldGlobal == 1) call tune(6, "Temperature/ion corrections")


    if (writeoutput) &
         write(*,*) "Calculating ionization and thermal equilibria"

!    if (myrankWorldGlobal == 1) then
!       i = 0
!       j = 0 
!       call testforZero(grid%octreeRoot, i, j)
!       write(*,*) "photoion(subcell,1) is zero for ", 100.d0*real(i)/real(j), "% of subcells"
!    endif


!    call writeVtkFile(grid, "beforethermalcalc.vtk", &
!         valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  ", &
!         "hydrovelocity","sourceCont   ","pressure     ", &
!         "crossings    "/))!


    if (loadbalancing) call  setLoadBalancingThreadsByCells(grid)

    if (allocated(octalArray)) deallocate(octalArray)
    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif


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
          maxChild = grid%octreeRoot%maxChildren
          allocate(tempCell(1:nOctal,1:maxChild))
          tempCell = 0.
          allocate(tempIon(1:nOctal,1:maxChild,1:grid%nIon))
          tempIon = 0.
       endif

       if (loadBalancing) then
          np = nLoadBalanceList(copyOfThread)
          
          call MPI_COMM_RANK(loadBalanceCommunicator(copyOfThread), localRank, ierr)

          n_rmdr = MOD(nOctal, np)
          m = nOctal/np
            
          if (localRank .lt. n_rmdr ) then
             ioctal_beg = (m+1)*localRank + 1
             ioctal_end = ioctal_beg + m
          else
             ioctal_beg = m*localRank + 1 + n_rmdr
             ioctal_end = ioctal_beg + m - 1
          end if
  
          maxChild = grid%octreeRoot%maxChildren
          allocate(tempCell(1:nOctal,1:maxChild))
          tempCell = 0.
          allocate(tempIon(1:nOctal,1:maxChild,1:grid%nIon))
          tempIon = 0.
       endif

       !$OMP PARALLEL DEFAULT(NONE) &
       !$OMP PRIVATE(iOctal, thisOctal, subcell, v, kappap, i, endTime,startTime) &
       !$OMP PRIVATE(dustHeating, tempcell, oldf, oldt, iter, fract, fracf, converged, tempion) &
       !$OMP SHARED(iOctal_beg, iOctal_end, dustOnly, octalArray, grid, epsOverDeltaT, uv_vector) &
       !$OMP SHARED(timedep, quickThermal, deltaTime, tminGlobal, myrankGlobal, nhydrosetsglobal) &
       !$OMP SHARED(augerArray, firstwarning, radpressuretest, loadbalancing, copyOfThread)

       call walltime(startTime)

       !$OMP DO SCHEDULE(DYNAMIC,2)
       do iOctal =  iOctal_beg, iOctal_end
          
          thisOctal => octalArray(iOctal)%content

          if (.not.radpressuretest) then
          
          if (dustOnly) then
             do subcell = 1, thisOctal%maxChildren
                if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
                if (.not.thisOctal%hasChild(subcell)) then
                   v = cellVolume(thisOctal, subcell)
                   call returnKappa(grid, thisOctal, subcell, kappap=kappap)
                   dustHeating = (epsOverDeltaT / (v * 1.d30))*thisOctal%distanceGrid(subcell) ! equation 14 of Lucy 1999
                   kappap = max(1.d-30,kappap)
                   thisOctal%temperature(subcell) = max(tMinGlobal,real((pi/stefanBoltz) * dustHeating / (fourPi * kappaP))**0.25e0)
                endif
             enddo

          else             
!             if(.not. uv_vector) then

             if(timeDep) then
                do i = 1, nTimes
                   call calculateIonizationBalanceTimeDep(grid,thisOctal, epsOverDeltaT, deltaTime/dble(nTimes))
                   if (quickThermal) then
                      call quickThermalCalc(thisOctal)
                   else
                      call calculateThermalBalanceTimeDep(grid, thisOctal, epsOverDeltaT, deltaTime/dble(nTimes))
                   end if
                enddo
                
             else
                oldT(1:thisOctal%maxChildren) = thisOctal%temperature(1:thisOctal%maxChildren)
                if (.not.associated(thisOctal%ionFrac)) then
                   write(*,*) "not associated ",myrankGlobal,copyOfThread
                endif
                oldF(1:thisOctal%maxChildren) = thisOctal%ionFrac(1:thisOctal%maxChildren,1)
                converged = .false.
                iter = 0
                do while(.not.converged)
                   iter = iter + 1
                   call calculateIonizationBalance(grid,thisOctal, epsOverDeltaT, augerArray)
                   if (quickThermal) then
                      call quickThermalCalc(thisOctal)
                   else
                      call calculateThermalBalanceNew(grid, thisOctal, epsOverDeltaT)
                   endif
                   fracT(1:thisOctal%maxChildren) = abs(thisOctal%temperature(1:thisOctal%maxChildren) - &
                        oldT(1:thisOctal%maxChildren))/oldT(1:thisOctal%maxChildren)
                   fracF(1:thisOctal%maxChildren) = abs(thisOctal%ionFrac(1:thisOctal%maxChildren,1) &
                        - oldF(1:thisOctal%maxChildren))/oldF(1:thisOctal%maxChildren)
                   oldT(1:thisOctal%maxChildren) = thisOctal%temperature(1:thisOctal%maxChildren)
                   oldF(1:thisOctal%maxChildren) = thisOctal%ionFrac(1:thisOctal%maxChildren,1)
                   converged = (ALL(fracT(1:thisOctal%maxChildren) < 0.01d0)).and. (ALL(FracF(1:thisOctal%maxChildren) < 0.01d0))
                   !                      if (inSubcell(thisOctal, 1, VECTOR(0.d0, 0.d0, 0.d0)).and.(octalOnThread(thisOctal,1,myRankWorldGlobal))) then
                   !                         write(*,*) i, "temp, frac ",thisOctal%temperature(1),thisOctal%ionFrac(1,1)
                   !                      endif
                   if (iter > 100) then
                      if (firstWarning) write(*,*) "temperature/ionization fraction not converged after 100 iterations "!,thisOctal%ghostCell(subcell)
                      firstWarning = .false.
                      converged = .true.
                   endif
                   !                   write(*,*) "ion/thermal converged after ",iter
                enddo
             end if

             
          endif
       endif
       !end if
          if (nHydroSetsGlobal > 1) tempCell(iOctal,1:thisOctal%maxChildren) = thisOctal%temperature(1:thisOctal%maxChildren)
          if (nHydroSetsGlobal > 1) tempIon(iOctal,1:thisOctal%maxChildren,:) = real(thisOctal%ionFrac(1:thisOctal%maxChildren,:))

          if (loadBalancing) tempCell(iOctal,1:thisOctal%maxChildren) = thisOctal%temperature(1:thisOctal%maxChildren)
          if (loadBalancing) tempIon(iOctal,1:thisOctal%maxChildren,:) = real(thisOctal%ionFrac(1:thisOctal%maxChildren,:))


       enddo
       !$OMP END DO
       !$OMP END PARALLEL

       if (nHydroSetsGlobal > 1) then
          allocate(temp1(1:(nOctal*maxChild)), temp2(1:nOctal*maxChild))
          temp1 = RESHAPE(tempCell, (/SIZE(temp1)/))
          call MPI_ALLREDUCE(temp1, temp2, SIZE(temp1), MPI_REAL, MPI_SUM, amrParallelCommunicator(myRankGlobal),ierr)
          tempCell = RESHAPE(temp2,SHAPE(tempCell))
          do iOctal = 1, nOctal
             octalArray(iOctal)%content%temperature(1:maxChild) = tempCell(iOctal,1:maxChild)
          enddo
          deallocate(temp1, temp2, tempCell)
          allocate(temp1(1:(nOctal*maxChild*grid%nIon)), temp2(1:nOctal*maxChild*grid%nion))
          temp1 = RESHAPE(tempIon, (/SIZE(temp1)/))
          call MPI_ALLREDUCE(temp1, temp2, SIZE(temp1), MPI_REAL, MPI_SUM, amrParallelCommunicator(myRankGlobal),ierr)
          tempIon = RESHAPE(temp2,SHAPE(tempIon))
          do iOctal = 1, nOctal
             octalArray(iOctal)%content%ionFrac(1:maxChild,:) = tempIon(iOctal,1:maxChild,:)
          enddo
          deallocate(temp1, temp2, tempIon)
          call resetNe(grid%octreeRoot)
       endif

       if (loadBalancing) then
          allocate(temp1(1:(nOctal*maxChild)), temp2(1:nOctal*maxChild))
          temp1 = RESHAPE(tempCell, (/SIZE(temp1)/))
          call MPI_ALLREDUCE(temp1, temp2, SIZE(temp1), MPI_REAL, MPI_SUM, loadBalanceCommunicator(copyOfThread),ierr)
          tempCell = RESHAPE(temp2,SHAPE(tempCell))
          do iOctal = 1, nOctal
             octalArray(iOctal)%content%temperature(1:maxChild) = tempCell(iOctal,1:maxChild)
          enddo
          deallocate(temp1, temp2, tempCell)
          allocate(temp1(1:(nOctal*maxChild*grid%nIon)), temp2(1:nOctal*maxChild*grid%nion))
          temp1 = RESHAPE(tempIon, (/SIZE(temp1)/))
          call MPI_ALLREDUCE(temp1, temp2, SIZE(temp1), MPI_REAL, MPI_SUM, loadBalanceCommunicator(copyOfThread),ierr)
          tempIon = RESHAPE(temp2,SHAPE(tempIon))
          do iOctal = 1, nOctal
             octalArray(iOctal)%content%ionFrac(1:maxChild,:) = tempIon(iOctal,1:maxChild,:)
          enddo
          deallocate(temp1, temp2, tempIon)
          call resetNe(grid%octreeRoot)
       endif
       call walltime(endTime)
!       write(*,*) myrankGlobal, " did thermal balance in ",endTime - startTime
    endif

    !deallocate(octalArray)
    call adjustChiline(grid%octreeRoot)
    
    call calculateKappaTimesFlux(grid%octreeRoot, epsOverDeltaT)
    if(uv_vector) then
       call calculateUVfluxVec(grid%octreeroot, epsoverdeltaT)
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    if (writeoutput) &
         write(*,*) "Finished calculating ionization and thermal equilibria"
    if (myrankWorldGlobal == 1) call tune(6, "Temperature/ion corrections")


    if(grid%geometry == "lexington" .or. grid%geometry == "lexpdr") then
      if (.not.loadBalancingThreadGlobal) call dumpLexingtonMPI(grid, epsoverdeltat, niter)
   end if

    if(grid%geometry == "lexington" .and. 0 == 1) then
      call dumpLexingtonMPI(grid, epsoverdeltat, niter)

      if(myrankGlobal /= 0) then
         call getHbetaLuminosity(grid%octreeRoot, grid, luminosity1)         
!         print *, "rank ", myRankWorldGlobal, "sent ", luminosity1/1.e37, "to 0"
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "N II", 1.22d6, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "N II", 6584.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "N II", 6548.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "N III", 5.73d5, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O I", 6300.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O I", 6363.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O II", 7320.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O II", 7330.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O II", 3726.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O II", 3729.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O III", 5007.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O III", 4959.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O III", 4363.d0, luminosity3)
         call MPI_SEND(luminosity3, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O III", 518145.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "O III", 883562.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "Ne II", 1.28d5, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "Ne III", 1.56d5, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "Ne III", 3869.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "Ne III", 3968.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "S II", 6716.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "S II", 6731.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "S II", 4068.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "S II", 4076.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "S III", 1.87d5, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "S III", 9532.d0, luminosity1)
         call MPI_SEND(luminosity1, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
         call getForbiddenLineLuminosity(grid, "S III", 9069.d0, luminosity2)
         call MPI_SEND(luminosity2, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator,  ierr)
      else if (0 == 1) then
         zTemp = 0.d0
         do recCounter = 1, 27
            do threadCounter = 1, nHydroThreadsGlobal
               call MPI_RECV(thisLum, 1, MPI_DOUBLE_PRECISION, threadCounter, tag, localWorldCommunicator, status, ierr)
               zTemp = zTemp + thisLum
!               print *, "0 got ", thisLum/1.e37, "from ", threadCounter
!               print *, "total so far is ", zTemp/1.e37
            end do

            select case(recCounter)
               case(1)
                  luminosity1 = zTemp
                  write(*,'(a20,2f12.4)') "H beta :",luminosity1/1.e37,luminosity1/2.05e37         
                  fac = luminosity1
               case(2)
                  luminosity1 = zTemp
                  write(*,'(a20,2f12.4)') "N II (122 um):",(luminosity1)/fac,(luminosity1)/(0.034*2.05e37)
               case(3)
                  luminosity1 = zTemp
               case(4)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "N II (6584+6548):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(0.730*2.05e37)
               case(5)
                  luminosity1 = zTemp
                  write(*,'(a20,2f12.4)') "N III (57.3 um):",(luminosity1)/fac,(luminosity1+luminosity2)/(0.292*2.05e37)
               case(6)
                  luminosity1 = zTemp
               case(7)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "O I (6300+6363):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(0.0086*2.05e37)
               case(8)
                  luminosity1 = zTemp
               case(9)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "O II (7320+7330):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(0.029*2.05e37)
               case(10)
                  luminosity1 = zTemp
               case(11)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "O II (3726+3729):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(2.03*2.05e37)
               case(12)
                  luminosity1 = zTemp
               case(13)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "O III (5007+4959):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(2.18*2.05e37)
               case(14)
                  luminosity3 = zTemp
                  write(*,'(a20,2f12.4)') "O III (4363):",(luminosity3)/1.e37,luminosity3/(0.0037*2.05e37)
               case(15)
                  luminosity1 = zTemp
               case(16)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "O III (52+88um):,",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/((1.06+1.22)*2.05e37)
               case(17)
                  luminosity1 = zTemp
                  write(*,'(a20,2f12.4)') "Ne II (12.8um):",(luminosity1)/fac,luminosity1/(0.195*2.05e37)

               case(18)
                  luminosity1 = zTemp
                  write(*,'(a20,2f12.4)') "Ne III (15.5um):",(luminosity1)/fac,luminosity1/(0.322*2.05e37)
               case(19)
                  luminosity1 = zTemp
               case(20)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "Ne III (3869+3968):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(0.085*2.05e37)
               case(21)
                  luminosity1 = zTemp
               case(22)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "S II (6716+6731):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(0.147*2.05e37)
               case(23)
                  luminosity1 = zTemp
               case(24)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "S II (4068+4076):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(0.008*2.05e37)
               case(25)
                  luminosity1 = zTemp
                  write(*,'(a20,2f12.4)') "S III (18.7um):",(luminosity1)/fac,luminosity1/(0.577*2.05e37)
               case(26)
                  luminosity1 = zTemp
               case(27)
                  luminosity2 = zTemp
                  write(*,'(a20,2f12.4)') "S III (9532+9069):",(luminosity1+luminosity2)/fac,&
                       (luminosity1+luminosity2)/(1.22*2.05e37)
               end select

               zTemp = 0.d0
            end do
         end if
      end if

    if(grid%geometry == "point") then
       write(mpiFilename,'(a, i4.4, a)') "point.grid"
       call writeAmrGrid(mpiFilename, .false., grid)
    end if


    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    thisThreadConverged = .false.    
    failed = .false.

     anyUndersampled = .false.
!     if(grid%geometry == "hii_test") then
!        minCrossings = 100
!     else if(grid%geometry == "lexington") then
!        if(grid%octreeroot%oned) then
!           minCrossings = 50000
!        else
!           minCrossings = 30000
!        end if
!        !     else if(singleMegaPhoto) then
!        !        minCrossings = 50000 
!     else if(.not. hOnly) then
!        minCrossings = 3000
!     else
!        !        minCrossings = 5000
!        minCrossings = 10
!     end if

   !Thaw - auto convergence testing I. Temperature, will shortly make into a subroutine
     maxDeltaT = -1.d30
     undersampled = .false.
     if (myRankGlobal /=0 ) then
        call findMaxFracTempChangeAndUndersampled(grid%octreeRoot, maxDeltaT, undersampled)
     endif
     call MPI_ALLREDUCE(maxDeltaT, tempDouble, 1, MPI_DOUBLE_PRECISION, MPI_MAX, localWorldCommunicator, ierr)
     maxDeltaT = tempDouble
     call MPI_ALLREDUCE(undersampled, tempLogical, 1, MPI_LOGICAL, MPI_LOR, localWorldCommunicator, ierr)
     undersampled = tempLogical

     undersampled = .false.

     converged = .false.
     if (myrankWorldGlobal == 0) write(*,*) "Maximum fractional change in T is ",maxDeltaT
     if (maxDeltaT < 0.05d0) converged = .true.
     if ((maxiter > 1).and.(niter <= 8).and.variableDustSublimation) converged = .false.
     if (nIter >= maxIter) then
        if (myRankWorldGlobal == 0) write(*,*) "Exceeded maxiter iterations, forcing convergence"
        converged = .true.
     else if (undersampled) then
        converged = .false.
        nTotalMonte = nTotalMonte * 2
        if (myrankWorldGlobal == 0) write(*,*) "Undersampled cells found. Increasing nMonte to ",nTotalMonte
     endif

     if(uv_vector) then
!        print *, "REDUCING UV_VEC NOISE"
        if(converged .and. firstConverged) then
!           print *, "forcing another iteration to reduce UV vector noise", nIter
!           converged = .false.
           firstConverged = .false.
           nTotalMonte = nTotalMonte * 10
        elseif(.not. firstConverged) then
           print *, "Allowing end to calculation at", nIter
           converged = .true.
        endif
     endif

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)
     if (.not.loadBalancingThreadGlobal) deallocate(octalArray)    
     
     if(uv_vector .or. singlemegaphoto .or. dumpRegularVTUS) then
        write(mpiFilename,'(a, i4.4, a)') "photo_", grid%iDump,".grid"
        call writeAmrGrid(mpiFilename, .false., grid)
        write(mpiFilename,'(a, i4.4, a)') "photo", nIter,".vtk"!
        call writeVtkFile(grid, mpiFilename, &
             valueTypeString=(/"rho          ", "HI           " , "temperature  ", "uvvec        ", &
             "crossings    "/))
     end if

!, &
!             "OI           ","HeI          ","HeII         ", "OI           ", "OII          "/))

!          "hydrovelocity","sourceCont   ","pressure     ", &
!          "crossings    ", &
!          "chiline      ", &
!          "dust         ", &
!          "diff         "/))

!     if(singleMegaPhoto) then



!        write(mpiFilename,'(a, i4.4, a)') "photo", nIter,".vtk"!

!     if(hydrodynamics) then
 !       call writeVtkFile(grid, mpiFilename, &
  !           valueTypeString=(/"rho          ","logRho       ", "HI           " , "temperature  "/))
!             "OI           ","HeI          ","HeII         ", "OI           ", "OII          ", &
!             "OIII         ","dust         ","dust2        " /))
!     end if
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

     if ((myRankGlobal /= 0).and.(.not.loadBalancingThreadGlobal)) then
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
        if(doPhotoRefine .and. 0 == 1) then!
           dprCounter = dprCounter + 1
           if(dprCounter == 4) then
              call setAllUnchanged(grid%octreeRoot)
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
!     call  dumpValuesAlongLine(grid, "radial_converged.dat", VECTOR(1.d0,0.d0,0.0d0), &
!          VECTOR(grid%octreeRoot%subcellSize, 0.d0, 0.d0),1000)
  
!        call writeVtkFile(grid, "temp.vtk", &
!             valueTypeString=(/"rho          ", "HI           " , "temperature  "/))


     call torus_mpi_barrier
     call MPI_BUFFER_DETACH(buffer,bufferSize, ierr)
     if(optimizeStack) then
        deallocate(photonPacketStack)     
        deallocate(currentStack)     
        deallocate(toSendStack)     
        deallocate(buffer)
     endif
  enddo
!  write(*,*) "set ",myHydroSetGlobal, " has converged ",niter,maxiter, maxdeltaT

 call MPI_BARRIER(MPI_COMM_WORLD, ierr)

 call MPI_TYPE_FREE(mpi_vector, ierr)
 call MPI_TYPE_FREE(mpi_photon_stack, ierr)

 call freeOctalSpectrum(grid%octreeRoot)
 deallocate(nSaved)
 deallocate(nEscapedArray)
 if(allocated(buffer)) then
    deallocate(buffer)
 end if
666 continue
end subroutine photoIonizationloopAMR


!!$recursive subroutine defineDiffusionZone(grid, maxRadius)
!!$  use mpi
!!$  use inputs_mod,only : resetDiffusion
!!$  type(GRIDTYPE) :: grid
!!$  integer :: nDiff
!!$  real(double) :: maxRadius, temp
!!$  integer :: ierr
!!$  resetDiffusion = .true.
!!$  ndiff = 0
!!$  call defineDiffusionOnKappap(grid, grid%octreeRoot, 1., nDiff)
!!$  maxRadius = 0.d0
!!$  call findMaxRadiusToDiffusionCell(grid%octreeRoot, maxRadius)
!!$  call MPI_ALLREDUCE(maxRadius, temp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
!!$       amrCommunicator, ierr)
!!$  maxRadius = temp
!!$  if (myrankWorldGlobal == 1) write(*,*) "Maximum diffusion radius ",maxRadius
!!$  call setDiffusionZoneOnRadius(grid%octreeRoot, VECTOR(0.d0, 0.d0, 0.d0), maxRadius)
!!$end subroutine defineDiffusionZone

recursive subroutine findMaxRadiusToDiffusionCell(thisOctal, maxRadius)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  real(double) :: maxRadius
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call findMaxRadiusToDiffusionCell(child, maxRadius)
              exit
           end if
        end do
     else
        if (thisOctal%diffusionApprox(subcell)) then
           maxRadius = max(maxRadius, modulus(subcellCentre(thisOctal, subcell) &
                - globalSourceArray(1)%position))
        endif
     end if
  end do
end subroutine findMaxRadiusToDiffusionCell

recursive subroutine  setDiffusionZoneOnRadius(thisOctal, position, maxRadius)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  type(VECTOR) :: position
  integer :: i, subcell
  real(double) :: maxRadius
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call setDiffusionZoneOnRadius(child, position, maxRadius)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
        if (modulus(subcellCentre(thisOctal, subcell) - position) < maxRadius) then
           thisOctal%diffusionApprox(subcell) = .true.
        endif
     end if
  end do
end subroutine setDiffusionZoneOnRadius

recursive subroutine  setPhotoionIsothermal(thisOctal)
  use inputs_mod, only : tMinGlobal
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call setPhotoionIsoThermal(child)
              exit
           end if
        end do
     else
        thisOctal%temperature(subcell) = tMinGlobal
        if (.not.associated(thisOctal%kappaTimesFlux)) then
           allocate(thisOctal%kappaTimesFlux(1:thisOctal%maxChildren))
        endif
!        thisOctal%kappaTimesFlux(subcell) = VECTOR(0.d0,0.d0,0.d0)
     end if
  end do
end subroutine setPhotoionIsothermal

recursive subroutine  zeroNfreq(thisOctal)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call zeroNfreq(child)
              exit
           end if
        end do
     else
        if (.not.associated(thisOctal%nfreq)) allocate(thisOctal%nFreq(1:thisOctal%maxChildren))
        thisOctal%nFreq(subcell) = 0
     end if
  end do
end subroutine zeroNfreq

recursive subroutine  freeOctalSpectrum(thisOctal)
  TYPE(OCTAL),pointer :: thisOctal
  TYPE(OCTAL),pointer :: child
  integer :: i, subcell
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call freeOctalSpectrum(child)
              exit
           end if
        end do
     else
        if (associated(thisOctal%spectrum)) then
           deallocate(thisOctal%spectrum)
           thisOctal%spectrum => null()
        endif
     end if
  end do
end subroutine freeOctalSpectrum


SUBROUTINE toNextEventPhoto(grid, rVec, uHat,  escaped,  thisFreq, nLambda, lamArray, photonPacketWeight, epsOverDeltaT, &
     nfreq, freq, dfreq, tPhoton, tLimit, crossedMPIboundary, newThread, sourcePhoton, crossedPeriodic, MiePhase, nMuMie, &
     movedCells)

  use inputs_mod, only : periodicX, periodicY, periodicZ, amrgridcentrey
  use inputs_mod, only : amrgridcentrez, radpressuretest, nDensity, timeDependentRT
  use mpi

   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec,uHat, octVec,thisOctVec, tvec, oldRvec, upperBound, lowerBound, iniVec
   type(PHASEMATRIX) :: miePhase(:,:,:)
   integer :: nMuMie
   logical :: movedCells
   type(OCTAL), pointer :: thisOctal
   type(OCTAL),pointer :: oldOctal
   real(double) :: tPhoton, tLimit, epsOverDeltaT
   real(double) :: photonMomentum
   integer, intent(out) :: newThread
   real(double) :: photonPacketWeight
   integer :: nFreq
   real(double) :: freq(:), dfreq(:), tempDouble
   integer :: subcell
   real(oct) :: tval, tau, r
   real :: lamArray(:)
   integer :: nLambda
   logical :: stillinGrid
   logical :: escaped
   real(double) :: kappaScaDb, kappaAbsDb, kappaRoss, wallDist
   real(double), parameter ::  gamma = 3.d0
   real(oct) :: thisTau
   real(oct) :: thisFreq
   real(oct) :: thisLam
   integer :: iLam
   logical ::inFlow
   real(double) :: kappap
!   real :: lambda
!   integer :: ilambda
   logical, intent(out) :: crossedMPIboundary
   type(OCTAL), pointer :: nextOctal
   integer :: nextSubcell, i, nStillInGrid
   logical :: outofTime, neverInteracted
   logical :: sourcePhoton, crossedPeriodic
   logical :: usedMRW

    Photonmomentum = epsOverDeltaT / cSpeed

    stillinGrid = .true.
    escaped = .false.
    
    movedCells = .false.
    crossedMPIboundary = .false.
    outOfTime = .false.

    neverInteracted = .true.

    thisLam = (cSpeed / thisFreq) * 1.e8

    call locate(lamArray, nLambda, real(thisLam), iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam,thisLam
    endif

! select an initial random tau and find distance to next cell

    call findSubcellTD(rVec, grid%octreeRoot,thisOctal, subcell)

    call randomNumberGenerator(getDouble=r)
    tau = -log(1.0-r)

!    write(*,*) "calling distancetocellboundary with ",cylindricalHydro
!    write(*,*) "rVec ",rVec
!    write(*,*) "uHat ",uHat

!    if (rVec%z < 0.d0) then
!       write(*,*) "rVec ",rVec
!       write(*,*) "uHat ",uHat
!    endif
    call distanceToCellBoundary(grid, rVec, uHat, tval, thisOctal, subcell)
!    write(*,*) "tval ",tval
    
    octVec = rVec
    thisOctVec = rVec

    call locate(lamArray, nLambda, real(thisLam), iLam)

    call amrGridValues(grid%octreeRoot, octVec,  iLambda=iLam, lambda=real(thisLam), startOctal=thisOctal, &
         actualSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, rosselandKappa=kappaRoss, &
         grid=grid, inFlow=inFlow)
    kappaP = thisOctal%kappap(subcell)
    oldOctal => thisOctal


    tempDouble = dfreq(1)
    usedMRW = .false.
    i = 0
    call distanceToNearestWall(rVec, wallDist, thisOctal, subcell)

!    write(*,*) "ross tau to wall ",(thisOctal%rho(subcell)*kappaRoss*1.d10*wallDist)
    if (.not.radpressuretest) then
       if (wallDist > gamma/(thisOctal%rho(subcell)*kappaRoss*1.d10) .and. .not. cart2d) then
          call modifiedRandomWalk(grid, thisOctal, subcell, rVec, uHat, &
               freq, dfreq, nfreq, lamArray, nlambda, thisFreq, photonPacketWeight)
          usedMRW = .true.
          thisLam = (cspeed/thisFreq)*1.d8
          call locate(lamArray, nLambda, real(thisLam), iLam)
          call amrGridValues(grid%octreeRoot, octVec,  iLambda=iLam, lambda=real(thisLam), startOctal=thisOctal, &
               actualSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
               grid=grid, inFlow=inFlow)
          call distanceToCellBoundary(grid, rVec, uHat, tval, thisOctal, subcell)
!          write(*,*) "used MRW ",rVec
       endif
    endif


    if (radpressureTest) then
       if (thisOctal%rho(subcell)<(1d-3*nDensity*mHydrogen)) then
          kappaAbsDb = 0.d0
          kappaScaDb = 0.d0
       else
          if (thisLam < 2.d5) then
             kappaScaDb = 0.d0
             kappaAbsDb = (100.d0/thisOctal%subcellSize)
          endif
       endif
    endif
    if (inFlow) then
       thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
    else
       thisTau = 1.0e-28
    end if

! if tau > thisTau then the photon has traversed a cell with no interactions

    nStillInGrid = 0
   do while(stillinGrid .and. (tau > thisTau) .and. .not. outOfTime) 
       nStillInGrid = nStillInGrid + 1
!       if (nStillInGrid > 10000) then
!          write(*,*) myRankWorldGlobal," nStillInGrid ",nStillInGrid
!          write(*,*) "rVec ",rVec
!          write(*,*) "uHat ",uHat
!          write(*,*) "tVal ",tval
!          write(*,*) "xmin,xmax ",thisOctal%xmin, thisOctal%xmax
!          write(*,*) "thisTau ",thisTau
!          write(*,*) "tau ",tau
!       endif
! add on the distance to the next cell

       oldrVec = rVec
       rVec = rVec + (tVal+1.d-5*grid%halfSmallestSubcell) * uHat ! rvec is now at face of thisOctal, subcell
       movedCells = .true.

       tPhoton = tPhoton + (tVal * 1.d10) / cSpeed
       if ((tPhoton > tLimit).and.timeDependentRT) then
          escaped = .true.
!          write(*,*) myrankGlobal, &
!               " big photon packet escaped due to tlimit. bug"
          outOfTime = .true.
       endif


       if(cart2d) then
          if(grid%octreeroot%oned) then
             if(abs(rVec%y - amrgridcentrey) > grid%halfsmallestsubcell &
                  .or. abs(rVec%z - amrgridcentrez) > grid%halfsmallestsubcell) then
                stillInGrid = .false.
                escaped = .true.
             end if
          else
             if(abs(rVec%y - amrgridcentrey) > (grid%halfsmallestsubcell)) then
!             if(abs(rVec%y - amrgridcentrey) > (grid%octreeroot%subcellsize)) then
                stillInGrid = .false.
                escaped = .true.
             end if
          end if
       end if
! check whether the photon has escaped from the grid
       if (inOctal(grid%octreeRoot,rVec) .and. stillInGrid) then
          nextOctal => thisOctal
          nextSubcell = subcell
          call findSubcellLocal(rVec, nextOctal, nextSubcell)
       else

          if(.not. crossedPeriodic .and. neverInteracted .and. .not. cart2d) then
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
          stillInGrid = .false.
          escaped = .true. 
          
       endif
       
       octVec = rVec
       
! check whether the photon has  moved across an MPI boundary
        
      if (stillInGrid) then

          if (.not.octalOnThread(nextOctal, nextsubcell, myRankGlobal)) then
             stillInGrid = .false.
             escaped = .true.
             crossedMPIboundary = .true.
             movedCells = .true.
             newThread = nextOctal%mpiThread(nextSubcell)
          endif
       endif

! update the distance grid

       if (.not.outOfTime)  then
          tVec = octVec - tval*uHat
          call updateGrid(grid, thisOctal, subcell, thisFreq, tVal, photonPacketWeight, ilam, nfreq, freq, &
               sourcePhoton, uHat, tVec,"first", miePhase, nMuMie)!octVec)
       endif

         
       if (stillinGrid) then
          
          thisOctal => nextOctal
          subcell = nextSubcell


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

!          if (rVec%z < 0.d0) then
!             write(*,*) "rVec ",rVec
!             write(*,*) "uHat ",uHat
!          endif
          call distanceToCellBoundary(grid, rVec, uHat, tval, thisOctal, subcell)



          if (radpressureTest) then
             if (thisOctal%rho(subcell)<(1d-3*nDensity*mHydrogen)) then
                kappaAbsDb = 0.d0
                kappaScaDb = 0.d0
             else
                if (thisLam < 2.d5) then
                   kappaScaDb = 0.d0
                   kappaAbsDb = (100.d0/thisOctal%subcellSize)
                endif
             endif
          endif


          octVec = rVec

! calculate the optical depth to the next cell boundary

          if (inFlow) then
             thisTau = dble(tVal) * (kappaAbsdb + kappaScadb)
          else
             thisTau = 1.0e-28
          end if

          if (tVal == 0.0d0) then
             escaped = .true.
             write(*,*) "tval zero error!"
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



       if (.not.inFlow) kappaAbsdb =0.0d0


! update the distance grid

       if (thisTau > 0.d0) then

!          lambda = cSpeed*1.e8/thisFreq
!          call locate(lamArray, nLambda, lambda, ilambda)
          if (.not.outOfTime) then 
             call updateGrid(grid, thisOctal, subcell, thisFreq, &
                  dble(tval)*dble(tau)/thisTau, photonPacketWeight, ilam, nfreq, freq, sourcePhoton, uHat, rVec, &
                  "sec", miePhase, nMuMie)
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
       if ((tPhoton > tLimit).and.(timeDependentRT)) then
          escaped = .true.
          outOfTime = .true.        
       endif

       if (.not.inOctal(grid%octreeRoot, rVec)) then  ! this is only needed due to floating point boundary issues
          escaped = .true.
          goto 666
       endif

       octVec = rVec
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

  recursive subroutine adjustChiline(thisOctal)
    TYPE(OCTAL),pointer :: thisOctal
    TYPE(OCTAL),pointer :: child
    integer :: i, subcell

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call adjustChiline(child)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          if (thisOctal%ncrossings(subcell) > 0) &
          thisOctal%chiline(subcell) = thisOctal%chiline(subcell) / dble(thisOctal%nCrossings(subcell))
       end if
    end do
  end subroutine adjustChiline

  recursive subroutine findMaxFracTempChangeAndUndersampled(thisOctal, frac, undersampled)
    use inputs_mod, only : minCrossings
    TYPE(OCTAL),pointer :: thisOctal
    TYPE(OCTAL),pointer :: child
    real(double) :: frac
    real         :: thisFrac
    integer :: i, subcell
    logical :: undersampled

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findMAxFracTempChangeAndUndersampled(child, frac, undersampled)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) cycle
          thisFrac = thisOctal%temperature(subcell) - thisOctal%TlastIter(subcell) / thisOctal%temperature(subcell)
          frac = max(abs(dble(thisFrac)), frac)
          thisOctal%tLastIter(subcell) = thisOctal%temperature(subcell)
          if (thisOctal%nCrossings(subcell) < minCrossings) undersampled = .true.
       endif
    end do
  end subroutine findMaxFracTempChangeAndUndersampled




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
    use inputs_mod, only : nDustType, grainFrac
    type(GRIDTYPE) :: grid
    TYPE(OCTAL), POINTER  :: thisOctal 
    TYPE(OCTAL), POINTER  :: child
    INTEGER :: i, j

    IF ( thisOctal%nChildren > 0 ) THEN
       ! call this subroutine recursively on each of its children
       DO i = 1, thisOctal%nChildren, 1
          child => thisOctal%child(i)
          CALL allocatePhotoIonAttributes(child, grid)
       END DO
    endif

      call allocateAttribute(thisOctal%oldFrac, thisOctal%maxChildren)
      thisOctal%oldFrac = 1.e-30
      call allocateAttribute(thisOctal%dustType, thisOctal%maxChildren)
      thisOctal%dustType = 1
      call allocateAttribute(thisOctal%dustTypeFraction, thisOctal%maxChildren, nDustType)
      thisOctal%dustTypeFraction = 0.d0
      do j = 1, ndustType
         thisOctal%dustTypeFraction(:,j) = grainFrac(j)
      enddo

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

      call allocateAttribute(thisOctal%corner,thisOctal%maxchildren)
      call allocateAttribute(thisOctal%boundaryPartner,thisOctal%maxChildren)
      call allocateAttribute(thisOctal%boundaryCondition,thisOctal%maxchildren)
      call allocateAttribute(thisOctal%edgeCell,thisOctal%maxchildren)
      call allocateAttribute(thisOctal%correction,thisOctal%maxchildren)

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
          thisOctal%chiline(subcell) = 0.d0
          thisOctal%radiationMomentum(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
          if (.not.associated(thisOctal%kappaTimesFlux)) then
             allocate(thisOctal%kappaTimesFlux(1:thisOctal%maxChildren))
          endif
          if (.not.associated(thisOctal%UVvector)) then
             allocate(thisOctal%UVvector(1:thisOctal%maxChildren))
          endif
          if (.not.associated(thisOctal%UVvectorPlus)) then
             allocate(thisOctal%UVvectorPlus(1:thisOctal%maxChildren))
          endif
          if (.not.associated(thisOctal%UVvectorMinus)) then
             allocate(thisOctal%UVvectorMinus(1:thisOctal%maxChildren))
          endif
          thisOctal%kappaTimesFlux(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
          thisOctal%UVvector(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
          thisOctal%UVvectorPlus(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
          thisOctal%UVvectorMinus(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
       endif
    enddo
  end subroutine zeroDistanceGrid


  recursive subroutine enforceMinRho(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call enforceMinRho(child)
                exit
             end if
          end do
       else
          if(thisOctal%rho(subcell)/mhydrogen < 1.d0) then
             thisOctal%rho(subcell) = mhydrogen
          endif
       endif
    enddo
  end subroutine enforceMinRho

  recursive subroutine enforceMinRhodisc(thisOctal)
    use inputs_mod, only : sourcepos, maxdepthamr, amrgridsize
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  type(vector) :: rvec
  real(double) :: dx
  logical, save :: message=.true.
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call enforceMinRho(child)
                exit
             end if
          end do
       else
          if(thisOctal%rho(subcell)/mhydrogen < 1.d0) then
             thisOctal%rho(subcell) = mhydrogen
          endif
          dx = amrgridsize/(2.d0**maxdepthamr)
          rvec = subcellCentre(thisoctal, subcell)
          if(modulus(rvec-sourcepos(1)) < 2.d0*dx) then
             if(message) then
                print *, "accreting"
                message = .false.
             endif
             thisOctal%rho(subcell) = 1.d2*mhydrogen
          endif


       endif
    enddo
  end subroutine enforceMinRhodisc

  recursive subroutine enforceMinRhofont(thisOctal)
    use inputs_mod, only : maxdepthamr, amrgridsize, sourcemass
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  real(double) :: dx, alpha, rg, cs
  type(vector) :: rvec

  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call enforceMinRhofont(child)
                exit
             end if
          end do
       else
          dx = amrgridsize/2.d0**maxdepthamr
          rvec = subcellcentre(thisoctal, subcell)
          if(abs(rvec%z) < dx) then
             cs = sqrt(kerg*10.0/mhydrogen)
             rg = bigG*sourceMass(1)/cs**2
             alpha = 3.d0/2.d0
             if(thisOctal%rho(subcell) < 1.d8*mhydrogen*(rvec%x/rg)**(-alpha)) then
                
                thisOctal%rho(subcell) = 1.d8*mhydrogen*(rvec%x/rg)**(-alpha)
             endif
          endif
       endif
    enddo
  end subroutine enforceMinRhofont

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
          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then             
!             v = 3.d0*cellVolume(thisOctal, subcell)/fourpi 
             v = cellVolume(thisOctal, subcell)
!             if(thisOctal%oneD) then
!                v = thisOctal%subcellSize**3 * 1.d30
!             else if(thisOctal%twoD) then
!                v = thisOctal%subcellSize**2 * 1.d30
!             else
!                v = 
!             end if
             hbeta = (10.d0**(-0.870d0*log10(thisOctal%temperature(subcell))+3.57d0)) * &
                  thisOctal%ne(subcell) * thisOctal%ionFrac(subcell, 2) * &
                  thisOctal%nh(subcell)*grid%ion(1)%abundance*1.d-25
             luminosity = luminosity + hbeta * (v*1.d30)
!             print *, "position ", subcellCentre(thisOctal, subcell)
!             print *, "luminosity ", luminosity
!             print *, "v ", v
!             print *, "hbeta ", hbeta
!             print *, "thisOctal%ne(subcell) ", thisOctal%nh(subcell)
!             print *, "thisOctal%nh(subcell) ", thisOctal%ne(subcell)
!             print *, "grid%ion(1)%abundance", grid%ion(1)%abundance
!             print *, "thisOctal%temperature(subcell)", thisOctal%temperature(subcell)
!             stop
          endif
       end if
    enddo
  end subroutine getHbetaluminosity

  recursive subroutine getRecombinationTime(thisOctal, tRecomb)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: tRecomb,thisT
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getRecombinationTime(child, tRecomb)
                exit
             end if
          end do
       else
          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then             
             thisT = 1.d0/(thisOctal%nh(subcell) * recombrate(globalIonArray(1), thisOctal%temperature(subcell)))
             tRecomb = min(thisT, tRecomb)


          endif
       end if
    enddo
  end subroutine getRecombinationTime


  recursive subroutine getIonizationTime(thisOctal, grid, tIonization, epsOverDeltaT)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: tIonization,thisT,photoionRate, recombinationRate
    integer :: k, iStart, iEnd, nIonizationStages
    real(double) :: chargeExchangeIonization, chargeExchangeRecombination
    integer :: subcell, i, iIon
    real(double) :: v,  epsOverDeltaT
    real :: temperature
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getIonizationTime(child, grid, tIonization, epsOverDeltaT)
                exit
             end if
          end do
       else
          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then             

             if (.not.thisOctal%ghostCell(subcell)) then
             v = cellVolume(thisOctal, subcell)
             k = 1
             temperature = thisOctal%temperature(subcell)
             

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


             end do
                thisT = 1.d0/max(abs(recombinationRate - photoIonRate),1.d-20)
                tIonization = min(thisT, tIonization)
                thisOctal%biasline3d(subcell) = thisT
             endif
          endif
       end if
    enddo
  end subroutine getIonizationTime

  recursive subroutine getThermalTime(thisOctal, grid, tThermal, epsOverDeltaT)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real(double) :: tThermal,thisT
    integer :: subcell, i
    real(double) :: epsOverDeltaT
    real(double) :: hHeating, heHeating, dustHeating, totalHeating, totalCooling
    real(double) :: mu, energy
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getThermalTime(child, grid, tThermal, epsOverDeltaT)
                exit
             end if
          end do
       else
          if (octalOnThread(thisOctal, subcell, myRankGlobal)) then             

             if (.not.thisOctal%ghostCell(subcell)) then

                call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)

                mu = returnMu(thisOctal, subcell, grid%ion, grid%nion)
                energy = 1.5d0*(thisOctal%rho(subcell)/(mu*mhydrogen))*kerg*thisOctal%temperature(subcell)
                totalCooling = HHecooling(grid, thisOctal, subcell, thisOctal%temperature(subcell)) 

                thisT =  (energy / max(abs(totalCooling-totalHeating),1.d-20))
                tThermal = min(thisT, tThermal)
             endif

          endif
       end if
    enddo
  end subroutine getThermalTime


  subroutine calculateIonizationBalance(grid, thisOctal, epsOverDeltaT, augerArray)
    use inputs_mod, only : xraycalc, useionparam
    type(gridtype) :: grid
    type(AUGER) :: augerArray(5, 5, 10)
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT
    integer :: subcell
    
    do subcell = 1, thisOctal%maxChildren
       
 !      if (hydrodynamics) then
!          if (thisOctal%ghostCell(subcell)) cycle
  !     endif

       if (.not.thisOctal%hasChild(subcell)) then

          if (thisOctal%inflow(subcell)) then
             if (.not.thisOctal%undersampled(subcell)) then
                if(xraycalc .and. .not. useionparam) then
                   call solveIonizationBalance_xray(grid, thisOctal, subcell, thisOctal%temperature(subcell), &
                        epsOverdeltaT, augerArray)
                else
                   call solveIonizationBalance(grid, thisOctal, subcell, thisOctal%temperature(subcell), epsOverdeltaT)
                end if
             else

                thisOctal%ionFrac(subcell, 1) = 1.d0
                thisOctal%ionFrac(subcell, 2) = 1.d-30
!                if(thisOctal%nCrossings(subcell) /= 0) then
!                   write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
!                end if
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
!             if (.not.thisOctal%undersampled(subcell)) then
                call solveIonizationBalanceTimeDep(grid, thisOctal, subcell, thisOctal%temperature(subcell), epsOverdeltaT, deltaT)
                
!             else
!                thisOctal%ionFrac(subcell, 1) = 1.d0
!                thisOctal%ionFrac(subcell, 2) = 1.d-30
!                if(thisOctal%nCrossings(subcell) /= 0) then
!                 !write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
!                end if
!             endif
          endif
       endif
    enddo
  end subroutine calculateIonizationBalanceTimeDep

  subroutine quickThermalCalc(thisOctal)
    use mpi
    use inputs_mod, only : tMinGlobal
    type(OCTAL), pointer :: thisOctal
    integer :: subcell


    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then
          if(octalOnThread(thisOCtal, subcell, myRankGlobal)) then
             thisOctal%temperature(subcell) = real(dble(tminGlobal) + (10000.d0-dble(tminGlobal)) * thisOctal%ionFrac(subcell,2))
          end if
       endif
       
    enddo
  end subroutine quickThermalCalc

!!$  subroutine calculateThermalBalance(grid, thisOctal, epsOverDeltaT)
!!$    use mpi
!!$    use inputs_mod, only : tMinGlobal, hydrodynamics
!!$    type(gridtype) :: grid
!!$    type(octal), pointer   :: thisOctal
!!$    real(double) :: epsOverDeltaT
!!$    real(double) :: totalHeating
!!$    integer :: subcell
!!$    logical :: converged, found
!!$    real :: t1, t2, tm
!!$    real(double) :: y1, y2, ym, Hheating, Heheating, dustHeating
!!$    real :: deltaT
!!$    real :: underCorrection = 1.
!!$    integer :: nIter
!!$    
!!$
!!$    do subcell = 1, thisOctal%maxChildren
!!$
!!$       if (.not.thisOctal%hasChild(subcell)) then
!!$
!!$ !THAW - trying to root out mpi things
!!$          if(octalOnThread(thisOctal, subcell, myRankGlobal)) then
!!$          
!!$             if (hydrodynamics) then
!!$                if (thisOctal%ghostCell(subcell)) cycle
!!$             endif
!!$
!!$
!!$
!!$             if (thisOctal%inflow(subcell)) then
!!$                
!!$!             write(*,*) thisOctal%nCrossings(subcell),thisOctal%undersampled(subcell)
!!$                
!!$                if (.not.thisOctal%undersampled(subcell)) then
!!$                   
!!$                   
!!$                   call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
!!$!                if (thisOctal%ionFrac(subcell,1) > 0.9d0) write(*,*) "total heating ",totalheating
!!$                   if (totalHeating < 1.d-30) then
!!$                      thisOctal%temperature(subcell) = tminGlobal
!!$                   else
!!$                      nIter = 0
!!$                      converged = .false.
!!$                      
!!$                      t1 = tMinGlobal
!!$                      t2 = 30000.
!!$!
!!$                      
!!$                      found = .true.
!!$                      
!!$                      if (found) then
!!$                         y1 = (HHecooling(grid, thisOctal, subcell, t1) &
!!$                         - totalHeating)
!!$                         y2 = (HHecooling(grid, thisOctal, subcell, t2) &
!!$                         - totalHeating)
!!$                         if (y1*y2 > 0.d0) then
!!$                            if (HHecooling(grid, thisOctal, subcell, t1) > totalHeating) then
!!$                               tm = t1
!!$!                               write(*,*) "Cell set to cold. Cooling is ", &
!!$!                                 HHecooling(grid, thisOctal, subcell, t1), " heating is ", totalHeating
!!$                            else
!!$                               tm  = t2
!!$                            write(*,*) "Cell set to hot. Cooling is ", &
!!$                                 HHecooling(grid, thisOctal, subcell, t1), " heating is ", totalHeating
!!$                            endif
!!$                            converged = .true.
!!$                         endif
!!$                         if (totalHeating < 1.d-30) then
!!$                            write(*,*) "no heating, set to cold ",totalheating
!!$                            tm = t1
!!$                            converged = .true.
!!$                         endif
!!$                         
!!$! Find root of heating and cooling by bisection
!!$                         
!!$                         do while(.not.converged)
!!$                            tm = 0.5*(t1+t2)
!!$                            y1 = (HHecooling(grid, thisOctal, subcell, t1) &
!!$                            - totalheating)
!!$                            y2 = (HHecooling(grid, thisOctal, subcell, t2) &
!!$                            - totalheating)
!!$                            ym = (HHecooling(grid, thisOctal, subcell, tm) &
!!$                            - totalheating)
!!$                            
!!$                            if (y1*ym < 0.d0) then
!!$                               t1 = t1
!!$                               t2 = tm
!!$                            else if (y2*ym < 0.d0) then
!!$                               t1 = tm
!!$                               t2 = t2
!!$                            else
!!$                               converged = .true.
!!$                               tm = 0.5*(t1+t2)
!!$                               write(*,*) t1, t2, y1,y2,ym
!!$                            endif
!!$                            
!!$                            if (abs((t1-t2)/t1) .le. 1.e-2) then
!!$                               converged = .true.
!!$                            endif
!!$                            niter = niter + 1
!!$!                         if (myrankGlobal == 1) write(*,*) niter,tm, abs(t1-t2)/t1
!!$                         enddo
!!$                      endif
!!$                      deltaT = tm - thisOctal%temperature(subcell)
!!$                      thisOctal%temperature(subcell) = &
!!$                      max(thisOctal%temperature(subcell) + underCorrection * deltaT,tMinGlobal)
!!$!                                write(*,*) thisOctal%temperature(subcell), niter
!!$                   endif
!!$                else
!!$                !                write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
!!$                endif
!!$             endif
!!$          end if
!!$       endif
!!$      enddo
!!$  end subroutine calculateThermalBalance

  subroutine calculateThermalBalanceNew(grid, thisOctal, epsOverDeltaT)
    use mpi
    use inputs_mod, only : tMinGlobal
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT
    real(double) :: totalHeating
    integer :: subcell
    logical :: converged
    real(double) :: a, b, c, d, e, fa, fb, fc, r, s, p, q, tol, tol1, xm
    integer, parameter :: itmax = 100
    integer :: iter
    real(double), parameter :: eps = 3.d-8
    real :: t1, t2
    real(double) :: Hheating, Heheating, dustHeating, newT, deltaT
    real :: underCorrection = 0.5
    integer :: nIter
    tol = 1.d-2

    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then

 !THAW - trying to root out mpi things
          if(octalOnThread(thisOctal, subcell, myRankGlobal)) then
          
             if (thisOctal%inflow(subcell)) then
                
                
                if (.not.thisOctal%undersampled(subcell)) then
                   

                   
                   call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)

                   iter = 0
                   if ((totalHeating < HHecooling(grid, thisOctal, subcell, tMinGlobal)).or.(totalHeating < 1.d-30)) then
                      newT = tminGlobal
                   else if ((totalHeating > HHecooling(grid, thisOctal, subcell, 20000.))) then
                      newT = 20000.
                   else
                      nIter = 0
                      converged = .false.
                      
                      !try nearby temperatures to start the bracket
                      t1 = thisOctal%temperature(subcell) * 0.8
                      t2 = thisOctal%temperature(subcell) * 1.2
                      FA = (totalHeating - HHecooling(grid, thisOctal, subcell, t1)) !FUNC(A)
                      FB = (totalHeating - HHecooling(grid, thisOctal, subcell, t2)) !FUNC(B)

                      IF (FB*FA.GT.0.d0) then
                         t1 = tMinGlobal
                         t2 = 20000.
!                         write(*,*) "failed using nearby temps ",thisOctal%temperature(subcell)
                      else
!                         write(*,*) "success using nearby temps ",thisOctal%temperature(subcell)
                      endif

                      A = t1
                      B = t2
                      FA = (totalHeating - HHecooling(grid, thisOctal, subcell, real(a))) !FUNC(A)
                      FB = (totalHeating - HHecooling(grid, thisOctal, subcell, real(b))) !FUNC(B)
                      IF (FB*FA.GT.0.d0) write(*,*) 'Root must be bracketed for ZBRENT.',totalHeating,fa,fb
                      FC = FB
                      DO ITER=1,ITMAX
                         IF(FB*FC.GT.0.d0) THEN
                            C=A
                            FC=FA
                            D=B-A
                            E=D
                         ENDIF
                         IF(ABS(FC).LT.ABS(FB)) THEN
                            A=B
                            B=C
                            C=A
                            FA=FB
                            FB=FC
                            FC=FA
                         ENDIF
                         TOL1=2.d0*EPS*ABS(B)+0.5d0*TOL
                         XM=0.5d0*(C-B)
                         IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.d0)THEN
                            newT = B
                            exit
                         ENDIF
                         IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
                            S=FB/FA
                            IF(A.EQ.C) THEN
                               P=2.*XM*S
                               Q=1.-S
                            ELSE
                               Q=FA/FC
                               R=FB/FC
                               P=S*(2.d0*XM*Q*(Q-R)-(B-A)*(R-1.d0))
                               Q=(Q-1.d0)*(R-1.d0)*(S-1.d0)
                            ENDIF
                            IF(P.GT.0.d0) Q=-Q
                            P=ABS(P)
                            IF(2.d0*P .LT. MIN(3.d0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                               E=D
                               D=P/Q
                            ELSE
                               D=XM
                               E=D
                            ENDIF
                         ELSE
                            D=XM
                            E=D
                         ENDIF
                         A=B
                         FA=FB
                         IF(ABS(D) .GT. TOL1) THEN
                            B=B+D
                         ELSE
                            B=B+SIGN(TOL1,XM)
                         ENDIF
                         FB = (totalHeating - HHecooling(grid, thisOctal, subcell, real(b))) !FUNC(B)
                      enddo
!                      write(*,*) "finished after ",iter, " iterations"
                      if (ABS(XM).gt.TOL1) write(*,*) 'ZBRENT exceeding maximum iterations.', xm, tol1
                      newT = B


                   END if
                   deltaT = newT- dble(thisOctal%temperature(subcell))
                   thisOctal%temperature(subcell) = thisOctal%temperature(subcell) + underCorrection * real(deltaT)
                   !                   write(*,*) "thermal balance new converged after ",iter
                endif
             endif
          end if
       endif
    enddo
  end subroutine calculateThermalBalanceNew

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

!             if (.not.thisOctal%undersampled(subcell)) then
             if (.not.thisOctal%ghostCell(subcell)) then
                
                call calculateEquilibriumTemperature(grid, thisOctal, subcell, epsOverDeltaT, Teq)

                call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)

                mu = returnMu(thisOctal, subcell, grid%ion, grid%nion)
                energy = 1.5d0*(thisOctal%rho(subcell)/(mu*mhydrogen))*kerg*thisOctal%temperature(subcell)
                totalCooling = HHecooling(grid, thisOctal, subcell, thisOctal%temperature(subcell)) 

                thermalTime = (energy / abs(totalCooling))

!                if (deltaTime > thermalTime) then
!                   thisOctal%temperature(subcell) = real(Teq)
!                else

                   nTimes = 10
                
                   smallDeltaT = deltaTime /dble(nTimes)
                   
                   do i = 1, nTimes
                      energy = energy + (totalHeating - totalcooling)*smallDeltaT
                      thisOctal%temperature(subcell) = real((2.d0/3.d0)*energy*(mu*mhydrogen)/(thisOctal%rho(subcell)*kerg))
                      totalCooling = HHecooling(grid, thisOctal, subcell, thisOctal%temperature(subcell)) 
                   enddo
!                endif

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
    real(double) :: kappap
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

!    call returnKappa(grid, thisOctal, subcell, kappap=kappap, atthistemperature=temperature)
    kappaP = returnKappaP(thisOctal, subcell, dble(temperature))
    dustCooling = fourPi * kappaP * (stefanBoltz/pi) * temperature**4

    coolingRate = coolingRate + dustCooling

  end function HHeCooling

  subroutine updateGrid(grid, thisOctal, subcell, thisFreq, distance, &
       photonPacketWeight, ilambda, nfreq, freq, sourcePhoton, uHat, rVec, flag, miePhase, nMumie)
    use inputs_mod,only : dustOnly, radPressureTest, UV_vector, UV_high, UV_low, nDensity
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    type(VECTOR) :: uHat, rVec, uHatDash, rHat, zHat, rVecTemp
    integer :: nDist
    real(double) :: dS
    integer :: nFreq, iFreq
    real(double) :: freq(:)
    integer :: subcell
    real(double) :: thisFreq, distance, kappaAbs,kappaAbsDust, kappaSca, kappaExt
    integer :: ilambda
    real(double) :: photonPacketWeight
    integer :: i 
    real(double) :: fac, xSec,gfac
    logical :: sourcePhoton
    character(len=*) :: flag
    integer :: nMuMie
    type(PHASEMATRIX) :: miePhase(:,:,:)

    
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

    call returnkappa(grid, thisoctal, subcell, ilambda=ilambda, kappaabsdust=kappaabsdust, kappaabs=kappaabs, kappaSca=kappaSca)
    kappaExt = kappaAbs + kappaSca


    thisoctal%distancegrid(subcell) = thisoctal%distancegrid(subcell) &
         + dble(distance) * dble(kappaabsdust) * photonPacketWeight

    uHatDash = uHat
    if (thisOctal%twoD .and. .not. cart2d) then
       rHat = VECTOR(rVec%x, rVec%y, 0.d0)
       call normalize(rHat)
       zHat = VECTOR(0.d0, 0.d0, 1.d0)
       uHatDash = VECTOR(rHat.dot.uHat, 0.d0, zHat.dot.uHat)
    endif

    if (thisOctal%oneD) then
       rHat = VECTOR(rVec%x, rVec%y, rVec%z)
       call normalize(rHat)
       uHatDash = VECTOR(rHat.dot.uHat, 0.d0, 0.d0)
    endif

    if (radpressureTest) then
       if (thisOctal%rho(subcell)<(1d-3*nDensity*mHydrogen)) then
          kappaext = 0.d0
       else
          if (grid%lamArray(ilambda) < 2.d5) then
             kappaExt = (100.d0/thisOctal%subcellSize)
          endif
       endif
    endif

! Kludge to avoid compiler warning
    if (.false.) write(*,*) "insubcell ",inSubcell(thisOctal, subcell, rvec), flag, nmumie

    rVecTemp = rVec
    nDist = 50
    dS = distance / dble(nDist)
    gFac = miePhase(1,iLambda, 1)%gfac
    gFac = 0.d0
    do i = 1, nDist
       rVecTemp = rVec + (dble(i-1)* ds) * uHat
       rHat = VECTOR(rVecTemp%x, rVecTemp%y, rVecTemp%z)
       call normalize(rHat)
       uHatDash = uHat
       if (thisOctal%twoD .and. .not. cart2d) then
          zHat = VECTOR(0.d0, 0.d0, 1.d0)
          uHatDash = VECTOR(rHat.dot.uHat, 0.d0, zHat.dot.uHat)
       endif

       if (thisOctal%oneD) then
          rHat = VECTOR(rVec%x, rVec%y, rVec%z)
          uHatDash = VECTOR(rHat.dot.uHat, 0.d0, 0.d0)
       endif
       thisoctal%kappaTimesFlux(subcell) = thisoctal%kappaTimesFlux(subcell) &
            + (dS * dble(kappaAbs + kappaSca*(1.d0-gfac)) * photonPacketWeight)*uHatDash
    enddo

       


!     thisoctal%kappaTimesFlux(subcell) = thisoctal%kappaTimesFlux(subcell) &
!         + (dble(distance) * dble(kappaExt) * photonPacketWeight)*uHatDash


    if(uv_vector) then
!       if(thisFreq > (2.99792458d8/100.d-9) .and. thisFreq < (2.99792458d8/10.d-9)) then
       if(thisFreq*hcgs > (UV_low*evToErg) .and. thisFreq*hcgs < (UV_high*evToErg)) then
          thisOctal%UVvector(subcell) = thisOctal%UVvector(subcell)&
               + (dble(distance) * photonpacketweight*uHatdash)

!          if(uhatdash 
          if(uhatdash%x > 0.d0) then
             thisOctal%UVvectorPlus(subcell)%x = thisOctal%UVvectorPlus(subcell)%x&
                  + (dble(distance) * photonpacketweight*uHatdash%x)
          else
             thisOctal%UVvectorMinus(subcell)%x = thisOctal%UVvectorMinus(subcell)%x&
                  + (dble(distance) * photonpacketweight*uHatdash%x)
          endif

          if(uhatdash%y > 0.d0) then
             thisOctal%UVvectorPlus(subcell)%y = thisOctal%UVvectorPlus(subcell)%y&
                  + (dble(distance) * photonpacketweight*uHatdash%y)
          else
             thisOctal%UVvectorMinus(subcell)%y = thisOctal%UVvectorMinus(subcell)%y&
                  + (dble(distance) * photonpacketweight*uHatdash%y)
          endif

          if(uhatdash%x > 0.d0) then
             thisOctal%UVvectorPlus(subcell)%z = thisOctal%UVvectorPlus(subcell)%z&
                  + (dble(distance) * photonpacketweight*uHatdash%z)
          else
             thisOctal%UVvectorMinus(subcell)%z = thisOctal%UVvectorMinus(subcell)%z&
                  + (dble(distance) * photonpacketweight*uHatdash%z)
          endif


          !            + (dble(distance) *dble(kappaExt)* photonPacketWeight)*uHatDash
       end if
    end if
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
  real(double), allocatable :: newFrac(:)
  real(double), allocatable :: xplus1overx(:)
  real(double), parameter :: underCorrection = 0.6d0

  v = cellVolume(thisOctal, subcell)
  k = 1
  allocate(newFrac(1:SIZE(thisOctal%ionFrac,2)))
  newFrac = thisOctal%ionFrac(subcell,:)
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
!     thisOctal%ionFrac(subcell, iStart:iEnd) = 1.d0
     newFrac(iStart:iEnd) = 1.d0
     do i = 1, nIonizationStages - 1
        iIon = iStart+i-1
!        thisOctal%ionFrac(subcell,iIon+1) = thisOctal%ionFrac(subcell,iIon) * xplus1overx(i)
        newFrac(iIon+1) = newFrac(iIon) * xplus1overx(i)
!        if ((myRankGlobal==1).and.grid%ion(iion)%species(1:2) =="H ") write(*,*) i,thisOctal%ionFrac(subcell,iIon)
     enddo
     
!     if (SUM(thisOctal%ionFrac(subcell,iStart:iEnd)) /= 0.d0) then
!        thisOctal%ionFrac(subcell,iStart:iEnd) = &
!             max(1.d-50,thisOctal%ionFrac(subcell,iStart:iEnd))/SUM(thisOctal%ionFrac(subcell,iStart:iEnd))
!     else
!        thisOctal%ionFrac(subcell,iStart:iEnd) = 1.d-50
!     endif

     if (SUM(newFrac(iStart:iEnd)) /= 0.d0) then
        newFrac(iStart:iEnd) = &
             max(1.d-50,newFrac(iStart:iEnd))/SUM(newFrac(iStart:iEnd))
     else
        newFrac(iStart:iEnd) = 1.d-50
     endif

        deallocate(xplus1overx)

        k = iEnd + 1
     end do

     thisOctal%ionFrac(subcell,:) = thisOctal%ionFrac(subcell,:) + underCorrection * (newFrac - thisOctal%ionFrac(subcell,:))
     deallocate(newFrac)
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
  real(double) :: photoIonRate, recombinationRate
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

        thisOctal%ionFrac(subcell,iion) = thisOctal%ionFrac(subcell,iion) + deltaT * (recombinationRate - photoIonRate)


        if(thisOctal%ionFrac(subcell, iIon) <= 0.d0) then
           thisOctal%ionFrac(subcell, iIon) = 1.d-50
        end if

     enddo

!     if (SUM(thisOctal%ionFrac(subcell,iStart:iEnd)) /= 0.d0) then
!        thisOctal%ionFrac(subcell,iStart:iEnd) = &
!             max(1.d-50,thisOctal%ionFrac(subcell,iStart:iEnd))/SUM(thisOctal%ionFrac(subcell,iStart:iEnd))
!     else
!        thisOctal%ionFrac(subcell,iStart:iEnd) = 1.d-50
!     endif
     


     thisOctal%ionFrac(subcell, iEnd) = 1.d0 - SUM(thisOctal%ionFrac(subcell,iStart:iEnd-1))
!     thisOctal%ionFrac(subcell, iend) = min(max(thisOctal%ionFrac(subcell,iend),ionFrac(iend)),1.d0)

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




subroutine dumpWhalenNormanTest(grid)
  use mpi
  use inputs_mod, only : quickThermal, readgrid
  type(gridtype) :: grid
  integer :: ier, ierr, i, j
  logical, save :: firstTime=.true.
  real(double) :: tempStorage(24), rhoBins(15), rhoDist(15)
  integer :: status(MPI_STATUS_SIZE)
  integer :: tag = 50
  real(double) :: minR, meanR, maxR, mIo, mNeu, Movd, KE, mom
  integer :: nRadii


  if(myRankGlobal == 0) then
     if(firstTime .and. .not. readGrid) then
        if(quickThermal) then
           open(unit=65, file="WN08simple_IFradius_torus_haworth.txt", &
                status="unknown", form="formatted",iostat=ier)
           open(unit=66, file="WN08simple_MassFn_torus_haworth.txt", &
                status="unknown", form="formatted",iostat=ier)
        else
           open(unit=65, file="WN08native_IFradius_torus_haworth.txt", &
                status="unknown", form="formatted",iostat=ier)
           open(unit=66, file="WN08native_MassFn_torus_haworth.txt", &
                status="unknown", form="formatted", iostat=ier)
        end if
        firstTime = .false.
        
        write (65, '(9(a7, 3x))') "time", "r_min", "r_mean", "r_max", "m_i", "m_n", "m_dens", "K.E.", "mom"
        write(66, '(16(a12, 3x))') "time", "0-10", "10-20", "20-40", "40-80", "80-160", "160-320", "320-640", "640-1280", &
              "1280-2560", "2560-5120","5120-10240", "10240-20480", "20480-40960" ,"40960-81920","81920+" 
     else
        if(quickThermal) then
           open(unit=65, file="WN08simple_IFradius_torus_haworth.txt", &
                status="unknown", position="append", form="formatted", iostat=ier)
           open(unit=66, file="WN08simple_MassFn_torus_haworth.txt", &
                status="unknown", position="append", form="formatted", iostat=ier)
        else
           open(unit=65, file="WN08native_IFradius_torus_haworth.txt", &
                status="unknown", position="append", form="formatted",iostat=ier)
           open(unit=66, file="WN08simple_MassFn_torus_haworth.txt", &
                status="unknown", position="append", form="formatted", iostat=ier)
        end if
     end if
  end if

  call  setupDensityBins(rhoBins)
  minR = 1.d30
  maxR = 0.d0
  meanR = 0.d0
  mIo = 0.d0
  mNeu = 0.d0
  Movd = 0.d0
  KE = 0.d0
  mom = 0.d0
  nRadii = 0
  rhoDist = 0.d0
  if(myRankGlobal /= 0) then
     call getTestValues(grid%octreeRoot, minR, meanR, maxR, mIo, mNeu, &
          Movd, KE, mom, nRadii, rhoBins, rhoDist)
     
     do i = 1, nhydrothreadsglobal
        if(i == myRankGlobal) then
           tempStorage(1) = minR
           tempStorage(2) = meanR
           tempStorage(3) = maxR
           tempStorage(4) = mIo
           tempStorage(5) = mNeu
           tempStorage(6) = Movd
           tempStorage(7) = KE
           tempStorage(8) = mom
           tempStorage(9) = nRadii
           tempStorage(10) = rhoDist(1)
           tempStorage(11) = rhoDist(2)
           tempStorage(12) = rhoDist(3)
           tempStorage(13) = rhoDist(4)
           tempStorage(14) = rhoDist(5)
           tempStorage(15) = rhoDist(6)
           tempStorage(16) = rhoDist(7)
           tempStorage(17) = rhoDist(8)
           tempStorage(18) = rhoDist(9)
           tempStorage(19) = rhoDist(10)
           tempStorage(20) = rhoDist(11)
           tempStorage(21) = rhoDist(12)
           tempStorage(22) = rhoDist(13)
           tempStorage(23) = rhoDist(14)
           tempStorage(24) = rhoDist(15)
           call MPI_SEND(tempStorage, 24, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)        
        end if
     end do
  end if

  if(myRankGlobal == 0) then
     do i = 1, nhydrothreadsglobal
        call MPI_RECV(tempStorage, 24, MPI_DOUBLE_PRECISION, i, &
             tag, localWorldCommunicator, status, ierr)
        if(tempStorage(1) < minR) then
           minR = tempStorage(1)
        end if
        if(tempStorage(3) > maxR) then
           maxR = tempStorage(3)
        end if
        mIo = mIo + tempStorage(4)
        mNeu = mNeu + tempStorage(5)
        movd = movd + tempStorage(6)
        KE = KE + tempStorage(7)
        mom = mom + tempStorage(8)
        meanR = meanR + tempStorage(2)
        nRadii = nRadii + int(tempStorage(9))
        do j = 1, 15
           rhoDist(j) = rhoDist(j) + tempStorage(j+9)
        end do
     end do
     if(nRadii /= 0) then
        meanR = meanR /dble(nRadii)
     else
        meanR = 0.d0
     end if
     write(65, '(9(e12.3, 3x))') grid%currentTime, minR, meanR, maxR, mIo, mNeu, Movd, KE, mom
     close(65)
     write(66, '(16(e12.3, 3x))') grid%currentTime, rhoDist(1), rhoDist(2), rhoDist(3), rhoDist(4), rhoDist(5), &
          rhoDist(6), rhoDist(7), rhoDist(8), rhoDist(9), rhoDist(10), rhoDist(11), rhoDist(12), &
          rhoDist(13), rhoDist(14), rhoDist(15)
     close(66)
  end if
end subroutine dumpWhalenNormanTest

subroutine dumpIfrontTest(grid)
  use mpi
  use inputs_mod, only : readgrid
  type(gridtype) :: grid
  integer :: ier, ierr, i, j
  logical, save :: firstTime=.true.
  real(double) :: tempStorage(24), rhoBins(15), rhoDist(15)
  integer :: stat(MPI_STATUS_SIZE)
  integer :: tag = 50
  real(double) :: minR, meanR, maxR, mIo, mNeu, Movd, KE, mom
  integer :: nRadii


  if(myRankGlobal == 0) then
     if(firstTime .and. .not. readGrid) then
        open(unit=25, file="planarIfront_testA_torus.txt", &
             status="replace", form="formatted",iostat=ier)
        firstTime = .false.
        write (25, '(3(a10, 3x))') "time", "m_ionized", "m_neutral"
     else
        open(unit=25, file="planarIfront_testA_torus.txt", &
             status="old", position="append", iostat=ier)
     end if
  end if
  
  minR = 1.d30
  maxR = 0.d0
  meanR = 0.d0
  mIo = 0.d0
  mNeu = 0.d0
  Movd = 0.d0
  KE = 0.d0
  mom = 0.d0
  nRadii = 0
  rhoDist = 0.d0
  if(myRankGlobal /= 0) then
     call getTestValues(grid%octreeRoot, minR, meanR, maxR, mIo, mNeu, &
          Movd, KE, mom, nRadii, rhoBins, rhoDist)
     
     do i = 1, nhydrothreadsglobal
        if(i == myRankGlobal) then
           tempStorage(1) = minR
           tempStorage(2) = meanR
           tempStorage(3) = maxR
           tempStorage(4) = mIo
           tempStorage(5) = mNeu
           tempStorage(6) = Movd
           tempStorage(7) = KE
           tempStorage(8) = mom
           tempStorage(9) = nRadii
           tempStorage(10) = rhoDist(1)
           tempStorage(11) = rhoDist(2)
           tempStorage(12) = rhoDist(3)
           tempStorage(13) = rhoDist(4)
           tempStorage(14) = rhoDist(5)
           tempStorage(15) = rhoDist(6)
           tempStorage(16) = rhoDist(7)
           tempStorage(17) = rhoDist(8)
           tempStorage(18) = rhoDist(9)
           tempStorage(19) = rhoDist(10)
           tempStorage(20) = rhoDist(11)
           tempStorage(21) = rhoDist(12)
           tempStorage(22) = rhoDist(13)
           tempStorage(23) = rhoDist(14)
           tempStorage(24) = rhoDist(15)
           call MPI_SEND(tempStorage, 24, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)        
        end if
     end do
  end if

  if(myRankGlobal == 0) then
     do i = 1, nhydrothreadsglobal
        call MPI_RECV(tempStorage, 24, MPI_DOUBLE_PRECISION, i, &
             tag, localWorldCommunicator, stat, ierr)
        if(tempStorage(1) < minR) then
           minR = tempStorage(1)
        end if
        if(tempStorage(3) > maxR) then
           maxR = tempStorage(3)
        end if
        mIo = mIo + tempStorage(4)
        mNeu = mNeu + tempStorage(5)
        movd = movd + tempStorage(6)
        KE = KE + tempStorage(7)
        mom = mom + tempStorage(8)
        meanR = meanR + tempStorage(2)
        nRadii = nRadii + int(tempStorage(9))
        do j = 1, 15
           rhoDist(j) = rhoDist(j) + tempStorage(j+9)
        end do
     end do
     if(nRadii /= 0) then
        meanR = meanR /dble(nRadii)
     else
        meanR = 0.d0
     end if
     write(25, '(3(e12.3, 3x))') grid%currentTime, mIo, mNeu
     close(25)

  end if
end subroutine dumpIfrontTest

subroutine setupDensityBins(rhoBins)
  real(double) :: rhoBins(15)
  integer :: i

  rhoBins(1) = 0
  rhoBins(2) = 10
  do i = 3, 15
     rhoBins(i) = rhoBins(i-1)*2.d0
  end do

end subroutine

subroutine updateBins(rhoBins, rhoDist, rho, dv)
  real(double) :: rho, rhoBins(15),  rhoDist(15), dv
  integer :: i
  logical :: found

  found = .false.  

  do i = 1, 14
     if (rho > rhoBins(i) .and. rho < rhoBins(i+1)) then
        rhoDist(i) = rhoDist(i) + rho*dV*mHydrogen
        found = .true.
     end if
  end do
  if(.not. found) then
     rhoDist(15) =  rhoDist(15) + rho*dV*mHydrogen
     found = .true.
  end if

end subroutine

recursive subroutine getTestValues(thisOctal, minR, meanR, maxR, mIo, &
     mNeu, Movd, KE, mom, nRadii, rhoBins, rhoDist)
  use mpi
  type(octal), pointer :: thisOCtal, child
  integer :: subcell
  real(double) :: minR, meanR, maxR, mIo, mNeu, Movd, KE, mom
  integer :: nRadii, i
  real(double) :: u2, eKinetic, dv, thisR, rho
  type(VECTOR) :: thisRVec, cen
  real(double) :: rhoBins(15), rhoDist(15)

  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call getTestValues(child, minR, meanR, maxR, mIo, mNeu, &
                   Movd, KE, mom, nRadii, rhoBins, rhoDist)
              exit
           end if
        end do
     else
        if (octalOnThread(thisOctal, subcell, myrankGlobal)) then

           dv = cellVolume(thisOctal, subcell)*1.d30

           if(thisOctal%ionFrac(subcell,2) < 0.9 .and. &
                thisOctal%ionFrac(subcell,2) > 0.1) then
              !found the ionization front             
              thisRVec = subcellCentre(thisOctal, subcell)

              if(thisOctal%threeD) then
                 cen = VECTOR(3.160d8, 3.160d8, 3.160d8)
              else if (thisOctal%twoD) then
                 cen = VECTOR(3.160d8, 0.d0, 3.160d8)
              end if              

              thisR = modulus(thisRVec - cen)

              if(thisR < minR) then
                 !update minimum I front radius
                 minR = thisR
              end if

              if(thisR > maxR) then
                 !update maximum I front radius
                 maxR = thisR
              end if
              
              meanR = meanR + thisR
              nRadii = nRadii + 1
           end if
           
           rho = thisOctal%rho(subcell)/mHydrogen
           call updateBins(rhoBins, rhoDist, rho, dv)
           
           if(thisOctal%ionFrac(subcell,2) > 0.9) then           
              !found ionized gas
              mIo = mIo + (thisOctal%rho(subcell)*dv)

           else if (thisOctal%ionFrac(subcell, 2) < 0.1) then
              !found neutral gas
              mNeu = mNeu + (thisOctal%rho(subcell)*dv)
           end if

           if(thisOctal%rho(subcell) > 2.338e-20) then
              !found an overdense cell
              Movd = Movd + (thisOctal%rho(subcell)*dv)
           end if
           

           u2 = (thisOctal%rhou(subcell)**2 + thisOctal%rhov(subcell)**2 + thisOctal%rhow(subcell)**2)/&
                thisOctal%rho(subcell)**2
           
           eKinetic = u2 / 2.d0

           KE = KE + (eKinetic*thisOctal%rho(subcell)*dv)

           mom = mom + (thisOctal%rho(subcell)*dv*sqrt(u2))
          
        end if
     end if
  end do
end subroutine getTestValues



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

  if(grid%octreeroot%twod) then
     startPoint = vector(4.4d9, 0.d0, 4.4d9)
     endPoint = vector(8.8d9, 0.d0, 4.4d9)
  else
     startPoint = vector(0.d0, 0.d0, 0.d0)
     endPoint = vector(4.4d9, 0.d0, 0.d0)
  end if
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

           
           if(grid%octreeroot%twod) then
              r = 4.4d9 + ((1.+7.d0*dble(i-1)/499.d0)*pctocm/1.e10)
              position = vector(r, 0.d0, 4.4d9)
           else
              r = (1.+7.d0*dble(i-1)/499.d0)*pctocm/1.e10
              position = vector(r, 0.d0, 0.d0)
           end if
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

           if(grid%octreeroot%twod) then
              write(21,'(f6.3,1p,6e12.3,0p)') (r*1.e10/pctocm)-4.4e19/pctocm,heating,cooling,oirate,oiirate,oiiirate,oivrate
              write(20,'(f6.3,f9.1,  14f8.3)') &
                   (r*1.e10/pctocm)-4.4e19/pctocm,t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
              write(22,*) (r*1.e10/pctocm)-4.4e19/pctocm,netot
           else
              write(21,'(f6.3,1p,6e12.3,0p)') (r*1.e10/pctocm),heating,cooling,oirate,oiirate,oiiirate,oivrate
              write(20,'(f6.3,f9.1,  14f8.3)') &
                   (r*1.e10/pctocm),t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
              write(22,*) r*1.e10/pctocm,netot
           end if
!           if(grid%octreeroot%twod) then
!              write(20,'(f6.3,f9.1,  14f8.3)') &
!                   (r*1.e10/pctocm)-4.4e19/pctocm,t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
!              write(22,*) r*1.e10/pctocm,netot
!           else
!              write(20,'(f6.3,f9.1,  14f8.3)') &
!                   (r*1.e10/pctocm),t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
!              write(22,*) r*1.e10/pctocm,netot
!           end if

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

        if(grid%octreeroot%oned) then
           octVec = VECTOR(r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta))
        else
           octVec = VECTOR(r, 0.d0, 4.4d9)
        end if

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
          if (octalOnThread(thisOctal, subcell, myrankGlobal)) then
             rVec = subcellCentre(thisOctal,subcell)
!             v = 3.d0*cellVolume(thisOctal, subcell)/fourpi
             v = cellVolume(thisOctal, subcell)
             
             call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell))
             
             rate =  pops(grid%ion(iion)%transition(iTrans)%j) * grid%ion(iion)%transition(itrans)%energy * &
                  grid%ion(iion)%transition(itrans)%a/ergtoev
             rate = rate * grid%ion(iion)%abundance * thisOctal%nh(subcell) * thisOctal%ionFrac(subcell, iion)
             luminosity = luminosity + rate * v * 1.d30
          end if
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
  
!!$subroutine getSahaMilneFreq(table,temperature, thisFreq)
!!$  type(SAHAMILNETABLE) :: table
!!$  real(double) :: temperature, thisfreq, r, t, fac
!!$  integer :: i, j
!!$
!!$  t = max(5000.d0, min(20000.d0, temperature))
!!$  call locate(table%temp, table%nTemp, t, i)
!!$  call randomNumberGenerator(getDouble=r)
!!$  call locate(table%Clyc(i,1:table%nfreq), table%nFreq, r, j)
!!$  fac = (r - table%Clyc(i,j))/(table%Clyc(i,j+1)-table%cLyc(i,j))
!!$  thisFreq = table%freq(j) + fac * (table%freq(j+1)-table%freq(j))
!!$end subroutine getSahaMilneFreq

!!$subroutine twoPhotonContinuum(thisFreq)
!!$
!!$! based on table ii of drake, victor, dalgarno, 1969, PhyRev Vol 180, pg 25
!!$
!!$  real(double) :: thisFreq
!!$  real :: y(21) = (/ 0., 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300, &
!!$       0.325, 0.350, 0.375, 0.400, 0.425, 0.450, 0.475, 0.500 /)
!!$  real :: hei(21) = (/ 0., 7.77e0, 2.52e1, 4.35e1, 5.99e1, 7.42e1, 8.64e1, 9.69e1, 1.06e2, 1.13e2, 1.20e2, 1.25e2, &
!!$       1.30e2, 1.34e2, 1.37e2, 1.40e2, 1.42e2, 1.43e2, 1.45e2, 1.45e2, 1.45e2 /)
!!$  real :: freq = 3.86e15, fac, r
!!$  real :: prob(21)
!!$  integer :: i
!!$  prob(1) = 0.
!!$  do i = 2, 21
!!$     prob(i) = prob(i-1) + (y(i)-y(i-1)) * hei(i)
!!$  enddo
!!$  prob(1:21) = prob(1:21)/prob(21)
!!$  thisFreq = 0.
!!$  do while((thisFreq*hcgs*ergtoev) < 13.6)
!!$     call randomNumberGenerator(getReal=r)
!!$     call locate(prob, 21, r, i)
!!$     fac = y(i) + ((r - prob(i))/(prob(i+1)-prob(i)))*(y(i+1)-y(i))
!!$     thisFreq = (1.-fac)*freq
!!$  enddo
!!$end subroutine twoPhotonContinuum

subroutine getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  integer :: subcell
  real(double) :: v, epsOverDeltaT
  real(double), intent(out) :: hHeating, heHeating, totalHeating, dustHeating

  dustHeating  = 0.d0
  heHeating = 0.d0
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

  logT = log10(max(table%temp(1),min(temp,table%temp(table%nTemp))))
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


!!$subroutine findForbiddenLine(lambda, grid, iIon, iTrans)
!!$  integer :: iIon, iTrans
!!$  type(GRIDTYPE) :: grid
!!$  real(double) :: lambda, lineLambda
!!$  logical :: foundLine
!!$  integer :: i, j
!!$
!!$  foundLine = .false.
!!$  do i = 3, grid%nIon
!!$     do j = 1, grid%ion(iIon)%nTransitions
!!$        lineLambda  =  cspeed * hcgs / (grid%ion(i)%transition(j)%energy / ergtoEV) / (angstromToCm)
!!$        if (abs((lineLambda -lambda)/lambda) < 0.01d0) then
!!$           foundLine = .true.
!!$           iIon = i
!!$           iTrans = j
!!$         endif
!!$     enddo
!!$  enddo
!!$  if (.not.foundLine) then
!!$     write(*,*) "Error finding forbidden line at: ",lambda
!!$     stop
!!$  endif
!!$end subroutine findForbiddenLine
!!$
!!$       
!!$subroutine addHeRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
!!$
!!$  integer :: nFreq
!!$  real(double) :: spectrum(:), freq(:)
!!$  type(OCTAL) :: thisOctal
!!$  integer :: subcell
!!$  type(GRIDTYPE) :: grid
!!$  integer :: i, j, k
!!$  real :: fac, t, aj ,bj, cj
!!$  real(double) :: lineFreq !, lambda
!!$  real(double) :: emissivity
!!$!  real :: heII4686
!!$!  integer :: ilow, iup
!!$  integer,parameter :: nHeIILyman = 4
!!$!  real(double) :: heIILyman(4)
!!$!  real(double) :: freqheIILyman(4) = (/ 3.839530, 3.749542, 3.555121, 2.99963 /)
!!$
!!$  ! HeI lines 
!!$
!!$  call locate(heIrecombinationNe, 3, real(log10(thisOctal%ne(subcell))), i)
!!$  fac = real ( &
!!$       (log10(thisOctal%ne(subcell)) - heIrecombinationNe(i))/(heIrecombinationNe(i+1)-heIrecombinationNe(i)) )
!!$
!!$  do j = 1, 32
!!$     aj = heIrecombinationFit(j,i,1) + fac*(heIrecombinationfit(j,i+1,1)-heIrecombinationfit(j,i,1))
!!$     bj = heIrecombinationFit(j,i,2) + fac*(heIrecombinationfit(j,i+1,2)-heIrecombinationfit(j,i,2))
!!$     cj = heIrecombinationFit(j,i,3) + fac*(heIrecombinationfit(j,i+1,3)-heIrecombinationfit(j,i,3))
!!$     t = thisOctal%temperature(subcell)/1.e4
!!$     emissivity = aj * (t**bj) * exp(cj / t) ! Benjamin et al. 1999 ApJ 514 307
!!$     emissivity = emissivity * thisOctal%ne(subcell) * thisOctal%nh(subcell) * &
!!$          thisOctal%ionFrac(subcell, 3) * grid%ion(3)%abundance
!!$
!!$     lineFreq = cspeed / (heiRecombinationLambda(j)*1.e-8)
!!$     call locate(freq, nFreq, lineFreq, k)
!!$     k = k + 1
!!$     spectrum(k) = spectrum(k) + emissivity
!!$  enddo
!!$
!!$end subroutine addHeRecombinationLines

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

recursive subroutine packvalues(thisOctal,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec, &
     kappaTimesFlux, uvVec)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: nCrossings(:)
  real(double) :: photoIonCoeff(:,:)
  real(double) :: hHeating(:)
  real(double) :: heHeating(:)
  real(double) :: distanceGrid(:)
  type(VECTOR) :: radMomVec(:)
  type(VECTOR) :: kappaTimesFlux(:)
  type(VECTOR) :: uvVec(:)
  
  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call packvalues(child, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec, &
                   kappaTimesFlux, uvVec)
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
        kappaTimesFlux(nIndex) = thisOctal%kappaTimesFlux(subcell)
        uvVec(nIndex) = thisOctal%uvVector(subcell)
     endif
  enddo
end subroutine packvalues

recursive subroutine unpackvalues(thisOctal,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec, &
     kappaTimesFlux, uvVec)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: nCrossings(:)
  real(double) :: photoIonCoeff(:,:)
  real(double) :: hHeating(:)
  real(double) :: heHeating(:)
  real(double) :: distanceGrid(:)
  type(VECTOR) :: radMomVec(:)
  type(VECTOR) :: kappaTimesFlux(:)
  type(VECTOR) :: uvvec(:)

  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unpackvalues(child, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec, &
                     kappaTimesFlux, uvvec)
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
          thisOctal%kappaTimesFlux(subcell) = kappaTimesFlux(nIndex)
          thisOctal%uvvector(subcell) = uvvec(nIndex)
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

  recursive subroutine  zeroCornerEmissivity(grid, thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i

  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroCornerEmissivity(grid, child)
                exit
             end if
          end do
       else
          if(octalonthread(thisOctal, subcell, myrankglobal)) then
             if (thisOctal%edgecell(subcell)) then
                thisOctal%etaCont(subcell) = 0.d0
             endif
          end if
       endif
    enddo

  end subroutine zeroCornerEmissivity


!!$  function getRandomWavelengthPhotoion(grid, thisOctal, subcell, lamArray, nLambda) result(thisLambda)
!!$
!!$    type(GRIDTYPE) :: grid
!!$    type(OCTAL), pointer :: thisOctal
!!$    integer :: subcell
!!$    real :: thisLambda
!!$    integer :: nFreq
!!$    real(double), allocatable :: freq(:), spectrum(:), tspec(:),lamspec(:), dfreq(:)
!!$    real(double) :: nuStart, nuEnd, r, fac
!!$    integer :: nLambda
!!$    real :: lamArray(:)
!!$    integer :: i
!!$
!!$    nFreq = 1000
!!$
!!$    allocate(freq(1:nFreq), spectrum(1:nFreq), lamSpec(1:nFreq), dFreq(1:nFreq))
!!$    nuStart = cSpeed / (1000.d4 * 1.d-8)
!!$    nuEnd = cSpeed / (10.d0 * 1.d-8)
!!$
!!$    do i = 1, nFreq
!!$       freq(i) = log10(nuStart) + dble(i-1)/dble(nFreq-1) * (log10(nuEnd)-log10(nuStart))
!!$       freq(i) = 10.d0**freq(i)
!!$    enddo
!!$    do i = 2, nFreq-1
!!$       dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
!!$    enddo
!!$    dfreq(1) = (freq(2)-freq(1))
!!$    dfreq(nfreq) = (freq(nfreq)-freq(nfreq-1))
!!$
!!$
!!$    spectrum = 1.d-30
!!$
!!$    call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
!!$    call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
!!$    call addHydrogenRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
!!$    !                        call addHeRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
!!$    call addForbiddenLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
!!$    call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
!!$    
!!$
!!$    do i = 1, nFreq
!!$       lamSpec(i) = cspeed/freq(i)
!!$       spectrum(i) = spectrum(i) * cspeed/(lamspec(i)**2)
!!$    enddo
!!$    lamSpec = lamSpec * 1.d8
!!$
!!$    allocate(tSpec(1:nFreq))
!!$  
!!$    tSpec(1:nFreq) = spectrum(1:nFreq)
!!$
!!$    do i = 2, nFreq
!!$       tSpec(i) = tSpec(i) + tSpec(i-1)
!!$    enddo
!!$    tSpec(1:nFreq) = tSpec(1:nFreq) - tSpec(1)
!!$    if (tSpec(nFreq) > 0.d0) then
!!$       tSpec(1:nFreq) = tSpec(1:nFreq) / tSpec(nFreq)
!!$       call randomNumberGenerator(getDouble=r)
!!$       call locate(tSpec, nFreq, r, i)
!!$       fac = (r - tSpec(i)) / (tSpec(i+1)-tSpec(i))
!!$       thisLambda = real(lamspec(i) + fac * (lamspec(i+1)-lamspec(i)))
!!$    else
!!$       thisLambda = 1000.e4
!!$  endif
!!$    
!!$  deallocate(freq, spectrum, lamSpec, dfreq)
!!$
!!$  end function getRandomWavelengthPhotoion


  recursive subroutine calculateEnergyFromTemperature(thisOctal)
    use inputs_mod, only : hOnly, simplemu
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
          if(hOnly .and. simplemu) mu = returnMuSimple(thisOctal, subcell)

             
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

  recursive subroutine calculateKappaTimesFlux(thisOctal, epsOverDeltaT)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: v, epsOverDeltaT

    do subcell = 1, thisOctal%maxChildren
       if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateKappaTimesFlux(child, epsOverDeltaT)
                exit
             end if
          end do
       else
          V = cellVolume(thisOctal, subcell)*1.d30
          thisOctal%kappaTimesFlux(subcell) = epsOverDeltaT * thisOctal%kappaTimesFlux(subcell) / V

          thisOctal%radiationMomentum(subcell) = thisOctal%radiationMomentum(subcell) / V


!          write(*,*) " k times f ",thisOctal%kappaTimesFlux(subcell)
       endif
    enddo
  end subroutine calculateKappaTimesFlux


  recursive subroutine calculateUVfluxVec(thisOctal, epsOverDeltaT)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: v, epsOverDeltaT

    do subcell = 1, thisOctal%maxChildren
       if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call calculateUVfluxVec(child, epsOverDeltaT)
                exit
             end if
          end do
       else
          V = cellVolume(thisOctal, subcell)*1.d30
          thisOctal%UVvector(subcell) = epsOverDeltaT * thisOctal%Uvvector(subcell) / V
!          write(*,*) " k times f ",thisOctal%kappaTimesFlux(subcell)
       endif
    enddo
  end subroutine calculateUVfluxVec


  subroutine createImageSplitGrid(grid, nSource, source, imageNum)
    use inputs_mod, only: nPhotImage, gridDistance
    use hydrodynamics_mod, only: setupEdges
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
    integer :: np(100)
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

    if(myRankGlobal /= 0.d0) then
       call setupEdges(grid%octreeRoot, grid)
       
       call zeroCornerEmissivity(grid, grid%octreeRoot)
    end if
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

    powerPerPhoton = (lCore + totalEmission) / (dble(nPhotImage)*1.d20)
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
       mainloop: do iPhoton = 1, nPhotImage

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
          if (nDone == nPhotImage) photonsStillProcessing = .false.
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
     endif


     tArray = 0.d0
     call MPI_ALLREDUCE(totalEmissionArray, tArray, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, &
          zeroPlusAmrCommunicator, ierr)
     totalEmissionArray = tArray

     tArray = 0.d0
     call MPI_ALLREDUCE(totalProbArray, tArray, nHydroThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, &
          zeroPlusAmrCommunicator, ierr)
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

!!$  subroutine  checkSetsHaveSameNumberOfOctals(grid, iThread)
!!$    use mpi
!!$    type(GRIDTYPE) :: grid
!!$    integer :: iThread
!!$    integer, allocatable :: nVoxels(:), temp(:)
!!$    integer :: ierr, i
!!$    logical :: same
!!$    allocate(nVoxels(1:NHydroSetsGlobal), temp(1:NHydroSetsGlobal))
!!$    nVoxels = 0
!!$    temp = 0
!!$    call countVoxelsOnThread(grid%octreeRoot, nVoxels(myHydroSetGlobal+1))
!!$    call mpi_allreduce(nVoxels, temp, nHydroSetsGlobal, &
!!$         MPI_INTEGER, MPI_SUM, amrParallelCommunicator(myrankGlobal), &
!!$         ierr)
!!$    nVoxels = temp
!!$    same = .true.
!!$    do i = 1, size(nVoxels)
!!$       if (any(nVoxels(i) /= nVoxels)) same = .false.
!!$    enddo
!!$    if (.not.same) then
!!$       write(*,*) "ERROR!! Octal numbers differ for thread ",ithread, " sets ",nVoxels
!!$    endif
!!$    deallocate(nVoxels)
!!$  end subroutine checkSetsHaveSameNumberOfOctals
    
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
    type(VECTOR), allocatable :: kappaTimesFlux(:)
    type(VECTOR), allocatable :: uvvec(:)
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
    allocate(kappaTimesFlux(1:nVoxels))
    allocate(uvvec(1:nVoxels))

    nIndex = 0
    call packValues(grid%octreeRoot,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec, &
         kappaTimesFlux, uvvec)

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

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(kappaTimesFlux(1:nVoxels)%x,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    kappaTimesFlux(1:nVoxels)%x = tempDoubleArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(kappaTimesFlux(1:nVoxels)%y,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    kappaTimesFlux(1:nVoxels)%y = tempDoubleArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(kappaTimesFlux(1:nVoxels)%z,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    kappaTimesFlux(1:nVoxels)%z = tempDoubleArray 


    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(uvvec(1:nVoxels)%x,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    uvvec(1:nVoxels)%x = tempDoubleArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(uvvec(1:nVoxels)%y,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    uvvec(1:nVoxels)%y = tempDoubleArray 

    tempDoubleArray = 0.0
    call MPI_ALLREDUCE(uvvec(1:nVoxels)%z,tempDoubleArray,nVoxels,MPI_DOUBLE_PRECISION,&
         MPI_SUM, amrParComm, ierr)
    uvvec(1:nVoxels)%z = tempDoubleArray 

    deallocate(tempRealArray, tempDoubleArray)
     
    call MPI_BARRIER(amrParComm, ierr) 
    
    nIndex = 0
    call unpackValues(grid%octreeRoot, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid, radMomVec, &
         kappaTimesFlux, uvvec)
    deallocate(nCrossings, photoIonCoeff, hHeating, heHeating, distanceGrid, radMomVec, kappaTimesFlux, uvvec)

  end subroutine updateGridMPIphoto


  subroutine modifiedRandomWalk(grid, thisOctal, subcell, rVec, uHat, &
       freq, dfreq, nfreq, lamArray, nlambda, thisFreq, photonPacketWeight)
!    use inputs_mod, only : smallestCellSize
    type(GRIDTYPE) :: grid
    real(double) :: spectrum(2000)
    real(double) :: freq(:), dfreq(:), thisFreq, photonPacketWeight
    real :: lamArray(:)
    integer :: nFreq, nlambda
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(VECTOR) :: rVec, uHat, uHatDash, rHat, zHat
    real(double) :: r0, zeta, kappaRoss, diffCoeff, mrwDist
    real(double) :: kappaP
    logical, save :: firstTime = .true.
    integer, parameter :: ny = 100
    real(double) :: y(ny), prob(ny), thisY
    integer :: i, n, iLoop
    real(double), parameter :: gamma = 3.d0

    if (firstTime) then
       do i = 1, ny
          y(i) = dble(i-1)/dble(ny-1)
          prob(i) = 0.d0
          do n = 1, 1000
             prob(i) = prob(i) + 2.d0 * (-1.d0)**(n+1) * y(i)**(n**2)
          enddo
       enddo
       prob(ny) = 1.d0
       firstTime = .false.
!       do i = 1, ny
!          write(*,*) i, y(i), prob(i)
!       enddo
    endif
       
       

    call returnKappa(grid, thisOctal, subcell,  rosselandKappa = kappaRoss)
    kappap = thisOctal%kappap(subcell)
    diffCoeff = 1.d0 / (3.d0*thisOctal%rho(subcell)*kappaRoss)
    iLoop = 0
    call distanceToNearestWall(rVec, r0, thisOctal, subcell)
    do
       iLoop = iLoop + 1
       if (iLoop > 100000) then
          write(*,*) "modified random walk loop exitted early"
          exit
       endif
       call randomNumberGenerator(getDouble=zeta)
       call locate(prob, ny, zeta, i)
       thisY = y(i) + (y(i+1)-y(i))*(zeta-prob(i))/(prob(i+1)-prob(i))
       mrwDist = -log(thisy) * (r0/pi)**2 * (1.d0 / diffCoeff)
       thisOctal%distanceGrid(subcell) = thisOctal%distanceGrid(subcell) + mrwDist * kappap * photonPacketWeight
       uHatDash = uHat
       if (thisOctal%twoD .and. .not. cart2d) then
          rHat = VECTOR(rVec%x, rVec%y, 0.d0)
          call normalize(rHat)
          zHat = VECTOR(0.d0, 0.d0, 1.d0)
          uHatDash = VECTOR(rHat.dot.uHat, 0.d0, zHat.dot.uHat)
       endif
       if (thisOctal%oneD) then
          rHat = VECTOR(rVec%x, rVec%y, rVec%z)
          call normalize(rHat)
          uHatDash = VECTOR(rHat.dot.uHat, 0.d0, 0.d0)
       endif

!       thisOctal%kappaTimesFlux(subcell) = thisOctal%kappaTimesFlux(subcell) + (mrwDist * kappap * photonPacketWeight)*uHatDash
       thisOctal%UVvector(subcell) = thisOctal%UVvector(subcell) + (mrwDist * photonPacketWeight)*uHatDash
       rVec = rVec + uHat * r0
       uHat = randomUnitVector()
       call distanceToNearestWall(rVec, r0, thisOctal, subcell)
       if (r0 < gamma/(thisOctal%rho(subcell)*kappaRoss*1.d10)) exit
    enddo
    spectrum = 1.d-50
    call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
    thisFreq =  getPhotonFreq(nfreq, freq, spectrum)
!    write(*,*) "exited modified random walk after ",i
  end subroutine modifiedRandomWalk
       
  subroutine setupKappaPArrays(grid)
    use inputs_mod, only : nDustType
    type(GRIDTYPE) :: grid
    integer :: i, j, itemp, iDust
    real(double) :: temperature, freq, dfreq, kappap, norm
    do iTemp = 1, 2000
       do iDust = 1, nDustType
          temperature = dble(itemp)
          kappaP = 0.d0
          norm = 0.d0
          do i = 2, grid%nLambda
             freq = cSpeed / (grid%lamArray(i)*1.d-8)
             dfreq = cSpeed / (grid%lamArray(i-1)*1.d-8) - cSpeed / (grid%lamArray(i)*1.d-8)
             do j = 1, nDustType
                kappaP = kappaP + dble(grid%oneKappaAbs(iDust,i)) * &
                     bnu(freq,temperature)  * dfreq
             enddo
             norm = norm + bnu(freq,temperature)  * dfreq
          enddo
          if (norm /= 0.d0) then
             kappaP = ((kappaP / norm) /1.d10)
          else
             kappaP = tiny(kappap)
          endif
          kappapArray(iDust, iTemp) = kappaP
       enddo
    enddo
  end subroutine setupKappaPArrays

  function returnKappaP(thisOctal, subcell, temperature) result(kappap)
    use inputs_mod, only : nDustType
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: temperature
    real(double) :: kappaP, frac
    integer :: itemp, i, itemp1

    itemp = min(int(temperature),2000)
    itemp1 = min(itemp+1,2000)
    frac = min(temperature,2000.d0) - dble(itemp)

    kappaP = 0.d0
    do i = 1, nDustType
       kappaP = kappaP + thisOctal%dustTypeFraction(subcell,i) * &
            ((1.d0-frac) * kappaPArray(i, itemp) + frac*kappaPArray(i, itemp1))
    enddo
    kappaP = kappaP * thisOctal%rho(subcell)
  end function returnKappaP

  recursive subroutine setKappaP(thisOctal, grid)
    use gas_opacity_mod, only: returnGasKappaValue
    use atom_mod, only : bnu
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, k
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do k = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(k) == subcell) then
                child => thisOctal%child(k)
                call setKappaP(child, grid)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%kappap)) then
             allocate(thisOctal%kappap(1:ThisOctal%maxChildren))
          endif
          thisOctal%kappaP(subcell) = returnKappaP(thisOctal, subcell, dble(thisOctal%temperature(subcell)))
          
       endif
    enddo
  end subroutine setKappaP




  recursive subroutine radpressureTimeStep(thisoctal, dt)
    use mpi
    use inputs_mod, only : gridDistanceScale
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: dt, dx, acc, dv
    



    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call radpressureTimestep(child, dt)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankGlobal)) cycle
          if (.not.thisoctal%ghostcell(subcell)) then

             dv = cellVolume(thisOctal, subcell)*1.d30
             dx = thisoctal%subcellsize * griddistancescale
             if (.not.associated(thisOctal%kappaTimesFlux)) then
                allocate(thisOctal%kappaTimesFlux(1:thisOctal%maxChildren))
                thisOctal%kappaTimesFlux = VECTOR(0.d0,0.d0,0.d0)
             endif
             if (.not.associated(thisOctal%radiationMomentum)) then
                allocate(thisOctal%radiationMomentum(1:thisOctal%maxChildren))
                thisOctal%kappaTimesFlux = VECTOR(0.d0,0.d0,0.d0)
             endif
!             acc =  (1.d0/thisOctal%rho(subcell)) * modulus(thisOctal%kappaTimesFlux(subcell)/cSpeed)
             acc =  (1.d0/thisOctal%rho(subcell)) * modulus(thisOctal%radiationMomentum(subcell))
             acc = max(1.d-30, acc)
             dt = min(dt, sqrt(2.d0*dx/acc))
          endif
       endif
    enddo
  end subroutine radpressureTimeStep

subroutine putStarsInGridAccordingToDensity(grid, nSource, source)
  use mpi
  use source_mod, only : sourcetype
  use math_mod, only : computeProbDist
  use inputs_mod, only : rhoFloor
  type(GRIDTYPE) :: grid
  integer :: nSource
  type(SOURCETYPE) :: source(:)
  real(double) :: mass, temp(3), massInCell, starMass
  type(OCTAL), pointer :: sourceOctal
  integer :: sourceSubcell
  integer :: i, j
  logical :: looping
  real(double), allocatable :: threadProbArray(:)
  real(double) :: r, dv, rhoLeft
  integer :: iThread, iter
  integer, parameter :: tag = 66
  integer :: ierr, signal
  integer :: status(MPI_STATUS_SIZE)
  type(VECTOR) :: rVec

  allocate(threadProbArray(1:nHydroThreadsGlobal))
  ! calculate the prob distribution based on density
  call putDensityInEtaCont(grid%octreeRoot)
  call computeProbDistAMRMpi(grid, mass, threadProbArray)
 
!  if (writeoutput) write(*,*) "Total mass for put stars ",mass*1.d30/msol

  ! zeroth thread 
  if (myrankGlobal == 0) then
     do i = 1, nSource
        looping = .true.
        do while (looping)
           call randomNumberGenerator(getDouble=r)
           if (r < threadProbArray(1)) then
              iThread = 1
           else
              call locate(threadProbArray, SIZE(threadProbArray), r, iThread)
              iThread = iThread + 1
           endif
           signal = -99
           call mpi_send(signal, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
           call mpi_send(source(i)%mass, 1, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
           ! receives x,y,z position of source from other threads
           call MPI_RECV(temp, 3, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, status, ierr)
           looping = .false.
           source(i)%position%x = temp(1)
           ! if threads report back that not enough mass is present in the cell, keep looping
           if (temp(1) < -1.d20) then
              looping = .true.
              cycle
           endif
           source(i)%position%y = temp(2)
           source(i)%position%z = temp(3)
           looping = .false.
           if (i > 1) then
              do  j = 1, i-1
                 ! don't put source within 0.1 of accretion radius of another source
!                if (modulus(source(j)%position - source(i)%position) < 0.1d0*source(j)%accretionRadius/1.d10) then
                 if (modulus(source(j)%position - source(i)%position) < source(j)%accretionRadius/1.d10) then
                    looping = .true.
                 endif
              enddo
           endif
        enddo
     enddo
     do iThread = 1, nHydroThreadsGlobal
        signal = 0
        call mpi_send(signal, 1, MPI_INTEGER, iThread, tag, localWorldCommunicator, ierr)
     enddo
  ! other threads
  else
     do 
        call MPI_RECV(signal, 1, MPI_INTEGER, 0, tag, localWorldCommunicator, status, ierr)
        ! checks if enough mass is present in cell to put star there
        if (signal == -99) then
           call MPI_RECV(starMass, 1, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
           massInCell = 0.d0
           iter = 0
           do while(massInCell < starMass)
              iter = iter + 1
              call randomNumberGenerator(getDouble=r)
              sourceOctal => grid%octreeRoot
              call locateContProbAMR(r, sourceOctal,sourcesubcell)
              massInCell = sourceOctal%rho(sourcesubcell) * cellVolume(sourceOctal, sourceSubcell) * 1.d30
!              write(*,*) iter, myRankGlobal, " mass in cell ",massinCell/msol, "star mass ",starMass/msol
              if (iter > 1000) exit
           enddo
           write(*,'(i4,i3,a15,f12.5,a12,f12.5)') iter, myRankGlobal, " mass in cell: ",massinCell/msol, "star mass: ",starMass/msol
           rVec = randomPositionInCell(sourceOctal, sourceSubcell)
           temp(1) = rVec%x
           ! assumes the check is unsuccessful if iterations exceed 1000
! (put the star in anyway, even if iterations exceed 1000) 
!!!           if (iter > 1000) temp(1) = -1.d21
           temp(2) = rVec%y
           temp(3) = rVec%z
           call mpi_send(temp, 3, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, ierr)
        else
           exit
        endif
     enddo
  endif
  ! send the source positions to other threads
  if (myrankGlobal == 0) then
     do ithread = 1, nHydroThreadsGlobal
        call mpi_send(source(1:nSource)%position%x, nSource, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
        call mpi_send(source(1:nSource)%position%y, nSource, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
        call mpi_send(source(1:nSource)%position%z, nSource, MPI_DOUBLE_PRECISION, iThread, tag, localWorldCommunicator, ierr)
     enddo
  else
     call MPI_RECV(source(1:nSource)%position%x, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
     call MPI_RECV(source(1:nSource)%position%y, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
     call MPI_RECV(source(1:nSource)%position%z, nSource, MPI_DOUBLE_PRECISION, 0, tag, localWorldCommunicator, status, ierr)
     
  endif

  if (myrankGlobal /= 0) then
     sourceOctal => grid%octreeRoot
     do i = 1, nSource
        call findSubcellLocal(source(i)%position, sourceOctal, sourcesubcell)
        if (octalOnThread(sourceOctal, sourceSubcell, myRankGlobal)) then
           dv = cellVolume(sourceOctal, sourceSubcell) * 1.d30
           
           source(i)%velocity%x = sourceOctal%rhou(sourcesubcell)/sourceoctal%rho(sourcesubcell)
           source(i)%velocity%y = sourceOctal%rhov(sourcesubcell)/sourceoctal%rho(sourcesubcell)
           source(i)%velocity%z = sourceOctal%rhow(sourcesubcell)/sourceoctal%rho(sourcesubcell)
           ! remove source mass from the cell 
           rhoLeft = max(rhoFloor, sourceOctal%rho(sourcesubcell)  - source(i)%mass/dv)
           sourceOctal%rho(sourcesubcell) = rhoLeft
           sourceOctal%rhou(sourcesubcell) = rhoLeft * source(i)%velocity%x
           sourceOctal%rhov(sourcesubcell) = rhoLeft * source(i)%velocity%y
           sourceOctal%rhow(sourcesubcell) = rhoLeft * source(i)%velocity%z
        endif
     enddo
  end if


end subroutine putStarsInGridAccordingToDensity

  recursive subroutine putDensityInEtaCont(thisOctal) 
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call putDensityInEtaCont(child)
                exit
             end if
          end do
       else 
          if (.not.associated(thisOctal%etaCont)) allocate(thisOctal%etaCont(1:thisOctal%maxChildren))
          if (.not.associated(thisOctal%etaLine)) allocate(thisOctal%etaLine(1:thisOctal%maxChildren))
          if (.not.associated(thisOctal%biasCont3D)) allocate(thisOctal%biasCont3D(1:thisOctal%maxChildren))
          thisOctal%biasCont3d = 1.d0
          thisOctal%etaCont(1:thisOctal%maxChildren) = thisOctal%rho(1:thisOctal%maxChildren)**1.5d0
       endif
    enddo
  end subroutine putDensityInEtaCont


#endif


end module photoionAMR_mod


#endif
