!Photoionization module - started on October 4th 2005 by th


module photoionAMR_mod

#ifdef MPI
use constants_mod
use messages_mod

use hydrodynamics_mod
use parallel_mod
use lucy_mod, only : calcContinuumEmissivityLucyMono
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
public :: photoIonizationloopAMR, radiationHydro, createImagesplitgrid, ionizeGrid, &
     neutralGrid, resizePhotoionCoeff, resetNH

type SAHAMILNETABLE
   integer :: nFreq 
   integer :: nTemp 
   real(double), pointer :: temp(:)
   real(double), pointer :: freq(:)
   real(double), pointer :: Clyc(:, :)
   real(double), pointer :: emissivity(:)
end type SAHAMILNETABLE

type RECOMBTABLE
   integer :: nrho
   integer :: nTemp 
   real(double), pointer :: rho(:)
   real(double),  pointer :: temp(:)
   real(double), pointer :: emissivity(:, :)
end type RECOMBTABLE


type GAMMATABLE
   integer :: nGamma
   integer :: nFreq
   integer :: nTemp
   real(double), pointer :: freq(:)
   real(double), pointer :: temp(:)
   real(double), pointer :: gamma(:,:)
end type GAMMATABLE



  real :: heIRecombinationFit(32, 3, 3)
  real :: heIRecombinationLambda(32)
  real :: heIRecombinationNe(3)
  real :: heIIrecombinationLines(3:30, 2:16)
  type(GAMMATABLE) :: gammaTableArray(3) ! H, HeI, HeII
  integer :: nEscapedGlobal=0

contains


  subroutine radiationHydro(grid, source, nSource, nLambda, lamArray)
    use ion_mod, only : returnMu
    use input_variables, only : iDump, doselfgrav !, hOnly
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: nLambda
    real :: lamArray(:)
    character(len=80) :: mpiFilename, datFilename
    real(double) :: dt, tc(65), temptc(65),cfl, gamma, mu
    integer :: iUnrefine
    integer :: myRank, ierr
    real(double) :: tDump, nextDumpTime, tEnd
    type(VECTOR) :: direction, viewVec
    logical :: gridConverged
    integer :: thread1(200), thread2(200), nBound(1000), nPairs
    integer :: group(1000), nGroup
!    logical :: globalConverged(64), tConverged(64)
    integer :: nHydroThreads
    logical :: dumpThisTime
    real(double) :: deltaTforDump, timeOfNextDump, loopLimitTime
    integer :: iRefine, nUnrefine
    logical :: startFromNeutral


    nHydroThreads = nThreadsGlobal-1
    dumpThisTime = .false.



    direction = VECTOR(1.d0, 0.d0, 0.d0)
    gamma = 7.d0 / 5.d0
    cfl = 0.3d0
    
    mu = 1.d0

    viewVec = VECTOR(-1.d0,0.d0,0.d0)
    viewVec = rotateZ(viewVec, 40.d0*degtorad)
    viewVec = rotateY(viewVec, 30.d0*degtorad)

    if (writeoutput) write(*,*) "Ionizing photons per second ",ionizingFlux(source(1))

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    grid%currentTime = 0.d0
    grid%iDump = 0
    tDump = 0.005d0
    deltaTforDump = 3.14d10 !1kyr
    nextDumpTime = deltaTforDump
    if (grid%geometry == "hii_test") deltaTforDump = 1.d8
    if(grid%geometry == "bonnor") deltaTforDump = 1.57d11 !5kyr

    iunrefine = 0
    startFromNeutral = .false.
!    if (grid%geometry == "bonnor") startFromNeutral = .true.


    if (readlucy) then
       write(mpiFilename,'(a, i4.4, a)') "dump_", iDump, ".grid"
       write(*,*) "Reading from: ",trim(mpiFilename)
       call readAmrGrid(mpiFilename,.false.,grid)
    endif

    timeofNextDump = 0.d0
    !timeofNextDump = deltaTforDump
    

    if (myRank == 1) write(*,*) "CFL set to ", cfl


    if (myrank /= 0) then
       call returnBoundaryPairs(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       call writeInfo("Calling exchange across boundary", TRIVIAL)
       if (myRank == 1) call tune(6,"Exchange across boundary")
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       if (myRank == 1) call tune(6,"Exchange across boundary")
       call writeInfo("Done", TRIVIAL)

       direction = VECTOR(1.d0, 0.d0, 0.d0)
       call calculateRhoU(grid%octreeRoot, direction)
       direction = VECTOR(0.d0, 1.d0, 0.d0)
       call calculateRhoV(grid%octreeRoot, direction)
       direction = VECTOR(0.d0, 0.d0, 1.d0)
       call calculateRhoW(grid%octreeRoot, direction)

       call calculateEnergyFromTemperature(grid%octreeRoot)
       
       call calculateRhoE(grid%octreeRoot, direction)



       call writeInfo("Refining individual subgrids", TRIVIAL)
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
          call evenUpGridMPI(grid, .false., .true.)
       endif


       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       call evenUpGridMPI(grid,.false.,.true.)      
       call refineGridGeneric(grid, 1.d-2)
       call writeInfo("Evening up grid", TRIVIAL)    
       call evenUpGridMPI(grid, .false.,.true.)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


    endif

    call neutralGrid(grid%octreeRoot)

    looplimitTime = deltaTForDump
    !looplimitTime = 0.1375d10
    do irefine = 1, 1
       if(grid%geometry == "hii_test") then
          loopLimitTime = deltaTForDump
       else 
          loopLimitTime = 2.d11       
       end if
       if (irefine == 1) then
          call writeInfo("Calling photoionization loop",TRIVIAL)
          call setupNeighbourPointers(grid, grid%octreeRoot)
          call photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, 60, loopLimitTime, looplimittime, .True.,&
               .false.)
             call photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, 20, loopLimitTime, looplimittime, .True.,&
                  .true.)

          call writeInfo("Done",TRIVIAL)
       else
          call writeInfo("Calling photoionization loop",TRIVIAL)
          call setupNeighbourPointers(grid, grid%octreeRoot)
          call photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, 20, loopLimitTime, looplimittime, .true., .true.)
          call writeInfo("Done",TRIVIAL)
       endif


!Track the evolution of the ionization front with time
    !      write(datFilename, '(a, i4.4, a)') "Ifront.dat"
    !      call dumpStromgrenRadius(grid, datFileName, VECTOR(-1.5d9,  -1.5d9, 1.5d9), &
     !          VECTOR(1.5d9, 1.5d9, 1.5d9), 1000)



       call writeInfo("Dumping post-photoionization data", TRIVIAL)
       call writeVtkFile(grid, "start.vtk", &
         valueTypeString=(/"rho        ","HI         " ,"temperature" /))


       !       call testIonFront(grid%octreeRoot, grid%currentTime)

       if (myrank /= 0) then
          call calculateEnergyFromTemperature(grid%octreeRoot)
          call calculateRhoE(grid%octreeRoot, direction)
       endif


       if (myrank/=0) then

          call writeInfo("Refining individual subgrids", TRIVIAL)
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
             call evenUpGridMPI(grid, .false., .true.)
          endif



       call evenUpGridMPI(grid,.false.,.true.)      
       call refineGridGeneric(grid, 1.d-2)
       call writeInfo("Evening up grid", TRIVIAL)    
       call evenUpGridMPI(grid, .false.,.true.)
       call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)


    endif
    enddo


    if (myrank /= 0) then

       call calculateEnergyFromTemperature(grid%octreeRoot)

       call calculateRhoE(grid%octreeRoot, direction)

       !       directi0on = VECTOR(1.d0, 0.d0, 0.d0)
       !       call calculateRhoU(grid%octreeRoot, direction)
       !       direction = VECTOR(0.d0, 1.d0, 0.d0)
       !       call calculateRhoV(grid%octreeRoot, direction)
       !       direction = VECTOR(0.d0, 0.d0, 1.d0)
       !       call calculateRhoW(grid%octreeRoot, direction)

    endif

    tEnd = 200.d0*3.14d10 !200kyr 

    do while(grid%currentTime < tEnd)
       tc = 0.d0
       tc(myrank+1) = 1.d30
       call computeCourantTime(grid, grid%octreeRoot, tc(myRank+1))
       call MPI_ALLREDUCE(tc, tempTc, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierr)
       tc = tempTc
       dt = MINVAL(tc(2:nThreadsGlobal)) * cfl

       

       !Dump at every step for now
       !if(grid%geometry == "hii_test") then
       !   dt = 1.d11
	!  !dumpThisTime = .true.
      ! e!nd if       

       if (myrank == 1) then
          write(*,*) tc(1:9)
          write(*,*) "courantTime", dt
       endif
       dumpThisTime = .false.

       if ((grid%currentTime + dt) >= timeOfNextDump) then
          dt =  timeofNextDump - grid%currentTime
          dumpThisTime = .true.
       endif


       if (myrank == 1) write(*,*) "dump ",dumpThisTime, " current ", &
            grid%currentTime, " deltaTfordump ",deltaTforDump, " dt ", dt

       if (myrank == 1) write(*,*) "Time step", dt
              
       


       if (myRank == 1) call tune(6,"Hydrodynamics step")
       call writeInfo("calling hydro step",TRIVIAL)

       if (myrank /= 0) call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       if (myrank /= 0) call hydroStep3d(grid, dt, nPairs, thread1, thread2, nBound, group, nGroup,doSelfGrav=doselfGrav)
       if (myRank == 1) call tune(6,"Hydrodynamics step")
       if (myrank /= 0) call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)
       if (myrank /= 0) call resetNh(grid%octreeRoot)


       if (myrank /= 0) then
          iUnrefine = iUnrefine + 1
          if (iUnrefine == 5) then
             if (myrankglobal == 1) call tune(6, "Unrefine grid")
             call unrefineCells(grid%octreeRoot, grid, nUnrefine,1.d-3)
             if (myrankglobal == 1) call tune(6, "Unrefine grid")
             iUnrefine = 0
          endif
          call evenUpGridMPI(grid, .true., .true.)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       endif

 

       call writeInfo("Calling photoionization loop",TRIVIAL)
       !       call testIonFront(grid%octreeRoot, grid%currentTime)
	       !call neutralGrid(grid%octreeRoot)

       if(dt /= 0.d0) then
          loopLimitTime = grid%currentTime + dt
       !if(looplimittime < 1.d5) then
       !  print *, "Running first photoionization sweep"
       ! loopLimitTime = 1.d15
          call setupNeighbourPointers(grid, grid%octreeRoot)
          call photoIonizationloopAMR(grid, source, nSource, nLambda,lamArray, 10, loopLimitTime, loopLimitTime, .True., .true.)

       else
          call setupNeighbourPointers(grid, grid%octreeRoot)
          call photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, 1, loopLimitTime, dt, .True., .true.)
       end if
          call writeInfo("Done",TRIVIAL)
       !end if
  
       
       if (myrank /= 0) then
          call calculateEnergyFromTemperature(grid%octreeRoot)
          call calculateRhoE(grid%octreeRoot, direction)
       endif

  

       if (myRank /= 0) then

          call evenUpGridMPI(grid,.false.,.true.)
       
          call refineGridGeneric(grid, 1.d-2)

          call writeInfo("Evening up grid", TRIVIAL)    
          call evenUpGridMPI(grid, .false.,.true.)
          call exchangeAcrossMPIboundary(grid, nPairs, thread1, thread2, nBound, group, nGroup)

       endif

!       call writeVtkFile(grid, "current.vtk", &
!            valueTypeString=(/"rho        ","HI        " ,"temperature" /))

       grid%currentTime = grid%currentTime + dt
       

       if (myRank == 1) write(*,*) "Current time: ",grid%currentTime

!Track the evolution of the ionization front with time
!          write(datFilename, '(a, i4.4, a)') "Ifront.dat"
!          call dumpStromgrenRadius(grid, datFileName, VECTOR(-3.3d8,  -3.3d8, 3.3d8), &
!               VECTOR(3.3d8, 3.3d8, 3.3d8), 1000)


       if (dumpThisTime) then

          !Thaw, to match time dumps of other codes
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
          !if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d12)) then
          !   deltaTForDump = 1.d12
	  !end if
          !if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d13)) then
          !   deltaTForDump = 1.d13
	  !end if
!          if(grid%geometry == "hii_test" .and. grid%currentTime >= (1.d14))then
!             deltaTForDump = 1.d14
 !         end if


          timeOfNextDump = timeOfNextDump + deltaTForDump
          grid%iDump = grid%iDump + 1

!          if(grid%geometry /= "hii_test") then
             write(mpiFilename,'(a, i4.4, a)') "dump_", grid%iDump,".grid"
             call writeAmrGrid(mpiFilename, .false., grid)
             write(mpiFilename,'(a, i4.4, a)') "dump_", grid%iDump,".vtk"
             call writeVtkFile(grid, mpiFilename, &
                  valueTypeString=(/"rho          ","HI           " , "temperature  ", &
                  "hydrovelocity","dust1        ","pressure     "/))
!          end if




!Track the evolution of the ionization front with time
       if(grid%geometry == "hii_test") then
          write(datFilename, '(a, i4.4, a)') "hii_test",grid%iDump,".dat"
          call dumpValuesAlongLine(grid, datFileName, VECTOR(-1.5d9,  -1.5d9, 1.5d9), &
               VECTOR(1.5d9, 1.5d9, 1.5d9), 1000)

          write(datFilename, '(a, i4.4, a)') "Ifront.dat"     
          call dumpStromgrenRadius(grid, datFileName, VECTOR(-1.5d9,  -1.5d9, 1.5d9), &
               VECTOR(1.5d9, 1.5d9, 1.5d9), 1000)
       end if


       if(grid%geometry == "radHydroTest") then
          write(datFileName,'(a,i4.4,a)') "radHydro",grid%iDump,".dat"
          call  dumpValuesAlongLine(grid, datFileName, VECTOR(0.d0,0.d0,0.0d0), &
          VECTOR(3.0d9, 0.d0, 0.0d0), 1000)

          write(datFilename, '(a, i4.4, a)') "Ifront.dat"
          call dumpStromgrenRadius(grid, datFileName, VECTOR(0.0d0,  0.0d0, 0.0d0), &
               VECTOR(3.0d9, 0.0d0, 0.0d0), 1000)

       end if


       if(grid%geometry == "bonnor") then
          write(datFilename, '(a, i4.4, a)') "bonnor",grid%iDump,".dat"
          call dumpValuesAlongLine(grid, datFileName, VECTOR(0.d0,  0.d0,0.d0), &
               VECTOR(-1.5d9, -0.d0, 0.d0), 1000)

       end if
	 
       endif
       
    enddo
  end subroutine radiationHydro


  subroutine photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, maxIter, tLimit, deltaTime, timeDep, monteCheck, &
       sublimate)
    use input_variables, only : quickThermal, inputnMonte, noDiffuseField, minDepthAMR, maxDepthAMR
    implicit none
    include 'mpif.h'
    integer :: myRank, ierr
    integer :: nMonte
    logical, optional :: sublimate
    logical :: doSublimate, timeDep, monteCheck
    real(double) :: tLimit
    real(double) :: deltaTime
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
!    integer :: nCellsInDiffusion
!    character(len=80) :: message
!    integer :: tempSubcell
    integer :: nlambda
    real :: lamArray(:)
    integer :: nSource
    type(SOURCETYPE) :: source(:), thisSource
    integer :: iSource
    type(VECTOR) :: rVec, uHat, rHat
    real(double) :: lCore
    integer :: iMonte
    integer :: subcell
    integer :: i, j
    logical :: escaped
    real(double) :: wavelength, thisFreq
    real :: thisLam
    type(VECTOR) :: octVec
    real(double) :: r
    integer :: ilam
    integer :: nInf
    real(double) :: kappaScadb, kappaAbsdb
    real(double) :: epsOverDeltaT
    integer :: nIter
    logical :: converged, thisThreadConverged, failed
!    real :: temp
    real(double) :: luminosity1, luminosity2, luminosity3
    real(double) :: photonPacketWeight
    real(double) :: fac
    real(double) :: tPhoton
    real(double) :: albedo
    integer :: maxIter
    type(SAHAMILNETABLE),save :: hTable, heTable
    type(RECOMBTABLE),save :: Hrecombtable

    real(double) :: freq(1000), dfreq(1000), spectrum(1000), nuStart, nuEnd
    real(double) :: r1, kappaAbsGas, kappaAbsDust, escat
    integer, parameter :: nTemp = 9
    real(double) :: tempstorage(nTemp)
    real(double) :: v, dustHeating
    real :: kappaP
    integer, parameter :: nFreq = 1000
    logical, save :: firsttime = .true.
    integer :: iMonte_beg, iMonte_end, nSCat
    type(octalWrapper), allocatable :: octalArray(:) ! array containing pointers to octals
    integer :: np, iOctal, iOctal_beg, iOctal_end, nOctal
    integer :: iThread, nThreads
    logical :: endLoop
    logical :: crossedMPIboundary
    integer :: tag = 41
    integer :: nTotScat, nPhot
    integer :: newThread
    logical, save :: firstTimeTables = .true.
    integer :: nEscaped, iSignal
    integer, parameter :: nTimes = 1
    logical :: photonsStillProcessing
    integer, allocatable :: nEscapedArray(:)
    integer :: status(MPI_STATUS_SIZE)

!================TOMS VARIABLES=======================
    real(double) :: deltaT  
!    integer :: nFailedRank, nFailedTotal, nFailedLast
!    integer :: deltaFails
!    real :: unconTLastIter = 0.0
!    real :: unconTThisIter, unconTThisIterRank
!    integer :: failCount
    logical :: anyUndersampled, undersampledTOT
    character(len=80) :: vtkFilename
    logical :: underSamFailed

    !optimisation variables
    !1integer :: nSaved, pSend
    !integer, parameter :: maxStore=1000
    !type(VECTOR) :: rVecStore(maxStore)
    !type(VECTOR) :: uHatStore(maxStore)
    !real(double) :: thisFreqStore(maxStore) 
    !real(double) :: tPhotonStore(maxStore) 
    !real(double) :: ppwStore(maxStore)
    !integer :: newThreadStore(maxStore) 
    !integer :: photonStackSize
    !type(VECTOR) :: photonStackrVec(maxStore)
    !type(VECTOR) :: photonStackuHat(maxStore)
    !real(double) :: photonStackFreq(maxStore)
    !real(double) :: photonStacktPhot(maxStore)
    !real(double) :: photonStackppw(maxStore)
    !integer :: photonStackDestination(maxStore)
    !logical :: readyToSend

!====================================================

    doSublimate = .true.
    if (PRESENT(sublimate)) doSublimate = sublimate

    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)

    allocate(nEscapedArray(1:nThreads-1))
    
!    write(*,*) "abundances ",grid%ion(1:5)%abundance

       nuStart = cSpeed / (1000.e4 * 1.d-8)
       nuEnd =  cSpeed / (10. * 1.d-8) ! 2.d0*maxval(grid%ion(1:grid%nIon)%nuThresh)

    do i = 1, nFreq
       freq(i) = log10(nuStart) + dble(i-1)/dble(nFreq-1) * (log10(nuEnd)-log10(nuStart))
       freq(i) = 10.d0**freq(i)
    enddo
    do i = 2, nFreq-1
       dfreq(i) = (freq(i+1)-freq(i-1))/2.d0
    enddo
    dfreq(1) = 2.d0*(freq(2)-freq(1))
    dfreq(nfreq) = 2.d0*(freq(nfreq)-freq(nfreq-1))


    if (firstTimeTables) then

       do i = 1, grid%nIon
          call addxSectionArray(grid%ion(i), nfreq, freq)
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
       lCore = lCore + source(i)%luminosity
    enddo



    if (myrankglobal == 1) write(*,'(a,1pe12.5)') "Total source luminosity (lsol): ",lCore/lSol

    if (writeoutput) then
       write(*,'(a,1pe12.1)') "Ionizing photons per cm^2 ", ionizingFlux(source(1))
    endif

!    call writeVtkFile(grid, "start.vtk", &
!         valueTypeString=(/"rho        ", "HI         ","temperature", "dust1      " /))


    !nmonte selector: Only works for fixed grids at present
    if (inputnMonte == 0) then
       
       if(.not. monteCheck) then
          !waymaker photoionization loop
          nmonte = 3000000
       else
          if(minDepthAMR == maxDepthAMR) then
             if(grid%octreeRoot%twoD) then
                nMonte = 10000.d0 * (4.d0**(maxDepthAMR))
             else if(grid%octreeRoot%threeD) then
                nMonte = 10.d0 * (8.d0**(maxDepthAMR))
		!nmonte = 100000
             else
                !nMonte = 100.d0 * 2**(maxDepthAMR)
                nmonte = 300000
             end if
          else
             call writeInfo("Non uniform grid, setting arbitrary nMonte", TRIVIAL)
             nMonte = 10000000
             write(*,*) "nMonte = ", nMonte
          end if
       end if
       
    else
       nMonte = inputnMonte
    endif

    nIter = 0
    
    if(grid%geometry == "Lexington") then
       maxIter = 20
    end if

    converged = .false.
    if (nSource > 1) &
         call randomSource(source, nSource, iSource, photonPacketWeight, lamArray, nLambda, initialize=.true.)
    do while(.not.converged)
       nIter = nIter + 1
       nInf=0
       

!       nCellsInDiffusion = 0
!      call defineDiffusionOnRosseland(grid,grid%octreeRoot, ndiff=nCellsInDiffusion)
!      write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
!      call writeInfo(message,IMPORTANT)


       epsoverdeltat = lcore/dble(nMonte)
       
       if (myrank /= 0) call zeroDistanceGrid(grid%octreeRoot)

       if (myrank == 1) write(*,*) "Running photoionAMR loop with ",nmonte," photons. Iteration: ",niter, maxIter

       if (myrank == 1) call tune(6, "One photoionization itr")  ! start a stopwatch

       iMonte_beg = 1
       iMonte_end = nMonte

       nTotScat = 0
       nPhot = 0


       call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!        photonStackppw = 0.d0
!       nSaved = 0
          if (myRank == 0) then
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
!                write(*,*) inOctal(grid%octreeRoot, rVec), "rvec ",rVec," dir ",uhat
 !               call amrGridValues(grid%octreeRoot, rVec, foundOctal=tempOctal, &
 !                    foundSubcell=tempsubcell)
 !               
 !               if (tempOctal%diffusionApprox(tempsubcell)) then
 !                  call randomWalk(grid, tempOctal, tempSubcell, thisOctal, Subcell, temp)
 !                  rVec = subcellCentre(thisOctal, subcell)
 !               endif
                
                
                call getWavelength(thisSource%spectrum, wavelength)                
                thisFreq = cSpeed/(wavelength / 1.e8)
                call findSubcellTD(rVec, grid%octreeRoot,thisOctal, subcell)

                iThread = thisOctal%mpiThread(subcell)

		!readyToSend=.true.
		!nSaved = 1
		!call MPI_SEND(nSaved , 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, ierr)		
		!call MPI_SEND(readyToSend , 1, MPI_LOGICAL, iThread, tag, MPI_COMM_WORLD, ierr)		
                call sendMPIPhoton(rVec, uHat, thisFreq, tPhoton,photonPacketWeight, iThread)
		!readyToSend=.false.
		!call MPI_SEND(readyToSend , 1, MPI_LOGICAL, iThread, tag, MPI_COMM_WORLD, ierr)
		!write(*,*) "Rank 0 sending photon with rvec ",rvec," to ", iThread
                nInf = nInf + 1
		     

             end do mainloop

             photonsStillProcessing = .true.
             do while(photonsStillProcessing)                   

                do iThread = 1, nThreads - 1
                   tempStorage(1) = HUGE(tempStorage(1))
                   tempStorage(2) = 0.d0
                   call MPI_SEND(tempStorage, nTemp, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD,  ierr)
                   call MPI_RECV(nEscapedArray(iThread), 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
                enddo
                nEscaped = SUM(nEscapedArray(1:nThreads-1))
                nEscapedGlobal = nEscaped
                if (nEscaped == nMonte) photonsStillProcessing = .false.


             end do

             do iThread = 1, nThreads - 1
                tempStorage(1) = HUGE(tempStorage(1))
                tempStorage(2) = HUGE(tempStorage(2))
                call MPI_SEND(tempStorage, nTemp, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD,  ierr)
             enddo

          else
             endLoop = .false.
             nEscaped = 0
	     
	     !needNewPhotonArray = .true.	     
             do while(.not.endLoop)
                crossedMPIboundary = .false.
		!readyToSend = .false.
                !HERE:CHECK STACK FOR PHOTONS - IF NONE THEN GET MORE
                !THAW
		
                !photonStackSize = SUM(photonStackppw)
		!print *, photonStackSize
		!i = 0
!		!stop
                !if(SUM(photonStackppw) == 0.d0) then
		!   print *, "A "
		!   
		!   do iThread = 0, nThreads-1
		!      print *, "Doing ", iThread, "of ", nThreads 
		!      if(iThread /= myRank .and. myRank /= 0) then
		
                !       call MPI_RECV(readyToSend ,1, MPI_LOGICAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
		!       
		
		!       if(readyToSend) then
		
                !        do While(readyToSend)
                !           call getNewMPIPhoton(rVec, uHat, thisFreq, tPhoton, photonPacketWeight, iSignal)
		!           i = i + 1
                !           photonStackrVec(i) = rVec
                !           photonStackuHat(i) = uHat
                !           photonStackFreq(i) = thisFreq
                !           photonStackppw(i) = photonPacketWeight
                !           photonStackDestination(i) = iSignal
                !           call MPI_RECV(readyToSend, 1, MPI_LOGICAL, iThread, tag, MPI_COMM_WORLD, status, ierr)                           
		!				   			   
		!
		!	end do
		!	print *, "sad face :(", myrank
                !       endif
		!       photonStackSize = i
		!      end if
		!      print *, "Problem? ^_^"
		!   end do
		   
                !else
		
		   !stop
                !   rVec = photonStackrVec(photonStackSize-1)
                !   uHat = photonStackuHat(photonStackSize-1)
                !   thisFreq = photonStackFreq(photonStackSize-1)
                !   tPhoton = photonStacktPhot(photonStackSize-1)
                !   photonPacketWeight = photonStackppw(photonStackSize-1)
                !   iSignal = photonStackDestination(photonStackSize-1)
		!   photonStackSize = photonStackSize - 1
		!   !Clear the array for another run
		!   if (photonStackSize == 0) then
		!      photonStackppw = 0.d0
		!   end if 
                !end if
                call getNewMPIPhoton(rVec, uHat, thisFreq, tPhoton, photonPacketWeight, iSignal)

		


                if (iSignal == 0) then
                   endLoop = .true.
                   goto 777
                endif
                if (iSignal == 1) then
                   call MPI_SEND(nEscaped, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD,  ierr)
                   goto 777
                endif

                ! here a new photon has been received, and we add its momentum

                call findSubcellTD(rVec, grid%octreeRoot,thisOctal, subcell)
                thisOctal%radiationMomentum(subcell) = thisOctal%radiationMomentum(subcell) + uHat * (epsOverDeltaT/cSpeed)
!                if (myrankGlobal == 1) write(*,*) "mom add new ",thisOctal%radiationMomentum(subcell)

                if (.not.endLoop) then

                   nScat = 0
                   escaped = .false.
                   do while(.not.escaped)
                      call toNextEventPhoto(grid, rVec, uHat, escaped, thisFreq, nLambda, lamArray, &
                           photonPacketWeight, epsOverDeltaT, nfreq, freq, tPhoton, tLimit, &
                           crossedMPIboundary, newThread)

                      
                      if (crossedMPIBoundary) then
                            !OLD STUFF
                         call sendMPIPhoton(rVec, uHat, thisFreq,tPhoton,photonPacketWeight, newThread)
                         goto 777

                         !HERE: SEND PHOTONS AS BEFORE
!THAW =============================================================================================
			 !Store the photon details
			 !nSaved = nSaved + 1
			 !rVecStore(nSaved) = rVec
			 !uHatStore(nSaved) = uHat
           		 !thisFreqStore(nSaved) = thisFreq
			 !tPhotonStore(nSaved) = tPhoton
			 !ppwStore(nSaved) = photonPacketWeight
			 !newThreadStore(nSaved) = newThread
			 !
			!!
                         !!once storage limit is reached, send photons
			 !if(nSaved /= 0) then
			 !  if(nSaved == maxStore .or. (nMonte-nEscapedGlobal) < (nThreads*maxStore)) then
			 !    readyToSend = .true.
			 !    !call MPI_SEND(nSaved , 1, MPI_INTEGER, newThreadStore(pSend), tag, MPI_COMM_WORLD, ierr)
			 !    do pSend=1, (nSaved)
			 !     call MPI_SEND(readyToSend , 1, MPI_LOGICAL, newThreadStore(pSend), tag, MPI_COMM_WORLD, ierr)
                         !     call sendMPIPhoton(rVecStore(pSend),uHatStore(pSend), thisFreqStore(pSend), & 
			 !     tPhotonStore(pSend), ppwStore(pSend), newThreadStore(pSend))
!
!!!			    !print *, "rank ", myRank, "sends photon to thread", newThread
			 !   end do
                         !   readyToSend = .false.
			 !   call MPI_SEND(readyToSend , 1, MPI_LOGICAL, newThreadStore(pSend), tag, MPI_COMM_WORLD, ierr)
	  		 !   nSaved = 0
			 !   goto 777
           		 !  else
			 !   readyToSend = .false.
	           	 !   call MPI_SEND(readyToSend , 1, MPI_LOGICAL, newThreadStore(pSend), tag, MPI_COMM_WORLD, ierr)
			 !  end if
			 !end if
	!		 !print *, nSaved
                        !g!oto 777			 
! =============================================================================================
                      endif

                      if (noDiffuseField) escaped = .true.
		
                      if (escaped) nEscaped = nEscaped + 1
		      
                      if (.not. escaped) then
                         
                         thisLam = (cSpeed / thisFreq) * 1.e8
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
                            
!                            write(*,*) "calling return kappa ",myrank
                            call returnKappa(grid, thisOctal, subcell, ilambda=ilam, &
                                 kappaAbsDust=kappaAbsDust, kappaAbsGas=kappaAbsGas, &
                                 kappaSca=kappaScadb, kappaAbs=kappaAbsdb, kappaScaGas=escat)
!                            write(*,*) "done return kappa ", myrank
                            
                            if ((thisFreq*hcgs*ergtoev) > 13.6) then ! ionizing photon
                               call randomNumberGenerator(getDouble=r1)

                               if (r1 < (kappaAbsGas / max(1.d-30,(kappaAbsGas + kappaAbsDust)))) then  ! absorbed by gas rather than dust
                                  call addLymanContinua(nFreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
                                  call addHigherContinua(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, GammaTableArray)
                                  call addHydrogenRecombinationLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
                                  !                        call addHeRecombinationLines(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid)
                                  call addForbiddenLines(nfreq, freq, spectrum, thisOctal, subcell, grid)
                               else
                                  call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
                               endif
                            else ! non-ionizing photon must be absorbed by dust
                               call addDustContinuum(nfreq, freq, dfreq, spectrum, thisOctal, subcell, grid, nlambda, lamArray)
                            endif
                            if (firsttime.and.(myrankglobal==1)) then
                               firsttime = .false.
                               open(67,file="pdf.dat",status="unknown",form="formatted")
                               do i = 1, nfreq
                                  write(67,*) freq(i), spectrum(i)
                               enddo
                               close(67)
                            endif
                            thisFreq =  getPhotonFreq(nfreq, freq, spectrum)
                            uHat = randomUnitVector() ! isotropic emission
!                            if (myrank == 1) write(*,*) "new uhat emission ",uhat

                            nScat = nScat + 1
                            
                            
!                            if ((thisFreq*hcgs*ergtoev) < 13.6) escaped = .true.


                         endif
                

                      endif
                      if (nScat > 100000) then
                         write(*,*) "Nscat exceeded 10000, forcing escape"
                         write(*,*) 1.e8*cspeed/thisFreq
                         write(*,*) albedo, kappaScaDb, kappaAbsdb,escat
                         
                         thisLam = (cSpeed / thisFreq) * 1.e8
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
                   enddo
                   nTotScat = nTotScat + nScat
                   nPhot = nPhot + 1
                endif
777             continue
             enddo
          endif

       if (myrank == 1) call tune(6, "One photoionization itr")  ! stop a stopwatch
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       epsOverDeltaT = (lCore) / dble(nMonte)
       if (myrank == 1) then
          if (nPhot > 0) then
             write(*,*) "Thread 1 had ",real(nTotScat)/real(nPhot), " scatters per photon"
          endif
       endif

       call  identifyUndersampled(grid%octreeRoot)

    np = 1
    firstTime = .true.

    allocate(octalArray(grid%nOctals))
    nOctal = 0
    call getOctalArray(grid%octreeRoot,octalArray, nOctal)
    if (nOctal /= grid%nOctals) then
       write(*,*) "Screw up in get octal array", nOctal,grid%nOctals
       stop
    endif

    ! default loop indices
    ioctal_beg = 1
    ioctal_end = nOctal


    if (myRank == 0) write(*,*) "Ninf ",ninf
    if (myrank == 1) call tune(6, "Temperature/ion corrections")

    if (writeoutput) &
         write(*,*) "Calculating ionization and thermal equilibria"

    if (myrank == 1) then
       i = 0
       j = 0 
       call testforZero(grid%octreeRoot, i, j)
       write(*,*) "photoion(subcell,1) is zero for ", 100.d0*real(i)/real(j), "% of subcells"
    endif

    if (myRank /= 0) then
       do iOctal =  iOctal_beg, iOctal_end
          
          thisOctal => octalArray(iOctal)%content
          
          if (dustOnly) then
             do subcell = 1, thisOctal%maxChildren
                if (.not.thisOctal%hasChild(subcell)) then
                   v = cellVolume(thisOctal, subcell)
                   call returnKappa(grid, thisOctal, subcell, kappap=kappap)
                   dustHeating = (epsOverDeltaT / (v * 1.d30))*thisOctal%distanceGrid(subcell) ! equation 14 of Lucy 1999
                   thisOctal%temperature(subcell) = ((pi/stefanBoltz) * dustHeating / (fourPi * kappaP))**0.25d0
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
       enddo


    endif



    !deallocate(octalArray)

    if (writeoutput) &
         write(*,*) "Finished calculating ionization and thermal equilibria"
    if (myrank == 1) call tune(6, "Temperature/ion corrections")

!       call defineDiffusionOnRosseland(grid,grid%octreeRoot)
!       call tune(6, "Gauss-Seidel sweeps")


!       nCellsInDiffusion = 0
!       call defineDiffusionOnUndersampled(grid%octreeroot, nDiff=nCellsInDiffusion)


!       call solveArbitraryDiffusionZones(grid)
!       call defneDiffusionOnRosseland(grid,grid%octreeRoot, nDiff=nCellsInDiffusion)
!       write(message,*) "Number of cells in diffusion zone: ", nCellsInDiffusion
!       call writeInfo(message,IMPORTANT)

!       call tune(6, "Gauss-Seidel sweeps")

!    call writeVtkFile(grid, "current.vtk", &
!         valueTypeString=(/"rho        ","HI         " ,"temperature" /))

    if(grid%geometry == "lexington") then
       if (myRank == 1) call dumpLexington(grid, epsoverdeltat)
    end if
if (.false.) then
       call dumpLexington(grid, epsoverdeltat)
       fac = 2.06e37

       luminosity1 = 0.d0
       call getHbetaLuminosity(grid%octreeRoot, grid, luminosity1)
       write(*,'(a,2f12.4)') "H beta :",luminosity1/1.e37,luminosity1/2.05e37
       fac = luminosity1
       

       call getForbiddenLineLuminosity(grid, "N II", 1.22d6, luminosity1)
       write(*,'(a,2f12.4)') "N II (122 um):",(luminosity1)/fac,(luminosity1)/(0.034*2.05e37)

       call getForbiddenLineLuminosity(grid, "N II", 6584.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "N II", 6548.d0, luminosity2)
       write(*,'(a,2f12.4)') "N II (6584+6548):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.730*2.05e37)

       call getForbiddenLineLuminosity(grid, "N III", 5.73d5, luminosity1)
       write(*,'(a,2f12.4)') "N III (57.3 um):",(luminosity1)/fac,(luminosity1+luminosity2)/(0.292*2.05e37)


       call getForbiddenLineLuminosity(grid, "O I", 6300.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O I", 6363.d0, luminosity2)
       write(*,'(a,2f12.4)') "O I (6300+6363):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.0086*2.05e37)

       call getForbiddenLineLuminosity(grid, "O II", 7320.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O II", 7330.d0, luminosity2)
       write(*,'(a,2f12.4)') "O II (7320+7330):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.029*2.05e37)

       call getForbiddenLineLuminosity(grid, "O II", 3726.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O II", 3729.d0, luminosity2)
       write(*,'(a,2f12.4)') "O II (3726+3729):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(2.03*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 5007.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O III", 4959.d0, luminosity2)
       write(*,'(a,2f12.4)') "O III (5007+4959):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(2.18*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 4363.d0, luminosity3)
       write(*,'(a,2f12.4)') "O III (4363):",(luminosity3)/1.e37,luminosity3/(0.0037*2.05e37)

       call getForbiddenLineLuminosity(grid, "O III", 518145.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "O III", 883562.d0, luminosity2)
       write(*,'(a,2f12.4)') "O III (52+88um):,",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/((1.06+1.22)*2.05e37)


       call getForbiddenLineLuminosity(grid, "Ne II", 1.28d5, luminosity1)
       write(*,'(a,2f12.4)') "Ne II (12.8um):",(luminosity1)/fac,luminosity1/(0.195*2.05e37)

       call getForbiddenLineLuminosity(grid, "Ne III", 1.56d5, luminosity1)
       write(*,'(a,2f12.4)') "Ne III (15.5um):",(luminosity1)/fac,luminosity1/(0.322*2.05e37)

       call getForbiddenLineLuminosity(grid, "Ne III", 3869.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "Ne III", 3968.d0, luminosity2)
       write(*,'(a,2f12.4)') "Ne III (3869+3968):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.085*2.05e37)


       call getForbiddenLineLuminosity(grid, "S II", 6716.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "S II", 6731.d0, luminosity2)
       write(*,'(a,2f12.4)') "S II (6716+6731):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.147*2.05e37)

       call getForbiddenLineLuminosity(grid, "S II", 4068.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "S II", 4076.d0, luminosity2)
       write(*,'(a,2f12.4)') "S II (4068+4076):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(0.008*2.05e37)

       call getForbiddenLineLuminosity(grid, "S III", 1.87d5, luminosity1)
       write(*,'(a,2f12.4)') "S III (18.7um):",(luminosity1)/fac,luminosity1/(0.577*2.05e37)

       call getForbiddenLineLuminosity(grid, "S III", 9532.d0, luminosity1)
       call getForbiddenLineLuminosity(grid, "S III", 9069.d0, luminosity2)
       write(*,'(a,2f12.4)') "S III (9532+9069):",(luminosity1+luminosity2)/fac,(luminosity1+luminosity2)/(1.22*2.05e37)
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    thisThreadConverged = .false.    
    failed = .false.

     anyUndersampled = .false.
     if(grid%geometry == "hii_test") then
        minCrossings = 10000
     else if(grid%geometry == "lexington") then
        minCrossings = 10000
     else
        minCrossings = 10000
     end if
   !Thaw - auto convergence testing I. Temperature, will shortly make into a subroutine
       if (myRank /= 0) then
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
                      
                      !    print *, "deltaT = ", deltaT
                      
                      if(deltaT > 1.5d-2) then
                         if (thisOctal%nCrossings(subcell) /= 0 .and. thisOctal%nCrossings(subcell) < minCrossings) then
                            anyUndersampled = .true.
                         endif
                      end if
                      
                      if(deltaT < 1.5d-2 .and. .not. failed) then
                         thisThreadConverged = .true.
                      else 
                         thisThreadConverged = .false.  
			 
                         if(deltaT /= 0.d0 .and. .not. failed) then
                            print *, "deltaT = ", deltaT
                            print *, "thisOctal%temperature(subcell) ", thisOctal%temperature(subcell)
                            print *, "thisOctal%TLastIter(subcell) ", thisOctal%TLastIter(subcell)
                            print *, "cell center ", subcellCentre(thisOctal,subcell)
                            print *, "nCrossings ", thisOctal%nCrossings(subcell)
                         end if
                         failed = .true.
                      end if
                      thisOctal%TLastIter(subcell) = thisOctal%temperature(subcell)
                   end if
                end if
             end do
             !Send result to master rank 
          end do

!       print *, "RANK ", myRank, " SENDS  DATA TO RANK 0"
!       print *, "thisThreadConverged = ", thisThreadConverged
       !Send converged information
       call MPI_SEND(thisThreadConverged , 1, MPI_LOGICAL, 0, tag, MPI_COMM_WORLD, ierr)
       call MPI_SEND(anyUndersampled, 1, MPI_LOGICAL, 0, tag, MPI_COMM_WORLD, ierr)

    else
        failed = .false.
        underSamFailed = .false.
        underSampledTOT = .false.

	!!!Rank 0, collate results and decide if converged 
        converged = .false.
        do iThread = 1 , (nThreadsGlobal-1)
              call MPI_RECV(thisThreadConverged,1, MPI_LOGICAL, iThread, tag, MPI_COMM_WORLD, status, ierr )
              call MPI_RECV(anyUndersampled, 1, MPI_LOGICAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
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

           do iThread = 1, (nThreadsGlobal-1)
              call MPI_SEND(converged, 1, MPI_LOGICAL, iThread, tag, MPI_COMM_WORLD, ierr) 
              call MPI_SEND(underSampledTOT, 1, MPI_LOGICAL, iThread, tag, MPI_COMM_WORLD, ierr)   
           end do
    end if

    if(myRank /= 0) then
       call MPI_RECV(converged, 1, MPI_LOGICAL, 0, tag, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(underSampledTOT, 1, MPI_LOGICAL, 0, tag, MPI_COMM_WORLD, status, ierr)

    end if

     call MPI_BARRIER(MPI_COMM_WORLD, ierr)

     deallocate(octalArray)    

!     converged = .false.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (niter >= maxIter) converged = .true. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    converged = .true.
! call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if(converged) then
     if(myRank == 0) then
      print *, "photoionization loop converged at iteration ", niter
     end if
  else if(underSampledTOT .and. nMonte < 2.5d9 .and. monteCheck) then
     if(myRank == 0) then
        print *, "Undersampled cell, increasing nMonte"
     end if 
      nMonte = nMonte *2.d0
  end if
  
!  call quickSublimate(grid%octreeRoot) ! do dust sublimation
  
    call torus_mpi_barrier

! if(tlimit /= 1.d20) then
!   write(vtkFilename,'(a,i2.2,a)') "photo",niter,".vtk"
!  call writeVtkFile(grid, vtkFilename, &
!      valueTypeString=(/"rho          ","HI           " , "temperature  ", &
!      "dust1        ","radmom       "/))
! call writeAMRgrid("tmp.grid",.false.,grid)
! end if
enddo


! thisOctal => grid%octreeRoot
! call writeInfo("Calculating continuum emissivities...",TRIVIAL)
!  call  calcContinuumEmissivity(grid, thisOctal, nlambda, lamArray)
! call writeInfo("Done.",TRIVIAL)

 call MPI_BARRIER(MPI_COMM_WORLD, ierr)


! if (writelucy) then
!    call writeAmrGrid(lucyfileout,.false.,grid)
! endif
end subroutine photoIonizationloopAMR


SUBROUTINE toNextEventPhoto(grid, rVec, uHat,  escaped,  thisFreq, nLambda, lamArray, photonPacketWeight, epsOverDeltaT, &
     nfreq, freq, tPhoton, tLimit, crossedMPIboundary, newThread)
  include 'mpif.h'
  integer :: myRank, ierr
   type(GRIDTYPE) :: grid
   type(VECTOR) :: rVec,uHat, octVec,thisOctVec, tvec, oldRvec
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
   logical :: ok, outofTime
   
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

    stillinGrid = .true.
    escaped = .false.
    
    crossedMPIboundary = .false.
    outOfTime = .false.

    photonMomentum = epsOverDeltaT / cSpeed
    thisLam = (cSpeed / thisFreq) * 1.e8
    call locate(lamArray, nLambda, real(thisLam), iLam)
    if ((ilam < 1).or.(ilam > nlambda)) then
       write(*,*) "ilam errro",ilam,thisLam
    endif

! select an initial random tau and find distance to next cell

    call findSubcellTD(rVec, grid%octreeRoot,thisOctal, subcell)

    call randomNumberGenerator(getDouble=r)
    tau = -log(1.0-r)

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
!   write(*,*) "thisTau ",thisTau,thisLam

! if tau > thisTau then the photon has traversed a cell with no interactions

!    write(*,*) "tau, thisTau ", tau, thistau

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
          stillingrid = .false.
          escaped = .true.
       endif
          
       octVec = rVec

       
! check whether the photon has  moved across an MPI boundary
        
      if (stillInGrid) then

          if (.not.octalOnThread(nextOctal, nextsubcell, myRank)) then
             stillInGrid = .false.
             escaped = .true.
             crossedMPIboundary = .true.
             newThread = nextOctal%mpiThread(nextSubcell)
             !print *, "next on thread? ", octalOnThread(nextOctal, nextSubcell, newThread)
          endif
       endif


       !if(crossedMPIBoundary) then
       !   stop
       !end if

! update the distance grid

       if (.not.outOfTime) then
          call updateGrid(grid, thisOctal, subcell, thisFreq, tVal, photonPacketWeight, ilam, nfreq, freq)
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
          call randomNumberGenerator(getDouble=r)
          tau = -log(1.0-r)
          call locate(lamArray, nLambda, real(thisLam), iLam)
          call amrGridValues(grid%octreeRoot, octVec, iLambda=iLam,lambda=real(thisLam), startOctal=thisOctal, &
               actualSubcell=subcell, kappaSca=kappaScadb, kappaAbs=kappaAbsdb, &
               grid=grid, inFlow=inFlow)
          oldOctal => thisOctal
          thisOctVec = octVec

! calculate the distance to the next cell


!          write(*,*) "disttocell rvec ", inSubcell(thisOCtal, subcell, rVec)
!          write(*,*) "disttocell octvec ", inSubcell(thisOCtal, subcell, octVec)


          !if (.not. inSubcell(thisOctal, subcell, rVec)) then
          !   write(*,*) "rvec ",rVec
          !   write(*,*) "octvec ", octvec
          !   write(*,*) "xmin xmax ", thisOctal%xmin, thisOctal%xmax
          !   write(*,*) "ymin ymax ", thisOctal%ymin, thisOctal%ymax
          !   write(*,*) "zmin zmax ", thisOctal%zmin, thisOctal%zmax
          !endif
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
                  dble(tval)*dble(tau)/thisTau, photonPacketWeight, ilam, nfreq, freq)
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

!    if (crossedMpiBoundary) write(*,*) "leaving tonexteventphoto with rvec ",rvec

 end subroutine toNextEventPhoto


  subroutine intersectCube(grid, posVec, i1,i2,i3,direction, tval)
   use vector_mod
   use grid_mod
   implicit none
   type(GRIDTYPE) :: grid
   type(VECTOR) :: direction
   type(VECTOR) :: posVec, norm(6), p3(6)
   real(oct) :: t(6),tval,denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)
   integer :: i1, i2, i3

   ok = .true.

   norm(1) = VECTOR(1., 0., 0.)
   norm(2) = VECTOR(0., 1., 0.)
   norm(3) = VECTOR(0., 0., 1.)
   norm(4) = VECTOR(-1., 0., 0.)
   norm(5) = VECTOR(0., -1., 0.)
   norm(6) = VECTOR(0., 0., -1.)

   p3(1) = VECTOR(grid%xAxis(i1+1), 0., 0.)
   p3(2) = VECTOR(0.,grid%yAxis(i2+1),0.)
   p3(3) = VECTOR(0.,0.,grid%zAxis(i3+1))
   p3(4) = VECTOR(grid%xAxis(i1), 0., 0.)
   p3(5) = VECTOR(0.,grid%yAxis(i2),0.)
   p3(6) = VECTOR(0.,0.,grid%zAxis(i3))

   thisOk = .true.
   
   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.
      endif
      if (t(i) < 0.) thisOk(i) = .false.
!      if (denom > 0.) thisOK(i) = .false.
 enddo




  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo

  if (j == 0) ok = .false.
   
  if (.not.ok) then
     write(*,*) i1, i2, i3
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble((grid%xAxis(2)-grid%xAxis(1))/1000.))


  if (tval == 0.) then
     write(*,*) i1, i2, i3,tval
     write(*,*) posVec
     write(*,*) grid%xAxis(i1),grid%yAxis(i2),grid%zAxis(i3)
     write(*,*) grid%xAxis(i1+1),grid%yAxis(i2+1),grid%zAxis(i3+1)
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > 2.*(grid%xAxis(2)-grid%xAxis(1))) then
!     write(*,*) "tval too big",tval,i1,i2,i3,posvec
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif


  end subroutine intersectCube 
  subroutine intersectCubeAMR(grid, posVec, direction, tval)
   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(in) :: direction
   real(oct), intent(out) :: tval
   !
   type(VECTOR) :: norm(6), p3(6)
   type(OCTAL),pointer :: thisOctal
   type(VECTOR) :: subcen, point
   integer :: subcell
   
   real(oct) :: t(6),denom(6)
   integer :: i,j
   logical :: ok, thisOk(6)


   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.

   norm(1) = VECTOR(1.0d0, 0.d0, 0.0d0)
   norm(2) = VECTOR(0.0d0, 1.0d0, 0.0d0)
   norm(3) = VECTOR(0.0d0, 0.0d0, 1.0d0)
   norm(4) = VECTOR(-1.0d0, 0.0d0, 0.0d0)
   norm(5) = VECTOR(0.0d0, -1.0d0, 0.0d0)
   norm(6) = VECTOR(0.0d0, 0.0d0, -1.0d0)

   p3(1) = VECTOR(subcen%x+thisOctal%subcellsize/2.0d0, subcen%y, subcen%z)
   p3(2) = VECTOR(subcen%x, subcen%y+thisOctal%subcellsize/2.0d0 ,subcen%z)
   p3(3) = VECTOR(subcen%x,subcen%y,subcen%z+thisOctal%subcellsize/2.0d0)
   p3(4) = VECTOR(subcen%x-thisOctal%subcellsize/2.0d0, subcen%y,  subcen%z)
   p3(5) = VECTOR(subcen%x,subcen%y-thisOctal%subcellsize/2.0d0, subcen%z)
   p3(6) = VECTOR(subcen%x,subcen%y,subcen%z-thisOctal%subcellsize/2.0d0)

   thisOk = .true.

   do i = 1, 6

      denom(i) = norm(i) .dot. direction
      if (denom(i) /= 0.0d0) then
         t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
      else
         thisOk(i) = .false.
         t(i) = 0.0d0
      endif
      if (t(i) < 0.) thisOk(i) = .false.
!      if (denom > 0.) thisOK(i) = .false.
 enddo

  
  j = 0
  do i = 1, 6
    if (thisOk(i)) j=j+1
  enddo

  if (j == 0) ok = .false.
   
  if (.not.ok) then
     write(*,*) "Error: j=0 (no intersection???) in lucy_mod::intersectCubeAMR. "
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  tval = minval(t, mask=thisOk)
  tval = max(tval * 1.001d0,dble(thisOctal%subCellSize/1000.))


  if (tval == 0.) then
     write(*,*) posVec
     write(*,*) direction%x,direction%y,direction%z
     write(*,*) t(1:6)
     stop
  endif

  if (tval > sqrt(3.)*thisOctal%subcellsize) then
!     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
!     write(*,*) "direction",direction
!     write(*,*) t(1:6)
!     write(*,*) denom(1:6)
  endif


  end subroutine intersectCubeAMR

  subroutine intersectCubeAMR2D(grid, posVec, direction, tval)

! this is to find a cell intersection for a 2D AMR grid
! which is essentially a 2D-grid that projects into
! a 3D cylindrical coordinate system


   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(inout) :: direction
   real(oct), intent(out) :: tval
   !
   type(OCTAL),pointer :: thisOctal
   type(VECTOR) :: subcen, point
   integer :: subcell
   real(double) :: compZ, currentZ
   real(double) :: distToZBoundary, distToXboundary
   real(oct) :: r1,r2,d,cosmu,x1,x2,distTor1,distTor2, theta, mu
   logical :: ok
   type(VECTOR) :: xHat, zHAt

   point = posVec

   call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   subcen =  subcellCentre(thisOctal,subcell)
   ok = .true.


   r1 = subcen%x - thisOctal%subcellSize/2.d0
   r2 = subcen%x + thisOctal%subcellSize/2.d0
   d = sqrt(point%x**2+point%y**2)
   xHat = VECTOR(point%x, point%y,0.d0)
   call normalize(xHat)

   cosmu =((-1.d0)*xHat).dot.direction
   call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
   if (.not.ok) then
      write(*,*) "Quad solver failed in intersectcubeamr2d"
      direction = randomUnitVector()
      x1 = thisoctal%subcellSize/2.d0
      x2 = 0.d0
   endif
   distTor2 = max(x1,x2)

   theta = asin(max(-1.d0,min(1.d0,r1 / d)))
   cosmu = xHat.dot.direction
   mu = acos(max(-1.d0,min(1.d0,cosmu)))
   distTor1 = 1.e30
   if (mu  < theta ) then
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in intersectcubeamr2d"
         direction = randomUnitVector()
         x1 = thisoctal%subcellSize/2.d0
         x2 = 0.d0
      endif
      distTor1 = max(x1,x2)
   endif
         
   distToXboundary = min(distTor1, distTor2)


   zHat = VECTOR(0.d0, 0.d0, 1.d0)
   compZ = zHat.dot.direction
   currentZ = point%z

   if (compZ /= 0.d0 ) then
      if (compZ > 0.d0) then
         distToZboundary = (subcen%z + thisOctal%subcellsize/2.d0 - currentZ ) / compZ
      else
         distToZboundary = abs((subcen%z - thisOctal%subcellsize/2.d0 - currentZ ) / compZ)
      endif
   else
      disttoZboundary = 1.e30
   endif

   tVal = min(distToZboundary, distToXboundary) +0.0001d0*grid%halfsmallestsubcell
   if (tVal > 1.e29) then
      write(*,*) tVal,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "z",currentZ
   endif
   if (tval < 0.) then
      write(*,*) tVal,compZ, distToZboundary,disttoxboundary
      write(*,*) "subcen",subcen
      write(*,*) "z",currentZ
   endif

  end subroutine intersectCubeAMR2D

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
                   write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
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
                 write(*,*) "Undersampled cell",thisOctal%ncrossings(subcell)
                end if
             endif
          endif
       endif
    enddo
  end subroutine calculateIonizationBalanceTimeDep

  subroutine quickThermalCalc(thisOctal)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell

    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then
!          thisOctal%temperature(subcell) = 100.d0 + (20000.d0-100.d0) * thisOctal%ionFrac(subcell,2)
          thisOctal%temperature(subcell) = 10.d0 + (10000.d0-10.d0) * thisOctal%ionFrac(subcell,2)

          !if(thisOctal%ionFrac(subcell,2) /= 0.d0 .and. thisOctal%ionFrac(subcell, 2) /= 1.d-30) then
          !   print *, thisOctal%ionFrac(subcell,2)
          !end if

       endif

    enddo
  end subroutine quickThermalCalc

  subroutine calculateThermalBalance(grid, thisOctal, epsOverDeltaT)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    real(double) :: epsOverDeltaT
    real(double) :: totalHeating
    integer :: subcell
    logical :: converged, found
    real :: t1, t2, tm
    real, parameter :: Tlow = 100.
    real(double) :: y1, y2, ym, Hheating, Heheating, dustHeating
    real :: deltaT
    real :: underCorrection = 1.
    integer :: nIter
    

    do subcell = 1, thisOctal%maxChildren

       if (.not.thisOctal%hasChild(subcell)) then
          
          if (thisOctal%inflow(subcell)) then



             if (.not.thisOctal%undersampled(subcell)) then
                
                
                call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
!                if (thisOctal%ionFrac(subcell,1) > 0.9d0) write(*,*) "total heating ",totalheating
                if (totalHeating < 1.d-30) then
                   thisOctal%temperature(subcell) = tLow
                else
                   nIter = 0
                   converged = .false.
                   

!                   if (thisOctal%dustTypeFraction(subcell,1) > 0.9) then
!                      t1 = tlow
!                      t2 = 2000.
!                   else
                      t1 = 100.
                      t2 = 30000.
 !                  endif
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
    real, parameter :: Tlow = 100.
    real(double) :: y1, y2, ym, Hheating, Heheating, dustHeating
    integer :: nIter
    real :: t1, t2, tm
                
    call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDeltaT)
    !                if (thisOctal%ionFrac(subcell,1) > 0.9d0) write(*,*) "total heating ",totalheating
    nIter = 0
    converged = .false.
    
    
    t1 = 100.
    t2 = 30000.
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
                   thisOctal%temperature(subcell) = Teq
                else

                   nTimes = 10
                
                   smallDeltaT = deltaTime /dble(nTimes)
                   
                   do i = 1, nTimes
                      energy = energy + (totalHeating - totalcooling)*smallDeltaT
                      thisOctal%temperature(subcell) = (2.d0/3.d0)*energy*(mu*mhydrogen)/(thisOctal%rho(subcell)*kerg)
                      totalCooling = HHecooling(grid, thisOctal, subcell, thisOctal%temperature(subcell)) 
                   enddo
                endif

             endif
          endif
       endif
    enddo
  end subroutine calculateThermalBalanceTimeDep


  function hRecombination(temperature) result (rate)
    real :: temperature
    real(double) :: rate, t0, t1, b

    t0 = sqrt(temperature/3.148d0)
    t1 = sqrt(temperature/7.036d5)
    b = 0.7480

    rate = 7.982d-11 / (t0*(1.d0+t0)**(1.d0-b) * (1.d0 + t1)**(1.d0+b))
  end function hRecombination

  function recombToGround(temperature) result (alpha1)
    real :: temperature
    real(double) :: alpha1

    alpha1 = 1.58d-13 * (temperature/1.d4)**(-0.53d0)  ! kenny's photo paper equation 24
  end function recombToGround

  
  function HHeCooling(grid, thisOctal, subcell, temperature, debug) result (coolingRate)
    use input_variables, only : dustOnly, hOnly
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
          
          coolColl = (ch12*exp(-th12/TeUsed)*Te4**ex12 + &
               & ch13*exp(-th13/TeUsed)*Te4**ex13) * &
               & hcRyd*(nh-nhii)*Ne
          
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
          print *, "=========================================="
          print *, "n ", n
          print *, "rootTbetaH(n) ", rootTbetaH(n)
          print *, "rootTbetaH(n+1) ", rootTbetaH(n+1)
          print *, "fac ", fac
          print *, "thisLogT ", thisLogT
          print *, "logT(n) ", logT(n)
          print *, "logT(n+1) ", logT(n+1)
          print *, "=========================================="
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
       
       if(coolingRate < 0.d0) then
          print *, "coolingRate 2 ", coolingRate
          !stop
       end if
 !betaH is negative in my model!
       
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
       !    write(*,*) coolingRate,crate,coolingRate/(coolingrate+crate)
       coolingRate = coolingRate + crate
       
    endif

    call returnKappa(grid, thisOctal, subcell, kappap=kappap)
       
    dustCooling = fourPi * kappaP * (stefanBoltz/pi) * temperature**4

    coolingRate = coolingRate + dustCooling

  end function HHeCooling

  subroutine updateGrid(grid, thisOctal, subcell, thisFreq, distance, photonPacketWeight, ilambda, nfreq, freq)
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
    kappaAbs = 0.d0; kappaAbsDust = 0.d0
    thisOctal%nCrossings(subcell) = thisOctal%nCrossings(subcell) + 1


    fac = distance * photonPacketWeight / (hCgs * thisFreq)

    call locate(freq, nFreq, thisFreq, iFreq)

    do i = 1, grid%nIon

       xSec = returnxSec(grid%ion(i), thisFreq, iFreq=iFreq)
!       call phfit2(grid%ion(i)%z, grid%ion(i)%n, grid%ion(i)%outerShell , e , xsec)
       if (xSec > 0.d0) then
          thisOctal%photoIonCoeff(subcell,i) = thisOctal%photoIonCoeff(subcell,i) &
               + fac * xSec

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
        
        
!        write(*,*) "Photoionization rate ", (epsOverDeltaT / (v*1.d30))*thisOctal%photoIonCoeff(subcell,iIon)
!        write(*,*) "Recombination rate ",recombRate(grid%ion(iion), temperature) * thisOctal%ne(subcell)
!        write(*,*) " " 
        
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
        
!        write(*,*) "Photoionization rate ", (epsOverDeltaT / (v*1.d30))*thisOctal%photoIonCoeff(subcell,iIon)
!        write(*,*) "Recombination rate ",recombRate(grid%ion(iion), temperature) * thisOctal%ne(subcell)
!        write(*,*) " "

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

!There is something wrong in this subroutine

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

!THAW=================

     if (SUM(thisOctal%ionFrac(subcell,iStart:iEnd)) /= 0.d0) then
        thisOctal%ionFrac(subcell,iStart:iEnd) = &
             max(1.d-50,thisOctal%ionFrac(subcell,iStart:iEnd))/SUM(thisOctal%ionFrac(subcell,iStart:iEnd))
     else
        thisOctal%ionFrac(subcell,iStart:iEnd) = 1.d-50
     endif

!====================
     
        thisOctal%ionFrac(subcell, iIon) = min(max(thisOctal%ionFrac(subcell,iIon),ionFrac(iIon)),1.d0)
        if(thisOctal%ionFrac(subcell, iIon) == 0.d0) then
           print *, ionFrac(iend)
           stop
        end if

        if(thisOctal%ionFrac(subcell, iIon) == 0.d0) then
           thisOctal%ionFrac(subcell, iIon) = 1.d-50
        end if

     enddo

!Thaw
     thisOctal%ionFrac(subcell, iEnd) = 1.d0 - SUM(thisOctal%ionFrac(subcell,iStart:iEnd-1))

     thisOctal%ionFrac(subcell, iend) = min(max(thisOctal%ionFrac(subcell,iend),ionFrac(iend)),1.d0)

     if(thisOctal%ionFrac(subcell, iend) == 0.d0) then
        thisOctal%ionFrac(subcell, iend) = 1.d-50
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

real function vernerFerland(a, t, t0, t1, b)
  real(double):: a, t, t0, t1, b

! based on Verner and Ferland (1996) ApJS 103 467

  vernerFerland = a / (sqrt(t/t0) * (1.d0+sqrt(t/t0))**(1.d0-b) * (1.d0+sqrt(t/t1))**(1.d0+b))
end function vernerFerland



function nussbaumerStorey1983(t, a, b, c, d, f) result(rate)

! based on Nussbaumer & Storey 1983 AA 126 75
  real(double) :: t, a, b, c, d, f, rate

  if (t < 0.05d0) then
     rate = tiny(rate)
  else
     rate = (1.d-12) * ((a/t) + b + (c*t) + d*(t**2))*(t**(-1.5d0))*exp(-f/t)
  endif

end function nussbaumerStorey1983

function ppb1991(z, a, b, c, d, temperature) result(rate)

! radiative recombination rates
! based on Pequignot, Petitjean & Boisson (1991) A&A 251 680

  real(double) :: z, a, b, c, d, temperature, t, rate

  t = 1.d-4 * temperature / z**2

  rate = 1.d-13 * z * (a*t**b)/(1.d0+c*t**d)

end function ppb1991

function svs1982(t, alpharad, xrad) result (rate)

! radiative recombination rates based on
! Shull and Van Steenberg, 1992, ApJS, 48, 95

  real(double) :: t, alpharad, xrad, rate

  rate = alpharad * (t /1.d4)**(-xrad)
end function svs1982


subroutine dumpLexington(grid, epsoverdt)
  type(GRIDTYPE) :: grid
  type(OCTAL), pointer :: thisOctal
  integer :: subcell
  integer :: i, j
  real(double) :: r, theta, phi
  real :: t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
  real(double) :: oirate, oiirate, oiiirate, oivrate
  real(double) :: v, epsoverdt
  type(VECTOR) :: octVec
  real :: fac
  real(double) :: hHeating, heHeating, totalHeating, heating, nh, nhii, nheii, ne
  real(double) :: cooling, dustHeating
  real :: netot
  open(20,file="lexington.dat",form="formatted",status="unknown")
  open(21,file="orates.dat",form="formatted",status="unknown")
  open(22,file="ne.dat",form="formatted",status="unknown")

  do i = 1, 50
     r = (1.+7.d0*dble(i-1)/49.d0)*pctocm/1.e10

     t=0;hi=0; hei=0;oii=0;oiii=0;cii=0;ciii=0;civ=0;nii=0;niii=0;niv=0;nei=0;neii=0;neiii=0;neiv=0;ne=0.

     oirate = 0; oiirate = 0; oiiirate = 0; oivrate = 0
     heating = 0.d0; cooling = 0.d0
     do j = 1, 1000
661     continue
        call randomNumberGenerator(getDouble=theta)
        theta = theta * Pi
        call randomNumberGenerator(getDouble=phi)
        phi = phi * twoPi
        
        octVec = VECTOR(r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta))
        
        
        call amrgridvalues(grid%octreeRoot, octVec,  foundOctal=thisOctal, foundsubcell=subcell)
!        if (thisOctal%mpiThread(subcell) /= myRankGlobal) goto 661
        if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) goto 661

        nHii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,2) * grid%ion(2)%abundance
        nHeii = thisOctal%nh(subcell) * thisOctal%ionFrac(subcell,4) * grid%ion(4)%abundance
        nh = thisOctal%nh(subcell)
        ne = thisOctal%ne(subcell)

!        cooling = cooling + HHeCooling(grid, thisOctal, subcell, thisOctal%temperature(subcell))

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
        NIV = NIV + thisOctal%ionfrac(subcell,returnIonNumber("N IV", grid%ion, grid%nIon))
        NeI = NeI + thisOctal%ionfrac(subcell,returnIonNumber("Ne I", grid%ion, grid%nIon))
        NeII = NeII + thisOctal%ionfrac(subcell,returnIonNumber("Ne II", grid%ion, grid%nIon))
        NeIII = NeIII + thisOctal%ionfrac(subcell,returnIonNumber("Ne III", grid%ion, grid%nIon))
        NeIV = NeIV + thisOctal%ionfrac(subcell,returnIonNumber("Ne IV", grid%ion, grid%nIon))
        netot = netot + thisOctal%ne(subcell)
        call getHeating(grid, thisOctal, subcell, hHeating, heHeating, dustHeating, totalHeating, epsOverDT)
        heating = heating + totalHeating
!        fac = thisOctal%nh(subcell) * returnAbundance(8) !* thisOctal%ionfrac(subcell,returnIonNumber("O I", grid%ion, grid%nIon))
        fac = 1.
        oirate = oirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O I", grid%ion, grid%nIon)))
        oiirate = oiirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O II", grid%ion, grid%nIon)))
        oiiirate = oiiirate + &
             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O III", grid%ion, grid%nIon)))
!        oivrate = oivrate + &
!             fac*((epsOverDT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell,returnIonNumber("O IV", grid%ion, grid%nIon)))

        t  = t + thisOctal%temperature(subcell)
     enddo


     hi = hi / 1000.; hei = hei/1000.; oii = oii/1000.; oiii = oiii/1000.; cii=cii/1000.
     ciii = ciii/1000; civ=civ/1000.; nii =nii/1000.; niii=niii/1000.; niv=niv/1000.
     nei=nei/1000.;neii=neii/1000.; neiii=neiii/1000.; neiv=neiv/1000.;t=t/1000.
     netot = netot / 1000.

     oirate = oirate / 1000.
     oiirate = oiirate / 1000.
     oiiirate = oiiirate / 1000.
     oivrate = oivrate / 1000.
     heating = heating / 1000.
     cooling = cooling / 1000.

     hi = log10(max(hi, 1e-10))
     hei = log10(max(hei, 1e-10))
     oii = log10(max(oii, 1e-10))
     oiii = log10(max(oiii, 1e-10))
     cii = log10(max(cii, 1e-10))
     ciii = log10(max(ciii, 1e-10))
     civ = log10(max(civ, 1e-10))
     nii = log10(max(nii, 1e-10))
     niii = log10(max(niii, 1e-10))
     niv= log10(max(niv, 1e-10))
     nei = log10(max(nei, 1e-10))
     neii = log10(max(neii, 1e-10))
     neiii = log10(max(neiii, 1e-10))
     neiv = log10(max(neiv, 1e-10))
     ne = log10(max(ne,1.d-10))


     write(21,'(f5.3,1p,6e12.3,0p)') r*1.e10/pctocm,heating,cooling,oirate,oiirate,oiiirate,oivrate

     write(20,'(f5.3,f9.1,  14f8.3)') &
          r*1.e10/pctocm,t,hi,hei,oii,oiii,cii,ciii,civ,nii,niii,niv,nei,neii,neiii,neiv
     write(22,*) r*1.e10/pctocm,netot
  enddo
  close(20)
  close(21)
  close(22)
end subroutine dumpLexington

subroutine getChargeExchangeRecomb(parentIon, temperature, nhi, recombRate)
  type(IONTYPE) :: parentIon
  real(double) :: nhi
  real(double), intent(out) :: recombRate
  real :: t4, a, b, c, d
  real :: temperature

  recombRate  = 0.d0

  select case(parentIon%z)
     case(7)
        select case(parentIon%n)
           case(4) ! N IV
              t4 = temperature / 1.e4
              a = 3.05e-10
              b = 0.60
              c = 2.65
              d = -0.93
              recombRate  = kingdonFerland96(t4, a, b, c, d)
        end select
     case(8)
        select case(parentIon%n)
           case(7) ! O II
              t4 = temperature / 1.e4
              a = 1.04e-9
              b = 3.15e-2
              c = -0.61
              d = -9.73
              recombRate  = kingdonFerland96(t4, a, b, c, d)
           case(6) ! OIII
              t4 = temperature / 1.e4
              a = 1.04e-9
              b = 0.27
              c = 2.02
              d = -5.92
              recombRate  = kingdonFerland96(t4, a, b, c, d)
       end select
   end select

  recombRate = recombRate * nhi

end subroutine getChargeExchangeRecomb

subroutine getChargeExchangeIon(parentIon, temperature,  nHII, IonRate)
  type(IONTYPE) :: parentIon
  real(double) :: nHII
  real(double), intent(out) :: ionRate
  real :: temperature
  real :: t4, a, b, c, d

  IonRate = 0.d0

  select case(parentIon%z)
     case(7)
        select case(parentIon%n)
           case(7) ! N I charge exchange Ionization
              t4 = temperature / 1.e4
              a = 4.55e-12
              b = -0.29
              c = -0.92
              d = -8.38
              ionRate  = kingdonFerland96(t4, a, b, c, d)
        end select
     case(8)
        select case(parentIon%n)
           case(8) ! O I charge exchange ionization
              t4 = temperature / 1.e4
              a = 7.40e-11
              b = 0.47
              c = 24.37
              d = -0.74
              ionRate  = kingdonFerland96(t4, a, b, c, d)
       end select
   end select

  ionRate = ionRate * nhii
end subroutine getChargeExchangeIon


function kingdonFerland96(t4, a, b, c, d) result (alpha)
  real :: alpha
  real :: t4, a, b, c, d
  alpha = a*(t4**b)*(1.+c*exp(d*t4))
end function kingdonFerland96
  

subroutine getCollisionalRates(thisIon, iTransition, temperature, excitation, deexcitation)
  real, intent(out) :: excitation, deexcitation
  type(IONTYPE) :: thisIon
  integer :: iTransition, i
  real :: temperature
  real :: thisGamma
  real :: t , fac
  real :: boltzFac

  t = max(min(real(thisIon%transition(iTransition)%t(thisIon%transition(iTransition)%ngamma)),temperature), &
       real(thisIon%transition(iTransition)%t(1)))
  call locate(thisIon%transition(iTransition)%t, thisIon%transition(iTransition)%ngamma, t, i)
  fac = (t - thisIon%transition(iTransition)%t(i))/(thisIon%transition(iTransition)%t(i+1) - thisIon%transition(iTransition)%t(i))
  thisGamma = thisIon%transition(iTransition)%gamma(i) + &
       fac * (thisIon%transition(iTransition)%gamma(i+1) - thisIon%transition(iTransition)%gamma(i))

  boltzFac =  exp(-thisIon%transition(iTransition)%energy / (kev*temperature))

  fac = (8.63e-6 / sqrt(temperature)) * thisGamma

  deexcitation =  fac / thisIon%level(thisIon%transition(iTransition)%j)%g

  excitation =  fac / thisIon%level(thisIon%transition(iTransition)%i)%g * boltzFac

!  excitation = fac * boltzFac
!  deexcitation = fac * thisIon%level(thisIon%transition(iTransition)%j)%g &
!       / thisIon%level(thisIon%transition(iTransition)%i)%g

!  if (thisIon%species == "O III") then
!     do i = 1, thisIon%transition(1)%ngamma
!        write(*,*) i,thisIon%transition(1)%t(i),thisIon%transition(1)%gamma(i)
!     enddo
!     if (iTransition ==1) then
!        write(*,*) thisGamma, boltzFac, fac, excitation, deexcitation
!     endif
!  endif


end subroutine getCollisionalRates


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

recursive subroutine quickSublimate(thisOctal)
  use input_variables, only : hOnly
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call quickSublimate(child)
                exit
             end if
          end do
       else

          if(.not. hOnly) then
             if ((thisOctal%ionFrac(subcell,1) < 0.1).or.(thisOctal%temperature(subcell) > 1500.)) then
                thisOctal%dustTypeFraction(subcell,:) = 0.d0
             else
                thisOctal%dustTypeFraction(subcell,:) = 1.d0
                thisOctal%ne(subcell) = tiny(thisOctal%ne(subcell))
                thisOctal%ionFrac(subcell,1) = 1.d0
                thisOctal%ionFrac(subcell,2) = tiny(thisOctal%ionFrac(subcell,2))
             endif

          else
             thisOctal%dustTypeFraction(subcell,:) = 0.d0
          end if
       endif
    enddo
  end subroutine quickSublimate


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
  

subroutine solvePops(thisIon, pops, ne, temperature, debug)
  type(IONTYPE) :: thisIon
  real(double) :: ne
  real, intent(out) :: pops(:)
  real :: temperature
  real(double), allocatable :: matrixA(:,:), MatrixB(:), tempMatrix(:,:), qeff(:,:),  rates(:)
  integer :: n, iTrans, i, j
  real :: excitation, deexcitation, arateji
  logical :: ok
  logical, optional :: debug

  n = thisIon%nLevels
  allocate(matrixA(1:n, 1:n), matrixB(1:n), tempMatrix(1:n,1:n), qeff(1:n,1:n), rates(1:n))

  matrixA = 1.d-30
  matrixB = 0.d0


  call getRecombs(rates, thision, dble(temperature))

  matrixA(1,:) = 1.d0
  matrixB(1) = 1.d0

  matrixB(2:n) = 0.d0 !rates(2:n) * ne * (nh * ionFrac * thisIon%abundance)
  do iTrans = 1, thisIon%nTransitions
     i = thision%transition(itrans)%i
     j = thision%Transition(itrans)%j
     call getCollisionalRates(thisIon, iTrans, temperature, excitation, deexcitation)
     qeff(i,j) = excitation
     qeff(j,i) = deexcitation
  enddo

  do i = 2, n
     do j = 1, n

        do iTrans = 1, thisIon%nTransitions
           if (((i == thision%transition(itrans)%i).and.(j == thision%Transition(itrans)%j)).or. &
               ((i == thision%transition(itrans)%j).and.(j == thision%Transition(itrans)%i))) then
              arateji =  thisIon%transition(iTrans)%a
              
              matrixA(i,j) = matrixA(i,j) + ne * qeff(j, i)
              matrixA(i,i) = matrixA(i,i) - ne * qeff(i, j)
              if (j > i) then
                 matrixA(i,j) = matrixA(i,j) + arateji
              else
                 matrixA(i,i) = matrixA(i,i) - arateji
              endif

           endif
        enddo

     enddo
  enddo

  tempMatrix = matrixA

  if (PRESENT(debug)) then
     if (debug) then
        do i = 1, n
           write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
        enddo
     endif
  endif
  
  call luSlv(matrixA, matrixB)

  matrixB(1:n) = matrixB(1:n) / SUM(matrixB(1:n))
  
  ok = .true.

  if (.not.ok) then
     write(*,*) "Population solver failed for: ",thisIon%species
     write(*,*) matrixB(1:n)
     write(*,*) "nlevels",thisIon%nLevels,"ntrans",thisIon%nTransitions 
     write(*,*) "temp",temperature,"ne",ne
     
     do i = 1, n
        write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
     enddo
     matrixB = 0.d0
     matrixB(1) = 1.d0
     write(*,*) "Setting pops to ground state"
  endif

  do i = 1, n
     pops(i) = max(1.d-30,matrixB(i))
  enddo

  deallocate(matrixA, matrixB, tempMatrix, qeff, rates)

end subroutine solvePops

subroutine calcHeRecombs(te, alpha1, alpha21s, alpha21p, alpha23s)
  real :: te, t
  real(double) :: alpha1, alpha21s, alpha21p, alpha23s

  t = te / 1.e4

  alpha1 = 1.54e-13 * t**(-0.486)

  alpha21s = 2.06e-14 * t**(-0.676)

  alpha21p = 4.17e-14 * t**(-0.861)
  
  alpha23s = 2.10e-13 * t**(-0.778)
end subroutine calcHeRecombs


subroutine createSahaMilneTables(hTable, heTable)
  type(SAHAMILNETABLE), intent(out) :: hTable, heTable
  real(double) :: nu0_h, nu0_he, nufinal_h, nufinal_he
  integer :: nFreq, nTemp
  integer :: i, j
  real :: e
  real(double) :: t1, t2
  real :: hxsec, hexsec
  real(double) :: dfreq, jnu

  nFreq = 10000
  nTemp = 100
  hTable%nFreq = nFreq
  heTable%nFreq = nFreq
  hTable%nTemp = nTemp
  heTable%nTemp = nTemp

  allocate(hTable%freq(1:nFreq), heTable%freq(1:nFreq))
  allocate(hTable%emissivity(1:nTemp), heTable%emissivity(1:ntemp))
  allocate(hTable%temp(1:ntemp), heTable%temp(1:ntemp))
  allocate(hTable%Clyc(1:nTemp,1:nFreq))
  allocate(heTable%Clyc(1:nTemp,1:nFreq))

  nu0_h = 13.6d0/ergtoev/hcgs
  nu0_he = 24.59d0/ergtoev/hcgs

  nufinal_h = 2.d0*nu0_h
  nufinal_he = 2.d0*nu0_he

  do i = 1, nFreq
     hTable%freq(i) = log10(nu0_h) + (log10(nuFinal_h)-log10(nu0_h))*dble(i-1)/dble(nFreq-1)
     heTable%freq(i) = log10(nu0_he) + (log10(nuFinal_he)-log10(nu0_he))*dble(i-1)/dble(nFreq-1)
  enddo
  hTable%freq(1:hTable%nFreq) = 10.d0**hTable%freq(1:hTable%nfreq)
  heTable%freq(1:hTable%nFreq) = 10.d0**heTable%freq(1:hTable%nfreq)

  t1 = 5000.d0
  t2 = 20000.d0

  do i = 1, nTemp
     hTable%temp(i) = t1 + (t2-t1)*dble(i-1)/dble(nTemp-1)
     heTable%temp(i) = t1 + (t2-t1)*dble(i-1)/dble(nTemp-1)
  enddo

  do i = 1, nTemp
     hTable%Clyc(i,1) = 0.d0
     heTable%Clyc(i,1) = 0.d0
     do j = 2, nFreq


        e = hTable%freq(j) * hcgs* ergtoev
        call phfit2(1, 1, 1 , e , hxsec)

        dFreq = hTable%freq(j)-hTable%freq(j-1)
        jnu = ((hcgs*hTable%freq(j)**3)/(cSpeed**2)) * ((hcgs**2) /(twoPi*mElectron*Kerg*hTable%temp(i)))**(1.5d0) * &
             dble(hxsec/1.d10) *  exp(-hcgs*(hTable%freq(j)-nu0_h)/(kerg*hTable%temp(i)))
        hTable%Clyc(i,j) = hTable%Clyc(i,j-1) + jnu * dfreq


        e = heTable%freq(j) * hcgs * ergtoev
        call phfit2(2, 2, 1 , e , hexsec)

        dFreq = heTable%freq(j)-heTable%freq(j-1)
        jnu = 2.d0*((hcgs*heTable%freq(j)**3)/(cSpeed**2)) * ((hCgs**2) / (twoPi*mElectron*kerg*heTable%temp(i)))**(1.5d0) * &
             dble(hexsec/1.d10) * exp(-hcgs*(heTable%freq(j)-nu0_he)/(kerg*heTable%temp(i)))
        heTable%Clyc(i,j) = heTable%Clyc(i,j-1) + jnu * dfreq
     enddo
     hTable%emissivity(i) = hTable%Clyc(i,hTable%nFreq)
     heTable%emissivity(i) = heTable%Clyc(i,heTable%nFreq)
     hTable%Clyc(i,1:hTable%nFreq) = hTable%Clyc(i,1:hTable%nFreq) / hTable%Clyc(i,hTable%nFreq)
     heTable%Clyc(i,1:heTable%nFreq) = heTable%Clyc(i,1:heTable%nFreq) / heTable%Clyc(i,heTable%nFreq)
  end do

!  do j = 1, nFreq
!     write(99  ,*) htable%freq(j),htable%clyc(50,j)
!  enddo


end subroutine createSahaMilneTables

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

subroutine getRecombs(rates, thision, temperature)  
  real(double) :: rates(:)
  type(IONTYPE) :: thisIon
  real(double) :: temperature

  select case(thisIon%species)
     case("O II")
        rates(1) = 0.d0
        rates(2) = rateFit(7.218d0, -0.575d0, temperature)
        rates(3) = rateFit(4.812d0, -0.575d0, temperature)
        rates(4) = rateFit(3.581d0, -0.495d0, temperature)
        rates(5) = rateFit(1.790d0, -0.495d0, temperature)
     case DEFAULT
        rates = 0.d0
   end select
 end subroutine getRecombs


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
!  dustHeating = (epsOverDeltaT / (v * 1.d30))*thisOctal%distanceGrid(subcell) ! equation 14 of Lucy 1999
  totalHeating = (Hheating + HeHeating + dustHeating)
  
end subroutine getHeating

function rateFit(a,b,t) result (rate)
  real(double) :: a, b, t, rate

  rate = 1.d-13 * a * (t/1.d4)**b
end function rateFit

subroutine createRecombTable(table, tablefilename)
  type(RECOMBTABLE), intent(out) :: table
  character(len=*) :: tablefilename
  character(len=200) :: filename, datadirectory
  integer :: ia, ib, ne, ncut, i
  real :: e(1000)

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//tablefilename

  open(20,file=filename,status="old",form="formatted")
  read(20,*) table%nTemp, table%nrho

  allocate(table%rho(1:table%nRho))
  allocate(table%temp(1:table%ntemp))
  allocate(table%emissivity(1:table%ntemp,1:table%nrho))

  do ia=1,table%ntemp
     do ib=1,table%nrho
        read(20,'(1x,e10.3,5x,e10.3,13x,i2)') table%rho(ib),table%temp(ia),ncut
        ne=ncut*(ncut-1)/2
        read(20,'(8e10.3)') e(1:ne)
        table%emissivity(ia,ib) = SUM(e(1:ne-1))
     enddo
  enddo
  close(20)
end subroutine createRecombTable



  
subroutine createGammaTable(table, thisfilename)

! Ferland 1980 PASP 92 596

  type(GAMMATABLE), intent(out) :: table
  character(len=*) :: thisfilename
  character(len=200) :: dataDirectory, filename
  integer :: i

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//thisfilename

  open(40, file=filename, form="formatted", status="old")
  read(40,*) table%nTemp, table%nFreq

  allocate(table%freq(1:table%nFreq))
  allocate(table%temp(1:table%nTemp))
  allocate(table%gamma(table%nFreq, table%nTemp))

  table%temp = 0.d0
  read(40,*)  table%temp(1:table%nTemp)

  do i = 1, table%nFreq 
     read(40,*) table%freq(i), table%gamma(i,1:table%nTemp)
  enddo

  do i = 1, table%nFreq
     table%freq(i) = table%freq(i) * nuHydrogen
  enddo

  where (table%gamma(1:table%nFreq,1:table%nTemp) == 0.d0)
     table%gamma(1:table%nFreq,1:table%nTemp) = 1.d-10
  end where
  table%gamma(1:table%nFreq,1:table%nTemp) =log10(table%gamma(1:table%nFreq,1:table%nTemp))
  table%freq(1:table%nFreq) = log10(table%freq(1:table%nFreq))
  table%temp(1:table%nTemp) = log10(table%temp(1:table%nTemp))
  close(40)
end subroutine createGammaTable

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
  type(GRIDTYPE) :: grid
  TYPE(OCTAL) :: thisOctal
  integer :: subcell
  integer :: nFreq
  integer :: i, k, iIon, n1, n2
  real(double) :: freq(:), spectrum(:), dfreq(:)
  real(double) :: jnu
  real :: e, hxsec
  real(double), parameter :: statisticalWeight(3) = (/ 2.d0, 0.5d0, 2.d0 /)


  ! do Saha-Milne continua for H, HeI and HeII

  do k = 1 , 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (k == 1) iIon = 1
     if (k == 2 .and. .not. hOnly) iIon = 3
     if (k == 3 .and. .not. hOnly) iIon = 4


!     print *, "k = ", k
     

     call locate(freq, nfreq, grid%ion(iIon)%nuThresh, n1)
     n2 = nFreq
     if (iIon == 3) then
        call locate(freq, nfreq, grid%ion(iIon+1)%nuThresh, n2)
     endif
     do i = n1, n2


!     print *, "n1 = ", n1
!     print *, "n2 = ", n2
     
        
        e = freq(i) * hcgs* ergtoev

!       print *, "====DEBUG===="
!1       print *, " i ", i
!       print *, " e ", e
      ! print *, " grid%ion(iIon) ", grid%ion(iIon)
!       print *, " freq(i) ", freq(i)

       hxSec = returnxSec(grid%ion(iIon), freq(i), iFreq=i)

!        call phfit2(grid%ion(iIon)%z, grid%ion(iIon)%n, grid%ion(iIon)%outerShell , e , hxsec)

        jnu = tiny(jnu)

        if (hxSec > 0.) then
           if (thisOctal%temperature(subcell) > 100.) then
              jnu = statisticalWeight(k) * ((hcgs*freq(i)**3)/(cSpeed**2)) * &
                   ((hcgs**2) /(twoPi*mElectron*Kerg*thisOctal%temperature(subcell)))**(1.5d0) * &
                   dble(hxsec/1.d10) *  &
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
  real :: emissivity
!  real :: heII4686
!  integer :: ilow, iup
  integer,parameter :: nHeIILyman = 4
!  real(double) :: heIILyman(4)
!  real(double) :: freqheIILyman(4) = (/ 3.839530, 3.749542, 3.555121, 2.99963 /)



  ! HeI lines 

  call locate(heIrecombinationNe, 3, real(log10(thisOctal%ne(subcell))), i)
  fac = (log10(thisOctal%ne(subcell)) - heIrecombinationNe(i))/(heIrecombinationNe(i+1)-heIrecombinationNe(i))

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
!

!  HeII4686 = 10.d0**(-0.997d0*log10(thisOctal%temperature(subcell))+5.16d0)
!  HeII4686 = HeII4686*thisOctal%ne(subcell)*thisOctal%nh(subcell)*thisOctal%ionFrac(subcell,5)*grid%ion(4)%abundance
!  
!  ! calculate emission due to HeII recombination lines [e-25 ergs/s/cm^3]                                                     
!  do iup = 30, 3, -1
!     do ilow = 2, min0(16, iup-1)
!        emissivity= HeIIrecombinationLines(iup, ilow)*HeII4686*1.d-25
!
!
!        lambda = 227.838 / (1./real(ilow**2)  - 1./real(iup**2))!!!!!!!!!!!!!!!!!!!
!        lineFreq = cSpeed/(lambda * 1.d-8)
!
!     call locate(freq, nFreq, lineFreq, k)
!     spectrum(k) = spectrum(k) + emissivity
!     end do
!  end do
!

  ! He II Lyman series

!  heIILyman(1:4) = (/ 0.0334, 0.0682, 0.1849, 1. /)

!  ! calculate Lyman alpha first
!  HeIILyman(4) = 10.d0**(-0.792d0*log10(thisOctal%temperature(subcell))+6.01d0)
!  HeIILyman(4) = HeIILyman(4)*thisOctal%ne(subcell)*thisOctal%nh(subcell) * &
!       grid%ion(3)%abundance * thisOctal%ionFrac(subcell, 5) * 1.d-25
!
!  do i = 1, NHeIILyman-1
!     HeIILyman(i) = HeIILyman(i)*HeIILyman(4)
!  end do
!  do i = 1, nHeIILyman
!     lineFreq = freqHeIILyman(i) * nuHydrogen
!     call locate(freq, nFreq, lineFreq, k)
!     spectrum(k) = spectrum(k) + HeIILyman(i)
!  enddo



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
  real(double), allocatable :: kabsArray(:)


  allocate(kAbsArray(1:nlambda))

  call returnKappa(grid, thisOctal, subcell, kappaAbsArray=kAbsArray)


  do i = 1, nFreq
     thisLam = (cSpeed / freq(i)) * 1.e8
     if ((thisLam >= lamArray(1)).and.(thisLam <= lamArray(nlambda))) then
        call hunt(lamArray, nLambda, real(thisLam), iLam)
        spectrum(i) = spectrum(i) + bnu(freq(i), dble(thisOctal%temperature(subcell))) * &
             kAbsArray(iLam) *1.d-10* dFreq(i) * fourPi
     endif
  enddo

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


  recursive subroutine packvalues(thisOctal,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: nCrossings(:)
  real(double) :: photoIonCoeff(:,:)
  real(double) :: hHeating(:)
  real(double) :: heHeating(:)
  real(double) :: distanceGrid(:)

  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call packvalues(child, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          nCrossings(nIndex) = real(thisOctal%nCrossings(subcell))
          photoIonCoeff(nIndex, :) = thisOctal%photoIonCoeff(subcell, :)
          hHeating(nIndex) = thisOctal%hHeating(subcell)
          heHeating(nIndex) = thisOctal%heHeating(subcell)
          distanceGrid(nIndex) = thisOctal%distanceGrid(subcell)
       endif
    enddo
  end subroutine packvalues

  recursive subroutine unpackvalues(thisOctal,nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: nCrossings(:)
  real(double) :: photoIonCoeff(:,:)
  real(double) :: hHeating(:)
  real(double) :: heHeating(:)
  real(double) :: distanceGrid(:)

  integer :: nIndex
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unpackvalues(child, nIndex,nCrossings, photoIonCoeff, hHeating, HeHeating, distanceGrid)
                exit
             end if
          end do
       else
          nIndex = nIndex + 1
          thisOctal%nCrossings(subcell) = int(nCrossings(nIndex))
          thisOctal%photoIonCoeff(subcell, :) = photoIonCoeff(nIndex, :)
          thisOctal%hHeating(subcell) = hHeating(nIndex) 
          thisOctal%heHeating(subcell) = heHeating(nIndex) 
          thisOctal%distanceGrid(subcell) = distanceGrid(nIndex) 
       endif
    enddo
  end subroutine unpackvalues

  recursive subroutine  identifyUndersampled(thisOctal)
    use input_variables, only : minCrossings
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
  
    minCrossings = 100

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
       thisLambda = lamspec(i) + fac * (lamspec(i+1)-lamspec(i))
    else
       thisLambda = 1000.e4
  endif
    
  deallocate(freq, spectrum, lamSpec, dfreq)

  end function getRandomWavelengthPhotoion


  subroutine refineLambdaArray(lamArray, nLambda, grid)
    type(GRIDTYPE) :: grid
    real :: lamArray(:)
    integer :: nLambda
    integer :: i, j
    integer :: iup, ilow

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


  do iup = 15, 3, -1
     do ilow = 2, min0(8, iup-1)
        call insertBin(lamArray, nLambda, real(lambdaTrans(iup, ilow)*1.e8), 1.)
    enddo
 enddo
 call insertBin(lamArray, nLambda, 1215.67, 1.)


! refine for forbidden line transitions

    do i = 1 , grid%nIon
       do j = 1, grid%ion(i)%nTransitions
          call insertBin(lamArray, nLambda, &
               real(grid%ion(i)%transition(j)%lambda), 1.)
       enddo
    enddo

     
  end subroutine refineLambdaArray


  subroutine getNewMPIPhoton(position, direction, frequency, tPhoton, photonPacketWeight, iSignal)
    include 'mpif.h'
    integer :: ierr
    type(VECTOR) :: position, direction
    real(double) :: frequency, tPhoton, photonPacketWeight
    integer, parameter :: nTemp = 9
    real(double) :: tempstorage(nTemp)
    integer :: status(MPI_STATUS_SIZE)
    integer :: tag = 41
    integer, intent(out) :: iSignal



    iSignal = -1

    call MPI_RECV(tempStorage, nTemp, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, ierr)
    if (tempStorage(1) > 1.d30) then
       if (tempStorage(2) > 1.d30) then
          iSignal = 0
       else
          iSignal = 1
       endif
       goto 666
    else
       position%x = tempStorage(1)
       position%y = tempStorage(2)
       position%z = tempStorage(3)
       direction%x = tempStorage(4)
       direction%y = tempStorage(5)
       direction%z = tempStorage(6)
       frequency = tempStorage(7)
       tPhoton = tempStorage(8)
       photonPacketWeight = tempStorage(9)
    endif
666 continue
  end subroutine getNewMPIPhoton

  subroutine sendMPIPhoton(position, direction, frequency, tPhoton, photonPacketWeight, iThread)
    include 'mpif.h'
    integer :: ierr
    type(VECTOR) :: position, direction
    real(double) :: frequency, tPhoton, photonPacketWeight
    integer, parameter :: nTemp = 9
    real(double) :: tempstorage(nTemp)
    integer :: iThread
    integer :: tag = 41


    tempStorage(1) =        position%x  
    tempStorage(2) =        position%y   
    tempStorage(3) =        position%z   
    tempStorage(4) =        direction%x  
    tempStorage(5) =        direction%y  
    tempStorage(6) =        direction%z  
    tempStorage(7) =        frequency    
    tempStorage(8) =        tPhoton
    tempStorage(9) =        photonPacketWeight
    if (iThread == myRankGlobal) then
       write(*,*) "sending to self bug ", ithread
       stop
    endif
    call MPI_SEND(tempStorage, nTemp, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD,  ierr)
  end subroutine sendMPIPhoton
       


  recursive subroutine calculateEnergyFromTemperature(thisOctal)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: eThermal
    real(double) :: mu

    do subcell = 1, thisOctal%maxChildren
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
          eThermal = 1.5d0 * thisOctal%temperature(subcell)*kerg/(mu*mHydrogen)
          thisOctal%energy(subcell) =  eThermal
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)

       endif
    enddo
  end subroutine calculateEnergyFromTemperature

  subroutine createImageSplitGrid(grid, nSource, source, observerDirection, imageFilename, &
       lambdaImage, outputimageType, nPixels)
    use photoion_mod, only : addForbiddenEmissionLine, addRecombinationEmissionLine
    use input_variables, only : freeFreeImage, nPhotons
    real :: lambdaImage
    include 'mpif.h'
    type(GRIDTYPE) :: grid
    integer :: npixels
    character(len=*) :: imageFilename, outputImageType
    integer :: nSource
    type(SOURCETYPE) :: source(:), thisSource
    type(PHOTON) :: thisPhoton, observerPhoton
    type(OCTAL), pointer :: thisOctal
    real(double) :: totalFlux
    real(double), allocatable :: totalFluxArray(:), tempTotalFlux(:)
    logical :: directFromSource
    integer :: subcell
    integer :: iPhoton
    integer :: iSource
    integer :: iThread
    type(VECTOR) ::  rHat, observerDirection
    logical :: endLoop, addToiMage
    integer :: newthread
    type(IMAGETYPE) :: thisImage
    logical :: escaped, absorbed, crossedBoundary, photonsStillProcessing, stillSCattering
    real(double) :: totalEmission
    integer :: iLambdaPhoton, nInf
    real(double) :: lCore, probsource, r
    real(double), allocatable :: threadProbArray(:)
    integer :: np(10)
    integer :: nDone
    integer :: tag = 41
    integer, allocatable :: nDoneArray(:)
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    integer :: isignal
    real(double) :: powerPerPhoton, photonPacketWeight

    call randomNumberGenerator(randomSeed=.true.)

    absorbed = .false.
    escaped = .false.

    allocate(nDoneArray(1:nThreadsGlobal-1))
    allocate(totalFluxArray(1:nThreadsGlobal-1))
    allocate(tempTotalFlux(1:nThreadsGlobal-1))

    totalFluxArray = 0.d0



    call zeroEtaCont(grid%octreeRoot)
    
    call quickSublimate(grid%octreeRoot) ! do dust sublimation
    
    call torus_mpi_barrier


    thisImage = initImage(npixels, npixels, real(2.*grid%octreeRoot%subcellSize), &
         real(2.*grid%octreeRoot%subcellSize), 0., 0.)

    allocate(threadProbArray(1:nThreadsGlobal-1))


    select case (outputimageType)
       case("freefree")
          freefreeImage = .true.
          ilambdaPhoton = grid%nLambda
          call  addRadioContinuumEmissivity(grid%octreeRoot)
          lcore = tiny(lcore)

       case("forbidden")

          call locate(grid%lamArray, grid%nlambda, real(lambdaImage), ilambdaPhoton)
          call calcContinuumEmissivityLucyMono(grid, grid%octreeRoot, nlambda, grid%lamArray, lambdaImage, iLambdaPhoton)
          call addForbiddenEmissionLine(grid, 1.d0, dble(lambdaImage))

          if (nSource > 0) then              
             lCore = sumSourceLuminosityMonochromatic(source, nsource, dble(grid%lamArray(iLambdaPhoton)))
          else
             lcore = tiny(lcore)
          endif

       case("recombination")

          call locate(grid%lamArray, grid%nlambda, real(lambdaImage), ilambdaPhoton)
          call calcContinuumEmissivityLucyMono(grid, grid%octreeRoot, nlambda, grid%lamArray, lambdaImage, iLambdaPhoton)
          call addRecombinationEmissionLine(grid, 1.d0, dble(lambdaImage))

          if (nSource > 0) then              
             lCore = sumSourceLuminosityMonochromatic(source, nsource, dble(grid%lamArray(iLambdaPhoton)))
          else
             lcore = tiny(lcore)
          endif

       case("dustonly")

          call locate(grid%lamArray, grid%nlambda, real(lambdaImage), ilambdaPhoton)
          call calcContinuumEmissivityLucyMono(grid, grid%octreeRoot, nlambda, grid%lamArray, lambdaImage,iLambdaPhoton)

          if (nSource > 0) then              
             lCore = sumSourceLuminosityMonochromatic(source, nsource, dble(grid%lamArray(iLambdaPhoton)))
          else
             lcore = tiny(lcore)
          endif
       case DEFAULT
          call writeFatal("Imagetype "//trim(outputimageType)//" not recognised")
          stop

    end select


    call computeProbDistAMRMpi(grid, totalEmission, threadProbArray)

    if (myrankglobal == 0) write(*,*) "prob array ", threadProbArray(1:nThreadsGlobal-1)
    totalEmission = totalEmission * 1.d30



    probSource = lCore / (lCore + totalEmission)

    if (myRankGlobal == 0) then
       write(*,*) "Probability of photon from sources: ", probSource
    endif

    Ninf = 0

    powerPerPhoton = (lCore + totalEmission) / dble(nPhotons)
    if (Writeoutput) write(*,*) "power per photon ",powerperphoton

    if (myRankGlobal == 0) then
       np = 0
       mainloop: do iPhoton = 1, nPhotons

          thisPhoton%stokes = STOKESVECTOR(1.d0, 0.d0, 0.d0, 0.d0)
          thisPhoton%iLam = iLambdaPhoton
          thisPhoton%lambda = grid%lamArray(iLambdaPhoton)
          thisPhoton%observerPhoton = .false.
          call randomNumberGenerator(getDouble=r)


          if (r < probSource) then
             call randomSource(source, nSource, iSource, photonPacketWeight)
             thisSource = source(iSource)
             call getPhotonPositionDirection(thisSource, thisPhoton%position, thisPhoton%direction, rHat,grid)


!             if (thisSource%outsideGrid) then
!                thisPhoton%stokes%i = thisPhoton%stokes%i * (2.d0*grid%octreeRoot%subcellSize*1.d10)**2 * &
!                   (thisSource%radius*1.d10)**2 / (thisSource%distance**2)
!             endif

             call findSubcellTD(thisPhoton%position, grid%octreeRoot,thisOctal, subcell)
             iThread = thisOctal%mpiThread(subcell)
             call sendPhoton(thisPhoton, iThread, endloop = .false.) 
             directFromSource = .true.
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
          do iThread = 1, nThreadsGlobal-1
             call sendPhoton(thisPhoton, iThread, endLoop = .false., report=.true.)
             call MPI_RECV(nDoneArray(iThread), 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
          enddo
          nDone = SUM(nDoneArray(1:nThreadsGlobal-1))
          write(*,*) myrankglobal, " thinks ", nDone, " photons have completed "
          if (nDone == nPhotons) photonsStillProcessing = .false.
       enddo
       do iThread = 1, nThreadsGlobal-1
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
             call MPI_SEND(nDone, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD,  ierr)
             goto 777
          end if


          stillScattering = .true.
          ninf = ninf + 1
          do while ((.not.endLoop).and.stillScattering)

             if (thisPhoton%observerPhoton) then
                newThread = -1
                if (.not.freeFreeImage) then
                   call propagateObserverPhoton(grid, thisPhoton, addToImage, newThread)
                else
                   addtoImage = .true.
                endif
                if (addToImage) then
                   call addPhotonToImageLocal(observerDirection, thisImage, thisPhoton, totalFluxArray(myRankGlobal))
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
                      call addPhotonToImageLocal(observerDirection, thisImage, observerPhoton, totalFluxArray(myRankGlobal))
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

     call MPI_ALLREDUCE(totalFluxArray, tempTotalFlux, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     totalFluxArray = tempTotalFlux
     totalFlux = SUM(totalFluxArray(1:nThreadsGlobal-1))
#ifdef USECFITSIO
    if (myrankGlobal == 0) then
       call writeFitsImage(thisimage, imageFilename, 1.d0, "intensity")
    endif
#else
    call writeInfo("FITS not enabled, not writing "//trim(imageFilename),FORINFO)
#endif
    call freeImage(thisImage)
  end subroutine createImageSplitGrid


  subroutine sendPhoton(thisPhoton, iThread, endLoop, report, getPosition)
    include 'mpif.h'
    type(PHOTON) :: thisPhoton
    integer :: iThread
    logical :: endLoop
    logical, optional :: report, getPosition
    integer, parameter  :: nTemp = 15
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

    if (iThread == myRankGlobal) then
       write(*,*) "sending to self bug ", ithread
       stop
    endif

    call MPI_SEND(temp, nTemp, MPI_DOUBLE_PRECISION, iThread, tag, MPI_COMM_WORLD,  ierr)
    
  end subroutine sendPhoton

  subroutine receivePhoton(thisPhoton, iSignal)
    include 'mpif.h'
    type(PHOTON) :: thisPhoton
    integer, parameter :: nTemp = 15
    real(double) :: temp(15)
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    integer :: tag = 42
    integer, intent(out) :: iSignal

    call MPI_RECV(temp, nTemp, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, ierr)
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

    thisPhoton%lambda =     temp(11) 

    thisPhoton%tau =     temp(13) 
    thisPhoton%iLam =     nint(temp(14) )

    if (temp(15) > 0.d0) then
       thisPhoton%observerPhoton = .true.
    else
       thisPhoton%observerPhoton = .false.
    endif
    
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
       thisPhoton%tau = thisPhoton%tau + tval * kappaExt
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


  subroutine collateImages(thisImage)
    include 'mpif.h'
    type(IMAGETYPE) :: thisImage
    real, allocatable :: tempRealArray(:), tempRealArray2(:)
    real(double), allocatable :: tempDoubleArray(:), tempDoubleArray2(:)
    integer :: ierr

     allocate(tempRealArray(SIZE(thisImage%pixel)))
     allocate(tempRealArray2(SIZE(thisImage%pixel)))
     allocate(tempDoubleArray(SIZE(thisImage%pixel)))
     allocate(tempDoubleArray2(SIZE(thisImage%pixel)))
     tempRealArray = 0.0
     tempRealArray2 = 0.0
     tempDoubleArray = 0.0_db
     tempDoubleArray2 = 0.0_db

     if (myrankGlobal == 1) write(*,*) "Collating images..."
     tempDoubleArray = reshape(thisImage%pixel%i,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%i = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%i))

     tempDoubleArray = reshape(thisImage%pixel%q,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%q = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%q))

     tempDoubleArray = reshape(thisImage%pixel%u,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%u = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%u))

     tempDoubleArray = reshape(thisImage%pixel%v,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%v = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%v))


     tempRealArray = reshape(thisImage%vel,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%vel = reshape(tempRealArray2,SHAPE(thisImage%vel))

     tempRealArray = reshape(thisImage%totWeight,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%totWeight = reshape(tempRealArray2,SHAPE(thisImage%totWeight))

     if (myrankGlobal == 1) write(*,*) "Done."
     deallocate(tempRealArray)
     deallocate(tempRealArray2)
     deallocate(tempDoubleArray)
     deallocate(tempDoubleArray2)

  end subroutine collateImages

   subroutine addPhotonToImageLocal(observerDirection, thisImage, thisPhoton, totalFlux)
     
     type(IMAGETYPE), intent(inout) :: thisImage
     type(PHOTON) :: thisPhoton
     type(VECTOR) :: observerDirection,  xProj, yProj
     real :: xDist, yDist
     integer :: xPix, yPix
     real(double) :: totalFlux

     type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.d0)

     xPix = 0; yPix = 0


     xProj =  zAxis .cross. observerDirection
     call normalize(xProj)
     yProj = observerDirection .cross. xProj
     call normalize(yProj)
     xDist = (thisPhoton%position) .dot. xProj
     yDist = (thisPhoton%position) .dot. yProj
           

     call pixelLocate(thisImage, xDist, yDist, xPix, yPix)

     if ((xPix >= 1) .and. &
          (yPix >= 1) .and. &
          (xPix <= thisImage%nx) .and. &
          (yPix <= thisImage%ny)) then

              
        thisImage%pixel(xPix, yPix) = thisImage%pixel(xPix, yPix)  &
             + thisPhoton%stokes * oneOnFourPi * exp(-thisPhoton%tau)

!        write(*,*) xpix,ypix, thisImage%pixel(xpix,ypix)%i, &
!             thisPhoton%stokes%i,  thisPhoton%stokes * oneOnFourPi * exp(-thisPhoton%tau)
     endif
     totalFlux = totalFlux + thisPhoton%stokes%i * oneOnFourPi * exp(-thisPhoton%tau)
           
   end subroutine addPhotonToImageLocal

   subroutine  computeProbDistAMRMpi(grid, totalEmission, threadProbArray)
     include 'mpif.h'
     type(GRIDTYPE) :: grid
     real(double) :: totalEmission, totalProb, biasCorrection
     real(double) :: threadProbArray(:)
     real(double), allocatable :: totalEmissionArray(:), totalProbArray(:), tArray(:)
     integer :: ierr, i

     totalEmission = 0.d0

     allocate(totalEmissionArray(1:nThreadsGlobal), totalProbArray(1:nThreadsGlobal), &
          tarray(1:nThreadsGlobal))
     totalEmissionArray = 0.d0
     totalProbArray = 0.d0

     call computeProbDist2AMRMpi(grid%octreeRoot,totalEmissionArray(myRankGlobal+1), totalProbArray(myRankGlobal+1))
     write(*,*) myrankGlobal, " total emission  ", totalEmissionArray(myrankGlobal+1)

     tArray = 0.d0
     call MPI_ALLREDUCE(totalEmissionArray, tArray, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     totalEmissionArray = tArray

     tArray = 0.d0
     call MPI_ALLREDUCE(totalProbArray, tArray, nThreadsGlobal, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     totalProb = SUM(tArray(2:nThreadsGlobal))
     
     threadProbArray = tArray(2:nThreadsGlobal)
     do i = 2, SIZE(threadProbArray)
        threadProbArray(i) = threadProbArray(i) + threadProbArray(i-1)
     enddo
     if (threadProbArray(SIZE(threadPRobArray)) /= 0.d0) then
        threadProbArray = threadProbArray / threadProbArray(SIZE(threadProbArray))
     else
        threadProbArray = 1.d0
     endif


     totalEmission = SUM(totalEmissionArray(2:nThreadsGlobal))

     if (totalProb /= 0.0) then
        biasCorrection = totalEmission / totalProb
     else
        biasCorrection = 1.0
     end if
     
     call computeProbDist3AMRMpi(grid%octreeRoot, biasCorrection, totalProbArray(myRankGlobal+1))

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

  recursive subroutine zeroEtaCont(thisOctal)

    implicit none
    type(octal), pointer              :: thisOctal
    integer :: subcell
    type(octal), pointer  :: child 

    if (thisOctal%nChildren > 0) then
       ! call this subroutine recursively on each of its children
       do subcell = 1, thisOctal%nChildren, 1 
          child => thisOctal%child(subcell)
          call zeroEtaCont(child)
       end do 
    end if
    thisOctal%biasCont3D = 1.d0
    thisOctal%etaCont = 0.d0
  end subroutine zeroEtaCont

  recursive subroutine addRadioContinuumEmissivity(thisOctal)

    use stateq_mod, only : alpkk
    type(octal), pointer                 :: thisOctal
    type(octal), pointer  :: child 
    integer               :: subcell
    integer               :: i
    real(double) :: eta, freq
    
    
    do subcell = 1, thisOctal%maxChildren, 1

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call addRadioContinuumEmissivity(child)
                exit
             end if
          end do
            
       else

          freq = cspeed / (20.d0) ! 20 cm radio free-free
          eta =  thisOctal%Ne(subcell)**2 * &
               alpkk(freq,real(thisOctal%temperature(subcell),kind=db))* &
               exp(-(hcgs*freq)/(kerg*thisOctal%temperature(subcell)))
          
          eta=eta*real((2.0*dble(hcgs)*dble(freq)**3)/(dble(cspeed)**2))
             
          thisOctal%etaCont(subcell) = eta * 1.d10
          thisOctal%biasCont3d(subcell) = 1.d0
       end if

    end do

  end subroutine addRadioContinuumEmissivity


#else

contains

! Dummy subroutines for non-MPI case
  SUBROUTINE resizePhotoionCoeff(thisOctal,grid)

    use grid_mod
    implicit none
    type(GRIDTYPE) :: grid
    TYPE(OCTAL), POINTER  :: thisOctal 
    integer :: i
    i = thisOctal%nChildren
    i = grid%octreeRoot%nChildren
  END SUBROUTINE resizePhotoionCoeff

!  subroutine photoIonizationloopAMR(grid, source, nSource, nLambda, lamArray, readlucy, writelucy, &
!       lucyfileout, lucyfilein, maxIter, tLimit, sublimate)
!
!    use grid_mod
!    use source_mod
!
!    implicit none
!
!    type(GRIDTYPE) :: grid
!    logical, optional :: sublimate
!    character(len=*) :: lucyfileout, lucyfilein    
!    integer :: nSource
!    type(SOURCETYPE) :: source(:)
!    integer :: nlambda
!    real :: lamArray(:)
!    logical :: readLucy, writeLucy
!    integer :: maxIter
!    real(double) :: tlimit
!
!  end subroutine photoIonizationloopAMR
!
!  subroutine ionizeGrid(thisOctal)
!    use octal_mod
!    type(octal) :: thisOctal
!  end subroutine ionizeGrid
!
!
!  subroutine createImageSplitGrid(grid, nSource, source, observerDirection, imageFilename, lambdaImage, outputtype, npix)
!    use gridtype_mod
!    use vector_mod
!    use source_mod
!    type(GRIDTYPE) :: grid
!    integer :: nSource
!    type(SOURCETYPE) :: source(:)
!    type(VECTOR) :: observerDirection
!    real :: lambdaImage
!    character(len=*) :: imageFilename, outputType
!    integer :: npix
!  end subroutine createImageSplitGrid

#endif    
end module photoionAMR_mod

