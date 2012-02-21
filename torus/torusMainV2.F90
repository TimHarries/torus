!
!            T   O   R   U  S
!
! torus is a general purpose monte-carlo radiative transfer code used to
! produce polarized spectra of hot and cool stellar winds
!

! written by tjh
! broken by chris

! v1.0 on 16/09/99

! raman scattering stuff added 3/3/2000

! OMP parallelization calls added 1/7/2001

! TJH: 19/7/02  Jorick's spotty star stuff added
! NHS: 02/9/02  adaptive mesh code merged
! RK : 12/02/03 Parallelized major loops in lucyStateEquilibriumAMR, stateqAMR, torusMain.
! RK : 07/07/05 Added the observed flux solver in oberver's frame. 
! TJH: 10/09/08 pgplot calls removed...
program torus
  use torus_version_mod
  use inputs_mod         ! variables filled by inputs subroutine
  use dimensionality_mod
  use constants_mod
  use messages_mod
  use mpi_global_mod
  use zlib_mod
  use utils_mod
  use inputs_mod
  use timing
  use grid_mod
  use gridio_mod
  use setupamr_mod
  use physics_mod
  use outputs_mod
  use source_mod
  use random_mod
  use memory_mod
  use zlib_mod
#ifdef MPI
  use mpi
#endif

  implicit none

  character(len=80) :: message
  type(GRIDTYPE) :: grid

#ifdef MPI
  ! For MPI implementations =====================================================
  integer ::   ierr           ! error flag
#endif
#ifdef _OPENMP
  call omp_set_dynamic(.false.)
#endif


#ifdef MPI
  ! FOR MPI IMPLEMENTATION=======================================================
  !  initialize the system for running MPI
  call MPI_INIT(ierr) 
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRankGlobal, ierr)
  !===============================================================================
#endif

  writeoutput    = .true.
  doTuning       = .true.
  outputwarnings = .true.
  outputinfo     = .true.
  myRankIsZero   = .true.

  if (TorusMpi) then 
     if (myRankGlobal/=1) writeoutput  = .false.
     if (myRankGlobal/=1) doTuning     = .false.
     if (myRankGlobal/=0) myRankIsZero = .false.
  end if

  call initializeCodeUnits()

  call writeTorusBanner()

  call setVersion("V2.0")
  grid%version = torusVersion
  verbosityLevel = 5
  call writeBanner("TORUS "//trim(torusVersion)//" model","-",IMPORTANT)
  call report_parallel_type

  ! For time statistics
  if (doTuning) call tune(6, "Torus Main") ! start a stopwatch  


  !  call testSuiteRandom()  
  call inputs()

#ifdef MPI
  ! Set up amrCOMMUNICATOR and global mpi groups
  call setupAMRCOMMUNICATOR
#endif

  ! set up a random seed
  call run_random_test_suite

  smallestCellSize = amrGridSize / dble(2**maxDepthAMR)

  do iModel = nModelStart, nModelEnd
     if (multimodels) then
        write(message,'(a,i6.6)') "Performing calculation for model number ", iModel
        call writeBanner(message,"-",TRIVIAL)
     endif

#ifdef USEZLIB

#else
     useBinaryXMLVTKfiles = .false.
#endif

     call setupMicrophysics(grid)



     call  setupamrgrid(grid)


     !  call checkAMRgrid(grid, .false.)

     !  call writeAMRgrid("test.dat",.false.,grid)
     !  call writeVtkFile(grid, "mpi.vtk",  valueTypeString=(/"mpithread"/))
     !  call torus_mpi_barrier
     !  goto 666

     !  call writeVtkFile(grid, "rho.vtk")

     call setupGlobalSources(grid)


     call writeBanner("Run-time messages","+",TRIVIAL)

     call randomNumberGenerator(randomSeed=.true.)

     call doPhysics(grid)

     call doOutputs(grid)

     if (doTuning) call tune(6, "Torus Main") ! stop a stopwatch  


     call torus_mpi_barrier
     call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
     call freeGrid(grid)
     call freeGlobalSourceArray()

  enddo

#ifdef MPI
  call torus_mpi_barrier
  if (hydrodynamics) call freeAMRCOMMUNICATOR
  call MPI_FINALIZE(ierr)
#endif
  !666 continue
  call writeBanner("Torus completed","o",TRIVIAL)

end program torus
