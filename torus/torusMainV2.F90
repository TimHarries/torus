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
  use input_variables         ! variables filled by inputs subroutine
  use constants_mod
  use messages_mod
  use mpi_global_mod
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
  type(GRIDTYPE) :: grid
#ifdef MPI
   include 'mpif.h'  
#endif
#ifdef MPI
  ! For MPI implementations =====================================================
  integer ::   ierr           ! error flag
#endif

#ifdef MPI
  ! FOR MPI IMPLEMENTATION=======================================================
  !  initialize the system for running MPI
  call MPI_INIT(ierr) 

  ! Set up amrCOMMUNICATOR and global mpi groups
  call setupAMRCOMMUNICATOR

  !===============================================================================
#endif

  globalMemoryChecking = .false.
#ifdef MEMCHECK
  globalMemoryChecking = .true.
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


  call writeTorusBanner()

  call setVersion("V2.0")
  grid%version = torusVersion
  verbosityLevel = 5
  call writeBanner("TORUS ("//trim(torusVersion)//") model","-",IMPORTANT)
  call report_parallel_type

  ! For time statistics
  if (doTuning) call tune(6, "Torus Main") ! start a stopwatch  

  ! set up a random seed
  if (torusOpenMP) call run_random_test_suite

!  call testSuiteRandom()  
  call inputs()

  call setupMicrophysics(grid)

  call  setupamrgrid(grid)


!  call checkAMRgrid(grid, .false.)

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
#ifdef MPI
  call torus_mpi_barrier
  call freeAMRCOMMUNICATOR
  call MPI_FINALIZE(ierr)
#endif

  call writeBanner("Torus completed","o",TRIVIAL)

end program torus
