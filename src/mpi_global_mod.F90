module mpi_global_mod

  implicit none

  integer :: amrCOMMUNICATOR
  integer :: myFitterSetGlobal
  integer :: localWorldCommunicator
  integer :: myRankGlobal
  integer :: myRankWorldGlobal
  integer :: nThreadsGlobal
  integer :: nHydroThreadsGlobal
  integer :: nLoadBalancingThreadsGlobal
  integer :: nHydroSetsGlobal
  integer :: myHydroSetGlobal
  integer :: zeroPlusAMRCommunicator
  integer, allocatable :: amrParallelCommunicator(:)
  integer, allocatable :: loadBalanceCommunicator(:)
  integer, allocatable :: hydroCommunicator(:)
  integer :: zeroThreadCommunicator
  integer :: allDomainsCommunicator
  logical :: loadBalancingThreadGlobal
  integer :: copyOfThread

! Set logical variables to say how Torus is parallelised
! For hybrid MPI/OpenMP configurations all three are set to true
#ifdef MPI
  logical, parameter :: TorusMpi=.true.
#else
  logical, parameter :: TorusMpi=.false.
#endif

#ifdef _OPENMP
  logical, parameter :: TorusOpenmp=.true.
#else
  logical, parameter :: TorusOpenmp=.false.
#endif

  logical, parameter :: TorusHybrid= TorusOpenmp .and. TorusMpi
  logical, parameter :: TorusSerial= (.not. TorusOpenmp) .and. (.not. TorusMpi)

end module mpi_global_mod
