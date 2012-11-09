module mpi_global_mod

  implicit none

  integer :: amrCOMMUNICATOR
  integer :: localWorldCommunicator
  integer :: myRankGlobal
  integer :: myRankWorldGlobal
  integer :: nThreadsGlobal
  integer :: nHydroThreadsGlobal
  integer :: nHydroSetsGlobal
  integer :: myHydroSetGlobal
  integer, allocatable :: amrParallelCommunicator(:)
  integer, allocatable :: hydroCommunicator(:)
  integer :: zeroThreadCommunicator

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

contains

  subroutine report_parallel_type
    use messages_mod
    implicit none
    
    if (TorusHybrid) then 
       call  writeInfo('Using hybrid MPI/OpenMP', FORINFO)
    elseif (TorusMpi) then 
       call  writeInfo('Using MPI', FORINFO)
    elseif (TorusOpenmp) then 
       call  writeInfo('Using OpenMP', FORINFO)
    elseif (TorusSerial) then
       call  writeInfo('No parallelism is in use', FORINFO)
    end if

  end subroutine report_parallel_type

end module mpi_global_mod
