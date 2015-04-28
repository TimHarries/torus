module mpi_global_mod

  implicit none

  integer :: amrCOMMUNICATOR
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

contains

! Report information about how Torus was built
  subroutine report_parallel_type
    use messages_mod
    implicit none
    
    call writeBanner("Build options","%",TRIVIAL)


    if (TorusHybrid) then 
       call  writeInfo('Using hybrid MPI/OpenMP', TRIVIAL)
    elseif (TorusMpi) then 
       call  writeInfo('Using MPI', TRIVIAL)
    elseif (TorusOpenmp) then 
       call  writeInfo('Using OpenMP', TRIVIAL)
    elseif (TorusSerial) then
       call  writeInfo('No parallelism is in use', TRIVIAL)
    end if

    call  writeInfo('', TRIVIAL)
    call writeInfo("Options used to build Torus were:", TRIVIAL)
#ifdef HYDRO
     call writeInfo("hydro=yes", TRIVIAL)
#else
     call writeInfo("hydro=no", TRIVIAL)
#endif
#ifdef PHOTOION
     call writeInfo("photoion=yes", TRIVIAL)
#else
     call writeInfo("photoion=no", TRIVIAL)
#endif
#ifdef MOLECULAR
     call writeInfo("molecular=yes", TRIVIAL)
#else
     call writeInfo("molecular=no", TRIVIAL)
#endif
#ifdef CMFATOM
     call writeInfo("atomic=yes", TRIVIAL)
#else
     call writeInfo("atomic=no", TRIVIAL)
#endif
#ifdef SPH
     call writeInfo("sph=yes", TRIVIAL)
#else
     call writeInfo("sph=no", TRIVIAL)
#endif
#ifdef USECFITSIO
     call writeInfo("cfitsio=yes", TRIVIAL)
#else
     call writeInfo("cfitsio=no", TRIVIAL)
#endif
    call  writeInfo('', TRIVIAL)

  end subroutine report_parallel_type

end module mpi_global_mod
