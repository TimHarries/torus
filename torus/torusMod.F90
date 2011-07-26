!
!            T   O   R   U  S
!
! torus is a general purpose monte-carlo radiative transfer code used to
! produce polarized spectra of hot and cool stellar winds
!

module torus_mod

  implicit none 

  public :: torus

  private

contains

subroutine torus(b_idim,  b_npart,       b_nptmass,  b_num_gas,   &  
                   b_xyzmh, b_rho,         b_iphase,              &
                   b_udist, b_umass,       b_utime,               &
                   b_time,  b_temp,        b_totalgasmass,        &
                   file_tag              )

  use torus_version_mod
  use inputs_mod         ! variables filled by inputs subroutine
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
  use sph_data_class 
  use phasematrix_mod, only: resetNewDirectionMie
  use amr_mod, only: returnKappa
  type(GRIDTYPE) :: grid
  type(VECTOR) :: dummy

! Variables used when linking to sph code
  integer, intent(in)   :: b_idim, b_npart, b_nptmass
  integer*1, intent(in) :: b_iphase(b_idim)
  real*8, intent(in)    :: b_xyzmh(5,b_idim)
  real*4, intent(in)    :: b_rho(b_idim)
  real*8, intent(in)    :: b_udist, b_umass, b_utime, b_time
  integer, intent(in)   :: b_num_gas           ! Number of gas particles
  real*8, intent(inout) :: b_temp(b_idim)   ! Temperature of gas particles
  real*8, intent(in)    :: b_totalgasmass      ! Total gas mass for this MPI process
  character(len=11), intent(in) :: file_tag

#ifdef MPI
  ! FOR MPI IMPLEMENTATION=======================================================

  ! Set up amrCOMMUNICATOR and global mpi groups
  call setupAMRCOMMUNICATOR

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
  call writeBanner("TORUS ("//trim(torusVersion)//") model","-",IMPORTANT)
  call report_parallel_type
  call writeInfo("Using SPH step "//file_tag, FORINFO)

  ! For time statistics
  if (doTuning) call tune(6, "Torus Main") ! start a stopwatch  

  ! set up a random seed
  call run_random_test_suite

!  call testSuiteRandom()  
  call inputs()

  call setupMicrophysics(grid)

  call init_sphtorus(b_udist, b_umass, b_utime, b_num_gas, b_time, b_nptmass, &
       b_npart, b_idim, b_iphase, b_xyzmh, b_rho, b_temp, b_totalgasmass)
  call gather_sph_data

  call setupamrgrid(grid)

  call setupGlobalSources(grid)

  call writeBanner("Run-time messages","+",TRIVIAL)

  call randomNumberGenerator(randomSeed=.true.)

  call doPhysics(grid)

  call doOutputs(grid)

  call update_sph_temperature (b_idim, b_npart, b_iphase, b_xyzmh, grid, b_temp)

  if (doTuning) call tune(6, "Torus Main") ! stop a stopwatch  


  call torus_mpi_barrier
  call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
  call freeGrid(grid)
  call freeGlobalSourceArray()
  dummy = clusterparameter(VECTOR(0.d0,0.d0,0.d0),grid%octreeroot, subcell = 1, isdone = .true.)
  call kill() ! Free SPH data type
  call resetNewDirectionMie()
  call returnKappa(grid, grid%octreeRoot, 1, reset_Kappa=.true.)
#ifdef MPI
  call torus_mpi_barrier
  call freeAMRCOMMUNICATOR
#endif

  call writeBanner("Torus completed","o",TRIVIAL)

end subroutine torus

!-----------------------------------------------------------------------------------------------------------------------

! Name:    update_sph_temperature
! Purpose: Update the temperatures of the SPH particles using the torus temperature field
!          Each process works on its own subset of the total number of particles.
! Author:  D. Acreman, November 2007

  subroutine update_sph_temperature (b_idim, b_npart, b_iphase, b_xyzmh, grid, b_temp)

    USE vector_mod, only:     vector
    USE amr_mod, only:        amrGridValues
    USE gridtype_mod, only:   gridType
    USE sph_data_class, only: sphData, get_udist
    USE messages_mod
    USE kind_mod
#ifdef MPI
    USE mpi
#endif

    implicit none

#ifdef MPI
! MPI specific variables
    real    :: mpi_max_deltaT, mpi_sum_deltaT
    integer :: ierr, mpi_iiigas
#endif

! Arguments 
    integer, intent(in)   :: b_idim, b_npart
    integer*1, intent(in) :: b_iphase(b_idim)
    real*8, intent(in)    :: b_xyzmh(5,b_idim)
    real*8, intent(inout) :: b_temp(b_idim)
    type(GRIDTYPE), intent(in) :: grid

! Local variables
    character(len=80) :: message
    integer      :: iiigas, i
    real*8       :: xgas, ygas, zgas
    real(double) :: sphDistFac
    real         :: tgas
    real         :: deltaT, sum_deltaT, mean_deltaT, max_deltaT
    type(VECTOR) :: positionVec

! Begin executable statements

! 1. Calculate conversion from sph distance units to torus distance units
       sphDistfac  = get_udist() ! [cm]
       sphDistfac = sphDistfac / 1.e10  ! to torus units

! 2. Update particle temperatures and calculate some statistics 
       iiigas = 0
       sum_deltaT = 0.0
       max_deltaT = 0.0
       do i=1, b_npart
          if (b_iphase(i) == 0) then
             iiigas = iiigas + 1 
! 2.1 Determine position of gas particle on the amr grid and extract the temperature value
             xgas = b_xyzmh(1,i) * sphDistFac
             ygas = b_xyzmh(2,i) * sphDistFac
             zgas = b_xyzmh(3,i) * sphDistFac
             positionVec = VECTOR(xgas,ygas,zgas)
             call amrGridValues(grid%octreeRoot, positionVec, temperature=tgas, grid=grid)
! 2.2 Calculate statistics of temperature change
             deltaT     = tgas - b_temp(i)
             sum_deltaT = sum_deltaT + deltaT
             max_deltaT = MAX(max_deltaT, deltaT)
! 2.3 Update the gas particle temperature to pass back to sph code
             b_temp(i) = tgas
          endif
       enddo

! 3. Calculate mean temperature change and perform MPI communication if required
!    MPI processes which did not do the SPH step will not have any particles to update
#ifdef MPI
! Calculate global values
       call MPI_ALLREDUCE( iiigas,     mpi_iiigas,     1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
       call MPI_ALLREDUCE( max_deltaT, mpi_max_deltaT, 1, MPI_REAL,    MPI_MAX, MPI_COMM_WORLD, ierr )
       call MPI_ALLREDUCE( sum_deltaT, mpi_sum_deltaT, 1, MPI_REAL,    MPI_SUM, MPI_COMM_WORLD, ierr )
! Update values to be output
       iiigas      = mpi_iiigas
       max_deltaT  = mpi_max_deltaT
       mean_deltaT = mpi_sum_deltaT / real(mpi_iiigas)
#else
       mean_deltaT = sum_deltaT / real(iiigas)
#endif

! 4. Write out the temperature change statistics
       write(message, *) "Number of particles= ", iiigas
       call writeInfo(message, FORINFO)
       write(message, *) "Maximum temperature change= ", max_deltaT
       call writeInfo(message, FORINFO)
       write(message, *) "Mean temperature change=  ", mean_deltaT
       call writeInfo(message, FORINFO)

  end subroutine update_sph_temperature

!-------------------------------------------------------------------------------

#ifdef MPI 

! Name:    gather_sph_data
! Purpose: Perform the mpi communication required to ensure all processes have
!          access to the full list of sph particle data
! Author:  D. Acreman October 2007

  subroutine gather_sph_data

    USE kind_mod
    USE sph_data_class, only: sphData, npart
    USE mpi

    implicit none

    integer :: my_rank, n_proc, ierr, i 
    integer :: npart_all, nptmass_all
    real(double)    :: totalgasmass_all
    integer, allocatable :: npart_arr(:),   nptmass_arr(:)
    integer, allocatable :: npart_displ(:), nptmass_displ(:)
    real(double), allocatable :: xn_tmp(:), yn_tmp(:), zn_tmp(:), rhon_tmp(:), temperature_tmp(:)
    real(double), allocatable :: gasmass_tmp(:), hn_tmp(:)
    real(double), allocatable :: x_tmp(:),  y_tmp(:),  z_tmp(:),  ptmass_tmp(:)
    character(len=3) :: char_my_rank
    logical, parameter :: ll_testwrite = .false.

! Begin executable statements -------

! 0. Preliminaries
  !  Get my process rank  
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  
  ! Find the total number of processes being used in this run
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr)

! 1.1 Get total number of gas particles
  call MPI_ALLREDUCE(sphData%npart, npart_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

! 1.2 Get number of gas particles for each process
  ALLOCATE( npart_arr(n_proc) )
  CALL MPI_GATHER(sphData%npart, 1, MPI_INTEGER, npart_arr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

! 1.3 Get total gas mass
  call MPI_ALLREDUCE(sphData%totalgasmass, totalgasmass_all, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr) 
  sphData%totalgasmass = totalgasmass_all

! 2.1 Get total number of point masses
  call MPI_ALLREDUCE(sphData%nptmass, nptmass_all, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

! 2.2 Get total number of point masses for each proccess
  ALLOCATE( nptmass_arr(n_proc) )
  CALL MPI_GATHER(sphData%nptmass, 1, MPI_INTEGER, nptmass_arr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

! 3. Communicate arrays xn, yn, zn, rhon, x, y, z, ptmass and keep in tmp storage

! 3.1 Gas particles
  ALLOCATE ( npart_displ(n_proc) )
  npart_displ(1)=0
  DO i=2, n_proc
     npart_displ(i) = npart_displ(i-1) + npart_arr(i-1)
  END DO

  ALLOCATE ( xn_tmp(npart_all)          )
  ALLOCATE ( yn_tmp(npart_all)          )
  ALLOCATE ( zn_tmp(npart_all)          )
  ALLOCATE ( rhon_tmp(npart_all)        )
  ALLOCATE ( temperature_tmp(npart_all) )
  ALLOCATE ( gasmass_tmp(npart_all)     )
  ALLOCATE ( hn_tmp(npart_all)          )

! Gather the particle data on process 0 then broadcast to all. For some reason MPI_ALLGATHERV 
! doesn't work  for doing this. 
  CALL MPI_GATHERV(sphData%xn, sphData%npart, MPI_DOUBLE_PRECISION, xn_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(xn_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%yn, sphData%npart, MPI_DOUBLE_PRECISION, yn_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(yn_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%zn, sphData%npart, MPI_DOUBLE_PRECISION, zn_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(zn_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%rhon, sphData%npart, MPI_DOUBLE_PRECISION, rhon_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(rhon_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%temperature, sphData%npart, MPI_DOUBLE_PRECISION, temperature_tmp(:), npart_arr(:), &
       npart_displ(:), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(temperature_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%gasmass, sphData%npart, MPI_DOUBLE_PRECISION, gasmass_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(gasmass_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  CALL MPI_GATHERV(sphData%hn, sphData%npart, MPI_DOUBLE_PRECISION, hn_tmp(:), npart_arr(:), npart_displ(:), &
                      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST(hn_tmp(:), npart_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

! 3.2 Point masses
  IF ( nptmass_all > 0 ) THEN

     ALLOCATE ( nptmass_displ(n_proc) )
     nptmass_displ(1)=0
     DO i=2, n_proc
        nptmass_displ(i) = nptmass_displ(i-1) + nptmass_arr(i-1)
     END DO

     ALLOCATE ( x_tmp(npart_all)      )
     ALLOCATE ( y_tmp(npart_all)      )
     ALLOCATE ( z_tmp(npart_all)      )
     ALLOCATE ( ptmass_tmp(npart_all) )

     CALL MPI_GATHERV(sphData%x, sphData%nptmass, MPI_DOUBLE_PRECISION, x_tmp(:), nptmass_arr(:), nptmass_displ(:), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST(x_tmp(:), nptmass_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     
     CALL MPI_GATHERV(sphData%y, sphData%nptmass, MPI_DOUBLE_PRECISION, y_tmp(:), nptmass_arr(:), nptmass_displ(:), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST(y_tmp(:), nptmass_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

     CALL MPI_GATHERV(sphData%z, sphData%nptmass, MPI_DOUBLE_PRECISION, z_tmp(:), nptmass_arr(:), nptmass_displ(:), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST(z_tmp(:), nptmass_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

     CALL MPI_GATHERV(sphData%ptmass, sphData%nptmass, MPI_DOUBLE_PRECISION, ptmass_tmp(:), nptmass_arr(:), nptmass_displ(:), &
          MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
     CALL MPI_BCAST(ptmass_tmp(:), nptmass_all, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  END IF

! Write out particle list for tests
  IF ( ll_testwrite ) THEN
     write(char_my_rank, '(i3)') my_rank
     open (unit=60, status='replace', file='mpi_test_tmp_'//TRIM(ADJUSTL(char_my_rank))//'.txt')
     do i=1, npart_all
        write(60,*) xn_tmp(i), yn_tmp(i), zn_tmp(i), rhon_tmp(i), temperature_tmp(i)
     end do
     IF ( nptmass_all > 0 ) THEN
        DO i=1, nptmass_all
           write(60, *) x_tmp(i), y_tmp(i), z_tmp(i), ptmass_tmp(i)
        END DO 
     END IF
     close(60)
  END IF

! 4. Deallocate and reAllocate arrays xn, yn, zn, rhon, x, y, z, ptmass

! 4.1 Gas particles
  DEALLOCATE ( sphData%xn          )
  DEALLOCATE ( sphData%yn          )
  DEALLOCATE ( sphData%zn          )
  DEALLOCATE ( sphData%rhon        )
  DEALLOCATE ( sphData%temperature )
  DEALLOCATE ( sphData%gasmass     )
  DEALLOCATE ( sphData%hn          ) 

  ALLOCATE   ( sphData%xn(npart_all)          )
  ALLOCATE   ( sphData%yn(npart_all)          )
  ALLOCATE   ( sphData%zn(npart_all)          )
  ALLOCATE   ( sphData%rhon(npart_all)        )
  ALLOCATE   ( sphData%temperature(npart_all) )
  ALLOCATE   ( sphData%gasmass(npart_all)     )
  ALLOCATE   ( sphData%hn(npart_all)          ) 

! 4.2 Point masses
  IF  ( nptmass_all > 0 ) THEN
     DEALLOCATE ( sphData%x   )
     DEALLOCATE ( sphData%y   )
     DEALLOCATE ( sphData%z   )
     DEALLOCATE ( sphData%ptmass )
     ALLOCATE   ( sphData%x(nptmass_all)      )
     ALLOCATE   ( sphData%y(nptmass_all)      )
     ALLOCATE   ( sphData%z(nptmass_all)      )
     ALLOCATE   ( sphData%ptmass(nptmass_all) )
  END IF

! 5. Populate new arrays with values from the tmp storage

! 5.1 Gas particles
  sphData%xn(:)          = xn_tmp(:)
  sphData%yn(:)          = yn_tmp(:)
  sphData%zn(:)          = zn_tmp(:)
  sphData%rhon(:)        = rhon_tmp(:)
  sphData%temperature(:) = temperature_tmp(:)
  sphData%npart          = npart_all
  npart                  = npart_all
  sphData%gasmass(:)     = gasmass_tmp(:)
  sphData%hn(:)          = hn_tmp(:)

! 5.2 Point masses
  IF  ( nptmass_all > 0 ) THEN
     sphData%x(:)      = x_tmp(:)
     sphData%y(:)      = y_tmp(:)
     sphData%z(:)      = z_tmp(:)
     sphData%ptmass(:) = ptmass_tmp(:)
     sphData%nptmass   = nptmass_all
  END IF

! Write out particle list for tests
  If ( ll_testwrite ) THEN
     open (unit=60, status='replace', file='mpi_test_'//TRIM(ADJUSTL(char_my_rank))//'.txt')
     do i=1, sphData%npart
        write(60,*) sphData%xn(i), sphData%yn(i), sphData%zn(i), sphData%rhon(i), sphData%temperature(i)
     end do
     IF ( nptmass_all > 0 ) THEN
        DO i=1, sphData%nptmass
           write(60, *) sphData%x(i), sphData%y(i), sphData%z(i), sphData%ptmass(i)
        END DO 
     END IF
     close(60)
     open (unit=60, status='replace', file='smoothing_length_'//TRIM(ADJUSTL(char_my_rank))//'.txt')
     do i=1, sphData%npart
        write(60,*) ( sphData%gasmass(i) / sphData%rhon(i) ) ** (1.0/3.0),  sphData%hn(i) 
     end do
     close(60)
  END IF

! 6. Deallocate temporary storage
  DEALLOCATE ( rhon_tmp        )
  DEALLOCATE ( temperature_tmp )
  DEALLOCATE ( zn_tmp          )
  DEALLOCATE ( yn_tmp          )
  DEALLOCATE ( xn_tmp          )
  DEALLOCATE ( hn_tmp          )
  DEALLOCATE ( gasmass_tmp     ) 
  IF  ( nptmass_all > 0 ) THEN
     DEALLOCATE ( x_tmp      )
     DEALLOCATE ( y_tmp      )
     DEALLOCATE ( z_tmp      )
     DEALLOCATE ( ptmass_tmp )
  END IF
  DEALLOCATE ( npart_displ )
  DEALLOCATE ( nptmass_arr )
  DEALLOCATE ( npart_arr   )

! 7. Communicate data to processes which did not run the sph step
  CALL MPI_BCAST(sphData%udist,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  CALL MPI_BCAST(sphData%utime,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  CALL MPI_BCAST(sphData%umass,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 
  CALL MPI_BCAST(sphData%time,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr) 

  end subroutine gather_sph_data

#else

! Dummy routines  for use in non-mpi case
  subroutine gather_sph_data

    USE sph_data_class
    implicit none

    return
  
  end subroutine gather_sph_data

#endif

end module torus_mod

! End of file ------------------------------------------------------------------

