!
!            T   O   R   U  S
!
! torus is a general purpose monte-carlo radiative transfer code used to
! produce polarized spectra of hot and cool stellar winds
!

! written by tjh
! Extracted from torusMain by D. Acreman (March 2008)
 
module torus_mod

  implicit none 

  public :: torus

  private

contains

  subroutine torus(b_idim,  b_npart,       b_nactive, b_nptmass, b_num_gas, &
                   b_xyzmh, b_rho,         b_iphase,                        &
                   b_udist, b_umass,       b_utime,   b_time,    b_temp,    &
                   temp_min, b_totalgasmass, file_tag )

  use kind_mod               ! variable type KIND parameters
  use vector_mod             ! vector math
  use photon_mod             ! photon manipulation
  use gridtype_mod           ! type definition for the 3-d grid
  use grid_mod               ! opacity grid routines
  use phasematrix_mod        ! phase matrices and stokes vectors
  use math_mod               ! misc maths subroutines
  use utils_mod
  use stateq_mod
  use inputs_mod
  use cmfgen_class
  use unix_mod
  use path_integral
  use puls_mod
  use input_variables         ! variables filled by inputs subroutine
  use lucy_mod
  use amr_mod
  use dust_mod
  use source_mod, only: SOURCETYPE
  use spectrum_mod
  use sph_data_class
  use cluster_class
  use timing
  use isochrone_class
  use cluster_utils
  use surface_mod, only: SURFACETYPE
  use disc_class
  use gas_opacity_mod
  use diffusion_mod
  use messages_mod
  use formal_solutions
  use modelatom_mod
  use cmf_mod
  use mpi_global_mod
  use parallel_mod
  use vtk_mod
  use phaseloop_mod, only: do_phaseloop
  use TTauri_mod, only: alpha_disc

  implicit none
#ifdef MPI
   include 'mpif.h'  
#endif

  integer :: nSource
  type(SOURCETYPE), allocatable :: source(:)
  type(SOURCETYPE) a_star

  ! variables for the grid

  type(GRIDTYPE) :: grid

  integer :: isize
  integer, allocatable :: iseed(:)

  ! variables to do with dust
  integer :: itestlam, ismoothlam
  integer, parameter :: nMuMie = 1000
  type(PHASEMATRIX), allocatable :: miePhase(:,:,:)

  ! filenames
  character(len=80) :: phasePopFilename

  ! model flags
  logical :: flatSpec
  logical :: ok
  logical :: greyContinuum

  ! model parameters
  integer :: i, j
  character(len=80) :: tempChar
  
  real :: theta1, theta2
  type(SURFACETYPE) :: starSurface

  ! adaptive grid stuff
  type(OCTALVECTOR) :: amrGridCentre ! central coordinates of grid
  integer           :: nOctals       ! number of octals in grid
  integer           :: nVoxels       ! number of unique voxels in grid
                                     !   (i.e. the number of childless subcells)
  logical :: gridConverged           ! true when adaptive grid structure has 
                                     !   been finalised
  character(len=80) :: newContFluxFile ! modified flux file (i.e. with accretion)

  ! For romanova geometry case
  type(romanova) :: romData ! parameters and data for romanova geometry

  type(alpha_disc)  :: ttauri_disc       ! parameters for ttauri disc
  real(double) :: objectDistance
  real :: inclination
  integer, parameter :: maxTau = 20000
  type(BLOBTYPE), allocatable :: blobs(:)
  real dTime
  !
  ! SPH data of Matthew
  type(sph_data) :: sphData

  ! Used for multiple sources (when geometry=cluster)
  type(cluster)   :: young_cluster

  ! Name of the file to output various message from torus
  character(len=80) :: message

! Variables used when linking to sph code
  type(isochrone)       :: isochrone_data
  integer, intent(in)   :: b_idim,b_npart,b_nactive,b_nptmass
  integer*1, intent(in) :: b_iphase(b_idim)
  real*8, intent(in)    :: b_xyzmh(5,b_idim)
  real*4, intent(in)    :: b_rho(b_idim)
  real*8, intent(in)    :: b_udist, b_umass, b_utime, b_time
  integer, intent(in)   :: b_num_gas           ! Number of gas particles
  real*8, intent(inout) :: b_temp(b_num_gas)   ! Temperature of gas particles
  real*8, intent(in)    :: b_totalgasmass      ! Total gas mass for this MPI process
  real(kind=8), intent(in)      :: temp_min
  character(len=11), intent(in) :: file_tag
  character(len=11), save       :: prev_file_tag="none"
  integer, save :: num_calls = 0
  character(len=4) :: char_num_calls
  character(len=100) :: filename

#ifdef MPI
  ! For MPI implementations =====================================================
  integer ::   ierr           ! error flag
  integer ::   tempInt        !
  
! Begin executable statements --------------------------------------------------

  ! FOR MPI IMPLEMENTATION=======================================================

  ! Set up amrCOMMUNICATOR and global mpi groups
  call setupAMRCOMMUNICATOR

  call unixGetHostname(tempChar, tempInt) 
  print *, 'Process ', myRankGlobal,' running on host ',TRIM(ADJUSTL(tempChar))
  print *, 'Process ', myRankGlobal, 'b_npart=', b_npart, 'b_nactive=', b_nactive, &
           'b_nptmass=', b_nptmass, 'b_num_gas=', b_num_gas

  !===============================================================================

#endif

  writeoutput    = .true.
  outputwarnings = .true.
  outputinfo     = .true.
  doTuning       = .true.
  myRankIsZero   = .true.
#ifdef MPI
  if (myRankGlobal/=1) writeoutput  = .false.
  if (myRankGlobal/=1) doTuning     = .false.
  if (myRankGlobal/=0) myRankIsZero = .false.
#endif
  
  ! For time statistics
  if (doTuning) call tune(6, "Torus Main") ! start a stopwatch  

  ! set up a random seed
  
  call random_seed

  call random_seed(size=iSize)
  allocate(iSeed(1:iSize))
  call random_seed(get=iSeed)

  call torus_mpi_barrier
#ifdef MPI
  call MPI_BCAST(iSeed, iSize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call random_seed(put=iseed)
#endif
  deallocate(iSeed)


  ! initialize

  useNdf = .true.
  sed = .false.
  jansky = .false.
  SIsed = .false.
  movie = .false.
  thinLine = .false.
  flatspec = .false.
  greyContinuum = .false.
  secondSource = .false.
  doRaman = .false.
  enhance = .false.

  contFluxFile = "none"
  intProFilename = "none"

  inputKappaSca = 0.
  inputKappaAbs = 0.

  lucyRadiativeEq = .false. ! this has to be initialized here

! Used to  set up output file names
! Reset the counter each time the SPH dump file name changes 
  if ( file_tag /= prev_file_tag ) num_calls = 0
  prev_file_tag = file_tag
  num_calls = num_calls + 1
  write(char_num_calls,'(i4.4)') num_calls

  ! get the model parameters

  call inputs() ! variables are passed using the input_variables module
  if (.not.inputOK) goto 666

! check that the minimum temperature used in the SPH code matches that used in Torus
  if ( temp_min /= TminGlobal ) then
     write(*,*) "Minimum temperature in the SPH code does not match that in Torus"
     write(*,*) "Torus: TminGlobal=", TminGlobal
     write(*,*) "SPH: temp_min=", temp_min
     STOP
  endif

  if (geometry /= "cluster" ) then
     print *, "Error: cluster geometry required"
     goto 666
  endif

  if (cmf) then
     print *, "Error: not compatible with cmf"
     goto 666
  endif

  if (.not. gridUsesAMR) then
     print *, "Error: AMR grid required"
     STOP
  end if

  grid%photoionization = photoionization
  ! time elipsed since the beginning of the model
  grid%timeNow = phaseTime * REAL(nStartPhase-1)
  

  grid%resonanceLine = resonanceLine


  amrGridCentre = OCTALVECTOR(amrGridCentreX, amrGridCentreY, amrGridCentreZ)


  if (doRaman) screened = .true.

  scale = scale * rSol

  ! switches for line emission

  if (lineEmission) then
     flatspec = .true.
     greyContinuum = .true.
  endif

  if (allocated(inclinations)) then
     inclination = inclinations(1)
  else
     inclination = firstInclination
  end if
  if (inclination .lt. 1.e-6) then
     call writeWarning("inclination less than 1.e-6, now set to 1.e-6.")
     call writeWarning("Did you specify an inclination of 0?")
  end if
  inclination = max(inclination, 1.e-4)

  !=====================================================================
  ! INIIALIZATIONS BEFORE CALLING INITAMR ROUTINE SHOULD BE HERE
  !=====================================================================

  ! The total number of gas particles is the total number of active particles less the number of point masses.
  call init_sph_data2(sphData, b_udist, b_umass, b_utime, b_num_gas, b_time, b_nptmass, &
       b_npart, b_idim, b_iphase, b_xyzmh, b_rho, b_temp, b_totalgasmass)
  ! Communicate particle data. Non-mpi case has a stubbed routine. 
  call gather_sph_data(sphData)

  ! Writing basic info of this data
  if (myRankIsZero) call info(sphData, "info_sph_"//trim(adjustl(file_tag))//"."//TRIM(ADJUSTL(char_num_calls))//".dat")

  ! reading in the isochrone data needed to build an cluster object.
  call new(isochrone_data, "dam98_0225")   
  call read_isochrone_data(isochrone_data)
     
  ! making a cluster object
  call new(young_cluster, sphData, dble(amrGridSize), disc_on)
  call build_cluster(young_cluster, sphData, dble(lamstart), dble(lamend), isochrone_data)
    
  ! Wrting the stellar catalog readble for a human
  filename="catalogue_"//trim(adjustl(file_tag))//"."//TRIM(ADJUSTL(char_num_calls))//".dat"
  if (myRankIsZero) call write_catalog(young_cluster, sphData, trim(filename) )

  ! Finding the inclinations of discs seen from +z directions...
  !     if (myRankIsZero) call find_inclinations(sphData, 0.0d0, 0.0d0, 1.0d0, "inclinations_z.dat")

  !==========================================================================
  !==========================================================================

  ! allocate the grid - this might crash out through memory problems

     ! any AMR allocation stuff goes here
  call initAMRGrid(greyContinuum, &
       newContFluxFile,flatspec,grid,ok,theta1,theta2)
  if (.not.ok) goto 666

  grid%splitOverMpi = .false.


  grid%resonanceLine = resonanceLine
  grid%lambda2 = lamline

  ! Note: the first index should be either lambda or mu
  !       in order to speedup the array operations!!!  (RK) 
  allocate(miePhase(1:nDustType,1:nLambda,1:nMumie)) 

  grid%doRaman = doRaman

  if (photoionization) call addIons(grid%ion, grid%nion)
  

! Set up the wavelength arrays with either linear or log spaced values.
  call set_up_lambda_array

  !
  ! Setting the number of opacity (kappa) arrays in grid.
  !
  if (flatspec) then
     grid%nopacity = 1
  else
     grid%nopacity = nLambda
  end if

  call  createDustCrossSectionPhaseMatrix(grid, grid%lamArray, nLambda, miePhase, nMuMie)

  if (noScattering) then
     if (writeoutput) write(*,*) "! WARNING: Scattering opacity turned off in model"
     grid%oneKappaSca(1:nDustType,1:nLambda) = TINY(grid%oneKappaSca)
  endif

  if (includeGasOpacity) then

     ! This routine may cause trouble with MPI routine! 
     ! Propably multiple nodes are trying to write to a same file at a same time.

     call readTsujiPPTable()
     call readTsujiKPTable()
     call createAllMolecularTables(20, 1.e-6, 20000., grid%nLambda, grid%lamArray)

  end if

  if (mie) then
     call locate(grid%lamArray, nLambda,lambdaTau,itestlam)
     call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)
     grid%itestlam = iTestLam
     write(message,*) "Test wavelength index: ",itestlam,ismoothlam
     call writeInfo(message, TRIVIAL)
  end if

  ! Set up the AMR grid
  call amr_grid_setup

  ! Write the information on the grid to file using a routine in grid_mod.f90
  filename = trim ("info_grid_"//trim(adjustl(file_tag))//"."//TRIM(ADJUSTL(char_num_calls))//".dat")
  if ( myRankIsZero ) call grid_info(grid, filename)

  ! set up the sources
  call set_up_sources

  call checkSphTotalMass(grid, sphData, 2.0, ok)
  if ( .not. ok ) goto 666

  call torus_mpi_barrier

  call random_seed

  filename = trim ( "torus_in_"//trim(adjustl(file_tag))//'.'//char_num_calls//".vtk" )
  if (myRankIsZero) call  writeVtkFile(grid, filename, "vtk.txt")

  if (doTuning) call tune(6, "LUCY Radiative Equilbrium")  ! start a stopwatch
  
  call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, & 
       nLambda, grid%lamArray, source, nSource, nLucy, massEnvelope, tthresh, &
       lucy_undersampled, .false., IterLucy )

  if (myRankIsZero .and. writeLucy) call writeAMRgrid(lucyFilenameOut,writeFileFormatted,grid)

  if (doTuning) call tune(6, "LUCY Radiative Equilbrium")  ! stop a stopwatch

  filename = trim ( "torus_out_"//trim(adjustl(file_tag))//'.'//char_num_calls//".vtk" )
  if (myRankIsZero) call  writeVtkFile(grid, filename, "vtk.txt")

  call update_sph_temperature (b_idim, b_npart, b_iphase, b_xyzmh, sphData, grid, b_temp)

  if ( nInclination > 0 ) then

     call writeInfo ("Calling phaseloop", FORINFO)
     objectDistance = griddistance * pctocm
     call do_phaseloop(grid, .false., 0.0, 0.0, 0.0,                                            &
          theta1, theta2, VECTOR(0.0,0.0,0.0), 0.0_db, 0.0, 0.0, 0.0, 0.0_db,                     &
          starsurface, newContFluxFile, 0.0, 0.0, 0.0, ttauri_disc, (/VECTOR(0.0,0.0,0.0)/), 1,   &
          0.0, 10000, flatspec, objectDistance, inclination, maxTau,                              &
          miePhase, nsource, source, blobs, nmumie, dTime)

  else
     call writeInfo ("Not calling phaseloop", FORINFO)
  end if

! Tidy up and finish the run 

666 continue

if (doTuning) call tune(6, "Torus Main") ! stop a stopwatch  

call writeInfo("TORUS exiting", FORINFO)

call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
call freeGrid(grid)

deallocate(miePhase) 


call torus_mpi_barrier

CONTAINS
!-----------------------------------------------------------------------------------------------------------------------

  subroutine set_up_lambda_array

    real :: deltaLambda
    real :: loglamStart, logLamEnd
    real :: xarray(nLambda)

    if (lamLinear) then
       deltaLambda = (lamEnd - lamStart) / real(nLambda)
     
       xArray(1) = lamStart + deltaLambda/2.
       do i = 2, nLambda
          xArray(i) = xArray(i-1) + deltaLambda
       enddo

    else

       logLamStart = log10(lamStart)
       logLamEnd   = log10(lamEnd)
       
       do i = 1, nLambda
          xArray(i) = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
          xArray(i) = 10.**xArray(i)
       enddo

    endif


       if (lamFile) then
          call writeInfo("Reading wavelength points from file.", TRIVIAL)
          open(77, file=lamfilename, status="old", form="formatted")
          nLambda = 0
333       continue
          nLambda = nLambda + 1
          read(77,*,end=334) xArray(nLambda)
          if (writeoutput) write(*,*) nlambda,xArray(nlambda)
          goto 333
334       continue
          nLambda = nLambda - 1
          close(77)
       endif

    !
    ! Copying the wavelength array to the grid
    do i = 1, nLambda
       grid%lamArray(i) = xArray(i)
    enddo
    grid%nLambda = nLambda

  end subroutine set_up_lambda_array

!-----------------------------------------------------------------------------------------------------------------------

  subroutine amr_grid_setup

    use constants_mod, only: mSol

    real(double) :: removedMass

    if (doTuning) call tune(6, "AMR grid construction.")  ! start a stopwatch

    if (readPops .or. readPhasePops .or. readLucy) then 

       print *, "Error: TorusMod not configured with readPops, readPhasePops, readLucy" 
       STOP

    else  ! not reading a population file

       amrGridCentre = octalVector(amrGridCentreX,amrGridCentreY,amrGridCentreZ)
       call writeInfo("Starting initial set up of adaptive grid...", TRIVIAL)
       
       call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, sphData, young_cluster, nDustType)
       call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, sphData, young_cluster)
       call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
       call estimateRhoOfEmpty(grid, grid%octreeRoot, sphData)	
       !Removing the cells within 10^14 cm from the stars.
       removedMass = 0.0
       call remove_too_close_cells(young_cluster,grid%octreeRoot,real(grid%rCore, kind=db), removedMass, amr_min_rho)
       write(message,*) "Mass removed by remove_too_close_cells= ", removedMass / mSol
       call writeInfo(message, TRIVIAL)
   
        nOctals = 0
        nVoxels = 0
        call countVoxels(grid%octreeRoot,nOctals,nVoxels)
        write(message,*) "Adaptive grid contains: ",nOctals," octals"
        call writeInfo(message, TRIVIAL)
        write(message,*) "                      : ",nVoxels," unique voxels"
        call writeInfo(message, TRIVIAL)
        grid%nOctals = nOctals

        call writeInfo("Calling routines to finalize the grid variables...",TRIVIAL)
        gridConverged = .false.
     
        do
           call finishGrid(grid%octreeRoot,grid,gridConverged,romData=romData)
           if (gridConverged) exit
        end do        

        call writeInfo("...final adaptive grid configuration complete",TRIVIAL)

        if (writeoutput) call writeAMRgrid("after_creation.grid",.false.,grid)


        if (lineEmission.and.(.not.cmf)) then
           !  calculate the statistical equilibrium (and hence the emissivities 
           !  and the opacities) for all of the subcells in an
           !  adaptive octal grid.
           !  Using a routine in stateq_mod module.
           if (writeoutput) write(*,*) "Calling statistical equilibrium routines..."
           if (doTuning) call tune(6, "amrStateq") ! start a stopwatch  

           call amrStateq(grid, newContFluxFile, lte, nLower, nUpper, &
                      starSurface, recalcPrevious=.false., ion_name=ion_name, ion_frac=ion_frac)

           call torus_mpi_barrier('waiting for other amrStatEq calls to return...')

           if (doTuning) call tune(6, "amrStateq") ! stop a stopwatch  
           if (writeoutput) write(*,*) "... statistical equilibrium routines complete"

        end if ! (lineEmission)

        !
        ! cleaning up unused memoryusing the routine in sph_data_class
        call kill(sphData)

        if (myRankIsZero) call delete_particle_lists(grid%octreeRoot)

     end if ! (readPops .or. readPhasePops)


  if (myRankIsZero)  then
     ! writing pop files after statEq routines
     if (writePhasePops) then
        write(tempChar,'(i3.3)') nStartPhase
        phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
        call writeAMRgrid(phasePopFilename,writeFileFormatted,grid)
     end if
     if (writePops) then
        call writeAMRgrid(popFilename,writeFileFormatted,grid)
     end if

  end if
  call torus_mpi_barrier ! sync here
  
  if (doTuning) call tune(6, "AMR grid construction.") ! stop a stopwatch

end subroutine amr_grid_setup

!-----------------------------------------------------------------------------------------------------------------------
subroutine set_up_sources

  integer      :: nstar

  ! set up the sources
  nSource = 0
 
  ! Extract some info from cluster object.
  nstar = get_nstar(young_cluster)  ! number of stars in the cluster
  nSource =  n_stars_in_octal(young_cluster, grid%octreeRoot)
       
  ! copy the star over in the array.
  ! This is ugly. Maybe lucyRadiativeEquilibriumAMR should be changed to take
  ! an cluster_class object as an input variable in future.
  ALLOCATE(source(nSource))
  
  ! Restricting the source to be within the root cell (in case the root cell is 
  ! is smaller than the sph model space!
  j = 0 
  do i = 1, nstar
     a_star = get_a_star(young_cluster, i)
     ! using a function in source_mod
     if ( source_within_octal(a_star, grid%octreeRoot) ) then
        j = j+1
        source(j) = a_star
     end if
  end do

  if (nSource > 0) then
     call randomSource(source, nSource, i, grid%lamArray, nLambda, initialize=.true.)  
     call writeInfo("Sources set up.",TRIVIAL)
  endif


end subroutine set_up_sources

!-----------------------------------------------------------------------------------------------------------------------

end subroutine torus

end module torus_mod

!-----------------------------------------------------------------------------------------------------------------------

! Name:    update_sph_temperature
! Purpose: Update the temperatures of the SPH particles using the torus temperature field
!          Each process works on its own subset of the total number of particles.
! Author:  D. Acreman, November 2007

  subroutine update_sph_temperature (b_idim, b_npart, b_iphase, b_xyzmh, sphData, grid, b_temp)

    USE vector_mod, only:     octalVector
    USE amr_mod, only:        amrGridValues
    USE gridtype_mod, only:   gridType
    USE sph_data_class, only: sph_data, get_udist
    USE messages_mod

    implicit none

#ifdef MPI
! MPI specific variables
    include 'mpif.h'
    real    :: mpi_max_deltaT, mpi_sum_deltaT
    integer :: ierr, mpi_iiigas
#endif

! Arguments 
    integer, intent(in)   :: b_idim, b_npart
    integer*1, intent(in) :: b_iphase(b_idim)
    real*8, intent(in)    :: b_xyzmh(5,b_idim)
    real*8, intent(inout) :: b_temp(b_idim)
    type(sph_data), intent(in) :: sphData
    type(GRIDTYPE), intent(in) :: grid

! Local variables
    character(len=80) :: message
    integer      :: iiigas, i
    real*8       :: xgas, ygas, zgas
    real(double) :: sphDistFac
    real         :: tgas
    real         :: deltaT, sum_deltaT, mean_deltaT, max_deltaT
    type(OCTALVECTOR) :: octVec

! Begin executable statements

! 1. Calculate conversion from sph distance units to torus distance units
       sphDistfac  = get_udist(sphData) ! [cm]
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
             octVec = OCTALVECTOR(xgas,ygas,zgas)
             call amrGridValues(grid%octreeRoot, octVec, temperature=tgas, grid=grid)
! 2.2 Calculate statistics of temperature change
             deltaT     = tgas - b_temp(iiigas)
             sum_deltaT = sum_deltaT + deltaT
             max_deltaT = MAX(max_deltaT, deltaT)
! 2.3 Update the gas particle temperature to pass back to sph code
             b_temp(iiigas) = tgas
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
! Check the total mass on the grid against the total mass component of the SPH data structure.
! 
! D. Acreman, July 2008
   
  subroutine checkSphTotalMass(grid, sphData, threshold, ok)

    use gridtype_mod, only: gridtype
    use sph_data_class, only : sph_data
    use constants_mod, only: mSol
    use messages_mod
    use amr_mod

    implicit none

    type(GRIDTYPE), intent(in) :: grid
    type(sph_data), intent(in) :: sphData
    real, intent(in)           :: threshold
    logical, intent(out)       :: ok

    real(double) :: totalMass
    real(double) :: diff

    character(len=100) :: message

    totalMass = 0.d0
    call findTotalMass(grid%octreeRoot, totalMass)
    totalMass = totalMass / mSol

    write(message,*) "Mass of envelope: ",totalMass, " solar masses"
    call writeInfo(message, FORINFO)
    write(message,*) "Mass from SPH data structure: ", sphData%totalgasmass
    call writeInfo(message, FORINFO)

    ok = .true. 
    diff = ( abs(totalMass - sphData%totalgasmass) ) / sphData%totalgasmass
    if ( diff > threshold ) then
       ok = .false. 
       write(message,*) "Difference between mass on grid and mass from particles exceeds allowable limit."
       call writeFatal(message)
    end if

  end subroutine checkSphTotalMass

! End of file ------------------------------------------------------------------

