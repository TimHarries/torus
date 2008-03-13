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

  subroutine torus(b_idim,b_npart,b_nactive,b_xyzmh,b_rho,b_iphase, &
                  b_nptmass,b_listpm,b_udist,b_umass, &
                  b_utime,b_time,b_gaspartmass, b_num_gas, b_temp)


  use constants_mod          ! physical constants
  use kind_mod               ! variable type KIND parameters
  use vector_mod             ! vector math
  use photon_mod             ! photon manipulation
  use gridtype_mod           ! type definition for the 3-d grid
  use grid_mod               ! opacity grid routines
  use phasematrix_mod        ! phase matrices and stokes vectors
  use math_mod               ! misc maths subroutines
  use blob_mod               ! clumps initialization and movement
  use distortion_mod         ! for distorting opacity grid
  use image_mod              ! stokes image
  use utils_mod
  use stateq_mod
  use inputs_mod
  use TTauri_mod
  use luc_cir3d_class
  use cmfgen_class
  use romanova_class
  use unix_mod
  use path_integral
  use puls_mod
  use input_variables         ! variables filled by inputs subroutine
  use lucy_mod
  use amr_mod
  use jets_mod
  use dust_mod
  use source_mod
  use spectrum_mod
  use wr104_mod
  use sph_data_class
  use cluster_class
  use timing
  use isochrone_class
  use filter_set_class
  use cluster_utils
  use surface_mod
  use disc_hydro_mod
  use disc_class
  use discwind_class
  use jet_class
  use gas_opacity_mod
  use diffusion_mod
  use messages_mod
  use photoion_mod
  use photoionAMR_mod
  use formal_solutions
  use starburst_mod
  use molecular_mod
  use modelatom_mod
  use cmf_mod
  use mpi_global_mod
  use parallel_mod

  implicit none
#ifdef MPI
   include 'mpif.h'  
#endif

  integer :: nOuterLoop = 10
  integer :: nSource
  type(SOURCETYPE), allocatable :: source(:)
  type(SOURCETYPE) a_star
  real ::  weightDust=1.0, weightPhoto=1.0
  real :: inclination
  real :: meanDustParticleMass

  character(len=80) :: thisdevice

  real(double) :: objectDistance

  ! variables for the grid

  type(GRIDTYPE) :: grid

  real(double) :: cOrecontinuumflux
  real :: probLinePhoton 

  integer :: isize
  integer, allocatable :: iseed(:)

  ! optical depth variables
  integer, parameter :: maxTau = 20000

  real :: sigmaExt0, sigmaAbs0, sigmaSca0  ! cross section at the line centre

  ! variables to do with dust

  integer :: itestlam, ismoothlam
  integer, parameter :: nXmie = 20, nMuMie = 1000
  type(PHASEMATRIX), allocatable :: miePhase(:,:, :)

  ! torus images

  type(PVIMAGETYPE), allocatable :: pvimage(:)

  ! intrinsic profile variables

  integer, parameter :: maxIntPro = 1000

  real, allocatable :: statArray(:)
  real, allocatable :: sourceSpectrum(:)
  real, allocatable :: sourceSpectrum2(:)
!  logical :: resonanceLine


  ! variables for clumped wind models
  
  integer, parameter :: maxBlobs = 10000
  real, parameter :: blobTime = 1000.
  real :: timeEnd = 24.*60.*60.
  real :: timeStart = 0.

  ! filenames

  character(len=80) :: filename
  character(len=80) :: phasePopFilename

  ! vectors

  type(VECTOR), parameter :: zAxis = VECTOR(0.,0.,1.)
  type(VECTOR), parameter :: yAxis = VECTOR(0.,1.,0.)
  type(VECTOR), parameter :: xAxis = VECTOR(1.,0.,0.)
  type(VECTOR) :: viewVec, outVec, originalViewVec
  type(VECTOR) :: rotationAxis, normToRotation
  type(VECTOR) :: zeroVec

  ! output arrays

  integer :: iLambda
  type(STOKESVECTOR), allocatable :: yArray(:)
  type(STOKESVECTOR), allocatable :: errorArray(:,:)
  real, allocatable :: xArray(:)

  ! model flags

  logical :: flatSpec
  logical :: ok
  logical :: greyContinuum
  logical :: firstPlot

  ! model parameters

  real :: foreground = 0., background = 0.  ! plotting intensities

  logical :: velocitySpace

  ! misc
  real :: meanr_line = 0., meanr_cont = 0.
  real :: wtot_line =0., wtot_cont = 0.
  real :: meanr0_line = 0., meanr0_cont = 0.
  real :: wtot0_line =0., wtot0_cont = 0.

  real(oct) :: t
  integer :: i, j
  real :: nu
  real(double) :: fac
  real :: rStar
  character(len=80) :: tempChar
  logical :: lineResAbs    ! T if you want to include absorption
  !                        !   of line by line resonance zones.

  !
  real(double) :: totalMass

  real :: theta1, theta2
  type(SURFACETYPE) :: starSurface

  ! Spot stuff

  real :: chanceDust = 0.

  ! adaptive grid stuff

  type(OCTALVECTOR) :: amrGridCentre ! central coordinates of grid
  type(OCTALVECTOR) :: octVec
  real :: ang
  integer           :: nOctals       ! number of octals in grid
  integer           :: nVoxels       ! number of unique voxels in grid
                                     !   (i.e. the number of childless subcells)
  logical :: gridConverged           ! true when adaptive grid structure has 
                                     !   been finalised
  integer           :: tooFewSamples ! number of errors from integratePathAMR
  integer           :: boundaryProbs ! number of errors from integratePathAMR
  integer           :: negativeOpacity ! number of errors from integratePathAMR
  character(len=80) :: newContFluxFile ! modified flux file (i.e. with accretion)
  logical :: alreadyDoneInfall = .false. ! whether we have already done an infall calculation
  type(alpha_disc)  :: ttauri_disc       ! parameters for ttauri disc

  ! For romanova geometry case
  type(romanova) :: romData ! parameters and data for romanova geometry

  !
  ! SPH data of Matthew
  type(sph_data) :: sphData

  ! Used for multiple sources (when geometry=cluster)
  type(cluster)   :: young_cluster

  ! Used in "plot_AMR_planes" and "plot_AMR_values"
  integer :: nmarker           ! number of markers
  real, allocatable    :: xmarker(:)  ! position of x
  real, allocatable    :: ymarker(:)  ! position of y
  real, allocatable    :: zmarker(:)  ! position of z
  real    :: width_3rd_dim         ! Use this to restrict the markers to be plotted..  
  real val_3rd_dim

  ! Name of the file to output various message from torus
  character(LEN=7), parameter :: messageFile = "MESSAGE"
  character(len=80) :: message

! Variables used when linking to sph code
  type(isochrone)       :: isochrone_data
  integer, intent(in)   :: b_idim,b_npart,b_nactive,b_nptmass
  integer*1, intent(in) :: b_iphase(b_idim)
  integer, intent(in)   :: b_listpm(b_nptmass)
  real*8, intent(in)    :: b_xyzmh(5,b_idim)
  real*4, intent(in)    :: b_rho(b_idim)
  real*8, intent(in)    :: b_udist, b_umass, b_utime, b_time, b_gaspartmass
  integer, intent(in)   :: b_num_gas           ! Number of gas particles
  real*8, intent(inout) :: b_temp(b_num_gas)   ! Temperature of gas particles
  integer :: ngaspart
  integer, save :: num_calls = 0

  type(modelatom), allocatable :: thisAtom(:)
  
  integer :: nRBBTrans
  integer :: indexRBBTrans(1000), indexAtom(1000)

  real(double) :: tempArray(10)
#ifdef MPI
  ! For MPI implementations =====================================================
  integer ::   ierr           ! error flag
  integer ::   tempInt        !
  integer, parameter :: tag = 0
  
! Begin executable statements --------------------------------------------------

  ! FOR MPI IMPLEMENTATION=======================================================

  ! Set up amrCOMMUNICATOR and global mpi groups
  call setupAMRCOMMUNICATOR

  call unixGetHostname(tempChar, tempInt) 
  print *, 'Process ', myRankGlobal,' running on host ',TRIM(ADJUSTL(tempChar))

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
! write(*,*) myRankGlobal,"random seed",iseed(1:isize)
#endif
  deallocate(iSeed)


  ! initialize

  useNdf = .true.
  sed = .false.
  jansky = .false.
  SIsed = .false.
  velocitySpace = .false.
  movie = .false.
  thinLine = .false.
  flatspec = .false.
  greyContinuum = .false.
  secondSource = .false.
  firstplot = .true.
  doRaman = .false.
  enhance = .false.
  lineResAbs = .false.  

  contFluxFile = "none"
  intProFilename = "none"

  inputKappaSca = 0.
  inputKappaAbs = 0.

  zeroVec = VECTOR(0.,0.,0.)

  nMarker = 0
  tooFewSamples = 0 
  boundaryProbs = 0 
  negativeOpacity = 0 
  lucyRadiativeEq = .false. ! this has to be initialized here
  num_calls = num_calls + 1

  ! get the model parameters

  call inputs() ! variables are passed using the input_variables module
  if (.not.inputOK) goto 666

  if (geometry /= "cluster" ) then
     print *, "Error: cluster geometry required"
     goto 666
  endif

  if (cmf) then
     allocate(thisAtom(1:nAtom))
     do i = 1, nAtom
        call readAtom(thisAtom(i),atomFilename(i))
        call stripAtomLevels(thisAtom(i), 5)
     enddo
    call createRBBarrays(nAtom, thisAtom, nRBBtrans, indexAtom, indexRBBTrans)
  endif

  objectDistance = griddistance * pctocm

  grid%photoionization = photoionization
  ! time elipsed since the beginning of the model
  grid%timeNow = phaseTime * REAL(nStartPhase-1)
  

  grid%resonanceLine = resonanceLine


  amrGridCentre = OCTALVECTOR(amrGridCentreX, amrGridCentreY, amrGridCentreZ)

  
  ! For plotting routines used later in this program
  if (plane_for_plot == "x-y") then
     val_3rd_dim = amrGridCentre%z
  else if (plane_for_plot == "y-z") then
     val_3rd_dim = amrGridCentre%x
  else if (plane_for_plot == "z-x") then
     val_3rd_dim = amrGridCentre%y
  else if (plane_for_plot == "x-z") then
     val_3rd_dim = amrGridCentre%y
  else
     val_3rd_dim = 0.0d0        
  end if


  if (doRaman) screened = .true.

  if (dopvimage) then
     allocate(pvimage(1:nSlit))
  endif

  probLinePhoton = 1. - probContPhoton

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

  scale = scale * rSol

  rStar = 100.*rSol

  ! switches for line emission

  if (lineEmission) then
     flatspec = .true.
     greyContinuum = .true.
     velocitySpace = .true.
  endif

  ! the observer's viewing direction

  rotationAxis = VECTOR(0., 0., 1.)

  normToRotation = rotationAxis .cross. yAxis
  call normalize(normToRotation)

  originalViewVec%x = 0.
  originalViewVec%y = -sin(inclination)
  originalViewVec%z = -cos(inclination)
  viewVec = originalViewVec

  !=====================================================================
  ! INIIALIZATIONS BEFORE CALLING INITAMR ROUTINE SHOULD BE HERE
  !=====================================================================

  !
  ! Special case
  !

     ! The total number of gas particles is the total number of active particles less the number of point masses.
  ngaspart = b_nactive-b_nptmass
  call init_sph_data2(sphData, b_udist, b_umass, b_utime, ngaspart, b_time, b_nptmass, &
       b_gaspartmass, b_npart, b_idim, b_iphase, b_xyzmh, b_rho, b_temp)
  ! Communicate particle data. Non-mpi case has a stubbed routine. 
  call gather_sph_data(sphData)

  ! Writing basic info of this data
  if (myRankIsZero) call info(sphData, "info_sph.dat")

  ! reading in the isochrone data needed to build an cluster object.
  call new(isochrone_data, "dam98_0225")   
  call read_isochrone_data(isochrone_data)
     
  ! making a cluster object
  call new(young_cluster, sphData, dble(amrGridSize), disc_on)
  call build_cluster(young_cluster, sphData, dble(lamstart), dble(lamend), isochrone_data)
    
  ! Wrting the stellar catalog readble for a human
  if (myRankIsZero) call write_catalog(young_cluster, sphData)

  ! Finding the inclinations of discs seen from +z directions...
  !     if (myRankIsZero) call find_inclinations(sphData, 0.0d0, 0.0d0, 1.0d0, "inclinations_z.dat")

  !==========================================================================
  !==========================================================================

  ! allocate the grid - this might crash out through memory problems

  if (gridUsesAMR) then
     ! any AMR allocation stuff goes here
     call initAMRGrid(greyContinuum, &
                        newContFluxFile,flatspec,grid,ok,theta1,theta2)
        if (.not.ok) goto 666

  else
     print *, "Error: AMR grid required"
     STOP
  end if

  grid%splitOverMpi = .false.


  grid%resonanceLine = resonanceLine
  grid%lambda2 = lamline

  ! Note: the first index should be either lambda or mu
  !       in order to speedup the array operations!!!  (RK) 
  allocate(miePhase(1:nDustType,1:nLambda,1:nMumie)) 
  allocate(statArray(1:nLambda))
  allocate(sourceSpectrum(1:nLambda))
  allocate(sourceSpectrum2(1:nLambda))

  grid%doRaman = doRaman

  statArray = 0.
  sourceSpectrum = 1.
  sourceSpectrum2 = 1.


  ! allocate the output arrays

  if (mie) then
     nOuterLoop = nLambda
  endif

  allocate(xArray(1:nLambda))
  allocate(yArray(1:nLambda))
  allocate(errorArray(1:nOuterLoop,1:nLambda))


  if (photoionization) then
     call addIons(grid%ion, grid%nion)
  endif

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


  errorArray(1:nOuterLoop,1:nLambda) = STOKESVECTOR(0.,0.,0.,0.)

  if (doRaman) then
     yarray(1:nLambda)%i = 1.e-20
     errorArray(1:nOuterLoop,1:nLambda)%i = 1.e-20
  endif


  call  createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)

  if (noScattering) then
     if (writeoutput) write(*,*) "! WARNING: Scattering opacity turned off in model"
     grid%oneKappaSca(1:nDustType,1:nLambda) = TINY(grid%oneKappaSca)
  endif

  if (writeoutput) then
     open(76, file="phasematrix.dat",status="unknown",form="formatted")
     ilambda = findIlambda(10000.0, xArray, nLambda, ok)
     do i = 1, nMuMie
	ang =  pi*real(i-1)/real(nMuMie-1)
        write(76,*) 180.*ang/pi, miephase(1, ilambda, i)%element(1,1)
     enddo
     close(76)
     tempArray(1) = 0.
     do i = 1, nMumie
        tempArray(1) = tempArray(1) + miePhase(1, ilambda, i)%element(1,1)/real(nMuMie)
     enddo
     write(*,*) "phase normalization ",tempArray(1)


  endif

  if (includeGasOpacity) then

     ! This routine may cause trouble with MPI routine! 
     ! Propably multiple nodes are trying to write to a same file at a same time.

     call readTsujiPPTable()
     call readTsujiKPTable()
     call createAllMolecularTables(20, 1.e-6, 20000., grid%nLambda, xArray)

  end if

  if (mie) then
     call locate(grid%lamArray, nLambda,lambdaTau,itestlam)
     call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)
     grid%itestlam = iTestLam
     write(message,*) "Test wavelength index: ",itestlam,ismoothlam
     call writeInfo(message, TRIVIAL)
  end if

  ! chris
  if (hydroWarp) then
     print *, "Creating and reading hydroGrid..."
     allocate(grid%hydroGrid)
     call initAMRGrid(greyContinuum, &
                        newContFluxFile,flatspec,grid%hydroGrid,ok,theta1,theta2)
     call readAMRgrid("hydrogrid.dat",readFileFormatted,grid%hydroGrid)
     print *, "hydroGrid read..."
     call hydroWarpFitSplines(grid%hydroGrid)
     print *, "splines fitted..."
     print *, "hydroGrid setup done."
  end if


  !===============================================================
  if (gridUsesAMR) then  !========================================
  !===============================================================


! below is really yucky code - please change me

     if (molecular .and. readmol .and. (.not. lucyRadiativeEq)) then
        if(.not. openlucy) then
           call readAMRgrid(molfilenamein,.false.,grid)
        else
           call readAMRgrid(lucyFileNamein,.false.,grid)
        endif
        goto 667
     endif

     if (cmf.and.readLucy) then
        call readAMRgrid("atom_tmp.grid",.false.,grid)
     else if (photoionization.and.readlucy) then
        continue
     else
        ! Set up the AMR grid
        call amr_grid_setup

        ! Write the information on the grid to file using a routine in grid_mod.f90
        if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")

        ! Plotting the various values stored in the AMR grid.
        call do_amr_plots
     end if

  !=================================================================
  else ! grid is not adaptive ======================================     
  !=================================================================


     print *, "ERROR: adaptive grid required when using torusMod"
     goto 666

  !=================================================================
  endif ! (gridusesAMR)
  !=================================================================


  outVec = (-1.)* originalViewVec

  ! set up the sources
  call set_up_sources

     if (associated(grid%octreeRoot)) then
        totalMass =0.d0
        call findTotalMass(grid%octreeRoot, totalMass)
        write(message,*) "Mass of envelope: ",totalMass/mSol, " solar masses"
        call writeInfo(message, IMPORTANT)
     endif

     call torus_mpi_barrier

667 continue
     call random_seed

     if (cmf) then
        if (movie) call plotAMRthreeDMovie(grid, source, nsource)
        if (.not.readlucy) call atomLoop(grid, nAtom, thisAtom, nsource, source)
!        call atomLoop(grid, nAtom, thisAtom, nsource, source)

           if (writeoutput) then
              call plot_AMR_values(grid, "n=3", plane_for_plot, val_3rd_dim,  &
                   "etaline.ps/vcps", .true., .false., & 
                   nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim)
              !
              
              fac = 120.*degtorad
              do i = 1, 1
                 ang = real(i-1)/40. * twoPi
                 octVec = OCTALVECTOR(sin(fac)*cos(ang), sin(fac)*sin(ang), cos(fac))
                 write(thisdevice,'(a,i3.3,a)') "col",i,".ps/vcps"
!                 call columnDensityPlot(grid, source, nsource, octVec, thisDevice)
              enddo
           endif

        do i = 1, 50
           ang = twoPi * real(i-1)/51. 
!           ang = 45.*degtorad
           t = 120.d0*degtorad
           call calculateAtomSpectrum(grid, thisAtom, nAtom, iTransAtom, iTransLine, &
                OCTALVECTOR(sin(t)*cos(ang), sin(t)*sin(ang), cos(t)), 100.d0*pctocm/1.d10, &
                source, nsource,i)
        enddo
        stop
     endif

  if (lucyRadiativeEq) then
     if (doTuning) call tune(6, "LUCY Radiative Equilbrium")  ! start a stopwatch
 
     if (.not.grid%adaptive) then
        call lucyRadiativeEquilibrium(grid, miePhase, nDustType, nMuMie, nLambda, xArray, dble(teff), nLucy)
     else
        if (readLucy .and. .not. redoLucy) then
           continue
        else

           if (solveVerticalHydro) then
              call verticalHydrostatic(grid, mCore, sigma0, rInner, miePhase, nDustType, nMuMie, nLambda, xArray, &
                   source, nSource, nLucy, massEnvelope, tThresh, .false., mDisc)
           else
              call lucyRadiativeEquilibriumAMR(grid, miePhase, nDustType, nMuMie, & 
                   nLambda, xArray, source, nSource, nLucy, massEnvelope, tthresh, &
                   lucy_undersampled, .false., IterLucy, plot_i=num_calls)
           endif

        endif

        if (myRankIsZero .and. writeLucy) call writeAMRgrid(lucyFilenameOut,writeFileFormatted,grid)
     endif

     if (doTuning) call tune(6, "LUCY Radiative Equilbrium")  ! stop a stopwatch

     call update_sph_temperature (b_idim, b_npart, b_iphase, b_xyzmh, sphData, grid, b_temp)

  endif

! Tidy up and finish the run 

666 continue

if (doTuning) call tune(6, "Torus Main") ! stop a stopwatch  

call writeInfo("TORUS exiting", FORINFO)

call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
call freeGrid(grid)

deallocate(miePhase) 
deallocate(statArray)
deallocate(sourceSpectrum)
deallocate(sourceSpectrum2)
deallocate(xArray)
deallocate(yArray)
deallocate(errorArray)

call torus_mpi_barrier

CONTAINS
!-----------------------------------------------------------------------------------------------------------------------

  subroutine set_up_lambda_array

    real :: deltaLambda
    real :: loglamStart, logLamEnd

    if (lamLinear) then
       deltaLambda = (lamEnd - lamStart) / real(nLambda)
     
       xArray(1) = lamStart + deltaLambda/2.
       yArray(1)%i = 0.
       yArray(1)%q = 0.
       yArray(1)%u = 0.
       yArray(1)%v = 0.
       do i = 2, nLambda
          xArray(i) = xArray(i-1) + deltaLambda
          yArray(i)%i = 0.
          yArray(i)%q = 0.
          yArray(i)%u = 0.
          yArray(i)%v = 0.
       enddo

    else

       logLamStart = log10(lamStart)
       logLamEnd   = log10(lamEnd)
       
       do i = 1, nLambda
          xArray(i) = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
          xArray(i) = 10.**xArray(i)
          yArray(i) = STOKESVECTOR(0.,0.,0.,0.)
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

    if (doTuning) call tune(6, "AMR grid construction.")  ! start a stopwatch

    if (readPops .or. readPhasePops .or. readLucy) then 

       print *, "Error: TorusMod not configures with readPops, readPhasePops, readLucy" 
       STOP

    else  ! not reading a population file

       amrGridCentre = octalVector(amrGridCentreX,amrGridCentreY,amrGridCentreZ)
       call writeInfo("Starting initial set up of adaptive grid...", TRIVIAL)
       
       call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, sphData, young_cluster, nDustType)
       call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, sphData, young_cluster)
       call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
       call estimateRhoOfEmpty(grid, grid%octreeRoot, sphData)	
       !Removing the cells within 10^14 cm from the stars.
       call remove_too_close_cells(young_cluster,grid%octreeRoot,1.0d4)

   
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

        if (myRankIsZero) call grid_info(grid, "info_grid.dat")

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

! Plot values of variables on AMR grid. Extracted from torusMain by D. Acreman 

subroutine do_amr_plots

  ! Do some preparation for the arrays used in plot_AMR_* which will be used later
  if (grid%geometry(1:7) == "cluster") then
     nmarker = get_nstar(young_cluster)
     allocate(xMarker(1:nMarker))
     allocate(yMarker(1:nMarker))
     allocate(zMarker(1:nMarker))
     do i = 1, nmarker
        a_star = get_a_star(young_cluster, i)
        xmarker(i)= a_star%position%x
        ymarker(i)= a_star%position%y
        zmarker(i)= a_star%position%z
     end do
  else
     nmarker = 0
     ALLOCATE(xmarker(nmarker), ymarker(nmarker), zmarker(nmarker))
  end if
  width_3rd_dim = amrGridSize

  if ( plot_maps .and. myRankIsZero ) then

     call plot_AMR_values(grid, "rho", plane_for_plot, val_3rd_dim, &
          "rho_grid",.true., .true., nmarker, xmarker, ymarker, zmarker, &
          width_3rd_dim, show_value_3rd_dim, suffix="default", index=num_calls, &
          fixValMin=sph_rho_min, fixValMax=sph_rho_max, useFixedRange=.true. )
     call plot_AMR_values(grid, "rho", plane_for_plot, val_3rd_dim, &
          "rho_zoom",.true., .true., nmarker, xmarker, ymarker, zmarker, &
          width_3rd_dim, show_value_3rd_dim, boxfac=zoomFactor, suffix="default", index=num_calls, &
          fixValMin=sph_rho_min, fixValMax=sph_rho_max, useFixedRange=.true. )
     if ((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry=="warpeddisc")) then
        call plot_AMR_values(grid, "rho", plane_for_plot, val_3rd_dim, &
             "rho_ultrazoom",.true., .true., nmarker, xmarker, ymarker, zmarker, &
             width_3rd_dim, show_value_3rd_dim, boxfac=0.005, suffix="default", index=num_calls)
     end if
     call plot_AMR_planes(grid, "rho", plane_for_plot, 3, "rho", .true., .false., &
          nmarker, xmarker, ymarker, zmarker, show_value_3rd_dim)
     if (grid%lineEmission) then
        call plot_AMR_values(grid, "Vx", plane_for_plot, val_3rd_dim,  &
             "Vx", .false., .false.,  &
             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim, suffix="default")
        call plot_AMR_values(grid, "Vy", plane_for_plot, val_3rd_dim, &
             "Vy", .false., .false., &
             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim, suffix="default")
        call plot_AMR_values(grid, "Vz", plane_for_plot, val_3rd_dim, &
             "Vz", .false., .false., &
             nmarker, xmarker, ymarker, zmarker, width_3rd_dim, show_value_3rd_dim, suffix="default")
     end if

     call plot_AMR_values(grid, "temperature", plane_for_plot, val_3rd_dim, &
          "temperature", .true., .false., nmarker, xmarker, ymarker, zmarker, &
          width_3rd_dim, show_value_3rd_dim, suffix="default", index=num_calls, &
          fixValMin=sph_tem_min, fixValMax=sph_tem_max, useFixedRange=.false. )
  
     if (mie) then
        do i = 1, nDustType
           write(message,'(a,i1)') "dusttype",i
           write(filename,'(a,i1)') "dusttype",i
           call plot_AMR_values(grid, message, "x-z", 0., &
                trim(filename),.false., .false., nmarker, xmarker, ymarker, zmarker, &
                width_3rd_dim, show_value_3rd_dim, boxfac=zoomfactor, suffix="default", index=num_calls)
        enddo
     endif

  endif

end subroutine do_amr_plots

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
     call randomSource(source, nSource, i, xArray, nLambda, initialize=.true.)  
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

