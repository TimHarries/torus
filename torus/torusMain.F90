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
  use utils_mod
  use input_variables         ! variables filled by inputs subroutine
  use constants_mod

  use amr_mod, only: amrGridValues, deleteOctreeBranch, hydroWarpFitSplines, polardump, pathtest, &
       setupNeighbourPointers, findtotalmass

  use gridtype_mod, only: gridType         ! type definition for the 3-d grid
  use grid_mod, only: initCartesianGrid, initPolarGrid, plezModel, initAMRgrid, freegrid, grid_info
  use phasematrix_mod, only: phasematrix        ! phase matrices
  use blob_mod, only: blobType
  use inputs_mod, only: inputs
  use TTauri_mod, only: infallenhancment, initInfallEnhancement 
  use romanova_class, only: romanova
  use dust_mod, only: createDustCrossSectionPhaseMatrix, MieCrossSection 
  use source_mod, only: sourceType, buildSphere 
  use sph_data_class, only: read_sph_data, kill, sphdata, clusterparameter
  use cluster_class, only: cluster
  use surface_mod, only: surfaceType
  use disc_class, only: alpha_disc, new, add_alpha_disc, finish_grid, turn_off_disc
  use discwind_class, only: discwind, new, add_discwind
  use jet_class, only: jet, new, add_jet, finish_grid_jet, turn_off_jet
  use photoion_mod, only: photoIonizationloop
  use molecular_mod, only: moleculetype, calculatemoleculespectrum, molecularloop, readmolecule, make_h21cm_image
  use modelatom_mod, only: modelAtom, createrbbarrays, readatom, stripatomlevels
  use cmf_mod, only: atomLoop, calculateAtomSpectrum
  use vtk_mod, only: writeVtkfile
  use mpi_global_mod, only: myRankGlobal
  use parallel_mod, only: torus_mpi_barrier
  use gridio_mod, only: writeAMRgrid, readamrgrid
  use bitstring_mod, only: constructBitStrings
  use gas_opacity_mod, only: createAllMolecularTables, readTsujiKPTable, readTsujiPPTable
  use density_mod, only: calcPlanetMass
  use phaseloop_mod, only: do_phaseloop
  use timing, only: tune
  use ion_mod, only: addions
  use angularImage, only: make_angular_image
#ifdef MPI
  use photoionAMR_mod, only: radiationhydro, createImageSplitGrid
  use hydrodynamics_mod, only: doHydrodynamics1d, doHydrodynamics2d, doHydrodynamics3d, readAMRgridMpiALL 
  use mpi_amr_mod, only: setupAMRCOMMUNICATOR, findMassOverAllThreads, grid_info_mpi
  use unix_mod, only: unixGetHostname
  use parallel_mod, only: sync_random_seed
#endif

  implicit none
#ifdef MPI
   include 'mpif.h'  
#endif

  integer :: nSource
  type(SOURCETYPE), allocatable :: source(:)
  real :: inclination

  ! variables for the grid

  type(GRIDTYPE) :: grid

  real(double) :: cOrecontinuumflux

  ! optical depth variables
  !
  ! for monte calro with high optical depth models
  !  integer, parameter :: maxTau = 2000000
  ! for monte carlo with optically thin model.  
  integer, parameter :: maxTau = 20000
  ! for direct integration
  !  integer, parameter :: maxTau = 8000
  !
  ! for Chris
  !  integer, parameter :: maxTau = 600000

  real :: sigmaAbs0, sigmaSca0  ! cross section at the line centre

  ! variables to do with dust

  integer :: itestlam, ismoothlam
  integer, parameter :: nMuMie = 1800
  type(PHASEMATRIX), pointer :: miePhase(:,:,:) => null()

  ! variables for clumped wind models
  
  integer, parameter :: maxBlobs = 10000
  integer :: nCurrent
  type(BLOBTYPE), allocatable :: blobs(:)
  real :: timeEnd = 24.*60.*60.
  real :: timeStart = 0.
  real :: dTime

  ! vectors

  type(VECTOR) :: originalViewVec
  ! output arrays

  integer :: iLambda
  real, allocatable :: xArray(:)

  real :: lamSmoothArray(5)

  ! model flags

  logical :: flatSpec
  logical :: ok
  logical :: greyContinuum

  ! model parameters

  real :: vel

  ! single dust blob parameters (WR137 type model)

  real :: meanDustParticleMass
  real(oct) :: t1, t2

  ! raman scattering model parameters
  type(VECTOR) :: hotSourcePosition, coolStarPosition


  ! misc
  real :: junk
  real :: r
  real :: r1
  real(oct) :: t
  integer :: i, j
  real :: nu
  real(double) :: fac
  real :: rStar
  character(len=80) :: tempChar

  !
  real(double) :: Laccretion
  real :: Taccretion, fAccretion, sAccretion
  real :: theta1, theta2 
  type(SURFACETYPE) :: starSurface



  ! binary parameters

  integer :: nVec
  type(VECTOR), allocatable :: distortionVec(:)

  real :: ang
  integer :: nt
  character(len=80) :: newContFluxFile ! modified flux file (i.e. with accretion)
  real :: infallParticleMass         ! for T Tauri infall models
  logical :: alreadyDoneInfall = .false. ! whether we have already done an infall calculation
  type(alpha_disc)  :: ttauri_disc       ! parameters for ttauri disc
  type(discwind)    :: ttauri_discwind   ! parameters for ttauri disc wind
  type(jet)         :: ttauri_jet        ! parameters for ttauri jets

  ! For romanova geometry case
  type(romanova) :: romData ! parameters and data for romanova geometry

  ! Used for multiple sources (when geometry=cluster)
  type(cluster)   :: young_cluster

  ! Name of the file to output various message from torus
  character(len=80) :: message
  real :: h 

! molecular line stuff
  type(MOLECULETYPE) :: co

  type(modelatom), allocatable :: thisAtom(:)
  integer :: nRBBTrans
  integer :: indexRBBTrans(1000), indexAtom(1000)

  real(double) :: totalmass, totalmasstrap, maxRho, minRho


#ifdef MPI
  ! For MPI implementations =====================================================
  integer ::   ierr           ! error flag
  integer ::   tempInt        !
  type(VECTOR) :: viewVec, outVec
#endif
  
! Begin executable statements --------------------------------------------------

  ! set the version number HERE!!!!!!

  call setVersion("V1.2")

  myRankGlobal = 0 ! set rank to zero for single processor job

#ifdef MPI


  ! FOR MPI IMPLEMENTATION=======================================================
  !  initialize the system for running MPI
  call MPI_INIT(ierr) 

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
  
  call init_random_seed()


  lamSmoothArray = (/5500., 1.e4, 2.e4, 5.e4, 10.e4/)
  co%molecule = " "
  indexAtom = 0
  indexRBBTrans = 0; nRBBTrans = 0
  laccretion = 0.; meanDustParticlemass = 0.
  taccretion = 0.; saccretion  = 0.; vel = 0.
  

  allocate(distortionVec(1:1))

#ifdef MPI
  call sync_random_seed()
#endif

  ! initialize

  useNdf = .true.
  sed = .false.
  jansky = .false.
  SIsed = .false.
  thinLine = .false.
  flatspec = .false.
  greyContinuum = .false.
  secondSource = .false.
  doRaman = .false.
  enhance = .false.

  contFluxFile = "none"
  intProFilename = "none"

  lucyRadiativeEq = .false. ! this has to be initialized here

  hydrodynamics = .false.
  storeScattered = .false.

  ! get the model parameters

  call inputs() ! variables are passed using the input_variables module

!  call test_profiles()  ! Testing Lorentz profile with Voigt profile

  nLambda = nLambdaInput
  sobolev = lineEmission .and. (.not.cmf)

  if(geometry == "planetgap") then
     ! takes the gapWidth and calculates what mPlanet should be
     if (planetgap) then
        call calcPlanetMass
        write (*,*) "mPlanet set to ", mPlanet
     end if
  end if

  if (molecular) then
     call readMolecule(co, moleculefile)
  endif

  if (cmf) then
     allocate(thisAtom(1:nAtom))
     do i = 1, nAtom
        call readAtom(thisAtom(i),atomFilename(i))
        call stripAtomLevels(thisAtom(i), 5)
     enddo
    call createRBBarrays(nAtom, thisAtom, nRBBtrans, indexAtom, indexRBBTrans)
  endif


  grid%photoionization = photoionization
  ! time elipsed since the beginning of the model
  grid%timeNow = phaseTime * REAL(nStartPhase-1)
  

  grid%resonanceLine = resonanceLine


  if (geometry.eq."raman") then
     hotSourcePosition = 0.75d0*secondSourcePosition
     coolStarPosition = (-0.25d0)*secondSourcePosition
     secondSourcePosition = hotSourcePosition
  endif

  if (doRaman) screened = .true.

!  if (mie .or. (geometry == "ttauri" .and. ttau_disc_on)) then
!     meanDustParticleMass = getMeanMass2(amin(1), amax(1), a0(1), qDist(1), pdist(1), graintype(1))
!     if (meanDustParticleMass <=0.0) then 
!        write(message,*) "Error: meanDustParticleMass <=0 in torusMain."
!        call writeFatal(message)
!        stop
!     end if
!  endif

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
  endif

  if (geometry == "planet") then
     flatspec = .true.
     greyContinuum = .true.
  endif

  if (geometry == "rolf") then
!     flatspec = .true.
     greyContinuum = .true.
  endif
     
  
  ! the observer's viewing direction

  originalViewVec%x = 0.
  originalViewVec%y = -sin(inclination)
  originalViewVec%z = -cos(inclination)

! Carry out geometry specific initialisation before setting up the AMR grid. 
! Applies to: cluster, molcluster, wr104, ttauri, magstream, luc_cir3d, cmfgen, romanova.
  call pre_initAMRGrid

  ! allocate the grid - this might crash out through memory problems

  if (gridUsesAMR) then
     ! any AMR allocation stuff goes here
     call initAMRGrid(newContFluxFile,flatspec,grid,ok,theta1,theta2)
     if (.not.ok) goto 666

  else
     if (gridcoords /= "polar") then
        if ((.not.doRaman).and.(geometry /= "binary"))  then
           grid = initCartesianGrid(nx, ny, nz, nLambda, lamStart, lamEnd, &
                greyContinuum, ok)
           else
           grid = initCartesianGrid(nx, ny, nz, 1, lamStart, lamEnd, &
                greyContinuum, ok)
        endif
         if (.not.ok) goto 666
     else   if (plezModelOn) then
        grid = plezModel("model.plez", lamStart, lamEnd, nLambda, kfac)
     else
        grid = initPolarGrid(nr, nmu, nphi, nlambda, lamStart, &
                             lamEnd, greyContinuum, ok)
        if (.not.ok) goto 666
     endif
  end if

!  if (mie) sed = .true.

  grid%splitOverMpi = .false.

#ifdef MPI
 if (hydrodynamics) grid%splitOverMpi = .true.
#endif

  grid%resonanceLine = resonanceLine
  grid%lambda2 = lamline

  grid%doRaman = doRaman

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


  call  createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)
  
  if (noScattering) then
     if (writeoutput) call writeWarning("Scattering opacity turned off in model")
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

  endif

  if (includeGasOpacity) then

     ! This routine may cause trouble with MPI routine! 
     ! Propably multiple nodes are trying to write to a same file at a same time.

     call readTsujiPPTable()
     call readTsujiKPTable()
     call createAllMolecularTables(20, 1.e-6, 20000., grid%nLambda, xArray)

  end if

  ! if the grid uses an adaptive mesh, create it
     

!  if (gridUsesAMR) then
!     if ((mie) .or. (geometry == "ttauri" .and. ttau_disc_on)) then 
!        call setKappaTest(grid, scale, aMin, aMax, a0, qDist, pDist, grainType, &
!             ngrain, X_grain, grainname, lambdaTau)
!     end if
!  end if


  if (mie) then
     call locate(grid%lamArray, nLambda,lambdaTau,itestlam)
     call locate(grid%lamArray, nLambda,lambdasmooth,ismoothlam)
     grid%itestlam = iTestLam
  end if

  ! chris
  if (hydroWarp) then
     print *, "Creating and reading hydroGrid..."
     allocate(grid%hydroGrid)
     call initAMRGrid(newContFluxFile,flatspec,grid%hydroGrid,ok,theta1,theta2)
     call readAMRgrid("hydrogrid.dat",readFileFormatted,grid%hydroGrid)
     print *, "hydroGrid read..."
     call hydroWarpFitSplines(grid%hydroGrid)
     print *, "splines fitted..."
     print *, "hydroGrid setup done."
  end if


  !===============================================================
  if (gridUsesAMR) then  !========================================
  !===============================================================

     if (cmf.and.readLucy) then
        call readAMRgrid("atom_tmp.grid",.false.,grid)
     else if (photoionization.and.readlucy) then
        continue
     elseif (molecular .and. (readmol .or. addnewmoldata)) then
        continue
     else if (.not. restart ) then 
        ! Set up the AMR grid
        call amr_grid_setup
        if(molecular .and. writeoutput) call writeAMRgrid('notmolecular.grid',writeFileFormatted,grid) 

        if (geometry(1:7) .eq. "testamr") call testamr
        
        if (grid%geometry(1:8) == 'windtest') call windtest

#ifdef MPI
        if (grid%splitOverMPI) then
           call grid_info_mpi(grid, "info_grid.dat")
        else
           if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
        endif
#else
        if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
#endif

        ! Plotting the various values stored in the AMR grid.
        if ( plot_maps .and. myRankIsZero ) call writeVtkFile(grid, "rho.vtk")

        ! Output H 21cm emissivity and opacity 
        if ( h21cm ) call writeVtkFile(grid, "h21cm.vtk", valueTypeString=(/"etaline","chiline"/) )

     end if

  !=================================================================
  else ! grid is not adaptive ======================================     
  !=================================================================

     ! fill up the grid with the appropriate opacities
     call fill_opacities_noamr( ok )
     if ( .not. ok ) goto 666

     ! Distort grid if required
     call distort_noamr

     if ( enhance ) then
        nVec = 100
        nPhi = 360
	if (allocated(distortionVec)) deallocate(distortionVec)
        allocate(distortionVec(nVec))
        timeStart = 0.
        timeEnd = 12.*3600.
        dTime = (timeEnd - timeStart)/real(nPhase)
        call initInfallEnhancement(distortionVec, nVec, nPhi, infallParticleMass)
        do i = 1, nStartPhase-1
           call infallEnhancment(grid, distortionVec, nVec, nPhi, dTime, &
                .false., infallParticleMass, alreadyDoneInfall)
        enddo
     endif

     ! Write the information on the grid to file using a routine in grid_mod.f90
#ifdef MPI
        if (grid%splitOverMPI) then
           call grid_info_mpi(grid, "info_grid.dat")
        else
           if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
        endif
#else
        if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
#endif

  !=================================================================
  endif ! (gridusesAMR)
  !=================================================================

  ! set up the sources
  if(.not. (restart .or. addnewmoldata)) call set_up_sources

! Any tasks to be performed after the AMR grid is set up should be done here.
! Calculates the mass on thr grid and also runs some geometry specific code. 
  call post_initAMRgrid

     call init_random_seed()
     if (cmf) then
        if (.not.readlucy) call atomLoop(grid, nAtom, thisAtom, nsource, source)

        do i = 1, 1
           ang = twoPi * real(i-1)/51. 
!           ang = 45.*degtorad
           t = 120.d0*degtorad
           call calculateAtomSpectrum(grid, thisAtom, nAtom, iTransAtom, iTransLine, &
                VECTOR(sin(t)*cos(ang), sin(t)*sin(ang), cos(t)), 100.d0*pctocm/1.d10, &
                source, nsource,i)
        enddo
        stop
     endif

#ifdef MPI
     if (photoIonization.and.hydrodynamics) then
        grid%splitOverMPI = .true.
        grid%photoionization = .true.

       inclination = 77.5d0 * degToRad
       viewVec%x = 0.
       viewVec%y = -sin(inclination)
       viewVec%z = -cos(inclination)
       viewVec = rotateZ(viewVec, 30.d0*degtorad)
       outVec = (-1.d0)*viewVec

       call createImageSplitGrid(grid, nSource, source, outVec, i, fac)
       call torus_mpi_barrier
       stop

           call radiationHydro(grid, source, nSource, nLambda, xArray, readlucy, writelucy, &
              lucyfilenameout, lucyfilenamein)
           call torus_mpi_barrier
     endif

     if (hydrodynamics) then

        grid%splitOverMPI = .true.
	
        if (idump /= 1) then
           call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
           write(message,'(a,i4.4,a)') "dump",idump,".grid"
           write(*,*) myrankglobal, " calling read ",trim(message)
           call readAMRgrid(message,.false.,grid)
           write(*,*) myrankglobal, " done reading ",trim(message)
        endif


        if (grid%octreeRoot%twoD) then
           call doHydrodynamics2d(grid)
           goto 666
        else if (grid%octreeRoot%oneD) then
           call doHydrodynamics1d(grid)
           goto 666
        else if (grid%octreeRoot%threeD) then
           call doHydrodynamics3d(grid)
           goto 666
        endif

     endif
#endif

  if (photoionization) then 
        call photoIonizationloop(grid, source, nSource, nLambda, xArray, readlucy, writelucy, &
             lucyfileNameout, lucyfileNamein)
        goto 666
  end if

  if (lucyRadiativeEq) call do_lucyRadiativeEq

  if (molecular) then
     ! Plotting the various values stored in the AMR grid.
     if ( geometry .eq. "molcluster" .and. plot_maps .and. myRankIsZero .and. .not. (restart .or. addnewmoldata .or. readmol)) &
          call writeVtkFile(grid, "rho.vtk", "test.txt")

     if (writemol) call molecularLoop(grid, co)
     
     if(.not. writemol) then
        if (molecular .and. readmol .and. .not. (lucyRadiativeEq)) then
           if(openlucy) then
              call readAMRgrid(lucyFileNamein,.false.,grid)
           else
              if(lte) then
                 write(molfilenamein,*) trim(co%molecule),"_lte.grid"
              else
                 write(molfilenamein,*) trim(co%molecule),"_grid.grid"
              endif
              call readAMRgrid(molfilenamein,.false.,grid)
           endif
        endif
     endif

     if (readmol) then 
        call findTotalMass(grid%octreeRoot, totalMass, totalmasstrap = totalmasstrap, maxrho=maxrho, minrho=minrho)
        write(message,*) "Mass of envelope: ",totalMass/mSol, " solar masses"
           call writeInfo(message, TRIVIAL)
        if(geometry .eq. 'molcluster') then
           write(message,*) "Mass of envelope (TRAP): ",totalMasstrap/mSol, " solar masses"
           call writeInfo(message, TRIVIAL)
           write(message,*) "Maximum Density: ",maxrho, " g/cm^3"
           call writeInfo(message, TRIVIAL)
           write(message,*) "Minimum Density: ",minrho, " g/cm^3"
           call writeInfo(message, TRIVIAL)
        endif
        call calculateMoleculeSpectrum(grid, co)
     endif

     goto 666

  endif

! Generate H 21cm image
  if ( h21cm ) then 

     if (restart) then 
        call readAMRgrid("h21cm_grid",.false.,grid)
     else
        call writeAMRgrid("h21cm_grid",.false.,grid)
     endif

     if ( internalView ) then 
        call make_angular_image(grid)
     else
        call make_h21cm_image(grid)
     end if

     goto 666
  end if

!  if (grid%geometry == "shakara") call defineDiffusionZone(grid, .false., .false.)

  if (useBias .and. .not. formalsol) call set_emission_bias

  if (myRankIsZero .and. geometry == "shakara") call polardump(grid)
  
  ! initialize the blobs if required

  if (nBlobs > 0) call initialize_blobs


  if (nPhase /= 1) then
     dTime = (timeEnd - timeStart) / real(nPhase-1)
  endif


  lucyRadiativeEq = .false.
  call set_up_lambda_array
  
  if (geometry=="pathtest") then
     call setupNeighbourPointers(grid, grid%octreeRoot)
     call pathTest(grid)
     goto 666
  endif

  if ( nInclination > 0 ) then
     call do_phaseloop(grid, alreadyDoneInfall, meanDustParticleMass, rstar, vel, &
          theta1, theta2, coolstarposition, Laccretion, Taccretion, fAccretion, sAccretion, corecontinuumflux, &
          starsurface, newContFluxFile, sigmaAbs0, sigmaSca0, ttauri_disc, distortionVec, nvec,       &
          infallParticleMass, maxBlobs, flatspec, inclination, maxTau, &
          miePhase, nsource, source, blobs, nmumie, dTime)
  end if

! Tidy up and finish the run 

666 continue

#ifdef MPI
write(*,*) myrankGlobal , " is waiting to finish"
call torus_mpi_barrier
#endif

if (doTuning) call tune(6, "Torus Main") ! stop a stopwatch  

call writeInfo("TORUS exiting", FORINFO)

call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
call freeGrid(grid)

deallocate(miePhase) 
deallocate(xArray)

call torus_mpi_barrier
#ifdef MPI
call MPI_FINALIZE(ierr)
#endif

CONTAINS
!-----------------------------------------------------------------------------------------------------------------------

  subroutine pre_initAMRGrid

    use input_variables, only: sphdatafilename
    use isochrone_class, only: isochrone, read_isochrone_data, new
    use cluster_class, only: write_catalog, build_cluster, new
    use romanova_class, only: get_dble_parameter, new
    use sph_data_class, only: find_inclinations, new_read_sph_data, read_stellar_disc_data, info, read_galaxy_sph_data
    use wr104_mod, only: readwr104particles
    use cmfgen_class, only: read_cmfgen_data, put_cmfgen_Rmin, put_cmfgen_Rmax
    use luc_cir3d_class, only: new

    type(isochrone) :: isochrone_data
    real(double)    :: objectDistance

  !=====================================================================
  ! INIIALIZATIONS BEFORE CALLING INITAMR ROUTINE SHOULD BE HERE
  !=====================================================================

  if (geometry == "cluster") then

     ! HSC
     ! read in the sph data from a file
     call read_sph_data(sphdata,"sph.dat")
     call read_stellar_disc_data(sphdata,"stellar_disc.dat")

        ! Writing basic info of this data
     if (myRankIsZero) call info("info_sph.dat")

     ! reading in the isochrone data needed to build an cluster object.
     call new(isochrone_data, "dam98_0225")   
     call read_isochrone_data(isochrone_data)
     
     ! making a cluster object
     call new(young_cluster, dble(amrGridSize), disc_on)
     call build_cluster(young_cluster, dble(lamstart), dble(lamend), iso_data=isochrone_data)
    
     ! Wrting the stellar catalog readble for a human
     if (myRankIsZero) call write_catalog(young_cluster)

      ! Finding the inclinations of discs seen from +z directions...
     if (myRankIsZero) call find_inclinations(0.0d0, 0.0d0, 1.0d0, "inclinations_z.dat")

  elseif (geometry .eq. "molcluster" .and. .not. readmol .and. .not. restart .and. .not. addnewmoldata) then

     ! read in the sph data from a file

     call new_read_sph_data(sphdatafilename)

!     call read_stellar_disc_data(sphData, "stellar_disc.dat")

     ! Writing basic info of this data
     if (myRankIsZero) call info("info_sph.dat")

!     ! reading in the isochrone data needed to build an cluster object.
     call new(isochrone_data, "dam98_0225")   
     call read_isochrone_data(isochrone_data)
     
     ! making a cluster object

     call new(young_cluster, dble(amrGridSize), disc_on)

     call build_cluster(young_cluster, dble(lamstart), dble(lamend), iso_data=isochrone_data)

     ! Wrting the stellar catalog readble for a human
     if (myRankIsZero) call write_catalog(young_cluster)

      ! Finding the inclinations of discs seen from +z directions...
!     if (myRankIsZero) call find_inclinations(sphData, 0.0d0, 0.0d0, 1.0d0, "inclinations_z.dat")

  elseif (geometry == "wr104") then
     if (.not.(readPops.or.readlucy)) then
        objectDistance = griddistance * pctocm
        call readWR104Particles("harries_wr104.txt",sphdata , objectDistance)
        call info("*")
     endif

  elseif (geometry == "ttauri".or.(geometry == "magstream")) then
     onekappa=.false.
     if (ttau_disc_on)  &
          call new(ttauri_disc, dble(TTauriDiskRin*TTauriRstar/1.0e10), 1.5d3*100.d0, &
          dble(TTauriMstar/100.0), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, dble(TTauriMstar/mSol))
     if (ttau_discwind_on) &
          call new(ttauri_discwind, DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
          DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, &
          dble(TTauriMstar/mSol),  dble(TTauriDiskHeight/1.0e10) )
     if (ttau_jet_on) &
!          call new(ttauri_jet,  dble(TTauriRouter/1.0e10), dble(amrgridsize), &
          call new(ttauri_jet,  JET_Rmin, dble(amrgridsize), &
          JET_theta_j, dble(TTauriMstar/mSol), JET_Mdot, JET_a_param, JET_b_param,  &
          JET_Vbase, JET_Vinf, JET_beta, JET_gamma, limitscalar2, JET_T)

  elseif (geometry == "luc_cir3d") then
     onekappa=.false.
     !                    Rmin [10^10cm]         
     call new(CIR_Rstar*Rsol*1.0e-10, CIR_Twind, CIR_Mdot_scale,  &
     "zeus_rgrid.dat", "zeus_rho_vel.dat")     

  elseif (geometry == "cmfgen") then
     onekappa=.false.
     call read_cmfgen_data("OPACITY_DATA")
     call put_cmfgen_Rmin(CMFGEN_Rmin)  ! in [10^10cm]
     call put_cmfgen_Rmax(dble(amrGridSize)/2.0d0)  ! in [10^10cm]

  elseif (geometry == "romanova") then
     onekappa=.false.
     ! create the romanova data object and readin the data at the same time.
     ! After this call, the obejct will be ready to use. 
     call new(romData, ROM_Rs*DBLE(Rsol)*1.0d-10, 70.0d0, &
          ROM_Mass*DBLE(Msol), ROM_isoT, ROM_T_flow, ROM_r_ref, &
          ROM_rho_ref, ROM_T_ref, ROM_v_ref, ROM_tilt, ROM_period, ROM_datafile)

     ! this is a temporary solution... It should be called something 
     ! else .
     TTauriRouter = REAL(get_dble_parameter(romData, "Rmax")*1.0d10)           

  elseif (geometry == "theGalaxy" ) then
     call read_galaxy_sph_data(sphdatafilename)

  end if

end subroutine pre_initAMRGrid

!-----------------------------------------------------------------------------------------------------------------------

  subroutine set_up_lambda_array

    use photoion_mod, only: refineLambdaArray

    real :: deltaLambda
    real :: loglamStart, logLamEnd
    real, allocatable :: tArray(:)
    real :: testlam

    if (allocated(xArray)) then
       deallocate(xArray)
    endif

     if (lucyradiativeEq) then
        call writeInfo("Doing radiative equilibrium so setting own wavelength arrays", TRIVIAL)
        nLambda = 200
        allocate(xArray(1:nLambda))
        logLamStart = log10(lamstart)
        logLamEnd   = log10(lamend)
        do i = 1, nLambda
           xArray(i) = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
           xArray(i) = 10.**xArray(i)
        enddo
        goto 777
     endif



     if(nLambdaInput == 0)  then
        if (allocated(source)) then
           call writeInfo("Basing SED wavelength grid on input photospheric spectrum",TRIVIAL)
           nLambda = SIZE(source(1)%spectrum%lambda)
           allocate(xArray(1:nLambda))
           xArray = source(1)%spectrum%lambda
           
           
           lamStart = 1200.
           lamEnd = 2.e7
           logLamStart = log10(lamStart)
           logLamEnd   = log10(lamEnd)
           
           do i = 1, 200
              testLam = logLamStart + real(i-1)/real(200-1)*(logLamEnd - logLamStart)
              if (testLam < xArray(1)) then
                 allocate(tArray(1:nLambda))
                 tArray = xArray
                 nLambda = nLambda + 1
                 deallocate(xArray)
                 allocate(xArray(1:nLambda))
                 xArray(1) = testLam
                 xArray(2:nLambda) = tArray
                 deallocate(tArray)
              endif
              if (testLam > xArray(nLambda)) then
                 allocate(tArray(1:nLambda))
                 tArray = xArray
                 nLambda = nLambda + 1
                 deallocate(xArray)
                 allocate(xArray(1:nLambda))
                 xArray(nLambda) = testLam
                 xArray(1:nLambda-1) = tArray
                 deallocate(tArray)
              endif
           enddo
        else
           nLambda = 200
           allocate(xArray(1:nLambda))
           lamStart = 1200.
           lamEnd = 2.e7
           logLamStart = log10(lamStart)
           logLamEnd   = log10(lamEnd)
           do i = 1, nLambda
              xArray(i) = logLamStart + real(i-1)/real(nLambda-1)*(logLamEnd - logLamStart)
              xArray(i) = 10.**xArray(i)
           enddo
        endif
     else
        nLambda = nLambdaInput
        allocate(xArray(1:nLambda))
        
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
           
           if (photoionization) then
              xArray(1) = lamStart
              xArray(2) = lamEnd
              nCurrent = 2
              call refineLambdaArray(xArray, nCurrent, grid)
              nt = nLambda - nCurrent
              do i = 1, nt
                 fac = logLamStart + real(i)/real(nt+1)*(logLamEnd - logLamStart)
                 fac = 10.**fac
                 nCurrent=nCurrent + 1
                 xArray(nCurrent) = fac
                 call sort(nCurrent, xArray)
              enddo
           endif

        endif


!       if (mie) then
!          if ((lambdaTau > Xarray(1)).and.(lambdaTau < xArray(nLambda))) then
!             call locate(xArray, nLambda, lambdaTau, i)
!             t1 = (lambdaTau - xArray(i))/(xArray(i+1)-xArray(i))
!             if (t1 > 0.5) then
!                write(message,*) "Replacing ",xArray(i+1), " wavelength step with ",lambdaTau
!                call writeInfo(message, TRIVIAL)
!                xArray(i+1) = lambdaTau
!             else
!                write(message,*) "Replacing ",xArray(i), " wavelength step with ",lambdaTau
!                call writeInfo(message, TRIVIAL)
!                xArray(i) = lambdaTau
!             endif
!          endif
!       endif

     endif


       if (lamFile) then
          call writeInfo("Reading wavelength points from file.", TRIVIAL)
          open(77, file=lamfilename, status="old", form="formatted")
          nLambda = 1
333       continue
          read(77,*,end=334) junk
          xArray(nLambda) = junk
          nLambda = nLambda + 1
          goto 333
334       continue
          nlambda = nlambda - 1
          close(77)
       endif

777 continue
    !
    ! Copying the wavelength array to the grid
       if (associated(grid%lamArray)) deallocate(grid%lamArray)
       allocate(grid%lamArray(1:nLambda))
    do i = 1, nLambda
       grid%lamArray(i) = xArray(i)
    enddo
    grid%nLambda = nLambda
  end subroutine set_up_lambda_array

!-----------------------------------------------------------------------------------------------------------------------

  subroutine windtest

    use stateq_mod, only: BoltzSaha
    use gridtype_mod, only: statEqMaxLevels 

    type(VECTOR) :: octVec
    integer, parameter :: nrGrid = 1000
    real :: rGrid(nrGrid)
    real,dimension(statEqMAxLevels) :: meanDepart ! for testing
    real :: treal
!    type(VECTOR) :: rVec 
    real(double)  :: Ne         ! for testing
    real(double),dimension(statEqMAxLevels) :: levelPops  ! for testing
    integer                         :: level      ! for testing

    do i = 1, 1000
       r = log10(grid%rInner) + (log10(grid%rOuter)-log10(grid%rInner))* real(i-1)/999.
       r = 10.**r
       r1 = log10(grid%rInner) + (log10(grid%rOuter)-log10(grid%rInner))* real(i)/999.
       r1 = 10.**r1
       rGrid(i) = r
    enddo
    
    open(21,file="rDepart.dat",status="unknown",form="formatted")

    do i = 1, 1000
       meanDepart = 0.
       nt = 0
!       outVec = VECTOR(1., 0., 0.)
!       rVec = dble(rGrid(i)) * outVec 
       r = rgrid(i)
!       call integratePathAMR(10.e4,  lamLine, VECTOR(1.,1.,1.), &
!            rVec, outVec, grid, lambda, &
!            tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
!            .false., lamStart, lamEnd, nLambda, contTau, &
!            hitCore, thinLine,.false.,  &
!            .false., nUpper, nLower, 0., 0., 0., junk,sampleFreq,intPathError)

       do j = 1, 100
          ang = twoPi * real(j-1)/100.
          octVec = VECTOR(r*cos(ang), r*sin(ang),0.)
          call amrGridValues(grid%octreeRoot, octVec, temperature=tReal, &
              ilambda=1, N=levelPops, Ne=Ne)
          t1 = dble(treal)
          do level = 1, statEqMaxLevels, 1
             meanDepart(level) = meanDepart(level) + real(levelPops(level))/boltzSaha(level,Ne,dble(treal))
          end do
          nt = nt + 1
       enddo
       if (nt > 0) then
          if (writeoutput) write(21,*) r/grid%rInner,meanDepart/real(nt)
       endif
    Enddo
    close(21)

  end subroutine windtest

!-----------------------------------------------------------------------------------------------------------------------

  subroutine testamr

    type(VECTOR) :: octVec
    real :: meant, meaneta
    integer, parameter :: nrGrid = 1000
    real :: rGrid(nrGrid), drGrid(nrgrid)
    real :: treal
!    type(VECTOR) :: rVec 
    real(double) :: kabs, eta

    kabs = 0.

    do i = 1, 1000
       meant = 0.
       meaneta =0.
       nt = 0
       r = log10(grid%rInner) + (log10(grid%rOuter)-log10(grid%rInner))* real(i-1)/999.
       r = 10.**r
       r1 = log10(grid%rInner) + (log10(grid%rOuter)-log10(grid%rInner))* real(i)/999.
       r1 = 10.**r1
       rGrid(i) = r
       drGrid(i) = r1 - r
    enddo

    open(21,file="r.dat",status="unknown",form="formatted")
    call locate(xarray, nLambda, 10.e4, ilambda)

    do i = 1, 1000
       meant = 0.
       meaneta =0.
       nt = 0
!       outVec = VECTOR(1., 0., 0.)
 !      rVec = dble(rGrid(i)) * outVec 
       r = rgrid(i)
!       call integratePathAMR(10.e4,  lamLine, VECTOR(1.,1.,1.), &
!            rVec, outVec, grid, lambda, &
!            tauExt, tauAbs, tauSca, maxTau, nTau, opaqueCore, escProb, &
!            .false., lamStart, lamEnd, nLambda, contTau, &
!            hitCore, thinLine,.false.,  &
!            .false., nUpper, nLower, 0., 0., 0., junk,sampleFreq,intPathError)



       do j = 1, 100
          ang = twoPi * real(j-1)/100.
          octVec = VECTOR(r*cos(ang), r*sin(ang),0.)
          call amrGridValues(grid%octreeRoot, octVec, temperature=treal, &
              ilambda=1, grid=grid ,etacont = eta, kappaAbs=kabs)
          if (treal > 1.) then
             t1 = dble(treal)
             meant = meant + t1
             meaneta = meaneta + eta
             nt = nt + 1
          endif
       enddo
       if (nt > 0) then
          if (writeoutput) write(21,*) r/grid%rInner,meant/real(nt), meaneta * r**2 * drgrid(i)/ real(nt)
       endif
    Enddo
    close(21)

  end subroutine testamr

!-----------------------------------------------------------------------------------------------------------------------

  subroutine fill_opacities_noamr( ok )

    use puls_mod, only: fillGridPuls
    use grid_mod
    use TTauri_mod, only: fillGridFlaredDisk, fillgridmagneticaccretion, fillgridttauriwind
    use stateq_mod, only: initgridstateq

    logical, intent(out) :: ok

         ! fill up the grid with the appropriate opacities
    ok = .true.
    if (.not.plezModelOn) then
       select case(geometry)
       case("torus")
          call fillGridTorus(grid, rho, rTorus, rOuter)
       case("sphere")
          call fillGridSpheriod(grid, rho, radius, kFac)
       case("flared")
          call fillGridFlaredDisk(grid, meanDustParticleMass)
       case("ellipse")
          call fillGridEllipse(Grid,rho, rMin, rMaj, rinner, teff)
       case("disk")
          call fillGridDisk(grid, rho, rCore, rInner, rOuter, height, mCore, diskTemp)
       case("star")
          call fillGridStar(grid, radius, mdot, vel, kfac, scale)
       case("spiral")
          call fillGridSpiral(grid, radius, mdot, vel, scale)
       case("shell")
          call fillGridShell(grid, radius, shellFrac, rho, kfac)
       case("stateq")
          call fillGridStateq(grid, opacityDataFile, kfac, scaleDensity)
       case("bipolar")
          call fillGridBipolar(grid, rho,  30.)
       case("collide")
          call fillGridCollide(grid, rho, momRatio, binarySep, mie, meanDustParticleMass, logMassLossRate)
       case("dustblob")
          call fillGridDustBlob(grid, dustBlobdistance, rho, phiDustBlob, &
               xDustBlobSize, yDustBlobSize, zDustBlobSize)
       case("raman")
          call  fillGridRaman(grid, 10.*rStar, mdot, vterm, rStar, coolStarPosition, beta)
       case("binary")
          call fillGridBinary(grid, radius1, radius2, mass1, mass2, &
               temp1, temp2, vNought1, vNought2, vTerm1, vTerm2, beta1, beta2, period, &
               mdot1, mdot2, deflectionAngle, shockWidth, shockFac)
          
          call initGridStateq(grid, contFluxFile, contFluxFile2, &
               popFilename, readPops, writePops, lte, nLower, nUpper)
        
       case("rolf")
          call fillGridRolf(grid, mdot, vterm, o6width)
       case("wr137")
          call fillGridWR137(grid, rCore, mDot, vTerm, beta, temp1)
       case("planet")
          call fillGridPlanet(grid)
        case("hourglass")
           call fillGridHourglass(grid)
        case("ttauri")
           call fillGridMagneticAccretion(grid,contfluxfile, popFileName, &
                readPops, writePops, lte,  lamLine, Laccretion, Taccretion, sAccretion, &
                curtains, dipoleOffset, nLower, nUpper, theta1, theta2)
        case("ttwind")
           call fillGridTTauriWind(grid,contfluxfile, popFileName, &
                readPops, writePops, lte, nLower, nUpper)
        case("betacep")
           call fillGridBetaCep(grid)
        case("donati")
           resonanceLine = .false.
           if (lamLine < 2000.) resonanceLine = .true.
           call fillGridDonati2(grid, resonanceLine, misc)
           if (.not.resonanceLine) then
              call initGridStateq(grid, contFluxFile, contFluxFile2, popFilename, &
                   readPops, writePops, lte, nLower, nUpper)
           endif
           grid%kappaAbs = 1.e-30

        case("puls")
           call fillGridPuls(grid, mDot, rcore, tEff, v0, vterm, beta, xfac, blobs, maxBlobs, .false., vContrast)
   !        call initGridStateq(grid, contFluxFile, contFluxFile2, popFilename, &
   !                            readPops, writePops, lte, nLower, nUpper)
        case("wind")
!           call fillGridWind(grid, mDot, rStar, tEff, v0, vterm, beta, &
!           lte, contFluxFile, writePops, readPops, popFilename, nLower, nUpper)

           call initgridstateq(grid, contfluxFile, " ", popFileName, &
                readPops, writePops, lte, nLower, nUpper)
           grid%etaCont = 1.e-30
           grid%kappaAbs = 1.e-30
        
        case("resonance")
           call fillGridResonance(grid, rCore, mDot, vTerm, beta, temp1)

        case DEFAULT
           if (writeoutput) write(*,*) "! Unrecognised grid geometry in fill_opacities_noamr: ",trim(geometry)
           ok = .false.
           return
        end select
     endif ! (.not.plezModelOn)

     if (fillTio) then
        call fillGridTio(grid, scale)
     endif

     if (fillRayleighOpacity) then
        call fillGridRayleigh(grid,scale)
     endif 

     if (fillThomson) then
        call fillGridThomson(grid)
     endif


  end subroutine fill_opacities_noamr

!-----------------------------------------------------------------------------------------------------------------------

  subroutine distort_noamr

    use distortion_mod, only: distortgridspiral, distortgridtest, distortrotation, distortstrom, &
         distortwindcollision, distortwrdisk

  ! the distortion types

     select case(distortionType)
     case("spiral")
        call distortGridSpiral(grid, vRot, nSpiral)
     case("rotation")
        call distortRotation(grid, vRot)
        if (writeoutput) write(*,'(a,f5.1,a)') "Grid distorted by a rotational velocity of ", &
                               vRot/1.e5," km/s"
     case("test")
        call distortGridTest(grid)
!     case("raman")
!        call distortRaman(grid)
     case("raman")
        call  distortStrom(grid, secondSourcePosition, .true., .true., 0.5*rStar, coolStarPosition, ramanDist)

     case("wrdisk")
        call distortWRdisk(grid)
     case("windwind")
        call distortWindCollision(grid, momRatio, binarySep)

     end select

  end subroutine distort_noamr

!-----------------------------------------------------------------------------------------------------------------------

  subroutine amr_grid_setup

    use amr_mod
    use spectrum_mod, only: fillSpectrumBB, normalizedSpectrum
    use stateq_mod, only: map_cmfgen_opacities, amrStateq, generateOpacitiesAMR
    use dust_mod, only: filldustshakara, filldustuniform, filldustwhitney
    use cluster_class, only: remove_too_close_cells
    use surface_mod, only: createTTauriSurface, createTTauriSurface2, createSurface
    use luc_cir3d_class, only: deallocate_zeus_data
    use lucy_mod, only: lucyRadiativeEquilibrium, lucyRadiativeEquilibriumAMR, allocateMemoryForLucy, putTau

    type(VECTOR) :: amrGridCentre ! central coordinates of grid
    real(double) :: mass_scale, mass_accretion_old, mass_accretion_new
    real(double) :: removedMass
    real         :: sigmaExt0
    logical      :: gridConverged  ! true when adaptive grid structure has 
                                   !   been finalised
    character(len=80) :: phasePopFilename

    type(STREAMTYPE)  :: thisStream(2000), bigStream
    integer           :: nStreams

  ! adaptive grid stuff
    integer           :: nOctals       ! number of octals in grid
    integer           :: nVoxels       ! number of unique voxels in grid
                                       ! (i.e. the number of childless subcells)
    type(VECTOR) :: dummy

    if (doTuning) call tune(6, "AMR grid construction.")  ! start a stopwatch

    nStreams = 0
    
    if (readPops .or. readPhasePops .or. readLucy) then 


       if (readLucy) call readAMRgrid(lucyFilenameIn,readFileFormatted,grid)

       if (readPhasePops) then ! need to get the right file for the phase
          write(tempChar,'(i3.3)') nStartPhase
          phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
          call readAMRgrid(phasePopFilename,readFileFormatted,grid)
       else ! just read the normal pops file
          if (readpops) call readAMRgrid(popFilename,readFileFormatted,grid)
       end if

       if (forceLineChange) then
          write (message,'(A,I1,A,I1,A)') 'Recalculating for n =',nUpper,' /',nLower,'levels...'
          call writeInfo(message, IMPORTANT)
          call generateOpacitiesAMR(grid, nLower, nUpper)
          do i = 1, nLambda
             grid%lamArray(i) = xArray(i)
          enddo
          print *, '...level change done.'
       end if

       ! for some geometries, it is valid to alter the 'dipoleOffset parameter
       !   for a grid we've already created. 
       if (geometry == "ttauri".or.(geometry == "magstream")) then 
          ! if (ABS(grid%dipoleOffset/dipoleoffset-1.0) > 0.01) then 
          write(message,'(a,f5.2,a)') 'Using new dipole offset value (',&
               dipoleOffset*radToDeg,'deg)'
          call writeInfo(message, IMPORTANT)
          grid%dipoleOffset = dipoleOffset                         
          grid%diskNormal = VECTOR(0.,0.,1.)
          grid%diskNormal = rotateX(grid%diskNormal,dble(grid%dipoleOffSet))
!          end if
       elseif (geometry == "romanova") then
!           write(*,'(a,f5.2,a)') 'Using new dipole offset value (',&
!                dipoleOffset*radToDeg,'deg)'
!           grid%dipoleOffset = dipoleOffset                         
          grid%diskNormal = VECTOR(0.,0.,1.)
!           grid%diskNormal = rotateY(grid%diskNormal,-grid%dipoleOffSet)
!           grid%diskNormal = VECTOR(0.,0.,1.)
       end if

        ! In case, a user is using a new nlambda (number of wavelength bins)
        ! we have to update the wavelength array in the grid.
        ! This may not work for dust calculation...
        ! Copying the wavelength array to the grid
!        if (grid%nlambda /= nlambda) then
!           if (mie .or. (grid%geometry == "ttauri" .and. ttau_disc_on)) then              
!              write(*,*) "Error:: The number of the wavelength bins in the data file "
!              write(*,*) "        does not macth with that specifed in your parameter file."
!              write(*,*) "nlambda(old) = ", grid%nlambda
!              write(*,*) "nlambda(new) = ", nlambda
!              write(*,*) " "
!              write(*,*) "Make the new nlambda same as the old one, otherwise you cannot use"
!              write(*,*) "this data file. This should be changed in future...."
!!              stop
!           else 
!              write(*,*) "Warning:: The number of the wavelength bins in the data file "
!              write(*,*) "      does not macth with that specifed in your parameter file."
!              write(*,*) "nlambda(old) = ", grid%nlambda
!              write(*,*) "nlambda(new) = ", nlambda
!              write(*,*) "==> We recompute the wavelength array with new nlambda value, and continue."
              ! update the wavelength array.
              ! It over writes whatever the wavelength you had in your file!
       deallocate(grid%lamArray)
       allocate(grid%lamArray(nlambda))
       do i = 1, nLambda
          grid%lamArray(i) = xArray(i)
       enddo
       grid%nLambda = nLambda
!           end if
!        end if

        !
        ! If the grid read from file contains ttauri disc or jet, you can turn it
        ! off by setting ttau_trun_off_disc = .true. and/or ttaur_turn_off_jet
        ! in your parameter file.
       if (grid%geometry=="ttauri" .and. ttau_turn_off_disc) then
          write(*,*) " "
          write(*,*) "Turning off the alpha disc read in from file, but keeping"
          write(*,*) "the magnetorsphere alive."
          if (.not.ttau_disc_on)   then
              ! we need to create the disc parameter object 
             call new(ttauri_disc, dble(TTauriDiskRin*TTauriRstar/1.0e10), 1.5d3*100.d0, &
                  dble(TTauriMstar/100.0), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, dble(TTauriMstar/mSol))
          end if
          call turn_off_disc(grid%octreeroot, grid, ttauri_disc)
       end if

       if (grid%geometry=="ttauri" .and. ttau_turn_off_jet) then
          write(*,*) " "
          write(*,*) "Turning off the jet read in from file, but keeping"
          write(*,*) "the magnetorsphere alive."
          if (.not.ttau_jet_on) then
              ! the parameter file has to created here.
             call new(ttauri_jet,  dble(TTauriRouter/1.0e10), dble(amrgridsize), &
                  JET_theta_j, dble(TTauriMstar/mSol), JET_Mdot, JET_a_param, JET_b_param,  &
                  JET_Vbase, JET_Vinf, JET_beta, JET_gamma, limitscalar, JET_T)
          end if
          call turn_off_jet(grid%octreeroot, grid, ttauri_jet)
       end if

       if (grid%geometry=="ttauri" .and. ttau_turn_off_acc) then
          write(*,*) " "
          write(*,*) "Turning off the magnetosphere read in from file."
          write(*,*) " "
           ! using the routine in amr_mod.f90
          call turn_off_magnetosphere(grid%octreeroot, grid, dble(TTauriRouter/1.0e10))
       end if

       if (myRankIsZero) then   !-----------------------------------------------
          if (writePhasePops) then
             write(tempChar,'(i3.3)') nStartPhase
             phasePopFilename = trim(popFilename)//'_phase'//TRIM(tempChar)
             call writeAMRgrid(phasePopFilename,writeFileFormatted,grid)
          else if (writePops) then
             call writeAMRgrid(popFilename,writeFileFormatted,grid)
          end if
       end if  ! (myRankIsZero) --------------------------------------------------------

       call torus_mpi_barrier

    else  ! not reading a population file

       amrGridCentre = VECTOR(amrGridCentreX,amrGridCentreY,amrGridCentreZ)
       call writeInfo("Starting initial set up of adaptive grid...", TRIVIAL)
       
       select case (geometry)
       case("cluster", "theGalaxy")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, young_cluster, nDustType)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, young_cluster)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
!          call estimateRhoOfEmpty(grid, grid%octreeRoot)	
           !Removing the cells within 10^14 cm from the stars.
          removedMass = 0.0
          call remove_too_close_cells(young_cluster,grid%octreeRoot,1.0d4, removedMass, amr_min_rho, 's')
          write(message,*) "Mass removed by remove_too_close_cells= ", removedMass / mSol
          call writeInfo(message, TRIVIAL)

       case("molcluster")
          if(.not. readmol) then
             call writeInfo("Initialising adaptive grid...", TRIVIAL)
             call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, young_cluster, nDustType)
             call writeInfo("Done. Splitting grid...", TRIVIAL)
             call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, young_cluster)
             call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
!          call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
!          do
!             gridConverged = .true.
!             call myScaleSmooth(smoothFactor, grid%octreeRoot, grid, &
!                  gridConverged,  inheritProps = .false., interpProps = .false., &
!                  sphData=sphData, stellar_cluster=young_cluster, romData=romData)
!             if (gridConverged) exit
!          end do
!          call writeInfo("...grid smoothing complete", TRIVIAL)
          endif
       case("wr104")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, young_cluster, nDustType)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
          call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
          do
             gridConverged = .true.
             call myScaleSmooth(smoothFactor, grid, &
                  gridConverged,  inheritProps = .false., interpProps = .false.)
             if (gridConverged) exit
          end do
          call writeInfo("...grid smoothing complete", TRIVIAL)
          
       case("starburst")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, young_cluster, nDustType)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid)
          gridconverged = .false.
!          do while(.not.gridconverged) 
!             call splitGridFractal(grid%octreeRoot, real(100.*mHydrogen), 0.1, grid, gridconverged)
!          enddo
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

       case("magstream")
          nSource = 1
          allocate(source(1:1))
          source(1)%luminosity = grid%lCore
          source(1)%radius = ttaurirStar/1.d10
          source(1)%teff = 4000.
          source(1)%position = VECTOR(0.,0.,0.)
          call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(100.), dble(2.e8), 200)
          call normalizedSpectrum(source(1)%spectrum)
          call buildSphere(grid%starPos1, dble(grid%rCore), source(1)%surface, 400, contFluxFile)
          nu =1.d15
!           call createMagStreamSurface(source(1)%surface, grid, nu, coreContinuumFlux, fAccretion)
!           call testSurface(source(1)%surface)

          call readStreams(thisStream,nStreams,"stream.dat")

          call allocateStream(bigStream, nStreams*200)
          do i = 1, nStreams
             do j = 1, thisStream(i)%nsamples
                bigStream%nSamples = bigStream%nSamples + 1
                bigStream%rho(bigStream%nSamples) = thisStream(i)%rho(j)
                bigStream%temperature(bigStream%nSamples) = thisStream(i)%temperature(j)
                bigStream%position(bigStream%nSamples) = thisStream(i)%position(j)
                bigStream%velocity(bigStream%nSamples) = thisStream(i)%velocity(j)
             enddo
          enddo

          write(*,*) "stream read . now refining"
          write(*,*) "splitting on stream"
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, young_cluster, nDustType, &
               stream=bigStream)

          if (doTuning) call tune(6, "Magstream grid construction.") ! start a stopwatch


          write(*,*) "Bigstream has ",bigStream%nsamples, " samples" 
          call splitGridOnStream3(grid%octreeRoot,  grid, bigstream) 
          call freeStream(bigStream)

          nOctals = 0
          nVoxels = 0
          call countVoxels(grid%octreeRoot,nOctals,nVoxels)
          write(message,*) "Adaptive grid currently contains: ",nOctals," octals"
          call writeInfo(message, TRIVIAL)
          write(message,*) "                                : ",nVoxels," unique voxels"
          call writeInfo(message, TRIVIAL)
           


          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
          if (doTuning) call tune(6, "Magstream grid construction.") ! stop a stopwatch

       case DEFAULT
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, young_cluster, nDustType, romData=romData) 
          call writeInfo("First octal initialized.", TRIVIAL)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,romData=romData)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

           ! This section is getting rather long. Maybe this should be done in 
           ! wrapper subroutine in amr_mod.f90.
          if (geometry=="ttauri") then 
              ! Finding the total mass in the accretion flow
             mass_accretion_old = 0.0d0
             call TTauri_accretion_mass(grid%octreeRoot, grid, mass_accretion_old)
             write(message,*) "Total mass in accretion flow is ",  mass_accretion_old, "[g]"
             call writeInfo(message,FORINFO)

             if (ttau_fuzzy_edge) then
                call writeInfo("Fuzzy edge for TTauri accreation stream will be used.",FORINFO)
                call ttauri_fuzzy_edge(grid%octreeRoot)
                
                 ! We now have to correct the density of the accretion flow 
                 ! due to the fuzziness introduced above.
                mass_accretion_new = 0.0d0
                call TTauri_accretion_mass(grid%octreeRoot, grid, mass_accretion_new)
                write(message,*) "Total mass in accretion flow (after fuzzy edge) is ",  mass_accretion_new, "[g]"
                call writeInfo(message, FORINFO)
                mass_scale = mass_accretion_old/mass_accretion_new
                write(message,*) "Scaling the density by ", mass_scale , "[g]"
                call writeInfo(message, FORINFO)

                call TTauri_accretion_scale_density(grid%octreeRoot, grid, mass_scale) 
                 ! Check to see if it has scaled correctly.
                mass_accretion_new = 0.0d0
                call TTauri_accretion_mass(grid%octreeRoot, grid, mass_accretion_new)                 
                if (writeoutput) write(*,*) "Total mass in accretion flow after rescaling: ", mass_accretion_new, "[g]"
                if (writeoutput) write(*,*) " "  
             end if

             if (ttau_discwind_on) then
                if (writeoutput) write(*,*) " "
                if (writeoutput) write(*,*) "Adding TTauri disc wind to magnetosphere model... "
                call add_discwind(grid%octreeRoot, grid, ttauri_discwind, limitscalar2)
                if (writeoutput) write(*,*) " "
             elseif (ttau_jet_on) then
                if (writeoutput) write(*,*) " "
                if (writeoutput) write(*,*) "Adding TTauri jets to magnetosphere model... "
                call add_jet(grid%octreeRoot, grid, ttauri_jet)
                call finish_grid_jet(grid%octreeroot, ttauri_jet)
                if (writeoutput) write(*,*) " "
             end if


             if (grid%geometry=="ttauri" .and. ttau_turn_off_disc) then
                write(*,*) " "
                write(*,*) "Turning off the alpha disc read in from file, but keeping"
                write(*,*) "the magnetorsphere alive."
                if (.not.ttau_disc_on)   then
                   ! we need to create the disc parameter object 
                    call new(ttauri_disc, dble(TTauriDiskRin*TTauriRstar/1.0e10), 1.5d3*100.d0, &
                         dble(TTauriMstar/100.0), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, dble(TTauriMstar/mSol))
                 end if
                 call turn_off_disc(grid%octreeroot, grid, ttauri_disc)
              end if
              
              if (grid%geometry=="ttauri" .and. ttau_turn_off_jet) then
                 write(*,*) " "
                 write(*,*) "Turning off the jet read in from file, but keeping"
                 write(*,*) "the magnetorsphere alive."
                 if (.not.ttau_jet_on) then
                    ! the parameter file has to created here.
                    call new(ttauri_jet,  dble(TTauriRouter/1.0e10), dble(amrgridsize), &
                         JET_theta_j, dble(TTauriMstar/mSol), JET_Mdot, JET_a_param, JET_b_param,  &
                         JET_Vbase, JET_Vinf, JET_beta, JET_gamma, limitscalar, JET_T)
                 end if
                 call turn_off_jet(grid%octreeroot, grid, ttauri_jet)
              end if


           end if  ! gemoetry == "ttaruri"

           ! 
           if (doSmoothGrid) then
	      call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
              do
                 gridConverged = .true.
                 ! The following is Tim's replacement for soomthAMRgrid.
                 call myScaleSmooth(smoothFactor, grid, &
                      gridConverged,  inheritProps = .false., &
		      interpProps = .false.,  &
                      stellar_cluster=young_cluster, romData=romData)
                 if (gridConverged) exit
              end do
              call writeInfo("...grid smoothing complete", TRIVIAL)
           endif

           ! Smooth the grid with respect to optical depth, if requested
           if (doSmoothGridTau.and.mie) then
              call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
              call locate(grid%lamArray, nLambda,lambdaSmooth,ismoothlam)
              do j = iSmoothLam,  nLambda, 2
                 write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
                 call writeInfo(message, TRIVIAL)
                 do
                    gridConverged = .true.
                    call putTau(grid, grid%lamArray(j))
                    call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                         inheritProps = .false., interpProps = .false.)!, photosphereSplit = .not.variableDustSublimation)
                    if (gridConverged) exit
                 end do
              enddo
              call writeInfo("...grid smoothing complete", TRIVIAL)

              call writeVtkFile(grid, "aftersmooth.vtk")
              ! The tau smoothing may result in large differences in the size
              ! of neighbouring octals, so we smooth the grid again.
              
              if (doSmoothgrid) then
                 call writeInfo("Smoothing adaptive grid structure (again)...", TRIVIAL)
                 do
                    gridConverged = .true.
                 ! The following is Tim's replacement for soomthAMRgrid.
                    call myScaleSmooth(smoothFactor, grid, &
                         gridConverged,  inheritProps = .false., interpProps = .false., &
                         stellar_cluster=young_cluster, romData=romData)
                    if (gridConverged) exit
                 end do
                 call writeInfo("...grid smoothing complete", TRIVIAL)
              endif
           end if

     end select
   
!        call writeInfo("Unrefining optically thin cells...", TRIVIAL)
!        gridconverged = .false.
!        do while(.not.gridconverged)
!           gridconverged = .true.
!           call unrefineThinCells(grid%octreeRoot, grid, ismoothlam, gridconverged)
!        end do
!        call writeInfo("Done.", TRIVIAL)

!        call writeInfo("Unrefining optically thick cells...", TRIVIAL)
!        gridconverged = .false.
!        do while(.not.gridconverged)
!           gridconverged = .true.
!           call unrefineThickCells(grid%octreeRoot, grid, ismoothlam, gridconverged)
!        end do
!        call writeInfo("Done.", TRIVIAL)


        nOctals = 0
        nVoxels = 0
        call countVoxels(grid%octreeRoot,nOctals,nVoxels)
        write(message,*) "Adaptive grid contains: ",nOctals," octals"
        call writeInfo(message, TRIVIAL)
        write(message,*) "                      : ",nVoxels," unique voxels"
        call writeInfo(message, TRIVIAL)
        grid%nOctals = nOctals
        call howmanysplits()

        call writeInfo("Calling routines to finalize the grid variables...",TRIVIAL)
        call finishGrid(grid%octreeRoot, grid, romData=romData)
!        call writeInfo("...cleaning up dynamic memory",TRIVIAL)
!        call cleanupAMRgrid(grid%octreeRoot)

        call writeInfo("Constructing bitstrings...",TRIVIAL)
        call constructBitstrings(grid%octreeRoot)
        call writeInfo("Done.",TRIVIAL)


        call writeInfo("...final adaptive grid configuration complete",TRIVIAL)
        call howmanysplits()

#ifdef MPI
        if (grid%splitOverMPI) then
           call grid_info_mpi(grid, "info_grid.dat")
        else
           if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
        endif
#else
        if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
#endif

	if (lineEmission) then
           nu = cSpeed / (lamLine * angstromtocm)
           call contread(contFluxFile, nu, coreContinuumFlux)
           call buildSphere(grid%starPos1, dble(grid%rCore), starSurface, 400, contFluxFile)
           if (geometry == "ttauri") then
              call createTTauriSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion) 
           elseif (geometry == "magstream") then
              
           elseif (geometry == "romanova") then
              call createTTauriSurface2(starSurface, grid, romData, nu, coreContinuumFlux,fAccretion) 
           else
              call createSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion) 
           end if
        endif

        if ((geometry == "shakara").and.(nDustType>1)) then
           if ((nDustType ==2).and.(aMax(1) <  aMax(2))) then
              call writeInfo("Filling dust with large dust in midplane", FORINFO)
              call fillDustShakara(grid, grid%octreeRoot)
           else
              call writeInfo("Filling disc with uniform dust fractions", FORINFO)
              call fillDustUniform(grid, grid%octreeRoot)
           endif
        endif

        if ((geometry == "whitney").and.(nDustType ==4)) then
           call	fillDustWhitney(grid, grid%octreeRoot)	
        endif

        if (((geometry == "ppdisk").or.(geometry == "planetgap").or.(geometry=="warpeddisc")).and.(nDustType>1)) then
           call fillDustUniform(grid, grid%octreeRoot)
        endif
        
!        if (variableDustSublimation) then
!           call sublimateDust(grid, grid%octreeRoot, totFrac, nFrac)
!        endif



! ttauri source  / central star creation code moved from here to set_up_sources by th (18/7/08)
        if (sobolev) then
           !  calculate the statistical equilibrium (and hence the emissivities 
           !  and the opacities) for all of the subcells in an
           !  adaptive octal grid.
           !  Using a routine in stateq_mod module.
           if (writeoutput) write(*,*) "Calling statistical equilibrium routines..."
           if (doTuning) call tune(6, "amrStateq") ! start a stopwatch  
           if (grid%geometry=="ttauri" .or. grid%geometry(1:8)=="windtest" .or. &
                grid%geometry(1:9)=="luc_cir3d" .or. geometry == "romanova") then
              call amrStateq(grid, newContFluxFile,lte, nLower, nUpper, &
                             starSurface, recalcPrevious=.false.)

              call torus_mpi_barrier('waiting for other amrStatEq calls to return...')

              ! the 'lte' setting is intended so that a model with infall
              ! enhancement is first set up with LTE values. Only after we have 
              ! done the first time-dependent change so we calculate non-LTE 
              ! values.

              
              ! RK added the followings ------------------------------------------
              !  adding alpha disc to the grid
              if (ttau_disc_on) then
                 if (writeoutput) write(*,*) "Adding an accretion disc ... "
                 call add_alpha_disc(grid%octreeroot, grid, ttauri_disc)
                 ! finding the photon-dust interaction x-secions
                 call mieCrossSection(sigmaExt0, sigmaAbs0, sigmaSca0,  &
                      aMin(1), aMax(1), a0(1), qDist(1), pDist(1), grainType(1), &
                      ngrain, X_grain, grainname, lamLine)
                 if (writeoutput) write(*,*) " "
                 if (writeoutput) write(*,*) "Photon-dust cross section at lambda = ", lamLine, " [A]"
                 if (writeoutput) write(*,*) "    sigma(tot) = ", sigmaExt0, " [cm^2]"
                 if (writeoutput) write(*,*) "    sigma(abs) = ", sigmaAbs0, " [cm^2]"
                 if (writeoutput) write(*,*) "    sigma(sca) = ", sigmaSca0, " [cm^2]"
                 if (writeoutput) write(*,*) " "
                 call finish_grid(grid%octreeroot, grid, ttauri_disc, 1.0, &
                      sigmaAbs0, sigmaSca0, meanDustParticleMass)
                 !
                 !
!                 allocate(lambda(1:maxTau),tauExt(1:maxTau),tauAbs(1:maxTau),tauSca(1:maxTau),&
!                          contTau(1:maxTau,1:nLambda), linePhotonalbedo(1:maxtau))
!                 call integratePath(gridUsesAMR, VoigtProf, &
!                      lambdatau,  lamLine, VECTOR(1.0d-20,1.0d-20,1.0d-20), VECTOR(150.0,0,-200.), &
!                      VECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
!                      tauSca, linePhotonAlbedo, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .false. , &
!                      lamStart, lamEnd, nLambda, contTau, hitCore, thinLine, lineResAbs, .false., &
!                      .false., nUpper, nLower, 0., 0., 0., junk,&
!                      sampleFreq,intPathError, useInterp, grid%Rstar1, coolStarPosition)              
!                 write(*,*) " ----  Vertical optical depth at r = 0.1 AU is ", tauExt(ntau)
!                 write(*,*) " ----                           at lambda =  ", lambdatau
!                 write(*,*) "Rescaling scattering and absorption coefficients ... "
!                 call finish_grid(grid%octreeroot, grid, ttauri_disc, tauExt(ntau)/100.0)
!                 call integratePath(gridUsesAMR, VoigtProf, &
!                      lambdatau,  lamLine, VECTOR(1.,1.,1.), VECTOR(75.0,0,-50.), &
!                      VECTOR(0.,0.,1.), grid, lambda, tauExt, tauAbs, &
!                      tauSca, maxTau, nTau, thin_disc_on, opaqueCore, escProb, .false. , &
!                      lamStart, lamEnd, nLambda, contTau, hitCore, thinLine, lineResAbs, .false., &
!                      .false., nUpper, nLower, 0., 0., 0., junk,&
!                      sampleFreq,intPathError, useInterp, grid%Rstar1, coolStarPosition)              
!                 write(*,*) " ---- New Optical depth  at r = 0.05 AU    is ", tauExt(ntau)
!                 write(*,*) " ----                           at lambda =  ", lambdatau
!                 deallocate(lambda,tauExt,tauAbs,tauSca,contTau,linePhotonalbedo)
!                 write(*,*) ".. Finished rescaling scattering and absorption coefficients ... "
                 if (writeoutput) write(*,*) ".. Finished adding an acreation disc ... "
              end if
              !-------------------------------------------------------------------

           elseif (geometry == "cmfgen") then
              ! simply map CMFGEN opacity data to the AMR grid
              call map_cmfgen_opacities(grid)

           else

              call amrStateq(grid, newContFluxFile, lte, nLower, nUpper, &
                      starSurface, recalcPrevious=.false., ion_name=ion_name, ion_frac=ion_frac)
!              if (ttau_disc_on) then
!                 ! amrStateq will have messed up the disc, so we reset those cells
!                 call finish_grid(grid%octreeroot, grid, ttauri_disc, 1.0, sigmaAbs0, sigmaSca0)
!              end if

              call torus_mpi_barrier('waiting for other amrStatEq calls to return...')

           end if ! (grid%geometry=="ttaur" .or. ...)

           if (doTuning) call tune(6, "amrStateq") ! stop a stopwatch  
           if (writeoutput) write(*,*) "... statistical equilibrium routines complete"



        end if ! (lineEmission)

        !
        ! cleaning up unused memory here ....
        if ((geometry(1:7) == "cluster").or.(geometry(1:5)=="wr104") .or. geometry == "molcluster") then
           ! using the routine in sph_data_class
           call kill()
           ! using the routine in amr_mod.f90

          if (myRankIsZero) call delete_particle_lists(grid%octreeRoot)
          dummy = clusterparameter(VECTOR(0.d0,0.d0,0.d0),grid%octreeroot, subcell = 1, isdone = .true.)

        elseif (geometry == "luc_cir3d") then
           call deallocate_zeus_data()

! *** This part is commented out for now since romData is needed 
! *** for creating the surface elements later..... 
! *** This can be changed later (RK)
!
!        elseif (geometry == "romanova") then
!           call kill(romData)
!
        end if


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


!  if (geometry(1:6) == "ttauri" .and. myRankIsZero) then
!     call writeHartmannValues(grid,'hartmann_logNH')
!     call writeHartmannValues(grid,'hartmann_logNe')
!     call writeHartmannValues(grid,'hartmann_temperature')
!     call writeHartmannValues(grid,'hartmann_velPol')
!     call writeHartmannValues(grid,'hartmann_velAz')
!     call writeHartmannValues(grid,'hartmann_line')
!     !call writeHartmannValues(grid,'hartmann_Nelectron')
!     call writeHartmannValues(grid,'hartmann_Nlevel2')
!     call writeHartmannValues(grid,'hartmann_NH')
!     !call writeHartmannValues(grid,'hartmann_departCoeff')
!     call writeHartmannValues(grid,'hartmann_N')
!  end if
 
end subroutine amr_grid_setup


!-----------------------------------------------------------------------------------------------------------------------
subroutine set_up_sources

  use amr_mod, only: genericAccretionSurface
  use spectrum_mod, only: fillSpectrumBB, readSpectrum, normalizedSpectrum
  use starburst_mod, only: createsources
  use source_mod, only: source_within_octal, randomSource
  use cluster_class, only: get_nstar, n_stars_in_octal, get_a_star
  use surface_mod, only: createTTauriSurface, createTTauriSurface2, sumSurface, createSurface, testSurface

  integer :: isize
  integer, allocatable :: iseed(:)

  type(SOURCETYPE) :: a_star
  integer      :: nstar
  real(double) :: d1, d2, massRatio
  real         :: tmp
!  real :: sTot

  ! The source spectrum is normally a black body
!  if (.not.grid%lineEmission) then
!     stot = 0.
!     do i = 1, nLambda
!        sourceSpectrum(i) = bLambda(dble(xArray(i)), dble(teff))
!        stot = stot + sourceSpectrum(i)
!     enddo
!     sourceSpectrum  = sourceSpectrum / stot
!  endif


  ! set up the sources
  call writeInfo("Setting up sources", TRIVIAL)
  nSource = 0

  select case(geometry)

     case("ttauri", "luc_cir3d", "cmfgen", "romanova")

        if (.not.cmf) then
           nu = cSpeed / (lamLine * angstromtocm)
           call contread(contFluxFile, nu, coreContinuumFlux)
           call buildSphere(grid%starPos1, dble(grid%rCore), starSurface, 400, contFluxFile)
           if (geometry == "ttauri") then
              call createTTauriSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion) 
           elseif (geometry == "magstream") then
              
           elseif (geometry == "romanova") then
              call createTTauriSurface2(starSurface, grid, romData, nu, coreContinuumFlux,fAccretion) 
           else
              call createSurface(starSurface, grid, nu, coreContinuumFlux,fAccretion) 
           end if
        else
           nu = cSpeed / (lamLine * angstromtocm)

           nSource = 1
           allocate(source(1:1))
           source(:)%outsideGrid = .false.
           source(1)%luminosity = grid%lCore
           source(1)%radius = ttaurirStar/1.d10
           source(1)%teff = 4000.
           source(1)%position = VECTOR(0.,0.,0.)
           call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(100.), dble(2.e8), 200)
           call normalizedSpectrum(source(1)%spectrum)
           call buildSphere(grid%starPos1, dble(grid%rCore), source(1)%surface, 400, contFluxFile)
           call createTTauriSurface(source(1)%surface, grid, nu, coreContinuumFlux,fAccretion) 


        endif

        
    case("testamr","benchmark")
       nSource = 1
       allocate(source(1:1))
       source(:)%outsideGrid = .false.
       source(1)%luminosity = grid%lCore
       source(1)%radius = grid%rCore
       source(1)%teff = teff
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(lamstart), dble(lamEnd), nLambda)
       call normalizedSpectrum(source(1)%spectrum)
       call buildSphere(source(1)%position, source(1)%radius, source(1)%surface, 400, "blackbody")
       call sumSurface(source(1)%surface)
       call testSurface(source(1)%surface)

    case("magstream")
       nSource = 1
       if (.not.allocated(source)) then
          allocate(source(1:1))
          source(:)%outsideGrid = .false.
          source(1)%luminosity = grid%lCore
          source(1)%radius = ttaurirStar/1.d10
          source(1)%teff = 4000.
          source(1)%position = VECTOR(0.,0.,0.)
          call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(100.), dble(2.e8), 200)
          call normalizedSpectrum(source(1)%spectrum)
          call buildSphere(grid%starPos1, dble(grid%rCore), source(1)%surface, 400, contFluxFile)
       endif
       nu =1.d15
       call genericAccretionSurface(source(1)%surface, grid, nu, coreContinuumFlux, fAccretion)
       call testSurface(source(1)%surface)

    case("melvin")
       nSource = 1
       teff = 30000.
       allocate(source(1:1))
       source(:)%outsideGrid = .false.
       source(1)%luminosity = fourPi * (10.*rsol)*(10.*rsol) * stefanBoltz * teff**4
       source(1)%radius = 10.*rSol / 1.e10
       source(1)%teff = teff
       source(1)%position = VECTOR(0.,0.,0.)
       fac = 1.e8*cspeed/5.d16
       call fillSpectrumBB(source(1)%spectrum, dble(source(1)%teff), fac, 1000.d4,1000)
       call normalizedSpectrum(source(1)%spectrum)
       rstar = source(1)%radius

    case("whitney")
       nSource = 1
       allocate(source(1:1))
       source(:)%outsideGrid = .false.
       source(1)%luminosity = fourPi * rStellar**2 * stefanBoltz * teff**4
       source(1)%radius = rStellar / 1.e10
       source(1)%teff = teff
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(teff),  dble(lamStart), dble(lamEnd),nLambda)
       call normalizedSpectrum(source(1)%spectrum)
       rstar = source(1)%radius

    case("toruslogo")
       nSource = 1
       allocate(source(1:1))
       source(:)%outsideGrid = .false.
       source(1)%luminosity = lsol
       source(1)%radius = rSol / 1.e10
       source(1)%teff = 6000.
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(6000),  dble(lamStart), dble(lamEnd),nLambda)
       call normalizedSpectrum(source(1)%spectrum)
       rstar = source(1)%radius


    case("starburst")
       call random_seed(size=iSize)
       allocate(iSeed(1:iSize))
       iseed = 2343245
       call random_seed(put=iSeed)
       deallocate(iSeed)

       allocate(source(1:10000))
       source(:)%outsideGrid = .false.
       nSource = 0
       call createSources(nSource, source, "instantaneous", 1.d6, 1.d3, 1.d0)
       call writeInfo(message, TRIVIAL)
       call init_random_seed()

    case("wr104")
       nSource = 2

       allocate(source(1:nSource)) 
       source(:)%outsideGrid = .false.
       source(1)%teff = 30000.  ! o star
       source(1)%radius = 10. * rSol / 1.e10
       source(1)%position = (VECTOR(1.,0.,0.)*dble(autocm))/1.d10
       source(1)%luminosity = fourPi * stefanBoltz * (10.*rSol)**2.0 * (source(1)%teff)**4
       call readSpectrum(source(1)%spectrum, "ostar.flx", ok)
       call normalizedSpectrum(source(1)%spectrum)

       source(2)%teff = 40000.             ! wr star 
       source(2)%radius = 20. * rSol / 1.e10
       source(2)%position = (VECTOR(-1.,0.,0.)*dble(autocm))/1.d10
       source(2)%luminosity = 0.5 * source(1)%luminosity
       call readSpectrum(source(2)%spectrum, "wr.flx", ok)
       call normalizedSpectrum(source(2)%spectrum)

    case("lexington","fractal")
       nSource = 1
       allocate(source(1:nSource)) 
       source(:)%outsideGrid = .false.
       source(1)%teff = 40000.  ! o star
       source(1)%radius = 18.67 * rSol / 1.e10
       source(1)%position = VECTOR(0.,0.,0.)
       source(1)%luminosity = fourPi * stefanBoltz * (source(1)%radius*1.e10)**2.0 * (source(1)%teff)**4
       if (writeoutput) write(*,*) "Lexington source: ",source(1)%luminosity/1.e37
       fac = 1.e8*cspeed/5.d16
       call fillSpectrumBB(source(1)%spectrum, dble(source(1)%teff), fac, 1000.d4,1000)
       call normalizedSpectrum(source(1)%spectrum)

    case("bonnor")
       nSource = 1
       allocate(source(1:nSource)) 
       source(1)%teff = 40000.  ! o star
       source(1)%radius = 18.67 * rSol / 1.e10
       source(1)%position = VECTOR(0.,0.,0.)
       source(1)%luminosity = fourPi * stefanBoltz * (source(1)%radius*1.e10)**2.0 * (source(1)%teff)**4

       source(1)%luminosity = source(1)%luminosity * (2.d0*1.5d9*1.d10)**2 / &
            (fourPi*(50.d0* pctocm)**2)
       if (writeoutput) write(*,*) "Lexington source: ",source(1)%luminosity/1.e37
       fac = 1.e8*cspeed/5.d16
       call fillSpectrumBB(source(1)%spectrum, dble(source(1)%teff), fac, 1000.d4,1000)
       call normalizedSpectrum(source(1)%spectrum)
       source%outsideGrid = .true.

    case("symbiotic")
       nSource = 2
       allocate(source(1:nSource)) 
       source(:)%outsideGrid = .false.
       source(1)%teff = 2500.  
       source(1)%radius = 100. * rSol / 1.e10
       source(1)%position = VECTOR(-250.*rSol/1.e10,0.,0.)
       source(1)%luminosity = fourPi * stefanBoltz * (source(1)%radius*1.e10)**2.0 * (source(1)%teff)**4
       call fillSpectrumBB(source(1)%spectrum, dble(source(1)%teff), 10.d0, 1000.d4,1000)
       call normalizedSpectrum(source(1)%spectrum)

       source(2)%teff = 150000.  
       source(2)%radius = 5000.e5 / 1.e10 ! 5000 km
       source(2)%position = VECTOR(250.*rSol/1.e10,0.,0.)
       source(2)%luminosity = fourPi * stefanBoltz * (source(2)%radius*1.e10)**2.0 * (source(2)%teff)**4
       call fillSpectrumBB(source(2)%spectrum, dble(source(2)%teff), 10.d0, 1000.d4,1000)
       call normalizedSpectrum(source(2)%spectrum)
       write(*,*) "!!!!!!! light ratio",source(1)%luminosity/source(2)%luminosity



    case ("cluster", "molcluster")
       ! Extract some info from cluster object.
       nstar = get_nstar(young_cluster)  ! number of stars in the cluster
       nSource =  n_stars_in_octal(young_cluster, grid%octreeRoot)
       
       ! copy the star over in the array.
       ! This is ugly. Maybe lucyRadiativeEquilibriumAMR should be changed to take
       ! an cluster_class object as an input variable in future.
       ALLOCATE(source(nSource))
       source(:)%outsideGrid = .false.
       
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

    case("shakara","clumpydisc","wrshell","warpeddisc","iras04158")
       nSource = 1
       allocate(source(1:1))
       source(:)%outsideGrid = .false.
       source(1)%radius = grid%rCore
       source(1)%teff = teff   
       source(1)%position = VECTOR(0.,0.,0.)
       if (contFluxfile .eq. "blackbody") then
          call fillSpectrumBB(source(1)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd),nLambda, lamArray=xArray)
          call normalizedSpectrum(source(1)%spectrum)
          tmp = source(1)%radius * 1.e10  ! [cm]
          source(1)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(1)%teff)**4
       else
          call buildSphere(source(1)%position, source(1)%radius, source(1)%surface, 400, contFluxFile)
          call readSpectrum(source(1)%spectrum, contfluxfile, ok)
          call normalizedSpectrum(source(1)%spectrum)
          call sumSurface(source(1)%surface, source(1)%luminosity)
       endif

       tmp = source(1)%radius * 1.e10  ! [cm]
       fac = fourPi * stefanBoltz * (tmp*tmp) * (source(1)%teff)**4
       if (abs(fac-source(1)%luminosity)/source(1)%luminosity > 0.01d0) then
          if (writeoutput) then
             write(*,*) "WARNING: luminosity from effective temperature and that from the SED differ by >1%"
             write(*,*) "Lum from Teff (lSol): ",fac/lSol
             write(*,*) "Lum from SED  (lSol): ",source(1)%luminosity/lSol
             write(*,*) "Implied source of accretion luminosity of: ",(source(1)%luminosity-fac)/lSol
             fac = (source(1)%luminosity-fac) * tmp / (bigG * mCore)
             write(*,*) "Mass accretion rate could be ",fac/mSol/ secstoYears, " solar masses/year"
          endif
       endif


    case("circumbin")
       nSource = 2
       allocate(source(1:2))
       source(:)%outsideGrid = .false.
       massRatio = mstar2/mstar1
       source(1)%radius = rStar1
       source(1)%teff = teff1  
       source(1)%position = VECTOR(dble(1.-eccentricity)*binarySep/(1.+massRatio), 0., 0.)
       tmp = rstar1 * 1.e10  ! [cm]
       source(1)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(1)%teff)**4
       if (contFluxfile1 .eq. "blackbody") then
          call fillSpectrumBB(source(1)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd),nLambda, lamArray=xArray)
       else
          call buildSphere(source(1)%position,source(1)%radius, source(1)%surface, 400, contFluxFile1)
          call readSpectrum(source(1)%spectrum, contfluxfile1, ok)
       endif
       call normalizedSpectrum(source(1)%spectrum)

       source(2)%radius = rStar2
       source(2)%teff = teff2  
       source(2)%position = VECTOR(-dble(1.-eccentricity)*binarySep*(1.-1./(1.+massratio)), 0., 0.)
       tmp = rstar2 * 1.e10  ! [cm]
       source(2)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(2)%teff)**4
       if (contFluxfile2 .eq. "blackbody") then
          call fillSpectrumBB(source(2)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd),nLambda, lamArray=xArray)
       else
          call buildSphere(source(2)%position, source(2)%radius, source(2)%surface, 400, contFluxFile2)
          call readSpectrum(source(2)%spectrum, contfluxfile2, ok)
       endif
       call normalizedSpectrum(source(2)%spectrum)

    case("gammavel")
       nSource = 2
       allocate(source(1:2))
       source(:)%outsideGrid = .false.


       massRatio = mass1/mass2

       d1 = binarySep * (1./(massRatio+1.))
       d2 = binarySep - d1

! source 1 is the WR star
       source(1)%radius = rstar1
       source(1)%teff = teff1
       source(1)%position = VECTOR(0.,0.,-d1)
       tmp = source(1)%radius * 1.e10  ! [cm]
       source(1)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(1)%teff)**4

       call buildSphere(source(1)%position,source(1)%radius, source(1)%surface, 400, contFluxFile1)


       if (contFluxfile1 .eq. "blackbody") then
          call fillSpectrumBB(source(1)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd), nLambda, lamArray=xArray)
       else
          call readSpectrum(source(1)%spectrum, contfluxfile1, ok)
       endif
       call normalizedSpectrum(source(1)%spectrum)

! source 2 is the WR star
       source(2)%radius = rstar2
       source(2)%teff = teff2
       source(2)%position = VECTOR(0.,0.,d2)
       tmp = source(2)%radius * 1.e10  ! [cm]
       source(2)%luminosity = fourPi * stefanBoltz * (tmp*tmp) * (source(2)%teff)**4

       call buildSphere(source(2)%position, source(2)%radius, source(2)%surface, 400, contFluxFile2)

       if (contFluxfile2 .eq. "blackbody") then
          call fillSpectrumBB(source(2)%spectrum, dble(teff), &
               dble(lamStart), dble(lamEnd), nLambda)
       else
          call readSpectrum(source(2)%spectrum, contfluxfile2, ok)
       endif
       call normalizedSpectrum(source(2)%spectrum)

    ! chris (26/05/04)
    case ("ppdisk")
       nSource = 1
       allocate(source(1:nSource))
       source(:)%outsideGrid = .false.
!       source(1)%teff = 5780.
       source(1)%teff = Teff
!       source(1)%luminosity = 3.83d+32 
!       source(1)%luminosity = 10.**((4.75-Mbol)/2.5) * lSol
! Baraffe, et al. use 4.64 for solar bolometric magnitude
       source(1)%luminosity = 10.**((4.64-Mbol)/2.5) * lSol
!       source(1)%radius = 6.96d+0
       source(1)%radius = sqrt(source(1)%luminosity/(4.*pi*stefanBoltz*Teff**4))/1.d10
       source(1)%position = VECTOR(0.,0.,0.)

       call fillSpectrumBB(source(1)%spectrum, dble(teff), dble(lamStart), dble(lamEnd), nLambda)
       call normalizedSpectrum(source(1)%spectrum)

    case ("planetgap")
       nSource = 1
       allocate(source(1:nSource))
       source(:)%outsideGrid = .false.
       source(1)%teff = Teff
       source(1)%luminosity = fourPi*(rcore*1.e10)**2 * stefanBoltz*teff**4
       source(1)%radius = rCore
       source(1)%position = VECTOR(0.,0.,0.)
       call fillSpectrumBB(source(1)%spectrum, dble(teff), dble(lamStart), dble(lamEnd), nLambda)
       call normalizedSpectrum(source(1)%spectrum)


    case default
       ! Allocating the source array with size =0 to avoid, non-allocated array passed problem
       ! in subroutine initPhoton...
       if (.not. allocated(source)) ALLOCATE(source(0))
       nSource = 0  ! This must be zero!
       
    end select

  if (nSource > 0) then
     call randomSource(source, nSource, i, xArray, nLambda, initialize=.true.)  
     call writeInfo("Sources set up.",TRIVIAL)
  endif


end subroutine set_up_sources

!-----------------------------------------------------------------------------------------------------------------------

subroutine post_initAMRgrid 

  use amr_mod, only: find_average_temperature, findTotalMass, scaleDensityAMR

  real         :: scaleFac
  real(double) :: totalMass, totalmasstrap, maxrho, minrho
  real(double) :: T_ave   ! average temperature of cluster
  real(double) :: T_mass  ! mass weighted temperature

  if (geometry == "wr104") then

     totalMass = 0.d0
     call findTotalMass(grid%octreeRoot, totalMass)
     scaleFac = massEnvelope / totalMass
     if (writeoutput) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
     call scaleDensityAMR(grid%octreeRoot, scaleFac)

  else if (geometry == "cluster") then
!        call find_average_temperature(grid, T_ave, T_mass, totalMass)
!        if (writeoutput) write(*,*) " "
!        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Total dust mass in cluster  : ", TotalMass, " [g]"
!        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Ave. temperature of cluster : ", T_ave,    " [K]"
!        if (writeoutput) write(*,'(a,1pe12.4, a5)') "Mass weighted ave. temperature of cluster : ",T_mass,  " [K]"
!        if (writeoutput) write(*,*) " "
!        ! computing the max, min and average tau at 2 micron (of all
!        ! active cells)
!        ilambda = findIlambda(20000.0, xArray, nLambda, ok)
!        call find_max_min_ave_tau(grid, ilambda, tau_max, tau_min,  tau_ave)
!        if (writeoutput) write(*,*) " "
!        if (writeoutput) write(*,'(a,1pe12.4, a11)') "Max tau(1 micron) = ", tau_max, " [contiuum]"
!        if (writeoutput) write(*,'(a,1pe12.4, a11)') "Min tau(1 micron) = ", tau_min, " [contiuum]"
!        if (writeoutput) write(*,'(a,1pe12.4, a11)') "Ave tau(1 micron) = ", tau_ave, " [contiuum]"
!        if (writeoutput) write(*,*) " "

  else if (geometry == "ttauri") then
     call find_average_temperature(grid, T_ave, T_mass, totalMass)
     if (writeoutput) write(*,*) " "
     if (writeoutput) write(*,'(a,1pe12.4, a5)') "Total mass in computational domain : ", TotalMass, " [g]"
     if (writeoutput) write(*,'(a,1pe12.4, a5)') "Ave. temperature in computational domain : ", T_ave,    " [K]"
     if (writeoutput) write(*,'(a,1pe12.4, a5)') "Mass weighted ave. temperature in computational domain : ",T_mass,  " [K]"
     if (writeoutput) write(*,*) " "
  end if


  if (associated(grid%octreeRoot)) then
     totalMass =0.d0
     if (.not.grid%splitOverMpi) then
        if(geometry .eq. 'molcluster') then 
           call findTotalMass(grid%octreeRoot, totalMass, totalmasstrap = totalmasstrap, minrho = minrho, maxrho = maxrho)
        else
           call findTotalMass(grid%octreeRoot, totalMass, minrho = minrho, maxrho = maxrho)
        endif
     else
#ifdef MPI
        if (myrankglobal /= 0) then
           call findMassOverAllThreads(grid, totalMass)
        endif
#endif
     endif
     write(message,*) "Mass of envelope: ",totalMass/mSol, " solar masses"
     call writeInfo(message, TRIVIAL)
     if(geometry .eq. 'molcluster') then
        write(message,*) "Mass of envelope (TRAP): ",totalMasstrap/mSol, " solar masses"
        call writeInfo(message, TRIVIAL)
        write(message,*) "Maximum Density: ",maxrho, " g/cm^3"
        call writeInfo(message, TRIVIAL)
        write(message,*) "Minimum Density: ",minrho, " g/cm^3"
        call writeInfo(message, TRIVIAL)
     endif
  endif
        

  if (mie) then
     if (geometry == "shakara" .or. geometry .eq. 'iras04158' .or. geometry == "circumbin") then

!        sigma0 = totalMass / (twoPi*(rOuter*1.e10-rInner*1.e10)*1.*real(autocm)) ! defined at 1AU

!           sigma0 = totalMass * real(betaDisc-alphaDisc+2) / & 
!                (twoPi*((rOuter*1.e10)**(betaDisc-alphaDisc+2) - (rInner*1.e10)**(betaDisc-alphaDisc+2))*1.*auToCm)

!           sigma0 = rho0 * (rinner*1.e10/(1.*autocm))**alphaDisc * &
!           (height*1.e10) * (autocm / (100.d0*autocm))**betaDisc * sqrt(twopi)



        sigma0 = rho0 * (height*1.e10) * ((rinner*1.e10) / (100.d0*autocm))**betaDisc * sqrt(twopi)

        if (writeoutput) write(*,*) "Sigma0: ",sigma0
     endif

     if (geometry == "ppdisk") then
        sigma0 = totalMass / ((amrgridsize**2 - twoPi*0.4)*1.e10*1.*autocm) ! defined at 1AU
        if (writeoutput) write(*,*) "Sigma0: ",sigma0
     endif

     if (geometry == "planetgap") then
        rho0  = mDisc *(betaDisc-alphaDisc+2.) / ( twoPi**1.5 * height * (rCore*1.e10) &
             * (rCore*1.e10)**(alphaDisc-betaDisc) * &
             (((rOuter*1.e10)**(betaDisc-alphaDisc+2.)-(rInner*1.e10)**(betaDisc-alphaDisc+2.))) )
        h = height * rCore * 1.e10
        sigma0 = rho0 * sqrt(twoPi) * h
        
        if (writeoutput) write(*,*) "Sigma0: ",sigma0
     endif

     if (geometry == "warpeddisc") then
        rho0  = mDisc *(betadisc-alphadisc+2.) / ( twoPi**1.5 * (height*1.e10)  &
             * (rOuter*1.d10)**(alphadisc-betadisc) * ( &
             ((router*1.d10)**(betadisc-alphadisc+2.)-(rInner*1.d10)**(betadisc-alphadisc+2.))) )
        sigma0 = rho0 * sqrt(twoPi) * height * 1.e10
     endif

  endif

  call torus_mpi_barrier

end subroutine post_initAMRgrid

!-----------------------------------------------------------------------------------------------------------------------

subroutine do_lucyRadiativeEq

  use benchmark_mod, only: check_benchmark_values
  use lucy_mod, only: lucyRadiativeEquilibrium, lucyRadiativeEquilibriumAMR, allocateMemoryForLucy
  use disc_hydro_mod, only: verticalHydrostatic
  use cluster_utils, only: analyze_cluster
  use cluster_class, only: reassign_10k_temperature, restrict, kill_all

  type(VECTOR) :: outVec
  outVec = (-1.d0)* originalViewVec

     if (doTuning) call tune(6, "LUCY Radiative Equilbrium")  ! start a stopwatch
 
     call allocateMemoryForLucy(grid%octreeRoot)

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
                   lucy_undersampled, IterLucy )
           endif

        endif

        

        if (myRankIsZero .and. writeLucy) call writeAMRgrid(lucyFilenameOut,writeFileFormatted,grid)
     endif

     if (doTuning) call tune(6, "LUCY Radiative Equilbrium")  ! stop a stopwatch

     if (grid%geometry(1:7) == "cluster") then

        if (myRankIsZero) then
           ! Finding apperant magnitudes and colors of the stars in cluster
           write (*,*) " "
           write (*,*) "Computing the magnitudes and colors of stars ..."
           ! -- note grid distance here is in [pc]
           call analyze_cluster(young_cluster,outVec,dble(gridDistance),grid)
        end if

        call torus_mpi_barrier

        ! restricting the sed and images calculations to the star 
        ! with ID number = idx_restrict_star. 
        ! By default idx_restrict_star =0 which means all stars the cluster
        ! are included in the calculations.  Specifty a value in your
        ! parameter file to restrict to one star.  (10 AU cylinder
        ! zone is used to restrict the effective computational domain.
        call restrict(grid%octreeroot, idx_restrict_star, nsource, &
                      source,  outVec, 7.5d4) ! the last value is 50 AU (in 10^10cm)
	!                      source,  s2o(outVec), 1.5d4) ! the last value is 10 AU (in 10^10cm)

        call reassign_10K_temperature(grid%octreeroot)
        ! delete the cluster object since it won't be used any more.
        call kill_all(young_cluster)
     endif


     call torus_mpi_barrier

! Check benchmark results against values read in from file
! If the file does not exist then the routine will return cleanly
     if (myRankIsZero .and. grid%geometry == 'benchmark') call check_benchmark_values(grid, "part_out.dat")

end subroutine do_lucyRadiativeEq

!-----------------------------------------------------------------------------------------------------------------------

subroutine set_emission_bias

  use cluster_class, only: assign_emission_bias
  use amr_mod, only: set_bias_shakara, set_bias_cmfgen, set_bias_rosseland, &
       set_bias_ttauri, set_bias_whitney, setbiasamr

  type(VECTOR) :: outVec
  outVec = (-1.d0)* originalViewVec

  !
  ! Setting the emission bias.
  !
     select case(grid%geometry)
        case("testamr")
           call setBiasAMR(grid%octreeRoot, grid)
        case("cluster")
           ! Computing the emission bias based on the optical depth at
           ! 2 microns. 
           call assign_emission_bias(grid%octreeroot, grid, 20000.0, xArray, nLambda)
        case("ttauri")
           call set_bias_ttauri(grid%octreeRoot, grid, lamline, outVec)
        case("cmfgen")
           call set_bias_cmfgen(grid%octreeRoot, grid, lamline)
        case("shakara","ppdisk","planetgap")
           if (forcedWavelength) then
              call locate(grid%lamArray, nLambda, usePhotonWavelength, itestlam)
              call set_bias_shakara(grid%octreeRoot, grid, ilam=itestlam, ross=.false.)
           else
              call locate(grid%lamArray, nLambda, usePhotonWavelength, itestlam)
              call set_bias_shakara(grid%octreeRoot, grid, ilam=1, ross=.true.)
!              call setBiasOnTau(grid)
           endif
!        case("ppdisk")
!           call locate(grid%lamArray, nLambda, usePhotonWavelength, itestlam)
!           call setDiscPhotosphereBias(grid, itestlam)
        case("whitney")
           call set_bias_whitney(grid%octreeRoot, grid)
        case("warpeddisc")
           call set_bias_rosseland(grid%octreeRoot, grid)
!	   call set_bias_radius(grid%octreeRoot, 3)
        case DEFAULT
           continue
     end select

end subroutine set_emission_bias

!-----------------------------------------------------------------------------------------------------------------------

subroutine initialize_blobs
  
  use blob_mod, only: addnewblobs, moveblobs, writeBlobs

  character(len=80) :: filename, specFile
  real, parameter   :: blobTime = 1000.

  allocate(blobs(1:maxBlobs))
  if (freshBlobs) then
     do i = 1 , maxBlobs
        blobs(i)%inUse = .false.
     enddo
     if (writeoutput) write(*,'(a)') "Running blobs for five days..."
   
     ! now we run a few of days worth of blobs
  
     t1 = 0.
     t2 = 240. * 60. * 60.
     dTime = (t2 - t1) / 100.
     do i = 1, 100
        call addNewBlobs(grid, maxBlobs, blobs, blobTime, dTime, &
                              nCurrent, blobContrast)
        call moveBlobs(maxBlobs, blobs, 0., dTime, grid)
        if (writeoutput) write(*,*) i,nCurrent
     enddo
  
     j = 0 
     do i = 1, maxBlobs
        if (blobs(i)%inUse) j = j + 1
     enddo
  
     if (writeoutput) write(*,'(a,i3)') "Number of blobs after 5 days: ",j
  
  
     dTime = 0.
     if (nPhase /= 1) then
        dTime = (timeEnd - timeStart) / real(nPhase-1)
     endif
  
     open(50,file="files.lis",status="unknown",form="formatted")
  
     if (writeoutput) write(*,'(a)') "Running blobs and writing configuration files"
     do i = 1, nPhase
  
        grid%timeNow = timeStart + (timeEnd-timeStart)*real(i-1)/real(nPhase-1)	    
        write(specFile,'(a,i3.3,a)') trim(outfile),i,".dat"
        write(50,*) specFile(1:30), grid%timeNow
  
        call addNewBlobs(grid, maxBlobs, blobs, blobTime, dTime, nCurrent, blobContrast)
        call moveBlobs(maxBlobs, blobs, 0., dTime, grid)
        write(filename,"(a,i3.3,a)") "run",i,".blob"
        call writeBlobs(filename, maxBlobs, blobs)
  
     enddo
  
     close(50)
  
  endif

end subroutine initialize_blobs

!-----------------------------------------------------------------------------------------------------------------------

end program torus

