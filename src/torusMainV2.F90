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
!
program torus
  use torus_version_mod
  use inputs_mod         ! variables filled by inputs subroutine
  use citations_mod
  use dimensionality_mod
!  use intensity_storage_mod
  use constants_mod
  use messages_mod
  use mpi_global_mod
#ifdef USEZLIB
  use zlib_mod
#endif
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
  use turbulence_mod
  use loadbalance_mod
  use biophysics_mod
  use photoionAMR_mod

#ifdef CHEMISTRY
  use chemistry_mod
#endif

!  use pah_mod
!  use zlib_mod
  use sph_data_class, only: deallocate_sph
#ifdef MPI
  use mpi
  use mpi_fitter_mod
#endif



  implicit none

  character(len=80) :: message, dataDirectory
!  character(len=10) :: stringArray(10)
!  integer :: i
!  type(PAHtabletype) :: PAHtable
  type(GRIDTYPE) :: grid
!  type(VECTOR) :: box(64,64,64)
#ifdef MPI
  ! For MPI implementations =====================================================
  integer ::   ierr           ! error flag
#endif
  !#ifdef _OPENMP
!  call omp_set_dynamic(.false.)
!#endif


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
  loadBalancingThreadGlobal = .false.

!  call createBox(box,64)

  
  if (TorusMpi) then 
     if (myRankGlobal/=1) writeoutput  = .false.
     if (myRankGlobal/=1) doTuning     = .false.
     if (myRankGlobal/=0) myRankIsZero = .false.
  end if

!  call readPAHemissivityTable(PAHtable, "MW3.1_00")

  call initializeCodeUnits()

  call writeTorusBanner()

  call unixGetEnv("TORUS_DATA", dataDirectory)
  if (len(trim(dataDirectory)) == 0) then
    call writeFatal("You need to set the environment variable TORUS_DATA")
    stop
  endif



  call initBibCode()
  call setVersion("V4.0")
  grid%version = torusVersion
  currentlyDoingHydroStep = .false.
  verbosityLevel = 5
  call writeBanner("TORUS "//trim(torusVersion)//" model","-",IMPORTANT)

  ! For time statistics
  if (doTuning) call tune(6, "Torus Main") ! start a stopwatch  

  call inputs()

  call setUpBibcodesOnParametersFile()
!  if (storeScattered) call initHealpix(1)

#ifdef MPI
  ! Set up amrCOMMUNICATOR and global mpi groups
  call setupAMRCOMMUNICATOR
  write(*,*) "nthreadsglobal ",nThreadsGlobal
!  call setupFitterCommunicators
!  do i = 0, nThreadsGlobal-1
!     if (myrankGlobal == i) then
!        write(*,*) "rank ",myrankGlobal, " fitter ",myFitterSetGlobal
!     endif
!  enddo
     
#endif

! Report build options 
  call report_parallel_type
  call write_job_info_file

  ! set up a random seed
  call run_random_test_suite

  smallestCellSize = amrGridSize / dble(2**maxDepthAMR)

  do iModel = nModelStart, nModelEnd, modelStep
     if (multimodels) then
        write(message,'(a,i6.6)') "Performing calculation for model number ", iModel
        call writeBanner(message,"-",TRIVIAL)
     endif

#ifdef USEZLIB

#else
     useBinaryXMLVTKfiles = .false.
#endif

     call setupMicrophysics(grid)

!     if (trim(geometry) == "envelope") then
!        call calculategrid()
!        goto 666
!     endif

     if (doSetupAMRgrid) then

        call  setupamrgrid(grid)



!        call sanityCheckGrid(grid)

     !  call checkAMRgrid(grid, .false.)

     !  call writeAMRgrid("test.dat",.false.,grid)
     !  call writeVtkFile(grid, "mpi.vtk",  valueTypeString=(/"mpithread"/))
     !  call torus_mpi_barrier
     !  goto 666
!        write(*,*) "OMP THREAD NUMBER ",omp_get_thread_num()

       if (.not.readGrid) then
          call writeVtkFile(grid, "rho.vtk")
!          call writeVtkFile(grid, "rho2.vtk",valuetypestring=(/"velcall"/))
       endif

       if (dustPhysics) then
          if (nDustType >= 1) then
!             do i = 1, nDustType
!                write(stringArray(i),'(a,i1.1)') "dust",i
!             enddo
!             call writeVTKfile(grid,"dust.vtk",valueTypeString=stringArray(1:nDustType))
             call writeVTKfile(grid, "dust.vtk",valueTypestring=(/"dust"/))
          endif
       endif

!       call writeVtkFile(grid, "mpi.vtk",valueTypeString=(/"rho        ","mpithread        "/))
    endif
     call setupGlobalSources(grid)

#ifdef MPI
!     call testBranchCopying(grid)
#endif

! SPH data may be required for setting up sources so only free memeory after this is done

#ifdef SPH
        call deallocate_sph
#endif

     call torus_mpi_barrier

     call writeBanner("Run-time messages","+",TRIVIAL)

     call randomNumberGenerator(randomSeed=.true.)

     call doPhysics(grid)

     call doOutputs(grid)

     if (doTuning) call tune(6, "Torus Main") ! stop a stopwatch  

     if (associated(grid%octreeRoot)) then
        call writeInfo("Before deallocation",TRIVIAL)
        call resetGlobalMemory(grid)
     endif
     call torus_mpi_barrier
     if (associated(grid%octreeRoot)) then
        call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
     else
        call writewarning("Want to delete grid but octreeroot not associated")
     endif


     call freeGrid(grid)
     call freeGlobalSourceArray()
  enddo

#ifdef MPI
  call torus_mpi_barrier
  if (hydrodynamics) call freeAMRCOMMUNICATOR
  call MPI_FINALIZE(ierr)
#endif

  if (Writeoutput) call writeBibtex("torus.bib")


  call writeBanner("Torus completed","o",TRIVIAL)

contains 

!  subroutine calculategrid()
!    implicit none
!    integer :: iMass, iRinner, iRadius, iLum, iAcc, iRouter, iMassStar
!    integer :: nMass, nRinner, nRadius, nLum, nAcc, nRouter, nMassStar
!
!    real(double) :: massStart, massEnd
!    real(double) :: massStarStart, massStarEnd
!    real(double) :: AccStart, accEnd
!    real(double) :: rInnerStart, rInnerEnd
!    real(double) :: lumStart, lumEnd
!    real(double) :: radiusStart, radiusEnd
!    real(double) :: rOuterStart, rOuterEnd
!    real(double) :: massStar, lum, lacc, tacc
!    
!    massStart = 1.d0; massEnd = 1.d0; nMass = 1
!    massStarStart  = 1.d0*mSol; massStarEnd = 1.d0*mSol; nMassStar = 1
!    radiusStart = 3.58d0 * rSol; radiusEnd = 3.58d0 * rSol; nRadius = 1
!    lumStart = 11.6d0*lSol; lumEnd = 11.6d0 * lSol; nLum = 1
!    accStart = 1.d-20; accEnd = 1.d-20; nAcc = 1
!    rOuterStart = 4.d4; rOuterEnd = 4.d4; nrouter = 1
!    rInnerStart = 400.d0; RinnerEnd = 400.d0; nRinner = 1
!    
!    
!    do iMass = 1, nMass
!       do iRinner = 1, nRinner
!          do iRadius = 1, nRadius
!             do iLum = 1, nLum
!                do iAcc = 1, nAcc
!                   do iRouter = 1, nRouter
!                      do iMassStar = 1, nMassStar
!
!                      massEnvelope =  real(getGridValue(iMass, massStart, massEnd, nMass, .false.)*mSol)
!                      rInner =  real(getGridValue(iRinner, rinnerStart, rinnerEnd, nrInner, .true.)*autocm/1.d10)
!                      rOuter = real(getGridValue(irOuter, rOuterStart, rOuterEnd, nrOuter, .false.)*autocm/1.d10)
!!
!                      massStar =  getGridValue(iMassStar, massStarStart, massStarEnd, nMassStar, .false.)
!                      lum = getGridValue(iLum, lumStart, lumEnd, nlum, .true.)
!                      radius = real(getGridValue(iRadius, radiusStart, radiusEnd, nRadius, .false.))
!                      mDot = real(getGridValue(iAcc, accStart, accEnd, nAcc, .true.)*msol / (365.25d0*24.d0*3600.d0))
!                      teff = real((lum / (fourPi * radius**2 * stefanBoltz))**0.25d0)
!
!                      lAcc = bigG * mDot * massStar / radius
!                      tAcc = (lAcc / (fourPi * radius**2 * stefanBoltz))**0.25d0
!                      
!
!                      if (associated(globalsourceArray)) then
!                         deallocate(globalSourceArray)
!                      endif
!
!                      globalNSource = 1
!                      allocate(globalsourceArray(1:1))
!
!
!                      call fillSpectrumBB(globalsourceArray(1)%spectrum, dble(teff), 120.d0, 1.d7, 100)
!                      call addToSpectrumBB(globalsourceArray(1)%spectrum, tAcc, 1.d0)
!                      globalsourceArray(1)%stellar = .true.
!                      globalsourceArray(1)%mass = massStar
!                      globalsourceArray(1)%radius = radius/1.d10
!                      globalsourceArray(1)%luminosity = lum + lAcc
!                      globalsourceArray(1)%teff = teff
!                      globalSourceArray(1)%position = VECTOR(0.d0, 0.d0, 0.d0)
!                      call testRandomSource(globalsourceArray, globalnsource)
!                      call writeSourceList(globalSourceArray, globalnSource)
!                      
!                      amrGridSize = rOuter
!                      amrgridCentreX = rOuter/2.d0
!                      call  setupamrgrid(grid)
!
!                      call writeVtkFile(grid, "rho.vtk")
!
!                      call randomNumberGenerator(randomSeed=.true.)
!
!                      call doPhysics(grid)
!
!                      write(radialFilename, &
!                           '(a,sp,f6.3,a,f6.3,a,f6.3,a,f6.3,a,f6.3,a,f6.3,a,f8.3,a)') &
!                           "menv_",log10(massEnvelope/msol), &
!                           "_rinner_",log10(rinner*1.d10/autocm), &
!                           "_router_",log10(router*1.d10/autocm), &
!                           "_mstar_",log10(massStar/msol), &
!                           "_rstar_",log10(radius/rsol), &
!                           "_lum_",log10(lum/lsol), &
!                           "_mdot_",log10(mdot*365.25d0*24.d0*3600.d0/mSol), &
!                           ".dat"
!
!                          
!
!                      call doOutputs(grid)
!
!
!                      call torus_mpi_barrier
!                      if (associated(grid%octreeRoot)) then
!                         call deleteOctreeBranch(grid%octreeRoot,onlyChildren=.false., adjustParent=.false.)
!                      else
!                         call writewarning("Want to delete grid but octreeroot not associated")
!                      endif
!
!
!                      call freeGrid(grid)
!                      call freeGlobalSourceArray()
!
!                   enddo
!                enddo
!             enddo
!          enddo
!       enddo
!    enddo
! end do
!
! 
!end subroutine calculategrid







end program torus
