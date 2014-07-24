module setupamr_mod

  Use cmfgen_class
  use amr_mod
  use vtk_mod, only: writeVTKFile, writeVTKFileSource
  use vector_mod
  use magnetic_mod
  use discwind_class
  use messages_mod
  USE constants_mod
  USE octal_mod, only: OCTAL, wrapperArray, octalWrapper, subcellCentre, cellVolume, &
       allocateattribute, copyattribute, deallocateattribute
  use gridtype_mod, only:   gridtype
  USE parallel_mod, ONLY:   torus_abort, torus_mpi_barrier
  use mpi_global_mod
  use TTauri_mod, only: TTauri_accretion_mass

  implicit none

contains
    
  subroutine setupamrgrid(grid)

    use gridio_mod
    use amr_mod
    use lucy_mod
    use grid_mod
    use inputs_mod, only : readgrid, gridinputfilename, geometry!, mdot
    use inputs_mod, only : amrGridCentreX, amrGridCentreY, amrGridCentreZ
    use inputs_mod, only : amr1d, amr2d, amr3d, splitOverMPI, atomicPhysics, molecularPhysics, variableDustSublimation
    use inputs_mod, only : amrGridSize, doSmoothGrid, ttauriMagnetosphere, discWind
    use inputs_mod, only : ttauriRstar, mDotparameter1, ttauriWind, ttauriDisc, ttauriWarp
    use inputs_mod, only : limitScalar, limitScalar2, smoothFactor, onekappa
    use inputs_mod, only : CMFGEN_rmin, CMFGEN_rmax, intextFilename
    use inputs_mod, only : rCore, rInner, rOuter, lamline,gridDistance, massEnvelope
    use inputs_mod, only : gridShuffle, minDepthAMR, maxDepthAMR, logspacegrid, nmag, dospiral
    use disc_class, only:  new
#ifdef ATOMIC
    use cmf_mod, only : checkVelocityInterp
#endif
    use inputs_mod, only : xplusbound, yplusbound, zplusbound
    use inputs_mod, only : xminusbound, yminusbound, zminusbound
#ifdef STATEQ
    use stateq_mod, only: map_cmfgen_opacities
#endif
#ifdef SPH
    use cluster_class
    use sph_data_class, only: read_sph_data_wrapper, deallocate_sph
#endif
#ifdef MPI 
    use mpi_amr_mod
    use inputs_mod, only : photoionPhysics, rho0
#ifdef PHOTOION
    use photoionAMR_mod, only : ionizeGrid, resetNh, resizePhotoionCoeff, &
         hasPhotoionAllocations, allocatePhotoionAttributes
#endif
#ifdef HYDRO
    use inputs_mod, only : hydrodynamics
#endif
#endif
    use vh1_mod
    use memory_mod
#ifdef USECFITSIO
    use gridFromFitsFile
#endif
    use gridFromFlash

    implicit none

    ! For romanova geometry case
    type(romanova) :: romData ! parameters and data for romanova geometry
    type(VECTOR) :: amrGridCentre
    type(GRIDTYPE) :: grid
    logical :: gridConverged
    real(double) :: astar, mass_accretion_old, totalMass
    real(double) :: minRCubedRhoSquared
    character(len=80) :: message
    integer :: nVoxels, nOctals
    integer :: counter, nTimes
!    integer :: nUnrefine
    real :: scalefac
#ifdef SPH
    type(cluster) :: young_cluster
    real(double)  ::  removedMass
    real(double) :: objectDistance
#endif
#ifdef MPI 
    integer :: i
#endif


#ifdef MPI
    call randomNumberGenerator(randomSeed=.true.)
    call randomNumberGenerator(syncIseed=.true.)
#endif

    call writeBanner("Setting up AMR grid","-",TRIVIAL)

    totalmass = 0.

    if ((xplusbound==6).or. &
         (yplusbound==6).or. &
         (zplusbound==6).or. &
         (xminusbound==6).or. &
         (yminusbound==6).or. &
         (zminusbound==6)) call setupInflowParameters()


    if (readgrid) then
       grid%splitOverMPI = splitOverMPI
       call readAMRgrid(gridInputfilename, .false., grid)
       grid%splitOverMPI = splitOverMPI
#ifdef MPI

       do i = 1, 8
          grid%octreeRoot%mpiThread(i) = i
       enddo

       call distributeMPIthreadLabels(grid%octreeRoot)

#ifdef PHOTOION

       if (photoionPhysics) then
          if (.not.hasPhotoionAllocations(grid)) then
             call allocatePhotoionAttributes(grid%octreeRoot, grid)
          endif
       endif

       if (photoIonPhysics) call resizePhotoionCoeff(grid%octreeRoot, grid)
#endif
#endif

       if (ttauriwind) call addGlobalDiscWind

       call findTotalMemory(grid, globalMemoryFootprint)

    else

       grid%lineEmission = atomicPhysics.or.molecularPhysics

       grid%splitOverMPI = splitOverMPI
       globalMemoryFootprint = 0

       select case (geometry)
          case("cmfgen")
             onekappa=.false.
             call read_cmfgen_data("OPACITY_DATA")
             call put_cmfgen_Rmin(CMFGEN_Rmin)  ! in [10^10cm]
             call put_cmfgen_Rmax(CMFGEN_Rmax)  ! in [10^10cm]
          case DEFAULT
       end select

       call initAMRGrid(grid)
       grid%splitOverMPI = splitOverMPI

       amrGridCentre = VECTOR(amrGridCentreX,amrGridCentreY,amrGridCentreZ)
       call writeInfo("Starting initial set up of adaptive grid...", TRIVIAL)
       if (writeoutput) write(*,*) "amr grid centre ", amrgridcentre

       select case (geometry)

       case("fogel")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call setupFogel(grid, intextFilename, "HCN")
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

       case("kengo")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call readgridKengo(grid)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

         
       case("turbulence")
#ifdef MPI
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .false.)
          call readGridTurbulence(grid)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
#endif

#ifdef SPH
       case("cluster")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .false.)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

! General case of setting up a grid from an SPH file
       case("molcluster","theGalaxy","sphfile", "dale")
          if (.not.grid%splitOverMPI) then
             call read_sph_data_wrapper
          else
#ifdef MPI
             if (myrankGlobal == 0) then
                call read_sph_data_wrapper
             endif
             call distributeSphDataOverMPI()
#endif
          endif
          call writeInfo("Initialising adaptive grid...", TRIVIAL)
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .false.)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)




       case("wr104")
          objectDistance = griddistance * pctocm
          call read_sph_data_wrapper
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
!          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .false.)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
          call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
          do
             gridConverged = .true.
             call myScaleSmooth(smoothFactor, grid, &
                  gridConverged,  inheritProps = .false., interpProps = .false.)
             if (gridConverged) exit
          end do
          call writeInfo("...grid smoothing complete", TRIVIAL)
#endif

#ifdef USECFITSIO
       case("fitsfile")
          call read_fits_file_for_grid()
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .false.)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
#endif

       case("clumpyagb")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
!          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .false.)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
          call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
          call writeInfo("...grid smoothing complete", TRIVIAL)

       case("starburst")
#ifdef MPI
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
!          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid)
          gridconverged = .false.
          grid%octreeRoot%rho = 100.*mhydrogen
          grid%octreeRoot%temperature = 10.
          call randomNumberGenerator(randomSeed = .true.)
          do while(.not.gridconverged) 
             gridConverged = .true.
             call splitGridFractal(grid%octreeRoot, rho0, 0.3, grid, gridconverged)
          enddo
          write(*,*) myrankglobal, " conv ",gridConverged
          call randomNumberGenerator(syncIseed=.true.)
          do i = 1, 8
             grid%octreeRoot%mpiThread(i) = i
          enddo
          call distributeMPIthreadLabels(grid%octreeRoot)

#ifdef PHOTOION
          if (photoionPhysics) then
             call ionizeGrid(grid%octreeRoot)
             call resetNH(grid%octreeRoot)
          endif
#endif

          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
#ifdef MPI
          call findMassOverAllThreads(grid, totalmass)
          write(message,'(a,1pe12.5,a)') "Total mass in fractal cloud (solar masses): ",totalMass/lsol
          call writeInfo(message,TRIVIAL)
#endif




#endif

       case("runaway")

          if (vh1FileRequired())   call read_vh1
          if (flashFileRequired()) call read_flash_hdf

          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d,  romData=romData) 
          call writeInfo("First octal initialized.", TRIVIAL)
!          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,romData=romData)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .false., romData=romData)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

          if (doSmoothGrid) then
             call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
             do
                gridConverged = .true.
                ! The following is Tim's replacement for soomthAMRgrid.
                call myScaleSmooth(smoothfactor, grid, &
                     gridConverged,  inheritProps = .false., &
                     interpProps = .false.)
                if (gridConverged) exit
             end do
             call writeInfo("...grid smoothing complete", TRIVIAL)
          endif

       case("ttauri")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, romData=romData) 
          call writeInfo("First octal initialized.", TRIVIAL)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .false., romData=romData)
          if (doSmoothGrid) then
          call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
          do
             gridConverged = .true.
             call myScaleSmooth(3., grid, &
                  gridConverged,  inheritProps = .false., interpProps = .false.)
             if (gridConverged) exit
          end do
          call writeInfo("...grid smoothing complete", TRIVIAL)
          endif
          call fixParentPointers(grid%octreeRoot)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

           ! This section is getting rather long. Maybe this should be done in 
           ! wrapper subroutine in amr_mod.f90.
          call zeroDensity(grid%octreeRoot)
          astar = accretingAreaMahdavi()
          if (writeoutput) write(*,*) "accreting area (%) ",100.*astar/(fourpi*ttauriRstar**2)
          minRCubedRhoSquared = 1.d30
          if (ttauriMagnetosphere) then
             call assignDensitiesMahdavi(grid, grid%octreeRoot, astar, mDotparameter1*mSol/(365.25d0*24.d0*3600.d0), &
                  minRCubedRhoSquared)
             call assignTemperaturesMahdavi(grid, grid%octreeRoot, astar, mDotparameter1*mSol/(365.25d0*24.d0*3600.d0), &
                     minRCubedRhoSquared)
          endif
!          if (ttauriwind) call assignDensitiesBlandfordPayne(grid, grid%octreeRoot)
          if (ttauriwind)  call addDiscWind(grid)
             
          call checkAMRgrid(grid, .false.)
#ifdef ATOMIC
          call writeInfo("Checking velocity interpolation...",TRIVIAL)
          call checkVelocityInterp(grid%octreeRoot, grid)
          call writeInfo("Done.",TRIVIAL)
#endif
          if (ttauridisc) call assignDensitiesAlphaDisc(grid, grid%octreeRoot)
          if (ttauriwarp) call addWarpedDisc(grid%octreeRoot)
          ! Finding the total mass in the accretion flow
          mass_accretion_old = 0.0d0
          call TTauri_accretion_mass(grid%octreeRoot, grid, mass_accretion_old)
          write(message,*) "Total mass in accretion flow is ",  mass_accretion_old, "[g]"
          call writeInfo(message,FORINFO)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)



       case DEFAULT
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, romData=romData) 
          call writeInfo("First octal initialized.", TRIVIAL)
!          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,romData=romData)
          call writeInfo("Doing first grid split.", TRIVIAL)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,wvars=.false., romData=romData)
          call writeInfo("grid split.", TRIVIAL)
          call fixParentPointers(grid%octreeRoot)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

          if (doSmoothGrid) then
             call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
             if (writeoutput) write(*,*) "Scale factor ",smoothFactor
             do
                gridConverged = .true.
                ! The following is Tim's replacement for soomthAMRgrid.
                call myScaleSmooth(smoothFactor, grid, &
                     gridConverged,  inheritProps = .false., &
                     interpProps = .false.)
                if (gridConverged) exit
             end do
             call writeInfo("...grid smoothing complete", TRIVIAL)
          endif

       end select




        call countVoxels(grid%octreeRoot,nOctals,nVoxels)
        grid%nOctals = nOctals
        call howmanysplits()

        call writeInfo("Calling routines to finalize the grid variables...",TRIVIAL)
        call fixParentPointers(grid%octreeRoot)
        call finishGrid(grid%octreeRoot, grid, romData=romData)
        call writeInfo("...final adaptive grid configuration complete",TRIVIAL)

        if(gridShuffle) then
           call writeInfo("Shuffling the grid for more accurate initial configuration...", TRIVIAL)
           nTimes = maxDepthAMR - minDepthAMR + 2
           if(logspacegrid) nTimes = nmag*10
           do counter = 1, nTimes
              call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, .true., romData=romData)     
              call fixParentPointers(grid%octreeRoot)
           end do
           
           if (doSmoothGrid) then
              call writeInfo("Smoothing adaptive grid structure 2...", TRIVIAL)
              do
                 gridConverged = .true.
                 ! The following is Tim's replacement for soomthAMRgrid.
                 call myScaleSmooth(smoothFactor, grid, &
                      gridConverged,  inheritProps = .false., &
                      interpProps = .false.)
                 if (gridConverged) exit
              end do
              call writeInfo("...grid smoothing complete", TRIVIAL)
           endif
           
           call countVoxels(grid%octreeRoot,nOctals,nVoxels)
           grid%nOctals = nOctals
           call howmanysplits()
           call writeInfo("Grid shuffle phase of initiation completed", TRIVIAL)          
           call finishGrid(grid%octreeRoot, grid, romData=romData)
#ifdef MPI
#ifdef HYDRO
           if(hydrodynamics) then
              call findMassOverAllThreads(grid, totalmass)
              write(message,'(a,1pe12.5,a)') "Total mass in grid (solar masses): ",totalMass/msol
              call writeInfo(message,TRIVIAL)
           endif
#endif
#endif

     end if
        


       select case (geometry)


          case("theGalaxy","fitsfile")
             if (variableDustsublimation) then
                call writeInfo("Calling quickSublimate",FORINFO)
                call quickSublimate(grid%octreeRoot)
             endif

          case("shakara")
             if (discWind) call fillgridSafier(grid)
             if (doSpiral) call addSpiralWake(grid)
          case("ttauri")
             if (ttauriMagnetosphere) then
                call assignDensitiesMahdavi(grid, grid%octreeRoot, astar, mDotparameter1*mSol/(365.25d0*24.d0*3600.d0), &
                     minRcubedRHoSquared)
                call assignTemperaturesMahdavi(grid, grid%octreeRoot, astar, mDotparameter1*mSol/(365.25d0*24.d0*3600.d0), &
                     minRCubedRhoSquared)
             endif
             call testVelocity(grid%octreeRoot,grid)
       
             if (ttauriwind)  call addDiscWind(grid)
             if (ttauridisc) call assignDensitiesAlphaDisc(grid, grid%octreeRoot)
             if (ttauriwarp) call addWarpedDisc(grid%octreeRoot)

             call writeVtkFile(grid, "deriv.vtk",  valueTypeString=(/"rho        ", &
!                                                                     "dust1      ", &
                                                                     "etaline    ", &
                                                                     "chiline    ", &
                                                                     "inflow     ", &
                                                                     "cornervel  ", &
                                                                     "velcall    ", &
                                                                     "vely       ", &
                                                                     "temperature", &
                                                                     "fixedtemp  "/))
#ifdef STATEQ
          case("cmfgen")
              call map_cmfgen_opacities(grid)
              call distort_cmfgen(grid%octreeRoot, grid)
              call set_bias_cmfgen(grid%octreeRoot, grid, lamline)
              call writeVtkFile(grid, "cmfgen.vtk",  valueTypeString=(/"etaline ", &
                                                                       "chiline ", &
                                                                       "ne      ", &
                                                                       "velocity"/))
#endif

          case("wrshell")
             grid%geometry = "wrshell"
             grid%lineEmission = .false.
             grid%rInner = rInner
             grid%rOuter = rOuter
             grid%rCore = rCore

#ifdef SPH
          case("cluster", "dale")
             removedMass = 0.0
             call new(young_cluster, .false.)
             call remove_too_close_cells(young_cluster, grid%octreeRoot, real(rCore,db), removedMass, 1.0d-37, 'z')
             call kill_all(young_cluster)
             write(message,*) "Mass removed by remove_too_close_cells: ", removedMass / mSol, " solar masses"
             call writeInfo(message, TRIVIAL)


          case("wr104")
             totalMass = 0.d0
             call findTotalMass(grid%octreeRoot, totalMass)
             write(*,*) "mass envelope ",massEnvelope
             write(*,*) "total mass ",totalmass
             scaleFac = real(massEnvelope / totalMass)
             if (writeoutput) write(*,'(a,1pe12.5)') "Density scale factor: ",scaleFac
             call scaleDensityAMR(grid%octreeRoot, dble(scaleFac))
#endif
          case("envelope")
             totalMass = 0.d0
             call findTotalMass(grid%octreeRoot, totalMass)
             if (writeoutput) write(*,*) "Mass of envelope: ",totalMass/msol

          case("spiral")
             totalMass = 0.d0
             call findTotalMass(grid%octreeRoot, totalMass)
             scaleFac = real(massEnvelope / totalMass)
             call scaleDensityAMR(grid%octreeRoot, dble(scaleFac))
             totalMass = 0.d0
             call findTotalMass(grid%octreeRoot, totalMass)
             if (writeoutput) write(*,*) "Total mass in spiral (solar masses): ",totalMass/msol
          case("turbbox")
             call turbulentVelocityField(grid, 1.d0)
#ifdef MPI
          case("unisphere","gravtest")
             if (hydrodynamics) then
             totalMass = 0.d0
                call findmassoverallthreads(grid, totalmass)
                if (writeoutput) write(*,*) "Total mass in sphere before (solar masses): ",totalMass/msol
                call  torus_mpi_barrier
                if (myrankGlobal /= 0) call scaleDensityAMR(grid%octreeRoot, 2.d0*msol/totalMass)
                call  torus_mpi_barrier
                totalMass = 0.d0
                call findmassoverallthreads(grid, totalmass)
                if (writeoutput) write(*,*) "Total mass in sphere (solar masses): ",totalMass/msol
             endif
#endif
          case DEFAULT
       end select

       call fixParentPointers(grid%octreeRoot)
       call postSetupChecks(grid)
    endif
#ifdef MPI
        if (grid%splitOverMPI) then
           call grid_info_mpi(grid, "info_grid.dat")
           !Check an appropriate no. of MPI threads is being used
           call checkThreadNumber(grid)

           if(readGrid) then
              do i = 1, 8
                 grid%octreeRoot%mpiThread(i) = i
              enddo
              
              !label each cell with its appropriate MPI thread
!              write(*,*) "Distributing MPI Labels"
              call distributeMPIthreadLabels(grid%octreeRoot)
!              write(*,*) "Label Distribution Completed"
           end if
        else
           if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
        endif
        call torus_mpi_barrier

#else
        if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
#endif
               grid%adaptive = .true.


!        pause
!       gridConverged = .false.
!       do while(.not.gridconverged)
!          gridconverged = .true.
!          nUnrefine = 0
!          call unrefineBack(grid%octreeRoot, grid, 1.125d0, 100.d0*autocm/1.d10, dble(grid%rinner), &
!               nUnrefine, gridconverged)
!          if (writeoutput) write(*,*) "Unrefined ",nUnrefine, " cells on this pass"
!       end do
!       pause
!        call bigarraytest(grid%octreeRoot)
!    call findTotalMemory(grid, i)
!    call reportMemory(i)
!        if (myrankGlobal == 1) then 
!           call freeGrid(grid)
!           do;enddo
!           else
!              do; enddo
!              endif

! Free allocatable arrays used in the grid set up which are no longer required. 
! These subroutines should check that the arrays they are deallocating have actually 
! been allocated, that way it is always safe to call the subroutine.
#ifdef USECFITSIO
        call deallocate_gridFromFitsFile
#endif
        call deallocate_gridFromFlash
        call deallocate_vh1
        call delete_particle_lists(grid%octreeRoot)

  end subroutine setupamrgrid

  subroutine doSmoothOnTau(grid)

    use inputs_mod, only: doSmoothGridTau, dustPhysics, lambdaSmooth!, cylindrical
    use inputs_mod, only: dosmoothgrid, smoothfactor !,  photoionPhysics, variableDustSublimation, 
    use utils_mod, only: locate
    use lucy_mod, only: putTau
    use grid_mod, only: grid_info

    implicit none

    type(GRIDTYPE) :: grid
    integer :: j, ismoothlam, iter
    logical :: gridConverged, doPhotoSphereSplit, smoothConverged, convergedFirstTime
    character(len=80) :: message, thisFile


    ! Smooth the grid with respect to optical depth, if requested
    if (doSmoothGridTau.and.dustPhysics) then

       call writeVtkFile(grid, "beforesmooth.vtk")

       dophotosphereSplit = .false.
!       if (grid%geometry == "shakara") doPhotosphereSplit = .true.
!       doPhotoSphereSplit =  ((.not.cylindrical).and.(.not.variableDustSublimation).and.(.not.photoionPhysics))


      if ( doPhotoSphereSplit ) then
          call writeInfo("Doing photsphere split", TRIVIAL)
       else
          call writeInfo("Not doing photsphere split", TRIVIAL)
       end if

       call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
       call locate(grid%lamArray, grid%nLambda,lambdaSmooth,ismoothlam)
       iter = 1
       do j = iSmoothLam,  grid%nLambda, 5
          write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
          call writeInfo(message, TRIVIAL)
          do
             gridConverged = .true.
             ! The following is Tim's replacement for soomthAMRgrid.
             call myScaleSmooth(smoothFactor, grid, &
                  gridConverged,  inheritProps = .false., interpProps = .false.)
             if (gridConverged) exit
          end do

          do
             gridConverged = .true.
             if (dophotospheresplit) call putTau(grid, grid%lamArray(j))
             call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                  inheritProps = .false., interpProps = .false., &
                  photosphereSplit = doPhotoSphereSplit )
             write(thisFile,'(a,i5.5,a)') "smooth",iter,".vtk"
             iter = iter + 1
             call writeVtkFile(grid, thisFile)
             convergedFirstTime = .true.
             do
                smoothConverged = .true.
                ! The following is Tim's replacement for soomthAMRgrid.
                call myScaleSmooth(smoothFactor, grid, &
                     smoothConverged,  inheritProps = .false., interpProps = .false.)
                if (.not.smoothConverged) convergedFirstTime = .false.
                if (smoothConverged) exit
             end do
             if (gridConverged.and.convergedFirstTime) exit
          enddo
       enddo
       call writeInfo("...grid smoothing complete", TRIVIAL)

       ! The tau smoothing may result in large differences in the size
       ! of neighbouring octals, so we smooth the grid again.

       if (doSmoothgrid) then
          call writeInfo("Smoothing adaptive grid structure (again)...", TRIVIAL)
          do
             gridConverged = .true.
             ! The following is Tim's replacement for soomthAMRgrid.
             call myScaleSmooth(smoothFactor, grid, &
                  gridConverged,  inheritProps = .false., interpProps = .false.)
             if (gridConverged) exit
          end do
          call writeInfo("...grid smoothing complete", TRIVIAL)
       endif
       call writeVtkFile(grid, "aftersmooth.vtk")
             
       if ( myRankIsZero ) call grid_info(grid, "info_grid_aftersmooth.dat")

    end if

  end subroutine doSmoothOnTau

  subroutine setupFogel(grid, filename, speciesName)
    use inputs_mod, only : rinner, rOuter, molecularPhysics
    type(GRIDTYPE) :: grid
    character(len=*) :: filename, speciesName
    integer, parameter :: maxR = 200, maxZ = 200
    integer :: nr
    integer, allocatable :: nz(:)
    real(double), allocatable :: r(:), z(:,:), rho(:,:), t(:,:)
    real(double), allocatable :: abundance(:,:)
    integer :: iSpecies
    character(len=80) :: cJunk
    real(double) :: junk(5)
    real(double), allocatable :: rev(:)
    integer :: i, j
    logical :: gridConverged

    rinner = real(5.*autocm/1.d10)
    grid%rinner = rinner
    router = real(400.*autocm/1.d10)
    allocate(r(maxR), nz(maxR), z(maxR, maxZ), rho(maxR, maxZ), &
         t(maxR,maxZ), abundance(maxR,maxZ))

    select case (speciesName)
       case("CN")
          iSpecies = 1
       case("HCN")
          iSpecies = 2
       case("C")
          iSpecies = 3
       case("C+")
          iSpecies = 4
       case("H2CO")
          iSpecies = 5
       case DEFAULT
          call writeFatal("setupFogel: unknown species")
    end select
    

    open(20, file=filename, status="old", form="formatted")
    read(20,'(A)') cJunk
    read(20,'(A)') cJunk
    
    nr = 1
    nz(1) = 0

10  continue
    read(20,'(a)',end=30) cJunk
    if (len(trim(cJunk)) < 10) then ! blank line
       read(20,'(a)',end=30) cJunk
       read(20,'(a)',end=30) cJunk
       nr = nr + 1
       nz(nr) = 0
       goto 10
    endif

    nz(nr) = nz(nr) + 1
    read(cJunk,*) r(nr), z(nr,nz(nr)), rho(nr,nz(nr)), t(nr,nz(nr)), junk(1:5)
    abundance(nr, nz(nr)) = junk(iSpecies)
    goto 10
30  continue
    close(20)
    nr = nr - 1
    r = r * autocm/1.d10
    z = z * autocm/1.d10

    do i = 1, nr
       allocate(rev(1:nz(i)))
       do j = 1, nz(i)
          rev(j) = z(i,nz(i) - j + 1)
       enddo
       z(i,1:nz(i)) = rev(1:nz(i))

       do j = 1, nz(i)
          rev(j) = rho(i,nz(i) - j + 1)
       enddo
       rho(i,1:nz(i)) = rev(1:nz(i))

       do j = 1, nz(i)
          rev(j) = t(i,nz(i) - j + 1)
       enddo
       t(i,1:nz(i)) = rev(1:nz(i))

       do j = 1, nz(i)
          rev(j) = abundance(i,nz(i) - j + 1)
       enddo
       abundance(i,1:nz(i)) = rev(1:nz(i))
       deallocate(rev)
    enddo



    call splitGridFogel(grid%octreeRoot, grid, r, z, nr, nz, rho, t, abundance)
    do
       gridConverged = .true.
       call myScaleSmooth(3., grid, &
            gridConverged,  inheritProps = .true., interpProps = .false.)
       if (gridConverged) exit
    end do
!    call writeWarning("Need to check whether abundances are really relative to N(H) rather than N(H_2)")
    call fillGridFogel(grid%octreeRoot, grid, r, z, nr, nz, rho, t, abundance)

    if (molecularPhysics) then
       call writeVtkFile(grid, "fogel.vtk",  valueTypeString=(/"rho         ",&
            "temperature ","molabundance","microturb   ","velocity    "/))
    else
       call writeVtkFile(grid, "fogel.vtk",  valueTypeString=(/"rho         ",&
         "temperature ","velocity    "/))
    endif
  end subroutine setupFogel

  subroutine writeFogel(grid, infilename, outfilename)
    type(GRIDTYPE) :: grid
    character(len=*) :: infilename, outfilename
    integer, parameter :: maxR = 200, maxZ = 200
    integer :: nr
    integer, allocatable :: nz(:)
    real(double), allocatable :: r(:), z(:,:), rho(:,:), t(:,:)
    real(double), allocatable :: abundance(:,:)
    character(len=120) :: cJunk
    real(double) :: junk(5)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: rTemp, zTemp
    type(VECTOR) :: rVec

    allocate(r(maxR), nz(maxR), z(maxR, maxZ), rho(maxR, maxZ), &
         t(maxR,maxZ), abundance(maxR,maxZ))
    thisOctal => grid%octreeRoot

    open(20, file=infilename, status="old", form="formatted")
    open(21, file=outfilename, status="unknown", form="formatted")
    read(20,'(A)') cJunk
    write(21,'(A)') trim(cjunk)
    read(20,'(A)') cJunk
    write(21,'(A)') trim(cjunk)
    
    nr = 1
    nz(1) = 0

10  continue
    read(20,'(a)',end=30) cJunk
    if (len(trim(cJunk)) < 10) then ! blank line
       write(21,'(A)') trim(cjunk)
       read(20,'(a)',end=30) cJunk
       write(21,'(A)') trim(cjunk)
       read(20,'(a)',end=30) cJunk
       write(21,'(A)') trim(cjunk)
       nr = nr + 1
       nz(nr) = 0
       goto 10
    endif

    nz(nr) = nz(nr) + 1
    read(cJunk,*) r(nr), z(nr,nz(nr)), rho(nr,nz(nr)), t(nr,nz(nr)), junk(1:5)
    rtemp = r(nr)*autocm/1.d10
    ztemp = z(nr, nz(nr))*autocm/1.d10
    rVec = VECTOR(rTemp, 0.d0, zTemp)
    call findSubcellLocal(rvec,thisOctal,subcell)
    t(nr, nz(nr)) = thisOctal%temperature(subcell)
    write(21, '(F5.1, f7.2, 1p, e11.3, 0p, f7.1, 1p, 5e11.3)')  r(nr), z(nr,nz(nr)), rho(nr,nz(nr)), &
         t(nr,nz(nr)), junk(1:5)


    goto 10
30  continue
    close(20)
    close(22)


  end subroutine writeFogel

  subroutine unrefineGridOnTau(grid)
    type(GRIDTYPE) :: grid
    logical :: converged
    integer :: i
    i = 10
    call writeInfo("Unrefining optically thin cells...")
    do
       converged = .true.
       call unrefineThinCells(grid%octreeRoot, grid, i, converged)
       if (converged) exit
    end do
    call writeInfo("Done.")
  end subroutine unrefineGridOnTau


  recursive subroutine splitGridFogel(thisOctal, grid, r, z, nr, nz, rho, t, abundance)
    use inputs_mod, only : minDepthAMR, maxDepthAMR
    type(GRIDTYPE) :: grid
    type(OCTAL),pointer :: thisOctal !TJH 9 JULY
    type(OCTAL), pointer :: childPointer
!    type(OCTAL) :: thisOctal
    integer :: iIndex
    real(double) :: r(:), z(:,:), rho(:,:), t(:,:), abundance(:,:)
    integer :: nr, nz(:)
    logical :: splitInAzimuth
    logical :: split
    integer :: iSubcell, i,j,k
    type(VECTOR) :: rVec
    real(double) :: thisR, thisZ,s
    logical :: outsideGrid

    splitInAzimuth = .false.


    DO iSubcell = 1, thisOctal%maxChildren

       split = .false.
       rVec = subcellCentre(thisOctal, iSubcell)
       s = thisOctal%subcellSize/2.d0
       thisR = sqrt(rVec%x**2 + rVec%y**2)
       thisZ = abs(rVec%z)
       outsideGrid = .false.
       if ((thisR+s < r(1)).or.(thisR-s > r(nr))) outsideGrid = .true.
       if (.not.outsideGrid) then
          call locate(r, nr, thisR, i)
          if (thisZ-s > z(i,nz(i))) outsidegrid = .true.
       endif
       if (.not.OutsideGrid) then
          call locate(z(i,1:nz(i)), nz(i), thisZ, j)
          if (2.d0*thisOctal%subcellSize > (r(i+1)-r(i))) split = .true.
          if (thisOctal%subcellSize > (z(i,j+1)-z(i,j))) then
             split = .true.
          endif
       endif
       if (thisOctal%nDepth >= maxDepthAMR) split = .false.
       if (thisOctal%nDepth <= minDepthAMR) split = .true.
       IF (split) then

          CALL addNewChild(thisOctal, iSubcell, grid, adjustGridInfo=.TRUE., &
               inherit = .true., &
               splitAzimuthally=splitInAzimuth)

          if (.not.thisOctal%hasChild(isubcell)) then
             write(*,*) "add child failed in splitGrid"
             do
             enddo
          endif



       END IF

    END DO

    do i = 1, thisOctal%nChildren
       if (.not.thisOctal%hasChild(thisOctal%indexchild(i))) then
          write(*,*) "octal children messed up"
          do ; enddo
          endif
       enddo

       do i = 1, thisOctal%maxChildren
          k = -99
          if (thisOctal%hasChild(i)) then
             do j = 1, thisOctal%nChildren
                if (thisOctal%indexChild(j) == i) then
                   k = j
                   exit
                endif
             enddo
             if (k==-99) then
                write(*,*) "subcell screwup"
                do
                enddo
             endif
          endif
       enddo

       if (any(thisOctal%haschild(1:thisOctal%maxChildren)).and.(thisOctal%nChildren==0)) then
          write(*,*) "nchildren screw up"
          do;enddo
          endif

          DO iIndex = 1, thisOctal%nChildren
             childPointer => thisOctal%child(iIndex) ! TJH 9 JULY
             call  splitGridFogel(childPointer, grid, r, z, nr, nz, rho, t, abundance)
!             call  splitGridFogel(thisOctal%child(iIndex), grid, r, z, nr, nz, rho, t, abundance)

          END DO

        end subroutine splitGridFogel



        recursive subroutine fillGridFogel(thisOctal, grid, r, z, nr, nz, rho, t, abundance)
          use inputs_mod, only : mcore, vturb, atomicPhysics, molecularPhysics,sourcemass
          type(GRIDTYPE) :: grid
          type(octal), pointer   :: thisOctal
          type(octal), pointer  :: child 
          type(VECTOR) :: rVec
          real(double) :: thisR, thisZ,fac1,fac2,fac3,s
          integer :: subcell, i, j, k1,k2
          logical :: outsideGrid
          real(double) :: r(:), z(:,:), rho(:,:), t(:,:), abundance(:,:)
          integer :: nr, nz(:)

          do subcell = 1, thisOctal%maxChildren
             if (thisOctal%hasChild(subcell)) then
                ! find the child
                do i = 1, thisOctal%nChildren, 1
                   if (thisOctal%indexChild(i) == subcell) then
                      child => thisOctal%child(i)
                      call fillGridFogel(child, grid, r, z, nr, nz, rho, t, abundance)
                      exit
                   end if
                end do
             else

                Mcore = real(sourcemass(1))
                rVec = subcellCentre(thisOctal, subcell)
                thisR = sqrt(rVec%x**2 + rVec%y**2)
                thisZ = abs(rVec%z)
                s = thisOctal%subcellsize/2.d0
                thisOctal%velocity(subcell) = keplerianVelocity(rvec)

                if (atomicPhysics.or.molecularPhysics) CALL fillVelocityCorners(thisOctal,keplerianVelocity)
                outsideGrid = .false.

                if ((thisR+s < r(1)).or.(thisR-s > r(nr))) outsideGrid = .true.
                if (.not.outsideGrid) then
                   call locate(r, nr, thisR, j)
                   if ((thisZ < z(j,1)).or.(thisZ > z(j+1,nz(j+1)))) outsideGrid = .true.
                endif

                thisOctal%rho(subcell) = 1.d-25
                thisOctal%temperature(subcell) = 3.d0
                if (molecularPhysics) then
                   thisOctal%molAbundance(subcell) = 1.e-20
                   thisOctal%microTurb(subcell) = 0.d0
                endif
                if (.not.outsideGrid) then
                   call locate(z(j,1:nz(j)), nz(j), thisZ, k1)
                   call locate(z(j+1,1:nz(j+1)), nz(j+1), thisZ, k2)
                   fac1 = (thisR - r(j))/(r(j+1)-r(j))
                   fac2 = (thisZ - z(j,k1))/(z(j,k1+1)-z(j,k1))
                   fac3 = (thisZ - z(j+1,k2))/(z(j+1,k2+1)-z(j+1,k2))

                   thisOctal%rho(subcell) = (1.d0-fac1)*(1.d0-fac2)* rho(j,k1) + &
                        (     fac1)*(1.d0-fac3)* rho(j+1,k2) + &
                        (1.d0-fac1)*(     fac2)* rho(j,k1+1) + &
                        (     fac1)*(     fac3)* rho(j+1,k2+1)  


                   thisOctal%temperature(subcell) = real(max(3.d0,(1.d0-fac1)*(1.d0-fac2)* t(j,k1) + &
                        (     fac1)*(1.d0-fac3)* t(j+1,k2) + &
                        (1.d0-fac1)*(     fac2)* t(j,k1+1) + &
                        (     fac1)*(     fac3)* t(j+1,k2+1) ) )


                   if (molecularPhysics) then
                      thisOctal%molAbundance(subcell) = real((1.d0-fac1)*(1.d0-fac2)* abundance(j,k1) + &
                           (     fac1)*(1.d0-fac3)* abundance(j+1,k2) + &
                           (1.d0-fac1)*(     fac2)* abundance(j,k1+1) + &
                           (     fac1)*(     fac3)* abundance(j+1,k2+1)  )

                      thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.d0 * mHydrogen) 

                      thisOctal%molAbundance(subcell) = real(thisOctal%molAbundance(subcell) * 2.d0)
                   endif
                endif
                if (molecularPhysics) then
                   thisOctal%microturb(subcell) = max(vturb/(1.d-5*cspeed),sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) &
                        / (28.0 * amu)) + vturb**2) &
                        / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence
                   
                   thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.d0 * mHydrogen) 
                endif
             endif
          enddo
        end subroutine fillGridFogel


#ifdef MPI
!routine to read in the velocity field produced by the generating code of M. R. Bate.
      subroutine readgridTurbulence(grid)
      use mpi
      use inputs_mod, only: turbvelfilex, turbvelfiley, turbvelfilez, readTurb, nturblines
      type(GRIDTYPE) :: grid
      type(OCTAL), pointer :: thisOctal
      integer :: subcell
      integer :: i, ierr
      real(double), allocatable :: u(:), v(:), w(:)
      real(double), allocatable :: x(:), y(:), z(:)
      real(double) ::  dx, range,  max, min
      type(VECTOR) :: locator, cen, minusBound
      integer :: myRank


      call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
      if(readTurb) then

         allocate(u(nturblines**3))
         allocate(v(nturblines**3))
         allocate(w(nturblines**3))
         allocate(x(nturblines**3))
         allocate(y(nturblines**3))
         allocate(z(nturblines**3))
         
         dx = 2.d0*grid%octreeRoot%subcellSize/dble(nturblines)

         cen = grid%octreeRoot%centre

         minusBound%x = cen%x - grid%octreeRoot%subcellSize + dx/2.d0!grid%halfsmallestsubcell
         minusBound%y = cen%y - grid%octreeRoot%subcellSize + dx/2.d0!grid%halfsmallestsubcell
         Minusbound%z = cen%z - grid%octreeRoot%subcellSize + dx/2.d0!grid%halfsmallestsubcell

!         print *, "dx ", dx
!         print *, "cen ", cen
!         print *, "minusBound%x", minusBound%x
!         print *, "minusBound%y", minusBound%y
!         print *, "minusBound%z", minusBound%z

!         print *, "importing turbulent velocity field", dx
         open(1, file=turbvelfilex, status="old", iostat=ierr)
         if(ierr /= 0) then
            print *, "Trouble opening " //turbvelfilex
         end if
         open(2, file=turbvelfiley, status="old", iostat=ierr)
         if(ierr /= 0) then
            print *, "Trouble opening " //turbvelfiley
         end if
         open(3, file=turbvelfilez, status="old", iostat=ierr)
         if(ierr /= 0) then
            print *, "Trouble opening " //turbvelfilez
         end if

 !        print *, "reading lines", nTurblines
         thisOctal => grid%octreeRoot
         
         do i=1, (nTurblines**3)
            read(1, *, iostat=ierr) u(i), x(i), y(i), z(i)
            if(ierr /= 0) then
               print *, "error reading from " //turbvelfilex
            end if
            read(2, *, iostat=ierr) v(i), x(i), y(i), z(i)
            if(ierr /= 0) then
               print *, "error reading from " //turbvelfiley
            end if
            read(3, *, iostat=ierr) w(i), x(i), y(i), z(i)
            if(ierr /= 0) then
               print *, "error reading from " //turbvelfilez
            end if
        end do

        
        if(maxval(abs(u)) > maxval(abs(v))) then
           if(maxval(abs(u)) > maxval(abs(w))) then
              max = maxval(abs(u))
           else
              max = maxval(abs(w))
           end if
        else
           if(maxval(abs(v)) > maxval(abs(w))) then
              max = maxval(abs(v))
           else
              max = maxval(abs(w))
           end if
        end if


        if(minval(abs(u)) < minval(abs(v))) then
           if(minval(abs(u)) < minval(abs(w))) then
              min = minval(abs(u))
           else
              min = minval(abs(w))
           end if
        else
           if(minval(abs(v)) < minval(abs(w))) then
              min = minval(abs(v))
           else
              min = minval(abs(w))
           end if
        end if

        range = max-min

        do i = 1, nturblines**3

           locator = VECTOR(minusBound%x + dble(x(i)-1)*dx, minusBound%y + &
           dble(y(i)-1)*dx, minusBound%z + dble(z(i)-1)*dx)
          
           call findSubcellLocal(locator, thisOctal, subcell)
           if(octalOnThread(thisOctal, subcell, myRank)) then
                     
              thisOctal%rhou(subcell) = thisOctal%rho(subcell)*(u(i)/max)*10.d0*sqrt(thisOctal%rho(subcell)* &
                 thisOctal%energy(subcell))
              thisOctal%rhov(subcell) = thisOctal%rho(subcell)*(v(i)/max)*10.d0*sqrt(thisOctal%rho(subcell)* &
                 thisOctal%energy(subcell))
              thisOctal%rhow(subcell) = thisOctal%rho(subcell)*(w(i)/max)*10.d0*sqrt(thisOCtal%rho(subcell)* &
                 thisOctal%energy(subcell))
              thisOctal%velocity(subcell) = VECTOR((u(i)/max)*10.d0*sqrt(thisOctal%rho(subcell)*thisOctal%energy(subcell)), &
              (v(i)/max)*10.d0*sqrt(thisOctal%rho(subcell)*thisOctal%energy(subcell)), &
              (w(i)/max)*10.d0*sqrt(thisOctal%rho(subcell)*thisOctal%energy(subcell)))
              
                           
              thisOctal%energy(subcell) = thisOctal%energy(subcell)+ 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
              
              thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
           end if
        end do

         deallocate(u)
         deallocate(v)
         deallocate(w)
         deallocate(x)
         deallocate(y)
         deallocate(z)

      end if


      end subroutine readgridTurbulence
#endif


        subroutine readgridKengo(grid)
          type(GRIDTYPE) :: grid
! *** data read routine
! *** input: step, levmin, levmax
! *** output: x, y, z, rho
!      subroutine readd(step,levmin,levmax,x,y,z,rho)

! *** parameters for grid configuration
! *** imax,jmax,kmax: grid size, lmax: max of levels
      integer, parameter :: imax=64, jmax=64, kmax=64, lmax=10
! *** ibl,jbl,kbl: offset between coarse and fine levels
! *** relation between coarse and fine levels is:
! *** rho(ibl+i,jbl+j,kbl+k,l-1)=
! *** (rho(2*i-1, 2*j-1, 2*k-1, l)+...+rho(2*i, 2*j, 2*k, l))/8d0
!      integer, parameter :: ibl=16, ibr=48, jbl=16, jbr=48, kbl=16, kbr=48

! *** conversion factors from computational units to cgs units.
      real(double), parameter :: rho0=1.104d-17, l0=6.3728748d15 

! *** input/output variables
      integer levmin,levmax
      real(double) ::  x(imax), y(jmax) ,z(kmax)
      real(double) ::  rho(imax,jmax,kmax)


! *** internal variables
      real(double) ::  time(lmax)
      integer(kind=8) ::  lstep(lmax)
      real :: version, arrtmp(imax,jmax,kmax)
      real :: xtmp(imax),ytmp(jmax),ztmp(kmax)
      integer :: i, j, k, l, level
      integer :: subcell
      character(len=80) :: fn
      logical :: converged
      type(VECTOR) :: point
      type(OCTAL), pointer :: thisOctal
      levmin = 5
      levmax = 9
      do l=levmin,levmax
         write(fn,'(a,i1,a)') "st798128.",l,".d"
         open(11,file=fn,form='unformatted')
         read(11) version,i,j,k,level,lstep(l),time(l)
         write(*,*) version, i, j, k,level,lstep(l),time(l)
         read(11) xtmp
         if (Writeoutput)write(*,*) xtmp
         read(11) ytmp
         read(11) ztmp
         do i=1,imax
            x(i)=xtmp(i)*l0/1.d10
         enddo
         do j=1,jmax
            y(j)=ytmp(j)*l0/1.d10
         enddo
         do k=1,kmax
            z(k)=ztmp(k)*l0/1.d10
         enddo
         if (Writeoutput) write(*,*) "Grid size is ",(x(64)-x(1))+(x(2)-x(1)),2.d0*(x(64)-x(1))+(x(2)-x(1))
         read(11) arrtmp
         do k=1,kmax
            do j=1,jmax
                do i=1,imax
                  rho(i,j,k)=arrtmp(i,j,k)*rho0
               enddo
            enddo
         enddo
         close(11)
         converged = .false.
         thisOctal => grid%octreeRoot
         do while (.not.converged)
            converged = .true.
            do k = 1, kmax
               do j = 1, jmax
                  do i = 1, imax
                     point = VECTOR( dble(x(i)), dble(y(j)), dble(z(k)) )
                     call findSubcellLocal(point,thisOctal,subcell)
                     if (thisOctal%subcellSize > (x(2)-x(1))) then
                        converged = .false.
                        call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                             inherit=.true., interp=.false.)
                     endif
                  enddo
               enddo
            enddo
         enddo


         do i = 1, imax
            do j = 1, jmax
               do k = 1, kmax
                  point = VECTOR( dble(x(i)), dble(y(j)), dble(z(k)) )
                  call findSubcellLocal(point,thisOctal,subcell)
                  thisOctal%rho(subcell) = rho(i,j,k)
               enddo
            enddo
         enddo


      enddo
    end subroutine readgridKengo


        subroutine writegridKengo(grid)
          type(GRIDTYPE) :: grid
! *** data read routine
! *** input: step, levmin, levmax
! *** output: x, y, z, rho
!      subroutine readd(step,levmin,levmax,x,y,z,rho)

! *** parameters for grid configuration
! *** imax,jmax,kmax: grid size, lmax: max of levels
      integer, parameter :: imax=64, jmax=64, kmax=32, lmax=20
! *** ibl,jbl,kbl: offset between coarse and fine levels
! *** relation between coarse and fine levels is:
! *** rho(ibl+i,jbl+j,kbl+k,l-1)=
! *** (rho(2*i-1, 2*j-1, 2*k-1, l)+...+rho(2*i, 2*j, 2*k, l))/8d0
!      integer, parameter :: ibl=16, ibr=48, jbl=16, jbr=48, kbl=16, kbr=48

! *** conversion factors from computational units to cgs units.
      real(double), parameter :: l0=6.3728748d15

! *** input/output variables
      integer levmin,levmax
      real(double) ::  x(imax), y(jmax) ,z(kmax)


! *** internal variables
      real(double) ::  time(lmax)
      integer(kind=8) ::  lstep(lmax)
      real :: version, arrtmp(imax,jmax,kmax)
      real :: xtmp(imax),ytmp(jmax),ztmp(kmax)
      integer :: i, j, k, l, level
      integer :: subcell
      character(len=80) :: fn
      type(VECTOR) :: point
      type(OCTAL), pointer :: thisOctal
      levmin = 5
      levmax = 9
      do l=levmin,levmax
         write(fn,'(a,i1,a)') "st798128.",l,".d"
         open(11,file=fn,form='unformatted')
         write(fn,'(a,i1,a)') "out798128.",l,".d"
         open(12,file=fn,form='unformatted')
         read(11) version,i,j,k,level,lstep(l),time(l)
         write(12) version,i,j,k,level,lstep(l),time(l)
         write(*,*) version, i, j, k,level,lstep(l),time(l)
         read(11) xtmp
         read(11) ytmp
         read(11) ztmp
         write(12) xtmp
         write(12) ytmp
         write(12) ztmp
         do i=1,imax
            x(i)=xtmp(i)*l0/1.d10
         enddo
         do j=1,jmax
            y(j)=ytmp(j)*l0/1.d10
         enddo
         do k=1,kmax
            z(k)=ztmp(k)*l0/1.d10
         enddo
         if (Writeoutput) write(*,*) "Grid size is ",(x(64)-x(1))+(x(2)-x(1)),2.d0*(x(64)-x(1))+(x(2)-x(1))
         read(11) arrtmp
         write(12) arrtmp

         thisOctal => grid%octreeRoot
         subcell = 1
         do i = 1, imax
            do j = 1, jmax
               do k = 1, kmax
                  point = VECTOR( dble(x(i)), dble(y(j)), dble(z(k)) )
                  call findSubcellLocal(point,thisOctal,subcell)
                  arrtmp(i,j,k) = thisOctal%temperature(subcell)
               enddo
            enddo
         enddo


         write(12) arrtmp

         close(11)
         close(12)

      enddo
    end subroutine writegridKengo




  recursive subroutine splitGridFractal(thisOctal, rho, aFac, grid, converged)
    use inputs_mod, only : photoionPhysics, hydrodynamics, inputNSource, &
         sourcepos, smallestCellSize, maxDepthAmr
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    type(VECTOR) :: centre
    integer :: isource
    real :: rho, aFac
    integer :: subcell, i, j
    real, allocatable :: r(:), s(:)
    real(double) :: ethermal, r1
    real :: rmin, rmax, tot, fac, mean
    logical :: converged, split

    if (.not.converged) return
    converged = .true.

! based on method 2 of Hetem and Lepine 1993 A&A 270 451

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call splitGridFractal(child, rho, aFac, grid, converged)
                exit
             end if
          end do
       else
!          if (thisOctal%rho(subcell)*cellVolume(thisOctal, subcell) > limitScalar) then
          split = .true.
          if (thisOctal%threeD.and.((nThreadsGlobal-1)==8)) then
             if (thisOctal%nDepth == 1) then
                split = .true.
             else if (thisOctal%nDepth > 1) then
                if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                   split = .false.
                endif
             endif
          endif

          if (thisOctal%threed.and.((nThreadsGlobal-1)==64)) then
             if (thisOctal%nDepth <= 2) then
                split = .true.
             else
                if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                   split = .false.
                endif
             endif
          endif

          if (thisOctal%threed.and.((nThreadsGlobal-1)==512)) then
             if (thisOctal%nDepth <= 3) then
                split = .true.
             else
                if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                   split = .false.
                endif
             endif
          endif


          if (inputnSource>0) then
             centre = subcellCentre(thisOctal, subcell)
             do iSource = 1, inputNsource
                r1 = modulus(centre - sourcePos(iSource))/smallestCellSize
                if ((r1 < 12.d0) .and. (thisOctal%nDepth < maxDepthAMR)) then
                   split = .true.
                endif
                if ((r1 < 24.d0) .and. (thisOctal%nDepth < maxDepthAMR-1)) then
                   split = .true.
                endif
             enddo
          endif

          if (inputnSource > 0) then
             do iSource = 1, inputnSource
                if (inSubcell(thisOctal,subcell, sourcePos(isource)) &
                     .and. (thisOctal%nDepth < maxDepthAMR)) then
                   split = .true.
                   exit
                endif
             enddo
          endif

!          if (thisOctal%nDepth == minDepthAMR) split = .false.
!          if (myrankGlobal==1) then
!             write(*,*) "depth ",thisOctal%nDepth, " subcell ",subcell, " mpithread ",thisOctal%mpithread(subcell), " split ",split
!          endif
          if (split) then
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=.true., interp=.false.)
             ! find the child
             do j = 1, thisOctal%nChildren, 1
                if (thisOctal%indexChild(j) == subcell) then
                   child => thisOctal%child(j)
                   exit
                endif
             enddo

             allocate(r(1:thisOctal%maxChildren), s(1:thisOctal%maxChildren))
             call randomNumberGenerator(getRealArray=r)
             tot = sum(r)
             mean = tot / real(thisOctal%maxChildren)
             r = r / mean
             rmin = minval(r)
             rmax = maxval(r)
             fac = (afac*rmin-rmax)/(1.-afac)
             s = r + fac
             tot=SUM(s)
             s = s / tot
             do j = 1, thisOctal%maxChildren
                child%rho(j) = s(j) * thisOctal%rho(subcell) * &
                     cellVolume(thisOctal, subcell)/cellVolume(child,j)

                child%dustTypeFraction(j,:) = 1.d-20
                child%temperature(j) = 10.d0
                child%velocity(j) = VECTOR(0.d0, 0.d0, 0.d0)
                !Thaw - will probably want to change this to use returnMu
                
                ethermal = (1.d0/(mHydrogen))*kerg*child%temperature(j)

                if (hydrodynamics) then
                   child%pressure_i(j) = child%rho(j)*ethermal
                   child%energy(j) = ethermal + 0.5d0*(cspeed*modulus(child%velocity(j)))**2
                   child%rhoe(j) = child%rho(j) * child%energy(j)
                   child%gamma(j) = 1.0
                   child%iEquationOfState(j) = 1
                endif

                child%inFlow(j) = .true.


                if (photoionPhysics) then
                   child%nh(j) = child%rho(j) / mHydrogen
                   child%ne(j) = child%nh(j)
                   child%nhi(j) = 1.e-5
                   child%nhii(j) = child%ne(j)
                   child%nHeI(j) = 0.d0 !0.1d0 *  child%nH(j)    
                   child%ionFrac(j,1) = 1.               !HI
                   child%ionFrac(j,2) = 1.e-10           !HII
                endif



             enddo
             deallocate(r, s)
             converged = .false.
             exit
          endif
          
       endif
    end do
  end subroutine splitGridFractal



  subroutine postSetupChecks(grid)
#ifdef SPH
    use sph_data_class, only: sph_mass_within_grid, info_sph
#endif
    use inputs_mod, only : mDisc, geometry, sphWithChem
    use memory_mod, only : findTotalMemory, reportMemory
    type(GRIDTYPE) :: grid
    integer(kind=bigInt) :: i
#ifdef SPH
    character(len=80) :: message
    real(double) :: minRho, maxRho, totalmasstrap, totalmass, totalMassMol
#endif

    call findTotalMemory(grid, i)
    call reportMemory(i)


    select case (geometry)
    case ("shakara")
       call testAMRmass(grid, dble(mdisc))

#ifdef SPH
    case("molcluster", "theGalaxy", "cluster","sphfile", "dale")
       totalmasstrap = 0.0; maxrho=0.0; minrho=1.0e30; totalmass=0.0; totalmassMol=0.0
       call findTotalMass(grid%octreeRoot, totalMass, totalmasstrap = totalmasstrap, totalmassMol=totalmassMol, &
            maxrho=maxrho, minrho=minrho)
! This write statement is checked by the test suite so update checkSphToGrid.pl if it changes.
       write(message,*) "Mass of envelope: ",totalMass/mSol, " solar masses"
       call writeInfo(message, TRIVIAL)
       write(message,*) "Mass of envelope (TRAP): ",totalMasstrap/mSol, " solar masses"
       call writeInfo(message, TRIVIAL)
       if (sphwithchem) then 
          write(message,*) "Molecular mass of envelope: ",totalMassMol/mSol, " solar masses"
          call writeInfo(message, TRIVIAL)
       endif
       write(message,*) "Mass of SPH particles within grid: ", sph_mass_within_grid(grid), " solar masses"
       call writeInfo(message, TRIVIAL)
       write(message,*) "Maximum Density: ",maxrho, " g/cm^3"
       call writeInfo(message, TRIVIAL)
       write(message,*) "Minimum Density: ",minrho, " g/cm^3"
       call writeInfo(message, TRIVIAL)

! Write a VTK file so we can check the SPH to grid conversion
       if ( sphWithChem ) then 
          call writeVTKfile(grid, "gridFromSph.vtk", valueTypeString=(/"rho         ",&
         "temperature ", "velocity    ", "molabundance", "nh2         "/))
       else
          call writeVTKfile(grid, "gridFromSph.vtk")
       end if

! Write SPH information (e.g. code units) to a file
       call info_sph("info_sph.dat")

! N.B. This subroutine changes the geometry parameters
!       call fitAlphaDisc(grid)

#endif

    case("fitsfile")
       call writeVTKfile(grid, "gridFromFitsFile.vtk",  &
            valueTypeString=(/"rho        ", "temperature", "dust1      ","velocity   "/))

    case DEFAULT
    end select
  end subroutine postSetupChecks


  subroutine testAMRmass(grid, massWanted)
    type(GRIDTYPE) :: grid
    real(double) :: massWanted, actualMass, fac
    character(len=80) :: message
    
    actualMass = 0.d0
    call findTotalMass(grid%octreeRoot, actualMass)

    fac = abs(actualMass-massWanted)/massWanted
    write(message,'(a,f9.5,a)') "Grid mass is ",actualMass/msol, " solar masses"
    call writeINFO(message,TRIVIAL)

    if (fac > 0.01d0) then
       write(message,'(a,f7.1,a)') "Grid mass differs from required mass by: ",100.d0*fac, " %"
       call writeWarning(message)
       write(message,'(a,f8.5,a)') "Mass requested is: ",massWanted/msol, " solar masses"
       call writeWarning(message)
       write(message,'(a,f8.5,a)') "Grid mass is: ",actualMass/mSol, " solar masses"
       call writeWarning(message)
    endif
  end subroutine testAMRmass

  recursive subroutine myFreeGrid(thisOctal)
    type(octal), pointer :: thisOctal, child
    integer :: i

    call deallocateOctalDynamicAttributes(thisOctal)
    if (thisOctal%nChildren > 0) then
       do i = 1 , thisOctal%nChildren
          child => thisOctal%child(i)
          call myFreeGrid(child)
       enddo
       deallocate(thisOctal%child)
    endif

       
  end subroutine myFreeGrid

  recursive subroutine bigArraytest(thisOctal)
    type(octal), pointer :: thisOctal, child
    integer :: i

    deallocate(thisOctal%distanceGrid)
    allocate(thisOctal%distanceGrid(1:100000))

    if (thisOctal%nChildren > 0) then
       do i = 1 , thisOctal%nChildren
          child => thisOctal%child(i)
          call bigarraytest(child)
       enddo
    endif
  end subroutine bigarraytest


  recursive subroutine testVelocity(thisOctal, grid)
    type(octal), pointer :: thisOctal, child
    type(GRIDTYPE) :: grid
    integer :: subcell , i, j, k
    real(double) :: t1, t2, t3, r1, r2, phi1, phi2, inc
    real(Double) :: z1, z2,z3,r3,phi3
    type(VECTOR) :: centre, corner(8)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call testVelocity(child,grid )
                exit
             end if
          end do
       else

          thisOctal%etaline(subcell) = 0.d0
          if (thisOctal%inflow(subcell)) then
             thisOctal%chiLine(subcell) =  0.d0
             do j = 1, 27
                thisOCtal%chiline(subcell) = thisOctal%chiline(subcell) + modulus(thisOctal%cornerVelocity(j))
             enddo
             thisOctal%etaLine(subcell) =  VECTOR(1.d0, 0.d0, 0.d0).dot.&
                  amrGridVelocity(grid%octreeRoot,subcellCentre(thisOctal,subcell), foundOctal=thisOctal,&
                  actualSubcell=subcell)
             if (thisOctal%etaLine(subcell) > 1.d30) then
                write(*,*) "Vel bug ",thisOctal%etaline(subcell),thisOctal%splitAzimuthally, subcell

               centre = subcellCentre(thisOctal,subcell)
               r1 = sqrt(centre%x**2 + centre%y**2)
               r2 = sqrt(centre%x**2 + centre%y**2)
               phi1 = atan2(centre%y, centre%x)
               if (phi1 < 0.d0) phi1 = phi1 + twoPi
               phi2 = atan2(centre%y, centre%x)
               if (phi2 < 0.d0) phi2 = phi2 + twoPi

               inc = thisOctal%subcellSize / 2.0               

               t1 = MAX(0.0_oc, (r1 - (r2 - inc)) / thisOctal%subcellSize)
               t2 = (phi1 - (phi2 - thisOctal%dPhi/4.d0))/(thisOctal%dPhi/2.d0)
               t3 = MAX(0.0_oc, (centre%z - (centre%z - inc)) / thisOctal%subcellSize)
               write(*,*) "t1, t2, t3 ",t1,t2,t3
                do k = 1, 8
                   write(*,*) "velocities " ,k,(cspeed/1.d5)* thisOctal%velocity(k)
                enddo


             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phiMin
             phi2 = thisOctal%phi 
             phi3 = thisOctal%phiMax
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize

             write(*,*) "phi1 phi2 phi3 ",phi1*radtodeg, phi2*radtodeg, phi3*radtodeg

               select case(subcell)
            CASE(1)
               write(*,*) thisOctal%cornerVelocity( 1)
               write(*,*) "corner ",vector(r1*cos(phi1),r1*sin(phi1),z1), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi1),r1*sin(phi1),z1)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity( 2)
               write(*,*) "corner ",vector(r1*cos(phi2),r1*sin(phi2),z1), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi2),r1*sin(phi2),z1)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity( 4)
               write(*,*) "corner ",vector(r2*cos(phi1),r2*sin(phi1),z1), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi1),r2*sin(phi1),z1)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity( 5)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z1), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z1)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(10)
               write(*,*) "corner ",vector(r1*cos(phi1),r1*sin(phi1),z2), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi1),r1*sin(phi1),z2)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(11)
               write(*,*) "corner ",vector(r1*cos(phi2),r1*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi2),r1*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(13)
               write(*,*) "corner ",vector(r2*cos(phi1),r2*sin(phi1),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi1),r2*sin(phi1),z2)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(14)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z2)), phi2*radtodeg
               
            CASE(2)
               write(*,*) thisOctal%cornerVelocity( 2)
               write(*,*) "corner ",vector(r1*cos(phi2),r1*sin(phi2),z1), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi2),r1*sin(phi2),z1)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity( 3)
               write(*,*) "corner ",vector(r1*cos(phi3),r1*sin(phi3),z1), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi3),r1*sin(phi3),z1)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity( 5)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z1), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z1)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity( 6)
               write(*,*) "corner ",vector(r2*cos(phi3),r2*sin(phi3),z1), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi3),r2*sin(phi3),z1)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(11)
               write(*,*) "corner ",vector(r1*cos(phi2),r1*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi2),r1*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(12)
               write(*,*) "corner ",vector(r1*cos(phi3),r1*sin(phi3),z2), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi3),r1*sin(phi3),z2)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(14)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(15)
               write(*,*) "corner ",vector(r2*cos(phi3),r2*sin(phi3),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi3),r2*sin(phi3),z2)), phi3*radtodeg
               
            CASE(3)
               write(*,*) thisOctal%cornerVelocity( 4)
               write(*,*) "corner ",vector(r2*cos(phi1),r2*sin(phi1),z1), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi1),r2*sin(phi1),z1)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity( 5)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z1), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z1)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity( 7)
               write(*,*) "corner ",vector(r3*cos(phi1),r3*sin(phi1),z1), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi1),r3*sin(phi1),z1)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity( 8)
               write(*,*) "corner ",vector(r3*cos(phi2),r3*sin(phi2),z1), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi2),r3*sin(phi2),z1)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(13)
               write(*,*) "corner ",vector(r2*cos(phi1),r2*sin(phi1),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi1),r2*sin(phi1),z2)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(14)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(16)
               write(*,*) "corner ",vector(r3*cos(phi1),r3*sin(phi1),z2), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi1),r3*sin(phi1),z2)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(17)
               write(*,*) "corner ",vector(r3*cos(phi2),r3*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi2),r3*sin(phi2),z2)), phi2*radtodeg
               
            CASE(4)
               write(*,*) thisOctal%cornerVelocity( 5)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z1), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z1)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity( 6)
               write(*,*) "corner ",vector(r2*cos(phi3),r2*sin(phi3),z1), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi3),r2*sin(phi3),z1)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity( 8)
               write(*,*) "corner ",vector(r3*cos(phi2),r3*sin(phi2),z1), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi2),r3*sin(phi2),z1)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity( 9)
               write(*,*) "corner ",vector(r3*cos(phi3),r3*sin(phi3),z1), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi3),r3*sin(phi3),z1)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(14)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(15)
               write(*,*) "corner ",vector(r2*cos(phi3),r2*sin(phi3),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi3),r2*sin(phi3),z2)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(17)
               write(*,*) "corner ",vector(r3*cos(phi2),r3*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi2),r3*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(18)
               write(*,*) "corner ",vector(r3*cos(phi3),r3*sin(phi3),z2), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi3),r3*sin(phi3),z2)), phi3*radtodeg
               
            CASE(5)
               write(*,*) thisOctal%cornerVelocity(10)
               write(*,*) "corner ",vector(r1*cos(phi1),r1*sin(phi1),z2), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi1),r1*sin(phi1),z2)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(11)
               write(*,*) "corner ",vector(r1*cos(phi2),r1*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi2),r1*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(13)
               write(*,*) "corner ",vector(r2*cos(phi1),r2*sin(phi1),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi1),r2*sin(phi1),z2)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(14)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(19)
               write(*,*) "corner ",vector(r1*cos(phi1),r1*sin(phi1),z3), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi1),r1*sin(phi1),z3)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(20)
               write(*,*) "corner ",vector(r1*cos(phi2),r1*sin(phi2),z3), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi2),r1*sin(phi2),z3)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(22)
               write(*,*) "corner ",vector(r2*cos(phi1),r2*sin(phi1),z3), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi1),r2*sin(phi1),z3)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(23)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z3), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z3)), phi2*radtodeg
               
            CASE(6)
               write(*,*) thisOctal%cornerVelocity(11)
               write(*,*) "corner ",vector(r1*cos(phi2),r1*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi2),r1*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(12)
               write(*,*) "corner ",vector(r1*cos(phi3),r1*sin(phi3),z2), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi3),r1*sin(phi3),z2)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(14)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(15)
               write(*,*) "corner ",vector(r2*cos(phi3),r2*sin(phi3),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi3),r2*sin(phi3),z2)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(20)
               write(*,*) "corner ",vector(r1*cos(phi2),r1*sin(phi2),z3), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi2),r1*sin(phi2),z3)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(21)
               write(*,*) "corner ",vector(r1*cos(phi3),r1*sin(phi3),z3), &
                    inflowMahdavi(1.d10*vector(r1*cos(phi3),r1*sin(phi3),z3)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(23)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z3), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z3)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(24)
               write(*,*) "corner ",vector(r2*cos(phi3),r2*sin(phi3),z3), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi3),r2*sin(phi3),z3)), phi3*radtodeg
            
            CASE(7)
               write(*,*) thisOctal%cornerVelocity(13)
               write(*,*) "corner ",vector(r2*cos(phi1),r2*sin(phi1),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi1),r2*sin(phi1),z2)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(14)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(16)
               write(*,*) "corner ",vector(r3*cos(phi1),r3*sin(phi1),z2), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi1),r3*sin(phi1),z2)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(17)
               write(*,*) "corner ",vector(r3*cos(phi2),r3*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi2),r3*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(22)
               write(*,*) "corner ",vector(r2*cos(phi1),r2*sin(phi1),z3), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi1),r2*sin(phi1),z3)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(23)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z3), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z3)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(25)
               write(*,*) "corner ",vector(r3*cos(phi1),r3*sin(phi1),z3), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi1),r3*sin(phi1),z3)), phi1*radtodeg
               write(*,*) thisOctal%cornerVelocity(26)
               write(*,*) "corner ",vector(r3*cos(phi2),r3*sin(phi2),z3), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi2),r3*sin(phi2),z3)), phi2*radtodeg
               
            CASE(8)
               write(*,*) thisOctal%cornerVelocity(14)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(15)
               write(*,*) "corner ",vector(r2*cos(phi3),r2*sin(phi3),z2), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi3),r2*sin(phi3),z2)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(17)
               write(*,*) "corner ",vector(r3*cos(phi2),r3*sin(phi2),z2), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi2),r3*sin(phi2),z2)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(18)
               write(*,*) "corner ",vector(r3*cos(phi3),r3*sin(phi3),z2), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi3),r3*sin(phi3),z2)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(23)
               write(*,*) "corner ",vector(r2*cos(phi2),r2*sin(phi2),z3), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi2),r2*sin(phi2),z3)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(24)
               write(*,*) "corner ",vector(r2*cos(phi3),r2*sin(phi3),z3), &
                    inflowMahdavi(1.d10*vector(r2*cos(phi3),r2*sin(phi3),z3)), phi3*radtodeg
               write(*,*) thisOctal%cornerVelocity(26)
               write(*,*) "corner ",vector(r3*cos(phi2),r3*sin(phi2),z3), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi2),r3*sin(phi2),z3)), phi2*radtodeg
               write(*,*) thisOctal%cornerVelocity(27)
               write(*,*) "corner ",vector(r3*cos(phi3),r3*sin(phi3),z3), &
                    inflowMahdavi(1.d10*vector(r3*cos(phi3),r3*sin(phi3),z3)), phi3*radtodeg

            CASE DEFAULT
               PRINT *, 'Invalid subcell in amrGridVelocity'
               end select

             call subcellCorners(thisOctal, subcell, corner)
             write(*,*) "centre ",subcellCentre(thisOctal,subcell), &
                  inflowMahdavi(1.d10*subcellCentre(thisOCtal,subcell))
             do j = 1, 8
                write(*,*) "subcell corner ",j,corner(j)
                write(*,*) "corner test ",j,inflowMahdavi(1.d10*corner(j))
                phi1 = atan2(corner(j)%y,corner(j)%x)*radtodeg
                if (phi1 < 0.d0) phi1 = phi1 + 360.d0
                write(*,*) "phi ", phi1
             enddo



             endif






          endif


       endif
    enddo
  end subroutine testVelocity



  recursive subroutine set_bias_cmfgen(thisOctal, grid, lambda0)
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real, intent(in)          :: lambda0                ! rest wavelength of line
  integer :: subcell, i
  real(double) :: d, dV, r, tauSob, escProb
  type(vector)  :: rvec, rhat, direction
  real(double):: nu0, dr, dA, tau
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call set_bias_cmfgen(child, grid, lambda0)
                exit
             end if
          end do

       else
          if (thisOctal%inflow(subcell)) then
             d = thisOctal%subcellsize 

             rVec = subcellCentre(thisOctal,subcell)
             r = modulus(rvec)
             if (r /= 0.0d0) then
                rhat = rvec/r
             else
                rhat = VECTOR(0.0d0, 0.0d0, 1.0d0)
             end if

             if (thisOctal%threed) then
                dV = d*d*d
             else
                dr = d
                dA = 2.0_db*pi*(SQRT(rVec%x*rVec%x  + rVec%y*rVec%y))*dr
                dV = d*dA
             endif

             nu0  = cSpeed_dbl / dble(lambda0*angstromtocm)
              ! in a radial direction
             tauSob = thisOctal%chiline(subcell)  / nu0
             tauSob = tauSob / amrGridDirectionalDeriv(grid, rvec, rhat, &
                  startOctal=thisOctal)
           
             if (tauSob < 0.01) then
                escProb = 1.0d0-tauSob*0.5d0*(1.0d0 -   &
                     tauSob/3.0d0*(1. - tauSob*0.25d0*(1.0d0 - 0.20d0*tauSob)))
             else if (tauSob < 15.) then
                escProb = (1.0d0-exp(-tauSob))/tauSob
             else
                escProb = 1.d0/tauSob
             end if
             escProb = max(escProb, 1.d-10)

             rVec = subcellCentre(thisoctal, subcell)
             direction = rVec
             call normalize(direction)
             call tauAlongPathCMF(grid, rVec, direction, tau, tauMax=1.d3, &
                  startOctal=thisOctal, startSubcell=subcell)
             thisOctal%biasCont3d(subcell) = max(exp(-tau),1.d-6)
             thisOctal%biasLine3D(subcell) = max(escProb*thisOctal%biasCont3D(subcell),1.d-6)


          else  ! this subcell is not "inFlow"
             thisOctal%biasCont3D(subcell) = 1.0d-150
             thisOctal%biasLine3D(subcell) = 1.0d-150
          end if
          ! just in case ....
          thisOctal%biasCont3D(subcell) = MAX(thisOctal%biasCont3D(subcell), 1.0d-150)
          thisOctal%biasLine3D(subcell) = MAX(thisOctal%biasLine3D(subcell), 1.0d-150)
       endif ! if (thisOctal%hasChild(subcell)) then
    enddo

  end subroutine set_bias_cmfgen

  subroutine turbulentVelocityField(grid, vDispersion)
    use inputs_mod, only : maxDepthAMR
    use utils_mod, only: gasdev
    type(GRIDTYPE) :: grid
    real(double) :: vDispersion
    real(double) :: deltaV, wavenumber, vmag
    type(VECTOR) :: vel
    real(double) :: xmax, xmin, ymax, ymin, zmax, zmin, sigma, fac
    integer :: i, j, k, idepth, iAcross

    do iDepth = 1, maxDepthAmr
       waveNumber = twoPi * dble(2**(iDepth-1))
       deltaV  = wavenumber**(-4.d0)
       vMag = deltaV * gasdev()
       iAcross = 2**(iDepth)
       do i = 1, iAcross
          do j = 1, iAcross
             do k = 1, iAcross

                vel = vMag * randomUnitVector()
                xmin = 2.d0 * grid%octreeRoot%subcellSize * dble(i-1)/dble(iAcross) - grid%octreeRoot%subcellSize
                xmax = 2.d0 * grid%octreeRoot%subcellSize * dble(i)/dble(iAcross) - grid%octreeRoot%subcellSize
                
                ymin = 2.d0 * grid%octreeRoot%subcellSize * dble(j-1)/dble(iAcross) - grid%octreeRoot%subcellSize
                ymax = 2.d0 * grid%octreeRoot%subcellSize * dble(j)/dble(iAcross) - grid%octreeRoot%subcellSize
                
                zmin = 2.d0 * grid%octreeRoot%subcellSize * dble(k-1)/dble(iAcross) - grid%octreeRoot%subcellSize
                zmax = 2.d0 * grid%octreeRoot%subcellSize * dble(k)/dble(iAcross) - grid%octreeRoot%subcellSize
                call applyVOverRange(grid%octreeRoot, vel, xmin, xmax, ymin, ymax, zmin, zmax)
             enddo
          enddo
       enddo
    enddo
    vel = VECTOR(0.d0,0.d0,0.d0)
    i = 0
    call sumVelocity(grid%octreeRoot, vel, i)
    vel = vel / dble(i)
    call subtractVelocity(grid%octreeRoot,vel)
    sigma = 0.d0
    call findSigma(grid%octreeRoot, sigma)
    sigma = sqrt(sigma / dble(i))
    fac = vDispersion/sigma
    call scaleVelocity(grid%octreeRoot, fac)
  end subroutine turbulentVelocityField

  recursive subroutine applyVoverRange(thisOctal,  vel, xmin, xmax, ymin, ymax, zmin, zmax)
    type(octal), pointer :: thisOctal, child
    type(VECTOR) :: vel,cvec
    real(double) :: xmin, xmax, ymin, ymax, zmin, zmax
    integer :: subcell , i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call applyVoverRange(child,  vel, xmin, xmax, ymin, ymax, zmin, zmax)
                exit
             end if
          end do
       else
          cVec = subcellCentre(thisOctal, subcell)
          if ( (cVec%x >= xmin).and.(cVec%x <= xmax).and. &
               (cVec%y >= ymin).and.(cVec%y <= ymax).and. &
               (cVec%z >= zmin).and.(cVec%z <= zmax)) then
             thisOctal%velocity(subcell) = thisOctal%velocity(subcell) + vel
          endif
       endif
    enddo
  end subroutine applyVoverRange



  recursive subroutine sumVelocity(thisOctal,vel,n)
    type(octal), pointer :: thisOctal, child
    integer :: n
    type(VECTOR) :: vel
    integer :: subcell , i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sumVelocity(child, vel,n)
                exit
             end if
          end do
       else
          n = n + 1
          vel = vel + thisOctal%velocity(subcell)
       endif
    enddo
  end subroutine sumVelocity

  recursive subroutine subtractVelocity(thisOctal,vel)
    type(octal), pointer :: thisOctal, child
    type(VECTOR) :: vel
    integer :: subcell , i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call subtractVelocity(child, vel)
                exit
             end if
          end do
       else
          thisOctal%velocity(subcell) = thisOctal%velocity(subcell) - vel
       endif
    enddo
  end subroutine subtractVelocity

  recursive subroutine scaleVelocity(thisOctal, fac)
    type(octal), pointer :: thisOctal, child
    real(double) :: fac
    integer :: subcell , i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call scaleVelocity(child, fac)
                exit
             end if
          end do
       else
          thisOctal%velocity(subcell) = thisOctal%velocity(subcell) * fac
       endif
    enddo
  end subroutine scaleVelocity

  recursive subroutine findSigma(thisOctal, sigma)
    type(octal), pointer :: thisOctal, child
    real(double) :: sigma
    integer :: subcell , i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findSigma(child, sigma)
                exit
             end if
          end do
       else
          sigma = sigma + thisOctal%velocity(subcell)%x**2
       endif
    enddo
  end subroutine findSigma

recursive subroutine quickSublimate(thisOctal)
  use inputs_mod, only: grainFrac, nDustType

  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call quickSublimate(child)
                exit
             end if
          end do
       else

          if (associated(thisOctal%dustTypeFraction)) then 
             if (thisOctal%temperature(subcell) > 1500.) then
                thisOctal%dustTypeFraction(subcell,:) = 1.d-20
             else
                thisOctal%dustTypeFraction(subcell,nDustType) = grainFrac(nDustType)
             endif
          end if

    endif
    enddo
  end subroutine quickSublimate


  subroutine SanityCheckGrid(grid)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, boundVec
    type(OCTAL), pointer :: thisOctal, testOctal
    integer :: subcell, i
    if (myrankWorldGlobal == 1) then
       do i = 1, 100
          write(*,*) "Running AMR sanity check ",i

10        continue
          rVec = randomUnitVector()
          rVec = grid%OctreeRoot%centre + rVec * grid%octreeRoot%subcellSize

          testOctal => grid%octreeRoot
          thisOctal => grid%octreeRoot

          call findSubcellLocal(rVec, testOctal, subcell)
          if (.not.octalOnThread(thisOctal, subcell, myrankGlobal)) goto 10
          write(*,*) "testOctal found OK"

          write(*,*) "grid centre ",grid%octreeRoot%centre
          write(*,*) "grid x min.max ",grid%octreeRoot%xMin, grid%octreeRoot%xMax
          write(*,*) "grid z min.max ",grid%octreeRoot%zMin, grid%octreeRoot%zMax


          boundVec = VECTOR(testOctal%xMin, 0.d0, rVec%z)
          call findSubcellTD(boundVec, grid%octreeRoot, thisOctal, subcell)
          write(*,*) "xMin test passed"

          call findSubcellLocal(boundVec, thisOctal, subcell)
          write(*,*) "xMin test passed with findsubcelllocal"


          boundvec = VECTOR(testOctal%xMax, 0.d0, rVec%z)
          call findSubcellTD(boundVec, grid%octreeRoot, thisOctal, subcell)
          write(*,*) "xMax test passed"

          call findSubcellLocal(boundVec, thisOctal, subcell)
          write(*,*) "xMax test passed with findsubcelllocal"

          boundVec = VECTOR(rVec%x, 0.d0, testOctal%zMin)
          call findSubcellTD(boundVec, grid%octreeRoot, thisOctal, subcell)
          write(*,*) "zMin test passed"

          call findSubcellLocal(boundVec, thisOctal, subcell)
          write(*,*) "zMin test passed with findsubcelllocal"

          boundVec = VECTOR(rVec%x, 0.d0, testOctal%zMax)
          call findSubcellTD(boundVec, grid%octreeRoot, thisOctal, subcell)
          write(*,*) "zMax test passed"

          call findSubcellLocal(boundVec, thisOctal, subcell)
          write(*,*) "zMax test passed with findsubcelllocal"

       enddo
    endif
    call torus_mpi_barrier
  end subroutine SanityCheckGrid

!!$subroutine fitAlphaDisc(grid)
!!$  use inputs_mod, only : alphaDisc, betaDisc, height, rinner, router, rho0, geometry, heightsplitfac
!!$  type(GRIDTYPE) :: grid
!!$  real(double) :: alpha, beta, heightDouble !, rhoDouble
!!$  real(double) :: rInnerDouble, rOuterDouble!, chisq, minChisq
!!$!  integer :: n, i1,i2,i3,i4
!!$  
!!$  rInnerDouble = 750d0
!!$  rOuterDouble = 23d3
!!$
!!$  alpha = 0.5
!!$  beta = 1.125
!!$  rho0 = 1.e-7
!!$  heightDouble = 10.d0*autocm/1.d10
!!$
!!$!  minChisq = 1.d30
!!$!  do i1 = 1, 10
!!$!     do i2 = 1, 10
!!$!        do i3 = 1,10
!!$!           do i4 = 1, 10
!!$!              alpha = 2.d0*dble(i1-1)/10.d0
!!$!              beta = 2.d0*dble(i2-1)/10.d0
!!$!              heightDouble = (1.d0+dble(i3-1)/10.d0 * 9.d0)*autocm/1.d10
!!$!              rhodouble = -11.d0 + 2.d0 * dble(i4-1)/10.d0
!!$!              rhodouble = 10.d0**rhodouble
!!$!
!!$!              chisq = 0.
!!$!              n = 0
!!$!              call chisqAlphaDisc(grid%octreeRoot,  alpha, beta, rhodouble, heightDouble, rInnerDouble, rOuterDouble, chisq, n)
!!$!              chisq = chisq / dble(n-4)
!!$!              if (chisq < minChisq) then
!!$!                 if (writeoutput) write(*,'(1p,5e12.3)') alpha, beta, heightDouble, rho0, chisq
!!$!                 minChisq = chisq
!!$!              endif
!!$!           enddo
!!$!        enddo
!!$!     enddo
!!$!  enddo
!!$  rinner = 10. * real(rSol   / 1.d10 )
!!$  router = 100.* real(autocm / 1.d10 )
!!$  height = 5.*   real(autocm / 1.d10 )
!!$  alphaDisc = 0.5
!!$  betaDisc = 1.
!!$  rho0 = 1.e-10
!!$  heightsplitFac = 1.
!!$  grid%geometry = "adddisc"
!!$  geometry = "adddisc"
!!$  call addDisc(grid%octreeRoot, grid)
!!$
!!$
!!$end subroutine fitAlphaDisc


recursive subroutine chisqAlphaDisc(thisOctal, alpha, beta, rho0, height, rInner, rOuter, chisq, n)
  real(double) :: alpha, beta, rho0, rinner, router, chisq, r, h, height, fac, rho
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  type(VECTOR) :: rVec
  integer :: subcell, i, n
  
  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call chisqAlphaDisc(child, alpha, beta, rho0, height, rInner, rOuter, chisq, n)
              exit
           end if
        end do
     else 
        rVec = subcellCentre(thisOctal, subcell)
        r = sqrt(rVec%x**2 + rVec%y**2)
        if ((r > rInner).and.(r < rOuter).and.(thisOctal%rho(subcell) < 1.d-10).and.&
             (thisOctal%rho(subcell) > 1.d-15)) then
           h = height * (r / (100.d0*autocm/1.d10))**beta
           fac = -0.5d0 * (rVec%z/h)**2
           fac = max(-50.d0,fac)
           rho = dble(rho0) * (dble(rInner/r))**dble(alpha) * exp(fac)
           chisq = chisq + (thisOctal%rho(subcell)-rho)**2 / rho**2
           n = n + 1
        endif
     endif
  enddo
end subroutine chisqAlphaDisc


subroutine addSpiralWake(grid)

  use inputs_mod, only : sourcePos, rInner, rOuter, rho0, betaDisc, alphaDisc, height, sourceMass

  type(GRIDTYPE) :: grid
  real(double) :: r, r0, epsilon, phi, rSpiral, hillRadius
  type(VECTOR) :: spiralVec, thisVec
  integer :: i, k, m
  integer :: nr, nphi, ndr, nz, subcell
  real(double) :: dphi, dr, rwidth, phiwidth, fac, rhoSpiral, h, z, largescalefac
  type(OCTAL), pointer :: thisOctal
  thisOctal => grid%octreeRoot
  r0 = modulus(sourcePos(2))
  nr = 100000
  nphi = 10
  ndr = 100
  nz = 201
  epsilon = 0.1
  hillRadius = modulus(sourcePos(2)-sourcePos(1))*sqrt(sourceMass(2)/(3.d0*sourceMass(1)))
  if (Writeoutput) write(*,*) "Adding spiral wake..."
  do i = 1, nr
     if (Writeoutput) write(*,*) i,nr

     r  = rinner + (rOuter - rInner)*dble(i-1)/dble(nr-1)

     if (abs(r-r0) > 0.3*hillRadius) then
     
     r = r/r0
     if (r > 1.d0) then
        phi = (2.d0/(3.d0*epsilon)) * ( r**1.5d0 - 1.5d0*log(r) - 1.d0)
     else
        phi = -(2.d0/(3.d0*epsilon)) * ( r**1.5d0 - 1.5d0*log(r) - 1.d0)
     endif
     phi = phi + pi
     rSpiral = r * r0
     spiralVec = VECTOR(rSpiral*cos(phi), rSpiral*sin(phi), 0.d0)

     phi = atan2(spiralVec%y, spiralVec%x)
     dr = abs(r-rSpiral)

     rwidth = r * fourpi * epsilon**2
     phiWidth = fourPi * epsilon

     h = height * (rspiral / (100.d0*autocm/1.d10))**betaDisc
     do m = 1, nz
        z = -5.d0*h  + 10.d0*h*dble(m-1)/dble(nz-1)
        fac = -0.5d0 * (dble(z)/h)**2
        fac = max(-50.d0,fac)
        rhoSpiral = dble(rho0) * (dble(rInner/rSpiral))**dble(alphaDisc) * exp(fac)

        do k = 1, ndr
           dphi = 0.d0
           dr = -rWidth/2.d0 + dble(k-1)/dble(ndr-1) * rWidth
           thisVec = VECTOR((rSpiral+dr)*cos(phi+dphi),(rSpiral+dr)*sin(phi+dphi),z)
           call findSubcellLocal(thisVec, thisOctal, subcell)

           fac = (abs(dr)/rwidth) + dphi/phiWidth
           largescaleFac = abs((rSpiral-r0)/(2.d0*r0))

           if ((dphi < phiWidth).and.(abs(dr) < rWidth)) then
              thisOctal%rho(subcell) = max(thisOctal%rho(subcell), rhoSpiral * 10.d0 * exp(-fac)*exp(-largescalefac))
           endif
        enddo
     enddo
  endif
  enddo



  if (Writeoutput) write(*,*) "Done."
end subroutine addSpiralWake




end module setupamr_mod
