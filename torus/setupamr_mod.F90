module setupamr_mod

  use cmfgen_class
  use stateq_mod
  use amr_mod
  use vtk_mod
  use vector_mod
  use magnetic_mod
  use messages_mod
  USE constants_mod
  USE octal_mod, only: OCTAL, wrapperArray, octalWrapper, subcellCentre, cellVolume, &
       allocateattribute, copyattribute, deallocateattribute
  use gridtype_mod, only:   gridtype
  USE parallel_mod, ONLY:   torus_abort
  use mpi_global_mod

  implicit none

contains
    

  subroutine setupamrgrid(grid)
    use cluster_class
    use gridio_mod
    use amr_mod
    use lucy_mod
    use grid_mod
    use photoionAMR_mod, only : resizePhotoionCoeff
    use input_variables, only : readgrid, gridinputfilename, geometry, mdot, constantAbundance
    use input_variables, only : amrGridCentreX, amrGridCentreY, amrGridCentreZ
    use input_variables, only : amr1d, amr2d, amr3d, splitOverMPI
    use input_variables, only : amrGridSize, doSmoothGrid, photoionPhysics, dustPhysics
    use input_variables, only : ttauriRstar, mDotparameter1, ttauriWind, ttauriDisc, ttauriWarp
    use input_variables, only : limitScalar, limitScalar2, smoothFactor, onekappa
    use input_variables, only : CMFGEN_rmin, CMFGEN_rmax, textFilename, sphDataFilename, inputFileFormat
    use input_variables, only : rCore, rInner, rOuter, lamline
    use disc_class, only:  new
    use discwind_class, only:  new
    use sph_data_class, only: new_read_sph_data, read_galaxy_sph_data
#ifdef MPI 
    use mpi_amr_mod
    use photoionAMR_mod, only : ionizeGrid, resetNh
#endif
    use vh1_mod, only: read_vh1
    use memory_mod

    implicit none

    ! For romanova geometry case
    type(romanova) :: romData ! parameters and data for romanova geometry
    type(cluster)   :: young_cluster
    type(VECTOR) :: amrGridCentre
    character(len=80) :: newContFluxFile
    real :: theta1, theta2
    logical :: ok, flatspec
    type(GRIDTYPE) :: grid
    logical :: gridConverged
    real(double) :: astar, mass_accretion_old, totalMass
    real(double) :: minRho, maxRho, totalmasstrap
    character(len=80) :: message
    integer :: nVoxels, nOctals
!    integer(bigInt) :: i
!    integer :: nUnrefine
    constantAbundance = .true.

#ifdef MPI
    call randomNumberGenerator(randomSeed=.true.)
    call randomNumberGenerator(syncIseed=.true.)
#endif

    call writeBanner("Setting up AMR grid","-",TRIVIAL)

    totalmass = 0.
    if ( dustPhysics ) then 
       flatspec = .false.
    else
       flatspec = .true.
    end if
    call new(young_cluster, .false.)

    if (readgrid) then
       grid%splitOverMPI = splitOverMPI
       call readAMRgrid(gridInputfilename, .false., grid)
       grid%splitOverMPI = splitOverMPI
       if (photoIonPhysics) call resizePhotoionCoeff(grid%octreeRoot, grid)

       call findTotalMemory(grid, globalMemoryFootprint)

    else

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

       call initAMRGrid(newContFluxFile,flatspec,grid,ok,theta1,theta2)
       grid%splitOverMPI = splitOverMPI

       amrGridCentre = VECTOR(amrGridCentreX,amrGridCentreY,amrGridCentreZ)
       call writeInfo("Starting initial set up of adaptive grid...", TRIVIAL)

       select case (geometry)

       case("fogel")
          constantAbundance = .false.
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call setupFogel(grid, textFilename, "HCN")
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)


       case("cluster")

          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

       case("theGalaxy")

          select case (inputFileFormat)
          case("binary")
             call read_galaxy_sph_data(sphdatafilename)

          case("ascii")
             call new_read_sph_data(sphdatafilename)

          case default
             call writeFatal("Unrecognised file format "//inputFileFormat)

          end select

          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, young_cluster)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

       case("molcluster")
          if (.not.grid%splitOverMPI) then
             call new_read_sph_data(sphdatafilename)
          else
#ifdef MPI
             if (myrankGlobal == 0) then
                call new_read_sph_data(sphdatafilename)
             endif
             call distributeSphDataOverMPI()
#endif
          endif
          call writeInfo("Initialising adaptive grid...", TRIVIAL)
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid, young_cluster)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

       case("wr104")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
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

       case("clumpyagb")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid)
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
         do while(.not.gridconverged) 
             call splitGridFractal(grid%octreeRoot, real(100.*mHydrogen), 0.1, grid, gridconverged)
          enddo
          call ionizeGrid(grid%octreeRoot)
          call resetNH(grid%octreeRoot)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)
          call findMassOverAllThreads(grid, totalmass)
          write(message,'(a,1pe12.5,a)') "Total mass in fractal cloud (solar masses): ",totalMass/lsol
          call writeInfo(message,TRIVIAL)
#endif

       case("runaway")
          call read_vh1
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d,  romData=romData) 
          call writeInfo("First octal initialized.", TRIVIAL)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,romData=romData)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

          if (doSmoothGrid) then
             call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
             do
                gridConverged = .true.
                ! The following is Tim's replacement for soomthAMRgrid.
                call myScaleSmooth(smoothfactor, grid, &
                     gridConverged,  inheritProps = .false., &
                     interpProps = .false.,  &
                     stellar_cluster=young_cluster, romData=romData)
                if (gridConverged) exit
             end do
             call writeInfo("...grid smoothing complete", TRIVIAL)
          endif

       case("ttauri")
          call initFirstOctal(grid,amrGridCentre,amrGridSize, amr1d, amr2d, amr3d, romData=romData) 
          call writeInfo("First octal initialized.", TRIVIAL)
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,romData=romData)
          call writeInfo("Smoothing adaptive grid structure...", TRIVIAL)
          do
             gridConverged = .true.
             call myScaleSmooth(3., grid, &
                  gridConverged,  inheritProps = .false., interpProps = .false.)
             if (gridConverged) exit
          end do
          call writeInfo("...grid smoothing complete", TRIVIAL)
          
          call fixParentPointers(grid%octreeRoot)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

           ! This section is getting rather long. Maybe this should be done in 
           ! wrapper subroutine in amr_mod.f90.
          mdot = 2.d-8 * msol * secstoyears
          call zeroDensity(grid%octreeRoot)
          astar = accretingAreaMahdavi(grid)
          ttauriwind = .false.
          ttauridisc = .false.
          if (writeoutput) write(*,*) "accreting area (%) ",100.*astar/(fourpi*ttauriRstar**2)
          call assignDensitiesMahdavi(grid, grid%octreeRoot, astar, mDotparameter1*mSol/(365.25d0*24.d0*3600.d0))
          if (ttauriwind) call assignDensitiesBlandfordPayne(grid, grid%octreeRoot)
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
          call splitGrid(grid%octreeRoot,limitScalar,limitScalar2,grid,romData=romData)
          call fixParentPointers(grid%octreeRoot)
          call writeInfo("...initial adaptive grid configuration complete", TRIVIAL)

          ! This section is getting rather long. Maybe this should be done in 
          ! wrapper subroutine in amr_mod.f90.
          if (geometry=="ttauri") then 
             mdot = 2.d-8 * msol * secstoyears
             !             do i = 1, 2
             !                call zeroDensity(grid%octreeRoot)
             !                call assignDensitiesMahdavi(grid, dble(mdot))
             !                gridconverged = .false.
             !                do while (.not.gridconverged)
             !                   gridConverged = .true.
             !                   call massSplit(grid%octreeRoot, grid, gridconverged, inheritProps=.true., interpProps=.false.)
             !                   if (gridConverged) exit
             !                enddo
             !             enddo
             call zeroDensity(grid%octreeRoot)
             astar = accretingAreaMahdavi(grid)
             if (writeoutput) write(*,*) "accreting area (%) ",100.*astar/(fourpi*ttauriRstar**2)
             call assignDensitiesMahdavi(grid, grid%octreeRoot, astar, mDotparameter1*mSol/(365.25d0*24.d0*3600.d0))
             if (ttauriwind) call assignDensitiesBlandfordPayne(grid, grid%octreeRoot)
             if (ttauridisc) call assignDensitiesAlphaDisc(grid, grid%octreeRoot)
             if (ttauriwarp) call addWarpedDisc(grid%octreeRoot)
             ! Finding the total mass in the accretion flow
             mass_accretion_old = 0.0d0
             call TTauri_accretion_mass(grid%octreeRoot, grid, mass_accretion_old)
             write(message,*) "Total mass in accretion flow is ",  mass_accretion_old, "[g]"
             call writeInfo(message,FORINFO)

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

       end select


        call countVoxels(grid%octreeRoot,nOctals,nVoxels)
        grid%nOctals = nOctals
        call howmanysplits()

        call writeInfo("Calling routines to finalize the grid variables...",TRIVIAL)
        call finishGrid(grid%octreeRoot, grid, romData=romData)

          if( geometry .eq. 'molcluster'.or. geometry == "theGalaxy" ) then
             totalmasstrap = 0.0; maxrho=0.0; minrho=1.0e30
             call findTotalMass(grid%octreeRoot, totalMass, totalmasstrap = totalmasstrap, maxrho=maxrho, minrho=minrho)
             write(message,*) "Mass of envelope: ",totalMass/mSol, " solar masses"
             call writeInfo(message, TRIVIAL)
             write(message,*) "Mass of envelope (TRAP): ",totalMasstrap/mSol, " solar masses"
             call writeInfo(message, TRIVIAL)
             write(message,*) "Maximum Density: ",maxrho, " g/cm^3"
             call writeInfo(message, TRIVIAL)
             write(message,*) "Minimum Density: ",minrho, " g/cm^3"
             call writeInfo(message, TRIVIAL)
          endif


        call writeInfo("...final adaptive grid configuration complete",TRIVIAL)
        call howmanysplits()

       select case (geometry)
          case("cmfgen")
              call map_cmfgen_opacities(grid)
              call distort_cmfgen(grid%octreeRoot, grid)
              call set_bias_cmfgen(grid%octreeRoot, grid, lamline)
              call writeVtkFile(grid, "cmfgen.vtk",  valueTypeString=(/"etaline ","chiline ","ne      ", &
                   "velocity"/))

          case("wrshell")
             grid%geometry = "wrshell"
             grid%lineEmission = .false.
             grid%rInner = rInner
             grid%rOuter = rOuter
             grid%rCore = rCore

          case DEFAULT
       end select

       call fixParentPointers(grid%octreeRoot)
       call postSetupChecks(grid)
       call writeVTKfile(grid, "rho.vtk")
    endif
#ifdef MPI
        if (grid%splitOverMPI) then
           call grid_info_mpi(grid, "info_grid.dat")
        else
           if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
        endif
        call torus_mpi_barrier
#else
        if ( myRankIsZero ) call grid_info(grid, "info_grid.dat")
#endif
        

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
        call kill_all(young_cluster)
  end subroutine setupamrgrid

  subroutine doSmoothOnTau(grid)

    use input_variables, only: doSmoothGridTau, dustPhysics, lambdaSmooth
    use input_variables, only: photoionPhysics, variableDustSublimation, dosmoothgrid, smoothfactor
    use utils_mod, only: locate
    use lucy_mod, only: putTau
    use grid_mod, only: grid_info

    implicit none

    type(GRIDTYPE) :: grid
    integer :: j, ismoothlam
    logical :: gridConverged, doPhotoSphereSplit
    character(len=80) :: message


    ! Smooth the grid with respect to optical depth, if requested
    if (doSmoothGridTau.and.dustPhysics) then

       call writeVtkFile(grid, "beforesmooth.vtk")

       doPhotoSphereSplit =  ((.not.variableDustSublimation).and.(.not.photoionPhysics))
       if ( doPhotoSphereSplit ) then
          call writeInfo("Doing photsphere split", TRIVIAL)
       else
          call writeInfo("Not doing photsphere split", TRIVIAL)
       end if

       call writeInfo("Smoothing adaptive grid structure for optical depth...", TRIVIAL)
       call locate(grid%lamArray, grid%nLambda,lambdaSmooth,ismoothlam)
       do j = iSmoothLam,  grid%nLambda, 20
          write(message,*) "Smoothing at lam = ",grid%lamArray(j), " angs"
          call writeInfo(message, TRIVIAL)
          do
             gridConverged = .true.
             call putTau(grid, grid%lamArray(j))
             call myTauSmooth(grid%octreeRoot, grid, j, gridConverged, &
                  inheritProps = .false., interpProps = .false., &
                  photosphereSplit = doPhotoSphereSplit )
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
                  gridConverged,  inheritProps = .false., interpProps = .false.)
             if (gridConverged) exit
          end do
          call writeInfo("...grid smoothing complete", TRIVIAL)
       endif
             
       if ( myRankIsZero ) call grid_info(grid, "info_grid_aftersmooth.dat")

    end if

  end subroutine doSmoothOnTau

  subroutine setupFogel(grid, filename, speciesName)
    use input_variables, only : rinner, rOuter
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

    rinner = 5.*autocm/1.d10
    grid%rinner = rinner
    router = 400.*autocm/1.d10
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
    call writeWarning("Need to check whether abundances are really relative to N(H) rather than N(H_2)")
    call fillGridFogel(grid%octreeRoot, grid, r, z, nr, nz, rho, t, abundance)

    call writeVtkFile(grid, "fogel.vtk",  valueTypeString=(/"rho         ",&
         "temperature ","molabundance","microturb   ","velocity    "/))

  end subroutine setupFogel

  recursive subroutine splitGridFogel(thisOctal, grid, r, z, nr, nz, rho, t, abundance)
    use input_variables, only : minDepthAMR, maxDepthAMR
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
          if (thisOctal%subcellSize > (r(i+1)-r(i))) split = .true.
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
          use input_variables, only : mcore, vturb
          type(GRIDTYPE) :: grid
          type(octal), pointer   :: thisOctal
          type(octal), pointer  :: child 
          type(VECTOR) :: rVec
          real(double) :: thisR, thisZ,fac1,fac2,fac3
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

                Mcore = 0.5d0 * mSol
                rVec = subcellCentre(thisOctal, subcell)
                thisR = sqrt(rVec%x**2 + rVec%y**2)
                thisZ = abs(rVec%z)

                thisOctal%velocity(subcell) = keplerianVelocity(rvec)
                CALL fillVelocityCorners(thisOctal,keplerianVelocity)
                outsideGrid = .false.

                if ((thisR < r(1)).or.(thisR > r(nr))) outsideGrid = .true.
                if (.not.outsideGrid) then
                   call locate(r, nr, thisR, j)
                   if ((thisZ < z(j,1)).or.(thisZ > z(j,nz(j)))) outsideGrid = .true.
                endif

                thisOctal%rho(subcell) = 1.d-25
                thisOctal%temperature(subcell) = 3.d0
                thisOctal%molAbundance(subcell) = 1.d-20
                thisOctal%microTurb(subcell) = 0.d0
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


                   thisOctal%temperature(subcell) = max(3.d0,(1.d0-fac1)*(1.d0-fac2)* t(j,k1) + &
                        (     fac1)*(1.d0-fac3)* t(j+1,k2) + &
                        (1.d0-fac1)*(     fac2)* t(j,k1+1) + &
                        (     fac1)*(     fac3)* t(j+1,k2+1) ) 


                   thisOctal%molAbundance(subcell) = (1.d0-fac1)*(1.d0-fac2)* abundance(j,k1) + &
                        (     fac1)*(1.d0-fac3)* abundance(j+1,k2) + &
                        (1.d0-fac1)*(     fac2)* abundance(j,k1+1) + &
                        (     fac1)*(     fac3)* abundance(j+1,k2+1)  

                   thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.d0 * mHydrogen) 

                   thisOctal%molAbundance(subcell) = thisOctal%molAbundance(subcell) * 2.d0
                endif
                thisOctal%microturb(subcell) = max(vturb/(1.d-5*cspeed),sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) &
                        / (28.0 * amu)) + vturb**2) &
                        / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence

                thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.d0 * mHydrogen) 
             endif
          enddo
        end subroutine fillGridFogel


  recursive subroutine splitGridFractal(thisOctal, rho, aFac, grid, converged)
    use input_variables, only : maxDepthAMR
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    real :: rho, aFac
    integer :: subcell, i, j
    real, allocatable :: r(:), s(:)
    real :: rmin, rmax, tot, fac, mean
    logical :: converged, split

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
          split = .false.
          if (octalonThread(thisOctal, subcell, myrankGlobal)) split = .true.
          if (thisOctal%nDepth == maxDepthAMR) split = .false.

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
             enddo
             deallocate(r, s)
             converged = .false.
             exit
          endif
          
       endif
    end do
  end subroutine splitGridFractal



  subroutine postSetupChecks(grid)
    use input_variables, only : mDisc, geometry
    use memory_mod, only : findTotalMemory, reportMemory
    type(GRIDTYPE) :: grid
    integer(kind=bigInt) :: i
    
    call findTotalMemory(grid, i)
    call reportMemory(i)


    select case (geometry)
    case ("shakara")
       call testAMRmass(grid, dble(mdisc))
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
    if (fac > 0.01d0) then
       write(message,'(a,f7.1,a)') "Grid mass differs from required mass by: ",100.d0*fac, " %"
       call writeFatal(message)
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



end module setupamr_mod
