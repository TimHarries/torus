MODULE amr_mod
 
! 21 nov
  ! routines for adaptive mesh refinement. nhs
  ! twod stuff added by tjh started 25/08/04

  USE unix_mod
  USE octal_mod
  use utils_mod
  USE constants_mod
  use density_mod, only:    density, TTauriInFlow
  use romanova_class, only: romanova
  use gridtype_mod, only:   gridtype, hydrospline
  use surface_mod, ONLY:    SURFACETYPE
  USE cluster_class, only:  cluster
  USE parallel_mod, ONLY:   torus_abort
  use mpi_global_mod

  IMPLICIT NONE

  type STREAMTYPE
     integer :: nSamples
     type(VECTOR),pointer :: position(:) => null()
     type(VECTOR), pointer :: direction(:) => null()
     type(VECTOR), pointer :: velocity(:) => null()
     real, pointer :: speed(:) => null()
     real(double),pointer :: distanceAlongStream(:)   => null()
     real(double), pointer :: rho(:)  => null()
     real, pointer :: temperature(:)  => null()
     real(double),pointer :: streamRadius(:)  => null()
  end type STREAMTYPE

  integer :: somecounter

  type(STREAMTYPE),save :: globalStream(2000)
  integer,save :: globalnStream

  type curtaintype
     integer :: nr, ntheta
     type(VECTOR), pointer :: position(:,:) => null()
     type(VECTOR), pointer :: velocity(:,:) => null()
     real(double), pointer :: density(:,:) => null()
     real(double), pointer :: temperature(:,:) => null()
  end type curtaintype

  real(double), parameter :: amr_min_rho = 1.0e-30_db 

  integer :: mass_split, mass_split2, density_split, velocity_split, both_split, maxdensity_split
  integer :: scaleheighta_count, scaleheightb_count, scaleheightc_count
  TYPE(octal), POINTER :: recentOctal

CONTAINS

  SUBROUTINE calcValuesAMR(thisOctal,subcell,grid, stellar_cluster, inherit, interp, &
       romData, stream)
    ! calculates the variables describing one subcell of an octal.
    ! each geometry that can be used with AMR should be described here, 
    !   otherwise the program will print a warning and exit.
   
    use cluster_class, only:  assign_density
    use luc_cir3d_class, only: calc_cir3d_mass_velocity
    use cmfgen_class, only: cmfgen_mass_velocity
    use jets_mod, only: calcJetsMassVelocity
    use romanova_class, only: calc_romanova_mass_velocity

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT)    :: thisOctal   ! the octal being changed
    INTEGER, INTENT(IN)           :: subcell   ! the subcell being changed
    TYPE(gridtype), INTENT(INOUT) :: grid      ! the grid
    !
    type(STREAMTYPE), optional :: stream
    TYPE(cluster), optional, INTENT(IN)   :: stellar_cluster
    LOGICAL, OPTIONAL :: inherit               ! inherit densities, temp, etc of parent
    LOGICAL, OPTIONAL :: interp                ! interpolate densities, temp, etc of parent
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
    !
    TYPE(octal), POINTER :: parentOctal
!    real(double) :: rhoDouble, r
    INTEGER :: parentSubcell
!    type(VECTOR) :: rVec
    LOGICAL :: inheritProps
    LOGICAL :: interpolate

    inheritProps = .false.
    if (present(inherit)) then
       inheritProps = inherit
    endif

    interpolate = .false.
    if (present(interp)) then
       interpolate = interp
    endif
       parentOctal => thisOctal%parent
       parentSubcell = thisOctal%parentSubcell

       if (associated (thisOctal%boundaryCondition)) &
            thisOCtal%boundaryCondition(subcell) = parentOctal%boundaryCondition(parentSubcell)
       if (associated(thisOctal%gamma)) thisOctal%gamma(subcell)  = parentOctal%gamma(parentSubcell)
       if (associated(thisOctal%iEquationOfState)) &
            thisOctal%iEquationOfState(subcell)  = parentOctal%iEquationOfState(parentSubcell)
    if (inheritProps) then
       thisOctal%rho(subcell) = parentOctal%rho(parentSubcell)
       thisOctal%temperature(subcell) = parentOctal%temperature(parentSubcell)
       if (associated(thisOctal%etaCont)) thisOctal%etaCont(subcell) = parentOctal%etaCont(parentSubcell)
       thisOctal%inFlow(subcell) = parentOctal%inFlow(parentSubcell)
       thisOctal%velocity(subcell) = parentOctal%velocity(parentSubcell)
       if (associated(thisOctal%biasCont3d)) thisOctal%biasCont3D(subcell) = parentOctal%biasCont3D(parentSubcell)
       if (associated(thisOctal%etaLine)) thisOctal%etaLine(subcell) = parentOctal%etaLine(parentSubcell)
       if (associated(thisOctal%chiLine)) thisOctal%chiLine(subcell) = parentOctal%chiLine(parentSubcell)
       if (associated(thisOctal%boundaryCondition)) thisOCtal%boundaryCondition(subcell) &
            = parentOctal%boundaryCondition(parentSubcell)
       if (associated(thisOctal%dustTypeFraction)) then
          thisOctal%dustTypeFraction(subcell,:) = parentOctal%dustTypeFraction(parentSubcell,:)
       endif
       if (associated(thisOctal%oldFrac)) thisOctal%oldFrac(subcell) = parentOctal%oldFrac(parentSubcell)
       if (associated(thisOctal%ionFrac)) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif
       if (associated(thisOctal%nh)) thisOctal%nh(subcell) = parentOctal%nh(parentsubcell)
       if (associated(thisOctal%ne)) thisOctal%ne(subcell) = parentOctal%ne(parentsubcell)

       if (associated(thisOctal%rhou)) thisOctal%rhou(subcell) = parentOctal%rhou(parentSubcell)
       if (associated(thisOctal%rhov)) thisOctal%rhov(subcell) = parentOctal%rhov(parentSubcell)
       if (associated(thisOctal%rhow)) thisOctal%rhow(subcell) = parentOctal%rhow(parentSubcell)
       if (associated(thisOctal%rhoe)) thisOctal%rhoe(subcell) = parentOctal%rhoe(parentSubcell)
       if (associated(thisOctal%phi_i)) thisOctal%phi_i(subcell) = parentOctal%phi_i(parentSubcell)
       if (associated(thisOctal%energy)) thisOctal%energy(subcell) = parentOctal%energy(parentSubcell)

    else if (interpolate) then
       if (associated(thisOctal%boundaryCondition)) &
            thisOCtal%boundaryCondition(subcell) = parentOctal%boundaryCondition(parentSubcell)
       Thisoctal%etacont(subcell) = parentOctal%etaCont(parentSubcell)
       thisOctal%inFlow(subcell) = parentOctal%inFlow(parentSubcell)
       thisOctal%velocity(subcell) = parentOctal%velocity(parentSubcell)
       thisOctal%biasCont3D(subcell) = parentOctal%biasCont3D(parentSubcell)
       thisOctal%etaLine(subcell) = parentOctal%etaLine(parentSubcell)
       thisOctal%chiLine(subcell) = parentOctal%chiLine(parentSubcell)
!       thisOctal%dustTypeFraction(subcell,:) = parentOctal%dustTypeFraction(parentSubcell,:)
       thisOctal%oldFrac(subcell) = parentOctal%oldFrac(parentSubcell)
       if (associated(thisOctal%ionFrac)) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif
       if (associated(thisOctal%nh)) &
            thisOctal%nh(subcell) = parentOctal%nh(parentsubcell)
       if (associated(thisOctal%ne)) &
            thisOctal%ne(subcell) = parentOctal%ne(parentsubcell)

       parentOctal%hasChild(parentsubcell) = .false.
       call interpFromParent(subcellCentre(thisOctal, subcell), thisOctal%subcellSize, grid, &
            thisOctal%temperature(subcell), thisOctal%rho(subcell), thisOctal%dusttypeFraction(subcell, :), &
            thisOctal%etaLine(subcell))
       parentOctal%hasChild(parentsubcell) = .true.

    else

    SELECT CASE (grid%geometry)

    CASE("pathtest")
       call calcPathTestDensity(thisOctal,subcell,grid)

    CASE ("ttauri")
      CALL calcTTauriMassVelocity(thisOctal,subcell,grid)
      
    CASE ("jets")
      CALL calcJetsMassVelocity(thisOctal,subcell,grid)

    CASE ("luc_cir3d")
      CALL calc_cir3d_mass_velocity(thisOctal,subcell)

    CASE ("cmfgen")
      CALL cmfgen_mass_velocity(thisOctal,subcell)

    CASE ("romanova")
      CALL calc_romanova_mass_velocity(romData, thisOctal,subcell)

    CASE ("testamr")
       CALL calcTestDensity(thisOctal,subcell,grid)

    CASE("lexington")
       CALL calcLexington(thisOctal, subcell, grid)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("starburst")
       CALL calcStarburst(thisOctal, subcell, grid)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("symbiotic")
       CALL calcSymbiotic(thisOctal, subcell, grid)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("proto")
       CALL calcProtoDensity(thisOctal,subcell,grid)

    CASE("wrshell")
       CALL calcWRShellDensity(thisOctal,subcell,grid)

    CASE("hydro1d")
       call calcHydro1DDensity(thisOctal, subcell, grid)

    CASE("kelvin")
       call calcKelvinDensity(thisOctal, subcell, grid)

    CASE("bonnor")
       call calcBonnorEbertDensity(thisOctal, subcell, grid)

    CASE("unisphere")
       call calcUniformsphere(thisOctal, subcell, grid)

    CASE("sedov")
       call calcSedovDensity(thisOctal, subcell, grid)

    CASE("protobin")
       call calcProtoBinDensity(thisOctal, subcell, grid)


    CASE("gammavel")
       CALL calcGammaVel(thisOctal,subcell,grid)

    CASE ("spiralwind")
       CALL spiralWindSubcell(thisOctal, subcell ,grid)
       
    CASE("cluster","molcluster","theGalaxy")
       ! using a routine in cluster_class.f90
!       call assign_density(thisOctal,subcell, grid%geometry, stellar_cluster)

       thisoctal%rho(subcell) = -9.9d99
       thisoctal%cornervelocity = VECTOR(-9.9d99,-9.9d99,-9.9d99)

    CASE("wr104")
       ! using a routine in cluster_class.f90
       call assign_density(thisOctal,subcell, grid%geometry)

    CASE ("benchmark")
       CALL benchmarkDisk(thisOctal, subcell ,grid)

    CASE ("molebench")
       thisoctal%rho(subcell) = -9.9d99
       CALL molecularBenchmark(thisOctal, subcell ,grid)

    CASE ("h2obench1")
       CALL WaterBenchmark1(thisOctal, subcell)

    CASE ("h2obench2")
       CALL WaterBenchmark2(thisOctal, subcell, grid)

    CASE ("agbstar")
       CALL AGBStarBenchmark(thisOctal, subcell ,grid)

    CASE ("molefil")
       CALL molecularFilamentFill(thisOctal, subcell ,grid)

    CASE ("ggtau")
       CALL ggtauFill(thisOctal, subcell ,grid)

    CASE ("shakara","aksco","circumbin")
       CALL shakaraDisk(thisOctal, subcell ,grid)

    CASE ("iras04158")
       CALL iras04158(thisOctal, subcell, grid)

    CASE ("warpeddisc")
       CALL warpedDisk(thisOctal, subcell ,grid)

    CASE("ppdisk")
       CALL calcPPDiskDensity(thisOctal,subcell,grid)

    CASE("melvin")
       CALL assign_melvin(thisOctal,subcell,grid)

    CASE("whitney")
       CALL assign_whitney(thisOctal,subcell,grid)

    CASE("planetgap")
       CALL assign_planetgap(thisOctal,subcell,grid)

    CASE("toruslogo")
       CALL assign_toruslogo(thisOctal,subcell,grid)

    CASE("clumpydisc")
       CALL assign_clumpydisc(thisOctal, subcell, grid)

    CASE ("windtest")
      CALL calcWindTestValues(thisOctal,subcell,grid)

   CASE ("fractal")
      thisOctal%rho = 100.d0 * mHydrogen
      thisOctal%temperature = 8000.


    CASE ("magstream")

!      thisOctal%rho(subcell) = 1.d-30
!      thisOctal%temperature(subcell) = 6000.
!
!       rVec = subcellCentre(thisOctal, subcell)
!       if (iGlobalSample /= 0) then
!          r = modulus(rVec - globalStream%position(iGlobalSample))
!          if (r < globalStream%streamradius(iGlobalSample)) then
!             thisOctal%rho(subcell) = globalStream%rho(iGlobalSample)
!             write(*,*) thisOctal%rho(subcell),iglobalsample
!             thisOctal%temperature(subcell) = 7500.
!          endif
!       endif


!      CALL getMagStreamValues2(point=subcellCentre(thisOctal,subcell),&
!                              grid=grid,                            &
!                              rho=rhoDouble,                        &
!                              temperature=thisOctal%temperature(subcell),&
!                              velocity=thisOctal%velocity(subcell),  &
!                              inFlow=thisOctal%inFlow(subcell))

      call getMagStreamValues3(thisOctal, subcell, stream, grid)


!      thisOctal%rho(subcell) = REAL(rhoDouble)
!      IF (subcell == thisOctal%maxChildren) CALL fillVelocityCorners(thisOctal,grid,magStreamVelocity, .true.)
!       thisOctal%microturb(subcell) = (20.d5/cSpeed) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

    CASE DEFAULT
      WRITE(*,*) "! Unrecognised grid geometry in calcValuesAMR: ",TRIM(grid%geometry)
      STOP

    END SELECT
 
!    CALL fillGridDummyValues(thisOctal,subcell, grid)
 
    end if
 
  END SUBROUTINE calcValuesAMR

  SUBROUTINE initFirstOctal(grid, centre, size, oned, twod, threed, stellar_cluster, nDustType, romData ,&
       stream)
    ! creates the first octal of a new grid (the root of the tree).
    ! this should only be used once; use addNewChild for subsequent
    !  additions.
    use input_variables, only : cylindrical

    IMPLICIT NONE
    type(OCTAL), pointer :: thisOctal
    TYPE(gridtype), INTENT(INOUT)    :: grid 
    TYPE(vector), INTENT(IN)    :: centre ! coordinates of the grid centre
    REAL, INTENT(IN)                 :: size 
      ! 'size' should be the vertex length of the cube that contains the whole
      !   of the simulation space, *not* the size of a subcell.
    TYPE(cluster), optional, INTENT(IN)    :: stellar_cluster
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
    type(streamtype), optional :: stream
    !
    LOGICAL :: oned, twod, threed  ! true if this is a twoD amr grid
    INTEGER :: subcell ! loop counter 
    INTEGER :: nDustType ! number of different dust types
    integer :: i

    character(len=100) :: message

    ALLOCATE(grid%octreeRoot)

    do i = 1, 8
       grid%octreeRoot%mpiThread(i) = i 
    enddo
    
    ! allocate any variables that need to be 
    if (.not.grid%oneKappa) then
       ALLOCATE(grid%octreeRoot%kappaAbs(8,grid%nOpacity))
       ALLOCATE(grid%octreeRoot%kappaSca(8,grid%nOpacity))
       grid%octreeRoot%kappaAbs = 1.e-30
       grid%octreeRoot%kappaSca = 1.e-30
    endif


    if (oneD) then
       grid%octreeRoot%oneD = .true.
       grid%octreeRoot%twoD = .false.
       grid%octreeRoot%threeD = .false.
       grid%octreeRoot%maxChildren = 2
    endif

    if (twoD) then
       grid%octreeRoot%oneD = .false.
       grid%octreeRoot%twoD = .true.
       grid%octreeRoot%threeD = .false.
       grid%octreeRoot%maxChildren = 4
    else if (threed) then
       grid%octreeRoot%cylindrical = .false.
       grid%octreeRoot%oneD = .false.
       grid%octreeRoot%twoD = .false.
       grid%octreeRoot%threeD = .true.
       grid%octreeRoot%maxChildren = 8
       if (cylindrical) then
          grid%octreeRoot%oneD = .false.
          grid%octreeRoot%twoD = .false.
          grid%octreeRoot%threeD = .true.
          grid%octreeRoot%cylindrical = .true.
          grid%octreeRoot%splitAzimuthally = .false.
          grid%octreeRoot%maxChildren = 4
          grid%octreeRoot%phi = pi
          grid%octreeRoot%dPhi = twoPi
       endif
    endif


    grid%octreeRoot%nDepth = 1
    grid%octreeRoot%nChildren = 0
    grid%octreeRoot%hasChild = .FALSE.
    grid%octreeRoot%subcellSize = size/2.0_oc
    grid%octreeRoot%centre = centre
    if (cylindrical) then
       grid%octreeRoot%centre%x = -size / 2.d0
       grid%octreeRoot%r        = size / 2.d0
       grid%octreeRoot%centre%y = 0.d0
       grid%octreeRoot%centre%z = 0.d0
    endif
    grid%octreeRoot%xMin = grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize
    grid%octreeRoot%yMin = grid%octreeRoot%centre%y - grid%octreeRoot%subcellSize
    grid%octreeRoot%zMin = grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize
    grid%octreeRoot%xMax = grid%octreeRoot%centre%x + grid%octreeRoot%subcellSize
    grid%octreeRoot%yMax = grid%octreeRoot%centre%y + grid%octreeRoot%subcellSize
    grid%octreeRoot%zMax = grid%octreeRoot%centre%z + grid%octreeRoot%subcellSize
    grid%octreeRoot%indexChild = -999 ! values are undefined
    NULLIFY(grid%octreeRoot%parent)   ! tree root does not have a parent
    NULLIFY(grid%octreeRoot%child)    ! tree root does not yet have children
    ! initialize some values
    grid%octreeRoot%rho = 1.e-30

    thisOctal => grid%octreeRoot
    call allocateOctalAttributes(grid, thisOctal)

    select case (grid%geometry)
       case("cluster","molcluster","theGalaxy")
          ! Initially we copy the idecies of particles (in SPH data)
          ! to the root node. The indecies will copy over to
          ! to the subcells if the particles are in the subcells.
          ! This will allow us to work with the subsets of gas particle
          ! list hence reduces the computation time when we are 
          ! splitting/constructing the octree data structure. 
          !
          ! Using the routine in grid_mod.f90

          write(message,*) "Copying SPH index to root"
          call writeinfo(message, TRIVIAL)
          call copy_sph_index_to_root(grid)
          write(message,*) "Initializing OctreeRoot"
          call writeinfo(message, TRIVIAL)
          !
          DO subcell = 1, grid%octreeRoot%maxChildren
            ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, stellar_cluster)
            ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO

          write(message,*) "Done."
          call writeinfo(message, TRIVIAL)

       case("wr104")
          ! Using the routine in grid_mod.f90
          call copy_sph_index_to_root(grid)
          !
          DO subcell = 1, grid%octreeRoot%maxChildren
             ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, stellar_cluster)
             ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO
       case("romanova")
          ! just checking..
          if (.not. PRESENT(romData)) then
             print *, "Error:: romData is not present in initFirstOctal!"
             stop
          end if
          DO subcell = 1, grid%octreeRoot%maxChildren
             ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, romData=romData)
             ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO

       case("magstream")
          ! just checking..
          if (.not. PRESENT(stream)) then
             print *, "Error:: stream is not present in initFirstOctal!"
             stop
          end if
          DO subcell = 1, grid%octreeRoot%maxChildren
             ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, stream=stream)
             ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO
       case DEFAULT

!!!!!!!!!! EDITTED OUT BY TJH
!          DO subcell = 1, grid%octreeRoot%maxChildren
            ! calculate the values at the centre of each of the subcells
!             CALL calcValuesAMR(grid%octreeRoot,subcell,grid)
            ! label the subcells
!             grid%octreeRoot%label(subcell) = subcell
!          END DO
       end select
    


    ! we keep track of the maximum depth of the grid...
    grid%maxDepth = 1
    ! ...and the size of the smallest subcells in the grid.
    ! we will actually store the value which is half the size of the 
    !   smallest subcell because this is more useful for later
    !   calculations.
    grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
                                2.0_oc**REAL(grid%maxDepth,kind=oct)

  END SUBROUTINE initFirstOctal


  SUBROUTINE addNewChild(parent, iChild, grid, adjustGridInfo, &
                         stellar_cluster, inherit, interp, splitAzimuthally, romData, &
                         isample, stream)
    ! adds one new child to an octal

    USE input_variables, ONLY : cylindrical
    USE cluster_class, only: update_particle_list
    USE sph_data_class, only: isAlive
    use mpi_global_mod, only: nThreadsGlobal

    IMPLICIT NONE
    
    TYPE(octal), TARGET, INTENT(INOUT) :: parent ! the parent octal 
    type(octal), pointer :: thisOctal
    INTEGER, INTENT(IN)  :: iChild     ! the label (1-8) of the subcell gaining the child 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that
                                          !   calculate the variables stored in
                                          !   the tree.
    type(STREAMTYPE), optional :: stream
    integer, optional :: isample
    LOGICAL, INTENT(IN) :: adjustGridInfo
    LOGICAL, optional :: splitAzimuthally
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell    
    TYPE(cluster), OPTIONAL, INTENT(IN)   :: stellar_cluster
    LOGICAL, OPTIONAL :: inherit       ! inherit densities, temps, etc from parent
    LOGICAL, OPTIONAL :: interp        ! interpolate densities, temps, etc from parent    
    
    INTEGER       :: subcell           ! loop counter
    INTEGER, SAVE :: counter = 9       ! keeps track of the current subcell label
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: nChildren         ! number of children the parent octal has
    INTEGER       :: newChildIndex     ! the storage location for the new child
    integer :: i
    logical :: inheritProps, interpolate
    type(VECTOR) :: rVec
    ! array of octals that may be needed for temporarily storing child octals

    ! For "romanova" geometry
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
    




    inheritProps = .false.
    if (present(inherit)) then
       inheritProps = inherit
    endif

    interpolate = .false.
    if (present(interp)) then
       interpolate = interp
    endif


    ! store the number of children that already exist
    nChildren = parent%nChildren

    ! safety checks of child array
    IF ( ASSOCIATED(parent%child) ) THEN
      IF ( ( nChildren == 0 ) .OR.                  &
           ( nChildren /= SIZE(parent%child) ) ) THEN
        PRINT *, 'Panic: in addNewChild, %child array wrong size'
        PRINT *, 'nChildren:',nChildren,' SIZE %child:', SIZE(parent%child)
        STOP
      END IF
    END IF
    IF ( (.NOT. ASSOCIATED(parent%child)) .AND. (nChildren > 0) ) THEN
      PRINT *, 'Panic: in addNewChild, %child array wrong size'
      PRINT *, 'nChildren:',nChildren,' ASSOCIATED %child:', ASSOCIATED(parent%child)
      STOP
    END IF

    ! check that new child does not already exist
    IF ( parent%hasChild(iChild) .EQV. .TRUE. ) THEN
      PRINT *, 'Panic: in addNewChild, attempted to add a child ',&
               '       that already exists'
      STOP
    ENDIF

!    call interpFromParent(subcellCentre(parent, iChild, parent%subcellSize, &
!         grid, temperature, density, dusttypeFraction)

    CALL growChildArray(parent, nNewChildren=1, grid=grid )

    ! update the bookkeeping
    nChildren = nChildren + 1
    newChildIndex = nChildren

    ! update the parent octal
    parent%nChildren = nChildren
    parent%hasChild(iChild) = .TRUE.
    parent%indexChild(newChildIndex) = iChild

    ! allocate any variables that need to be  
    IF (.NOT.grid%oneKappa) THEN
       ! The kappa arrays should be allocated with grid%nopacity instead of grid%nlambda
       ! because for line calculation, there is only one kappa needed.
       ! (but grid%nlambda is not 1). If you allocate the arrays with grid%nlambda,
       ! it will be a huge waste of RAM. ---  (RK) 
       ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nopacity))
       ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nopacity))
       ! ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nlambda))
       ! ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nlambda))
       parent%child(newChildIndex)%kappaAbs = 1.e-30
       parent%child(newChildIndex)%kappaSca = 1.e-30
    ENDIF
    NULLIFY(parent%child(newChildIndex)%child)


    parent%child(newChildIndex)%nDepth = parent%nDepth + 1

! setup mpiThread values



    if ( ((parent%twoD)  .and.((nThreadsGlobal - 1) == 4)) .or. &
         ((parent%threed).and.((nThreadsGlobal - 1) == 8)).or. &
         ((parent%oneD)  .and.((nThreadsGlobal - 1) == 2)) ) then
       parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
    else

       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          if (parent%oneD) then
             do i = 1, 2
                parent%child(newChildIndex)%mpiThread(i) = 2 * (parent%mpiThread(iChild) - 1) + i
             enddo
          else if (parent%twoD) then
             do i = 1, 4
                parent%child(newChildIndex)%mpiThread(i) = 4 * (parent%mpiThread(iChild) - 1) + i
             enddo
          else if (parent%Threed) then
             do i = 1, 8
                parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
             enddo
          endif
       endif
    endif

    if ((parent%threed).and.(nThreadsGlobal - 1) == 64) then
       if (parent%child(newChildIndex)%nDepth > 2) then
          parent%child(newChildIndex)%mpiThread = parent%mpiThread(iChild)
       else
          do i = 1, 8
             parent%child(newChildIndex)%mpiThread(i) = 8 * (parent%mpiThread(iChild) - 1) + i
          enddo
       endif
    endif

    ! set up the new child's variables
    parent%child(newChildIndex)%threeD = parent%threeD
    parent%child(newChildIndex)%twoD = parent%twoD
    parent%child(newChildIndex)%oneD = parent%oneD
    parent%child(newChildIndex)%maxChildren = parent%maxChildren
    parent%child(newChildIndex)%cylindrical = parent%cylindrical


    
    ! if splitAzimuthally is not present then we assume we are not

    if (cylindrical) then  
       if (parent%splitAzimuthally) then
          rVec =  subcellCentre(parent,iChild)
          parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
          if (parent%child(newChildIndex)%phi < 0.d0) then
              parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
          endif
          parent%child(newChildIndex)%dphi = parent%dphi/2.d0
       else
          parent%child(newChildIndex)%phi = parent%phi
          parent%child(newChildIndex)%dphi = parent%dphi
       endif
       parent%child(newChildIndex)%splitAzimuthally = .false.
       parent%child(newChildIndex)%maxChildren = 4

       if (PRESENT(splitAzimuthally)) then
          if (splitAzimuthally) then
             parent%child(newChildIndex)%splitAzimuthally = .true.
             parent%child(newChildIndex)%maxChildren = 8
             rVec =  subcellCentre(parent,iChild)
             parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
             if (parent%child(newChildIndex)%phi < 0.d0) then
                parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
             endif
             if (parent%splitAzimuthally) then
                parent%child(newChildIndex)%dphi = parent%dphi/2.d0
             else
                parent%child(newChildIndex)%dphi = parent%dphi
             endif
          else
             if (parent%splitAzimuthally) then
                rVec =  subcellCentre(parent,iChild)
                parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
                if (parent%child(newChildIndex)%phi < 0.d0) then
                   parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
                endif
                parent%child(newChildIndex)%dphi = parent%dphi/2.d0
             else
                parent%child(newChildIndex)%phi = parent%phi
                parent%child(newChildIndex)%dphi = parent%dphi
             endif
             parent%child(newChildIndex)%splitAzimuthally = .false.
             parent%child(newChildIndex)%maxChildren = 4
          endif
       endif
    endif


    parent%child(newChildIndex)%inFlow = parent%inFlow
    parent%child(newChildIndex)%parent => parent
    parent%child(newChildIndex)%parentSubcell = iChild
    parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
    parent%child(newChildIndex)%hasChild = .false.
    parent%child(newChildIndex)%nChildren = 0
    parent%child(newChildIndex)%indexChild = -999 ! values are undefined
    parent%child(newChildIndex)%centre = subcellCentre(parent,iChild)
    if (parent%cylindrical) then
       parent%child(newChildIndex)%r = subcellRadius(parent,iChild)
    endif

    parent%child(newChildIndex)%xMin = parent%child(newChildIndex)%centre%x - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMin = parent%child(newChildIndex)%centre%y - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMin = parent%child(newChildIndex)%centre%z - parent%child(newChildIndex)%subcellSize

    parent%child(newChildIndex)%xMax = parent%child(newChildIndex)%centre%x + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMax = parent%child(newChildIndex)%centre%y + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMax = parent%child(newChildIndex)%centre%z + parent%child(newChildIndex)%subcellSize


    thisOctal => parent%child(newChildIndex)
    call allocateOctalAttributes(grid, thisOctal)


    if (isAlive()) then
       ! updates the sph particle list.           
       CALL update_particle_list(parent, iChild, newChildIndex)
    endif
    
    ! put some data in the four/eight subcells of the new child
    DO subcell = 1, parent%child(newChildIndex)%maxChildren
      CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid, &
                         stellar_cluster=stellar_cluster, inherit=inheritProps, &
                         interp=interpolate, romData=romData, stream=stream)
      parent%child(newChildIndex)%label(subcell) = counter
      counter = counter + 1
    END DO

    IF ( adjustGridInfo ) THEN
      grid%nOctals = grid%nOctals + 1

      ! check for a new maximum depth 
      IF (parent%child(newChildIndex)%nDepth > grid%maxDepth) THEN
        grid%maxDepth = parent%child(newChildIndex)%nDepth
        CALL setSmallestSubcell(grid)
      END IF
    END IF

    IF (( COUNT(parent%hasChild(:)) /= parent%nChildren ) .OR. &     
        (  SIZE(parent%child(:))    /= parent%nChildren )) THEN
        PRINT *, "Problem in addNewChild"
        PRINT *, "  nchildren: ",parent%nChildren
        PRINT *, "  haschild: ",parent%hasChild
        PRINT *, "  indexchild: ",parent%indexchild
        PRINT *, "  ASSOCIATED(child): ",ASSOCIATED(parent%child)
        PRINT *, "  SIZE(child): ",SIZE(parent%child)
        do;enddo
    END IF

    !CALL checkAMRgrid(grid,checkNoctals=.FALSE.)

  END SUBROUTINE addNewChild
  
  SUBROUTINE growChildArray(parent, nNewChildren, grid)
    ! adds storage space for new children to an octal.
    ! 
    ! some of the bookkeeping must be done by the routine that calls
    !   this one (e.g. making sure the new children are indexed).

    IMPLICIT NONE
    
    TYPE(octal), TARGET, INTENT(INOUT) :: parent ! the parent octal 
    INTEGER, INTENT(IN) :: nNewChildren ! number of children to add
    TYPE(GRIDTYPE), INTENT(INOUT) :: grid
    
    TYPE(wrapperArray):: tempChildStorage 
      ! holder for existing children, while we shuffle them around to 
      !   make room for new ones.
                                       
    INTEGER       :: iChild            ! loop counter
    INTEGER       :: nChildren         ! number of children the parent octal has
    INTEGER       :: error
   
    ! store the number of children that already exist
    nChildren = parent%nChildren

    IF ( nChildren == 0 ) THEN
      ! if there are no existing children, we can just allocate
      ! the 'child' array with the new size
      ALLOCATE(parent%child(nNewChildren), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed in growChildArray.(A)'
        STOP
      END IF

    ELSE ! there are existing children
     
      ! check that there is not a full quota of children
      IF ( (nChildren + nNewChildren) > parent%maxChildren ) THEN
        PRINT *, 'Panic: in growChildArray, attempted to have too ',&
                 '       many children: ',(nChildren + nNewChildren) 
        STOP
      ENDIF
      
      ! if there are existing children, we must enlarge the allocated   
      ! array. we need to use temporary octals[1] and copy the  
      ! existing children into them[2]; then increase the 'child' array size 
      ! [3]; then copy the children back in[4]. 

      ! [1]
      ALLOCATE(tempChildStorage%wrappers(nChildren), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed in growChildArray. (B)'
        STOP
      END IF     

      ! [2]
      DO iChild = 1, nChildren, 1
        
        ALLOCATE(tempChildStorage%wrappers(iChild)%content, STAT=error)
        IF ( error /= 0 ) THEN
          PRINT *, 'Panic: allocation failed in growChildArray. (C)'
          STOP
        END IF           
        tempChildStorage%wrappers(iChild)%inUse = .TRUE.
               
        CALL deleteOctreeBranch(parent%child(iChild),                   &
               onlyChildren=.FALSE.,                                    &
               deletedBranch=tempChildStorage%wrappers(iChild)%content, &
               adjustParent=.FALSE., grid=grid, adjustGridInfo=.FALSE. )
      END DO

      ! [3]
      IF ( ASSOCIATED(parent%child) ) THEN
        DEALLOCATE(parent%child)
        NULLIFY(parent%child)
        ALLOCATE(parent%child( nChildren + nNewChildren ), STAT=error)
        IF ( error /= 0 ) THEN
          PRINT *, 'Panic: allocation failed in growChildArray. (D)'
          STOP
        END IF
      ELSE
        PRINT *, "Error in growChildArray:"
        PRINT *, "  parent%child not associated."
        STOP
      END IF

      ! [4]
      DO iChild = 1, nChildren, 1
        
        CALL insertOctreeBranch(parent%child(iChild),               &
               branch=tempChildStorage%wrappers(iChild)%content,    &
               onlyChildren=.FALSE.)

        parent%child(iChild)%parent => parent       
               
        if (associated(tempChildStorage%wrappers(iChild)%content)) then
           DEALLOCATE(tempChildStorage%wrappers(iChild)%content, STAT=error)
           if (error /= 0) then
              write(*,*) "error",error
              write(*,*) ichild
              write(*,*) tempChildStorage%wrappers(ichild)%content%ndepth
              parent%child(iChild)%rho = 1.d30
!              do;enddo
           endif
        endif

        NULLIFY(tempChildStorage%wrappers(iChild)%content)
        tempChildStorage%wrappers(iChild)%inUse = .FALSE.

      END DO
      
      ! can now clean up the temporary storage
      DEALLOCATE(tempChildStorage%wrappers)
      NULLIFY(tempChildStorage%wrappers)
      
    END IF ! ( nChildren == 0 )

  END SUBROUTINE growChildArray


  RECURSIVE SUBROUTINE splitGrid(thisOctal,amrLimitScalar,amrLimitScalar2,grid,&
 stellar_cluster, setChanged, romData)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    IMPLICIT NONE


    TYPE(OCTAL), intent(inout) :: thisOctal
!    TYPE(OCTAL), POINTER :: thisOctal
    real(double), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
                                          !   routines that this subroutine calls

    ! Object containg the output from the (Mattew's) SPH code.
    ! This will be just passed to decideSplit routine.
    TYPE(cluster), OPTIONAL,  intent(in) :: stellar_cluster
    LOGICAL, INTENT(IN), OPTIONAL :: setChanged
    !
    INTEGER              :: iSubcell, iIndex ! loop counters
    INTEGER              :: i, j, k
    logical :: splitInAzimuth

    !
    ! For "romanova" geometry
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry



    DO iSubcell = 1, thisOctal%maxChildren
      
      IF (decideSplit(thisOctal,iSubcell,amrLimitScalar,amrLimitScalar2,grid,&
            splitInAzimuth, &
            stellar_cluster=stellar_cluster, romData=romData)) THEN

        CALL addNewChild(thisOctal, iSubcell, grid, adjustGridInfo=.TRUE., &
                         stellar_cluster=stellar_cluster, &
                         splitAzimuthally=splitInAzimuth, romData=romData)
       
        if (.not.thisOctal%hasChild(isubcell)) then
          write(*,*) "add child failed in splitGrid"
          do
          enddo
       endif
        
       IF (PRESENT(setChanged)) THEN
          IF (setChanged) THEN
             
             DO iIndex = 1, thisOctal%nChildren, 1
                IF (thisOctal%indexChild(iIndex) == iSubcell) &
                     thisOctal%child(iIndex)%changed = .TRUE.
             END DO
             
          END IF
       END IF
        
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
      
      CALL splitGrid(thisOctal%child(iIndex),amrLimitScalar,amrLimitScalar2,grid,&
                     stellar_cluster, setChanged, romData=romData)
      
   END DO

!    if (associated(thisOctal%gas_particle_list)) then
!       DEALLOCATE(thisOctal%gas_particle_list)
!    endif


  END SUBROUTINE splitGrid
  
  
  RECURSIVE SUBROUTINE getOctalArray(thisOctal,array,counter,maxOctals) 
    ! returns an array of pointers to all of the subcells in the grid.
    ! NB because fortran cannot create arrays of pointers, the output
    !   array is actually of a derived type which *contains* the 
    !   pointer to an octal.
    ! counter should be set to 0 before this routine is called

    IMPLICIT NONE

    TYPE(octal), POINTER                            :: thisOctal
    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    INTEGER, INTENT(INOUT)                          :: counter 
    INTEGER, INTENT(IN), OPTIONAL                   :: maxOctals
  
    INTEGER              :: i
    TYPE(octal), POINTER :: child

    ! if this is the root of the tree, we initialize the counter
    IF (.NOT. ASSOCIATED(thisOctal%parent)) counter = 0
    
    counter = counter + 1 
    array(counter)%content => thisOctal
    !array(counter)%inUse = .TRUE. 
    array(counter)%inUse = .NOT. thisOctal%hasChild 
    IF (PRESENT(maxOctals)) THEN
      IF (counter >= maxOctals) RETURN
    END IF
    
    IF ( thisOctal%nChildren > 0 ) THEN
      DO i = 1, thisOctal%nChildren, 1
        
        ! call this subroutine recursively on each of its children
        child => thisOctal%child(i)
        CALL getOctalArray(child,array,counter)
        
      END DO
    END IF

  END SUBROUTINE getOctalArray

  RECURSIVE SUBROUTINE getOctalArrayLevel(thisOctal,array,counter,ndepth,maxOctals) 
    ! returns an array of pointers to all of the subcells in the grid.
    ! NB because fortran cannot create arrays of pointers, the output
    !   array is actually of a derived type which *contains* the 
    !   pointer to an octal.
    ! counter should be set to 0 before this routine is called

    IMPLICIT NONE

    TYPE(octal), POINTER                            :: thisOctal
    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    INTEGER, INTENT(INOUT)                          :: counter 
    INTEGER, INTENT(IN), OPTIONAL                   :: maxOctals
    integer :: nDepth
  
    INTEGER              :: i
    TYPE(octal), POINTER :: child

    ! if this is the root of the tree, we initialize the counter
    IF (.NOT. ASSOCIATED(thisOctal%parent)) counter = 0
    
    if (thisOctal%nDepth == nDepth) then
       counter = counter + 1 
       array(counter)%content => thisOctal
       array(counter)%inUse = .true.
    endif

    IF (PRESENT(maxOctals)) THEN
      IF (counter >= maxOctals) RETURN
    END IF
    
    IF ( thisOctal%nChildren > 0 .and.(thisOctal%nDepth < nDepth)) THEN
      DO i = 1, thisOctal%nChildren, 1
        
        ! call this subroutine recursively on each of its children
        child => thisOctal%child(i)
        CALL getOctalArrayLevel(child,array,counter,nDepth)
        
      END DO
    END IF

  END SUBROUTINE getOctalArrayLevel

  
  SUBROUTINE sortOctalArray(array,grid)
    ! sorts by radius using HeapSort
    
    IMPLICIT NONE

    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    TYPE(gridType), INTENT(IN) :: grid
    
    TYPE(octal), POINTER  :: tempContent
    LOGICAL(KIND=logic), DIMENSION(8) :: tempInUse
    
    INTEGER               :: i, n

    n = SIZE(array) 
    
    DO i = n/2, 1, -1 
      CALL sift_down(i,n)
    END DO
    
    DO i = n, 2, -1
      ! swap the array elements
      tempContent      => array(1)%content
      tempInUse        =  array(1)%inUse
      array(1)%content => array(i)%content 
      array(1)%inUse   =  array(i)%inUse 
      array(i)%content => tempContent 
      array(i)%inUse   =  tempInUse 
      CALL sift_down(1,i-1)
    END DO

    CONTAINS

    SUBROUTINE sift_down(l,r)
    
      INTEGER, INTENT(IN) :: l, r
      INTEGER :: j, jold 
      real(oct) :: valueX
      real(oct) :: valueY
      real(oct) :: valueA
      real(oct) :: valueJ
      TYPE(octal), POINTER :: Acontent
      LOGICAL(KIND=logic), DIMENSION(8) :: AinUse
      TYPE(vector)    :: starPos
      
      starPos = grid%starPos1
      valueA = modulus(array(l)%content%centre - starPos)
      Acontent => array(l)%content
      AinUse   =  array(l)%inUse
      
      jold = l
      j = l + l 
      DO 
        IF (j > r) EXIT 
        
        IF (j < r) THEN
          valueX = modulus(array(j)%content%centre   - starPos)
          valueY = modulus(array(j+1)%content%centre - starPos)
          IF (valueX < valueY) j = j + 1
        END IF 
        
        valueJ = modulus(array(j)%content%centre - starPos)
        IF (valueA >= valueJ) EXIT

        array(jold)%content => array(j)%content
        array(jold)%inUse   =  array(j)%inUse
        
        jold = j
        j = j + j
      END DO
      array(jold)%content => Acontent
      array(jold)%inUse   =  AinUse

    END SUBROUTINE sift_down 

  END SUBROUTINE sortOctalArray    
    
  
  RECURSIVE SUBROUTINE finishGrid(thisOctal, grid, romData)
    ! takes the octree grid that has been created using 'splitGrid'
    !   and calculates all the other variables in the model.
    ! this should be called once the structure of the grid is complete.
    
    USE input_variables, ONLY : useHartmannTemp
    USE cluster_class, ONLY:    assign_grid_values 
    USE luc_cir3d_class, ONLY:  calc_cir3d_temperature
    USE cmfgen_class, ONLY:     calc_cmfgen_temperature
    USE jets_mod, ONLY:         calcJetsTemperature
    USE romanova_class, ONLY:   calc_romanova_temperature

    IMPLICIT NONE

    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry

    TYPE(octal), POINTER   :: child
  
    INTEGER :: subcell, iChild
    REAL T1, T2, TSTART
    logical, save :: firsttime = .true.
    integer, save :: counter = 0
     
    ! all of the work that must be done recursively goes here: 
    DO subcell = 1, thisOctal%maxChildren
   
       SELECT CASE (grid%geometry)

       CASE ("ttauri")
          IF (.NOT. useHartmannTemp) &
               CALL calcTTauriTemperature(thisOctal,subcell)
          
       CASE ("jets")
          CALL calcJetsTemperature(thisOctal,subcell, grid)
          
       CASE ("luc_cir3d")
          CALL calc_cir3d_temperature(thisOctal,subcell)
        
       CASE ("cmfgen")
         CALL calc_cmfgen_temperature(thisOctal,subcell)

       CASE ("romanova")
          CALL calc_romanova_temperature(romData, thisOctal,subcell)
        
       CASE ("cluster","wr104")
          call assign_grid_values(thisOctal,subcell)
          
       CASE ("molcluster", "theGalaxy")
          if(thisoctal%haschild(subcell)) then 
             continue
          else
             if(firsttime) then
                call CPU_TIME(tstart)
                firsttime = .false.
             endif
             CALL CPU_TIME(T1)
            
             if(thisoctal%cornervelocity(14)%x .eq. -9.9d99) then
                recentoctal => thisoctal
                CALL fillDensityCorners(thisOctal,grid,clusterdensity, clustervelocity, thisOctal%threed)
                thisOctal%velocity = thisoctal%cornervelocity(14)
             endif

             call assign_grid_values(thisOctal,subcell,grid)

             counter = counter + 1
             CALL CPU_TIME(T2)
             write(112, *) counter, t2 - t1, t2 - tstart
          endif
          
          CASE("molebench")
             CALL fillDensityCorners(thisOctal,grid, molebenchdensity, molebenchvelocity, thisOctal%threed)

       CASE DEFAULT
          ! Nothing to be done for this geometry so just return. 
          goto 666
          
       END SELECT
      
    END DO
   
    DO iChild = 1, thisOctal%nChildren, 1
       child => thisOctal%child(iChild)
       CALL finishGrid(child, grid, romData=romData)
    END DO

666 continue
    
  END SUBROUTINE finishGrid
 
  
  SUBROUTINE countVoxels(thisOctal,nOctals,nVoxels)  
    ! count the number of octals in the current section of the grid.
    ! also counts the number of unique volume elements (voxels) i.e.
    !   those subcells that are not subdivided

    IMPLICIT NONE

    TYPE(OCTAL), POINTER  :: thisOctal 
    INTEGER,INTENT(INOUT) :: nOctals   ! number of octals
    INTEGER,INTENT(INOUT) :: nVoxels   ! number of childless subcells
    
    nOctals = 0 
    nVoxels = 0
    CALL countVoxelsPrivate(thisOctal)
    
    CONTAINS
    
      RECURSIVE SUBROUTINE countVoxelsPrivate(thisOctal)
      
        TYPE(OCTAL), POINTER  :: thisOctal 
        TYPE(OCTAL), POINTER  :: child
        INTEGER :: i
        
        nOctals = nOctals + 1
        
        IF ( thisOctal%nChildren > 0 ) THEN
          ! call this subroutine recursively on each of its children
          DO i = 1, thisOctal%nChildren, 1
            child => thisOctal%child(i)
            CALL countVoxelsPrivate(child)
          END DO
        END IF

        ! increment the counter once for each of its childless subcells
        nVoxels = nVoxels + (thisOctal%maxChildren - thisOctal%nchildren)
        
      END SUBROUTINE countVoxelsPrivate

  END SUBROUTINE countVoxels


  SUBROUTINE startReturnSamples (startPoint,direction,grid,          &
             sampleFreq,nSamples,maxSamples,thin_disc_on, opaqueCore,hitCore,      &
             usePops,iLambda,error,lambda,kappaAbs,kappaSca,velocity,&
             velocityDeriv,chiLine,levelPop,rho, temperature, Ne, inflow, &
             etaCont, etaLine)
    ! samples the grid at points along the path.
    ! this should be called by the program, instead of calling 
    !   returnSamples directly, because this checks the start and finish
    !   points are within the grid bounds. returnSamples *assumes* this 
    !   criterion is met. also, the path of the photon is checked for 
    !   intersection with the stellar surface(s), disc etc. 
 
    IMPLICIT NONE

    TYPE(vector), INTENT(IN)      :: startPoint ! photon start point
    TYPE(vector), INTENT(IN)      :: direction  ! photon direction 
    TYPE(gridtype), INTENT(IN)         :: grid       ! the entire grid structure
    REAL, INTENT(IN)                   :: sampleFreq ! the maximum number of
!    real(oct), INTENT(IN)   :: sampleFreq ! the maximum number of 
                       ! samples that will be made in any subcell of the octree. 
                       
    INTEGER, INTENT(OUT)             :: nSamples   ! number of samples made
    INTEGER, INTENT(IN)                :: maxSamples ! size of sample arrays 
    logical, intent(in)                :: thin_disc_on   ! T to include thin disc
    LOGICAL, INTENT(IN)                :: opaqueCore ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: hitCore    ! true if the core is opaque
    LOGICAL, INTENT(IN)                :: usePops    ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda    ! wavelength index
    INTEGER, INTENT(INOUT)             :: error      ! error code
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda     ! path distances of sample locations

    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs   ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca   ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: chiLine    ! line opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: rho        ! density at sample points
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL:: levelPop ! level populations
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: temperature! temperature [K] at sample points
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: Ne ! electron density
    logical,DIMENSION(:),INTENT(INOUT),OPTIONAL  ::inFlow   ! flag to tell if the point is inflow
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity

    TYPE(vector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    TYPE(vector)       :: currentPoint ! current position of ray 
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(octal)             :: octree       ! the octree structure within 'grid'
    TYPE(vector)       :: directionNormalized
    real(oct)    :: distanceLimit ! max length of ray before aborting
    ! margin is the size of the region around the edge of a subcell
    !   where numerical inaccuracies may cause problems.
    real(oct)    :: margin
 
    ! variables for testing special cases (stellar intersections etc.)
    TYPE(vector)       :: starPosition       ! position vector of stellar centre
    TYPE(vector)       :: diskNormal         ! disk normal vector
    real(oct)    :: rStar              ! stellar radius
    real(oct)    :: endLength          ! max path length of photon
    TYPE(vector)       :: endPoint           ! where photon leaves grid
    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
    real(oct)    :: distanceFromOrigin ! closest distance of plane from origin
    TYPE(vector)       :: diskIntersection   ! point of photon intersection with disk
    real(oct)    :: starIntersectionDistance1 ! distance to first  intersection
    real(oct)    :: starIntersectionDistance2 ! distance to second intersection
    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
    real(oct)    :: intersectionRadius ! disk radius when photon intersects
    real(oct)    :: diskDistance       ! distance to disk intersection
    real(oct)    :: distanceThroughStar! distance of chord through star
    TYPE(vector)       :: dummyStartPoint    ! modified start point 
!    real(oct), PARAMETER :: fudgefactor = 1.00001 ! overestimates stellar size
    real(oct), PARAMETER :: fudgefactor = 1.000001 ! overestimates stellar size
    
    
    ! we will abort tracking a photon just before it reaches the edge of the
    !   simulation space. This is the fraction of the total distance to use:
    real(oct), PARAMETER :: distanceFraction = 0.999_oc 
    
    ! we will abort tracking any rays which are too close to 
    !   to a cell wall. The distance to use is defined by:
    margin = 6.0_oc * REAL(grid%maxDepth,kind=oct) * EPSILON(1.0_oc)
    ! some more experimentation is required to find out the best value for
    !   margin.

    
    ! set up some variables
    octree = grid%octreeRoot  
    abortRay = .FALSE.
    hitCore = .FALSE.
    directionNormalized = direction
    CALL normalize(directionNormalized)
    distanceLimit = HUGE(distanceLimit)
    locator = startPoint
    nSamples = 0

    currentPoint = startPoint
    
    IF (.NOT. inOctal(octree,startPoint)) THEN
      PRINT *, 'Attempting to find path between point(s) outwith the grid.'
      PRINT *, ' in amr_mod::startReturnSamples.'
      PRINT *, ' ==> StartPoint = (', startPoint, ')'
      error = -30
      return
!      startpoint%x = min(max(startpoint%x,octree%centre%x - octree%subcellSize),octree%centre%x + octree%subcellSize)
!      startpoint%y = min(max(startpoint%y,octree%centre%y - octree%subcellSize),octree%centre%y + octree%subcellSize)
!      startpoint%z = min(max(startpoint%z,octree%centre%z - octree%subcellSize),octree%centre%z + octree%subcellSize)

    ENDIF
   
   
    ! geometry-specific tests should go here
    IF (grid%geometry(1:6) == "ttauri" .OR. grid%geometry(1:4) == "jets" .or. &
         grid%geometry(1:9) == "luc_cir3d" .or. grid%geometry(1:6) == "cmfgen" .or. &
         grid%geometry(1:8) == "romanova") THEN
      
       ! need to test for both star and disc intersections 
      
       ! we will find out when and where the photon leaves the simulation space 

       call getExitPoint2(currentPoint, directionNormalized ,locator, abortRay, error, &
       grid%halfSmallestSubcell, endPoint, endLength, grid, octree, edgeofgrid=.true.)


!       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
!                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
!                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 
       endLength = endLength * distanceFraction
       
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startPoint
       
       ! for the moment, we will just assume a 2-d disc.
       absorbPhoton = .FALSE.
       diskNormal = grid%diskNormal
       starPosition = grid%starPos1
                                  
       ! find the (geometrycally thin) disk intersection point
       if (grid%geometry == "luc_cir3d" .or. grid%geometry == "cmfgen") then
          intersectionFound = .false.
       else
          if (thin_disc_on) then
             distanceFromOrigin = modulus(grid%starPos1)
             diskIntersection = intersectionLinePlane(startPoint, directionNormalized,&
                  diskNormal, distanceFromOrigin, intersectionFound)
          else
             intersectionFound = .false.
          end if
       end if

       
       IF (intersectionFound) THEN
       
         ! we need to check whether the photon passes through the disk's
         !   central hole, or the outside of the outer radius of accretion disc.
         intersectionRadius =  modulus(diskIntersection - starPosition)
         IF (intersectionRadius > grid%diskRadius .and.  &
              intersectionRadius < 1.5d5) THEN  ! assuming 100 AU disc size
!              intersectionRadius < grid%octreeRoot%subcellsize*2.0) THEN
           absorbPhoton = .TRUE.

           ! we need to check whether the intersection occurs within our
           !   simulation space.
           diskDistance = modulus(diskIntersection-startPoint)
           IF (diskDistance < endLength) THEN
           
             endLength = diskDistance
             
           END IF
         END IF
       END IF
         
         
       ! now we check for intersections with the star

       rStar = grid%rStar1 
       CALL intersectionLineSphere(startPoint,directionNormalized,endLength,starPosition, &
                                   rStar,starIntersectionFound1,starIntersectionFound2,   &
                                   starIntersectionDistance1,starIntersectionDistance2)
       ! by passing a line segment to intersectionLineSphere, we ensure that we
       !   do not find intersections with the star that take place after the 
       !   photon has been absorbed by the disk.
       IF (starIntersectionFound1) THEN
       
         endLength = starIntersectionDistance1
         IF (opaqueCore) absorbPhoton = .TRUE.
         
       END IF

       ! we trace the photon path until it encounters the disk, the star,
       !   or the edge of the simulation space.
       error = 0
       CALL returnSamples(currentPoint,startPoint,locator,directionNormalized,&
                     octree,grid,sampleFreq,nSamples,maxSamples,abortRay,     &
                     lambda,usePops,iLambda,error,margin,endLength,           &
                     kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,   &
                     velocityDeriv=velocityDeriv,chiLine=chiLine,             &
                     levelPop=levelPop,rho=rho, temperature=temperature,      &
                     Ne=Ne, inFlow=inFlow, etaCont=etaCont, etaLine=etaLine)
      
       IF (error < 0)  RETURN

       IF (.NOT. opaqueCore .AND. starIntersectionFound2) THEN 
         ! the photon passes through the star and continues on the other side.

         ! we have to adjust the arguments to returnSamples to make the 
         !   output arrays correct.
         distanceThroughStar = fudgeFactor * &
                      (starIntersectionDistance2 - starIntersectionDistance1)
         currentPoint = currentPoint + distanceThroughStar * directionNormalized
         locator = currentPoint
         dummystartPoint = startPoint + distanceThroughStar * directionNormalized
         distanceLimit = endLength - distanceThroughStar
           
        CALL returnSamples(currentPoint,dummyStartPoint,locator,             &
                    directionNormalized,octree,grid,sampleFreq,nSamples,     &
                    maxSamples,abortRay,lambda,usePops,iLambda,error,margin, &
                    distanceLimit,kappaAbs=kappaAbs,kappaSca=kappaSca,       &
                    velocity=velocity,velocityDeriv=velocityDeriv,           &
                    chiLine=chiLine,levelPop=levelPop,rho=rho,               &
                    temperature=temperature, Ne=Ne, inFlow=inFlow,           &
                    etaCont=etaCont, etaLine=etaLine)

       END IF

       ! if the photon ends up in the disk or the star, we make sure it absorbed.
       IF (absorbPhoton) THEN
       
         nSamples = nSamples + 1
         IF (nSamples > maxSamples) THEN
           PRINT *, "Error:: nSamples > maxSamples in startReturnSamples subroutine"
           PRINT *, "        nSamples   = ", nSamples
           PRINT *, "        maxSamples = ", maxSamples
           STOP
         END IF
         
         lambda(nSamples) = lambda(nSamples-1)
         hitCore = .TRUE.
         IF (PRESENT(kappaSca))      kappaSca(nSamples) = 0.
         IF (PRESENT(kappaAbs))      kappaAbs(nSamples) = 1.e20
         IF (PRESENT(velocity))      velocity(nSamples) = vector(0.,0.,0.)
         IF (PRESENT(velocityDeriv)) velocityDeriv(nSamples) = 1.0
         IF (PRESENT(chiLine))       chiLine(nSamples) = 1.e20
         IF (PRESENT(levelPop))      levelPop(nSamples,:) = 0.0
         IF (PRESENT(Ne))            Ne(nSamples) = 0.0
         IF (PRESENT(rho))           rho(nSamples) = 0.0
         IF (PRESENT(temperature))   temperature(nSamples) = 0.0
         IF (PRESENT(inFlow))        inFlow(nSamples) = .true.
         IF (PRESENT(etaCont))       etaCont(nSamples) = 1.e-20
         IF (PRESENT(etaLine))       etaLine(nSamples) = 1.e-20
       END IF
       
    ELSE
            

       call getExitPoint2(currentPoint, directionNormalized ,locator, abortRay, error, &
       grid%halfSmallestSubcell, endPoint, endLength, grid, octree, edgeofgrid=.true.)


!       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
!                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
!                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 
       distanceLimit = endLength * distanceFraction
       
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startPoint
       
       CALL returnSamples(currentPoint,startPoint,locator,directionNormalized,&
                   octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,&
                   usePops,iLambda,error,margin,distanceLimit,                &
                   kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,     &
                   velocityDeriv=velocityDeriv,chiLine=chiLine,               &
                   levelPop=levelPop,rho=rho,temperature=temperature,         &
                   Ne=Ne,inFlow=inFlow,etaCont=etaCont,etaLine=etaLine)
    END IF
      
  END SUBROUTINE startReturnSamples

  
  RECURSIVE SUBROUTINE returnSamples (currentPoint,startPoint,locator,direction, &
             octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,usePops, &
             iLambda,error,margin,distanceLimit,kappaAbs,kappaSca,velocity,      &
             velocityDeriv,chiLine,levelPop,rho,temperature, Ne, inFlow,         &
             etaCont,etaLine)
    ! this uses a recursive ray traversal algorithm to sample the octal
    !   grid at points along the path of the ray. 
    ! no checks are made that the ray lies within the boundaries of the
    !  grid, so this subroutine not be called directly. use  
    !  'startReturnSamples' instead.

    IMPLICIT NONE

    ! note: in the code comments, the terms 'subcell' and 'cell' are
    !  mostly interchangeable.

    TYPE(vector), INTENT(INOUT)    :: currentPoint ! current ray position
    TYPE(vector), INTENT(IN)       :: startPoint   ! initial ray position
    TYPE(vector)                   :: locator, rotloc
                  ! 'locator' is used to indicate a point that lies within the  
                  !   *next* cell of the octree that the ray will interesect.
                  !   initially this will be the same as the currentPoint
      
    TYPE(vector), INTENT(IN)       :: direction    ! ray's direction vector 
    TYPE(octal), INTENT(IN)             :: octree       ! an octal grid
    TYPE(gridtype), INTENT(IN)          :: grid         ! grid structure
    REAL, INTENT(IN)                    :: sampleFreq
!    real(oct), INTENT(IN)    :: sampleFreq
                  ! 'sampleFreq' is the maximum number of samples that will be made
                  !   in any subcell of the octree. 
    INTEGER, INTENT(INOUT)              :: nSamples     ! number of samples made
    INTEGER, INTENT(IN)                 :: maxSamples   ! size of sample arrays 
    LOGICAL, INTENT(INOUT)              :: abortRay     ! flag to signal completion
    REAL, DIMENSION(:), INTENT(INOUT)   :: lambda       ! distance travelled by photon
    INTEGER, INTENT(IN)                 :: iLambda      ! wavelength index
    INTEGER, INTENT(INOUT)              :: error        ! error code
    real(oct)                :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    real(oct), INTENT(IN)    :: distanceLimit ! max length of ray before aborting
    LOGICAL, INTENT(IN)                 :: usePops      ! whether to use level populations
    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs     ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca     ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL :: chiLine       ! line opacities
    REAL(double),DIMENSION(:),OPTIONAL               :: rho           ! density at sample points
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL   :: temperature! in [K]
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL   :: Ne         ! electron density
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL :: levelPop   ! level populations 
    logical, DIMENSION(:),INTENT(INOUT),OPTIONAL       :: inFlow     ! inflow flag
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity



    TYPE(vector)      :: exitPoint      ! where ray leaves current cell
    TYPE(vector)      :: centre         ! centre of current subcell
    real(oct)   :: subcellSize    ! size of current subcell
    real(oct)   :: minWallDistance ! distance to *nearest* wall
    
    real(oct)   :: length         ! distance from start of the ray's path.
    real(oct)   :: sampleLength   ! distance interval between samples 
    real(oct)   :: trialLength    ! trial distance to a possible next sample
    TYPE(vector)      :: trialPoint     ! trial location for a next sample 
    INTEGER                :: trial          ! loop counter for trial points 
    INTEGER                :: subcell        ! current subcell 
    INTEGER                :: subIndex       ! octal's child index for current subcell 
    INTEGER                :: i
    real(oct) :: dL


    ! find which of the subcells the point lies in
    DO 
       if (octree%twod) then
          rotloc = projecttoxz(locator)
       else
          rotloc = locator
       endif
       if (octree%oneD) then
          rotloc  = VECTOR(modulus(locator),0.d0,0.d0)
       endif



      subcell = whichSubcell(octree,rotloc)

      ! the IF statement below is for debugging
      !IF ( looseInOctal(octree,locator) .EQV. .FALSE. ) THEN
      !  PRINT *, "Panic: looseInOctal failed in returnSamples"
      !  STOP
      !ENDIF
 
      ! if the subcell has a child, we will find a pointer to that
      !   child and use recursion to sample it.

      IF (octree%hasChild(subcell)) THEN

        ! find pointer to child
        subindex = -99
        DO i = 1, octree%nChildren, 1
          IF ( octree%indexChild(i) == subcell ) THEN
            subIndex = i  
            EXIT
          ENDIF
        ENDDO
        IF (subindex == -99) THEN
          PRINT *, ' Panic: subindex not found'
          STOP
        ENDIF

        CALL returnSamples(currentPoint,startPoint,locator,direction,        &
                    octree%child(subIndex),grid,sampleFreq,nSamples,         &
                    maxSamples,abortRay,lambda,usePops,iLambda,error,margin, &
                    distanceLimit,kappaAbs=kappaAbs,kappaSca=kappaSca,       &
                    velocity=velocity,velocityDeriv=velocityDeriv,           &
                    chiLine=chiLine,levelPop=levelPop,rho=rho,               &
                    temperature=temperature, Ne=Ne, inFlow=inFlow,           &
                    etaCont=etaCont,etaLine=etaLine)
        
        ! after returning from the recursive subroutine, we may have 
        !   finished tracing the ray's path
        IF ( abortRay .EQV. .TRUE. ) RETURN

        ! otherwise, we check whether the ray is still within the
        !   boundaries of the current subsection of the octree.
        IF (.NOT. inOctal(octree,locator)) RETURN
      
      
      ELSE
        ! we now consider the case of the current subcell being childless

        ! find the subcell centre
        centre = subcellCentre(octree,subcell)
        subcellSize = octree%subcellSize

        ! we find the exit point from the current subcell, 
        !  and also 'locator' - a point that lies in the *next* subcell 
        CALL getExitPoint(currentPoint,direction,locator,abortRay,error, &
                          grid%halfSmallestSubcell,exitPoint,centre,subcellSize,&
                          minWallDistance,margin,grid%octreeRoot%threed)

!       call getExitPoint2(currentPoint, direction ,locator, abortRay, error, &
!       grid%halfSmallestSubcell, exitPoint, minWallDistance, grid, octree)

        IF ( abortRay .EQV. .TRUE. ) RETURN

        ! we now decide where we are going to sample the quantities
        ! we will define an approximate rate that is a fraction of
        !   the current subcell size
        sampleLength = octree%subcellSize / sampleFreq

        ! if there are no previously sampled points, we definitely have to take a
        !   sample here
!        IF ( nSamples == 0 ) then
!          length = modulus(currentPoint - startPoint)
!          CALL takeSample(currentPoint,length,direction,grid,octree,subcell,    &
!                          nSamples,maxSamples,usePops,iLambda,error,lambda,     &
!                          kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,&
!                          velocityDeriv=velocityDeriv,chiLine=chiLine,          &
!                          levelPop=levelPop,rho=rho, Ne=Ne)
!        ENDIF

        
        ! we check whether we should take a sample at currentPoint
        !   (which will usually be the entry point to the subcell).

        length = modulus(currentPoint - startPoint)

! force sampling here        
        !IF ( modulus(currentPoint - (startPoint + (direction *         & 
        !                      REAL(lambda(nSamples),kind=oct)))) &
        !                                       > sampleLength ) THEN
          CALL takeSample(currentPoint,length,direction,grid,octree,subcell,    &
                          nSamples,maxSamples,usePops,iLambda,error,lambda,     &
                          kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,&
                          velocityDeriv=velocityDeriv,chiLine=chiLine,          &
                          levelPop=levelPop,rho=rho,temperature=temperature,    &
                          Ne=Ne, inFlow=inFlow, etaCont=etaCont, etaLine=etaLine)
        !END IF

        ! we add sampleLength to the distance from the last location
        !   that was sampled, and decide whether to take a new sample.
        trial = 1
        DO 
          trialPoint = currentPoint + direction * &
                       ((REAL(lambda(nSamples),kind=oct) + &
                       (REAL(trial,kind=oct) * sampleLength)) - length)
                       
          ! we only want to take a sample if we are still within the subcell
          dL = modulus(trialPoint - currentPoint)
          IF (  dL < minWallDistance ) THEN
          
            trialLength = length + dL
            IF (trialLength > distanceLimit) THEN
              abortRay = .TRUE.
              RETURN
            END IF
              
            IF (trialLength >= 0.0_oc) THEN
              CALL takeSample(trialPoint,trialLength,direction,grid,octree,subcell, &
                         nSamples,maxSamples,usePops,iLambda,error,lambda,     &
                         kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,&
                         velocityDeriv=velocityDeriv,chiLine=chiLine,          &
                         levelPop=levelPop,rho=rho,temperature=temperature,    &
                         Ne=Ne, inFlow=inFlow, etaCont=etaCont, etaLine=etaLine)
            ELSE
              EXIT
            END IF
            
          ELSE 
            EXIT
          END IF
          trial = trial + 1
        END DO    

        ! adjust some variables
        length = modulus(exitPoint - startPoint)
        currentPoint = exitPoint
        
      ! if we have left the boundaries of the simulation space, we are finished
      IF (.NOT. inOctal(octree,locator) .OR. &
               length > distanceLimit ) RETURN 
      
    ENDIF
    
  ENDDO


  END SUBROUTINE returnSamples


  subroutine getExitPoint2(currentPoint, direction ,locator, abortRay, error, &
       halfSmallestSubcell, exitPoint, minWallDistance, grid, thisOctal, edgeOfGrid)


    implicit none
    type(GRIDTYPE) :: grid
    type(VECTOR), intent(in) :: currentPoint
    type(VECTOR), intent(in) :: direction
    type(VECTOR), intent(out) :: locator    
    logical, intent(inout) :: abortRay
    type(OCTAL), target  :: thisOctal
    integer, intent(inout) :: error
    real(oct), INTENT(IN)    :: halfSmallestSubcell
    type(VECTOR), intent(out) :: exitPoint
    real(oct), intent(out) :: minWallDistance
    type(OCTAL), pointer :: sOctal
    logical, optional :: edgeOfGrid
    
    error = 0
    abortRay = .false.

    if (present(edgeofGrid)) then
       call distanceToGridEdge(grid, currentPoint, direction, minWallDistance)
    else
       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPoint, direction, minWallDistance, sOctal)
    endif

    exitPoint = currentPoint + minWallDistance * direction
    locator = currentPoint + (minWallDistance + 0.001d0*halfSmallestSubcell) * direction

  end subroutine getExitPoint2

  
  SUBROUTINE getExitPoint(currentPoint,direction,locator,abortRay,error,    &
                          halfSmallestSubcell,exitPoint,centre,subcellSize, &
                          minWallDistance,margin,threed)
    ! this subroutine finds the closest subcell wall in the direction the
    !   photon is travelling. 

    IMPLICIT NONE
          
    TYPE(vector), INTENT(IN)       :: currentPoint ! current ray position
    TYPE(vector), INTENT(IN)       :: direction    ! ray's direction vector 
    TYPE(vector), INTENT(OUT)      :: locator       
                  ! 'locator' is used to indicate a point that lies within the  
                  !   *next* cell of the octree that the ray will interesect.
    LOGICAL, INTENT(INOUT)              :: abortRay     ! flag to signal completion
    INTEGER, INTENT(INOUT)              :: error        ! error code
    real(oct), INTENT(IN)    :: halfSmallestSubcell
    real(oct), INTENT(IN)    :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    TYPE(vector), INTENT(OUT)      :: exitPoint    ! where ray leaves current cell
    TYPE(vector), INTENT(IN)       :: centre       ! centre of current subcell
    real(oct), INTENT(IN)         :: subcellSize  ! size of current subcell
    real(oct), INTENT(OUT)   :: minWallDistance ! distance to *nearest* wall
    
    
    real(oct)   :: wallDistanceX  ! distance to next x-wall-plane 
    real(oct)   :: wallDistanceY  ! distance to next y-wall-plane
    real(oct)   :: wallDistanceZ  ! distance to next z-wall-plane
                  ! 'wallDistanceX,Y,Z' are the distances between the current 
                  !   position and the intersections with the cell walls (along the 
                  !   'direction' vector)
    LOGICAL                :: found          ! status flag  
    real(oct)   :: wallFromOrigin ! distance of cell wall from the origin
    TYPE(vector)      :: wallNormal     ! normal to plane of cell wall

    INTEGER, parameter :: max_num_err = 10;
    INTEGER, save :: num_err = 0    

    LOGICAL :: threed
    REAL(oct) :: r2, d, cosMu, disttor2, disttoxboundary,disttozboundary
    REAL(oct) :: compZ, currentZ, tval, x1, x2
    REAL(oct) :: r1, theta, mu, disttor1

    TYPE(vector) :: xDir, zDir
    logical :: ok
    ! Specify the ratio of extra length to give it for "locater" to the 
    ! "halfSmallestSubcell" size.
!    REAL(oct), parameter :: frac =1.e-6_oc  
    REAL(oct), parameter :: frac =1.e-2_oc  

    if (threed) then
       
    ! this code is ugly, but should work OK.   
                  

    ! there are six subcell walls - and each wall is part of a plane    
    !   parallel with an axis. we use the direction of the photon to 
    !   find the three planes that the photon will intersect.
    
    IF ( direction%x > 0.0_oc ) THEN
      walldistanceX = (centre%x - currentPoint%x + subcellSize / 2.0_oc  ) / ABS(direction%x) 
    ELSE IF ( direction%x < 0.0_oc ) THEN
      wallDistanceX = (currentPoint%x - centre%x + subcellSize / 2.0_oc ) / ABS(direction%x) 
    ELSE 
      wallDistanceX = HUGE(wallDistanceX) 
    END IF

    IF ( direction%y > 0.0_oc ) THEN
      wallDistanceY = (centre%y - currentPoint%y + subcellSize / 2.0_oc ) / ABS(direction%y) 
    ELSE IF ( direction%y < 0.0_oc ) THEN
      wallDistanceY = (currentPoint%y - centre%y + subcellSize / 2.0_oc ) / ABS(direction%y) 
    ELSE 
      wallDistanceY = HUGE(wallDistanceX) 
    END IF
        
    IF ( direction%z > 0.0_oc ) THEN
      wallDistanceZ = (centre%z - currentPoint%z + subcellSize / 2.0_oc ) / ABS(direction%z) 
    ELSE IF ( direction%z < 0.0_oc ) THEN
      wallDistanceZ = (currentPoint%z - centre%z + subcellSize / 2.0_oc ) / ABS(direction%z) 
    ELSE 
      wallDistanceZ = HUGE(wallDistanceX) 
    END IF

    minWallDistance = MIN(wallDistanceX,wallDistanceY,wallDistanceZ)

    ! we may have problems if "currentPoint" actually lies slightly outside
    !  the current subcell (probably due to numerical accuracy problems). 

    ! if any of the wallDistances are negative or very small, we will 
    !  abandon this photon.

    IF ( wallDistanceX < margin .OR. &
         wallDistanceY < margin .OR. &
         wallDistanceZ < margin      ) THEN 
      abortRay = .TRUE.
      if (num_err < max_num_err) then
         PRINT *, 'In getExitPoint, aborting ray because distance is less than ''margin'''
         PRINT *, '  wallDistances=',wallDistanceX,wallDistanceY,wallDistanceZ
      elseif (num_err ==  max_num_err) then
         PRINT *, ' '
         PRINT *, 'Surpressing the further message from subroutine getExitPoint....'
         PRINT *, ' '
      else
         continue
      end if
      
      error = -20
      num_err = num_err + 1
      
      RETURN      
    END IF
       
       
    IF ( wallDistanceX < wallDistanceY .AND. &
         wallDistanceX < wallDistanceZ ) THEN
      wallNormal =  xHat
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          PRINT *, 'wallDistances = ',wallDistanceX,&
          wallDistanceY,wallDistanceZ, margin 
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHat)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*xHat)
      ENDIF

    ELSEIF ( wallDistanceY < wallDistanceX .AND. &
             wallDistanceY < wallDistanceZ ) THEN
      wallNormal =  yHat
      IF ( direction%y > 0.0_oc ) THEN
        wallFromOrigin = centre%y + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*yHat)
      ELSE
        wallFromOrigin = centre%y - (subcellSize / 2.0_oc) 
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*yHat)
      ENDIF
      
    ELSEIF ( wallDistanceZ < wallDistanceX .AND. &
             wallDistanceZ < wallDistanceY ) THEN
      wallNormal =  zHat
      IF ( direction%z > 0.0_oc ) THEN
       wallFromOrigin = centre%z + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*zHat)
      ELSE
        wallFromOrigin = centre%z - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*zHat)
      ENDIF

    ! we now consider the case where the ray is leaving through one
    !   of the vertices of the cell.
    ! we could integrate this into the code above, but it is very
    !   unlikely that it will be executed, so we avoid evaluating a 
    !   few IFs by keeping it separate.
    ELSEIF ( wallDistanceX == wallDistanceY .AND. &
             wallDistanceX /= wallDistanceZ ) THEN
      wallNormal =  xHat
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHat)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*xHat)
      ENDIF
      IF ( direction%y > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*yHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*yHat)
      ENDIF  
    ELSEIF ( wallDistanceX == wallDistanceZ .AND. &
             wallDistanceX /= wallDistanceY ) THEN
      wallNormal =  xHat
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHat)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
      IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*xHat)
      ENDIF
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHat)
      ENDIF  

    ELSEIF ( wallDistanceY == wallDistanceZ  .AND. &
             wallDistanceX /= wallDistanceZ ) THEN
      wallNormal =  yHat
      IF ( direction%y > 0.0_oc ) THEN
        wallFromOrigin = centre%y + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*yHat)
      ELSE
        wallFromOrigin = centre%y - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*yHat)
      ENDIF
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHat)
      ENDIF  
      
    ! we now consider the case where the ray is leaving through one
    !   of the corners of the cell.
         
    ELSEIF ( wallDistanceX == wallDistanceY  .AND. &
             wallDistanceY == wallDistanceZ ) THEN
      wallNormal =  xHat
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                   (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHat)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*xHat)
      ENDIF
      IF ( direction%y > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*yHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*yHat)
      ENDIF  
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHat)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHat)
      ENDIF  

    ! we should now have run out of possibilities
    ELSE
      PRINT *, 'Panic: Fell through direction finder!'
      stop
!      DO ; END DO ; STOP ! sit in loop for debugging purposes
    END IF

    else       

       ! two-d case written by TJH on 27/1/05
       ! this code isn't ugly and it does work !!!!

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

       r1 = centre%x - subcellSize/2.d0
       r2 = centre%x + subcellSize/2.d0
       d = sqrt(currentpoint%x**2+currentpoint%y**2)
       xDir = VECTOR(currentpoint%x, currentpoint%y,0.d0)


       if ((direction%x**2+direction%y**2) /= 0.0) then
          if (modulus(xDir) /= 0.d0) then
             call normalize(xDir)       
             cosmu =((-1.d0)*xDir).dot.direction
             call solveQuadDble(1.d0, -2.d0*d*cosmu, d*d-r2*r2, x1, x2, ok)
             if (.not.ok) then
                write(*,*) "Error:: Quad solver failed in [amr_mod::getExitPoint] (1)"
                write(*,*) " Correctionn forced for now. If problem insistet, you need to fix it."
                write(*,*) "           d = ", d
                write(*,*) "      cosmud = ", cosmu
                write(*,*) "          r2 = ", r2
                write(*,*) "currentPoint = ", currentPoint
                x1 = 1.0d-5; x2=1.0d-5
!                stop
!                do;enddo
             endif
             distTor2 = max(x1,x2)      
         
         
             theta = asin(max(-1.d0,min(1.d0,r1 / d)))
             cosmu = xDir.dot.direction
             mu = acos(max(-1.d0,min(1.d0,cosmu)))
             distTor1 = 1.e30
             if (mu  < theta ) then
                call solveQuadDble(1.d0, -2.d0*d*cosmu, d*d-r1*r1, x1, x2, ok)
                if (.not.ok) then
                   write(*,*) "Quad solver failed in  [amr_mod::getExitPoint] (2)"
                   write(*,*) "           d = ", d
                   write(*,*) "      cosmud = ", cosmu
                   write(*,*) "          r2 = ", r2
                   write(*,*) "currentPoint = ", currentPoint
                   stop
!                   do;enddo
                endif
             endif
             distTor1 = max(x1,x2)
             distToXboundary = min(distTor1, distTor2)
          else
             distToXboundary = 1.e30
          endif
       else
          distToXboundary = 1.e30
       endif

       zDir = VECTOR(0.d0, 0.d0, 1.d0)
       compZ = zDir.dot.direction
       currentZ = currentpoint%z

       if (compZ /= 0.d0 ) then
          if (compZ > 0.d0) then
             distToZboundary = (centre%z + subcellsize/2.d0 - currentZ ) / compZ
          else
             distToZboundary = abs((centre%z - subcellsize/2.d0 - currentZ ) / compZ)
          endif
       else
          disttoZboundary = 1.e30
       endif

   
       tVal = min(distToZboundary, distToXboundary)
       if (tval <= 0.) tval = 1.0e-5
       if (tVal > 1.e29) then
          write(*,*) "Error :: tVal > 1.e29 [amr_mod:getExitPoint]."
          write(*,*) "tVal,compZ, distToZboundary,disttoxboundary = "
          write(*,*) tVal,compZ, distToZboundary,disttoxboundary
          write(*,*) "z = ", currentZ
          stop 
       elseif (tval <= 0.) then
          write(*,*) "Error :: tVal <= 0  [amr_mod:getExitPoint]."
          write(*,*) "tVal, compZ, distToZboundary,disttoxboundary = "
          write(*,*) tVal,compZ, distToZboundary,disttoxboundary
          write(*,*) "z = ", currentZ
          stop
       endif
   
       minWallDistance = tVal
   
       exitPoint = currentPoint + tval * direction
       locator = exitPoint + frac*halfSmallestSubcell * direction

!       write(*,*) tval, exitpoint

    endif

  END SUBROUTINE getExitPoint 



  SUBROUTINE takeSample(point,length,direction,grid,thisOctal,subcell,nSamples,&
                        maxSamples,usePops,iLambda,error,lambda,kappaAbs,      &
                        kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,  &
                        temperature,Ne,inFlow,etaCont,etaLine) 
  
    use jets_mod, only: JetsVelocity, get_jets_parameter, dV_dn_jets    

    IMPLICIT NONE
    
    TYPE(vector), INTENT(IN)      :: point     ! place to make sample
    real(oct), INTENT(IN)   :: length    ! 
    TYPE(vector), INTENT(IN)      :: direction ! direction vector
    TYPE(gridtype), INTENT(IN)         :: grid      ! grid structure
    TYPE(octal), TARGET, INTENT(IN)    :: thisOctal ! grid structure
    INTEGER, INTENT(IN)                :: subcell   ! subcell containing 'point'
    INTEGER, INTENT(INOUT)             :: nSamples  ! number of samples so far
    INTEGER, INTENT(IN)                :: maxSamples! size of sample arrays 
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda    ! distances of samples from startPoint 
    LOGICAL, INTENT(IN)                :: usePops   ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda   ! wavelength index
    INTEGER, INTENT(INOUT)             :: error     ! error code
    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL          :: kappaAbs  ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL          :: kappaSca  ! scattering opacities 
    TYPE(vector), DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: velocityDeriv   ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: chiLine   ! line opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL          :: rho       ! density at sample points
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL :: levelPop  ! level populations
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL   :: Ne        ! electron density
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL           :: temperature ! in [K]
    logical,DIMENSION(:),INTENT(INOUT),OPTIONAL        :: inFlow     ! indicates if the cell is in use.
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity


    TYPE(vector)                       :: directionReal ! direction vector (REAL values)
    TYPE(octal), POINTER               :: localPointer  ! pointer to the current octal
 
    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    REAL :: Vr, r, Rs
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn)

    nSamples = nSamples + 1
    IF (nSamples > maxSamples) THEN
      PRINT *, "Erro:: nSamples > maxSamples in takeSample subroutine"
      PRINT *, "nSamples   = ", nSamples
      PRINT *, "maxSamples = ", maxSamples      
      write(*,*) "lambda(nSamples-10:nSamples) = ", lambda(nSamples-10:nSamples)
      STOP
    END IF

    
    lambda(nSamples) = length
    directionReal = direction

    localPointer => thisOctal

    IF (usePops) THEN
      CALL amrGridValues(grid%octreeRoot,point,startOctal=localPointer,&
                         actualSubcell=subcell,                        &
                         iLambda=iLambda,                              &
                         direction=directionReal,                      &
                         velocity=velocity(nSamples),                  &
                         velocityDeriv=velocityDeriv(nSamples),        &
                         kappaAbs=kappaAbs(nSamples),                  &
                         kappaSca=kappaSca(nSamples),                  &
                         rho=rho(nSamples),                            &
                         N=levelPop(nSamples,:),                       &
                         grid=grid,                                    &
                         temperature=temperature(nSamples),            &
                         Ne=Ne(nSamples),                              &
                         inFlow=inFlow(nSamples)                       &
                         )
    ELSE
      CALL amrGridValues(grid%octreeRoot,point,startOctal=localPointer,&
                         actualSubcell=subcell,                        &
                         iLambda=iLambda,                              &
                         direction=directionReal,                      &
                         velocity=velocity(nSamples),                  &
                         velocityDeriv=velocityDeriv(nSamples),        &
                         kappaAbs=kappaAbs(nSamples),                  &
                         kappaSca=kappaSca(nSamples),                  &
                         rho=rho(nSamples),                            &
                         chiLine=chiLine(nSamples),                    &
                         grid=grid,                                    &
                         temperature=temperature(nSamples),            &
                         Ne=Ne(nSamples),                              &
                         inFlow=inFlow(nSamples),                      &
                         etaCont=etaCont(nSamples),                    &
                         etaLine=etaLine(nSamples)                     &                         
                         )
    END IF

    
    ! special case if the geometry is "jets" (for accuracy)
    If (grid%geometry(1:4) == "jets" ) then
       ! using routines in jets_mod.f90
       rho(nSamples) = Density(pointVec, grid)                 ! in [g/cm^3]
       Vr = JetsVelocity(pointVec, grid)/(cSpeed/1.0e5)        ! in [c]
       r = modulus(pointVec)
       Rs = get_jets_parameter("Rmin")
       if (r < Rs) r =Rs
       velocity(nSamples) = REAL(Vr, kind=oct)*(pointVec/REAL(r,kind=oct))
       velocityDeriv(nSamples) = dV_dn_jets(VECTOR(pointVec%x, pointVec%y, pointVec%z), &
            VECTOR(direction%x, direction%y, direction%z))    ! [1/sec]
    end If

    ! use variables to silence compiler warnings
    error = error

  END SUBROUTINE takeSample

  
  SUBROUTINE amrGridValues(octalTree,point,startOctal,foundOctal,&
                           foundSubcell,actualSubcell,iLambda,lambda,direction,&
                           velocity,velocityDeriv,temperature,kappaAbs,&
                           kappaSca,rho,chiLine,etaLine,etaCont, &
                           probDistLine,probDistCont,N,Ne,nTot,inflow,grid, &
                           interp, departCoeff,kappaAbsArray,kappaScaArray, dusttypeFraction, rosselandKappa, kappap, &
                           atthistemperature)

    USE input_variables, only: nLambda, hydrodynamics

    ! POINT, direction --> should be in unrotated coordinates for 2D case (not projected onto x-z plane!)
    !

    ! returns one or more physical variables at a given point in the grid.
    ! optional arguments should be specified for all of the variables that
    !   are wanted.
    
    ! if the foundOctal argument is supplied, it is made to point to 
    !   the octal containing 'point'.
    ! if the foundSubcell argument is supplied, it is made to point to 
    !   the subcell containing 'point'.
    ! if the startOctal argument is supplied, the function uses a 
    !   local search for the correct octal starting at that octal.
    ! if actualSubcell and startSubcell are both supplied, these 
    !   locations are assumed to be correct and no search is performed.

    IMPLICIT NONE

    TYPE(octal), POINTER              :: octalTree
    TYPE(vector), INTENT(IN)     :: point
    TYPE(vector)                 :: point2 ! may be projected point
    TYPE(octal), OPTIONAL, POINTER    :: startOctal
    TYPE(octal), OPTIONAL, POINTER    :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL    :: foundSubcell
    INTEGER, INTENT(IN), OPTIONAL     :: actualSubcell 
    INTEGER, INTENT(IN), OPTIONAL     :: iLambda       ! wavelength index
    TYPE(vector),INTENT(IN),OPTIONAL  :: direction     
    TYPE(gridtype),INTENT(IN),OPTIONAL:: grid          
    LOGICAL, INTENT(IN), OPTIONAL     :: interp        ! use interpolation
                                      !  ^^^^^^ ! not implemented yet???
    REAL,INTENT(IN),OPTIONAL          :: lambda
    TYPE(vector),INTENT(OUT),OPTIONAL :: velocity
    REAL,INTENT(OUT),OPTIONAL         :: velocityDeriv
    REAL,INTENT(OUT),OPTIONAL         :: temperature
    REAL(double),OPTIONAL         :: kappaAbs
    REAL(double),OPTIONAL         :: kappaSca
    REAL(double),OPTIONAL         :: kappaAbsArray(nLambda)
    REAL(double),OPTIONAL         :: kappaScaArray(nLambda)
    REAL(double),OPTIONAL        :: rosselandKappa
    REAL, OPTIONAL        :: kappap
    REAL,INTENT(IN), OPTIONAL         :: atThisTemperature
    REAL(double),INTENT(OUT),OPTIONAL         :: rho
    REAL,INTENT(OUT),OPTIONAL         :: chiLine
    REAL(double),INTENT(OUT),OPTIONAL         :: etaLine
    REAL(double),INTENT(OUT),OPTIONAL         :: etaCont
    real(double), dimension(:), intent(out), optional :: dusttypeFraction
    real(double),INTENT(OUT),OPTIONAL :: probDistLine
    real(double),INTENT(OUT),OPTIONAL :: probDistCont
    real(double),DIMENSION(:),INTENT(OUT),OPTIONAL :: N
    real(double),INTENT(OUT),OPTIONAL :: Ne
    real(double),INTENT(OUT),OPTIONAL :: nTot
    REAL,DIMENSION(:),INTENT(OUT),OPTIONAL     :: departCoeff
    LOGICAL,INTENT(OUT),optional      :: inFlow
    
    TYPE(octal), POINTER              :: resultOctal
    INTEGER                           :: subcell
    LOGICAL                           :: interpolate
    LOGICAL                           :: boundaryProblem = .false.

!    type(vector), save :: lastpoint, lastpoint2

! this for possibility of twoD AMR grid


    point2 = point
!    if(point .eq. lastpoint) then
!       point2 = lastpoint2
!    else
       if (octaltree%twoD) then
          if (.not.hydrodynamics) then
             point2 = projectToXZ(point)
          else
             point2 = point
          endif
       elseif(octaltree%oneD) then
          point2  = VECTOR(modulus(point),0.d0,0.d0)
       endif
!    endif
!    lastpoint = point
!    lastpoint2 = point2

    IF (PRESENT(interp)) THEN
      interpolate = interp
    ELSE
      interpolate = .FALSE. 
    END IF
    
    ! first we find the correct subcell (if we need to)
    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
    ! called with rotated position if 2d octal
        CALL findSubcellLocal(point2,startOctal,subcell, boundaryProblem)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
        if (boundaryProblem) then
          boundaryProblem = .false.
          CALL findSubcellTD(point2,octalTree,resultOctal,subcell)
          IF (PRESENT(foundOctal))   foundOctal   => resultOctal
          IF (PRESENT(foundSubcell)) foundSubcell =  subcell
        end if
      END IF
      resultOctal => startOctal
      
    ELSE
      CALL findSubcellTD(point2,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal))   foundOctal   => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

    END IF
 
    IF (interpolate) THEN
       PRINT *, 'Interpolation not implemented!' 
       STOP
    ELSE


!      IF (PRESENT(velocity))         velocity = resultOctal%velocity(subcell)

       ! unrotated poistion should be passed here!

       IF (PRESENT(velocity))  velocity = amrGridVelocity(octalTree,point,startOctal=resultOctal,&
            actualSubcell=subcell) 



      IF (PRESENT(temperature))   temperature = resultOctal%temperature(subcell)
      IF (PRESENT(rho))                   rho = resultOctal%rho(subcell)
      IF (PRESENT(chiLine))           chiLine = resultOctal%chiLine(subcell)

!      IF (PRESENT(chiLine)) then
!         CALL interpAMR(resultOctal, point, chiLine)
!         write(*,*) point,resultOCTAL%chiLine(subcell), chiline
!      ENDIF

      IF (PRESENT(etaLine))           etaLine = resultOctal%etaLine(subcell)
      IF (PRESENT(etaCont))           etaCont = resultOctal%etaCont(subcell)
      IF (PRESENT(probDistLine)) probDistLine = resultOctal%probDistLine(subcell)
      IF (PRESENT(probDistCont)) probDistCont = resultOctal%probDistCont(subcell)
      IF (PRESENT(Ne).and.associated(resultOctal%ne))                     Ne = resultOctal%Ne(subcell)
      IF (PRESENT(N))                       N = resultOctal%N(subcell,:)
      IF (PRESENT(nTot))                 nTot = resultOctal%nTot(subcell)
      IF (PRESENT(departCoeff))   departCoeff = resultOctal%departCoeff(subcell,:)
      IF (PRESENT(inFlow))           inFlow = resultOctal%inFlow(subcell)     

      IF (PRESENT(dusttypeFraction))           dusttypeFraction = resultOctal%dusttypeFraction(subcell,:)     

      IF (PRESENT(kappaAbsArray)) THEN
         call returnKappa(grid, resultOctal, subcell, kappaAbsArray=kappaAbsArray)
      ENDIF
      IF (PRESENT(kappaScaArray)) THEN
         call returnKappa(grid, resultOctal, subcell, kappaScaArray=kappaScaArray)
      ENDIF

      IF (PRESENT(rosselandKappa)) THEN
         if (PRESENT(atthistemperature)) then
            call returnKappa(grid, resultOctal, subcell, rosselandKappa=rosselandKappa, &
                 atthisTemperature=atthisTemperature)
         else
            call returnKappa(grid, resultOctal, subcell, rosselandKappa=rosselandKappa)
         endif
      ENDIF

      IF (PRESENT(kappap)) THEN
         if (PRESENT(atthistemperature)) then
            call returnKappa(grid, resultOctal, subcell, kappap=kappap, &
                 atthistemperature=atthistemperature)
         else
            call returnKappa(grid, resultOctal, subcell, kappap=kappap)
         endif
      ENDIF

      IF (PRESENT(kappaAbs)) THEN 
        IF (PRESENT(iLambda)) THEN
!           if (.not.grid%oneKappa) then
!              kappaAbs = resultOctal%kappaAbs(subcell,iLambda)
!           else
              call returnKappa(grid, resultOctal, subcell, ilambda,&
              lambda=lambda, kappaAbs=kappaAbs)
!              if (resultOctal%gasOpacity) then
!                 call returnKappaValue(resultOctal%temperature(subcell), lambda, kappaAbs)
!              else
!                 IF (.NOT.PRESENT(lambda)) THEN
!                    kappaAbs = grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell)
!                 ELSE
!                    kappaAbs = logint(lambda, grid%lamArray(ilambda), grid%lamArray(ilambda+1), &
!                         grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell), &
!                         grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda+1)*resultOctal%rho(subcell))
!                 ENDIF
!           endif
!        endif
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''kappaAbs'' without',&
                   ' ''iLambda''.'
          stop
        END IF
      END IF
      
      IF (PRESENT(kappaSca)) THEN 
        IF (PRESENT(iLambda)) THEN
!           if (.not.grid%oneKappa) then
!              kappaSca = resultOctal%kappaSca(subcell,iLambda)
!           else
              call returnKappa(grid, resultOctal, subcell, &
              ilambda,lambda=lambda, kappaSca=kappaSca)
!              if (resultOctal%gasOpacity) then
!                 call returnKappaValue(resultOctal%temperature(subcell), lambda, kappaSca)
!              else
!                 IF (.NOT.PRESENT(lambda)) THEN
!                    kappaSca = grid%oneKappaSca(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell)
!                 ELSE
!                    kappaSca = logint(lambda, grid%lamArray(ilambda), grid%lamArray(ilambda+1), &
!                         grid%oneKappaSca(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell), &
!                         grid%oneKappaSca(resultOctal%dustType(subcell),iLambda+1)*resultOctal%rho(subcell))
!           ENDIF
!           endif
!        endif
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''kappaSca'' without',&
                   ' ''iLambda''.'
          STOP
       END IF
     END IF

      IF (PRESENT(velocityDeriv)) THEN
        IF (PRESENT(direction) .AND. PRESENT(grid)) THEN
           ! unrotated poistion should be passed here!
          velocityDeriv = amrGridDirectionalDeriv(grid,point,direction,&
                                                  startOctal=resultOctal)
!          velocityDeriv = amrGridDirectionalDeriv(grid,point2,direction,&
!                                                  startOctal=resultOctal)
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''velocityDeriv'' without',&
                   ' ''direction'' and ''grid''.'
          STOP
        END IF
      END IF
    END IF
      

  END SUBROUTINE amrGridValues


  FUNCTION amrGridVelocity(octalTree,point,startOctal,foundOctal,&
                                        foundSubcell,actualSubcell, linearinterp) 
    ! POINT --> should be in unrotated coordinates for 2D case (not projected onto x-z plane!)
    !

    ! returns the velocity at a given point in the grid.
    ! this function can be called with just the first two arguments.
    !   and it will start at the root of the octal tree to locate
    !   the correct octal.
    ! if the foundOctal argument is supplied, it is made to point to 
    !   the octal containing 'point'.
    ! if the foundSubcell argument is supplied, it is made to point to 
    !   the subcell containing 'point'.
    ! if the startOctal argument is supplied, the function uses a 
    !   local search for the correct octal starting at that octal.
    ! if actualSubcell and startSubcell are both supplied, these 
    !   locations are assumed to be correct and no search is performed.

    IMPLICIT NONE

    TYPE(vector)                   :: amrGridVelocity
    TYPE(octal), POINTER           :: octalTree
    TYPE(vector), INTENT(IN)  :: point
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell
    INTEGER, INTENT(IN),  OPTIONAL :: actualSubcell

    TYPE(octal), POINTER           :: resultOctal
    INTEGER                        :: subcell

    TYPE(vector)              :: centre
    real(oct)           :: fac, inc
    real(oct)           :: t1, t2, t3
    real(double)        :: dt1, dt2, dt3
    real(oct)           :: r1, r2, phi1, phi2
    real(double) :: phi
    type(vector) :: newvec
    TYPE(vector) :: point_local, vvec, rHat

    real(double) :: weights(27)
    logical, optional :: linearinterp
    logical :: linear

    if(present(linearinterp)) then
       linear = linearinterp
    else
       linear = .true.
    endif
       
    if (octalTree%threeD) then
       point_local = point
    elseif (octalTree%twoD) then
       ! roate "point" back to z-x plane!
       point_local = projectToXZ(point)
    else 
       ! assume it's threeD for now
       point_local = point
    end if
    if (octalTree%oneD) then
       point_local = VECTOR(modulus(point), 0.d0, 0.d0)
    endif

    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
         ! called with rotated point if necessary for 2d
        CALL findSubcellLocal(point_local,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      resultOctal => startOctal
      
   ELSE
         ! called with rotated point if necessary for 2d
      CALL findSubcellTD(point_local,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal))   foundOctal   => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

   END IF

      inc = 0.5 * resultOctal%subcellSize
      centre = subcellCentre(resultOctal,subcell)
      fac = 1. / resultOctal%subcellsize
      
      t1 = MAX(0.0_oc, fac * (point_local%x - (centre%x - inc)))
      t2 = MAX(0.0_oc, fac * (point_local%y - (centre%y - inc)))
      t3 = MAX(0.0_oc, fac * (point_local%z - (centre%z - inc)))

      if (resultOctal%oneD) then

         select case(subcell)
         case(1)
            vvec = (1.d0-t1) * resultOctal%cornerVelocity(1) + &
                 (   t1) * resultOctal%cornerVelocity(2)
         case(2)
            vvec = (1.d0-t1) * resultOctal%cornerVelocity(2) + &
                 (   t1) * resultOctal%cornerVelocity(3)
         end select
         rHat = point
         call normalize(rHat)
         if (vvec%x >= 0.d0) then
            amrGridVelocity = modulus(vvec) * rHat
         else
            amrGridVelocity = modulus(vvec) * ((-1.d0)*rHat)
         endif
         goto 666
      endif

      if (resultOctal%threed) then
         if (.not.resultOctal%cylindrical) then
            if(linear) then

            SELECT CASE(subcell)
               
            CASE(1)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 1) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 2) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 4) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(10) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(11) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(14)
               
            CASE(2)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 2) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 3) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 6) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(11) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(12) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(14) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(15)
               
            CASE(3)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 4) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 7) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 8) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(14) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(16) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(17)
               
            CASE(4)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 6) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 8) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 9) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(14) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(15) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(17) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(18)
               
            CASE(5)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(10) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(11) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(14) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(19) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(20) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(22) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(23)
               
            CASE(6)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(11) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(12) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(14) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(15) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(20) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(21) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(23) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(24)
            
            CASE(7)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(14) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(16) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(17) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(22) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(23) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(25) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(26)
               
            CASE(8)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(14) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(15) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(17) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(18) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(23) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(24) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(26) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(27)
            CASE DEFAULT
               PRINT *, 'Invalid subcell in amrGridVelocity'
               
            end SELECT
         else

            inc = resultOctal%subcellSize
            centre = resultoctal%centre
            fac = 0.5d0 / resultOctal%subcellsize
      
            dt1 = fac * (point_local%x - (centre%x - inc))
            dt2 = fac * (point_local%y - (centre%y - inc))
            dt3 = fac * (point_local%z - (centre%z - inc))

            call quadint(dt1,dt2,dt3,weights)
            amrgridvelocity = &
            VECTOR(sum(weights(:) * resultoctal%cornervelocity(:)%x), &
                   sum(weights(:) * resultoctal%cornervelocity(:)%y), &
                   sum(weights(:) * resultoctal%cornervelocity(:)%z))
         endif

         else ! cylindrical
            if (resultOctal%splitAzimuthally) then

               centre = subcellCentre(resultOctal,subcell)
               r1 = sqrt(point_local%x**2 + point_local%y**2)
               r2 = sqrt(centre%x**2 + centre%y**2)
               phi1 = atan2(point_local%y, point_local%x)
               if (phi1 < 0.d0) phi1 = phi1 + twoPi
               phi2 = atan2(centre%y, centre%x)
               if (phi2 < 0.d0) phi2 = phi2 + twoPi

               inc = resultOctal%subcellSize / 2.0               

               t1 = MAX(0.0_oc, (r1 - (r2 - inc)) / resultOctal%subcellSize)
               t2 = (phi1 - (phi2 - resultOctal%dPhi/4.d0))/(resultOctal%dPhi/2.d0)
               t3 = MAX(0.0_oc, (point_local%z - (centre%z - inc)) / resultOctal%subcellSize)

               select case(subcell)
            CASE(1)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 1) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 2) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 4) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(10) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(11) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(14)
               
            CASE(2)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 2) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 3) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 6) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(11) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(12) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(14) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(15)
               
            CASE(3)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 4) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 7) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 8) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(14) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(16) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(17)
               
            CASE(4)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 6) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 8) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 9) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(14) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(15) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(17) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(18)
               
            CASE(5)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(10) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(11) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(14) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(19) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(20) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(22) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(23)
               
            CASE(6)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(11) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(12) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(14) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(15) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(20) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(21) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(23) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(24)
            
            CASE(7)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(14) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(16) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(17) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(22) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(23) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(25) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(26)
               
            CASE(8)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(14) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(15) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(17) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(18) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(23) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(24) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(26) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(27)

            CASE DEFAULT
               PRINT *, 'Invalid subcell in amrGridVelocity'
               end select
            else

               centre = subcellCentre(resultOctal,subcell)
               r1 = sqrt(point_local%x**2 + point_local%y**2)
               r2 = sqrt(centre%x**2 + centre%y**2)
               phi1 = atan2(point_local%y, point_local%x)
               if (phi1 < 0.d0) phi1 = phi1 + twoPi
               phi2 = atan2(centre%y, centre%x)
               if (phi2 < 0.d0) phi2 = phi2 + twoPi

               inc = resultOctal%subcellSize * 0.5_oc               

               t1 = MAX(0.0_oc, (r1 - (r2 - inc)) / resultOctal%subcellSize)
               t2 = (phi1 - (phi2 - resultOctal%dPhi*0.5d0))/(resultOctal%dPhi)
               t3 = MAX(0.0_oc, (point_local%z - (centre%z - inc)) / resultOctal%subcellSize)

               select case(subcell)
            CASE(1)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 1) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 2) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 3) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 4) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(7) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(8) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(9) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(10)
               
            CASE(2)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 3) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 4) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 6) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity( 9) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(10) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(11) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(12)
               
            CASE(3)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 7) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 8) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 9) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(10) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(13) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(14) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(15) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(16)
               
            CASE(4)
               amrGridVelocity = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity( 9) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(10) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(11) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(12) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(15) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(16) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(17) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(18)

            CASE DEFAULT
               PRINT *, 'Invalid subcell in amrGridVelocity'
               end select

            endif

         endif

   else
      select case(subcell)
         CASE(1)
            amrGridVelocity = &
                 ((1.d0-t1) * (1.d0-t3)) * resultOctal%cornerVelocity( 1) + &
                 ((     t1) * (1.d0-t3)) * resultOctal%cornerVelocity( 2) + &
                 ((1.d0-t1) * (     t3)) * resultOctal%cornerVelocity( 4) + &
                 ((     t1) * (     t3)) * resultOctal%cornerVelocity( 5)

         CASE(2)
            amrGridVelocity = &
                 ((1.d0-t1) * (1.d0-t3)) * resultOctal%cornerVelocity( 2) + &
                 ((     t1) * (1.d0-t3)) * resultOctal%cornerVelocity( 3) + &
                 ((1.d0-t1) * (     t3)) * resultOctal%cornerVelocity( 5) + &
                 ((     t1) * (     t3)) * resultOctal%cornerVelocity( 6)
         CASE(3)
            amrGridVelocity = &
                 ((1.d0-t1) * (1.d0-t3)) * resultOctal%cornerVelocity( 4) + &
                 ((     t1) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                 ((1.d0-t1) * (     t3)) * resultOctal%cornerVelocity( 7) + &
                 ((     t1) * (     t3)) * resultOctal%cornerVelocity( 8)

         CASE(4)
            amrGridVelocity = &
                 ((1.d0-t1) * (1.d0-t3)) * resultOctal%cornerVelocity( 5) + &
                 ((     t1) * (1.d0-t3)) * resultOctal%cornerVelocity( 6) + &
                 ((1.d0-t1) * (     t3)) * resultOctal%cornerVelocity( 8) + &
                 ((     t1) * (     t3)) * resultOctal%cornerVelocity( 9)
         end select

         phi = atan2(point%y, point%x)
         newVec = rotateZ(amrGridVelocity, -phi)
         amrGridVelocity = newVec
      endif
666   continue

!    endif
  END FUNCTION amrGridVelocity

  FUNCTION amrGridDirectionalDeriv(grid,position,direction,startOctal,&
                                   foundOctal,foundSubcell) 
    !
    ! POINT --> should be in unrotated coordinates for 2D case (not projected onto x-z plane!)
    !

    ! returns the directional derivative of velocity at a given point
    !   in the grid.
    ! this function can be called with just the first three arguments
    !   and it will start at the root of the octal tree to locate
    !   the correct octal.
    ! if the foundOctal argument is supplied, it is made to point to 
    !   the octal containing 'position'.
    ! if the foundSubcell argument is supplied, it is made to point to 
    !   the subcell containing 'position'.
    ! if the startOctal argument is supplied, the function uses a 
    !   local search for the correct octal starting at that octal.
    !   NOTE that startOctal may be *changed* by this function! 

    ! Note that the results of this function DO have a 1^10 factor in them. 
    !   (they are in cSpeed, but the distance is calculated 1e10cm)

    ! THIS SHOULD WORK ALSO IN TWO-D CASE as long as position and direction
    ! are both are not projected on z-x plane (un-rotated coodinates).

    IMPLICIT NONE

    TYPE(gridtype), INTENT(IN)     :: grid 
    REAL                           :: amrGridDirectionalDeriv
    TYPE(vector), INTENT(IN)  :: position
    TYPE(vector), INTENT(IN)       :: direction
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    TYPE(octal),  POINTER         :: thisOctal !, currentOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell

    TYPE(vector)              :: octalDirection
    TYPE(octal), POINTER           :: firstOctal
    real(oct)           :: dr, dx, dphi
    real(oct)           :: r
    real(oct)           :: phi1, phi2
    TYPE(vector)              :: position1
    TYPE(vector)              :: position2
    INTEGER                        :: subcell


    if (.not.grid%lineEmission) then
       amrGridDirectionalDeriv = 1.e30
       goto 666
    endif
    
    octalDirection = direction

    call amrGridValues(grid%octreeRoot,position,foundOctal=thisOctal,&
                           foundSubcell=subcell)
    
!    currentOctal => grid%octreeRoot
!    thisOctal => grid%octreeRoot
!    CALL findSubcellTD(position, currentOctal, thisOctal, subcell)
    ! dr is a small increment of distance
    dr = thisOctal%subcellSize * 1.d-1

    ! get a new position a little way back from current position
    position1 = position - (dr * octalDirection)
    
    ! this might be inside core or outside grid - in which case
    !   just use the current position as the first point

    r = modulus(position1)
    IF (.NOT. inOctal(grid%octreeRoot,position1) .OR. (r < grid%rCore)) THEN
      position1 = position
    END IF
      
    ! first line of sight velocity
    phi1 = direction .dot. (AMRgridVelocity(grid%octreeRoot,position1,&
                                            startOctal=startOctal,&
                                            foundOctal=firstOctal))
    
    ! now go forward a bit from current position
    position2 = position + (dr * octalDirection)

    ! check we're still inside grid
    r = modulus(position2)
    IF (.NOT. inOctal(grid%octreeRoot,position2) .OR. (r < grid%rCore)) THEN
      position2 = position
    END IF

    ! the second position l.o.s. velocity
    phi2 = direction .dot. (AMRgridVelocity(grid%octreeRoot,position2,&
                                            startOctal=firstOctal,&
                                            foundOctal=foundOctal,&
                                            foundSubcell=subcell))

    IF (PRESENT(foundSubcell)) foundSubcell = subcell                                        

    dx = modulus(position2 - position1)

    dphi = phi2 - phi1

    ! the line of sight velocity gradient

    IF (ABS(dx) >= 1.e-20) THEN
       amrGridDirectionalDeriv = (abs(dphi / dx))
    ELSE
       amrGridDirectionalDeriv = 1.e-20
    ENDIF
    IF (.NOT. ABS(amrGridDirectionalDeriv) >= 1.e-10) &
       amrGridDirectionalDeriv = 1.e-10
666 continue
  END FUNCTION amrGridDirectionalDeriv


  FUNCTION amrGridTemperature(octalTree,point,startOctal,foundOctal,& 
                                        foundSubcell,actualSubcell) 
    ! POINT --> can be in both in roteated or unrotated coordinates for 2D case 
    !

    ! returns the temperature at a given point in the grid.
    ! this function can be called with just the first two arguments
    !   and it will start at the root of the octal tree to locate
    !   the correct octal.
    ! if the foundOctal argument is supplied, it is made to point to 
    !   the octal containing 'point'.
    ! if the startOctal argument is supplied, the function uses a 
    !   local search for the correct octal starting at that octal.
    ! if actualSubcell and startSubcell are both supplied, these 
    !   locations are assumed to be correct and no search is performed.

    IMPLICIT NONE

    REAL                           :: amrGridTemperature
    TYPE(octal), POINTER           :: octalTree
    TYPE(vector), INTENT(IN)  :: point
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell
    INTEGER, INTENT(IN),  OPTIONAL :: actualSubcell

    TYPE(octal), POINTER           :: resultOctal
    INTEGER                        :: subcell
    TYPE(vector) :: point_local

    if (octalTree%threeD) then
       point_local = point
    elseif (octalTree%twoD) then
       ! roate "point" back to z-x plane!
       point_local = projectToXZ(point)
    else 
       ! assume it's threeD for now
       point_local = point
    end if
    if (octalTree%oneD) then
       point_local = VECTOR(modulus(point), 0.d0, 0.d0)
    endif

    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point_local,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      
      ! should add in some interpolation routine
      amrGridTemperature = startOctal%temperature(subcell)
      
    ELSE
       ! rotated point for 2d case
      CALL findSubcellTD(point_local,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal)) foundOctal => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

      ! should add in some interpolation routine
      amrGridTemperature = resultOctal%temperature(subcell)

    END IF

  END FUNCTION amrGridTemperature


  RECURSIVE SUBROUTINE locateLineProbAMR(probability,thisOctal,subcell) 
    ! finds the subcell that contains a given value of 'probability'.
    ! each subcell of the tree's octals has a value for line emission 
    !   probability which is an upper bound of the subcell's value in the cumulative probability
    !   distribution for the 

    IMPLICIT NONE

    real(double), INTENT(IN) :: probability
    TYPE(octal), POINTER              :: thisOctal
    INTEGER, INTENT(OUT)              :: subcell

    INTEGER              :: i, j 
   
   
    ! we need to treat the first subcell as a special case
    
    IF (probability < thisOctal%probDistLine(1)) THEN

      IF (thisOctal%hasChild(1)) THEN
      
        ! find the child
        DO j = 1, thisOctal%nChildren, 1
          IF (thisOctal%indexChild(j) == 1) THEN
            thisOctal => thisOctal%child(j)
            CALL locateLineProbAMR(probability,thisOctal,subcell)
            RETURN
          END IF
        END DO
        
      ELSE 
        subcell = 1
        RETURN
        
      END IF
    END IF   
  
   
    DO i = 2, thisOctal%maxChildren, 1
      
      IF (probability > thisOctal%probDistLine(i-1) .AND. &
          probability < thisOctal%probDistLine(i)) THEN
      
        IF (thisOctal%hasChild(i)) THEN
          
          ! find the child
          DO j = 1, thisOctal%nChildren, 1
            IF (thisOctal%indexChild(j) == i) THEN
              thisOctal => thisOctal%child(j)
              CALL locateLineProbAMR(probability,thisOctal,subcell)
              RETURN
            END IF
          END DO
          
        ELSE 
          subcell = i
          RETURN
        
        END IF
      END IF
    END DO

      
  END SUBROUTINE locateLineProbAMR


  RECURSIVE SUBROUTINE locateContProbAMR(probability,thisOctal,subcell) 
    ! finds the subcell that contains a given value of 'probability'.
    ! each subcell of the tree's octals has a value for continuous emission 
    !   probability which is an upper bound of the subcell's value in the cumulative probability
    !   distribution for the 

    IMPLICIT NONE

    real(double), INTENT(IN) :: probability
    TYPE(octal), POINTER              :: thisOctal
    INTEGER, INTENT(OUT)              :: subcell

    INTEGER              :: i, j 
   
!    write(*,'(9f9.5)') probability,thisOctal%probDistCont(1:8)
    ! we need to treat the first subcell as a special case
    
    IF (probability < thisOctal%probDistCont(1)) THEN
    
      IF (thisOctal%hasChild(1)) THEN
      
        ! find the child
        DO j = 1, thisOctal%nChildren, 1
          IF (thisOctal%indexChild(j) == 1) THEN
            thisOctal => thisOctal%child(j)
            CALL locateContProbAMR(probability,thisOctal,subcell)
            RETURN
          END IF
        END DO
        
      ELSE 
        subcell = 1
        RETURN
        
      END IF
    END IF   
  
   
    DO i = 2, thisOctal%maxChildren, 1
    
      IF (probability > thisOctal%probDistCont(i-1) .AND. &
          probability < thisOctal%probDistCont(i)) THEN
      
        IF (thisOctal%hasChild(i)) THEN
          
          ! find the child
          DO j = 1, thisOctal%nChildren, 1
            IF (thisOctal%indexChild(j) == i) THEN
               thisOctal => thisOctal%child(j)
              CALL locateContProbAMR(probability,thisOctal,subcell)
              RETURN
            END IF
          END DO
          
        ELSE 
          subcell = i
          RETURN
          
        END IF
      END IF
    END DO

      
  END SUBROUTINE locateContProbAMR


  FUNCTION whichSubcell(thisOctal,point) RESULT (subcell)
    ! POINT --> MUST be pre-rotated for 2d case!!!!!!!!!!!!!!
    !
    ! returns the identification number (1-8) of the subcell of the 
    ! current octal which contains a given point
    ! NB this does NOT check that the point lies within the bounds of the octal!



    IMPLICIT NONE

    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(vector), INTENT(IN) :: point
    INTEGER                       :: subcell
    real(double) :: r, phi

    if (thisOctal%oneD) then
       if (point%x <= thisOctal%centre%x) then
          subcell = 1
       else
          subcell = 2
       endif
       goto 666
    endif


    if (thisOctal%threed) then ! threed case 

       if (.not.thisOctal%cylindrical) then ! cartesian case

          IF ( point%x < thisOctal%centre%x ) THEN
             IF ( point%y < thisOctal%centre%y ) THEN
                IF ( point%z < thisOctal%centre%z ) THEN
                   subcell = 1
                ELSE 
                   subcell = 5
                   
                ENDIF
             ELSE 
                IF (point%z < thisOctal%centre%z) THEN
                   subcell = 3
                ELSE 
                   subcell = 7
                ENDIF
             END IF
          ELSE
             IF (point%y < thisOctal%centre%y) THEN
                IF (point%z < thisOctal%centre%z) THEN
                   subcell = 2
                ELSE 
                   subcell = 6
                ENDIF
             ELSE 
                IF (point%z < thisOctal%centre%z) THEN
                   subcell = 4
                ELSE 
                   subcell = 8
                ENDIF
             END IF
          ENDIF ! cartesian case

       else ! cylindrical case

          r = sqrt(point%x**2+point%y**2)
          phi = atan2(point%y, point%x)
          if (phi < 0.d0) phi = phi + twoPi

          if (thisOctal%splitAzimuthally) then ! azimuthal split case

             if (phi <= thisOctal%phi) then
                IF ( r <= thisOctal%r) THEN
                   IF ( point%z <= thisOctal%centre%z ) THEN
                      subcell = 1
                   ELSE 
                      subcell = 3
                   ENDIF
                ELSE
                   IF (point%z <= thisOctal%centre%z) THEN
                      subcell = 2
                   ELSE 
                      subcell = 4
                   ENDIF
                END IF
             else
                IF ( r <= thisOctal%r ) THEN
                   IF ( point%z <= thisOctal%centre%z ) THEN
                      subcell = 5
                   ELSE 
                      subcell = 7
                   ENDIF
                ELSE
                   IF (point%z <= thisOctal%centre%z) THEN
                      subcell = 6
                   ELSE 
                      subcell = 8
                   ENDIF
                END IF
             endif
          else

             IF ( r <= thisOctal%r ) THEN
                IF ( point%z <= thisOctal%centre%z ) THEN
                   subcell = 1
                ELSE 
                   subcell = 3
                ENDIF
                ELSE
                   IF (point%z <= thisOctal%centre%z) THEN
                      subcell = 2
                   ELSE 
                      subcell = 4
                   ENDIF
                END IF
             endif ! azi case
          endif ! cylindrical

    else ! twoD case
       
       IF ( point%x <= thisOctal%centre%x ) THEN
          IF ( point%z <= thisOctal%centre%z ) THEN
             subcell = 1
          ELSE 
             subcell = 3
          ENDIF
       ELSE
          IF (point%z <= thisOctal%centre%z) THEN
             subcell = 2
          ELSE 
             subcell = 4
          ENDIF
       END IF
    endif

666 continue

  END FUNCTION whichSubcell    


  FUNCTION inOctal(thisOctal,point,alreadyRotated) 
    ! true if the point lies within the boundaries of the current octal
  
    use input_variables, only : hydrodynamics
    IMPLICIT NONE
    LOGICAL                       :: inOctal
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(vector), INTENT(IN) :: point
    TYPE(vector)             :: octVec2D
    real(double)                  :: r, phi, dphi, eps
    logical, optional :: alreadyRotated
    logical :: doRotate

    doRotate = .true.
    if (PRESENT(alreadyRotated)) doRotate = .not.alreadyRotated

    if (thisOctal%threeD) then
       if (.not.thisOctal%cylindrical) then
          IF (point%x < thisOctal%xMin) THEN ; inOctal = .FALSE. ; goto 101
          ELSEIF (point%x > thisOctal%xMax) THEN ; inOctal = .FALSE.; goto 101
          ELSEIF (point%y < thisOctal%yMin) THEN ; inOctal = .FALSE.; goto 101
          ELSEIF (point%y > thisOctal%yMax) THEN ; inOctal = .FALSE.; goto 101
          ELSEIF (point%z < thisOctal%zMin) THEN ; inOctal = .FALSE.; goto 101
          ELSEIF (point%z > thisOctal%zMax) THEN ; inOctal = .FALSE.; goto 101
          ELSE  
             inOctal = .TRUE.
          ENDIF
101       continue
       else
          phi = atan2(point%y,point%x)
          if (phi < 0.d0) phi = phi + twoPi
          dphi = abs(phi - thisOctal%phi)
          r = sqrt(point%x**2 + point%y**2)
          eps =  0.d0 ! sqrt(epsilon(thisOctal%subcellsize))
          IF     (r < thisOctal%r - thisOctal%subcellSize - eps) THEN ; inOctal = .FALSE. 
          ELSEIF (r > thisOctal%r + thisOctal%subcellSize + eps) THEN ; inOctal = .FALSE.
          ELSEIF (dphi > thisOctal%dphi/2.d0) THEN ; inOctal = .FALSE.
          ELSEIF (point%z < thisOctal%zMin) THEN ; inOctal = .FALSE.
          ELSEIF (point%z > thisOctal%zMax) THEN ; inOctal = .FALSE.
          ELSE  
             inOctal = .TRUE.
          ENDIF
       endif
    else ! twoD case
       if (.not.hydrodynamics) then
          if (doRotate) then
             octVec2D = projectToXZ(point)
          else
             octVec2D = point   
          endif
       else
          octVec2D = point
       endif
       IF (octVec2D%x < thisOctal%xMin) THEN ; inOctal = .FALSE. 
       ELSEIF (octVec2D%x > thisOctal%xMax) THEN ; inOctal = .FALSE.
       ELSEIF (octVec2D%z < thisOctal%zMin) THEN ; inOctal = .FALSE.
       ELSEIF (octVec2D%z > thisOctal%zMax) THEN ; inOctal = .FALSE.
       ELSE  
          inOctal = .TRUE.
       ENDIF
    endif

    if (thisOctal%oneD) then
       r = modulus(point)
       if ( r < thisOctal%centre%x  - thisOctal%subcellSize) then ; inoctal = .false.
       else if (r > thisOctal%centre%x + thisOctal%subcellSize) then; inOctal = .false.
       else
          inOctal = .true.
       endif
       goto 666
    endif


666 continue
  END FUNCTION inOctal

  FUNCTION inSubcell(thisOctal,thisSubcell,point) 
    ! true if the point lies within the boundaries of the current octal
  
    IMPLICIT NONE
 
    LOGICAL                       :: inSubcell
    TYPE(octal), INTENT(IN)       :: thisOctal
    INTEGER, INTENT(IN)           :: thisSubcell
    TYPE(vector), INTENT(IN) :: point

    IF (inOctal(thisOctal,point)) THEN
      inSubcell = whichSubcell(thisOctal,point) == thisSubcell
    ELSE
      inSubcell = .FALSE.
    END IF
  
  END FUNCTION inSubcell

  FUNCTION looseInOctal(thisOctal,point) 
    ! true if the point lies 'loosely' in the current octal
    ! ( a 10% margin of error is allowed )
    ! this is useful only for testing purposes!
  
    IMPLICIT NONE
 
    LOGICAL                       :: looseInOctal
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(vector), INTENT(IN) :: point

    if (thisOctal%threeD) then
       IF ((point%x <= thisOctal%centre%x - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%x >= thisOctal%centre%x + 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%y <= thisOctal%centre%y - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%y >= thisOctal%centre%y + 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%z <= thisOctal%centre%z - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%z >= thisOctal%centre%z + 1.1_oc * thisOctal%subcellSize )) THEN
          looseInOctal = .FALSE.
       ELSE  
          looseInOctal = .TRUE.
       ENDIF
    else
       IF ((point%x <= thisOctal%centre%x - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%x >= thisOctal%centre%x + 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%z <= thisOctal%centre%z - 1.1_oc * thisOctal%subcellSize ) .OR. &
            (point%z >= thisOctal%centre%z + 1.1_oc * thisOctal%subcellSize )) THEN
          looseInOctal = .FALSE.
       ELSE  
          looseInOctal = .TRUE.
       ENDIF
    endif
  END FUNCTION looseInOctal

  SUBROUTINE smoothAMRgrid(grid,factor, stellar_cluster, inheritProps, interpProps, &
       romData)
    ! checks whether each octal's neighbours are much bigger than it, 
    !   if so, makes the neighbours smaller.

    TYPE(gridtype), INTENT(INOUT) :: grid 
    REAL, INTENT(IN)              :: factor
    TYPE(cluster), optional, intent(in)  :: stellar_cluster
    LOGICAL, INTENT(IN),optional  :: inheritProps
    logical, intent(in), optional :: interpProps 
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
    
    LOGICAL :: gridConverged

    CALL setAllUnchanged(grid%octreeRoot)

    DO 
      gridConverged = .TRUE.
      CALL smoothAMRgridPrivate(grid%octreeRoot,grid,gridConverged,romData)
      IF ( gridConverged ) EXIT 
    END DO
    
  CONTAINS
    
    RECURSIVE SUBROUTINE smoothAMRgridPrivate(thisOctal,grid,gridConverged, &
         romData)

!      TYPE(octal), INTENT(INOUT), TARGET :: thisOctal
      TYPE(octal), INTENT(INOUT) :: thisOctal
      TYPE(gridtype), INTENT(INOUT   ) :: grid 
      LOGICAL, INTENT(INOUT)               :: gridConverged
      TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry

      INTEGER              :: i
      REAL(oct) :: halfSmallestSubcell
      REAL(oct) :: offset
      TYPE(octal), POINTER :: neighbour
      TYPE(vector), ALLOCATABLE, DIMENSION(:) :: locator
      TYPE(vector) :: aHat
      INTEGER              :: subcell
      INTEGER              :: nLocator ! number of locators (4 for twoD, 6 for threed)

      ! we will find the coordinates of a point that lies outside the current
      !   octal. we then compare the size of the cell that contains that point
      !   with the size of the current cell, if it is bigger be more than a 
      !   factor of 'factor', we subdivide the neighbouring cell.
      ! we do this in each of six directions

      ! we do not have to test the other subcells in the current octal because
      !   they can be smaller than the any of the other subcells, but they
      !   cannot be *bigger*. this saves some time.

      if (thisOctal%threed) then
         if (.not.thisOctal%cylindrical) then
            nlocator = 6
         else
            nlocator = 4
         endif
      else
        nlocator = 4
      endif

      ALLOCATE(locator(nlocator))

      ! we find points which are outside the current octal by a distance
      !   equivalent to half the size of the tree's smallest subcell.
      halfSmallestSubcell = grid%halfSmallestSubcell

      ! we also add a slight offset to our test positions to avoid testing 
      !   at cell boundaries.
      offset = halfSmallestSubcell / 2.0_oc

      IF ( thisOctal%threed ) THEN
        locator(:) = thisOctal%centre + ( offset * vector(1.0_oc,1.0_oc,1.0_oc) )
      ELSE
        locator(:) = thisOctal%centre + ( offset * vector(1.0_oc,0.0_oc,1.0_oc) )
      END IF

      IF ( thisOctal%threeD ) THEN
         if (.not.thisOctal%cylindrical) then
            locator(1)%x = thisOctal%centre%x + thisOctal%subcellSize + halfSmallestSubcell
            locator(2)%y = thisOctal%centre%y + thisOctal%subcellSize + halfSmallestSubcell
            locator(3)%z = thisOctal%centre%z + thisOctal%subcellSize + halfSmallestSubcell
            locator(4)%x = thisOctal%centre%x - thisOctal%subcellSize - halfSmallestSubcell
            locator(5)%y = thisOctal%centre%y - thisOctal%subcellSize - halfSmallestSubcell
            locator(6)%z = thisOctal%centre%z - thisOctal%subcellSize - halfSmallestSubcell
         else
            locator(:) = thisOctal%centre
            locator(1) = locator(1) + (thisOctal%subcellSize + halfSmallestSubcell) * zHat
            locator(2) = locator(2) - (thisOctal%subcellSize + halfSmallestSubcell) * zHat
            aHat = VECTOR(thisOctal%centre%x,thisOctal%centre%y,0.d0)
            call normalize(aHat)
            locator(3) = locator(3) + (thisOctal%subcellSize + halfSmallestSubcell) * aHat
            locator(4) = locator(4) - (thisOctal%subcellSize + halfSmallestSubcell) * aHat
         endif
      ELSE
        locator(1)%x = thisOctal%centre%x + thisOctal%subcellSize + halfSmallestSubcell
        locator(2)%z = thisOctal%centre%z + thisOctal%subcellSize + halfSmallestSubcell
        locator(3)%x = thisOctal%centre%x - thisOctal%subcellSize - halfSmallestSubcell
        locator(4)%z = thisOctal%centre%z - thisOctal%subcellSize - halfSmallestSubcell
      ENDIF
      
      DO i = 1, nLocator, 1
        IF ( inOctal(grid%octreeRoot,locator(i)) ) THEN
          CALL findSubcellTD(locator(i),grid%octreeRoot,neighbour,subcell)
          !neighbour => thisOctal
          !CALL findSubcellLocal(locator(i),neighbour,subcell)

          IF ( neighbour%subcellSize > (factor * thisOctal%subcellSize) ) THEN
            IF ( neighbour%hasChild(subcell) ) THEN 
              PRINT *, "neighbour already has child. (B)"
              !STOP
              do ;end do
            END IF
              neighbour%changed(1:neighbour%maxChildren) = .TRUE.
              call addNewChild(neighbour, subcell, grid, adjustGridInfo=.TRUE.,  &
                               stellar_cluster=stellar_cluster, &
                               inherit=inheritProps, interp=interpProps, romData=romData)
              gridConverged = .FALSE.
          ENDIF
        END IF
          
! force return until algorithm is fixed
IF ( .NOT. gridConverged ) RETURN
          
      END DO

      ! call this subroutine recursively on each of any children.
      IF ( (.NOT. ANY(thisOctal%Changed)) .AND. thisOctal%nChildren > 0 ) THEN
        DO i = 1, thisOctal%nChildren, 1 
          CALL smoothAMRgridPrivate(thisOctal%child(i),grid,gridConverged,romData=romData)
          !CALL checkAMRgrid(grid,checkNoctals=.FALSE.)                                  
        
! force return until algorithm is fixed
IF ( .NOT. gridConverged ) RETURN
        END DO
      END IF

      DEALLOCATE(locator)
          
    END SUBROUTINE smoothAMRgridPrivate
  
  END SUBROUTINE smoothAMRgrid

  RECURSIVE SUBROUTINE findSubcellTD(point,currentOctal,resultOctal,subcell)
  ! finds the octal (and that octal's subcell) containing a point.
  !   only searches in downwards direction (TD = top-down) , so
  !   probably best to start from root of tree

    IMPLICIT NONE

    TYPE(vector), INTENT(IN) :: point
    TYPE(octal), POINTER :: currentOctal
    TYPE(octal), POINTER :: resultOctal
    INTEGER, INTENT(OUT) :: subcell
    TYPE(octal), POINTER :: child

    INTEGER :: i

    resultOctal => currentOctal
    subcell = whichSubcell(currentOctal,point)

    IF ( currentOctal%hasChild(subcell) ) THEN
      ! search the index to see where it is stored
      DO i = 1, currentOctal%maxChildren, 1
        IF ( currentOctal%indexChild(i) == subcell ) THEN
          child => currentOctal%child(i)
          CALL findSubcellTD(point,child,resultOctal,subcell)
          EXIT
        END IF
      END DO
    END IF

  END SUBROUTINE findSubcellTD

  RECURSIVE SUBROUTINE findSubcellTDLevel(point,currentOctal,resultOctal,subcell,nDepth)
  ! finds the octal (and that octal's subcell) containing a point.
  !   only searches in downwards direction (TD = top-down) , so
  !   probably best to start from root of tree

    IMPLICIT NONE

    TYPE(vector), INTENT(IN) :: point
    TYPE(octal), POINTER :: currentOctal
    TYPE(octal), POINTER :: resultOctal
    INTEGER, INTENT(OUT) :: subcell
    integer :: nDepth
    TYPE(octal), POINTER :: child

    INTEGER :: i

    resultOctal => currentOctal
    subcell = whichSubcell(currentOctal,point)

    IF ( currentOctal%hasChild(subcell).and.(currentOctal%nDepth < nDepth) ) THEN
      ! search the index to see where it is stored
      DO i = 1, currentOctal%maxChildren, 1
        IF ( currentOctal%indexChild(i) == subcell ) THEN
          child => currentOctal%child(i)
          CALL findSubcellTDLevel(point,child,resultOctal,subcell,ndepth)
          EXIT
        END IF
      END DO
    END IF

  END SUBROUTINE findSubcellTDLevel


  SUBROUTINE findSubcellLocal(point,thisOctal,subcell,  prob)
    ! finds the octal (and that octal's subcell) containing a point.
    !   starts searching from the current octal, and goes up and down the
    !   tree as needed to find the correct octal.
    use input_variables, only : hydrodynamics, suppresswarnings
    IMPLICIT NONE
    TYPE(vector), INTENT(IN) :: point
    TYPE(vector) :: point_local
    TYPE(octal),POINTER    :: thisOctal
    INTEGER, INTENT(OUT)   :: subcell
    LOGICAL, INTENT(OUT),optional   :: prob
    
    LOGICAL                :: haveDescended    ! see comments below
    LOGICAL                :: boundaryProblem  ! see comments below
    
                             
    haveDescended = .FALSE.   ! if the 'point' lies very close to an 
    boundaryProblem = .FALSE. !   boundary, the program may go into 
                              !   a loop going up and down the tree.
                              ! we will keep track of the progress of
                              !   the search using these flags.
                             
    if (thisOctal%twoD) then
       if (.not.hydrodynamics) then
          point_local = projectToXZ(point)
       else
          point_local = point
       endif
    else
       point_local = point
    endif
    if (thisOctal%oneD) then
       if (.not.hydrodynamics) then
          point_local = VECTOR(modulus(point), 0.d0, 0.d0)
       else
          point_local = point
       endif
    endif


    CALL findSubcellLocalPrivate(point_local,thisOctal,subcell,&
                                 haveDescended,boundaryProblem)

    if (present(prob)) then
      prob = boundaryProblem
    else
      if (boundaryProblem .and. .not. suppresswarnings) then
        call torus_abort("Torus aborting due to panic in findSubcellLocal")
     else
        CALL findSubcellTD(point_local,thisoctal,thisOctal,subcell)
     endif
    endif
                                 
  CONTAINS

    RECURSIVE SUBROUTINE findSubcellLocalPrivate(point,thisOctal,subcell,&
                                                 haveDescended,boundaryProblem)
      use input_variables, only : suppressWarnings
      TYPE(vector), INTENT(IN) :: point
      TYPE(octal),POINTER    :: thisOctal
      INTEGER, INTENT(OUT)   :: subcell
      LOGICAL, INTENT(INOUT) :: haveDescended
      LOGICAL, INTENT(INOUT) :: boundaryProblem
!      type(vector) :: rVec
      INTEGER :: i
      
      IF ( inOctal(thisOctal,point,alreadyRotated=.true.) ) THEN

        haveDescended = .TRUE. ! record that we have gone down the tree.
      
        ! if the point lies within the current octal, we identify the
        !   subcell
        subcell = whichSubcell(thisOctal,point)

        ! if a problem has been detected, this is where we complete the search
        IF (boundaryProblem) RETURN 
      
        ! if the subcell has a child, we look in the child for the point
        IF ( thisOctal%hasChild(subcell) ) THEN
                
          ! search the index to see where it is stored
          DO i = 1, thisOctal%nChildren, 1
            IF ( thisOctal%indexChild(i) == subcell ) THEN
                    
              thisOctal => thisOctal%child(i)


              CALL findSubcellLocalPrivate(point,thisOctal,subcell,haveDescended,boundaryProblem)
              RETURN
              
            END IF
          END DO
          
        ELSE 
          RETURN
        END IF

      ELSE
        ! if the point is outside the current octal, we look in its
        !   parent octal

        ! first check that we are not outside the grid
        IF ( thisOctal%nDepth == 1 ) THEN
           if(.not. suppresswarnings) then
           write(*,*) "octal", thisOctal%ndepth
           write(*,*) "inoctal min", thisOctal%xMin,thisOctal%yMin,thisOctal%zMin
           write(*,*) "inoctal max", thisOctal%xMax,thisOctal%yMax,thisOctal%zMax

          PRINT *, 'Panic: In findSubcellLocalPrivate, point is outside the grid'
          write(*,*) point
          write(*,*) sqrt(point%x**2+point%y**2)
          write(*,*) atan2(point%y,point%x)*radtodeg
          write(*,*) " "
          write(*,*) thisOctal%centre
          write(*,*) thisOctal%subcellSize
!          write(*,*) thisOctal%phi*radtodeg,thisOctal%dphi*radtodeg
          write(*,*) sqrt(thisOctal%centre%x**2+thisOctal%centre%y**2)
       endif
          if(.not. suppresswarnings) then
                STOP
          endif
          boundaryProblem = .TRUE.
          RETURN
       END IF
     
        ! if we have previously gone down the tree, and are now going back up, there
        !   must be a problem.
        IF (haveDescended) then
           boundaryProblem = .TRUE.
           PRINT *, 'Panic: In findSubcellLocalPrivate, have descended and are now going back up'
           write(*,*) "rank ",myrankglobal, thisOctal%mpithread(1:8)
           write(*,*) "split az ",thisOctal%splitAzimuthally
           write(*,*) point
           write(*,*) atan2(point%y,point%x)*radtodeg
           write(*,*) sqrt(point%x**2 + point%y**2)
           write(*,*) " "
           write(*,*) thisOctal%nDepth
           write(*,*) thisOctal%centre
           write(*,*) thisOctal%subcellSize
           write(*,*) thisOctal%phi*radtodeg,thisOctal%dphi*radtodeg
           write(*,*) sqrt(thisOctal%centre%x**2+thisOctal%centre%y**2)
           write(*,*) atan2(thisOctal%centre%y,thisOctal%centre%x)*radtodeg
           write(*,*) " x min/max, z min max ",thisOctal%xMin, thisOctal%xMax, thisOctal%zMin, thisOctal%zMax
           write(*,*) "parent x min/max, z min max ",thisOctal%parent%xMin, thisOctal%parent%xMax, thisOctal%parent%zMin, &
                thisOctal%parent%zMax
           write(*,*) "cen ",thisOctal%centre
           write(*,*) "size ",thisOctal%subcellsize
!           rVec = subcellCentre(thisOctal,subcell)
!           write(*,*) rVec%x+thisOctal%subcellSize/2.
!           write(*,*) rVec%x-thisOctal%subcellSize/2.
!           write(*,*) rVec%y+thisOctal%subcellSize/2.
!           write(*,*) rVec%y-thisOctal%subcellSize/2.
!           write(*,*) rVec%z+thisOctal%subcellSize/2.
!           write(*,*) rVec%z-thisOctal%subcellSize/2.
!           do ; enddo
!           STOP
           return
        endif
        
        IF ( thisOctal%nDepth /= 1 ) THEN
           thisOctal => thisOctal%parent
        ENDIF

        CALL findSubcellLocalPrivate(point,thisOctal,subcell,haveDescended,boundaryProblem)
       
     END IF
    
    END SUBROUTINE findSubcellLocalPrivate

  END SUBROUTINE findSubcellLocal

  SUBROUTINE findSubcellLocalLevel(point,thisOctal,subcell, nDepth, prob)
    ! finds the octal (and that octal's subcell) containing a point.
    !   starts searching from the current octal, and goes up and down the
    !   tree as needed to find the correct octal.
    use input_variables, only : hydrodynamics
    IMPLICIT NONE
    integer :: nDepth
    TYPE(vector), INTENT(IN) :: point
    TYPE(vector) :: point_local
    TYPE(octal),POINTER    :: thisOctal
    INTEGER, INTENT(OUT)   :: subcell
    LOGICAL, INTENT(OUT),optional   :: prob
    
    LOGICAL                :: haveDescended    ! see comments below
    LOGICAL                :: boundaryProblem  ! see comments below
    
                             
    haveDescended = .FALSE.   ! if the 'point' lies very close to an 
    boundaryProblem = .FALSE. !   boundary, the program may go into 
                              !   a loop going up and down the tree.
                              ! we will keep track of the progress of
                              !   the search using these flags.
                             
    if (thisOctal%twoD) then
       if (.not.hydrodynamics) then
          point_local = projectToXZ(point)
       else
          point_local = point
       endif
    else
       point_local = point
    endif
    if (thisOctal%oneD) then
       if (.not.hydrodynamics) then
          point_local = VECTOR(modulus(point), 0.d0, 0.d0)
       else
          point_local = point
       endif
    endif


    CALL findSubcellLocalPrivateLevel(point_local,thisOctal,subcell,&
                                 haveDescended,boundaryProblem, nDepth)
    if (present(prob)) then
      prob = boundaryProblem
    else
      if (boundaryProblem) then
        stop 1
      endif
    endif
                                 
  CONTAINS

    RECURSIVE SUBROUTINE findSubcellLocalPrivateLevel(point,thisOctal,subcell,&
                                                 haveDescended,boundaryProblem, nDepth)
      TYPE(vector), INTENT(IN) :: point
      TYPE(octal),POINTER    :: thisOctal
      INTEGER, INTENT(OUT)   :: subcell
      integer :: nDepth
      LOGICAL, INTENT(INOUT) :: haveDescended
      LOGICAL, INTENT(INOUT) :: boundaryProblem
!      type(vector) :: rVec
      INTEGER :: i
      
      IF ( inOctal(thisOctal,point,alreadyRotated=.true.) ) THEN

        haveDescended = .TRUE. ! record that we have gone down the tree.
      
        ! if the point lies within the current octal, we identify the
        !   subcell
        subcell = whichSubcell(thisOctal,point)

        ! if a problem has been detected, this is where we complete the search
        IF (boundaryProblem) RETURN 
      
        ! if the subcell has a child, we look in the child for the point
        IF ( thisOctal%hasChild(subcell).and.(thisOctal%ndepth < nDepth) ) THEN
                
          ! search the index to see where it is stored
          DO i = 1, thisOctal%nChildren, 1
            IF ( thisOctal%indexChild(i) == subcell ) THEN
                    
              thisOctal => thisOctal%child(i)


              CALL findSubcellLocalPrivateLevel(point,thisOctal,subcell,haveDescended,boundaryProblem,nDepth)
              RETURN
              
            END IF
          END DO
          
        ELSE 
          RETURN
        END IF

      ELSE
        ! if the point is outside the current octal, we look in its
        !   parent octal

        ! first check that we are not outside the grid
        IF ( thisOctal%nDepth == 1 ) THEN
          PRINT *, 'Panic: In findSubcellLocal, point is outside the grid'
          write(*,*) point
          write(*,*) sqrt(point%x**2+point%y**2)
          write(*,*) atan2(point%y,point%x)*radtodeg
          write(*,*) " "
          write(*,*) thisOctal%centre
          write(*,*) thisOctal%subcellSize
          write(*,*) thisOctal%phi*radtodeg,thisOctal%dphi*radtodeg
          write(*,*) sqrt(thisOctal%centre%x**2+thisOctal%centre%y**2)
           DO ; END DO
          STOP
          boundaryProblem = .TRUE.
          RETURN
        END IF
     
        ! if we have previously gone down the tree, and are now going back up, there
        !   must be a problem.
        IF (haveDescended) then
           boundaryProblem = .TRUE.
           PRINT *, 'Panic: In findSubcellLocal, have descended and are now going back up'
           write(*,*) point
           write(*,*) atan2(point%y,point%x)*radtodeg
           write(*,*) sqrt(point%x**2 + point%y**2)
           write(*,*) " "
           write(*,*) thisOctal%nDepth
           write(*,*) thisOctal%centre
           write(*,*) thisOctal%subcellSize
           write(*,*) thisOctal%phi*radtodeg,thisOctal%dphi*radtodeg
           write(*,*) sqrt(thisOctal%centre%x**2+thisOctal%centre%y**2)
           
!           rVec = subcellCentre(thisOctal,subcell)
!           write(*,*) rVec%x+thisOctal%subcellSize/2.
!           write(*,*) rVec%x-thisOctal%subcellSize/2.
!           write(*,*) rVec%y+thisOctal%subcellSize/2.
!           write(*,*) rVec%y-thisOctal%subcellSize/2.
!           write(*,*) rVec%z+thisOctal%subcellSize/2.
!           write(*,*) rVec%z-thisOctal%subcellSize/2.
!           do ; enddo
!           STOP
           return
        endif
        
        IF ( thisOctal%nDepth /= 1 ) THEN
           thisOctal => thisOctal%parent
        ENDIF

        CALL findSubcellLocalPrivateLevel(point,thisOctal,subcell,haveDescended,boundaryProblem,nDepth)
       
      END IF    
    
    END SUBROUTINE findSubcellLocalPrivateLevel

  END SUBROUTINE findSubcellLocalLevel
   




  FUNCTION decideSplit(thisOctal,subcell,amrLimitScalar,amrLimitScalar2,grid, splitInAzimuth,&
       stellar_cluster,romData) RESULT(split)
    ! returns true if the current voxel is to be subdivided. 
    ! decision is made by comparing 'amrLimitScalar' to some value
    !   derived from information in the current cell  

    use input_variables, only: height, betadisc, rheight, flaringpower, rinner, router, hydrodynamics
    use input_variables, only: drInner, drOuter, rStellar, cavangle, erInner, erOuter, rCore
    use input_variables, only: warpFracHeight, warpRadius, warpSigma, warpAngle
    use input_variables, only: solveVerticalHydro, hydroWarp, rsmooth
    use input_variables, only: rGap, gapWidth, rStar1, rStar2, mass1, mass2, binarysep, mindepthamr, maxdepthamr, vturbmultiplier
    use input_variables, only: planetgap, heightSplitFac, refineCentre, h21cm
    use luc_cir3d_class, only: get_dble_param, cir3d_data
    use cmfgen_class,    only: get_cmfgen_data_array, get_cmfgen_nd, get_cmfgen_Rmin
    use romanova_class, only:  romanova_density
    USE cluster_class, only:   find_n_particle_in_subcell
    use mpi_global_mod, only:  nThreadsGlobal, myRankGlobal

    IMPLICIT NONE

    TYPE(octal), intent(inout) :: thisOctal

!    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    LOGICAL, INTENT(INOUT) :: splitInAzimuth
    real(double), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 ! used for split decision
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(cluster), OPTIONAL, intent(in) :: stellar_cluster
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
    !
    LOGICAL                    :: split        
    real(double) :: massratio, d1, d2
    real(oct)  :: cellSize
    TYPE(vector)     :: searchPoint, rVec
    TYPE(vector)     :: cellCentre
    REAL                  :: x, y, z
    REAL(double) :: hr, rd, fac, warpHeight, phi
    INTEGER               :: i
    real(double)      :: total_mass
    real(double), save :: rgrid(1000)
    real(double)      :: ave_density,  r, dr
    INTEGER               :: nr, nr1, nr2
    real(double)          :: minDensity, maxDensity
    INTEGER               :: nsample = 400
    INTEGER               :: nparticle, limit, npt_subcell
!    real(double) :: timenow
    real(double) :: dummyDouble
    real(double),save  :: R_tmp(204)  ! [10^10cm]
    real(double),allocatable, save  :: R_cmfgen(:)  ! [10^10cm]
    real(double),save  :: Rmin_cmfgen  ! [10^10cm]
    real(double) :: rho
    real(double) :: OstrikerRho(2), r0
    logical, save :: first_time=.true.
    logical :: close_to_star
    real(double)      :: thisScale
    real(double) :: h0
    logical,save  :: firstTime = .true.

    type(VECTOR) :: minV, maxV
    real(double) :: vgradx, vgrady, vgradz, vgrad
    real(double) :: T, vturb

    splitInAzimuth = .false.
    split = .false.
    
   select case(grid%geometry)

    case("melvin")

      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      if (thisOctal%nDepth < 4) split = .true.
      if ((cellCentre%x  < 2.e6).and.(cellSize > 1.e5)) split = .true.

    case("whitney")

       cellSize = thisOctal%subcellSize * 1.d10
       cellCentre = 1.d10 * subcellCentre(thisOctal,subCell)
       nr1 = 50
       nr2 = 10
       nr = nr1 + nr2

      do i = 1, nr1
         rgrid(i) = log10(0.5*drInner)+dble(i)*(log10(drOuter)-log10(0.5*drInner))/dble(nr1)
      end do
      do i = 1, nr2
         rgrid(nr1+i) = log10(drOuter)+dble(i)*(log10(erOuter)-log10(drInner))/dble(nr2)
      end do
      rgrid(1:nr) = 10.d0**rgrid(1:nr)
      r = modulus(cellcentre)
      if (thisOctal%nDepth < 5) split = .true.
      if ((r < rGrid(nr)).and.(r > rGrid(1))) then
         call locate(rGrid, nr, r, i)      
         if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
      endif
       r = sqrt(cellcentre%x**2 + cellcentre%y**2)
       if ((r > drInner*0.9).and.(r < drOuter)) then
          hr = 0.01 * rStellar * (r / rStellar)**1.25
          if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 1.)) split = .true.
          if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 1.)) split = .true.
       endif
       dr = tan(cavAngle/2.) * abs(cellCentre%z)
       if ( ((abs(cellCentre%x) - cellsize/2.) < dr).and.(cellSize > dr/4.) .and.(abs(cellCentre%z)>erInner)) then
          split = .true.
       endif


    case("ttauri","jets","spiralwind","romanova")
      nsample = 100
      ! the density is only sampled at the centre of the grid
      ! we will search in each subcell to see if any point exceeds the 
      ! threshold density

      ! get the size and centre of the current cell
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
    
      ! check the density of random points in the current cell - 
      !   if any of them exceed the critical density, set the flag
      !   indicating the cell is to be split, and return from the function

      ave_density = 0.0_db
      minDensity = 1.e30
      maxDensity = -1.e30
      DO i = 1, nsample
        searchPoint = cellCentre
        CALL RANDOM_NUMBER(x)
        CALL RANDOM_NUMBER(z)
        if (thisOctal%threed) then
           CALL RANDOM_NUMBER(y)
           searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=oct)
           searchPoint%y = searchPoint%y - (cellSize / 2.0_oc) + cellSize*REAL(y,KIND=oct) 
           searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=oct) 
        else
           searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=oct) 
           searchPoint%y = 0.d0
           searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=oct) 
        endif
        if (grid%geometry=="romanova") then
           rho = romanova_density(romData, searchPoint)
        else
           ! using a generic density function in density_mod.f90
           rho = density(searchPoint,grid)
        end if
        ave_density  =  rho + ave_density
        minDensity = min(minDensity, rho)
        maxDensity = max(maxDensity, rho)
      END DO

      ave_density = ave_density / REAL(nsample,KIND=double)


      if (thisOctal%threed) then
         total_mass = maxDensity * (cellSize*1.e10_db)**3
      else
         total_mass = maxDensity * pi * ((searchPoint%x+cellsize/2.)**2-(searchPoint%x-cellsize/2.)**2)*cellsize*1.d30
      endif


      ! weigting toward the smaller radial positions
      if (grid%geometry=="romanova") then
         thisScale = grid%rstar1/modulus(cellCentre)
         total_mass = total_mass * thisScale**3
      elseif (grid%geometry=="ttauri") then
         thisScale = grid%rstar1/modulus(cellCentre)
         total_mass = total_mass * thisScale**4
      end if


      IF (total_mass > amrLimitScalar) then
         split = .TRUE.
      END IF



      ! If 2D and uses large root cell special care must be taken
      if (grid%geometry == "ttauri" .and. .not. thisOctal%threeD) then         
!         if (cellSize > 40.d0  .and.  &
!              (subcell == 1 .or. subcell == 3) )  split=.true.
         close_to_star = .false.
         r = modulus(cellcentre)
         if (r < 300.d0) close_to_star =.true.
         if (close_to_star) then
            if (cellSize > 20.d0) split=.true.
         else
            ! Splits the two cells closer to the origin.
            if (r < cellSize*2.0d0)  split = .true.
         end if

         !
         ! We use a slightly higher resolution near the edge of the accretion flow
!         if (ttau_fuzzy_edge) then
            if (in_fuzzy_edge(cellCentre) ) then
               if (total_mass*4.0d0 > amrLimitScalar) then
               !               if (total_mass*20.0d0 > amrLimitScalar) then  ! for Halpha106
               !               if (total_mass*40.0d0 > amrLimitScalar) then  ! for Halpha107
                  !           ^^^^^ Note the extra factor here
                  split = .TRUE.
               end if
            end if
!         end if

      end if

      if (grid%geometry == "romanova") then         
         if (cellSize > 100.d0) split=.true.
      end if


   case("lexington")
      if (thisOctal%nDepth < mindepthamr) then
         split = .true.
      else
         split = .false.
      endif
!      if (modulus(subcellCentre(thisoctal,subcell)) < 60.*grid%rinner) then
!         if (thisOctal%nDepth < 7) split = .true.
!      endif
!      if (modulus(subcellCentre(thisoctal,subcell)) < 6.*grid%rinner) then
!         if (thisOctal%nDepth < 7) split = .true.
!      endif

   case("starburst")
      if (thisOctal%nDepth < mindepthamr) then
         split = .true.
      else
         split = .false.
      endif

   case("symbiotic")
      if (thisOctal%nDepth < 5) then
         split = .true.
      else
         split = .false.
      endif

   case("gammavel")

      if (thisOctal%nDepth < 6) then
         split = .true.
      else
         split = .false.
      endif
      massRatio = mass1/mass2
      
      d1 = binarySep * (1./(massRatio+1.))
      d2 = binarySep - d1

      rVec = subcellCentre(thisOctal,subcell)

      do i = 1, 100
         rgrid(i) = log10(0.9d0) + (log10(100.d0) - log10(0.9d0))*dble(i-1)/99.d0
      enddo
      rgrid(1:100) = 10.d0**rgrid(1:100)

      r = modulus(rVec - VECTOR(0.d0,0.d0,-d1))/rstar1
      if ((r > rgrid(1)).and.(r<rgrid(99))) then
         call locate(rgrid,100,r,i)
         if (thisOctal%subcellsize/rstar1 > (rgrid(i+1)-rgrid(i))) split = .true.
      endif

      r = modulus(rVec - VECTOR(0.d0,0.d0,d2))/rstar2
      if ((r > rgrid(1)).and.(r<rgrid(99))) then
         call locate(rgrid,100,r,i)
         if (thisOctal%subcellsize/rstar2 > (rgrid(i+1)-rgrid(i))) split = .true.
      endif



   case ("testamr","proto")
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)

      nr1 = 8
      nr2 = 100
      rgrid(1) = 0.8
      rgrid(2) = 0.9
      rgrid(3) = 0.999
      rGrid(4) = 1.
      rGrid(5) = 1.001
      rGrid(6) = 1.002
      rGrid(7) = 1.004
      rGrid(8) = 1.008
      rGrid(1:nr1) = log10(rGrid(1:nr1)*grid%rInner)
      nr = nr1 + nr2
!      do i = 1, nr1
!         rgrid(i) = log10(grid%rInner)+dble(i-1)*(log10(4.*grid%rInner)-log10(grid%rInner))/dble(nr1-1)
!      end do
      do i = 1, nr2
         rgrid(nr1+i) = log10(1.01*grid%rInner)+dble(i)*(log10(grid%rOuter)-log10(1.01*grid%rInner))/dble(nr2)
      end do
      rgrid(1:nr) = 10.d0**rgrid(1:nr)
      r = modulus(cellcentre)
      if (thisOctal%nDepth < 4) split = .true.
      if ((r-cellsize/2.) < grid%rOuter) then
         call locate(rGrid, nr, r, i)      
         if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
      endif

!      if (thisOctal%nDepth > 3) split = .false.
   case ("wrshell")
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)

      nr1 = 8
      nr2 = 100
      rgrid(1) = 0.8
      rgrid(2) = 0.9
      rgrid(3) = 0.999
      rGrid(4) = 1.
      rGrid(5) = 1.001
      rGrid(6) = 1.002
      rGrid(7) = 1.004
      rGrid(8) = 1.008
      rGrid(1:nr1) = log10(rGrid(1:nr1)*grid%rInner)
      nr = nr1 + nr2
      do i = 1, nr2
         rgrid(nr1+i) = log10(1.01*grid%rInner)+dble(i)*(log10(grid%rOuter)-log10(1.01*grid%rInner))/dble(nr2)
      end do
      rgrid(1:nr) = 10.d0**rgrid(1:nr)
      r = modulus(cellcentre)-cellsize/2.
      if (thisOctal%nDepth < 4) split = .true.
      call locate(rGrid, nr, r, i)      
      if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
      if ( (r+cellsize) < rgrid(1)) split = .false.

   case("hydro1d",  "bonnor", "unisphere")

      rVec = subcellCentre(thisOctal, subcell)
      if (thisOctal%nDepth < minDepthAMR) split = .true.
!      rVec = subcellCentre(thisOctal, subcell)
!      if ( (modulus(rVec)< 0.1d0).and.(thisOctal%nDepth < maxDepthAMR) ) split = .true.
!      if ((rVec%x > 0.d0).and.(thisOctal%nDepth < maxDepthAMR) ) split = .true.

   case("kelvin")
      rVec = subcellCentre(thisOctal, subcell)

      if ((abs(rVec%x) < 0.3d0).and.(abs(rvec%z) < 0.3d0)) then
         if (thisOctal%nDepth < (maxDepthAMR)) then
            split = .true.
         endif
      endif

   case("sedov")

      ! Split is decided using mindepthAMR defined globally

   case("protobin")

      ! Split is decided using mindepthAMR defined globally

!      if ((rVec%z > 0).and.(thisOctal%nDepth < maxDepthAMR))  split = .true.
!      if ( (modulus(rVec)< 0.5d5).and.(thisOctal%nDepth < 10) ) split = .true.

   case("benchmark")

      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
!      if (thisOctal%nDepth < 6) split = .true.
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      rd = grid%rOuter / 2. 
      hr = height * (r/rd)**1.125

!      if (.not.thisOctal%cylindrical) then
         if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.2)) split = .true.
         if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 0.2)) split = .true.
!      else
!         if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.2)) split = .true.
!      endif

      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 91.)) then
         splitInAzimuth = .true.
      endif
      if ((r+cellsize/2.d0) < grid%rinner) split = .false.
      if ((r-cellsize/2.d0) > rOuter) split = .false.

   case("molebench")
      cellCentre = subcellCentre(thisOctal, subcell)

      nr = 50
      if (firsttime) then
         open(31, file="model_1.dat", status="old", form="formatted")
         do i = nr,1,-1
            read(31,*) rgrid(i)
         enddo
         rgrid = rgrid * 1.e-10
         close(31)
         firsttime = .false.
      endif
      rd = modulus(cellCentre)
      call locate(rgrid, nr, rd, i)

      if (thisOctal%subcellSize > (rgrid(i+1)-rgrid(i))) split = .true.
      if (rd+0.5d0*thisOctal%subcellSize < rgrid(1)) split = .false.
      if (rd-0.5d0*thisOctal%subcellSize > rgrid(nr)) split = .false.

   case("molefil")

      r0 = (1d3/sqrt(4.d0*pi*6.67d-11*5d-18)) / 1d8 ! 1d3 sigma in m/s, 6.67d-11 is G, 5d-18 is RHOc, 1d8 metres to TORUS

      cellCentre = subcellCentre(thisOctal, subcell)
      rd = modulus(VECTOR(cellCentre%x, cellCentre%y, 0.d0))
         ! change the parameter
      rd = rd+thisOctal%subcellSize/2.d0
      OstrikerRho(1) = 5d-18/(1 + rd**2/(8.d0*r0**2))**2
      rd = rd-thisOctal%subcellSize
      OstrikerRho(2) = 5d-18/(1 + rd**2/(8.d0*r0**2))**2
      if(abs((OstrikerRho(1) - OstrikerRho(2))/5d-18) .ge. 0.02) split = .true.

   case("h2obench1")
      cellCentre = subcellCentre(thisOctal, subcell)
      ! mindepthamr replaces value of 8 in line below
      if (thisOctal%nDepth < mindepthamr) split = .true.
      nr = 200
      if (firsttime) then
         open(31, file="grid.dat", status="old", form="formatted")
         do i = 1,nr
            read(31,*) rgrid(i)
         enddo
         rgrid = rgrid * 3.08568025e8
         close(31)
         firsttime = .false.
      endif
      rd = modulus(cellCentre)
      call locate(rgrid, nr, rd, i)
!      if (thisOctal%subcellSize > (rgrid(i+1)-rgrid(i))) split = .true.
      if (rd+0.5d0*thisOctal%subcellSize < rgrid(1)) split = .false.
      if (rd-0.5d0*thisOctal%subcellSize > rgrid(nr)) split = .false.

   case("h2obench2")
      cellCentre = subcellCentre(thisOctal, subcell)
      ! mindepthamr replaces value of 8 in line below
      if (thisOctal%nDepth < mindepthamr) split = .true.
      nr = 200
      if (firsttime) then
         open(31, file="grid.dat", status="old", form="formatted")
         do i = 1,nr
            read(31,*) rgrid(i)
         enddo
         rgrid = rgrid * 3.08568025e8
         close(31)
         firsttime = .false.
      endif
      rd = modulus(cellCentre)
      call locate(rgrid, nr, rd, i)
!      if (thisOctal%subcellSize > (rgrid(i+1)-rgrid(i))) split = .true.
      if (rd+0.5d0*thisOctal%subcellSize < rgrid(1)) split = .false.
      if (rd-0.5d0*thisOctal%subcellSize > rgrid(nr)) split = .false.

   case("agbstar")
      cellCentre = subcellCentre(thisOctal, subcell)
      ! mindepthamr replaces 4
      if (thisOctal%nDepth < mindepthamr) split = .true.
      nr = 100
      if (firsttime) then
         open(31, file="mc_100.dat", status="old", form="formatted")
         do i = nr,1,-1
            read(31,*) rgrid(i)
         enddo
         close(31)
         rgrid = rgrid * 1e-10
         firsttime = .false.
      endif

      rd = modulus(cellCentre)
      call locate(rgrid, nr, rd, i)
      if (2. * thisOctal%subcellSize > (rgrid(i+1)-rgrid(i))) split = .true.
      if (rd+0.5d0*thisOctal%subcellSize < rgrid(1)) split = .false.
      if (rd-0.5d0*thisOctal%subcellSize > rgrid(nr)) split = .false.
                                                          
   case("luc_cir3d") 
      if (first_time) then
         open(unit=77, file ="zeus_rgrid.dat", status="old")
         do i = 1, 204
            read(77,*) R_tmp(i)
         end do
         close(77)
         R_tmp(:) = R_tmp(:) * get_dble_param(cir3d_data,"Rs") ! [10^10cm]
         first_time=.false.         
      end if
      cellSize = thisOctal%subcellSize
      cellCentre = subcellCentre(thisOctal,subCell)
      r = modulus(cellcentre)
      nr=204
      call locate(R_tmp, nr, r, i)      
      if (i == 0) i = nr-1
      if (i == nr) i = nr-1
      if (cellsize*amrlimitscalar > (R_tmp(i+1)-R_tmp(i))) then
         split = .true.
      else
         split = .false.
      end if

      if (cellSize > 100.0d0)  split=.true.

   case("cmfgen") 
      nr = get_cmfgen_nd()
      if (first_time) then
         ! retriving the r grid in CMFGEN data.
         ALLOCATE(R_cmfgen(nr))
         call get_cmfgen_data_array("R", R_cmfgen) ! [10^10cm]
         Rmin_cmfgen = get_cmfgen_Rmin()
         first_time=.false.         
      end if
      cellSize = thisOctal%subcellSize
      cellCentre = subcellCentre(thisOctal,subCell)
      r = modulus(cellcentre)

      if (r > R_cmfgen(nr) .or. r < Rmin_cmfgen ) then
         split = .false.
      else
         call locate(R_cmfgen, nr, r, i)      
         if (i == 0) i = nr-1
         if (i == nr) i = nr-1

         dR = ABS(R_cmfgen(i+1)-R_cmfgen(i))
         thisScale = cellsize/amrlimitscalar
         
         if ( thisScale > dR) then
            split = .true.
         else
            split = .false.
         end if
      end if

      if (cellSize > ABS(R_cmfgen(nr)-Rmin_cmfgen)/4.0d0)  split=.true.

   case ("cluster","molcluster","theGalaxy")

      call find_n_particle_in_subcell(nparticle, ave_density, &
           thisOctal, subcell, rho_min=minDensity, rho_max=maxDensity, n_pt=npt_subcell, v_min=minV, v_max=maxV)

      total_mass = cellVolume(thisOctal, subcell)  * 1.d30

      thisOctal%rho(subcell) = ((maxdensity * total_mass) / mindensity)**(1.d0/3.d0) ! placeholder for maximum expected smoothing length (not rho!) 

      total_mass = ave_density * total_mass

      if (total_mass > amrlimitscalar) then
         split = .true.
         mass_split = mass_split + 1
      endif

      ! Split in order to capture density gradients. 

! Use sign of amrlimitscalar2 to determine which condition to use.
! This allows testing of the grid generation method without recompiling the code

      if ( amrlimitscalar2 > 0.0 ) then 
         if ( ( (maxDensity-minDensity) / (maxDensity+minDensity) )  > amrlimitscalar2 ) split =.true.
      elseif (amrlimitscalar2 < 0.0 ) then
         if (  (maxDensity / minDensity) > abs(amrlimitscalar2) ) then 
            if(split) then
               both_split = both_split + 1
            else
               density_split = density_split + 1
               split = .true.
            endif
         end if
      endif

!      if(maxdensity .gt. 1d-13 .and. nparticle .gt. 1) then
!         split = .true.
!         maxdensity_split = maxdensity_split + 1
!      endif


      if(.not. split .and. (nparticle .ge. 2) .and. (.not. h21cm) ) then
         if(ave_density .gt. 1d-13) then
            T = 10.d0 * (ave_density * 1d13)**(0.4d0)
         Vturb = 5.d0 * sqrt(2d-10 * kerg * T / (28.d0 * amu) + 0.3**2) / (cspeed * 1d-5) ! 5 is fudge factor to make sure condition isn't too strigent ! 28 is mass of CO
            Vturb = sqrt(5.938e-4 * T + 0.09) / (cspeed * 1d-5)
         else
            Vturb = 1.03246E-06 ! above calculation with T = 10
         endif

         vgradx = maxV%x - minV%x
         vgrady = maxV%y - minV%y
         vgradz = maxV%z - minV%z
         
         vgrad = max(vgradx,max(vgrady,vgradz))

         if(vgrad .gt. vturbmultiplier * vturb) then
            velocity_split = velocity_split + 1
            split = .true.
         endif
      endif


      !Jeans mass condition
!      if(.not. split .and. )then
!         split = .true.
!         mass_split2 = mass_split2 + 1
!      endif
      
! Additional refinement at the grid centre used for SPH-Torus discs. 
      if ( refineCentre ) then 

         cellSize   = thisOctal%subcellSize
         cellCentre = subcellCentre(thisOctal,subCell)

         if ( abs(cellCentre%x) < 10.0*grid%rCore .and.  abs(cellCentre%y) < 10.0*grid%rCore .and.  abs(cellCentre%z) &
              < 40.0*grid%rCore ) then 
            if ( cellSize > grid%rCore ) split = .true. 
         end if

         if ( abs(cellCentre%x) < 2.5*grid%rCore .and.  abs(cellCentre%y) < 2.5*grid%rCore .and.  & 
              abs(cellCentre%z) < 20.0*grid%rCore ) then 
            if ( cellSize > 0.5*grid%rCore ) split = .true. 
         end if

      end if

! The stellar disc code is retained in case this functionality needs to be used in future 
!!$
!!$      if (include_disc(stellar_cluster)) then
!!$      
!!$         ! If the subcell intersects with the stellar disk.. then do additional check.
!!$         ! This time, the cells will be split when the mass of the cell exceeds a 
!!$         ! critical mass specified by "amrlimitscalar.
!!$
!!$         if (stellar_disc_exists(sphData) .and.  &       
!!$              disc_intersects_subcell(stellar_cluster, sphData, thisOctal, subcell) ) then
!!$         rho_disc_ave = average_disc_density_fast(sphData, thisOctal, &
!!$              subcell, stellar_cluster, scale_length)
!!$
!!$!!         rho_disc_ave = average_disc_density(sphData, thisOctal, &
!!$!!              subcell, stellar_cluster, scale_length)
!!$!            rho_disc_ave = max_disc_density_from_sample(sphData, thisOctal, &
!!$!                 subcell, stellar_cluster, scale_length)
!!$
!!$            total_mass = total_mass + rho_disc_ave * (cellSize*1.e10_db)**3  !  [g]
!!$
!!$            if (cellSize > 1.0e4  .or.  &
!!$                 (cellsize > scale_length) ) then
!!$!                 (cellsize > scale_length .and. rho_disc_ave > 1.0e-20) ) then
!!$               split = .true.  ! units in 10^10cm
!!$            end if
!!$            
!!$         end if
!!$
!!$      end if

   case ("wr104")
      ! Splits if the number of particle is more than a critical value (~3).
      limit = nint(amrLimitScalar)
      
      call find_n_particle_in_subcell(nparticle, dummyDouble, thisOctal, subcell)

      rVec = subcellCentre(thisOctal,subcell)
      r = sqrt(rVec%x**2 + rVec%y**2)

      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 89.)) then
         splitInAzimuth = .true.
         limit = nint(amrLimitScalar)*1000
      endif

      if ((r < 2.e6).and.(thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 31.)) then
         limit = nint(amrLimitScalar)*100
         splitInAzimuth = .true.
      endif

      if ((r < 2.e5).and.(thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 15)) then
         splitInAzimuth = .true.
         limit = nint(amrLimitScalar)*1

      endif



      !
      if (nparticle > limit ) then 
	 split = .TRUE.
      else
         split = .FALSE.
      end if
      cellCentre = subcellCentre(thisOctal,subCell)



!      write(*,*) nparticle,thisOctal%nDepth,subcell

   case ("windtest")

      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)

      nr = 25
      r  = 1. 
      dr = 1.e-3
      do i = 1, nr, 1
        rgrid(i) = r 
        r = r + dr
        dr = dr * 1.25
      end do
      do i = 2, nr, 1
        rgrid(i) = grid%rInner + ((grid%rOuter - grid%rInner) * (rgrid(i)-rgrid(1))/(rgrid(nr) - rGrid(1)))
      end do
      rgrid(1) = grid%rInner
 !    print *, rgrid(1), rgrid(nr) 
!      nr1 = 30
!      nr2 = 0
!      nr = nr1 + nr2
!      do i = 1, nr1
!         rgrid(i) = log10(grid%rInner)+dble(i-1)*(log10(grid%rOuter)-log10(grid%rInner))/dble(nr1-1)
!      end do
!      do i = 1, nr2
!         rgrid(nr1+i) = log10(0.01*grid%rOuter)+dble(i)*(log10(grid%rOuter)-log10(0.01*grid%rOuter))/dble(nr2)
!      end do
      
!      rgrid(1:nr) = 10.d0**rgrid(1:nr)
      r = modulus(cellcentre)
      if ((r < grid%rOuter).and.(r > grid%rinner)) then
         call locate(rGrid, nr, r, i)      
         if (i == 0) i = nr-1
         if (i == nr) i = nr-1
         if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
      endif
      if (thisOctal%nDepth < 4) split = .true.

   case("ggtau")
! used to be 5
      if (thisOctal%ndepth < mindepthamr) split = .true.

      cellCentre = subcellCentre(thisOctal,subCell)

      r = sqrt(cellcentre%x**2+cellcentre%y**2)
      r = r * 6.68458134e-06! (torus units to 100's of AUs)

      cellcentre%z = cellcentre%z * 6.68458134e-04
      cellSize = thisOctal%subcellSize * 6.68458134e-4

      H0 = 14.55
      hr = H0 * (r ** (1.25))

      if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 0.1)) then
         split = .true.
         scaleheighta_count = scaleheighta_count + 1
         goto 1001
      endif

      if ((abs(cellcentre%z)/hr > 5.) .and. (abs(cellcentre%z/cellsize) < 5.)) then
         split = .true.
          scaleheightb_count =  scaleheightb_count + 1
         goto 1001
      endif

      if ((r .gt. 1.8*0.99).and.(r .lt. 1.01*1.8)) then
         if ((abs(cellcentre%z)/hr < 1.)) then
            if (cellsize > 0.1*1.8) then
               split = .true.
               scaleheightc_count = scaleheightc_count + 1
               goto 1001
            endif
         endif
      endif

1001  if ((r+(cellsize/(2.d0*100.))) < 1.8) split = .false.
      if ((r-(cellsize/(2.d0*100.))) > 8.) split = .false.

   case("shakara","aksco")
 ! used to be 5
      if (thisOctal%ndepth  < mindepthamr) split = .true.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      hr = height * (r / (100.d0*autocm/1.d10))**betadisc

!      if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 1.)) split = .true.

      if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > heightSplitFac)) split = .true.

      if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.

      if (((r-cellsize/2.d0) < grid%rinner).and. ((r+cellsize/2.d0) > grid%rInner) .and. &
           (thisOctal%nDepth < maxdepthamr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.

      if (((r-cellsize/2.d0) < rOuter).and. ((r+cellsize/2.d0) > rOuter) .and. &
           (thisOctal%subcellSize/rOuter > 0.01) .and. (abs(cellCentre%z/hr) < 7.d0) ) split=.true.

      if ((r+cellsize/2.d0) < grid%rinner*1.) split = .false.
      if ((r-cellsize/2.d0) > grid%router*1.) split = .false.

!      if ((r > grid%rinner).and.(r < 1.01d0*grid%rinner)) then
!         if ((abs(cellcentre%z)/hr < 1.)) then
!            if (cellsize > 5.d-1*grid%rinner) split = .true.
!         endif
!      endif
!
!      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 31.)) then
!         splitInAzimuth = .true.
!         split = .true.
!      endif
!
      if ((r > rOuter*1.1d0).and.(thisOctal%nDepth > 4)) then
         split = .false.
         splitInAzimuth = .false.
      endif

      if (hydrodynamics) then
         split = .false.
         if (thisOctal%nDepth < minDepthAMR) split = .true.
      endif

!      if ((r > grid%rinner).and.(r < 1.001d0*grid%rinner)) then
!         if ((abs(cellcentre%z)/hr < 2.)) then
!            split = .true.
!         endif
!      endif

   case("circumbin")

      if (thisOctal%ndepth  < 5) split = .true.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      hr = height * (r / (100.d0*autocm/1.d10))**betadisc

!      if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 1.)) split = .true.
      if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 0.2)) split = .true.

      if ((abs(cellcentre%z)/hr > 2.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.

      if ((r+cellsize/2.d0) < grid%rinner*1.) split = .false.
      if ((r-cellsize/2.d0) > grid%router*1.) split = .false.

      if (((r-cellsize/2.d0) < grid%rinner).and. ((r+cellsize/2.d0) > grid%rInner) .and. &
           (thisOctal%nDepth < maxDepthAmr) .and. (abs(cellCentre%z/hr) < 3.d0) ) split=.true.

!      splitInAzimuth = .false.
!      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 181.)) then
!         splitInAzimuth = .true.
!         split = .true.
!      endif

!      if (abs(cellCentre%z) < rinner/2.) then
!         if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 11.).and.(r < rInner)) then
!            splitInAzimuth = .true.
!            split = .true.
!         endif
!      endif
!      if (abs(cellCentre%z) < rinner/5.) then
!         if ((r < rinner).and.(thisOctal%subcellSize > (0.05*rinner))) split = .true.
!      endif

      if ((r > rOuter*1.1d0).and.(thisOctal%nDepth > 4)) then
         split = .false.
         splitInAzimuth = .false.
      endif


   case("iras04158")

      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      hr = height * (r / (50.d0*autocm*1.d-10))**betadisc

!! This splitting law replaces the commented out code below. It does the same stuff but further out
!! and more gradually. I hope the splitting in r has become redundant because it should be linked to
!! the increased resolution in z.

      if(abs(cellcentre%z) .lt. 6.d0 * hr) then ! 6 is a reasonable number for this law
!         if(cellsize .gt. 2.d0**(abs(cellcentre%z)/hr - 3.d0) * hr) split = .true.
         if(8.d0 * cellsize .gt. 2.d0**(abs(cellcentre%z)/hr) * hr) split = .true.
      endif

      if((abs(cellcentre%z) .gt. 3.d0 * hr).and.(abs(cellcentre%z) < 4.d0 * cellsize)) then
         split = .true.
      endif

!      if ((abs(cellcentre%z) < 1.d0 * hr) .and. (cellsize > 0.25d0 * hr)) then
!         split = .true.
!      elseif ((abs(cellcentre%z) < 2.d0 * hr) .and. (cellsize > 0.5d0 * hr)) then
!         split = .true.
!      elseif((r > grid%rInner).and.(r < grid%rInner * 1.02)) then
!         if ((abs(cellcentre%z) < 3. * hr)) then        
!            if (cellsize > 0.1d0 * grid%rInner) split = .true.
!         endif
!      endif

!      if((r > grid%rInner).and.(r < grid%rInner * 1.1)) then
!         if ((abs(cellcentre%z) < 1.d0 * hr)) then        
!            if (cellsize > grid%rInner) split = .true.
!         endif
!      endif

      if ((r-cellsize*0.5d0) > grid%router*1.005) then
         split = .false.
         goto 555
      elseif ((r+cellsize*0.5d0) < grid%rinner*0.995) then
         split = .false.
         goto 555
      endif

      555 continue

   case("oldshakara")

!      if (thisOctal%ndepth  < 6) split = .true.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      hr = height * (r / (100.d0*autocm/1.d10))**betadisc

      if ((abs(cellcentre%z)/hr < 7.) .and. (cellsize/hr > 0.2)) split = .true.

!      if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 0.2)) split = .true.

      if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.

!      if ((r > grid%rInner).and.(r < grid%rInner * 1.01)) then
!         if ((abs(cellcentre%z)/hr < 5.)) then
!            if (cellsize > 1.e-3 * grid%rInner) split = .true.
!         endif
!      endif



      if ((r+cellsize/2.d0) < grid%rinner*0.9) split = .false.

!      if ((r > grid%rinner).and.(r < 1.01d0*grid%rinner)) then
!         if ((abs(cellcentre%z)/hr < 2.)) then
!            if (cellsize > 1.d-3*grid%rinner) split = .true.
!         endif
!      endif

!      if ((r > grid%rinner).and.(r < 1.001d0*grid%rinner)) then
!         if ((abs(cellcentre%z)/hr < 2.)) then
!            split = .true.
!         endif
!      endif


   case("ppdisk")

      cellSize = thisOctal%subcellSize
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
! The 100 is because Tim's code uses the scale height at 100AU.
!      hr = height * (r / (100.d0*autocm/1.d10))**1.25
      hr = (auTocm * height * rHeight * (r / (rHeight*autocm/1.e10))**flaringPower)/1.e10
      if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 0.5)) split = .true.
      if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
      if ((r+cellsize/2.d0) < rsmooth*autocm/1.e10) split = .false.

   case("planetgap")

      cellSize = thisOctal%subcellSize
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      hr = height * rCore * (r/rCore)**betaDisc
      if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
      if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.5)) split = .true.
      ! The next few lines are used in the Varniere models instead of the line
      ! above. The final line tries to prevent splitting in the denser
      ! midplane of the disc.
!      if (r > 6000.) then
!        if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.5)) split = .true.
!      else
!        if ((abs(cellcentre%z)/hr < 12.) .and. (cellsize/hr > 0.5)) split = .true.
!      end if
!      if (abs(cellcentre%z)/hr < 0.9) split = .false.
      ! Added in to get an extra level of refinement at the outer edge of the
      ! gap. Cells at radii between that of the gap and (gap+gapWidth) must
      ! have a resolution of a quarter of a scale height instead of a half...
      if (planetgap) then
         if ((r > rGap*autocm/1e10) .and. (r < (rGap+1.5*gapWidth)*autocm/1e10) .and.  &
             (abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.25)) split = .true.
        ! The alternative version below was used for the models where we want a
        ! high resolution gap (hi-res-gap models).
!        if ((r > (rGap-gapWidth)*autocm/1e10) .and. (r < (rGap+1.5*gapWidth)*autocm/1e10) .and.  &
!            (abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.1)) split = .true.
      end if
      if ((r+cellsize/2.d0) < 0.9*grid%rinner) split = .false.



   case("clumpydisc")

      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
!      if (thisOctal%nDepth < 4) split = .true.
!      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
!      phi = atan2(cellcentre%y,cellcentre%x)
!      warpHeight = 0. !cos(phi) * grid%rInner * sin(30.*degtorad) * sqrt(grid%rinner / r)
      do i = 1, grid%ng
         if (inOctal(thisOctal,grid%gArray(i)%centre).and.(cellSize > grid%gArray(i)%sigma)) then
            split = .true.
         endif
      enddo
      
!      if (r > grid%rInner*0.8) then
!         hr = height * (r / (100.*autocm/1.e10))**1.25
!         fac = cellsize/hr
!         if (abs((cellCentre%z-warpHeight)/hr) < 10.) then
!            if (fac > 2.) split = .true.
!         endif
!      endif

   case("pathtest")
      if (thisOctal%nDepth < minDepthAMR) split = .true.
      if ((modulus(subcellCentre(thisOctal, subcell)) < 0.2).and. &
           (thisOctal%nDepth<maxDepthAMR)) split = .true.
   case("toruslogo")
!used to be 6
      if (thisOctal%nDepth  < 8) split = .true.

      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 10.)) then
         split = .true.
         splitInAzimuth = .true.
      endif
      cellCentre = subcellCentre(thisOctal, subcell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2 + cellcentre%z**2)
      if ((r + thisOctal%subcellSize/2.d0)  .lt. 1.5*rsol/1.e10) then
         split = .false.
         splitinAzimuth = .false.
      endif
      if ((r - thisOctal%subcellSize/2.d0)  .gt. 2.*rsol/1.e10) then
         split = .false.
         splitinAzimuth = .false.
      endif

!      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 90.)) then
!         split = .true.
!         splitInAzimuth = .true.
!      endif


   case("warpeddisc")

      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      phi = atan2(cellcentre%y,cellcentre%x)
      warpheight  = cos(phi+warpAngle) * warpFracHeight * warpradius * exp(-0.5d0*((r - warpRadius)/warpSigma)**2)


      if ((r+cellsize/2.d0) > rInner*0.8) then
         hr = height * (r / rOuter)**betaDisc
         fac = cellsize/hr
         if (solveVerticalHydro) then
           ! split the grid in a way more similar to planetgap, to give hydro the best chance possible
           if ((abs((cellCentre%z-warpHeight)/hr) < 10.).and.(fac .gt. 0.5)) split = .true.
           if ((abs(cellcentre%z-warpheight)/hr > 5.).and.(abs((cellcentre%z-warpheight)/cellsize) < 2.)) split = .true.
         else
           if (hydroWarp) then
             fac=hr
             hr = getHydroWarpScaleHeight(r, grid)
!print *, fac, hr
             if (.not. (hr > 0.)) then
               hr = fac
!               if (r < rOuter*1.1d0) then
!                 print *, "Error: scaleheight returned by getHydroWarpScaleHeight = 0."
!                 print *, r, hr
!                 stop
!               end if
             end if
             fac = cellsize/hr
             ! o The first criterion of this if (cell is less than 5 scaleheights from midplane) is fine because
             !   it overestimates the height above the midplane that needs to be split.
             ! o The second criterion will not split some cells that should be split (the new scaleheight is smaller,
             !   so fac will be smaller than it should be and may not be > 1). To compensate, we change 1. to 0.75.
             ! If we were to do this properly, we would find out what the new scale height of the disc is...
           end if
           if ((abs((cellCentre%z-warpHeight)/hr) < 5.).and.(fac .gt. 1.)) split = .true.
           if ((abs(cellcentre%z-warpheight)/hr > 5.).and.(abs((cellcentre%z-warpheight)/cellsize) < 1.)) split = .true.
           if ((abs(r - warpradius) < warpsigma).and.(abs((cellCentre%z-warpHeight)/hr) < 8.).and.(fac .gt. 0.5)) split = .true.
         end if
      endif

      if ((abs(r-rinner) < 0.9*rinner).and.(cellSize > 0.02*rInner).and.(abs(cellCentre%z) < 2.*hr)) then
         split = .true.
      endif

      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 30.)) then
         split = .true.
         splitInAzimuth = .true.
      endif
      
      if (r > rOuter*1.1d0) then
         split = .false.
         splitInAzimuth = .false.
      endif

   case("magstream")

      
      IF (thisOctal%nDepth < 5) THEN 
        split = .TRUE.
     endif
!      ELSE
!
!        nsample = 1
!        ! the density is only sampled at the centre of the grid
!        ! we will search in each subcell to see if any point exceeds the 
!        ! threshold density
!
!        ! get the size and centre of the current cell
!        cellSize = thisOctal%subcellSize 
!        cellCentre = subcellCentre(thisOctal,subCell)
!    
!        ! check the density of random points in the current cell - 
!        !   if any of them exceed the critical density, set the flag
!        !   indicating the cell is to be split, and return from the function
!        split = .FALSE.
!        ave_density = 0.0_db
!        DO i = 1, nsample
!          searchPoint = cellCentre
!          CALL RANDOM_NUMBER(x)
!          CALL RANDOM_NUMBER(y)
!          CALL RANDOM_NUMBER(z)
!          searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=oct) 
!          searchPoint%y = searchPoint%y - (cellSize / 2.0_oc) + cellSize*REAL(y,KIND=oct) 
!          searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=oct)
!
!          ! using a generic density function in density_mod.f90
!          ave_density  = density(searchPoint,grid,timeNow) + ave_density
!
!        END DO
!
!        ave_density = ave_density / REAL(nsample,KIND=double)
!        total_mass = ave_density * (cellSize*1.e10_db)**3.0_db
!        IF (total_mass > amrLimitScalar) then
!           split = .TRUE.
!        END IF
!
!      END IF


   case DEFAULT
      PRINT *, 'Invalid grid geometry option passed to amr_mod::decideSplit'
      PRINT *, 'grid%geometry ==', TRIM(grid%geometry)
      PRINT *, 'Exiting the program .... '
      STOP
   end select

   if (thisOctal%nDepth == maxDepthAmr) then
      split = .false.
      if (firstTime) then
         call writeWarning("AMR cell depth capped")
         firstTime = .false.
      endif
   endif

   if (thisOctal%nDepth < minDepthAmr) then
      split = .true.
   endif

   if (grid%splitOverMPI) then
      if ((thisOctal%twoD.and.(nThreadsGlobal-1)==4).or. &
           (thisOctal%threed.and.(nThreadsGlobal-1)==8).or. &
           (thisOctal%oned  .and.(nThreadsGlobal-1)==2) ) then
         if ((thisOctal%nDepth == 1).and.(subcell /= myRankGlobal)) then
            split = .false.
         endif
         if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
            split = .false.
         endif
      else
         if (thisOctal%oneD) then
            if (thisOctal%nDepth == 1) then
               split = .true.
            else if (thisOctal%nDepth > 1) then
               if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                  split = .false.
               endif
            endif
         endif
         if (thisOctal%twoD) then
            if (thisOctal%nDepth == 1) then
               split = .true.
            else if (thisOctal%nDepth > 1) then
               if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                  split = .false.
               endif
            endif
         endif
         if (thisOctal%threeD) then
            if (thisOctal%nDepth == 1) then
               split = .true.
            else if (thisOctal%nDepth > 1) then
               if (thisOctal%mpiThread(subcell) /= myRankGlobal) then
                  split = .false.
               endif
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
   endif

  END FUNCTION decideSplit


  PURE SUBROUTINE fillGridDummyValues(thisOctal,subcell, grid) 
    ! for testing, just put some generic values in an octal

    IMPLICIT NONE
    
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(gridtype), intent(in) :: grid

    if (.not.grid%oneKappa) then
       thisOctal%kappaAbs(subcell,:) = 1.0e-20
       thisOctal%kappaSca(subcell,:) = 1.0e-20 
    endif
    thisOctal%chiLine(subcell)    = 1.0e-20
    thisOctal%etaLine(subcell)    = 1.0e-20 
    thisOctal%etaCont(subcell)    = 1.0e-20 
    thisOctal%N(subcell,:)        = 1.0e-20 
    thisOctal%Ne(subcell)         = 1.0e-20 
    thisOctal%nTot(subcell)       = 1.0e-20
    thisOctal%biasLine3D(subcell) = 1.0 
    thisOctal%biasCont3D(subcell) = 1.0 

  END SUBROUTINE fillGridDummyValues
  
  
  SUBROUTINE fillVelocityCorners(thisOctal,grid,velocityFunc, threed)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
  
    TYPE(octal), INTENT(INOUT) :: thisOctal
    TYPE(gridtype), INTENT(IN) :: grid
    LOGICAL, INTENT(IN) :: threed
    real(oct)      :: r1, r2, r3
    real(oct)      :: phi1, phi2, phi3
    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    
    INTERFACE 
      TYPE(vector) FUNCTION velocityFunc(point,grid)
        USE vector_mod
        USE gridtype_mod
        TYPE(vector), INTENT(IN) :: point
        TYPE(gridtype), INTENT(IN)    :: grid
      END FUNCTION velocityFunc
    END INTERFACE

    if (thisOctal%oneD) then
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       y1 = 0.d0
       z1 = 0.d0
       thisOctal%cornerVelocity(1) = velocityFunc(vector(x1,y1,z1),grid)
       thisOctal%cornerVelocity(2) = velocityFunc(vector(x2,y1,z1),grid)
       thisOctal%cornerVelocity(3) = velocityFunc(vector(x3,y1,z1),grid)
       goto 666
    endif


    if (thisOctal%threed) then
       if (.not.thisOctal%cylindrical) then ! 3d cartesian case
          ! we first store the values we use to assemble the position vectors
          
          x1 = thisOctal%centre%x - thisOctal%subcellSize
          x2 = thisOctal%centre%x
          x3 = thisOctal%centre%x + thisOctal%subcellSize
          
          y1 = thisOctal%centre%y - thisOctal%subcellSize
          y2 = thisOctal%centre%y
          y3 = thisOctal%centre%y + thisOctal%subcellSize
          
          z1 = thisOctal%centre%z - thisOctal%subcellSize
          z2 = thisOctal%centre%z
          z3 = thisOctal%centre%z + thisOctal%subcellSize
                    
          ! now store the 'base level' values
          
          thisOctal%cornerVelocity(1) = velocityFunc(vector(x1,y1,z1),grid)
          thisOctal%cornerVelocity(2) = velocityFunc(vector(x2,y1,z1),grid)
          thisOctal%cornerVelocity(3) = velocityFunc(vector(x3,y1,z1),grid)
          thisOctal%cornerVelocity(4) = velocityFunc(vector(x1,y2,z1),grid)
          thisOctal%cornerVelocity(5) = velocityFunc(vector(x2,y2,z1),grid)
          thisOctal%cornerVelocity(6) = velocityFunc(vector(x3,y2,z1),grid)
          thisOctal%cornerVelocity(7) = velocityFunc(vector(x1,y3,z1),grid)
          thisOctal%cornerVelocity(8) = velocityFunc(vector(x2,y3,z1),grid)
          thisOctal%cornerVelocity(9) = velocityFunc(vector(x3,y3,z1),grid)
          
          ! middle level
          
          thisOctal%cornerVelocity(10) = velocityFunc(vector(x1,y1,z2),grid)
          thisOctal%cornerVelocity(11) = velocityFunc(vector(x2,y1,z2),grid)
          thisOctal%cornerVelocity(12) = velocityFunc(vector(x3,y1,z2),grid)
          thisOctal%cornerVelocity(13) = velocityFunc(vector(x1,y2,z2),grid)
          thisOctal%cornerVelocity(14) = velocityFunc(vector(x2,y2,z2),grid)
          thisOctal%cornerVelocity(15) = velocityFunc(vector(x3,y2,z2),grid)
          thisOctal%cornerVelocity(16) = velocityFunc(vector(x1,y3,z2),grid)
          thisOctal%cornerVelocity(17) = velocityFunc(vector(x2,y3,z2),grid)
          thisOctal%cornerVelocity(18) = velocityFunc(vector(x3,y3,z2),grid)
          
          ! top level
          
          thisOctal%cornerVelocity(19) = velocityFunc(vector(x1,y1,z3),grid)
          thisOctal%cornerVelocity(20) = velocityFunc(vector(x2,y1,z3),grid)
          thisOctal%cornerVelocity(21) = velocityFunc(vector(x3,y1,z3),grid)
          thisOctal%cornerVelocity(22) = velocityFunc(vector(x1,y2,z3),grid)
          thisOctal%cornerVelocity(23) = velocityFunc(vector(x2,y2,z3),grid)
          thisOctal%cornerVelocity(24) = velocityFunc(vector(x3,y2,z3),grid)
          thisOctal%cornerVelocity(25) = velocityFunc(vector(x1,y3,z3),grid)
          thisOctal%cornerVelocity(26) = velocityFunc(vector(x2,y3,z3),grid)
          thisOctal%cornerVelocity(27) = velocityFunc(vector(x3,y3,z3),grid)

       else ! cylindrical 
          if (thisOctal%splitAzimuthally) then
             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi 
             phi3 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize

             ! bottom level

             thisOctal%cornerVelocity(1) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1),grid)
             thisOctal%cornerVelocity(2) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z1),grid)
             thisOctal%cornerVelocity(3) = velocityFunc(vector(r1*cos(phi3),r1*sin(phi3),z1),grid)
             thisOctal%cornerVelocity(4) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z1),grid)
             thisOctal%cornerVelocity(5) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z1),grid)
             thisOctal%cornerVelocity(6) = velocityFunc(vector(r2*cos(phi3),r2*sin(phi3),z1),grid)
             thisOctal%cornerVelocity(7) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z1),grid)
             thisOctal%cornerVelocity(8) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1),grid)
             thisOctal%cornerVelocity(9) = velocityFunc(vector(r3*cos(phi3),r3*sin(phi3),z1),grid)

             ! middle level

             thisOctal%cornerVelocity(10) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2),grid)
             thisOctal%cornerVelocity(11) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z2),grid)
             thisOctal%cornerVelocity(12) = velocityFunc(vector(r1*cos(phi3),r1*sin(phi3),z2),grid)
             thisOctal%cornerVelocity(13) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z2),grid)
             thisOctal%cornerVelocity(14) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z2),grid)
             thisOctal%cornerVelocity(15) = velocityFunc(vector(r2*cos(phi3),r2*sin(phi3),z2),grid)
             thisOctal%cornerVelocity(16) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z2),grid)
             thisOctal%cornerVelocity(17) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2),grid)
             thisOctal%cornerVelocity(18) = velocityFunc(vector(r3*cos(phi3),r3*sin(phi3),z2),grid)

             ! top level

             thisOctal%cornerVelocity(10) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3),grid)
             thisOctal%cornerVelocity(11) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z3),grid)
             thisOctal%cornerVelocity(12) = velocityFunc(vector(r1*cos(phi3),r1*sin(phi3),z3),grid)
             thisOctal%cornerVelocity(13) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z3),grid)
             thisOctal%cornerVelocity(14) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z3),grid)
             thisOctal%cornerVelocity(15) = velocityFunc(vector(r2*cos(phi3),r2*sin(phi3),z3),grid)
             thisOctal%cornerVelocity(16) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z3),grid)
             thisOctal%cornerVelocity(17) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3),grid)
             thisOctal%cornerVelocity(18) = velocityFunc(vector(r3*cos(phi3),r3*sin(phi3),z3),grid)

          else

             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize


             ! bottom level

             thisOctal%cornerVelocity(1) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1),grid)
             thisOctal%cornerVelocity(2) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z1),grid)
             thisOctal%cornerVelocity(3) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z1),grid)
             thisOctal%cornerVelocity(4) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z1),grid)
             thisOctal%cornerVelocity(5) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z1),grid)
             thisOctal%cornerVelocity(6) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1),grid)

             ! middle level

             thisOctal%cornerVelocity(7) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2),grid)
             thisOctal%cornerVelocity(8) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z2),grid)
             thisOctal%cornerVelocity(9) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z2),grid)
             thisOctal%cornerVelocity(10) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z2),grid)
             thisOctal%cornerVelocity(11) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z2),grid)
             thisOctal%cornerVelocity(12) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2),grid)

             ! top level

             thisOctal%cornerVelocity(13) = velocityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3),grid)
             thisOctal%cornerVelocity(14) = velocityFunc(vector(r1*cos(phi2),r1*sin(phi2),z3),grid)
             thisOctal%cornerVelocity(15) = velocityFunc(vector(r2*cos(phi1),r2*sin(phi1),z3),grid)
             thisOctal%cornerVelocity(16) = velocityFunc(vector(r2*cos(phi2),r2*sin(phi2),z3),grid)
             thisOctal%cornerVelocity(17) = velocityFunc(vector(r3*cos(phi1),r3*sin(phi1),z3),grid)
             thisOctal%cornerVelocity(18) = velocityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3),grid)

          endif
       endif
    else
       
       
    ! we first store the values we use to assemble the position vectors
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values
       
       thisOctal%cornerVelocity(1) = velocityFunc(vector(x1,0.d0,z1),grid)
       thisOctal%cornerVelocity(2) = velocityFunc(vector(x2,0.d0,z1),grid)
       thisOctal%cornerVelocity(3) = velocityFunc(vector(x3,0.d0,z1),grid)
       thisOctal%cornerVelocity(4) = velocityFunc(vector(x1,0.d0,z2),grid)
       thisOctal%cornerVelocity(5) = velocityFunc(vector(x2,0.d0,z2),grid)
       thisOctal%cornerVelocity(6) = velocityFunc(vector(x3,0.d0,z2),grid)
       thisOctal%cornerVelocity(7) = velocityFunc(vector(x1,0.d0,z3),grid)
       thisOctal%cornerVelocity(8) = velocityFunc(vector(x2,0.d0,z3),grid)
       thisOctal%cornerVelocity(9) = velocityFunc(vector(x3,0.d0,z3),grid)
    endif
666 continue

!    if(isnan(thisOctal%cornerVelocity(1)%x)) then
!          write(*,*) "cornervel",thisOctal%cornerVelocity(1)
!          write(*,*) velocityFunc(vector(x1,0.d0,z1),grid)
!          write(*,*) x1,z1
!          write(*,*) x2,z2
!          write(*,*) x3,z3
!       enddo
!    endif
    
  END SUBROUTINE fillVelocityCorners

  

  TYPE(vector) FUNCTION TTauriVelocity(point,grid)
    ! calculates the velocity vector at a given point for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    use input_variables, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               TTauriMstar, TTauriDiskHeight
                               
    IMPLICIT NONE

    TYPE(vector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    TYPE(vector)      :: vP
    REAL(double)  :: modVp
    REAL(double)  :: phi
    REAL(double)  :: r, rM, theta, y

    starPosn = grid%starPos1
    pointVec = (point - starPosn) * 1.e10_oc
    r = modulus( pointVec ) 
    if (r /= 0.d0) then
       theta = ACOS( max(-1.d0,min(1.d0,pointVec%z / dble(r) )))
    else
       theta = pi/2.d0
    endif
    if (theta == 0.d0) theta=1.d-20
    rM  = r / SIN(theta)**2
    y = SIN(theta)**2 

    ! test if the point lies within the star
    IF ( r < TTauriRstar ) THEN
      TTauriVelocity = vector(1.e-15,1.e-15,1.e-15)
      
    ! test if the point lies too close to the disk
    ELSE IF ( ABS(pointVec%z) < TTauriDiskHeight) THEN
      TTauriVelocity = vector(1.e-14,1.e-14,1.e-14)
    ! test if the point lies outside the accretion stream
    ELSE IF ((rM > TTauriRinner) .AND. (rM < TTauriRouter )) THEN
  
      vP = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
                  0.0, &
                 (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))
      modVp = SQRT((2.0 * bigG * TTauriMstar / TTauriRstar) * &
                     (TTauriRstar/r - TTauriRstar/rM))
      vP = (-1.0 * (modVp/cSpeed)) * vP
      phi = ATAN2(pointVec%y,pointVec%x)
      vP = rotateZ(vP,-phi)
      
      IF (pointVec%z < 0.0) vP%z = -vP%z
      TTauriVelocity = vP

    ELSE
      TTauriVelocity = vector(1.e-12,1.e-12,1.e-12)
    END IF

  END FUNCTION TTauriVelocity

  SUBROUTINE calcTTauriMassVelocity(thisOctal,subcell,grid) 
    ! calculates some of the variables at a given point for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    USE input_variables, ONLY : useHartmannTemp, maxHartTemp, TTauriRstar,&
                                TTauriRinner, TTauriRouter, ttau_acc_on

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(INOUT) :: grid

    TYPE(vector) :: point

    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec

    REAL :: r, rM, theta, bigR, thetaStar
    REAL :: bigRstar, thetaStarHartmann, bigRstarHartmann, rMnorm
    REAL :: tmp

    starPosn = grid%starPos1

    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn) * 1.e10_oc

    if (.not. ttau_acc_on) then
       ! we don't include the magnetopsherical accretion
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%rho(subcell) = 1.e-19
       IF (useHartmannTemp) thisOctal%temperature(subcell) = 6500.0
    else
    
       IF (TTauriInFlow(point,grid)) THEN
!      thisOctal%rho(subcell) = TTauriDensity(point,grid)
          thisOctal%rho(subcell) = Density(point,grid)
          thisOctal%inFlow(subcell) = .TRUE.
          IF (useHartmannTemp) THEN
             ! need to calculate the flow point AS IF the magnetosphere
             !   was the same size as the one in the Hartmann et al paper
             r = modulus( pointVec ) 
             theta = acos( pointVec%z  / r )
        
             rM  = r / SIN(theta)**2
        
             ! work out the intersection angle (at the stellar surface) for
             !   the current flow line  
             thetaStar = ASIN(SQRT(TTauriRstar/rM))
             bigRstar = TTauriRstar * SIN(thetaStar) 

             ! normalize the magnetic radius (0=inside, 1=outside)
             rMnorm = (rM-TTauriRinner) / (TTauriRouter-TTauriRinner)
        
             ! convert rM to Hartmann units
             rM = 2.2*(2.0*rSol) + rMnorm*(0.8*(2.0*rSol))
        
             bigR = SQRT(pointVec%x**2+pointVec%y**2)
             bigR = (bigR-bigRstar) / (TTauriRouter-bigRstar) ! so that we can rescale it
             bigR = min(bigR,TTauriRouter)
             bigR = max(bigR,0.0)

             ! get the equivalent bigR for the Hartmann geometry
             ! first work out the intersection angle (at the stellar surface) for
             !   the current flow line  
             thetaStarHartmann = ASIN(SQRT((2.0*rSol)/rM))
             ! work out bigR for this point
             bigRstarHartmann = (2.0*rSol) * SIN(thetaStarHartmann) 
             
             bigR = bigRstarHartmann + (bigR * (3.0*(2.0*rSol) - bigRstarHartmann))
             
             r = (rM * bigR**2)**(1./3.) ! get r in the Hartmann setup
             
             ! get theta in the Hartmann setup
             tmp = r/rM
             ! quick fix for out of range argments (R. Kurosawa)
             if (tmp > 1.0) tmp=1.0
             if (tmp < -1.0) tmp = -1.0
             theta = ASIN(SQRT(tmp))
             theta = MIN(theta,REAL(pi/2.0-0.0001))
             theta = MAX(theta,0.0001)
        
             rMnorm = MAX(rMnorm, 0.0001)
             rMnorm = MIN(rMnorm, 0.9999)
             !        thisOctal%temperature(subcell) = hartmannTemp(rMnorm, theta, maxHartTemp)
             thisOctal%temperature(subcell) = hartmannTemp2(rMnorm, theta, maxHartTemp)
          END IF
       ELSE  ! not inflow
          thisOctal%inFlow(subcell) = .FALSE.
          !     thisOctal%rho(subcell) = 1.e-25
          thisOctal%rho(subcell) = 1.e-19
          IF (useHartmannTemp) thisOctal%temperature(subcell) = 6500.0
       END IF
    end if ! ttau_acc_on

  thisOctal%velocity(subcell) = TTauriVelocity(point,grid)
  !thisOctal%velocity(subcell) = TTauriRotation(point,grid)    
     
  if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = 50.d5/cspeed!!!!!!!!!!!!!!!!!!!!
  IF ((thisoctal%threed).and.(subcell == 8)) &
       CALL fillVelocityCorners(thisOctal,grid,TTauriVelocity, .true.)
  
  IF ((thisoctal%twod).and.(subcell == 4)) &
       CALL fillVelocityCorners(thisOctal,grid,TTauriVelocity, .false.)

  

  END SUBROUTINE calcTTauriMassVelocity


  !
  !  Given a octal object, this routine will "reassign" the values of density at the outer
  !  edge of the T Taur accretion stream. THIS ROUTINE SHOULD BE CALLED AS A POST PROCESS (AFTER THE 
  !  NORMAL DENSITY OF ACCRETION FLOW HAS BEEN ASSIGNED). It has to be done this way since 
  !  if we were to assined very small density, the cells at the edge will not split (because they
  !  are split by the total mass, c.f. decideSplit routine in this module).
  !  
  RECURSIVE SUBROUTINE TTauri_fuzzy_edge(thisOctal)
    use input_variables, only: TTauriRinner, TTauriRouter

    IMPLICIT NONE
    TYPE(OCTAL), POINTER :: thisOctal    
    !
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    TYPE(vector)     :: rvec
    integer :: n
    real(double) :: r, theta, rM, w, p, rho
    real(double) :: rM_fuzzy_in, rM_fuzzy_out  ! beginning of the fuzzy edges
    !
    real(double), parameter :: scale = 7.0d0  ! a scale in exponential decay of density.
!    real(double), parameter :: scale = 10.0d0  ! a scale in exponential decay of density.
    !
    
    if (thisOctal%threeD) then
       n = 8
    else
       n = 4
    endif

    do subcell = 1, n    
       if (thisOctal%hasChild(subcell)) then
          ! find the index of the new child and decend the tree
          do i = 1, n
             if (thisOctal%indexChild(i) == subcell) then
                childPointer => thisOctal%child(i)
                CALL TTauri_fuzzy_edge(childPointer)
                exit
             end if
          end do
       else 
          ! This is a leaf, so modefiy the density value if required.
          if (thisOctal%inFlow(subcell)) then
             ! get the size and the position of the centre of the current cell
             rVec = subcellCentre(thisOctal,subcell)
             r = modulus(rvec)
             
             if (r/=0.0d0) then
                theta = ACOS( MIN(ABS(rVec%z/r),0.998_oc) )
             else
                theta=0.01
             end if
             
             if (theta == 0.0d0) theta = 0.01
             rM  = r / SIN(theta)**2
                          
             ! The fuzzy density starts from a 5-th of the thickness (2h) 
             ! below the surface.
             w = get_fuzzy_width() ! [10^10cm] using a fucntion in this module

             rM_fuzzy_in  = TTauriRinner*1.0d-10 + w   ! [10^10cm]
             rM_fuzzy_out = TTauriRouter*1.0d-10 - w   ! [10^10cm]
             
             ! If the point is close to the edge, we make it fuzzy
             if ( rM < rM_fuzzy_in) then
!                thisOctal%rho(subcell) =  1.0d-13  ! limit the minimum value...
                p = (rM_fuzzy_in-rM)
                rho = thisOctal%rho(subcell)* EXP(-scale*p/w)
                thisOctal%rho(subcell) =  MAX(rho, 1d-25)  ! limit the minimum value...
             elseif (rM > rM_fuzzy_out) then
                p = (rM - rM_fuzzy_out)
                rho = thisOctal%rho(subcell)* EXP(-scale*p/w)
                thisOctal%rho(subcell) =  MAX(rho, 1d-25)  ! limit the minimum value...
             else
                ! don't touch the density, and just continue
                continue
             end if

          end if  ! if inflow

       end IF

    end do

  END SUBROUTINE TTauri_fuzzy_edge


  !
  ! The function to check if a given point (rvec) is in "fuzzy" edge region
  logical function  in_fuzzy_edge(rvec) 
    use input_variables, only: TTauriRinner, TTauriRouter

    IMPLICIT NONE    
    TYPE(vector), intent(in):: rvec
    real(double) :: r, theta, rM,  w
    real(double) :: rM_fuzzy_in, rM_fuzzy_out  ! beginning of the fuzzy edges
    !
    !
    r = modulus(rvec)
    
    if (r/=0.0d0) then
       theta = ACOS( MIN(ABS(rVec%z/r),0.998_oc) )
    else
       theta=0.01
    end if
    
    if (theta == 0.0d0) theta = 0.01
    rM  = r / SIN(theta)**2
    
    
    ! The fuzzy density starts from a 5-th of the thickness (2h) 
    ! below the surface.
    w = get_fuzzy_width() ! [10^10cm] using a fucntion in this module

    rM_fuzzy_in  = TTauriRinner*1.0d-10 + w   ! [10^10cm]
    rM_fuzzy_out = TTauriRouter*1.0d-10 - w   ! [10^10cm]
             
    ! If the point is close to the edge, we make it fuzzy
    if ( rM < rM_fuzzy_in) then
       in_fuzzy_edge = .true.
    elseif (rM > rM_fuzzy_out) then
       in_fuzzy_edge = .true.
    else
       in_fuzzy_edge = .false.
    end if

  end function in_fuzzy_edge



  !
  ! The function to check if a given point (rvec) is in "fuzzy" edge region
  function  get_fuzzy_width()  RESULT(w)
    use input_variables, only: TTauriRinner, TTauriRouter

    IMPLICIT NONE    
    real(double) :: w  ! width in  [10^10cm]
    real(double) :: rM_center, h
    !
    ! bias more toward the edges of the accreation stream
    h = 0.5d0*(TTauriRouter - TTauriRinner)*1.0d-10  ! [10^10cm] distance from the centere to edge
    rM_center = 0.5d0*(TTauriRouter + TTauriRinner)*1.0d-10 ! [10^10cm]   
    
    ! The fuzzy density starts from a 5-th of the thickness (2h) 
    ! below the surface.
             w = 0.1d0*h;  ! Halpha081 model
!    w = 0.2d0*h;  ! a fifth for now  ! Halpha067 model
!             w = 0.25d0*h;  ! Halpha083 model
!             w = 0.3d0*h;  ! Halpha080 model


  end function get_fuzzy_width



  ! find the total mass in accretion flow
  ! Before the initial call the total mass must be set to zero
  ! The output should be in grams.
  RECURSIVE SUBROUTINE TTauri_accretion_mass(thisOctal, grid, total_mass)

    IMPLICIT NONE
    TYPE(OCTAL), POINTER :: thisOctal    
    TYPE(gridtype), INTENT(INOUT)    :: grid
    real(double), intent(inout) :: total_mass
    !
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    TYPE(vector)     :: point
    integer :: n
    real(double) :: d, dV
    !
    
    if (thisOctal%threeD) then
       n = 8
    else
       n = 4
    endif

    do subcell = 1, n    
       if (thisOctal%hasChild(subcell)) then
          ! find the index of the new child and decend the tree
          do i = 1, n
             if (thisOctal%indexChild(i) == subcell) then
                childPointer => thisOctal%child(i)
                CALL TTauri_accretion_mass(childPointer, grid, total_mass)
                exit
             end if
          end do
       else 
          point = subcellCentre(thisOctal,subcell)
          ! This is a leaf, so modefiy the density value if required.
          if (thisOctal%inFlow(subcell) .and. TTauriInFlow(point,grid)) then
             d = thisOctal%subcellsize  ! [10^10cm]
             if (thisOctal%threed) then
                dV = d*d*d   ! [10^30cm]
             else
                dv = 2.0_db*pi*d*d*SQRT(point%x*point%x + point%y*point%y)  ! [10^30cm]
             endif
             total_mass = total_mass + thisOctal%rho(subcell) * dV*1.0d30  ! [g]

             total_mass = total_mass + thisOctal%rho(subcell) * dV*1.0d30  ! [g]
          end if  ! if inflow

       end IF

    end do

  END SUBROUTINE TTauri_accretion_mass


  ! Scale the density in the accretion flow by the scale factor passed on this routine.
  ! Before the initial call the total mass must be set to zero
  RECURSIVE SUBROUTINE TTauri_accretion_scale_density(thisOctal, grid, scale)

    IMPLICIT NONE
    TYPE(OCTAL), POINTER :: thisOctal    
    TYPE(gridtype), INTENT(INOUT)    :: grid
    real(double), intent(in) :: scale
    !
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    TYPE(vector)     :: point
    integer :: n
    real(double) :: rho

    !
    
    if (thisOctal%threeD) then
       n = 8
    else
       n = 4
    endif

    do subcell = 1, n    
       if (thisOctal%hasChild(subcell)) then
          ! find the index of the new child and decend the tree
          do i = 1, n
             if (thisOctal%indexChild(i) == subcell) then
                childPointer => thisOctal%child(i)
                CALL TTauri_accretion_scale_density(childPointer, grid, scale)
                exit
             end if
          end do
       else 
          point = subcellCentre(thisOctal,subcell)
          ! This is a leaf, so modefiy the density value if required.
          if (thisOctal%inFlow(subcell) .and. TTauriInFlow(point,grid)) then
             rho = thisOctal%rho(subcell) * scale
             thisOctal%rho(subcell) = rho
          end if  ! if inflow
       end IF

    end do

  END SUBROUTINE TTauri_accretion_scale_density




  subroutine infallEnhancmentAMR(grid, distortionVector, nVec, timeStep, doDistortion, &
                                 particleMass, alreadyDoneInfall)
    ! creates an increase in the accretion rate for a T Tauri geometry.
    ! currently, the setup must be simple (no disc inclination etc.)
                                 
    type(GRIDTYPE), intent(inout) :: grid
    integer, intent(in) :: nVec
    logical, intent(in) :: doDistortion
    type(VECTOR), intent(inout) :: distortionVector(nVec)
    real, intent(in) :: timeStep
    real, intent(in) :: particleMass
    logical, intent(inout) :: alreadyDoneInfall

    real(oct) :: dTime
    real, parameter :: etaFac = 9
    real, parameter :: chiFac = 1
    type(vector) :: distortionVec(nVec)
    type(vector) :: thisVec
    type(vector) :: thisVel
    type(vector) :: thisVelOc
    logical :: setAllChanged
    
    integer, parameter :: nTimes = 1000
    integer :: i, j
    type(vector) :: starPos

    ! if we are not undoing a previous infall, we will assume that we need to
    ! run stateq on ALL subcells
    setAllChanged = .not. alreadyDoneInfall
    
    do i = 1, nVec
       distortionVec(i) = distortionVector(i)
    end do

    if (alreadyDoneInfall) then 
      ! we first have to undo the previous infall phase's changes
      call infallEnhancePrivate(grid%octreeRoot,distortionVec,undoPrevious=.true.,&
                                setAllChanged=.false.)
    end if
    
    write(*,*) "Time stepping vectors..."
    dTime = timeStep/real(nTimes,kind=db)
    do j = 1, nVec
      do i = 1, nTimes
        thisVec = distortionVec(j) / 1.e10_oc
        call amrGridValues(grid%octreeRoot,thisVec,velocity=thisVel)
        thisVelOc = thisVel
        thisVelOc = cSpeed * thisVel
        distortionVec(j) = distortionVec(j) + (dTime * thisVelOc)
      enddo
    enddo
    write(*,*) "done."
    
    do i = 1, nVec
       distortionVector(i) = distortionVec(i)
    end do
    
    if (doDistortion) then
      write(*,*) "Distorting grid..."    
      starPos = grid%starPos1
      call infallEnhancePrivate(grid%octreeRoot,distortionVec, undoPrevious=.false.,&
                                setAllChanged=(.not.alreadyDoneInfall))
      alreadyDoneInfall = .true.
      write(*,*) "done."
    end if
    
  contains

    recursive subroutine infallEnhancePrivate(thisOctal, distortionVec, undoPrevious,&
                                              setAllChanged)

      type(octal), pointer :: thisOctal
      type(vector), intent(in) :: distortionVec(nVec)
      logical, intent(in) :: undoPrevious
      logical, intent(in) :: setAllChanged
      
      integer, parameter :: nPhi = 360
      real(oct) :: phi
      integer :: containedParticles
      real :: addedDensity
      real :: deltaRho
      integer :: iSubcell, iChild
      type(octal), pointer :: child
      type(vector) :: trueVector
      type(vector) :: thisVec, thisVec2
 
      do iSubcell = 1, thisOctal%maxChildren, 1
        if (thisOctal%hasChild(iSubcell)) then
          
          ! descend to a child
          
          do iChild = 1, thisOctal%nChildren, 1
            if (thisOctal%indexChild(iChild) == iSubcell) exit
          end do 

          child => thisOctal%child(iChild)
          call infallEnhancePrivate(child, distortionVec, undoPrevious,&
                                    setAllChanged)
          
        else

          if (setAllChanged) thisOctal%changed(:) = .true.
          
          ! if we are out of the accretion flow, we can ignore this cell
          if (thisOctal%rho(iSubcell) < 1.e-24) cycle
          
          ! check if any of the 'particles' lie in the current subcell

          containedParticles = 0 
          do i = 1, nVec
            do j = 1, nPhi

              phi = twoPi * real(j-1) / real(nPhi-1)
              thisVec = rotateZ(distortionVec(i), phi)
              thisVec = thisVec * 1.0e-10_oc
              thisVec = thisVec + starPos

              if (inOctal(thisOctal,thisVec)) then
                 if (thisOctal%twod) then
                    thisVec2 = projecttoXZ(thisVec)
                 else
                    thisVec2 = thisVec
                 endif
                if (whichSubcell(thisOctal,thisVec2) == iSubcell) &
                  containedParticles = containedParticles + 1
                  
              else
                trueVector = subcellCentre(thisOctal,iSubcell) - starPos
                
                if (trueVector%z < 0.0_oc) then 

                  ! mirror the particles about the disc, and try again 
                  thisVec%z = starPos%z - (thisVec%z - starPos%z)
                  
                  if (inOctal(thisOctal,thisVec)) then
                 if (thisOctal%twod) then
                    thisVec2 = projecttoXZ(thisVec)
                 else
                    thisVec2 = thisVec
                 endif
                    if (whichSubcell(thisOctal,thisVec2) == iSubcell) &
                      containedParticles = containedParticles + 1
                  end if
                end if

              end if

            end do
          end do

          if (containedParticles > 0) then 
            
            addedDensity = (real(containedParticles) * particleMass) / &
                            (real(thisOctal%subcellSize,kind=db)*1.e10_db)**3
            ! calculate the density change in the subcell
            deltaRho = addedDensity / thisOctal%rho(iSubcell)

            if (.not. undoPrevious) then
              thisOctal%rho(iSubcell) = thisOctal%rho(iSubcell) + addedDensity
              thisOctal%etaLine(iSubcell) = thisOctal%etaLine(iSubcell) * (1.0+deltaRho)**2
              thisOctal%chiLine(iSubcell) = thisOctal%chiLine(iSubcell) * (1.0+deltaRho)**2
            else  
              ! we are reversing the changes made on a previous run of this
              !   subroutine.
              thisOctal%rho(iSubcell) = thisOctal%rho(iSubcell) - addedDensity
              thisOctal%etaLine(iSubcell) = thisOctal%etaLine(iSubcell) / (1.0+deltaRho)**2
              thisOctal%chiLine(iSubcell) = thisOctal%chiLine(iSubcell) / (1.0+deltaRho)**2
            end if

            thisOctal%changed(iSubcell) = .true.
              
          endif
        end if
      end do
    
    end subroutine infallEnhancePrivate

  end subroutine infallEnhancmentAMR
  

  
  PURE SUBROUTINE calcTTauriTemperature(thisOctal,subcell) 
    ! calculates the temperature in an octal for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    USE input_variables, ONLY: isoTherm, isoThermTemp

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell

    REAL :: rho

   ! TTauri star and disk
   ! we may need to track the minimum and maximum densities of the
   !   accretion flow 
! RK COMMENTED OUT THIS. IBM XFL compiler does not allow intrinsic function 
! in declearations.
!   real, save      :: TTauriMinRho = huge(TTauriMinRho)
!   real, save      :: TTauriMaxRho = tiny(TTauriMaxRho)
    real, parameter :: TTauriMinRho = 1.0e25
    real, parameter :: TTauriMaxRho = 1.0e-25

    IF (isoTherm) THEN
      thisOctal%temperature(subcell) = isoThermTemp
      RETURN
    END IF
      
    IF ( thisOctal%inFlow(subcell) ) THEN
      rho = thisOctal%rho(subcell)
      thisOctal%temperature(subcell) = MAX(5000.0, &
        7000.0 - ((2500.0 * rho/real(mHydrogen) - TTauriMinRho) / (TTauriMaxRho-TTauriMinRho)))
      ! we will initialise the bias distribution
      thisOctal%biasLine3D(subcell) = 1.0d0
      thisOctal%biasCont3D(subcell) = 1.0d0
    ELSE
      thisOctal%temperature(subcell) = 6500.0
      thisOctal%biasLine3D(subcell) = 1.0d-150
      thisOctal%biasCont3D(subcell) = 1.0d-150
    END IF


  
  END SUBROUTINE calcTTauriTemperature

  
  FUNCTION hartmannTemp(rM,inTheta,maxHartTemp)
    ! the paper: Hartman, Hewett & Calvet 1994ApJ...426..669H 
    !   quotes data from points along the magnetic field lines marked 
    !   on one of the figures. this function returns the temperature
    !   at any point within the accretion flow, by interpolating
    !   between the values in the published figure 6.

    USE unix_mod, only: unixgetenv

    REAL             :: hartmannTemp
    REAL, INTENT(IN) :: rM ! normalized (0=inside, 1=outside)
    REAL, INTENT(IN) :: inTheta
    REAL, INTENT(IN) :: maxHartTemp
    
    
    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: temperatures
    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: angles
    INTEGER, DIMENSION(5), SAVE             :: nSamples
    
    LOGICAL, SAVE:: alreadyLoaded = .FALSE.
    REAL :: radius
    REAL :: theta
    CHARACTER(LEN=80) :: dataDirectory 
    INTEGER :: errNo, i
    
    INTEGER :: subPos
    INTEGER :: fieldline
    INTEGER :: lowerLine, upperLine
    INTEGER :: iSample
    INTEGER :: maxSamples
    INTEGER :: lowerSample, upperSample
    REAL :: lowerTemp, upperTemp
    INTEGER :: iStat
    logical, save :: warned_already_01 = .false.
    logical, save :: warned_already_02 = .false.

    IF (.NOT. alreadyLoaded) THEN
      ! if this is the first time the function has been called, need
      !   to load in the data.
      
      call unixGetenv("TORUS_DATA",dataDirectory,i,errno)
      dataDirectory = trim(dataDirectory)//"/hartmann/"
      OPEN(31,FILE=TRIM(dataDirectory)//"tempProfile1.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
      OPEN(32,FILE=TRIM(dataDirectory)//"tempProfile2.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
      OPEN(33,FILE=TRIM(dataDirectory)//"tempProfile3.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF  
      OPEN(34,FILE=TRIM(dataDirectory)//"tempProfile4.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF  
      OPEN(35,FILE=TRIM(dataDirectory)//"tempProfile5.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 

      DO fieldline = 1, 5, 1
        READ(UNIT=(30+fieldline),FMT=*) nSamples(fieldLine)
      END DO

      maxSamples = MAXVAL(nSamples)

      ALLOCATE(temperatures(5,maxSamples))
      ALLOCATE(angles(5,maxSamples))
      ! make sure that there are no uninitialized variables
      temperatures = 0.0
      angles = 0.0
      
      DO fieldline = 1, 5, 1
        DO iSample = 1,  nSamples(fieldLine), 1
          READ(UNIT=(30+fieldline),FMT=*) angles(fieldLine,iSample),     &
                                          temperatures(fieldLine,iSample)
        END DO
        CLOSE(UNIT=(30+fieldline))
      END DO


      ! might need to rescale the temperature distribution to have a new maximum
      temperatures = temperatures * maxHartTemp/MAXVAL(temperatures)
      
      alreadyLoaded = .TRUE.
    END IF

    ! find which fieldlines bracket the point
   
    radius = rM * 0.8 ! because the Hartmann magnetic radii span 0.8 R_star
    radius = radius / 0.2 ! divide by the field line spacing

    IF (radius < 0.0) THEN 
       if (.not. warned_already_01) then
          PRINT *, 'In hartmannTemp, rM value is out of range! (',rM,')'
          PRINT *, '  assuming fieldlines 1-2'
          PRINT *, '  ==> Further warning surpressed.'
          warned_already_01 = .true.
       end if
       lowerLine = 1
       upperLine = 2 
       subPos    = 0.0
    ELSE IF (radius > 5.0) THEN
       if (.not. warned_already_02) then
          PRINT *, 'In hartmannTemp, rM value is out of range! (',rM,')'
          PRINT *, '  assuming fieldlines 1-2'
          PRINT *, '  ==> Further warning surpressed.'
          warned_already_02 = .true.
       end if
       lowerLine = 4
       upperLine = 5 
    ELSE 
       lowerLine = INT(radius) + 1
       upperLine = lowerLine + 1 
       subPos = radius - REAL(lowerline)
    END IF
      
    IF (inTheta > pi) THEN
      theta = inTheta - pi
    ELSE
      theta = inTheta
    END IF
    
    IF (inTheta > pi/2.0) theta = pi - theta
    
    ! find the two temperatures at the fieldlines 
    CALL locate(angles(lowerLine,1:nSamples(lowerLine)), &
                nSamples(lowerLine), ABS(theta), lowerSample)
    CALL locate(angles(upperLine,1:nSamples(upperLine)), &
                nSamples(upperLine), ABS(theta), upperSample)


    ! Forcing the values to be within the range...
    ! ... quick fix to avoid out of range problem... (R. Kurosawa)
    if (lowerSample <= 0) lowerSample =1
    if (upperSample <= 0) upperSample =1
    if (lowerSample >= nSamples(lowerLine)) lowerSample = nSamples(lowerLine)-1
    if (upperSample >= nSamples(upperLine)) upperSample = nSamples(upperLine)-1

                
    IF ((lowerSample == 0) .OR. (lowerSample >= nSamples(lowerLine)) .OR. &
       (upperSample == 0) .OR. (upperSample >= nSamples(upperLine))) THEN       
       PRINT *, 'In hartmannTemp, theta value is out of range! (',ABS(theta),')'
       STOP
    END IF
                
    ! interpolate along each fieldline
    lowerTemp = (temperatures(lowerLine,lowerSample) +        &
                 temperatures(lowerLine,lowerSample+1)) / 2.0
    upperTemp = (temperatures(upperLine,upperSample) +        &
                 temperatures(upperLine,upperSample+1)) / 2.0
   
    ! interpolate between the fieldLines
    
    hartmannTemp = (lowerTemp * (1.0-subPos)) + (upperTemp * subPos)
    
  END FUNCTION hartmannTemp


  !
  ! This is based on hartmannTemp. The temperature is interpolarted
  ! by using all five field lines instead of just two to avoid 
  ! the jumps in the temperature across the field lines. (RK)
  FUNCTION hartmannTemp2(rM,inTheta,maxHartTemp)
    ! the paper: Hartman, Hewett & Calvet 1994ApJ...426..669H 
    !   quotes data from points along the magnetic field lines marked 
    !   on one of the figures. this function returns the temperature
    !   at any point within the accretion flow, by interpolating
    !   between the values in the published figure 6.

    USE unix_mod, only: unixGetEnv

    REAL             :: hartmannTemp2
    REAL, INTENT(IN) :: rM ! normalized (0=inside, 1=outside)
    REAL, INTENT(IN) :: inTheta
    REAL, INTENT(IN) :: maxHartTemp
    
    
    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: temperatures
    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: angles
    INTEGER, DIMENSION(5), SAVE             :: nSamples
    
    LOGICAL, SAVE:: alreadyLoaded = .FALSE.
    REAL :: radius
    REAL :: theta
    CHARACTER(LEN=80) :: dataDirectory 
    INTEGER :: errNo, i
    
    INTEGER :: subPos
    INTEGER :: fieldline
    INTEGER :: lowerLine, upperLine
    INTEGER :: iSample
    INTEGER :: maxSamples
    INTEGER :: iStat
    logical, save :: warned_already_01 = .false.
    logical, save :: warned_already_02 = .false.
    real:: T(5), R(5), dT, R_ip
    integer :: indx_T(5)

    dt = 0.; hartmanntemp2 = 0.

    IF (.NOT. alreadyLoaded) THEN
      ! if this is the first time the function has been called, need
      !   to load in the data.
      
      call unixGetenv("TORUS_DATA",dataDirectory,i,errno)
      dataDirectory = trim(dataDirectory)//"/hartmann/"
      OPEN(31,FILE=TRIM(dataDirectory)//"tempProfile1.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
      OPEN(32,FILE=TRIM(dataDirectory)//"tempProfile2.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
      OPEN(33,FILE=TRIM(dataDirectory)//"tempProfile3.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF  
      OPEN(34,FILE=TRIM(dataDirectory)//"tempProfile4.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF  
      OPEN(35,FILE=TRIM(dataDirectory)//"tempProfile5.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 

      DO fieldline = 1, 5, 1
        READ(UNIT=(30+fieldline),FMT=*) nSamples(fieldLine)
      END DO

      maxSamples = MAXVAL(nSamples)

      ALLOCATE(temperatures(5,maxSamples))
      ALLOCATE(angles(5,maxSamples))
      ! make sure that there are no uninitialized variables
      temperatures = 0.0
      angles = 0.0
      
      DO fieldline = 1, 5, 1
        DO iSample = 1,  nSamples(fieldLine), 1
          READ(UNIT=(30+fieldline),FMT=*) angles(fieldLine,iSample),     &
                                          temperatures(fieldLine,iSample)
        END DO
        CLOSE(UNIT=(30+fieldline))
      END DO


      ! might need to rescale the temperature distribution to have a new maximum
      temperatures = temperatures * maxHartTemp/MAXVAL(temperatures)
      
      alreadyLoaded = .TRUE.
    END IF

    ! find which fieldlines bracket the point
   
    radius = rM * 0.8 ! because the Hartmann magnetic radii span 0.8 R_star
    radius = radius / 0.2 ! divide by the field line spacing

    IF (radius < 0.0) THEN 
       if (.not. warned_already_01) then
          PRINT *, 'In hartmannTemp, rM value is out of range! (',rM,')'
          PRINT *, '  assuming fieldlines 1-2'
          PRINT *, '  ==> Further warning surpressed.'
          warned_already_01 = .true.
       end if
       lowerLine = 1
       upperLine = 2 
       subPos    = 0.0
    ELSE IF (radius > 5.0) THEN
       if (.not. warned_already_02) then
          PRINT *, 'In hartmannTemp, rM value is out of range! (',rM,')'
          PRINT *, '  assuming fieldlines 1-2'
          PRINT *, '  ==> Further warning surpressed.'
          warned_already_02 = .true.
       end if
       lowerLine = 4
       upperLine = 5 
    ELSE 
       lowerLine = INT(radius) + 1
       upperLine = lowerLine + 1 
!       subPos = radius - REAL(lowerline)
       subPos = radius - 1.0
    END IF
      
    IF (inTheta > pi) THEN
      theta = inTheta - pi
    ELSE
      theta = inTheta
    END IF
    
    IF (inTheta > pi/2.0) theta = pi - theta
    
    ! find the index position in each field line 
    do i = 1, 5 
       CALL locate(angles(i,1:nSamples(i)), &
            nSamples(i), ABS(theta), indx_T(i))
       ! Forcing the values to be within the range...
       ! ... quick fix to avoid out of range problem... (R. Kurosawa)
       if (indx_T(i) <=0) indx_T(i)=1
       if (indx_T(i) >=nSamples(i)) indx_T(i)=nSamples(i)-1

       ! get the cooresponding temperature from each lines.       
       T(i) = temperatures(i, indx_T(i))
       R(i) = REAL(i)
    end do
                
    ! interpolate between the fieldLines 
    ! polynomial interpolation 
    ! using a routine in utils_mod.f90
    ! The output here is hartmannTemp2
    R_ip = radius
    if (R_ip<=1.0) R_ip=1.0
    if (R_ip>=5.0) R_ip=5.0
    call polint(R, T, 5, R_ip, hartmannTemp2, dT)

    
  END FUNCTION hartmannTemp2


  SUBROUTINE hartmannLines(fieldline,RAD,phi,grid,point,azVec,polVec,ok)
    ! the paper: Hartman, Hewett & Calvet 1994ApJ...426..669H 
    !   quotes data from points along the magnetic field lines marked 
    !   on one of the figures. this function returns the coordinates 
    !   of points along these lines, for obtaining comparable data. 

    use input_variables, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               TTauriDiskHeight
                               
    INTEGER, INTENT(IN)              :: fieldline
      ! 1 is the inner fieldline, 5 is the outer fieldline
    real(oct), INTENT(IN) :: RAD 
      ! RAD is projection of radius onto z=0 plane (in units of R_star)
    REAL(double), INTENT(IN)                 :: phi    ! azimuth angle (radians)
    TYPE(gridtype), INTENT(IN)       :: grid
    TYPE(vector), INTENT(OUT)   :: point  ! coordinates returned
    TYPE(vector), INTENT(OUT)        :: azVec  ! azimuth vector
    TYPE(vector), INTENT(OUT)        :: polVec ! poloidal vector
    LOGICAL, INTENT(OUT)             :: ok     ! coordinates are valid
    
    TYPE(vector) :: starPosn
    REAL              :: rM, theta, radius
    REAL              :: swap, y

    starPosn = grid%starPos1

    ok = .FALSE.
    
    ! test if the radius is sensible
    IF ( RAD  > (TTauriRouter/TTauriRstar) ) RETURN 

    rM = TTauriRinner + REAL(fieldline-1)*((TTauriRouter - TTauriRinner) / 4.0)
    rM = rM / TTauriRstar
    
    ! test if the radius is greater than this fieldline's max radius
    IF ( RAD  > rM ) RETURN 
    
    radius = (rM * RAD**2)**(1./3.)  

    ! test if the point lies inside the star
    IF ( radius < 1.0 ) RETURN
    
    y = radius / rM
    theta = ASIN(SQRT(y))

    ! set up vector in x-z plane
    point = vector(RAD, 0., radius*cos(theta))
    azVec = yHat

    ! test if the point lies too close to the disk
    IF ( ABS(point%z) < (TTauriDiskHeight/TTauriRstar)) RETURN

    ! calculate the poloidal vector
    polVec = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
                     0.0, &
                    (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))

    polVec = rotateZ(polVec,-phi)

    CALL RANDOM_NUMBER(swap)
    IF ( swap > 0.5 ) THEN 
      point%z  = -1.0_oc * point%z
      polVec%z = -1.0_oc * polVec%z
    END IF

    ! rotate about z-axis
    point  = rotateZ(point,REAL(-phi,kind=oct))
    polVec = rotateZ(polVec,-phi)
    azVec  = rotateZ(azVec,-phi)
    
    CALL normalize(polVec)

    ! convert to grid units (1.e10 cm)
    point = (TTauriRstar * 1.e-10_oc) * point
    
    ! offset to match centre of star
    point = point + starPosn
    
    ok = .TRUE.

  END SUBROUTINE hartmannLines

  SUBROUTINE writeHartmannValues(grid,variable)
    ! writes out values from the grid at places defined by the
    !   hartmannLines function
    
    use input_variables, only: TTauriRouter, TTauriRstar, &
                               MdotType, curtainsPhi1s, curtainsPhi1e, &
                               curtainsPhi2s, curtainsPhi2e
    
                               
    TYPE(gridtype), INTENT(IN)   :: grid
    CHARACTER(LEN=*), INTENT(IN) :: variable

    INTEGER, PARAMETER :: nBins = 100
    INTEGER, PARAMETER :: radiiPerBin = 1
    INTEGER, PARAMETER :: samplesPerRadii = 10
    INTEGER, PARAMETER :: samplesPerBin = radiiPerBin  * samplesPerRadii
    INTEGER, PARAMETER :: nFieldLines = 5
    
    real(double), DIMENSION(nFieldLines,nBins) :: outputArray
    real(double), DIMENSION(nFieldLines,nBins,grid%maxLevels) :: output2dArray
    REAL, DIMENSION(nBins)             :: radii
    TYPE(vector)                       :: azVec  ! azimuth vector
    TYPE(vector)                       :: polVec ! poloidal vector
    
    INTEGER           :: iBin, iRadius, iSample, iFieldLine
    REAL              :: rand
    REAL(double)              :: phi
    real(oct) :: radius
    LOGICAL           :: ok
    TYPE(vector) :: point 
    REAL              :: value
    real(double) :: valueDouble
    TYPE(vector)      :: velocityValue
    REAL(double)      :: etaLine
    REAL              :: chiLine
    REAL,DIMENSION(grid%maxLevels+1) :: departCoeff
    real(double),DIMENSION(grid%maxLevels) :: valueArrayDouble
    LOGICAL           :: scalar, logarithmic
    CHARACTER(80)     :: tempChar
    INTEGER           :: iLevel
    LOGICAL, SAVE     :: firstRun = .TRUE.
    REAL              :: curtainsFill
    
    scalar = .TRUE.
    logarithmic = .FALSE.
    curtainsFill = 1.
    IF (variable == 'hartmann_departCoeff' .or. variable == 'hartmann_N') scalar = .false.
   
    IF (scalar) THEN
      OPEN(UNIT=22, FILE=TRIM(variable)//'.csv', STATUS='REPLACE',FORM='FORMATTED')
    ELSE
      DO iLevel = 1, grid%maxLevels, 1
        write(tempChar,'(i2.2)') iLevel
        OPEN(UNIT=30+iLevel, FILE=TRIM(variable)//TRIM(tempChar)//'.csv', STATUS='REPLACE',FORM='FORMATTED')
      END DO
    END IF
    
    outputArray = 0.0
    output2dArray = 0.0
    
    DO iBin = 1, nBins, 1
    
      DO iFieldLine = 1, nFieldLines, 1
      
        DO iRadius = 1, radiiPerBin

          radius = 0.5*TTauriRstar + ((REAL(iBin-1)*REAL(radiiPerBin)+REAL(iRadius)) * &
                     ((TTauriRouter - (0.1*TTauriRstar)) / REAL(nBins*radiiPerBin)))
          radius = radius / TTauriRstar
          radii(iBin) = radius
        
          DO iSample = 1, samplesPerRadii
            IF (grid%amr2dOnly) THEN
              phi = pi/4.0 
            ELSE
              IF (TRIM(mDotType) == "constantcurtains") THEN
                IF (firstRun) THEN
                  PRINT *, 'Using only region in curtains for ''writeHartmannValues'''
                  curtainsFill = curtainsPhi1e-curtainsPhi1s 
                  IF (curtainsPhi2e-curtainsPhi2s /= curtainsFill) &
                    PRINT *, 'WARNING: curtains not equal size. Hartman values will be wrong'
                  curtainsFill = curtainsFill / 180.
                  PRINT *, '  curtains filling-factor: ',curtainsFill
                    
                  firstRun = .FALSE.
                END IF
                ok = .FALSE.
                DO 
                  CALL RANDOM_NUMBER(rand)
                  phi = (rand-0.5) * twoPi
                  IF (((phi > curtainsPhi1s+1.0).AND.(phi < curtainsPhi1e-1.0)) .OR.  &
                      ((phi > curtainsPhi2s+1.0).AND.(phi < curtainsPhi2e-1.0))) ok = .TRUE.
                  IF (ok) EXIT
                END DO
              ELSE
                CALL RANDOM_NUMBER(rand)
                phi = (rand-0.5) * twoPi
              END IF
            END IF

            CALL hartmannLines(iFieldline,radius,phi,grid,point,azVec,polVec,ok)
            IF (.NOT. ok) EXIT

            SELECT CASE (variable)

            CASE ('hartmann_logNH')
              CALL amrGridValues(grid%octreeRoot,point,rho=valueDouble) 
              valueDouble = (valueDouble/(mHydrogen*curtainsFill))
              logarithmic = .TRUE.
            
            CASE ('hartmann_temperature')
              CALL amrGridValues(grid%octreeRoot,point,temperature=value) 
              valueDouble = value / curtainsFill
              
            CASE ('hartmann_velPol')
              CALL amrGridValues(grid%octreeRoot,point,velocity=velocityValue) 
              value = velocityValue .dot. polVec
              value = (value * cSpeed) * 1.e-5 ! km/s
              valueDouble = abs(value) / curtainsFill
              
            CASE ('hartmann_velAz')
              CALL amrGridValues(grid%octreeRoot,point,velocity=velocityValue) 
              value = velocityValue .dot. azVec
              value = (value * cSpeed) * 1.e-5 ! km/s
              valueDouble = abs(value) / curtainsFill

            CASE ('hartmann_line')
              CALL amrGridValues(grid%octreeRoot,point,chiLine=chiLine,etaLine=etaLine)
              if (chiLine*curtainsFill /= 0.) then
                 valueDouble = etaLine / (chiLine*curtainsFill*10.) 
              else
                 valueDouble = tiny(valueDouble)
              endif
              logarithmic = .TRUE.
              
            CASE ('hartmann_Nelectron')
              CALL amrGridValues(grid%octreeRoot,point,Ne=valueDouble)
              valueDouble = valueDouble/curtainsFill
              logarithmic = .TRUE.
              
            CASE ('hartmann_Nlevel2')
              CALL amrGridValues(grid%octreeRoot,point,N=valueArrayDouble)
              valueDouble = valueArrayDouble(2)/curtainsFill
              logarithmic = .TRUE.
              
            CASE ('hartmann_logNe')
              CALL amrGridValues(grid%octreeRoot,point,Ne=valueDouble)
              valueDouble = valueDouble/curtainsFill
              logarithmic = .TRUE.
              
            CASE ('hartmann_NH')
              CALL amrGridValues(grid%octreeRoot,point,rho=valueDouble)
              valueDouble = valueDouble/(mHydrogen*curtainsFill)
              logarithmic = .TRUE.
              
            CASE ('hartmann_departCoeff')
              scalar = .false.
              CALL amrGridValues(grid%octreeRoot,point,departCoeff=departCoeff)
              valueArrayDouble = departCoeff(1:grid%maxLevels) / curtainsFill
              
            CASE ('hartmann_N')
              scalar = .false.
              CALL amrGridValues(grid%octreeRoot,point,N=valueArrayDouble)
              valueArrayDouble = valueArrayDouble / curtainsFill
              !logarithmic = .TRUE.
              
            CASE DEFAULT
              PRINT *, '''variable'' not recognized in writeHartmannValues'
              STOP
              
            END SELECT

            IF (scalar) THEN
              outputArray(iFieldLine,iBin) = outputArray(iFieldLine,iBin) + valueDouble 
            ELSE
              output2dArray(iFieldLine,iBin,:) = output2dArray(iFieldLine,iBin,:) + valueArrayDouble 
            END IF

          END DO ! iSample
          
          IF (.NOT. ok) EXIT
        
        END DO ! iRadius
        
        IF (ok) THEN
          IF (scalar) THEN
            outputArray(iFieldLine,iBin) = outputArray(iFieldLine,iBin) / REAL(samplesPerBin,KIND=db)
            IF (logarithmic) outputArray = LOG10(MAX(1.0e-10_db,outputArray))
          ELSE
            output2dArray(iFieldLine,iBin,:) = output2dArray(iFieldLine,iBin,:) / REAL(samplesPerBin,KIND=db)
            IF (logarithmic) output2dArray = LOG10(MAX(1.0e-10_db,output2dArray))
          END IF
        ELSE 
          IF (scalar) THEN
            outputArray(iFieldLine,iBin) = 0.0
            !outputArray(iFieldLine,iBin) = -1.0 * HUGE(outputArray(1,iBin))
          ELSE
            output2dArray(iFieldLine,iBin,:) = 0.0
            !output2dArray(iFieldLine,iBin,:) = -1.0 * HUGE(output2dArray(1,iBin,1))
          END IF
        END IF

      END DO ! iFieldLine
      
      IF (scalar) THEN
        WRITE(22,'(F11.4,3(" ",E11.4))') radii(iBin), outputArray(2:4,iBin)
        !                ^                                  ^^^   
        ! these constants might change!
      ELSE
        DO iLevel = 1, grid%maxLevels, 1
          WRITE(30+iLevel,'(F11.4,3(" ",E11.4))') radii(iBin), REAL(output2dArray(2:4,iBin,iLevel))
        END DO
      END IF
    
    END DO ! iBin
    
    IF (scalar) THEN
      CLOSE(UNIT=22)
    ELSE
      DO iLevel = 1, grid%maxLevels, 1
        CLOSE(UNIT=30+iLevel)
      END DO
    END IF

  END SUBROUTINE writeHartmannValues    
    
  subroutine initTTauriAMR(grid,theta1,theta2)

    use constants_mod
    use vector_mod

    use input_variables, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               dipoleOffset, contFluxFile, &
                               magStreamFile, thinDiskRin !, TTauriMstar
    USE magField, only: loadMagField

    implicit none
    
    type(GRIDTYPE), intent(inout)      :: grid                
    real, intent(out) :: theta1, theta2
   
    integer :: nNu 
    real :: nuArray(3000),fnu(3000)
    real :: tot
    integer :: i 
    real :: rStar



    rStar  = TTauriRstar / 1.e10
    
    grid%rCore = rStar
    grid%rStar1 = rStar
    grid%rStar2 = 0.
    grid%dipoleOffset = dipoleOffset
!    grid%diskRadius = TTauriRInner * 1.e-10
!    grid%diskRadius = TTauriDiskRin*rStar  ! [10^10cm]
    grid%diskRadius = ThinDiskRin*rStar  ! [10^10cm]
    grid%diskNormal = VECTOR(0.,0.,1.)
    grid%diskNormal = rotateX(grid%diskNormal,dble(grid%dipoleOffSet))
    grid%starPos1 = VECTOR(0.,0.,0.) ! in units of 1.e-10 cm
    grid%starPos2 = VECTOR(9.e9,9.e9,9.e9) ! in units of 1.e-10 cm


    theta1 = asin(sqrt(TTauriRstar/TTauriRouter))
    theta2 = asin(sqrt(TTauriRstar/TTauriRinner))

    !Laccretion = (REAL(bigG,KIND=db)* &
    !              REAL(TTauriMstar,KIND=db)* &
    !              !REAL(TTauriMdot,KIND=db)/ &
    !              REAL(10e-8*Msol*secstoyears,KIND=db)/ &
    !              REAL(TTauriRstar,kind=db))* &
    !        REAL((1.0_db-(2.0_db*TTauriRstar/(TTauriRouter+TTauriRinner))),KIND=db)

    !TaccretionDouble = Laccretion / REAL(((fourPi * TTauriRstar**2)*stefanBoltz* &
    !                                  abs(cos(theta1)-cos(theta2))),kind=db)

!    sAccretion = (fourPi * TTauriRstar**2)*abs(cos(theta1)-cos(theta2))!1.e20
!    tAccretionDouble = max(1.d0, tAccretionDouble)
!    Taccretion = TaccretionDouble**0.25

!    write(*,*) contfluxfile
    open(20,file=contFluxFile,status="old",form="formatted")
    nnu = 1
10  continue
    read(20,*,end=20) nuarray(nnu), fnu(nnu)
    nnu = nnu + 1
    goto 10
20  continue
    nnu = nnu  - 1
    close(20)
    tot = 0.
    do i = 1, nnu-1
       tot = tot + 0.5*(nuArray(i+1)-nuArray(i))*(fnu(i+1)+fnu(i))
    enddo
!    write(*,*) (fourPi*TTauriRstar**2)*tot/lSol," solar luminosities"   
!
!    write(*,*) (fourPi*TTauriRstar**2)*tot/ &
!               (fourPi*TTauriRstar**2*stefanBoltz*4000.**4)

    ! add the accretion luminosity spectrum to the stellar spectrum,
    ! write it out and pass it to the stateq routine.
    

    IF ( grid%geometry == "magstream" ) THEN           
      CALL loadMagField(fileName=magStreamFile,starPosn=grid%starPos1,Rstar=grid%rStar1)           
    END IF

    !open(22,file="star_plus_acc.dat",form="formatted",status="unknown")
    !do i = 1, nNu
    !   fNu(i) = fNu(i) + pi*blackbody(tAccretion, 1.e8*cSpeed/ nuArray(i))* &
    !            ((1.e20*sAccretion)/(fourPi*TTauriRstar**2))
    !   write(22,*) nuArray(i), fNu(i)
    !enddo
    !tot = 0.
    !do i = 1, nnu-1
    !   tot = tot + 0.5*(nuArray(i+1)-nuArray(i))*(fnu(i+1)+fnu(i))
    !enddo
    !write(*,*) "Final star+accretion luminosity (in solar)",(fourPi*TTauriRstar**2)*tot/lSol
    !close(22)
    !newContfile="star_plus_acc.dat"


  end subroutine initTTauriAMR


  subroutine calcTestDensity(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 30.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.

    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       thisOctal%rho(subcell) = rho * (grid%rInner / r)**2 
       thisOctal%temperature(subcell) = 30.
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaCont(subcell) = 0.
    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

! for test of diffusion zone

!    thisOctal%diffusionApprox(subcell) = .false.
!    if ((r > grid%rinner).and.(r < grid%rinner*2.)) then
!       thisOctal%temperature(subcell) = 1000.
!    endif
!    if ((r > grid%rOuter*0.5).and.(r < grid%rOuter)) then
!       thisOctal%temperature(subcell) = 100.
!    endif
!    if ((r > grid%rinner*2.).and.(r < 0.5*grid%rOuter)) then
!       thisOctal%diffusionApprox(subcell) = .true.
!    endif
    

  end subroutine calcTestDensity

  subroutine calcPathTestDensity(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = 1.d0

  end subroutine calcPathTestDensity

  subroutine calcLexington(thisOctal,subcell,grid)

    use input_variables, only : hydrodynamics
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(vector) :: rVec
    real(double) :: ethermal, gamma
    integer,save :: it
    real :: fac
    real,save :: radius(400),temp(400)
    integer :: i
    logical,save :: firsttime = .true.
    ethermal  = 0.;
    gamma = 5.d0/3.d0

!    if (firsttime) then
!       open(20,file="overview.txt", form="formatted",status="old")
!       it = 1
!77     continue
!       read(20,*,end=88) radius(it),temp(it)
!       radius(it) = (radius(it)+30.e17) / 1.e10
!       temp(it) = 10.**temp(it) 
!       it = it + 1
!       goto 77
!88     continue
!       it = it - 1
!       close(20)
!       firsttime = .false.
!    endif


    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.d-30
    thisOctal%temperature(subcell) = 10000.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.


    if (r > grid%rinner) then
       thisOctal%inFlow(subcell) = .true.
       thisOctal%rho(subcell) = 100.*mHydrogen
       thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
       thisOctal%ne(subcell) = thisOctal%nh(subcell)
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
       thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)

       thisOctal%ionFrac(subcell,1) = 1.e-10
       thisOctal%ionFrac(subcell,2) = 1.
       thisOctal%ionFrac(subcell,3) = 1.e-10
       thisOctal%ionFrac(subcell,4) = 1.       
       thisOctal%etaCont(subcell) = 0.
       thisOctal%temperature(subcell) = 10000.
!       if ((r > radius(1)).and.(r < radius(it))) then
!          call locate(radius, it, r, i)
!          fac = (r-radius(i))/(radius(i+1)-radius(i))
!          thisOctal%temperature(subcell) = temp(i) + fac * (temp(i+1)-temp(i))
!       endif

    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
!    thisOctal%dustTypeFraction(subcell,1)=0.d0
!    if (r < radius(it)*2.) then
!       thisOctal%dusttypeFraction(subcell, 1) = 0.d0
!    else
!       thisOctal%dusttypeFraction(subcell, 1) = 1.d0
!       thisOctal%rho(subcell) = 1.d3*mHydrogen
!       thisOctal%temperature(subcell) = 10.
!       thisOctal%nh(subcell) = thisOctal%rho(subcell)/mHydrogen
!       thisOctal%ne(subcell) = thisOctal%rho(subcell)/mHydrogen
!    endif


    if (hydrodynamics) then
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
       thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
       thisOctal%pressure_i(subcell) = (gamma-1.d0)* thisOctal%rho(subcell)*ethermal
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
       thisOctal%boundaryCondition(subcell) = 4
    endif

  end subroutine calcLexington

  subroutine calcStarburst(thisOctal,subcell,grid)

    use input_variables, only : hydrodynamics
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    real(double) :: ethermal, gamma = 5.d0/3.d0

    rVec = subcellCentre(thisOctal,subcell)

    thisOctal%rho(subcell) = tiny(thisOctal%rho(subcell))
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.
    thisOctal%rho(subcell) = 1.e2*mHydrogen
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-5
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    thisOctal%ionFrac(subcell,1) = 1.e-10
    thisOctal%ionFrac(subcell,2) = 1.
    thisOctal%ionFrac(subcell,3) = 1.e-10
    thisOctal%ionFrac(subcell,4) = 1.       
    thisOctal%etaCont(subcell) = 0.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
    thisOctal%dustTypeFraction(subcell,1)=1.d0

    if (hydrodynamics) then
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
       thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
       thisOctal%pressure_i(subcell) = (gamma-1.d0)* thisOctal%rho(subcell)*ethermal
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
       thisOctal%boundaryCondition(subcell) = 4
    endif
  end subroutine calcStarburst

  subroutine calcSymbiotic(thisOctal,subcell,grid)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    real(double) :: r, v, mdot


    mDot = 1.e-8*mSol/(365.25*24.*3600.)

    rVec = subcellCentre(thisOctal, subcell)

    r = modulus(rVec - VECTOR(-250.d0*rSol/1.d10,0.d0,0.d0))
    v = 10.e5 ! 10 km/s
    thisOctal%rho(subcell) = mDot/(fourPi * r**2 * 1.d20 * v)

    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .true.

    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
    thisOctal%dustTypeFraction(subcell,1)=0.d0

  end subroutine calcSymbiotic
  
  
  subroutine initWindTestAMR(grid)

    use constants_mod
    use vector_mod
    use input_variables

    implicit none
    
    type(GRIDTYPE), intent(inout)      :: grid                
   
    real :: rStar
    
    rStar  = (rSol * 20.) / 1.e10
    
    grid%rCore = rStar
    grid%rStar1 = rStar
    grid%rStar2 = 0.
    grid%rInner = rStar
    grid%rOuter = rStar * 500.
    grid%starPos1 = VECTOR(0.,0.,0.) ! in units of 1.e-10 cm
    grid%starPos2 = VECTOR(9.e9,9.e9,9.e9) ! in units of 1.e-10 cm

  end subroutine initWindTestAMR


  subroutine calcProtoDensity(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 1.e-3
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .false.

    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       thisOctal%rho(subcell) = grid%densityScaleFac * rho * (grid%rInner / r)**rPower 
       thisOctal%temperature(subcell) = 10.
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaCont(subcell) = 0.
    endif

! put different sort of dust within 4 rinner

    if (r < (grid%rInner*4.)) then
       thisOctal%dustType(subcell) = 1
    else
       thisOctal%dustType(subcell) = 2
    endif

    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine calcProtoDensity


  subroutine calcGammaVel(thisOctal, subcell, grid)
    use input_variables
    type(OCTAL) :: thisOctal
    type(GRIDTYPE) :: grid
    integer :: subcell
    integer, parameter :: nsteps = 1000
    real :: ydist(nSteps), xDist(nSteps)
    real :: d1, d2, curlyR, dx, dybydx, ddash1, ddash2
    type(VECTOR) :: direction(nSteps), stagVec, shockdirection
    type(VECTOR) :: rvec,direction2
    real :: stagpoint
    real :: v, r
    integer :: i, j

    thisOctal%inflow(subcell) = .true.
    thisOctal%temperature(subcell) = teff1
    massRatio = mass1/mass2

    d1 = binarySep * (1./(massRatio+1.))
    d2 = binarySep - d1

    momRatio = (mdot1 * vterm1) / (mdot2 * vterm2)

    curlyR = sqrt(momRatio)         ! = d1/d2 (Equ 1.) 
    ! Stevens, Blondin & Pollock 1992
    ! ApJ 386 265

    stagPoint = binarySep * sqrt(momRatio) / (1. + sqrt(momRatio))

    stagVec = VECTOR(0.,0.,-d1) + VECTOR(0., 0., stagPoint)

    shockDirection = VECTOR(0.,0.,0.)


    dx = (grid%octreeRoot%subcellSize-stagVec%z)/real(nsteps)
    yDist(1) = dx*10.
    xDist(1) = stagVec%z 
    direction(1) = VECTOR(0.,1.,0.)
    do j = 2,nsteps
       xDist(j) = stagVec%z + real(j-1)*(grid%octreeroot%subcellsize-stagVec%z)/real(nSteps-1)
       ddash1 = sqrt((xDist(j-1) + d1)**2 + yDist(j-1)**2)
       ddash2 = sqrt((xDist(j-1) - d2)**2 + yDist(j-1)**2)

       dybydx = ((curlyR*ddash2**2 + ddash1**2)*yDist(j-1)) / &
            (curlyR *ddash2**2*(xDist(j-1)+d1) + ddash1**2*(xDist(j-1)-d2))

       yDist(j) = yDist(j-1) + dx * dyBydx

       direction(j) = VECTOR(dx,dybydx,0.)
       call normalize(direction(j))
    enddo

    rVec = subcellCentre(thisOctal, subcell)
    

    if (rVec%z > xDist(1)) then
       call locate(xDist, nSteps, real(rVec%z), i)
       if (rvec%x > ydist(i)) then
          direction2 = (rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
          call normalize(direction2)
          r = modulus(rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
          if (r > rStar1) then
             v = vNought1 + (vterm1-vNought1)*(1.d0 - rstar1/r)**beta1
             thisOctal%rho(subcell) = mdot1 / (fourPi * (r*1.d10)**2 * v)
             thisOctal%temperature(subcell) = 0.9d0 * teff1
             thisOctal%velocity(subcell) = dble(v/cspeed) * direction2
!             thisOctal%atomAbundance(subcell, 1) =  1.d-5 / mHydrogen
             thisOctal%atomAbundance(subcell, :) =  1.d0 / (4.d0*mHydrogen)
!             thisOctal%atomAbundance(subcell, 2) =  1.d0 / (4.d0*mHydrogen)
          endif
       else
          direction2 = (rVec - VECTOR(0.d0, 0.d0, dble(d2)))
          call normalize(direction2)
          r = modulus(rVec - VECTOR(0.d0, 0.d0, dble(d2)))
          if (r > rStar2) then
             v = vNought2 + (vterm2-vNought2)*(1.d0 - rstar2/r)**beta2
             thisOctal%rho(subcell) = mdot2 / (fourPi * (r*1.d10)**2 * v)
             thisOctal%temperature(subcell) = 0.9d0 * teff2
             thisOctal%velocity(subcell) = dble(v/cspeed) * direction2
!             thisOctal%atomAbundance(subcell, 1) =  0.71d0 / mHydrogen
             thisOctal%atomAbundance(subcell, :) =  0.27d0 / (4.d0*mHydrogen)
!             thisOctal%atomAbundance(subcell, 2) =  0.27d0 / (4.d0*mHydrogen)
          endif
       endif
    else
       direction2 = (rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
       call normalize(direction2)
       r = modulus(rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
       if (r > rStar1) then
          v = vNought1 + (vterm1-vNought1)*(1.d0 - rstar1/r)**beta1
          thisOctal%rho(subcell) = mdot1 / (fourPi * (r*1.d10)**2 * v)
          thisOctal%temperature(subcell) = 0.9d0 * teff1
          thisOctal%velocity(subcell) = dble(v/cspeed) * direction2
!          thisOctal%atomAbundance(subcell, 1) =  1.d-5 / mHydrogen
          thisOctal%atomAbundance(subcell, :) =  1.d0 / (4.d0*mHydrogen)
!          thisOctal%atomAbundance(subcell, 2) =  1.d0 / (4.d0*mHydrogen)
       endif
    endif


    thisOctal%microturb(subcell) = sqrt((2.d0*kErg*thisOctal%temperature(subcell))/mhydrogen)/cspeed

  end subroutine calcGammaVel

    
  TYPE(vector)  function gammaVelVelocity(point, grid)
    use input_variables, only : vnought1, vnought2, vterm1, vterm2, beta1, beta2
    use input_variables, only : rstar1, rstar2, mass1, mass2, mdot1, mdot2, binarySep

    type(vector), intent(in) :: point
    type(vector) :: rvec
    type(GRIDTYPE), intent(in) :: grid
    real(double) :: v, r
    real :: massRatio
    integer, parameter :: nsteps = 1000
    real :: ydist(nSteps), xDist(nSteps)
    real :: d1, d2, curlyR, dx, dybydx, ddash1, ddash2
    real(double) :: momRatio
    type(VECTOR) :: direction(nSteps), stagVec, shockdirection
    type(VECTOR) :: direction2
    real :: stagpoint
    integer :: i, j


    rVec = point

    gammaVelVelocity = VECTOR(0.d0, 0.d0, 0.d0)
    r = modulus(rVec)

    massRatio = mass1/mass2

    d1 = binarySep * (1./(massRatio+1.))
    d2 = binarySep - d1

    momRatio = (mdot1 * vterm1) / (mdot2 * vterm2)

    curlyR = sqrt(momRatio)         ! = d1/d2 (Equ 1.) 
    ! Stevens, Blondin & Pollock 1992
    ! ApJ 386 265

    stagPoint = binarySep * sqrt(momRatio) / (1. + sqrt(momRatio))

    stagVec = VECTOR(0.,0.,-d1) + VECTOR(0., 0., stagPoint)

    shockDirection = VECTOR(0.,0.,0.)


    dx = (grid%octreeRoot%subcellSize-stagVec%z)/real(nsteps)
    yDist(1) = dx*10.
    xDist(1) = stagVec%z 
    direction(1) = VECTOR(0.,1.,0.)
    do j = 2,nsteps
       xDist(j) = stagVec%z + real(j-1)*(grid%octreeroot%subcellsize-stagVec%z)/real(nSteps-1)
       ddash1 = sqrt((xDist(j-1) + d1)**2 + yDist(j-1)**2)
       ddash2 = sqrt((xDist(j-1) - d2)**2 + yDist(j-1)**2)

       dybydx = ((curlyR*ddash2**2 + ddash1**2)*yDist(j-1)) / &
            (curlyR *ddash2**2*(xDist(j-1)+d1) + ddash1**2*(xDist(j-1)-d2))

       yDist(j) = yDist(j-1) + dx * dyBydx

       direction(j) = VECTOR(dx,dybydx,0.)
       call normalize(direction(j))
    enddo

   
    if (rVec%z > xDist(1)) then
       call locate(xDist, nSteps, real(rVec%z), i)
       if (rvec%x > ydist(i)) then
          direction2 = (rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
          call normalize(direction2)
          r = modulus(rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
          if (r > rStar1) then
             v = vNought1 + (vterm1-vNought1)*(1.d0 - rstar1/r)**beta1
             gammaVelVelocity = dble(v/cspeed) * direction2
          endif
       else
          direction2 = (rVec - VECTOR(0.d0, 0.d0, dble(d2)))
          call normalize(direction2)
          r = modulus(rVec - VECTOR(0.d0, 0.d0, dble(d2)))
          if (r > rStar2) then
             v = vNought2 + (vterm2-vNought2)*(1.d0 - rstar2/r)**beta2
             gammaVelVelocity = dble(v/cspeed) * direction2
          endif
       endif
    else
       direction2 = (rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
       call normalize(direction2)
       r = modulus(rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
       if (r > rStar1) then
          v = vNought1 + (vterm1-vNought1)*(1.d0 - rstar1/r)**beta1
          gammaVelVelocity = dble(v/cspeed) * direction2
       endif
    endif


  end function gammaVelVelocity




  subroutine calcWRShellDensity(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(vector) :: rVec
    real(double) :: v
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 0.9*teff
    thisOctal%etaCont(subcell) = 1.e-30
    thisOctal%inFlow(subcell) = .false.
    thisOctal%velocity(subcell) = VECTOR(0.,0.,0.)

    if (((r-thisOctal%subcellSize/2.d0) > grid%rInner)) then !.and.(r < grid%rOuter)) then
       thisOctal%rho(subcell) = density(rVec, grid)
       thisOctal%temperature(subcell) = 30000.d0 + (100000.d0-30000.d0)*(rcore/r)**4
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaCont(subcell) = 0.
       thisOctal%ne(subcell) = thisOctal%rho(subcell)/mHydrogen
       v = 0.001d5+(vterm-0.001d5)*(1.d0 - grid%rinner/r)**beta
       thisOctal%microturb(subcell) = 20.d5/cspeed
       thisOctal%velocity(subcell) = rVec
       thisOctal%inFlow(subcell) = .true.
       call normalize(thisOctal%velocity(subcell))
       thisOctal%velocity(subcell) = thisOctal%velocity(subcell) * v/cSpeed

    endif
    CALL fillVelocityCorners(thisOctal,grid,wrshellVelocity,thisOctal%threed)


    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

!    thisOctal%atomAbundance(subcell, 1) =  1.d0 / (mHydrogen)

    thisOctal%atomAbundance(subcell, 1) =  0.9d0 * 1.d0 / (4.d0*mHydrogen)
    if (natom>1) &
    thisOctal%atomAbundance(subcell, 2) =  0.9d0 * 1.d0 / (4.d0*mHydrogen)
    

  end subroutine calcWRShellDensity

  subroutine calcHydro1DDensity(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    type(VECTOR) :: rVec
    real(double) :: gd, xmid, x, z, r , zprime


    rVec = subcellCentre(thisOctal, subcell)
    xmid = (x1 + x2)/2.d0
    x = rVec%x
    z = rVec%z
    gd = 0.05d0 * (x2 - x1)
    if (thisOctal%twod) then
       r = modulus(rVec - VECTOR(0.5d0, 0.d0, 0.5d0))
    else
       r = modulus(rVec - VECTOR(0.5d0, 0.d0, 0.0d0))
    endif
!    thisOctal%rho(subcell) = 1.d0 + 0.3d0 * exp(-r**2/gd**2)
    thisOctal%energy(subcell) = 2.5d0
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
    thisOctal%pressure_i(subcell) = 1.d0

    thisOctal%rho(subcell) = 1.d0

!    if (thisOctal%twod) then
!       xMid = 0.d0 - z
!    endif

!    if (thisOctal%threed) then
    xMid = -0.75d0 - z
!    endif

    if  (thisOctal%threed) then
       zprime = -1.d0 - (rVec%x+rVec%y)
    else
!       zprime = -0.5d0 - rVec%x
       zprime = -rVec%x
!       zprime = sqrt(0.5d0**2 - rVec%x**2) - 0.5d0
    endif

!    if (rvec%z < zprime) then
    if (rvec%x < 0.5d0) then
       thisOctal%rho(subcell) = 1.d0
       thisOctal%energy(subcell) = 2.5d0
       thisOctal%pressure_i(subcell) = 1.d0
       thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    else
       thisOctal%rho(subcell) = 0.125d0
       thisOctal%energy(subcell) = 2.d0
       thisOctal%pressure_i(subcell) = 0.1d0
       thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    endif

    if (thisOctal%oneD) then
       if (rvec%x < 0.5d0) then
          thisOctal%rho(subcell) = 1.d0
          thisOctal%energy(subcell) = 2.5d0
          thisOctal%pressure_i(subcell) = 1.d0
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       else
          thisOctal%rho(subcell) = 0.125d0
          thisOctal%energy(subcell) = 2.d0
          thisOctal%pressure_i(subcell) = 0.1d0
          thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
       endif
    endif

    thisOctal%phi_i(subcell) = 0.d0
    thisOctal%boundaryCondition(subcell) = 1
    thisOctal%gamma(subcell) = 7.d0/5.d0
    thisOctal%iEquationOfState(subcell) = 0

  end subroutine calcHydro1DDensity

  subroutine calcKelvinDensity(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    type(VECTOR) :: rVec
    real :: u1 !, u2
    real(double), parameter :: gamma = 5.d0/3.d0


    rVec = subcellCentre(thisOctal, subcell)
    
    if (abs(rVec%z) > 0.25d0) then
       thisOctal%velocity(subcell) = VECTOR(-0.50, 0.d0, 0.d0)
       thisOctal%rho(subcell) = 1.d0
    else
       thisOctal%velocity(subcell) = VECTOR(0.50, 0.d0, 0.d0)
       thisOctal%rho(subcell) = 2.d0
    endif

!    call random_number(u1)
!    call random_number(u2)
!    u1 = (2.*u1-1.) * 0.005
!    u2 = (2.*u2-1.) * 0.005

    u1 = 0.d0
!    if (abs(rVec%z-0.25d0) < 0.025d0) then
!       u1 = 0.025d0 * sin ( -twoPi*(rvec%x+0.5d0)/(1.d0/6.d0) )
!    else if (abs(rVec%z+0.25d0) < 0.025d0) then
!       u1 = 0.025d0 * sin (  twoPi*(rvec%x+0.5d0)/(1.d0/6.d0) )
!    endif
    u1 = -0.025d0
    
    thisOctal%velocity(subcell) = VECTOR(0.5d0, 0.d0, 0.d0)/cspeed
    if (modulus(rvec) < 0.2d0) thisOctal%rho(subcell) = 0.5d0
!    if (inSubcell(thisOctal, subcell, VECTOR(0.295d0,0.d0, -0.005d0))) thisOctal%rho(subcell) = 0.5d0
    thisOctal%pressure_i(subcell) = 2.5d0
!    thisOctal%velocity(subcell) = (thisOctal%velocity(subcell) + VECTOR(0.d0, 0.d0, u1))/cSpeed

    thisOctal%energy(subcell) = thisOctal%pressure_i(subcell) /( (gamma-1.d0) * thisOctal%rho(subcell))
    thisOctal%boundaryCondition(subcell) = 2
    thisOctal%gamma(subcell) = 7.d0/5.d0
    thisOctal%iEquationOfState(subcell) = 1
  end subroutine calcKelvinDensity

  subroutine calcBonnorEbertDensity(thisOctal,subcell,grid)

    use input_variables, only : photoionization
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    type(VECTOR) :: rVec
    real(double), parameter :: gamma = 5.d0/3.d0
    real(double) :: eThermal, rMod, fac
    logical, save :: firstTime = .true.
    integer, parameter :: nr = 1000
    real(double), save :: r(nr), rho(nr)
    integer :: i

    if (firstTime) then
       firstTime = .false.
       call bonnorEbertRun(10.d0, 2.d0, 1000.d0*2.d0*mhydrogen, 6.d0,  nr, r, rho)
       r = r / 1.d10
       if (myrankGlobal==1) then
          do i =1 , nr
             write(55, *) r(i)*1.d10/autocm, rho(i)
          enddo
       endif
    endif

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec)
    if (rMod < r(nr)) then
       call locate(r, nr, rMod, i)
       fac = (rMod-r(i))/(r(i+1)-r(i))
       thisOctal%rho(subcell) = rho(i) + fac*(rho(i+1)-rho(i))
       thisOctal%temperature(subcell) = 10.d0
    else
       thisOctal%rho(subcell) = rho(nr)
       thisOctal%temperature(subcell) = 10.d0
    endif

!    thisOctal%rho(subcell) = rho(nr)

    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    ethermal = (1.d0/(2.d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = thisOctal%rho(subcell)*ethermal
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)
    thisOctal%boundaryCondition(subcell) = 4 ! outflow no inflow
    thisOctal%phi_i(subcell) = -bigG * 6.d0 * mSol / (modulus(rVec)*1.d10)
    thisOctal%gamma(subcell) = 1.d0
    thisOctal%iEquationOfState(subcell) = 1
      
    if (photoionization) then
       thisOctal%inFlow(subcell) = .true.
       thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
       thisOctal%ne(subcell) = thisOctal%nh(subcell)
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
       thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    
       thisOctal%ionFrac(subcell,1) = 1.e-10
       thisOctal%ionFrac(subcell,2) = 1.
       thisOctal%ionFrac(subcell,3) = 1.e-10
       thisOctal%ionFrac(subcell,4) = 1.       
       thisOctal%etaCont(subcell) = 0.
    endif
  end subroutine calcBonnorEbertDensity

  subroutine calcUniformSphere(thisOctal,subcell,grid)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    type(VECTOR) :: rVec
    real(double) :: eThermal, rMod

    rVec = subcellCentre(thisOctal, subcell)
    rMod = modulus(rVec)
    if (rMod < (pctocm/1.d10)) then
       thisOctal%rho(subcell) = 1000.d0*2.d0*mHydrogen
       thisOctal%temperature(subcell) = 10.d0
    else
       thisOctal%rho(subcell) = 1.d-24
       thisOctal%temperature(subcell) = 10.d0
    endif
    thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
    ethermal = (1.d0/(2.d0*mHydrogen))*kerg*thisOctal%temperature(subcell)
    thisOctal%pressure_i(subcell) = ((4.1317d9*1.d-20)/(1.d-20**2)) * thisOctal%rho(subcell)**2
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%rhoe(subcell) = thisOctal%rho(subcell) * thisOctal%energy(subcell)

    thisOctal%boundaryCondition(subcell) = 2

    thisOctal%phi_i(subcell) = -bigG * 100.d0 * mSol / (modulus(rVec)*1.d10)
    thisOctal%gamma(subcell) = 2.d0
    thisOctal%iEquationOfState(subcell) = 3
  end subroutine calcUniformSphere



  subroutine calcSedovDensity(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    type(VECTOR) :: rVec, cVec
    real(double) :: gamma, ethermal
    logical :: blast
    type(VECTOR) :: lVec, offset

    if (thisOctal%threed) then
       lVec = VECTOR(0., 0., 2.)
    else
       lVec = VECTOR(0., 2., 0.)
    endif

    offset = vector(0., 0., 0.)

    gamma = 7.d0/3.d0
    rVec = subcellCentre(thisOctal, subcell)
    cVec = VECTOR(0.d0, 0.0d0, -0.d0)
    blast = .false.
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)

    if (modulus(rvec-cVec) < 0.45d0) blast = .true.

    if (blast) then
       thisOctal%rho(subcell) = 0.1
       ethermal = 1.d-4
       thisOctal%velocity(subcell) = (1./cspeed) * ((rVec .cross. lVec)+offset)
       thisOctal%energy(subcell) = ethermal + 0.5 * (cspeed*modulus(thisOctal%velocity(subcell)))**2
       thisOctal%pressure_i(subcell) = (gamma-1.d0)*thisOctal%rho(subcell)*ethermal
    else
       thisOctal%rho(subcell) = 1.d-4
       ethermal = 0.1d0
       thisOctal%velocity(subcell) = (1./cspeed) * ((rVec .cross. lVec)+offset)
       thisOctal%energy(subcell) = ethermal + 0.5 * (cspeed*modulus(thisOctal%velocity(subcell)))**2
       thisOctal%pressure_i(subcell) = (gamma-1.d0)*thisOctal%rho(subcell)*ethermal
    endif
    thisOctal%energy(subcell) = 0.1d0

    thisOctal%boundaryCondition(subcell) = 4
    Thisoctal%phi_i(subcell) = 0.d0!-1.d-1 / modulus(rVec-cVec)
  end subroutine calcSedovDensity

  subroutine calcProtoBinDensity(thisOctal,subcell,grid)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    type(VECTOR) :: rVec
    real(double) :: gamma, ethermal
    real(double) :: rho0, r0, n, soundSpeed, omega, a, r, v, phi
    real(double) :: inertia, beta, rCloud, mCloud, eGrav

    gamma = 7.d0/4.d0
    rho0 = 1.0d0
    r0 = 0.4d0
    soundSpeed = 0.01d0
    omega = 0.8d0
    a = 0.1d0
    n = 2.d0
    ethermal = 0.1d0

    rVec = subcellCentre(thisOctal, subcell)
    r = modulus(rVec)
    thisOctal%rho(subcell) = rho0 / (1.d0 + (r/r0)**n)**(4.d0/n)
    if (thisOctal%twoD) then
       phi = atan2(rVec%z, rVec%x)
    else
       phi = atan2(rVec%y, rVec%x)
    endif

    soundSpeed = sqrt(gamma*(gamma-1.d0)*eThermal)

    thisOctal%rho(subcell) = thisOctal%rho(subcell) * (1.d0 + a*cos(2.d0*phi))
    if (thisOctal%twoD) then
       r = modulus(rVec)
       v = (4.d0*omega*r / (r0 + sqrt(r0**2 + r**2))) * soundSpeed
       thisOctal%velocity(subcell) = (1./cspeed)*VECTOR(v*cos(phi+pi/2.d0), 0., v*sin(phi+pi/2.d0))
    else
       r = sqrt(rVec%x**2 + rVec%y**2)
       v = (4.d0*omega*r / (r0 + sqrt(r0**2 + r**2))) * soundSpeed
       thisOctal%velocity(subcell) = (1./cspeed)*VECTOR(v*cos(phi),  v*sin(phi), 0.)
    endif

    thisOctal%pressure_i(subcell) = (gamma-1.d0)*thisOctal%rho(subcell)*ethermal
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%boundaryCondition(subcell) = 4


    mCloud = 1.d0 * msol
    rCloud = 7.d15
    inertia = (2.d0/5.d0)*mCloud*rCloud**2
    eGrav = 3.d0/5.d0 * bigG * mCloud**2 / rCloud
    beta = 0.3d0
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)

    omega = sqrt(2.d0 * beta * eGrav / inertia)
    thisOctal%phi_i(subcell) = -bigG * mCloud / (modulus(rVec)*1.d10)
    thisOctal%velocity(subcell) = vector(0., 0., 0.)

    if (thisOctal%twoD) then
       phi = atan2(rVec%z, rVec%x)
    else
       phi = atan2(rVec%y, rVec%x)
    endif

    if (thisOctal%twoD) then
       r = modulus(rVec)*1.d10
       v = r * omega
       thisOctal%velocity(subcell) = (1./cspeed)*VECTOR(v*cos(phi+pi/2.d0), 0., v*sin(phi+pi/2.d0))
    else
       r = sqrt(rVec%x**2 + rVec%y**2)*1.d10
       v = r * omega
       thisOctal%velocity(subcell) = (1./cspeed)*VECTOR(v*cos(phi+pi/2.d0),  v*sin(phi+pi/2.d0), 0.0)
    endif

    if (modulus(rVec) < (rCloud/1.d10)) then
       thisOctal%rho(subcell) = mCloud / (4.d0/3.d0*pi*rCloud**3)
       thisOctal%rho(subcell) = thisOctal%rho(subcell) * (1.d0 + 0.5d0 * cos(2.d0 * phi)) !m=2 dens perturbation
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
    else
       thisOctal%rho(subcell) = 1.d-2 * mCloud / (4.d0/3.d0*pi*rCloud**3)
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * 10.d0
    endif


    thisOctal%pressure_i(subcell) = (gamma-1.d0)* thisOctal%rho(subcell)*ethermal
    thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
    thisOctal%boundaryCondition(subcell) = 4
!    thisOctal%phi_i(subcell) = 0.d0

  end subroutine calcProtoBinDensity
    


  subroutine benchmarkDisk(thisOctal,subcell,grid)

    use input_variables, ONLY : rInner, rOuter, height, rho
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, hr, rd
    real(double), parameter :: min_rho = 1.0d-33 ! minimum density
    TYPE(vector) :: rVec

    real :: rInnerGap, rOuterGap
    logical :: gap

    gap = .false.

    rInnerGap = 2. * auToCm / 1.e10
    rOuterGap = 3. * auToCm / 1.e10
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = min_rho
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    rd = rOuter / 2.
    r = sqrt(rVec%x**2 + rVec%y**2)
!    if (gap.and.((r < rInnerGap).or.(r > rOuterGap))) then
       if ((r > rInner).and.(r < rOuter)) then
          hr = height * (r/rd)**1.125
! Calculate density and check the exponential won't underflow
          if ( rVec%z/hr < 20.0 ) THEN
             thisOctal%rho(subcell) = rho * ((r / rd)**(-1.))*exp(-pi/4.*(rVec%z/hr)**2)
          endif
          thisOctal%rho(subcell) = max(thisOctal%rho(subcell), min_rho)
          thisOctal%temperature(subcell) = 100.
          thisOctal%inFlow(subcell) = .true.
          thisOctal%etaCont(subcell) = 0.
       endif
!    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine benchmarkDisk

  subroutine molecularBenchmark(thisOctal,subcell,grid)

    use input_variables, only : molAbundance !, amr2d

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 50
    real(double),save :: r(nr), nh2(nr), junk,t(nr), v(nr) , mu(nr)
    real(double) :: mu1, r1, t1, t2
    real(double) :: v1 !, vDopp
    integer :: i

    type(VECTOR) :: vel

    if (firsttime) then
       open(31, file="model_1.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
       do i = nr,1,-1                                             
          read(31,*) r(i), nh2(i), junk,t(i), v(i) , mu(i)
       enddo
       r = r * 1.e-10
       close(31)
       firsttime = .false.
    endif

    r1 = modulus(subcellCentre(thisOctal,subcell))
    thisOctal%temperature(subcell) = tcbr
    thisOctal%microTurb(subcell) = 1d-8 !0.159e5/cspeed
    thisOctal%velocity(subcell) = VECTOR(-1d4,-1d4,-1d4)
    thisOctal%molAbundance(subcell) = molAbundance

    if(r1 > r(nr) .or. r1 < r(1)) then 
       thisOctal%nh2(subcell) = 1.e-20
       thisOctal%rho(subcell) = 1.e-20 * 2. * mhydrogen
    else

       call locate(r, nr, r1, i)
       t2 = (r1 - r(i))/(r(i+1)-r(i)) ! linear but know its a power law so use better interpolation

       t1 = log(r1/r(i))/log(r(i+1)/r(i))

       thisOctal%nh2(subcell) = exp((1.d0 - t1) * log(nh2(i))  +  t1 * log(nh2(i+1)))
       
       thisOctal%rho(subcell) = thisOctal%nh2(subcell)*2.d0*mhydrogen

!       thisOctal%temperature(subcell) = t(i) + t1 * (t(i+1)-t(i))
!       thisOctal%nh2(subcell) = nh2(i) + t1 * (nh2(i+1)-nh2(i))
!       thisOctal%rho(subcell) = (nh2(i) + t1 * (nh2(i+1)-nh2(i)))*2.*mhydrogen

       thisOctal%temperature(subcell) = real((t(i) + t2 * (t(i+1)-t(i))))
!       thisOctal%temperature(subcell) = 15. ! debug 30/8/08
!       thisOctal%nh2(subcell) = nh2(i) + t1 * (nh2(i+1)-nh2(i))
!       thisOctal%rho(subcell) = (nh2(i) + t1 * (nh2(i+1)-nh2(i)))*2.*mhydrogen


!       if(amr2d) then

!          thisOctal%nh2(subcell) = thisOctal%nh2(subcell) * costheta**2
!          thisOctal%rho(subcell) = thisOctal%rho(subcell) * costheta**2

!       endif

       v1 = (v(i) + t2 * (v(i+1)-v(i)))*1.d5
       vel = subcellCentre(thisOctal, subcell)
       call normalize(vel)
       thisOctal%velocity(subcell) = (v1 * vel)/cspeed
!       thisOctal%velocity(subcell) = vector(0.d0,0.d0,0.d0) ! debug 17/1/08
       mu1 = mu(i) + t2 * (mu(i+1)-mu(i))

      thisOctal%microturb(subcell) = max(1d-8,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / (29.0 * amu)) + mu1**2) &
                                     / (cspeed * 1e-5)) ! mu is subsonic turbulence
   endif
  end subroutine molecularBenchmark

  subroutine WaterBenchmark1(thisOctal, subcell)

    use input_variables, only : molAbundance !, amr2d

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 200
    real,save :: r(nr)
    real :: r1
    integer :: i

    if (firsttime) then
       open(31, file="grid.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
       do i = 1,nr                                             
          read(31,*) r(i)
       enddo
       r = r * 3.08568025e8 ! (pc -> 10^10cm)
       close(31)
       firsttime = .false.
    endif

    r1 = modulus(subcellCentre(thisOctal,subcell))
    thisOctal%temperature(subcell) = 1d-20
    if(r1 > r(nr) .or. r1 < r(1)) then 
       thisOctal%nh2(subcell) = 1.e-20
       thisOctal%rho(subcell) = 1.e-20 * 2. * mhydrogen
    endif

    thisOctal%microTurb(subcell) = 1d-20
    thisOctal%molAbundance(subcell) = molabundance

    if ((r1 > r(1)).and.(r1 < r(nr))) then

       thisOctal%nh2(subcell) = 1e4
       thisOctal%rho(subcell) = thisOctal%nh2(subcell)*2.*mhydrogen
       thisOctal%temperature(subcell) = 40.
       thisOctal%velocity(subcell) = vector(0d-20,0d-20,0d-20)
       thisOctal%microturb(subcell) = max(1d-9,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) &
                                      / (18.0 * amu)) + 10.d0**2) / (cspeed * 1e-5)) ! mu is subsonic turbulence
    endif
  end subroutine WaterBenchmark1

  subroutine WaterBenchmark2(thisOctal,subcell, grid)

    use input_variables, only : molAbundance !, amr2d

    TYPE(octal), INTENT(INOUT) :: thisOctal
    TYPE(gridtype), INTENT(IN) :: grid
    INTEGER, INTENT(IN) :: subcell
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 200
    real,save :: r(nr)
    real :: r1
    integer :: i

    if (firsttime) then
       open(31, file="grid.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
       do i = 1,nr                                             
          read(31,*) r(i)
       enddo
       r = r * 3.08568025e8 ! (pc -> 10^10cm)
       close(31)
       firsttime = .false.
    endif

    r1 = modulus(subcellCentre(thisOctal,subcell))

    if(r1 > r(nr) .or. r1 < r(1)) then 
       thisOctal%nh2(subcell) = 1.e-20
       thisOctal%rho(subcell) = 1.e-20 * 2. * mhydrogen
    endif

    thisOctal%temperature(subcell) = 40.
    thisOctal%microTurb(subcell) = 1d-10
    thisOctal%molAbundance(subcell) = 1e-20
    thisOctal%velocity(subcell) = vector(1d-20,1d-20,1d-20)

    if ((r1 > r(1)).and.(r1 < r(nr))) then
       thisOctal%molAbundance(subcell) = molabundance
       thisOctal%nh2(subcell) = 1e4
       thisOctal%rho(subcell) = thisOctal%nh2(subcell)*2.*mhydrogen
       thisOctal%velocity(subcell) = WaterBenchmarkvelocity(subcellCentre(thisOctal,subcell),grid)
       thisOctal%microturb(subcell) = max(1d-9,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) &
                                      / (18.0 * amu)) + 10.d0**2) / (cspeed * 1e-5)) ! mu is subsonic turbulence
    endif

   CALL fillVelocityCorners(thisOctal,grid,WaterBenchmarkVelocity,thisOctal%threed)
  end subroutine WaterBenchmark2

  subroutine AGBStarBenchmark(thisOctal,subcell,grid)

    use input_variables, only : molAbundance

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 100
    real,save :: r(nr), nh2(nr), td(nr), tg(nr), v(nr), mu(nr)
    real :: r1!, mu1 , t1, t2
    real(double) :: v1 !, vDopp
    integer :: i

    type(VECTOR) :: vel

    if (firsttime) then
       open(31, file="mc_100.dat", status="old", form="formatted")
       do i = nr,1,-1                                             
          read(31,*) r(i), nh2(i), tg(i), td(i), v(i), mu(i)
       enddo
       r = r * 1e-10
       mu = mu / 2.35482 ! FWHM -> sigma
       close(31)
       firsttime = .false.
    endif

    r1 = modulus(subcellCentre(thisOctal,subcell))

    if(r1 > r(nr)) then
       thisOctal%nh2(subcell) = 1.e-20
       thisOctal%rho(subcell) = 1.e-20 * 2. * mhydrogen
       thisOctal%temperaturegas(subcell) = tcbr
       thisOctal%temperaturedust(subcell) = tcbr
       thisOctal%microturb(subcell) = 1e-8
    elseif(r1 < r(1)) then
       thisOctal%nh2(subcell) = nh2(1)
       thisOctal%rho(subcell) = thisOctal%nh2(subcell) * 2. * mhydrogen
       thisOctal%temperaturegas(subcell) = 3000.
       thisOctal%temperaturedust(subcell) = 3000.
       thisOctal%microturb(subcell) = sqrt((2.d-10 * kerg * thisOctal%temperaturegas(subcell) / (18.0 * amu)) + mu(1)**2) &
                                     / (cspeed * 1e-5)
    endif

    thisOctal%molAbundance(subcell) = molabundance
    thisOctal%velocity(subcell) = vector(1d-20,1d-20,1d-20)

    if ((r1 > r(1)).and.(r1 < r(nr))) then

       thisOctal%nh2(subcell) = 2.4652426912e19*r1**(-2.5)
       thisOctal%rho(subcell) = thisOctal%nh2(subcell)*2.*mhydrogen

       thisOctal%temperaturegas(subcell) = 4.93117864e8 / sqrt(r1*r1*r1) ! 5e8
       thisOctal%temperaturedust(subcell) = 9.001104e6 / r1 ! 9e6 

       if(r1 .gt. (6e0 * 1.49598e3) ) then ! if r gt than 6 AU out this is outward velocity
          v1 = 8.0673550611e-3 * r1**(0.65) * 1d5 ! km -> cm
          thisOctal%microturb(subcell) = sqrt((2.d-10 * kerg * thisOctal%temperaturegas(subcell) / (18.0 * amu)) + 1.d0**2) &
                                         / (cspeed * 1e-5) ! mu is subsonic turbulence
       else
          v1 = 1d-20
          thisOctal%microturb(subcell) = sqrt((2.d-10 * kerg * thisOctal%temperaturegas(subcell) / (18.0 * amu)) + 3.d0**2) &
                                         / (cspeed * 1e-5) ! mu is subsonic turbulence
       endif

       vel = subcellCentre(thisOctal, subcell)
       call normalize(vel)
       thisOctal%velocity(subcell) = (v1 * vel)/cspeed
    endif

    thisOctal%temperature(subcell) = thisOctal%temperaturegas(subcell)

   CALL fillVelocityCorners(thisOctal,grid,AGBStarVelocity,thisOctal%threed)
  end subroutine AGBStarBenchmark

  subroutine iras04158(thisOctal, subcell ,grid)

    use input_variables, only : router, rinner
    use density_mod, only: iras04158Disc
    
    TYPE(octal), INTENT(INOUT) :: thisOctal
    type(GRIDTYPE), intent(in) :: grid
    INTEGER, INTENT(IN) :: subcell
    real :: r
    TYPE(vector) :: rVec
    character(len=10) :: out

    rvec = subcellCentre(thisOctal,subcell)
    r = sqrt(rVec%x**2 + rVec%y**2)

    thisOctal%temperature(subcell) = tcbr 
    thisOctal%nh2(subcell) = 1d-20
    thisOctal%rho(subcell) = thisOctal%nh2(subcell) * 2.*mhydrogen
    thisOctal%molabundance(subcell) = 1d-20
    thisOctal%temperaturedust(subcell) = tcbr
    thisOctal%temperaturegas(subcell) = tcbr
    thisOctal%microturb(subcell) = 1d-10
    thisOctal%Velocity(subcell) = VECTOR(1d-20,1d-20,1d-20)

    if ((r .lt. rOuter) .and. (r .gt. rinner)) then
      
       out ='nh2'
       thisOctal%rho(subcell) = iras04158Disc(rvec)
       thisOctal%nh2(subcell) = thisOctal%rho(subcell) / (2.*mhydrogen)
       
       out ='abundance'
       thisOctal%molabundance(subcell) = readparameterfrom2dmap(rvec,out,.true.)

       out='td'
!       thisOctal%temperaturedust(subcell) = readparameterfrom2dmap(rvec,out,.true.)
       thisOctal%temperaturedust(subcell) = 10.
       thisOctal%temperaturegas(subcell) = thisOctal%temperaturedust(subcell)
       thisOctal%temperature(subcell) = thisOctal%temperaturedust(subcell)
       
       thisOctal%microturb(subcell) = max(1d-8,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) / (28.0 * amu)) + 0.3**2) &
                                      / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence
       thisOctal%velocity(subcell) = keplerianVelocity(rvec, grid)
    else
       thisOctal%velocity = VECTOR(1d-20,1d-20,1d-20)
    endif
!    write(*,*) thisOctal%temperature(subcell), thisOctal%rho(subcell), modulus(thisOctal%velocity(subcell))
   CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
 end subroutine iras04158
 
 function readparameterfrom2dmap(point,output,dologint) result(out)
   
   TYPE(vector), INTENT(IN) :: point
   logical, save :: firsttime = .true.
   integer, parameter :: nr = 140, nz = 60
   real,save :: rgrid(nr), zgrid(nr,nz), nh2(nr,nz), td(nr,nz), abundance(nr,nz)
   real :: junk, zarray(nz,2)
   integer :: r1, r2, z11, z12, z21, z22
   integer :: ri, zi
   real :: rmax, rmin, z1max, z1min, z2min, z2max, r, z, u, v
   real :: zminint, zmaxint, gradmin, gradmax, d
   real(double) :: out
   character(len=10) :: output

   logical, optional :: dologint
   logical :: logint

   if(present(dologint)) then
      logint = dologint
   else
      logint = .false.
   endif

   if (firsttime) then
      open(31, file="molecularmap.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
      
      do ri = 1, nr
         read(31,*) rgrid(ri), zgrid(ri,1), td(ri,1), nh2(ri,1), abundance(ri,1)
         do zi = 2,nz
            read(31,*) junk, zgrid(ri,zi), td(ri,zi), nh2(ri,zi), abundance(ri,zi)
         enddo
      enddo
             
      rgrid = rgrid * 1e-10
      zgrid = zgrid * 1e-10
      close(31)
      firsttime = .false.
   endif

    r = sqrt(point%x**2+point%y**2)
    z = abs(point%z)

    if(r .gt. rgrid(1) .and. r .lt. rgrid(nr)) then
       call locate(rgrid,size(rgrid),r,r1)
    elseif(r .lt. rgrid(1)) then
       r1 = 1
    else
       r1 = nr
    endif
       
    if( (r .gt. rgrid(r1))) then
       r2 = min(r1+1,nr)
    else
       r2 = r1
       r1 = max(1,r2-1)
    endif

    zarray(:,1) = zgrid(r1,:)
    zarray(:,2) = zgrid(r2,:)

    if((z .gt. zarray(1,1) .and. z .lt. zarray(nz,1))) then
       call locate(zarray(:,1),size(zarray(:,1)),z,z11)
    elseif(z .lt. zarray(1,1)) then
       z11 = 1
    else
       z11 = nz
    endif

    if(z .gt. zarray(z11,1)) then
       z12 = min(z11 + 1,nz)
    else
       z12 = z11
       z11 = max(1,z12 - 1)
    endif

    if((z .gt. zarray(1,2)) .and. (z .lt. zarray(nz,2))) then
       call locate(zarray(:,2),size(zarray(:,2)),z,z21)
    elseif(z .lt. zarray(1,2)) then
       z21 = 1
    else
       z21 = nz
    endif

    if( (z .gt. zarray(z21,2))) then
       z22 = min(z21 + 1,nz)
    else
       z22 = z21
       z21 = max(1,z22 - 1)
    endif

     rmax = rgrid(r2)
     rmin = rgrid(r1)
     z1min = zarray(z11,1)
     z1max = zarray(z12,1)
     z2min = zarray(z21,2)
     z2max = zarray(z22,2)

     if(logint) then
        r = log(r)
        rmax= log(rmax) 
        rmin= log(rmin)
        z = log(z)
        z1max= log(z1max) 
        z2max= log(z2max) 
        z1min= log(z1min)
        z2min= log(z2min)
     endif

     if(r1 .ne. r2) then
!        d = sqrt((rmax-rmin)**2+(z2min-z1min)**2)
        u = (r - rmin) / (rmax - rmin)
     else
        u = 0
     endif
     
!     if((zarray(z12,1) .ne. zarray(z11,1)) .and. (zarray(z22,2) .ne. zarray(z21,2)) .and. (r1 .ne. r2)) then
     if((z12 .ne. z11) .and. (z22 .ne. z21) .and. (r1 .ne. r2)) then
        if((z2min-z)*(z2max-z) .gt. 0.d0 .or. (z1min-z)*(z1max-z) .gt. 0.d0) stop
        gradmin = (z2min - z1min) / (rmax - rmin)
        gradmax = (z2max - z1max) / (rmax - rmin)
        zminint = z1min + (r - rmin) * gradmin
        zmaxint = z1max + (r - rmin) * gradmax
        d = zmaxint - zminint
        v = (z - zminint) / d

     elseif((z1max .ne. z1min) .and. (z2max .ne. z2min) .and. (r1 .eq. r2)) then
        zminint = z1min
        zmaxint = z1max
        d = zmaxint - zminint
        v = (z - zminint) / d

     elseif(z .gt. z1max .or. z .lt. z2min) then ! possible for v to be greater than 1 so limit it
        if (rmax /= rMin) then
           gradmin = (z2min - z1min) / (rmax - rmin)
           gradmax = (z2max - z1max) / (rmax - rmin)
        else
           gradmin = 0.
           gradmax = 0.
        endif
        zminint = z1min + (r - rmin) * gradmin
        zmaxint = z1max + (r - rmin) * gradmax
        d = zmaxint - zminint
        if (d /= 0.) then
           v = max(min((z - zminint) / d,2.d0),-1.d0)
        else
           v = 0.d0
        endif

     else

        v = 0

     endif

     if(output .eq. 'abundance') then
        if(z .gt. 1.1 * z1max .and. z .gt. 1.1 * z2max) then
           out = 1e-13
        elseif( .not. logint) then
           out = (1-v) * ((1-u) * abundance(r1,z11) + u * abundance(r2,z21)) +  &
                     v * ((1-u) * abundance(r1,z12) + u * abundance(r2,z22))
        else
           out = 10.d0**((1-v) * ((1-u) * log10(abundance(r1,z11)) + u * log10(abundance(r2,z21))) +  &
                             v * ((1-u) * log10(abundance(r1,z12)) + u * log10(abundance(r2,z22))))
        endif
!        out = v
     
       
     elseif(output .eq. 'nh2') then
        if(z .gt. 1.1 * z1max .and. z .gt. 1.1 * z2max) then
           out = 1e-5
        elseif( .not. logint) then
           out = (1-v) * ((1-u) * nh2(r1,z11) + u * nh2(r2,z21)) + &
                     v * ((1-u) * nh2(r1,z12) + u * nh2(r2,z22))        
        else
           out = 10.d0**((1-v) * ((1-u) * log10(nh2(r1,z11)) + u * log10(nh2(r2,z21))) +  &
                             v * ((1-u) * log10(nh2(r1,z12)) + u * log10(nh2(r2,z22))))
        endif

     elseif(output .eq. 'td') then
        if(z .gt. 1.1 * z1max .and. z .gt. 1.1 * z2max) then
           out = tcbr
        elseif( .not. logint) then
           out = 10.d0**((1-v) * ((1-u) * log10(td(r1,z11)) + u * log10(td(r2,z21))) +  &
                             v * ((1-u) * log10(td(r1,z12)) + u * log10(td(r2,z22))))
        else
           out = (1-v) * ((1-u) * td(r1,z11) + u * td(r2,z21)) + &
                     v * ((1-u) * td(r1,z12) + u * td(r2,z22))        
        endif

     else
        write(*,*) "Undefined parameter"
        out = 1d-20
     endif
end function readparameterfrom2dmap

  subroutine molecularFilamentFill(thisOctal,subcell,grid)

    use input_variables, only : molAbundance!,rhoC,r0

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real(double) :: r0, r1
    TYPE(VECTOR) :: cellCentre

    r0 = (1d3/sqrt(4.d0*pi*6.67d-11*5d-18)) / 1d8 ! 1d3 sigma in m/s, 6.67d-11 is G, 5d-18 is RHOc, 1d8 metres to TORUS

    cellCentre = subcellCentre(thisOctal, subcell)
    r1 = modulus(VECTOR(cellCentre%x, cellCentre%y, 0.d0))

    thisOctal%molAbundance(subcell) = molAbundance
    thisOctal%microTurb(subcell) = 1.e5/cspeed
    thisOctal%temperature(subcell) = 10.d0

    thisOctal%rho(subcell) = 5d-18/(1 + r1**2/(8.d0*r0**2))**2
    thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.*mhydrogen)

    thisOctal%velocity(subcell) = VECTOR(0.d0,0.d0,0.d0)

   CALL fillVelocityCorners(thisOctal,grid,molebenchVelocity,thisOctal%threed)
  end subroutine molecularFilamentFill

  subroutine ggtauFill(thisOctal,subcell,grid)

    use input_variables, only : molAbundance 

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real(double) :: r,t0,nh20,H0,H,z, rmin, rmax
    TYPE(VECTOR) :: cellCentre

    CellCentre = subcellCentre(thisOctal, subcell)
!    r = modulus(cellCentre) ! spherical
    r = sqrt(cellcentre%x**2+cellcentre%y**2) ! cylindrical
    rmin = sqrt((cellcentre%x-thisoctal%subcellsize)**2+(cellcentre%y-thisoctal%subcellsize)**2)
    rmax = sqrt((cellcentre%x+thisoctal%subcellsize)**2+(cellcentre%y+thisoctal%subcellsize)**2)
    r = r * 6.68458134e-06! (torus units to 100's of AUs)
    rmin = rmin * 6.68458134e-06! (torus units to 100's of AUs)
    rmax = rmax * 6.68458134e-06! (torus units to 100's of AUs)

    z = 6.68458134e-04 * Cellcentre%z !in AU

    t0 = 30.
    nh20 = 6.3e9
    H0 = 14.55

    H = H0 * r * sqrt(sqrt(r))  !H0*r^(5/4) r in 100s of AU and H in AU

    thisOctal%molabundance = molAbundance

    if(rmin .gt. 1.8 .and. rmax .lt. 8.) then
       thisOctal%nh2(subcell) = nh20 * (r **(-2.75)) * exp(-((z/H)**2))!(in cm-3) 
       thisOctal%rho(subcell) = thisOctal%nh2(subcell) * 2. * mHydrogen

       thisOctal%temperature(subcell) = t0 / sqrt(r) 
    else
       thisOctal%nh2(subcell) = 1.d-60
       thisOctal%rho(subcell) = 1.d-60 * 2. * mhydrogen
       
       thisOctal%temperature(subcell) = tcbr
    endif

    CALL fillVelocityCorners(thisOctal,grid,ggtauVelocity,thisOctal%threed)

  end subroutine ggtauFill


  TYPE(vector)  function wrshellVelocity(point, grid)
    use input_variables, only : vterm, beta
    type(vector), intent(in) :: point
    type(vector) :: rvec
    type(GRIDTYPE), intent(in) :: grid
    real(double) :: v, r
    rVec = point

    wrShellVelocity = VECTOR(0.d0, 0.d0, 0.d0)
    r = modulus(rVec)
    If ((r > grid%rInner)) then !.and.(r < grid%rOuter)) then
       v = 0.001d5+(vterm-0.001d5)*(1.d0 - grid%rinner/r)**beta
       call normalize(rvec)
       wrshellvelocity = rvec * (v/cSpeed)
    endif

  end function wrshellVelocity

  TYPE(vector) FUNCTION moleBenchVelocity(point, grid)

    type(VECTOR), intent(in) :: point
    type(GRIDTYPE), intent(in) :: grid
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 50
    integer :: i
    real(double), save :: r(nr), nh2(nr), junk,t(nr), v(nr) , mu(nr), OneOverrdiff(nr)
    real(double) :: v1, t1
    real(double) :: r1, pos
    real(double), save :: nstartOverlogdiff, noverlogdiff

    if (firsttime) then
       open(31, file="model_1.dat", status="old", form="formatted")
       do i = nr,1,-1
          read(31,*) r(i), nh2(i), junk,t(i), v(i) , mu(i)
       enddo
       r = r * 1d-10
       close(31)
       firsttime = .false.

       OneOverrdiff(nr) = 1d30
       do i = 1, nr-1
          OneOverrdiff(i) = 1.d0 / (r(i+1)-r(i))
       enddo

       nOverlogdiff = dble(nr-1) / (log(r(nr)) - log(r(1)))
       nstartOverlogdiff = log(r(1)) * dble(nr-1) / (log(r(nr)) - log(r(1)))
    endif

    moleBenchVelocity = VECTOR(0d4,0d4,0d4)
      
    r1 = modulus(point)
    
    if ((r1 > r(1)).and.(r1 < r(nr))) then

!       call locate(r, nr, r1, i)
       
!       t1 = (r1 - r(i)) * OneOverrdiff(i)

       pos = (nOverlogdiff * log(r1) - nstartOverlogdiff)
       
       t1 = fraction(pos)
       i = int(pos) + 1
       v1 = (v(i) + t1 * (v(i+1)-v(i)))*1.d5
       moleBenchVelocity = (v1 * OneOvercSpeed / r1) * point
    endif

  end FUNCTION moleBenchVelocity

  TYPE(vector) FUNCTION WaterBenchmarkVelocity(point, grid)

    type(VECTOR), intent(in) :: point
    type(GRIDTYPE), intent(in) :: grid
    real :: v1, r1
    type(VECTOR) :: vel

    r1 = modulus(point) * 3.24077649e-9 ! 10^10cm -> pc
    WaterBenchmarkVelocity = VECTOR(1d-20,1d-20,1d-20)

    if ((r1 > 0.001) .and. (r1 < 0.1)) then
       v1 = 100. * r1
       vel = point
       call normalize(vel)
       WaterBenchmarkVelocity = (v1/cSpeed) * vel
    endif

  end FUNCTION WaterBenchmarkVelocity

  TYPE(vector) FUNCTION AGBStarVelocity(point,grid)

    type(VECTOR), intent(in) :: point
    type(GRIDTYPE), intent(in) :: grid
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 100
    integer :: i
    real,save :: r(nr)
    real :: v1, r1
    type(VECTOR) :: vel

    if (firsttime) then
       open(31, file="mc_100.dat", status="old", form="formatted")
       do i = nr,1,-1                                             
          read(31,*) r(i)
       enddo
       r = r * 1e-10
       close(31)
       firsttime = .false.
    endif

    r1 = modulus(point)
    AGBStarVelocity = VECTOR(1d-20,1d-20,1d-20)

    if ((r1 > r(1)).and.(r1 < r(nr))) then
       v1 = 8.0673550611E-03 * r1**(0.65)
       vel = point
       call normalize(vel)
       AGBStarVelocity = (v1/cSpeed) * vel
    endif

  end FUNCTION AGBStarVelocity

  TYPE(vector)  function keplerianVelocity(point, grid)
    use input_variables, only : mcore
    type(vector), intent(in) :: point
    type(GRIDTYPE), intent(in) :: grid
    type(vector) :: rvec
    real(double) :: v, r
    rVec = point

    keplerianVelocity = VECTOR(1.d-20, 1.d-20, 1.d-20)

    r = sqrt(point%x**2+point%y**2)

!    write(*,*) "r", r
!    write(*,*) "rvec", rvec, modulus(rvec)

    if ((r > grid%rInner)) then !.and.(r < grid%rOuter)) then
       v = sqrt(2.d0*6.672d-8*mcore/(r*1d10)) ! G in cgs and M in g (from Msun)
!    write(*,*) "v", v
 
       keplerianvelocity = VECTOR(0.d0,0.d0,1.d0) .cross. rvec
      call normalize(keplerianvelocity) 
      keplerianvelocity = (v/cSpeed) * keplerianvelocity        

    endif
  end function keplerianVelocity

  TYPE(vector)  function ggtauVelocity(point, grid)

    type(vector), intent(in) :: point
    TYPE(gridtype), INTENT(IN) :: grid
    type(vector) :: rvec
    real(double) :: v, r,v0

    rvec = point

    r = sqrt(point%x**2+point%y**2)
    r = r * 6.68458134d-06! (torus units to 100's of AUs)

    ggtauvelocity = VECTOR(1.d-20,1.d-20,1.d-20)
    v0 = 3.3d5

    if(r .gt. 1.6d0 .and. r .lt.  10.d0) then
       v = v0 / sqrt(r)
       call normalize(rvec)
       ggtauvelocity = VECTOR(0.d0,0.d0,1.d0) .cross. rvec
       call normalize(ggtauvelocity)
       ggtauvelocity = (v/cSpeed) * ggtauvelocity        
    else
       ggtauvelocity = VECTOR(v0/cspeed,v0/cspeed,v0/cspeed)
    endif    

  end function ggtauVelocity
 
 TYPE(vector) FUNCTION ClusterVelocity(point, grid)

   use sph_data_class, only : clusterparameter

    type(VECTOR), intent(in) :: point
    type(GRIDTYPE), intent(in) :: grid
    integer :: subcell

    call findSubcellLocal(point, recentOctal,subcell)

    clustervelocity = Clusterparameter(point, recentoctal, subcell, theparam = 1) ! use switch for storing velocity
  end FUNCTION ClusterVelocity

  real(double) FUNCTION ClusterDensity(point, grid)
    
    use sph_data_class, only : clusterparameter

    type(VECTOR) :: out
    type(VECTOR), intent(in) :: point
    type(GRIDTYPE), intent(in) :: grid
    integer :: subcell

    call findSubcellLocal(point, recentOctal,subcell)
    
    out = Clusterparameter(point, recentoctal, subcell, theparam = 2) ! use switch for storing velocity
    clusterdensity = out%x
  end FUNCTION ClusterDensity

  subroutine assign_melvin(thisOctal,subcell,grid)

    use density_mod, only: melvinDensity
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = melvinDensity(rVec, grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)

    thisOctal%ionFrac(subcell,1) = 1.e-10
    thisOctal%ionFrac(subcell,2) = 1.
    thisOctal%ionFrac(subcell,3) = 1.e-10
    thisOctal%ionFrac(subcell,4) = 1.       
    thisOctal%etaCont(subcell) = 0.
    if (thisOctal%rho(subcell) < 110.d0*mHydrogen)  then ! in cavity
       thisOctal%dustTypeFraction(subcell, :) = 1.d-20 ! no dust in cavity
       thisOctal%temperature(subcell) = 1.d4 ! 10,000K in cavity
    endif


  end subroutine assign_melvin

  subroutine assign_whitney(thisOctal,subcell,grid)

    use density_mod, only: whitneyDensity
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = whitneyDensity(rVec, grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine assign_whitney

  subroutine assign_planetgap(thisOctal,subcell,grid)
    use density_mod, only: planetgapDensity

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = planetgapDensity(rVec, grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine assign_planetgap

  subroutine assign_toruslogo(thisOctal,subcell,grid)

    use density_mod, only: toruslogoDensity
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = toruslogoDensity(rVec)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine assign_toruslogo


  subroutine shakaraDisk(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, h, rd, ethermal
    TYPE(vector) :: rVec
    
    type(VECTOR),save :: velocitysum
    logical,save :: firsttime = .true.

    if(firsttime) then
       velocitysum = VECTOR(1d-20,1d-20,1d-20)
       firsttime = .false.
    endif

    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)
    thisOctal%inflow(subcell) = .true.
    thisOctal%temperature(subcell) = 10. 
    thisOctal%etaCont(subcell) = 0.
    rd = rOuter / 2.

    thisOctal%rho(subcell) = 1.d-30
    if (associated(thisOctal%nh)) &
         thisOctal%nh(subcell) =  thisOctal%rho(subcell) / mHydrogen
    if (associated(thisOctal%ne)) &
         thisOctal%ne(subcell) = 1.d-5 !thisOctal%nh(subcell)
    if (photoionization) then
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
       thisOctal%nHeI(subcell) = 0.d0 !0.1d0 *  thisOctal%nH(subcell)
    endif
    if (associated(thisOctal%ionFrac)) then
       thisOctal%ionFrac(subcell,1) = 1.
       thisOctal%ionFrac(subcell,2) = 1.e-5
       thisOctal%ionFrac(subcell,3) = 1.
       thisOctal%ionFrac(subcell,4) = 1.e-5
       thisOctal%etaCont(subcell) = 0.
    endif

    if (photoionization) then
       thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
       thisOctal%ne(subcell) = 1.d-5!thisOctal%nh(subcell)
       thisOctal%nhi(subcell) = 1.e-5
       thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    endif

    r = sqrt(rVec%x**2 + rVec%y**2)
    if (r < rOuter) then
       thisOctal%rho(subcell) = density(rVec, grid)
       thisOctal%temperature(subcell) = 10. ! tinkered from 10K - I figured the cooler bits will gently drop but a large no. of cells are close to this temp.
       thisOctal%etaCont(subcell) = 0.
       thisOctal%inflow(subcell) = .true.

       if (photoionization) then
          thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
          thisOctal%ne(subcell) = 1.e-5!thisOctal%nh(subcell)
          thisOctal%nhi(subcell) = 1.e-5
          thisOctal%nhii(subcell) = thisOctal%ne(subcell)
       endif

       h = height * (r / (100.d0*autocm/1.d10))**betaDisc
    
    endif
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    if(molecular) then
 !      if(modulus(rvec) .lt. 1000.) then
          thisOctal%velocity(subcell) = keplerianVelocity(rvec,grid)
          CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
 !      else
 !         thisOctal%velocity(subcell) = vector(0.,0.,0.)
 !      endif
    else
       thisOctal%velocity = VECTOR(0.,0.,0.)
    endif

    if (molecular) then
       thisOctal%microturb(subcell) = sqrt((2.d0*kErg*thisOctal%temperature(subcell))/(29.0 * amu))/cspeed
       thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.*mhydrogen)
    endif

    if (hydrodynamics) then
       thisOctal%rho(subcell) = max(thisOCtal%rho(subcell), 1.e-20)
       thisOctal%temperature(subcell) = (1.d-15*10.d0)/thisOCtal%rho(subcell)
       thisOctal%velocity(subcell) = keplerianVelocity(rvec,grid)
       thisOctal%boundaryCondition(subcell) = 4
       thisOctal%iEquationOfState(subcell) = 1
       thisOctal%phi_i(subcell) = -bigG * mCore / (modulus(rVec)*1.d10)
       if (modulus(rvec) < rInner) then
          thisOctal%phi_i(subcell) = -bigG * mCore / (rInner*1.d10)
          thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
       endif
       thisOctal%gamma(subcell) = 7.d0/5.d0
       ethermal = 1.5d0 * (1.d0/(2.d0*mHydrogen)) * kerg * thisOCtal%temperature(subcell)
       thisOctal%energy(subcell) = ethermal + 0.5d0*(cspeed*modulus(thisOctal%velocity(subcell)))**2
       thisOctal%rhoe(subcell) = thisOctal%energy(subcell) * thisOctal%rho(subcell)
    endif


  end subroutine shakaraDisk

  subroutine warpedDisk(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, h, rd
    TYPE(vector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)
    thisOctal%inflow(subcell) = .true.
    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    rd = rOuter / 2.
    r = sqrt(rVec%x**2 + rVec%y**2)
    if ((r > rInner).and.(r < rOuter)) then
       if (hydroWarp) then
         thisOctal%rho(subcell) = hydroWarpDensity(rVec, grid)
         thisOctal%temperature(subcell) = hydroWarpTemperature(rVec, grid)
       else
         thisOctal%rho(subcell) = density(rVec, grid)
         thisOctal%temperature(subcell) = 20.
       end if
       thisOctal%etaCont(subcell) = 0.
       thisOctal%inflow(subcell) = .true.
       h = height * (r / (100.d0*autocm/1.d10))**betaDisc
    endif


    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine warpedDisk

  ! chris (26/05/04)
  subroutine calcPPDiskDensity(thisOctal, subcell, grid)

    use input_variables
    type(octal), intent(inout) :: thisOctal
    integer, intent(in) :: subcell
    type(gridtype), intent(in) :: grid

    type(vector) :: rVec

    rVec = subcellCentre(thisOctal,subcell)

    thisOctal%rho(subcell) = density(rVec,grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 0.
    thisOctal%chiLine = 0.

  end subroutine calcPPDiskDensity



  subroutine assign_clumpydisc(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector) :: rVec
    real :: r
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.e-33
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .false.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

    thisOctal%rho(subcell) = density(rVec, grid)
   if (thisOctal%rho(subcell) > 1.e-33) thisOctal%inflow(subcell) = .true.
  end subroutine assign_clumpydisc


  !
  ! Recursively deletes the sph_particle list in the grid.
  !
  recursive subroutine delete_particle_lists(thisoctal)
    implicit none
    
    type(octal), pointer :: thisoctal
    type(octal), pointer :: pChild  
    integer              :: i             ! loop counters
    

    IF (ASSOCIATED(thisOctal%gas_particle_list)) then
       DEALLOCATE(thisOctal%gas_particle_list)
    endif
    
    if ( thisOctal%nChildren > 0) then
       do i = 1, thisOctal%nChildren
          pChild => thisOctal%child(i)
          call delete_particle_lists(pChild)
       end do
    end if
    
  end subroutine delete_particle_lists
  
  


  !
  ! Assign the indecies to root node in the grid
  !
  
  subroutine copy_sph_index_to_root(grid)
    use sph_data_class, only: get_npart
    implicit none
    type(gridtype), intent(inout) :: grid    
    !
    integer :: i, n

    ! extracting the number of gas particles in sph
    n = get_npart()

    ! initialize the list
    ALLOCATE(grid%octreeRoot%gas_particle_list(n))

    ! Actually the initial indecies are just 1, 2, 3 ..n
    do i = 1, n
       ! assigin the values to the root
       grid%octreeRoot%gas_particle_list(i) = i
    end do
       
  end subroutine copy_sph_index_to_root

    
  recursive subroutine scaleDensityAMR(thisOctal, scaleFac)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  real :: scaleFac
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call scaleDensityAMR(child, scaleFac)
                exit
             end if
          end do
       else
          thisOctal%rho(subcell) = thisOctal%rho(subcell) * scaleFac
       endif
    enddo
  end subroutine scaleDensityAMR

  recursive subroutine findTotalMass(thisOctal, totalMass, totalMassTrap, minRho, maxRho)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: totalMass
  real(double),optional :: totalMassTrap
  real(double),optional :: minRho, maxRho
  real(double) :: dv
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findtotalMass(child, totalMass, totalmasstrap, minRho, maxRho)
                exit
             end if
          end do
       else

          dv = cellVolume(thisOctal, subcell)
          totalMass = totalMass + (1.d30)*thisOctal%rho(subcell) * dv
          if(present(totalmasstrap)) totalMassTrap = totalMassTrap + (1.d30) * averagerhofromoctal(thisoctal,subcell) * dv

          if (PRESENT(minRho)) then
             if (thisOctal%rho(subcell) > 1.d-40) then
                minRho = min(dble(thisOctal%rho(subcell)), minRho)
             endif
          endif
          if (PRESENT(maxRho)) maxRho = max(dble(thisOctal%rho(subcell)), maxRho)
       endif
    enddo
  end subroutine findTotalMass

!
! Computes average temperature (T_ave) and mass-weighted temperature (T_mass)
!
  subroutine find_average_temperature(grid, T_ave, T_mass, TotalMass)
    implicit none
    type(gridtype), intent(in)  :: grid
    real(double), intent(out) :: T_ave   ! average temperature in [k}
    real(double), intent(out) :: T_mass  ! mass weighted temperature in [K]
    real(double), intent(out) :: TotalMass  ! total mass [g]
    real(double) :: sum_T    ! Sum of temperatures in all cells [K]
    real(double) :: sum_M    ! Sum of mass in all cells [g]
    real(double) :: sum_TM   ! Sum of temperatures*Mass in all cells [K*g]
    integer :: ncell

    ! initialize some values
    sum_T=0.0d0; sum_TM = 0.0d0; ncell=0; sum_M=0.0d0
    
    ! calling the routine in this module to do sum of T and T*M
    call sum_temp_mass(grid%octreeroot, sum_T, sum_TM, sum_M, ncell)

    ! outputs
    T_ave = sum_T/ncell  ![K]
    T_mass = sum_TM/sum_M
    TotalMass = sum_M

  end subroutine find_average_temperature


  ! find the sum of temperature
  recursive subroutine sum_temp_mass(thisOctal, totalTemperature, &
       totalMassTemperature, TotalMass, ncell)
    implicit none
    type(octal), pointer   :: thisOctal
    real(double), intent(inout) :: totalTemperature      ! simple total mass
    real(double), intent(inout) :: totalMassTemperature
    real(double), intent(inout) :: totalMass
    integer, intent(inout) :: ncell ! number of cells used for this calculation
    !
    type(octal), pointer  :: child
    real(double) :: M, T
    integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sum_temp_mass(child, totalTemperature, totalMassTemperature, &
                     totalMass, ncell)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell)) then
             M = (1.d30)*thisOctal%rho(subcell) * thisOctal%subcellSize**3
             T = thisOctal%temperature(subcell)
             totalTemperature = totalTemperature + T
             totalMassTemperature = totalMassTemperature  + M*T
             totalMass = totalMass  + M
             ncell = ncell + 1
          endif
       endif
    enddo
  end subroutine sum_temp_mass



  !
  ! Computes the average, maximum, minimum, contiuum optical depth across each cell, 
  subroutine find_max_min_ave_tau(grid, ilambda, tau_max, tau_min, tau_ave)
    implicit none
    type(gridtype), intent(in)  :: grid
    integer, intent(in)         :: ilambda  ! wavelength index
    real(double), intent(out) :: tau_max
    real(double), intent(out) :: tau_min
    real(double), intent(out) :: tau_ave
    ! 
    integer :: ncell   ! number of cells

    ! initialize some values
    tau_max = -1.0d10; tau_min=1.0d10; tau_ave=0.0d0; ncell=0    

    ! calling the routine in this module to do sum of T and T*M
    call max_min_ave_tau(grid%octreeroot, grid, ilambda, tau_max, tau_min, tau_ave, ncell)


  end subroutine find_max_min_ave_tau


  ! find the average, maximum, minimum, contiuum optical depth across each cell,
  recursive subroutine max_min_ave_tau(thisOctal, grid, ilambda, tau_max, tau_min, tau_ave, ncell)
    implicit none
    type(octal), pointer   :: thisOctal
    type(gridtype), intent(in)  :: grid    
    integer, intent(in)  :: ilambda  ! wavelength index at which tau is computed.
    real(double), intent(inout) :: tau_max
    real(double), intent(inout) :: tau_min
    real(double), intent(inout) :: tau_ave
    integer, intent(inout) :: ncell ! number of cells used for this calculation
    !
    type(octal), pointer  :: child
    real(double) :: chi, tau
    integer :: subcell, i
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call max_min_ave_tau(child, grid, ilambda,  &
                     tau_max, tau_min, tau_ave, ncell)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell)) then
             if (.not.grid%oneKappa) then
                chi = thisOctal%kappaSca(subcell,iLambda) &
                     + thisOctal%kappaAbs(subcell,iLambda)
             else
                chi = ( grid%oneKappaSca(thisOctal%dustType(subcell),iLambda) + &
                     grid%oneKappaAbs(thisOctal%dustType(subcell),iLambda) ) &
                     *thisOctal%rho(subcell)
             endif
             ! chi is normal chi times 10^10cm, and sucellsize is in 10^10 cm.
             ! so tau should be in correct scale here...
             tau = chi*thisOctal%subcellSize
             
             ! now update the max, min and average tau values.
             tau_max = MAX(tau_max, tau)
             tau_min = MIN(tau_min, tau)
             
             ncell = ncell + 1
             
             tau_ave = (tau_ave*dble(ncell-1) + tau)/dble(ncell)
          else
             continue
          endif
       endif
    enddo
  end subroutine max_min_ave_tau





  subroutine spiralWindSubcell(thisOctal, subcell, grid)
    use ostar_mod, only: returnSpiralFactor
    use input_variables
    real(double) :: v, r, rhoOut
    type(VECTOR) :: rVec, rHat
    type(VECTOR) :: octVec
    integer :: subcell
    type(OCTAL) :: thisOctal
    type(GRIDTYPE) :: grid
    character(len=80) :: opacityFilename
    real, allocatable :: rgrid(:), sigma(:), t(:), vgrid(:), eta(:), chi(:), etal(:), chil(:), escat(:)
    integer :: i, j


    rVec = subcellCentre(thisOctal,subcell)
    opacityFilename = "opacity.dat"


    open(20,file=opacityFilename, status="old", form="formatted")
    read(20,*) nr

    allocate(rgrid(1:nr))
    allocate(sigma(1:nr))
    allocate(t(1:nr))
    allocate(vgrid(1:nr))
    allocate(eta(1:nr))
    allocate(etal(1:nr))
    allocate(chi(1:nr))
    allocate(chil(1:nr))
    allocate(escat(1:nr))


!    write(*,*) "Reading opacities from: ",trim(opacityFilename)
    do j = 1, nr
       i = nr - j + 1
       read(20,*) rgrid(i), t(i), sigma(i), vgrid(i), eta(i), chi(i), escat(i), &
            etal(i), chil(i)
    enddo
    close(20)




    r = modulus(rVec)

    thisOctal%inFlow(subcell) = .false.
    thisOctal%velocity(subcell) = VECTOR(0., 0., 0.)
    if (r > grid%rCore) then
       r = modulus(rVec)
       v = v0 + (vTerm - v0) * (1. - grid%rCore/r)**beta
       rhoOut = mdot / (fourPi * (r*1.e10)**2 * v)
       rHat = rVec / r
       octVec = rVec
       thisOctal%rho(subcell) = rhoOut* returnSpiralFactor(octVec, 10.*grid%rCore/real(twoPi), grid)
       thisOctal%temperature(subcell) = Teff
       thisOctal%velocity(subcell)  = (v/cSpeed) * rHat
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaLine(subcell) =  logInterp(etal, nr, rgrid, real(r))
       thisOctal%etaCont(subcell) =  1.e-30
       thisOctal%chiLine(subcell) =  logInterp(chil, nr, rgrid, real(r))
       thisOctal%kappaSca(subcell,1) = logInterp(etal, nr, rgrid,  real(r))/1.e10
       thisOctal%kappaAbs(subcell,1) = 1.e-30
    endif

    CALL fillVelocityCorners(thisOctal,grid,ostarVelocity,thisOctal%threed)

  end subroutine spiralWindSubcell


  TYPE(vector) FUNCTION ostarVelocity(point,grid)
    real(double) :: r1
    type(VECTOR), intent(in) :: point
!    type(VECTOR) rHat
    type(GRIDTYPE), intent(in) :: grid
    real :: v, v0, vterm, beta
    v0 = 100.e5
    vterm = 2500.e5
    beta = 1.


    r1 = modulus(point)

    if (r1 > grid%rCore) then
!       rHat = (1.d0/r1) * point 
       v = v0 + (vTerm - v0) * (1. - grid%rCore/r1)**beta
!       ostarVelocity  = (dble(v/cSpeed)) * rHat
       ostarvelocity = vector(0.,0.,0.)
    endif
  end FUNCTION ostarVelocity


  SUBROUTINE deleteOctreeBranch(thisOctal,onlyChildren,deletedBranch,&
                                adjustParent, grid, adjustGridInfo)
    ! recursively deletes an octal's contents (and any children it has). 
    ! if 'deletedBranch' pointer is supplied, the branch is copied there as it
    !   is deleted. the %parent variable is NOT set here - you must do it 
    !   elsewhere.
    ! if onlyChildren is set, the octal's variables are not changed, only
    !   its children.
    ! *** WARNINGS:                                                   ***
    ! *** This does not delete the top octal itself -                 ***
    ! ***   you must do this yourself after you call this subroutine! ***

    ! HOWTO use this subroutine:
    ! 1. to "unrefine" part of the grid: call without a 'deletedBranch'
    !    and the children of 'thisOctal' will all be removed.
    ! 2. ...

    IMPLICIT NONE

    TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal ! top of branch to be deleted
    LOGICAL, INTENT(IN)            :: onlyChildren ! only delete this octal's *children*
    TYPE(octal), INTENT(INOUT), TARGET, OPTIONAL :: deletedBranch ! optional copy of deleted branch
    LOGICAL, INTENT(IN) :: adjustParent 
      ! whether the physical parameters stored in the parent's subcells
      !   should be filled with data derived from the children being deleted.
    TYPE(gridtype), INTENT(INOUT), OPTIONAL :: grid 
    LOGICAL, INTENT(IN), OPTIONAL :: adjustGridInfo
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell    
      
    INTEGER :: iChild ! loop counter
    
    IF ( PRESENT(deletedBranch) ) THEN

       IF ( PRESENT(adjustGridInfo) ) THEN
         IF ( adjustGridInfo ) THEN
           PRINT *, "Sorry, adjustGridInfo is not supported with"
           PRINT *, "  a deletedBranch in deleteOctreeBranch"
           STOP
         END IF
       END IF

       IF (.NOT. onlyChildren) THEN
          ! need to copy variables from thisOctal to deletedBranch
          CALL copyOctalComponents(source=thisOctal, &
                                   dest=deletedBranch)
          NULLIFY(deletedBranch%parent)
          NULLIFY(deletedBranch%child)
       END IF
       
       ! now need to copy any children from thisOctal to deletedBranch
       IF ( thisOctal%nChildren > 0 ) THEN 
          deletedBranch%child => thisOctal%child
          
          DO iChild = 1, thisOctal%nChildren, 1
             IF ( onlyChildren ) THEN
                NULLIFY(deletedBranch%child(iChild)%parent)
             ELSE
                deletedBranch%child(iChild)%parent => deletedBranch
             END IF
          END DO
          
       END IF ! octal has children
       
       deletedBranch%nChildren = thisOctal%nChildren
       deletedBranch%indexChild = thisOctal%indexChild
       deletedBranch%hasChild = thisOctal%hasChild
       deletedBranch%maxChildren = thisOctal%maxChildren
       
       CALL deleteOctal(thisOctal, deleteChildren=.FALSE.,    &
                        adjustParent=adjustParent, grid=grid, & 
                        adjustGridInfo=.FALSE. )
       
    ELSE ! no deletedBranch
       
       ! can delete everything
       CALL deleteOctal(thisOctal, deleteChildren=.TRUE.,     &
                        adjustParent=adjustParent, grid=grid, & 
                        adjustGridInfo=adjustGridInfo )
       
    END IF ! PRESENT(deletedBranch)
    
    thisOctal%nChildren = 0
    thisOctal%indexChild(:) = -999
    thisOctal%hasChild(:) = .FALSE.
    
  END SUBROUTINE deleteOctreeBranch
  
  SUBROUTINE deleteOctal(thisOctal, deleteChildren,          &
                         adjustParent, grid, adjustGridInfo, &
                         newMaxDepth )
    ! deallocates the variables in an octal.
    ! optionally deallocates the children of the octal.
    
    TYPE(octal), INTENT(INOUT) :: thisOctal
    LOGICAL, INTENT(IN) :: deleteChildren 
    LOGICAL, INTENT(IN) :: adjustParent 
      ! whether the physical parameters stored in the parent's subcells
      !   should be filled with data derived from the children being deleted.
    TYPE(gridtype), INTENT(INOUT), OPTIONAL :: grid 
    LOGICAL, INTENT(IN), OPTIONAL :: adjustGridInfo
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell
    LOGICAL, INTENT(OUT), OPTIONAL :: newMaxDepth ! true if grid depth has changed
      
    LOGICAL :: doAdjustGridInfo
    INTEGER :: maxDeletionDepth
      ! used for tracking the depth in the tree that has been altered   
    INTEGER, PARAMETER :: hugeInt = HUGE(hugeInt)

    doAdjustGridInfo = .FALSE.
    
    IF ( PRESENT(adjustGridInfo) ) THEN
      IF ( adjustGridInfo .AND. (.NOT. (PRESENT(grid))) ) THEN
        PRINT *, "Panic: you must supply a grid if you"
        PRINT *, "       want to 'adjustGridInfo' in deleteOctal"
        STOP
      END IF
      doAdjustGridInfo = adjustGridInfo
    END IF
    
    maxDeletionDepth = -99

    CALL deleteOctalPrivate(thisOctal, deleteChildren,        &
                            adjustParent,                     &
                            adjustGridInfo=doAdjustGridInfo,  &
                            grid=grid, maxDeletionDepth=maxDeletionDepth )

    IF ( doAdjustGridInfo ) THEN
      ! we should see whether the maximum depth of the grid has shrunk

      IF ( .NOT. ( maxDeletionDepth < grid%maxDepth ) ) THEN
        CALL updateMaxDepth(grid,searchLimit=hugeInt, &
                            changeMade=newMaxDepth)
        IF ( newMaxDepth ) CALL setSmallestSubcell(grid)
      END IF
    END IF

  CONTAINS

    RECURSIVE SUBROUTINE deleteOctalPrivate(thisOctal, deleteChildren,      &
                                            adjustParent, adjustGridInfo,   &
                                            grid, maxDeletionDepth )
    
      TYPE(octal), INTENT(INOUT) :: thisOctal
      LOGICAL, INTENT(IN) :: deleteChildren 
      LOGICAL, INTENT(IN) :: adjustParent 
      LOGICAL, INTENT(IN) :: adjustGridInfo
      TYPE(gridtype), INTENT(INOUT), OPTIONAL :: grid 
      INTEGER, INTENT(INOUT), OPTIONAL :: maxDeletionDepth

      INTEGER :: iChild
      INTEGER :: error

      error = 0

      maxDeletionDepth = MAX( thisOctal%nDepth, maxDeletionDepth )

      IF (deleteChildren) THEN

        DO iChild = 1, thisOctal%nChildren, 1
          CALL deleteOctalPrivate(thisOctal%child(iChild),         &
                           deleteChildren=deleteChildren,          &
                           adjustParent=adjustParent, grid=grid,   &
                           adjustGridInfo=adjustGridInfo,          &
                           maxDeletionDepth=maxDeletionDepth )
        END DO
        IF (ASSOCIATED(thisOctal%child)) DEALLOCATE(thisOctal%child)
        IF ( error /= 0 ) CALL deallocationError(error,location=1) 
        NULLIFY(thisOctal%child)

      END IF ! (deleteChildren)

      IF ( adjustParent ) CALL updateParentFromChild(childOctal=thisOctal)

      IF ( adjustGridInfo ) grid%nOctals = grid%nOctals - 1


      IF (ASSOCIATED(thisOctal%ionFrac)) DEALLOCATE(thisOctal%ionFrac,STAT=error)
      IF ( error /= 0 ) CALL deallocationError(error,location=2) 
      NULLIFY(thisOctal%ionFrac)

      IF (ASSOCIATED(thisOctal%photoIonCoeff)) DEALLOCATE(thisOctal%photoIonCoeff,STAT=error)
      IF ( error /= 0 ) CALL deallocationError(error,location=2) 
      NULLIFY(thisOctal%photoIonCoeff)


      IF (ASSOCIATED(thisOctal%kappaAbs)) DEALLOCATE(thisOctal%kappaAbs,STAT=error)
      IF ( error /= 0 ) CALL deallocationError(error,location=2) 
      NULLIFY(thisOctal%kappaAbs)

      IF (ASSOCIATED(thisOctal%kappaSca)) DEALLOCATE(thisOctal%kappaSca,STAT=error)
      IF ( error /= 0 ) CALL deallocationError(error,location=3) 
      NULLIFY(thisOctal%kappaSca)
    
      IF (ASSOCIATED(thisOctal%N)) DEALLOCATE(thisOctal%N,STAT=error)
      IF ( error /= 0 ) CALL deallocationError(error,location=4) 
      NULLIFY(thisOctal%N)

      IF (ASSOCIATED(thisOctal%departCoeff)) DEALLOCATE(thisOctal%departCoeff,STAT=error)
      IF ( error /= 0 ) CALL deallocationError(error,location=5) 
      NULLIFY(thisOctal%departCoeff)

      IF (ASSOCIATED(thisOctal%gas_particle_list)) DEALLOCATE(thisOctal%gas_particle_list,STAT=error)
      IF ( error /= 0 ) CALL deallocationError(error,location=6) 
      NULLIFY(thisOctal%gas_particle_list)

      IF (ASSOCIATED(thisOctal%dustTypeFraction)) DEALLOCATE(thisOctal%dustTypeFraction,STAT=error)
      IF ( error /= 0 ) CALL deallocationError(error,location=7) 
      NULLIFY(thisOctal%dustTypeFraction)

      call deallocateOctalDynamicAttributes(thisOctal)

    END SUBROUTINE deleteOctalPrivate

    SUBROUTINE deallocationError(error,location)
      INTEGER, INTENT(IN) :: error, location
      PRINT *, "DEALLOCATE error in deleteOctal"
      PRINT *, "Error:", error, " Location:", location
      STOP
    END SUBROUTINE deallocationError
    
  END SUBROUTINE deleteOctal


  
  SUBROUTINE insertOctreeBranch(thisOctal,branch,onlyChildren)
    ! adds one octree into another.
    ! use with care, it's not sensible to insert a tree anywhere - this 
    !   subroutine is meant to be used to replace a tree after it
    !   was temporarily removed from the same location.
    ! the %parent variable of thisOctal is not set - you must do this 
    !   elsewhere
    ! if onlyChildren is set, the octal's variables are not changed, only
    !   its children.
    ! the contents of 'branch' are deleted. you should delete 'branch'
    !   itself after calling this subroutine.
    
    ! *** WARNINGS:                                                   ***
    ! *** This does not delete the top octal of the branch itself -   ***
    ! ***   you must do this yourself after you call this subroutine! ***
    ! *** If the octree being changed is part of a 'grid', you'll     ***
    ! ***   have to fix grid%nOctals etc.                             ***

    IMPLICIT NONE

    TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal ! octal where branch is to
                                                    !   be inserted
    TYPE(octal), INTENT(INOUT) :: branch ! branch being inserted
    LOGICAL, INTENT(IN) :: onlyChildren ! only insert the *children* on the
                                        !   branch, leaving the other variables
                                        !   of thisOctal unchanged
    INTEGER :: iChild ! loop counter
    
    IF (ASSOCIATED(thisOctal%child)) THEN
      WRITE(*,*) "Error in insertOctreeBranch, attempt to overwrite existing children"
      STOP
    END IF

    IF (.NOT. onlyChildren) CALL copyOctalComponents(source=branch,dest=thisOctal)

    IF (ASSOCIATED(branch%child)) THEN
      thisOctal%child => branch%child

      DO iChild = 1, SIZE(thisOctal%child), 1
        thisOctal%child(iChild)%parent => thisOctal
      END DO
      thisOctal%nChildren = branch%nChildren
      thisOctal%indexChild = branch%indexChild
      thisOctal%hasChild = branch%hasChild
      thisOctal%maxChildren = branch%maxChildren
      branch%nChildren = 0 
      branch%hasChild(:) = .FALSE. 
    END IF

    CALL deleteOctal(branch, deleteChildren=.FALSE., adjustParent=.FALSE.)
        
  END SUBROUTINE insertOctreeBranch
  
  SUBROUTINE copyOctalComponents(source,dest)
    ! copy the components within an octal variable to a new octal variable
    !
    ! WARNING: this does not change the parent and child variables - you
    !   must update those yourself elsewhere.
 
    TYPE(octal), INTENT(IN) :: source
    TYPE(octal), INTENT(INOUT) :: dest

    ! first make sure that 'dest' is empty
    IF ( ASSOCIATED(dest%child) ) THEN
      PRINT *, "Problem in copyOctalComponents:"
      PRINT *, "destination seems to have some children"
      do;enddo
    END IF
    CALL deleteOctal(dest, deleteChildren=.FALSE., adjustParent=.FALSE. )

    dest%nDepth =  source%nDepth
    dest%mpiThread =  source%mpiThread
    dest%nChildren = source%nChildren
    dest%indexChild = source%indexChild
    dest%threeD =   source%threeD 
    dest%twod = source%twoD  
    dest%oneD = source%oneD   
    dest%iXbitString = source%iXbitString
    dest%iYbitString = source%iYbitString
    dest%iZbitString = source%iZbitString
    dest%cylindrical = source%cylindrical
    dest%splitAzimuthally = source%splitAzimuthally
    dest%maxChildren =  source%maxChildren
    dest%hasChild = source%hasChild
    dest%centre =  source%centre
    dest%rho =  source%rho
    dest%phi =  source%phi
    dest%dphi = source%dphi
    dest%r =   source%r
    dest%velocity = source%velocity
    dest%temperature = source%temperature
    dest%inFlow = source%inFlow
    dest%label =   source%label
    dest%parentSubcell =   source%parentSubcell
    dest%subcellSize =  source%subcellSize
    dest%gasOpacity =  source%gasOpacity 
    dest%cornerVelocity =  source%cornerVelocity

    dest%xMax = source%xMax
    dest%yMax = source%yMax
    dest%zMax = source%zMax


    dest%xMin = source%xMin
    dest%yMin = source%yMin
    dest%zMin = source%zMin

    dest%inStar = source%inStar

    call copyAttribute(dest%nCrossings, source%nCrossings)
    call copyAttribute(dest%chiLine, source%chiLine)
    call copyAttribute(dest%etaLine, source%etaLine)
    call copyAttribute(dest%etaCont, source%etaCont)
    call copyAttribute(dest%biasCont3d, source%biasCont3d)
    call copyAttribute(dest%biasLine3d, source%biasLine3d)
    call copyAttribute(dest%distanceGrid, source%distanceGrid)

    call copyAttribute(dest%probDistLine, source%probDistLine)
    call copyAttribute(dest%probDistCont,  source%probDistCont)
    call copyAttribute(dest%Ne,  source%Ne)
    call copyAttribute(dest%Ntot,  source%Ntot)
    call copyAttribute(dest%NH,  source%NH)
    call copyAttribute(dest%NHI,  source%NHI)
    call copyAttribute(dest%boundaryCondition,  source%boundaryCondition)
    call copyAttribute(dest%boundaryPartner,  source%boundaryPartner)
    call copyAttribute(dest%NHII,  source%NHII)
    call copyAttribute(dest%NHeI,  source%NHeI)
    call copyAttribute(dest%NHeII,  source%NHeII)
    call copyAttribute(dest%Hheating,  source%Hheating)
    call copyAttribute(dest%Heheating,  source%Heheating)
    call copyAttribute(dest%changed,  source%changed)


    call copyAttribute(dest%scatteredIntensity, source%scatteredIntensity)
    call copyAttribute(dest%temperaturedust, source%temperatureDust)
    call copyAttribute(dest%temperaturegas, source%temperaturegas)
    call copyAttribute(dest%dustType, source%dustType)
    call copyAttribute(dest%oldFrac, source%oldFrac)
    call copyAttribute(dest%diffusionApprox, source%diffusionApprox)
    call copyAttribute(dest%diffusionCoeff, source%diffusionCoeff)
    call copyAttribute(dest%undersampled, source%undersampled)
    call copyAttribute(dest%nDiffusion, source%nDiffusion)
    call copyAttribute(dest%nDirectPhotons, source%nDirectPhotons)
    call copyAttribute(dest%oldtemperature, source%oldtemperature)
    call copyAttribute(dest%eDens, source%eDens)
    call copyAttribute(dest%oldeDens, source%oldeDens)
    call copyAttribute(dest%kappaRoss, source%kappaRoss)
 
    call copyAttribute(dest%nh2, source%nh2)
    call copyAttribute(dest%microturb, source%microturb)
    call copyAttribute(dest%molmicroturb, source%molmicroturb)

    call copyAttribute(dest%rho_i_minus_1, source%rho_i_minus_1)

    call copyAttribute(dest%rho_i_plus_1, source%rho_i_plus_1)


    call copyAttribute(dest%x_i_minus_1, source%x_i_minus_1)

    call copyAttribute(dest%x_i, source%x_i)
    call copyAttribute(dest%x_i_plus_1, source%x_i_plus_1)

    call copyAttribute(dest%q_i_minus_2, source%q_i_minus_2)
    call copyAttribute(dest%q_i_minus_1, source%q_i_minus_1)
    call copyAttribute(dest%q_i, source%q_i)
    call copyAttribute(dest%q_i_plus_1, source%q_i_plus_1)

    call copyAttribute(dest%flux_i_minus_1, source%flux_i_minus_1)
    call copyAttribute(dest%flux_i, source%flux_i)
    call copyAttribute(dest%flux_i_plus_1, source%flux_i_plus_1)

    call copyAttribute(dest%pressure_i_minus_1, source%pressure_i_minus_1)
    call copyAttribute(dest%pressure_i, source%pressure_i)
    call copyAttribute(dest%phi_i, source%phi_i)
    call copyAttribute(dest%phi_i_plus_1, source%phi_i_plus_1)
    call copyAttribute(dest%phi_i_minus_1, source%phi_i_minus_1)
    call copyAttribute(dest%pressure_i_plus_1, source%pressure_i_plus_1)

    call copyAttribute(dest%u_interface, source%u_interface)
    call copyAttribute(dest%u_i_plus_1, source%u_i_plus_1)
    call copyAttribute(dest%u_i_minus_1, source%u_i_minus_1)
    call copyAttribute(dest%rLimit, source%rLimit)
    call copyAttribute(dest%phiLimit, source%phiLimit)
    call copyAttribute(dest%ghostCell, source%ghostCell)
    call copyAttribute(dest%edgeCell, source%edgeCell)
    call copyAttribute(dest%feederCell, source%feederCell)
    call copyAttribute(dest%energy, source%energy)
    call copyAttribute(dest%rhou, source%rhou)
    call copyAttribute(dest%rhov, source%rhov)
    call copyAttribute(dest%rhow, source%rhow)
    call copyAttribute(dest%rhoe, source%rhoe)
    call copyAttribute(dest%refinedLastTime, source%refinedLastTime)

    call copyAttribute(dest%photoIonCoeff, source%photoIonCoeff)
    call copyAttribute(dest%molecularLevel, source%molecularLevel)

    call copyAttribute(dest%molAbundance, source%molAbundance)

    call copyAttribute(dest%atomAbundance, source%atomAbundance)

    call copyAttribute(dest%atomLevel, source%atomLevel)

    call copyAttribute(dest%newatomLevel, source%newatomLevel)

    call copyAttribute(dest%ionFrac, source%ionFrac)

    call copyAttribute(dest%gamma, source%gamma)
    call copyAttribute(dest%iEquationOfState, source%iEquationOfState)


    IF (ASSOCIATED(source%mpiboundaryStorage)) THEN                   
      ALLOCATE(dest%mpiboundaryStorage( SIZE(source%mpiboundaryStorage,1),       &
                              SIZE(source%mpiboundaryStorage,2), &
                              SIZE(source%mpiboundaryStorage,3)))
      dest%mpiboundaryStorage = source%mpiboundaryStorage
    END IF  

    IF (ASSOCIATED(source%kappaAbs)) THEN                   
      ALLOCATE(dest%kappaAbs( SIZE(source%kappaAbs,1),       &
                              SIZE(source%kappaAbs,2)))
      dest%kappaAbs = source%kappaAbs
    END IF  

    IF (ASSOCIATED(source%kappaSca)) THEN
      ALLOCATE(dest%kappaSca( SIZE(source%kappaSca,1),       &
                              SIZE(source%kappaSca,2)))
      dest%kappaSca = source%kappaSca
    END IF  

    IF (ASSOCIATED(source%N)) THEN 
      ALLOCATE(dest%N( SIZE(source%N,1),              &
                       SIZE(source%N,2)))
      dest%N = source%N
    END IF  

    IF (ASSOCIATED(source%departCoeff)) THEN
      ALLOCATE(dest%departCoeff( SIZE(source%departCoeff,1),    &
                                 SIZE(source%departCoeff,2)))
      dest%departCoeff = source%departCoeff
    END IF  

    IF (ASSOCIATED(source%gas_particle_list)) THEN
      ALLOCATE(dest%gas_particle_list(SIZE(source%gas_particle_list)))
      dest%gas_particle_list = source%gas_particle_list
    END IF  

    IF (ASSOCIATED(source%dustTypeFraction)) THEN                   
      ALLOCATE(dest%dustTypeFraction( SIZE(source%dustTypeFraction,1), &
                                      SIZE(source%dustTypeFraction,2)))
      dest%dustTypeFraction = source%dustTypeFraction
    END IF  

  END SUBROUTINE copyOctalComponents

  SUBROUTINE checkAMRgrid(grid,checkNoctals)
    ! does some checking that the cells in an AMR grid are
    !   set up and linked to correctly.

    TYPE(gridType), INTENT(IN) :: grid
    LOGICAL, INTENT(IN) :: checkNoctals ! whether to confirm grid%nOctals
    
    TYPE(OCTAL), POINTER :: rootOctal
    INTEGER :: nOctals

    nOctals = 1
    rootOctal => grid%octreeRoot
    
    CALL checkAMRgridPrivate(grid=grid,              &
                             thisOctal=rootOctal,    &
                             thisDepth=1,            &
                             thisParent=NULL(),      &
                             thisParentSubcell=-999, &
                             nOctals=nOctals)

    IF ( checkNoctals .AND. nOctals /= grid%nOctals ) THEN
      PRINT *, "In checkAMRgrid, nOctals mismatch:"
      PRINT *, "  nOctals = ", nOctals
      PRINT *, "  grid%nOctals = ", grid%nOctals
      PRINT *, "Ignoring..."
    END IF

  CONTAINS
  
    RECURSIVE SUBROUTINE checkAMRgridPrivate(grid,thisOctal, thisDepth ,&
                                             thisParent,thisParentSubcell,  &
                                             nOctals)
      TYPE(gridType), INTENT(IN) :: grid
      TYPE(OCTAL), INTENT(IN), TARGET :: thisOctal
      INTEGER, INTENT(IN) :: thisDepth
      TYPE(OCTAL), POINTER :: thisParent
      INTEGER, INTENT(IN) :: thisParentSubcell
      INTEGER, INTENT(INOUT) :: nOctals

      INTEGER :: iSubcell, iIndex
      TYPE(OCTAL), POINTER :: thisOctalPointer
      REAL(oct) :: sizeRatio

      nOctals = nOctals + 1
      
      IF ( thisOctal%nDepth /= thisDepth ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, depth mismatch"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
      END IF
      
      ! check parent (except at root of tree)
      IF ( thisDepth /= 1 ) THEN
        IF ( .NOT. ASSOCIATED(thisOctal%parent,thisParent) ) THEN
          PRINT *, "Error: In checkAMRgridPrivate, mismatch with parent argument"
          CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
        END IF
      END IF
       
      IF ( thisDepth /= 1 ) THEN
         IF ( thisOctal%parentSubcell /= thisParentSubcell ) THEN
            PRINT *, "Error: In checkAMRgridPrivate, parentSubcell mismatch"
            CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
         END IF
      END IF
      
      IF ( thisOctal%nDepth > grid%maxDepth ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, maxDepth exceeded"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
      END IF
      
      IF ( thisOctal%subcellSize < grid%halfSmallestSubcell ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, halfSmallestSubcell exceeded"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
      END IF
       
      IF ( thisOctal%nChildren > thisOctal%maxChildren ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, maxChildren exceeded"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
      END IF

      IF ( ASSOCIATED(thisOctal%child) .AND. (thisOctal%nChildren == 0) ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, thisOctal%child shouldn't be associated"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
      END IF
      
      IF ( .NOT. ASSOCIATED(thisOctal%child) .AND. (thisOctal%nChildren > 0) ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, thisOctal%child should be associated"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
      END IF

      ! check that invalid children are not set:
      IF ( ANY( thisOctal%hasChild(thisOctal%maxChildren+1:SIZE(thisOctal%hasChild)) )) THEN
        PRINT *, "Error: In checkAMRgridPrivate, invalid hasChild variables set"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
      END IF
      
      ! check that the number of children agree 
      IF ( COUNT(thisOctal%hasChild) /= thisOctal%nChildren ) THEN
        PRINT *, "Error: In checkAMRgridPrivate, nChildren does not match hasChild"
        CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
      END IF

      IF ( thisOctal%nChildren > 0 ) THEN

        ! see if %child is sized correctly
        IF ( SIZE(thisOctal%child) /= thisOctal%nChildren ) THEN
          PRINT *, "Error: In checkAMRgridPrivate, %child has wrong size"
          CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
        END IF

        DO iSubcell = 1, thisOctal%maxChildren

          IF ( .NOT. thisOctal%hasChild(iSubcell) ) CYCLE
          
          ! find the correct index of %child for this child
          DO iIndex = 1, thisOctal%nChildren
            IF ( thisOctal%indexChild(iIndex) == iSubcell ) EXIT
            
            IF ( iIndex == thisOctal%nChildren ) THEN
              ! shouldn't get here
              PRINT *, "Error: In checkAMRgridPrivate, indexChild not found"
              PRINT *, "       for iSubcell = ",iSubcell
              CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals) 
            END IF
          END DO

          ! now we know the correct %child variable, let's do some tests

          IF ( thisOctal%child(iIndex)%parentSubcell /= iSubcell ) THEN
            PRINT *, "Error: In checkAMRgridPrivate, child's parentSubcell doesn't match:"
            PRINT *, "       thisOctal%child(iIndex)%parentSubcell = ",thisOctal%child(iIndex)%parentSubcell
            PRINT *, "       iIndex = ", iIndex
            PRINT *, "       iSubcell = ", iSubcell
            CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
          END IF

          ! see if the child's coordinates really lie in the parent subcell
          IF ( .NOT. inSubcell(thisOctal,iSubcell,point=thisOctal%child(iIndex)%centre) ) THEN
            PRINT *, "Error: In checkAMRgridPrivate, child isn't in parentSubcell"
            PRINT *, "       thisOctal%centre = ",thisOctal%centre 
            PRINT *, "       iSubcell = ", iSubcell
            PRINT *, "       subcellCentre(thisOctal,iSubcell) = ",subcellCentre(thisOctal,iSubcell)
            PRINT *, "       thisOctal%child(iIndex)%centre = ",thisOctal%child(iIndex)%centre 
            PRINT *, "       iIndex = ", iIndex
            CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
          END IF

          ! see if the child is of the correct size
          sizeRatio = thisOctal%subcellSize / thisOctal%child(iIndex)%subcellSize
          sizeRatio = sizeRatio / 2.0_oc
          IF ( ABS(sizeRatio-1.0_oc) > 0.1 ) THEN
            PRINT *, "Error: In checkAMRgridPrivate, size of child wrong:"
            PRINT *, "       thisOctal%subcellSize = ", thisOctal%subcellSize
            PRINT *, "       thisOctal%child(iIndex)%subcellSize = ", thisOctal%child(iIndex)%subcellSize
          END IF

          ! check the child's parent pointer
          IF ( .NOT. ASSOCIATED(thisOctal%child(iIndex)%parent,thisOctal) ) THEN
            PRINT *, "Error: In checkAMRgridPrivate, child has wrong %parent"
            CALL printErrorPrivate(grid,thisOctal,thisDepth,nOctals)
          END IF

          ! call recursively on this child
          thisOctalPointer => thisOctal
          CALL checkAMRgridPrivate(grid,                              &
                                   thisOctal=thisOctal%child(iIndex), &
                                   thisDepth=thisDepth+1,             &
                                   thisParent=thisOctalPointer,       &
                                   thisParentSubcell=iSubcell,        &
                                   nOctals=nOctals)
                                   
        END DO                           
                                   
      END IF ! ( thisOctal%nChildren > 0 )
      
    END SUBROUTINE checkAMRgridPrivate
    
    SUBROUTINE printErrorPrivate(grid,thisOctal,thisDepth,nOctals)  
        
      TYPE(gridType), INTENT(IN) :: grid
      TYPE(OCTAL), INTENT(IN) :: thisOctal
      INTEGER, INTENT(IN) :: thisDepth
      INTEGER, INTENT(IN) :: nOctals

      PRINT *, "  thisOctal%nDepth = ", thisOctal%nDepth
      PRINT *, "  thisDepth = ", thisDepth
      PRINT *, "  grid%maxDepth = ", grid%maxDepth
      PRINT *, "  grid%halfSmallestSubcell = ", grid%halfSmallestSubcell
      PRINT *, "  thisOctal%subcellSize = ", thisOctal%subcellSize
      PRINT *, "  thisOctal%nChildren = ", thisOctal%nChildren
      PRINT *, "  thisOctal%maxChildren = ", thisOctal%maxChildren
      PRINT *, "  thisOctal%twoD = ", thisOctal%twoD
      PRINT *, "  thisOctal%threeD = ", thisOctal%threeD
      PRINT *, "  thisOctal%hasChild = ", thisOctal%hasChild
      PRINT *, " ASSOCIATED(thisOctal%child) = ", ASSOCIATED(thisOctal%child) 
      IF ( ASSOCIATED(thisOctal%child) ) THEN
        PRINT *, "  SIZE(thisOctal%child) = ", SIZE(thisOctal%child)
      END IF
      PRINT *, "  grid%nOctals = ", grid%nOctals
      PRINT *, "  nOctals (counted so far) = ", nOctals
      !STOP
      PRINT *, "Entering infinite loop..."
      DO ; END DO
      
    END SUBROUTINE printErrorPrivate       
    
  END SUBROUTINE checkAMRgrid

  ! chris (19/05/04)
  ! Smooths the AMR grid to ensure 'adequate' resolution at the boundary
  ! between optically thin and optically thick regions. This routine should be
  ! called repeatedly until gridConverged returns true (it should be set to
  ! true in the initial call.
  ! The value of kappa is the opacity of the grid at a test wavelength (set in
  ! the parameter file).

  RECURSIVE SUBROUTINE smoothAMRgridTau(thisOctal,grid,gridConverged, ilam, &
                                        stellar_cluster, &
                                        inheritProps, interpProps, romData)
    USE input_variables, ONLY : tauSmoothMax,tauSmoothMin
    IMPLICIT NONE

    TYPE(octal), POINTER             :: thisOctal
    TYPE(gridtype), INTENT(INOUT   ) :: grid 
    LOGICAL, INTENT(INOUT)               :: gridConverged
    TYPE(cluster), optional, intent(in)  :: stellar_cluster
    LOGICAL, INTENT(IN) :: inheritProps
    LOGICAL, INTENT(IN), optional :: interpProps
    TYPE(romanova), optional, INTENT(IN)   :: romData  ! used for "romanova" geometry

    INTEGER              :: i, ilam
    TYPE(vector)    :: thisSubcellCentre
    REAL(oct) :: dSubcellCentre
    real(double) :: kappaAbs, kappaSca
    TYPE(octal), POINTER :: child
    TYPE(octal), POINTER :: neighbour
    TYPE(vector), ALLOCATABLE, DIMENSION(:) :: locator
    INTEGER              :: subcell
    INTEGER              :: thisSubcell
    REAL                 :: thisTau, thatTau !, tauDiff
    INTEGER              :: nLocator
    LOGICAL              :: prob

    kappaAbs = 0.d0; kappaSca = 0.d0
    ! For each subcell, we find the coordinates of at least one point in every
    ! neighbouring subcell. The optical depths are compared and if one cell is
    ! optically thick while the other is optically thin with an optical depth
    ! less than some value, both cells are split. There may be up to 26
    ! neighbouring subcells.

    if (thisOctal%threed) then
       nLocator = 26
    else
       nlocator = 8
    endif

    ALLOCATE(locator(nLocator))

    thisSubcell = 1
    do
      if (thisOctal%hasChild(thisSubcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == thisSubcell) then
              child => thisOctal%child(i)
              call smoothAMRgridTau(child,grid,gridConverged,  ilam, stellar_cluster, &
                   inheritProps, interpProps, romData)
              ! force return until algorithm is fixed
              IF ( .NOT. gridConverged ) RETURN
              exit
           end if
        end do
      else
        thisSubcellCentre = subcellCentre(thisOctal,thisSubcell)

        ! don't split outer edge of disc

        if ((grid%geometry == "shakara").or.(grid%geometry == "warpeddisc").or.(grid%geometry == "iras04158").or.&
             (grid%geometry=="circumbin")) then
           if (sqrt(thissubcellcentre%x**2 + thissubcellcentre%y**2) > grid%router*0.9) goto 666
        endif
           

        FORALL (i = 1:nLocator)
          locator(i) = thisSubcellCentre
        END FORALL

        ! Moving this distance in any direction will take us into a
        ! neighbouring subcell.
        dSubcellCentre = (0.5 * thisOctal%subcellSize) + grid%halfSmallestSubcell

        if (thisOctal%threed) then 

           ! faces
           locator(1)%x = thisSubcellCentre%x + dSubcellCentre
           locator(2)%y = thisSubcellCentre%y + dSubcellCentre
           locator(3)%z = thisSubcellCentre%z + dSubcellCentre
           locator(4)%x = thisSubcellCentre%x - dSubcellCentre
           locator(5)%y = thisSubcellCentre%y - dSubcellCentre
           locator(6)%z = thisSubcellCentre%z - dSubcellCentre
           ! x-edges
           locator(7)%y = thisSubcellCentre%y + dSubcellCentre
           locator(7)%z = thisSubcellCentre%z + dSubcellCentre
           locator(8)%y = thisSubcellCentre%y + dSubcellCentre
           locator(8)%z = thisSubcellCentre%z - dSubcellCentre
           locator(9)%y = thisSubcellCentre%y - dSubcellCentre
           locator(9)%z = thisSubcellCentre%z + dSubcellCentre
           locator(10)%y = thisSubcellCentre%y - dSubcellCentre
           locator(10)%z = thisSubcellCentre%z - dSubcellCentre
           ! y-edges
           locator(11)%x = thisSubcellCentre%x + dSubcellCentre
           locator(11)%z = thisSubcellCentre%z + dSubcellCentre
           locator(12)%x = thisSubcellCentre%x + dSubcellCentre
           locator(12)%z = thisSubcellCentre%z - dSubcellCentre
           locator(13)%x = thisSubcellCentre%x - dSubcellCentre
           locator(13)%z = thisSubcellCentre%z + dSubcellCentre
           locator(14)%x = thisSubcellCentre%x - dSubcellCentre
           locator(14)%z = thisSubcellCentre%z - dSubcellCentre
           ! z-edges
           locator(15)%x = thisSubcellCentre%x + dSubcellCentre
           locator(15)%y = thisSubcellCentre%y + dSubcellCentre
           locator(16)%x = thisSubcellCentre%x + dSubcellCentre
           locator(16)%y = thisSubcellCentre%y - dSubcellCentre
           locator(17)%x = thisSubcellCentre%x - dSubcellCentre
           locator(17)%y = thisSubcellCentre%y + dSubcellCentre
           locator(18)%x = thisSubcellCentre%x - dSubcellCentre
           locator(18)%y = thisSubcellCentre%y - dSubcellCentre
           ! corners
           locator(19)%x = thisSubcellCentre%x + dSubcellCentre
           locator(19)%y = thisSubcellCentre%y + dSubcellCentre
           locator(19)%z = thisSubcellCentre%z + dSubcellCentre
           locator(20)%x = thisSubcellCentre%x + dSubcellCentre
           locator(20)%y = thisSubcellCentre%y + dSubcellCentre
           locator(20)%z = thisSubcellCentre%z - dSubcellCentre
           locator(21)%x = thisSubcellCentre%x + dSubcellCentre
           locator(21)%y = thisSubcellCentre%y - dSubcellCentre
           locator(21)%z = thisSubcellCentre%z + dSubcellCentre
           locator(22)%x = thisSubcellCentre%x + dSubcellCentre
           locator(22)%y = thisSubcellCentre%y - dSubcellCentre
           locator(22)%z = thisSubcellCentre%z - dSubcellCentre
           locator(23)%x = thisSubcellCentre%x - dSubcellCentre
           locator(23)%y = thisSubcellCentre%y + dSubcellCentre
           locator(23)%z = thisSubcellCentre%z + dSubcellCentre
           locator(24)%x = thisSubcellCentre%x - dSubcellCentre
           locator(24)%y = thisSubcellCentre%y + dSubcellCentre
           locator(24)%z = thisSubcellCentre%z - dSubcellCentre
           locator(25)%x = thisSubcellCentre%x - dSubcellCentre
           locator(25)%y = thisSubcellCentre%y - dSubcellCentre
           locator(25)%z = thisSubcellCentre%z + dSubcellCentre
           locator(26)%x = thisSubcellCentre%x - dSubcellCentre
           locator(26)%y = thisSubcellCentre%y - dSubcellCentre
           locator(26)%z = thisSubcellCentre%z - dSubcellCentre
        else

           ! edges

           locator(1)%x = thisSubcellCentre%x + dSubcellCentre
           locator(2)%z = thisSubcellCentre%z + dSubcellCentre
           locator(3)%x = thisSubcellCentre%x - dSubcellCentre
           locator(4)%z = thisSubcellCentre%z - dSubcellCentre
           ! corners
           locator(5)%x = thisSubcellCentre%x + dSubcellCentre
           locator(5)%z = thisSubcellCentre%z + dSubcellCentre

           locator(6)%x = thisSubcellCentre%x + dSubcellCentre
           locator(6)%z = thisSubcellCentre%z - dSubcellCentre

           locator(7)%x = thisSubcellCentre%x - dSubcellCentre
           locator(7)%z = thisSubcellCentre%z + dSubcellCentre

           locator(8)%x = thisSubcellCentre%x - dSubcellCentre
           locator(8)%z = thisSubcellCentre%z - dSubcellCentre


        endif

        call returnKappa(grid, thisOctal, thissubcell, ilambda=ilam,&
             kappaSca=kappaSca, kappaAbs=kappaAbs)

        thisTau = (kappaAbs + kappaSca) * thisOctal%subcellSize !* thisOctal%rho(thisSubcell)

        DO i = 1, nLocator, 1
          IF ( inOctal(grid%octreeRoot,locator(i)) ) THEN
            neighbour => thisOctal
            CALL findSubcellLocal(locator(i),neighbour,subcell, prob=prob)
                IF ( neighbour%hasChild(subcell) ) THEN 
                  PRINT *, "neighbour already has child (a)",prob
                  do;enddo
                END IF
            
            ! The neighbouring subcell must be larger than the current subcell
            ! otherwise our locators won't cover all neighbouring subcells
            ! (and we'll hit cell boundaries, which is not good).
            IF ( neighbour%subcellSize >= thisOctal%subcellSize) THEN

               call returnKappa(grid, neighbour, subcell, &
                    ilambda=ilam,  kappaSca=kappaSca, kappaAbs=kappaAbs)
               thatTau = (kappaSca + kappaAbs) * neighbour%subcellSize! * neighbour%rho(subcell)

              ! Original critera was to split if:
              ! ((o1 - o2)/(o1+o2) > 0.5 .and. (o1 - o2) > 1)
!              tauDiff = abs(thisTau - thatTau)
!!              if ((tauDiff.gt.1000000.0).and. &
!              if ((tauDiff.gt.1.).and. &
!!                  .not.((thisTau.gt.5.).and.(thatTau.gt.5.)).and. &
!                  ((tauDiff/(thisTau + thatTau)).gt.0.5)) then

              ! Critera for cell splitting is that one cell is optically think
              ! (tau > 1) and the other is optically thin with an optical
              ! depth less than some value (which should really be set in the
              ! parameter file).
!              if (max(thisTau, thatTau).gt.0.5.and.min(thisTau, thatTau).lt.0.01) then
!              if ((max(thisTau, thatTau).gt.tauSmoothMax.and.min(thisTau, thatTau).lt.tauSmoothMin).or.&
!                   ((max(thisTau, thatTau)<5.).and.(abs(thistau-thatTau)> 2.))) then

               if ((max(thisTau, thatTau).gt.tauSmoothMax.and.min(thisTau, thatTau).lt.tauSmoothMin)) then


!write (*,*) thisSubcellCentre%x, thisSubcellCentre%y, thisSubcellCentre%z, thisTau, thatTau, i
                ! Because addNewChild is broken, we must use addNewChildren below,
                ! which makes these tests on the current and neighbouring subcells
                ! insufficient. We perform them anyway.
                IF ( neighbour%hasChild(subcell) ) THEN 
                  PRINT *, "neighbour already has child"
                  STOP
                END IF
                IF ( thisOctal%hasChild(thisSubcell) ) THEN 
                  PRINT *, "thisOctal already has child"
                  STOP
                END IF
                ! Note that the version of addNewChild called here does not
                ! support sphData, etc.
!                call addNewChild(neighbour,subcell,grid)
!                call addNewChild(thisOctal,thisSubcell,grid)

              call addNewChild(neighbour, subcell, grid, adjustGridInfo=.TRUE.,  &
                               stellar_cluster=stellar_cluster, &
                               inherit=inheritProps, interp=interpProps, romData=romData)
!              call addNewChild(thisOctal, thissubcell, grid, adjustGridInfo=.TRUE.,  &
!                               sphData=sphData, stellar_cluster=stellar_cluster, &
!                               inherit=inheritProps, interp=interpProps, romData=romData)
!              write(*,*) "added new child",thisOctal%hasChild,thisOctal%nchildren


!                call addNewChildren(thisOctal, grid, sphData, stellar_cluster, inheritProps, interpProps, romData)
                ! If we are splitting two subcells in the same octal, we don't
                ! need to addNewChildren to the neighbour. If we switch to
                ! addNewChild this test becomes redundant.
!                if (.not. associated(thisOctal, neighbour)) then
!                  call addNewChildren(neighbour, grid, sphData, stellar_cluster, inheritProps, interpProps, romData)
!                end if 

                gridConverged = .FALSE.
                if (.not.gridConverged) return!!!!
                exit ! loop through locators, move onto next subcell
              end if ! grid must be refined
            ENDIF ! neighbour subcell is larger
          END IF ! in grid
        END DO ! loop through locators

      end if ! thisOctal%hasChild(thisSubcell)
      thisSubcell = thisSubcell + 1 ! next subcell
      if (thisSubcell > thisOctal%maxChildren) exit ! loop through subcells in an octal
    end do ! loop through subcells in an octal
666 continue
    deallocate(locator)
  END SUBROUTINE smoothAMRgridTau

  subroutine polarDump(grid)
    implicit none
    type(GRIDTYPE) :: grid
    integer :: nt = 200, nr = 400, i, j
    real :: theta, r, t
    real(double) :: dens
    type(octal), pointer   :: thisOctal
    integer :: subcell
    type(VECTOR) :: rVec

    open(30, file="polardump.dat", status="unknown", form="formatted")
    do j = 1, nt
       do i = 1, nr
          r = grid%rInner + (grid%rOuter-grid%rInner)*(real(i-1)/real(nr-1))**3
          theta = pi*real(j-1)/real(nt-1)
          rVec = VECTOR(dble(r*sin(theta)),0.d0,dble(r*cos(theta)))
          call amrGridValues(grid%octreeRoot, rVec,temperature=t, rho=dens, foundoctal=thisOctal, foundsubcell=subcell)
          if (thisOCtal%inFlow(subcell)) then
             write(30,*) r/1.5e3, t, dens
          endif
       enddo
    enddo
    close(30)

    open(30, file="vert1.dat", status="unknown", form="formatted")
    call verticalDump(grid%octreeRoot, grid%rInner*2.)
    close(30)

    open(30, file="vert2.dat", status="unknown", form="formatted")
    call verticalDump(grid%octreeRoot, grid%rOuter/2.)
    close(30)
    open(30,file="radial.dat",status="unknown",form="formatted")
    do i = 1, nr
       r = grid%rInner + (grid%rOuter-grid%rInner)*(real(i-1)/real(nr-1))**3
       rVec = VECTOR(dble(r), 0.d0, 0.d0)
       call amrGridValues(grid%octreeRoot, rVec,temperature=t, rho=dens, foundoctal=thisOctal, foundsubcell=subcell)
       write(30,*) r/1496., t, dens
    enddo
    close(30)


  end subroutine polarDump



  recursive subroutine verticalDump(thisOctal, xPos)
    real :: xPos
    type(VECTOR) :: rVec
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call verticalDump(child, xPos)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          if ( (xPos >= rVec%x-thisOctal%subcellsize/2.).and. &
               (xPos <= rVec%x+thisOctal%subcellsize/2.)) then
             write(30,*) rVec%z/1496., thisOctal%temperature(subcell),thisOctal%rho(subcell)
          endif
       endif
          
    enddo

  end subroutine verticalDump


  recursive subroutine unrefineThickCells(thisOctal, grid, ilambda, converged)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: ilambda
    real(double) :: kappaAbs, kappaSca, tau
    integer :: subcell, i
    logical :: unrefine, converged
    kappaAbs =0.d0; kappaSca = 0.d0
    unrefine = .true.

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unrefineThickCells(child, grid, ilambda, converged)
                exit
             end if
          end do
       else
          if (.not.ASSOCIATED(thisOctal%dustTypeFraction)) then
             write(*,*) "unalloc dusttypefraction!!"
          endif
          call returnKappa(grid, thisOctal, subcell, ilambda, kappaAbs=kappaAbs,kappaSca=kappaSca)
          tau = thisOctal%subcellSize*(kappaAbs+kappaSca)
          if (tau < 1.e4) unrefine = .false.
       endif
    enddo

    if ((thisOctal%nChildren == 0).and.unrefine.and.converged) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
       converged = .false.
    endif

  end subroutine unrefineThickCells

  SUBROUTINE shrinkChildArray(parent, childrenToDelete, adjustParent )
    ! removes children from an octal.
    ! you probably don't want to call this directly - use the 'deleteChild' wrapper instead.
    ! NB this subroutine doesn't update grid%nOctals etc.

    IMPLICIT NONE
    
    TYPE(octal), TARGET, INTENT(INOUT) :: parent ! the parent octal 
    LOGICAL, INTENT(IN), DIMENSION(parent%maxChildren) :: childrenToDelete
      ! mask defining which children to get rid of 
      ! NB childrenToDelete does not map to the index number in the 
      !   %child array, it is the "real" number of the children (the 
      !   number of the subcell that was refined when the child was
      !   created).
    LOGICAL, INTENT(IN) :: adjustParent 
      ! whether the physical parameters stored in the parent's subcells
      !   should be filled with data derived from the children being deleted.
    
    TYPE(wrapperArray) :: tempChildStorage  
      ! holder for remaining children, while we shrink the %child array 
                                       
    LOGICAL, DIMENSION(parent%maxChildren) :: checkMask 
      ! used for testing the validity of the 'childrenToDelete' input
      
    INTEGER, DIMENSION(SIZE(parent%indexChild)) :: temporaryIndexChild
      ! used for assembling a valid indexChild array, to be used after
      !   the children have been deleted
      
    LOGICAL, DIMENSION(parent%maxChildren) :: deleteMask
      ! which of the elements of %child will be deleted
        
    INTEGER :: iChild ! loop counter
    INTEGER :: nChildren ! number of children the parent octal has
    INTEGER :: error
    INTEGER :: nChildrenToDelete
    TYPE(octal), POINTER :: thisChild ! convenient alias to current child 
    LOGICAL :: deleteALLchildren ! are we getting rid of ALL the children
    INTEGER :: nChildrenStay ! how many children will be left when done
    INTEGER :: insertLocation ! the next location to use in tempChildStorage
    
    NULLIFY(thisChild)
    temporaryIndexChild = -999
    
    ! setup some useful accounting variables
    nChildren = parent%nChildren
    nChildrenToDelete = COUNT(childrenToDelete)
    nChildrenStay = ( nChildren - nChildrenToDelete )
    deleteALLchildren = ( nChildrenStay == 0 )
    
    ! some safety checks
    error = 0
    
    IF ( nChildrenStay < 0 ) error = -1

    ! the following lines check that all the children to be deleted have 
    !   their %hasChild flag set.
    checkMask = childrenToDelete .AND. parent%hasChild(1:SIZE(childrenToDelete))
    checkMask = checkMask .NEQV. childrenToDelete ! (exclusive OR operation)
    IF ( ANY(checkMask) ) error = -2

    IF (error /= 0) THEN
      PRINT *, "In shrinkChildArray, attempting to delete a "
      PRINT *, "child that doesn't exist."
      PRINT *, error, childrenToDelete
      STOP
    END IF
 
    ! we transform 'childrenToDelete' which is of SIZE(1:maxChildren) into
    !   an array of SIZE(1:nChildren), so that it matches the %child
    !   array
    deleteMask(:) = .FALSE.
    FORALL ( iChild = 1:parent%nChildren ) &
      deleteMask(iChild) = childrenToDelete(parent%indexChild(iChild))
    
    IF ( .NOT. deleteALLchildren ) THEN
      ! we need to allocate some temporary storage for the children that 
      !   are not going to be deleted so that we can reposition them
      !   in the %child array.
      ALLOCATE(tempChildStorage%wrappers(nChildrenStay), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed in shrinkChildArray. (A)'
        STOP
      END IF
      
      FORALL ( iChild = 1:nChildrenStay ) &
        tempChildStorage%wrappers(iChild)%inUse = .FALSE.
        
      insertlocation = 1
        
    END IF  
      
    DO iChild = 1, parent%nChildren
      thisChild => parent%child(iChild)

      IF ( deleteMask(iChild) ) THEN
        ! we want to delete this octal from the %child array

        CALL deleteOctal(thisChild, deleteChildren=.TRUE.,     &
                         adjustParent=adjustParent )

     ELSE 
        ! we do not want to delete this child. 
        ! instead we move it to the temporary storage array.

        ALLOCATE(tempChildStorage%wrappers(insertlocation)%content, STAT=error)
        IF ( error /= 0 ) THEN
          PRINT *, 'Panic: allocation failed in shrinkChildArray. (B)'
          STOP
        END IF           
        tempChildStorage%wrappers(insertlocation)%inUse = .TRUE.
               
        CALL deleteOctreeBranch(thisOctal=thisChild,                           &
               onlyChildren=.FALSE.,                                           &
               deletedBranch=tempChildStorage%wrappers(insertlocation)%content,&
               adjustParent=.FALSE.)
               
        temporaryIndexChild(insertLocation) = parent%indexChild(iChild)
        
        insertlocation = insertlocation + 1

      END IF
      
    END DO
   
    ! all the unwanted children have now been deleted.
    ! we now get rid of the old %child array.
    DEALLOCATE(parent%child)
    NULLIFY(parent%child) 
    
    ! if there are any remaining children, we need to create a new
    !   %child array, and copy them back there.
    IF ( .NOT. deleteALLchildren ) THEN    
      
      ALLOCATE(parent%child(nChildrenStay), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed in shrinkChildArray. (C)'
        STOP
      END IF

      DO iChild = 1, nChildrenStay, 1
        
        CALL insertOctreeBranch(parent%child(iChild),               &
               branch=tempChildStorage%wrappers(iChild)%content,    &
               onlyChildren=.FALSE.)                              
               
        parent%child(iChild)%parent => parent
               
        DEALLOCATE(tempChildStorage%wrappers(iChild)%content)
        NULLIFY(tempChildStorage%wrappers(iChild)%content)
        tempChildStorage%wrappers(iChild)%inUse = .FALSE.

      END DO
      
      ! can now clean up the temporary storage
      DEALLOCATE(tempChildStorage%wrappers)
      NULLIFY(tempChildStorage%wrappers)

    END IF

    ! some bookkeeping
    parent%nChildren = nChildrenStay
    parent%indexChild(:) = temporaryIndexChild(:)

    parent%hasChild( 1:SIZE(childrenToDelete) ) =  &
         parent%hasChild(1:SIZE(childrenToDelete)) .NEQV. childrenToDelete

  END SUBROUTINE shrinkChildArray
  
    
  SUBROUTINE deleteChild(parent, childToDelete, adjustParent, &
                         grid, adjustGridInfo)
    ! removes one child from an octal.
    ! can be used to "unrefine" part of the grid.
    ! this is a simpler wrapper for the shrinkChildArray subroutine

    IMPLICIT NONE
    
    TYPE(octal), INTENT(INOUT) :: parent ! the parent octal 
    INTEGER :: childToDelete ! number of the child to delete
      ! NB childToDelete does not map to the index number in the 
      !   %child array, it is the "real" number of the child (the 
      !   number of the subcell that was refined when the child was
      !   created).
    LOGICAL, INTENT(IN) :: adjustParent 
      ! whether the physical parameters stored in the parent's subcells
      !   should be filled with data derived from the children being deleted.
    TYPE(gridtype), INTENT(INOUT) :: grid 
    LOGICAL, INTENT(IN) :: adjustGridInfo
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell

    LOGICAL, DIMENSION(parent%maxChildren) :: deletionMask
    INTEGER, PARAMETER :: hugeInt = HUGE(hugeInt)
    LOGICAL :: newMaxDepth ! true if grid depth has changed
    
    IF ( (childToDelete > parent%maxChildren) .OR.   &
         (childToDelete < 0)                      ) THEN
      PRINT *, "Invalid child number passed to deleteChild: ", childToDelete 
      STOP
    END IF
    
    deletionMask(:) = .FALSE.
    deletionMask(childToDelete) = .TRUE.
    
    CALL shrinkChildArray(parent=parent,                 &
                          childrenToDelete=deletionMask, &
                          adjustParent=adjustParent)

    IF ( adjustGridInfo ) THEN
      grid%nOctals = grid%nOctals - 1

      IF ( ( parent%nChildren == 0 ) .AND.              &
           ( grid%maxDepth == (parent%nDepth+1) ) ) THEN
        ! we should see whether the maximum depth of the grid has shrunk
        CALL updateMaxDepth(grid,searchLimit=hugeInt, &
                            changeMade=newMaxDepth)
        IF ( newMaxDepth ) CALL setSmallestSubcell(grid)
      END IF
    END IF  

    ! only for debugging - comment out later:
    CALL checkAMRgrid(grid,checkNoctals=.FALSE.)
        
      
  END SUBROUTINE deleteChild

  SUBROUTINE updateParentFromChild(childOctal)
    ! uses the 4 or 8 subcells in an octal to compute the physical
    !   variables in the appropriate subcell of that octal's parent.

    TYPE(OCTAL), INTENT(INOUT) :: childOctal 
    
    TYPE(OCTAL), POINTER :: parentOctal
    INTEGER :: parentSubcell
    INTEGER :: nVals, i
    REAL(double) :: nValsREAL

    IF ( childOctal%nDepth == 1 ) THEN
      ! we're at the root of the tree
      PRINT *, "Warning: Attempted to update parent, but already at tree root."
      PRINT *, "         Ignoring this and continuing..."
      RETURN
    END IF
    
    parentOctal => childOctal%parent
    parentSubcell = childOctal%parentSubcell
    
    if (childOctal%oneD) then
       nVals = 2
    else if (childOctal%twoD) then
       nVals = 4
    else if (childOctal%threeD) then
       nVals = 8
    endif

    nValsREAL = REAL( nVals, KIND=double )
    
    parentOctal%rho(parentSubcell) =                    &
    SUM(childOctal%rho(1:nVals)) / nValsREAL

    if (associated(parentOctal%nh)) parentOctal%nh(parentSubcell) =                    &
    SUM(childOctal%nh(1:nVals)) / nValsREAL

    if (associated(parentOctal%ne)) parentOctal%ne(parentSubcell) =                    &
    SUM(childOctal%ne(1:nVals)) / nValsREAL

    if (associated(parentOctal%rhoe)) parentOctal%rhoe(parentSubcell) =  &
    SUM(childOctal%rhoe(1:nVals)) / nValsREAL

    if (associated(parentOctal%rhou))     parentOctal%rhou(parentSubcell) = &
    SUM(childOctal%rhou(1:nVals)) / nValsREAL
 
    if (associated(parentOctal%rhov))     parentOctal%rhov(parentSubcell) = &
    SUM(childOctal%rhov(1:nVals)) / nValsREAL

    if (associated(parentOctal%rhow))     parentOctal%rhow(parentSubcell) = &
    SUM(childOctal%rhow(1:nVals)) / nValsREAL
   
    parentOctal%temperature(parentSubcell) = &
     SUM(childOctal%temperature(1:nVals)) / nValsREAL


    if (associated(parentOctal%ionFrac)) then
       do i = 1, SIZE(parentOctal%ionFrac,2)
          parentOctal%ionFrac(parentSubcell,i) =            &
               SUM(childOctal%ionFrac(1:nVals,i)) / nValsREAL
       enddo
    endif
    
    if (associated(parentOctal%dustTypeFraction)) then
       do i = 1, SIZE(childOCtal%dustTypeFraction,2)
          parentOctal%dustTypeFraction(parentSubcell,i) =     &
               SUM(childOctal%dustTypeFraction(1:nVals,i)) / nValsREAL
       enddo
    endif

  END SUBROUTINE updateParentFromChild
  
  SUBROUTINE updateMaxDepth(grid,searchLimit,changeMade)
    ! if octals have been deleted from the grid, we may have to 
    !   check whether the grid%maxDepth variable should be updated
  
    TYPE(gridtype), INTENT(INOUT)  :: grid
    INTEGER, INTENT(IN), OPTIONAL :: searchLimit 
      ! the depth at which to abandon the search.
      ! e.g. if octals have been deleted from levels 2 and 3, pass 
      !   "3" in here, and the search can be abandoned if an octal
      !   of depth 4 is encountered anywhere
    LOGICAL, INTENT(OUT), OPTIONAL :: changeMade
      ! true if grid%maxDepth is actually changed

    TYPE(octal), POINTER  :: thisOctal
    INTEGER :: maxDepthFound
    INTEGER :: oldMaxDepth
    LOGICAL :: searchLimitReached

    oldMaxDepth = grid%maxDepth
    thisOctal => grid%octreeRoot
    maxDepthFound = 0
    searchLimitReached = .FALSE.
    
    CALL updateMaxDepthPrivate(thisOctal,maxDepthFound,       &
                               searchLimitReached,searchLimit)
    IF ( searchLimitReached .OR. ( maxDepthFound == oldMaxDepth ) ) THEN
      IF ( PRESENT(changeMade) ) changeMade = .FALSE.
    ELSE
      IF ( PRESENT(changeMade) ) changeMade = .TRUE.
      grid%maxDepth = maxDepthFound
    END IF
    
    CONTAINS
    
      RECURSIVE SUBROUTINE updateMaxDepthPrivate(thisOctal,maxDepthFound,&
                                          searchLimitReached,searchLimit)
      
        TYPE(OCTAL), INTENT(IN) :: thisOctal 
        INTEGER, INTENT(INOUT) :: maxDepthFound
        LOGICAL, INTENT(INOUT) :: searchLimitReached
        INTEGER, INTENT(IN), OPTIONAL :: searchLimit 
        INTEGER :: iChild
        
        maxDepthFound = MAX( maxDepthFound, thisOctal%nDepth )
        
        IF ( thisOctal%nChildren > 0 ) THEN

          IF ( PRESENT(searchLimit) ) THEN
            IF ( (thisOctal%nDepth+1) > searchLimit ) THEN
              searchLimitReached = .TRUE.
              RETURN
            END IF
          END IF
          
          ! call this subroutine recursively on each of its children
          DO iChild = 1, thisOctal%nChildren, 1
            CALL updateMaxDepthPrivate(thisOctal%child(iChild),maxDepthFound,&
                                          searchLimitReached,searchLimit)
                                          
            IF ( searchLimitReached ) RETURN
          END DO
        END IF

      END SUBROUTINE updateMaxDepthPrivate

  END SUBROUTINE updateMaxDepth

  
  SUBROUTINE setSmallestSubcell(grid)
    ! calculates and stores the grid%halfSmallestSubcell value.
    ! grid%maxDepth must already be set correctly.
  
    TYPE(gridtype), INTENT(INOUT) :: grid 
    
    ! we actually store the value which is half the size of the 
    !   smallest subcell because this is more useful for later
    !   calculations.
    grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
                                2.0_oc**REAL(grid%maxDepth,kind=oct)

  END SUBROUTINE setSmallestSubcell

  RECURSIVE SUBROUTINE setAllUnchanged(thisOctal)
    ! goes through an octree and sets all the %changed variables
    !   to .FALSE.
  
    TYPE(octal), INTENT(INOUT) :: thisOctal 
    INTEGER :: iChild
    
    thisOctal%changed(:) = .FALSE.
    
    DO iChild = 1, thisOctal%nChildren
      CALL setAllUnchanged(thisOctal%child(iChild))
    END DO

  END SUBROUTINE setAllUnchanged


  recursive subroutine biasCentreVolume(thisOctal, boxSize)
    real :: boxSize
    type(VECTOR) :: rVec
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call biasCentreVolume(child, boxSize)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          if ( ((rVec%x < boxSize/2.).and.(rVec%x > -boxsize/2.)) .and. &
               ((rVec%y < boxSize/2.).and.(rVec%y > -boxsize/2.)) .and. &
               ((rVec%z < boxSize/2.).and.(rVec%z > -boxsize/2.)) ) then
             thisOctal%biasCont3D(subcell) = 100.
          endif
       endif
          
    enddo

  end subroutine biasCentreVolume


  recursive subroutine setbiasAMR(thisOctal, grid)
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  real :: r
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setBiasAMR(child, grid)
                exit
             end if
          end do
       else
          r = modulus(subcellcentre(thisOctal, subcell)) / grid%rInner
          thisOctal%biasCont3D(subcell) = r
       endif
    enddo
  end subroutine setBiasAMR

  recursive subroutine set_bias_shakara(thisOctal, grid, ilam, ross)
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: tau, kappaAbs, kappaSca !, thisTau
!  type(vector) :: rVec, direction
  logical :: ross
  integer :: subcell, i, ilam !, ntau
!  real :: r
  kappaAbs = 0.d0; kappaSca = 0.d0
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call set_bias_shakara(child, grid, ilam, ross)
                exit
             end if
          end do
       else
          if (.not.ross) then
             call returnKappa(grid, thisOctal, subcell, ilam, kappaAbs=kappaAbs, kappaSca=kappaSca)
             tau = thisOctal%subcellSize*(kappaAbs + kappaSca)
             thisOctal%biasCont3D(subcell) = exp(-tau)
          else
!             rVec = subcellCentre(thisOctal, subcell)
!             tau = 1.d30
!             ntau = 20
!             do i = 1, nTau
!                direction = randomUnitVector()
!                call tauAlongPath(grid, rVec, direction, thistau, 100.d0)
!                tau = min(tau, thisTau)
!             enddo
!             write(*,*) tau
             thisOctal%biasCont3D(subcell) = 1.d0/(cellVolume(thisOCtal,subcell)*thisOctal%etaCont(subcell))
          endif
!          r = rVec%x
!          thisOCtal%biasCont3d(subcell) = thisOctal%biasCont3d(subcell) * r
       endif

    enddo
  end subroutine set_bias_shakara

  recursive subroutine set_bias_rosseland(thisOctal, grid)
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: tau, kappaAbs
  integer :: subcell, i, ilam
  ilam = 0; kappaAbs = 0.d0
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call set_bias_rosseland(child, grid)
                exit
             end if
          end do
       else
          call returnKappa(grid, thisOctal, subcell, ilam, rosselandKappa = kappaAbs)
          tau = thisOctal%subcellSize*kappaAbs*thisOctal%rho(subcell)*1.d10
          thisOctal%biasCont3D(subcell) = MAX(exp(-tau),1.d-8)
       endif
    enddo
  end subroutine set_bias_rosseland

  recursive subroutine set_bias_radius(thisOctal, p)
  type(octal), pointer  :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  integer :: p
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call set_bias_radius(child, p)
                exit
             end if
          end do
       else
          thisOctal%biasCont3D(subcell) = thisOctal%biasCont3D(subcell) &
               * modulus(subcellCentre(thisOctal, subcell))**p
       endif
    enddo
  end subroutine set_bias_radius

  recursive subroutine set_bias_whitney(thisOctal, grid)
  use input_variables, only : erOuter
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: tau, kappaAbs
  integer :: subcell, i, ilam
  real(double) :: r
  kappaAbs = 0.d0; ilam = 0
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call set_bias_whitney(child, grid)
                exit
             end if
          end do
       else
          r = modulus(subcellCentre(thisoctal,subcell))
          call returnKappa(grid, thisOctal, subcell, ilam, rosselandKappa = kappaAbs)
             tau = thisOctal%subcellSize*kappaAbs*thisOctal%rho(subcell)*1.d10
             thisOctal%biasCont3D(subcell) = MAX(exp(-tau),1.d-4) * (1.d10*r/erOuter)
!             thisOctal%biasCont3D(subcell) = 1.d0/(thisOctal%etacont(subcell)*cellVolume(thisOctal, subcell))
       endif
    enddo
  end subroutine set_bias_whitney

  recursive subroutine set_bias_ttauri(thisOctal, grid, lambda0, outVec)
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  real, intent(in)       :: lambda0                ! rest wavelength of line
  type(octal), pointer  :: child 
  integer :: subcell, i
  real(double) :: d, dV, r, nu0, tauSob, escProb
!  real(double) :: tauSob_tmp, tauSob_sum, theta, phi, sin_theta
  type(vector)  :: rvec, rhat !, rvec_sample
!  integer :: nsample
  type(VECTOR) :: outVec
  real(double):: dtau_cont, dtau_line, rho
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call set_bias_ttauri(child, grid, lambda0, outVec)
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
                dv = 2.0_db*pi*d*d*rVec%x
             endif

             nu0  = cSpeed_dbl / dble(lambda0*angstromtocm)

             tauSob = thisOctal%chiline(subcell)  / nu0

!             ! Find the average tauSob over all azimuth positions
!             nsample = 40
!             tauSob_sum = 0.0d0
!             do i = 1, nsample
!                phi = 2.0d0*PIdouble*dble(i-1)/dble(nsample-1)
!                theta = ACOS(rhat%z); sin_theta = SIN(theta)
!                rvec_sample = VECTOR(r*cos(phi)*sin_theta, r*sin(phi)*sin_theta,  rvec%z)
!                tauSob_tmp = tauSob / amrGridDirectionalDeriv(grid, rvec_sample, outVec, &
!                  startOctal=thisOctal)
!                tauSob_sum = tauSob_sum + tauSob_tmp
!             end do
!             tauSob = tauSob_sum/dble(nsample)

             ! in obs direction, but works only for i=0 case.
!             tauSob = tauSob / amrGridDirectionalDeriv(grid, rvec, outVec, &
!                  startOctal=thisOctal)
             ! in a radial direction
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
             escProb = max(escProb, 1.d-5)
!             escProb = min(1.d-2,max(escProb, 1.d-5))


             dtau_cont = d*(thisOctal%kappaAbs(subcell,1) + thisOctal%kappaSca(subcell,1))
             dtau_line = d*(thisOctal%chiline(subcell))  / nu0
             !
             rho = thisOctal%rho(subcell)
!             thisOctal%biasCont3D(subcell) = MAX(dV, 1.d-7) ! Limits the minimum value
!             thisOctal%biasLine3D(subcell) = dV*thisOctal%biasCont3D(subcell)
             thisOctal%biasCont3D(subcell) = MAX(EXP(-dtau_cont), 1.d-7) ! Limits the minimum value
!             thisOctal%biasCont3D(subcell) = 1.0d0
!             thisOctal%biasLine3D(subcell) = thisOctal%rho(subcell)*thisOctal%biasCont3D(subcell)
!             thisOctal%biasLine3D(subcell) = escProb*thisOctal%biasCont3D(subcell)
             thisOctal%biasLine3D(subcell) = 1.0d0/SQRT(rho) * thisOctal%biasCont3D(subcell)
!             thisOctal%biasLine3D(subcell) = 1.0d0/(rho**0.25)
!             thisOctal%biasLine3D(subcell) = EXP(-dtau_line)*thisOctal%biasCont3D(subcell)
!             thisOctal%biasLine3D(subcell) = EXP(-dtau_line*d*rVec%x)*thisOctal%biasCont3D(subcell)

          else  ! this subcell is not "inFlow"
             thisOctal%biasCont3D(subcell) = 1.0d-150
             thisOctal%biasLine3D(subcell) = 1.0d-150
          end if
          ! just in case ....
          thisOctal%biasCont3D(subcell) = MAX(thisOctal%biasCont3D(subcell), 1.0d-150)
          thisOctal%biasLine3D(subcell) = MAX(thisOctal%biasLine3D(subcell), 1.0d-150)       
       endif ! if (thisOctal%hasChild(subcell)) then
    enddo
  end subroutine set_bias_ttauri


  recursive subroutine set_bias_ttauri2(thisOctal, grid)
    use input_variables, only: TTauriRinner, TTauriRouter
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    real(double) :: r
    type(vector)  :: rvec
    real(double)::  rM, theta,  h, rM_center
  
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call set_bias_ttauri2(child, grid)
                exit
             end if
          end do
       else
          if (thisOctal%inflow(subcell)) then
             
             rVec = subcellCentre(thisOctal,subcell)
             r = modulus(rvec)
             
             if (r/=0.0d0) then
                theta = ACOS(MIN(ABS(rVec%z/r),0.998_oc))
             else
                theta=0.01
             end if

             if (theta == 0.0d0) theta = 0.01
             rM  = r / SIN(theta)**2

             ! bias more toward the edges of the accreation stream
             h = 0.5d0*(TTauriRouter - TTauriRinner)*1.0d-10           ! [10^10cm]   
             rM_center = 0.5d0*(TTauriRouter + TTauriRinner)*1.0d-10   ![10^10cm]   
             
             thisOctal%biasCont3D(subcell) = 1.0d0  ! no bias for contiuum
!             thisOctal%biasLine3D(subcell) = 1.0d0/(dtau_line*d*rvec%x)
!             thisOctal%biasLine3D(subcell) = 1.0d0/thisOctal%rho(subcell)
!             thisOctal%biasLine3D(subcell) = EXP(10.0d0*ABS(rM - rM_center)/h)
!             if (ABS(rM - rM_center)/h > 0.80d0) then
             if (ABS(rM - rM_center)/h > 0.40d0) then
                thisOctal%biasLine3D(subcell) = 1.0d5
             else
                thisOctal%biasLine3D(subcell) = 1.0d-150
             end if
          else  ! this subcell is not "inFlow"
             thisOctal%biasCont3D(subcell) = 1.0d-150
             thisOctal%biasLine3D(subcell) = 1.0d-150
          end if
          ! just in case ....
          thisOctal%biasCont3D(subcell) = MAX(thisOctal%biasCont3D(subcell), 1.0d-150)
          thisOctal%biasLine3D(subcell) = MAX(thisOctal%biasLine3D(subcell), 1.0d-150)       
       endif ! if (thisOctal%hasChild(subcell)) then
    enddo
  end subroutine set_bias_ttauri2


  recursive subroutine set_bias_cmfgen(thisOctal, grid, lambda0)
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real, intent(in)          :: lambda0                ! rest wavelength of line
  integer :: subcell, i
  real(double) :: d, dV, r, tauSob, escProb
  type(vector)  :: rvec, rhat
  real(double):: dtau_cont, dtau_line, nu0, xc, xl, nr, dr, dA
  
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
             escProb = max(escProb, 1.d-7)
             
             dtau_cont = d*(thisOctal%kappaAbs(subcell,1) + thisOctal%kappaSca(subcell,1))
             dtau_line = d*(thisOctal%chiline(subcell))  / nu0
             xc = dv*(thisOctal%kappaAbs(subcell,1) + thisOctal%kappaSca(subcell,1))
             xl = dv*(thisOctal%chiline(subcell))  / nu0
             nr=dA/d**2  ! Approximate number of cells (for spherical case) in dr stripe
             !
!             thisOctal%biasCont3D(subcell) = dV
!             thisOctal%biasLine3D(subcell) = dV*dV
             thisOctal%biasCont3D(subcell) = MAX( (EXP(-dtau_cont)), 1.d-7) ! Limits the minimum value
!             thisOctal%biasLine3D(subcell) = (escProb*nr)*thisOctal%biasCont3D(subcell)
             thisOctal%biasLine3D(subcell) = EXP(-dtau_line)*thisOctal%biasCont3D(subcell)

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




  recursive subroutine setBiasPpdisk(thisOctal, grid)
!  use input_variables, only: rInner
  type(gridtype) :: grid
  type(VECTOR) :: rVec
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
!  real(double) :: r
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setBiasPpdisk(child, grid)
                exit
             end if
          end do
       else
!          r = modulus(subcellcentre(thisOctal, subcell)) / rInner
!          thisOctal%biasCont3D(subcell) =  sqrt(r)
!!          if ((r > 1.).and.(r < 1.05)) then
!          if ((r > 1.).and.(r < 2.)) then
!             thisOctal%biasCont3D(subcell) =  thisOctal%biasCont3D(subcell) * 10.d0
!          endif

          rVec = subcellcentre(thisOctal, subcell)
          thisOctal%biasCont3D(subcell) = rVec%x**2.


          if (thisOctal%diffusionApprox(subcell)) then
             thisOctal%biasCont3D(subcell) = thisOctal%biasCont3D(subcell) * 1.e-2
          endif
       endif
    enddo
  end subroutine setBiasPpdisk


  recursive subroutine massInAnn(thisOctal, r1, r2, mass)
  type(VECTOR) :: rVec
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: mass, r1, r2, r
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call massInAnn(child,  r1, r2, mass)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          r  = sqrt(rVec%x**2 + rVec%y**2)
          if ( (r > r1).and.(r < r2) ) then
             mass = mass + thisOctal%rho(subcell)*thisOctal%subCellSize**3 
          endif
       endif
    enddo
  end subroutine massInAnn


  SUBROUTINE calcWindTestValues(thisOctal,subcell,grid) 
    ! calculates some of the variables for a spherical wind flow

    USE constants_mod

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(INOUT) :: grid

    TYPE(vector) :: point

    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    TYPE(vector) :: pointVecNorm
    TYPE(Vector)      :: pointVecNormSingle

    REAL :: r, rStar
    REAL :: velocity, rho
    REAL, PARAMETER :: v0         = 100.e5
    REAL, PARAMETER :: vTerminal  = 2000.e5
    REAL, PARAMETER :: tEff       = 30.e3
    REAL, PARAMETER :: mDot       = 1.e-7 * mSol * secsToYears 
!    REAL, PARAMETER :: mDot       = 1.e-2 * mSol * secsToYears 
    REAL, PARAMETER :: vBeta      = 1.0
    

    rStar = grid%rStar1
    starPosn = grid%starPos1

    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn) 
    
    r = modulus( pointVec ) 
    pointVecNorm = pointVec 
    CALL normalize(pointVecNorm)

    IF (r > grid%rInner .AND. r < grid%rOuter) THEN

       ! calculate the velocity
       velocity = REAL((v0 + (vTerminal - v0) * (1. - rStar / r)**vBeta),kind=oct)

       ! calculate the density
       rho = mDot / ( fourPi * (r*1.e10)**2 * velocity)

       ! store the data 
       thisOctal%inFlow(subcell) = .TRUE.
       thisOctal%temperature(subcell) = 0.8 * tEff
       pointVecNormSingle = pointVecNorm
       thisOctal%velocity(subcell) = (velocity / cSpeed) * pointVecNormSingle
       thisOctal%rho(subcell) = rho

    ELSE
       thisOctal%inFlow(subcell) = .FALSE.
       thisOctal%rho(subcell) = 1.e-20
       thisOctal%temperature(subcell) = 0.8 * tEff
       thisOctal%velocity(subcell) = vector(1.e-25,1.e-25,1.e-25)
    END IF
    
    IF (subcell == 8) CALL fillVelocityCorners(thisOctal,grid,windTestVelocity,thisOctal%threed)
      
  END SUBROUTINE calcWindTestValues

  
  TYPE(vector) FUNCTION windTestVelocity(point,grid)

    IMPLICIT NONE

    TYPE(vector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    
    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    TYPE(vector) :: pointVecNorm
    TYPE(Vector)      :: pointVecNormSingle

    REAL :: r, rStar
    REAL :: velocity
    REAL, PARAMETER :: v0         = 100.e5
    REAL, PARAMETER :: vTerminal  = 2000.e5
    REAL, PARAMETER :: vBeta      = 1.0
    
    rStar = grid%rStar1
    starPosn = grid%starPos1

    pointVec = (point - starPosn) 
    
    r = modulus( pointVec ) 
    pointVecNorm = pointVec 
    CALL normalize(pointVecNorm)
    pointVecNormSingle = pointVecNorm
    
    IF (r > grid%rInner .AND. r < grid%rOuter) THEN
       velocity = REAL((v0 + (vTerminal - v0) * (1. - rStar / r)**vBeta),kind=oct)
       windTestVelocity = (velocity / cSpeed) * pointVecNormSingle
    ELSE
       windTestVelocity = vector(1.e-25,1.e-25,1.e-25)
    END IF
    
  END FUNCTION windTestVelocity

  

  RECURSIVE SUBROUTINE getIntersectedOctals(thisOctal,listHead,grid,nOctals,onlyChanged)  
    ! returns a linked list of all the octals which are intersected by a
    !   plane. we assume that the z-axis lies in the plane and the plane
    !   is at 45deg from the x and y axes. we only need to consider half
    !   the plane (one quadrant of the simulation space).
    !   
    ! if onlyChanged is True, we ignore cells without the Changed flag

    USE vector_mod
    
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT), TARGET :: thisOctal
    TYPE(octalListElement), POINTER :: listHead
    TYPE(gridType), INTENT(IN)      :: grid
    INTEGER,INTENT(INOUT)           :: nOctals   ! number of octals
    LOGICAL, INTENT(IN), OPTIONAL   :: onlyChanged
    
    LOGICAL               :: octalAdded
    TYPE(vector)     :: centrePoint
    TYPE(vector)     :: lineOrigin
    TYPE(vector)     :: starPos
    TYPE(vector)     :: nearestPlanePoint
    TYPE(octal), POINTER  :: child
    INTEGER               :: i, iSubcell, subIndex
    TYPE(vector)     :: planeVec
    TYPE(octalListElement), POINTER :: newElement
    
    octalAdded = .FALSE.
    starPos = grid%starPos1
    lineOrigin = starPos
    planeVec = vector(1.0_oc/SQRT(2.0_oc),1.0_oc/SQRT(2.0_oc),0.0_oc)

    DO iSubcell = 1, thisOctal%maxChildren, 1 
    
      centrePoint = subcellCentre(thisOctal,iSubcell) - starPos
      lineOrigin%z = centrePoint%z
      centrePoint%z = 0.0_oc
      
      IF (centrePoint%x < 0.0 .OR. centrePoint%y < 0.0) THEN
        ! we are in the wrong quadrant.
        ! this makes some assumptions about the numerical accuracy of the
        !   grid, but it should be mostly OK.
        CYCLE 
      END IF

      nearestPlanePoint = lineOrigin + (planeVec * (centrePoint .dot. planeVec))
    
      IF (inSubcell(thisOctal,iSubcell,nearestPlanePoint)) THEN
                            
        IF (thisOctal%hasChild(iSubcell)) THEN
           ! call the recursive subroutine on its child
           
          subIndex = -99
          DO i = 1, thisOctal%nChildren, 1
            IF ( thisOctal%indexChild(i) == iSubcell ) THEN
              subIndex = i  
              EXIT
            ENDIF
          ENDDO
          IF (subIndex == -99) THEN
            PRINT *, ' Panic: subindex not found'
            STOP
          ENDIF

          child => thisOctal%child(subIndex)
           
          CALL getIntersectedOctals(child,listHead,grid,nOctals) 
                            
        ELSE ! childless
          
          IF (PRESENT(onlyChanged)) THEN
            IF (onlyChanged .AND. thisOctal%changed(iSubcell)) CYCLE
          END IF
          
          IF (.NOT. octalAdded) THEN
            nOctals = nOctals + 1
            ALLOCATE(newElement)
            newElement%content => thisOctal
            newElement%next    => listHead
            newElement%inUse   = .FALSE.
            newElement%inUse(iSubcell) = .TRUE.
            listHead => newElement
            octalAdded = .TRUE.
          ELSE
            newElement%inUse(iSubcell) = .TRUE.
          END IF

        END IF
      END IF
    END DO
        
  END SUBROUTINE getIntersectedOctals

  !
  ! Based on getIntersectedOctals, but for all the octals in the first 
  ! octant (x>0, y>0 and z>0 domain).
  ! 
  ! Returns a linked list of all the octals which are in x>0, y>0 and z>0 
  ! domain. 
  RECURSIVE SUBROUTINE getOctalsInFirstOctant(thisOctal,listHead,grid,nOctals,onlyChanged)  
    !   
    ! if onlyChanged is True, we ignore cells without the Changed flag

    USE vector_mod
    
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT), TARGET :: thisOctal
    TYPE(octalListElement), POINTER :: listHead
    TYPE(gridType), INTENT(IN)      :: grid
    INTEGER,INTENT(INOUT)           :: nOctals   ! number of octals
    LOGICAL, INTENT(IN), OPTIONAL   :: onlyChanged
    
    LOGICAL               :: octalAdded
    TYPE(vector)     :: centrePoint
    TYPE(vector)     :: starPos
    TYPE(octal), POINTER  :: child
    INTEGER               :: i, iSubcell, subIndex
    TYPE(octalListElement), POINTER :: newElement
    
    octalAdded = .FALSE.
    starPos = grid%starPos1

    MAINLOOP: DO iSubcell = 1, thisOctal%maxChildren, 1 
    
      centrePoint = subcellCentre(thisOctal,iSubcell) - starPos
      
      IF (centrePoint%x < 0.0 .OR. centrePoint%y < 0.0 .OR. centrePoint%z < 0.0) THEN
         ! we are in the wrong octant.
         ! this makes some assumptions about the numerical accuracy of the
         !   grid, but it should be mostly OK.
         CYCLE  MAINLOOP

      ELSE
         IF (thisOctal%hasChild(iSubcell)) THEN
            ! call the recursive subroutine on its child           
            subIndex = -99
            DO i = 1, thisOctal%nChildren, 1
               IF ( thisOctal%indexChild(i) == iSubcell ) THEN
                  subIndex = i  
                  EXIT
               ENDIF
            ENDDO
            IF (subIndex == -99) THEN
               PRINT *, ' Panic: subindex not found'
               STOP
            ENDIF

            child => thisOctal%child(subIndex)
           
            CALL getOctalsInFirstOctant(child,listHead,grid,nOctals) 
                            
         ELSE ! childless
          
            IF (PRESENT(onlyChanged)) THEN
               IF (onlyChanged .AND. thisOctal%changed(iSubcell)) CYCLE
            END IF
            
            IF (.NOT. octalAdded) THEN
               nOctals = nOctals + 1
               ALLOCATE(newElement)
               newElement%content => thisOctal
               newElement%next    => listHead
               newElement%inUse   = .FALSE.
               newElement%inUse(iSubcell) = .TRUE.
               listHead => newElement
               octalAdded = .TRUE.
            ELSE
               newElement%inUse(iSubcell) = .TRUE.
            END IF
            
         END IF

      END IF

   END DO MAINLOOP
        
 END SUBROUTINE getOctalsInFirstOctant
  

  SUBROUTINE moveOctalListToArray(listHead,octalArray)
    ! copies all the pointers (to octals) from a linked list into an array.
    ! assumes that the SIZE of octalArray is the number of list elements.
  
    IMPLICIT NONE 
    
    TYPE(octalListElement), POINTER :: listHead
    TYPE(octalWrapper),DIMENSION(:) :: octalArray

    TYPE(octalListElement), POINTER :: oldHead
    INTEGER :: element
    
    DO element = 1, SIZE(octalArray), 1

      octalArray(element)%content => listHead%content
      octalArray(element)%inUse   =  listHead%inUse
      
      oldHead => listHead
      listHead => listHead%next
      DEALLOCATE(oldHead)

    END DO
    NULLIFY(listHead)
  
  END SUBROUTINE moveOctalListToArray


  SUBROUTINE amrUpdateGrid(amrLimitScalar,amrLimitScalar2,grid)
    ! checks whether each octal has changed significantly since the last phase
    !   of the simulation, and splits/deletes it as necessary.

    IMPLICIT NONE

    real(double), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
                                          !   routines that this subroutine calls
    
    TYPE(OCTAL), POINTER :: thisOctal

    IF ((.NOT. grid%adaptive) .OR. (.NOT.grid%geometry(1:6)=="ttauri")) THEN
      print *, 'amrUpdateGrid probably doesn''t work with this setup'
      stop
    END IF
  
    ! first we check whether we are going to delete any child octals
   
    !print *, 'Deleting cells no longer needed...'
    !thisOctal => grid%octreeRoot
    !CALL amrUpdateGridDelete(thisOctal)
   
    ! then we check if we are adding any new octals
    
    print *, 'Adding new cells to grid...'
    thisOctal => grid%octreeRoot
    CALL amrUpdateGridAdd(thisOctal)

    ! finally, check whether any subcells have changed significantly
    print *, 'Flagging updated cells...'
    thisOctal => grid%octreeRoot
    CALL amrUpdateGridChanged(thisOctal)
    
  CONTAINS

    RECURSIVE SUBROUTINE amrUpdateGridDelete(thisOctal)
      ! if a child should no longer exist, we delete it.

      TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal
      TYPE(octal), POINTER :: thisChild
      logical :: splitInAzimuth
      INTEGER :: iChild, iSubcell

      DO iChild = 1, thisOctal%nChildren, 1
        
        thisChild => thisOctal%child(iChild)
        
        DO iSubcell = 1, thisChild%maxChildren, 1
          IF (thisOctal%indexChild(iSubcell) == iChild) EXIT
        END DO 
        
        IF (.NOT. decideSplit(thisOctal,iSubcell,amrLimitScalar,amrLimitScalar2,grid,splitInAzimuth)) THEN
          PRINT *, 'Deleting unneeded children'
          CALL deleteChild(parent=thisOctal, childToDelete=iSubcell, &
                           adjustParent=.FALSE., grid=grid,          &
                           adjustGridInfo=.TRUE.)
          ! setting adjustParent and adjustGridInfo to be true might
          !   not be sensible here?

          thisOctal%changed(iSubcell) = .TRUE.
        ELSE
          CALL amrUpdateGridDelete(thisChild)
        END IF

      END DO  
    
    END SUBROUTINE amrUpdateGridDelete
    
    RECURSIVE SUBROUTINE amrUpdateGridAdd(thisOctal)
      ! subdivide any octals that now exceed the threshold

      TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal
      TYPE(octal), POINTER :: thisChild
      INTEGER :: iSubcell, j
      logical :: splitInAzimuth

      DO iSubcell = 1, thisOctal%maxChildren, 1
        IF (thisOctal%hasChild(iSubcell)) THEN 
          ! find the child
          DO j = 1, thisOctal%nChildren, 1
            IF (thisOctal%indexChild(j) == iSubcell) THEN
              thisChild => thisOctal%child(j)
              CALL amrUpdateGridAdd(thisChild)
              EXIT
            END IF
          END DO
        ELSE 
          IF (decideSplit(thisOctal,iSubcell,amrLimitScalar,amrLimitScalar2,grid,splitinAzimuth)) THEN
            PRINT *, 'Adding new child'
            CALL addNewChild(thisOctal,iSubcell,grid,adjustGridInfo=.TRUE.,splitAzimuthally=splitinAzimuth)
            thisOctal%child(iSubcell)%changed = .TRUE.
            grid%nOctals = grid%nOctals + 1
          END IF
        END IF 
          
      END DO 
        
    END SUBROUTINE amrUpdateGridAdd

    RECURSIVE SUBROUTINE amrUpdateGridChanged(thisOctal)
      ! update the octals in the grid
      ! flag any octals that have changed "significantly" since the last phase

      ! NB This is only for T Tauri models. It should be made more elegant and general!
      
      TYPE(octal), TARGET, INTENT(INOUT) :: thisOctal
      TYPE(octal), POINTER :: thisChild
      INTEGER :: iSubcell
      REAL :: newDensity
      INTEGER :: iChild

      DO iSubcell = 1, thisOctal%maxChildren, 1

!        newDensity = TTauriDensity(subcellCentre(thisOctal,iSubcell),grid)
         newDensity = Density(subcellCentre(thisOctal,iSubcell),grid)
          IF ( ABS((newDensity/(MAX(thisOctal%rho(iSubcell),1.d-25))-1.0)) > 0.1 ) &
            thisOctal%changed(iSubcell) = .TRUE.
         
        CALL calcTTauriMassVelocity(thisOctal,iSubcell,grid)

      END DO

      DO iChild = 1, thisOctal%nChildren, 1
        thisChild => thisOctal%child(iChild)
        CALL amrUpdateGridChanged(thisChild)
      END DO 
      
    END SUBROUTINE amrUpdateGridChanged

  END SUBROUTINE amrUpdateGrid

  subroutine returnKappa(grid, thisOctal, subcell, ilambda, lambda, kappaSca, kappaAbs, kappaAbsArray, kappaScaArray, &
       rosselandKappa, kappap, atthistemperature, kappaAbsDust, kappaAbsGas, kappaScaDust, kappaScaGas)
    use input_variables, only: nDustType, photoionization !, mie, includeGasOpacity
    use atom_mod, only: bnu
    implicit none
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer, optional :: ilambda
    real, optional :: lambda
    real(double), optional :: kappaSca, kappaAbs, kappaAbsArray(grid%nLambda), kappaScaArray(grid%nLambda)
    real(double), optional :: rosselandKappa
    real(double), optional :: kappaAbsDust, kappaScaDust, kappaAbsGas, kappaScaGas
    real, optional :: kappap
    real, optional :: atthistemperature
    real :: temperature
    real :: frac
!    real :: tlambda
    real, parameter :: sublimationTemp = 1500., subRange = 100.
!    real :: tArray(1000)
    real(double) :: freq, dfreq, norm !,  bnutot
    integer :: i,j,m
    real :: fac
    real :: e, h0, he0
    real(double) :: kappaH, kappaHe

    real(double), allocatable, save :: oneKappaAbsT(:,:)
    real(double), allocatable, save :: oneKappaScaT(:,:)

    logical,save :: firsttime = .true. 
    integer(double),save :: nlambda
    

    !! Commented out by DAR 14/10/08 as I don't think this is doing anything here. Moved to
    !! a more relevant place - rosselandkappa and others
!    temperature = thisOctal%temperature(subcell)


!    if (temperature < sublimationTemp) frac = 1.
!    if (temperature > (sublimationTemp+subRange)) frac = 0.
!    
!    if ((temperature > sublimationTemp).and.(temperature < (sublimationTemp+subRange))) then
!       frac = 1.-(temperature-sublimationTemp)/subRange
!    endif

!    frac = max(1.e-20,frac)

    frac = 1.

    if(firsttime) then

       allocate(oneKappaAbsT(grid%nlambda, ndusttype))
       allocate(oneKappaScaT(grid%nlambda, ndusttype))

       do i = 1, ndusttype
          do m = 1, grid%nLambda
             oneKappaAbsT(m, i) = grid%oneKappaAbs(i,m)
             oneKappaScaT(m, i) = grid%oneKappaSca(i,m)
          enddo
       enddo

       nlambda = grid%nlambda
       firsttime = .false.
    endif

    if (PRESENT(kappaAbsArray)) then

      if(grid%onekappa) then

         if (nDustType .eq. 1) then
            kappaAbsArray(1:nLambda) = thisOctal%rho(subcell) * oneKappaAbsT(1:nlambda,1)
         else
            kappaAbsArray(1:grid%nLambda) = 0.            
            do i = 1, nDustType
               kappaAbsArray(1:nLambda) = kappaAbsArray(1:nLambda) + thisOctal%dustTypeFraction(subcell, i) * & 
                    oneKappaAbsT(1:nLambda,i)*thisOctal%rho(subcell) * frac
            enddo
         endif
      else
         kappaAbsArray(1:nlambda) = thisOctal%kappaAbs(subcell,:)
      endif
!       if (includeGasOpacity) then
!          call returnGasKappaValue(temperature, thisOctal%rho(subcell),  kappaAbsArray=tarray)
!          kappaAbsArray(1:grid%nLambda) = kappaAbsArray(1:grid%nLambda) + tarray(1:grid%nLambda)*thisOctal%rho(subcell)
!       endif
!       write(*,*) nDustType,thisOctal%dusttypeFraction(subcell,1), grid%oneKappaAbs(1,1:grid%nLambda)

    endif
       
    if (PRESENT(kappaScaArray)) then
       if(grid%onekappa) then
 
          if (nDustType .eq. 1) then
             kappaScaArray(1:nLambda) = thisOctal%rho(subcell) * oneKappaScaT(1:nlambda,1)
          else
             kappaScaArray(1:grid%nLambda) = 0.   
             do i = 1, nDustType
                kappaScaArray(1:nLambda) = kappaScaArray(1:nLambda) + thisOctal%dustTypeFraction(subcell, i) * & 
                     oneKappaScaT(1:nLambda,i)*thisOctal%rho(subcell) * frac
             enddo
          endif
       else
          kappaScaArray(1:nlambda) = thisOctal%kappaSca(subcell,:)
       endif
!       if (includeGasOpacity) then
!          call returnGasKappaValue(temperature, thisOctal%rho(subcell),  kappaScaArray=tarray)
!          kappaScaArray(1:grid%nLambda) = kappaScaArray(1:grid%nLambda) + tarray(1:grid%nLambda)*thisOctal%rho(subcell)
!       endif
    endif


    if (PRESENT(kappaSca)) then
       kappaSca = 0.
       IF (.NOT.PRESENT(lambda)) THEN
          if (grid%oneKappa) then
             kappaSca = 0.
             
             if(ndusttype .eq. 1) then
                kappaSca = OneKappaScaT(iLambda,1) * thisOctal%rho(subcell) 
             else
                do i = 1, nDustType
                   kappaSca = kappaSca + thisOctal%dustTypeFraction(subcell, i) * grid%oneKappaSca(i,iLambda)*thisOctal%rho(subcell)
                enddo
             endif
          else 
             ! For line computation (without dust).
             ! Needs modification for a model which include gas and dust here.
             ! This is a temporarily solution
             kappaSca = thisOctal%kappaSca(subcell,iLambda)
          end if
       else
          kappaSca = 0.
          do i = 1, nDustType
             kappaSca = kappaSca + thisOctal%dustTypeFraction(subcell, i) * &
                  logint(dble(lambda), dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), &
                  grid%oneKappaSca(i,iLambda)*thisOctal%rho(subcell), &
                  grid%oneKappaSca(i,iLambda+1)*thisOctal%rho(subcell))
          enddo
       endif
       kappaSca = kappaSca * frac
    endif
    if (PRESENT(kappaScaDust)) kappaScaDust = kappaSca

    if (PRESENT(kappaAbs)) then
       kappaAbs = 0.
       
       IF (.NOT.PRESENT(lambda)) THEN
          if (grid%oneKappa) then
             kappaAbs = 0.
             if(ndustType .eq. 1) then
                kappaAbs = oneKappaAbsT(iLambda,1)*thisOctal%rho(subcell)
             else
                do i = 1, nDustType
                   kappaAbs = kappaAbs + thisOctal%dustTypeFraction(subcell, i) * grid%oneKappaAbs(i,iLambda)*thisOctal%rho(subcell)
                enddo
             endif
          else
             ! For line computation (without dust).
             ! Needs modification for a model which include gas and dust here.
             ! This is a temporarily solution
             kappaAbs = thisOctal%kappaAbs(subcell,iLambda)
          end if
       else
          kappaAbs = 0.
          if (ndusttype .eq. 1) then
             kappaAbs = logint(dble(lambda), dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), &
                  oneKappaAbsT(ilambda,1)*thisoctal%rho(subcell), &
                  oneKappaAbsT(ilambda+1,1)*thisoctal%rho(subcell))

          else
          
             do i = 1, nDustType
             !write(*,*) grid%oneKappaAbs(i,iLambda),grid%oneKappaAbs(i,iLambda+1),thisOctal%rho(subcell), &
             !     dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), ilambda, dble(lambda)
             !write(*,*) logint(dble(lambda), dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), &
             !     grid%oneKappaAbs(i,iLambda)*thisOctal%rho(subcell), &
             !     grid%oneKappaAbs(i,iLambda+1)*thisOctal%rho(subcell))
             !write(*,*) i, subcell
             !write(*,*) nDusttype,kappaAbs 
             !write(*,*) thisOctal%dustTypeFraction(subcell, i)
                kappaAbs = kappaAbs + thisOctal%dustTypeFraction(subcell, i) * &
                     logint(dble(lambda), dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), &
                     grid%oneKappaAbs(i,iLambda)*thisOctal%rho(subcell), &
                     grid%oneKappaAbs(i,iLambda+1)*thisOctal%rho(subcell))

             enddo
          endif
       endif
!       write(*,*) nDustType,thisOctal%dusttypeFraction(subcell,1), grid%oneKappaAbs(1,1:grid%nLambda)
       
       kappaAbs = kappaAbs * frac

   endif
   if (PRESENT(kappaAbsDust)) kappaAbsDust = kappaAbs


!   if (includeGasOpacity) then
!      if (.not.PRESENT(lambda)) then
!         tlambda = grid%lamArray(iLambda)
!      else
!         tlambda = lambda
!      endif
!
!      if (PRESENT(kappaAbs)) then
!         call returnGasKappaValue(temperature, thisOctal%rho(subcell), tlambda, kappaAbs=kappaAbsGas)
!         kappaAbs = kappaAbs + kappaAbsGas*thisOctal%rho(subcell)
!      endif
!      
!      if (PRESENT(kappaSca)) then
!         call returnGasKappaValue(temperature, thisOctal%rho(subcell), tlambda, kappaSca=kappaScaGas)
!         kappaSca = kappaSca + kappaScaGas*thisOctal%rho(subcell)
!      endif
!   endif

   if (PRESENT(rosselandKappa)) then
      temperature = thisOctal%temperature(subcell)

      if (PRESENT(atthistemperature)) then
         temperature = atthistemperature
      endif
      if (temperature < grid%tempRossArray(grid%nTempRossArray)) then
         call locate(grid%tempRossArray, grid%nTempRossArray, temperature, m)
         fac = (temperature - grid%tempRossArray(m))/(grid%tempRossArray(m+1)-grid%tempRossArray(m))
      else
         m = grid%nTempRossArray-1
         fac = 1.
      endif
      rosselandKappa = 0.
      
      do i = 1, nDustType
         rosselandKappa = rosselandKappa + (grid%kappaRossArray(i, m) + &
              fac*(grid%kappaRossArray(i,m+1)-grid%kappaRossArray(i,m))) * &
              max(1.d-30,thisOctal%dustTypeFraction(subcell, i))
      enddo

   endif
   

   if (PRESENT(kappap)) then
      temperature = thisOctal%temperature(subcell)

      if (PRESENT(atthistemperature)) then
         temperature = atthistemperature
      endif
      kappaP = 0.d0
      norm = 0.d0
      do i = grid%nLambda,2,-1
         freq = cSpeed / (grid%lamArray(i)*1.e-8)
         dfreq = cSpeed / (grid%lamArray(i)*1.e-8) - cSpeed / (grid%lamArray(i-1)*1.e-8)
         do j = 1, nDustType
            kappaP = kappaP + thisOctal%dustTypeFraction(subcell, j) * dble(grid%oneKappaAbs(j,i)) * &
                 thisOctal%rho(subcell) *&
                 dble(bnu(dble(freq),dble(temperature)))  * dfreq
         enddo
         norm = norm + dble(bnu(dble(freq),dble(temperature)))  * dfreq
      enddo
      if (norm /= 0.d0) then
         kappaP = (kappaP / norm) /1.d10
      else
         kappaP = 1.d-30
      endif
   endif
   

   if (photoionization) then

      if (PRESENT(kappaAbs)) then
         if (present(lambda)) then
            e = (hCgs * (cSpeed / (lambda * 1.e-8))) * ergtoev
         else
            e = (hCgs * (cSpeed / (grid%lamArray(iLambda) * 1.e-8))) * ergtoev
            !            write(*,*) "! using rough grid"
         endif
         call phfit2(1, 1, 1 , e , h0)
         call phfit2(2, 2, 1 , e , he0)
         kappaH =  thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1) * h0
         kappaHe = thisOctal%nh(subcell)*grid%ion(3)%abundance*thisOctal%ionFrac(subcell,3) * he0
         kappaAbs = kappaAbs + (kappaH + kappaHe)
      endif
      if (PRESENT(kappaAbsGas)) kappaAbsGas = (kappaH + kappaHe)
      if (PRESENT(kappaSca)) then
         kappaSca = kappaSca + thisOctal%ne(subcell) * sigmaE * 1.e10
      endif
      if (PRESENT(kappaScaGas)) kappaScaGas = thisOctal%ne(subcell) * sigmaE * 1.e10 
   endif
   
  end subroutine returnKappa
   


   subroutine averageofNearbyCells(grid, thisOctal, thisSubcell, meantemp, meanrho)
     type(GRIDTYPE) :: grid
     type(OCTAL),pointer :: thisOctal
     type(OCTAL),pointer :: nearbyOctal
     real(double) :: rho, meanRho
     integer :: thisSubcell
     integer :: subcell
     integer :: i
     real :: temp, meanTemp
     integer :: nTemp
     type(VECTOR), allocatable :: locator(:)
     real(double) :: dSubcellcentre
     type(VECTOR) :: thisSubcellCentre
     integer :: nlocator, j
     meanTemp = 0.
     meanRho = 0.d0
     nTemp = 0
     
    if (thisOctal%threed) then
       nLocator = 26
    else
       nlocator = 8
    endif

    ALLOCATE(locator(nLocator))
    thisSubcellCentre = subcellCentre(thisOctal,thisSubcell)
    FORALL (i = 1:nLocator)
       locator(i) = thisSubcellCentre
    END FORALL
     
     ! Moving this distance in any direction will take us into a
     ! neighbouring subcell.
     dSubcellCentre = (0.5 * thisOctal%subcellSize) + grid%halfSmallestSubcell
     
        if (thisOctal%threed) then 

           ! faces
           locator(1)%x = thisSubcellCentre%x + dSubcellCentre
           locator(2)%y = thisSubcellCentre%y + dSubcellCentre
           locator(3)%z = thisSubcellCentre%z + dSubcellCentre
           locator(4)%x = thisSubcellCentre%x - dSubcellCentre
           locator(5)%y = thisSubcellCentre%y - dSubcellCentre
           locator(6)%z = thisSubcellCentre%z - dSubcellCentre
           ! x-edges
           locator(7)%y = thisSubcellCentre%y + dSubcellCentre
           locator(7)%z = thisSubcellCentre%z + dSubcellCentre
           locator(8)%y = thisSubcellCentre%y + dSubcellCentre
           locator(8)%z = thisSubcellCentre%z - dSubcellCentre
           locator(9)%y = thisSubcellCentre%y - dSubcellCentre
           locator(9)%z = thisSubcellCentre%z + dSubcellCentre
           locator(10)%y = thisSubcellCentre%y - dSubcellCentre
           locator(10)%z = thisSubcellCentre%z - dSubcellCentre
           ! y-edges
           locator(11)%x = thisSubcellCentre%x + dSubcellCentre
           locator(11)%z = thisSubcellCentre%z + dSubcellCentre
           locator(12)%x = thisSubcellCentre%x + dSubcellCentre
           locator(12)%z = thisSubcellCentre%z - dSubcellCentre
           locator(13)%x = thisSubcellCentre%x - dSubcellCentre
           locator(13)%z = thisSubcellCentre%z + dSubcellCentre
           locator(14)%x = thisSubcellCentre%x - dSubcellCentre
           locator(14)%z = thisSubcellCentre%z - dSubcellCentre
           ! z-edges
           locator(15)%x = thisSubcellCentre%x + dSubcellCentre
           locator(15)%y = thisSubcellCentre%y + dSubcellCentre
           locator(16)%x = thisSubcellCentre%x + dSubcellCentre
           locator(16)%y = thisSubcellCentre%y - dSubcellCentre
           locator(17)%x = thisSubcellCentre%x - dSubcellCentre
           locator(17)%y = thisSubcellCentre%y + dSubcellCentre
           locator(18)%x = thisSubcellCentre%x - dSubcellCentre
           locator(18)%y = thisSubcellCentre%y - dSubcellCentre
           ! corners
           locator(19)%x = thisSubcellCentre%x + dSubcellCentre
           locator(19)%y = thisSubcellCentre%y + dSubcellCentre
           locator(19)%z = thisSubcellCentre%z + dSubcellCentre
           locator(20)%x = thisSubcellCentre%x + dSubcellCentre
           locator(20)%y = thisSubcellCentre%y + dSubcellCentre
           locator(20)%z = thisSubcellCentre%z - dSubcellCentre
           locator(21)%x = thisSubcellCentre%x + dSubcellCentre
           locator(21)%y = thisSubcellCentre%y - dSubcellCentre
           locator(21)%z = thisSubcellCentre%z + dSubcellCentre
           locator(22)%x = thisSubcellCentre%x + dSubcellCentre
           locator(22)%y = thisSubcellCentre%y - dSubcellCentre
           locator(22)%z = thisSubcellCentre%z - dSubcellCentre
           locator(23)%x = thisSubcellCentre%x - dSubcellCentre
           locator(23)%y = thisSubcellCentre%y + dSubcellCentre
           locator(23)%z = thisSubcellCentre%z + dSubcellCentre
           locator(24)%x = thisSubcellCentre%x - dSubcellCentre
           locator(24)%y = thisSubcellCentre%y + dSubcellCentre
           locator(24)%z = thisSubcellCentre%z - dSubcellCentre
           locator(25)%x = thisSubcellCentre%x - dSubcellCentre
           locator(25)%y = thisSubcellCentre%y - dSubcellCentre
           locator(25)%z = thisSubcellCentre%z + dSubcellCentre
           locator(26)%x = thisSubcellCentre%x - dSubcellCentre
           locator(26)%y = thisSubcellCentre%y - dSubcellCentre
           locator(26)%z = thisSubcellCentre%z - dSubcellCentre
        else

           ! edges

           locator(1)%x = thisSubcellCentre%x + dSubcellCentre
           locator(2)%z = thisSubcellCentre%z + dSubcellCentre
           locator(3)%x = thisSubcellCentre%x - dSubcellCentre
           locator(4)%z = thisSubcellCentre%z - dSubcellCentre
           ! corners
           locator(5)%x = thisSubcellCentre%x + dSubcellCentre
           locator(5)%z = thisSubcellCentre%z + dSubcellCentre

           locator(6)%x = thisSubcellCentre%x + dSubcellCentre
           locator(6)%z = thisSubcellCentre%z - dSubcellCentre

           locator(7)%x = thisSubcellCentre%x - dSubcellCentre
           locator(7)%z = thisSubcellCentre%z + dSubcellCentre

           locator(8)%x = thisSubcellCentre%x - dSubcellCentre
           locator(8)%z = thisSubcellCentre%z - dSubcellCentre


        endif


     do j = 1, nLocator
        call amrGridValues(grid%octreeRoot,locator(j),foundOctal=nearbyOctal, &
             foundsubcell=subcell, temperature=temp, rho=rho)
        if ( rho > amr_min_rho ) then
           meanTemp = meanTemp + temp
           meanRho = meanRho + log10(rho)
           ntemp = ntemp + 1
        endif
     end do
     if (nTemp > 0) then
        meanTemp = meanTemp/ real(nTemp)
        meanRho = meanRho /  dble(nTemp)
     else
        meanTemp = thisOctal%temperature(thisSubcell)
        meanRho = log10(thisOctal%rho(thisSubcell))
     endif
     meanRho = 10.d0**meanRho
     deallocate(locator)
   end subroutine averageofNearbyCells

  recursive subroutine estimateRhoOfEmpty(grid, thisOctal)
    USE sph_data_class, only: sphData, get_udist
    USE cluster_class, only:  find_n_particle_in_subcell

    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    real :: temp
    real(double) :: rho, rho_tmp
    integer :: subcell, i, np
    type(VECTOR) :: rVec
    real(double), allocatable, save :: recip_sm(:)
    real(double)  :: udist
    logical, save :: first_time = .true.

    if (first_time) then 
       allocate( recip_sm(sphData%npart) )
       
       udist = get_udist()
       do i=1, sphData%npart 

          recip_sm(i) = 1.0_db / (sphData%hn(i) * udist * 1.0e-10_db)
       end do
       recip_sm(:) = recip_sm(:) * recip_sm(:)
       first_time = .false.
    end if

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call estimateRhoOfEmpty(grid, child)
                exit
             end if
          end do
       else
          call find_n_particle_in_subcell(np, rho_tmp, &
               thisOctal, subcell)
          if (np < 1) then
                rVec = subcellCentre(thisOctal,subcell)
               call averageofNearbyCells(grid, thisOctal, subcell, temp, rho)
                thisOctal%rho(subcell) = rho
             thisOctal%temperature(subcell) = temp
          endif
       endif
    enddo


  end subroutine estimateRhoOfEmpty


  subroutine setDiscPhotosphereBias(grid, iLam)
    use input_variables, only: rgap
    type(gridtype) :: grid
    integer :: iLam
    real :: xAxis(1:100000)
    real :: tau(1:100000)
    integer :: nx, i, i1, i2, iGap
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) kappaSca, kappaAbs
    type(VECTOR) :: octVec
    kappaAbs = 0.d0; kappaSca = 0.d0
    write(*,*) "Setting disc photosphere bias..."
    nx = 0
    xAxis = 0
    call getxValuesAMR(grid%octreeRoot, nx, xAxis)
    call stripSimilarValues(xAxis,nx,real(1.d-5*grid%halfSmallestSubcell))

    tau(1) = 0.
    do i = 2, nx
       octVec = VECTOR(xAxis(i), 0.d0, 0.d0)
       CALL findSubcellTD(octVec,grid%octreeRoot,thisOctal,subcell)
       call returnKappa(grid, thisOctal, subcell, ilambda=ilam,&
            kappaSca=kappaSca, kappaAbs=kappaAbs)
       tau(i) = tau(i-1) + (xAxis(i)-xAxis(i-1)) * (kappaAbs + kappaSca)
       call setVerticalBias(grid, xAxis(i), 0., ilam, thisTau=tau(i))
       if (tau(i) > 5.) then
          i1 = i
          exit
       endif
    enddo
       

    call locate(xAxis, nx, rGap*real(autocm)/1.e10, iGap)

    tau(i+1) = 0.
    do i = iGap, 1, -1
       octVec = VECTOR(xAxis(i), 0.d0, 0.d0)
       CALL findSubcellTD(octVec,grid%octreeRoot,thisOctal,subcell)
       call returnKappa(grid, thisOctal, subcell, ilambda=ilam,&
            kappaSca=kappaSca, kappaAbs=kappaAbs)
       tau(i) = tau(i+1) + (xAxis(i+1)-xAxis(i)) * (kappaAbs + kappaSca)
       call setVerticalBias(grid, xAxis(i), 0., ilam, thisTau=tau(i))
       if (tau(i) > 5.) then
          i2 = i
          exit
       endif
    enddo
       

    write(*,*) "range",xAxis(i1),xAxis(i2)
    do i = i1, i2
       call setVerticalBias(grid, xAxis(i), 0., ilam)
    enddo

    tau(i-1) = 0.
    do i = iGap, nx
       octVec = VECTOR(xAxis(i), 0.d0, 0.d0)
       CALL findSubcellTD(octVec,grid%octreeRoot,thisOctal,subcell)
       tau(i) = tau(i-1) + (xAxis(i)-xAxis(i-1)) * (kappaAbs + kappaSca)
       call setVerticalBias(grid, xAxis(i), 0., ilam, thisTau=tau(i))
       if (tau(i) > 5.) then
          i1 = i
          exit
       endif
    enddo
       
    write(*,*) "from ",xAxis(i1)
    do i = i1 , nx
       call setVerticalBias(grid, xAxis(i), 0., ilam)
    enddo








    write(*,*) "Done."

  end subroutine setDiscPhotosphereBias

    

  recursive subroutine getxValuesAMR(thisOctal, nx, xAxis)

    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(vector) :: rVec
    integer :: nx, subcell, i
    real :: xAxis(:)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getxValuesamr(child, nx, xAxis)
                exit
             end if
          end do
       else

          rVec = subcellCentre(thisOctal, subcell)
          nx = nx + 1
          xAxis(nx) = rVec%x
       end if
    end do

  end subroutine getxValuesAMR

  subroutine setVerticalBias(grid, xPos, yPos, iLam, thisTau)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    real, optional :: thisTau
    integer :: nz
    real :: zAxis(100000)
    real :: tau
    real :: xPos, yPos
    integer :: subcell
    real(double) :: rhotemp
    real :: temptemp
    integer :: iLam
    type(VECTOR) :: currentPos, temp
    real :: halfSmallestSubcell
    real(double) :: kappaSca, kappaAbs
    kappaSca = 0.d0; kappaAbs = 0.d0
    nz = 0
    tau = 0.
    halfSmallestSubcell = grid%halfSmallestSubcell

    ! bottom up

    currentPos = VECTOR(xpos, yPos, -1.*grid%octreeRoot%subcellsize)
    do while(currentPos%z < 0.)
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp, temperature=temptemp, grid=grid)
       nz = nz + 1
       temp = subCellCentre(thisOctal, subcell)
       zAxis(nz) = temp%z

       call returnKappa(grid, thisOctal, subcell, ilambda=ilam,&
            kappaSca=kappaSca, kappaAbs=kappaAbs)


       if (nz > 1) then
          tau = tau + thisOctal%subcellsize*(kappaSca+kappaAbs)
       endif
       if (.not.present(thisTau)) then
          thisOctal%biasCont3D(subcell) = max(1.e-6,exp(-tau))
       else
          thisOctal%biasCont3D(subcell) = max(1.e-6,exp(-thistau))
       endif
             
       currentPos = VECTOR(xpos, yPos, zAxis(nz)+0.5*thisOctal%subcellsize+halfSmallestSubcell)
    end do
    
    tau = 0.
    nz = 0
    currentPos = VECTOR(xpos, yPos, grid%octreeRoot%subcellsize)
    do while(currentPos%z > 0.)
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp, temperature=temptemp, grid=grid)
       nz = nz + 1
       temp = subCellCentre(thisOctal, subcell)
       zAxis(nz) = temp%z
       call returnKappa(grid, thisOctal, subcell, ilambda=ilam,&
            kappaSca=kappaSca, kappaAbs=kappaAbs)

       if (nz > 1) then
          tau = tau + thisOctal%subcellsize*(kappaAbs+kappaSca)
       endif
       thisOctal%biasCont3D(subcell) = max(1.e-6,exp(-tau))
       currentPos = VECTOR(xpos, yPos, zAxis(nz)-0.5*thisOctal%subcellsize-halfSmallestSubcell)
    end do
    
  end subroutine setVerticalBias

  subroutine interpFromParent(centre, cellSize, grid, temperature, density, dusttypeFraction, thisEtaLine)
    type(GRIDTYPE) :: grid
    real(double) :: cellSize
    real :: temperature
    real(double) :: density, thisEtaLine
    real(double) :: dusttypeFraction(:)
    real(double), allocatable :: tdusttype(:,:)
    real(double), allocatable :: rho(:), etaline(:)
    real, allocatable :: temp(:)
    type(VECTOR) :: centre, octVec
    real(double) :: r
    integer :: j, i

    if (grid%octreeRoot%oneD) then
       write(*,*) "interp from parent not implemented for one-d case"
       stop
    endif

    r = cellSize/2.d0 + 0.1d0 * grid%halfSmallestSubcell
    if (grid%octreeRoot%twod) then
       j = 4
       allocate(tDustType(1:j,1:SIZE(dusttypeFraction)))
       allocate(rho(1:j), temp(1:j), etaline(1:j))

       octVec = centre + r * VECTOR(1.d0, 0.d0, 0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(1), rho=rho(1), dusttypeFraction=tdusttype(1,:), &
            etaline=etaLine(1))
       
       octVec = centre + r * VECTOR(-1.d0,0.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(2), rho=rho(2), dusttypeFraction=tdusttype(2,:), &
            etaline=etaline(2))
       
       octVec = centre + r * VECTOR(0.d0,0.d0,-1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(3), rho=rho(3), dusttypeFraction=tdusttype(3,:), &
            etaline=etaLine(3))
       
       octVec = centre + r * VECTOR(0.d0,0.d0,+1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(4), rho=rho(4), dusttypeFraction=tdusttype(4,:), &
            etaline=etaline(4))
    else
       j = 6
       allocate(tDustType(1:j,1:SIZE(dusttypeFraction)))
       allocate(rho(1:j), temp(1:j),etaline(1:j))
       octVec = centre + r * VECTOR(1.d0, 0.d0, 0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(1), rho=rho(1), dusttypeFraction=tdusttype(1,:), &
            etaline=etaline(1))
       
       octVec = centre + r * VECTOR(-1.d0,0.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(2), rho=rho(2), dusttypeFraction=tdusttype(2,:), &
            etaline=etaline(2))
       
       octVec = centre + r * VECTOR(0.d0,0.d0,-1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(3), rho=rho(3), dusttypeFraction=tdusttype(3,:), &
            etaline=etaline(3))
       
       octVec = centre + r * VECTOR(0.d0,0.d0,+1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(4), rho=rho(4), dusttypeFraction=tdusttype(4,:), &
            etaline=etaline(4))

       octVec = centre + r * VECTOR(0.d0,+1.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(5), rho=rho(5), dusttypeFraction=tdusttype(5,:), &
            etaline=etaline(5))

       octVec = centre + r * VECTOR(0.d0,-1.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(6), rho=rho(6), dusttypeFraction=tdusttype(6,:), &
            etaline=etaline(6))

    endif
  

    density = sum(rho) / dble(j)

    thisetaline = sum(etaline) / dble(j)

    temperature = sum(temp) / real(j)

    do i = 1, SIZE(dustTypefraction)
       dusttypeFraction(i) = SUM(tdusttype(1:j,i)) / dble(j)
    enddo
  end subroutine interpFromParent

  subroutine interpHydroProperties(grid, thisOctal, subcell)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, parentOctal, neighbourOctal
    type(OCTAL),pointer :: neighbourOctalMinus, neighbourOctalPlus
    integer :: subcell, parentSubcell, neighbourSubcell
    type(VECTOR) :: direction, cellCentre, rVec, parentSubcellCentre
    real(double) :: xMinus, xMid, xPlus
    integer :: neighbourSubcellPlus, neighbourSubcellMinus

    parentOctal => thisOctal%parent
    parentSubcell = thisOctal%parentSubcell
    cellCentre = subcellCentre(thisOctal, subcell)
    xMid = cellCentre%x
    parentSubcellCentre = subcellCentre(parentOctal, parentSubcell)
    direction = VECTOR(1.d0, 0.d0, 0.d0)

    neighbourOctalMinus => thisOctal
    rVec = parentSubcellCentre - direction*(parentoctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
    call findSubcellLocal(rVec, neighbourOctalMinus, neighbourSubcellMinus)
    rVec = subcellCentre(neighbourOctalminus, neighbourSubcellMinus)
    xMinus = rvec%x
    neighbourOctalPlus => thisOctal
    rVec = parentSubcellCentre - direction*(parentoctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
    call findSubcellLocal(rVec, neighbourOctalPlus, neighbourSubcellPlus)
    rVec = subcellCentre(neighbourOctalPlus, neighbourSubcellPlus)
    xPlus = rvec%x
    

  end subroutine interpHydroProperties



  recursive subroutine myTauSmooth(thisOctal, grid, ilambda, converged, inheritProps, interpProps, photosphereSplit)
    use input_variables, only : tauSmoothMin, tauSmoothMax, erOuter, router, maxDepthAmr, rinner
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    logical, optional :: inheritProps, interpProps, photoSphereSplit
    integer :: subcell, i, ilambda
    logical :: converged
    real(double) :: kabs, ksca, r, fac
    type(VECTOR) :: dirVec(6), centre, octVec, aHat, rVec
    real :: thisTau, neighbourTau
    integer :: neighbourSubcell, j, nDir
    logical :: split
    logical, save :: firsttime = .true.
    character(len=30) :: message
    logical :: converged_tmp
    kabs = 0.d0; ksca = 0.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call myTauSmooth(child, grid, ilambda, converged, inheritProps, interpProps, photosphereSplit)
                exit
             end if
          end do
       else

          split = .true.

          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda,&
               kappaSca=ksca, kappaAbs=kabs)

          thisTau  = thisOctal%subcellSize * (ksca + kabs)

          r = thisOctal%subcellSize/2.d0 + grid%halfSmallestSubcell * 0.001d0
          centre = subcellCentre(thisOctal, subcell)
          if (thisOctal%threed) then
             if (.not.thisOctal%cylindrical) then
                nDir = 6
                dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
                dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
                dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
                dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
                dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
                dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
             else
                nDir = 4
                aHat = VECTOR(centre%x, centre%y, 0.d0)
                call normalize(aHat)
                dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
                dirVec(2) = aHat
                dirVec(3) = (-1.d0)*aHat
                dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
             endif
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)  
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          endif

          do j = 1, nDir
             octVec = centre + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then
                startOctal => thisOctal
                call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                     foundOctal=neighbourOctal, foundsubcell=neighbourSubcell, kappaSca=ksca, kappaAbs=kabs, ilambda=ilambda)
                neighbourTau = neighbourOctal%subcellSize * (ksca + kabs)

                if ((grid%geometry.eq."whitney").and.&
                     (modulus(subcellCentre(thisOctal,subcell)) > 0.9*erouter/1.e10)) split = .false.

                rVec = subcellCentre(thisOctal,subcell)

                if ((grid%geometry.eq."warpeddisc").and.&
                     (sqrt(rVec%x**2+rVec%y**2)) > 0.9*router) split = .false.


                if ((grid%geometry.eq."shakara").and.&
                     ((rVec%x+thisOctal%subcellsize/2.d0) < rinner)) split = .false.



                if (thisOctal%nDepth == maxDepthamr) then
                   split = .false.
                   if (firstTime) then
                      write(message,'(a,i3)') "AMR cell depth capped at: ",maxDepthamr
                      call writeWarning(message)
                      firstTime = .false.
                   endif
                endif
                   
                
                if ((min(thisTau, neighbourTau) < tauSmoothMin).and.(max(thisTau, neighbourTau) > tauSmoothMax).and.split) then
                   if (thisTau > neighbourTau) then
                      call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                           inherit=inheritProps, interp=interpProps)
                      converged = .false.
                      return
                   endif
                endif

                fac = abs(thisOctal%temperature(subcell)-neighbourOctal%temperature(neighbourSubcell)) / &
                     thisOctal%temperature(subcell)

                if (split.and.(thisTau > neighbourTau).and.(fac > 0.05d0).and. &
                     (thisOctal%temperature(subcell)>100.d0).and.(thisTau > 1.d-4)) then
                   call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                        inherit=inheritProps, interp=interpProps)
                   converged = .false.
                   return
                endif

                if (PRESENT(photosphereSplit)) then
                   if (photosphereSplit) then
                      if ((thisOctal%etaLine(subcell) /= 0.d0).and.(neighbourOctal%etaLine(neighbourSubcell)/=0.d0)) then
                         if ((j==3).or.(j==4)) then
                            fac = abs(neighbourOctal%etaLine(neighbourSubcell)-thisOctal%etaLine(subcell))
                            if (split.and.(fac > 0.2d0).and.(thisOctal%etaLine(subcell) > 0.1d0).and. & 
                                 (thisOctal%etaLine(subcell) < 1.d0).and. &
                                 (thisOctal%etaLine(subcell)>neighbourOctal%etaLine(neighbourSubcell))) then
                               !                               write(*,*) " tau split ", fac, " eta ",thisOctal%etaline(subcell), "depth ",thisOctal%ndepth
                               call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                                    inherit=inheritProps, interp=interpProps)
                               converged = .false.
                               return
                            endif
                         endif
                      endif
                   endif
                endif



             endif
          enddo
       endif
    end do

  end subroutine myTauSmooth


  recursive subroutine myTauSplit(thisOctal, grid, converged, inheritProps, interpProps)
    use input_variables, only : tauSmoothMin, tauSmoothMax, erOuter, router, maxDepthAmr!, rinner
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    logical, optional :: inheritProps, interpProps
    integer :: subcell, i
    logical :: converged
    real(double) :: kabs, ksca, r
    type(VECTOR) :: dirVec(6), centre, octVec, aHat, rVec
    real :: thisTau
    integer :: neighbourSubcell, j, nDir
    logical :: split
    logical, save :: firsttime = .true.
    character(len=30) :: message
    kabs = 0.d0; ksca = 0.d0

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call myTauSplit(child, grid, converged, inheritProps, interpProps)
                exit
             end if
          end do
       else

          split = .true.

          call returnKappa(grid, thisOctal, subcell, rosselandKappa=kabs)
          thisTau = thisOctal%subcellsize * kabs * thisOctal%rho(subcell) * 1.e10

          if (thisOctal%nDepth == maxDepthamr) then
             split = .false.
             if (firstTime) then
                write(message,'(a,i3)') "AMR cell depth capped at: ",maxDepthamr
                call writeWarning(message)
                firstTime = .false.
             endif
          endif



          if (split.and.(thisTau > 2.).and.(thisOctal%chiLine(subcell) < 100.)) then
             write(*,*) "splitting cell with ",thisTau, " at depth ",thisOctal%chiline(subcell)
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=inheritProps, interp=interpProps)
             converged = .false.
             return
          endif
       endif
    enddo

  end subroutine myTauSplit


  recursive subroutine myScaleSmooth(factor, grid, converged, &
       inheritProps, interpProps, stellar_cluster, romData)
    type(gridtype) :: grid
    real :: factor
    integer :: nTagged
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    logical, optional :: inheritProps, interpProps
    !
    TYPE(cluster), optional, intent(in)  :: stellar_cluster
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry    
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
    real(double) :: r
    type(VECTOR) :: dirVec(6), centre, octVec, aHat
    integer :: neighbourSubcell, j, nDir

    call zeroChiLineLocal(grid%octreeRoot)
    nTagged = 0
    call tagScaleSmooth(nTagged, factor, grid%octreeRoot, grid,  converged, &
         inheritProps, interpProps, stellar_cluster, romData)
    call splitTagged(grid%octreeRoot, grid, inheritProps, interpProps, stellar_cluster, romData)

  end subroutine myScaleSmooth

  recursive subroutine zeroChiLineLocal(thisOctal)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call zeroChiLineLocal(child)
                exit
             end if
          end do
       else

          call allocateAttribute(thisOctal%chiLine,thisOctal%maxChildren)
          thisOctal%chiLine = 0.d0

       endif
    enddo
  end subroutine zeroChiLineLocal

  recursive subroutine splitTagged(thisOctal, grid, inheritProps, interpProps, stellar_cluster, romData)
    type(GRIDTYPE) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i
    logical, optional :: inheritProps, interpProps
    !
    TYPE(cluster), optional, intent(in)  :: stellar_cluster
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry    
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call splitTagged(child, grid, inheritProps, interpProps, stellar_cluster, romData)
                exit
             end if
          end do
       else

          if (thisOctal%chiline(subcell) > 0.d0) then
             thisOctal%chiLine(subcell) = 0.d0
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=inheritProps, interp=interpProps)
             return
          endif
       endif
    enddo
  end subroutine splitTagged


  recursive subroutine tagScaleSmooth(ntagged, factor, thisOctal, grid,  converged, &
       inheritProps, interpProps, stellar_cluster, romData)
    type(gridtype) :: grid
    integer :: nTagged
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    logical, optional :: inheritProps, interpProps
    !
    TYPE(cluster), optional, intent(in)  :: stellar_cluster
    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry    
    !
    integer :: subcell, i
    logical :: converged, converged_tmp
    real(double) :: r
    type(VECTOR) :: dirVec(6), centre, octVec, aHat
    integer :: neighbourSubcell, j, nDir



    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call tagScaleSmooth(nTagged, factor, child, grid, converged, inheritProps, interpProps, &
                     stellar_cluster=stellar_cluster, romData=romData)
                exit
             end if
          end do
       else

          r = thisOctal%subcellSize/2.d0 + 0.0001d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell) 
          if (thisOctal%threed) then
             if (.not.thisOctal%cylindrical) then
                nDir = 6
                dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
                dirVec(2) = VECTOR( 0.d0,+1.d0,  0.d0)
                dirVec(3) = VECTOR(+1.d0, 0.d0,  0.d0)
                dirVec(4) = VECTOR(-1.d0, 0.d0,  0.d0)
                dirVec(5) = VECTOR( 0.d0,-1.d0,  0.d0)
                dirVec(6) = VECTOR( 0.d0, 0.d0, -1.d0)
             else
                nDir = 4
                aHat = VECTOR(centre%x, centre%y, 0.d0)
                call normalize(aHat)
                dirVec(1) = VECTOR( 0.d0, 0.d0, +1.d0)
                dirVec(2) = aHat
                dirVec(3) = (-1.d0)*aHat
                dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
             endif
          else if (thisOctal%twod) then
             nDir = 4
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)  
             dirVec(2) = VECTOR(-1.d0,0.d0, 0.d0)
             dirVec(3) = VECTOR( 0.d0, 0.d0,  1.d0)
             dirVec(4) = VECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 2
             dirVec(1) = VECTOR( 1.d0, 0.d0, 0.d0)
             dirVec(2) = VECTOR(-1.d0, 0.d0, 0.d0)
          endif
          do j = 1, nDir
             octVec = centre + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then

                startOctal => thisOctal
                call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                     foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)

                if ((thisOctal%subcellSize/neighbourOctal%subcellSize) > factor) then
                   thisOctal%chiLine(subcell) = 1.d0
                   converged = .false.
                   nTagged = nTagged + 1
                endif

                if ((neighbourOctal%subcellSize/thisOctal%subcellSize) > factor) then
                   neighbourOctal%chiLine(neighboursubcell) = 1.d0
                   converged = .false.
                   nTagged = nTagged + 1
                endif

             endif
          enddo
       endif
    end do

  end subroutine tagScaleSmooth

  subroutine distanceToCellBoundary(grid, posVec, direction, tVal, sOctal, sSubcell)

    implicit none
    type(GRIDTYPE), intent(in)    :: grid
    type(VECTOR), intent(in) :: posVec
    type(VECTOR), intent(in) :: direction
    type(OCTAL), pointer, optional :: sOctal
    integer, optional :: sSubcell
    real(oct), intent(out) :: tval
    !
    type(VECTOR), parameter :: norm(6) = (/ VECTOR( 1.0d0,  0.d0,   0.0d0),  &
                                            VECTOR( 0.0d0,  1.0d0,  0.0d0),  &
                                            VECTOR( 0.0d0,  0.0d0,  1.0d0),  &
                                            VECTOR(-1.0d0,  0.0d0,  0.0d0),  &
                                            VECTOR( 0.0d0, -1.0d0,  0.0d0),  &
                                            VECTOR( 0.0d0,  0.0d0, -1.0d0) /)
    type(VECTOR) :: rDirection
    type(OCTAL),pointer :: thisOctal
    real(double) :: distTor1, distTor2, theta, mu
    real(double) :: distToRboundary, compz,currentZ
    real(double) :: phi, distToZboundary, ang1, ang2
    type(VECTOR) :: subcen, point, xHat, rVec, rplane, rnorm, xVec
    integer :: subcell
    real(double) :: distToSide1, distToSide2, distToSide
    real(double) ::  compx,disttoxBoundary, halfCellSize, d2, fac
    real(oct) :: t(6),denom(6), r, r1, r2, d, cosmu,x1,x2, halfSubCellsize
    real(double) :: a, b, c 
    logical :: ok, thisOk(6)
!    integer :: jarray(6)

    type(VECTOR) :: normdiff

    point = posVec

    if (PRESENT(sOctal)) then
       if (PRESENT(sSubcell)) then
          subcell = sSubcell
          thisOctal => sOctal
       else
          call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid, startOctal=sOctal)
       endif
    else
       call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
    endif
    subcen =  subcellCentre(thisOctal,subcell)

    if (thisOctal%oneD) then

       distToR1 = 1.d30
       distToR2 = 1.d30

       rVec = posVec
       call normalize(rVec)
       cosmu = ((-1.d0)*direction).dot.rVec
       d = modulus(posVec)

       ! distance to outer radius

       r2 = subcen%x + thisOctal%subcellSize/2.d0
       call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
       distToR2 = max(x1,x2)
!             write(*,*) "r2",x1,x2,disttor2

       !   inner radius

       r1 = subcen%x - thisOctal%subcellSize/2.d0
       theta = asin(max(-1.d0,min(1.d0,r1 / d)))
       cosmu =((-1.d0)*rVec).dot.direction
       mu = acos(max(-1.d0,min(1.d0,cosmu)))
       distTor1 = 1.e30
       if (mu  < theta ) then
          call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
          distTor1 = min(x1,x2)
       endif
!             write(*,*) "r1",x1,x2,disttor1,mu,theta

       tval = min(distTor1, distTor2)
       goto 666
    endif

    if (thisOctal%threed) then

       if (.not.thisOctal%cylindrical) then
          ok = .true.

          halfSubCellsize = thisOctal%subcellsize * 0.5d0

          if(direction%x .ne. 0.d0) then            
             denom(1) = 1.d0 / direction%x
          else
             denom(1) = 0.d0
          endif
          denom(4) = -denom(1)

          if(direction%y .ne. 0.d0) then            
             denom(2) = 1.d0 / direction%y
          else
             denom(2) = 0.d0
          endif
          denom(5) = -denom(2)
          
          if(direction%z .ne. 0.d0) then            
             denom(3) = 1.d0 / direction%z
          else
             denom(3) = 0.d0
          endif
          denom(6) = -denom(3)

          normdiff = subcen - posvec
          
          thisOK = .false.
          ok = .false.
          t(1) =  (normdiff%x + halfsubcellsize) * denom(1)
          if (t(1) > 0.d0) then
             thisOK(1) = .true.
             ok = .true.
          endif
          t(2) =  (normdiff%y + halfsubcellsize) * denom(2)
          if (t(2) > 0.d0) then
             thisOK(2) = .true.
             ok = .true.
          endif
          t(3) =  (normdiff%z + halfsubcellsize) * denom(3)
          if (t(3) > 0.d0) then
             thisOK(3) = .true.
             ok = .true.
          endif
          t(4) =  (normdiff%x - halfsubcellsize) * denom(1)
          if (t(4) > 0.d0) then
             thisOK(4) = .true.
             ok = .true.
          endif
          t(5) =  (normdiff%y - halfsubcellsize) * denom(2)
          if (t(5) > 0.d0) then
             thisOK(5) = .true.
             ok = .true.
          endif
          t(6) =  (normdiff%z - halfsubcellsize) * denom(3)
          if (t(6) > 0.d0) then
             thisOK(6) = .true.
             ok = .true.
          endif

!          where(t > 0.d0)
!             jarray = 1
!             thisOk = .true.
!          elsewhere
!             jarray = 0
!             thisOk = .false.
!          end where
!          
!          j = sum(jarray)

!          if (j .eq. 0) ok = .false.
 
          if (.not.ok) then
             write(*,*) "Error: j=0 (no intersection???) in amr_mod::distanceToCellBoundary. "
             write(*,*) direction%x,direction%y,direction%z
             write(*,*) t(1:6)
             call torus_abort
          endif
          
          tval = minval(t, mask=thisOk)
          
! Commented out by Dave Acreman, October 2008
! tval == 0 is handled at the end of this subroutine    
!          if (tval == 0.) then
!             write(*,*) posVec
!             write(*,*) direction%x,direction%y,direction%z
!             write(*,*) t(1:6)
!             call torus_abort("tval==0 in distanceToCellBoundary")
!          endif

          !if (tval > sqrt(3.)*thisOctal%subcellsize) then
             !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
             !     write(*,*) "direction",direction
             !     write(*,*) t(1:6)
             !     write(*,*) denom(1:6)
          !endif

       else

          ! now look at the cylindrical case

          halfCellSize = thisOctal%subcellSize/2.d0
          rVec = subcellCentre(thisOctal,subcell)
          r = sqrt(rVec%x**2 + rVec%y**2)
          r1 = r - halfCellSize
          r2 = r + halfCellSize

          distToR1 = 1.d30
          distToR2 = 1.d30
          d = sqrt(point%x**2+point%y**2)
          xHat = VECTOR(point%x, point%y,0.d0)
          if (modulus(xhat)/=0.d0)    call normalize(xHat)
          rDirection = VECTOR(direction%x, direction%y,0.d0)
          compX = modulus(rDirection)
          if (modulus(rDirection) /= 0.d0) call normalize(rDirection)
          if (compX /= 0.d0) then
             cosmu =((-1.d0)*xHat).dot.rdirection
             call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
             if (.not.ok) then
                write(*,*) "Quad solver failed in intersectcubeamr2d I",d,cosmu,r2
                write(*,*) "xhat",xhat
                write(*,*) "dir",direction
                write(*,*) "point",point
                do;enddo
                endif
                if ((x1.lt.0.d0).and.(x2.lt.0.d0)) then
                   write(*,*) "x1, x2 ",x1,x2
                   write(*,*) "rdirection ",rdirection
                   write(*,*) "xhat ",xhat
                   write(*,*) "compx ",compx
                   write(*,*) "cosmu ",cosmu
                endif
                distTor2 = max(x1,x2)/compX

                if ((d .ne. 0.).and.(r1 > 0.1d0*grid%halfSmallestSubcell)) then
                   theta = asin(max(-1.d0,min(1.d0,r1 / d)))
                   cosmu = ((-1.d0)*xHat).dot.rdirection
                   mu = acos(max(-1.d0,min(1.d0,cosmu)))
                   distTor1 = 1.e30
                   if (mu  < theta ) then
                      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
                      if (.not.ok) then
                         write(*,*) "Quad solver failed in intersectcubeamr2d II",d,cosmu,r1,x1,x2
                         write(*,*) "coeff b",-2.d0*d*cosmu, "coeff c", d**2-r2**2
                         write(*,*) "direction ",direction
                         write(*,*) "mu ",mu*radtodeg, "theta ",theta*radtodeg
                         x1 = thisoctal%subcellSize/2.d0
                         x2 = 0.d0
                      endif
                      distTor1 = min(x1,x2)/compX
                   endif
                else
                   distTor1 = 1.d30
                end if
             endif
             distToRboundary = min(distTor1, distTor2)

             ! now do the upper and lower (z axis) surfaces

             compZ = zHat.dot.direction
             currentZ = point%z

             if (compZ /= 0.d0 ) then
                if (compZ > 0.d0) then
                   distToZboundary = (subcen%z + halfCellSize - currentZ ) / compZ
                else
                   distToZboundary = abs((subcen%z - halfCellSize - currentZ ) / compZ)
                endif
             else
                disttoZboundary = 1.e30
             endif

             ! ok now we have to tackle the two angled sides...
         thisOk = .true.
         
! Commented out by Dave Acreman, October 2008
! tval == 0 is handled at the end of this subroutine         
!         if (tval == 0.) then
!            write(*,*) posVec
!            write(*,*) direction%x,direction%y,direction%z
!            write(*,*) t(1:6)
!            call torus_abort("tval==0 in distanceToCellBoundary")
!         endif
         
         !if (tval > sqrt(3.)*thisOctal%subcellsize) then
            !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
            !     write(*,*) "direction",direction
            !     write(*,*) t(1:6)
            !     write(*,*) denom(1:6)
         !endif

             rVec = subcellCentre(thisOctal, subcell)
             phi = atan2(rVec%y, rVec%x)
             if (phi < 0.d0) phi = phi + twoPi

             ang1 = phi - returndPhi(thisOctal)
             rPlane = VECTOR(cos(ang1),sin(ang1),0.d0)
             rnorm = rplane .cross. VECTOR(0.d0, 0.d0, 1.d0)
             call normalize(rnorm)
             distToSide1 = 1.d30
             if ((rnorm .dot. direction) /= 0.d0) then
                distToSide1 = (rnorm.dot.(rPlane-posVec))/(rnorm.dot.direction)
                if (distToSide1 < 0.d0) distToSide1 = 1.d30
             endif

             ang2 = phi + returndPhi(thisOctal)
             rPlane = VECTOR(cos(ang2),sin(ang2),0.d0)
             rnorm = rplane .cross.  VECTOR(0.d0, 0.d0, 1.d0)
             call normalize(rnorm)
             distToSide2 = 1.d30
             if ((rnorm .dot. direction) /= 0.d0) then
                distToSide2 = (rnorm.dot.(rPlane-posVec))/(rnorm.dot.direction)
                if (distToSide2 < 0.d0) distToSide2 = 1.d30
             endif

             distToSide = min(distToSide1, distToSide2)


             tVal = min(distToZboundary, distToRboundary, distToSide)
             if (tVal > 1.e29) then
                write(*,*) "Cylindrical ",tval
                write(*,*) tVal,compX,compZ, distToZboundary,disttorboundary, disttoside
                write(*,*) "subcen",subcen
                write(*,*) "z", currentZ
             endif
             if (tval < 0.) then
                write(*,*) "Cylindrical ",tval
                write(*,*) tVal,distToZboundary,disttorboundary, disttoside
                write(*,*) "subcen",subcen
                write(*,*) "z", currentZ
                write(*,*) "disttor1, disttor2 ",disttor1,disttor2
             endif

          endif

       else ! two-d grid case below

          halfCellSize = thisOctal%subcellSize/2.d0
          r1 = subcen%x - halfCellSize
          r2 = subcen%x + halfCellSize

          distToR1 = 1.d30
          distToR2 = 1.d30
          d2 = point%x**2+point%y**2
          d = sqrt(d2)
          xVec = VECTOR(point%x, point%y,0.d0)
          xHat = xVec
          call normalize(xHat)
          rDirection = VECTOR(direction%x, direction%y,0.d0)
          compX = modulus(rDirection)
          call normalize(rDirection)

          if (compX /= 0.d0) then
             cosmu =((-1.d0)*xHat).dot.rdirection
             call solveQuadDble(1.d0, -2.d0*d*cosmu, d2-r2**2, x1, x2, ok)
             if (.not.ok) then
                write(*,*) "Quad solver failed in intersectcubeamr2d I",d,cosmu,r2
                write(*,*) "xhat",xhat
                write(*,*) "dir",direction
                write(*,*) "point",point
                do
                enddo
                x1 = thisoctal%subcellSize/2.d0
                x2 = 0.d0
             endif
             distTor2 = max(x1,x2)/compX

             distTor1 = 1.e30

             if (cosmu > 0.d0) then
                a = 1.d0
                b = -2.d0*d*cosmu
                c = d2-r1**2
                fac = b*b-4.d0*a*c
                if (fac > 0.d0) then
                   call solveQuadDble(a, b, c, x1, x2, ok)
                   !               if(ok) then
                   !                  write(*,*) "All good",d,cosmu,r1,x1,x2
                   !                  write(*,*) "coeff b",-2.d0*d*cosmu, "coeff c", d**2-r2**2
                   !               endif
                   
                   if (.not.ok) then
                      write(*,*) "mu", mu, "theta", theta
                      write(*,*) "Quad solver failed in intersectcubeamr2d IIb",d,cosmu,r1,x1,x2
                      write(*,*) "coeff b",-2.d0*d*cosmu, "coeff c", d**2-r1**2
                      write(*,*) "xhat",xhat
                      write(*,*) "dir",direction
                      write(*,*) "point",point
                      
                      x1 = thisoctal%subcellSize/2.d0
                      x2 = 0.d0
                   endif
                   distTor1 = min(x1,x2)/compX
                endif
             endif
          endif
          distToXboundary = min(distTor1, distTor2)


          compZ = zHat.dot.direction
          currentZ = point%z

          if (compZ /= 0.d0 ) then
             if (compZ > 0.d0) then
                distToZboundary = (subcen%z + halfCellSize - currentZ ) / compZ
             else
                distToZboundary = abs((subcen%z - halfCellSize - currentZ ) / compZ)
             endif
          else
             disttoZboundary = 1.e30
          endif

          tVal = min(distToZboundary, distToXboundary)
          if (tVal > 1.e29) then
             write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
             write(*,*) "subcen",subcen
             write(*,*) "z", currentZ
             write(*,*) "TVAL", tval
             write(*,*) "direction", direction
             call torus_abort
          endif

!          if (tval < 0.) then
!             write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
!             write(*,*) "subcen",subcen
             !         write(*,*) "x,z",currentX,currentZ
!          endif
      
   endif

666    continue

!       tVal = max(tVal, 0.001d0*grid%halfSmallestSubcell) ! avoid sticking on a cell boundary


     end subroutine distanceToCellBoundary

  subroutine distanceToGridEdge(grid, posVec, direction, tVal)

   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(VECTOR), intent(in) :: posVec
   type(VECTOR), intent(in) :: direction
   real(oct), intent(out) :: tval
   !
   real(double) :: distTor1, distTor2, theta, mu
   real(double) :: distToRboundary, compz,currentZ
   real(double) :: distToZboundary !, ang1, ang2 , phi
   type(VECTOR) :: subcen, point, xHat, zHat !, rVec
   real(double) :: distToSide  !, distToSide1, distToSide2
   real(double) :: disttoxBoundary, subcellsize, halfsubcellsize
   real(oct) :: t(6),denom(6), r, r1, r2, d, cosmu,x1,x2
   integer :: i,j
   logical :: ok, thisOk(6)

   type(VECTOR) :: normdiff

   point = posVec

   subcen =  grid%octreeRoot%centre
   subcellsize = grid%octreeRoot%subcellSize * 2.d0
   halfsubcellsize = grid%octreeRoot%subcellSize

   if (grid%octreeRoot%threed) then

      if (.not.grid%octreeRoot%cylindrical) then
         ok = .true.
         
         if(direction%x .ne. 0.d0) then            
            denom(1) = 1.d0 / direction%x
         else
            denom(1) = 0.d0
         endif
         
         if(direction%y .ne. 0.d0) then            
            denom(2) = 1.d0 / direction%y
         else
            denom(2) = 0.d0
         endif
         
         if(direction%z .ne. 0.d0) then            
            denom(3) = 1.d0 / direction%z
         else
            denom(3) = 0.d0
         endif
                     
         normdiff = subcen - posvec

         t(1) =  (normdiff%x + halfsubcellsize) * denom(1)
         t(2) =  (normdiff%y + halfsubcellsize) * denom(2)
         t(3) =  (normdiff%z + halfsubcellsize) * denom(3)
         t(4) =  (normdiff%x - halfsubcellsize) * denom(1)
         t(5) =  (normdiff%y - halfsubcellsize) * denom(2)
         t(6) =  (normdiff%z - halfsubcellsize) * denom(3)

         thisOk = .true.
         
         do i = 1, 6
            
            if (denom(i) .eq. 0.0d0) then
               thisOk(i) = .false.
            endif

            if (t(i) < 0.) thisOk(i) = .false.
         enddo
                  
         j = 0
         do i = 1, 6
            if (thisOk(i)) j=j+1
         enddo
         
         if (j == 0) ok = .false.
         
         if (.not.ok) then
            write(*,*) "Error: j=0 (no intersection???) in amr_mod::distanceToGridEdge. "
            write(*,*) direction%x,direction%y,direction%z
            write(*,*) t(1:6)
            stop
         endif
         
         tval = minval(t, mask=thisOk)

         if (tval == 0.) then
            write(*,*) posVec
            write(*,*) direction%x,direction%y,direction%z
            write(*,*) t(1:6)
            stop
         endif
         
         if (tval > sqrt(3.)*subcellsize) then
            !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
            !     write(*,*) "direction",direction
            !     write(*,*) t(1:6)
            !     write(*,*) denom(1:6)
         endif

      else

! now look at the cylindrical case

         ! first do the inside and outside curved surfaces
         r = sqrt(subcen%x**2 + subcen%y**2)
         r1 = r -subcellSize/2.d0
         r2 = r +subcellSize/2.d0
         d = sqrt(point%x**2+point%y**2)
         xHat = VECTOR(point%x, point%y,0.d0)
         call normalize(xHat)
      
         cosmu =((-1.d0)*xHat).dot.direction
         call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
         if (.not.ok) then
            write(*,*) "Quad solver failed in intersectcubeamr2d"
            x1 = subcellSize/2.d0
            x2 = 0.d0
         endif
         distTor2 = max(x1,x2)
         
         theta = asin(max(-1.d0,min(1.d0,r1 / d)))
         cosmu = xHat.dot.direction
         mu = acos(max(-1.d0,min(1.d0,cosmu)))
         distTor1 = 1.e30
         if (mu  < theta ) then
            call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
            if (.not.ok) then
               write(*,*) "Quad solver failed in intersectcubeamr2d"
               x1 = subcellSize/2.d0
               x2 = 0.d0
            endif
            distTor1 = max(x1,x2)
         endif
      
         distToRboundary = min(distTor1, distTor2)

         ! now do the upper and lower (z axis) surfaces
      
         zHat = VECTOR(0.d0, 0.d0, 1.d0)
         compZ = zHat.dot.direction
         currentZ = point%z
      
         if (compZ /= 0.d0 ) then
            if (compZ > 0.d0) then
               distToZboundary = (subcen%z + subcellsize/2.d0 - currentZ ) / compZ
            else
               distToZboundary = abs((subcen%z - subcellsize/2.d0 - currentZ ) / compZ)
            endif
         else
            disttoZboundary = 1.e30
         endif
      
        
!         ! ok now we have to tackle the two angled sides...
!
!         ! find posvec to surface centre
!
!         phi = atan2(posVec%y,posVec%x)
!         if (phi < 0.d0) phi = phi + twoPi
!
!         rVec = VECTOR(r, 0.d0, 0.d0)
!         if (grid%octreeRoot%splitAzimuthally) then
!            if (phi < grid%octreeRoot%phi) then
!               ang1 = grid%octreeRoot%phi - grid%octreeRoot%dPhi/2.d0
!               ang2 = grid%octreeRoot%phi
!            else
!               ang1 = grid%octreeRoot%phi
!               ang2 = grid%octreeRoot%phi + grid%octreeRoot%dPhi/2.d0
!            endif
!         else
!            ang1 = grid%octreeRoot%phi - grid%octreeRoot%dPhi/2.d0
!            ang2 = grid%octreeRoot%phi + grid%octreeRoot%dPhi/2.d0
!         endif
!
!         rVec = VECTOR(r, 0.d0, 0.d0)
!         rVec = rotateZ(rVec, -ang1)
!         thisnorm = rVec .cross. zHat
!         call normalize(thisnorm)
!         if ((thisnorm.dot.direction) /= 0.d0) then
!            distToSide1 = (thisnorm.dot.(rVec-posVec))/(thisnorm.dot.direction)
!            if (distToSide1 < 0.d0) distToSide1 = 1.d30
!         endif
!
!         rVec = VECTOR(r, 0.d0, 0.d0)
!         rVec = rotateZ(rVec, -ang2)
!         thisnorm = rVec .cross. zHat
!         call normalize(thisnorm)
!         if ((thisnorm.dot.direction) /= 0.d0) then
!            distToSide2 = (thisnorm.dot.(rVec-posVec))/(thisnorm.dot.direction)
!            if (distToSide2 < 0.d0) distToSide2 = 1.d30
!         endif
!
!         distToSide = min(distToSide1, distToside2)
         distToSide = 1.d30

         tVal = min(distToZboundary, distToRboundary, distToSide)

         write(*,*) disttoside,disttoZboundary,disttoRboundary
         if (tVal > 1.e29) then
            write(*,*) "Cylindrical"
            write(*,*) tVal,compZ, distToZboundary,disttorboundary, disttoside
            write(*,*) "subcen",subcen
            write(*,*) "z", currentZ
         endif
         if (tval < 0.) then
            write(*,*) "Cylindrical"
            write(*,*) tVal,distToZboundary,disttorboundary, disttoside
            write(*,*) "subcen",subcen
            write(*,*) "z", currentZ
         endif

      endif

   else ! two-d grid case below

      r1 = subcen%x - subcellSize/2.d0
      r2 = subcen%x + subcellSize/2.d0
      d = sqrt(point%x**2+point%y**2)
      xHat = VECTOR(point%x, point%y,0.d0)
      call normalize(xHat)
      
      cosmu =((-1.d0)*xHat).dot.direction
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in intersectcubeamr2d"
         x1 = subcellSize/2.d0
         x2 = 0.d0
      endif
      distTor2 = max(x1,x2)
      
      theta = asin(max(-1.d0,min(1.d0,r1 / d)))
      cosmu = xHat.dot.direction
      mu = acos(max(-1.d0,min(1.d0,cosmu)))
      distTor1 = 1.e30
      if (mu  < theta ) then
         call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
         if (.not.ok) then
            write(*,*) "Quad solver failed in intersectcubeamr2d"
            x1 = subcellSize/2.d0
            x2 = 0.d0
         endif
         distTor1 = max(x1,x2)
      endif
      
      distToXboundary = min(distTor1, distTor2)
      
      
      zHat = VECTOR(0.d0, 0.d0, 1.d0)
      compZ = zHat.dot.direction
      currentZ = point%z
      
      if (compZ /= 0.d0 ) then
         if (compZ > 0.d0) then
            distToZboundary = (subcen%z + subcellsize/2.d0 - currentZ ) / compZ
         else
            distToZboundary = abs((subcen%z - subcellsize/2.d0 - currentZ ) / compZ)
         endif
      else
         disttoZboundary = 1.e30
      endif
      
      tVal = min(distToZboundary, distToXboundary)
      if (tVal > 1.e29) then
         write(*,*) tVal,compZ, distToZboundary,disttoxboundary
         write(*,*) "subcen",subcen
         write(*,*) "z",currentZ
      endif
      if (tval < 0.) then
         write(*,*) tVal,compZ, distToZboundary,disttoxboundary
         write(*,*) "subcen",subcen
         write(*,*) "z", currentZ
      endif
      
   endif

   tVal = max(tVal, 1.d-4*grid%halfSmallestSubcell) ! avoid sticking on a cell boundary

 end subroutine distanceToGridEdge


  type(VECTOR) function randomPositionInCell(thisOctal, subcell)
    type(OCTAL) :: thisOctal
    integer :: subcell
    type(VECTOR) :: octalCentre
    real(double) :: r1, r2, r3, r
    real(double) :: xOctal, yOctal, zOctal
    real(double) :: ang, ang1, ang2, phi


    octalCentre = subcellCentre(thisOctal,subcell)
    
!!! we will just choose a random point within the subcell.
!!! this *should* be done in a better way.


    if (thisOctal%oneD) then
       call random_number(r1)
       r1 = r1 - 0.5d0
       r1 = r1 * 0.9999
       r = r1 * thisOctal%subcellSize + octalCentre%x
       randomPositionInCell = r * randomUnitVector()
       goto 666
    endif


    if (thisOctal%threed) then

       if (.not.thisOctal%cylindrical) then

          call random_number(r1)
          r1 = r1 - 0.5  ! shift value mean value to zero
          r1 = r1 * 0.9999 ! to avoid any numerical accuracy problems
          xOctal = r1 * thisOctal%subcellSize + octalCentre%x
       
          call random_number(r2)
          r2 = r2 - 0.5                                  
          r2 = r2 * 0.9999                                          
          yOctal = r2 * thisOctal%subcellSize + octalCentre%y
          
          
          call random_number(r3)
          r3 = r3 - 0.5                                  
          r3 = r3 * 0.9999                                          
          zOctal = r3 * thisOctal%subcellSize + octalCentre%z
          
          randomPositionInCell = VECTOR(xOctal,yOctal,zOctal)
          
       else


          call random_number(r1)
          r1 = r1 - 0.5  ! shift value mean value to zero
          r1 = r1 * 0.9999 ! to avoid any numerical accuracy problems
          xOctal = r1 * thisOctal%subcellSize + sqrt(octalCentre%x**2+octalCentre%y**2)


          call random_number(r3)
          r3 = r3 - 0.5                                  
          r3 = r3 * 0.9999                                          
          zOctal = r3 * thisOctal%subcellSize + octalCentre%z

          randomPositionInCell = VECTOR(xOctal,0.,zOctal)
          
          call random_number(r2)
          phi = atan2(octalCentre%y, octalCentre%x)
          if (phi < 0.d0) phi = phi + twoPi

          ang1 = phi - returndPhi(thisOctal)
          ang2 = phi + returndPhi(thisOctal)
          
!          if (thisOctal%splitAzimuthally) then
!             ang1 = thisOctal%phi - thisOctal%dphi/4.
!             ang2 = thisOctal%phi + thisOctal%dphi/4.
!          else
!             ang1 = thisOctal%phi - thisOctal%dphi/2.
!             ang2 = thisOctal%phi + thisOctal%dphi/2.
!          endif
          ang = ang1 + r2 * (ang2 - ang1)
          randomPositionInCell = rotateZ(randomPositionInCell, -ang)

       endif
    else

       call random_number(r1)
       r1 = r1 - 0.5  ! shift value mean value to zero
       r1 = r1 * 0.9999 ! to avoid any numerical accuracy problems
       xOctal = r1 * thisOctal%subcellSize + octalCentre%x
       
          
       call random_number(r3)
       r3 = r3 - 0.5                                  
       r3 = r3 * 0.9999                                          
       zOctal = r3 * thisOctal%subcellSize + octalCentre%z
       
       randomPositionInCell = VECTOR(xOctal,0.d0,zOctal)

       if (thisOctal%twod) then
          call random_number(ang)
          ang = ang * twoPi
          randomPositionInCell = rotateZ(randomPositionInCell, ang)
       endif
    endif
666 continue
  end function randomPositionInCell


  !-----------------------------------------------------------------------------------
  ! Base on startReturnSample of NHS. Optimized for a solving formal solution.
  !    --- (R. Kurosawa)
  ! Returns the amr grid values along a ray in the direction (direction) originated
  ! from "startPoint".
  ! Main differences from startReturnSample is:
  !   1. If a ray encontour "thin disc", starting point will be shifted to the point
  !      of insertersection. 
  !------------------------------------------------------------------------------------
  SUBROUTINE amr_values_along_ray (startPoint,direction,grid,          &
             sampleFreq,nSamples,maxSamples,thin_disc_on, opaqueCore,hitCore, fromDisc, &
             usePops,iLambda,error,lambda,kappaAbs,kappaSca,velocity,&
             velocityDeriv,chiLine,levelPop,rho, temperature, Ne, inflow, &
             etaCont, etaLine)
    ! samples the grid at points along the path.
    ! this should be called by the program, instead of calling 
    !   returnSamples directly, because this checks the start and finish
    !   points are within the grid bounds. returnSamples *assumes* this 
    !   criterion is met. also, the path of the photon is checked for 
    !   intersection with the stellar surface(s), disc etc. 
 
    IMPLICIT NONE

    TYPE(vector), INTENT(INOUT)   :: startPoint ! photon start point
    TYPE(vector), INTENT(IN)      :: direction  ! photon direction 
    TYPE(gridtype), INTENT(IN)         :: grid       ! the entire grid structure
    REAL, INTENT(IN)                   :: sampleFreq ! the maximum number of
!    real(oct), INTENT(IN)   :: sampleFreq ! the maximum number of 
                       ! samples that will be made in any subcell of the octree. 
                       
    INTEGER, INTENT(OUT)               :: nSamples   ! number of samples made
    INTEGER, INTENT(IN)                :: maxSamples ! size of sample arrays 
    logical, intent(in)                :: thin_disc_on   ! T to include thin disc
    LOGICAL, INTENT(IN)                :: opaqueCore ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: hitCore    ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: fromDisc   ! starting point the ray was adjusted to the disc
    LOGICAL, INTENT(IN)                :: usePops    ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda    ! wavelength index
    INTEGER, INTENT(INOUT)             :: error      ! error code
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda     ! path distances of sample locations

    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs   ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca   ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: chiLine    ! line opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: rho        ! density at sample points
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL:: levelPop ! level populations
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: temperature! temperature [K] at sample points
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: Ne ! electron density
    logical,DIMENSION(:),INTENT(INOUT),OPTIONAL  ::inFlow   ! flag to tell if the point is inflow
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity

    TYPE(vector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    TYPE(vector)       :: currentPoint ! current position of ray 
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(octal)             :: octree       ! the octree structure within 'grid'
    TYPE(vector)       :: directionNormalized
    real(oct)    :: distanceLimit ! max length of ray before aborting
    ! margin is the size of the region around the edge of a subcell
    !   where numerical inaccuracies may cause problems.
    real(oct)    :: margin
 
    ! variables for testing special cases (stellar intersections etc.)
    TYPE(vector)       :: starPosition       ! position vector of stellar centre
    TYPE(vector)       :: diskNormal         ! disk normal vector
    real(oct)    :: rStar              ! stellar radius
    real(oct)    :: endLength          ! max path length of photon
    TYPE(vector)       :: endPoint           ! where photon leaves grid
    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
    real(oct)    :: distanceFromOrigin ! closest distance of plane from origin
    TYPE(vector)       :: diskIntersection   ! point of photon intersection with disk
    real(oct)    :: starIntersectionDistance1 ! distance to first  intersection
    real(oct)    :: starIntersectionDistance2 ! distance to second intersection
    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
    real(oct)    :: intersectionRadius ! disk radius when photon intersects
    real(oct)    :: diskDistance       ! distance to disk intersection
    real(oct)    :: distanceThroughStar! distance of chord through star
    TYPE(vector)       :: dummyStartPoint    ! modified start point 
!    real(oct), PARAMETER :: fudgefactor = 1.00001 ! overestimates stellar size
    real(oct), PARAMETER :: fudgefactor = 1.000001 ! overestimates stellar size
    
    
    ! we will abort tracking a photon just before it reaches the edge of the
    !   simulation space. This is the fraction of the total distance to use:
    real(oct), PARAMETER :: distanceFraction = 0.999_oc 

    ! used internally
    TYPE(vector)       ::startPointNew       ! 
!    TYPE(vector)       ::entryPoint       ! 
 
    
    ! we will abort tracking any rays which are too close to 
    !   to a cell wall. The distance to use is defined by:
    margin = 6.0_oc * REAL(grid%maxDepth,kind=oct) * EPSILON(1.0_oc)
    ! some more experimentation is required to find out the best value for
    !   margin.
    
    
    ! set up some variables
    octree = grid%octreeRoot  
    abortRay = .FALSE.
    hitCore = .FALSE.
    fromDisc = .false.
    directionNormalized = direction
    CALL normalize(directionNormalized)
    distanceLimit = HUGE(distanceLimit)
    locator = startPoint
    nSamples = 0
    startPointnew = startPoint

    currentPoint = startPointNew
    
    IF (.NOT. inOctal(octree,startpointnew)) THEN
       
      PRINT *, 'Attempting to find path between point(s) outwith the grid.'
      PRINT *, ' in [amr_mod::amr_values_along_ray].'
      PRINT *, ' ==> StartPoint = (', startPoint, ')'
      error = -30

      !
      print *, "Set the size of the box should be more than the distance between "
      print *, "the center of the star and the edge of the density field."

      stop

    ENDIF
   
   
    ! geometry-specific tests should go here
    IF (grid%geometry(1:6) == "ttauri" .OR. grid%geometry(1:4) == "jets" .or. &
         grid%geometry(1:9) == "luc_cir3d" .or. grid%geometry(1:6) == "cmfgen" .or. &
         grid%geometry(1:8) == "romanova") THEN
      
       ! need to test for both star and disc intersections 
      
       ! we will find out when and where the photon leaves the simulation space 
       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 

       endLength = endLength * distanceFraction
       
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startpointnew
       
       ! for the moment, we will just assume a 2-d disc.
       absorbPhoton = .FALSE.
       diskNormal = grid%diskNormal
       starPosition = grid%starPos1
                                  
       ! find the (geometrycally thin) disk intersection point
       if (grid%geometry == "luc_cir3d" .or. grid%geometry == "cmfgen" ) then
          intersectionFound = .false.
       else
          if (thin_disc_on) then
             distanceFromOrigin = modulus(grid%starPos1)
             diskIntersection = intersectionLinePlane(startpointnew, directionNormalized,&
                  diskNormal, distanceFromOrigin, intersectionFound)
          else
             intersectionFound = .false.
          end if
       end if

       
       IF (intersectionFound) THEN
       
         ! we need to check whether the photon passes through the disk's
         !   central hole, or the outside of the outer radius of accretion disc.
         intersectionRadius =  modulus(diskIntersection - starPosition)
         IF (intersectionRadius > grid%diskRadius .and.  &
              intersectionRadius < 1.5d5) THEN  ! assuming 100 AU disc size

           ! we need to check whether the intersection occurs within our
           !   simulation space.
           diskDistance = modulus(diskIntersection-startpointnew)
           IF (diskDistance < endLength) THEN
              ! Ajdust the current poistion to the insersecting point
              abortRay = .FALSE.
              fromDisc = .true.
              currentPoint = diskIntersection + (fudgefactor*direction)
              startpointnew = currentPoint
              endLength = fudgefactor*(endLength - diskDistance)
              locator = currentPoint                  

           else
              ! no need to do any calculation for this ray.
              abortRay = .TRUE.

           END IF
         END IF
       END IF
         
         
       ! now we check for intersections with the star

       rStar = grid%rStar1 
       CALL intersectionLineSphere(startpointnew,directionNormalized,endLength,starPosition, &
                                   rStar,starIntersectionFound1,starIntersectionFound2,   &
                                   starIntersectionDistance1,starIntersectionDistance2)
       ! by passing a line segment to intersectionLineSphere, we ensure that we
       !   do not find intersections with the star that take place after the 
       !   photon has been absorbed by the disk.
       IF (starIntersectionFound1) THEN
       
         endLength = starIntersectionDistance1
         IF (opaqueCore) absorbPhoton = .TRUE.
         
       END IF

       ! we trace the photon path until it encounters the disk, the star,
       !   or the edge of the simulation space.
       error = 0
       CALL returnSamples(currentPoint,startpointnew,locator,directionNormalized,&
                     octree,grid,sampleFreq,nSamples,maxSamples,abortRay,     &
                     lambda,usePops,iLambda,error,margin,endLength,           &
                     kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,   &
                     velocityDeriv=velocityDeriv,chiLine=chiLine,             &
                     levelPop=levelPop,rho=rho, temperature=temperature,      &
                     Ne=Ne, inFlow=inFlow, etaCont=etaCont, etaLine=etaLine)
      
       IF (error < 0)  RETURN

       IF (.NOT. opaqueCore .AND. starIntersectionFound2) THEN 
         ! the photon passes through the star and continues on the other side.

         ! we have to adjust the arguments to returnSamples to make the 
         !   output arrays correct.
         distanceThroughStar = fudgeFactor * &
                      (starIntersectionDistance2 - starIntersectionDistance1)
         currentPoint = currentPoint + distanceThroughStar * directionNormalized
         locator = currentPoint
         dummystartpoint = startpointnew + distanceThroughStar * directionNormalized
         distanceLimit = endLength - distanceThroughStar
           
        CALL returnSamples(currentPoint,dummyStartPoint,locator,             &
                    directionNormalized,octree,grid,sampleFreq,nSamples,     &
                    maxSamples,abortRay,lambda,usePops,iLambda,error,margin, &
                    distanceLimit,kappaAbs=kappaAbs,kappaSca=kappaSca,       &
                    velocity=velocity,velocityDeriv=velocityDeriv,           &
                    chiLine=chiLine,levelPop=levelPop,rho=rho,               &
                    temperature=temperature, Ne=Ne, inFlow=inFlow,           &
                    etaCont=etaCont, etaLine=etaLine)

       END IF

       ! if the photon ends up in the star, we make sure it absorbed.
       IF (absorbPhoton) THEN
       
         nSamples = nSamples + 1
         IF (nSamples > maxSamples) THEN
           PRINT *, "Error:: nSamples > maxSamples in [amr_mod::amr_values_along_ray] "
           PRINT *, "        nSamples   = ", nSamples
           PRINT *, "        maxSamples = ", maxSamples
           STOP
         END IF
         
         lambda(nSamples) = lambda(nSamples-1)
         hitCore = .TRUE.
         IF (PRESENT(kappaSca))      kappaSca(nSamples) = 0.
         IF (PRESENT(kappaAbs))      kappaAbs(1:nSamples) = 1.e20
         IF (PRESENT(velocity))      velocity(nSamples) = vector(0.,0.,0.)
         IF (PRESENT(velocityDeriv)) velocityDeriv(nSamples) = 1.0
         IF (PRESENT(chiLine))       chiLine(1:nSamples) = 1.e20
         IF (PRESENT(levelPop))      levelPop(nSamples,:) = 0.0
         IF (PRESENT(Ne))            Ne(nSamples) = 0.0
         IF (PRESENT(rho))           rho(nSamples) = 0.0
         IF (PRESENT(temperature))   temperature(nSamples) = 0.0
         IF (PRESENT(inFlow))        inFlow(nSamples) = .true.
         IF (PRESENT(etaCont))       etaCont(nSamples) = 0.0
         IF (PRESENT(etaLine))       etaLine(nSamples) = 0.0
       END IF
       
    ELSE
            
       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 
       distanceLimit = endLength * distanceFraction
       
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startpointnew
       
       CALL returnSamples(currentPoint,startpointnew,locator,directionNormalized,&
                   octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,&
                   usePops,iLambda,error,margin,distanceLimit,                &
                   kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,     &
                   velocityDeriv=velocityDeriv,chiLine=chiLine,               &
                   levelPop=levelPop,rho=rho,temperature=temperature,         &
                   Ne=Ne,inFlow=inFlow,etaCont=etaCont,etaLine=etaLine)
    END IF

    ! The new starting point will be returned to a parent routine.
    startpoint = startpointnew
      
  end SUBROUTINE amr_values_along_ray


  !
  ! What is the difference between this and the original startReturnSamples routine??????
  ! 
  SUBROUTINE startReturnSamples2(startPoint,direction,grid,          &
             sampleFreq,nSamples,maxSamples,thin_disc_on, opaqueCore,hitCore,      &
             usePops,iLambda,error,lambda,nSource,source,kappaAbs,kappaSca,velocity,&
             velocityDeriv,chiLine,levelPop,rho, temperature, Ne, inflow, etaCont, etaLine, startOctal, startSubcell)
    ! samples the grid at points along the path.
    ! this should be called by the program, instead of calling 
    !   returnSamples directly, because this checks the start and finish
    !   points are within the grid bounds. returnSamples *assumes* this 
    !   criterion is met. also, the path of the photon is checked for 
    !   intersection with the stellar surface(s), disc etc. 
 
    USE magField, only: innerDiskData_Radii, innerDiskData_Phi, maxSizeMagFieldGrid
    USE source_mod, only: SOURCETYPE, distanceToSource

    IMPLICIT NONE

    type(octal), pointer, optional :: startOctal
    integer, optional :: startSubcell
    TYPE(vector), INTENT(IN)      :: startPoint ! photon start point
    TYPE(vector), INTENT(IN)      :: direction  ! photon direction 
    TYPE(gridtype), INTENT(IN)         :: grid       ! the entire grid structure
    REAL, INTENT(IN)                   :: sampleFreq ! the maximum number of
!    real(oct), INTENT(IN)   :: sampleFreq ! the maximum number of 
                       ! samples that will be made in any subcell of the octree. 
                       
    INTEGER, INTENT(OUT)             :: nSamples   ! number of samples made
    INTEGER, INTENT(IN)                :: maxSamples ! size of sample arrays 
    logical, intent(in)                :: thin_disc_on   ! T to include thin disc
    LOGICAL, INTENT(IN)                :: opaqueCore ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: hitCore    ! true if the core is opaque
    LOGICAL, INTENT(IN)                :: usePops    ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda    ! wavelength index
    INTEGER, INTENT(INOUT)             :: error      ! error code
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda     ! path distances of sample locations
    integer, intent(in), optional      :: nSource
    type(SOURCETYPE),intent(in), optional :: source(:)

    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs   ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca   ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: chiLine    ! line opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: rho        ! density at sample points
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL:: levelPop ! level populations
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: temperature! temperature [K] at sample points
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: Ne ! electron density
    logical,DIMENSION(:),INTENT(INOUT),OPTIONAL  ::inFlow   ! flag to tell if the point is inflow
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity

    TYPE(vector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(vector)       :: directionNormalized
 
    ! variables for testing special cases (stellar intersections etc.)
    TYPE(vector)       :: starPosition       ! position vector of stellar centre
    TYPE(vector)       :: diskNormal         ! disk normal vector
    real(oct)    :: rStar              ! stellar radius
    real(oct)    :: endLength          ! max path length of photon
    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
    real(oct)    :: distanceFromOrigin ! closest distance of plane from origin
    TYPE(vector)       :: diskIntersection   ! point of photon intersection with disk
    real(oct)    :: starIntersectionDistance1 ! distance to first  intersection
    real(oct)    :: starIntersectionDistance2 ! distance to second intersection
    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
    real(oct)    :: intersectionRadius ! disk radius when photon intersects
    real(oct)    :: diskDistance       ! distance to disk intersection
    real(oct), PARAMETER :: fudgefactor = 1.00001 ! overestimates stellar size
    
    type(vector) :: currentPosition
    type(octal), pointer :: sOctal, thisOctal
    integer :: subcell
    real(double) :: length, distToNextCell
    ! we will abort tracking a photon just before it reaches the edge of the
    !   simulation space. This is the fraction of the total distance to use:
    real(oct), PARAMETER :: distanceFraction = 0.999_oc 
    real(double) :: fudgeFac = 0.1d0
    integer :: nSource_local
    TYPE(vector)       :: vectorToIntersectionXYZ ! intersection with disk (cartesian)
    TYPE(spVector)    :: vectorToIntersectionSP ! intersection with disk (sph. polar)
    REAL(KIND=oct)    :: intersectionPhi ! disk azimuth angle where photon intersects
    INTEGER                 :: lowerBound ! index returned from 'hunt' with lower bound of phi bin
    INTEGER                 :: upperBound ! usually lowerBound+1
    REAL(KIND=oct)    :: localDiskRadius ! disk inner radius at current azimuth (1.e10cm)
    REAL(KIND=oct)    :: lowerDistance ! angular distance to lower phi bin
    REAL(KIND=oct)    :: upperDistance ! angular distance to upper phi bin
    LOGICAL           :: boundaryProblem = .false.



    currentPosition = startPoint
    locator = startPoint
    nSamples = 0
    length = 0.d0
    endLength = 1.d30
    directionNormalized = direction
    CALL normalize(directionNormalized)

    absorbPhoton = .false.
    hitcore = .false.

    ! Since not all gemoetry uses a "source" object, you need this 
    ! statement (RK).
    if (PRESENT(nSource)) then
       nSource_local = nSource
    else
       nSource_local = 0
    end if

    if (nSource_local == 0) then
       ! first we check for intersections with the star
       starPosition = grid%starPos1
       rStar = grid%rStar1 
       CALL intersectionLineSphere(startPoint,direction,1.d30,starPosition, &
            rStar,starIntersectionFound1,starIntersectionFound2,   &
            starIntersectionDistance1,starIntersectionDistance2)
       ! by passing a line segment to intersectionLineSphere, we ensure that we
       !   do not find intersections with the star that take place after the 
       !   photon has been absorbed by the disk.
       IF (starIntersectionFound1) THEN
          endLength = starIntersectionDistance1
          IF (opaqueCore) absorbPhoton = .TRUE.
       END IF
       IF (.NOT. opaqueCore .AND. starIntersectionFound2) then
          absorbPhoton = .false.
          endLength = 1.d30
       endif
    else
       call distanceToSource(source, nSource_local, startPoint, direction, absorbPhoton, endLength)
    endif


    !
    ! geometry-specific tests should go here
    !
    IF (grid%geometry(1:6) == "ttauri" .OR. grid%geometry(1:4) == "jets" .or. &
         grid%geometry(1:8) == "romanova") THEN
      
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startPoint
       
       ! for the moment, we will just assume a 2-d disc.
       absorbPhoton = .FALSE.
       diskNormal = grid%diskNormal
       starPosition = grid%starPos1
                                  
       ! find the (geometrycally thin) disk intersection point
       if (thin_disc_on) then
          distanceFromOrigin = modulus(grid%starPos1)
          diskIntersection = intersectionLinePlane(startPoint, directionNormalized,&
               diskNormal, distanceFromOrigin, intersectionFound)
       else
          intersectionFound = .false.
       end if
       
       IF (intersectionFound) THEN
       
         ! we need to check whether the photon passes through the disk's
         !   central hole, or the outside of the outer radius of accretion disc.
         intersectionRadius =  modulus(diskIntersection - starPosition)

         IF (.NOT. grid%geometry == "magstream") THEN
           localDiskRadius = grid%diskRadius
         ELSE ! grid%geometry == "magstream"
           vectorToIntersectionXYZ = diskIntersection - starPosition
           ! convert to spherical polars
            call XYZtoSPvector(vectorToIntersectionSP, vectorToIntersectionXYZ)
           intersectionPhi = vectorToIntersectionSP%phi
           
           ! need to find the appropriate disk truncation radius at this azimuth
           
           IF (.NOT. ALLOCATED(innerDiskData_Phi) .OR. &
               .NOT. ALLOCATED(innerDiskData_Radii)) THEN
             PRINT *, "Panic in startReturnSamples: innerDiskData not ALLOCATED"  
             STOP
           END IF

           CALL locate(innerDiskData_Phi,SIZE(innerDiskData_Phi),intersectionPhi,lowerBound)
           
           ! need to check for special cases - close to 0 / 2pi
           IF (lowerBound == 0) THEN
             upperBound = 1
             lowerBound = SIZE(innerDiskData_Phi)
             ! find out whether phi is closer to the lower or upper value
             lowerDistance = intersectionPhi - (innerDiskData_Phi(lowerBound) - twoPi)
             upperDistance = innerDiskData_Phi(upperBound) - intersectionPhi
           ELSE IF (lowerBound == SIZE(innerDiskData_Phi)) THEN
             upperBound = 1
             lowerBound = SIZE(innerDiskData_Phi)
             lowerDistance = intersectionPhi - innerDiskData_Phi(lowerBound)
             upperDistance = (innerDiskData_Phi(upperBound) + twoPi) - intersectionPhi
           ELSE
             ! not a special case
             upperBound = lowerBound + 1
             lowerDistance = intersectionPhi - innerDiskData_Phi(lowerBound)
             upperDistance = innerDiskData_Phi(upperBound) - intersectionPhi
           END IF

           ! if there is an angular distance of more than 5 degrees to the 
           !   nearest sample point, we will assume we are in a region
           !   with no accretion streams, and adopt the maximum disk radius
           IF ( MIN(lowerDistance,upperDistance) > (5.0_oct * degToRad) ) THEN
             localDiskRadius = maxSizeMagFieldGrid
           ! else, we have a close sample
           ELSE IF (lowerDistance <= upperDistance) THEN
             localDiskRadius = innerDiskData_Radii(lowerBound)
           ELSE
             localDiskRadius = innerDiskData_Radii(upperBound)
           END IF
           
         END IF





         IF (intersectionRadius > localDiskRadius) THEN  ! assuming 100 AU disc size
!              intersectionRadius < grid%octreeRoot%subcellsize*2.0) THEN
           absorbPhoton = .TRUE.

           ! we need to check whether the intersection occurs within our
           !   simulation space.
           diskDistance = modulus(diskIntersection-startPoint)
           ! now compare this distance with the endLength computed
           ! earlier. 
           endLength = min(endLength, diskDistance)
         END IF
       END IF
       
    end IF


    if (present(startOctal)) then
       thisOctal => startOctal
       subcell = startSubcell
    else
       CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
    endif

    do while (inOctal(grid%octreeRoot, currentPosition).and.(length < endLength))

       call findSubcellLocal(currentPosition,thisOctal,subcell,boundaryProblem)
       if (boundaryProblem) then
         boundaryProblem = .false.
         call findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
       end if

       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)

       if (nSamples < maxSamples) then
          call takeSample(currentPosition,length,direction,grid,thisOctal,subcell,nSamples,&
               maxSamples,usePops,iLambda,error,lambda,&
               kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity, velocityDeriv=velocityDeriv, &
               chiLine=chiLine,levelPop=levelPop,rho=rho,  &
               temperature=temperature,Ne=Ne,inFlow=inFlow, etaCont=etaCont, etaLine=etaLine) 
       else
          call writeWarning("Reached maxSamples limit in ray trace")
          exit
       endif
  
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       length = length + distToNextCell+fudgeFac*grid%halfSmallestSubcell
    end do


    ! if the photon ends up in the disk or the star, we make sure it absorbed.
    IF (absorbPhoton) THEN
       
       nSamples = nSamples + 1
       IF (nSamples > maxSamples) THEN
          PRINT *, "Error:: nSamples > maxSamples in startReturnSamples subroutine"
          PRINT *, "        nSamples   = ", nSamples
          PRINT *, "        maxSamples = ", maxSamples
          STOP
       END IF
       
       lambda(nSamples) = endLength
       hitCore = .TRUE.
       IF (PRESENT(kappaSca))      kappaSca(nSamples) = 0.
       IF (PRESENT(kappaAbs))      kappaAbs(nSamples) = 1.e20
       IF (PRESENT(velocity))      velocity(nSamples) = vector(0.,0.,0.)
       IF (PRESENT(velocityDeriv)) velocityDeriv(nSamples) = 1.0
       IF (PRESENT(chiLine))       chiLine(nSamples) = 1.e20
       IF (PRESENT(levelPop))      levelPop(nSamples,:) = 0.0
       IF (PRESENT(Ne))            Ne(nSamples) = 0.0
       IF (PRESENT(rho))           rho(nSamples) = 0.0
       IF (PRESENT(temperature))   temperature(nSamples) = 0.0
       IF (PRESENT(inFlow))        inFlow(nSamples) = .true.
       IF (PRESENT(etaLine))       etaLine(nSamples) = 0.0d0
       IF (PRESENT(etaCont))       etaCont(nSamples) = 0.0d0

    END IF



  end SUBROUTINE startReturnSamples2

  recursive subroutine getxValues(thisOctal, nx, xAxis)

    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: rVec
    integer :: nx, subcell, i
    real(double) :: xAxis(:)

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call getxValues(child, nx, xAxis)
                exit
             end if
          end do
       else

          rVec = subcellCentre(thisOctal, subcell)
          nx = nx + 1
          xAxis(nx) = rVec%x
       end if
    end do

  end subroutine getxValues

  recursive subroutine splitGridFractal(thisOctal, rho, aFac, grid, converged)
!    use input_variables, only : limitScalar
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    real :: rho, aFac
    integer :: subcell, i, j
!    real, allocatable :: r(:), s(:)
!    real :: rmin, rmax, tot, fac, mean
    logical :: converged

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
          if (thisOctal%nDepth < 5) then
!          if (thisOctal%rho(subcell)*cellVolume(thisOctal, subcell) > limitScalar) then
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                  inherit=.true., interp=.false.)
             ! find the child
             do j = 1, thisOctal%nChildren, 1
                if (thisOctal%indexChild(j) == subcell) then
                   child => thisOctal%child(j)
                   exit
                endif
             enddo

!             allocate(r(1:thisOctal%maxChildren), s(1:thisOctal%maxChildren))
!             call random_number(r)
!             tot = sum(r)
!             mean = tot / real(thisOctal%maxChildren)
!             r = r / mean
!             rmin = minval(r)
!             rmax = maxval(r)
!             fac = (afac*rmin-rmax)/(1.-afac)
!             s = r + fac
!             tot=SUM(s)
!             s = s / tot
!             do j = 1, thisOctal%maxChildren
!                child%rho(j) = s(j) * thisOctal%rho(subcell) * &
!                     cellVolume(thisOctal, subcell)/cellVolume(child,j)
!             enddo
!             deallocate(r, s)
             converged = .false.
             exit
          endif
          
          
       endif
    end do
  end subroutine splitGridFractal

  !
  ! Recursively turn off the magnetosphere
  !
  RECURSIVE SUBROUTINE turn_off_magnetosphere(thisOctal,grid, Rmax)    
    
    IMPLICIT NONE
    
    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    ! The outer radius of the magnetosphere.
    real(double), intent(in) :: Rmax  ! [10^10cm] 
    
    TYPE(octal), POINTER   :: pChild
    INTEGER :: subcell, n

    if (thisOctal%threeD) then
       n = 8
    else
       n =4
    end if
    

    do subcell = 1, n
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_off_magnetosphere(pChild,grid,Rmax)
       else
          ! turnning it off 
!          if (TTauriInFlow(thisOctal%centre, grid)) then          
          if ( (modulus(thisOctal%centre) <= Rmax .and. abs(thisOctal%centre%z)>5.0d0)  &
               .or.  &
               (modulus(thisOctal%centre) <= Rmax*1.01 .and. abs(thisOctal%centre%z)<5.0d0) ) then 
             thisOctal%inFlow(subcell) = .false.
             thisOctal%kappaAbs(subcell,:) = 1.0e-30
             thisOctal%kappaSca(subcell,:) = 1.0e-30          
             thisOctal%temperature(subcell) = 6500.0
             thisOctal%biasCont3D(subcell) = 1.0e-30  
             thisOctal%biasLine3D(subcell) = 1.0e-30  
             thisOctal%etaLine(subcell) = 1.e-30
             thisOctal%etaCont(subcell) = 1.e-30
             thisOctal%cornerVelocity = vector(1.e-30,1.e-30,1.e-30)
             thisOctal%velocity = vector(1.e-30,1.e-30,1.e-30)
             thisOctal%rho(subcell) = 1.0d-19
          end if
       end if
    end do
    
  END SUBROUTINE turn_off_magnetosphere


  function distanceToGridFromOutside(grid, posVec, direction) result (tval)
    use input_variables, only : suppressWarnings
    type(GRIDTYPE) :: grid
    type(VECTOR) :: subcen, direction, posVec, point, hitVec, rdirection, xhat
    type(OCTAL), pointer :: thisOctal
    real(double) :: tval
    type(VECTOR),save :: norm(6)
    logical, save :: firsttime = .true.
    real(double) :: distTor1, theta, mu
    real(double) :: distToRboundary, compz,currentZ
    real(double) :: distToZboundary
    type(VECTOR) ::  zHat, rhat
    real(double) ::  compx, gridRadius
    real(oct) :: t(6),denom(6), r, r1, d, cosmu,x1,x2
    integer :: i,j
    logical :: ok, thisOk(6)
    logical :: debug
    
    real(double) :: subcellsize
    type(VECTOR) :: normdiff

    if(firsttime) then
       norm(1) = VECTOR(1.0d0, 0.d0, 0.0d0)
       norm(2) = VECTOR(0.0d0, 1.0d0, 0.0d0)
       norm(3) = VECTOR(0.0d0, 0.0d0, 1.0d0)
       norm(4) = VECTOR(-1.0d0, 0.0d0, 0.0d0)
       norm(5) = VECTOR(0.0d0, -1.0d0, 0.0d0)
       norm(6) = VECTOR(0.0d0, 0.0d0, -1.0d0)
       firsttime = .false.
    endif

   tval = HUGE(tval)

   point = posVec

   subcen = grid%OctreeRoot%centre

   thisOctal => grid%octreeRoot

   if (thisOctal%oneD) then
   
      r1 = thisOctal%subcellSize*2.d0
      d = modulus(point)
      rHat = (-1.d0)*posVec
      call normalize(rhat)
      theta = asin(max(-1.d0,min(1.d0,r1 / d)))
      cosmu = rHat.dot.direction
      mu = acos(max(-1.d0,min(1.d0,cosmu)))
      distTor1 = 1.e30
      if (mu  < theta ) then
         call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
         if (.not.ok) then
            write(*,*) "Quad solver failed in intersectcubeamr2d"
            x1 = thisoctal%subcellSize
            x2 = 0.d0
         endif
         tval = min(x1,x2)
      endif
      goto 666
   endif

    if (thisOctal%threed.and.(.not.thisOctal%cylindrical)) then

          ! cube

         ok = .true.

         subcellsize = thisOctal%subcellSize

         if(direction%x .ne. 0.d0) then            
            denom(1) = 1.d0 / direction%x
         else
            denom(1) = 0.d0
         endif
         denom(4) = -denom(1)

         if(direction%y .ne. 0.d0) then            
            denom(2) = 1.d0 / direction%y
         else
            denom(2) = 0.d0
         endif
         denom(5) = -denom(2)

         if(direction%z .ne. 0.d0) then            
            denom(3) = 1.d0 / direction%z
         else
            denom(3) = 0.d0
         endif
         denom(6) = -denom(3)

         normdiff = subcen - posvec

         t(1) =  (normdiff%x + subcellsize) * denom(1)
         t(2) =  (normdiff%y + subcellsize) * denom(2)
         t(3) =  (normdiff%z + subcellsize) * denom(3)
         t(4) =  (normdiff%x - subcellsize) * denom(1)
         t(5) =  (normdiff%y - subcellsize) * denom(2)
         t(6) =  (normdiff%z - subcellsize) * denom(3)

         thisOk = .true.

         do i = 1, 6
            if (denom(i) .ge. 0.d0) thisOK(i) = .false.
            if (t(i) < 0.) thisOk(i) = .false.
         enddo
         
         j = 0
         do i = 1, 6
            if (thisOk(i)) j=j+1
         enddo
         
         if (j == 0) ok = .false.
         
         if (.not.ok) then
            write(*,*) "Error: j=0 (no intersection???) in lucy_mod::distancetoGridFromOutside. "
            write(*,*) direction%x,direction%y,direction%z
            write(*,*) t(1:6)
         endif
         
         tval = maxval(t, mask=thisOk)
         
!         write(*,*) t(1:6),thisOK(1:6)

         if (tval == 0.) then
            write(*,*) " tval=0 (no intersection???) in lucy_mod::intersectCubeAMR. "
            write(*,*) posVec
            write(*,*) direction%x,direction%y,direction%z
            write(*,*) t(1:6)
            stop
         endif
         

      else


! now look at the cylindrical case

         ! first do the inside and outside curved surfaces

         if (thisOctal%twoD) then
            r = sqrt(subcen%x**2 + subcen%y**2)
            r1 = r + thisOctal%subcellSize
            gridRadius = r1
         else
            gridRadius = thisOctal%subcellSize
            r1 = gridRadius
         endif


         d = sqrt(point%x**2+point%y**2)
         xHat = (-1.d0)*VECTOR(point%x, point%y,0.d0)
         call normalize(xHat)
         rDirection = VECTOR(direction%x, direction%y, 0.d0)
         compX = modulus(rDirection)
         call normalize(rDirection)
               
         theta = asin(max(-1.d0,min(1.d0,r1 / d)))
         cosmu = xHat.dot.rdirection
         mu = acos(max(-1.d0,min(1.d0,cosmu)))
         distTor1 = 1.e30
         if (compx /= 0.d0) then
            if (mu  < theta ) then
               call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
               if (.not.ok) then
                  write(*,*) "Quad solver failed in intersectcubeamr2d"
                  x1 = thisoctal%subcellSize
                  x2 = 0.d0
               endif
               distTor1 = min(x1,x2)/CompX
               hitVec = posVec + distToR1 * direction
               if (abs(hitVec%z) > thisOctal%subcellSize) distToR1 = 1.d30
            endif
         endif
         
         distToRboundary = distTor1
         if (distToRboundary < 0.d0) then
            distToRboundary = 1.e30
         endif

         ! now do the upper and lower (z axis) surfaces

         zHat = VECTOR(0.d0, 0.d0, 1.d0)
         compZ = zHat.dot.direction
         currentZ = point%z

         debug = .false.
         if(debug) then
            write(*,*) "subcen",subcen
            write(*,*) "thisoctal%subcellsize",thisoctal%subcellsize
            write(*,*) "point", point
            write(*,*) "posvec", posvec
            write(*,*) "direction",direction
         endif

         if (compZ /= 0.d0 ) then
            if (compZ > 0.d0) then
               distToZboundary= (subcen%z - thisOctal%subcellsize - currentZ ) / compZ                              
               hitVec = posvec + disttoZboundary * direction
               if (sqrt(hitVec%x**2 + hitVec%y**2) > gridRadius) distToZboundary = 1.d30
            else
               distToZboundary = abs((subcen%z + thisOctal%subcellsize - currentZ ) / compZ)
               hitVec = posvec + disttoZboundary * direction
               if (sqrt(hitVec%x**2 + hitVec%y**2) > gridRadius) distToZboundary = 1.d30
            endif
         else
            disttoZboundary = 1.e30
         endif
      
         tVal = min(distToZboundary, distToRboundary)
         if(.not. suppressWarnings) then

            if (tVal > 1.e29) then
               write(*,*) "Cylindrical"
               write(*,*) "posVec",posvec
               write(*,*) "direction",direction
               write(*,*) tVal,compX,compZ, distToZboundary,disttorboundary
               write(*,*) "subcen",subcen
               write(*,*) "z",currentZ
               write(*,*) "x1,x2",x1,x2
            endif

            if (tval < 0.) then
               write(*,*) "Cylindrical"
               write(*,*) tVal,distToZboundary,disttorboundary
               write(*,*) "subcen",subcen
               write(*,*) "z",currentZ
            endif
         endif
      endif
      
666   continue
    end function distanceToGridFromOutside

  !
  ! hydroWarp stuff
  !
  function hydroWarpTemperature(point, grid) result (tempOut)
    use input_variables, only : warpFracHeight, warpRadius, warpSigma, warpAngle
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector), INTENT(IN) :: point
    real(double) :: tempOut
    TYPE(vector) :: rVec
    real(double) :: r, phi, warpHeight
    real :: thisTemp

    tempOut = 20.

    rVec = point
    r = sqrt(rVec%x**2 + rVec%y**2)
    phi = atan2(rVec%y,rVec%x)
    warpheight  = cos(phi+warpAngle) * warpFracHeight * warpradius * exp(-0.5d0*((r - warpRadius)/warpSigma)**2)
    rVec%z = rVec%z - warpheight

    if (inOctal(grid%hydroGrid%octreeRoot, rVec)) then
      call amrGridValues(grid%hydroGrid%octreeRoot, rVec, temperature=thisTemp)
      tempOut = thisTemp
    end if
  end function hydroWarpTemperature

  function hydroWarpDensity(point, grid) result (rhoOut)
    use input_variables, only : warpFracHeight, warpRadius, warpSigma, warpAngle
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(vector), INTENT(IN) :: point
    real(double) :: rhoOut
    TYPE(vector) :: rVec
    real(double) :: r, phi, warpHeight

    integer :: nx, thisX
    type(hydroSpline) :: thisSpline

    rhoOut = 1.d-30

    rVec = point
    r = sqrt(rVec%x**2 + rVec%y**2)
    phi = atan2(rVec%y,rVec%x)
    warpheight  = cos(phi+warpAngle) * warpFracHeight * warpradius * exp(-0.5d0*((r - warpRadius)/warpSigma)**2)
    rVec%z = rVec%z - warpheight

    ! Determine which x-value most closely matches our r-value
    nx = size(grid%hydroGrid%hydroSplines)
    call locate_hydroSpline(grid%hydroGrid,nx,r,thisX)
    if (thisX < 1) then
      thisX = 1
    else if (thisX < nx) then
      if (abs(grid%hydroGrid%hydroSplines(thisX+1)%x-r) .lt. abs(r-grid%hydroGrid%hydroSplines(thisX)%x)) thisX = thisX + 1
    end if

    ! Interpolate over the appropriate spline to get the log(rho) at the point
    thisSpline = grid%hydroGrid%hydroSplines(thisX)
    call splint(thisSpline%z, thisSpline%rho, thisSpline%rhoDD, thisSpline%nz, rVec%z, rhoOut)

    ! The spline function stores and interpolates on the log of the density, so we
    ! must unlog the densities to give our result.
    rhoOut = 10**rhoOut
  end function hydroWarpDensity

  function getHydroWarpScaleHeight(r, grid) result (scaleheight)
    type(gridtype), intent(in) :: grid
    real(double), intent(in) :: r
    real(double) :: scaleheight

    integer :: nx, thisX

    ! Determine which x-value most closely matches our r-value
    nx = size(grid%hydroGrid%hydroSplines)
    call locate_hydroSpline(grid%hydroGrid,nx,r,thisX)
    if (thisX < 1) then
      thisX = 1
    else if (thisX < nx) then
      if (abs(grid%hydroGrid%hydroSplines(thisX+1)%x-r) .lt. abs(r-grid%hydroGrid%hydroSplines(thisX)%x)) thisX = thisX + 1
    end if

    scaleheight = grid%hydroGrid%hydroSplines(thisX)%scaleheight
  end function getHydroWarpScaleHeight

  subroutine hydroWarpFitSplines(grid)
    use input_variables, only: rOuter

    type(GRIDTYPE) :: grid
    real :: xAxis(1000000)
    real(double) :: zAxis(10000), rho(10000), subcellsize(10000)
    integer :: nx, nz
    real :: xPos, yPos
    integer :: i, j

    real(double) :: newrho(10000), newzAxis(10000), newrhodd(10000)
    integer :: newnz
    zAxis = 0
    nx = 0
    nz = 0
    newrhodd = 0.d0 ; rho = 0.d0; subcellSize = 0.d0
    call getxValuesAMR(grid%octreeRoot, nx, xAxis)
    call stripSimilarValues(xAxis,nx,real(1.d-5*grid%halfSmallestSubcell))
    do while (xAxis(nx) > rOuter)
      nx = nx - 1
    end do
    xAxis(1:nx) = xAxis(1:nx) + 1.d-5*grid%halfSmallestSubcell

    allocate(grid%hydroSplines(nx))

    yPos = 0.
    do i = 1, nx
       xPos = xAxis(i)
!       if ((xPos > grid%rInner).and.(xPos < grid%rOuter)) then
          ! Read in a complete density run for this x position
          call getDensityRun(grid, zAxis, subcellsize, rho, xPos, yPos, nz, -1.)
          do j = 1, nz
            newRho(j) = rho(nz+1-j)
            newzAxis(j) = -1.*zAxis(nz+1-j)
          end do
          newnz = nz
          call getDensityRun(grid, zAxis, subcellsize, rho, xPos, yPos, nz, +1.)
          newRho(newnz+1:newnz+nz) = rho(1:nz)
          newzAxis(newnz+1:newnz+nz) = zAxis(1:nz)
          newnz = newnz + nz

          ! We use the log of the densities to avoid numerical funnies in the
          ! spline fitting.
          do j = 1, newnz
            newRho(j) = log10(newRho(j))
          enddo

          ! Calculate the spline function
          call spline(newzAxis, newRho, newnz, 1.d30, 1.d30, newRhoDD)

          ! Store the spline in the appropriate array
          allocate(grid%hydroSplines(i)%z(newnz))
          allocate(grid%hydroSplines(i)%rho(newnz))
          allocate(grid%hydroSplines(i)%rhoDD(newnz))
          grid%hydroSplines(i)%x = xPos
          grid%hydroSplines(i)%nz = newnz
          grid%hydroSplines(i)%z(1:newnz) = newzAxis(1:newnz)
          grid%hydroSplines(i)%rho(1:newnz) = newRho(1:newnz)
          grid%hydroSplines(i)%rhoDD(1:newnz) = newRhoDD(1:newnz)
          ! Set the scaleheight
          call hydroWarpScaleHeight(grid%hydroSplines(i))
!       endif
    enddo
  end subroutine hydroWarpFitSplines

  subroutine getDensityRun(grid, zAxis, subcellsize, rho, xPos, yPos, nz, direction)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer :: nz
    real(double) :: rho(:)
    real(double) :: zAxis(:), subcellsize(:)
    real :: xPos, yPos
    integer :: subcell
    real(double) :: rhotemp 
    real :: direction
    type(VECTOR) :: currentPos, temp
    real :: halfSmallestSubcell

    nz = 0
    halfSmallestSubcell = grid%halfSmallestSubcell

    currentPos = VECTOR(xPos, yPos, direction*halfSmallestSubcell)

    do while(abs(currentPos%z) < grid%ocTreeRoot%subcellsize)
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp)
!       if (thisOctal%inFlow(subcell)) then
          nz = nz + 1
          rho(nz) = rhotemp
          temp = subCellCentre(thisOctal, subcell)
          zAxis(nz) = temp%z
          subcellsize(nz) = thisOctal%subcellsize
!       endif
          currentPos = VECTOR(xPos, yPos, zAxis(nz)+0.5*direction*thisOctal%subcellsize+direction*halfSmallestSubcell)
    end do
    zAxis(1:nz) = abs(zAxis(1:nz)) !* 1.d10  ! convert to cm
  end subroutine getDensityRun

    SUBROUTINE LOCATE_hydroSpline(grid,N,X,J)
    type(gridtype), intent(in) :: grid
    integer, intent(in)              :: n
    real(double), intent(in) :: x
    integer,intent(out)              :: j
    integer :: jl, ju,jm
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((grid%hydroSplines(N)%x.GT.grid%hydroSplines(1)%x).EQV.(X.GE.grid%hydroSplines(JM)%x))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF

      ! Will force to be between 1 and the array size -1.
      if(x <= grid%hydroSplines(1)%x)then
        j=1
      else if(x>=grid%hydroSplines(n)%x)then
        j=n-1
      else
        j=jl
      end if

      return

    END SUBROUTINE LOCATE_hydroSpline

  subroutine hydroWarpScaleHeight(spline)
    type(hydroSpline) :: spline
    real(double) :: max_rho, max_rho_z
    integer :: max_rho_i
    real(double) :: target_rho

    integer :: i

    integer :: sv1_i, sv2_i
    real(double) :: sv1, sv2, step1, step2
    real(double) :: sh1, sh2

    call hydroWarpSplinePeak(spline, max_rho, max_rho_z, max_rho_i)
    target_rho = max_rho - log10(2.718281828)

! improve this conditional
    if (max_rho .lt. log10(1.d-29)) then
      ! We have a straight line
      spline%scaleHeight = 0.
      return
    endif

    ! assume the spline is fairly stable around one scaleheight from the peak
    ! find the cell before which we pass through the target density
    sv1_i = 0
    sv2_i = 0
    do i=1, max_rho_i
      if (spline%rho(i) .gt. target_rho) then
        if (i == 1) then
          print *, "Error finding scaleheight: sv1_i = 1!"
          stop
        end if
        sv1_i = i-1
        sv1 = spline%z(sv1_i)
        step1 = abs(spline%z(sv1_i+1)-sv1)/2.
        exit
      end if
    end do
    if (sv1_i == 0) then
      sv1_i = max_rho_i
      sv1 = spline%z(sv1_i)
      step1 = abs(max_rho_z-sv1)/2.
    end if

    do i=spline%nz, max_rho_i+1, -1
      if (spline%rho(i) .gt. target_rho) then
        if (i == spline%nz) then
          print *, "Error finding scaleheight: sv2_i = nz!"
          stop
        end if
        sv2_i = i+1
        sv2 = spline%z(sv2_i)
        step2 = -1.*abs(sv2-spline%z(sv2_i-1))/2.
        exit
      end if
    end do
    if (sv2_i == 0) then
      sv2_i = max_rho_i+1
      sv2 = spline%z(sv2_i)
      step2 = -1.*abs(sv2-max_rho_z)/2.
    end if

    sh1 = hydroWarpSplinty(spline, target_rho, sv1, step1)
    sh2 = hydroWarpSplinty(spline, target_rho, sv2, step2)

    spline%scaleHeight = abs(sh2-sh1)/2.

    ! perform a quick sanity check
    if (abs(sh1+sh2-2.*max_rho_z)/spline%scaleHeight .gt. 0.05) then
      print *, "Warning: Scaleheights either side of midplane differ by"
      print *, "         more than 5%:", sh1, sh2
    end if
  end subroutine hydroWarpScaleHeight

  real(double) function hydroWarpSplinty(spline, target, startval, step)
    implicit none

    type(hydroSpline) :: spline
    real(double) :: target    ! target density value
    real(double) :: startval  ! initial value of z to try
    real(double) :: step      ! initial stepping in value of z

    ! initial values
    real(double) :: tol = 0.0001   ! +/- tol * target
    real(double) :: reduxFac = 0.5 ! step is reduced by this factor when homing in
    integer :: maxIter = 100       ! maximum number of iterations before we give up
    integer :: i                   ! iteration counter
   
    real(double) :: xNew, yNew, xOld, yOld

    xNew = startval
    call splint(spline%z, spline%rho, spline%rhoDD, spline%nz, xNew, yNew)
    i = 0
    do
      i = i + 1
      if (i > maxIter) then
         write (*,*) "Numerical solver: exceeded maximum iterations allowed (", maxIter, ")"
         stop
      end if

      xOld = xNew
      yOld = yNew
      xNew = xNew + step
      call splint(spline%z, spline%rho, spline%rhoDD, spline%nz, xNew, yNew)
      if (abs(yNew-target) < abs(tol * target)) then
         exit
      else if (yNew < target) then
         if (yNew < yOld) then
            if (yOld < target) then
               step = -1. * step
               xNew = xOld
               yNew = yOld
            else ! (yOld > target)
               step = -1. * (reduxFac * step)
            end if
         end if
      else ! (yNew > target)
         if (yNew > yOld) then
            if (yOld < target) then
               step = -1. * (reduxFac * step)
            else ! (yOld > target)
               step = -1. * step
               xNew = xOld
               yNew = yOld
            end if
         end if
      end if
    end do

    hydroWarpSplinty = xNew
  end function hydroWarpSplinty

  ! returns, the max density, z value of max density, and index of element
  ! immediately preceeding the max density location
  subroutine hydroWarpSplinePeak(spline, rho_max, rho_max_z, rho_max_i)
    type(hydroSpline), intent(in) :: spline
    real(double), intent(out) :: rho_max, rho_max_z
    integer, intent(out) :: rho_max_i

    integer :: i

    real(double) :: step, z_test, rho_test
    integer :: nsteps
    nsteps = 100

    rho_test = 0.d0
    ! desnity is log(10) density, so we need a very low value here
    rho_max = -350.
    rho_max_i = 0

    ! find cell with the max density
    do i=1, spline%nz
      if (spline%rho(i) .gt. rho_max) then
        rho_max = spline%rho(i)
        rho_max_i = i
      end if
    end do
    rho_max_z = spline%z(rho_max_i)
    
    ! the peak in the spline should be somewhere between this max cell value
    ! and the two neighbouring cells
    if (rho_max_i == 1) then
      step = (spline%z(2) - spline%z(1)) / (nsteps-1)
      z_test = spline%z(1)
    else if (rho_max_i == spline%nz) then
      step = (spline%z(rho_max_i) - spline%z(rho_max_i-1)) / (nsteps-1)
      z_test = spline%z(rho_max_i-1)
    else
      step = (spline%z(rho_max_i+1) - spline%z(rho_max_i-1)) / (nsteps-1)
      z_test = spline%z(rho_max_i-1)
    end if

    do i=1, nsteps
        call splint(spline%z, spline%rho, spline%rhoDD, spline%nz, z_test, rho_test)
        if (rho_test .gt. rho_max) then
          rho_max = rho_test
          rho_max_z = z_test
        end if
        z_test = z_test + step
    end do

    ! set the index for the density peak to the index of the element
    ! immediately preceeding the peak value
    if (rho_max_z < spline%z(rho_max_i)) rho_max_i = rho_max_i - 1
  end subroutine hydroWarpSplinePeak
  !
  ! end hydroWarp stuff
  !

  recursive subroutine splitGridOnStream(thisOctal, thisStream, grid, converged)

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    type(VECTOR) :: rVec !, corner
    integer :: i, j, subcell
    logical :: converged, converged_tmp
    logical :: split
    type(STREAMTYPE), pointer :: thisStream
    real(double), parameter :: fac = 2.d0
    converged = .true.
    converged_tmp=.true.

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          children : do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call splitGridOnStream(child, thisStream, grid, converged)
                if (.not.converged) converged_tmp = converged
                exit children
             end if
          end do children
          if (.not.converged_tmp) converged=converged_tmp
       else

          rVec = subcellCentre(thisOctal, subcell)

          stream : do j = 1, thisStream%nSamples
             
             
             split = .false.
             !corner = closestCorner(thisOctal, subcell, thisStream%position(j))
             
             if ((inSubcell(thisOctal, subcell, thisStream%position(j))).and. &
                  (thisOctal%subcellSize > thisStream%streamRadius(j)/fac) ) then 
                split = .true.
             else
                
                if ((dist_from_closestCorner(thisOctal, subcell, thisStream%position(j)) < &
                     thisStream%StreamRadius(j)).and. &
                     (thisOctal%subcellSize > thisStream%streamRadius(j)/fac) ) split = .true.
             endif
             
             if (split) then

!                write(*,*) "adding child",thisStream%streamRadius(j),thisOctal%subcellSize
                call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                     inherit=.false., interp=.false.)
                
                converged = .false.
                exit stream
             endif
          enddo stream 
          if (.not.converged) exit 
       end if
    enddo

    return
    
  end subroutine splitGridOnStream

  recursive subroutine splitGridOnStream2(thisOctal, thisStream, grid, childrenAdded)

    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    type(VECTOR) :: rVec
    integer :: i, j, subcell
    logical :: childrenAdded
    logical :: split
    type(STREAMTYPE) :: thisStream
    real(double), parameter :: fac = 2.d0

    subcell = 1
    do while (subcell <= thisOctal%maxChildren)
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          children : do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call splitGridOnStream2(child, thisStream, grid, childrenAdded)
                exit children
             end if
          end do children
       else


          childrenAdded = .false.
          rVec = subcellCentre(thisOctal, subcell)

          stream : do j = 1, thisStream%nSamples
             
             
             split = .false.
             


             if ((inSubcell(thisOctal, subcell, thisStream%position(j))).and. &
                  (thisOctal%subcellSize > thisStream%streamRadius(j)/fac) ) then 
                split = .true.
             else
                
                if ((dist_from_closestEdge(thisOctal, subcell, thisStream%position(j)) < &
                     thisStream%StreamRadius(j)).and. &
                     (thisOctal%subcellSize > thisStream%streamRadius(j)/fac) ) split = .true.
             endif
             
             if (split) then

                call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                     inherit=.false., interp=.false.)
                childrenAdded = .true.
                subcell = subcell - 1
                exit stream
             endif
          enddo stream
       end if
       subcell = subcell + 1
    enddo

  end subroutine splitGridOnStream2

  recursive subroutine splitGridOnStream3(thisOctal, grid, thisStream)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    integer :: i, subcell
    logical :: split
    type(STREAMTYPE) :: thisStream
    type(STREAMTYPE) :: octalStream, subcellStream
!    real(double) :: side, rStar
    logical :: splitAz

    call createOctalStream(thisOctal, thisStream, octalStream)

    subcell = 1
    do while (subcell <= thisOctal%maxChildren)
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          children : do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call splitGridOnStream3(child, grid, octalStream)
                exit children
             end if
          end do children
       else

          split = .false.
          call createSubcellStream(thisOctal, subcell, octalStream, subcellStream)
!          if (subcellStream%nSamples > 1) split = .true.
 

          split = .true.
!          rStar = 2.d0*rSol/1.d10
!          side = 0.1d0 * (thisOctal%r/rStar)**3
!          if (thisOctal%subcellSize > side) split = .true.

          if (thisOctal%dphi*radtodeg > 15.) then
             splitAz = .true.
          else
             splitAz = .false.
          endif

          if (thisOctal%nDepth > 7) split = .false.

          if (subcellstream%nsamples == 0) split = .false.

          if (split) then
             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                     inherit=.false., interp=.false., stream=subcellStream, splitAzimuthally=splitAz)
             subcell = subcell - 1
          endif
       end if
       subcell = subcell + 1
       call freeStream(subcellStream)
    enddo
 

    call freeStream(octalStream)
  end subroutine splitGridOnStream3




  subroutine allocateStream(thisStream, nSamples)
    type(STREAMTYPE) :: thisStream
    integer :: nSamples
!    call freeStream(thisStream)
    allocate(thisStream%position(1:nSamples))
    allocate(thisStream%direction(1:nSamples))
    allocate(thisStream%velocity(1:nSamples))
    allocate(thisStream%rho(1:nSamples))
    allocate(thisStream%temperature(1:nSamples))
    allocate(thisStream%speed(1:nSamples))
    allocate(thisStream%distanceAlongStream(1:nSamples))
    allocate(thisStream%streamRadius(1:nSamples))


  end subroutine allocateStream

  subroutine freeStream(thisStream)
    type(STREAMTYPE) :: thisStream
    if (associated(thisStream%position)) deallocate(thisStream%position)
    if (associated(thisStream%velocity)) deallocate(thisStream%velocity)
    if (associated(thisStream%rho)) deallocate(thisStream%rho)
    if (associated(thisStream%temperature)) deallocate(thisStream%temperature)

    if (associated(thisStream%direction)) deallocate(thisStream%direction)
    if (associated(thisStream%speed)) deallocate(thisStream%speed)
    if (associated(thisStream%streamRadius)) deallocate(thisStream%streamRadius)
    if (associated(thisStream%distanceAlongStream)) deallocate(thisStream%distanceAlongStream)

  end subroutine freeStream

  subroutine createOctalStream(thisOctal,  thisStream, octalStream)
    type(OCTAL),pointer :: thisOctal
    type(STREAMTYPE) :: thisStream, octalStream
    integer :: i

    call allocateStream(octalStream, thisStream%nSamples)

    octalStream%nSamples = 0
    do i = 1, thisStream%nSamples
       if (inOctal(thisOctal, thisStream%position(i))) then
          octalStream%nSamples = octalStream%nSamples + 1
          octalStream%rho(octalStream%nSamples) = thisStream%rho(i)
          octalStream%position(octalStream%nSamples) = thisStream%position(i)
          octalStream%velocity(octalStream%nSamples) = thisStream%velocity(i)
          octalStream%temperature(octalStream%nSamples) = thisStream%temperature(i)
       endif
    end do
  end subroutine createOctalStream

  subroutine createSubcellStream(thisOctal,  subcell, thisStream, octalStream)
    integer :: subcell
    type(OCTAL),pointer :: thisOctal
    type(STREAMTYPE) :: thisStream, octalStream
    integer :: i
    call allocateStream(octalStream, thisStream%nSamples)
    octalStream%nSamples = 0 
    do i = 1, thisStream%nSamples
       if (inSubcell(thisOctal, subcell, thisStream%position(i))) then
          octalStream%nSamples = octalStream%nSamples + 1
          octalStream%rho(octalStream%nSamples) = thisStream%rho(i)
          octalStream%position(octalStream%nSamples) = thisStream%position(i)
          octalStream%velocity(octalStream%nSamples) = thisStream%velocity(i)
          octalStream%temperature(octalStream%nSamples) = thisStream%temperature(i)
       endif
    end do
  end subroutine createSubcellStream


  subroutine readStreams(thisStream, nStream, filename)
    use input_variables, only : scaleFlowRho, flowrhoScale
    type(STREAMTYPE) :: thisStream(:)
    character(len=*) :: filename
    integer :: nstream
    integer :: i, j
    real(double) :: r, theta, phi, v, rho
    real(double) :: tot
    type(VECTOR) :: rVec
    open(20, file=filename, status="old", form="formatted")
    nStream = 0


10 continue
       read(20,*,end=20) r, theta, phi, v, rho
       if (r < 1.00001d0) then
          if (nStream /= 0) then
             if (thisStream(nStream)%nSamples == 1) then
                write(*,*) "Stream ",nstream, " only has one sample!!!!"
                nStream = nStream - 1
             endif
          endif
          nStream = nStream + 1
          thisStream(nStream)%nSamples = 0
          call allocateStream(thisStream(nStream), 200)
       endif

          theta = theta !* degtorad
          phi = phi ! * degtorad
          rho = rho * 1.d-3
          r = r * 2.d0 * rsol / 1.d10
          rVec = r * VECTOR(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
          
          thisStream(nStream)%nSamples = thisStream(nStream)%nSamples + 1
          
          thisStream(nStream)%position(thisStream(nStream)%nSamples) = rVec
          thisStream(nStream)%speed(thisStream(nStream)%nSamples) = v * 1.e5/cspeed
          if (scaleFlowRho) then
             thisStream(nStream)%rho(thisStream(nStream)%nSamples) = rho * flowRhoScale
          else
             thisStream(nStream)%rho(thisStream(nStream)%nSamples) = rho
          endif
          thisStream(nStream)%temperature(thisStream(nStream)%nSamples) = 7500.
          thisStream(nStream)%streamRadius(thisStream(nStream)%nSamples)  = 0.01d0
    goto 10
20  continue
    close(20)
    write(*,*) "Read in ",nStream, " streams."
    do i = 1, nStream
       do j = 1, thisStream(i)%nSamples
          if (j == thisStream(i)%nSamples) then
             thisStream(i)%direction(j) = thisStream(i)%position(j) - thisStream(i)%position(j-1)
             call normalize(thisStream(i)%direction(j))
          
          else
             thisStream(i)%direction(j) = thisStream(i)%position(j+1) - thisStream(i)%position(j)
             call normalize(thisStream(i)%direction(j))
          endif
             if (j == 1) then
                thisStream(i)%distanceAlongStream(j) = 0.d0
             else
                thisStream(i)%distanceAlongStream(j) = thisStream(i)%distanceAlongStream(j-1) + &
                     modulus(thisStream(i)%position(j) - thisStream(i)%position(j-1))
             endif
          ! direction is outward so need to use -speed...
          thisStream(i)%velocity(j) = (-1.d0*thisStream(i)%speed(j)) * thisStream(i)%direction(j)
       enddo
    enddo

! scale temperature according to distance along the stream
    do i = 1, nStream
       tot = MAXVAL(thisStream(i)%rho(1:thisStream(i)%nSamples))
       do j = 1, thisStream(i)%nSamples
          tot = 1.d0-(thisStream(i)%distanceAlongStream(j) / &
               thisStream(i)%distanceAlongStream(thisStream(i)%nSamples))
          thisStream(i)%temperature(j) = 4000.d0 + 4000.d0 * tot

!          thisStream(i)%temperature(j) = 6000.d0
   
          ! Behaviour similar to Martin96 fig 1b
          ! Arbitrary !!!!
!          thisStream(i)%temperature(j) = min(2850.d0  + (tot/0.8)**3  * 3000.d0, 6000.d0)
       enddo
    enddo

    ! Analytical description of the T along the stream 
    ! reproduces behaviour calculated by Martin 1996 fig 1b
!    do i = 1, nStream
!       tot = MAXVAL(thisStream(i)%rho(1:thisStream(i)%nSamples))
!       do j = 1, thisStream(i)%nSamples
!          ! evolution is adiabatic up to 6000 K and then isothermal.
!          ! adiabatic exponent is gamma = 5/3 for monoatomic gas (T.V^(gamma-1) = cst)
!          thisStream(i)%temperature(j) = min(3000.d0 * (thisStream(i)%rho(j) / thisStream(i)%rho(1))**(2./3.), 6000.d0)
!       enddo
!    enddo
    

    do i = 1, nStream
       call allocateStream(globalStream(i),thisStream(i)%nSamples)
       globalStream(i) = thisStream(i)
    enddo
    globalnStream = nStream

  end subroutine readStreams

  subroutine readCurtain(thisCurtain,filename)

    implicit none

    type(curtaintype) :: thisCurtain
    character(len=*) :: filename

    integer :: i,j

    real(double), parameter :: r0 = 4.e9, rho0 = 6.9e-12, v0 = 1.63e5

    open(unit=1,file=filename,status="old")
    read(1,*) thisCurtain%nr
    read(1,*) thisCurtain%ntheta
    read(1,*)
    do i=1, thisCurtain%nr
       do j=1, thisCurtain%ntheta
          read(1,*) thisCurtain%position(i,j)%x, thisCurtain%position(i,j)%z, &
               thisCurtain%density(i,j), thisCurtain%velocity(i,j)%x, &
               thisCurtain%velocity(i,j)%z, thisCurtain%velocity(i,j)%y, &
               thisCurtain%temperature(i,j)
       enddo !j
    enddo !i

    ! SI system
    thisCurtain%position(:,:)%x = thisCurtain%position(:,:)%x * r0
    thisCurtain%position(:,:)%z = thisCurtain%position(:,:)%z * r0
    thisCurtain%velocity(:,:)%x = thisCurtain%velocity(:,:)%x * v0
    thisCurtain%velocity(:,:)%y = thisCurtain%velocity(:,:)%y * v0
    thisCurtain%velocity(:,:)%z = thisCurtain%velocity(:,:)%z * v0
    thisCurtain%temperature(:,:) = thisCurtain%temperature(:,:) / thisCurtain%density(:,:) 
    thisCurtain%density(:,:) = thisCurtain%density(:,:) * rho0

    close(unit=1)

    return

  end subroutine readCurtain


  function  closestCorner(thisOctal, subcell, posVec) result (closeCorner)
    type(VECTOR) :: closeCorner, posVec, cen, thisCorner
    type(OCTAL), pointer :: thisOctal
    real(double) :: minDist, r, dist
    integer :: subcell
    
    minDist =  1.d30
    r = thisOctal%subcellSize/2.d0
    cen = subcellCentre(thisOctal, subcell)

    thisCorner = cen + r * VECTOR(1.d0, 1.d0, 1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
       closeCorner = thisCorner
    endif
    thisCorner = cen + r * VECTOR(1.d0, 1.d0, -1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
       closeCorner = thisCorner
    endif
    thisCorner = cen + r * VECTOR(1.d0, -1.d0, -1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
       closeCorner = thisCorner
    endif
    thisCorner = cen + r * VECTOR(1.d0, -1.d0, 1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
       closeCorner = thisCorner
    endif
    thisCorner = cen + r * VECTOR(-1.d0, -1.d0,  -1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
       closeCorner = thisCorner
    endif
    thisCorner = cen + r * VECTOR(-1.d0, 1.d0, -1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
       closeCorner = thisCorner
    endif
    thisCorner = cen + r * VECTOR(-1.d0, -1.d0, 1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
       closeCorner = thisCorner
    endif
    thisCorner = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
       closeCorner = thisCorner
    endif
  end function closestCorner

  function  dist_from_closestCorner(thisOctal, subcell, posVec) result (minDist)
    type(VECTOR) :: posVec, cen, thisCorner
    type(OCTAL), pointer :: thisOctal
    real(double) :: minDist, r, dist
    integer :: subcell
    
    minDist =  1.d30
    r = thisOctal%subcellSize/2.d0
    cen = subcellCentre(thisOctal, subcell)

    thisCorner = cen + r * VECTOR(1.d0, 1.d0, 1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
    endif
    thisCorner = cen + r * VECTOR(1.d0, 1.d0, -1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
    endif
    thisCorner = cen + r * VECTOR(1.d0, -1.d0, -1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
    endif
    thisCorner = cen + r * VECTOR(1.d0, -1.d0, 1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
    endif
    thisCorner = cen + r * VECTOR(-1.d0, -1.d0,  -1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
    endif
    thisCorner = cen + r * VECTOR(-1.d0, 1.d0, -1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
    endif
    thisCorner = cen + r * VECTOR(-1.d0, -1.d0, 1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
    endif
    thisCorner = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
    dist  =  modulus(posVec - thisCorner)
    if (dist < minDist) then
       minDist = dist
    endif
  end function dist_from_closestCorner

  function  dist_from_closestEdge(thisOctal, subcell, posVec) result (minDist)
    type(VECTOR) :: posVec, cen, corner1, corner2
    type(OCTAL), pointer :: thisOctal
    real(double) :: minDist, r, dist
    integer :: subcell
    
    minDist =  1.d30
    r = thisOctal%subcellSize/2.d0
    cen = subcellCentre(thisOctal, subcell)

!1
    corner1 = cen + r * VECTOR(-1.d0,-1.d0, 1.d0)
    corner2 = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!2
    corner1 = cen + r * VECTOR(-1.d0,-1.d0, 1.d0)
    corner2 = cen + r * VECTOR( 1.d0,-1.d0, 1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!3
    corner1 = cen + r * VECTOR(-1.d0,-1.d0, 1.d0)
    corner2 = cen + r * VECTOR(-1.d0,-1.d0,-1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!4
    corner1 = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
    corner2 = cen + r * VECTOR( 1.d0, 1.d0, 1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!5
    corner1 = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
    corner2 = cen + r * VECTOR(-1.d0, 1.d0,-1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!6
    corner1 = cen + r * VECTOR( 1.d0,-1.d0, 1.d0)
    corner2 = cen + r * VECTOR( 1.d0, 1.d0, 1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!7
    corner1 = cen + r * VECTOR( 1.d0,-1.d0, 1.d0)
    corner2 = cen + r * VECTOR( 1.d0,-1.d0,-1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!8
    corner1 = cen + r * VECTOR( 1.d0, 1.d0, 1.d0)
    corner2 = cen + r * VECTOR( 1.d0, 1.d0,-1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!9
    corner1 = cen + r * VECTOR(-1.d0, 1.d0,-1.d0)
    corner2 = cen + r * VECTOR( 1.d0, 1.d0,-1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!10
    corner1 = cen + r * VECTOR(-1.d0, 1.d0,-1.d0)
    corner2 = cen + r * VECTOR(-1.d0,-1.d0,-1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!11
    corner1 = cen + r * VECTOR(-1.d0,-1.d0,-1.d0)
    corner2 = cen + r * VECTOR( 1.d0,-1.d0,-1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif
!12
    corner1 = cen + r * VECTOR( 1.d0,-1.d0,-1.d0)
    corner2 = cen + r * VECTOR( 1.d0, 1.d0,-1.d0)
    dist  =  distancePointLineSegment(corner1, corner2, posVec)
    if (dist < minDist) then
       minDist = dist
    endif


  end function dist_from_closestEdge


  subroutine tauAlongPath(ilambda, grid, rVec, direction, tau, tauMax, ross, startOctal, startSubcell, nTau, xArray, tauArray)
    use input_variables, only : rGap
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, direction, currentPosition, beforeVec, afterVec
    real(double), optional :: xArray(:), tauArray(:)
    integer, optional :: nTau
    integer :: iLambda
    real(double) :: tau, distToNextCell
    real(double), optional :: tauMax
    type(OCTAL), pointer :: thisOctal, sOctal
    type(OCTAL), pointer, optional :: startOctal
    integer, optional :: startSubcell
    real(double) :: fudgeFac = 1.d-1
    real(double) :: kappaSca, kappaAbs, kappaExt
    integer :: subcell
    logical, optional :: ross
    logical :: planetGap
    kappaAbs = 0.d0; kappaSca = 0.d0
    tau = 0.d0
    currentPosition = rVec
    if (PRESENT(nTau)) then
       xArray = 0.d0
       tauArray = 0.d0
       ntau = 1
    endif
    planetGap  = .false.
    if (grid%geometry == "planetgap") planetgap = .true.

    if (PRESENT(startOctal)) then
       thisOctal => startOctal
       subcell = startSubcell
!       if (.not.inOctal(thisOctal, currentPosition)) then
!          write(*,*) "bug in startreturnsamples2"
!          stop
!       endif
    else
       CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
    endif

    do while (inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition, thisOctal,subcell)
       if (.not.PRESENT(ross)) then
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
          kappaExt = kappaAbs + kappaSca
       else
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, rosselandKappa=kappaExt)
          kappaExt = kappaExt * thisOctal%rho(subcell) * 1.d10
       endif
       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
       
       beforeVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       afterVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       
       if (thisOctal%twod.and.(direction%x < 0.d0).and.(currentPosition%x < 0.d0)) exit

       if ((thisOctal%cylindrical).and.((beforeVec.dot.afterVec).lt.0.d0)) then
          exit
       endif
          
       if (planetGap) then
          if ((direction%x < 0.d0).and.(rVec%x > rGap*autocm/1.d10) &
               .and.(currentPosition%x < rGap*autocm/1.d10)) exit
          if ((direction%x > 0.d0).and.(rVec%x < rGap*autocm/1.d10) &
               .and.(currentPosition%x > rGap*autocm/1.d10)) exit
       endif

       tau = tau + distToNextCell*kappaExt
       if (PRESENT(nTau)) then
          nTau = nTau + 1
          xArray(nTau) = xArray(nTau-1) + distToNextCell
          tauArray(nTau) = tau
       endif
       if (PRESENT(tauMax)) then
          if (tau > tauMax) exit
       endif
    end do
  end subroutine tauAlongPath

  subroutine tauAlongPath2(ilambda, grid, rVec, direction, tau, tauMax, ross, startOctal, startSubcell)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, direction, currentPosition, beforeVec, afterVec
    integer :: iLambda
    real(double) :: tau, distToNextCell
    real(double), optional :: tauMax
    type(OCTAL), pointer :: thisOctal, sOctal
    type(OCTAL), pointer, optional :: startOctal
    integer, optional :: startSubcell
    real(double) :: fudgeFac = 1.d-1
    real(double) :: kappaSca, kappaAbs, kappaExt
    integer :: subcell
    logical, optional :: ross
    kappaAbs = 0.d0; kappaSca = 0.d0
    tau = 0.d0
    currentPosition = rVec
    

    if (PRESENT(startOctal)) then
       thisOctal => startOctal
       subcell = startSubcell
    else
       CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
    endif

    do while (inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition, thisOctal,subcell)
       if (.not.PRESENT(ross)) then
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)
          kappaExt = kappaAbs + kappaSca
       else
          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda, rosselandKappa=kappaExt)
          kappaExt = kappaExt * thisOctal%rho(subcell) * 1.d10
       endif
       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)
       
       beforeVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*direction
       afterVec = VECTOR(currentPosition%x, currentPosition%y,0.d0)
       

       tau = tau + distToNextCell*kappaExt
       if (PRESENT(tauMax)) then
          if (tau > tauMax) exit
       endif
    end do
  end subroutine tauAlongPath2


  subroutine getMagStreamValues2(point, grid,  rho, temperature, &
       velocity, inflow)
    type(STREAMTYPE) :: thisStream
    type(GRIDTYPE) :: grid
    integer :: iSample
    type(VECTOR) :: point
    real(double) :: t, thisR
    real(double),optional :: rho
    real,optional :: temperature
    type(VECTOR), optional :: velocity
    logical(kind=logic), optional :: inFlow
    logical :: outsideStream
    outsideStream = .false.
    iSample = 0
    t = 0.d0

    call findNearestSample(thisStream, point, iSample, t, outsideStream)


    thisR =  thisStream%streamradius(isample) + &
         t * (thisStream%streamradius(iSample+1)-thisStream%streamradius(iSample))

!    outsideStream = .false.
!    if (modulus(thisStream%position(isample)-point) > thisR) then
!       outSideStream = .true.
!    endif
!
!    if (iSample == thisStream%nSamples) then
!       if (((point-thisStream%position(iSample)).dot.thisStream%position(iSample)) < 0.d0) then
!          outsideStream = .true.
!       endif
!    endif

    if (PRESENT(rho)) then
       if (outSideStream) then
          rho = 1.d-25
       else
          rho = thisStream%rho(isample) + t * (thisStream%rho(iSample+1)-thisStream%rho(iSample))
       endif
    endif

    if (PRESENT(temperature)) then
       if (outSideStream) then
          temperature = 6000.
       else
          temperature = thisStream%temperature(isample) + &
               t * (thisStream%temperature(iSample+1)-thisStream%temperature(iSample))
       endif
    endif

    if (PRESENT(velocity)) then
       if (outsideStream) then
          velocity = VECTOR(1.e-30, 1.e-30, 1.e-30)
       else
          velocity = thisStream%velocity(isample) + &
               t * (thisStream%velocity(iSample+1)-thisStream%velocity(iSample))
       endif
    endif

    if (present(inflow)) then
       inflow = .not.outsideStream
    endif
  end subroutine getMagStreamValues2

  subroutine getMagStreamValues3(thisOctal, subcell, stream, grid)
    type(GRIDTYPE) :: grid
    type(OCTAL) :: thisOctal
    integer :: subcell
    type(STREAMTYPE) :: stream
    real(double) :: rho
    real  :: temperature
    type(VECTOR) :: vel
    integer :: n, i

    temperature = 0.
    rho = 0.d0
    vel = VECTOR(0.d0, 0.d0, 0.d0)
    
    n = 0
    do i = 1, stream%nSamples
       if (inSubcell(thisOctal, subcell, stream%position(i))) then
          n = n + 1
          rho = rho + stream%rho(i)
          temperature = temperature + stream%temperature(i)
          vel = vel + stream%velocity(i)
       endif
    enddo
    if (n ==0 ) then
       thisOctal%inFlow(subcell) = .false.
       thisOctal%rho(subcell) = 1.d-25
       thisOctal%temperature(subcell) = 6000.
       thisOctal%velocity(subcell) = VECTOR(0.,0.,0.)
    else
       thisOctal%inFlow(subcell) = .true.
       thisOctal%rho(subcell) = rho / dble(n)
       thisOctal%temperature(subcell) = temperature/real(n)
       thisOctal%velocity(subcell) = vel / dble(n)
       thisOctal%microturb(subcell) = 50.d5/cspeed!!!!!!!!!!!!!!!!!!!!

       if (subcell == thisOctal%maxchildren) then
          call fillVelocityCorners(thisOctal, grid, magstreamvelocity, .true.)
!          write(*,*) " "
!          write(*,*) thisOctal%velocity(subcell)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(14)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(15)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(17)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(18)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(23)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(24)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(26)*real(cspeed/1.e5)
!          write(*,*) thisOctal%cornervelocity(27)*real(cspeed/1.e5)
       endif
    endif

  end subroutine getMagStreamValues3

  type(VECTOR)  function magStreamVelocity(point, grid) result(vel)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: point
    real(double) :: minDist, dist, dist1, dist2, t
    integer :: i, j, isample, istream
    minDist = 1.d30

    do iStream = 1, globalnStream
       do iSample = 2, globalStream(iStream)%nSamples
          dist = distancePointLineSegment(globalStream(iStream)%position(iSample-1), &
               globalStream(iStream)%position(iSample), point)
          if (dist < minDist) then
             minDist = dist
             i = iStream
             j = iSample-1
          endif
       enddo
    enddo
    dist1 = modulus(point - globalStream(i)%position(j))
    dist2 = modulus(point - globalStream(i)%position(j+1))
    t = (minDist - dist1)/ (dist2 - dist1)
    vel = globalStream(i)%velocity(j) + t * (globalStream(i)%velocity(j+1) - globalStream(i)%velocity(j))
!    write(*,*) modulus(vel)*cspeed/1.d5
  end function magStreamVelocity


  subroutine findNearestSample(thisStream, position, iSample, t, outsideStream)
    type(STREAMTYPE) :: thisStream
    type(VECTOR) :: rVec, position
    real(double) :: minDist, dist
    integer :: iSample
    real(double) :: t
    integer :: i
    logical :: outsideStream

    outsideStream = .true.
    minDist = 1.d30
    iSample = 0

    rVec = thisStream%position(thisStream%nSamples) - position
    if ((rVec .dot. thisStream%direction(thisStream%nSamples)) < 0.d0) then
       outsideStream = .true.
       iSample = thisStream%nSamples
       t = 0.d0
       goto 666
    endif


    if (modulus(position) > modulus(thisStream%position(thisStream%nSamples))) then
       iSample = thisStream%nSamples
       t = 0.d0
       goto 666
    endif
    
    if (modulus(position) < modulus(thisStream%position(1))) then
       iSample = 1
       t = 0.d0
       goto 666
    endif

    do i = 1, thisStream%nSamples
       rVec = thisStream%position(i) - position
       dist = modulus(rVec)
       if (dist < thisStream%streamRadius(i)) outsideStream = .false.
       if ( (rVec .dot. thisStream%direction(i)) >= 0.d0 ) then
          if (dist < minDist) then
             minDist = dist
             iSample = i
          endif
       endif

    enddo
    if ((iSample == 0).and.(.not.outsideStream)) then
       write(*,*) "error in find nearest sample"
       write(*,*) "position",position
       write(*,*) modulus(position)/modulus(thisStream%position(1))
       iSample = 1
    endif
    if (outSideStream) iSample = 1
    rVec = thisStream%position(iSample) - position
    if (modulus(thisStream%position(iSample+1)-thisStream%position(iSample)) /= 0.d0) then
       t = (rVec.dot.thisStream%direction(iSample)) / &
            modulus(thisStream%position(iSample+1)-thisStream%position(iSample))
    else
       t = 0.d0
    endif
666 continue

  end  subroutine findNearestSample


  function  getVel(grid, thisOctal, subcell, posVec, direction) result (outVec)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer subcell
    real(double) :: distToNextCell, r, x1
    type(VECTOR) :: direction, rVec, sVec, posVec, thisDirection
    type(VECTOR) :: v1,  outVec
    integer :: neighbourSubcell

    rVec = subcellCentre(thisOctal, subcell)
    r = ((posVec - rVec).dot.direction)
    if (r  > 0.d0 ) then

       thisDirection = direction
       rVec = subcellCentre(thisOctal, subcell)
       call distanceToCellBoundary(grid, rVec, thisdirection, DisttoNextCell)
       v1 = thisOctal%velocity(subcell)
       x1 = distToNextCell
       distToNextCell = distToNextcell + 0.01d0*grid%halfSmallestSubcell
       sVec = rVec + distToNextCell * thisdirection
       if (inOctal(grid%octreeRoot, sVec)) then
          neighbourOctal => thisOctal
          CALL findSubcellLocal(sVec,neighbourOctal,neighboursubcell)
          if (neighbourOctal%inFlow(neighbourSubcell)) then
             v1 = neighbourOctal%velocity(neighbourSubcell)
             x1 = disttoNextcell + neighbourOctal%subcellSize/2.d0
          endif
       endif
       outVec = thisOctal%velocity(subcell) + (r/x1) * (v1 - thisOctal%velocity(subcell))
    else
       thisDirection = (-1.d0)*direction
       rVec = subcellCentre(thisOctal, subcell)
       call distanceToCellBoundary(grid, rVec, thisdirection, DisttoNextCell)
       v1 = thisOctal%velocity(subcell)
       x1 = distToNextCell
       distToNextCell = distToNextcell + 0.01d0*grid%halfSmallestSubcell
       sVec = rVec + distToNextCell * thisdirection
       if (inOctal(grid%octreeRoot, sVec)) then
          neighbourOctal => thisOctal
          CALL findSubcellLocal(sVec,neighbourOctal,neighboursubcell)
          if (neighbourOctal%inFlow(neighbourSubcell)) then
             v1 = neighbourOctal%velocity(neighbourSubcell)
             x1 = disttoNextcell + neighbourOctal%subcellSize/2.d0
          endif
       endif
       outVec = thisOctal%velocity(subcell) + (abs(r)/x1) * (v1 - thisOctal%velocity(subcell))
    endif

  end function getVel


  subroutine genericAccretionSurface(surface, grid, lineFreq,coreContFlux,fAccretion)

    USE surface_mod, only: createProbs, sumSurface

    type(SURFACETYPE) :: surface
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    type(VECTOR) :: rVec
    integer :: subcell
    real(double) :: v, area, T, flux, power, totalArea, accretingArea, mdot, totalMdot
    integer :: i
    REAL(double), INTENT(IN) :: coreContFlux
    REAL, INTENT(IN) :: lineFreq
    REAL, INTENT(OUT) :: fAccretion ! erg s^-1 Hz^-1

    write(*,*) "calculating generic accretion surface"
    accretingArea = 0.d0
    totalArea = 0.d0
    totalmdot = 0.d0

    do i = 1, surface%nElements
       rVec = surface%element(i)%position + grid%halfSmallestSubcell * surface%element(i)%norm
       CALL findSubcellTD(rVec,grid%octreeRoot,thisOctal,subcell)


       v = modulus(thisOctal%velocity(subcell))*cspeed
       area = (surface%element(i)%area*1.d20)
       mdot = thisOctal%rho(subcell) * v * area
       totalMdot = totalMdot + mdot
       power = 0.5d0 * mdot * v**2
       if (area /= 0.d0) then
          flux = power / area
       else
          flux = 0.d0
       endif


       T = (flux/stefanBoltz)**0.25d0

       totalArea = totalArea + area
       
       if (T > 4000.0d0) then 
          ! in futire this should be somehow compared to the effective
          ! temperature of photosphere. 
          !       if (T > T_eff) then
          surface%element(i)%hot = .true.
          allocate(surface%element(i)%hotFlux(surface%nNuHotFlux))
          
          surface%element(i)%hotFlux(:) = &
               pi*blackbody(REAL(T), 1.e8*REAL(cSpeed)/surface%nuArray(:))
          surface%element(i)%temperature = T
          accretingArea = accretingArea + area
       else 
          surface%element(i)%hot = .false.
       end if
    enddo


    write(*,*) "Spot fraction is: ",100.d0*accretingArea/totalArea, "%"
    write(*,*) "Mass accretion rate is: ",(totalMdot/mSol)*(365.25d0*24.d0*3600.d0), "solar masses/year"
    CALL createProbs(surface,lineFreq,coreContFlux,fAccretion)
    CALL sumSurface(surface)

  end subroutine genericAccretionSurface


  subroutine allocateOctalAttributes(grid, thisOctal)
    use input_variables, only : mie, cmf, nAtom, nDustType, molecular, TminGlobal, &
         photoionization, hydrodynamics, sobolev, storeScattered, h21cm
    use gridtype_mod, only: statEqMaxLevels
    type(OCTAL), pointer :: thisOctal
    type(GRIDTYPE) :: grid
    integer, parameter :: nTheta = 10 , nphi = 10

    if ( h21cm ) then 
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
    end if

    if (mie) then
       call allocateAttribute(thisOctal%oldFrac, thisOctal%maxChildren)
       thisOctal%oldFrac = 1.d-30
       call allocateAttribute(thisOctal%dustType, thisOctal%maxChildren)
       thisOctal%dustType = 1
       ALLOCATE(thisOctal%dusttypefraction(thisOctal%maxchildren,  nDustType))
       thisOctal%dustTypeFraction = 0.d0
       thisOctal%dustTypeFraction(:,:) = 1.d0
       thisOctal%inflow = .true.

       if (storescattered) allocate(thisOctal%scatteredIntensity(thisOctal%maxChildren, ntheta, nPhi))

       call allocateAttribute(thisOctal%diffusionApprox, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nDiffusion, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%eDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%diffusionCoeff, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldeDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nDirectPhotons, thisOctal%maxChildren)
       
       call allocateAttribute(thisOctal%underSampled, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldTemperature, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%kappaRoss, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%distanceGrid, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nCrossings, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaCont, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasCont3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasLine3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%probDistLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%probDistCont, thisOctal%maxChildren)


    endif

    if (cmf) then
       ALLOCATE(thisOctal%N(8,grid%maxLevels))

       allocate(thisOctal%atomAbundance(8, 1:nAtom))
       thisOctal%atomAbundance(:, 1) = 0.71d0 / mHydrogen
       if (nAtom > 1) then
          thisOctal%atomAbundance(:, 2:nAtom) =  0.27d0 / (4.d0*mHydrogen) !assume higher atoms are helium
       endif
       call allocateAttribute(thisOctal%microturb, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%ne, thisOctal%maxChildren)

    endif

    if (sobolev) then
       call allocateAttribute(thisOctal%biasCont3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasLine3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaCont, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%ne, thisOctal%maxChildren)
       ALLOCATE(thisOctal%N(8,1:stateqMaxLevels))
       ALLOCATE(thisOctal%kappaAbs(8,1))
    endif

    if (molecular) then
       call allocateAttribute(thisOctal%molAbundance, thisOctal%maxChildren)
       thisOctal%molAbundance(:) = 1.e-30
       call allocateAttribute(thisOctal%temperatureGas, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%temperatureDust, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nh2, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%microturb, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%molmicroturb, thisOctal%maxChildren)
    endif

    thisOctal%rho = amr_min_rho
    thisOctal%gasOpacity = .false.
    thisOctal%temperature = TMinGlobal


    if (photoionization) then
       call allocateAttribute(thisOctal%oldFrac, thisOctal%maxChildren)
       thisOctal%oldFrac = 1.d-30
       call allocateAttribute(thisOctal%dustType, thisOctal%maxChildren)
       thisOctal%dustType = 1
       ALLOCATE(thisOctal%dusttypefraction(thisOctal%maxchildren,  nDustType))
       thisOctal%dustTypeFraction = 0.d0
       thisOctal%dustTypeFraction(:,:) = 1.d0

       call allocateAttribute(thisOctal%diffusionApprox, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nDiffusion, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%eDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%diffusionCoeff, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldeDens, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nDirectPhotons, thisOctal%maxChildren)
       
       call allocateAttribute(thisOctal%underSampled, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%oldTemperature, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%kappaRoss, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%distanceGrid, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nCrossings, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nTot, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%chiLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasLine3D, thisOctal%maxChildren)

       call allocateAttribute(thisOctal%etaCont, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nh, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%ne, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nhi, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nhei, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%nhii, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%biasCont3D, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%etaLine, thisOctal%maxChildren)

       call allocateAttribute(thisOctal%HHeating, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%HeHeating, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%probDistLine, thisOctal%maxChildren)
       call allocateAttribute(thisOctal%probDistCont, thisOctal%maxChildren)

       allocate(thisOctal%ionFrac(1:thisOctal%maxchildren, 1:grid%nIon))
       allocate(thisOctal%photoionCoeff(1:thisOctal%maxchildren, 1:grid%nIon))
    endif

    if (associated(thisOctal%ionFrac)) thisOctal%ionFrac = 1.e-30
    
    if (associated(thisOctal%photoIonCoeff)) then
       thisOctal%photoIonCoeff = 0.d0
    endif
       
    if (hydrodynamics) then



       call allocateAttribute(thisOctal%q_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_minus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%q_i_minus_2,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%x_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%x_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%x_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%u_interface,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%u_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%u_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%flux_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%flux_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%flux_i_minus_1,thisOctal%maxchildren)


       call allocateAttribute(thisOctal%phiLimit,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%ghostCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%feederCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%edgeCell,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%refinedLastTime,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%pressure_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%pressure_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%pressure_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rhou,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhov,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rhow,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rhoe,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%energy,thisOctal%maxchildren)


       call allocateAttribute(thisOctal%phi_i,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%phi_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%rho_i_plus_1,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%rho_i_minus_1,thisOctal%maxchildren)

       call allocateAttribute(thisOctal%boundaryCondition,thisOctal%maxchildren)
       call allocateAttribute(thisOctal%boundaryPartner,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%changed,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%rLimit,thisOctal%maxChildren)

       call allocateAttribute(thisOctal%iEquationOfState,thisOctal%maxChildren)
       call allocateAttribute(thisOctal%gamma,thisOctal%maxChildren)

    endif
  end  subroutine allocateOctalAttributes


    subroutine deallocateOctalDynamicAttributes(thisOctal)
      type(OCTAL):: thisOctal

       call deallocateAttribute(thisOctal%oldFrac)
       call deallocateAttribute(thisOctal%dustType)
       call deallocateAttribute(thisOctal%dusttypefraction)
       call deallocateAttribute(thisOctal%diffusionApprox)
       call deallocateAttribute(thisOctal%scatteredIntensity)
       call deallocateAttribute(thisOctal%nDiffusion)
       call deallocateAttribute(thisOctal%eDens)
       call deallocateAttribute(thisOctal%diffusionCoeff)
       call deallocateAttribute(thisOctal%oldeDens)
       call deallocateAttribute(thisOctal%nDirectPhotons)
       
       call deallocateAttribute(thisOctal%underSampled)
       call deallocateAttribute(thisOctal%oldTemperature)
       call deallocateAttribute(thisOctal%kappaRoss)
       call deallocateAttribute(thisOctal%distanceGrid)
       call deallocateAttribute(thisOctal%nCrossings)
       call deallocateAttribute(thisOctal%nTot)
       call deallocateAttribute(thisOctal%etaCont)
       call deallocateAttribute(thisOctal%etaLine)
       call deallocateAttribute(thisOctal%chiLine)
       call deallocateAttribute(thisOctal%biasCont3D)
       call deallocateAttribute(thisOctal%biasLine3D)

       call deallocateAttribute(thisOctal%probDistLine)
       call deallocateAttribute(thisOctal%probDistCont)



       call deAllocateAttribute(thisOctal%N)
       call deAllocateAttribute(thisOctal%atomAbundance)

       call deAllocateAttribute(thisOctal%molAbundance)
       call deAllocateAttribute(thisOctal%temperatureGas)
       call deAllocateAttribute(thisOctal%temperatureDust)
       call deAllocateAttribute(thisOctal%nh2)
       call deAllocateAttribute(thisOctal%microturb)
       call deAllocateAttribute(thisOctal%molmicroturb)

       call deAllocateAttribute(thisOctal%ionFrac)
       call deAllocateAttribute(thisOctal%PhotoIonCoeff)
       call deAllocateAttribute(thisOctal%scatteredIntensity)
       

       call deAllocateAttribute(thisOctal%q_i)
       call deAllocateAttribute(thisOctal%q_i_plus_1)
       call deAllocateAttribute(thisOctal%q_i_minus_1)
       call deAllocateAttribute(thisOctal%q_i_minus_2)

       call deAllocateAttribute(thisOctal%x_i)
       call deAllocateAttribute(thisOctal%x_i_plus_1)
       call deAllocateAttribute(thisOctal%x_i_minus_1)

       call deAllocateAttribute(thisOctal%u_interface)
       call deAllocateAttribute(thisOctal%u_i_plus_1)
       call deAllocateAttribute(thisOctal%u_i_minus_1)

       call deAllocateAttribute(thisOctal%flux_i)
       call deAllocateAttribute(thisOctal%flux_i_plus_1)
       call deAllocateAttribute(thisOctal%flux_i_minus_1)

       call deAllocateAttribute(thisOctal%phiLimit)

       call deAllocateAttribute(thisOctal%ghostCell)
       call deAllocateAttribute(thisOctal%feederCell)
       call deAllocateAttribute(thisOctal%edgeCell)
       call deAllocateAttribute(thisOctal%refinedLastTime)

       call deAllocateAttribute(thisOctal%pressure_i)
       call deAllocateAttribute(thisOctal%pressure_i_plus_1)
       call deAllocateAttribute(thisOctal%pressure_i_minus_1)

       call deAllocateAttribute(thisOctal%rhou)
       call deAllocateAttribute(thisOctal%rhov)
       call deAllocateAttribute(thisOctal%rhow)

       call deAllocateAttribute(thisOctal%rhoe)
       call deAllocateAttribute(thisOctal%energy)

       call deAllocateAttribute(thisOctal%iEquationOfState)
       call deAllocateAttribute(thisOctal%gamma)


       call deAllocateAttribute(thisOctal%phi_i)
       call deAllocateAttribute(thisOctal%phi_i_plus_1)
       call deAllocateAttribute(thisOctal%phi_i_minus_1)

       call deAllocateAttribute(thisOctal%rho_i_plus_1)
       call deAllocateAttribute(thisOctal%rho_i_minus_1)

       call deAllocateAttribute(thisOctal%boundaryCondition)
       call deAllocateAttribute(thisOctal%boundaryPartner)
       call deAllocateAttribute(thisOctal%changed)
       call deAllocateAttribute(thisOctal%rLimit)

       call deallocateAttribute(thisOctal%neighbourOctal)
       call deallocateAttribute(thisOctal%neighbourSubcell)

       call deallocateAttribute(thisOctal%mpiBoundaryStorage)
     end subroutine deallocateOctalDynamicAttributes


  RECURSIVE SUBROUTINE cleanupAMRgrid(thisOctal)
      
        TYPE(OCTAL), POINTER  :: thisOctal 
        TYPE(OCTAL), POINTER  :: child
        INTEGER :: i
        
        if (thisOctal%nChildren == thisOctal%maxChildren) then
           call deallocateOctalDynamicAttributes(thisOctal)
        endif

        IF ( thisOctal%nChildren > 0 ) THEN
           ! call this subroutine recursively on each of its children
           DO i = 1, thisOctal%nChildren, 1
              child => thisOctal%child(i)
              CALL cleanupAMRgrid(child)
           END DO
        END IF
        
  END SUBROUTINE cleanupAMRgrid


  recursive subroutine setupNeighbourPointers(grid, thisOctal)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer :: neighbourOctal
    integer :: neighbourSubcell
    type(octal), pointer  :: child 
    type(VECTOR) :: centre, rvec
    real(double) :: d
    integer :: subcell, i
    type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.d0)
    type(VECTOR), parameter :: yAxis = VECTOR(0.d0, 1.d0, 0.d0)
    type(VECTOR), parameter :: xAxis = VECTOR(1.d0, 0.d0, 0.d0)
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  setupNeighbourPointers(grid, child)
                exit
             end if
          end do
       else

          if (.not.associated(thisOctal%neighbourOctal)) then
             allocate(thisOctal%neighbourOctal(thisOctal%maxChildren, 6, 4))
          endif
          if (.not.associated(thisOctal%neighbourSubcell)) then
             allocate(thisOctal%neighbourSubcell(thisOctal%maxChildren, 6, 4))
             thisOctal%neighbourSubcell = -666
          endif

          centre = subcellCentre(thisOctal, subcell)
          d = thisOctal%subcellSize/2.d0

          ! top face (plus z)

          rVec = centre + zAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) - yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 1, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 1, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 1, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 1, 1) = -1
          endif

          if (inOctal(grid%octreeRoot, rVec)) then
             rVec = centre + zAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) - yAxis*(d/2.d0)
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 1, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 1, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 1, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 1, 2) = -1
          endif


          rVec = centre + zAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) + yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 1, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 1, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 1, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 1, 3) = -1
          endif

          rVec = centre + zAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) + yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 1, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 1, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 1, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 1, 4) = -1
          endif

          ! bottom face (minus z)

          rVec = centre - zAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) - yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 2, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 2, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 2, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 2, 1) = -1
          endif

          rVec = centre - zAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) - yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 2, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 2, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 2, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 2, 2) = -1
          endif

          rVec = centre - zAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) + yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 2, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 2, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 2, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 2, 3) = -1
          endif

          rVec = centre - zAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) + yAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 2, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 2, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 2, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 2, 4) = -1
          endif


          ! left face (minus x)

          rVec = centre - xAxis * (d+grid%halfSmallestSubcell*00.1d0) - yAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 3, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 3, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 3, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 3, 1) = -1
          endif

          rVec = centre - xAxis * (d+grid%halfSmallestSubcell*00.1d0) + yAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 3, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 3, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 3, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 3, 2) = -1
          endif

          rVec = centre - xAxis * (d+grid%halfSmallestSubcell*00.1d0) - yAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 3, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 3, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 3, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 3, 3) = -1
          endif

          rVec = centre - xAxis * (d+grid%halfSmallestSubcell*00.1d0) + yAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 3, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 3, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 3, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 3, 4) = -1
          endif

          ! right face (plus x)

          rVec = centre + xAxis * (d+grid%halfSmallestSubcell*00.1d0) - yAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 4, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 4, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 4, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 4, 1) = -1
          endif

          rVec = centre + xAxis * (d+grid%halfSmallestSubcell*00.1d0) + yAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 4, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 4, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 4, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 4, 2) = -1
          endif

          rVec = centre + xAxis * (d+grid%halfSmallestSubcell*00.1d0) - yAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 4, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 4, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 4, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 4, 3) = -1
          endif

          rVec = centre + xAxis * (d+grid%halfSmallestSubcell*00.1d0) + yAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 4, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 4, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 4, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 4, 4) = -1
          endif


          ! front face (plus y)

          rVec = centre + yAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 5, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 5, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 5, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 5, 1) = -1
          endif

          rVec = centre + yAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 5, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 5, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 5, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 5, 2) = -1
          endif

          rVec = centre + yAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 5, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 5, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 5, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 5, 3) = -1
          endif

          rVec = centre + yAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 5, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 5, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 5, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 5, 4) = -1
          endif

          ! back face (minus y)

          rVec = centre - yAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 6, 1)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 6, 1) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 6, 1)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 6, 1) = -1
          endif

          rVec = centre - yAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) - zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 6, 2)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 6, 2) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 6, 2)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 6, 2) = -1
          endif

          rVec = centre - yAxis * (d+grid%halfSmallestSubcell*00.1d0) - xAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 6, 3)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 6, 3) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 6, 3)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 6, 3) = -1
          endif

          rVec = centre - yAxis * (d+grid%halfSmallestSubcell*00.1d0) + xAxis*(d/2.d0) + zAxis*(d/2.d0)
          if (inOctal(grid%octreeRoot, rVec)) then
             neighbourOctal => thisOctal
             neighbourSubcell = subcell
             call findSubcellLocal(rVec, neighbourOctal, neighbourSubcell)
             thisOctal%neighbourOctal(subcell, 6, 4)%pointer => neighbourOctal
             thisOctal%neighbourSubcell(subcell, 6, 4) = neighbourSubcell
          else
             thisOctal%neighbourOctal(subcell, 6, 4)%pointer => null()
             thisOctal%neighbourSubcell(subcell, 6, 4) = -1
          endif
!          write(*,*) "subcell ",subcell
!          write(*,*) "top ",thisOctal%neighbourSubcell(subcell, 1, 1:4)
!          write(*,*) "bot ",thisOctal%neighbourSubcell(subcell, 2, 1:4)
!          write(*,*) "left ",thisOctal%neighbourSubcell(subcell, 3, 1:4)
!          write(*,*) "right ",thisOctal%neighbourSubcell(subcell, 4, 1:4)
!          write(*,*) "front ",thisOctal%neighbourSubcell(subcell, 5, 1:4)
!          write(*,*) "back ",thisOctal%neighbourSubcell(subcell, 6, 1:4)
       endif
    enddo
  end subroutine setupNeighbourPointers

  subroutine getNeighbourFromPointOnFace(rVec, uHat, thisOctal, subcell, neighbourOctal, neighbourSubcell)
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: neighbourSubcell, subcell
    type(VECTOR) :: rVec, centre, normVec, uHat
    logical :: xPos, xNeg, yPos, yNeg, zPos, zNeg

    centre = subcellCentre(thisOctal, subcell)
    normVec = 2.d0*((rVec - centre)/thisOctal%subcellSize)

!    write(*,*) "normVec ",normvec
    xPos = normVec%x >= 0.d0
    xNeg = .not.xPos

    yPos = normVec%y >= 0.d0
    yNeg = .not.yPos

    zPos = normVec%z >= 0.d0
    zNeg = .not.zPos


!    write(*,*) "xpos,xneg ",xpos,xneg
!    write(*,*) "ypos,yneg ",ypos,yneg
!    write(*,*) "zpos,zneg ",zpos,zneg

!    write(*,*) "current subcell", subcell
    if ((normVec%z > 0.9999999d0).and.(uHat%z > 0.d0)) then ! top face
       if ((xNeg.and.yNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 1, 1)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 1, 1)
       else if ((xPos.and.yNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 1, 2)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 1, 2)
       else if ((xNeg.and.yPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 1, 3)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 1, 3)
       else if ((xPos.and.yPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 1, 4)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 1, 4)
       endif
    else if ((normVec%z < -0.9999999d0).and.(uHat%z < 0.d0)) then ! bottom face
       if ((xNeg.and.yNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 2, 1)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 2, 1)
       else if ((xPos.and.yNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 2, 2)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 2, 2)
       else if ((xNeg.and.yPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 2, 3)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 2, 3)
       else if ((xPos.and.yPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 2, 4)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 2, 4)
       endif
    else if ((normVec%x < -0.9999999d0).and.(uHat%x < 0.d0)) then ! left face
       if ((yNeg.and.zNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 3, 1)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 3, 1)
       else if ((yPos.and.zNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 3, 2)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 3, 2)
       else if ((yNeg.and.zPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 3, 3)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 3, 3)
       else if ((yPos.and.zPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 3, 4)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 3, 4)
       endif
    else if ((normVec%x > 0.9999999d0).and.(uHat%x > 0.d0)) then ! right face
       if ((yNeg.and.zNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 4, 1)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 4, 1)
       else if ((yPos.and.zNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 4, 2)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 4, 2)
       else if ((yNeg.and.zPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 4, 3)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 4, 3)
       else if ((yPos.and.zPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 4, 4)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 4, 4)
       endif
    else if ((normVec%y > 0.9999999d0).and.(uHat%y > 0.d0)) then ! front face
       if ((xNeg.and.zNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 5, 1)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 5, 1)
       else if ((xPos.and.zNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 5, 2)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 5, 2)
       else if ((xNeg.and.zPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 5, 3)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 5, 3)
       else if ((xPos.and.zPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 5, 4)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 5, 4)
       endif
    else if ((normVec%y < -0.9999999d0).and.(uHat%y < 0.d0)) then ! back
       if ((xNeg.and.zNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 6, 1)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 6, 1)
       else if ((xPos.and.zNeg)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 6, 2)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 6, 2)
       else if ((xNeg.and.zPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 6, 3)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 6, 3)
       else if ((xPos.and.zPos)) then
          neighbourOctal => thisOctal%neighbourOctal(subcell, 6, 4)%pointer
          neighbourSubcell = thisOctal%neighbourSubcell(subcell, 6, 4)
       endif
    endif
  end subroutine getNeighbourFromPointOnFace


  subroutine howmanysplits()
    implicit none
    character(len = 100) :: message    

    if(mass_split .ne. 0) then
       write(message, *) "There were ", mass_split, " splits by mass"
       call writeinfo(message, TRIVIAL)
    endif

    if(mass_split2 .ne. 0) then
       write(message, *) "There were ", mass_split2, " splits by mass2"
       call writeinfo(message, TRIVIAL)
    endif
    
    if(density_split .ne. 0) then
       write(message, *) "There were ", density_split, " splits by density"
       call writeinfo(message, TRIVIAL)
    endif

    if(both_split .ne. 0) then
       write(message, *) "There were ", both_split, " splits by both"
       call writeinfo(message, TRIVIAL)
    endif

    if(maxdensity_split .ne. 0) then
       write(message, *) "There were ", maxdensity_split, " splits by maximum density"
       call writeinfo(message, TRIVIAL)
    endif

    if(velocity_split .ne. 0) then
       write(message, *) "There were ", velocity_split, " splits by velocity"
       call writeinfo(message, TRIVIAL)
    endif

    if(scaleheighta_count .ne. 0) then
       write(message, *) "There were ", scaleheighta_count, " splits by scaleheight condition A"
       call writeinfo(message, TRIVIAL)
    endif

    if(scaleheightb_count .ne. 0) then
       write(message, *) "There were ", scaleheightb_count, " splits by scaleheight condition B"
       call writeinfo(message, TRIVIAL)
    endif

    if(scaleheightc_count .ne. 0) then
       write(message, *) "There were ", scaleheightc_count, " splits by scaleheight condition C"
       call writeinfo(message, TRIVIAL)
    endif

  end subroutine howmanysplits

  function returnScatteredIntensity(position, thisOctal, subcell, uHat) result(intensity)
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(VECTOR) :: uHat, thisVec, position
    real(double) :: intensity, thisTheta, thisPhi, ang
    integer :: nTheta, nPhi
    integer :: iTheta, iPhi
    
    nTheta = SIZE(thisOctal%scatteredIntensity,2)
    nPhi = SIZE(thisOctal%scatteredIntensity,3)

    ang = atan2(position%y, position%x)

    thisVec = uHat
    
    if (thisOctal%twoD) then
       thisVec = rotateZ(uHat, -ang)
    endif

    thisTheta = acos(thisvec%z)
    thisPhi = atan2(thisVec%y,thisVec%x)

    if (thisPhi < 0.d0) thisPhi = thisPhi + twoPi 
    iTheta = nint((thisTheta / pi) * dble(nTheta-1))+1
    iphi = nint((thisPhi / twoPi) * dble(nPhi-1))+1

    intensity = thisOctal%scatteredIntensity(subcell,iTheta, iPhi) 
!    intensity = SUM(thisOctal%scatteredIntensity(subcell,:,:))/100.d0
  end function returnScatteredIntensity
  
  subroutine quadint(t1,t2,t3,weights)

    real(double) :: t1, t2, t3
    real(double), save :: t1old, t2old, t3old
    real(double) :: xu,xv,xw,yu,yv,yw,zu,zv,zw
    real(double) :: xbasis1, xbasis2, ybasis1, ybasis2, zbasis1, zbasis2
    real(double) :: xyweights(9)
    real(double) :: weights(27)
    real(double), save :: oldweights(27)
      
    if(t1 .eq. t1old) then
       if(t2 .eq. t2old) then
          if(t3 .eq. t3old) then
             weights = oldweights
             return
          endif
       endif
    endif

    xbasis1 = t1 - 1.d0
    xbasis2 = 2.d0 * t1 - 1.d0
    ybasis1 = t2 - 1.d0
    ybasis2 = 2.d0 * t2 - 1.d0
    zbasis1 = t3 - 1.d0
    zbasis2 = 2.d0 * t3 - 1.d0
    
    xu = xbasis1 * xbasis2
    xv = -4.d0 * xbasis1 * t1
    xw = xbasis2 * t1
            
    yu = ybasis1 * ybasis2
    yv = -4.d0 * ybasis1 * t2
    yw = ybasis2 * t2

    zu = zbasis1 * zbasis2
    zv = -4.d0 * zbasis1 * t3
    zw = zbasis2 * t3

    ! base levels
            
    xyweights(:) = (/ xu*yu, xv*yu, xw*yu, xu*yv, xv*yv, xw*yv, xu*yw, xv*yw, xw*yw /) 
            
    weights(1:9) = zu * xyweights
    weights(10:18) = zv * xyweights
    weights(19:27) = zw * xyweights
    
    oldweights = weights
    t1old = t1
    t2old = t2
    t3old = t3

  end subroutine quadint

  real(double) function amrgriddensity(position, grid, linearinterp) result(rho)

    type(GRIDTYPE) :: grid
    type(OCTAL),save, pointer :: resultoctal
    integer :: subcell
    type(VECTOR) :: position, centre, point_local
    real(double) :: t1, t2, t3
    real(double) :: inc, fac
    real(double) :: weights(27)
    real(double) :: rhoa(1)

    logical, save :: firsttime = .true.
    logical, optional :: linearinterp
    logical :: linear

    if(present(linearinterp)) then
       linear = linearinterp
    else
       linear = .true.
    endif

    if(firsttime) then
       CALL findSubcellTD(position,grid%octreeroot, resultOctal,subcell)
       firsttime = .false.
    endif

 !   if(.not. inoctal(resultoctal,position)) then
!       CALL findSubcellTD(position,grid%octreeroot, resultOctal,subcell)
       CALL findSubcellLocal(position,resultOctal,subcell)
!    endif

    if (resultoctal%threeD) then
       point_local = position
    elseif (resultoctal%twoD) then
       point_local = projectToXZ(position)
    else !oneD
       point_local = VECTOR(modulus(position), 0.d0, 0.d0)
    end if

    inc = resultOctal%subcellSize
    centre = resultoctal%centre
    fac = 0.5d0 / resultOctal%subcellsize

    t1 = fac * (point_local%x - (centre%x - inc))
    t2 = fac * (point_local%y - (centre%y - inc))
    t3 = fac * (point_local%z - (centre%z - inc))

    if(resultoctal%oneD) then
       centre = VECTOR(modulus(resultoctal%centre), 0.d0, 0.d0)
       t1 = fac * (point_local%x - (centre%x - inc))

       rho = (2.d0 * t1**2 - 3.d0 * t1 + 1.d0) * resultOctal%cornerrho(1) + &
            (-4.d0 * t1 * (t1 - 1.d0))         * resultOctal%cornerrho(2) + &
            (t1 * (2.d0 * t1 - 1.d0))          * resultOctal%cornerrho(3)

       return
    endif
    

    if(linear) then

       inc = 0.5 * resultOctal%subcellSize
       centre = subcellCentre(resultOctal,subcell)
       fac = 1. / resultOctal%subcellsize
      
       t1 = MAX(0.0_oc, fac * (point_local%x - (centre%x - inc)))
       t2 = MAX(0.0_oc, fac * (point_local%y - (centre%y - inc)))
       t3 = MAX(0.0_oc, fac * (point_local%z - (centre%z - inc)))

       SELECT CASE(subcell)
               
       CASE(1)
          rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 1) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 2) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 4) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 5) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(10) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(11) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(13) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(14)
               
            CASE(2)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 2) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 3) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 5) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 6) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(11) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(12) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(14) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(15)
               
            CASE(3)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 4) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 5) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 7) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 8) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(13) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(14) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(16) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(17)
               
            CASE(4)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 5) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho( 6) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 8) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho( 9) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(14) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(15) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(17) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(18)
               
            CASE(5)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(10) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(11) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(13) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(14) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(19) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(20) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(22) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(23)
               
            CASE(6)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(11) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(12) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(14) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(15) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(20) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(21) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(23) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(24)
            
            CASE(7)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(13) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(14) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(16) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(17) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(22) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(23) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(25) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(26)
               
            CASE(8)
               rho = &
                    ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(14) + &
                    ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerrho(15) + &
                    ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(17) + &
                    ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerrho(18) + &
                    ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(23) + &
                    ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerrho(24) + &
                    ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerrho(26) + &
                    ((     t1) * (     t2) * (     t3)) * resultOctal%cornerrho(27)
            CASE DEFAULT
               PRINT *, 'Invalid subcell in amrGridDensity'
               
            end SELECT
            
         else
            call quadint(t1,t2,t3,weights)
            rho = sum(weights(:) * resultoctal%cornerrho(:))
            if(rho .lt. 0.d0) then
               rhoa = resultoctal%cornerrho(maxloc(weights))
               rho = rhoa(1)
            endif
         endif
    
       end function amrgriddensity

  SUBROUTINE fillDensityCorners(thisOctal,grid,densityFunc, velocityfunc,threed)
    ! store the density values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
  
    TYPE(octal), INTENT(INOUT) :: thisOctal
    TYPE(gridtype), INTENT(IN) :: grid
    LOGICAL, INTENT(IN) :: threed
    real(oct)      :: r1, r2, r3
    real(oct)      :: phi1, phi2, phi3
    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    
    INTERFACE 
      real(double) FUNCTION densityFunc(point,grid)
        USE vector_mod
        USE gridtype_mod
        TYPE(vector), INTENT(IN) :: point
        TYPE(gridtype), INTENT(IN)    :: grid
      END FUNCTION densityFunc
    END INTERFACE

    INTERFACE
      type(VECTOR) FUNCTION velocityFunc(point,grid)
        USE vector_mod
        USE gridtype_mod
        TYPE(vector), INTENT(IN) :: point
        TYPE(gridtype), INTENT(IN)    :: grid
      END FUNCTION velocityFunc
    END INTERFACE

    if (thisOctal%oneD) then
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       y1 = 0.d0
       z1 = 0.d0
       thisOctal%cornerrho(1) = densityFunc(vector(x1,y1,z1),grid)
       thisOctal%cornerrho(2) = densityFunc(vector(x2,y1,z1),grid)
       thisOctal%cornerrho(3) = densityFunc(vector(x3,y1,z1),grid)
       goto 667
    endif


    if (thisOctal%threed) then
       if (.not.thisOctal%cylindrical) then ! 3d cartesian case
          ! we first store the values we use to assemble the position vectors
          
          x1 = thisOctal%centre%x - thisOctal%subcellSize
          x2 = thisOctal%centre%x
          x3 = thisOctal%centre%x + thisOctal%subcellSize
          
          y1 = thisOctal%centre%y - thisOctal%subcellSize
          y2 = thisOctal%centre%y
          y3 = thisOctal%centre%y + thisOctal%subcellSize
          
          z1 = thisOctal%centre%z - thisOctal%subcellSize
          z2 = thisOctal%centre%z
          z3 = thisOctal%centre%z + thisOctal%subcellSize
                    
          ! now store the 'base level' values

          thisOctal%cornerrho(14) = densityFunc(vector(x2,y2,z2),grid)          
          thisOctal%cornerVelocity(14) = velocityFunc(vector(x2,y2,z2),grid)          
                              
          thisOctal%cornerrho(1) = densityFunc(vector(x1,y1,z1),grid)
          thisOctal%cornerVelocity(1) = velocityFunc(vector(x1,y1,z1),grid)
          thisOctal%cornerrho(2) = densityFunc(vector(x2,y1,z1),grid)
          thisOctal%cornerVelocity(2) = velocityFunc(vector(x2,y1,z1),grid)
          thisOctal%cornerrho(3) = densityFunc(vector(x3,y1,z1),grid)
          thisOctal%cornerVelocity(3) = velocityFunc(vector(x3,y1,z1),grid)
          thisOctal%cornerrho(4) = densityFunc(vector(x1,y2,z1),grid)
          thisOctal%cornerVelocity(4) = velocityFunc(vector(x1,y2,z1),grid)
          thisOctal%cornerrho(5) = densityFunc(vector(x2,y2,z1),grid)
          thisOctal%cornerVelocity(5) = velocityFunc(vector(x2,y2,z1),grid)
          thisOctal%cornerrho(6) = densityFunc(vector(x3,y2,z1),grid)
          thisOctal%cornerVelocity(6) = velocityFunc(vector(x3,y2,z1),grid)
          thisOctal%cornerrho(7) = densityFunc(vector(x1,y3,z1),grid)
          thisOctal%cornerVelocity(7) = velocityFunc(vector(x1,y3,z1),grid)
          thisOctal%cornerrho(8) = densityFunc(vector(x2,y3,z1),grid)
          thisOctal%cornerVelocity(8) = velocityFunc(vector(x2,y3,z1),grid)
          thisOctal%cornerrho(9) = densityFunc(vector(x3,y3,z1),grid)
          thisOctal%cornerVelocity(9) = velocityFunc(vector(x3,y3,z1),grid)

          ! middle level
          
          thisOctal%cornerrho(10) = densityFunc(vector(x1,y1,z2),grid)
          thisOctal%cornerVelocity(10) = velocityFunc(vector(x1,y1,z2),grid)          
          thisOctal%cornerrho(11) = densityFunc(vector(x2,y1,z2),grid)
          thisOctal%cornerVelocity(11) = velocityFunc(vector(x2,y1,z2),grid)
          thisOctal%cornerrho(12) = densityFunc(vector(x3,y1,z2),grid)
          thisOctal%cornerVelocity(12) = velocityFunc(vector(x3,y1,z2),grid)
          thisOctal%cornerrho(13) = densityFunc(vector(x1,y2,z2),grid)
          thisOctal%cornerVelocity(13) = velocityFunc(vector(x1,y2,z2),grid)

          thisOctal%cornerrho(15) = densityFunc(vector(x3,y2,z2),grid)
          thisOctal%cornerVelocity(15) = velocityFunc(vector(x3,y2,z2),grid)
          thisOctal%cornerrho(16) = densityFunc(vector(x1,y3,z2),grid)
          thisOctal%cornerVelocity(16) = velocityFunc(vector(x1,y3,z2),grid)
          thisOctal%cornerrho(17) = densityFunc(vector(x2,y3,z2),grid)
          thisOctal%cornerVelocity(17) = velocityFunc(vector(x2,y3,z2),grid)
          thisOctal%cornerrho(18) = densityFunc(vector(x3,y3,z2),grid)
          thisOctal%cornerVelocity(18) = velocityFunc(vector(x3,y3,z2),grid)          
          ! top level         

          thisOctal%cornerrho(19) = densityFunc(vector(x1,y1,z3),grid)
          thisOctal%cornerVelocity(19) = velocityFunc(vector(x1,y1,z3),grid)
          thisOctal%cornerrho(20) = densityFunc(vector(x2,y1,z3),grid)
          thisOctal%cornerVelocity(20) = velocityFunc(vector(x2,y1,z3),grid)
          thisOctal%cornerrho(21) = densityFunc(vector(x3,y1,z3),grid)
          thisOctal%cornerVelocity(21) = velocityFunc(vector(x3,y1,z3),grid)
          thisOctal%cornerrho(22) = densityFunc(vector(x1,y2,z3),grid)
          thisOctal%cornerVelocity(22) = velocityFunc(vector(x1,y2,z3),grid)
          thisOctal%cornerrho(23) = densityFunc(vector(x2,y2,z3),grid)
          thisOctal%cornerVelocity(23) = velocityFunc(vector(x2,y2,z3),grid)
          thisOctal%cornerrho(24) = densityFunc(vector(x3,y2,z3),grid)
          thisOctal%cornerVelocity(24) = velocityFunc(vector(x3,y2,z3),grid)
          thisOctal%cornerrho(25) = densityFunc(vector(x1,y3,z3),grid)
          thisOctal%cornerVelocity(25) = velocityFunc(vector(x1,y3,z3),grid)
          thisOctal%cornerrho(26) = densityFunc(vector(x2,y3,z3),grid)
          thisOctal%cornerVelocity(26) = velocityFunc(vector(x2,y3,z3),grid)
          thisOctal%cornerrho(27) = densityFunc(vector(x3,y3,z3),grid)
          thisOctal%cornerVelocity(27) = velocityFunc(vector(x3,y3,z3),grid)

       else ! cylindrical 
          if (thisOctal%splitAzimuthally) then
             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi 
             phi3 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize

             ! bottom level

             thisOctal%cornerrho(1) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1),grid)
             thisOctal%cornerrho(2) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z1),grid)
             thisOctal%cornerrho(3) = densityFunc(vector(r1*cos(phi3),r1*sin(phi3),z1),grid)
             thisOctal%cornerrho(4) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z1),grid)
             thisOctal%cornerrho(5) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z1),grid)
             thisOctal%cornerrho(6) = densityFunc(vector(r2*cos(phi3),r2*sin(phi3),z1),grid)
             thisOctal%cornerrho(7) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z1),grid)
             thisOctal%cornerrho(8) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1),grid)
             thisOctal%cornerrho(9) = densityFunc(vector(r3*cos(phi3),r3*sin(phi3),z1),grid)

             ! middle level

             thisOctal%cornerrho(10) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2),grid)
             thisOctal%cornerrho(11) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z2),grid)
             thisOctal%cornerrho(12) = densityFunc(vector(r1*cos(phi3),r1*sin(phi3),z2),grid)
             thisOctal%cornerrho(13) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z2),grid)
             thisOctal%cornerrho(14) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z2),grid)
             thisOctal%cornerrho(15) = densityFunc(vector(r2*cos(phi3),r2*sin(phi3),z2),grid)
             thisOctal%cornerrho(16) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z2),grid)
             thisOctal%cornerrho(17) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2),grid)
             thisOctal%cornerrho(18) = densityFunc(vector(r3*cos(phi3),r3*sin(phi3),z2),grid)

             ! top level

             thisOctal%cornerrho(10) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3),grid)
             thisOctal%cornerrho(11) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z3),grid)
             thisOctal%cornerrho(12) = densityFunc(vector(r1*cos(phi3),r1*sin(phi3),z3),grid)
             thisOctal%cornerrho(13) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z3),grid)
             thisOctal%cornerrho(14) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z3),grid)
             thisOctal%cornerrho(15) = densityFunc(vector(r2*cos(phi3),r2*sin(phi3),z3),grid)
             thisOctal%cornerrho(16) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z3),grid)
             thisOctal%cornerrho(17) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3),grid)
             thisOctal%cornerrho(18) = densityFunc(vector(r3*cos(phi3),r3*sin(phi3),z3),grid)

          else

             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize


             ! bottom level

             thisOctal%cornerrho(1) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z1),grid)
             thisOctal%cornerrho(2) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z1),grid)
             thisOctal%cornerrho(3) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z1),grid)
             thisOctal%cornerrho(4) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z1),grid)
             thisOctal%cornerrho(5) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z1),grid)
             thisOctal%cornerrho(6) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z1),grid)

             ! middle level

             thisOctal%cornerrho(7) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z2),grid)
             thisOctal%cornerrho(8) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z2),grid)
             thisOctal%cornerrho(9) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z2),grid)
             thisOctal%cornerrho(10) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z2),grid)
             thisOctal%cornerrho(11) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z2),grid)
             thisOctal%cornerrho(12) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z2),grid)

             ! top level

             thisOctal%cornerrho(13) = densityFunc(vector(r1*cos(phi1),r1*sin(phi1),z3),grid)
             thisOctal%cornerrho(14) = densityFunc(vector(r1*cos(phi2),r1*sin(phi2),z3),grid)
             thisOctal%cornerrho(15) = densityFunc(vector(r2*cos(phi1),r2*sin(phi1),z3),grid)
             thisOctal%cornerrho(16) = densityFunc(vector(r2*cos(phi2),r2*sin(phi2),z3),grid)
             thisOctal%cornerrho(17) = densityFunc(vector(r3*cos(phi1),r3*sin(phi1),z3),grid)
             thisOctal%cornerrho(18) = densityFunc(vector(r3*cos(phi2),r3*sin(phi2),z3),grid)

          endif
       endif
    else
       
       
    ! we first store the values we use to assemble the position vectors
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values
       
       thisOctal%cornerrho(1) = densityFunc(vector(x1,0.d0,z1),grid)
       thisOctal%cornerrho(2) = densityFunc(vector(x2,0.d0,z1),grid)
       thisOctal%cornerrho(3) = densityFunc(vector(x3,0.d0,z1),grid)
       thisOctal%cornerrho(4) = densityFunc(vector(x1,0.d0,z2),grid)
       thisOctal%cornerrho(5) = densityFunc(vector(x2,0.d0,z2),grid)
       thisOctal%cornerrho(6) = densityFunc(vector(x3,0.d0,z2),grid)
       thisOctal%cornerrho(7) = densityFunc(vector(x1,0.d0,z3),grid)
       thisOctal%cornerrho(8) = densityFunc(vector(x2,0.d0,z3),grid)
       thisOctal%cornerrho(9) = densityFunc(vector(x3,0.d0,z3),grid)
    endif
667 continue
    
  END SUBROUTINE fillDensityCorners
           
  real(double) function molebenchDensity(position,grid) result(rho)

    type(VECTOR) :: position
    TYPE(gridtype), INTENT(IN) :: grid
    logical, save :: firsttime = .true.
    integer, parameter :: nr = 50
    real(double),save :: r(nr), nh2(nr), junk
    real(double) :: r1, t1, t2, nh2out
    integer :: i

    if (firsttime) then
       open(31, file="model_1.dat", status="old", form="formatted") ! Model 2 in the Hogerheijde 2000 paper. 
       do i = nr,1,-1                                             
          read(31,*) r(i), nh2(i), junk
       enddo
       r = r * 1.d-10
       close(31)
       firsttime = .false.
    endif

    r1 = modulus(position)

    if(r1 > r(nr) .or. r1 < r(1)) then 
       nh2out = 1.d-20
       rho = 1.d-20 * 2. * mhydrogen
!       nh2out = nh2(1)
!       rho = nh2(1) * 2.d0 * mhydrogen
    endif

    if ((r1 > r(1)).and.(r1 < r(nr))) then
       call locate(r, nr, r1, i)
!       t2 = (r1 - r(i))/(r(i+1)-r(i)) ! linear but know its a power law so use better interpolation
       
       t1 = log(r1/r(i))/log(r(i+1)/r(i))
       
       nh2out = exp((1.d0 - t1) * log(nh2(i))  +  t1 * log(nh2(i+1)))
       rho = nh2out * 2.d0 * mhydrogen
       
    endif
  end function molebenchDensity

  real(double) function averagerhofromoctal(thisoctal,subcell) result(rho)
    
    type(OCTAL) :: thisoctal
    integer :: subcell

       SELECT CASE(subcell)
               
       CASE(1)
          rho =( &
               thisoctal%cornerrho( 1) + &
               thisoctal%cornerrho( 2) + &
               thisoctal%cornerrho( 4) + &
               thisoctal%cornerrho( 5) + &
               thisoctal%cornerrho(10) + &
               thisoctal%cornerrho(11) + &
               thisoctal%cornerrho(13) + &
               thisoctal%cornerrho(14)) * 1.25d-1
          
       CASE(2)
          rho =( &
               thisoctal%cornerrho( 2) + &
               thisoctal%cornerrho( 3) + &
               thisoctal%cornerrho( 5) + &
               thisoctal%cornerrho( 6) + &
               thisoctal%cornerrho(11) + &
               thisoctal%cornerrho(12) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(15)) * 1.25d-1
          
       CASE(3)
          rho =( &
               thisoctal%cornerrho( 4) + &
               thisoctal%cornerrho( 5) + &
               thisoctal%cornerrho( 7) + &
               thisoctal%cornerrho( 8) + &
               thisoctal%cornerrho(13) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(16) + &
               thisoctal%cornerrho(17)) * 1.25d-1
          
       CASE(4)
          rho =( &
               thisoctal%cornerrho( 5) + &
               thisoctal%cornerrho( 6) + &
               thisoctal%cornerrho( 8) + &
               thisoctal%cornerrho( 9) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(15) + &
               thisoctal%cornerrho(17) + &
               thisoctal%cornerrho(18)) * 1.25d-1
          
       CASE(5)
          rho =( &
               thisoctal%cornerrho(10) + &
               thisoctal%cornerrho(11) + &
               thisoctal%cornerrho(13) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(19) + &
               thisoctal%cornerrho(20) + &
               thisoctal%cornerrho(22) + &
               thisoctal%cornerrho(23)) * 1.25d-1
          
       CASE(6)
          rho =( &
               thisoctal%cornerrho(11) + &
               thisoctal%cornerrho(12) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(15) + &
               thisoctal%cornerrho(20) + &
               thisoctal%cornerrho(21) + &
               thisoctal%cornerrho(23) + &
               thisoctal%cornerrho(24)) * 1.25d-1
          
       CASE(7)
          rho =( &
               thisoctal%cornerrho(13) + &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(16) + &
               thisoctal%cornerrho(17) + &
               thisoctal%cornerrho(22) + &
               thisoctal%cornerrho(23) + &
               thisoctal%cornerrho(25) + &
               thisoctal%cornerrho(26)) * 1.25d-1
          
       CASE(8)
          rho =( &
               thisoctal%cornerrho(14) + &
               thisoctal%cornerrho(15) + &
               thisoctal%cornerrho(17) + &
               thisoctal%cornerrho(18) + &
               thisoctal%cornerrho(23) + &
               thisoctal%cornerrho(24) + &
               thisoctal%cornerrho(26) + &
               thisoctal%cornerrho(27)) * 1.25d-1
          
       end select

  end function averagerhofromoctal

  subroutine pathTest(grid)
    type(GRIDTYPE) :: grid
    integer :: i, nPaths
    type(VECTOR) :: rVec, uHat
    real(double) :: r, s, tot
    integer :: cpuTime, startTime, endTime

    nPaths = 100000
    write(*,*) "Doing path test..."
    tot = 0.d0
    call unixTimes(cpuTime, startTime)
    do i = 1, nPaths
       call random_Number(r)
       rVec%x = grid%octreeRoot%centre%x + (2.d0*r-1.d0) * grid%octreeRoot%subcellSize
       call random_Number(r)
       rVec%y = grid%octreeRoot%centre%y + (2.d0*r-1.d0) * grid%octreeRoot%subcellSize
       call random_Number(r)
       rVec%z = grid%octreeRoot%centre%z + (2.d0*r-1.d0) * grid%octreeRoot%subcellSize
       uHat = randomUnitVector()
       call testToBoundary(grid, rVec, uHat, s)
       tot = tot+s
    enddo
    call unixTimes(cpuTime, endTime)
    write(*,*) "Path integrals per second: ",real(nPaths)/real(endTime - startTime)
    write(*,*) "Check: ",tot/real(nPaths)

    tot = 0.d0
    call unixTimes(cpuTime, startTime)
    do i = 1, nPaths
       call random_Number(r)
       rVec%x = grid%octreeRoot%centre%x + (2.d0*r-1.d0) * grid%octreeRoot%subcellSize
       call random_Number(r)
       rVec%y = grid%octreeRoot%centre%y + (2.d0*r-1.d0) * grid%octreeRoot%subcellSize
       call random_Number(r)
       rVec%z = grid%octreeRoot%centre%z + (2.d0*r-1.d0) * grid%octreeRoot%subcellSize
       uHat = randomUnitVector()
       call testToBoundary2(grid, rVec, uHat, s)
       tot = tot + s
    enddo
    call unixTimes(cpuTime, endTime)
    write(*,*) "Path integrals per second with neighbour pointers: ",real(nPaths)/real(endTime - startTime)
    write(*,*) "Check: ",tot/real(nPaths)
  end subroutine pathTest
   
  subroutine testToBoundary(grid, rVec, uHat, totDist)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, uHat, currentPosition
    type(OCTAL), pointer :: thisOctal, sOctal
    integer :: subcell
    real(double) :: distToNextCell, totDist
    real(double) :: fudgeFac = 0.01d0

    totDist = 0.d0
    currentPosition = rVec
    CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)


    do while (inOctal(grid%octreeRoot, currentPosition))

       call findSubcellLocal(currentPosition,thisOctal,subcell)

       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, uHat, DisttoNextCell, sOctal)
  
       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*uHat
       totDist = totDist + distToNextCell

    end do
  end subroutine testToBoundary

  subroutine testToBoundary2(grid, rVec, uHat, totDist)
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec, uHat, currentPosition
    type(OCTAL), pointer :: thisOctal, sOctal, tempOctal
    integer :: subcell, sSubcell
    real(double) :: distToNextCell, totDist
    real(double) :: fudgeFac = 0.000001d0
    integer :: tempSubcell

    totDist = 0.d0
    currentPosition = rVec
    CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)

    sOctal => thisOctal
    sSubcell = subcell

    do while (sSubcell > 0)

       call distanceToCellBoundary(grid, currentPosition, uHat, DisttoNextCell, thisOctal, subcell)

       currentPosition = currentPosition + distToNextCell*uHat
       call getNeighbourFromPointOnFace(currentPosition, uHat, thisOctal, subcell, sOctal, sSubcell)
       thisOctal => sOctal
       subcell = sSubcell

!       currentPosition = currentPosition + (fudgeFac*grid%halfSmallestSubcell)*uHat



       totDist = totDist + distToNextCell

    end do
  end subroutine testToBoundary2
END MODULE amr_mod

