MODULE amr_mod
  ! routines for adaptive mesh refinement. nhs
  ! twod stuff added by tjh started 25/08/04


  USE messages_mod
  USE vector_mod          ! vector maths routines
  USE kind_mod            ! variable kind parameters    
  USE octal_mod           ! type definition for the grid elements
  USE gridtype_mod        ! type definition for the 3-d grid 
  USE parameters_mod      ! parameters for specific geometries
  USE jets_mod            ! 
  USE luc_cir3d_class 
  USE cmfgen_class
  USE constants_mod, only: cSpeed, pi
  USE sph_data_class
  USE cluster_class
  USE density_mod
  USE wr104_mod
  USE utils_mod
  USE gas_opacity_mod
  USE math_mod2!, only : mnewt  
  USE atom_mod
  USE source_mod

  IMPLICIT NONE


CONTAINS

  SUBROUTINE calcValuesAMR(thisOctal,subcell,grid, sphData, stellar_cluster, inherit, interp)
    ! calculates the variables describing one subcell of an octal.
    ! each geometry that can be used with AMR should be described here, 
    !   otherwise the program will print a warning and exit.
   
    IMPLICIT NONE

    TYPE(octal)    :: thisOctal   ! the octal being changed
    INTEGER, INTENT(IN)           :: subcell   ! the subcell being changed
    TYPE(gridtype), INTENT(INOUT) :: grid      ! the grid
    !
    TYPE(sph_data), optional, INTENT(IN)  :: sphData   ! Matthew's SPH data.
    TYPE(cluster), optional, INTENT(IN)   :: stellar_cluster
    LOGICAL, OPTIONAL :: inherit               ! inherit densities, temp, etc of parent
    LOGICAL, OPTIONAL :: interp                ! interpolate densities, temp, etc of parent
    TYPE(octal), POINTER :: parentOctal
    INTEGER :: parentSubcell
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


    if (inheritProps) then
       thisOctal%rho(subcell) = parentOctal%rho(parentSubcell)
       thisOctal%temperature(subcell) = parentOctal%temperature(parentSubcell)
       thisOctal%etaCont(subcell) = parentOctal%etaCont(parentSubcell)
       thisOctal%inFlow(subcell) = parentOctal%inFlow(parentSubcell)
       thisOctal%velocity(subcell) = parentOctal%velocity(parentSubcell)
       thisOctal%biasCont3D(subcell) = parentOctal%biasCont3D(parentSubcell)
       thisOctal%etaLine(subcell) = parentOctal%etaLine(parentSubcell)
       thisOctal%chiLine(subcell) = parentOctal%chiLine(parentSubcell)
       thisOctal%dustTypeFraction(subcell,:) = parentOctal%dustTypeFraction(parentSubcell,:)
       thisOctal%oldFrac(subcell) = parentOctal%oldFrac(parentSubcell)
       if (associated(thisOctal%ionFrac)) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif
       thisOctal%nh(subcell) = parentOctal%nh(parentsubcell)
       thisOctal%ne(subcell) = parentOctal%ne(parentsubcell)
    else if (interpolate) then
       thisOctal%etaCont(subcell) = parentOctal%etaCont(parentSubcell)
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
       thisOctal%nh(subcell) = parentOctal%nh(parentsubcell)
       thisOctal%ne(subcell) = parentOctal%ne(parentsubcell)

       parentOctal%hasChild(parentsubcell) = .false.
       call interpFromParent(subcellCentre(thisOctal, subcell), thisOctal%subcellSize, grid, &
            thisOctal%temperature(subcell), thisOctal%rho(subcell), thisOctal%dusttypeFraction(subcell, :))
       parentOctal%hasChild(parentsubcell) = .true.

    else

    SELECT CASE (grid%geometry)

    CASE ("ttauri")
      CALL calcTTauriMassVelocity(thisOctal,subcell,grid)
      
    CASE ("jets")
      CALL calcJetsMassVelocity(thisOctal,subcell,grid)

    CASE ("luc_cir3d")
      CALL calc_cir3d_mass_velocity(thisOctal,subcell)

    CASE ("cmfgen")
      CALL cmfgen_mass_velocity(thisOctal,subcell)

    CASE ("testamr")
       CALL calcTestDensity(thisOctal,subcell,grid)

    CASE("lexington")
       CALL calcLexington(thisOctal, subcell, grid)
       if (thisOctal%nDepth > 1) then
          thisOctal%ionFrac(subcell,:) = parentOctal%ionFrac(parentsubcell,:)
       endif

    CASE("proto")
       CALL calcProtoDensity(thisOctal,subcell,grid)

    CASE("wrshell")
       CALL calcWRShellDensity(thisOctal,subcell,grid)

    CASE ("spiralwind")
       CALL spiralWindSubcell(thisOctal, subcell ,grid)
       

    CASE("cluster")
       ! using a routine in cluster_class.f90
       call assign_density(thisOctal,subcell, sphData, grid%geometry, stellar_cluster)
       
    CASE("wr104")
       ! using a routine in cluster_class.f90
       call assign_density(thisOctal,subcell, sphData, grid%geometry)

    CASE ("benchmark")
       CALL benchmarkDisk(thisOctal, subcell ,grid)

    CASE ("shakara","aksco")
       CALL shakaraDisk(thisOctal, subcell ,grid)

    CASE ("warpeddisc")
       CALL warpedDisk(thisOctal, subcell ,grid)

    CASE("ppdisk")
       CALL calcPPDiskDensity(thisOctal,subcell,grid)

    CASE("melvin")
       CALL assign_melvin(thisOctal,subcell,grid)

    CASE("whitney")
       CALL assign_whitney(thisOctal,subcell,grid)

    CASE("toruslogo")
       CALL assign_toruslogo(thisOctal,subcell,grid)

    CASE("clumpydisc")
       CALL assign_clumpydisc(thisOctal, subcell, grid)

    CASE ("windtest")
      CALL calcWindTestValues(thisOctal,subcell,grid)
       
    CASE DEFAULT
      WRITE(*,*) "! Unrecognised grid geometry: ",TRIM(grid%geometry)
      STOP

    END SELECT
 
!    CALL fillGridDummyValues(thisOctal,subcell, grid)
 
    end if
 
  END SUBROUTINE calcValuesAMR

  SUBROUTINE initFirstOctal(grid, centre, size, twod, sphData, stellar_cluster, nDustType)
    ! creates the first octal of a new grid (the root of the tree).
    ! this should only be used once; use addNewChild for subsequent
    !  additions.
    use input_variables, only : cylindrical, photoionization

    IMPLICIT NONE
    
    TYPE(gridtype), INTENT(INOUT)    :: grid 
    TYPE(octalVector), INTENT(IN)    :: centre ! coordinates of the grid centre
    REAL, INTENT(IN)                 :: size 
      ! 'size' should be the vertex length of the cube that contains the whole
      !   of the simulation space, *not* the size of a subcell.
    TYPE(sph_data), optional, intent(in)   :: sphData   ! Matthew's SPH model data
    TYPE(cluster), optional, INTENT(IN)   :: stellar_cluster
    LOGICAL :: twod  ! true if this is a twoD amr grid
    INTEGER :: subcell ! loop counter 
    INTEGER :: nDustType ! number of different dust types


    ALLOCATE(grid%octreeRoot)
    
    ! allocate any variables that need to be 
    if (.not.grid%oneKappa) then
       ALLOCATE(grid%octreeRoot%kappaAbs(8,grid%nOpacity))
       ALLOCATE(grid%octreeRoot%kappaSca(8,grid%nOpacity))
       grid%octreeRoot%kappaAbs = 1.e-30
       grid%octreeRoot%kappaSca = 1.e-30
    endif

    ALLOCATE(grid%octreeRoot%N(8,grid%maxLevels))

    if (twoD) then
       grid%octreeRoot%twoD = .true.
       grid%octreeRoot%threeD = .false.
       grid%octreeRoot%maxChildren = 4
    else
       grid%octreeRoot%cylindrical = .false.
       grid%octreeRoot%twoD = .false.
       grid%octreeRoot%threeD = .true.
       grid%octreeRoot%maxChildren = 8
       if (cylindrical) then
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
    grid%octreeRoot%indexChild = -999 ! values are undefined
    grid%octreeRoot%probDistLine = 0.0
    grid%octreeRoot%probDistCont = 0.0
    NULLIFY(grid%octreeRoot%parent)   ! tree root does not have a parent
    NULLIFY(grid%octreeRoot%child)    ! tree root does not yet have children
    ! initialize some values
    grid%octreeRoot%biasLine3D = 1.0 
    grid%octreeRoot%biasCont3D = 1.0 
    grid%octreeRoot%velocity = vector(1.e-30,1.e-30,1.e-30)
    grid%octreeRoot%cornerVelocity= vector(1.e-30,1.e-30,1.e-30)
    grid%octreeRoot%chiLine = -9.9e9
    grid%octreeRoot%etaLine = -9.9e9
    grid%octreeRoot%etaCont = -9.9e9 
    grid%octreeRoot%nDirectPhotons = -1
    grid%octreeRoot%rho = 1.e-30
    grid%octreeRoot%N = -9.9e9
    grid%octreeRoot%Ne = -9.9e9
    grid%octreeRoot%nTot = -9.9e9
    grid%octreeRoot%changed = .false.
    grid%octreeRoot%dustType = 1
    ALLOCATE(grid%octreeRoot%dusttypefraction(8,  nDustType))
    grid%octreeroot%dustTypeFraction(1:8,1:nDustType) = 0.d0
    grid%octreeroot%dustTypeFraction(1:8,1) = 1.d0
    grid%octreeRoot%inflow = .true.
    grid%octreeRoot%gasOpacity = .false.
    grid%octreeRoot%diffusionApprox = .false.
    grid%octreeRoot%nDiffusion = 0.
    grid%octreeRoot%oldFrac = 1.e-5
    if (photoionization) then
       allocate(grid%octreeRoot%ionFrac(1:grid%octreeRoot%maxChildren, 1:grid%nIon))
       allocate(grid%octreeRoot%photoionCoeff(1:grid%octreeRoot%maxChildren, 1:grid%nIon))
    endif

    if (associated(grid%octreeRoot%ionFrac)) grid%octreeRoot%ionFrac = 1.e-30
    
    if (associated(grid%octreeRoot%photoIonCoeff)) then
       grid%octreeRoot%photoIonCoeff = 0.d0
    endif

    select case (grid%geometry)
       case("cluster")
          ! Initially we copy the idecies of particles (in SPH data)
          ! to the root node. The indecies will copy over to
           ! to the subcells if the particles are in the subcells.
          ! This will allow us to work with the subsets of gas particle
          ! list hence reduces the computation time when we are 
          ! splitting/constructing the octree data structure. 
          !
          ! Using the routine in grid_mod.f90
          call copy_sph_index_to_root(grid, sphData)
          !
          DO subcell = 1, grid%octreeRoot%maxChildren
             ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, sphData, stellar_cluster)
             ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO
       case("wr104")
          ! Using the routine in grid_mod.f90
          call copy_sph_index_to_root(grid, sphData)
          !
          DO subcell = 1, grid%octreeRoot%maxChildren
             ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, sphData, stellar_cluster)
             ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO
       case DEFAULT
          DO subcell = 1, grid%octreeRoot%maxChildren
             ! calculate the values at the centre of each of the subcells
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid)
             ! label the subcells
             grid%octreeRoot%label(subcell) = subcell
          END DO
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


  SUBROUTINE addNewChild(parent, iChild, grid, adjustGridInfo, sphData, &
                         stellar_cluster, inherit, interp, splitAzimuthally)
    ! adds one new child to an octal

    USE input_variables, ONLY : nDustType, cylindrical, photoionization
    IMPLICIT NONE
    
    TYPE(octal), TARGET, INTENT(INOUT) :: parent ! the parent octal 
    INTEGER, INTENT(IN)  :: iChild     ! the label (1-8) of the subcell gaining the child 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that
                                          !   calculate the variables stored in
                                          !   the tree.
    LOGICAL, INTENT(IN) :: adjustGridInfo
    LOGICAL, optional :: splitAzimuthally
      ! whether these variables should be updated: 
      !   grid%nOctals, grid%maxDepth, grid%halfSmallestSubcell    
    ! For only cluster geometry ...
    TYPE(sph_data), OPTIONAL, INTENT(IN) :: sphData 
    TYPE(cluster), OPTIONAL, INTENT(IN)   :: stellar_cluster
    LOGICAL, OPTIONAL :: inherit       ! inherit densities, temps, etc from parent
    LOGICAL, OPTIONAL :: interp        ! interpolate densities, temps, etc from parent    
    
    INTEGER       :: subcell           ! loop counter
    INTEGER, SAVE :: counter = 9       ! keeps track of the current subcell label
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: nChildren         ! number of children the parent octal has
    INTEGER       :: newChildIndex     ! the storage location for the new child
    logical :: inheritProps, interpolate
    type(OCTALVECTOR) :: rVec
    ! array of octals that may be needed for temporarily storing child octals
    

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

    ALLOCATE(parent%child(newChildIndex)%N(8,grid%maxLevels))
    ! set up the new child's variables
    parent%child(newChildIndex)%threeD = parent%threeD
    parent%child(newChildIndex)%twoD = parent%twoD
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
    parent%child(newChildIndex)%nDirectPhotons = -1
    parent%child(newChildIndex)%indexChild = -999 ! values are undefined
    parent%child(newChildIndex)%nDepth = parent%nDepth + 1
    parent%child(newChildIndex)%centre = subcellCentre(parent,iChild)
    if (parent%cylindrical) then
       parent%child(newChildIndex)%r = subcellRadius(parent,iChild)
    endif
    parent%child(newChildIndex)%probDistLine = 0.0
    parent%child(newChildIndex)%probDistCont = 0.0
    parent%child(newChildIndex)%biasLine3D = 1.0 
    parent%child(newChildIndex)%biasCont3D = 1.0 
    parent%child(newChildIndex)%velocity = vector(1.e-30,1.e-30,1.e-30)
    parent%child(newChildIndex)%cornerVelocity = vector(1.e-30,1.e-30,1.e-30)
    parent%child(newChildIndex)%chiLine = 1.e-30
    parent%child(newChildIndex)%etaLine = 1.e-30
    parent%child(newChildIndex)%etaCont = 1.e-30
    parent%child(newChildIndex)%rho = 1.e-30
    parent%child(newChildIndex)%N = 1.e-30
    parent%child(newChildIndex)%dusttype = 1
    ALLOCATE(parent%child(newChildIndex)%dusttypefraction(8,  nDustType))
    parent%child(newChildIndex)%dustTypeFraction(1:8,1:nDustType) = 0.d0
    parent%child(newChildIndex)%dustTypeFraction(1:8,1) = 1.d0
    parent%child(newChildIndex)%gasOpacity = .false.
    parent%child(newChildIndex)%Ne = 1.e-30
    parent%child(newChildIndex)%temperature = 3.0
    parent%child(newChildIndex)%nTot = 1.e-30
    parent%child(newChildIndex)%eDens = 1.d-10
    parent%child(newChildIndex)%changed = .false.
    parent%child(newChildIndex)%diffusionApprox = .false.
    parent%child(newChildindex)%nDiffusion  = 0.
    parent%child(newChildindex)%oldFrac = 1.
    parent%child(newChildindex)%ncrossings = 10000
    if (photoionization) then
       allocate(parent%child(newChildIndex)%ionFrac(1:parent%maxChildren, 1:grid%nIon))
       allocate(parent%child(newChildIndex)%photoionCoeff(1:parent%maxChildren, 1:grid%nIon))
    endif

    if (associated(parent%child(newChildindex)%ionFrac)) parent%child(newChildIndex)%ionFrac = 1.e-30
    
    if (associated(parent%child(newChildIndex)%photoIonCoeff)) then
       parent%child(newChildIndex)%photoIonCoeff = 0.d0
    endif
       

    IF (PRESENT(sphData)) THEN
      ! updates the sph particle list.           
!      CALL update_particle_list(parent, newChildIndex, newChildIndex, sphData)
      CALL update_particle_list(parent, iChild, newChildIndex, sphData)
!      write(*,*) "update_particle_list broken"
    END IF 
    
    ! put some data in the four/eight subcells of the new child
    DO subcell = 1, parent%child(newChildIndex)%maxChildren
      CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid, sphData, &
                         stellar_cluster, inherit=inheritProps, interp=interpolate)
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


  RECURSIVE SUBROUTINE splitGrid(thisOctal,amrLimitScalar,amrLimitScalar2,grid, &
       sphData, stellar_cluster, setChanged)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    IMPLICIT NONE

    TYPE(OCTAL), INTENT(INOUT) :: thisOctal
    real(double), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
                                          !   routines that this subroutine calls

    ! Object containg the output from the (Mattew's) SPH code.
    ! This will be just passed to decideSplit routine.
    TYPE(sph_data), OPTIONAL, intent(in) :: sphData
    TYPE(cluster), OPTIONAL,  intent(in) :: stellar_cluster
    LOGICAL, INTENT(IN), OPTIONAL :: setChanged
    !
    INTEGER              :: iSubcell, iIndex ! loop counters
    INTEGER              :: i, j, k
    logical :: splitInAzimuth


    DO iSubcell = 1, thisOctal%maxChildren, 1
      
      IF (decideSplit(thisOctal,iSubcell,amrLimitScalar,amrLimitScalar2,grid,&
            splitInAzimuth, &
            sphData, stellar_cluster)) THEN

        CALL addNewChild(thisOctal, iSubcell, grid, adjustGridInfo=.TRUE., &
                         sphData=sphData, stellar_cluster=stellar_cluster, splitAzimuthally=splitInAzimuth)
       
        if (.not.thisOctal%hasChild(isubcell)) then
          write(*,*) "add child failed in splitGrid"
          do; enddo
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
          do;enddo
        endif
      endif
    enddo   

    if (any(thisOctal%haschild(1:thisOctal%maxChildren)).and.(thisOctal%nChildren==0)) then
      write(*,*) "nchildren screw up"
      do;enddo
    endif

    
    DO iIndex = 1, thisOctal%nChildren, 1
      
      CALL splitGrid(thisOctal%child(iIndex),amrLimitScalar,amrLimitScalar2,grid,&
                     sphData, stellar_cluster, setChanged)
      
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
      REAL    :: a
      real(oct) :: valueX
      real(oct) :: valueY
      real(oct) :: valueA
      real(oct) :: valueJ
      TYPE(octal), POINTER :: Acontent
      LOGICAL(KIND=logic), DIMENSION(8) :: AinUse
      TYPE(octalVector)    :: starPos
      
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
    
  
  RECURSIVE SUBROUTINE finishGrid(thisOctal,grid,gridConverged)
    ! takes the octree grid that has been created using 'splitGrid'
    !   and calculates all the other variables in the model.
    ! this should be called once the structure of the grid is complete.
    ! the gridConverged variable should be set .TRUE. when the entire
    !   grid has been finished. Until that happens, this subroutine 
    !   will be called repeatedly.
    
    USE input_variables, ONLY : useHartmannTemp
    
    IMPLICIT NONE

    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    LOGICAL, INTENT(INOUT) :: gridConverged
    
    TYPE(octal), POINTER   :: child
  
    INTEGER :: subcell, iChild
     
    ! all of the work that must be done recursively goes here: 
    DO subcell = 1, thisOctal%maxChildren, 1
   
      SELECT CASE (grid%geometry)

      CASE ("ttauri")
        IF (.NOT. useHartmannTemp) &
          CALL calcTTauriTemperature(thisOctal,subcell)
        gridConverged = .TRUE.
        
      CASE ("spiralwind")
        gridConverged = .TRUE.

      CASE ("jets")
        CALL calcJetsTemperature(thisOctal,subcell, grid)
        gridConverged = .TRUE.
        
      CASE ("luc_cir3d")
        CALL calc_cir3d_temperature(thisOctal,subcell)
        gridConverged = .TRUE.

      CASE ("cmfgen")
        CALL calc_cmfgen_temperature(thisOctal,subcell)
        gridConverged = .TRUE.
        
      CASE ("testamr","proto")
        gridConverged = .TRUE.

      CASE ("windtest")
        gridConverged = .TRUE.

      CASE("benchmark","shakara","aksco", "melvin","clumpydisc", &
           "lexington", "warpeddisc", "whitney")
         gridConverged = .TRUE.

      CASE ("cluster","wr104")
        call assign_grid_values(thisOctal,subcell, grid)
        gridConverged = .TRUE.

      CASE ("ppdisk","wrshell","toruslogo")
        gridConverged = .TRUE.


      CASE DEFAULT
        WRITE(*,*) "! Unrecognised grid geometry: ",trim(grid%geometry)
        STOP
      
      END SELECT
      
    END DO
   
    DO iChild = 1, thisOctal%nChildren, 1
      child => thisOctal%child(iChild)
      CALL finishGrid(child,grid,gridConverged)
    END DO

    ! any stuff that gets done *after* the finishGrid recursion goes here:
    IF (.NOT. ASSOCIATED(thisOctal%parent)) THEN
     
      ! nothing at the moment

    END IF


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
             velocityDeriv,chiLine,levelPop,rho, temperature, Ne, inflow)
    ! samples the grid at points along the path.
    ! this should be called by the program, instead of calling 
    !   returnSamples directly, because this checks the start and finish
    !   points are within the grid bounds. returnSamples *assumes* this 
    !   criterion is met. also, the path of the photon is checked for 
    !   intersection with the stellar surface(s), disc etc. 
 
    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN)      :: startPoint ! photon start point
    TYPE(octalVector), INTENT(IN)      :: direction  ! photon direction 
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

    TYPE(octalVector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    TYPE(octalVector)       :: currentPoint ! current position of ray 
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(octal)             :: octree       ! the octree structure within 'grid'
    TYPE(octalVector)       :: directionNormalized
    real(oct)    :: distanceLimit ! max length of ray before aborting
    ! margin is the size of the region around the edge of a subcell
    !   where numerical inaccuracies may cause problems.
    real(oct)    :: margin
 
    ! variables for testing special cases (stellar intersections etc.)
    TYPE(octalVector)       :: starPosition       ! position vector of stellar centre
    TYPE(octalVector)       :: diskNormal         ! disk normal vector
    real(oct)    :: rStar              ! stellar radius
    real(oct)    :: endLength          ! max path length of photon
    TYPE(octalVector)       :: endPoint           ! where photon leaves grid
    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
    real(oct)    :: distanceFromOrigin ! closest distance of plane from origin
    TYPE(octalVector)       :: diskIntersection   ! point of photon intersection with disk
    real(oct)    :: starIntersectionDistance1 ! distance to first  intersection
    real(oct)    :: starIntersectionDistance2 ! distance to second intersection
    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
    real(oct)    :: intersectionRadius ! disk radius when photon intersects
    real(oct)    :: diskDistance       ! distance to disk intersection
    real(oct)    :: distanceThroughStar! distance of chord through star
    TYPE(octalVector)       :: dummyStartPoint    ! modified start point 
    real(oct), PARAMETER :: fudgefactor = 1.00001 ! overestimates stellar size
    
    
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
         grid%geometry(1:9) == "luc_cir3d" .or. grid%geometry(1:6) == "cmfgen" ) THEN
      
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
                     Ne=Ne, inFlow=inFlow)
      
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
                    temperature=temperature, Ne=Ne, inFlow=inFlow)

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
                   kappaAbs=kappaAbs,kappaSca=kappaSca, velocity=velocity,     &
                   velocityDeriv=velocityDeriv,chiLine=chiLine,               &
                   levelPop=levelPop,rho=rho,temperature=temperature,         &
                   Ne=Ne, inFlow=inFlow)
    END IF
      
  END SUBROUTINE startReturnSamples

  
  RECURSIVE SUBROUTINE returnSamples (currentPoint,startPoint,locator,direction, &
             octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,usePops, &
             iLambda,error,margin,distanceLimit,kappaAbs,kappaSca,velocity,      &
             velocityDeriv,chiLine,levelPop,rho,temperature, Ne, inFlow)
    ! this uses a recursive ray traversal algorithm to sample the octal
    !   grid at points along the path of the ray. 
    ! no checks are made that the ray lies within the boundaries of the
    !  grid, so this subroutine not be called directly. use  
    !  'startReturnSamples' instead.

    IMPLICIT NONE

    ! note: in the code comments, the terms 'subcell' and 'cell' are
    !  mostly interchangeable.

    TYPE(octalVector), INTENT(INOUT)    :: currentPoint ! current ray position
    TYPE(octalVector), INTENT(IN)       :: startPoint   ! initial ray position
    TYPE(octalVector)                   :: locator, rotloc
                  ! 'locator' is used to indicate a point that lies within the  
                  !   *next* cell of the octree that the ray will interesect.
                  !   initially this will be the same as the currentPoint
      
    TYPE(octalVector), INTENT(IN)       :: direction    ! ray's direction vector 
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
    
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL :: kappaAbs      ! continuous absorption opacities
    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca     ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL :: chiLine       ! line opacities
    REAL(double),DIMENSION(:),OPTIONAL               :: rho           ! density at sample points
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL   :: temperature! in [K]
    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL   :: Ne         ! electron density
    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL :: levelPop   ! level populations 
    logical, DIMENSION(:),INTENT(INOUT),OPTIONAL       :: inFlow     ! inflow flag
    


    TYPE(octalVector)      :: exitPoint      ! where ray leaves current cell
    TYPE(octalVector)      :: centre         ! centre of current subcell
    real(oct)   :: subcellSize    ! size of current subcell
    real(oct)   :: minWallDistance ! distance to *nearest* wall
    
    real(oct)   :: length         ! distance from start of the ray's path.
    real(oct)   :: sampleLength   ! distance interval between samples 
    real(oct)   :: trialLength    ! trial distance to a possible next sample
    TYPE(octalVector)      :: trialPoint     ! trial location for a next sample 
    INTEGER                :: trial          ! loop counter for trial points 
    INTEGER                :: subcell        ! current subcell 
    INTEGER                :: subIndex       ! octal's child index for current subcell 
    INTEGER                :: i


    ! find which of the subcells the point lies in
    DO 
       if (octree%twod) then
          rotloc = projecttoxz(locator)
       else
          rotloc = locator
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
          do; enddo
        ENDIF

        CALL returnSamples(currentPoint,startPoint,locator,direction,        &
                    octree%child(subIndex),grid,sampleFreq,nSamples,         &
                    maxSamples,abortRay,lambda,usePops,iLambda,error,margin, &
                    distanceLimit,kappaAbs=kappaAbs,kappaSca=kappaSca,       &
                    velocity=velocity,velocityDeriv=velocityDeriv,           &
                    chiLine=chiLine,levelPop=levelPop,rho=rho,               &
                    temperature=temperature, Ne=Ne, inFlow=inFlow)
        
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


       call getExitPoint2(currentPoint, direction ,locator, abortRay, error, &
       grid%halfSmallestSubcell, exitPoint, minWallDistance, grid, octree)


!        CALL getExitPoint(currentPoint,direction,locator,abortRay,error, &
!                          grid%halfSmallestSubcell,exitPoint,centre,subcellSize,&
!                          minWallDistance,margin,grid%octreeRoot%threed, octree, subcell)
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
                          kappaAbs=kappaAbs,kappaSca=kappaSca, velocity=velocity,&
                          velocityDeriv=velocityDeriv,chiLine=chiLine,          &
                          levelPop=levelPop,rho=rho,temperature=temperature,    &
                          Ne=Ne, inFlow=inFlow)
        !END IF

        ! we add sampleLength to the distance from the last location
        !   that was sampled, and decide whether to take a new sample.
        trial = 1
        DO 
          trialPoint = currentPoint + direction * &
                       ((REAL(lambda(nSamples),kind=oct) + &
                       (REAL(trial,kind=oct) * sampleLength)) - length)
                       
          ! we only want to take a sample if we are still within the subcell
          IF ( modulus(trialPoint - currentPoint) < minWallDistance ) THEN
          
            trialLength = length + modulus(trialPoint - currentPoint)
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
                         Ne=Ne, inFlow=inFlow)
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
    type(OCTALVECTOR), intent(in) :: currentPoint
    type(OCTALVECTOR), intent(in) :: direction
    type(OCTALVECTOR), intent(out) :: locator    
    logical, intent(inout) :: abortRay
    type(OCTAL), target  :: thisOctal
    integer, intent(inout) :: error
    real(oct), INTENT(IN)    :: halfSmallestSubcell
    type(OCTALVECTOR), intent(out) :: exitPoint
    real(oct), intent(out) :: minWallDistance
    type(OCTAL), pointer :: sOctal
    logical, optional :: edgeOfGrid
    
    integer :: subcell

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
          
    TYPE(octalVector), INTENT(IN)       :: currentPoint ! current ray position
    TYPE(octalVector), INTENT(IN)       :: direction    ! ray's direction vector 
    TYPE(octalVector), INTENT(OUT)      :: locator       
                  ! 'locator' is used to indicate a point that lies within the  
                  !   *next* cell of the octree that the ray will interesect.
    LOGICAL, INTENT(INOUT)              :: abortRay     ! flag to signal completion
    INTEGER, INTENT(INOUT)              :: error        ! error code
    real(oct), INTENT(IN)    :: halfSmallestSubcell
    real(oct), INTENT(IN)    :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    TYPE(octalVector), INTENT(OUT)      :: exitPoint    ! where ray leaves current cell
    TYPE(octalVector), INTENT(IN)       :: centre       ! centre of current subcell
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
    TYPE(octalVector)      :: wallNormal     ! normal to plane of cell wall

    INTEGER, parameter :: max_num_err = 10;
    INTEGER, save :: num_err = 0    

    LOGICAL :: threed
    REAL(oct) :: r2, d, cosMu, disttor2, disttoxboundary,disttozboundary
    REAL(oct) :: compZ, currentZ, tval, x1, x2, currentx, compx
    REAL(oct) :: r1, theta, mu, disttor1

    TYPE(octalvector) :: xDir, zDir
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
      wallNormal =  xHatOctal
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
        locator = exitpoint + (halfSmallestSubcell*xHatOctal)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*xHatOctal)
      ENDIF

    ELSEIF ( wallDistanceY < wallDistanceX .AND. &
             wallDistanceY < wallDistanceZ ) THEN
      wallNormal =  yHatOctal
      IF ( direction%y > 0.0_oc ) THEN
        wallFromOrigin = centre%y + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*yHatOctal)
      ELSE
        wallFromOrigin = centre%y - (subcellSize / 2.0_oc) 
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*yHatOctal)
      ENDIF
      
    ELSEIF ( wallDistanceZ < wallDistanceX .AND. &
             wallDistanceZ < wallDistanceY ) THEN
      wallNormal =  zHatOctal
      IF ( direction%z > 0.0_oc ) THEN
       wallFromOrigin = centre%z + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*zHatOctal)
      ELSE
        wallFromOrigin = centre%z - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*zHatOctal)
      ENDIF

    ! we now consider the case where the ray is leaving through one
    !   of the vertices of the cell.
    ! we could integrate this into the code above, but it is very
    !   unlikely that it will be executed, so we avoid evaluating a 
    !   few IFs by keeping it separate.
    ELSEIF ( wallDistanceX == wallDistanceY .AND. &
             wallDistanceX /= wallDistanceZ ) THEN
      wallNormal =  xHatOctal
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHatOctal)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*xHatOctal)
      ENDIF
      IF ( direction%y > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*yHatOctal)
      ELSE  
        locator = locator - (halfSmallestSubcell*yHatOctal)
      ENDIF  
    ELSEIF ( wallDistanceX == wallDistanceZ .AND. &
             wallDistanceX /= wallDistanceY ) THEN
      wallNormal =  xHatOctal
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHatOctal)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
      IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*xHatOctal)
      ENDIF
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHatOctal)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHatOctal)
      ENDIF  

    ELSEIF ( wallDistanceY == wallDistanceZ  .AND. &
             wallDistanceX /= wallDistanceZ ) THEN
      wallNormal =  yHatOctal
      IF ( direction%y > 0.0_oc ) THEN
        wallFromOrigin = centre%y + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*yHatOctal)
      ELSE
        wallFromOrigin = centre%y - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*yHatOctal)
      ENDIF
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHatOctal)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHatOctal)
      ENDIF  
      
    ! we now consider the case where the ray is leaving through one
    !   of the corners of the cell.
         
    ELSEIF ( wallDistanceX == wallDistanceY  .AND. &
             wallDistanceY == wallDistanceZ ) THEN
      wallNormal =  xHatOctal
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                   (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHatOctal)
      ELSE
        wallFromOrigin = centre%x - (subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint - (halfSmallestSubcell*xHatOctal)
      ENDIF
      IF ( direction%y > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*yHatOctal)
      ELSE  
        locator = locator - (halfSmallestSubcell*yHatOctal)
      ENDIF  
      IF ( direction%z > 0.0_oc ) THEN             
        locator = locator + (halfSmallestSubcell*zHatOctal)
      ELSE  
        locator = locator - (halfSmallestSubcell*zHatOctal)
      ENDIF  

    ! we should now have run out of possibilities
    ELSE
      PRINT *, 'Panic: Fell through direction finder!'
      DO ; END DO ; STOP ! sit in loop for debugging purposes
    END IF

    else       

       ! two-d case written by TJH on 27/1/05
       ! this code isn't ugly and it does work !!!!

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

       r1 = centre%x - subcellSize/2.d0
       r2 = centre%x + subcellSize/2.d0
       d = sqrt(currentpoint%x**2+currentpoint%y**2)
       xDir = OCTALVECTOR(currentpoint%x, currentpoint%y,0.d0)


       if ((direction%x**2+direction%y**2) /= 0.0) then
          if (modulus(xDir) /= 0.d0) then
             call normalize(xDir)       
             cosmu =((-1.d0)*xDir).dot.direction
             call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
             if (.not.ok) then
                write(*,*) "Error:: Quad solver failed in [amr_mod::intersectcubeamr2d] (1)"
                write(*,*) " Correctionn forced for now. If problem insistet, you need to fix it."
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
                call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
                if (.not.ok) then
                   write(*,*) "Quad solver failed in  [amr_mod::intersectcubeamr2d] (2)"
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

       zDir = OCTALVECTOR(0.d0, 0.d0, 1.d0)
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
          write(*,*) "tVal,compX,compZ, distToZboundary,disttoxboundary = "
          write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
          write(*,*) "x,z = ",currentX,currentZ
          stop 
       elseif (tval <= 0.) then
          write(*,*) "Error :: tVal <= 0  [amr_mod:getExitPoint]."
          write(*,*) "tVal,compX, compZ, distToZboundary,disttoxboundary = "
          write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
          write(*,*) "x,z = ",currentX,currentZ
          stop
       endif
   
       minWallDistance = tVal
   
       exitPoint = currentPoint + tval * direction
       locator = exitPoint + frac*halfSmallestSubcell * direction

!       write(*,*) tval, exitpoint

    endif

  END SUBROUTINE getExitPoint 



  SUBROUTINE takeSample(point,length,direction,grid,thisOctal,subcell,nSamples,&
                        maxSamples,usePops,iLambda,error,lambda,kappaAbs,     &
                        kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,  &
                        temperature,Ne,inFlow) 
  
    
    IMPLICIT NONE
    
    TYPE(octalVector), INTENT(IN)      :: point     ! place to make sample
    real(oct), INTENT(IN)   :: length    ! 
    TYPE(octalVector), INTENT(IN)      :: direction ! direction vector
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

    TYPE(vector)                       :: directionReal ! direction vector (REAL values)
    TYPE(octal), POINTER               :: localPointer  ! pointer to the current octal
 
    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    REAL :: Vr, r, Rs
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn)

    nSamples = nSamples + 1
    IF (nSamples > maxSamples) THEN
      PRINT *, "Erro:: nSamples > maxSamples in takeSample subroutine"
      PRINT *, "nSamples   = ", nSamples
      PRINT *, "maxSamples = ", maxSamples      
      write(*,*) lambda(nSamples-10:nSamples)
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
                         inFlow=inFlow(nSamples)                       &
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
    TYPE(octalVector), INTENT(IN)     :: point
    TYPE(octalVector)                 :: point2 ! may be projected point
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
    REAL(double),INTENT(OUT),OPTIONAL         :: kappaAbs
    REAL(double),INTENT(OUT),OPTIONAL         :: kappaSca
    REAL(double),INTENT(OUT),OPTIONAL         :: kappaAbsArray(:)
    REAL(double),INTENT(OUT),OPTIONAL         :: kappaScaArray(:)
    REAL(double),INTENT(OUT), OPTIONAL        :: rosselandKappa
    REAL,INTENT(OUT), OPTIONAL        :: kappap
    REAL,INTENT(IN), OPTIONAL         :: atThisTemperature
    REAL(double),INTENT(OUT),OPTIONAL         :: rho
    REAL,INTENT(OUT),OPTIONAL         :: chiLine
    REAL,INTENT(OUT),OPTIONAL         :: etaLine
    REAL,INTENT(OUT),OPTIONAL         :: etaCont
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

! this for possibility of twoD AMR grid

    point2 = point

    if (octaltree%twoD) then
       point2 = projectToXZ(point)
    endif


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
        CALL findSubcellLocal(point2,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
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
      IF (PRESENT(velocity))  velocity =amrGridVelocity(octalTree,point,foundOctal=resultOctal,&
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
      IF (PRESENT(Ne))                     Ne = resultOctal%Ne(subcell)
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
                                        foundSubcell,actualSubcell) 
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
    TYPE(octalVector), INTENT(IN)  :: point
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell
    INTEGER, INTENT(IN),  OPTIONAL :: actualSubcell

    TYPE(octal), POINTER           :: resultOctal
    INTEGER                        :: subcell

    TYPE(octalVector)              :: centre
    real(oct)           :: inc
    real(oct)           :: t1, t2, t3
    real :: phi
    type(vector) :: newvec
    TYPE(octalVector) :: point_local

    if (octalTree%threeD) then
       point_local = point
    elseif (octalTree%twoD) then
       ! roate "point" back to z-x plane!
       point_local = projectToXZ(point)
    else 
       ! assume it's threeD for now
       point_local = point
    end if

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


      inc = resultOctal%subcellSize / 2.0
      centre = subcellCentre(resultOctal,subcell)
      
      t1 = MAX(0.0_oc, (point_local%x - (centre%x - inc)) / resultOctal%subcellSize)
      t2 = MAX(0.0_oc, (point_local%y - (centre%y - inc)) / resultOctal%subcellSize)
      t3 = MAX(0.0_oc, (point_local%z - (centre%z - inc)) / resultOctal%subcellSize)


      if (resultOctal%threed) then
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
    TYPE(octalVector), INTENT(IN)  :: position
    TYPE(vector), INTENT(IN)       :: direction
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    TYPE(octal),  POINTER         :: thisOctal, currentOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell

    TYPE(octalVector)              :: octalDirection
    TYPE(octal), POINTER           :: firstOctal
    real(oct)           :: dr, dx, dphi
    real(oct)           :: r
    real(oct)           :: phi1, phi2
    TYPE(octalVector)              :: position1
    TYPE(octalVector)              :: position2
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

    
  FUNCTION amrGridKappaAbs(octalTree,point,iLambda,startOctal,&
                           foundOctal,foundSubcell,actualSubcell) 
    ! POINT --> can be in both in roteated or unrotated coordinates for 2D case 
    !
    ! returns the continuous absorption opacity at a given point in the grid.
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

    REAL                           :: amrGridKappaAbs
    TYPE(octal), POINTER           :: octalTree
    TYPE(octalVector), INTENT(IN)  :: point
    INTEGER, INTENT(IN)            :: iLambda
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell
    INTEGER, INTENT(IN),  OPTIONAL :: actualSubcell

    TYPE(octal), POINTER           :: resultOctal
    INTEGER                        :: subcell
    TYPE(octalVector) :: point_local

    if (octalTree%threeD) then
       point_local = point
    elseif (octalTree%twoD) then
       ! roate "point" back to z-x plane!
       point_local = projectToXZ(point)
    else 
       ! assume it's threeD for now
       point_local = point
    end if

    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point_local,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      
      ! should add in some interpolation routine
      amrGridKappaAbs = startOctal%kappaAbs(subcell,iLambda)
      
    ELSE
       ! rotated point used for 2d case
      CALL findSubcellTD(point_local,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal))   foundOctal   => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

      ! should add in some interpolation routine
      amrGridKappaAbs = resultOctal%kappaAbs(subcell,iLambda)

    END IF
  
  END FUNCTION amrGridKappaAbs

  FUNCTION amrGridKappaSca(octalTree,point,iLambda,startOctal,&
                           foundOctal,foundSubcell,actualSubcell) 
    ! POINT --> can be in both in roteated or unrotated coordinates for 2D case 
    !

    ! returns the scattering opacity at a given point in the grid.
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

    REAL                           :: amrGridKappaSca
    TYPE(octal), POINTER           :: octalTree
    INTEGER, INTENT(IN)            :: iLambda
    TYPE(octalVector), INTENT(IN)  :: point
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell
    INTEGER, INTENT(IN),  OPTIONAL :: actualSubcell

    TYPE(octal), POINTER           :: resultOctal
    INTEGER                        :: subcell
    TYPE(octalVector) :: point_local

    if (octalTree%threeD) then
       point_local = point
    elseif (octalTree%twoD) then
       ! roate "point" back to z-x plane!
       point_local = projectToXZ(point)
    else 
       ! assume it's threeD for now
       point_local = point
    end if

    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point_local,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      
      ! should add in some interpolation routine
      amrGridKappaSca = startOctal%kappaSca(subcell,iLambda)
      
    ELSE
      ! rotated point for 2d case
      CALL findSubcellTD(point_local,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal))   foundOctal   => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

      ! should add in some interpolation routine
      amrGridKappaSca = resultOctal%kappaSca(subcell,iLambda)

    END IF
  
  END FUNCTION amrGridKappaSca

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
    TYPE(octalVector), INTENT(IN)  :: point
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell
    INTEGER, INTENT(IN),  OPTIONAL :: actualSubcell

    TYPE(octal), POINTER           :: resultOctal
    INTEGER                        :: subcell
    TYPE(octalVector) :: point_local

    if (octalTree%threeD) then
       point_local = point
    elseif (octalTree%twoD) then
       ! roate "point" back to z-x plane!
       point_local = projectToXZ(point)
    else 
       ! assume it's threeD for now
       point_local = point
    end if

    
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

  FUNCTION amrGridDensity(octalTree,point,startOctal,foundOctal,& 
                                        foundSubcell,actualSubcell) 
    ! POINT --> can be in both in roteated or unrotated coordinates for 2D case 
    !

    ! returns the density at a given point in the grid.
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

    REAL                           :: amrGridDensity
    TYPE(octal), POINTER           :: octalTree
    TYPE(octalVector), INTENT(IN)  :: point
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell
    INTEGER, INTENT(IN),  OPTIONAL :: actualSubcell

    TYPE(octal), POINTER           :: resultOctal
    INTEGER                        :: subcell
    TYPE(octalVector) :: point_local

    if (octalTree%threeD) then
       point_local = point
    elseif (octalTree%twoD) then
       ! roate "point" back to z-x plane!
       point_local = projectToXZ(point)
    else 
       ! assume it's threeD for now
       point_local = point
    end if

    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point_local,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      
      ! should add in some interpolation routine
      amrGridDensity = startOctal%rho(subcell)
      
    ELSE
      ! rotated point for 2d case
      CALL findSubcellTD(point_local,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal)) foundOctal => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

      ! should add in some interpolation routine
      amrGridDensity = resultOctal%rho(subcell)

    END IF

  END FUNCTION amrGridDensity

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
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(octalVector)             :: rotpoint
    INTEGER                       :: subcell
    real(double) :: r, phi
    integer :: i

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

             if (phi < thisOctal%phi) then
                IF ( r < thisOctal%r) THEN
                   IF ( point%z <= thisOctal%centre%z ) THEN
                      subcell = 1
                   ELSE 
                      subcell = 3
                   ENDIF
                ELSE
                   IF (point%z < thisOctal%centre%z) THEN
                      subcell = 2
                   ELSE 
                      subcell = 4
                   ENDIF
                END IF
             else
                IF ( r < thisOctal%r ) THEN
                   IF ( point%z < thisOctal%centre%z ) THEN
                      subcell = 5
                   ELSE 
                      subcell = 7
                   ENDIF
                ELSE
                   IF (point%z < thisOctal%centre%z) THEN
                      subcell = 6
                   ELSE 
                      subcell = 8
                   ENDIF
                END IF
             endif
          else

             IF ( r < thisOctal%r ) THEN
                IF ( point%z < thisOctal%centre%z ) THEN
                   subcell = 1
                ELSE 
                   subcell = 3
                ENDIF
                ELSE
                   IF (point%z < thisOctal%centre%z) THEN
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

  END FUNCTION whichSubcell    


  FUNCTION inOctal(thisOctal,point) 
    ! true if the point lies within the boundaries of the current octal
  
    IMPLICIT NONE
 
    LOGICAL                       :: inOctal
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(octalVector)             :: octVec2D
    real(double)                  :: r, phi, r1, dphi, eps
    

    if (thisOctal%threeD) then
       if (.not.thisOctal%cylindrical) then
          IF (point%x < thisOctal%centre%x - thisOctal%subcellSize ) THEN ; inOctal = .FALSE. 
          ELSEIF (point%x > thisOctal%centre%x + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
          ELSEIF (point%y < thisOctal%centre%y - thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
          ELSEIF (point%y > thisOctal%centre%y + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
          ELSEIF (point%z < thisOctal%centre%z - thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
          ELSEIF (point%z > thisOctal%centre%z + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
          ELSE  
             inOctal = .TRUE.
          ENDIF
       else
          phi = atan2(point%y,point%x)
          if (phi < 0.d0) phi = phi + twoPi
          dphi = abs(phi - thisOctal%phi)
          r = sqrt(point%x**2 + point%y**2)
          eps =  0.d0 ! sqrt(epsilon(thisOctal%subcellsize))
          IF     (r < thisOctal%r - thisOctal%subcellSize - eps) THEN ; inOctal = .FALSE. 
          ELSEIF (r > thisOctal%r + thisOctal%subcellSize + eps) THEN ; inOctal = .FALSE.
          ELSEIF (dphi > thisOctal%dphi/2.d0) THEN ; inOctal = .FALSE.
          ELSEIF (point%z < thisOctal%centre%z - thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
          ELSEIF (point%z > thisOctal%centre%z + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
          ELSE  
             inOctal = .TRUE.
          ENDIF
       endif
    else ! twoD case
       octVec2D = projectToXZ(point)
           IF (octVec2D%x <  thisOctal%centre%x - thisOctal%subcellSize ) THEN ; inOctal = .FALSE. 
       ELSEIF (octVec2D%x > thisOctal%centre%x + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSEIF (octVec2D%z < thisOctal%centre%z - thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSEIF (octVec2D%z > thisOctal%centre%z + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSE  
          inOctal = .TRUE.
       ENDIF
    endif
  END FUNCTION inOctal

  FUNCTION inSubcell(thisOctal,thisSubcell,point) 
    ! true if the point lies within the boundaries of the current octal
  
    IMPLICIT NONE
 
    LOGICAL                       :: inSubcell
    TYPE(octal), INTENT(IN)       :: thisOctal
    INTEGER, INTENT(IN)           :: thisSubcell
    TYPE(octalVector), INTENT(IN) :: point

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
    TYPE(octalVector), INTENT(IN) :: point

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

  SUBROUTINE smoothAMRgrid(grid,factor,sphData, stellar_cluster, inheritProps, interpProps)
    ! checks whether each octal's neighbours are much bigger than it, 
    !   if so, makes the neighbours smaller.

    TYPE(gridtype), INTENT(INOUT) :: grid 
    REAL, INTENT(IN)              :: factor
    TYPE(sph_data), optional, INTENT(IN) :: sphData   ! Matthew's SPH data.
    TYPE(cluster), optional, intent(in)  :: stellar_cluster
    LOGICAL, INTENT(IN),optional  :: inheritProps
    logical, intent(in), optional :: interpProps 
    
    LOGICAL :: gridConverged

    CALL setAllUnchanged(grid%octreeRoot)

    DO 
      gridConverged = .TRUE.
      CALL smoothAMRgridPrivate(grid%octreeRoot,grid,gridConverged)
      IF ( gridConverged ) EXIT 
    END DO
    
  CONTAINS
    
    RECURSIVE SUBROUTINE smoothAMRgridPrivate(thisOctal,grid,gridConverged)

      TYPE(octal), INTENT(INOUT), TARGET :: thisOctal
      TYPE(gridtype), INTENT(INOUT   ) :: grid 
      LOGICAL, INTENT(INOUT)               :: gridConverged

      INTEGER              :: i
      REAL(oct) :: halfSmallestSubcell
      REAL(oct) :: offset
      TYPE(octal), POINTER :: neighbour
      TYPE(octalVector), ALLOCATABLE, DIMENSION(:) :: locator
      TYPE(octalVector) :: aHat
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
        locator(:) = thisOctal%centre + ( offset * octalVector(1.0_oc,1.0_oc,1.0_oc) )
      ELSE
        locator(:) = thisOctal%centre + ( offset * octalVector(1.0_oc,0.0_oc,1.0_oc) )
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
            locator(1) = locator(1) + (thisOctal%subcellSize + halfSmallestSubcell) * zHatOctal
            locator(2) = locator(2) - (thisOctal%subcellSize + halfSmallestSubcell) * zHatOctal
            aHat = OCTALVECTOR(thisOctal%centre%x,thisOctal%centre%y,0.d0)
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
                               sphData=sphData, stellar_cluster=stellar_cluster, &
                               inherit=inheritProps, interp=interpProps)
              gridConverged = .FALSE.
          ENDIF
        END IF
          
! force return until algorithm is fixed
IF ( .NOT. gridConverged ) RETURN
          
      END DO

      ! call this subroutine recursively on each of any children.
      IF ( (.NOT. ANY(thisOctal%Changed)) .AND. thisOctal%nChildren > 0 ) THEN
        DO i = 1, thisOctal%nChildren, 1 
          CALL smoothAMRgridPrivate(thisOctal%child(i),grid,gridConverged)
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

    TYPE(octalVector), INTENT(IN) :: point
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


  SUBROUTINE findSubcellLocal(point,thisOctal,subcell, prob)
    ! finds the octal (and that octal's subcell) containing a point.
    !   starts searching from the current octal, and goes up and down the
    !   tree as needed to find the correct octal.

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(octalVector) :: point_local
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
       point_local = projectToXZ(point)
    else
       point_local = point
    endif


    CALL findSubcellLocalPrivate(point_local,thisOctal,subcell,&
                                 haveDescended,boundaryProblem)
    if (present(prob)) prob = boundaryProblem
                                 
  CONTAINS

    RECURSIVE SUBROUTINE findSubcellLocalPrivate(point,thisOctal,subcell,&
                                                 haveDescended,boundaryProblem)
      TYPE(octalVector), INTENT(IN) :: point
      TYPE(octal),POINTER    :: thisOctal
      INTEGER, INTENT(OUT)   :: subcell
      LOGICAL, INTENT(INOUT) :: haveDescended
      LOGICAL, INTENT(INOUT) :: boundaryProblem
      type(octalvector) :: rVec
      INTEGER :: i
      
      IF ( inOctal(thisOctal,point) ) THEN

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
           do ; enddo
        endif
        
        IF ( thisOctal%nDepth /= 1 ) THEN
           thisOctal => thisOctal%parent
        ENDIF

        CALL findSubcellLocalPrivate(point,thisOctal,subcell,haveDescended,boundaryProblem)
       
      END IF    
    
    END SUBROUTINE findSubcellLocalPrivate

  END SUBROUTINE findSubcellLocal
   

  recursive subroutine amrPlotGrid(thisOctal)
    ! plots the outline of each of the cells in the grid

    IMPLICIT NONE
    type(OCTAL), POINTER :: thisOctal
    type(OCTAL), POINTER :: child
    real :: x1, y1
    integer :: j


       x1 = REAL(thisOctal%centre%x-thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%z-thisOctal%subcellSize)
       call pgmove(x1,y1)
       x1 = REAL(thisOctal%centre%x+thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%z-thisOctal%subcellSize)
       call pgdraw(x1,y1)
       x1 = REAL(thisOctal%centre%x+thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%z+thisOctal%subcellSize)
       call pgdraw(x1,y1)
       x1 = REAL(thisOctal%centre%x-thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%z+thisOctal%subcellSize)
       call pgdraw(x1,y1)
       x1 = REAL(thisOctal%centre%x-thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%z-thisOctal%subcellSize)
       call pgdraw(x1,y1)

       x1 = REAL(thisOctal%centre%x)
       y1 = REAL(thisOctal%centre%z-thisOctal%subcellSize)
       call pgmove(x1,y1)
       x1 = REAL(thisOctal%centre%x)
       y1 = REAL(thisOctal%centre%z+thisOctal%subcellSize)
       call pgdraw(x1,y1)
       
       x1 = REAL(thisOctal%centre%x-thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%z)
       call pgmove(x1,y1)
       x1 = REAL(thisOctal%centre%x+thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%z)
       call pgdraw(x1,y1)
       
    do j = 1, thisOctal%nChildren
       child => thisOctal%child(j)
       call amrPlotGrid(child)
    enddo


  end subroutine amrPlotGrid

  RECURSIVE SUBROUTINE plotPointGridXY(thisOctal)
    ! recursively plots a point at the centre of each subcell which
    !   does not have a child

    IMPLICIT NONE
    TYPE(octal), POINTER :: thisOctal
    TYPE(octal), POINTER :: child
    INTEGER :: sub, j
    TYPE(octalVector) :: centre

    DO sub = 1, thisOctal%maxChildren, 1
      IF ( thisOctal%hasChild(sub) ) THEN
        ! if the current subcell has a child, call this subroutine on
        !   each of its children
        DO j = 1, thisOctal%nChildren
          IF ( thisOctal%indexChild(j) == sub ) THEN
            child => thisOctal%child(j)
            CALL plotPointGridXY(child)
            EXIT
          END IF
        END DO
      ELSE
        ! if the current subcell does not have a child, plot a point
        centre = subcellCentre(thisOctal,sub)
        CALL pgpt1(REAL(centre%x),REAL(centre%y),-1)
      END IF
    END DO

  END SUBROUTINE plotPointGridXY
  
  RECURSIVE SUBROUTINE plotPointGridXZ(thisOctal)
    ! recursively plots a point at the centre of each subcell which
    !   does not have a child

    IMPLICIT NONE
    TYPE(octal), POINTER :: thisOctal
    TYPE(octal), POINTER :: child
    INTEGER :: sub, j
    TYPE(octalVector) :: centre


    DO sub = 1, thisOctal%maxChildren, 1
      IF ( thisOctal%hasChild(sub) ) THEN
        ! if the current subcell has a child, call this subroutine on
        !   each of its children
        DO j = 1, thisOctal%nChildren
          IF ( thisOctal%indexChild(j) == sub ) THEN
            child => thisOctal%child(j)
            CALL plotPointGridXZ(child)
            EXIT
          END IF
        END DO
      ELSE
        ! if the current subcell does not have a child, plot a point
        centre = subcellCentre(thisOctal,sub)
        CALL pgpt1(REAL(centre%x),REAL(centre%z),-1)
      END IF
    END DO

  END SUBROUTINE plotPointGridXZ

  RECURSIVE SUBROUTINE plotMiddlePointGridXZ(thisOctal,thickness)
    ! recursively plots a point at the centre of each subcell which
    !   does not have a child
    ! ignores any cell which does not lie close to the central plane of
    !   the simulation space. use this to look at a crosscut through 
    !   the middle of the simulation

    IMPLICIT NONE
    TYPE(octal), POINTER :: thisOctal
    TYPE(octal), POINTER :: child
    REAL :: thickness ! maximum distance from centre to plot points
    INTEGER :: sub, j
    TYPE(octalVector) :: centre

    DO sub = 1, thisOctal%maxChildren, 1
      IF ( thisOctal%hasChild(sub) ) THEN
        DO j = 1, thisOctal%nChildren
          IF ( thisOctal%indexChild(j) == sub ) THEN
            child => thisOctal%child(j)
            CALL plotMiddlePointGridXZ(child,thickness)
            EXIT
          END IF
        END DO
      ELSE
        centre = subcellCentre(thisOctal,sub)
        IF ( ABS(centre%y) < thickness ) THEN
          CALL pgpt1(REAL(centre%x),REAL(centre%z),-1)
        END IF
      END IF
    END DO

  END SUBROUTINE plotMiddlePointGridXZ

  RECURSIVE SUBROUTINE plotMiddlePointGridYZ(thisOctal,thickness)
    ! recursively plots a point at the centre of each subcell which
    !   does not have a child
    ! ignores any cell which does not lie close to the central plane of
    !   the simulation space. use this to look at a crosscut through 
    !   the middle of the simulation

    IMPLICIT NONE
    TYPE(octal), POINTER :: thisOctal
    TYPE(octal), POINTER :: child
    REAL :: thickness ! maximum distance from centre to plot points
    INTEGER :: sub, j
    TYPE(octalVector) :: centre

    DO sub = 1, thisOctal%maxChildren, 1
      IF ( thisOctal%hasChild(sub) ) THEN
        DO j = 1, thisOctal%nChildren
          IF ( thisOctal%indexChild(j) == sub ) THEN
            child => thisOctal%child(j)
            CALL plotMiddlePointGridYZ(child,thickness)
            EXIT
          END IF
        END DO
      ELSE
        centre = subcellCentre(thisOctal,sub)
        IF ( ABS(centre%x) < thickness ) THEN
          CALL pgpt1(REAL(centre%y),REAL(centre%z),-1)
        END IF
      END IF
    END DO

  END SUBROUTINE plotMiddlePointGridYZ

  subroutine plotAMRthreeDMovie(device, thisType, grid)
    character(len=*) :: device,thisType
    character(len=30) :: filename
    type(GRIDTYPE) :: grid
    integer :: nFrames, maxFrames
    type(VECTOR),allocatable :: viewVec(:), posVec(:)
    type(VECTOR) :: rHat
    real :: defDist, fac
    type(OCTALVECTOR) :: octViewVec,octPosVec
    real :: ang
    real(double) :: totalMass, minRho, maxRho
    integer maxPoints
    integer :: nPoints
    type(VECTOR),allocatable :: points(:), lineGrid(:,:)
    real :: x, y, dr
    integer :: nGrid
    integer, allocatable :: nColour(:)
    real, allocatable :: depth(:)
    integer :: i
    real :: dist, r1, r2
    type(VECTOR) :: newPoint, oldPoint, thisPoint, direction,campoint(4),startpoint
    integer :: nScene = 100,  j
   
    defDist = 4.*grid%octreeRoot%subcellSize
    
    minRho = 1.e30
    maxRho = -1.e30
    nGrid = 0
    allocate(lineGrid(1:10000,1:4))
    r1 = 2.*grid%ocTreeRoot%subcellSize
    dr = r1/9.
    do i = 1, 10
       do j = 1, 10
          x = (real(i-1)/9.-0.5)*r1
          y = (real(j-1)/9.-0.5)*r1
          nGrid = nGrid + 1
          lineGrid(nGrid,1) = VECTOR(x - dr, y - dr, 0.)
          lineGrid(nGrid,2) = VECTOR(x + dr, y - dr, 0.)
          lineGrid(nGrid,3) = VECTOR(x + dr, y + dr, 0.)
          lineGrid(nGrid,4) = VECTOR(x - dr, y + dr, 0.)
       enddo
    enddo
    nGrid = 100
    allocate(lineGrid(1:ngrid,1:2))
    r1 = 100.*grid%ocTreeRoot%subcellSize
    do i = 1, nGrid/2
       r2 = r1 * (2.*real(i-1)/real(nGrid/2-1)-1.)
       lineGrid(i,1) = VECTOR(-r1,r2,0.)
       lineGrid(i,2) = VECTOR(+r1,r2,0.)
    enddo
    do i = ngrid/2+1, nGrid
       r2 = r1 * (2.*real(i-ngrid/2-1)/real(nGrid/2-1)-1.)
       lineGrid(i,1) = VECTOR(r2,-r1,0.)
       lineGrid(i,2) = VECTOR(r2,+r1,0.)
    enddo
    maxFrames = 1000

    allocate(viewVec(1:maxFrames))
    allocate(posVec(1:maxFrames))
    allocate(depth(1:maxFrames))
    nFrames = 0

    oldPoint = randomUnitVector()
    oldPoint = 2.*defDist*oldPoint
    newPoint = VECTOR(0.0001, 0.0001, 0.01 * defDist)

    dist = modulus(newPoint - oldPoint)
    direction = (newPoint - oldPoint) / dist
    rHat = randomUnitVector()
    nFrames = 0
    nScene = 100

    campoint(1) = defDist*randomUnitVector()
    campoint(2) = defDist*randomUnitVector()
    campoint(3) = defDist*randomUnitVector()
    campoint(4) = defDist*randomUnitVector()
    call cameraPositions(campoint, 4, posVec, viewVec, nScene, nFrames)
    campoint(1) = camPoint(4)
    campoint(2) = defDist*randomUnitVector()
    campoint(3) = defDist*randomUnitVector()
    campoint(4) = defDist*randomUnitVector()
    call cameraPositions(campoint, 4, posVec, viewVec, nScene, nFrames)
    campoint(1) = camPoint(4)
    campoint(2) = defDist*randomUnitVector()
    campoint(3) = defDist*randomUnitVector()
    campoint(4) = defDist*randomUnitVector()
    call cameraPositions(campoint, 4, posVec, viewVec, nScene, nFrames)
    campoint(1) = camPoint(4)
    campoint(2) = defDist*randomUnitVector()
    campoint(3) = defDist*randomUnitVector()
    campoint(4) = defDist*randomUnitVector()
    call cameraPositions(campoint, 4, posVec, viewVec, nScene, nFrames)
    

    totalMass =0.
    minrho=1.e30
    maxrho=1.e-30
    call findTotalMass(grid%octreeRoot, totalMass, minRho, maxRho)
!    minrho = maxrho/10000.
    maxPoints = 1000000
    write(*,*) log10(minrho),"->",log10(maxrho)
    allocate(points(1:maxPoints))
    allocate(nColour(1:maxPoints))

    write(*,*) "creating point list"
    npoints = 0

    write(filename,'(a,i3.3,a)') "movie",1,".gif/gif"
    write(*,*) "Writing frame to: ",trim(filename)
    call pgbegin(0,filename,1,1)
    call createPointList(grid%octreeRoot, nPoints, points, nColour, maxPoints, totalMass, minRho, maxRho)
    write(*,*) "done",npoints
    call pgend

    do i = 1, nFrames
       write(filename,'(a,i3.3,a)') "movie",i,".gif/gif"
       write(*,*) "Writing frame to: ",trim(filename)
       call pgbegin(0,filename,1,1)
       call pgpaper(4.,1.)
       call pgwnad(REAL(grid%octreeRoot%centre%x - grid%octreeRoot%subcellSize),&
            REAL(grid%octreeRoot%centre%x + grid%octreeRoot%subcellSize),&
            REAL(grid%octreeRoot%centre%z - grid%octreeRoot%subcellSize),&
            REAL(grid%octreeRoot%centre%z + grid%octreeRoot%subcellSize))
       octViewVec = viewVec(i)
       octPosVec = posVec(i)
       

!       call plotThreed(grid%octreeRoot, octPosVec, octViewVec, depth(i), defDist, points, nPoints)

       call plotPoints(posVec(i),viewVec(i),defDist, points, lineGrid, nGrid, nColour, nPoints)

       call pgsci(1)
       call pgbox('bcnst',0.0,0,'bcnst',0.0,0)
       call pgend

    enddo
    deallocate(points)
  end subroutine plotAMRthreeDMovie


  RECURSIVE SUBROUTINE plotThreeD(thisOctal, posVec,viewVec, thisDepth, defDist, points, nPoints)

    IMPLICIT NONE

    integer :: j
    type(OCTAL), INTENT(IN) :: thisOctal
    type(octalVector) :: viewVec, posVec
    real :: defDist
    real(double) :: totalMass
    integer :: thisDepth
    integer :: nPoints
    type(VECTOR) :: points(:)

    do j = 1, thisOctal%nChildren
       if ((thisOctal%nDepth == thisDepth))then
!          call plotCube(viewVec,subcellCentre(thisOctal,j), &
!                        thisOctal%subcellsize,thisOctal%nDepth)
       endif
    enddo

    do j = 1, thisOctal%maxChildren
!       if ((thisOctal%nDepth == thisDepth) .AND. &
!               (.NOT. (thisOctal%hasChild(j))) )then
       if (.not.thisOctal%hasChild(j).and.(thisOctal%nDepth >= thisDepth).and.(thisOctal%inFlow(j))) then
!          call plotCube(posVec,viewVec,subcellCentre(thisOctal,j), &
!                        thisOctal%subcellsize,thisOctal%nDepth, defDist)


!          call plotPoints(posVec,viewVec,subcellCentre(thisOctal,j), &
!                        thisOctal%subcellsize,thisOctal%nDepth, defDist, points, nPoints)

       endif
    enddo

    do j = 1, thisOctal%nChildren
          call plotThreed(thisOctal%child(j), posVec, viewVec, thisDepth, defDist, points, nPoints)
    enddo


  END SUBROUTINE plotThreeD

  SUBROUTINE plotCube(posVecOctal, viewVecOctal, cubeCentreOctal, cubeSizeOctal, nCol, defDist)
    
    IMPLICIT NONE
    
    type(octalVector) :: viewVecOctal, cubeCentreOctal, posVecOctal
    type(VECTOR) :: cubeCentre,viewVec, posVec
    real(oct) :: cubeSizeOctal
    real :: cubeSize, s, defDist
    type(VECTOR) :: xAxis, yAxis, zAxis, xProj, yProj
    real :: xCube(8), yCube(8)
    type(VECTOR) :: cubeCorner(8)
    integer :: i, nCol
    real :: distFac
    call pgsci(nCol)

    xAxis = VECTOR(1., 0., 0.)
    yAxis = VECTOR(0., 1., 0.)
    zAxis = VECTOR(0., 0., 1.)

    distFac = (cubeCentreOctal - posVecOctal).dot.ViewVecOctal


    if (distFac > 0.) then
       

       cubesize = cubeSizeOctal
       
       viewVec%x = viewVecOctal%x
       viewVec%y = viewVecOctal%y
       viewVec%z = viewVecOctal%z
       
       cubeCentre%x = cubeCentreOctal%x
       cubeCentre%y = cubeCentreOctal%y
       cubeCentre%z = cubeCentreOctal%z
       
       posVec = posVecOctal

       xProj =  zAxis .cross. viewVec
!       xProj = xAxis
       call normalize(xProj)
       yProj = viewVec .cross. xProj
       call normalize(yProj)
       
       s = cubeSize / 2.
       
       cubeCorner(1) = cubeCentre + (s * xAxis) + (s * yAxis) + (s * zAxis)
       cubeCorner(2) = cubeCentre + (s * xAxis) - (s * yAxis) + (s * zAxis)
       cubeCorner(3) = cubeCentre - (s * xAxis) - (s * yAxis) + (s * zAxis)
       cubeCorner(4) = cubeCentre - (s * xAxis) + (s * yAxis) + (s * zAxis)
       cubeCorner(5) = cubeCentre + (s * xAxis) + (s * yAxis) - (s * zAxis)
       cubeCorner(6) = cubeCentre + (s * xAxis) - (s * yAxis) - (s * zAxis)
       cubeCorner(7) = cubeCentre - (s * xAxis) - (s * yAxis) - (s * zAxis)
       cubeCorner(8) = cubeCentre - (s * xAxis) + (s * yAxis) - (s * zAxis)

       do i = 1, 8
          distFac = modulus(cubeCorner(i) - posVec)/defDist
          xCube(i) = ((cubeCorner(i) - posVec).dot. xProj)/distFac
          yCube(i) = ((cubeCorner(i) - posVec).dot. yProj)/distFac
       end do
       
       call pgmove(xCube(1), yCube(1))
       call pgdraw(xCube(2), yCube(2))
       call pgdraw(xCube(3), yCube(3))
       call pgdraw(xCube(4), yCube(4))
       call pgdraw(xCube(1), yCube(1))
       
       
       call pgmove(xCube(5), yCube(5))
       call pgdraw(xCube(6), yCube(6))
       call pgdraw(xCube(7), yCube(7))
       call pgdraw(xCube(8), yCube(8))
       call pgdraw(xCube(5), yCube(5))
       
       call pgmove(xCube(1), yCube(1))
       call pgdraw(xCube(5), yCube(5))
       
       call pgmove(xCube(2), yCube(2))
       call pgdraw(xCube(6), yCube(6))
       
       call pgmove(xCube(3), yCube(3))
       call pgdraw(xCube(7), yCube(7))
       
       call pgmove(xCube(4), yCube(4))
       call pgdraw(xCube(8), yCube(8))
    endif
  END SUBROUTINE plotCube

  SUBROUTINE plotPoints(posVec, viewVec, defDist, points, lineGrid, nGrid, nColour, nPoints)
    
    IMPLICIT NONE
    
    type(VECTOR) :: viewVec, posVec, up
    real :: s, defDist
    type(VECTOR) :: xAxis, yAxis, zAxis, xProj, yProj
    type(VECTOR) :: points(:)
    real, allocatable :: zDepth(:)
    integer, allocatable :: indx(:)
    integer :: nColour(:)
    integer :: i, nPoints, j
    real :: distFac,r1,r2,r3,xs,ys
    integer :: nGrid
    type(VECTOR) :: lineGrid(:,:)

    allocate(zDepth(1:nPoints))
    allocate(indx(1:nPoints))

    do i = 1, nPoints
       zDepth(i) = -((points(i) - posVec).dot.viewVec)
    enddo
    call indexx(nPoints,zDepth,indx)
    

    call palette(3)

    xAxis = VECTOR(1., 0., 0.)
    yAxis = VECTOR(0., 1., 0.)
    zAxis = VECTOR(0., 0., 1.)


    up = zAxis
    yProj = up - (up .dot. viewVec)*viewVec
!    xProj =  zAxis .cross. viewVec
    call normalize(yProj)
    xProj = viewVec .cross. yProj
    call normalize(xProj)
       

    do i = 1, nPoints
       j = indx(i)
       distFac = modulus(points(j) - posVec)/defDist
       if ((viewVec.dot.(points(j) - posVec)) > 0.d0) then
          if ((points(j)%z * viewVec%z) > 0.) then
             xs = ((points(j) - posVec).dot. xProj)/distFac
             ys = ((points(j) - posVec).dot. yProj)/distFac
             call pgsci(nColour(j))
             call pgpoint(1,xs,ys,-1)
          endif
       endif
    end do
    do i = 1, nGrid
       if ((viewVec.dot.(lineGrid(i,1) - posVec)) > 0.d0) then
             distFac = modulus(lineGrid(i,1) - posVec)/defDist
             xs = ((lineGrid(i,1) - posVec).dot. xProj)/distFac
             ys = ((lineGrid(i,1) - posVec).dot. yProj)/distFac
             call pgsci(1)
             call pgmove(xs,ys)
             distFac = modulus(lineGrid(i,2) - posVec)/defDist
             xs = ((lineGrid(i,2) - posVec).dot. xProj)/distFac
             ys = ((lineGrid(i,2) - posVec).dot. yProj)/distFac
             call pgdraw(xs,ys)
             distFac = modulus(lineGrid(i,3) - posVec)/defDist
             xs = ((lineGrid(i,3) - posVec).dot. xProj)/distFac
             ys = ((lineGrid(i,3) - posVec).dot. yProj)/distFac
             call pgdraw(xs,ys)
             distFac = modulus(lineGrid(i,4) - posVec)/defDist
             xs = ((lineGrid(i,4) - posVec).dot. xProj)/distFac
             ys = ((lineGrid(i,4) - posVec).dot. yProj)/distFac
             call pgdraw(xs,ys)
             distFac = modulus(lineGrid(i,1) - posVec)/defDist
             xs = ((lineGrid(i,1) - posVec).dot. xProj)/distFac
             ys = ((lineGrid(i,1) - posVec).dot. yProj)/distFac
             call pgdraw(xs,ys)
       endif
    enddo
    do i = 1, nPoints
       j = indx(i)
       distFac = modulus(points(j) - posVec)/defDist
       if ((viewVec.dot.(points(j) - posVec)) > 0.d0) then
          if ((points(j)%z * viewVec%z) < 0.) then
             xs = ((points(j) - posVec).dot. xProj)/distFac
             ys = ((points(j) - posVec).dot. yProj)/distFac
             call pgsci(nColour(j))
             call pgpoint(1,xs,ys,-1)
          endif
       endif
    end do
    call pgend

  END SUBROUTINE plotPoints



  FUNCTION decideSplit(thisOctal,subcell,amrLimitScalar,amrLimitScalar2,grid, splitInAzimuth,&
       sphData, stellar_cluster) RESULT(split)
    ! returns true if the current voxel is to be subdivided. 
    ! decision is made by comparing 'amrLimitScalar' to some value
    !   derived from information in the current cell  

    use input_variables, only: height, betadisc, rheight, flaringpower, rinner, router
    use input_variables, only: drInner, drOuter, rStellar, cavangle, erInner, erOuter, rCore
    use input_variables, only: warpFracHeight, warpRadius, warpSigma
    IMPLICIT NONE
    TYPE(octal)       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    LOGICAL, INTENT(INOUT) :: splitInAzimuth
    real(double), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 ! used for split decision
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(sph_data), OPTIONAL, intent(in) :: sphData
    TYPE(cluster), OPTIONAL, intent(in) :: stellar_cluster
    LOGICAL                    :: split          
    real(oct)  :: cellSize
    TYPE(octalVector)     :: searchPoint, rVec
    TYPE(octalVector)     :: cellCentre
    REAL                  :: x, y, z
    REAL(double) :: hr, rd, fac, warpHeight, phi
    INTEGER               :: i
    real(double) :: thisTau, kappaSca, kappaAbs
    real(double)      :: total_mass
    real(double)      :: ave_density, rGrid(1000), r, dr
    INTEGER               :: nr, nr1, nr2
    real(double)      :: total_opacity, minDensity, maxDensity, thisDensity
    INTEGER               :: nsample = 400
    LOGICAL               :: inUse
    INTEGER               :: nparticle, limit
    real(double) :: dummyDouble
    real(double) :: fac1,fac2
    real(double) :: rho_disc_ave, scale_length
    real(double),save  :: R_tmp(204)  ! [10^10cm]
    real(double),allocatable, save  :: R_cmfgen(:)  ! [10^10cm]
    real(double),save  :: Rmin_cmfgen  ! [10^10cm]
    real(double) :: rho
    logical, save :: first_time=.true.
    logical :: close_to_star
    real(double)      :: thisScale
    logical,save  :: firstTime = .true.

    splitInAzimuth = .false.
    
   select case(grid%geometry)

    case("melvin")
       split = .false.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      if (thisOctal%nDepth < 4) split = .true.
      if ((cellCentre%x  < 2.e6).and.(cellSize > 1.e5)) split = .true.

    case("whitney")
       split = .false.
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


    case("ttauri","jets","spiralwind")
      nsample = 40
      ! the density is only sampled at the centre of the grid
      ! we will search in each subcell to see if any point exceeds the 
      ! threshold density

      ! get the size and centre of the current cell
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
    
      ! check the density of random points in the current cell - 
      !   if any of them exceed the critical density, set the flag
      !   indicating the cell is to be split, and return from the function
      split = .FALSE.
      ave_density = 0.0_db
      minDensity = 1.e30
      maxDensity = -1.e30
      DO i = 1, nsample
        searchPoint = cellCentre
        CALL RANDOM_NUMBER(x)
        CALL RANDOM_NUMBER(y)
        CALL RANDOM_NUMBER(z)
        if (thisOctal%threed) then
           searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=oct)
           searchPoint%y = searchPoint%y - (cellSize / 2.0_oc) + cellSize*REAL(y,KIND=oct) 
           searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=oct) 
        else
           searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=oct) 
           searchPoint%y = 0.d0
           searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=oct) 
        endif
        ! using a generic density function in density_mod.f90
        rho = density(searchPoint,grid)
        ave_density  =  rho + ave_density
        minDensity = min(minDensity, rho)
        maxDensity = max(maxDensity, rho)
      END DO

      ave_density = ave_density / REAL(nsample,KIND=double)


      if (thisOctal%threed) then
         total_mass = maxDensity * (cellSize*1.e10_db)**3.0_db
      else
         total_mass = maxDensity * pi * ((searchPoint%x+cellsize/2.)**2-(searchPoint%x-cellsize/2.)**2)*cellsize*1.d30
      endif

      if (grid%geometry == "ttauri") then
         thisScale = grid%rstar1/modulus(cellCentre)
         ! weigting toward the smaller radial positions
         total_mass = total_mass * thisScale**4
      endif

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
      end if

   case("lexington")
      if (thisOctal%nDepth < 5) then
         split = .true.
      else
         split = .false.
      endif
      if (modulus(subcellCentre(thisoctal,subcell)) < 6.*grid%rinner) then
         if (thisOctal%nDepth < 6) split = .true.
      endif
      
   case ("testamr","proto","wrshell")
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      split = .FALSE.
      nr1 = 5
      nr2 = 50
      rGrid(1) = 1.
      rGrid(2) = 1.04
      rGrid(3) = 1.08
      rGrid(4) = 1.050
      rGrid(5) = 1.1
      rGrid(1:nr1) = log10(rGrid(1:nr1)*grid%rInner)
      nr = nr1 + nr2
!      do i = 1, nr1
!         rgrid(i) = log10(grid%rInner)+dble(i-1)*(log10(4.*grid%rInner)-log10(grid%rInner))/dble(nr1-1)
!      end do
      do i = 1, nr2
         rgrid(nr1+i) = log10(1.1*grid%rInner)+dble(i)*(log10(grid%rOuter)-log10(1.1*grid%rInner))/dble(nr2)
      end do
      rgrid(1:nr) = 10.d0**rgrid(1:nr)
      r = modulus(cellcentre)
      if (thisOctal%nDepth < 6) split = .true.
      if ((r < grid%rOuter).and.(r > grid%rinner*0.9)) then
         call locate(rGrid, nr, r, i)      
         if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
      endif

   case("benchmark")
      split = .false.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
!      if (thisOctal%nDepth < 6) split = .true.
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      rd = grid%rOuter / 2. 
      hr = height * (r/rd)**1.125

      if ((abs(cellcentre%z)/hr < 10.) .and. (cellsize/hr > 0.2)) split = .true.
      if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 0.2)) split = .true.
      if ((r+cellsize/2.d0) < grid%rinner) split = .false.
      


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

   case ("cluster")
      ! Splits if the number of particle is more than a critical mass.
      ! using the function in amr_mod.f90 
      call find_n_particle_in_subcell(nparticle, ave_density, sphData, &
           thisOctal, subcell)


      ! get the size and centre of the current cell
      cellSize = thisOctal%subcellSize

      !
      total_mass = ave_density * (cellSize*1.e10_db)**3  ! should be in [g]


      if (total_mass > amrlimitscalar .and. nparticle > 0) then
         split = .true.
      else
         split = .false.
      end if

      ! Extra check
      ! if # of SPH particle is greater than 5 then splits...
      if (nparticle> 5) then
         split = .true.
      else
         split = .false.
      end if
!      write(*,*) nparticle,total_mass,amrlimitscalar

      if (include_disc(stellar_cluster)) then
      
         ! If the subcell intersects with the stellar disk.. then do additional check.
         ! This time, the cells will be split when the mass of the cell exceeds a 
         ! critical mass specified by "amrlimitscalar.

         if (stellar_disc_exists(sphData) .and.  &       
              disc_intersects_subcell(stellar_cluster, sphData, thisOctal, subcell) ) then
         rho_disc_ave = average_disc_density_fast(sphData, thisOctal, &
              subcell, stellar_cluster, scale_length)

!!         rho_disc_ave = average_disc_density(sphData, thisOctal, &
!!              subcell, stellar_cluster, scale_length)
!            rho_disc_ave = max_disc_density_from_sample(sphData, thisOctal, &
!                 subcell, stellar_cluster, scale_length)

            total_mass = total_mass + rho_disc_ave * (cellSize*1.e10_db)**3  !  [g]


!            if (cellSize > 1.0e2  .or.  &
            if (cellSize > 1.0e4  .or.  &
                 (cellsize > scale_length) ) then
!                 (cellsize > scale_length .and. rho_disc_ave > 1.0e-20) ) then




               split = .true.  ! units in 10^10cm
            else
               split = .false.
            end if
            
         end if

      end if
            
            
   case ("wr104")
      ! Splits if the number of particle is more than a critical value (~3).
      limit = nint(amrLimitScalar)
      
      call find_n_particle_in_subcell(nparticle, dummyDouble, sphData, thisOctal, subcell)

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

      if ((r < 2.e5).and.(thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 15.)) then
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
      split = .FALSE.

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


   case("shakara","aksco")
      split = .false.
!      if (thisOctal%ndepth  < 6) split = .true.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      hr = height * (r / (100.d0*autocm/1.d10))**betadisc
      if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 1.)) split = .true.
!      if (r < 2.*grid%rInner) then
!         if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 0.2)) split = .true.
!      endif

      if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.




!      if ((r > grid%rInner).and.(r < grid%rInner * 1.01)) then
!         if ((abs(cellcentre%z)/hr < 5.)) then
!            if (cellsize > 1.e-3 * grid%rInner) split = .true.
!         endif
!      endif



      if ((r+cellsize/2.d0) < grid%rinner) split = .false.


   case("ppdisk")
      split = .false.
      cellSize = thisOctal%subcellSize
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
! The 100 is because Tim's code uses the scale height at 100AU.
!      hr = height * (r / (100.d0*autocm/1.d10))**1.25
      hr = (auTocm * height * rHeight * (r / (rHeight*autocm/1.e10))**flaringPower)/1.e10
      if ((abs(cellcentre%z)/hr < 5.) .and. (cellsize/hr > 0.5)) split = .true.
      if ((abs(cellcentre%z)/hr > 5.).and.(abs(cellcentre%z/cellsize) < 2.)) split = .true.
      if ((r+cellsize/2.d0) < grid%rinner) split = .false.



   case("clumpydisc")
      split = .false.
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

   case("toruslogo")
      split = .false.
      if (thisOctal%nDepth  < 6) split = .true.

      if ((thisOctal%cylindrical).and.(thisOctal%dPhi*radtodeg > 30.)) then
         split = .true.
         splitInAzimuth = .true.
      endif


   case("warpeddisc")
      split = .false.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      phi = atan2(cellcentre%y,cellcentre%x)
!      warpheight = 0.3 * rOuter * (r / rOuter)**2 * cos(phi)
      warpheight  = cos(phi) * warpFracHeight * warpradius * exp(-0.5d0*((r - warpRadius)/warpSigma)**2)


      if (r > rInner*0.8) then
         hr = height * (r / rOuter)**betaDisc
         fac = cellsize/hr
         if (abs((cellCentre%z-warpHeight)/hr) < 5.) then
            if (fac > 1.) split = .true.
         endif
         if ((abs(cellcentre%z-warpheight)/hr > 5.).and.(abs((cellcentre%z-warpheight)/cellsize) < 1.)) split = .true.
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

   case DEFAULT
      PRINT *, 'Invalid grid geometry option passed to amr_mod::decideSplit'
      PRINT *, 'grid%geometry ==', TRIM(grid%geometry)
      PRINT *, 'Exiting the program .... '
      STOP
   end select
   

   if (thisOctal%nDepth == 28) then
      split = .false.
      if (firstTime) then
         call writeWarning("AMR cell depth capped at 28")
         firstTime = .false.
      endif
   endif

  END FUNCTION decideSplit


  FUNCTION columnDensity(grid,startPoint,direction,sampleFreq) RESULT(rho)
                 
    ! integrates the density along a line through an amr grid.

    IMPLICIT NONE

    TYPE(gridtype), INTENT(IN)        :: grid
    TYPE(octalVector), INTENT(IN)     :: startPoint
    TYPE(octalVector), INTENT(IN)     :: direction
    REAL, INTENT(IN)  :: sampleFreq

    REAL(double)                              :: rho

    INTEGER, PARAMETER                :: maxSamples = 10000
    REAL, DIMENSION(maxSamples)       :: distances
    real(double), dimension(maxSamples) :: densities
    REAL, DIMENSION(maxSamples)       :: dummy
    REAL(double), DIMENSION(maxSamples)       :: ddummy
    real(double), DIMENSION(maxSamples,1) :: dummyPops
    TYPE(vector),DIMENSION(maxSamples):: dummyVel
    LOGICAL                           :: hitCore
    
    INTEGER :: nSamples, iSample

    INTEGER :: error

    error = 0 

    nSamples = 0
                
    CALL startReturnSamples(startPoint,direction,grid,sampleFreq,     &
                 nSamples,maxSamples,.false.,.false.,hitCore,.false.,1,error, &
                 distances,kappaAbs=ddummy,kappaSca=ddummy,velocity=dummyVel,&
                 velocityDeriv=dummy,chiLine=dummy,                &
                 levelPop=dummyPops,rho=densities)
 
    IF (nSamples <= 1 .OR. error /= 0)  THEN  
      rho = 0.0
    ELSE
      rho = 0.0
      DO iSample = 2, nSamples, 1
         rho = rho + ( (densities(iSample) + densities(iSample-1)) / 2.0 ) * &
                    ABS(distances(iSample-1) - distances(iSample) * 1.e10  )                    
      END DO
    END IF

  END FUNCTION columnDensity



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

    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    
    INTERFACE 
      TYPE(vector) FUNCTION velocityFunc(point,grid)
        USE vector_mod
        USE gridtype_mod
        TYPE(octalVector), INTENT(IN) :: point
        TYPE(gridtype), INTENT(IN)    :: grid
      END FUNCTION velocityFunc
    END INTERFACE

    if (threed) then

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
       
       thisOctal%cornerVelocity(1) = velocityFunc(octalVector(x1,y1,z1),grid)
       thisOctal%cornerVelocity(2) = velocityFunc(octalVector(x2,y1,z1),grid)
       thisOctal%cornerVelocity(3) = velocityFunc(octalVector(x3,y1,z1),grid)
       thisOctal%cornerVelocity(4) = velocityFunc(octalVector(x1,y2,z1),grid)
       thisOctal%cornerVelocity(5) = velocityFunc(octalVector(x2,y2,z1),grid)
       thisOctal%cornerVelocity(6) = velocityFunc(octalVector(x3,y2,z1),grid)
       thisOctal%cornerVelocity(7) = velocityFunc(octalVector(x1,y3,z1),grid)
       thisOctal%cornerVelocity(8) = velocityFunc(octalVector(x2,y3,z1),grid)
       thisOctal%cornerVelocity(9) = velocityFunc(octalVector(x3,y3,z1),grid)
       
       ! middle level
       
       thisOctal%cornerVelocity(10) = velocityFunc(octalVector(x1,y1,z2),grid)
       thisOctal%cornerVelocity(11) = velocityFunc(octalVector(x2,y1,z2),grid)
       thisOctal%cornerVelocity(12) = velocityFunc(octalVector(x3,y1,z2),grid)
       thisOctal%cornerVelocity(13) = velocityFunc(octalVector(x1,y2,z2),grid)
       thisOctal%cornerVelocity(14) = velocityFunc(octalVector(x2,y2,z2),grid)
       thisOctal%cornerVelocity(15) = velocityFunc(octalVector(x3,y2,z2),grid)
       thisOctal%cornerVelocity(16) = velocityFunc(octalVector(x1,y3,z2),grid)
       thisOctal%cornerVelocity(17) = velocityFunc(octalVector(x2,y3,z2),grid)
       thisOctal%cornerVelocity(18) = velocityFunc(octalVector(x3,y3,z2),grid)
       
       ! top level
       
       thisOctal%cornerVelocity(19) = velocityFunc(octalVector(x1,y1,z3),grid)
       thisOctal%cornerVelocity(20) = velocityFunc(octalVector(x2,y1,z3),grid)
       thisOctal%cornerVelocity(21) = velocityFunc(octalVector(x3,y1,z3),grid)
       thisOctal%cornerVelocity(22) = velocityFunc(octalVector(x1,y2,z3),grid)
       thisOctal%cornerVelocity(23) = velocityFunc(octalVector(x2,y2,z3),grid)
       thisOctal%cornerVelocity(24) = velocityFunc(octalVector(x3,y2,z3),grid)
       thisOctal%cornerVelocity(25) = velocityFunc(octalVector(x1,y3,z3),grid)
       thisOctal%cornerVelocity(26) = velocityFunc(octalVector(x2,y3,z3),grid)
       thisOctal%cornerVelocity(27) = velocityFunc(octalVector(x3,y3,z3),grid)
       
       
    else
       
       
    ! we first store the values we use to assemble the position vectors
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values
       
       thisOctal%cornerVelocity(1) = velocityFunc(octalVector(x1,0.d0,z1),grid)
       thisOctal%cornerVelocity(2) = velocityFunc(octalVector(x2,0.d0,z1),grid)
       thisOctal%cornerVelocity(3) = velocityFunc(octalVector(x3,0.d0,z1),grid)
       thisOctal%cornerVelocity(4) = velocityFunc(octalVector(x1,0.d0,z2),grid)
       thisOctal%cornerVelocity(5) = velocityFunc(octalVector(x2,0.d0,z2),grid)
       thisOctal%cornerVelocity(6) = velocityFunc(octalVector(x3,0.d0,z2),grid)
       thisOctal%cornerVelocity(7) = velocityFunc(octalVector(x1,0.d0,z3),grid)
       thisOctal%cornerVelocity(8) = velocityFunc(octalVector(x2,0.d0,z3),grid)
       thisOctal%cornerVelocity(9) = velocityFunc(octalVector(x3,0.d0,z3),grid)
    endif

  END SUBROUTINE fillVelocityCorners

  
  TYPE(vector) FUNCTION TTauriVelocity(point,grid)
    ! calculates the velocity vector at a given point for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    use input_variables, only: TTauriRinner, TTauriRouter, TTauriRstar, &
                               TTauriMstar, TTauriDiskHeight
                               
    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    TYPE(octalvector)      :: vP
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
    ELSE IF ((rM > TTauriRinner*0.99d0) .AND. (rM < TTauriRouter*1.01d0 )) THEN
  
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

  SUBROUTINE equationA7(m,func,deriv)
  
    USE rotation_variables
    use input_variables, only: TTauriMstar
    real(double), DIMENSION(:), INTENT(IN) :: m
    real(double), DIMENSION(:), INTENT(OUT) :: func
    real(double), DIMENSION(:,:), INTENT(OUT) :: deriv

    !print *, 'eta = ',eta
    !print *, 'm = ',m
    !print *, 'Bp = ',Bp
    !print *, 'Einitial = ',Einitial
    !print *, '2ndterm = ',-(Bp**2 * m**2)/(8.0*pi*eta)
    !print *, '3rdterm = ',bigG*TTauriMstar/r
    

    !print *,(Bp**2/(8.0_db*pi*eta)) * m**4                                    
    !print *, - (Einitial + bigG*TTauriMstar/r - l**2/(2.0_db*bigR**2) - (Bp**2/(8.0_db*pi*eta))) * m**2 
    !print *,        + (bigR**2 * omegaStar * omega) * m                              
    !print *,        - (0.5_db * bigR**2 * omega**2) + Einitial + (bigG*TTauriMstar/r)

    func = (Bp**2/(8.0_db*pi*eta)) * m**4                                    &
            - (Einitial + bigG*TTauriMstar/r + l**2/(2.0_db*bigR**2) + (Bp**2/(8.0_db*pi*eta))) * m**2 &
            + (omegaStar * l) * m                                            &
            - (0.5_db * bigR**2 * omegaStar**2) + Einitial + (bigG*TTauriMstar/r)

    deriv = spread(dim=1,ncopies=1,source=                                        &
              (4.0 * Bp**2/(8.0_db*pi*eta)) * m**3                                &
              + 2.0 * (-Einitial - bigG*TTauriMstar/r - l**2/(2.0_db*bigR**2) - (Bp**2/(8.0_db*pi*eta))) * m &
              + (omegaStar * l))

    
  END SUBROUTINE equationA7
  
  real(double) PURE FUNCTION funcA7(m)

    USE rotation_variables
    use input_variables, only: TTauriMstar
    real(double), INTENT(IN) :: m
!!!!!!!!!!!!!! INCORRECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    funcA7 = (Bp**2/8.0_db*pi*eta) * m**4                                      &
              - (Einitial + bigG*TTauriMstar/r - l**2/(2.0_db*bigR**2)) * m**2 &
              + (bigR**2 * omegaStar * omega) * m                              &
              - (0.5_db * bigR**2 * omega**2)

  END FUNCTION funcA7

  SUBROUTINE calcTTauriMassVelocity(thisOctal,subcell,grid) 
    ! calculates some of the variables at a given point for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    USE input_variables, ONLY : useHartmannTemp, maxHartTemp, TTauriRstar,&
                                TTauriRinner, TTauriRouter

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(INOUT) :: grid

    TYPE(octalVector) :: point
    REAL :: rho

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    TYPE(vector) :: vP

    REAL :: modVp
    REAL :: phi
    REAL :: r, rM, theta, y, ang, bigR, thetaStar
    REAL :: bigRstar, thetaStarHartmann, bigRstarHartmann, rMnorm
    REAL :: tmp, sum
    integer :: i

    starPosn = grid%starPos1

    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn) * 1.e10_oc
    thisOctal%inflow(subcell) = .false.
    
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
        theta = MIN(theta,real(pi/2.0)-0.0001)
        theta = MAX(theta,0.0001)
        
        rMnorm = MAX(rMnorm, 0.0001)
        rMnorm = MIN(rMnorm, 0.9999)
!        thisOctal%temperature(subcell) = hartmannTemp(rMnorm, theta, maxHartTemp)
        thisOctal%temperature(subcell) = hartmannTemp2(rMnorm, theta, maxHartTemp)
        thisOctal%inFlow(subcell) = .true.
     END IF
  ELSE
     thisOctal%inFlow(subcell) = .FALSE.
!     thisOctal%rho(subcell) = 1.e-25
     thisOctal%rho(subcell) = 1.e-19
     IF (useHartmannTemp) thisOctal%temperature(subcell) = 6500.0
  END IF

  thisOctal%velocity(subcell) = TTauriVelocity(point,grid)
  !thisOctal%velocity(subcell) = TTauriRotation(point,grid)    
     

  IF ((thisoctal%threed).and.(subcell == 8)) &
       CALL fillVelocityCorners(thisOctal,grid,TTauriVelocity, .true.)
  
  IF ((thisoctal%twod).and.(subcell == 4)) &
       CALL fillVelocityCorners(thisOctal,grid,TTauriVelocity, .false.)
  
!  if (subcell ==4) then
!     sum = 0.
!     do i = 1, 9
!        sum = sum + modulus(thisOctal%cornerVelocity(i))
!     enddo
!     if (sum < 1.e-10) then
!        if (thisOctal%inflow(1).or.thisOctal%inflow(2).or.thisOctal%inflow(3).or.thisOctal%inflow(4)) then
!           write(*,*) "velocity test",thisOctal%inflow(1:4)
!           write(*,*) thisOCtal%cornervelocity(1:9)
!        endif
!     endif
!  endif

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
    TYPE(octalVector)     :: rvec
    integer :: n
    real(double) :: r, theta, rM, rM_center, h, w, p, rho
    real(double) :: rM_fuzzy_in, rM_fuzzy_out  ! beginning of the fuzzy edges
    !
!    real(double), parameter :: scale = 7.0d0  ! a scale in exponential decay of density.
    real(double), parameter :: scale = 10.0d0  ! a scale in exponential decay of density.
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
             
             ! bias more toward the edges of the accreation stream
             h = 0.5d0*(TTauriRouter - TTauriRinner)*1.0d-10  ! [10^10cm] distance from the centere to edge
             rM_center = 0.5d0*(TTauriRouter + TTauriRinner)*1.0d-10 ! [10^10cm]   
             
             ! The fuzzy density starts from a 5-th of the thickness (2h) 
             ! below the surface.
!             w = 0.1d0*h;  ! Halpha081 model
!             w = 0.2d0*h;  ! a fifth for now  ! Halpha067 model
!             w = 0.25d0*h;  ! Halpha083 model
             w = 0.3d0*h;  ! Halpha080 model

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


  ! find the total mass in accretion flow
  ! Before the initial call the total mass must be set to zero
  ! The output should be in grams.
  RECURSIVE SUBROUTINE TTauri_accretion_mass(thisOctal, grid, total_mass)
    use input_variables, only: TTauriRinner, TTauriRouter

    IMPLICIT NONE
    TYPE(OCTAL), POINTER :: thisOctal    
    TYPE(gridtype), INTENT(INOUT)    :: grid
    real(double), intent(inout) :: total_mass
    !
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    TYPE(octalVector)     :: point
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
    use input_variables, only: TTauriRinner, TTauriRouter

    IMPLICIT NONE
    TYPE(OCTAL), POINTER :: thisOctal    
    TYPE(gridtype), INTENT(INOUT)    :: grid
    real(double), intent(in) :: scale
    !
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    TYPE(octalVector)     :: point
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
    type(octalVector) :: distortionVec(nVec)
    type(octalVector) :: thisVec
    type(vector) :: thisVel
    type(octalVector) :: thisVelOc
    logical :: setAllChanged
    
    integer, parameter :: nTimes = 1000
    integer :: i, j
    type(octalVector) :: starPos

    ! if we are not undoing a previous infall, we will assume that we need to
    ! run stateq on ALL subcells
    setAllChanged = .not. alreadyDoneInfall
    
    do i = 1, nVec
       distortionVec(i) = s2o(distortionVector(i))
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
       distortionVector(i) = o2s(distortionVec(i))
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
      type(octalVector), intent(in) :: distortionVec(nVec)
      logical, intent(in) :: undoPrevious
      logical, intent(in) :: setAllChanged
      integer, parameter :: nPhi = 360
      real(oct) :: phi
      integer :: containedParticles
      real :: addedDensity
      real :: deltaRho
      integer :: iSubcell, iChild
      type(octal), pointer :: child
      type(octalVector) :: trueVector
      type(octalVector) :: thisVec, thisVec2
 
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

    USE parameters_mod
    USE input_variables, ONLY: isoTherm, isoThermTemp

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell

    REAL :: rho


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
    
    INTEGER :: nLines, subPos
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
    
    INTEGER :: nLines, subPos
    INTEGER :: fieldline
    INTEGER :: lowerLine, upperLine
    INTEGER :: iSample
    INTEGER :: maxSamples
    REAL :: lowerTemp, upperTemp
    INTEGER :: iStat
    logical, save :: warned_already_01 = .false.
    logical, save :: warned_already_02 = .false.
    real:: T(5), R(5), dT, R_ip
    integer :: indx_T(5)
    

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
    REAL, INTENT(IN)                 :: phi    ! azimuth angle (radians)
    TYPE(gridtype), INTENT(IN)       :: grid
    TYPE(octalVector), INTENT(OUT)   :: point  ! coordinates returned
    TYPE(vector), INTENT(OUT)        :: azVec  ! azimuth vector
    TYPE(vector), INTENT(OUT)        :: polVec ! poloidal vector
    LOGICAL, INTENT(OUT)             :: ok     ! coordinates are valid
    
    TYPE(octalVector) :: starPosn
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
    point = octalVector(RAD, 0., radius*cos(theta))
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
    REAL, DIMENSION(samplesPerBin)     :: samples
    TYPE(vector)                       :: azVec  ! azimuth vector
    TYPE(vector)                       :: polVec ! poloidal vector
    
    INTEGER           :: iBin, iRadius, iSample, iFieldLine
    REAL              :: rand
    REAL              :: phi
    real(oct) :: radius
    LOGICAL           :: ok
    TYPE(octalVector) :: point 
    REAL              :: value
    real(double) :: valueDouble
    TYPE(vector)      :: velocityValue
    REAL              :: etaLine, chiLine
    REAL,DIMENSION(grid%maxLevels) :: valueArraySingle
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
              valueDouble = etaLine / (chiLine*curtainsFill*10.) 
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
                               TTauriMstar, dipoleOffset, contFluxFile

    implicit none
    
    type(GRIDTYPE), intent(inout)      :: grid                
    real, intent(out) :: theta1, theta2
   
    real(double) :: TaccretionDouble
    integer :: nNu 
    real :: nuArray(3000),fnu(3000)
    real :: tot
    integer :: i 
    real :: rStar, saccretion, taccretion

    rStar  = TTauriRstar / 1.e10
    
    grid%rCore = rStar
    grid%rStar1 = rStar
    grid%rStar2 = 0.
    grid%dipoleOffset = dipoleOffset
    grid%diskRadius = TTauriRInner * 1.e-10
    grid%diskNormal = VECTOR(0.,0.,1.)
    grid%diskNormal = rotateX(grid%diskNormal,grid%dipoleOffSet)
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

    sAccretion = (fourPi * TTauriRstar**2)*abs(cos(theta1)-cos(theta2))!1.e20
    Taccretion = TaccretionDouble**0.25

    write(*,*) contfluxfile
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
    write(*,*) (fourPi*TTauriRstar**2)*tot/lSol," solar luminosities"   

    write(*,*) (fourPi*TTauriRstar**2)*tot/ &
               (fourPi*TTauriRstar**2*stefanBoltz*4000.**4)

    ! add the accretion luminosity spectrum to the stellar spectrum,
    ! write it out and pass it to the stateq routine.
    
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
    TYPE(octalVector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.e-20
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

    thisOctal%diffusionApprox(subcell) = .false.
    if ((r > grid%rinner).and.(r < grid%rinner*2.)) then
       thisOctal%temperature(subcell) = 1000.
    endif
    if ((r > grid%rOuter*0.5).and.(r < grid%rOuter)) then
       thisOctal%temperature(subcell) = 100.
    endif
    if ((r > grid%rinner*2.).and.(r < 0.5*grid%rOuter)) then
       thisOctal%diffusionApprox(subcell) = .true.
    endif
    

  end subroutine calcTestDensity

  subroutine calcLexington(thisOctal,subcell,grid)

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(octalVector) :: rVec
    integer,save :: it
    real :: fac
    real,save :: radius(400),temp(400)
    integer :: i
    logical,save :: firsttime = .true.

    if (firsttime) then
       open(20,file="overview.txt", form="formatted",status="old")
       it = 1
77     continue
       read(20,*,end=88) radius(it),temp(it)
       radius(it) = (radius(it)+30.e17) / 1.e10
       temp(it) = 10.**temp(it) 
       it = it + 1
       goto 77
88     continue
       it = it - 1
       close(20)
       firsttime = .false.
    endif


    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = tiny(thisOctal%rho(subcell))
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
    thisOctal%ne(subcell) = thisOctal%nh(subcell)
    thisOctal%nhi(subcell) = 1.e-8
    thisOctal%nhii(subcell) = thisOctal%ne(subcell)
    thisOctal%inFlow(subcell) = .false.


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
       thisOctal%temperature(subcell) = 8000.
       if ((r > radius(1)).and.(r < radius(it))) then
          call locate(radius, it, r, i)
          fac = (r-radius(i))/(radius(i+1)-radius(i))
          thisOctal%temperature(subcell) = temp(i) + fac * (temp(i+1)-temp(i))
       endif

    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30

  end subroutine calcLexington
  
  
  subroutine initWindTestAMR(grid)

    use constants_mod
    use vector_mod
    use input_variables
    use parameters_mod

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
    TYPE(octalVector) :: rVec
    
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

  subroutine calcWRShellDensity(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r
    TYPE(octalVector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 1.e-3
    thisOctal%etaCont(subcell) = 1.e-30
    thisOctal%inFlow(subcell) = .true.

    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       thisOctal%rho(subcell) = density(rVec, grid)
       thisOctal%temperature(subcell) = 10.
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaCont(subcell) = 0.
    endif

    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine calcWRShellDensity

  subroutine benchmarkDisk(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, hr, rd, r1
    TYPE(octalVector) :: rVec

    real :: rInnerGap, rOuterGap
    logical :: gap

    gap = .false.

    rInnerGap = 2. * auToCm / 1.e10
    rOuterGap = 3. * auToCm / 1.e10
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.e-33
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .false.
    rd = rOuter / 2.
    r = sqrt(rVec%x**2 + rVec%y**2)
!    if (gap.and.((r < rInnerGap).or.(r > rOuterGap))) then
       if ((r > rInner).and.(r < rOuter)) then
          hr = height * (r/rd)**1.125
          thisOctal%rho(subcell) = rho * ((r / rd)**(-1.))*exp(-pi/4.*(rVec%z/hr)**2)
          thisOctal%rho(subcell) = max(thisOctal%rho(subcell), 1.d-33)
          thisOctal%temperature(subcell) = 100.
          thisOctal%inFlow(subcell) = .true.
          thisOctal%etaCont(subcell) = 0.
       endif
!    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine benchmarkDisk

  subroutine assign_melvin(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, hr, rd, r1
    TYPE(octalVector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = melvinDensity(rVec, grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .false.
    if (thisOctal%rho(subcell) > 1.e-20) then
       thisOctal%inFlow(subcell) = .true.
    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine assign_melvin

  subroutine assign_whitney(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, hr, rd, r1
    TYPE(octalVector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    thisOctal%rho(subcell) = whitneyDensity(rVec, grid)
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .true.
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine assign_whitney

  subroutine assign_toruslogo(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, hr, rd, r1
    TYPE(octalVector) :: rVec
    
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
    real :: r, h, rd, r1
    TYPE(octalVector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)
    thisOctal%inflow(subcell) = .true.
    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    rd = rOuter / 2.
    r = sqrt(rVec%x**2 + rVec%y**2)
    if ((r > rInner).and.(r < rOuter)) then
       thisOctal%rho(subcell) = density(rVec, grid)
       thisOctal%temperature(subcell) = 10.
       thisOctal%etaCont(subcell) = 0.
       thisOctal%inflow(subcell) = .true.


       h = height * (r / (100.d0*autocm/1.d10))**betaDisc
    
    endif


    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine shakaraDisk

  subroutine warpedDisk(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, h, rd, r1
    TYPE(octalVector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)
    thisOctal%inflow(subcell) = .true.
    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    rd = rOuter / 2.
    r = sqrt(rVec%x**2 + rVec%y**2)
    if ((r > rInner).and.(r < rOuter)) then
       thisOctal%rho(subcell) = density(rVec, grid)
       thisOctal%temperature(subcell) = 20.
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

    type(octalVector) :: rVec

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
    TYPE(octalVector) :: rVec
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

  
  SUBROUTINE addNewChildren(parent, grid, sphData, stellar_cluster, inherit, interp)
    ! adds all eight new children to an octal
    use input_variables, only : nDustType, mie, photoionization
    IMPLICIT NONE
    
    TYPE(octal), POINTER :: parent     ! pointer to the parent octal 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that
                                          !   calculate the variables stored in
                                          !   the tree.
    TYPE(cluster), optional, INTENT(IN)   :: stellar_cluster
    INTEGER       :: subcell           ! loop counter
    INTEGER, SAVE :: counter = 9       ! keeps track of the current subcell label
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: newChildIndex     ! the storage location for the new child
    INTEGER       :: error
    LOGICAL, OPTIONAL :: inherit       ! inherit densities, temps, etc from parent
    LOGICAL, OPTIONAL :: interp        ! interpolate densities, temps, etc from parent
    LOGICAL       :: inheritProps
    LOGICAL       :: interpolate

    ! For only cluster geometry ...
    TYPE(sph_data), optional, intent(in) :: sphData
    
    inheritProps = .false.
    if (present(inherit)) then
       inheritProps = inherit
    endif

    interpolate = .false.
    if (present(interp)) then
       interpolate = interp
    endif

    
    if (any(parent%hasChild(1:parent%maxChildren))) write(*,*) "parent already has a child"

    if (parent%threeD) then
       parent%nChildren = 8
       parent%maxChildren = 8
    else
       parent%nChildren = 4
       parent%maxChildren = 4
    endif

    if (parent%threed) then
       ALLOCATE(parent%child(1:8), STAT=error)
    else
       ALLOCATE(parent%child(1:4), STAT=error)
    endif
    IF ( error /= 0 ) THEN
       PRINT *, 'Panic: allocation failed.'
       STOP
    END IF

    ! update the parent octal
    
    parent%hasChild(1:parent%maxChildren) = .TRUE.
    parent%indexChild(1) = 1
    parent%indexChild(2) = 2
    parent%indexChild(3) = 3
    parent%indexChild(4) = 4
    if (parent%threeD) then
       parent%indexChild(5) = 5
       parent%indexChild(6) = 6
       parent%indexChild(7) = 7
       parent%indexChild(8) = 8
    endif



    do newChildIndex = 1, parent%maxChildren

       ! allocate any variables that need to be  
       if (.not.grid%oneKappa) then
          ! The kappa arrays should be allocated with grid%nopacity instead of grid%nlambda
          ! because for line calculation, there is only one kappa needed.
          ! (but grid%nlambda is not 1). If you allocate the arrays with grid%nlambda,
          ! it will be a huge waste of RAM. ---  (RK) 
          ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nopacity))
          ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nopacity))
!          ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nlambda))
!          ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nlambda))
          parent%child(newChildIndex)%kappaAbs = 1.e-30
          parent%child(newChildIndex)%kappaSca = 1.e-30
       endif
       NULLIFY(parent%child(newChildIndex)%child)

       if (.not.mie) then
          ALLOCATE(parent%child(newChildIndex)%N(8,grid%maxLevels))
       endif

       ! set up the new child's variables
       parent%child(newChildIndex)%threed = parent%threed
       parent%child(newChildIndex)%twoD = parent%twod
       parent%child(newChildIndex)%maxChildren = parent%maxChildren
       parent%child(newChildIndex)%twoD = parent%twod
       parent%child(newChildIndex)%inFlow = parent%inFlow
       parent%child(newChildIndex)%parent => parent
       parent%child(newChildIndex)%parentSubcell = newChildIndex
       parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
       parent%child(newChildIndex)%hasChild = .false.
       parent%child(newChildIndex)%nChildren = 0
       parent%child(newChildIndex)%indexChild = -999 ! values are undefined
       parent%child(newChildIndex)%nDepth = parent%nDepth + 1
       parent%child(newChildIndex)%centre = subcellCentre(parent,newChildIndex)
       parent%child(newChildIndex)%probDistLine = 0.0
       parent%child(newChildIndex)%probDistCont = 0.0
       parent%child(newChildIndex)%biasLine3D = 1.0 
       parent%child(newChildIndex)%biasCont3D = 1.0 
       parent%child(newChildIndex)%velocity = vector(1.e-30,1.e-30,1.e-30)
       parent%child(newChildIndex)%cornerVelocity = vector(1.e-30,1.e-30,1.e-30)
       parent%child(newChildIndex)%chiLine = 1.e-30
       parent%child(newChildIndex)%etaLine = 1.e-30
       parent%child(newChildIndex)%etaCont = 1.e-30
       parent%child(newChildIndex)%rho = 1.e-30
       parent%child(newChildIndex)%N = 1.e-30
       parent%child(newChildIndex)%dusttype = 1
       ALLOCATE(parent%child(newChildIndex)%dusttypefraction(8,  nDustType))
       parent%child(newChildIndex)%dustTypeFraction(1:8,1:nDustType) = 0.d0
       parent%child(newChildIndex)%dustTypeFraction(1:8,1) = 1.d0
       parent%child(newChildIndex)%gasOpacity = .false.
       parent%child(newChildIndex)%Ne = 1.e-30
       parent%child(newChildIndex)%temperature = 3.0
       parent%child(newChildIndex)%nTot = 1.e-30
       parent%child(newChildIndex)%changed = .false.
       parent%child(newChildIndex)%diffusionApprox = .false.
       parent%child(newChildindex)%nDiffusion  = 0.
       parent%child(newChildindex)%oldFrac = 1.d-5
       if (photoionization) then
          allocate(parent%child(newChildIndex)%ionFrac(1:parent%maxChildren, 1:grid%nIon))
          allocate(parent%child(newChildIndex)%photoionCoeff(1:parent%maxChildren, 1:grid%nIon))
       endif

       if (associated(parent%child(newChildindex)%ionFrac)) parent%child(newChildIndex)%ionFrac = 1.e-30

       if (associated(parent%child(newChildIndex)%photoIonCoeff)) then
          parent%child(newChildIndex)%photoIonCoeff = 0.d0
       endif

       if (present(sphData)) then
          ! updates the sph particle list.           
          call update_particle_list(parent, newChildIndex, newChildIndex, sphData)
       end if
       
       ! put some data in the eight subcells of the new child
       DO subcell = 1, parent%maxChildren
          CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid, sphData, stellar_cluster,  &
               inheritProps, interpolate)
          parent%child(newChildIndex)%label(subcell) = counter
          counter = counter + 1
       END DO
 
    enddo

    ! check for a new maximum depth 
    IF (parent%child(1)%nDepth > grid%maxDepth) THEN
       grid%maxDepth = parent%child(1)%nDepth
       grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
            2.0_oc**REAL(grid%maxDepth,KIND=oct)
       ! we store the value which is half the size of the 
       !   smallest subcell because this is more useful for later
       !   calculations.
    END IF
    
  END SUBROUTINE addNewChildren
  !
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
  
  subroutine copy_sph_index_to_root(grid, sphData)
    implicit none
    type(gridtype), intent(inout) :: grid    
    type(sph_data), intent(in) :: sphData
    !
    integer :: i, n

    ! extracting the number of gas particles in sph
    n = get_npart(sphData)

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

  recursive subroutine findTotalMass(thisOctal, totalMass, minRho, maxRho)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: totalMass
  real(double),optional :: minRho, maxRho
  real(double) :: r1,r2,dv
  type(OCTALVECTOR) :: rVec
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findtotalMass(child, totalMass, minRho, maxRho)
                exit
             end if
          end do
       else

          dv = cellVolume(thisOctal, subcell)
          totalMass = totalMass + (1.d30)*thisOctal%rho(subcell) * dv
                
          if (PRESENT(minRho)) minRho = min(dble(thisOctal%rho(subcell)), minRho)
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
    use input_variables
    real :: v, r, rhoOut
    type(VECTOR) :: rVec, rHat
    type(OCTALVECTOR) :: octVec
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
    thisOctal%velocity(subcell) = OCTALVECTOR(0., 0., 0.)
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
       thisOctal%etaLine(subcell) =  logInterp(etal, nr, rgrid, r)
       thisOctal%etaCont(subcell) =  1.e-30
       thisOctal%chiLine(subcell) =  logInterp(chil, nr, rgrid, r)
       thisOctal%kappaSca(subcell,1) = logInterp(etal, nr, rgrid, r)/1.e10
       thisOctal%kappaAbs(subcell,1) = 1.e-30
    endif

    CALL fillVelocityCorners(thisOctal,grid,ostarVelocity,thisOctal%threed)

  end subroutine spiralWindSubcell


  recursive subroutine  createPointList(thisOctal, nPoints, points, nColour, maxPoints, totalMass, minRho, maxRho)
    type(OCTAL), pointer :: thisOctal
    type(OCTAL), pointer :: child
    type(VECTOR) :: points(:), xAxis, yAxis, zAxis, a
    integer :: nColour(:)
    integer :: maxpoints, nPoints
    integer :: subcell, i, j, m
    real(double) :: totalMass, minRho, maxRho
    real :: s, r1, r2, r3
    integer :: iclow, ichi
    xAxis = VECTOR(1., 0., 0.)
    yAxis = VECTOR(0., 1., 0.)
    zAxis = VECTOR(0., 0., 1.)
    call pgqcir(iclow,ichi)
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call createPointList(child, nPoints, points, nColour, maxPoints, totalMass, minrho, maxRho)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell)) then             
             m = (maxPoints-1000) * dble((1.d30*thisOctal%rho(subcell)*thisOctal%subcellSize**3)/totalMass)
             do j = 1 , m
                nPoints = nPoints + 1             
                call random_number(r1)
                call random_number(r2)
                call random_number(r3)
                r1 = 2.*r1 - 1.
                r2 = 2.*r2 - 1.
                r3 = 2.*r3 - 1.
                s = thisOctal%subCellSize / 2.
                a = subcellCentre(thisOctal, subcell)
                points(nPoints) = a + (r1 * s)*xAxis + (r2 * s)*yAxis + (r3 * s)*zAxis
                nColour(nPoints) = max(iclow,iclow+int(real(ichi-iclow) * &
                     (log10(thisOctal%rho(subcell))-log10(minRho))/(log10(maxRho)-log10(minRho))))
             end do
          endif
       endif
    enddo
  end subroutine createPointList

  TYPE(vector) FUNCTION ostarVelocity(point,grid)
    real(double) :: r1
    type(OCTALVECTOR), intent(in) :: point
    type(OCTALVECTOR) rHat
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

  subroutine cameraPositions(point, nPoints, cameraPos, cameraDirection, nScene, nFrames)
    integer :: nPoints, nScene, nFrames
    type(VECTOR) point(nPoints), cameraPos(:), cameraDirection(:)
    integer :: i
    real :: t
    do i = 1, nScene
       t = real(i-1)/real(nScene-1)
       nFrames = nFrames + 1
       call bezier(point(1), point(2), point(3), point(4), t, cameraPos(nFrames))
       cameraDirection(nFrames) = (-1.)*cameraPos(nFrames)
       call normalize(cameraDirection(nFrames))
    enddo
  end subroutine cameraPositions

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

    dest%nDepth           = source%nDepth
    dest%nChildren        = source%nChildren
    dest%indexChild       = source%indexChild
    dest%threeD           = source%threeD 
    dest%cylindrical      = source%cylindrical
    dest%splitAzimuthally= source%splitAzimuthally
    dest%twoD             = source%twoD   
    dest%maxChildren      = source%maxChildren
    dest%hasChild         = source%hasChild
    dest%centre           = source%centre
    dest%rho              = source%rho
    dest%phi              = source%phi
    dest%dphi             = source%dphi
    dest%r                = source%r
    dest%velocity         = source%velocity
    dest%nDirectPhotons = source%nDirectPhotons
    dest%temperature      = source%temperature
    dest%oldTemperature   = source%oldTemperature
    dest%distanceGrid     = source%distanceGrid
    dest%nCrossings       = source%nCrossings
    dest%chiLine          = source%chiLine
    dest%etaLine          = source%etaLine
    dest%etaCont          = source%etaCont
    dest%biasLine3D       = source%biasLine3D
    dest%biasCont3D       = source%biasCont3D
    dest%probDistLine     = source%probDistLine
    dest%probDistCont     = source%probDistCont
    dest%Ne               = source%Ne
    dest%NH               = source%NH
    dest%NHI              = source%NHI
    dest%NHII             = source%NHII
    dest%NHeI             = source%NHeI
    dest%NHeII            = source%NHeII
    dest%Hheating         = source%Hheating
    dest%Heheating        = source%Heheating
    dest%oldFrac          = source%oldFrac
    dest%nTot             = source%nTot
    dest%inStar           = source%inStar
    dest%inFlow           = source%inFlow
    dest%label            = source%label
    dest%subcellSize      = source%subcellSize
    dest%changed          = source%changed
    dest%dustType         = source%dustType         
    dest%parentSubcell    = source%parentSubcell    
    dest%gasOpacity       = source%gasOpacity       
    dest%diffusionApprox  = source%diffusionApprox  
    dest%underSampled     = source%underSampled
    dest%nDiffusion       = source%nDiffusion
    dest%cornerVelocity   = source%cornerVelocity


    if (associated(source%photoIonCoeff)) then
       allocate(dest%photoIonCoeff(SIZE(source%photoIonCoeff,1),SIZE(source%photoIonCoeff,2)))
       dest%photoIonCoeff          = source%photoIonCoeff
    endif

    if (associated(source%ionFrac)) then
       allocate(dest%ionFrac(SIZE(source%ionFrac,1),SIZE(source%ionFrac,2)))
       dest%ionFrac          = source%ionFrac
    endif


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

      INTEGER :: error 
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
                                        sphData, stellar_cluster, inheritProps, interpProps)
    USE input_variables, ONLY : tauSmoothMax,tauSmoothMin
    IMPLICIT NONE

    TYPE(octal), POINTER             :: thisOctal
    TYPE(gridtype), INTENT(INOUT   ) :: grid 
    LOGICAL, INTENT(INOUT)               :: gridConverged
    TYPE(sph_data), optional, INTENT(IN) :: sphData   ! Matthew's SPH data.
    TYPE(cluster), optional, intent(in)  :: stellar_cluster
    LOGICAL, INTENT(IN) :: inheritProps
    LOGICAL, INTENT(IN), optional :: interpProps

    INTEGER              :: i, ilam
    TYPE(octalVector)    :: thisSubcellCentre
    REAL(oct) :: dSubcellCentre
    real(double) :: kappaAbs, kappaSca
    TYPE(octal), POINTER :: child
    TYPE(octal), POINTER :: neighbour
    TYPE(octalVector), ALLOCATABLE, DIMENSION(:) :: locator
    INTEGER              :: subcell
    INTEGER              :: thisSubcell
    REAL                 :: thisTau, thatTau, tauDiff
    INTEGER              :: nLocator
    LOGICAL              :: prob
    
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
              call smoothAMRgridTau(child,grid,gridConverged,  ilam, sphData, stellar_cluster, inheritProps, interpProps)
              ! force return until algorithm is fixed
              IF ( .NOT. gridConverged ) RETURN
              exit
           end if
        end do
      else
        thisSubcellCentre = subcellCentre(thisOctal,thisSubcell)

        ! don't split outer edge of disc

        if ((grid%geometry == "shakara").or.(grid%geometry == "warpeddisc")) then
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
                  do;enddo
                END IF
                IF ( thisOctal%hasChild(thisSubcell) ) THEN 
                  PRINT *, "thisOctal already has child"
                  do;enddo
                END IF
                ! Note that the version of addNewChild called here does not
                ! support sphData, etc.
!                call addNewChild(neighbour,subcell,grid)
!                call addNewChild(thisOctal,thisSubcell,grid)

              call addNewChild(neighbour, subcell, grid, adjustGridInfo=.TRUE.,  &
                               sphData=sphData, stellar_cluster=stellar_cluster, &
                               inherit=inheritProps, interp=interpProps)
!              call addNewChild(thisOctal, thissubcell, grid, adjustGridInfo=.TRUE.,  &
!                               sphData=sphData, stellar_cluster=stellar_cluster, &
!                               inherit=inheritProps, interp=interpProps)
!              write(*,*) "added new child",thisOctal%hasChild,thisOctal%nchildren


!                call addNewChildren(thisOctal, grid, sphData, stellar_cluster, inheritProps, interpProps)
                ! If we are splitting two subcells in the same octal, we don't
                ! need to addNewChildren to the neighbour. If we switch to
                ! addNewChild this test becomes redundant.
!                if (.not. associated(thisOctal, neighbour)) then
!                  call addNewChildren(neighbour, grid, sphData, stellar_cluster, inheritProps, interpProps)
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
    type(OCTALVECTOR) :: rVec

    open(30, file="polardump.dat", status="unknown", form="formatted")
    do j = 1, nt
       do i = 1, nr
          r = grid%rInner + (grid%rOuter-grid%rInner)*(real(i-1)/real(nr-1))**3
          theta = pi*real(j-1)/real(nt-1)
          rVec = OCTALVECTOR(dble(r*sin(theta)),0.d0,dble(r*cos(theta)))
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


  end subroutine polarDump



  recursive subroutine verticalDump(thisOctal, xPos)
    real :: xPos
    type(OCTALVECTOR) :: rVec
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
             write(30,*) rVec%z, thisOctal%temperature(subcell),thisOctal%rho(subcell)
          endif
       endif
          
    enddo

  end subroutine verticalDump


  recursive subroutine unrefineThinCells(thisOctal, grid, ilambda, converged)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: ilambda
    real(double) :: kappaAbs, kappaSca, tau
    integer :: subcell, i, j
    logical :: unrefine, ok, converged
    integer :: nTau
    integer :: nVals, nc

    unrefine = .true.
    nc = 0
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call unrefineThinCells(child, grid, ilambda, converged)
                exit
             end if
          end do
       else
          if (.not.ASSOCIATED(thisOctal%dustTypeFraction)) then
             write(*,*) "unalloc dusttypefraction!!"
          endif
          call returnKappa(grid, thisOctal, subcell, ilambda, kappaAbs=kappaAbs,kappaSca=kappaSca)
          tau = thisOctal%subcellSize*(kappaAbs+kappaSca)
          nc = nc + 1
          if (tau > 1.e-10) then
             unrefine = .false.
             exit
          endif
       endif
    enddo
!    write(*,*) thisOctal%nChildren,unrefine,tau, nc,thisOctal%maxChildren

    if ((thisOctal%nChildren == 0).and.unrefine) then
       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
            grid = grid, adjustGridInfo = .true.)
       converged = .false.
    endif

  end subroutine unrefineThinCells

  recursive subroutine unrefineThickCells(thisOctal, grid, ilambda, converged)
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: ilambda
    real(double) :: kappaAbs, kappaSca, tau
    integer :: subcell, i, j
    logical :: unrefine, ok, converged
    integer :: nTau
    integer :: nVals

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
    
    TYPE(octal), TARGET, INTENT(INOUT) :: parent ! the parent octal 
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
    INTEGER :: nVals
    REAL(double) :: nValsREAL

    IF ( childOctal%nDepth == 1 ) THEN
      ! we're at the root of the tree
      PRINT *, "Warning: Attempted to update parent, but already at tree root."
      PRINT *, "         Ignoring this and continuing..."
      RETURN
    END IF
    
    parentOctal => childOctal%parent
    parentSubcell = childOctal%parentSubcell
    
    if (childOctal%twoD) then
       nVals = 4
    else
       nVals = 8
    endif

    nValsREAL = REAL( nVals, KIND=double )
    
    parentOctal%rho(parentSubcell) =                    &
    SUM(childOctal%rho(1:nVals)) / nValsREAL
    
    parentOctal%temperature(parentSubcell) =            &
    SUM(childOctal%temperature(1:nVals)) / nValsREAL
    
    parentOctal%dustTypeFraction(parentSubcell,:) =     &
    SUM(childOctal%dustTypeFraction(1:nVals,:)) / nValsREAL
    
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
      
        TYPE(OCTAL), TARGET, INTENT(IN) :: thisOctal 
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
    type(OCTALVECTOR) :: rVec
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
  use input_variables, only : rinner
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: kap
  real(double) :: tau, kappaAbs, kappaSca
  type(octalvector) :: rVec
  logical :: ross
  integer :: subcell, i, ilam
  real :: r
  
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
             thisOctal%biasCont3D(subcell) = MAX(exp(-tau),1.d-4)
          else
             call returnKappa(grid, thisOctal, subcell, ilam, rosselandKappa = kappaAbs)
             tau = thisOctal%subcellSize*(kappaAbs + kappaSca)*thisOctal%rho(subcell)*1.d10
             thisOctal%biasCont3D(subcell) = MAX(exp(-tau),1.d-4)
          endif
       endif

    enddo
  end subroutine set_bias_shakara

  recursive subroutine set_bias_rosseland(thisOctal, grid)
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: kap
  real(double) :: tau, kappaAbs, kappaSca
  type(octalvector) :: rVec
  logical :: ross
  integer :: subcell, i, ilam
  real :: r
  
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
  type(gridtype) :: grid
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real :: kap
  real(double) :: tau, kappaAbs, kappaSca
  type(octalvector) :: rVec
  logical :: ross
  integer :: subcell, i, ilam
  real :: r
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
  real :: kap
  real(double) :: tau, kappaAbs, kappaSca
  type(octalvector) :: rVec
  logical :: ross
  integer :: subcell, i, ilam
  real(double) :: r
  
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
             thisOctal%biasCont3D(subcell) = MAX(exp(-tau),1.d-15) * (1.d10*r/erOuter)
!             thisOctal%biasCont3D(subcell) = 1.d0/thisOctal%etacont(subcell)
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
  type(octalvector)  :: rvec, rhat
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
                rhat = OCTALVECTOR(0.0d0, 0.0d0, 1.0d0)
             end if

             if (thisOctal%threed) then
                dV = d*d*d
             else
                dv = 2.0_db*pi*d*d*rVec%x
             endif

             nu0  = cSpeed_dbl / dble(lambda0*angstromtocm)
              ! in a radial direction
             tauSob = thisOctal%chiline(subcell)  / nu0
             tauSob = tauSob / amrGridDirectionalDeriv(grid, rvec, outVec, &
                  startOctal=thisOctal)
           
             if (tauSob < 0.01) then
                escProb = 1.0d0-tauSob*0.5d0*(1.0d0 -   &
                     tauSob/3.0d0*(1. - tauSob*0.25d0*(1.0d0 - 0.20d0*tauSob)))
             else if (tauSob < 15.) then
                escProb = (1.0d0-exp(-tauSob))/tauSob
             else
                escProb = 1.d0/tauSob
             end if
             escProb = min(1.d-2,max(escProb, 1.d-5))


             dtau_cont = d*(thisOctal%kappaAbs(subcell,1) + thisOctal%kappaSca(subcell,1))
             dtau_line = d*(thisOctal%chiline(subcell))  / nu0
             !
             rho = thisOctal%rho(subcell)
!             thisOctal%biasCont3D(subcell) = MAX(dV, 1.d-7) ! Limits the minimum value
!             thisOctal%biasLine3D(subcell) = dV*thisOctal%biasCont3D(subcell)
!             thisOctal%biasCont3D(subcell) = MAX(EXP(-dtau_cont), 1.d-7) ! Limits the minimum value
             thisOctal%biasCont3D(subcell) = 1.0d0
!             thisOctal%biasLine3D(subcell) = thisOctal%rho(subcell)*thisOctal%biasCont3D(subcell)
             thisOctal%biasLine3D(subcell) = escProb
!             thisOctal%biasLine3D(subcell) = 1.0d0/(rho*rho)
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
    type(octalvector)  :: rvec
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
  type(octalvector)  :: rvec, rhat
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
                rhat = OCTALVECTOR(0.0d0, 0.0d0, 1.0d0)
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
             tauSob = tauSob / amrGridDirectionalDeriv(grid, rvec, o2s(rhat), &
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
  use input_variables, only: rInner
  type(gridtype) :: grid
  type(OCTALVECTOR) :: rVec
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(double) :: r
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
  type(OCTALVECTOR) :: rVec
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

  subroutine dumpSmoothedSurfaceDensity(filename, grid)
    type(GRIDTYPE) :: grid
    real :: mass, area, r1, r2, r
    integer :: i
    integer, parameter :: npts = 100
    real :: rAxis(npts),sigma
    character(len=*) :: filename


    write(*,*) "Calculating smoothed surface density profile..."
    open(20,file=filename,status="unknown",form="formatted")

    do i = 1, npts
       rAxis(i) = log10(grid%rInner)+real(i-1)/real(npts-1) * (log10(grid%rOuter)-log10(grid%rInner))
       rAxis(i) = 10.**rAxis(i)
    enddo
    do i = 2, npts
       r1 = rAxis(i-1)
       r2 = rAxis(i)
       r = (r1+r2)/2.
       mass = 0.
       area = pi*(rAxis(i)**2-rAxis(i-1)**2)
       call massInAnn(grid%octreeRoot, r1, r2, mass)
       if (area > 0.) then
          sigma = (mass/area)*1.e10
       else
          sigma = 0.
       endif
       write(20,*) r,sigma
    enddo
    close(20)
  end subroutine dumpSmoothedSurfaceDensity


  SUBROUTINE calcWindTestValues(thisOctal,subcell,grid) 
    ! calculates some of the variables for a spherical wind flow

    USE parameters_mod
    USE constants_mod

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(INOUT) :: grid

    TYPE(octalVector) :: point

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    TYPE(octalVector) :: pointVecNorm
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

    USE parameters_mod

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    
    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    TYPE(octalVector) :: pointVecNorm
    TYPE(Vector)      :: pointVecNormSingle

    REAL :: r, rStar
    REAL :: velocity, rho
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
    
    LOGICAL               :: intersectionFound
    LOGICAL               :: octalAdded
    TYPE(octalVector)     :: centrePoint
    TYPE(octalVector)     :: lineOrigin
    TYPE(octalVector)     :: starPos
    TYPE(octalVector)     :: nearestPlanePoint
    TYPE(octal), POINTER  :: child
    INTEGER               :: i, iSubcell, subIndex
    TYPE(octalVector)     :: planeVec
    TYPE(octalListElement), POINTER :: newElement
    
    octalAdded = .FALSE.
    starPos = grid%starPos1
    lineOrigin = starPos
    planeVec = octalVector(1.0_oc/SQRT(2.0_oc),1.0_oc/SQRT(2.0_oc),0.0_oc)

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
            do; enddo
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
    
    LOGICAL               :: intersectionFound
    LOGICAL               :: octalAdded
    TYPE(octalVector)     :: centrePoint
    TYPE(octalVector)     :: starPos
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
               do;enddo
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
       rosselandKappa, kappap, atthistemperature)
    use input_variables, only: includeGasOpacity, nDustType, photoionization
    implicit none
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    integer, optional :: ilambda
    real, optional :: lambda
    real(double), optional :: kappaSca, kappaAbs, kappaAbsArray(:), kappaScaArray(:)
    real(double), optional :: rosselandKappa
    real, optional :: kappap
    real, optional :: atthistemperature
    real :: temperature
    real :: kappaAbsGas, kappaScaGas
    real :: frac
    real :: tlambda
    real, parameter :: sublimationTemp = 1500., subRange = 100.
    real :: tArray(1000)
    real(double) :: freq, dfreq, bnutot, norm
    integer :: i,j,m
    real :: fac
    real :: e, h0, he0
    real(double) :: kappaH, kappaHe

    temperature = thisOctal%temperature(subcell)
!    if (temperature < sublimationTemp) frac = 1.
!    if (temperature > (sublimationTemp+subRange)) frac = 0.
!    
!    if ((temperature > sublimationTemp).and.(temperature < (sublimationTemp+subRange))) then
!       frac = 1.-(temperature-sublimationTemp)/subRange
!    endif

!    frac = max(1.e-20,frac)

    frac = 1.


    if (PRESENT(kappaAbsArray)) then

       kappaAbsArray(1:grid%nLambda) = 0.
       do i = 1, nDustType
          kappaAbsArray(1:grid%nLambda) = kappaAbsArray(1:grid%nLambda) + & 
               thisOctal%dustTypeFraction(subcell, i) * grid%oneKappaAbs(i,1:grid%nLambda)*thisOctal%rho(subcell) * frac
       enddo
       if (includeGasOpacity) then
          call returnGasKappaValue(temperature, thisOctal%rho(subcell),  kappaAbsArray=tarray)
          kappaAbsArray(1:grid%nLambda) = kappaAbsArray(1:grid%nLambda) + tarray(1:grid%nLambda)*thisOctal%rho(subcell)
       endif
!       write(*,*) nDustType,thisOctal%dusttypeFraction(subcell,1), grid%oneKappaAbs(1,1:grid%nLambda)

    endif
       
    if (PRESENT(kappaScaArray)) then

       kappaScaArray(1:grid%nLambda) = 0.
       do i = 1, nDustType
          kappaScaArray(1:grid%nLambda) = kappaScaArray(1:grid%nLambda) + & 
               thisOctal%dustTypeFraction(subcell, i) * grid%oneKappaSca(i,1:grid%nLambda)*thisOctal%rho(subcell) * frac
       enddo
       if (includeGasOpacity) then
          call returnGasKappaValue(temperature, thisOctal%rho(subcell),  kappaScaArray=tarray)
          kappaScaArray(1:grid%nLambda) = kappaScaArray(1:grid%nLambda) + tarray(1:grid%nLambda)*thisOctal%rho(subcell)
       endif
    endif


    if (PRESENT(kappaSca)) then
       kappaSca = 0.
       if (.not.PRESENT(lambda)) then
          tlambda = grid%lamArray(iLambda)
       else
          tlambda = lambda
       endif
       IF (.NOT.PRESENT(lambda)) THEN
          kappaSca = 0
          do i = 1, nDustType
             kappaSca = kappaSca + thisOctal%dustTypeFraction(subcell, i) * grid%oneKappaSca(i,iLambda)*thisOctal%rho(subcell) 
          enddo
       else
          kappaSca = 0
          do i = 1, nDustType
             kappaSca = kappaSca + thisOctal%dustTypeFraction(subcell, i) * &
                  logint(dble(lambda), dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), &
                  grid%oneKappaSca(i,iLambda)*thisOctal%rho(subcell), &
                  grid%oneKappaSca(i,iLambda+1)*thisOctal%rho(subcell))
          enddo
       endif
       kappaSca = kappaSca * frac
    endif

    if (PRESENT(kappaAbs)) then
       kappaAbs = 0.
       if (.not.PRESENT(lambda)) then
          tlambda = grid%lamArray(iLambda)
       else
          tlambda = lambda
       endif
       IF (.NOT.PRESENT(lambda)) THEN
            kappaAbs = 0
          do i = 1, nDustType
             kappaAbs = kappaAbs + thisOctal%dustTypeFraction(subcell, i) * grid%oneKappaAbs(i,iLambda)*thisOctal%rho(subcell)
          enddo
       else
          kappaAbs = 0
          do i = 1, nDustType
             kappaAbs = kappaAbs + thisOctal%dustTypeFraction(subcell, i) * &
                  logint(dble(lambda), dble(grid%lamArray(ilambda)), dble(grid%lamArray(ilambda+1)), &
                  grid%oneKappaAbs(i,iLambda)*thisOctal%rho(subcell), &
                  grid%oneKappaAbs(i,iLambda+1)*thisOctal%rho(subcell))
          enddo
       endif
!       write(*,*) nDustType,thisOctal%dusttypeFraction(subcell,1), grid%oneKappaAbs(1,1:grid%nLambda)
      kappaAbs = kappaAbs * frac
   endif


   if (includeGasOpacity) then
      if (PRESENT(kappaAbs)) then
         if (.not.PRESENT(lambda)) then
            tlambda = grid%lamArray(iLambda)
         else
            tlambda = lambda
         endif
         call returnGasKappaValue(temperature, thisOctal%rho(subcell), tlambda, kappaAbs=kappaAbsGas)
         kappaAbs = kappaAbs + kappaAbsGas*thisOctal%rho(subcell)
      endif
      
      if (PRESENT(kappaSca)) then
         call returnGasKappaValue(temperature, thisOctal%rho(subcell), tlambda, kappaSca=kappaScaGas)
         kappaSca = kappaSca + kappaScaGas*thisOctal%rho(subcell)
      endif
   endif

   if (PRESENT(rosselandKappa)) then
      if (PRESENT(atthistemperature)) then
         temperature = atthistemperature
      endif
      call locate(grid%tempRossArray, grid%nTempRossArray, temperature, m)
      fac = (temperature - grid%tempRossArray(m))/(grid%tempRossArray(m+1)-grid%tempRossArray(m))
      rosselandKappa = 0.
      
      do i = 1, nDustType
         rosselandKappa = rosselandKappa + (grid%kappaRossArray(i, m) + &
              fac*(grid%kappaRossArray(i,m+1)-grid%kappaRossArray(i,m))) * &
              max(1.d-30,thisOctal%dustTypeFraction(subcell, i))
      enddo

!      rosselandKappa = 0.
!      Bnutot = 0.
!      if (thisOctal%inFlow(subcell)) then
!         do i =  grid%nLambda,2,-1
!            freq = cSpeed / (grid%lamArray(i)*1.e-8)
!            dfreq = cSpeed / (grid%lamArray(i)*1.e-8) - cSpeed / (grid%lamArray(i-1)*1.e-8)
!            do j = 1, nDustType
!               rosselandKappa = rosselandKappa + bnu(freq, dble(temperature)) * dFreq / &
!                 ((grid%oneKappaabs(j,i)+grid%oneKappaSca(j,i))*max(1.d-30,thisOctal%dustTypeFraction(subcell, j)))
!            enddo
!            bnutot = bnutot + bnu(freq, dble(temperature)) * dfreq
!         enddo
!         if (rosselandkappa /= 0.) then
!            rosselandKappa = (bnutot / rosselandKappa)/1.d10
!         endif
!      endif
   endif
   

   if (PRESENT(kappap)) then
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
                 dble(bnu(dble(freq),dble(temperature)))  * dfreq
         enddo
         norm = norm + dble(bnu(dble(freq),dble(temperature)))  * dfreq
      enddo
      kappaP = kappaP / norm /1.d10
   endif

   if (photoionization) then
      if (PRESENT(kappaAbs)) then
         kappaAbs = 0. !!!!!!!!!!!!!!!!!!!!!!! need to change when dust included!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
         kappaAbs = kappaH + kappaHe
      endif
      if (PRESENT(kappaSca)) then
         kappaSca = kappaSca + thisOctal%ne(subcell) * sigmaE * 1.e10
      endif
   endif
   
  end subroutine returnKappa
   


   subroutine averageofNearbyCells(grid, thisOctal, thisSubcell, meantemp, meanrho)
     type(GRIDTYPE) :: grid
     type(OCTAL),pointer :: thisOctal
     type(OCTAL),pointer :: nearbyOctal
     real(double) :: rho, meanRho
     integer :: thisSubcell
     integer :: found
     integer :: subcell
     real(double) :: r
     type(OCTALVECTOR) :: rVec, aVec
     integer :: i, nMonte
     real :: temp, meanTemp
     integer :: nTemp
     logical :: differentOctal
     type(OCTALVECTOR), allocatable :: locator(:)
     real(double) :: dSubcellcentre
     type(OCTALVECTOR) :: thisSubcellCentre
     integer :: nlocator, j
     nMonte = 10000
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
        if (rho /= 0.d0) then
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


  recursive subroutine estimateRhoOfEmpty(grid, thisOctal, sphData)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    type(sph_data), intent(in) :: sphData
    real :: temp
    real(double) :: rho, rho_tmp
    integer :: subcell, i, np
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call estimateRhoOfEmpty(grid, child, sphData)
                exit
             end if
          end do
       else
          call find_n_particle_in_subcell(np, rho_tmp, sphData, &
               thisOctal, subcell)
          if (np < 1) then
             call averageofNearbyCells(grid, thisOctal, subcell, temp, rho)
             thisOctal%rho(subcell) = rho
          endif
       endif
    enddo
  end subroutine estimateRhoOfEmpty




  recursive subroutine updateTemps(grid, thisOctal)
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call updateTemps(grid, child)
                exit
             end if
          end do
       else
          if (thisOctal%undersampled(subcell)) then
             thisOctal%temperature(subcell) = thisOctal%chiLine(subcell)
          endif
       endif
    enddo
  end subroutine updateTemps

  subroutine amrGridDirectionalDeriv_sub(grid,position,direction,startOctal,&
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
    TYPE(octalVector), INTENT(IN)  :: position
    TYPE(vector), INTENT(IN)       :: direction
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    TYPE(octal),  POINTER         :: thisOctal, currentOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell

    TYPE(octalVector)              :: octalDirection
    TYPE(octal), POINTER           :: firstOctal
    real(oct)           :: dr, dx, dphi
    real(oct)           :: r
    real(oct)           :: phi1, phi2
    TYPE(octalVector)              :: position1
    TYPE(octalVector)              :: position2
    INTEGER                        :: subcell

    
    octalDirection = direction

    call amrGridValues(grid%octreeRoot,position,foundOctal=thisOctal,&
                           foundSubcell=subcell)
    
!    currentOctal => grid%octreeRoot
!    thisOctal => grid%octreeRoot
!    CALL findSubcellTD(position, currentOctal, thisOctal, subcell)
    ! dr is a small increment of distance
    dr = thisOctal%subcellSize 

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

    write(*,*) "phi1 phi2",phi1,phi2
    write(*,*) "position1",position1
    write(*,*) "position2",position2
    write(*,*) "dx",dx

    ! the line of sight velocity gradient

    IF (ABS(dx) >= 1.e-20) THEN
       amrGridDirectionalDeriv = (abs(dphi / dx))
    ELSE
       amrGridDirectionalDeriv = 1.e-20
    ENDIF
    IF (.NOT. ABS(amrGridDirectionalDeriv) >= 1.e-10) &
       amrGridDirectionalDeriv = 1.e-10

  END subroutine amrGridDirectionalDeriv_sub



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
    type(OCTALVECTOR) :: octVec

    write(*,*) "Setting disc photosphere bias..."
    nx = 0
    call getxValuesAMR(grid%octreeRoot, nx, xAxis)
    call stripSimilarValues(xAxis,nx,real(1.d-5*grid%halfSmallestSubcell))

    tau(1) = 0.
    do i = 2, nx
       octVec = OCTALVECTOR(xAxis(i), 0.d0, 0.d0)
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
       octVec = OCTALVECTOR(xAxis(i), 0.d0, 0.d0)
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
       octVec = OCTALVECTOR(xAxis(i), 0.d0, 0.d0)
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
    type(octalvector) :: rVec
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
    type(OCTALVECTOR) :: currentPos, temp
    real :: halfSmallestSubcell
    real(double) :: kappaSca, kappaAbs

    nz = 0
    tau = 0.
    halfSmallestSubcell = grid%halfSmallestSubcell

    ! bottom up

    currentPos = OCTALVECTOR(xpos, yPos, -1.*grid%octreeRoot%subcellsize)
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
             
       currentPos = OCTALVECTOR(xpos, yPos, zAxis(nz)+0.5*thisOctal%subcellsize+halfSmallestSubcell)
    end do
    
    tau = 0.
    nz = 0
    currentPos = OCTALVECTOR(xpos, yPos, grid%octreeRoot%subcellsize)
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
       currentPos = OCTALVECTOR(xpos, yPos, zAxis(nz)-0.5*thisOctal%subcellsize-halfSmallestSubcell)
    end do
    
  end subroutine setVerticalBias

  subroutine interpFromParent(centre, cellSize, grid, temperature, density, dusttypeFraction)
    type(GRIDTYPE) :: grid
    real(double) :: cellSize
    real :: temperature
    real(double) :: density
    real(double) :: dusttypeFraction(:)
    real(double), allocatable :: tdusttype(:,:)
    real(double), allocatable :: rho(:)
    real, allocatable :: temp(:)
    type(OCTALVECTOR) :: centre, octVec
    real(double) :: r, t1, t2
    integer :: j

    r = cellSize/2.d0 + 0.1d0 * grid%halfSmallestSubcell
    if (grid%octreeRoot%twod) then
       j = 4
       allocate(tDustType(1:j,1:SIZE(dusttypeFraction)))
       allocate(rho(1:j), temp(1:j))

       octVec = centre + r * OCTALVECTOR(1.d0, 0.d0, 0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(1), rho=rho(1), dusttypeFraction=tdusttype(1,:))
       
       octVec = centre + r * OCTALVECTOR(-1.d0,0.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(2), rho=rho(2), dusttypeFraction=tdusttype(2,:))
       
       octVec = centre + r * OCTALVECTOR(0.d0,0.d0,-1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(3), rho=rho(3), dusttypeFraction=tdusttype(3,:))
       
       octVec = centre + r * OCTALVECTOR(0.d0,0.d0,+1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(4), rho=rho(4), dusttypeFraction=tdusttype(4,:))
    else
       j = 6
       allocate(tDustType(1:j,1:SIZE(dusttypeFraction)))
       allocate(rho(1:j), temp(1:j))
       octVec = centre + r * OCTALVECTOR(1.d0, 0.d0, 0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(1), rho=rho(1), dusttypeFraction=tdusttype(1,:))
       
       octVec = centre + r * OCTALVECTOR(-1.d0,0.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(2), rho=rho(2), dusttypeFraction=tdusttype(2,:))
       
       octVec = centre + r * OCTALVECTOR(0.d0,0.d0,-1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(3), rho=rho(3), dusttypeFraction=tdusttype(3,:))
       
       octVec = centre + r * OCTALVECTOR(0.d0,0.d0,+1.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(4), rho=rho(4), dusttypeFraction=tdusttype(4,:))

       octVec = centre + r * OCTALVECTOR(0.d0,+1.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(5), rho=rho(5), dusttypeFraction=tdusttype(5,:))

       octVec = centre + r * OCTALVECTOR(0.d0,-1.d0,0.d0)
       call amrGridValues(grid%octreeRoot, octVec, temperature=temp(6), rho=rho(6), dusttypeFraction=tdusttype(6,:))

    endif
    tdusttype=log10(max(1.d-30,tdusttype))
  
    rho = log10(rho)
    temp = log10(temp)

    density = sum(rho) / dble(j)
    density = 10.d0**density

    temperature = sum(temp) / real(j)
    temperature = 10.d0**temperature

    dusttypeFraction=SUM(tdusttype(1:j,:)) / dble(j)
    dusttypeFraction = 10.d0**dusttypeFraction

  end subroutine interpFromParent


  RECURSIVE SUBROUTINE smoothAMRgridIonization(thisOctal,grid,gridConverged, ilam, &
                                        sphData, stellar_cluster, inheritProps, interpProps)
    USE input_variables, ONLY : tauSmoothMax,tauSmoothMin
    IMPLICIT NONE

    TYPE(octal), POINTER             :: thisOctal
    TYPE(gridtype), INTENT(INOUT   ) :: grid 
    LOGICAL, INTENT(INOUT)               :: gridConverged
    TYPE(sph_data), optional, INTENT(IN) :: sphData   ! Matthew's SPH data.
    TYPE(cluster), optional, intent(in)  :: stellar_cluster
    LOGICAL, INTENT(IN) :: inheritProps
    LOGICAL, INTENT(IN), optional :: interpProps

    INTEGER              :: i, ilam
    TYPE(octalVector)    :: thisSubcellCentre
    REAL(oct) :: dSubcellCentre
    real(double) :: kappaAbs, kappaSca
    TYPE(octal), POINTER :: child
    TYPE(octal), POINTER :: neighbour
    TYPE(octalVector), ALLOCATABLE, DIMENSION(:) :: locator
    INTEGER              :: subcell
    INTEGER              :: thisSubcell
    REAL                 :: thisTau, thatTau, tauDiff
    INTEGER              :: nLocator
    
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
              call smoothAMRgridIonization(child,grid,gridConverged,  ilam, sphData, stellar_cluster, inheritProps, interpProps)
              exit
           end if
        end do
      else
        thisSubcellCentre = subcellCentre(thisOctal,thisSubcell)

        ! don't split outer edge of disc

        if (grid%geometry == "shakara") then
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


        thisTau = thisOctal%ionFrac(thisSubcell,1)

        DO i = 1, nLocator, 1
          IF ( inOctal(grid%octreeRoot,locator(i)) ) THEN
            neighbour => thisOctal
            CALL findSubcellLocal(locator(i),neighbour,subcell)

            ! The neighbouring subcell must be larger than the current subcell
            ! otherwise our locators won't cover all neighbouring subcells
            ! (and we'll hit cell boundaries, which is not good).
            IF ( neighbour%subcellSize >= thisOctal%subcellSize) THEN

               thatTau = neighbour%ionFrac(subcell,1)

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

               if ((abs(thisTau-thatTau) > 0.5).and. &
                    (thisOctal%photoIonCoeff(subcell,1) /= 0.d0)) then


!write (*,*) thisSubcellCentre%x, thisSubcellCentre%y, thisSubcellCentre%z, thisTau, thatTau, i
                ! Because addNewChild is broken, we must use addNewChildren below,
                ! which makes these tests on the current and neighbouring subcells
                ! insufficient. We perform them anyway.
                IF ( neighbour%hasChild(subcell) ) THEN 
                  PRINT *, "neighbour already has child"
                  do;enddo
                END IF
                IF ( thisOctal%hasChild(thisSubcell) ) THEN 
                  PRINT *, "thisOctal already has child"
                  do;enddo
                END IF
                ! Note that the version of addNewChild called here does not
                ! support sphData, etc.
!                call addNewChild(neighbour,subcell,grid)
!                call addNewChild(thisOctal,thisSubcell,grid)

                call addNewChildren(thisOctal, grid, sphData, stellar_cluster, inheritProps, interpProps)
                ! If we are splitting two subcells in the same octal, we don't
                ! need to addNewChildren to the neighbour. If we switch to
                ! addNewChild this test becomes redundant.
                if (.not. associated(thisOctal, neighbour)) then
                  call addNewChildren(neighbour, grid, sphData, stellar_cluster, inheritProps, interpProps)
                end if 

                gridConverged = .FALSE.
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
  END SUBROUTINE smoothAMRgridIonization
  
  recursive subroutine dumpdensitytemperature(thisOctal, unit)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i, unit

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call dumpdensitytemperature(child, unit)
                exit
             end if
          end do
       else
          write(unit) thisOctal%rho(subcell),thisOctal%temperature(subcell)
       end if
    end do

  end subroutine dumpdensitytemperature

  recursive subroutine dumpdiffusion(thisOctal, unit)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i, unit

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call dumpdiffusion(child, unit)
                exit
             end if
          end do
       else
          write(unit) thisOctal%diffusionApprox(subcell)
       end if
    end do

  end subroutine dumpdiffusion

  recursive subroutine myTauSmooth(thisOctal, grid, ilambda, converged, inheritProps, interpProps)
    use input_variables, only : tauSmoothMin, tauSmoothMax, erOuter, drouter, router, rinner
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    logical, optional :: inheritProps, interpProps
    integer :: subcell, i, ilambda
    logical :: converged
    real(double) :: kabs, ksca, r
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat, rVec
    real :: thisTau, neighbourTau
    integer :: neighbourSubcell, j, nDir
    logical :: split
    logical, save :: firsttime = .true.
    integer, parameter :: maxDepth = 28
    character(len=30) :: message

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call myTauSmooth(child, grid, ilambda, converged, inheritProps, interpProps)
                exit
             end if
          end do
       else

          split = .true.

          call returnKappa(grid, thisOctal, subcell, ilambda=ilambda,&
               kappaSca=ksca, kappaAbs=kabs)

          thisTau  = thisOctal%subcellSize * (ksca + kabs)

          r = thisOctal%subcellSize/2. + grid%halfSmallestSubcell * 0.1d0
          centre = subcellCentre(thisOctal, subcell) + &
               (0.01d0*grid%halfSmallestsubcell)*OCTALVECTOR(+1.d0,+1.d0,+1.d0)
          if (.not.thisOctal%cylindrical) then
             nDir = 6
             dirVec(1) = OCTALVECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = OCTALVECTOR( 0.d0,+1.d0,  0.d0)
             dirVec(3) = OCTALVECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(4) = OCTALVECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(5) = OCTALVECTOR( 0.d0,-1.d0,  0.d0)
             dirVec(6) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 4
             aHat = OCTALVECTOR(centre%x, centre%y, 0.d0)
             call normalize(aHat)

             dirVec(1) = OCTALVECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = aHat
             dirVec(3) = (-1.d0)*aHat
             dirVec(4) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
          endif
          do j = 1, nDir
             octVec = centre + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then
                startOctal => thisOctal
                call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                     foundOctal=neighbourOctal, foundsubcell=neighbourSubcell, kappaSca=ksca, kappaAbs=kabs, ilambda=ilambda)
                neighbourTau = thisOctal%subcellSize * (ksca + kabs)

                if ((grid%geometry.eq."whitney").and.&
                     (modulus(subcellCentre(thisOctal,subcell)) > 0.9*erouter/1.e10)) split = .false.

                rVec = subcellCentre(thisOctal,subcell)

                if ((grid%geometry.eq."warpeddisc").and.&
                     (sqrt(rVec%x**2+rVec%y**2)) > 0.9*router) split = .false.


                if (thisOctal%nDepth == maxDepth) then
                   split = .false.
                   if (firstTime) then
                      write(message,'(a,i3)') "AMR cell depth capped at: ",maxDepth
                      call writeWarning(message)
                      firstTime = .false.
                   endif
                endif



                if ((min(thisTau, neighbourTau) < tauSmoothMin).and.(max(thisTau, neighbourTau) > tauSmoothMax).and.split) then
!                   write(*,*) modulus(subcellCentre(neighbourOctal, neighboursubcell)),thisOctal%ndepth,neighbourOctal%ndepth
                   if (thisTau > neighbourTau) then
                      call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                           inherit=inheritProps, interp=interpProps)
!                      write(*,*) thisOctal%nDepth,thisTau, neighbourTau,sqrt(rVec%x**2+rVec%y**2)/rinner,&
!                           sqrt(rVec%x**2+rVec%y**2)/router
                      converged = .false.
                      return
                   endif
                endif
             endif
          enddo
       endif
    end do

  end subroutine myTauSmooth

  recursive subroutine myScaleSmooth(factor, thisOctal, grid,  converged, inheritProps, interpProps)
    type(gridtype) :: grid
    real :: factor
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child, neighbourOctal, startOctal
    logical, optional :: inheritProps, interpProps
    integer :: subcell, i, ilambda
    logical :: converged
    real(double) :: kabs, ksca, r
    type(OCTALVECTOR) :: dirVec(6), centre, octVec, aHat
    real :: thisTau, neighbourTau
    integer :: neighbourSubcell, j, nDir



    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call myScaleSmooth(factor, child, grid, converged, inheritProps, interpProps)
                exit
             end if
          end do
       else

          r = thisOctal%subcellSize/2.d0 + 0.1d0*grid%halfSmallestSubcell
          centre = subcellCentre(thisOctal, subcell) + &
               (0.01d0*grid%halfSmallestsubcell)*OCTALVECTOR(+1.d0,+1.d0,+1.d0)
          if (.not.thisOctal%cylindrical) then
             nDir = 6
             dirVec(1) = OCTALVECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = OCTALVECTOR( 0.d0,+1.d0,  0.d0)
             dirVec(3) = OCTALVECTOR(+1.d0, 0.d0,  0.d0)
             dirVec(4) = OCTALVECTOR(-1.d0, 0.d0,  0.d0)
             dirVec(5) = OCTALVECTOR( 0.d0,-1.d0,  0.d0)
             dirVec(6) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
          else
             nDir = 4
             aHat = OCTALVECTOR(centre%x, centre%y, 0.d0)
             call normalize(aHat)

             dirVec(1) = OCTALVECTOR( 0.d0, 0.d0, +1.d0)
             dirVec(2) = aHat
             dirVec(3) = (-1.d0)*aHat
             dirVec(4) = OCTALVECTOR( 0.d0, 0.d0, -1.d0)
          endif
          do j = 1, nDir
             octVec = centre + r * dirvec(j)
             if (inOctal(grid%octreeRoot, octVec)) then

                startOctal => thisOctal
                call amrGridValues(grid%octreeRoot, octVec, grid=grid, startOctal=startOctal, &
                     foundOctal=neighbourOctal, foundsubcell=neighbourSubcell)

!                if ((thisOctal%nDepth < (neighbourOctal%ndepth-1))) then
                if ((thisOctal%subcellSize/neighbourOctal%subcellSize) > factor) then
                      call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
                           inherit=inheritProps, interp=interpProps)
                      converged = .false.
                      return
                endif

             endif
          enddo
       endif
    end do

  end subroutine myScaleSmooth


  
  subroutine distanceToCellBoundary(grid, posVec, direction, tVal, sOctal)


   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTALVECTOR), intent(in) :: posVec
   type(OCTALVECTOR), intent(in) :: direction
   type(OCTAL), pointer, optional :: sOctal
   real(oct), intent(out) :: tval
   !
   type(OCTALVECTOR) :: norm(6), p3(6), thisNorm
   type(OCTAL),pointer :: thisOctal
   real(double) :: distTor1, distTor2, theta, mu
   real(double) :: distToRboundary, compz,currentZ
   real(double) :: phi, distToZboundary, ang1, ang2
   type(OCTALVECTOR) :: subcen, point, xHat, zHat, rVec
   integer :: subcell
   real(double) :: distToSide1, distToSide2, distToSide
   real(double) ::  compx,disttoxBoundary, currentX
   real(oct) :: t(6),denom(6), r, r1, r2, d, cosmu,x1,x2
   integer :: i,j
   logical :: ok, thisOk(6)


   point = posVec

   if (PRESENT(sOctal)) then
      call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid, startOctal=sOctal)
   else
      call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
   endif
   subcen =  subcellCentre(thisOctal,subcell)

   if (thisOctal%threed) then

      if (.not.thisOctal%cylindrical) then
         ok = .true.
         
         norm(1) = OCTALVECTOR(1.0d0, 0.d0, 0.0d0)
         norm(2) = OCTALVECTOR(0.0d0, 1.0d0, 0.0d0)
         norm(3) = OCTALVECTOR(0.0d0, 0.0d0, 1.0d0)
         norm(4) = OCTALVECTOR(-1.0d0, 0.0d0, 0.0d0)
         norm(5) = OCTALVECTOR(0.0d0, -1.0d0, 0.0d0)
         norm(6) = OCTALVECTOR(0.0d0, 0.0d0, -1.0d0)
         
         p3(1) = OCTALVECTOR(subcen%x+thisOctal%subcellsize/2.0d0, subcen%y, subcen%z)
         p3(2) = OCTALVECTOR(subcen%x, subcen%y+thisOctal%subcellsize/2.0d0 ,subcen%z)
         p3(3) = OCTALVECTOR(subcen%x,subcen%y,subcen%z+thisOctal%subcellsize/2.0d0)
         p3(4) = OCTALVECTOR(subcen%x-thisOctal%subcellsize/2.0d0, subcen%y,  subcen%z)
         p3(5) = OCTALVECTOR(subcen%x,subcen%y-thisOctal%subcellsize/2.0d0, subcen%z)
         p3(6) = OCTALVECTOR(subcen%x,subcen%y,subcen%z-thisOctal%subcellsize/2.0d0)

         thisOk = .true.
         
         do i = 1, 6
            
            denom(i) = norm(i) .dot. direction
            if (denom(i) /= 0.0d0) then
               t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
            else
               thisOk(i) = .false.
               t(i) = 0.0d0
            endif
            if (t(i) < 0.) thisOk(i) = .false.
            !      if (denom > 0.) thisOK(i) = .false.
         enddo
         
         
         j = 0
         do i = 1, 6
            if (thisOk(i)) j=j+1
         enddo
         
         if (j == 0) ok = .false.
         
         if (.not.ok) then
            write(*,*) "Error: j=0 (no intersection???) in lucy_mod::intersectCubeAMR. "
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
         
         if (tval > sqrt(3.)*thisOctal%subcellsize) then
            !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
            !     write(*,*) "direction",direction
            !     write(*,*) t(1:6)
            !     write(*,*) denom(1:6)
         endif

      else

! now look at the cylindrical case

         ! first do the inside and outside curved surfaces
         r = sqrt(subcen%x**2 + subcen%y**2)
         r1 = r - thisOctal%subcellSize/2.d0
         r2 = r + thisOctal%subcellSize/2.d0
         d = sqrt(point%x**2+point%y**2)
         xHat = VECTOR(point%x, point%y,0.d0)
         call normalize(xHat)
      
         cosmu =((-1.d0)*xHat).dot.direction
         call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
         if (.not.ok) then
            write(*,*) "Quad solver failed in intersectcubeamr2d"
            x1 = thisoctal%subcellSize/2.d0
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
               x1 = thisoctal%subcellSize/2.d0
               x2 = 0.d0
            endif
            distTor1 = max(x1,x2)
         endif
      
         distToRboundary = min(distTor1, distTor2)
         if (distToRboundary < 0.d0) then
            distToRboundary = 1.e30
         endif

         ! now do the upper and lower (z axis) surfaces
      

         zHat = VECTOR(0.d0, 0.d0, 1.d0)
         compZ = zHat.dot.direction
         currentZ = point%z
      
         if (compZ /= 0.d0 ) then
            if (compZ > 0.d0) then
               distToZboundary = (subcen%z + thisOctal%subcellsize/2.d0 - currentZ ) / compZ
            else
               distToZboundary = abs((subcen%z - thisOctal%subcellsize/2.d0 - currentZ ) / compZ)
            endif
         else
            disttoZboundary = 1.e30
         endif
      
        
         ! ok now we have to tackle the two angled sides...

         ! find posvec to surface centre

         phi = atan2(posVec%y,posVec%x)
         if (phi < 0.d0) phi = phi + twoPi

         rVec = OCTALVECTOR(r, 0.d0, 0.d0)

         ang1 = thisOctal%phi - returndPhi(thisOctal)
         ang2 = thisOctal%phi + returndPhi(thisOctal)
!         if (thisOctal%splitAzimuthally) then
!            if (phi < thisOctal%phi) then
!               ang1 = thisOctal%phi - thisOctal%dPhi/2.d0
!               ang2 = thisOctal%phi
!            else
!               ang1 =  thisOctal%phi
!               ang2 = thisOctal%phi + thisOctal%dPhi/2.d0
!            endif
!         else
!            ang1 = thisOctal%phi - thisOctal%dPhi/2.d0
!            ang2 = thisOctal%phi + thisOctal%dPhi/2.d0
!         endif

         distToSide1 = 1.e30
         rVec = OCTALVECTOR(r, 0.d0, 0.d0)
         rVec = rotateZ(rVec, -ang1)
         thisnorm = rVec .cross. zHat
         call normalize(thisnorm)
         if ((thisnorm.dot.direction) /= 0.d0) then
            distToSide1 = (thisnorm.dot.(rVec-posVec))/(thisnorm.dot.direction)
            if (distToSide1 < 0.d0) distToSide1 = 1.d30
         endif

         distToSide2 = 1.e30
         rVec = OCTALVECTOR(r, 0.d0, 0.d0)
         rVec = rotateZ(rVec, -ang2)
         thisnorm = rVec .cross. zHat
         call normalize(thisnorm)
         if ((thisnorm.dot.direction) /= 0.d0) then
            distToSide2 = (thisnorm.dot.(rVec-posVec))/(thisnorm.dot.direction)
            if (distToSide2 < 0.d0) distToSide2 = 1.d30
         endif

         distToSide = min(distToSide1, distToside2)

         tVal = min(distToZboundary, distToRboundary, distToSide)
         if (tVal > 1.e29) then
            write(*,*) "Cylindrical"
            write(*,*) tVal,compX,compZ, distToZboundary,disttorboundary, disttoside
            write(*,*) "subcen",subcen
            write(*,*) "x,z",currentX,currentZ
         endif
         if (tval < 0.) then
            write(*,*) "Cylindrical"
            write(*,*) tVal,distToZboundary,disttorboundary, disttoside
            write(*,*) "subcen",subcen
            write(*,*) "x,z",currentX,currentZ
         endif

      endif

   else ! two-d grid case below

      r1 = subcen%x - thisOctal%subcellSize/2.d0
      r2 = subcen%x + thisOctal%subcellSize/2.d0
      d = sqrt(point%x**2+point%y**2)
      xHat = VECTOR(point%x, point%y,0.d0)
      call normalize(xHat)
      
      cosmu =((-1.d0)*xHat).dot.direction
      call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
      if (.not.ok) then
         write(*,*) "Quad solver failed in intersectcubeamr2d"
         x1 = thisoctal%subcellSize/2.d0
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
            x1 = thisoctal%subcellSize/2.d0
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
            distToZboundary = (subcen%z + thisOctal%subcellsize/2.d0 - currentZ ) / compZ
         else
            distToZboundary = abs((subcen%z - thisOctal%subcellsize/2.d0 - currentZ ) / compZ)
         endif
      else
         disttoZboundary = 1.e30
      endif
      
      tVal = min(distToZboundary, distToXboundary)
      if (tVal > 1.e29) then
         write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
         write(*,*) "subcen",subcen
         write(*,*) "x,z",currentX,currentZ
      endif
      if (tval < 0.) then
         write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
         write(*,*) "subcen",subcen
         write(*,*) "x,z",currentX,currentZ
      endif
      
   endif

   tVal = max(tVal, 0.001d0*grid%halfSmallestSubcell) ! avoid sticking on a cell boundary


 end subroutine distanceToCellBoundary

  subroutine distanceToGridEdge(grid, posVec, direction, tVal)

   implicit none
   type(GRIDTYPE), intent(in)    :: grid
   type(OCTALVECTOR), intent(in) :: posVec
   type(OCTALVECTOR), intent(in) :: direction
   real(oct), intent(out) :: tval
   !
   type(OCTALVECTOR) :: norm(6), p3(6), thisNorm
   real(double) :: distTor1, distTor2, theta, mu
   real(double) :: distToRboundary, compz,currentZ
   real(double) :: phi, distToZboundary, ang1, ang2
   type(OCTALVECTOR) :: subcen, point, xHat, zHat, rVec
   integer :: subcell
   real(double) :: distToSide1, distToSide2, distToSide
   real(double) ::  compx,disttoxBoundary, currentX, subcellsize
   real(oct) :: t(6),denom(6), r, r1, r2, d, cosmu,x1,x2
   integer :: i,j
   logical :: ok, thisOk(6)


   point = posVec

   subcen =  grid%octreeRoot%centre
   subcellsize = grid%octreeRoot%subcellSize * 2.d0

   if (grid%octreeRoot%threed) then

      if (.not.grid%octreeRoot%cylindrical) then
         ok = .true.
         
         norm(1) = OCTALVECTOR(1.0d0, 0.d0, 0.0d0)
         norm(2) = OCTALVECTOR(0.0d0, 1.0d0, 0.0d0)
         norm(3) = OCTALVECTOR(0.0d0, 0.0d0, 1.0d0)
         norm(4) = OCTALVECTOR(-1.0d0, 0.0d0, 0.0d0)
         norm(5) = OCTALVECTOR(0.0d0, -1.0d0, 0.0d0)
         norm(6) = OCTALVECTOR(0.0d0, 0.0d0, -1.0d0)
         
         p3(1) = OCTALVECTOR(subcen%x+subcellsize/2.0d0, subcen%y, subcen%z)
         p3(2) = OCTALVECTOR(subcen%x, subcen%y+subcellsize/2.0d0 ,subcen%z)
         p3(3) = OCTALVECTOR(subcen%x,subcen%y,subcen%z+subcellsize/2.0d0)
         p3(4) = OCTALVECTOR(subcen%x-subcellsize/2.0d0, subcen%y,  subcen%z)
         p3(5) = OCTALVECTOR(subcen%x,subcen%y-subcellsize/2.0d0, subcen%z)
         p3(6) = OCTALVECTOR(subcen%x,subcen%y,subcen%z-subcellsize/2.0d0)

         thisOk = .true.
         
         do i = 1, 6
            
            denom(i) = norm(i) .dot. direction
            if (denom(i) /= 0.0d0) then
               t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
            else
               thisOk(i) = .false.
               t(i) = 0.0d0
            endif
            if (t(i) < 0.) thisOk(i) = .false.
            !      if (denom > 0.) thisOK(i) = .false.
         enddo
         
         
         j = 0
         do i = 1, 6
            if (thisOk(i)) j=j+1
         enddo
         
         if (j == 0) ok = .false.
         
         if (.not.ok) then
            write(*,*) "Error: j=0 (no intersection???) in lucy_mod::intersectCubeAMR. "
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
!         rVec = OCTALVECTOR(r, 0.d0, 0.d0)
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
!         rVec = OCTALVECTOR(r, 0.d0, 0.d0)
!         rVec = rotateZ(rVec, -ang1)
!         thisnorm = rVec .cross. zHat
!         call normalize(thisnorm)
!         if ((thisnorm.dot.direction) /= 0.d0) then
!            distToSide1 = (thisnorm.dot.(rVec-posVec))/(thisnorm.dot.direction)
!            if (distToSide1 < 0.d0) distToSide1 = 1.d30
!         endif
!
!         rVec = OCTALVECTOR(r, 0.d0, 0.d0)
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
            write(*,*) tVal,compX,compZ, distToZboundary,disttorboundary, disttoside
            write(*,*) "subcen",subcen
            write(*,*) "x,z",currentX,currentZ
         endif
         if (tval < 0.) then
            write(*,*) "Cylindrical"
            write(*,*) tVal,distToZboundary,disttorboundary, disttoside
            write(*,*) "subcen",subcen
            write(*,*) "x,z",currentX,currentZ
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
         write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
         write(*,*) "subcen",subcen
         write(*,*) "x,z",currentX,currentZ
      endif
      if (tval < 0.) then
         write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
         write(*,*) "subcen",subcen
         write(*,*) "x,z",currentX,currentZ
      endif
      
   endif

   tVal = max(tVal, 1.d-4*grid%halfSmallestSubcell) ! avoid sticking on a cell boundary

 end subroutine distanceToGridEdge


  subroutine writeCylindrical(thisOctal, subcell)
    type(OCTAL) :: thisOctal
    integer :: subcell
    type(OCTALVECTOR) :: octVec
    real(double) :: r, z, phi

    octVec = subcellCentre(thisoctal,subcell)

    phi = atan2(octVec%y, octVec%x)
    phi = phi + twoPi
    r = sqrt(octVec%x**2 + octVec%y**2)

    write(*,*) "r",r,thisOctal%r
    write(*,*) "phi",phi*radtodeg
    write(*,*) "z",octVec%z

  end subroutine writeCylindrical


  type(OCTALVECTOR) function randomPositionInCell(thisOctal, subcell)
    type(OCTAL) :: thisOctal
    integer :: subcell
    type(OCTALVECTOR) :: octalCentre
    real(double) :: r1, r2, r3
    real(double) :: xOctal, yOctal, zOctal
    real(double) :: ang, ang1, ang2, phi


    octalCentre = subcellCentre(thisOctal,subcell)
    
!!! we will just choose a random point within the subcell.
!!! this *should* be done in a better way.

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
  end function randomPositionInCell




  SUBROUTINE startReturnSamples2(startPoint,direction,grid,          &
             sampleFreq,nSamples,maxSamples,thin_disc_on, opaqueCore,hitCore,      &
             usePops,iLambda,error,lambda,nSource,source,kappaAbs,kappaSca,velocity,&
             velocityDeriv,chiLine,levelPop,rho, temperature, Ne, inflow)
    ! samples the grid at points along the path.
    ! this should be called by the program, instead of calling 
    !   returnSamples directly, because this checks the start and finish
    !   points are within the grid bounds. returnSamples *assumes* this 
    !   criterion is met. also, the path of the photon is checked for 
    !   intersection with the stellar surface(s), disc etc. 
 
    IMPLICIT NONE
    integer :: nSource
    type(SOURCETYPE) :: source(:)

    TYPE(octalVector), INTENT(IN)      :: startPoint ! photon start point
    TYPE(octalVector), INTENT(IN)      :: direction  ! photon direction 
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

    TYPE(octalVector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    TYPE(octalVector)       :: currentPoint ! current position of ray 
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(octal)             :: octree       ! the octree structure within 'grid'
    TYPE(octalVector)       :: directionNormalized
    real(oct)    :: distanceLimit ! max length of ray before aborting
    ! margin is the size of the region around the edge of a subcell
    !   where numerical inaccuracies may cause problems.
    real(oct)    :: margin
 
    ! variables for testing special cases (stellar intersections etc.)
    TYPE(octalVector)       :: starPosition       ! position vector of stellar centre
    TYPE(octalVector)       :: diskNormal         ! disk normal vector
    real(oct)    :: rStar              ! stellar radius
    real(oct)    :: endLength          ! max path length of photon
    TYPE(octalVector)       :: endPoint           ! where photon leaves grid
    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
    real(oct)    :: distanceFromOrigin ! closest distance of plane from origin
    TYPE(octalVector)       :: diskIntersection   ! point of photon intersection with disk
    real(oct)    :: starIntersectionDistance1 ! distance to first  intersection
    real(oct)    :: starIntersectionDistance2 ! distance to second intersection
    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
    real(oct)    :: intersectionRadius ! disk radius when photon intersects
    real(oct)    :: diskDistance       ! distance to disk intersection
    real(oct)    :: distanceThroughStar! distance of chord through star
    TYPE(octalVector)       :: dummyStartPoint    ! modified start point 
    real(oct), PARAMETER :: fudgefactor = 1.00001 ! overestimates stellar size
    
    type(octalvector) :: currentPosition
    type(octal), pointer :: sOctal, thisOctal
    integer :: subcell
    real(double) :: length, distToNextCell
    ! we will abort tracking a photon just before it reaches the edge of the
    !   simulation space. This is the fraction of the total distance to use:
    real(oct), PARAMETER :: distanceFraction = 0.999_oc 
    real(double) :: fudgeFac = 0.001d0


    currentPosition = startPoint
    locator = startPoint
    nSamples = 0
    length = 0.d0
    endLength = 1.d30

    absorbPhoton = .false.
    hitcore = .false.

    if (nSource == 0) then
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
       call distanceToSource(source, nSource, startPoint, direction, absorbPhoton, endLength)
    endif




    CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)

    do while (inOctal(grid%octreeRoot, currentPosition).and.(length < endLength))

       call findSubcellLocal(currentPosition,thisOctal,subcell)

       sOctal => thisOctal
       call distanceToCellBoundary(grid, currentPosition, direction, DisttoNextCell, sOctal)

       if (nSamples < maxSamples) then
          call takeSample(currentPosition,length,direction,grid,thisOctal,subcell,nSamples,&
               maxSamples,usePops,iLambda,error,lambda,kappaAbs,     &
               kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,  &
               temperature,Ne,inFlow) 
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
    END IF



  end SUBROUTINE startReturnSamples2



END MODULE amr_mod
