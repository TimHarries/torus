MODULE amr_mod
  ! routines for adaptive mesh refinement. nhs
  ! twod stuff added by tjh started 25/08/04


  USE vector_mod          ! vector maths routines
  USE kind_mod            ! variable kind parameters    
  USE octal_mod           ! type definition for the grid elements
  USE gridtype_mod        ! type definition for the 3-d grid 
  USE parameters_mod      ! parameters for specific geometries
  USE jets_mod            ! 
  USE constants_mod, only: cSpeed
  USE sph_data_class
  USE cluster_class
  USE density_mod
  USE wr104_mod
  USE utils_mod

  IMPLICIT NONE


CONTAINS


  SUBROUTINE calcValuesAMR(thisOctal,subcell,grid, sphData, stellar_cluster)
    ! calculates the variables describing one subcell of an octal.
    ! each geometry that can be used with AMR should be described here, 
    !   otherwise the program will print a warning and exit.
   
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT)    :: thisOctal ! the octal being changed
    INTEGER, INTENT(IN)           :: subcell   ! the subcell being changed
    TYPE(gridtype), INTENT(INOUT) :: grid      ! the grid
    !
    TYPE(sph_data), optional, INTENT(IN)  :: sphData   ! Matthew's SPH data.
    TYPE(cluster), optional, INTENT(IN)   :: stellar_cluster

    SELECT CASE (grid%geometry)

    CASE ("ttauri")
      CALL calcTTauriMassVelocity(thisOctal,subcell,grid)
      
    CASE ("jets")
      CALL calcJetsMassVelocity(thisOctal,subcell,grid)

    CASE ("testamr")
       CALL calcTestDensity(thisOctal,subcell,grid)

    CASE("proto")
       CALL calcProtoDensity(thisOctal,subcell,grid)

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

    CASE("melvin")
       CALL assign_melvin(thisOctal,subcell,grid)

    CASE("clumpydisc")
       CALL assign_clumpydisc(thisOctal, subcell, grid)
       
    CASE DEFAULT
      WRITE(*,*) "! Unrecognised grid geometry: ",TRIM(grid%geometry)
      STOP

    END SELECT
 
!    CALL fillGridDummyValues(thisOctal,subcell, grid)
 
 
  END SUBROUTINE calcValuesAMR


  SUBROUTINE initFirstOctal(grid, centre, size, twod, sphData, stellar_cluster, nDustType)
    ! creates the first octal of a new grid (the root of the tree).
    ! this should only be used once; use addNewChild for subsequent
    !  additions.

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
       ALLOCATE(grid%octreeRoot%kappaAbs(8,grid%nLambda))
       ALLOCATE(grid%octreeRoot%kappaSca(8,grid%nLambda))
       grid%octreeRoot%kappaAbs = 1.e-30
       grid%octreeRoot%kappaSca = 1.e-30
!    ELSE
!       ALLOCATE(grid%oneKappaAbs(nDustType,grid%nLambda))
!       ALLOCATE(grid%oneKappaSca(nDustTYpe,grid%nLambda))
!       grid%oneKappaAbs = 1.e-30
!       grid%oneKappaSca = 1.e-30
    endif

    ALLOCATE(grid%octreeRoot%N(8,grid%maxLevels))

    if (twoD) then
       grid%octreeRoot%twoD = .true.
       grid%octreeRoot%threeD = .false.
       grid%octreeRoot%maxChildren = 4
    else
       grid%octreeRoot%twoD = .false.
       grid%octreeRoot%threeD = .true.
       grid%octreeRoot%maxChildren = 8
    endif
    grid%octreeRoot%nDepth = 1
    grid%octreeRoot%nChildren = 0
    grid%octreeRoot%hasChild = .FALSE.
    grid%octreeRoot%subcellSize = size/2.0_oc
    grid%octreeRoot%centre = centre
    grid%octreeRoot%indexChild = -999 ! values are undefined
    grid%octreeRoot%probDistLine = 0.0
    grid%octreeRoot%probDistCont = 0.0
    NULLIFY(grid%octreeRoot%parent)   ! tree root does not have a parent
    NULLIFY(grid%octreeRoot%child)    ! tree root does not yet have children
    ! initialize some values
    grid%octreeRoot%biasLine3D = 1.0 
    grid%octreeRoot%biasCont3D = 1.0 
    grid%octreeRoot%velocity = vector(1.e-30,1.e-30,1.e-30)
    grid%octreeRoot%cornerVelocity = vector(1.e-30,1.e-30,1.e-30)
    grid%octreeRoot%chiLine = 1.e-30
    grid%octreeRoot%etaLine = 1.e-30
    grid%octreeRoot%etaCont = 1.e-30
    grid%octreeRoot%N = 1.e-30

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
             CALL calcValuesAMR(grid%octreeRoot,subcell,grid, sphData)
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
                                2.0_oc**REAL(grid%maxDepth,KIND=octalKind)
      
  END SUBROUTINE initFirstOctal


! THIS ROUTINE IS NOW REDUNDANT AND HAS BEEN REPLACED BY ADDNEWCHILDEN!!!!!!

  SUBROUTINE addNewChild(parent, nChild, grid, sphData, stellar_cluster)
    ! adds one new child to an octal

    IMPLICIT NONE
    
    TYPE(octal), POINTER :: parent     ! pointer to the parent octal 
    INTEGER, INTENT(IN)  :: nChild     ! the label (1-8) of the subcell gaining the child 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that
                                          !   calculate the variables stored in
                                          !   the tree.
    TYPE(octal)   :: tempChildStorage  ! holder for existing children, while we
                                       !   shuffle them around to make room for 
                                       !   the new child.
    INTEGER       :: subcell           ! loop counter
    INTEGER, SAVE :: counter = 9       ! keeps track of the current subcell label
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: nChildren         ! number of children the parent octal has
    INTEGER       :: newChildIndex     ! the storage location for the new child
    INTEGER       :: error
    


    ! For only cluster geometry ...
    TYPE(sph_data), optional, intent(in) :: sphData 
    TYPE(cluster), optional, INTENT(IN)   :: stellar_cluster

    ! store the number of children that already exist
    nChildren = parent%nChildren

    IF ( nChildren == 0 ) THEN
      ! if there are no existing children, we can just allocate
      ! the 'child' array with one item
      ALLOCATE(parent%child(1), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed.'
        STOP
      END IF

    ELSE 
      ! check that child does not already exist
      IF ( parent%hasChild(nChild) .EQV. .TRUE. ) THEN
        PRINT *, 'Panic: in addNewChild, attempted to add a child ',&
                 '       that already exists'
        STOP
      ENDIF
      
      ! if there are existing children, we must enlarge the allocated   
      ! array. we need to use a temporary octal structure[1] and copy the  
      ! existing children into it[2]; then increase the 'child' array size 
      ! by one[3]; then copy the children back in[4]. 

      ! [1]
      ALLOCATE(tempChildStorage%child(1:nChildren), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed.'
        STOP
      END IF
      
      ![2]
      tempChildStorage%child(1:nChildren) = parent%child(1:nChildren)

      
!      ![3] 
!
      ALLOCATE(parent%child(nChildren+1), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed.'
        STOP
      END IF

      ![4]
      parent%child(1:nChildren) = tempChildStorage%child(1:nChildren)   

      
      ! we should now release the temporary storage
      DEALLOCATE(tempChildStorage%child, STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: deallocation failed.'
        STOP
      END IF
      NULLIFY(tempChildStorage%child)  ! For safety
    ENDIF


    ! update the bookkeeping
    nChildren = nChildren + 1
    newChildIndex = nChildren

    ! update the parent octal
    parent%nChildren = nChildren
    parent%hasChild(nChild) = .TRUE.
    parent%indexChild(newChildIndex) = nChild

    ! allocate any variables that need to be  
    if (.not.grid%oneKappa) then
       ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nLambda))
       parent%child(newChildIndex)%kappaAbs = 1.e-30
       ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nLambda))
       parent%child(newChildIndex)%kappaSca = 1.e-30
    endif
    ALLOCATE(parent%child(newChildIndex)%N(8,grid%maxLevels))
    NULLIFY(parent%child(newChildIndex)%child)

    ! set up the new child's variables
    parent%child(newChildIndex)%parent => parent
    parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
    parent%child(newChildIndex)%hasChild = .false.
    parent%child(newChildIndex)%nChildren = 0
    parent%child(newChildIndex)%indexChild = -999 ! values are undefined
    parent%child(newChildIndex)%nDepth = parent%nDepth + 1
    parent%child(newChildIndex)%centre = subcellCentre(parent,nChild)
    parent%child(newChildIndex)%probDistLine = 0.0
    parent%child(newChildIndex)%probDistCont = 0.0
    parent%child(newChildIndex)%biasLine3D = 1.0 
    parent%child(newChildIndex)%biasCont3D = 1.0 
    parent%child(newChildIndex)%velocity = vector(1.e-30,1.e-30,1.e-30)
    parent%child(newChildIndex)%cornerVelocity = vector(1.e-30,1.e-30,1.e-30)
    parent%child(newChildIndex)%chiLine = 1.e-30
    parent%child(newChildIndex)%etaLine = 1.e-30
    parent%child(newChildIndex)%etaCont = 1.e-30
    parent%child(newChildIndex)%N = 1.e-30


    select case (grid%geometry)
       case("cluster","wr104")
          if (present(sphData))  then
             ! updates the sph particle linked list.           
             call update_particle_list(parent, nChild, newChildIndex, sphData)
             
             ! put some data in the eight subcells of the new child
             DO subcell = 1, parent%maxChildren 
                CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid, &
                     sphData, stellar_cluster)
                parent%child(newChildIndex)%label(subcell) = counter
                counter = counter + 1
             END DO
          endif
       case DEFAULT
          ! put some data in the eight subcells of the new child
          DO subcell = 1, parent%maxChildren
             CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid)
             parent%child(newChildIndex)%label(subcell) = counter
             counter = counter + 1
          END DO
    end select
    
 
    ! check for a new maximum depth 
    IF (parent%child(newChildIndex)%nDepth > grid%maxDepth) THEN
      grid%maxDepth = parent%child(newChildIndex)%nDepth
      grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
                                2.0_oc**REAL(grid%maxDepth,KIND=octalKind)
      ! we store the value which is half the size of the 
      !   smallest subcell because this is more useful for later
      !   calculations.
    END IF

  END SUBROUTINE addNewChild


  

  RECURSIVE SUBROUTINE splitGrid(thisOctal,amrLimitScalar,amrLimitScalar2,grid, &
       sphData, stellar_cluster)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    IMPLICIT NONE

    TYPE(OCTAL), POINTER :: thisOctal
    REAL(KIND=doubleKind), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
                                          !   routines that this subroutine calls

    ! Object containg the output from the (Mattew's) SPH code.
    ! This will be just passed to decideSplit routine.
    TYPE(sph_data), OPTIONAL, intent(in) :: sphData
    TYPE(cluster), OPTIONAL,  intent(in) :: stellar_cluster
    !
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    logical :: splitThis

    splitThis = .false.
    subcell = 1
    do while ((.not.splitThis).and.(subcell <= thisOctal%maxChildren))
       IF (decideSplit(thisOctal,subcell,amrLimitScalar,amrLimitScalar2,grid,&
            sphData, stellar_cluster))splitThis = .true.
       subcell = subcell + 1
    enddo

         ! tjh changed from
!        CALL addNewChild(thisOctal,subcell,grid)
    IF (splitThis) then
       CALL addNewChildren(thisOctal, grid, sphData, stellar_cluster)
       ! find the index of the new child and call splitGrid on it
        DO i = 1, thisOctal%maxChildren, 1
           childPointer => thisOctal%child(i)
           CALL splitGrid(childPointer,amrLimitScalar,amrLimitScalar2,grid,&
                sphData, stellar_cluster )
        END DO
     ENDIF
      

  END SUBROUTINE splitGrid

  
  
  RECURSIVE SUBROUTINE getOctalArray(thisOctal,array,counter) 
    ! returns an array of pointers to all of the subcells in the grid.
    ! NB because fortran cannot create arrays of pointers, the output
    !   array is actually of a derived type which *contains* the 
    !   pointer to an octal.
    ! counter should be set to 0 before this routine is called

    IMPLICIT NONE

    TYPE(octal), POINTER                            :: thisOctal
    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    INTEGER, INTENT(INOUT)                          :: counter 
  
    INTEGER              :: i
    TYPE(octal), POINTER :: child

    ! if this is the root of the tree, we initialize the counter
    IF (.NOT. ASSOCIATED(thisOctal%parent)) counter = 0
    
    counter = counter + 1 
    array(counter)%content => thisOctal
    !array(counter)%inUse = .TRUE. 
    array(counter)%inUse = .NOT. thisOctal%hasChild 
    
    IF ( thisOctal%nChildren > 0 ) THEN
      DO i = 1, thisOctal%nChildren, 1
        
        ! call this subroutine recursively on each of its children
        child => thisOctal%child(i)
        CALL getOctalArray(child,array,counter)
        
      END DO
    END IF

  END SUBROUTINE getOctalArray
    
  

  RECURSIVE SUBROUTINE finishGrid(thisOctal,grid,gridConverged)
    ! takes the octree grid that has been created using 'splitGrid'
    !   and calculates all the other variables in the model.
    ! this should be called once the structure of the grid is complete.
    ! the gridConverged variable should be set .TRUE. when the entire
    !   grid has been finished. Until that happens, this subroutine 
    !   will be called repeatedly.
    
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
        CALL calcTTauriTemperature(thisOctal,subcell)
        gridConverged = .TRUE.
        
      CASE ("spiralwind")
        gridConverged = .TRUE.

      CASE ("jets")
        CALL calcJetsTemperature(thisOctal,subcell, grid)
        gridConverged = .TRUE.
        
      CASE ("testamr","proto")
        gridConverged = .TRUE.

      CASE("benchmark","shakara","aksco", "melvin","clumpydisc")
         gridConverged = .TRUE.

      CASE ("cluster","wr104")
        call assign_grid_values(thisOctal,subcell, grid)
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
             sampleFreq,nSamples,maxSamples,opaqueCore,hitCore,      &
             usePops,iLambda,error,lambda,kappaAbs,kappaSca,velocity,&
             velocityDeriv,chiLine,levelPop,rho, temperature)
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
!    REAL(KIND=octalKind), INTENT(IN)   :: sampleFreq ! the maximum number of 
                       ! samples that will be made in any subcell of the octree. 
                       
    INTEGER, INTENT(OUT)             :: nSamples   ! number of samples made
    INTEGER, INTENT(IN)                :: maxSamples ! size of sample arrays 
    LOGICAL, INTENT(IN)                :: opaqueCore ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: hitCore    ! true if the core is opaque
    LOGICAL, INTENT(IN)                :: usePops    ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda    ! wavelength index
    INTEGER, INTENT(INOUT)             :: error      ! error code
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda     ! path distances of sample locations
    
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs   ! continuous absorption opacities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca   ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: chiLine    ! line opacities
    REAL,DIMENSION(:,:),INTENT(INOUT),OPTIONAL:: levelPop   ! level populations
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: rho        ! density at sample points
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: temperature        ! density at sample points

    TYPE(octalVector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    TYPE(octalVector)       :: currentPoint ! current position of ray 
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(octal)             :: octree       ! the octree structure within 'grid'
    TYPE(octalVector)       :: directionNormalized
    REAL(KIND=octalKind)    :: distanceLimit ! max length of ray before aborting
    ! margin is the size of the region around the edge of a subcell
    !   where numerical inaccuracies may cause problems.
    REAL(KIND=octalKind)    :: margin
 
    ! variables for testing special cases (stellar intersections etc.)
    TYPE(octalVector)       :: starPosition       ! position vector of stellar centre
    TYPE(octalVector)       :: diskNormal         ! disk normal vector
    REAL(KIND=octalKind)    :: rStar              ! stellar radius
    REAL(KIND=octalKind)    :: endLength          ! max path length of photon
    TYPE(octalVector)       :: endPoint           ! where photon leaves grid
    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
    REAL(KIND=octalKind)    :: distanceFromOrigin ! closest distance of plane from origin
    TYPE(octalVector)       :: diskIntersection   ! point of photon intersection with disk
    REAL(KIND=octalKind)    :: starIntersectionDistance1 ! distance to first  intersection
    REAL(KIND=octalKind)    :: starIntersectionDistance2 ! distance to second intersection
    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
    REAL(KIND=octalKind)    :: intersectionRadius ! disk radius when photon intersects
    REAL(KIND=octalKind)    :: diskDistance       ! distance to disk intersection
    REAL(KIND=octalKind)    :: distanceThroughStar! distance of chord through star
    TYPE(octalVector)       :: dummyStartPoint    ! modified start point 
    REAL(KIND=octalKind), PARAMETER :: fudgefactor = 1.00001 ! overestimates stellar size
    
    
    ! we will abort tracking a photon just before it reaches the edge of the
    !   simulation space. This is the fraction of the total distance to use:
    REAL(KIND=octalKind), PARAMETER :: distanceFraction = 0.999_oc 
    
    ! we will abort tracking any rays which are too close to 
    !   to a cell wall. The distance to use is defined by:
    margin = 6.0_oc * REAL(grid%maxDepth,KIND=octalKind) * EPSILON(1.0_oc)
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

      write(*,*) "size",octree%subcellsize
      write(*,*) "centre",octree%centre
      return
!      startpoint%x = min(max(startpoint%x,octree%centre%x - octree%subcellSize),octree%centre%x + octree%subcellSize)
!      startpoint%y = min(max(startpoint%y,octree%centre%y - octree%subcellSize),octree%centre%y + octree%subcellSize)
!      startpoint%z = min(max(startpoint%z,octree%centre%z - octree%subcellSize),octree%centre%z + octree%subcellSize)

    ENDIF
   
   
    ! geometry-specific tests should go here
    IF (grid%geometry(1:6) == "ttauri" .OR. grid%geometry(1:4) == "jets") THEN
      
       ! need to test for both star and disc intersections 
      
       ! we will find out when and where the photon leaves the simulation space 
       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 
       endLength = endLength * distanceFraction
       
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startPoint
       
       ! for the moment, we will just assume a 2-d disc.
       absorbPhoton = .FALSE.
       diskNormal = grid%diskNormal
       starPosition = grid%starPos1
                                  
       ! find the disk intersection point
       distanceFromOrigin = modulus(grid%starPos1)
       diskIntersection = intersectionLinePlane(startPoint, directionNormalized,&
                                 diskNormal, distanceFromOrigin, intersectionFound)
       
       IF (intersectionFound) THEN
       
         ! we need to check whether the photon passes through the disk's
         !   central hole.
         intersectionRadius =  modulus(diskIntersection - starPosition)
         IF (intersectionRadius > grid%diskRadius) THEN
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
                     levelPop=levelPop,rho=rho, temperature=temperature)
      
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
                    chiLine=chiLine,levelPop=levelPop,rho=rho,temperature=temperature)

       END IF

       ! if the photon ends up in the disk or the star, we make sure it absorbed.
       IF (absorbPhoton) THEN
       
         nSamples = nSamples + 1
         IF (nSamples > maxSamples) THEN
           PRINT *, "nSamples > maxSamples in takeSample subroutine"
           STOP
         END IF
         
         kappaSca(nSamples) = 0.
         kappaAbs(nSamples) = 1.e20
         lambda(nSamples) = lambda(nSamples-1)
         hitCore = .TRUE.
         velocity(nSamples) = vector(0.,0.,0.)
         velocityDeriv(nSamples) = 1.0
         chiLine(nSamples) = 1.e20
         levelPop(nSamples,:) = 0.0
         rho(nSamples) = 0.0
       END IF
       
    ELSE
            
       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 
       distanceLimit = endLength * distanceFraction
       
       ! need to reset some of the variables
       abortRay = .FALSE.
       locator = startPoint
       
       CALL returnSamples(currentPoint,startPoint,locator,directionNormalized,&
                   octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,&
                   usePops,iLambda,error,margin,distanceLimit,                &
                   kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,     &
                   velocityDeriv=velocityDeriv,chiLine=chiLine,               &
                   levelPop=levelPop,rho=rho,temperature=temperature)
    END IF
      
  END SUBROUTINE startReturnSamples

  
  RECURSIVE SUBROUTINE returnSamples (currentPoint,startPoint,locator,direction, &
             octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,usePops, &
             iLambda,error,margin,distanceLimit,kappaAbs,kappaSca,velocity,      &
             velocityDeriv,chiLine,levelPop,rho,temperature)
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
    TYPE(octalVector)                   :: locator       
                  ! 'locator' is used to indicate a point that lies within the  
                  !   *next* cell of the octree that the ray will interesect.
                  !   initially this will be the same as the currentPoint
      
    TYPE(octalVector), INTENT(IN)       :: direction    ! ray's direction vector 
    TYPE(octal), INTENT(IN)             :: octree       ! an octal grid
    TYPE(gridtype), INTENT(IN)          :: grid         ! grid structure
    REAL, INTENT(IN)                    :: sampleFreq
!    REAL(KIND=octalKind), INTENT(IN)    :: sampleFreq
                  ! 'sampleFreq' is the maximum number of samples that will be made
                  !   in any subcell of the octree. 
    INTEGER, INTENT(INOUT)              :: nSamples     ! number of samples made
    INTEGER, INTENT(IN)                 :: maxSamples   ! size of sample arrays 
    LOGICAL, INTENT(INOUT)              :: abortRay     ! flag to signal completion
    REAL, DIMENSION(:), INTENT(INOUT)   :: lambda       ! distance travelled by photon
    INTEGER, INTENT(IN)                 :: iLambda      ! wavelength index
    INTEGER, INTENT(INOUT)              :: error        ! error code
    REAL(KIND=octalKind)                :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    REAL(KIND=octalKind), INTENT(IN)    :: distanceLimit ! max length of ray before aborting
    LOGICAL, INTENT(IN)                 :: usePops      ! whether to use level populations
    
    REAL,DIMENSION(:), INTENT(INOUT)   :: kappaAbs     ! continuous absorption opacities
    REAL,DIMENSION(:), INTENT(INOUT)   :: kappaSca     ! scattering opacities
    TYPE(vector), DIMENSION(:), INTENT(INOUT) :: velocity ! sampled velocities
    REAL,DIMENSION(:), INTENT(INOUT)   :: velocityDeriv ! sampled velocity derivatives
    REAL,DIMENSION(:), INTENT(INOUT)   :: chiLine      ! line opacities
    REAL,DIMENSION(:,:), INTENT(INOUT) :: levelPop     ! level populations 
    REAL,DIMENSION(:)                  :: rho          ! density at sample points
    REAL,DIMENSION(:),INTENT(INOUT)    :: temperature  ! density at sample points



    TYPE(octalVector)      :: exitPoint      ! where ray leaves current cell
    TYPE(octalVector)      :: centre         ! centre of current subcell
    REAL(KIND=octalKind)   :: subcellSize    ! size of current subcell
    REAL(KIND=octalKind)   :: minWallDistance ! distance to *nearest* wall
    
    REAL(KIND=octalKind)   :: length         ! distance from start of the ray's path.
    REAL(KIND=octalKind)   :: sampleLength   ! distance interval between samples 
    REAL(KIND=octalKind)   :: trialLength    ! trial distance to a possible next sample
    TYPE(octalVector)      :: trialPoint     ! trial location for a next sample 
    INTEGER                :: trial          ! loop counter for trial points 
    INTEGER                :: subcell        ! current subcell 
    INTEGER                :: subIndex       ! octal's child index for current subcell 
    INTEGER                :: i


    ! find which of the subcells the point lies in
    DO 
      subcell = whichSubcell(octree,locator)

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
                    chiLine=chiLine,levelPop=levelPop,rho=rho,temperature=temperature)
        
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

!        write(*,*) "centre",centre
!        write(*,*) "subcellsize",subcellsize
!        write(*,*) "current",currentpoint

        ! we find the exit point from the current subcell, 
        !  and also 'locator' - a point that lies in the *next* subcell 
        CALL getExitPoint(currentPoint,direction,locator,abortRay,error, &
                          grid%halfSmallestSubcell,exitPoint,centre,subcellSize,&
                          minWallDistance,margin,grid%octreeRoot%threed)
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
!                          levelPop=levelPop,rho=rho)
!        ENDIF

        
        ! we check whether we should take a sample at currentPoint
        !   (which will usually be the entry point to the subcell).

        length = modulus(currentPoint - startPoint)

! force sampling here        
        !IF ( modulus(currentPoint - (startPoint + (direction *         & 
        !                      REAL(lambda(nSamples),KIND=octalKind)))) &
        !                                       > sampleLength ) THEN
          CALL takeSample(currentPoint,length,direction,grid,octree,subcell,    &
                          nSamples,maxSamples,usePops,iLambda,error,lambda,     &
                          kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,&
                          velocityDeriv=velocityDeriv,chiLine=chiLine,          &
                          levelPop=levelPop,rho=rho,temperature=temperature)
        !END IF

        ! we add sampleLength to the distance from the last location
        !   that was sampled, and decide whether to take a new sample.
        trial = 1
        DO 
          trialPoint = currentPoint + direction * &
                       ((REAL(lambda(nSamples),KIND=octalKind) + &
                       (REAL(trial,KIND=octalKind) * sampleLength)) - length)
                       
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
                         levelPop=levelPop,rho=rho,temperature=temperature)
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
    REAL(KIND=octalKind), INTENT(IN)    :: halfSmallestSubcell
    REAL(KIND=octalKind), INTENT(IN)    :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    TYPE(octalVector), INTENT(OUT)      :: exitPoint    ! where ray leaves current cell
    TYPE(octalVector), INTENT(IN)       :: centre       ! centre of current subcell
    REAL(octalKind), INTENT(IN)         :: subcellSize  ! size of current subcell
    REAL(KIND=octalKind), INTENT(OUT)   :: minWallDistance ! distance to *nearest* wall
    
    
    REAL(KIND=octalKind)   :: wallDistanceX  ! distance to next x-wall-plane 
    REAL(KIND=octalKind)   :: wallDistanceY  ! distance to next y-wall-plane
    REAL(KIND=octalKind)   :: wallDistanceZ  ! distance to next z-wall-plane
                  ! 'wallDistanceX,Y,Z' are the distances between the current 
                  !   position and the intersections with the cell walls (along the 
                  !   'direction' vector)
    LOGICAL                :: found          ! status flag  
    REAL(KIND=octalKind)   :: wallFromOrigin ! distance of cell wall from the origin
    TYPE(octalVector)      :: wallNormal     ! normal to plane of cell wall

    INTEGER, parameter :: max_num_err = 10;
    INTEGER, save :: num_err = 0    

    LOGICAL :: threed
    REAL(KIND=octalKind) :: r2, d, cosMu, disttor2, disttoxboundary,disttozboundary
    REAL(KIND=octalKind) :: compZ, currentZ, tval, x1, x2, currentx, compx
    REAL(KIND=octalKind) :: r1, theta, mu, disttor1

    TYPE(octalvector) :: xDir, zDir
    logical :: ok
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

       r1 = centre%x - subcellSize/2.
       r2 = centre%x + subcellSize/2.
       d = sqrt(currentpoint%x**2+currentpoint%y**2)
       xDir = VECTOR(currentpoint%x, currentpoint%y,0.d0)

       if (modulus(xDir) /= 0.d0) then
          call normalize(xDir)       
          cosmu =((-1.d0)*xDir).dot.direction
          call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
          if (.not.ok) then
             write(*,*) "Quad solver failed in intersectcubeamr2d"
             do;enddo
          endif
          distTor2 = max(x1,x2)      
          
          
          theta = asin(max(-1.d0,min(1.d0,r1 / d)))
          cosmu = xDir.dot.direction
          mu = acos(max(-1.d0,min(1.d0,cosmu)))
          distTor1 = 1.e30
          if (mu  < theta ) then
             call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
             if (.not.ok) then
                write(*,*) "Quad solver failed in intersectcubeamr2d"
                do;enddo
             endif
          endif
          distTor1 = max(x1,x2)
          distToXboundary = min(distTor1, distTor2)
       else
          distToXboundary = 1.e30
       endif

       zDir = VECTOR(0.d0, 0.d0, 1.d0)
       compZ = zDir.dot.direction
       currentZ = currentpoint%z

       if (compZ /= 0.d0 ) then
          if (compZ > 0.d0) then
             distToZboundary = (centre%z + subcellsize/2. - currentZ ) / compZ
          else
             distToZboundary = abs((centre%z - subcellsize/2. - currentZ ) / compZ)
          endif
       else
          disttoZboundary = 1.e30
       endif

       tVal = min(distToZboundary, distToXboundary)
       if (tVal > 1.e29) then
          write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
          write(*,*) "x,z",currentX,currentZ
          do;enddo
       endif
       if (tval < 0.) then
          write(*,*) tVal,compX,compZ, distToZboundary,disttoxboundary
          write(*,*) "x,z",currentX,currentZ
          do;enddo
       endif
       
       minWallDistance = tVal
       
       exitPoint = currentPoint + tval * direction
       locator = exitPoint + halfSmallestSubcell * direction

    endif





  END SUBROUTINE getExitPoint 


  SUBROUTINE takeSample(point,length,direction,grid,thisOctal,subcell,nSamples,&
                        maxSamples,usePops,iLambda,error,lambda,kappaAbs,      &
                        kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,temperature) 
  
    
    IMPLICIT NONE
    
    TYPE(octalVector), INTENT(IN)      :: point     ! place to make sample
    REAL(KIND=octalKind), INTENT(IN)   :: length    ! 
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
    
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: kappaAbs  ! continuous absorption opacities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: kappaSca  ! scattering opacities 
    TYPE(vector), DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: velocityDeriv   ! sampled velocity derivatives
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: chiLine   ! line opacities
    REAL,DIMENSION(:,:),INTENT(INOUT),OPTIONAL        :: levelPop  ! level populations
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: rho       ! density at sample points
    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL          :: temperature       ! temp at sample points

    TYPE(vector)                       :: directionReal ! direction vector (REAL values)
    TYPE(octal), POINTER               :: localPointer  ! pointer to the current octal
 
    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    REAL :: Vr, r, Rs
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn)

    nSamples = nSamples + 1
    IF (nSamples > maxSamples) THEN
      PRINT *, "nSamples > maxSamples in takeSample subroutine"
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
                         grid=grid, temperature=temperature(nSamples))
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
                         grid=grid,temperature=temperature(nSamples))
    END IF

    
    ! special case if the geometry is "jets" (for accuracy)
    If (grid%geometry(1:4) == "jets" ) then
       ! using routines in jets_mod.f90
       rho(nSamples) = Density(pointVec, grid)                 ! in [g/cm^3]
       Vr = JetsVelocity(pointVec, grid)/(cSpeed/1.0e5)        ! in [c]
       r = modulus(pointVec)
       Rs = get_jets_parameter("Rmin")
       if (r < Rs) r =Rs
       velocity(nSamples) = REAL(Vr, kind=octalkind)*(pointVec/REAL(r,kind=octalkind))
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
                           probDistLine,probDistCont,N,Ne,nTot,inflow,grid,interp)
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
    REAL,INTENT(OUT),OPTIONAL         :: kappaAbs
    REAL,INTENT(OUT),OPTIONAL         :: kappaSca
    REAL,INTENT(OUT),OPTIONAL         :: rho
    REAL,INTENT(OUT),OPTIONAL         :: chiLine
    REAL,INTENT(OUT),OPTIONAL         :: etaLine
    REAL,INTENT(OUT),OPTIONAL         :: etaCont
    REAL,INTENT(OUT),OPTIONAL         :: probDistLine
    REAL,INTENT(OUT),OPTIONAL         :: probDistCont
    REAL,DIMENSION(:),INTENT(OUT),OPTIONAL :: N
    REAL,INTENT(OUT),OPTIONAL         :: Ne
    REAL,INTENT(OUT),OPTIONAL         :: nTot
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
      IF (PRESENT(velocity))         velocity = resultOctal%velocity(subcell)
      IF (PRESENT(temperature))   temperature = resultOctal%temperature(subcell)
      IF (PRESENT(rho))                   rho = resultOctal%rho(subcell)
      IF (PRESENT(chiLine))           chiLine = resultOctal%chiLine(subcell)
      IF (PRESENT(etaLine))           etaLine = resultOctal%etaLine(subcell)
      IF (PRESENT(etaCont))           etaCont = resultOctal%etaCont(subcell)
      IF (PRESENT(probDistLine)) probDistLine = resultOctal%probDistLine(subcell)
      IF (PRESENT(probDistCont)) probDistCont = resultOctal%probDistCont(subcell)
      IF (PRESENT(Ne))                     Ne = resultOctal%Ne(subcell)
      IF (PRESENT(N))                       N = resultOctal%N(subcell,:)
      IF (PRESENT(nTot))                 nTot = resultOctal%nTot(subcell)
      IF (PRESENT(inFlow))           inFlow = resultOctal%inFlow(subcell)     

      IF (PRESENT(kappaAbs)) THEN 
        IF (PRESENT(iLambda)) THEN
           if (.not.grid%oneKappa) then
              kappaAbs = resultOctal%kappaAbs(subcell,iLambda)
           else
              IF (.NOT.PRESENT(lambda)) THEN
                 kappaAbs = grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell)
              ELSE
                 kappaAbs = logint(lambda, grid%lamArray(ilambda), grid%lamArray(ilambda+1), &
                 grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell), &
                 grid%oneKappaAbs(resultOctal%dustType(subcell),iLambda+1)*resultOctal%rho(subcell))
              ENDIF
           endif
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''kappaAbs'' without',&
                   ' ''iLambda''.'
          do;enddo
        END IF
      END IF
      
      IF (PRESENT(kappaSca)) THEN 
        IF (PRESENT(iLambda)) THEN
           if (.not.grid%oneKappa) then
              kappaSca = resultOctal%kappaSca(subcell,iLambda)
           else
              IF (.NOT.PRESENT(lambda)) THEN
                 kappaSca = grid%oneKappaSca(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell)
              ELSE
                 kappaSca = logint(lambda, grid%lamArray(ilambda), grid%lamArray(ilambda+1), &
                 grid%oneKappaSca(resultOctal%dustType(subcell),iLambda)*resultOctal%rho(subcell), &
                     grid%oneKappaSca(resultOctal%dustType(subcell),iLambda+1)*resultOctal%rho(subcell))
              ENDIF
           endif
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''kappaSca'' without',&
                   ' ''iLambda''.'
          STOP
        END IF
      END IF

      IF (PRESENT(velocityDeriv)) THEN
        IF (PRESENT(direction) .AND. PRESENT(grid)) THEN
          velocityDeriv = amrGridDirectionalDeriv(grid,point2,direction,&
                                                  startOctal=resultOctal)
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
    type(vector) :: newVec
    REAL(KIND=octalKind)           :: inc
    real :: phi
    REAL(KIND=octalKind)           :: t1, t2, t3

    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      resultOctal => startOctal
      
    ELSE
      CALL findSubcellTD(point,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal))   foundOctal   => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

    END IF

      inc = resultOctal%subcellSize / 2.0
      centre = subcellCentre(resultOctal,subcell)
      
      t1 = MAX(0.0_oc, (point%x - (centre%x - inc)) / resultOctal%subcellSize)
      t2 = MAX(0.0_oc, (point%y - (centre%y - inc)) / resultOctal%subcellSize)
      t3 = MAX(0.0_oc, (point%z - (centre%z - inc)) / resultOctal%subcellSize)
    

!      amrGridVelocity = vector(0.,0.,0.)  ! TJH did this 
!      amrGridVelocity = resultOctal%velocity(subcell) ! RK did this.

!       if (.false.) then

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
            
      END SELECT
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

    IMPLICIT NONE

    TYPE(gridtype), INTENT(IN)     :: grid 
    REAL                           :: amrGridDirectionalDeriv
    TYPE(octalVector), INTENT(IN)  :: position
    TYPE(vector), INTENT(IN)       :: direction
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL :: foundSubcell

    TYPE(octalVector)              :: octalDirection
    TYPE(octal), POINTER           :: firstOctal
    REAL(KIND=octalKind)           :: dr, dx, dphi
    REAL(KIND=octalKind)           :: r
    REAL(KIND=octalKind)           :: phi1, phi2
    TYPE(octalVector)              :: position1
    TYPE(octalVector)              :: position2
    INTEGER                        :: subcell
    
    octalDirection = direction
    
    ! dr is a small increment of distance
    dr = grid%halfSmallestSubcell * 2.0 

    ! get a new position a little way back from current position
    position1 = position - (dr * octalDirection)
    
    ! this might be inside core or outside grid - in which case
    !   just use the current position as the first point

    r = modulus(position1)
    IF (.NOT. inOctal(grid%octreeRoot,position1) .OR. (r < grid%rCore)) THEN
      position1 = position
    END IF
      
    ! first line of sight velocity
    phi1 = direction .dot. AMRgridVelocity(grid%octreeRoot,position1,&
                                        startOctal=startOctal,&
                                        foundOctal=firstOctal)
                                        
    
    ! now go forward a bit from current position
    position2 = position + (dr * octalDirection)

    ! check we're still inside grid
    r = modulus(position2)
    IF (.NOT. inOctal(grid%octreeRoot,position2) .OR. (r < grid%rCore)) THEN
      position2 = position
    END IF

    ! the second position l.o.s. velocity
    phi2 = direction .dot. AMRgridVelocity(grid%octreeRoot,position2,&
                                        startOctal=firstOctal,&
                                        foundOctal=foundOctal,&
                                        foundSubcell=subcell)

    IF (PRESENT(foundSubcell)) foundSubcell = subcell                                        

    dx = modulus(position2 - position1)

    dphi = phi2 - phi1

    ! the line of sight velocity gradient

    IF (ABS(dx) >= 1.e-10) THEN
       amrGridDirectionalDeriv = (abs(dphi / dx))
    ELSE
       amrGridDirectionalDeriv = 1.e-10
    ENDIF
    
    IF (.NOT. ABS(amrGridDirectionalDeriv) >= 1.e-10) &
       amrGridDirectionalDeriv = 1.e-10

  END FUNCTION amrGridDirectionalDeriv

    
  FUNCTION amrGridKappaAbs(octalTree,point,iLambda,startOctal,&
                           foundOctal,foundSubcell,actualSubcell) 
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
    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      
      ! should add in some interpolation routine
      amrGridKappaAbs = startOctal%kappaAbs(subcell,iLambda)
      
    ELSE
      CALL findSubcellTD(point,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal))   foundOctal   => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

      ! should add in some interpolation routine
      amrGridKappaAbs = resultOctal%kappaAbs(subcell,iLambda)

    END IF
  
  END FUNCTION amrGridKappaAbs

  FUNCTION amrGridKappaSca(octalTree,point,iLambda,startOctal,&
                           foundOctal,foundSubcell,actualSubcell) 
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
    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      
      ! should add in some interpolation routine
      amrGridKappaSca = startOctal%kappaSca(subcell,iLambda)
      
    ELSE
      CALL findSubcellTD(point,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal))   foundOctal   => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

      ! should add in some interpolation routine
      amrGridKappaSca = resultOctal%kappaSca(subcell,iLambda)

    END IF
  
  END FUNCTION amrGridKappaSca

  FUNCTION amrGridTemperature(octalTree,point,startOctal,foundOctal,& 
                                        foundSubcell,actualSubcell) 
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
    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      
      ! should add in some interpolation routine
      amrGridTemperature = startOctal%temperature(subcell)
      
    ELSE
      CALL findSubcellTD(point,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal)) foundOctal => resultOctal
      IF (PRESENT(foundSubcell)) foundSubcell =  subcell

      ! should add in some interpolation routine
      amrGridTemperature = resultOctal%temperature(subcell)

    END IF

  END FUNCTION amrGridTemperature

  FUNCTION amrGridDensity(octalTree,point,startOctal,foundOctal,& 
                                        foundSubcell,actualSubcell) 
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
    
    IF (PRESENT(startOctal)) THEN
      IF (PRESENT(actualSubcell)) THEN
        subcell = actualSubcell
      ELSE 
        CALL findSubcellLocal(point,startOctal,subcell)
        IF (PRESENT(foundOctal))   foundOctal   => startOctal
        IF (PRESENT(foundSubcell)) foundSubcell =  subcell
      END IF
      
      ! should add in some interpolation routine
      amrGridDensity = startOctal%rho(subcell)
      
    ELSE
      CALL findSubcellTD(point,octalTree,resultOctal,subcell)
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

    REAL(KIND=doubleKind), INTENT(IN) :: probability
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

    REAL(KIND=doubleKind), INTENT(IN) :: probability
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
    ! returns the identification number (1-8) of the subcell of the 
    ! current octal which contains a given point
    ! NB this does NOT check that the point lies within the bounds of the octal!
  
    IMPLICIT NONE
    
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(octalVector)             :: rotpoint
    INTEGER                       :: subcell

    if (thisOctal%threed) then
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
       ENDIF
    else ! twoD case
       rotpoint = projecttoxz(point)
       IF ( rotpoint%x <= thisOctal%centre%x ) THEN
          IF ( rotpoint%z <= thisOctal%centre%z ) THEN
             subcell = 1
          ELSE 
             subcell = 3
          ENDIF
       ELSE
          IF (rotpoint%z <= thisOctal%centre%z) THEN
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


    

    if (thisOctal%threeD) then
           IF (point%x <= thisOctal%centre%x - thisOctal%subcellSize ) THEN ; inOctal = .FALSE. 
       ELSEIF (point%x >= thisOctal%centre%x + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSEIF (point%y <= thisOctal%centre%y - thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSEIF (point%y >= thisOctal%centre%y + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSEIF (point%z <= thisOctal%centre%z - thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSEIF (point%z >= thisOctal%centre%z + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSE  
          inOctal = .TRUE.
       ENDIF
    else ! twoD case
       octVec2D = projectToXZ(point)
           IF (octVec2D%x <  thisOctal%centre%x - thisOctal%subcellSize ) THEN ; inOctal = .FALSE. 
       ELSEIF (octVec2D%x >= thisOctal%centre%x + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSEIF (octVec2D%z <= thisOctal%centre%z - thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSEIF (octVec2D%z >= thisOctal%centre%z + thisOctal%subcellSize ) THEN ; inOctal = .FALSE.
       ELSE  
          inOctal = .TRUE.
       ENDIF
    endif
  END FUNCTION inOctal


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


  RECURSIVE SUBROUTINE smoothAMRgrid(thisOctal,grid,factor,gridConverged, sphData, stellar_cluster)
    ! checks whether each octal's neighbours are much bigger than it, 
    !   if so, makes the neighbours smaller.
    ! the 'needRestart' flag will be set if a change is made to the octree.
    !   smoothAMRgrid should be called repeatedly until the flag is no 
    !   longer set - this should ensure that the octree eventually stabilizes
  
    IMPLICIT NONE

    TYPE(octal), POINTER             :: thisOctal
    TYPE(gridtype), INTENT(INOUT   ) :: grid 
     REAL, INTENT(IN)                 :: factor
    LOGICAL, INTENT(INOUT)               :: gridConverged
    TYPE(sph_data), optional, INTENT(IN) :: sphData   ! Matthew's SPH data.
    TYPE(cluster), optional, intent(in)  :: stellar_cluster

    INTEGER              :: i
    REAL(KIND=octalKind) :: halfSmallestSubcell
    TYPE(octal), POINTER :: child
    TYPE(octal), POINTER :: neighbour
    TYPE(octalVector), ALLOCATABLE, DIMENSION(:) :: locator
    INTEGER              :: subcell
    INTEGER              :: nLocator ! number of locators (4 for twoD, 6 for threed)
    
    gridConverged = .TRUE.

    ! we will find the coordinates of a point that lies outside the current
    !   octal. we then compare the size of the cell that contains that point
    !   with the size of the current cell, if it is bigger be more than a 
    !   factor of 'factor', we subdivide the neighbouring cell.
    ! we do this in each of six directions

    ! we do not have to test the other subcells in the current octal because
    !   they can be smaller than the any of the other subcells, but they
    !   cannot be *bigger*. this saves some time.


    if (thisOctal%threed) then
       nlocator = 6
    else
       nlocator = 4
    endif

    ALLOCATE(locator(nlocator))

    FORALL (i = 1:nLocator)
      locator(i) = thisOctal%centre
    END FORALL

    ! we find points which are outside the current octal by a distance
    !   equivalent to half the size of the tree's smallest subcell.
    halfSmallestSubcell = grid%halfSmallestSubcell

    if (thisOctal%threed) then
       locator(1)%x = thisOctal%centre%x + thisOctal%subcellSize + halfSmallestSubcell
       locator(2)%y = thisOctal%centre%y + thisOctal%subcellSize + halfSmallestSubcell
       locator(3)%z = thisOctal%centre%z + thisOctal%subcellSize + halfSmallestSubcell
       locator(4)%x = thisOctal%centre%x - thisOctal%subcellSize - halfSmallestSubcell
       locator(5)%y = thisOctal%centre%y - thisOctal%subcellSize - halfSmallestSubcell
       locator(6)%z = thisOctal%centre%z - thisOctal%subcellSize - halfSmallestSubcell
    else
       locator(1)%x = thisOctal%centre%x + thisOctal%subcellSize + halfSmallestSubcell
       locator(2)%z = thisOctal%centre%z + thisOctal%subcellSize + halfSmallestSubcell
       locator(3)%x = thisOctal%centre%x - thisOctal%subcellSize - halfSmallestSubcell
       locator(4)%z = thisOctal%centre%z - thisOctal%subcellSize - halfSmallestSubcell
    endif
    DO i = 1, nLocator, 1
      IF ( inOctal(grid%octreeRoot,locator(i)) ) THEN
        CALL findSubcellTD(locator(i),grid%octreeRoot,neighbour,subcell)
        neighbour => thisOctal
        CALL findSubcellLocal(locator(i),neighbour,subcell)

        IF ( neighbour%subcellSize > (factor * thisOctal%subcellSize) ) THEN
          IF ( neighbour%hasChild(subcell) ) THEN 
            PRINT *, "neighbour already has child"
            STOP
          END IF
          ! tjh changed from
!          CALL addNewChild(neighbour,subcell,grid)
          call addNewChildren(neighbour, grid, sphData, stellar_cluster)
          gridConverged = .FALSE.
        ENDIF
      END IF
    END DO

    ! call this subroutine recursively on each of any children.
    IF ( thisOctal%nChildren > 0 ) THEN
       DO i = 1, thisOctal%nChildren, 1 
          child => thisOctal%child(i)
          CALL smoothAMRgrid(child,grid,factor,gridConverged, sphData, stellar_cluster)
       END DO
    END IF

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


  SUBROUTINE findSubcellLocal(point,thisOctal,subcell)
    ! finds the octal (and that octal's subcell) containing a point.
    !   starts searching from the current octal, and goes up and down the
    !   tree as needed to find the correct octal.

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(octal),POINTER    :: thisOctal
    INTEGER, INTENT(OUT)   :: subcell
    
    LOGICAL                :: haveDescended    ! see comments below
    LOGICAL                :: boundaryProblem  ! see comments below
    
                             
    haveDescended = .FALSE.   ! if the 'point' lies very close to an 
    boundaryProblem = .FALSE. !   boundary, the program may go into 
                              !   a loop going up and down the tree.
                              ! we will keep track of the progress of
                              !   the search using these flags.
                              
    CALL findSubcellLocalPrivate(point,thisOctal,subcell,&
                                 haveDescended,boundaryProblem)
                                 
  CONTAINS

    RECURSIVE SUBROUTINE findSubcellLocalPrivate(point,thisOctal,subcell,&
                                                 haveDescended,boundaryProblem)
      TYPE(octalVector), INTENT(IN) :: point
      TYPE(octal),POINTER    :: thisOctal
      INTEGER, INTENT(OUT)   :: subcell
      LOGICAL, INTENT(INOUT) :: haveDescended
      LOGICAL, INTENT(INOUT) :: boundaryProblem
      
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
!          STOP
!           DO ; END DO
          boundaryProblem = .TRUE.
          RETURN
        END IF
     
        ! if we have previously gone down the tree, and are now going back up, there
        !   must be a problem.
        IF (haveDescended) boundaryProblem = .TRUE.
        
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
    real(kind=doubleKind) :: totalMass, minRho, maxRho
    integer maxPoints
    integer :: nPoints
    type(VECTOR),allocatable :: points(:), lineGrid(:,:)
    real :: x, y, dr
    integer,allocatable :: depth(:)
    integer :: nGrid
    integer, allocatable :: nColour(:)
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
    real(kind=doubleKind) :: totalMass
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
    real(KIND=octalKind) :: cubeSizeOctal
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



  FUNCTION decideSplit(thisOctal,subcell,amrLimitScalar,amrLimitScalar2,grid, &
       sphData, stellar_cluster) RESULT(split)
    ! returns true if the current voxel is to be subdivided. 
    ! decision is made by comparing 'amrLimitScalar' to some value
    !   derived from information in the current cell  

    use input_variables, only: height
    IMPLICIT NONE
    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    REAL(KIND=doubleKind), INTENT(IN) :: amrLimitScalar, amrLimitScalar2 ! used for split decision
    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(sph_data), OPTIONAL, intent(in) :: sphData
    TYPE(cluster), OPTIONAL, intent(in) :: stellar_cluster
    LOGICAL                    :: split          
    REAL(KIND=octalKind)  :: cellSize
    TYPE(octalVector)     :: searchPoint, octVec
    TYPE(octalVector)     :: cellCentre
    REAL                  :: x, y, z
    REAL :: hr, rd, fac, warpHeight, phi
    INTEGER               :: i
    DOUBLE PRECISION      :: total_mass
    DOUBLE PRECISION      :: ave_density, rGrid(1000), r
    INTEGER               :: nr, nr1, nr2
    DOUBLE PRECISION      :: total_opacity, minDensity, maxDensity, thisDensity, rho
    INTEGER               :: nsample = 400
    LOGICAL               :: inUse
    INTEGER               :: nparticle
    REAL(KIND=doubleKind) :: dummyDouble
    REAL(KIND=doubleKind) :: rho_disc_ave, scale_length


   select case(grid%geometry)

    case("ttauri","jets","spiralwind","melvin")
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
           searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=octalKind) 
           searchPoint%y = searchPoint%y - (cellSize / 2.0_oc) + cellSize*REAL(y,KIND=octalKind) 
           searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=octalKind) 
        else
           searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=octalKind) 
           searchPoint%y = 1.e-30
           searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=octalKind) 
        endif
        ! using a generic density function in density_mod.f90
        rho = density(searchPoint,grid)
        ave_density  =  rho + ave_density
        minDensity = min(minDensity, rho)
        maxDensity = max(maxDensity, rho)
      END DO

      ave_density = ave_density / REAL(nsample,KIND=doubleKind)
      total_mass = maxDensity * (cellSize*1.e10_db)**3.0_db
      IF (total_mass > amrLimitScalar) then
         split = .TRUE.
      END IF

   case ("testamr","proto")
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      split = .FALSE.
      nr1 = 5
      nr2 = 30
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
         rgrid(nr1+i) = log10(1.2*grid%rInner)+dble(i)*(log10(grid%rOuter)-log10(1.2*grid%rInner))/dble(nr2)
      end do
      rgrid(1:nr) = 10.d0**rgrid(1:nr)
      r = modulus(cellcentre)
      if (thisOctal%nDepth < 6) split = .true.
      if ((r < grid%rOuter).and.(r > grid%rinner)) then
         call locate(rGrid, nr, r, i)      
         if (cellsize > (rGrid(i+1)-rGrid(i))) split = .true.
      endif

   case("benchmark")
      split = .false.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
      if (thisOctal%nDepth < 6) split = .true.
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      if (grid%geometry == "benchmark") then
         if (r > grid%rInner) then
            rd = grid%rOuter / 2. 
            hr = height * (r/rd)**1.125
            if (abs(cellCentre%z)/hr < 10.) then
               if (cellsize/hr > 1.) split = .true.
            endif
            if (sqrt(cellcentre%x**2 + cellcentre%y**2)/grid%rInner < 4.) then
               if (abs(cellCentre%z)/hr < 10.) then
                  if (cellsize/hr > 0.4) split = .true.
               endif
            endif
         endif
      endif

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
      if (nparticle> 6) then
         split = .true.
!      else
!         split = .false.
      end if
      

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
      endif

   case("shakara","aksco")
      split = .false.
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
!      if ((thisOctal%nDepth < 8).and.(cellCentre%x < grid%rOuter)) split = .true.
      r = sqrt(cellcentre%x**2 + cellcentre%y**2)
      phi = atan2(cellcentre%y,cellcentre%x)
      warpHeight = 0. !cos(phi) * grid%rInner * sin(30.*degtorad) * sqrt(grid%rinner / r)

      
      if (r > grid%rInner*0.8) then
         hr = height * (r / (100.*autocm/1.e10))**1.25
         fac = cellsize/hr
         if (abs((cellCentre%z-warpHeight)/hr) < 10.) then
            if (fac > 2.) split = .true.
         endif
      endif
!      if ((r < 100.*rSol/1.e10).and.(abs(cellCentre%z)< 20.*rSol/1.e10)) then
!         if (thisOctal%subcellSize > 0.5*rSol/1.e10) split = .true.
!      endif

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

     
   case ("wr104")
      ! Splits if the number of particle is more than a critical value (~3).
      
      call find_n_particle_in_subcell(nparticle, dummyDouble, sphData, thisOctal, subcell)

      !
      if (nparticle > nint(amrLimitScalar) ) then 
	 split = .TRUE.
      else
         split = .FALSE.
      end if



      
   case DEFAULT
      PRINT *, 'Invalid grid geometry option passed to amr_mod::decideSplit'
      PRINT *, 'grid%geometry ==', TRIM(grid%geometry)
      PRINT *, 'Exiting the program .... '
      STOP
   end select
   

  END FUNCTION decideSplit


  FUNCTION columnDensity(grid,startPoint,direction,sampleFreq) RESULT(rho)
                 
    ! integrates the density along a line through an amr grid.

    IMPLICIT NONE

    TYPE(gridtype), INTENT(IN)        :: grid
    TYPE(octalVector), INTENT(IN)     :: startPoint
    TYPE(octalVector), INTENT(IN)     :: direction
    REAL, INTENT(IN)  :: sampleFreq

    REAL                              :: rho

    INTEGER, PARAMETER                :: maxSamples = 10000
    REAL, DIMENSION(maxSamples)       :: distances, densities
    REAL, DIMENSION(maxSamples)       :: dummy
    REAL, DIMENSION(maxSamples,1)     :: dummyPops
    TYPE(vector),DIMENSION(maxSamples):: dummyVel
    LOGICAL                           :: hitCore
    
    INTEGER :: nSamples, iSample

    INTEGER :: error

    error = 0 

    nSamples = 0
                
    CALL startReturnSamples(startPoint,direction,grid,sampleFreq,     &
                 nSamples,maxSamples,.false.,hitCore,.false.,1,error, &
                 distances,kappaAbs=dummy,kappaSca=dummy,velocity=dummyVel,&
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

    REAL(octalKind)      :: x1, x2, x3
    REAL(octalKind)      :: y1, y2, y3
    REAL(octalKind)      :: z1, z2, z3
    
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

    USE parameters_mod

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    TYPE(octalvector)      :: vP
    REAL(kind=doubleKind)  :: modVp
    REAL(kind=doubleKind)  :: phi
    REAL(kind=doubleKind)  :: r, rM, theta, y

    starPosn = grid%starPos1
    pointVec = (point - starPosn) * 1.e10_oc
    r = modulus( pointVec ) 
    if (r /= 0.d0) then
       theta = ACOS( max(-1.d0,min(1.d0,pointVec%z / dble(r) )))
    else
       theta = 90.d0
    endif
    if (theta == 0.d0) theta=1.d-20
    rM  = r / SIN(theta)**2
    y = SIN(theta)**2 

    ! test if the point lies within the star
    IF ( r < TTauriRstar ) THEN
      TTauriVelocity = vector(1.e-25,1.e-25,1.e-25)
      
    ! test if the point lies too close to the disk
    ELSE IF ( ABS(pointVec%z) < 4.0e-2 * TTauriRstar) THEN
      TTauriVelocity = vector(1.e-25,1.e-25,1.e-25)
   
    ! test if the point lies outside the accretion stream
    ELSE IF ((rM > TTauriRinner) .AND. (rM < TTauriRouter )) THEN
  
      vP = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
                  0.0, &
                 (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))
      modVp = SQRT((2.0 * bigG * TTauriMstar / TTauriRstar) * &
                     (TTauriRstar/r - TTauriRstar/rM))
      vP = (-1.0 * (modVp/cSpeed)) * vP
      phi = ATAN2(pointVec%y,pointVec%x)
      vP = rotateZ(vP,phi) 
      IF (theta > pi/2.0) vP%z = -vP%z
      TTauriVelocity = vP

    ELSE
      TTauriVelocity = vector(1.e-25,1.e-25,1.e-25)
    END IF

  END FUNCTION TTauriVelocity


  SUBROUTINE calcTTauriMassVelocity(thisOctal,subcell,grid) 
    ! calculates some of the variables at a given point for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    USE parameters_mod

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
    REAL :: r, rM, theta, y

    starPosn = grid%starPos1

      point = subcellCentre(thisOctal,subcell)
      pointVec = (point - starPosn) * 1.e10_oc
      r = modulus( pointVec ) 
      theta = ACOS( pointVec%z / r )
      rM  = r / SIN(theta)**2
      y = SIN(theta)**2 

      ! test if the point lies within the star
      IF ( r < TTauriRstar ) THEN
        thisOctal%rho(subcell) = 1.e-25
        thisOctal%velocity(subcell) = vector(1.e-25,1.e-25,1.e-25)
        thisOctal%inFlow(subcell) = .FALSE.
      
      ! test if the point lies too close to the disk
      ELSE IF ( ABS(pointVec%z) < 4.0e-2 * TTauriRstar) THEN
        thisOctal%rho(subcell) = 1.e-25
        thisOctal%velocity(subcell) = vector(1.e-25,1.e-25,1.e-25)
        thisOctal%inFlow(subcell) = .FALSE.
   
      ! test if the point lies outside the accretion stream
      ELSE IF ((rM > TTauriRinner) .AND. (rM < TTauriRouter )) THEN
        thisOctal%inFlow(subcell) = .TRUE.
  
        rho = (TTauriMdot * TTauriRstar) / (4.0 * pi * &
                 (TTauriRstar/TTauriRinner - TTauriRstar/TTauriRouter))
        rho = rho * r**(-5.0/2.0) / SQRT( 2.0 * bigG * TTauriMstar ) 
        rho = rho * SQRT( 4.0 - 3.0*y) / SQRT( 1.0 - y) 
        TTauriMinRho = MIN(TTauriMinRho, rho/mHydrogen)
        TTauriMaxRho = MAX(TTauriMaxRho, rho/mHydrogen)
        thisOctal%rho(subcell) = rho

        vP = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
                    0.0, &
                   (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))
        modVp = SQRT((2.0 * bigG * TTauriMstar / TTauriRstar) * &
                       (TTauriRstar/r - TTauriRstar/rM))
        vP = (-1.0 * (modVp/cSpeed)) * vP
        phi = ATAN2(pointVec%y,pointVec%x)
        vP = rotateZ(vP,phi) 
        IF (theta > pi/2.0) vP%z = -vP%z
        thisOctal%velocity(subcell) = vP

      ELSE
        thisOctal%inFlow(subcell) = .FALSE.
        thisOctal%rho(subcell) = 1.e-25
        thisOctal%velocity(subcell) = vector(1.e-25,1.e-25,1.e-25)
      END IF
      

      IF ((thisoctal%threed).and.(subcell == 8)) &
        CALL fillVelocityCorners(thisOctal,grid,TTauriVelocity, .true.)

      IF ((thisoctal%twod).and.(subcell == 4)) &
        CALL fillVelocityCorners(thisOctal,grid,TTauriVelocity, .false.)


  END SUBROUTINE calcTTauriMassVelocity

  
  PURE SUBROUTINE calcTTauriTemperature(thisOctal,subcell) 
    ! calculates the temperature in an octal for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    USE parameters_mod

    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell

    REAL :: rho

    IF ( thisOctal%inFlow(subcell) ) THEN
      rho = thisOctal%rho(subcell)
      thisOctal%temperature(subcell) = MAX(5000.0, &
        7000.0 - ((2500.0 * rho/mHydrogen - TTauriMinRho) / (TTauriMaxRho-TTauriMinRho)))
    ELSE
      thisOctal%temperature(subcell) = 6500.0
    END IF

    ! we will initialise the bias distribution
    thisOctal%biasLine3D(subcell) = 1000.0
  
  END SUBROUTINE calcTTauriTemperature


  subroutine initTTauriAMR(grid,Laccretion,Taccretion,&
                           sAccretion,newContFile,theta1,theta2)

    use constants_mod
    use vector_mod
    use input_variables
    use parameters_mod

    implicit none
    
    type(GRIDTYPE), intent(inout)      :: grid                
    real(kind=doublekind), intent(out) :: Laccretion 
    real, intent(out)                  :: Taccretion 
    real, intent(out)                  :: sAccretion
    character(len=80), intent(out)     :: newContFile
    real, intent(out) :: theta1, theta2
   
    real(kind=doublekind) :: TaccretionDouble
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
    grid%diskRadius = TTauriRInner * 1.e-10
    grid%diskNormal = VECTOR(0.,0.,1.)
    grid%diskNormal = rotateX(grid%diskNormal,grid%dipoleOffSet)
    grid%starPos1 = VECTOR(0.,0.,0.) ! in units of 1.e-10 cm
    grid%starPos2 = VECTOR(9.e9,9.e9,9.e9) ! in units of 1.e-10 cm


    theta1 = asin(sqrt(TTauriRstar/TTauriRouter))
    theta2 = asin(sqrt(TTauriRstar/TTauriRinner))


    Laccretion = (REAL(bigG,KIND=db)* &
                  REAL(TTauriMstar,KIND=db)* &
                  REAL(TTauriMdot,KIND=db)/ &
                  REAL(TTauriRstar,kind=db))* &
            REAL((1.0_db-(2.0_db*TTauriRstar/(TTauriRouter+TTauriRinner))),KIND=db)

    TaccretionDouble = Laccretion / REAL(((fourPi * TTauriRstar**2)*stefanBoltz* &
                                      abs(cos(theta1)-cos(theta2))),kind=db)

    sAccretion = (fourPi * TTauriRstar**2)*abs(cos(theta1)-cos(theta2))!1.e20
    Taccretion = TaccretionDouble**0.25

    write(*,*) "accretion lum/temp",Laccretion/Lsol, Taccretion

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
    
    open(22,file="star_plus_acc.dat",form="formatted",status="unknown")
    do i = 1, nNu
       fNu(i) = fNu(i) + blackbody(tAccretion, 1.e8*cSpeed/ nuArray(i))* &
                (sAccretion/(fourPi*TTauriRstar**2))
       write(22,*) nuArray(i), fNu(i)
    enddo
    close(22)
    newContfile="star_plus_acc.dat"


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

    thisOctal%rho(subcell) = 1.e-30
    thisOctal%temperature(subcell) = 1.e-3
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .false.

    if ((r > grid%rInner).and.(r < grid%rOuter)) then
       thisOctal%rho(subcell) = rho * (grid%rInner / r)**2 
       thisOctal%temperature(subcell) = 10.
       thisOctal%inFlow(subcell) = .true.
       thisOctal%etaCont(subcell) = 0.
    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine calcTestDensity

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
          thisOctal%rho(subcell) = max(thisOctal%rho(subcell), 1.e-33)
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





  subroutine shakaraDisk(thisOctal,subcell,grid)

    use input_variables
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    real :: r, hr, rd, r1
    TYPE(octalVector) :: rVec
    
    rVec = subcellCentre(thisOctal,subcell)
    r = modulus(rVec)

    thisOctal%rho(subcell) = 1.e-33
    thisOctal%temperature(subcell) = 10.
    thisOctal%etaCont(subcell) = 0.
    thisOctal%inFlow(subcell) = .false.
    rd = rOuter / 2.
    r = sqrt(rVec%x**2 + rVec%y**2)
    if ((r > rInner).and.(r < rOuter)) then
       thisOctal%rho(subcell) = density(rVec, grid)
       thisOctal%rho(subcell) = max(thisOctal%rho(subcell), 1.e-30)
       thisOctal%temperature(subcell) = 10.
       thisOctal%etaCont(subcell) = 0.
    endif
    if (thisOctal%rho(subcell) > 1.e-30) then
       thisOctal%inFlow(subcell) = .true.
    endif
    thisOctal%velocity = VECTOR(0.,0.,0.)
    thisOctal%biasCont3D = 1.
    thisOctal%etaLine = 1.e-30
  end subroutine shakaraDisk

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



  SUBROUTINE addNewChildren(parent, grid, sphData, stellar_cluster)
    ! adds all eight new children to an octal

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


    
    ! For only cluster geometry ...
    TYPE(sph_data), optional, intent(in) :: sphData
    
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
          ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%nLambda))
          ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%nLambda))
          parent%child(newChildIndex)%kappaAbs = 1.e-30
          parent%child(newChildIndex)%kappaSca = 1.e-30
       endif
       ALLOCATE(parent%child(newChildIndex)%N(8,grid%maxLevels))
       NULLIFY(parent%child(newChildIndex)%child)

       ! set up the new child's variables
       parent%child(newChildIndex)%threed = parent%threed
       parent%child(newChildIndex)%twoD = parent%twod
       parent%child(newChildIndex)%maxChildren = parent%maxChildren
       parent%child(newChildIndex)%twoD = parent%twod
       parent%child(newChildIndex)%parent => parent
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
       parent%child(newChildIndex)%N = 1.e-30
       parent%child(newChildIndex)%dusttype = 1

       if (present(sphData)) then
          ! updates the sph particle list.           
          call update_particle_list(parent, newChildIndex, newChildIndex, sphData)
       end if
       
       ! put some data in the eight subcells of the new child
       DO subcell = 1, parent%maxChildren
          CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid, sphData, stellar_cluster)
          parent%child(newChildIndex)%label(subcell) = counter
          counter = counter + 1
       END DO
 
    enddo

    ! check for a new maximum depth 
    IF (parent%child(1)%nDepth > grid%maxDepth) THEN
       grid%maxDepth = parent%child(1)%nDepth
       grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
            2.0_oc**REAL(grid%maxDepth,KIND=octalKind)
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
    

    DEALLOCATE(thisOctal%gas_particle_list)
    
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
          thisOctal%rho(subcell) = max(1.e-30,thisOctal%rho(subcell) * scaleFac)
       endif
    enddo
  end subroutine scaleDensityAMR

  recursive subroutine findTotalMass(thisOctal, totalMass, minRho, maxRho)
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  real(kind=doubleKind) :: totalMass
  real(kind=doubleKind),optional :: minRho, maxRho
  real(kind=doubleKind) :: r1,r2,dv
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
          if (thisOctal%inFlow(subcell)) then
             if (thisOctal%threed) then
                totalMass = totalMass + (1.d30)*thisOctal%rho(subcell) * thisOctal%subcellSize**3
             else
                rVec = subcellCentre(thisOctal,subcell)
                r1 = rVec%x - thisOctal%subcellSize/2.
                r2 = rVec%x + thisOctal%subcellSize/2.
                dv = pi * (r2**2-r1**2) * thisOctal%subcellSize
                totalMass = totalMass + (1.d30)*thisOctal%rho(subcell) * dv
             endif
                
             if (PRESENT(minRho)) minRho = min(dble(thisOctal%rho(subcell)), minRho)
             if (PRESENT(maxRho)) maxRho = max(dble(thisOctal%rho(subcell)), maxRho)
          endif
       endif
    enddo
  end subroutine findTotalMass

!
! Computes average temperature (T_ave) and mass-weighted temperature (T_mass)
!
  subroutine find_average_temperature(grid, T_ave, T_mass, TotalMass)
    implicit none
    type(gridtype), intent(in)  :: grid
    real(kind=doubleKind), intent(out) :: T_ave   ! average temperature in [k}
    real(kind=doubleKind), intent(out) :: T_mass  ! mass weighted temperature in [K]
    real(kind=doubleKind), intent(out) :: TotalMass  ! total mass [g]
    ! 
    real(kind=doublekind) :: sum_T    ! Sum of temperatures in all cells [K]
    real(kind=doublekind) :: sum_M    ! Sum of mass in all cells [g]
    real(kind=doublekind) :: sum_TM   ! Sum of temperatures*Mass in all cells [K*g]
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
    real(kind=doubleKind), intent(inout) :: totalTemperature      ! simple total mass
    real(kind=doubleKind), intent(inout) :: totalMassTemperature
    real(kind=doubleKind), intent(inout) :: totalMass
    integer, intent(inout) :: ncell ! number of cells used for this calculation
    !
    type(octal), pointer  :: child
    real(kind=doubleKind) :: M, T
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
    real(kind=doubleKind), intent(out) :: tau_max
    real(kind=doubleKind), intent(out) :: tau_min
    real(kind=doubleKind), intent(out) :: tau_ave
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
    real(kind=doubleKind), intent(inout) :: tau_max
    real(kind=doubleKind), intent(inout) :: tau_min
    real(kind=doubleKind), intent(inout) :: tau_ave
    integer, intent(inout) :: ncell ! number of cells used for this calculation
    !
    type(octal), pointer  :: child
    real(kind=doubleKind) :: chi, tau
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
       thisOctal%rho(subcell) = rhoOut* returnSpiralFactor(octVec, 10.*grid%rCore/twoPi, grid)
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
    real(kind=doubleKind) :: totalMass, minRho, maxRho
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
    real(kind=doubleKind) :: r1
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


  !
  !
  !

  SUBROUTINE deleteOctreeBranch(topOctal,deletedBranch)
    ! recursively deletes an octal (and any children it has). If the
    !   'deletedBranch' pointer is supplied, the branch is copied 
    !   there as it is deleted.
    ! NB This does not delete the 'topOctal' of the deleted tree - 
    !   you must do this yourself after you call this subroutine!
  
    IMPLICIT NONE

    TYPE(octal), POINTER           :: topOctal      ! top of branch to be deleted
    TYPE(octal), POINTER, OPTIONAL :: deletedBranch ! optional copy of deleted branch
    
    TYPE(octal), POINTER :: deletedBranchcopy
    TYPE(octal), POINTER :: topOctalPointer

    NULLIFY(deletedBranchcopy)
    NULLIFY(topOctalPointer)
    
    ! if required, we copy the octal at the top of the branch
    IF (PRESENT(deletedBranch)) THEN
            
      ALLOCATE(deletedBranch)
      
      ! the octal structure contains allocatable array pointers. we must
      !   create similarly sized arrays in the new octal before copying
      !   the contents from the existing octal.
      
      IF (ASSOCIATED(topOctal%kappaAbs)) THEN                  
        ALLOCATE(deletedBranch%kappaAbs((SIZE(topOctal%kappaAbs,1)),   &
                                        (SIZE(topOctal%kappaAbs,2))))
      ELSE
        NULLIFY(deletedBranch%kappaAbs)
      END IF
          
      IF (ASSOCIATED(topOctal%kappaSca)) THEN                  
        ALLOCATE(deletedBranch%kappaSca((SIZE(topOctal%kappaSca,1)),   &
                                        (SIZE(topOctal%kappaSca,2))))
      ELSE
        NULLIFY(deletedBranch%kappaSca)
      END IF
          
      IF (ASSOCIATED(topOctal%N)) THEN                  
        ALLOCATE(deletedBranch%N((SIZE(topOctal%N,1)),(SIZE(topOctal%N,2))))
      ELSE
        NULLIFY(deletedBranch%N)
      END IF
     
!      IF (ASSOCIATED(topOctal%departCoeff)) THEN                  
!        ALLOCATE(deletedBranch%departCoeff((SIZE(topOctal%departCoeff,1)),(SIZE(topOctal%departCoeff,2))))
!      ELSE
!        NULLIFY(deletedBranch%departCoeff)
!      END IF
     
      ! now that the original and destination octals have identical layouts, we
      !   can copy the data from one to the other.
      deletedBranch = topOctal  
      
      NULLIFY(deletedBranch%parent) 
      NULLIFY(deletedBranch%child)  ! for now, assume no children

    END IF

    ! if the original octal has any children, we must recursively delete them
    !   from the tree (and optionally make copies of them in the
    !   'deletedBranch')
    IF (topOctal%nChildren > 0) THEN
            
      deletedBranchCopy => deletedBranch
      topOctalPointer => topOctal
      
      CALL deleteOctreeBranchPrivate(topOctalPointer,deletedBranch=deletedBranchCopy)
      
      ! we now dellocate the children
      DEALLOCATE(topOctal%child) 
      
    END IF
    
  CONTAINS
    
    RECURSIVE SUBROUTINE deleteOctreeBranchPrivate(thisOctal,deletedBranch)
      ! delete (and optionally make copies) of an octal's children.

      TYPE(octal), POINTER           :: thisOctal
      TYPE(octal), POINTER, OPTIONAL :: deletedBranch

      TYPE(octal), POINTER  :: thisChild          
      TYPE(octal), POINTER  :: thisChildPointer   
      TYPE(octal), POINTER  :: deletedBranchChild 
      TYPE(octal), POINTER  :: deletedBranchChildPointer
      INTEGER               :: iChild
      INTEGER               :: status

      TYPE(octal), POINTER  :: parentPointer1 
 
      NULLIFY(thisChild)
      NULLIFY(thisChildPointer)
      NULLIFY(deletedBranchChild)
      NULLIFY(deletedBranchChildPointer)

      ! if we want to copy the children before deletion, we must first
      !   allocate space.
      IF (PRESENT(deletedBranch)) THEN 
              
        ALLOCATE(deletedBranch%child(thisOctal%nChildren),STAT=status)
        IF (status /= 0) print *, 'allocation status /= 0 in deleteOctreeBranchPrivate'
        
      END IF

      ! now loop over each child
      DO iChild = 1, thisOctal%nChildren, 1
      
        thisChild => thisOctal%child(iChild)
         
        ! again, have to allocate arrays of same size before copying contents.
        IF (PRESENT(deletedBranch)) THEN
          deletedBranchChild => deletedBranch%child(iChild)

          ! 
          IF (ASSOCIATED(thisChild%kappaAbs)) THEN                  
            ALLOCATE(deletedBranchChild%kappaAbs(      &
                                (SIZE(thisChild%kappaAbs,1)),   &
                                (SIZE(thisChild%kappaAbs,2))))
          ELSE 
            NULLIFY(deletedBranchChild%kappaAbs)
          END IF
         
          IF (ASSOCIATED(thisChild%kappaSca)) THEN                  
            ALLOCATE(deletedBranchChild%kappaSca(      &
                                (SIZE(thisChild%kappaSca,1)),   &
                                (SIZE(thisChild%kappaSca,2))))
          ELSE
            NULLIFY(deletedBranchChild%kappaSca)
          END IF
            
          IF (ASSOCIATED(thisChild%N)) THEN                  
            ALLOCATE(deletedBranchChild%N(             &
                                (SIZE(thisChild%N,1)),          &
                                (SIZE(thisChild%N,2))))
          ELSE
            NULLIFY(deletedBranchChild%N)
          END IF
     
!          IF (ASSOCIATED(thisChild%departCoeff)) THEN                  
!            ALLOCATE(deletedBranchChild%departCoeff(             &
!                                (SIZE(thisChild%departCoeff,1)),          &
!                                (SIZE(thisChild%departCoeff,2))))
!          ELSE
!            NULLIFY(deletedBranchChild%departCoeff)
!          END IF
     
          ! can now copy contents
          deletedBranchChild = thisChild  

          ! and set the parent pointer
          deletedBranchChild%parent => deletedBranch
          NULLIFY(deletedBranchChild%child)

            
        END IF



        ! if this child has children of its own, we delete/(copy) those too
        IF (thisChild%nChildren > 0) THEN
                
          IF (PRESENT(deletedBranch)) THEN 
            deletedBranchChildPointer => deletedBranch%child(iChild)
            thisChildPointer => thisOctal%child(iChild)
            CALL deleteOctreeBranchPrivate(thisChildPointer,deletedBranch=deletedBranchChildPointer)
          ELSE
            CALL deleteOctreeBranchPrivate(thisChildPointer)
          END IF
          
          ! can now free-up the memory occupied by this child's child
          IF (ASSOCIATED(thisOctal%child(iChild)%child)) DEALLOCATE(thisOctal%child(iChild)%child)
        END IF
          
      END DO

    END SUBROUTINE deleteOctreeBranchPrivate
  
  END SUBROUTINE deleteOctreeBranch
 

  SUBROUTINE insertOctreeBranch(branch,insertLocation,delete)
    ! copies an orphaned branch back to a location within a tree. If 
    !   'delete' is TRUE, the branch is deleted as it is copied.
    ! NB This subroutine does not set the 'parent' pointer - you must
    !   do this after you call it.

    IMPLICIT NONE
 
    TYPE(octal), POINTER :: branch ! the branch to be inserted
    TYPE(octal), POINTER :: insertLocation ! the octal where the top of the
                                           !   branch will be placed
    LOGICAL, INTENT(IN)  :: delete ! TRUE for 'branch' to be deleted afterwards
    
    TYPE(octal), POINTER :: insertLocationPointer
    
    NULLIFY(insertLocationPointer)
    
    ! the octal structure contains allocatable array pointers. we must
    !   create similarly sized arrays in the insertLocation octal before
    !   copying the contents from the branch octal.
    
    IF (ASSOCIATED(branch%kappaAbs)) THEN                  
      ALLOCATE(insertLocation%kappaAbs(SIZE(branch%kappaAbs,1),SIZE(branch%kappaAbs,2)))
    ELSE
      NULLIFY(insertLocation%kappaAbs)
    END IF
         
    IF (ASSOCIATED(branch%kappaSca)) THEN                  
      ALLOCATE(insertLocation%kappaSca(SIZE(branch%kappaSca,1),SIZE(branch%kappaSca,2)))
    ELSE
      NULLIFY(insertLocation%kappaSca)
    END IF
          
    IF (ASSOCIATED(branch%N)) THEN                  
      ALLOCATE(insertLocation%N(SIZE(branch%N,1),SIZE(branch%N,2)))
    ELSE
      NULLIFY(insertLocation%N)
    END IF
    
!    IF (ASSOCIATED(branch%departCoeff)) THEN                  
!      ALLOCATE(insertLocation%departCoeff(SIZE(branch%departCoeff,1),SIZE(branch%departCoeff,2)))
!    ELSE
!      NULLIFY(insertLocation%departCoeff)
!    END IF
    
    ! now we can copy the data over.
    insertLocation = branch 
    
    NULLIFY(insertLocation%child) ! assume no children
    NULLIFY(insertLocation%parent) ! set this elsewhere
   
    ! recursively copy any children
    IF (branch%nChildren > 0) THEN
            
      insertLocationPointer => insertLocation
      CALL insertOctreeBranchPrivate(branch,insertLocationPointer,delete)
      
    ELSE IF (delete) THEN
      ! optionally delete the octal      
      DEALLOCATE(branch)
    END IF
      
  CONTAINS
    
    RECURSIVE SUBROUTINE insertOctreeBranchPrivate(thisOctal,insertLocation,delete)
      ! insert the octal's children into another location

      TYPE(octal), INTENT(INOUT) :: thisOctal
      TYPE(octal), POINTER       :: insertLocation
      LOGICAL, INTENT(IN)        :: delete
      
      TYPE(octal), POINTER :: thisChild           
      TYPE(octal), POINTER :: insertLocationChild 
      INTEGER              :: iChild
      NULLIFY(thisChild)
      NULLIFY(insertLocationChild)
      ALLOCATE(insertLocation%child(thisOctal%nChildren))
      
      DO iChild = 1, thisOctal%nChildren, 1
      
        thisChild => thisOctal%child(iChild)
        insertLocationChild => insertLocation%child(iChild)
        
        IF (ASSOCIATED(thisChild%kappaAbs)) THEN                  
          ALLOCATE(insertLocationChild%kappaAbs(                    &
                     SIZE(thisChild%kappaAbs,1),SIZE(thisChild%kappaAbs,2)))
        ELSE
          NULLIFY(insertLocationChild%kappaAbs)
        END IF
         
        IF (ASSOCIATED(thisChild%kappaSca)) THEN                  
          ALLOCATE(insertLocationChild%kappaSca(                    &
                     SIZE(thisChild%kappaSca,1),SIZE(thisChild%kappaSca,2)))
        ELSE
          NULLIFY(insertLocationChild%kappaSca)
        END IF
           
        IF (ASSOCIATED(thisChild%N)) THEN                  
          ALLOCATE(insertLocationChild%N(                           &
                     SIZE(thisChild%N,1),SIZE(thisChild%N,2)))
        ELSE
          NULLIFY(insertLocationChild%N)
        END IF
     
!        IF (ASSOCIATED(thisChild%departCoeff)) THEN                  
!          ALLOCATE(insertLocationChild%departCoeff(                           &
!                     SIZE(thisChild%departCoeff,1),SIZE(thisChild%departCoeff,2)))
!        ELSE
!          NULLIFY(insertLocationChild%departCoeff)
!        END IF
     
        insertLocationChild = thisChild  
        insertLocationChild = thisChild  
        insertLocationChild%parent => insertLocation
        
        NULLIFY(insertLocationChild%child)
           
        IF (thisChild%nChildren > 0) THEN
          CALL insertOctreeBranchPrivate(thisChild,insertLocationChild,delete)
        END IF
          
      END DO
      
      IF (delete) DEALLOCATE(thisOctal%child)

    END SUBROUTINE insertOctreeBranchPrivate
  
  END SUBROUTINE insertOctreeBranch

  
  ! chris (19/05/04)
  ! Smooths the AMR grid to ensure 'adequate' resolution at the boundary
  ! between optically thin and optically thick regions. This routine should be
  ! called repeatedly until gridConverged returns true (it should be set to
  ! true in the initial call.
  ! The value of kappa is the opacity of the grid at a test wavelength (set in
  ! the parameter file).

  RECURSIVE SUBROUTINE smoothAMRgridTau(thisOctal,grid,gridConverged, ilam, &
                                        sphData, stellar_cluster)
    use input_variables, only: tauSmoothMax, tauSmoothMin
    IMPLICIT NONE
    TYPE(octal), POINTER             :: thisOctal
    TYPE(gridtype), INTENT(INOUT   ) :: grid 
    LOGICAL, INTENT(INOUT)               :: gridConverged
    TYPE(sph_data), optional, INTENT(IN) :: sphData   ! Matthew's SPH data.
    TYPE(cluster), optional, intent(in)  :: stellar_cluster

    INTEGER              :: i, ilam
    TYPE(octalVector)    :: thisSubcellCentre
    REAL(KIND=octalKind) :: dSubcellCentre
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
        child => thisOctal%child(thisOctal%indexChild(thisSubcell))
        call smoothAMRgridTau(child,grid,gridConverged,  ilam, sphData, stellar_cluster)
      else
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

        thisTau = (grid%oneKappaAbs(thisOctal%dusttype(thissubcell),ilam) + &
                   grid%oneKappaSca(thisOctal%dusttype(thissubcell),ilam)) * &
                   thisOctal%subcellSize * thisOctal%rho(thisSubcell)
        DO i = 1, nLocator, 1
          IF ( inOctal(grid%octreeRoot,locator(i)) ) THEN
            neighbour => thisOctal
            CALL findSubcellLocal(locator(i),neighbour,subcell)

            ! The neighbouring subcell must be larger than the current subcell
            ! otherwise our locators won't cover all neighbouring subcells
            ! (and we'll hit cell boundaries, which is not good).
            IF ( neighbour%subcellSize >= thisOctal%subcellSize) THEN

               thatTau = (grid%oneKappaAbs(neighbour%dusttype(subcell),ilam) + &
                    grid%oneKappaSCa(neighbour%dusttype(subcell),ilam)) * &
                    neighbour%subcellSize * neighbour%rho(subcell)
!               write(*,*) thisTau,thatTau
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
              if (max(thisTau, thatTau).gt.tauSmoothMax.and.min(thisTau, thatTau).lt.tauSmoothMin) then

!               if (abs(thisTau-thatTau) .gt. 1.) then
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

                call addNewChildren(thisOctal, grid, sphData, stellar_cluster)
                ! If we are splitting two subcells in the same octal, we don't
                ! need to addNewChildren to the neighbour. If we switch to
                ! addNewChild this test becomes redundant.
                if (.not. associated(thisOctal, neighbour)) then
                  call addNewChildren(neighbour, grid, sphData, stellar_cluster)
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

  END SUBROUTINE smoothAMRgridTau

  subroutine polarDump(grid)
    implicit none
    type(GRIDTYPE) :: grid
    integer :: nt = 200, nr = 400, i, j
    real :: theta, r, t, dens
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


END MODULE amr_mod

