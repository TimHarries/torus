MODULE amr_mod
  ! routines for adaptive mesh refinement. nhs


  USE vector_mod          ! vector maths routines
  USE kind_mod            ! variable kind parameters    
  USE octal_mod           ! type definition for the grid elements
  USE gridtype_mod        ! type definition for the 3-d grid 
  USE parameters_mod      ! parameters for specific geometries


  IMPLICIT NONE


CONTAINS


  SUBROUTINE calcValuesAMR(thisOctal,subcell,grid)
    ! calculates the variables describing one subcell of an octal.
    ! each geometry that can be used with AMR should be described here, 
    !   otherwise the program will print a warning and exit.
   
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT)    :: thisOctal ! the octal being changed
    INTEGER, INTENT(IN)           :: subcell   ! the subcell being changed
    TYPE(gridtype), INTENT(INOUT) :: grid      ! the grid 


    SELECT CASE (grid%geometry)

    CASE ("ttauri")
      CALL calcTTauriMassVelocity(thisOctal,subcell,grid)
      
    CASE DEFAULT
      WRITE(*,*) "! Unrecognised grid geometry: ",TRIM(grid%geometry)
      STOP

    END SELECT
 
    CALL fillGridDummyValues(thisOctal,subcell)
 
 
  END SUBROUTINE calcValuesAMR


  SUBROUTINE initFirstOctal(grid, centre, size)
    ! creates the first octal of a new grid (the root of the tree).
    ! this should only be used once; use addNewChild for subsequent
    !  additions.

    IMPLICIT NONE
    
    TYPE(gridtype), INTENT(INOUT)    :: grid 
    TYPE(octalVector), INTENT(IN)    :: centre ! coordinates of the grid centre
    REAL, INTENT(IN)                 :: size 
      ! 'size' should be the vertex length of the cube that contains the whole
      !   of the simulation space, *not* the size of a subcell.
                                         
    INTEGER :: subcell ! loop counter 


    ALLOCATE(grid%octreeRoot)
    
    ! allocate any variables that need to be 
    ALLOCATE(grid%octreeRoot%kappaAbs(8,grid%maxLevels))
    ALLOCATE(grid%octreeRoot%kappaSca(8,grid%maxLevels))
    ALLOCATE(grid%octreeRoot%N(8,grid%maxLevels))
    
    grid%octreeRoot%nDepth = 1
    grid%octreeRoot%nChildren = 0
    grid%octreeRoot%hasChild = .FALSE.
    grid%octreeRoot%subcellSize = size/2.0_oc
    grid%octreeRoot%centre = centre
    grid%octreeRoot%indexChild = -999 ! values are undefined
    NULLIFY(grid%octreeRoot%parent)   ! tree root does not have a parent
    NULLIFY(grid%octreeRoot%child)    ! tree root does not yet have children
    
    DO subcell = 1, 8
      ! calculate the values at the centre of each of the subcells
      CALL calcValuesAMR(grid%octreeRoot,subcell,grid)
      ! label the subcells
      grid%octreeRoot%label(subcell) = subcell
    END DO

    ! we keep track of the maximum depth of the grid...
    grid%maxDepth = 1
    ! ...and the size of the smallest subcells in the grid.
    ! we will actually store the value which is half the size of the 
    !   smallest subcell because this is more useful for later
    !   calculations.
    grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
                                2.0_oc**REAL(grid%maxDepth,KIND=octalKind)
      
  END SUBROUTINE initFirstOctal


  SUBROUTINE addNewChild(parent, nChild, grid)
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
       
      ![3] 
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

    ENDIF


    ! update the bookkeeping
    nChildren = nChildren + 1
    newChildIndex = nChildren

    ! update the parent octal
    parent%nChildren = nChildren
    parent%hasChild(nChild) = .TRUE.
    parent%indexChild(newChildIndex) = nChild

    ! allocate any variables that need to be  
    ALLOCATE(parent%child(newChildIndex)%kappaAbs(8,grid%maxLevels))
    ALLOCATE(parent%child(newChildIndex)%kappaSca(8,grid%maxLevels))
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

    ! put some data in the eight subcells of the new child
    DO subcell = 1, 8
      CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid)
      parent%child(newChildIndex)%label(subcell) = counter
      counter = counter + 1
    END DO

 
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


  RECURSIVE SUBROUTINE splitGrid(thisOctal,amrLimitScalar,grid)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    IMPLICIT NONE

    TYPE(OCTAL), POINTER :: thisOctal
    REAL, INTENT(IN)     :: amrLimitScalar 
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
                                          !   routines that this subroutine calls
    
    TYPE(OCTAL), POINTER :: childPointer  
    INTEGER              :: subcell, i    ! loop counters
    
    DO subcell = 1, 8, 1
    
      IF (decideSplit(thisOctal,subcell,amrLimitScalar,grid)) THEN
      
        CALL addNewChild(thisOctal,subcell,grid)

        ! find the index of the new child and call splitGrid on it
        DO i = 1, 8, 1
          IF ( thisOctal%indexChild(i) == subcell ) THEN 
            childPointer => thisOctal%child(i)
            CALL splitGrid(childPointer,amrLimitScalar,grid)
            EXIT
          END IF
        END DO
        
      END IF
      
    END DO

  END SUBROUTINE splitGrid
  
  
  RECURSIVE SUBROUTINE getOctalArray(thisOctal,array,counter) 
  ! returns an array of pointers to all of the subcells in the grid.
    ! NB because fortran cannot create arrays of pointers, the output
    !   array is actually of a derived type which *contains* the 
    !   pointer to an octal.

    IMPLICIT NONE

    TYPE(octal), POINTER                            :: thisOctal
    TYPE(octalWrapper), DIMENSION(:), INTENT(INOUT) :: array 
    INTEGER, INTENT(INOUT)                          :: counter 
  
    INTEGER              :: i
    TYPE(octal), POINTER :: child


    counter = counter + 1 
    array(counter)%content => thisOctal
    array(counter)%inUse = .TRUE. 
    
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
    TYPE(octalWrapper), DIMENSION(:), ALLOCATABLE :: octalArray
  
    INTEGER :: subcell, i 
    INTEGER :: counter ! used by getOctalArray
     
    ! all of the work that must be done recursively goes here: 
    DO subcell = 1, 8, 1
   
      SELECT CASE (grid%geometry)

      CASE ("ttauri")
        CALL calcTTauriTemperature(thisOctal,subcell)
        gridConverged = .TRUE.
        
      CASE DEFAULT
        WRITE(*,*) "! Unrecognised grid geometry: ",trim(grid%geometry)
        STOP
      
      END SELECT
      
      IF ( thisOctal%hasChild(subcell) )  THEN
      ! find the index of the child and call finishGrid on it
        DO i = 1, 8, 1
          IF ( thisOctal%indexChild(i) == subcell ) THEN 
            child => thisOctal%child(i)
            CALL finishGrid(child,grid,gridConverged)
            EXIT 
          END IF
        END DO
      END IF
      
    END DO

    ! any stuff that gets done *after* the finishGrid recursion goes here:
    IF (.NOT. ASSOCIATED(thisOctal%parent)) THEN

       ! we create an array of pointers to the all of the array's octals.
       ALLOCATE(octalArray(grid%nOctals))
       counter = 0 
       CALL getOctalArray(grid%octreeRoot,octalArray,counter) 
       !
       ! [do the statistical equilibrium setup here]
       !
       DEALLOCATE(octalArray)

    END IF


  END SUBROUTINE finishGrid
 
  
  RECURSIVE SUBROUTINE countVoxels(thisOctal,nOctals,nVoxels)  
    ! count the number of octals in the current section of the grid.
    ! also counts the number of unique volume elements (voxels) i.e.
    !   those subcells that are not subdivided

    IMPLICIT NONE

    TYPE(OCTAL), POINTER  :: thisOctal 
    INTEGER,INTENT(INOUT) :: nOctals   ! number of octals
    INTEGER,INTENT(INOUT) :: nVoxels   ! number of childless subcells
    
    TYPE(OCTAL), POINTER  :: child

    INTEGER :: i

    nOctals = nOctals + 1

    IF ( thisOctal%nChildren > 0 ) THEN
      ! call this subroutine recursively on each of its children
      DO i = 1, thisOctal%nChildren, 1
      
        child => thisOctal%child(i)
        CALL countVoxels(child,nOctals,nVoxels)
      END DO
    END IF

    ! increment the counter once for each of its childless subcells
    nVoxels = nVoxels + (8 - thisOctal%nchildren)

  END SUBROUTINE countVoxels


  SUBROUTINE startReturnSamples (startPoint,direction,grid, &
             sampleFreq,nSamples,maxSamples,lambda,kappaAbs,&
             kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,&
             opaqueCore,hitCore,usePops,iLambda,error)
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
                       ! samples that will be made in any subcell of the octree. 
                       
    INTEGER, INTENT(INOUT)             :: nSamples   ! number of samples made
    INTEGER, INTENT(IN)                :: maxSamples ! size of sample arrays 
    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda     ! path distances of sample locations
    REAL, DIMENSION(:), INTENT(INOUT)  :: kappaAbs   ! continuous absorption opacities
    REAL, DIMENSION(:), INTENT(INOUT)  :: kappaSca   ! scattering opacities
    TYPE(vector),DIMENSION(:),INTENT(INOUT) :: velocity ! sampled velocities
    REAL, DIMENSION(:), INTENT(INOUT)  :: velocityDeriv ! sampled velocity derivatives
    REAL, DIMENSION(:), INTENT(INOUT)  :: chiLine    ! line opacities
    REAL, DIMENSION(:,:), INTENT(INOUT):: levelPop   ! level populations
    REAL, DIMENSION(:)                 :: rho        ! density at sample points
    LOGICAL, INTENT(IN)                :: opaqueCore ! true if the core is opaque
    LOGICAL, INTENT(OUT)               :: hitCore    ! true if the core is opaque
    LOGICAL, INTENT(IN)                :: usePops    ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda    ! wavelength index
    INTEGER, INTENT(INOUT)             :: error      ! error code

    TYPE(octalVector)       :: locator 
                       ! 'locator' is used to indicate a point that lies within the  
                       !   *next* cell of the octree that the ray will interesect.
                       !   initially this will be the same as the startPoint
      
    TYPE(octalVector)       :: currentPoint ! current position of ray 
    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
    TYPE(octal)             :: octree       ! the octree structure within 'grid'
    TYPE(octalVector)       :: directionNormalized
    REAL(KIND=octalKind)    :: margin
                       ! margin is the size of the region around the edge of a subcell
                       !   where numerical inaccuracies may cause problems.
    REAL(KIND=octalKind)    :: distanceLimit ! max length of ray before aborting
 
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
    
    
    ! set up some variables
    octree = grid%octreeRoot  
    abortRay = .FALSE.
    hitCore = .FALSE.
    directionNormalized = direction
    CALL normalize(directionNormalized)
    distanceLimit = HUGE(distanceLimit)
    locator = startPoint
    
    ! we will abort tracking any rays which are too close to 
    !   to a cell wall. The distance to use is defined by:
    margin = 8.0_oc * REAL(grid%maxDepth,KIND=octalKind) * EPSILON(1.0_oc)
    ! some more experimentation is required to find out the best value for
    !   margin
    
    currentPoint = startPoint
    
    IF (.NOT. inOctal(octree,startPoint)) THEN
      PRINT *, 'Attempting to find path between point(s) outwith the grid'
      STOP
    ENDIF
   
   
    ! geometry-specific tests should go here
    SELECT CASE (grid%geometry)
      
    CASE ("ttauri")
       ! need to test for both star and disc intersections 
      
       ! we will find out when and where the photon leaves the simulation space 
       CALL getExitPoint(currentPoint,directionNormalized,locator,octree,&
                         abortRay,error,grid%halfSmallestSubcell,endPoint,&
                         octree%centre,endLength,margin) 
       
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
       CALL returnSamples(currentPoint,startPoint,locator,directionNormalized,octree, &
                 grid,sampleFreq,nSamples,maxSamples,&
                 abortRay,lambda,kappaAbs,kappaSca,velocity,&
                 velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,error,&
                 margin,endLength)
      
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
           
         CALL returnSamples(currentPoint,dummystartPoint,locator,directionNormalized,&
                   octree,grid,sampleFreq,nSamples,maxSamples,&
                   abortRay,lambda,kappaAbs,kappaSca,velocity,&
                   velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,error,&
                   margin,distanceLimit)

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
       
    CASE DEFAULT       
       CALL returnSamples(currentPoint,startPoint,locator,directionNormalized,octree, &
                 grid,sampleFreq,nSamples,maxSamples,&
                 abortRay,lambda,kappaAbs,kappaSca,velocity,&
                 velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,error,&
                 margin,distanceLimit)
    END SELECT
      
   


  END SUBROUTINE startReturnSamples

  
  RECURSIVE SUBROUTINE returnSamples (currentPoint,startPoint,locator,direction, &
             octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,kappaAbs, &
             kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,&
             usePops,iLambda,error,margin,distanceLimit)
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
                  ! 'sampleFreq' is the maximum number of samples that will be made
                  !   in any subcell of the octree. 
    INTEGER, INTENT(INOUT)              :: nSamples     ! number of samples made
    INTEGER, INTENT(IN)                 :: maxSamples   ! size of sample arrays 
    LOGICAL, INTENT(INOUT)              :: abortRay     ! flag to signal completion
    REAL, DIMENSION(:), INTENT(INOUT)   :: lambda       ! distance travelled by photon
    REAL, DIMENSION(:), INTENT(INOUT)   :: kappaAbs     ! continuous absorption opacities
    REAL, DIMENSION(:), INTENT(INOUT)   :: kappaSca     ! scattering opacities
    TYPE(vector), DIMENSION(:), INTENT(INOUT) :: velocity ! sampled velocities
    REAL, DIMENSION(:), INTENT(INOUT)   :: velocityDeriv ! sampled velocity derivatives
    REAL, DIMENSION(:), INTENT(INOUT)   :: chiLine      ! line opacities
    REAL, DIMENSION(:,:), INTENT(INOUT) :: levelPop     ! level populations 
    REAL, DIMENSION(:)                  :: rho          ! density at sample points
    LOGICAL, INTENT(IN)                 :: usePops      ! whether to use level populations
    INTEGER, INTENT(IN)                 :: iLambda      ! wavelength index
    INTEGER, INTENT(INOUT)              :: error        ! error code
    REAL(KIND=octalKind)                :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    REAL(KIND=octalKind), INTENT(IN)    :: distanceLimit ! max length of ray before aborting


    TYPE(octalVector)      :: exitPoint      ! where ray leaves current cell
    TYPE(octalVector)      :: centre         ! centre of current subcell
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

        CALL returnSamples(currentPoint,startPoint,locator,direction, &
                    octree%child(subIndex),grid,sampleFreq,nSamples,&
                    maxSamples,abortRay,lambda,&
                    kappaAbs,kappaSca,velocity,velocityDeriv,chiLine,&
                    levelPop,rho,usePops,iLambda,error,margin,distanceLimit)
        
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

        ! we find the exit point from the current subcell, 
        !  and also 'locator' - a point that lies in the *next* subcell 
        CALL getExitPoint(currentPoint,direction,locator,octree,abortRay,error,&
                          grid%halfSmallestSubcell,exitPoint,centre,minWallDistance,&
                          margin)
        IF ( abortRay .EQV. .TRUE. ) RETURN

        ! we now decide where we are going to sample the quantities
        ! we will define an approximate rate that is a fraction of
        !   the current subcell size
        sampleLength = octree%subcellSize / sampleFreq

        ! if there are no previously sampled points, we definitely have to take a
        !   sample here
        IF ( nSamples == 0 ) &
          CALL takeSample(currentPoint,length,direction,grid,octree,subcell,&
                          nSamples,maxSamples,lambda,kappaAbs,kappaSca,velocity,&
                          velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,&
                          error) 

        
        ! we check whether we should take a sample at currentPoint
        !   (which will usually be the entry point to the subcell).

        length = modulus(currentPoint - startPoint)

        IF ( modulus(currentPoint - (startPoint + (direction *         & 
                              REAL(lambda(nSamples),KIND=octalKind)))) &
                                               > sampleLength ) THEN
          CALL takeSample(currentPoint,length,direction,grid,octree,subcell,&
                          nSamples,maxSamples,lambda,kappaAbs,kappaSca,velocity,&
                          velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,error) 
        END IF

        ! we add sampleLength to the distance from the last location
        !   that was sampled, and decide whether to take a new sample.
        trial = 1
        DO 
          trialPoint = currentPoint + direction * &
                       (length - (REAL(lambda(nSamples),KIND=octalKind) + &
                       (REAL(trial,KIND=octalKind) * sampleLength)))

          ! we only want to take a sample if we are still within the subcell
          IF ( modulus(trialPoint - currentPoint) < minWallDistance ) THEN
          
            trialLength = length + modulus(trialPoint - currentPoint)
            IF (trialLength > distanceLimit) THEN
              abortRay = .TRUE.
              RETURN
            END IF
              
            IF (trialLength >= 0.0_oc) THEN
              CALL takeSample(trialPoint,trialLength,direction,grid, &
                       octree,subcell,nSamples,maxSamples,lambda,kappaAbs,&
                       kappaSca,velocity,velocityDeriv,chiLine,levelPop,&
                       rho,usePops,iLambda,error) 
            ELSE
              EXIT
            END IF
            
          ELSE 
            EXIT
          END IF
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

  
  SUBROUTINE getExitPoint(currentPoint,direction,locator,octree,abortRay,error,&
                          halfSmallestSubcell,exitPoint,centre,minWallDistance,&
                          margin)
    ! this subroutine finds the closest subcell wall in the direction the
    !   photon is travelling. 

    IMPLICIT NONE
          
    TYPE(octalVector), INTENT(IN)       :: currentPoint ! current ray position
    TYPE(octalVector), INTENT(IN)       :: direction    ! ray's direction vector 
    TYPE(octalVector), INTENT(OUT)      :: locator       
                  ! 'locator' is used to indicate a point that lies within the  
                  !   *next* cell of the octree that the ray will interesect.
    TYPE(octal), INTENT(IN)             :: octree       ! an octal grid
    LOGICAL, INTENT(INOUT)              :: abortRay     ! flag to signal completion
    INTEGER, INTENT(INOUT)              :: error        ! error code
    REAL(KIND=octalKind), INTENT(IN)    :: halfSmallestSubcell
    REAL(KIND=octalKind), INTENT(IN)    :: margin
                  ! margin is the size of the region around the edge of a subcell
                  !   where numerical inaccuracies may cause problems.
    TYPE(octalVector), INTENT(OUT)      :: exitPoint    ! where ray leaves current cell
    TYPE(octalVector), INTENT(IN)       :: centre       ! centre of current subcell
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

       
    ! this code is ugly, but should work OK.   
                  

    ! there are six subcell walls - and each wall is part of a plane    
    !   parallel with an axis. we use the direction of the photon to 
    !   find the three planes that the photon will intersect.
    
    IF ( direction%x > 0.0_oc ) THEN
      walldistanceX = (centre%x - currentPoint%x + octree%subcellSize / 2.0_oc  ) / ABS(direction%x) 
    ELSE IF ( direction%x < 0.0_oc ) THEN
      wallDistanceX = (currentPoint%x - centre%x + octree%subcellSize / 2.0_oc ) / ABS(direction%x) 
    ELSE 
      wallDistanceX = HUGE(wallDistanceX) 
    END IF

    IF ( direction%y > 0.0_oc ) THEN
      wallDistanceY = (centre%y - currentPoint%y + octree%subcellSize / 2.0_oc ) / ABS(direction%y) 
    ELSE IF ( direction%y < 0.0_oc ) THEN
      wallDistanceY = (currentPoint%y - centre%y + octree%subcellSize / 2.0_oc ) / ABS(direction%y) 
    ELSE 
      wallDistanceY = HUGE(wallDistanceX) 
    END IF
        
    IF ( direction%z > 0.0_oc ) THEN
      wallDistanceZ = (centre%z - currentPoint%z + octree%subcellSize / 2.0_oc ) / ABS(direction%z) 
    ELSE IF ( direction%z < 0.0_oc ) THEN
      wallDistanceZ = (currentPoint%z - centre%z + octree%subcellSize / 2.0_oc ) / ABS(direction%z) 
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
      PRINT *, 'In getExitPoint, aborting ray because distance is less than ''margin'''
      PRINT *, '  wallDistances=',wallDistanceX,wallDistanceY,wallDistanceZ 
      error = -20
      RETURN
    END IF
       
       
    IF ( wallDistanceX < wallDistanceY .AND. &
         wallDistanceX < wallDistanceZ ) THEN
      wallNormal =  xHatOctal
      IF ( direction%x > 0.0_oc ) THEN
        wallFromOrigin = centre%x + (octree%subcellSize / 2.0_oc)
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
        wallFromOrigin = centre%x - (octree%subcellSize / 2.0_oc)
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
        wallFromOrigin = centre%y + (octree%subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*yHatOctal)
      ELSE
        wallFromOrigin = centre%y - (octree%subcellSize / 2.0_oc) 
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
       wallFromOrigin = centre%z + (octree%subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*zHatOctal)
      ELSE
        wallFromOrigin = centre%z - (octree%subcellSize / 2.0_oc)
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
        wallFromOrigin = centre%x + (octree%subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHatOctal)
      ELSE
        wallFromOrigin = centre%x - (octree%subcellSize / 2.0_oc)
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
        wallFromOrigin = centre%x + (octree%subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHatOctal)
      ELSE
        wallFromOrigin = centre%x - (octree%subcellSize / 2.0_oc)
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
        wallFromOrigin = centre%y + (octree%subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                     (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*yHatOctal)
      ELSE
        wallFromOrigin = centre%y - (octree%subcellSize / 2.0_oc)
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
        wallFromOrigin = centre%x + (octree%subcellSize / 2.0_oc)
        exitPoint = intersectionLinePlane &
                   (currentPoint,direction,wallNormal,wallFromOrigin,found)
        IF ( found .EQV. .FALSE. ) THEN 
          PRINT *, "Panic: wall intersection not found"
          DO ; END DO ; STOP ! sit in loop for debugging purposes
        END IF
        locator = exitpoint + (halfSmallestSubcell*xHatOctal)
      ELSE
        wallFromOrigin = centre%x - (octree%subcellSize / 2.0_oc)
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

  END SUBROUTINE getExitPoint 


  SUBROUTINE takeSample(point,length,direction,grid,thisOctal,subcell,nSamples,&
                        maxSamples,lambda,kappaAbs,kappaSca,velocity,&
                        velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,error) 
  
    
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
    REAL, DIMENSION(:), INTENT(INOUT)  :: kappaAbs  ! continuous absorption opacities
    REAL, DIMENSION(:), INTENT(INOUT)  :: kappaSca  ! scattering opacities 
    TYPE(vector), DIMENSION(:), INTENT(INOUT) :: velocity ! sampled velocities
    REAL, DIMENSION(:), INTENT(INOUT)  :: velocityDeriv   ! sampled velocity derivatives
    REAL, DIMENSION(:), INTENT(INOUT)  :: chiLine   ! line opacities
    REAL, DIMENSION(:,:), INTENT(INOUT):: levelPop  ! level populations
    REAL, DIMENSION(:), INTENT(INOUT)  :: rho       ! density at sample points
    LOGICAL, INTENT(IN)                :: usePops   ! whether to use level populations
    INTEGER, INTENT(IN)                :: iLambda   ! wavelength index
    INTEGER, INTENT(INOUT)             :: error     ! error code

    TYPE(vector)                       :: directionReal ! direction vector (REAL values)
    TYPE(octal), POINTER               :: localPointer  ! pointer to the current octal

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
                         grid=grid)
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
                         grid=grid)
    END IF

    ! use variables to silence compiler warnings
    error = error

  END SUBROUTINE takeSample

  
  SUBROUTINE amrGridValues(octalTree,point,startOctal,foundOctal,&
                           foundSubcell,actualSubcell,iLambda,direction,&
                           velocity,velocityDeriv,temperature,kappaAbs,&
                           kappaSca,rho,chiLine,etaLine,etaCont, &
                           probDistLine,probDistCont,N,Ne,nTot,grid,interp)
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
    TYPE(octal), OPTIONAL, POINTER    :: startOctal
    TYPE(octal), OPTIONAL, POINTER    :: foundOctal
    INTEGER, INTENT(OUT), OPTIONAL    :: foundSubcell
    INTEGER, INTENT(IN), OPTIONAL     :: actualSubcell 
    INTEGER, INTENT(IN), OPTIONAL     :: iLambda       ! wavelength index
    TYPE(vector),INTENT(IN),OPTIONAL  :: direction     
    TYPE(gridtype),INTENT(IN),OPTIONAL:: grid          
    LOGICAL, INTENT(IN), OPTIONAL     :: interp        ! use interpolation
                                      !  ^^^^^^ ! not implemented yet???

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
    
    
    TYPE(octal), POINTER              :: resultOctal
    INTEGER                           :: subcell
    LOGICAL                           :: interpolate

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
     
      IF (PRESENT(kappaAbs)) THEN 
        IF (PRESENT(iLambda)) THEN
          kappaAbs = resultOctal%kappaAbs(subcell,1)
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''kappaAbs'' without',&
                   ' ''iLambda''.'
          STOP
        END IF
      END IF
      
      IF (PRESENT(kappaSca)) THEN 
        IF (PRESENT(iLambda)) THEN
          kappaSca = resultOctal%kappaSca(subcell,iLambda)
        ELSE
          PRINT *, 'In amrGridValues, can''t evaluate ''kappaSca'' without',&
                   ' ''iLambda''.'
          STOP
        END IF
      END IF

      IF (PRESENT(velocityDeriv)) THEN
        IF (PRESENT(direction) .AND. PRESENT(grid)) THEN
          velocityDeriv = amrGridDirectionalDeriv(grid,point,direction,&
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
    REAL(KIND=octalKind)           :: inc
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
      
      t1 = point%x - (centre%x - inc)
      t2 = point%y - (centre%y - inc)
      t3 = point%z - (centre%z - inc)
      
      amrGridVelocity = &
        ((1.d0-t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(subcell,1) + &
        ((     t1) * (1.d0-t2) * (1.d0-t3)) * resultOctal%cornerVelocity(subcell,2) + &
        ((1.d0-t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(subcell,3) + &
        ((     t1) * (     t2) * (1.d0-t3)) * resultOctal%cornerVelocity(subcell,4) + &
        ((1.d0-t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(subcell,5) + &
        ((     t1) * (1.d0-t2) * (     t3)) * resultOctal%cornerVelocity(subcell,6) + &
        ((1.d0-t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(subcell,7) + &
        ((     t1) * (     t2) * (     t3)) * resultOctal%cornerVelocity(subcell,8)
      

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

    IF (dx /=0.) THEN
       amrGridDirectionalDeriv = (abs(dphi / dx))
    ELSE
       amrGridDirectionalDeriv = 1.e-10
    ENDIF

    IF (amrGridDirectionalDeriv == 0.) amrGridDirectionalDeriv = 1.e-10


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
      amrGridTemperature = resultOctal%temperature(subcell)
      
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
      amrGridDensity = resultOctal%rho(subcell)
      
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

    TYPE(octal), POINTER :: child
    INTEGER              :: i, j 
   
   
    ! we need to treat the first subcell as a special case
    
    IF (probability < thisOctal%probDistLine(1)) THEN

      IF (thisOctal%hasChild(1)) THEN
      
        ! find the child
        DO j = 1, thisOctal%nChildren, 1
          IF (thisOctal%indexChild(j) == 1) THEN
            child => thisOctal%child(j)
            CALL locateLineProbAMR(probability,child,subcell)
            RETURN
          END IF
        END DO
        
      ELSE 
        subcell = 1
        RETURN
        
      END IF
    END IF   
  
   
    DO i = 2, 8, 1
      
      IF (probability > thisOctal%probDistLine(i-1) .AND. &
          probability < thisOctal%probDistLine(i)) THEN
      
        IF (thisOctal%hasChild(i)) THEN
          
          ! find the child
          DO j = 1, thisOctal%nChildren, 1
            IF (thisOctal%indexChild(j) == i) THEN
              child => thisOctal%child(j)
              CALL locateLineProbAMR(probability,child,subcell)
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

    TYPE(octal), POINTER :: child
    INTEGER              :: i, j 
   

   
    ! we need to treat the first subcell as a special case
    
    IF (probability < thisOctal%probDistCont(1)) THEN
    
      IF (thisOctal%hasChild(1)) THEN
      
        ! find the child
        DO j = 1, thisOctal%nChildren, 1
          IF (thisOctal%indexChild(j) == 1) THEN
            child => thisOctal%child(j)
            CALL locateContProbAMR(probability,child,subcell)
            RETURN
          END IF
        END DO
        
      ELSE 
        subcell = 1
        RETURN
        
      END IF
    END IF   
  
   
    DO i = 2, 8, 1
    
      IF (probability > thisOctal%probDistCont(i-1) .AND. &
          probability < thisOctal%probDistCont(i)) THEN
      
        IF (thisOctal%hasChild(i)) THEN
          
          ! find the child
          DO j = 1, thisOctal%nChildren, 1
            IF (thisOctal%indexChild(j) == i) THEN
              child => thisOctal%child(j)
              CALL locateContProbAMR(probability,child,subcell)
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


  PURE FUNCTION whichSubcell(thisOctal,point) RESULT (subcell)
    ! returns the identification number (1-8) of the subcell of the 
    ! current octal which contains a given point
    ! NB this does NOT check that the point lies within the bounds of the octal!
  
    IMPLICIT NONE
    
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(octalVector), INTENT(IN) :: point
    INTEGER                       :: subcell

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

  END FUNCTION whichSubcell    


  PURE FUNCTION inOctal(thisOctal,point) 
    ! true if the point lies within the boundaries of the current octal
  
    IMPLICIT NONE
 
    LOGICAL                       :: inOctal
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(octalVector), INTENT(IN) :: point

    IF ((point%x <= thisOctal%centre%x - thisOctal%subcellSize ) .OR. &
        (point%x >= thisOctal%centre%x + thisOctal%subcellSize ) .OR. &
        (point%y <= thisOctal%centre%y - thisOctal%subcellSize ) .OR. &
        (point%y >= thisOctal%centre%y + thisOctal%subcellSize ) .OR. &
        (point%z <= thisOctal%centre%z - thisOctal%subcellSize ) .OR. &
        (point%z >= thisOctal%centre%z + thisOctal%subcellSize )) THEN
      inOctal = .FALSE.
    ELSE  
      inOctal = .TRUE.
    ENDIF
  
  END FUNCTION inOctal


  FUNCTION looseInOctal(thisOctal,point) 
    ! true if the point lies 'loosely' in the current octal
    ! ( a 10% margin of error is allowed )
    ! this is useful only for testing purposes!
  
    IMPLICIT NONE
 
    LOGICAL                       :: looseInOctal
    TYPE(octal), INTENT(IN)       :: thisOctal
    TYPE(octalVector), INTENT(IN) :: point

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

  END FUNCTION looseInOctal


  RECURSIVE SUBROUTINE smoothAMRgrid(thisOctal,grid,factor,gridConverged)
    ! checks whether each octal's neighbours are much bigger than it, 
    !   if so, makes the neighbours smaller.
    ! the 'needRestart' flag will be set if a change is made to the octree.
    !   smoothAMRgrid should be called repeatedly until the flag is no 
    !   longer set - this should ensure that the octree eventually stabilizes
  
    IMPLICIT NONE

    TYPE(octal), POINTER             :: thisOctal
    TYPE(gridtype), INTENT(INOUT   ) :: grid 
    REAL, INTENT(IN)                 :: factor
    LOGICAL, INTENT(INOUT)           :: gridConverged
    
    INTEGER              :: i
    REAL(KIND=octalKind) :: halfSmallestSubcell
    TYPE(octal), POINTER :: child
    TYPE(octal), POINTER :: neighbour
    TYPE(octalVector), ALLOCATABLE, DIMENSION(:) :: locator
    INTEGER              :: subcell
    
    gridConverged = .TRUE.

    ! we will find the coordinates of a point that lies outside the current
    !   octal. we then compare the size of the cell that contains that point
    !   with the size of the current cell, if it is bigger be more than a 
    !   factor of 'factor', we subdivide the neighbouring cell.
    ! we do this in each of six directions

    ! we do not have to test the other subcells in the current octal because
    !   they can be smaller than the any of the other subcells, but they
    !   cannot be *bigger*. this saves some time.

    ALLOCATE(locator(6))
    
    FORALL (i = 1:6)
      locator(i) = thisOctal%centre
    END FORALL

    ! we find points which are outside the current octal by a distance
    !   equivalent to half the size of the tree's smallest subcell.
    halfSmallestSubcell = grid%halfSmallestSubcell

    locator(1)%x = thisOctal%centre%x + thisOctal%subcellSize + halfSmallestSubcell
    locator(2)%y = thisOctal%centre%y + thisOctal%subcellSize + halfSmallestSubcell
    locator(3)%z = thisOctal%centre%z + thisOctal%subcellSize + halfSmallestSubcell
    locator(4)%x = thisOctal%centre%x - thisOctal%subcellSize - halfSmallestSubcell
    locator(5)%y = thisOctal%centre%y - thisOctal%subcellSize - halfSmallestSubcell
    locator(6)%z = thisOctal%centre%z - thisOctal%subcellSize - halfSmallestSubcell

    DO i = 1, 6, 1
      IF ( inOctal(grid%octreeRoot,locator(i)) ) THEN
        CALL findSubcellTD(locator(i),grid%octreeRoot,neighbour,subcell)
        IF ( neighbour%subcellSize > (factor * thisOctal%subcellSize) ) THEN
          IF ( neighbour%hasChild(subcell) ) THEN 
            PRINT *, "neighbour already has child"
            STOP
          END IF
          CALL addNewChild(neighbour,subcell,grid)
          gridConverged = .FALSE.
        ENDIF
      END IF
    END DO

    ! call this subroutine recursively on each of any children.
    IF ( thisOctal%nChildren > 0 ) THEN
      DO i = 1, thisOctal%nChildren, 1 
        child => thisOctal%child(i)
        CALL smoothAMRgrid(child,grid,factor,gridConverged)
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
      DO i = 1, 8, 1
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
           DO ; END DO
        END IF
     
        ! if we have previously gone down the tree, and are now going back up, there
        !   must be a problem.
        IF (haveDescended) boundaryProblem = .TRUE.
        
        thisOctal => thisOctal%parent
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
       y1 = REAL(thisOctal%centre%y-thisOctal%subcellSize)
       call pgmove(x1,y1)
       x1 = REAL(thisOctal%centre%x+thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%y-thisOctal%subcellSize)
       call pgdraw(x1,y1)
       x1 = REAL(thisOctal%centre%x+thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%y+thisOctal%subcellSize)
       call pgdraw(x1,y1)
       x1 = REAL(thisOctal%centre%x-thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%y+thisOctal%subcellSize)
       call pgdraw(x1,y1)
       x1 = REAL(thisOctal%centre%x-thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%y-thisOctal%subcellSize)
       call pgdraw(x1,y1)

       x1 = REAL(thisOctal%centre%x)
       y1 = REAL(thisOctal%centre%y-thisOctal%subcellSize)
       call pgmove(x1,y1)
       x1 = REAL(thisOctal%centre%x)
       y1 = REAL(thisOctal%centre%y+thisOctal%subcellSize)
       call pgdraw(x1,y1)
       
       x1 = REAL(thisOctal%centre%x-thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%y)
       call pgmove(x1,y1)
       x1 = REAL(thisOctal%centre%x+thisOctal%subcellSize)
       y1 = REAL(thisOctal%centre%y)
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

    DO sub = 1, 8, 1
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


    DO sub = 1, 8, 1
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

    DO sub = 1, 8, 1
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

    DO sub = 1, 8, 1
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

  RECURSIVE SUBROUTINE plotThreeD(thisOctal, viewVec, thisDepth)

    IMPLICIT NONE

    integer :: j
    type(OCTAL), INTENT(IN) :: thisOctal
    type(octalVector) :: viewVec
    integer :: thisDepth

    do j = 1, thisOctal%nChildren
       if ((thisOctal%nDepth == thisDepth))then
  !        call plotCube(viewVec,subcellCentre(thisOctal,j), &
  !                      thisOctal%subcellsize,thisOctal%nDepth)
       endif
    enddo

    do j = 1, 8
       if ((thisOctal%nDepth == thisDepth) .AND. &
               (.NOT. (thisOctal%hasChild(j))) )then
          call plotCube(viewVec,subcellCentre(thisOctal,j), &
                        thisOctal%subcellsize,thisOctal%nDepth)
       endif
    enddo

    do j = 1, thisOctal%nChildren
          call plotThreed(thisOctal%child(j), viewVec, thisDepth)
    enddo


  END SUBROUTINE plotThreeD

  SUBROUTINE plotCube(viewVecOctal, cubeCentreOctal, cubeSizeOctal, nCol)
    
    IMPLICIT NONE
    
    type(octalVector) :: viewVecOctal, cubeCentreOctal
    type(Vector) :: viewVec, cubeCentre
    real(KIND=octalKind) :: cubeSizeOctal
    real :: cubeSize, s
    type(VECTOR) :: xAxis, yAxis, zAxis, xProj, yProj
    real :: xCube(8), yCube(8)
    type(VECTOR) :: cubeCorner(8)
    integer :: i, nCol
    call pgsci(nCol)

    xAxis = VECTOR(1., 0., 0.)
    yAxis = VECTOR(0., 1., 0.)
    zAxis = VECTOR(0., 0., 1.)

    cubesize = cubeSizeOctal

    viewVec%x = viewVecOctal%x
    viewVec%y = viewVecOctal%y
    viewVec%z = viewVecOctal%z
    
    cubeCentre%x = cubeCentreOctal%x
    cubeCentre%y = cubeCentreOctal%y
    cubeCentre%z = cubeCentreOctal%z
    
    xProj =  zAxis .cross. viewVec
    call normalize(xProj)
    yProj = viewVec .cross. xProj
    call normalize(yProj)

    s = cubeSize 

    cubeCorner(1) = cubeCentre + (s * xAxis) + (s * yAxis) + (s * zAxis)
    cubeCorner(2) = cubeCentre + (s * xAxis) - (s * yAxis) + (s * zAxis)
    cubeCorner(3) = cubeCentre - (s * xAxis) - (s * yAxis) + (s * zAxis)
    cubeCorner(4) = cubeCentre - (s * xAxis) + (s * yAxis) + (s * zAxis)
    cubeCorner(5) = cubeCentre + (s * xAxis) + (s * yAxis) - (s * zAxis)
    cubeCorner(6) = cubeCentre + (s * xAxis) - (s * yAxis) - (s * zAxis)
    cubeCorner(7) = cubeCentre - (s * xAxis) - (s * yAxis) - (s * zAxis)
    cubeCorner(8) = cubeCentre - (s * xAxis) + (s * yAxis) - (s * zAxis)
    forall (i = 1:8)
       xCube(i) = cubeCorner(i) .dot. xProj
       yCube(i) = cubeCorner(i) .dot. yProj
    end forall

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
  END SUBROUTINE plotCube


  FUNCTION decideSplit(thisOctal,subcell,amrLimitScalar,grid) RESULT(split)
    ! returns true if the current voxel is to be subdivided. 
    ! decision is made by comparing 'amrLimitScalar' to some value
    !   derived from information in the current cell  

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    REAL, INTENT(IN)           :: amrLimitScalar ! used for split decision
    TYPE(gridtype), INTENT(IN) :: grid
    LOGICAL                    :: split          

    REAL                  :: criticalValue
    REAL(KIND=doubleKind) :: criticalValueDouble
    REAL(KIND=octalKind)  :: cellSize
    TYPE(octalVector)     :: searchPoint
    TYPE(octalVector)     :: cellCentre
    REAL                  :: x, y, z
    INTEGER               :: i
   

    IF ( grid%geometry == "ttauri" ) THEN
    
      ! the density is only sampled at the centre of the grid
      ! we will search in each subcell to see if any point exceeds the 
      ! threshold density

      ! get the size and centre of the current cell
      cellSize = thisOctal%subcellSize 
      cellCentre = subcellCentre(thisOctal,subCell)
    
      ! calculate the threshold density for the subcell
      criticalValueDouble = REAL(amrLimitScalar,KIND=doubleKind) / cellSize**3.0_oc
      IF ( criticalValueDouble >= HUGE(amrLimitScalar)) THEN
        PRINT *, 'In decideSplit, criticalValue exceeds floating point limit'
        STOP
      END IF
      criticalValue = amrLimitScalar / cellSize**3.0

      ! check the density of random points in the current cell - 
      !   if any of them exceed the critical density, set the flag
      !   indicating the cell is to be split, and return from the function
      split = .FALSE.
      DO i = 1, 400, 1
        searchPoint = cellCentre
        CALL RANDOM_NUMBER(x)
        CALL RANDOM_NUMBER(y)
        CALL RANDOM_NUMBER(z)
        searchPoint%x = searchPoint%x - (cellSize / 2.0_oc) + cellSize*REAL(x,KIND=octalKind) 
        searchPoint%y = searchPoint%y - (cellSize / 2.0_oc) + cellSize*REAL(y,KIND=octalKind) 
        searchPoint%z = searchPoint%z - (cellSize / 2.0_oc) + cellSize*REAL(z,KIND=octalKind) 
!print *, TTauriDensity(searchPoint,grid), criticalValue
        IF ( TTauriDensity(searchPoint,grid) > criticalValue ) THEN
          split = .TRUE.
          RETURN
        END IF
      END DO  

    ELSE
      PRINT *,'In decideSplit, there is no procedure for handling this geometry'
      STOP
    END IF

  END FUNCTION decideSplit


  FUNCTION columnDensity(grid,startPoint,direction,sampleFreq) RESULT(rho)
                 
    ! integrates the density along a line through an amr grid.

    IMPLICIT NONE

    TYPE(gridtype), INTENT(IN)        :: grid
    TYPE(octalVector), INTENT(IN)     :: startPoint
    TYPE(octalVector), INTENT(IN)     :: direction
    REAL(KIND=octalKind), INTENT(IN)  :: sampleFreq

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
                
    CALL startReturnSamples(startPoint,direction,grid,sampleFreq, &
                 nSamples,maxSamples,distances,dummy,dummy, &
                 dummyVel,dummy,dummy,dummyPops,densities,.false.,&
                 hitCore,.false.,1,error)

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



  PURE SUBROUTINE fillGridDummyValues(thisOctal,subcell) 
    ! for testing, just put some generic values in an octal

    IMPLICIT NONE
    
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN)        :: subcell

    thisOctal%kappaAbs(subcell,:) = 1.0e-2 
    thisOctal%kappaSca(subcell,:) = 1.0e-2 
    thisOctal%chiLine(subcell)    = 1.0e-2 
    thisOctal%etaLine(subcell)    = 1.0e-2 
    thisOctal%etaCont(subcell)    = 1.0e-2 
    thisOctal%N(subcell,:)        = 1.0e-2 
    thisOctal%Ne(subcell)         = 1.0e-2 
    thisOctal%nTot(subcell)       = 1.0e-2 
    thisOctal%biasLine3D(subcell) = 1.0 
    thisOctal%biasCont3D(subcell) = 1.0 

  END SUBROUTINE fillGridDummyValues

  
  SUBROUTINE fillVelocityCorners(thisOctal,grid,velocityFunc)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.
    ! this is currently very inefficient, many of the points are identical.
    ! it should probably be changed later.

    IMPLICIT NONE
  
    TYPE(octal), INTENT(INOUT) :: thisOctal
    TYPE(gridtype), INTENT(IN) :: grid

    INTEGER              :: subcell
    TYPE(octalVector)    :: centre
    REAL(KIND=octalKind) :: inc
    
    INTERFACE 
      TYPE(vector) FUNCTION velocityFunc(point,grid)
        USE vector_mod
        USE gridtype_mod
        TYPE(octalVector), INTENT(IN) :: point
        TYPE(gridtype), INTENT(IN)    :: grid
      END FUNCTION velocityFunc
    END INTERFACE

    inc = thisOctal%subcellSize / 2.0
        
    DO subcell = 1, 8, 1
    
      centre = subcellCentre(thisOctal,subcell)

      thisOctal%cornerVelocity(subcell,1) = velocityFunc(centre + &
         octalVector( (-1.0*inc), (-1.0*inc), (-1.0*inc)),grid)
      thisOctal%cornerVelocity(subcell,2) = velocityFunc(centre + &
         octalVector( ( 1.0*inc), (-1.0*inc), (-1.0*inc)),grid)
      thisOctal%cornerVelocity(subcell,3) = velocityFunc(centre + &
         octalVector( (-1.0*inc), ( 1.0*inc), (-1.0*inc)),grid)
      thisOctal%cornerVelocity(subcell,4) = velocityFunc(centre + &
         octalVector( ( 1.0*inc), ( 1.0*inc), (-1.0*inc)),grid)
      thisOctal%cornerVelocity(subcell,5) = velocityFunc(centre + &
         octalVector( (-1.0*inc), (-1.0*inc), ( 1.0*inc)),grid)
      thisOctal%cornerVelocity(subcell,6) = velocityFunc(centre + &
         octalVector( ( 1.0*inc), (-1.0*inc), ( 1.0*inc)),grid)
      thisOctal%cornerVelocity(subcell,7) = velocityFunc(centre + &
         octalVector( (-1.0*inc), ( 1.0*inc), ( 1.0*inc)),grid)
      thisOctal%cornerVelocity(subcell,8) = velocityFunc(centre + &
         octalVector( ( 1.0*inc), ( 1.0*inc), ( 1.0*inc)),grid)
    END DO
      
  END SUBROUTINE fillVelocityCorners


  PURE FUNCTION TTauriDensity(point,grid) RESULT(rho)
    ! calculates the density at a given point for a model of a T Tauri 
    !   star with magnetospheric accretion
    !   see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    use parameters_mod

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    REAL                          :: rho

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec

    REAL :: r, rM, theta, y
    REAL :: rStar
    REAL :: diskRinner
    REAL :: diskRouter

    starPosn = grid%starPos1
    pointVec = (point - starPosn)
    r = modulus( pointVec ) 
    rStar = TTauriRstar / 1.e10
    diskRinner = TTauriRinner / 1.e10
    diskRouter = TTauriRouter / 1.e10

    ! test if the point lies within the star
    IF ( r < rStar ) THEN
      rho = 0.0
      RETURN
    END IF
    
    ! test if the point lies too close to the disk
    IF ( ABS(pointVec%z) < 4.0e-2 * rStar) THEN
      ! we will fake the coordinates so that the density
      !  remains constant within this region 
      pointVec%z = SIGN( 4.0e-2 * rStar, REAL(pointVec%z) )
      r = modulus( pointVec ) 
      RETURN
    END IF
   
    theta = ACOS( pointVec%z  / r )
    IF (ABS(MODULO(theta,pi)) > 1.e-10 ) THEN 
      rM  = r / SIN(theta)**2
    ELSE
      rM = HUGE(rM)
    END IF
     
    ! test if the point lies outside the accretion stream
    IF  ((rM > diskRinner) .AND. (rM < diskRouter )) THEN
      y = SIN(theta)**2 
  
      rho = (TTauriMdot * rStar) / (4.0 * pi * &
              (rStar/diskRinner - rStar/diskRouter))
      rho = rho * r**(-5.0/2.0) / SQRT( 2.0 * bigG * TTauriMstar ) 
      rho = rho * SQRT( 4.0 - 3.0*y) / SQRT( 1.0 - y) 
    ELSE
      rho = 0.0
    END IF
    
  END FUNCTION TTauriDensity

  TYPE(vector) PURE FUNCTION TTauriVelocity(point,grid)
    ! calculates the velocity vector at a given point for a model
    !   of a T Tauri star with magnetospheric accretion
    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    USE parameters_mod

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    TYPE(vector)      :: vP
    REAL              :: modVp
    REAL              :: rStar
    REAL              :: diskRinner
    REAL              :: diskRouter
    REAL              :: phi
    REAL              :: r, rM, theta, y

    starPosn = grid%starPos1
    rStar  = TTauriRstar / 1.e10
    diskRinner = TTauriRinner / 1.e10
    diskRouter = TTauriRouter / 1.e10
 
    pointVec = (point - starPosn)
    r = modulus( pointVec ) 
    theta = ACOS( pointVec%z / r )
    rM  = r / SIN(theta)**2
    y = SIN(theta)**2 

    ! test if the point lies within the star
    IF ( r < rStar ) THEN
      TTauriVelocity = vector(0.0,0.0,0.0)
      
    ! test if the point lies too close to the disk
    ELSE IF ( ABS(pointVec%z) < 4.0e-2 * rStar) THEN
      TTauriVelocity = vector(0.0,0.0,0.0)
   
    ! test if the point lies outside the accretion stream
    ELSE IF ((rM > diskRinner) .AND. (rM < diskRouter )) THEN
  
      vP = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
                  0.0, &
                 (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))
      modVp = SQRT((2.0 * bigG * TTauriMstar / rStar) * &
                     (rStar/r - rStar/rM))
      vP = (-1.0 * (modVp/cSpeed)) * vP
      phi = ATAN2(pointVec%y,pointVec%x)
      vP = rotateZ(vP,phi) 
      IF (theta > pi/2.0) vP%z = -vP%z
      TTauriVelocity = vP

    ELSE
      TTauriVelocity = vector(0.0,0.0,0.0)
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
    REAL :: rStar
    REAL :: diskRinner
    REAL :: diskRouter
    REAL :: phi
    REAL :: r, rM, theta, y

    starPosn = grid%starPos1
    rStar  = TTauriRstar / 1.e10
    diskRinner = TTauriRinner / 1.e10
    diskRouter = TTauriRouter / 1.e10

      point = subcellCentre(thisOctal,subcell)
      pointVec = (point - starPosn)
      r = modulus( pointVec ) 
      theta = ACOS( pointVec%z / r )
      rM  = r / SIN(theta)**2
      y = SIN(theta)**2 

      ! test if the point lies within the star
      IF ( r < rStar ) THEN
        thisOctal%rho(subcell) = 0.0
        thisOctal%velocity(subcell) = vector(0.0,0.0,0.0)
        thisOctal%inFlow(subcell) = .FALSE.
      
      ! test if the point lies too close to the disk
      ELSE IF ( ABS(pointVec%z) < 4.0e-2 * rStar) THEN
        thisOctal%rho(subcell) = 0.0
        thisOctal%velocity(subcell) = vector(0.0,0.0,0.0)
        thisOctal%inFlow(subcell) = .FALSE.
   
      ! test if the point lies outside the accretion stream
      ELSE IF ((rM > diskRinner) .AND. (rM < diskRouter )) THEN
        thisOctal%inFlow(subcell) = .TRUE.
  
        rho = (TTauriMdot * rStar) / (4.0 * pi * &
                 (rStar/diskRinner - rStar/diskRouter))
        rho = rho * r**(-5.0/2.0) / SQRT( 2.0 * bigG * TTauriMstar ) 
        rho = rho * SQRT( 4.0 - 3.0*y) / SQRT( 1.0 - y) 
        TTauriMinRho = MIN(TTauriMinRho, rho/mHydrogen)
        TTauriMaxRho = MAX(TTauriMaxRho, rho/mHydrogen)
        thisOctal%rho(subcell) = rho

        vP = vector(3.0 * SQRT(y) * SQRT(1.0-y) / SQRT(4.0 - (3.0*y)), &
                    0.0, &
                   (2.0 - 3.0 * y) / SQRT(4.0 - 3.0 * y))
        modVp = SQRT((2.0 * bigG * TTauriMstar / rStar) * &
                       (rStar/r - rStar/rM))
        vP = (-1.0 * (modVp/cSpeed)) * vP
        phi = ATAN2(pointVec%y,pointVec%x)
        vP = rotateZ(vP,phi) 
        IF (theta > pi/2.0) vP%z = -vP%z
        thisOctal%velocity(subcell) = vP

      ELSE
        thisOctal%inFlow(subcell) = .FALSE.
        thisOctal%rho(subcell) = 0.0
        thisOctal%velocity(subcell) = vector(0.0,0.0,0.0)
      END IF

      IF (subcell == 8) CALL fillVelocityCorners(thisOctal,grid,TTauriVelocity)

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
        7000.0 - ((2500.0 * rho/mHydrogen - TTauriMinRho) / (TTauriMaxRho/TTauriMinRho)))
    ELSE
      thisOctal%temperature(subcell) = 0.0
    END IF

    ! we will initialise the bias distribution
    thisOctal%biasLine3D(subcell) = 1.0
  
  END SUBROUTINE calcTTauriTemperature


  subroutine initTTauriAMR(grid,Laccretion,Taccretion,&
                           sAccretion,newContFile)

    use constants_mod
    use vector_mod
    use input_variables
    use parameters_mod

    implicit none
    
    type(GRIDTYPE), intent(out)        :: grid                
    real(kind=doublekind), intent(out) :: Laccretion 
    real, intent(out)                  :: Taccretion 
    real, intent(out)                  :: sAccretion
    character(len=80), intent(out)     :: newContFile
   
    real(kind=doublekind) :: TaccretionDouble
    integer :: nNu 
    real :: nuArray(3000),fnu(3000)
    real :: tot
    integer :: i 
    real :: theta1, theta2
    real :: rStar
    real :: diskRinner
    real :: diskRouter
    
    
    rStar  = TTauriRstar / 1.e10
    diskRouter = TTauriRouter /1.e10
    diskRInner = TTauriRinner /1.e10
    
    grid%rCore = rStar
    grid%rStar1 = rStar
    grid%dipoleOffset = dipoleOffset
    grid%diskRadius = diskRInner
    grid%diskNormal = VECTOR(0.,0.,1.)
    grid%diskNormal = rotateX(grid%diskNormal,grid%dipoleOffSet)
    grid%starPos1 = VECTOR(0.,0.,0.)


    theta1 = asin(sqrt(TTauriRstar/TTauriRouter))
    theta2 = asin(sqrt(TTauriRstar/TTauriRinner))


    Laccretion = (REAL(bigG,KIND=db)* &
                  REAL(TTauriMstar,KIND=db)* &
                  REAL(TTauriMdot,KIND=db)/ &
                  REAL(TTauriRstar,kind=db))* &
            REAL((1.0_db-(2.0_db*TTauriRstar/(TTauriRouter+TTauriRinner))),KIND=db)

    TaccretionDouble = Laccretion / REAL(((fourPi * TTauriRstar**2)*stefanBoltz* &
                                      abs(cos(theta1)-cos(theta2))),kind=db)

    sAccretion = (fourPi * TTauriRstar**2)*abs(cos(theta1)-cos(theta2))!/1.e20
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
                (1.e20*sAccretion/(fourPi*TTauriRstar**2))
       write(22,*) nuArray(i), fNu(i)
    enddo
    close(22)
    newContfile="star_plus_acc.dat"


  end subroutine initTTauriAMR


END MODULE amr_mod

