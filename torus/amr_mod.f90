MODULE amr_mod
  ! routines for adaptive mesh refinement. nhs


  USE vector_mod
  USE kind_mod
  USE octal_mod
  USE gridtype_mod            ! type definition for the 3-d grid 
  !use phasematrix_mod         ! needed for photon_mod stuff
  !use utils_mod               ! needed for photon_mod stuff
  !use math_mod                ! needed for photon_mod stuff



  IMPLICIT NONE

  ! we need to track the minimum and maximum densities of the
  !   accretion flow
  REAL, SAVE, PRIVATE :: TTauriMinRho = HUGE(1.0)
  REAL, SAVE, PRIVATE :: TTauriMaxRho = TINY(1.0)


CONTAINS


  SUBROUTINE calcValuesAMR(thisOctal,subcell,grid)
    ! calculates the variables describing one subcell of an octal
  
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal     ! the octal of being considered
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(INOUT) :: grid

    SELECT CASE (grid%geometry)

    CASE ("ttauri")
      CALL calcTTauriMassVelocity(thisOctal,subcell,grid)
      
    CASE DEFAULT
      WRITE(*,*) "! Unrecognised grid geometry: ",trim(grid%geometry)
      STOP

    END SELECT
 
    CALL fillGridDummyValues(thisOctal,subcell)
 
  END SUBROUTINE calcValuesAMR


  SUBROUTINE initFirstOctal(grid, centre, size)
    ! creates the first octal of a new grid (the root of the tree)

    IMPLICIT NONE
    
    TYPE(gridtype), INTENT(INOUT) :: grid 
    TYPE(octalVector), INTENT(IN) :: centre
    REAL(KIND=octalKind), INTENT(IN) :: size 
      ! 'size' should be the vertex length of the cube that contains the whole
      !   of the simulation space, *not* the size of a subcell.
                                         
    
    INTEGER :: subcell ! counter 

    ALLOCATE(grid%octreeRoot)
    grid%octreeRoot%nDepth = 1
    grid%octreeRoot%nChildren = 0
    grid%octreeRoot%hasChild = .FALSE.
    grid%octreeRoot%subcellSize = size/2.0_oc
    grid%octreeRoot%centre = centre
    grid%octreeRoot%indexChild = -999 ! values are undefined
    NULLIFY(grid%octreeRoot%parent) ! because tree root does not have a parent

    ! allocate any variables that need to be 
    ALLOCATE(grid%octreeRoot%kappaAbs(8,grid%maxLevels))
    ALLOCATE(grid%octreeRoot%kappaSca(8,grid%maxLevels))
    ALLOCATE(grid%octreeRoot%N(8,grid%maxLevels))
    
    DO subcell = 1, 8
      ! calculate the values at the centre of each of the subcells
      CALL calcValuesAMR(grid%octreeRoot,subcell,grid)
      ! label the subcells
      grid%octreeRoot%label(subcell) = subcell
    END DO

    grid%maxDepth = 1
    grid%halfSmallestSubcell = grid%octreeRoot%subcellSize / &
                                2.0_oc**REAL(grid%maxDepth,KIND=octalKind)
      
  END SUBROUTINE initFirstOctal


  SUBROUTINE addNewChild(parent, nChild, grid)
    ! adds one child to the parent octal

    IMPLICIT NONE
    
    TYPE(octal), POINTER :: parent     ! pointer to the parent octal 
    INTEGER, INTENT(IN) :: nChild      ! the label (1-8) of the subcell gaining the child 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that calculate
                                          !   the variables stored in the tree
    
    TYPE(octal) :: tempChildStorage    ! holder for existing children, while we make changes
    INTEGER :: subcell
    INTEGER, SAVE :: counter = 9       ! keeps track of the current subcell label
                                       ! this isn't very clever. might change it. 
    INTEGER :: nChildren               ! number of children the parent octal has
    INTEGER :: newChildIndex           ! the storage location for the new child
    INTEGER :: error 


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
      ! array. we need to use a temporary octal structure and copy the
      ! existing children into it; then increase the 'child' array size
      ! by one; then copy the children back in. 

      ALLOCATE(tempChildStorage%child(1:nChildren), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed.'
        STOP
      END IF
      
      tempChildStorage%child(1:nChildren) = parent%child(1:nChildren)
       
      ALLOCATE(parent%child(nChildren+1), STAT=error)
      IF ( error /= 0 ) THEN
        PRINT *, 'Panic: allocation failed.'
        STOP
      END IF

      parent%child(1:nChildren) = tempChildStorage%child(1:nChildren)   
      
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

    ! set up the new child's variables
    parent%child(newChildIndex)%parent => parent
    parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
    parent%child(newChildIndex)%hasChild = .false.
    parent%child(newChildIndex)%nChildren = 0
    parent%child(newChildIndex)%indexChild = -999 ! values are undefined
    parent%child(newChildIndex)%nDepth = parent%nDepth + 1
    parent%child(newChildIndex)%centre = &
                    subcellCentre(parent,nChild)

    ! put some data in the eight subcells of the new child
    DO subcell = 1, 8
      CALL calcValuesAMR(parent%child(newChildIndex),subcell,grid)
      parent%child(newChildIndex)%label(subcell) = counter
      counter = counter + 1
    END DO

        ! VVVVVVVVVV GET RID OF THIS 
        !! if this child's subcells are the smallest ones in the grid,
        !! update the halfSmallestSubcell variable
        !IF ( (parent%child(newChildIndex)%subcellSize / 2.0_oc) &
        !                < halfSmallestSubcell )               &
        !halfSmallestSubcell = (parent%child(newChildIndex)%subcellSize / 2.0_oc)

    ! check for a new maximum depth 
    IF (parent%child(newChildIndex)%nDepth > grid%maxDepth) THEN
      grid%maxDepth = parent%child(newChildIndex)%nDepth
    END IF

  END SUBROUTINE addNewChild


  RECURSIVE SUBROUTINE splitGrid(thisOctal,amrLimitScalar,grid)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    IMPLICIT NONE

    TYPE(OCTAL), POINTER :: thisOctal
    REAL, INTENT(IN) :: amrLimitScalar 
      ! 'limitScalar' is the value the decideSplit function uses to
      !   decide whether or not to split cell.
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid through to the 
                                          !   routines that this subroutine calls
    
    TYPE(OCTAL), POINTER :: childPointer
    
    INTEGER :: subcell, i
    
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
  
  
  RECURSIVE SUBROUTINE finishGrid(thisOctal,grid,gridConverged)
    ! takes the octree grid that has been created according to some
    !   criterion and calculates all the other variables in the model.
    ! this should be called once the grid is complete.
    
    IMPLICIT NONE

    TYPE(octal), POINTER :: thisOctal
    TYPE(gridtype) :: grid
    LOGICAL, INTENT(INOUT) :: gridConverged
    
    TYPE(octal), POINTER :: child
    
    INTEGER :: subcell, i 
     
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

  END SUBROUTINE finishGrid
 
  
  RECURSIVE SUBROUTINE countVoxels(thisVoxel,nOctals,nVoxels)  
    ! count the number of octals in the current section of the grid
    ! also counts the number of unique volume elements (voxels) i.e.
    !   those that are not subdivided

    IMPLICIT NONE

    TYPE(OCTAL), POINTER :: thisVoxel
    INTEGER,INTENT(INOUT) :: nOctals ! number of octals
    INTEGER,INTENT(INOUT) :: nVoxels ! number of childless subcells
    
    TYPE(OCTAL), POINTER :: child

    INTEGER :: i

    nOctals = nOctals + 1

    IF ( thisVoxel%nChildren > 0 ) THEN
      ! call this subroutine recursively on each of its children
      DO i = 1, thisVoxel%nChildren, 1
        child => thisVoxel%child(i)
        CALL countVoxels(child,nOctals,nVoxels)
      END DO
    END IF

    ! increment the counter once for each of its childless subcells
    DO i = (thisVoxel%nchildren+1), 8, 1
        nVoxels = nVoxels + 1
    END DO

  END SUBROUTINE countVoxels


  SUBROUTINE startReturnSamples (startPoint,direction,grid, &
             sampleFreq,nSamples,maxSamples,lambda,kappaAbs,&
             kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,&
             usePops,iLambda,error)
    ! samples the grid at points along the path.
    ! this should be called by the program, instead of calling 
    !   returnSamples directly, because this checks the start and finish
    !   points are within the grid bounds. returnSamples *assumes* this 
    !   criterion is met.
 
    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: startPoint
    TYPE(octalVector), INTENT(IN) :: direction
    TYPE(gridtype), INTENT(IN) :: grid 
    REAL(KIND=octalKind), INTENT(IN) :: sampleFreq
      ! 'sampleFreq' is the maximum number of samples that will be made
      !   in any subcell of the octree. 
    INTEGER, INTENT(INOUT) :: nSamples ! number of samples made
    INTEGER, INTENT(IN) :: maxSamples  ! size of sample arrays 
    REAL, DIMENSION(:), INTENT(INOUT) :: lambda ! distance travelled by photon
    REAL, DIMENSION(:), INTENT(INOUT) :: kappaAbs
    REAL, DIMENSION(:), INTENT(INOUT) :: kappaSca
    TYPE(vector), DIMENSION(:), INTENT(INOUT) :: velocity
    REAL, DIMENSION(:), INTENT(INOUT) :: velocityDeriv
    REAL, DIMENSION(:), INTENT(INOUT) :: chiLine
    REAL, DIMENSION(:,:), INTENT(INOUT) :: levelPop 
    REAL, DIMENSION(:) :: rho
    LOGICAL, INTENT(IN) :: usePops 
    INTEGER, INTENT(IN) :: iLambda 
    INTEGER, INTENT(INOUT) :: error 

    TYPE(octalVector) :: locator 
      ! 'locator' is used to indicate a point that lies within the  
      !   *next* cell of the octree that the ray will interesect.
      !   initially this will be the same as the startPoint
      
    TYPE(octalVector) :: currentPoint  ! current position of ray 
    LOGICAL :: abortRay ! flag to signal completion of ray trace
    TYPE(octal) :: octree 
    TYPE(octalVector) :: directionNormalized

    octree = grid%octreeRoot
    abortRay = .FALSE.
    directionNormalized = direction
    CALL normalize(directionNormalized)
    
    currentPoint = octalVector(startPoint%x, &
                               startPoint%y, &
                               startPoint%z)
    !currentPoint = octalVector(REAL(startPoint%x,KIND=octalKind), &
    !                           REAL(startPoint%y,KIND=octalKind), &
    !                           REAL(startPoint%z,KIND=octalKind))
    
    locator = startPoint


    IF ( inOctal(octree,startPoint) ) THEN
      ! we call the recursive subroutine 'returnSamples'
      CALL returnSamples(currentPoint,startPoint,locator,directionNormalized,octree, &
                   grid,sampleFreq,nSamples,maxSamples,&
                   abortRay,lambda,kappaAbs,kappaSca,velocity,&
                   velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,error)
    ELSE
      PRINT *, 'Attempting to find path between point(s) outwith the grid'
      STOP
    ENDIF

  END SUBROUTINE startReturnSamples

  
  RECURSIVE SUBROUTINE returnSamples (currentPoint,startPoint,locator,direction, &
             octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,kappaAbs, &
             kappaSca,velocity,velocityDeriv,chiLine,levelPop,rho,&
             usePops,iLambda,error)
    ! this uses a recursive ray traversal algorithm to sample the octal
    !   grid at points along the path of the ray. 
    ! no checks are made that the ray lies within the boundaries of the
    !  grid, so this subroutine not be called directly. use  
    !  'startReturnSamples' instead.

    IMPLICIT NONE

    ! note: in the code comments, the terms 'subcell' and 'cell' are
    !  mostly interchangeable.

    TYPE(octalVector), INTENT(INOUT) :: currentPoint ! current ray position
    TYPE(octalVector), INTENT(IN) :: startPoint ! initial ray position
    TYPE(octalVector) :: locator 
      ! 'locator' is used to indicate a point that lies within the  
      !   *next* cell of the octree that the ray will interesect.
      !   initially this will be the same as the currentPoint
    TYPE(octalVector), INTENT(IN) :: direction ! ray's direction vector 
    TYPE(octal), INTENT(IN) :: octree ! an octal grid
    TYPE(gridtype), INTENT(IN) :: grid 
    REAL(KIND=octalKind), INTENT(IN) :: sampleFreq
      ! 'sampleFreq' is the maximum number of samples that will be made
      !   in any subcell of the octree. 
    INTEGER, INTENT(INOUT) :: nSamples ! number of samples made
    INTEGER, INTENT(IN) :: maxSamples  ! size of sample arrays 
    LOGICAL, INTENT(INOUT) :: abortRay ! flag to signal completion
    REAL, DIMENSION(:), INTENT(INOUT) :: lambda ! distance travelled by photon
    REAL, DIMENSION(:), INTENT(INOUT) :: kappaAbs
    REAL, DIMENSION(:), INTENT(INOUT) :: kappaSca
    TYPE(vector), DIMENSION(:), INTENT(INOUT) :: velocity
    REAL, DIMENSION(:), INTENT(INOUT) :: velocityDeriv
    REAL, DIMENSION(:), INTENT(INOUT) :: chiLine
    REAL, DIMENSION(:,:), INTENT(INOUT) :: levelPop 
    REAL, DIMENSION(:) :: rho
    LOGICAL, INTENT(IN) :: usePops 
    INTEGER, INTENT(IN) :: iLambda 
    INTEGER, INTENT(INOUT) :: error 


    REAL(KIND=octalKind) :: halfSmallestSubcell
      ! 'halfSmallestSubcell' is half the vertex length of the grid's smallest
      !    subcells.
    INTEGER :: maxDepth ! maximum depth of octree
    TYPE(octalVector) :: wallNormal ! normal to plane of cell wall
    TYPE(octalVector) :: exitPoint ! where ray leaves current cell
    TYPE(octalVector) :: centre ! centre of current subcell
    REAL(KIND=octalKind) :: wallDistanceX,wallDistanceY,wallDistanceZ
      ! 'wallDistanceX,Y,Z' are the distances between the current 
      !   position and the intersections with the cell walls (along the 
      !   'direction' vector)
    REAL(KIND=octalKind) :: minWallDistance
      ! distance to the nearest cell wall
    
    REAL(KIND=octalKind) :: length 
      ! 'length' is the distance from the start of the ray's path. 
    REAL(KIND=octalKind) :: sampleLength
    REAL(KIND=octalKind) :: trialLength
    TYPE(octalVector) :: trialPoint
    INTEGER :: trial
    REAL(KIND=octalKind) :: wallFromOrigin 
      ! 'wallFromOrigin' is the distance of the cell wall from the origin
    REAL(KIND=octalKind) :: margin
      ! margin is the size of the region around the edge of a subcell
      !   where numerical inaccuracies may cause problems.
    LOGICAL :: found
    INTEGER :: sub, subIndex, i

    halfSmallestSubcell = grid%halfSmallestSubcell
    maxDepth = grid%maxDepth


    DO 
      ! find which of the subcells the point lies in
      sub = whichSubcell(octree,locator)
      IF ( looseInOctal(octree,locator) .EQV. .FALSE. ) THEN
        PRINT *, "looseInOctal failed"
!        PRINT *, "  label = ",octree%label(sub)
        !PRINT *, "  parent labels = ",grid%parent%label
        STOP
      ENDIF
 
      IF ( octree%hasChild(sub) .EQV. .TRUE. ) THEN
        
        ! if the subcell has a child, we will find a pointer to that
        !   child and use recursion to sample it.

        ! find pointer to child
        subindex = -99
        DO i = 1, octree%nChildren, 1
          IF ( octree%indexChild(i) == sub ) THEN
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
                    levelPop,rho,usePops,iLambda,error)
        
        ! after returning from the recursive subroutine, we may have 
        !   finished tracing the ray's path
        IF ( abortRay .EQV. .TRUE. ) RETURN

        ! otherwise, we check whether the ray is still within the
        !   boundaries of the current subsection of the octree.
        IF ( inOctal(octree,locator) .EQV. .FALSE. ) RETURN
      
      
      ELSE
        ! we now consider the case of the current subcell being childless


        centre = subcellCentre(octree,sub)

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
        !  the current subcell. this may not be the case if the ray
        !  entered the subcell very close to a vertex. in that case,
        !  we can try searching the octree for a more appropriate
        !  neighbouring subcell  

        ! if any of the wallDistances are negative, we will jump to
        !  another subroutine that will correct the problem

        IF ( wallDistanceX <= 0.0_oc .OR. &
             wallDistanceY <= 0.0_oc .OR. &
             wallDistanceZ <= 0.0_oc      ) THEN 
           abortRay = .TRUE.
           PRINT *, 'Should have called fixBoundary - negative distance'
           PRINT *, 'wallDistances=',wallDistanceX,wallDistanceY,wallDistanceZ 
           error = -20
           RETURN
        END IF
          
        ! we will abort tracking any rays which are too close to 
        ! to a cell wall
        margin = 8.0_oc * REAL(maxDepth,KIND=octalKind) * EPSILON(1.0_oc)
          
        IF ( wallDistanceX < margin .OR. &
             wallDistanceY < margin .OR. &
             wallDistanceZ < margin ) THEN 
           abortRay = .TRUE.
           PRINT *, 'Should have called fixBoundary, less than margin'
           PRINT *, 'wallDistances=',wallDistanceX,wallDistanceY,wallDistanceZ 
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
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint + (halfSmallestSubcell*xHatOctal)
          ELSE
            wallFromOrigin = centre%x - (octree%subcellSize / 2.0_oc)
            exitPoint = intersectionLinePlane &
                         (currentPoint,direction,wallNormal,wallFromOrigin,found)
            IF ( found .EQV. .FALSE. ) THEN 
              PRINT *, "Panic: wall intersection not found"
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
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
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint + (halfSmallestSubcell*yHatOctal)
          ELSE
            wallFromOrigin = centre%y - (octree%subcellSize / 2.0_oc) 
            exitPoint = intersectionLinePlane &
                         (currentPoint,direction,wallNormal,wallFromOrigin,found)
            IF ( found .EQV. .FALSE. ) THEN 
              PRINT *, "Panic: wall intersection not found"
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
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
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint + (halfSmallestSubcell*zHatOctal)
          ELSE
            wallFromOrigin = centre%z - (octree%subcellSize / 2.0_oc)
            exitPoint = intersectionLinePlane &
                         (currentPoint,direction,wallNormal,wallFromOrigin,found)
            IF ( found .EQV. .FALSE. ) THEN 
              PRINT *, "Panic: wall intersection not found"
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
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
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint + (halfSmallestSubcell*xHatOctal)
          ELSE
            wallFromOrigin = centre%x - (octree%subcellSize / 2.0_oc)
            exitPoint = intersectionLinePlane &
                         (currentPoint,direction,wallNormal,wallFromOrigin,found)
            IF ( found .EQV. .FALSE. ) THEN 
              PRINT *, "Panic: wall intersection not found"
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
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
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint + (halfSmallestSubcell*xHatOctal)
          ELSE
            wallFromOrigin = centre%x - (octree%subcellSize / 2.0_oc)
            exitPoint = intersectionLinePlane &
                         (currentPoint,direction,wallNormal,wallFromOrigin,found)
            IF ( found .EQV. .FALSE. ) THEN 
              PRINT *, "Panic: wall intersection not found"
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint - (halfSmallestSubcell*xHatOctal)
          ENDIF
          IF ( direction%z > 0.0_oc ) THEN             
            locator = locator + (halfSmallestSubcell*zHatOctal)
          ELSE  
            locator = locator - (halfSmallestSubcell*zHatOctal)
          ENDIF  

!PRINT *, '  x=z     locator = ',locator,' exitPoint = ',exitpoint
         ELSEIF ( wallDistanceY == wallDistanceZ  .AND. &
                  wallDistanceX /= wallDistanceZ ) THEN
          wallNormal =  yHatOctal
          IF ( direction%y > 0.0_oc ) THEN
            wallFromOrigin = centre%y + (octree%subcellSize / 2.0_oc)
            exitPoint = intersectionLinePlane &
                         (currentPoint,direction,wallNormal,wallFromOrigin,found)
            IF ( found .EQV. .FALSE. ) THEN 
              PRINT *, "Panic: wall intersection not found"
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint + (halfSmallestSubcell*yHatOctal)
          ELSE
            wallFromOrigin = centre%y - (octree%subcellSize / 2.0_oc)
            exitPoint = intersectionLinePlane &
                         (currentPoint,direction,wallNormal,wallFromOrigin,found)
            IF ( found .EQV. .FALSE. ) THEN 
              PRINT *, "Panic: wall intersection not found"
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint - (halfSmallestSubcell*yHatOctal)
          ENDIF
          IF ( direction%z > 0.0_oc ) THEN             
            locator = locator + (halfSmallestSubcell*zHatOctal)
          ELSE  
            locator = locator - (halfSmallestSubcell*zHatOctal)
          ENDIF  
       
!PRINT *, '  y=z      locator = ',locator,' exitPoint = ',exitpoint
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
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
            END IF
            locator = exitpoint + (halfSmallestSubcell*xHatOctal)
          ELSE
            wallFromOrigin = centre%x - (octree%subcellSize / 2.0_oc)
            exitPoint = intersectionLinePlane &
                         (currentPoint,direction,wallNormal,wallFromOrigin,found)
            IF ( found .EQV. .FALSE. ) THEN 
              PRINT *, "Panic: wall intersection not found"
PRINT *, 'wllDstnce ',wallDistanceX,wallDistanceY,wallDistanceZ 
              DO ! sit in loop for debugging purposes
              END DO
              STOP
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
          PRINT *, 'Fell through direction finder!'
              DO ! sit in loop for debugging purposes
              END DO
          STOP
        END IF

        ! we now decide where we are going to sample the quantities
        ! we will define an approximate rate that is a fraction of
        !   the current subcell size
!print *, 'octree%subcellSize = ', octree%subcellSize       
!print *, 'samplefreq = ', samplefreq       
        sampleLength = octree%subcellSize / sampleFreq

        ! if there are no previously sampled points, we definitely have to take a
        !   sample here
        IF ( nSamples == 0 ) &
          CALL takeSample(currentPoint,length,direction,octree,grid,sub,&
                          nSamples,maxSamples,lambda,kappaAbs,kappaSca,velocity,&
                          velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,error) 

        
        ! we check whether we should take a sample at currentPoint
        !   (which will usually be the entry point to the subcell).

!print *, 'currentPoint = ',currentPoint 
!print *, 'startPoint = ',startPoint 
        length = modulus(currentPoint - startPoint)

        IF ( modulus(currentPoint - (startPoint + (direction *         & 
                              REAL(lambda(nSamples),KIND=octalKind)))) &
                                               > sampleLength ) THEN
          CALL takeSample(currentPoint,length,direction,octree,grid,sub,&
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
          IF ( modulus(trialPoint - currentPoint) < minWallDistance ) THEN 
            trialLength = length + modulus(trialPoint - currentPoint)
            IF ( trialLength >= 0.0_oc ) THEN
              CALL takeSample(trialPoint,trialLength,direction,octree,&
                       grid,sub,nSamples,maxSamples,lambda,kappaAbs,&
                       kappaSca,velocity,velocityDeriv,chiLine,levelPop,&
                       rho,usePops,iLambda,error) 
            ELSE
              EXIT
            END IF
          ELSE 
            EXIT
          END IF
        END DO    

        length = modulus(exitPoint - startPoint)
        currentPoint = exitPoint
        
      ! if we have left the boundaries of the simulation space, we are finished
      IF (inOctal(octree,locator) .EQV. .FALSE.) RETURN 
      
    ENDIF
    
  ENDDO


  END SUBROUTINE returnSamples
   
   
  SUBROUTINE takeSample(point,length,direction,octree,grid,sub,nSamples,&
                        maxSamples,lambda,kappaAbs,kappaSca,velocity,&
                        velocityDeriv,chiLine,levelPop,rho,usePops,iLambda,error) 
  
    
    IMPLICIT NONE
    
    TYPE(octalVector), INTENT(IN) :: point
    REAL(KIND=octalKind), INTENT(IN) :: length 
    TYPE(octalVector), INTENT(IN) :: direction
    TYPE(octal), INTENT(IN) :: octree 
    TYPE(gridtype), INTENT(IN) :: grid 
    INTEGER, INTENT(IN) :: sub 
    INTEGER, INTENT(INOUT) :: nSamples
    INTEGER, INTENT(IN) :: maxSamples  ! size of sample arrays 
    REAL, DIMENSION(:), INTENT(INOUT) :: lambda 
    REAL, DIMENSION(:), INTENT(INOUT) :: kappaAbs
    REAL, DIMENSION(:), INTENT(INOUT) :: kappaSca
    TYPE(vector), DIMENSION(:), INTENT(INOUT) :: velocity
    REAL, DIMENSION(:), INTENT(INOUT) :: velocityDeriv
    REAL, DIMENSION(:), INTENT(INOUT) :: chiLine
    REAL, DIMENSION(:,:), INTENT(INOUT) :: levelPop 
    REAL, DIMENSION(:), INTENT(INOUT) :: rho 
    LOGICAL, INTENT(IN) :: usePops 
    INTEGER, INTENT(IN) :: iLambda 
    INTEGER, INTENT(INOUT) :: error
    
!print *, 'in takeSample, point = ',point
    nSamples = nSamples + 1
    IF (nSamples > maxSamples) THEN
      PRINT *, "nSamples > maxSamples in takeSample subroutine"
      STOP
    END IF
    
    lambda(nSamples) = length
 
!!! Need to add in the routines to interpolate these quantities
 
    kappaAbs(nSamples) = octree%kappaAbs(sub,iLambda) 
    rho(nSamples) = octree%rho(sub) 
    kappaSca(nSamples) = octree%kappaSca(sub,iLambda) 
    velocity(nSamples) = octree%velocity(sub)
    velocityDeriv(nSamples) = 0.001


    IF (usePops) THEN
      levelPop(nSamples,:) = octree%n(sub,:)
    ELSE
      chiLine(nSamples) = octree%chiLine(sub)
    END IF
    

  END SUBROUTINE takeSample


  FUNCTION amrGridVelocity(octalTree,point,startOctal,foundOctal) 
    ! returns the velocity at a given point in the grid.
    ! this function can be called with just the first two arguments.
    !   and it will start at the root of the octal tree to locate
    !   the correct octal.
    ! if the foundOctal argument is supplied, it is made to point to 
    !   the octal containing 'point'.
    ! if the startOctal argument is supplied, the function uses a 
    !   local search for the correct octal starting at that octal.

    IMPLICIT NONE

    TYPE(vector) :: amrGridVelocity
    TYPE(octal), POINTER :: octalTree
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal

    TYPE(octal), POINTER :: resultOctal
    INTEGER :: subcell
    
    IF (PRESENT(startOctal)) THEN
      CALL findSubcellLocal(point,startOctal,subcell)
      IF (PRESENT(foundOctal)) foundOctal => startOctal
      
      ! should add in some interpolation routine
      amrGridVelocity = startOctal%velocity(subcell)
      
    ELSE
      CALL findSubcellTD(point,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal)) foundOctal => resultOctal

      ! should add in some interpolation routine
      amrGridVelocity = resultOctal%velocity(subcell)

    END IF
  
  END FUNCTION amrGridVelocity


  FUNCTION amrGridTemperature(octalTree,point,startOctal,foundOctal) 
    ! returns the temperature at a given point in the grid.
    ! this function can be called with just the first two arguments
    !   and it will start at the root of the octal tree to locate
    !   the correct octal.
    ! if the foundOctal argument is supplied, it is made to point to 
    !   the octal containing 'point'.
    ! if the startOctal argument is supplied, the function uses a 
    !   local search for the correct octal starting at that octal.

    IMPLICIT NONE

    REAL :: amrGridTemperature
    TYPE(octal), POINTER :: octalTree
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(octal), OPTIONAL, POINTER :: startOctal
    TYPE(octal), OPTIONAL, POINTER :: foundOctal

    TYPE(octal), POINTER :: resultOctal
    INTEGER :: subcell
    
    IF (PRESENT(startOctal)) THEN
      CALL findSubcellLocal(point,startOctal,subcell)
      IF (PRESENT(foundOctal)) foundOctal => startOctal
      
      ! should add in some interpolation routine
      amrGridTemperature = resultOctal%temperature(subcell)
      
    ELSE
      CALL findSubcellTD(point,octalTree,resultOctal,subcell)
      IF (PRESENT(foundOctal)) foundOctal => resultOctal

      ! should add in some interpolation routine
      amrGridTemperature = resultOctal%temperature(subcell)

    END IF

  END FUNCTION amrGridTemperature


  RECURSIVE SUBROUTINE locateLineProbAMR(probability,thisOctal,subcell) 
    ! finds the subcell that contains 'probability'.

    IMPLICIT NONE

    REAL(KIND=doubleKind), INTENT(IN) :: probability
    TYPE(octal), POINTER :: thisOctal
    INTEGER, INTENT(OUT) :: subcell

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
    ! finds the subcell that contains 'probability'.

    IMPLICIT NONE

    REAL(KIND=doubleKind), INTENT(IN) :: probability
    TYPE(octal), POINTER :: thisOctal
    INTEGER, INTENT(OUT) :: subcell

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


  FUNCTION whichSubcell(thisOctal,point) RESULT (sub)
    ! returns the identification number (1-8) of the subcell of the 
    ! current octal which contains a given point
    ! NB this does NOT check that the point lies within the bounds of the octal!
  
    IMPLICIT NONE
    
    TYPE(octal), INTENT(IN) :: thisOctal
    TYPE(octalVector), INTENT(IN) :: point
    INTEGER :: sub

    IF ( point%x < thisOctal%centre%x ) THEN
      IF ( point%y < thisOctal%centre%y ) THEN
        IF ( point%z < thisOctal%centre%z ) THEN
          sub = 1
        ELSE 
          sub = 5

        ENDIF
      ELSE 
        IF (point%z < thisOctal%centre%z) THEN
          sub = 3
        ELSE 
          sub = 7
        ENDIF
      END IF
    ELSE
      IF (point%y < thisOctal%centre%y) THEN
        IF (point%z < thisOctal%centre%z) THEN
          sub = 2
        ELSE 
          sub = 6
        ENDIF
      ELSE 
        IF (point%z < thisOctal%centre%z) THEN
          sub = 4
        ELSE 
          sub = 8
        ENDIF
      END IF
    ENDIF

  END FUNCTION whichSubcell    


  FUNCTION inOctal(thisOctal,point) RESULT (found)
    ! true if the point lies in the current octal
  
    IMPLICIT NONE
 
    TYPE(octal), INTENT(IN) :: thisOctal
    TYPE(octalVector), INTENT(IN) :: point
    LOGICAL :: found

    IF ((point%x <= thisOctal%centre%x - thisOctal%subcellSize ) .OR. &
        (point%x >= thisOctal%centre%x + thisOctal%subcellSize ) .OR. &
        (point%y <= thisOctal%centre%y - thisOctal%subcellSize ) .OR. &
        (point%y >= thisOctal%centre%y + thisOctal%subcellSize ) .OR. &
        (point%z <= thisOctal%centre%z - thisOctal%subcellSize ) .OR. &
        (point%z >= thisOctal%centre%z + thisOctal%subcellSize )) THEN
      found = .FALSE.
    ELSE  
      found = .TRUE.
    ENDIF
  
  END FUNCTION inOctal


  FUNCTION looseInOctal(thisOctal,point) RESULT (found)
    ! true if the point lies 'loosely' in the current octal
    ! ( a 10% margin of error is allowed )
    ! this is useful only for testing purposes!
  
    IMPLICIT NONE
 
    TYPE(octal), INTENT(IN) :: thisOctal
    TYPE(octalVector), INTENT(IN) :: point
    LOGICAL :: found

    IF ((point%x <= thisOctal%centre%x - 1.1_oc * thisOctal%subcellSize ) .OR. &
        (point%x >= thisOctal%centre%x + 1.1_oc * thisOctal%subcellSize ) .OR. &
        (point%y <= thisOctal%centre%y - 1.1_oc * thisOctal%subcellSize ) .OR. &
        (point%y >= thisOctal%centre%y + 1.1_oc * thisOctal%subcellSize ) .OR. &
        (point%z <= thisOctal%centre%z - 1.1_oc * thisOctal%subcellSize ) .OR. &
        (point%z >= thisOctal%centre%z + 1.1_oc * thisOctal%subcellSize )) THEN
      found = .FALSE.
    ELSE  
      found = .TRUE.
    ENDIF

  END FUNCTION looseInOctal


  RECURSIVE SUBROUTINE smoothAMRgrid(thisOctal,grid,factor,gridConverged)
    ! checks whether each octal's neighbours are much bigger than it, 
    !   if so, makes the neighbours smaller.
    ! the 'needRestart' flag will be set if a change is made to the octree.
    !   smoothAMRgrid should be called repeatedly until the flag is no 
    !   longer set - this should ensure that the octree eventually stabilizes
  
    IMPLICIT NONE

    TYPE(octal), POINTER :: thisOctal
    TYPE(gridtype), INTENT(INOUT) :: grid 
    REAL(KIND=octalKind), INTENT(IN) :: factor
    LOGICAL, INTENT(INOUT) :: gridConverged
    
    INTEGER :: i
    REAL(KIND=octalKind) :: halfSmallestSubcell
    TYPE(octal), POINTER :: child
    TYPE(octal), POINTER :: neighbour
    TYPE(octalVector), ALLOCATABLE, DIMENSION(:) :: locator
    INTEGER :: subcell
    
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
    
    DO i = 1, 6, 1
      locator(i) = thisOctal%centre
    END DO

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

  RECURSIVE SUBROUTINE findSubcellLocal(point,thisOctal,subcell)
  ! finds the octal (and that octal's subcell) containing a point.
  !   starts searching from the current octal, and goes up and down the
  !   tree as needed to find the correct octal.

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(OUT) :: subcell
    
    TYPE(octal), POINTER :: nextOctal 

    INTEGER :: i

    IF ( inOctal(thisOctal,point) ) THEN
      
      ! if the point lies within the current octal, we identify the
      !   subcell
      subcell = whichSubcell(thisOctal,point)

      ! if the subcell has a child, we look in the child for the point
      IF ( thisOctal%hasChild(subcell) ) THEN
        ! search the index to see where it is stored
        DO i = 1, 8, 1
          IF ( thisOctal%indexChild(i) == subcell ) THEN
            nextOctal => thisOctal%child(i)
            CALL findSubcellLocal(point,nextOctal,subcell)
            EXIT
          END IF
        END DO
      END IF

    ELSE
      ! if the point is outside the current octal, we look in its
      !   parent octal

      ! first check that we are not outside the grid
      IF ( thisOctal%nDepth == 1 ) THEN
        PRINT *, 'Panic: In findSubcellLocal, point is outside the grid'
        STOP
      END IF

      nextOctal => thisOctal%parent
      CALL findSubcellLocal(point,nextOctal,subcell)
      
    END IF    

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
    type(octalVectoR) :: viewVec
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
    do i = 1, 8
       xCube(i) = cubeCorner(i) .dot. xProj
       yCube(i) = cubeCorner(i) .dot. yProj
    enddo

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

    TYPE(octal), POINTER :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    REAL, INTENT(IN) :: amrLimitScalar
    TYPE(gridtype), INTENT(IN) :: grid
    LOGICAL :: split

    REAL :: criticalValue
    REAL(KIND=doubleKind) :: criticalValueDouble
    REAL(KIND=octalKind) :: cellSize
    TYPE(octalVector) :: searchPoint
    TYPE(octalVector) :: cellCentre
    REAL :: x, y, z
    INTEGER :: i
   

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

  SUBROUTINE columnDensity(grid,angle1,angle2,resolution,sampleFreq,densities)

    IMPLICIT NONE

    TYPE(gridtype), INTENT(IN) :: grid
    TYPE(octalVector) :: viewVec
    REAL(KIND=octalKind), INTENT(INOUT) :: angle1, angle2 
    INTEGER, INTENT(IN) :: resolution
    REAL(KIND=octalKind), INTENT(IN) :: sampleFreq
    REAL, DIMENSION(:,:), INTENT(OUT) :: densities

    INTEGER :: u, v, wall
    REAL(KIND=octalKind) :: uReal, vReal

    REAL(KIND=octalKind), DIMENSION(6) :: dist
    REAL(KIND=octalKind), DIMENSION(6) :: Walldist
    REAL(KIND=octalKind) :: WalldistCurrent
    TYPE(octalVector), DIMENSION(6) :: wallNormal
    TYPE(octalVector) :: wallNormalCurrent
    TYPE(octalVector) :: reverseViewVec

    REAL(KIND=octalKind) :: minusOne
    REAL, DIMENSION(:), ALLOCATABLE :: distances, values
    REAL, DIMENSION(:), ALLOCATABLE :: dummy
    REAL, DIMENSION(:,:), ALLOCATABLE :: dummyPops
    TYPE(vector), DIMENSION(:), ALLOCATABLE :: dummyVel
    
    INTEGER :: nSamples, i
    INTEGER :: halfResolution

    TYPE(octalVector) :: startPoint, endPoint
    TYPE(octalVector) :: viewPoint
    LOGICAL :: found
    REAL(KIND=octalKind) :: cubeSize
    REAL :: rho
    LOGICAL :: abandon
    INTEGER :: maxSamples
    INTEGER :: error

    ALLOCATE(distances(maxSamples),values(maxSamples))
    ALLOCATE(dummy(maxSamples),dummyVel(maxSamples))
    ALLOCATE(dummyPops(maxSamples,1))
    cubeSize = grid%octreeRoot%subcellSize * 0.99_oc 
            ! might want to make this slightly smaller?
            
    minusOne = -1.0_oc
    error = 0 

    ! define the planes of the grid walls
    wallNormal(1) = xHatOctal  
    wallNormal(2) = minusOne * xHatOctal  
    wallNormal(3) = yHatOctal  
    wallNormal(4) = minusOne * yHatOctal  
    wallNormal(5) = zHatOctal  
    wallNormal(6) = minusOne * zHatOctal  
    wallDist(1) = 0.99 * cubesize
    wallDist(2) = 0.99 * cubesize
    wallDist(3) = 0.99 * cubesize
    wallDist(4) = 0.99 * cubesize
    wallDist(5) = 0.99 * cubesize
    wallDist(6) = 0.99 * cubesize

    halfResolution = resolution / 2  -1


    DO u = -halfResolution, halfResolution,  1
      uReal = 1.0_oc * REAL(u,KIND=octalKind) * &
                       cubeSize / REAL(halfResolution,KIND=octalKind) +&
                       0.01_oc * cubeSize / REAL(halfResolution,KIND=octalKind)
                       
      DO v = -halfResolution, halfResolution,  1
      vReal = 1.0_oc * REAL(v,KIND=octalKind) * &
                       cubeSize / REAL(halfResolution,KIND=octalKind) + &  
                       0.01_oc * cubeSize / REAL(halfResolution,KIND=octalKind)
      abandon = .FALSE.

      ! we start tracing the ray at the point where it hits the nearest
      !   wall
      viewVec = octalVector(0.0,1.0,0.0)
      viewPoint = octalVector(uReal,0.0_oc,vReal)
      viewPoint = rotateX(viewPoint,angle1)
      viewPoint = rotateZ(viewPoint,angle2)
      viewVec = rotateX(viewVec,angle1)
      viewVec = rotateZ(viewVec,angle2)
     
      DO wall = 1, 6
      
        dist(wall) = modulus( viewPoint - intersectionLinePlane( &
              viewPoint,viewVec,wallNormal(wall),wallDist(wall),found))

        IF ( found .EQV. .FALSE. ) dist(wall) = HUGE(1.0)
        IF ( dist(wall) < 0.0_oc ) dist(wall) = HUGE(1.0) 

      END DO

      IF ( MINVAL(dist) > (2.0_oc * cubeSize)) THEN
        abandon = .TRUE.
      END IF

      wallNormalCurrent = wallNormal(TRANSFER(MINLOC(dist),1))
      wallDistCurrent =  wallDist(TRANSFER(MINLOC(dist),1))
      startPoint = intersectionLinePlane( viewPoint,viewVec, &
                wallNormalCurrent,wallDistCurrent,found)
      
      ! we stop tracing the ray at the point where it hits the nearest
      !   wall in the negative direction
      viewPoint = octalVector(uReal,0.0_oc,vReal)
      viewPoint = rotateX(viewPoint,angle1)
      viewPoint = rotateZ(viewPoint,angle2)
      reverseViewVec = minusOne * viewVec 
       
      DO wall = 1, 6
      
        dist(wall) = modulus( viewPoint - intersectionLinePlane( &
              viewPoint,reverseViewVec,wallNormal(wall),wallDist(wall),found))

        IF ( found .EQV. .FALSE. ) dist(wall) = HUGE(1.0)
        IF ( dist(wall) < 0.0_oc ) dist(wall) = HUGE(1.0) 

      END DO

      IF ( MINVAL(dist) > (2.0_oc * cubeSize)) THEN
        abandon = .TRUE.
      END IF

      wallNormalCurrent = wallNormal(TRANSFER(MINLOC(dist),1))
      wallDistCurrent =  wallDist(TRANSFER(MINLOC(dist),1))
      endPoint = intersectionLinePlane( viewPoint,reverseViewVec, &
                wallNormalCurrent,wallDistCurrent,found)
    
      IF (( .NOT. abandon ) .AND. ( inOctal(grid%octreeRoot,startPoint) ) &
           .AND. ( inOctal(grid%octreeRoot,endPoint) )) THEN
      
        
        nSamples = 0
        
        CALL startReturnSamples(startPoint,viewVec,grid,sampleFreq, &
                     nSamples,maxSamples,distances,dummy,dummy, &
                     dummyVel,dummy,dummy,dummyPops,values,.false.,1,error)
 
        IF (nSamples == 1) THEN  
          rho = 0.0
        ELSE
          rho = 0.0
          DO i = 2, nSamples, 1
             rho = rho + ( (values(i) + values(i-1)) / 2.0 ) * &
                        (distances(i-1) - distances(i)  )
          END DO
        END IF

      ELSE
        rho = 0.0
      END IF
    densities(u+halfResolution+1,v+halfResolution+1) = rho

    END DO
  END DO

  END SUBROUTINE columnDensity




  SUBROUTINE fillGridDummyValues(thisOctal,subcell) 
    ! for testing, just put some generic values in an octal

    IMPLICIT NONE
    
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell

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


  FUNCTION TTauriDensity(point,grid) RESULT(rho)
    ! calculates the density at a given point for a model of a T Tauri 
    !   star with magnetospheric accretion
    !   see Hartman, Hewett & Calvet 1994ApJ...426..669H 

    use parameters_mod

    IMPLICIT NONE

    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN) :: grid
    REAL :: rho

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec

    REAL :: r, rM, theta, y
    REAL :: rStar
    REAL :: diskRinner
    REAL :: diskRouter

    starPosn = octalVector(grid%starPos1%x, &
                           grid%starPos1%y, &
                           grid%starPos1%z)

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
    rM  = r / SIN(theta)**2
     
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
    
    starPosn = octalVector(grid%starPos1%x, &
                           grid%starPos1%y, &
                           grid%starPos1%z)

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
        phi = atan2(pointVec%y,pointVec%z)
        vP = rotateZ(vP,phi) 
        IF ( theta > pi/2.0 ) vP%z = -vP%z
        thisOctal%velocity(subcell) = vP

      ELSE
        thisOctal%inFlow(subcell) = .FALSE.
        thisOctal%rho(subcell) = 0.0
        thisOctal%velocity(subcell) = vector(0.0,0.0,0.0)
      END IF


  END SUBROUTINE calcTTauriMassVelocity
    
  SUBROUTINE calcTTauriTemperature(thisOctal,subcell) 
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
    
    type(GRIDTYPE)    :: grid                
    real(kind=doublekind), intent(out) :: Laccretion 
    real, intent(out) :: Taccretion 
    real, intent(out) :: sAccretion
    character(len=80), intent(out) :: newContFile
   
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
    grid%diskRadius = diskRouter
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
print *, 'Laccretion = ',laccretion    

    TaccretionDouble = Laccretion / REAL(((fourPi * TTauriRstar**2)*stefanBoltz* &
                                      abs(cos(theta1)-cos(theta2))),kind=db)

    sAccretion = (fourPi * TTauriRstar**2)*abs(cos(theta1)-cos(theta2))!/1.e20
print *, 'Saccretion = ',Saccretion    
    Taccretion = TaccretionDouble**0.25
print *, 'Taccretion = ',Taccretion    

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

