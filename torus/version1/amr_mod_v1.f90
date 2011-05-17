!!$  SUBROUTINE startReturnSamples (startPoint,direction,grid,          &
!!$             sampleFreq,nSamples,maxSamples,thin_disc_on, opaqueCore,hitCore,      &
!!$             usePops,iLambda,error,lambda,kappaAbs,kappaSca,velocity,&
!!$             velocityDeriv,chiLine,levelPop,rho, temperature, Ne, inflow, &
!!$             etaCont, etaLine)
!!$    ! samples the grid at points along the path.
!!$    ! this should be called by the program, instead of calling 
!!$    !   returnSamples directly, because this checks the start and finish
!!$    !   points are within the grid bounds. returnSamples *assumes* this 
!!$    !   criterion is met. also, the path of the photon is checked for 
!!$    !   intersection with the stellar surface(s), disc etc. 
!!$ 
!!$    IMPLICIT NONE
!!$
!!$    TYPE(vector), INTENT(IN)      :: startPoint ! photon start point
!!$    TYPE(vector), INTENT(IN)      :: direction  ! photon direction 
!!$    TYPE(gridtype), INTENT(IN)         :: grid       ! the entire grid structure
!!$    REAL, INTENT(IN)                   :: sampleFreq ! the maximum number of
!!$!    real(oct), INTENT(IN)   :: sampleFreq ! the maximum number of 
!!$                       ! samples that will be made in any subcell of the octree. 
!!$                       
!!$    INTEGER, INTENT(OUT)             :: nSamples   ! number of samples made
!!$    INTEGER, INTENT(IN)                :: maxSamples ! size of sample arrays 
!!$    logical, intent(in)                :: thin_disc_on   ! T to include thin disc
!!$    LOGICAL, INTENT(IN)                :: opaqueCore ! true if the core is opaque
!!$    LOGICAL, INTENT(OUT)               :: hitCore    ! true if the core is opaque
!!$    LOGICAL, INTENT(IN)                :: usePops    ! whether to use level populations
!!$    INTEGER, INTENT(IN)                :: iLambda    ! wavelength index
!!$    INTEGER, INTENT(INOUT)             :: error      ! error code
!!$    REAL, DIMENSION(:), INTENT(INOUT)  :: lambda     ! path distances of sample locations
!!$
!!$    
!!$    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaAbs   ! continuous absorption opacities
!!$    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: kappaSca   ! scattering opacities
!!$    TYPE(vector),DIMENSION(:),INTENT(INOUT),OPTIONAL :: velocity ! sampled velocities
!!$    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: velocityDeriv ! sampled velocity derivatives
!!$    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: chiLine    ! line opacities
!!$    REAL(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: rho        ! density at sample points
!!$    real(double),DIMENSION(:,:),INTENT(INOUT),OPTIONAL:: levelPop ! level populations
!!$    REAL,DIMENSION(:),INTENT(INOUT),OPTIONAL  :: temperature! temperature [K] at sample points
!!$    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: Ne ! electron density
!!$    logical,DIMENSION(:),INTENT(INOUT),OPTIONAL  ::inFlow   ! flag to tell if the point is inflow
!!$    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaCont ! contiuum emissivity
!!$    real(double),DIMENSION(:),INTENT(INOUT),OPTIONAL  :: etaLine ! line emissivity
!!$
!!$    TYPE(vector)       :: locator 
!!$                       ! 'locator' is used to indicate a point that lies within the  
!!$                       !   *next* cell of the octree that the ray will interesect.
!!$                       !   initially this will be the same as the startPoint
!!$      
!!$    TYPE(vector)       :: currentPoint ! current position of ray 
!!$    LOGICAL                 :: abortRay     ! flag to signal completion of ray trace
!!$    TYPE(octal)             :: octree       ! the octree structure within 'grid'
!!$    TYPE(vector)       :: directionNormalized
!!$    real(oct)    :: distanceLimit ! max length of ray before aborting
!!$    ! margin is the size of the region around the edge of a subcell
!!$    !   where numerical inaccuracies may cause problems.
!!$    real(oct)    :: margin
!!$ 
!!$    ! variables for testing special cases (stellar intersections etc.)
!!$    TYPE(vector)       :: starPosition       ! position vector of stellar centre
!!$    TYPE(vector)       :: diskNormal         ! disk normal vector
!!$    real(oct)    :: rStar              ! stellar radius
!!$    real(oct)    :: endLength          ! max path length of photon
!!$    TYPE(vector)       :: endPoint           ! where photon leaves grid
!!$    LOGICAL                 :: absorbPhoton       ! photon will be absorbed
!!$    real(oct)    :: distanceFromOrigin ! closest distance of plane from origin
!!$    TYPE(vector)       :: diskIntersection   ! point of photon intersection with disk
!!$    real(oct)    :: starIntersectionDistance1 ! distance to first  intersection
!!$    real(oct)    :: starIntersectionDistance2 ! distance to second intersection
!!$    LOGICAL                 :: starIntersectionFound1 ! 1st star intersection takes place      
!!$    LOGICAL                 :: starIntersectionFound2 ! 2nd star intersection takes place      
!!$    LOGICAL                 :: intersectionFound  ! true when intersection takes place      
!!$    real(oct)    :: intersectionRadius ! disk radius when photon intersects
!!$    real(oct)    :: diskDistance       ! distance to disk intersection
!!$    real(oct)    :: distanceThroughStar! distance of chord through star
!!$    TYPE(vector)       :: dummyStartPoint    ! modified start point 
!!$!    real(oct), PARAMETER :: fudgefactor = 1.00001 ! overestimates stellar size
!!$    real(oct), PARAMETER :: fudgefactor = 1.000001 ! overestimates stellar size
!!$    
!!$    
!!$    ! we will abort tracking a photon just before it reaches the edge of the
!!$    !   simulation space. This is the fraction of the total distance to use:
!!$    real(oct), PARAMETER :: distanceFraction = 0.999_oc 
!!$    
!!$    ! we will abort tracking any rays which are too close to 
!!$    !   to a cell wall. The distance to use is defined by:
!!$    margin = 6.0_oc * REAL(grid%maxDepth,kind=oct) * EPSILON(1.0_oc)
!!$    ! some more experimentation is required to find out the best value for
!!$    !   margin.
!!$
!!$    
!!$    ! set up some variables
!!$    octree = grid%octreeRoot  
!!$    abortRay = .FALSE.
!!$    hitCore = .FALSE.
!!$    directionNormalized = direction
!!$    CALL normalize(directionNormalized)
!!$    distanceLimit = HUGE(distanceLimit)
!!$    locator = startPoint
!!$    nSamples = 0
!!$
!!$    currentPoint = startPoint
!!$    
!!$    IF (.NOT. inOctal(octree,startPoint)) THEN
!!$      PRINT *, 'Attempting to find path between point(s) outwith the grid.'
!!$      PRINT *, ' in amr_mod::startReturnSamples.'
!!$      PRINT *, ' ==> StartPoint = (', startPoint, ')'
!!$      error = -30
!!$      return
!!$!      startpoint%x = min(max(startpoint%x,octree%centre%x - octree%subcellSize),octree%centre%x + octree%subcellSize)
!!$!      startpoint%y = min(max(startpoint%y,octree%centre%y - octree%subcellSize),octree%centre%y + octree%subcellSize)
!!$!      startpoint%z = min(max(startpoint%z,octree%centre%z - octree%subcellSize),octree%centre%z + octree%subcellSize)
!!$
!!$    ENDIF
!!$   
!!$   
!!$    ! geometry-specific tests should go here
!!$    IF (grid%geometry(1:6) == "ttauri" .OR. grid%geometry(1:4) == "jets" .or. &
!!$         grid%geometry(1:9) == "luc_cir3d" .or. grid%geometry(1:6) == "cmfgen" .or. &
!!$         grid%geometry(1:8) == "romanova") THEN
!!$      
!!$       ! need to test for both star and disc intersections 
!!$      
!!$       ! we will find out when and where the photon leaves the simulation space 
!!$
!!$       call getExitPoint2(currentPoint, directionNormalized ,locator, abortRay, error, &
!!$       grid%halfSmallestSubcell, endPoint, endLength, grid, octree, edgeofgrid=.true.)
!!$
!!$
!!$!       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
!!$!                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
!!$!                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 
!!$       endLength = endLength * distanceFraction
!!$       
!!$       ! need to reset some of the variables
!!$       abortRay = .FALSE.
!!$       locator = startPoint
!!$       
!!$       ! for the moment, we will just assume a 2-d disc.
!!$       absorbPhoton = .FALSE.
!!$       diskNormal = grid%diskNormal
!!$       starPosition = grid%starPos1
!!$                                  
!!$       ! find the (geometrycally thin) disk intersection point
!!$       if (grid%geometry == "luc_cir3d" .or. grid%geometry == "cmfgen") then
!!$          intersectionFound = .false.
!!$       else
!!$          if (thin_disc_on) then
!!$             distanceFromOrigin = modulus(grid%starPos1)
!!$             diskIntersection = intersectionLinePlane(startPoint, directionNormalized,&
!!$                  diskNormal, distanceFromOrigin, intersectionFound)
!!$          else
!!$             intersectionFound = .false.
!!$          end if
!!$       end if
!!$
!!$       
!!$       IF (intersectionFound) THEN
!!$       
!!$         ! we need to check whether the photon passes through the disk's
!!$         !   central hole, or the outside of the outer radius of accretion disc.
!!$         intersectionRadius =  modulus(diskIntersection - starPosition)
!!$         IF (intersectionRadius > grid%diskRadius .and.  &
!!$              intersectionRadius < 1.5d5) THEN  ! assuming 100 AU disc size
!!$!              intersectionRadius < grid%octreeRoot%subcellsize*2.0) THEN
!!$           absorbPhoton = .TRUE.
!!$
!!$           ! we need to check whether the intersection occurs within our
!!$           !   simulation space.
!!$           diskDistance = modulus(diskIntersection-startPoint)
!!$           IF (diskDistance < endLength) THEN
!!$           
!!$             endLength = diskDistance
!!$             
!!$           END IF
!!$         END IF
!!$       END IF
!!$         
!!$         
!!$       ! now we check for intersections with the star
!!$
!!$       rStar = grid%rStar1 
!!$       CALL intersectionLineSphere(startPoint,directionNormalized,endLength,starPosition, &
!!$                                   rStar,starIntersectionFound1,starIntersectionFound2,   &
!!$                                   starIntersectionDistance1,starIntersectionDistance2)
!!$       ! by passing a line segment to intersectionLineSphere, we ensure that we
!!$       !   do not find intersections with the star that take place after the 
!!$       !   photon has been absorbed by the disk.
!!$       IF (starIntersectionFound1) THEN
!!$       
!!$         endLength = starIntersectionDistance1
!!$         IF (opaqueCore) absorbPhoton = .TRUE.
!!$         
!!$       END IF
!!$
!!$       ! we trace the photon path until it encounters the disk, the star,
!!$       !   or the edge of the simulation space.
!!$       error = 0
!!$       CALL returnSamples(currentPoint,startPoint,locator,directionNormalized,&
!!$                     octree,grid,sampleFreq,nSamples,maxSamples,abortRay,     &
!!$                     lambda,usePops,iLambda,error,margin,endLength,           &
!!$                     kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,   &
!!$                     velocityDeriv=velocityDeriv,chiLine=chiLine,             &
!!$                     levelPop=levelPop,rho=rho, temperature=temperature,      &
!!$                     Ne=Ne, inFlow=inFlow, etaCont=etaCont, etaLine=etaLine)
!!$      
!!$       IF (error < 0)  RETURN
!!$
!!$       IF (.NOT. opaqueCore .AND. starIntersectionFound2) THEN 
!!$         ! the photon passes through the star and continues on the other side.
!!$
!!$         ! we have to adjust the arguments to returnSamples to make the 
!!$         !   output arrays correct.
!!$         distanceThroughStar = fudgeFactor * &
!!$                      (starIntersectionDistance2 - starIntersectionDistance1)
!!$         currentPoint = currentPoint + distanceThroughStar * directionNormalized
!!$         locator = currentPoint
!!$         dummystartPoint = startPoint + distanceThroughStar * directionNormalized
!!$         distanceLimit = endLength - distanceThroughStar
!!$           
!!$        CALL returnSamples(currentPoint,dummyStartPoint,locator,             &
!!$                    directionNormalized,octree,grid,sampleFreq,nSamples,     &
!!$                    maxSamples,abortRay,lambda,usePops,iLambda,error,margin, &
!!$                    distanceLimit,kappaAbs=kappaAbs,kappaSca=kappaSca,       &
!!$                    velocity=velocity,velocityDeriv=velocityDeriv,           &
!!$                    chiLine=chiLine,levelPop=levelPop,rho=rho,               &
!!$                    temperature=temperature, Ne=Ne, inFlow=inFlow,           &
!!$                    etaCont=etaCont, etaLine=etaLine)
!!$
!!$       END IF
!!$
!!$       ! if the photon ends up in the disk or the star, we make sure it absorbed.
!!$       IF (absorbPhoton) THEN
!!$       
!!$         nSamples = nSamples + 1
!!$         IF (nSamples > maxSamples) THEN
!!$           PRINT *, "Error:: nSamples > maxSamples in startReturnSamples subroutine"
!!$           PRINT *, "        nSamples   = ", nSamples
!!$           PRINT *, "        maxSamples = ", maxSamples
!!$           STOP
!!$         END IF
!!$         
!!$         lambda(nSamples) = lambda(nSamples-1)
!!$         hitCore = .TRUE.
!!$         IF (PRESENT(kappaSca))      kappaSca(nSamples) = 0.
!!$         IF (PRESENT(kappaAbs))      kappaAbs(nSamples) = 1.e20
!!$         IF (PRESENT(velocity))      velocity(nSamples) = vector(0.,0.,0.)
!!$         IF (PRESENT(velocityDeriv)) velocityDeriv(nSamples) = 1.0
!!$         IF (PRESENT(chiLine))       chiLine(nSamples) = 1.e20
!!$         IF (PRESENT(levelPop))      levelPop(nSamples,:) = 0.0
!!$         IF (PRESENT(Ne))            Ne(nSamples) = 0.0
!!$         IF (PRESENT(rho))           rho(nSamples) = 0.0
!!$         IF (PRESENT(temperature))   temperature(nSamples) = 0.0
!!$         IF (PRESENT(inFlow))        inFlow(nSamples) = .true.
!!$         IF (PRESENT(etaCont))       etaCont(nSamples) = 1.e-20
!!$         IF (PRESENT(etaLine))       etaLine(nSamples) = 1.e-20
!!$       END IF
!!$       
!!$    ELSE
!!$            
!!$
!!$       call getExitPoint2(currentPoint, directionNormalized ,locator, abortRay, error, &
!!$       grid%halfSmallestSubcell, endPoint, endLength, grid, octree, edgeofgrid=.true.)
!!$
!!$
!!$!       CALL getExitPoint(currentPoint,directionNormalized,locator,abortRay,error,&
!!$!                         grid%halfSmallestSubcell,endPoint,octree%centre,        &
!!$!                         (octree%subcellSize*2.0_oc),endLength,margin,grid%octreeRoot%threed) 
!!$       distanceLimit = endLength * distanceFraction
!!$       
!!$       ! need to reset some of the variables
!!$       abortRay = .FALSE.
!!$       locator = startPoint
!!$       
!!$       CALL returnSamples(currentPoint,startPoint,locator,directionNormalized,&
!!$                   octree,grid,sampleFreq,nSamples,maxSamples,abortRay,lambda,&
!!$                   usePops,iLambda,error,margin,distanceLimit,                &
!!$                   kappaAbs=kappaAbs,kappaSca=kappaSca,velocity=velocity,     &
!!$                   velocityDeriv=velocityDeriv,chiLine=chiLine,               &
!!$                   levelPop=levelPop,rho=rho,temperature=temperature,         &
!!$                   Ne=Ne,inFlow=inFlow,etaCont=etaCont,etaLine=etaLine)
!!$    END IF
!!$      
!!$  END SUBROUTINE startReturnSamples


!!$  subroutine getExitPoint2(currentPoint, direction ,locator, abortRay, error, &
!!$       halfSmallestSubcell, exitPoint, minWallDistance, grid, thisOctal, edgeOfGrid)
!!$
!!$
!!$    implicit none
!!$    type(GRIDTYPE) :: grid
!!$    type(VECTOR), intent(in) :: currentPoint
!!$    type(VECTOR), intent(in) :: direction
!!$    type(VECTOR), intent(out) :: locator    
!!$    logical, intent(inout) :: abortRay
!!$    type(OCTAL), target  :: thisOctal
!!$    integer, intent(inout) :: error
!!$    real(oct), INTENT(IN)    :: halfSmallestSubcell
!!$    type(VECTOR), intent(out) :: exitPoint
!!$    real(oct), intent(out) :: minWallDistance
!!$    type(OCTAL), pointer :: sOctal
!!$    logical, optional :: edgeOfGrid
!!$    
!!$    error = 0
!!$    abortRay = .false.
!!$
!!$    if (present(edgeofGrid)) then
!!$       call distanceToGridEdge(grid, currentPoint, direction, minWallDistance)
!!$    else
!!$       sOctal => thisOctal
!!$       call distanceToCellBoundary(grid, currentPoint, direction, minWallDistance, sOctal)
!!$    endif
!!$
!!$    exitPoint = currentPoint + minWallDistance * direction
!!$    locator = currentPoint + (minWallDistance + 0.001d0*halfSmallestSubcell) * direction
!!$
!!$  end subroutine getExitPoint2

!!$  PURE SUBROUTINE fillGridDummyValues(thisOctal,subcell, grid) 
!!$    ! for testing, just put some generic values in an octal
!!$
!!$    IMPLICIT NONE
!!$    
!!$    TYPE(octal), INTENT(INOUT) :: thisOctal
!!$    INTEGER, INTENT(IN)        :: subcell
!!$    type(gridtype), intent(in) :: grid
!!$
!!$    if (.not.grid%oneKappa) then
!!$       thisOctal%kappaAbs(subcell,:) = 1.0e-20
!!$       thisOctal%kappaSca(subcell,:) = 1.0e-20 
!!$    endif
!!$    thisOctal%chiLine(subcell)    = 1.0e-20
!!$    thisOctal%etaLine(subcell)    = 1.0e-20 
!!$    thisOctal%etaCont(subcell)    = 1.0e-20 
!!$    thisOctal%N(subcell,:)        = 1.0e-20 
!!$    thisOctal%Ne(subcell)         = 1.0e-20 
!!$    thisOctal%nTot(subcell)       = 1.0e-20
!!$    thisOctal%biasLine3D(subcell) = 1.0 
!!$    thisOctal%biasCont3D(subcell) = 1.0 
!!$
!!$  END SUBROUTINE fillGridDummyValues


  !
  !  Given a octal object, this routine will "reassign" the values of density at the outer
  !  edge of the T Taur accretion stream. THIS ROUTINE SHOULD BE CALLED AS A POST PROCESS (AFTER THE 
  !  NORMAL DENSITY OF ACCRETION FLOW HAS BEEN ASSIGNED). It has to be done this way since 
  !  if we were to assined very small density, the cells at the edge will not split (because they
  !  are split by the total mass, c.f. decideSplit routine in this module).
  !  
!!$  RECURSIVE SUBROUTINE TTauri_fuzzy_edge(thisOctal)
!!$    use input_variables, only: TTauriRinner, TTauriRouter
!!$
!!$    IMPLICIT NONE
!!$    TYPE(OCTAL), POINTER :: thisOctal    
!!$    !
!!$    TYPE(OCTAL), POINTER :: childPointer  
!!$    INTEGER              :: subcell, i    ! loop counters
!!$    TYPE(vector)     :: rvec
!!$    integer :: n
!!$    real(double) :: r, theta, rM, w, p, rho
!!$    real(double) :: rM_fuzzy_in, rM_fuzzy_out  ! beginning of the fuzzy edges
!!$    !
!!$    real(double), parameter :: scale = 7.0d0  ! a scale in exponential decay of density.
!!$!    real(double), parameter :: scale = 10.0d0  ! a scale in exponential decay of density.
!!$    !
!!$    
!!$    if (thisOctal%threeD) then
!!$       n = 8
!!$    else
!!$       n = 4
!!$    endif
!!$
!!$    do subcell = 1, n    
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the index of the new child and decend the tree
!!$          do i = 1, n
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                childPointer => thisOctal%child(i)
!!$                CALL TTauri_fuzzy_edge(childPointer)
!!$                exit
!!$             end if
!!$          end do
!!$       else 
!!$          ! This is a leaf, so modefiy the density value if required.
!!$          if (thisOctal%inFlow(subcell)) then
!!$             ! get the size and the position of the centre of the current cell
!!$             rVec = subcellCentre(thisOctal,subcell)
!!$             r = modulus(rvec)
!!$             
!!$             if (r/=0.0d0) then
!!$                theta = ACOS( MIN(ABS(rVec%z/r),0.998_oc) )
!!$             else
!!$                theta=0.01
!!$             end if
!!$             
!!$             if (theta == 0.0d0) theta = 0.01
!!$             rM  = r / SIN(theta)**2
!!$                          
!!$             ! The fuzzy density starts from a 5-th of the thickness (2h) 
!!$             ! below the surface.
!!$             w = get_fuzzy_width() ! [10^10cm] using a fucntion in this module
!!$
!!$             rM_fuzzy_in  = TTauriRinner*1.0d-10 + w   ! [10^10cm]
!!$             rM_fuzzy_out = TTauriRouter*1.0d-10 - w   ! [10^10cm]
!!$             
!!$             ! If the point is close to the edge, we make it fuzzy
!!$             if ( rM < rM_fuzzy_in) then
!!$!                thisOctal%rho(subcell) =  1.0d-13  ! limit the minimum value...
!!$                p = (rM_fuzzy_in-rM)
!!$                rho = thisOctal%rho(subcell)* EXP(-scale*p/w)
!!$                thisOctal%rho(subcell) =  MAX(rho, 1d-25)  ! limit the minimum value...
!!$             elseif (rM > rM_fuzzy_out) then
!!$                p = (rM - rM_fuzzy_out)
!!$                rho = thisOctal%rho(subcell)* EXP(-scale*p/w)
!!$                thisOctal%rho(subcell) =  MAX(rho, 1d-25)  ! limit the minimum value...
!!$             else
!!$                ! don't touch the density, and just continue
!!$                continue
!!$             end if
!!$
!!$          end if  ! if inflow
!!$
!!$       end IF
!!$
!!$    end do
!!$
!!$  END SUBROUTINE TTauri_fuzzy_edge
!!$

  ! Scale the density in the accretion flow by the scale factor passed on this routine.
  ! Before the initial call the total mass must be set to zero
!!$  RECURSIVE SUBROUTINE TTauri_accretion_scale_density(thisOctal, grid, scale)
!!$
!!$    IMPLICIT NONE
!!$    TYPE(OCTAL), POINTER :: thisOctal    
!!$    TYPE(gridtype), INTENT(INOUT)    :: grid
!!$    real(double), intent(in) :: scale
!!$    !
!!$    TYPE(OCTAL), POINTER :: childPointer  
!!$    INTEGER              :: subcell, i    ! loop counters
!!$    TYPE(vector)     :: point
!!$    integer :: n
!!$    real(double) :: rho
!!$
!!$    !
!!$    
!!$    if (thisOctal%threeD) then
!!$       n = 8
!!$    else
!!$       n = 4
!!$    endif
!!$
!!$    do subcell = 1, n    
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the index of the new child and decend the tree
!!$          do i = 1, n
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                childPointer => thisOctal%child(i)
!!$                CALL TTauri_accretion_scale_density(childPointer, grid, scale)
!!$                exit
!!$             end if
!!$          end do
!!$       else 
!!$          point = subcellCentre(thisOctal,subcell)
!!$          ! This is a leaf, so modefiy the density value if required.
!!$          if (thisOctal%inFlow(subcell) .and. TTauriInFlow(point,grid)) then
!!$             rho = thisOctal%rho(subcell) * scale
!!$             thisOctal%rho(subcell) = rho
!!$          end if  ! if inflow
!!$       end IF
!!$
!!$    end do
!!$
!!$  END SUBROUTINE TTauri_accretion_scale_density

!!$  PURE SUBROUTINE calcTTauriTemperature(thisOctal,subcell) 
!!$    ! calculates the temperature in an octal for a model
!!$    !   of a T Tauri star with magnetospheric accretion
!!$    ! see Hartman, Hewett & Calvet 1994ApJ...426..669H 
!!$
!!$    USE input_variables, ONLY: isoTherm, isoThermTemp
!!$
!!$    IMPLICIT NONE
!!$
!!$    TYPE(octal), INTENT(INOUT) :: thisOctal
!!$    INTEGER, INTENT(IN) :: subcell
!!$
!!$    REAL :: rho
!!$
!!$   ! TTauri star and disk
!!$   ! we may need to track the minimum and maximum densities of the
!!$   !   accretion flow 
!!$! RK COMMENTED OUT THIS. IBM XFL compiler does not allow intrinsic function 
!!$! in declearations.
!!$!   real, save      :: TTauriMinRho = huge(TTauriMinRho)
!!$!   real, save      :: TTauriMaxRho = tiny(TTauriMaxRho)
!!$    real, parameter :: TTauriMinRho = 1.0e25
!!$    real, parameter :: TTauriMaxRho = 1.0e-25
!!$
!!$    IF (isoTherm) THEN
!!$      thisOctal%temperature(subcell) = isoThermTemp
!!$      RETURN
!!$    END IF
!!$      
!!$    IF ( thisOctal%inFlow(subcell) ) THEN
!!$      rho = thisOctal%rho(subcell)
!!$      thisOctal%temperature(subcell) = MAX(5000.0, &
!!$        7000.0 - ((2500.0 * rho/real(mHydrogen) - TTauriMinRho) / (TTauriMaxRho-TTauriMinRho)))
!!$      ! we will initialise the bias distribution
!!$!      thisOctal%biasLine3D(subcell) = 1.0d0
!!$!      thisOctal%biasCont3D(subcell) = 1.0d0
!!$    ELSE
!!$      thisOctal%temperature(subcell) = 6500.0
!!$!      thisOctal%biasLine3D(subcell) = 1.0d-150
!!$!      thisOctal%biasCont3D(subcell) = 1.0d-150
!!$    END IF
!!$
!!$
!!$  
!!$  END SUBROUTINE calcTTauriTemperature

!!$  FUNCTION hartmannTemp(rM,inTheta,maxHartTemp)
!!$    ! the paper: Hartman, Hewett & Calvet 1994ApJ...426..669H 
!!$    !   quotes data from points along the magnetic field lines marked 
!!$    !   on one of the figures. this function returns the temperature
!!$    !   at any point within the accretion flow, by interpolating
!!$    !   between the values in the published figure 6.
!!$
!!$    USE unix_mod, only: unixgetenv
!!$
!!$    REAL             :: hartmannTemp
!!$    REAL, INTENT(IN) :: rM ! normalized (0=inside, 1=outside)
!!$    REAL, INTENT(IN) :: inTheta
!!$    REAL, INTENT(IN) :: maxHartTemp
!!$    
!!$    
!!$    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: temperatures
!!$    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: angles
!!$    INTEGER, DIMENSION(5), SAVE             :: nSamples
!!$    
!!$    LOGICAL, SAVE:: alreadyLoaded = .FALSE.
!!$    REAL :: radius
!!$    REAL :: theta
!!$    CHARACTER(LEN=80) :: dataDirectory 
!!$    INTEGER :: errNo, i
!!$    
!!$    INTEGER :: subPos
!!$    INTEGER :: fieldline
!!$    INTEGER :: lowerLine, upperLine
!!$    INTEGER :: iSample
!!$    INTEGER :: maxSamples
!!$    INTEGER :: lowerSample, upperSample
!!$    REAL :: lowerTemp, upperTemp
!!$    INTEGER :: iStat
!!$    logical, save :: warned_already_01 = .false.
!!$    logical, save :: warned_already_02 = .false.
!!$
!!$    IF (.NOT. alreadyLoaded) THEN
!!$      ! if this is the first time the function has been called, need
!!$      !   to load in the data.
!!$      
!!$      call unixGetenv("TORUS_DATA",dataDirectory,i,errno)
!!$      dataDirectory = trim(dataDirectory)//"/hartmann/"
!!$      OPEN(31,FILE=TRIM(dataDirectory)//"tempProfile1.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
!!$        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
!!$      OPEN(32,FILE=TRIM(dataDirectory)//"tempProfile2.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
!!$        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
!!$      OPEN(33,FILE=TRIM(dataDirectory)//"tempProfile3.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
!!$        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF  
!!$      OPEN(34,FILE=TRIM(dataDirectory)//"tempProfile4.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
!!$        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF  
!!$      OPEN(35,FILE=TRIM(dataDirectory)//"tempProfile5.csv",STATUS="OLD",FORM="FORMATTED",IOSTAT=iStat)
!!$        IF (iStat/=0) THEN; PRINT *, 'File open error in hartmannTemp: ',iStat; STOP ; END IF 
!!$
!!$      DO fieldline = 1, 5, 1
!!$        READ(UNIT=(30+fieldline),FMT=*) nSamples(fieldLine)
!!$      END DO
!!$
!!$      maxSamples = MAXVAL(nSamples)
!!$
!!$      ALLOCATE(temperatures(5,maxSamples))
!!$      ALLOCATE(angles(5,maxSamples))
!!$      ! make sure that there are no uninitialized variables
!!$      temperatures = 0.0
!!$      angles = 0.0
!!$      
!!$      DO fieldline = 1, 5, 1
!!$        DO iSample = 1,  nSamples(fieldLine), 1
!!$          READ(UNIT=(30+fieldline),FMT=*) angles(fieldLine,iSample),     &
!!$                                          temperatures(fieldLine,iSample)
!!$        END DO
!!$        CLOSE(UNIT=(30+fieldline))
!!$      END DO
!!$
!!$
!!$      ! might need to rescale the temperature distribution to have a new maximum
!!$      temperatures = temperatures * maxHartTemp/MAXVAL(temperatures)
!!$      
!!$      alreadyLoaded = .TRUE.
!!$    END IF
!!$
!!$    ! find which fieldlines bracket the point
!!$   
!!$    radius = rM * 0.8 ! because the Hartmann magnetic radii span 0.8 R_star
!!$    radius = radius / 0.2 ! divide by the field line spacing
!!$
!!$    IF (radius < 0.0) THEN 
!!$       if (.not. warned_already_01) then
!!$          PRINT *, 'In hartmannTemp, rM value is out of range! (',rM,')'
!!$          PRINT *, '  assuming fieldlines 1-2'
!!$          PRINT *, '  ==> Further warning surpressed.'
!!$          warned_already_01 = .true.
!!$       end if
!!$       lowerLine = 1
!!$       upperLine = 2 
!!$       subPos    = 0.0
!!$    ELSE IF (radius > 5.0) THEN
!!$       if (.not. warned_already_02) then
!!$          PRINT *, 'In hartmannTemp, rM value is out of range! (',rM,')'
!!$          PRINT *, '  assuming fieldlines 1-2'
!!$          PRINT *, '  ==> Further warning surpressed.'
!!$          warned_already_02 = .true.
!!$       end if
!!$       lowerLine = 4
!!$       upperLine = 5 
!!$    ELSE 
!!$       lowerLine = INT(radius) + 1
!!$       upperLine = lowerLine + 1 
!!$       subPos = radius - REAL(lowerline)
!!$    END IF
!!$      
!!$    IF (inTheta > pi) THEN
!!$      theta = inTheta - pi
!!$    ELSE
!!$      theta = inTheta
!!$    END IF
!!$    
!!$    IF (inTheta > pi/2.0) theta = pi - theta
!!$    
!!$    ! find the two temperatures at the fieldlines 
!!$    CALL locate(angles(lowerLine,1:nSamples(lowerLine)), &
!!$                nSamples(lowerLine), ABS(theta), lowerSample)
!!$    CALL locate(angles(upperLine,1:nSamples(upperLine)), &
!!$                nSamples(upperLine), ABS(theta), upperSample)
!!$
!!$
!!$    ! Forcing the values to be within the range...
!!$    ! ... quick fix to avoid out of range problem... (R. Kurosawa)
!!$    if (lowerSample <= 0) lowerSample =1
!!$    if (upperSample <= 0) upperSample =1
!!$    if (lowerSample >= nSamples(lowerLine)) lowerSample = nSamples(lowerLine)-1
!!$    if (upperSample >= nSamples(upperLine)) upperSample = nSamples(upperLine)-1
!!$
!!$                
!!$    IF ((lowerSample == 0) .OR. (lowerSample >= nSamples(lowerLine)) .OR. &
!!$       (upperSample == 0) .OR. (upperSample >= nSamples(upperLine))) THEN       
!!$       PRINT *, 'In hartmannTemp, theta value is out of range! (',ABS(theta),')'
!!$       STOP
!!$    END IF
!!$                
!!$    ! interpolate along each fieldline
!!$    lowerTemp = (temperatures(lowerLine,lowerSample) +        &
!!$                 temperatures(lowerLine,lowerSample+1)) / 2.0
!!$    upperTemp = (temperatures(upperLine,upperSample) +        &
!!$                 temperatures(upperLine,upperSample+1)) / 2.0
!!$   
!!$    ! interpolate between the fieldLines
!!$    
!!$    hartmannTemp = (lowerTemp * (1.0-subPos)) + (upperTemp * subPos)
!!$    
!!$  END FUNCTION hartmannTemp

!!$  SUBROUTINE writeHartmannValues(grid,variable)
!!$    ! writes out values from the grid at places defined by the
!!$    !   hartmannLines function
!!$    
!!$    use input_variables, only: TTauriRouter, TTauriRstar, &
!!$                               MdotType, curtainsPhi1s, curtainsPhi1e, &
!!$                               curtainsPhi2s, curtainsPhi2e
!!$    
!!$                               
!!$    TYPE(gridtype), INTENT(IN)   :: grid
!!$    CHARACTER(LEN=*), INTENT(IN) :: variable
!!$
!!$    INTEGER, PARAMETER :: nBins = 100
!!$    INTEGER, PARAMETER :: radiiPerBin = 1
!!$    INTEGER, PARAMETER :: samplesPerRadii = 10
!!$    INTEGER, PARAMETER :: samplesPerBin = radiiPerBin  * samplesPerRadii
!!$    INTEGER, PARAMETER :: nFieldLines = 5
!!$    
!!$    real(double), DIMENSION(nFieldLines,nBins) :: outputArray
!!$    real(double), DIMENSION(nFieldLines,nBins,grid%maxLevels) :: output2dArray
!!$    REAL, DIMENSION(nBins)             :: radii
!!$    TYPE(vector)                       :: azVec  ! azimuth vector
!!$    TYPE(vector)                       :: polVec ! poloidal vector
!!$    
!!$    INTEGER           :: iBin, iRadius, iSample, iFieldLine
!!$    REAL              :: rand
!!$    REAL(double)              :: phi
!!$    real(oct) :: radius
!!$    LOGICAL           :: ok
!!$    TYPE(vector) :: point 
!!$    REAL              :: value
!!$    real(double) :: valueDouble
!!$    TYPE(vector)      :: velocityValue
!!$    REAL(double)      :: etaLine
!!$    REAL              :: chiLine
!!$    REAL,DIMENSION(grid%maxLevels+1) :: departCoeff
!!$    real(double),DIMENSION(grid%maxLevels) :: valueArrayDouble
!!$    LOGICAL           :: scalar, logarithmic
!!$    CHARACTER(80)     :: tempChar
!!$    INTEGER           :: iLevel
!!$    LOGICAL, SAVE     :: firstRun = .TRUE.
!!$    REAL              :: curtainsFill
!!$    
!!$    scalar = .TRUE.
!!$    logarithmic = .FALSE.
!!$    curtainsFill = 1.
!!$    IF (variable == 'hartmann_departCoeff' .or. variable == 'hartmann_N') scalar = .false.
!!$   
!!$    IF (scalar) THEN
!!$      OPEN(UNIT=22, FILE=TRIM(variable)//'.csv', STATUS='REPLACE',FORM='FORMATTED')
!!$    ELSE
!!$      DO iLevel = 1, grid%maxLevels, 1
!!$        write(tempChar,'(i2.2)') iLevel
!!$        OPEN(UNIT=30+iLevel, FILE=TRIM(variable)//TRIM(tempChar)//'.csv', STATUS='REPLACE',FORM='FORMATTED')
!!$      END DO
!!$    END IF
!!$    
!!$    outputArray = 0.0
!!$    output2dArray = 0.0
!!$    
!!$    DO iBin = 1, nBins, 1
!!$    
!!$      DO iFieldLine = 1, nFieldLines, 1
!!$      
!!$        DO iRadius = 1, radiiPerBin
!!$
!!$          radius = 0.5*TTauriRstar + ((REAL(iBin-1)*REAL(radiiPerBin)+REAL(iRadius)) * &
!!$                     ((TTauriRouter - (0.1*TTauriRstar)) / REAL(nBins*radiiPerBin)))
!!$          radius = radius / TTauriRstar
!!$          radii(iBin) = radius
!!$        
!!$          DO iSample = 1, samplesPerRadii
!!$            IF (grid%amr2dOnly) THEN
!!$              phi = pi/4.0 
!!$            ELSE
!!$              IF (TRIM(mDotType) == "constantcurtains") THEN
!!$                IF (firstRun) THEN
!!$                  PRINT *, 'Using only region in curtains for ''writeHartmannValues'''
!!$                  curtainsFill = curtainsPhi1e-curtainsPhi1s 
!!$                  IF (curtainsPhi2e-curtainsPhi2s /= curtainsFill) &
!!$                    PRINT *, 'WARNING: curtains not equal size. Hartman values will be wrong'
!!$                  curtainsFill = curtainsFill / 180.
!!$                  PRINT *, '  curtains filling-factor: ',curtainsFill
!!$                    
!!$                  firstRun = .FALSE.
!!$                END IF
!!$                ok = .FALSE.
!!$                DO 
!!$                   call randomNumberGenerator(getReal=rand)
!!$                  phi = (rand-0.5) * twoPi
!!$                  IF (((phi > curtainsPhi1s+1.0).AND.(phi < curtainsPhi1e-1.0)) .OR.  &
!!$                      ((phi > curtainsPhi2s+1.0).AND.(phi < curtainsPhi2e-1.0))) ok = .TRUE.
!!$                  IF (ok) EXIT
!!$                END DO
!!$              ELSE
!!$                 call randomNumberGenerator(getReal=rand)
!!$                phi = (rand-0.5) * twoPi
!!$              END IF
!!$            END IF
!!$
!!$            CALL hartmannLines(iFieldline,radius,phi,grid,point,azVec,polVec,ok)
!!$            IF (.NOT. ok) EXIT
!!$
!!$            SELECT CASE (variable)
!!$
!!$            CASE ('hartmann_logNH')
!!$              CALL amrGridValues(grid%octreeRoot,point,rho=valueDouble) 
!!$              valueDouble = (valueDouble/(mHydrogen*curtainsFill))
!!$              logarithmic = .TRUE.
!!$            
!!$            CASE ('hartmann_temperature')
!!$              CALL amrGridValues(grid%octreeRoot,point,temperature=value) 
!!$              valueDouble = value / curtainsFill
!!$              
!!$            CASE ('hartmann_velPol')
!!$              CALL amrGridValues(grid%octreeRoot,point,velocity=velocityValue) 
!!$              value = velocityValue .dot. polVec
!!$              value = (value * cSpeed) * 1.e-5 ! km/s
!!$              valueDouble = abs(value) / curtainsFill
!!$              
!!$            CASE ('hartmann_velAz')
!!$              CALL amrGridValues(grid%octreeRoot,point,velocity=velocityValue) 
!!$              value = velocityValue .dot. azVec
!!$              value = (value * cSpeed) * 1.e-5 ! km/s
!!$              valueDouble = abs(value) / curtainsFill
!!$
!!$            CASE ('hartmann_line')
!!$              CALL amrGridValues(grid%octreeRoot,point,chiLine=chiLine,etaLine=etaLine)
!!$              if (chiLine*curtainsFill /= 0.) then
!!$                 valueDouble = etaLine / (chiLine*curtainsFill*10.) 
!!$              else
!!$                 valueDouble = tiny(valueDouble)
!!$              endif
!!$              logarithmic = .TRUE.
!!$              
!!$            CASE ('hartmann_Nelectron')
!!$              CALL amrGridValues(grid%octreeRoot,point,Ne=valueDouble)
!!$              valueDouble = valueDouble/curtainsFill
!!$              logarithmic = .TRUE.
!!$              
!!$            CASE ('hartmann_Nlevel2')
!!$              CALL amrGridValues(grid%octreeRoot,point,N=valueArrayDouble)
!!$              valueDouble = valueArrayDouble(2)/curtainsFill
!!$              logarithmic = .TRUE.
!!$              
!!$            CASE ('hartmann_logNe')
!!$              CALL amrGridValues(grid%octreeRoot,point,Ne=valueDouble)
!!$              valueDouble = valueDouble/curtainsFill
!!$              logarithmic = .TRUE.
!!$              
!!$            CASE ('hartmann_NH')
!!$              CALL amrGridValues(grid%octreeRoot,point,rho=valueDouble)
!!$              valueDouble = valueDouble/(mHydrogen*curtainsFill)
!!$              logarithmic = .TRUE.
!!$              
!!$            CASE ('hartmann_departCoeff')
!!$              scalar = .false.
!!$              CALL amrGridValues(grid%octreeRoot,point,departCoeff=departCoeff)
!!$              valueArrayDouble = departCoeff(1:grid%maxLevels) / curtainsFill
!!$              
!!$            CASE ('hartmann_N')
!!$              scalar = .false.
!!$              CALL amrGridValues(grid%octreeRoot,point,N=valueArrayDouble)
!!$              valueArrayDouble = valueArrayDouble / curtainsFill
!!$              !logarithmic = .TRUE.
!!$              
!!$            CASE DEFAULT
!!$              PRINT *, '''variable'' not recognized in writeHartmannValues'
!!$              STOP
!!$              
!!$            END SELECT
!!$
!!$            IF (scalar) THEN
!!$              outputArray(iFieldLine,iBin) = outputArray(iFieldLine,iBin) + valueDouble 
!!$            ELSE
!!$              output2dArray(iFieldLine,iBin,:) = output2dArray(iFieldLine,iBin,:) + valueArrayDouble 
!!$            END IF
!!$
!!$          END DO ! iSample
!!$          
!!$          IF (.NOT. ok) EXIT
!!$        
!!$        END DO ! iRadius
!!$        
!!$        IF (ok) THEN
!!$          IF (scalar) THEN
!!$            outputArray(iFieldLine,iBin) = outputArray(iFieldLine,iBin) / REAL(samplesPerBin,KIND=db)
!!$            IF (logarithmic) outputArray = LOG10(MAX(1.0e-10_db,outputArray))
!!$          ELSE
!!$            output2dArray(iFieldLine,iBin,:) = output2dArray(iFieldLine,iBin,:) / REAL(samplesPerBin,KIND=db)
!!$            IF (logarithmic) output2dArray = LOG10(MAX(1.0e-10_db,output2dArray))
!!$          END IF
!!$        ELSE 
!!$          IF (scalar) THEN
!!$            outputArray(iFieldLine,iBin) = 0.0
!!$            !outputArray(iFieldLine,iBin) = -1.0 * HUGE(outputArray(1,iBin))
!!$          ELSE
!!$            output2dArray(iFieldLine,iBin,:) = 0.0
!!$            !output2dArray(iFieldLine,iBin,:) = -1.0 * HUGE(output2dArray(1,iBin,1))
!!$          END IF
!!$        END IF
!!$
!!$      END DO ! iFieldLine
!!$      
!!$      IF (scalar) THEN
!!$        WRITE(22,'(F11.4,3(" ",E11.4))') radii(iBin), outputArray(2:4,iBin)
!!$        !                ^                                  ^^^   
!!$        ! these constants might change!
!!$      ELSE
!!$        DO iLevel = 1, grid%maxLevels, 1
!!$          WRITE(30+iLevel,'(F11.4,3(" ",E11.4))') radii(iBin), REAL(output2dArray(2:4,iBin,iLevel))
!!$        END DO
!!$      END IF
!!$    
!!$    END DO ! iBin
!!$    
!!$    IF (scalar) THEN
!!$      CLOSE(UNIT=22)
!!$    ELSE
!!$      DO iLevel = 1, grid%maxLevels, 1
!!$        CLOSE(UNIT=30+iLevel)
!!$      END DO
!!$    END IF
!!$
!!$  END SUBROUTINE writeHartmannValues    

!!$  subroutine checkMassLossRate(grid)
!!$    use input_variables, only : CMFGEN_Rmin
!!$    type(GRIDTYPE) :: grid
!!$    real(double) :: mu, dmu, r, sintheta, tot
!!$    type(OCTAL), pointer :: thisOctal
!!$    integer :: subcell
!!$    integer :: i, j
!!$    type(VECTOR) :: rVec
!!$
!!$    do j = 1, 100
!!$       r = 1.01d0 * CMFGEN_Rmin * dble(j)
!!$       tot = 0.d0
!!$       do i = 1, 1000
!!$          mu = -1.d0 + 2.d0*dble(i-1)/999.d0
!!$          dmu = 2.d0/999.d0
!!$          sinTheta = sqrt(1.d0-mu**2)
!!$          rVec = VECTOR(r*sinTheta, 0.d0, r*mu)
!!$          CALL findSubcellTD(rVec, grid%octreeRoot, thisOctal, subcell)
!!$          tot = tot + thisOctal%rho(subcell)*modulus(thisOctal%velocity(subcell))*cSpeed * r**2 * dmu
!!$       enddo
!!$       tot = tot * 1.d20
!!$       tot = tot * twoPi
!!$       tot = (tot / mSol) / secstoyears
!!$       write(*,*) j
!!$       write(*,'(a,1pe12.3,a)') "Mass-loss rate (assuming fully ionized pure helium): ", tot, " solar masses/year"
!!$    enddo
!!$  end subroutine checkMassLossRate

    
!!$  TYPE(vector)  function gammaVelVelocity(point, grid)
!!$    use input_variables, only : vnought1, vnought2, vterm1, vterm2, beta1, beta2
!!$    use input_variables, only : rstar1, rstar2, mass1, mass2, mdot1, mdot2, binarySep
!!$
!!$    type(vector), intent(in) :: point
!!$    type(vector) :: rvec
!!$    type(GRIDTYPE), intent(in) :: grid
!!$    real(double) :: v, r
!!$    real :: massRatio
!!$    integer, parameter :: nsteps = 1000
!!$    real :: ydist(nSteps), xDist(nSteps)
!!$    real :: d1, d2, curlyR, dx, dybydx, ddash1, ddash2
!!$    real(double) :: momRatio
!!$    type(VECTOR) :: direction(nSteps), stagVec, shockdirection
!!$    type(VECTOR) :: direction2
!!$    real :: stagpoint
!!$    integer :: i, j
!!$
!!$
!!$    rVec = point
!!$
!!$    gammaVelVelocity = VECTOR(0.d0, 0.d0, 0.d0)
!!$    r = modulus(rVec)
!!$
!!$    massRatio = mass1/mass2
!!$
!!$    d1 = binarySep * (1./(massRatio+1.))
!!$    d2 = binarySep - d1
!!$
!!$    momRatio = (mdot1 * vterm1) / (mdot2 * vterm2)
!!$
!!$    curlyR = sqrt(momRatio)         ! = d1/d2 (Equ 1.) 
!!$    ! Stevens, Blondin & Pollock 1992
!!$    ! ApJ 386 265
!!$
!!$    stagPoint = binarySep * sqrt(momRatio) / (1. + sqrt(momRatio))
!!$
!!$    stagVec = VECTOR(0.,0.,-d1) + VECTOR(0., 0., stagPoint)
!!$
!!$    shockDirection = VECTOR(0.,0.,0.)
!!$
!!$
!!$    dx = (grid%octreeRoot%subcellSize-stagVec%z)/real(nsteps)
!!$    yDist(1) = dx*10.
!!$    xDist(1) = stagVec%z 
!!$    direction(1) = VECTOR(0.,1.,0.)
!!$    do j = 2,nsteps
!!$       xDist(j) = stagVec%z + real(j-1)*(grid%octreeroot%subcellsize-stagVec%z)/real(nSteps-1)
!!$       ddash1 = sqrt((xDist(j-1) + d1)**2 + yDist(j-1)**2)
!!$       ddash2 = sqrt((xDist(j-1) - d2)**2 + yDist(j-1)**2)
!!$
!!$       dybydx = ((curlyR*ddash2**2 + ddash1**2)*yDist(j-1)) / &
!!$            (curlyR *ddash2**2*(xDist(j-1)+d1) + ddash1**2*(xDist(j-1)-d2))
!!$
!!$       yDist(j) = yDist(j-1) + dx * dyBydx
!!$
!!$       direction(j) = VECTOR(dx,dybydx,0.)
!!$       call normalize(direction(j))
!!$    enddo
!!$
!!$   
!!$    if (rVec%z > xDist(1)) then
!!$       call locate(xDist, nSteps, real(rVec%z), i)
!!$       if (rvec%x > ydist(i)) then
!!$          direction2 = (rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
!!$          call normalize(direction2)
!!$          r = modulus(rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
!!$          if (r > rStar1) then
!!$             v = vNought1 + (vterm1-vNought1)*(1.d0 - rstar1/r)**beta1
!!$             gammaVelVelocity = dble(v/cspeed) * direction2
!!$          endif
!!$       else
!!$          direction2 = (rVec - VECTOR(0.d0, 0.d0, dble(d2)))
!!$          call normalize(direction2)
!!$          r = modulus(rVec - VECTOR(0.d0, 0.d0, dble(d2)))
!!$          if (r > rStar2) then
!!$             v = vNought2 + (vterm2-vNought2)*(1.d0 - rstar2/r)**beta2
!!$             gammaVelVelocity = dble(v/cspeed) * direction2
!!$          endif
!!$       endif
!!$    else
!!$       direction2 = (rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
!!$       call normalize(direction2)
!!$       r = modulus(rVec - VECTOR(0.d0, 0.d0, -dble(d1)))
!!$       if (r > rStar1) then
!!$          v = vNought1 + (vterm1-vNought1)*(1.d0 - rstar1/r)**beta1
!!$          gammaVelVelocity = dble(v/cspeed) * direction2
!!$       endif
!!$    endif
!!$
!!$
!!$  end function gammaVelVelocity
!!$
!!$

!
! Computes average temperature (T_ave) and mass-weighted temperature (T_mass)
!
!!$  subroutine find_average_temperature(grid, T_ave, T_mass, TotalMass)
!!$    implicit none
!!$    type(gridtype), intent(in)  :: grid
!!$    real(double), intent(out) :: T_ave   ! average temperature in [k}
!!$    real(double), intent(out) :: T_mass  ! mass weighted temperature in [K]
!!$    real(double), intent(out) :: TotalMass  ! total mass [g]
!!$    real(double) :: sum_T    ! Sum of temperatures in all cells [K]
!!$    real(double) :: sum_M    ! Sum of mass in all cells [g]
!!$    real(double) :: sum_TM   ! Sum of temperatures*Mass in all cells [K*g]
!!$    integer :: ncell
!!$
!!$    ! initialize some values
!!$    sum_T=0.0d0; sum_TM = 0.0d0; ncell=0; sum_M=0.0d0
!!$    
!!$    ! calling the routine in this module to do sum of T and T*M
!!$    call sum_temp_mass(grid%octreeroot, sum_T, sum_TM, sum_M, ncell)
!!$
!!$    ! outputs
!!$    T_ave = sum_T/ncell  ![K]
!!$    T_mass = sum_TM/sum_M
!!$    TotalMass = sum_M
!!$
!!$  end subroutine find_average_temperature
!!$



  !
  ! Computes the average, maximum, minimum, contiuum optical depth across each cell, 
!!$  subroutine find_max_min_ave_tau(grid, ilambda, tau_max, tau_min, tau_ave)
!!$    implicit none
!!$    type(gridtype), intent(in)  :: grid
!!$    integer, intent(in)         :: ilambda  ! wavelength index
!!$    real(double), intent(out) :: tau_max
!!$    real(double), intent(out) :: tau_min
!!$    real(double), intent(out) :: tau_ave
!!$    ! 
!!$    integer :: ncell   ! number of cells
!!$
!!$    ! initialize some values
!!$    tau_max = -1.0d10; tau_min=1.0d10; tau_ave=0.0d0; ncell=0    
!!$
!!$    ! calling the routine in this module to do sum of T and T*M
!!$    call max_min_ave_tau(grid%octreeroot, grid, ilambda, tau_max, tau_min, tau_ave, ncell)
!!$
!!$
!!$  end subroutine find_max_min_ave_tau

  ! find the average, maximum, minimum, contiuum optical depth across each cell,
!!$  recursive subroutine max_min_ave_tau(thisOctal, grid, ilambda, tau_max, tau_min, tau_ave, ncell)
!!$    implicit none
!!$    type(octal), pointer   :: thisOctal
!!$    type(gridtype), intent(in)  :: grid    
!!$    integer, intent(in)  :: ilambda  ! wavelength index at which tau is computed.
!!$    real(double), intent(inout) :: tau_max
!!$    real(double), intent(inout) :: tau_min
!!$    real(double), intent(inout) :: tau_ave
!!$    integer, intent(inout) :: ncell ! number of cells used for this calculation
!!$    !
!!$    type(octal), pointer  :: child
!!$    real(double) :: chi, tau
!!$    integer :: subcell, i
!!$  
!!$    do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call max_min_ave_tau(child, grid, ilambda,  &
!!$                     tau_max, tau_min, tau_ave, ncell)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          if (thisOctal%inFlow(subcell)) then
!!$             if (.not.grid%oneKappa) then
!!$                chi = thisOctal%kappaSca(subcell,iLambda) &
!!$                     + thisOctal%kappaAbs(subcell,iLambda)
!!$             else
!!$                chi = ( grid%oneKappaSca(thisOctal%dustType(subcell),iLambda) + &
!!$                     grid%oneKappaAbs(thisOctal%dustType(subcell),iLambda) ) &
!!$                     *thisOctal%rho(subcell)
!!$             endif
!!$             ! chi is normal chi times 10^10cm, and sucellsize is in 10^10 cm.
!!$             ! so tau should be in correct scale here...
!!$             tau = chi*thisOctal%subcellSize
!!$             
!!$             ! now update the max, min and average tau values.
!!$             tau_max = MAX(tau_max, tau)
!!$             tau_min = MIN(tau_min, tau)
!!$             
!!$             ncell = ncell + 1
!!$             
!!$             tau_ave = (tau_ave*dble(ncell-1) + tau)/dble(ncell)
!!$          else
!!$             continue
!!$          endif
!!$       endif
!!$    enddo
!!$  end subroutine max_min_ave_tau

  ! chris (19/05/04)
  ! Smooths the AMR grid to ensure 'adequate' resolution at the boundary
  ! between optically thin and optically thick regions. This routine should be
  ! called repeatedly until gridConverged returns true (it should be set to
  ! true in the initial call.
  ! The value of kappa is the opacity of the grid at a test wavelength (set in
  ! the parameter file).

!!$  RECURSIVE SUBROUTINE smoothAMRgridTau(thisOctal,grid,gridConverged, ilam, &
!!$                                        stellar_cluster, &
!!$                                        inheritProps, interpProps, romData)
!!$    USE input_variables, ONLY : tauSmoothMax,tauSmoothMin
!!$    IMPLICIT NONE
!!$
!!$    TYPE(octal), POINTER             :: thisOctal
!!$    TYPE(gridtype), INTENT(INOUT   ) :: grid 
!!$    LOGICAL, INTENT(INOUT)               :: gridConverged
!!$    TYPE(cluster), optional, intent(in)  :: stellar_cluster
!!$    LOGICAL, INTENT(IN) :: inheritProps
!!$    LOGICAL, INTENT(IN), optional :: interpProps
!!$    TYPE(romanova), optional, INTENT(IN)   :: romData  ! used for "romanova" geometry
!!$
!!$    INTEGER              :: i, ilam
!!$    TYPE(vector)    :: thisSubcellCentre
!!$    REAL(oct) :: dSubcellCentre
!!$    real(double) :: kappaAbs, kappaSca
!!$    TYPE(octal), POINTER :: child
!!$    TYPE(octal), POINTER :: neighbour
!!$    TYPE(vector), ALLOCATABLE, DIMENSION(:) :: locator
!!$    INTEGER              :: subcell
!!$    INTEGER              :: thisSubcell
!!$    REAL                 :: thisTau, thatTau !, tauDiff
!!$    INTEGER              :: nLocator
!!$    LOGICAL              :: prob
!!$
!!$    kappaAbs = 0.d0; kappaSca = 0.d0
!!$    ! For each subcell, we find the coordinates of at least one point in every
!!$    ! neighbouring subcell. The optical depths are compared and if one cell is
!!$    ! optically thick while the other is optically thin with an optical depth
!!$    ! less than some value, both cells are split. There may be up to 26
!!$    ! neighbouring subcells.
!!$
!!$    if (thisOctal%threed) then
!!$       nLocator = 26
!!$    else
!!$       nlocator = 8
!!$    endif
!!$
!!$    ALLOCATE(locator(nLocator))
!!$
!!$    thisSubcell = 1
!!$    do
!!$      if (thisOctal%hasChild(thisSubcell)) then
!!$        ! find the child
!!$        do i = 1, thisOctal%nChildren, 1
!!$           if (thisOctal%indexChild(i) == thisSubcell) then
!!$              child => thisOctal%child(i)
!!$              call smoothAMRgridTau(child,grid,gridConverged,  ilam, stellar_cluster, &
!!$                   inheritProps, interpProps, romData)
!!$              ! force return until algorithm is fixed
!!$              IF ( .NOT. gridConverged ) RETURN
!!$              exit
!!$           end if
!!$        end do
!!$      else
!!$        thisSubcellCentre = subcellCentre(thisOctal,thisSubcell)
!!$
!!$        ! don't split outer edge of disc
!!$
!!$        if ((grid%geometry == "shakara").or.(grid%geometry == "warpeddisc").or.(grid%geometry == "iras04158").or.&
!!$             (grid%geometry=="circumbin")) then
!!$           if (sqrt(thissubcellcentre%x**2 + thissubcellcentre%y**2) > grid%router*0.9) goto 666
!!$        endif
!!$           
!!$
!!$        FORALL (i = 1:nLocator)
!!$          locator(i) = thisSubcellCentre
!!$        END FORALL
!!$
!!$        ! Moving this distance in any direction will take us into a
!!$        ! neighbouring subcell.
!!$        dSubcellCentre = (0.5 * thisOctal%subcellSize) + grid%halfSmallestSubcell
!!$
!!$        if (thisOctal%threed) then 
!!$
!!$           ! faces
!!$           locator(1)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(2)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(3)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(4)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(5)%y = thisSubcellCentre%y - dSubcellCentre
!!$           locator(6)%z = thisSubcellCentre%z - dSubcellCentre
!!$           ! x-edges
!!$           locator(7)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(7)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(8)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(8)%z = thisSubcellCentre%z - dSubcellCentre
!!$           locator(9)%y = thisSubcellCentre%y - dSubcellCentre
!!$           locator(9)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(10)%y = thisSubcellCentre%y - dSubcellCentre
!!$           locator(10)%z = thisSubcellCentre%z - dSubcellCentre
!!$           ! y-edges
!!$           locator(11)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(11)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(12)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(12)%z = thisSubcellCentre%z - dSubcellCentre
!!$           locator(13)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(13)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(14)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(14)%z = thisSubcellCentre%z - dSubcellCentre
!!$           ! z-edges
!!$           locator(15)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(15)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(16)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(16)%y = thisSubcellCentre%y - dSubcellCentre
!!$           locator(17)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(17)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(18)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(18)%y = thisSubcellCentre%y - dSubcellCentre
!!$           ! corners
!!$           locator(19)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(19)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(19)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(20)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(20)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(20)%z = thisSubcellCentre%z - dSubcellCentre
!!$           locator(21)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(21)%y = thisSubcellCentre%y - dSubcellCentre
!!$           locator(21)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(22)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(22)%y = thisSubcellCentre%y - dSubcellCentre
!!$           locator(22)%z = thisSubcellCentre%z - dSubcellCentre
!!$           locator(23)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(23)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(23)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(24)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(24)%y = thisSubcellCentre%y + dSubcellCentre
!!$           locator(24)%z = thisSubcellCentre%z - dSubcellCentre
!!$           locator(25)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(25)%y = thisSubcellCentre%y - dSubcellCentre
!!$           locator(25)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(26)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(26)%y = thisSubcellCentre%y - dSubcellCentre
!!$           locator(26)%z = thisSubcellCentre%z - dSubcellCentre
!!$        else
!!$
!!$           ! edges
!!$
!!$           locator(1)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(2)%z = thisSubcellCentre%z + dSubcellCentre
!!$           locator(3)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(4)%z = thisSubcellCentre%z - dSubcellCentre
!!$           ! corners
!!$           locator(5)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(5)%z = thisSubcellCentre%z + dSubcellCentre
!!$
!!$           locator(6)%x = thisSubcellCentre%x + dSubcellCentre
!!$           locator(6)%z = thisSubcellCentre%z - dSubcellCentre
!!$
!!$           locator(7)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(7)%z = thisSubcellCentre%z + dSubcellCentre
!!$
!!$           locator(8)%x = thisSubcellCentre%x - dSubcellCentre
!!$           locator(8)%z = thisSubcellCentre%z - dSubcellCentre
!!$
!!$
!!$        endif
!!$
!!$        call returnKappa(grid, thisOctal, thissubcell, ilambda=ilam,&
!!$             kappaSca=kappaSca, kappaAbs=kappaAbs)
!!$
!!$        thisTau = (kappaAbs + kappaSca) * thisOctal%subcellSize !* thisOctal%rho(thisSubcell)
!!$
!!$        DO i = 1, nLocator, 1
!!$          IF ( inOctal(grid%octreeRoot,locator(i)) ) THEN
!!$            neighbour => thisOctal
!!$            CALL findSubcellLocal(locator(i),neighbour,subcell, prob=prob)
!!$                IF ( neighbour%hasChild(subcell) ) THEN 
!!$                  PRINT *, "neighbour already has child (a)",prob
!!$                  do;enddo
!!$                END IF
!!$            
!!$            ! The neighbouring subcell must be larger than the current subcell
!!$            ! otherwise our locators won't cover all neighbouring subcells
!!$            ! (and we'll hit cell boundaries, which is not good).
!!$            IF ( neighbour%subcellSize >= thisOctal%subcellSize) THEN
!!$
!!$               call returnKappa(grid, neighbour, subcell, &
!!$                    ilambda=ilam,  kappaSca=kappaSca, kappaAbs=kappaAbs)
!!$               thatTau = (kappaSca + kappaAbs) * neighbour%subcellSize! * neighbour%rho(subcell)
!!$
!!$              ! Original critera was to split if:
!!$              ! ((o1 - o2)/(o1+o2) > 0.5 .and. (o1 - o2) > 1)
!!$!              tauDiff = abs(thisTau - thatTau)
!!$!!              if ((tauDiff.gt.1000000.0).and. &
!!$!              if ((tauDiff.gt.1.).and. &
!!$!!                  .not.((thisTau.gt.5.).and.(thatTau.gt.5.)).and. &
!!$!                  ((tauDiff/(thisTau + thatTau)).gt.0.5)) then
!!$
!!$              ! Critera for cell splitting is that one cell is optically think
!!$              ! (tau > 1) and the other is optically thin with an optical
!!$              ! depth less than some value (which should really be set in the
!!$              ! parameter file).
!!$!              if (max(thisTau, thatTau).gt.0.5.and.min(thisTau, thatTau).lt.0.01) then
!!$!              if ((max(thisTau, thatTau).gt.tauSmoothMax.and.min(thisTau, thatTau).lt.tauSmoothMin).or.&
!!$!                   ((max(thisTau, thatTau)<5.).and.(abs(thistau-thatTau)> 2.))) then
!!$
!!$               if ((max(thisTau, thatTau).gt.tauSmoothMax.and.min(thisTau, thatTau).lt.tauSmoothMin)) then
!!$
!!$
!!$!write (*,*) thisSubcellCentre%x, thisSubcellCentre%y, thisSubcellCentre%z, thisTau, thatTau, i
!!$                ! Because addNewChild is broken, we must use addNewChildren below,
!!$                ! which makes these tests on the current and neighbouring subcells
!!$                ! insufficient. We perform them anyway.
!!$                IF ( neighbour%hasChild(subcell) ) THEN 
!!$                  PRINT *, "neighbour already has child"
!!$                  STOP
!!$                END IF
!!$                IF ( thisOctal%hasChild(thisSubcell) ) THEN 
!!$                  PRINT *, "thisOctal already has child"
!!$                  STOP
!!$                END IF
!!$                ! Note that the version of addNewChild called here does not
!!$                ! support sphData, etc.
!!$!                call addNewChild(neighbour,subcell,grid)
!!$!                call addNewChild(thisOctal,thisSubcell,grid)
!!$
!!$              call addNewChild(neighbour, subcell, grid, adjustGridInfo=.TRUE.,  &
!!$                               inherit=inheritProps, interp=interpProps, romData=romData)
!!$!              call addNewChild(thisOctal, thissubcell, grid, adjustGridInfo=.TRUE.,  &
!!$!                               sphData=sphData, stellar_cluster=stellar_cluster, &
!!$!                               inherit=inheritProps, interp=interpProps, romData=romData)
!!$!              write(*,*) "added new child",thisOctal%hasChild,thisOctal%nchildren
!!$
!!$
!!$!                call addNewChildren(thisOctal, grid, sphData, stellar_cluster, inheritProps, interpProps, romData)
!!$                ! If we are splitting two subcells in the same octal, we don't
!!$                ! need to addNewChildren to the neighbour. If we switch to
!!$                ! addNewChild this test becomes redundant.
!!$!                if (.not. associated(thisOctal, neighbour)) then
!!$!                  call addNewChildren(neighbour, grid, sphData, stellar_cluster, inheritProps, interpProps, romData)
!!$!                end if 
!!$
!!$                gridConverged = .FALSE.
!!$                if (.not.gridConverged) return!!!!
!!$                exit ! loop through locators, move onto next subcell
!!$              end if ! grid must be refined
!!$            ENDIF ! neighbour subcell is larger
!!$          END IF ! in grid
!!$        END DO ! loop through locators
!!$
!!$      end if ! thisOctal%hasChild(thisSubcell)
!!$      thisSubcell = thisSubcell + 1 ! next subcell
!!$      if (thisSubcell > thisOctal%maxChildren) exit ! loop through subcells in an octal
!!$    end do ! loop through subcells in an octal
!!$666 continue
!!$    deallocate(locator)
!!$  END SUBROUTINE smoothAMRgridTau
!!$
!!$  subroutine polarDump(grid)
!!$    implicit none
!!$    type(GRIDTYPE) :: grid
!!$    integer :: nt = 200, nr = 400, i, j
!!$    real :: theta, r, t
!!$    real(double) :: dens
!!$    type(octal), pointer   :: thisOctal
!!$    integer :: subcell
!!$    type(VECTOR) :: rVec
!!$
!!$    open(30, file="polardump.dat", status="unknown", form="formatted")
!!$    do j = 1, nt
!!$       do i = 1, nr
!!$          r = grid%rInner + (grid%rOuter-grid%rInner)*(real(i-1)/real(nr-1))**3
!!$          theta = pi*real(j-1)/real(nt-1)
!!$          rVec = VECTOR(dble(r*sin(theta)),0.d0,dble(r*cos(theta)))
!!$          call amrGridValues(grid%octreeRoot, rVec,temperature=t, rho=dens, foundoctal=thisOctal, foundsubcell=subcell)
!!$          if (thisOCtal%inFlow(subcell)) then
!!$             write(30,*) r/1.5e3, t, dens
!!$          endif
!!$       enddo
!!$    enddo
!!$    close(30)
!!$
!!$    open(30, file="vert1.dat", status="unknown", form="formatted")
!!$    call verticalDump(grid%octreeRoot, grid%rInner*2.)
!!$    close(30)
!!$
!!$    open(30, file="vert2.dat", status="unknown", form="formatted")
!!$    call verticalDump(grid%octreeRoot, grid%rOuter/2.)
!!$    close(30)
!!$    open(30,file="radial.dat",status="unknown",form="formatted")
!!$    do i = 1, nr
!!$       r = grid%rInner + (grid%rOuter-grid%rInner)*(real(i-1)/real(nr-1))**3
!!$       rVec = VECTOR(dble(r), 0.d0, 0.d0)
!!$       call amrGridValues(grid%octreeRoot, rVec,temperature=t, rho=dens, foundoctal=thisOctal, foundsubcell=subcell)
!!$       write(30,*) r/1496., t, dens
!!$    enddo
!!$    close(30)
!!$
!!$
!!$  end subroutine polarDump
!!$

!!$  recursive subroutine unrefineThickCells(thisOctal, grid, ilambda, converged)
!!$    type(GRIDTYPE) :: grid
!!$    type(octal), pointer   :: thisOctal
!!$    type(octal), pointer  :: child
!!$    integer :: ilambda
!!$    real(double) :: kappaAbs, kappaSca, tau
!!$    integer :: subcell, i
!!$    logical :: unrefine, converged
!!$    kappaAbs =0.d0; kappaSca = 0.d0
!!$    unrefine = .true.
!!$
!!$    do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call unrefineThickCells(child, grid, ilambda, converged)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          if (.not.ASSOCIATED(thisOctal%dustTypeFraction)) then
!!$             write(*,*) "unalloc dusttypefraction!!"
!!$          endif
!!$          call returnKappa(grid, thisOctal, subcell, ilambda, kappaAbs=kappaAbs,kappaSca=kappaSca)
!!$          tau = thisOctal%subcellSize*(kappaAbs+kappaSca)
!!$          if (tau < 1.e4) unrefine = .false.
!!$       endif
!!$    enddo
!!$
!!$    if ((thisOctal%nChildren == 0).and.unrefine.and.converged) then
!!$       call deleteChild(thisOctal%parent, thisOctal%parentSubcell, adjustParent = .true., &
!!$            grid = grid, adjustGridInfo = .true.)
!!$       converged = .false.
!!$    endif
!!$
!!$  end subroutine unrefineThickCells


!!$  recursive subroutine biasCentreVolume(thisOctal, boxSize)
!!$    real :: boxSize
!!$    type(VECTOR) :: rVec
!!$    type(octal), pointer   :: thisOctal
!!$    type(octal), pointer  :: child 
!!$    integer :: subcell, i
!!$    
!!$    do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call biasCentreVolume(child, boxSize)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          rVec = subcellCentre(thisOctal, subcell)
!!$          if ( ((rVec%x < boxSize/2.).and.(rVec%x > -boxsize/2.)) .and. &
!!$               ((rVec%y < boxSize/2.).and.(rVec%y > -boxsize/2.)) .and. &
!!$               ((rVec%z < boxSize/2.).and.(rVec%z > -boxsize/2.)) ) then
!!$             thisOctal%biasCont3D(subcell) = 100.
!!$          endif
!!$       endif
!!$          
!!$    enddo
!!$
!!$  end subroutine biasCentreVolume
!!$
!!$
!!$  recursive subroutine setbiasAMR(thisOctal, grid)
!!$  type(gridtype) :: grid
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  integer :: subcell, i
!!$  real :: r
!!$  
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call setBiasAMR(child, grid)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          r = modulus(subcellcentre(thisOctal, subcell)) / grid%rInner
!!$          thisOctal%biasCont3D(subcell) = r
!!$       endif
!!$    enddo
!!$  end subroutine setBiasAMR
!!$
!!$  recursive subroutine set_bias_shakara(thisOctal, grid, ilam, ross)
!!$  type(gridtype) :: grid
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  real(double) :: tau, kappaAbs, kappaSca !, thisTau
!!$!  type(vector) :: rVec, direction
!!$  logical :: ross
!!$  integer :: subcell, i, ilam !, ntau
!!$!  real :: r
!!$  kappaAbs = 0.d0; kappaSca = 0.d0
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call set_bias_shakara(child, grid, ilam, ross)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          if (.not.ross) then
!!$             call returnKappa(grid, thisOctal, subcell, ilam, kappaAbs=kappaAbs, kappaSca=kappaSca)
!!$             tau = thisOctal%subcellSize*(kappaAbs + kappaSca)
!!$             thisOctal%biasCont3D(subcell) = exp(-tau)
!!$          else
!!$!             rVec = subcellCentre(thisOctal, subcell)
!!$!             tau = 1.d30
!!$!             ntau = 20
!!$!             do i = 1, nTau
!!$!                direction = randomUnitVector()
!!$!                call tauAlongPath(grid, rVec, direction, thistau, 100.d0)
!!$!                tau = min(tau, thisTau)
!!$!             enddo
!!$!             write(*,*) tau
!!$             thisOctal%biasCont3D(subcell) = 1.d0/(cellVolume(thisOCtal,subcell)*thisOctal%etaCont(subcell))
!!$          endif
!!$!          r = rVec%x
!!$!          thisOCtal%biasCont3d(subcell) = thisOctal%biasCont3d(subcell) * r
!!$       endif
!!$
!!$    enddo
!!$  end subroutine set_bias_shakara
!!$
!!$  recursive subroutine set_bias_rosseland(thisOctal, grid)
!!$  type(gridtype) :: grid
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  real(double) :: tau, kappaAbs
!!$  integer :: subcell, i, ilam
!!$  ilam = 0; kappaAbs = 0.d0
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call set_bias_rosseland(child, grid)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          call returnKappa(grid, thisOctal, subcell, ilam, rosselandKappa = kappaAbs)
!!$          tau = thisOctal%subcellSize*kappaAbs*thisOctal%rho(subcell)*1.d10
!!$          thisOctal%biasCont3D(subcell) = MAX(exp(-tau),1.d-8)
!!$       endif
!!$    enddo
!!$  end subroutine set_bias_rosseland
!!$
!!$  recursive subroutine set_bias_radius(thisOctal, p)
!!$  type(octal), pointer  :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  integer :: subcell, i
!!$  integer :: p
!!$  
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call set_bias_radius(child, p)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          thisOctal%biasCont3D(subcell) = thisOctal%biasCont3D(subcell) &
!!$               * modulus(subcellCentre(thisOctal, subcell))**p
!!$       endif
!!$    enddo
!!$  end subroutine set_bias_radius
!!$
!!$  recursive subroutine set_bias_whitney(thisOctal, grid)
!!$  use input_variables, only : erOuter
!!$  type(gridtype) :: grid
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  real(double) :: tau, kappaAbs
!!$  integer :: subcell, i, ilam
!!$  real(double) :: r
!!$  kappaAbs = 0.d0; ilam = 0
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call set_bias_whitney(child, grid)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          r = modulus(subcellCentre(thisoctal,subcell))
!!$          call returnKappa(grid, thisOctal, subcell, ilam, rosselandKappa = kappaAbs)
!!$             tau = thisOctal%subcellSize*kappaAbs*thisOctal%rho(subcell)*1.d10
!!$             thisOctal%biasCont3D(subcell) = MAX(exp(-tau),1.d-4) * (1.d10*r/erOuter)
!!$!             thisOctal%biasCont3D(subcell) = 1.d0/(thisOctal%etacont(subcell)*cellVolume(thisOctal, subcell))
!!$       endif
!!$    enddo
!!$  end subroutine set_bias_whitney
!!$
!!$  recursive subroutine set_bias_ttauri(thisOctal, grid, lambda0, outVec)
!!$  type(gridtype) :: grid
!!$  type(octal), pointer   :: thisOctal
!!$  real, intent(in)       :: lambda0                ! rest wavelength of line
!!$  type(octal), pointer  :: child 
!!$  integer :: subcell, i
!!$  real(double) :: d, dV, r, nu0, tauSob, escProb
!!$!  real(double) :: tauSob_tmp, tauSob_sum, theta, phi, sin_theta
!!$  type(vector)  :: rvec, rhat !, rvec_sample
!!$!  integer :: nsample
!!$  type(VECTOR) :: outVec
!!$  real(double):: dtau_cont, dtau_line, rho
!!$  
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call set_bias_ttauri(child, grid, lambda0, outVec)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          if (thisOctal%inflow(subcell)) then
!!$             d = thisOctal%subcellsize 
!!$
!!$             rVec = subcellCentre(thisOctal,subcell)
!!$             r = modulus(rvec)
!!$             if (r /= 0.0d0) then
!!$                rhat = rvec/r
!!$             else
!!$                rhat = VECTOR(0.0d0, 0.0d0, 1.0d0)
!!$             end if
!!$
!!$             if (thisOctal%threed) then
!!$                dV = d*d*d
!!$             else
!!$                dv = 2.0_db*pi*d*d*rVec%x
!!$             endif
!!$
!!$             nu0  = cSpeed_dbl / dble(lambda0*angstromtocm)
!!$
!!$             tauSob = thisOctal%chiline(subcell)  / nu0
!!$
!!$!             ! Find the average tauSob over all azimuth positions
!!$!             nsample = 40
!!$!             tauSob_sum = 0.0d0
!!$!             do i = 1, nsample
!!$!                phi = 2.0d0*PIdouble*dble(i-1)/dble(nsample-1)
!!$!                theta = ACOS(rhat%z); sin_theta = SIN(theta)
!!$!                rvec_sample = VECTOR(r*cos(phi)*sin_theta, r*sin(phi)*sin_theta,  rvec%z)
!!$!                tauSob_tmp = tauSob / amrGridDirectionalDeriv(grid, rvec_sample, outVec, &
!!$!                  startOctal=thisOctal)
!!$!                tauSob_sum = tauSob_sum + tauSob_tmp
!!$!             end do
!!$!             tauSob = tauSob_sum/dble(nsample)
!!$
!!$             ! in obs direction, but works only for i=0 case.
!!$!             tauSob = tauSob / amrGridDirectionalDeriv(grid, rvec, outVec, &
!!$!                  startOctal=thisOctal)
!!$             ! in a radial direction
!!$             tauSob = tauSob / amrGridDirectionalDeriv(grid, rvec, rhat, &
!!$                  startOctal=thisOctal)
!!$           
!!$             if (tauSob < 0.01) then
!!$                escProb = 1.0d0-tauSob*0.5d0*(1.0d0 -   &
!!$                     tauSob/3.0d0*(1. - tauSob*0.25d0*(1.0d0 - 0.20d0*tauSob)))
!!$             else if (tauSob < 15.) then
!!$                escProb = (1.0d0-exp(-tauSob))/tauSob
!!$             else
!!$                escProb = 1.d0/tauSob
!!$             end if
!!$             escProb = max(escProb, 1.d-5)
!!$!             escProb = min(1.d-2,max(escProb, 1.d-5))
!!$
!!$
!!$             dtau_cont = d*(thisOctal%kappaAbs(subcell,1) + thisOctal%kappaSca(subcell,1))
!!$             dtau_line = d*(thisOctal%chiline(subcell))  / nu0
!!$             !
!!$             rho = thisOctal%rho(subcell)
!!$!             thisOctal%biasCont3D(subcell) = MAX(dV, 1.d-7) ! Limits the minimum value
!!$!             thisOctal%biasLine3D(subcell) = dV*thisOctal%biasCont3D(subcell)
!!$             thisOctal%biasCont3D(subcell) = MAX(EXP(-dtau_cont), 1.d-7) ! Limits the minimum value
!!$!             thisOctal%biasCont3D(subcell) = 1.0d0
!!$!             thisOctal%biasLine3D(subcell) = thisOctal%rho(subcell)*thisOctal%biasCont3D(subcell)
!!$!             thisOctal%biasLine3D(subcell) = escProb*thisOctal%biasCont3D(subcell)
!!$             thisOctal%biasLine3D(subcell) = 1.0d0/SQRT(rho) * thisOctal%biasCont3D(subcell)
!!$!             thisOctal%biasLine3D(subcell) = 1.0d0/(rho**0.25)
!!$!             thisOctal%biasLine3D(subcell) = EXP(-dtau_line)*thisOctal%biasCont3D(subcell)
!!$!             thisOctal%biasLine3D(subcell) = EXP(-dtau_line*d*rVec%x)*thisOctal%biasCont3D(subcell)
!!$
!!$          else  ! this subcell is not "inFlow"
!!$             thisOctal%biasCont3D(subcell) = 1.0d-150
!!$             thisOctal%biasLine3D(subcell) = 1.0d-150
!!$          end if
!!$          ! just in case ....
!!$          thisOctal%biasCont3D(subcell) = MAX(thisOctal%biasCont3D(subcell), 1.0d-150)
!!$          thisOctal%biasLine3D(subcell) = MAX(thisOctal%biasLine3D(subcell), 1.0d-150)       
!!$       endif ! if (thisOctal%hasChild(subcell)) then
!!$    enddo
!!$  end subroutine set_bias_ttauri
!!$
!!$
!!$  recursive subroutine set_bias_ttauri2(thisOctal, grid)
!!$    use input_variables, only: TTauriRinner, TTauriRouter
!!$    type(gridtype) :: grid
!!$    type(octal), pointer   :: thisOctal
!!$    type(octal), pointer  :: child 
!!$    integer :: subcell, i
!!$    real(double) :: r
!!$    type(vector)  :: rvec
!!$    real(double)::  rM, theta,  h, rM_center
!!$  
!!$    do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call set_bias_ttauri2(child, grid)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          if (thisOctal%inflow(subcell)) then
!!$             
!!$             rVec = subcellCentre(thisOctal,subcell)
!!$             r = modulus(rvec)
!!$             
!!$             if (r/=0.0d0) then
!!$                theta = ACOS(MIN(ABS(rVec%z/r),0.998_oc))
!!$             else
!!$                theta=0.01
!!$             end if
!!$
!!$             if (theta == 0.0d0) theta = 0.01
!!$             rM  = r / SIN(theta)**2
!!$
!!$             ! bias more toward the edges of the accreation stream
!!$             h = 0.5d0*(TTauriRouter - TTauriRinner)*1.0d-10           ! [10^10cm]   
!!$             rM_center = 0.5d0*(TTauriRouter + TTauriRinner)*1.0d-10   ![10^10cm]   
!!$             
!!$             thisOctal%biasCont3D(subcell) = 1.0d0  ! no bias for contiuum
!!$!             thisOctal%biasLine3D(subcell) = 1.0d0/(dtau_line*d*rvec%x)
!!$!             thisOctal%biasLine3D(subcell) = 1.0d0/thisOctal%rho(subcell)
!!$!             thisOctal%biasLine3D(subcell) = EXP(10.0d0*ABS(rM - rM_center)/h)
!!$!             if (ABS(rM - rM_center)/h > 0.80d0) then
!!$             if (ABS(rM - rM_center)/h > 0.40d0) then
!!$                thisOctal%biasLine3D(subcell) = 1.0d5
!!$             else
!!$                thisOctal%biasLine3D(subcell) = 1.0d-150
!!$             end if
!!$          else  ! this subcell is not "inFlow"
!!$             thisOctal%biasCont3D(subcell) = 1.0d-150
!!$             thisOctal%biasLine3D(subcell) = 1.0d-150
!!$          end if
!!$          ! just in case ....
!!$          thisOctal%biasCont3D(subcell) = MAX(thisOctal%biasCont3D(subcell), 1.0d-150)
!!$          thisOctal%biasLine3D(subcell) = MAX(thisOctal%biasLine3D(subcell), 1.0d-150)       
!!$       endif ! if (thisOctal%hasChild(subcell)) then
!!$    enddo
!!$  end subroutine set_bias_ttauri2
!!$
!!$
!!$
!!$
!!$  recursive subroutine setBiasPpdisk(thisOctal, grid)
!!$!  use input_variables, only: rInner
!!$  type(gridtype) :: grid
!!$  type(VECTOR) :: rVec
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$!  real(double) :: r
!!$  integer :: subcell, i
!!$  
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call setBiasPpdisk(child, grid)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$!          r = modulus(subcellcentre(thisOctal, subcell)) / rInner
!!$!          thisOctal%biasCont3D(subcell) =  sqrt(r)
!!$!!          if ((r > 1.).and.(r < 1.05)) then
!!$!          if ((r > 1.).and.(r < 2.)) then
!!$!             thisOctal%biasCont3D(subcell) =  thisOctal%biasCont3D(subcell) * 10.d0
!!$!          endif
!!$
!!$          rVec = subcellcentre(thisOctal, subcell)
!!$          thisOctal%biasCont3D(subcell) = rVec%x**2.
!!$
!!$
!!$          if (thisOctal%diffusionApprox(subcell)) then
!!$             thisOctal%biasCont3D(subcell) = thisOctal%biasCont3D(subcell) * 1.e-2
!!$          endif
!!$       endif
!!$    enddo
!!$  end subroutine setBiasPpdisk
!!$
!!$
!!$  recursive subroutine massInAnn(thisOctal, r1, r2, mass)
!!$  type(VECTOR) :: rVec
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  real :: mass, r1, r2, r
!!$  integer :: subcell, i
!!$  
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call massInAnn(child,  r1, r2, mass)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          rVec = subcellCentre(thisOctal, subcell)
!!$          r  = sqrt(rVec%x**2 + rVec%y**2)
!!$          if ( (r > r1).and.(r < r2) ) then
!!$             mass = mass + thisOctal%rho(subcell)*thisOctal%subCellSize**3 
!!$          endif
!!$       endif
!!$    enddo
!!$  end subroutine massInAnn
!!$

!!$
!!$  subroutine setDiscPhotosphereBias(grid, iLam)
!!$    use input_variables, only: rgap
!!$    type(gridtype) :: grid
!!$    integer :: iLam
!!$    real :: xAxis(1:100000)
!!$    real :: tau(1:100000)
!!$    integer :: nx, i, i1, i2, iGap
!!$    type(OCTAL), pointer :: thisOctal
!!$    integer :: subcell
!!$    real(double) kappaSca, kappaAbs
!!$    type(VECTOR) :: octVec
!!$    kappaAbs = 0.d0; kappaSca = 0.d0
!!$    write(*,*) "Setting disc photosphere bias..."
!!$    nx = 0
!!$    xAxis = 0
!!$    call getxValuesAMR(grid%octreeRoot, nx, xAxis)
!!$    call stripSimilarValues(xAxis,nx,real(1.d-5*grid%halfSmallestSubcell))
!!$
!!$    tau(1) = 0.
!!$    do i = 2, nx
!!$       octVec = VECTOR(xAxis(i), 0.d0, 0.d0)
!!$       CALL findSubcellTD(octVec,grid%octreeRoot,thisOctal,subcell)
!!$       call returnKappa(grid, thisOctal, subcell, ilambda=ilam,&
!!$            kappaSca=kappaSca, kappaAbs=kappaAbs)
!!$       tau(i) = tau(i-1) + (xAxis(i)-xAxis(i-1)) * (kappaAbs + kappaSca)
!!$       call setVerticalBias(grid, xAxis(i), 0., ilam, thisTau=tau(i))
!!$       if (tau(i) > 5.) then
!!$          i1 = i
!!$          exit
!!$       endif
!!$    enddo
!!$       
!!$
!!$    call locate(xAxis, nx, rGap*real(autocm)/1.e10, iGap)
!!$
!!$    tau(i+1) = 0.
!!$    do i = iGap, 1, -1
!!$       octVec = VECTOR(xAxis(i), 0.d0, 0.d0)
!!$       CALL findSubcellTD(octVec,grid%octreeRoot,thisOctal,subcell)
!!$       call returnKappa(grid, thisOctal, subcell, ilambda=ilam,&
!!$            kappaSca=kappaSca, kappaAbs=kappaAbs)
!!$       tau(i) = tau(i+1) + (xAxis(i+1)-xAxis(i)) * (kappaAbs + kappaSca)
!!$       call setVerticalBias(grid, xAxis(i), 0., ilam, thisTau=tau(i))
!!$       if (tau(i) > 5.) then
!!$          i2 = i
!!$          exit
!!$       endif
!!$    enddo
!!$       
!!$
!!$    write(*,*) "range",xAxis(i1),xAxis(i2)
!!$    do i = i1, i2
!!$       call setVerticalBias(grid, xAxis(i), 0., ilam)
!!$    enddo
!!$
!!$    tau(i-1) = 0.
!!$    do i = iGap, nx
!!$       octVec = VECTOR(xAxis(i), 0.d0, 0.d0)
!!$       CALL findSubcellTD(octVec,grid%octreeRoot,thisOctal,subcell)
!!$       tau(i) = tau(i-1) + (xAxis(i)-xAxis(i-1)) * (kappaAbs + kappaSca)
!!$       call setVerticalBias(grid, xAxis(i), 0., ilam, thisTau=tau(i))
!!$       if (tau(i) > 5.) then
!!$          i1 = i
!!$          exit
!!$       endif
!!$    enddo
!!$       
!!$    write(*,*) "from ",xAxis(i1)
!!$    do i = i1 , nx
!!$       call setVerticalBias(grid, xAxis(i), 0., ilam)
!!$    enddo
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$    write(*,*) "Done."
!!$
!!$  end subroutine setDiscPhotosphereBias
!!$
!!$    

!!$  subroutine interpHydroProperties(grid, thisOctal, subcell)
!!$    type(GRIDTYPE) :: grid
!!$    type(OCTAL), pointer :: thisOctal, parentOctal
!!$    type(OCTAL),pointer :: neighbourOctalMinus, neighbourOctalPlus
!!$    integer :: subcell, parentSubcell
!!$    type(VECTOR) :: direction, cellCentre, rVec, parentSubcellCentre
!!$    real(double) :: xMinus, xMid, xPlus
!!$    integer :: neighbourSubcellPlus, neighbourSubcellMinus
!!$
!!$    parentOctal => thisOctal%parent
!!$    parentSubcell = thisOctal%parentSubcell
!!$    cellCentre = subcellCentre(thisOctal, subcell)
!!$    xMid = cellCentre%x
!!$    parentSubcellCentre = subcellCentre(parentOctal, parentSubcell)
!!$    direction = VECTOR(1.d0, 0.d0, 0.d0)
!!$
!!$    neighbourOctalMinus => thisOctal
!!$    rVec = parentSubcellCentre - direction*(parentoctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
!!$    call findSubcellLocal(rVec, neighbourOctalMinus, neighbourSubcellMinus)
!!$    rVec = subcellCentre(neighbourOctalminus, neighbourSubcellMinus)
!!$    xMinus = rvec%x
!!$    neighbourOctalPlus => thisOctal
!!$    rVec = parentSubcellCentre - direction*(parentoctal%subcellSize/2.d0 + 0.01d0*grid%halfSmallestSubcell)
!!$    call findSubcellLocal(rVec, neighbourOctalPlus, neighbourSubcellPlus)
!!$    rVec = subcellCentre(neighbourOctalPlus, neighbourSubcellPlus)
!!$    xPlus = rvec%x
!!$    
!!$
!!$  end subroutine interpHydroProperties

!!$
!!$  recursive subroutine myTauSplit(thisOctal, grid, converged, inheritProps, interpProps)
!!$    use input_variables, only : maxDepthAMR
!!$    type(gridtype) :: grid
!!$    type(octal), pointer   :: thisOctal
!!$    type(octal), pointer  :: child
!!$    logical, optional :: inheritProps, interpProps
!!$    integer :: subcell, i
!!$    logical :: converged
!!$    real(double) :: kabs, ksca
!!$    real :: thisTau
!!$    logical :: split
!!$    logical, save :: firsttime = .true.
!!$    character(len=30) :: message
!!$    kabs = 0.d0; ksca = 0.d0
!!$
!!$    do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call myTauSplit(child, grid, converged, inheritProps, interpProps)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$
!!$          split = .true.
!!$
!!$          call returnKappa(grid, thisOctal, subcell, rosselandKappa=kabs)
!!$          thisTau = thisOctal%subcellsize * kabs * thisOctal%rho(subcell) * 1.e10
!!$
!!$          if (thisOctal%nDepth == maxDepthamr) then
!!$             split = .false.
!!$             if (firstTime) then
!!$                write(message,'(a,i3)') "AMR cell depth capped at: ",maxDepthamr
!!$                call writeWarning(message)
!!$                firstTime = .false.
!!$             endif
!!$          endif
!!$
!!$
!!$
!!$          if (split.and.(thisTau > 2.).and.(thisOctal%chiLine(subcell) < 100.)) then
!!$             write(*,*) "splitting cell with ",thisTau, " at depth ",thisOctal%chiline(subcell)
!!$             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
!!$                  inherit=inheritProps, interp=interpProps)
!!$             converged = .false.
!!$             return
!!$          endif
!!$       endif
!!$    enddo
!!$
!!$  end subroutine myTauSplit
!!$
!!$  recursive subroutine massSplit(thisOctal, grid, converged, inheritProps, interpProps, rLimit)
!!$    use input_variables, only : maxDepthAMR, limitScalar
!!$    type(gridtype) :: grid
!!$    type(octal), pointer   :: thisOctal
!!$    type(octal), pointer  :: child
!!$    real(double), optional :: rLimit
!!$    logical, optional :: inheritProps, interpProps
!!$    integer :: subcell, i
!!$    real(double) :: cellMass, dv
!!$    logical :: converged
!!$    logical :: split
!!$    logical, save :: firsttime = .true.
!!$    character(len=30) :: message
!!$
!!$    do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call massSplit(child, grid, converged, inheritProps, interpProps, rLimit)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$
!!$          split = .false.
!!$
!!$          dv = cellVolume(thisOctal, subcell)*1.d30
!!$          cellMass = dv * thisOctal%rho(subcell) 
!!$
!!$          if (thisOctal%nDepth == maxDepthamr) then
!!$             split = .false.
!!$             if (firstTime) then
!!$                write(message,'(a,i3)') "massSplit: AMR cell depth capped at: ",maxDepthamr
!!$                call writeWarning(message)
!!$                firstTime = .false.
!!$             endif
!!$          endif
!!$
!!$          if (cellMass > limitScalar) split = .true.
!!$          if (PRESENT(rlimit)) then
!!$             if (modulus(subcellCentre(thisOctal,subcell)) > rLimit) split = .false.
!!$          endif
!!$          if (split) then
!!$!             write(*,*) "splitting cell with mass ",cellmass
!!$             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
!!$                  inherit=inheritProps, interp=interpProps)
!!$             converged = .false.
!!$             return
!!$          endif
!!$       endif
!!$    enddo
!!$
!!$  end subroutine massSplit
!!$
!!$

!!$  recursive subroutine stripMahdavi(thisOctal)
!!$    use magnetic_mod, only : inflowMahdavi
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  integer :: subcell, i
!!$  
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call stripMahdavi(child)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$          thisOctal%inflow(subcell) = .true.
!!$          if (.not.inflowMahdavi(1.d10*subcellCentre(thisOctal,subcell))) then
!!$             thisOctal%rho(subcell) = 1.e-25
!!$             thisOctal%inflow(subcell) = .false.
!!$          endif
!!$       endif
!!$    enddo
!!$  end subroutine stripMahdavi
!!$

!!$  recursive subroutine fillVelocityCornersMahdavi(thisOctal,grid)
!!$    use magnetic_mod, only : velocityMahdavi
!!$    use input_variables, only : vturb
!!$    type(GRIDTYPE) :: grid
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  integer :: subcell, i
!!$  
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call fillVelocityCornersMahdavi(child, grid)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$
!!$          if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = vturb
!!$
!!$          IF ((thisoctal%threed).and.(subcell == 8)) &
!!$               CALL fillVelocityCorners(thisOctal, velocityMahdavi)
!!$          
!!$          IF ((thisoctal%twod).and.(subcell == 4)) &
!!$               CALL fillVelocityCorners(thisOctal, Velocitymahdavi)
!!$       endif
!!$    enddo
!!$  end subroutine fillVelocityCornersMahdavi
!!$
!!$  recursive subroutine convertToDensity(thisOctal)
!!$  type(octal), pointer   :: thisOctal
!!$  type(octal), pointer  :: child 
!!$  integer :: subcell, i
!!$  
!!$  do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call convertToDensity(child)
!!$                exit
!!$             end if
!!$          end do
!!$       else
!!$
!!$          thisOctal%rho(subcell) = thisOctal%rho(subcell) &
!!$               / (1.d30*cellVolume(thisOctal,subcell))
!!$
!!$       endif
!!$    enddo
!!$  end subroutine convertToDensity

!!$  SUBROUTINE smoothAMRgrid(grid,factor, inheritProps, interpProps, &
!!$       romData)
!!$    ! checks whether each octal's neighbours are much bigger than it, 
!!$    !   if so, makes the neighbours smaller.
!!$
!!$    TYPE(gridtype), INTENT(INOUT) :: grid 
!!$    REAL, INTENT(IN)              :: factor
!!$    LOGICAL, INTENT(IN),optional  :: inheritProps
!!$    logical, intent(in), optional :: interpProps 
!!$    TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
!!$    
!!$    LOGICAL :: gridConverged
!!$
!!$    CALL setAllUnchanged(grid%octreeRoot)
!!$
!!$    DO 
!!$      gridConverged = .TRUE.
!!$      CALL smoothAMRgridPrivate(grid%octreeRoot,grid,gridConverged,romData)
!!$      IF ( gridConverged ) EXIT 
!!$    END DO
!!$    
!!$  CONTAINS
!!$    
!!$    RECURSIVE SUBROUTINE smoothAMRgridPrivate(thisOctal,grid,gridConverged, &
!!$         romData)
!!$
!!$!      TYPE(octal), INTENT(INOUT), TARGET :: thisOctal
!!$      TYPE(octal), INTENT(INOUT) :: thisOctal
!!$      TYPE(gridtype), INTENT(INOUT   ) :: grid 
!!$      LOGICAL, INTENT(INOUT)               :: gridConverged
!!$      TYPE(romanova), optional, INTENT(IN)   :: romDATA  ! used for "romanova" geometry
!!$
!!$      INTEGER              :: i
!!$      REAL(oct) :: halfSmallestSubcell
!!$      REAL(oct) :: offset
!!$      TYPE(octal), POINTER :: neighbour
!!$      TYPE(vector), ALLOCATABLE, DIMENSION(:) :: locator
!!$      TYPE(vector) :: aHat
!!$      INTEGER              :: subcell
!!$      INTEGER              :: nLocator ! number of locators (4 for twoD, 6 for threed)
!!$
!!$      ! we will find the coordinates of a point that lies outside the current
!!$      !   octal. we then compare the size of the cell that contains that point
!!$      !   with the size of the current cell, if it is bigger be more than a 
!!$      !   factor of 'factor', we subdivide the neighbouring cell.
!!$      ! we do this in each of six directions
!!$
!!$      ! we do not have to test the other subcells in the current octal because
!!$      !   they can be smaller than the any of the other subcells, but they
!!$      !   cannot be *bigger*. this saves some time.
!!$
!!$      if (thisOctal%threed) then
!!$         if (.not.thisOctal%cylindrical) then
!!$            nlocator = 6
!!$         else
!!$            nlocator = 4
!!$         endif
!!$      else
!!$        nlocator = 4
!!$      endif
!!$
!!$      ALLOCATE(locator(nlocator))
!!$
!!$      ! we find points which are outside the current octal by a distance
!!$      !   equivalent to half the size of the tree's smallest subcell.
!!$      halfSmallestSubcell = grid%halfSmallestSubcell
!!$
!!$      ! we also add a slight offset to our test positions to avoid testing 
!!$      !   at cell boundaries.
!!$      offset = halfSmallestSubcell / 2.0_oc
!!$
!!$      IF ( thisOctal%threed ) THEN
!!$        locator(:) = thisOctal%centre + ( offset * vector(1.0_oc,1.0_oc,1.0_oc) )
!!$      ELSE
!!$        locator(:) = thisOctal%centre + ( offset * vector(1.0_oc,0.0_oc,1.0_oc) )
!!$      END IF
!!$
!!$      IF ( thisOctal%threeD ) THEN
!!$         if (.not.thisOctal%cylindrical) then
!!$            locator(1)%x = thisOctal%centre%x + thisOctal%subcellSize + halfSmallestSubcell
!!$            locator(2)%y = thisOctal%centre%y + thisOctal%subcellSize + halfSmallestSubcell
!!$            locator(3)%z = thisOctal%centre%z + thisOctal%subcellSize + halfSmallestSubcell
!!$            locator(4)%x = thisOctal%centre%x - thisOctal%subcellSize - halfSmallestSubcell
!!$            locator(5)%y = thisOctal%centre%y - thisOctal%subcellSize - halfSmallestSubcell
!!$            locator(6)%z = thisOctal%centre%z - thisOctal%subcellSize - halfSmallestSubcell
!!$         else
!!$            locator(:) = thisOctal%centre
!!$            locator(1) = locator(1) + (thisOctal%subcellSize + halfSmallestSubcell) * zHat
!!$            locator(2) = locator(2) - (thisOctal%subcellSize + halfSmallestSubcell) * zHat
!!$            aHat = VECTOR(thisOctal%centre%x,thisOctal%centre%y,0.d0)
!!$            call normalize(aHat)
!!$            locator(3) = locator(3) + (thisOctal%subcellSize + halfSmallestSubcell) * aHat
!!$            locator(4) = locator(4) - (thisOctal%subcellSize + halfSmallestSubcell) * aHat
!!$         endif
!!$      ELSE
!!$        locator(1)%x = thisOctal%centre%x + thisOctal%subcellSize + halfSmallestSubcell
!!$        locator(2)%z = thisOctal%centre%z + thisOctal%subcellSize + halfSmallestSubcell
!!$        locator(3)%x = thisOctal%centre%x - thisOctal%subcellSize - halfSmallestSubcell
!!$        locator(4)%z = thisOctal%centre%z - thisOctal%subcellSize - halfSmallestSubcell
!!$      ENDIF
!!$      
!!$      DO i = 1, nLocator, 1
!!$        IF ( inOctal(grid%octreeRoot,locator(i)) ) THEN
!!$          CALL findSubcellTD(locator(i),grid%octreeRoot,neighbour,subcell)
!!$          !neighbour => thisOctal
!!$          !CALL findSubcellLocal(locator(i),neighbour,subcell)
!!$
!!$          IF ( neighbour%subcellSize > (factor * thisOctal%subcellSize) ) THEN
!!$            IF ( neighbour%hasChild(subcell) ) THEN 
!!$              PRINT *, "neighbour already has child. (B)"
!!$              !STOP
!!$              do ;end do
!!$            END IF
!!$              neighbour%changed(1:neighbour%maxChildren) = .TRUE.
!!$              call addNewChild(neighbour, subcell, grid, adjustGridInfo=.TRUE.,  &
!!$                               inherit=inheritProps, interp=interpProps, romData=romData)
!!$              gridConverged = .FALSE.
!!$          ENDIF
!!$        END IF
!!$          
!!$! force return until algorithm is fixed
!!$IF ( .NOT. gridConverged ) RETURN
!!$          
!!$      END DO
!!$
!!$      ! call this subroutine recursively on each of any children.
!!$      IF ( (.NOT. ANY(thisOctal%Changed)) .AND. thisOctal%nChildren > 0 ) THEN
!!$        DO i = 1, thisOctal%nChildren, 1 
!!$          CALL smoothAMRgridPrivate(thisOctal%child(i),grid,gridConverged,romData=romData)
!!$          !CALL checkAMRgrid(grid,checkNoctals=.FALSE.)                                  
!!$        
!!$! force return until algorithm is fixed
!!$IF ( .NOT. gridConverged ) RETURN
!!$        END DO
!!$      END IF
!!$
!!$      DEALLOCATE(locator)
!!$          
!!$    END SUBROUTINE smoothAMRgridPrivate
!!$  
!!$  END SUBROUTINE smoothAMRgrid

  !
  ! Recursively turn off the magnetosphere
  !
!!$  RECURSIVE SUBROUTINE turn_off_magnetosphere(thisOctal,grid, Rmax)    
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    TYPE(octal), POINTER   :: thisOctal
!!$    TYPE(gridtype)         :: grid
!!$    ! The outer radius of the magnetosphere.
!!$    real(double), intent(in) :: Rmax  ! [10^10cm] 
!!$    
!!$    TYPE(octal), POINTER   :: pChild
!!$    INTEGER :: subcell, n
!!$
!!$    if (thisOctal%threeD) then
!!$       n = 8
!!$    else
!!$       n =4
!!$    end if
!!$    
!!$
!!$    do subcell = 1, n
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! just decdend the tree branch
!!$          pChild => thisOctal%child(subcell)
!!$          CALL turn_off_magnetosphere(pChild,grid,Rmax)
!!$       else
!!$          ! turnning it off 
!!$!          if (TTauriInFlow(thisOctal%centre, grid)) then          
!!$          if ( (modulus(thisOctal%centre) <= Rmax .and. abs(thisOctal%centre%z)>5.0d0)  &
!!$               .or.  &
!!$               (modulus(thisOctal%centre) <= Rmax*1.01 .and. abs(thisOctal%centre%z)<5.0d0) ) then 
!!$             thisOctal%inFlow(subcell) = .false.
!!$             thisOctal%kappaAbs(subcell,:) = 1.0e-30
!!$             thisOctal%kappaSca(subcell,:) = 1.0e-30          
!!$             thisOctal%temperature(subcell) = 6500.0
!!$             thisOctal%biasCont3D(subcell) = 1.0e-30  
!!$             thisOctal%biasLine3D(subcell) = 1.0e-30  
!!$             thisOctal%etaLine(subcell) = 1.e-30
!!$             thisOctal%etaCont(subcell) = 1.e-30
!!$             thisOctal%cornerVelocity = vector(1.e-30,1.e-30,1.e-30)
!!$             thisOctal%velocity = vector(1.e-30,1.e-30,1.e-30)
!!$             thisOctal%rho(subcell) = 1.0d-19
!!$          end if
!!$       end if
!!$    end do
!!$    
!!$  END SUBROUTINE turn_off_magnetosphere
!!$

!!$  subroutine hydroWarpFitSplines(grid)
!!$    use input_variables, only: rOuter
!!$
!!$    type(GRIDTYPE) :: grid
!!$    real :: xAxis(1000000)
!!$    real(double) :: zAxis(10000), rho(10000), subcellsize(10000)
!!$    integer :: nx, nz
!!$    real :: xPos, yPos
!!$    integer :: i, j
!!$
!!$    real(double) :: newrho(10000), newzAxis(10000), newrhodd(10000)
!!$    integer :: newnz
!!$    zAxis = 0
!!$    nx = 0
!!$    nz = 0
!!$    newrhodd = 0.d0 ; rho = 0.d0; subcellSize = 0.d0
!!$    call getxValuesAMR(grid%octreeRoot, nx, xAxis)
!!$    call stripSimilarValues(xAxis,nx,real(1.d-5*grid%halfSmallestSubcell))
!!$    do while (xAxis(nx) > rOuter)
!!$      nx = nx - 1
!!$    end do
!!$    xAxis(1:nx) = xAxis(1:nx) + 1.d-5*grid%halfSmallestSubcell
!!$
!!$    allocate(grid%hydroSplines(nx))
!!$
!!$    yPos = 0.
!!$    do i = 1, nx
!!$       xPos = xAxis(i)
!!$!       if ((xPos > grid%rInner).and.(xPos < grid%rOuter)) then
!!$          ! Read in a complete density run for this x position
!!$          call getDensityRun(grid, zAxis, subcellsize, rho, xPos, yPos, nz, -1.)
!!$          do j = 1, nz
!!$            newRho(j) = rho(nz+1-j)
!!$            newzAxis(j) = -1.*zAxis(nz+1-j)
!!$          end do
!!$          newnz = nz
!!$          call getDensityRun(grid, zAxis, subcellsize, rho, xPos, yPos, nz, +1.)
!!$          newRho(newnz+1:newnz+nz) = rho(1:nz)
!!$          newzAxis(newnz+1:newnz+nz) = zAxis(1:nz)
!!$          newnz = newnz + nz
!!$
!!$          ! We use the log of the densities to avoid numerical funnies in the
!!$          ! spline fitting.
!!$          do j = 1, newnz
!!$            newRho(j) = log10(newRho(j))
!!$          enddo
!!$
!!$          ! Calculate the spline function
!!$          call spline(newzAxis, newRho, newnz, 1.d30, 1.d30, newRhoDD)
!!$
!!$          ! Store the spline in the appropriate array
!!$          allocate(grid%hydroSplines(i)%z(newnz))
!!$          allocate(grid%hydroSplines(i)%rho(newnz))
!!$          allocate(grid%hydroSplines(i)%rhoDD(newnz))
!!$          grid%hydroSplines(i)%x = xPos
!!$          grid%hydroSplines(i)%nz = newnz
!!$          grid%hydroSplines(i)%z(1:newnz) = newzAxis(1:newnz)
!!$          grid%hydroSplines(i)%rho(1:newnz) = newRho(1:newnz)
!!$          grid%hydroSplines(i)%rhoDD(1:newnz) = newRhoDD(1:newnz)
!!$          ! Set the scaleheight
!!$          call hydroWarpScaleHeight(grid%hydroSplines(i))
!!$!       endif
!!$    enddo
!!$  end subroutine hydroWarpFitSplines

!!$  recursive subroutine splitGridOnStream(thisOctal, thisStream, grid, converged)
!!$
!!$    type(GRIDTYPE) :: grid
!!$    type(OCTAL), pointer :: thisOctal, child
!!$    type(VECTOR) :: rVec !, corner
!!$    integer :: i, j, subcell
!!$    logical :: converged, converged_tmp
!!$    logical :: split
!!$    type(STREAMTYPE), pointer :: thisStream
!!$    real(double), parameter :: fac = 2.d0
!!$    converged = .true.
!!$    converged_tmp=.true.
!!$
!!$    do subcell = 1, thisOctal%maxChildren
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          children : do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call splitGridOnStream(child, thisStream, grid, converged)
!!$                if (.not.converged) converged_tmp = converged
!!$                exit children
!!$             end if
!!$          end do children
!!$          if (.not.converged_tmp) converged=converged_tmp
!!$       else
!!$
!!$          rVec = subcellCentre(thisOctal, subcell)
!!$
!!$          stream : do j = 1, thisStream%nSamples
!!$             
!!$             
!!$             split = .false.
!!$             !corner = closestCorner(thisOctal, subcell, thisStream%position(j))
!!$             
!!$             if ((inSubcell(thisOctal, subcell, thisStream%position(j))).and. &
!!$                  (thisOctal%subcellSize > thisStream%streamRadius(j)/fac) ) then 
!!$                split = .true.
!!$             else
!!$                
!!$                if ((dist_from_closestCorner(thisOctal, subcell, thisStream%position(j)) < &
!!$                     thisStream%StreamRadius(j)).and. &
!!$                     (thisOctal%subcellSize > thisStream%streamRadius(j)/fac) ) split = .true.
!!$             endif
!!$             
!!$             if (split) then
!!$
!!$!                write(*,*) "adding child",thisStream%streamRadius(j),thisOctal%subcellSize
!!$                call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
!!$                     inherit=.false., interp=.false.)
!!$                
!!$                converged = .false.
!!$                exit stream
!!$             endif
!!$          enddo stream 
!!$          if (.not.converged) exit 
!!$       end if
!!$    enddo
!!$
!!$    return
!!$    
!!$  end subroutine splitGridOnStream
!!$
!!$  recursive subroutine splitGridOnStream2(thisOctal, thisStream, grid, childrenAdded)
!!$
!!$    type(GRIDTYPE) :: grid
!!$    type(OCTAL), pointer :: thisOctal, child
!!$    type(VECTOR) :: rVec
!!$    integer :: i, j, subcell
!!$    logical :: childrenAdded
!!$    logical :: split
!!$    type(STREAMTYPE) :: thisStream
!!$    real(double), parameter :: fac = 2.d0
!!$
!!$    subcell = 1
!!$    do while (subcell <= thisOctal%maxChildren)
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          children : do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call splitGridOnStream2(child, thisStream, grid, childrenAdded)
!!$                exit children
!!$             end if
!!$          end do children
!!$       else
!!$
!!$
!!$          childrenAdded = .false.
!!$          rVec = subcellCentre(thisOctal, subcell)
!!$
!!$          stream : do j = 1, thisStream%nSamples
!!$             
!!$             
!!$             split = .false.
!!$             
!!$
!!$
!!$             if ((inSubcell(thisOctal, subcell, thisStream%position(j))).and. &
!!$                  (thisOctal%subcellSize > thisStream%streamRadius(j)/fac) ) then 
!!$                split = .true.
!!$             else
!!$                
!!$                if ((dist_from_closestEdge(thisOctal, subcell, thisStream%position(j)) < &
!!$                     thisStream%StreamRadius(j)).and. &
!!$                     (thisOctal%subcellSize > thisStream%streamRadius(j)/fac) ) split = .true.
!!$             endif
!!$             
!!$             if (split) then
!!$
!!$                call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
!!$                     inherit=.false., interp=.false.)
!!$                childrenAdded = .true.
!!$                subcell = subcell - 1
!!$                exit stream
!!$             endif
!!$          enddo stream
!!$       end if
!!$       subcell = subcell + 1
!!$    enddo
!!$
!!$  end subroutine splitGridOnStream2
!!$
!!$  recursive subroutine splitGridOnStream3(thisOctal, grid, thisStream)
!!$    type(GRIDTYPE) :: grid
!!$    type(OCTAL), pointer :: thisOctal, child
!!$    integer :: i, subcell
!!$    logical :: split
!!$    type(STREAMTYPE) :: thisStream
!!$    type(STREAMTYPE) :: octalStream, subcellStream
!!$!    real(double) :: side, rStar
!!$    logical :: splitAz
!!$
!!$    call createOctalStream(thisOctal, thisStream, octalStream)
!!$
!!$    subcell = 1
!!$    do while (subcell <= thisOctal%maxChildren)
!!$       if (thisOctal%hasChild(subcell)) then
!!$          ! find the child
!!$          children : do i = 1, thisOctal%nChildren, 1
!!$             if (thisOctal%indexChild(i) == subcell) then
!!$                child => thisOctal%child(i)
!!$                call splitGridOnStream3(child, grid, octalStream)
!!$                exit children
!!$             end if
!!$          end do children
!!$       else
!!$
!!$          split = .false.
!!$          call createSubcellStream(thisOctal, subcell, octalStream, subcellStream)
!!$!          if (subcellStream%nSamples > 1) split = .true.
!!$ 
!!$
!!$          split = .true.
!!$!          rStar = 2.d0*rSol/1.d10
!!$!          side = 0.1d0 * (thisOctal%r/rStar)**3
!!$!          if (thisOctal%subcellSize > side) split = .true.
!!$
!!$          if (thisOctal%dphi*radtodeg > 15.) then
!!$             splitAz = .true.
!!$          else
!!$             splitAz = .false.
!!$          endif
!!$
!!$          if (thisOctal%nDepth > 7) split = .false.
!!$
!!$          if (subcellstream%nsamples == 0) split = .false.
!!$
!!$          if (split) then
!!$             call addNewChild(thisOctal,subcell,grid,adjustGridInfo=.TRUE., &
!!$                     inherit=.false., interp=.false., stream=subcellStream, splitAzimuthally=splitAz)
!!$             subcell = subcell - 1
!!$          endif
!!$       end if
!!$       subcell = subcell + 1
!!$       call freeStream(subcellStream)
!!$    enddo
!!$ 
!!$
!!$    call freeStream(octalStream)
!!$  end subroutine splitGridOnStream3

!!$  subroutine readStreams(thisStream, nStream, filename)
!!$    use input_variables, only : scaleFlowRho, flowrhoScale, rcore
!!$    type(STREAMTYPE) :: thisStream(:)
!!$    character(len=*) :: filename
!!$    integer :: nstream
!!$    integer :: i, j
!!$    real(double) :: r, theta, phi, v, rho
!!$    real(double) :: tot
!!$    type(VECTOR) :: rVec
!!$    open(20, file=filename, status="old", form="formatted")
!!$    nStream = 0
!!$
!!$
!!$10 continue
!!$       read(20,*,end=20) r, theta, phi, v, rho
!!$       if (r < 1.00001d0) then
!!$          if (nStream /= 0) then
!!$             if (thisStream(nStream)%nSamples == 1) then
!!$                write(*,*) "Stream ",nstream, " only has one sample!!!!"
!!$                nStream = nStream - 1
!!$             endif
!!$          endif
!!$          nStream = nStream + 1
!!$          thisStream(nStream)%nSamples = 0
!!$          call allocateStream(thisStream(nStream), 200)
!!$       endif
!!$
!!$          theta = theta !* degtorad
!!$          phi = phi ! * degtorad
!!$          rho = rho * 1.d-27
!!$          r = r * rcore
!!$          rVec = r * VECTOR(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
!!$          
!!$          thisStream(nStream)%nSamples = thisStream(nStream)%nSamples + 1
!!$          
!!$          thisStream(nStream)%position(thisStream(nStream)%nSamples) = rVec
!!$          thisStream(nStream)%speed(thisStream(nStream)%nSamples) = v * 1.e5/cspeed
!!$          if (scaleFlowRho) then
!!$             thisStream(nStream)%rho(thisStream(nStream)%nSamples) = rho * flowRhoScale
!!$          else
!!$             thisStream(nStream)%rho(thisStream(nStream)%nSamples) = rho
!!$          endif
!!$          thisStream(nStream)%temperature(thisStream(nStream)%nSamples) = 7500.
!!$          thisStream(nStream)%streamRadius(thisStream(nStream)%nSamples)  = 0.01d0
!!$    goto 10
!!$20  continue
!!$    close(20)
!!$    write(*,*) "Read in ",nStream, " streams."
!!$    do i = 1, nStream
!!$       do j = 1, thisStream(i)%nSamples
!!$          if (j == thisStream(i)%nSamples) then
!!$             thisStream(i)%direction(j) = thisStream(i)%position(j) - thisStream(i)%position(j-1)
!!$             call normalize(thisStream(i)%direction(j))
!!$          
!!$          else
!!$             thisStream(i)%direction(j) = thisStream(i)%position(j+1) - thisStream(i)%position(j)
!!$             call normalize(thisStream(i)%direction(j))
!!$          endif
!!$             if (j == 1) then
!!$                thisStream(i)%distanceAlongStream(j) = 0.d0
!!$             else
!!$                thisStream(i)%distanceAlongStream(j) = thisStream(i)%distanceAlongStream(j-1) + &
!!$                     modulus(thisStream(i)%position(j) - thisStream(i)%position(j-1))
!!$             endif
!!$          ! direction is outward so need to use -speed...
!!$          thisStream(i)%velocity(j) = (-1.d0*thisStream(i)%speed(j)) * thisStream(i)%direction(j)
!!$       enddo
!!$    enddo
!!$
!!$! scale temperature according to distance along the stream
!!$    do i = 1, nStream
!!$       tot = MAXVAL(thisStream(i)%rho(1:thisStream(i)%nSamples))
!!$       do j = 1, thisStream(i)%nSamples
!!$          tot = 1.d0-(thisStream(i)%distanceAlongStream(j) / &
!!$               thisStream(i)%distanceAlongStream(thisStream(i)%nSamples))
!!$          thisStream(i)%temperature(j) = 4000.d0 + 4000.d0 * tot
!!$
!!$!          thisStream(i)%temperature(j) = 6000.d0
!!$   
!!$          ! Behaviour similar to Martin96 fig 1b
!!$          ! Arbitrary !!!!
!!$!          thisStream(i)%temperature(j) = min(2850.d0  + (tot/0.8)**3  * 3000.d0, 6000.d0)
!!$       enddo
!!$    enddo
!!$
!!$    ! Analytical description of the T along the stream 
!!$    ! reproduces behaviour calculated by Martin 1996 fig 1b
!!$!    do i = 1, nStream
!!$!       tot = MAXVAL(thisStream(i)%rho(1:thisStream(i)%nSamples))
!!$!       do j = 1, thisStream(i)%nSamples
!!$!          ! evolution is adiabatic up to 6000 K and then isothermal.
!!$!          ! adiabatic exponent is gamma = 5/3 for monoatomic gas (T.V^(gamma-1) = cst)
!!$!          thisStream(i)%temperature(j) = min(3000.d0 * (thisStream(i)%rho(j) / thisStream(i)%rho(1))**(2./3.), 6000.d0)
!!$!       enddo
!!$!    enddo
!!$    
!!$
!!$    do i = 1, nStream
!!$       call allocateStream(globalStream(i),thisStream(i)%nSamples)
!!$       globalStream(i) = thisStream(i)
!!$    enddo
!!$    globalnStream = nStream
!!$
!!$  end subroutine readStreams
!!$
!!$  subroutine readCurtain(thisCurtain,filename)
!!$
!!$    implicit none
!!$
!!$    type(curtaintype) :: thisCurtain
!!$    character(len=*) :: filename
!!$
!!$    integer :: i,j
!!$
!!$    real(double), parameter :: r0 = 4.e9, rho0 = 6.9e-12, v0 = 1.63e5
!!$
!!$    open(unit=1,file=filename,status="old")
!!$    read(1,*) thisCurtain%nr
!!$    read(1,*) thisCurtain%ntheta
!!$    read(1,*)
!!$    do i=1, thisCurtain%nr
!!$       do j=1, thisCurtain%ntheta
!!$          read(1,*) thisCurtain%position(i,j)%x, thisCurtain%position(i,j)%z, &
!!$               thisCurtain%density(i,j), thisCurtain%velocity(i,j)%x, &
!!$               thisCurtain%velocity(i,j)%z, thisCurtain%velocity(i,j)%y, &
!!$               thisCurtain%temperature(i,j)
!!$       enddo !j
!!$    enddo !i
!!$
!!$    ! SI system
!!$    thisCurtain%position(:,:)%x = thisCurtain%position(:,:)%x * r0
!!$    thisCurtain%position(:,:)%z = thisCurtain%position(:,:)%z * r0
!!$    thisCurtain%velocity(:,:)%x = thisCurtain%velocity(:,:)%x * v0
!!$    thisCurtain%velocity(:,:)%y = thisCurtain%velocity(:,:)%y * v0
!!$    thisCurtain%velocity(:,:)%z = thisCurtain%velocity(:,:)%z * v0
!!$    thisCurtain%temperature(:,:) = thisCurtain%temperature(:,:) / thisCurtain%density(:,:) 
!!$    thisCurtain%density(:,:) = thisCurtain%density(:,:) * rho0
!!$
!!$    close(unit=1)
!!$
!!$    return
!!$
!!$  end subroutine readCurtain
!!$
!!$
!!$  function  closestCorner(thisOctal, subcell, posVec) result (closeCorner)
!!$    type(VECTOR) :: closeCorner, posVec, cen, thisCorner
!!$    type(OCTAL), pointer :: thisOctal
!!$    real(double) :: minDist, r, dist
!!$    integer :: subcell
!!$    
!!$    minDist =  1.d30
!!$    r = thisOctal%subcellSize/2.d0
!!$    cen = subcellCentre(thisOctal, subcell)
!!$
!!$    thisCorner = cen + r * VECTOR(1.d0, 1.d0, 1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$       closeCorner = thisCorner
!!$    endif
!!$    thisCorner = cen + r * VECTOR(1.d0, 1.d0, -1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$       closeCorner = thisCorner
!!$    endif
!!$    thisCorner = cen + r * VECTOR(1.d0, -1.d0, -1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$       closeCorner = thisCorner
!!$    endif
!!$    thisCorner = cen + r * VECTOR(1.d0, -1.d0, 1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$       closeCorner = thisCorner
!!$    endif
!!$    thisCorner = cen + r * VECTOR(-1.d0, -1.d0,  -1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$       closeCorner = thisCorner
!!$    endif
!!$    thisCorner = cen + r * VECTOR(-1.d0, 1.d0, -1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$       closeCorner = thisCorner
!!$    endif
!!$    thisCorner = cen + r * VECTOR(-1.d0, -1.d0, 1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$       closeCorner = thisCorner
!!$    endif
!!$    thisCorner = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$       closeCorner = thisCorner
!!$    endif
!!$  end function closestCorner
!!$
!!$  function  dist_from_closestCorner(thisOctal, subcell, posVec) result (minDist)
!!$    type(VECTOR) :: posVec, cen, thisCorner
!!$    type(OCTAL), pointer :: thisOctal
!!$    real(double) :: minDist, r, dist
!!$    integer :: subcell
!!$    
!!$    minDist =  1.d30
!!$    r = thisOctal%subcellSize/2.d0
!!$    cen = subcellCentre(thisOctal, subcell)
!!$
!!$    thisCorner = cen + r * VECTOR(1.d0, 1.d0, 1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$    thisCorner = cen + r * VECTOR(1.d0, 1.d0, -1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$    thisCorner = cen + r * VECTOR(1.d0, -1.d0, -1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$    thisCorner = cen + r * VECTOR(1.d0, -1.d0, 1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$    thisCorner = cen + r * VECTOR(-1.d0, -1.d0,  -1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$    thisCorner = cen + r * VECTOR(-1.d0, 1.d0, -1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$    thisCorner = cen + r * VECTOR(-1.d0, -1.d0, 1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$    thisCorner = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
!!$    dist  =  modulus(posVec - thisCorner)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$  end function dist_from_closestCorner
!!$
!!$  function  dist_from_closestEdge(thisOctal, subcell, posVec) result (minDist)
!!$    type(VECTOR) :: posVec, cen, corner1, corner2
!!$    type(OCTAL), pointer :: thisOctal
!!$    real(double) :: minDist, r, dist
!!$    integer :: subcell
!!$    
!!$    minDist =  1.d30
!!$    r = thisOctal%subcellSize/2.d0
!!$    cen = subcellCentre(thisOctal, subcell)
!!$
!!$!1
!!$    corner1 = cen + r * VECTOR(-1.d0,-1.d0, 1.d0)
!!$    corner2 = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!2
!!$    corner1 = cen + r * VECTOR(-1.d0,-1.d0, 1.d0)
!!$    corner2 = cen + r * VECTOR( 1.d0,-1.d0, 1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!3
!!$    corner1 = cen + r * VECTOR(-1.d0,-1.d0, 1.d0)
!!$    corner2 = cen + r * VECTOR(-1.d0,-1.d0,-1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!4
!!$    corner1 = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
!!$    corner2 = cen + r * VECTOR( 1.d0, 1.d0, 1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!5
!!$    corner1 = cen + r * VECTOR(-1.d0, 1.d0, 1.d0)
!!$    corner2 = cen + r * VECTOR(-1.d0, 1.d0,-1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!6
!!$    corner1 = cen + r * VECTOR( 1.d0,-1.d0, 1.d0)
!!$    corner2 = cen + r * VECTOR( 1.d0, 1.d0, 1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!7
!!$    corner1 = cen + r * VECTOR( 1.d0,-1.d0, 1.d0)
!!$    corner2 = cen + r * VECTOR( 1.d0,-1.d0,-1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!8
!!$    corner1 = cen + r * VECTOR( 1.d0, 1.d0, 1.d0)
!!$    corner2 = cen + r * VECTOR( 1.d0, 1.d0,-1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!9
!!$    corner1 = cen + r * VECTOR(-1.d0, 1.d0,-1.d0)
!!$    corner2 = cen + r * VECTOR( 1.d0, 1.d0,-1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!10
!!$    corner1 = cen + r * VECTOR(-1.d0, 1.d0,-1.d0)
!!$    corner2 = cen + r * VECTOR(-1.d0,-1.d0,-1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!11
!!$    corner1 = cen + r * VECTOR(-1.d0,-1.d0,-1.d0)
!!$    corner2 = cen + r * VECTOR( 1.d0,-1.d0,-1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$!12
!!$    corner1 = cen + r * VECTOR( 1.d0,-1.d0,-1.d0)
!!$    corner2 = cen + r * VECTOR( 1.d0, 1.d0,-1.d0)
!!$    dist  =  distancePointLineSegment(corner1, corner2, posVec)
!!$    if (dist < minDist) then
!!$       minDist = dist
!!$    endif
!!$
!!$
!!$  end function dist_from_closestEdge
!!$

!!$  subroutine getMagStreamValues2(point,  rho, temperature, &
!!$       velocity, inflow)
!!$    type(STREAMTYPE) :: thisStream
!!$    integer :: iSample
!!$    type(VECTOR) :: point
!!$    real(double) :: t, thisR
!!$    real(double),optional :: rho
!!$    real,optional :: temperature
!!$    type(VECTOR), optional :: velocity
!!$    logical, optional :: inFlow
!!$    logical :: outsideStream
!!$    outsideStream = .false.
!!$    iSample = 0
!!$    t = 0.d0
!!$
!!$    call findNearestSample(thisStream, point, iSample, t, outsideStream)
!!$
!!$
!!$    thisR =  thisStream%streamradius(isample) + &
!!$         t * (thisStream%streamradius(iSample+1)-thisStream%streamradius(iSample))
!!$
!!$!    outsideStream = .false.
!!$!    if (modulus(thisStream%position(isample)-point) > thisR) then
!!$!       outSideStream = .true.
!!$!    endif
!!$!
!!$!    if (iSample == thisStream%nSamples) then
!!$!       if (((point-thisStream%position(iSample)).dot.thisStream%position(iSample)) < 0.d0) then
!!$!          outsideStream = .true.
!!$!       endif
!!$!    endif
!!$
!!$    if (PRESENT(rho)) then
!!$       if (outSideStream) then
!!$          rho = 1.d-25
!!$       else
!!$          rho = thisStream%rho(isample) + t * (thisStream%rho(iSample+1)-thisStream%rho(iSample))
!!$       endif
!!$    endif
!!$
!!$    if (PRESENT(temperature)) then
!!$       if (outSideStream) then
!!$          temperature = 6000.
!!$       else
!!$          temperature = thisStream%temperature(isample) + &
!!$               t * (thisStream%temperature(iSample+1)-thisStream%temperature(iSample))
!!$       endif
!!$    endif
!!$
!!$    if (PRESENT(velocity)) then
!!$       if (outsideStream) then
!!$          velocity = VECTOR(1.e-30, 1.e-30, 1.e-30)
!!$       else
!!$          velocity = thisStream%velocity(isample) + &
!!$               t * (thisStream%velocity(iSample+1)-thisStream%velocity(iSample))
!!$       endif
!!$    endif
!!$
!!$    if (present(inflow)) then
!!$       inflow = .not.outsideStream
!!$    endif
!!$  end subroutine getMagStreamValues2
!!$

!!$
!!$  function  getVel(grid, thisOctal, subcell, posVec, direction) result (outVec)
!!$    type(GRIDTYPE) :: grid
!!$    type(OCTAL), pointer :: thisOctal, neighbourOctal
!!$    integer subcell
!!$    real(double) :: distToNextCell, r, x1
!!$    type(VECTOR) :: direction, rVec, sVec, posVec, thisDirection
!!$    type(VECTOR) :: v1,  outVec
!!$    integer :: neighbourSubcell
!!$
!!$    rVec = subcellCentre(thisOctal, subcell)
!!$    r = ((posVec - rVec).dot.direction)
!!$    if (r  > 0.d0 ) then
!!$
!!$       thisDirection = direction
!!$       rVec = subcellCentre(thisOctal, subcell)
!!$       call distanceToCellBoundary(grid, rVec, thisdirection, DisttoNextCell)
!!$       v1 = thisOctal%velocity(subcell)
!!$       x1 = distToNextCell
!!$       distToNextCell = distToNextcell + 0.01d0*grid%halfSmallestSubcell
!!$       sVec = rVec + distToNextCell * thisdirection
!!$       if (inOctal(grid%octreeRoot, sVec)) then
!!$          neighbourOctal => thisOctal
!!$          CALL findSubcellLocal(sVec,neighbourOctal,neighboursubcell)
!!$          if (neighbourOctal%inFlow(neighbourSubcell)) then
!!$             v1 = neighbourOctal%velocity(neighbourSubcell)
!!$             x1 = disttoNextcell + neighbourOctal%subcellSize/2.d0
!!$          endif
!!$       endif
!!$       outVec = thisOctal%velocity(subcell) + (r/x1) * (v1 - thisOctal%velocity(subcell))
!!$    else
!!$       thisDirection = (-1.d0)*direction
!!$       rVec = subcellCentre(thisOctal, subcell)
!!$       call distanceToCellBoundary(grid, rVec, thisdirection, DisttoNextCell)
!!$       v1 = thisOctal%velocity(subcell)
!!$       x1 = distToNextCell
!!$       distToNextCell = distToNextcell + 0.01d0*grid%halfSmallestSubcell
!!$       sVec = rVec + distToNextCell * thisdirection
!!$       if (inOctal(grid%octreeRoot, sVec)) then
!!$          neighbourOctal => thisOctal
!!$          CALL findSubcellLocal(sVec,neighbourOctal,neighboursubcell)
!!$          if (neighbourOctal%inFlow(neighbourSubcell)) then
!!$             v1 = neighbourOctal%velocity(neighbourSubcell)
!!$             x1 = disttoNextcell + neighbourOctal%subcellSize/2.d0
!!$          endif
!!$       endif
!!$       outVec = thisOctal%velocity(subcell) + (abs(r)/x1) * (v1 - thisOctal%velocity(subcell))
!!$    endif
!!$
!!$  end function getVel

!!$  RECURSIVE SUBROUTINE cleanupAMRgrid(thisOctal)
!!$      
!!$        TYPE(OCTAL), POINTER  :: thisOctal 
!!$        TYPE(OCTAL), POINTER  :: child
!!$        INTEGER :: i
!!$        
!!$        if (thisOctal%nChildren == thisOctal%maxChildren) then
!!$           call deallocateOctalDynamicAttributes(thisOctal)
!!$        endif
!!$
!!$        IF ( thisOctal%nChildren > 0 ) THEN
!!$           ! call this subroutine recursively on each of its children
!!$           DO i = 1, thisOctal%nChildren, 1
!!$              child => thisOctal%child(i)
!!$              CALL cleanupAMRgrid(child)
!!$           END DO
!!$        END IF
!!$        
!!$  END SUBROUTINE cleanupAMRgrid

!!$   
!!$  subroutine testToBoundary(grid, rVec, uHat, totDist)
!!$    type(GRIDTYPE) :: grid
!!$    type(VECTOR) :: rVec, uHat, currentPosition
!!$    type(OCTAL), pointer :: thisOctal, sOctal
!!$    integer :: subcell
!!$    real(double) :: distToNextCell
!!$    real(double), intent(out) :: totDist
!!$    real(double) :: fudgeFac = 0.01d0
!!$
!!$    totDist = 0.d0
!!$    currentPosition = rVec
!!$    CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
!!$
!!$
!!$    do while (inOctal(grid%octreeRoot, currentPosition))
!!$
!!$       call findSubcellLocal(currentPosition,thisOctal,subcell)
!!$
!!$       sOctal => thisOctal
!!$       call distanceToCellBoundary(grid, currentPosition, uHat, DisttoNextCell, sOctal)
!!$  
!!$       currentPosition = currentPosition + (distToNextCell+fudgeFac*grid%halfSmallestSubcell)*uHat
!!$       totDist = totDist + distToNextCell
!!$
!!$    end do
!!$  end subroutine testToBoundary
!!$
!!$  subroutine testToBoundary2(grid, rVec, uHat, totDist)
!!$    type(GRIDTYPE) :: grid
!!$    type(VECTOR) :: rVec, uHat, currentPosition
!!$    type(OCTAL), pointer :: thisOctal, sOctal
!!$    integer :: subcell, sSubcell
!!$    real(double) :: distToNextCell, totDist
!!$!    real(double) :: fudgeFac = 0.000001d0
!!$
!!$    totDist = 0.d0
!!$    currentPosition = rVec
!!$    CALL findSubcellTD(currentPosition,grid%octreeRoot,thisOctal,subcell)
!!$
!!$    sOctal => thisOctal
!!$    sSubcell = subcell
!!$
!!$    do while (sSubcell > 0)
!!$
!!$       call distanceToCellBoundary(grid, currentPosition, uHat, DisttoNextCell, thisOctal, subcell)
!!$
!!$       currentPosition = currentPosition + distToNextCell*uHat
!!$       call getNeighbourFromPointOnFace(currentPosition, uHat, thisOctal, subcell, sOctal, sSubcell)
!!$       thisOctal => sOctal
!!$       subcell = sSubcell
!!$
!!$!       currentPosition = currentPosition + (fudgeFac*grid%halfSmallestSubcell)*uHat
!!$
!!$
!!$
!!$       totDist = totDist + distToNextCell
!!$
!!$    end do
!!$  end subroutine testToBoundary2
!!$

!  subroutine assignDensitiesMahdavi(grid)
!    use input_variables, only : ttauriRstar, dipoleOffset, isothermTemp, vturb, mdotParameter1
!    use magnetic_mod, only : inflowMahdavi, velocityMahdavi
!    type(GRIDTYPE) :: grid
!    integer :: i, j
!    type(VECTOR) :: rVec, rVecDash, thisRvec, rVec1, cellcentre
!    real(double) :: mdot
!    integer :: nLines, iline
!    type(OCTAL), pointer :: thisOctal
!    integer :: subcell, nr
!    real(double) :: accretingArea, vstar, astar, beta
!    real(double) :: thetaDash, phiDash, rDash, cosThetaDash, mDotFraction
!    real(double) :: thisrMax, sin2theta0dash, thisr,thisrho, thisPhi, thisTheta,ds, thisv
!
!
!    mdot = mDotparameter1*mSol/(365.25d0*24.d0*3600.d0)
!    beta = dipoleOffset
!    nLines = 1000000
!    j = 0
!    do i = 1, nLines
!       rVec = (ttauriRstar + 1.d10 *grid%halfSmallestSubcell)*randomUnitVector() 
!       if (inFlowMahdavi(rVec)) j = j + 1
!    enddo
!    accretingArea = fourPi * ttauriRstar**2 * dble(j)/dble(nLines)
!    if (Writeoutput) then
!       write(*,*) "Fractional area of accretion (%): ",100.d0 *  dble(j)/dble(nLines)
!       write(*,*) "Using mdot of ",(mdot/msol)*(365.25d0*24d0*3600d0)
!    endif
!    nLines = 100000
!    aStar = accretingArea / dble(nLines)
!    mDotFraction = mDot / dble(nLines)
!    do iLine = 1, nLines
!       rVec = (ttauriRstar + 1.d10 *grid%halfSmallestSubcell)*randomUnitVector() 
!       do while(.not.inFlowMahdavi(rVec))
!          rVec = (ttauriRstar + 1.d10 *grid%halfSmallestSubcell)*randomUnitVector() 
!       enddo
!       call findSubcellTD(rVec/1.d10, grid%octreeRoot, thisOctal, subcell)
!       thisOctal%velocity(subcell) = velocityMahdavi(rVec/1.d10, grid)
!       vStar = modulus(thisOctal%velocity(subcell)*cspeed)
!       rVecDash = rotateY(rVec, -beta)
!       rDash = modulus(rVecDash)
!       cosThetaDash = rVecDash%z / rDash
!       thetaDash = acos(rVecDash%z/rDash)
!       phiDash = atan2(rVecDash%y, rVecDash%x)
!       sin2theta0dash = (1.d0 + tan(beta)**2 * cos(phiDash)**2)**(-1.d0)
!       thisrMax = rDash / sin(thetaDash)**2
!
!       nr = 10000
!       thisOctal => grid%octreeRoot
!       rVec1 = rVec
!       do i = 1, nr
!          thisr = rDash + (thisRmax-rDash) * dble(i-1)/dble(nr-1)
!          thisphi = phiDash
!          thisTheta = asin(sqrt(thisr/thisRmax))
!          if (cosThetaDash  < 0.d0) thisTheta = pi - thisTheta
!          thisRVec = VECTOR(thisR*cos(thisPhi)*sin(thisTheta),thisR*sin(thisPhi)*sin(thisTheta), &
!               thisR * cos(thisTheta))
!          thisRvec = rotateY(thisRvec, beta)
!          if (inflowMahdavi(thisRVec)) then
!             if (inOctal(grid%octreeRoot, thisRvec/1.d10)) then
!                
!                ds = modulus(thisRvec - rVec1)
!                rVec1 = thisRvec
!                call findSubcellLocal(thisRvec/1.d10, thisOctal, subcell)
!                cellCentre = subcellCentre(thisOctal, subcell)
!!                if (inflowMahdavi(1.d10*cellCentre)) then
!                   thisOctal%velocity(subcell) = velocityMahdavi(thisrVec/1.d10, grid)
!                   thisV = modulus(thisOctal%velocity(subcell))*cSpeed
!                   if (thisV /= 0.d0) then
!                      thisRho =  mdotFraction /(aStar * thisV)  * (ttauriRstar/thisR)**3 
!                      thisOctal%rho(Subcell) = thisRho
!                      thisOctal%inflow(subcell) = .true.
!                      thisOctal%temperature(subcell) = isothermTemp
!                      if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = vturb
!                   endif
!
!                endif
!!             endif
!          endif
!       enddo
!    end do
!!    call convertToDensity(grid%octreeRoot)
!  end subroutine assignDensitiesMahdavi
!!$
!!$  subroutine createTTauriSurfaceMahdavi(surface, grid, lineFreq, &
!!$                                 coreContFlux,fAccretion)
!!$    use magnetic_mod, only : inflowMahdavi
!!$    use surface_mod, only : surfaceType, createProbs, sumSurface
!!$    type(SURFACETYPE),intent(inout) :: surface
!!$    type(gridType), intent(in) :: grid
!!$    type(OCTAL), pointer :: thisOctal
!!$    integer :: subcell
!!$    real(double) :: v, rho,mdot,power
!!$    real(double), intent(in) :: coreContFlux
!!$    real, intent(in) :: lineFreq
!!$    real, intent(out) :: fAccretion ! erg s^-1 Hz^-1
!!$    integer :: iElement
!!$    type(VECTOR) :: aboveSurface
!!$    real(double) :: Taccretion
!!$    
!!$    call writeInfo("Creating T Tauri stellar surface",TRIVIAL)
!!$    
!!$    do iElement = 1, SIZE(surface%element)
!!$      aboveSurface = surface%element(iElement)%position - surface%centre
!!$      aboveSurface = aboveSurface * 1.01_oc 
!!$
!!$
!!$      surface%element(iElement)%hot = .false.
!!$
!!$      if (inFlowMahdavi(aboveSurface*1.d10)) then
!!$         call findSubcellTD(aboveSurface, grid%octreeRoot, thisOctal, subcell)
!!$         v = modulus(thisOctal%velocity(subcell))*cspeed
!!$         rho = thisOctal%rho(subcell)
!!$         mdot = v * rho * surface%element(iElement)%area*1.d20
!!$         power = 0.5d0*v*mdot**2
!!$         taccretion = (power/(stefanBoltz*surface%element(iElement)%area*1.d20))**0.25d0
!!$
!!$        surface%element(iElement)%hot = .true.
!!$        allocate(surface%element(iElement)%hotFlux(surface%nNuHotFlux))
!!$        surface%element(iElement)%hotFlux(:) = &
!!$           pi*blackbody(REAL(tAccretion), 1.e8*real(cSpeed)/surface%nuArray(:)) 
!!$        surface%element(iElement)%temperature = Taccretion
!!$     endif
!!$    end do
!!$    call createProbs(surface,lineFreq,coreContFlux,fAccretion)
!!$    call sumSurface(surface)
!!$    
!!$  end subroutine createTTauriSurfaceMahdavi


