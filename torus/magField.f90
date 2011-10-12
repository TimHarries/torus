MODULE magField

USE vector_mod

IMPLICIT NONE

TYPE gridsample
  TYPE(vector) :: position ! same units as AMR grid  
  REAL :: velocity ! velocity magnitude (cSpeed)
  TYPE(vector) :: flowVector 
    ! local direction of flow.
    ! should be normalised, and (usually) directed inwards towards star.
  TYPE(vector) :: prevFlowVector ! flowVectors for the next and
  TYPE(vector) :: nextFlowVector !   previous samples
  REAL :: nextVelocity ! velocity magnitude of next sample(cSpeed)
  REAL :: prevVelocity ! velocity magnitude of previous sample (cSpeed)
  REAL :: rho ! g cm-3
  REAL(oct) :: radius ! radius of (circular) stream (1.e10cm)
  REAL :: temperature ! K
  REAL(oct) :: prevDistance ! used to interpolate between adjacent samples
                            !   in a stream. this is distance of this sample
                            !   from the previous one. note that the "previous"
                            !   sample will actually be in the forward direction
                            !   (in the direction defined by flowvector)
  REAL(oct) :: nextDistance ! distance to next sample (1.e10 cm)
  TYPE(vector) :: prevDirection ! directions to previous/next samples in flow 
  TYPE(vector) :: nextDirection
  REAL(oct) :: distanceUpperLimit 
    ! maximum possible distance that a point could be from this gridSample, 
    !   and still be in the flow defined by this gridSample. used to quickly
    !   reject points in the density finding routine
  LOGICAL :: endOfStream ! .true. if sample is last one in a stream
END TYPE gridSample

TYPE hotSpotVariable
  TYPE(spVector) :: position ! position, r should be 1.0  
  REAL :: speed ! (cm s-1)
  REAL :: rho ! g cm-3
  REAL(oct) :: radius ! radius of (circular) stream (radians)
  REAL :: temperature ! K

  ! variables below are for an alternative implementation, for testing
  TYPE(vector) :: XYZposition ! centre of hot spot (usual grid coords + units) 
  TYPE(vector) :: subSurfacePosn ! a position directly below the centre of the
      !   hot spot; if a surface element has a (cartesian) distance that is less than 
      !   "radius_1e10" from this point, it is within the spot
      !   (usual grid unit; coordinates are relative to stellar centre) 
  REAL(oct) :: radius_1e10 ! radius of (circular) stream (1.e10cm)
  
END TYPE hotSpotVariable

TYPE(gridSample), ALLOCATABLE, TARGET, SAVE :: magFieldGrid(:)
TYPE(hotSpotVariable), ALLOCATABLE, SAVE :: magFieldHotspots(:)
INTEGER, ALLOCATABLE :: nSamplesInStream(:)
REAL, SAVE :: maxSizeMagFieldGrid = -HUGE(1.0) ! maximum radius of accretion (1.e10cm)

TYPE gridSampleResults
  ! this type stores the results from finding whether a particular
  !   'gridSample' is close to a point in the AMR grid
  INTEGER :: iSample
  REAL(oct) :: sampleDistance 
    ! distance between point and gridSample point (1.e10cm)
END TYPE  

TYPE(gridSampleResults), ALLOCATABLE :: sampleResults(:)

  ! the following two arrays store information about how far inwards
  !   the accretion disk extends at different azimuth angles
REAL(oct), DIMENSION(:), ALLOCATABLE :: innerDiskData_Phi ! radians
REAL(oct), DIMENSION(:), ALLOCATABLE :: innerDiskData_Radii ! 1.e10 cm

REAL(oct), PARAMETER :: minusOneOctal = -1.0_oc

CONTAINS
 
  SUBROUTINE loadMagField(fileName,starPosn,Rstar)
    ! loads in data defining magnetically controlled accretion
    !   flows around a star, from the models of Gregory,
    !   Jardine et al.

    USE inputs_mod, ONLY: isothermStream, isothermTemp, &
                               flowRhoScale, scaleFlowRho,   &
                               magStreamFileDegrees, limitSpotTemp, &
                               maxSpotTemp
    USE constants_mod
   
    CHARACTER(LEN=*), INTENT(IN) :: fileName ! input file
    TYPE(vector), INTENT(IN) :: starPosn 
      ! position of star in grid
    REAL, INTENT(IN) :: Rstar ! (1.e10 cm)
      ! radius of star, used to de-normalize units
    
    INTEGER :: nSamples ! total num. of samples in file
    INTEGER :: nStreams ! total num. of accretion streams
    REAL :: rValue     ! r, from file  
    REAL :: thetaValue ! theta, from file  
    REAL :: phiValue   ! phi, from file  
    REAL :: velValue   ! velocity, from file  
    REAL :: rhoValue   ! rho, from file  
    REAL(oct) :: areaValue  ! cross-sectional area, from file  
    INTEGER :: fileRead ! loop counter
    TYPE(spVector) :: spPosition ! position, spherical polar coords.
    TYPE(vector) :: prevVector
    TYPE(gridSample), POINTER :: thisSample
    TYPE(gridSample), POINTER :: prevSample
    INTEGER :: iStream ! index of current stream
    INTEGER :: iSample ! index of current sample 
    INTEGER :: iSampleInStream ! index of sample in current stream
    INTEGER :: returnVal
    REAL(double) :: localAccretionRate
    REAL(double) :: totalAccretionRate
    REAL :: distance
    TYPE(spVector) :: spSubSurfacePosn ! sub-surface position (sph. pol.) 

    IF (scaleFlowRho) THEN
      PRINT *, "Will rescale densities by factor of: ",flowRhoScale
    ELSE
      PRINT *, "Will not rescale densities"
    END IF
    
    ! file format is:
    
    ! Col 1: r-coordinate in units of R_*
    ! Col 2: theta-coordinate in radians 
    !   (0 at North pole, pi at south pole)
    ! *OR* Col 2: theta-coordinate in degrees 
    !        (0 at North pole, 180 at south pole)
    ! Col 3: phi-coordinate in radians (0 - 2pi)
    ! *OR* Col 3: phi-coordinate in degrees (0 - 360)
    ! Col 4: velocity in km/s
    ! Col 5: mass density in kg/m^3
    ! Col 6: cross sectional area (m^2)
    totalAccretionRate = 0.0_db
    
    fileReadLoop: &
    DO fileRead = 1, 2
      ! we will read all the data from the file twice. the first
      !   time we find out how big we need to make the arrays. the second
      !   time we store the data properly, and compute
      !   the local flow vector.

      SELECT CASE ( fileRead )
      CASE (1)
        PRINT *, "Reading magnetic field data: phase 1, getting sizes"
      CASE (2)
        PRINT *, "Reading magnetic field data: phase 2, storing data"
      END SELECT
      
      OPEN(UNIT=27,FILE=fileName,STATUS="old",FORM="formatted")
      nStreams = 0
      iStream  = 0
      IF (fileRead == 1) nSamples = 0
      iSample  = 0
      iSampleInStream = 1
    
      readLineLoop: &
      DO
        
        READ(UNIT=27,FMT=*,IOSTAT=returnVal) &
               rValue, thetaValue, phiValue, velValue, rhoValue, areaValue
          IF (magStreamFileDegrees) THEN
            ! need to convert to radians
            thetaValue = thetaValue * real(DegToRad)
            phiValue = phiValue * real(DegToRad)
          END IF
        IF (returnVal /= 0) THEN
          EXIT 
        ELSE
          iSample = iSample + 1
        END IF

        ! we assume we have reached a new accretion stream if we
        !   find an r = 1 value, or r < (previousR)
        
        IF ( ABS(rValue-1.0) < 5*EPSILON(1.) ) THEN
          ! new stream
          iStream = iStream + 1
          iSampleInStream = 1 
          nStreams = MAX(iStream, nStreams)

          IF (fileRead == 2) THEN
            ! safety check of previous stream
            IF ( iStream > 1 ) THEN
              IF ( nSamplesInStream(iStream-1) < 2 ) THEN
                PRINT *, "Problem in loadMagField: this stream only has one sample:",(iStream-1)
              END IF
            END IF
          END IF
        ELSE
          iSampleInStream = iSampleInStream + 1
        END IF  

        
        IF (fileRead == 2) nSamplesInStream(iStream) = nSamplesInStream(iStream) + 1
        IF (fileRead == 1) nSamples = MAX(nSamples, iSample)
        
        IF ( fileRead == 2 ) THEN
          ! we want to store the data in our arrays
          
          spPosition = spVector(rValue, thetaValue, phiValue)
          ! transform to grid coordinates
          spPosition%r = spPosition%r * Rstar
          magFieldGrid(iSample)%position = spPosition ! (includes conversion to
                                              ! cartesian coordinates)
          magFieldGrid(iSample)%position =                           &
            magFieldGrid(iSample)%position + &
              vector(starPosn%x,starPosn%y,starPosn%z)
                                              
          magFieldGrid(iSample)%velocity = real((velValue * 1.e5) / cSpeed )
!          magFieldGrid(iSample)%velocity = (velValue * 1.e-5) 
          IF (scaleFlowRho) THEN
            rhoValue = rhoValue * flowRhoScale
          END IF
          magFieldGrid(iSample)%rho = rhoValue * 1.e-3 ! convert to g cm-3
          areaValue = areaValue * 1.e4 ! m^2 to cm^2 !!!!!! changed from neils 1.e5
          areaValue = SQRT(areaValue / pi) ! area to radius
          magFieldGrid(iSample)%radius = areaValue / 1.e10_oct ! to 1.e10 cm

          maxSizeMagFieldGrid = MAX(maxSizeMagFieldGrid, real(spPosition%r) )

          ! store inner disk radius, for this stream.
          ! we don't care if this is not the last sample in this stream,
          !   because we will overwrite it later with the correct value
          innerDiskData_Phi(iStream) = phiValue
          innerDiskData_Radii(iStream) = rValue * rStar
          
          IF ( .NOT. isothermStream ) THEN
            PRINT *, "Error: accretion streams must (for now) be isothermal"
            STOP
          ELSE
            magFieldGrid(iSample)%temperature = isothermTemp
          END IF
          
          ! if we are in the middle of an accretion stream, we find the
          !   local path vector of the stream by averaging the line 
          !   segments to the previous and the next samples. if we are 
          !   at an end of the stream, we can only get one line segment.
          ! (remember that we are working outwards from the star, so the 
          !  direction of the vectors is reversed)

          IF ( iSampleInStream == 1 ) THEN
              
            thisSample => magFieldGrid(iSample)
            thisSample%flowVector = zeroVector

            ! store data about a surface hotspot
            magFieldHotspots(iStream)%position = spPosition

            magFieldHotspots(iStream)%XYZposition = thisSample%position
            magFieldHotspots(iStream)%XYZposition = &
              magFieldHotspots(iStream)%XYZposition + starPosn
            magFieldHotspots(iStream)%radius_1e10 = thisSample%radius

            magFieldHotspots(iStream)%position%r = 1.0
            magFieldHotspots(iStream)%speed = real(thisSample%velocity * cSpeed)
            ! approximation, radius will not be quite correct 
            magFieldHotspots(iStream)%radius = &
              thisSample%radius / Rstar ! 1.e10 cm to radians
            magFieldHotspots(iStream)%rho = thisSample%rho

            ! for temperature, assume uniform circular cross-section
            !   and all kinetic energy radiated as a blackbody
            magFieldHotspots(iStream)%temperature =          & 
              real(((0.5 * magFieldHotspots(iStream)%rho *        &
                (magFieldHotSpots(iStream)%speed)**3 )         &
                / stefanBoltz )**0.25 )
                
            IF (limitSpotTemp) CALL rescaleHotSpot(magFieldHotspots(iStream), maxSpotTemp)

            ! work out how far beneath the surface "subSurfacePosn" should lie
            distance = real(COS(magFieldHotspots(iStream)%radius)) ! (Rstar)
            spSubSurfacePosn = magFieldHotspots(iStream)%position
            spSubSurfacePosn%r = distance * Rstar
            magFieldHotspots(iStream)%subSurfacePosn = spSubSurfacePosn ! convert to XYZ
              ! note that we haven't added on "starPosn" 
            

            localAccretionRate = magFieldHotspots(iStream)%speed * &
            pi * (REAL(thisSample%radius*1.e10,KIND=double))**2 * &
                                   magFieldHotspots(iStream)%rho
            totalAccretionRate = localAccretionRate + totalAccretionRate                       
            
          ELSE IF ( iSampleInStream == 2 ) THEN
            ! store a line segment for the first sample in this stream
            prevSample => magFieldGrid(iSample-1)
            thisSample => magFieldGrid(iSample)
            
            prevSample%flowVector = prevSample%position - thisSample%position
            prevSample%nextDistance = modulus(prevSample%flowVector)
            prevSample%prevDistance = prevSample%nextDistance
            CALL normalize(prevSample%flowVector)  
            prevSample%prevDirection = prevSample%flowVector
            prevSample%nextDirection = minusOneOctal * prevSample%prevDirection
            prevSample%prevFlowVector = prevSample%prevFlowVector
            
            thisSample%prevDirection = prevSample%prevDirection
            thisSample%prevDistance = prevSample%nextDistance
            thisSample%prevFlowVector = prevSample%flowVector
            thisSample%prevVelocity = prevSample%velocity
            prevSample%nextFlowVector = thisSample%prevFlowVector
            
            prevSample%prevVelocity = prevSample%velocity
            prevSample%nextVelocity = thisSample%velocity
            
            ! in case this stream only has 2 samples, we should fill in this
            !   samples's values
            thisSample%flowVector = thisSample%prevDirection 
            CALL normalize(thisSample%flowVector)  
            thisSample%nextDistance = thisSample%prevDistance
            thisSample%nextFlowVector = thisSample%flowVector

          ELSE IF ( iSampleInStream > 2 ) THEN
            prevSample => magFieldGrid(iSample-1)
            thisSample => magFieldGrid(iSample)
              
            ! the general case:
            prevVector = prevSample%position - thisSample%position
            thisSample%prevDistance = modulus(prevVector)
            prevSample%nextDistance = thisSample%prevDistance
            CALL normalize(prevVector)
            thisSample%prevDirection = prevVector
            prevSample%nextDirection = minusOneOctal * prevVector
            
            prevSample%flowVector = & ! get weighted average of vectors
              ( (1.0_oc/prevSample%prevDistance) * prevSample%prevDirection) + &
              ( (1.0_oc/prevSample%nextDistance) * prevSample%nextDirection)
            CALL normalize(prevSample%flowVector)

            prevSample%nextVelocity = thisSample%velocity
            thisSample%prevVelocity = prevSample%velocity
            
            ! in case we are at the end of the stream, fill all this sample's
            !   variables:
            thisSample%flowVector = thisSample%prevDirection 
            CALL normalize(thisSample%flowVector)  
            thisSample%nextDistance = thisSample%prevDistance
            thisSample%nextFlowVector = thisSample%flowVector
            thisSample%nextVelocity = thisSample%velocity

          END IF ! ( iSampleInStream == 1 )
            
        END IF ! ( fileRead == 2 )
      
      END DO readLineLoop
      
      CLOSE(UNIT=27)
      
      IF ( fileRead == 1 ) THEN
        PRINT *, "   ",nSamples," samples found"
        PRINT *, "   ",nStreams," streams found"
        ALLOCATE(nSamplesInStream(nStreams))
        nSamplesInStream = 0
        ALLOCATE(magFieldHotspots(nStreams))
        ALLOCATE(magFieldGrid(nSamples))
        ALLOCATE(innerDiskData_Phi(nStreams))
        ALLOCATE(innerDiskData_Radii(nStreams))
        iStream = 0
      END IF
     
      SELECT CASE ( fileRead )
      CASE (1)
        PRINT *, "                             phase 1 complete"
      CASE (2)
        PRINT *, "                             phase 2 complete"
      END SELECT
      
    END DO fileReadLoop    
    
    ! store distances between samples, for interpolating quantities later
    PRINT *, "Reading magnetic field data: phase 3, storing distances"
    
    iStream = 1 
    iSampleInStream = 1
   
    DO iSample = 1, nSamples
      magFieldGrid(iSample)%endOfStream = .FALSE.
      
      ! store distance to previous sample
      IF ( iSampleInStream == 1 ) THEN
        magFieldGrid(iSample)%prevDistance = 0.0_oc
      ELSE 
        magFieldGrid(iSample)%prevDistance =  &
          modulus( magFieldGrid(iSample)%position - magFieldGrid(iSample-1)%position )
      END IF
        
      ! store distance to next sample

      IF ( iSampleInStream == nSamplesInStream(iStream) ) THEN
        magFieldGrid(iSample)%nextDistance = 0.0_oc
      ELSE
        magFieldGrid(iSample)%nextDistance = &
          modulus( magFieldGrid(iSample+1)%position - magFieldGrid(iSample)%position )
      END IF
      
      IF ( iSampleInStream < nSamplesInStream(iStream) ) THEN
        iSampleInStream = iSampleInStream + 1
      ELSE
        magFieldGrid(iSample)%endOfStream = .TRUE.
        iStream = iStream + 1 
        iSampleInStream = 1
      END IF

      ! to speed up searching later, we find an upper limit for the distance
      !   that a point can be from a gridSample, and still has the possibility
      !   of being in the flow defined by that sample
      magFieldGrid(iSample)%distanceUpperLimit =                                    &
        MAX(magFieldGrid(iSample)%nextDistance, magFieldGrid(iSample)%prevDistance) &
        + magFieldGrid(iSample)%radius  
            
    END DO
    PRINT *, "                             phase 3 complete"
    PRINT *, " maximum size of magnetic field = ",maxSizeMagFieldGrid/rStar," stellar radii"
    PRINT *, " total accretion rate = ", &
                  REAL(totalAccretionRate/mSol * (1.0_db / secsToYears))
    
    DEALLOCATE(nSamplesInStream)
    
    !   for the nearest sample to a point in the simulation
    ALLOCATE( sampleResults(SIZE(magFieldGrid)) )

    CALL sortInnerDiskData(innerDiskData_Phi,innerDiskData_Radii)
  
  END SUBROUTINE loadMagField

  SUBROUTINE sortInnerDiskData(innerDiskData_Phi,innerDiskData_Radii)
    ! this moves the data in the innerDiskData arrays so that
    !   the smallest phi angle is at index 1
    USE utils_mod
    REAL(oct), DIMENSION(:), INTENT(INOUT) :: innerDiskData_Phi
    REAL(oct), DIMENSION(:), INTENT(INOUT) :: innerDiskData_Radii
      ! radius of star, used for debugging
      
    REAL(oct), DIMENSION(SIZE(innerDiskData_Phi)) :: Phi_temp
    REAL(oct), DIMENSION(SIZE(innerDiskData_Phi)) :: Radii_temp
    INTEGER, DIMENSION(SIZE(innerDiskData_Phi)) :: phiIndex
    INTEGER :: i

    if (SIZE(innerDiskdata_phi,1) > 1) then
       CALL INDEXX(n=SIZE(innerDiskData_Phi,1),ARRIN=innerDiskData_Phi,INDX=phiIndex)
    else
       phiIndex=1
    endif

    DO i = 1, SIZE(innerDiskData_Phi)
      Phi_temp(i) = innerDiskData_Phi(phiIndex(i))
      Radii_temp(i) = innerDiskData_Radii(phiIndex(i))
    END DO

    innerDiskData_Phi = Phi_temp 
    innerDiskData_Radii = Radii_temp

  END SUBROUTINE sortInnerDiskData

  SUBROUTINE rescaleHotSpot(hotspot, maxSpotTemp)
    ! makes some stellar hot-spots bigger, so that they are cooler

    TYPE(hotSpotVariable), INTENT(INOUT) :: hotspot
    REAL, INTENT(IN) :: maxSpotTemp ! (K)
    REAL :: TscaleFactor
    REAL :: RscaleFactor
    REAL :: RhoScaleFactor
    
    IF ( hotspot%temperature > maxSpotTemp ) THEN
      TscaleFactor = maxSpotTemp / hotspot%temperature
      RscaleFactor  = TscaleFactor**(0.25)
      RhoScaleFactor = RscaleFactor**(-2.0)
      hotspot%radius = hotspot%radius * RscaleFactor
      hotspot%radius_1e10 = hotspot%radius_1e10 * RscaleFactor
      hotspot%temperature = maxSpotTemp
      hotspot%rho = hotspot%rho * RhoScaleFactor
    END IF

  END SUBROUTINE rescaleHotSpot

END MODULE magField
