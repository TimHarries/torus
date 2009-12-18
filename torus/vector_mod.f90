module vector_mod
!
! This is a module for vector maths. This includes the vector
! type definitions, as well as simple operations such as add,
! substract, multiple and divide. The two vector functions
! dot product and cross product are also defined.
!

! written by tjh


! v1.0 on 13/08/99
! real(double) and 'octal' versions of some routines added by nhs

  use kind_mod

  implicit none

  public:: &
       modulus, &
       normalize, &
       rotateX, rotateY, rotateZ, &
       intersectionLinePlane, &
       intersectionLineSphere, &
       arbitraryRotate, &
       getPolar, &
       fromPhotosphereVector, &
       distancePointLineSegment, &
       greatCircleDistance       
       
!  private :: &
!       rmultSingle, rmultDouble, rmultOctal, rmultSingleReversed, rmultDoubleReversed, &
!       rmultOctalReversed, dmult, divideVecSingle, divideVecDouble, divideVecOctal, &
!       addSingle, addDouble, addOctal, addOctalDouble,  addDoubleOctal, addOctalSingle, &
!       addSingleOctal, subtractSingle, subtractDouble, subtractOctal, &
!       singleToDoubleVector, doubleToSingleVector, singleToOctalVector, octalToSingleVector, &
!       doubleToOctalVector, octalToDoubleVector, dotProdSingle, dotProdDouble, dotProdOctal, &
!       dotProdOctalSingle, dotProdSingleOctal, crossProdSingle, crossProdDouble, crossProdOctal, &
!       crossProdOctalSingle, crossProdSingleOctal, modulusSingle, modulusDouble, modulusOctal, &
!       normalizeSingle, normalizeDouble, normalizeOctal, rotateXSingle, rotateXDouble, rotateXOctal, &
!       rotateYSingle, rotateYDouble, rotateYOctal, rotateZSingle, rotateZDouble, rotateZOctal, &
!       intersectionLinePlaneSingle, intersectionLinePlaneDouble, intersectionLinePlaneOctal, &
!       intersectionLineSphereSingle, intersectionLineSphereDouble, intersectionLineSphereOctal, &
!       arbitraryRotateSingle, arbitraryRotateDouble, arbitraryRotateOctal, &
!       getPolarSingle, getPolarDouble, getPolarOctal, fromPhotosphereSingleVector, &
!       fromPhotosphereDoubleVector, fromPhotosphereOctalVector, distancePointLineSegmentSingle, &
!       distancePointLineSegmentOctal, distancePointLineSegmentDouble, greatCircleDistanceSingle, &
!       greatCircleDistanceDouble, greatCircleDistanceOctal
  
  

  ! The definition of the vector type

  type VECTOR
     real(double) :: x
     real(double) :: y
     real(double) :: z
  end type VECTOR

  type spVector
    ! spherical polar vector
    real(double) :: r
    real(double) :: theta
    real(double) :: phi
  end type spVector

  ! Define multiply

  interface operator(*)
     module procedure rmult
     module procedure rmultReversed
  end interface

  ! and divide

  interface operator(/)
     module procedure divide
  end interface

  ! add

  interface operator(+)
     module procedure add
  end interface

  ! subtract

  interface operator(-)
     module procedure subtract
  end interface

  ! assignment

  interface assignment(=)
     module procedure SPtoXYZvector
  end interface

  ! equivalence

  interface operator(==)
     module procedure VectorEquivalence
  end interface

  ! dot product

  interface operator(.dot.)
     module procedure dotProd
  end interface

  ! cross product

  interface operator(.cross.)
     module procedure crossProd
  end interface
  
  type(VECTOR), parameter :: xHat = VECTOR(1.d0, 0.d0, 0.d0)
  type(VECTOR), parameter :: yHat = VECTOR(0.d0, 1.d0, 0.d0)
  type(VECTOR), parameter :: zHat = VECTOR(0.d0, 0.d0, 1.d0)
  
  type(VECTOR), parameter :: zeroVector = VECTOR(0.d0, 0.d0, 0.d0)
  
contains

  ! the dot product function

  real(double) pure function dotProd(a , b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    dotProd = a%x*b%x + a%y*b%y + a%z*b%z

  end function dotProd
  
  ! the cross product function

  type(VECTOR) pure function crossProd(a ,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    crossProd%x =  (a%y*b%z - a%z*b%y)
    crossProd%y = -(a%x*b%z - a%z*b%x)
    crossProd%z =  (a%x*b%y - a%y*b%x)
  end function crossProd
  
  ! normalization subroutine - checks for zero vector

  subroutine normalize(a)
    type(VECTOR), intent(inout) :: a
    type(VECTOR) :: ZeroVec = VECTOR(0.d0,0.d0,0.d0)
    real(double) :: m, oneOverM
    logical, save :: firstTime = .true.

    if (a .eq. ZeroVec) then
       if (firstTime) then
          write(*,'(a)') "! Attempt to normalize the zero vector"
          firstTime = .false.
       endif
       a = VECTOR(1.d0,0.d0,0.d0)
    else
       m = modulus(a)
       OneOverm = 1.d0 / m

       a = a * OneOverM
    endif

  end subroutine normalize
  
  ! find the modulus of a vector

  real(double) pure function modulus(a)
    type(VECTOR), intent(in) :: a

    modulus = a%x*a%x + a%y*a%y + a%z*a%z
    modulus = sqrt(modulus)

  end function modulus

  ! multiply function

  type(VECTOR) pure function rmult(a,b)
    real(double), intent(in) :: a
    type(VECTOR), intent(in) :: b

    rmult%x = a * b%x
    rmult%y = a * b%y
    rmult%z = a * b%z

  end function rmult

  type(VECTOR) pure function rmultReversed(b,a)
    real(double), intent(in) :: a
    type(VECTOR), intent(in) :: b

    rmultReversed%x = a * b%x
    rmultReversed%y = a * b%y
    rmultReversed%z = a * b%z

  end function rmultReversed
  ! divide vector by a scalar

  type(VECTOR) pure function divide(a,b)
    type(VECTOR), intent(in) :: a
    real(double), intent(in) :: b

    divide%x = a%x / b
    divide%y = a%y / b
    divide%z = a%z / b

  end function divide
  ! add two vectors
  
  type(VECTOR) pure function add(a,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    add%x = a%x + b%x
    add%y = a%y + b%y
    add%z = a%z + b%z

  end function add
  
  ! subtract two vectors

  type(VECTOR) pure function subtract(a,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    subtract%x = a%x - b%x
    subtract%y = a%y - b%y
    subtract%z = a%z - b%z

  end function subtract

  pure subroutine getPolar(vec, r, theta, phi)
    use constants_mod, only: twoPi
    implicit none
    type(vector), intent(in) :: vec
    real(double), intent(out) :: r, theta, phi
    real(double) :: cosTheta 

    r = modulus(vec)
    if ((vec%y == 0.d0) .and. (vec%x == 0.d0)) then
       phi = 0.d0
    else
       phi = atan2(vec%y, vec%x)
    endif
    if (phi < 0.d0) phi = phi + twoPi
    if (r /= 0.d0) then
       cosTheta = vec%z/r
    else
       cosTheta = 0.d0
    endif
       
    theta = acos(cosTheta)
  end subroutine getPolar
  
  ! rotate a vector "a" about the z-axis by angle b - CLOCKWISE!

  type(VECTOR) pure function rotateZ(a,b)
    type(VECTOR), intent(in) :: a
    real(double), intent(in) :: b   ! angle in radians
    real(double) :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateZ%x = cosb * a%x + sinb * a%y
    rotateZ%y =-sinb * a%x + cosb * a%y
    rotateZ%z = a%z

  end function rotateZ
  
  ! rotate a vector "a" about the x-axis by angle b - CLOCKWISE!

  type(VECTOR) pure function rotateX(a,b)
    type(VECTOR), intent(in) :: a
    real(double), intent(in) :: b   ! angle in radians
    real(double) :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateX%x = a%x 
    rotateX%y = cosb * a%y + sinb * a%z
    rotateX%z =-sinb * a%y + cosb * a%z

  end function rotateX
  
  ! rotate a vector "a" about the y-axis by angle b - ANTI-CLOCKWISE!

  type(VECTOR) pure function rotateY(a,b)
    type(VECTOR), intent(in) :: a
    real(double), intent(in) :: b   ! angle in radians
    real(double) :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateY%x = cosb * a%x + sinb * a%z
    rotateY%y = a%y
    rotateY%z =-sinb * a%x + cosb * a%z

  end function rotateY

!thap
  subroutine melvinUnitVector(rVec, photonPacketWeight)
    use constants_mod, only: pi, degToRad
    real(double) :: r, u, v, w, t, ang, openingAngle = 20.*degToRad, percentPhoCav
    real(double), intent(out) :: photonPacketWeight
    type(VECTOR), intent(out) :: rVec

    percentPhoCav = 0.9d0

    call random_number(r)
    if (r < percentPhoCav) then   !% of photons in cavities 
       call random_number(r)
       if(r < 0.5) then
          call random_number(r)
          w = cos(openingAngle) + (r * (1.0- cos(openingAngle)))
       else
          call random_number(r)
          w = -cos(openingAngle) + (r*(-1.0 + cos(openingAngle)))
       end if
       photonPacketWeight = (1.0 - cos(openingAngle)) / percentPhoCav          
    else  !rest towards disc
       call random_number(r)
       w = -cos(openingAngle) + (r * 2.0 * cos(openingAngle))
       photonPacketWeight = cos(openingAngle) / (1.0 - percentPhoCav)
   end if

    t = sqrt(1.0-w*w)
    call random_number(r)
    ang = pi*(2.d0*r-1.d0)
    u = t*cos(ang)
    v = t*sin(ang)

    rVec = VECTOR(u,v,w)
  end subroutine melvinUnitVector


  type(VECTOR) function randomUnitVector()
    use constants_mod, only: pi
    real(double) :: r1, r2, u, v, w, t, ang
    call random_number(r1)
    w = 2.*r1 - 1.
    t = sqrt(1.0-w*w)
    call random_number(r2)
    ang = pi*(2.d0*r2-1.d0)
    u = t*cos(ang)
    v = t*sin(ang)

    randomUnitVector = VECTOR(u,v,w)
  end function randomUnitVector

  type(VECTOR) function specificUnitVector(r1,r2)
    use constants_mod, only: pi
    real(double) :: r1, r2, u, v, w, t, ang
!    type(VECTOR) :: specificunitvector
    w = 2.d0*r1 - 1.
    t = sqrt(1.d0-w*w)
    ang = pi*(2.d0*r2-1.d0)
    u = t*cos(ang)
    v = t*sin(ang)

    specificUnitVector = VECTOR(u,v,w)
  end function specificUnitVector

  type (VECTOR) function intersectionLinePlane(r0, rHat, nHat, d, ok)

! finds the intersection between a line and a plane

    implicit none

    type(VECTOR), intent(in) :: r0, rHat ! equation of line
    type(VECTOR), intent(in) :: nHat     ! the normal to the plane
    real(double), intent(in) :: d                ! minimum distance of plane from origin
    logical, intent(out) :: ok            ! is there in intersection?
    real(double) :: fac

    ok = .false.
    if ((nHat.dot.rHat) /= 0.) then
       fac = ((d - (nHat.dot.r0))/(nHat.dot.rHat))
       if (fac > 0.d0) then
          intersectionLinePlane = r0 + (fac  * rHat )
          ok = .true.
       else
          intersectionLinePlane = VECTOR(0.d0,0.d0,0.d0)
          ok = .false.
       endif
    else
       intersectionLinePlane = VECTOR(0.d0,0.d0,0.d0)
       ok = .false.
    endif

  end function intersectionLinePlane
  
  pure subroutine intersectionLineSphere(r0, rHat, length, s0, sR, found1,&
                              found2,intersectionDistance1,intersectionDistance2)
    ! finds the intersection(s) between a line segment and a sphere
    
    implicit none

    type(vector), intent(in)      :: r0, rHat ! equation of line
    real(double), intent(in)              :: length   ! length of line segment 
    type(vector), intent(in)      :: s0       ! centre of sphere
    real(double), intent(in)              :: sR       ! radius of sphere
    logical,intent(out)           :: found1   ! there is one intersection
    logical,intent(out)           :: found2   ! there is a second intersection
    real(double), intent(out)             :: intersectionDistance1 ! distances along line 
    real(double), intent(out)             :: intersectionDistance2 !   to intersections
    real(double)  :: a, b, c       ! quadratic formula variables
    real(double)  :: discriminant  ! quadratic formula variable
    real(double)  :: solution1     ! quadratic formula variable
    real(double)  :: solution2     ! quadratic formula variable
    
    type(Vector) :: r0Double, rHatDouble ! equation of line (real(double))
    type(Vector) :: s0Double             ! centre of sphere (real(double))
    real(double)   :: sRDouble             ! radius of sphere (real(double))
    
    r0Double = Vector(r0%x,r0%y,r0%z)
    rHatDouble = Vector(rHat%x,rHat%y,rHat%z)
    s0Double = Vector(s0%x,s0%y,s0%z)
    sRDouble = sR
    
    a = rHatDouble%x**2.0_db + rHatDouble%y**2.0_db + rHatDouble%z**2.0_db
    b = 2.0_db * ( rHatDouble%x * (r0Double%x - s0Double%x) + rHatDouble%y * &
                (r0Double%y - s0Double%y) + rHatDouble%z * (r0Double%z - s0Double%z)) 
    c = (r0Double%x - s0Double%x)**2.0_db + (r0Double%y - s0Double%y)**2.0_db + &
                (r0Double%z - s0Double%z)**2.0_db - sRDouble**2.0_db 
 
    discriminant = b**2.0_db - (4.0_db * a * c)

    intersectionDistance1 = -9.99e9
    intersectionDistance2 = -9.99e9
    found1 = .false.
    found2 = .false.
    
    if (discriminant >= 0.0_db) then
      
      ! we only need to find the smallest positive solution of the quadratic
      
      solution1 = (-b - (b**2.0_db - 4.0_db*a*c)**0.5) / (2.0_db * a) 
      solution2 = (-b + (b**2.0_db - 4.0_db*a*c)**0.5) / (2.0_db * a) 
      
      if (solution1 >= 0.0_db) then 
        
        if (solution1 < length) then
          found1 = .true.
          intersectionDistance1 = solution1
        end if
        
        if (solution2 < length) then
          if (found1) then      
            found2 = .true.
            intersectionDistance2 = solution2
          else
            found1 = .true.
            intersectionDistance1 = solution2
          end if
        end if
      end if
    end if
    
  end subroutine intersectionLineSphere
  
  type(VECTOR) pure function arbitraryRotate(p, theta, r)
    type(VECTOR),intent(in) :: p      ! position vector
    real(double),intent(in)         :: theta  ! angle in radians
    type(VECTOR),intent(in) :: r      ! the arbitrary axis
    
    real(double) :: cosTheta, sinTheta
    type(VECTOR)            :: q    

    costheta = cos(theta)
    sintheta = sin(theta)
    
    q%x = (costheta + (1. - costheta) * r%x * r%x) * p%x
    q%x = q%x + ((1. - costheta) * r%x * r%y - r%z * sintheta) * p%y
    q%x = q%x + ((1. - costheta) * r%x * r%z + r%y * sintheta) * p%z
    
    q%y = ((1. - costheta) * r%x * r%y + r%z * sintheta) * p%x
    q%y = q%y + (costheta + (1. - costheta) * r%y * r%y) * p%y
    q%y = q%y + ((1. - costheta) * r%y * r%z - r%x * sintheta) * p%z
    
    q%z = ((1. - costheta) * r%x * r%z - r%y * sintheta) * p%x
    q%z = q%z + ((1. - costheta) * r%y * r%z + r%x * sintheta) * p%y
    q%z = q%z + (costheta + (1. - costheta) * r%z * r%z) * p%z
    
    arbitraryRotate = q
    
  end function arbitraryRotate
  
  
  type(VECTOR) function fromPhotosphereVector(rVec)
    use constants_mod, only: twoPi
    real(double) :: ang, z
    type(VECTOR) :: norm, zAxis, v, rVec, n
    
    zAxis = VECTOR(0.,0.,1.)

    do while (ABS(rVec .dot. zAxis) == 1.0d0)
       ! choose another one
       rVec = randomUnitVector()
    end do
    
    norm = rVec
    call normalize(norm)

    call random_number(z)
    z = sqrt(z)

    ang = acos(z)
    n = norm .cross. zAxis
    call normalize(n)
    v = norm

    v = arbitraryRotate(v, ang, n)
    call random_number(ang)
    ang = ang * twoPi
    v = arbitraryRotate(v, ang, norm)

    fromPhotosphereVector = v
  end function fromPhotosphereVector

  function projectToXZ(rVec) result (out)
    implicit none
    type(VECTOR) :: OUT
    type(VECTOR), save :: outprev, rvecprev
    type(VECTOR),intent(in) :: rVec

    if(rvec .eq. rvecprev) then
       out = outprev
    else
       out%x = sqrt(rVec%x**2 + rVec%y**2)
       out%y = 0.d0
       out%z = rVec%z

       outprev = out
       rvecprev = rvec
    endif
       
  end function projectToXZ

  elemental subroutine SPtoXYZvector(out,in)

    TYPE(vector), INTENT(OUT) :: out
    TYPE(spVector), INTENT(IN)  :: in

    out%x = in%r * SIN(in%theta) * COS(in%phi)
    out%y = in%r * SIN(in%theta) * SIN(in%phi)
    out%z = in%r * COS(in%theta)
    
  end subroutine SPtoXYZvector

  ELEMENTAL SUBROUTINE XYZtoSPvector(out,in)
    use constants_mod, only: twoPi, pi
    TYPE(SPvector), INTENT(OUT) :: out
    TYPE(vector),  INTENT(IN) :: in

    out%r = modulus(in)
    
    IF ( ABS(out%r) > TINY(1.0) ) THEN
      out%theta = ACOS( in%z / out%r )
    ELSE   
      out%theta = pi / 2.0
    END IF
      
    IF ( ( ABS(in%x) > TINY(1.0) ) .AND. ( ABS(in%y) > TINY(1.0) )) THEN
      out%phi = ATAN2(in%y, in%x) 
    ELSE
      out%phi = 0.0
    ENDIF
    ! transform range from (-pi to pi) to (0 to 2pi)
    IF (out%phi < 0.0) out%phi = out%phi + twoPi
    
  END SUBROUTINE XYZtoSPvector

  elemental logical function VectorEquivalence(a, b)
    TYPE(vector),  INTENT(IN) :: a,b

    VectorEquivalence = .false.

    if(a%x .ne. b%x) then
       return
    else
       if(a%y .ne. b%y) then
          return
       else
          if(a%z .eq. b%z) then
             VectorEquivalence = .true.
             return
          endif
       endif
    endif
    
  end function VectorEquivalence
  
  FUNCTION greatCircleDistance(pointA, pointB) RESULT(distance)
    ! returns the great circle distance between two points
    !   on a spherical surface
 
    REAL :: distance ! (radians)
    TYPE(spVector), INTENT(IN) :: pointA, pointB
      ! the "r" components of the vectors are ignored
    
    REAL(double) :: dLon, dLat, a
    
    dLon = pointB%phi - pointA%phi
    dLat = pointB%theta - pointA%theta
    a = (SIN(dLat/2.))**2 + COS(pointA%theta) * COS(pointB%theta) * SIN(dLon/2.)**2
    distance = 2. * ASIN(MIN(1.,SQRT(a)))
  
  END FUNCTION greatCircleDistance
  
  
  FUNCTION distancePointLineSegment(linePoint1, linePoint2, testPoint) &
                                          RESULT(distance)
    ! given a line segment, defined by the two points at either end, 
    !   this returns the distance to the line segment from another point
   
    REAL(double) :: distance
    TYPE(vector), INTENT(IN) :: linePoint1, linePoint2 ! points defining line segment
    TYPE(vector), INTENT(IN) :: testPoint

    TYPE(vector) :: line ! line segment
    TYPE(vector) :: externalLine ! line to external point
    REAL(double) :: a, b, c
    REAL(double) :: distanceA, distanceB
    TYPE(vector) :: extendedPoint

    line = linePoint2 - linePoint1
    externalLine = testPoint - linePoint1

    a = line .dot. externalLine
    distanceA = modulus( testPoint - linePoint1 )
    IF ( a <= 0.0 ) THEN
      distance = distanceA
      RETURN
    END IF
    
    b = line .dot. line
    distanceB = modulus( testPoint - linePoint2 )
    IF ( b <= a ) THEN
      distance = distanceB
      RETURN
    END IF

    c = a / b
    extendedPoint = linePoint1 + (c * line)
    distance = modulus( testPoint - extendedPoint )
    
  END FUNCTION distancePointLineSegment

  function minimumDistanceFromPointToLine(point, lineStart, lineEnd) result (distance)
    real(double) :: distance
    type(VECTOR) :: point, lineStart,lineEnd, a, b

    a = lineEnd - lineStart
    b = lineStart - point

    distance = modulus(a .cross. b) / modulus(a)
  end function minimumDistanceFromPointToLine
    


end module vector_mod

