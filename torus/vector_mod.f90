!
! This is a module for vector maths. This includes the vector
! type definitions, as well as simple operations such as add,
! substract, multiple and divide. The two vector functions
! dot product and cross product are also defined.
!

! written by tjh


! v1.0 on 13/08/99
! double precision and 'octal' versions of some routines added by nhs

module vector_mod


  use kind_mod
  use constants_mod

  implicit none

  public

  ! The definition of the vector type

  type VECTOR
     real(KIND=singleKind) :: x
     real(KIND=singleKind) :: y
     real(KIND=singleKind) :: z
  end type VECTOR

  TYPE doubleVector
    REAL(KIND=doubleKind) :: x
    REAL(KIND=doubleKind) :: y
    REAL(KIND=doubleKind) :: z
  END TYPE doubleVector
  
  TYPE octalVector
    REAL(KIND=octalKind) :: x
    REAL(KIND=octalKind) :: y
    REAL(KIND=octalKind) :: z
  END TYPE octalVector

  ! Define multiply

  interface operator(*)
     module procedure rmultSingle
     module procedure rmultDouble
     module procedure rmultOctal
     module procedure rmultSingleReversed
     module procedure rmultDoubleReversed
     module procedure rmultOctalReversed
     module procedure dmult ! see dmult comments
  end interface

  ! and divide

  interface operator(/)
     module procedure divideVecSingle
     module procedure divideVecDouble
     module procedure divideVecOctal
  end interface

  ! add

  interface operator(+)
     module procedure addSingle
     module procedure addDouble
     module procedure addOctal
  end interface

  ! subtract

  interface operator(-)
     module procedure subtractSingle
     module procedure subtractDouble
     module procedure subtractOctal
  end interface

  ! dot product

  interface operator(.dot.)
     module procedure dotProdSingle
     module procedure dotProdDouble
     module procedure dotProdOctal
  end interface

  ! cross product

  interface operator(.cross.)
     module procedure crossProdSingle
     module procedure crossProdDouble
     module procedure crossProdOctal
  end interface
  
  INTERFACE modulus
    MODULE PROCEDURE modulusSingle
    MODULE PROCEDURE modulusDouble
    MODULE PROCEDURE modulusOctal
  END INTERFACE

  INTERFACE normalize
    MODULE PROCEDURE normalizeSingle
    MODULE PROCEDURE normalizeDouble
    MODULE PROCEDURE normalizeOctal
  END INTERFACE

  INTERFACE rotateX
    MODULE PROCEDURE rotateXSingle
    MODULE PROCEDURE rotateXDouble
    MODULE PROCEDURE rotateXOctal
  END INTERFACE

  INTERFACE rotateY
    MODULE PROCEDURE rotateYSingle
    MODULE PROCEDURE rotateYDouble
    MODULE PROCEDURE rotateYOctal
  END INTERFACE

  INTERFACE rotateZ
    MODULE PROCEDURE rotateZSingle
    MODULE PROCEDURE rotateZDouble
    MODULE PROCEDURE rotateZOctal
  END INTERFACE

  INTERFACE intersectionLinePlane
    MODULE PROCEDURE intersectionLinePlaneSingle
    MODULE PROCEDURE intersectionLinePlaneDouble
    MODULE PROCEDURE intersectionLinePlaneOctal
  END INTERFACE
  
  INTERFACE arbitraryRotate
    MODULE PROCEDURE arbitraryRotateSingle
    MODULE PROCEDURE arbitraryRotateDouble
    MODULE PROCEDURE arbitraryRotateOctal
  END INTERFACE
  
  INTERFACE getPolar
    MODULE PROCEDURE getPolarSingle
    MODULE PROCEDURE getPolarDouble
    MODULE PROCEDURE getPolarOctal
  END INTERFACE
  
  type(VECTOR), parameter :: xHat = VECTOR(1., 0., 0.)
  type(VECTOR), parameter :: yHat = VECTOR(0., 1., 0.)
  type(VECTOR), parameter :: zHat = VECTOR(0., 0., 1.)
  
  TYPE(doubleVector),PARAMETER :: xHatDouble=doubleVector(1.0_db,0.0_db,0.0_db)
  TYPE(doubleVector),PARAMETER :: yHatDouble=doubleVector(0.0_db,1.0_db,0.0_db)
  TYPE(doubleVector),PARAMETER :: zHatDouble=doubleVector(0.0_db,0.0_db,1.0_db)

  TYPE(octalVector),PARAMETER :: xHatOctal=octalVector(1.0_oc,0.0_oc,0.0_oc)
  TYPE(octalVector),PARAMETER :: yHatOctal=octalVector(0.0_oc,1.0_oc,0.0_oc)
  TYPE(octalVector),PARAMETER :: zHatOctal=octalVector(0.0_oc,0.0_oc,1.0_oc)

contains

  ! the dot product function

  real function dotProdSingle(a , b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    dotProdSingle = a%x*b%x + a%y*b%y + a%z*b%z

  end function dotProdSingle
  
  REAL(KIND=doubleKind) FUNCTION dotProdDouble(a , b)
    TYPE(doubleVector), INTENT(IN) :: a
    TYPE(doubleVector), INTENT(IN) :: b

    dotProdDouble = a%x*b%x + a%y*b%y + a%z*b%z

  END FUNCTION dotProdDouble
  
  REAL(KIND=octalKind) FUNCTION dotProdOctal(a , b)
    TYPE(octalVector), INTENT(IN) :: a
    TYPE(octalVector), INTENT(IN) :: b

    dotProdOctal = a%x*b%x + a%y*b%y + a%z*b%z

  END FUNCTION dotProdOctal


  ! the cross product function

  type(VECTOR) function crossProdSingle(a ,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    crossProdSingle%x =  (a%y*b%z - a%z*b%y)
    crossProdSingle%y = -(a%x*b%z - a%z*b%x)
    crossProdSingle%z =  (a%x*b%y - a%y*b%x)
  end function crossProdSingle
  
  TYPE(doubleVector) FUNCTION crossProdDouble(a,b)
    TYPE(doubleVector), INTENT(IN) :: a
    TYPE(doubleVector), INTENT(IN) :: b

    crossProdDouble%x =  (a%y*b%z - a%z*b%y)
    crossProdDouble%y = -(a%x*b%z - a%z*b%x)
    crossProdDouble%z =  (a%x*b%y - a%y*b%x)
  END FUNCTION crossProdDouble

  TYPE(octalVector) FUNCTION crossProdOctal(a,b)
    TYPE(octalVector), INTENT(IN) :: a
    TYPE(octalVector), INTENT(IN) :: b

    crossProdOctal%x =  (a%y*b%z - a%z*b%y)
    crossProdOctal%y = -(a%x*b%z - a%z*b%x)
    crossProdOctal%z =  (a%x*b%y - a%y*b%x)
  END FUNCTION crossProdOctal


  ! normalization subroutine - checks for zero vector

  subroutine normalizeSingle(a)
    type(VECTOR), intent(inout) :: a
    real :: m

    m = modulus(a)

    if (m == 0.) then
       write(*,'(a)') "! Attempt to normalize the zero vector"
       m = sqrt(m - 1.d0)
       a = VECTOR(1.,0.,0.)
    endif

    a%x = a%x / m
    a%y = a%y / m
    a%z = a%z / m

  end subroutine normalizeSingle
  
  SUBROUTINE normalizeDouble(a)
    TYPE(doubleVector), INTENT(INOUT) :: a
    REAL(KIND=doubleKind) :: m

    m = modulus(a)

    IF (m == 0.0_db) THEN
       WRITE(*,'(a)') "! Attempt to normalize the zero vector"
       m = SQRT(m - 1.0_db)
       a = doubleVector(1.0_db,0.0_db,0.0_db)
    endif

    a%x = a%x / m
    a%y = a%y / m
    a%z = a%z / m

  END SUBROUTINE normalizeDouble

  SUBROUTINE normalizeOctal(a)
    TYPE(octalVector), INTENT(INOUT) :: a
    REAL(KIND=octalKind) :: m

    m = modulus(a)

    IF (m == 0.0_oc) THEN
       WRITE(*,'(a)') "! Attempt to normalize the zero vector"
       m = SQRT(m - 1.0_oc)
       a = octalVector(1.0_oc,0.0_oc,0.0_oc)
    endif

    a%x = a%x / m
    a%y = a%y / m
    a%z = a%z / m

  END SUBROUTINE normalizeOctal


  ! find the modulus of a vector

  real function modulusSingle(a)
    type(VECTOR) :: a

    modulusSingle = a%x*a%x + a%y*a%y + a%z*a%z
    modulusSingle = sqrt(modulusSingle)

  end function modulusSingle
  
  REAL(KIND=doubleKind) FUNCTION modulusDouble(a)
    TYPE(doubleVector) :: a

    modulusDouble = a%x*a%x + a%y*a%y + a%z*a%z
    modulusDouble = SQRT(modulusDouble)

  END FUNCTION modulusDouble
  
  REAL(KIND=octalKind) FUNCTION modulusOctal(a)
    TYPE(octalVector) :: a

    modulusOctal = a%x*a%x + a%y*a%y + a%z*a%z
    modulusOctal = SQRT(modulusOctal)

  END FUNCTION modulusOctal

  ! multiply function

  type(VECTOR) function rmultSingle(a,b)
    real, intent(in) :: a
    type(VECTOR), intent(in) :: b

    rmultSingle%x = a * b%x
    rmultSingle%y = a * b%y
    rmultSingle%z = a * b%z

  end function rmultSingle

  TYPE(doubleVector) function rmultDouble(a,b)
    real(kind=doubleKind), intent(in) :: a
    TYPE(doubleVector), intent(in) :: b

    rmultDouble%x = a * b%x
    rmultDouble%y = a * b%y
    rmultDouble%z = a * b%z

  end function rmultDouble

  TYPE(octalVector) function rmultOctal(a,b)
    real(kind=octalKind), intent(in) :: a
    TYPE(octalVector), intent(in) :: b

    rmultOctal%x = a * b%x
    rmultOctal%y = a * b%y
    rmultOctal%z = a * b%z

  end function rmultOctal

  type(VECTOR) function rmultSingleReversed(b,a)
    real, intent(in) :: a
    type(VECTOR), intent(in) :: b

    rmultSingleReversed%x = a * b%x
    rmultSingleReversed%y = a * b%y
    rmultSingleReversed%z = a * b%z

  end function rmultSingleReversed

  TYPE(doubleVector) function rmultDoubleReversed(b,a)
    real(kind=doubleKind), intent(in) :: a
    TYPE(doubleVector), intent(in) :: b

    rmultDoubleReversed%x = a * b%x
    rmultDoubleReversed%y = a * b%y
    rmultDoubleReversed%z = a * b%z

  end function rmultDoubleReversed

  TYPE(octalVector) function rmultOctalReversed(b,a)
    real(kind=octalKind), intent(in) :: a
    TYPE(octalVector), intent(in) :: b

    rmultOctalReversed%x = a * b%x
    rmultOctalReversed%y = a * b%y
    rmultOctalReversed%z = a * b%z

  end function rmultOctalReversed

 
  type(VECTOR) function dmult(a,b)
    ! this has one single precision and one double precision 
    !   argument, which doesn't really fit with the scheme used
    !   elsewhere in this module, but is required by some of the
    !   other torus routines.
    real(kind=doubleKind), intent(in) :: a
    type(VECTOR), intent(in) :: b

    dmult%x = a * b%x
    dmult%y = a * b%y
    dmult%z = a * b%z

  end function dmult


  ! divide vector by a scalar

  type(VECTOR) function divideVecSingle(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b

    divideVecSingle%x = a%x / b
    divideVecSingle%y = a%y / b
    divideVecSingle%z = a%z / b

  end function divideVecSingle
  
  TYPE(doubleVector) FUNCTION divideVecDouble(a,b)
    TYPE(doubleVector), INTENT(IN) :: a
    REAL(KIND=doubleKind), INTENT(IN) :: b

    divideVecDouble%x = a%x / b
    divideVecDouble%y = a%y / b
    divideVecDouble%z = a%z / b

  END FUNCTION divideVecDouble
  
  TYPE(octalVector) FUNCTION divideVecOctal(a,b)
    TYPE(octalVector), INTENT(IN) :: a
    REAL(KIND=octalKind), INTENT(IN) :: b

    divideVecOctal%x = a%x / b
    divideVecOctal%y = a%y / b
    divideVecOctal%z = a%z / b

  END FUNCTION divideVecOctal


  ! add two vectors
  
  type(VECTOR) function addSingle(a,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    addSingle%x = a%x + b%x
    addSingle%y = a%y + b%y
    addSingle%z = a%z + b%z

  end function addSingle
  
  TYPE(doubleVector) FUNCTION addDouble(a,b)
    TYPE(doubleVector), INTENT(IN) :: a
    TYPE(doubleVector), INTENT(IN) :: b

    addDouble%x = a%x + b%x
    addDouble%y = a%y + b%y
    addDouble%z = a%z + b%z

  END FUNCTION addDouble
  
  TYPE(octalVector) FUNCTION addOctal(a,b)
    TYPE(octalVector), INTENT(IN) :: a
    TYPE(octalVector), INTENT(IN) :: b

    addOctal%x = a%x + b%x
    addOctal%y = a%y + b%y
    addOctal%z = a%z + b%z

  END FUNCTION addOctal


  ! subtract two vectors

  type(VECTOR) function subtractSingle(a,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    subtractSingle%x = a%x - b%x
    subtractSingle%y = a%y - b%y
    subtractSingle%z = a%z - b%z

  end function subtractSingle
  
  TYPE(doubleVector) FUNCTION subtractDouble(a,b)
    TYPE(doubleVector), INTENT(IN) :: a
    TYPE(doubleVector), INTENT(IN) :: b

    subtractDouble%x = a%x - b%x
    subtractDouble%y = a%y - b%y
    subtractDouble%z = a%z - b%z

  END FUNCTION subtractDouble

  TYPE(octalVector) FUNCTION subtractOctal(a,b)
    TYPE(octalVector), INTENT(IN) :: a
    TYPE(octalVector), INTENT(IN) :: b

    subtractOctal%x = a%x - b%x
    subtractOctal%y = a%y - b%y
    subtractOctal%z = a%z - b%z

  END FUNCTION subtractOctal


  ! get polar form of a cartesian vector
  
  subroutine getPolarSingle(vec, r, theta, phi)

    IMPLICIT NONE
    TYPE(vector), INTENT(IN) :: vec
    REAL, INTENT(OUT) :: r, theta, phi
    REAL :: cosTheta 

    r = modulus(vec)
    if ((vec%y == 0.) .and. (vec%x == 0.)) then
       phi = 0.
    else
       phi = atan2(vec%y, vec%x)
    endif
    if (phi < 0.) phi = phi + twoPi
    if (r /= 0.) then
       cosTheta = vec%z/r
    else
       cosTheta = 0.
    endif
       
    theta = acos(cosTheta)
  end subroutine getPolarSingle
  
  SUBROUTINE getPolarDouble(vec, r, theta, phi)

    IMPLICIT NONE
    TYPE(doubleVector), INTENT(IN) :: vec
    REAL(KIND=doubleKind),INTENT(OUT) :: r, theta, phi
    REAL(KIND=doubleKind) :: cosTheta 

    r = modulus(vec)
    if ((vec%y == 0.0_db) .and. (vec%x == 0.0_db)) then
       phi = 0.0_db
    else
       phi = atan2(vec%y, vec%x)
    endif
    if (phi < 0.0_db) phi = phi + twoPi
    if (r /= 0.0_db) then
       cosTheta = vec%z/r
    else
       cosTheta = 0.0_db
    endif
       
    theta = acos(cosTheta)
  end subroutine getPolarDouble
  
  SUBROUTINE getPolarOctal(vec, r, theta, phi)

    IMPLICIT NONE
    TYPE(octalVector), INTENT(IN) :: vec
    REAL(KIND=octalKind),INTENT(OUT) :: r, theta, phi
    REAL(KIND=octalKind) :: cosTheta 

    r = modulus(vec)
    if ((vec%y == 0.0_oc) .and. (vec%x == 0.0_oc)) then
       phi = 0.0_oc
    else
       phi = atan2(vec%y, vec%x)
    endif
    if (phi < 0.0_oc) phi = phi + twoPi
    if (r /= 0.0_oc) then
       cosTheta = vec%z/r
    else
       cosTheta = 0.0_oc
    endif
       
    theta = acos(cosTheta)
  end subroutine getPolarOctal


  ! rotate a vector "a" about the z-axis by angle b

  type(VECTOR) function rotateZSingle(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateZSingle%x = cosb * a%x + sinb * a%y
    rotateZSingle%y =-sinb * a%x + cosb * a%y
    rotateZSingle%z = a%z

  end function rotateZSingle
  
  TYPE(doubleVector) FUNCTION rotateZDouble(a,b)
    TYPE(doubleVector), INTENT(IN) :: a
    REAL(KIND=doubleKind), INTENT(IN) :: b   ! angle in radians
    REAL(KIND=doubleKind) :: cosb, sinb

    cosb = COS(b)
    sinb = SIN(b)

    rotateZDouble%x = cosb * a%x + sinb * a%y
    rotateZDouble%y =-sinb * a%x + cosb * a%y
    rotateZDouble%z = a%z

  END FUNCTION rotateZDouble
  
  TYPE(octalVector) FUNCTION rotateZOctal(a,b)
    TYPE(octalVector), INTENT(IN) :: a
    REAL(KIND=octalKind), INTENT(IN) :: b   ! angle in radians
    REAL(KIND=octalKind) :: cosb, sinb

    cosb = COS(b)
    sinb = SIN(b)

    rotateZOctal%x = cosb * a%x + sinb * a%y
    rotateZOctal%y =-sinb * a%x + cosb * a%y
    rotateZOctal%z = a%z

  END FUNCTION rotateZOctal


  type(VECTOR) function rotateXSingle(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateXSingle%x = a%x 
    rotateXSingle%y = cosb * a%y + sinb * a%z
    rotateXSingle%z =-sinb * a%y + cosb * a%z

  end function rotateXSingle
  
  TYPE(doubleVector) FUNCTION rotateXDouble(a,b)
    TYPE(doubleVector), INTENT(IN) :: a
    REAL(KIND=doubleKind), INTENT(IN) :: b   ! angle in radians
    REAL(KIND=doubleKind)  :: cosb, sinb

    cosb = COS(b)
    sinb = SIN(b)

    rotateXDouble%x = a%x 
    rotateXDouble%y = cosb * a%y + sinb * a%z
    rotateXDouble%z =-sinb * a%y + cosb * a%z

  END FUNCTION rotateXDouble
  
  TYPE(octalVector) FUNCTION rotateXOctal(a,b)
    TYPE(octalVector), INTENT(IN) :: a
    REAL(KIND=octalKind), INTENT(IN) :: b   ! angle in radians
    REAL(KIND=octalKind)  :: cosb, sinb

    cosb = COS(b)
    sinb = SIN(b)

    rotateXOctal%x = a%x 
    rotateXOctal%y = cosb * a%y + sinb * a%z
    rotateXOctal%z =-sinb * a%y + cosb * a%z

  END FUNCTION rotateXOctal


  type(VECTOR) function rotateYSingle(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateYSingle%x = cosb * a%x + sinb * a%z
    rotateYSingle%y = a%y
    rotateYSingle%z =-sinb * a%x + cosb * a%z

  end function rotateYSingle
  
  TYPE(doubleVector) FUNCTION rotateYDouble(a,b)
    TYPE(doubleVector), INTENT(IN) :: a
    REAL(KIND=doubleKind), INTENT(IN) :: b   ! angle in radians
    REAL(KIND=doubleKind) :: cosb, sinb

    cosb = COS(b)
    sinb = SIN(b)

    rotateYDouble%x = cosb * a%x + sinb * a%z
    rotateYDouble%y = a%y
    rotateYDouble%z =-sinb * a%x + cosb * a%z

  END FUNCTION rotateYDouble
  
  TYPE(octalVector) FUNCTION rotateYOctal(a,b)
    TYPE(octalVector), INTENT(IN) :: a
    REAL(KIND=octalKind), INTENT(IN) :: b   ! angle in radians
    REAL(KIND=octalKind) :: cosb, sinb

    cosb = COS(b)
    sinb = SIN(b)

    rotateYOctal%x = cosb * a%x + sinb * a%z
    rotateYOctal%y = a%y
    rotateYOctal%z =-sinb * a%x + cosb * a%z

  END FUNCTION rotateYOctal


  type(VECTOR) function randomUnitVector()
    real :: r1, r2, u, v, w, t, ang
    call random_number(r1)
    w = 2.*r1 - 1.
    t = sqrt(1.d0-w*w)
    call random_number(r2)
    ang = pi*(2.*r2-1.)
    u = t*cos(ang)
    v = t*sin(ang)

    randomUnitVector = VECTOR(u,v,w)
  end function randomUnitVector


  type (VECTOR) function intersectionLinePlaneSingle(r0, rHat, nHat, d, ok)

! finds the intersection between a line and a plane

    type(VECTOR) :: r0, rHat    ! equation of line
    type(VECTOR) :: nHat        ! the normal to the plane
    real :: d                   ! minimum distance of plane from origin
    logical :: ok               ! is there in intersection?
    real :: fac

    ok = .false.
    if ((nHat.dot.rHat) /= 0.) then
       fac = ((d - (nHat.dot.r0))/(nHat.dot.rHat))
       if (fac > 0.) then
          intersectionLinePlaneSingle = r0 + (fac  * rHat )
          ok = .true.
       else
          intersectionLinePlaneSingle = VECTOR(0.,0.,0.)
          ok = .false.
       endif
    else
       intersectionLinePlaneSingle = VECTOR(0.,0.,0.)
       ok = .false.
    endif

  end function intersectionLinePlaneSingle
  
  TYPE (doubleVector) FUNCTION intersectionLinePlaneDouble &
                                           (r0, rHat, nHat, d, ok)
    ! finds the intersection between a line and a plane

    TYPE(doubleVector), INTENT(IN) :: r0, rHat    ! equation of line
    TYPE(doubleVector), INTENT(IN) :: nHat        ! the normal to the plane
    REAL(KIND=doubleKind) , INTENT(IN) :: d
                         ! minimum distance of plane from origin
    LOGICAL,INTENT(OUT) :: ok               ! is there in intersection?
    REAL(KIND=doubleKind) :: fac

    ok = .FALSE.
    IF ((nHat.dot.rHat) /= 0.0_db) THEN
       fac = ((d - (nHat.dot.r0))/(nHat.dot.rHat))
       IF (fac > 0.0_db) THEN
          intersectionLinePlaneDouble = r0 + (fac * rHat)
          ok = .TRUE.
       ELSE
          intersectionLinePlaneDouble = doubleVector(0.0_db, 0.0_db, 0.0_db)
          ok = .FALSE.
       ENDIF
    ELSE
       intersectionLinePlaneDouble = doubleVector(0.0_db ,0.0_db ,0.0_db)
       ok = .false.
    ENDIF

  END FUNCTION intersectionLinePlaneDouble
  
  TYPE (octalVector) FUNCTION intersectionLinePlaneOctal &
                                           (r0, rHat, nHat, d, ok)
    ! finds the intersection between a line and a plane

    TYPE(octalVector), INTENT(IN) :: r0, rHat    ! equation of line
    TYPE(octalVector), INTENT(IN) :: nHat        ! the normal to the plane
    REAL(KIND=octalKind) , INTENT(IN) :: d
                         ! minimum distance of plane from origin
    LOGICAL,INTENT(OUT) :: ok               ! is there in intersection?
    REAL(KIND=octalKind) :: fac

    ok = .FALSE.
    IF ((nHat.dot.rHat) /= 0.0_oc) THEN
       fac = ((d - (nHat.dot.r0))/(nHat.dot.rHat))
       IF (fac > 0.0_oc) THEN
          intersectionLinePlaneOctal = r0 + (fac * rHat)
          ok = .TRUE.
       ELSE
          intersectionLinePlaneOctal = octalVector(0.0_oc, 0.0_oc, 0.0_oc)
          ok = .FALSE.
       ENDIF
    ELSE
       intersectionLinePlaneOctal = octalVector(0.0_oc ,0.0_oc ,0.0_oc)
       ok = .false.
    ENDIF

  END FUNCTION intersectionLinePlaneOctal


  type(VECTOR) function arbitraryRotateSingle(p, theta, r)
    type(VECTOR) :: p,q    ! position vector
    real :: theta          ! angle in radians
    type(VECTOR) :: r      ! the arbitrary axis
    real :: cosTheta, sinTheta


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
    
    arbitraryRotateSingle = q
    
  end function arbitraryRotateSingle
  
  TYPE(doubleVector) FUNCTION arbitraryRotateDouble(p, theta, r)
    TYPE(doubleVector) :: p,q    ! position vector
    REAL :: theta          ! angle in radians
    TYPE(doubleVector) :: r      ! the arbitrary axis
    REAL :: cosTheta, sinTheta


    costheta = cos(theta)
    sintheta = sin(theta)
    
    q%x = (costheta + (1.0_db - costheta) * r%x * r%x) * p%x
    q%x = q%x + ((1.0_db - costheta) * r%x * r%y - r%z * sintheta) * p%y
    q%x = q%x + ((1.0_db - costheta) * r%x * r%z + r%y * sintheta) * p%z
    
    q%y = ((1.0_db - costheta) * r%x * r%y + r%z * sintheta) * p%x
    q%y = q%y + (costheta + (1.0_db - costheta) * r%y * r%y) * p%y
    q%y = q%y + ((1.0_db - costheta) * r%y * r%z - r%x * sintheta) * p%z
    
    q%z = ((1.0_db - costheta) * r%x * r%z - r%y * sintheta) * p%x
    q%z = q%z + ((1.0_db - costheta) * r%y * r%z + r%x * sintheta) * p%y
    q%z = q%z + (costheta + (1.0_db - costheta) * r%z * r%z) * p%z
    
    arbitraryRotateDouble = q
    
  end function arbitraryRotateDouble
  
  TYPE(octalVector) FUNCTION arbitraryRotateOctal(p, theta, r)
    TYPE(octalVector) :: p,q    ! position vector
    REAL :: theta          ! angle in radians
    TYPE(octalVector) :: r      ! the arbitrary axis
    REAL :: cosTheta, sinTheta


    costheta = cos(theta)
    sintheta = sin(theta)
    
    q%x = (costheta + (1.0_oc - costheta) * r%x * r%x) * p%x
    q%x = q%x + ((1.0_oc - costheta) * r%x * r%y - r%z * sintheta) * p%y
    q%x = q%x + ((1.0_oc - costheta) * r%x * r%z + r%y * sintheta) * p%z
    
    q%y = ((1.0_oc - costheta) * r%x * r%y + r%z * sintheta) * p%x
    q%y = q%y + (costheta + (1.0_oc - costheta) * r%y * r%y) * p%y
    q%y = q%y + ((1.0_oc - costheta) * r%y * r%z - r%x * sintheta) * p%z
    
    q%z = ((1.0_oc - costheta) * r%x * r%z - r%y * sintheta) * p%x
    q%z = q%z + ((1.0_oc - costheta) * r%y * r%z + r%x * sintheta) * p%y
    q%z = q%z + (costheta + (1.0_oc - costheta) * r%z * r%z) * p%z
    
    arbitraryRotateOctal = q
    
  end function arbitraryRotateOctal


end module vector_mod

