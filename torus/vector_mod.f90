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
     real(kind=singleKind) :: x
     real(kind=singleKind) :: y
     real(kind=singleKind) :: z
  end type VECTOR

  type doubleVector
    real(kind=doubleKind) :: x
    real(kind=doubleKind) :: y
    real(kind=doubleKind) :: z
  end type doubleVector
  
  type octalVector
    real(kind=octalKind) :: x
    real(kind=octalKind) :: y
    real(kind=octalKind) :: z
  end type octalVector

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

  ! assignment

  interface assignment(=)
     module procedure singleToDoubleVector
     module procedure doubleToSingleVector
     module procedure singleToOctalVector
     module procedure octalToSingleVector
     module procedure doubleToOctalVector
     module procedure octalToDoubleVector
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
  
  interface modulus
    module procedure modulusSingle
    module procedure modulusDouble
    module procedure modulusOctal
  end interface

  interface normalize
    module procedure normalizeSingle
    module procedure normalizeDouble
    module procedure normalizeOctal
  end interface

  interface rotateX
    module procedure rotateXSingle
    module procedure rotateXDouble
    module procedure rotateXOctal
  end interface

  interface rotateY
    module procedure rotateYSingle
    module procedure rotateYDouble
    module procedure rotateYOctal
  end interface

  interface rotateZ
    module procedure rotateZSingle
    module procedure rotateZDouble
    module procedure rotateZOctal
  end interface

  interface intersectionLinePlane
    module procedure intersectionLinePlaneSingle
    module procedure intersectionLinePlaneDouble
    module procedure intersectionLinePlaneOctal
  end interface

  interface intersectionLineSphere
    module procedure intersectionLineSphereSingle
    module procedure intersectionLineSphereDouble
    module procedure intersectionLineSphereOctal
  end interface

  interface arbitraryRotate
    module procedure arbitraryRotateSingle
    module procedure arbitraryRotateDouble
    module procedure arbitraryRotateOctal
  end interface
  
  interface getPolar
    module procedure getPolarSingle
    module procedure getPolarDouble
    module procedure getPolarOctal
  end interface
  
  type(VECTOR), parameter :: xHat = VECTOR(1., 0., 0.)
  type(VECTOR), parameter :: yHat = VECTOR(0., 1., 0.)
  type(VECTOR), parameter :: zHat = VECTOR(0., 0., 1.)
  
  type(doubleVector),parameter :: xHatDouble=doubleVector(1.0_db,0.0_db,0.0_db)
  type(doubleVector),parameter :: yHatDouble=doubleVector(0.0_db,1.0_db,0.0_db)
  type(doubleVector),parameter :: zHatDouble=doubleVector(0.0_db,0.0_db,1.0_db)

  type(octalVector),parameter :: xHatOctal=octalVector(1.0_oc,0.0_oc,0.0_oc)
  type(octalVector),parameter :: yHatOctal=octalVector(0.0_oc,1.0_oc,0.0_oc)
  type(octalVector),parameter :: zHatOctal=octalVector(0.0_oc,0.0_oc,1.0_oc)

contains

  ! the dot product function

  real pure function dotProdSingle(a , b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    dotProdSingle = a%x*b%x + a%y*b%y + a%z*b%z

  end function dotProdSingle
  
  real(kind=doubleKind) pure function dotProdDouble(a , b)
    type(doubleVector), intent(in) :: a
    type(doubleVector), intent(in) :: b

    dotProdDouble = a%x*b%x + a%y*b%y + a%z*b%z

  end function dotProdDouble
  
  real(kind=octalKind) pure function dotProdOctal(a , b)
    type(octalVector), intent(in) :: a
    type(octalVector), intent(in) :: b

    dotProdOctal = a%x*b%x + a%y*b%y + a%z*b%z

  end function dotProdOctal


  ! the cross product function

  type(VECTOR) pure function crossProdSingle(a ,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    crossProdSingle%x =  (a%y*b%z - a%z*b%y)
    crossProdSingle%y = -(a%x*b%z - a%z*b%x)
    crossProdSingle%z =  (a%x*b%y - a%y*b%x)
  end function crossProdSingle
  
  type(doubleVector) pure function crossProdDouble(a,b)
    type(doubleVector), intent(in) :: a
    type(doubleVector), intent(in) :: b

    crossProdDouble%x =  (a%y*b%z - a%z*b%y)
    crossProdDouble%y = -(a%x*b%z - a%z*b%x)
    crossProdDouble%z =  (a%x*b%y - a%y*b%x)
  end function crossProdDouble

  type(octalVector) pure function crossProdOctal(a,b)
    type(octalVector), intent(in) :: a
    type(octalVector), intent(in) :: b

    crossProdOctal%x =  (a%y*b%z - a%z*b%y)
    crossProdOctal%y = -(a%x*b%z - a%z*b%x)
    crossProdOctal%z =  (a%x*b%y - a%y*b%x)
  end function crossProdOctal


  ! normalization subroutine - checks for zero vector

  subroutine normalizeSingle(a)
    type(VECTOR), intent(inout) :: a
    real :: m

    m = modulus(a)

    if (m == 0.) then
       write(*,'(a)') "! Attempt to normalize the zero vector"
       a = VECTOR(1.,0.,0.)
    else
       a%x = a%x / m
       a%y = a%y / m
       a%z = a%z / m
    endif

  end subroutine normalizeSingle
  
  subroutine normalizeDouble(a)
    type(doubleVector), intent(inout) :: a
    real(kind=doubleKind) :: m

    m = modulus(a)

    if (m == 0.0_db) then
       write(*,'(a)') "! Attempt to normalize the zero vector"
       m = SQRT(m - 1.0_db)
       a = doubleVector(1.0_db,0.0_db,0.0_db)
    endif

    a%x = a%x / m
    a%y = a%y / m
    a%z = a%z / m

  end subroutine normalizeDouble

  subroutine normalizeOctal(a)
    type(octalVector), intent(inout) :: a
    real(kind=octalKind) :: m

    m = modulus(a)

    if (m == 0.0_oc) then
       WRITE(*,'(a)') "! Attempt to normalize the zero vector"
       m = SQRT(m - 1.0_oc)
       a = octalVector(1.0_oc,0.0_oc,0.0_oc)
    endif

    a%x = a%x / m
    a%y = a%y / m
    a%z = a%z / m

  end subroutine normalizeOctal


  ! find the modulus of a vector

  real pure function modulusSingle(a)
    type(VECTOR), intent(in) :: a

    modulusSingle = a%x*a%x + a%y*a%y + a%z*a%z
    modulusSingle = sqrt(modulusSingle)

  end function modulusSingle
  
  real(kind=doubleKind) pure function modulusDouble(a)
    type(doubleVector), intent(in) :: a

    modulusDouble = a%x*a%x + a%y*a%y + a%z*a%z
    modulusDouble = SQRT(modulusDouble)

  end function modulusDouble
  
  real(kind=octalKind) pure function modulusOctal(a)
    type(octalVector), intent(in) :: a

    modulusOctal = a%x*a%x + a%y*a%y + a%z*a%z
    modulusOctal = SQRT(modulusOctal)

  end function modulusOctal

  ! multiply function

  type(VECTOR) pure function rmultSingle(a,b)
    real, intent(in) :: a
    type(VECTOR), intent(in) :: b

    rmultSingle%x = a * b%x
    rmultSingle%y = a * b%y
    rmultSingle%z = a * b%z

  end function rmultSingle

  type(doubleVector) pure function rmultDouble(a,b)
    real(kind=doubleKind), intent(in) :: a
    type(doubleVector), intent(in) :: b

    rmultDouble%x = a * b%x
    rmultDouble%y = a * b%y
    rmultDouble%z = a * b%z

  end function rmultDouble

  type(octalVector) pure function rmultOctal(a,b)
    real(kind=octalKind), intent(in) :: a
    type(octalVector), intent(in) :: b

    rmultOctal%x = a * b%x
    rmultOctal%y = a * b%y
    rmultOctal%z = a * b%z

  end function rmultOctal

  type(VECTOR) pure function rmultSingleReversed(b,a)
    real, intent(in) :: a
    type(VECTOR), intent(in) :: b

    rmultSingleReversed%x = a * b%x
    rmultSingleReversed%y = a * b%y
    rmultSingleReversed%z = a * b%z

  end function rmultSingleReversed

  type(doubleVector) pure function rmultDoubleReversed(b,a)
    real(kind=doubleKind), intent(in) :: a
    type(doubleVector), intent(in) :: b

    rmultDoubleReversed%x = a * b%x
    rmultDoubleReversed%y = a * b%y
    rmultDoubleReversed%z = a * b%z

  end function rmultDoubleReversed

  type(octalVector) pure function rmultOctalReversed(b,a)
    real(kind=octalKind), intent(in) :: a
    type(octalVector), intent(in) :: b

    rmultOctalReversed%x = a * b%x
    rmultOctalReversed%y = a * b%y
    rmultOctalReversed%z = a * b%z

  end function rmultOctalReversed

 
  type(VECTOR) pure function dmult(a,b)
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

  type(VECTOR) pure function divideVecSingle(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b

    divideVecSingle%x = a%x / b
    divideVecSingle%y = a%y / b
    divideVecSingle%z = a%z / b

  end function divideVecSingle
  
  type(doubleVector) pure function divideVecDouble(a,b)
    type(doubleVector), intent(in) :: a
    real(kind=doubleKind), intent(in) :: b

    divideVecDouble%x = a%x / b
    divideVecDouble%y = a%y / b
    divideVecDouble%z = a%z / b

  end function divideVecDouble
  
  type(octalVector) pure function divideVecOctal(a,b)
    type(octalVector), intent(in) :: a
    real(kind=octalKind), intent(in) :: b

    divideVecOctal%x = a%x / b
    divideVecOctal%y = a%y / b
    divideVecOctal%z = a%z / b

  end function divideVecOctal


  ! add two vectors
  
  type(VECTOR) pure function addSingle(a,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    addSingle%x = a%x + b%x
    addSingle%y = a%y + b%y
    addSingle%z = a%z + b%z

  end function addSingle
  
  type(doubleVector) pure function addDouble(a,b)
    type(doubleVector), intent(in) :: a
    type(doubleVector), intent(in) :: b

    addDouble%x = a%x + b%x
    addDouble%y = a%y + b%y
    addDouble%z = a%z + b%z

  end function addDouble
  
  type(octalVector) pure function addOctal(a,b)
    type(octalVector), intent(in) :: a
    type(octalVector), intent(in) :: b

    addOctal%x = a%x + b%x
    addOctal%y = a%y + b%y
    addOctal%z = a%z + b%z

  end function addOctal


  ! subtract two vectors

  type(VECTOR) pure function subtractSingle(a,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    subtractSingle%x = a%x - b%x
    subtractSingle%y = a%y - b%y
    subtractSingle%z = a%z - b%z

  end function subtractSingle
  
  type(doubleVector) pure function subtractDouble(a,b)
    type(doubleVector), intent(in) :: a
    type(doubleVector), intent(in) :: b

    subtractDouble%x = a%x - b%x
    subtractDouble%y = a%y - b%y
    subtractDouble%z = a%z - b%z

  end function subtractDouble

  type(octalVector) pure function subtractOctal(a,b)
    type(octalVector), intent(in) :: a
    type(octalVector), intent(in) :: b

    subtractOctal%x = a%x - b%x
    subtractOctal%y = a%y - b%y
    subtractOctal%z = a%z - b%z

  end function subtractOctal
  
  pure subroutine doubleToSingleVector(a,b)

    type(vector),       intent(out) :: a
    type(doubleVector), intent(in)  :: b

    a%x = b%x 
    a%y = b%y
    a%z = b%z
    
  end subroutine doubleToSingleVector 
  
  pure subroutine singleToDoubleVector(a,b)

    type(doubleVector), intent(out) :: a
    type(vector),       intent(in)  :: b

    a%x = b%x 
    a%y = b%y
    a%z = b%z
    
  end subroutine singleToDoubleVector 
  
  pure subroutine octalToSingleVector(a,b)

    type(vector),       intent(out) :: a
    type(octalVector), intent(in)   :: b

    a%x = b%x 
    a%y = b%y
    a%z = b%z
    
  end subroutine octalToSingleVector 
  
  pure subroutine singleToOctalVector(a,b)

    type(octalVector), intent(out)  :: a
    type(vector),       intent(in)  :: b

    a%x = b%x 
    a%y = b%y
    a%z = b%z
    
  end subroutine singleToOctalVector 

  pure subroutine octalToDoubleVector(a,b)

    type(doubleVector), intent(out) :: a
    type(octalVector), intent(in)   :: b

    a%x = b%x 
    a%y = b%y
    a%z = b%z
    
  end subroutine octalToDoubleVector 
  
  pure subroutine doubleToOctalVector(a,b)

    type(octalVector), intent(out)  :: a
    type(doubleVector), intent(in)  :: b

    a%x = b%x 
    a%y = b%y
    a%z = b%z
    
  end subroutine doubleToOctalVector 

  ! get polar form of a cartesian vector
  
  pure subroutine getPolarSingle(vec, r, theta, phi)

    implicit none
    type(vector), intent(in) :: vec
    real, intent(out) :: r, theta, phi
    real :: cosTheta 

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
  
  pure subroutine getPolarDouble(vec, r, theta, phi)

    implicit none
    type(doubleVector), intent(in) :: vec
    real(kind=doubleKind),intent(out) :: r, theta, phi
    real(kind=doubleKind) :: cosTheta 

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
  
  pure subroutine getPolarOctal(vec, r, theta, phi)

    implicit none
    type(octalVector), intent(in) :: vec
    real(kind=octalKind),intent(out) :: r, theta, phi
    real(kind=octalKind) :: cosTheta 

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

  type(VECTOR) pure function rotateZSingle(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateZSingle%x = cosb * a%x + sinb * a%y
    rotateZSingle%y =-sinb * a%x + cosb * a%y
    rotateZSingle%z = a%z

  end function rotateZSingle
  
  type(doubleVector) pure function rotateZDouble(a,b)
    type(doubleVector), intent(in) :: a
    real(kind=doubleKind), intent(in) :: b   ! angle in radians
    real(kind=doubleKind) :: cosb, sinb

    cosb = COS(b)
    sinb = Sin(b)

    rotateZDouble%x = cosb * a%x + sinb * a%y
    rotateZDouble%y =-sinb * a%x + cosb * a%y
    rotateZDouble%z = a%z

  end function rotateZDouble
  
  type(octalVector) pure function rotateZOctal(a,b)
    type(octalVector), intent(in) :: a
    real(kind=octalKind), intent(in) :: b   ! angle in radians
    real(kind=octalKind) :: cosb, sinb

    cosb = COS(b)
    sinb = Sin(b)

    rotateZOctal%x = cosb * a%x + sinb * a%y
    rotateZOctal%y =-sinb * a%x + cosb * a%y
    rotateZOctal%z = a%z

  end function rotateZOctal


  type(VECTOR) pure function rotateXSingle(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateXSingle%x = a%x 
    rotateXSingle%y = cosb * a%y + sinb * a%z
    rotateXSingle%z =-sinb * a%y + cosb * a%z

  end function rotateXSingle
  
  type(doubleVector) pure function rotateXDouble(a,b)
    type(doubleVector), intent(in) :: a
    real(kind=doubleKind), intent(in) :: b   ! angle in radians
    real(kind=doubleKind)  :: cosb, sinb

    cosb = COS(b)
    sinb = Sin(b)

    rotateXDouble%x = a%x 
    rotateXDouble%y = cosb * a%y + sinb * a%z
    rotateXDouble%z =-sinb * a%y + cosb * a%z

  end function rotateXDouble
  
  type(octalVector) pure function rotateXOctal(a,b)
    type(octalVector), intent(in) :: a
    real(kind=octalKind), intent(in) :: b   ! angle in radians
    real(kind=octalKind)  :: cosb, sinb

    cosb = COS(b)
    sinb = Sin(b)

    rotateXOctal%x = a%x 
    rotateXOctal%y = cosb * a%y + sinb * a%z
    rotateXOctal%z =-sinb * a%y + cosb * a%z

  end function rotateXOctal


  type(VECTOR) pure function rotateYSingle(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateYSingle%x = cosb * a%x + sinb * a%z
    rotateYSingle%y = a%y
    rotateYSingle%z =-sinb * a%x + cosb * a%z

  end function rotateYSingle
  
  type(doubleVector) pure function rotateYDouble(a,b)
    type(doubleVector), intent(in) :: a
    real(kind=doubleKind), intent(in) :: b   ! angle in radians
    real(kind=doubleKind) :: cosb, sinb

    cosb = COS(b)
    sinb = Sin(b)

    rotateYDouble%x = cosb * a%x + sinb * a%z
    rotateYDouble%y = a%y
    rotateYDouble%z =-sinb * a%x + cosb * a%z

  end function rotateYDouble
  
  type(octalVector) pure function rotateYOctal(a,b)
    type(octalVector), intent(in) :: a
    real(kind=octalKind), intent(in) :: b   ! angle in radians
    real(kind=octalKind) :: cosb, sinb

    cosb = COS(b)
    sinb = Sin(b)

    rotateYOctal%x = cosb * a%x + sinb * a%z
    rotateYOctal%y = a%y
    rotateYOctal%z =-sinb * a%x + cosb * a%z

  end function rotateYOctal


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

    implicit none

    type(VECTOR), intent(in) :: r0, rHat ! equation of line
    type(VECTOR), intent(in) :: nHat     ! the normal to the plane
    real, intent(in) :: d                ! minimum distance of plane from origin
    logical, intent(out) :: ok           ! is there in intersection?
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
  
  type (doubleVector) function intersectionLinePlaneDouble &
                                           (r0, rHat, nHat, d, ok)
    ! finds the intersection between a line and a plane

    implicit none

    type(doubleVector), intent(in) :: r0, rHat    ! equation of line
    type(doubleVector), intent(in) :: nHat        ! the normal to the plane
    real(kind=doubleKind) , intent(in) :: d
                         ! minimum distance of plane from origin
    logical,intent(out) :: ok               ! is there in intersection?
    real(kind=doubleKind) :: fac

    ok = .false.
    if ((nHat.dot.rHat) /= 0.0_db) then
       fac = ((d - (nHat.dot.r0))/(nHat.dot.rHat))
       if (fac > 0.0_db) then
          intersectionLinePlaneDouble = r0 + (fac * rHat)
          ok = .true.
       else
          intersectionLinePlaneDouble = doubleVector(0.0_db, 0.0_db, 0.0_db)
          ok = .false.
       endif
    else
       intersectionLinePlaneDouble = doubleVector(0.0_db ,0.0_db ,0.0_db)
       ok = .false.
    endif

  end function intersectionLinePlaneDouble
  
  type (octalVector) function intersectionLinePlaneOctal &
                                           (r0, rHat, nHat, d, ok)
    ! finds the intersection between a line and a plane

    implicit none

    type(octalVector), intent(in) :: r0, rHat    ! equation of line
    type(octalVector), intent(in) :: nHat        ! the normal to the plane
    real(kind=octalKind) , intent(in) :: d
                         ! minimum distance of plane from origin
    logical,intent(out) :: ok               ! is there in intersection?
    real(kind=octalKind) :: fac

    ok = .false.
    if ((nHat.dot.rHat) /= 0.0_oc) then
       fac = ((d - (nHat.dot.r0))/(nHat.dot.rHat))
       if (fac > 0.0_oc) then
          intersectionLinePlaneOctal = r0 + (fac * rHat)
          ok = .true.
       else
          intersectionLinePlaneOctal = octalVector(0.0_oc, 0.0_oc, 0.0_oc)
          ok = .false.
       endif
    else
       intersectionLinePlaneOctal = octalVector(0.0_oc ,0.0_oc ,0.0_oc)
       ok = .false.
    endif

  end function intersectionLinePlaneOctal

  
  pure subroutine intersectionLineSphereSingle(r0, rHat, length, s0, sR, found1,&
                              found2,intersectionDistance1,intersectionDistance2)
    ! finds the intersection(s) between a line segment and a sphere
    
    implicit none

    type(vector), intent(in)      :: r0, rHat ! equation of line
    real, intent(in)              :: length   ! length of line segment 
    type(vector), intent(in)      :: s0       ! centre of sphere
    real, intent(in)              :: sR       ! radius of sphere
    logical,intent(out)           :: found1   ! there is one intersection
    logical,intent(out)           :: found2   ! there is a second intersection
    real, intent(out)             :: intersectionDistance1 ! distances along line 
    real, intent(out)             :: intersectionDistance2 !   to intersections
    real(doubleKind)  :: a, b, c       ! quadratic formula variables
    real(doubleKind)  :: discriminant  ! quadratic formula variable
    real(doubleKind)  :: solution1     ! quadratic formula variable
    real(doubleKind)  :: solution2     ! quadratic formula variable
    
    type(doubleVector) :: r0Double, rHatDouble ! equation of line (double precision)
    type(doubleVector) :: s0Double             ! centre of sphere (double precision)
    real(doubleKind)   :: sRDouble             ! radius of sphere (double precision)
    
    r0Double = doubleVector(r0%x,r0%y,r0%z)
    rHatDouble = doubleVector(rHat%x,rHat%y,rHat%z)
    s0Double = doubleVector(s0%x,s0%y,s0%z)
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
    
  end subroutine intersectionLineSphereSingle
  
  pure subroutine intersectionLineSphereDouble(r0, rHat, length, s0, sR, found1,&
                              found2,intersectionDistance1,intersectionDistance2)
    ! finds the intersection between a line segment and a sphere
    
    implicit none

    type(doubleVector), intent(in) :: r0, rHat ! equation of line
    real(doubleKind), intent(in)   :: length   ! length of line segment 
    type(doubleVector), intent(in) :: s0       ! centre of sphere
    real(doubleKind), intent(in)   :: sR       ! radius of sphere
    logical,intent(out)            :: found1   ! there is one intersection
    logical,intent(out)            :: found2   ! there is a second intersection
    real(doubleKind), intent(out)  :: intersectionDistance1 ! distances along line 
    real(doubleKind), intent(out)  :: intersectionDistance2 !   to intersections
    real(doubleKind)  :: a, b, c       ! quadratic formula variables
    real(doubleKind)  :: discriminant  ! quadratic formula variable
    real(doubleKind)  :: solution1     ! quadratic formula variable
    real(doubleKind)  :: solution2     ! quadratic formula variable
    
    a = rHat%x**2.0_db + rHat%y**2.0_db + rHat%z**2.0_db
    b = 2.0_db * ( rHat%x * (r0%x - s0%x) + rHat%y * (r0%y - s0%y) + rHat%z * (r0%z - s0%z)) 
    c = (r0%x - s0%x)**2.0_db + (r0%y - s0%y)**2.0_db + (r0%z - s0%z)**2.0_db - sR**2.0_db 
 
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
    
  end subroutine intersectionLineSphereDouble

  pure subroutine intersectionLineSphereOctal(r0, rHat, length, s0, sR, found1, &
                               found2,intersectionDistance1,intersectionDistance2)
    ! finds the intersection between a line segment and a sphere
    
    implicit none

    type(octalVector), intent(in) :: r0, rHat ! equation of line
    real(octalKind), intent(in)   :: length   ! length of line segment 
    type(octalVector), intent(in) :: s0       ! centre of sphere
    real(octalKind), intent(in)   :: sR       ! radius of sphere
    logical,intent(out)           :: found1   ! there is one intersection
    logical,intent(out)           :: found2   ! there is a second intersection
    real(octalKind), intent(out)  :: intersectionDistance1 ! distances along line 
    real(octalKind), intent(out)  :: intersectionDistance2 !   to intersections
    real(doubleKind)  :: a, b, c       ! quadratic formula variables
    real(doubleKind)  :: discriminant  ! quadratic formula variable
    real(doubleKind)  :: solution1     ! quadratic formula variable
    real(doubleKind)  :: solution2     ! quadratic formula variable
    
    type(doubleVector) :: r0Double, rHatDouble ! equation of line (double precision)
    type(doubleVector) :: s0Double             ! centre of sphere (double precision)
    real(doubleKind)   :: sRDouble             ! radius of sphere (double precision)
    
    r0Double = doubleVector(r0%x,r0%y,r0%z)
    rHatDouble = doubleVector(rHat%x,rHat%y,rHat%z)
    s0Double = doubleVector(s0%x,s0%y,s0%z)
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
    
  end subroutine intersectionLineSphereOctal


  type(VECTOR) pure function arbitraryRotateSingle(p, theta, r)
    type(VECTOR),intent(in) :: p      ! position vector
    real,intent(in)         :: theta  ! angle in radians
    type(VECTOR),intent(in) :: r      ! the arbitrary axis
    
    real :: cosTheta, sinTheta
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
    
    arbitraryRotateSingle = q
    
  end function arbitraryRotateSingle
  
  type(doubleVector) pure function arbitraryRotateDouble(p, theta, r)
    type(doubleVector),intent(in) :: p     ! position vector
    real,intent(in)               :: theta ! angle in radians
    type(doubleVector),intent(in) :: r     ! the arbitrary axis
    
    real :: cosTheta, sinTheta
    type(doubleVector)            :: q     ! position vector


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
  
  type(octalVector) pure function arbitraryRotateOctal(p, theta, r)
    type(octalVector),intent(in) :: p     ! position vector
    real,intent(in)              :: theta ! angle in radians
    type(octalVector),intent(in) :: r     ! the arbitrary axis
    
    real :: cosTheta, sinTheta
    type(octalVector) :: q    ! position vector


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

  type(VECTOR) function fromPhotosphereVector(rVec)
    real :: ang, z
    type(VECTOR) :: norm, zAxis, v, rVec, n
    
    zAxis = VECTOR(0.,0.,1.)

    norm = rVec
    call normalize(norm)

    call random_number(z)
    z = sqrt(z)

    ang = acos(z)
    n = norm .cross. zAxis
    call normalize(n)
    v = norm

    v = arbitraryRotateSingle(v, ang, n)
    call random_number(ang)
    ang = ang * twoPi
    v = arbitraryRotateSingle(v, ang, norm)

    fromPhotosphereVector = v
  end function fromPhotosphereVector



  
end module vector_mod

