!
! This is a module for vector maths. This includes the vector
! type definitions, as well as simple operations such as add,
! substract, multiple and divide. The two vector functions
! dot product and cross product are also defined.
!

! written by tjh


! v1.0 on 13/08/99

module vector_mod


  use kind_mod
  use constants_mod

  implicit none

  public

  ! The definition of the vector type

  type VECTOR
     real :: x
     real :: y
     real :: z
  end type VECTOR

  ! Define multiply

  interface operator(*)
     module procedure rmult, dmult
  end interface

  ! and divide

  interface operator(/)
     module procedure divideVec
  end interface

  ! add

  interface operator(+)
     module procedure add
  end interface

  ! subtract

  interface operator(-)
     module procedure subtract
  end interface

  ! dot product

  interface operator(.dot.)
     module procedure dotProd
  end interface

  ! cross product

  interface operator(.cross.)
     module procedure crossProd
  end interface

contains

  ! the dot product function

  real function dotProd(a , b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    dotProd = a%x*b%x + a%y*b%y + a%z*b%z

  end function dotProd

  ! the cross product function

  type(VECTOR) function crossProd(a ,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    crossProd%x =  (a%y*b%z - a%z*b%y)
    crossProd%y = -(a%x*b%z - a%z*b%x)
    crossProd%z =  (a%x*b%y - a%y*b%x)
  end function crossProd

  ! normalization subroutine - checks for zero vector

  subroutine normalize(a)
    type(VECTOR) :: a
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

  end subroutine normalize

  ! find the modulus of a vector

  real function modulus(a)
    type(VECTOR) :: a

    modulus = a%x*a%x + a%y*a%y + a%z*a%z
    modulus = sqrt(modulus)

  end function modulus

  ! multiply function

  type(VECTOR) function rmult(a,b)
    real, intent(in) :: a
    type(VECTOR), intent(in) :: b

    rmult%x = a * b%x
    rmult%y = a * b%y
    rmult%z = a * b%z

  end function rmult

  ! multiply function

  type(VECTOR) function dmult(a,b)
    real(kind=doubleKind), intent(in) :: a
    type(VECTOR), intent(in) :: b

    dmult%x = a * b%x
    dmult%y = a * b%y
    dmult%z = a * b%z

  end function dmult


  ! divide vector by a scalar

  type(VECTOR) function divideVec(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b

    divideVec%x = a%x / b
    divideVec%y = a%y / b
    divideVec%z = a%z / b

  end function divideVec

  ! add two vectors
  
  type(VECTOR) function add(a,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    add%x = a%x + b%x
    add%y = a%y + b%y
    add%z = a%z + b%z

  end function add

  ! subtract two vectors

  type(VECTOR) function subtract(a,b)
    type(VECTOR), intent(in) :: a
    type(VECTOR), intent(in) :: b

    subtract%x = a%x - b%x
    subtract%y = a%y - b%y
    subtract%z = a%z - b%z

  end function subtract

  ! get polar form of a cartesian vector
  
  subroutine getPolar(vec, r, theta, phi)

    implicit none
    type(VECTOR) :: vec
    real :: r, theta, phi, cosTheta 

    r = modulus(vec)
    if ((vec%y == 0.) .and. (vec%x == 0)) then
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
  end subroutine getPolar

  ! rotate a vector "a" about the z-axis by angle b

  type(VECTOR) function rotateZ(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateZ%x = cosb * a%x + sinb * a%y
    rotateZ%y =-sinb * a%x + cosb * a%y
    rotateZ%z = a%z

  end function rotateZ

  type(VECTOR) function rotateX(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateX%x = a%x 
    rotateX%y = cosb * a%y + sinb * a%z
    rotateX%z =-sinb * a%y + cosb * a%z

  end function rotateX

  type(VECTOR) function rotateY(a,b)
    type(VECTOR), intent(in) :: a
    real, intent(in) :: b   ! angle in radians
    real :: cosb, sinb

    cosb = cos(b)
    sinb = sin(b)

    rotateY%x = cosb * a%x + sinb * a%z
    rotateY%y = a%y
    rotateY%z =-sinb * a%x + cosb * a%z

  end function rotateY

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

  type (VECTOR) function intersectionLinePlane(r0, rHat, nHat, d, ok)

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
          intersectionLinePlane = r0 + (fac  * rHat )
          ok = .true.
       else
          intersectionLinePlane = VECTOR(0.,0.,0.)
          ok = .false.
       endif
    else
       intersectionLinePlane = VECTOR(0.,0.,0.)
       ok = .false.
    endif

  end function intersectionLinePlane

  type(VECTOR) function arbitraryRotate(p, theta, r)
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
    
    arbitraryRotate = q
    
  end function arbitraryRotate

end module vector_mod

