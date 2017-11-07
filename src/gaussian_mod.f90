module gaussian_mod

  use vector_mod
  use kind_mod

  implicit none

  public :: findFactor, GAUSSIAN
  private :: eval

  ! The definition of the vector type

  type GAUSSIAN
     real(double) :: sigma
     type(VECTOR) :: centre
     real(double) :: amplitude
  end type GAUSSIAN

contains

  subroutine findFactor(fac, position, gArray, ng)
    real :: fac
    integer :: ng, i
    type(VECTOR) :: position
    type(GAUSSIAN) :: gArray(:)

    fac = 0.
    do i = 1, ng
       fac = fac + eval(gArray(i),position)
    end do
  end subroutine findFactor

  real function eval(thisgaussian, position)
    type(GAUSSIAN) :: thisGaussian
    type(VECTOR) :: position
    real(oct) :: x

    x = modulus(thisGaussian%centre - position)
    eval = real(thisGaussian%amplitude * exp(-x**2 / (2.*thisGaussian%sigma**2)))
  end function eval

end module gaussian_mod
