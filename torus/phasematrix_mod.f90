!
! this module sets up the stokes vector type and the corresponding
! phase matrix type. It includes subroutines to create and
! apply phase matrices

! written by tjh

! v1.0 on 16/08/99


module phasematrix_mod

  use utils_mod
  implicit none

  public

  ! the phase matrix is 4x4

  type PHASEMATRIX
     real :: element(4,4)
  end type PHASEMATRIX


  ! the stokes vector includes each of the stokes intensities

  type STOKESVECTOR
     real :: i
     real :: q
     real :: u
     real :: v
  end type STOKESVECTOR

  ! + signifies adding two stokes vectors

  interface operator(+)
     module procedure addStokes
  end interface

  ! * multiplies a stokes vector by a constant

  interface operator(*)
     module procedure multStokes
  end interface

contains

  ! this function applies a phase matrix to a stokes vector
  

  type(STOKESVECTOR) function apply(a, b)

    type(PHASEMATRIX), intent(in) :: a
    type(STOKESVECTOR), intent(in) :: b
    type(STOKESVECTOR) :: c

    ! perform the matrix multiplication
    
    c%i = a%element(1,1) * b%i + &
         a%element(1,2) * b%q + &
         a%element(1,3) * b%u + &
         a%element(1,4) * b%v


    c%q = a%element(2,1) * b%i + &
         a%element(2,2) * b%q + &
         a%element(2,3) * b%u + &
         a%element(2,4) * b%v


    c%u = a%element(3,1) * b%i + &
         a%element(3,2) * b%q + &
         a%element(3,3) * b%u + &
         a%element(3,4) * b%v


    c%v = a%element(4,1) * b%i + &
         a%element(4,2) * b%q + &
         a%element(4,3) * b%u + &
         a%element(4,4) * b%v 

    apply = c


  end function apply


  ! this sets up a Rayleigh phase matrix

  type(PHASEMATRIX) function fillRayleigh(costheta)


    type(PHASEMATRIX) :: b
    real :: costheta
    real :: cos2t

    cos2t = costheta*costheta

    b%element = 0.
    b%element(1,1) = 0.75*(1.d0+cos2t)
    b%element(1,2) = 0.75*(cos2t-1.d0)
    b%element(2,1) = 0.75*(cos2t-1.d0)
    b%element(2,2) = 0.75*(1.d0+cos2t)
    b%element(3,3) = 0.75*2.*costheta
    b%element(4,4) = 0.75*2.*costheta

    fillRayleigh = b

  end function fillRayleigh

  ! this writes out a 4x4 phase matrix - used for debugging

  subroutine writePhaseMatrix(a)

    type(PHASEMATRIX) :: a
    integer :: i

    write(*,*) " "
    write(*,'(a)') "----------------------------"
    do i = 1, 4
       write(*,'(f6.3,1x,f6.3,1x,f6.3,1x,f6.3)') a%element(i,1),  &
            a%element(i,2),  a%element(i,3),  a%element(i,4)
    enddo
    write(*,'(a)') "----------------------------"
  end subroutine writePhaseMatrix

  ! function to add to stokes vectors

  type(STOKESVECTOR) function addStokes(a , b)
    type(STOKESVECTOR), intent(in) :: a
    type(STOKESVECTOR), intent(in) :: b

    addStokes%i = a%i + b%i
    addStokes%q = a%q + b%q
    addStokes%u = a%u + b%u
    addStokes%v = a%v + b%v

  end function addStokes

  ! multiply a stokes vector by a constant

  type(STOKESVECTOR) function multStokes(a , b)
    type(STOKESVECTOR), intent(in) :: a
    real, intent(in) :: b

    multStokes%i = a%i * b
    multStokes%q = a%q * b
    multStokes%u = a%u * b
    multStokes%v = a%v * b

  end function multStokes

  type(VECTOR) function newDirectionMie(oldDirection, wavelength, lamArray, nLambda, miePhase, nMuMie)
    type(VECTOR) :: oldDirection
    real :: wavelength
    integer :: nMuMie, nLambda, i, j
    real :: costheta, theta, r, phi
    real :: lamArray(*)
    real, allocatable :: prob(:)
    real, allocatable :: cosArray(:)
    type(VECTOR) :: tVec, perpVec, newVec
    
    type(PHASEMATRIX) :: miePhase(nMuMie, nLambda)
    

    allocate(prob(1:nMuMie))
    allocate(cosArray(1:nMuMie))
    call locate(lamArray, nLambda, wavelength, j)
    if (j < 1) j = 1
    if (j > nLambda) j = nLambda

    prob = 0.
    do i = 1, nMuMie
       cosArray(i) = -1. + 2.*real(i-1)/real(nMuMie-1)
       prob(i) = miePhase(j,i)%element(1,1)
    enddo
    do i = 2, nMuMie
       prob(i) = prob(i) + prob(i-1)
    enddo
    prob(1:nMuMie) = prob(1:nMuMie)/prob(nMuMie)
    call random_number(r)
    call locate(prob, nMuMie, r, j)
    cosTheta = cosArray(j) + &
         (cosArray(j+1)-cosArray(j))*(r - prob(j))/(prob(j+1)-prob(j))
    call random_number(r)
    phi = twoPi * r

    tVec = oldDirection

    perpVec%x =  tVec%y
    perpVec%y = -tVec%x
    perpVec%z = 0.
    call normalize(perpVec)

    theta = acos(min(1.,max(-1.,costheta)))
    newVec = arbitraryRotate(oldDirection, theta, perpVec)
    newVec = arbitraryRotate(newVec, phi, oldDirection)

    Deallocate(cosArray)
    deallocate(prob)

    newDirectionMie = newVec
    

  end function newDirectionMie

end module phasematrix_mod


