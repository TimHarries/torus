
! this module sets up the stokes vector type and the corresponding
! phase matrix type. It includes subroutines to create and
! apply phase matrices

! written by tjh

! v1.0 on 16/08/99


module phasematrix_mod

  use constants_mod
  use vector_mod
  use messages_mod 
  use random_mod
  implicit none

  public

  ! the phase matrix is 4x4

  type PHASEMATRIX
     real(double) :: element(4,4)
     real(double) :: gfac
  end type PHASEMATRIX


  ! the stokes vector includes each of the stokes intensities

  type STOKESVECTOR
     real(double) :: i
     real(double) :: q
     real(double) :: u
     real(double) :: v
  end type STOKESVECTOR

  ! + signifies adding two stokes vectors

  interface operator(+)
     module procedure addStokes
  end interface

  ! - signifies subtracting two stokes vectors

  interface operator(-)
     module procedure subStokes
  end interface

  ! * multiplies a stokes vector by a constant

  interface operator(*)
     module procedure multStokes
     module procedure multStokes_dble
  end interface

  interface operator(+)
     module procedure addPhase
  end interface

  interface operator(-)
     module procedure subPhase
  end interface


  interface operator(*)
     module procedure multPhase
  end interface




contains

  ! this function applies a phase matrix to a stokes vector
  
  


  pure function apply(a, b) result(c)

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

  end function apply

  function applyMean(a, f,  n, b) result(c)

    integer, intent(in) :: n
    type(PHASEMATRIX), intent(in) :: a(n)
    real(double), intent(in) :: f(:)
    type(STOKESVECTOR), intent(in) :: b
    type(STOKESVECTOR) :: c
    integer :: i
    real(double), allocatable :: normf(:)
    type(PHASEMATRIX)  :: mean
    allocate(normf(1:n))
    mean%element = 0.
    normf(1:n) = f(1:n)/SUM(f)
    if (ANY(normf > 1.d0)) then
       write(*,*) "f ",f(1:n)
       write(*,*) "normf ",normf(1:n)
       write(*,*) "sum f ",sum(f)
    endif
    do i = 1 , n
       mean%element = real(mean%element + normf(i) * a(i)%element)
    enddo

    ! perform the matrix multiplication
    
    c%i = mean%element(1,1) * b%i + &
         mean%element(1,2) * b%q + &
         mean%element(1,3) * b%u + &
         mean%element(1,4) * b%v


    c%q = mean%element(2,1) * b%i + &
         mean%element(2,2) * b%q + &
         mean%element(2,3) * b%u + &
         mean%element(2,4) * b%v


    c%u = mean%element(3,1) * b%i + &
         mean%element(3,2) * b%q + &
         mean%element(3,3) * b%u + &
         mean%element(3,4) * b%v


    c%v = mean%element(4,1) * b%i + &
         mean%element(4,2) * b%q + &
         mean%element(4,3) * b%u + &
         mean%element(4,4) * b%v 
    deallocate(normf)
  end function applyMean


  ! this sets up a Rayleigh phase matrix

  pure function fillRayleigh(costheta) result(b)

    type(PHASEMATRIX) :: b
    real,intent(in) :: costheta
    real :: cos2t

    cos2t = costheta*costheta

    b%element = 0.
    b%element(1,1) = 0.75*(1.0+cos2t)
    b%element(1,2) = 0.75*(cos2t-1.0)
    b%element(2,1) = 0.75*(cos2t-1.0)
    b%element(2,2) = 0.75*(1.0+cos2t)
    b%element(3,3) = 0.75*2.*costheta
    b%element(4,4) = 0.75*2.*costheta

  end function fillRayleigh


  pure function fillHenyey(costheta, g) result(b)

    type(PHASEMATRIX) :: b
    real(double),intent(in) :: costheta, g

    b%element = 0.
    b%element(1,1) = oneOnFourPi * (1.d0-g**2) / (1.d0 + g**2 - 2.d0*g*costheta)**1.5
    b%gFac = g

  end function fillHenyey

  pure function fillIsotropic() result(b)

    type(PHASEMATRIX) :: b

    b%element = 0.
    b%element(1,1) = 1.

  end function fillIsotropic

  ! this writes out a 4x4 phase matrix - used for debugging

  subroutine writePhaseMatrix(a)

    type(PHASEMATRIX), intent(in) :: a
    integer :: i

    write(*,*) " "
    write(*,'(a)') "----------------------------"
    do i = 1, 4
       write(*,'(f6.3,1x,f6.3,1x,f6.3,1x,f6.3)') a%element(i,1),  &
            a%element(i,2),  a%element(i,3),  a%element(i,4)
    enddo
    write(*,'(a)') "----------------------------"
  end subroutine writePhaseMatrix


  type(PHASEMATRIX) pure function multPhase(a , b)
    real, intent(in) :: a
    type(PHASEMATRIX), intent(in) :: b

    multPhase%element(1:4,1:4) = a * b%element(1:4,1:4)
  end function multPhase

  type(PHASEMATRIX) pure function addPhase(a , b)
    type(PHASEMATRIX), intent(in) :: a
    type(PHASEMATRIX), intent(in) :: b

    addPhase%element(1:4,1:4) = a%element(1:4,1:4) + b%element(1:4,1:4)
  end function addPhase

  type(PHASEMATRIX) pure function subPhase(a , b)
    type(PHASEMATRIX), intent(in) :: a
    type(PHASEMATRIX), intent(in) :: b

    subPhase%element(1:4,1:4) = a%element(1:4,1:4) - b%element(1:4,1:4)
  end function subPhase


  ! function to add to stokes vectors

  type(STOKESVECTOR) pure function addStokes(a , b)
    type(STOKESVECTOR), intent(in) :: a
    type(STOKESVECTOR), intent(in) :: b

    addStokes%i = a%i + b%i
    addStokes%q = a%q + b%q
    addStokes%u = a%u + b%u
    addStokes%v = a%v + b%v

  end function addStokes

  ! function to subtract to stokes vectors

  type(STOKESVECTOR) pure function subStokes(a , b)
    type(STOKESVECTOR), intent(in) :: a
    type(STOKESVECTOR), intent(in) :: b

    subStokes%i = a%i - b%i
    subStokes%q = a%q - b%q
    subStokes%u = a%u - b%u
    subStokes%v = a%v - b%v

  end function subStokes

  ! multiply a stokes vector by a constant

  type(STOKESVECTOR) pure function multStokes(a , b)
    type(STOKESVECTOR), intent(in) :: a
    real, intent(in) :: b

    multStokes%i = a%i * b
    multStokes%q = a%q * b
    multStokes%u = a%u * b
    multStokes%v = a%v * b

  end function multStokes


  type(STOKESVECTOR) pure function multStokes_dble(a , b)
    type(STOKESVECTOR), intent(in) :: a
    real(double), intent(in) :: b

    multStokes_dble%i = a%i * b
    multStokes_dble%q = a%q * b
    multStokes_dble%u = a%u * b
    multStokes_dble%v = a%v * b

  end function multStokes_dble

  subroutine testMiePhase(wavelength, lamArray, nLambda, miePhase, nDustType, nMuMie, dustTypeFraction)
    real :: wavelength
    real :: lamArray(:)
    integer :: nLambda, nDustType, nMuMie
    real(double) :: dustTypeFraction(:), meanscat
    type(PHASEMATRIX) :: miePhase(nDustType, nLambda, nMuMie)
    type(VECTOR) :: oldDirection, newDirection, tot
    integer :: i, ntest

!    miePhase(:,:,:)%element(1,1) = 1.
    ntest = 100000
    tot=vector(0. ,0., 0)
    meanscat = 0
    do i = 1, nTest
       oldDirection = VECTOR(1., 0., 0.)
!       newDirection =  newDirectionMie(oldDirection, wavelength, lamArray, &
!            nLambda, miePhase, nDustType, nMuMie, dustTypeFraction)
!       write(*,*) "Scattering angle: ",acos(oldDirection.dot.newDirection)*180./pi
       meanscat=meanscat+acos(oldDirection.dot.newDirection)*180./pi
       tot=tot+newDirection
    enddo
    tot = tot / dble(nTest)
    meanscat = meanscat / real(ntest)
    write(*,*) "Final vector: ",tot
    write(*,*) "mean scattering angle ",meanscat
  end subroutine testMiePhase


  type(VECTOR) function newDirectionMie(grid, currentOctal, currentSubcell, &
       oldDirection, wavelength, &
       lamArray, nLambda, miePhase, nDustType, nMuMie, dustTypeFraction, weight)
    use amr_mod, only : returnKappa
    use gridtype_mod
    use inputs_mod, only : inputGfac, henyeyGreensteinPhaseFunction
    use utils_mod, only: locate
    type(GRIDTYPE) ::  grid
    type(OCTAL), pointer :: currentOctal
    integer :: currentSubcell
    type(VECTOR), intent(in) :: oldDirection
    real, intent(in) :: wavelength
    real(double) :: dustTypeFraction(:)
    integer :: nDustType
    real, optional :: weight
    integer, intent(in) :: nMuMie, nLambda
    integer :: i, j, k, m, ilam
    real(double) :: theta, phi
    real(double) :: r
    real, intent(in) :: lamArray(:)
    real(double) :: normfac(100)
    real(double) :: allSca(10)
    type(VECTOR) :: tVec, perpVec, newVec
    type(PHASEMATRIX) :: miePhase(:,:,:)
    integer :: nMuTemp, refineAT
    logical :: refine
    real(double), allocatable :: prob(:)
    real(double), allocatable :: cosArray(:)

    call locate(lamArray, nLambda, wavelength, ilam)
    if (ilam < 1) ilam = 1
    if (ilam > nLambda) ilam = nLambda

    allocate(prob(1:nMuMie))
    allocate(cosArray(1:nMuMie))
    do i = 1, nMuMie
       cosArray(i) = -1.d0 + 2.d0*dble(i-1)/dble(nMuMie-1)
    enddo
       
    call returnKappa(grid, currentOctal, currentSubcell, ilambda=ilam, allSca=allSca)
       
    k = randomIndex(allSca(1:nDustType), nDustType)
    
    prob = 0.d0
    do i = 2, nMuMie
       prob(i) = prob(i) + miePhase(k,ilam,i)%element(1,1)*(cosArray(i)-cosArray(i-1))
    enddo
    do i = 2, nMuMie
       prob(i) = prob(i) + prob(i-1)
    enddo
    prob(1:nMuMie) = prob(1:nMuMie)/prob(nMuMie)
 

    call randomNumberGenerator(getDouble=r)
    call locate(prob(1:nMuMie), nMuMie, r, j)
    theta = cosArray(j) + &
         (cosArray(j+1)-cosArray(j))*(r - prob(j))/(prob(j+1)-prob(j))

    if (present(weight)) then
       weight = real(miePhase(k, iLam, j)%element(1,1) + & 
            (miePhase(k, iLam, j+1)%element(1,1) - &
             miePhase(k, iLam, j)%element(1,1) ) * &
            (r - prob(j))/(prob(j+1)-prob(j)))
       if (weight <= 0.d0) then
          write(*,*) "weight ", weight, " ilam ", ilam, " j ",j
       endif
    endif


    theta = acos(max(-1.0_db,min(1.0_db,theta)))
    call randomNumberGenerator(getDouble=r)
    phi = twoPi * r

    tVec = oldDirection

    perpVec%x =  tVec%y
    perpVec%y = -tVec%x
    perpVec%z = 1.d-20
    call normalize(perpVec)


    newVec = arbitraryRotate(oldDirection, theta, perpVec)
    newVec = arbitraryRotate(newVec, phi, oldDirection)
    call normalize(newvec)

    newDirectionMie = newVec
666 continue
  end function newDirectionMie


  subroutine fixMiePhase(miePhase, nDustType, nLambda, nMuMie)
    implicit none
    integer :: nDustType, nLambda, nMuMie
    type(PHASEMATRIX) :: miePhase(:,:,:)
    integer :: i, j
    logical, save :: firstTime = .true.

!$OMP THREADPRIVATE(firstTime)

    do i = 1, nDustType
       do j = 1, nLambda
          if ((miePhase(i,j,nMuMie)%element(1,1)/max(1.e-30,miePhase(i,j,nMuMie-1)%element(1,1))) > 10.d0) then
             miePhase(i,j,nMuMie)%element  = real(10.d0 * miePhase(i,j,nMuMie-1)%element)
             if (writeoutput.and.firstTime) then
                call writeInfo("Undersampeld miephase fixed (near 180)",TRIVIAL)
                firstTime = .false.
             endif
          endif

          if ((miePhase(i,j,1)%element(1,1)/max(1.e-30,miePhase(i,j,2)%element(1,1))) > 10.d0) then
             miePhase(i,j,1)%element  = real(10.d0 * miePhase(i,j,2)%element)
             if (writeoutput.and.firstTime) then
                call writeInfo("Undersampeld miephase fixed (near 0)",TRIVIAL)
                firstTime = .false.
             endif
          endif
       enddo
    enddo
  end subroutine fixMiePhase

end module phasematrix_mod


