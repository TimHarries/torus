
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

  type(PHASEMATRIX) function fillIsotropic(costheta)


    type(PHASEMATRIX) :: b
    real :: costheta
    real :: cos2t

    cos2t = costheta*costheta

    b%element = 0.
    b%element(1,1) = 1.
    b%element(1,2) = 0.
    b%element(2,1) = 0.
    b%element(2,2) = 0.
    b%element(3,3) = 0.
    b%element(4,4) = 0.

    fillIsotropic = b

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

  ! function to add to stokes vectors

  type(STOKESVECTOR) pure function addStokes(a , b)
    type(STOKESVECTOR), intent(in) :: a
    type(STOKESVECTOR), intent(in) :: b

    addStokes%i = a%i + b%i
    addStokes%q = a%q + b%q
    addStokes%u = a%u + b%u
    addStokes%v = a%v + b%v

  end function addStokes

  ! multiply a stokes vector by a constant

  type(STOKESVECTOR) pure function multStokes(a , b)
    type(STOKESVECTOR), intent(in) :: a
    real, intent(in) :: b

    multStokes%i = a%i * b
    multStokes%q = a%q * b
    multStokes%u = a%u * b
    multStokes%v = a%v * b

  end function multStokes

  type(VECTOR) function newDirectionMie(oldDirection, wavelength, lamArray, nLambda, miePhase, nMuMie)
    type(VECTOR), intent(in) :: oldDirection
    real, intent(in) :: wavelength
    integer, intent(in) :: nMuMie, nLambda
    integer :: i, j
    real :: costheta, theta, r, phi
    real, intent(in) :: lamArray(:)
    real, allocatable :: prob(:)
    real, allocatable :: cosArray(:)
    type(VECTOR) :: tVec, perpVec, newVec
    
    type(PHASEMATRIX) :: miePhase(nLambda, nMuMie)
    

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

    call normalize(newvec)
    newDirectionMie = newVec
    

  end function newDirectionMie


subroutine writeSpectrum(outFile,  nLambda, xArray, yArray,  errorArray, nOuterLoop, &
     normalizeSpectrum, useNdf, sed, objectDistance, jansky, SI, velocitySpace, lamLine)

  implicit none
  integer, intent(in) :: nLambda
  character(len=*), intent(in) :: outFile
  real, intent(in) :: xArray(nLambda)
  logical, intent(in) :: useNdf
  logical, intent(in) :: jansky
  integer, intent(in) :: nOuterLoop
  logical, intent(in) :: normalizeSpectrum, sed, velocitySpace
  real(double), intent(in) :: objectDistance
  real(double) :: area
  real, intent(in) :: lamLine
  type(STOKESVECTOR), intent(in) :: yArray(nLambda), errorArray(nOuterloop,nLambda)
  type(STOKESVECTOR), allocatable :: ytmpArray(:), tmpErrorArray(:,:)
  real, allocatable :: meanQ(:), meanU(:), sigQ(:), sigU(:)
  real(double), allocatable :: stokes_i(:), stokes_q(:), stokes_qv(:)
  real(double), allocatable :: stokes_u(:), stokes_uv(:), dlam(:)
  real, allocatable, dimension(:) :: tmpXarray
  real :: tot
  real :: x
  integer :: i
  logical :: SI

  allocate(ytmpArray(1:nLambda))
  allocate(tmpXarray(1:nLambda))
  allocate(tmpErrorArray(nOuterloop,nLambda))

  allocate(meanQ(1:nLambda))
  allocate(meanU(1:nLambda))
  allocate(sigQ(1:nLambda))
  allocate(sigU(1:nLambda))

  allocate(dlam(1:nLambda))
  allocate(stokes_i(1:nLambda))
  allocate(stokes_q(1:nLambda))
  allocate(stokes_qv(1:nLambda))
  allocate(stokes_u(1:nLambda))
  allocate(stokes_uv(1:nLambda))

!  x = SUM(yArray(1:min(10,nLambda))%i)/real(min(10,nLambda))

!  x = 1./x

!  x = 
  write(*,*) "starting to write spectrum"

  if (normalizeSpectrum) then
     if (yArray(1)%i /= 0.) then
        x = 1.d0/yArray(1)%i
!        x = 1.d0/yArray(nLambda)%i
     else
       x  = 1.d0
     endif
  else
     x = 1.d0
  endif

  write(*,*) "scaling by ",x
!  
  do i = 1, nLambda
     ytmpArray(i) = yArray(i) * x
  enddo


  where(errorArray%i /= 0.)
     tmpErrorArray%q = errorArray%q / errorArray%i
     tmpErrorArray%u = errorArray%u / errorArray%i
  elsewhere 
     tmpErrorArray%q = errorArray%q 
     tmpErrorArray%u = errorArray%u
  end where

  do i = 1, nLambda
     meanQ(i) = sum(tmpErrorArray(1:nOuterLoop,i)%q) / real(nOuterLoop)
     meanU(i) = sum(tmpErrorArray(1:nOuterLoop,i)%u) / real(nOuterLoop)
  enddo

  do i = 1, nLambda
     sigQ(i) = sqrt(sum((tmpErrorArray(1:nOuterLoop,i)%q-meanQ(i))**2)/real(nOuterLoop-1))
     sigU(i) = sqrt(sum((tmpErrorArray(1:nOuterLoop,i)%u-meanU(i))**2)/real(nOuterLoop-1))
  enddo

  write(*,*) "setting up stokes"

  stokes_i = ytmpArray%i
  stokes_q = ytmpArray%q
  stokes_u = ytmpArray%u
  stokes_qv = (ytmpArray%i * sigQ)**2
  stokes_uv = (ytmpArray%i * sigU)**2

  dlam(1) = (xArray(2)-xArray(1))
  dlam(nlambda) = (xArray(nLambda)-xArray(nLambda-1))
  do i = 2, nLambda-1
     dlam(i) = 0.5*((xArray(i+1)+xArray(i))-(xArray(i)+xArray(i-1)))
  enddo

  if (.not.normalizeSpectrum) then
     ! convert from erg/s to erg/s/A
     stokes_i(1:nLambda) = stokes_i(1:nLambda) / dlam(1:nLambda)
     stokes_q(1:nLambda) = stokes_q(1:nLambda) / dlam(1:nLambda)
     stokes_u(1:nLambda) = stokes_u(1:nLambda) / dlam(1:nLambda)
     stokes_qv(1:nLambda) = stokes_qv(1:nLambda) / dlam(1:nLambda)**2
     stokes_uv(1:nLambda) = stokes_uv(1:nLambda) / dlam(1:nLambda)**2
     
     ! convert to erg/s/A to erg/s/cm^2/A
     
     area =  objectDistance**2  ! (nb flux is alread per sterad)
     !
     stokes_i(1:nLambda) = stokes_i(1:nLambda) / area
     stokes_q(1:nLambda) = stokes_q(1:nLambda) / area
     stokes_u(1:nLambda) = stokes_u(1:nLambda) / area
     stokes_qv(1:nLambda) = stokes_qv(1:nLambda) / area**2
     stokes_uv(1:nLambda) = stokes_uv(1:nLambda) / area**2
  

     stokes_i(1:nLambda) = stokes_i(1:nLambda) * 1.d20
     stokes_q(1:nLambda) = stokes_q(1:nLambda)  * 1.d20
     stokes_u(1:nLambda) = stokes_u(1:nLambda)  * 1.d20
     stokes_qv(1:nLambda) = stokes_qv(1:nLambda)  * 1.d40
     stokes_uv(1:nLambda) = stokes_uv(1:nLambda)  * 1.d40
     
  if (jansky) then
     do i = 1, nLambda
        stokes_i(i) = convertToJanskies(dble(stokes_i(i)), dble(xArray(i)))
        stokes_q(i) = convertToJanskies(dble(stokes_i(i)), dble(xArray(i)))
        stokes_u(i) = convertToJanskies(dble(stokes_i(i)), dble(xArray(i)))
        stokes_qv(i) = convertToJanskies(sqrt(dble(stokes_qv(i))), dble(xArray(i)))**2
        stokes_uv(i) = convertToJanskies(sqrt(dble(stokes_uv(i))), dble(xArray(i)))**2
     enddo
  endif
endif

  if (sed) then
     write(*,'(a)') "Writing spectrum as lambda F_lambda"


!     tot = 0.
!     do i = 1, nLambda
!        tot = tot + stokes_i(i)* dlam(i)
!     enddo

!     stokes_i = stokes_i / tot
!     stokes_q = stokes_q / tot
!     stokes_u = stokes_u / tot
!     stokes_qv = stokes_qv / tot**2
!     stokes_uv = stokes_uv / tot**2
     
     stokes_i(1:nLambda) = stokes_i(1:nLambda) * xArray(1:nLambda)
     stokes_q(1:nLambda) = stokes_q(1:nLambda) * xArray(1:nLambda)
     stokes_u(1:nLambda) = stokes_u(1:nLambda) * xArray(1:nLambda)
     stokes_qv(1:nLambda) = stokes_qv(1:nLambda) * xArray(1:nLambda)**2
     stokes_uv(1:nLambda) = stokes_uv(1:nLambda) * xArray(1:nLambda)**2
  endif

  if (SI) then
     write(*,'(a)') "Writing spectrum as lambda (microns) vs lambda F_lambda (microns * W/m^2)"

     tmpXarray(1:nLambda) = xArray(1:nLambda) / 1.e4
     
     stokes_i(1:nLambda) = stokes_i(1:nLambda) * tmpxArray(1:nLambda) * 10.
     stokes_q(1:nLambda) = stokes_q(1:nLambda) * tmpxArray(1:nLambda) * 10.
     stokes_u(1:nLambda) = stokes_u(1:nLambda) * tmpxArray(1:nLambda) * 10.
     stokes_qv(1:nLambda) = stokes_qv(1:nLambda) * tmpxArray(1:nLambda)**2  * 100.
     stokes_uv(1:nLambda) = stokes_uv(1:nLambda) * tmpxArray(1:nLambda)**2 * 100.
  else 
    tmpXarray = xArray
  endif

  if (velocitySpace) then
     tmpXarray(1:nLambda) = cSpeed*((xArray(1:nLambda)-lamLine)/lamLine)/1.e5 ! wavelength to km/s
  endif

  if (useNdf) then
     write(*,*) "Writing spectrum to ",trim(outfile),".sdf"
     call wrtsp(nLambda,real(stokes_i),real(stokes_q),real(stokes_qv), &
          real(stokes_u), &
          real(stokes_uv),tmpXarray,outFile)
  else
     write(*,*) "Writing spectrum to ",trim(outfile),".dat"
     open(20,file=trim(outFile)//".dat",status="unknown",form="formatted")
33   format(6(1x, 1PE14.5))
     do i = 1, nLambda
        write(20,33) tmpXarray(i),stokes_i(i), stokes_q(i), stokes_qv(i), &
             stokes_u(i), stokes_uv(i)
     enddo
     close(20)
  endif

  deallocate(ytmpArray,meanQ,meanU,sigQ,sigU,tmpErrorArray)
  deallocate(dlam,stokes_i,stokes_q,stokes_qv,stokes_u,stokes_uv)

end subroutine writeSpectrum


end module phasematrix_mod


