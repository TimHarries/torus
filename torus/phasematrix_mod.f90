
! this module sets up the stokes vector type and the corresponding
! phase matrix type. It includes subroutines to create and
! apply phase matrices

! written by tjh

! v1.0 on 16/08/99


module phasematrix_mod

  use utils_mod
  use kind_mod
  implicit none

  public

  ! the phase matrix is 4x4

  type PHASEMATRIX
     real :: element(4,4)
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
    type(PHASEMATRIX)  :: mean

    mean%element = 0.
    do i = 1 , n
       mean%element = mean%element + f(i) * a(i)%element
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
       newDirection =  newDirectionMie(oldDirection, wavelength, lamArray, &
            nLambda, miePhase, nDustType, nMuMie, dustTypeFraction)
!       write(*,*) "Scattering angle: ",acos(oldDirection.dot.newDirection)*180./pi
       meanscat=meanscat+acos(oldDirection.dot.newDirection)*180./pi
       tot=tot+newDirection
    enddo
    tot = tot / real(nTest)
    meanscat = meanscat / real(ntest)
    write(*,*) "Final vector: ",tot
    write(*,*) "mean scattering angle ",meanscat
  end subroutine testMiePhase


  type(VECTOR) function newDirectionMie(oldDirection, wavelength, &
       lamArray, nLambda, miePhase, nDustType, nMuMie, dustTypeFraction, weight)
    type(VECTOR), intent(in) :: oldDirection
    real, intent(in) :: wavelength
    real(double) :: dustTypeFraction(:)
    integer :: nDustType
    real, optional :: weight
    integer, intent(in) :: nMuMie, nLambda
    integer :: i, j, k, m, ilam
    real :: costheta, theta, r, phi
    real, intent(in) :: lamArray(:)
    real, allocatable, save :: prob(:,:)
    real, allocatable, save :: cosArray(:)
    type(VECTOR) :: tVec, perpVec, newVec
    logical,save ::firsttime = .true.
    type(PHASEMATRIX) :: miePhase(nDustType, nLambda, nMuMie)
    

    call locate(lamArray, nLambda, wavelength, ilam)
    if (ilam < 1) ilam = 1
    if (ilam > nLambda) ilam = nLambda

    if (firstTime) then
       allocate(cosArray(1:nMuMie))
       allocate(prob(1:nLambda,1:nMuMie))
       do i = 1, nMuMie
          cosArray(i) = -1. + 2.*real(i-1)/real(nMuMie-1)
       enddo
       prob = 0.
       do m = 1, nLambda
          do i = 2, nMuMie
             do k = 1, nDustType
                prob(m,i) = prob(m,i) + max(1.d-20,dustTypeFraction(k))*miePhase(k,m,i)%element(1,1)
             enddo
          enddo
          do i = 2, nMuMie
             prob(m,i) = prob(m,i) + prob(m,i-1)
          enddo
          prob(m,1:nMuMie) = prob(m,1:nMuMie)/prob(m,nMuMie)
       enddo
       firstTime = .false.
    endif

    call random_number(r)
    call locate(prob(ilam,1:nMuMie), nMuMie, r, j)
    theta = cosArray(j) + &
         (cosArray(j+1)-cosArray(j))*(r - prob(ilam,j))/(prob(ilam,j+1)-prob(ilam,j))

    if (present(weight)) then
       weight = SUM(miePhase(1:nDustType, iLam, j)%element(1,1)*dustTypeFraction(1:nDustType))/dble(nDustType) + & 
            (SUM(miePhase(1:nDustType, iLam, j+1)%element(1,1)*dustTypeFraction(1:nDustType)) - &
            SUM(miePhase(1:nDustType, iLam, j)%element(1,1)*dustTypeFraction(1:nDustType)))/dble(nDustType) * &
            (r - prob(ilam,j))/(prob(ilam,j+1)-prob(ilam,j))
    endif


    theta = acos(max(-1.,min(1.,theta)))
    call random_number(r)
    phi = twoPi * r

    tVec = oldDirection

    perpVec%x =  tVec%y
    perpVec%y = -tVec%x
    perpVec%z = 0.
    call normalize(perpVec)

    newVec = arbitraryRotate(oldDirection, theta, perpVec)
    newVec = arbitraryRotate(newVec, phi, oldDirection)


    call normalize(newvec)
    newDirectionMie = newVec
    

  end function newDirectionMie


subroutine writeSpectrum(outFile,  nLambda, xArray, yArray,  errorArray, nOuterLoop, &
     normalizeSpectrum, useNdf, sed, objectDistance, jansky, SI, velocitySpace, lamLine)

  implicit none
  integer, intent(in) :: nLambda
  character(len=*), intent(in) :: outFile
!  character(len=80) :: tfile
  real, intent(in) :: xArray(nLambda)
  logical, intent(in) :: useNdf
  logical, intent(in) :: jansky
  integer, intent(in) :: nOuterLoop
  logical, intent(in) :: normalizeSpectrum, sed, velocitySpace
  real(double), intent(in) :: objectDistance
  real(double) :: area
  real, intent(in) :: lamLine
  type(STOKESVECTOR), intent(in) :: yArray(nLambda), errorArray(nOuterloop,nLambda)
  type(STOKESVECTOR), allocatable :: ytmpArray(:), tmpErrorArray(:,:), ymedian(:)
  real, allocatable :: meanQ(:), meanU(:), sigQ(:), sigU(:)
  real(double), allocatable :: stokes_i(:), stokes_q(:), stokes_qv(:)
  real(double), allocatable :: stokes_u(:), stokes_uv(:), dlam(:), tArray(:)
  real, allocatable, dimension(:) :: tmpXarray
!  real :: tot
  real :: x
  integer :: i,j
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

  allocate(yMedian(1:nLambda))

!  do j = 1, nOuterloop
!     write(tfile,'(a,i2.2,a)') "errorarray",j,".dat"
!     open(55,file=tfile,status="unknown",form="formatted")
!     do i = 1, nLambda
!        write(55,*) xarray(i),errorArray(j,i)%i
!     enddo
!    close(55)
!  enddo

  allocate(tArray(1:nOuterLoop))
  do i = 1,nLambda
     do j = 1, nOuterLoop
        tarray(j) = errorArray(j,i)%i
     enddo
     yMedian(i)%i = median(tArray, nOuterLoop)
     do j = 1, nOuterLoop
        tarray(j) = errorArray(j,i)%q
     enddo
     yMedian(i)%q = median(tArray, nOuterLoop)
     do j = 1, nOuterLoop
        tarray(j) = errorArray(j,i)%u
     enddo
     yMedian(i)%u = median(tArray, nOuterLoop)
     do j = 1, nOuterLoop
        tarray(j) = errorArray(j,i)%v
     enddo
     yMedian(i)%v = median(tArray, nOuterLoop)
!
     yMedian(i)%i = SUM(errorArray(1:nOuterloop,i)%i)/dble(nOuterloop)
     yMedian(i)%q = SUM(errorArray(1:nOuterloop,i)%q)/dble(nOuterloop)
     yMedian(i)%u = SUM(errorArray(1:nOuterloop,i)%u)/dble(nOuterloop)
     yMedian(i)%v = SUM(errorArray(1:nOuterloop,i)%v)/dble(nOuterloop)

  enddo
  deallocate(tarray)
  do i = 1, nLambda
     write(79,*) xArray(i), ymedian(i)
  enddo


  do i = 1, nLambda
     yMedian(i) = yMedian(i) * dble(nouterloop)
  enddo

!  x = SUM(yArray(1:min(10,nLambda))%i)/real(min(10,nLambda))

!  x = 1./x

!  x = 
  write(*,*) "starting to write spectrum"

  if (normalizeSpectrum) then
     if (yMedian(1)%i /= 0.) then
        x = 1.d0/yMedian(1)%i
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
     ytmpArray(i) = yMedian(i) * x !!!!!!!!!!!!!
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
!     ! convert from erg/s to erg/s/A
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
     write(*,'(a)') "Writing spectrum as lambda (microns) vs lambda F_lambda (W/m^2)"

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
!34   format(a1, a14, 5(a6))
     ! You should always put a header!
!     write(20, "(a)") "# Writteng by writeSpectrum."
!     write(20,34)  "#", "xaxis", "I", "Q", "QV", "U", "UV"
     do i = 1, nLambda
        write(20,33) tmpXarray(i),stokes_i(i), stokes_q(i), stokes_qv(i), &
             stokes_u(i), stokes_uv(i)
     enddo
     close(20)
  endif

  deallocate(ytmpArray,meanQ,meanU,sigQ,sigU,tmpErrorArray)
  deallocate(dlam,stokes_i,stokes_q,stokes_qv,stokes_u,stokes_uv)

end subroutine writeSpectrum


subroutine plotspec(xArray, yArray, nLambda)
  implicit none
  real :: xArray(*)
  type(STOKESVECTOR) :: yArray(*)
  integer :: nLambda
  real, save :: x(1000), y(1000)
  logical, save :: firsttime = .true.

  if (firstTime) then
     firstTime = .false.
     x(1:nLambda) = xArray(1:nLambda)
     y(1:nLambda) = yArray(1:nLambda)%i
     call pgbegin(0,"/xs",1,1)
     call pgenv(x(1),x(nlambda),0.,12.,0,0)
  endif

  if (.not.firstTime) then
     call pgsci(0)
     if (y(1) /= 0.) then
        call pgbin(nlambda,x,y,.true.)
     endif
     if (yArray(1)%i /= 0.) then
        x(1:nLambda) = xArray(1:nLambda)
        y(1:nLambda) = yArray(1:nLambda)%i/yArray(1)%i
        call pgsci(1)
        call pgbin(nlambda,x,y,.true.)
     endif
  endif
end subroutine plotspec

end module phasematrix_mod


