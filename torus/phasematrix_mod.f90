
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

  interface operator(+)
     module procedure addPhase
  end interface

  interface operator(-)
     module procedure subPhase
  end interface


  interface operator(*)
     module procedure multPhase
  end interface


! Private variables
  real, allocatable, private, save :: prob(:,:), probtemp(:), cosArray(:)
  logical, save, private :: setupMieDir = .true.

!$OMP THREADPRIVATE (prob, probtemp, cosArray, setupMieDir)

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
!    write(*,*) "mean element(1,1): ",mean%element(1,1)

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
       newDirection =  newDirectionMie(oldDirection, wavelength, lamArray, &
            nLambda, miePhase, nDustType, nMuMie, dustTypeFraction)
!       write(*,*) "Scattering angle: ",acos(oldDirection.dot.newDirection)*180./pi
       meanscat=meanscat+acos(oldDirection.dot.newDirection)*180./pi
       tot=tot+newDirection
    enddo
    tot = tot / dble(nTest)
    meanscat = meanscat / real(ntest)
    write(*,*) "Final vector: ",tot
    write(*,*) "mean scattering angle ",meanscat
  end subroutine testMiePhase

  subroutine resetNewDirectionMie

    if ( allocated(cosArray) ) deallocate ( cosArray )
    if ( allocated(prob)     ) deallocate ( prob     )
    if ( allocated(probtemp) ) deallocate ( probtemp )
    setupMieDir = .true.

  end subroutine resetNewDirectionMie

  type(VECTOR) function newDirectionMie(oldDirection, wavelength, &
       lamArray, nLambda, miePhase, nDustType, nMuMie, dustTypeFraction, weight)
    use utils_mod, only: locate
    type(VECTOR), intent(in) :: oldDirection
    real, intent(in) :: wavelength
    real(double) :: dustTypeFraction(:)
    integer :: nDustType
    real, optional :: weight
    integer, intent(in) :: nMuMie, nLambda
    integer :: i, j, k, m, ilam
    real(double) :: theta, phi
    real :: r
    real, intent(in) :: lamArray(:)
    type(VECTOR) :: tVec, perpVec, newVec
    type(PHASEMATRIX) :: miePhase(nDustType, nLambda, nMuMie)
    

    call locate(lamArray, nLambda, wavelength, ilam)
    if (ilam < 1) ilam = 1
    if (ilam > nLambda) ilam = nLambda

    if (setupMieDir) then
       allocate(cosArray(1:nMuMie))
       allocate(prob(1:nMuMie,1:nlambda))
       allocate(probtemp(1:nMuMie))
       do i = 1, nMuMie
          cosArray(i) = -1. + 2.*real(i-1)/real(nMuMie-1)
       enddo
       prob = 0.
       do m = 1, nLambda
          do i = 2, nMuMie
             do k = 1, nDustType
                prob(i,m) = prob(i,m) + max(1.d-20,dustTypeFraction(k))*miePhase(k,m,i)%element(1,1)
             enddo
          enddo
          do i = 2, nMuMie
             prob(i,m) = prob(i,m) + prob(i-1,m)
          enddo
          prob(1:nMuMie,m) = prob(1:nMuMie,m)/prob(nMuMie,m)
       enddo
       setupMieDir = .false.
    endif

    call randomNumberGenerator(getReal=r)

    probtemp(:) = prob(1:nMuMie,ilam)
!    call locate(prob(ilam,1:nMuMie), nMuMie, r, j)
    call locate(probtemp(:), nMuMie, r, j)
    theta = cosArray(j) + &
         (cosArray(j+1)-cosArray(j))*(r - prob(j,ilam))/(prob(j+1,ilam)-prob(j,ilam))

    if (present(weight)) then
       weight = SUM(miePhase(1:nDustType, iLam, j)%element(1,1)*dustTypeFraction(1:nDustType)) + & 
            (SUM(miePhase(1:nDustType, iLam, j+1)%element(1,1)*dustTypeFraction(1:nDustType)) - &
             SUM(miePhase(1:nDustType, iLam, j)%element(1,1)*dustTypeFraction(1:nDustType)) ) * &
            (r - prob(j,ilam))/(prob(j+1,ilam)-prob(j,ilam))
       if (weight < 0.d0) then
          write(*,*) "weight ", weight, " ilam ", ilam, " j ",j
       endif
    endif


    theta = acos(max(-1.0_db,min(1.0_db,theta)))
    call randomNumberGenerator(getReal=r)
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

  end function newDirectionMie


  subroutine fixMiePhase(miePhase, nDustType, nLambda, nMuMie)
    implicit none
    integer :: nDustType, nLambda, nMuMie
    type(PHASEMATRIX) :: miePhase(nDustType, nLambda, nMuMie)
    integer :: i, j
    logical, save :: firstTime = .true.

!$OMP THREADPRIVATE(firstTime)

    do i = 1, nDustType
       do j = 1, nLambda
          if ((miePhase(i,j,nMuMie)%element(1,1)/miePhase(i,j,nMuMie-1)%element(1,1)) > 10.d0) then
             miePhase(i,j,nMuMie)%element  = 10.d0 * miePhase(i,j,nMuMie-1)%element
             if (writeoutput.and.firstTime) then
                write(*,*) "! Undersampeld miephase fixed (near 180)"
                firstTime = .false.
             endif
          endif

          if ((miePhase(i,j,1)%element(1,1)/miePhase(i,j,2)%element(1,1)) > 10.d0) then
             miePhase(i,j,1)%element  = 10.d0 * miePhase(i,j,2)%element
             if (writeoutput.and.firstTime) then
                write(*,*) "! Undersampeld miephase fixed (near 0)"
                firstTime = .false.
             endif
          endif
       enddo
    enddo
  end subroutine fixMiePhase

  subroutine writeSpectrum(outFile,  nLambda, xArray, yArray, varianceArray,&
       normalizeSpectrum, sed, objectDistance, jansky, SI, velocitySpace, lamLine)
    use input_variables, only: useNdf
    use utils_mod, only: convertToJanskies

    implicit none
    integer, intent(in) :: nLambda
    character(len=*), intent(in) :: outFile
    !  character(len=80) :: tfile
    real, intent(in) :: xArray(nLambda)
    logical, intent(in) :: jansky
    logical, intent(in) :: normalizeSpectrum, sed, velocitySpace
    real(double), intent(in) :: objectDistance
    real(double) :: area
    real, intent(in) :: lamLine
    type(STOKESVECTOR), intent(in) :: yArray(nLambda), varianceArray(nLambda)
    type(STOKESVECTOR), allocatable :: ytmpArray(:),  ymedian(:)
    real(double), allocatable :: stokes_i(:), stokes_q(:), stokes_qv(:)
    real(double), allocatable :: stokes_u(:), stokes_uv(:), dlam(:)
    real, allocatable, dimension(:) :: tmpXarray
    !  real :: tot
    real :: x
    integer :: i
    logical :: SI
    character(len=80) :: message

    allocate(ytmpArray(1:nLambda))
    allocate(tmpXarray(1:nLambda))

    allocate(dlam(1:nLambda))
    allocate(stokes_i(1:nLambda))
    allocate(stokes_q(1:nLambda))
    allocate(stokes_qv(1:nLambda))
    allocate(stokes_u(1:nLambda))
    allocate(stokes_uv(1:nLambda))

    allocate(yMedian(1:nLambda))

    yMedian = yArray


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

    !  
    do i = 1, nLambda
       ytmpArray(i) = yMedian(i) * x !!!!!!!!!!!!!
    enddo



    stokes_i = ytmpArray%i
    stokes_q = ytmpArray%q
    stokes_u = ytmpArray%u
    stokes_qv = varianceArray%q
    stokes_uv = varianceArray%u

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
             stokes_qv(i) = convertToJanskies(sqrt(dble(stokes_qv(i))), dble(xArray(i)))
             stokes_uv(i) = convertToJanskies(sqrt(dble(stokes_uv(i))), dble(xArray(i)))
          enddo
       endif
    endif

    if (sed) then
       write(message,'(a)') "Writing spectrum as lambda F_lambda"
       call writeInfo(message, TRIVIAL)


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
       write(message,'(a)') "Writing spectrum as lambda (microns) vs lambda F_lambda (W/m^2)"
       call writeInfo(message, TRIVIAL)

       tmpXarray(1:nLambda) = xArray(1:nLambda) / 1.e4

       stokes_i(1:nLambda) = stokes_i(1:nLambda) * tmpxArray(1:nLambda) * 10.
       stokes_q(1:nLambda) = stokes_q(1:nLambda) * tmpxArray(1:nLambda) * 10.
       stokes_u(1:nLambda) = stokes_u(1:nLambda) * tmpxArray(1:nLambda) * 10.
       stokes_qv(1:nLambda) = stokes_qv(1:nLambda) * tmpxArray(1:nLambda) * 10.
       stokes_uv(1:nLambda) = stokes_uv(1:nLambda) * tmpxArray(1:nLambda) * 10.
    else 
       tmpXarray = xArray
    endif

    if (velocitySpace) then
       tmpXarray(1:nLambda) = cSpeed*((xArray(1:nLambda)-lamLine)/lamLine)/1.e5 ! wavelength to km/s
    endif

    if (useNdf) then
       write(message,*) "Writing spectrum to ",trim(outfile),".sdf"
       call writeInfo(message, TRIVIAL)
       call wrtsp(nLambda,real(stokes_i),real(stokes_q),real(stokes_qv), &
            real(stokes_u), &
            real(stokes_uv),tmpXarray,outFile)
    else
       write(message,*) "Writing spectrum to ",trim(outfile),".dat"
       call writeInfo(message, TRIVIAL)

       open(20,file=trim(outFile)//".dat",status="unknown",form="formatted")
       if (sed) then
          write(20,*) '# Columns are: Lambda (Angstroms) and Flux (Flux * lambda) (ergs/s/cm^2)'
       else if (SI) then
          write(20,*) '# Columns are: Lambda (Microns) and Flux (W/m^2)'
       else
          write(20,*) '# Columns are: Lambda (Angstroms) and Flux (ergs/s/cm^2/Ang)'
       end if

33     format(6(1x, 1PE14.5))
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

    deallocate(ytmpArray)
    deallocate(dlam,stokes_i,stokes_q,stokes_qv,stokes_u,stokes_uv)

  end subroutine writeSpectrum

end module phasematrix_mod


