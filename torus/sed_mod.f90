module sed_mod

  implicit none

  logical, private, save :: SedIsInitialised=.false.

! How to write the SED
  logical, private, save :: SedInJansky
  logical, private, save :: SedInSiUnits
  logical, private, save :: SedInLambdaFLambda

! Which inclinations to use
  integer, private, save :: SedNInc 
  real, private, save, allocatable :: SedInclinations(:)

! Root for name of SED files
  character(len=80), save :: SedFileName

! Parameters of wavelength array used for SEDs
  real, save    :: SEDlamMin, SEDlamMax
  logical, save :: SEDwavLin
  integer, save :: SEDnumLam

  public :: setSedParameters, getSedInc, getNumSedInc, getSedViewVec, writeSpectrum

contains

  subroutine setSedParameters(fileName,jansky,SIsed,sed,nInclination,firstInc,LastInc,cosSpacing,incList)
    use kind_mod
    use messages_mod
    use constants_mod

    logical, intent(in) :: jansky, SIsed, sed
    integer, intent(in), optional :: nInclination
    real, intent(in), optional    :: firstInc, LastInc
    logical, optional             :: cosSpacing
    real, intent(in), optional    :: incList(:)
    character(len=80), intent(in) :: fileName

    real(double) :: cos_inc_first, cos_inc_last, d_cos_inc, cos_inc
    integer :: i

! If this is not the first call then deallocate allocatable arrays
    if (SedIsInitialised) then
       deallocate(SedInclinations)
    end if

    SedIsInitialised   = .true.
    SedInJansky        = jansky
    SedInSiUnits       = SIsed
    SedInLambdaFLambda = sed

! Set up array of inclinations, either using the supplied list of values 
! or calculate from first, last and number of values. 
    if (present(incList)) then 
       call writeInfo("Using supplied list of SED inclinations", TRIVIAL)
       SedNInc = size(incList)
       allocate(SedInclinations(SedNInc))
       SedInclinations(:) = incList(:)
    else

! Check the required arguments are present
       if (present(firstInc).and.present(lastInc).and.present(nInclination).and.present(cosSpacing)) then
          call writeInfo("Calculating SED inclinations", TRIVIAL)
       else
          call WriteFatal("Error in call to setSedParameters")
       end if

! Check input values and issue a warning if anything is odd  
       if (lastInc <= firstInc) call writeWarning("lastInc is not greater than firstInc")
       if (firstInc > piBy2)    call writeWarning("firstInc > pi/2")
       if (lastInc  > piBy2)    call writeWarning("LastInc > pi/2")
       
       SedNInc = nInclination
       allocate(SedInclinations(SedNInc))

       SedInclinations(1) = firstInc
       if ( nInclination > 1 ) then
          do i=2, SedNInc
             if (cosSpacing) then 
                cos_inc_first = COS(firstInc)
                cos_inc_last = COS(lastInc)
                d_cos_inc = (cos_inc_first - cos_inc_last)/ REAL(nInclination-1)
                cos_inc = cos_inc_first - d_cos_inc * REAL(i-1)
                SedInclinations(i) = real(max(ACOS(cos_inc),1.e-4_db))
             else
                SedInclinations(i) = firstInc + REAL(i-1) * &
                     (lastInc-firstInc)/REAL(nInclination-1)
             endif
          end do
       end if
    end if
    sedFilename = filename

  end subroutine setSedParameters

! Return the number of inclinations to be used
  integer function getNumSedInc()
    getNumSedInc = SedNInc
  end function getNumSedInc

! Return inclination i from the list
  real function getSedInc(i)
    integer, intent(in) :: i
    getSedInc = SedInclinations(i)
  end function getSedInc

! Return the viewing vector for the ith SED
  function getSedViewVec(i)
    use vector_mod, only: vector
    integer, intent(in) :: i
    TYPE(vector) :: getSedViewVec

    getSedViewVec%x = -sin(SedInclinations(i))
    getSedViewVec%y = 0.0
    getSedViewVec%z = -cos(SedInclinations(i))
    
  end function getSedViewVec

  subroutine writeSpectrum(outFile,  nLambda, xArray, yArray, varianceArray,&
       normalizeSpectrum, objectDistance, velocitySpace, lamLine)
    use kind_mod
    use messages_mod
    use constants_mod
    use phasematrix_mod
    use utils_mod, only: convertToJanskies

    implicit none
    integer, intent(in) :: nLambda
    character(len=*), intent(in) :: outFile
    real, intent(in) :: xArray(:)
    logical, intent(in) :: normalizeSpectrum, velocitySpace
    real(double), intent(in) :: objectDistance
    real(double) :: area
    real, intent(in) :: lamLine
    type(STOKESVECTOR), intent(in) :: yArray(:), varianceArray(:)
    type(STOKESVECTOR), allocatable :: ytmpArray(:),  ymedian(:)
    real(double), allocatable :: stokes_i(:), stokes_q(:), stokes_qv(:)
    real(double), allocatable :: stokes_u(:), stokes_uv(:), dlam(:)
    real, allocatable, dimension(:) :: tmpXarray
    !  real :: tot
    real :: x
    integer :: i
    character(len=80) :: message

    if (.not. SedIsInitialised) then 
       call writeWarning("SED parameters have not been initialised")
    end if

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
          x = real(1.d0/yMedian(1)%i)
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

       if (SedInJansky) then
          do i = 1, nLambda
             stokes_i(i) = convertToJanskies(dble(stokes_i(i)), dble(xArray(i)))
             stokes_q(i) = convertToJanskies(dble(stokes_i(i)), dble(xArray(i)))
             stokes_u(i) = convertToJanskies(dble(stokes_i(i)), dble(xArray(i)))
             stokes_qv(i) = convertToJanskies(sqrt(dble(stokes_qv(i))), dble(xArray(i)))
             stokes_uv(i) = convertToJanskies(sqrt(dble(stokes_uv(i))), dble(xArray(i)))
          enddo
       endif
    endif

    if (SedInLambdaFlambda) then
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

    if (SedInSiUnits) then
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
       tmpXarray(1:nLambda) = real(cSpeed*((xArray(1:nLambda)-lamLine)/lamLine)/1.e5) ! wavelength to km/s
    endif

    write(message,*) "Writing spectrum to ",trim(outfile),".dat"
    call writeInfo(message, TRIVIAL)

    open(20,file=trim(outFile)//".dat",status="unknown",form="formatted")
    if (SedInLambdaFLambda) then
       write(20,*) '# Columns are: Lambda (Angstroms) and Flux (Flux * lambda) (ergs/s/cm^2)'
    else if (SedInSiUnits) then
       write(20,*) '# Columns are: Lambda (Microns) and Flux (W/m^2)'
    else if (SedInJansky) then
       write(20,*) '# Columns are: Lambda (Microns) and Flux (janskies)'
    else
       write(20,*) '# Columns are: Lambda (Angstroms) and Flux (ergs/s/cm^2/Hz)'
    end if

33  format(6(1x, 1PE14.5))
       !34   format(a1, a14, 5(a6))
       ! You should always put a header!
       !     write(20, "(a)") "# Writteng by writeSpectrum."
       !     write(20,34)  "#", "xaxis", "I", "Q", "QV", "U", "UV"
    do i = 1, nLambda
       write(20,33) tmpXarray(i),stokes_i(i), stokes_q(i), stokes_qv(i), &
            stokes_u(i), stokes_uv(i)
    enddo
    close(20)

    deallocate(ytmpArray)
    deallocate(dlam,stokes_i,stokes_q,stokes_qv,stokes_u,stokes_uv)

  end subroutine writeSpectrum

end module sed_mod
