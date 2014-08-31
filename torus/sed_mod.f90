module sed_mod

  implicit none

  logical, public, save :: SedIsInitialised=.false.

! How to write the SED
  logical, public, save :: SedInJansky
  logical, public, save :: SedInSiUnits
  logical, public, save :: SedInLambdaFLambda

! Which inclinations to use
  integer, private, save :: SedNInc 
  real, private, save, allocatable :: SedInclinations(:)
  real, private, save, allocatable :: SedPositionAngles(:)

! Root for name of SED files
  character(len=80), save :: SedFileName

! Parameters of wavelength array used for SEDs
  real, save    :: SEDlamMin, SEDlamMax
  logical, save :: SEDwavLin
  integer, save :: SEDnumLam

  public :: setSedParameters, getSedInc, getSedPA, getNumSedInc, getSedViewVec

contains

  subroutine setSedParameters(fileName,jansky,SIsed,sed, firstPA, lastPA, &
       nInclination,firstInc,LastInc,cosSpacing,incList,PAlist,thisInclination, thisPA)
    use kind_mod
    use messages_mod
    use constants_mod

    logical, intent(in) :: jansky, SIsed, sed
    integer, intent(in), optional :: nInclination
    real, intent(in), optional    :: firstInc, LastInc
    real, intent(in)              :: firstPA, lastPA ! Don't need to be optional, default to zero if not required
    logical, optional             :: cosSpacing
    real, optional :: thisInclination, thisPA
    real, intent(in), optional    :: incList(:), PAlist(:)
    character(len=80), intent(in) :: fileName

    real(double) :: cos_inc_first, cos_inc_last, d_cos_inc, cos_inc
    integer :: i

! If this is not the first call then deallocate allocatable arrays
    if (SedIsInitialised) then
       deallocate(SedInclinations)
       deallocate(SedPositionAngles)
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
       allocate(SedPositionAngles(SedNInc))
       SedInclinations(:)   = incList(:)
       SedPositionAngles(:) = PAlist(:)

    else

! Check the required arguments are present
       if (nInclination > 1) then
          if (present(firstInc).and.present(lastInc).and.present(nInclination).and.present(cosSpacing)) then
             call writeInfo("Calculating SED inclinations", TRIVIAL)
          else
             call WriteFatal("Error in call to setSedParameters")
          end if
       endif

! Check input values and issue a warning if anything is odd  
       if (lastInc < firstInc) call writeWarning("lastInc is less than firstInc")
       if (lastPA < firstPA)   call writeWarning("lastPA is less than firstPA")

       SedNInc = nInclination
       allocate(SedInclinations(SedNInc))
       allocate(SedPositionAngles(SedNInc))

       SedInclinations(1)   = firstInc
       SedPositionAngles(1) = firstPA

       if (PRESENT(thisInclination)) then
          SedInclinations(1)   = thisInclination
          SedPositionAngles(1) = thisPA
       endif

       if ( nInclination > 1 ) then
          do i=2, SedNInc

! Inclinations, either linear spaced or cos spaced in angle
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

! Position angles, always linearly spaced
             SedPositionAngles(i) = firstPA + REAL(i-1) * &
                  (lastPA-firstPA)/REAL(nInclination-1)

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

! Return position angle i from the list
  real function getSedPA(i)
    integer, intent(in) :: i
    getSedPA = sedPositionAngles(i)
  end function getSedPA

! Return the viewing vector for the ith SED
  function getSedViewVec(i)
    use vector_mod, only: vector
    integer, intent(in) :: i
    TYPE(vector) :: getSedViewVec

    getSedViewVec%x = -sin(SedInclinations(i))
    getSedViewVec%y = 0.0
    getSedViewVec%z = -cos(SedInclinations(i))
    
  end function getSedViewVec


end module sed_mod
