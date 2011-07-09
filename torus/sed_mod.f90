module sed_mod

  implicit none

  logical, save :: SedIsInitialised=.false.

! How to write the SED
  logical, save :: SedInJansky
  logical, save :: SedInSiUnits
  logical, save :: SedInLambdaFLambda

! Which inclinations to use
  integer, private, save :: SedNInc 
  real, private, save, allocatable :: SedInclinations(:)

  public :: setSedParameters, getSedInc, getNumSedInc

contains

  subroutine setSedParameters(jansky,SIsed,sed,nInclination,firstInc,LastInc,cosSpacing,incList)
    use kind_mod
    use messages_mod
    use constants_mod

    logical, intent(in) :: jansky, SIsed, sed
    integer, intent(in), optional :: nInclination
    real, intent(in), optional    :: firstInc, LastInc
    logical, optional             :: cosSpacing
    real, intent(in), optional    :: incList(:)

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
                SedInclinations(i) = max(ACOS(cos_inc),1.e-4_db)
             else
                SedInclinations(i) = firstInc + REAL(i-1) * &
                     (lastInc-firstInc)/REAL(nInclination-1)
             endif
          end do
       end if
    end if

  end subroutine setSedParameters

  integer function getNumSedInc()
    getNumSedInc = SedNInc
  end function getNumSedInc

  real function getSedInc(i)
    integer, intent(in) :: i
    getSedInc = SedInclinations(i)
  end function getSedInc

end module sed_mod
