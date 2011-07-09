module sed_mod

  implicit none

  logical, save :: SedIsInitialised=.false.

! SED parameters

  logical, save :: SedInJansky
  logical, save :: SedInSiUnits
  logical, save :: SedInLambdaFLambda

contains

  subroutine setSedParameters(jansky,SIsed,sed)

    logical, intent(in) :: jansky, SIsed, sed

    SedIsInitialised   = .true.
    SedInJansky        = jansky
    SedInSiUnits       = SIsed
    SedInLambdaFLambda = sed

  end subroutine setSedParameters

end module sed_mod
