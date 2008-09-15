!
! written by tjh


module messages_mod

  use kind_mod
  use unix_mod
  
  implicit none

  interface writeFormatted
     module procedure writeFormattedReal
     module procedure writeFormattedInteger
  end interface

  public

  logical :: writeoutput
  logical :: outputWarnings
  logical :: outputInfo
  logical :: doTuning
  logical :: myRankIsZero
  integer :: verbosityLevel
  integer, parameter :: TRIVIAL = 3
  integer, parameter :: FORINFO = 2
  integer, parameter :: IMPORTANT = 1
  
contains

  subroutine writeFormattedReal(formatString, message, value, level)
    character(len=*) :: formatString, message
    real :: value
    integer :: level
    logical :: thisOutputInfo
    
    thisOutputInfo = outputInfo
    if (verbosityLevel .lt. Level) then
       thisoutputinfo = .false.
    endif

    if (writeoutput.and.thisoutputInfo) then
       write(*,formatString) "! "//trim(message), value
    endif
  end subroutine writeFormattedReal

  subroutine writeFormattedInteger(formatString, message, value, level)
    character(len=*) :: formatString, message
    integer :: value
    integer :: level
    logical :: thisOutputInfo
    
    thisOutputInfo = outputInfo
    if (verbosityLevel .lt. Level) then
       thisoutputinfo = .false.
    endif

    if (writeoutput.and.thisoutputInfo) then
       write(*,formatString) "! "//trim(message), value
    endif
  end subroutine writeFormattedInteger


  subroutine writewarning(wstring)
    character(len=*) :: wstring
    logical :: doOutput
    if (outputWarnings.and.(verbosityLevel .ge. IMPORTANT)) dooutput = .true.
    if (writeoutput.and.dooutput) then
       write(*,'(a,a,a)') "WARNING: ", trim(wstring)," !"
    endif
  end subroutine writewarning

  subroutine writefatal(wstring)
    character(len=*) :: wstring
    if (writeoutput) then
       write(*,'(a,a,a)') "FATAL ERROR: ", trim(wstring)
    endif
  end subroutine writefatal

  subroutine writeInfo(wstring, level)
    character(len=*) :: wstring
    integer, optional :: level
    logical :: thisoutputinfo

    thisoutputinfo = outputinfo
    if (present(level)) then
       if (verbositylevel .lt. Level) then
          thisoutputinfo = .false.
       endif
    endif


    if (writeoutput.and.thisoutputInfo) then
       write(*,'(a,a)') "! ",trim(wstring)
    endif
  end subroutine writeInfo


end module messages_mod




