!
! written by tjh


module messages_mod

  use kind_mod
  use unix_mod
  
  implicit none

  public

  logical :: writeoutput
  logical :: outputWarnings
  logical :: outputInfo
  logical :: doTuning
  logical :: myRankIsZero
  integer :: verbosityLevel
  integer, parameter :: TRIVIAL = 1
  integer, parameter :: FORINFO = 2
  integer, parameter :: IMPORTANT = 3
  
contains

  subroutine writewarning(wstring)
    character(len=*) :: wstring
    if (writeoutput.and.outputWarnings) then
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
       if (level < verbosityLevel) then
          thisoutputinfo = .false.
       endif
    endif


    if (writeoutput.and.thisoutputInfo) then
       write(*,'(a,a)') "! ",trim(wstring)
    endif
  end subroutine writeInfo


end module messages_mod




