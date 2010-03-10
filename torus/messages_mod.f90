!
! written by tjh


module messages_mod
  
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
!    if (writeoutput) then
       write(*,'(a,a,a)') "FATAL ERROR: ", trim(wstring)
!    endif
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

  subroutine writeBanner(message, cLine, level)
    integer, optional :: level
    character(len=*) :: message
    character(len=1) :: cline
    character(len=80) :: banner
    logical :: thisoutputinfo
    integer :: i, n


    thisoutputinfo = outputinfo
    if (present(level)) then
       if (verbositylevel .lt. Level) then
          thisoutputinfo = .false.
       endif
    endif

    if (writeoutput.and.thisoutputInfo) then
       n = LEN(TRIM(message))
       banner = " "
       do i = 1 , n 
          banner(i:i) = cline
       enddo
       
       write(*,'(a)') " "
       write(*,'(a)') trim(banner)
       write(*,'(a)') trim(message)
       write(*,'(a)') trim(banner)
       write(*,'(a)') " "

    endif

  end subroutine writeBanner

end module messages_mod




