module timing

  use kind_mod

  public :: tune, mySleep


contains


  subroutine tune(lu_error,ID)

    !********************************************************************
    !
    ! Unix verision of routine to collect timing statistics for 
    ! a program (major loops only).  
    !
    ! This routine calculate the elaped CPU time between two points in a parent 
    ! program.  The call to this routine has to be always 
    ! in pair.  Time here is the time between the two calls to this routine.
    ! A user must supply string to specify the tuning points, and file unit #
    ! to keep the error messages.  The output will be written in file called 
    ! "tune.dat".
    !
    ! lu_error is the logical unit for outputing the error messages.
    !
    ! NB.  
    !      SUN fortran routine ETIME gives elapsed CPU time.
    !      These library functions will return times in seconds. (See Man pages
    !      for detail.)
    !
    ! Modefied 21-Apr-2004 by Ryuichi Kurosawa :: Now uses intrinsic function (fortran95) CPU_TIME
    !             instead of ETIME. 
    ! Modefied 11-Nov-1998 by Ryuichi Kurosawa
    ! Created 02-Jan-1992 by D. J. Hillier 
    !
    implicit none
    
    integer, parameter :: max_id=500      ! max # of IDs
    integer :: lu_error                   ! file unit # for error messages
    character(len=*) :: id                   ! ID name 
    character(len=30),save  :: id_stored(max_id)    ! names of id stored
    integer,save :: N_th_call(max_id)          ! stores how many time this routine
                                          !  was called with a given ID
    real,save :: cpu_start(max_id)             ! time at the first tuning pt
    real :: cputime                       ! elapsed CPU time
    real :: clockMax                      ! Maximum time from system_clock in seconds 
    logical :: found
    logical,save :: very_first_time            ! true if this routine is called
                                          !  very first time

!    real :: time_array(2)                 ! dummy array to get cpu time
    data very_first_time/.true./

    integer :: i,id_num, dum_i
    integer :: minutes
    real ::  seconds
    

!    real :: etime
!    external etime  ! libarary function


!    CALL setrteopts ('cpu_time_type=usertime')
    ! Initiallizing arrays if used for the fisrt time
    if (very_first_time) then
       do i=1,max_id
          N_th_call(i) = 0
          cpu_start(i) = 0.0
          id_stored(i) = ' '
       end do

       open(unit=1,file='tune.dat',status='replace')
       !  -- If first time to write, an old file will be overwritten.
       close(1)
       very_first_time = .false.
    end if
    
    ! Checks if ID is already exists in stored array (id_stored)
    found = .false.
    i=0
    do while (.not.found)
       i = i + 1
       if (i.gt.max_id) then
          write(lu_error,*) 'Error: The # of tuning point &
               &exceeds the max value set (200) -- [tune.f]'
          return
       elseif (id.eq.id_stored(i)) then
          id_num=i
          N_th_call(id_num)=N_th_call(id_num) + 1
          found = .true.
       elseif (id_stored(i).eq.' ') then ! reached the end of 
                                         ! the list and 
                                         ! could not find the matching ID.
          id_num = i                    
          id_stored(id_num)=id
          N_th_call(id_num)=1
          found = .true.
       end if
    end do

    ! Cheking if the call is the first time of the pair.
    dum_i = MOD(N_th_call(id_num),2)
    
    if (dum_i.eq.1) then ! the first call in pair 
       !cpu_start(id_num) = ETIME(time_array)
!       call CPU_TIME(cpu_start(id_num))
       call wallTime(cpu_start(id_num))
       
    elseif (dum_i.eq.0) then ! the second call in pair
       !cpu_time = ETIME(time_array) - cpu_start(id_num)
!       call CPU_TIME(cputime)
       call wallTime(cpuTime, clockMax)
       cputime = cputime - cpu_start(id_num)
! Handle cases where the system clock reaches the 32bit limit and wraps.
       if (cputime < 0.0) cputime = cputime+clockMax
!       minutes = NINT(cputime)/60 ! this is wrong! (RK)
       minutes = INT(cputime/60.)
       seconds = MOD(cputime,60.)
       
! Trap the negative timing bug. DMA 15/4/14
! This should not be triggered any more but leave it in just in case ...
       if (minutes<0) then
          write(*,*) "Diagnsotic messages: minutes<0 in tune"
          write(*,*) "cputime= ", cputime
          write(*,*) "cputime/60.= ", cputime/60.
          write(*,*) "minutes= ", minutes
          write(*,*) "cpu_start=", cpu_start(id_num)
          write(*,*) "End of diagnsotic messages"
       endif

40     format(a12,1x,a30,1x,i6,1x, a5,1x,f9.4,1x,a5)

       open(unit=1,file='tune.dat',status='old', position = 'append')
       write(1,40) 'Tune ID:: ',id_stored(id_num),minutes, &
            'min',seconds,'sec'
       close(1)
    end if
    
    return
  end subroutine tune


  subroutine wallTime(rSec, clockMax)
    integer ::  count, countRate, countMax
    real, optional, intent(out) :: clockMax
    real, intent(out) :: rSec

    if (present(clockMax)) then
       call system_clock(count=count,count_rate=countRate, count_max=countMax)
       clockMax = real(countMax)/real(countRate)
    else
       call system_clock(count=count,count_rate=countRate)
    endif
    rSec = real(count)/real(countRate)

  end subroutine wallTime

  subroutine mySleep(delay)
    use messages_mod
    implicit none
    integer(bigInt) ::  count,countRate
    real(double) :: startTime, endTime, thisTime, delay
    character(len=80) :: message
    call system_clock(count=count,count_rate=countRate)


    startTime = dble(count)/dble(countRate)
    endTime = startTime + delay
    thisTime = startTime

    write(message,'(a,f7.1,a)') "Sleeping for ",delay, " seconds."
    call writeInfo(message,TRIVIAL)
    do while(thisTime < endTime)
       call system_clock(count=count,count_rate=countRate)
       thisTime = dble(count)/dble(countRate)
    end do
    call writeInfo("Done.",TRIVIAL)
  end subroutine mySleep
end module timing

