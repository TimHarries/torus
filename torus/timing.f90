module timing


  public :: tune


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
    ! Modefied 11-Nov-1998 by Ryuichi Kurosawa
    ! Created 02-Jan-1992 by D. J. Hillier 
    !
    implicit none
    
    integer :: max_id                     ! max # of IDs
    parameter (max_id=40)
    integer :: lu_error                   ! file unit # for error messages
    character*(*) :: id                   ! ID name 
    character*30  :: id_stored(max_id)    ! names of id stored
    integer :: N_th_call(max_id)          ! stores how many time this routine
                                          !  was called with a given ID
    real :: cpu_start(max_id)             ! time at the first tuning pt
    real :: cpu_time                      ! elapsed CPU time
    logical :: found
    logical :: very_first_time            ! true if this routine is called
                                          !  very first time

    real :: time_array(2)                 ! dummy array to get cpu time
    data very_first_time/.true./

    integer :: i,id_num, dum_i
    integer :: minutes
    real ::  seconds
    
    save id_stored,very_first_time,cpu_start, N_th_call

    real :: etime
    external etime  ! libarary function

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

    else 
       continue
    end if
    
    ! Checks if ID is already exists in stored array (id_stored)
    found = .false.
    i=0
    do while (.not.found)
       i = i + 1
       if (i.gt.max_id) then
          write(lu_error,*) 'Error: The # of tuning point &
               &exceeds the max value set (40) -- [tune.f]'
          return
       elseif (id.eq.id_stored(i)) then
          id_num=i
          N_th_call(id_num)=N_th_call(id_num) + 1
          found = .true.
       elseif (id_stored(i).eq.' ') then ! reached the end of 
                                         ! the list and 
                                         ! could not found the matching ID.
          id_num = i                    
          id_stored(id_num)=id
          N_th_call(id_num)=1
          found = .true.
       else 
          continue
       end if
    end do

    ! Cheking if the call is the first time of the pair.
    dum_i = MOD(N_th_call(id_num),2)
    
    if (dum_i.eq.1) then ! the first call in pair 
       cpu_start(id_num) = ETIME(time_array)
       
    elseif (dum_i.eq.0) then ! the second call in pair
       cpu_time = ETIME(time_array) - cpu_start(id_num)
       minutes = NINT(cpu_time)/60
       seconds = MOD(cpu_time,60.)
       
39     continue
40     format(a12,1x,a30,1x,i6,1x, a5,1x,f9.4,1x,a5)

       open(unit=1,file='tune.dat',status='old', position = 'append')
       write(1,40) 'Tune ID:: ',id_stored(id_num),minutes, &
            'min',seconds,'sec'
       close(1)
    else
       continue
    end if
    
    return
  end subroutine tune


end module timing

