module filter_set_class

  
  !  Choices of the filter sets avilable so far are:
  !  1. "sdss"           : Jim Gunn's Sloan Digital Sky Survey (SDSS) filters.
  !                        Data source: http://home.fnal.gov/~dtucker/ugriz/Filters/response.html
  !  2. "step_functions" : 1+3 different square filters at 1(+/-0.1), 10(+/-1) and 
  !                        100 (+/-10) microns.
  !
  !                        c.f. http://www.stsci.edu/hst/nicmos/design/filters
  !  3. "NIC1          " : 1+3 different square filters using similar 
  !                        bands as in HST NIC1 Medium Band Filters.
  !                        c.f. http://www.stsci.edu/hst/nicmos/design/filters


  public::  &
       make_filter_set, &
       pass_a_filter, &
       get_filter_name, &
       get_nfilter 
  

  
  private::  &
       init_filter, &
       response, &
       make_step_functions, &
       make_sdss_filters, &       
       make_nic1_filters       

       
  type filter_set
     private 
     character(LEN=30) :: name   ! For now options are "sdss" or "step_functions"     
     integer :: nfilter                ! number of filters
     type(filter), pointer  :: filters(:)  ! filter set
  end type filter_set

  
  type filter
     private
     character(LEN=30) :: name   ! the name of the filter
     ! Filter functions table
     integer ::  nlambda         ! # of elements in the table  
     ! The followin should be allocated by nf.
     double precision, pointer :: response_function(:)    !  The values should be between 0-1.
     double precision, pointer :: wavelength(:)           !  in Angstrome.
  end type filter



contains
 
  ! iniialize a new filter
  subroutine init_filter(this_filter, name, nlambda, lambda_min, lambda_max)
    implicit none
    type(filter), intent(inout)     :: this_filter
    character(LEN=*), intent(in) :: name         ! the name of the filter
    integer, intent(in)          ::  nlambda     ! # of elements in the table  
    double precision, intent(in) ::  lambda_min  !  [A]. The minimum wavelength of the filter.
    double precision, intent(in) ::  lambda_max  !  [A]. The maximum wavelength of the filter. 
    !
    integer :: i
    double precision :: del
    
    ! save some inputs
    this_filter%name = name
    this_filter%nlambda = nlambda

    ! allocate the memory for arrays.
    allocate(this_filter%response_function(nlambda))
    allocate(this_filter%wavelength(nlambda))
    
    ! initialize the elements
    this_filter%response_function(:) = 0.0d0
    this_filter%wavelength(:) = 0.0d0

    ! asssign the wavelength
    del = (lambda_max-lambda_min)/dble(nlambda-1)
    do i = 1, nlambda
       this_filter%wavelength(i) = lambda_min + del*dble(i-1)
    end do

  end subroutine init_filter


  !
  !  For a given wavelength [A] and a filter, it returns the 
  !  response of this function.
  !  Note: If the wavelength is out of range, it will 
  !        simply returns 0.
  function response(this_filter, wavelength) result(out)
    implicit none 
    double precision :: out 
    type(filter), intent(in) :: this_filter
    double precision, intent(in) :: wavelength ! in [A]
    ! 
    integer :: n, i, j 
    double precision :: tmp

    n = SIZE(this_filter%wavelength)
    if (wavelength < this_filter%wavelength(1) &
         .or.                                  &
         wavelength > this_filter%wavelength(n) ) then
       out = 0.0d0
    else
       ! interpolate the value.
       
       ! finding the index first
       j = 0
       do i = 1, n-1
          if (wavelength >= this_filter%wavelength(i) &
               .or.                                  &
               wavelength < this_filter%wavelength(i+1) ) then
             ! found the index
             j = i
             exit 
          end if
       end do

       if (j==0) then ! something went wrong.
          write(*,*) " "
          write(*,*) "Error:: j = 0 in filter_set_class::response."
          stop
       end if

       ! Doing a simple linear interpolations!       
       tmp =  (this_filter%response_function(i+1)-this_filter%response_function(i))  &
            &                                    /                                   &
            &      (this_filter%wavelength(i+1) - this_filter%wavelength(i))        
       
            
       out = tmp * (wavelength - this_filter%wavelength(i))  &
            &        + this_filter%response_function(i) 
    end if
    
  end function response
    

  ! 
  ! Initialize the filter set with a given name
  !

  subroutine make_filter_set(this_set, name)
    implicit none

    type(filter_set), intent(inout) :: this_set
    character(LEN=*),  intent(in)  :: name
    
    select case (name)
    case ("step_functions")
       call make_step_functions(this_set, name)
    case ("sdss") 
       call make_sdss_filters(this_set, name)
    case ("nic1") 
       call make_nic1_filters(this_set, name)
    case default
       write(*,*) " "
       write(*,*) "Error:: Unknown filter name passed to filter_set_class::init_filter_set."
       write(*,*) "The name was : ", TRIM(name)
       write(*,*) "Exiting the program ... "
       stop
    end select
       
  end subroutine make_filter_set


  !
  ! making the filter sets for step_functions options
  !
  !
  subroutine make_step_functions(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=4  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: natural        ! One filter to pass everything    
    type(filter) :: F001
    type(filter) :: F010
    type(filter) :: F100

    !
    ! Setting up the natural filter
    !
    ! 100 A  to  0.36 mm
    call init_filter(natural, "natural", nlam, 100.0d0, 3.6d7)
    natural%response_function(:) = 1.0d0  ! set everything to 1.0
    
    ! 
    ! F001: 0.9 to 1.1 micron
    call init_filter(F001, "F001",  nlam, 0.9d4, 1.1d4)
    F001%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! F010: 9 to 11 micron
    call init_filter(F010, "F010",  nlam, 9.0d4, 11.0d4)
    F010%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! F100: 90 to 110 micron
    call init_filter(F100, "F100",  nlam, 90.0d4, 110.0d4)
    F001%response_function(:) = 1.0d0  ! set everything to 1.0
    

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = natural
    this_set%filters(2) = F001
    this_set%filters(3) = F010
    this_set%filters(4) = F100

    ! finished

    
  end subroutine make_step_functions
  



  !
  ! Making Jim Gunn's Sloan Digital Sky Survey (SDSS) filters
  ! Data source: http://home.fnal.gov/~dtucker/ugriz/Filters/response.html
  subroutine make_sdss_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=4  ! number of filters
    type(filter) :: natural           ! One filter to pass everything    
    type(filter) :: u_filter     
    type(filter) :: g_filter
    type(filter) :: r_filter
    type(filter) :: i_filter
    type(filter) :: z_filter

    ! data arrays
    integer, parameter :: nu = 47
    real, parameter :: u(nu) =  &
         (/  &
         0.00006, 0.00012, 0.00078, 0.00328, 0.01061, &
         0.02313, 0.04150, 0.06555, 0.08696, 0.10517, &
         0.12773, 0.14634, 0.16120, 0.17775, 0.18779, &
         0.19577, 0.20279, 0.21102, 0.21718, 0.21109, &
         0.21180, 0.22420, 0.22251, 0.21329, 0.21809, &
         0.22169, 0.20859, 0.20229, 0.20748, 0.20295, &
         0.18390, 0.16873, 0.15936, 0.14399, 0.11957, &
         0.09243, 0.06573, 0.04097, 0.02178, 0.00999, & 
         0.00389, 0.00129, 0.00034, 0.00006, 0.00001, &
         0.00000, 0.00000 &
         /)
    
    integer, parameter :: ng = 45
    real, parameter :: g(ng) = &
         (/ &
         0.00000, 0.00000, 0.00000, 0.00048, 0.00893, &
         0.05384, 0.11165, 0.17106, 0.21967, 0.29012, &
         0.31400, 0.34319, 0.35647, 0.34381, 0.38717, &
         0.38128, 0.40641, 0.40627, 0.41968, 0.42600, &
         0.43143, 0.43407, 0.43410, 0.43146, 0.44100, &
         0.44549, 0.44841, 0.45192, 0.44444, 0.44881, &
         0.45457, 0.45156, 0.45005, 0.44660, 0.44057, &
         0.42125, 0.31373, 0.13353, 0.04888, 0.01955, &
         0.00731, 0.00011, 0.00000, 0.00000, 0.00000  &
         /)
   
    integer, parameter :: nr = 41
    real, parameter :: r(nr) =  &
         (/ &
         0.00000, 0.00000, 0.00001, 0.00129, 0.02186, &
         0.09726, 0.22022, 0.30848, 0.37540, 0.41713, &
         0.43044, 0.45827, 0.46276, 0.46808, 0.47567, &
         0.45958, 0.45903, 0.47507, 0.48164, 0.47961, &
         0.47151, 0.46629, 0.47223, 0.47491, 0.47146, &
         0.47198, 0.47159, 0.46501, 0.46112, 0.46193, &
         0.46033, 0.45216, 0.41676, 0.30167, 0.15345, &
         0.06400, 0.02572, 0.01069, 0.00472, 0.00224, &
         0.00112 &
         /)
         
         
    integer, parameter :: ni = 45
    real, parameter :: i(ni) =  &
         (/ &
         0.00282, 0.00281, 0.00281, 0.00280, 0.00280, &
         0.00279, 0.00693, 0.02552, 0.07164, 0.14991, &
         0.24042, 0.31680, 0.36739, 0.39594, 0.40939, &
         0.41383, 0.41305, 0.40966, 0.40575, 0.40138, &
         0.39677, 0.39004, 0.38112, 0.37065, 0.36043, &
         0.35212, 0.34597, 0.34006, 0.33312, 0.32441, &
         0.31461, 0.30492, 0.29616, 0.28844, 0.28121, &
         0.27398, 0.26657, 0.25727, 0.24744, 0.22753, &
         0.17354, 0.09048, 0.03446, 0.01226, 0.00466  &
         /)


    integer, parameter :: nz = 45
    real, parameter :: z(nz) =  &
         (/  &
         0.00003, 0.00003, 0.00013, 0.00062, 0.00248, &
         0.00845, 0.02389, 0.05252, 0.09052, 0.12920, &
         0.16093, 0.18175, 0.19355, 0.19778, 0.19711, &
         0.19312, 0.18697, 0.17949, 0.17091, 0.16157, &
         0.15156, 0.14133, 0.13084, 0.12025, 0.10949, &
         0.09879, 0.08814, 0.07773, 0.06775, 0.05800, &
         0.04841, 0.03913, 0.03062, 0.02334, 0.01735, &
         0.01233, 0.00800, 0.00443, 0.00186, 0.00049, &
         0.00003, 0.00000, 0.00002, 0.00001, 0.00000  &
         /)


    integer :: j
    
    !
    ! Setting up the natural filter
    !
    ! 100 A  to  0.36 mm
    call init_filter(natural, "natural", 5, 100.0d0, 3.6d7)
    natural%response_function(:) = 1.0d0  ! set everything to 1.0

    !
    ! u_filter
    !
    call init_filter(u_filter, "u", nu, 2980.0d0, 4130.d0)
    do j = 1, nu
       u_filter%response_function(j) = dble(u(j))
    end do

    !
    ! g_filter
    !
    call init_filter(g_filter, "g", ng, 3630.d0, 5830.d0)
    do j = 1, ng
       g_filter%response_function(j) = dble(g(j))
    end do


    !
    ! r_filter
    !
    call init_filter(r_filter, "r", nr, 5230.d0, 7230.d0)
    do j = 1, nr
       r_filter%response_function(j) = dble(r(j))
    end do

    
    !
    ! i_filter 
    !
    call init_filter(i_filter, "i", ni, 6430.d0, 8630.d0)
    do j = 1, ni
       i_filter%response_function(j) = dble(i(j))
    end do


    !
    ! z_filter 
    !
    call init_filter(z_filter, "z", nz, 7730.d0, 11230.d0)
    do j = 1, nz
       z_filter%response_function(j) = dble(z(j))
    end do


    !
    ! Now save them as a set.
    ! 
    this_set%name = name
    this_set%nfilter = 6 
    
    ALLOCATE(this_set%filters(nfilter))
    
    this_set%filters(1) = natural
    this_set%filters(2) = u_filter
    this_set%filters(3) = g_filter
    this_set%filters(4) = r_filter
    this_set%filters(5) = i_filter
    this_set%filters(6) = z_filter

    ! finished

  end subroutine make_sdss_filters



  !
  ! making the filter sets for nic1 options
  !
  !
  subroutine make_nic1_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=4  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: natural        ! One filter to pass everything    
    ! Followings are similar to HST NIC1 Medium Band Filters.
    ! See http://www.stsci.edu/hst/nicmos/design/filters
    type(filter) :: F110M
    type(filter) :: F145M
    type(filter) :: F170M

    !
    ! Setting up the natural filter
    !
    ! 100 A  to  0.36 mm
    call init_filter(natural, "natural", nlam, 100.0d0, 3.6d7)
    natural%response_function(:) = 1.0d0  ! set everything to 1.0
    
    ! F110M
    ! 1.0-1.2 microns
    call init_filter(F110M, "F110M",  nlam, 1.0d4, 1.2d4)
    F110M%response_function(:) = 1.0d0  ! set everything to 1.0

    ! F145M
    ! 1.35-1.55 microns
    call init_filter(F145M, "F145M",  nlam, 1.35d4, 1.55d4)
    F145M%response_function(:) = 1.0d0  ! set everything to 1.0

    ! F170M
    ! 1.6-1.8 microns
    call init_filter(F170M, "F170M",  nlam, 1.6d4, 1.8d4)
    F170M%response_function(:) = 1.0d0  ! set everything to 1.0

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = natural
    this_set%filters(2) = F110M
    this_set%filters(3) = F145M
    this_set%filters(4) = F170M        

    ! finished

    
  end subroutine make_nic1_filters
  




  !
  ! Given a photon, this function will return a photon after passing through 
  ! the i th filter in a filter_set object. The name of the filter set must be 
  ! also supplied as an input parameter.
  ! 
  function pass_a_filter(this_set, i, lambda) result(out)
    implicit none
    
    double precision :: out  ! betweem 0 and 1 
    !
    type(filter_set), intent(in)  :: this_set
    integer, intent(in) :: i                  ! index of the filters in this set.
    double precision, intent(in)  :: lambda    ! wavelength in Angstrome
    
    ! using a function in this module
    out = response(this_set%filters(i), lambda) 

  end function pass_a_filter



  !
  ! Returns the name of the i-th filter in a set
  !
  
  function get_filter_name(this_set, i) result(out)
    implicit none
    character(LEN=30) :: out 
    type(filter_set), intent(in) :: this_set
    integer, intent(in) :: i

    out = this_set%filters(i)%name

  end function get_filter_name
    


  ! 
  ! Returns the number of filters in a set
  !
  function get_nfilter(this_set) result(out)
    implicit none
    integer :: out 
    type(filter_set), intent(in) :: this_set
    
    out = this_set%nfilter

  end function get_nfilter
    
  
end module filter_set_class
