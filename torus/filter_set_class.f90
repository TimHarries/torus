module filter_set_class

  
  !  Choices of the filter sets avilable so far are:
  !  1. "sdss"           : Jim Gunn's Sloan Digital Sky Survey (SDSS) filters.
  !                        Data source: http://home.fnal.gov/~dtucker/ugriz/Filters/response.html
  !  2. "step_functions" : 1+4 different square filters at 5000(+/-500) Angstrome, 
  !                        1 (+/-0.1), 10(+/-1) and 100 (+/-10) microns.
  !
  !                        c.f. http://www.stsci.edu/hst/nicmos/design/filters
  !  3. "nic1"           : 1+3 different square filters using similar 
  !                        bands as in HST NIC1 Medium Band Filters.
  !                        c.f. http://www.stsci.edu/hst/nicmos/design/filters
  !  4. "wr104"          : WR104 filters used by monnier et al.
  !                        http://www.astro.caltech.edu/mirror/keck/inst/nirc/manual/fw_macros.html
  !  5. "mips"           : Detector arrays of MIPS on SIRTF.
  !                        c.f. http://sirtf.caltech.edu/SSC/mips/documents/pocketguide.pdf
  !  6. "irac"           : IR array camera (IRAC) on SIRTF.
  !                        c.f. http://sirtf.caltech.edu/SSC/irac/spectral_response.html
  !  7. "ukirt"          : J, H, L', M' and K filters at UKIRT.
  !                        c.f. http://www.ast.cam.ac.uk:81/JACpublic/UKIRT/instruments/ircam/filters.html
  !  8. "combo"          : ukirt+irac+Vband
  !  9. "raman"          : Raman-scattering filters for the optical (6830) and UV (1032) lines
  use utils_mod

  public::  &
       make_filter_set, &
       pass_a_filter, &
       get_set_name, &
       get_filter_name, &
       get_nfilter, &
       FWHM_filters, &
       lambda_eff_filters
  

  
  private::  &
       init_filter, &
       response, &
       make_step_functions, &
       make_sdss_filters, &       
       make_nic1_filters, &
       make_wr104_filters, &
       make_mips_filters, &
       make_irac_filters, &
       make_ukirt_filters, &
       make_combo_filters

       
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
 
  ! initialize a new filter
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
    integer :: n,  i
    double precision :: tmp


    n = this_filter%nLambda

    if (wavelength < this_filter%wavelength(1) &
         .or.                                  &
         wavelength > this_filter%wavelength(n) ) then
       out = 0.0d0
    else
       ! interpolate the value.
       
       ! using a routine in utils_mod.f90
       call locate(this_filter%wavelength, n, wavelength, i)

       if (i==0) then ! something went wrong.
          write(*,*) " "
          write(*,*) "Error:: i = 0 in filter_set_class::response."
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
    case ("wr104") 
       call make_wr104_filters(this_set, name)
    case ("mips") 
       call make_mips_filters(this_set, name)
    case ("irac") 
       call make_irac_filters(this_set, name)
    case ("ukirt") 
       call make_ukirt_filters(this_set, name)
    case ("combo") 
       call make_combo_filters(this_set, name)
    case ("raman") 
       call make_raman_filters(this_set, name)
    case ("midi") 
       call make_midi_filters(this_set, name)
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
    integer,  parameter :: nfilter=5  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: natural        ! One filter to pass everything    
    type(filter) :: F5000A         ! optical (5000 Angstrome)
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
    ! F5000A: 4500 to 5500 Angstrome
    call init_filter(F5000A, "F5000A",  nlam, 4500.d0, 5500.d0)
    F5000A%response_function(:) = 1.0d0  ! set everything to 1.0

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
    F100%response_function(:) = 1.0d0  ! set everything to 1.0
    

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = natural
    this_set%filters(2) = F5000A
    this_set%filters(3) = F001
    this_set%filters(4) = F010
    this_set%filters(5) = F100

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
    integer,  parameter :: nfilter=6  ! number of filters
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
  ! making the filter sets for step_functions options
  !
  !
  subroutine make_wr104_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=4  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: natural        ! One filter to pass everything    
    type(filter) :: H
    type(filter) :: ch4
    type(filter) :: pahcs

    !
    ! Setting up the natural filter
    !
    ! 100 A  to  0.36 mm
    call init_filter(natural, "natural", nlam, 100.0d0, 3.6d7)
    natural%response_function(:) = 1.0d0  ! set everything to 1.0
    
    ! 
    ! H
    call init_filter(H, "H",  nlam, 1.5d4, 1.82d4)
    h%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! ch4
    call init_filter(ch4, "ch4",  nlam, 2.19d4, 2.34d4)
    ch4%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! pahcs
    call init_filter(pahcs, "pahcs",  nlam, 3.035d4, 3.125d4)
    pahcs%response_function(:) = 1.0d0  ! set everything to 1.0
    

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = natural
    this_set%filters(2) = H
    this_set%filters(3) = ch4
    this_set%filters(4) = pahcs

    ! finished

    
  end subroutine make_wr104_filters
  

  !
  ! making the filter sets for mips options
  !
  !
  subroutine make_mips_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=4  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: natural        ! One filter to pass everything    
    type(filter) :: F024
    type(filter) :: F070
    type(filter) :: F160

    !
    ! Setting up the natural filter
    !
    ! 100 A  to  0.36 mm
    call init_filter(natural, "natural", nlam, 100.0d0, 3.6d7)
    natural%response_function(:) = 1.0d0  ! set everything to 1.0
    
    ! 
    ! F024: center=24.0microns, width=4.7microns
    call init_filter(F024, "F024",  nlam, 21.65d4, 26.35d4)
    F024%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! F070: center=70.0microns, width=19.0microns
    call init_filter(F070, "F070",  nlam, 60.5d4, 79.5d4)
    F070%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! F160: center=160.0microns, width=35.0microns
    call init_filter(F160, "F160",  nlam, 142.5d4, 177.5d4)
    F160%response_function(:) = 1.0d0  ! set everything to 1.0



    

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = natural
    this_set%filters(2) = F024
    this_set%filters(3) = F070
    this_set%filters(4) = F160

    ! finished

    
  end subroutine make_mips_filters





  !
  ! making the filter sets for mips options
  !
  !
  subroutine make_irac_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=5  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: natural        ! One filter to pass everything    
    type(filter) :: F3_6
    type(filter) :: F4_5
    type(filter) :: F5_8
    type(filter) :: F8_0
    
    !
    ! Setting up the natural filter
    !
    ! 100 A  to  0.36 mm
    call init_filter(natural, "natural", nlam, 100.0d0, 3.6d7)
    natural%response_function(:) = 1.0d0  ! set everything to 1.0
    
    ! 
    ! F3_6: center=3.6 microns, width=0.7 microns
    call init_filter(F3_6, "F3_6",  nlam, 3.2d4, 3.9d4)
    F3_6%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! F4_5: center=4.5 microns, width=1.0 microns
    call init_filter(F4_5, "F4_5",  nlam, 4.0d4, 5.0d4)
    F4_5%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! F5_8: center=5.8 microns, width=1.5microns
    call init_filter(F5_8, "F5_8",  nlam, 5.0d4, 6.5d4)
    F5_8%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    ! F8_0: center=8.0 microns, width=3.0microns
    call init_filter(F8_0, "F8_0",  nlam, 6.5d4, 9.5d4)
    F8_0%response_function(:) = 1.0d0  ! set everything to 1.0


    

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = natural
    this_set%filters(2) = F3_6
    this_set%filters(3) = F4_5
    this_set%filters(4) = F5_8
    this_set%filters(5) = F8_0


    ! finished

    
  end subroutine make_irac_filters






  !
  ! Making J, H, L', M'  and K filters at UKIRT.
  ! Data source: http://www.ast.cam.ac.uk:81/JACpublic/UKIRT/instruments/ircam/filters.html/
  subroutine make_ukirt_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=6  ! number of filters
    type(filter) :: natural           ! One filter to pass everything    
    type(filter) :: J_filter     
    type(filter) :: H_filter
    type(filter) :: L_filter
    type(filter) :: M_filter
    type(filter) :: K_filter

    ! data arrays
    integer, parameter :: nj = 28
    real, parameter :: J(nj) =  &
         (/  &
         0.000000,  2.260830,  6.151254, 15.786206,  39.438770, &
         79.162857, 83.766869, 80.143471, 85.265762, 85.986000, &
         79.677597, 79.252739, 79.977997, 90.649467, 91.078468, &
         90.763840, 88.868996, 89.621101, 90.863686, 89.763725, &
         75.340927, 74.761604, 67.121780, 11.667536,  4.714067, &
         1.847847,  0.418457,  0.000000 &
         /)  ! this is in percentage!    
  
    integer, parameter :: nh = 36
    real, parameter :: H(nh) = &
         (/ &
         0.000000,  1.130662,  3.382868,  8.159169,  24.612747, &
         77.004608, 82.073036, 82.314682, 83.788208, 84.677483, &
         88.575111, 88.291664, 85.996643, 85.995949, 86.520782, &
         87.283333, 87.399727, 89.966148, 81.513557, 81.709839, &
         87.528709, 87.424072, 86.596039, 86.861832, 85.054771, &
         78.558357, 78.208626, 79.980522, 81.229370, 74.677246, &
         68.908279, 19.099644, 6.466859,  2.021544,   0.566167, &
         0.000000 &
        /)
   
    integer, parameter :: nl = 45
    real, parameter :: L(nl) =  &
         (/ &
         0.000000000, 0.036086982, 0.039898091, 0.031725455,  0.32689798,  &
         0.12034394,   1.4595498,   2.2350506,   4.2311155,   11.738147,   & 
         43.798082,     84.2407,   88.254423,   91.013578,   88.282997,    &
         80.101474,   85.225245,   84.099969,   87.512064,   90.767686,    &
         92.089577,   88.850656,     90.6195,   88.105127,   87.829018,    &
         87.980626,   86.376512,   87.137095,   87.484141,   86.581636,    &
         79.749162,   84.529178,   82.452241,   84.058625,   81.798227,    &
         78.115858,   76.038293,   72.831403,   72.510816,    64.65833,    &
         44.332483,   16.565997,   3.5336631,  0.39929267,  0.00000000 &
         /)
         
         
    integer, parameter :: nm = 45
    real, parameter :: M(nm) =  &
         (/ &
         0.3237915,   0.3064842,  0.57274036,   1.0353759,   2.4187509,  &
         5.9438863,   12.813169,   27.896501,    43.20953,    61.38293,  &
         72.896656,   74.012586,   75.064453,   76.547569,   77.755741,  &
         80.213367,   80.859521,   82.894614,   82.588527,   79.546416,  &
         72.752209,   75.563961,   66.493175,   59.948437,   61.247354,  &
         63.087041,    75.17032,   66.150479,   70.397813,   67.104888,  &
         64.003926,   56.702409,   51.507717,   44.577544,   27.779548,  & 
         20.47704,   14.452045,   8.9165691,   4.6622539,   2.7718573,   &
         2.3491976,    1.597735,  0.85916467,  0.57745782,  0.46482488   &
         /)


    integer, parameter :: nk = 32
    real, parameter :: K(nk) =  &
         (/  &
         0.000000,   0.787641,   1.958967,   4.547612,  10.947986, &
         58.313885,  68.592232,  77.532722,  64.455177,  72.375771,&
         86.587502,  88.620926,  92.280205,  92.903427,  88.789703,&
         88.196770,  91.462334,  91.947006,  94.770317,  93.438721,&
         95.959778,  93.058548,  94.249901,  93.387123,  89.938362,&
         67.043968,  18.038059,   9.573532,   3.788256,   1.881710,&
         0.487001,   0.000000 &
         /)

    integer :: i
    
    !
    ! Setting up the natural filter
    !
    ! 100 A  to  0.36 mm
    call init_filter(natural, "natural", 5, 100.0d0, 3.6d7)
    natural%response_function(:) = 1.0d0  ! set everything to 1.0

    !
    ! J_filter
    !
    call init_filter(J_filter, "J", nj,1.140300d04 , 1.362900d04)
    do i = 1, nj
       J_filter%response_function(i) = dble(J(i))*0.01d0
    end do

    !
    ! H_filter
    !
    call init_filter(H_filter, "H", nh, 1.4433d04, 1.8222d04)
    do i = 1, nh
       H_filter%response_function(i) = dble(H(i))*0.01d0
    end do


    !
    ! L_filter
    !
    call init_filter(L_filter, "L", nl, 3.2d04, 4.2d04)
    do i = 1, nl
       L_filter%response_function(i) = dble(L(i))*0.01d0
    end do

    
    !
    ! M_filter 
    !
    call init_filter(M_filter, "M", nm, 4.5d04, 4.9d04)
    do i = 1, nm
       M_filter%response_function(i) = dble(M(i))*0.01d0
    end do


    !
    ! K_filter 
    !
    call init_filter(K_filter, "K", nk, 1.9625d04, 2.4138d4)
    do i = 1, nk
       K_filter%response_function(i) = dble(K(i))*0.01d0
    end do


    !
    ! Now save them as a set.
    ! 
    this_set%name = name
    this_set%nfilter = 6 
    
    ALLOCATE(this_set%filters(nfilter))
    
    this_set%filters(1) = natural
    this_set%filters(2) = J_filter
    this_set%filters(3) = H_filter
    this_set%filters(4) = L_filter
    this_set%filters(5) = M_filter
    this_set%filters(6) = K_filter

    ! finished

  end subroutine make_ukirt_filters





  !
  ! making the filter sets for combo (V+ukirt+irac)
  !
  !
  subroutine make_combo_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=11  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: V
    ! subset of filters
    type(filter_set) :: ukirt
    type(filter_set) :: irac


    call init_filter(V, "V",  nlam, 5.055d3, 5.945d3)
    V%response_function(:) = 1.0d0  ! set everything to 1.0
    call make_ukirt_filters(ukirt, "ukirt")
    call make_irac_filters(irac, "irac")

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = ukirt%filters(1) ! natural
    this_set%filters(2) = V                ! V
    this_set%filters(3) = ukirt%filters(2) ! J
    this_set%filters(4) = ukirt%filters(3) ! H
    this_set%filters(5) = ukirt%filters(4) ! L
    this_set%filters(6) = ukirt%filters(5) ! M
    this_set%filters(7) = ukirt%filters(6) ! K
    this_set%filters(8) = irac%filters(2)  ! 3.6 microns
    this_set%filters(9) = irac%filters(3)  ! 4.5 microns
    this_set%filters(10) = irac%filters(4) ! 5.8 microns
    this_set%filters(11) = irac%filters(5) ! 8.0 microns


    ! finished

    
  end subroutine make_combo_filters


  !
  ! making the filter sets for raman options
  !
  !
  subroutine make_raman_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=2  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: filt6830        ! One filter to pass 6830 line
    type(filter) :: filt1032        ! UV (OVI 1032 Angstrome)

    !
    ! Setting up the 6830 filter
    !
    call init_filter(filt6830, "6830", nlam, 6800.d0, 6900.d0)
    filt6830%response_function(:) = 1.0d0  ! set everything to 1.0
    
    ! 
    ! UV filter (1032)
    call init_filter(filt1032, "1032",  nlam, 1000.d0, 1050.d0)
    filt1032%response_function(:) = 1.0d0  ! set everything to 1.0

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = filt6830
    this_set%filters(2) = filt1032

    ! finished

    
  end subroutine make_raman_filters

  !
  ! making the filter sets for raman options
  !
  !
  subroutine make_midi_filters(this_set, name)
    implicit none 
    type(filter_set), intent(inout) :: this_set
    character(LEN=*), intent(in) :: name
    !
    !
    integer,  parameter :: nfilter=4  ! number of filters
    integer,  parameter :: nlam=5  ! number of wavelegth samples
    type(filter) :: filt08        ! 8 micron passband
    type(filter) :: filt10        ! 10 micron passband
    type(filter) :: filt13        ! 13 micron passband
    type(filter) :: filtN             ! 13 micron passband

    !
    ! Setting up the 8 micron filter
    !
    call init_filter(filt08, "08", nlam, 76000.d0, 84000.d0)
    filt08%response_function(:) = 1.0d0  ! set everything to 1.0
    
    ! 
    call init_filter(filt10, "10",  nlam, 95000.d0, 105000.d0)
    filt10%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    call init_filter(filt13, "13",  nlam, 123500.d0, 136500.d0)
    filt13%response_function(:) = 1.0d0  ! set everything to 1.0

    ! 
    call init_filter(filtN, "N",  nlam, 75000.d0, 135000.d0)
    filtN%response_function(:) = 1.0d0  ! set everything to 1.0

    !
    ! Now store them as a set    
    this_set%name = name
    this_set%nfilter = nfilter 
    
    ALLOCATE(this_set%filters(nfilter))

    this_set%filters(1) = filt08
    this_set%filters(2) = filt10
    this_set%filters(3) = filt13
    this_set%filters(4) = filtN

    ! finished

    
  end subroutine make_midi_filters






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
  ! Returns the name of the filter set
  !
  
  function get_filter_set_name(this_set) result(out)
    implicit none
    character(LEN=30) :: out 
    type(filter_set), intent(in) :: this_set

    out = this_set%name

  end function get_filter_set_name
    


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


  !
  ! creates a bias and weight array based on the given filterset
  !
  subroutine createWeightArrays(thisSet, lamArray, dlam, nLambda, bias, weight, lams, lame)
    integer :: nLambda
    real(kind=doubleKind) :: lamArray(:), bias(:), weight(:), resp, dlam(:)
    integer :: i 
    real(kind=doubleKind) :: lams, lame, biasCorrection
    type(filter_set) :: thisSet

    nLambda = 0
    do i = 1, get_nfilter(thisSet)
       do j = 1, thisSet%filters(i)%nlambda
          nLambda = nLambda + 1
          lamArray(nLambda) = thisSet%filters(i)%wavelength(j)
       enddo
    enddo
    do i = 1, 2000
       nLambda = nLambda + 1
       lamArray(nLambda) = 10.**(log10(lams) + (log10(lame)-log10(lams))*real(i-1)/1999.)
    enddo
    call sort(nLambda, lamArray)
    do i = 2, nLambda-1
       dlam(i) = 0.5*((lamarray(i+1)+lamarray(i))-(lamarray(i)+lamarray(i-1)))
    enddo
    dlam(1) = lamarray(2)-lamarray(1)
    dlam(nLambda) = lamarray(nlambda)-lamarray(nLambda-1)

    do i = 1, nLambda
       resp = -1.e30
       do j = 1, get_nfilter(thisSet)
          resp = max(pass_a_filter(thisset, j, lamArray(i)), resp)
       enddo
       bias(i) = resp
    enddo
    bias(1:nLambda) = bias(1:nLambda) + 0.001
    biasCorrection =  SUM(bias(1:nLambda))
    weight(1:nLambda) = biasCorrection * (dlam(1:nLambda) / (bias(1:nLambda)))

  end subroutine createWeightArrays
    



  !
  ! Given a filter_set object and an index of an filter, this function returns
  ! the full width half maximum (FWHM) of the filter in Ansgtrom [A]
  !

  function FWHM_filters(this_set, i) RESULT(out)
    implicit none
    real(kind=doublekind) :: out   ! width in [A]
    !
    type(filter_set), intent(in) :: this_set
    integer, intent(in) :: i  ! index of the filter
    
    ! 
    real(kind=doublekind) ::  area   ! are under response curve [unit less] = 1/lambda * lambda   
    real(kind=doublekind) ::  R_max  ! Maximum value of the response curve [1/A]
   
    integer :: j
    real(kind=doublekind) :: dLambda, dA

    area = 0.0
    R_max = this_set%filters(i)%response_function(1)

    ! Integration by Trapezoidal rule
    do j = 2, this_set%filters(i)%nlambda

       dLambda = this_set%filters(i)%wavelength(j) - this_set%filters(i)%wavelength(j-1)
       dA = 0.5 * dLambda * &
            ( this_set%filters(i)%response_function(j) + this_set%filters(i)%response_function(j-1) )

       R_max = MAX(R_max, this_set%filters(i)%response_function(j))
       area = area + dA
    end do
  
    !
    ! Setting Area = 0.5*R_max*width => width = 2.0*area/R_max
    out = 2.0*area/R_max    

  end function FWHM_filters


  !
  ! Given a filter set, and the index of a filter in this set, this
  ! function will compute, the effective wavelength of the filter.
  !
  !                Integrate[L*S(L), {0,infnity}]
  ! lambda_eff =  ---------------------------------
  !                   Integrate[S(L}, {0,infinity}]
  !
  ! where L is the wavelength, and S(L) is the response function of 
  ! a filter.
  !
  function lambda_eff_filters(this_set, i) RESULT(out)
    implicit none
    real(kind=doublekind) :: out   ! width in [A]
    !
    type(filter_set), intent(in) :: this_set
    integer, intent(in) :: i  ! index of the filter
    
    ! 
    real(kind=doublekind) ::  area1, area2
    real(kind=doublekind) ::  R_max  ! Maximum value of the response curve [1/A]
   
    integer :: j
    real(kind=doublekind) :: dLambda, dA

    area1 = 0.0; area2 = 0.0
    R_max = this_set%filters(i)%response_function(1)

    ! Integration by Trapezoidal rule
    do j = 2, this_set%filters(i)%nlambda

       dLambda = this_set%filters(i)%wavelength(j) - this_set%filters(i)%wavelength(j-1)
       dA = 0.5 * dLambda * &
            ( this_set%filters(i)%response_function(j) + this_set%filters(i)%response_function(j-1) )

       R_max = MAX(R_max, this_set%filters(i)%response_function(j))
       area1 = area1 + dA
       area2 = area2 + dA*0.5*(this_set%filters(i)%wavelength(j) + this_set%filters(i)%wavelength(j-1))
    end do
  
    if (area1 > 0.0d0) then
       out = area2/area1
    else
       write(*,*) "Error:: area1 <= 0.0d0 in [filter_set_class::lambda_eff_filter]!"
       write(*,*) "        area1  = ", area1 
       write(*,*) "Exiting the program... "
       stop       
    end if


  end function lambda_eff_filters
  



  
end module filter_set_class
