module isochrone_class

  !
  ! DATA from http://www.mporzio.astro.it/~dantona/prems.html
  !
  ! All the data used here are from 1998 updated version.
  !
  ! See REF:
  ! 1. F.D'Antona and I. Mazzitelli, 1997 , `Evolution of low mass stars'
  !   in ``Cool stars in Clusters and Associations'', eds. 
  !   G. Micela and R. Pallavicini, Mem. S.A.It., 68, 807
  !
  ! 2. F. D'Antona, I. Mazzitelli 1998 , `A role for superadiabatic 
  !    convection in low mass structures?',
  !    in ``Brown Dwarfs and Extrasolar Planets", ASP Conference Series, 
  !    eds. R. Rebolo, E. Martin, M.R. Zapatero Osorio, p. 442
  !
  ! 3. Censori, C. & D'Antona, F. 1998 , ``Brown Dwarfs and the luminosity
  !    function of young stellar populations", in ``Brown Dwarfs and Extrasolar Planets"
  !    ASP Conference Series, eds. R. Rebolo, E. Martin, M.R. Zapatero Osorio, p. 518 
  !


  public ::                            &
       &     new,                      &
       &     read_isochrone_data,      &
       &     get_log_T,                &
       &     get_log_L,                &
       &     allocate_arrays,          &
       &     interpol_log_L_log_T,     &
       &     mass_age_to_Temp_Rad_Lum
  
  private :: init_isochrone



  !
  ! Data
  type isochrone
     private
     ! Name of the data file.
     ! So far the choices are : 1. dam98_0215, dam98_0225, dam98_0245
     character(LEN=30) :: data_type
     !
     integer :: nmass ! The number of masses in the data
     integer :: nage  ! The number of ages in the data
     ! data array
     ! Should be allocated with size of nmass x nage.
     real, pointer :: log_T(:, :)  ! T in Kelvins
     real, pointer :: log_L(:, :)  ! L in L_sun
     ! grids
     real, pointer :: age(:)    ! [yr]  and   size = nage
     real, pointer :: mass(:)   ! [M_sun] and size = nmass     
  end type isochrone 



  !
  ! Interface for your convenience
  interface new
     module procedure init_isochrone
  end interface



contains

  
  !===================================================
  ! Constructor
  !==================================================

  ! Nullify the pointers (arrays) and save the data_type
  subroutine init_isochrone(this, data_type) 
    implicit none 
    type(isochrone), intent(inout) :: this
    character(LEN=*), intent(in) :: data_type
    !
    this%data_type = data_type

    NULLIFY(this%log_T)
    NULLIFY(this%log_L)
    
    NULLIFY(this%age)
    NULLIFY(this%mass)

  end subroutine init_isochrone
  

  ! Allocate arrays
  subroutine allocate_arrays(this, nmass, nage) 
    implicit none 
    type(isochrone), intent(inout) :: this
    integer, intent(in) :: nmass 
    integer, intent(in) :: nage
     
    !
    this%nmass = nmass
    this%nage = nage

    ALLOCATE(this%log_T(nmass, nage))
    ALLOCATE(this%log_L(nmass, nage))
    
    ALLOCATE(this%age(nage))
    ALLOCATE(this%mass(nmass))

  end subroutine allocate_arrays
  

  !====================================================
  ! Accessors
  !====================================================
  
  !
  ! Get an element of 2D array
  function get_log_T(this, i, j) RESULT(out)
    implicit none 
    real :: out 
    
    type(isochrone), intent(in) :: this
    integer, intent(in) :: i  ! the first index
    integer, intent(in) :: j  ! the second index
    
    out = this%log_T(i,j) 
  end function get_log_T

  !
  ! Get an element of 2D array
  function get_log_L(this, i, j) RESULT(out)
    implicit none 
    real :: out 
    
    type(isochrone), intent(in) :: this
    integer, intent(in) :: i  ! the first index
    integer, intent(in) :: j  ! the second index
    
    out = this%log_L(i,j) 
  end function get_log_L



  !=================================================
  ! TOOLS
  !=================================================

  ! Reads in the data from a file. 
  ! This routine must be used after the data arrays has 
  ! been allocated by subrotine new or init_isochrone
  subroutine read_isochrone_data(this) 
    implicit none
    type(isochrone), intent(inout) :: this

    integer :: FIN=33
    character(LEN=30) :: filename
    character(LEN=4)  :: dum_a
    real :: logT, logL, mass
    integer :: nage , nmass, i, j
    !
    filename(1:30) = " "
    select case (this%data_type)
    case ("dam98_0215")
       filename(1:13) = "iso0215.tot98"
    case ("dam98_0225") 
       filename(1:13) = "iso0225.tot98"
    case ("dam98_0245")
       filename(1:13) = "iso0245.tot98"
    case default
       write(*,*) "Error::Unknown data_type passed to isochrone::read_isochrone_data."
       write(*,*) "  Exiting the program ..."
       stop
    end select
       

    !
    open(unit=FIN, file=TRIM(filename), status ="old")

    
    !
    ! For the data from  F.D'Antona and I. Mazzitelli...
    
    if (this%data_type(1:5) == "dam98") then
       if (this%data_type(7:10) == "0215") then
          nmass = 25
          nage = 16
       elseif (this%data_type(7:10) == "0225") then
          nmass = 27
          nage = 16
       elseif (this%data_type(7:10) == "0245") then
          nmass = 25
          nage = 16
       else
          write(*,*) "Error:: Strange data_type in read_isochrone_data."
          stop
       end if
       
       ! allocate arrays
       ! --using a routine in this module
       call allocate_arrays(this, nmass, nage)
       
       ! read data

       do j = 1,nage
          read(FIN, *)  dum_a
          read(FIN, *)  dum_a, dum_a, this%age(j)  ! Hope this works.
          read(FIN, *)  dum_a
          read(FIN, *)  dum_a
          do i = 1, nmass
             read(FIN, *) mass,  logL, LogT
             this%mass(i) = mass     ! M_sun
             this%log_L(i,j) = logL  ! Kelvins
             this%log_T(i,j) = logT  ! L_sun
          end do
       end do

    else
       
       continue

    end if

  end subroutine read_isochrone_data


  !
  ! Find the corresponding effective temperature, radius and luminosity of a star
  subroutine mass_age_to_Temp_Rad_Lum(this, mass, age, temperature, radius, luminosity)
    implicit none
    type(isochrone),  intent(in) :: this
    double precision, intent(in) :: mass   ! should be in Solar masses
    double precision, intent(in) :: age    ! shoudle be in years
    double precision, intent(out) :: temperature  ! in  kelvins
    double precision, intent(out) :: radius       ! in  10^10 cm
    double precision, intent(out) :: luminosity   ! in  erg/sec
    !
    real :: log_L, log_T
    double precision, parameter  :: T_sun = 5778.d0  ! Kelvins
    double precision, parameter :: L_sun = 3.85d33          !erg/s
    double precision, parameter :: R_sun = 6.96d10          !cm

    !
    !  --- Using the routine in this module.
    call interpol_log_L_log_T(this, real(mass), real(age), log_L, log_T)

    temperature = (10.0**log_T)                   ! Kelvins
    
    radius =  SQRT(luminosity/L_sun) * (T_sun/temperature)**2  * R_sun ! [cm]
    radius = radius/1.0d10 ! [10^10cm]

    luminosity = (10.0**log_L) * L_sun  ! [erg/s]
    
  end subroutine mass_age_to_Temp_Rad_Lum
  



  !
  !
  ! For given mass [M_sun], age [yrs] and the isochrone data, this routine 
  ! interpolates the value of log(T[K]) and log(L[L_sun]).
  ! 
  !
  !  f(i,j) ---------- f(i, j+1)
  !     !                   !
  !     !                   !
  !     !                   !   Perform simple linear interpolations
  !     !                   !   using these four points.
  !     !                   !
  !     !                   !
  !  f(i+1,j) -------- f(i+1, j+1)
  !
  subroutine interpol_log_L_log_T(this, mass, age, log_L, log_T)
    implicit none 
    type(isochrone), intent(in) :: this  ! data
    real, intent(in) :: mass  ! [M_sun]
    real, intent(in) :: age   ! [yr]
    real, intent(out) :: log_L  ! L in [L_sun]
    real, intent(out)  :: log_T  ! T in Kelvins

    integer :: i, j, im, ia, nmass, nage
    real :: tmp, f1, f2
    
    nmass = this%nmass; nage = this%nage

    ! do simple linear interpolarions.
    ! assuming the mass and age array are in increasing order. 

    ! Fining the indecies    
    if (mass <=this%mass(1)) then
       im = 1
    elseif (mass>=this%mass(nmass)) then
       im = nmass - 1 
    else
       ! scan through the indecies
       do i = 1, nmass-1
          if (this%mass(i)< mass .and. mass <= this%mass(i+1)) then
             im = i
             exit 
          end if
       end do
    end if


    if (age <=this%age(1)) then
       ia = 1
    elseif (age>=this%age(nage)) then
       ia = nage - 1 
    else
       ! scan through the indecies
       do i = 1, nage-1
          if (this%age(i)< age .and. age <= this%age(i+1)) then
             ia = i
             exit 
          end if
       end do
       
    end if



    ! Interpolation of log(T)
    ! 
    ! interpolartes in mass first
    tmp = ( this%log_T(im+1, ia) - this%log_T(im, ia) )   &
         &                     /                          &
         &       (this%mass(im+1) - this%mass(im))
    
    f1 = tmp * (mass - this%mass(im)) + this%log_T(im, ia)


    tmp = ( this%log_T(im+1, ia+1) - this%log_T(im, ia+1) )   &
         &                     /                              &
         &       (this%mass(im+1) - this%mass(im))
    
    f2 = tmp * (mass - this%mass(im)) + this%log_T(im, ia+1)

    ! then interpolate in age
    tmp = ( f2 - f1 ) /  (this%age(ia+1) - this%age(ia))    
    log_T = tmp * (age - this%age(ia)) + f1
    


    ! Interpolation of log(L)
    ! 
    ! interpolartes in mass first
    tmp = ( this%log_L(im+1, ia) - this%log_L(im, ia) )   &
         &                     /                          &
         &       (this%mass(im+1) - this%mass(im))
    
    f1 = tmp * (mass - this%mass(im)) + this%log_L(im, ia)


    tmp = ( this%log_L(im+1, ia+1) - this%log_L(im, ia+1) )   &
         &                     /                              &
         &       (this%mass(im+1) - this%mass(im))
    
    f2 = tmp * (mass - this%mass(im)) + this%log_L(im, ia+1)

    ! then interpolate in age
    tmp = ( f2 - f1 ) /  (this%age(ia+1) - this%age(ia))    
    log_L = tmp * (age - this%age(ia)) + f1
    
  end subroutine interpol_log_L_log_T




end module isochrone_class
