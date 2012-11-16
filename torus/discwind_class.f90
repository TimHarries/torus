module discwind_class
   
  use kind_mod
  use vector_mod
  use constants_mod
  use octal_mod, only: OCTAL, subcellCentre
  use gridtype_mod, only: GRIDTYPE
  implicit none

  !  
  ! Class definition for a simple disc wind model of 
  ! Knigge et al, 1997, ApJ, 476, 291.
  !
  ! Created: june-23-2003 (Ryuich Kurosawa) 
  !
 
  public::  &
       &   new, &
       &   adddiscwind, &
       &   get_parameters, &
       &   component, &
       &   in_discwind, &
       &   all_in_discwind, &
       &   discwind_density, &
       &   discwind_Vr, &
       &   discwind_Vphi, &
       &   ave_discwind_density, &
       &   add_discwind, &
       &   turn_on_discwind, &
       &   turn_off_discwind

  private::  &
       &   int_discwind, &
       &   need_to_split, &
       &   need_to_split2, &
       &   fill_velocity_corners, &
       &   speed_of_sound, &
       &   add_new_children_discwind, &
       &   escape_velocity,&
       &   temperature_disc, &
       &   mdot_on_disc,&
       &   ave_discwind_density_slow, &
       &   ave_discwind_density_fast


  !----------------------------------------------------------
  type discwind
     private  ! DO NOT MAKE THE COMPONENTS PUBLIC!!!!!
     ! disc wind conditions
     real(double) :: d           ![10^10cm] displacement of souce point from the center of star
     real(double) :: Rmin        ! the inner most radius of the disc [10^10cm]
     real(double) :: Rmax        ! outer limit of the disc [10^10cm]
     !
     ! Temperature : T(R)= Tmax*(R/Rmin)^gamma where R is the distance from the center
     !               along the disc along the disc
     real(double) :: Tmax        ! [K] Temperature at the inner edge of the disc
     real(double) :: gamma       ! exponet in the temperature power low: 
     !
     ! mass loss rate per unit area (from the disc)
     !                                       Mdot*R^delta
     !     mdot_per_area =   ---------------------------------------------
     !                         4Pi*(Rmax ^(delta+2) - Rmin^(delta+2))
     ! where delta = 4*alpha*gamma
     !
     real(double) :: Mdot        ! [Msun/yr] total mass-loss rate 
     real(double) :: alpha       ! [-] exponent in the mass-loss rate per unit area
     !
     ! modefied  beta-velocity low
     !                                                  Rs
     !  V(r) = Cs(R) + ( f*Vesc(R) - Cs(R) ) * ( 1 - ------- )^beta
     !                                                s - Rs
     ! 
     !  Cs -- speed of sound
     !  f  -- scaling of the asymptotic terminal velocity
     !  Vesc(R) -- escape velocity from R.
     !  s -- distance from the disc along a stream line. Note this l in the paper.
     !  Rs -- constant  
     real(double) :: beta  ! [-]
     real(double) :: Rs    ! [10^10 cm]  usually 50 times of Rmin
     real(double) :: f     ! [-]  usually 2.0
     !
     ! temperature of the disc wind 
     !  -- set to be isothermal for now.
     real(double) :: Twind     ! [Kelvin] Isothermal temperature of the wind
     !
     real(double) :: Mstar ! [M_sun]  mass of the central object
     real(double) :: Hdisc ! [10^10cm]  disc height 

  end type discwind
  


  ! overload functions to make interface here ----------
  interface new
     module procedure int_discwind_default
     module procedure int_discwind
  end interface

  interface ave_discwind_density
     module procedure ave_discwind_density_max
!     module procedure ave_discwind_density_fast
!     module procedure ave_discwind_density_slow
  end interface
  

contains
  ! =====================================================================
  ! constructors
  !======================================================================

  ! initializing the discwind with default parameters
  subroutine int_discwind_default(this)
    implicit none
    type(discwind), intent(inout) :: this 
    this%d = 6.96d0*10.d0    ! [10^10cm] 10 solar radius
    this%Rmin=6.96d0*5.d0    ! the inner most radius of the disc [10^10cm]
    this%Rmax=10d0*this%Rmin ! outer limit of the disc [10^10cm]
    this%Tmax=2000.d0        ! [K] Temperature at the inner edge of the disc
    this%gamma=-0.5          ! exponet in the temperature power low: 
    this%Mdot=1.0d-8         ! [Msun/yr] total mass-loss rate 
    this%alpha=0.5d0         ! [-] exponent in the mass-loss rate per unit area
    this%beta=0.5d0          ! [-] exponent in the beta velocity law
    this%Rs=50.d0*this%Rmin  ! [10^10 cm]  usually 50 times of Rmin
    this%f=2.0d0             ! [-]  usually 2.0
    this%Twind=5000.d0       ! [Kelvin] Isothermal temperature of the wind
    this%Mstar=1.d0          ! [M_sun] mass of the central object
    this%Hdisc=0.06d0*this%Rmin   ! [10^10^cm] disc heigh (cosntant)
  end subroutine int_discwind_default
  
  
  ! inilializing with prarameters
  subroutine int_discwind(this, d, Rmin, Rmax, Tmax, gamma, &
       Mdot, alpha, beta, Rs, f, Twind, Mstar, Hdisc)
    implicit none 
    
    type(discwind), intent(inout) :: this 
    !
    ! disc wind conditions
    real(double),intent(in) :: d       ![10^10cm] displacement of souce point from the center of star
    real(double),intent(in) :: Rmin    ! the inner most radius of the disc [10^10cm]
    real(double),intent(in) :: Rmax    ! outer limit of the disc [10^10cm]
    real(double),intent(in) :: Tmax    ! [K] Temperature at the inner edge of the disc
    real(double),intent(in) :: gamma   ! exponet in the temperature power low: 
    real(double),intent(in) :: Mdot    ! [Msun/yr] total mass-loss rate 
    real(double),intent(in) :: alpha   ! [-] exponent in the mass-loss rate per unit area
    real(double),intent(in) :: beta    ! [-]
    real(double),intent(in) :: Rs      ! [10^10 cm]  usually 50 times of Rmin
    real(double),intent(in) :: f       ! [-]  usually 2.0
    real(double),intent(in) :: Twind   ! [Kelvin] Isothermal temperature of the wind
    real(double),intent(in) :: Mstar   ! [M_sun]  mass of the central object
    real(double),intent(in) :: Hdisc   ! [10^10cm]  height of the disc

    this%d = d
    this%Rmin = Rmin
    this%Rmax = Rmax
    this%Tmax = Tmax
    this%gamma = gamma
    this%Mdot = Mdot 
    this%alpha = alpha
    this%beta = beta
    this%Rs = Rs
    this%f = f
    this%Twind = Twind
    this%Mstar = Mstar
    this%Hdisc = Hdisc
    
  end subroutine int_discwind


  !==========================================================================
  ! Accessors
  !==========================================================================

  ! Given a discwind object, it returns all components of the objects.
  subroutine get_parameters(this, d, Rmin, Rmax, Tmax, gamma, &
       Mdot, alpha, beta, Rs, f, Twind, Mstar, Hdisc)
    implicit none 
    
    type(discwind), intent(in) :: this
    !
    real(double),intent(out) :: d       ![10^10cm] displacement of souce point from the center of star
    real(double),intent(out) :: Rmin    ! the inner most radius of the disc [10^10cm]
    real(double),intent(out) :: Rmax    ! outer limit of the disc [10^10cm]
    real(double),intent(out) :: Tmax    ! [K] Temperature at the inner edge of the disc
    real(double),intent(out) :: gamma   ! exponet in the temperature power low: 
    real(double),intent(out) :: Mdot    ! [Msun/yr] total mass-loss rate 
    real(double),intent(out) :: alpha   ! [-] exponent in the mass-loss rate per unit area
    real(double),intent(out) :: beta    ! [-]
    real(double),intent(out) :: Rs      ! [10^10 cm]  usually 50 times of Rmin
    real(double),intent(out) :: f       ! [-]  usually 2.0
    real(double),intent(out) :: Twind   ! [Kelvin] Isothermal temperature of the wind
    real(double),intent(out) :: Mstar ! [M_sun]  mass of the central object
    real(double),intent(out) :: Hdisc   ! [10^10cm]  height of the disc

    d      =      this%d 
    Rmin   =   this%Rmin 
    Rmax   =   this%Rmax 
    Tmax   =   this%Tmax 
    gamma  =  this%gamma 
    Mdot   =   this%Mdot 
    alpha  =  this%alpha 
    beta   =   this%beta 
    Rs     =     this%Rs 
    f      =      this%f 
    Twind  =  this%Twind
    Mstar  =  this%Mstar
    Hdisc  =  this%Hdisc


  end subroutine get_parameters


  !
  ! Given a discwind object and the name of its component, this function 
  ! returns the component
  real(double) function component(this, name) 
    implicit none 
    
    type(discwind), intent(in) :: this 
    character(LEN=*), intent(in) :: name 
    
    select case (name)
    case ("d")
       component =      this%d 
    case ("Rmin")
       component   =   this%Rmin 
    case ("Rmax")
       component   =   this%Rmax 
    case ("Tmax")
       component   =   this%Tmax 
    case ("gamma") 
       component  =  this%gamma 
    case ("Mdot")
       component   =   this%Mdot 
    case ("alpha")       
       component  =  this%alpha 
    case ("beta") 
       component   =   this%beta 
    case ("Rs")
       component     =     this%Rs 
    case ("f")
       component      =      this%f 
    case ("Twind")
       component  =  this%Twind        
    case ("Mstar")
       component  =  this%Mstar
    case ("Hdisc")
       component  =  this%Hdisc
    case default
       write(*,*) "Error:: Unknown NAME was passed to [discwind_class::component]."
       write(*,*) "Exiting the program because of this error!"
       stop
    end select

  end function component


  !=============================================================================
  ! Tools
  !=============================================================================

  !
  ! disc temperature (at z=0) power law  [K]
  ! Temperature : T(R)= Tmax*(R/Rmin)^gamma where R is the distance from the center
  ! (private)
  function temperature_disc(this, R) RESULT(out)
    implicit none
    real(double) :: out
    type(discwind), intent(in) :: this
    real(double), intent(in) :: R ! [10^10cm] of distance from z-axis
    
    out = this%Tmax * (r/this%Rmin)**this%gamma  ! [K]

    if (out < 3.0d0) out = 3.0d0  ! setting the minimum temperature

  end function temperature_disc


  !
  ! Approximated speed of sound in [cm/s]
  ! (private)
  real(double) function speed_of_sound(this, r)  
    implicit none
    type(discwind), intent(in) :: this
    real(double), intent(in) :: R ! [10^10cm] of distance from z-axis
        
    ! using the function in this module
    speed_of_sound = 10.0d0*SQRT(temperature_disc(this, R)/1.0d4) ! [km/s]
    speed_of_sound = speed_of_sound*1.d5  ! [cm/s]

  end function speed_of_sound

  !
  ! Vesc(R) = SQRT(2GM/R)
  ! (private) 
  real(double) function escape_velocity(this, r) 
    implicit none
    type(discwind), intent(in) :: this
    real(double), intent(in) :: r ! [10^10cm] of distance from the center of the star
    real(double) :: Radius
    real(double) :: Mass
    
    Radius = r*1.0d8  ! [m]
    Mass = this%Mstar*1.98892d30   ! [kg]

    escape_velocity = SQRT(2.0d0*6.673d-11*Mass/Radius) * 1.d2 ! [cm/s]
    

  end function escape_velocity




  !
  ! Modefied beta-velocity law: 
  !                                                  Rs
  !  V(r) = Cs(R) + ( f*Vesc(R) - Cs(R) ) * ( 1 - ------- )^beta
  !                                                s - Rs
  !   r -- the distace measured from the source point located at the distnace d below
  !        or above the origin along z axis.
  !  Cs -- speed of sound
  !  f  -- scaling of the asymptotic terminal velocity
  !  Vesc(R) -- escape velocity from R.
  !  s -- distance from the disc along a stream line. Note this l in the paper.
  !  Rs -- constant    
  !
  ! Will rerurn the value in [cm/s].
  !
  real(double) function discwind_Vr(this, x, y, z)
    implicit none
    !
    type(discwind), intent(in) :: this
    real(double), intent(in) :: x, y, z      ! [10^10 cm] coordinates of a point
    !
    real(double) :: r  ! [10^10cm]  distance from the souce point
    real(double) :: s  ! [10^10cm]  distance from the disc
    real(double) :: d, rho  ! [10^10cm]  
    real(double) :: rho_sq  ! [10^20 cm]
    real(double) :: Cs       ! speed of sound in [cm/s]
    real(double) :: fac, tmp

    d = this%d  ! [10^10cm]
    rho_sq = x*x + y*y; rho = SQRT(rho_sq)
    !
    r = SQRT(rho_sq + (ABS(z)+d)**2)  ! [10^10cm]
    
    ! distance from the disc to the field point along 
    ! the stream line (straight line connection the source point 
    ! and the field line.)
    ! -- using a proportionality in similar triangles
    s = r*(1.0d0 - d/(ABS(z)+d))  ! [10^10cm]

    ! 
    Cs = speed_of_sound(this, rho)  ! [cm/s]    
    fac = 1.0d0 -  this%Rs/(s+this%Rs)  ! [-]
    fac = fac**this%beta

    tmp = Cs + (this%f*escape_velocity(this, rho) - Cs)*fac  ! [cm/s]
    discwind_Vr = tmp
 
  end function discwind_Vr


  !
  ! Assuming that angular momentum around z-axis conserves...
  ! The output is in cm/s.
  real(double) function discwind_Vphi(this, x, y, z)
    implicit none
    !
    type(discwind), intent(in) :: this
    real(double), intent(in) :: x, y, z      ! [10^10 cm] coordinates of a point
    !
    real(double) :: r  ! [10^10cm]  distance from the souce point
    real(double) :: d, rho  ! [10^10cm]  
    real(double) :: rho_sq  ! [10^20 cm]
    real(double) :: rho_i, z_dist, Vk
    real(double), parameter :: G = 6.67259d-8 ! in cgs
    real(double) :: Mstar

    
    d = this%d  ! [10^10cm]
    rho_sq = x*x + y*y; rho = SQRT(rho_sq)
    !
    r = SQRT(rho_sq + (ABS(z)+d)**2)  ! [10^10cm]


    ! finding the polar radius of the 
    ! initial wind insection point on the disc by using geometry.
    ! -- using similar triangles.
    z_dist = d+ABS(z)
    rho_i = (d/z_dist) * rho

    ! Now compute the Keplerian orbital speed at this radius


    Mstar = this%Mstar * 1.9891d33  ! converting from [Msun] to [g]

    if (rho_i/=0) then
       Vk = SQRT( G*Mstar / (rho_i*1.0d10) )  ! [cm/s]
    else
       Vk = 0.0d0
    end if
    
    ! Now using the conservation of angular momentum around z...   
    if (rho /= 0.0d0) then
       discwind_Vphi = Vk * (rho_i/rho)   ! [10^10cm]
    else
       discwind_Vphi = 0.0d0  ! [10^10cm]
    end if
 
  end function discwind_Vphi
    

  !
  ! Given a  position vector  [10^10cm] and this objects, this routine returns the  
  ! the velocity vector. componets are in [c] speed of lights
  !
  type(vector) function discwind_velocity(this, r_vec) 
    implicit none
    type(discwind), intent(in) :: this
    type(vector), intent(in) :: r_vec
    !
    !
    real(double), parameter :: C =2.99792458d10   ! speed of light in [cm/s]
    !
    real(double) :: x,y,z   ! should be in [10^10cm]
    real(double) :: rho     ! cylindical radius [10^10cm]
    real(double) :: rho_sq  ! cylindical radius squared [10^20cm]
    real(double) :: rp      ! distance from a souce point to a field point [10^10cm]
    type(vector) :: rp_vec  ! position vector orinined from a source points
    type(vector) :: dir_vec ! direction vector (from source)
    real(double) :: d       ! distance from the z=0 plane to the source
    real(double) :: Vr,Vphi ! 
    real(double) :: phi, theta
    real(double) :: Vx, Vy, Vz
    
    
    ! position
    x=r_vec%x; y=r_vec%y; z=r_vec%z   ! [10^10cm]
    rho_sq = x*x + y*y; rho = SQRT(rho_sq)
    d = this%d  ! [10^10cm]

    if (in_discwind(this, x, y, z) ) then 
       if (z>0) then
          rp_vec = VECTOR(x, y, z+d)
       else
          rp_vec = VECTOR(x, y, z-d)
       end if

       rp = modulus(rp_vec)  ! normalizetion
       if (rp /=0.0d0) then
          dir_vec = rp_vec/rp
       else
          ! this could happen at the poles
          dir_vec = VECTOR(0.0, 0.0, 1.0)
       end if

       Vr = discwind_Vr(this,x, y, z) /c   ! in [c] the unit of the speed of light    
       Vphi = discwind_Vphi(this,x, y, z) /c   ! in [c] the unit of the speed of light    
       
       phi = ATAN2(y,x)
       theta = ACOS( (ABS(z)+d) / rp )  ! origin shifted
       
       Vx = Vr*SIN(theta)*COS(phi) - Vphi*SIN(phi)      ! [c]
       Vy = Vr*SIN(theta)*SIN(phi) + Vphi*COS(phi)      ! [c]
       Vz = Vr*COS(theta)                               ! [c]
       if (z <0 ) Vz = -Vz

       discwind_velocity= VECTOR(Vx, Vy, Vz)  ! [c]

    else

       discwind_velocity = VECTOR(1.e-30,1.e-30,1.e-30)  ! (magnitude x direction)       

    end if

  end function discwind_velocity




  ! 
  ! Function to check if a given point (xpos, ypos, zpos) are in the discwind zone
  ! 
  logical function in_discwind(this, xpos,ypos,zpos, subcellsize)
    implicit none
    type(discwind), intent(in) :: this
    real(double), intent(in) :: xpos,ypos,zpos  ! should be in [10^10cm]
    real(double), optional, intent(in) :: subcellsize
    !
    real(double) :: x, y, z     ! [10^10cm]  position with respect to the center of star
    real(double) :: R           ! [10^10cm]  cylindical radius 
    real(double) :: s1, s2      ! slope of the inner and outer edge of the wind
    real(double) :: R1, R2      ! The range of the R allowed for the win
    real(double) :: dr
    ! position with respect to the center of the star
    x = xpos
    y = ypos
    z = zpos

    ! distance from the center of the star on the disc (z=0)
    R = SQRT(x*x + y*y)

    ! slopes of the two boundaries
    s1 = this%d/this%Rmin
    s2 = this%d/this%Rmax
    if (s1 ==0)  then
       write(*,*) "Error:: s1 = 0 in [discwind_class::in_discwid]."
       stop
    end if
    if (s2 ==0)  then
       write(*,*) "Error:: s2 = 0 in [discwind_class::in_discwid]."
       stop
    end if

    ! the boundaries
    R1 = (ABS(z)+this%d)/s1
    R2 = (ABS(z)+this%d)/s2


    ! Now checks if in a vaild zone
    if (PRESENT(subcellsize)) then
       ! this will allow to split
       ! the cell at the edge
       dr = 2.0d0*subcellsize
    else
       dr = 0.0d0
    end if
    if ( (R+dr) >R1 .and. (R-dr)<R2 ) then
       in_discwind = .true.
    else
       in_discwind = .false.
    end if


!    ! set it false if in the disc
!    if (ABS(z) < this%Hdisc/2.0d0) in_discwind = .false.

  end function in_discwind


  logical function all_in_discwind(thisOctal, subcell, this)
    type(DISCWIND) :: this
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: x, y, z, d, xcen
    type(VECTOR) :: cVec

    all_in_discwind = .true.
    cVec = subcellCentre(thisOctal,subcell)
    d = thisOctal%subcellSize/2.d0 * 0.999d0
    xcen = sqrt(cVec%x**2 + cVec%y**2)

    x = xcen+d; y = 1.d-20; z = cvec%z+d
    all_in_discwind = all_in_discwind.and.in_discwind(this, x, y, z)

    x = xcen+d; y = 1.d-20; z = cvec%z-d
    all_in_discwind = all_in_discwind.and.in_discwind(this, x, y, z)

    x = xcen-d; y = 1.d-20; z = cvec%z+d
    all_in_discwind = all_in_discwind.and.in_discwind(this, x, y, z)

    x = xcen-d; y = 1.d-20; z = cvec%z-d
    all_in_discwind = all_in_discwind.and.in_discwind(this, x, y, z)

  end function all_in_discwind

  !
  !
  ! Mass loss rate per unit are on the disc surface [g/s/cm^2]
  function mdot_on_disc(this, rho) RESULT(out)
    implicit none
    real(double) :: out
    !
    type(discwind), intent(in) :: this
    real(double), intent(in) :: rho ! [10^10cm]  cylindical diatance from the center
    !
    real(double) :: Mdot, Rmax, Rmin, R, delta, fac, p
    real(double), parameter ::  pi = 3.1415927d0
    
    delta =4.0d0*this%alpha*this%gamma
    p = delta + 2.0d0
    ! converting everthing in cgs
    R = rho*1.d10  !  [cm]
    Rmin = this%Rmin * 1.0d10 ! converting from [10^10cm] to [cm]
    Rmax = this%Rmax * 1.0d10 ! converting from [10^10cm] to [cm]
    Mdot = this%Mdot*1.989d33/3.1536d07   ! Converting from [Msun/yr] to [g/s]

    fac = p / ( Rmax**p - Rmin**p )

    ! just checking
    if (fac < 0.0d0) then
       write(*,*) "Error:: fac <0.0 in [discwind_class::mdot_on_disc]."
       stop
    end if

    out = Mdot*(R**delta)*fac/(4.0d0*pi)   ! [g/s/cm^2]

  end function mdot_on_disc


          
  ! Given a position (x, y and z) in 10^10 cm, this routine
  ! will return the density [g/cm^3] at the point.
  !
  function discwind_density(this, xpos, ypos, zpos) RESULT(out)
    implicit none
    real(double) :: out  ! in [g/cm^3]
    !
    type(discwind), intent(in) :: this
    real(double), intent(in) :: xpos,ypos,zpos  ! should be in [10^10cm]
    !
    real(double) :: x, y, z, d     ! [cm]  
    real(double) :: rho    ! cylyndical distance [cm]
    real(double) :: r      ! distance from a souce displaced [cm]
    real(double), parameter :: density_min = tiny(density_min) ! [g/cm^3]
    real(double), parameter :: density_max = 1.0d10            ! [g/cm^3]
    real(double) :: rho_sq, s, rp, fac

    if ( in_discwind(this, xpos, ypos, zpos) ) then
       ! position with respect to the center of the star
       x = xpos*1.0d10  ! [cm]
       y = ypos*1.0d10  ! [cm]
       z = zpos*1.0d10   ! [cm]
       d = this%d*1.0d10 ! [cm]

       ! distance from the z-axis
       rho_sq = x*x + y*y
       rho = SQRT(rho_sq)  ! [cm]

       ! distance from the source
       r  = SQRT(rho_sq + (ABS(z)+d)**2)  ! [cm]

       ! distance from the disc to the field point along
       ! the stream line (straight line connection the source point
       ! and the field line.)
       ! -- using a proportionality in similar triangles
       s = r*(1.0d0 - d/(ABS(z)+d))  ! [cm]

       rp = r-s
       fac = rp/r
       fac = fac*fac
       ! using the functions in this module
       out = mdot_on_disc(this, rho/1.0d10)*fac  &
            / discwind_Vr(this, xpos, ypos, zpos)   ! [g/cm^3]

    else
       out = density_min    
    end if

    if (out > density_max) out=density_max

    
  end function discwind_density



  

  !
  !  Given a grid object, this routine will add extra grid and assign density (and etc) 
  !  to the grid using the density function defined in this module.
  !  


  subroutine addDiscWind(grid)
    use inputs_mod, only : DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
       DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, limitscalar, ttauriMstar
    real(double) :: DW_Hdisc
    type(GRIDTYPE) :: grid
    type(DISCWIND) :: myDiscWind

    call new(mydiscWind, DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
       DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, dble(ttauriMstar)/msol, DW_Hdisc)
    call add_discwind(grid%octreeRoot, grid, myDiscWind, limitscalar)


  end subroutine addDiscWind



  RECURSIVE SUBROUTINE add_discwind(thisOctal,grid,this, limitscalar)
    ! uses an external function to decide whether to split a subcell of
    !   the current octal. 

    use inputs_mod, only : maxDepthAMR
    IMPLICIT NONE
    TYPE(OCTAL), POINTER :: thisOctal
    TYPE(gridtype), INTENT(INOUT) :: grid   ! need to pass the grid through to the 
    type(discwind), intent(in)  :: this     
    ! the value the decideSplit function uses to  decide whether or not to split cell.
    real(double), intent(in) :: limitscalar 
    !
    !
    TYPE(OCTAL), POINTER :: child
    INTEGER              :: i, j, k    ! loop counters
    integer :: iindex
    integer :: isubcell


    DO iSubcell = 1, thisOctal%maxChildren

       if (thisOctal%hasChild(isubcell)) cycle
       IF (need_to_split2(thisOctal,isubcell,this).and.(thisOctal%nDepth < maxDepthAMR)) then

          CALL add_new_children_discwind(thisOctal, isubcell, grid, this)

          if (.not.thisOctal%hasChild(isubcell)) then
             write(*,*) "add child failed in splitGrid"
             do
             enddo
          endif
          
       END IF
       
    END DO

  do i = 1, thisOctal%nChildren
     if (.not.thisOctal%hasChild(thisOctal%indexchild(i))) then
        write(*,*) "octal children messed up"
        do ; enddo
        endif
     enddo

    do i = 1, thisOctal%maxChildren
      k = -99
      if (thisOctal%hasChild(i)) then
        do j = 1, thisOctal%nChildren
          if (thisOctal%indexChild(j) == i) then
            k = j
            exit
          endif
       enddo
       if (k==-99) then
          write(*,*) "subcell screwup"
          do
          enddo
       endif
    endif
 enddo

   if (any(thisOctal%haschild(1:thisOctal%maxChildren)).and.(thisOctal%nChildren==0)) then
      write(*,*) "nchildren screw up"
      do;enddo
      endif
    
    DO iIndex = 1, thisOctal%nChildren
       child => thisOctal%child(iIndex)
       CALL add_discwind(child,grid,this,limitscalar)      
   END DO
666 continue    
  END SUBROUTINE add_discwind


  !
  !
  ! This is based decideSplit in amr_mod.f90
  logical function need_to_split(thisOctal,subcell, this, limitscalar)
    
    IMPLICIT NONE
    
    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(discwind), intent(in) :: this
    ! the value the decideSplit function uses to  decide whether or not to split cell.
    real(double), intent(in) :: limitscalar 
    
    
    real(oct)  :: cellSize
    real(double) :: rho_discwind,  mass
    TYPE(VECTOR)     :: cellCentre 

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize*2.0d0
    cellCentre = subcellCentre(thisOctal, subcell)
    
    rho_discwind = ave_discwind_density(thisOctal, subcell, this)  ! [g/cm^3]
    ! the mass of the cell
    mass = rho_discwind * (cellsize*1.0d10)**3          ! [g]
    
    if (mass > limitscalar) then
       need_to_split = .true.  
    else
       need_to_split = .false.
    end if
               

  end function need_to_split


  !
  ! Split using the log-scaled radial grid.
  logical function need_to_split2(thisOctal,subcell, this)
    use utils_mod, only: locate 

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(discwind), intent(in) :: this
    
    real(oct)  :: cellSize, d
    TYPE(VECTOR)     :: cellCentre
!    integer, parameter :: nr = 150  ! normal resolution
    integer, parameter :: nr = 180  ! normal resolution
!    integer, parameter :: nr = 40  ! low resolution

    real(double) :: r
    integer :: i
    !
    logical, save :: first_time = .true.
    real(double) , save:: rGrid(nr)
    real(double) :: Rmax = 1.5d5  ! [10^10cm] = 100 AU
    TYPE(VECTOR)     :: VecInnerEdge
    real(double) :: wi

    need_to_split2 = .false.

    if (first_time) then  ! setup ref. r grid
       do i = 1, nr
          ! this should not depend on the size of the model boundary box.
          rGrid(i) = log10(this%Rmin)+dble(i-1)/dble(nr-1)*(log10(Rmax)-log10(this%Rmin))
!          rGrid(i) = log10(this%Rmin)+dble(i-1)/dble(nr-1)*(log10(this%Rmax)-log10(this%Rmin))
       enddo
       do i = 1, nr
          rGrid(i) = 10.d0**rGrid(i)
       end do
       first_time = .false.
    end if

    cellSize = (thisOctal%subcellSize)*2.0d0
    cellCentre = subcellCentre(thisOctal, subcell)


!    if (.not.in_jet_flow(this, cellCentre) ) then
    if (.not.in_discwind(this, cellCentre%x, cellCentre%y, cellCentre%z, thisOctal%subcellSize) ) then
       need_to_split2 = .false.
    else
       ! get the size and the position of the centre of the current cell
!       r = modulus(cellCentre)  ! used in the paper
       wi = sqrt(cellCentre%x*cellCentre%x + cellCentre%y*cellCentre%y)
       VecInnerEdge = VECTOR(cellCentre%x, cellCentre%y, 0.0d0)* (this%Rmin/wi)
       r = modulus(cellCentre-VecInnerEdge)  ! shift it to the inner edge of the disc
!       r = ABS(wi-this%Rmin)  ! just a cylindical radius
       call locate(rGrid,nr,r,i)
       if (i > (nr-1)) i = nr-1
       d = rGrid(i+1) - rGrid(i)
       if (cellSize > d ) then
          need_to_split2 = .true.
       end if
    endif
  end function need_to_split2


  

  real(double) function ave_discwind_density_slow(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(discwind), intent(in) :: this
    
    real(oct)  :: cellSize, x, y, z
    TYPE(VECTOR)     :: cellCentre 
    real(double) :: xc, yc, zc, d, r, rho_sum, rho_sample
    integer, parameter :: nsample = 200
    integer :: i

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize
    cellCentre = subcellCentre(thisOctal, subcell)

    xc = cellcentre%x
    yc = cellcentre%y
    zc = cellcentre%z
    d = cellsize/2.0d0

    rho_sum=0.0d0
    if (thisOctal%threeD)  then
       do i = 1, nsample
          ! picking a random position in the subcell
          call randomNumberGenerator(getDouble=r)
          x =  xc - d + r*2.0d0*d
          call randomNumberGenerator(getDouble=r)
          y =  yc - d + r*2.0d0*d
          call randomNumberGenerator(getDouble=r)
          z =  zc - d + r*2.0d0*d               
          ! evaluate the disc density density
          ! Taking the max value samples
          rho_sample =  discwind_density(this, x, y, z)
          rho_sum = rho_sum + rho_sample
       end do
    else  ! 2D amr
       do i = 1, nsample
          ! picking a random position in the subcell
          call randomNumberGenerator(getDouble=r)
          x =  xc - d + r*2.0d0*d
          y =  yc 
          call randomNumberGenerator(getDouble=r)
          z =  zc - d + r*2.0d0*d               
          ! evaluate the disc density density
          ! Taking the max value samples
          rho_sample =  discwind_density(this, x, y, z)
          rho_sum = rho_sum + rho_sample
       end do
    end if

    ave_discwind_density_slow = rho_sum/dble(nsample)

  end function ave_discwind_density_slow


  !
  !
  !
  real(double) function ave_discwind_density_fast(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(discwind), intent(in) :: this
    
    real(oct)  :: cellSize, x, y, z
    TYPE(VECTOR)     :: cellCentre 
    real(double) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, d

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize
    cellCentre = subcellCentre(thisOctal, subcell)
    x = cellcentre%x
    y = cellcentre%y
    z = cellcentre%z
    d = cellsize/3.00
    if (.not. thisOctal%threeD)  y = 0.0d0  ! 2D case

    ! assigning density 
    if (thisOctal%threeD)  then
       c0 =  discwind_density(this, x, y, z)    
       c1 =  discwind_density(this, x+d, y+d, z+d)
       c2 =  discwind_density(this, x+d, y+d, z-d)
       c3 =  discwind_density(this, x+d, y-d, z+d)
       c4 =  discwind_density(this, x+d, y-d, z-d)
       c5 =  discwind_density(this, x-d, y+d, z+d)
       c6 =  discwind_density(this, x-d, y+d, z-d)
       c7 =  discwind_density(this, x-d, y-d, z+d)
       c8 =  discwind_density(this, x-d, y-d, z-d)
       
       ave_discwind_density_fast =  &
            (c0 + c1+ c2 + c3 + c4 + c5 + c6 + c7 + c8)/9.0d0
    else  ! 2D case
       c0 =  discwind_density(this, x, y, z)    
       c1 =  discwind_density(this, x+d, y, z+d)
       c2 =  discwind_density(this, x+d, y, z-d)
       c3 =  discwind_density(this, x-d, y, z+d)
       c4 =  discwind_density(this, x-d, y, z-d)
       
       ave_discwind_density_fast =  &
            (c0 + c1+ c2 + c3 + c4)/5.0d0
    end if

  end function ave_discwind_density_fast




  !
  !
  !
  real(double) function ave_discwind_density_max(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(discwind), intent(in) :: this
    
    real(oct)  :: cellSize, x, y, z
    TYPE(VECTOR)     :: cellCentre 
    real(double) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, d

    ! get the size and the position of the centre of the current cell
    cellSize = thisOctal%subcellSize
    cellCentre = subcellCentre(thisOctal, subcell)
    x = cellcentre%x
    y = cellcentre%y
    z = cellcentre%z
    d = cellsize/3.00
    if (.not. thisOctal%threeD)  y = 0.0d0  ! 2D case

    ! assigning density 
    if (thisOctal%threeD)  then
       c0 =  discwind_density(this, x, y, z)    
       c1 =  discwind_density(this, x+d, y+d, z+d)
       c2 =  discwind_density(this, x+d, y+d, z-d)
       c3 =  discwind_density(this, x+d, y-d, z+d)
       c4 =  discwind_density(this, x+d, y-d, z-d)
       c5 =  discwind_density(this, x-d, y+d, z+d)
       c6 =  discwind_density(this, x-d, y+d, z-d)
       c7 =  discwind_density(this, x-d, y-d, z+d)
       c8 =  discwind_density(this, x-d, y-d, z-d)
       
       ave_discwind_density_max =  &
            (c0 + c1+ c2 + c3 + c4 + c5 + c6 + c7 + c8)/9.0d0
    else  ! 2D case
       c0 =  discwind_density(this, x, y, z)    
       c1 =  discwind_density(this, x+d, y, z+d)
       c2 =  discwind_density(this, x+d, y, z-d)
       c3 =  discwind_density(this, x-d, y, z+d)
       c4 =  discwind_density(this, x-d, y, z-d)
       
       ave_discwind_density_max =  max(c0, c1, c2, c3, c4)
    end if

  end function ave_discwind_density_max


    

  !
  ! Recursively turn off the inflow flag to false if grid > Rh.
  !
  RECURSIVE SUBROUTINE turn_off_discwind(thisOctal,grid, this)    
    
    IMPLICIT NONE
    
    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(discwind), INTENT(IN)  :: this
    
    TYPE(octal), POINTER   :: pChild
    
    INTEGER :: subcell, n
    
    if (thisOctal%threeD) then
       n = 8
    else
       n =4
    end if
    
    do subcell = 1, n
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_off_discwind(pChild,grid, this)
       else
          ! turnning it off 
          if (modulus(thisOctal%centre) > this%Rmin) thisOctal%inFlow(subcell) = .false.
       end if
    end do
    
  END SUBROUTINE turn_off_discwind


  !
  !
  ! Recursively turn off the inflow flag to false if grid > Rh.
  !
  RECURSIVE SUBROUTINE turn_on_discwind(thisOctal,grid, this)  
    
    IMPLICIT NONE

    TYPE(octal), POINTER   :: thisOctal
    TYPE(gridtype)         :: grid
    TYPE(discwind), INTENT(IN)  :: this

    TYPE(octal), POINTER   :: pChild
  
    INTEGER :: subcell, n

    if (thisOctal%threeD) then
       n = 8
    else
       n =4
    end if
     
    do subcell = 1, n
       if (thisOctal%hasChild(subcell)) then
          ! just decdend the tree branch
          pChild => thisOctal%child(subcell)
          CALL turn_on_discwind(pChild,grid, this)
       else
          ! turnning it on
          if (modulus(thisOctal%centre) > this%Rmin) thisOctal%inFlow(subcell) = .true.        
       end if
    end do
    
  END SUBROUTINE turn_on_discwind


  SUBROUTINE fill_Velocity_Corners(this,thisOctal, debug)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
    type(discwind) :: this
    TYPE(octal), INTENT(INOUT) :: thisOctal
    logical, optional :: debug
    logical :: writedebug
    real(oct)      :: r1, r2, r3
    real(oct)      :: phi1, phi2, phi3
    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    

    writedebug = .false.
    if (present(debug)) writedebug=debug

    if (thisOctal%oneD) then
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       y1 = 0.d0
       z1 = 0.d0
       thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(x1,y1,z1))
       thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(x2,y1,z1))
       thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(x3,y1,z1))
       goto 666
    endif


    if (thisOctal%threed) then
       if (.not.thisOctal%cylindrical) then ! 3d cartesian case
          ! we first store the values we use to assemble the position vectors
          
          x1 = thisOctal%centre%x - thisOctal%subcellSize
          x2 = thisOctal%centre%x
          x3 = thisOctal%centre%x + thisOctal%subcellSize
          
          y1 = thisOctal%centre%y - thisOctal%subcellSize
          y2 = thisOctal%centre%y
          y3 = thisOctal%centre%y + thisOctal%subcellSize
          
          z1 = thisOctal%centre%z - thisOctal%subcellSize
          z2 = thisOctal%centre%z
          z3 = thisOctal%centre%z + thisOctal%subcellSize
                    
           ! now store the 'base level' values
          
          thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(x1,y1,z1))
          thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(x2,y1,z1))
          thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(x3,y1,z1))
          thisOctal%cornerVelocity(4) = discwind_velocity(this,vector(x1,y2,z1))
          thisOctal%cornerVelocity(5) = discwind_velocity(this,vector(x2,y2,z1))
          thisOctal%cornerVelocity(6) = discwind_velocity(this,vector(x3,y2,z1))
          thisOctal%cornerVelocity(7) = discwind_velocity(this,vector(x1,y3,z1))
          thisOctal%cornerVelocity(8) = discwind_velocity(this,vector(x2,y3,z1))
          thisOctal%cornerVelocity(9) = discwind_velocity(this,vector(x3,y3,z1))
          
          ! middle level
          
          thisOctal%cornerVelocity(10) = discwind_velocity(this,vector(x1,y1,z2))
          thisOctal%cornerVelocity(11) = discwind_velocity(this,vector(x2,y1,z2))
          thisOctal%cornerVelocity(12) = discwind_velocity(this,vector(x3,y1,z2))
          thisOctal%cornerVelocity(13) = discwind_velocity(this,vector(x1,y2,z2))
          thisOctal%cornerVelocity(14) = discwind_velocity(this,vector(x2,y2,z2))
          thisOctal%cornerVelocity(15) = discwind_velocity(this,vector(x3,y2,z2))
          thisOctal%cornerVelocity(16) = discwind_velocity(this,vector(x1,y3,z2))
          thisOctal%cornerVelocity(17) = discwind_velocity(this,vector(x2,y3,z2))
          thisOctal%cornerVelocity(18) = discwind_velocity(this,vector(x3,y3,z2))
          
          ! top level
          
          thisOctal%cornerVelocity(19) = discwind_velocity(this,vector(x1,y1,z3))
          thisOctal%cornerVelocity(20) = discwind_velocity(this,vector(x2,y1,z3))
          thisOctal%cornerVelocity(21) = discwind_velocity(this,vector(x3,y1,z3))
          thisOctal%cornerVelocity(22) = discwind_velocity(this,vector(x1,y2,z3))
          thisOctal%cornerVelocity(23) = discwind_velocity(this,vector(x2,y2,z3))
          thisOctal%cornerVelocity(24) = discwind_velocity(this,vector(x3,y2,z3))
          thisOctal%cornerVelocity(25) = discwind_velocity(this,vector(x1,y3,z3))
          thisOctal%cornerVelocity(26) = discwind_velocity(this,vector(x2,y3,z3))
          thisOctal%cornerVelocity(27) = discwind_velocity(this,vector(x3,y3,z3))

       else ! cylindrical 
          if (thisOctal%splitAzimuthally) then
             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phiMin
             phi2 = thisOctal%phi 
             phi3 = thisOctal%phiMax
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize

             ! bottom level

             thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z1))
             thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(r3*cos(phi1),r3*sin(phi1),z1))
             thisOctal%cornerVelocity(4) = discwind_velocity(this,vector(r1*cos(phi2),r1*sin(phi2),z1))
             thisOctal%cornerVelocity(5) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z1))
             thisOctal%cornerVelocity(6) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z1))
             thisOctal%cornerVelocity(7) = discwind_velocity(this,vector(r1*cos(phi3),r1*sin(phi3),z1))
             thisOctal%cornerVelocity(8) = discwind_velocity(this,vector(r2*cos(phi3),r2*sin(phi3),z1))
             thisOctal%cornerVelocity(9) = discwind_velocity(this,vector(r3*cos(phi3),r3*sin(phi3),z1))

             thisOctal%cornerVelocity(10) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(11) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z2))
             thisOctal%cornerVelocity(12) = discwind_velocity(this,vector(r3*cos(phi1),r3*sin(phi1),z2))
             thisOctal%cornerVelocity(13) = discwind_velocity(this,vector(r1*cos(phi2),r1*sin(phi2),z2))
             thisOctal%cornerVelocity(14) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z2))
             thisOctal%cornerVelocity(15) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z2))
             thisOctal%cornerVelocity(16) = discwind_velocity(this,vector(r1*cos(phi3),r1*sin(phi3),z2))
             thisOctal%cornerVelocity(17) = discwind_velocity(this,vector(r2*cos(phi3),r2*sin(phi3),z2))
             thisOctal%cornerVelocity(18) = discwind_velocity(this,vector(r3*cos(phi3),r3*sin(phi3),z2))

             thisOctal%cornerVelocity(19) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(20) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z3))
             thisOctal%cornerVelocity(21) = discwind_velocity(this,vector(r3*cos(phi1),r3*sin(phi1),z3))
             thisOctal%cornerVelocity(22) = discwind_velocity(this,vector(r1*cos(phi2),r1*sin(phi2),z3))
             thisOctal%cornerVelocity(23) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z3))
             thisOctal%cornerVelocity(24) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z3))
             thisOctal%cornerVelocity(25) = discwind_velocity(this,vector(r1*cos(phi3),r1*sin(phi3),z3))
             thisOctal%cornerVelocity(26) = discwind_velocity(this,vector(r2*cos(phi3),r2*sin(phi3),z3))
             thisOctal%cornerVelocity(27) = discwind_velocity(this,vector(r3*cos(phi3),r3*sin(phi3),z3))


          else

             z1 = thisOctal%centre%z - thisOctal%subcellSize
             z2 = thisOctal%centre%z
             z3 = thisOctal%centre%z + thisOctal%subcellSize
             phi1 = thisOctal%phi - thisOctal%dPhi/2.d0
             phi2 = thisOctal%phi + thisOctal%dPhi/2.d0
             r1 = thisOctal%r - thisOctal%subcellSize
             r2 = thisOctal%r
             r3 = thisOctal%r + thisOctal%subcellSize


             ! bottom level

             thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z1))
             thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z1))
             thisOctal%cornerVelocity(4) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z1))
             thisOctal%cornerVelocity(5) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z1))
             thisOctal%cornerVelocity(6) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z1))

             ! middle level

             thisOctal%cornerVelocity(7) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(8) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z2))
             thisOctal%cornerVelocity(9) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z2))
             thisOctal%cornerVelocity(10) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z2))
             thisOctal%cornerVelocity(11) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z2))
             thisOctal%cornerVelocity(12) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z2))

             ! top level

             thisOctal%cornerVelocity(13) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(14) = discwind_velocity(this,vector(r1*cos(phi1),r1*sin(phi1),z3))
             thisOctal%cornerVelocity(15) = discwind_velocity(this,vector(r2*cos(phi1),r2*sin(phi1),z3))
             thisOctal%cornerVelocity(16) = discwind_velocity(this,vector(r2*cos(phi2),r2*sin(phi2),z3))
             thisOctal%cornerVelocity(17) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z3))
             thisOctal%cornerVelocity(18) = discwind_velocity(this,vector(r3*cos(phi2),r3*sin(phi2),z3))

          endif
       endif
    else       
       
    ! we first store the values we use to assemble the position vectors
       
       x1 = thisOctal%centre%x - thisOctal%subcellSize
       x2 = thisOctal%centre%x
       x3 = thisOctal%centre%x + thisOctal%subcellSize
       
       z1 = thisOctal%centre%z - thisOctal%subcellSize
       z2 = thisOctal%centre%z
       z3 = thisOctal%centre%z + thisOctal%subcellSize
       
       ! now store the 'base level' values
       
       thisOctal%cornerVelocity(1) = discwind_velocity(this,vector(x1,0.d0,z1))
       thisOctal%cornerVelocity(2) = discwind_velocity(this,vector(x2,0.d0,z1))
       thisOctal%cornerVelocity(3) = discwind_velocity(this,vector(x3,0.d0,z1))
       thisOctal%cornerVelocity(4) = discwind_velocity(this,vector(x1,0.d0,z2))
       thisOctal%cornerVelocity(5) = discwind_velocity(this,vector(x2,0.d0,z2))
       thisOctal%cornerVelocity(6) = discwind_velocity(this,vector(x3,0.d0,z2))
       thisOctal%cornerVelocity(7) = discwind_velocity(this,vector(x1,0.d0,z3))
       thisOctal%cornerVelocity(8) = discwind_velocity(this,vector(x2,0.d0,z3))
       thisOctal%cornerVelocity(9) = discwind_velocity(this,vector(x3,0.d0,z3))
    endif
666 continue

!    if(isnan(thisOctal%cornerVelocity(1)%x)) then
!          write(*,*) "cornervel",thisOctal%cornerVelocity(1)
!          write(*,*) discwind_velocity(this,vector(x1,0.d0,z1),grid)
!          write(*,*) x1,z1
!          write(*,*) x2,z2
!          write(*,*) x3,z3
!       enddo
!    endif
    
  END SUBROUTINE fill_Velocity_Corners

  SUBROUTINE add_new_children_discwind(parent, ichild, grid, this, splitAzimuthally)
    use memory_mod
    use inputs_mod, only : cylindrical
    use amr_mod, only : growChildArray, allocateOctalAttributes, setSmallestSubcell
    ! adds all eight new children to an octal
    IMPLICIT NONE
    logical, optional :: splitAzimuthally
    type(VECTOR) :: rVec
    TYPE(octal), POINTER :: parent     ! pointer to the parent octal 
    TYPE(gridtype), INTENT(INOUT) :: grid ! need to pass the grid to routines that
                                          !   calculate the variables stored in
                                          !   the tree.
    TYPE(discwind), INTENT(IN)  :: this
    INTEGER       :: subcell           ! loop counter
                                       ! - this isn't very clever. might change it. 
    INTEGER       :: newChildIndex     ! the storage location for the new child
    INTEGER       :: iChild
    integer :: parentSubcell
    real(double) :: x, y, z
    integer :: nChildren
    TYPE(OCTAL), POINTER :: thisOctal

    nChildren = parent%nChildren

    parentSubcell = iChild

    ! safety checks of child array
    IF ( ASSOCIATED(parent%child) ) THEN
      IF ( ( nChildren == 0 ) .OR.                  &
           ( nChildren /= SIZE(parent%child) ) ) THEN
        PRINT *, 'Panic: in addNewChild, %child array wrong size'
        PRINT *, 'nChildren:',nChildren,' SIZE %child:', SIZE(parent%child)
        STOP
      END IF
    END IF
    IF ( (.NOT. ASSOCIATED(parent%child)) .AND. (nChildren > 0) ) THEN
      PRINT *, 'Panic: in addNewChild, %child array wrong size'
      PRINT *, 'nChildren:',nChildren,' ASSOCIATED %child:', ASSOCIATED(parent%child)
      STOP
    END IF

    ! check that new child does not already exist
    IF ( parent%hasChild(iChild) .EQV. .TRUE. ) THEN
      PRINT *, 'Panic: in addNewChild, attempted to add a child ',&
               '       that already exists'
      STOP
    ENDIF

!    call interpFromParent(subcellCentre(parent, iChild, parent%subcellSize, &
!         grid, temperature, density, dusttypeFraction)

    CALL growChildArray(parent, nNewChildren=1, grid=grid )

    ! update the bookkeeping
    nChildren = nChildren + 1
    newChildIndex = nChildren

    ! update the parent octal
    parent%nChildren = nChildren
    parent%hasChild(iChild) = .TRUE.
    parent%indexChild(newChildIndex) = iChild

    NULLIFY(parent%child(newChildIndex)%child)

    parent%child(newChildIndex)%nDepth = parent%nDepth + 1
    if (parent%child(newChildIndex)%nDepth  > grid%maxDepth) then
       grid%maxDepth = grid%maxDepth + 1
       CALL setSmallestSubcell(grid)
    endif
    ! set up the new child's variables
    parent%child(newChildIndex)%threeD = parent%threeD
    parent%child(newChildIndex)%twoD = parent%twoD
    parent%child(newChildIndex)%oneD = parent%oneD
    parent%child(newChildIndex)%maxChildren = parent%maxChildren
    parent%child(newChildIndex)%cylindrical = parent%cylindrical

    if (cylindrical) then  
       if (parent%splitAzimuthally) then
          rVec =  subcellCentre(parent,iChild)
          parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
          if (parent%child(newChildIndex)%phi < 0.d0) then
              parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
          endif
          parent%child(newChildIndex)%dphi = parent%dphi/2.d0
          parent%child(newChildIndex)%phimin = parent%child(newChildIndex)%phi - parent%dPhi/4.d0
          parent%child(newChildIndex)%phimax = parent%child(newChildIndex)%phi + parent%dPhi/4.d0
       else
          parent%child(newChildIndex)%phi = parent%phi
          parent%child(newChildIndex)%dphi = parent%dphi
          parent%child(newChildIndex)%phimin = parent%phimin
          parent%child(newChildIndex)%phimax = parent%phimax
       endif
       if (parent%child(newChildIndex)%phimin < 1.d-10) parent%child(newChildIndex)%phimin = 0.d0 ! fixed!!!!
       parent%child(newChildIndex)%splitAzimuthally = .false.
       parent%child(newChildIndex)%maxChildren = 4

       if (PRESENT(splitAzimuthally)) then
          if (splitAzimuthally) then
             parent%child(newChildIndex)%splitAzimuthally = .true.
             parent%child(newChildIndex)%maxChildren = 8
             rVec =  subcellCentre(parent,iChild)
             parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
             if (parent%child(newChildIndex)%phi < 0.d0) then
                parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
             endif
             if (parent%splitAzimuthally) then
                parent%child(newChildIndex)%dphi = parent%dphi/2.d0
             else
                parent%child(newChildIndex)%dphi = parent%dphi
             endif
          else
             if (parent%splitAzimuthally) then
                rVec =  subcellCentre(parent,iChild)
                parent%child(newChildIndex)%phi = atan2(rvec%y, rVec%x)
                if (parent%child(newChildIndex)%phi < 0.d0) then
                   parent%child(newChildIndex)%phi = parent%child(newChildIndex)%phi + twoPi 
                endif
                parent%child(newChildIndex)%dphi = parent%dphi/2.d0
             else
                parent%child(newChildIndex)%phi = parent%phi
                parent%child(newChildIndex)%dphi = parent%dphi
             endif
             parent%child(newChildIndex)%splitAzimuthally = .false.
             parent%child(newChildIndex)%maxChildren = 4
          endif
       endif
    endif
    

    parent%child(newChildIndex)%inFlow = parent%inFlow
    parent%child(newChildIndex)%parent => parent
    parent%child(newChildIndex)%parentSubcell = iChild
    parent%child(newChildIndex)%subcellSize = parent%subcellSize / 2.0_oc
    parent%child(newChildIndex)%hasChild = .false.
    parent%child(newChildIndex)%nChildren = 0
    parent%child(newChildIndex)%indexChild = -999 ! values are undefined
    parent%child(newChildIndex)%centre = subcellCentre(parent,iChild)
    if (parent%cylindrical) then
       parent%child(newChildIndex)%r = subcellRadius(parent,iChild)
    endif


    parent%child(newChildIndex)%xMin = parent%child(newChildIndex)%centre%x - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMin = parent%child(newChildIndex)%centre%y - parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMin = parent%child(newChildIndex)%centre%z - parent%child(newChildIndex)%subcellSize

    parent%child(newChildIndex)%xMax = parent%child(newChildIndex)%centre%x + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%yMax = parent%child(newChildIndex)%centre%y + parent%child(newChildIndex)%subcellSize
    parent%child(newChildIndex)%zMax = parent%child(newChildIndex)%centre%z + parent%child(newChildIndex)%subcellSize


    thisOctal => parent%child(newChildIndex)
    call allocateOctalAttributes(grid, thisOctal)

    globalMemoryFootprint = globalMemoryFootprint + octalMemory(thisOctal)

    call fill_velocity_corners(this, thisOctal)


    do subcell = 1, thisOctal%maxChildren
       thisOctal%temperature(subcell) = parent%temperature(parentSubcell)
       thisOctal%rho(subcell) = parent%rho(parentSubcell)
       thisOctal%velocity(subcell) = parent%velocity(parentSubcell)

       rVec = subcellCentre(thisOctal, subcell)
       x = rVec%x
       y = rVec%y
       z = rVec%z
       thisOctal%inFlow(subcell) = .true. !all_in_discwind(thisOctal, subcell, this)
       if (thisOctal%inflow(subcell)) then
          thisOctal%temperature(subcell) = real(this%Twind)
          thisOctal%rho(subcell) = ave_discwind_density(thisOctal, subcell, this)
          thisOctal%velocity(subcell) = discwind_velocity(this, vector(x,y,z))
       endif
    enddo
  end SUBROUTINE add_new_children_discwind





end module discwind_class
