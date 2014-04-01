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
       &   get_parameters, &
       &   component, &
       &   in_discwind, &
       &   all_in_discwind, &
       &   discwind_density, &
       &   discwind_Vr, &
       &   discwind_Vphi, &
       &   discwind_velocity, &
       &   ave_discwind_density, &
       &   turn_on_discwind, &
       &   turn_off_discwind

  private::  &
       &   int_discwind, &
       &   need_to_split, &
       &   speed_of_sound, &
       &   escape_velocity,&
       &   temperature_disc, &
       &   mdot_on_disc,&
       &   ave_discwind_density_slow, &
       &   ave_discwind_density_fast


  !----------------------------------------------------------
  type discwind_type
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

  end type discwind_type
  


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

  type(DISCWIND_TYPE) :: globalDiscWind

contains
  ! =====================================================================
  ! constructors
  !======================================================================

  ! initializing the discwind with default parameters
  subroutine int_discwind_default(this)
    implicit none
    type(discwind_type), intent(inout) :: this 
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
    
    type(discwind_type), intent(inout) :: this 
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

! Set up the globalDiscWind stored in this module using values from the parameter file
  subroutine addGlobalDiscWind
!#    use inputs_mod, only : DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
!       DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, limitscalar, ttauriMstar

    use inputs_mod, only : DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
       DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, ttauriMstar
    real(double) :: DW_Hdisc

    call new(globalDiscWind, DW_d, DW_Rmin, DW_Rmax, DW_Tmax, DW_gamma, &
       DW_Mdot, DW_alpha, DW_beta, DW_Rs, DW_f, DW_Twind, dble(ttauriMstar)/msol, DW_Hdisc)

  end subroutine addGlobalDiscWind


  !==========================================================================
  ! Accessors
  !==========================================================================

  ! Given a discwind object, it returns all components of the objects.
  subroutine get_parameters(this, d, Rmin, Rmax, Tmax, gamma, &
       Mdot, alpha, beta, Rs, f, Twind, Mstar, Hdisc)
    implicit none 
    
    type(discwind_type), intent(in) :: this
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
    
    type(discwind_type), intent(in) :: this 
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
    type(discwind_type), intent(in) :: this
    real(double), intent(in) :: R ! [10^10cm] of distance from z-axis
    
    out = this%Tmax * (r/this%Rmin)**this%gamma  ! [K]

    if (out < 3.0d0) out = 3.0d0  ! setting the minimum temperature

  end function temperature_disc


  !
  ! Approximated speed of sound in [cm/s]
  ! (private)
  real(double) function speed_of_sound(this, r)  
    implicit none
    type(discwind_type), intent(in) :: this
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
    type(discwind_type), intent(in) :: this
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
    type(discwind_type), intent(in) :: this
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
    type(discwind_type), intent(in) :: this
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
    type(discwind_type), intent(in) :: this
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
       if (phi < 0.d0) phi = phi + twoPi
       theta = ACOS( (ABS(z)+d) / rp )  ! origin shifted
       
       Vx = Vr*SIN(theta)*COS(phi) - Vphi*SIN(phi)      ! [c]
       Vy = -1.d0*(Vr*SIN(theta)*SIN(phi) + Vphi*COS(phi))      ! [c]
       Vz = Vr*COS(theta)                               ! [c]
       if (z < 0.d0) Vz = -Vz

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
    type(discwind_type), intent(in) :: this
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
    type(DISCWIND_TYPE) :: this
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
    type(discwind_type), intent(in) :: this
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
    type(discwind_type), intent(in) :: this
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
  !
  ! This is based decideSplit in amr_mod.f90
  logical function need_to_split(thisOctal,subcell, this, limitscalar)
    
    IMPLICIT NONE
    
    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(discwind_type), intent(in) :: this
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

  

  real(double) function ave_discwind_density_slow(thisOctal,subcell, this)

    IMPLICIT NONE

    TYPE(octal), POINTER       :: thisOctal
    INTEGER, INTENT(IN)        :: subcell
    type(discwind_type), intent(in) :: this
    
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
    type(discwind_type), intent(in) :: this
    
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
    type(discwind_type), intent(in) :: this
    
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
    TYPE(discwind_type), INTENT(IN)  :: this
    
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
    TYPE(discwind_type), INTENT(IN)  :: this

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




  recursive subroutine assignDensitiesDiscwind(grid, thisOctal, thisWind)
    use inputs_mod, only : vturb
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, child
    TYPE(discwind_type), INTENT(IN)  :: thisWind
    type(VECTOR) :: rVec
    real(double) :: x, y, z
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child =>thisOctal%child(i)
                call assignDensitiesDiscwind(grid, child, thisWind)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal, subcell)
          x = rVec%x
          y = rVec%y
          z = rVec%z
          if (all_in_discwind(thisOctal, subcell, thisWind)) then
             thisOctal%inFlow(subcell) = .true.
             thisOctal%iAnalyticalVelocity(subcell) = 1
             thisOctal%temperature(subcell) = real(thisWind%Twind)
             thisOctal%rho(subcell) = ave_discwind_density(thisOctal, subcell, thisWind)
             thisOctal%velocity(subcell) = discwind_velocity(thisWind, vector(x,y,z))
             thisOctal%fixedTemperature(subcell) = .true.
             if (associated(thisOctal%microturb)) thisOctal%microturb(subcell) = vturb
          endif
       endif
    enddo
  end subroutine assignDensitiesDiscwind


end module discwind_class
