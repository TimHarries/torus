module jets_mod
  ! Current version handle the following geometry case:
  !
  !  1. bipolar jets
  !  
  ! --------------------------------------------------------------------------------
  !  Created: OCT-03-2002  (R.Kurosawa)
  !
  !---------------------------------------------------------------------------------

  USE vector_mod           ! vector maths routines
  USE kind_mod               ! variable kind parameters
  USE octal_mod              ! type definition for the grid elements
  USE gridtype_mod         ! type definition for the 3-d grid
  USE constants_mod

    
  public :: &
       set_jets_parameters, &
       get_jets_parameter, &
       JetsDensity, &
       JetsVelocity,&
       JetsVelocity_in_C,&
       JetsTemperature,&
       calcJetsMassVelocity, &
       calcJetsTemperature, &
       initJetsAMR, &
       dV_dn_jets

  private :: &
       VelocityCorners
  
  
  !--------------------------------------------------------------------------------
  ! Bipolar Jet flow
  ! Flow into cones with openning angle theta_o around
  !  z-axis as the symetry axis.
  !
  !  The velocity is assumed to follow the classical
  !  beta velocity law with an offset:
  !
  !  In this formulation:
  !
  !                    Mdot
  !  rho(r) = ------------------------------
  !            4Pi V(r) r^2 (1-Cos(theta_o))
  !
  !  and
  !                         r*
  !   V(r) =  Vinf * ( 1 - --- )^beta + Vo
  !                         r
  !
  !   where Vinf is the terminal speed and
  !         Vo is the offeset speed.
  ! theta_o  is "half" opening angle.
  !
  ! For temperature structure of the bipolar is assumed as:
  !                Rmin
  !    T = Tcore*(------)^e6
  !                  R
  !
  !===============================================================
  type jets_parameters
     private
     real :: Rmin    ! radius of central star  [10^10cm]
     real :: Rmax    ! cutoff radius in [10^10 cm]
     real :: Vo      ! a small offset in V in  [km/s]
     
     ! For jets
     real :: Vinf_jets    ! terminal velocity in  [kms]
     real :: beta_jets    ! beta in the beta-belocity law [-]
     real :: Mdot_jets    ! mass loss rate in jets [M_solar/yr]
     real :: theta_o_jets ! half opening angle [degrees]     
     real :: Tcore_jets   ! temperature at the core at Rmin in [10^4 K]
     real :: e6_jets      ! an exponet of in tempereture eq. [-]
     
     ! For disk wind
     real :: Vinf_disk    ! terminal velocity in  [kms]
     real :: beta_disk    ! beta in the beta-belocity law [-]
     real :: Mdot_disk    ! mass loss rate in jets [M_solar/yr]
     real :: theta_o_disk ! half opening angle [degrees]     
     real :: Tcore_disk   ! temperature at the core at Rmin in [10^4 K]
     real :: e6_disk      ! an exponet of in tempereture eq. [-]

     ! For equatorial disk
     real :: Rdisk_min  ! The minimum radius of the disk. in 10^10 cm.
     real :: Rdisk_max  ! The minimum  of the disk in 10^10 cm.
     real :: h_disk     ! Thickness of the disk in 10^10 cm.
     real :: rho_scale  ! The density in the units of rho max for disk wind.     
     
  end type jets_parameters

  
  type(jets_parameters), private :: jets  ! One instance to be used ONLY IN THIS MODULE.

  !
  !
  !
  
contains

  !
  ! set the jets_parameters.
  !
  subroutine set_jets_parameters(Rmin, Rmax, Vo, &
       Vinf_jets, beta_jets, Mdot_jets, theta_o_jets, Tcore_jets, e6_jets, & 
       Vinf_disk, beta_disk, Mdot_disk, theta_o_disk, Tcore_disk, e6_disk, &
       & Rdisk_min, Rdisk_max, h_disk, rho_scale)
    
    implicit none
    real,intent(in) :: Rmin    ! radius of central star  [10^10cm]
    real,intent(in) :: Rmax    ! cutoff radius in [10^10 cm]
    real,intent(in) :: Vo      ! a small offset in V in  [km/s]

    ! For jets
    real, intent(in)  :: Vinf_jets    ! terminal velocity in  [kms]
    real, intent(in)  :: beta_jets    ! beta in the beta-belocity law [-]
    real, intent(in)  :: Mdot_jets    ! mass loss rate in jets [M_solar/yr]
    real, intent(in)  :: theta_o_jets ! half opening angle [degrees]     
    real, intent(in)  :: Tcore_jets   ! temperature at the core at Rmin in [10^4 K]
    real, intent(in)  :: e6_jets      ! an exponet of in tempereture eq. [-]
    
    ! For disk wind
    real, intent(in)  :: Vinf_disk    ! terminal velocity in  [kms]
    real, intent(in)  :: beta_disk    ! beta in the beta-belocity law [-]
    real, intent(in)  :: Mdot_disk    ! mass loss rate in jets [M_solar/yr]
    real, intent(in)  :: theta_o_disk ! half opening angle [degrees]     
    real, intent(in)  :: Tcore_disk   ! temperature at the core at Rmin in [10^4 K]
    real, intent(in)  :: e6_disk      ! an exponet of in tempereture eq. [-]
    
    ! For equatorial disk
    real, intent(in)  :: Rdisk_min  ! The minimum radius of the disk. in 10^10 cm.
    real, intent(in)  :: Rdisk_max  ! The minimum  of the disk in 10^10 cm.
    real, intent(in)  :: h_disk     ! Thickness of the disk in 10^10 cm.
    real, intent(in)  :: rho_scale  ! The density in the units of rho max for disk wind.     
    
    jets%Rmin     =     Rmin   
    jets%Rmax     =     Rmax   
    jets%Vo       =     Vo

    jets%Vinf_jets     =     Vinf_jets 
    jets%beta_jets     =     beta_jets   
    jets%Mdot_jets     =     Mdot_jets
    jets%theta_o_jets  =     theta_o_jets
    jets%Tcore_jets    =     Tcore_jets
    jets%e6_jets       =     e6_jets

    jets%Vinf_disk     =     Vinf_disk
    jets%beta_disk     =     beta_disk
    jets%Mdot_disk     =     Mdot_disk
    jets%theta_o_disk  =     theta_o_disk
    jets%Tcore_disk    =     Tcore_disk
    jets%e6_disk       =     e6_disk

    jets%Rdisk_min =  Rdisk_min 
    jets%Rdisk_max =  Rdisk_max 
    jets%h_disk    =  h_disk    
    jets%rho_scale =  rho_scale 

  end subroutine set_jets_parameters


  !
  ! access a parameter in jets_parameters
  !
  function get_jets_parameter(name)  RESULT(out)
    implicit none
    real :: out
    character(LEN=*) :: name
    
    select case(name)
    case ("Rmin")
       out = jets%Rmin
    case ("Rmax")
       out = jets%Rmax
    case ("Vo_jets")
       out = jets%Vo
    case ("Vinf_jets")
       out = jets%Vinf_jets
    case ("beta_jets")
       out = jets%beta_jets
    case ("Mdot_jets")
       out = jets%Mdot_jets
    case ("theta_o_jets")
       out = jets%theta_o_jets
    case ("Tcore_jets")
       out = jets%Tcore_jets
    case ("e6_jets")
       out = jets%e6_jets
    case ("Vinf_disk")
       out = jets%Vinf_disk
    case ("beta_disk")
       out = jets%beta_disk
    case ("Mdot_disk")
       out = jets%Mdot_disk
    case ("theta_o_disk")
       out = jets%theta_o_disk
    case ("Tcore_disk")
       out = jets%Tcore_disk
    case ("e6_disk")
       out = jets%e6_disk
    case default
       print *, 'Unknown name passed rho_vel_temp::get_jets_parameters.'
       print *, 'name =', TRIM(name)
       print *, 'Exiting the program ...'
       stop
    end select

  end function get_jets_parameter

  
  !
  !
  !    
  FUNCTION JetsDensity(point,grid) RESULT(out)
    ! Given a point and grid, this function returns
    ! the density at this point in [kg/m^3].
    ! 
    !  Jets flow. Flow in to cones with openning angle theta_o
    !  The velocity is assumed to follow the classical
    !  beta velocity law with an offset:
    !
    !  In this formulation:
    !
    !                     Mdot
    !  rho(r) = ------------------------------
    !            4Pi V(r) r^2 (1-Cos(theta_o))
    !
    !  and
    !                         r*
    !   V(r) =  Vinf * ( 1 - --- )^beta + Vo
    !                         r
    !
    !   where Vinf is the terminal speed and
    !         Vo is the offeset speed.
    ! theta_o  is "half" opening angle
    !
    
    !    USE inputs_mod  ! inputs needed for this function are defined in here.
    USE constants_mod
    
    IMPLICIT NONE

    REAL                          :: out
    TYPE(vector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec

    !
    LOGICAL, SAVE  :: first_time = .true. 
    real(double) ::  r, rn  ! distance from the center of central star
    real(double) :: rp, rho_jets, rho_disk
    real(double) :: vr
    real(double) :: cos_theta
    real(double) :: theta_o_jets, theta_o_disk
    !
    real(double), SAVE :: cos_theta_o_jets    
    real(double), SAVE :: Mdot_jets             ! [kg/sec]
    real(double), SAVE :: cos_theta_o_disk    
    real(double), SAVE :: Mdot_disk             ! [kg/sec]
    !
    real(double), parameter :: M_sun=1.9891d30 ! in [kg]

    real(double), save :: rho_max_disk_wind, rho_min_disk_wind

    IF (first_time) THEN
       ! Converting the half open angle to radians
       theta_o_jets = jets%theta_o_jets * Pi/180.0d0 
       cos_theta_o_jets = COS(theta_o_jets)
       theta_o_disk = jets%theta_o_disk * Pi/180.0d0 
       cos_theta_o_disk = COS(theta_o_disk)

       ! check the validity if theta_o_jets
       if (theta_o_jets > Pi/2.0d0 .or. theta_o_disk > Pi/2.0d0) then
          write(*,*) "Error:: theta_o_jets > Pi/2.0d0 .or. theta_o_disk > Pi/2.0d0."
          write(*,*) "        in [jets_mod::JetsDensity]."
          stop
       end if
       ! Converting Mdot to kg/s
       Mdot_jets = jets%Mdot_jets*M_sun/(60.0d0*60.0d0*24.0d0*365.0d0) ! [kg/sec]
       Mdot_disk = jets%Mdot_disk*M_sun/(60.0d0*60.0d0*24.0d0*365.0d0) ! [kg/sec]

       ! Computing the maximum desnity of the disk wind
       Rp = 1.0d0    
       Vr = jets%Vo     ! [km/s]
       ! density in [kg/m^3]     
       !-- convert r in [m] was in [10^10cm]
       r = jets%Rmin  * 1.0d8
       ! convert [km/s] to [m/s]
       Vr = Vr*10.0d3  ! m/s    
       
       ! this value is save for the successive calls
       rho_max_disk_wind  =       Mdot_disk &
            &                         /                 &
            ( 4.0*Pi*r*r*Vr*(1.0-cos_theta_o_disk) )


       ! this value is save for the successive calls
       Rp = jets%Rmin/jets%Rmax
       Vr = jets%Vinf_disk * ( (1.0 - Rp)**jets%beta_disk ) + jets%Vo     ! [km/s]
       ! convert [km/s] to [m/s]
       Vr = Vr*10.0d3  ! m/s           
       r = jets%Rmax  * 1.0d8
       rho_min_disk_wind  =       Mdot_disk &
            &                         /                 &
            ( 4.0*Pi*r*r*Vr*(1.0-cos_theta_o_disk) )
  

       ! Quick check
       if (theta_o_disk < theta_o_jets) then
          print *, 'Error:: theta_o_disk must be greater than theta_o_jets.'
          print *, '  Exiting the program ...'
          stop
       end if


       first_time = .false.
    END IF
    !
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn)
    r = modulus( pointVec ) 
    !-- convert r in [m] was in [10^10cm]
    rn = r * 1.0d8
        
    
    IF (pointVec%z /= 0.0) THEN
       cos_theta = pointVec%z/r
    ELSE
       cos_theta = 0.0d0
    END IF



    if ( r >= jets%Rmin ) then
       Rp = jets%Rmin/r
    else
       Rp = 1.0d0
    end if
    
    !
    ! computing the density contribution from the jets
    !
    Vr =  jets%Vinf_jets * (1.0d0 - Rp)**jets%beta_jets  + jets%Vo     ! [km/s]
    ! density in [kg/m^3]     
    
    ! convert [km/s] to [m/s]
    Vr = Vr*10.0d3  ! m/s
    
    rho_jets  =       Mdot_jets               &
         /                 &
         ( 4.0d0*Pi*rn*rn*Vr*(1.0d0-cos_theta_o_jets) )  ! [kg/m^3]


    !
    ! Contribution from the disk wind (outer cones)
    !

    Vr =  jets%Vinf_disk * (1.0d0 - Rp)**jets%beta_disk  + jets%Vo     ! [km/s]
    ! density in [kg/m^3]     
    
    ! convert [km/s] to [m/s]
    Vr = Vr*10.0d3  ! m/s
    
    rho_disk  =       Mdot_disk               &
         /                 &
         ( 4.0d0*Pi*rn*rn*Vr*(1.0d0-cos_theta_o_disk) )  ! [kg/m^3]
        


    !
!    if (r >= jets%Rmin .and. r < jets%Rmax &
    if (r >= jets%Rdisk_min .and. r < jets%Rmax &
       .and. abs(cos_theta) > cos_theta_o_jets) then

       out = rho_jets + rho_disk

    elseif (r >= jets%Rdisk_min .and. r < jets%Rmax &
         .and. abs(cos_theta) > cos_theta_o_disk .and. &
         abs(cos_theta) < cos_theta_o_jets ) then

       out = rho_disk

    else if (r >= jets%Rmin .and.  r < jets%Rdisk_min) then

       ! Spherical wind up to the disk radius
       out = rho_disk

    elseif (abs(pointVec%z)<jets%h_disk/2.0 .and.  &
         (r>jets%Rdisk_min .and. r<jets%Rdisk_max)) then

       out = jets%rho_scale*rho_max_disk_wind

    else
       
       out = rho_min_disk_wind

    end if
                   
  END FUNCTION JetsDensity

  !===================================================================
  ! Given a point and grid, this function returns
  ! the radial component of the jed velocity (other components at this point
  ! are assumed to be zero for now).  The returned values is in [km/s]
  ! 
  !  Jets flow. Flow into the cones with openning angle theta_o.
  !  The symmetry axis is z.
  !
  !                         r*
  !   V(r) =  Vinf * ( 1 - --- )^beta + Vo  if within the cones
  !                         r
  !
  !        = 0  otherwise
  !
  !   where Vinf is the terminal speed and
  !         Vo is the offeset speed.
  !  NB: theta_o  is "half" opening angle
  !====================================================================
  FUNCTION JetsVelocity(point,grid) RESULT(out)
!    USE inputs_mod  ! inputs needed for this function are defined in here.
    USE constants_mod
    
    IMPLICIT NONE

    REAL                          :: out
    TYPE(vector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec

    !
    LOGICAL, SAVE  :: first_time = .true. 
    real(double) ::  r  ! distance from the center of central star
    real(double) :: Vr
    real(double) :: Rp, Vinf, beta
    real(double) :: cos_theta 
    real(double) :: theta_o_jets, theta_o_disk
    !
    real(double), SAVE :: cos_theta_o_jets
    real(double), SAVE :: cos_theta_o_disk    

    IF (first_time) THEN
       ! Converting the half open angle to radians
       theta_o_jets = jets%theta_o_jets * Pi/180.0d0 
       cos_theta_o_jets = COS(theta_o_jets)

       theta_o_disk = jets%theta_o_disk * Pi/180.0d0 
       cos_theta_o_disk = COS(theta_o_disk)

       ! Quick checks
       if (theta_o_jets > Pi/2.0d0 .or. theta_o_disk > Pi/2.0d0) then
          write(*,*) "Error:: theta_o_jets > Pi/2.0d0 .or. theta_o_disk > Pi/2.0d0."
          write(*,*) "        in [jets_mod::JetVelocity]."
          stop
       end if

       if (theta_o_disk < theta_o_jets) then
          print *, 'Error:: theta_o_disk must be greater than theta_o_jets.'
          print *, ' Exiting the program ...'
          stop
       end if
       first_time = .false.
    END IF
    !
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn)
    r = modulus( pointVec )   
    
    
    if (pointVec%z /= 0.0d0)then
       cos_theta = pointVec%z/r
    else
       cos_theta = 0.0d0
    end if
    !

    !
    if (r >= jets%Rmin .and. r < jets%Rmax &
         .and. abs(cos_theta) > cos_theta_o_jets) then
       Vinf = jets%Vinf_jets
       beta = jets%beta_jets
    elseif (r >= jets%Rmin .and. r < jets%Rmax  &
         .and. abs(cos_theta) > cos_theta_o_disk .and. &
         abs(cos_theta) < cos_theta_o_jets ) then
       Vinf = jets%Vinf_disk
       beta = jets%beta_disk
    else if (r >= jets%Rmin .and.  r < jets%Rdisk_min .and. &
         &  abs(cos_theta) < cos_theta_o_disk) then
       ! Spherical wind up to the disk radius
       Vinf = jets%Vinf_disk
       beta = jets%beta_disk
    else
       Vinf = jets%Vinf_disk
       beta = jets%beta_disk
       r = jets%Rmin + jets%Rmin/100.d0
    end if

    
    Rp = jets%Rmin/r
        
    
    Vr = Vinf * (1.0d0 - Rp)**beta  + jets%Vo  ! [km/s]

    out = Vr ! [km/s]

    
  END FUNCTION JetsVelocity



  ! Same as JetsVelocity function defined above except that
  ! this one returns a vector with it's components in the unit of speed of light. 
  FUNCTION JetsVelocity_in_C(point,grid) RESULT(out)
    USE constants_mod
    
    IMPLICIT NONE

    TYPE(vector)                  :: out  ! componets in [c]
    TYPE(vector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid
    real(double)                          :: Vr
    TYPE(vector)                  :: rhat

    ! radial speed in [c]
    Vr = JetsVelocity(point,grid)*1.0d5/cSpeed  
    
    rhat = (point)/modulus(point)
    
    out = Vr * rhat

  end FUNCTION JetsVelocity_in_C



  !
  !
  !
  !===============================================================
  !  Returns jet temperature in [10^4K].
  !
  ! For temperature structure of the jets is assumed as:
  !                Rmin
  !    T = Tcore*(------)^e6        [10^4 K]
  !                  R
  !
  !   Tcore  and e6 are to be assinged in the main parameter file.
  !
  !===============================================================  
  FUNCTION JetsTemperature(point,grid)  RESULT(out)
    
    IMPLICIT NONE
    
    REAL  ::  out 
    TYPE(vector), INTENT(IN) :: point 
    TYPE(gridtype), INTENT(IN) :: grid
    
    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    
    REAL :: Rp, cos_theta
    real :: r
    
    real(double) :: theta_o_jets, theta_o_disk
    real(double),save :: cos_theta_o_jets, cos_theta_o_disk
    logical, save :: first_time = .true.
    real(double) :: Tcore, e6

    starPosn = grid%starPos1
    pointVec = (point - starPosn)
    r = modulus( pointVec )



    IF (first_time) THEN
       ! Converting the half open angle to radians
       theta_o_jets = jets%theta_o_jets * Pi/180.0d0 
       cos_theta_o_jets = COS(theta_o_jets)

       theta_o_disk = jets%theta_o_disk * Pi/180.0d0 
       cos_theta_o_disk = COS(theta_o_disk)

       ! Quick check
       if (theta_o_jets > Pi/2.0d0 .or. theta_o_disk > Pi/2.0d0) then
          write(*,*) "Error:: theta_o_jets > Pi/2.0d0 .or. theta_o_disk > Pi/2.0d0."
          write(*,*) "        in [jets_mod::JetsTemperature]."
          stop
       end if
       if (theta_o_disk < theta_o_jets) then
          print *, 'Error:: theta_o_disk must be greater than theta_o_jets.'
          print *, ' Exiting the program ...'
          stop
       end if
       first_time = .false.
    END IF



    if (pointVec%z /= 0.0d0)then
       cos_theta = pointVec%z/r
    else
       cos_theta = 0.0d0
    end if
    !
    !
    if (r >= jets%Rmin .and. r < jets%Rmax &
         .and. abs(cos_theta) > cos_theta_o_jets) then
       Tcore = jets%Tcore_jets
       e6 = jets%e6_jets
    elseif (r >= jets%Rmin .and. r < jets%Rmax &
         .and. abs(cos_theta) > cos_theta_o_disk .and. &
         abs(cos_theta) < cos_theta_o_jets ) then
       Tcore = jets%Tcore_disk
       e6 = jets%e6_disk       
    else if (r >= jets%Rmin .and.  r < jets%Rdisk_min .and. &
         &  abs(cos_theta) < cos_theta_o_disk) then
       ! Spherical wind up to the disk radius
       Tcore = jets%Tcore_disk
       e6 = jets%e6_disk
    else
       Tcore = jets%Tcore_disk
       e6 = jets%e6_disk       
       r = jets%Rmax
    end if

    !
    Rp = jets%Rmin/r          
    out = Tcore*(Rp**e6) ! [10^4 K]                      
    

  END FUNCTION JetsTemperature

  
  
  SUBROUTINE calcJetsMassVelocity(thisOctal,subcell,grid)
    !  Assigings the density and velocity values to thisOctal at a given
    !  subcell using the density and the velocity functions as descrubed
    !  below.

    !    use inputs_mod
        
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    
    TYPE(vector) :: point
    
    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    TYPE(vector) :: vP

    REAL :: r    
    REAL :: Vr
        

    starPosn = grid%starPos1

    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn)
    r = modulus( pointVec ) 

                      
    ! test if the point lies within the star
    IF ( r >= jets%Rmin .AND. r < jets%Rmax) THEN
       thisOctal%inFlow(subcell) = .TRUE.       
    ELSE 
       thisOctal%inFlow(subcell) = .FALSE.
    END IF

    thisOctal%rho(subcell) = JetsDensity(pointVec, grid)/1.0d3
    ! [g/cm^3]

!!===========================================================
!! TEST FUNCTION!
!    thisOctal%rho(subcell) = density_inverse_sq(pointVec)
!    thisOctal%inFlow(subcell) = .TRUE.       
!!===========================================================   

    Vr = JetsVelocity(pointVec, grid)/ dble(cSpeed/1.0e5)
    ! in [c]
    
    vp = pointVec*real((Vr/r),kind=oct)  ! vector operation done here
    
    thisOctal%velocity(subcell) = vP
       
    
    IF (subcell == 8) CALL VelocityCorners(thisOctal,grid)
    
    

  END SUBROUTINE calcJetsMassVelocity
  

  SUBROUTINE calcJetsTemperature(thisOctal,subcell,grid)
  !===============================================================
  ! Assigns temperatures
  !
  ! For temperature structure of the jets is assumed as:
  !                Rmin
  !    T = Tcore*(------)^e6        [10^4 K]
  !                  R
  !
  !   Tcore  and e6 are to be assinged in the main parameter file.
  !
  ! NB: The temperature assigned to the cell is in [K].
  !===============================================================  
  
!  use inputs_mod    
  
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    
    TYPE(vector) :: point    
    TYPE(vector) :: starPosn
    TYPE(vector) :: pointVec
    
    starPosn = grid%starPos1
   
    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn)

    ! using the function in this module
    thisOctal%temperature(subcell) = JetsTemperature(pointVec, grid)*1.0e4  ! [K]
       
    ! we will initialise the bias distribution
    thisOctal%biasLine3D(subcell) = 1.0
  
  END SUBROUTINE calcJetsTemperature


  !
  !
  !
  SUBROUTINE initJetsAMR(grid)    
    use constants_mod
    use vector_mod
!    use inputs_mod

    implicit none
    
    type(GRIDTYPE), intent(inout)  :: grid
   
    real :: rStar
    
    !rStar  = jets%Rmin / 1.e10
    rStar  = jets%Rmin   ! this is already in [10^10cm]
    
    grid%rCore = rStar
    grid%rStar1 = rStar
    grid%starPos1 = vector(0.,0.,0.)


  end subroutine initJetsAMR


  ! Given a position in the model space and the direction specified
  ! by a directional vector (n), this function returns the directional
  ! derivative of the velocity field V(r) in this function. 
  function dV_dn_jets(r, n) RESULT(out)
    implicit none
    real :: out 
    !
    type(vector), intent(in) :: n  ! direction
    type(vector), intent(in) :: r  ! position

    real(double) :: rval ! size of r  [10^10cm]
    real(double) :: nval ! size of n
    real(double) :: Vinf, beta

    real(double) :: theta_o_jets, theta_o_disk
    real(double),save :: cos_theta_o_jets, cos_theta_o_disk
    logical, save :: first_time = .true.
    real(double) :: cos_theta


    IF (first_time) THEN
       ! Converting the half open angle to radians
       theta_o_jets = jets%theta_o_jets * Pi/180.0d0 
       cos_theta_o_jets = COS(theta_o_jets)

       theta_o_disk = jets%theta_o_disk * Pi/180.0d0 
       cos_theta_o_disk = COS(theta_o_disk)

       if (theta_o_jets > Pi/2.0d0 .or. theta_o_disk > Pi/2.0d0) then
          write(*,*) "Error:: theta_o_jets > Pi/2.0d0 .or. theta_o_disk > Pi/2.0d0."
          write(*,*) "        in [jets_mod::dv_dn_jets]."
          stop
       end if
       ! Quick check
       if (theta_o_disk < theta_o_jets) then
          print *, 'Error:: theta_o_disk must be greater than theta_o_jets.'
          print *, ' Exiting the program ...'
          stop
       end if
       first_time = .false.
    END IF



    ! Setting everything in consistent units.
    rval = dble(modulus(r))  ! [10^10cm]
    nval=  dble(modulus(n))



    if (r%z /= 0.0d0) then
       cos_theta = r%z/rval
    else
       cos_theta = 0.0d0
    end if
    !
    !
    if (rval >= jets%Rmin .and. rval < jets%Rmax &
       .and. abs(cos_theta) > cos_theta_o_jets) then

       Vinf = jets%Vinf_jets / 1.0d5  ! [10^10cm/s]
       beta = jets%beta_jets

    elseif (rval >= jets%Rmin .and. rval < jets%Rmax &
         .and. abs(cos_theta) > cos_theta_o_disk .and. &
         abs(cos_theta) < cos_theta_o_jets ) then

       Vinf = jets%Vinf_disk  / 1.0d5 ! [10^10cm/s]
       beta = jets%beta_disk

    else
       Vinf = jets%Vinf_disk       
       beta = jets%beta_disk

    end if

    

    if (rval <= jets%Rmin) then
       rval = jets%Rmin
       out =  ( jets%Rmin/(rval*rval) ) * Vinf * beta *(r.dot.n)/rval/nval
    else
       out =  ( jets%Rmin/(rval*rval) ) * Vinf * beta &
            * (1.0-jets%Rmin/rval)**(beta-1.0)  &
            * (r.dot.n)/rval/nval
    end if

    out = abs(out)      ! The units of the output should be sec^-1.

    if (out < 1.0e-10) out = 1.0e-10

    
  end function dV_dn_jets



  SUBROUTINE VelocityCorners(thisOctal,grid)
    ! store the velocity values at the subcell corners of an octal so
    !   that they can be used for interpolation.

    IMPLICIT NONE
  
    TYPE(octal), INTENT(INOUT) :: thisOctal
    TYPE(gridtype), INTENT(IN) :: grid

    real(oct)      :: x1, x2, x3
    real(oct)      :: y1, y2, y3
    real(oct)      :: z1, z2, z3
    

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
    
    thisOctal%cornerVelocity(1) = JetsVelocity_in_C(vector(x1,y1,z1),grid)
    thisOctal%cornerVelocity(2) = JetsVelocity_in_C(vector(x2,y1,z1),grid)
    thisOctal%cornerVelocity(3) = JetsVelocity_in_C(vector(x3,y1,z1),grid)
    thisOctal%cornerVelocity(4) = JetsVelocity_in_C(vector(x1,y2,z1),grid)
    thisOctal%cornerVelocity(5) = JetsVelocity_in_C(vector(x2,y2,z1),grid)
    thisOctal%cornerVelocity(6) = JetsVelocity_in_C(vector(x3,y2,z1),grid)
    thisOctal%cornerVelocity(7) = JetsVelocity_in_C(vector(x1,y3,z1),grid)
    thisOctal%cornerVelocity(8) = JetsVelocity_in_C(vector(x2,y3,z1),grid)
    thisOctal%cornerVelocity(9) = JetsVelocity_in_C(vector(x3,y3,z1),grid)

    ! middle level
  
    thisOctal%cornerVelocity(10) = JetsVelocity_in_C(vector(x1,y1,z2),grid)
    thisOctal%cornerVelocity(11) = JetsVelocity_in_C(vector(x2,y1,z2),grid)
    thisOctal%cornerVelocity(12) = JetsVelocity_in_C(vector(x3,y1,z2),grid)
    thisOctal%cornerVelocity(13) = JetsVelocity_in_C(vector(x1,y2,z2),grid)
    thisOctal%cornerVelocity(14) = JetsVelocity_in_C(vector(x2,y2,z2),grid)
    thisOctal%cornerVelocity(15) = JetsVelocity_in_C(vector(x3,y2,z2),grid)
    thisOctal%cornerVelocity(16) = JetsVelocity_in_C(vector(x1,y3,z2),grid)
    thisOctal%cornerVelocity(17) = JetsVelocity_in_C(vector(x2,y3,z2),grid)
    thisOctal%cornerVelocity(18) = JetsVelocity_in_C(vector(x3,y3,z2),grid)

    ! top level
    
    thisOctal%cornerVelocity(19) = JetsVelocity_in_C(vector(x1,y1,z3),grid)
    thisOctal%cornerVelocity(20) = JetsVelocity_in_C(vector(x2,y1,z3),grid)
    thisOctal%cornerVelocity(21) = JetsVelocity_in_C(vector(x3,y1,z3),grid)
    thisOctal%cornerVelocity(22) = JetsVelocity_in_C(vector(x1,y2,z3),grid)
    thisOctal%cornerVelocity(23) = JetsVelocity_in_C(vector(x2,y2,z3),grid)
    thisOctal%cornerVelocity(24) = JetsVelocity_in_C(vector(x3,y2,z3),grid)
    thisOctal%cornerVelocity(25) = JetsVelocity_in_C(vector(x1,y3,z3),grid)
    thisOctal%cornerVelocity(26) = JetsVelocity_in_C(vector(x2,y3,z3),grid)
    thisOctal%cornerVelocity(27) = JetsVelocity_in_C(vector(x3,y3,z3),grid)
    
  END SUBROUTINE VelocityCorners


  ! Use this for quick test function.
  function density_inverse_sq(point) RESULT(out)
    implicit none
    real                          :: out
    type(vector), intent(in) :: point
!    type(gridtype), intent(in)    :: grid

    real(oct)        :: rr                     !  [10^10 cm]
    real(oct)        :: r                      !  [-]
    
    real(oct), parameter  :: rmin = 1.0_oc     !  [-]
    real(oct), parameter  :: rmax = 101.0_oc   !  [-]
!    real(double), parameter  :: rmax = 100.0_oc      !  [-]
    real(oct), parameter  :: rho_0=1.0_oc      ! density at rmin

    rr  = modulus(point)  ! [10^10cm]
    r = rr/rmin          ! dimensionless radius
    if ( r >= rmin*0.5_oc  .and. r < rmax ) then
!    if ( r >= rmin  .and. r < rmax ) then
       out = rho_0/(r**2)
    else
       out = (rho_0/ (rmax**2))*1.0e-9
    end if

  end function density_inverse_sq



  
  
end module jets_mod

