module jets_mod
  ! Current version handle the following geometry case:
  !
  !  1. bipolar jets
  !  
  ! --------------------------------------------------------------------------------
  !  Created: OCT-03-2002  (R.Kurosawa)
  !
  !---------------------------------------------------------------------------------

  USE vector_mod          ! vector maths routines
  USE kind_mod            ! variable kind parameters
  USE octal_mod           ! type definition for the grid elements
  USE gridtype_mod        ! type definition for the 3-d grid
  USE constants_mod

    
  public :: set_jets_parameters, get_jets_parameter, &
       &    JetsDensity, JetsVelocity, &
       &    calcJetsMassVelocity, calcJetsTemperature, initJetsAMR, &
       &    dV_dn_jets

  !  privare :: 
  
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
     real :: Vinf    ! terminal velocity in  [kms]
     real :: beta    ! beta in the beta-belocity law [-]
     real :: Vo      ! a small offset in V in  [km/s]
     real :: Mdot    ! mass loss rate in jets [M_solar/yr]
     real :: theta_o ! half opening angle [degrees]
     
     real :: Tcore   ! temperature at the core at Rmin in [10^4 K]
     real :: e6      ! an exponet of in tempereture eq. [-]
  end type jets_parameters
   
  type(jets_parameters) :: jets  ! One incidence to be used ONLY IN THIS MODULE.

  !
  !
  !
  
contains

  ! set the jets_parameters.
  subroutine set_jets_parameters(Rmin, Rmax, Vinf, beta, &
       & Vo, Mdot, theta_o, Tcore, e6)
    
    implicit none
    real,intent(in) :: Rmin    ! radius of central star  [10^10cm]
    real,intent(in) :: Rmax    ! cutoff radius in [10^10 cm]
    real,intent(in) :: Vinf    ! terminal velocity in  [kms]
    real,intent(in) :: beta    ! beta in the beta-belocity law [-]
    real,intent(in) :: Vo      ! a small offset in V in  [km/s]
    real,intent(in) :: Mdot    ! mass loss rate in jets [M_solar/yr]
    real,intent(in) :: theta_o ! half opening angle [degrees]
    
    real,intent(in) :: Tcore   ! temperature at the core at Rmin in [10^4 K]
    real,intent(in) :: e6      ! an exponet of in tempereture eq. [-]

    
    jets%Rmin     =     Rmin   
    jets%Rmax     =     Rmax   
    jets%Vinf     =     Vinf   
    jets%beta     =     beta   
    jets%Vo       =     Vo     
    jets%Mdot     =     Mdot   
    jets%theta_o  =  theta_o
    jets%Tcore    =    Tcore  
    jets%e6       =    e6

  end subroutine set_jets_parameters

  ! access a parameter in jets_parameters
  function get_jets_parameter(name)  RESULT(out)
    implicit none
    real :: out
    character(LEN=*) :: name
    
    select case(name)
    case ("Rmin")
       out = jets%Rmin
    case ("Rmax")
       out = jets%Rmax
    case ("Vinf")
       out = jets%Vinf
    case ("beta")
       out = jets%beta
    case ("Vo")
       out = jets%Vo
    case ("Mdot")
       out = jets%Mdot
    case ("theta_o")
       out = jets%theta_o
    case ("Tcore")
       out = jets%Tcore
    case ("e6")
       out = jets%e6
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
    
    !    USE input_variables  ! inputs needed for this function are defined in here.
    USE constants_mod
    
    IMPLICIT NONE

    REAL                          :: out
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec

    !
    LOGICAL, SAVE  :: first_time = .true. 
    DOUBLE PRECISION ::  r  ! distance from the center of central star
    DOUBLE PRECISION :: rp
    DOUBLE PRECISION :: vr
    DOUBLE PRECISION  :: cos_theta, theta_o
    !
    DOUBLE PRECISION, SAVE :: cos_theta_o    
    DOUBLE PRECISION, SAVE :: Mdot_jet             ! [kg/sec]
    DOUBLE PRECISION, parameter :: M_sun=1.9891d30 ! in [kg]

    IF (first_time) THEN
       ! Converting the half open angle to radians
       theta_o = jets%theta_o * Pi/180.0d0 
       cos_theta_o = COS(theta_o)
       ! Converting Mdot to kg/s
       Mdot_jet = jets%Mdot*M_sun/(60.0d0*60.0d0*24.0d0*365.0d0) ! [kg/sec]
       first_time = .false.
    END IF
    !
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn)
    r = modulus( pointVec ) 
    
    
    IF (pointVec%z /= 0.0) THEN
       cos_theta = pointVec%z/r
    ELSE
       cos_theta = 0.0d0
    END IF
    IF (r >= jets%Rmin .AND. r < jets%Rmax .AND. ABS(cos_theta) > cos_theta_o) THEN
       Rp = jets%Rmin/r  
       Vr = jets%Vinf * (1.0d0 - Rp)**(jets%beta)  + jets%Vo     ! [km/s]
       ! density in [kg/m^3]
       !-- convert r into [m] ... r was in 10^10cm
       r = r * 1.0d8
       
       ! convert [km/s] to [m/s]
       Vr = Vr*10.0d3  ! m/s
       out  =                    Mdot_jet               &
	    &                     /                     &
	    &           ( 4.0d0*Pi*r*r*Vr*(1.0d0-cos_theta_o) )  ! [kg/m^3]      
    ELSE
       Rp = jets%Rmin/jets%Rmax
       Vr = jets%Vinf * (1.0d0 - Rp)**jets%beta +jets%Vo    ! [km/s]
       ! density in [kg/m^3]
       !-- convert r into [m] was in [10^10cm]
       r = jets%Rmax * 1.0d8
       
       ! convert [km/s] to [m/s]
       Vr = Vr*10.0d3  ! m/s
       
       out =                      Mdot_jet             &
	    &                       /                  &
	    &        ( 4.0d0*Pi*r*r*Vr*(1.0d0-cos_theta_o)  )
       
       out =  out*1.0d-3
       
    END IF

    
  END FUNCTION JetsDensity


  FUNCTION JetsVelocity(point,grid) RESULT(out)
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
    !
    
!    USE input_variables  ! inputs needed for this function are defined in here.
    USE constants_mod
    
    IMPLICIT NONE

    REAL                          :: out
    TYPE(octalVector), INTENT(IN) :: point
    TYPE(gridtype), INTENT(IN)    :: grid

    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec

    !
    LOGICAL, SAVE  :: first_time = .true. 
    DOUBLE PRECISION ::  r  ! distance from the center of central star
    DOUBLE PRECISION :: Vr
    DOUBLE PRECISION :: Rp
    DOUBLE PRECISION  :: cos_theta, theta_o
    !
    DOUBLE PRECISION, SAVE :: cos_theta_o    
    DOUBLE PRECISION, parameter :: M_sun=1.9891d30 ! in [kg]

    IF (first_time) THEN
       ! Converting the half open angle to radians
       theta_o = jets%theta_o * Pi/180.0d0 
       cos_theta_o = COS(theta_o)
       first_time = .false.
    END IF
    !
    
    starPosn = grid%starPos1
    pointVec = (point - starPosn)
    r = modulus( pointVec )   
    
    IF (pointVec%z /= 0.0) THEN
       cos_theta = pointVec%z/r
    ELSE
       cos_theta = 0.0d0
    END IF
    IF (r >= jets%Rmin .AND. r < jets%Rmax .AND. ABS(cos_theta) > cos_theta_o) THEN
       Rp = jets%Rmin/r  
    ELSE
       Rp = jets%Rmin/jets%Rmax
    END IF

    Vr = jets%Vinf * (1.0d0 - Rp)**(jets%beta)  + jets%Vo     ! [km/s]
    out = Vr ! [km/s]

    
  END FUNCTION JetsVelocity


  
  
  SUBROUTINE calcJetsMassVelocity(thisOctal,subcell,grid)
    !  Assigings the density and velocity values to thisOctal at a given
    !  subcell using the density and the velocity functions as descrubed
    !  below.

!    use input_variables
	
    IMPLICIT NONE
    
    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    
    TYPE(octalVector) :: point
    
    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    TYPE(vector) :: vP

    REAL :: r    
    REAL :: Vr
        
    starPosn = grid%starPos1

    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn)
    r = modulus( pointVec ) 

		      
    ! test if the point lies within the star
    IF ( r >= jets%Rmin .AND. r < jets%Rmax) THEN
       thisOctal%inFlow(subcell) = .FALSE.       
    ELSE 
       thisOctal%inFlow(subcell) = .TRUE.
    END IF

    thisOctal%rho(subcell) = JetsDensity(pointVec, grid)/1.0d3
    ! [g/cm^3]
    
    Vr = JetsVelocity(pointVec, grid)/ dble(cSpeed/1.0e5)
    ! in [c]
    
    vp = pointVec*(Vr/r)  ! vector operation done here
    
    thisOctal%velocity(subcell) = vP
       

    

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
  
!  use input_variables    
  
    IMPLICIT NONE

    TYPE(octal), INTENT(INOUT) :: thisOctal
    INTEGER, INTENT(IN) :: subcell
    TYPE(gridtype), INTENT(IN) :: grid
    
    TYPE(octalVector) :: point    
    TYPE(octalVector) :: starPosn
    TYPE(octalVector) :: pointVec
    
    REAL :: Rp, cos_theta
    real :: r
    
    starPosn = grid%starPos1
   
    point = subcellCentre(thisOctal,subcell)
    pointVec = (point - starPosn)
    r = modulus( pointVec )
    
    IF (pointVec%z /= 0.0) THEN
       cos_theta = pointVec%z/r
    ELSE
       cos_theta = 0.0d0
    END IF

    
    IF (r >= jets%Rmin .AND. r < jets%Rmax .and. &
	             &  ABS(cos_theta) > Cos(jets%theta_o)) THEN
       Rp = jets%Rmin/r
    ELSE
       Rp = jets%Rmin/jets%Rmax
    END IF

    thisOctal%temperature(subcell) = jets%Tcore*1.0e4 &
	 &  * Rp**(jets%e6)    ! [K]
       
    
    ! we will initialise the bias distribution
    thisOctal%biasLine3D(subcell) = 1.0
  
  END SUBROUTINE calcJetsTemperature

  
  SUBROUTINE initJetsAMR(grid)    
    use constants_mod
    use vector_mod
!    use input_variables
    use parameters_mod

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

    double precision :: rval ! size of r  [10^10cm]
    double precision :: nval ! size of n
    double precision :: Vinf
    double precision  :: Rmin

    ! Setting everything in consistent units.
    rval = dble(modulus(r))  ! [10^10cm]
    nval=  dble(modulus(n))
    Vinf = jets%Vinf / 1.0d5 ! [10^10cm/s]
    Rmin = jets%Rmin 

    if (rval <= Rmin) rval = Rmin
    
    out =  ( Rmin/(rval*rval) ) * Vinf * jets%beta * (1.0-Rmin/rval)**(jets%beta-1.0)  &
	 &  *abs(r.dot.n)/rval/nval

    ! The units of the output should be sec^-1.
    
  end function dV_dn_jets
  
  
end module jets_mod



