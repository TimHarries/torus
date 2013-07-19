! This module relates to angularImage_mod and contains things which are 
! needed before amr_mod is built
!
! D. Acreman, July 2013

module angularImage_utils

  use kind_mod

  implicit none

  public :: calcAngImgTest

! Is the angular image module in use? 
  logical, public, save :: internalView=.false.

  logical, public, save :: obsVelFromGrid ! Is observer velocity taken from Torus grid?
  logical, public, save :: thermalLineWidth ! Use thermal line width? 
  logical, public, save :: SplitCubes ! Split cube into +ve and -ve contributions? 
  logical, public, save :: refineQ2Only ! Limit grid refinement to 2nd quadrant
  real(double), public, save :: intPosX, intPosY, intPosZ ! Position of internal observer
  real(double), public, save :: intDeltaVx, intDeltaVy, intDeltaVz ! Additional velocity applied to internal observer
  real(double), public, save :: galaxyInclination, galaxyPositionAngle
  integer, public, save :: dssMinSample ! Minimum number of samples per cell when density subsample is used

  contains

! Set up octals for the angular image test case
    subroutine calcAngImgTest(thisOctal, subcell)

      use octal_mod

! Arguments
      TYPE(OCTAL), intent(inout) :: thisOctal
      integer, intent(in) :: subcell

! Local variables
      TYPE(vector) :: thisCentre
      real(double) :: r, z, theta
      real(double) :: vx, vy
! Radius of model "galaxy" in kpc
      real(double), parameter :: r_gal=10.0
! Height of model "galaxy" in kpc
      real(double), parameter :: z_gal=1.0
! Density inside the galaxy
      real(double), parameter :: rho_in=1.0e-23
! Density outside the galaxy
      real(double), parameter :: rho_out=1.0e-33
! Circular velocity in km/s
      real(double), parameter :: v_circ = 220.0 
! Local velocity
      TYPE(VECTOR) :: thisVel

! Calculate cylindrical polar co-ordinates of this cell
      thisCentre = subcellcentre(thisOctal, subcell)
      r = sqrt( thisCentre%x**2 + thisCentre%y**2 )     ! distance from grid centre in torus units
      theta = asin(thisCentre%y/r)
      r = r * 1.0e10_db / kpcToCm ! convert to kpc
      z = abs( thisCentre%z  * 1.0e10_db / kpcToCm )
      

      if (.not. associated(thisoctal%cornervelocity) ) allocate(thisoctal%cornervelocity(27))
      if (.not. associated(thisoctal%cornerrho)      ) allocate(thisoctal%cornerrho(27))

      if ( r < r_gal .and. z < z_gal ) then 
         thisOctal%rho(subcell) = rho_in
         thisoctal%cornerrho(:) = rho_in
         thisOctal%temperature(subcell) = 30.0
      else
         thisOctal%rho(subcell) = rho_out
         thisoctal%cornerrho(:) = rho_out
         thisOctal%temperature(subcell) = 1.0e4_db
      end if

      if (thisCentre%x > 0.0 ) then
         vx =        (v_circ*1.0e5/cspeed) * sin(theta)
         vy = -1.0 * (v_circ*1.0e5/cspeed) * cos(theta)
      else
         vx =        (v_circ*1.0e5/cspeed) * sin(theta)
         vy =        (v_circ*1.0e5/cspeed) * cos(theta)
      endif

      thisVel = VECTOR(vx,vy,0.0)
      thisOctal%velocity= thisVel

! Set corner values from first subcell. This test case uses uniform values so they don't vary across an octal.
      if (subcell == 1) then 
         thisoctal%cornervelocity(:) = thisoctal%velocity(1)
         thisoctal%cornerrho(:)      = thisoctal%rho(1)
      endif

    end subroutine calcAngImgTest

  end module angularImage_utils

