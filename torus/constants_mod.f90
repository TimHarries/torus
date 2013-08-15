!
! A module containing physical constants, mostly c.g.s.
!

module constants_mod



  
  use kind_mod

  public
  character(len=20) :: torusVersion
  
  real, parameter :: reallySmall = 1.e-30

  ! pi stuff
  
  real(double), parameter :: pi = 3.1415926535897932_db
  real(double), parameter :: piDouble = 3.1415926535897932_db
  real(double), parameter :: twoPi = 2.d0*pi
  real(double), parameter :: fourPi = 4.d0*pi
  real(double), parameter :: piBy4 = pi/4.d0
  real(double), parameter :: piBy2 = pi/2.d0
  real(double), parameter :: oneOnFourPi = 1.d0/(4.d0*pi)
  real(double), parameter :: oneOnTwoPi = 1.d0/(2.d0*pi)
  real(double), parameter :: degToRad = pi/180.d0
  real(double), parameter :: radToDeg = 180.d0/pi
  real(double), parameter :: oneByRootTwo = 1.d0/1.4142136d0
  real(double), parameter :: radiansToArcSec = radToDeg * 60.d0 * 60.d0
  real(double), parameter :: ArcSecsToRadians = pi/6.48d5 ! 6.48d5 = 180*60*60
  real(double), parameter :: arcsec = 1.d0/3600.d0 * degtorad
  real(double), parameter :: rootTwo = 1.4142136d0
  real(double), parameter :: oneOverrootTwo = 1.d0/1.4142136d0
  real(double), parameter :: sqrtPi = 1.77245385090551602730
  real(double), parameter :: OneOversqrtPi = 1.d0/1.77245385090551602730
  real(double), parameter :: OneOverPi = 1.d0 / Pi
  ! times

  real(double), parameter :: secsToYears = 1.d0/(365.25d0*24.d0*3600.d0)

  real(double), parameter :: barn = 1.d-24
  real(double), parameter :: megabarn = 1.d6 * barn

  ! lengths

  real(double), parameter :: angstrom = 1.d-10
  real(double), parameter :: micron = 1.d-6
  real(double), parameter :: angsToMicrons = angstrom/micron
  real(double), parameter :: micronsToAngs = micron/angstrom
  real(double), parameter :: centimeter = 1.d-2
  real(double), parameter :: micronTocm = micron/centimeter
  real(double), parameter :: angstromToCm = angstrom/centimeter
  real(double), parameter :: rSol = 6.96d10              ! cm
  real(double), parameter :: auToCm = 1.495979d13
  real(double), parameter :: pcToCm = 3.0856776d18
  real(double), parameter :: kpcToCm = pcToCm * 1000.0_db

  ! atomic

  real(double), parameter :: hConst = 6.626176d-34   ! Js
  real(double), parameter :: hCgs = 6.626205d-27
  real(double), parameter :: twoTimeshCgs = 2.d0 * hCgs  
  real(double), parameter :: hCgsOverfourPi = hCgs / fourPi
  real(double), parameter :: kConst = 1.380662d-23   ! J/k
  real(double), parameter :: kErg = 1.380626d-16     ! erg/k
  real(double), parameter :: hCgsOverKerg = hCgs/kErg
  real(double), parameter :: kev = 8.6171d-5         ! erg/k
  real(double), parameter :: sigmaE = 6.6525d-25       ! cm^2
  real(double), parameter :: eCharge = 4.803242384d-10 ! 
  real(double), parameter :: hydE0eV = 13.598433       ! eV
  real(double), parameter :: rydbergtoEv = 13.605692312d0 ! ev/ry
  real(double), parameter :: rydbergtoCm = 1.0973731568539d5 !cm^-1
  real(double), parameter :: hydE0eVdb = 13.598433_db! eV

  real(double), parameter :: ergToEv = 6.24145d11   
  real(double), parameter :: evtoerg = 1.d0/6.24145d11   

  real(double), parameter :: nuHydrogen = (hydE0eVdb/(ergToEv)) / (hCgs)

  real(double), parameter :: stefanBoltz = 5.6705d-5 ! erg/cm^2 s
  real(double), parameter :: lSol = 3.85d33          !erg/s
  ! speeds

  real(double), parameter :: cSpeed = 2.99792458d10  ! cm/s
  real(double), parameter :: cSpeedSquared = cSpeed**2
  real(double), parameter :: OneOvercSpeed = 3.33564095d-11  ! s/cm
  real, parameter :: cSpeed_sgl = 2.99792458e10  ! cm/s
  real(double), parameter :: cSpeed_dbl = 2.99792458d10  ! cm/s

  ! masses

  real(double), parameter :: mHydrogen = 1.6733d-24      ! g
  real(double), parameter :: amu = 1.6605402d-24         ! g
  real(double), parameter :: mSol = 1.9891d33            ! g
  real(double), parameter :: mEarth = 5.98d27            ! g
  real(double), parameter :: mMoon = 7.349d25            ! g
  real(double), parameter :: bigG = 6.67259d-8           ! cgs
  real(double), parameter :: mElectron = 9.109565D-28    ! g

  real(double), parameter :: aRad = 4.d0 * stefanBoltz / cSpeed
  real(double), parameter :: OneOveraRad = cSpeed / (4.d0 * stefanBoltz)

  real(double), parameter :: Tcbr = 2.728d0 !K
  real(double), parameter :: Navogadro = 6.0221415d23   ! per mol
  real(double), parameter :: Rgas = Navogadro * kerg    

  ! densities

  real(double), parameter :: msolpercAUtogpercc = msol / (autocm**3)    

  ! unit conversions
  real(double), parameter :: kmsToC = 1.d5/cSpeed

  ! refractive indices

  real(double) :: nAir = 1.000277

  !UV stuff
  real(double), parameter :: Habing = 1.6d-3 !erg cm-2 s-1
  real(double), parameter :: draine = 1.71*Habing
  real(double), parameter :: AV_fac = 6.289d-22
  real(double), parameter :: UV_fac = 3.02
  real(double), parameter :: zeta = 5.d-17/1.3d-17 !numerator is hardwired (would be user defined)

end module constants_mod


