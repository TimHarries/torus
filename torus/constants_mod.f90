!
! A module containing physical constants, mostly c.g.s.
!

module constants_mod

  public
  
  real, parameter :: reallySmall = 1.e-30

  ! pi stuff
  
  real, parameter :: pi = 3.141592654
  real, parameter :: twoPi = 2.*pi
  real, parameter :: fourPi = 4.*pi
  real, parameter :: piBy4 = pi/4.
  real, parameter :: oneOnFourPi = 1./(4.*pi)
  real, parameter :: oneOnTwoPi = 1./(2.*pi)
  real, parameter :: degToRad = pi/180.
  real, parameter :: radToDeg = 180./pi
  real, parameter :: oneByRootTwo = 1./1.4142136
  real, parameter :: radiansToArcSec = radToDeg * 60. * 60.
  ! times

  real, parameter :: secsToYears = 1./(365.25*24.*3600.)
  

  ! lengths

  real, parameter :: angstrom = 1.e-10
  real, parameter :: micron = 1.e-6
  real, parameter :: angsToMicrons = angstrom/micron
  real, parameter :: centimeter = 1.e-2
  real, parameter :: micronTocm = micron/centimeter
  real, parameter :: angstromToCm = angstrom/centimeter
  real, parameter :: rSol = 6.96e10              ! cm
  real, parameter :: auToCm = 1.495979e13
  real, parameter :: pcToCm = 3.0856776e18

  ! atomic

  real, parameter :: hConst = 6.626176e-34   ! Js
  real, parameter :: hCgs = 6.626205e-27
  real, parameter :: kConst = 1.380662e-23   ! J/k
  real, parameter :: kErg = 1.380626e-16     ! erg/k
  real, parameter :: kev = 8.6171e-5         ! erg/k
  real, parameter :: sigmaE = 6.65e-25       ! cm^2
  real, parameter :: eCharge = 4.803242384E-10 ! 

  real, parameter :: ergToEv = 6.24145e11   

  real, parameter :: stefanBoltz = 5.6705e-5 ! erg/cm^2 s
  real, parameter :: lSol = 3.85e33          !erg/s
  ! speeds

  real, parameter :: cSpeed = 2.99792458e10  ! cm/s

  ! masses

  real, parameter :: mHydrogen = 1.67e-24        ! g
  real, parameter :: mSol = 1.9891e33            ! g
  real, parameter :: mEarth = 5.98e27            ! g
  real, parameter :: bigG = 6.67259e-8           ! cgs
  real, parameter :: mElectron = 9.109565D-28    ! g

end module constants_mod


