real function getMeanMass(aMin, aMax, qDist)

  use constants_mod

  implicit none
  real :: aMin, aMax, qDist
  real :: tot, a1, a2, a, vol, fac, mass
  real :: ah, al, da, meanA
  integer :: i
  integer, parameter :: n = 1000
  real, parameter :: density = 2.   ! mean density of graphite 2 g /cm^3
     a1 = log10(aMin)
     a2 = log10(aMax)
     getMeanMass = 0.
     meanA = 0.
     tot = 0.
     do i = 1, n-1
        al = 10.**(a1 + (a2 - a1) * real(i-1)/real(n-1))
        ah = 10.**(a1 + (a2 - a1) * real(i)/real(n-1))
        a = (ah+al)/2.
        da = ah-al
        vol = (4./3.)* pi * a**3
        fac = da*0.5*(al**qDist + ah*qDist)
        mass = vol * density
        getMeanMass = getMeanMass + mass*fac
        meanA = meanA + a * fac
        tot = tot + fac
     enddo
     getMeanMass = getMeanMass / tot

end function getMeanMass






