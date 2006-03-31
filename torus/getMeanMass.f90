real function getMeanMass(aMin, aMax, a0, qDist, pDist)

  use constants_mod

  implicit none
  real :: aMin, aMax, a0, qDist, pDist
  real :: tot, a1, a2, a,vol, fac, mass
  real :: ah, al, da, meanA
  integer :: i
  integer, parameter :: n = 1000
  real :: density
  character(len=80) :: grainType

     select case(grainType)
       case("sil_dl")
          density = 3.6
       case("amc_hn","amc_zb")
          density = 2.
       case DEFAULT
          write(*,*) "!!! Unknown grain type in getMeanMass"
          density = 3.6
     end select
     a1 = log10(aMin)
     a2 = log10(aMax)
     getMeanMass = 0.
     meanA = 0.
     tot = 0.
     write(*,*) "qdist",qdist
     do i = 1, n-1
        al = 10.**(a1 + (a2 - a1) * real(i-1)/real(n-1))
        ah = 10.**(a1 + (a2 - a1) * real(i)/real(n-1))
        a = (ah+al)/2.
        da = ah-al
        vol = (4./3.)* pi * (a*microntocm)**3
        fac = da*0.5*(al**(-qDist)*exp(-(al/a0)**pDist) &
             + ah**(-qDist)* exp(-(ah/a0)**pDist))
        mass = vol * density
        getMeanMass = getMeanMass + mass*fac
        meanA = meanA + a * fac
        tot = tot + fac
     enddo
     getMeanMass = getMeanMass / tot
     meanA = meanA / tot
     write(*,*) "mean mass by number: ",getMeanMass
     write(*,*) "mean mass by radius",(4./3.)* pi * (meana*microntocm)**3 * density
end function getMeanMass

!
! This is similar to getMeanMass, but it uses more accurate integration method.
real function getMeanMass2(aMin, aMax, a0, qDist, pDist, graintype)  

  use constants_mod

  implicit none
  real, intent(in) :: aMin, aMax, a0, qDist, pDist
  real :: a1, a2, vol, fac
  integer :: i
  integer, parameter :: n = 1000
  real :: a(n)     ! grain sizes (log spaced)
  real :: f(n)     ! distribution function (normalized)
  real :: mass(n)  ! 
  real :: normFac
  character(len=*) :: grainType
  real :: density

  select case(grainType)
  case("sil_dl")
     density = 3.6
  case("amc_hn","amc_zb")
     density = 2.   ! mean density of graphite 2 g /cm^3
  case DEFAULT
     density = 3.6
     write(*,*) "==== WARNING ==== WARNING ==== WARNING ====="
     write(*,*) "Unknown grain type in getMeanMass2."
     write(*,*) "       grainType =", grainType
     write(*,*) "  Assuming the density of grain to be ", density, "[g/cm^3] "
     write(*,*) "       and continuing .... "
     write(*,*) "====================================== ====="
  end select


  a1 = log(aMin)
  a2 = log(aMax)
  !
  ! setting up the grain sizes
  do i = 1, n
     a(i) = (a1 + (a2 - a1) * real(i-1)/real(n-1))
     a(i) = exp(a(i))
     f(i) = a(i)**(-qDist) * exp(-(a(i)/a0)**pDist) 
  end do
  
  !
  ! normalize the dist function
  call PowerInt(n, 1, n, a, f, normFac)
  f(:) = f(:)/normFac

  !
  ! Finding the mean mass now.
  do i = 1, n
     vol = (4./3.)* pi * (a(i)*microntocm)**3
     mass(i) = vol * density * f(i)    ! weighted by dist function
  end do
  
  call PowerInt(n, 1, n, a, mass, fac)

  getMeanMass2 = fac

  !  !
  !  !  For debug
  !  !
  !  getMeanMass2 = getMeanMass2*1000.0
  
end function getMeanMass2



real function getMeanRadius(aMin, aMax, a0, qDist, pDist)

  use constants_mod

  implicit none
  real :: aMin, aMax, a0, qDist, pDist
  real :: tot, a1, a2, a, fac
  real :: ah, al, da, radius
  integer :: i
  integer, parameter :: n = 1000
     a1 = log10(aMin)
     a2 = log10(aMax)
     getMeanRadius = 0.
     tot = 0.
     do i = 1, n-1
        al = 10.**(a1 + (a2 - a1) * real(i-1)/real(n-1))
        ah = 10.**(a1 + (a2 - a1) * real(i)/real(n-1))
        a = (ah+al)/2.
        da = ah-al
        radius = a * microntocm
        fac = da*0.5*(al**(-qDist)*exp(-(al/a0)**pDist)  &
             + ah**(-qDist)*exp(-(ah/a0)**pDist))
        getMeanRadius = getMeanRadius + a*fac
        tot = tot + fac
     enddo
     getMeanRadius = getMeanRadius / tot

end function getMeanRadius




!
! This is similar to getMeanRaius, but it uses more accurate integration method.
real function getMeanRadius2(aMin, aMax, a0, qDist, pDist)  

  use constants_mod

  implicit none
  real, intent(in) :: aMin, aMax, a0, qDist, pDist
  real :: a1, a2,  vol, fac
  integer :: i
  integer, parameter :: n = 1000
  real :: a(n)     ! grain sizes (log spaced)
  real :: f(n)     ! distribution function (normalized)
  real :: radius(n)  ! 

  real :: normFac
  real, parameter :: density = 2.   ! mean density of graphite 2 g /cm^3
  a1 = log(aMin)
  a2 = log(aMax)
  !
  ! setting up the grain sizes
  do i = 1, n
     a(i) = (a1 + (a2 - a1) * real(i-1)/real(n-1))
     a(i) = exp(a(i))
     f(i) = a(i)**(-qDist) * exp(-(a(i)/a0)**pDist) 
  end do
  
  !
  ! normalize the dist function
  call PowerInt(n, 1, n, a, f, normFac)
  f(:) = f(:)/normFac

  !
  ! Finding the mean radius now.
  do i = 1, n
     radius(i) = a(i)*f(i)  ! weighted by dist function
  end do
  
  call PowerInt(n, 1, n, a, radius, fac)

  getMeanRadius2 = fac
  
end function getMeanRadius2




