program check

  implicit none
  logical :: passed
  real(kind=8) :: r(1000), phi_stars(1000),phi_gas(1000)
  real(kind=8) :: phiAnalytical
  real(kind=8), parameter :: bigG = 6.67259d-8
  real(kind=8), parameter :: msol = 1.9891d33
  real(kind=8), parameter :: pc = 3.0856776d18
  real(kind=8), parameter :: pi = 3.141592654d0
  real(kind=8) :: rhoMean, rToStar
  real(kind=8) :: junk1, junk2, junk3, junk4
  real(kind=8) :: maxfrac_gas,frac_gas
  real(kind=8) :: maxfrac_tot,frac_tot
  integer :: n, i

  rhoMean = mSol / ((4.d0/3.d0)*pi*(0.5d0*pc)**3)
  open(20,file="radial0001.dat", status="old",form="formatted")
  
  n = 1
  do 

     read(20,*,end=30) r(n), junk1, junk2, junk3, junk4,  phi_stars(n), phi_gas(n)
     n = n + 1
  enddo
  30 continue
  n = n - 1
  close(20)

  if ( n == 0 ) then
     write(*,*) "Torus test failed: file radial0001.dat is empty"
     STOP
  endif

  r(1:n) = r(1:n)*1.d10
  maxFrac_Gas = 0.d0
  maxFrac_tot = 0.d0
  do i = 1, n

     if (r(i) > 0.5d0*pc) then
        phiAnalytical = -bigG*mSol/(r(i))
    else
        phiAnalytical = (2.d0/3.d0)*pi*bigG*rhoMean*(r(i)**2 - 3.d0*(0.5d0*pc)**2)
     endif
     phiAnalytical = phiAnalytical-bigG*mSol/(sqrt(r(i)**2+5.d18**2))

     frac_gas = abs((phi_gas(i) - phiAnalytical)/phiAnalytical)
     maxFrac_Gas = max(maxFrac_gas, frac_gas)

     rToStar = sqrt(r(i)**2 + 3e18**2)
     phiAnalytical = phiAnalytical - bigG*0.5d0*msol/rToStar

     frac_tot = abs(((phi_gas(i)+phi_stars(i)) - phiAnalytical)/phiAnalytical)
     maxFrac_tot = max(maxFrac_tot, frac_tot)

     write(*,'(1p,5e12.4)') r(i), frac_gas ,frac_tot,phi_gas(i),phiAnalytical
  enddo
  passed = .true.
  if (maxFrac_gas > 0.01d0) then
     write(*,*) "Torus test failed on phi_gas"
     passed = .false.
  else
     write(*,*) "Torus test passed on phi_gas"
  endif

  if (maxFrac_tot > 0.01d0) then
     write(*,*) "Torus test failed on phi_i"
     passed = .false.
  else
     write(*,*) "Torus test passed on phi_i"
  endif

  if (passed) then
     write(*,*) "Torus gravity solver test successful."
  else
     write(*,*) "Torus gravity test failed"
  endif
end program check
