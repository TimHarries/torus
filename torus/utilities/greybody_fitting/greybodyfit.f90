!greybody fits to SEDs - SI version
!Author: Thomas Haworth (haworth@astro.ex.ac.uk)
!Created: 08/02/2012
!Last updated: 14/02/2012

!- based on what I read in Thompson et al. 2004 A&A 414, 1017
!- also uses Hildebrand, 1983, QJRAS, 24, 267 - the much better paper

!- not very pretty at present but I think it works

program greybodyfit
implicit none

double precision, parameter :: pi = 3.1415926535897932
double precision, parameter :: h = 6.626176d-34   ! Js 
double precision, parameter :: c = 2.99792458d8   ! ms^-1
double precision, parameter :: kB = 1.380662d-23
double precision, parameter :: um = 1.d-6
double precision, parameter :: degToRad = pi/180.d0
double precision, parameter :: radToDeg = 180.d0/pi 
double precision, parameter :: radiansToArcSec = radToDeg * 60.d0 * 60.d0
double precision, parameter :: ArcSecsToRadians = pi/6.48d5 
double precision, parameter :: arcsec = 1.d0/3600.d0 * degtorad


character(len=80) :: sedFilename
integer :: ier
integer :: numLines
integer :: i, j, k
integer :: numJobs
double precision :: lambda, flux, x1, x2, x3, x4
double precision, allocatable :: lambdaArray(:), fluxArray(:)
double precision :: Omega, T_dust
double precision :: Tmin, Tmax
double precision :: dT
double precision :: tau
double precision :: B
double precision :: C_nu
double precision :: F_nu
double precision :: minResidual, zerothResidual, minbeta
double precision :: distance, size, theta
double precision :: mass, totalMass, flux_ref, cutoff
integer :: minResidualID
double precision :: nu_ref, beamSize, lam_ref
double precision :: tau_ref, grad
double precision :: beta, percentDiff
double precision :: upperLimit
double precision, allocatable :: totalResidual(:), residualB(:)
logical :: ok, debug
integer :: choice


debug = .false.

!!write(*,*) "!- Enter thermal emission sed filename: "
!read(*,*) sedFilename
sedFilename="test_inc013_thermal_direct.dat"
!sedFilename="2Med_SED_conv_inc000_thermal_direct.dat"
!sedFilename="2Low_SED_conv_inc000.dat"


!write(*,*) "!- Enter number of lines in file: "
!read(*,*) numLines
numLines= 200

!write(*,*) "!- Enter starting search temperature (recommend 1 K): "
!read(*,*) Tmin
Tmin = 0.1d0

!write(*,*) "!- Enter final search temperature (recommend 1500K): "
!read(*,*) Tmax

Tmax = 1500.d0

!write(*,*) "!- Enter temperature search interval (recommend 0.1 K): "
!read(*,*) dT

dT = 0.1d0

!write(*,*) "!- Enter distance to object (m): "
!read(*,*) distance


distance = 3.0857d18
!distance = 3.08568025d18

!write(*,*) "!- Enter image size (m): "
!read(*,*) size

size = 4.0d14
!size = 1.5d17

beta = 2.d0

!calculate the solid angle subtended by the image

!theta = asin(size/(2.d0*distance))

theta = atan(size/distance)

Omega = 2.d0*pi*(1.d0 - cos(theta))

print *, "Solid angle subtended is: ", Omega

!lam_ref = 450.*um
!nu_ref = c/450.*um

lam_ref = 850.*Um
Nu_ref = c/(850.*um)
!cutoff = 450.*um
!cutoff = 8.*um
!cutoff = 10.d0*um
cutoff = 0.25*um
!beamSize = 17.3!*ArcSecsToRadians
beamSize = 22.9*ArcSecsToRadians
!beamSize = 25.*ArcSecsToRadians

print *, "nu_ref ", nu_ref

sedFilename = trim(sedFilename)


open(1, file=sedFilename, status="old", iostat=ier)
if(ier /= 0) then
   write(*,*) "!- Failed to open "//sedFilename
   stop
end if


allocate(fluxArray(numLines))
allocate(lambdaArray(numLines))

!read in the data from sed file
ok = .false.
do i = 1, numLines
   read(1,*) lambda, flux, x1, x2, x3, x4
   lambdaArray(i) = lambda*um
   fluxArray(i) = flux
   if(i /= 1) then
      if(lambdaArray(i) > (lam_ref) .and. .not. ok) then
         ok = .true.
         grad = (fluxArray(i) - fluxArray(i-1))/(lambdaArray(i) - lambdaArray(i-1))
         flux_ref = fluxArray(i-1) + (grad*((lam_ref) - lambdaArray(i-1)))
!         upperLimit = lambdaArray(i)
         upperLimit = 1.d10
         print *, "flux_ref", flux_ref
      end if
   end if
end do

close(1)

numJobs = int((Tmax-Tmin)/dT)

allocate(totalResidual(numJobs))
allocate(residualB(numJobs))


T_dust = Tmin
minResidual = 1.d20
totalResidual = 0.d0

!run check
call run_tests(Omega, lam_ref, flux_ref)

do i = 1, numJobs  
   totalResidual(i) = 0.d0
   zerothResidual = 0.d0
   do j = 1, numLines
      if(lambdaArray(j) > cutoff .and. lambdaArray(j) <= upperLimit) then
         if(fluxArray(j) /= 0.d0) then
            
            B = calcPlanck(T_dust, lam_ref, debug)
!            tau_ref = flux_ref/(pi*(beamSize**2)*B)
            tau_ref = flux_ref/(Omega*B)
            tau = tau_ref*((c/lambdaArray(j))/(nu_ref))**beta
            B = calcPlanck(T_dust, lambdaArray(j),debug)  
            F_nu = Omega*B*(1.d0 - exp(-tau))    
!            F_nu = Omega*B*tau    
            totalResidual(i) = totalResidual(i) + (abs(F_nu-fluxArray(j))/fluxArray(j))
            
            residualB(i) = residualB(i) +  (abs(fluxArray(j)-F_nu)/F_nu)
            zerothResidual = zerothResidual + fluxArray(j)
         end if
      end if
   end do
   t_dust = t_dust + dt
   if(totalResidual(i) < minResidual .and. totalResidual(i) /= zerothResidual) then
      percentDiff = (1.d0-(totalResidual(i)/minresidual))
!      print *, "temperature ", t_dust
!      print *, "improved by ", (1.d0-(totalResidual(i)/minresidual))*100.d0, "%"
      minResidual = totalResidual(i)
      minResidualID = i
!      minbeta = beta
!      print *, "minResidual ", minResidual
!      if(percentDiff < 1.d-6) then
!      if(percentDiff < 1.d-5) then
!         print *, "converged "
!         exit
!      end if
   end if
end do

!beta = minbeta

T_dust = Tmin + dble(minResidualID)*dT

write(*,*) "!- Derived dust temperature is: ", T_dust

!write(*,*) "!- Enter mass conversion factor (C_nu)"
!read(*,*) C_nu

!C_nu = 50.d0
C_nu = (50.d0)*1.d-7

!now to calculate the dust mass

B = calcPlanck(T_dust, lam_ref, debug)
if(B == 0.d0) then
   B = calcPlanck(T_dust, lam_ref, debug)
end if
Mass = (distance**2. * flux_ref*C_nu)/B

write(*,*) "Dust mass:"
write(*,'(5e14.5)') mass

B = calcPlanck(10.d0, lam_ref, debug)
if(B == 0.d0) then
   B = calcPlanck(10.d0, lam_ref, debug)
end if
write(*,*) "Mass at 10 K:"
write(*,'(5e14.5)') mass
Mass = (distance**2. * flux_ref*C_nu)/B

B = calcPlanck(50.d0, lam_ref, debug)
if(B == 0.d0) then
   B = calcPlanck(50.d0, lam_ref, debug)
end if
write(*,*) "Mass at 50 K:"
write(*,'(5e14.5)') mass
Mass = (distance**2. * flux_ref*C_nu)/B

!print *, "flux_ref ", flux_ref
!print *, "lam_ref", lam_ref
!print *, "B ", B
!print *, "distance ", distance



open(2, file="matched_spectrum.dat", status="unknown", iostat=ier)
open(3, file="residuals.dat", status="unknown", iostat=ier)
if(ier /= 0) then
   print *, "! ! ! Problem opening matched_spectrum.dat"
   stop
end if

print *, "Proceeding with temperature ", T_dust
print *, "Proceeding with wavelength ", lam_ref
print *, "Proceeding with tau_ref ", tau_ref

do i = 1, numLines
   B = calcPlanck(T_dust, lam_ref, debug)
!   tau_ref = flux_ref/(pi*(beamSize**2)*B)
   tau_ref = flux_ref/(Omega*B)
   if(lambdaArray(i) > cutoff .and. lambdaArray(i) <= upperLimit) then
      tau = tau_ref*((c/lambdaArray(i))/(nu_ref))**beta
      B = calcPlanck(T_dust, lambdaArray(i), debug)
      F_nu = Omega*B*(1.d0 - exp(-tau))
!      F_nu = Omega*B*tau       
      write(2,'(5e14.5)') lambdaArray(i)/um, F_nu
   end if
end do

do i = 1, numjobs
   write(3, *) i, totalResidual(i)
enddo

close(2)
close(3)

print *, "All tasks completed."

deallocate(lambdaArray)
deallocate(fluxArray)
deallocate(totalResidual)

contains


subroutine run_tests(Omega, lam_ref, flux_ref)
  implicit none
  double precision, parameter :: b = 2.8977685d-3
  double precision, parameter :: um = 1.d-6
  double precision, parameter :: pi = 3.1415926535897932
  double precision :: Omega
  double precision :: lam_ref, peak_lam, lam_max
  double precision :: flux_ref
  double precision :: B1, B2
  double precision :: beamSize
  double precision :: tau, lambda
  double precision :: tau_ref
  double precision :: F_nu, max_I, this_I, upperbound, lowerbound
  double precision :: lam_start, lam_end, dlam
  integer :: nJobs, i, ier
  logical :: debug
  debug = .false.

  B1 = calcPlanck(10.d0, lam_ref, debug)
!  tau_ref = flux_ref/(pi*(beamSize**2)*B1)
  tau_ref = (flux_ref/(Omega*B1))
  tau = tau_ref
  B2 = calcPlanck(10.d0, lam_ref,debug)
  F_nu = Omega*B2*(1.d0 - exp(-tau))    
!  F_nu = Omega*B2*tau    

  print *, " "
  print *, " "
  print *, " "
  print *, "*****************************"
  print *, "*     TEST RESULTS          *"
  print *, "*****************************"
  print *, " "
  print *, " "
  print *, "*****************************"
  print *, "*    FLUX INTERSECT TEST    *"
  print *, "*****************************"
  print *, "F_nu", F_nu
  print *, "flux_ref", flux_ref
  print *, "fractional difference: ", (F_nu-flux_ref)/F_nu
  print *, "tau_ref ", tau_ref
  print *, "tau ", tau
  print *, "Omega ", omega
  print *, "lam_ref ", lam_ref
  print *, "B1 ", B1
  print *, "B2 ", B2
  print *, "(1.d0 - exp(-tau))", (1.d0 - exp(-tau))
  print *, " "
  print *, " "


  if(abs((F_nu-flux_ref)/F_nu) > 1.d-2) then
     print *, "FAILED TEST"
     print *, "*****************************"
!     stop
  else
     print *, "Test passed."
  end if
  print *, "*****************************"
  print *, " "
  print *, " "

  print *, "*****************************"
  print *, "*     WIENS LAW TEST        *"
  print *, "*****************************"
  
  lam_start = 0.01*um
  lam_end = 500*um
  dlam = 0.01*um
  nJobs = (lam_end-lam_start)/dlam
  max_I = 0.d0
  lambda = lam_start

  open(9, file="spectrum.dat", status="unknown", iostat=ier)
  if(ier/=0) then
     print *, "failed to open spectrum.dat in tests"
     stop
  end if

  do i = 1, nJobs
     this_I =  calcPlanck(5778.d0, lambda, debug)
     if(this_I > max_I) then
        max_I = this_I
        peak_lam = lambda
     end if
     write(9,*) lambda, this_I
     lambda = lambda + dlam     
  end do !
  close(9)
  lam_max = b/5778.d0
  
  print *, "Wiens law predicts: ", lam_max
  print *, "my model gets: ", peak_lam
  print *, "fractional difference is: ", (peak_lam - lam_max)/lam_max
  upperbound = peak_lam + (dlam/2.d0)
  print *, "upper bound frac diff is: ", (upperbound - lam_max)/lam_max
  lowerbound = peak_lam - (dlam/2.d0)
  print *, "lower bound frac diff is: ", (lowerbound - lam_max)/lam_max

  if(abs((peak_lam - lam_max)/lam_max) > 1.d-2 .and. abs((peak_lam - lam_max)/lam_max) > 1.d-2 &
       .and. abs((peak_lam - lam_max)/lam_max) > 1.d-2) then
     print *, " "
     print *, " "    
     print *, "FAILED TEST"
     print *, " "
     print *, "lam_start ", lam_start
     print *, "lam_end ", lam_end
     print *, "nJobs ", nJobs
     print *, "max_I ", max_I    
     print *, "*****************************"
!     stop
  else
     print *, "Test passed."
  end if
  print *, "*****************************"
  print *, " "
  print *, " "

end subroutine


double precision function calcPlanck(temperature, lambda, debug)

  double precision :: temperature
  double precision, intent(in) :: lambda
  double precision, parameter :: h = 6.626176d-34   ! Js 
  double precision, parameter :: c = 2.99792458d8   ! ms^-1
  double precision, parameter :: kB = 1.380662d-23
  double precision :: fac1, fac2, nu
  logical :: debug, rjean

  nu = (c/lambda)


!nu
!  fac1 = ((2.d0*h*nu**3.d0)/(c**2.d0))*(1.d0/((exp(h*nu/(kB*temperature))-1.d0)))

!rayleigh-jeans
!     fac1 = 2.d0*(nu**2)*kB*Temperature/(c**2)

!lam
  fac1 = ((2.d0*h*c**2.d0)/(lambda**5.d0))*(1.d0/((exp((h*c)/(lambda*kB*temperature))-1.d0)))

  calcPlanck = fac1

  if(debug) then
     print *, "fac1 ", fac1
  end if

  
end function calcPlanck


end program greybodyfit
