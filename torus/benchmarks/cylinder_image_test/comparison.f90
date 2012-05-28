program comparison

implicit none

real, parameter :: gridSize = 1.e19 
real, parameter :: r_c = 3.e8  !cylinder radius (10^10cm)
integer, parameter :: nLines = 201
integer,parameter :: tolerance = 5
real, parameter :: n_sigma=3.0 
real :: A_p(nLines)         !Torus pixel value
real :: E_p(nLines)         !Torus sampling Error
real :: A_a(nLines)         !Analytical pixel value
real :: analytical_error_plus(nLines)    !Upper limit on analytical pixel value
real :: analytical_error_minus(nLines)   !Lower limit on analytical pixel value
real :: analytical_error(nLines)         !Total pixel analytical error
real :: comparisonValue                  
real :: pixelValue, r, upperLimit, lowerLimit
integer :: i, nPhotons
integer :: ierr
integer :: nPassed, nZero


character(len=*), parameter :: pixelFile = "pixelFile.dat"
character(len=*), parameter :: analyticalFile = "analytical.dat"


!open files
open(1, file=pixelFile, status="old", iostat = ierr)

if(ierr /= 0) then
   print *, "Error opening pixel file: " //pixelFile
   stop
end if

open(2, file=analyticalFile, status="old", iostat = ierr)
if(ierr /= 0) then
   print *, "Error opening analytical solution file: " //analyticalFile
   stop
end if

A_a = 0.0
A_p = 0.0
E_p = 0.0
analytical_error = 0.0
analytical_error_plus = 0.0
analytical_error_minus = 0.0


!read torus values
do i = 1, nLines
   read(1,*, iostat=ierr) r, pixelValue, nPhotons
   if(ierr == 0) then
      if(abs(r) > r_c) then
         if(pixelValue == 0.0 .and. nPhotons == 0) then
            A_p(i) = 1.0
            E_p(i) = 1.0
         else
            A_p(i) = pixelValue
            E_p(i) = n_sigma/(sqrt(real(nPhotons)))
         end if
      else
         A_p(i) = pixelValue
         E_p(i) = n_sigma/(sqrt(real(nPhotons)))
      end if
   else
      print *, "Error reading from :" //pixelFile   
      stop
   end if
end do 

!Read analytical values
do i = 1, nLines
   read(2,*, iostat=ierr) r, pixelValue, upperLimit, lowerLimit
   if(ierr == 0) then
      A_a(i) = pixelValue
      if(upperLimit /= 0.0 .and. lowerLimit /= 0.0) then
         analytical_error_plus(i) = upperLimit
         analytical_error_minus(i) = lowerLimit
         analytical_error(i) = abs(analytical_error_plus(i) - &
              analytical_error_minus(i))/analytical_error_minus(i)
      end if
   else
      print *, "Error reading from :" //analyticalFile
      stop
   end if
end do

nPassed = 0
nZero   = 0

!Check pixels lie in acceptable range
do i = 1, nLines
   if ( A_p(i) /= 0.0 ) then 
      comparisonValue = abs(A_a(i) - A_p(i))/A_p(i)
      if(comparisonValue <= sqrt(E_p(i)**2 + analytical_error(i)**2)) then
         nPassed = nPassed + 1
      end if
   else
      nZero = nZero + 1
   end if
end do

print *, "nPassed : ", nPassed
print *, "nZero : ", nZero
print *, "nFailed : ", nLines - nPassed
print *, "Tolerance: ", tolerance
if((nLines-nPassed) <= tolerance) then
   print *, "Test Successful"
else
   print *, "Test Failed"
end if

close(1)
close(2)

end program comparison
