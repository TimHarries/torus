program comparison

implicit none

real, parameter :: gridSize = 1.e19 
real, parameter :: r_c = 3.e8  !cylinder radius (10^10cm)
integer, parameter :: nLines = 201
integer,parameter :: tolerance = 5
real :: A_p(nLines)         !Torus pixel value
real :: E_p(nLines)         !Torus sampling Error
real :: A_a(nLines)         !Analytical pixel value
real :: analytical_error_plus(nLines)    !Upper limit on analytical pixel value
real :: analytical_error_minus(nLines)   !Lower limit on analytical pixel value
real :: analytical_error(nLines)         !Total pixel analytical error
real :: comparisonValue                  
real :: pixelValue, r, upperLimit, lowerLimit
real :: dx
integer :: i, nPhotons
integer :: ierr
integer :: nPassed


character(len=*), parameter :: pixelFile = "pixelFile.dat"
character(len=*), parameter :: analyticalFile = "analytical.dat"

dx = gridSize/real(nLines)

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
            A_p(i-1) = 1.0
            E_p(i-1) = 1.0
         else
            A_p(i-1) = pixelValue
            E_p(i-1) = 3./(sqrt(real(nPhotons)))
         end if
      else
         A_p(i-1) = pixelValue
         E_p(i-1) = 3./(sqrt(real(nPhotons)))
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
      A_a(i-1) = pixelValue
      if(upperLimit /= 0.0 .and. lowerLimit /= 0.0) then
         analytical_error_plus(i-1) = upperLimit
         analytical_error_minus(i-1) = lowerLimit
         analytical_error(i-1) = abs(analytical_error_plus(i-1) - &
              analytical_error_minus(i-1))/analytical_error_minus(i-1)
      end if
   else
      print *, "Error reading from :" //analyticalFile
      stop
   end if
end do

nPassed = 0

!Check pixels lie in acceptable range
do i = 1, nLines
   comparisonValue = abs(A_a(i-1) - A_p(i-1))/A_p(i-1)
   if(comparisonValue <= sqrt(E_p(i-1)**2 + analytical_error(i-1)**2)) then
      nPassed = nPassed + 1
   end if
end do

print *, "nPassed : ", nPassed
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
