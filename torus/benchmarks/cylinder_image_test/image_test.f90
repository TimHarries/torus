!Thaw - 16/07/2011

program compare_program

implicit none

real, parameter :: pi=3.14159
real, parameter :: angsToHz
real, parameter :: scaleFac=1.e20
real, parameter :: MJyConversionFactor=1.e17
real :: j_to           !Total emission across grid
real :: D              !Distance of observer from center of cylinder
real :: R              !Radius of cylinder
real :: F_tot_p        !Total predicting flux
real, allocatable :: F_s(:)            !Single cell flux
real, allocatable :: A(:)              !Pixel value
real :: steradPerPixel !Number of steradians in a pixel
real :: F_tot_sum
real :: comparison
real :: tolerance=1.0

integer :: i, j, line
integer :: iostat

character(len=*) :: pixelFile
character(len=*) :: valuesFile

pixelFile = "pixelData.dat"
valuesFile = "valuesData.dat"

!Open files
open(1, file=pixelFile, status="old", iostat = ierr)

if(ierr /= 0) then
   print *, "Error opening pixel file: " //pixelFile
   stop
end if

open(2, file=valuesFile, status="old", iostat = ierr)
if(ierr /= 0) then
   print *, "Error opening values file: " //valuesFile
   stop
end if

!Read in pixels
!Get nLines

line = 0
do
   read(1,*, iostat=ierr) pixelValue
   if(ierr < 0) then
      !eof
      exit
   else  if(ierr == 0) then
      line = line + 1
   else
      print *, "Error counting lines of : " //pixelFile
   end if
end do

allocate(A(Line))
allocate(F_s(Line))

A = 0.0
F_s = 0.0

do i = 1, line
   read(1,*, iostat=ierr) pixelValue
   if(ierr /= 0) then
      A(i-1) = pixelValue
   else
      print *, "Error reading from :" //pixelFile   
      stop
   end if
end do 

!Read in values
read(2, *, iostat=ierr) j_tot, D, R, steradPerPixel
if(ierr /= 0) then
   print *, "Error reading from :" //valuesFile
   stop
end if

!Expected result

F_tot_p = j_tot/(2.*(pi**2)*(D**0.5)*(R**(3./2.)))

!Convert pixel values 

do i = 1, nLines
   F_s(i-1) = (A(i-1) * steradPerPixel)/(MJyConvertsionFactor*angsToHz*scaleFac)
end do

!Comparisons
F_tot_sum = sum(F_s)

comparison = (1.0 - (F_tot_sum/F_tot_p))*100.0

print *, "Total image fluxes differ by: ", comparison,"%"
print *, "Tolerance is : ", tolerance, "%"
if(comparsion < tolerance) then
   print *, "Test Successful"
else
   print *, "Test Failed"
end if

close(1)
close(2)
deallocate(A)
deallocate(F_s)

end program compare_program
