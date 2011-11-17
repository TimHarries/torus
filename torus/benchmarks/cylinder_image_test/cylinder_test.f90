program cylinder_test

implicit none

double precision,parameter :: pi = 3.1415926535897932d0
integer, parameter :: nPixels = 201
double precision :: theta
double precision :: r_cylinder= 3.e18
double precision :: gridSize=1.e19
double precision :: r = 0.0
double precision :: dx
double precision :: I=0.0
double precision :: I_plus=0.0
double precision :: I_minus=0.0
double precision :: j 
double precision :: j_o = 5.830988926108036e32
double precision :: L = 3.e18
integer :: counter, ierr

dx = gridSize/201.

!Get average emission
j = j_o /(pi*L*(r_cylinder**2))

r = -(gridSize/2.0) + (dx/2.0)

open (1, file="analytical.dat", status="unknown", iostat=ierr)
if (ierr /= 0) then
   print *, "Trouble opening file: analytical.dat"
   stop
end if


do counter = 1, nPixels
!Pixel in cylinder
   if(abs(r) <= r_cylinder) then
     theta = asin(r/r_cylinder)
!Central 
     I = (j/(4.*pi))*(r_cylinder*cos(theta)*dx**2)
!Upper limit
     I_plus = (j/(4.*pi))*(dx*abs(tan(theta))+r_cylinder*cos(theta))*dx**2
!lower limit
     I_minus = (j/(4.*pi))*(r_cylinder*cos(theta)- dx*abs(tan(theta)))*dx**2      
   else
!Pixel outside of cylinder
      I = 0.0
      I_plus = 0.0
      I_minus = 0.0
   end if
! Output distance in Torus length units
   write(1, *) r/1.0e10, I, I_plus, I_minus
   r = r + dx
end do

close(1)

end program cylinder_test
