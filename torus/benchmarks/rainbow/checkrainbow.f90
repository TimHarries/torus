program checkrainbow
implicit none
integer :: i, j, k(1)
character(len=80) :: junk
integer :: maxvalue, nx, ny
real :: angle1, angle2,angle3,angle4,s,r1
integer, allocatable :: r(:,:), g(:,:), b(:,:), slice(:)
open(20,file="detector.ppm",status="old",form="formatted")

read(20,'(a)') junk
read(20,'(a)') junk
read(20, '( i5, 1x, i5 )') nx,ny
read(20,'(i5)') maxvalue

allocate(r(nx,ny))
allocate(g(nx,ny))
allocate(b(nx,ny))
allocate(slice(1:nx))
do j = 1, ny
   do i = 1, nx
      read(20,'( 3(i5,1x) )') R(i,j),G(i,j),B(i,j)
   enddo
enddo
close(20)

slice = 0
do i = 1,nx
   slice(i) = r(i,ny/2)+g(i,ny/2)+b(i,ny/2)
enddo
do i = 1,ny
   slice(i) = slice(i) + r(nx/2,i)+g(nx/2,i)+b(nx/2,i)
enddo


open(21,file="slice.dat",status="unknown",form="formatted")
do i = 1, nx
   write(21,*) i, slice(i)
enddo
close(21)
s = real(nx)/2.
k = maxloc(slice(nx/2:nx)) 
r1 = real(k(1)) - 0.5
angle1 = atan2(r1/s*1.5,1.0)* 180./3.141592654

k = maxloc(slice(1:nx/2))
r1 = real(k(1)) - 0.5
angle2 = atan2((s - r1)/s*1.5,1.0)* 180./3.141592654


k = 5*nx/6 + maxloc(slice((5*nx/6):nx)) 
r1 = real(k(1)) - 0.5
angle3 = atan2(real(r1-s)/s*1.5,1.0)* 180./3.141592654

k = maxloc(slice(1:(nx/6)))
r1 = real(k(1)) - 0.5
angle4 = atan2((s - r1)/s*1.5,1.0)* 180./3.141592654


   write(*,*) "primary rainbow at angles ",angle1,angle2
   write(*,*) "secondary rainbow at angles ",angle3,angle4


if ((angle1 > 40.89).and.(angle1 < 43.).and.(angle2 > 40.89).and.(angle2 < 43.).and. &
     (angle3 > 50.).and.(angle3 < 52.0).and.(angle4>50.).and.(angle4 < 52.)) then
   write(*,*) "TORUS: Test successful"
else
   write(*,*) "Rainbow benchmark failed"
endif


end program checkrainbow

 
