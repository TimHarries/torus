program checkrainbow
implicit none
integer :: i, j, k(1)
character(len=80) :: junk
integer :: maxvalue, nx, ny
real :: angle1, angle2,s,r1
integer, allocatable :: r(:,:), g(:,:), b(:,:), slice(:)
open(20,file="detector.ppm",status="old",form="formatted")

read(20,'(a)') junk
read(20,'(a)') junk
read(20, '( i5, 1x, i5 )') nx,ny
read(20,'(i5)') maxvalue

allocate(r(nx,ny))
allocate(g(nx,ny))
allocate(b(nx,ny))
allocate(slice(0:(nx*2)))
do j = 1, ny
   do i = 1, nx
      read(20,'( 3(i5,1x) )') R(i,j),G(i,j),B(i,j)
   enddo
enddo
close(20)

slice = 0
do i = 1,nx
   do j = 1,ny
      k = nint(sqrt(real(i-nx/2)**2 + real(j-ny/2)**2))
      slice(k) = slice(k) + r(i,j)+g(i,j)+b(i,j)
   enddo
enddo

open(21,file="slice.dat",status="unknown",form="formatted")
do i = 1, nx*2
   write(21,*) i, slice(i)
enddo
close(21)
s = real(nx)/2.
k = maxloc(slice(1:nx)) 
r1 = real(k(1)) - 0.5
angle1 = atan2(r1/s*1.5,1.0)* 180./3.141592654

k = maxloc(slice(150:200))
k = k + 150
r1 = real(k(1)) - 0.5
angle2 = atan2(r1/s*1.5,1.0)* 180./3.141592654

   write(*,*) "primary rainbow at angles ",angle1
   write(*,*) "secondary rainbow at angles ",angle2


if ((angle1 > 40.89).and.(angle1 < 43.).and.(angle2 > 50.).and.(angle2 < 52.0)) then
   write(*,*) "TORUS: Test successful"
else
   write(*,*) "Rainbow benchmark failed"
endif
deallocate(r,g,b,slice)
end program checkrainbow

 
