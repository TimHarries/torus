
subroutine fillGridTorus(Grid, xAxis, yAxis, zAxis, nx, ny, nz, rho, scale, &
                         rTorus, rOuter)

use vector_mod
use constants_mod

 implicit none
 integer :: nx, ny, nz
 integer, parameter :: nRad = 10
 real :: thisRad
 real :: Grid(nx, ny, nz)
 real :: xAxis(nx)
 real :: yAxis(nx)
 real :: zAxis(nx)
 real :: theta, phi
 real :: rTorus
 real :: rOuter
 type(VECTOR) :: r0Vec, rVec, aVec, bVec, cVec, rHat
 type(VECTOR) :: torusAxis, normVec
 integer, parameter :: nAng = 10
 integer :: i, j, k
 integer :: i1, i2, i3
 real :: rho, scale


 do i = 1, nx
    xAxis(i) = 2.*real(i-1)/real(nx-1) - 1.
 enddo


 do i = 1, ny
    yAxis(i) = 2.*real(i-1)/real(ny-1) - 1.
 enddo


 do i = 1, nz
    zAxis(i) = 2.*real(i-1)/real(nz-1) - 1.
 enddo

 xAxis = xAxis * scale
 yAxis = yAxis * scale
 zAxis = zAxis * scale

 rTorus = rTorus * scale
 rOuter = rOuter * scale

 torusAxis%x = 0.
 torusAxis%y = 0.
 torusAxis%z = 1.

 do i = 1, nAng

	theta  = real(i-1)/real(nAng-1)*twoPi
  	r0Vec%x = rTorus * cos(theta)
	r0Vec%y = rTorus * sin(theta)
	r0Vec%z = 0.
	normVec = crossProd(r0Vec, torusAxis)
	call normalize(normVec)

	rHat = r0Vec
	call normalize(rHat)

      do k = 1 , nRad
	thisRad = rOuter*real(k-1)/real(nRad-1)
        do j = 1, nAng
	  phi  = real(j-1)/real(nAng-1)*twoPi        
          aVec = mult(thisRad,mult(cos(phi),rHat))
	  bVec = mult(thisRad,mult(sin(phi),torusAxis))
	  cVec = add(aVec, bVec)
 	  rVec = add(r0Vec, cVec)
	  call locate(xAxis,nx,rVec%x, i1)
	  call locate(yAxis,ny,rVec%y, i2)
	  call locate(zAxis,nz,rVec%z, i3)
	  Grid(i1,i2,i3) = rho
        enddo
      enddo
 enddo

end subroutine fillGridTorus
		
