subroutine fillGridRayleigh(grid,  scale)

use gridtype_mod
use constants_mod
implicit none
type(GRIDTYPE) :: grid
real :: scale
real :: abundance
real :: sigmaRayleigh
integer :: i, j, k, m

abundance = 1.

if (grid%cartesian) then
write(*,*) "Filling grid with Rayleigh scattering opacity"
do i = 1, grid%nx
  do j = 1, grid%ny
    do k = 1, grid%nz
     do m = 1, grid%nLambda
       grid%kappaSca(i,j,k,m) = real(grid%kappaSca(i,j,k,m) + sigmaRayleigh(grid%lamArray(m))*&
            abundance*grid%rho(i,j,k)/mHydrogen*scale)
     enddo
    enddo
  enddo
enddo

else

do i = 1, grid%nr
  do j = 1, grid%nmu
    do k = 1, grid%nPhi
     do m = 1, grid%nLambda
       grid%kappaSca(i,j,k,m)=real(grid%kappaSca(i,j,k,m)+sigmaRayleigh(grid%lamArray(m))*abundance*grid%rho(i,j,k)/mHydrogen*&
           scale)
     enddo
    enddo
  enddo
enddo

endif

end subroutine fillGridRayleigh

real function sigmaRayleigh(lambda)
implicit none
real lambda

sigmaRayleigh = 0.66520e-24*(1026./lambda)**4

end function sigmaRayleigh


