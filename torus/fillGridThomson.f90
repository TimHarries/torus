subroutine fillGridThomson(grid)

use gridtype_mod
use grid_mod
use constants_mod
implicit none
type(GRIDTYPE) :: grid
real, parameter :: sigmaThomson = 6.65e-25
integer :: i, j, k, m


write(*,*) "Filling grid with electron scattering opacity"

if (grid%cartesian) then
do i = 1, grid%nx
  do j = 1, grid%ny
    do k = 1, grid%nz
     do m = 1, grid%nLambda
       grid%kappaSca(i,j,k,m) = grid%rho(i,j,k) * sigmaThomson 

     enddo
    enddo
  enddo
enddo

else

do i = 1, grid%nr
  do j = 1, grid%nmu
    do k = 1, grid%nPhi
     do m = 1, grid%nLambda
       grid%kappaSca(i,j,k,m) = grid%rho(i,j,k) * sigmaThomson
    enddo
  enddo
enddo
enddo
endif

end subroutine fillGridThomson


