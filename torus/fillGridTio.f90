subroutine fillGridTio(grid,  scale)

use grid_mod
use constants_mod
implicit none
real :: scale
type(GRIDTYPE) :: grid
real :: abundance
real :: massTio
real, allocatable :: xArray(:), sigmaArray(:)
real :: r, mu
integer :: i, j, k, m

massTio = 15.9994*mHydrogen + 47.867*mHydrogen
abundance = 1.e-7


allocate(xArray(1:grid%nlambda))
allocate(sigmaArray(1:grid%nlambda))

xArray = grid%lamarray


call readTioCrossSection(xArray, grid%nLambda, sigmaArray)

write(*,*) "Filling grid with TiO opacity"

if (.not.grid%cartesian) then
do i = 1, grid%nr
  do j = 1, grid%nmu
    do k = 1, grid%nPhi
     do m = 1, grid%nLambda
       grid%kappaAbs(i,j,k,m) = grid%kappaAbs(i,j,k,m) + &
            sigmaArray(m)*abundance*grid%rho(i,j,k)/massTio * scale
     enddo
    enddo
  enddo
enddo
else
do i = 1, grid%nx
  do j = 1, grid%ny
    do k = 1, grid%nz
     do m = 1, grid%nLambda
       r = grid%xAxis(i)**2 + grid%yAxis(j)**2 + grid%zAxis(k)**2
       mu = grid%zAxis(k)/r
       grid%kappaAbs(i,j,k,m) = grid%kappaAbs(i,j,k,m) + &
            sigmaArray(m)*abundance*grid%rho(i,j,k)/massTio * scale
     enddo
    enddo
  enddo
enddo
endif
end subroutine fillGridTio


