subroutine fillGridMie(grid, scale, aMin, aMax, qDist, grainType)

  use grid_mod
  use constants_mod
  implicit none
  type(GRIDTYPE) :: grid
  real :: aMin, aMax, qDist
  real, allocatable :: sigmaAbs(:), sigmaSca(:)
  real :: abundance
  real :: scale
  character(len=*) :: grainType
  integer :: i, j, k

  scale = 1.

  allocate(sigmaAbs(1:grid%nLambda))
  allocate(sigmaSca(1:grid%nLambda))

  abundance = 1.

  write(*,'(a)') "Filling grid with mie cross-sections..."
  call draineCrossSection(grid%lamArray, grid%nLambda, &
       aMin, aMax, qDist, &
       sigmaAbs, sigmaSca, grainType)

  !write(*,*) amin,amax,qdist
  !do i = 1, grid%nLambda
  !   write(*,*) grid%lamArray(i),sigmaAbs(i),sigmaSca(i)
  !enddo

  if (grid%cartesian) then

     do i = 1, grid%nx
        do j = 1, grid%ny
           do k = 1, grid%nz


              grid%kappaAbs(i,j,k,1:grid%nLambda) = sigmaAbs  * grid%rho(i,j,k)
              grid%kappaSca(i,j,k,1:grid%nLambda) = sigmaSca  * grid%rho(i,j,k)

              !      write(*,*) grid%kappaAbs(i,j,k,1:grid%nLambda),grid%kappaSca(i,j,k,1:grid%nLambda), grid%rho(i,j,k),scale
           enddo
        enddo
     enddo
  else

     do i = 1, grid%nr
        do j = 1, grid%nmu
           do k = 1, grid%nPhi

              grid%kappaAbs(i,j,k,1:grid%nLambda) = sigmaAbs * grid%rho(i,j,k)
              grid%kappaSca(i,j,k,1:grid%nLambda) = sigmaSca * grid%rho(i,j,k)

           enddo
        enddo
     enddo

  endif

  where(grid%kappaAbs < 1.e-25) grid%kappaAbs = 1.e-25
  where(grid%kappaSca < 1.e-25) grid%kappaSca = 1.e-25

  write(*,'(a)') "mie cross-sections done."
end subroutine fillGridMie

