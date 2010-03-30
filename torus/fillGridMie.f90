
subroutine fillGridMie2(grid, scale, aMin, aMax, a0, qDist, pDist, grainType)

  use gridtype_mod
  use constants_mod
  implicit none
  type(GRIDTYPE) :: grid
  real :: aMin, aMax, a0, qDist, pDist
  real, allocatable :: sigmaAbs(:), sigmaSca(:)
  real :: scale
  character(len=*) :: grainType
  integer :: i, j, k

  scale = 1.

  allocate(sigmaAbs(1:grid%nLambda))
  allocate(sigmaSca(1:grid%nLambda))

  write(*,'(a)') "Filling grid with mie cross-sections..."
  call draineCrossSection(grid%lamArray, grid%nLambda, &
       aMin, aMax, a0, qDist, pDist, &
       sigmaAbs, sigmaSca, grainType)

  !write(*,*) amin,amax,qdist
  !do i = 1, grid%nLambda
  !   write(*,*) grid%lamArray(i),sigmaAbs(i),sigmaSca(i)
  !enddo

  if (grid%cartesian) then

     do i = 1, grid%nx
        do j = 1, grid%ny
           do k = 1, grid%nz


              if (grid%inUse(i,j,k)) then
                 grid%kappaAbs(i,j,k,1:grid%nLambda) = sigmaAbs  * grid%rho(i,j,k)
                 grid%kappaSca(i,j,k,1:grid%nLambda) = sigmaSca  * grid%rho(i,j,k)
              endif

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


  grid%kappaAbs = grid%kappaAbs * 1.e10
  grid%kappaSca = grid%kappaSca * 1.e10

  write(*,'(a)') "mie cross-sections done. Note 10^10 factor"
end subroutine fillGridMie2

