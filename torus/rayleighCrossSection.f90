subroutine rayleighCrossSection(lamArray, kappaSca, kappaExt, nLambda)

 implicit none
 integer :: nLambda
 real :: sigma
 real :: lamArray(1:nLambda)
 real :: kappaSca(1:nLambda)
 real :: kappaExt(1:nLambda)
 integer :: i

 do i = 1, nLambda
   sigma = 0.66250e-24 * (1026./lamArray(i))
   kappaSca(i) = sigma
   kappaExt(i) = kappaExt(i) + kappaSca(i)
 enddo

end subroutine rayleighCrossSection