module ostar_mod
  use constants_mod
  use vector_mod
  use gridtype_mod
 
  implicit none
  

contains

  function spiralWindDensity(rVec, grid) result (rhoOut)
    use inputs_mod
    real :: v, r, rhoOut
    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec
    
    r = real(modulus(rVec))
    if (r < grid%rCore) then
       rhoOut = 1.e-30
    else
       v = v0 + (vTerm - v0) * (1. - grid%rCore/r)**beta
       rhoOut = real(mdot / (fourPi * (r*1.e10)**2 * v) * returnSpiralFactor(rVec, 10.*grid%rCore/real(twoPi)))
    endif

  end function spiralWindDensity


  function returnSpiralFactor(rVec, alpha) result (fac)
    real :: alpha, fac, r, rSpiral1, rSpiral2, theta
    real ::  mu, x
    integer :: n
    type(VECTOR) :: rVec
    fac  = 1.

    r = real(modulus(rVec))
    theta = real(atan2(rVec%y, rVec%x))
    if (theta < 0.) theta = theta + real(twoPi)
    mu = real(rVec%z/r)
    if (abs(mu) < 0.5) then
       n = 0
       rSpiral2 = 0.
       do while (rSpiral2 < r)
          n = n + 1
          rSpiral2 = real(alpha*(theta+(real(n)*twoPi)))
       enddo
       rSpiral1= real(alpha*(theta+(real(n-1)*twoPi)))
       x = (r - rSpiral1)/(rSpiral2 - rSpiral1)
       
       if (((x > 0.).and.(x < 0.1)).or.(x > 0.9)) then
          fac = 100.
       endif
    endif
  end function returnSpiralFactor

    
    

    

end module ostar_mod
