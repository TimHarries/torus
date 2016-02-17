module viscosity_mod
#ifdef MPI
  
  use messages_mod
  use vector_mod
  use amr_mod
  use mpi_amr_mod
  use mpi_global_mod
  implicit none
  
contains
  
  real(double) function dudx(thisOctal, subcell, dir_u, dir_x, grid, divR)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: u_i_minus_1, u_i_plus_1, u_i
    type(VECTOR) :: cen, dx, cen_i_plus_1, cen_i_minus_1
    type(VECTOR) :: dir_u, dir_x, locator, tmpV
    real(double) :: q, rho_m1,rho_p1, rhoe, rhou, rhov, rhow, x, xnext, qnext, pressure, flux, phi, phigas, px, py, pz, qViscosity(3,3)
    real(double) :: rm1, um1, pm1, r, correction
    integer :: nd, nc
    logical ::divideByR
    logical, optional :: divR

    divideByR=.false.
    if (present(divR)) divideByR=divR
    
    
    cen=subcellCentre(thisOctal,subcell)
    locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir_x, q, rho_m1, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (abs(rhov) < 1.d-30) rhov = sign(1.d-30,rhov)
    cen_i_minus_1 = VECTOR(px,py,pz)
    if (cylindricalHydro) then
       r = sqrt(px**2 + py**2)*gridDistanceScale
       u_i_minus_1 = 0.d0
       if (rho_m1 /= 0.d0) then
          tmpV = VECTOR(rhou,rhov/r,rhow)
       else
          tmpV = VECTOR(rhou,0.0,rhow)
       endif
!       r = sqrt(cen%x**2 + cen%y**2)*gridDistanceScale
!       if (r/=0) then
!          tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
!                               thisOctal%rhov(subcell)/r,&
!                               thisOctal%rhow(subcell))
!       else
!          tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
!                               0.0,&
!                               thisOctal%rhow(subcell))
!       endif
!       u_i_minus_1 = (tmpV / (rho_m1+thisOctal%rho(subcell))).dot.dir_u
       u_i_minus_1 = (tmpV / rho_m1).dot.dir_u
       if (divideByR) u_i_minus_1 = u_i_minus_1 /r
    else
       u_i_minus_1 = 0.d0
       if (rho_m1 /= 0.d0) &
            tmpV = VECTOR(rhou,rhov,rhow)
       tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
            thisOctal%rhov(subcell),&
            thisOctal%rhow(subcell))
       u_i_minus_1 = (tmpV / (rho_m1+thisOctal%rho(subcell))).dot.dir_u
       if (divideByR) u_i_minus_1 = u_i_minus_1 /r
    endif
    
    
    locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir_x, q, rho_p1, rhoe, & 
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc,xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (abs(rhov) < 1.d-30) rhov = sign(1.d-30,rhov)
    cen_i_plus_1 = VECTOR(px,py,pz)
    if (cylindricalHydro) then
       r = sqrt(px**2 + py**2)*gridDistanceScale
       u_i_plus_1 = 0.d0
       if (rho_p1 /= 0.d0) then
          tmpV = VECTOR(rhou,rhov/r,rhow)
       else
          tmpV = VECTOR(rhou,0.0,rhow)
       endif
!       r = sqrt(cen%x**2 + cen%y**2)*gridDistanceScale
!       if (r/=0) then
!          tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
!               thisOctal%rhov(subcell)/r,&
!               thisOctal%rhow(subcell))
!       else
!          tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
!               0.0,&
!               thisOctal%rhow(subcell))
!       endif
!       u_i_plus_1 = (tmpV / (rho_p1+thisOctal%rho(subcell))).dot.dir_u
       u_i_plus_1 = (tmpV / rho_p1).dot.dir_u
       if (divideByR) u_i_plus_1 = u_i_plus_1 /r
    else
       u_i_plus_1 = 0.d0
       if (rho_p1 /= 0.d0) &
            tmpV = VECTOR(rhou,rhov,rhow)
       tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
            thisOctal%rhov(subcell),&
            thisOctal%rhow(subcell))
       u_i_plus_1 = (tmpV / (rho_p1+thisOctal%rho(subcell))).dot.dir_u
       if (divideByR) u_i_plus_1 = u_i_plus_1 /r
    endif

    
    if (cylindricalHydro) then
       tmpV = VECTOR(thisOctal%rhou(subcell),&
            thisOctal%rhov(subcell)/r,&
            thisOctal%rhow(subcell))
    else
       tmpV = VECTOR(thisOctal%rhou(subcell),&
            thisOctal%rhov(subcell),&
            thisOctal%rhow(subcell))
    endif
    
    tmpV = tmpV * (1/thisOctal%rho(subcell))
    u_i = tmpV .dot. dir_u
    
    
!    dx = cen_i_plus_1 - cen
!    dudx =        rho_p1  * (u_i_plus_1 - u_i)  / abs((dx.dot.dir_x)*gridDistanceScale/2) ! /2 as we are now using face based velocities for the cell so the distace is ((cen_i+1 + cen_i)/2 - (cen_i-1 + cen_i)/2) = (cen_i+1 - cen_i-1)/2
!    dx = cen - cen_i_minus_1
!    dudx = dudx + rho_m1  * (u_i - u_i_minus_1) / abs((dx.dot.dir_x)*gridDistanceScale/2) ! /2 as we are now using face based velocities for the cell so the distace is ((cen_i+1 + cen_i)/2 - (cen_i-1 + cen_i)/2) = (cen_i+1 - cen_i-1)/2
!
!    dudx = dudx / (rho_p1+rho_m1)

    dx=cen_i_plus_1 - cen_i_minus_1
    dudx = (u_i_plus_1 -u_i_minus_1) / abs(dx.dot.dir_x * gridDistanceScale)
    
  end function dudx
  
  real(double) function ddudxx(thisOctal, subcell, dir_u, dir_x, grid, divR)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: u_i_minus_1, u_i_plus_1, u_i
    type(VECTOR) :: cen, dx, cen_i_plus_1, cen_i_minus_1
    type(VECTOR) :: dir_u, dir_x, locator, tmpV
    real(double) :: q, rho_m1,rho_p1, rhoe, rhou, rhov, rhow, x, xnext, qnext, pressure, flux, phi, phigas, px, py, pz, qViscosity(3,3)
    real(double) :: rm1, um1, pm1, r, correction, q1,q2
    integer :: nd, nc
    logical :: divideByR
    logical, optional :: divR

    divideByR=.false.
    if (present(divR)) divideByR=divR    
    
    cen = subcellCentre(thisOctal,subcell)
    r=sqrt(cen%x**2+cen%y**2)*GridDistanceScale
    
    if (cylindricalHydro) then
       tmpV = VECTOR(thisOctal%rhou(subcell),&
                     thisOctal%rhov(subcell)/r,&
                     thisOctal%rhow(subcell))
    else
       tmpV = VECTOR(thisOctal%rhou(subcell),&
                     thisOctal%rhov(subcell),&
                     thisOctal%rhow(subcell))
    endif
    
    tmpV = tmpV / thisOctal%rho(subcell)
    u_i = tmpV .dot. dir_u
    if (divideByR) u_i = u_i /r

    locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir_x, q, rho_m1, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (abs(rhov) < 1.d-30) rhov = sign(1.d-30,rhov)
    cen_i_minus_1 = VECTOR(px,py,pz)
    if (cylindricalHydro) then
       u_i_minus_1 = 0.d0         
       if (rho_m1 /= 0.d0) then
          r = sqrt(px**2 + py**2)*gridDistanceScale
          tmpV = VECTOR(rhou,rhov/r,rhow)
          u_i_minus_1 = (tmpV / rho_m1).dot.dir_u
          if (divideByR) u_i_minus_1 = u_i_minus_1 /r
       endif
    else
       u_i_minus_1 = 0.d0
       if (rho_m1 /= 0.d0) then
          tmpV = VECTOR(rhou,rhov,rhow)
          u_i_minus_1 = (tmpV / rho_m1).dot.dir_u
          if (divideByR) u_i_minus_1 = u_i_minus_1 /r
       endif
    endif
    
    q1=qViscosity(1,1)
    
    locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir_x, q, rho_p1, rhoe, & 
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc,xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (abs(rhov) < 1.d-30) rhov = sign(1.d-30,rhov)
    cen_i_plus_1 = VECTOR(px,py,pz)
    if (cylindricalHydro) then
       u_i_plus_1 = 0.d0
       if (rho_p1 /= 0.d0) then
          r = sqrt(px**2 + py**2)*gridDistanceScale
          tmpV = VECTOR(rhou,rhov/r,rhow)
          u_i_plus_1 = (tmpV / rho_p1).dot.dir_u
          if (divideByR) u_i_plus_1 = u_i_plus_1 /r
       endif
    else
       u_i_plus_1 = 0.d0
       if (rho_p1 /= 0.d0) then
          tmpV = VECTOR(rhou,rhov,rhow)
          u_i_plus_1 = (tmpV / rho_p1).dot.dir_u
          if (divideByR) u_i_plus_1 = u_i_plus_1 /r
       endif
    endif

    q2=qViscosity(1,1)

    
    dx = (cen_i_plus_1 - cen_i_minus_1)
!    print *, u_i_minus_1 , u_i, u_i_plus_1
!    print *, q1, thisOctal%qViscosity(subcell,1,1), q2
    ddudxx =(u_i_plus_1 + u_i_minus_1 - 2.0*u_i) / ((dx.dot.dir_x)*gridDistanceScale/2)**2 ! /2 as dx . dir_x  is twice the distance between 2 cells
 !   print *, ddudxx
  end function ddudxx
  
  real(double) function ddudxdz(thisOctal, subcell, dir_u, grid)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: dudx_m1, dudx_p1
    type(VECTOR) :: cen,  cen_i_plus_1, cen_i_minus_1, locator, dir_u
    logical      :: problem
    
    cen = subcellCentre(thisOctal,subcell)
    if ((cen%x - (3*thisOctal%subcellSize/2.d0+0.01d0*smallestCellSize)) > 0.d0) then ! not on axis
       locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
       neighbouroctal => thisoctal
       call findsubcelllocal(locator, neighbouroctal, neighboursubcell, problem)
       if (problem) then
          dudx_m1 = 0
       else
          dudx_m1 = dudx(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0, 0.0), grid)
       endif
       cen_i_minus_1= subcellCentre(neighbourOctal, neighbourSubcell)
       
       locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
       neighbouroctal => thisoctal
       call findsubcelllocal(locator, neighbouroctal, neighboursubcell, problem)
       if (problem) then
          dudx_p1 = 0
       else
          dudx_p1 = dudx(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0,0.0), grid)
       endif
       cen_i_plus_1= subcellCentre(neighbourOctal, neighbourSubcell)
       
       ddudxdz = (dudx_p1 - dudx_m1) / abs((cen_i_plus_1-cen_i_minus_1).dot.vector(0.0, 0.0, 1.0)*gridDistanceScale)
    else
       ddudxdz = 0.0
    endif
  end function ddudxdz
  
  real(double) function div_u_cylindrical(thisOctal, subcell, grid)
    use inputs_mod, only : gridDistanceScale
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    type(VECTOR) :: cen
    real(double) :: r
    
    cen = subcellCentre(thisOctal,subcell)
    r= sqrt(cen%x**2 +cen%y**2)* gridDistanceScale
    div_u_cylindrical= 0.d0
    div_u_cylindrical = div_u_cylindrical + dudx(thisOctal, subcell, vector(1.0, 0.0, 0.0), vector(1.0, 0.0, 0.0), grid)
    div_u_cylindrical = div_u_cylindrical + thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*r)
    div_u_cylindrical = div_u_cylindrical + dudx(thisOctal, subcell, vector(0.0, 0.0, 1.0), vector(0.0, 0.0, 1.0), grid)
    print *, "div_u: ", div_u_cylindrical
  end function div_u_cylindrical

  real(double) function diff_div_u_cylindrical(thisOctal, subcell, dir_x, grid)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: divu_m1, divu_p1
    type(VECTOR) :: cen,  cen_i_plus_1, cen_i_minus_1, locator, dir_x
    
    cen = subcellCentre(thisOctal,subcell)
    locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    divU_m1 = div_u_cylindrical(neighbourOctal,neighbourSubcell, grid)
    cen_i_minus_1= subcellCentre(neighbourOctal, neighbourSubcell)
    
    locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    divU_p1 = div_u_cylindrical(neighbourOctal,neighbourSubcell, grid)
    cen_i_plus_1= subcellCentre(neighbourOctal, neighbourSubcell)
    
    diff_div_u_cylindrical = (divU_p1 - divU_m1) / abs((cen_i_plus_1-cen_i_minus_1).dot.dir_x*gridDistanceScale)
    print *,divU_p1, divU_m1
    print *, diff_div_u_cylindrical, ((cen_i_plus_1-cen_i_minus_1).dot.dir_x*gridDistanceScale)
    
  end function diff_div_u_cylindrical
  
  real(double) function div_u(thisOctal, subcell, grid)
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal
    integer :: subcell, i, nDim
    type(VECTOR) :: dir(3)
    if (thisOctal%threed) then
       ndim = 3
       dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
       dir(2) = VECTOR(0.d0, 1.d0, 0.d0)
       dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
    else
       ndim = 2
       dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
       dir(2) = VECTOR(0.d0, 0.d0, 1.d0)
    endif
    
    div_u = 0.d0
    do i = 1, nDim
       div_u = div_u + dudx(thisOctal, subcell, dir(i), dir(i), grid)
    enddo
  end function div_u
  
  function symmetricVelocityGradientTensor(thisOctal, subcell, grid) result(out) !strain rate tensor
    type(OCTAL), pointer :: thisOctal
    integer :: subcell
    real(double) :: out(3,3), temp(3,3)
    type(GRIDTYPE) :: grid
    integer :: i, j, nDim
    type(VECTOR) :: dir(3)
    if (thisOctal%threed) then
       ndim = 3
       dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
       dir(2) = VECTOR(0.d0, 1.d0, 0.d0)
       dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
    else
       ndim = 2
       dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
       dir(2) = VECTOR(0.d0, 0.d0, 1.d0)
    endif
    
    
    
    do i = 1, nDim
       do j = 1, nDim
          temp(i,j)  = dudx(thisOctal, subcell, dir(i), dir(j), grid)
       enddo
    enddo
    do i = 1, nDim
       do j = 1, nDim
          out(i,j)  = 0.5d0*(temp(i,j) + temp(j,i))
       enddo
    enddo
  end function symmetricVelocityGradientTensor
  
  function viscosityQ(thisOctal, subcell, grid) result(out)
    use inputs_mod, only : etaViscosity, smallestCellSize, gridDistanceScale
    type(OCTAL), pointer :: thisoctal
    integer :: subcell
    type(GRIDTYPE) :: grid
    real(double) :: out(3,3), symVelGrad(3,3), temp(3,3), divV, lengthScale
    integer :: i,nDim
    
    if (thisOctal%threed) then
       ndim = 3
    else
       nDim = 2
    endif
    
    lengthScale = etaViscosity * smallestCellSize * gridDistanceScale
    out = 0.d0
    divV = div_u(thisOctal, subcell, grid)
    if (divV < 0.d0) then
       symVelGrad = symmetricVelocityGradientTensor(thisOctal, subcell, grid)
       temp = symVelGrad
       do i = 1, nDim
          temp(i,i) = temp(i,i) - (1.d0/3.d0)*divV
       enddo
       out = (lengthScale**2 * thisOctal%rho(subcell) * divV)*temp
    endif
  end function viscosityQ
  
  function qColonDelV(thisOctal, subcell, grid) result(out)
    use inputs_mod, only : smallestCellsize,etaViscosity, gridDistanceScale
    type(octal), pointer :: thisOctal
    integer :: subcell
    type(GRIDTYPE) :: grid
    real(double) :: out, symVelGrad(3,3), lengthScale, divV
    
    lengthScale = etaViscosity * smallestCellSize * gridDistanceScale
    out = 0.d0
    divV = div_u(thisOctal, subcell, grid)
    symVelGrad = symmetricVelocityGradientTensor(thisOctal, subcell, grid)
    
    out = 0.d0
    if (thisOctal%threed) then
       out = lengthScale**2 * thisOctal%rho(subcell) * divV * &
            ((symVelGrad(1,1)-symVelGrad(2,2))**2 + &
            (symVelGrad(1,1)-symVelGrad(3,3))**2 + &
            (symVelGrad(2,2)-symVelGrad(3,3))**2)/3.d0
    else
       out = lengthScale**2 * thisOctal%rho(subcell) * divV * &
            ((symVelGrad(1,1)-symVelGrad(2,2))**2)/2.d0
    endif
    
  end function qColonDelV
  
  
  
  function divQ(thisOctal, subcell,  grid) result(out)
    use inputs_mod, only : smallestCellSize,gridDistanceScale
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    type(VECTOR) :: out, cen_i_plus_1, cen_i_minus_1
    real(double) :: q, rho, rhoe, rhou,rhov,rhow, x, qnext, pressure, flux, phi, phigas,xnext,px,py,pz,dx
    real(double) :: qViscosity1(3,3), qViscosity2(3,3),correction
    integer :: subcell, neighbourSubcell
    type(VECTOR) :: dir(3), cen2, locator
    real(double) :: tmp(3,3), tmp2(3)
    real(double) :: rm1, um1, pm1, r
    integer :: nd, iDir, jDir,nc
    dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
    dir(2) = VECTOR(0.d0, 1.d0, 0.d0)
    dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
    
    
    out = VECTOR(0.d0,0.d0,0.d0)
    cen2 = subcellCentre(thisOctal,subcell)
    locator = cen2 - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir(1)
    
    r = sqrt(cen2%x**2 + cen2%y**2) * gridDistanceScale
    
    do iDir = 1, 3
       do jDir = 1, 3
          
          
          
          if (jDir /= 2) then
             locator = cen2 + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir(jDir)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir(jDir), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, &
                  qViscosity1)
             cen_i_plus_1 = VECTOR(px, py,pz)
             
             locator = cen2 - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir(jDir)
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir(jDir), q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, &
                  qViscosity2)
             cen_i_minus_1 = VECTOR(px, py,pz)
             
             dx = (cen_i_plus_1 - cen_i_minus_1).dot.dir(jDir)
             tmp(iDir, jDir) = (qviscosity1(iDir, jDir) - qviscosity2(iDir, jDir))/(dx*gridDistanceScale)
          else
             
             select case(idir)
             case(1)
                tmp(iDir, jDir) = (1.d0/r)*(thisOctal%qViscosity(subcell,1,1)-thisOctal%qViscosity(subcell,2,2))
             case(2)
                tmp(iDir, jDir) = (1.d0/r)*(thisOctal%qViscosity(subcell,1,2)+thisOctal%qViscosity(subcell,2,1))
             case(3)
                tmp(iDir, jDir) = (1.d0/r)*(thisOctal%qViscosity(subcell,1,3))
             end select
          endif
          
       enddo
    enddo
    do iDir = 1, 3
       tmp2(iDir) = SUM(tmp(iDir,1:3))
    enddo
    
    out%x = tmp2(1)
    out%y = tmp2(2)
    out%z = tmp2(3)
    !666 continue
  end function divQ
  
  function newdivQ(thisOctal, subcell,  grid) result(out)
    use inputs_mod, only : smallestCellSize,gridDistanceScale
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    type(VECTOR) :: out
    integer :: subcell, neighbourSubcell
    type(VECTOR) :: dir, cen2, cen, locator
    real(double) :: r, oldq
    logical :: flag
    
    out = VECTOR(0.d0, 0.d0, 0.d0)
    
    dir = VECTOR(1.d0, 0.d0, 0.d0)
    cen2 = subcellCentre(thisOctal,subcell)
    
    r = cen2%x * gridDistanceScale
    
    
    if ( (cen2%x - (thisOctal%subcellSize/2.d0+0.01d0*smallestCellSize)) > 0.d0) then ! not on axis
       out%x = out%x + dqdx(thisOctal, subcell, grid, 1, 1, VECTOR(1.d0, 0.d0, 0.d0))
       out%x = out%x + thisOctal%qViscosity(subcell,1,1) / r
       out%x = out%x - thisOctal%qViscosity(subcell,2,2) / r
       out%x = out%x + dqdx(thisOctal, subcell, grid, 1, 3, VECTOR(0.d0, 0.d0, 1.d0))
       
       out%y = out%y  + dqdx(thisOctal, subcell, grid, 1, 2, VECTOR(1.d0, 0.d0, 0.d0))
       out%y = out%y + 2.d0 * thisOctal%qViscosity(subcell,1,2)/r
       out%y = out%y  + dqdx(thisOctal, subcell, grid, 2, 3, VECTOR(0.d0, 0.d0, 1.d0))
       
       out%z = out%z + dqdx(thisOctal, subcell, grid, 1, 3, VECTOR(1.d0, 0.d0, 0.d0))
       out%z = out%z + thisOctal%qViscosity(subcell,1,3)/r
       out%z = out%z + dqdx(thisOctal, subcell, grid, 3, 3, VECTOR(0.d0, 0.d0, 1.d0))
    endif

!    flag= associated(thisOctal%mpiBoundaryStorage)
!    if (.not. flag) then
!       cen=subcellCentre(thisOctal,subcell)
!       locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(1,0,0)
!       neighbouroctal => thisoctal
!       call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
!       flag =  associated(neighbourOctal%mpiBoundaryStorage)
!       if (.not. flag) then
!          cen=subcellCentre(thisOctal,subcell)
!          locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(1,0,0)
!          neighbouroctal => thisoctal
!          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
!          flag =  associated(neighbourOctal%mpiBoundaryStorage)
!          if (.not. flag) then
!             cen=subcellCentre(thisOctal,subcell)
!             locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0,0,1)
!             neighbouroctal => thisoctal
!             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
!             flag =  associated(neighbourOctal%mpiBoundaryStorage)
!             if (.not. flag) then
!                cen=subcellCentre(thisOctal,subcell)
!                locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0,0,1)
!                neighbouroctal => thisoctal
!                call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
!                flag =  associated(neighbourOctal%mpiBoundaryStorage)
!             endif
!          endif
!       endif
!    endif
    
    if ((cen2%z >3.0d6 ))then !.or. cen2%z<-3.0d6) ) then !(abs(cen2%x-8.0d7)>25.0e6 .and.(cen2%x>3.0e6 .or. cen2%x<1.5e8))) then
       oldq=out%y
       out%x=thisOctal%etaline(subcell)*(4.0/3*ddudxx(thisOctal, subcell, vector(1.0,0,0),vector(1.0,0,0),grid)+&
                                         1.0/3*ddudxx(thisOctal, subcell, vector(0,0,1.0),vector(0,0,1.0),grid)+&
                                         ddudxx(thisOctal, subcell, vector(1.0,0,0),vector(0,0,1.0),grid)+&
                                         4.0/3/r*dudx(thisOctal, subcell, vector(1.0,0,0), vector(1.0,0,0),grid)-&
                                         4.0*3/r/r *thisOctal%rhou(subcell)/thisOctal%rho(subcell))
       out%y=thisOctal%etaline(subcell)*(ddudxx(thisOctal, subcell, vector(0,1.0,0),vector(1.0,0,0),grid)+&
                                         ddudxx(thisOctal, subcell, vector(0,1.0,0),vector(0,0,1.0),grid)+&
                                         dudx(thisOctal,subcell,vector(0,1.0,0),vector(1.0,0,0),grid)/r-&
                                         thisOctal%rhov(subcell)/thisOctal%rho(subcell)/r**3)
       out%z=thisOctal%etaline(subcell)*(ddudxx(thisOctal, subcell, vector(0,0,1.0),vector(1.0,0,0),grid)+&
                                         1.0/3*ddudxdz(thisOctal, subcell, vector(1.0,0,0),grid)+&
                                         4.0/3*ddudxx(thisOctal, subcell, vector(0,0,1.0), vector(0,0,1.0), grid)+&
                                         1.0/r*dudx(thisOctal, subcell, vector(0,0,1.0), vector(1.0,0,0), grid)+&
                                         1.0/3/r*dudx(thisOctal, subcell, vector(1.0,0,0), vector(0,0,1.0), grid))

       thisOctal%qViscosity(subcell,1,1)=ddudxx(thisOctal, subcell, vector(0,1.0,0),vector(1.0,0,0),grid)
       thisOctal%qViscosity(subcell,2,2)=thisOctal%rhov(subcell)/thisOctal%rho(subcell)/r**3
       
       
!       if (1==myRankGlobal) then
!          print *, "^^^^^^^^^^"
!          print *, "old:", oldq
!          !       print *, "old comp:", thisOctal%qViscosity(subcell,1,1)/r-&
!          !                             thisOctal%qViscosity(subcell,2,2)/r,&
!          !                             dqdx(thisOctal, subcell, grid, 1, 1, VECTOR(1.d0, 0.d0, 0.d0)),&
!          !                             dqdx(thisOctal, subcell, grid, 1, 3, VECTOR(0.d0, 0.d0, 1.d0))
!          print *, "new:", out
!          print *, "eta:", thisOctal%etaline(subcell)
!          print *, "new comp (theta):", r*ddudxx(thisOctal, subcell, vector(0,1.0,0),vector(1.0,0,0),grid, .true.),&
!                                        ddudxx(thisOctal, subcell, vector(0,0,1.0),vector(0,0,1.0),grid),&
!                                        3*dudx(thisOctal, subcell, vector(0,1.0,0),vector(1.0,0,0),grid, .true.)
!          
!          !      thisOctal%etaline(subcell)*2*ddudxx(thisOctal, subcell, vector(1,0,0),vector(1,0,0),grid)
!          !print *, "new comp:",thisOctal%qViscosity(subcell,1,1)/r-&
!          !       !                     thisOctal%qViscosity(subcell,2,2)/r,&
!          !       !                     thisOctal%etaline(subcell)*(2*ddudxx(thisOctal, subcell, vector(1,0,0),vector(1,0,0),grid)),&
!          !       !                     thisOctal%etaline(subcell)*ddudxdz(thisOctal, subcell, vector(0,0,1),grid),&
!          !       !                     thisOctal%etaline(subcell)*ddudxx(thisOctal, subcell, vector(1,0,0),vector(0,0,1),grid),&
!          !       !                     thisOctal%etaline(subcell)*2/3*diff_div_u_cylindrical(thisOctal, subcell, vector(1,0,0),grid)
!          !       if  (oldq /=0) print *, "ratio:", out%y/oldq
!          !       print *, "vvvvvvvvvv"
!       endif
    endif
    
  end function newdivQ
  
  function dQdx(thisOctal, subcell,  grid, i, j, dir, multbyr) 
    use inputs_mod, only : smallestCellSize,gridDistanceScale
    type(gridtype) :: grid
    real(double) :: dqdx
    logical, optional :: multbyr
    logical :: doMultByR
    type(octal), pointer   :: thisoctal, neighbourOctal
    type(VECTOR) :: cen_i_plus_1, cen_i_minus_1
    real(double) :: q, rho, rhoe, rhou,rhov,rhow, x, qnext, pressure, flux, phi, phigas,xnext,px,py,pz,dx
    real(double) :: qViscosityM1(3,3), qViscosityM2(3,3), qViscosityP1(3,3), qViscosityP2(3,3), correction
    integer :: subcell, neighbourSubcell
    type(VECTOR) :: dir, cen2, locator
    real(double) :: rm1, um1, pm1, r,r1 ,r2
    integer :: nd, i, j, nc
    
    doMultByR = .false.
    if (present(multbyr)) doMultByR = multByR
    
    dqdx = 0.d0
    
    cen2 = subcellCentre(thisOctal,subcell)
    r = sqrt(cen2%x**2 + cen2%y**2) * gridDistanceScale
    
    
    locator = cen2 + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, qViscosityP1)
    cen_i_plus_1 = VECTOR(px, py,pz)
    r2 = abs(px)*gridDistanceScale
    
    locator = cen_i_plus_1 + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, qViscosityP2)
    
    locator = cen2 - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, qViscosityM1)
    cen_i_minus_1 = VECTOR(px, py,pz)
    r1 = abs(px)*gridDistanceScale
    
    locator = cen_i_minus_1 - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, qViscosityM2)
    
    dx = (cen_i_plus_1 - cen_i_minus_1).dot.dir
    if (doMultByR) then
       dqdx = (r2*qviscosityP1(i, j) - r1*qviscosityM1(i, j))/(dx*gridDistanceScale)
    else
       dqdx = (qViscosityM2(i,j)/6  - qViscosityM1(i,j)*4/3 +&
            qViscosityP1(i,j)*4/3 - qViscosityP2(i,j)/6)/(dx*gridDistanceScale)
    endif !changed dqdx to be 4th order accurate rather than 2nd, should stop checkerboarding of viscous forces 
  end function dQdx
  
  function divV(thisOctal, subcell, grid) result(out)
    use inputs_mod, only : smallestCellSize,gridDistanceScale
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    real(double) :: out
    real(double) :: q, rho, rhoe, rhou,rhov,rhow, x, qnext, pressure, flux, phi, phigas,xnext,px,py,pz,qViscosity(3,3)
    real(double) :: speedPlus, speedMinus, correction
    integer :: subcell, neighbourSubcell
    type(VECTOR) :: dir(3), cen, locator
    integer :: iDir, nDir, nd, nc
    real(double) :: rm1, um1, pm1
    if (thisOctal%threed) then
       ndir = 3
       dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
       dir(2) = VECTOR(0.d0, 1.d0, 0.d0)
       dir(3) = VECTOR(0.d0, 0.d0, 1.d0)
    else if (thisOctal%twoD) then
       ndir = 2
       dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
       dir(2) = VECTOR(0.d0, 0.d0, 1.d0)
    else
       ndir = 1
       dir(1) = VECTOR(1.d0, 0.d0, 0.d0)
    endif
    
    cen = subcellCentre(thisOctal,subcell)
    out = 0.d0
    
    if (.not.thisOctal%edgeCell(subcell)) then
       do iDir = 1, nDir
          
          locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir(iDir)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir(iDir), q, rho, rhoe, &
               rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, &
               qViscosity)
          if (thisOctal%threed) then
             select case (iDir)
             case(1)
                speedplus = rhou/rho
             case(2)
                speedplus = rhov/rho
             case(3)
                speedplus = rhow/rho
             end select
          else if (thisOctal%twoD) then
             select case (iDir)
             case(1)
                speedplus = rhou/rho
             case(2)
                speedplus = rhow/rho
             end select
          else
             speedPlus = rhou/rho
          endif
          
          locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir(iDir)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir(iDir), q, rho, rhoe, &
               rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, px, py, pz, rm1, um1, pm1, &
               qViscosity)
          if (thisOctal%threed) then
             select case (iDir)
             case(1)
                speedminus = rhou/rho
             case(2)
                speedminus = rhov/rho
             case(3)
                speedminus = rhow/rho
             end select
          else if (thisOctal%twoD) then
             select case (iDir)
             case(1)
                speedminus = rhou/rho
             case(2)
                speedminus = rhow/rho
             end select
          else
             speedminus = rhou/rho
          endif
          
          out = out + (speedplus - speedminus) / (thisOctal%subcellSize*gridDistanceScale)
       end do
    endif
  end function divV
  
  
  recursive subroutine setupViscosity(thisoctal, grid)
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupViscosity(child, grid)
                exit
             end if
          end do
       else
          
          if (.not.associated(thisOctal%chiline)) then
             allocate(thisOctal%chiline(1:thisOctal%maxChildren))
             thisOctal%chiline = 1.d-30
          endif
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle
          thisoctal%qViscosity(subcell,1:3,1:3) = 0.d0
          
          if (.not.thisOctal%edgeCell(subcell)) then
             thisoctal%qViscosity(subcell,1:3,1:3) = viscosityQ(thisOctal, subcell, grid)
             thisOctal%chiLine(subcell) = qColonDelV(thisOctal, subcell, grid)
          endif
       endif
    enddo
  end subroutine setupViscosity
  
  recursive subroutine setupCylindricalViscosity(thisoctal, grid)
    use inputs_mod, only : gridDistanceScale!, smallestCellSize
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child
    real(double) :: divV, r, vTheta, fac
    integer :: subcell, i
    type(VECTOR) :: rVec
    
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupCylindricalViscosity(child, grid)
                exit
             end if
          end do
       else
          
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle
          thisoctal%qViscosity(subcell,1:3,1:3) = 0.d0
          
          if (.not.thisOctal%edgeCell(subcell)) then
             
             ! first calculate div V
             
             rVec = subcellCentre(thisOctal, subcell)
             r = sqrt(rVec%x**2+rVec%y**2)*gridDistanceScale
             
             divV = dudx(thisOctal, subcell, VECTOR(1.d0, 0.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid) + &
                  thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*r) +  & 
                  dudx(thisOctal, subcell, VECTOR(0.d0, 0.d0, 1.d0), VECTOR(0.d0, 0.d0, 1.d0), grid)
             
             !             if (thisOctal%rho(subcell) > 1.d-15) then
             !                write(*,*) "dvphi/dr ", &
             !                     dudx(thisOctal, subcell, VECTOR(0.d0, 1.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid), &
             !                     -0.5d0*sqrt(bigG * mSol) * r**(-1.5d0)
             !             endif
             
             
             
             ! now tau_rr
             
             thisOctal%qViscosity(subcell,1,1) = thisOctal%etaline(subcell) * &
                  (2.d0 * dudx(thisOctal, subcell, VECTOR(1.d0, 0.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid) - 0.6666666666d0 * divV)
             
             ! now tau_thetatheta
             thisOctal%qViscosity(subcell,2,2) =  thisOctal%etaline(subcell) * &
                  (2.d0 * thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*r) - 0.6666666666d0 * divV)
             
             ! now tau_zz
             
             thisOctal%qViscosity(subcell,3,3) = thisOctal%etaline(subcell) * &
                  (2.d0 * dudx(thisOctal, subcell, VECTOR(0.d0, 0.d0, 1.d0), VECTOR(0.d0, 0.d0, 1.d0), grid) - 0.6666666666d0 * divV)
             
             ! now tau_rtheta
             
             vTheta = thisOctal%rhov(subcell) / thisOctal%rho(subcell) / r
             thisOctal%qViscosity(subcell,1,2) = thisOctal%etaline(subcell) * &
                  (dudx(thisOctal, subcell, VECTOR(0.d0, 1.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid) - (vTheta/r))
             
             ! now tau_thetar
             
             thisOctal%qViscosity(subcell,2,1) = thisOctal%qViscosity(subcell,1,2) 
             
             ! now tau_thetaz
             
             thisOctal%qViscosity(subcell,2,3) = thisOctal%etaline(subcell) * &
                  dudx(thisOctal, subcell, VECTOR(0.d0, 1.d0, 0.d0), VECTOR(0.d0, 0.d0, 1.d0), grid)
             
             
             ! now tau_thetaz
             
             thisOctal%qViscosity(subcell,3,2) = thisOctal%qViscosity(subcell,2,3)
             
             
             ! now tau_rz
             
             thisOctal%qViscosity(subcell,1,3) = thisOctal%etaline(subcell) * &
                  (dudx(thisOctal, subcell, VECTOR(0.d0, 0.d0, 1.d0), VECTOR(1.d0, 0.d0, 0.d0), grid) + &
                  dudx(thisOctal, subcell, VECTOR(1.d0, 0.d0, 0.d0), VECTOR(0.d0, 0.d0, 1.d0), grid))
             
             ! now tau_zr
             
             thisOctal%qViscosity(subcell,3,1) = thisOctal%qViscosity(subcell,1,3) 
             
             
             
             fac = thisOctal%etaLine(subcell) * r  * &
                  (dudx(thisOctal,subcell, VECTOR(0.d0, 1.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid)/r &
                  - (thisOctal%rhov(subcell)/(thisOctal%rho(subcell)*r**3)) )
             
             if (maxval(abs(thisOctal%qViscosity(subcell,:,:)))> 1) then
                write(*,*) "Q ",thisOctal%qViscosity(subcell,1,2), " mu r domegabydr ", fac, " ratio ", thisOctal%qViscosity(subcell,1,2)/fac
                write(*,*) "kep ",sqrt(bigG * msol) * (-3.d0/2.d0)*r**(-5.d0/2.d0), &
                     " model ",dudx(thisOctal,subcell, VECTOR(0.d0, 1.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid)/r
                fac = sqrt(bigG *mSol /r)
                write(*,*) "kep v ",fac, " model ",thisOctal%rhoV(subcell)/(thisOctal%rho(subcell)*r)
                write(*,*)  " "
                write(*,'(a,1p,3e9.1)') "qvisc ",thisOctal%qViscosity(subcell,1,1:3)
                write(*,'(a,1p,3e9.1)') "      ",thisOctal%qViscosity(subcell,2,1:3)
                write(*,'(a,1p,3e9.1)') "      ",thisOctal%qViscosity(subcell,3,1:3)
                write(*,*)  "Trace (%): ",(thisOctal%qViscosity(subcell,1,1) +  &
                     thisOctal%qViscosity(subcell,2,2) +  thisOctal%qViscosity(subcell,3,3) ) &
                     / maxval(abs(thisOctal%qViscosity(subcell,:,:)))
             endif
          endif
       endif
    enddo
  end subroutine setupCylindricalViscosity
  
  
  !"
  
  
  recursive subroutine viscousTimescale(thisoctal, grid, dt)
    use inputs_mod, only : etaViscosity, smallestCellsize, gridDistanceScale
    real(double) :: lengthScale
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    real(double) :: dt, thisTime, divV
    integer :: subcell, i
    
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call viscousTimescale(child, grid, dt)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle
          if (.not.thisOctal%ghostCell(subcell)) then
             lengthScale = max(1.d-20,etaViscosity * smallestCellSize * gridDistanceScale)
             divV = div_u(thisOctal, subcell, grid)
             thisTime = (smallestCellSize*gridDistanceScale)**2 /(4.d0* max(abs(divV),1.d-20) * lengthScale**2)
             dt = min(thisTime, dt)
          endif
       endif
    enddo
  end subroutine viscousTimescale
  
  recursive subroutine viscousTimescaleCylindrical(thisoctal, grid, dt)
    use inputs_mod, only : smallestCellsize, gridDistanceScale
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    real(double) :: dt, thisTime, acc, r, torque, cs
    type(VECTOR) :: fVisc,rvec
    integer :: subcell, i
    real(double) :: a,b,c
    
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call viscousTimescaleCylindrical(child, grid, dt)
                exit
             end if
          end do
       else
          if (.not.octalonthread(thisoctal, subcell, myrankglobal)) cycle
          if (.not.thisOctal%ghostCell(subcell)) then
             
             fVisc =  newdivQ(thisOctal, subcell,  grid)
             acc = max(abs(fVisc%x), abs(fVisc%y), abs(fvisc%z))/ thisOctal%rho(subcell)
             acc = max(acc, 1.d-30)
             
             thisTime = sqrt(smallestCellSize*gridDistanceScale/acc)*0.5
             dt = min(thisTime, dt)
             
!             rVec = subcellCentre(thisOCtal, subcell)            
!             r = max(sqrt(rVec%x**2+rVec%y**2) * gridDistanceScale, gridDistanceScale/2)
!             torque = abs(fVisc%y) * r
!             thisTime = max(0.5*abs(thisOctal%rhov(subcell))**2/(thisOctal%rho(subcell)*r**2),&
!                  thisOctal%rhoe(subcell)/thisOctal%rho(subcell))&
!                  / max(abs(torque),1.d-60)
!             cs=sqrt(thisOctal%gamma(subcell)*kErg*thisOctal%temperature(subcell)/mHydrogen)
!             a=thisTime
             b=dt
!             c=(abs(thisOctal%rhov(subcell))+cs*thisOctal%rho(subcell)*r/10)/&
!                  max(abs(torque), 1d-60)
!             dt=min(thisTime,dt)
!             dt = min(thisTime,&
!                  dt,&
!                  (abs(thisOctal%rhov(subcell))+cs*thisOctal%rho(subcell)*r/10)/&
!                  max(abs(torque), 1d-60)) !angular speed change is limited to current angular speed + r*soundSpeed/10
          endif
       
          if (dt<1) then
             print *, thisOctal%rhou(subcell)/thisOctal%rho(subcell),&
                  thisOctal%rhoV(subcell)/r/thisOctal%rho(subcell),&
                  thisOctal%rhow(subcell)/thisOctal%rho(subcell)
             PRINT *, "1", sqrt(smallestCellSize*gridDistanceScale/acc)*0.5
             print *, "2", a
             print *, "3", b
!             print *, "4", c
             
             print *, fVisc
             print *, torque
             print *, dt
             print *, "                     "
          endif
       endif
    enddo
  end subroutine viscousTimescaleCylindrical
  
#endif
end module viscosity_mod
