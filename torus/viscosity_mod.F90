module viscosity_mod
#ifdef MPI

  use messages_mod
  use vector_mod
  use amr_mod
  use mpi_amr_mod
  use mpi_global_mod
  implicit none

contains
  
  real(double) function dudx(thisOctal, subcell, dir_u, dir_x, grid)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: u_i_minus_1, u_i_plus_1
    type(VECTOR) :: cen
    type(VECTOR) :: dir_u, dir_x, locator
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, xnext, qnext, pressure, flux, phi, phigas, px, py, pz, qViscosity(3,3)
    real(double) :: rm1, um1, pm1, r
    integer :: nd
    
    cen = subcellCentre(thisOctal, subcell)    
    locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir_x, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1, um1, pm1, qViscosity)

    if (cylindricalHydro) then
       r = sqrt(px**2 + py**2)*gridDistanceScale
       u_i_minus_1 = (VECTOR(rhou,rhov/r,rhow).dot.dir_u)/rho
    else
       u_i_minus_1 = (VECTOR(rhou,rhov,rhow).dot.dir_u)/rho
    endif

    locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir_x, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (cylindricalHydro) then
       r = sqrt(px**2 + py**2)*gridDistanceScale
       u_i_plus_1 = (VECTOR(rhou,rhov/r,rhow).dot.dir_u)/rho
    else
       u_i_plus_1 = (VECTOR(rhou,rhov,rhow).dot.dir_u)/rho
    endif


    dudx = (u_i_plus_1 - u_i_minus_1) / (2.d0*thisOctal%subcellSize*gridDistanceScale)
  end function dudx

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

  function symmetricVelocityGradientTensor(thisOctal, subcell, grid) result(out)
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
    


  function divQ(thisOctal, subcell, iDir, grid) result(out)
    use inputs_mod, only : smallestCellSize,gridDistanceScale
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    real(double) :: out
    real(double) :: q, rho, rhoe, rhou,rhov,rhow, x, qnext, pressure, flux, phi, phigas,xnext,px,py,pz,qViscosity(3,3)
    integer :: subcell, neighbourSubcell
    type(VECTOR) :: dir(3), cen2, locator
    real(double) :: qa, qb
    real(double) :: rm1, um1, pm1
    integer :: nd, iDir, nDim
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

    out = 0.d0

    cen2 = subcellCentre(thisOctal,subcell)

    locator = cen2 + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir(iDir)
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir(iDir), q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1, um1, pm1, qViscosity)

    select case(iDir)
    case(1)
       qa = qViscosity(1,1)
    case(2)
       qa = qViscosity(2,2)
    case(3)
       qa = qViscosity(3,3)
    end select

    locator = cen2 - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir(iDir)
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir(iDir), q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1, um1, pm1, qViscosity)

    select case(iDir)
    case(1)
       qb = qViscosity(1,1)
    case(2)
       qb = qViscosity(2,2)
    case(3)
       qb = qViscosity(3,3)
    end select

    out = (qa - qb) / (2.d0*thisOctal%subcellSize*gridDistanceScale)
  end function divQ

  function divV(thisOctal, subcell, grid) result(out)
    use inputs_mod, only : smallestCellSize,gridDistanceScale
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal, neighbourOctal
    real(double) :: out
    real(double) :: q, rho, rhoe, rhou,rhov,rhow, x, qnext, pressure, flux, phi, phigas,xnext,px,py,pz,qViscosity(3,3)
    real(double) :: speedPlus, speedMinus
    integer :: subcell, neighbourSubcell
    type(VECTOR) :: dir(3), cen, locator
    integer :: iDir, nDir, nd
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
            rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1, um1, pm1, qViscosity)
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
            rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1, um1, pm1, qViscosity)
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
    use inputs_mod, only : gridDistanceScale, smallestCellSize
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child, neighbourOctal
    real(double) :: divV, r, vTheta
    integer :: subcell, i, neighbourSubcell, nd
    type(VECTOR) :: rVec, locator, cen, dir_x
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, xnext
    real(double) :: px, py, pz, q11, q22, q33, rm1, um1, pm1, u_i_minus_1, u_i_plus_1, drvrdr, qViscosity(3,3)
  
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

          if (.not.thisOctal%ghostCell(subcell)) then

! first calculate div V

             divV = -9999.d0
             rVec = subcellCentre(thisOctal, subcell)
             r = sqrt(rVec%x**2+rVec%y**2)*gridDistanceScale
	     dir_x = VECTOR(1.d0, 0.d0, 0.d0)
             cen = subcellCentre(thisOctal, subcell)    
             locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir_x, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1, um1, pm1, qViscosity)

             u_i_minus_1 = sqrt(px**2 + py**2) * rhou / rho

             locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
             neighbouroctal => thisoctal
             call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
             call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir_x, q, rho, rhoe, &
                  rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, nd, xnext, px, py, pz, rm1, um1, pm1, qViscosity)

             u_i_plus_1 = sqrt(px**2 + py**2) * rhou / rho


             drvrdr = (u_i_plus_1 - u_i_minus_1) / (2.d0*thisOctal%subcellSize*gridDistanceScale)


             divV = (1.d0/r) * drvrdr +  dudx(thisOctal, subcell, VECTOR(0.d0, 0.d0, 1.d0), VECTOR(0.d0, 0.d0, 1.d0), grid)


! now tau_rr

             thisOctal%qViscosity(subcell,1,1) = thisOctal%etaline(subcell) * thisOctal%rho(subcell) * &
             (2.d0 * dudx(thisOctal, subcell, VECTOR(1.d0, 0.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid) - 0.6666666666d0 * divV)

! now tau_thetatheta
             thisOctal%qViscosity(subcell,2,2) =  thisOctal%qViscosity(subcell,1,1) 

! now tau_zz

             thisOctal%qViscosity(subcell,3,3) = thisOctal%etaline(subcell) * thisOctal%rho(subcell) * &
             (2.d0 * dudx(thisOctal, subcell, VECTOR(0.d0, 0.d0, 1.d0), VECTOR(0.d0, 0.d0, 1.d0), grid) - 0.6666666666d0 * divV)

! now tau_rtheta

             vTheta = thisOctal%rhov(subcell) / (thisOctal%rho(subcell)*r)
             thisOctal%qViscosity(subcell,1,2) = thisOctal%etaline(subcell) * thisOctal%rho(subcell) * &
                  (dudx(thisOctal, subcell, VECTOR(0.d0, 1.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid) - (vTheta/r))

! now tau_thetar

             thisOctal%qViscosity(subcell,2,1) = thisOctal%qViscosity(subcell,1,2) 

! now tau_thetaz

             thisOctal%qViscosity(subcell,2,3) = thisOctal%etaline(subcell) * thisOctal%rho(subcell) * &
             (2.d0 * dudx(thisOctal, subcell, VECTOR(0.d0, 1.d0, 0.d0), VECTOR(0.d0, 0.d0, 1.d0), grid) - 0.6666666666d0 * divV)


! now tau_thetaz
             
             thisOctal%qViscosity(subcell,3,2) = thisOctal%qViscosity(subcell,2,3)


! now tau_rz

             thisOctal%qViscosity(subcell,1,3) = thisOctal%etaline(subcell) * thisOctal%rho(subcell) * &
             (dudx(thisOctal, subcell, VECTOR(0.d0, 0.d0, 1.d0), VECTOR(1.d0, 0.d0, 0.d0), grid) + &
              dudx(thisOctal, subcell, VECTOR(1.d0, 0.d0, 0.d0), VECTOR(0.d0, 0.d0, 1.d0), grid))

! now tau_zr

             thisOctal%qViscosity(subcell,3,1) = thisOctal%qViscosity(subcell,1,3) 





          endif
       endif
    enddo
  end subroutine setupCylindricalViscosity





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

#endif
end module viscosity_mod
