module viscosity_mod
#ifdef MPI
  
  use messages_mod
  use vector_mod
  use amr_mod
  use mpi_amr_mod
  use mpi_global_mod
  implicit none
  
contains

  logical function floatEqual(a,b)
    real(double) :: a,b

    if (abs(a-b) <1.0d-6 * abs(a)) then
       floatEqual=.true.
    else
       floatEqual = .false.
    endif
    
  end function floatEqual

!  logical function inviscid(thisOctal, subcell) ! if in contact with a ghost cell or an axis cell then cell should be inviscid
!    use inputs_mod, only : smallestCellsize
!    type(octal) :: thisOctal
!    integer :: subcell
!    type(vector) :: cen

!    if (thisOctal%ghostcell(subcell)) then
!       inviscid =.true.
!    else
!       cen = subcellCentre(thisOctal,subcell)
!       if (cen%x <thisOctal%subcellSize) then
!          inviscid = .true.
!       else
!          inviscid = .false.
!       endif
!    endif
!  end function inviscid
  
  real(double) function dudx(thisOctal, subcell, dir_u, dir_x, grid)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: u_i_minus_1, u_i_plus_1, u_i
    type(VECTOR) :: cen, dx, cen_i_plus_1, cen_i_minus_1
    type(VECTOR) :: dir_u, dir_x, locator, tmpV
    real(double) :: q, rho_m1,rho_p1, rhoe, rhou, rhov, rhow, x, qViscosity(3,3)
    real(double) :: xnext, qnext, pressure, flux, phi, phigas, px, py, pz
    real(double) :: rm1, um1, pm1, r, correction, dx_m,dx_p,frac
    integer :: nd, nc

    cen=subcellCentre(thisOctal,subcell)
    r=sqrt(cen%x**2+cen%y**2)*gridDistanceScale
    
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
          tmpV = VECTOR(0.0,0.0,0.0)
       endif
       u_i_minus_1 = (tmpV / rho_m1).dot.dir_u
       if (u_i_minus_1 .ne. u_i_minus_1) then
          write(*,*) neighbourOctal%rho(neighboursubcell)
          write(*,*) cen_i_minus_1 * (1.0/neighbourOctal%subcellsize)
          write(*,*) neighbourOctal%rhou(neighboursubcell),neighbourOctal%rhov(neighboursubcell)/r,&
               neighbourOctal%rhow(neighboursubcell)
       endif
    else
       u_i_minus_1 = 0.d0
       if (rho_m1 /= 0.d0) &
            tmpV = VECTOR(rhou,rhov,rhow)
       tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
            thisOctal%rhov(subcell),&
            thisOctal%rhow(subcell))
       u_i_minus_1 = (tmpV / (rho_m1+thisOctal%rho(subcell))).dot.dir_u
    endif
!    if (inviscid(neighbourOctal, neighbourSubcell)) then
!       dudx =0.0
!       return
!       u_i_minus_1 = u_i
!       cen_i_minus_1=cen
!    endif
    
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
          tmpV = VECTOR(0.0,0.0,0.0)
       endif
       u_i_plus_1 = (tmpV / rho_p1).dot.dir_u
    else
       u_i_plus_1 = 0.d0
       if (rho_p1 /= 0.d0) &
            tmpV = VECTOR(rhou,rhov,rhow)
       tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
            thisOctal%rhov(subcell),&
            thisOctal%rhow(subcell))
       u_i_plus_1 = (tmpV / (rho_p1+thisOctal%rho(subcell))).dot.dir_u
    endif
!        if (inviscid(neighbourOctal, neighbourSubcell)) then
!       dudx=0.0
!       return
!       !       u_i_plus_1 = u_i
!!       cen_i_plus_1=cen
!    endif

    dx=cen_i_plus_1 - cen_i_minus_1
    dx_m = abs((cen - cen_i_minus_1).dot.dir_x)
    dx_p = abs((cen_i_plus_1 - cen ).dot.dir_x)

    if (floatEqual(dx_m,dx_p)) then !plus and minus dxs are equidistant from central cell
       dudx =(u_i_plus_1 - u_i_minus_1) / ((dx.dot.dir_x)*gridDistanceScale/2) ! /2 as dx . dir_x  is twice the distance between 2 cells
    else if (dx_m .eq. 0) then
       dudx =(u_i_plus_1 - u_i) / (dx_p*gridDistanceScale)
    else if (dx_p .eq. 0) then
       dudx =(u_i - u_i_minus_1) / (dx_m*gridDistanceScale)
    elseif (dx_m > dx_p) then  ! if they're not equidistant lerp along the largest cell to find a velocity at the point opposite the centre of the smaller cell and use that
       frac=dx_p/dx_m
       u_i_minus_1=u_i_minus_1*frac + u_i*(1.0-frac)
       dx=(cen_i_plus_1 - cen)
       dudx =(u_i_plus_1 - u_i_minus_1) / ((dx.dot.dir_x)*gridDistanceScale) 
    else
       frac=dx_m/dx_p
       u_i_plus_1=u_i_plus_1*frac + u_i*(1.0-frac)
       dx=(cen - cen_i_minus_1)
       dudx =(u_i_plus_1 - u_i_minus_1 ) / ((dx.dot.dir_x)*gridDistanceScale)
!       if (myrankglobal .eq. 1) then
!          write (*,*) "frac", frac
!          write (*,*) "dx_m/p", dx_m, dx_p
!          write (*,*) dx
!          write (*,*) u_i_minus_1, u_i_plus_1
!          write (*,*) dudx
!          endif
    endif

    if (dudx .ne. dudx .and. myrankglobal .eq.1) then
       write(*,*) "nan in dudx", floatEqual(dx_m,dx_p)
       write(*,*) "dir_x",dir_x
       write(*,*) "dir_u",dir_u
       write(*,*) cen_i_minus_1, cen, cen_i_plus_1
       write(*,*) cen - (3*thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
       write(*,*) locator
       write(*,*) thisOctal%subcellsize, smallestCellSize
       write(*,*) dx_m, dx_p, floatEqual(dx_m,dx_p)
       write(*,*) dx
       write(*,*) "u_is",u_i_minus_1, u_i, u_i_plus_1
       write(*,*) rho_m1, thisOctal%rho(subcell), rho_p1
       write(*,*) dudx
       write(*,*) "-------------"
       stop 1
       dudx = 0.0
    endif    
    
!    dudx = (u_i_plus_1 - u_i_minus_1) / abs((dx.dot.dir_x) * gridDistanceScale)
 !   if (myrankglobal .eq. 1) write(*,*) "dudx", u_i_minus_1, u_i_plus_1, ((dx.dot.dir_x)*gridDistanceScale)
  end function dudx

  real(double) function dudx_onesided(thisOctal, subcell, dir_u, dir_x, grid, backwards)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: u_i_plus_1, u_i_plus_2, u_i
    type(VECTOR) :: cen, dx, cen_i_plus_1, cen_i_plus_2
    type(VECTOR) :: dir_u, dir_x, locator, tmpV
    real(double) :: q, rho_p1,rho_p2, rhoe, rhou, rhov, rhow, x, xnext, qnext, pressure
    real(double) :: flux, phi, phigas, px, py, pz, qViscosity(3,3)
    real(double) :: rm1, um1, pm1, r, correction
    integer :: nd, nc, backDeriv
    logical, optional :: backwards !are we doing a backwards deriv i.e. using values from the lhs of the cell
    real(double) :: coeffs(3)

    backDeriv=1
    if (present(backwards))then !backwards defaults to false
       if (backwards) then
          backDeriv=-1
       endif
    endif
    
    coeffs=(/-3.0/2, 2.0, -1.0/2/)* backDeriv
    
    cen=subcellCentre(thisOctal,subcell)
    r=sqrt(cen%x**2+cen%y**2)*gridDistanceScale
    
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
    
    
    locator = cen + backDeriv*(thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
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
          tmpV = VECTOR(rhou/rho_p1,rhov/r/rho_p1,rhow/rho_p1)
       else
          tmpV = VECTOR(0,0.0,0)
       endif
       u_i_plus_1 = tmpV.dot.dir_u
    else
       u_i_plus_1 = 0.d0
       if (rho_p1 /= 0.d0) &
            tmpV = VECTOR(rhou,rhov,rhow)
       tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
            thisOctal%rhov(subcell),&
            thisOctal%rhow(subcell))
       u_i_plus_1 = (tmpV / (rho_p1+thisOctal%rho(subcell))).dot.dir_u
    endif
!    if (inviscid(neighbourOctal, neighbourSubcell)) then
!       dudx_onesided=0.0
!       return
!       !       u_i_plus_1 = u_i
!!       cen_i_plus_1=cen
!    endif
    
    locator = cen + backDeriv*(3*thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir_x, q, rho_p2, rhoe, & 
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc,xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (abs(rhov) < 1.d-30) rhov = sign(1.d-30,rhov)
    cen_i_plus_2 = VECTOR(px,py,pz)
    if (cylindricalHydro) then
       r = sqrt(px**2 + py**2)*gridDistanceScale
       u_i_plus_2 = 0.d0
       if (rho_p2 /= 0.d0) then
          tmpV = VECTOR(rhou/rho_p2,rhov/r/rho_p2,rhow/rho_p2)
       else
          tmpV = VECTOR(0,0.0,0)
       endif
       u_i_plus_2 = tmpV.dot.dir_u
    else
       u_i_plus_2 = 0.d0
       if (rho_p2 /= 0.d0) tmpV = VECTOR(rhou,rhov,rhow)
       tmpV = tmpV + VECTOR(thisOctal%rhou(subcell),&
            thisOctal%rhov(subcell),&
            thisOctal%rhow(subcell))
       u_i_plus_2 = (tmpV / (rho_p2+thisOctal%rho(subcell))).dot.dir_u
    endif
!    if (inviscid(neighbourOctal, neighbourSubcell)) then
!       dudx_onesided=0.0
!       return
!!       u_i_plus_2 = u_i_plus_1
!!       cen_i_plus_2=cen_i_plus_1
!    endif
    
    dx=cen_i_plus_1 - cen
    dudx_onesided = (u_i*coeffs(1) + u_i_plus_1*coeffs(2) + u_i_plus_2*coeffs(3)) /&
                    (dx.dot.dir_x * gridDistanceScale)

!    write(*,*) "dudx oneside ",dx, u_i, u_i_plus_1, u_i_plus_2,coeffs,dudx_onesided
    if (dudx_onesided .ne. dudx_onesided .and. myrankglobal .eq. 1) then
       write(*,*) "nan in dudx_onesided"
       write(*,*) "dir_x",dir_x
       write(*,*) "dir_u",dir_u
       write(*,*) cen
       write(*,*) dx
       write(*,*) u_i,u_i_plus_1, u_i_plus_2
       write(*,*) dudx_onesided
       write(*,*) "-------------"
    endif
    
  end function dudx_onesided
  
  real(double) function ddudxx(thisOctal, subcell, dir_u, dir_x, grid)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: u_i_minus_1, u_i_plus_1, u_i
    type(VECTOR) :: cen, dx, cen_i_plus_1, cen_i_minus_1
    type(VECTOR) :: dir_u, dir_x, locator, tmpV
    real(double) :: q, rho_m1,rho_p1, rhoe, rhou, rhov, rhow, x, xnext, qnext, pressure
    real(double) :: flux, phi, phigas, px, py, pz, qViscosity(3,3)
    real(double) :: rm1, um1, pm1, r, correction
    integer :: nd, nc
    real(double) :: dx_m, dx_p, frac
    
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
       else
          u_i_minus_1 = 0.0
       endif
    else
       u_i_minus_1 = 0.d0
       if (rho_m1 /= 0.d0) then
          tmpV = VECTOR(rhou,rhov,rhow)
          u_i_minus_1 = (tmpV / rho_m1).dot.dir_u
       endif
    endif
!    if (inviscid(neighbourOctal, neighbourSubcell)) then
!       ddudxx=0.0
!       return
!       !       u_i_minus_1 = u_i
!!       cen_i_minus_1=cen
!    endif
    
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
       else
          u_i_plus_1=0.0
       endif
    else
       u_i_plus_1 = 0.d0
       if (rho_p1 /= 0.d0) then
          tmpV = VECTOR(rhou,rhov,rhow)
          u_i_plus_1 = (tmpV / rho_p1).dot.dir_u
       endif
    endif
!    if (inviscid(neighbourOctal, neighbourSubcell)) then
!       ddudxx=0.0
!       return
!       !       u_i_plus_1 = u_i
!!       cen_i_plus_1=cen
!    endif
    
    dx = (cen_i_plus_1 - cen_i_minus_1)
    dx_m = abs((cen - cen_i_minus_1).dot.dir_x)
    dx_p = abs((cen_i_plus_1 - cen ).dot.dir_x)
    
    if (floatEqual(dx_m,dx_p)) then !plus and minus dxs are equidistant from central cell
       ddudxx =(u_i_plus_1 + u_i_minus_1 - 2.0*u_i) / ((dx.dot.dir_x)*gridDistanceScale/2)**2 ! /2 as dx . dir_x  is twice the distance between 2 cells
    elseif (dx_m > dx_p) then  ! if they're not equidistant lerp along the largest cell to find a velocity at the point opposite the centre of the smaller cell and use that
       frac=dx_p/dx_m
       u_i_minus_1=u_i_minus_1*frac + u_i*(1.0-frac)
       dx=(cen_i_plus_1 - cen)
       ddudxx =(u_i_plus_1 + u_i_minus_1 - 2.0*u_i) / ((dx.dot.dir_x)*gridDistanceScale)**2 
    else
       frac=dx_m/dx_p
       u_i_plus_1=u_i_plus_1*frac + u_i*(1.0-frac)
       dx=(cen - cen_i_minus_1)
       ddudxx =(u_i_plus_1 + u_i_minus_1 - 2.0*u_i) / ((dx.dot.dir_x)*gridDistanceScale)**2
    endif
!    if (myrankglobal .eq. 1) write(*,*) "ddudxx", u_i_minus_1, u_i_plus_1, -2*u_i, ((dx.dot.dir_x)*gridDistanceScale)**2

!    if (ddudxx .ne. ddudxx .and. myrankglobal .eq. 1) then
!       write(*,*) "nan in ddudxx"
!       write(*,*) cen
!       write(*,*) dx
!       write(*,*) u_i_minus_1, u_i,u_i_plus_1
!       write(*,*) ddudxx
!       write(*,*) "-------------"
!    endif

  end function ddudxx

    real(double) function ddudxx_onesided(thisOctal, subcell, dir_u, dir_x, grid, backwards)
    use inputs_mod, only : smallestCellSize, gridDistanceScale, cylindricalHydro
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: u_i, u_i_plus_1, u_i_plus_2, u_i_plus_3
    type(VECTOR) :: cen, dx, cen_i_plus_1, cen_i_plus_2, cen_i_plus_3
    type(VECTOR) :: dir_u, dir_x, locator, tmpV
    real(double) :: q, rho, rhoe, rhou, rhov, rhow, x, xnext, qnext, pressure, flux, phi, phigas, px, py, pz, qViscosity(3,3)
    real(double) :: rm1, um1, pm1, r, correction
    integer :: nd, nc, backDeriv
    logical, optional :: backwards !are we doing a backwards deriv i.e. using values from the lhs of the cell
    real(double) :: coeffs(4)

    backDeriv=1
    if (present(backwards))then !backwards defaults to false
       if (backwards) then
          backDeriv=-1
       endif
    endif
    
    coeffs=(/2.0, -5.0, 4.0, -1.0 /)*backDeriv
    
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

    locator = cen + backDeriv*(thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir_x, q, rho, rhoe, & 
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc,xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (abs(rhov) < 1.d-30) rhov = sign(1.d-30,rhov)
    cen_i_plus_1 = VECTOR(px,py,pz)
    if (cylindricalHydro) then
       u_i_plus_1 = 0.d0
       if (rho /= 0.d0) then
          r = sqrt(px**2 + py**2)*gridDistanceScale
          tmpV = VECTOR(rhou,rhov/r,rhow)
          u_i_plus_1 = (tmpV / rho).dot.dir_u
       else
          u_i_plus_1 =0.0
       endif
    else
       u_i_plus_1 = 0.d0
       if (rho /= 0.d0) then
          tmpV = VECTOR(rhou,rhov,rhow)
          u_i_plus_1 = (tmpV / rho).dot.dir_u
       endif
    endif
!    if (inviscid(neighbourOctal, neighbourSubcell)) then
!       ddudxx_onesided=0.0
!       return
!       !       u_i_plus_1 = u_i
!!       cen_i_plus_1=cen
!    endif
    
    
    locator = cen + backDeriv*(3*thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir_x, q, rho, rhoe, & 
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc,xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (abs(rhov) < 1.d-30) rhov = sign(1.d-30,rhov)
    cen_i_plus_2 = VECTOR(px,py,pz)
    if (cylindricalHydro) then
       u_i_plus_2 = 0.d0
       if (rho /= 0.d0) then
          r = sqrt(px**2 + py**2)*gridDistanceScale
          tmpV = VECTOR(rhou,rhov/r,rhow)
          u_i_plus_2 = (tmpV / rho).dot.dir_u
       else
          u_i_plus_2 =0.0
       endif
    else
       u_i_plus_1 = 0.d0
       if (rho /= 0.d0) then
          tmpV = VECTOR(rhou,rhov,rhow)
          u_i_plus_2 = (tmpV / rho).dot.dir_u
       endif
    endif
!    if (inviscid(neighbourOctal, neighbourSubcell)) then
!       ddudxx_onesided=0.0
!       !       u_i_plus_2 = u_i_plus_1
!!    cen_i_plus_2=cen_i_plus_1
!    endif
    
    locator = cen + backDeriv*(5*thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir_x
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir_x, q, rho, rhoe, & 
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc,xnext, px, py, pz, rm1, um1, pm1, qViscosity)
    if (abs(rhov) < 1.d-30) rhov = sign(1.d-30,rhov)
    cen_i_plus_3 = VECTOR(px,py,pz)
    if (cylindricalHydro) then
       u_i_plus_3 = 0.d0
       if (rho /= 0.d0) then
          r = sqrt(px**2 + py**2)*gridDistanceScale
          tmpV = VECTOR(rhou,rhov/r,rhow)
          u_i_plus_3 = (tmpV / rho).dot.dir_u
       endif
    else
       u_i_plus_1 = 0.d0
       if (rho /= 0.d0) then
          tmpV = VECTOR(rhou,rhov,rhow)
          u_i_plus_3 = (tmpV / rho).dot.dir_u
       endif
    endif
!    if (inviscid(neighbourOctal, neighbourSubcell)) then
!       ddudxx_onesided=0.0
!       return
!       !       u_i_plus_3 = u_i_plus_2
!!       cen_i_plus_3=cen_i_plus_2
!    endif

    dx=(cen - cen_i_plus_3)
    ddudxx_onesided = (u_i*coeffs(1) + u_i_plus_1*coeffs(2) + u_i_plus_2*coeffs(3) + u_i_plus_3*coeffs(4)) /&
                      ((dx.dot.dir_x) * gridDistanceScale/3.0)**2

!    if (myrankglobal .eq. 1) write(*,*) "ddudxx oneside ",dx, u_i, u_i_plus_1, u_i_plus_2
    
  end function ddudxx_oneSided
  
  real(double) function ddudxdz(thisOctal, subcell, dir_u, grid, oneSideIn)
    use inputs_mod, only : smallestCellSize, gridDistanceScale
    type(GRIDTYPE) :: grid
    type(OCTAL), pointer :: thisOctal, neighbourOctal
    integer :: subcell, neighbourSubcell
    real(double) :: dudx_m1, dudx_p1, dudx_m2, dudx_0
    type(VECTOR) :: cen,  cen_i_plus_1, cen_i_minus_1, locator, dir_u
    logical      :: backx
    integer :: oneSide(2)
    integer, optional :: oneSideIn(2)

    if (.not. present(oneSideIn)) then
       oneSide=(/ 0,0 /)
    else
       oneSide=oneSideIn
    endif
    ddudxdz=1.0d-100
    cen = subcellCentre(thisOctal,subcell)
 !   if ((cen%x - (3.0/2*thisOctal%subcellSize+0.01d0*smallestCellSize)) > 0.d0) then ! not on axis
    if (oneSide(2).eq.0) then
       if (oneSide(1).eq.0) then
          locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          dudx_m1 = dudx(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0, 0.0), grid)
          cen_i_minus_1= subcellCentre(neighbourOctal, neighbourSubcell)
          locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          dudx_p1 = dudx(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0,0.0), grid)
          cen_i_plus_1= subcellCentre(neighbourOctal, neighbourSubcell)
       else
          backx= oneSide(1).eq.-1
          locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          dudx_m1 = dudx_oneSided(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0, 0.0), grid, backx)
          cen_i_minus_1= subcellCentre(neighbourOctal, neighbourSubcell)
          locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          dudx_p1 = dudx_oneSided(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0,0.0), grid,backx)
          cen_i_plus_1= subcellCentre(neighbourOctal, neighbourSubcell)
       endif
       ddudxdz = (dudx_p1 - dudx_m1) /&
             abs((cen_i_plus_1-cen_i_minus_1).dot.vector(0.0, 0.0, 1.0)*gridDistanceScale)
    else
       if (oneSide(1).eq.0) then
          dudx_0=dudx_onesided(thisOctal,Subcell, dir_u, vector(1.0, 0.0, 0.0), grid)*(-3.0/2.0)*oneSide(1)
          locator = cen + oneside(2)*1.0*(thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          dudx_m1 = dudx(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0, 0.0), grid)*2.0*oneSide(1)
          locator = cen + oneside(2)*1.0*(thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          dudx_m2 = dudx(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0, 0.0), grid)*(-1.0/2.0)*oneSide(1)
          cen_i_minus_1= subcellCentre(neighbourOctal, neighbourSubcell)
             
       else
             
          dudx_0=dudx_onesided(thisOctal,Subcell, dir_u, vector(1.0, 0.0, 0.0), grid,&
                               backwards=backx)*(-3.0/2.0)*oneSide(1)
             
          locator = cen + oneside(2)*1.0*(thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          dudx_m1 = dudx_onesided(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0, 0.0), grid,&
                                  backwards=backx)*2.0*oneSide(1)
             
             
          locator = cen + oneside(2)*1.0*(thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*vector(0.0, 0.0, 1.0)
          neighbouroctal => thisoctal
          call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
          dudx_m2 = dudx_onesided(neighbourOctal,neighbourSubcell, dir_u, vector(1.0, 0.0, 0.0), grid,&
                                  backwards=backx)*(-1.0/2.0)*oneSide(1)
          cen_i_minus_1= subcellCentre(neighbourOctal, neighbourSubcell)
             
       endif
!          if (myrankglobal .eq. 1) write(*,*) "dudx's ", dudx_0,dudx_m1,dudx_m2
          ddudxdz = (dudx_0 + dudx_m1 + dudx_m2) /&
               abs((cen-cen_i_minus_1).dot.vector(0.0, 0.0, 1.0)*(gridDistanceScale*0.5))
    endif
!    else
!       if (myrankglobal .eq. 1) write(*,*) "on axis ", cen
!       ddudxdz=0
!    endif

    if (ddudxdz .ne. ddudxdz) then
       write(*,*) "nan in ddudxdz"
       write(*,*) cen
       write(*,*) (cen-cen_i_minus_1)
       write(*,*) dudx_m1,dudx_0, dudx_p1
       write(*,*) ddudxdz
       write(*,*) "-------------"
    endif
    
  end function ddudxdz
  
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
    type(octal), pointer   :: thisoctal
    type(VECTOR) :: out
    integer :: subcell 
    type(VECTOR) :: dir, cen
    real(double) :: r, invR
    real(double) :: drdrr,drdrz,drdzz, dzdrr,dzdrz,dzdzz, dqdrr,dqdzz
    real(double) :: drdr,drdz,dqdr,dzdr

    out = VECTOR(0.d0, 0.d0, 0.d0)
    
    dir = VECTOR(1.d0, 0.d0, 0.d0)
    if (thisOctal%ghostCell(subcell)) return
    cen = subcellCentre(thisOctal,subcell)
    
!    if (amrGridSize/2+amrGridCentreX < cen%x .or.&
!        amrGridSize/2+amrGridCentreZ < cen%z .or.&
!        amrGridCentreX-amrGridSize/2 > cen%x .or.&
!        amrGridCentreZ-amrGridSize/2 > cen%z) return
    
    r = cen%x * gridDistanceScale
    invR=1.0/r

    if (cen%x > 0.d0) then
       drdr =dudx(thisOctal, subcell,   vector(1.0 ,0.0 ,0.0), vector(1.0, 0.0, 0.0),grid)
       dqdr =dudx(thisOctal, subcell,   vector(0.0 ,1.0, 0.0), vector(1.0, 0.0, 0.0),grid)
       dzdr =dudx(thisOctal, subcell,   vector(0.0 ,0.0 ,1.0), vector(1.0, 0.0, 0.0),grid)
       drdz =dudx(thisOctal, subcell,   vector(1.0, 0.0, 0.0), vector(0.0, 0.0, 1.0),grid)
       drdrr=ddudxx(thisOctal, subcell, vector(1.0, 0.0 ,0.0), vector(1.0, 0.0, 0.0),grid)
       dqdrr=ddudxx(thisOctal, subcell, vector(0.0, 1.0, 0.0), vector(1.0, 0.0, 0.0),grid)
       dqdzz=ddudxx(thisOctal, subcell, vector(0.0, 1.0, 0.0), vector(0.0, 0.0, 1.0),grid)
       dzdrr=ddudxx(thisOctal, subcell, vector(0.0 ,0.0 ,1.0), vector(1.0, 0.0, 0.0),grid)
       drdzz=ddudxx(thisOctal, subcell, vector(1.0, 0.0, 0.0), vector(0.0, 0.0, 1.0),grid)
       dzdzz=ddudxx(thisOctal, subcell, vector(0.0, 0.0, 1.0), vector(0.0, 0.0, 1.0),grid)
       
       if (all(thisOctal%isOnBoundary(subcell,:).eq.0)) then
          drdrz=ddudxdz(thisOctal, subcell, vector(1.0,0,0),grid, (/0,0/))
          dzdrz=ddudxdz(thisOctal, subcell, vector(0,0,1.0),grid, (/0,0/))
!          if (myrankglobal .eq. 1) write(*,*) "case0",drdrz,dzdrz
       else if ((thisOctal%isOnBoundary(subcell,1).ne.0 .and.&
            thisOctal%isOnBoundary(subcell,2).ne.0).or.&
         (thisOctal%isOnBoundary(subcell,3).ne.0 .and.&
         thisOctal%isOnBoundary(subcell,4).ne.0)) then
          write(*,*) "trying to calc derivs with boundaries on both sides(ghost cell?) at: "
          write(*,*) cen
          write(*,*) "boundaries: "
          write(*,*) thisOctal%isOnBoundary(subcell,:)
          drdrz=0.0
          dzdrz=0.0
       else
          drdrz=0.0
          dzdrz=0.0
          !isOnBondary is an integer array set to 1 if the cell is on an mpi in the directions specified
          !1: x-right,   2: x-left,   3: z-down,   4: z-up
          !if (myrankglobal .eq. 1) write(*,*) "on a boundary ", thisOctal%isOnBoundary(subcell,:)
          if (thisOctal%isOnBoundary(subcell,3).eq.0 .and.&
               thisOctal%isOnBoundary(subcell,4).eq.0) then         
             drdrz=ddudxdz(thisOctal, subcell, vector(1.0,0,0), grid)
             dzdrz=ddudxdz(thisOctal, subcell, vector(0,0,1.0), grid)
            ! if (myrankglobal .eq. 1) write(*,*) "case1",drdrz,dzdrz
          else if (thisOctal%isOnBoundary(subcell,3).ne.0) then
             drdrz=ddudxdz(thisOctal, subcell, vector(1.0,0,0),grid, (/0,-1/))
             dzdrz=ddudxdz(thisOctal, subcell, vector(0,0,1.0),grid, (/0,-1/))
    !         drdzz=ddudxx_onesided(thisOctal, subcell, vector(1.0,0,0),vector(0,0,1.0),grid, .true.)
    !         dzdzz=ddudxx_onesided(thisOctal, subcell, vector(0,0,1.0),vector(0,0,1.0),grid, .true.)
            ! if (myrankglobal .eq. 1) write(*,*) "case2",drdrz,dzdrz
          else 
             drdrz=ddudxdz(thisOctal, subcell, vector(1.0,0,0),grid, (/0,1/))
             dzdrz=ddudxdz(thisOctal, subcell, vector(0,0,1.0),grid, (/0,1/))
     !        drdzz=ddudxx_onesided(thisOctal, subcell, vector(1.0,0,0),vector(0,0,1.0),grid, .false.)
     !        dzdzz=ddudxx_onesided(thisOctal, subcell, vector(0,0,1.0),vector(0,0,1.0),grid, .false.)
          endif
       endif
       
       out%x=thisOctal%etaline(subcell)*(4.0/3*drdrr&
                                        +drdzz&
                                        +1.0/3*dzdrz&
                                        +4.0*invR/3*drdr&
                                        -4.0*invR**2/3*thisOctal%rhou(subcell)/thisOctal%rho(subcell))
       if (abs(r)< smallestCellSize) then
          out%y=0.0
       else
          out%y=thisOctal%etaline(subcell)*(dqdrr&
                                           +dqdzz&
                                           +invR*dqdr&
                                           -thisOctal%rhov(subcell)/thisOctal%rho(subcell)/r**3)
       endif
       out%z=thisOctal%etaline(subcell)*(1.0/3*drdrz&
                                        +dzdrr&
                                        +4.0/3*dzdzz&
                                        +invR/3*drdz&
                                        +invR*dzdr)
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
    type(VECTOR) :: dir, cen, locator
    real(double) :: rm1, um1, pm1, r,r1 ,r2
    integer :: nd, i, j, nc
    
    doMultByR = .false.
    if (present(multbyr)) doMultByR = multByR
    
    dqdx = 0.d0
    
    cen = subcellCentre(thisOctal,subcell)
    r = sqrt(cen%x**2 + cen%y**2) * gridDistanceScale
    
    
    locator = cen + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, &
         px, py, pz, rm1, um1, pm1, qViscosityP1)
    cen_i_plus_1 = VECTOR(px, py,pz)
    r2 = abs(px)*gridDistanceScale
    
    locator = cen_i_plus_1 + (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, dir, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, &
         px, py, pz, rm1, um1, pm1, qViscosityP2)
    
    locator = cen - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, &
         px, py, pz, rm1, um1, pm1, qViscosityM1)
    cen_i_minus_1 = VECTOR(px, py,pz)
    r1 = abs(px)*gridDistanceScale
    
    locator = cen_i_minus_1 - (thisOctal%subcellSize/2.d0 + 0.1d0*smallestCellSize)*dir
    neighbouroctal => thisoctal
    call findsubcelllocal(locator, neighbouroctal, neighboursubcell)
    call getneighbourvalues(grid, thisoctal, subcell, neighbouroctal, neighboursubcell, (-1.d0)*dir, q, rho, rhoe, &
         rhou, rhov, rhow, x, qnext, pressure, flux, phi, phigas, correction, nd, nc, xnext, &
         px, py, pz, rm1, um1, pm1, qViscosityM2)
    
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
                  (2.d0 * dudx(thisOctal, subcell, VECTOR(1.d0, 0.d0, 0.d0), VECTOR(1.d0, 0.d0, 0.d0), grid)&
                  - 0.6666666666d0 * divV)
             
             ! now tau_thetatheta
             thisOctal%qViscosity(subcell,2,2) =  thisOctal%etaline(subcell) * &
                  (2.d0 * thisOctal%rhou(subcell)/(thisOctal%rho(subcell)*r)&
                  - 0.6666666666d0 * divV)
             
             ! now tau_zz
             
             thisOctal%qViscosity(subcell,3,3) = thisOctal%etaline(subcell) * &
                  (2.d0 * dudx(thisOctal, subcell, VECTOR(0.d0, 0.d0, 1.d0), VECTOR(0.d0, 0.d0, 1.d0), grid)&
                  - 0.6666666666d0 * divV)
             
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
             
             if (.false.) then !maxval(abs(thisOctal%qViscosity(subcell,:,:)))> 1) then
                write(*,*) "Q ",thisOctal%qViscosity(subcell,1,2), " mu r domegabydr ", fac, " ratio ",&
                     thisOctal%qViscosity(subcell,1,2)/fac
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
    use inputs_mod, only : gridDistanceScale
    type(gridtype) :: grid
    type(octal), pointer   :: thisoctal
    type(octal), pointer  :: child 
    real(double) :: dt, thisTime, acc, r, rho
    type(VECTOR) :: fVisc, vel, cen
    integer :: subcell, i
    
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
             cen = subcellCentre(thisOctal,subcell)
             rho=thisOctal%rho(subcell)
             r=sqrt(cen%x**2+cen%y**2)*GridDistanceScale
             fVisc =  newdivQ(thisOctal, subcell,  grid)
             acc = sqrt(fVisc%x**2  + fvisc%z**2)/rho
             acc = max(acc, 1.d-60)
             vel=vector(thisOctal%rhou(subcell),&
                        thisOctal%rhov(subcell)/r,&
                        thisOctal%rhow(subcell))/rho
!             write(*,*) "vel ", vel
             if (acc .ne. 0) then
                thisTime = sqrt(thisOctal%subcellSize*gridDistanceScale/acc)*0.5!*cflnumber
             endif
             thisTime = min(thisTime, &
                            abs(max(abs(vel%y),1.0d4)/&
                                   (fvisc%y/rho/r)/10.0)) !limit change in angular speed to 0.1km/s or rotation speed per step
             thisTime = min(thisTime, &
                            max(1.0d5, sqrt(vel%x**2+vel%z**2)/&
                                       sqrt(fvisc%x**2+fvisc%z**2)/rho)) !limit change in velocity max(1km/s,currentVel)  per step
             dt = min(thisTime, dt)
          endif
       endif
    enddo
  end subroutine viscousTimescaleCylindrical
  
#endif
end module viscosity_mod



   !if ((cen%x - (thisOctal%subcellSize/2.d0+0.1d0*smallestCellSize)) > 0.d0) then ! not on axis
    !       out%x = out%x + dqdx(thisOctal, subcell, grid, 1, 1, VECTOR(1.d0, 0.d0, 0.d0))
    !       out%x = out%x + thisOctal%qViscosity(subcell,1,1) * invR
    !       out%x = out%x - thisOctal%qViscosity(subcell,2,2) * invR
    !       out%x = out%x + dqdx(thisOctal, subcell, grid, 1, 3, VECTOR(0.d0, 0.d0, 1.d0))
    !       
    !       out%y = out%y  + dqdx(thisOctal, subcell, grid, 1, 2, VECTOR(1.d0, 0.d0, 0.d0))
    !       out%y = out%y + 2.d0 * thisOctal%qViscosity(subcell,1,2) * invR
    !       out%y = out%y  + dqdx(thisOctal, subcell, grid, 2, 3, VECTOR(0.d0, 0.d0, 1.d0))
    !       
    !       out%z = out%z + dqdx(thisOctal, subcell, grid, 1, 3, VECTOR(1.d0, 0.d0, 0.d0))
    !       out%z = out%z + thisOctal%qViscosity(subcell,1,3) * invR
    !       out%z = out%z + dqdx(thisOctal, subcell, grid, 3, 3, VECTOR(0.d0, 0.d0, 1.d0))
    !    endif
!old divQ method
