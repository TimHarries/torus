!Photoionization utilities module
!Started 01/04/2011

!This module consolidates established functions and routines that are
!shared by both photoion_mod and photoionAMR_mod

! Structure Guide (quickSearchCode)
! 1. Recombination Rates (pu01)
! 2. Intersect Cube Routines(pm02)
! X. Unused routines (pu0X)

!NOTE: Until this module is well established, routines will have either a P,
!      a U or an F next to their leading comments. These stand for:
!      P - passed: photoionAMR_mod & photoion_mod match for this routine
!      F - failed: Pass criterion unfulfilled (rm or modify this routine!)
!      U - unchecked: yet to verify P or F


module photoion_utils_mod

use constants_mod
use messages_mod
use gridio_mod
use source_mod
use timing
use image_mod
use grid_mod
use amr_mod
use diffusion_mod
use mpi_amr_mod
use mpi_global_mod
use photon_mod
use phasematrix_mod

implicit none

contains

! 1.Recombination Rates (pu01)


! Radiative recombination rate for HI, HeI, HeII and CIV (P)
  real function vernerFerland(a, t, t0, t1, b)
    real(double):: a, t, t0, t1, b
    
    ! based on Verner and Ferland (1996) ApJS 103 467
    vernerFerland = a / (sqrt(t/t0) * (1.d0+sqrt(t/t0))**(1.d0-b) * (1.d0+sqrt(t/t1))**(1.d0+b))
    
  end function vernerFerland


! Radiative recombination rate contribution for CI, CII, CIII,  NI, (P)
! NII, NIII, N IV, OI, OII, OIII,  NeII, NeIII and SI
  function nussbaumerStorey1983(t, a, b, c, d, f) result(rate)
    
    ! based on Nussbaumer & Storey 1983 AA 126 75
    real(double) :: t, a, b, c, d, f, rate
    
    if (t < 0.05d0) then
       rate = tiny(rate)
    else
       rate = (1.d-12) * ((a/t) + b + (c*t) + d*(t**2))*(t**(-1.5d0))*exp(-f/t)
    endif
    
  end function nussbaumerStorey1983


! Radiative recombination rate contribution for CI, CII, CIII, NI, (P)
! NII, NIII, NIV, OI, OII, OIII,  NeII, NeIII
  function ppb1991(z, a, b, c, d, temperature) result(rate)

    ! radiative recombination rates
    ! based on Pequignot, Petitjean & Boisson (1991) A&A 251 680
    
    real(double) :: z, a, b, c, d, temperature, t, rate
    
    t = 1.d-4 * temperature / z**2
    
    rate = 1.d-13 * z * (a*t**b)/(1.d0+c*t**d)
    
  end function ppb1991

! Radiative recombination rate for SI, SII, SIII, SIV (P)
  function svs1982(t, alpharad, xrad) result (rate)

    ! radiative recombination rates based on
    ! Shull and Van Steenberg, 1992, ApJS, 48, 95
    
    real(double) :: t, alpharad, xrad, rate
    
    rate = alpharad * (t /1.d4)**(-xrad)
  end function svs1982

!---awaiting comments                                  (F/P)
!Photoion_mod passes some unnecessary variables so although they 
!don't agree, it should use this one
  subroutine getRecombs(rates, thision, temperature)
    real(double) :: rates(:)
    type(IONTYPE) :: thisIon
    real(double) :: temperature
    
    select case(thisIon%species)
    case("O II")
       rates(1) = 0.d0
       rates(2) = rateFit(7.218d0, -0.575d0, temperature)
       rates(3) = rateFit(4.812d0, -0.575d0, temperature)
       rates(4) = rateFit(3.581d0, -0.495d0, temperature)
       rates(5) = rateFit(1.790d0, -0.495d0, temperature)
    case DEFAULT
       rates = 0.d0
    end select
  end subroutine getRecombs


!-- Awaiting comments
function rateFit(a,b,t) result (rate)
  real(double) :: a, b, t, rate

  rate = 1.d-13 * a * (t/1.d4)**b
end function rateFit


!   End of pu01

!   2.Intersect Cube Routines (pu02)

!return tVal                                                 (P)
  subroutine intersectCube(grid, posVec, i1,i2,i3,direction, tval)
    use vector_mod
    use grid_mod
    implicit none
    type(GRIDTYPE) :: grid
    type(VECTOR) :: direction
    type(VECTOR) :: posVec, norm(6), p3(6)
    real(oct) :: t(6),tval,denom(6)
    integer :: i,j
    logical :: ok, thisOk(6)
    integer :: i1, i2, i3

    ok = .true.
    
    norm(1) = VECTOR(1., 0., 0.)
    norm(2) = VECTOR(0., 1., 0.)
    norm(3) = VECTOR(0., 0., 1.)
    norm(4) = VECTOR(-1., 0., 0.)
    norm(5) = VECTOR(0., -1., 0.)
    norm(6) = VECTOR(0., 0., -1.)
    
    p3(1) = VECTOR(grid%xAxis(i1+1), 0., 0.)
    p3(2) = VECTOR(0.,grid%yAxis(i2+1),0.)
    p3(3) = VECTOR(0.,0.,grid%zAxis(i3+1))
    p3(4) = VECTOR(grid%xAxis(i1), 0., 0.)
    p3(5) = VECTOR(0.,grid%yAxis(i2),0.)
    p3(6) = VECTOR(0.,0.,grid%zAxis(i3))
    
    thisOk = .true.
    
    do i = 1, 6
       
       denom(i) = norm(i) .dot. direction
       if (denom(i) /= 0.) then
          t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
       else
          thisOk(i) = .false.
          t(i) = 0.
       endif
       if (t(i) < 0.) thisOk(i) = .false.
       !      if (denom > 0.) thisOK(i) = .false.
    enddo
    
    j = 0
    do i = 1, 6
       if (thisOk(i)) j=j+1
    enddo
    
    if (j == 0) ok = .false.
    
    if (.not.ok) then
       write(*,*) i1, i2, i3
       write(*,*) direction%x,direction%y,direction%z
       write(*,*) t(1:6)
       stop
    endif
    
    tval = minval(t, mask=thisOk)
    tval = max(tval * 1.001d0,dble((grid%xAxis(2)-grid%xAxis(1))/1000.))
    
    if (tval == 0.) then
       write(*,*) i1, i2, i3,tval
       write(*,*) posVec
       write(*,*) grid%xAxis(i1),grid%yAxis(i2),grid%zAxis(i3)
       write(*,*) grid%xAxis(i1+1),grid%yAxis(i2+1),grid%zAxis(i3+1)
       write(*,*) direction%x,direction%y,direction%z
       write(*,*) t(1:6)
       stop
    endif
    
    if (tval > 2.*(grid%xAxis(2)-grid%xAxis(1))) then
       !     write(*,*) "tval too big",tval,i1,i2,i3,posvec
       !     write(*,*) "direction",direction
       !     write(*,*) t(1:6)
       !     write(*,*) denom(1:6)
    endif

  end subroutine intersectCube
  
!Get tVal for a position on an amr grid                    (P)
  subroutine intersectCubeAMR(grid, posVec, direction, tval)
    implicit none
    type(GRIDTYPE), intent(in)    :: grid
    type(VECTOR), intent(in) :: posVec
    type(VECTOR), intent(in) :: direction
    real(oct), intent(out) :: tval
    !
    type(VECTOR) :: norm(6), p3(6)
    type(OCTAL),pointer :: thisOctal
    type(VECTOR) :: subcen, point
    integer :: subcell
    
    real(oct) :: t(6),denom(6)
    integer :: i,j
    logical :: ok, thisOk(6)
    
    
    point = posVec
    
    call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
    subcen =  subcellCentre(thisOctal,subcell)
    ok = .true.
    
    norm(1) = VECTOR(1.0d0, 0.d0, 0.0d0)
    norm(2) = VECTOR(0.0d0, 1.0d0, 0.0d0)
    norm(3) = VECTOR(0.0d0, 0.0d0, 1.0d0)
    norm(4) = VECTOR(-1.0d0, 0.0d0, 0.0d0)
    norm(5) = VECTOR(0.0d0, -1.0d0, 0.0d0)
    norm(6) = VECTOR(0.0d0, 0.0d0, -1.0d0)
    
    p3(1) = VECTOR(subcen%x+thisOctal%subcellsize/2.0d0, subcen%y, subcen%z)
    p3(2) = VECTOR(subcen%x, subcen%y+thisOctal%subcellsize/2.0d0 ,subcen%z)
    p3(3) = VECTOR(subcen%x,subcen%y,subcen%z+thisOctal%subcellsize/2.0d0)
    p3(4) = VECTOR(subcen%x-thisOctal%subcellsize/2.0d0, subcen%y,  subcen%z)
    p3(5) = VECTOR(subcen%x,subcen%y-thisOctal%subcellsize/2.0d0, subcen%z)
    p3(6) = VECTOR(subcen%x,subcen%y,subcen%z-thisOctal%subcellsize/2.0d0)
    
    thisOk = .true.
    
    do i = 1, 6
       
       denom(i) = norm(i) .dot. direction
       if (denom(i) /= 0.0d0) then
          t(i) = (norm(i) .dot. (p3(i)-posVec))/denom(i)
       else
          thisOk(i) = .false.
          t(i) = 0.0d0
       endif
       if (t(i) < 0.) thisOk(i) = .false.
       !      if (denom > 0.) thisOK(i) = .false.
    enddo
    
    
    j = 0
    do i = 1, 6
       if (thisOk(i)) j=j+1
    enddo
    
    if (j == 0) ok = .false.
    
    if (.not.ok) then
       write(*,*) "Error: j=0 (no intersection???) in lucy_mod::intersectCubeAMR. "
       write(*,*) direction%x,direction%y,direction%z
       write(*,*) t(1:6)
       stop
    endif
    
    tval = minval(t, mask=thisOk)
    tval = max(tval * 1.001d0,dble(thisOctal%subCellSize/1000.))
    
    
    if (tval == 0.) then
       write(*,*) posVec
       write(*,*) direction%x,direction%y,direction%z
       write(*,*) t(1:6)
       stop
    endif
    
    if (tval > sqrt(3.)*thisOctal%subcellsize) then
       !     write(*,*) "tval too big",tval/(sqrt(3.)*thisOctal%subcellSize)
       !     write(*,*) "direction",direction
       !     write(*,*) t(1:6)
       !     write(*,*) denom(1:6)
    endif
    
    
  end subroutine intersectCubeAMR
  

! this is to find a cell intersection for a 2D AMR grid            (P)
! which is essentially a 2D-grid that projects into
! a 3D cylindrical coordinate system
  subroutine intersectCubeAMR2D(grid, posVec, direction, tval)
    
    implicit none
    type(GRIDTYPE), intent(in)    :: grid
    type(VECTOR), intent(in) :: posVec
    type(VECTOR), intent(inout) :: direction
    real(oct), intent(out) :: tval
    !
    type(OCTAL),pointer :: thisOctal
    type(VECTOR) :: subcen, point
    integer :: subcell
    real(double) :: compZ, currentZ
    real(double) :: distToZBoundary, distToXboundary
    real(oct) :: r1,r2,d,cosmu,x1,x2,distTor1,distTor2, theta, mu
    logical :: ok
    type(VECTOR) :: xHat, zHAt
    
    point = posVec
    
    call amrGridValues(grid%octreeRoot, point, foundOctal=thisOctal, foundSubcell=subcell, grid=grid)
    subcen =  subcellCentre(thisOctal,subcell)
    ok = .true.
    
    
    r1 = subcen%x - thisOctal%subcellSize/2.d0
    r2 = subcen%x + thisOctal%subcellSize/2.d0
    d = sqrt(point%x**2+point%y**2)
    xHat = VECTOR(point%x, point%y,0.d0)
    call normalize(xHat)
    
    cosmu =((-1.d0)*xHat).dot.direction
    call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r2**2, x1, x2, ok)
    if (.not.ok) then
       write(*,*) "Quad solver failed in intersectcubeamr2d"
       direction = randomUnitVector()
       x1 = thisoctal%subcellSize/2.d0
       x2 = 0.d0
    endif
    distTor2 = max(x1,x2)
    
    theta = asin(max(-1.d0,min(1.d0,r1 / d)))
    cosmu = xHat.dot.direction
    mu = acos(max(-1.d0,min(1.d0,cosmu)))
    distTor1 = 1.e30
    if (mu  < theta ) then
       call solveQuadDble(1.d0, -2.d0*d*cosmu, d**2-r1**2, x1, x2, ok)
       if (.not.ok) then
          write(*,*) "Quad solver failed in intersectcubeamr2d"
          direction = randomUnitVector()
          x1 = thisoctal%subcellSize/2.d0
          x2 = 0.d0
       endif
       distTor1 = max(x1,x2)
    endif
    
    distToXboundary = min(distTor1, distTor2)
    
    zHat = VECTOR(0.d0, 0.d0, 1.d0)
    compZ = zHat.dot.direction
    currentZ = point%z
    
    if (compZ /= 0.d0 ) then
       if (compZ > 0.d0) then
          distToZboundary = (subcen%z + thisOctal%subcellsize/2.d0 - currentZ ) / compZ
       else
          distToZboundary = abs((subcen%z - thisOctal%subcellsize/2.d0 - currentZ ) / compZ)
       endif
    else
       disttoZboundary = 1.e30
    endif
    
    tVal = min(distToZboundary, distToXboundary) +0.0001d0*grid%halfsmallestsubcell
    if (tVal > 1.e29) then
       write(*,*) tVal,compZ, distToZboundary,disttoxboundary
       write(*,*) "subcen",subcen
       write(*,*) "z",currentZ
    endif
    if (tval < 0.) then
       write(*,*) tVal,compZ, distToZboundary,disttoxboundary
       write(*,*) "subcen",subcen
       write(*,*) "z",currentZ
    endif
    
  end subroutine intersectCubeAMR2D
  
!   End of pu02  
  


!   pu0X Unused routines

! Currently never called. HI recombination rates (P)
  function hRecombination(temperature) result (rate)
    real :: temperature
    real(double) :: rate, t0, t1, b

    t0 = sqrt(temperature/3.148d0)
    t1 = sqrt(temperature/7.036d5)
    b = 0.7480

    rate = 7.982d-11 / (t0*(1.d0+t0)**(1.d0-b) * (1.d0 + t1)**(1.d0+b))
  end function hRecombination

! Currently never called. He recombination rates (P)
  subroutine calcHeRecombs(te, alpha1, alpha21s, alpha21p, alpha23s)
    real :: te, t
    real(double) :: alpha1, alpha21s, alpha21p, alpha23s
    
    t = te / 1.e4
    
    alpha1 = 1.54e-13 * t**(-0.486)

    alpha21s = 2.06e-14 * t**(-0.676)

    alpha21p = 4.17e-14 * t**(-0.861)

    alpha23s = 2.10e-13 * t**(-0.778)
  end subroutine calcHeRecombs


! Currently never called. Recomb rate into ground state. (P)
  function recombToGround(temperature) result (alpha1)
    real :: temperature
    real(double) :: alpha1

    alpha1 = 1.58d-13 * (temperature/1.d4)**(-0.53d0)  ! kenny's photo paper equation 24
  end function recombToGround
  



!   End of pu0X
  
end module photoion_utils_mod
