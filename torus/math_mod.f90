! 
! this module holds various maths subroutines and NR subroutines
! such as searching, sorting and interpolation
!

! written by tjh

! v1.0 on 16/09/99

module math_mod


  use grid_mod            ! opacity grid
  use vector_mod          ! vector maths
  use constants_mod       ! physical constants
  use utils_mod
  use kind_mod

  implicit none

  public



contains


  ! interpolated in 3d for grid in chiline

  real function interpGridChil(grid,  i1, i2, i3, t1, t2, t3)
    implicit none
    type(GRIDTYPE) :: grid
    integer i1,i2,i3
    real :: t1, t2, t3

!    weight = 0.
!    if (grid%inUse(i1  , i2  , i3  )) weight(1) =  ((1.-t1)  * (1.-t2) * (1.-t3))
!    if (grid%inUse(i1+1, i2  , i3  )) weight(2) =  ((t1   )  * (1.-t2) * (1.-t3))
!    if (grid%inUse(i1+1, i2+1, i3  )) weight(3) =  ((t1   )  * (t2   ) * (1.-t3))
!    if (grid%inUse(i1  , i2+1, i3+1)) weight(4) =  ((1.-t1)  * (t2   ) * (t3   ))
!    if (grid%inUse(i1  , i2+1, i3  )) weight(5) =  ((1.-t1)  * (t2   ) * (1.-t3))
!    if (grid%inUse(i1+1, i2  , i3+1)) weight(6) =  ((t1   )  * (1.-t2) * (t3   ))
!    if (grid%inUse(i1  , i2  , i3+1)) weight(7) =  ((1.-t1)  * (1.-t2) * (t3   ))
!    if (grid%inUse(i1+1, i2+1, i3+1)) weight(8) =  ((t1   )  * (t2   ) * (t3   ))


    interpGridChil = &
         ((1.d0-t1)  * (1.d0-t2) * (1.d0-t3))* grid%chiline(i1  , i2   , i3   ) + &
         ((t1   )  * (1.d0-t2) * (1.d0-t3))* grid%chiline(i1+1, i2   , i3   ) + &
         ((t1   )  * (t2   ) * (1.d0-t3))* grid%chiline(i1+1, i2+1 , i3   ) + &
         ((1.d0-t1)  * (t2   ) * (t3   ))* grid%chiline(i1  , i2+1 , i3+1 ) + &
         ((1.d0-t1)  * (t2   ) * (1.d0-t3))* grid%chiline(i1  , i2+1 , i3   ) + &
         ((t1   )  * (1.d0-t2) * (t3   ))* grid%chiline(i1+1, i2   , i3+1 ) + &
         ((1.d0-t1)  * (1.d0-t2) * (t3   ))* grid%chiline(i1  , i2   , i3+1 ) + &
         ((t1   )  * (t2   ) * (t3   ))* grid%chiline(i1+1, i2+1 , i3+1 )

!    interpGridChil = &
!          grid%chiline(i1  , i2   , i3   ) * weight(1)+ &
!          grid%chiline(i1+1, i2   , i3   ) * weight(2)+ &
!          grid%chiline(i1+1, i2+1 , i3   ) * weight(3)+ &
!          grid%chiline(i1  , i2+1 , i3+1 ) * weight(4)+ &
!          grid%chiline(i1  , i2+1 , i3   ) * weight(5)+ &
!          grid%chiline(i1+1, i2   , i3+1 ) * weight(6)+ &
!          grid%chiline(i1  , i2   , i3+1 ) * weight(7)+ &
!          grid%chiline(i1+1, i2+1 , i3+1 ) * weight(8)

!    totWeight = SUM(weight)
!    if (totWeight /= 0.) then
!       interpGridChil = interpGridChil / totWeight
!    else
!       interpGridChil = 1.e10
!    endif

 
 end  function interpGridChil

  real function interpGridKappaSca(grid,  i1, i2, i3, iLambda, t1, t2, t3)
    implicit none
    type(GRIDTYPE) :: grid
    integer i1,i2,i3, iLambda
    real :: t1, t2, t3



    interpGridKappaSca = &
         ((1.e0-t1)  * (1.e0-t2) * (1.e0-t3))* log10(grid%kappaSca(i1  , i2   , i3   , ilambda)) + &
         ((t1   )  *   (1.e0-t2) * (1.e0-t3))* log10(grid%kappaSca(i1+1, i2   , i3   , ilambda)) + &
         ((t1   )  *   (t2   ) *   (1.e0-t3))* log10(grid%kappaSca(i1+1, i2+1 , i3   , ilambda)) + &
         ((1.e0-t1)  * (t2   ) *     (t3   ))* log10(grid%kappaSca(i1  , i2+1 , i3+1 , ilambda)) + &
         ((1.e0-t1)  * (t2   )   * (1.e0-t3))* log10(grid%kappaSca(i1  , i2+1 , i3   , ilambda)) + &
         ((t1   )  *   (1.e0-t2) *   (t3   ))* log10(grid%kappaSca(i1+1, i2   , i3+1 , ilambda)) + &
         ((1.e0-t1)  * (1.e0-t2) *   (t3   ))* log10(grid%kappaSca(i1  , i2   , i3+1 , ilambda)) + &
         ((t1   )  *   (t2   ) *     (t3   ))* log10(grid%kappaSca(i1+1, i2+1 , i3+1 , ilambda))

    
    interpGridKappaSca = 10.e0**interpGridKappaSca

  end  function interpGridKappaSca

  real function interpGridKappaScaRed(grid,  i1, i2, i3, iLambda, t1, t2, t3)
    implicit none
    type(GRIDTYPE) :: grid
    integer i1,i2,i3, iLambda
    real :: t1, t2, t3

    interpGridKappaScaRed = &
         ((1.d0-t1)  * (1.d0-t2) * (1.d0-t3))* log10(grid%kappaScaRed(i1  , i2   , i3   , ilambda)) + &
         ((t1   )  * (1.d0-t2) * (1.d0-t3))* log10(grid%kappaScaRed(i1+1, i2   , i3   , ilambda)) + &
         ((t1   )  * (t2   ) * (1.d0-t3))* log10(grid%kappaScaRed(i1+1, i2+1 , i3   , ilambda)) + &
         ((1.d0-t1)  * (t2   ) * (t3   ))* log10(grid%kappaScaRed(i1  , i2+1 , i3+1 , ilambda)) + &
         ((1.d0-t1)  * (t2   ) * (1.d0-t3))* log10(grid%kappaScaRed(i1  , i2+1 , i3   , ilambda)) + &
         ((t1   )  * (1.d0-t2) * (t3   ))* log10(grid%kappaScaRed(i1+1, i2   , i3+1 , ilambda)) + &
         ((1.d0-t1)  * (1.d0-t2) * (t3   ))* log10(grid%kappaScaRed(i1  , i2   , i3+1 , ilambda)) + &
         ((t1   )  * (t2   ) * (t3   ))* log10(grid%kappaScaRed(i1+1, i2+1 , i3+1 , ilambda))

    interpGridKappaScaRed = 10.e0**interpGridKappaScaRed

  end  function interpGridKappaScaRed


  real function interpGridKappaAbs(grid,i1, i2, i3, iLambda,t1, t2, t3)
    implicit none
    type(GRIDTYPE) :: grid
    integer i1,i2,i3, iLambda
    real :: t1, t2, t3



    interpGridKappaAbs = &
         ((1.e0-t1)  * (1.e0-t2) * (1.e0-t3)) * log10((grid%kappaAbs(i1  , i2   , i3   , ilambda))) + &
         ((t1   )  * (1.e0-t2) * (1.e0-t3))   * log10((grid%kappaAbs(i1+1, i2   , i3   , ilambda))) + &
         ((t1   )  * (t2   ) * (1.e0-t3))     * log10((grid%kappaAbs(i1+1, i2+1 , i3   , ilambda))) + &
         ((1.e0-t1)  * (t2   ) * (t3   ))     * log10((grid%kappaAbs(i1  , i2+1 , i3+1 , ilambda))) + &
         ((1.d0-t1)  * (t2   ) * (1.e0-t3))   * log10((grid%kappaAbs(i1  , i2+1 , i3   , ilambda))) + &
         ((t1   )  * (1.e0-t2) * (t3   ))     * log10((grid%kappaAbs(i1+1, i2   , i3+1 , ilambda))) + &
         ((1.e0-t1)  * (1.e0-t2) * (t3   ))   * log10((grid%kappaAbs(i1  , i2   , i3+1 , ilambda))) + &
         ((t1   )  * (t2   ) * (t3   ))       * log10((grid%kappaAbs(i1+1, i2+1 , i3+1 , ilambda)))

    if (interpGridKappaAbs  > 20.) then
       write(*,*) grid%kappaAbs(i1-1:i1+1,i2-1:i2+1,i3-1:i3+1,1)
       write(*,*) i1,i2,i3,t1,t2,t3
       write(*,*) interpGridKappaAbs
    endif

    interpGridKappaAbs = 10.e0**interpGridKappaAbs


  end function interpGridKappaAbs

  real function interpGridKappaAbsRed(grid,i1, i2, i3, iLambda,t1, t2, t3)
    implicit none
    type(GRIDTYPE) :: grid
    integer i1,i2,i3, iLambda
    real :: t1, t2, t3



    interpGridKappaAbsRed = &
         ((1.d0-t1)  * (1.d0-t2) * (1.d0-t3))* log10((grid%kappaAbsRed(i1  , i2   , i3   , ilambda))) + &
         ((t1   )  * (1.d0-t2) * (1.d0-t3))* log10((grid%kappaAbsRed(i1+1, i2   , i3   , ilambda))) + &
         ((t1   )  * (t2   ) * (1.d0-t3))* log10((grid%kappaAbsRed(i1+1, i2+1 , i3   , ilambda))) + &
         ((1.d0-t1)  * (t2   ) * (t3   ))* log10((grid%kappaAbsRed(i1  , i2+1 , i3+1 , ilambda))) + &
         ((1.d0-t1)  * (t2   ) * (1.d0-t3))* log10((grid%kappaAbsRed(i1  , i2+1 , i3   , ilambda))) + &
         ((t1   )  * (1.d0-t2) * (t3   ))* log10((grid%kappaAbsRed(i1+1, i2   , i3+1 , ilambda))) + &
         ((1.d0-t1)  * (1.d0-t2) * (t3   ))* log10((grid%kappaAbsRed(i1  , i2   , i3+1 , ilambda))) + &
         ((t1   )  * (t2   ) * (t3   ))* log10((grid%kappaAbsRed(i1+1, i2+1 , i3+1 , ilambda)))

    if (interpGridKappaAbsRed  > 20.) then
       write(*,*) grid%kappaAbsRed(i1-1:i1+1,i2-1:i2+1,i3-1:i3+1,1)
       write(*,*) "red",i1,i2,i3,t1,t2,t3
       write(*,*) interpGridKappaAbsRed
    endif

    interpGridKappaAbsRed = 10.e0**interpGridKappaAbsRed

    if (abs(t1) > 1.d0) then
       write(*,*) "t1",t1,i1,i2,i3
       stop
    endif
    if (abs(t2) > 1.d0) then
       write(*,*) "t2",t2,i1,i2,i3
       stop
    endif
    if (abs(t3) > 1.d0) then
       write(*,*) "t3",t3,i1,i2,i3
       stop
    endif

  end function interpGridKappaAbsRed

  real function interpGridScalar(scalarGrid,nx,ny,nz,nLambda,i1, i2, i3, iLambda,t1, t2, t3)
    integer :: nx,ny,nz,nLambda
    real :: scalarGrid(1:nx,1:ny,1:nz,1:nLambda)
    integer i1,i2,i3, iLambda
    real :: t1, t2, t3

    interpGridScalar = &
         ((1.d0-t1)  * (1.d0-t2) * (1.d0-t3))* scalargrid(i1  , i2   , i3   , ilambda) + &
         ((t1   )  * (1.d0-t2) * (1.d0-t3))* scalargrid(i1+1, i2   , i3   , ilambda) + &
         ((t1   )  * (t2   ) * (1.d0-t3))* scalargrid(i1+1, i2+1 , i3   , ilambda) + &
         ((1.d0-t1)  * (t2   ) * (t3   ))* scalargrid(i1  , i2+1 , i3+1 , ilambda) + &
         ((1.d0-t1)  * (t2   ) * (1.d0-t3))* scalargrid(i1  , i2+1 , i3   , ilambda) + &
         ((t1   )  * (1.d0-t2) * (t3   ))* scalargrid(i1+1, i2   , i3+1 , ilambda) + &
         ((1.d0-t1)  * (1.d0-t2) * (t3   ))* scalargrid(i1  , i2   , i3+1 , ilambda) + &
         ((t1   )  * (t2   ) * (t3   ))* scalargrid(i1+1, i2+1 , i3+1 , ilambda)

  end  function interpGridScalar

  real function interpGridScalar4(scalarGrid,nx,ny,nz,nLambda,i1, i2, i3, iLambda,t1, t2, t3)
    integer nx,ny,nz,nLambda
    real :: scalarGrid(1:nx,1:ny,1:nz,1:nLambda)
    integer i1,i2,i3, iLambda
    real :: t1, t2, t3


    interpGridScalar4 = &
         ((1.d0-t1)  * (1.d0-t2) * (1.d0-t3))* log10(scalargrid(i1  , i2   , i3   , ilambda)) + &
         ((t1   )  * (1.d0-t2) * (1.d0-t3))* log10(scalargrid(i1+1, i2   , i3   , ilambda)) + &
         ((t1   )  * (t2   ) * (1.d0-t3))* log10(scalargrid(i1+1, i2+1 , i3   , ilambda)) + &
         ((1.d0-t1)  * (t2   ) * (t3   ))* log10(scalargrid(i1  , i2+1 , i3+1 , ilambda)) + &
         ((1.d0-t1)  * (t2   ) * (1.d0-t3))* log10(scalargrid(i1  , i2+1 , i3   , ilambda)) + &
         ((t1   )  * (1.d0-t2) * (t3   ))* log10(scalargrid(i1+1, i2   , i3+1 , ilambda)) + &
         ((1.d0-t1)  * (1.d0-t2) * (t3   ))* log10(scalargrid(i1  , i2   , i3+1 , ilambda)) + &
         ((t1   )  * (t2   ) * (t3   ))* log10(scalargrid(i1+1, i2+1 , i3+1 , ilambda))

    interpGridScalar4 = 10.e0**interpGridScalar4

  end  function interpGridScalar4

  real function interpGridScalar2(scalarGrid,nx,ny,nz,i1, i2, i3, t1, t2, t3)
    integer nx,ny,nz
    real :: scalarGrid(1:nx,1:ny,1:nz)
    integer i1,i2,i3
    real :: t1, t2, t3


    interpGridScalar2 = &
         ((1.d0-t1)  * (1.d0-t2) * (1.d0-t3))* scalargrid(i1  , i2   , i3   ) + &
         ((t1   )  * (1.d0-t2) * (1.d0-t3))* scalargrid(i1+1, i2   , i3   ) + &
         ((t1   )  * (t2   ) * (1.d0-t3))* scalargrid(i1+1, i2+1 , i3   ) + &
         ((1.d0-t1)  * (t2   ) * (t3   ))* scalargrid(i1  , i2+1 , i3+1 ) + &
         ((1.d0-t1)  * (t2   ) * (1.d0-t3))* scalargrid(i1  , i2+1 , i3   ) + &
         ((t1   )  * (1.d0-t2) * (t3   ))* scalargrid(i1+1, i2   , i3+1 ) + &
         ((1.d0-t1)  * (1.d0-t2) * (t3   ))* scalargrid(i1  , i2   , i3+1 ) + &
         ((t1   )  * (t2   ) * (t3   ))* scalargrid(i1+1, i2+1 , i3+1 )

  end  function interpGridScalar2

  real function interpGridScalar3(scalarGrid,nx,ny,nz,i1, i2, i3, t1, t2, t3)
    integer nx,ny,nz
    real :: scalarGrid(1:nx,1:ny,1:nz)
    integer i1,i2,i3
    real :: t1, t2, t3


    interpGridScalar3 = &
         ((1.d0-t1)  * (1.d0-t2) * (1.d0-t3))* log10(scalargrid(i1  , i2   , i3   )) + &
         ((t1   )  * (1.d0-t2) * (1.d0-t3))* log10(scalargrid(i1+1, i2   , i3   )) + &
         ((t1   )  * (t2   ) * (1.d0-t3))* log10(scalargrid(i1+1, i2+1 , i3   )) + &
         ((1.d0-t1)  * (t2   ) * (t3   ))* log10(scalargrid(i1  , i2+1 , i3+1 )) + &
         ((1.d0-t1)  * (t2   ) * (1.d0-t3))* log10(scalargrid(i1  , i2+1 , i3   )) + &
         ((t1   )  * (1.d0-t2) * (t3   ))* log10(scalargrid(i1+1, i2   , i3+1 )) + &
         ((1.d0-t1)  * (1.d0-t2) * (t3   ))* log10(scalargrid(i1  , i2   , i3+1 )) + &
         ((t1   )  * (t2   ) * (t3   ))* log10(scalargrid(i1+1, i2+1 , i3)+1 )


    interpGridScalar3 = 10.e0**interpGridScalar3
  end  function interpGridScalar3




  type (VECTOR) function interpGridVelocity(grid,  i1, i2, i3, t1, t2, t3)
    implicit none
    type(GRIDTYPE) :: grid
    integer i1,i2,i3
    integer j1,j2,j3
    real :: t1, t2, t3


    if (grid%cartesian) then
       j1 = min(i1+1,grid%nx)
       j2 = min(i2+1,grid%ny)
       j3 = min(i3+1,grid%nz)
    else
       j1 = min(i1+1,grid%nr)
       j2 = min(i2+1,grid%nmu)
       j3 = min(i3+1,grid%nphi)
    endif

    interpGridVelocity = &
         ((1.d0-t1)  * (1.d0-t2) * (1.d0-t3)) * grid%velocity(i1  , i2   , i3   ) + &
         ((t1   )  * (1.d0-t2) * (1.d0-t3))* grid%velocity(j1, i2   , i3   ) + &
         ((t1   )  * (t2   ) * (1.d0-t3))* grid%velocity(j1, j2 , i3   ) + &
         ((1.d0-t1)  * (t2   ) * (t3   ))* grid%velocity(i1  , j2 , j3 ) + &
         ((1.d0-t1)  * (t2   ) * (1.d0-t3))* grid%velocity(i1  , j2 , i3   ) + &
         ((t1   )  * (1.d0-t2) * (t3   ))* grid%velocity(j1, i2   , j3 ) + &
         ((1.d0-t1)  * (1.d0-t2) * (t3   ))* grid%velocity(i1  , i2   , j3 ) + &
         ((t1   )  * (t2   ) * (t3   ))* grid%velocity(j1, j2 , j3 )



  end  function interpGridVelocity


  ! spline interpolation
  ! this subroutine computes the line of sight directional derivative
  ! numerical - ie completely generalized

  real function directionalDeriv(grid, position, i1, i2, i3 , direction)

    implicit none
    type(GRIDTYPE) :: grid                       ! the opacity grid
    integer :: i1, i2, i3, j1, j2, j3            ! indices
    real :: t1, t2, t3                           ! multipliers
    type(VECTOR) :: direction, position          ! vectors
    type(VECTOR) :: position1, position2
    type(VECTOR) :: rHat
    real :: r, theta, mu, phi                    ! spherical polar coords
    real, parameter  :: h = 1.                  ! factor
    logical :: hit
    real :: dr
    real :: phi1, phi2, dphi, dx
    real :: dy, dz, dtheta

    hit = .false.

    ! dr is a small increment of distance

    if (grid%cartesian) then
       dx = grid%xAxis(i1+1)-grid%xAxis(i1)
       dy = grid%yAxis(i2+1)-grid%yAxis(i2)
       dz = grid%zAxis(i3+1)-grid%zAxis(i3)
       dr = h * min(dx,dy,dz)
    else
       r = grid%rAxis(i1)
       if (i1 /= grid%nR) then
          dr = grid%rAxis(i1+1) - grid%rAxis(i1)
       else
          dr = grid%rAxis(grid%nr) - grid%rAxis(grid%nr-1)
       endif
       if (i2 /= grid%nMu) then
          dtheta = abs(acos(grid%muAxis(i2+1))-acos(grid%muAxis(i2)))
       else
          dtheta = abs(acos(grid%muAxis(grid%nMu))-acos(grid%muAxis(grid%nMu-1)))
       endif
       dr = h * min(dr, r*dtheta)
    endif



    if (modulus(position) /= 0.) then
       rHat = position/modulus(position)
    else
       rHat = direction
       
    endif

    j1=i1
    j2=i2
    j3=i3

    call getPolar(position, r, theta, phi)

    ! get a new position a little way back from current position

    position1 = position - dr * direction

    ! this might be inside core or outside grid- 
    ! in which case just use the current
    ! position as the first point

    r = modulus(position1)
    if (.not.grid%cartesian) then
       if ((r < grid%rAxis(1)).or.(r > grid%rAxis(grid%nr))) then
          position1 = position
          hit = .true.
       endif
       call getPolar(position1, r, theta, phi)
       mu = position1%z/r
    else 
       if (outsideGrid(position1, grid)) then
          position1 = position
          hit = .true.
       endif
    endif

    call getIndices(grid, position1, j1, j2, j3, t1, t2, t3)

    ! first line of sight velocity

    phi1 = direction .dot. interpGridVelocity(grid,j1,j2,j3,t1,t2,t3)

    ! now go forward a bit from current position


    position2 = position + dr * direction
    r = modulus(position2)


    ! check we're still inside grid

    if (.not.grid%cartesian) then
       if ((r < grid%rAxis(1)).or.(r > grid%rAxis(grid%nr))) then
          position2 = position
          hit = .true.
       endif
       call getPolar(position2,r, theta, phi)
       mu = position2%z/r
    else
       if (outsideGrid( position2, grid)) then
          position2 = position
          hit = .true.
       endif
    endif

    call getIndices(grid, position2, j1, j2, j3, t1, t2, t3)


    ! the second position l.o.s. velocity

    phi2 = direction .dot. interpGridVelocity(grid,j1,j2,j3,t1,t2,t3)

    dx = modulus(position2 - position1)

    dphi = phi2-phi1

    ! the line of sight velocity gradient

    if (dx /=0. ) then
       directionalDeriv = abs(dphi / dx)
    else
       directionalDeriv = 1.e-10
    endif

    if (directionalDeriv == 0.) directionalDeriv = 1.e-10

  end function directionalDeriv


  
  subroutine computeCoreEmissionProfile(xArray, sourceSpectrum, nLambda, &
       lamLine, velFWHM, relInt)
    integer :: nLambda, i
    real :: fac1, fac2
    real :: xArray(nLambda)
    real :: sourceSpectrum(nLambda)
    real :: lamLine, velFWHM, relInt, sigma
    real :: thisVel

    sigma = velFWHM / 2.35
    
    do i = 1, nLambda
       thisVel = cSpeed * (xArray(i) - lamLine)/lamLine
       fac1 = relInt - 1.d0
       fac2 = -(thisVel**2 / (2.*sigma**2))
       sourceSpectrum(i) = 1.d0 + fac1 * exp(fac2)
    enddo
  end subroutine computeCoreEmissionProfile


  subroutine computeProbDist(grid, totalLineEmission, totalContEmission, lambda0, useBias)

    type(GRIDTYPE) :: grid
    real :: totalLineEmission, totalContEmission, lambda0
    integer :: i,j,k, ierr
    real, allocatable :: chi(:,:,:)
    logical :: useBias

    
    if (grid%cartesian) then
       allocate(chi(1:grid%nx,1:grid%ny,1:grid%nz), stat=ierr)
       if (ierr /=0) then
         write(*,'(a)') "! Cannot allocate tmp chi memory"
         stop
       endif
       chi = grid%chiLine
    else
       allocate(chi(1:grid%nr,1:grid%nmu,1:grid%nphi), stat=ierr)
       if (ierr /=0) then
         write(*,'(a)') "! Cannot allocate tmp chi memory"
         stop
       endif
       chi = grid%chiLine
    endif


    call computeProbDist2(grid, grid%cartesian, grid%xProbDistLine, grid%xAxis, grid%nx, grid%yProbDistLine, &
         grid%yAxis, grid%ny, grid%zProbDistLine, grid%zAxis, grid%nz, grid%rProbDistLine, grid%rAxis, &
         grid%nr, grid%muProbDistLine, grid%muAxis, grid%nMu, &
         grid%phiProbDistLine, grid%phiAxis, grid%nphi, &
         grid%etaLine, chi, grid%biasLine, totalLineEmission, .true.,lambda0, useBias)


    if (grid%cartesian) then
       do i = 1, grid%nx
          do j = 1, grid%ny
             do k = 1, grid%nz
                chi(i,j,k) = grid%kappaAbs(i,j,k,1)
             enddo
          enddo
       enddo
    else
       do i = 1, grid%nr
          do j = 1, grid%nmu
             do k = 1, grid%nphi
                chi(i,j,k) = grid%kappaAbs(i,j,k,1)
             enddo
          enddo
       enddo
    endif

    call computeProbDist2(grid, grid%cartesian, grid%xProbDistCont, grid%xAxis, grid%nx, grid%yProbDistCont, &
         grid%yAxis, grid%ny, grid%zProbDistCont, grid%zAxis, grid%nz, grid%rProbDistCont, grid%rAxis, &
         grid%nr, grid%muProbDistCont, grid%muAxis, grid%nMu, grid%phiProbDistCont, grid%phiAxis, &
         grid%nphi, grid%etaCont, chi, grid%biasCont,totalContEmission, .false., lambda0, useBias)


    deallocate(chi)


  end subroutine computeProbDist

  subroutine computeProbDist2(grid, cartesian, xProbDist, xAxis, nx, yProbDist, yAxis, ny, zProbDist, zAxis, nz, &
       rProbDist, rAxis, nr, muProbDist, muAxis, nMu, phiProbDist, phiAxis, nPhi, eta, chi,  &
       bias, totalEmission, lineEmission, lambda0, useBias)

    type(GRIDTYPE) :: grid
    type(VECTOR) :: rVec
    real :: tot, tot2
    real :: totalEmission
    logical :: lineEmission, useBias
    integer :: i, j, k
    real :: t1, t2, t3
    integer :: i1,i2,i3
    real :: fac, gridStep
    real :: dR, dtheta, dphi, sintheta
    real :: r,phi
    real :: mu
    real :: dx, dy, dz
    real :: eta(:,:,:)
    real :: chi(:,:,:)
    real :: bias(:)
    integer :: nx ,ny, nz, nr, nmu, nphi
    logical :: cartesian
    real :: scaleFac
    real :: tau1,  nu, escProb, lambda0, tausob
    real :: xProbDist(:), yProbDist(:,:), zProbDist(:,:,:)
    real :: rProbDist(:), muProbDist(:,:), phiProbDist(:,:,:)
    real :: xAxis(:), yAxis(:), zAxis(:)
    real :: rAxis(:), muAxis(:), phiAxis(:)
    real, allocatable :: tmp(:),tmp2(:)
    real(kind=doubleKind) :: dV

    bias = 1.d0
    totalEmission = 0.
    gridStep = xAxis(2) - xAxis(1)

    write(*,*) "Computing probability distributions..."


    if (useBias) then
       if (.not.cartesian) then
          bias = 1.d0
          if (lineEmission) then
             do i = nr, 1, -1
                nu = cSpeed / (lambda0*angstromtocm)
                tauSob = chi(i,1,1)*rAxis(i)/nu/modulus(grid%velocity(i,1,1))
                if (tauSob < 0.1) then
                   escProb=1.d0-tauSob*0.5*(  1.d0 - tauSob/3.0* ( 1.d0 - tauSob*0.25*(1.d0 -0.20e0*tauSob)))
                else if (tauSob < 15.e0) then
                   escProb = (1.d0 - exp(-tauSob))/tauSob
                else
                   escProb = 1.d0/tauSob
                endif
                bias(i) = escProb
                bias(i) = max(1.e-7,bias(i))
             enddo
             
!             allocate(tmp(1:nr))
!             allocate(tmp2(1:nr))
!             tmp = 0.e0
!             do i = 1, nr
!                tmp(i) = sqrt((grid%kappaAbs(i,1,1,1)+grid%kappaSca(i,1,1,1))*grid%kappaAbs(i,1,1,1))
!             enddo
!             tau1 = 0.
!             tmp2(nr) = 1.d0
!             do i = nr-1,1,-1
!                tau1 = tau1 + 0.5*(rAxis(i+1)-rAxis(i))*(tmp(i+1)+tmp(i))
!                tmp2(i) = exp(-tau1)
!                tmp2(i) = max(1.e-7,tmp2(i))
!             enddo
!             do i = 1, nr
!                bias(i) = bias(i) * tmp2(i)
!             enddo
!             deallocate(tmp)
!             deallocate(tmp2)

          else

             allocate(tmp(1:nr))
             tmp = 0.
             do i = 1, nr
                tmp(i) = sqrt((grid%kappaAbs(i,1,1,1)+grid%kappaSca(i,1,1,1))*grid%kappaAbs(i,1,1,1))
             enddo
             do i = 1, nr
                tau1 = 1.e-20
                do j = i, nr-1
                   tau1 = tau1 + 0.5*(rAxis(j+1)-rAxis(j))*(tmp(j+1)+tmp(j))
                enddo
                bias(i) = exp(-tau1)
                bias(i) = max(1.e-7,bias(i))
             enddo
             deallocate(tmp)

             bias(1:nr) = 1.

          endif
          
!          bias(1) = bias(2)
!          do i = nr,2,-1
!             bias(i) = 0.5*(bias(i-1)+bias(i))
!          enddo
!          bias(1) = bias(2)
          

       endif

       write(*,*) "done biases..."
    endif

    if (cartesian) then


       xProbDist = 0.
       tot = 0.
       tot2 = 0.
       do i = 2, nx
             if (i > 1) then
                dx = abs(xAxis(i)-xAxis(i-1))
             else
                dx = abs(xAxis(2)-xAxis(1))
             endif
          do j = 1, ny
             if (j > 1) then
                dy = abs(yAxis(j)-yAxis(j-1))
             else
                dy = abs(yAxis(2)-yAxis(1))
             endif
             do k = 1,nz
                if (k > 1) then
                   dz = abs(zAxis(k)-zAxis(k-1))
                else
                   dz = abs(zAxis(2)-zAxis(1))
                endif
                dV = dble(dx)*dble(dy)*dble(dz)
                if (.not.grid%inStar(i,j,k).and.grid%inUse(i,j,k)) then
                   if (lineEmission) then
                      fac = grid%biasLine3D(i,j,k)
                   else
                      fac = grid%biasCont3D(i,j,k)
                   endif
                   tot = tot + dble(eta(i,j,k))*dV*fac
                endif
             enddo
          enddo
          xProbDist(i) =  tot
       enddo

       scaleFac = xProbDist(grid%nx)


!       xProbDist = xProbDist - xProbDist(1)
       if (xProbDist(nx) /= 0.) then
          do i = 1, nx
             xProbDist(i) = xProbDist(i) / xProbDist(nx)
          enddo
       endif



       yProbDist = 0.
       do i = 1, nx
          if (i > 1) then
             dx = abs(xAxis(i)-xAxis(i-1))
          else
             dx = abs(xAxis(2) - xAxis(1))
          endif
          tot = 0.
          do j = 2, ny
             if (j > 1) then
                dy = abs(yAxis(j)-yAxis(j-1))
             else
                dy = abs(yAxis(2)-yAxis(1))
             endif
             do k = 1, nz
                if (k > 1) then
                   dz = abs(zAxis(k)-zAxis(k-1))
                else
                   dz = abs(zAxis(2)-zAxis(1))
                endif
                dV = dble(dx)*dble(dy)*dble(dz)
                if (.not.grid%inStar(i,j,k).and.grid%inUse(i,j,k)) then
                   if (lineEmission) then
                      fac = grid%biasLine3D(i,j,k)
                   else
                      fac = grid%biasCont3D(i,j,k)
                   endif
                   tot = tot + dble(eta(i,j,k))*dV*fac
                endif
             enddo
             yProbDist(i,j) =  tot
          enddo
       enddo
!       do i = 1, nx
!          yProbDist(i,1:ny) = yProbDist(i,1:ny) - yProbDist(i,1)
!       enddo


       do i = 1, nx
          if (yProbDist(i,ny) /= 0.) then
             do j = 1, ny
                yProbDist(i,j) = yProbDist(i,j)/yProbDist(i,ny)
             enddo
          endif
       enddo

       zProbDist = 0.
       do i = 1, nx
          if (i > 1) then
             dx = abs(xAxis(i)-xAxis(i-1))
          else
             dx = abs(xAxis(2) - xAxis(1))
          endif
          do j = 1,ny
             if (j > 1) then
                dy = abs(yAxis(j)-yAxis(j-1))
             else
                dy = abs(yAxis(2)-yAxis(1))
             endif
             tot2 = 0.
             do k = 2, nz
                if (k > 1) then
                   dz = abs(zAxis(k)-zAxis(k-1))
                else
                   dz = abs(zAxis(2)-zAxis(1))
                endif
                tot = 0.
                dV = dble(dx)*dble(dy)*dble(dz)
                if (.not.grid%inStar(i,j,k).and.grid%inUse(i,j,k)) then
                   tot = tot + dble(eta(i,j,k))*dv
                endif
                if (lineEmission) then
                   fac = grid%biasLine3D(i,j,k)
                else
                   fac = grid%biasCont3D(i,j,k)
                endif

                tot2 = tot2 + tot*fac
                zProbDist(i,j,k) = tot2
                totalEmission = totalEmission + tot
             enddo
          enddo
       enddo
!       do i = 1, nx
!          do j = 1, ny
!             zProbDist(i,j,1:nz) = zProbDist(i,j,1:nz)  - zProbDist(i,j,1)
!          enddo
!       enddo
       
       do i = 1, nx
          do j = 1, ny
             do k = 1, nz
                if (zProbDist(i,j,nz) /= 0.) then
                   zProbDist(i,j,k) = zProbDist(i,j,k) / zProbDist(i,j,nz)
                endif
             enddo
          enddo
       enddo
       if (scaleFac /= 0.) then
          write(*,*) "Bias correction: ",totalEmission/scaleFac
          if (lineEmission) then
             where(grid%inUse)
                grid%biasLine3d = grid%biasLine3d * totalEmission/scaleFac
             end where
          else
             where(grid%inUse)
                grid%biasCont3D = grid%biasCont3d * totalEmission/scaleFac
             end where
          endif
       endif

    else


       write(*,*) "starting r"
       rProbDist = 0.


       do i = 2,nr
          tot = 0.
          do j = 1, nMu-1
             do k = 1, nPhi-1
                dTheta = acos(muAxis(j+1))-acos(muAxis(j))
                mu = 0.5*(muAxis(j+1) + muAxis(j))
                dPhi = phiAxis(k+1) - phiAxis(k)
                sinTheta = sqrt(1.d0 - mu**2)
!                r = 0.5*(rAxis(i) + rAxis(i-1))
                r = rAxis(i)
                phi = 0.5*(phiAxis(k+1) + phiAxis(k))
                dr = rAxis(i) - rAxis(i-1)
                dV = abs(sinTheta * dr * dTheta * dPhi)*r**2
                rVec = VECTOR(r*sinTheta*cos(phi),r*sinTheta*sin(phi),r*mu)

                call getIndices(grid, rVec, i1, i2, i3, t1, t2, t3)

                fac = interpGridScalar3(eta,nr,nmu,nphi,i1, i2, i3, t1, t2, t3)
                tot = tot +  fac * dV * &
                     logint(r,rAxis(i1),rAxis(i1+1),bias(i1),bias(i1+1))

                totalEmission = totalEmission  + fac*dV

             enddo
          enddo
          rProbDist(i) = rProbDist(i-1) + tot 
       enddo

       scaleFac = rProbDist(nr)

       if (rProbDist(nr) > 0.) then
          rProbDist = rProbDist / rProbDist(nr)
       endif
       write(*,*) "r done..."

       muProbDist = 0.

       do i = 1, nr
             tot = 0.
          do j = 2, nMu
             if (j > 1) then
                dTheta = acos(muAxis(j))-acos(muAxis(j-1))
             else
                dTheta = acos(muAxis(2))-acos(muAxis(1))
             endif
             mu = muAxis(j)
             sinTheta = sqrt(1.d0 - mu**2)
             dr = rAxis(2) - rAxis(1)
             r = rAxis(i)


             do k = 1, nPhi-1
                dPhi = phiAxis(k+1) - phiAxis(k)
                phi = 0.5*(phiAxis(k+1) - phiAxis(k))
                if (i > 1) then
                   dr = rAxis(i) - rAxis(i-1)
                else
                   dr = rAxis(2) - rAxis(1)
                endif
                dV = abs(sinTheta * dr * dTheta * dPhi)*r**2
                rVec = VECTOR(r*sinTheta*cos(phi),r*sinTheta*sin(phi),r*mu)

                call getIndices(grid, rVec, i1, i2, i3, t1, t1, t3)


!                tot = tot + interpGridScalar3(eta,nr,nmu,nphi,i1, i2, i3, t1, t2, t3)*dV
                tot = tot + eta(i,j,k)

             enddo
             muProbDist(i,j) = tot
          enddo 
!         do j = nMu, 2, -1
!             muProbDist(i,j) = muProbDist(i,j-1)
!          enddo

       enddo



       do i = 1, nr
          if (muProbDist(i,nMu) > 0.) then
             muProbDist(i,1:nMu) = &
                  muProbDist(i,1:nMu)/muProbDist(i,nMu)
          endif
       enddo
       write(*,*) "mu done..."


       phiProbDist = 0.

       do i = 1, nr
          r = rAxis(i)
          do j = 1, nMu
             if (j > 1) then
                dTheta = abs(acos(muAxis(j))-acos(muAxis(j-1)))
             else
                dTheta = abs(acos(muAxis(2))-acos(muAxis(1)))
             endif
             mu = muAxis(j)
             sinTheta = max(1.d-20,sqrt(1.d0 - muAxis(j)**2))
             do k = 2, nPhi
                dPhi = phiAxis(k) - phiAxis(k-1)
                if (i > 1) then
                   dr = rAxis(i) - rAxis(i-1)
                else
                   dr = rAxis(2) - rAxis(1)
                endif
                phi = 0.5*(phiAxis(k) + phiAxis(k-1))
                dV = abs( sinTheta * dr * dTheta * dPhi)*r**2
                rVec = VECTOR(r*sinTheta*cos(phi),r*sinTheta*sin(phi),r*mu)

                call getIndices(grid, rVec, i1, i2, i3, t1, t1, t3)

                phiProbDist(i,j,k) = phiProbDist(i,j,k-1)  &
                     + interpGridScalar3(eta,nr,nmu,nphi,i1, i2, i3, t1, t2, t3)*dV

             enddo
          enddo
       enddo


       do i = 1, nr
          do j = 1, nMu
             if (phiProbDist(i,j,nPhi) > 0.) then
                phiProbDist(i,j,1:nPhi) = &
                     PhiProbDist(i,j,1:nPhi)/phiProbDist(i,j,nPhi)
             endif
          enddo
       enddo



       write(*,*) "phi done..."


       if (useBias) then
          do i = 1, nr
             bias(i) = bias(i) * totalEmission/scaleFac
          enddo
          if (lineEmission) then
             open(21,file="prob.dat",form="formatted",status="unknown")
             do i = 1, nr
                write(21,*) i,rAxis(i),rProbDist(i),bias(i)*scaleFac/totalEmission
             enddo
             close(21)
          endif
       endif


    endif
    write(*,*) "Distributions done"
  end subroutine computeProbDist2

  type(VECTOR) function thermalElectronVelocity(temperature)
    real :: temperature
    type(VECTOR) :: rHat
    real :: vel
    real :: sigmaVel

!    rHat = randomUnitVector()    
!    sigmaVel = 100.e0*sqrt( kConst * temperature / (mElectron/1000.))
!    vel = gasdev()*sigmaVel/cSpeed


    thermalElectronVelocity = maxwellianVelocity(mElectron,temperature)/cSpeed


  end function thermalElectronVelocity

  type(VECTOR) function thermalHydrogenVelocity(temperature)
    real :: temperature
    type(VECTOR) :: rHat
    real :: vel
    real :: sigmaVel
 
    thermalHydrogenVelocity = maxwellianVelocity(mHydrogen, temperature)/cSpeed


  end function thermalHydrogenVelocity





end module math_mod




