module puls_mod

  use gridtype_mod
  use grid_mod
  use stateq_mod

contains

  subroutine fillGridPuls(grid, mDot, rStar, tEff, v0, vterm, beta, &
       xfac)

    implicit none
    real :: mDot, rStar, tEff, v0, vTerm, beta
    integer :: i, j, k
    real :: sinTheta, fac
    logical :: ok
    type(VECTOR) :: rHat, rVec
    real :: vel
    integer, parameter :: maxLevels = 9
    real :: x, dx
    real(kind=doubleKind) :: phiT, ne1, ne2, ntot
    integer ::m
    real :: v, b2, b3, xfac, pfac, kfac
    type(GRIDTYPE) :: grid

!
! Departure coeffs from Puls et al 1996, A&A, 305, 171
!

    pfac = 2.e0
    kfac = 1.0e0/(1.0e0-xFac/(1.0e0+pFac))


    grid%geometry = "puls"

    do i = 1, grid%nmu
       grid%muAxis(i) = 2.e0*real(i-1)/real(grid%nmu-1)-1.
    enddo

    do i = 1, grid%nPhi
       grid%phiAxis(i) = twoPi*real(i-1)/real(grid%nphi-1)
    enddo

    x = 1.
    dx = 1.e-3
    do i = 1, grid%nr
       grid%rAxis(i) = x * rStar / 1.e10
       x = x + dx
       dx = dx * 1.25
    enddo

    grid%rCore = grid%rAxis(1)
    grid%rStar1 = grid%rCore
    grid%rStar2 = 0.
    grid%starPos1 = VECTOR(0.,0.,0.)
    grid%starPos2 = VECTOR(1.e30,0.,0.)
    grid%temperature = 0.75 * tEff
    grid%lineEmission = .true.
    
    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             sinTheta = sqrt(1. - grid%muAxis(j)**2)
             rVec = VECTOR(grid%rAxis(i)*cos(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*sin(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*grid%muAxis(j))
             rHat = rVec
             call normalize(rHat)
             vel = v0 + (vTerm - v0) * (1. - grid%rAxis(1)/grid%rAxis(i))**beta
             grid%velocity(i,j,k) = (vel / cSpeed) * rHat
             grid%rho(i,j,k) = mDot / (fourPi * vel * grid%rAxis(i)**2 * 1.e20)
          enddo
       enddo
    enddo


    do i = 1,grid%nr-1
       write(*,*) i , grid%rAxis(i)/grid%rAxis(1), grid%rAxis(i+1)-grid%rAxis(i),&
            modulus(grid%velocity(i,1,1))*cSpeed/1.e5
    enddo


    allocate(grid%n(1:grid%na1, 1:grid%na2, 1:grid%na3,1:maxlevels))
    allocate(grid%ne(1:grid%na1, 1:grid%na2, 1:grid%na3))
    allocate(grid%nTot(1:grid%na1, 1:grid%na2, 1:grid%na3))

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi

             nTot = dble(grid%rho(i,j,k) / mhydrogen)

             if (grid%temperature(i,j,k) < 1.e5) then
                phiT = 1.5d0*log10(2.d0*pi*mElectron*kErg*grid%temperature(i,j,k)) - &
                     3.d0*log10(hCgs) + log10(exp(-13.6d0/(kEV*grid%temperature(i,j,k))))

                phiT = 10.d0**phiT
                call solveQuadDble(1.d0, phiT, &
                     -1.d0*phiT*dble(nTot), ne1, ne2, ok)
                grid%ne(i,j,k) = min(max(ne1,ne2),nTot)
                grid%ne(i,j,k) = max(grid%ne(i,j,k),1.d0)
             else
                grid%ne(i,j,k) = dble(ntot)
             endif

             grid%nTot(i,j,k) = nTot 
             do m = 1, maxlevels
                grid%n(i,j,k,m) = boltzsaha(m, grid%ne(i,j,k), &
                     dble(grid%temperature(i,j,k)))
             enddo

             v = (modulus(grid%velocity(i,j,k))*cSpeed)/vTerm

             if ((v >= 0.).and.(v < 0.01)) then
                b2 = 1.0 + ((1.5 - 1.)/0.01) * v
                b3 = 0.9 + ((1.1 - 0.9)/0.01) * v
             endif

             if ((v >= 0.01).and.(v < 0.1)) then
                b2 = 1.5 + ((1.2 - 1.5)/0.09)*(v - 0.01)
                b3 = 1.1 + ((1.1 - 1.1)/0.09)*(v - 0.01)
             endif

             if ((v >= 0.1).and.(v < 1.)) then
                b2 = 1.2 + ((1.3 - 1.2)/0.9)*(v - 0.1)
                b3 = 1.1 + ((1.1 - 1.1)/0.9)*(v - 0.1)
             endif

             grid%n(i,j,k,2) = grid%n(i,j,k,2) * b2
             grid%n(i,j,k,3) = grid%n(i,j,k,3) * b3


          enddo
       enddo
    enddo

    call generateOpacities(grid, 2, 3)


    deallocate(grid%ne)
    deallocate(grid%n)
    deallocate(grid%nTot)
    
    grid%etaCont = 1.e-30

    do i = 1, grid%nr
       do j = 1, grid%nMu
          do k = 1, grid%nPhi
             fac = kFac * (1.e0-xfac*abs(grid%muAxis(j))**pFac)
             grid%chiLine(i,j,k) = grid%chiLine(i,j,k)*fac**2
             grid%etaLine(i,j,k) = grid%etaLine(i,j,k)*fac**2
             grid%kappaSca(i,j,k,1) = grid%kappaSca(i,j,k,1)*fac
          enddo
       enddo
    enddo


  end subroutine fillGridPuls

end module puls_mod
