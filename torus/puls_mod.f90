module puls_mod

  use gridtype_mod
  use grid_mod
  use stateq_mod
  use blob_mod

contains

  subroutine fillGridPuls(grid, mDot, rStar, tEff, v0, vterm, beta, &
       xfac, blobs, maxBlobs, doBlobs, vContrast)

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
    real :: w, qfrac,beta2,vext
    integer ::m
    real :: v, b2, b3, xfac, pfac, kfac, beta1, v1
    integer :: maxBlobs
    logical :: doBlobs
    real :: mdotcone, mdotDisk, vTermCone, vTermDisk, betaCone, betaDisk
    real :: coneAng, diskAng
    real :: vContrast
    integer :: muCone1, muCone2, muDisk1, muDisk2
    type(GRIDTYPE) :: grid
    type(BLOBTYPE)  :: blobs(:)

!
! Departure coeffs from Puls et al 1996, A&A, 305, 171
!

    pfac = 2.
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
    
!    mdotCone = 4.e-5
!    mdotCone = mdotCone * msol / (365.25 * 24. * 60. * 60.)
!    mdotDisk = 2.e-6
!    mdotDisk = mdotDisk * msol / (365.25 * 24. * 60. * 60.)
!    vTermCone = 1000. * 1.e5
!    vTermDisk = 100.* 1.e5
!    betaCone = 2.
!    betaDisk = 2.
!    coneAng = 30.*degtorad
!    diskAng = 40.* degtoRad
!    call locate(grid%muAxis, grid%nMu, cos(pi-coneAng/2.), muCone1)
!    call locate(grid%muAxis, grid%nMu, cos(coneAng/2.), muCone2)
!    call locate(grid%muAxis, grid%nMu, cos(pi/2.+diskAng/2.), muDisk1)
!    call locate(grid%muAxis, grid%nMu, cos(pi/2.-diskAng/2.), muDisk2)
    
    grid%etaLine = 1.e-30


    do i = 1, grid%nr
       do k = 1, grid%nPhi
          do j = 1, grid%nMu
             sinTheta = sqrt(1. - grid%muAxis(j)**2)
             rVec = VECTOR(grid%rAxis(i)*cos(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*sin(grid%phiAxis(k))*sinTheta, &
                  grid%rAxis(i)*grid%muAxis(j))
             rHat = rVec
             call normalize(rHat)
!             v1 = vTerm * (1.+(vContrast-1.)*abs(grid%muAxis(j)))/vContrast
             vel = v0 + (vTerm - v0) * (1. - grid%rAxis(1)/grid%rAxis(i))**beta
             grid%velocity(i,j,k) = (vel / cSpeed) * rHat
             grid%rho(i,j,k) = mDot / (fourPi * vel * grid%rAxis(i)**2 * 1.e20)

             grid%rho(i,j,k) = grid%rho(i,j,k)* kfac*(1.-xfac*abs(grid%muAxis(j))**pfac)

             grid%inUse(i,j,k) = .true.
          enddo
       enddo
    enddo

    if (doBlobs) then
       call distortGridWithBlobs(grid, maxBlobs, blobs)
    endif


    if (.not.grid%resonanceLine) then

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
                grid%n(i,j,k,1) = grid%n(i,j,k,1) * 1./(0.5*(1.-sqrt(max(0.,(1.-grid%rAxis(1)**2/grid%rAxis(i)**2)))))

                v = modulus(grid%velocity(i,j,k))/modulus(grid%velocity(grid%nr,j,k))
                
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
    endif

    if (grid%resonanceLine) then
       do i = 1, grid%nr
          do j = 1, grid%nMu
             do k = 1, grid%nPhi
                grid%chiLine(i,j,k) = grid%rho(i,j,k)*10.e11
                grid%etaLine(i,j,k) = 1.e-20
                grid%kappaSca(i,j,k,1) = 1.e-20
                grid%kappaAbs(i,j,k,1) = 1.e-20
                w = 0.5 * sqrt(1. - (1. - (grid%rAxis(1)/grid%rAxis(i))**2))
                qfrac = (w/grid%rho(i,j,k))*(grid%rAxis(1)/grid%rAxis(i))
                grid%chiLine(i,j,k) = grid%chiLine(i,j,k) * qfrac
             enddo
          enddo
       enddo
    endif


  end subroutine fillGridPuls

end module puls_mod
