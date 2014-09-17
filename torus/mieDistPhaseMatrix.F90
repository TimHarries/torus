module mieDistPhaseMatrix_mod

   type beshEntry
      real :: x
      integer :: nci
      complex, allocatable :: hankel(:)
   end type beshEntry

   type sphereEntry
      complex, allocatable :: f(:)
      complex, allocatable :: g(:)
      real, allocatable :: cnrm(:)
   end type sphereEntry

contains

   subroutine mieDistPhaseMatrixWrapper(nDustType, nLambda, nMuMie, xArray, mReal, mImg, miePhase)
      use phasematrix_mod
#ifdef MPI
      use mpi
#endif
      use inputs_mod, only : aMin, aMax, a0, qDist, pDist
      implicit none

      integer, intent(in) :: nDustType, nLambda, nMuMie
      real, intent(in) :: xArray(:)
      real, intent(in) :: mReal(:,:), mImg(:, :)
      type(PHASEMATRIX), intent(out) :: miePhase(:,:,:)!nDustType, nLambda, nMuMie)

      integer :: i, j, k

      ! Dust distribution parameters
      integer, parameter :: nDist = 100
      real :: a(nDist-1), da(nDist-1), dist(nDist-1)
      real :: aFac

      ! Variables used in calculating various precomputed tables
      integer :: n
      type(beshEntry) :: beshTable(nLambda, nDist-1)
      type(sphereEntry) :: sphereTable(nDist-1)
      real :: x
      integer :: nc, nci
      integer, save :: max_nci = 0
      integer :: ilam_beg, ilam_end
      real :: cosTheta
      real, allocatable :: pnmllg(:,:)
#ifdef MPI
    real, allocatable :: temp(:,:,:,:), tempArray(:), tempArray2(:)
    integer :: np, n_rmdr, m, ierr, i1, i2
#endif



    max_nci = 0
      do i = 1, nDustType
         ! First set up precomputed tables
         call mapDustDistribution(aMin(i),aMax(i),a0(i),qDist(i),pDist(i), nDist, &
                 a, da, dist, aFac)

         do j = 1, nLambda
            do n = 1, nDist-1
               ! convert a(um) and xArray(A) to cm
               x = max(1.e-5, 2.*real(pi)*(a(n)*1.e-4)/(xArray(j)*1.e-8))
               x = min(2.*real(pi)*1000./0.1, x)
               ! xc = x+4.05*x**.3333+2.0                         !               eq 4.16
               ! nci = int(xc)+1
               nci = int(x+4.05*x**.3333+2.0 )+1
               if (nci > max_nci) max_nci = nci

               beshTable(j,n)%x = x
               beshTable(j,n)%nci = nci

               allocate(beshTable(j,n)%hankel(nci))
               call besh(x, nci, beshTable(j,n)%hankel)
            end do   ! n
         end do   ! j
         ! Calculate the associated Legendre functions
         allocate(pnmllg(nMuMie, max_nci))
         do k = 1, nMuMie
            cosTheta = -1. + 2.*real(k-1)/real(nMumie-1)
!            call genlgp(cosTheta, max_nci, pnmllg(k, 1:max_nci))
            ! In order to match the results of the old code when optimisations
            ! are turned off, we must convert to theta and back again
            call genlgpOld(acos(cosTheta), pnmllg(k, 1:max_nci), max_nci)
         end do

    ilam_beg = 1
    ilam_end = nLambda
#ifdef MPI
    ! Set the range of index for a photon loop used later.     
    np = nThreadsGlobal
    n_rmdr = MOD(nLambda,np)
    m = nLambda/np
          
    if (myRankGlobal .lt. n_rmdr ) then
       ilam_beg = (m+1)*myRankGlobal + 1
       ilam_end = ilam_beg + m
    else
       ilam_beg = m*myRankGlobal + 1 + n_rmdr
       ilam_end = ilam_beg + m -1
    end if
#endif

    do j = ilam_beg, ilam_end
            ! Now set up a precomputed table for sphere for fixed values of
            ! dustType and lambda and varying values of a
            do n = 1, nDist-1
               nc = beshTable(j,n)%nci - 1
               allocate(sphereTable(n)%f(nc))
               allocate(sphereTable(n)%g(nc))
               allocate(sphereTable(n)%cnrm(nc))
               call sphere(beshTable(j,n)%x, nc, cmplx(mReal(i,j), mImg(i,j)), &
                       beshTable(j,n)%hankel, sphereTable(n)%f, sphereTable(n)%g, sphereTable(n)%cnrm)
            end do   ! n

            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP PRIVATE (k, cosTheta) &
            !$OMP SHARED (nMuMie,  beshTable, sphereTable, pnmllg, a, da, dist, afac, miePhase, j, i, max_nci, xArray, mreal, mimg)
            !$OMP DO SCHEDULE(DYNAMIC)
            do k = 1, nMumie
               cosTheta = -1. + 2.*real(k-1)/real(nMumie-1)
               call mieDistPhaseMatrix(cosTheta, nDist, beshTable(j,1:nDist-1), sphereTable, &
                       pnmllg(k,1:max_nci), a, da, dist, aFac, miePhase(i,j,k), mReal(i,j), mImg(i,j), &
                       xArray(j))
!               call mieDistPhaseMatrixOld(aMin(i), aMax(i), a0(i), qDist(i), pDist(i), xArray(j), &
!                       cosTheta, miePhase(i,j,k), mReal(i,j), mImg(i,j))
            enddo
            !$OMP END DO
            !$OMP END PARALLEL


            call normalizeMiePhase(miePhase(i,j,1:nMuMie), nMuMie)

            do n = 1, nDist-1
               deallocate(sphereTable(n)%f)
               deallocate(sphereTable(n)%g)
               deallocate(sphereTable(n)%cnrm)
            end do   ! n
         end do   ! j

#ifdef MPI                
                allocate(temp(1:nlambda,1:nMuMie,1:4,1:4))
                temp = 0.
                do j = iLam_beg, iLam_end
                   do k = 1, nMuMie
                      do i1 = 1, 4
                         do i2 = 1, 4
                            temp(j,k,i1,i2)= miePhase(i,j,k)%element(i1,i2)
                         enddo
                      enddo
                   enddo
                enddo
                allocate(tempArray(1:(nLambda*nMuMie*4*4)))
                allocate(tempArray2(1:(nLambda*nMuMie*4*4)))
                tempArray = reshape(temp, shape(tempArray))

                call MPI_ALLREDUCE(tempArray, tempArray2, size(tempArray), MPI_REAL,&
                     MPI_SUM, MPI_COMM_WORLD, ierr)
                temp(:,:,:,:) = reshape(tempArray2, shape(temp))
                do j = 1, nLambda
                   do k = 1, nMuMie
                      do i1 = 1, 4
                         do i2 = 1, 4
                            miePhase(i,j,k)%element(i1,i2) = temp(j,k,i1,i2)
                         enddo
                      enddo
                   enddo
                enddo
                deallocate(temp, temparray, temparray2)
#endif



         deallocate(pnmllg)
         do j = 1, nLambda
            do n = 1, nDist-1
               deallocate(beshTable(j,n)%hankel)
            end do
         end do
         if (writeoutput) write(*,*) "Dust type ",i," done."
      end do   ! i

   end subroutine mieDistPhaseMatrixWrapper

   subroutine mapDustDistribution(aMin, aMax, a0, qDist, pDist, nDist, a, da, dist, aFac)
      use constants_mod
      implicit none

      real, intent(in) :: aMin, aMax, a0, qDist, pDist
      integer, intent(in) :: nDist
      real, intent(out) :: a(nDist-1), da(nDist-1), dist(nDist-1)
      real, intent(out) :: aFac

      integer :: i
      real :: logAmin, logAmax
      real :: loga1, loga2
      real :: a1, a2
      real :: p

      if (aMin == aMax) then
         logAmin = log(aMin)
         logAmax = log(aMax*1.01)
      else
         logAmin = log(aMin)
         logAmax = log(aMax)
      endif
      aFac = 0.d0
      p = (logAmax-logAmin)/real(nDist-1)
      do i = 1, nDist-1
!         a1 = exp(logAmin + p*real(i-1))
!         a2 = exp(logAmin + p*real(i))
         ! The a1 and a2 calculations have to be done in two steps in order to
         ! give identical results to the original code when optimisations are off
         loga1 = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
         loga2 = logAmin+(logAmax-logAmin)*real(i)/real(nDist-1)
         a1 = exp(loga1)
         a2 = exp(loga2)
         a(i) = 0.5*(a1+a2)
         da(i) = a2 - a1
         ! a is in microns - convert to cm by multiplying by 1.e-4
         ! gfac was only used to calculate t, which is not itself used
!         gfac(i) = pi*(a(i) * 1.e-4)**2
         dist(i) = a(i)**(-qDist)*exp(-(a(i)/a0)**pDist)

         aFac = aFac + 0.5*(a1**(-qDist)*exp(-(a1/a0)**pDist) &
                + a2**(-qDist)*exp(-(a2/a0)**pDist))*(a2-a1)
      enddo
      aFac = 1./aFac
   end subroutine mapDustDistribution

   subroutine mieDistPhaseMatrix(cosTheta, nDist, beshTable, sphereTable, pnmllg, &
                 aArray, da, dist, aFac, miePhase, mReal, mImg, lambda)
      use phasematrix_mod
      use bhmie_mod, only: MXNANG, bhmie
      implicit none
      real, parameter :: micronsToCm=1.e-4
      real, intent(in) :: cosTheta
      real :: mReal, mImg, lambda, qExt, qSca, qBack, gSca
      integer, intent(in) :: nDist
      type(beshEntry), intent(in) :: beshTable(:) !1:nDist-1)
      type(sphereEntry), intent(in) :: sphereTable(:)!nDist-1)
      real, intent(in) :: pnmllg(:)
      real, intent(in) :: aArray(nDist-1), da(nDist-1), dist(nDist-1)
      real, intent(in) :: aFac
      type(PHASEMATRIX), intent(out) :: miePhase

!     ..................................................................
!     .  Light Scattering by Particles: Computational Methods          .
!     .  by P.W. Barber and S.C. Hill                                  .
!     .  copyright (c) 1990 by World Scientific Publishing Co Pte Ltd  .
!     .                                                                .
!     .  equation numbers in columns 73-80 are references to the text  .
!     ..................................................................
!     ...........................................................
!     .  calculate the elements of the scattering matrix        .
!     .    for a sphere                                         .
!     .  inputs: x = size parameter (ka)                        .
!     .          cm = complex index of refraction, (real,imag)  .
!     .               (imag is positive for absorption)         .
!     .          dlt = increment in scattering angle, degrees   .
!     .                                                         .
!     .  dimension of arrays f(*), g(*), amat(*) and cnrm(*):   .
!     .    nc = int(x+4.05*x**.3333+2.0), e.g., for x = 200,    .
!     .    nc = 225                                             .
!     .                                                         .
!     .  dimension of arrays bj(*), by(*), hkl(*), and          .
!     .    pnmllg(*): nc+1, e.g., for x = 200, nc+1 = 226       .
!     .                                                         .
!     .  arrays are set for a maximum size parameter of 200     .
!     ...........................................................

      integer :: i, n, n1
      complex :: cim
      real :: p1, p2
      complex :: s1, s2
      real :: p11, pl, p33p11, p34p11
      real :: p11tot, pltot, p33p11tot, p34p11tot
      real :: a, x, fac,totfac
      complex :: refrel
      complex,save :: Ci = (0.0,1.0)
      complex s12(2*MXNANG-1),s22(2*MXNANG-1)
      integer :: nang=2
!$OMP THREADPRIVATE (ci)

      p11tot = 0.
      pltot = 0.
      p33p11tot = 0.
      p34p11tot = 0.
      totfac = 0.
      miePhase%element = 0.
     do i = 1, nDist-1




          a = aArray(i) 
          x = 2.*pi*(a * micronsTocm)/(lambda*1.e-8)
          x = max(x, 1.e-5)
          refrel = cmplx(mReal, mimg)
!          call BHMIE(X,REFREL,NANG ,S12,S22,QEXT,QSCA,QBACK, GSCA)

!          qSca = qSca * pi * (a * micronsToCm)**2







!        .............................................................
!        .  calculate t = Qsca*x**2                                  .
!        .              = (scattering cross section)*x**2/(pi*a**2)  .
!        .............................................................
!        ! t is not used in any expression, so is no longer calculated --chris
!         t = 0.0
!         do  n = 1,nc
!           t = t+(abs(f(n))**2+abs(g(n))**2)*cnrm(n)
!         enddo
!         t = gfac*t/x**2
         s1 = 0.0
         s2 = 0.0
         do n = 1, beshTable(i)%nci-1
            n1 = n+1
            cim = ci**(-n1)
            p1 = real(n)*costheta*pnmllg(n1)-(real(n)+1.0)*pnmllg(n)
            p2 = pnmllg(n1)
!           ........................................
!           .  calculate parallel field amplitude  .
!           .    for parallel incident             .                      ! eq 3.50b
!           ........................................                          and
            s2 = s2+cim*(-ci*p2*sphereTable(i)%f(n) &                     ! eq 4.10a
                         +p1*sphereTable(i)%g(n))*sphereTable(i)%cnrm(n)
!           .............................................
!           .  calculate perpendicular field amplitude  .
!           .    for perpendicular incident             .
!           .  use coefficients for parallel incident   .
!           .    f(e1n) = -f(o1n),  g(o1n) = g(ein)     .                 ! eq 3.50a
!           .............................................                     and
            s1 = s1+cim*(ci*p1*(-sphereTable(i)%f(n)) &                   ! eq 4.11b
                         +p2*sphereTable(i)%g(n))*sphereTable(i)%cnrm(n)
         enddo
!       ....................................................
!       .  calculate P11, PL = -P12/P11, P33/P11, P34/P11  .              ! eq 4.25
!       ....................................................                  and
         p11 = 2.0*(abs(s1)**2+abs(s2)**2)
         pl = -2.0*(abs(s2)**2-abs(s1)**2)/(p11)
         p33p11 = 4.0*real(s1*conjg(s2))/(p11)
         p34p11 = 4.0*aimag(s2*conjg(s1))/(p11)

         fac = afac * dist(i) * da(i) !* qsca/1.d-15
!         write(*,*) "Fac ",fac,p11,pl,p33p11,p34p11
         p11tot = p11tot + fac*p11
         pltot = pltot + fac*pl
         p33p11tot = p33p11tot + fac*p33p11
         p34p11tot = p34p11tot + fac*p34p11
         totfac = totfac + fac

         miePhase%element(1,1) = miePhase%element(1,1) + p11 * fac
         miePhase%element(1,2) = miePhase%element(1,2) - pl * p11 * fac 
         miePhase%element(2,1) = miePhase%element(2,1) - pl * p11 *fac
         miePhase%element(2,2) = miePhase%element(2,2) + p11 * fac
         miePhase%element(3,3) = miePhase%element(3,3) + p33p11 * p11 * fac
         miePhase%element(3,4) = miePhase%element(3,4) + p34p11 * p11 * fac 
         miePhase%element(4,3) = miePhase%element(4,3) - p34p11 * p11 * fac
         miePhase%element(4,4) = miePhase%element(4,4) + p33p11 * p11 * fac


      enddo
      p11tot = p11tot / totfac
      pltot = pltot / totfac
      p33p11tot = p33p11tot / totfac
      p34p11tot = p34p11tot / totfac

      miePhase%element = miePhase%element / totFac

!      miePhase%element(1,1) = p11tot
!      miePhase%element(1,2) = -pltot*p11tot
!      miePhase%element(2,1) = -pltot*p11tot
!      miePhase%element(2,2) = p11tot

!      miePhase%element(3,3) = p33p11tot*p11tot
!      miePhase%element(3,4) = p34p11tot*p11tot
!      miePhase%element(4,3) = -p34p11tot*p11tot
!      miePhase%element(4,4) = p33p11tot*p11tot
   end subroutine mieDistPhaseMatrix

   subroutine sphere(x, nc, cm, hkl, f, g, cnrm)
      implicit none
!     ..............................................................
!     .  calculate the scattered field f(n) and g(n) coefficients  .
!     .    the f(n) and g(n) for theta incident polarization are,  .
!     .    within an n-dependent factor, the same as the b(n) and  .
!     .    a(n) coefficients, respectively, defined in C.F.        .
!     .    Bohren and D.R. Huffman, Absorption and Scattering of   .
!     .    Light by Small Particles (Wiley- Interscience,New       .
!     .    York,1983), p.100                                       .
!     ..............................................................
      real, intent(in) :: x
      integer, intent(in) :: nc
      complex, intent(in) :: cm
      complex, intent(in) :: hkl(nc+1)
      complex, intent(out) :: f(nc), g(nc)
      real, intent(out) :: cnrm(nc)

      complex, allocatable :: amat(:)
      complex :: b,z,ci,an
      real :: xc, bj, bjm, rf, rn
      integer :: nci, n, nmx

      ci = (0.0,1.0)
!     ......................................................
!     .  set the number of terms required for convergence  .
!     ......................................................
      xc = x+4.05*x**.3333+2.0                                            ! eq 4.16
!     nc = int(xc)
      nci = nc+1
      z = cm*x
!     ..................................................
!     .  logarithmic derivative calculation - set the  .
!     .    starting order for downward recursion       .
!     ..................................................
      nmx = int(max(xc,abs(z)))+15                                        ! eq 4.20
      allocate(amat(nc))

      an = 0.0
      do n = nmx,nc+1,-1
         rn = real(n)
         an = rn/z-1.0/(an+rn/z)
      end do
      amat(nc) = an
      do n = nc,2,-1
         rn = real(n)
         amat(n-1) = rn/z-1.0/(amat(n)+rn/z)                              ! eq 4.19
      end do
!     ...................................................
!     .  calculate the Bessel functions - the order is  .
!     .    incremented by one in the hkl(*) array       .
!     ...................................................
      ! We no longer need to calculate the Bessel functions here because we've
      ! precomputed them and passed them in as an array
!      call besh(x,nci,hkl)
      bj = real(hkl(1))
!     ................................
!     .  calculate the coefficients  .
!     ................................
      do n = 1,nc
         rn = real(n)
         rf = 2.0*rn*(rn+1.0)
         bjm = bj
         bj = real(hkl(n+1))
!        .......................................................
!        .  scattering coefficients for theta                  .
!        .    (parallel) incident polarization                 .
!        .    f(n) = -ci**n*rf*(Bohren and Huffman's b(n))     .
!        .    g(n) = ci**(n+1)*rf*(Bohren and Huffman's a(n))  .
!        .......................................................
         b = cm*amat(n)+rn/x
         f(n) = -ci**n*rf*(b*bj-bjm)/(b*hkl(n+1)-hkl(n))                  ! eq 4.18a
         b = amat(n)/cm+rn/x
         g(n) = ci**(n+1)*rf*(b*bj-bjm)/(b*hkl(n+1)-hkl(n))               ! eq 4.18b
!        ........................................
!        .  calculate the normalization factor  .
!        .    (used in main program)            .
!        ........................................
         cnrm(n) = (2.0*rn+1.0)/(rf*rn*(rn+1.0))
      end do

      deallocate(amat)
   end subroutine sphere

   subroutine besh(x, nc, hankel)
      implicit none
!     ...................................................
!     .  calculate Hankel functions                     .
!     .  bj = Bessel function of the first kind         .
!     .  by = Bessel function of the second kind        .
!     .  x = real argument                              . 
!     .  nc = number of orders (0 to nc-1)              .
!     .  the order of the functions is incremented by   .
!     .    one in the bj(*),by(*) and hankel(*) arrays  .
!     .                                                 .
!     .  arrays are set for nc = 226 maximum            .
!     ...................................................
      real, intent(in) :: x
      integer, intent(in) :: nc
      complex, intent(out) :: hankel(nc)

      real :: a, alpha
      integer :: n, nst
      real, allocatable :: bj(:), by(:)
      real :: t(3)

      allocate(bj(nc))
      allocate(by(nc))
!     ................................................
!     .  by(*) calculation - obtain the zeroeth and  .
!     .                      first order functions   .
!     ................................................
      a = sin(x)/x                                                        ! eq 4.68
      by(1) = -cos(x)/x                                                   ! eq 4.69a
      by(2) = by(1)/x-a                                                   ! eq 4.69b
!     ...........................................................
!     .  obtain the higher order functions by upward recursion  .
!     ...........................................................
      do n = 3,nc
         by(n) = (2.0*real(n-2)+1.0)*by(n-1)/x-by(n-2)
      end do
!     ................................................
!     .  bj(*) calculation - set the starting order  .
!     .                      for downward recursion  .
!     ................................................
      nst = nc+int((101.0+x)**.5)                                         ! eq 4.21
!     ....................................................
!     .  the t(*) array is used to recur down to the     .
!     .    two highest order functions that are needed   .
!     .  set starting values for the two highest orders  .
!     .    nst and nst-1                                 .
!     ....................................................
      t(3) = 0.0
      t(2) = 1.0e-35
!     ...................................................
!     .  recur downward to obtain orders nc-1 and nc-2  .
!     ...................................................
      do n = nst-1,nc-1,-1
         t(1) = (2.0*real(n)+1.0)*t(2)/x-t(3)
         t(3) = t(2)
         t(2) = t(1)
      end do
!     ...............................................
!     .  continue downward recursion to order zero  .
!     ...............................................
      bj(nc) = t(3)
      bj(nc-1) = t(2)
      do n = nc-2,1,-1
         bj(n) = (2.0*real(n)+1.0)*bj(n+1)/x-bj(n+2)
      end do                
!     ..................................................
!     .  calculate the scale factor and the functions  .
!     ..................................................
      alpha = a/bj(1)
      do n = 1,nc
          hankel(n) = cmplx(bj(n)*alpha,by(n))
      end do

      deallocate(bj)
      deallocate(by)
   end subroutine besh

   subroutine genlgp(cosTheta, nc, pnmllg)
      implicit none
!     ........................................................
!     .  calculate associated Legendre functions (argument   .
!     .    cos(theta)) divided by sin(theta) for m = 1       .
!     .  generate first two orders by formula and remaining  .
!     .    orders by recursion                               .
!     .                                                      .
!     .  pnmllg = associated Legendre function/sin(theta)    .
!     .  nc = number of orders (0 to nc-1)                   .
!     .  the order of the associated Legendre functions is   .
!     .    incremented by one in the pnmllg(*) array         .
!     ........................................................
      real, intent(in) :: cosTheta
      integer, intent(in) :: nc
      real, intent(out) :: pnmllg(:)

      integer :: n
      real :: rn

!     ..............................
!     .  calculate orders 0 and 1  .
!     ..............................
      pnmllg(1) = 0.0                                                     ! eq 4.70a
      pnmllg(2) = 1.0                                                     ! eq 4.70b

!     .................................................
!     .  recur upward to obtain all remaining orders  .
!     .................................................
      do n = 3, nc 
         rn = real(n-1)
         pnmllg(n) = ((2.0*rn-1.0)*cosTheta*pnmllg(n-1) - rn*pnmllg(n-2)) / (rn-1.0) ! eq 4.71
      end do 
   end subroutine genlgp

   subroutine normalizeMiePhase(miePhase, nMuMie)
      use phasematrix_mod

      integer, intent(in) :: nMuMie
      type(PHASEMATRIX), intent(inout) :: miePhase(:)
      integer :: i, m ,n
      real(double) :: normfac
  
      normFac = 0.d0
      do i = 2, nMuMie
         normFac = normFac + miePhase(i)%element(1,1)
      enddo
      normfac = normFac / dble(nMuMie-1)
  
      do i = 1, nMuMie
         do m = 1, 4
            do n = 1, 4
               if (normFac /= 0.d0) then
                  miePhase(i)%element(m,n) = real(miePhase(i)%element(m,n) / normFac)
               endif
            enddo
         enddo
      enddo
   end subroutine normalizeMiePhase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The subroutines below here are old and obsolete
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mieDistPhaseMatrixOld(aMin, aMax, a0, qDist, pDist, lambda, &
                                costheta, mieMatrix, cmr, cmi)

      use constants_mod
      use phasematrix_mod

      implicit none


      type(PHASEMATRIX) :: mieMatrix
!      real :: normFac, dMu
      real :: cosTheta
!      integer :: nMu = 20
      real :: aMin, aMax, a0, qDist, pDist, lambda, aFac
      integer :: nDist
      real ::  x, cmr, cmi,  t, theta, costh, rn, p1, p2
      real :: logamax, logamin
      integer :: nc,nci,n, i, n1 !,j

      real :: p11, pl, p33p11, p34p11
      
!     ..................................................................
!     .  Light Scattering by Particles: Computational Methods          .
!     .  by P.W. Barber and S.C. Hill                                  .
!     .  copyright (c) 1990 by World Scientific Publishing Co Pte Ltd  .
!     .                                                                .
!     .  equation numbers in columns 73-80 are references to the text  .
!     ..................................................................
!     ...........................................................
!     .  calculate the elements of the scattering matrix        .
!     .    for a sphere                                         .
!     .  inputs: x = size parameter (ka)                        .
!     .          cm = complex index of refraction, (real,imag)  .
!     .               (imag is positive for absorption)         .
!     .          dlt = increment in scattering angle, degrees   .
!     .                                                         .
!     .  dimension of arrays f(*), g(*), amat(*) and cnrm(*):   .
!     .    nc = int(x+4.05*x**.3333+2.0), e.g., for x = 200,    .
!     .    nc = 225                                             .
!     .                                                         .
!     .  dimension of arrays bj(*), by(*), hkl(*), and          .
!     .    pnmllg(*): nc+1, e.g., for x = 200, nc+1 = 226       .
!     .                                                         .
!     .  arrays are set for a maximum size parameter of 200     .
!     ...........................................................
!      complex cm,ci,cim,f(1000),g(1000),s1,s2
      complex cm,ci,cim,f(1000000),g(1000000),s1,s2 !tjh
      common /cfcom/ f,g,cnrm
!      real :: pnmllg(1001),cnrm(1000)
      real :: pnmllg(1000001),cnrm(1000000)
      real :: a, da, gfac, dist
      real :: a1,a2,loga1,loga2
      real :: tot
      real :: p11tot, pltot, p33p11tot, p34p11tot
      real :: micronsTocm

      nc = 0
      pnmllg = 0

      ci = (0.0,1.0)

      ! for amorphous carbon grains we assume we're in the rayleigh regime

      cm = cmplx(cmr,cmi)


      micronsTocm = 1.e-4

      p11tot = 0.
      pltot = 0.
      p33p11tot = 0.
      p34p11tot = 0.

      nDist = 100
      tot = 0.
      logamin = log(aMin)
      logamax = log(aMax)
      do i = 1, nDist-1
         loga1 = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
         loga2 = logAmin+(logAmax-logAmin)*real(i)/real(nDist-1)
         a1 = exp(loga1)
         a2 = exp(loga2)
         tot = tot + 0.5*(a1**(-qDist)*exp(-(a1/a0)**pDist) &
              + a2**(-qDist)*exp(-(a2/a0)**pDist))*(a2-a1)
      enddo
      aFac = 1./tot

      do i = 1 , nDist-1
        loga1 = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
        loga2 = logAmin+(logAmax-logAmin)*real(i)/real(nDist-1)
        a1 = exp(loga1)
        a2 = exp(loga2)
        da = a2 - a1
        a = 0.5*(a1+a2)
        x = real(2.*pi*(a * micronsToCm)/(lambda*1.e-8))
        x = max(1.e-5, x)
        x = min(2.*real(pi)*1000./0.1, x)
        gfac = real(pi*(a * micronsToCm)**2)

!       .........................................
!       .  set the complex index of refraction  .
!       .    for an exp(-iwt) time variation    .
!       .........................................
        call sphereOld(x,cm,nc)
        nci = nc+1
        t = 0.0
!       .............................................................
!       .  calculate t = Qsca*x**2                                  .
!       .              = (scattering cross section)*x**2/(pi*a**2)  .
!       .............................................................
        do  n = 1,nc
          t = t+(abs(f(n))**2+abs(g(n))**2)*cnrm(n)
        enddo
        t = gfac*t/x**2
        costh = costheta
        theta = acos(costheta)
        call genlgpOld(theta,pnmllg,nci)
        s1 = 0.0
        s2 = 0.0
        do n = 1,nc
          n1 = n+1
          cim = ci**(-n1)
          rn = real(n)
          p1 = rn*costh*pnmllg(n1)-(rn+1.0)*pnmllg(n)
          p2 = pnmllg(n1)
!         ........................................
!         .  calculate parallel field amplitude  .
!         .    for parallel incident             .     !                  eq 3.50b
!         ........................................                          and
          s2 = s2+cim*(-ci*p2*f(n)+p1*g(n))*cnrm(n) !                  eq 4.10a
!         .............................................
!         .  calculate perpendicular field amplitude  .
!         .    for perpendicular incident             .
!         .  use coefficients for parallel incident   .
!         .    f(e1n) = -f(o1n),  g(o1n) = g(ein)     .  !                 eq 3.50a
!         .............................................                     and
          s1 = s1+cim*(ci*p1*(-f(n))+p2*g(n))*cnrm(n) !                eq 4.11b
        enddo
!       ....................................................
!       .  calculate P11, PL = -P12/P11, P33/P11, P34/P11  .              eq 4.25
!       ....................................................              and
         p11 = 2.0*(abs(s1)**2+abs(s2)**2)
         pl = -2.0*(abs(s2)**2-abs(s1)**2)/(p11)
         p33p11 = 4.0*real(s1*conjg(s2))/(p11)
         p34p11 = 4.0*aimag(s2*conjg(s1))/(p11)
         dist = a**(-qDist)*exp(-(a/a0)**pDist)

         p11tot = p11tot + aFac*dist*da*p11
         pltot = pltot-aFac*dist*da*pl
         p33p11tot = p33p11tot + aFac*dist*da*p33p11
         p34p11tot = p34p11tot + aFac*dist*da*p34p11
       enddo

       mieMatrix%element = 0.
       mieMatrix%element(1,1) = p11tot
       mieMatrix%element(1,2) = pltot*p11tot
       mieMatrix%element(2,1) = pltot*p11tot
       mieMatrix%element(2,2) = p11tot

       mieMatrix%element(3,3) = p33p11tot*p11tot
       mieMatrix%element(3,4) = p34p11tot*p11tot
       mieMatrix%element(4,3) = -p34p11tot*p11tot
       mieMatrix%element(4,4) = p33p11tot*p11tot

!
!
!     dMu = 2./real(nMu)
!     normFac = 0.
!     do j = 1, nMu
!
!      costh = (2.*real(j-1)/real(nMu-1))-1.
!      theta = acos(costh)
!
!      p11tot = 0.
!      pltot = 0.
!      p33p11tot = 0.
!      p34p11tot = 0.
!      do i = 1 , nDist-1
!        loga1 = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
!        loga2 = logAmin+(logAmax-logAmin)*real(i)/real(nDist-1)
!        a1 = exp(loga1)
!        a2 = exp(loga2)
!        da = a2 - a1
!        a = 0.5*(a1+a2)
!        x = 2.*pi*(a * micronsToCm)/(lambda*1.e-8)
!        x = max(1.e-5,x)
!        gfac = pi*(a * micronsTocm)**2
!      
!
!!     .........................................
!!     .  set the complex index of refraction  .
!!     .    for an exp(-iwt) time variation    .
!!     .........................................
!      cm = cmplx(cmr,cmi)
!      call sphereOld(x,cm,nc)
!      nci = nc+1
!      t = 0.0
!!     .............................................................
!!     .  calculate t = Qsca*x**2                                  .
!!     .              = (scattering cross section)*x**2/(pi*a**2)  .
!!     .............................................................
!      do  n = 1,nc
!        t = t+(abs(f(n))**2+abs(g(n))**2)*cnrm(n)
!      enddo
!      t = gfac*t/x**2
!      call genlgpOld(theta,pnmllg,nci)
!      s1 = 0.0
!      s2 = 0.0
!      do n = 1,nc
!          n1 = n+1
!          cim = ci**(-n1)
!          rn = real(n)
!          p1 = rn*costh*pnmllg(n1)-(rn+1.0)*pnmllg(n)
!          p2 = pnmllg(n1)
!!     ........................................
!!     .  calculate parallel field amplitude  .
!!     .    for parallel incident             .     !                  eq 3.50b
!!     ........................................                          and
!          s2 = s2+cim*(-ci*p2*f(n)+p1*g(n))*cnrm(n) !                  eq 4.10a
!!     .............................................
!!     .  calculate perpendicular field amplitude  .
!!     .    for perpendicular incident             .
!!     .  use coefficients for parallel incident   .
!!     .    f(e1n) = -f(o1n),  g(o1n) = g(ein)     .  !                 eq 3.50a
!!     .............................................                     and
!          s1 = s1+cim*(ci*p1*(-f(n))+p2*g(n))*cnrm(n) !                eq 4.11b
!        enddo
!!     ....................................................
!!     .  calculate P11, PL = -P12/P11, P33/P11, P34/P11  .              eq 4.25
!!     ....................................................              and
!       p11 = 2.0*(abs(s1)**2+abs(s2)**2)
!       pl = -2.0*(abs(s2)**2-abs(s1)**2)/(p11)
!       p33p11 = 4.0*real(s1*conjg(s2))/(p11)
!       p34p11 = 4.0*aimag(s2*conjg(s1))/(p11)
!       dist = a**(-qDist)*exp(-(a/a0)**pDist)
!       p11tot = p11tot + aFac*dist*da*p11
!       pltot = pltot-aFac*dist*da*pl
!       p33p11tot = p33p11tot + aFac*dist*da*p33p11
!       p34p11tot = p34p11tot + aFac*dist*da*p34p11
!
!    enddo
!
!          normFac = normFac + p11tot*dmu
!       enddo
!
!
!      normFac = normFac * 0.5
!
!      do i = 1, 4
!       do j = 1 , 4
!         mieMatrix%element(i,j) = MieMatrix%element(i,j)/normFac
!       enddo
!      enddo

      end subroutine mieDistPhaseMatrixOld

      subroutine sphereOld(x,cm,nc)
        implicit none
!     ..............................................................
!     .  calculate the scattered field f(n) and g(n) coefficients  .
!     .    the f(n) and g(n) for theta incident polarization are,  .
!     .    within an n-dependent factor, the same as the b(n) and  .
!     .    a(n) coefficients, respectively, defined in C.F.        .
!     .    Bohren and D.R. Huffman, Absorption and Scattering of   .
!     .    Light by Small Particles (Wiley- Interscience,New       .
!     .    York,1983), p.100                                       .
!     ..............................................................
!      complex b,z,cm,ci,hkl(1001),an,amat(1000),f(1000),g(1000)
      complex b,z,cm,ci,hkl(1000001),an,amat(1000000),f(1000000),g(1000000) !tjh
      common /cfcom/ f,g,cnrm
!      dimension cnrm(1000)
      real :: cnrm(1000000) ! tjh, real declaration dma
      real :: x, xc, bj, bjm, rf, rn
      integer :: nc, nci, n, nmx
      hkl = 0
      ci = (0.0,1.0)
!     ......................................................
!     .  set the number of terms required for convergence  .
!     ......................................................
      xc = x+4.05*x**.3333+2.0                         !               eq 4.16
      nc = int(xc)
      nci = nc+1
      z = cm*x
!     ..................................................
!     .  logarithmic derivative calculation - set the  .
!     .    starting order for downward recursion       .
!     ..................................................
      nmx = int(max(xc,abs(z)))+15                      !              eq 4.20
      an = 0.0
      do n = nmx,nc+1,-1
        rn = real(n)
        an = rn/z-1.0/(an+rn/z)
      end do
      amat(nc) = an
      do n = nc,2,-1
        rn = real(n)
        amat(n-1) = rn/z-1.0/(amat(n)+rn/z)              !             eq 4.19
      end do
!     ...................................................
!     .  calculate the Bessel functions - the order is  .
!     .    incremented by one in the hkl(*) array       .
!     ...................................................
      call beshOld(x,hkl,nci)
      bj = real(hkl(1))
!     ................................
!     .  calculate the coefficients  .
!     ................................
      do n = 1,nc
        rn = real(n)
        rf = 2.0*rn*(rn+1.0)
        bjm = bj
        bj = real(hkl(n+1))
!     .......................................................
!     .  scattering coefficients for theta                  .
!     .    (parallel) incident polarization                 .
!     .    f(n) = -ci**n*rf*(Bohren and Huffman's b(n))     .
!     .    g(n) = ci**(n+1)*rf*(Bohren and Huffman's a(n))  .
!     .......................................................
        b = cm*amat(n)+rn/x
        f(n) = -ci**n*rf*(b*bj-bjm)/(b*hkl(n+1)-hkl(n))   !            eq 4.18a
        b = amat(n)/cm+rn/x
        g(n) = ci**(n+1)*rf*(b*bj-bjm)/(b*hkl(n+1)-hkl(n)) !           eq 4.18b
!     ........................................
!     .  calculate the normalization factor  .
!     .    (used in main program)            .
!     ........................................
        cnrm(n) = (2.0*rn+1.0)/(rf*rn*(rn+1.0))
      end do
      end subroutine sphereOld

      subroutine beshOld(x,hankel,nc)
        implicit none
!     ...................................................
!     .  calculate Hankel functions                     .
!     .  bj = Bessel function of the first kind         .
!     .  by = Bessel function of the second kind        .
!     .  x = real argument                              . 
!     .  nc = number of orders (0 to nc-1)              .
!     .  the order of the functions is incremented by   .
!     .    one in the bj(*),by(*) and hankel(*) arrays  .
!     .                                                 .
!     .  arrays are set for nc = 226 maximum            .
!     ...................................................
      real :: x, a, alpha, by, bj, rn, ri, t   
      integer :: nc, i, k, nst, n
      complex hankel(nc)
      integer, parameter :: nmax = 1000001
      dimension bj(nmax),by(nmax),t(3)
!     ................................................
!     .  by(*) calculation - obtain the zeroeth and  .
!     .                      first order functions   .
!     ................................................
      a = sin(x)/x                                        !             eq 4.68
      by(1) = -cos(x)/x                                    !           eq 4.69a
      by(2) = by(1)/x-a                                     !          eq 4.69b
!     ...........................................................
!     .  obtain the higher order functions by upward recursion  .
!     ...........................................................
        do n = 3,nc
        rn = real(n-2)
        by(n) = (2.0*rn+1.0)*by(n-1)/x-by(n-2)
        end do
!     ................................................
!     .  bj(*) calculation - set the starting order  .
!     .                      for downward recursion  .
!     ................................................
      nst = nc+int((101.0+x)**.5)                            !         eq 4.21
!     ....................................................
!     .  the t(*) array is used to recur down to the     .
!     .    two highest order functions that are needed   .
!     .  set starting values for the two highest orders  .
!     .    nst and nst-1                                 .
!     ....................................................
      t(3) = 0.0
      t(2) = 1.0e-35
!     ...................................................
!     .  recur downward to obtain orders nc-1 and nc-2  .
!     ...................................................
        do i = nst-1,nc-1,-1
        ri = real(i)
        t(1) = (2.0*ri+1.0)*t(2)/x-t(3)
        t(3) = t(2)
        t(2) = t(1)
        end do
!     ...............................................
!     .  continue downward recursion to order zero  .
!     ...............................................
      bj(nc) = t(3)
      bj(nc-1) = t(2)
      do i = nc-2,1,-1
        ri = real(i)
        bj(i) = (2.0*ri+1.0)*bj(i+1)/x-bj(i+2)
      end do                
!     ..................................................
!     .  calculate the scale factor and the functions  .
!     ..................................................
      alpha = a/bj(1)
      do k = 1,nc
        hankel(k) = cmplx(bj(k)*alpha,by(k))
      end do
      end subroutine beshOld

      subroutine genlgpOld(theta,pnmllg,nc)
!     ........................................................
!     .  calculate associated Legendre functions (argument   .
!     .    cos(theta)) divided by sin(theta) for m = 1       .
!     .  generate first two orders by formula and remaining  .
!     .    orders by recursion                               .
!     .                                                      .
!     .  pnmllg = associated Legendre function/sin(theta)    .
!     .  nc = number of orders (0 to nc-1)                   .
!     .  the order of the associated Legendre functions is   .
!     .    incremented by one in the pnmllg(*) array         .
!     ........................................................
      real    :: theta, costh, rn
      integer :: nc, n
      real pnmllg(:)
      costh = cos(theta)
!     ..............................
!     .  calculate orders 0 and 1  .
!     ..............................
      pnmllg(1) = 0.0                                                
      pnmllg(2) = 1.0                                                
!     .................................................
!     .  recur upward to obtain all remaining orders  .
!     .................................................
      do n = 3,nc 
      rn = real(n-1)
      pnmllg(n) = ((2.0*rn-1.0)*costh*pnmllg(n-1) &
                  -rn*pnmllg(n-2))/(rn-1.0)     
      enddo
      end subroutine genlgpOld

end module mieDistPhaseMatrix_mod
