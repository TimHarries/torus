
      subroutine mieDistPhaseMatrix(aMin, aMax, a0, qDist, pDist, lambda, &
                                costheta, mieMatrix, cmr, cmi)

      use constants_mod
      use phasematrix_mod

      implicit none


      type(PHASEMATRIX) :: mieMatrix
      real :: normFac, dMu
      real :: cosTheta
      integer :: nMu = 200
      real :: aMin, aMax, a0, qDist, pDist, lambda, aFac
      integer :: nDist
      real ::  x, cmr, cmi,  t, theta, costh, rn, p1, p2
      real :: logamax, logamin
      integer :: nc,nci,n, i, n1,j

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
      ci = (0.0,1.0)

      ! for amorphous carbon grains we assume we're in the rayleigh regime

      cm = cmplx(cmr,cmi)


      micronsTocm = 1.e-4

      p11tot = 0.
      pltot = 0.
      p33p11tot = 0.
      p34p11tot = 0.

      nDist = 50
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
        x = 2.*pi*(a * micronsToCm)/(lambda*1.e-8)
        x = max(1.e-5, x)
        gfac = pi*(a * micronsToCm)**2

!     .........................................
!     .  set the complex index of refraction  .
!     .    for an exp(-iwt) time variation    .
!     .........................................
      call sphere(x,cm,nc)
      nci = nc+1
      t = 0.0
!     .............................................................
!     .  calculate t = Qsca*x**2                                  .
!     .              = (scattering cross section)*x**2/(pi*a**2)  .
!     .............................................................
      do  n = 1,nc
        t = t+(abs(f(n))**2+abs(g(n))**2)*cnrm(n)
      enddo
      t = gfac*t/x**2
      costh = costheta
      theta = acos(costheta)
      call genlgp2(theta,pnmllg,nci)
      s1 = 0.0
      s2 = 0.0
      do n = 1,nc
          n1 = n+1
          cim = ci**(-n1)
          rn = real(n)
          p1 = rn*costh*pnmllg(n1)-(rn+1.0)*pnmllg(n)
          p2 = pnmllg(n1)
!     ........................................
!     .  calculate parallel field amplitude  .
!     .    for parallel incident             .     !                  eq 3.50b
!     ........................................                          and
          s2 = s2+cim*(-ci*p2*f(n)+p1*g(n))*cnrm(n) !                  eq 4.10a
!     .............................................
!     .  calculate perpendicular field amplitude  .
!     .    for perpendicular incident             .
!     .  use coefficients for parallel incident   .
!     .    f(e1n) = -f(o1n),  g(o1n) = g(ein)     .  !                 eq 3.50a
!     .............................................                     and
          s1 = s1+cim*(ci*p1*(-f(n))+p2*g(n))*cnrm(n) !                eq 4.11b
        enddo
!     ....................................................
!     .  calculate P11, PL = -P12/P11, P33/P11, P34/P11  .              eq 4.25
!     ....................................................              and
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
!      call sphere(x,cm,nc)
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
!      call genlgp2(theta,pnmllg,nci)
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
!       enddo
!
!       normFac = normFac + p11tot*dmu
!    enddo
!      
!
!      normFac = normFac * 0.5
!
!      do i = 1, 4
!       do j = 1 , 4
!         mieMatrix%element(i,j) = MieMatrix%element(i,j)/normFac
!       enddo
!      enddo
!
     end subroutine mieDistPhaseMatrix


      subroutine sphere(x,cm,nc)
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
      do 10 n = nmx,nc+1,-1
        rn = real(n)
        an = rn/z-1.0/(an+rn/z)
10    continue
      amat(nc) = an
      do 20 n = nc,2,-1
        rn = real(n)
        amat(n-1) = rn/z-1.0/(amat(n)+rn/z)              !             eq 4.19
20    continue
!     ...................................................
!     .  calculate the Bessel functions - the order is  .
!     .    incremented by one in the hkl(*) array       .
!     ...................................................
      call besh(x,hkl,nci)
      bj = real(hkl(1))
!     ................................
!     .  calculate the coefficients  .
!     ................................
      do 30 n = 1,nc
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
30    continue
      return
      end
      subroutine besh(x,hankel,nc)
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
        do 10 n = 3,nc
        rn = real(n-2)
        by(n) = (2.0*rn+1.0)*by(n-1)/x-by(n-2)
10      continue
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
        do 20 i = nst-1,nc-1,-1
        ri = real(i)
        t(1) = (2.0*ri+1.0)*t(2)/x-t(3)
        t(3) = t(2)
        t(2) = t(1)
20      continue
!     ...............................................
!     .  continue downward recursion to order zero  .
!     ...............................................
      bj(nc) = t(3)
      bj(nc-1) = t(2)
        do 30 i = nc-2,1,-1
        ri = real(i)
        bj(i) = (2.0*ri+1.0)*bj(i+1)/x-bj(i+2)
30      continue                
!     ..................................................
!     .  calculate the scale factor and the functions  .
!     ..................................................
      alpha = a/bj(1)
        do 40 k = 1,nc
        hankel(k) = cmplx(bj(k)*alpha,by(k))
40      continue
      return
      end

      subroutine genlgp(theta,pnmllg,nc)
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
        real    :: theta, pnmllg, costh, rn
        integer :: nc, n
        dimension pnmllg(nc)

      costh = cos(theta)
!     ..............................
!     .  calculate orders 0 and 1  .
!     ..............................
      pnmllg(1) = 0.0                                          !       eq 4.70a
      pnmllg(2) = 1.0                                            !     eq 4.70b
!     .................................................
!     .  recur upward to obtain all remaining orders  .
!     .................................................
      do 10 n = 3,nc 
      rn = real(n-1)
      pnmllg(n) = ((2.0*rn-1.0)*costh*pnmllg(n-1) &
            -rn*pnmllg(n-2))/(rn-1.0)                      !     eq 4.71
10    continue 
      return 
      end 


      subroutine genlgp2(theta,pnmllg,nc)
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
      real pnmllg(nc)
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
      end 
