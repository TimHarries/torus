      subroutine mieDistCrossSection(aMin, aMax, a0, 
     & qDist, pDist, lambda, cmr, cmi,
     & kappaExt, kappaSca, kappaAbs, gscadist)

      use bhmie_mod, only: MXNANG

      implicit none
      real  cmr, cmi, x, kappaExt, kappaAbs, kappaSca
      real aMin, aMax, a0, qDist, pDist, lambda
      real logAmin, logAmax
      integer nDist
      real pi
      real qsca(100), qext(100), qback(100), gsca(100)
      integer ic, ip, i
      real micronsToCm
      complex cm,ci,f(1000),g(1000)
      common /cfcom2/ f,g,cnrm
      real cnrm(1000)
      real normFac, a
      real gscadist
      real nsd(100), aDist(100)
      complex refrel
      complex s1(2*MXNANG-1),s2(2*MXNANG-1)
      integer :: nang=2
      
      pi = 3.14159265358979
      ci = (0.0,1.0)
      ic = 1
      ip = 1

      micronsToCm = 1.e-4
      
c     .........................................
c     .  set the complex index of refraction  .
c     .    for an exp(-iwt) time variation    .
c     .........................................
      cm = cmplx(cmr,cmi)
        
        nDist = 100
        logamin = log(aMin)
        logamax = log(aMax)

        do i = 1, nDist

          aDist(i) = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
          aDist(i) = exp(aDist(i))
          nsd(i) = aDist(i)**(-qDist) * exp(-(aDist(i)/a0)**pDist)
c     For safty (R. Kurosawa fixed this.)
          if (nsd(i) .lt. 1.0e-30) nsd(i) = 1.0e-30

        enddo


        call PowerInt(nDist,1,nDist,aDist,nsd,normFac)
!        open(unit=54, file='g_scat.dat', status='unknown')

!        refrel = cmplx(0.d0, 0.d0)
!        s1 = cmplx(0.d0, 0.d0)
!        s2 = cmplx(0.d0, 0.d0)
!        QEXT = 0.d0
!        QSCA = 0.d0
!        QBACK =0.d0
!        GSCA = 0.d0

        do i = 1, nDist

          a = aDist(i) 
          x = 2.*pi*(a * micronsTocm)/(lambda*1.e-8)
          x = max(x, 1.e-5)
          refrel = cmplx(cmr, cmi)

          call BHMIE(X,REFREL,NANG ,S1,S2,QEXT(i),QSCA(i),QBACK(i),
     &       GSCA(i))

          qExt(i) = qExt(i) * nsd(i) * pi * (a * micronsToCm)**2
          qSca(i) = qSca(i) * nsd(i) * pi * (a * micronsToCm)**2
!          write(54, *) x, gsca(i)
        enddo

!        close(54)

         call powerInt(nDist, 1, nDist, aDist, qExt, kappaExt)
         kappaExt = kappaExt / normFac
         call powerInt(nDist, 1, nDist, aDist, qSca, kappaSca)
         kappaSca = kappaSca / normFac
         gScaDist = gSca(ndist)
c     ..........................................
c     .   calculate the absorption efficiency  .
c     ..........................................

        kappaabs = kappaext-kappasca                     
      END

      SUBROUTINE PowerInt(N,N1,N2,x,y,integral)
c =======================================================================       
c This subroutine calculates integral I(y(x)*dx). Both y and x are              
c 1D arrays, y(i), x(i) with i=1,N (declared with NN). Lower and upper          
c integration limits are x(N1) and x(N2), respectively. The method used         
c is a power-law approximation for y(x) between any two points .                
c                                                      [Z.I., Mar. 1996]        
c -----------------------------------------------------------------------
c Modefied on 24-apr-2003: Now computations are done internally 
c                          using double precision.  The inputs and 
c                          the outputs remains single precision.
c                                    --- (R. Kurosawa)
c                          
c
c =======================================================================       
      IMPLICIT none
      INTEGER i, N, N1, N2
      REAL x(N), y(N), integral
!      REAL x(N), y(N), integral, pow, C, delint                     
      DOUBLE PRECISION x1, x2, y1, y2, pow, delint, C, sum
c -----------------------------------------------------------------------       
c     set integral to 0 and accumulate result in the loop                       
      sum = 0.0d0 
      integral=0.0
c     calculate weight, wgth, and integrate in the same loop                    
      IF (N2.GT.N1) THEN 
        DO i = N1, N2-1 
           x1 = dble(x(i))
           x2 = dble(x(i+1))
           y1 = dble(y(i))
           y2 = dble(y(i+1))
c
          pow = log(y2/y1) / log(x2/x1)                           
          C = y1 / x1**pow                      
          delint = (x2**(pow+1.0d0)-x1**(pow+1.0d0))*C/(pow+1.d0)
c         add contribution to the integral                                      
          sum = sum + delint                                          
        END DO                
        ELSE              
c        integral = 0.0                                                         
        sum = y(1)          
      END IF                   
c -----------------------------------------------------------------------       
c Hope this is ok....
      integral = real(sum)

      RETURN                                
      END                           
c ***********************************************************************       
        
