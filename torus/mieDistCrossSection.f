      subroutine mieDistCrossSection(aMin, aMax,
     & qDist, lambda, cmr, cmi,
     & kappaExt, kappaSca, kappaAbs)
      implicit none
      real  cmr, cmi, x, kappaExt, kappaAbs, kappaSca
      real aMin, aMax, qDist, lambda
      real logAmin, logAmax
      real tot
      integer nDist
      real pi, rtd
      real qsca(100), qext(100), qback(100), gsca
      integer ic, ip, id, npnts,i,n,n1, nc,nci,nang
      real x1,xn, dk, x1i,xni,dwl,rn,rf,ri,xi,theta,costh
      real snorm, p1,p2,tp,dlt
      real micronsToCm
      complex cm,ci,cim,f(1000),g(1000),t
      common /cfcom2/ f,g,cnrm
      real pnmllg(1001),cnrm(1000)
      real normFac, a, da
      real gfac
      real loga1, loga2, a1, a2
      real nsd(100), aDist(100)
      complex refrel
      complex s1(1000),s2(1000)
      
      pi = 3.14159265358979
      rtd = 180.0/pi
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
          nsd(i) = aDist(i)**(-qDist)
        enddo
        call PowerInt(nDist,1,nDist,aDist,nsd,normFac)
        
        do i = 1, nDist
          a = aDist(i) 
          x = 2.*pi*(a * micronsTocm)/(lambda*1.e-8)
          x = max(x, 1.e-5)
          refrel = cmplx(cmr, cmi)
          call BHMIE(X,REFREL,2 ,S1,S2,QEXT(i),QSCA(i),QBACK(i),GSCA)
          
          qExt(i) = qExt(i) * nsd(i) * pi * (a * micronsToCm)**2
          qSca(i) = qSca(i) * nsd(i) * pi * (a * micronsToCm)**2
        enddo
        
         call powerInt(nDist, 1, nDist, aDist, qExt, kappaExt)
         kappaExt = kappaExt / normFac
         call powerInt(nDist, 1, nDist, aDist, qSca, kappaSca)
         kappaSca = kappaSca / normFac
        
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
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER i, N, N1, N2                                                      
      REAL x(N), y(N), integral, pow, C, delint                     
c -----------------------------------------------------------------------       
c     set integral to 0 and accumulate result in the loop                       
      integral = 0.0                                                            
c     calculate weight, wgth, and integrate in the same loop                    
      IF (N2.GT.N1) THEN                                                        
        DO i = N1, N2-1                                                         
          pow = log(y(i+1)/y(i)) / log(x(i+1)/x(i))                           
          C = y(i) / x(i)**pow                                                  
          delint = (x(i+1)**(pow+1)-x(i)**(pow+1.))*C/(pow+1.)                  
c         add contribution to the integral                                      
          integral = integral + delint                                          
        END DO                                                                  
        ELSE                                                                    
c        integral = 0.0                                                         
        integral = y(1)                                                         
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
        
