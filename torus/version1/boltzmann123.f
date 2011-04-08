       program boltzmann123

       implicit none
       
       integer npts, ncut, n, i
       parameter (npts = 1001)

       double precision k, eV, Z, part
       double precision g1,g2,g3,E1,E2,E3,T, Tincr, Tmin, Tmax
       double precision n21(npts),n31(npts),n32(npts)
       double precision n1rat(npts), n2rat(npts), n3rat(npts)


       print *
       write(6,*)'compute H-like rel.level pops (1,2,3) over range in T'
       print *

       write(6,*)'enter number of protons in H-like specie; Z=1 for H'
       read(5,*) Z

* some constants

       k = 1.380658d-16     ! Boltzmann's constant
       eV = 1.60217733d-12  ! ergs/eV
       Z = 1.d+00
       E1 = -13.59844d+00 * Z * Z
       E2 = E1 / (2.d+00*2.d+00)
       E3 = E1 / (3.d+00*3.d+00)
       g1 = 2.d+00
       g2 = 8.d+00
       g3 = 18.d+00

       open(16,file='bolt123.out',status='unknown')


       print *
       write(6,*)'also compute N_n/Ntot n = 1,2,3, with less accuracy'
       write(6,*)'because H-like partition function blows up at large n'
       print *

       write(6,*)'enter max n to consider for H-like partition function'
       read(5,*) ncut

       print *
       write(6,*)'enter min, max T(K)'
       read(5,*) Tmin, Tmax

       write(16,96) Z
96     format('# relative level populations for H-like Z = ', f3.0)
       write(16,97)
97     format('# current form of H-like part.func uncertain at high T')
       write(16,98) ncut
98     format('# max. n level considered in partition function: ', i3)
       write(16,99)
99     format('# T(K)   N21    N31    N32   N1rat    N2rat   N3rat   Z')

       Tincr = (Tmax - Tmin)/(npts - 1)
       T = Tmin - Tincr

       do i = 1, npts

       T = T + Tincr

       n21(i) = (g2/g1) * dexp( -(E2 - E1)*eV/(k*T) )
       n31(i) = (g3/g1) * dexp( -(E3 - E1)*eV/(k*T) )
       n32(i) = n31(i)/n21(i)

! compute the partition function, cutoff at n < 100
       part = 0.
       do n = 1, ncut
        part = part + 2.d+00*n*n * 
     -  dexp( E1*eV*( 1.d+00 - 1.d+00/(n*n) )/(k*T) )
       enddo

       n1rat(i) = (g1/part) * 
     - dexp( E1*eV*( 1.d+00 - 1.d+00/(1.d+00*1.d+00) )/(k*T) )
       n2rat(i) = (g2/part) * 
     - dexp( E1*eV*( 1.d+00 - 1.d+00/(2.d+00*2.d+00) )/(k*T) )
       n3rat(i) = (g3/part) * 
     - dexp( E1*eV*( 1.d+00 - 1.d+00/(3.d+00*3.d+00) )/(k*T) )

       write(16,100)T,n21(i),n31(i),n32(i),n1rat(i),n2rat(i),n3rat(i),
     -              part
100    format(1p, 7(e9.3, 2x), e10.4 )

       enddo
       print *
       close(16)

       write(6,*)' output sent to bolt123.out'
       print *

       stop
       end
