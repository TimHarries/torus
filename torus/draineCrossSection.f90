subroutine draineCrossSection(lambda, nLambda, aMin, aMax, qDist, &
     kappaAbs, kappaSca, grainType)
  use constants_mod
  use utils_mod
  use unix_mod
  implicit none
  character(len=80) :: dataDirectory, filename
  real :: lambda(*)
  real :: aMin, aMax, qDist, normFac
  character(len=*) :: grainType
  integer, parameter :: nDist = 100
  integer :: nLambda
  real :: t
  real,allocatable  :: kappaExt(:)
  real :: kappaSca(*)
  real :: kappaAbs(*)
  integer, parameter :: maxpts = 300
  integer :: npts
  real :: xArray1(maxpts), qAbs1(maxpts), omega1(maxpts), qSca1(maxpts)
  real :: xArray2(maxpts), qAbs2(maxpts), omega2(maxpts)
  real :: nindex(maxpts), kindex(maxpts)
  complex :: m(maxpts), thisM
  real :: junk1, junk2
  integer :: errNo
  integer :: i, j, k
  real :: logamin, logamax, loga1, loga2, a1, a2, da, a, tot
  real :: q1, om1, q2, om2, q, om

  real :: x, qSca, qExt, qAbs, fac



  allocate(kappaExt(1:nLambda))
  kappaSca(1:nLambda) = 0.
  kappaAbs(1:nLambda) = 0.
  kappaExt(1:nLambda) = 0.

  call unixGetenv("TORUS_DATA",dataDirectory,i,errno)
  dataDirectory = trim(dataDirectory)//"/"
  select case(grainType)


  case("draine")

     npts = 241
     filename = trim(dataDirectory)//"draine-graphite-0.01.dat"
     open(20,file=filename, status="old", form="formatted")
     do i = npts, 1, -1
        read(20,*) xArray1(i), qAbs1(i), qSca1(i)
     enddo
     close(20)

     do j = 1, nLambda

        if ((lambda(j)*angsToMicrons < xArray1(1)).or.(lambda(j)*angsToMicrons > xArray1(npts))) then
           kappaAbs(j) = 0.
           kappaSca(j) = 0.
           kappaExt(j) = 0.
        else

           call locate(xArray1, npts, lambda(j)*angsToMicrons, k)
           t   = (lambda(j)*angsToMicrons - xArray1(k))/(xArray1(k+1)-xArray1(k))
           kappaAbs(j)= (qAbs1(k)+ t * (qAbs1(k+1)-qAbs1(k))) *  pi * (0.01 * microntocm)**2
           kappaSca(j)= (qSca1(k)+ t * (qSca1(k+1)-qSca1(k))) *  pi * (0.01 * microntocm)**2
           kappaExt(j) = kappaAbs(j) + kappaSca(j)
           write(*,*) kappaAbs(j), kappaExt(j), kappaSca(j)
        endif

     enddo




  case("graphite")

     npts = 53
     filename = trim(dataDirectory)//"graphite0.01.dat"
     open(20,file=filename, status="old", form="formatted")
     do i = 1, npts
        read(20,*) xArray1(i), qAbs1(i), omega1(i)
     enddo
     close(20)
     filename = trim(dataDirectory)//"graphite0.1.dat"
     open(20,file=filename, status="old", form="formatted")
     do i = 1, npts
        read(20,*) xArray2(i), qAbs2(i), omega2(i)
     enddo
     close(20)

     logamin = log(aMin)
     logamax = log(aMax)


     do j = 1, nLambda

        if ((lambda(j)*angsToMicrons < xArray1(1)).or.(lambda(j)*angsToMicrons > xArray1(npts))) then
           kappaAbs(j) = 0.
           kappaSca(j) = 0.
           kappaExt(j) = 0.
        else

           call locate(xArray1, npts, lambda(j)*angsToMicrons, k)
           t   = (lambda(j)*angsToMicrons - xArray1(k))/(xArray1(k+1)-xArray1(k))
           q1 = qAbs1(k)+ t * (qAbs1(k+1)-qAbs1(k))
           om1 = omega1(k) + t * (omega1(k+1)-omega1(k))

           call locate(xArray2, npts, lambda(j)*angsToMicrons, k)
           t   = (lambda(j)*angsToMicrons - xArray2(k))/(xArray2(k+1)-xArray2(k))
           q2 = qAbs2(k)+ t * (qAbs2(k+1)-qAbs2(k))
           om2 = omega2(k) + t * (omega2(k+1)-omega2(k))

           tot = 0.
           do i = 1, nDist-1
              loga1 = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
              loga2 = logAmin+(logAmax-logAmin)*real(i)/real(nDist-1)
              a1 = exp(loga1)
              a2 = exp(loga2)
              tot = tot + 0.5*(a1**(qDist) + a2**(qDist))*(a2-a1)
           enddo
           normFac = 1./tot


           kappaAbs(j) = 0.
           kappaSca(j) = 0.
           kappaExt(j) = 0.
           do i = 1, nDist-1
              loga1 = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
              loga2 = logAmin+(logAmax-logAmin)*real(i)/real(nDist-1)
              a1 = exp(loga1)
              a2 = exp(loga2)

              da = a2 - a1
              a = 0.5*(a1+a2)

              if ((a > 0.01e-4).and.(a < 0.1e-4)) then
                 q = q1 + q2 *(a - 0.01e-4)/(0.1e-4 - 0.01e-4)
                 om = om1 + om2 *(a - 0.01e-4)/(0.1e-4 - 0.01e-4)
              endif
              if (a < 0.01e-4) then
                 q = q1
                 om = om1
              endif
              if (a > 0.1e-4) then
                 q = q2
                 om = om2
              endif

              kappaAbs(j) = kappaAbs(j) + normFac*a**(qDist)*da * q * pi * a**2
              kappaExt(j) = kappaExt(j) + normFac*a**(qDist)*da * q/(1.-om) * pi * a**2
              kappaSca(j) = kappaSca(j) + normFac*a**(qDist)*da * om*q/(1.-om) * pi * a**2
           enddo

           write(*,*) kappaAbs(j), kappaExt(j), kappaSca(j)
        endif
     enddo

  case("silicate")


     npts = 30
     filename = trim(dataDirectory)//"silicate.dat"
     open(20,file=filename, status="old", form="formatted")
     do i = 1, npts
        read(20,*) xArray1(i), junk1, junk2, qAbs1(i), omega1(i),qAbs2(i), omega2(i)
     enddo
     close(20)


     xArray2 = xArray1

     logamin = log(aMin)
     logamax = log(aMax)


     do j = 1, nLambda

        if ((lambda(j)*angsToMicrons < xArray1(1)).or.(lambda(j)*angsToMicrons > xArray1(npts))) then
           kappaAbs(j) = 0.
           kappaSca(j) = 0.
           kappaExt(j) = 0.
        else
           call locate(xArray1, npts, lambda(j)*angsToMicrons, k)
           t   = (lambda(j)*angsToMicrons - xArray1(k))/(xArray1(k+1)-xArray1(k))
           q1 = qAbs1(k)+ t * (qAbs1(k+1)-qAbs1(k))
           om1 = omega1(k) + t * (omega1(k+1)-omega1(k))

           call locate(xArray2, npts, lambda(j)*angsToMicrons, k)
           t   = (lambda(j)*angsToMicrons - xArray2(k))/(xArray2(k+1)-xArray2(k))
           q2 = qAbs2(k)+ t * (qAbs2(k+1)-qAbs2(k))
           om2 = omega2(k) + t * (omega2(k+1)-omega2(k))
           tot = 0.

           do i = 1, nDist-1
              loga1 = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
              loga2 = logAmin+(logAmax-logAmin)*real(i)/real(nDist-1)
              a1 = exp(loga1)
              a2 = exp(loga2)
              tot = tot + 0.5*(a1**(qDist) + a2**(qDist)) *(a2-a1)
           enddo
           normFac = 1./tot


           kappaAbs(j) = 0.
           kappaSca(j) = 0.
           kappaExt(j) = 0.
           do i = 1, nDist-1
              loga1 = logAmin+(logAmax-logAmin)*real(i-1)/real(nDist-1)
              loga2 = logAmin+(logAmax-logAmin)*real(i)/real(nDist-1)
              a1 = exp(loga1)
              a2 = exp(loga2)

              da = a2 - a1
              a = 0.5*(a1+a2)

              q = q2
              om = om2

              if (a < 0.05*1e-4) then
                 q = q1
                 om = om1
              endif
              kappaAbs(j) = kappaAbs(j) + normFac*a**(qDist)*da * q * pi * a**2
              kappaExt(j) = kappaExt(j) + normFac*a**(qDist)*da * q/(1.-om)* pi * a**2
              kappaSca(j) = kappaSca(j) + normFac*a**(qDist)*da * om*q/(1.-om)* pi * a**2

           enddo
        endif
     enddo

  case("amorphous")

     ! amorphous carbon refractive indices from W.W.Duley, 1984, ApJ, 694

     filename=trim(dataDirectory)//"amorphous.dat"
     open(20,file=filename,status="old",form="formatted")
     npts = 14
     do i = 1, npts
        read(20,*) xArray1(i), nindex(i), kindex(i)
        xArray1(i) = 1./xArray1(i)
        m(i) = cmplx(nindex(i), -kindex(i))
     enddo
     close(20)
     do i = 1, nLambda

        if ((lambda(j)*angsToMicrons < xArray1(1)).or.(lambda(j)*angsToMicrons > xArray1(npts))) then
           kappaAbs(j) = 0.
           kappaSca(j) = 0.
           kappaExt(j) = 0.
        else
           call locate(xArray1, npts, lambda(i)*angsToMicrons, k)
           t   = (lambda(i)*angsToMicrons - xArray1(k))/(xArray1(k+1)-xArray1(k))
           thisM = m(k) + t * (m(k+1)-m(k))

           tot = 0.
           do j = 1, 1000
              a1 = 10.**(log10(amin)+(log10(amax)-log10(amin))*real(j-1)/1000.)
              a2 = 10.**(log10(amin)+(log10(amax)-log10(amin))*real(j)/1000.)
              a = 0.5*(a1 + a2)
              x = (2. * pi * a)/(lambda(i)*angstromTocm)
              qSca = (8./3.)*(x**4)* abs( (thism**2 - 1)/(thism**2 + 2))**2
              qAbs = -4. * x * aimag( (thism**2 - 1)/(thism**2 + 2) )
              qExt = qSca + qAbs
              fac = a**(qDist)*(a2-a1)
              kappaSca(i) = kappaSca(i) + qSca * pi * a**2 * fac
              kappaAbs(i) = kappaAbs(i) + qAbs * pi * a**2 * fac
              kappaExt(i) = kappaExt(i) + qExt * pi * a**2 * fac
              tot = tot + fac
           enddo
           kappaSca(i) = kappaSca(i) / tot
           kappaExt(i) = kappaExt(i) / tot
           kappaAbs(i) = kappaAbs(i) / tot
           write(*,*) kappaAbs(i), kappaExt(i), kappaSca(i)
        endif

     enddo



  case DEFAULT

     write(*,'(a,a)') "! Grain type not recognised: ",grainType
     stop

  end select



end subroutine draineCrossSection

  
  
