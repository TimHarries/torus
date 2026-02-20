module dust_mod

  use constants_mod
  use citations_mod
  use messages_mod
  use vector_mod
  use gridtype_mod, only: GRIDTYPE
  use utils_mod, only: locate
  use octal_mod, only: OCTAL, subcellCentre, cellVolume
  use amr_mod, only: amrGridValues, returnKappa, octalOnThread, findTotalMass
  use mpi_global_mod
  use mieDistCrossSection_mod, only: mieDistCrossSection

  implicit none
  public
  private :: fillGridMie, fillAMRgridMie, dustPropertiesfromFile, &
       parseGrainType, &
       getMeanMass2 !, getTemperatureDensityRundust, returnScaleHeight,  rtnewtdust, Equation2dust

contains

   

  subroutine effectiveMedium(epsilonEffective, epsilonMatrix, epsilonInclusion, fillingFactor)
    real :: fillingFactor
    complex :: epsilonEffective, epsilonMatrix, epsilonInclusion, epsilonEffectiveMG

! maxwell garnett formula (see Maron & Maron, 2005, MNRAS, 357, 873)

    epsilonEffectiveMG = epsilonMatrix + 3.0 * fillingFactor * epsilonMatrix * &
         (epsilonInclusion - epsilonMatrix)/(epsilonInclusion + 2.0*epsilonMatrix - &
         fillingFactor * (epsilonInclusion - epsilonMatrix))

!    write(*,*) "eps maxwell ",epsilonEffective
! bruggeman formula (same paper)

    call searchBisection(epsilonEffective,epsilonMatrix, epsilonInclusion, fillingFactor, &
         0.5*epsilonEffectiveMG, 2.*epsilonEffectiveMG)
!    write(*,*) "eps bruggeman ",epsilonEffective

  end subroutine effectiveMedium

  subroutine searchBisection(eps, epsMatrix, epsInclusion, f, eps1, eps2)
    complex :: eps, epsMatrix, epsInclusion,eps1,eps2
    real :: f
    complex :: a, b, c
    real :: fa, fb, fc
    logical :: converged
    integer :: i
    a = eps1
    b = eps2
    fa = bruggeman(a, epsMatrix, epsInclusion, f)
    fb = bruggeman(b, epsMatrix, epsInclusion, f)
    if (fa*fb > 0.d0) then
       if (writeoutput) then
          write(*,*) i,a,b,fa,fb," bisection failed"
          stop
       endif
    endif
    i = 0
    converged = .false.
    do while (.not.converged)
       i = i + 1
       fa = bruggeman(a, epsMatrix, epsInclusion, f)
       fb = bruggeman(b, epsMatrix, epsInclusion, f)
!       if (writeoutput) write(*,*) i,a,b,fa,fb
       c = 0.5 * (a + b)
       fc = bruggeman(c, epsMatrix, epsInclusion, f)
       if (fa*fc < 0.d0) then
          b = c
       else
          a = c
       endif
       if (abs(a-b)/abs(a+b) < 1.d-6) then
          converged = .true.
       endif
    end do
    eps = 0.5*(a+b)
  end subroutine searchBisection


    real function bruggeman(eps, epsMatrix, epsInclusion, f)
      complex :: eps, epsMatrix, epsInclusion
      real :: f

      bruggeman = real(f * (epsInclusion - eps)/(epsInclusion+2.*eps) + &
           (1.-f)*(epsMatrix - eps)/(epsMatrix + 2.*eps))
    end function bruggeman

  subroutine refractiveIndexToPermittivity(n, k, epsilonReal, epsilonImg)
    real :: n, k, epsilonReal, epsilonImg
    epsilonReal = n**2 - k**2
    epsilonImg = 2.0 * n * k
  end subroutine refractiveIndexToPermittivity

  subroutine permittivityToRefractiveIndex(epsilonReal, epsilonImg, n, k)
    real :: n, k, epsilonReal, epsilonImg

    n = sqrt( (sqrt(epsilonreal**2 + epsilonImg**2) + epsilonReal)/2.0 )
    k = sqrt( (sqrt(epsilonreal**2 + epsilonImg**2) - epsilonReal)/2.0 )
  end subroutine permittivityToRefractiveIndex

  subroutine dumpPolarizability(miePhase, nMuMie, lambda, nLambda)
    use inputs_mod, only : polarWavelength, polarFilename
    use phasematrix_mod, only : phasematrix
    type(PHASEMATRIX) :: miePhase(:,:,:)
    integer :: nMuMie
    real :: lambda(:)
    integer :: nLambda
    real :: ang, mu
    integer :: i, j, k


    if (writeoutput) then
    do i = 1, 1
       call locate(lambda, nlambda, real(polarWavelength), k)
       open(23, file=polarFileName, status="unknown", form="formatted")
       write(23,'(a9,a9)') "# angle","-S21/S11"
       do j = nMuMie, 1, -1
          mu = 2.*real(j-1)/real(nMumie-1)-1.
          ang = acos(mu) * real(radtoDeg)
          write(23,'(f9.2,4f9.3)') ang, -miePhase(1,k,j)%element(1,2)/miePhase(1,k,j)%element(1,1)
!               -miePhase(2,k,j)%element(1,2)/miePhase(2,k,j)%element(1,1), &
!               -miePhase(3,k,j)%element(1,2)/miePhase(3,k,j)%element(1,1), &
!               -miePhase(4,k,j)%element(1,2)/miePhase(4,k,j)%element(1,1)
       enddo
       close(23)
    enddo
    endif
  end subroutine dumpPolarizability

       

  subroutine getRefractiveIndex(lambda, nLambda, graintype, mReal, mImg, porousFillingFactor)
    
    use unix_mod, only: unixGetenv
    real :: lambda(:)
    real :: porousFillingFactor
    integer :: nLambda, nRef
    real, allocatable :: tempIm(:), tempReal(:), lamRef(:)
    real(double) :: dydx
    integer :: i, n, j
    character(len=*), intent(in) :: graintype
    character(len=200) :: filename, dataDirectory
    character(len=100) :: textline, cjunk
    real :: junk1, junk2
    real :: mReal(:), mImg(:), t,  x(5000), y1(5000,3), y2(5000,3)
    complex :: epsilonEffective, epsilonMatrix, epsilonInclusion
    real :: epsilonImg, epsilonReal
    logical :: firstTime = .true.
    dataDirectory = " "
    select case(graintype)

    case("amc_zb")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"amC-zb2.nk"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       nRef = 1245
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       do i = 1, nRef
          read(20,*) lamRef(i), tempReal(i), tempIm(i)
       enddo
       close(20)

    case("draine_sil")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"callindex.out_silD03.txt"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       read(20,*) cjunk
       read(20,*) cjunk
       read(20,*) cjunk
       read(20,*) cjunk
       read(20,*) cjunk
       nRef = 837
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       do i = nRef, 1, -1
          read(20,*) lamRef(i), junk1, junk2, tempReal(i), tempIm(i)
          tempReal(i) = tempReal(i) + 1.
       enddo
       close(20)


    case("sil_dl")
       call addBibcode("1984ApJ...285...89D","Silicite grain optical properties")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"eps_Sil.txt"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       read(20,*) cjunk
       read(20,*) cjunk
       read(20,*) cjunk
       read(20,*) cjunk
       read(20,*) cjunk
       read(20,*) cjunk
       nRef = 1201
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       do i = nRef, 1, -1
          read(20,*) lamRef(i), junk1, junk2, tempReal(i), tempIm(i)
          tempReal(i) = tempReal(i) + 1.
       enddo
       close(20)


    case("am_olivine")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"am_olivine.txt"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       n = 0
31     continue
       read(20,'(a)',end=32) textline
       if (textLine(1:1) == "#") goto 31
       n = n + 1
       read(textline,*) x(n),y1(n,1), y2(n,1)
       goto 31
32     continue
       nRef = n
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = y1(1:n,1)
       tempIm(1:nRef) = y2(1:n,1)
       close(20)


    case("am_pyroxene")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"am_pyroxene.txt"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       n = 0
41     continue
       read(20,'(a)',end=42) textline
       if (textLine(1:1) == "#") goto 41
       n = n + 1
       read(textline,*) x(n),y1(n,1), y2(n,1)
       goto 41
42     continue
       nRef = n
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = y1(1:n,1)
       tempIm(1:nRef) = y2(1:n,1)
       close(20)

    case("forsterite")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       do j = 1, 3
          write(filename,'(a,i1.1,a)') trim(dataDirectory)//"/"//"forsterite",j,".txt"
          if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
          open(20,file=filename,status="old",form="formatted")
          n = 0
51        continue
          read(20,'(a)',end=52) textline
          if (textLine(1:1) == "#") goto 51
          n = n + 1
          read(textline,*) x(n),y1(n,j), y2(n,j)
          goto 51
52        continue
          nRef = n
          close(20)
       enddo
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = (y1(1:n,1)+y1(1:n,2)+y1(1:n,3))/3.0
       tempIm(1:nRef) =  (y2(1:n,1)+y2(1:n,2)+y2(1:n,3))/3.0

    case("enstatite")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       do j = 1, 3
          write(filename,'(a,i1.1,a)') trim(dataDirectory)//"/"//"enstatite",j,".txt"
          if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
          open(20,file=filename,status="old",form="formatted")
          n = 0
61        continue
          read(20,'(a)',end=62) textline
          if (textLine(1:1) == "#") goto 61
          n = n + 1
          read(textline,*) x(n),y1(n,j), y2(n,j)
          goto 61
62        continue
          nRef = n
          close(20)
       enddo
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = (y1(1:n,1)+y1(1:n,2)+y1(1:n,3))/3.0
       tempIm(1:nRef) =  (y2(1:n,1)+y2(1:n,2)+y2(1:n,3))/3.0

    case("sio2")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"sio2.txt"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       n = 0
71     continue
       read(20,'(a)',end=72) textline
       if (textLine(1:1) == "#") goto 71
       n = n + 1
       read(textline,*) x(n),y1(n,1), y2(n,1)
       goto 71
72     continue
       nRef = n
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = y1(1:n,1)
       tempIm(1:nRef) = y2(1:n,1)
       close(20)

    case("iron")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"iron.dat"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       n = 0
81     continue
       read(20,'(a)',end=82) textline
       if (textLine(1:1) == "#") goto 81
       n = n + 1
       read(textline,*) x(n),y1(n,1), y2(n,1)
       goto 81
82     continue
       nRef = n
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = y1(1:n,1)
       tempIm(1:nRef) = y2(1:n,1)
       close(20)

    case("alumina")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"alumina.dat"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       n = 0
91     continue
       read(20,'(a)',end=92) textline
       if (textLine(1:1) == "#") goto 91
       n = n + 1
       read(textline,*) x(n),y1(n,1), y2(n,1)
       goto 91
92     continue
       nRef = n
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = y1(1:n,1)
       tempIm(1:nRef) = y2(1:n,1)
       close(20)

    case("glassy_pyroxene")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"glassy_pyroxene.dat"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       n = 0
93     continue
       read(20,'(a)',end=94) textline
       if (textLine(1:1) == "#") goto 93
       n = n + 1
       read(textline,*) x(n),y1(n,1), y2(n,1)
       goto 93
94     continue
       nRef = n
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = y1(1:n,1)
       tempIm(1:nRef) = y2(1:n,1)
       close(20)

    case("olivine_glass")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"olivine_glass.dat"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       n = 0
95     continue
       read(20,'(a)',end=96) textline
       if (textLine(1:1) == "#") goto 95
       n = n + 1
       read(textline,*) x(n),y1(n,1), y2(n,1)
       goto 95
96     continue
       nRef = n
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = y1(1:n,1)
       tempIm(1:nRef) = y2(1:n,1)
       close(20)


    case("pinteISM")
       call unixGetenv("TORUS_DATA", dataDirectory, i)
       filename = trim(dataDirectory)//"/"//"pinteISM.dust"
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       nRef = 38
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       do i = nRef, 1, -1
          read(20,*) lamRef(i), tempReal(i), tempIm(i)
       enddo
       close(20)

    case DEFAULT
       filename = trim(graintype)
       if (writeoutput) write(*,'(a,a)') "Reading grain properties from: ",trim(filename)
       open(20,file=filename,status="old",form="formatted")
       n = 0
501     continue
       read(20,'(a)',end=502) textline
       if (textLine(1:1) == "#") goto 501
       n = n + 1
       read(textline,*) x(n),y1(n,1), y2(n,1)
       goto 501
502     continue
       nRef = n
       allocate(lamRef(1:nRef))
       allocate(tempIm(1:nRef))
       allocate(tempReal(1:nRef))
       lamRef(1:nRef) = x(1:n)
       tempReal(1:nRef) = y1(1:n,1)
       tempIm(1:nRef) = y2(1:n,1)
       close(20)

    end select

    firstTime = .true.
    do i = 1, nLambda 
       if (lambda(i)*real(angsToMicrons) < lamRef(1)) then
          mReal(i) = tempReal(1)
          mImg(i) = tempIm(1)
       elseif (((lambda(i)*real(angsToMicrons)) >= lamRef(1)) .and. &
            (lambda(i)*real(angsToMicrons) <= lamRef(nRef))) then


          call locate(lamRef, nRef, lambda(i)*real(angsToMicrons), j)
          t = real((lambda(i)*angsToMicrons - lamRef(j))/(lamRef(j+1) - lamRef(j)))
          mReal(i) = tempReal(j) + t * (tempReal(j+1) - tempReal(j))
          mImg(i) = tempIm(j) + t * (tempIm(j+1) - tempIm(j))         
       else
          if (firstTime) then
             call writeWarning("Extrapolating grain properties")
             firstTime = .false.
          endif
          dydx = (log10(tempReal(nref)) - log10(tempReal(nRef-1))) / &
               (log10(lamRef(nref))-log10(lamRef(nRef-1)))
          mReal(i) = real(log10(tempReal(nref)) + dydx * &
               (log10(lambda(i)*angsToMicrons) - log10(lamRef(nRef))))
          mReal(i) = real(10.d0**mreal(i))
          dydx = (log10(tempIm(nref)) - log10(tempIm(nRef-1))) / &
               (log10(lamRef(nref))-log10(lamRef(nRef-1)))
          mImg(i) = real(log10(tempIm(nref)) + dydx * &
               (log10(lambda(i)*angsToMicrons) - log10(lamRef(nRef))))
          mImg(i) = real(10.d0**mImg(i))
       endif


       if (porousFillingFactor > 0.d0) then
          if (writeoutput) write(*,*) "before ",mreal(i),mimg(i)
          call refractiveIndexToPermittivity(mReal(i), mImg(i), epsilonReal, epsilonImg)
          epsilonMatrix = cmplx(epsilonReal, epsilonImg)
          epsilonInclusion = cmplx(1.e0, 0.e0)
          call effectiveMedium(epsilonEffective, epsilonMatrix, epsilonInclusion, porousFillingFactor)
          epsilonReal = real(epsilonEffective)
          epsilonImg = aimag(epsilonEffective)
          call permittivityToRefractiveIndex(epsilonReal, epsilonImg, mreal(i), mImg(i))
          if (writeoutput) write(*,*) "after ",mreal(i),mimg(i)
       endif

    enddo



  end subroutine getRefractiveIndex

  subroutine fillGridMie(grid, aMin, aMax, a0, qDist, pDist, fillingFactor,&
       ngrain, abundance, grainname, thisDust)
!DEC$ NOOPTIMIZE
    use inputs_mod, only : grainFrac, grainDensity
! This compiler directive disables optimisation in this subroutine, as  
! ifort 12 was incorrectly setting up grid%oneKappaAbs and grid%oneKappaSca. 
    use mieDistCrossSection_mod, only: mieDistCrossSection, mieSingleCrossSection
#ifdef MPI
    use mpi
#endif


    implicit none
    type(GRIDTYPE) :: grid
    integer :: thisDust
    real :: aMin, aMax,a0, qDist, pDist
    real :: fillingFactor
    real, allocatable :: sigmaAbs(:), sigmaSca(:), sigmaExt(:)
    real, allocatable :: mReal(:), mImg(:)          ! size = nlamda
    real, allocatable :: mReal2D(:,:), mImg2D(:,:)  ! size = ngrain x nlambda
    real :: meanParticleMass
    real :: rayleigh, gsca
    integer, intent(in) :: ngrain  ! number of grain types
    real, intent(in) :: abundance(:)   ! relative abundance of grains
    character(len=*) :: grainname(:)   ! names of grains available
    real :: sig_ext, sig_scat, sig_abs
    real :: total_abundance
    integer :: ilam_beg, ilam_end
    character(len=80) :: albedoFilename
    integer :: i, j, k

#ifdef MPI
    real, allocatable :: tempArray(:)
    integer :: np, n_rmdr, m, ierr
#endif



    allocate(sigmaAbs(1:grid%nLambda))
    allocate(sigmaSca(1:grid%nLambda))
    allocate(sigmaExt(1:grid%nLambda))


    call writeInfo("NEW: Filling grid with mie cross-sections...", TRIVIAL)
    call writeInfo("Dust law: n(a) = const * a^-q * Exp( -(a/a0)^p )", TRIVIAL)
    call writeInfo("          where  amin < a < amax", TRIVIAL)
    call writeFormatted("(a,f10.3)","    amin  = ",  aMin, TRIVIAL)
    call writeFormatted("(a,f10.3)","    amax  = ",  aMax, TRIVIAL)
    call writeFormatted("(a,e12.3)","    a0    = ",  a0, TRIVIAL)
    call writeFormatted("(a,e12.3)","    qDist = ",  qDist, TRIVIAL)
    call writeFormatted("(a,e12.3)","    pDist = ",  pDist, TRIVIAL)
    call writeFormatted("(a,f10.3)"," Porosity = ",  fillingFactor, TRIVIAL)
    call writeFormatted("(a,f10.3)"," aMedian  = ", getMedianSize(aMin, aMax, a0, qDist, pDist), TRIVIAL)


    allocate(mReal(1:grid%nLambda))
    allocate(mImg(1:grid%nLambda))

    ! Synthetic grains

    ! quick test for zero total abundance.
    total_abundance = SUM(abundance)
    if ( total_abundance <= 0.0 ) then
       write(*,*) "Error:: total_abundance <= 0.0 in  grain_mod::fillGridMie."
       write(*,*) "  ==> You probably forgot to assign abundance in your parameter file!"
       write(*,*) "  ==> Exiting the program ... "
       stop 
    end if

    ! allocate mem for temp arrays
    allocate(mReal2D(1:ngrain, 1:grid%nLambda))
    allocate(mImg2D(1:ngrain, 1:grid%nLambda))
    ! initializing the values
    mReal2D(:,:) = 0.0; mImg2D(:,:) = 0.0

    ! Find the index of refractions for all types of grains available
    do j = 1, ngrain
       call getRefractiveIndex(grid%lamArray, grid%nLambda, grainname(j), mReal, mImg, FillingFactor)
       mReal2d(j,:) = mReal(:)  ! copying the values to a 2D maxtrix
       mImg2D(j,:)  = mImg(:)   ! copying the values to a 2D maxtrix            
    end do

    ! finding the cross sections
    sigmaExt(:) = 0.0; sigmaAbs(:)=0.0; sigmaSca(:)=0.0 ! initializing the values

!    if (writeoutput) open(20,file="albedo.dat",form="formatted",status="unknown")
!    if (writeoutput) open(21,file="gfactor.dat",form="formatted",status="unknown")

    ilam_beg = 1
    ilam_end = grid%nLambda
#ifdef MPI
    ! Set the range of index for a photon loop used later.     


    np = nThreadsGlobal

    if (np <= grid%nLambda) then
       
       n_rmdr = MOD(grid%nLambda,np)
       m = grid%nLambda/np
       
       if (myRankWorldGlobal .lt. n_rmdr ) then
          ilam_beg = (m+1)*myRankWorldGlobal + 1
          ilam_end = ilam_beg + m
       else
          ilam_beg = m*myRankWorldGlobal + 1 + n_rmdr
          ilam_end = ilam_beg + m -1
       end if
    else
       if (myrankWorldGlobal <= (grid%nLambda-1)) then
          ilam_beg = myRankWorldGlobal+1
          ilam_end = myRankWorldGlobal+1
       else
          ilam_beg = -1
          ilam_end = -1
       endif
    endif
       
#endif

    sigmaExt = 0.0
    sigmaAbs = 0.0
    sigmaSca = 0.0
    if (ilam_beg >= 1) then
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP SHARED(ilam_beg, ilam_end, nGrain, amin, amax, a0, qdist, pdist, grid) &
    !$OMP SHARED(mreal2d, mimg2d) &
    !$OMP SHARED(sigmaExt, sigmaAbs, sigmaSca, abundance, total_abundance) &
    !$OMP PRIVATE(i,j,sig_ext, sig_scat, sig_abs, gsca)
    !$OMP DO SCHEDULE(DYNAMIC)
       do i = ilam_beg, ilam_end
          do j = 1, ngrain
             if (aMin /= aMax) then
                call mieDistCrossSection(aMin, aMax, a0, qDist, pDist, grid%lamArray(i), &
                     mReal2D(j,i), mImg2D(j,i), sig_ext, sig_scat, sig_abs, gsca)
             else
                call mieSingleCrossSection(aMin, grid%lamArray(i), &
                     mReal2D(j,i), mImg2D(j,i), sig_ext, sig_scat, sig_abs, gsca)
             endif
             ! Weighting the cross section according to their abundance...            
             sigmaExt(i) = sig_ext*abundance(j)+ sigmaExt(i)
             sigmaAbs(i) = sig_abs*abundance(j)+ sigmaAbs(i)
             sigmaSca(i) = sig_scat*abundance(j)+ sigmaSca(i)
          end do
          sigmaExt(i) =    sigmaExt(i)/total_abundance 
          sigmaAbs(i) =    sigmaAbs(i)/total_abundance 
          sigmaSca(i) =    sigmaSca(i)/total_abundance 
       end do
    !$OMP END DO
    !$OMP END PARALLEL
    else
       ilam_beg = 1
       ilam_end = 1
    endif
#ifdef MPI
    allocate(tempArray(1:grid%nLambda))
    tempArray = 0.

    call MPI_ALLREDUCE(sigmaExt, tempArray, grid%nLambda, MPI_REAL,&
         MPI_SUM, MPI_COMM_WORLD, ierr)
    sigmaExt = tempArray

    tempArray = 0.
    call MPI_ALLREDUCE(sigmaAbs, tempArray, grid%nLambda, MPI_REAL,&
         MPI_SUM, MPI_COMM_WORLD, ierr)
    sigmaAbs = tempArray
    
    tempArray = 0.
    call MPI_ALLREDUCE(sigmaSca, tempArray, grid%nLambda, MPI_REAL,&
         MPI_SUM, MPI_COMM_WORLD, ierr)
    sigmaSca = tempArray
    deallocate(tempArray)
#endif
    if (.not.grid%oneKappa) then
       if (grid%adaptive) then
          if (writeoutput) write(*,'(a,i3)') "Filling AMR grid with mie cross sections...",grid%nLambda
          call fillAMRgridMie(grid%OctreeRoot, sigmaSca, sigmaAbs, grid%nLambda)
       endif

       if (grid%cartesian) then

          do i = 1, grid%nx
             do j = 1, grid%ny
                do k = 1, grid%nz


                   if (grid%inUse(i,j,k)) then
                      grid%kappaAbs(i,j,k,1:grid%nLambda) = sigmaAbs  * grid%rho(i,j,k)
                      grid%kappaSca(i,j,k,1:grid%nLambda) = sigmaSca  * grid%rho(i,j,k)
                   endif

                   !      write(*,*) grid%kappaAbs(i,j,k,1:grid%nLambda),grid%kappaSca(i,j,k,1:grid%nLambda), grid%rho(i,j,k)
                enddo
             enddo
          enddo
       endif

       if (grid%polar) then
          do i = 1, grid%nr
             do j = 1, grid%nmu
                do k = 1, grid%nPhi

                   grid%kappaAbs(i,j,k,1:grid%nLambda) = sigmaAbs * grid%rho(i,j,k)
                   grid%kappaSca(i,j,k,1:grid%nLambda) = sigmaSca * grid%rho(i,j,k)

                enddo
             enddo
          enddo

       endif

       where(grid%kappaAbs < 1.e-25) grid%kappaAbs = 1.e-25
       where(grid%kappaSca < 1.e-25) grid%kappaSca = 1.e-25


       grid%kappaAbs = grid%kappaAbs * 1.e10
       grid%kappaSca = grid%kappaSca * 1.e10
    else
       call writeFormatted("(a,i4)", "Filling the oneKappa arrays: ",grid%nLambda, TRIVIAL)

       meanParticleMass = 0.
       do i = 1, ngrain
          meanParticleMass = meanParticleMass + getMeanMass2(fillingFactor, aMin, aMax, a0, qDist, pDist, &
               grainname(i),graindensity(i))*abundance(i)
       enddo
       grid%oneKappaAbs(thisDust,1:grid%nLambda) = (sigmaAbs(1:grid%nLambda) * 1.e10)/meanParticleMass
       grid%oneKappaSca(thisDust,1:grid%nLambda) = (sigmaSca(1:grid%nLambda) * 1.e10)/meanParticleMass

       write(albedoFilename,'(a,i2.2,a)') "albedo",thisDust,".dat"
       if (writeoutput) open(20,file=albedoFilename,form="formatted",status="unknown")
       if (writeoutput) write(20,'(a120)') "# Columns are: wavelength (microns), kappa ext (cm^2 g^-1), &
           &    kappa abs (cm^2 g^-1), kappa sca (cm^2 g^-1), albedo"
       if (writeoutput) write(20,*) "# Note that the opacities are per gram of gas"

       do i = 1, grid%nLambda
          rayleigh = real((8.*pi**2)/(grid%lamArray(i)*angstromtocm)* &
               aimag((cmplx(mreal(i),mimg(i))**2-cmplx(1.,0.))/(cmplx(mreal(i),mimg(i))**2+cmplx(2.,0.)))*(amin*microntocm)**3)
          rayleigh = rayleigh / meanParticleMass
          if (writeoutput) &
               write(20,'(1p,5e13.4)') grid%lamArray(i)*angstomicrons,&
               grainFrac(thisDust)*(grid%oneKappaAbs(thisdust,i)+grid%oneKappaSca(thisdust,i))/1.e10, &
               grainFrac(thisDust)*grid%oneKappaAbs(thisdust,i)/1.e10,grainFrac(thisDust)*grid%oneKappaSca(thisdust,i)/1.e10, &
               grid%oneKappaSca(thisdust,i)/(grid%oneKappaAbs(thisdust,i)+grid%oneKappaSca(thisdust,i))
       enddo
       if (writeoutput) close(20)




    endif
    deallocate(sigmaAbs, sigmaSca)
    call writeInfo("mie cross-sections done. Note 10^10 factor", TRIVIAL)
  end subroutine fillGridMie




  recursive subroutine fillAMRgridMie(thisOctal, sigmaSca, sigmaAbs, nLambda)
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: nLambda
    real :: sigmaSca(*), sigmaAbs(*)
    integer :: subcell, i
!    write(*,*) subcell,sigmasca(1),sigmaabs(1),nlambda

    do subcell = 1, 8
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillAMRgridMie(child, sigmaSca, sigmaAbs, nLambda)
                exit
             end if
          end do
       else
          if (thisOctal%inFlow(subcell)) then
             thisOctal%kappaAbs(subcell,1:nLambda)   &
                  = sigmaAbs(1:nLambda) * thisOctal%rho(subcell) * 1.e10
             thisOctal%kappaSca(subcell,1:nLambda)   &
                  = sigmaSca(1:nLambda) * thisOctal%rho(subcell) * 1.e10
          else
             thisOctal%kappaAbs(subcell,1:nLambda)   &
                  = 0.0
             thisOctal%kappaSca(subcell,1:nLambda)   &
                  = 0.0
          end if

       endif
    enddo
  end subroutine fillAMRgridMie

  subroutine dustPropertiesfromFile(filename, nlambda, lambda, kappaAbs, kappaSca, gfac)
    use utils_mod, only: logInt
    implicit none
    character(len=*) :: filename
    integer :: nlambda
    real :: lambda(:)
    real :: kappaAbs(:), kappaSca(:)
    real :: sigmaExt(2000),sigmaSca(2000), kappa(2000), albedo(2000), tlam(2000)
    real :: tSca(2000), tAbs(2000), sigmaAbs(2000), junk
    real(double) :: gfac(2000), tgfac(2000), rjunk
    character(len=40) :: filetype
    character(len=80) :: message, junkchar
    integer :: npts, i, j

    write(message,'(a,a)') "Reading dust properties from: ",trim(filename)
    call writeInfo(message, TRIVIAL)
    open(20, file=filename, status="old", form="formatted")
    read(20,'(a)') filetype

    select case (filetype)


    case("MW")
       do i = 1, 80
          read(20,'(a80)') junkchar
       enddo
       do npts = 1066,1,-1
          read(20,*) tlam(npts),albedo(npts),tgfac(npts),rjunk,kappa(npts)
!          if (writeoutput) write(*,*) tlam(npts),albedo(npts),tgfac(npts),rjunk,kappa(npts)
       enddo
       close(20)
       npts = 1066
       tabs(1:npts) = kappa(1:npts)
       tsca(1:npts) = tabs(1:npts)*albedo(1:npts)/(1.-albedo(1:npts))
       tlam(1:npts) = tlam(1:npts)*1e4
       
    case("kenny")
       npts = 1
10     read(20,*,end=20) tlam(npts), sigmaExt(npts),sigmaSca(npts),kappa(npts)
       npts = npts + 1
       goto 10
20     continue
       npts = npts - 1
       close(20)
       albedo(1:npts) = sigmaSca(1:npts) / sigmaExt(1:npts)

       tAbs(1:npts) = (1.-albedo(1:npts))*kappa(1:npts)
       tSca(1:npts) = albedo(1:npts)*kappa(1:npts)
       
       tlam(1:npts) = tlam(1:npts) * 1.e4 ! microns to angstrom
    case("jenny")
       npts = 1
30     read(20,*,end=40) tlam(npts), sigmaAbs(npts),sigmaSca(npts)
       npts = npts + 1
       goto 30
40     continue
       npts = npts - 1
       close(20)
       tAbs(1:npts) = sigmaAbs(1:npts)
       tSca(1:npts) = sigmaSca(1:npts)
       tlam(1:npts) = tlam(1:npts) * 1.e4 ! microns to angstrom

    case("trust")
       npts = 1201
       read(20,*) junkchar
       read(20,*) junkchar
       read(20,*) junkchar
       read(20,*) junkchar

       do i = 1, npts
          read(20,*) tlam(i), junk, junk, junk, kappa(i), albedo(i), tgfac(i)
       enddo

       close(20)
       tAbs(1:npts) = (1.0-albedo(1:npts)) * kappa(1:npts)
       tSca(1:npts) = albedo(1:npts) * kappa(1:npts)
       tlam(1:npts) = tlam(1:npts) * 1.e4 ! microns to angstrom

    case DEFAULT
       write(*,'(a)') "! Dust properties file has unknown type",trim(filetype)
       stop
       
    end select

    do i = 1, nLambda
       call locate(tlam,npts,lambda(i),j)
       kappaAbs(i) = logint(lambda(i), tlam(j), tlam(j+1), tAbs(j), tAbs(j+1))*1.e10
       kappaSca(i) = logint(lambda(i), tlam(j), tlam(j+1), tSca(j), tSca(j+1))*1.e10


!       kappaAbs(i) = 1.d10*(tabs(j) + (lambda(i) - tlam(j)) * (tabs(j+1) - tabs(j))/ ( tlam(j+1) - tlam(j)))
!       kappaSca(i) = 1.d10*(tsca(j) + (lambda(i) - tlam(j)) * (tsca(j+1) - tsca(j))/ ( tlam(j+1) - tlam(j)))
       
       gfac(i) = tgfac(j) + (lambda(i)-tlam(j)) * (tgfac(j+1)-tgfac(j)) / (tlam(j+1)-tlam(j))
!       if (writeoutput) write(*,*) lambda(i), kappaAbs(i)+ kappaSca(i), kappaAbs(i), kappaSca(i)
    enddo

  end subroutine dustPropertiesfromFile


  recursive subroutine fillDustUniform(grid, thisOctal)

    use inputs_mod, only : nDustType, grainFrac, rSublimation
    use inputs_mod, only : ttauriDisc
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: rVec
    real(double) :: r, fac
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustUniform(grid, child)
                exit
             end if
          end do
       else
          fac = 1.

!          if (((grid%geometry == "shakara").and.(.not.variableDustSublimation)).or. &
           if ((grid%geometry == "ttauri").and.ttauriDisc) then
              rVec = subCellCentre(thisOctal, Subcell)
              r = sqrt(rVec%x**2 + rVec%y**2)
              if (r < 1.01d0*rSublimation) then
                 fac = (1.01d0*rSublimation - r)/(0.01d0*rSublimation)
                 fac = exp(-fac)
              endif
           endif

          thisOctal%dustTypeFraction(subcell,1:nDustType) = grainFrac(1:nDustType) * fac
       end if
    end do

  end subroutine fillDustUniform

  recursive subroutine setupOrigDustFraction(thisOctal)
    use inputs_mod, only : nDustType
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call setupOrigDustFraction(child)
                exit
             end if
          end do
       else
          if (.not.associated(thisOctal%origDustTypeFraction)) &
               allocate(thisOctal%origDustTypeFraction(1:thisOctal%maxChildren,1:nDustType))
          thisOctal%origdustTypeFraction(subcell,:) =  thisOctal%dustTypeFraction(subcell,:) 
       end if
    end do

  end subroutine setupOrigDustFraction

  recursive subroutine fillDustShakara(grid, thisOctal, dustmass)

    use inputs_mod, only : rSublimation, nDustType, curvedInnerEdge, grainFrac, usemultidust, router, xidust, amin
    use inputs_mod, only : betaDisc, height, rinner
    use octal_mod, only : cellVolume
    use density_mod

    
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: rVec
    real(double) :: r, z, hr
    real(double) :: dustMass, cellMass, thisHeight
    real(double) :: fac, rhoFid, rho, thisRsub, gasheight
    integer :: subcell, i, idust

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustShakara(grid, child, dustMass)
                exit
             end if
          end do
       else
         
          rVec = VECTOR(rSublimation*1.001d0, 0.d0, 0.d0)
          rhoFid = shakaraSunyaevDisc(rVec, grid)

          thisOctal%DustTypeFraction(subcell,:) = 1.d-10
          rVec = subcellCentre(thisOctal, subcell)
          rho = shakaraSunyaevDisc(rVec, grid)
          thisRsub = 1.01d0 * rSublimation * max(1.d0,(1.d0/(rho/rhoFid)**0.0195)**2)
          r = sqrt(rVec%x**2+rVec%y**2)
          z = rVec%z
          thisOctal%dustTypeFraction(subcell,1:nDustType) = grainFrac(1:nDustType)
!          write(*,*) r/1496.,thisRsub/1496.d0
!          if (modulus(rVec) < rSublimation) then
!             thisOctal%dustTypeFraction(subcell,1:nDustType) = 1.d-25
!          endif

          if (curvedInnerEdge.and.(r < thisRsub).and.(modulus(rVec) < 2.d0*rsublimation)) then
             fac = (thisRsub-r)/(0.002d0*rSublimation)
             thisOctal%dustTypeFraction(subcell,1:nDustType) = max(1.d-20,grainFrac(1:nDustType)*exp(-fac))
          endif


          if (usemultidust) then

             thisOctal%dustTypeFraction(subcell,:) = 1.e-30
             
             if ((r > rinner).and.(r < router)) then
             
                hr =  height*(r/(100.d0*autocm/1.d10))**betaDisc
                if (abs(z/hr) < 7.0) then
                   thisOctal%dustTypeFraction(subcell,nDustType) = 1.e-30
!                   sigma = rho0 * (r/rinner)**(-alphaDisc) * (hr * 1.d10) * sqrt(2.d0*pi)
                   do iDust = 1, 10
!                      grainDensity(idust) = getCompositeGrainDensity(graintype(idust))
                      gasHeight =  height*(r/(100.d0*autocm/1.d10))**betaDisc
!                      f = alphaViscosity * sigma / (sqrt(6.*pi) * (amid(idust)*microntocm) * grainDensity(idust))
!                      !                write(*,*) "f ",f, sqrt(f/(f+1.d0))
                      thisHeight = gasHeight * (amin(idust)/amin(1))**(-xidust)
                      fac = exp(-0.5d0*(z/thisHeight)**2 + 0.5d0*(z/gasHeight)**2)
                      thisOctal%dustTypeFraction(subcell,idust) = max(fac,1.d-30)
                   enddo
                endif
                thisOctal%dustTypeFraction(subcell,1:10) = max(1.d-30,thisOctal%dustTypeFraction(subcell,1:10) * &
                  grainFrac(1:10))
             endif
          endif
             
          cellMass = cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell)
          dustMass = dustMass + SUM(thisOctal%dustTypeFraction(subcell,1:nDustType))*cellMass

       end if
    end do

  end subroutine fillDustShakara


  recursive subroutine fillDustJaehan(grid, thisOctal, dustmass)
    use inputs_mod, only : grainFrac, nDustType
    use octal_mod, only : cellVolume
    use density_mod

    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: cellMass, dustMass

    integer :: subcell, i
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustJaehan(grid, child, dustMass)
                exit
             end if
          end do
       else
         
          

          thisOctal%dustTypeFraction(subcell,1:nDustType) = grainFrac(1:nDustType)
          thisOctal%dustTypeFraction(subcell,3) = thisOctal%etaline(subcell)/(cellVolume(thisOctal,subcell)*1.d30)&
               /thisOctal%rho(subcell)
          cellMass = cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell)
          dustMass = dustMass + SUM(thisOctal%dustTypeFraction(subcell,1:nDustType))*cellMass

       end if
    end do

  end subroutine fillDustJaehan

  subroutine reportMasses(grid)
    use inputs_mod, only : nDustType
    type(GRIDTYPE) :: grid
    integer :: i
    real(double) :: gasMass, dustMass

    gasMass = 0.d0
    call findTotalMass(grid%octreeRoot, gasMass)

    if (writeoutput) then
       write(*,'(a,1pe12.5)') "Gas mass (solar masses): ",gasMass/msol
       do i = 1, nDustType
          dustMass = 0.d0
          call findDustMassSingle(grid, grid%octreeRoot, dustMass, i)
          write(*,'(a,i2,a,1pe12.5)') "Dust mass ",i, " (solar masses): ",dustMass/mSol
       enddo
    endif
  end  subroutine reportMasses

  subroutine getDustMasses(grid, dustMass)
    use inputs_mod, only : nDustType
    type(GRIDTYPE) :: grid
    integer :: i
    real(double) :: dustMass(:)

    do i = 1, nDustType
       dustMass(i) = 0.d0
       call findDustMassSingle(grid, grid%octreeRoot, dustMass(i), i)
    enddo
  end  subroutine getDustMasses

  recursive subroutine findDustMassSingle(grid, thisOctal, dustmass, iDust)

    use octal_mod, only : cellVolume
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: dustMass, cellMass
    integer :: subcell, i, iDust

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findDustMassSingle(grid, child, dustMass, iDust)
                exit
             end if
          end do
       else

          cellMass = cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell)
          dustMass = dustMass + thisOctal%dustTypeFraction(subcell,iDust)*cellMass
          
       endif

    end do

  end subroutine findDustMassSingle

  recursive subroutine findDustMass(grid, thisOctal, dustmass)

    use inputs_mod, only : nDustType
    use octal_mod, only : cellVolume
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: dustMass, cellMass
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findDustMass(grid, child, dustMass)
                exit
             end if
          end do
       else

          cellMass = cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell)
          dustMass = dustMass + SUM(thisOctal%dustTypeFraction(subcell,1:nDustType))*cellMass
          
       endif

    end do

  end subroutine findDustMass

    recursive subroutine findDustWeightedTemperature(grid, thisOctal, dustmass, tempsum, idust)

    use octal_mod, only : cellVolume
    type(gridtype) :: grid
    integer :: idust
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: dustMass, cellMass, tempSum
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call findDustWeightedTemperature(grid, child, dustMass, tempSum, idust)
                exit
             end if
          end do
       else
          cellMass = cellVolume(thisOctal, subcell) * 1.d30 * thisOctal%rho(subcell) * thisOctal%dustTypeFraction(subcell,idust)
          dustMass = dustMass + cellMass
          tempSum = tempSum + cellMass * thisOctal%temperature(subcell)
          
       endif

    end do

  end subroutine findDustWeightedTemperature

  subroutine reportDustMassWeightedTemperature(grid)
    use inputs_mod, only : nDustType
    type(GRIDTYPE) :: grid
    integer :: iDust
    real(double) :: thisTemp, tempSum, dustMass

    do iDust = 1, nDustType
       dustMass = 0.
       tempSum = 0.
       call findDustWeightedTemperature(grid, grid%octreeRoot, dustMass, tempSum, iDust)
       thisTemp = tempSum/dustMass
       if (writeoutput) write(*,*) "Dust mass weighted average temperature for ",idust, " is ",thisTemp
    enddo
  end subroutine reportDustMassWeightedTemperature

  
  recursive subroutine normalizeDustFractions(grid, thisOctal, thisDustMass, requiredDustMass, nDust)

    type(gridtype) :: grid
    integer :: ndust
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: thisDustMass(:), requiredDustMass(:)
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call normalizeDustFractions(grid, child, thisDustMass, requiredDustMass, ndust)
                exit
             end if
          end do
       else
          thisOctal%dustTypeFraction(subcell,1:nDust) = thisOctal%dustTypeFraction(subcell,1:nDust) * &
               requiredDustMass(1:nDust) / thisDustMass(1:nDust)
       end if
    end do

  end subroutine normalizeDustFractions

  recursive subroutine emptyDustCavity(thisOctal, position, radius)

    use inputs_mod, only : nDustType, grainFrac
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: position
    real(double) :: r, radius
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call emptyDustCavity(child, position, radius)
                exit
             end if
          end do
       else
          if (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) cycle
          r = modulus(subcellCentre(thisOctal, subcell) - position)
          if (r < radius) then
             thisOctal%dustTypeFraction(subcell,1:nDustType) = 1.d-6
          else
             thisOctal%dustTypeFraction(subcell,1:nDustType) = grainFrac(1:nDustType)
          endif
       end if
    end do

  end subroutine emptyDustCavity


  recursive subroutine sublimateDustWR104(thisOctal)
    use inputs_mod, only : tThresh, grainFrac
    type(OCTAL), pointer :: thisOctal, child
    integer :: subcell, i
    real(double) :: temperature
    real(double), parameter :: subRange = 1.d0
    real(double) :: frac
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sublimateDustWR104(child)
                exit
             end if
          end do
       else
          temperature = thisOctal%temperature(subcell)

          frac = 1.d0
          if (temperature > tThresh) then
             frac = exp(-dble((temperature-tThresh)/subRange))
          endif
          thisOctal%dustTypeFraction(subcell,1) = grainFrac(1)*max(1.d-20, frac)
       endif
    end do
  end subroutine sublimateDustWR104

  recursive subroutine sublimateDust(grid, thisOctal, totFrac, nFrac, tauMax, subTemp, minLevel)

    use inputs_mod, only : grainFrac, nDustType, tThresh, tSub, decoupleGasDustTemperature, tsubpower, subrange
    type(gridtype) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real :: totFrac
    real :: tauMax
    real(double), optional :: subTemp, minLevel
    real(double) :: smallVal
    integer :: nFrac
    real(double) :: frac, newFrac, deltaFrac, thistau
    real ::  temperature, sublimationTemp
    real :: underCorrect 
    integer :: ilambda
    real(double) :: kappaSca, kappaAbs
    integer :: subcell, i, j

    underCorrect = 1.

    kappaSca = 0.d0; kappaAbs = 0.d0
    if (present(minLevel)) then
       smallVal = minLevel
    else
       smallVal = 1.d-25
    endif

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sublimateDust(grid, child, totFrac,nFrac, tauMax, subTemp, minLevel)
                exit
             end if
          end do
       else
          do j = 1, nDustType
             if (decoupleGasDustTemperature) then
                temperature = real(thisOctal%tDust(subcell))
             else
                temperature = thisOctal%temperature(subcell)
             endif
             if (present(subTemp)) then
                sublimationTemp = real(subTemp)
             else
                sublimationTemp = real(max(700.d0,tSub(j) * thisOctal%rho(subcell)**(tsubpower(j))))
             endif

             if (tThresh /= 0.) sublimationTemp = tThresh
             if (temperature < sublimationTemp) newFrac = 1.

             
             if (temperature >= sublimationTemp) then
                newfrac = 0.5d0 * exp(-dble((temperature-sublimationtemp)/subRange))
             else
                newfrac = 1.d0 - 0.5d0 * exp(dble((temperature-sublimationtemp)/subRange))
             endif


             newfrac = max(newfrac,smallVal)
             
             deltaFrac = newFrac - thisOctal%oldFrac(subcell)
             
             frac = thisOctal%oldFrac(subcell) + underCorrect * deltaFrac
             
             frac = max(frac, smallVal)


             if (.not.associated(thisOctal%origDustTypeFraction)) then
                allocate(thisOctal%origDustTypeFraction(1:thisOctal%maxChildren,1:nDustType))
                do i = 1, thisOctal%maxChildren
                   thisOctal%origDustTypeFraction(i,1:nDustType) = grainFrac(1:nDustType)
                enddo
             endif
             thisOctal%dustTypeFraction(subcell,j) = thisOctal%origDustTypeFraction(subcell,j) * frac



             call locate(grid%lamArray, grid%nLambda, 5500., iLambda)
             call returnKappa(grid, thisOctal, subcell, ilambda, kappaSca=kappaSca, kappaAbs=kappaAbs)

             thisTau = (kappaAbs+kappaSca)*thisOctal%subcellSize
             if (thisTau > tauMax) then
                frac = tauMax / thisTau 
                thisOctal%dustTypeFraction(subcell,j) = thisOctal%dustTypeFraction(subcell,j) * frac
             endif

          enddo
!          where (thisOctal%dustTypeFraction(subcell,1:nDustType) < 1.d-25) 
!             thisOctal%dustTypeFraction(subcell,1:nDustType) = 1.d-25
!          end where
          
          if (deltaFrac /= 0.) then
             nfrac = nfrac + 1
             totFrac = totFrac + real(abs(deltaFrac))
          endif
          
          thisOctal%oldFrac(subcell) = real(frac)
       endif
777    continue
    end do

  end subroutine sublimateDust

!!$  subroutine returnScaleHeight(grid, x,  height)
!!$    type(gridtype) :: grid
!!$    integer :: nz
!!$    real :: x, z(10000), height
!!$    real :: subcellsize(10000), temperature(10000)
!!$    real(double) :: rho(10000)
!!$    real(double) :: rho_over_e
!!$    integer :: i, j
!!$    nz  = 0
!!$    temperature = 0.; subcellSize = 0.; rho = 0.0; z=0.d0
!!$    call getTemperatureDensityRundust(grid, z, subcellsize, rho, temperature, x, 0., nz, 1.)
!!$
!!$
!!$    if (rho(2) >= rho(1)) then
!!$       height = 1.e20
!!$       goto 666
!!$    endif
!!$
!!$    rho_over_e = rho(1) / exp(1.d0)
!!$
!!$    if (rho_over_e < rho(nz)) then
!!$       height = z(nz) /1.e10
!!$    else
!!$       j = 1
!!$       do i = 1, nz
!!$          if ((rho_over_e < rho(i)).and.(rho_over_e > rho(i+1))) then
!!$             j = i
!!$             exit
!!$          endif
!!$       enddo
!!$       height = real(z(j) + (z(j+1)-z(j))*(rho_over_e - rho(j))/(rho(j+1)-rho(i)))
!!$       height = height / 1.e10
!!$    endif
!!$666 continue
!!$  end subroutine returnScaleHeight
!!$
!!$  subroutine getTemperatureDensityRundust(grid, zAxis, subcellsize, rho, temperature, xPos, yPos, nz, direction)
!!$    type(GRIDTYPE) :: grid
!!$    type(octal), pointer   :: thisOctal
!!$    integer :: nz
!!$    real(double) :: rho(:)
!!$    real ::temperature(:), zAxis(:), subcellsize(:)
!!$    real :: xPos, yPos
!!$    integer :: subcell
!!$    real(double) :: rhotemp
!!$    real :: temptemp
!!$    real :: direction
!!$    type(VECTOR) :: currentPos, temp
!!$    real :: halfSmallestSubcell
!!$
!!$    nz = 0
!!$    halfSmallestSubcell = real(grid%halfSmallestSubcell)
!!$
!!$    currentPos = VECTOR(xPos, yPos, direction*halfSmallestSubcell)
!!$
!!$    do while(abs(currentPos%z) < grid%ocTreeRoot%subcellsize)
!!$       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
!!$            foundSubcell=subcell, rho=rhotemp, temperature=temptemp)
!!$       thisOctal%chiLine(subcell) = 1.e-30
!!$       !       if (thisOctal%inFlow(subcell)) then
!!$       nz = nz + 1
!!$       temperature(nz) = temptemp
!!$       rho(nz) = rhotemp
!!$       temp = subCellCentre(thisOctal, subcell)
!!$       zAxis(nz) = real(temp%z)
!!$       subcellsize(nz) = real(thisOctal%subcellsize)
!!$       !       endif
!!$       currentPos = VECTOR(xPos, yPos, zAxis(nz)+0.5*direction*thisOctal%subcellsize+direction*halfSmallestSubcell)
!!$       !       else
!!$       !          currentPos = VECTOR(xPos, yPos, grid%octreeRoot%subcellsize+halfSmallestSubcell)
!!$       !       endif
!!$    end do
!!$    zAxis(1:nz) = abs(zAxis(1:nz)) * 1.e10  ! convert to cm
!!$  end subroutine getTemperatureDensityRundust

  recursive subroutine stripDustAway(thisOctal, fac, radius, sourcePos)

    use inputs_mod, only : grainFrac, nDustType
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: fac
    real(double) :: radius, r
    integer :: subcell, i
    type(vector), optional :: sourcePos
    type(vector) :: thisSourcePos

! Use the source position if provided or default to the origin 
    if (present(sourcePos)) then 
       thisSourcePos = sourcePos
    else
       thisSourcePos = vector(0.0d0, 0.0d0, 0.0d0)
    endif

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                if (present(sourcePos)) then 
                   call stripDustAway(child, fac, radius, sourcePos)
                else
                   call stripDustAway(child, fac, radius)
                endif
                exit
             end if
          end do
       else
          r = modulus(subcellCentre(thisOctal, subcell) - thisSourcePos)
          if (r < radius) then
             if (.not.associated(thisOctal%origDustTypeFraction)) then
                allocate(thisOctal%origDustTypeFraction(1:thisOctal%maxChildren,1:nDustType))
                do i = 1, thisOctal%maxChildren
                   thisOctal%origDustTypeFraction(i,1:nDustType) = grainFrac(1:nDustType)
                enddo
             endif
             thisOctal%dustTypeFraction(subcell,1:nDustType) =  thisOctal%origdustTypeFraction(subcell,1:nDustType) * fac
             thisOctal%etaCont(subcell) = tiny(thisOctal%etaCont(subcell))
             if (.not.associated(thisOctal%oldFrac)) then
                allocate(thisOctal%oldFrac(1:thisOctal%maxChildren))
             endif
             thisOctal%oldFrac(subcell) = real(fac)
          endif
       end if
    end do

  end subroutine stripDustAway

  subroutine createRossArray(grid)
    use inputs_mod, only : nDustType
    use atom_mod, only: bnu
    type(GRIDTYPE) :: grid
    integer :: i, j, k
    real(double) :: bNuTot, rosselandKappa, temperature
    real(double) :: dFreq, Freq
    real :: maxTemp

    maxTemp = 3000.

    if (grid%geometry == "wr104") maxTemp = 30000.

    if (grid%nTempRossArray == 0) then
       call writeFatal("nTempRossArray is zero on call to createRossArray")
    endif

    do k = 1, grid%nTempRossArray
       temperature = 3. + (maxTemp-3.)*real(k-1)/real(grid%nTempRossArray-1)
       do j = 1, nDustType
       rosselandKappa = 0.
          Bnutot = 0.
          do i =  grid%nLambda,2,-1
             freq = cSpeed / (grid%lamArray(i)*1.e-8)
             dfreq = cSpeed / (grid%lamArray(i)*1.e-8) - cSpeed / (grid%lamArray(i-1)*1.e-8)
             if ((grid%oneKappaabs(j,i)+grid%oneKappaSca(j,i)) /= 0.) then
!                rosselandKappa = rosselandKappa + (bnu(freq, dble(temperature)) * dFreq / &
!                     (grid%oneKappaabs(j,i)+grid%oneKappaSca(j,i)))
!                bnutot = bnutot + bnu(freq, dble(temperature))*dfreq


                rosselandKappa = rosselandKappa + (dbnubydt(freq, dble(temperature)) * dFreq / &
                     (grid%oneKappaabs(j,i)+grid%oneKappaSca(j,i)))
                bnutot = bnutot + dbnubydt(freq, dble(temperature))*dfreq
             endif
          enddo
          if (rosselandkappa /= 0.) then
             rosselandKappa = (bnutot / rosselandKappa)/1.d10
          endif
          grid%kappaRossArray(j,k) = real(rosselandKappa)
       enddo
       grid%tempRossArray(k) = real(temperature)
    enddo
  end subroutine createRossArray

!!$  real function rtnewtdust(x1,x2,xacc, p1, p2) result(junk)
!!$
!!$    real :: x1, x2, xacc, p1, p2
!!$    integer :: jmax, j
!!$    real ::  dx, f, df
!!$    parameter (jmax=20)
!!$    df =0.; f = 0.; dx = 0.
!!$    junk = 0.5 * (x1+x2)
!!$    do j=1,jmax
!!$       call equation2dust(junk,f,df,p1,p2)
!!$       dx=f/df
!!$       junk=junk-dx
!!$       if((x1-junk)*(junk-x2).lt.0.) then
!!$          write(*,*) 'RTNEWT: jumped out of brackets',p1,p2,junk
!!$          stop
!!$       endif
!!$       if(abs(dx).lt.xacc) return
!!$    enddo
!!$    write(*,*) 'rtnewt exceeding maximum iterations'
!!$  end function rtnewtdust


!!$  subroutine Equation2dust(mu0, eq2, deq2, r, mu)
!!$    real :: r, mu, mu0
!!$    real :: eq2, deq2
!!$
!!$    eq2 = mu0**3 + (r-1.)*mu0 -r*mu
!!$    deq2 = 3.*mu0**2 + r - 1.
!!$
!!$  end subroutine Equation2dust


  real(double) function dbnubydt(nu, T) 
    real(double) :: nu, T, fac1, fac2, fac3
    real(double), parameter :: c1 = (2.d0*hCgs/cSpeed**2)
    real(double), parameter :: c2 = hCgs/kErg

    fac1 = (c1*c2*nu**4/T**2)
    fac3 = c2*nu/T
    if (fac3 < 20.d0) then
       fac2 = exp(fac3)/(exp(fac3)-1.d0)**2
    else
       fac2 = exp(-c2*nu/T)
    endif
    dBnuBydT = fac1 * fac2

  end function dbnubydt

  subroutine parseGrainType(grainString, nTypes, name, abundance)

    character(len=*) :: grainString
    integer :: nTypes
    character(len=*) :: name(:)
    character(len=80) :: tempString
    real :: abundance(:)
    integer :: i, j

    nTypes = 0

    i = index(grainString,":")

    ! not a mixed grain type

    if (i == 0) then
       nTypes = 1
       name(1) = trim(grainString)
       abundance(1) = 1.
       goto 999
    endif
    tempString = grainString
    do while (index(tempString,":") /=0)
       i = index(tempString,":")
       nTypes = nTypes + 1
       name(nTypes) = tempString(1:i-1)
       tempString(1:) = tempString(i+1:)
       if (index(tempString,":") /=0 ) then
          j = index(tempString,":")-1
          read(tempString(1:j), *, err=666) abundance(nTypes)
       else
          j = len(trim(tempString))   
          read(tempString(1:j), *,err=666) abundance(nTypes)
          goto 999
       endif
       tempString(1:) = tempString(j+2:)
    end do

666 continue
    write(*,*) "Error parsing graintype string"
    write(*,*) "Grain string: ",trim(grainString)
    write(*,*) "temp string: ",trim(tempString)
    stop
999 continue
  end subroutine parseGrainType

  subroutine writeDust(dustfile, iDustType, grid, xArray, nLambda, miePhase, nMuMie) 
    use phasematrix_mod, only : phasematrix
    character(len=*) :: dustFile
    type(GRIDTYPE) :: grid
    real :: xarray(:)
    integer :: nLambda, iDustType
    type(PHASEMATRIX) :: miePhase(:,:,:)
    integer :: nMuMie

    if (writeoutput) then
       open(22, file=dustfile,status="unknown", form="unformatted")
       write(22) nLambda
       write(22) xArray(1:nLambda)
       write(22) grid%oneKappaAbs(iDustType, 1:nLambda)
       write(22) grid%oneKappaSca(iDustType, 1:nLambda)
       write(22) miePhase(iDustType,1:nLambda, 1:nMuMie)
       close(22)
    endif

  end subroutine writeDust

  subroutine readDust(dustfile, iDustType, grid, xArray, nLambda, miePhase, nMuMie)
    use phasematrix_mod, only : phasematrix
    character(len=*) :: dustFile
    type(GRIDTYPE) :: grid
    real :: xarray(:)
    integer :: nLambda, iDustType
    type(PHASEMATRIX), pointer :: miePhase(:,:,:)
    integer :: nMuMie

    open(22, file=dustfile,status="old", form="unformatted")
    read(22) nLambda

    read(22) xArray(1:nLambda)
    read(22) grid%oneKappaAbs(iDustType, 1:nLambda)
    read(22) grid%oneKappaSca(iDustType, 1:nLambda)

    read(22) miePhase(iDustType,1:nLambda, 1:nMuMie)
    close(22)
!    grid%oneKappaAbs(iDustType, 1:nLambda) = grid%oneKappaAbs(iDustType, 1:nLambda) * dustTogas    
!    grid%oneKappaSca(iDustType, 1:nLambda) = grid%oneKappaSca(iDustType, 1:nLambda) * dustTogas    

  end subroutine readDust

  subroutine dustComparison(grid, miePhase, nMuMie)
    use mieDistPhaseMatrix_mod
    use phasematrix_mod, only: fillIsotropic, fixMiePhase, PHASEMATRIX, fillHenyey, &
         newDirectionMie
    use inputs_mod, only : nDustType, graintype, ngrain, &
         grainname, x_grain, amin, amax, a0, qdist, pdist, &
         grainFrac, porousFillingFactor
    type(PHASEMATRIX),pointer :: miePhase(:,:,:)
!    type(VECTOR) :: uHat, uNew, vec_tmp
!    real(double) :: cosang
    integer :: nMuMie
    real :: kAbs, kSca
    integer :: i, j
    real(double) :: pol1, pol2

    type(GRIDTYPE) :: grid
    real, pointer  :: xArray(:)
    real(double) ::  albedo
    integer :: nLambda

    if (associated(miePhase)) then
       deallocate(miePhase)
       nullify(miePhase)
    endif

    nLambda = 1
    if (associated(xarray)) deallocate(xarray)
    
    nLambda = 1
    allocate(xarray(1:nLambda))
    if (associated(grid%lamArray)) deallocate(grid%lamArray)
    allocate(grid%lamArray(1:1))
    allocate(grid%oneKappaAbs(1:nDustType,1:1))
    allocate(grid%oneKappaSca(1:nDustType,1:1))
    grid%nLambda = 1
    grid%oneKappa = .true.

    allocate(miePhase(1:nDustType,1:nLambda,1:nMumie)) 

    if (writeoutput) open(25,file="comparison.dat",status="unknown",form="formatted")
    do i = 2, 2

       do j = 1, 1000

          aMax(2) = (0.01 + 9.99 * real(j)/1000.)*1.001


          xArray(1) = 1.22e4
          grid%lamArray(1) = xArray(1)
          call parseGrainType(graintype(i), ngrain, grainname, x_grain)
          call fillGridMie(grid, aMin(i), aMax(i), a0(i), qDist(i), pDist(i), porousFillingFactor(i),&
               ngrain, X_grain, grainname, i)
          kAbs = SUM(grid%oneKappaAbs(1:nDustType,1)*grainFrac(1:nDustType))/1.e10
          kSca = SUM(grid%oneKappaSca(1:nDustType,1)*grainFrac(1:nDustType))/1.e10
          albedo = kSca / (kAbs + kSca)
          pol1 = miePhase(i,1,nMuMie/2)%element(1,2)/miePhase(i,1,nMuMie/2)%element(1,1) * albedo


          xArray(1) = 1.63e4
          grid%lamArray(1) = xArray(1)
          call parseGrainType(graintype(i), ngrain, grainname, x_grain)
          call fillGridMie(grid, aMin(i), aMax(i), a0(i), qDist(i), pDist(i), porousFillingFactor(i),&
               ngrain, X_grain, grainname, i)
          kAbs = SUM(grid%oneKappaAbs(1:nDustType,1)*grainFrac(1:nDustType))/1.e10
          kSca = SUM(grid%oneKappaSca(1:nDustType,1)*grainFrac(1:nDustType))/1.e10
          albedo = kSca / (kAbs + kSca)
          pol2 = miePhase(i,1,nMuMie/2)%element(1,2)/miePhase(i,1,nMuMie/2)%element(1,1) * albedo

          if (writeoutput) write(25,'(1p,5e13.4)') aMax(2), pol1, pol2, -2.5d0*log10(pol1/pol2)
          
       enddo
    enddo
       if (writeoutput) close(25)

  end subroutine dustComparison



  subroutine createDustCrossSectionPhaseMatrix(grid, xArray, nLambda, miePhase, nMuMie)
#ifdef MPI
    use mpi
#endif
    use mieDistPhaseMatrix_mod
    use phasematrix_mod, only: fillIsotropic, fixMiePhase, PHASEMATRIX, fillHenyey, &
         newDirectionMie
    use inputs_mod, only : mie, useDust, dustFile, nDustType, graintype, ngrain, &
         grainname, x_grain, amin, amax, a0, qdist, pdist, &
         kappafilename, isotropicScattering, readmiephase, writemiephase, &
         ttau_disc_on, grainFrac, henyeyGreensteinphaseFunction, porousFillingFactor, inputGfac
    real, allocatable :: mReal2D(:,:), mImg2D(:,:)
    character(len=80) :: miefile
    type(PHASEMATRIX),pointer :: miePhase(:,:,:)
!    type(VECTOR) :: uHat, uNew, vec_tmp
!    real(double) :: cosang
    real, allocatable :: mReal(:,:), mImg(:,:), tmimg(:), tmreal(:)
    integer :: nMuMie
    real(double) :: mu
    real :: total_dust_abundance
    real :: kAbs, kSca
    integer :: i, j, k


    type(GRIDTYPE) :: grid
    real :: xArray(:)
    real(double) :: gfac(2000)
    integer :: nLambda
!    integer :: ilam_beg, ilam_end

#ifdef MPI
!    real, allocatable :: temp(:,:,:,:), tempArray(:), tempArray2(:)
!    integer :: np, n_rmdr, m, ierr, i1,i2
#endif


    ! Note: the first index should be either lambda or mu
    !       in order to speedup the array operations!!!  (RK) 
    if (associated(miePhase)) then
       deallocate(miePhase)
       nullify(miePhase)
    endif

    allocate(miePhase(1:nDustType,1:nLambda,1:nMumie)) 



    if (mie .or. useDust) then
       if (associated(grid%onekappaAbs)) deallocate(grid%onekappaAbs)
       if (associated(grid%onekappaSca)) deallocate(grid%onekappaSca)
       allocate(grid%oneKappaAbs(1:nDustType, 1:nLambda), grid%oneKappaSca(1:nDustType, 1:nLambda))

       do i = 1, nDustType
          if (.not.dustfile(i)) then
             call parseGrainType(graintype(i), ngrain, grainname, x_grain)
             call fillGridMie(grid, aMin(i), aMax(i), a0(i), qDist(i), pDist(i), porousFillingFactor(i),&
                  ngrain, X_grain, grainname, i)
          else
             call dustPropertiesfromFile(kappafilename(i), grid%nlambda, xArray, &
                  grid%onekappaAbs(i,1:grid%nlambda), grid%onekappaSca(i,1:grid%nLambda), gfac)
          endif
       enddo
       if (writeoutput) then
          open(20, file="albedo.dat", status="unknown", form="formatted")
          write(20,'(a120)') "# Columns are: wavelength (microns), kappa ext (cm^2 g^-1), &
               &  kappa abs (cm^2 g^-1), kappa sca (cm^2 g^-1), albedo"
          write(20,*) "# Note that the opacities are per gram of gas"
          do i = 1, nLambda
             kAbs = SUM(grid%oneKappaAbs(1:nDustType,i)*grainFrac(1:nDustType))/1.e10
             kSca = SUM(grid%oneKappaSca(1:nDustType,i)*grainFrac(1:nDustType))/1.e10
             write(20,'(1p,5e13.4)') xArray(i)*angstomicrons, kAbs+kSca, kAbs, kSca, kSca/(kAbs+kSca)
          enddo
          close(20)
       endif



       call writeInfo("Creating Rosseland opacity lu table",TRIVIAL)
       call createRossArray(grid)
       call writeInfo("Done.",TRIVIAL)
    endif


    if (mie .or. (grid%geometry == "ttauri" .and. ttau_disc_on)) then
       ! construct the mie phase matrix
       call writeInfo("Computing Mie phase grid...",TRIVIAL)
       
       miePhase(:,:,:)%gfac = 0.d0

       if (isotropicScattering) then
          call writeInfo("Using isotropic scattering",FORINFO)
          miePhase = fillIsotropic()
          call writeInfo("Completed.",TRIVIAL)
          goto 666
       else if (henyeyGreensteinPhaseFunction) then
          call writeInfo("Using Henyey-Greenstein scattering",FORINFO)
          if (inputgFac /= 0.d0) gfac = inputgFac
          i = 1
          do j = 1, grid%nLambda
             do k = 1, nMumie
                mu = 2.d0*dble(k-1)/dble(nMumie-1)-1.d0
                miePhase(i,j,k) = fillHenyey(mu, gfac(j))
             enddo
          end do

!          cosang = 0.d0
!          do i = 1, 100000
!             uHat = randomUnitVector()
!             vec_tmp = uHat
!             uNew = newDirectionMie(vec_tmp, real(xArray(1)), xArray, nLambda, miePhase, nD!ustType, nMuMie, dble(grainFrac))
!             cosang = cosang + (uHat.dot.uNew)
!          enddo
!          write(*,*) "Average cos angle test ",cosang/100000.d0


          call writeInfo("Completed.",TRIVIAL)
          
       else

       if (readMiePhase.and.(nLambda>1)) then
          write(miefile,'(a,i3.3,a)') "miephasefile",nLambda,".dat"
          open(144, file=miefile, status="old", form="unformatted")
          read(unit=144) miePhase
          close(144)

       else

          allocate(mReal(1:nDusttype,1:nLambda))
          allocate(mImg(1:nDustType,1:nLambda))
          allocate(tmReal(1:nLambda))
          allocate(tmImg(1:nLambda))

          ! Set up refractive indices
          do i = 1, nDustType
             
             if (dustfile(i)) then
                do j = 1, grid%nLambda
                   do k = 1, nMumie
                      mu = 2.d0*dble(k-1)/dble(nMumie-1)-1.d0
                      miePhase(i,j,k) = fillHenyey(mu, gfac(j))
                   enddo
                end do
             else

                call parseGrainType(graintype(i), ngrain, grainname, x_grain)
                ! quick test for zero total dust abundance.
                total_dust_abundance = SUM(x_grain)
                if (total_dust_abundance <= 0.0) then
                   write(*,*) "Error:: total_dust_abundance <= 0.0 in torusMain."
                   write(*,*) "  ==> You probably forgot to assign dust abundance in your parameter file!"
                   write(*,*) "  ==> Exiting the program ... "
                   stop 
                end if
                
                ! allocate mem for temp arrays
                allocate(mReal2D(1:ngrain, 1:nLambda))
                allocate(mImg2D(1:ngrain, 1:nLambda))
                ! initializing the values
                mReal2D(:,:) = 0.0; mImg2D(:,:) = 0.0
                
                ! Find the index of refractions for all types of grains available
                do k = 1, ngrain
                   call getRefractiveIndex(xArray, nLambda, grainname(k), tmReal, tmImg, porousFillingFactor(k))
                   mReal2D(k,:) = tmReal(:)  ! copying the values to a 2D maxtrix
                   mImg2D(k,:)  = tmImg(:)   ! copying the values to a 2D maxtrix            
                end do
                
                ! Finding the weighted average of the refractive index.
                mReal(i,:) = 0.0; mImg(i,:) = 0.0
                do j = 1, nLambda
                   do k = 1, ngrain
                      mReal(i,j) = mReal(i,j) + mReal2D(k,j)*X_grain(k)
                      mImg(i,j)  = mImg(i,j) + mImg2D(k,j) *X_grain(k)
                   end do
                   mReal(i,j) = mReal(i,j) / total_dust_abundance
                   mImg(i,j)  = mImg(i,j)  / total_dust_abundance
                end do
                
                deallocate(mReal2D)
                deallocate(mImg2D)
             endif
          enddo
                ! You should use the new wrapper as it is much faster.
!                ! The old and new methods give exactly the same result when optimisations
!                ! are turned off.
!                ! When optimisations are turned on, the methods give different result,
!                ! and *neither* matches the result obtained when optimisations are off.
!                useOldMiePhaseCalc = .false.
!             
!                if (useOldMiePhaseCalc) then
!                      
!                      
!                      ilam_beg = 1
!                      ilam_end = grid%nLambda
!#ifdef MPI
!                      ! Set the range of index for a photon loop used later.     
!                      np = nThreadsGlobal
!                      n_rmdr = MOD(grid%nLambda,np)
!                      m = grid%nLambda/np
!                      
!                      if (myRankWorldGlobal .lt. n_rmdr ) then
!                         ilam_beg = (m+1)*myRankWorldGlobal + 1
!                         ilam_end = ilam_beg + m
!                      else
!                         ilam_beg = m*myRankWorldGlobal + 1 + n_rmdr
!                         ilam_end = ilam_beg + m -1
!                      end if
!#endif
!                      
!                      do j = ilam_beg, ilam_end
!                         do k = 1, nMumie
!                            mu = 2.*real(k-1)/real(nMumie-1)-1.
!                            call mieDistPhaseMatrixOld(aMin(i), aMax(i), a0(i), qDist(i), pDist(i), &
!                                 xArray(j), real(mu), miePhase(i,j,k), mReal(i,j), mImg(i,j))
!                         enddo
!                         call normalizeMiePhase(miePhase(i,j,1:nMuMie), nMuMie)
!                      end do
!#ifdef MPI                
!                      allocate(temp(1:grid%nlambda,1:nMuMie,1:4,1:4))
!                      temp = 0.
!!                      do j = 1, iLam_beg, iLam_end
 !                        do k = 1, nMuMie
 !                           do i1 = 1, 4
 !                              do i2 = 1, 4
 !                                 ! temp is real so need to explicitly convert from double precision to avoid compiler warning
 !                                 temp(j,k,i1,i2)= real(miePhase(i,j,k)%element(i1,i2))
 !                              enddo
 !                           enddo
 !                        enddo
 !                     enddo
 !                     allocate(tempArray(1:(grid%nLambda*nMuMie*4*4)))
 !                     allocate(tempArray2(1:(grid%nLambda*nMuMie*4*4)))
 !                     tempArray = reshape(temp, shape(tempArray))
                      
!                      call MPI_ALLREDUCE(tempArray, tempArray2, grid%nLambda, MPI_REAL,&
!                           MPI_SUM, MPI_COMM_WORLD, ierr)
!                      temp = reshape(tempArray2, shape(temp))
!                      do j = 1, grid%nLambda
!                         do k = 1, nMuMie
!                            do i1 = 1, 4
!                            do i2 = 1, 4
!                               miePhase(i,j,k)%element(i1,i2) = temp(j,k,i1,i2)
!                            enddo
!                         enddo
!                      enddo
!                   enddo
!                   deallocate(temp, temparray, temparray2)
!#endif
!                else
          call mieDistPhaseMatrixWrapper(nDustType, nLambda, nMuMie, xArray, mReal, mImg, miePhase)
!                end if
       endif

!       deallocate(mReal)
!       deallocate(mImg)
!       deallocate(tmReal)
!       deallocate(tmImg)
       !    endif
    endif
    if (writeMiePhase.and.writeoutput.and.(nlambda>1)) then
       write(miefile,'(a,i3.3,a)') "miephasefile",nLambda,".dat"
       open(144, file=mieFile, status="replace", form="unformatted")
       write(unit=144) miePhase
       close(144)
    end if
    
    
    call fixMiePhase(miePhase, nDustType, nLambda, nMuMie)
    do i = 1, nDustType
       do j = 1, nLambda
          call normalizeMiePhase(miePhase(i,j,1:nMuMie), nMuMie)
       enddo
    enddo
 
    
666 continue
    
    call returnKappa(grid, grid%OctreeRoot, 1, reset_kappa=.true.)
    call writeInfo("Completed.",TRIVIAL)
 endif
end subroutine createDustCrossSectionPhaseMatrix


  recursive subroutine allocateMemoryForDust(thisOctal)
    use inputs_mod, only : nDustType, grainFrac
    use octal_mod, only: allocateAttribute
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child 
    integer :: subcell, i
    
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call allocateMemoryForDust(child)
                exit
             end if
          end do
       else
          if (.not.associated(thisOctal%dustTypeFraction)) then
             call allocateAttribute(thisOctal%oldFrac, thisOctal%maxChildren)
             call allocateAttribute(thisOctal%dustType, thisOctal%maxChildren)
             call allocateAttribute(thisOctal%dustTypeFraction, thisOctal%maxChildren, nDustType)
             thisOctal%dustTypeFraction(subcell,1:nDustType) = grainfrac(1:nDusttype)
             thisOctal%oldfrac(subcell) = grainFrac(1)
          endif
       endif
    enddo
  end subroutine allocateMemoryForDust

  real function getCompositeGrainDensity(grainstring) result(density)
    implicit none
    integer :: nTypes, i
    character(len=80) :: name(100), grainstring
    real :: abundance(100)
    call parseGrainType(grainString, nTypes, name, abundance)
    density = 0.
    do i = 1, nTypes
       density = density + abundance(i)*getGrainDensity(name(i))
    enddo
  end function getCompositeGrainDensity
  
  real function getGrainDensity(graintype) result(density)
    character(len=*) :: graintype
    
    select case(grainType)
    case("am_olivine", "am_pyroxene")
       density = 3.71
    case("forsterite")
       density = 3.33
    case("enstatite")
       density = 2.8
    case("sio2")
       density = 2.21
    case("sil_dl")
       density = 3.6
    case("draine_sil")
       density = 3.5
    case("amc_zb")
       density = 2.0
    case("pinteISM")
       density = 0.5
    case DEFAULT
       call writeFatal("Unknown grain type in getGrainDensity: "//trim(graintype))
    end select
  end function getGrainDensity
  
real function getMeanMass2(porousFillingFactor, aMin, aMax, a0, qDist, pDist, graintype, grainDensity)  

  use constants_mod
  use mieDistCrossSection_mod, only: PowerInt

  implicit none
  real, intent(in) :: aMin, aMax, a0, qDist, pDist, porousFillingFactor
  real :: a1, a2, vol, fac
  integer :: i
  integer, parameter :: n = 1000
  real :: a(n)     ! grain sizes (log spaced)
  real :: f(n)     ! distribution function (normalized)
  real :: mass(n)  ! 
  real :: normFac
  character(len=*) :: grainType
  real :: grainDensity

  grainDensity = getGrainDensity(grainType)

  grainDensity = grainDensity * (1. - porousFillingFactor)

  if (aMin == aMax) then
     vol = real((4./3.)* pi * (aMin*microntocm)**3)
     getMeanMass2 = vol * graindensity 
  else
     a1 = log(aMin)
     a2 = log(aMax)
     !
     ! setting up the grain sizes
     do i = 1, n
        a(i) = (a1 + (a2 - a1) * real(i-1)/real(n-1))
        a(i) = exp(a(i))
        f(i) = a(i)**(-qDist) * exp(-(a(i)/a0)**pDist) 
     end do
     
     !
     ! normalize the dist function
     call PowerInt(n, 1, n, a, f, normFac)
     f(:) = f(:)/normFac
     
     !
     ! Finding the mean mass now.
     do i = 1, n
        vol = real((4./3.)* pi * (a(i)*microntocm)**3)
        mass(i) = vol * graindensity * f(i)    ! weighted by dist function
     end do
     
     call PowerInt(n, 1, n, a, mass, fac)
     
     getMeanMass2 = fac
  endif

  !  !
  !  !  For debug
  !  !
  !  getMeanMass2 = getMeanMass2*1000.0
  
end function getMeanMass2

subroutine readLambdaFile(lamFilename, lamArray, nLambda)
  character(len=*) :: lamFilename
  real(double) :: lamArray(:),junk
  integer :: nLambda, i


  call writeInfo("Reading wavelength points from file.", TRIVIAL)
  open(77, file=lamfilename, status="old", form="formatted")
  nLambda = 1
  ! Count the number of entries
333 continue
  read(77,*,end=334) junk
  nLambda = nLambda + 1
  goto 333
334 continue
  nlambda = nlambda - 1
  ! Rewind the file and read them in
  rewind(77)
  do i = 1, nLambda
     read(77,*) lamArray(i)
     lamArray(i) = lamArray(i) * 1e4 ! microns to angstrom
  enddo
  close(77)

end subroutine readLambdaFile

real(double) function dustToGasNumber(thisoctal, subcell, porousFillingFactor, aMin, aMax, &
     a0, qDist, pDist, graintype, grainDensity, nbins, binID)

  use inputs_mod, only : photoionPhysics, grainfrac
  use ion_mod, only : nGlobalIon, globalIonArray, returnMu

  type(octal) :: thisoctal
  integer :: subcell
  real(double) :: mu  !mean gas mass
  real(double) :: grainMass !mass of a single grain in the bin of interest
  real(double) :: massfrac  !fraction of the total dust mass in this bin

  !dust distribution properties
  real(double) :: porousFillingFactor, amin, amax, a0, qdist, pdist
  character(len=*) :: grainType
  real(double) :: grainDensity  

  integer :: binID, nbins
  

  !return the fraction of the dust to gas ratio in the bin of interest, as well
  !as the mass of a grain in that bin

  call getBinMassFraction(porousFillingFactor, aMin, aMax, a0, qDist, pDist, graintype, &
       grainDensity, nbins, binID, massfrac, grainMass)

  !mean gas mass
  if (photoionPhysics) then
     mu = returnMu(thisOctal, subcell, globalIonArray, nGlobalIon)
  else
     mu = 1.d0
  endif


  dustToGasNumber = (massFrac*grainfrac(1))*(mu*mHydrogen/grainmass)

  ! Md = total dust mass, md = mass of single grain
  !
  ! dust to gas number ratio = eta = nd / ng
  ! dust to gas mass ratio = delta = Md / Mg = eta * md/mu/mHydrogen
  ! therefore eta = delta*mu*mH/md
  !
end function dustToGasNumber



subroutine getBinMassFraction(porousFillingFactor, aMin, aMax, a0, qDist, pDist, graintype, &
     grainDensity, nbins, binID, massfrac, thisMass)  

  use constants_mod
  use mieDistCrossSection_mod, only: PowerInt

  implicit none
  real(double), intent(in) :: aMin, aMax, a0, qDist, pDist, porousFillingFactor
  real(double) :: a1, a2, vol
  integer, intent(in) :: nBins
  integer :: i, binID
!  integer, parameter :: n = 1000
  real(double) :: a(nBins)     ! grain sizes (log spaced)
  real(double) :: f(nBins)     ! distribution function (normalized)
  real(double) :: mass(nBins)  ! 
  character(len=*) :: grainType
  real(double) :: grainDensity
  real(double) :: density, massfrac, thismass
  real :: normfrac
  
  select case(grainType)
  case("am_olivine", "am_pyroxene")
     density = 3.71
  case("forsterite")
     density = 3.33
  case("enstatite")
     density = 2.8
  case("sio2")
     density = 2.21
  case("sil_dl")
     density = 3.6
  case("draine_sil")
     density = 3.5
  case("pinteISM")
     density = 0.5
  case DEFAULT
     density = grainDensity
  end select

  grainDensity = grainDensity * (1. - porousFillingFactor)

  if (aMin == aMax) then
     vol = real((4./3.)* pi * (aMin*microntocm)**3)
     Massfrac = 1.d0
     thisMass = vol * density 
  else
     a1 = log(aMin)
     a2 = log(aMax)
     !
     ! setting up the grain sizes
     do i = 1, nBins
        a(i) = (a1 + (a2 - a1) * real(i-1)/real(nBins-1))
        a(i) = exp(a(i))
        f(i) = a(i)**(-qDist) * exp(-(a(i)/a0)**pDist) 
     end do
     
     !
     ! normalize the dist function
     call PowerInt(nBins, 1, nBins, real(a), real(f), normFrac)
     f(:) = f(:)/normFrac
     
     !
     ! Finding the mean mass now.
     do i = 1, nBins
        vol = real((4./3.)* pi * (a(i)*microntocm)**3)
        mass(i) = vol * density * f(i)    ! weighted by dist function
     end do

     Massfrac = mass(binID)/sum(mass)
     thisMass = mass(binID)
    
  endif

end subroutine getBinMassFraction

real function getMedianSize(aMin, aMax, a0, qDist, pDist) result(aMedian)
  real, intent(in) :: aMin, aMax, a0, qDist, pDist
  real :: a1, a2
  integer :: i
  integer, parameter :: n = 1000
  real :: a(n)     ! grain sizes 
  real :: f(n)     
  real :: g(n)
  character(len=80) :: message

  if (aMin == aMax) then
     aMedian = aMin
  else
     a1 = log(aMin)
     a2 = log(aMax)
     !
     ! setting up the grain sizes
     do i = 1, n
        a(i) = (a1 + (a2 - a1) * real(i-1)/real(n-1))
        a(i) = exp(a(i))
        f(i) = a(i)**(-qDist) * exp(-(a(i)/a0)**pDist) 
     end do
     
     !
     ! find the median (in microns)
     g(:) = f(:)
     do i = 2, n
        g(i) = g(i) + g(i-1)
     enddo
! If a0 is small compared to the grain size the exponential term underflows and we cannot normalise g
     if (g(n) == 0.0) then
        write(message,'(3(a,f10.5),a,f4.1)') "aMin=", aMin, " aMax=", aMax, " a0=", a0, " pDist=", pDist
        call writeFatal("dust_mod::getMedianSize: cannot normalise grain size distribution "//message)
     end if
     g(:) = g(:)/g(n)
     call locate(g, n, 0.5, i)
     aMedian = a(i+1)
  endif
end function getMedianSize
  subroutine fillDustSettled(grid)
    type(GRIDTYPE) :: grid
    real(double) :: gasMass, dustMass(1:10)
    call writeInfo("Doing dust settling calculation...",TRIVIAL)
    call fillDustSettledRecur(grid, grid%octreeRoot)
    gasMass = 0.d0
    dustMass = 0.d0
    call sumDustMass(grid%octreeRoot, gasMass, dustMass)
    call normDustMass(grid%octreeRoot, gasMass, dustMass)
!    call testDustScaleheight(grid, grid%octreeRoot)
    call writeInfo("Done.",TRIVIAL)
  end subroutine fillDustSettled

  recursive subroutine fillDustSettledRecur(grid, thisOctal)

    use inputs_mod, only : rinner, router, nDustType, fracdustHeight,  rsublimation
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: rvec
    real(double) :: fitheight
    integer :: subcell, i, j, k
    real(double), allocatable :: zAxis(:), rho(:), subcellsize(:), normrho(:)
    real, allocatable :: temperature(:)
    integer :: nz
    integer, parameter :: m = 10000
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call fillDustSettledRecur(grid, child)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          if ((rVec%x > rinner).and.(rVec%x< rOuter)) then
             allocate(zAxis(m), rho(m), temperature(m), subcellsize(m), normrho(m))
             call getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, real(rVec%x), 0., nz, 1.)
             zAxis(1:nz) = zAxis(1:nz) / 1.d10
             zAxis(1:nz) = zAxis(1:nz)**2
             normrho(1:nz) = log(rho(1:nz)/rho(1))
             j = 1
             do while ((normrho(j) > -8.d0).and.(.not.(j > nz)))
                j = j + 1
             enddo
             nz  = j - 1
             
!             if (writeoutput) then
!                do j = 1,nz
!                   write(*,*) j,zaxis(j),normrho(j)
!                enddo
!                write(*,*) " "
!             endif
             
             if (all(normrho(1:nz) == 0.d0).or.(nz <= 2)) then
                fitHeight = 1.d30
             else
                call getLocalScaleheight(zAxis, normrho, nz, fitheight)
             endif


!             if (writeoutput) write(*,*) "expected height ",height*(rVec%x/(100.d0*autocm/1.d10))**betadisc, " fit ", fitheight

             do k = 1, nDustType

                if (abs(fracDustHeight(k)-1.d0)<1.d-3) then
                   thisOctal%dustTypeFraction(subcell,k) = 1.d0
                else
                   thisOctal%dustTypeFraction(subcell,k) = 1.d-30
                   if (thisOctal%rho(subcell) > 1d-30) then
                      thisOctal%dustTypeFraction(subcell, k) =  &
                           exp(-0.5d0 * (rVec%z**2)/((fracdustheight(k) * fitheight)**2))*rho(1)/thisOctal%rho(subcell)
                   endif
                endif

             enddo
             if (rVec%x < rSublimation) thisOctal%dustTypeFraction(subcell,:) = 1.d-30

!             thisOctal%dustTypeFraction(subcell, 2)  = grainfrac(2)
!             thisOctal%origDustTypeFraction(subcell,1:2) = thisOctal%dustTypeFraction(subcell,1:2)
             deallocate( zAxis, rho, temperature, subcellsize, normrho)

          endif
       end if
    end do

  end subroutine fillDustSettledRecur

  recursive subroutine testDustScaleheight(grid, thisOctal)

    use inputs_mod, only : rinner, router, nDustType
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    type(VECTOR) :: rvec
    integer :: subcell, i, j, k
    real(double), allocatable :: zAxis(:), rho(:), subcellsize(:), normrho(:)
    real, allocatable :: temperature(:)
    real(double) :: dustHeight(10), gasheight
    integer :: nz
    integer, parameter :: m = 10000
    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call testDustScaleheight(grid, child)
                exit
             end if
          end do
       else
          rVec = subcellCentre(thisOctal,subcell)
          if ((rVec%x > rinner).and.(rVec%x< rOuter)) then

             allocate(zAxis(m), rho(m), temperature(m), subcellsize(m), normrho(m))


             call getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, real(rVec%x), 0., nz, 1.)
             zAxis(1:nz) = zAxis(1:nz) / 1.d10
             zAxis(1:nz) = zAxis(1:nz)**2
             normrho(1:nz) = log(rho(1:nz)/rho(1))
             j = 1
             do while ((normrho(j) > -8.d0).and.(.not.(j > nz)))
                j = j + 1
             enddo
             nz  = j - 1
             call getLocalScaleheight(zAxis, normrho, nz, gasheight)

             do k = 1, nDustType
                call getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, real(rVec%x), 0., nz, 1., dustType=k)
                zAxis(1:nz) = zAxis(1:nz) / 1.d10
                zAxis(1:nz) = zAxis(1:nz)**2
                normrho(1:nz) = log(rho(1:nz)/rho(1))
                j = 1
                do while ((normrho(j) > -8.d0).and.(.not.(j > nz)))
                   j = j + 1
                enddo
                nz  = j - 1
                call getLocalScaleheight(zAxis, normrho, nz, dustheight(k))

             enddo
             if (writeoutput) write(*,*) "scaleheights ",dustheight(1:nDustType)/gasHeight


             deallocate( zAxis, rho, temperature, subcellsize, normrho)

          endif
       end if
    end do

  end subroutine testDustScaleheight

  recursive subroutine sumDustMass(thisOctal, gasMass, dustMass)

    use inputs_mod, only :  nDustType
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: dustMass(:), gasMass
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call sumDustMass(child, gasMass, dustMass)
                exit
             end if
          end do
       else

          dustMass(1:nDustType) = dustMass(1:nDustType) + &
               cellVolume(thisOctal,subcell) * thisOctal%dustTypeFraction(subcell,1:nDustType) * thisOctal%rho(subcell) * 1.d30
          gasMass = gasMass + cellVolume(thisOctal,subcell) * thisOctal%rho(subcell) * 1.d30
       end if
    end do

  end subroutine sumDustMass

  recursive subroutine normDustMass(thisOctal, gasMass, dustMass)
 
    use inputs_mod, only :  grainFrac, nDustType
    type(octal), pointer   :: thisOctal
    type(octal), pointer  :: child
    real(double) :: dustMass(:), gasMass, scaleFac(10)
    integer :: subcell, i

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call normDustMass(child, gasMass, dustMass)
                exit
             end if
          end do
       else

          scaleFac(1:nDustType) = GrainFrac(1:nDustType) / (dustMass(1:nDustType) / gasMass)

          thisOctal%dustTypeFraction(subcell,1:nDustType) = & 
               thisOctal%dustTypeFraction(subcell,1:nDustType) * scaleFac(1:nDustType)

          if (.not.associated(thisOctal%origDustTypeFraction)) &
               allocate(thisOctal%origDustTypeFraction(1:thisOctal%maxChildren, 1:nDustType))
          thisOctal%origDustTypeFraction(subcell,1:nDustType) = thisOctal%dustTypeFraction(subcell,1:nDustType)
       end if
    end do

  end subroutine normDustMass

  



  subroutine getLocalScaleheight(z, rho, nz, height)
    use utils_mod, only: linfit 

    use utils_mod, only : locate
    real(double) :: z(:), rho(:), height
    integer :: nz
    real(double) :: a, sigmaa, b, sigmab, rcoeff

    call LINFIT(z,rho,rho,nz, 0, A, SIGMAA, B, SIGMAB, Rcoeff)
    if (b < 0.d0) then
       height = sqrt(-1.d0/(2.d0*b))
    else
       height = 1.d30
    endif


  end subroutine getLocalScaleheight



  subroutine getTemperatureDensityRun(grid, zAxis, subcellsize, rho, temperature, xPos, yPos, nz, direction, dusttype)
    use amr_mod, only: amrGridValues
    use parallel_mod, only: torus_abort
    integer, optional :: dustType
    type(GRIDTYPE) :: grid
    type(octal), pointer   :: thisOctal
    integer, intent(out) :: nz
    real(double) :: rho(:)
    real :: temperature(:)
    real(double) :: zAxis(:), subcellsize(:)
    real :: xPos, yPos
    integer :: subcell
    real(double) :: rhotemp
    real :: temptemp
    real :: direction
    type(VECTOR) :: currentPos, temp
    real :: halfSmallestSubcell
    integer :: nzMax

    nzMax = SIZE(temperature)
    nz = 0
    halfSmallestSubcell = real(grid%halfSmallestSubcell)

    currentPos = VECTOR(xPos, yPos, direction*halfSmallestSubcell)

    do while(abs(currentPos%z) < grid%ocTreeRoot%subcellsize)
       call amrGridValues(grid%octreeRoot, currentPos, foundOctal=thisOctal, &
            foundSubcell=subcell, rho=rhotemp, temperature=temptemp)
       thisOctal%chiLine(subcell) = 1.e-30
!       if (thisOctal%inFlow(subcell)) then
          nz = nz + 1
          if (nz>nzmax) then
             call torus_abort("nz>nzMax in getTemperatureDensityRun. Aborting ...")
          endif
          temperature(nz) = temptemp
          rho(nz) = rhotemp
          if (present(dustType)) rho(nz) = thisOctal%rho(subcell) * thisOctal%dustTypeFraction(subcell,dustType)
          temp = subCellCentre(thisOctal, subcell)
          zAxis(nz) = temp%z
          subcellsize(nz) = thisOctal%subcellsize
!       endif
          currentPos = VECTOR(xPos, yPos, zAxis(nz)+0.5*direction*thisOctal%subcellsize+direction*halfSmallestSubcell)
!       else
!          currentPos = VECTOR(xPos, yPos, grid%octreeRoot%subcellsize+halfSmallestSubcell)
!       endif
    end do
    zAxis(1:nz) = abs(zAxis(1:nz)) * 1.d10  ! convert to cm
  end subroutine getTemperatureDensityRun
 
end module dust_mod

