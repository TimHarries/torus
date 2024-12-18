module pah_mod

  use citations_mod
  use constants_mod
  use utils_mod, only : locate
  use unix_mod
  use messages_mod
  use random_mod
  implicit none

  type PAHTABLETYPE
     integer :: nu
     real(double), allocatable :: u(:)        ! dimensionless scale fac
     real(double), allocatable :: adot(:)     ! erg cm-3 s-1 / rho
     real(double), allocatable :: freq(:)     ! Hz
     real(double), allocatable :: jnu(:,:)    ! erg cm-2 s-1 sr-1 H-1
     real(double), allocatable :: lamKappa(:) ! micron
     real(double), allocatable :: kappaSca(:) ! cm2 (g gas)-1
     real(double), allocatable :: kappaAbs(:) ! cm2 (g gas)-1
     real(double), allocatable :: gFac(:)
     real(double), allocatable :: pnu(:,:)
  end type PAHTABLETYPE

  real(double) :: a01, a02, sigma1, sigma2, b1, b2
  type(PAHTABLETYPE) :: PAHtable

contains

  real(double) function dnda(a, amin)
    real(double) :: a, amin, n01, n02, x1, x2
    real(double) :: am1, am2
    real(double), parameter :: mc = 12.d0 * mHydrogen
    real(double), parameter :: grainDensity = 3.5d0

    dnda = 0.d0

    am1 = a01 * exp(3.d0 * sigma1**2)
    am2 = a02 * exp(3.d0 * sigma2**2)

    write(*,*) "am1 am2 ",am1,am2
    x1 = log(am1 / amin)/(sqrt(2.d0)*sigma1)
    x2 = log(am2 / amin)/(sqrt(2.d0)*sigma2)

    write(*,*) " x1 x2 ",x1 ,x2, erf(x1),erf(x2)
    if (erf(x1) > -1.d0) then
       n01 = (3.d0)/(twoPi**1.5d0) * exp(4.5d0*sigma1**2) * mc / (1.d0 + erf(x1)) / (grainDensity * am1**3 * sigma1) * b1
       n02 = (3.d0)/(twoPi**1.5d0) * exp(4.5d0*sigma2**2) * mc / (1.d0 + erf(x2)) / (grainDensity * am2**3 * sigma2) * b2
       write(*,*) " n01, n02 ", n01, n02
       dnda = n01/a * exp(- (log(a/a01)**2)/(2.d0*sigma1**2))

       write(*,*) "dnda ", dnda
       dnda = dnda + n02/a * exp(- (log(a/a02)**2)/(2.d0*sigma2**2))
       write(*,*) "dnda ", dnda
    endif
  end function dnda

  subroutine readDraineOpacity()
    use inputs_mod, only : PAHscale
    character(len=80) :: filename, dataDirectory, cjunk
    integer :: i, n
    real(double) :: junk0, junk1, junk2, junk3, junk4, junk5, junk6
    real(double), parameter :: mu = 1.4d0

    call unixGetenv("TORUS_DATA", dataDirectory, i)
    filename = trim(dataDirectory)//"/comp_opacity_abs.out"

    n = 751
    allocate(PAHtable%kappaAbs(1:n), PAHtable%kappaSca(1:n), PAHtable%lamKappa(1:n))


    open(30, file=filename, status="old", form="formatted")
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    do i = 1,n
       read(30,*) PAHtable%lamKappa(i), junk1, junk2, junk3, junk4, junk5, junk6

       PAHtable%kappaAbs(i) = junk1 + junk2 + junk4 + junk5
    enddo
    close(30)

    filename = trim(dataDirectory)//"/comp_opacity_ext.out"

    open(30, file=filename, status="old", form="formatted")
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk


    do i = 1,n
       read(30,*) junk0, junk1, junk2, junk3, junk4, junk5, junk6
       PAHtable%kappaSca(i) =  (junk1 + junk2 + junk4 + junk5) - PAHtable%kappaAbs(i)
    enddo
    close(30)
    ! table reads xsec (um2/H) - convert to kappa (cm2/g_gas)
    PAHtable%kappaAbs = PAHscale * PAHtable%kappaAbs / (mu*mHydrogen) * 1.d-8
    PAHtable%kappaSca = PAHscale * PAHtable%kappaSca / (mu*mHydrogen) * 1.d-8

   if (writeoutput) then
      open(36, file='albedo_pah.dat', status="unknown", form="formatted")
      write(36,'(a)') "# Columns are: wavelength (microns), kappa ext (cm^2 g^-1), kappa abs (cm^2 g^-1), kappa sca (cm^2 g^-1)"
      write(36,*) "# Note that the opacities are per gram of gas"
      do i = 1, SIZE(PAHtable%lamKappa)
         write(36,*) PAHtable%lamKappa(i), (PAHtable%kappaAbs(i) + PAHtable%kappaSca(i)), &
         PAHtable%kappaAbs(i), PAHtable%kappaSca(i)
      enddo
      close(36)
   endif

  end subroutine readDraineOpacity

  subroutine readHensleyOpacity()
    use inputs_mod, only : PAHscale
    character(len=200) :: filename, dataDirectory
    character(len=80) :: cjunk
    integer :: i, n
    real(double) :: junk0, junk1, junk2, junk3, junk4, junk5, junk6, junk
    real(double), parameter :: mu = 1.4d0
    ! Hensley & Draine 2023 PAH opacities

    call unixGetenv("TORUS_DATA", dataDirectory, i)
    filename = trim(dataDirectory)//"/astrodust+PAH_opacities.dat"
    call addBibcode("2023ApJ...948...55H","Astrodust+PAH opacities")

    n = 1000
    allocate(PAHtable%kappaAbs(1:n), PAHtable%kappaSca(1:n), PAHtable%lamKappa(1:n))


    open(30, file=filename, status="old", form="formatted")
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    read(30,'(A)') cjunk
    ! AA: first entry in the file is entered by me to give a const Astrodust kappa in EUV
    ! But this is irrelevant for PAHs 
    read(30,*) junk, junk1, junk2, junk3, junk4, junk5, junk6, junk, junk

    ! read original data
    do i = 1,n
       read(30,*) PAHtable%lamKappa(i), junk1, junk2, junk3, junk4, junk5, junk6, junk, junk

       PAHtable%kappaAbs(i) = junk1
       PAHtable%kappaSca(i) = junk3
    enddo
    close(30)
    ! table reads xsec (cm2/H) - convert to kappa (cm2/g dust), then to per g of gas by setting PAHscale = desired torus d/g ratio
    PAHtable%kappaAbs = PAHscale * PAHtable%kappaAbs / mHydrogen / 0.00708d0
    PAHtable%kappaSca = PAHscale * PAHtable%kappaSca / mHydrogen / 0.00708d0

   if (writeoutput) then
      open(36, file='albedo_pah.dat', status="unknown", form="formatted")
      write(36,'(a)') "# Columns are: wavelength (microns), kappa ext (cm^2 g^-1), kappa abs (cm^2 g^-1), kappa sca (cm^2 g^-1)"
      write(36,*) "# Note that the opacities are per gram of gas"
      do i = 1, SIZE(PAHtable%lamKappa)
         write(36,*) PAHtable%lamKappa(i), (PAHtable%kappaAbs(i) + PAHtable%kappaSca(i)), &
         PAHtable%kappaAbs(i), PAHtable%kappaSca(i)
      enddo
      close(36)
   endif

  end subroutine readHensleyOpacity

  ! read DL07 table
  subroutine readPAHEmissivityTable()
    use inputs_mod, only : pahtype, pahscale, pahKappa
    integer :: i, j 
    character(len=10) :: cval
    character(len=120) :: filename, dataDirectory, cjunk
    real :: vjunk
    real(double) :: lambda

    call unixGetenv("TORUS_DATA", dataDirectory, i)

    PAHtable%nU = 30
    allocate(PAHtable%U(1:PAHtable%nU), PAHtable%freq(1001), PAHtable%jnu(1:PAHtable%nU,1:1001))
    PAHtable%U(01) = 0.10
    PAHtable%U(02) = 0.15
    PAHtable%U(03) = 0.20
    PAHtable%U(04) = 0.30
    PAHtable%U(05) = 0.40
    PAHtable%U(06) = 0.50
    PAHtable%U(07) = 0.70
    PAHtable%U(08) = 0.80
    PAHtable%U(09) = 1.00
    PAHtable%U(10) = 1.20
    PAHtable%U(11) = 1.50
    PAHtable%U(12) = 2.00
    PAHtable%U(13) = 2.50
    PAHtable%U(14) = 3.00
    PAHtable%U(15) = 4.00
    PAHtable%U(16) = 5.00
    PAHtable%U(17) = 7.00
    PAHtable%U(18) = 8.00
    PAHtable%U(19) = 12.0
    PAHtable%U(20) = 15.0
    PAHtable%U(21) = 20.0
    PAHtable%U(22) = 25.0
    PAHtable%U(23) = 1e2
    PAHtable%U(24) = 3e2
    PAHtable%U(25) = 1e3
    PAHtable%U(26) = 3e3
    PAHtable%U(27) = 1e4
    PAHtable%U(28) = 3e4
    PAHtable%U(29) = 1e5
    PAHtable%U(30) = 3e5
   
    do i = 1, PAHtable%nU

       if (PAHtable%u(i) <= 10.d0) then
          write(cval,'(f4.2)') PAHtable%u(i)          
       else if (PAHtable%u(i) <= 10.d0) then
          write(cval,'(f5.2)') PAHtable%u(i)
       else if ((PAHtable%u(i) > 10.d0).and.(PAHtable%u(i) < 100.d0)) then
          write(cval, '(f4.1)') PAHtable%u(i)
       else
          if (i == 23) cval="1e2"
          if (i == 24) cval="3e2"
          if (i == 25) cval="1e3"
          if (i == 26) cval="3e3"
          if (i == 27) cval="1e4"
          if (i == 28) cval="3e4"
          if (i == 29) cval="1e5"
          if (i == 30) cval="3e5"
       endif
       write(filename,'(a,a,a,a,a,a,a,a,a,a,a,a)') &
            trim(dataDirectory),"/PAH/", &
            "U",trim(cval),"/U",trim(cval),"_", &
            trim(cval),"_",trim(PAHtype),".txt"
       call writeInfo("Reading PAH emissivities from: "//trim(filename),TRIVIAL)

       open(20,file=filename,status="old",form="formatted")
       read(20,'(a)') cjunk
       read(20,'(a)') cjunk
       read(20,'(a)') cjunk
       read(20,*) a01, sigma1, b1
       read(20,*) a02, sigma2, b2
       a01 = a01 / microntocm
       a02 = a02 / microntocm
       do j = 1, 56
          read(20,'(a)') cjunk
       enddo
       do j = 1, 1001
          read(20,*) lambda, vjunk, PAHtable%jnu(i,j)

          ! jnu read in as Jy cm2 sr-1 H-1

          PAHtable%freq(j) = cSpeed/(lambda*micronToCm)
          PAHtable%jnu(i,j) = PAHscale * PAHtable%jnu(i,j) * 1.d-23 ! from jy to erg/cm^2/s
       enddo
       close(20)
    enddo



    call writeInfo("Reading PAH opacities for type: "//trim(pahKappa),TRIVIAL)
    select case(pahKappa)
       case("astrodust+PAH_opacities.dat")
          call readHensleyOpacity()
       case("DL07")
          call readDraineOpacity()
       case DEFAULT
          call writeWarning("Couldn't find PAH opacities, using DL07 instead")
          stop
          call readDraineOpacity()
    end select

    call createPAHprobs()
    call calculateAdots()
  end subroutine readPAHEmissivityTable

  ! read Draine et al. 2021 tables
  subroutine readPAHEmissivityTable2021()
    use inputs_mod, only : pahtype, pahscale, pahKappa
    integer :: i, j 
    integer, parameter :: nfreq = 1973
    character(len=10) :: cval
    character(len=120) :: filename, dataDirectory, cjunk, population, metallicity, populationAge
    real :: vjunk, ion, neutral
    real(double) :: lambda, logU

    call unixGetenv("TORUS_DATA", dataDirectory, i)
    write(population, '(a)') "bc03"    ! stellar population reference
    write(metallicity, '(a)') "z0.02"  ! Z = 0.02 == solar
    write(populationAge, '(a)') "3e6"  ! yr
    PAHtable%nU = 15

    allocate(PAHtable%U(1:PAHtable%nU), PAHtable%freq(nfreq), PAHtable%jnu(1:PAHtable%nU,1:nfreq))
    do i = 1, PAHtable%nU
       logU = dble(i-1)/2.d0
       PAHtable%U(i) = 10.d0**logU

       write(filename,'(a,a,a,a,a,a,f4.2,a)') &
            trim(dataDirectory),"/draine2021pah/", &
            "pahspec.out_", &
            population, "_", &
            metallicity, "_", &
            populationAge, "_" ,&
            logU, & 
            "_st_std"

       call writeInfo("Reading PAH emissivities from: "//trim(filename),TRIVIAL)

       open(20,file=filename,status="old",form="formatted")
       ! header
       do j = 1, 7
          read(20,'(a)') cjunk
       enddo
       ! variables not needed but DL07 routine uses them(?)
       a01 = 0.d0
       sigma1 = 0.d0
       b1 = 0.d0
       a02 = 0.d0
       sigma2 = 0.d0
       b2 = 0.d0
       ! data
       do j = 1, nfreq
          read(20,*) lambda, vjunk, vjunk, ion, neutral

          PAHtable%freq(j) = cSpeed/(lambda*micronToCm)  ! s-1

          ! nu*Pnu (erg/s/H) to jnu (erg/cm2/s)
          PAHtable%jnu(i,j) = dble(ion + neutral) / (PAHtable%freq(j) * 4.d0 * pi) * PAHscale
       enddo
       close(20)
    enddo



    call writeInfo("Reading PAH opacities for type: "//trim(pahKappa),TRIVIAL)
    select case(pahKappa)
       case("astrodust+PAH_opacities.dat")
          call readHensleyOpacity()
       case DEFAULT
          call writeFatal("Draine2021 PAH emission requires Hensley PAH opacity")
          stop
    end select

    call createPAHprobs()
    call calculateAdots()
  end subroutine readPAHEmissivityTable2021
 
  real(double) function getKappaAbsPAH(freq)
    use utils_mod
    real(double) :: freq, lambda
    integer :: i
    getKappaAbsPAH = tiny(getKappaAbsPAH)
    lambda = (cSpeed/freq)/microntocm

    if ((lambda >= PAHtable%lamKappa(1)).and.(lambda <= PAHtable%lamKappa(SIZE(PAHtable%lamKappa)))) then
       call locate(PAHtable%lamKappa, SIZE(PAHtable%lamKappa), lambda, i)
       getKappaAbsPAH = logint(lambda, PAHtable%lamKappa(i), PAHtable%lamKappa(i+1), PAHtable%kappaAbs(i), PAHtable%kappaAbs(i+1))
    endif
  end function getKappaAbsPAH

  real(double) function getKappaScaPAH(freq)
    use utils_mod
    real(double) :: freq, lambda
    integer :: i
    getKappaScaPAH = tiny(getKappaScaPAH)
    lambda = (cSpeed/freq)/microntocm

    if ((lambda >= PAHtable%lamKappa(1)).and.(lambda <= PAHtable%lamKappa(SIZE(PAHtable%lamKappa)))) then
       call locate(PAHtable%lamKappa, SIZE(PAHtable%lamKappa), lambda, i)
       getKappaScaPAH = logint(lambda, PAHtable%lamKappa(i), PAHtable%lamKappa(i+1), PAHtable%kappaSca(i), PAHtable%kappaSca(i+1))
    endif
  end function getKappaScaPAH


  real(double) function Jisrf(freq)
    use spectrum_mod
    character(len=80) :: dataDirectory, ifilename
    integer :: i
    real(double) :: freq
    real(double) :: lambda
    type(SPECTRUMTYPE),save :: spectrum
    logical, save :: firstTime = .true.
    logical :: ok
    


    if (firstTime) then
       firstTime = .false.

       call unixGetenv("TORUS_DATA", dataDirectory, i)
       ifilename = trim(dataDirectory)//"/"//"isrf.dat"
       
       call readSpectrum(spectrum, ifilename, ok)
       if (.not.ok) then
          call writeFatal("Cannot read isrf.dat")
          stop
       endif

       spectrum%lambda = spectrum%lambda * micronsToAngs ! A
       spectrum%dlambda = spectrum%dlambda * micronsToAngs ! A
       spectrum%flux = spectrum%flux / micronsToAngs ! per A
       
       spectrum%flux = spectrum%flux / (4.d0 * pi)  ! from fourpi Jlambda to Jlambda
    end if


    lambda = (cSpeed/freq) / AngstromToCm ! lambda in A
    ! Jlambda (per A) to Jnu (per Hz)
    Jisrf = getFlux(lambda, spectrum)*cspeed/freq**2 / angstromToCm
       
  end function Jisrf


  subroutine createPAHprobs()

    integer :: i,j 
    allocate(PAHtable%pnu(PAHtable%nu, 1001))

    PAHtable%pnu = 0.d0

    do i = 1, PAHtable%nu

       do j = 2, 1001
          PAHtable%pnu(i,j) = PAHtable%pnu(i,j-1) + PAHTable%jnu(i,j) * &
               (PAHtable%freq(j) - PAHtable%freq(j-1))
       enddo
       PAHtable%pnu(i,1:1001) = PAHtable%pnu(i,1:1001) / PAHtable%pnu(i,1001)
    enddo
  end subroutine createPAHprobs

  subroutine calculateAdots()
    integer :: i,j , n_nu
    real(double) :: g0(1:pahTable%nu), nu, nu_1, nu_n, dnu
    allocate(PAHtable%adot(PAHtable%nu))

    do i = 1, PAHtable%nu
       PAHtable%adot(i) = 0.d0
       do j = 2, 1001
          ! this is Adot/rho = 4pi sum(kappa J dnu); J = U * Jisrf (e.g. robitaille+ 2011)
          ! note MCRT Adot   = 4pi sum(kappa rho J dnu)
          ! cgs units 
          PAHtable%adot(i) = PAHtable%adot(i) + 4.d0 * pi * PAHtable%u(i) * Jisrf(PAHtable%freq(j)) &
               * getkappaAbsPAH(PAHtable%freq(j)) *  (PAHtable%freq(j) - PAHtable%freq(j-1))
       enddo
!       if (writeoutput) write(*,'(i6,a,es13.5,a,es13.5)') i, " U_i ", PAHtable%u(i), " adot_i ",PAHtable%adot(i)
    enddo

    ! test Jisrf is correct
    nu_1 = 6.d0*evtoErg/hcgs 
    nu_n = 13.6d0*evtoerg/hcgs
    n_nu = 100
    dnu = (nu_n - nu_1)/dble(n_nu)
    do i = 1, PAHtable%nu
       g0(i) = 0.d0
       do j = 2, n_nu
          nu = nu_1 + dble(j-1)*dnu
          if ((hcgs*nu > 6.d0*evtoerg) .and. (hcgs*nu < 13.6d0*evtoerg)) then
             g0(i) = g0(i) + 4.d0 * pi * PAHtable%u(i) * Jisrf(nu) * dnu
          endif
       enddo
       ! g0/habing/U ~ 1 is correct
       if (writeoutput) write(*,'(i6,a,es11.3,a,es11.3,2(a,es11.3))') i, " U_i ", PAHtable%u(i), " adot_i ",PAHtable%adot(i),&
          " g0_i/U_i ",g0(i)/habing/pahtable%u(i)!," u(6-13.6eV)_i ", g0(i)/cspeed
    enddo
  end subroutine calculateAdots



!  real(double) function PAHemissivity(lambda, u, rho)
!    real(double) :: lambda, u, thisjnu(1001)
!    real(double) :: freq, t1, rho, nH
!    integer :: i, j
!    freq = cspeed/(lambda * angstromTocm)
!
!
!    PAHemissivity = 0.d0
!    if ( (freq > PAHtable%freq(1)).and.(freq < PAHtable%freq(1001)) ) then
!       if ( (u > PAHtable%u(1)) .and. (u < PAHtable%u(PAHtable%nu)) ) then
!
!          nH = rho/mHydrogen
!
!          call locate(PAHtable%u, PAHtable%nu, u, i)
!          t1 = (u - PAHtable%u(i)) / (PAHtable%u(i+1) - PAHtable%u(i))
!          
!          thisJnu = PAHtable%jnu(i,:) + t1 * (PAHtable%jnu(i+1,:) + PAHtable%jnu(i,:))
!       
!          call locate(PAHtable%freq, 1001, freq, j)
!
!          t1 = (freq - PAHtable%freq(j)) / (PAHtable%freq(j+1) - PAHtable%freq(j))
!
!          PAHemissivity =  nh *  (thisJnu(j) + t1 * (thisJnu(j+1) - thisJnu(j)))
!       endif
!    endif
!  end function PAHemissivity

  real(double) function PAHemissivityFromAdot(lambda, adot, rho)
    real(double) :: lambda, adot, thisjnu(1001),thisAdot
    real(double) :: freq, t1, rho, nH, adotpermass
    integer :: i, j
    real(double), parameter :: mu = 1.4d0
    adotpermass = adot/rho

    freq = cspeed/(lambda * angstromTocm)

    ! input adot should be actual adot (erg/s/cm3)

    PAHemissivityfromAdot = 0.d0
    if ( (freq > PAHtable%freq(1)).and.(freq < PAHtable%freq(1001)) ) then
       if ( (adotpermass > PAHtable%adot(1)) .and. (adotpermass < PAHtable%adot(PAHtable%nu)) ) then

          nH = rho/(mu * mHydrogen)

          thisAdot = min(PAHtable%adot(PAHtable%nu),max(adotpermass, PAHtable%adot(1)))
          call locate(PAHtable%adot, PAHtable%nu, thisadot, i)

          t1 = (thisadot - PAHtable%adot(i)) / (PAHtable%adot(i+1) - PAHtable%adot(i))
          
          thisJnu = PAHtable%jnu(i,:) + t1 * (PAHtable%jnu(i+1,:) + PAHtable%jnu(i,:))
       
          call locate(PAHtable%freq, 1001, freq, j)

          t1 = (freq - PAHtable%freq(j)) / (PAHtable%freq(j+1) - PAHtable%freq(j))

          PAHemissivityfromAdot =  nh *  (thisJnu(j) + t1 * (thisJnu(j+1) - thisJnu(j)))
          ! jnu per Hz
       endif
    endif
  end function PAHemissivityFromAdot


  real(double) function  getPAHFreq(u)
    real(double) :: prob(1001), u 
    integer :: i, j
    real(double) :: t1, r

    call locate(PAHtable%u, PAHtable%nu, u, i)
    t1 = (u - PAHtable%u(i)) / (PAHtable%u(i+1) - PAHtable%u(i))

    prob = PAHtable%pnu(i,:) + t1 * (PAHtable%pnu(i+1,:) + PAHtable%pnu(i,:))
    call randomNumberGenerator(getDouble=r)

    call locate(prob, 1001, r, j)

    t1 = (r - prob(j))/(prob(j+1) - prob(j))

    getPAHfreq = PAHtable%freq(j) + t1 * (PAHtable%freq(j+1) - PAHtable%freq(j))
  end function getPAHFreq

  real(double) function  getPAHFreqfromAdot(adot,rho)
    real(double) :: prob(1001), adot, thisAdot,rho
    integer :: i, j
    real(double) :: t1, r, adotpermass
    adotpermass = adot/rho

    thisAdot = max(min(adotpermass, PAHtable%adot(PAHtable%nu)),PAHtable%adot(1))

    call locate(PAHtable%adot, PAHtable%nu, thisadot, i)
    t1 = (thisadot - PAHtable%adot(i)) / (PAHtable%adot(i+1) - PAHtable%adot(i))

    prob = PAHtable%pnu(i,:) + t1 * (PAHtable%pnu(i+1,:) + PAHtable%pnu(i,:))
    call randomNumberGenerator(getDouble=r)

    call locate(prob, 1001, r, j)

    t1 = (r - prob(j))/(prob(j+1) - prob(j))

    getPAHfreqfromAdot = PAHtable%freq(j) + t1 * (PAHtable%freq(j+1) - PAHtable%freq(j))
  end function getPAHFreqfromAdot


  subroutine testPAHtable(u,fname)
    integer, parameter :: n = 1000
    real(double) :: intensity(n), lambda(n)
    character(len=*) :: fname
    real(double) :: u, thisLambda
    integer :: i,j 

    do i = 1, n
       lambda(i) = 5.d0 + 19.d0*dble(i-1)/dble(n-1)
    enddo

    intensity = 0.d0

    do i = 1, 1000000
       thislambda = (cSpeed/getPAHfreq(u))/microntocm
       if ((thisLambda > lambda(1)).and.(thisLambda <= lambda(n))) then
          call locate(lambda, n, thisLambda, j)
          intensity(j) = intensity(j) + 1.d0
       endif
    enddo
    if (writeoutput) then
       open(43, file=fname, status="unknown", form="formatted")
       do i = 1, n
          write(43,*) lambda(i), intensity(i)
       enddo
       close(43)
    endif
          
  end subroutine testPAHtable

  ! not used?
  subroutine readPAHkappa()
    character(len=80) :: dataDirectory, ifilename, junk
    integer :: i, j, na, nLambda
    real(double) :: mass, dm, weight, totweight, aMin, aMax,  da
    integer :: istart, iend
    real(double) :: grainDensity
    real(double), allocatable :: a(:), lambda(:)
    real(double), allocatable :: qext(:,:), kappaSca(:), kappaAbs(:), kappaExt(:), gFac(:)
    real(double), allocatable :: qabs(:,:)
    real(double), allocatable :: qsca(:,:)
    real(double), allocatable :: g(:,:)

    call unixGetenv("TORUS_DATA", dataDirectory, i)
    ifilename = trim(dataDirectory)//"/PAH/PAHneu_30"

    open(33, file=ifilename, status="old", form="formatted")

    read(33,*) junk
    read(33,*) junk
    read(33,*) junk
    read(33,*) junk
    read(33,*) junk
    read(33,*) junk
    read(33,*) na
    read(33,*) nLambda
    if (writeoutput) write(*,*) "na ",na, " nalmbda ",nlambda

    allocate(a(1:na), lambda(1:nLambda))
    allocate(qext(1:na,1:nLambda))
    allocate(qabs(1:na,1:nLambda))
    allocate(qsca(1:na,1:nLambda))
    allocate(g(1:na,1:nLambda))



    allocate(PAHtable%kappaAbs(1:nLambda))
    allocate(PAHtable%kappaSca(1:nLambda))
    allocate(PAHtable%gFac(1:nLambda))
    allocate(PAHtable%lamKappa(1:nLambda))
    PAHtable%kappaAbs = 0.d0
    PAHtable%kappaSca = 0.d0
    PAHtable%gFac = 0.d0


    do i = 1, na
       read(33,*) a(i)
       read(33,*) junk
       do j = 1, nLambda
          read(33,'(a)') junk
          read(junk,'(4e10.3,e9.2)') lambda(nlambda+1-j), qext(i,nlambda+1-j), &
               qabs(i,nLambda+1-j), qsca(i,nLambda+1-j), g(i,nLambda+1-j)
       enddo
    enddo
    close (33)

    aMin = 1e-4
    amax = 1.e-2
    grainDensity = 3.5d0


    
    call locate(a, na, aMin, iStart)
    call locate(a, na, aMax, iEnd)

    mass = 0.d0
    totWeight = 0.d0

    allocate(kappaAbs(1:nLambda))
    allocate(kappaSca(1:nLambda))
    allocate(kappaExt(1:nLambda))
    allocate(gfac(1:nLambda))

    kappaAbs = 0.d0
    kappaSca = 0.d0
    kappaExt = 0.d0
    gfac = 0.d0

    do i = iStart+1,iEnd



       dm = (fourPi / 3.d0)*(a(i)*micronToCm)**3 * grainDensity - &
            (fourPi / 3.d0)*(a(i-1)*micronToCm)**3 * grainDensity

       mass = mass + dm
       da = a(i) - a(i-1)

       if (writeoutput) write(*,*) " a ",i,a(i)
       weight = dnda(a(i),amin) * da
       totWeight = totWeight + weight
      do j = 1, nLambda

         kappaExt(j) = kappaExt(j) + weight * qExt(i,j) * pi *(a(i)*micronToCm)**2 
         kappaAbs(j) = kappaAbs(j) + weight * qAbs(i,j) * pi *(a(i)*micronToCm)**2 
         kappaSca(j) = kappaSca(j) + weight * qSca(i,j) * pi *(a(i)*micronToCm)**2 
         gFac(j) = gFac(j) + weight * g(i,j) * pi *(a(i)*micronToCm)**2 
      enddo
   enddo
   kappaExt = kappaExt / totWeight / mass
   PAHtable%kappaSca = kappaSca / totWeight / mass
   PAHtable%kappaAbs = kappaAbs / totWeight / mass
   PAHtable%gFac = gFac / totWeight
   
   PAHtable%lamKappa = lambda


  end subroutine readPAHkappa
       

end module pah_mod
