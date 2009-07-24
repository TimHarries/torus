module gas_opacity_mod

use constants_mod
use utils_mod, only: locate
use unix_mod, only: unixGetenv

implicit none

type molecularKappaGrid
   character(len=80) :: filename
   character(len=80) :: source
   character(len=3) :: molecule
   real :: mu
   integer :: nLam
   integer :: nTemps
   real, pointer :: lamArray(:) => null()
   real, pointer :: tempArray(:) => null()
   real, pointer :: kapArray(:,:) => null()
end type molecularKappaGrid

type tsujiPPtable
   integer :: nAtoms
   integer :: nTemps
   integer :: nPressures
   real, pointer :: pp(:,:,:) => null()
   real, pointer :: temperature(:) => null()
   real, pointer :: logPressure(:) => null()
end type tsujiPPtable

type tsujiKPtable
   integer :: nMolecules 
   character(len=3), pointer :: molecule(:) => null()
   real, pointer :: a(:,:) => null()
end type tsujiKPtable


type(molecularKappaGrid), save :: TioLookuptable
type(molecularKappaGrid), save :: chLookuptable
type(molecularKappaGrid), save :: nhLookuptable
type(molecularKappaGrid), save :: ohLookuptable
type(molecularKappaGrid), save :: h2Lookuptable
type(molecularKappaGrid), save :: coLookuptable
type(molecularKappaGrid), save :: h2oLookuptable


type(tsujiPPtable), save :: tsujipplookuptable
type(tsujiKPtable), save :: tsujikplookuptable

contains

subroutine readTio(nLines,lambda,kappa,excitation,g)

  implicit none
  integer :: nLines
  real :: lambda(:),kappa(:),excitation(:),g(:)
  integer :: vUp, vLow, i
  integer, allocatable :: jLow(:)
  character(len=4) :: branch
  real, allocatable :: wavenumber(:),gf(:),alpha(:)
  real, parameter :: mHydrogen = 1.6733e-24      ! g
  real, parameter :: mElectron = 9.109565D-28    ! g
  real, parameter :: eCharge = 4.803242384E-10 !

  real, parameter :: cSpeed = 2.99792458e10  ! cm/s
  real, parameter :: pi = 3.141592654

  nLines = 11369860 


  allocate(gf(1:nLines))
  allocate(wavenumber(1:nLines))
  allocate(alpha(1:nLines))
  allocate(jlow(1:nLines))


  open(20,file="/home/th/plez/TiO.dat",status="old",form="formatted")

  write(*,*) "Reading all ",nLines," TiO transitions..."

  do i = 1, nLines

     if (mod(i,nLines/10) == 0) write(*,*) i 

!     read(20,'(4x,i2,i2,i3,1x,a4,e9.2,f10.3,f6.3)') &
!          vUp, Vlow, Jlow, branch, gf(i), wavenumber(i),excitation(i)

     read(20,'(4x,i3,i3,i3,2x,a4,4x,e9.4,1x,f10.3,f6.3)') &
          vUp, Vlow, Jlow(i),branch,gf(i),wavenumber(i),excitation(i)

!     write(*,*) vUp, Vlow, Jlow, branch, gf(i), wavenumber(i),excitation(i)

  enddo

  close(20)

  lambda(1:nLines) = (1./wavenumber(1:nLines)) * 1.e8 ! wavenumbers to angs

  g(1:nLines) = (2.*real(Jlow(1:nLines))+1.)

  alpha(1:nLines) = ((pi*eCharge**2)/(mElectron*cSpeed)) * gf(1:nLines) / g(1:nLines)

  kappa(1:nLines) = alpha(1:nLines) / (16.*mHydrogen + 48.*mHydrogen) ! divide by mass of Tio

  kappa(1:nLines) = kappa(1:nLines) / (4.*pi) ! ?!

  kappa(1:nLines) = kappa(1:nLines) * 1.e10 ! to code units

  write(*,*) "Done."

end subroutine readTio

subroutine createSample(temperature, lamArray, kapArray, nLam, nLines, lambda, &
     kappa, g, excitation, mu)
  real :: temperature
  real :: vTherm
  real(double):: partFunc
  integer :: i, j
  integer :: nLam
  real :: lam1, lam2
  integer :: i1, i2
  integer :: nLines
  real :: lamArray(:)
  real :: kapArray(:)
  real :: lambda(:)
  real :: kappa(:)
  real :: g(:)
  real :: excitation(:)
  real :: dv
  real :: gaussFac
  real :: fac
  real :: mu
  


  vTherm = sqrt(3.*kerg*temperature/(mu*mHydrogen))
  vTherm = sqrt(vTherm**2 + (2.e5)**2)

  partFunc = 0.
  do i = 1, nLines
     partFunc = partFunc + dble(g(i))*exp(-dble(excitation(i)/(kEv * temperature)))
 enddo

  kapArray(1:nLam)  = 0.

  do j = 1, nLam
     lam1 = (1.-10.*vTherm/cSpeed)*lamArray(j)
     lam2 = (1.+10.*vTherm/cSpeed)*lamArray(j)
     call locate(lambda, nLines, lam1, i1)
     call locate(lambda, nLines, lam2, i2)
     do i = i1,i2
        dv = cSpeed * (lambda(i)-lamArray(j))/lamArray(j)
        gaussFac = 1./(vTherm*sqrt(2.*pi))*exp(-(dv**2)/(2.*vTherm**2))
        if (partfunc /= 0.) then
           fac =  exp(-dble(excitation(i)/(kEv * temperature))) / partFunc * dble(g(i))
        else
           fac = 0.
        endif
        kapArray(j) = kapArray(j) + kappa(i)*fac*gaussFac
     enddo
  enddo
end subroutine createSample

subroutine createKappaGrid(lookuptable, nTemps, t1, t2, nLam, lamArray)
  implicit none
  real :: lamArray(:)
  type(molecularKappaGrid) :: lookuptable 
  integer :: nTemps, nLam
  real :: t1, t2
  integer :: i
  integer :: nLines
  real, allocatable :: excitation(:)
  real, allocatable  :: g(:)
  real, allocatable :: kappa(:)
  real, allocatable :: lambda(:)
  logical :: gridReadFromdisc
  gridReadFromDisc = .false.

  nLines = 0

  write(*,*) "Creating opacity look-up table..."

  write(*,*) "Attempting to read a previous look-up table: ",trim(lookuptable%filename)

  call readKappaGrid(lookuptable,nTemps, nLam, gridReadFromDisc)

  if (.not.gridReadFromDisc) then
     write(*,*) "Correct lookup table not available, creating..."
     LookupTable%nTemps = nTemps
     LookupTable%nLam = nLam
     
     allocate(LookupTable%tempArray(1:nTemps))
     allocate(LookupTable%lamArray(1:nLam))
     allocate(LookupTable%kapArray(1:nTemps, 1:nLam))
     
     LookupTable%lamArray(1:nLam) = lamArray(1:nLam)
     

     select case(Lookuptable%molecule)

        case("TiO")

           allocate(lambda(1:12000000))
           allocate(g(1:12000000))
           allocate(kappa(1:12000000))
           allocate(excitation(1:12000000))
           call readTio(nLines, lambda, kappa, excitation, g)

        case("H2O")
           allocate(lambda(1:65912356))
           allocate(g(1:65912356))
           allocate(kappa(1:65912356))
           allocate(excitation(1:65912356))
           call readh2o(nLines, lambda, kappa, excitation, g)

        case DEFAULT
           allocate(lambda(1:12000000))
           allocate(g(1:12000000))
           allocate(kappa(1:12000000))
           allocate(excitation(1:12000000))
           call readKurucz(lookuptable%source, nLines, lambda, kappa, excitation, g, lookuptable%mu)

     end select
     
     do i = 1, nTemps
        LookupTable%tempArray(i) = t1 + (t2-t1)*real(i-1)/real(nTemps-1)
     enddo
     do i = 1, nTemps
        call createSample(LookupTable%tempArray(i), lamArray, LookupTable%kapArray(i,1:nLam), nLam, &
             nLines, lambda, kappa, g, excitation, lookuptable%mu)
     enddo
     call writeKappaGrid(lookuptable)
     deallocate(lambda, g, kappa, excitation)
  endif


  write(*,*) lookuptable%molecule," opacity lookup table complete."
end subroutine createKappaGrid

subroutine readKappaGrid(lookuptable,  wantedntemps, wantednlam, gridreadfromdisc)
  type(molecularKappaGrid) :: lookuptable 
  integer :: wantednTemps, wantednlam
  integer :: nTemps, nLam
  logical :: gridreadFromDisc

  gridReadFromDisc = .false.

  open(21, file=lookuptable%filename, status="old", form="unformatted",err=666)
  read(21) lookuptable%molecule
  read(21) lookuptable%mu
  read(21) nTemps, nLam

  if ((nTemps /= wantedNtemps).or.(nLam /= wantednLam)) then
     write(*,*) "TiO lookup-table on disc is not of correct size"
     goto 666
  endif

  TioLookupTable%nTemps = nTemps
  TioLookupTable%nLam = nLam
     
  allocate(LookupTable%tempArray(1:nTemps))
  allocate(LookupTable%lamArray(1:nLam))
  allocate(LookupTable%kapArray(1:nTemps, 1:nLam))   
  read(21) LookupTable%tempArray(1:nTemps)
  read(21) LookupTable%lamArray(1:nLam)
  read(21) LookupTable%kapArray(1:nTemps, 1:nLam)
  close(21)
  gridReadFromDisc = .true.
666 continue
end subroutine readKappaGrid

subroutine writeKappaGrid(lookuptable)
  type(molecularKappaGrid) :: lookuptable 

  open(21, file=lookuptable%filename, status="unknown", form="unformatted")
  write(21) LookupTable%molecule
  write(21) lookupTable%mu
  write(21) LookupTable%nTemps, LookupTable%nLam
  write(21) LookupTable%tempArray(1:LookupTable%nTemps)
  write(21) LookupTable%lamArray(1:LookupTable%nLam)
  write(21) LookupTable%kapArray(1:LookupTable%nTemps, 1:LookupTable%nLam)
  close(21)

end subroutine writeKappaGrid


subroutine returnKappaArray(temperature, LookupTable, kappaAbs, KappaSca)
  type(molecularKappaGrid) :: lookuptable 
  real :: temperature
  real,optional :: kappaAbs(:),kappaSca(:)
  real :: fac
  integer :: i

  if ((temperature < LookupTable%tempArray(1)).or.(temperature > LookupTable%tempArray(LookupTable%nTemps))) then
     write(*,*) "! temperature is outside opacity array bounds"
     write(*,*) "temperature",temperature,LookupTable%temparray(1), LookupTable%temparray(LookupTable%nTemps)
     stop
  endif
  
  call locate(LookupTable%tempArray, LookupTable%nTemps, temperature, i)

  fac = (temperature-LookupTable%tempArray(i))/(LookupTable%tempArray(i+1)-LookupTable%tempArray(i))
  if (PRESENT(kappaAbs)) then 
     kappaAbs(1:LookupTable%nLam) = LookupTable%kapArray(i,1:LookupTable%nLam) + &
          fac*(LookupTable%kapArray(i+1,1:LookupTable%nLam)-LookupTable%kapArray(i,1:LookupTable%nLam))
  endif
  if (PRESENT(kappaSca)) then
     kappaSca(1:LookupTable%nLam) = 1.e-30
  endif
end subroutine returnKappaArray

subroutine returnGasKappaValue(temperature, rho, lambda, kappaAbs, kappaSca, kappaAbsArray, kappaScaArray)
  type(molecularKappaGrid) :: LookupTable
  real :: temperature
  real, optional :: lambda
  real(double) :: rho
  real, optional :: kappaAbs, kappaSca, kappaAbsArray(:), kappaScaArray(:)
  real :: t1, t2
  real :: pressure, mu
  integer :: i, j

  lookuptable = tioLookuptable

  mu = 2.46

  pressure = rho * kErg * temperature / (mu * mHydrogen)

  if ((temperature < LookupTable%tempArray(1)).or.(temperature > LookupTable%tempArray(LookupTable%nTemps))) then
     write(*,*) "! temperature is outside opacity array bounds"
     write(*,*) "temperature",temperature,LookupTable%temparray(1),LookupTable%temparray(LookupTable%nTemps)
     stop
  endif
  
  call locate(LookupTable%tempArray, LookupTable%nTemps, temperature, i)
  t1 = (temperature-LookupTable%tempArray(i))/(LookupTable%tempArray(i+1)-LookupTable%tempArray(i))

  if (present(lambda)) then
     call locate(LookupTable%lamArray, LookupTable%nLam, lambda, j)
     t2 = (lambda-LookupTable%lamArray(j))/(LookupTable%lamArray(j+1)-LookupTable%lamArray(j))
  endif

  if (PRESENT(kappaAbsArray)) then
     call returnKappaArray(temperature, LookupTable, kappaAbs=kappaAbsArray)
  endif
  if (PRESENT(kappaScaArray)) then
     call returnKappaArray(temperature, LookupTable, kappaSca=kappaScaArray)
  endif

  if (PRESENT(kappaAbs)) then
     kappaAbs = (1.-t1) * (1.-t2) * lookupTable%kapArray(i  , j  ) + &
                (   t1) * (1.-t2) * lookupTable%kapArray(i+1, j  ) + &
                (1.-t1) * (   t2) * lookupTable%kapArray(i  , j+1) + &
                (   t1) * (   t2) * lookupTable%kapArray(i+1, j+1) 
     kappaAbs = kappaAbs * fracMolecule(temperature, pressure, 12)
!     write(*,*) fracMolecule(temperature, pressure, 12)
  endif

  if (PRESENT(kappaSca)) then
     kappaSca = 0.
!     kappaSca = kappaSca + (hydrogenRayXsection(lambda)/mHydrogen)*1.e10
  endif

end subroutine returnGasKappaValue
  
real function hydrogenRayXsection(lambda)

  implicit none
  real :: lambda ! in angstrom
  real(double) :: freq, fac, freq0, gamma0
  real :: einsteinA(1:6), lambdaTrans(1:6)
  integer :: i

  data einsteinA  / 0.000D0 , 4.699D8 , 5.575D7 , 1.278D7 , 4.125D6 , 1.644D6 /
  data lambdaTrans / 0.00000D-8 , 1215.67D-8 , 1025.72D-8 , 992.537D-8 , 949.743D-8 , 937.803D-8 /

  hydrogenRayXsection = 0.
  do i = 2, 6

     freq = dble(cSpeed) / dble(lambda*angstromtoCm)
     freq0 = dble(cSpeed) / dble(lambdaTrans(i))
     gamma0 =  dble(einsteinA(i)/fourPi)

     fac = freq**4 / ((freq0**2 - freq**2)**2 + (gamma0*freq)**2)

     hydrogenRayXsection = hydrogenRayXsection + sigmaE * fac

  enddo

end function hydrogenRayXsection

subroutine readTsujiPPTable()
  implicit none
!  type(tsujiPPtable) :: tsujiPPlookuptable
  integer :: i, j
  character(len=80) :: filename, dataDirectory

  dataDirectory = " "
  tsujiPPlookuptable%nTemps = 16
  tsujiPPlookuptable%nAtoms = 5
  tsujiPPlookuptable%nPressures = 9
  allocate(tsujiPPlookuptable%temperature(1:16))
  allocate(tsujiPPlookuptable%logPressure(1:9))
  allocate(tsujiPPlookuptable%pp(1:9,1:5,1:16))

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"tsuji_pp.dat"

  open(20,file=filename,form="formatted",status="old")
  read(20,*) tsujiPPlookuptable%temperature(1:16)
  tsujiPPlookuptable%temperature(1:16) = 5040.0 / tsujiPPlookuptable%temperature(1:16)
  do i = 1, 9
     tsujiPPlookuptable%logPressure(i) = real(i-1)
     do j = 1, 5
        read(20,'(15x,16f7.2)') tsujiPPlookuptable%pp(i,j,1:16)
     enddo
  enddo
  close(20)
end subroutine readTsujiPPTable

subroutine readTsujiKPtable()
  implicit none
  character(len=80) :: junk, dataDirectory, filename
  integer :: i

  datadirectory = " "
  tsujikplookuptable%nMolecules = 12
  allocate(tsujikplookuptable%molecule(1:12))
  allocate(tsujikplookuptable%a(1:12,1:5))

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"tsuji_kp.dat"

  open(20, file=filename, form = "formatted", status="old")
  read(20,'(a)') junk
  do i = 1, 12
     read(20,'(a3, 1x, 5e12.4)') tsujikplookuptable%molecule(i), tsujikplookuptable%a(i,1:5)
  enddo
  close(20)
end subroutine readTsujiKPtable


real function logKp(temperature, nMol)
  implicit none
  real :: temperature
  integer :: nMol
  real :: theta, thisTemp

  logKp = 0.

  thisTemp  = min(max(temperature, tsujiPPlookuptable%temperature(1)), &
       tsujiPPlookuptable%temperature(tsujiPPlookuptable%nTemps))

  theta = 5040./thisTemp

  logKp = tsujikplookuptable%a(nMol,1) + &
       tsujikplookuptable%a(nMol,2) * theta + &
       tsujikplookuptable%a(nMol,3) * theta**2 + &
       tsujikplookuptable%a(nMol,4) * theta**3 + &
       tsujikplookuptable%a(nMol,5) * theta**4

end function logKp

real function PPatom(temperature, pressure, nAtom)
  implicit none
  real :: temperature, pressure
  integer :: nAtom
  real :: thisTemp
  real :: logP
  integer :: ip, it
  real :: t1, t2


  logP = log10(pressure)
  logP = min(max(logP, tsujiPPlookuptable%logPressure(1)), &
       tsujiPPlookuptable%logPressure(tsujiPPlookuptable%nPressures))
  thisTemp  = min(max(temperature, tsujiPPlookuptable%temperature(1)), &
       tsujiPPlookuptable%temperature(tsujiPPlookuptable%nTemps))

  call locate(tsujiPPlookuptable%logPressure, tsujiPPlookuptable%nPressures, &
       logP, ip)
  t1 = (logP - tsujiPPlookuptable%logPressure(ip)) / &
       (tsujiPPlookuptable%logPressure(ip+1) - tsujiPPlookuptable%logPressure(ip))
  
  call locate(tsujiPPlookuptable%temperature, tsujiPPlookuptable%nTemps, &
       thistemp, it)
  t2 = (thistemp - tsujiPPlookuptable%temperature(it)) / &
       (tsujiPPlookuptable%temperature(it+1) - tsujiPPlookuptable%temperature(it))
  PPatom = (1.-t1) * (1.-t2) * tsujipplookupTable%pp(ip  , nAtom, it  ) + &
             (   t1) * (1.-t2) * tsujipplookupTable%pp(ip+1, nAtom, it  ) + &
             (1.-t1) * (   t2) * tsujipplookupTable%pp(ip  , nAtom, it+1) + &
             (   t1) * (   t2) * tsujipplookupTable%pp(ip+1, nAtom, it+1) 


  PPatom = 10.**PPatom
end function PPatom
 
real function fracMolecule(temperature, pressure, nMol)
  implicit none
  real :: temperature, pressure
  integer :: nMol
  real :: kp
  real :: p1, p2
  real :: logP, ppMolecule

  logP = log10(pressure)
  logP = min(max(logP, tsujiPPlookuptable%logPressure(1)), &
       tsujiPPlookuptable%logPressure(tsujiPPlookuptable%nPressures))


  kp = 10.**(logKp(temperature, nMol))

  select case(nMol)

     case(1) ! H2
        p1 = PPatom(temperature, pressure, 1)
        ppMolecule = p1 * p1 / kp

     case(2) ! C2
        p1 = PPatom(temperature, pressure, 2)
        ppMolecule = p1 * p1 / kp

     case(3) ! CH
        p1 = PPatom(temperature, pressure, 1)
        p2 = PPatom(temperature, pressure, 2)
        ppMolecule = p1 * p2 / kp

     case(4) ! N2
        p1 = PPatom(temperature, pressure, 3)
        ppMolecule = p1 * p1 / kp

     case(5) ! NH
        p1 = PPatom(temperature, pressure, 1)
        p2 = PPatom(temperature, pressure, 3)
        ppMolecule = p1 * p2 / kp

     case(6) ! NO 
        p1 = PPatom(temperature, pressure, 3)
        p2 = PPatom(temperature, pressure, 4)
        ppMolecule = p1 * p2 / kp

     case(7) ! O2
        p1 = PPatom(temperature, pressure, 4)
        ppMolecule = p1 * p1 / kp
        
     case(8) ! OH
        p1 = PPatom(temperature, pressure, 1)
        p2 = PPatom(temperature, pressure, 4)
        ppMolecule = p1 * p2 / kp
        
     case(9) ! H2O
        p1 = PPatom(temperature, pressure, 1)
        p2 = PPatom(temperature, pressure, 4)
        ppMolecule = p1 * p1 * p2 / kp

     case(10) ! CO
        p1 = PPatom(temperature, pressure, 2)
        p2 = PPatom(temperature, pressure, 4)
        ppMolecule = p1 * p2 / kp

     case(11) ! CO2
        p1 = PPatom(temperature, pressure, 2)
        p2 = PPatom(temperature, pressure, 4)
        ppMolecule = p1 * p2 * p2 / kp

     case(12) ! TiO
        p1 = PPatom(temperature, pressure, 4)
        p2 = PPatom(temperature, pressure, 5)
        ppMolecule = p1*p2 / kp

     case DEFAULT
       write(*,*) "Molecule number ",nMol, " out of range."
       stop

   end select

   fracMolecule = ppMolecule/(10.**logP)

 end function fracMolecule

subroutine readKurucz(filename,nLines,lambda,kappa,excitation,g,mu)
  use utils_mod, only: wavenumbertoEv

  implicit none
  character(len=*) :: filename
  integer :: nLines
  real :: lambda(:),kappa(:),excitation(:),g(:)
  integer :: i
  real, allocatable :: jLow(:)
  real, allocatable :: gf(:),alpha(:)
  real, parameter :: mHydrogen = 1.6733e-24      ! g
  real, parameter :: mElectron = 9.109565D-28    ! g
  real, parameter :: eCharge = 4.803242384E-10 !

  real, parameter :: cSpeed = 2.99792458e10  ! cm/s
  real, parameter :: pi = 3.141592654
  real :: mu

  open(20,file=filename,status="old",form="formatted")

  read(20,*) nLines
  allocate(gf(1:nLines))
  allocate(alpha(1:nLines))
  allocate(jlow(1:nLines))



  write(*,*) "Reading from ",trim(filename)

  do i = 1, nLines

     if (mod(i,nLines/10) == 0) write(*,*) i 

     read(20,'(f10.4,f7.3,f5.1,f10.3)') &
          lambda(i), gf(i),  jLow(i), excitation(i)

  enddo

  close(20)

  lambda(1:nLines) = lambda(1:nLines) * 10.  ! nm to angstroms

  g(1:nLines) = (2.*real(Jlow(1:nLines))+1.)

  alpha(1:nLines) = ((pi*eCharge**2)/(mElectron*cSpeed)) * gf(1:nLines) / g(1:nLines)

  kappa(1:nLines) = alpha(1:nLines) / (mu * mHydrogen) ! divide by mass of molecule

  kappa(1:nLines) = kappa(1:nLines) * 1.e10 ! to code units
  
  excitation(1:Nlines) = abs(excitation(1:Nlines)) ! get rid of negative flags

  call wavenumberToEv(excitation, nLines)

  write(*,*) "Done."

end subroutine readKurucz

subroutine createAllMolecularTables(nTemps, t1, t2, nLam, lamArray)
  implicit none
  integer :: nTemps, nLam
  real :: t1, t2, lamArray(:)
  character(len=80) :: dataDirectory, filename
  integer :: i
  dataDirectory = " "
  ! h2o

  h2oLookupTable%molecule = "H2O"
  h2oLookupTable%mu = 1.+1.+16.
  h2oLookupTable%filename = "h2olookuptable.dat"
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"h2o.dat"
  h2oLookupTable%source = filename 
  call createKappaGrid(h2olookuptable, nTemps, t1, t2, nLam, lamArray)

  ! tio

  tioLookupTable%molecule = "TiO"
  tioLookupTable%mu = 48.+16.
  tioLookupTable%filename = "tiolookuptable.dat"
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"tio.dat"
  tioLookupTable%source = filename 
  call createKappaGrid(tiolookuptable, nTemps, t1, t2, nLam, lamArray)


  ! ch

  chLookupTable%molecule = "CH"
  chLookupTable%mu = 12.+1.
  chLookupTable%filename = "chlookuptable.dat"
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"ch.asc"
  chLookupTable%source = filename 
  call createKappaGrid(chlookuptable, nTemps, t1, t2, nLam, lamArray)

  ! nh

  nhLookupTable%molecule = "NH"
  nhLookupTable%mu = 14.+1.
  nhLookupTable%filename = "nhlookuptable.dat"
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"nh.asc"
  nhLookupTable%source = filename 
  call createKappaGrid(nhlookuptable, nTemps, t1, t2, nLam, lamArray)

  ! oh

  ohLookupTable%molecule = "OH"
  ohLookupTable%mu = 16.+1.
  ohLookupTable%filename = "ohlookuptable.dat"
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"oh.asc"
  ohLookupTable%source = filename 
  call createKappaGrid(ohlookuptable, nTemps, t1, t2, nLam, lamArray)

  ! h2

  h2LookupTable%molecule = "H2"
  h2LookupTable%mu = 2.
  h2LookupTable%filename = "h2lookuptable.dat"
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"h2.asc"
  h2LookupTable%source = filename 
  call createKappaGrid(h2lookuptable, nTemps, t1, t2, nLam, lamArray)

  ! co

  coLookupTable%molecule = "CO"
  coLookupTable%mu = 12.+16.
  coLookupTable%filename = "colookuptable.dat"
  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//"co.asc"
  coLookupTable%source = filename 
  call createKappaGrid(colookuptable, nTemps, t1, t2, nLam, lamArray)

end subroutine createAllMolecularTables

subroutine readH2O(nLines,lambda,kappa,excitation,g)
  use utils_mod, only: wavenumbertoEv, locate, convertbyte2, convertbyte4


  implicit none
  integer :: nLines
  integer :: i, j, k
  real :: lambda(:), kappa(:), excitation(:), g(:)
  real, allocatable :: wavenumber(:),gf(:)
  real, parameter :: mHydrogen = 1.6733e-24      ! g
  real, parameter :: mElectron = 9.109565D-28    ! g
  real, parameter :: eCharge = 4.803242384E-10 !

  real, parameter :: cSpeed = 2.99792458e10  ! cm/s
  real, parameter :: pi = 3.141592654
  integer(kind=single) :: ibyte(4), ibyte2(2), ibyte3(2)
  integer :: iwl
  integer(kind=double) :: igflog, ielo
  integer :: iso
  real(kind=double) :: ratiolog, wlvac
  real(kind=double) :: freq, congf, elo
  real :: TABLOG(32768)
  real :: xiso(4),  e(273201)
  INTEGER :: QUANTUMS(8,273201)
  REAL :: WT(273201)

  DATA XISO/  .9976,  .0004,  .0020, .00001/

  ielo = 0; iwl = 0
  ratiolog = LOG(1.D0+1.D0/2000000.D0)


  nLines = 65912356  / 10

  DO  I=1,32768
     TABLOG(I)=10.**((I-16384)*.001)
  enddo

  allocate(gf(1:nLines))
  allocate(wavenumber(1:nLines))


  open(20,file="/home/th/molecules/h2ofast.bin",status="old",form="unformatted", &
       recl=8, access="direct")




  do i = 1, nLines

     if (mod(i,nLines/10) == 0) write(*,*) i 

     READ(20,REC=i) ibyte(1:4) , ibyte2(1:2), ibyte3(1:2)

     call convertByte4(ibyte, iwl)
     WLVAC=EXP(IWL*RATIOLOG)
     FREQ=2.99792458E17/WLVAC
     call convertByte2(ibyte2, ielo)
     call convertByte2(ibyte3, igflog)
     ISO=1
     IF(IELO.GT.0.AND.IGFLOG.GT.0)GO TO 19
     ISO=2
     IF(IELO.GT.0)GO TO 19
     ISO=3
     IF(IGFLOG.GT.0)GO TO 19
     ISO=4

19 continue
      elo=  abs(ielo)
      igflog = abs(igflog)
     CONGF=.01502*TABLOG(IGFLOG) * XISO(ISO)
!     write(*,*) wlvac,congf,elo

     lambda(i) = wlvac*10. ! nm to angstroms
     gf(i) = congf
     wavenumber(i) = elo
     excitation(i) = elo
  enddo
  close(20)
  open(20,file="/h/th/molecules/eh2opartridge.asc",status="old", form="formatted")

  READ(20,1)
  READ(20,1)
  READ(20,1)
  READ(20,1)
  DO  I=1,273201
     READ(20,1) E(I),(QUANTUMS(K,I),K=1,8)
1    FORMAT(10X,F12.5,2I2,6I3)
     WT(I)=2.*QUANTUMS(3,I)+1.
     IF(QUANTUMS(2,I).EQ.1)WT(I)=WT(I)*3
  enddo
  close(20)

  write(*,*) "Associating statistical weights..."
  do i = 1, nLines
     if (mod(i,nLines/10) == 0) write(*,*) i 
     call locate(e, 273201, wavenumber(i),j)
     g(i) = wt(j)
  enddo

  kappa(1:nLines) = ((pi*eCharge**2)/(mElectron*cSpeed)) * gf(1:nLines) / g(1:nLines) &
       / (16.*mHydrogen + 2.*mHydrogen) ! divide by mass of h2o

  kappa(1:nLines) = kappa(1:nLines) * 1.e10 ! to code units

  call wavenumberToEv(excitation, nLines)

  write(*,*) "Done."

end subroutine readH2O

end module gas_opacity_mod
