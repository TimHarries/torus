module gas_opacity_mod

use constants_mod
use utils_mod

implicit none

type TioKappaGrid
   integer :: nLam
   integer :: nTemps
   real, pointer :: lamArray(:) => null()
   real, pointer :: tempArray(:) => null()
   real, pointer :: kapArray(:,:) => null()
end type TioKappaGrid

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


type(tioKappaGrid), save :: TioLookuptable
type(tsujiPPtable), save :: tsujipplookuptable
type(tsujiKPtable), save :: tsujikplookuptable

contains

subroutine readTio(nLines,lambda,kappa,excitation,g)

  implicit none
  integer :: nLines
  real :: lambda(*),kappa(*),excitation(*),g(*)
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

subroutine createTiOSample(temperature, lamArray, kapArray, nLam, nLines, lambda, &
     kappa, g, excitation)
  real :: temperature
  real :: vTherm
  real :: partFunc
  integer :: i, j
  integer :: nLam
  real :: lam1, lam2
  integer :: i1, i2
  integer :: nLines
  real :: lamArray(*)
  real :: kapArray(*)
  real :: lambda(*)
  real :: kappa(*)
  real :: g(*)
  real :: excitation(*)
  real :: dv
  real :: gaussFac
  real :: fac
  


  vTherm = sqrt(3.*kerg*temperature/(16.*mHydrogen + 48.*mHydrogen))
  vTherm = sqrt(vTherm**2 + (2.e5)**2)

  partFunc = 0.
  do i = 1, nLines
     partFunc = partFunc + g(i)*exp(-excitation(i)/(kEv * temperature))
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
        fac =  exp(-excitation(i)/(kEv * temperature)) / partFunc * g(i)
        kapArray(j) = kapArray(j) + kappa(i)*fac*gaussFac
     enddo
  enddo
end subroutine createTiOSample

subroutine createTioGrid(nTemps, t1, t2, nLam, lamArray)
  implicit none
  real :: lamArray(*)
  integer :: nTemps, nLam
  real :: t1, t2
  integer :: i
  integer :: nLines
  real :: excitation(12000000)
  real :: g(12000000)
  real :: kappa(12000000)
  real :: lambda(12000000)
  logical :: gridReadFromdisc

  write(*,*) "Creating TiO opacity look-up table..."

  write(*,*) "Attempting to read a previous look-up table."

  call readTioGrid(nTemps, nLam, gridReadFromDisc)

  if (.not.gridReadFromDisc) then
     write(*,*) "Correct lookup table not available, creating..."
     TioLookupTable%nTemps = nTemps
     TioLookupTable%nLam = nLam
     
     allocate(TioLookupTable%tempArray(1:nTemps))
     allocate(TioLookupTable%lamArray(1:nLam))
     allocate(TioLookupTable%kapArray(1:nTemps, 1:nLam))
     
     TioLookupTable%lamArray(1:nLam) = lamArray(1:nLam)
     
     call readTio(nLines, lambda, kappa, excitation, g)
     
     do i = 1, nTemps
        TioLookupTable%tempArray(i) = t1 + (t2-t1)*real(i-1)/real(nTemps-1)
     enddo
     do i = 1, nTemps
        call createTioSample(TioLookupTable%tempArray(i), lamArray, TioLookupTable%kapArray(i,1:nLam), nLam, &
             nLines, lambda, kappa, g, excitation)
     enddo
     call writeTioGrid()
  endif
  write(*,*) "TiO opacity lookup table complete."
end subroutine createTioGrid

subroutine readTioGrid(wantedntemps, wantednlam, gridreadfromdisc)
  integer :: wantednTemps, wantednlam
  integer :: nTemps, nLam
  logical :: gridreadFromDisc

  gridReadFromDisc = .false.

  open(21, file="tiolookuptable.dat", status="old", form="unformatted",err=666)
  read(21) nTemps, nLam

  if ((nTemps /= wantedNtemps).or.(nLam /= wantednLam)) then
     write(*,*) "TiO lookup-table on disc is not of correct size"
     goto 666
  endif

  TioLookupTable%nTemps = nTemps
  TioLookupTable%nLam = nLam
     
  allocate(TioLookupTable%tempArray(1:nTemps))
  allocate(TioLookupTable%lamArray(1:nLam))
  allocate(TioLookupTable%kapArray(1:nTemps, 1:nLam))   
  read(21) TioLookupTable%tempArray(1:nTemps)
  read(21) TioLookupTable%lamArray(1:nLam)
  read(21) TioLookupTable%kapArray(1:nTemps, 1:nLam)
  close(21)
  gridReadFromDisc = .true.
666 continue
end subroutine readTioGrid

subroutine writeTioGrid()

  open(21, file="tiolookuptable.dat", status="unknown", form="unformatted")
  write(21) tioLookupTable%nTemps, tioLookupTable%nLam
  write(21) TioLookupTable%tempArray(1:TioLookupTable%nTemps)
  write(21) TioLookupTable%lamArray(1:TioLookupTable%nLam)
  write(21) TioLookupTable%kapArray(1:TioLookupTable%nTemps, 1:TioLookupTable%nLam)
  close(21)

end subroutine writeTioGrid


subroutine returnKappaArray(temperature, TioLookupTable, kappaAbs, KappaSca)
  real :: temperature
  real :: kappaAbs(*),kappaSca(*)
  type(TioKappaGrid) :: TioLookupTable
  real :: fac
  integer :: i

  if ((temperature < TioLookupTable%tempArray(1)).or.(temperature > TioLookupTable%tempArray(TioLookupTable%nTemps))) then
     write(*,*) "! TiO temperature is outside opacity array bounds"
     write(*,*) "temperature",temperature,TioLookupTable%temparray(1),TioLookupTable%temparray(TioLookupTable%nTemps)
     stop
  endif
  
  call locate(TioLookupTable%tempArray, TioLookupTable%nTemps, temperature, i)

  fac = (temperature-TioLookupTable%tempArray(i))/(TioLookupTable%tempArray(i+1)-TioLookupTable%tempArray(i))
  kappaAbs(1:TioLookupTable%nLam) = TioLookupTable%kapArray(i,1:TioLookupTable%nLam) + &
       fac*(TioLookupTable%kapArray(i+1,1:TioLookupTable%nLam)-TioLookupTable%kapArray(i,1:TioLookupTable%nLam))
  kappaSca(1:TioLookupTable%nLam) = 1.e-30
end subroutine returnKappaArray

subroutine returnGasKappaValue(temperature, rho, lambda, kappaAbs, kappaSca)
  real :: temperature
  real :: lambda, rho
  real, optional :: kappaAbs, kappaSca
  real :: t1, t2
  real :: pressure, mu
  integer :: i, j

  mu = 2.46

  pressure = rho * kErg * temperature / (mu * mHydrogen)

  if ((temperature < TioLookupTable%tempArray(1)).or.(temperature > TioLookupTable%tempArray(TioLookupTable%nTemps))) then
     write(*,*) "! TiO temperature is outside opacity array bounds"
     write(*,*) "temperature",temperature,TioLookupTable%temparray(1),TioLookupTable%temparray(TioLookupTable%nTemps)
     stop
  endif
  
  call locate(TioLookupTable%tempArray, TioLookupTable%nTemps, temperature, i)
  call locate(TioLookupTable%lamArray, TioLookupTable%nLam, lambda, j)

  t1 = (temperature-TioLookupTable%tempArray(i))/(TioLookupTable%tempArray(i+1)-TioLookupTable%tempArray(i))
  t2 = (lambda-TioLookupTable%lamArray(j))/(TioLookupTable%lamArray(j+1)-TioLookupTable%lamArray(j))

  if (PRESENT(kappaAbs)) then
     kappaAbs = (1.-t1) * (1.-t2) * tiolookupTable%kapArray(i  , j  ) + &
                (   t1) * (1.-t2) * tiolookupTable%kapArray(i+1, j  ) + &
                (1.-t1) * (   t2) * tiolookupTable%kapArray(i  , j+1) + &
                (   t1) * (   t2) * tiolookupTable%kapArray(i+1, j+1) 
     kappaAbs = kappaAbs * ppMolecule(temperature, pressure, 12)/pressure
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
  character(Len=80) :: junk
  real :: pp(9,5,16)
  integer :: i, j
  character(len=80) :: filename, dataDirectory

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
  real :: logP, theta
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
 
real function PPmolecule(temperature, pressure, nMol)
  implicit none
  real :: temperature, pressure
  integer :: nMol
  real :: kp
  real :: p1, p2

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

 end function PPmolecule

end module gas_opacity_mod
