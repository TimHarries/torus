module gas_opacity_mod

use constants_mod
use utils_mod, only: locate
use unix_mod, only: unixGetenv
use gridtype_mod, only : GRIDTYPE
use messages_mod
use mpi_global_mod
use parallel_mod

#ifdef NETCDF
use netcdf
#endif

implicit none

type lineListType
   integer :: nLines
   real(double), pointer :: wavelength(:) => null()
   real(double), pointer :: logGF(:) => null()
   real(double), pointer :: Jlower(:) => null()
   real(double), pointer :: eLower(:) => null()
   real(double), pointer :: JUpper(:) => null()
   real(double), pointer :: eUpper(:) => null()
end type lineListType

type molecularKappaGrid
   character(len=80) :: filename
   character(len=80) :: source
   character(len=5) :: molecule
   real :: abundance ! relative to N(H)
   real :: mu
   integer :: nLam
   integer :: nTemps
   integer :: nPressures
   real(double), pointer :: lamArray(:) => null()
   real(double), pointer :: tempArray(:) => null()
   real(double), pointer :: pressureArray(:) => null()
   real(double), pointer :: kapArray(:,:) => null()
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


type(molecularKappaGrid), save :: Lookuptable(10)
integer :: nLookUpTables


type(tsujiPPtable), save :: tsujipplookuptable
type(tsujiKPtable), save :: tsujikplookuptable

contains

subroutine freeTable(lookuptable)
  type(molecularKappaGrid) :: lookupTable
  if (associated(lookuptable%lamArray)) then
     deallocate(lookuptable%lamArray)
  endif

  if (associated(lookuptable%tempArray)) then
     deallocate(lookuptable%tempArray)
  endif

  if (associated(lookuptable%pressureArray)) then
     deallocate(lookuptable%pressureArray)
  endif
  if (associated(lookuptable%kapArray)) then
     deallocate(lookuptable%kapArray)
  endif

end subroutine freeTable

subroutine readTio(nLines,lambda,kappa,excitation,g)

  implicit none
  integer :: nLines
  real :: lambda(:),kappa(:),excitation(:),g(:)
  integer :: vUp, vLow, i
  integer, allocatable :: jLow(:)
  character(len=4) :: branch
  real, allocatable :: wavenumber(:),gf(:),alpha(:)
  real, parameter :: mHydrogen = 1.6733e-24      ! g
  real, parameter :: mElectron = 9.109565e-28    ! g
  real, parameter :: eCharge = 4.803242384E-10 !

  real, parameter :: cSpeed = 2.99792458e10  ! cm/s
  real, parameter :: pi = 3.141592654

  nLines = 11369860 


  allocate(gf(1:nLines))
  allocate(wavenumber(1:nLines))
  allocate(alpha(1:nLines))
  allocate(jlow(1:nLines))


  open(20,file="/home/th/plez/TiO.dat",status="old",form="formatted")
  
  if (writeoutput) write(*,*) "Reading all ",nLines," TiO transitions..."

  do i = 1, nLines

     if (mod(i,nLines/10) == 0) write(*,*) i 

!     read(20,'(4x,i2,i2,i3,1x,a4,e9.2,f10.3,f6.3)') &
!          vUp, Vlow, Jlow, branch, gf(i), wavenumber(i),excitation(i)

     read(20,'(4x,i3,i3,i3,2x,a4,4x,e11.4,1x,f10.3,f6.3)') &
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

  if (writeoutput) write(*,*) "Done."

end subroutine readTio

subroutine createSample(temperature, lamArray, kapArray, nLam, nLines, lambda, &
     kappa, g, excitationLower, excitationUpper, mu, abundance)
  real :: temperature
  real :: vTherm
  real(double):: partFunc
  integer :: i, j
  integer :: nLam
  real :: lam1, lam2
  integer :: i1, i2
  integer :: nLines
  real :: lamArray(:)
  real(double) :: kapArray(:)
  real :: lambda(:)
  real :: kappa(:)
  real :: g(:)
  real :: excitationLower(:),excitationUpper(:)
  real :: dv
  real :: gaussFac
  real :: fac
  real :: mu
  real :: abundance
  


  vTherm = real(sqrt(3.*kerg*temperature/(mu*mHydrogen)))
  vTherm = sqrt(vTherm**2 + (2.e5)**2)

  partFunc = 0.
  do i = 1, nLines
     partFunc = partFunc + dble(g(i))*exp(-dble(excitationLower(i)/(kEv * temperature)))
 enddo

  kapArray(1:nLam)  = 0.

  do j = 1, nLam
     lam1 = real((1.-10.*vTherm/cSpeed)*lamArray(j))
     lam2 = real((1.+10.*vTherm/cSpeed)*lamArray(j))
     call locate(lambda, nLines, lam1, i1)
     call locate(lambda, nLines, lam2, i2)
     do i = i1,i2
        dv = real(cSpeed * (lambda(i)-lamArray(j))/lamArray(j))
        gaussFac = real(1./(vTherm*sqrt(2.*pi))*exp(-(dv**2)/(2.*vTherm**2)))
        if (partfunc /= 0.) then
           fac =  real(exp(-dble((excitationLower(i)/(kEv * temperature))) / partFunc))
           fac = fac * (1.e0 - exp(-(excitationUpper(i) - excitationLower(i))/(real(kEv) * temperature)))
        else
           fac = 0.
        endif
        kapArray(j) = kapArray(j) + kappa(i)*fac*gaussFac
     enddo
  enddo
  kapArray(1:nLam) = kapArray(1:nLam) * abundance 
end subroutine createSample

subroutine createKappaGrid(lookuptable, nLam, lamArray, reset)
  use timing, only : mySleep
  implicit none
  real(double) :: lamArray(:)
  type(molecularKappaGrid) :: lookuptable(:)
  integer :: nLam, i
  logical :: gridReadFromdisc, reset
  gridReadFromDisc = .false.


  do i = 1, nLookupTables
     if (writeoutput) write(*,*) "Creating opacity look-up table..."
     if (.not.reset) then
        if (Writeoutput) write(*,*) "Attempting to read a previous look-up table: ",trim(lookuptable(i)%filename)
        call readKappaGrid(lookuptable(i), nLam, lamArray, gridReadFromDisc)
        call torus_mpi_barrier()
     endif
     if ((.not.gridReadFromDisc).or.reset) then
        if (Writeoutput) write(*,*) "Correct lookup table not available, creating..."
        
        if (myrankGlobal == 0) then
#ifdef NETCDF
           call readAbsCoeffFileAmundsen(lookupTable(i)%source, lookupTable(i), nLam, lamArray)
#else
           write(*,*) "Code not linked with NETCDF. Aborting."
           stop
#endif
           call writeKappaGrid(lookuptable(i))
           call freeTable(lookuptable(i))
           write(*,*) "calling sleep"
           call mySleep(10.d0)
        endif
        call torus_mpi_barrier()
        call waitForFile(lookuptable(i)%filename)
        call readKappaGrid(lookuptable(i), nLam, lamArray, gridReadFromDisc)
        if (.not.gridReadFromDisc) then
           write(*,*) myrankGlobal, " had problem reading: ",trim(lookuptable(i)%filename)
           stop
        endif
     endif
     call torus_mpi_barrier()
     if (writeoutput) write(*,*) lookuptable(i)%molecule," opacity lookup table complete."
  enddo
end subroutine createKappaGrid

subroutine readKappaGrid(lookuptable,  wantednlam, wantedlamArray, gridreadfromdisc)
  type(molecularKappaGrid) :: lookuptable 
  integer :: wantednlam
  integer :: nTemps, nLam
  logical :: gridreadFromDisc
  real(double) :: wantedLamArray(:)
  logical :: resizeLookup
  
  gridReadFromDisc = .false.

  open(21, file=lookuptable%filename, status="old", form="unformatted",err=666)
  read(21) lookuptable%molecule
  read(21) lookuptable%mu
  read(21) lookuptable%abundance
  read(21) nTemps, nLam

  if (nLam == 0) then
     if (writeoutput) write(*,*) "Kappa file has zero wavelength array. Error. "
     goto 666
  endif
  resizeLookup = .false.
  if (nLam /= wantednLam) then
     if (writeoutput) write(*,*) "TiO lookup-table on disc is not of correct size ",nlam,wantedNlam
     if (writeoutput) write(*,*) "Attemping to re-size"
     resizeLookup = .true.
  endif

  LookupTable%nTemps = nTemps
  LookupTable%nLam = nLam
     
  allocate(LookupTable%tempArray(1:nTemps))
  allocate(LookupTable%lamArray(1:nLam))
  allocate(LookupTable%kapArray(1:nTemps, 1:nLam))   
  read(21) LookupTable%tempArray(1:nTemps)
  read(21) LookupTable%lamArray(1:nLam)
  read(21) LookupTable%kapArray(1:nTemps, 1:nLam)
  close(21)
  gridReadFromDisc = .true.

  if (resizeLookup) call resizeLookuptable(lookupTable, wantednLam, wantedLamArray)

666 continue
end subroutine readKappaGrid


subroutine resizeLookuptable(lookuptable, nLambda, lamArray)
  type(molecularKappaGrid) :: lookuptable 
  integer :: nLambda
  real(double) :: lamArray(:), fac
  integer :: i, j, k
  real(double), allocatable :: tempKapArray(:,:)

  allocate(tempKapArray(1:lookuptable%nTemps, 1:nLambda))

  do j = 1, lookupTable%nTemps
     do i = 1, nLambda
        call locate(lookupTable%lamArray, lookuptable%nLam, lamArray(i), k)
        fac = (lamArray(i)-lookupTable%lamArray(k))/(lookupTable%lamArray(k+1)-lookupTable%lamArray(k))
        tempKapArray(j,i) = lookupTable%kapArray(j,k) + fac * (lookupTable%kapArray(j,k+1)-lookupTable%kapArray(j,k))
     enddo
  enddo
  deallocate(lookupTable%lamArray, lookupTable%kapArray)
  allocate(lookupTable%lamArray(1:nLambda))
  lookupTable%nlam = nLambda
  lookupTable%lamArray = lamArray

  allocate(lookupTable%kapArray(1:lookupTable%nTemps,1:nLambda))
  lookupTable%kapArray = tempKapArray
  deallocate(tempKapArray)
end subroutine resizeLookuptable


subroutine writeKappaGrid(lookuptable)
  type(molecularKappaGrid) :: lookuptable 

  open(21, file=lookuptable%filename, status="unknown", form="unformatted")
  write(21) LookupTable%molecule
  write(21) lookupTable%mu
  write(21) lookupTable%abundance
  write(21) LookupTable%nTemps, LookupTable%nLam
  write(21) LookupTable%tempArray(1:LookupTable%nTemps)
  write(21) LookupTable%lamArray(1:LookupTable%nLam)
  write(21) LookupTable%kapArray(1:LookupTable%nTemps, 1:LookupTable%nLam)
  close(21)

end subroutine writeKappaGrid


subroutine returnKappaArray(temperature, LookupTable, kappaAbs, KappaSca)
  type(molecularKappaGrid) :: lookuptable 
  real :: temperature
  real(double) :: thisTemp
  real(double),optional :: kappaAbs(:),kappaSca(:)
  real(double) :: fac
  integer :: i

  thisTemp = max(LookupTable%tempArray(1), min(dble(temperature), LookupTable%tempArray(LookupTable%nTemps)))

!  if ((temperature < LookupTable%tempArray(1)).or.(temperature > LookupTable%tempArray(LookupTable%nTemps))) then
!     write(*,*) "! temperature is outside opacity array bounds"
!     write(*,*) "temperature",temperature,LookupTable%temparray(1), LookupTable%temparray(LookupTable%nTemps)
!     stop
!  endif
  
  call locate(LookupTable%tempArray, LookupTable%nTemps, thisTemp, i)

  fac = (thisTemp-LookupTable%tempArray(i))/(LookupTable%tempArray(i+1)-LookupTable%tempArray(i))
  if (PRESENT(kappaAbs)) then 
     kappaAbs(1:LookupTable%nLam) = LookupTable%kapArray(i,1:LookupTable%nLam) + &
          fac*(LookupTable%kapArray(i+1,1:LookupTable%nLam)-LookupTable%kapArray(i,1:LookupTable%nLam))
  endif
  if (PRESENT(kappaSca)) then
     kappaSca(1:LookupTable%nLam) = 1.e-30
  endif
end subroutine returnKappaArray

subroutine returnGasKappaValue(grid, temperature, rho, lambda, kappaAbs, kappaSca, kappaAbsArray, kappaScaArray, ilambda)
  type(GRIDTYPE) :: grid
  real :: temperature
  real, optional :: lambda
  real(double) :: rho
  logical, save :: firstTime = .true.
  integer,optional :: iLambda
  real, allocatable,save :: rayScatter(:)
  real(double), optional :: kappaAbs, kappaSca, kappaAbsArray(:), kappaScaArray(:)
  real(double),allocatable :: tArray(:)
  real(double) :: nhi !, nh, ne,freq
  integer :: i
  !$OMP THREADPRIVATE (firstTime, rayScatter)

!  Nh = rho/mHydrogen
!  ne = Ne_lte(nh, dble(temperature))
!  nHI = nh - ne
  allocate(tarray(1:LookupTable(1)%nlam))

  nhi = rho/mhydrogen
  if (PRESENT(kappaAbs)) then
     kappaAbs = 0.
     do i = 1, nLookupTables
        call returnKappaArray(temperature, LookupTable(i), kappaAbs=tArray)
        kappaAbs = kappaAbs + tArray(ilambda)
     enddo
  endif
     

  if (PRESENT(kappaAbsArray)) then
     kappaAbsArray = 0.d0
     do i = 1, nLookupTables
        call returnKappaArray(temperature, LookupTable(i), kappaAbs=tarray)        
        if (any(tarray > 1.d20)) then
           write(*,*) i,tarray
        endif
!        write(*,*) "i " ,i 
!        write(*,*) "kabs ",kappaabsArray(1:lookuptable(i)%nlam)
!        write(*,*) "tarr ",tArray(1:lookuptable(i)%nlam)
        kappaAbsArray(1:lookuptable(i)%nlam) = kappaAbsArray(1:lookuptable(i)%nlam) + tarray(1:lookuptable(i)%nlam)
     enddo
  endif

  if (PRESENT(kappaSca)) then
     kappaSca = atomhydrogenRayXsection(dble(lambda))*1.e10 / (2.33d0 * mHydrogen)
!     kappaSca = kappaSca + ne * sigmaE * 1.d10
!     kappaSca = 1.d-30
  endif
  if (PRESENT(kappaScaArray)) then
     if (firstTime) then
        allocate(rayScatter(1:grid%nLambda))
        do i = 1, grid%nLambda
           rayScatter(i) = real(atomhydrogenRayXsection(dble(grid%lamArray(i)))*1.d10)
        enddo
        
        firstTime = .false.
     endif
     kappaScaArray = rayScatter / (2.33d0 * mHydrogen)

  endif

  deallocate(tArray)

end subroutine returnGasKappaValue
  

!function Ne_LTE(nh, T) result (ne)
!  use stateq_mod, only : z_hi
!  real(double) :: nh, T, ne, phiT
!  real(double), parameter  :: CI = 2.07d-16   ! in cgs units

!  phiT = CI*Z_HI(10,T)*(T**(-1.5))*EXP(real(hydE0eV,kind=double)/(kev*T))
!
!  ! Solving for phi(T)*ne^2 + 2ne -nTot =0 and ne+N_H = nTot for ne where
  ! nTot is the number density of particles includeing all species.
  ! ==> phi(T)*ne^2 + ne - N_H =0
  ! Th physical solution  is chosen out of two ...  
  !    Ne = (sqrt(nTot*phiT+1.0_db) -1.0_db)/phiT
!  Ne = (sqrt(4.0_db*NH*phiT+1.0_db) -1.0_db)/(2.0_db*phiT)
!end function Ne_LTE

real(double) function atomhydrogenRayXsection(lambda) result(tot)


  implicit none

! Lee and Kim, 2004, MNRAS, 347, 802
  real(double) :: lambda
  real(double), parameter :: omega_l = twoPi * nuHydrogen
  real(double), parameter :: omega_21 = 0.75d0 * omega_l
  real(double) :: omega,deltaOmega
  integer :: i
  real(double) :: cp(0:9), fac

  data cp / 1.26537d0, 3.73766d0, 8.8127d0, 19.1515d0, 39.919d0, 81.1018d0, 161.896d0, 319.001d0, &
       622.229d0, 1203.82d0 /

  omega = cspeed / (lambda * angstromtocm)
  omega = twoPi * omega

  fac = omega/omega_l

  if (fac < 0.647d0) then

     tot = 0.d0
     do i = 0, 9
        tot = tot + cp(i)*fac**(2*i)
     enddo
     tot = tot * sigmaE * fac**4
  else
     deltaOmega = omega - omega_21
     fac = deltaOmega / omega_21
     tot = (0.0433056/fac**2) * (1.d0 - 1.792*fac - 23.637d0*fac**2 - &
          83.1393*fac**3 - 244.1453d0*fac**4 - 699.473d0*fac**5)
     tot = max(1.d-30,tot)
     tot = tot * sigmaE
  endif
end function atomhydrogenRayXsection

real(double) function molhydrogenRayXsection(lambda) result(tot)


  implicit none

! Dalgarno and Williams, 1962, ApJ, 136, 690

  real(double) :: lambda


  tot = 8.14d-13/lambda**4 + 1.28e-6/lambda**6 + 1.61d0/lambda**8

end function molhydrogenRayXsection




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

subroutine readKurucz(filename,nLines,lambda,kappa,excitationLower,excitationUpper,g)
  use utils_mod, only: wavenumbertoEv

  implicit none
  character(len=*) :: filename
  integer :: nLines
  real :: lambda(:),kappa(:),excitationLower(:),excitationUpper(:),g(:),junk
  integer :: i
  real, allocatable :: jLow(:)
  real, allocatable :: gf(:),alpha(:)
  real, parameter :: mHydrogen = 1.6733e-24      ! g
  real, parameter :: mElectron = 9.109565e-28    ! g
  real, parameter :: eCharge = 4.803242384E-10 !

  real, parameter :: cSpeed = 2.99792458e10  ! cm/s
  real, parameter :: pi = 3.141592654

  open(20,file=filename,status="old",form="formatted")

  read(20,*) nLines
  allocate(gf(1:nLines))
  allocate(alpha(1:nLines))
  allocate(jlow(1:nLines))



  if (writeoutput) write(*,*) "Reading from ",trim(filename)

  do i = 1, nLines

     if (writeoutput.and.(mod(i,nLines/10) == 0)) write(*,*) i 

     read(20,'(f10.4,f7.3,f5.1,f10.3,f5.1,f11.3)') &
          lambda(i), gf(i),  jLow(i), excitationLower(i), junk, excitationUpper(i)
  enddo

  close(20)

  lambda(1:nLines) = lambda(1:nLines) * 10.  ! nm to angstroms

  g(1:nLines) = (2.*real(Jlow(1:nLines))+1.)

  alpha(1:nLines) = ((pi*eCharge**2)/(mElectron*cSpeed)) * (10.e0**gf(1:nLines))

  kappa(1:nLines) = alpha(1:nLines) / (mHydrogen) ! divide by mass of H molecule

  kappa(1:nLines) = kappa(1:nLines) * 1.e10 ! to code units
  excitationUpper(1:Nlines) = abs(excitationUpper(1:Nlines)) ! get rid of negative flags
  excitationLower(1:Nlines) = abs(excitationLower(1:Nlines)) ! get rid of negative flags

  call wavenumberToEv(excitationLower, nLines)
  call wavenumberToEv(excitationUpper, nLines)

  write(*,*) "kappa ",kappa(1:nLines)
  if (writeoutput) write(*,*) "Done."

end subroutine readKurucz

subroutine createAllMolecularTables(nLam, lamArray, reset)
  implicit none
  integer :: nLam
  real(double) :: lamArray(:)
  character(len=80) :: dataDirectory
  logical, optional :: reset
  logical :: doReset

  doreset = .false.
  if (PRESENT(reset)) doReset = reset

  dataDirectory = " "


  nLookupTables = 8

  LookupTable(1)%molecule = "TiO"
  LookupTable(1)%source = "/data/dsa206/Abs_coeff/abs_coeff_TiO_Plez.nc"
  LookupTable(1)%filename = "tiolookuptable.dat"
  LookupTable(1)%abundance = 1.e-7
  lookupTable(1)%mu = 64.

  LookupTable(2)%molecule = "CH4"
  LookupTable(2)%source = "/data/dsa206/Abs_coeff/abs_coeff_CH4_YT10to10.nc"
  LookupTable(2)%filename = "ch4lookuptable.dat"
  LookupTable(2)%abundance = 1.e-8
  lookupTable(2)%mu = 16.

  LookupTable(3)%molecule = "H2-H2"
  LookupTable(3)%source = "/data/dsa206/Abs_coeff/abs_coeff_H2-H2_HITRAN.nc"
  LookupTable(3)%filename = "h2h2lookuptable.dat"
  LookupTable(3)%abundance = 1.
  lookupTable(3)%mu = 2.

  LookupTable(4)%molecule = "H2-He"
  LookupTable(4)%source = "/data/dsa206/Abs_coeff/abs_coeff_H2-He_HITRAN.nc"
  LookupTable(4)%filename = "h2helookuptable.dat"
  LookupTable(4)%abundance = 1.
  lookupTable(4)%mu = 2.

  LookupTable(5)%molecule = "H2O"
  LookupTable(5)%source = "/data/dsa206/Abs_coeff/abs_coeff_H2O_BT2.nc"
  LookupTable(5)%filename = "h2olookuptable.dat"
  LookupTable(5)%abundance = 1.e-8
  lookupTable(5)%mu = 18.

  LookupTable(6)%molecule = "NH3"
  LookupTable(6)%source = "/data/dsa206/Abs_coeff/abs_coeff_NH3_BYTe.nc"
  LookupTable(6)%filename = "nh3lookuptable.dat"
  LookupTable(6)%abundance = 1.e-8
  lookupTable(6)%mu = 17.

  LookupTable(7)%molecule = "CO"
  LookupTable(7)%source = "/data/dsa206/Abs_coeff/abs_coeff_CO_HITEMP.nc"
  LookupTable(7)%filename = "colookuptable.dat"
  LookupTable(7)%abundance = 1.e-5
  lookupTable(7)%mu = 28.

  LookupTable(8)%molecule = "VO"
  LookupTable(8)%source = "/data/dsa206/Abs_coeff/abs_coeff_VO_Plez.nc"
  LookupTable(8)%filename = "volookuptable.dat"
  LookupTable(8)%abundance = 1.e-8
  lookupTable(8)%mu = 67.

  call createKappaGrid(lookuptable, nLam, lamArray, doReset)
  

end subroutine createAllMolecularTables

subroutine readH2O(nLines,lambda,kappa,excitation,g)
  use utils_mod, only: wavenumbertoEv, locate, convertbyte2, convertbyte4


  implicit none
  integer :: nLines
  integer :: i, j, k
  real :: lambda(:), kappa(:), excitation(:), g(:)
  real, allocatable :: wavenumber(:),gf(:)
  real, parameter :: mHydrogen = 1.6733e-24      ! g
  real, parameter :: mElectron = 9.109565e-28    ! g
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

     if (writeoutput.and.(mod(i,nLines/10) == 0)) write(*,*) i 

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

     lambda(i) = real(wlvac*10.) ! nm to angstroms
     gf(i) = real(congf)
     wavenumber(i) = real(elo)
     excitation(i) = real(elo)
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

  if (Writeoutput) write(*,*) "Associating statistical weights..."
  do i = 1, nLines
     if ((mod(i,nLines/10) == 0).and.writeoutput) write(*,*) i 
     call locate(e, 273201, wavenumber(i),j)
     g(i) = wt(j)
  enddo

  kappa(1:nLines) = ((pi*eCharge**2)/(mElectron*cSpeed)) * gf(1:nLines) / g(1:nLines) &
       / (16.*mHydrogen + 2.*mHydrogen) ! divide by mass of h2o

  kappa(1:nLines) = kappa(1:nLines) * 1.e10 ! to code units

  call wavenumberToEv(excitation, nLines)

  if (writeoutput) write(*,*) "Done."

end subroutine readH2O

subroutine waitForFile(thisFile)
  character(len=*) :: thisFile
  integer :: i
  integer, parameter :: maxWait = 10000000

  i = 1
10 continue
  open(20,file=thisFile,status="old",form="unformatted",err=20)
  close(20)
  goto 666
20 continue
  i = i + 1
  if (i > maxWait) then
     write(*,*) "Tried to open file ",maxWait, " times"
     stop
  endif
  goto 10

666 continue
end subroutine waitForFile

#ifdef NETCDF
SUBROUTINE readAbsCoeffFileAmundsen(thisFile, thisKappaTable, nLambda, lamArray)

  type(MOLECULARKAPPAGRID) :: thisKappaTable
  integer :: nPTPairs, nnu
  integer :: ncidin
  integer, parameter :: ncvars = 30, ncdims = 30
  integer, dimension(ncvars) :: invarid
  integer, dimension(ncdims):: indimid
  real(double), allocatable :: nu(:), tArray(:), pArray(:)
  real(double), allocatable :: kAbs(:), tempkAbs(:)
  real(double) :: lamArray(:)
  real(double), allocatable :: lamKappa(:)
  integer :: nLambda
  character(len=*) :: thisFile
  character(len=256) :: dim_name
  integer :: i
  integer :: nTemps, it, iPT, i1, i2, j
  real(double) :: lam1, lam2, dLam, thisKap

  nTemps = 20
  iPT = 1
  call nf90_check(nf90_open(thisFile, NF90_NOWRITE, ncidin))
  call nf90_check(nf90_inq_dimid(ncidin, "pt_pair", indimid(1)))
  call nf90_check(nf90_inquire_dimension(ncidin, indimid(1), dim_name, nPTpairs))

  call nf90_check(nf90_inq_dimid(ncidin, "nu", indimid(2)))
  call nf90_check(nf90_inquire_dimension(ncidin, indimid(2), dim_name, nNu))
  allocate(pArray(1:nPTpairs), tArray(1:nPTpairs))
  allocate(nu(1:nNu))


  call nf90_check(nf90_inq_varid(ncidin, "nu", invarid(2)))
  call nf90_check(nf90_inq_varid(ncidin, "p_calc", invarid(3)))
  call nf90_check(nf90_inq_varid(ncidin, "t_calc", invarid(4)))
  call nf90_check(nf90_inq_varid(ncidin, "kabs", invarid(5)))



  call nf90_check(nf90_get_var(ncidin, invarid(2), nu))
  call nf90_check(nf90_get_var(ncidin, invarid(3), pArray))
  call nf90_check(nf90_get_var(ncidin, invarid(4), tArray))


  thisKappaTable%source = trim(thisFile)
  thisKappaTable%nLam = nLambda
  thisKappaTable%nTemps = nTemps
  allocate(lamKappa(1:nNu))
  do i = nNu, 1, -1
     lamKappa(i) = 1.d10/(nu(nNu-i+1))
  enddo
  allocate(thisKappaTable%tempArray(1:nTemps))
  allocate(thisKappaTable%lamArray(1:nLambda))
  allocate(thisKappaTable%kapArray(1:nTemps, 1:nLambda))
  thisKappaTable%lamArray(1:nLambda) = lamArray(1:nLambda)
  thisKappaTable%tempArray(1:nTemps) = tArray(1:nTemps)

  do it = 1, nTemps
     allocate(kAbs(1:nNu))
     allocate(tempkAbs(1:nNu))
     call nf90_check(nf90_get_var(ncidin, invarid(5), tempkAbs, start = (/ 1, it /), count=(/ nNu, 1 /)))
     do j = nNu, 1, -1
        kAbs(j) = tempKabs(nNu-j+1)
     enddo
     deallocate(tempKabs)
     do i = 1, nLambda
        if (i == 1) then 
           dLam = 0.5d0 * (lamArray(2)-lamArray(1))
           lam1 = lamArray(1) - dLam
           lam2 = lamArray(1) + dlam
        else if (i == nLambda) then
           dLam = 0.5d0*(lamArray(nLambda) - lamArray(nLambda-1))
           lam1 = lamArray(nLambda) - dlam
           lam2 = lamArray(nLambda) + dlam
        else
           lam1 = 0.5d0*(lamArray(i)+lamArray(i-1))
           lam2 = 0.5d0*(lamArray(i+1) + lamArray(i))
        endif
        call locate(lamKappa, nNu, lam1, i1)
        call locate(lamKappa, nNu, lam2, i2)
        thisKap = 0.d0
        do j = i1, i2
           thisKap = thisKap + kAbs(j)
        enddo
        thisKap = thisKap / dble(i2-i1+1)
        thisKappaTable%kapArray(it, i) = thisKap
     enddo
     deallocate(kAbs)
    enddo
    thisKappaTable%kapArray = thisKappaTable%kapArray * thisKappaTable%abundance * thisKappaTable%mu / 2.33d0
    thisKappaTable%kapArray = thisKappaTable%kapArray * 10.d0 * 1.d10 ! from m^2/kg to cm^2/g and then to code
  END SUBROUTINE readAbsCoeffFileAmundsen
#endif    

#ifdef NETCDF
  SUBROUTINE nf90_check(errcode)

    INTEGER,INTENT(IN) :: errcode

    IF (errcode == NF90_NOERR) THEN
       RETURN
    ELSE
       WRITE(*,*)'io_mod: NetCDF error ',nf90_strerror(errcode)
       STOP
    endif

  END SUBROUTINE nf90_check
#endif
  

end module gas_opacity_mod
