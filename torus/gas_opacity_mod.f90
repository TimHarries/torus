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


type(tioKappaGrid), save :: TioLookuptable

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

  kappa(1:nLines) = kappa(1:nLines) * 1.e-8 ! abundance
  

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

subroutine returnGasKappaValue(temperature, lambda, kappaAbs, kappaSca)
  real :: temperature
  real :: lambda
  real, optional :: kappaAbs, kappaSca
  real :: t1, t2
  integer :: i, j


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
  endif

  if (PRESENT(kappaSca)) then
     kappaSca = 1.e-30
  endif

end subroutine returnGasKappaValue
  
end module gas_opacity_mod
