#ifdef PHOTOION
!Photoionization utilities module
!Started 01/04/2011

!This module consolidates established functions and routines that are
!shared by both photoion_mod and photoionAMR_mod

! Structure Guide (quickSearchCode)
! 1. Recombination Rates (pu01)
! 2. Tables (pu02)
! 3. Charge Exchange Routines(pu03)
! 4. Population levels (??04)
! 5. Line emission (??05)
! 6. Continuum emission (??06)
! 7. Imaging (??07)
! 8. Wavelength array (??08)
! 9. Dust (??09)
! X. Unused routines (pu0X)


module photoion_utils_mod

use constants_mod
use messages_mod
use gridio_mod
use source_mod
use timing
use image_mod
use grid_mod
use amr_mod
use diffusion_mod
use mpi_amr_mod
use mpi_global_mod
use photon_mod
use phasematrix_mod
use ion_mod

implicit none

! - Shared type definitions

type SAHAMILNETABLE
   integer :: nFreq
   integer :: nTemp
   real(double), pointer :: temp(:)
   real(double), pointer :: freq(:)
   real(double), pointer :: Clyc(:, :)
   real(double), pointer :: emissivity(:)
end type SAHAMILNETABLE

type RECOMBTABLE
   integer :: nrho
   integer :: nTemp
   real(double), pointer :: rho(:)
   real(double),  pointer :: temp(:)
   real(double), pointer :: emissivity(:, :)
end type RECOMBTABLE

type GAMMATABLE
   integer :: nGamma
   integer :: nFreq
   integer :: nTemp
   real(double), pointer :: freq(:)
   real(double), pointer :: temp(:)
   real(double), pointer :: gamma(:,:)
end type GAMMATABLE

! If a subroutine/function needs to be used outside this module declare it as public here. 
public :: solvePops, refineLambdaArray, addRadioContinuumEmissivityMono, identifyForbiddenTransitionsInRange
private :: getCollisionalRates, identifyForbiddenTransition, addForbiddenToEmission, identifyRecombinationTransition, &
     addRecombinationToEmission, addForbiddenEmissionLine, addRecombinationEmissionLine, addRadioContinuumEmissivity
!     addRecombinationToEmission, addRadioContinuumEmissivity, addForbiddenEmissionLine, addRecombinationEmissionLine

contains

! 1.Recombination Rates (pu01)

! Radiative recombination rate for HI, HeI, HeII and CIV (P)
  real function vernerFerland(a, t, t0, t1, b)
    real(double):: a, t, t0, t1, b
    
    ! based on Verner and Ferland (1996) ApJS 103 467
    vernerFerland = real(a / (sqrt(t/t0) * (1.d0+sqrt(t/t0))**(1.d0-b) * (1.d0+sqrt(t/t1))**(1.d0+b)))
    
  end function vernerFerland


! Radiative recombination rate contribution for CI, CII, CIII,  NI, (P)
! NII, NIII, N IV, OI, OII, OIII,  NeII, NeIII and SI
  function nussbaumerStorey1983(t, a, b, c, d, f) result(rate)
    
    ! based on Nussbaumer & Storey 1983 AA 126 75
    real(double) :: t, a, b, c, d, f, rate
    
    if (t < 0.05d0) then
       rate = tiny(rate)
    else
       rate = (1.d-12) * ((a/t) + b + (c*t) + d*(t**2))*(t**(-1.5d0))*exp(-f/t)
    endif
    
  end function nussbaumerStorey1983


! Radiative recombination rate contribution for CI, CII, CIII, NI, (P)
! NII, NIII, NIV, OI, OII, OIII,  NeII, NeIII
  function ppb1991(z, a, b, c, d, temperature) result(rate)

    ! radiative recombination rates
    ! based on Pequignot, Petitjean & Boisson (1991) A&A 251 680
    
    real(double) :: z, a, b, c, d, temperature, t, rate
    
    t = 1.d-4 * temperature / z**2
    
    rate = 1.d-13 * z * (a*t**b)/(1.d0+c*t**d)
    
  end function ppb1991

! Radiative recombination rate for SI, SII, SIII, SIV (P)
  function svs1982(t, alpharad, xrad) result (rate)

    ! radiative recombination rates based on
    ! Shull and Van Steenberg, 1992, ApJS, 48, 95
    
    real(double) :: t, alpharad, xrad, rate
    
    rate = alpharad * (t /1.d4)**(-xrad)
  end function svs1982

!---awaiting comments                                  (F/P)
!Photoion_mod passes some unnecessary variables so although they 
!don't agree, it should use this one
  subroutine getRecombs(rates, thision, temperature)
    real(double) :: rates(:)
    type(IONTYPE) :: thisIon
    real(double) :: temperature
    
    select case(thisIon%species)
    case("O II")
       rates(1) = 0.d0
       rates(2) = rateFit(7.218d0, -0.575d0, temperature)
       rates(3) = rateFit(4.812d0, -0.575d0, temperature)
       rates(4) = rateFit(3.581d0, -0.495d0, temperature)
       rates(5) = rateFit(1.790d0, -0.495d0, temperature)
    case DEFAULT
       rates = 0.d0
    end select
  end subroutine getRecombs


!-- Awaiting comments
function rateFit(a,b,t) result (rate)
  real(double) :: a, b, t, rate

  rate = 1.d-13 * a * (t/1.d4)**b
end function rateFit


!   End of pu01

! 2. Tables (pu02)

!-- Awaiting comments                               (P)
subroutine createSahaMilneTables(hTable, heTable)
  use phfit_mod, only : phfit2
  type(SAHAMILNETABLE), intent(out) :: hTable, heTable
  real(double) :: nu0_h, nu0_he, nufinal_h, nufinal_he
  integer :: nFreq, nTemp
  integer :: i, j
  real :: e
  real(double) :: t1, t2
  real :: hxsec, hexsec
  real(double) :: dfreq, jnu

  nFreq = 10000
  nTemp = 100
  hTable%nFreq = nFreq
  heTable%nFreq = nFreq
  hTable%nTemp = nTemp
  heTable%nTemp = nTemp
  
  allocate(hTable%freq(1:nFreq), heTable%freq(1:nFreq))
  allocate(hTable%emissivity(1:nTemp), heTable%emissivity(1:ntemp))
  allocate(hTable%temp(1:ntemp), heTable%temp(1:ntemp))
  allocate(hTable%Clyc(1:nTemp,1:nFreq))
  allocate(heTable%Clyc(1:nTemp,1:nFreq))
  
  nu0_h = 13.6d0/ergtoev/hcgs
  nu0_he = 24.59d0/ergtoev/hcgs
  
  nufinal_h = 2.d0*nu0_h
  nufinal_he = 2.d0*nu0_he
  
  do i = 1, nFreq
     hTable%freq(i) = log10(nu0_h) + (log10(nuFinal_h)-log10(nu0_h))*dble(i-1)/dble(nFreq-1)
     heTable%freq(i) = log10(nu0_he) + (log10(nuFinal_he)-log10(nu0_he))*dble(i-1)/dble(nFreq-1)
  enddo
  hTable%freq(1:hTable%nFreq) = 10.d0**hTable%freq(1:hTable%nfreq)
  heTable%freq(1:hTable%nFreq) = 10.d0**heTable%freq(1:hTable%nfreq)

  t1 = 5000.d0
  t2 = 20000.d0

  do i = 1, nTemp
     hTable%temp(i) = t1 + (t2-t1)*dble(i-1)/dble(nTemp-1)
     heTable%temp(i) = t1 + (t2-t1)*dble(i-1)/dble(nTemp-1)
  enddo

  do i = 1, nTemp
     hTable%Clyc(i,1) = 0.d0
     heTable%Clyc(i,1) = 0.d0
     do j = 2, nFreq

        e = real(hTable%freq(j) * hcgs* ergtoev)
        call phfit2(1, 1, 1 , e , hxsec)

        dFreq = hTable%freq(j)-hTable%freq(j-1)
        jnu = ((hcgs*hTable%freq(j)**3)/(cSpeed**2)) * ((hcgs**2) /(twoPi*mElectron*Kerg*hTable%temp(i)))**(1.5d0) * &
             dble(hxsec/1.d10) *  exp(-hcgs*(hTable%freq(j)-nu0_h)/(kerg*hTable%temp(i)))
        hTable%Clyc(i,j) = hTable%Clyc(i,j-1) + jnu * dfreq
        
        e = real(heTable%freq(j) * hcgs * ergtoev)
        call phfit2(2, 2, 1 , e , hexsec)
        
        dFreq = heTable%freq(j)-heTable%freq(j-1)
        jnu = 2.d0*((hcgs*heTable%freq(j)**3)/(cSpeed**2)) * ((hCgs**2) / (twoPi*mElectron*kerg*heTable%temp(i)))**(1.5d0) * &
             dble(hexsec/1.d10) * exp(-hcgs*(heTable%freq(j)-nu0_he)/(kerg*heTable%temp(i)))
        heTable%Clyc(i,j) = heTable%Clyc(i,j-1) + jnu * dfreq
     enddo
     hTable%emissivity(i) = hTable%Clyc(i,hTable%nFreq)
     heTable%emissivity(i) = heTable%Clyc(i,heTable%nFreq)
     hTable%Clyc(i,1:hTable%nFreq) = hTable%Clyc(i,1:hTable%nFreq) / hTable%Clyc(i,hTable%nFreq)
     heTable%Clyc(i,1:heTable%nFreq) = heTable%Clyc(i,1:heTable%nFreq) / heTable%Clyc(i,heTable%nFreq)
  end do
  
  !  do j = 1, nFreq
  !     write(99  ,*) htable%freq(j),htable%clyc(50,j)
  !  enddo
  
end subroutine createSahaMilneTables


!-- Awaiting comments                                (P)
subroutine createRecombTable(table, tablefilename)
  type(RECOMBTABLE), intent(out) :: table
  character(len=*) :: tablefilename
  character(len=200) :: filename, datadirectory
  integer :: ia, ib, ne, ncut, i
  real :: e(1000)

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//tablefilename

  open(20,file=filename,status="old",form="formatted")
  read(20,*) table%nTemp, table%nrho

  allocate(table%rho(1:table%nRho))
  allocate(table%temp(1:table%ntemp))
  allocate(table%emissivity(1:table%ntemp,1:table%nrho))

  do ia=1,table%ntemp
     do ib=1,table%nrho
        read(20,'(1x,e10.3,5x,e10.3,13x,i2)') table%rho(ib),table%temp(ia),ncut
        ne=ncut*(ncut-1)/2
        read(20,'(8e10.3)') e(1:ne)
        table%emissivity(ia,ib) = SUM(e(1:ne-1))
     enddo
  enddo
  close(20)
end subroutine createRecombTable

!-- Awaiting comments                           (P)
subroutine createGammaTable(table, thisfilename)

! Ferland 1980 PASP 92 596

  type(GAMMATABLE), intent(out) :: table
  character(len=*) :: thisfilename
  character(len=200) :: dataDirectory, filename
  integer :: i

  call unixGetenv("TORUS_DATA", dataDirectory, i)
  filename = trim(dataDirectory)//"/"//thisfilename

  open(40, file=filename, form="formatted", status="old")
  read(40,*) table%nTemp, table%nFreq

  allocate(table%freq(1:table%nFreq))
  allocate(table%temp(1:table%nTemp))
  allocate(table%gamma(table%nFreq, table%nTemp))

  table%temp = 0.d0
  read(40,*)  table%temp(1:table%nTemp)

  do i = 1, table%nFreq
     read(40,*) table%freq(i), table%gamma(i,1:table%nTemp)
  enddo

  do i = 1, table%nFreq
     table%freq(i) = table%freq(i) * nuHydrogen
  enddo

  where (table%gamma(1:table%nFreq,1:table%nTemp) == 0.d0)
     table%gamma(1:table%nFreq,1:table%nTemp) = 1.d-10
  end where
  table%gamma(1:table%nFreq,1:table%nTemp) =log10(table%gamma(1:table%nFreq,1:table%nTemp))
  table%freq(1:table%nFreq) = log10(table%freq(1:table%nFreq))
  table%temp(1:table%nTemp) = log10(table%temp(1:table%nTemp))
  close(40)
end subroutine createGammaTable

!   End of pu02

! 3. Charge Exchange Routines (pu03)

!--Awaiting Comments                                         (U)
subroutine getChargeExchangeRecomb(parentIon, temperature, nhi, recombRate)
  type(IONTYPE) :: parentIon
  real(double) :: nhi
  real(double), intent(out) :: recombRate
  real :: t4, a, b, c, d
  real :: temperature
  
  recombRate  = 0.d0
  
  select case(parentIon%z)
  case(7)
     select case(parentIon%n)
     case(4) ! N IV
        t4 = temperature / 1.e4
        a = 3.05e-10
        b = 0.60
        c = 2.65
        d = -0.93
        recombRate  = kingdonFerland96(t4, a, b, c, d)
     end select
  case(8)
     select case(parentIon%n)
     case(7) ! O II
        t4 = temperature / 1.e4
        a = 1.04e-9
        b = 3.15e-2
        c = -0.61
        d = -9.73
        recombRate  = kingdonFerland96(t4, a, b, c, d)
     case(6) ! OIII
        t4 = temperature / 1.e4
        a = 1.04e-9
        b = 0.27
        c = 2.02
        d = -5.92
        recombRate  = kingdonFerland96(t4, a, b, c, d)
     end select
  end select
  
  recombRate = recombRate * nhi
  
end subroutine getChargeExchangeRecomb


!-- Awaiting Comments                                             (U)
subroutine getChargeExchangeIon(parentIon, temperature,  nHII, IonRate)
  type(IONTYPE) :: parentIon
  real(double) :: nHII
  real(double), intent(out) :: ionRate
  real :: temperature
  real :: t4, a, b, c, d

  IonRate = 0.d0

  select case(parentIon%z)
  case(7)
     select case(parentIon%n)
     case(7) ! N I charge exchange Ionization
        t4 = temperature / 1.e4
        a = 4.55e-12
        b = -0.29
        c = -0.92
        d = -8.38
        ionRate  = kingdonFerland96(t4, a, b, c, d)
     end select
  case(8)
     select case(parentIon%n)
     case(8) ! O I charge exchange ionization
        t4 = temperature / 1.e4
        a = 7.40e-11
        b = 0.47
        c = 24.37
        d = -0.74
        ionRate  = kingdonFerland96(t4, a, b, c, d)
     end select
  end select
  
  ionRate = ionRate * nhii
end subroutine getChargeExchangeIon
      

!--Awaiting comments                                (U)
function kingdonFerland96(t4, a, b, c, d) result (alpha)
  real :: alpha
  real :: t4, a, b, c, d
  alpha = a*(t4**b)*(1.+c*exp(d*t4))
end function kingdonFerland96

! End of pu03

! ??04 Moved by DMA. Awaiting comments. 

subroutine solvePops(thisIon, pops, ne, temperature, debug)
  type(IONTYPE) :: thisIon
  real(double) :: ne
  real, intent(out) :: pops(:)
  real :: temperature
  real(double), allocatable :: matrixA(:,:), MatrixB(:), tempMatrix(:,:), qeff(:,:),  rates(:)
  integer :: n, iTrans, i, j
  real :: excitation, deexcitation, arateji
  logical, optional :: debug
  excitation  = 0.; deexcitation = 0.
  n = thisIon%nLevels
  allocate(matrixA(1:n, 1:n), matrixB(1:n), tempMatrix(1:n,1:n), qeff(1:n,1:n), rates(1:n))

  matrixA = 1.d-30
  matrixB = 0.d0

  call getRecombs(rates, thision, dble(temperature))

  matrixA(1,:) = 1.d0
  matrixB(1) = 1.d0

  matrixB(2:n) = 0.d0 !rates(2:n) * ne * (nh * ionFrac * thisIon%abundance)
  do iTrans = 1, thisIon%nTransitions
     i = thision%transition(itrans)%i
     j = thision%Transition(itrans)%j
     call getCollisionalRates(thisIon, iTrans, temperature, excitation, deexcitation)
     qeff(i,j) = excitation
     qeff(j,i) = deexcitation
  enddo

  do i = 2, n
     do j = 1, n

        do iTrans = 1, thisIon%nTransitions
           if (((i == thision%transition(itrans)%i).and.(j == thision%Transition(itrans)%j)).or. &
               ((i == thision%transition(itrans)%j).and.(j == thision%Transition(itrans)%i))) then
              arateji =  thisIon%transition(iTrans)%a
              
              matrixA(i,j) = matrixA(i,j) + ne * qeff(j, i)
              matrixA(i,i) = matrixA(i,i) - ne * qeff(i, j)
              if (j > i) then
                 matrixA(i,j) = matrixA(i,j) + arateji
              else
                 matrixA(i,i) = matrixA(i,i) - arateji
              endif

           endif
        enddo

     enddo
  enddo

  tempMatrix = matrixA

  if (PRESENT(debug)) then
     if (debug) then
        do i = 1, n
           write(*,'(1p,9e12.3)') tempmatrix(i,1:n)
        enddo
     endif
  endif
  
  call luSlv(matrixA, matrixB)

  matrixB(1:n) = matrixB(1:n) / SUM(matrixB(1:n))

  do i = 1, n
     pops(i) = real(max(1.d-30,matrixB(i)))
  enddo

  deallocate(matrixA, matrixB, tempMatrix, qeff, rates)

end subroutine solvePops

subroutine getCollisionalRates(thisIon, iTransition, temperature, excitation, deexcitation)
  real, intent(out) :: excitation, deexcitation
  type(IONTYPE), intent(in) :: thisIon
  integer, intent(in) :: iTransition
  real, intent(in) :: temperature
  integer :: i
  real :: thisGamma
  real :: t , fac, maxTemp
  real :: boltzFac
  logical, save :: firstTime = .true.
  !$OMP THREADPRIVATE (firstTime)

 ! - thap -check to see if the temperature is too hot for the gamma table                                                
  maxTemp = maxval(thisIon%transition(iTransition)%t)

  if(temperature > maxTemp) then
    if(firstTime) then
        firstTime = .false.
        write(*,*) "Warning: Temperature Exceeds Scope Of Atomic Data"
        write(*,*) "         Proceeding With The Maximum Known Value..."
     end if

     t = maxTemp

  else
     t = max(min(real(thisIon%transition(iTransition)%t(thisIon%transition(iTransition)%ngamma)),temperature), &
          real(thisIon%transition(iTransition)%t(1)))
  endif

  call locate(thisIon%transition(iTransition)%t, thisIon%transition(iTransition)%ngamma, t, i)
  fac = (t - thisIon%transition(iTransition)%t(i))/(thisIon%transition(iTransition)%t(i+1) - thisIon%transition(iTransition)%t(i))
  thisGamma = thisIon%transition(iTransition)%gamma(i) + &
       fac * (thisIon%transition(iTransition)%gamma(i+1) - thisIon%transition(iTransition)%gamma(i))

  boltzFac =  real(exp(-thisIon%transition(iTransition)%energy / (kev*temperature)))

  fac = (8.63e-6 / sqrt(temperature)) * thisGamma

  deexcitation =  fac / thisIon%level(thisIon%transition(iTransition)%j)%g

  excitation =  fac / thisIon%level(thisIon%transition(iTransition)%i)%g * boltzFac

end subroutine getCollisionalRates

! ??05 Line emission routines

!
! Forbidden emission lines
!

subroutine addForbiddenEmissionLine(grid, dLambda, lineWavelength)
  type(GRIDTYPE) :: grid
  real(double), intent(in) :: lineWavelength, dLambda
  integer :: iIon, iTransition
  logical :: ok

  call identifyForbiddenTransition(grid, lineWavelength, iIon, iTransition, ok)
  if (ok) call addForbiddenToEmission(grid, grid%octreeRoot, dLambda, iIon, iTransition) 

end subroutine addForbiddenEmissionLine

!for SED/image generation over multiple lines
subroutine identifyForbiddenTransitionsInRange(grid, lamStart, lamEnd, transitionLams)

  type(GRIDTYPE),intent(in) :: grid
  real(double),intent(in) :: lamStart, lamEnd
  integer :: i, j, k
  real(double), intent(out) :: transitionLams(5000)

  transitionLams = 0.d0
  k = 1
  do i = 1, grid%nIon
     do j = 1, grid%ion(i)%nTransitions
        if (grid%ion(i)%transition(j)%lambda < lamEnd .and. &
             grid%ion(i)%transition(j)%lambda > lamStart) then
           
           transitionLams(k) = grid%ion(i)%transition(j)%lambda
           k = k + 1
           
        end if
     end do
  enddo

end subroutine identifyForbiddenTransitionsInRange

subroutine identifyForbiddenTransition(grid, lambda, iIon, iTransition, ok)
  type(GRIDTYPE),intent(in) :: grid
  real(double),intent(in) :: lambda
  integer, intent(out) :: iIon, iTransition
  logical, intent(inout) :: ok
  integer :: i, j
  character(len=80) :: message

  ok = .false.

  do i = 1, grid%nIon
     do j = 1, grid%ion(i)%nTransitions
        if (abs(lambda-grid%ion(i)%transition(j)%lambda) < 1.d0) then
           write(message, '(a, a, a)') "Transition of ", trim(grid%ion(i)%species), " identified."
           call writeInfo(message)
           iIon = i
           iTransition = j
           ok = .true.
        end if
     end do
  enddo
end subroutine identifyForbiddenTransition

recursive subroutine addForbiddenToEmission(grid, thisOctal, dLambda, iIon, iTransition)
  type(GRIDTYPE) :: grid
  real(double) :: dLambda, rate
  real :: pops(10)
  integer :: iIon, iTransition
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  integer :: subcell, i

  do subcell = 1, thisOctal%maxChildren
     if (thisOctal%hasChild(subcell)) then
        ! find the child
        do i = 1, thisOctal%nChildren, 1
           if (thisOctal%indexChild(i) == subcell) then
              child => thisOctal%child(i)
              call addForbiddenToEmission(grid, child, dLambda, iIon, iTransition)
              exit
           end if
        end do
     else

        call solvePops(grid%ion(iIon), pops, thisOctal%ne(subcell), thisOctal%temperature(subcell))
        rate =  pops(grid%ion(iion)%transition(iTransition)%j) * grid%ion(iion)%transition(itransition)%energy * &
             grid%ion(iion)%transition(itransition)%a/ergtoev
        rate = rate * grid%ion(iion)%abundance * thisOctal%nh(subcell) * thisOctal%ionFrac(subcell, iion)

        thisOctal%etaCont(subcell) =  thisOctal%etaCont(subcell) + rate/dLambda

     endif
  enddo
end subroutine addForbiddenToEmission

!
! Recombination emission lines
!

subroutine addRecombinationEmissionLine(grid, dLambda, lineWavelength)
  type(GRIDTYPE) :: grid
  real(double), intent(in) :: lineWavelength, dLambda
  integer :: ilower, iupper
  logical :: ok

  call identifyRecombinationTransition(lineWavelength, ilower, iupper, ok)
  if (ok) call addRecombinationToEmission(grid, grid%octreeRoot, dLambda, ilower, iupper)
   
end subroutine addRecombinationEmissionLine

  subroutine identifyRecombinationTransition(lambdaLine, ilower, iupper, ok)
    real(double),intent(in) :: lambdaLine
    integer, intent(out) :: ilower, iupper
    logical, intent(inout) :: ok
    integer :: iup, ilow
     real(double) :: lambdaTrans(20,20) = reshape( source=&
          (/ 000000.d-8,1215.67d-8,1025.72d-8,992.537d-8,949.743d-8,937.803d-8,930.748d-8,926.226d-8,923.150d-8,920.963d00,&
          919.352d-8,918.129d-8,917.181d-8,916.429d-8,915.824d-8,915.329d-8,914.919d-8,914.576d-8,914.286d-8,914.039d-8,&  
          0.00000d-8,0000000d-8,6562.80d-8,4861.32d-8,4340.36d-8,4101.73d-8,3970.07d-8,3889.05d-8,3835.38d-8,3797.90d-8,&
          3770.63d-8,3750.15d-8,3734.37d-8,3721.94d-8,3711.97d-8,3703.85d-8,3697.15d-8,3691.55d-8,3686.83d-8,3682.81d-8,&  
          0.00000d-8,0.00000d-8,0000000d-8,18751.0d-8,12818.1d-8,10938.1d-8,10049.4d-8,9545.98d-8,9229.02d-8,9014.91d-8,&
          8862.79d-8,8750.47d-8,8665.02d-8,8598.39d-8,8545.39d-8,8502.49d-8,8467.26d-8,8437.96d-8,8413.32d-8,8392.40d-8,&  
          0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,40512.0d-8,26252.0d-8,21655.0d-8,19445.6d-8,18174.1d-8,17362.1d-8,&
          16806.5d-8,16407.2d-8,16109.3d-8,15880.5d-8,15700.7d-8,15556.5d-8,15438.9d-8,15341.8d-8,15260.6d-8,15191.8d-8,&  
          0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,74578.0d-8,46525.0d-8,37395.0d-8,32961.0d-8,30384.0d-8,&
          28722.0d-8,27575.0d-8,26744.0d-8,26119.0d-8,25636.0d-8,25254.0d-8,24946.0d-8,24693.0d-8,24483.0d-8,24307.0d-8,&  
          0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,123680.d-8,75005.0d-8,59066.0d-8,51273.0d-8,&
          46712.0d-8,43753.0d-8,41697.0d-8,40198.0d-8,39065.0d-8,38184.0d-8,37484.0d-8,36916.0d-8,36449.0d-8,36060.0d-8,&  
          0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,190570.d-8,113060.d-8,87577.0d-8,&
          75061.0d-8,67701.0d-8,62902.0d-8,59552.0d-8,57099.0d-8,55237.0d-8,53783.0d-8,52622.5d-8,51679.0d-8,50899.0d-8,&  
          0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,0000000d-8,277960.d-8,162050.d-8,&
          123840.d-8,105010.d-8,93894.0d-8,86621.0d-8,81527.0d-8,77782.0d-8,74930.0d-8,72696.0d-8,70908.0d-8,69448.0d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,388590.d-8,&  
          223340.d-8,168760.d-8,141790.d-8,125840.d-8,115360.d-8,108010.d-8,102580.d-8,98443.0d-8,95191.0d-8,92579.0d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          525200.d-8,298310.d-8,223250.d-8,186100.d-8,164070.d-8,149580.d-8,139380.d-8,131840.d-8,126080.d-8,121530.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,690500.d-8,388320.d-8,288230.d-8,238620.d-8,209150.d-8,189730.d-8,176030.d-8,165900.d-8,158120.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,000000.d-8,887300.d-8,494740.d-8,364610.d-8,300020.d-8,261610.d-8,236260.d-8,218360.d-8,205090.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,111800.d-7,619000.d-8,453290.d-8,371000.d-8,322000.d-8,289640.d-8,266740.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,138600.d-7,762300.d-8,555200.d-8,452220.d-8,390880.d-8,350300.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,169400.d-7,926100.d-8,671200.d-8,544400.d-8,468760.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,204400.d-7,111200.d-7,802300.d-8,648200.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,243800.d-7,132100.d-7,949200.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,288200.d-7,155400.d-7,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
          000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,337400.d-7,&  
          0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /), shape=(/20,20/))

     ok = .false.
     do iup = 15, 3, -1
        do ilow = 2, min0(8, iup-1)
           if (abs(lambdaTrans(iup, ilow)*1.d8-lambdaLine)/lambdaLine < 0.01d0) then
              iLower = ilow
              iUpper = iUp
              ok = .true.
           endif
        enddo
     enddo
     if (.not.ok) then
        call writeFatal("Recombination line not identified")
        stop
     endif


  end subroutine identifyRecombinationTransition

   recursive subroutine addRecombinationToEmission(grid, thisOctal, dLambda, ilower, iupper)
     type(GRIDTYPE) :: grid
     real(double) :: dLambda
     integer :: iLower, iUpper
     type(octal), pointer   :: thisOctal
     type(octal), pointer  :: child 
     integer :: subcell, i
     real :: emissivity(3:15,2:8)
     real(double) :: hBeta, LineEmissivity

     emissivity(15, 2:8) = (/ 0.156E-01, 0.541E-02, 0.266E-02, 0.154E-02, 0.974E-03, 0.657E-03, 0.465E-03 /)
     emissivity(14, 2:8) = (/ 0.192E-01, 0.666E-02, 0.328E-02, 0.190E-02, 0.120E-02, 0.811E-03, 0.573E-03 /)
     emissivity(13, 2:8) = (/ 0.240E-01, 0.832E-02, 0.411E-02, 0.237E-02, 0.151E-02, 0.102E-02, 0.719E-03 /)
     emissivity(12, 2:8) = (/ 0.305E-01, 0.106E-01, 0.524E-02, 0.303E-02, 0.193E-02, 0.130E-02, 0.917E-03 /)
     emissivity(11, 2:8) = (/ 0.397E-01, 0.138E-01, 0.683E-02, 0.396E-02, 0.252E-02, 0.170E-02, 0.119E-02 /)
     emissivity(10, 2:8) = (/ 0.530E-01, 0.184E-01, 0.914E-02, 0.531E-02, 0.338E-02, 0.228E-02, 0.158E-02 /)
     emissivity(9, 2:8) = (/ 0.731E-01, 0.254E-01, 0.127E-01, 0.737E-02, 0.469E-02, 0.312E-02, 0.204E-02 /)
     emissivity(8,2:8) = (/ 0.105E+00, 0.366E-01, 0.183E-01, 0.107E-01, 0.673E-02, 0.425E-02, 0.000E+00 /)
     emissivity(7,2:8) = (/ 0.159E+00, 0.555E-01, 0.278E-01, 0.162E-01, 0.976E-02, 0.000E+00, 0.000E+00 /)
     emissivity(6,2:8) = (/ 0.259E+00, 0.904E-01, 0.455E-01, 0.255E-01, 0.000E+00, 0.000E+00, 0.000E+00 /)
     emissivity(5,2:8) = (/ 0.468E+00, 0.163E+00, 0.802E-01, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
     emissivity(4,2:8) = (/ 0.100E+01, 0.339E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)
     emissivity(3,2:8) = (/ 0.286E+01, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00 /)

     do subcell = 1, thisOctal%maxChildren
        if (thisOctal%hasChild(subcell)) then
           ! find the child
           do i = 1, thisOctal%nChildren, 1
              if (thisOctal%indexChild(i) == subcell) then
                 child => thisOctal%child(i)
                 call addRecombinationToEmission(grid, child, dLambda, iLower, iUpper)
                 exit
              end if
           end do
        else

           if (thisOctal%temperature(subcell) /= 0.d0) then

              Hbeta = 10.d0**(-0.870d0*log10(thisOctal%temperature(subcell)) + 3.57d0)
              Hbeta = Hbeta * thisOctal%ne(subcell) *  thisOctal%nh(subcell) * &
                   thisOctal%ionFrac(subcell, 2) * grid%ion(1)%abundance
              
              ! calculate emission due to HI recombination lines [e-25 ergs/s/cm^3]
              lineEmissivity = emissivity(iupper, ilower) * Hbeta * 1.e-25
              
              thisOctal%etaCont(subcell) =  thisOctal%etaCont(subcell) + lineEmissivity/dLambda
           endif 
        endif
     enddo
   end subroutine addRecombinationToEmission

! ??06 Continuum emission routines

  recursive subroutine addRadioContinuumEmissivity(thisOctal)

    use stateq_mod, only : alpkk
    type(octal), pointer  :: thisOctal
    type(octal), pointer  :: child 
    integer               :: subcell
    integer               :: i
    real(double) :: eta, freq
    
    do subcell = 1, thisOctal%maxChildren, 1

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call addRadioContinuumEmissivity(child)
                exit
             end if
          end do            
       else
          freq = cspeed / (20.d0) ! 20 cm radio free-free
          eta =  thisOctal%Ne(subcell)**2 * &
               alpkk(freq,real(thisOctal%temperature(subcell),kind=db))* &
               exp(-(hcgs*freq)/(kerg*thisOctal%temperature(subcell)))
          
          eta=eta*real((2.0*dble(hcgs)*dble(freq)**3)/(dble(cspeed)**2))
             
          thisOctal%etaCont(subcell) = eta !* 1.d10
          thisOctal%biasCont3d(subcell) = 1.d0
       end if

    end do

  end subroutine addRadioContinuumEmissivity

  recursive subroutine addRadioContinuumEmissivityMono(thisOctal, lambda, stack)

    use stateq_mod, only : alpkk
    type(octal), pointer  :: thisOctal
    type(octal), pointer  :: child 
    integer               :: subcell
    integer               :: i
    real(double) :: eta, freq
    real :: lambda
    logical, intent(in) :: stack
    
    do subcell = 1, thisOctal%maxChildren, 1

       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call addRadioContinuumEmissivityMono(child, lambda, stack)
                exit
             end if
          end do            
       else
          freq = cspeed / (lambda*angstromToCm)
          eta =  thisOctal%Ne(subcell)**2 * &
               alpkk(freq,real(thisOctal%temperature(subcell),kind=db))* &
               exp(-(hcgs*freq)/(kerg*thisOctal%temperature(subcell)))
          
          eta=eta*real((2.0*dble(hcgs)*dble(freq)**3)/(dble(cspeed)**2))
          
          if(stack) then
             thisOctal%etaCont(subcell) = thisOctal%etaCont(subcell) + (eta)!*1.d10)
          else
             thisOctal%etaCont(subcell) = eta! * 1.d10
          end if
          thisOctal%biasCont3d(subcell) = 1.d0
       end if

    end do

  end subroutine addRadioContinuumEmissivityMono

! ??07 Imaging
subroutine setupGridForImage(grid, outputimageType, lambdaImage, iLambdaPhoton, nsource, source, lcore)
  use lucy_mod, only: calcContinuumEmissivityLucyMono

  implicit none
  
  type(GRIDTYPE) :: grid
  character(len=*), intent(in) :: outputImageType
  real, intent(in)             :: lambdaImage
  integer, intent(out)         :: iLambdaPhoton
  integer, intent(in)          :: nsource
  type(SOURCETYPE), intent(in) :: source(:)
  real(double), intent(out)    :: lCore
  character(len=80) :: message

  write(message,*) "Setting up grid for "//trim(outputimageType)//" image."
  call writeInfo(message, FORINFO)

  select case (outputimageType)

  case("freefree")

     ilambdaPhoton = grid%nLambda
     call  addRadioContinuumEmissivity(grid%octreeRoot)
     lcore = tiny(lcore)

  case("forbidden")

     call locate(grid%lamArray, grid%nlambda, real(lambdaImage), ilambdaPhoton)
     call calcContinuumEmissivityLucyMono(grid, grid%octreeRoot, grid%lamArray, lambdaImage, iLambdaPhoton)
     call addForbiddenEmissionLine(grid, 1.d0, dble(lambdaImage))

     if (nSource > 0) then              
!THAW - need to add the grid to allow for off-grid sources?
        lCore = sumSourceLuminosityMonochromatic(grid, source, nsource, dble(grid%lamArray(iLambdaPhoton)))
     else
        lcore = tiny(lcore)
     endif

  case("recombination")

     call locate(grid%lamArray, grid%nlambda, real(lambdaImage), ilambdaPhoton)
     call calcContinuumEmissivityLucyMono(grid, grid%octreeRoot, grid%lamArray, lambdaImage, iLambdaPhoton)
     call addRecombinationEmissionLine(grid, 1.d0, dble(lambdaImage))
     
     if (nSource > 0) then              
        lCore = sumSourceLuminosityMonochromatic(grid, source, nsource, dble(grid%lamArray(iLambdaPhoton)))
     else
        lcore = tiny(lcore)
     endif

  case("dustonly")

     write(message,"(a,f8.2,a)") "Image wavelength: ", lambdaImage*angsToMicrons, " microns"
     call writeInfo(message, FORINFO)
     call locate(grid%lamArray, grid%nlambda, real(lambdaImage), ilambdaPhoton)
     call calcContinuumEmissivityLucyMono(grid, grid%octreeRoot, grid%lamArray, lambdaImage,iLambdaPhoton)
     
     if (nSource > 0) then              
        lCore = sumSourceLuminosityMonochromatic(grid, source, nsource, dble(grid%lamArray(iLambdaPhoton)))
     else
        lcore = tiny(lcore)
     endif

  case DEFAULT

     call writeFatal("Imagetype "//trim(outputimageType)//" not recognised")
     stop

  end select

end subroutine setupGridForImage

! 08 refineLambdaArray adds points to the wavelength array 

  subroutine refineLambdaArray(lamArray, nLambda, grid)
    use utils_mod, only: insertBin
    type(GRIDTYPE) :: grid
    real :: lamArray(:)
    integer :: nLambda
    integer :: i, j
    integer :: iup, ilow

  real(double) :: lambdaTrans(20,20) = reshape( source=&
    (/ 000000.d-8,1215.67d-8,1025.72d-8,992.537d-8,949.743d-8,937.803d-8,930.748d-8,926.226d-8,923.150d-8,920.963d00,&
       919.352d-8,918.129d-8,917.181d-8,916.429d-8,915.824d-8,915.329d-8,914.919d-8,914.576d-8,914.286d-8,914.039d-8,&  
       0.00000d-8,0000000d-8,6562.80d-8,4861.32d-8,4340.36d-8,4101.73d-8,3970.07d-8,3889.05d-8,3835.38d-8,3797.90d-8,&
       3770.63d-8,3750.15d-8,3734.37d-8,3721.94d-8,3711.97d-8,3703.85d-8,3697.15d-8,3691.55d-8,3686.83d-8,3682.81d-8,&  
       0.00000d-8,0.00000d-8,0000000d-8,18751.0d-8,12818.1d-8,10938.1d-8,10049.4d-8,9545.98d-8,9229.02d-8,9014.91d-8,&
       8862.79d-8,8750.47d-8,8665.02d-8,8598.39d-8,8545.39d-8,8502.49d-8,8467.26d-8,8437.96d-8,8413.32d-8,8392.40d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,40512.0d-8,26252.0d-8,21655.0d-8,19445.6d-8,18174.1d-8,17362.1d-8,&
       16806.5d-8,16407.2d-8,16109.3d-8,15880.5d-8,15700.7d-8,15556.5d-8,15438.9d-8,15341.8d-8,15260.6d-8,15191.8d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0000000d-8,74578.0d-8,46525.0d-8,37395.0d-8,32961.0d-8,30384.0d-8,&
       28722.0d-8,27575.0d-8,26744.0d-8,26119.0d-8,25636.0d-8,25254.0d-8,24946.0d-8,24693.0d-8,24483.0d-8,24307.0d-8,&  
       0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,0.00000d-8,123680.d-8,75005.0d-8,59066.0d-8,51273.0d-8,&
       46712.0d-8,43753.0d-8,41697.0d-8,40198.0d-8,39065.0d-8,38184.0d-8,37484.0d-8,36916.0d-8,36449.0d-8,36060.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,190570.d-8,113060.d-8,87577.0d-8,&
       75061.0d-8,67701.0d-8,62902.0d-8,59552.0d-8,57099.0d-8,55237.0d-8,53783.0d-8,52622.5d-8,51679.0d-8,50899.0d-8,&  
       0.00000d-8,0000000d-8,0000000d-8,0000000d-8,0000000d-8,0.00000d-8,0000000d-8,0000000d-8,277960.d-8,162050.d-8,&
       123840.d-8,105010.d-8,93894.0d-8,86621.0d-8,81527.0d-8,77782.0d-8,74930.0d-8,72696.0d-8,70908.0d-8,69448.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,388590.d-8,&  
       223340.d-8,168760.d-8,141790.d-8,125840.d-8,115360.d-8,108010.d-8,102580.d-8,98443.0d-8,95191.0d-8,92579.0d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       525200.d-8,298310.d-8,223250.d-8,186100.d-8,164070.d-8,149580.d-8,139380.d-8,131840.d-8,126080.d-8,121530.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,690500.d-8,388320.d-8,288230.d-8,238620.d-8,209150.d-8,189730.d-8,176030.d-8,165900.d-8,158120.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,887300.d-8,494740.d-8,364610.d-8,300020.d-8,261610.d-8,236260.d-8,218360.d-8,205090.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,111800.d-7,619000.d-8,453290.d-8,371000.d-8,322000.d-8,289640.d-8,266740.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,138600.d-7,762300.d-8,555200.d-8,452220.d-8,390880.d-8,350300.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,169400.d-7,926100.d-8,671200.d-8,544400.d-8,468760.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,204400.d-7,111200.d-7,802300.d-8,648200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,243800.d-7,132100.d-7,949200.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,288200.d-7,155400.d-7,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,&  
       000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,000000.d-8,337400.d-7,&  
       0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /), shape=(/20,20/))


  do iup = 15, 3, -1
     do ilow = 2, min0(8, iup-1)
        call insertBin(lamArray, nLambda, real(lambdaTrans(iup, ilow)*1.e8), 1.)
    enddo
 enddo
 call insertBin(lamArray, nLambda, 1215.67, 1.)


! refine for forbidden line transitions

    do i = 1 , grid%nIon
       do j = 1, grid%ion(i)%nTransitions
          call insertBin(lamArray, nLambda, &
               real(grid%ion(i)%transition(j)%lambda), 1.)
       enddo
    enddo

     
  end subroutine refineLambdaArray

!   pu09 Dust

recursive subroutine quickSublimate(thisOctal, fraction)
  use inputs_mod, only : hOnly
  type(octal), pointer   :: thisOctal
  type(octal), pointer  :: child 
  ! Where dust is present set dustTypeFraction to this value. 
  real, optional, intent(in) :: fraction
  integer :: subcell, i
  
  do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                if (present(fraction)) then 
                   call quickSublimate(child, fraction)
                else
                   call quickSublimate(child)
                endif
                exit
             end if
          end do
       else


          if (thisOctal%temperature(subcell) > 1500.) then
             thisOctal%dustTypeFraction(subcell,:) = 1.d-20
          else
             if (present(fraction)) then
                thisOctal%dustTypeFraction(subcell,:) = fraction
             else
                thisOctal%dustTypeFraction(subcell,:) = 1.d0
             endif
          endif
!
!          if(.not. hOnly) then
!             if ((thisOctal%ionFrac(subcell,1) < 0.1).or.(thisOctal%temperature(subcell) > 1500.)) then
!                thisOctal%dustTypeFraction(subcell,:) = 1.d-20
!             else
!                if (present(fraction)) then 
!                   thisOctal%dustTypeFraction(subcell,:) = fraction
!                else
!                   thisOctal%dustTypeFraction(subcell,:) = 1.d0
!                endif
!                thisOctal%ne(subcell) = tiny(thisOctal%ne(subcell))
!1                thisOctal%ionFrac(subcell,1) = 1.d0
!                thisOctal%ionFrac(subcell,2) = tiny(thisOctal%ionFrac(subcell,2))
!             endif
!          Else
!             thisOctal%dustTypeFraction(subcell,:) = 1.d-20
!          end if

    endif
    enddo
  end subroutine quickSublimate

!   pu0X Unused routines
      
! Currently never called. HI recombination rates (P)
function hRecombination(temperature) result (rate)
  real :: temperature
  real(double) :: rate, t0, t1, b
  
  t0 = sqrt(temperature/3.148d0)
  t1 = sqrt(temperature/7.036d5)
  b = 0.7480
  
  rate = 7.982d-11 / (t0*(1.d0+t0)**(1.d0-b) * (1.d0 + t1)**(1.d0+b))
end function hRecombination

! Currently never called. He recombination rates (P)
subroutine calcHeRecombs(te, alpha1, alpha21s, alpha21p, alpha23s)
  real :: te, t
  real(double) :: alpha1, alpha21s, alpha21p, alpha23s
  
  t = te / 1.e4
  
  alpha1 = 1.54e-13 * t**(-0.486)
  
  alpha21s = 2.06e-14 * t**(-0.676)
  
  alpha21p = 4.17e-14 * t**(-0.861)
  
  alpha23s = 2.10e-13 * t**(-0.778)
end subroutine calcHeRecombs
  
  
! Currently never called. Recomb rate into ground state. (P)
function recombToGround(temperature) result (alpha1)
  real :: temperature
  real(double) :: alpha1
  
  alpha1 = 1.58d-13 * (temperature/1.d4)**(-0.53d0)  ! kenny's photo paper equation 24
end function recombToGround




!   End of pu0X
  
end module photoion_utils_mod

#endif
