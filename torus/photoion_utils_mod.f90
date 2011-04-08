!Photoionization utilities module
!Started 01/04/2011

!This module consolidates established functions and routines that are
!shared by both photoion_mod and photoionAMR_mod

! Structure Guide (quickSearchCode)
! 1. Recombination Rates (pu01)
! 2. Tables (pu02)
! 3. Charge Exchange Routines(pu03)
! X. Unused routines (pu0X)

!NOTE: Until this module is well established, routines will have either a P,
!      a U or an F next to their leading comments. These stand for:
!      P - passed: photoionAMR_mod & photoion_mod match for this routine
!      F - failed: Pass criterion unfulfilled (rm or modify this routine!)
!      U - unchecked: yet to verify P or F


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


contains

! 1.Recombination Rates (pu01)

! Radiative recombination rate for HI, HeI, HeII and CIV (P)
  real function vernerFerland(a, t, t0, t1, b)
    real(double):: a, t, t0, t1, b
    
    ! based on Verner and Ferland (1996) ApJS 103 467
    vernerFerland = a / (sqrt(t/t0) * (1.d0+sqrt(t/t0))**(1.d0-b) * (1.d0+sqrt(t/t1))**(1.d0+b))
    
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

        e = hTable%freq(j) * hcgs* ergtoev
        call phfit2(1, 1, 1 , e , hxsec)

        dFreq = hTable%freq(j)-hTable%freq(j-1)
        jnu = ((hcgs*hTable%freq(j)**3)/(cSpeed**2)) * ((hcgs**2) /(twoPi*mElectron*Kerg*hTable%temp(i)))**(1.5d0) * &
             dble(hxsec/1.d10) *  exp(-hcgs*(hTable%freq(j)-nu0_h)/(kerg*hTable%temp(i)))
        hTable%Clyc(i,j) = hTable%Clyc(i,j-1) + jnu * dfreq
        
        e = heTable%freq(j) * hcgs * ergtoev
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
