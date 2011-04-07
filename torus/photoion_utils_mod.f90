!Photoionization utilities module
!Started 01/04/2011

!This module consolidates established functions and routines that are
!shared by both photoion_mod and photoionAMR_mod

! Structure Guide (quickSearchCode)
! 1. Recombination Rates (pu01)
! 2. Intersect Cube Routines(pu02)
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

!   2.Intersect Cube Routines (pu02)

! Removed by DMA as these routines are never called. 

!   End of pu02  


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
