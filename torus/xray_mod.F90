#ifdef XRAY
module xray_mod
!
!Started by Thaw - 22/07/2013
!
!
!Module for housing stuff to do with x-ray enabled torus
!

use constants_mod
use messages_mod
use parallel_mod
use amr_mod
use vtk_mod
use mpi_amr_mod
use mpi_global_mod
use grid_mod
use source_mod
use timing
use photoion_utils_mod

implicit none

type AUGER
   real(double) :: yield
   real(double) energy     
end type AUGER

contains


!Read in auger yield data
subroutine setUpAugerData()
  implicit none
  integer, parameter :: nAtoms=5 !C, O, Mg, Si and Fe
  integer, parameter :: maxShells=5 !m23 shell for iron
  character(len=200):: dataDirectory
  type(AUGER) :: augerArray(nAtoms,maxShells)
    
  augerArray%yield = 0.d0
  augerArray%energy = 0.d0

  call unixGetenv("TORUS_DATA", dataDirectory)
  open(1, file=trim(dataDirectory)//'auger_yields.dat', status="old", position="rewind")

  
  

end subroutine setUpAugerData


!Compton scattering stuff



!x-ray specific ionization balance (not yet working)
subroutine solveIonizationBalance_Xray(grid, thisOctal, subcell, temperature, epsOverdeltaT)
  implicit none
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  real(double) :: epsOverDeltaT, V
  real :: temperature
  integer :: subcell
  integer :: i, k
  integer :: nIonizationStages
  integer :: iStart, iEnd, iIon
  real(double) :: chargeExchangeIonization, chargeExchangeRecombination
  real(double), allocatable :: newFrac(:)
  real(double), allocatable :: xplus1overx(:)
  real(double), parameter :: underCorrection = 0.6d0

  v = cellVolume(thisOctal, subcell)

!currently just old photoion stuff, under development
    k = 1
  allocate(newFrac(1:SIZE(thisOctal%ionFrac,2)))
  newFrac = thisOctal%ionFrac(subcell,:)
  do while(k <= grid%nIon)

     iStart = k
     iEnd = k+1
     do while(grid%ion(istart)%z == grid%ion(iEnd)%z)
        iEnd = iEnd + 1
     enddo
     iEnd = iEnd - 1
     nIonizationStages = iEnd - iStart + 1
          
     allocate(xplus1overx(1:nIonizationStages-1))
     do i = 1, nIonizationStages-1
        iIon = iStart+i-1
        call getChargeExchangeRecomb(grid%ion(iion+1), temperature, &
             thisOctal%nh(subcell)*grid%ion(1)%abundance*thisOctal%ionFrac(subcell,1),  &
             chargeExchangeRecombination)
        
        call getChargeExchangeIon(grid%ion(iion), temperature, &
             thisOctal%nh(subcell)*grid%ion(2)%abundance*thisOctal%ionFrac(subcell,2),  &
             chargeExchangeIonization)
           
        xplus1overx(i) = ((epsOverDeltaT / (v * 1.d30))*thisOctal%photoIonCoeff(subcell, iIon) + chargeExchangeIonization) / &
             max(1.d-50,(recombRate(grid%ion(iIon),temperature) * thisOctal%ne(subcell) + chargeExchangeRecombination))
     enddo


     newFrac(iStart:iEnd) = 1.d0
     do i = 1, nIonizationStages - 1
        iIon = iStart+i-1
        newFrac(iIon+1) = newFrac(iIon) * xplus1overx(i)
     enddo
     
     if (SUM(newFrac(iStart:iEnd)) /= 0.d0) then
        newFrac(iStart:iEnd) = &
             max(1.d-50,newFrac(iStart:iEnd))/SUM(newFrac(iStart:iEnd))
     else
        newFrac(iStart:iEnd) = 1.d-50
     endif

        deallocate(xplus1overx)

        k = iEnd + 1
     end do

     thisOctal%ionFrac(subcell,:) = thisOctal%ionFrac(subcell,:) + underCorrection * (newFrac - thisOctal%ionFrac(subcell,:))
     deallocate(newFrac)
     thisOctal%ne(subcell) = returnNe(thisOctal, subcell, grid%ion, grid%nion)

end subroutine solveIonizationBalance_Xray

end module xray_mod
#endif
