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
   real(double) :: energy     
   real(double) :: ionizthresh
end type AUGER

contains


!
!we are going to be as good as cloudy but not as comprehensive as mocassin (at least for now)
!we will be lacking in some atomic species
!


!Read in auger yield data
subroutine setUpAugerData(augerArray)
  implicit none
  integer, parameter :: nAtoms=5 !C, O, Mg, Si and Fe
  integer, parameter :: maxShells=5 !m23 shell for iron
  integer, parameter :: nLines=49
  character(len=16), parameter :: file="auger_yields.dat"  
  character(len=200):: dataDirectory, path
  character(len=200) :: message
  logical :: infile
  type(AUGER), intent(out) :: augerArray(nAtoms, maxShells, 10)
  integer :: ier, i
  integer :: elem, shell, ishell
  real(double) :: ionizthresh, energy, yield


  augerArray%yield = 0.d0
  augerArray%energy = 0.d0
  augerArray%ionizthresh = 1.d30
  infile = .true.

  call unixGetenv("TORUS_DATA", dataDirectory)
  dataDirectory = trim(dataDirectory)  

  path = trim(dataDirectory)//"/"//file
  open(1, file=path, status="old", position="rewind", iostat=ier)
  if(ier /= 0) then
     message = trim("trouble opening file auger_yields.dat in xray_mod")
     call torus_abort(message)
  end if

  do i = 1, nLines
     read(1,*) elem, shell, ionizthresh, energy, yield, ishell
     augerArray(elem, shell, ishell)%ionizthresh = ionizthresh   !eV
     augerArray(elem, shell, ishell)%energy = energy  !eV
     augerArray(elem, shell, ishell)%yield = yield   !average number of electrons
  end do

  call writeInfo("Auger array filled", TRIVIAL)

end subroutine setUpAugerData


!
!Compton scattering stuff
!

!this is the same as is used by cloudy and mocassin. 
subroutine getComptonThomsonXsec(freq, nfreq)
  implicit none
  real(double), parameter :: sigma_thom = 6.65d-25
  real(double), parameter :: xraycutoff = 20.6d0*rydbergtoev
  real(double) :: freq
  real(double), allocatable :: CT_KN_Xsec(:), CT_heating(:), CT_cooling(:)
  real(double) :: x  !x is h*nu/(m*c^2)
  real(double) :: alpha, beta
  integer :: i, nfreq

  nfreq = 1000

  if(.not. allocated(CT_KN_Xsec)) allocate(CT_KN_Xsec(nfreq))
  if(.not. allocated(CT_heating)) allocate(CT_heating(nfreq))
  if(.not. allocated(CT_cooling)) allocate(CT_cooling(nfreq))

  do i = 1, nfreq    

     !Compton exchange factors from
     !C. B. Tarter, W. H. Tucker, and E. E. Salpeter. 
     !The Interaction of X-Ray Sources with Optically Thin Environments. 
     !ApJ, 156:943, June 1969. doi: 10.1086/150026.
     alpha = 1.d0/(1.d0+freq(i)*(1.1792d-4+7.084d-10*freq(i)))
     beta = (1.d0 - alpha*freq(i)*(1.1792d-4+2.d0*7.084d-10*freq(i))/4.d0)

     !heating rate
     CT_heating(i) = alpha*freq(i)*freq(i)*3.858d-25 
     !cooling rate
     CT_cooling(i) = alpha*beta*freq(i)*3.858d-25 
     
     if(freq(i)*hConst <= xraycutoff) then
        CT_KN_Xsec(i) = sigma_thom
     else

        !calculate compton cross section using Klein-Nishina formula
        !equation 7.5 of Rybicki and Lightman (2004), page 197
        CT_KN_Xsec(i) = sigma_thom * ((3.d0/4.d0)* ( &
             ((1.d0+x)/(x**3)) * ( &
             (2.d0*x*(1.d0+x)/(1.d0+2.d0*x)) - log(1.d0+2.d0*x)) + &
             (1.d0/2.d0*x)*log(1.d0+2.d0*x) - (1.d0-3.d0*x/(1.d0+2.d0*x)**2)))
     end if
  end do


end subroutine getComptonThomsonXsec


!x-ray specific ionization balance (not yet working)
subroutine solveIonizationBalance_Xray(grid, thisOctal, subcell, temperature, epsOverdeltaT, augerArray)
  implicit none
  type(GRIDTYPE) :: grid
  type(OCTAL) :: thisOctal
  type(AUGER) :: augerArray(5, 5, 10)
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
!  real(double) :: source, sink, bet


  v = cellVolume(thisOctal, subcell)
  print *, "auger array", augerArray

  stop
!currently just old photoion stuff, under development
    k = 1
  allocate(newFrac(1:SIZE(thisOctal%ionFrac,2)))
  newFrac = thisOctal%ionFrac(subcell,:)
  do while(k <= grid%nIon)

     !For hydrogen and helium there is no inner shell photoionization

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
