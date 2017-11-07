#ifdef XRAY
module xray_mod
#ifdef MPI
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
use mpi_amr_mod
use mpi_global_mod
use grid_mod
use source_mod
use timing
use photoion_utils_mod
use mpi


implicit none


type AUGER
   real(double) :: yield
   real(double) :: energy     
   real(double) :: ionizthresh
end type AUGER

contains

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
subroutine doComptonXsecs(freq, nfreq)
  implicit none
  real(double), parameter :: sigma_thom = 6.65d-25
  real(double), parameter :: xraycutoff = 20.6d0*rydbergtoev
  real(double) :: freq
  real(double), allocatable :: CT_KN_Xsec(:), CT_heating(:), CT_cooling(:)
  real(double) :: x  !x is h*nu/(m*c^2)
  real(double) :: alpha, beta
  integer :: i, nfreq

  nfreq = 1000
  x= 0.d0
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

     !heating xsec
     CT_heating(i) = alpha*freq(i)*freq(i)*3.858d-25 
     !cooling xsec
     CT_cooling(i) = alpha*beta*freq(i)*3.858d-25 
     
     if(freq(i)*hConst <= xraycutoff) then
        CT_KN_Xsec(i) = sigma_thom
     else
        !calculate compton cross section using Klein-Nishina formula
        !equation 7.5 of Rybicki and Lightman (2004), page 197
        !This is for recoil ionization
        x = hConst*freq(i)/(mElectron*cspeed**2)
        CT_KN_Xsec(i) = sigma_thom * ((3.d0/4.d0)* ( &
             ((1.d0+x)/(x**3)) * ( &
             (2.d0*x*(1.d0+x)/(1.d0+2.d0*x)) - log(1.d0+2.d0*x)) + &
             (1.d0/2.d0*x)*log(1.d0+2.d0*x) - (1.d0-3.d0*x/(1.d0+2.d0*x)**2)))
     end if
  end do

end subroutine doComptonXsecs

!set up stuff for determing photon properties post compton-scatter
!uses Klein Nishina formula
subroutine inititateComptonScatterArray(KN_sigmaT, KN_sigmaArray, KN_PDF, freq)
  implicit none
  
  real(double) :: Eini_over_Enew
  integer, parameter :: ntheta=180
  real(double), allocatable :: KN_sigmaT(:), KN_sigmaArray(:,:), KN_PDF(:,:)
  real(double) :: freq(1000), thisTheta
  integer :: nfreq, i, j

  nfreq = 1000
  
  allocate(KN_sigmaT(nfreq))
  allocate(KN_sigmaArray(nfreq, ntheta))
  allocate(KN_PDF(nfreq, ntheta))

  KN_sigmaT = 0.d0
  KN_sigmaArray = 0.d0
  KN_PDF = 0.d0

  thisTheta = 0.d0

  do i = 1, nfreq
     do j = 1, ntheta
        thisTheta = dble(j)*pi/180.d0
        Eini_over_Enew = calcComptonEnergyChange(freq(i), thisTheta)

!the ratio of initial to final photon energies
        KN_PDF(i, j) = Eini_over_Enew

!scattering angle probability distribution
        KN_sigmaArray(i, j) = calcComptonScatterAngleProb(Eini_over_Enew, thisTheta)

!The total probability
        KN_sigmaT(i) = KN_SigmaT(i) + KN_sigmaArray(i,j)*2.d0*pi*sin(thisTheta)*pi/180.d0
     end do
  end do

  !Complete the PDF
  do i = 1, nfreq
     KN_sigmaArray(i,1) = KN_sigmaArray(i, 1)*2.d0*pi*sin(pi/180.d0)*pi/180.d0
     do j = 2, ntheta
        thisTheta = dble(j)*pi/180.d0
        KN_sigmaArray(i, j) = KN_sigmaArray(i, j-1) + KN_sigmaArray(i, j)* &
             2.d0*pi*sin(thistheta)*pi/180.d0
     end do
     KN_sigmaArray(i,:) = KN_sigmaArray(i,:)/KN_sigmaT(i)
     KN_sigmaArray(i, nTheta) = 1.d0
  end do


end subroutine inititateComptonScatterArray



!calculate the ratio of incoming to outgoing photon energy in a compton scattering event
!This is equation 7.2 of Rybicki and Lightman 2004
real(double) function calcComptonEnergyChange(frequency, angle)
  implicit none
  real(double) :: frequency
  real(double) :: angle
  real(double), parameter :: constant=37557.66267

!const = m_electron / c^2 * erg-->rydbergs
  
  calcComptonEnergyChange = 1.d0/(1.d0+ (frequency/constant)*(1.d0 - cos(angle)))

end function calcComptonEnergyChange



!uses Klein-Nishina formula
! dsigma/dOmega = 0.5 r_e^2 (P(E,theta)- P(E,theta)^2 Sin(theta)^2 + P(E,theta)^3)
!P = Eini/Enew
real(double) function calcComptonScatterAngleProb(Erat, theta)
  implicit none

  real(double) :: Erat, theta
  
  real(double), parameter :: electronRadius= 2.8179403267d-13 !not 10^10cm
  
  calcComptonScatterAngleProb = (0.5d0*electronRadius**2)* &
       (Erat - (Erat**2)*sin(theta**2)+Erat**3)

end function calcComptonScatterAngleProb


!!!get new photon direction and frequency following compton scattering
subroutine scatterCompton(newDirection, newFreq, oldfreq, KN_sigmaArray, KN_PDF, freq)
  implicit none
  type(vector) :: newDirection
  real(double) :: u, v, w, t, r1, r2, newtheta, newfreq, oldfreq
  integer :: i, index
  integer :: nfreq, ifreq, thisTheta
  integer, parameter :: ntheta = 180
  real(double) :: KN_sigmaArray(:,:), KN_PDF(:,:), freq(1000)

  nfreq = 1000
  newfreq = 0.d0
  call locate(freq, nFreq, oldFreq, iFreq)
  !get a random scattering angle from Klein-Nishina PDF  -------------
  call random_number(r2)
  do i = 1, 10000
     if(r2 == 0.d0 .or. r2 == 1.d0) then
        call random_number(r2)
     else
        exit
     end if
  end do
  if(i >= 10000) then
     call torus_abort("random number selection broken in xray_mod")
  end if

!
!find the frequency bin that this photon resides in
!

  do thisTheta = 1, ntheta
     if (r2 >= KN_sigmaArray(iFreq ,thisTheta)) then
        index = thisTheta
     else
        exit
     end if
  end do

  newtheta = index * pi/180.d0
!---------------------------------------------------------------------


!usual random direction stuff-----------------------------------------
  call random_number(r1)
  w = 2.d0*r1 - 1.d0
  t = sqrt(1.d0 - w*w)
  u = t*cos(newTheta)
  v = t*sin(newTheta)

  newDirection = VECTOR(u, v, w)
  
!---------------------------------------------------------------------


!get the new frequency------------------------------------------------

  newFreq = oldfreq * KN_PDF(iFreq, index)
  call locate(freq, nFreq, newFreq, index)
!---------------------------------------------------------------------!

end subroutine

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


!SIMPLIFIED XRAY treatment using the ionization parameter scheme of Owen+2010
subroutine simpleXRay(grid, thisSource)!
  use inputs_mod, only : xheat
  implicit none

  type(gridtype) :: grid
  type(sourcetype) :: thisSource
  integer :: i, ier
  
  if(myrankglobal /= 0 .and. xheat) then
     do i = 1, nhydrothreadsglobal
!        print *, myrankglobal, "starting ", i
        if(myrankglobal == i) then
!           print *, "thread ", myrankglobal, "doing raytracing"
           call calculateColumnToStar(grid%octreeRoot, grid, thisSource%position,thisSource%luminosity/1.d6)
!           print *, "thread ", myrankglobal, "done tracing"
           call shutdownserversXRAY()
        else
!           print *, "thread ", myrankglobal, "serving"
           call raytracingserverXRAY(grid)
!           print *, "thread ", myrankglobal, "done"
           
        end if
        call MPI_BARRIER(amrCommunicator, ier)
     end do            
  end if
  call MPI_BARRIER(MPI_COMM_WORLD, ier)
  if(xheat) then
     call calcIonParamTemperature(grid%octreeRoot, thisSource%luminosity/1.d6,  thisSource%position)
  endif

end subroutine simpleXRay

recursive subroutine calcIonParamTemperature(thisOctal, xrayLuminosity, starPos)

  integer :: j
  type(octal), pointer :: thisOctal, child!, soctal
  integer :: subcell!, ssubcell
  type(vector) :: starPos, rVec
  real(double) :: ionparam, xrayLuminosity, column, starDist
!  type(gridtype) :: grid

  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call calcIonParamTemperature(child, xrayluminosity, starpos)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
        if(.not. thisOctal%ghostcell(subcell)) then
!           if(.not. thisOctal%ionfrac(subcell,2 ) > 0.9d0) then
              
              !        if(.not. thisOctal%temperature(subcell) > 50.d0) then
              !           thisOCtal%temperature(subcell) = 10.d0 ! until I get the zeta(T) relation from James
              rVec = subcellCentre(thisOctal, subcell)
              starDist = modulus(rVec - starpos)*1.d10
              column = thisOctal%columnRho(subcell)*1.d10/mHydrogen
!                      print *, "columnrho ", column
!                      print *, "xraylum ", xrayluminosity
!                      print *, "stardist ", stardist
              !        ionParam = xrayLuminosity / ((column*starDist**2))
              if(column /= 0.d0 .and. starDist /= 0.d0) then
                 ionParam = xrayLuminosity / ((column*starDist))
              else
                 ionParam = 100.d0
              endif

              
              thisOctal%temperature(subcell) = real(calcOwen(ionParam, thisOctal%dust_t(subcell)))
!              print *, "ionparam, temperature ", ionParam, thisOctal%temperature(subcell)
 !          end if
        end if
  end if
end do

end subroutine calcIonParamTemperature

!subroutine printIonParam()
!  real(double) :: xi, xi_max, temperature
!  integer :: ier!!
!
!  xi = 1.d-7
!  xi_max = 1.d-3
!  open (1, file="ionParam", status="unknown", iostat=ier)
!  do while (xi < xi_max)
!     temperature =  real(calcOwen(xi))
!     write(1, *) xi, temperature
!     xi = xi + (xi/10.d0)
!  end do
!  
!end subroutine printIonParam

real(double) function calcOwen(thisxi, dustt)
  real(double), intent(in) :: thisxi
  real(double) :: tHot, tCold, xi, dustt
!FITTING CONSTANTS
  real(double), parameter :: ao=8.9362527959248299D-03
  real(double), parameter :: bo=-4.0392424905367275d0
  real(double), parameter :: co=1.2870891083912458D+01
  real(double), parameter :: do=4.4233310301789743D+01
  real(double), parameter :: eo=4.3469496951396964d0
  real(double), parameter :: fo=3.15d0
  real(double), parameter :: go=23.9d0


!  calcOwen = 10.d0**(a*i + b*i**(-2))/
!  print *, "xi ", thisxi
  xi = log10(thisxi)
!  print *, "log10(xi) ", xi
  tHot = 0.D0
  if(xi /= 0.d0) tHot = 10.d0**((ao*xi+bo*xi**(-2.))/(1+co*xi**(-1.)+do*xi**(-2.))+eo)
  
  tCold = 10.d0**((fo*xi)+go)
!  tCold = log10(tCold)
!  print *, tHot, tCold
  tCold = max(dustt, tCold)
  calcOwen = min(tHot, tCold)
!  calcOwen = max(calcOwen, 10.d0)

end function calcOwen

!get the column density from the cell to the star
!for use with the ionization parameter scheme used
!by James Owen
recursive subroutine calculateColumnToStar(thisOctal, grid, starpos, lum)
  use inputs_mod, only : sourceRadius, sourceTeff
  integer :: j
  type(octal), pointer :: thisOctal, child!, soctal
  integer :: subcell!, ssubcell
  type(vector) :: testposition, startposition, uhat, starpos
  real(double) :: tval, rho, disttostar, Lum
  type(gridtype) :: grid
  logical :: foundstar 

  do subcell = 1, thisoctal%maxchildren
     if (thisoctal%haschild(subcell)) then
        ! find the child
        do j = 1, thisoctal%nchildren, 1
           if (thisoctal%indexchild(j) == subcell) then
              child => thisoctal%child(j)
              call calculateColumnToStar(child, grid, starpos, lum)
              exit
           end if
        end do
     else
        if (.not.octalOnThread(thisOctal,subcell,myrankGlobal)) cycle
        thisOctal%dust_t(subcell)= 10.d0
        if(.not. thisOctal%ghostcell(subcell)) then
!        if(.not. thisOctal%ionfrac(subcell, 2) > 0.9d0) then
           !           sOctal=> thisOctal           
           startPosition = subcellCentre(thisOctal, subcell)
           testPosition = startPosition
!           call findSubcellLocal(testPosition, sOctal, ssubcell)
           thisOctal%columnRho(subcell) = 0.d0
           
!           uhat = startposition - starpos
           uhat = starpos - startposition
           call normalize(uHat)
!           print *, "star", starpos
!           print *, "testpos", testposition
!           print *, " uHat ", uHat
!           uHat = VECTOR(-Uhat%x, -uHat%y, -uHat%z)
!           uhat = -uhat
!           print *, "tracing ray with ", uhat
           foundstar = .false.
           disttostar = modulus(testposition-starpos)
           thisOctal%dust_t(subcell) = sourceTeff(1)*((sourceRadius(1)*rsol)/(2.d0*disttostar*1.d10))**0.5d0
           
           do while (.not. foundstar)     
 !             print *, "went one"
!              call findSubcellLocal(testPosition, sOctal, ssubcell)
!              if(.not octalonthread(thisoctal, myrankglobal)) then
              tval = 0.d0

  !            print *, "got values  ", tval
!!              call distanceToCellBoundary(grid, testPosition, uHat, tVal, soctal, ssubcell

              disttostar = modulus(testposition-starpos)
              if (disttostar < 2.d0*grid%halfsmallestsubcell) then
!                 call getRayTracingValuesXRAY(grid, testposition, uHat, rho,  tval, getTval=.false.) 
!                 thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (rho*disttostar)
!                 testPosition = testPosition + (disttostar*uhat)
                 foundstar = .true.
              else
                 call getRayTracingValuesXRAY(grid, testposition, uHat, rho,  tval, getTval=.true.) 
                 thisOctal%columnRho(subcell) = thisOctal%columnRho(subcell) + (rho*tval)
                 testPosition = testPosition + ((tVal+1.d-10*grid%halfsmallestsubcell)*uhat)

                 if(Lum / ((thisOctal%columnrho(subcell)*disttostar*1.d20/mhydrogen)) < 1.d-6) then
                    thisOctal%columnrho(subcell) = 1.d50
                    foundstar = .true.
                 endif

              endif              
            
           end do
        end if
!     endif
     end if
  end do
  
     
end subroutine calculateColumnToStar


#endif
end module xray_mod
#endif
