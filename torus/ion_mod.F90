#ifdef PHOTOION
module ion_mod
!
!

! written by tjh
! updated 28/02/09 by tjhaworth and ajp (inc. carbon)


  use kind_mod
  use constants_mod
  use utils_mod
  use messages_mod

  implicit none

  type LEVELTYPE
     character(len=10) :: term
     real :: g
     real(double) :: energy
  end type LEVELTYPE

  type TRANSITIONTYPE
     integer :: i  ! lower level
     integer :: j  ! upper level
     real(double) :: lambda ! angstrom
     real(double) :: energy
     integer :: nGamma
     real :: gamma(20)  ! collisional gamma at T 
     real :: t(20) ! temperatures at which gammas are defined
     real :: A ! Einstein A coefficient
  end type TRANSITIONTYPE

  ! The definition of the ion type

  type IONTYPE
     integer :: z         ! nuclear charge
     integer :: n         ! number of electrons
     integer :: outerShell    ! outer shell
     real :: ipot         ! ionization potential
     real :: abundance    ! of element relative to H
     real(double) :: nuthresh ! freq equivalent of ion pot
     real(double) :: nucleons     ! number of nucleons, use to calculate mass
     integer :: nFreq  ! number of frequency bins in xsection array
     real(double), pointer :: freq(:)=>null(), xSec(:)=>null() ! frequency and ionization xsection arrays
     character(len=8) :: species  ! name of ion
     integer :: nTransitions ! number of bb transitions for this species
     integer :: nLevels ! number of energy levels
     type(LEVELTYPE) :: level(10) 
     type(TRANSITIONTYPE) :: transition(40) 
  end type IONTYPE

! Abundances 
  real, private, save :: h_abund, he_abund, c_abund, n_abund, o_abund, ne_abund, s_abund 

  type(IONTYPE), save :: globalIonArray(50)
  integer :: nGlobalIon

contains

  subroutine setAbundances(h_in, he_in, c_in, n_in, o_in, ne_in, s_in)
    real, intent(in) :: h_in, he_in, c_in, n_in, o_in, ne_in, s_in

    h_abund  = h_in
    he_abund = he_in
    c_abund  = c_in
    n_abund  = n_in 
    o_abund  = o_in
    ne_abund = ne_in
    s_abund  = s_in

  end subroutine setAbundances

  subroutine createIon(thisIon, z, n, nucleons, ipot)
    type(IONTYPE) :: thisIon
    integer :: z, n, i, nucleons
    real :: iPot
    character(len=2) :: element
    character(len=4) :: roman
    character(len=30) :: message

    element = " "; roman = " "
    thisIon%z = z
    thisIon%n = n
    thisIon%nucleons = nucleons
    i = z - n + 1
    call returnElement(z, element)
    call createRoman(i, roman)
    thisIon%species = trim(element)//" "//trim(roman)
    thisIon%abundance = returnAbundance(z)
    write(message, '(a10,f11.8)') thisIon%species, thisIon%abundance
    if (writeoutput) call writeInfo(message, TRIVIAL)
    thisIon%iPot = iPot
    thisIon%outerShell = returnOuterShell(n) 
    thisIon%nuThresh = (thisIon%iPot / dble(ergtoev)) / dble(hCgs)
    thisIon%nLevels = 0
    thisIon%nTransitions = 0
    call addLevels(thisIon)
    call addTransitions(thisIon)

  end subroutine createIon

  subroutine addxSectionArray(thisIon, nfreq, freq)
    use phfit_mod, only : phfit2
    type(IONTYPE) :: thisIon
    integer :: nFreq
    real(double) :: freq(:)
    real :: e, xSec
    integer :: i

    thision%nfreq = nFreq
    if (associated(thisIon%Freq)) deallocate(thisIon%freq)
    if (associated(thisIon%xsec)) deallocate(thisIon%xsec)
    allocate(thisIon%freq(1:nFreq), thisIon%xSec(1:nFreq))
    thision%freq(1:nFreq) = freq(1:nFreq)

    thisIon%Xsec = 0.d0
    do  i = 1, nFreq
       e = real(hCgs * freq(i) * ergtoEv)
!       print *, "e", e
 !      print *, "thisIon%iPot", thisIon%iPot
  !     print *, "i", i
   !    print *, " "
       if (e > thisIon%iPot) then
          call phfit2(thisIon%z, thisIon%n, thisIon%outerShell , e , xsec)
          thisIon%xSec(i) = dble(xSec)
       endif
    enddo
  end subroutine addxSectionArray

  function returnXsec(thisIon, freq, iFreq) result(xSec)
    type(IONTYPE) :: thisIon
    real(double) :: freq, xSec
    integer, optional :: iFreq
    integer :: i

       if (PRESENT(iFreq)) then
!          print *, "iFreq", iFreq
!          print *, "thisIon%xSec(iFreq)", thisIon%xSec(iFreq)
          xSec = thisIon%xSec(iFreq)
       else
          call locate(thisIon%freq, thisIon%nFreq, freq, i)
          xSec = thisIon%xSec(i)
       endif
  end function returnXsec

  subroutine addIons(ionArray, nIon, usemetals, hOnly)
    logical, intent(in) :: usemetals, hOnly
    integer :: nIon
    type(IONTYPE) :: ionArray(:)


    nIon = 1
    call createIon(ionArray(nIon), 1, 1, 1, 1.360e1) ! H I 1

    nIon = nIon + 1
    call createIon(ionArray(nIon), 1, 0, 1, 1.e-10) ! H II 2

    if (.not. hOnly) then

       nIon = nIon + 1
       call createIon(ionArray(nIon), 2, 2, 4, 2.459e1) ! He I 3

       nIon = nIon + 1
       call createIon(ionArray(nIon), 2, 1, 4, 5.442e1) ! He II 4

       nIon = nIon + 1
       call createIon(ionArray(nIon), 2, 0, 4, 1.e-10) ! He III (alpha particle!) 5

    end if

    if (usemetals) then

       nIon = nIon + 1
       call createIon(ionArray(nIon), 6, 6, 12, 1.126e1) ! C I

       nIon = nIon + 1
       call createIon(ionArray(nIon), 6, 5, 12, 2.438e1) ! C II

       nIon = nIon + 1
       call createIon(ionArray(nIon), 6, 4, 12, 4.789e1) ! C III

       nIon = nIon + 1
       call createIon(ionArray(nIon), 6, 3, 12, 6.449e1) ! C IV

       !    nIon = nIon + 1
       !    call createIon(ionArray(nIon), 6, 2, 12, 3.921e2) ! C V

       nIon = nIon + 1
       call createIon(ionArray(nIon), 7, 7, 14, 1.453e1) ! N I

       nIon = nIon + 1
       call createIon(ionArray(nIon), 7, 6, 14, 2.960e1) ! N II

       nIon = nIon + 1
       call createIon(ionArray(nIon), 7, 5, 14, 4.745e1) ! N III

!THAW - removed 03/03/2011 since no corresponding level data present.
!±       nIon = nIon + 1
!       call createIon(ionArray(nIon), 7, 4, 14, 7.747e1) ! N IV

       !    nIon = nIon + 1
       !    call createIon(ionArray(nIon), 7, 3, 9.789e1) ! N V

       nIon = nIon + 1
       call createIon(ionArray(nIon), 8, 8, 16, 1.362e1) ! O I

       nIon = nIon + 1
       call createIon(ionArray(nIon), 8, 7, 16, 3.512e1) ! O II

       nIon = nIon + 1
       call createIon(ionArray(nIon), 8, 6, 16, 5.494e1) ! O III

!THAW
!       nIon = nIon + 1
!       call createIon(ionArray(nIon), 8, 5, 16, 7.741e1) ! O IV

       !    nIon = nIon + 1
       !    call createIon(ionArray(nIon), 8, 4, 16, 1.139e2) ! O V

       nIon = nIon + 1
       call createIon(ionArray(nIon), 10, 10, 20, 2.156e1) ! Ne I

       nIon = nIon + 1
       call createIon(ionArray(nIon), 10, 9, 20, 4.096e1) ! Ne II

       nIon = nIon + 1
       call createIon(ionArray(nIon), 10, 8, 20, 6.342e1) ! Ne III

!       nIon = nIon + 1
!       call createIon(ionArray(nIon), 10, 7, 20, 9.712e1) ! Ne IV

       !    nIon = nIon + 1
       !    call createIon(ionArray(nIon), 10, 6, 20, 1.262e2) ! Ne V

!THAW
       nIon = nIon + 1
       call createIon(ionArray(nIon), 16, 16, 32, 1.036e1) ! S I

       nIon = nIon + 1
       call createIon(ionArray(nIon), 16, 15, 32, 2.334e1) ! S II

       nIon = nIon + 1
       call createIon(ionArray(nIon), 16, 14, 32, 3.479e1) ! S III

       nIon = nIon + 1
       call createIon(ionArray(nIon), 16, 13, 32, 4.722e1) ! S IV
    endif
    if (writeoutput) &
         write(*,*) "Added ",nion," species to photoionization calculation"

  end subroutine addIons
    
function returnOuterShell(n) result (ishell)
  integer :: n, ishell
! Shell numbers:
! 1 - 1s, 2 - 2s, 3 - 2p, 4 - 3s, 5 - 3p, 6 - 3d, 7 - 4s. 

  if (n <= 2) then
     ishell = 1
  else if (n > 2 .and. n <= 4) then
     ishell = 2
  else if (n > 4 .and. n <= 10) then
     ishell = 3
  else if (n > 10 .and. n <= 12) then
     ishell = 4
  else if (n > 12 .and. n <= 18) then
     ishell = 5
  else 
     write(*,*) "No configuration for ",n," electrons"
     ishell = 0 
  endif
end function returnOuterShell
  
function returnIonNumber(species, ionArray, nIon) result (iIon)
  type(IONTYPE) :: ionArray(:)
  integer :: nIon, i, iIon 
  character(len=*) :: species

  iIon = 0
  loop: do i = 1, nIon
     if (trim(ionArray(i)%species) == trim(species)) then
        iIon = i
        exit loop
     endif
  enddo loop

end function returnIonNumber


subroutine addLevels(thisIon)

  type(IONTYPE) :: thisIon

  select case(thisIon%species)
!Carbon (I, II, III, IV) energies taken from NIST atomic spectra database levels data 
  case("C I")
     call addLevel(thisIon,"3P_0", 0.d0, 0.)
     call addLevel(thisIon,"3P_1", 0.002033d0, 1.)
     call addLevel(thisIon,"3P_2", 0.005381d0, 2.)
     call addLevel(thisIon,"1D_2", 1.263725d0, 2.)
     call addLevel(thisIon,"1S_0", 2.684011d0, 0.)
     call addLevel(thisIon,"5S^0_2", 4.182631d0, 2.)

  case("C II")
     call addLevel(thisIon,"2P^0_1/2", 0.d0, 0.5)   
     call addLevel(thisIon,"2P^0_3/2", 0.007863d0, 1.5)
     call addLevel(thisIon,"4P_1/2", 5.33173d0, 0.5)
     call addLevel(thisIon,"4P_3/2", 5.33446d0, 1.5)
     call addLevel(thisIon,"4P_5/2", 5.33797d0, 2.5)
     !THAW
     call addLevel(thisIon,"2D_3/2", 9.290460d0 , 1.5)
     call addLevel(thisIon,"2D_5/2", 9.290148d0 , 2.5)
     call addLevel(thisIon,"2S_1/2", 11.96370d0, 0.5)
     call addLevel(thisIon,"2P_1/2", 13.715648d0, 0.5)
     call addLevel(thisIon,"2P_3/2", 13.720780d0, 1.5)

  case("C III")
     call addLevel(thisIon, "1S_0", 0.d0, 0.)
     call addLevel(thisIon, "3P^0_0", 6.492688d0, 0.)
     call addLevel(thisIon, "3P^0_1", 6.495625d0, 1.)
     call addLevel(thisIon, "3P^0_2", 6.502612d0, 2.)
     call addLevel(thisIon, "1P^0_1", 12.690035d0,1.)

  case("C IV")
     call addLevel(thisIon, "2S_1/2", 0.d0, 0.5)
     call addLevel(thisIon, "2P^0_1/2", 7.99500d0, 0.5)
     call addLevel(thisIon, "2P^0_3/2", 8.00835d0, 1.5)

  case("N I")
     call addLevel(thisIon, "4S^0_3/2", 0.d0, 1.5)
     call addLevel(thisIon, "2D^0_5/2", 2.383530d0, 2.5)
     call addLevel(thisIon, "2D^0_3/2", 2.384610d0, 1.5)
     call addLevel(thisIon, "2P^0_1/2", 3.575571d0, 0.5)
     call addLevel(thisIon, "2P^0_3/2", 3.575619d0, 1.5)

  case("N II")
     call addLevel(thisIon, "3P_0", 0.d0, 0.)
     call addLevel(thisIon, "3P_1", 0.00604d0, 1.)
     call addLevel(thisIon, "3P_2", 0.01622d0, 2.)
     call addLevel(thisIon, "1D_2", 1.89897d0, 2.)
     call addLevel(thisIon, "1S_0", 4.05290d0, 2.)
     call addLevel(thisIon, "5S^0_2", 5.80055d0, 2.)

  case("N III")
     call addLevel(thisIon, "2P^0_1/2", 0.d0, 0.5)
     call addLevel(thisIon, "2P^0_3/2", 0.02162d0, 1.5)
     call addLevel(thisIon, "4P_1/2", 7.09030d0, 0.5)
     call addLevel(thisIon, "4P_3/2", 7.09770d0, 1.5)
     call addLevel(thisIon, "4P_5/2", 7.10776d0, 2.5)

!  case("N IV")
!     call addLevel(thisIon, "1S_0", 0.d0, 0.)
!     call addLevel(thisIon, "3P^0_0", ???, 0.)
!     call addLevel(thisIon, "3P^0_1", ???, 1.)
!     call addLevel(thisIon, "3P^0_2", ???, 2.)
!     call addLevel(thisIon, "1P^0_1", ???, 1.)

  case("O I")
     call addLevel(thisIon, "3P_2", 0.d0, 2.)
     call addLevel(thisIon, "3P_1", 0.019622d0, 1.)
     call addLevel(thisIon, "3P_0", 0.028141d0, 0.)
     call addLevel(thisIon, "1D_2", 1.967364d0, 2.)
     call addLevel(thisIon, "1S_0", 4.189747d0, 0.)

  case("O II")
     call addLevel(thisIon, "4S^0_3/2", 0.d0, 1.5)
     call addLevel(thisIon, "2D^0_5/2", 3.324084d0, 2.5)
     call addLevel(thisIon, "2D^0_3/2", 3.326567d0, 1.5)
     call addLevel(thisIon, "2P^0_3/2", 5.017393d0, 1.5)
     call addLevel(thisIon, "2P^0_1/2", 5.017640d0, 0.5)


  case("O III")
     call addLevel(thisIon, "3P_0", 0.d0, 0.)
     call addLevel(thisIon, "3P_1", 0.014033d0, 1.)
     call addLevel(thisIon, "3P_2", 0.0379607d0, 2.)
     call addLevel(thisIon, "1D_2", 2.51365d0, 2.)
     call addLevel(thisIon, "1S_0", 5.354349d0, 0.)
     call addLevel(thisIon, "5S^0_2", 7.479320d0, 2.)
        
  case("Ne II")
     call addLevel(thisIon, "2P^0_3/2", 0.d0, 1.5)
     call addLevel(thisIon, "2P^0_1/2", 0.096750d0, 1.5)
!        call addLevel(thisIon, "2S_1/2", 26.91048d0, 1.5)
!        call addLevel(thisIon, "4P_5/2", 27.16876d0, 2.5)
!        call addLevel(thisIon, "4P_3/2", 27.2394d0, 1.5)
!        call addLevel(thisIon, "4P_1/2", 27.2700d0, 0.5)

  case("Ne III")
     call addLevel(thisIon, "3P_2", 0.d0, 2.)
     call addLevel(thisIon, "3P_1", 0.07971d0, 1.)
     call addLevel(thisIon, "3P_0", 0.11413d0, 0.)
     call addLevel(thisIon, "1D_2", 3.20385d0, 2.)
     call addLevel(thisIon, "1S_0", 6.91220d0, 0.)

  case("S II")
     call addLevel(thisIon, "4S^0_3/2", 0.d0, 1.5)
     call addLevel(thisIon, "2D^0_3/2", 1.841530d0, 1.5)
     call addLevel(thisIon, "2D^0_5/2", 1.845471d0, 2.5)
     call addLevel(thisIon, "2P^0_1/2", 3.040691d0, 0.5)
     call addLevel(thisIon, "2P^0_3/2", 3.046482d0, 1.5)
!        call addLevel(thisIon, "4P_5/2", 9.83774d0, 2.5)
!        call addLevel(thisIon, "4P_3/2", 9.880586d0, 1.5)
!        call addLevel(thisIon, "4P_1/2", 9.914099d0, 0.5)
!
  case("S III")
     call addLevel(thisIon, "3P_0", 0.d0, 0.)
     call addLevel(thisIon, "3P_1", 0.037033d0, 1.)
     call addLevel(thisIon, "3P_2", 0.10329d0, 2.)
     call addLevel(thisIon, "1D_2", 1.403841d0, 2.)
     call addLevel(thisIon, "1S_0", 3.36753d0, 0.)
     !   call addLevel(thisIon, "5S^0_2", 7.274391d0, 2.)

  case("S IV")
     call addLevel(thisIon, "2P^0_1/2", 0.d0, 0.5)
     call addLevel(thisIon, "2P^0_3/2", 0.117962d0, 1.5)

  case DEFAULT
  end select
end subroutine addLevels

subroutine addLevel(thisIon, term, energy, J)
  type(IONTYPE) :: thisIon
  character(len=*) :: term
  real(double) :: energy
  real :: j
  integer :: i

  i = thisIon%nLevels + 1
  thisIon%nLevels = i
  thisIon%level(i)%term = term
  thisIon%level(i)%energy = energy
  thisIon%level(i)%g = real(2.*J+1.)
end subroutine addLevel

subroutine addTransitions(thisIon)
  type(IONTYPE) :: thisIon
  real, allocatable :: t(:), gamma(:)

  select case(thisIon%species)

!THaw - freshly modding carbon to agree with MOCASSIN
  case("C I")
     allocate(t(6), gamma(6))
     t = (/ 500, 1000, 5000, 10000, 15000, 20000/)
     gamma = (/ 0.009,  0.022,  0.243,  0.371,  .4,  .4/) 
     call addTransition2(thisIon, "3P_1", "3P_0", 6.094E6, 7.93e-8, t, gamma, 6)
     gamma = (/ 0.006,  0.017,  0.182,  0.246,  .28,  .28/) 
     call addTransition2(thisIon, "3P_2", "3P_0", 2304147., 1.71e-14, t, gamma, 6)
     gamma = (/ 0.00694,  0.01389,  0.0670,  0.12667,  0.17778,  0.21778/) 
     call addTransition2(thisIon, "1D_2", "3P_0", 9811.03, 7.77e-8, t, gamma, 6)
     !1-5
     gamma = (/ 0.00191, 0.00377, 0.01656, 0.02800, 0.03556, 0.04056 /)
     call addTransition2(thisIon, "1S_0", "3P_0", 4619.43, 1.e-20, t, gamma, 6)
     !1-6
     gamma = (/ 0.01667, 0.02356, 0.05278, 0.07456, 0.09133, 0.10556 /)
     call addTransition2(thisIon, "5S^0_2", "3P_0", 2964.31, 1.e-20, t, gamma, 6)
     gamma = (/ 0.025,   0.065,   0.714,   1.020,   1.1,   1.1/) 
     call addTransition2(thisIon, "3P_2", "3P_1", 3704140., 2.65e-7, t, gamma, 6)
     gamma = (/ 0.02083,  0.0417,  0.201,  0.380,  0.533,  0.6533/) 
     call addTransition2(thisIon, "1D_2", "3P_1", 9824.12, 8.21e-5, t, gamma, 6)
     gamma = (/ 0.00573,  0.01130,  0.04967,  0.0840,  0.1067,  0.1217/) 
     call addTransition2(thisIon, "1S_0", "3P_1", 4621.57, 2.71e-3, t, gamma, 6)
     gamma = (/ 0.0500,  0.0707,  0.1583,  0.2237,  0.2740,  0.3167/) 
     call addTransition2(thisIon, "5S^0_2", "3P_1", 2965.70, 6.94, t, gamma, 6)
     gamma = (/ 0.03472,  0.06944,  0.3354,  0.6333,  0.8889,  1.0889/) 
     call addTransition2(thisIon, "1D_2", "3P_2", 9850.28, 2.44e-4, t, gamma, 6)
     gamma = (/ 0.00956,  0.01883,  0.08278,  0.1400,  0.17778,  0.20278/) 
     call addTransition2(thisIon, "1S_0", "3P_2", 4628.64, 2.00e-5, t, gamma, 6)
     gamma = (/ 0.0833,  0.1178,  0.2639,  0.3728,  0.4567,  0.5278/) 
     call addTransition2(thisIon, "5S^0_2", "3P_2", 2968.08, 1.56e1, t, gamma, 6)
     gamma = (/ 0.0620,  0.0877,  0.196,  0.277,  0.340,  0.392/) 
     call addTransition2(thisIon, "1S_0", "1D_2", 8727.18, 5.28e-1, t, gamma, 6)

     deallocate(t, gamma)
     
  case("C II") 

     allocate(t(12), gamma(12))
     t = (/ 2000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000 /)
     gamma = (/ 1.64, 1.82, 2.12, 2.30, 2.31, 2.27, 2.21, 2.16, 2.11, 2.05, 2.01, 1.96/)
     call addTransition2(thisIon, "2P^0_3/2", "2P^0_1/2", 1.5774e5, 2.29e-6,t, gamma, 12)
     gamma = (/ 0.291, 0.285, 0.280, 0.277, 0.272, 0.265, 0.258, 0.251, 0.245, 0.240, 0.235, 0.230/) 
     call addTransition2(thisIon, "4P_1/2", "2P^0_1/2", 2325., 5.53e1,t, gamma, 12)
     gamma = (/ 0.428, 0.420, 0.414, 0.411, 0.402, 0.393, 0.382, 0.373, 0.364, 0.356, 0.348, 0.342/) 
     call addTransition2(thisIon, "4P_3/2", "2P^0_1/2", 2324., 1.71,t, gamma, 12)
     gamma = (/ 0.256, 0.253, 0.254, 0.256, 0.253, 0.246, 0.240, 0.234, 0.228, 0.223, 0.219, 0.215/) 
     call addTransition2(thisIon, "4P_5/2", "2P^0_1/2", 2323., 1.e-20,t, gamma, 12)
     gamma = (/ 0.196, 0.194, 0.194, 0.195, 0.192, 0.187, 0.182, 0.178, 0.173, 0.170, 0.167, 0.164/) 
     call addTransition2(thisIon, "4P_1/2", "2P^0_3/2", 2329., 6.55e1,t, gamma, 12)
     gamma = (/ 0.547, 0.538, 0.535, 0.534, 0.524, 0.511, 0.498, 0.485, 0.473, 0.463, 0.454, 0.446/) 
     call addTransition2(thisIon, "4P_3/2", "2P^0_3/2", 2328., 5.24,t, gamma, 12)
     gamma = (/ 0.740, 0.739, 0.796, 0.881, 0.930, 0.968, 1.00, 1.03, 1.05, 1.07, 1.09, 1.11/) 
     call addTransition2(thisIon, "4P_3/2", "4P_1/2", 4.55e6, 2.39e-7,t, gamma, 12)
     gamma = (/ 1.21, 1.18, 1.17, 1.16, 1.14, 1.11, 1.08, 1.05, 1.03, 1.00, 0.982, 0.963/) 
     call addTransition2(thisIon, "4P_5/2", "2P^0_3/2", 2326., 4.32e1,t, gamma, 12)
     gamma = (/ 0.498, 0.517, 0.561, 0.621, 0.658, 0.688, 0.713, 0.733, 0.750, 0.764, 0.776, 0.787/) 
     call addTransition2(thisIon, "4P_5/2", "4P_1/2", 1.99e6, 1.e-20,t, gamma, 12)
     gamma = (/ 1.42, 1.44, 1.56, 1.73, 1.82, 1.90, 1.97, 2.02, 2.07, 2.11, 2.14, 2.17/)
     call addTransition2(thisIon, "4P_5/2", "4P_3/2", 3.53e6, 3.67e-7,t, gamma, 12)

     !Thaw
     !1-6
     gamma = (/ 1.10, 1.14, 1.20, 1.27, 1.31, 1.34, 1.35, 1.36, 1.36, 1.36, 1.36, 1.36/)
     call addTransition2(thisIon, "2D_3/2", "2P^0_1/2", 1334.55, 2.42E+8,t, gamma, 12)

     !1-7
     gamma = (/ 0.658, 0.633, 0.624, 0.613, 0.597, 0.581, 0.567, 0.553, 0.542, 0.532, 0.523, 0.515 /)
     call addTransition2(thisIon, "2D_5/2", "2P^0_1/2", 1334.59, 1.e-20,t, gamma, 12)

     !1-8
     gamma = (/ 0.802, 0.760, 0.736, 0.721, 0.721, 0.727, 0.734, 0.742, 0.749, 0.756, 0.762, 0.767/)
     call addTransition2(thisIon, "2S_1/2", "2P^0_1/2", 1036.35 , 7.74E+8,t, gamma, 12)

     !1-9
     gamma = (/ 1.05, 1.12, 1.21, 1.34, 1.43, 1.50, 1.56, 1.61, 1.66, 1.70, 1.74, 1.77 /)
     call addTransition2(thisIon, "2P_1/2", "2P^0_1/2", 903.98 , 2.74E+9,t, gamma, 12)

     !1-10
     gamma = (/ 0.603, 0.636, 0.674, 0.725, 0.759, 0.785, 0.806, 0.822, 0.837, 0.850, 0.857, 0.863 /)
     call addTransition2(thisIon, "2P_3/2", "2P^0_1/2", 903.64 , 6.86E+8,t, gamma, 12)

     !2-6
     gamma = (/ 1.01, 0.995, 0.998, 1.00, 0.996, 0.985, 0.973, 0.961, 0.950, 0.940, 0.930, 0.922 /)
     call addTransition2(thisIon, "2D_3/2", "2P^0_3/2", 1335.68, 4.78E+7,t, gamma, 12)

     !2-7
     gamma = (/ 1.77, 1.78, 1.81, 1.84, 1.85, 1.84, 1.82, 1.81, 1.79, 1.78, 1.77, 1.76 /)
     call addTransition2(thisIon, "2D_5/2", "2P^0_3/2", 1335.73, 2.88E+8,t, gamma, 12)

     !2-8
     gamma = (/ 1.61, 1.53, 1.48, 1.46, 1.47, 1.50, 1.52, 1.55, 1.58, 1.60, 1.62, 1.64 /)
     call addTransition2(thisIon, "2S_1/2", "2P^0_3/2", 1037.03 , 1.53E+9,t, gamma, 12)

     !2-9
     gamma = (/ 0.628, 0.667, 0.720, 0.794, 0.850, 0.895, 0.935, 0.969, 0.999, 1.03, 1.05, 1.07 /)
     call addTransition2(thisIon, "2P_1/2", "2P^0_3/2", 904.49 , 1.38E+9,t, gamma, 12)
     
     !2-10
     gamma = (/ 2.67, 2.83, 3.03, 3.30, 3.49, 3.64, 3.78, 3.89, 3.99, 4.07, 4.14, 4.21/)
     call addTransition2(thisIon, "2P_3/2", "2P^0_3/2", 904.16 , 3.43E+9,t, gamma, 12)

     !3-6
     gamma = (/ 0.560, 0.559, 0.562, 0.569, 0.573, 0.573, 0.577, 0.577, 0.576, 0.573, 0.571, 0.568/)
     call addTransition2(thisIon, "2D_3/2", "4P_1/2", 3131.97 , 1.e-20,t, gamma, 12)

     !3-7
     gamma = (/ 0.336, 0.337, 0.339, 0.343, 0.346, 0.349, 0.351, 0.352, 0.352, 0.352, 0.351, 0.350/)
     call addTransition2(thisIon, "2D_5/2", "4P_1/2", 3132.22 , 1.e-20,t, gamma, 12)

     !3-8
     gamma = (/ 0.092, 0.085, 0.0858, 0.087, 0.0878, 0.0888, 0.0897, 0.0906, 0.0914, 0.0921, 0.0926, 0.0932/)
     call addTransition2(thisIon, "2S_1/2", "4P_1/2", 1869.52 , 1.e-20,t, gamma, 12)

     !3-9
     gamma = (/ 0.057, 0.069, 0.100, 0.133, 0.146, 0.155, 0.159, 0.162, 0.167, 0.170, 0.174, 0.178/)
     call addTransition2(thisIon, "2P_1/2", "4P_1/2", 1478.86 , 1.e-20,t, gamma, 12)

     !3-10
     gamma = (/ 0.035, 0.048, 0.085, 0.122, 0.136, 0.143, 0.148, 0.152, 0.156, 0.160, 0.163, 0.166/)
     call addTransition2(thisIon, "2P_3/2", "4P_1/2", 1477.95 , 1.e-20,t, gamma, 12)

     !4-6
     gamma = (/ 0.877, 0.877, 0.882, 0.892, 0.900, 0.905, 0.907, 0.907, 0.906, 0.903, 0.900, 0.895/)
     call addTransition2(thisIon, "2D_3/2", "4P_3/2", 3134.13 , 1.e-20,t, gamma, 12)

     !4-7
     gamma = (/ 0.913, 0.915, 0.919, 0.930, 0.938, 0.944, 0.947, 0.947, 0.945, 0.943, 0.939, 0.934/)
     call addTransition2(thisIon, "2D_5/2", "4P_3/2", 3134.38 , 1.e-20,t, gamma, 12)

     !4-8
     gamma = (/ 0.184, 0.171, 0.172, 0.174, 0.176, 0.177, 0.179, 0.181, 0.183, 0.184, 0.185, 0.186/)
     call addTransition2(thisIon, "2S_1/2", "4P_3/2", 1870.30 , 1.e-20,t, gamma, 12)

     !4-9
     gamma = (/ 0.083, 0.102, 0.153, 0.209, 0.230, 0.242, 0.249, 0.256, 0.262, 0.268, 0.274, 0.280/)
     call addTransition2(thisIon, "2P_1/2", "4P_3/2", 1479.34 , 1.e-20,t, gamma, 12)

     !4-10
     gamma = (/ 0.102, 0.132, 0.215, 0.300, 0.331, 0.347, 0.359, 0.368, 0.375, 0.380, 0.390, 0.400/)
     call addTransition2(thisIon, "2P_3/2", "4P_3/2", 1478.43 , 1.e-20,t, gamma, 12)

     !5-6
     gamma = (/ 0.712, 0.714, 0.718, 0.727, 0.734, 0.739, 0.742, 0.744, 0.744, 0.742, 0.740, 0.737/)
     call addTransition2(thisIon, "2D_3/2", "4P_5/2", 3136.91 , 1.e-20,t, gamma, 12)

     !5-7
     gamma = (/ 1.97, 1.97, 1.98, 2.01, 2.02, 2.03, 2.03, 2.02, 2.02, 2.01, 1.99, 1.98/)
     call addTransition2(thisIon, "2D_5/2", "4P_5/2", 3137.16 , 1.e-20,t, gamma, 12)

     !5-8
     gamma = (/ 0.276, 0.256, 0.257, 0.261, 0.263, 0.266, 0.269, 0.272, 0.274, 0.276, 0.278, 0.279/)
     call addTransition2(thisIon, "2S_1/2", "4P_5/2", 1871.28 , 1.e-20,t, gamma, 12)

     !5-9
     gamma = (/ 0.045, 0.063, 0.116, 0.169, 0.188, 0.198, 0.206, 0.211, 0.217, 0.222, 0.227, 0.231/)
     call addTransition2(thisIon, "2P_1/2", "4P_5/2", 1479.96 , 1.e-20,t, gamma, 12)

     !5-10
     gamma = (/ 0.231, 0.288, 0.435, 0.582, 0.633, 0.657, 0.675, 0.690, 0.701, 0.710, 0.724, 0.738/)
     call addTransition2(thisIon, "2P_3/2", "4P_5/2", 1479.05 , 1.e-20,t, gamma, 12)

!     !7-6 
!     gamma = (/ 1.23, 1.15, 1.13, 1.16, 1.22, 1.27, 1.32, 1.37, 1.41, 1.45, 1.49, 1.52/)
!     call addTransition2(thisIon, "2D_3/2", "2D_5/2", 39739164.27 , 1.e-20,t, gamma, 12)
!
!     !6-8
!     gamma = (/ 0.680, 0.661, 0.660, 0.639, 0.624, 0.618, 0.615, 0.615, 0.616, 0.618, 0.620, 0.623/)
!     call addTransition2(thisIon, "2S_1/2", "2D_3/2", 4638.05 , 1.e-20,t, gamma, 12)

     deallocate(t, gamma)

  case("C III")
     allocate(t(4), gamma(4))
     t = (/ 5000, 10000, 15000, 20000/)
     gamma = (/ 0.611,  0.589,  0.578,  0.567/) 
     call addTransition2(thisIon, "3P^0_2", "1S_0", 1907., 5.15e-3, t, gamma, 4)
     gamma = (/ 0.367,  0.353,  0.347,  0.340/) 
     call addTransition2(thisIon, "3P^0_1", "1S_0", 1909., 104.0, t, gamma, 4)
     gamma = (/ 0.122,  0.118,  0.116,  0.113/) 
     call addTransition2(thisIon, "3P^0_0", "1S_0", 1909.6, 1.e-20, t, gamma, 4)
     gamma = (/ 3.58,  3.72,  3.81,  3.87/) 
     call addTransition2(thisIon, "1P^0_1", "1S_0", 977.02, 1.79e9, t, gamma, 4)
     gamma = (/ 2.23,  2.75,  3.05,  3.23/) 
     call addTransition2(thisIon, "3P^0_2", "3P^0_1", 1.774e6, 2.33e-6, t, gamma, 4)
     gamma = (/ 0.590,  0.695,  0.797,  0.855/) 
     call addTransition2(thisIon, "3P^0_2", "3P^0_0", 1.25e6, 1.37e-13, t, gamma, 4)
     gamma = (/ 0.800,  0.950,  1.01,  1.03/) 
     call addTransition2(thisIon, "3P^0_1", "3P^0_0", 4.22e6, 2.27e-7, t, gamma, 4)
     deallocate(t, gamma)


  case("C IV")
     allocate(t(2), gamma(2))

     t = (/ 10000, 20000/)
     gamma = (/ 2.96,  2.983/)
     call addTransition2(thisIon, "2P^0_1/2", "2S_1/2", 1550.8, 2.63E8, t, gamma, 2)
     gamma = (/ 5.92,  5.9667/)
     call addTransition2(thisIon, "2P^0_3/2", "2S_1/2", 1548.2, 2.65E8, t, gamma, 2)
     !	deallocate(t, gamma)


  case("N I") 
     allocate(t(15), gamma(15))
     
     t = (/ 5000, 6000, 7000, 8000, 9000, 10000, 12000, 13000, 14000, 15000, 16000, 17000, &
          18000, 19000, 20000/)
     gamma =  (/  0.1402 ,   .1732  ,  .2044  ,  .2338  ,  .2615  ,  .2877  ,  .3357  ,  .3576  ,  &
          .3784  ,  .3979  ,  .4164  ,  .4338  ,  .4503  ,  .4658  ,  .4804 /)
     call addTransition2(thisIon,"2D^0_5/2","4S^0_3/2", 5200.4, 7.27e-6, t, gamma, 15)
     gamma =  (/  0.0935 ,   .1155  ,  .1363  ,  .1559  ,  .1743  ,  .1918  ,  .2238  ,  .2384  ,  &
          .2523  ,  .2653  ,  .2776  ,  .2892  ,  .3002  ,  .3105  ,  .3203 /)
     call addTransition2(thisIon,"2D^0_3/2","4S^0_3/2", 5197.9, 2.02e-5, t, gamma, 15)
     gamma =  (/  0.1281 ,   .1591  ,  .1886  ,  .2167  ,  .2434  ,  .2609  ,  .3162  ,  .3382  ,  &
          .3591  ,  .3790  ,  .3979  ,  .4159  ,  .4331  ,  .4494  ,  .4650 /)
     call addTransition2(thisIon,"2D^0_3/2","2D^0_5/2", 1.1484e7, 1.27e-8, t, gamma, 15)
     gamma =  (/  0.0277 ,   .0342  ,  .0404  ,  .0463  ,  .0519  ,  .0573  ,  .0673  ,  .0720  ,  &
          .0765  ,  .0807  ,  .0848  ,  .0887  ,  .0924  ,  .0959  ,  .0993 /)
     call addTransition2(thisIon,"2P^0_1/2","4S^0_3/2", 3466.5, 2.71e-3, t, gamma, 15)
     gamma =  (/  0.0554 ,   .0683  ,  .0807  ,  .0925  ,  .1038  ,  .1146  ,  .1346  ,  .1440  ,  &
          0.1529  ,  .1615  ,  .1696  ,  .1774  ,  .1848  ,  .1919  ,  .1987 /)
     call addTransition2(thisIon,"2P^0_3/2","4S^0_3/2", 3466.5, 6.58e-3, t, gamma, 15)
     gamma =  (/  0.0626 ,   .0724  ,  .0819  ,  .0912  ,  .1002  ,  .1090  ,  .1261  , 0.1345  , &
          0.1427  , 0.1508  , 0.1588  , 0.1666  ,  .1744  , 0.1822  , 0.1898 /)
     call addTransition2(thisIon,"2P^0_1/2","2D^0_5/2", 10507.1, 3.45e-2, t, gamma, 15)
     gamma =  (/  0.1615 ,  0.1841  ,  .2058  ,  .2265  ,  .2466  ,  .2660  ,  .3033  ,  .3213  ,  &
          .3389  ,  .3562  ,  .3731  ,  .3898  ,  .4061  ,  .4223  ,  .4382 /)
     call addTransition2(thisIon,"2P^0_3/2","2D^0_5/2", 10397.7, 6.14e-2, t, gamma, 15)
     gamma =  (/  0.0601 ,   .0682  ,  .0758  ,  .0832  ,  .0902  ,  .0970  ,  .1100  ,  .1162  ,  &
          .1223  ,  .1283  ,  .1342  ,  .1399  ,  .1455  ,  .1510  ,  .1565 /)
     call addTransition2(thisIon,"2P^0_1/2","2D^0_3/2", 10407.6, 5.29e-2, t, gamma, 15)
     gamma =  (/  0.0856 ,   .0987  ,  .1113  ,  .1235  ,  .1354  ,  .1470  ,  .1695  ,  .1804  ,  &
          .1911  ,  .2017  ,  .2121  ,  .2224  ,  .2325  ,  .2425  ,  .2524 /)
     call addTransition2(thisIon,"2P^0_3/2","2D^0_3/2", 10407.2, 2.76e-2, t, gamma, 15)
     gamma =  (/  0.0329 ,   .0403  ,  .0478  ,  .0554  ,  .0632  ,  .0710  ,  .0869  ,  .0950  , &
          0.1031  , 0.1114  , 0.1196  , 0.1280  , 0.1363  , 0.1448  ,  .1533 /)
     call addTransition2(thisIon,"2P^0_3/2","2P^0_1/2", 2.59e8, 5.17e-13, t, gamma, 15)
     deallocate(t, gamma)
     
  case("N II")
     allocate( t(7), gamma(7))
     t = (/ 5000, 7500, 10000, 12500, 15000, 17500, 20000 /) 
     
     gamma = (/ .4097,  .4158,  .4232,  .4304,  .4371,  .4433,  .4491 /)
     call addTransition2(thisIon,"3P_1","3P_0", 2.055e6, 2.08e-6, t, gamma, 7)
     gamma = (/ .2439,  .2543,  .2657,  .2763,  .2859,  .2944,  .3019 /)
     call addTransition2(thisIon,"3P_2","3P_0", 7.65e5, 1.16e-12, t, gamma, 7)
     gamma =  (/ .3385,  .3357,  .3342,  .3333,  .3328,  .3324,  .3323 /)
     call addTransition2(thisIon,"1D_2","3P_0", 6529.0, 5.35e-7, t, gamma, 7)
     gamma = (/ .0424,  .0415,  .0410,  .0408,  .0406,  .0405,  .0405 /)
     call addTransition2(thisIon,"1S_0","3P_0", 3062.9, 1.e-20, t, gamma, 7)
     gamma = (/ .1282,  .1273,  .1274,  .1277,  .1280,  .1283,  .1285 /)
     call addTransition2(thisIon,"5S^0_2","3P_0", 2144., 1.e-20, t, gamma, 7)
     gamma = (/1.0608, 1.0918, 1.1268, 1.1598, 1.1898, 1.2167, 1.2408 /)
     call addTransition2(thisIon,"3P_2","3P_1", 1.22e6, 7.46e-6, t, gamma, 7)
     gamma = (/1.0157, 1.0069, 1.0024,  .9998,  .9982,  .9973,  .9967 /)
     call addTransition2(thisIon,"1D_2","3P_1", 6548.1, 1.01e-3, t, gamma, 7)
     gamma = (/ .1271,  .1245,  .1231,  .1223,  .1219,  .1216,  .1214 /)
     call addTransition2(thisIon,"1S_0","3P_1", 3062.9, 3.38e-2, t, gamma, 7)
     gamma =(/ .3849,  .3817,  .3819,  .3828,  .3839,  .3848,  .3856 /)
     call addTransition2(thisIon,"5S^0_2","3P_1", 2144., 1.00e2, t, gamma, 7)
     gamma = (/1.6929, 1.6782, 1.6707, 1.6663, 1.6637, 1.6621, 1.6612 /)
     call addTransition2(thisIon,"1D_2","3P_2", 6583.4, 2.99e-3, t, gamma, 7)
     gamma = (/ .2118,  .2074,  .2052,  .2039,  .2031,  .2026,  .2023 /)
     call addTransition2(thisIon,"1S_0","3P_2", 3071.4, 1.51e-4, t, gamma, 7)
     gamma = (/ .6412,  .6366,  .6369,  .6383,  .6399,  .6413,  .6424 /)
     call addTransition2(thisIon,"5S^0_2","3P_2", 2140., 1.00e2, t, gamma, 7)
     gamma = (/ .5027,  .5001,  .5013,  .5040,  .5071,  .5105,  .5140 /)
     call addTransition2(thisIon,"1S_0","1D_2", 5754.6, 1.12e0, t, gamma, 7)
     deallocate(t, gamma)

  case("N III")
     allocate(t(11), gamma(11))
     t = (/ 2000.,  4000.,   6000.,  8000.,  10000., 12000., 14000., 16000., 18000., 20000., 50000./)
     gamma = (/1.051, 1.011, 1.023, 1.050, 1.080, 1.105, 1.126, 1.144, 1.158, 1.168, 1.168 /)
     call addTransition2(thisIon,"2P^0_3/2","2P^0_1/2", 5.73e5, 4.77e-5, t, gamma, 11)
     gamma = (/.202, .202, .202, .202, .204, .206, .208, .210, .212, .214, .214 /)
     call addTransition2(thisIon,"4P_1/2","2P^0_1/2", 1748., 4.914e2, t, gamma, 11)
     gamma = (/.295, .295, .295, .295, .298, .301, .304, .307, .310, .312, .312 /)
     call addTransition2(thisIon,"4P_3/2","2P^0_1/2", 1747., 9.81, t, gamma, 11)
     gamma = (/.170, .170, .170, .170, .172, .173, .175, .176, .178, .179, .179 /)
     call addTransition2(thisIon,"4P_5/2","2P^0_1/2", 1744., 1.e-20, t, gamma, 11)
     gamma = (/.131, .131, .131, .131, .132, .133, .134, .135, .136, .138, .138 /)
     call addTransition2(thisIon,"4P_1/2","2P^0_3/2", 1754., 5.276e2, t, gamma, 11)
     gamma = (/.372, .372, .372, .372, .375, .378, .381, .385, .388, .392, .392 /)
     call addTransition2(thisIon,"4P_3/2","2P^0_3/2", 1752., 64.69, t, gamma, 11)
     gamma = (/.829, .829, .829, .829, .837, .844, .852, .859, .867, .875, .875 /)
     call addTransition2(thisIon,"4P_5/2","2P^0_3/2", 1747., 3.08e2, t, gamma, 11)
     gamma = (/.695, .695, .695, .695, .695, .695, .695, .695, .695, .695, .695 /)
     call addTransition2(thisIon,"4P_3/2","4P_1/2", 1.68e6, 1.e-20, t, gamma, 11)
     gamma = (/.397, .397, .397, .397, .397, .397, .397, .397, .397, .397, .397 /)
     call addTransition2(thisIon,"4P_5/2","4P_1/2", 7.10e5, 1.e-20, t, gamma, 11)
     gamma = (/1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26, 1.26 /)
     call addTransition2(thisIon,"4P_5/2","4P_3/2", 1.23e6, 1.e-20, t, gamma, 11)
     deallocate(t, gamma)

!  case("N IV")
!     allocate(t(4), gamma(4))
!     t = (/5000, 10000, 15000, 20000 /)
!     gamma = (/0.937, 0.905, 0.879, 0.858 /)
!     call addTransition2(thisIon,"3P^0_2","1S_0", 1483.3, 1.15e-2, t, gamma, 11)
!     deallocate(t, gamma)


  case("O I")
     allocate(t(5), gamma(5))
     t = (/ 500,  1000, 5000, 10000, 20000 /)
     gamma =  (/  5.69E-3   ,     9.72E-3   ,    4.52E-2   ,    0.106    ,     2.066E-1 /)
     call addTransition2(thisIon,"3P_1","3P_2", 6.32e5, 8.957e-5, t, gamma, 5)
     gamma =  (/  2.51E-3   ,     4.11E-3   ,    1.53E-2   ,    3.21E-2  ,     5.36E-2 /)
     call addTransition2(thisIon,"3P_0","3P_2", 4.41e5, 1.203e-10, t, gamma, 5)
     gamma =  (/  0.3222E-02,     0.8389E-02,    0.6889E-01,    0.1478E+00,    0.2783E+00  /)
     call addTransition2(thisIon,"1D_2","3P_2", 6300.3, 5.627e-3, t, gamma, 5)
     gamma =  (/  0.3611E-03,     0.1022E-02,    0.8500E-02,    0.1800E-01,    0.3372E-01  /)
     call addTransition2(thisIon,"1S_0","3P_2", 2959.2, 2.732e-4, t, gamma, 5)
     gamma =  (/  2.86E-4   ,     6.36E-4   ,    9.12E-3   ,    2.83E-2  ,     6.93E-2     /)
     call addTransition2(thisIon,"3P_0","3P_1", 1.46e6, 1.735e-5, t, gamma, 5)
     gamma =  (/  0.1933E-02,     0.5033E-02,    0.4133E-01,    0.8867E-01,    0.1670E+00  /)
     call addTransition2(thisIon,"1D_2","3P_1", 6363.9, 1.818e-3, t, gamma, 5)
     gamma =  (/  0.2167E-03,     0.6133E-03,    0.5100E-02,    0.1080E-01,    0.2023E-01  /)
     call addTransition2(thisIon,"1S_0","3P_1", 2972.3, 7.601e-2, t, gamma, 5)
     gamma =  (/  0.6444E-03,     0.1678E-02,    0.1378E-01,    0.2956E-01,    0.5567E-01  /)
     call addTransition2(thisIon,"1D_2","3P_0", 6393.5, 8.922e-7, t, gamma, 5)
     gamma =  (/  0.7222E-04,     0.2044E-03,    0.1700E-02,    0.3600E-02,    0.6744E-02  /)
     call addTransition2(thisIon,"1S_0","3P_0", 2979. , 1.e-20, t, gamma, 5)
     gamma =  (/  2.1E-2    ,     3.1E-2    ,    7.32E-2   ,   1.05E-1    , 1.48E-1      /)
     call addTransition2(thisIon,"1S_0","1D_2", 5577.3, 1.215e0, t, gamma, 5)
     deallocate(t, gamma)

  case("O II")
     allocate(t(11), gamma(11))
     t = (/5000, 6000, 7000, 8000, 9000, 10000, 12000, 14000, 16000, 18000, 20000/)
     gamma = (/0.7950, 0.7956, 0.7968, 0.7986, 0.7998, 0.8010, 0.8040, 0.8076, 0.8106, 0.8142, 0.8178/)
     call addTransition2(thisIon,"2D^0_5/2","4S^0_3/2", 3728.8, 3.82e-5, t, gamma, 11)
     gamma = (/0.5300, 0.5304, 0.5312, 0.5324, 0.5332, 0.5340, 0.5360, 0.5384, 0.5404, 0.5428, 0.5452/)
     call addTransition2(thisIon,"2D^0_3/2","4S^0_3/2", 3726.0, 1.65e-4, t, gamma, 11)
     gamma = (/0.26510  ,  .2657   ,  .2671   ,  .2684  ,   .2684  ,  .2704  ,  .2724  , .2744  , .2764  ,  .2784 ,   .2804 /)
     call addTransition2(thisIon,"2P^0_3/2", "4S^0_3/2", 2470.3, 5.640e-2, t, gamma, 11)
     gamma = (/  0.1325 ,    .1329 ,    .1335 ,    .1342,     .1345,    .1352,    .1362,   .1372,   .1382,    .1392,    .1402 /)
     call addTransition2(thisIon,"2P^0_1/2", "4S^0_3/2", 2470.2, 2.32e-2, t, gamma, 11)
     gamma = (/1.22     , 1.208    , 1.197    , 1.186   ,  1.177   , 1.168   , 1.152   ,1.139   ,1.128   , 1.119  ,  1.112  /)
     call addTransition2(thisIon,"2D^0_3/2", "2D^0_5/2", 4.97e6, 1.20e-7, t, gamma, 11)
     gamma = (/  .28    ,   .282   ,   .283   ,   .284  ,    .286  ,   .287  ,   .290  ,  .292  ,  .295  ,   .297 ,    .3    /)
     call addTransition2(thisIon,"2P^0_1/2", "2P^0_3/2", 5.00e7, 2.08e-11, t, gamma, 11)
     gamma = (/  .7178  ,   .7203  ,   .7229  ,   .7250 ,    .7276 ,   .7302 ,   .7353 ,  .7404 ,  .7451 ,   .7503,    .7554 /)
     call addTransition2(thisIon,"2P^0_3/2", "2D^0_5/2", 7319.9, 1.17e-1, t, gamma, 11)
     gamma = (/  .4008  ,   .4022  ,   .4037  ,   .4049 ,    .4063 ,   .4077 ,   .4106 ,  .4135 ,  .4161 ,   .4190,    .4218 /)
     call addTransition2(thisIon,"2P^0_3/2", "2D^0_3/2", 7330.7, 6.14e-2, t, gamma, 11)
     gamma = (/  .2901  ,   .2912  ,   .2922  ,   .2931 ,    .2941 ,   .2951 ,   .2972 ,  .2993 ,  .3012 ,   .3033,    .3053 /)
     call addTransition2(thisIon,"2P^0_1/2", "2D^0_5/2", 7321.8, 6.15e-2, t, gamma, 11)
     gamma = (/  .2700  ,   .2710  ,   .2719  ,   .2727 ,    .2737 ,   .2747 ,   .2766 ,  .2785 ,  .2803 ,   .2822,    .2842 /)
     call addTransition2(thisIon,"2P^0_1/2", "2D^0_3/2", 7329.6, 1.02e-1, t, gamma, 11)
     deallocate(t, gamma)

  case("O III")
     allocate(t(1:13), gamma(1:13))
     t = (/2500,   5000,   7500,  10000,  12500,  15000,  17500,  20000,  25000,  30000,  40000,  50000,  60000/) 
     gamma = (/0.5041, 0.5172, 0.5310, 0.5417, 0.5490, 0.5537, 0.5567, 0.5586, 0.5612, 0.5633, 0.5670, 0.5699, 0.5715/)
     call addTransition2(thisIon,"3P_1","3P_0", 883562., 2.62e-5, t, gamma, 13)
     gamma = (/0.2499, 0.2566, 0.2646, 0.2717, 0.2776, 0.2824, 0.2865, 0.2901, 0.2962, 0.3013, 0.3097, 0.3160, 0.3203/)
     call addTransition2(thisIon,"3P_2","3P_0", 326611., 3.02e-11, t, gamma, 13)
     gamma = (/0.2283, 0.2262, 0.2337, 0.2426, 0.2506, 0.2627, 0.2627, 0.2672, 0.2740, 0.2790, 0.2850, 0.2880, 0.2892/)
     call addTransition2(thisIon,"1D_2","3P_0", 4932.6, 2.74e-6, t, gamma, 13)
     gamma = (/0.0278, 0.0280, 0.0295, 0.0310, 0.0324, 0.0335, 0.0344, 0.0351, 0.0362, 0.0368, 0.0375, 0.0377, 0.0377/)
     call addTransition2(thisIon,"1S_0","3P_0", 2315.6, 1.e-20, t, gamma, 13)
     gamma = (/0.1104, 0.1223, 0.1334, 0.1398, 0.1434, 0.1454, 0.1465, 0.1472, 0.1474, 0.1467, 0.1439, 0.1402, 0.1363/)
     call addTransition2(thisIon,"5S^0_2","3P_0", 1657.7, 1.e-20, t, gamma, 13)
     gamma = (/1.1925, 1.2239, 1.2592, 1.2884, 1.3107, 1.3275, 1.3404, 1.3510, 1.3679, 1.3821, 1.4056, 1.4232, 1.4350/)
     call addTransition2(thisIon,"3P_2","3P_1", 518145., 9.76e-5, t, gamma, 13)
     gamma = (/0.6848, 0.6785, 0.7010, 0.7279, 0.7518, 0.7716, 0.7879, 0.8014, 0.8221, 0.8368, 0.8551, 0.8641, 0.8675/)
     call addTransition2(thisIon,"1D_2","3P_1", 4958.9, 6.74e-3, t, gamma, 13)
     gamma = (/0.0833, 0.0840, 0.0884, 0.0931, 0.0972, 0.1006, 0.1033, 0.1054, 0.1085, 0.1105, 0.1126, 0.1132, 0.1131/)
     call addTransition2(thisIon,"1S_0","3P_1", 2321.0, 2.23e-1, t, gamma, 13)
     gamma = (/0.3311, 0.3669, 0.4003, 0.4194, 0.4301, 0.4361, 0.4396, 0.4415, 0.4420, 0.4400, 0.4316, 0.4206, 0.4088/)
     call addTransition2(thisIon,"5S^0_2","3P_1", 1660.8, 2.12e2, t, gamma, 13)
     gamma = (/1.1413, 1.1308, 1.1683, 1.2131, 1.2529, 1.2860, 1.3132, 1.3357, 1.3702, 1.3947, 1.4252, 1.4402, 1.4458/)
     call addTransition2(thisIon,"1D_2","3P_2", 5006.7, 1.96e-2, t, gamma, 13)
     gamma = (/0.1388, 0.1401, 0.1473, 0.1552, 0.1620, 0.1676, 0.1721, 0.1757, 0.1809, 0.1842, 0.1876, 0.1886, 0.1884/)
     call addTransition2(thisIon,"1S_0","3P_2", 2332.0, 7.85e-4, t, gamma, 13)
     gamma = (/0.5519, 0.6114, 0.6672, 0.6991, 0.7168, 0.7268, 0.7326, 0.7357, 0.7367, 0.7333, 0.7194, 0.7010, 0.6814/)
     call addTransition2(thisIon,"5S^0_2","3P_2", 1666.1, 5.22e2, t, gamma, 13)
     gamma = (/0.4708, 0.5463, 0.6114, 0.6468, 0.6630, 0.6687, 0.6692, 0.6670, 0.6599, 0.6524, 0.6407, 0.6330, 0.6277/)
     call addTransition2(thisIon,"1S_0","1D_2", 4363.2, 1.78e0, t, gamma, 13)
     deallocate(t, gamma)

  case("Ne II")
     call addTransition(thisIon,"2P^0_1/2","2P^0_3/2", 1.28e5, 8.55e-3, 2.96e-1, 3.03e-1, 3.10e-1, 3.17e-1)
     
  case("Ne III")
     allocate(t(11), gamma(11))
     t = (/ 1000., 1585., 2512., 3981., 6310., 10000., 15849., 25119., 39811., 63096., 100000. /)
     gamma = (/ .154, .168, .194, .218, .235, .244, .247, .247, .246, .249, .256 /)
     call addTransition2(thisIon,"3P_0", "3P_1", 3.60e5, 1.15e-3, t, gamma, 11)
     gamma = (/ .128, .149, .174, .194, .204, .208, .208, .205, .203, .205, .211 /)
     call addTransition2(thisIon,"3P_0", "3P_2", 1.07e5, 2.18e-8, t, gamma, 11)
     gamma = (/ .481, .545, .634, .708, .752, .774, .778, .771, .764, .773, .794 /)
     call addTransition2(thisIon,"3P_1", "3P_2", 1.56e5, 5.97e-3, t, gamma, 11)
     gamma = (/ .150, .153, .154, .153, .152, .151, .150, .149, .150, .153, .156 /)
     call addTransition2(thisIon,"1D_2","3P_0", 4012.8, 8.51e-6, t, gamma, 11)
     gamma = (/ .450, .459, .462, .460, .456, .452, .449, .448, .449, .458, .469 /)
     call addTransition2(thisIon,"1D_2","3P_1", 3967.5, 5.42e-2, t, gamma, 11)
     gamma = (/ .749, .765, .771, .767, .760, .754, .749, .746, .749, .763, .782 /)
     call addTransition2(thisIon,"1D_2","3P_2", 3868.8, 1.71e-1, t, gamma, 11)
     gamma = (/ .017, .017, .017, .017, .017, .017, .017, .017, .018, .019, .020 /)
     call addTransition2(thisIon,"1S_0","3P_0", 1814.6, 1.e-20, t, gamma, 11)
     gamma = (/ .050, .050, .050, .050, .050, .050, .050, .051, .054, .057, .060 /)
     call addTransition2(thisIon,"1S_0","3P_1", 1814.6, 2.00e0,  t, gamma, 11)
     gamma = (/ .083, .083, .083, .083, .084, .084, .084, .086, .090, .095, .100 /)
     call addTransition2(thisIon,"1S_0","3P_2", 1793.7, 3.94e-3, t, gamma, 11)
     gamma = (/ .266, .266, .266, .266, .267, .269, .277, .292, .310, .325, .333 /)
     call addTransition2(thisIon,"1S_0","1D_2", 3342.5, 2.71e0, t, gamma, 11)
     deallocate(t, gamma)

  case("S II")
     allocate(t(16), gamma(16))
     t = (/5000., 6000., 7000., 8000., 9000.,10000.,11000.,12000.,13000.,14000.,15000.,16000.,17000.,18000.,19000.,20000. /)
     gamma = (/ 2.9228,  2.8978,  2.871,   2.8438,  2.8167,  2.7902,  2.7644,  2.7394,  2.7151,  2.6914,  2.6683,  &
          2.6456,  2.6234,  2.6014,  2.5796,  2.5581 /) 
     call addTransition2(thisIon,"2D^0_3/2", "4S^0_3/2", 6730.8, 8.82e-4, t, gamma, 16)
     gamma = (/   4.3842,  4.3467,  4.3065,  4.2656,  4.2250,  4.1853,  4.1466,  4.1091,  4.0726,  4.0371,  4.0024, &
          3.9684,  3.935 ,  3.9021,  3.8694,  3.8371 /) 
     call addTransition2(thisIon,"2D^0_5/2", "4S^0_3/2", 6716.5, 2.60e-4, t, gamma, 16)
     gamma = (/   8.1464,  8.0429,  7.9304,  7.8149,  7.7006,  7.5897,  7.4834,  7.3819,  7.2853,  7.1933,  7.1053, &
          7.0211,  6.9401,  6.862 ,  6.7864,  6.713  /) 
     call addTransition2(thisIon,"2D^0_5/2", "2D^0_3/2", 3.145e6, 3.35e-7, t, gamma, 16)
     gamma = (/    .7223,   .7329,   .7414,   .7484,   .7543,   .7592,   .7633,   .7666,   .7691,   .7709,   .7721, &
          .7726,   .7725,   .7718,   .7706,   .7688 /) 
     call addTransition2(thisIon,"2P^0_1/2", "4S^0_3/2", 4076.4, 9.06e-2, t, gamma, 16)
     gamma = (/   1.4446,  1.4658,  1.4828,  1.4968,  1.5086,  1.5184,  1.5266,  1.5331,  1.5382,  1.5419,  1.5442, &
          1.5452,  1.5450,  1.5436,  1.5411,  1.5376 /) 
     call addTransition2(thisIon,"2P^0_3/2", "4S^0_3/2", 4068.6, 2.25e-1, t, gamma, 16)
     gamma = (/   1.499 ,  1.5086,  1.515,   1.5192,  1.5219,  1.5232,  1.5234,  1.5226,  1.5208,  1.5180,  1.5144, &
          1.5099,  1.5046,  1.4985,  1.4918,  1.4844 /) 
     call addTransition2(thisIon,"2P^0_1/2", "2D^0_3/2", 10336.3, 1.63e-1, t, gamma, 16)
     gamma = (/   3.3615,  3.3773,  3.3850,  3.3872,  3.3851,  3.3798,  3.3716,  3.3610,  3.3482,  3.3335,  3.317 , &
          3.2988,  3.2793,  3.2584,  3.2363,  3.2132 /) 
     call addTransition2(thisIon,"2P^0_3/2", "2D^0_3/2", 10286.7, 1.33e-1, t, gamma, 16)
     gamma = (/   2.5514,  2.563 ,  2.5683,  2.5694,  2.5673,  2.5626,  2.5558,  2.5471,  2.5367,  2.5249,  2.5117, &
          2.4974,  2.482 ,  2.4656,  2.4483,  2.4303 /) 
     call addTransition2(thisIon,"2P^0_1/2", "2D^0_5/2", 10373.3, 7.79e-2, t, gamma, 16)
     gamma = (/   4.7394,  4.7658,  4.7816,  4.7902,  4.7932,  4.7919,  4.7868,  4.7783,  4.7668,  4.7523,  4.7352, &
          4.7157,  4.6938,  4.6698,  4.6439,  4.6162 /) 
     call addTransition2(thisIon,"2P^0_3/2", "2D^0_5/2", 10320.4, 1.79e-1, t, gamma, 16)
     gamma = (/   2.3915,  2.3948,  2.3937,  2.39000, 2.3846,  2.378 ,  2.3705,  2.3621,  2.3529,  2.3428,  2.3320, &
          2.3204,  2.3081,  2.2950,  2.2812,  2.2667 /) 
     call addTransition2(thisIon,"2P^0_3/2", "2P^0_1/2", 2.14e6, 1.03-6, t, gamma, 16)
     deallocate(t, gamma)



!      case("S III") !collision strength values from Tayal and Gupta 1999 (haworth)
!	allocate(t(6), gamma(6))
!	t = (/5000., 8000., 10000., 20000., 30000., 40000./)
!	gamma = (/  4.44  , 4.17  , 3.98  ,  3.24 , 2.85 , 2.61  /)!1 -2
!        call addTransition2(thisIon,"3P_1", "3P_0", 3.347e5, 4.72e-4, t, gamma, 4) !1 -2
!        gamma = (/  1.41  , 1.35  , 1.31  , 1.28 , 1.31 , 1.33  /)!1-3
!        call addTransition2(thisIon,"3P_2", "3P_0", 1.20e5, 4.61e-8, t, gamma, 4)!1-3
!        gamma = (/  0.802  , 0.778 , 0.773 , 0.786 , 0.792 , 0.785  /)!1-4
!        call addTransition2(thisIon,"1D_2", "3P_0", 8833.9, 5.82e-6, t, gamma, 4)!1-4
!        gamma = (/  0.129 , 0.130 , 0.131 ,  0.137 , 0.144 , 0.148  /)!1-5
!        call addTransition2(thisIon,"1S_0", "3P_0", 3681.8, 1.e-20, t, gamma, 4)!1-5
!        gamma = (/  8.72  , 8.20  , 7.87  ,  6.88 , 6.46 , 6.21  /)!2-3
!        call addTransition2(thisIon,"3P_2", "3P_1", 187129., 2.07e-3, t, gamma, 4)!2-3
!        gamma = (/  2.41  , 2.33  , 2.32  ,  2.36 , 2.38 , 2.35  /)!2-4
!        call addTransition2(thisIon,"1D_2", "3P_1", 9068.9, 2.21e-2, t, gamma, 4)
!        gamma = (/  0.388 , 0.391 , 0.393 ,  0.410, 0.432, 0.445  /)!2-5
!        call addTransition2(thisIon,"1S_0", "3P_1", 3721.7, 7.96e-1, t, gamma, 4)!2-5
!        gamma = (/  4.01  , 3.89  , 3.86  ,  3.93 , 3.96 , 3.92   /)!3-4
!        call addTransition2(thisIon,"1D_2", "3P_2", 9531.0, 5.76e-2, t, gamma, 4)!3-4
!        gamma = (/  0.646 , 0.651 , 0.655 ,  0.683 , 0.719 , 0.742 /)!3-5
!        call addTransition2(thisIon,"1S_0", "3P_2", 3797.8, 1.05e-2, t, gamma, 4)!3-5
!        gamma = (/  1.31  , 1.33 , 1.38  ,  1.59 , 1.72 , 1.79  /)!4-5
!        call addTransition2(thisIon,"1S_0", "1D_2", 6312.1, 2.22e0, t, gamma, 4)!4-5
!	deallocate(t, gamma)

!Old data
!     case("S III")
!        allocate(t(4), gamma(4))
!        t(1:4) = (/5000., 10000., 15000., 20000./)
!        gamma = (/  2.64  , 2.59  , 2.38  ,  2.20  /)
!        call addTransition2(thisIon,"3P_1", "3P_0", 3.347e5, 4.72e-4, t, gamma, 4)
!        gamma = (/  1.11  , 1.15  , 1.15  ,  1.14  /)
!        call addTransition2(thisIon,"3P_2", "3P_0", 1.20e5, 4.61e-8, t, gamma, 4)
!        gamma = (/  1.01  , 0.932 , 0.921 ,  0.911 /)
!        call addTransition2(thisIon,"1D_2", "3P_0", 8833.9, 5.82e-6, t, gamma, 4)
!        gamma = (/  0.129 , 0.132 , 0.134 ,  0.138 /)
!        call addTransition2(thisIon,"1S_0", "3P_0", 3681.8, 1.e-20, t, gamma, 4)
!        gamma = (/  5.79  , 5.81  , 5.56  ,  5.32  /)
!        call addTransition2(thisIon,"3P_2", "3P_1", 187129., 2.07e-3, t, gamma, 4)
!        gamma = (/  3.02  , 2.80  , 2.76  ,  2.73  /)
!        call addTransition2(thisIon,"1D_2", "3P_1", 9068.9, 2.21e-2, t, gamma, 4)
!        gamma = (/  0.387 , 0.397 , 0.403 ,  0.413 /)
!        call addTransition2(thisIon,"1S_0", "3P_1", 3721.7, 7.96e-1, t, gamma, 4)
!        gamma = (/  5.04  , 4.66  , 4.61  ,  4.56  /)
!        call addTransition2(thisIon,"1D_2", "3P_2", 9531.0, 5.76e-2, t, gamma, 4)
!        gamma = (/  0.644 , 0.661 , 0.672 ,  0.689 /)
!        call addTransition2(thisIon,"1S_0", "3P_2", 3797.8, 1.05e-2, t, gamma, 4)
!        gamma = (/  1.42  , 1.88  , 2.02  ,  2.08  /)
!        call addTransition2(thisIon,"1S_0", "1D_2", 6312.1, 2.22e0, t, gamma, 4)
!        deallocate(t, gamma)


!Matching Pradhan as closely as possible
!     case("S III")
!        allocate(t(4), gamma(4))
!        t(1:4) = (/5000., 10000., 15000., 20000./)
!        gamma = (/  2.64  , 2.59  , 2.38  ,  2.20  /)
!        call addTransition2(thisIon,"3P_1", "3P_0", 3.347e5, 4.72e-4, t, gamma, 4)
!        gamma = (/  1.11  , 1.15  , 1.15  ,  1.14  /)
!        call addTransition2(thisIon,"3P_2", "3P_0", 1.20e5, 4.61e-8, t, gamma, 4)
!        gamma = (/  9.07  , 8.39 , 8.29 ,  8.20 /)!!!!!!!!
!        call addTransition2(thisIon,"1D_2", "3P_0", 8833.9, 5.82e-6, t, gamma, 4)
!        gamma = (/  0.129 , 0.132 , 0.134 ,  0.138 /)
!        call addTransition2(thisIon,"1S_0", "3P_0", 3681.8, 1.e-20, t, gamma, 4)
!        gamma = (/  5.79  , 5.81  , 5.56  ,  5.32  /)
!        call addTransition2(thisIon,"3P_2", "3P_1", 187129., 2.07e-3, t, gamma, 4)
!        gamma = (/  3.02  , 2.80  , 2.76  ,  2.73  /)
!        call addTransition2(thisIon,"1D_2", "3P_1", 9068.9, 2.21e-2, t, gamma, 4)
!        gamma = (/  1.16 , 1.19 , 1.21 ,  1.24 /)!!!!!!!!!
!        call addTransition2(thisIon,"1S_0", "3P_1", 3721.7, 7.96e-1, t, gamma, 4)
!        gamma = (/  5.04  , 4.66  , 4.61  ,  4.56  /)
!        call addTransition2(thisIon,"1D_2", "3P_2", 9531.0, 5.76e-2, t, gamma, 4)
!        gamma = (/  0.644 , 0.661 , 0.672 ,  0.689 /)
!        call addTransition2(thisIon,"1S_0", "3P_2", 3797.8, 1.05e-2, t, gamma, 4)
!        gamma = (/  1.42  , 1.88  , 2.02  ,  2.08  /)
!        call addTransition2(thisIon,"1S_0", "1D_2", 6312.1, 2.22e0, t, gamma, 4)
!        deallocate(t, gamma)

!Matching Barbara Ercolano's S III data as closely as possible
  case("S III")
     allocate(t(16), gamma(16))
     t(1:16) = (/ 100., 158.5, 251.2, 398.1, 631.0, 1000.0, 1584.9, 2511.9, 3981.1, 6309.6, 10000.0, 15848.9, 25118.9, &
          39810.7, 63095.8, 100000.0 /)
     gamma = (/  4.78  , 4.78  , 4.77  ,  4.76 , 4.75 , 4.72 , 4.69 , 4.62 , 4.52 , 4.33 , 3.98 , 3.49 , 3.01 , &
          2.62 , 2.25 , 1.90 /)
     call addTransition2(thisIon,"3P_1", "3P_0", 3.347e5, 8.29e-4, t, gamma, 16) !1-2
     gamma = (/  1.61  , 1.60  , 1.60  ,  1.59 , 1.58 , 1.56 , 1.54 , 1.50 , 1.44 , 1.38 , 1.32 , 1.28 , 1.29 , &
          1.33 , 1.32 , 1.24 /)
     call addTransition2(thisIon,"3P_2", "3P_0", 1.20e5, 5.008e-8, t, gamma, 16) !1-3
     gamma = (/  0.901  , 0.900 , 0.897 ,  0.893 , 0.887 , 0.877 , 0.863 , 0.842 , 0.816 , 0.789 , 0.773 , 0.778 , &
          0.791 , 0.786 , 0.741 , 0.665 /)
     call addTransition2(thisIon,"1D_2", "3P_0", 8833.9, 2.570e-4, t, gamma, 16) !1-4
     gamma = (/  0.128 , 0.128 , 0.128 , 0.128 , 0.128 , 0.128 , 0.128 , 0.128 , 0.129 , 0.129 , 0.131 , 0.134 , &
          0.140 , 0.148 , 0.148 , 0.136 /)
     call addTransition2(thisIon,"1S_0", "3P_0", 3681.8, 1.e-20, t, gamma, 16) !1-5
     gamma = (/  9.31  , 9.31  , 9.30  ,  9.28 , 9.26 , 9.22 , 9.16 , 9.06 , 8.87 , 8.50 , 7.87 , 7.15 , 6.64 , &
          6.23 , 5.71 , 5.15 /)
     call addTransition2(thisIon,"3P_2", "3P_1", 187129., 1.877e-3, t, gamma, 16) !2-3
     gamma = (/  2.71  , 2.70  , 2.70  ,  2.68 , 2.67 , 2.64 , 2.59 , 2.53 , 2.45 , 2.37 , 2.32 , 2.33 , 2.38 , &
          2.36 , 2.22 , 1.99 /)
     call addTransition2(thisIon,"1D_2", "3P_1", 9068.9, 3.107e-2, t, gamma, 16) !2-4
     gamma = (/  0.371 , 0.371 , 0.372 ,  0.372 , 0.373 , 0.375 , 0.377 , 0.381 , 0.385 , 0.390 , 0.393 , 0.400 , &
          0.422 , 0.447 , 0.442 , 0.408 /)
     call addTransition2(thisIon,"1S_0", "3P_1", 3721.7, 1.016e0, t, gamma, 16) !2-5
     gamma = (/  4.39  , 4.38  , 4.37  ,  4.36 , 4.34 , 4.30 , 4.25 , 4.17 , 4.07 , 3.95 , 3.86 , 3.88 , 3.95 , &
          3.93 , 3.71 , 3.32  /)
     call addTransition2(thisIon,"1D_2", "3P_2", 9531.0, 1.856e-1, t, gamma, 16) !3-4
     gamma = (/  0.620 , 0.620 , 0.621 ,  0.622,  0.623,  0.626,  0.630,  0.635,  0.642,  0.649,  0.655,  0.666 , &
          0.702 ,  0.744 ,  0.737 ,  0.680 /)
     call addTransition2(thisIon,"1S_0", "3P_2", 3797.8, 1.408e-3, t, gamma, 16) !3-5
     gamma = (/  1.40  , 1.40  , 1.39  ,  1.39 ,  1.38 ,  1.37 ,  1.35 ,  1.33 ,  1.31 , 1.31 , 1.38 , 1.51 , &
          1.67 , 1.79 , 1.85 , 1.84 /)
     call addTransition2(thisIon,"1S_0", "1D_2", 6312.1, 3.045e0, t, gamma, 16) !4-5
     deallocate(t, gamma)
     
  case("S IV")
     allocate(t(15), gamma(15))
     t = (/2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000, 30000, 40000, 50000, 100000, 200000/)
     gamma = (/ 7.109, 7.85, 8.27, 8.47, 8.548, 8.55, 8.51, 8.44, 8.36, 8.275, 7.837, 7.467, 7.15, 5.918, 4.41 /)
     call addTransition2(thisIon,"2P^0_3/2", "2P^0_1/2", 105000., 7.74e-3, t, gamma, 15)
     deallocate(t, gamma)


  case DEFAULT
  end select

 end subroutine addTransitions

subroutine addTransition(thisIon,term1, term2, lambda, a, gamma1, gamma2, gamma3, gamma4)
  type(IONTYPE) :: thisIon
  character(len=*) :: term1, term2
  real :: lambda
  real :: a
  real :: gamma1, gamma2, gamma3, gamma4
  real :: gamma(4)
  real :: thisLam
  integer :: i, j, k, m
  real :: t(4) = (/5000., 10000., 15000., 20000./)
  character(len=80) :: message

  gamma(1) = gamma1
  gamma(2) = gamma2
  gamma(3) = gamma3
  gamma(4) = gamma4

  k = thisIon%nTransitions + 1
  thisIon%nTransitions = k

  j = returnLevel(thisIon,term1)
  i = returnLevel(thisIon,term2)

  thisIon%transition(k)%i = i
  thisIon%transition(k)%j = j
  thisLam = real(1.e8*cSpeed / ((thisIon%level(j)%energy - thisIon%level(i)%energy)/ergToEv/hCgs))

  thisIon%transition(k)%energy = (thisIon%level(j)%energy - thisIon%level(i)%energy)
  if (thisIon%transition(k)%energy < 0.) then
     write(*,*) "negative energy for ",thisIon%species,term1,term2
     stop
  endif


  if (abs(thisLam-lambda)/lambda > 0.05) then
     call writeWarning("WARNING! Given wavelength and  calculated wavelength differ by more than 5%")
     write(message,*) "Calculated: ",thisLam, "Given: ", lambda
     call writeWarning(message)
     write(message,*) "For transition: ",trim(term1)," ",trim(term2), " in ion ",thisIon%species
     call writeWarning(message)
  endif
  thisIon%transition(k)%lambda = lambda
  thisIon%transition(k)%a = a


  
  thisIon%transition(k)%nGamma = 0
  do m = 1, 4
     if (gamma(m) /= 0.) then
        thisIon%transition(k)%nGamma =  thisIon%transition(k)%nGamma + 1
        thisIon%transition(k)%gamma(thisIon%transition(k)%nGamma) = gamma(m)
        thisIon%transition(k)%t(thisIon%transition(k)%nGamma) = t(m)
     endif
  enddo
end subroutine addTransition

subroutine addTransition2(thisIon,term1, term2, lambda, a, t, gamma, nt)
  type(IONTYPE) :: thisIon
  character(len=*) :: term1, term2
  real :: lambda
  real :: a
  real :: gamma(:)
  real :: thisLam
  integer :: i, j, k, m
  real :: t(:)
  integer :: nt
  character(len=80) :: message

  k = thisIon%nTransitions + 1
  thisIon%nTransitions = k

  j = returnLevel(thisIon,term1)
  i = returnLevel(thisIon,term2)

  thisIon%transition(k)%i = i
  thisIon%transition(k)%j = j
  thisLam = real(1.e8*cSpeed / ((thisIon%level(j)%energy - thisIon%level(i)%energy)/ergToEv/hCgs))

  thisIon%transition(k)%energy = (thisIon%level(j)%energy - thisIon%level(i)%energy)
  if (thisIon%transition(k)%energy < 0.) then
     write(*,*) "negative energy for ",thisIon%species,term1,term2
     stop
  endif


  if (abs(thisLam-lambda)/lambda > 0.05) then
     call writeWarning("Given wavelength and  calculated wavelength differ by more than 5%")
     write(message,*) "Calculated: ",thisLam, "Given: ", lambda
     call writeWarning(message)
     write(message,*) "For transition: ",trim(term1)," ",trim(term2), " in ion ",thisIon%species
     call writeWarning(message)
  endif
  thisIon%transition(k)%lambda = lambda
  thisIon%transition(k)%a = a


  
  thisIon%transition(k)%nGamma = 0
  do m = 1, nt
     if (gamma(m) /= 0.) then
        thisIon%transition(k)%nGamma =  thisIon%transition(k)%nGamma + 1
        thisIon%transition(k)%gamma(thisIon%transition(k)%nGamma) = gamma(m)
        thisIon%transition(k)%t(thisIon%transition(k)%nGamma) = t(m)
     endif
  enddo
end subroutine addTransition2

        


function returnLevel(ion, term) result(level)
  integer :: level
  type(IONTYPE) :: ion
  character(len=*) :: term
  integer :: i

  level = 0
  loop: do i = 1, ion%nLevels
     if (trim(ion%level(i)%term) == trim(term)) then
        level = i
        exit loop
     endif
  enddo loop
  if (level == 0) then
     write(*,*) "no term: ",trim(term)," found in ion: ",trim(ion%species)
     stop
  endif
end function returnLevel

function returnAbundance(z) result(a)

  integer :: z
  real:: a

  select case(z)
     case(1)
        a = h_abund
     case(2)
        a =  he_abund
!        write(*,*) "NO HELIUM!!!!!!!!!!!!!!"
!        a = 1.e-10
     case(6)
!        write(*,*) "NO HEAVIES!!!!!!!!!!!!!!"
        a = c_abund
     case(7)
        a = n_abund
     case(8)
        a = o_abund
     case(10)
        a = ne_abund
     case(16)
        a = s_abund
     case DEFAULT
        write(*,*) "No abundance set for z=",z
        a = tiny(a)
  end select

end function returnAbundance


function returnMu(thisOctal, subcell, ionArray, nion) result (mu)

  real(double) :: mu, tot, ne, mA, mE
  integer :: subcell
  type(OCTAL) :: thisOctal
  type(IONTYPE) :: ionArray(:)
  integer :: nion, i

  tot = 0.d0
  ne = 0.d0
  mu = 0.d0
  mA = 0.d0

  do i = 1, nIon
     !sum number of ions
     tot = tot + ionArray(i)%abundance * thisOctal%nh(subcell) * &
          thisOctal%ionFrac(subcell, i)

     !sum ion masses
     mA = mA + ionArray(i)%abundance * thisOctal%nh(subcell) * &
          dble(ionArray(i)%nucleons)*thisOctal%ionFrac(subcell, i)
  enddo


  !get number of free electrons
  !do i = 1, nIon
  !   ne = ne + ionArray(i)%abundance * thisOctal%nh(subcell) * &
  !        thisOctal%ionFrac(subcell, i) * dble(ionArray(i)%z-ionArray(i)%n)
  !enddo
  mE = thisOctal%ne(subcell)*(mElectron/mHydrogen)

  tot = tot + thisOctal%ne(subcell)
  mu = (mA + mE) / tot
  !print *, "mu = ", mu


!  real(double) :: mu, tot
!  integer :: subcell
!  type(OCTAL) :: thisOctal
!  type(IONTYPE) :: ionArray(:)
!  integer :: nion, i
 
!  tot = 0.d0
!  do i = 1, nIon
!     tot = tot + ionArray(i)%abundance * thisOctal%nh(subcell) * &
!           dble(ionArray(i)%nucleons)*thisOctal%ionFrac(subcell, i)
! enddo
!  tot = tot + thisOctal%ne(subcell)*(mElectron/mHydrogen)
!  mu = tot

end function returnMu


function returnNe(thisOctal, subcell, ionArray, nion) result (ne)
  real(double) :: ne, tot
  integer :: subcell
  type(OCTAL) :: thisOctal
  type(IONTYPE) :: ionArray(:)
  integer :: nion, i

  tot = 0.d0 
  do i = 1, nIon
     tot = tot + ionArray(i)%abundance * thisOctal%nh(subcell) * &
          thisOctal%ionFrac(subcell, i) * dble(ionArray(i)%z-ionArray(i)%n)
  enddo
  ne = tot
end function returnNe

end module ion_mod

#endif
