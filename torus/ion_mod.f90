module ion_mod
!
!

! written by tjh


  use kind_mod
  use constants_mod
  use utils_mod


  implicit none

  type LEVELTYPE
     character(len=10) :: term
     integer :: g
     real(double) :: energy
  end type LEVELTYPE

  type TRANSITIONTYPE
     integer :: i  ! lower level
     integer :: j  ! upper level
     real(double) :: lambda ! angstrom
     real :: energy
     integer :: nGamma
     real :: gamma(4)  ! collisional gamma at T = 5000, 10000, 15000 and 20000K
     real :: t(4) ! temperatures at which gammas are defined
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
     character(len=8) :: species  ! name of ion
     integer :: nTransitions ! number of bb transitions for this species
     integer :: nLevels ! number of energy levels
     type(LEVELTYPE) :: level(10) 
     type(TRANSITIONTYPE) :: transition(20) 
  end type IONTYPE

contains

  subroutine createIon(thisIon, z, n, ipot)
    type(IONTYPE) :: thisIon
    integer :: z, n, i
    real :: iPot
    integer :: outerShell
    character(len=2) :: element
    character(len=4) :: roman

    thisIon%z = z
    thisIon%n = n
    i = z - n + 1
    call returnElement(z, element)
    call createRoman(i, roman)
    thisIon%species = trim(element)//" "//trim(roman)
    thisIon%abundance = returnAbundance(z)
    thisIon%iPot = iPot
    thisIon%outerShell = returnOuterShell(n) 
    thisIon%nuThresh = (thisIon%iPot / dble(ergtoev)) / dble(hCgs)
    thisIon%nLevels = 0
    thisIon%nTransitions = 0
    call addLevels(thisIon)
    call addTransitions(thisIon)

  end subroutine createIon


  subroutine addIons(ionArray, nIon)
    integer :: nIon
    type(IONTYPE) :: ionArray(:)


    nIon = 1
    call createIon(ionArray(nIon), 1, 1, 1.360e1) ! H I

    nIon = nIon + 1
    call createIon(ionArray(nIon), 1, 0, 999.e9) ! H II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 2, 2, 2.459e1) ! He I 

    nIon = nIon + 1
    call createIon(ionArray(nIon), 2, 1, 5.442e1) ! He II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 6, 6, 1.126e1) ! C I

    nIon = nIon + 1
    call createIon(ionArray(nIon), 6, 5, 2.438e1) ! C II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 6, 4, 4.789e1) ! C III

    nIon = nIon + 1
    call createIon(ionArray(nIon), 6, 3, 6.449e1) ! C IV

    nIon = nIon + 1
    call createIon(ionArray(nIon), 6, 3, 3.921e2) ! C V

    nIon = nIon + 1
    call createIon(ionArray(nIon), 7, 7, 1.453e1) ! N I

    nIon = nIon + 1
    call createIon(ionArray(nIon), 7, 6, 2.960e1) ! N II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 7, 5, 4.745e1) ! N III

    nIon = nIon + 1
    call createIon(ionArray(nIon), 7, 4, 7.747e1) ! N IV

    nIon = nIon + 1
    call createIon(ionArray(nIon), 7, 3, 9.789e1) ! N V

    nIon = nIon + 1
    call createIon(ionArray(nIon), 8, 8, 1.362e1) ! O I

    nIon = nIon + 1
    call createIon(ionArray(nIon), 8, 7, 3.512e1) ! O II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 8, 6, 5.494e1) ! O III

    nIon = nIon + 1
    call createIon(ionArray(nIon), 8, 5, 7.741e1) ! O IV

    nIon = nIon + 1
    call createIon(ionArray(nIon), 8, 4, 1.139e2) ! O V

    nIon = nIon + 1
    call createIon(ionArray(nIon), 10, 10, 2.156e1) ! Ne I

    nIon = nIon + 1
    call createIon(ionArray(nIon), 10, 9, 4.096e1) ! Ne II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 10, 8, 6.346e1) ! Ne III

    nIon = nIon + 1
    call createIon(ionArray(nIon), 10, 7, 9.712e1) ! Ne IV

    nIon = nIon + 1
    call createIon(ionArray(nIon), 10, 6, 1.262e2) ! Ne V

    nIon = nIon + 1
    call createIon(ionArray(nIon), 16, 16, 1.036e1) ! S I

    nIon = nIon + 1
    call createIon(ionArray(nIon), 16, 15, 2.333e1) ! S II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 16, 14, 3.483e1) ! S III

    nIon = nIon + 1
    call createIon(ionArray(nIon), 16, 13, 4.731e1) ! S IV


    write(*,*) "Added ",nion," species to photoionization calculation"

  end subroutine addIons
    
function returnOuterShell(n) result (ishell)
  integer :: z, n, ishell

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
        
!     case("Ne II")
!        call addLevel(thisIon, "2P^0_3/2", 0.d0, 1.5)
!        call addLevel(thisIon, "2P^0_1/2", 0.096750d0, 1.5)
!        call addLevel(thisIon, "2S_1/2", 26.91048d0, 1.5)
!        call addLevel(thisIon, "4P_5/2", 27.16876d0, 2.5)
!        call addLevel(thisIon, "4P_3/2", 27.2394d0, 1.5)
!        call addLevel(thisIon, "4P_1/2", 27.2700d0, 0.5)

!     case("Ne III")
!        call addLevel(thisIon, "3P_2", 0.d0, 2.)
!        call addLevel(thisIon, "3P_1", 0.07971d0, 1.)
!        call addLevel(thisIon, "3P_0", 0.11413d0, 0.)
!        call addLevel(thisIon, "1D_2", 3.20385d0, 2.)
!        call addLevel(thisIon, "1S_0", 6.91220d0, 0.)

     case("S II")
        call addLevel(thisIon, "4S^0_3/2", 0.d0, 1.5)
        call addLevel(thisIon, "2D^0_3/2", 1.841530d0, 1.5)
        call addLevel(thisIon, "2D^0_5/2", 1.845471d0, 2.5)
        call addLevel(thisIon, "2P^0_1/2", 3.040691d0, 0.5)
        call addLevel(thisIon, "2P^0_3/2", 3.046482d0, 1.5)
!        call addLevel(thisIon, "4P_5/2", 9.83774d0, 2.5)
!        call addLevel(thisIon, "4P_3/2", 9.880586d0, 1.5)
!        call addLevel(thisIon, "4P_1/2", 9.914099d0, 0.5)

     case("S III")
        call addLevel(thisIon, "3P_0", 0.d0, 0.)
        call addLevel(thisIon, "3P_1", 0.037033d0, 1.)
        call addLevel(thisIon, "3P_2", 0.10329d0, 2.)
        call addLevel(thisIon, "1D_2", 1.403841d0, 2.)
        call addLevel(thisIon, "1S_0", 3.36753d0, 0.)
        call addLevel(thisIon, "5S^0_2", 7.274391d0, 2.)
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
  thisIon%level(i)%g = int(real(2.*J+1.))
end subroutine addLevel

subroutine addTransitions(thisIon)
  type(IONTYPE) :: thisIon
  
  select case(thisIon%species)
     case("N I")
        call addTransition(thisIon,"2D^0_5/2","4S^0_3/2", 5200.4, 6.13e-6, 1.55e-1, 2.90e-1, 0., 4.76e-1)
        call addTransition(thisIon,"2D^0_3/2","4S^0_3/2", 5197.9, 2.28e-5, 1.30e-1, 1.93e-1, 0., 3.18e-1)
        call addTransition(thisIon,"2P^0_3/2","4S^0_3/2", 3466.5, 6.60e-3, 5.97e-2, 1.13e-1, 0., 1.89e-1)
        call addTransition(thisIon,"2D^0_3/2","2D^0_5/2", 1.1484e7, 1.24e-8, 1.28e-1, 2.69e-1, 0., 4.65e-1)
        call addTransition(thisIon,"2P^0_3/2","2P^0_1/2", 2.59e8, 5.17e-13, 3.29e-2, 7.10e-2, 0., 1.53e-1)
        call addTransition(thisIon,"2P^0_3/2","2D^0_5/2", 10397.7, 5.59e-2, 1.62e-1, 2.66e-1, 0., 4.38e-1)
        call addTransition(thisIon,"2P^0_3/2","2D^0_3/2", 10407.2, 2.52e-2, 8.56e-2, 1.47e-1, 0., 2.52e-1)
        call addTransition(thisIon,"2P^0_1/2","2D^0_5/2", 10507.1, 3.14e-2, 6.26e-2, 1.09e-1, 0., 1.90e-1)
        call addTransition(thisIon,"2P^0_1/2","2D^0_3/2", 10407.6, 4.80e-2, 6.01e-2, 9.70e-2, 0., 1.57e-1)
     case("N II")
        call addTransition(thisIon,"1D_2","3P_0", 6529.0, 5.35e-7, 2.57e0, 2.64e0, 2.70e0, 2.73e0)
        call addTransition(thisIon,"1D_2","3P_1", 6548.1, 1.01e-3, 2.57e0, 2.64e0, 2.70e0, 2.73e0)
        call addTransition(thisIon,"1D_2","3P_2", 6583.4, 2.99e-3, 2.57e0, 2.64e0, 2.70e0, 2.73e0)
        call addTransition(thisIon,"1S_0","3P_1", 3062.9, 3.38e-2, 2.87e-1, 2.93e-1, 3.00e-1, 3.05e-1)
        call addTransition(thisIon,"1S_0","3P_2", 3071.4, 1.51e-4, 2.87e-1, 2.93e-1, 3.00e-1, 3.05e-1)
        call addTransition(thisIon,"1S_0","1D_2", 5754.6, 1.12e0, 9.59e-1, 8.34e-1, 7.61e-1, 7.34e-1)
        call addTransition(thisIon,"3P_1","3P_0", 2.055e6, 2.08e-6, 3.71e-1, 4.08e-1, 4.29e-1, 4.43e-1)
        call addTransition(thisIon,"3P_2","3P_0", 7.65e5, 1.16e-12, 2.43e-1, 2.72e-1, 3.01e-1, 3.16e-1)
        call addTransition(thisIon,"3P_2","3P_1", 1.22e6, 7.46e-6, 1.01e0, 1.12e0, 1.21e0, 1.26e0)
        call addTransition(thisIon,"5S^0_2","3P_1", 2144., 4.80e1, 1.19e0, 1.19e0, 1.21e0, 1.21e0)
        call addTransition(thisIon,"5S^0_2","3P_2", 2140., 1.07e2, 1.19e0, 1.19e0, 1.21e0, 1.21e0)
     case("N III")
        call addTransition(thisIon,"2P^0_3/2","2P^0_1/2", 5.73e5, 4.77e-5, 1.32e0, 1.45e0, 1.55e0, 1.64e0)
        call addTransition(thisIon,"4P_1/2","2P^0_1/2", 1748., 3.39e2, 1.89e-1, 1.98e-1, 2.04e-1, 2.07e-1)
        call addTransition(thisIon,"4P_1/2","2P^0_3/2", 1754., 3.64e2, 1.35e-1, 1.51e-1, 1.62e-1, 1.68e-1)
        call addTransition(thisIon,"4P_3/2","2P^0_1/2", 1747., 8.95e2, 2.81e-1, 2.98e-1, 3.09e-1, 3.16e-1)
        call addTransition(thisIon,"4P_3/2","2P^0_3/2", 1752., 5.90e1, 3.67e-1, 3.99e-1, 4.23e-1, 4.35e-1)
        call addTransition(thisIon,"4P_3/2","4P_1/2", 1.68e6, 1.e-20, 1.01e0, 1.10e0, 1.14e0, 1.16e0)
        call addTransition(thisIon,"4P_5/2","2P^0_1/2", 1744., 1.e-20, 1.78e-1, 2.01e-1, 2.19e-1, 2.29e-1)
        call addTransition(thisIon,"4P_5/2","2P^0_3/2", 1747., 3.08e2, 7.93e-1, 8.44e-1, 8.80e-1, 8.98e-1)
        call addTransition(thisIon,"4P_5/2","4P_1/2", 7.10e5, 1.e-20, 6.12e-1, 6.67e-1, 6.95e-1, 7.11e-1)
        call addTransition(thisIon,"4P_5/2","4P_3/2", 1.23e6, 1.e-20, 1.88e0, 2.04e0, 2.12e0, 2.16e0)
     case("O I")
        call addTransition(thisIon,"1D_2","3P_0", 6393.5, 7.23e-7, 1.2e-1, 2.66e-1, 0., 5.01e-1)
        call addTransition(thisIon,"1D_2","3P_1", 6363.9, 2.11e-3, 1.2e-1, 2.66e-1, 0., 5.01e-1)
        call addTransition(thisIon,"1D_2","3P_2", 6300.3, 6.34e-3, 1.2e-1, 2.66e-1, 0., 5.01e-1)
        call addTransition(thisIon,"1S_0","3P_1", 2972.3, 7.32e-2, 1.53e-2, 3.24e-2, 0., 6.07e-2)
        call addTransition(thisIon,"1S_0","3P_2", 2959.2, 2.88e-4, 1.53e-2, 3.24e-2, 0., 6.07e-2)
        call addTransition(thisIon,"1S_0","1D_2", 5577.3, 1.22e0, 7.32e-2, 1.05e-1, 0., 1.48e-1)
        call addTransition(thisIon,"3P_0","3P_1", 1.46e6, 1.74e-5, 1.12e-2, 2.65e-2, 0., 6.93e-2)
        call addTransition(thisIon,"3P_0","3P_2", 4.41e5, 1.e-10, 1.48e-2, 2.92e-2, 0., 5.36e-2)
        call addTransition(thisIon,"3P_1","3P_2", 6.32e5, 8.92e-5, 4.74e-2, 9.87e-2, 0., 2.07e-1)
     case("O II")
        call addTransition(thisIon,"2D^0_5/2","4S^0_3/2", 3728.8, 3.50e-5, 7.95e-1, 8.01e-1, 8.10e-1, 8.18e-1)
        call addTransition(thisIon,"2D^0_3/2","4S^0_3/2", 3726.0, 1.79e-4, 5.30e-1, 5.34e-1, 5.41e-1, 5.45e-1)
        call addTransition(thisIon,"2P^0_3/2", "4S^0_3/2", 2470.3, 5.70e-2, 2.65e-1, 2.70e-1, 2.75e-1, 2.80e-1)
        call addTransition(thisIon,"2P^0_1/2", "4S^0_3/2", 2470.2, 2.34e-2, 1.33e-1, 1.35e-1, 1.37e-1, 1.40e-1)
        call addTransition(thisIon,"2D^0_3/2", "2D^0_5/2", 4.97e6, 1.30e-7, 1.22e0, 1.17e0, 1.14e0, 1.11e0)
        call addTransition(thisIon,"2P^0_1/2", "2P^0_3/2", 5.00e7, 2.08e-11, 2.80e-1, 2.87e-1, 2.93e-1, 3.00e-1)
        call addTransition(thisIon,"2P^0_3/2", "2D^0_5/2", 7319.9, 1.07e-1, 7.18e-1, 7.30e-1, 7.41e-1, 7.55e-1)
        call addTransition(thisIon,"2P^0_3/2", "2D^0_3/2", 7330.7, 5.78e-2, 4.01e-1, 4.08e-1, 4.14e-1, 4.22e-1)
        call addTransition(thisIon,"2P^0_1/2", "2D^0_5/2", 7321.8, 6.15e-2, 2.90e-1, 2.95e-1, 3.00e-1, 3.05e-1)
        call addTransition(thisIon,"2P^0_1/2", "2D^0_3/2", 7329.6, 1.02e-1, 2.70e-1, 2.75e-1, 2.80e-1, 2.84e-1)
     case("O III")
        call addTransition(thisIon,"1D_2","3P_0", 4932.6, 2.74e-6, 2.13e0, 2.29e0, 2.45e0, 2.52e0)
        call addTransition(thisIon,"1D_2","3P_1", 4958.9, 6.74e-3, 2.13e0, 2.29e0, 2.45e0, 2.52e0)
        call addTransition(thisIon,"1D_2","3P_2", 5006.7, 1.96e-2, 2.13e0, 2.29e0, 2.45e0, 2.52e0)
        call addTransition(thisIon,"1S_0","3P_1", 2321.0, 2.23e-1, 2.72e-1, 2.93e-1, 3.17e-1, 3.29e-1)
        call addTransition(thisIon,"1S_0","3P_2", 2332.0, 7.85e-4, 2.72e-1, 2.93e-1, 3.17e-1, 3.29e-1)
        call addTransition(thisIon,"1S_0","1D_2", 4363.2, 1.78e0, 4.94e-1, 5.82e-1, 6.10e-1, 6.10e-1)
        call addTransition(thisIon,"3P_1","3P_0", 883562., 2.62e-5, 5.24e-1, 5.45e-1, 5.59e-1, 5.63e-1)
        call addTransition(thisIon,"3P_2","3P_0", 326611., 3.02e-11, 2.58e-1, 2.71e-1, 2.83e-1, 2.89e-1)
        call addTransition(thisIon,"3P_2","3P_1", 518145., 9.76-5, 1.23e0, 1.29e0, 1.34e0, 1.35e0)
        call addTransition(thisIon,"5S^0_2","3P_1", 1660.8, 2.12e2, 1.07e0, 1.21e0, 1.25e0, 1.26e0)
        call addTransition(thisIon,"5S^0_2","3P_2", 1666.1, 5.22e2, 1.07e0, 1.21e0, 1.25e0, 1.26e0)
!     case("Ne II")
!        call addTransition(thisIon,"2P^0_1/2","2P^0_3/2", 1.28e5, 8.55e-3, 2.96e-1, 3.03e-1, 3.10e-1, 3.17e-1)
!     case("Ne III")
!        call addTransition(thisIon,"1D_2","3P_0", 4012.8, 8.51e-6, 1.63e0, 1.65e0, 1.65e0, 1.64e0)
!        call addTransition(thisIon,"1D_2","3P_1", 3967.5, 5.42e-2, 1.63e0, 1.65e0, 1.65e0, 1.64e0)
!        call addTransition(thisIon,"1D_2","3P_2", 3868.8, 1.71e-1, 1.63e0, 1.65e0, 1.65e0, 1.64e0)
!        call addTransition(thisIon,"1S_0","3P_1", 1814.6, 2.00e0, 1.51e-1, 1.69e-1, 1.75e-1, 1.79e-1)
!        call addTransition(thisIon,"1S_0","3P_2", 1793.7, 3.94e-2, 1.51e-1, 1.69e-1, 1.75e-1, 1.79e-1)
!        call addTransition(thisIon,"1S_0","1D_2", 3342.5, 2.71e0, 2.00e-1, 2.26e-1, 2.43e-1, 2.60e-1)
!        call addTransition(thisIon,"3P_0", "3P_1", 3.60e5, 1.15-3, 3.31e-1, 3.50e-1, 3.51e-1, 3.50e-1)
!        call addTransition(thisIon,"3P_0", "3P_2", 1.07e5, 2.18e-8, 3.00e-1, 3.07e-1, 3.03e-1, 2.98e-1)
!        call addTransition(thisIon,"3P_1", "3P_2", 1.56e5, 5.97e-3, 1.09e0, 1.65e0, 1.65e0, 1.64e0)
     case("S II")
        call addTransition(thisIon,"2D^0_5/2", "4S^0_3/2", 6716.5, 2.60e-4, 4.90e0, 4.66e0, 4.44e0, 4.26e0)
        call addTransition(thisIon,"2D^0_3/2", "4S^0_3/2", 6730.8, 8.82e-4, 3.27e0, 3.11e0, 2.97e0, 2.84e0)
        call addTransition(thisIon,"2P^0_3/2", "4S^0_3/2", 4068.6, 2.25e-1, 1.67e0, 2.07e0, 1.98e0, 2.07e0)
        call addTransition(thisIon,"2P^0_1/2", "4S^0_3/2", 4076.4, 9.06e-2, 8.31e-1, 8.97e-1, 9.87e-1, 1.03e0)
        call addTransition(thisIon,"2D^0_5/2", "2D^0_3/2", 3.145e6, 3.35e-7, 7.90e0, 7.46e0, 7.11e0, 8.65e0)
        call addTransition(thisIon,"2P^0_3/2", "2P^0_1/2", 2.14e6, 1.03-6, 2.02e0, 2.54e0, 2.13e0, 2.22e0)
        call addTransition(thisIon,"2P^0_3/2", "2D^0_5/2", 10320.4, 1.79e-1, 5.93e0, 4.77e0, 4.75e0, 4.68e0)
        call addTransition(thisIon,"2P^0_3/2", "2D^0_3/2", 10286.7, 1.33e-1, 3.41e0, 2.74e0, 2.74e0, 2.71e0)
        call addTransition(thisIon,"2P^0_1/2", "2D^0_5/2", 10373.3, 7.79e-2, 2.47e0, 1.99e0, 1.99e0, 1.97e0)
        call addTransition(thisIon,"2P^0_1/2", "2D^0_3/2", 10336.3, 1.63e-1, 2.20e0, 1.76e0, 1.76e0, 1.73e0)
     case("S III")
        call addTransition(thisIon,"1D_2", "3P_0", 8833.9, 5.82e-6, 9.07e0, 8.39e0, 8.29e0, 8.30e0)
        call addTransition(thisIon,"1D_2", "3P_1", 9068.9, 2.21e-2, 9.07e0, 8.39e0, 8.29e0, 8.30e0)
        call addTransition(thisIon,"1D_2", "3P_2", 9531.0, 5.76e-2, 9.07e0, 8.39e0, 8.29e0, 8.30e0)
        call addTransition(thisIon,"1S_0", "3P_1", 3721.7, 7.96e-1, 1.16e0, 1.19e0, 1.21e0, 1.24e0)
        call addTransition(thisIon,"1S_0", "3P_2", 3797.8, 1.05e-2, 1.16e0, 1.19e0, 1.21e0, 1.24e0)
        call addTransition(thisIon,"1S_0", "1D_2", 6312.1, 2.22e0, 1.42e0, 1.88e0, 2.02e0, 2.08e0)
        call addTransition(thisIon,"3P_1", "3P_0", 3.347e5, 4.72e-4, 2.64e0, 2.59e0, 2.38e0, 2.20e0)
        call addTransition(thisIon,"3P_2", "3P_0", 1.20e5, 4.61e-8, 1.11e0, 1.15e0, 1.15e0, 1.14e0)
        call addTransition(thisIon,"3P_2", "3P_1", 187129., 2.07e-3, 5.79e0, 5.81e0, 5.56e0, 5.32e0)
        call addTransition(thisIon,"5S^0_2", "3P_1", 1683.5, 6.22e3, 0., 3.8e0, 3.7e0, 3.6e0)
        call addTransition(thisIon,"5S^0_2", "3P_2", 1698.86, 1.70e4, 0., 3.8e0, 3.7e0, 3.6e0)
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
  thisLam = 1.e8*cSpeed / ((thisIon%level(j)%energy - thisIon%level(i)%energy)/ergToEv/hCgs)

  thisIon%transition(k)%energy = (thisIon%level(j)%energy - thisIon%level(i)%energy)


  if (abs(thisLam-lambda)/lambda > 0.05) then
     write(*,*) "WARNING! Given wavelength and  calculated wavelength differ by more than 5%"
     write(*,*) "Calculated: ",thisLam, "Given: ", lambda
     write(*,*) "For transition: ",trim(term1)," ",trim(term2), " in ion ",thisIon%species
  endif
  thisIon%transition(k)%lambda = lambda
  thisIon%transition(k)%a = a
  thisIon%transition(k)%nGamma = 0
  do m = 1, 4
     if (gamma(m) /= 0) then
        thisIon%transition(k)%nGamma =  thisIon%transition(k)%nGamma + 1
        thisIon%transition(k)%gamma(thisIon%transition(k)%nGamma) = gamma(m)
        thisIon%transition(k)%t(thisIon%transition(k)%nGamma) = t(m)
     endif
  enddo
end subroutine addTransition

        


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

end module ion_mod

