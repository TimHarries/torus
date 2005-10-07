module ion_mod
!
!

! written by tjh


  use kind_mod
  use constants_mod
  use utils_mod


  implicit none

  type LEVELTYPE
     character(len=6) :: term
     integer :: g
     real(double) :: energy
  end type LEVELTYPE

  type TRANSITIONTYPE
     integer :: i  ! lower level
     integer :: j  ! upper level
     real(double) :: lambda ! angstrom
     real :: gamma(4)  ! collisional gamma at T = 5000, 10000, 15000 and 20000K
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
        call addLevel(thisIon, "1S^0_2", 4.05290d0, 2.)
        call addLevel(thisIon, "3D^0_3", 11.43596d0, 3.)
        call addLevel(thisIon, "3D_2", 11.43758d0, 2.)
        call addLevel(thisIon, "3D_1", 11.43777d0, 1.)
        call addLevel(thisIon, "5S^0_2", 28.49641d0, 2.)

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
        call addLevel(thisIon, "2D^0_5/2", 3.324086d0, 2.5)
        call addLevel(thisIon, "2D^0_3/2", 3.326568d0, 1.5)
        call aadLevel(thisIon, "2P^0_3/2", 5.017396d0, 1.5)
        call addLevel(thisIon, "2P^0_1/2", 5.017642d0, 0.5)

     case("O III")
        call addLevel(thisIon, "3P_0", 0.d0, 0.)
        call addLevel(thisIon, "3P_1", 0.14032d0, 1.)
        call addLevel(thisIon, "3P_2", 0.037960d0, 2.)
        call addLevel(thisIon, "1D_2", 2.513566d0, 2.)
        call addLevel(thisIon, "1S_0", 5.354351d0, 0.)
        call addLevel(thisIon, "5S^0_2", 7.479323d0, 2.)
        
     case("Ne II")
        call addLevel(thisIon, "2P^0_3/2", 0.d0, 1.5)
        call addLevel(thisIon, "2P^0_1/2", 0.096750d0, 1.5)
        call addLevel(thisIon, "2S_1/2", 26.91048d0, 1.5)
        call addLevel(thisIon, "4P_5/2", 27.16876d0, 2.5)
        call addLevel(thisIon, "4P_3/2", 27.2394d0, 1.5)
        call addLevel(thisIon, "4P_1/2", 27.2700d0, 0.5)

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
        call addLevel(thisIon, "4P_5/2", 9.83774d0, 2.5)
        call addLevel(thisIon, "4P_3/2", 9.880586d0, 1.5)
        call addLevel(thisIon, "4P_1/2", 9.914099d0, 0.5)

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

  i = thisIon%nTransitions + 1
  thisIon%nTransitions = i
  thisIon%level(i)%term = term
  thisIon%level(i)%energy = energy
  thisIon%level(i)%g = int(real(2.*J+1.))
end subroutine addLevel

subroutine addTransitions(thisIon)
  type(IONTYPE) :: thisIon
  
  select case(thisIon%species)
     case("N I")
        call addTransition("
end module ion_mod

