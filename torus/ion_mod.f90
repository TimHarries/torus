module ion_mod
!
!

! written by tjh


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
    character(len=30) :: message


    nIon = 1
    call createIon(ionArray(nIon), 1, 1, 1.360e1) ! H I

    nIon = nIon + 1
    call createIon(ionArray(nIon), 1, 0, 1.e-10) ! H II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 2, 2, 2.459e1) ! He I 

    nIon = nIon + 1
    call createIon(ionArray(nIon), 2, 1, 5.442e1) ! He II

    nIon = nIon + 1
    call createIon(ionArray(nIon), 2, 0, 1.e-10) ! He III (alpha particle!)


    write(message,*) "Metals switched OFF"
    call writeWarning(message)

!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 6, 6, 1.126e1) ! C I
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 6, 5, 2.438e1) ! C II
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 6, 4, 4.789e1) ! C III
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 6, 3, 6.449e1) ! C IV
!
!!    nIon = nIon + 1
!!    call createIon(ionArray(nIon), 6, 2, 3.921e2) ! C V
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 7, 7, 1.453e1) ! N I
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 7, 6, 2.960e1) ! N II
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 7, 5, 4.745e1) ! N III
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 7, 4, 7.747e1) ! N IV
!
!!    nIon = nIon + 1
!!    call createIon(ionArray(nIon), 7, 3, 9.789e1) ! N V
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 8, 8, 1.362e1) ! O I
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 8, 7, 3.512e1) ! O II
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 8, 6, 5.494e1) ! O III
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 8, 5, 7.741e1) ! O IV
!
!!    nIon = nIon + 1
!!    call createIon(ionArray(nIon), 8, 4, 1.139e2) ! O V
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 10, 10, 2.156e1) ! Ne I
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 10, 9, 4.096e1) ! Ne II
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 10, 8, 6.346e1) ! Ne III
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 10, 7, 9.712e1) ! Ne IV
!
!!    nIon = nIon + 1
!!    call createIon(ionArray(nIon), 10, 6, 1.262e2) ! Ne V
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 16, 16, 1.036e1) ! S I
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 16, 15, 2.333e1) ! S II
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 16, 14, 3.483e1) ! S III
!
!    nIon = nIon + 1
!    call createIon(ionArray(nIon), 16, 13, 4.731e1) ! S IV
!
!
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
!        call addLevel(thisIon, "5S^0_2", 7.274391d0, 2.)

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
     case("N I")
        allocate(t(15), gamma(15))

        t = (/ 5000, 6000, 7000, 8000, 9000, 10000, 12000, 13000, 14000, 15000, 16000, 17000, &
             18000, 19000, 20000/)
        gamma =  (/  0.1402 ,   .1732  ,  .2044  ,  .2338  ,  .2615  ,  .2877  ,  .3357  ,  .3576  ,  &
             .3784  ,  .3979  ,  .4164  ,  .4338  ,  .4503  ,  .4658  ,  .4804 /)
        call addTransition2(thisIon,"2D^0_5/2","4S^0_3/2", 5200.4, 6.13e-6, t, gamma, 15)
        gamma =  (/  0.0935 ,   .1155  ,  .1363  ,  .1559  ,  .1743  ,  .1918  ,  .2238  ,  .2384  ,  &
             .2523  ,  .2653  ,  .2776  ,  .2892  ,  .3002  ,  .3105  ,  .3203 /)
        call addTransition2(thisIon,"2D^0_3/2","4S^0_3/2", 5197.9, 2.28e-5, t, gamma, 15)
        gamma =  (/  0.1281 ,   .1591  ,  .1886  ,  .2167  ,  .2434  ,  .2609  ,  .3162  ,  .3382  ,  &
             .3591  ,  .3790  ,  .3979  ,  .4159  ,  .4331  ,  .4494  ,  .4650 /)
        call addTransition2(thisIon,"2D^0_3/2","2D^0_5/2", 1.1484e7, 1.24e-8, t, gamma, 15)
        gamma =  (/  0.0277 ,   .0342  ,  .0404  ,  .0463  ,  .0519  ,  .0573  ,  .0673  ,  .0720  ,  &
             .0765  ,  .0807  ,  .0848  ,  .0887  ,  .0924  ,  .0959  ,  .0993 /)
        call addTransition2(thisIon,"2P^0_1/2","4S^0_3/2", 3466.5, 2.71e-3, t, gamma, 15)
        gamma =  (/  0.0554 ,   .0683  ,  .0807  ,  .0925  ,  .1038  ,  .1146  ,  .1346  ,  .1440  ,  &
             0.1529  ,  .1615  ,  .1696  ,  .1774  ,  .1848  ,  .1919  ,  .1987 /)
        call addTransition2(thisIon,"2P^0_3/2","4S^0_3/2", 3466.5, 6.60e-3, t, gamma, 15)
        gamma =  (/  0.0626 ,   .0724  ,  .0819  ,  .0912  ,  .1002  ,  .1090  ,  .1261  , 0.1345  , &
             0.1427  , 0.1508  , 0.1588  , 0.1666  ,  .1744  , 0.1822  , 0.1898 /)
        call addTransition2(thisIon,"2P^0_1/2","2D^0_5/2", 10507.1, 3.14e-2, t, gamma, 15)
        gamma =  (/  0.1615 ,  0.1841  ,  .2058  ,  .2265  ,  .2466  ,  .2660  ,  .3033  ,  .3213  ,  &
             .3389  ,  .3562  ,  .3731  ,  .3898  ,  .4061  ,  .4223  ,  .4382 /)
        call addTransition2(thisIon,"2P^0_3/2","2D^0_5/2", 10397.7, 5.59e-2, t, gamma, 15)
        gamma =  (/  0.0601 ,   .0682  ,  .0758  ,  .0832  ,  .0902  ,  .0970  ,  .1100  ,  .1162  ,  &
             .1223  ,  .1283  ,  .1342  ,  .1399  ,  .1455  ,  .1510  ,  .1565 /)
        call addTransition2(thisIon,"2P^0_1/2","2D^0_3/2", 10407.6, 4.80e-2, t, gamma, 15)
        gamma =  (/  0.0856 ,   .0987  ,  .1113  ,  .1235  ,  .1354  ,  .1470  ,  .1695  ,  .1804  ,  &
             .1911  ,  .2017  ,  .2121  ,  .2224  ,  .2325  ,  .2425  ,  .2524 /)
        call addTransition2(thisIon,"2P^0_3/2","2D^0_3/2", 10407.2, 2.52e-2, t, gamma, 15)
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
        call addTransition2(thisIon,"5S^0_2","3P_1", 2144., 4.80e1, t, gamma, 7)
        gamma = (/1.6929, 1.6782, 1.6707, 1.6663, 1.6637, 1.6621, 1.6612 /)
        call addTransition2(thisIon,"1D_2","3P_2", 6583.4, 2.99e-3, t, gamma, 7)
        gamma = (/ .2118,  .2074,  .2052,  .2039,  .2031,  .2026,  .2023 /)
        call addTransition2(thisIon,"1S_0","3P_2", 3071.4, 1.51e-4, t, gamma, 7)
        gamma = (/ .6412,  .6366,  .6369,  .6383,  .6399,  .6413,  .6424 /)
        call addTransition2(thisIon,"5S^0_2","3P_2", 2140., 1.07e2, t, gamma, 7)
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
     case("O I")
        allocate(t(5), gamma(5))
        t = (/ 500,  1000, 5000, 10000, 20000 /)
        gamma =  (/  5.69E-3   ,     9.72E-3   ,    4.52E-2   ,    0.106    ,     2.066-1 /)
        call addTransition2(thisIon,"3P_1","3P_2", 6.32e5, 8.92e-5, t, gamma, 5)
        gamma =  (/  2.51E-3   ,     4.11E-3   ,    1.53E-2   ,    3.21E-2  ,     5.36E-2 /)
        call addTransition2(thisIon,"3P_0","3P_2", 4.41e5, 1.e-10, t, gamma, 5)
        gamma =  (/  0.3222E-02,     0.8389E-02,    0.6889E-01,    0.1478E+00,    0.2783E+00  /)
        call addTransition2(thisIon,"1D_2","3P_2", 6300.3, 6.34e-3, t, gamma, 5)
        gamma =  (/  0.3611E-03,     0.1022E-02,    0.8500E-02,    0.1800E-01,    0.3372E-01  /)
        call addTransition2(thisIon,"1S_0","3P_2", 2959.2, 2.88e-4, t, gamma, 5)
        gamma =  (/  2.86E-4   ,     6.36E-4   ,    9.12E-3   ,    2.83E-2  ,     6.93E-2     /)
        call addTransition2(thisIon,"3P_0","3P_1", 1.46e6, 1.74e-5, t, gamma, 5)
        gamma =  (/  0.1933E-02,     0.5033E-02,    0.4133E-01,    0.8867E-01,    0.1670E+00  /)
        call addTransition2(thisIon,"1D_2","3P_1", 6363.9, 2.11e-3, t, gamma, 5)
        gamma =  (/  0.2167E-03,     0.6133E-03,    0.5100E-02,    0.1080E-01,    0.2023E-01  /)
        call addTransition2(thisIon,"1S_0","3P_1", 2972.3, 7.32e-2, t, gamma, 5)
        gamma =  (/  0.6444E-03,     0.1678E-02,    0.1378E-01,    0.2956E-01,    0.5567E-01  /)
        call addTransition2(thisIon,"1D_2","3P_0", 6393.5, 7.23e-7, t, gamma, 5)
        gamma =  (/  0.7222E-04,     0.2044E-03,    0.1700E-02,    0.3600E-02,    0.6744E-02  /)
        call addTransition2(thisIon,"1S_0","3P_0", 2979. , 1.e-20, t, gamma, 5)
        gamma =  (/  2.1E-2    ,     3.1E-2    ,    7.32E-2   ,   1.05E-1    , 1.48E-1      /)
        call addTransition2(thisIon,"1S_0","1D_2", 5577.3, 1.22e0, t, gamma, 5)
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
        gamma = (/  .4080  ,   .4022  ,   .4037  ,   .4049 ,    .4063 ,   .4077 ,   .4106 ,  .4135 ,  .4161 ,   .4190,    .4218 /)
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
        t = (/5000, 6000, 7000, 8000, 9000,10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000 /)
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

     case("S III")
        allocate(t(4), gamma(4))
        t(1:4) = (/5000., 10000., 15000., 20000./)
        gamma = (/  2.64  , 2.59  , 2.38  ,  2.20  /)
        call addTransition2(thisIon,"3P_1", "3P_0", 3.347e5, 4.72e-4, t, gamma, 4)
        gamma = (/  1.11  , 1.15  , 1.15  ,  1.14  /)
        call addTransition2(thisIon,"3P_2", "3P_0", 1.20e5, 4.61e-8, t, gamma, 4)
        gamma = (/  1.01  , 0.932 , 0.921 ,  0.911 /)
        call addTransition2(thisIon,"1D_2", "3P_0", 8833.9, 5.82e-6, t, gamma, 4)
        gamma = (/  0.129 , 0.132 , 0.134 ,  0.138 /)
        call addTransition2(thisIon,"1S_0", "3P_0", 3681.8, 1.e-20, t, gamma, 4)
        gamma = (/  5.79  , 5.81  , 5.56  ,  5.32  /)
        call addTransition2(thisIon,"3P_2", "3P_1", 187129., 2.07e-3, t, gamma, 4)
        gamma = (/  3.02  , 2.80  , 2.76  ,  2.73  /)
        call addTransition2(thisIon,"1D_2", "3P_1", 9068.9, 2.21e-2, t, gamma, 4)
        gamma = (/  0.387 , 0.397 , 0.403 ,  0.413 /)
        call addTransition2(thisIon,"1S_0", "3P_1", 3721.7, 7.96e-1, t, gamma, 4)
        gamma = (/  5.04  , 4.66  , 4.61  ,  4.56  /)
        call addTransition2(thisIon,"1D_2", "3P_2", 9531.0, 5.76e-2, t, gamma, 4)
        gamma = (/  0.644 , 0.661 , 0.672 ,  0.689 /)
        call addTransition2(thisIon,"1S_0", "3P_2", 3797.8, 1.05e-2, t, gamma, 4)
        gamma = (/  1.42  , 1.88  , 2.02  ,  2.08  /)
        call addTransition2(thisIon,"1S_0", "1D_2", 6312.1, 2.22e0, t, gamma, 4)
        deallocate(t, gamma)

     case("S IV")
        allocate(t(4), gamma(4))
        t(1:4) = (/5000., 10000., 15000., 20000./)
        gamma = (/  6.42, 6.42  , 6.41  , 6.40 /)
        call addTransition2(thisIon,"2P^0_3/2", "2P^0_1/2", 105000., 7.73e-3, t, gamma, 4)
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
  if (thisIon%transition(k)%energy < 0.) then
     write(*,*) "negative energy for ",thisIon%species,term1,term2
     stop
  endif


  if (abs(thisLam-lambda)/lambda > 0.05) then
     write(*,*) "WARNING! Given wavelength and  calculated wavelength differ by more than 5%"
     write(*,*) "Calculated: ",thisLam, "Given: ", lambda
     write(*,*) "For transition: ",trim(term1)," ",trim(term2), " in ion ",thisIon%species
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


  k = thisIon%nTransitions + 1
  thisIon%nTransitions = k

  j = returnLevel(thisIon,term1)
  i = returnLevel(thisIon,term2)

  thisIon%transition(k)%i = i
  thisIon%transition(k)%j = j
  thisLam = 1.e8*cSpeed / ((thisIon%level(j)%energy - thisIon%level(i)%energy)/ergToEv/hCgs)

  thisIon%transition(k)%energy = (thisIon%level(j)%energy - thisIon%level(i)%energy)
  if (thisIon%transition(k)%energy < 0.) then
     write(*,*) "negative energy for ",thisIon%species,term1,term2
     stop
  endif


  if (abs(thisLam-lambda)/lambda > 0.05) then
     write(*,*) "WARNING! Given wavelength and  calculated wavelength differ by more than 5%"
     write(*,*) "Calculated: ",thisLam, "Given: ", lambda
     write(*,*) "For transition: ",trim(term1)," ",trim(term2), " in ion ",thisIon%species
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

  use input_variables, only : h_abund, he_abund, c_abund, n_abund, &
       o_abund, ne_abund, s_abund
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


end module ion_mod

