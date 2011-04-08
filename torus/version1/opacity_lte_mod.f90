module opacity_lte_mod

  !
  ! Contains the function to return the LET oapcities.
  ! 
  
  use atom_mod
  
  implicit none

  public:: compute_LTE_opacities
  private:: solveForbidden


contains

  !
  !
  !
  ! For given ion name (ion_kind), ion fraction (ionFrac), electron number density (ne), 
  ! mass density (rho) and temperature (temp),  this routine computes and return
  ! the LTE opacities (etal, chil, eta, chi, esec). 
  !
  ! Currently, the following species are implemented in this rouine.
  !   
  ! ion_kind: "H I"                                -- (6 levels)
  !           "[N I]",    "[N II]",                -- (5 levels)
  !           "[O I]",    "[O II]",  "[O III]"     -- (5 levels)
  !           "[Ne III]", "[Ne IV]", "[Ne V]",     -- (5 levels)
  !           "[S II]",   "[S III]"                -- (5 levels)
  !           "[Cl II]",  "[Cl III]", "[Cl IV]"    -- (5 levels)
  !           "[Ar IV]",  "[Ar V]"                 -- (5 levels)
  !
  !
  !
  subroutine compute_LTE_opacities(ion_name, ion_frac, nLower, nUpper, ne, rho, temp, &
       eta_line, chi_line, eta_cont, chi_cont, chi_e_scat)
    implicit none

    character(LEN=*), intent(in)   :: ion_name  ! name the ion Must be one of the kinds listed above.
    real, intent(in)               :: ion_frac   ! ionization fraction of this specie.
    integer, intent(in)            :: nLower, nUpper ! Levels of transions.
    real, intent(in)               :: ne             ! electron nunmber density in [g/cm^3]
    real, intent(in)               :: rho            ! mass density in [g/cm^3]
    real, intent(in)               :: temp           ! temperature in [K]
    real, intent(out)              :: eta_line       ! line emissivity   in [cgs]
    real, intent(out)              :: chi_line       ! line opacity      in [cgs]
    real, intent(out)              :: eta_cont       ! thermal emissivity   in [cgs]
    real, intent(out)              :: chi_cont       ! thermal opacity   in [cgs]
    real, intent(out)              :: chi_e_scat     ! electron scattering opacity.  in [cgs]


    select case(ion_name)

    case("H I")
       ! uses the routine in atom_mod.f90
       call getLTEopacities(ne, temp, eta_line, chi_line,  eta_cont, &
            chi_cont, chi_e_scat, nLower, nUpper)
       
    case default
       ! uses the routine in this module. 
       call solveForbidden(temp, ne, rho, ion_name, ion_frac, nLower, nUpper, eta_line)

       ! we assume that all the other values (chi_line, eta_cont,chi_cont, chi_e_scat)
       ! to be insignificantlt small.
       chi_line   = 1.0e-21
       eta_cont   = 1.0e-21
       chi_cont   = 1.0e-21
       chi_e_scat = 1.0e-21       
       
    end select
    
    
  end subroutine compute_LTE_opacities



  !
  !
  !  Computes the line emissiovities of forbidden lines for certains species
  !
  subroutine solveForbidden(temperature, electronDensity, massDensity, ion, ionFrac, nLower, nUpper, etaLine)
    implicit none

    !        write(ttyout,*)'LICK 5-LEVEL ATOM-POPULATIONS CATEGORY'
    !        write(ttyout,*)'       *** AVAILABLE IONS ***'
    !        write(ttyout,*)'   7000: [N I];       7010: [N II];'
    !        write(ttyout,*)'   8000: [O I];       8010: [O II];'
    !        write(ttyout,*)'   8020: [O III];    10020: [Ne III];'
    !        write(ttyout,*)'  10030: [Ne IV];    10040: [Ne V];  '
    !        write(ttyout,*)'  16010: [S II];     16020: [S III]; '
    !        write(ttyout,*)'  17010: [Cl II];    17020: [Cl III];'
    !        write(ttyout,*)'  17030: [Cl IV];    18020: [Ar III];'
    !        write(ttyout,*)'  18030: [Ar IV];    18040: [Ar V];  '
    !  120   write(ttyout,20)
    !        read(kbin,*,err=120,end=120) ionum


    integer :: nLower, nUpper
    real :: e(5,5), q(5,5), a1(5,5), a2(5,5), jl(5,5)
    real :: pop(5), b1(5), b2(5), dtl(5), pop1(5), ncrit(5), kt, hc
    integer :: weight(5)
    integer :: filout=6
    real :: etaLine, nIon
    character(len=*) :: ion
    real :: electronDensity, den
    real :: temperature, tem
    integer :: i, j, k, l
    real :: t4
    real :: t(5)
    real :: cj(3)
    real :: a(5,5), c(5,5)
    real :: cq
    real :: ts
    real :: t1, t2, xl(5), xi(5)
    real :: asum, qsum
    real :: logAbund
    real :: nH
    real :: massDensity
    real :: ionFrac

    den = electronDensity
    tem = temperature




    do i=1,5
       do  j=1,5
          a(i,j)=0.0
          c(i,j)=0.0
       enddo
    enddo

    t(1)=0.0
    t4=1.0e-4*tem
    cj(1)=0.0
    cj(2)=0.0
    cj(3)=0.0


    select case(ion)
       !
       !  N I ion parameters
       !
    case("[N I]")
       a(2,1)=7.27e-6
       a(3,1)=2.02e-5
       a(4,1)=2.71e-3
       a(5,1)=6.58e-3
       a(3,2)=1.27e-8
       a(4,2)=3.45e-2
       a(5,2)=6.14e-2
       a(4,3)=5.29e-2
       a(5,3)=2.76e-2
       t(2)=1.92245e-4
       t(3)=1.92332e-4
       t(4)=2.88389e-4
       t(5)=2.88393e-4
       weight(1)=4
       weight(2)=6
       weight(3)=4
       weight(4)=2
       weight(5)=4
       !
       !  N I ion parameters
       !
       c(2,1)=-1.2124e-2+t4*(0.3616-t4*5.88e-2)
       c(3,1)=-8.386e-3+t4*(0.2419-t4*3.94e-2)
       c(4,1)=-2.601e-3+t4*(0.07003-t4*1.07e-2)
       c(5,1)=-5.074e-3+t4*(0.1395-t4*0.02126)
       c(3,2)=-4.166e-2+t4*(0.368-t4*0.05733)
       c(4,2)=1.227e-2+t4*(0.1046-t4*7.868e-3)
       c(5,2)=4.600e-2+t4*(0.244-t4*0.0240)
       c(4,3)=1.8601e-2+t4*(8.76e-2-t4*0.0092)
       c(5,3)=1.8265e-2+t4*(0.1406-t4*1.187e-2)
       c(5,4)=-3.265e-3+t4*(0.0704-t4*0.00387)
       logAbund = 8.05


       !
       !  N II ion parameters
       !
    case("[N II]")
       a(2,1)=2.08e-6
       a(3,1)=1.16e-12
       a(4,1)=5.35e-7
       a(3,2)=7.46e-6
       a(4,2)=1.01e-3
       a(5,2)=3.38e-2
       a(4,3)=2.99e-3
       a(5,3)=1.51e-4
       a(5,4)=1.12   
       t(2)=4.87e-7
       t(3)=1.308e-6
       t(4)=1.53162e-4
       t(5)=3.26888e-4
       weight(1)=1
       weight(2)=3
       weight(3)=5
       weight(4)=5
       weight(5)=1
       !
       !  N II ion parameters
       !
       cj(1)=2.577+t4*(0.137-t4*0.03)
       cj(2)=0.353+t4*(-0.005807+t4*0.006003)
       c(2,1)=0.401
       c(3,1)=0.279
       c(4,1)=cj(1)/9
       c(5,1)=cj(2)/9
       c(3,2)=1.128
       c(4,2)=cj(1)*3/9
       c(5,2)=cj(2)*3/9
       c(4,3)=cj(1)*5/9
       c(5,3)=cj(2)*5/9
       c(5,4)=0.3993+t4*(0.0109+t4*0.001)
       logAbund = 8.05


       !
       !  O I ion parameters
       !
    case("[O I]")
       a(2,1)=8.92e-5
       a(4,1)=6.34e-3
       a(5,1)=2.88e-4
       a(3,2)=1.74e-5
       a(4,2)=2.11e-3
       a(5,2)=7.32e-2
       a(4,3)=7.23e-7
       a(5,4)=1.22   
       t(2)=1.583e-6
       t(3)=2.27e-6
       t(4)=1.58679e-4
       t(5)=3.37926e-4
       weight(1)=5
       weight(2)=3
       weight(3)=1
       weight(4)=5
       weight(5)=1
       !
       !  O I ion parameters
       !
       cj(1)=0.016656+t4*(0.3004-t4*0.02067)
       cj(2)=0.00201+t4*(0.03686-t4*0.002742)
       c(2,1)=-0.00214+t4*(0.09715+t4*0.003707)
       c(3,1)=-0.00109+t4*(0.0332-t4*0.00295)
       c(4,1)=cj(1)*5/9
       c(5,1)=cj(2)*5/9
       c(3,2)=-0.000128+t4*(0.0186+t4*0.00807)
       c(4,2)=cj(1)*3/9
       c(5,2)=cj(2)*3/9
       c(4,3)=cj(1)/9
       c(5,3)=cj(2)/9
       c(5,4)=0.0218+t4*(0.1076-t4*0.02234)
       logAbund = 8.93


       !
       !  O II ion parameters
       !
    case("[O II]")
       a(2,1)=3.82e-5
       a(3,1)=1.65e-4
       a(4,1)=5.64e-2
       a(5,1)=2.32e-2
       a(3,2)=1.2e-7
       a(4,2)=1.17e-1
       a(5,2)=6.15e-2
       a(4,3)=6.14e-2
       a(5,3)=1.02e-1
       a(5,4)=2.08e-11
       t(2)=2.68107e-4
       t(3)=2.68302e-4
       t(4)=4.04675e-4
       t(5)=4.04686e-4
       weight(1)=4
       weight(2)=6
       weight(3)=4
       weight(4)=4
       weight(5)=2
       !
       !  O II ion parameters
       !
       c(2,1)=0.7894+t4*(0.0098+t4*0.002282)
       c(3,1)=0.5269+t4*(5.204e-3+t4*1.923e-3)
       c(4,1)=0.26+t4*(1.001e-2-t4*2.915e-6)
       c(5,1)=0.1311+t4*(3.445e-3+t4*4.99e-4)
       c(3,2)=1.283+t4*(-0.1389+t4*0.0262)
       c(4,2)=0.7061+t4*(0.0235+t4*4.9e-4)
       c(5,2)=0.285+t4*0.01
       c(4,3)=0.3948+t4*(0.0123+t4*6.262e-4)
       c(5,3)=0.2644+t4*(0.0117-t4*9.16e-4)
       c(5,4)=0.2731+t4*(0.014-t4*2.919e-4)
       logAbund = 8.93
       !
       !  O III ion parameters
       !
    case("[O III]")
       a(2,1)=2.62e-5
       a(3,1)=3.02e-11
       a(4,1)=2.74e-6
       a(3,2)=9.76e-5
       a(4,2)=6.74e-3
       a(5,2)=2.23e-1
       a(4,3)=1.96e-2
       a(5,3)=7.85e-4
       a(5,4)=1.78
       t(2)=1.132e-6
       t(3)=3.062e-6
       t(4)=2.02733e-4
       t(5)=4.318577e-4
       weight(1)=1
       weight(2)=3
       weight(3)=5
       weight(4)=5
       weight(5)=1
       !
       !  O III ion parameters
       !
       cj(1)=1.835+t4*(0.3981-t4*0.06)
       cj(2)=0.2127+t4*(0.0767-0.013*t4)
       c(2,1)=0.4825+t4*(0.0806-t4*0.022)
       c(3,1)=0.2397+t4*(0.0381-t4*0.007)
       c(4,1)=cj(1)/9
       c(5,1)=cj(2)/9
       c(3,2)=1.1325+t4*(0.203-t4*0.05)
       c(4,2)=cj(1)/3
       c(5,2)=cj(2)/3
       c(4,3)=cj(1)*5/9
       c(5,3)=cj(2)*5/9
       c(5,4)=0.3763+t4*(0.3375-t4*0.105)
       logAbund = 8.93


       !
       !  Ne III ion parameters
       !
    case("[Ne III]")
       a(2,1)=5.97e-3
       a(3,1)=2.18e-8
       a(4,1)=1.71e-1
       a(5,1)=3.94e-3
       a(3,2)=1.15e-3
       a(4,2)=5.42e-2
       a(5,2)=2.0
       a(4,3)=8.51e-6
       a(5,4)=2.71
       t(2)=6.429e-6
       t(3)=9.205e-6
       t(4)=2.58408e-4
       t(5)=5.57506e-4
       weight(1)=5
       weight(2)=3
       weight(3)=1
       weight(4)=5
       weight(5)=1
       !
       !  Ne III ion parameters
       !
       cj(1)=1.6028+t4*(0.07967-t4*0.03282)
       cj(2)=0.12956+t4*(0.05331-t4*0.01487)
       c(2,1)=1.0314+t4*(0.1544-t4*0.05638)
       c(3,1)=0.2910+t4*(0.02765-t4*0.012639)
       c(4,1)=cj(1)*5/9
       c(5,1)=cj(2)*5/9
       c(3,2)=0.3080+t4*(0.06120-t4*0.02096)
       c(4,2)=cj(1)/3
       c(5,2)=cj(2)/3
       c(4,3)=cj(1)/9
       c(5,3)=cj(2)/9
       c(5,4)=0.1739+t4*(0.05962-t4*0.00862)
       logAbund = 8.09

       !
       !  Ne IV ion parameters
       !
    case("[Ne IV]")
       a(2,1)=4.84e-4
       a(3,1)=5.54e-3
       a(4,1)=5.21e-1
       a(5,1)=1.27
       a(3,2)=1.48e-6
       a(4,2)=1.15e-1
       a(5,2)=4.0e-1
       a(4,3)=3.93e-1
       a(5,3)=4.37e-1
       a(5,4)=2.68e-9
       t(2)=4.12346e-4
       t(3)=4.12795e-4
       t(4)=6.24346e-4
       t(5)=6.24413e-4
       weight(1)=4
       weight(2)=6
       weight(3)=4
       weight(4)=2
       weight(5)=4
       !
       !  Ne IV ion parameters
       !
       c(2,1)=0.8473+t4*(-0.005832-t4*0.002886)
       c(3,1)=0.5652+t4*(-0.00452-t4*0.001535)
       c(4,1)=0.1496+t4*(9.539e-3-t4*3.437e-3)
       c(5,1)=0.2976+t4*(0.0232-t4*0.0088)
       c(3,2)=1.362+t4*(0.0217-t4*0.0186)
       c(4,2)=0.3057+t4*(0.0859-t4*0.026)
       c(5,2)=0.8014+t4*(0.1364-t4*0.0415)
       c(4,3)=0.308+t4*(0.0391-t4*0.0119)
       c(5,3)=0.4291+t4*(0.1103-t4*0.0336)
       c(5,4)=0.2883+t4*(0.0662-t4*0.0127)
       logAbund = 8.09


       !
       !  Ne V ion parameters
       !
    case("[Ne V]")
       a(2,1)=1.28e-3
       a(3,1)=5.08e-9
       a(4,1)=2.37e-5
       a(3,2)=4.59e-3
       a(4,2)=1.31e-1
       a(5,2)=4.21
       a(4,3)=3.65e-1
       a(5,3)=6.69e-3
       a(5,4)=2.85
       t(2)=4.124e-6
       t(3)=1.1101e-5
       t(4)=3.02915e-4
       t(5)=6.39136e-4
       weight(1)=1
       weight(2)=3
       weight(3)=5
       weight(4)=5
       weight(5)=1
       !
       !  Ne V ion parameters
       !
       cj(1)=1.6175+t4*(0.171-t4*0.01)
       cj(2)=0.3315+t4*(-0.1142+t4*0.034)
       c(2,1)=0.244
       c(3,1)=0.122
       c(4,1)=cj(1)/9
       c(5,1)=cj(2)/9
       c(3,2)=0.578
       c(4,2)=cj(1)/3
       c(5,2)=cj(2)/3
       c(4,3)=cj(1)*5/9
       c(5,3)=cj(2)*5/9
       c(5,4)=1.27-0.02*t4
       logAbund = 8.09


       !
       !  S II ion parameters
       !
    case("[S II]")
       a(2,1)=8.82e-4
       a(3,1)=2.60e-4
       a(4,1)=9.06e-2
       a(5,1)=2.25e-1
       a(3,2)=3.35e-7
       a(4,2)=1.63e-1
       a(5,2)=1.33e-1
       a(4,3)=7.79e-2
       a(5,3)=1.79e-1
       a(5,4)=1.03e-6
       t(2)=1.48530e-4
       t(3)=1.48848e-4
       t(4)=2.45249e-4
       t(5)=2.45718e-4
       weight(1)=4
       weight(2)=4
       weight(3)=6
       weight(4)=2
       weight(5)=4
       !
       !  S II ion parameters
       !
       c(2,1)=3.06+t4*(-0.29+t4*0.02)
       c(3,1)=4.592+t4*(-0.437+t4*0.03)
       c(4,1)=0.78+t4*(-0.0242+t4*0.01)
       c(5,1)=1.33+t4*(0.264-t4*0.08)
       c(3,2)=8.79+t4*(-1.36+t4*0.16)
       c(4,2)=1.4575+t4*(0.111-t4*0.05)
       c(5,2)=3.282+t4*(0.223-t4*0.13)
       c(4,3)=2.502+t4*(0.1431-t4*0.09)
       c(5,3)=4.613+t4*(0.343-t4*0.17)
       c(5,4)=2.383+t4*(0.043-t4*0.05)
       logAbund = 7.27


       !
       !  S III ion parameters
       !
    case("[S III]")
       a(2,1)=4.72e-4
       a(3,1)=4.61e-8
       a(4,1)=5.82e-6
       a(3,2)=2.07e-3
       a(4,2)=2.21e-2
       a(5,2)=7.96e-1
       a(4,3)=5.76e-2
       a(5,3)=1.05e-2
       a(5,4)=2.22
       t(2)=2.972e-6
       t(3)=8.325e-6
       t(4)=1.1320e-4
       t(5)=2.7163e-4
       weight(1)=1
       weight(2)=3
       weight(3)=5
       weight(4)=5
       weight(5)=1
       !
       !  S III ion parameters
       !
       cj(1)=9.903+t4*(-2.017+t4*0.59)
       cj(2)=1.135+t4*0.0522
       c(2,1)=2.6725+t4*(0.01897-t4*0.13)
       c(3,1)=1.053+t4*(0.143-t4*0.05)
       c(4,1)=cj(1)/9
       c(5,1)=cj(2)/9
       c(3,2)=5.71+t4*(0.3181-t4*0.26)
       c(4,2)=cj(1)/3
       c(5,2)=cj(2)/3
       c(4,3)=cj(1)*5/9
       c(5,3)=cj(2)*5/9
       c(5,4)=0.82+t4*(1.424-t4*0.4)
       logAbund = 7.27


       !
       !  Cl II ion parameters
       !
    case("[Cl II]")
       a(2,1)=7.57e-3
       a(3,1)=4.57e-7
       a(4,1)=1.04e-1
       a(5,1)=1.97e-2
       a(3,2)=1.46e-3
       a(4,2)=2.92e-2
       a(5,2)=1.31
       a(4,3)=9.82e-6
       a(5,4)=2.06
       t(2)=6.96e-6
       t(3)=9.965e-6
       t(4)=1.06571e-4
       t(5)=2.78780e-4
       weight(1)=5
       weight(2)=3
       weight(3)=1
       weight(4)=5
       weight(5)=1
       !
       !  Cl II ion parameters
       !
       c(2,1)=2.17
       c(3,1)=0.443
       c(4,1)=3.86*5/9
       c(5,1)=0.456*5/9
       c(3,2)=0.933
       c(4,2)=3.86/3
       c(5,2)=0.456/3
       c(4,3)=3.86/9
       c(5,3)=0.456/9
       c(5,4)=1.15
       logAbund = 5.27
       !
       !  Cl III ion parameters
       !
    case("[Cl III]")
       a(2,1)=4.83e-3
       a(3,1)=7.04e-4
       a(4,1)=3.05e-1
       a(5,1)=7.54e-1
       a(3,2)=3.22e-6
       a(4,2)=3.03e-1
       a(5,2)=3.23e-1
       a(4,3)=1.0e-1
       a(5,3)=3.16e-1
       a(5,4)=7.65e-6
       t(2)=1.8053e-4
       t(3)=1.81186e-4
       t(4)=2.9812e-4
       t(5)=2.9907e-4
       weight(1)=5
       weight(2)=3
       weight(3)=1
       weight(4)=5
       weight(5)=1
       !
       !  Cl III ion parameters
       !
       c(2,1)=1.1120+t4*(0.3837-t4*0.13750)
       c(3,1)=1.6660+t4*(0.5886-t4*0.21062)
       c(4,1)=0.3912+t4*(0.0085+t4*0.01505)
       c(5,1)=0.7810+t4*(0.0213+t4*0.02784)
       c(3,2)=4.0507+t4*(0.7741-t4*0.29264)
       c(4,2)=1.2051+t4*(0.6197-t4*0.18408)
       c(5,2)=1.8324+t4*(0.4803-t4*0.13531)
       c(4,3)=1.3373+t4*(0.2975-t4*0.08182)
       c(5,3)=3.2157+t4*(1.3672-t4*0.40935)
       c(5,4)=1.7478-t4*(0.0450-t4*0.05217)
       logAbund = 5.27

       !
       !  Cl IV ion parameters
       !
    case("[Cl IV]")
       a(2,1)=2.14e-3
       a(3,1)=2.70e-7
       a(4,1)=1.54e-5
       a(3,2)=8.25e-3
       a(4,2)=7.23e-2
       a(5,2)=2.47
       a(4,3)=1.79e-1
       a(5,3)=2.62e-2
       a(5,4)=2.80
       t(2)=4.92e-6
       t(3)=1.3419e-5
       t(4)=1.37676e-4
       t(5)=3.2547e-4
       weight(1)=1
       weight(2)=3
       weight(3)=5
       weight(4)=5
       weight(5)=1
       !
       !  Cl IV ion parameters
       !
       cj(1)=4.702+t4*(0.771-t4*0.01)
       cj(2)=1.712+t4*(0.791-t4*0.25)
       c(2,1)=0.475
       c(3,1)=0.4
       c(4,1)=cj(1)/9
       c(5,1)=cj(2)/9
       c(3,2)=1.5
       c(4,2)=cj(1)/3
       c(5,2)=cj(2)/3
       c(4,3)=cj(1)*5/9
       c(5,3)=cj(2)*5/9
       c(5,4)=0.3388+t4*(1.3214-t4*0.265)
       logAbund = 5.27


       !
       !  Ar III ion parameters
       !
    case("[Ar III]")
       a(2,1)=3.08e-2
       a(3,1)=2.37e-6
       a(4,1)=3.14e-1
       a(5,1)=4.17e-2
       a(3,2)=5.17e-3
       a(4,2)=8.32e-2
       a(5,2)=3.91
       a(4,3)=2.21e-5
       a(5,3)=0.0
       a(5,4)=2.59
       t(2)=1.1121e-5
       t(3)=1.5702e-5
       t(4)=1.40100e-4
       t(5)=3.32657e-4
       weight(1)=5
       weight(2)=3
       weight(3)=1
       weight(4)=5
       weight(5)=1
       !
       !  Ar III ion parameters
       !
       c(2,1)=2.24
       c(3,1)=0.531
       c(4,1)=4.74*5/9
       c(5,1)=0.68*5/9
       c(3,2)=1.18
       c(4,2)=4.74/3
       c(5,2)=0.68/3
       c(4,3)=4.74/9
       c(5,3)=0.68/9
       c(5,4)=0.823
       logAbund = 6.56

       !
       !  Ar IV ion parameters
       !
    case("[Ar IV]")
       a(2,1)=2.23e-2
       a(3,1)=1.77e-3
       a(4,1)=8.62e-1
       a(5,1)=2.11
       a(3,2)=2.3e-5
       a(4,2)=6.03e-1
       a(5,2)=7.89e-1
       a(4,3)=0.119
       a(5,3)=0.598
       a(5,4)=4.94e-5
       t(2)=2.10904e-4
       t(3)=2.12193e-4
       t(4)=3.48555e-4
       t(5)=3.50326e-4
       weight(1)=4
       weight(2)=4
       weight(3)=6
       weight(4)=2
       weight(5)=4
       !
       !  Ar IV ion parameters
       !
       c(2,1)=2.2661-t4*(1.2805-t4*0.32167)
       c(3,1)=3.3993-t4*(1.9217-t4*0.48293)
       c(4,1)=0.1637-t4*(0.0351-t4*0.01790)
       c(5,1)=0.3356-t4*(0.0817-t4*0.03930)
       c(3,2)=6.4696-t4*(0.3631-t4*0.04479)
       c(4,2)=1.4523+t4*(0.4098-t4*0.15965)
       c(5,2)=2.3424+t4*(0.2396-t4*0.11250)
       c(4,3)=1.7193+t4*(0.1391-t4*0.07235)
       c(5,3)=3.9465+t4*(0.7872-t4*0.30578)
       c(5,4)=2.1475+t4*(0.1047+t4*0.09440)
       logAbund = 6.56

       !
       !  Ar V ion parameters
       !
    case("[Ar V]")
       a(2,1)=7.99e-3
       a(3,1)=1.24e-6
       a(4,1)=3.50e-5
       a(3,2)=2.72e-2
       a(4,2)=0.204
       a(5,2)=6.55
       a(4,3)=0.476
       a(5,3)=5.69e-2
       a(5,4)=3.29
       t(2)=7.639e-6
       t(3)=2.0292e-5
       t(4)=1.62994e-4
       t(5)=3.79125e-4
       weight(1)=1
       weight(2)=3
       weight(3)=5
       weight(4)=5
       weight(5)=1
       !
       !  Ar V ion parameters
       !
       cj(1)=5.2075+t4*(-1.985+t4*0.55)
       cj(2)=1.133+t4*(0.127-t4*0.09)
       c(2,1)=0.257
       c(3,1)=0.32
       c(4,1)=cj(1)/9
       c(5,1)=cj(2)/9
       c(3,2)=1.04
       c(4,2)=cj(1)/3
       c(5,2)=cj(2)/3
       c(4,3)=cj(1)*5/9
       c(5,3)=cj(2)*5/9
       c(5,4)=1.27+t4*(-0.02+t4*2.4e-5)
       logAbund = 6.56


    case DEFAULT
       write(*,'(a,a10)') "! Ion not recognised: ",ion
       stop
    end select



    !
    !     _________________________________________________________________
    !     |                                                               |
    !     |                            SOLVE                              |
    !     |                                                               |
    !     |     TS:       square root of temperature                      |
    !     |     POP1:     storage array for level populations             |
    !     |     POP:      level population array                          |
    !     |     A:        radiative transition probabilities              |
    !     |     C:        collision strengths                             |
    !     |     Q:        collision transition probabilities              |
    !     |     E:        relative energy separation                      |
    !     |     KT,HC,CQ: "constants"                                     |
    !     |                                                               |
    !     -----------------------------------------------------------------
    !
    !
    kt=(1.3807e-16)*tem
    hc=1.9865e-08
    cq=8.629e-06
    ts=sqrt(tem)
    !
    !  Transition energy differences
    !
    do i=1,5
       do j=i+1,5
          e(j,i)=hc*(t(j)-t(i))
       enddo
    enddo
    !
    !  Collision transition probabilities (if i=j, Q=0)
    !
    do i=1,5
       do j=1,5
          if(j .gt. i) then
             q(i,j)=cq*c(j,i)*exp(-e(j,i)/kt)/(weight(i)*ts)
          else if(j .eq. i) then
             q(i,j)=0.0
          else
             q(i,j)=cq*c(i,j)/(weight(i)*ts)
          endif
       enddo
    enddo
    !
    !  Critical density calculation
    !
    do i=2,5
       asum=0.0
       qsum=0.0
       do j=1,i-1
          asum=asum+a(i,j)
          qsum=qsum+q(i,j)
       enddo
       ncrit(i)=asum/qsum
    enddo
    !
    !  Matrix manipulation
    !
    do  i=1,5
       do j=1,5
          if(j .eq. i) then
             a1(i,j)=0.0
          else if(j .gt. i) then
             a1(i,j)=den*q(i,j)
          else
             a1(i,j)=den*q(i,j)+a(i,j)
          endif
       enddo
    enddo
    do i=1,5
       b1(i)=0.0
       do j=1,5 
          b1(i)=b1(i)+a1(i,j)
       enddo
    enddo

    do i=1,5
       b2(i)=-a1(1,i)
       a1(i,i)=-b1(i)
    enddo

    do i=1,5
       do  j=1,5
          do  k=1,5
             a2(j,k)=a1(j,k)
          enddo
       enddo

       do j=1,5
          a2(i,j)=b2(j)
       enddo
       dtl(i)=1

       do j=2,4
          t1=a2(j,j)
          if(t1 .eq. 0.0) then
             write(filout,*)'Divide by 0 during T1 normalization in SOLVE'
             stop
          endif
          dtl(i)=dtl(i)*t1
          do k=j,5
             a2(j,k)=a2(j,k)/t1
          enddo
          do k=j+1,5
             t2=a2(k,j)
             do l=j,5 
                a2(k,l)=a2(k,l)-t2*a2(j,l)
             enddo
          enddo
       enddo

       dtl(i)=dtl(i)*a2(5,5)
    enddo

    pop1(1)=1
    do  i=2,5
       pop1(i)=dtl(i)/dtl(1)
    enddo

    t1=0
    do i=1,5
       t1=t1+pop1(i)
    enddo
    !
    !  Finally, level populations, emissivities
    !
    do i=1,5
       pop(i)=pop1(i)/t1
       do  j=1,5
          jl(i,j)=pop(i)*a(i,j)*e(i,j)
          if(jl(i,j) .lt. 1.0e-34) jl(i,j)=1.0e-34
       enddo
    enddo
    do i=2,5
       do j=1,4
          if(j .lt. i) then
             xl(j)=1./(t(i)-t(j))
             if(xl(j) .ge. 1.e5) xl(j)=xl(j)*1.e-4
             xi(j)=jl(i,j)/den
             if(xi(j) .lt. 1.e-35) xi(j)=0.
          endif
       enddo
       if(i .eq. 2) then
          write(filout,11) (xl(k),k=1,i-1),(i,k,k=1,i-1), (xi(k),k=1,i-1)
       elseif(i .eq. 3) then
          write(filout,12) (xl(k),k=1,i-1),(i,k,k=1,i-1), (xi(k),k=1,i-1)
       elseif(i .eq. 4) then
          write(filout,13) (xl(k),k=1,i-1),(i,k,k=1,i-1), (xi(k),k=1,i-1)
       elseif(i .eq. 5) then
          write(filout,14) (xl(k),k=1,i-1),(i,k,k=1,i-1), (xi(k),k=1,i-1)
       endif

       if (i == nUpper) etaLine = xi(nLower)
    enddo
!!$   10 format(10x,'Level Populations - Critical Densities (cm**-3):',/, &
!!$     10x,'Level 1: ',1pe8.2,/, &
!!$     10x,'Level 2: ',1pe8.2,13x,1pe8.2,/, &
!!$     10x,'Level 3: ',1pe8.2,13x,1pe8.2,/, &
!!$     10x,'Level 4: ',1pe8.2,13x,1pe8.2,/, &
!!$     10x,'Level 5: ',1pe8.2,13x,1pe8.2)
   11 format(5x,'Wavelength',5x,f11.2,/, &
            1x,'Upper--Lower Levels',5x,'(',i1,'--',i1,') ',/, &
            2x,'Volume Emissivity ',1pe11.3,/)
   12 format(20x,2f11.2,/,20x,2(5x,'(',i1,'--',i1,')'),/, &
            20x,2(1pe11.3),/)
   13 format(20x,3f11.2,/,20x,3(5x,'(',i1,'--',i1,')'),/, &
            20x,3(1pe11.3),/)
   14 format(20x,4f11.2,/,20x,4(5x,'(',i1,'--',i1,')'),/, &
            20x,4(1pe11.3),/)
!!$   15 format(1x,'H-beta volume emissivity: ',1pe9.2,' N(H+)Ne', &
!!$            ' erg/s',/)


    ! etaLine is per ion density per electron density * 4 pi

    nH = (0.707 * massDensity)/1.67e-24

    nIon = nH * 10.**(logAbund - 12.) * ionFrac
    etaLine = etaLine * (electronDensity * nIon) / (4.*3.141592654)


  end subroutine solveForbidden



  

    
end module opacity_lte_mod
