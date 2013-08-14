module setuppdr_mod
use healpix_types
implicit none



contains

  subroutine setupPDR()
    use unix_mod, only: unixGetenv

    implicit none
    character(len=200) :: dataDirectory, dataPath
    character(len=400) :: C12Oinput, CIIinput, CIinput, OIinput
    integer :: i
    integer :: cii_nlev, cii_ntemp, ci_nlev, ci_ntemp, oi_nlev, oi_ntemp
    integer :: c12o_nlev, c12o_ntemp



  real(kind=dp), allocatable :: COEFF(:)
  real(kind=dp), allocatable :: ENERGIES(:), WEIGHTS(:)
  real(kind=dp), allocatable :: A_COEFFS(:,:), B_COEFFS(:,:), C_COEFFS(:,:)
  real(kind=dp), allocatable :: FREQUENCIES(:,:), TEMPERATURES(:)
  real(kind=dp), allocatable :: H_COL(:,:,:)
  real(kind=dp), allocatable :: EL_COL(:,:,:)
  real(kind=dp), allocatable :: HE_COL(:,:,:)
  real(kind=dp), allocatable :: H2_COL(:,:,:)
  real(kind=dp), allocatable :: PH2_COL(:,:,:)
  real(kind=dp), allocatable :: OH2_COL(:,:,:)

  real(kind=dp), allocatable :: tau_ij(:)
  real(kind=dp), allocatable :: field(:,:)
  real(kind=dp), allocatable :: transition_CII(:,:),transition_CI(:,:)
  real(kind=dp), allocatable :: transition_OI(:,:),transition_C12O(:,:)
!CII cooling variables
  real(kind=dp) :: CII_COOLING
  real(kind=dp), allocatable :: CII_ENERGIES(:),CII_WEIGHTS(:)
  real(kind=dp), allocatable :: CII_A_COEFFS(:,:),CII_B_COEFFS(:,:),CII_C_COEFFS(:,:)
  real(kind=dp), allocatable :: CII_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: CII_TEMPERATURES(:,:), CII_HP(:,:,:)
  real(kind=dp), allocatable :: CII_H(:,:,:),CII_EL(:,:,:)
  real(kind=dp), allocatable :: CII_HE(:,:,:),CII_H2(:,:,:)
  real(kind=dp), allocatable :: CII_PH2(:,:,:),CII_OH2(:,:,:)
!  real(kind=dp), allocatable :: CII_LINE(:,:,:)!,CII_POP(:,:)
!CI cooling variables
  real(kind=dp) :: CI_COOLING
  real(kind=dp), allocatable :: CI_ENERGIES(:),CI_WEIGHTS(:)
  real(kind=dp), allocatable :: CI_A_COEFFS(:,:),CI_B_COEFFS(:,:),CI_C_COEFFS(:,:)
  real(kind=dp), allocatable :: CI_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: CI_TEMPERATURES(:,:),CI_HP(:,:,:)
  real(kind=dp), allocatable :: CI_H(:,:,:),CI_EL(:,:,:)
  real(kind=dp), allocatable :: CI_HE(:,:,:),CI_H2(:,:,:)
  real(kind=dp), allocatable :: CI_PH2(:,:,:),CI_OH2(:,:,:)
!  real(kind=dp), allocatable :: CI_LINE(:,:,:)!,CI_POP(:,:)
!OI cooling variables
  real(kind=dp) :: OI_COOLING
  real(kind=dp), allocatable :: OI_ENERGIES(:),OI_WEIGHTS(:)
  real(kind=dp), allocatable :: OI_A_COEFFS(:,:),OI_B_COEFFS(:,:),OI_C_COEFFS(:,:)
  real(kind=dp), allocatable :: OI_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: OI_TEMPERATURES(:,:),OI_HP(:,:,:)
  real(kind=dp), allocatable :: OI_H(:,:,:),OI_EL(:,:,:)
  real(kind=dp), allocatable :: OI_HE(:,:,:),OI_H2(:,:,:)
  real(kind=dp), allocatable :: OI_PH2(:,:,:),OI_OH2(:,:,:)
!  real(kind=dp), allocatable :: OI_LINE(:,:,:)!,OI_POP(:,:)
!12CO cooling variables
  real(kind=dp) :: C12O_COOLING
  real(kind=dp), allocatable :: C12O_ENERGIES(:),C12O_WEIGHTS(:)
  real(kind=dp), allocatable :: C12O_A_COEFFS(:,:),C12O_B_COEFFS(:,:),C12O_C_COEFFS(:,:)
  real(kind=dp), allocatable :: C12O_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: C12O_TEMPERATURES(:,:),C12O_HP(:,:,:)
  real(kind=dp), allocatable :: C12O_H(:,:,:),C12O_EL(:,:,:)
  real(kind=dp), allocatable :: C12O_HE(:,:,:),C12O_H2(:,:,:)
  real(kind=dp), allocatable :: C12O_PH2(:,:,:),C12O_OH2(:,:,:)
  real(kind=dp) :: SCO_GRID(1:8, 1:6)
!  real(kind=dp), allocatable :: C12O_LINE(:,:,:)!,C12O_POP(:,:)

  cii_nlev = 5
  cii_ntemp = 18
  ci_nlev = 5
  ci_ntemp = 29
  oi_nlev = 5
  oi_ntemp = 27
  c12o_nlev = 41
  c12o_ntemp = 25

!load SCO_GRID data [UCL_PDR]
SCO_GRID(1:8,1) = (/0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00/)
SCO_GRID(1:8,2) = (/-8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00/)
SCO_GRID(1:8,3) = (/-1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00/)
SCO_GRID(1:8,4) = (/-4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00/)
SCO_GRID(1:8,5) = (/-1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00/)
SCO_GRID(1:8,6) = (/-3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/)


ALLOCATE(CII_ENERGIES(1:CII_NLEV))
ALLOCATE(CII_WEIGHTS(1:CII_NLEV))
ALLOCATE(CII_A_COEFFS(1:CII_NLEV,1:CII_NLEV))
ALLOCATE(CII_B_COEFFS(1:CII_NLEV,1:CII_NLEV))
ALLOCATE(CII_FREQUENCIES(1:CII_NLEV,1:CII_NLEV))
ALLOCATE(CII_TEMPERATURES(1:7,1:CII_NTEMP))
ALLOCATE(CII_HP(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_H(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_EL(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_HE(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_H2(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_PH2(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
ALLOCATE(CII_OH2(1:CII_NLEV,1:CII_NLEV,1:CII_NTEMP))
!CI
ALLOCATE(CI_ENERGIES(1:CI_NLEV))
ALLOCATE(CI_WEIGHTS(1:CI_NLEV))
ALLOCATE(CI_A_COEFFS(1:CI_NLEV,1:CI_NLEV))
ALLOCATE(CI_B_COEFFS(1:CI_NLEV,1:CI_NLEV))
ALLOCATE(CI_FREQUENCIES(1:CI_NLEV,1:CI_NLEV))
ALLOCATE(CI_TEMPERATURES(1:7,1:CI_NTEMP))
ALLOCATE(CI_HP(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_H(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_EL(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_HE(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_H2(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_PH2(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
ALLOCATE(CI_OH2(1:CI_NLEV,1:CI_NLEV,1:CI_NTEMP))
!OI
ALLOCATE(OI_ENERGIES(1:OI_NLEV))
ALLOCATE(OI_WEIGHTS(1:OI_NLEV))
ALLOCATE(OI_A_COEFFS(1:OI_NLEV,1:OI_NLEV))
ALLOCATE(OI_B_COEFFS(1:OI_NLEV,1:OI_NLEV))
ALLOCATE(OI_FREQUENCIES(1:OI_NLEV,1:OI_NLEV))
ALLOCATE(OI_TEMPERATURES(1:7,1:OI_NTEMP))
ALLOCATE(OI_HP(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_H(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_EL(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_HE(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_H2(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_PH2(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
ALLOCATE(OI_OH2(1:OI_NLEV,1:OI_NLEV,1:OI_NTEMP))
!12CO
ALLOCATE(C12O_ENERGIES(1:C12O_NLEV))
ALLOCATE(C12O_WEIGHTS(1:C12O_NLEV))
ALLOCATE(C12O_A_COEFFS(1:C12O_NLEV,1:C12O_NLEV))
ALLOCATE(C12O_B_COEFFS(1:C12O_NLEV,1:C12O_NLEV))
ALLOCATE(C12O_FREQUENCIES(1:C12O_NLEV,1:C12O_NLEV))
ALLOCATE(C12O_TEMPERATURES(1:7,1:C12O_NTEMP))
ALLOCATE(C12O_HP(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_H(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_EL(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_HE(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_H2(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_PH2(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))
ALLOCATE(C12O_OH2(1:C12O_NLEV,1:C12O_NLEV,1:C12O_NTEMP))


!    C12Oinput = "12co.dat"; CIIinput = "12c+.dat" ; CIinput = "12c.dat" ; OIinput = "16o.dat"

    
    call unixGetenv("TORUS_DATA", dataPath, i)
    
    dataDirectory = trim(dataPath)

    C12Oinput = trim(dataDirectory) // "/pdrdata/12co.dat"
    CIIinput = trim(dataDirectory) // "/pdrdata/12c+.dat"
    CIinput = trim(dataDirectory) // "/pdrdata/12c.dat"
    OIinput = trim(dataDirectory) // "/pdrdata/16o.dat"
    
!    print *, C12Oinput
!    print *, CIIinput
!    print *, CIinput
!    print *, OIinput
    

    call readinput(C12Oinput,C12O_NLEV,C12O_NTEMP,C12O_ENERGIES,C12O_WEIGHTS,&
         &     C12O_A_COEFFS,C12O_B_COEFFS,C12O_FREQUENCIES,C12O_TEMPERATURES,&
         &     C12O_H,C12O_HP,C12O_EL,C12O_HE,C12O_H2,C12O_PH2,C12O_OH2)
    call readinput(CIIinput,CII_NLEV,CII_NTEMP,CII_ENERGIES,CII_WEIGHTS,&
         &     CII_A_COEFFS,CII_B_COEFFS,CII_FREQUENCIES,CII_TEMPERATURES,&
         &     CII_H,CII_HP,CII_EL,CII_HE,CII_H2,CII_PH2,CII_OH2)
    call readinput(CIinput,CI_NLEV,CI_NTEMP,CI_ENERGIES,CI_WEIGHTS,&
         &     CI_A_COEFFS,CI_B_COEFFS,CI_FREQUENCIES,CI_TEMPERATURES,&
         &     CI_H,CI_HP,CI_EL,CI_HE,CI_H2,CI_PH2,CI_OH2)
    call readinput(OIinput,OI_NLEV,OI_NTEMP,OI_ENERGIES,OI_WEIGHTS,&
         &     OI_A_COEFFS,OI_B_COEFFS,OI_FREQUENCIES,OI_TEMPERATURES,&
         &     OI_H,OI_HP,OI_EL,OI_HE,OI_H2,OI_PH2,OI_OH2)



!allocations..........
!ALLOCATE(COEFF(1:NTEMP))
!ALLOCATE(ENERGIES(1:NLEV))
!ALLOCATE(WEIGHTS(1:NLEV))
!ALLOCATE(A_COEFFS(1:NLEV,1:NLEV))
!ALLOCATE(B_COEFFS(1:NLEV,1:NLEV))
!ALLOCATE(FREQUENCIES(1:NLEV,1:NLEV))
!ALLOCATE(TEMPERATURES(1:NTEMP))
!ALLOCATE(H_COL(1:NLEV,1:NLEV,1:NTEMP))
!ALLOCATE(EL_COL(1:NLEV,1:NLEV,1:NTEMP))
!ALLOCATE(HE_COL(1:NLEV,1:NLEV,1:NTEMP))
!ALLOCATE(H2_COL(1:NLEV,1:NLEV,1:NTEMP))
!ALLOCATE(PH2_COL(1:NLEV,1:NLEV,1:NTEMP))
!ALLOCATE(OH2_COL(1:NLEV,1:NLEV,1:NTEMP))
!CII


!allocate(C_COEFFS(1:NLEV,1:NLEV)) !don't need to C_COEFFS=0.0d0 because it's done in subroutine





  end subroutine setupPDR



end module setuppdr_mod






