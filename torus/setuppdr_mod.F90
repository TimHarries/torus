module setuppdr_mod
#ifdef PDR
!use healpix_module
use healpix_guts
use grid_mod
implicit none


real(kind=dp), allocatable :: CII_ENERGIES(:),CII_WEIGHTS(:)
real(kind=dp), allocatable :: CII_A_COEFFS(:,:),CII_B_COEFFS(:,:)!,CII_C_COEFFS(:,:)
real(kind=dp), allocatable :: CII_FREQUENCIES(:,:)
real(kind=dp), allocatable :: CII_TEMPERATURES(:,:), CII_HP(:,:,:)
real(kind=dp), allocatable :: CII_H(:,:,:),CII_EL(:,:,:)
real(kind=dp), allocatable :: CII_HE(:,:,:),CII_H2(:,:,:)
real(kind=dp), allocatable :: CII_PH2(:,:,:),CII_OH2(:,:,:)

  real(kind=dp), allocatable :: CI_ENERGIES(:),CI_WEIGHTS(:)
  real(kind=dp), allocatable :: CI_A_COEFFS(:,:),CI_B_COEFFS(:,:)!,CI_C_COEFFS(:,:)
  real(kind=dp), allocatable :: CI_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: CI_TEMPERATURES(:,:),CI_HP(:,:,:)
  real(kind=dp), allocatable :: CI_H(:,:,:),CI_EL(:,:,:)
  real(kind=dp), allocatable :: CI_HE(:,:,:),CI_H2(:,:,:)
  real(kind=dp), allocatable :: CI_PH2(:,:,:),CI_OH2(:,:,:)


  real(kind=dp), allocatable :: OI_ENERGIES(:),OI_WEIGHTS(:)
  real(kind=dp), allocatable :: OI_A_COEFFS(:,:),OI_B_COEFFS(:,:)!,OI_C_COEFFS(:,:)
  real(kind=dp), allocatable :: OI_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: OI_TEMPERATURES(:,:),OI_HP(:,:,:)
  real(kind=dp), allocatable :: OI_H(:,:,:),OI_EL(:,:,:)
  real(kind=dp), allocatable :: OI_HE(:,:,:),OI_H2(:,:,:)
  real(kind=dp), allocatable :: OI_PH2(:,:,:),OI_OH2(:,:,:)

  real(kind=dp), allocatable :: C12O_ENERGIES(:),C12O_WEIGHTS(:)
  real(kind=dp), allocatable :: C12O_A_COEFFS(:,:),C12O_B_COEFFS(:,:)!,C12O_C_COEFFS(:,:)
  real(kind=dp), allocatable :: C12O_FREQUENCIES(:,:)
  real(kind=dp), allocatable :: C12O_TEMPERATURES(:,:),C12O_HP(:,:,:)
  real(kind=dp), allocatable :: C12O_H(:,:,:),C12O_EL(:,:,:)
  real(kind=dp), allocatable :: C12O_HE(:,:,:),C12O_H2(:,:,:)
  real(kind=dp), allocatable :: C12O_PH2(:,:,:),C12O_OH2(:,:,:)
  real(kind=dp) :: SCO_GRID(1:8, 1:6)

  integer, parameter :: cii_nlev=5
  integer, parameter :: cii_ntemp=18
  integer, parameter :: ci_nlev=5
  integer, parameter :: ci_ntemp=29
  integer, parameter :: oi_nlev=5
  integer, parameter :: oi_ntemp=27
  integer, parameter :: c12o_nlev=41
  integer, parameter :: c12o_ntemp=25

  INTEGER :: NH,ND,NH2,NHD,NPROTON,NH2O,NHe, &
       &        NMG,NMGx,NN,NFE,NFEx,NSI,NSIx,NCA,NCAx,NCAxx,NS,NSx,NCS, &
       &        NCL,NCLx,NH2x,NHEx,NOx,NNx,NNA,NNAx,NCH,NCH2,NOH,NO2, &
       &        NH3x, NH3Ox, NHCOx!, NOSH
  integer :: NCO, NC, NCX, NO

contains


!  subroutine 


  subroutine setupPDR(grid, reactant, product, alpha, beta, gamma, &
       rate, duplicate, rtmin, rtmax, nc12o, nci, ncii, noi, nelect)
    use unix_mod, only: unixGetenv
!    use healpix_guts
!    use read_input, only: readinput
    implicit none

!  integer, parameter :: cii_nlev=5
!  integer, parameter :: cii_ntemp=18
!  integer, parameter :: ci_nlev=5
!  integer, parameter :: ci_ntemp=29
!  integer, parameter :: oi_nlev=5
!  integer, parameter :: oi_ntemp=27
!  integer, parameter :: c12o_nlev=41
!  integer, parameter :: c12o_ntemp=25
    type(gridtype) :: grid
    character(len=200) :: dataDirectory, dataPath
    character(len=400) :: C12Oinput, CIIinput, CIinput, OIinput
    integer :: i
!    integer :: cii_nlev, cii_ntemp, ci_nlev, ci_ntemp, oi_nlev, oi_ntemp
!    integer :: c12o_nlev, c12o_ntemp
    integer :: nc12o, nci, ncii, noi, nelect
    

!  real(kind=dp), allocatable :: COEFF(:)
!  real(kind=dp), allocatable :: ENERGIES(:), WEIGHTS(:)
 ! real(kind=dp), allocatable :: A_COEFFS(:,:), B_COEFFS(:,:), C_COEFFS(:,:)
 ! real(kind=dp), allocatable :: FREQUENCIES(:,:), TEMPERATURES(:)
  !real(kind=dp), allocatable :: H_COL(:,:,:)
  !real(kind=dp), allocatable :: EL_COL(:,:,:)
  !real(kind=dp), allocatable :: HE_COL(:,:,:)
  !real(kind=dp), allocatable :: H2_COL(:,:,:)
  !real(kind=dp), allocatable :: PH2_COL(:,:,:)
  !real(kind=dp), allocatable :: OH2_COL(:,:,:)

!!  real(kind=dp), allocatable :: tau_ij(:)
 ! real(kind=dp), allocatable :: field(:,:)
!  real(kind=dp), allocatable :: transition_CII(:,:),transition_CI(:,:)
!  real(kind=dp), allocatable :: transition_OI(:,:),transition_C12O(:,:)
!CII cooling variables
!  real(kind=dp) :: CII_COOLING
!  real(kind=dp), allocatable :: CII_ENERGIES(:),CII_WEIGHTS(:)
!  real(kind=dp), allocatable :: CII_A_COEFFS(:,:),CII_B_COEFFS(:,:)!,CII_C_COEFFS(:,:)
! real(kind=dp), allocatable :: CII_FREQUENCIES(:,:)
!  real(kind=dp), allocatable :: CII_TEMPERATURES(:,:), CII_HP(:,:,:)
!  real(kind=dp), allocatable :: CII_H(:,:,:),CII_EL(:,:,:)
!  real(kind=dp), allocatable :: CII_HE(:,:,:),CII_H2(:,:,:)
!  real(kind=dp), allocatable :: CII_PH2(:,:,:),CII_OH2(:,:,:)
!  real(kind=dp), allocatable :: CII_LINE(:,:,:)!,CII_POP(:,:)
!CI cooling variables
!  real(kind=dp) :: CI_COOLING

!  real(kind=dp), allocatable :: CI_LINE(:,:,:)!,CI_POP(:,:)
!OI cooling variables
!  real(kind=dp) :: OI_COOLING

!  real(kind=dp), allocatable :: OI_LINE(:,:,:)!,OI_POP(:,:)
!12CO cooling variables
!  real(kind=dp) :: C12O_COOLING

!

  integer :: nspec
  real(double) :: dummyAbundance(1:33)
  
  integer :: nreac
  real(double), allocatable :: rate(:), alpha(:), beta(:), gamma(:)
  real(double), allocatable :: rtmin(:), rtmax(:)
  integer, allocatable :: duplicate(:)
  character(len=10), allocatable :: product(:,:), reactant(:,:)

!  real(kind=dp), allocatable :: C12O_LINE(:,:,:)!,C12O_POP(:,:)
  nreac = 329



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





!load SCO_GRID data [UCL_PDR]
SCO_GRID(1:8,1) = (/0.000D+00,-1.408D-02,-1.099D-01,-4.400D-01,-1.154D+00,-1.888D+00,-2.760D+00,-4.001D+00/)
SCO_GRID(1:8,2) = (/-8.539D-02,-1.015D-01,-2.104D-01,-5.608D-01,-1.272D+00,-1.973D+00,-2.818D+00,-4.055D+00/)
SCO_GRID(1:8,3) = (/-1.451D-01,-1.612D-01,-2.708D-01,-6.273D-01,-1.355D+00,-2.057D+00,-2.902D+00,-4.122D+00/)
SCO_GRID(1:8,4) = (/-4.559D-01,-4.666D-01,-5.432D-01,-8.665D-01,-1.602D+00,-2.303D+00,-3.146D+00,-4.421D+00/)
SCO_GRID(1:8,5) = (/-1.303D+00,-1.312D+00,-1.367D+00,-1.676D+00,-2.305D+00,-3.034D+00,-3.758D+00,-5.077D+00/)
SCO_GRID(1:8,6) = (/-3.883D+00,-3.888D+00,-3.936D+00,-4.197D+00,-4.739D+00,-5.165D+00,-5.441D+00,-6.446D+00/)



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


!This is hardwired for species_reduced ! ! ! ! ! ! ! 
    nspec = 33
!    allocate(dummyAbundance(1:nspec))
!    print *, "READING SPECIES"
    call read_species(nspec, dummyAbundance, nc12o, nci, ncii, noi, nelect)
    call setupDummyAbundances(grid%octreeroot, dummyabundance)
!    print *, dummyAbundance(1:10)
    
    allocate(reactant(1:nreac,1:3))
    allocate(product(1:nreac,1:4))
    allocate(alpha(1:nreac))
    allocate(beta(1:nreac))
    allocate(gamma(1:nreac))
    allocate(rate(1:nreac))
    allocate(rtmin(1:nreac))
    allocate(rtmax(1:nreac))
    allocate(duplicate(1:nreac))
    
!    print *, "READING RATES"
    call read_rates(nreac, reactant, product, alpha, beta, gamma, rate, duplicate, rtmin, rtmax)
!    print *, product(1:10, :)
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

  recursive subroutine setupDummyAbundances(thisOctal, dummyabundance)
    type(octal), pointer :: thisOctal, child
    real(double), intent(in) :: dummyabundance(1:33)
    integer :: subcell, i
    do subcell = 1, thisoctal%maxchildren
       if (thisoctal%haschild(subcell)) then
          ! find the child
          do i = 1, thisoctal%nchildren, 1
             if (thisoctal%indexchild(i) == subcell) then
                child => thisoctal%child(i)
                call setupDummyAbundances(child, dummyabundance)
                exit
             end if
          end do
       else
          thisOctal%abundance(subcell, :) = 0.d0
          if(.not. thisOctal%ionfrac(subcell, 2) > 0.99d0) then
             thisOctal%abundance(subcell,:) = dummyabundance(:)
          else
             thisOctal%abundance(subcell, :) = 0.d0
          end if

       end if

    end do
  end subroutine setupDummyAbundances

!C***********************************************************************
!C     Read in the chemical reaction rates and the species, masses and
!C     initial abundances (if specified). The  rates and species files
!C     are assumed to have comma separated values (CSV) format. This is
!C     in line with the Rate05 formatting, removing the need for file-
!C     dependent FORMAT statements.
!C***********************************************************************
!
      SUBROUTINE READ_RATES(NREAC,REAC,PROD,ALPHA,BETA,GAMMA,RATE,&
                           &DUPLICATE,RTMIN,RTMAX)
!T.Bell
!use definitions
!use healpix_types
!        use healpix_guts
!use pdr_mod
use unix_mod, only: unixGetenv
!use maincode_module, only : nreac, duplicate, alpha, beta, gamma, rate, &
!                           & rtmin, rtmax, reac, prod

      IMPLICIT NONE
      INTEGER, intent(in) :: NREAC
      integer, intent(out) :: DUPLICATE(1:nreac)
      real(double), intent(out) :: ALPHA(1:nreac),BETA(1:nreac),&
              &GAMMA(1:nreac),RATE(1:nreac),RTMIN(1:nreac),RTMAX(1:nreac)
      CHARACTER(len=10), intent(out) :: REAC(1:nreac,1:3),PROD(1:nreac,1:4)

      INTEGER :: I,J,N,RATEFILE
      CHARACTER(len=1) :: CLEM
      character(len=400) :: rateData
      character(len=200) :: dataDirectory, dataPath
!   allocate(reac(1:nreac,1:3))
!   allocate(prod(1:nreac,1:4))
!   allocate(rate(1:nreac))
!   allocate(alpha(1:nreac))
!   allocate(beta(1:nreac))
!   allocate(gamma(1:nreac))
!   allocate(rtmin(1:nreac))
!   allocate(rtmax(1:nreac))
!   allocate(duplicate(1:nreac))


      RATEFILE = 2

!C     Initialize the variables and read in the ratefile data. Check that
!C     the value of NREAC agrees with the number of reactions in the file
!C     and produce an error message if not.
      REAC="          "
      PROD="          "
      ALPHA=0.0D0
      BETA=0.0D0
      GAMMA=0.0D0
      RTMIN=0.0D0
      RTMAX=0.0D0
      DUPLICATE=0

      RATE=0.0D0

!#ifdef XRAYS
!      OPEN(RATEFILE,FILE="Xrates.d",STATUS="OLD")
!#else
!      OPEN(RATEFILE,FILE="rates.d",STATUS="OLD")
!#endif        
!#ifdef REDUCED
!#ifdef XRAYS
!      OPEN(RATEFILE,FILE="rates_Xreduced.d",STATUS="OLD")
!#else
      call unixGetenv("TORUS_DATA", dataPath, i)
      
      dataDirectory = trim(dataPath)
      
      rateData = trim(dataDirectory) // "/pdrdata/rates_reduced.d"

      OPEN(RATEFILE,FILE=rateData,STATUS="OLD")
!#endif
!#endif
!#ifdef FULL
!#ifdef XRAYS
!      OPEN(RATEFILE,FILE="rates_Xfull.d",STATUS="OLD")
!#else
!      OPEN(RATEFILE,FILE="rates_full.d",STATUS="OLD")
!#endif
!#endif
!#ifdef MYNETWORK
!#ifdef XRAYS
!      OPEN(RATEFILE,FILE="rates_Xmynetwork.d",STATUS="OLD")
!#else
!      OPEN(RATEFILE,FILE="rates_mynetwork.d",STATUS="OLD")
!#endif
!#endif
      REWIND(RATEFILE)
      DO I=1,NREAC
         READ(RATEFILE,*,END=1) N,(REAC(I,J),J=1,3),(PROD(I,J),J=1,4),&
     &                          ALPHA(I),BETA(I),GAMMA(I), &
     &                          CLEM,RTMIN(I),RTMAX(I)
         IF(CLEM.NE."") CLEM=""

!C     Check for duplicate reactions and set the DUPLICATE counter to the
!C     appropriate value. Adjust their minimum temperatures so that the
!C     temperature ranges are adjacent.
         IF(I.GT.1) THEN
           IF(REAC(I,1).EQ.REAC(I-1,1) .AND. &
     &     REAC(I,2).EQ.REAC(I-1,2) .AND. REAC(I,3).EQ.REAC(I-1,3) .AND. &
     &     PROD(I,1).EQ.PROD(I-1,1) .AND. PROD(I,2).EQ.PROD(I-1,2) .AND. &
     &     PROD(I,3).EQ.PROD(I-1,3) .AND. PROD(I,4).EQ.PROD(I-1,4)) THEN 
            IF(DUPLICATE(I-1).EQ.0) DUPLICATE(I-1)=1
            DUPLICATE(I)=DUPLICATE(I-1)+1
            RTMIN(I)=RTMAX(I-1)
           ELSE
            DUPLICATE(I)=0
           ENDIF
         ELSE
            DUPLICATE(I)=0
         ENDIF

!C     Check for negative gamma values as they could cause problems when
!C     calculating abundances. Produce a warning message if they occur.
         IF(GAMMA(I).LT.0.0D0 .and. myrankglobal == 1) THEN
          write(6,*) 'Negative gamma factor in rate',N
          WRITE(10,"('Negative gamma factor in rate',I5,' (',F8.1,')')")&
     &         N,GAMMA(I)
         ENDIF
      ENDDO
      I=I-1
      READ(RATEFILE,*,END=1)
      I=I+1
 1    IF(I.NE.NREAC .and. myrankglobal == 1) THEN
         write(6,*) 'ERROR! Number of reactions (NREAC) does not match ', &
     &           'the number of entries in the ratefile'
         STOP
      ENDIF

      CLOSE(RATEFILE)
      RETURN
      END SUBROUTINE
!C-----------------------------------------------------------------------


!C***********************************************************************
!C     Read in the chemical reaction rates and the species, masses and
!C     initial abundances (if specified). The  rates and species files
!C     are assumed to have comma separated values (CSV) format. This is
!C     in line with the Rate05 formatting, removing the need for file-
!C     dependent FORMAT statements.
!C***********************************************************************
!
!C-----------------------------------------------------------------------
!C     Read in the species data, including initial fractional abundances
!C     (3rd column) and their masses (4th column). Check that the value
!C     of NSPEC agrees with the number of species in the file and produce
!C     an error message if not.
!C-----------------------------------------------------------------------
!      SUBROUTINE READ_SPECIES(NSPEC,SPECIES,ABUNDANCE,MASS)
      SUBROUTINE READ_SPECIES(NSPEC,ABUNDANCE, nco, nc, ncx, no, nelect)

!T.Bell
!use definitions
!
!use healpix_types
!use healpix_guts
!use pdr_mod
use unix_mod, only: unixGetenv
!use maincode_module, only : nspec, abundance, mass, species

      IMPLICIT NONE
!      INTEGER(kind=i4b),intent(in) :: NSPEC
      INTEGER,intent(in) :: NSPEC
!      INTEGER(kind=i4b) :: NH,ND,NH2,NHD,NC,NCx,NCO,NO,NELECT,NPROTON,NH2O,NHe, &
!     &        NMG,NMGx,NN,NFE,NFEx,NSI,NSIx,NCA,NCAx,NCAxx,NS,NSx,NCS, &
!     &        NOSH,NCL,NCLx,NH2x,NHEx,NOx,NNx,NNA,NNAx,NCH,NCH2,NOH,NO2
      real(kind=dp), intent(out) :: ABUNDANCE(1:nspec)
      real(kind=dp) :: MASS(1:nspec) !intent(out)
      CHARACTER(len=10) :: SPECIES(1:nspec) !intent(out)
      character(len=400) :: speciesData
      character(len=200) :: dataDirectory, dataPath
      INTEGER(kind=i4b) :: I,INDEX,SPECIESFILE

!      INTEGER(kind=i4b) :: NH,ND,NH2,NHD,NC,NCx,NCO,NO,NELECT,NPROTON,NH2O,NHe, &
!           &        NMG,NMGx,NN,NFE,NFEx,NSI,NSIx,NCA,NCAx,NCAxx,NS,NSx,NCS, &
!           &        NCL,NCLx,NH2x,NHEx,NOx,NNx,NNA,NNAx,NCH,NCH2,NOH,NO2, &
!           &        NH3x, NH3Ox, NHCOx!, NOSH
      integer :: NC, NCx, NCO, NO, NELECT

!   allocate(species(1:nspec))
!   allocate(abundance(1:nspec))
!   allocate(mass(1:nspec))

      SPECIESFILE = 3

!C     Initialize the variables and read in the species data. Check that
!C     the value of NSPEC agrees with the number of species in the file
!C     and produce an error message if not.
      SPECIES="          "
      ABUNDANCE=0.0D0
      MASS=0.0D0

!C     Initialize all the species index labels. If they are not assigned
!C     subsequently, any attempt to access that species will generate an
!C     error and the code will crash. This is a useful bug catch.
      NH=0
      ND=0
      NH2=0
      NHD=0
      NH2x=0
      NPROTON=0
      NC=0
      NCx=0
      NO=0
      NOx=0
      NN=0
      NNx=0
      NS=0
      NSx=0
      NHE=0
      NHEx=0
      NNA=0
      NNAx=0
      NMG=0
      NMGx=0
      NSI=0
      NSIx=0
      NFE=0
      NFEx=0
      NCL=0
      NCLx=0
      NCA=0
      NCAx=0
      NCAxx=0
      NCO=0
      NCH=0
      NCH2=0
      NOH=0
      NO2=0
      NCS=0
      NH2O=0
      NELECT=0
      NH3x=0
      NH3Ox=0
      NHCOx=0
!#ifdef REDUCED
!#ifdef XRAYS
!      OPEN(SPECIESFILE,FILE="species_Xreduced.d",STATUS="OLD")
!#else

      call unixGetenv("TORUS_DATA", dataPath, i)
      
      dataDirectory = trim(dataPath)
      
      speciesData = trim(dataDirectory) // "/pdrdata/species_reduced.d"

      
      OPEN(SPECIESFILE,FILE=speciesData,STATUS="OLD")
!#endif
!#endif
!#ifdef FULL
!#ifdef XRAYS
!      OPEN(SPECIESFILE,FILE="species_Xfull.d",STATUS="OLD")
!#else
!      OPEN(SPECIESFILE,FILE="species_full.d",STATUS="OLD")
!#endif
!#endif
!#ifdef MYNETWORK
!#ifdef XRAYS
!      OPEN(SPECIESFILE,FILE="species_Xmynetwork.d",STATUS="OLD")
!#else
!      OPEN(SPECIESFILE,FILE="species_mynetwork.d",STATUS="OLD")
!#endif
!#endif
      REWIND(SPECIESFILE)
      DO I=1,NSPEC
         READ(SPECIESFILE,*,END=1) INDEX,SPECIES(I),ABUNDANCE(I),MASS(I)

!C        Assign the various index labels to their correct species.
         IF(SPECIES(I).EQ."H         ") then 
            print *, "NH IS", I
            NH      = I
         endif
         IF(SPECIES(I).EQ."D         ") ND      = I
         IF(SPECIES(I).EQ."H2        ") then
            NH2     = I
            print *, "NH2 IS ", NH2
         endif
         IF(SPECIES(I).EQ."HD        ") then
            print *, "NHD IS ", NHD
            NHD     = I
         endif
         IF(SPECIES(I).EQ."H2+       ") NH2x    = I
         IF(SPECIES(I).EQ."H3+       ") NH3x    = I
         IF(SPECIES(I).EQ."H+        ") then
            NPROTON = I
            print *, "nproton is", nproton
         endif
         IF(SPECIES(I).EQ."C         ") then
            print *, "NC is", I
            NC      = I
         endif
         IF(SPECIES(I).EQ."C+        ") then
            print *, "NCII is ", I
            NCx     = I
         endif
         IF(SPECIES(I).EQ."O         ") NO      = I
         IF(SPECIES(I).EQ."O+        ") NOx     = I
         IF(SPECIES(I).EQ."N         ") NN      = I
         IF(SPECIES(I).EQ."N+        ") NNx     = I
         IF(SPECIES(I).EQ."S         ") then
            NS      = I
            print *, "NS IS", I
         endif
         IF(SPECIES(I).EQ."S+        ") NSx     = I
         IF(SPECIES(I).EQ."He        ") then
            print *, "NHE IS", I
            NHE     = I
         endif
         IF(SPECIES(I).EQ."HE        ") NHE     = I
         IF(SPECIES(I).EQ."He+       ") then
            print *, "NHEx is ", I
            NHEx    = I
         endif
         IF(SPECIES(I).EQ."HE+       ") NHEx    = I
         IF(SPECIES(I).EQ."Na        ") NNA     = I
         IF(SPECIES(I).EQ."NA        ") NNA     = I
         IF(SPECIES(I).EQ."Na+       ") NNAx    = I
         IF(SPECIES(I).EQ."NA+       ") NNAx    = I
         IF(SPECIES(I).EQ."Mg        ") NMG     = I
         IF(SPECIES(I).EQ."MG        ") NMG     = I
         IF(SPECIES(I).EQ."Mg+       ") NMGx    = I
         IF(SPECIES(I).EQ."MG+       ") NMGx    = I
         IF(SPECIES(I).EQ."Si        ") NSI     = I
         IF(SPECIES(I).EQ."SI        ") NSI     = I
         IF(SPECIES(I).EQ."Si+       ") NSIx    = I
         IF(SPECIES(I).EQ."SI+       ") NSIx    = I
         IF(SPECIES(I).EQ."Fe        ") NFE     = I
         IF(SPECIES(I).EQ."FE        ") NFE     = I
         IF(SPECIES(I).EQ."Fe+       ") NFEx    = I
         IF(SPECIES(I).EQ."FE+       ") NFEx    = I
         IF(SPECIES(I).EQ."Cl        ") NCL     = I
         IF(SPECIES(I).EQ."CL        ") NCL     = I
         IF(SPECIES(I).EQ."Cl+       ") NCLx    = I
         IF(SPECIES(I).EQ."CL+       ") NCLx    = I
         IF(SPECIES(I).EQ."Ca        ") NCA     = I
         IF(SPECIES(I).EQ."CA        ") NCA     = I
         IF(SPECIES(I).EQ."Ca+       ") NCAx    = I
         IF(SPECIES(I).EQ."CA+       ") NCAx    = I
         IF(SPECIES(I).EQ."Ca++      ") NCAxx   = I
         IF(SPECIES(I).EQ."CA++      ") NCAxx   = I
         IF(SPECIES(I).EQ."CO        ") then
            print *, "NCO IS", I
            NCO     = I
         endif
         IF(SPECIES(I).EQ."CH        ") NCH     = I
         IF(SPECIES(I).EQ."CH2       ") NCH2    = I
         IF(SPECIES(I).EQ."OH        ") NOH     = I
         IF(SPECIES(I).EQ."O2        ") NO2     = I
         IF(SPECIES(I).EQ."CS        ") NCS     = I
         IF(SPECIES(I).EQ."H2O       ") NH2O    = I
         IF(SPECIES(I).EQ."H3O+      ") NH3Ox   = I
         IF(SPECIES(I).EQ."HCO+      ") NHCOx   = I
         IF(SPECIES(I).EQ."e-        ") then
            print *, "NELECT IS", I
            NELECT  = I
         endif
         IF(SPECIES(I).EQ."ELECTR    ") then
            print *, "NELECT IS B", I
            NELECT  = I
         endif
      ENDDO

      I=I-1
      READ(SPECIESFILE,*,END=1)
      I=I+1
 1    IF(I.NE.NSPEC .and. myrankglobal == 1) THEN
         write(6,*) 'ERROR! Number of species (NSPEC) does not match ',&
     &           'the number of entries in the species file'
         STOP
      ENDIF

!C     Check that the final species in the file is e-. Print a warning
!C     message to screen and logfile if not.
      IF(SPECIES(NSPEC).NE."e-") THEN
         write(6,*) 'WARNING! Last entry in species file is not e-'
         WRITE(10,*)'WARNING! Last entry in species file is not e-'
      ENDIF


!C     Check that the total hydrogen nuclei abundance adds up to 1.
!C     If not, modify the abundance of H2 (only consider H, H+ & H2)
      IF((ABUNDANCE(NH)+ABUNDANCE(NPROTON)+2.0D0*ABUNDANCE(NH2)).NE.1.0D0) THEN
         ABUNDANCE(NH2)=0.5D0*(1.0D0-ABUNDANCE(NH)-ABUNDANCE(NPROTON))
      ENDIF

!C     Calculate the intial electron abundance, if not 
!C     specified, as the sum of the metal ion abundances
      IF(ABUNDANCE(NELECT).LE.0.0D0) THEN
       ABUNDANCE(NELECT)=0.0D0
       IF(NCx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NCx)
       IF(NSx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NSx)
       IF(NNAx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NNAx)
       IF(NMGx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NMGx)
       IF(NSIx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NSIx)
       IF(NFEx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NFEx)
       IF(NCLx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NCLx)
       IF(NCAx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NCAx)
       IF(NCAxx.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+2.0D0*ABUNDANCE(NCAxx)
       IF(NPROTON.NE.0) ABUNDANCE(NELECT)=ABUNDANCE(NELECT)+ABUNDANCE(NPROTON)
      ENDIF

      CLOSE(SPECIESFILE)
      RETURN
      END SUBROUTINE



      SUBROUTINE READINPUT(FILENAME,NLEV,NTEMP,ENERGIES,WEIGHTS,&
           A_COEFFS,B_COEFFS,FREQUENCIES,TEMPERATURES,&
           H_COL,HP_COL,EL_COL,HE_COL,H2_COL,PH2_COL,OH2_COL)
        
        !T.Bell
        
!        use healpix_types
!        use healpix_guts
        implicit none
        INTEGER(kind=I4B), INTENT(IN)::NLEV
        INTEGER(kind=I4B), INTENT(IN)::NTEMP
        CHARACTER(len=*),INTENT(IN)::FILENAME
        
        real(kind=dp), intent(out)::ENERGIES(1:NLEV), WEIGHTS(1:NLEV)
        real(kind=dp), intent(out)::A_COEFFS(1:NLEV,1:NLEV), B_COEFFS(1:NLEV,1:NLEV)
        real(kind=dp), intent(out)::FREQUENCIES(1:NLEV,1:NLEV), TEMPERATURES(1:7,1:NTEMP)
        real(kind=dp), intent(out)::H_COL(1:NLEV,1:NLEV,1:NTEMP)
        real(kind=dp), intent(out)::HP_COL(1:NLEV,1:NLEV,1:NTEMP)
        real(kind=dp), intent(out)::EL_COL(1:NLEV,1:NLEV,1:NTEMP)
        real(kind=dp), intent(out)::HE_COL(1:NLEV,1:NLEV,1:NTEMP)
        real(kind=dp), intent(out)::H2_COL(1:NLEV,1:NLEV,1:NTEMP)
        real(kind=dp), intent(out)::PH2_COL(1:NLEV,1:NLEV,1:NTEMP)
        real(kind=dp), intent(out)::OH2_COL(1:NLEV,1:NLEV,1:NTEMP)
        
        INTEGER(kind=I4B)::NLIN,NPART,NCOL
        INTEGER(kind=I4B)::I,J,K,L,M,N,P
        real(kind=dp):: ENERGY,WEIGHT,EINSTEINA,FREQUENCY
        real(kind=dp)::COEFF(1:NTEMP)
        
        
        !     Initialize all variables to zero before reading in the data
        DO I=1,NLEV
           ENERGIES(I)=0.0D0
           WEIGHTS(I)=0.0D0
           DO J=1,NLEV
              A_COEFFS(I,J)=0.0D0
              B_COEFFS(I,J)=0.0D0
              FREQUENCIES(I,J)=0.0D0
              DO K=1,NTEMP
                 TEMPERATURES(:,K)=0.0D0
                 H_COL(I,J,K)=0.0D0
                 EL_COL(I,J,K)=0.0D0
               HE_COL(I,J,K)=0.0D0
               H2_COL(I,J,K)=0.0D0
               PH2_COL(I,J,K)=0.0D0
               OH2_COL(I,J,K)=0.0D0
            ENDDO
         ENDDO
      ENDDO

      OPEN(8,FILE=FILENAME,STATUS='OLD')
      READ(8,'(////)') !empty line
      READ(8,*) N !number of levels of the file
      IF(N.NE.NLEV) STOP "ERROR! Incorrect number of energy levels, N>NLEV"
      READ(8,*)
      N=0
      DO WHILE(N.LT.NLEV)
         READ(8,*) N,ENERGY,WEIGHT
         ENERGIES(N)=ENERGY*C*HP ! Convert from cm^-1 to erg
         WEIGHTS(N)=WEIGHT
      ENDDO
      READ(8,*)
      READ(8,*) NLIN
      READ(8,*)
      N=0
      DO WHILE(N.LT.NLIN)
         READ(8,*) N,I,J,EINSTEINA,FREQUENCY
         FREQUENCIES(I,J)=FREQUENCY*1.0D9 ! Convert from GHz to Hz
         FREQUENCIES(J,I)=FREQUENCIES(I,J)
         A_COEFFS(I,J)=EINSTEINA
!        Calculate the Einstein B coefficients using Bij = Aij/(2.h.nu^3/c^2)
         B_COEFFS(I,J)=A_COEFFS(I,J)&
     &                 /(2.0D0*HP*(FREQUENCIES(I,J)**3)/(C**2))
         B_COEFFS(J,I)=B_COEFFS(I,J)*(WEIGHTS(I)/WEIGHTS(J))
      ENDDO
!     Calculate the transition frequencies between all levels (even if forbidden)
      DO I=1,NLEV
         DO J=1,NLEV
            FREQUENCY=ABS(ENERGIES(I)-ENERGIES(J))/HP
            IF(FREQUENCIES(I,J).NE.0.0D0) THEN
!              Check if the calculated and measured frequencies differ by >1%
               IF(ABS(FREQUENCY-FREQUENCIES(I,J))&
     &                /FREQUENCIES(I,J).GT.1.0D-2) THEN
                  WRITE(6,*) 'ERROR! Calculated frequency differs by >1%:'
                  WRITE(6,*) FREQUENCY,' Hz vs',FREQUENCIES(I,J),' Hz'
                  STOP
               ENDIF
            ELSE
               FREQUENCIES(I,J)=FREQUENCY
            ENDIF
         ENDDO
      ENDDO
      READ(8,*)
      READ(8,*) NPART
!     Read the collisional rate coefficients (cm^3 s^-1) for each collision partner
      DO L=1,NPART
         READ(8,*)
         READ(8,*) P
         READ(8,*)
         READ(8,*) NCOL
         READ(8,*)
         READ(8,*) M
         IF(M.GT.NTEMP) THEN
            WRITE(6,*) 'ERROR! Too many temperature values (>NTEMP):',M
            STOP
         ENDIF
         IF(P.EQ.1) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  H2_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(H2_COL(I,J,K).NE.0.0D0 .AND. H2_COL(J,I,K).EQ.0.0D0) THEN
                     H2_COL(J,I,K)=H2_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.2) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  PH2_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(PH2_COL(I,J,K).NE.0.0D0 .AND. PH2_COL(J,I,K).EQ.0.0D0) THEN
                     PH2_COL(J,I,K)=PH2_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
!         ELSE IF(P.EQ.3) THEN
!            READ(8,*)
!            READ(8,*) (TEMPERATURES(P,K),K=1,M)
!            READ(8,*)
!            N=0
!            DO WHILE(N.LT.NCOL)
!               READ(8,*) N,I,J,(COEFF(K),K=1,M)
!               DO K=1,M
!                  OH2_COL(I,J,K)=COEFF(K)   
!!                 Calculate the reverse (excitation) rate coefficient
!!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
!                  IF(PH2_COL(I,J,K).NE.0.0D0 .AND. PH2_COL(J,I,K).EQ.0.0D0) THEN
!                     PH2_COL(J,I,K)=PH2_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
!     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
!                  ENDIF
!               ENDDO
!            ENDDO
         ELSE IF(P.EQ.3) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  OH2_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(OH2_COL(I,J,K).NE.0.0D0 .AND. OH2_COL(J,I,K).EQ.0.0D0) THEN
                     OH2_COL(J,I,K)=OH2_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.4) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  EL_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(EL_COL(I,J,K).NE.0.0D0 .AND. EL_COL(J,I,K).EQ.0.0D0) THEN
                     EL_COL(J,I,K)=EL_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.5) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  H_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(H_COL(I,J,K).NE.0.0D0 .AND. H_COL(J,I,K).EQ.0.0D0) THEN
                     H_COL(J,I,K)=H_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.6) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  HE_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(HE_COL(I,J,K).NE.0.0D0 .AND. HE_COL(J,I,K).EQ.0.0D0) THEN
                     HE_COL(J,I,K)=HE_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF(P.EQ.7) THEN
            READ(8,*)
            READ(8,*) (TEMPERATURES(P,K),K=1,M)
            READ(8,*)
            N=0
            DO WHILE(N.LT.NCOL)
               READ(8,*) N,I,J,(COEFF(K),K=1,M)
               DO K=1,M
                  HP_COL(I,J,K)=COEFF(K)
!                 Calculate the reverse (excitation) rate coefficient
!                 from detailed balance: Cji = Cij*gi/gj*exp(-(Ei-Ej)/kT)
                  IF(HP_COL(I,J,K).NE.0.0D0 .AND. HP_COL(J,I,K).EQ.0.0D0) THEN
                     HP_COL(J,I,K)=HP_COL(I,J,K)*(WEIGHTS(I)/WEIGHTS(J)) &
     &              *EXP(-(ENERGIES(I)-ENERGIES(J))/(KB*TEMPERATURES(P,K)))
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            WRITE(6,*) 'ERROR! Unrecognized collision partner ID:',P
            STOP
         ENDIF
      ENDDO

      CLOSE(8)
      if(myrankglobal == 1) then
         WRITE(6,*) 'Cooling datafile: ',FILENAME,' read successfully'
      end if
      RETURN

      END subroutine

#endif
end module setuppdr_mod






