#ifdef PDR
module pdr_utils_mod
!!
! use definitions
! use healpix_types
use healpix_guts
use constants_mod
use messages_mod
use parallel_mod
use gridio_mod
use source_mod
use timing
use grid_mod
use vtk_mod
use amr_mod
use mpi_amr_mod
use mpi_global_mod
use utils_mod
 implicit none


!UCLPDR TOM BELL



  logical :: start
  integer, save :: DIMH2=6
  integer, save :: DIMCO=8

  real(double), dimension(8), save :: NCO_GRID=(/&
    &12.0D0,13.0D0,14.0D0,15.0D0, 16.0D0,17.0D0,18.0D0,19.0D0/)
  real(double), dimension(6), save :: NH2_GRID=(/&
    &18.0D0,19.0D0,20.0D0,21.0D0,22.0D0,23.0D0/)
  real(double) :: SCO_GRID(1:8,1:6)
  integer, save :: N_GRID=30



  real(double) :: SH2_DERIV(1:105), SCO_DERIV(1:8,1:6), X_DERIV(1:30)


! interface
!    module procedure function H2PDRATE(K0,G0,AV,NH2)
!      use definitions	
!      use healpix_types		
!      real(double) :: H2PDRATE
!      real(double), intent(in) :: k0, g0, av
!      real(double), intent(in) :: nh2
!      real(double) :: lambda, scatter, h2shield2
!    end function H2PDRATE!
!
!!   function COPDRATE(K0,G0,AV,NCO,NH2)
 !    use definitions
!     use healpix_types
!     real(double) :: copdrate
!     real(double), intent(in) :: k0, g0, av, nco
!     real(double), intent(in) :: nh2
!!     real(double) :: lambda, lbar, coshield, scatter
!   end function copdrate
!
!   function CIPDRATE(K0,G0,AV,KAV,NCI,NH2,TGAS)
!     use definitions
!     use healpix_types
!     real(double) :: cipdrate
!     real(double), intent(in) :: K0,G0,AV,KAV,NCI,TGAS
!     real(double), intent(in) :: nh2
!     real(double) :: tauc
!   end function cipdrate
!
!   function SIPDRATE(K0,G0,AV,KAV,NSI)
!     use definitions
!     use healpix_types
!     real(double) :: sipdrate
!     real(double), intent(in) :: K0,G0,AV,KAV,NSI
!     real(double) :: taus
!   end function sipdrate
!
!   function H2SHIELD1(NH2,DOPW,RADW)
!     use definitions
!     use healpix_types
!     real(double) :: h2shield1
!     real(double), intent(in) :: nh2
!     real(double), intent(in) ::DOPW,RADW
!     real(double) :: FPARA, FOSC, TAUD, R, T, U, JD, JR
!   end function h2shield1
!
!   function h2shield2(nh2)
!     use definitions
!     use healpix_types
!!     use uclpdr_module, only : start, numh2, COL_GRID, SH2_GRID, SH2_DERIV
!     real(double) :: h2shield2
!     real(double), intent(in) :: nh2
!   end function h2shield2
!!
 !  function COSHIELD(NCO,NH2)  
 !    use definitions
 !    use healpix_types
!!     use uclpdr_module, only : start, NCO_GRID, NH2_GRID, SCO_GRID, SCO_DERIV
 !    real(double) :: COSHIELD
 !!    real(double) :: LOGNCO, LOGNH2
  !   real(double), intent(in) :: NCO, NH2
!   end function COSHIELD

!   function SCATTER(AV,LAMBDA)
!     use definitions
!     use healpix_types
!     real(double) :: scatter
!     real(double), intent(in) :: AV, LAMBDA
!     real(double), dimension(0:5), save :: A = (/&
!           &1.000D0,2.006D0,-1.438D0,0.7364D0,-0.5076D0,-0.0592D0/)
!     real(double), dimension(0:5), save :: K = (/&
!           &0.7514D0,0.8490D0,1.013D0,1.282D0,2.005D0,5.832D0/)
!     real(double) :: EXPONENT, XLAMBDA
!   end function scatter
!
!   function XLAMBDA(LAMBDA)
!     use definitions
!     use healpix_types
!!     use uclpdr_module, only : start, N_GRID, L_GRID, X_GRID, X_DERIV
!     real(double) :: xlambda
!     real(double), intent(in) :: lambda
!   end function xlambda!

!   function LBAR(NCO,NH2)
!     use definitions
!     use healpix_types
!     real(double) :: lbar
!!     real(double) :: U,W 
 !    real(double) :: NCO, NH2
 !  end function LBAR!!
!
!    function calculate_heating(density, gas_temperature, dust_temperature, UV_field, &
!             & v_turb, nspec, dummyabundance, nreac, rate)
!      use definitions
!      use healpix_types
!!      real(double) :: calculate_heating
 !     integer :: nspec, nreac
 !     real(double) :: density, gas_temperature, dust_temperature, UV_field, v_turb
 !     real(double) :: dummyabundance(1:nspec), rate(1:nreac)
 !   end function calculate_heating

!#ifdef H2FORM
!!
!      FUNCTION H2_FORMATION_RATE(GAS_TEMPERATURE,GRAIN_TEMPERATURE) RESULT(RATE)
!         USE DEFINITIONS
!         USE HEALPIX_TYPES
!         IMPLICIT NONE
!         REAL(double) :: RATE
!!         REAL(double), INTENT(IN) :: GAS_TEMPERATURE,GRAIN_TEMPERATURE
 !     END FUNCTION H2_FORMATION_RATE
!
!#endif
!
!!!
!
!!
! end interface

contains




integer function nray_func()
  use inputs_mod, only : hlevel
  integer :: nside
  nside = 2**hlevel
  nray_func = 12*nside**2
end function nray_func


!C=======================================================================
!C
!C     !Calculate the rate coefficient for each reaction at the specified
!C     temperature and visual extinction, Av. The photodissociation of H2
!C     and !CO and the photoionization of !CI and SI are treated separately
!C     in more detail (see the functions in shield.f). Multiple rates for
!C     the same reaction (duplicates) are allowed in the ratefile and are
!C     activated based on their minimum and maximum temperature specified
!C     in the file. Negative gamma factors are ignored below the minimum
!C     temperature at which the reaction rate is valid.
!C
!C-----------------------------------------------------------------------
!#ifdef XRAYS
!      SUBROUTINE CALCULATE_REACTION_RATES(TEMPERATURE,DUST_TEMPERATURE,NRAYS,RAD_SURFACE,XRAD_SURFACE,&
!                 &AV,COLUMN,NREAC, REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RATE,RTMIN,RTMAX,DUPLICATE,NSPEC,&
!                 &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
!#else
      SUBROUTINE CALCULATE_REACTION_RATES(TEMPERATURE,DUST_TEMPERATURE,RAD_SURFACE,AV,COLUMN, &
                 &REACTANT,PRODUCT,ALPHA,BETA,GAMMA,RATE,RTMIN,RTMAX,DUPLICATE,NSPEC,&
                 &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI, nreac, nrays, n12co, nci)
!#endif
 
!T.Bell

!        use definitions
!        use healpix_types
! use global_module
! use functions_module
! use maincode_module , only : zeta, AV_fac
      IMPLICIT NONE

      integer, intent(in) ::  NSPEC, nreac, nrays
      real(double),intent(in) :: TEMPERATURE, DUST_TEMPERATURE
      real(double),intent(in) :: RAD_SURFACE(1:nrays),AV(1:nrays),COLUMN(1:nrays,1:nspec)
!#ifdef XRAYS
!      real(double),intent(in) :: XRAD_SURFACE(0:nrays-1)
!#endif
      integer,intent(in) :: DUPLICATE(1:nreac)
      real(double),intent(in) :: ALPHA(1:nreac),BETA(1:nreac),GAMMA(1:nreac),RTMIN(1:nreac),RTMAX(1:nreac)
      real(double), intent(out) :: rate(1:nreac)
      integer,intent(out):: NRGR,NRH2,NRHD,NRCO,NRCI,NRSI
      CHARACTER(len=10),intent(in) ::  REACTANT(1:nreac,1:3),PRODUCT(1:nreac,1:4)

!      integer ::  NRGR,NRH2,NRHD,NRCO,NRCI,NRSI
!      real(double) :: H2PDRATE,COPDRATE,CIPDRATE,SIPDRATE
      real(double) :: CION,STICKING,FLUX,YIELD
      integer :: I,J,K!, nreac, nrays

      REAL(double), parameter :: OMEGA=0.42D0,GRAIN_RADIUS=1.0D-7,METALLICITY=1.0D0
!      real(double) :: NH2, NHD, NCO, NC, NS
      integer, intent(in) :: n12co,  nci
      integer ::NH2, NCO, NC, NS
      integer :: NHD
!#ifdef MOCASSIN
!      real(double) :: RADSURFTOT
!#endif
!C     Initialize the rate coefficients.
      
      NCO = n12co
      NC = nci
      NH2 = 31
      NS = 30

      RATE=0.0D0

!C     Initialize the stored reaction numbers. If they are not assigned
!C     subsequently, any attempt to access that reaction will generate an
!C     error and the code will crash. This is a useful bug catch.
      NRGR=0
      NRH2=0
      NRHD=0
      NRCO=0
      NRCI=0
      NRSI=0

!THIS NEEDS TO BE CHANGED SO THAT THE REACTANT AND PRODUCT ARE TRIMMED
!328 and 329, 0, 290, 269, 0 for NR...

      !print *, "NREAC ", NREAC
      DO I=1,NREAC
!         print *, "I is ", I
!C        Determine the type of reaction
     !    print *, len(REACTANT(I, 2))
     !    print *, "2 ", REACTANT(I, 2), "FRWFWE"
         IF(trim(REACTANT(I,2)).EQ."PHOTON") GOTO 1
         IF(trim(REACTANT(I,2)).EQ."CRP") GOTO 2
         IF(trim(REACTANT(I,2)).EQ."CRPHOT") GOTO 3
         IF(trim(REACTANT(I,2)).EQ."FREEZE") GOTO 4
         IF(trim(REACTANT(I,2)).EQ."ELFRZE") GOTO 5
         IF(trim(REACTANT(I,2)).EQ."CRH") GOTO 6
         IF(trim(REACTANT(I,2)).EQ."PHOTD") GOTO 7
         IF(trim(REACTANT(I,2)).EQ."THERM") GOTO 8
         IF(trim(REACTANT(I,2)(1:1)).EQ."G") GOTO 9
!#ifdef XRAYS
!!ADD IF REACTANT IS XRAY GOTO 11
!!         IF(REACTANT(I,2).EQ."XRAY  ") GOTO 11
!         IF(REACTANT(I,2).EQ."XRAY  ") GOTO 11
!         IF(REACTANT(I,2).EQ."XRSEC ") GOTO 11
!         IF(REACTANT(I,2).EQ."XRLYA ") GOTO 11
!         IF(REACTANT(I,2).EQ."XRPHOT") GOTO 11
!#else
         IF(trim(REACTANT(I,2)).EQ."XRAY") GOTO 110
         IF(trim(REACTANT(I,2)).EQ."XRSEC") GOTO 120
         IF(trim(REACTANT(I,2)).EQ."XRLYA") GOTO 130
         IF(trim(REACTANT(I,2)).EQ."XRPHOT") GOTO 140
!#endif

!C-----------------------------------------------------------------------

!C     Thermal reactions:

!C     Store the reaction number for H2 grain formation. The rate of H2
!C     formation on grains is calculated directly, using the formula of
!C     de Jong, 1977, A&A, 55, 137 (page 140, right column). The second
!C     dependence on temperature helps to prevent runaway H2 formation
!C     heating at high temperatures:
!C
!C        k_H2 = 3E-18 * T^0.5 * exp(-T/1000)   (cm3/s)
!C
!         print *, "1 ", REACTANT(I, 1), " TRWFE"
         IF(trim(REACTANT(I,1)).EQ."H  " .AND. trim(REACTANT(I,2)).EQ."H  "  .AND. &
         & (trim(REACTANT(I,3)).EQ."   " .OR.  trim(REACTANT(I,3)).EQ."#  ") .AND. &
         &  trim(PRODUCT(I,1)).EQ."H2"  .AND. &
         & (trim(PRODUCT(I,2)).EQ.""  .OR.  trim(PRODUCT(I,2)).EQ."#")) THEN

!old way to read reactants, omitting the #
!         IF(REACTANT(I,1).EQ."H  " .AND. REACTANT(I,2).EQ."H  " .AND. &
!           REACTANT(I,3).EQ."   " .AND. PRODUCT(I,1).EQ."H2 " .AND. &
!       &       PRODUCT(I,2).EQ."   ") THEN

!#ifdef H2FORM
!            print *, "PART 2"
          RATE(I)=H2_FORMATION_RATE(TEMPERATURE,DUST_TEMPERATURE)
!#else
!!            T.Bell,15/11/11 - Use the H2 grain formation rate from the
!!            Rollig et al. (2007) benchmark paper (no exponential term)
!            RATE(I)=3.0D-18*SQRT(TEMPERATURE)*EXP(-(TEMPERATURE/1.0D3))
!!            RATE(I)=3.0D-18*SQRT(TEMPERATURE)
!#endif
            NRGR=I
            GOTO 10
         ENDIF

!C     Check for large negative gamma values that might cause discrepant
!C     rates at low temperatures. Set these rates to zero when T < RTMIN.
!         print *, "PART 1"
         IF(DUPLICATE(I).EQ.0) THEN
            IF(GAMMA(I).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(I)) THEN
               RATE(I)=0.0D0
            ELSE
               RATE(I)=ALPHA(I)*(TEMPERATURE/300.0D0)**BETA(I)*EXP(-(GAMMA(I)/TEMPERATURE))
            ENDIF
         ELSE IF(DUPLICATE(I).EQ.1) THEN
            J=I
            DO
               IF(TEMPERATURE.LE.RTMAX(J)) THEN
                  IF(GAMMA(J).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(J)) THEN
                  RATE(J)=0.0D0
                  ELSE
                  RATE(J)=ALPHA(J)*(TEMPERATURE/300.0D0)**BETA(J)*EXP(-(GAMMA(J)/TEMPERATURE))
                  ENDIF
                  EXIT
               ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
                  IF(GAMMA(J).LT.-200.0D0 .AND. TEMPERATURE.LT.RTMIN(J)) THEN
                  RATE(J)=0.0D0
                  ELSE
                  RATE(J)=ALPHA(J)*(TEMPERATURE/300.0D0)**BETA(J)*EXP(-(GAMMA(J)/TEMPERATURE))
                  ENDIF
                  EXIT
               ELSE
                  RATE(J)=0.0D0
                  J=J+1
               ENDIF
            ENDDO
         ENDIF
         GOTO 10

!C-----------------------------------------------------------------------

!C     Photoreactions:

!C     Store the reaction number for H2 photodissociation. The rate itself
!C     is calculated separately by the function H2PDRATE (within shield.f)

 1       IF(trim(REACTANT(I,1)).EQ."H2" .AND. trim(REACTANT(I,3)).EQ."" .AND.&
     &      trim(PRODUCT(I,1)).EQ."H" .AND. trim(PRODUCT(I,2)).EQ."H") THEN
!C           Loop over all rays
!#ifdef MOCASSIN
!            RADSURFTOT = 0.0D0
!            DO K=0,NRAYS-1
!!               RATE(I)=RATE(I) + H2PDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NH2))
 !              RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
 !           ENDDO
!            RATE(I) = RATE(I)/RADSURFTOT
!#else
 !           if(myrankglobal == 1) then
!               print *, "DOING ", NRH2, "FOR ", NRAYS, "REAC ", NREAC
  !          end if
!            print *, "PART 3"
            DO K=1,NRAYS
!            DO K=0,NRAYS-1
               RATE(I)=RATE(I) + H2PDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NH2))
            ENDDO
!#endif
            NRH2=I
            GOTO 10
         ENDIF


!THAW NHD IS NOT INCLUDED IN THE SPECIES DATAFILE

!C     Store the reaction number for HD photodissociation. The rate itself
!C     is calculated separately by the function H2PDRATE (within shield.f)
!         IF(trim(REACTANT(I,1)).EQ."HD" .AND. trim(REACTANT(I,3)).EQ."" .AND.&
!     &    ((trim(PRODUCT(I,1)).EQ."H" .AND. trim(PRODUCT(I,2)).EQ."D") .OR.&
!     &     (trim(PRODUCT(I,1)).EQ."D" .AND. trim(PRODUCT(I,2)).EQ."H"))) THEN
!C           Loop over all rays
!#ifdef MOCASSIN
!            RADSURFTOT = 0.0D0
!            DO K=0,NRAYS-1
!               RATE(I)=RATE(I) + H2PDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NHD))
!               RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
!            ENDDO
!            RATE(I) = RATE(I)/RADSURFTOT
!#else
!            print *, "PART 4"
!            DO K=1,NRAYS
!!            DO K=0,NRAYS-1
!               RATE(I)=RATE(I) + H2PDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NHD))
!            ENDDO
!#endif
!            NRHD=I
!            GOTO 10
!         ENDIF

!C     Store the reaction number for !CO photodissociation. The rate itself
!C     is calculated separately by the function !COPDRATE (within shield.f)
         IF(trim(REACTANT(I,1)).EQ."CO" .AND. trim(REACTANT(I,3)).EQ."" .AND.&
     &    ((trim(PRODUCT(I,1)).EQ."C" .AND. trim(PRODUCT(I,2)).EQ."O") .OR.&
     &     (trim(PRODUCT(I,1)).EQ."O" .AND. trim(PRODUCT(I,2)).EQ."C"))) THEN
!C           Loop over all rays
!#ifdef MOCASSIN
!            RADSURFTOT = 0.0D0
!            DO K=0,NRAYS-1
!               RATE(I)=RATE(I) + COPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NCO),COLUMN(K,NH2))
!!               RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
 !           ENDDO
 !           RATE(I) = RATE(I)/RADSURFTOT
!#else
!            print *, "PART 5"
!            DO K=0,NRAYS-1
            DO K=1,NRAYS
!!               print *, "COMMENCING DUMP"
!               print *, RATE(I)
!               print *, ALPHA(I)
!               print *, RAD_SURFACE(K)
!               print *, AV(K) 
!               print *, COLUMN(K,NCO)
!               print *, COLUMN(K,NH2)
!               Print *, COPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NCO),COLUMN(K,NH2))
               RATE(I)=RATE(I) + COPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),COLUMN(K,NCO),COLUMN(K,NH2))
            ENDDO
!#endif
            NRCO=I
            GOTO 10
         ENDIF

!C     Store the reaction number for !CI photoionization. The rate itself
!C     is calculated separately by the function CIPDRATE (within shield.f)
         IF(trim(REACTANT(I,1)).EQ."C" .AND. trim(REACTANT(I,3)).EQ."" .AND.&
     &    ((trim(PRODUCT(I,1)).EQ."C+" .AND. trim(PRODUCT(I,2)).EQ."e-") .OR.&
     &     (trim(PRODUCT(I,1)).EQ."e-" .AND. trim(PRODUCT(I,2)).EQ."C+"))) THEN
!C           Loop over all rays
!#ifdef MOCASSIN
!            RADSURFTOT = 0.0D0
!            DO K=0,NRAYS-1
!               RATE(I)=RATE(I) + CIPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),GAMMA(I),COLUMN(K,NC),COLUMN(K,NH2),TEMPERATURE)
!               RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
!            ENDDO
!            RATE(I) = RATE(I)/RADSURFTOT
!#else
!            print *, "PART 6"
!            DO K=0,NRAYS-1
            DO K=1,NRAYS
               RATE(I)=RATE(I) + CIPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),GAMMA(I),COLUMN(K,NC),COLUMN(K,NH2),TEMPERATURE)
            ENDDO
!#endif
            NRCI=I
            GOTO 10
         ENDIF

!C     Store the reaction number for SI photoionization. The rate itself
!C     is calculated separately by the function SIPDRATE (within shield.f)
         IF(trim(REACTANT(I,1)).EQ."S" .AND. trim(REACTANT(I,3)).EQ."" .AND.&
     &    ((trim(PRODUCT(I,1)).EQ."S+" .AND. trim(PRODUCT(I,2)).EQ."e-") .OR.&
     &     (trim(PRODUCT(I,1)).EQ."e-" .AND. trim(PRODUCT(I,2)).EQ."S+"))) THEN
!C           Loop over all rays
!#ifdef MOCASSIN
!            RADSURFTOT = 0.0D0
!            DO K=0,NRAYS-1
!               RATE(I)=RATE(I) + SIPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),GAMMA(I),COLUMN(K,NS))
!               RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
!            ENDDO
!            RATE(I) = RATE(I)/RADSURFTOT
!#else
!            print *, "PART 7"
            DO K=1,NRAYS
!            DO K=0,NRAYS-1
               RATE(I)=RATE(I) + SIPDRATE(ALPHA(I),RAD_SURFACE(K),AV(K),GAMMA(I),COLUMN(K,NS))
            ENDDO
!#endif
            NRSI=I
            GOTO 10
         ENDIF

         IF(DUPLICATE(I).EQ.0) THEN
!C           Loop over all rays
!#ifdef MOCASSIN
!            RADSURFTOT = 0.0D0
!            DO K=0,NRAYS-1
!                RATE(I)=RATE(I) + ALPHA(I)*RAD_SURFACE(K)*EXP(-(GAMMA(I)*AV(K)))/2.0
!                RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
!!            ENDDO
 !           RATE(I) = RATE(I)/RADSURFTOT
!#else
!            print *, "part 8"
!            DO K=0,NRAYS-1
            DO K=1,NRAYS
 !              print *, ""
!               print *, "PART 8", RATE(I), ALPHA(I), RAD_SURFACE(K)
!               print *, GAMMA(I), AV(K)
!               print *, ""
                RATE(I)=RATE(I) + ALPHA(I)*RAD_SURFACE(K)*EXP(-(GAMMA(I)*AV(K)))/2.0
            ENDDO
!#endif
         ELSE IF(DUPLICATE(I).EQ.1) THEN
            J=I
            DO
               IF(TEMPERATURE.LE.RTMAX(J)) THEN
!C                 Loop over all rays
!#ifdef MOCASSIN
!                  RADSURFTOT = 0.0D0
!                  DO K=0,NRAYS-1
!                     RATE(J)=RATE(J) + ALPHA(J)*RAD_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))/2.0
!                     RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
!                  ENDDO
!                  RATE(J) = RATE(J)/RADSURFTOT
!#else
!                  print *, "PART 9"
!                  DO K=0,NRAYS-1
                  DO K=1,NRAYS
                     RATE(J)=RATE(J) + ALPHA(J)*RAD_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))/2.0
                  ENDDO
!#endif
                  EXIT
               ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
!C                 Loop over all rays
!#ifdef MOCASSIN
!                  RADSURFTOT = 0.0D0
!                  DO K=0,NRAYS-1
!                     RATE(J)=RATE(J) + ALPHA(J)*RAD_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))/2.0
!                     RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
!                  ENDDO
!                  RATE(J) = RATE(J)/RADSURFTOT
!#else                   
!                  print *, "PART 10"
!                  DO K=0,NRAYS-1
                  DO K=1,NRAYS
                     RATE(J)=RATE(J) + ALPHA(J)*RAD_SURFACE(K)*EXP(-(GAMMA(J)*AV(K)))/2.0
                  ENDDO
!#endif
                  EXIT
               ELSE
                  RATE(J)=0.0D0
                  J=J+1
               ENDIF
            ENDDO
         ENDIF
         GOTO 10

!C-----------------------------------------------------------------------

!C     Cosmic ray-induced ionization:
!         print *, "PART 11"
 2       IF(DUPLICATE(I).EQ.0) THEN
            RATE(I)=ALPHA(I)*ZETA
         ELSE IF(DUPLICATE(I).EQ.1) THEN
            J=I
            DO
               IF(TEMPERATURE.LE.RTMAX(J)) THEN
                  RATE(J)=ALPHA(J)*ZETA
                  EXIT
               ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
                  RATE(J)=ALPHA(J)*ZETA
                  EXIT
               ELSE
                  RATE(J)=0.0D0
                  J=J+1
               ENDIF
            ENDDO
         ENDIF
         GOTO 10

!#ifdef XRAYS
!
!11      IF(DUPLICATE(I).EQ.0) THEN
!            DO K=0,NRAYS-1
!!               RATE(I)=RATE(I)+ALPHA(I)*XRAD_SURFACE(K)*EXP(-sigma*AV(K)/AV_fac)
 !           ENDDO
 !        ELSE IF(DUPLICATE(I).EQ.1) THEN
 !           J=I
 !           DO
 !              IF(TEMPERATURE.LE.RTMAX(J)) THEN
 !                 DO K=0,NRAYS-1
 !                    RATE(I)=RATE(I)+ALPHA(I)*XRAD_SURFACE(K)*EXP(-sigma*AV(K)/AV_fac)
 !                 ENDDO
 !                 EXIT
 !              ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
 !                 DO K=0,NRAYS-1
 !                    RATE(I)=RATE(I)+ALPHA(I)*XRAD_SURFACE(K)*EXP(-sigma*AV(K)/AV_fac)
 !                 ENDDO
 !                 EXIT
 !              ELSE
 !                 RATE(J)=0.0D0
 !                 J=J+1
!               ENDIF
!            ENDDO
!         ENDIF
!         GOTO 10!

!#else
110     RATE(I)=0.0D0; GOTO 10
120     RATE(I)=0.0D0; GOTO 10
130     RATE(I)=0.0D0; GOTO 10
140     RATE(I)=0.0D0; GOTO 10
!#endif

!C-----------------------------------------------------------------------

!C     Photoreactions due to cosmic ray-induced secondary photons:
!         print *, "PART 12"
 3       IF(DUPLICATE(I).EQ.0) THEN
            RATE(I)=ALPHA(I)*ZETA*(TEMPERATURE/300.0D0)**BETA(I)&
     &           *GAMMA(I)/(1.0D0-OMEGA)
         ELSE IF(DUPLICATE(I).EQ.1) THEN
            J=I
            DO
               IF(TEMPERATURE.LE.RTMAX(J)) THEN
                  RATE(J)=ALPHA(J)*ZETA*(TEMPERATURE/300.0D0)**BETA(J)&
     &                 *GAMMA(J)/(1.0D0-OMEGA)
                  EXIT
               ELSE IF(DUPLICATE(J+1).LT.DUPLICATE(J)) THEN
                  RATE(J)=ALPHA(J)*ZETA*(TEMPERATURE/300.0D0)**BETA(J)&
     &                 *GAMMA(J)/(1.0D0-OMEGA)
                  EXIT
               ELSE
                  RATE(J)=0.0D0
                  J=J+1
               ENDIF
            ENDDO
         ENDIF
         GOTO 10

!C-----------------------------------------------------------------------

!C     Freeze-out of neutral species:

 4       IF(BETA(I).EQ.0.0D0) THEN
            CION=1.0D0
         ELSE IF(BETA(I).EQ.1.0D0) THEN
            CION=1.0D0+16.71D-4/(GRAIN_RADIUS*TEMPERATURE)
         ELSE
            CION=0.0D0
         ENDIF
         STICKING=0.3D0
!         print *, "PART 13"
         RATE(I)=ALPHA(I)*4.57D4*2.4D-22*SQRT(TEMPERATURE/GAMMA(I))*CION*STICKING
         GOTO 10

!C-----------------------------------------------------------------------

!C     Freeze-out of singly charged positive ions:

 5       IF(BETA(I).EQ.0.0D0) THEN
            CION=1.0D0
         ELSE IF(BETA(I).EQ.1.0D0) THEN
            CION=1.0D0+16.71D-4/(GRAIN_RADIUS*TEMPERATURE)
         ELSE
            CION=0.0D0
         ENDIF
         STICKING=0.3D0
!         print *, "PART 14"
         RATE(I)=ALPHA(I)*4.57D4*2.4D-22*SQRT(TEMPERATURE/GAMMA(I))*CION*STICKING
         GOTO 10

!C-----------------------------------------------------------------------

!C     Desorption due to cosmic ray heating:

!CC     Treatment of Hasegawa & Herbst (1993, MNRAS, 261, 83, Equation 15)
!C 6       RATE(I)=ALPHA(I)*ZETA

!C     Treatment of Roberts et al. (2007, MNRAS, 382, 773, Equation 3)
 6       IF(GAMMA(I).LE.1210.0D0) THEN
            YIELD=1.0D5 ! Number of adsorbed molecules released per cosmic ray impact
         ELSE
            YIELD=0.0D0
         ENDIF
         FLUX=2.06D-3 ! Flux of iron nuclei cosmic rays (in cm^-2 s^-1)
!         print *, "PART 15"
         RATE(I)=FLUX*ZETA*2.4D-22*YIELD
         GOTO 10

!C-----------------------------------------------------------------------

!C     Photodesorption:

 7       IF(TEMPERATURE.LT.50.0D0) THEN
            YIELD=3.5D-3
         ELSE IF(TEMPERATURE.LT.85.0D0) THEN
            YIELD=4.0D-3
         ELSE IF(TEMPERATURE.LT.100.0D0) THEN
            YIELD=5.5D-3
         ELSE
            YIELD=7.5D-3
         ENDIF
!C         FLUX=1.0D8 ! Flux of FUV photons in the unattenuated Habing field (in photons cm^-2 s^-1)
         FLUX=1.7D8 ! Flux of FUV photons in the unattenuated Draine field (in photons cm^-2 s^-1)
!C        Loop over all rays
!#ifdef MOCASSIN
!         RADSURFTOT = 0.0D0
!         DO K=0,NRAYS-1
!            RATE(I)=RATE(I) + FLUX*RAD_SURFACE(K)*EXP(-(1.8D0*AV(K)))*2.4D-22*YIELD
!            RADSURFTOT = RADSURFTOT + RAD_SURFACE(K)
!         ENDDO
!         RATE(I) = RATE(I)/RADSURFTOT
!#else       
!         print *, "PART 17"
         DO K=1,NRAYS
!         DO K=0,NRAYS-1
            RATE(I)=RATE(I) + FLUX*RAD_SURFACE(K)*EXP(-(1.8D0*AV(K)))*2.4D-22*YIELD
         ENDDO
!#endif
         GOTO 10

!C-----------------------------------------------------------------------

!C     Thermal desorption:

!C     Treatment of Hasegawa, Herbst & Leung (1992, ApJS, 82, 167, Equations 2 & 3)
!         print *, "PART 18"
 8       RATE(I)=SQRT(2.0D0*1.5D15*KB/(PI**2*AU)*ALPHA(I)/GAMMA(I))&
     &        *EXP(-(ALPHA(I)/DUST_TEMPERATURE))
         GOTO 10

!C-----------------------------------------------------------------------

!C     Grain mantle reactions:
!         print *, "PART 19"
 9       RATE(I)=ALPHA(I)
         GOTO 10

!C-----------------------------------------------------------------------

!C     Check that the rate is physical (0<RATE(I)<1) and produce an error
!C     message if not. Impose a lower cut-off on all rate coefficients to
!C     prevent the problem becoming too stiff. Rates less than 1E-99 are
!C     set to zero. Grain-surface reactions and desorption mechanisms are
!C     allowed rates greater than 1.
 10      IF(RATE(I).LT.0.0D0) THEN
           PRINT *,'ERROR! Negative rate for reaction',I
           STOP
         ENDIF
         IF(RATE(I).GT.1.0D0 .AND. REACTANT(I,1)(1:1).NE."G") THEN
           WRITE(10,*)'WARNING! Rate is too large for reaction',I
           WRITE(10,*)'RATE =',RATE(I)
           RATE(I)=1.0D0
         ENDIF
         IF(RATE(I).LT.1.0D-99) RATE(I)=0.0D0
!C     End of loop over rates
      ENDDO

      RETURN
      END SUBROUTINE CALCULATE_REACTION_RATES


!C=======================================================================
!C
!C     H2 photodissociation rate taking into account
!C     self-shielding and grain extinction
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0  = Unattenuated photodissociation rate (in cm^3/s)
!C     G0  = Incident FUV field (in Draine units)
!C     AV  = visual extinction (in magnitudes)
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     H2PDRATE = H2 photodissociation rate taking into
!C                account self-shielding and grain extinction
!C     DOPW     = Doppler linewidth (in Hz) of a typical transition
!C                (assuming turbulent broadening with b=3 km/s)
!C     RADW     = radiative linewidth (in Hz) of a typical transition
!C     LAMBDA   = wavelength (in Å) of a typical transition
!C
!C     Functions called:
!C     H2SHIELD = H2 self-shielding function
!C     S!CATTER  = attenuation due to scattering by dust
!C
!C-----------------------------------------------------------------------
      FUNCTION H2PDRATE(K0,G0,AV,NH2)
!      IMPLICIT NONE
!      real(double) :: K0,G0,AV,NH2
!      real(double) :: LAMBDA,SCATTER
!!C      DOUBLE PRECISION DOPW,RADW
!!C      DOUBLE PRECISION H2SHIELD1
!      real(double) :: H2SHIELD2
!     use definitions
!     use healpix_types
!     use maincode_module, only : v_turb
!     use global_module, only : nh2
    implicit none
     real(double) :: H2PDRATE
     real(double), intent(in) :: k0, g0, av
!     integer, intent(in) :: nh2
     real(double), intent(in) :: nh2
     real(double) :: lambda!, scatter!, h2shield2
     real(double) :: dopw, radw!,! !h2shield1!, v_turb
     real(double) :: h2shieldVAL, scatterVAL

     LAMBDA=1000.0D0
      DOPW=V_TURBPDR/(LAMBDA*1.0D-8)
      RADW=8.0D7

      h2shieldVAL = H2SHIELD1(NH2,DOPW,RADW)
      scatterVAL = SCATTER(AV,LAMBDA, START=.true.)
!     Calculate the H2 photodissociation rate (H2PDRATE)
      H2PDRATE=K0*G0*h2shieldVAL*scatterVAL/2.0
!      H2PDRATE=K0*G0*H2SHIELD2(NH2)*SCATTER(AV,LAMBDA)/2.0

      RETURN
    END FUNCTION H2PDRATE
!C=======================================================================

!C=======================================================================
!C
!C     !CO photodissociation rate taking into
!C     account shielding and grain extinction
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0  = Unattenuated photodissociation rate (in cm^3/s)
!C     G0  = Incident FUV field (in Draine units)
!C     AV  = visual extinction (in magnitudes)
!C     N!CO = !CO column density (in cm^-2)
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     !COPDRATE = !CO photodissociation rate taking into
!C                account self-shielding and grain extinction
!C     LAMBDA   = wavelength (in Å) of a typical transition
!C
!C     Functions called:
!C     LBAR     = function to determine the wavelength
!C     !COSHIELD = !CO shielding function
!C     S!CATTER  = attenuation due to scattering by dust
!C
!C-----------------------------------------------------------------------
      FUNCTION COPDRATE(K0,G0,AV,NCO,NH2)

!      IMPLICIT NONE
!      real(double) :: K0,G0,AV,NCO,NH2
!      real(double) :: LAMBDA,LBAR
!      real(double) :: COSHIELD,SCATTER

 !    use definitions
!     use healpix_types
!     use global_module, only : nh2
     implicit none
     real(double) :: copdrate
     real(double), intent(in) :: k0, g0, av, nco
!     integer, intent(in) :: nh2
     real(double), intent(in) :: nh2
     real(double) :: lambda, coshieldVAL, scatterVAL


      LAMBDA=LBAR(NCO,NH2)
      coshieldVAL = COSHIELD(NCO,NH2)
!      print *, "CALLING SCATTER WITH ", LAMBDA
      scatterVAL = SCATTER(AV,LAMBDA, start=.true.)
!C     Calculate the CO photodissociation rate (COPDRATE)
      COPDRATE=K0*G0*COSHIELDVAL*SCATTERVAL/2.0

      RETURN
    END FUNCTION COPDRATE
!C=======================================================================

!C=======================================================================
!C
!C     !CI photoionization rate taking into account grain extinction
!C     and shielding by !CI and H2 lines, adopting the treatment of
!C     Kamp & Bertoldi (2000, A&A, 353, 276, Equation 8)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0   = Unattenuated photoionization rate (in cm^3/s)
!C     G0   = Incident FUV field (in Draine units)
!C     AV   = visual extinction (in magnitudes)
!C     KAV  = tau(λ)/tau(V) correction factor
!C     N!CI  = !CI column density (in cm^-2)
!C     NH2  = H2 column density (in cm^-2)
!C     TGAS = gas temperature (in K)
!C
!C     Program variables:
!C     !CIPDRATE = !CI photoionization rate taking into
!C                account shielding and grain extinction
!C     TAU!C     = optical depth in the !CI absorption band
!C
!C-----------------------------------------------------------------------
       FUNCTION CIPDRATE(K0,G0,AV,KAV,NCI,NH2,TGAS)

!      IMPLICIT NONE
!      real(double) :: K0,G0,AV,KAV,NCI,NH2,TGAS
!      real(double) :: TAUC

!     use definitions
!     use healpix_types
!     use global_module, only : nh2
     implicit none
     real(double) :: cipdrate
     real(double), intent(in) :: K0,G0,AV,KAV,NCI,TGAS
!     integer, intent(in) :: nh2
     real(double), intent(in) :: nh2
     real(double) :: tauc


!C     !Calculate the optical depth in the !CI absorption band, accounting
!C     for grain extinction and shielding by !CI and overlapping H2 lines
      TAUC=KAV*AV+1.1D-17*NCI+(0.9D0*TGAS**0.27D0*(NH2/1.59D21)**0.45D0)

!C     Calculate the CI photoionization rate (CIPDRATE)
      CIPDRATE=K0*G0*EXP(-TAUC)/2.0

      RETURN
    END FUNCTION CIPDRATE
!C=======================================================================

!C=======================================================================
!C
!C     SI photoionization rate -- needs to be implemented!
!C     For now, use the standard expression for photorates
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     K0   = Unattenuated photoionization rate (in cm^3/s)
!C     G0   = Incident FUV field (in Draine units)
!C     AV   = visual extinction (in magnitudes)
!C     KAV  = tau(λ)/tau(V) correction factor
!C     NSI  = SI column density (in cm^-2)
!C
!C     Program variables:
!C     SIPDRATE = SI photoionization rate taking into
!C                account shielding and grain extinction
!C     TAUS     = optical depth in the SI absorption band
!C
!C-----------------------------------------------------------------------
      FUNCTION SIPDRATE(K0,G0,AV,KAV,NSI)

!      IMPLICIT NONE
!      real(double) K0,G0,AV,KAV,NSI
!      real(double) TAUS

 !    use definitions
!     use healpix_types
     implicit none
     real(double) :: sipdrate
     real(double), intent(in) :: G0,AV,KAV!,NSI
     real(double) :: taus, K0, NSI
!C     Calculate the optical depth in the SI absorption band, accounting
!C     for grain extinction and shielding by ???
     K0 = K0 !otherwise unused variable - thaw
     NSI = NSI
      TAUS=KAV*AV

!C     Calculate the SI photoionization rate (SIPDRATE)
      SIPDRATE=K0*G0*EXP(-TAUS)

      RETURN
    END FUNCTION SIPDRATE
!C=======================================================================

!C=======================================================================
!C
!C     H2 line self-shielding, adopting the treatment of
!C     Federman, Glassgold & Kwan (1979, ApJ, 227, 466)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     NH2  = H2 column density (in cm^-2)
!C     DOPW = Doppler linewidth (in Hz)
!C     RADW = radiative linewidth (in Hz)
!C
!C     Program variables:
!C     H2SHIELD1 = total self-shielding function containing
!C                 both Doppler and radiative contributions
!C     FPARA     = fraction of H2 in para state: 1/(1+o/p ratio)
!C     FOS!C      = oscillator strength of a typical transition
!C     TAUD      = parameter tauD (eq. A7) in Federman's paper
!C                 (optical depth at line centre)
!C     R         = parameter r  (eq. A2) in Federman's paper
!C     T         = parameter t1 (eq. A6) in Federman's paper
!C     U         = parameter u1 (eq. A6) in Federman's paper
!C     JD        = parameter JD (eq. A8) in Federman's paper
!C                 (Doppler contribution to self-shielding)
!C     JR        = parameter JR (eq. A9) in Federman's paper
!C                (radiative contribution to self-shielding)
!C
!C-----------------------------------------------------------------------
      FUNCTION H2SHIELD1(NH2,DOPW,RADW)

!      IMPLICIT NONE
!      real(double) ::  NH2,DOPW,RADW
!      real(double) ::  FPARA,FOSC,TAUD
!      real(double) ::  R,T,U,JD,JR

!     use definitions
!     use healpix_types
!     use global_module, only : nh2
     implicit none
     real(double) :: h2shield1
!     integer, intent(in) :: nh2
     real(double), intent(in) :: nh2
     real(double), intent(in) ::DOPW,RADW
     real(double) :: FPARA, FOSC, TAUD, R, T, U, JD, JR


!C     Calculate the optical depth at line centre = N(H2)*f_para*(πe^2/mc)*f/(√πß) ≈ N(H2)*f_para*(1.5E-2)*f/ß
      FPARA=0.5D0 ! (assume o/p ratio=1)
      FOSC=1.0D-2
      TAUD=NH2*FPARA*(1.497358985D-2)*FOSC/DOPW

!C     Calculate the Doppler core contribution to the self-shielding (JD)
      IF(TAUD.EQ.0.0D0) THEN
         JD=1.0D0
      ELSE IF(TAUD.LT.2.0D0) THEN
         JD=EXP(-(0.666666667D0*TAUD))
      ELSE IF(TAUD.LT.10.0D0) THEN
         JD=0.638D0*TAUD**(-1.25D0)
      ELSE IF(TAUD.LT.100.0D0) THEN
         JD=0.505D0*TAUD**(-1.15D0)
      ELSE
         JD=0.344D0*TAUD**(-1.0667D0)
      ENDIF

!C     Calculate the radiative wing contribution to self-shielding (JR)
      IF(RADW.EQ.0.0D0) THEN
         JR=0.0D0
      ELSE
         R=RADW/(1.772453851D0*DOPW)
         T=3.02D0*((R*1.0D3)**(-0.064D0))
         U=SQRT(TAUD*R)/T
         JR=R/(T*SQRT(0.785398163D0+U**2))
      ENDIF

!C     Calculate the total self-shielding function (H2SHIELD1)
      H2SHIELD1=JD+JR

      RETURN
    END FUNCTION H2SHIELD1
!C=======================================================================

!C=======================================================================
!C
!C     H2 line shielding, using the computed values listed in
!C     Lee et al. (1996, A&A, 311, 690, Table 10)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     H2SHIELD2 = total H2 shielding factor containing
!C                 contributions from both H2 and H lines
!C                 from spline interpolation over the grid
!C     SH2_GRID  = H2 shielding factors from Lee et al. (1996)
!C                 as a function of H2 column density
!C     SH2_DERIV = 2nd derivative of SH2_GRID values from SPLINE
!C     !COL_GRID  = H2 column densities (in cm^-2)
!C     NUMH2     = number of entries in the table
!C     START     = .TRUE. when H2SHIELD2 is first called
!C
!C     Functions called:
!C     SPLINE =
!C     SPLINT =
!C
!C-----------------------------------------------------------------------
      FUNCTION H2SHIELD2(NH2)

      IMPLICIT NONE
      LOGICAL :: START
      integer :: NUMH2
!      real(double) ::  NH2
        real(double) ::  SH2_DERIV(105)
  real(double), dimension(105) :: COL_GRID=(/&
     &              0.000D+00,3.690D+11,3.715D+12,3.948D+13,1.233D+14, &
     &              2.536D+14,4.342D+14,6.653D+14,6.689D+14,9.075D+14, &
     &              1.234D+15,1.631D+15,2.105D+15,2.363D+15,2.899D+15, &
     &              3.207D+15,3.848D+15,4.636D+15,5.547D+15,6.604D+15, &
     &              7.855D+15,9.368D+15,1.122D+16,1.352D+16,1.643D+16, &
     &              2.017D+16,2.515D+16,3.190D+16,4.128D+16,5.439D+16, &
     &              7.315D+16,1.009D+17,1.432D+17,2.092D+17,3.123D+17, &
     &              4.738D+17,5.388D+17,8.935D+17,1.381D+18,2.164D+18, &
     &              3.330D+18,5.024D+18,7.404D+18,9.029D+18,1.316D+19, &
     &              1.813D+19,2.453D+19,3.248D+19,4.216D+19,5.370D+19, &
     &              6.722D+19,8.277D+19,9.894D+19,1.186D+20,1.404D+20, &
     &              1.644D+20,1.908D+20,2.197D+20,2.510D+20,2.849D+20, &
     &              3.214D+20,3.604D+20,4.019D+20,4.456D+20,4.915D+20, &
     &              5.393D+20,5.886D+20,6.392D+20,6.909D+20,7.433D+20, &
     &              7.965D+20,8.505D+20,9.056D+20,9.627D+20,1.011D+21, &
     &              1.068D+21,1.125D+21,1.185D+21,1.250D+21,1.327D+21, &
     &              1.428D+21,1.578D+21,1.851D+21,2.128D+21,2.298D+21, &
     &              2.389D+21,2.459D+21,2.519D+21,2.571D+21,2.618D+21, &
     &              2.707D+21,2.790D+21,2.887D+21,3.001D+21,3.139D+21, &
     &              3.303D+21,3.497D+21,3.722D+21,3.983D+21,4.283D+21, &
     &              4.644D+21,5.127D+21,5.945D+21,8.205D+21,1.015D+22/)

  real(double), dimension(105) :: SH2_GRID=(/&
     &              1.000D+00,9.983D-01,9.853D-01,8.761D-01,7.199D-01, &
     &              5.728D-01,4.455D-01,3.431D-01,3.418D-01,2.732D-01, &
     &              2.110D-01,1.619D-01,1.236D-01,1.084D-01,8.447D-02, &
     &              7.410D-02,5.774D-02,4.416D-02,3.390D-02,2.625D-02, &
     &              2.048D-02,1.606D-02,1.264D-02,9.987D-03,7.937D-03, &
     &              6.343D-03,5.088D-03,4.089D-03,3.283D-03,2.640D-03, &
     &              2.130D-03,1.725D-03,1.397D-03,1.129D-03,9.097D-04, &
     &              7.340D-04,6.883D-04,5.377D-04,4.352D-04,3.475D-04, &
     &              2.771D-04,2.205D-04,1.753D-04,1.549D-04,1.210D-04, &
     &              9.666D-05,7.705D-05,6.148D-05,4.904D-05,3.909D-05, &
     &              3.112D-05,2.473D-05,1.997D-05,1.578D-05,1.244D-05, &
     &              9.769D-06,7.634D-06,5.932D-06,4.581D-06,3.515D-06, &
     &              2.679D-06,2.029D-06,1.527D-06,1.144D-06,8.523D-07, &
     &              6.332D-07,4.693D-07,3.475D-07,2.574D-07,1.907D-07, &
     &              1.413D-07,1.047D-07,7.739D-08,5.677D-08,4.386D-08, &
     &              3.227D-08,2.385D-08,1.750D-08,1.248D-08,8.389D-09, &
     &              5.026D-09,2.382D-09,6.259D-10,1.653D-10,7.399D-11, &
     &              4.824D-11,3.474D-11,2.633D-11,2.069D-11,1.663D-11, &
     &              1.099D-11,7.506D-12,4.825D-12,2.864D-12,1.534D-12, &
     &              7.324D-13,3.087D-13,1.135D-13,3.591D-14,9.689D-15, &
     &              2.045D-15,2.618D-16,8.918D-18,3.041D-21,1.739D-23/)


!      COMMON /STATUS/START
!      COMMON /H2GRID/SH2_GRID,SH2_DERIV,COL_GRID,NUMH2
!     use definitions
!     use healpix_types
!     use uclpdr_module, only : start, numh2, COL_GRID, SH2_GRID, SH2_DERIV
!     use global_module, only : nh2
!     implicit none
     real(double) :: h2shield2
!     integer, intent(in) :: nh2
     real(double), intent(inout) :: nh2
     START = .true.
     
      NUMH2=105
      IF(START) CALL SPLINE_PDR(COL_GRID,SH2_GRID,NUMH2, &
     &                      1.0D30,1.0D30,SH2_DERIV)
      IF(NH2.LT.COL_GRID(1))     NH2=COL_GRID(1)
      IF(NH2.GT.COL_GRID(NUMH2)) NH2=COL_GRID(NUMH2)
      CALL SPLINT_PDR(COL_GRID,SH2_GRID,SH2_DERIV,NUMH2,NH2,H2SHIELD2)
      IF(H2SHIELD2.LT.0.0D0) H2SHIELD2=0.0D0

      RETURN
    END FUNCTION H2SHIELD2
!C=======================================================================

!C=======================================================================
!C
!C     12!CO line shielding, using the computed values listed in
!C     van Dishoeck & Black (1988, ApJ, 334, 771, Table 5)
!C
!C     Appropriate shielding factors are determined by performing a
!C     2-dimensional spline interpolation over the values listed in
!C     Table 5 of van Dishoeck & Black, which include contributions
!C     from self-shielding and H2 screening
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     N!CO = !CO column density (in cm^-2)
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     !COSHIELD  = total 12!CO shielding factor containing
!C                 contributions from both H2 and !CO lines
!C                 from 2D spline interpolation over the grid
!C     S!CO_GRID  = log10 values of the 12!CO shielding factors
!C                 from van Dishoeck & Black (1988) as a function
!C                 of !CO column density (1st index) and H2 column
!C                 density (2nd index)
!C     S!CO_DERIV = 2nd derivative of S!CO_GRID values from SPLIE2
!C     N!CO_GRID  = log10 values of !CO column densities (in cm^-2)
!C     NH2_GRID  = log10 values of H2 column densities (in cm^-2)
!C     DIM!CO     = number of !CO column densities
!C     DIMH2     = number of H2 column densities
!C     START     = .TRUE. when !COSHIELD is first called
!C
!C     Functions called:
!C     SPLIE2 =
!C     SPLIN2 =
!C
!C-----------------------------------------------------------------------
    FUNCTION COSHIELD(NCO,NH2)
      
      IMPLICIT NONE
      LOGICAL, save :: START=.true.
      integer :: DIMCO,DIMH2
      real(double) ::  NCO,NH2
      !      real(double) ::  LOGNCO,LOGNH2
!      real(double) ::  NCO_GRID(8),NH2_GRID(6)
!      real(double) ::  SCO_GRID(8,6),SCO_DERIV(8,6)

      real(double), dimension(8) :: NCO_GRID=(/&
           &12.0D0,13.0D0,14.0D0,15.0D0, 16.0D0,17.0D0,18.0D0,19.0D0/)
      real(double), dimension(6) :: NH2_GRID=(/&
           &18.0D0,19.0D0,20.0D0,21.0D0,22.0D0,23.0D0/)
      real(double) :: SCO_GRID(1:8,1:6)
      real(double) :: SCO_DERIV(1:8,1:6)

      !      COMMON /STATUS/START
      !      COMMON /COGRID/SCO_GRID,SCO_DERIV,NCO_GRID,NH2_GRID,DIMCO,DIMH2
      
      !     use definitions
      !     use healpix_types
      !     use global_module, only : NCO, NH2
      !    use uclpdr_module, only : start, NCO_GRID, NH2_GRID, SCO_GRID, &
      !           & DIMCO, SCO_DERIV, DIMH2!, SCO_DERIV
      !     implicit none
     real(double) :: COSHIELD
     real(double) :: LOGNCO, LOGNH2
!     integer, intent(in) :: NCO, NH2
!     real(double), intent(in) :: nh2, nco
     DIMH2 = 6
     DIMCO=8


      IF(START) THEN
         CALL SPLIE2(NCO_GRID,NH2_GRID,SCO_GRID,DIMCO,DIMH2,SCO_DERIV)
         START=.FALSE.
      ENDIF

      LOGNCO=DLOG10(NCO+1.0D0)
      LOGNH2=DLOG10(NH2+1.0D0)

      IF(LOGNCO.LT.NCO_GRID(1)) LOGNCO=NCO_GRID(1)
      IF(LOGNH2.LT.NH2_GRID(1)) LOGNH2=NH2_GRID(1)
      IF(LOGNCO.GT.NCO_GRID(DIMCO)) LOGNCO=NCO_GRID(DIMCO)
      IF(LOGNH2.GT.NH2_GRID(DIMH2)) LOGNH2=NH2_GRID(DIMH2)

      CALL SPLIN2(NCO_GRID,NH2_GRID,SCO_GRID,SCO_DERIV,&
     &            DIMCO,DIMH2,LOGNCO,LOGNH2,COSHIELD)
      COSHIELD=10.0D0**COSHIELD

      RETURN
    END FUNCTION COSHIELD
!C=======================================================================

!C=======================================================================
!C
!C     Scattering by dust grains, adopting the treatment of
!C     Wagenblast & Hartquist (1989, MNRAS, 237, 1019) and
!C     Flannery, Roberge & Rybicki (1980, ApJ, 236, 598)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     AV     = visual extinction (in magnitudes)
!C     LAMBDA = wavelength (in Å) of incident radiation
!C
!C     Program variables:
!C     S!CATTER = attenuation factor describing the influence of
!C               grain scattering on the FUV flux, dependening
!C               on the total column density and wavelength of
!C               light (assuming albedo=0.3 gscat=0.8)
!C     TAUV    = optical depth at visual wavelength (λ=5500Å)
!C     TAUL    = optical depth at wavelength LAMBDA
!C     A(0)    = a(0)*exp(-k(0)*tau)
!C             = relative intensity decrease for 0 < tau < 1
!C     A(I)    = ∑ a(i)*exp(-k(i)*tau) for i=1,5
!C               relative intensity decrease for tau ≥ 1
!C     K(0)    = see A0
!C     K(I)    = see A(I)
!C
!C     Functions called:
!C     XLAMBDA = function to determine tau(λ)/tau(V)
!C
!C-----------------------------------------------------------------------
      FUNCTION SCATTER(AV,LAMBDA,START)

      IMPLICIT NONE
!      integer :: I
!      real(double) ::  AV,LAMBDA
!      real(double) ::  TAUV,TAUL
!      real(double) ::  A(0:5),K(0:5),EXPONENT
!      real(double) ::  XLAMBDA
!      DATA A/1.000D0,2.006D0,-1.438D0,0.7364D0,-0.5076D0,-0.0592D0/
!      DATA K/0.7514D0,0.8490D0,1.013D0,1.282D0,2.005D0,5.832D0/

!     use definitions
!     use healpix_types
!     implicit none
     real(double) :: scatter
     real(double), intent(in) :: AV
     real(double), dimension(0:5), save :: A = (/&
           &1.000D0,2.006D0,-1.438D0,0.7364D0,-0.5076D0,-0.0592D0/)
     real(double), dimension(0:5), save :: K = (/&
           &0.7514D0,0.8490D0,1.013D0,1.282D0,2.005D0,5.832D0/)
     real(double) :: EXPONENT, XLAMBDAVAL, LAMBDA
     integer :: i
     logical :: start
     real(double) :: TAUL, TAUV

!C     Calculate the optical depth at visual wavelength
      TAUV=AV/1.086D0

!C     Convert the optical depth to that at the desired wavelength
!      print *, "CALLING XLAMBDAVAL WITH ", LAMBDA
      XLAMBDAVAL = XLAMBDA(LAMBDA, START)
      TAUL=TAUV*XLAMBDAVAL

!C     Calculate the attenuation due to scattering by dust (SCATTER)
      SCATTER=0.0D0
      IF(TAUL.LT.1.0D0) THEN
         EXPONENT=K(0)*TAUL
         IF(EXPONENT.LT.100.0D0) THEN
            SCATTER=A(0)*EXP(-EXPONENT)
         ENDIF
      ELSE
         DO I=1,5
            EXPONENT=K(I)*TAUL
            IF(EXPONENT.LT.100.0D0) THEN
               SCATTER=SCATTER+A(I)*EXP(-EXPONENT)
            ENDIF
         ENDDO
      ENDIF

      RETURN
    END FUNCTION SCATTER
!C=======================================================================

!C=======================================================================
!C
!C     Determine the ratio of the optical depth at a given wavelength to
!C     that at visual wavelength (λ=5500Å) using the extinction curve of
!C     Savage & Mathis (1979, ARA&A, 17, 73, Table 2)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     LAMBDA  = wavelength (in Å)
!C
!C     Program variables:
!C     XLAMBDA = value of tau(λ)/tau(V) at the desired wavelength
!C               (by spline interpolation over a table of values)
!C     L_GRID  = wavelengths listed in Table 2 of Savage & Mathis
!C     X_GRID  = tau(λ)/tau(V) values, determined by dividing the
!C               Aλ/E(B-V) values in Table 2 by R=AV/E(B-V)=3.1
!C     X_DERIV = 2nd derivative of X_GRID values from SPLINE
!C     N_GRID  = number of wavelengths
!C
!C     Functions called:
!C     SPLIE =
!C     SPLIN =
!C
!C-----------------------------------------------------------------------
      FUNCTION XLAMBDA(LAMBDA, START)

      IMPLICIT NONE
      LOGICAL :: START
      integer ::  N_GRID
!      real(double) ::  LAMBDA
!      real(double) ::  X_GRID(30),X_DERIV(30)
  real(double), dimension(30) :: X_GRID=(/&
     &            5.76D0,5.18D0,4.65D0,4.16D0,3.73D0, &
     &            3.40D0,3.11D0,2.74D0,2.63D0,2.62D0, &
     &            2.54D0,2.50D0,2.58D0,2.78D0,3.01D0, &
     &            3.12D0,2.86D0,2.58D0,2.35D0,2.00D0, &
     &            1.58D0,1.42D0,1.32D0,1.00D0,0.75D0, &
     &            0.48D0,0.28D0,0.12D0,0.05D0,0.00D0/)

  real(double), dimension(30) :: L_GRID=(/&
     &             910.0D0, 950.0D0,1000.0D0,1050.0D0,1110.0D0, &
     &            1180.0D0,1250.0D0,1390.0D0,1490.0D0,1600.0D0, &
     &            1700.0D0,1800.0D0,1900.0D0,2000.0D0,2100.0D0, &
     &            2190.0D0,2300.0D0,2400.0D0,2500.0D0,2740.0D0, &
     &            3440.0D0,4000.0D0,4400.0D0,5500.0D0,7000.0D0, &
     &            9000.0D0,12500.0D0,22000.0D0,34000.0D0,1.0D9/)

!      COMMON /STATUS/START
!      COMMON /TAUGRID/L_GRID,X_GRID,X_DERIV,N_GRID

!     use definitions
!     use healpix_types
!     use uclpdr_module, only : start, N_GRID, L_GRID, X_GRID, X_DERIV
!     implicit none
     real(double) :: xlambda
     real(double), intent(inout) :: lambda

      N_GRID=30
!      print *, "L_GRID", L_GRID


!C     Find the appropriate value for XLAMBDA using spline interpolation
      IF(START) CALL SPLINE_PDR(L_GRID,X_GRID,N_GRID,1.0D30,1.0D30,X_DERIV)
!      print *, "XLAMBDAVAL HAS ", lambda, l_grid(n_grid)
      IF(LAMBDA.LT.L_GRID(1))      LAMBDA=L_GRID(1)
      IF(LAMBDA.GT.L_GRID(N_GRID)) LAMBDA=L_GRID(N_GRID)
!      print *, "XLAMBDAVAL HAS ", lambda
      CALL SPLINT_PDR(L_GRID,X_GRID,X_DERIV,N_GRID,LAMBDA,XLAMBDA)

      RETURN
    END FUNCTION XLAMBDA
!C=======================================================================

!C=======================================================================
!C
!C     !Calculate the mean wavelength (in Å) of the 33 dissociating bands,
!C     weighted by their fractional contribution to the total shielding
!C     van Dishoeck & Black (1988, ApJ, 334, 771, Equation 4)
!C
!C-----------------------------------------------------------------------
!C
!C     Input parameters:
!C     N!CO = !CO column density (in cm^-2)
!C     NH2 = H2 column density (in cm^-2)
!C
!C     Program variables:
!C     LBAR = mean wavelength (in Å)
!C     U    = log10(N!CO)
!C     W    = log10(NH2)
!C
!C-----------------------------------------------------------------------
      FUNCTION LBAR(NCO,NH2)

      IMPLICIT NONE
!      real(double) ::  NCO,NH2
!      real(double) ::  U,W


!     use definitions
!     use healpix_types
!     use global_module, only : NCO, NH2
!     implicit none
     real(double) :: lbar
     real(double) :: U,W
!
!     integer, intent(in) :: NCO, NH2
     real(double) :: nh2,nco

!     nh2d = dble(nh2)
!     ncod = dble(ncod)
!     print *, "nco pork",  nco
!     print *, "nh2 ham",  nh2

      U=DLOG10(NCO+1.0D0)
      W=DLOG10(NH2+1.0D0)

      LBAR=(5675.0D0 - 200.6D0*W) &
     &    - (571.6D0 - 24.09D0*W)*U &
     &   + (18.22D0 - 0.7664D0*W)*U**2

!C     LBAR cannot be larger than the wavelength of band 33 (1076.1Å)
!C     and cannot be smaller than the wavelength of band 1 (913.6Å)
      IF(LBAR.LT.913.6D0)  LBAR=913.6D0
      IF(LBAR.GT.1076.1D0) LBAR=1076.1D0
!      print *, "LBAR IS ", LBAR

      RETURN
    END FUNCTION LBAR
!C=======================================================================


!=======================================================================
!
!     Copied from Numerical Recipes
!
! Given a tabulated function YA (of size MxN) and tabulated independent
! variables X1A (M values) and X2A (N values), this routine constructs
! one-dimensional natural cubic splines of the rows of YA and returns
! the second derivatives in the array Y2A.
!
!-----------------------------------------------------------------------
      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)

!      USE DEFINITIONS
!      USE HEALPIX_TYPES

      IMPLICIT NONE

      integer, INTENT(IN) :: M,N
      REAL(double), INTENT(IN)     :: X2A(1:N),YA(1:M,1:N)
      REAL(double), INTENT(OUT)    :: Y2A(1:M,1:N)

      integer :: I,J
      REAL(double) :: YTMP(1:N),Y2TMP(1:N), X1A(1:M)

      X1A = X1A ! otherwise unused variable - thaw

      DO I=1,M
         DO J=1,N
            YTMP(J)=YA(I,J)
         ENDDO
!        Values of 1.0D30 indicate a natural spline
         CALL SPLINE_PDR(X2A,YTMP,N,1.0D30,1.0D30,Y2TMP)
         DO J=1,N
            Y2A(I,J)=Y2TMP(J)
         ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SPLIE2
!=======================================================================



!=======================================================================
!
!     Given X1A, X2A, YA, M, N (as described in SPLIE2) and Y2A (as
!     produced by that routine), and given a desired interpolating
!     point (X1,X2), this routine returns an interpolated function
!     value Y by performing a bicubic spline interpolation.
!
!-----------------------------------------------------------------------
      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)

!      USE DEFINITIONS
!      USE HEALPIX_TYPES

      IMPLICIT NONE

      integer, INTENT(IN) :: M,N
      REAL(double), INTENT(IN)     :: X1A(1:M),X2A(1:N),YA(1:M,1:N),Y2A(1:M,1:N)
      REAL(double), INTENT(IN)     :: X1,X2
      REAL(double), INTENT(OUT)    :: Y

      integer :: I,J
      REAL(double) :: YTMP(1:N),Y2TMP(1:N),YYTMP(1:M),YY2TMP(1:M)

!     Perform M evaluations of the row splines constructed by
!     SPLIE2 using the one-dimensional spline evaluator SPLINT
      DO I=1,M
         DO J=1,N
            YTMP(J)=YA(I,J)
            Y2TMP(J)=Y2A(I,J)
         ENDDO
         CALL SPLINT_PDR(X2A,YTMP,Y2TMP,N,X2,YYTMP(I))
      ENDDO

!     Construct the one-dimensional column spline and evaluate it
!     Values of 1.0D30 indicate a natural spline
      CALL SPLINE_PDR(X1A,YYTMP,M,1.0D30,1.0D30,YY2TMP)
      CALL SPLINT_PDR(X1A,YYTMP,YY2TMP,M,X1,Y)

      RETURN
    END SUBROUTINE SPLIN2
!=======================================================================



!=======================================================================
!
!  Calculate the rate of molecular hydrogen (H2) formation on grains
!  using the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222).
!
!-----------------------------------------------------------------------
FUNCTION H2_FORMATION_RATE(GAS_TEMPERATURE,GRAIN_TEMPERATURE) RESULT(RATE)

!   USE DEFINITIONS
!   USE HEALPIX_TYPES
   IMPLICIT NONE

   REAL(double) :: RATE
   REAL(double), INTENT(IN) :: GAS_TEMPERATURE,GRAIN_TEMPERATURE

   REAL(double) :: THERMAL_VELOCITY,STICKING_COEFFICIENT,TOTAL_CROSS_SECTION
   REAL(double) :: FLUX,FACTOR1,FACTOR2,EPSILON
   REAL(double) :: SILICATE_FORMATION_EFFICIENCY,GRAPHITE_FORMATION_EFFICIENCY
   REAL(double) :: SILICATE_CROSS_SECTION,SILICATE_MU,SILICATE_E_S,SILICATE_E_H2
   REAL(double) :: SILICATE_E_HP,SILICATE_E_HC,SILICATE_NU_H2,SILICATE_NU_HC
   REAL(double) :: GRAPHITE_CROSS_SECTION,GRAPHITE_MU,GRAPHITE_E_S,GRAPHITE_E_H2
   REAL(double) :: GRAPHITE_E_HP,GRAPHITE_E_HC,GRAPHITE_NU_H2,GRAPHITE_NU_HC

!  Mean thermal velocity of hydrogen atoms (cm s^-1)
   THERMAL_VELOCITY=1.45D5*SQRT(GAS_TEMPERATURE/1.0D2)

!  Calculate the thermally averaged sticking coefficient of hydrogen atoms on grains,
!  as given by Hollenbach & McKee (1979, ApJS, 41, 555, eqn 3.7)
   STICKING_COEFFICIENT=1.0D0/(1.0D0+0.04D0*SQRT(GAS_TEMPERATURE+GRAIN_TEMPERATURE) &
                    & + 0.2D0*(GAS_TEMPERATURE/1.0D2)+0.08D0*(GAS_TEMPERATURE/1.0D2)**2)

   FLUX=1.0D-10 ! Flux of H atoms in monolayers per second (mLy s^-1)

   TOTAL_CROSS_SECTION=6.273D-22 ! Total mixed grain cross section per H nucleus (cm^-2/nucleus)
   SILICATE_CROSS_SECTION=8.473D-22 ! Silicate grain cross section per H nucleus (cm^-2/nucleus)
   GRAPHITE_CROSS_SECTION=7.908D-22 ! Graphite grain cross section per H nucleus (cm^-2/nucleus)

!  Silicate grain properties
   SILICATE_MU=0.005D0   ! Fraction of newly formed H2 that stays on the grain surface
   SILICATE_E_S=110.0D0  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
   SILICATE_E_H2=320.0D0 ! Desorption energy of H2 molecules (K)
   SILICATE_E_HP=450.0D0 ! Desorption energy of physisorbed H atoms (K)
   SILICATE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
   SILICATE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
   SILICATE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)

   FACTOR1=SILICATE_MU*FLUX/(2*SILICATE_NU_H2*EXP(-SILICATE_E_H2/GRAIN_TEMPERATURE))

   FACTOR2=1.0D0*(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2 &
        & /4.0D0*EXP(-SILICATE_E_S/GRAIN_TEMPERATURE)

   EPSILON=1.0D0/(1.0D0+SILICATE_NU_HC/(2*FLUX)*EXP(-1.5*SILICATE_E_HC/GRAIN_TEMPERATURE) &
              & *(1.0D0+SQRT((SILICATE_E_HC-SILICATE_E_S)/(SILICATE_E_HP-SILICATE_E_S)))**2)

   SILICATE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON

!  Graphite grain properties
   GRAPHITE_MU=0.005D0   ! Fraction of newly formed H2 that stays on the grain surface
   GRAPHITE_E_S=260.0D0  ! Energy of the saddle point between a physisorbed and a chemisorbed site (K)
   GRAPHITE_E_H2=520.0D0 ! Desorption energy of H2 molecules (K)
   GRAPHITE_E_HP=800.0D0 ! Desorption energy of physisorbed H atoms (K)
   GRAPHITE_E_HC=3.0D4   ! Desorption energy of chemisorbed H atoms (K)
   GRAPHITE_NU_H2=3.0D12 ! Vibrational frequency of H2 molecules in surface sites (s^-1)
   GRAPHITE_NU_HC=1.3D13 ! Vibrational frequency of H atoms in their surface sites (s^-1)

   FACTOR1=GRAPHITE_MU*FLUX/(2*GRAPHITE_NU_H2*EXP(-GRAPHITE_E_H2/GRAIN_TEMPERATURE))

   FACTOR2=1.0D0*(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2 &
        & /4.0D0*EXP(-GRAPHITE_E_S/GRAIN_TEMPERATURE)

   EPSILON=1.0D0/(1.0D0+GRAPHITE_NU_HC/(2*FLUX)*EXP(-1.5*GRAPHITE_E_HC/GRAIN_TEMPERATURE) &
              & *(1.0D0+SQRT((GRAPHITE_E_HC-GRAPHITE_E_S)/(GRAPHITE_E_HP-GRAPHITE_E_S)))**2)

   GRAPHITE_FORMATION_EFFICIENCY=1.0D0/(1.0D0+FACTOR1+FACTOR2)*EPSILON

!!$!  Use the tradional rate, with a simple temperature dependence based on the
!!$!  thermal velocity of the H atoms in the gas and neglecting any temperature
!!$!  dependency of the formation and sticking efficiencies
!!$   RATE=3.0D-18*SQRT(GAS_TEMPERATURE)

!!$!  Use the treatment of de Jong (1977, A&A, 55, 137, p140 right column).
!!$!  The second exponential dependence on the gas temperature reduces the
!!$!  efficiency at high temperatures and so prevents runaway H2 formation
!!$!  heating at high temperatures:
!!$!
!!$!  k_H2 = 3E-18 * T^0.5 * exp(-T/1000)   [cm3/s]
!!$!
!!$   RATE=3.0D-18*SQRT(GAS_TEMPERATURE)*EXP(-(GAS_TEMPERATURE/1.0D3))

!!$!  Use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722, eqn 4)
!!$   RATE=0.5D0*THERMAL_VELOCITY*TOTAL_CROSS_SECTION*STICKING_COEFFICIENT

!  Use the treatment of Cazaux & Tielens (2002, ApJ, 575, L29) and
!  Cazaux & Tielens (2004, ApJ, 604, 222)
   RATE=0.5D0*THERMAL_VELOCITY*(SILICATE_CROSS_SECTION*SILICATE_FORMATION_EFFICIENCY &
    & + GRAPHITE_CROSS_SECTION*GRAPHITE_FORMATION_EFFICIENCY)*STICKING_COEFFICIENT

!!$!  Use the expression given by Markus Rollig during the February 2012 Leiden workshop
!!$   RATE=0.5D0*THERMAL_VELOCITY &
!!$      & *(SILICATE_CROSS_SECTION/((1.0D0 + 6.998D24/EXP(1.5*SILICATE_E_HC/GRAIN_TEMPERATURE)) &
!!$      & *(1.0D0 + 1.0D0/(EXP(SILICATE_E_HP/GRAIN_TEMPERATURE) &
!!$      & *(0.427D0/EXP((SILICATE_E_HP-SILICATE_E_S)/GRAIN_TEMPERATURE) + 2.5336D-14*SQRT(GRAIN_TEMPERATURE))))) &
!!$      & + GRAPHITE_CROSS_SECTION/((1.0D0 + 4.610D24/EXP(1.5*GRAPHITE_E_HC/GRAIN_TEMPERATURE)) &
!!$      & *(1.0D0 + 1.0D0/(EXP(GRAPHITE_E_HP/GRAIN_TEMPERATURE) &
!!$      & *(0.539D0/EXP((GRAPHITE_E_HP-GRAPHITE_E_S)/GRAIN_TEMPERATURE) + 5.6334D-14*SQRT(GRAIN_TEMPERATURE)))))) &
!!$      & *STICKING_COEFFICIENT

   RETURN
END FUNCTION H2_FORMATION_RATE
!=======================================================================


! Calculate the partition function for the given species
SUBROUTINE CALCULATE_PARTITION_FUNCTION(PARTITION_FUNCTION,NLEV,ENERGIES,WEIGHTS,TEMPERATURE)

 !     USE HEALPIX_TYPES
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NLEV
      REAL(double), INTENT(IN)     :: ENERGIES(1:NLEV),WEIGHTS(1:NLEV)
      REAL(double), INTENT(IN)     :: TEMPERATURE
      REAL(double), INTENT(OUT)    :: PARTITION_FUNCTION

      INTEGER :: ILEVEL

      PARTITION_FUNCTION=0.0D0
      DO ILEVEL=1,NLEV
         PARTITION_FUNCTION=PARTITION_FUNCTION + WEIGHTS(ILEVEL)*EXP(-ENERGIES(ILEVEL)/KB/TEMPERATURE)
      ENDDO

RETURN
END SUBROUTINE CALCULATE_PARTITION_FUNCTION

! Calculate the LTE level populations for the given species
SUBROUTINE CALCULATE_LTE_POPULATIONS(NLEV,LEVEL_POP,ENERGIES,WEIGHTS,PARTITION_FUNCTION,DENSITY,TEMPERATURE)

!      USE HEALPIX_TYPES
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NLEV
      REAL(double), INTENT(IN)     :: ENERGIES(1:NLEV),WEIGHTS(1:NLEV)
      REAL(double), INTENT(IN)     :: PARTITION_FUNCTION
      REAL(double), INTENT(IN)     :: DENSITY,TEMPERATURE
      REAL(double), INTENT(OUT)    :: LEVEL_POP(1:NLEV)

      INTEGER :: ILEVEL
      REAL(double) :: TOTAL_POP

      TOTAL_POP=0.0D0
      DO ILEVEL=1,NLEV
         LEVEL_POP(ILEVEL)=DENSITY*WEIGHTS(ILEVEL)*EXP(-ENERGIES(ILEVEL)/KB/TEMPERATURE)/PARTITION_FUNCTION
         TOTAL_POP=TOTAL_POP + LEVEL_POP(ILEVEL)
      ENDDO

      ! Check that the sum of the level populations adds up to the total density
      IF(ABS(TOTAL_POP-DENSITY)/DENSITY.GT.1.0D-3) THEN
         WRITE(6,*)'ERROR! Sum of LTE level populations differs from the total density by ', &
		 & INT(1.0D2*ABS(TOTAL_POP-DENSITY)/DENSITY),'%'
         STOP
      ENDIF

RETURN
END SUBROUTINE CALCULATE_LTE_POPULATIONS

     SUBROUTINE SPLINT_PDR(XA,YA,Y2A,N,X,Y)

 !     USE DEFINITIONS
!      USE HEALPIX_TYPES

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL(double), INTENT(IN)     :: XA(1:N),YA(1:N),Y2A(1:N)
      REAL(double), INTENT(IN)     :: X
      REAL(double), INTENT(OUT)    :: Y

      LOGICAL :: ASCND
      INTEGER :: JLO,JHI,JMID,INC
      REAL(double) :: A,B

      JLO=0
      JHI=0

!     ASCND is TRUE if the table values are in ascending order, FALSE otherwise
      ASCND=XA(N).GT.XA(1)

!     Find the interval XA(JLO) <= X <= XA(JLO+1) = XA(JHI)
      IF(JLO.LE.0 .OR. JLO.GT.N) THEN
!        Input guess not useful, go immediately to bisection
         JLO=0
         JHI=N+1
         GOTO 300
      ENDIF

!     Set the hunting increment
      INC=1

      IF(X.GE.XA(JLO) .EQV. ASCND) THEN
!        Hunt up:
 100     JHI=JLO+INC
         IF(JHI.GT.N) THEN
!           Done hunting, since off the end of the table
            JHI=N+1
         ELSE IF(X.GE.XA(JHI) .EQV. ASCND) THEN
!           Not done hunting...
            JLO=JHI
!           ...so double the increment...
            INC=INC+INC
!           ...and try again
            GOTO 100
         ENDIF
!     Done hunting, value bracketed
      ELSE
         JHI=JLO
!        Hunt down:
 200     JLO=JHI-INC
         IF(JLO.LT.1) THEN
!           Done hunting, since off the end of the table
            JLO=0
         ELSE IF(X.LT.XA(JLO) .EQV. ASCND) THEN
!           Not done hunting...
            JHI=JLO
!           ...so double the increment...
            INC=INC+INC
!           ...and try again
            GOTO 200
         ENDIF
!     Done hunting, value bracketed
      ENDIF

 300  IF((JHI-JLO).NE.1) THEN
!        Hunt is done, so begin the final bisection phase
         JMID=(JHI+JLO)/2
         IF(X.GT.XA(JMID) .EQV. ASCND) THEN
            JLO=JMID
         ELSE
            JHI=JMID
         ENDIF
         GOTO 300
      ENDIF

      IF(JLO.EQ.0) THEN
         JLO=1
         JHI=2
      ENDIF
      IF(JLO.EQ.N) THEN
         JLO=N-1
         JHI=N
      ENDIF

!     JLO and JHI now bracket the input value X
!     The cubic spline polynomial is now evaluated
!      print *, "JLO", JLO
!      print *, "JLI", JHI
!      print *, ""
!      print *, "XA ", XA
!      print *, ""
!      print *, "X ", X
!      print *, ""
!      print *, "A ", A
!      print *, ""
 !     print *, "(XA(JHI)-X)", (XA(JHI)-X)
 !!     print *, "TROLOL A(JHI)", XA(JHI)
 !     print *, "OMGWTFBBQ XA(JLO)", XA(JLO)
 !     print *, "(XA(JHI)-XA(JLO)", (XA(JHI)-XA(JLO))

      A=(XA(JHI)-X)/(XA(JHI)-XA(JLO))
      B=(X-XA(JLO))/(XA(JHI)-XA(JLO))
      Y=A*YA(JLO)+B*YA(JHI)+((A**3-A)*Y2A(JLO)+(B**3-B)*Y2A(JHI))&
     &  *((XA(JHI)-XA(JLO))**2)/6.0D0

      RETURN
    END SUBROUTINE SPLINT_PDR
!=======================================================================

      SUBROUTINE SPLINE_PDR(X,Y,N,YP1,YPN,Y2)

!      USE DEFINITIONS
!      USE HEALPIX_TYPES

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL(double), INTENT(IN)     :: X(1:N),Y(1:N)
      REAL(double), INTENT(IN)     :: YP1,YPN
      REAL(double), INTENT(OUT)    :: Y2(1:N)

      INTEGER :: I
      REAL(double) :: P,QN,SIG,U(1:N),UN

      IF(YP1.GE.1.0D30) THEN
!        The lower boundary condition is either set to be "natural"
         Y2(1)=0.0D0
         U(1)=0.0D0
      ELSE
!        or to have a specified first derivative
         Y2(1)=-0.5D0
         U(1)=(3.0D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

!     This is the decomposition loop of the tridiagonal algorithm
!     Y2 and U are used for temporary storage of the decomposed factors
      DO I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.0D0
         Y2(I)=(SIG-1.0D0)/P
         U(I)=(6.0D0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))&
     &              /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO

      IF(YPN.GE.1.0D30) THEN
!        The upper boundary condition is either set to be "natural"
         QN=0.0D0
         UN=0.0D0
      ELSE
!        or to have a specified first derivative
         QN=0.5D0
         UN=(3.0D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.0D0)
!     This is the back-substitution loop of the tridiagonal algorithm
      DO I=N-1,1,-1
         Y2(I)=Y2(I)*Y2(I+1)+U(I)
      ENDDO

      RETURN
    END SUBROUTINE SPLINE_PDR
!=======================================================================


end module pdr_utils_mod

#endif
