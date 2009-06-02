#ifdef USENDF

!
! ********************************************************************
!
      SUBROUTINE WRTSP(NPTS,STOKES_I,STOKES_Q,STOKES_QV,STOKES_U, &
                      STOKES_UV,LAMBDA,FILENAME)
!
! This routine reads in a polarization spectrum from a TSP format file
! and puts it into the 'current' arrays
!
      IMPLICIT NONE
      INTEGER OUT_LU
      INCLUDE '/star/include/DAT_PAR'
      INCLUDE '/star/include/SAE_PAR'
!
! The current arrays
!
      INTEGER NPTS
      REAL STOKES_I(*)
      REAL STOKES_Q(*)
      REAL STOKES_QV(*)
      REAL STOKES_U(*)
      REAL STOKES_UV(*)
      REAL LAMBDA(*)
!
! Misc.
!
      INTEGER SP
      INTEGER NDF
      CHARACTER*80 PATH,FILENAME
      CHARACTER*64 NAME
      CHARACTER*80 ERROR
      INTEGER UBND(1),LBND(1),NDFO,OAXISP,NDFQ,NDFU,QPTR
      INTEGER DIMS(10)
      INTEGER NDIMS
      INTEGER UPTR,QVPTR,UVPTR,PLACE
      INTEGER NELM,PTR1,NDF1,NDF2,NDF3,IPTR
      INTEGER PTRQ,PTRQV,PTRU,PTRUV,APTR
      CHARACTER*(DAT__SZLOC) PLOC,DLOC,QLOC,LOC      
      INTEGER STATUS

      STATUS = SAI__OK
      LBND(1) = 1
      UBND(1) = NPTS
!
! Begin the ndf and hds systems
!
      CALL NDF_BEGIN
      CALL HDS_START(STATUS)
!
! Create a new tsp ndf and map it
!
      CALL HDS_NEW(FILENAME,'OUTPUT','NDF',0,0,LOC,STATUS)
      CALL DAT_NEW(LOC,'DATA_ARRAY','_REAL',1,UBND,STATUS)
      CALL NDF_IMPRT(LOC,NDFO,STATUS)
      CALL NDF_ACRE(NDFO,STATUS)
      CALL NDF_ACPUT('Wavelength',NDFO,'Lab',1,STATUS)
      CALL NDF_ACPUT('Angstrom',NDFO,'Unit',1,STATUS)
      CALL NDF_AMAP(NDFO,'Centre',1,'_REAL','UPDATE',OAXISP,UBND,STATUS)

      CALL NDF_XNEW(NDFO,'POLARIMETRY','EXT',0,0,PLOC,STATUS)
      CALL NDF_PLACE(PLOC,'STOKES_Q',PLACE,STATUS)
      CALL NDF_NEW('_REAL',1,LBND,UBND,PLACE,NDFQ,STATUS)
      CALL NDF_PLACE(PLOC,'STOKES_U',PLACE,STATUS)
      CALL NDF_NEW('_REAL',1,LBND,UBND,PLACE,NDFU,STATUS)

      CALL NDF_MAP(NDFO,'DATA','_REAL','WRITE',IPTR,UBND,STATUS)
      CALL NDF_MAP(NDFQ,'DATA','_REAL','WRITE',QPTR,UBND,STATUS)
      CALL NDF_MAP(NDFU,'DATA','_REAL','WRITE',UPTR,UBND,STATUS)
      CALL NDF_MAP(NDFQ,'VARIANCE','_REAL','WRITE',QVPTR,UBND,STATUS)
      CALL NDF_MAP(NDFU,'VARIANCE','_REAL','WRITE',UVPTR,UBND,STATUS)
!
! If everything is oK then write out the arrays
!
      IF (STATUS.EQ.SAI__OK) THEN
      CALL WRITE_IT(NPTS,LAMBDA,STOKES_I,STOKES_Q,STOKES_QV,STOKES_U, &
                    STOKES_UV, &
      UBND(1),%VAL(IPTR),%VAL(QPTR),%VAL(QVPTR),%VAL(UPTR), &
      %VAL(UVPTR),%VAL(OAXISP))
      ENDIF
!
! Close down the ndf and hds
!
      CALL DAT_ANNUL(PLOC,STATUS)
      CALL NDF_END(STATUS)
      CALL HDS_CLOSE(LOC,STATUS)
      CALL HDS_STOP(STATUS)
    END SUBROUTINE WRTSP


!
! ********************************************************************
!

      SUBROUTINE WRITE_IT(NPTS,LAMBDA,STOKES_I,STOKES_Q,STOKES_QV, &
      STOKES_U,STOKES_UV,N,IAR,Q,QV,U,UV,WA)
!
! Writes out the array mapped previously
!
      INTEGER NPTS
      REAL STOKES_I(*)
      REAL STOKES_Q(*)
      REAL STOKES_QV(*)
      REAL STOKES_U(*)
      REAL STOKES_UV(*)
      REAL LAMBDA(*)
      INTEGER STATUS,N,I
      CHARACTER*80 INFILE
      REAL Q(N),U(N),QV(N),UV(N),WA(N),IAR(N)

      DO I = 1,N
       IAR(I) = STOKES_I(I)
       WA(I) = LAMBDA(I)
       Q(I) = STOKES_Q(I)
       QV(I) = STOKES_QV(I)
       U(I) = STOKES_U(I)
       UV(I) = STOKES_UV(I)
      ENDDO
    END SUBROUTINE WRITE_IT

#else

!
! Stub routine to use when NDF support is not required. 
!
      SUBROUTINE WRTSP(NPTS,STOKES_I,STOKES_Q,STOKES_QV,STOKES_U,STOKES_UV,LAMBDA,FILENAME)

      IMPLICIT NONE

      INTEGER NPTS
      REAL STOKES_I(*)
      REAL STOKES_Q(*)
      REAL STOKES_QV(*)
      REAL STOKES_U(*)
      REAL STOKES_UV(*)
      REAL LAMBDA(*)

      CHARACTER*80 FILENAME

    END SUBROUTINE WRTSP

#endif

