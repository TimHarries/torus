module pixplot_module
  ! based on  PROGRAM PGDEM4
  ! contains routines to plot the 2D pix map using  PGIMAG and etc.
  ! ====>  Need to be linked with PGPLOT library. <=====
  !
  !-----------------------------------------------------------------------
  ! Uses imaging routine PGIMAG and associated
  ! routines PGWEDG and PGCTAB.
  !-----------------------------------------------------------------------


  contains

    
    ! This is the main routine in this module.
    
    subroutine pixplot(F, MXI, MXJ, I1, I2, J1, J2, title)
            
      !
      !Arguments:
      ! F      (input), real     : the array to be plotted.
      ! MXI    (input), integer  : the first dimension of array A.
      ! MXJ    (input), integer  : the second dimension of array A.
      ! I1, I2 (input), integer  : the inclusive range of the first index
      !                            (I) to be plotted.
      ! J1, J2 (input), integer  : the inclusive range of the second
      !                            index (J) to be plotted.
      ! title  (input), character: the title of the plot

      implicit none
      INTEGER :: PGOPEN
      INTEGER ::  MXI, MXJ
      integer :: I1, I2, J1, J2
      INTEGER :: I, L, C1, C2, NC
      REAL :: F(MXI,MXJ)
      REAL :: FMIN,FMAX,TR(6), CONTRA, BRIGHT, ANGLE, C, S, ALEV(1)
      CHARACTER*16 :: VAL
      character (LEN=*) :: title
      real :: offset
      !
      ! Introduction.
      !
      WRITE(*,*)'PIXPLOT:: Plots 2D pixcel map.'
      WRITE(*,*)'This program requires a device with color capability.'
      WRITE(*,*)'On an interactive device, you can modify the color map'
      WRITE(*,*)'used for the image.'
      WRITE(*,*)
      !
      ! Open device for graphics.
      !
      IF (PGOPEN('?') .LT. 1) STOP
      CALL PGQINF('TYPE', VAL, L)
      WRITE (*,*) 'PGPLOT device type: ', VAL(1:L)
      CALL PGQCIR(C1, C2)
      NC = MAX(0, C2-C1+1)
      WRITE (*,*) 'Number of color indices used for image: ', NC
      IF (NC .LT.8) THEN 
         WRITE (*,*) 'Not enough colors available on this device'
         STOP
      ELSE
         WRITE (*,*)
      END IF
      !
      !-----------------------------------------------------------------------
      ! Example 1: simple transformation matrix
      !-----------------------------------------------------------------------
      !
      ! Set the coordinate transformation matrix: 
      ! world coordinate = pixel number.
      !
      TR(1) = 0.0
      TR(2) = 1.0
      TR(3) = 0.0
      TR(4) = 0.0
      TR(5) = 0.0
      TR(6) = 1.0
      !
      ! Clear the screen. Set up window and viewport.
      !
      CALL PGPAGE
      CALL SETVP
      CALL PGWNAD(0.0, 1.0+MXI, 0.0, 1.0+MXJ)
      !
      ! Set up the color map.
      !
      BRIGHT = 0.5
      CONTRA  = 1.0
      CALL PALETT(2, CONTRA, BRIGHT)
      !
      ! Draw the map with PGIMAG.  
      !
      
      offset = 0.0
      FMIN = MINVAL(F)
      FMAX = MAXVAL(F) - offset
      CALL PGIMAG(F,MXI,MXJ,I1,I2,J1,J2,FMIN,FMAX,TR)
      !
      ! Annotate the plot.
      !
      CALL PGMTXT('t',1.0,0.0,0.0,title)
      CALL PGSCH(0.6)
      CALL PGBOX('bcntsi',0.0,0,'bcntsiv',0.0,0)
      CALL PGMTXT('b',3.0,1.0,1.0,'pixel number')
      !
      ! Draw a wedge.
      !
      CALL PGWEDG('BI', 4.0, 5.0, FMIN, FMAX, 'pixel value')
      CALL PGSCH(1.0)
      !
      ! If the device has a cursor, allow user to fiddle with color table.
      !
      CALL PGQINF('CURSOR', VAL, L)
      IF (VAL(:L).EQ.'YES') THEN
         CALL FIDDLE
         CALL PGASK(.FALSE.)
      END IF
      !
      !-----------------------------------------------------------------------
      ! Example 2: rotation, overlay contours.
      !-----------------------------------------------------------------------
      !
      ! Compute the coordinate transformation matrix. The matrix is chosen
      ! to put array element (MXI/2, MXJ/2) at (X,Y)=(0,0), and map the
      ! entire array onto a square of side 2, rotated through angle ANGLE
      ! radians.
      !

      !ANGLE = 120.0/57.29578

      ANGLE = 0.0

      C = COS(ANGLE)
      S = SIN(ANGLE)
      TR(1) = -C - S
      TR(2) = 2.0*C/REAL(MXI)
      TR(3) = 2.0*S/REAL(MXJ)
      TR(4) = -C + S
      TR(5) = (-2.0)*S/REAL(MXI)
      TR(6) = 2.0*C/REAL(MXJ)
      !
      ! Clear the screen. Set up window and viewport.
      !
      CALL PGPAGE
      CALL SETVP
      CALL PGWNAD(-1.0, 1.0, -1.0, 1.0)
      CALL PGSCI(1)
      !
      ! Set up the color map.
      !
      BRIGHT = 0.5
      CONTRA  = 1.0
      CALL PALETT(2, CONTRA, BRIGHT)
      !
      ! Draw the map with PGIMAG.  
      !
      CALL PGIMAG(F,MXI,MXJ,I1,I2,J1,J2,FMIN,FMAX,TR)
      !
      ! Overlay contours in white.
      !
      CALL PGSCI(1)
      DO  I=1,21
         ALEV(1) = FMIN + (I-1)*(FMAX-FMIN)/20.0
         IF (MOD(I,5).EQ.0) THEN
            CALL PGSLW(3)
         ELSE
            CALL PGSLW(1)
         END IF
         IF (I.LT.10) THEN
            CALL PGSLS(2)
         ELSE
            CALL PGSLS(1)
         END IF
         CALL PGCONT(F,MXI,MXJ,I1,I2,J1,J2,ALEV,-1,TR)
      end do

      CALL PGSLS(1)
      CALL PGSLW(1)
      !         
      ! Annotate the plot.
      !
      CALL PGSCI(1)
      CALL OUTLIN(I1,I2,J1,J2,TR)
      CALL PGMTXT('t',1.0,0.0,0.0,title)
      CALL PGSCH(0.6)
      CALL PGBOX('bctsn',0.0,0,'bctsn',0.0,0)
      !
      ! Draw a wedge.
      !
      CALL PGWEDG('BI', 4.0, 5.0, FMIN, FMAX, 'pixel value')
      CALL PGSCH(1.0)
      !
      ! If the device has a cursor, allow user to fiddle with color table.
      !
      CALL PGQINF('CURSOR', VAL, L)
      IF (VAL(:L).EQ.'YES') THEN
         CALL FIDDLE
      END IF
      !
      ! Close the device and exit.
      !
      CALL PGEND
      !-----------------------------------------------------------------------
    END subroutine pixplot

    !
    !
    !

    SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
      !--------------------------------------------------------------------
      ! Set a "palette" of colors in the range of color indices used by
      ! PGIMAG.
      !--------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
      !
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
      !
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
      !
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
      !
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
      !
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
      !
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
           0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,   &
           0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,   &
           0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,   &
           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
      !
      IF (TYPE.EQ.1) THEN
         !        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
         !        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
         !        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
         !        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
         !        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
    END SUBROUTINE PALETT
    
    SUBROUTINE SETVP
      !-----------------------------------------------------------------------
      ! Set the viewport, allowing margins around the edge for annotation.
      ! (This is similar in effect to PGVSTD, but has different margins.)
      ! The routine determines the view-surface size and allocates margins
      ! as fractions of the minimum of width and height.
      !-----------------------------------------------------------------------
      REAL D, VPX1, VPX2, VPY1, VPY2
      !
      CALL PGSVP(0.0, 1.0, 0.0, 1.0)
      CALL PGQVP(1, VPX1, VPX2, VPY1, VPY2)
      D = MIN(VPX2-VPX1, VPY2-VPY1)/40.0
      VPX1 = VPX1 + 5.0*D
      VPX2 = VPX2 - 2.0*D
      VPY1 = VPY1 + 8.0*D
      VPY2 = VPY2 - 2.0*D
      CALL PGVSIZ(VPX1, VPX2, VPY1, VPY2)
    END SUBROUTINE SETVP

    SUBROUTINE FIDDLE
      !
      INTEGER P, IER, PGCURS
      REAL CONTRA, BRIGHT, X, Y, SIGN
      REAL X1, Y1, X2, Y2, B1, B2, C1, C2
      CHARACTER CH
      !
      WRITE (*,*) 'Use cursor to adjust color table:'
      WRITE (*,*) ' Keys 1,2,3,4,5 select different palettes'
      WRITE (*,*) ' Key P cycles through available palettes'
      WRITE (*,*) ' Key F adjusts contrast and brightness, with'
      WRITE (*,*) '  cursor x position setting brightness [0.0 - 1.0]'
      WRITE (*,*) '   and y position setting contrast [0.0 - 10.0]'
      WRITE (*,*) '  (Hold down F key while moving cursor to change'
      WRITE (*,*) '  contrast and brightness continuously)'
      WRITE (*,*) ' Key C resets contrast=1.0, brightness=0.5'
      WRITE (*,*) ' Key - reverses color palette'
      WRITE (*,*) ' Key X or right mouse button exits program' 
      !
      P = 2
      CONTRA = 1.0
      BRIGHT = 0.5
      X = 0.5
      Y = 1.0
      SIGN = +1.0
      !
      CALL PGQWIN(X1, X2, Y1, Y2)
      B1 = 0.0
      B2 = 1.0
      C1 = 0.0
      C2 = 10.0
      CALL PGSWIN(B1, B2, C1, C2)
10    IER = PGCURS(X, Y, CH)
      IF (CH.EQ.CHAR(0) .OR. CH.EQ.'x' .OR. CH.EQ.'X') THEN
         CALL PGSWIN(X1, X2, Y1, Y2)
         RETURN
      ELSE IF (CH.EQ.'F' .OR. CH.EQ.'f') THEN
         BRIGHT = MAX(B1, MIN(B2,X))
         CONTRA = MAX(C1, MIN(C2,Y))
      ELSE IF (CH.EQ.'C' .OR. CH.EQ.'c') THEN
         CONTRA = 1.0
         Y = 1.0
         BRIGHT = 0.5
         X = 0.5
      ELSE IF (CH.EQ.'-') THEN
         SIGN = -SIGN
      ELSE IF (CH.EQ.'1') THEN
         P = 1
      ELSE IF (CH.EQ.'2') THEN
         P = 2
      ELSE IF (CH.EQ.'3') THEN
         P = 3
      ELSE IF (CH.EQ.'4') THEN
         P = 4
      ELSE IF (CH.EQ.'5') THEN
         P = 5
      ELSE IF (CH.EQ.'P' .OR. CH.EQ.'p') THEN
         P = 1 + MOD(P,5)
      END IF
      CALL PALETT(P, SIGN*CONTRA, BRIGHT)
      GOTO 10
    END SUBROUTINE FIDDLE


    SUBROUTINE OUTLIN(I1,I2,J1,J2,TR)
      INTEGER I1,I2,J1,J2
      REAL TR(6)
      !-----------------------------------------------------------------------
      ! Draw the enclosing rectangle of the subarray to be contoured,
      ! applying the transformation TR.
      !
      ! For a contour map, the corners are (I1,J1) and (I2,J2); for
      ! a gray-scale map, they are (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).
      !-----------------------------------------------------------------------
      INTEGER K
      REAL XW(5), YW(5), T
      !
      XW(1) = I1
      YW(1) = J1
      XW(2) = I1
      YW(2) = J2
      XW(3) = I2
      YW(3) = J2
      XW(4) = I2
      YW(4) = J1
      XW(5) = I1
      YW(5) = J1
      DO  K=1,5
         T = XW(K)
         XW(K) = TR(1) + TR(2)*T + TR(3)*YW(K)
         YW(K) = TR(4) + TR(5)*T + TR(6)*YW(K)
      end DO
      CALL PGLINE(5,XW,YW)
    END SUBROUTINE OUTLIN



  end module pixplot_module


