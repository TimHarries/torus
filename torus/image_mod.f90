module image_mod


  use phasematrix_mod
  use vector_mod
  use constants_mod
  use photon_mod

  implicit none

  public
  
  type IMAGETYPE
    type(STOKESVECTOR), pointer :: pixel(:,:)
    real, pointer :: vel(:,:)
    real, pointer :: totWeight(:,:)
    integer :: nsize
    real :: scale
    real :: vMin
    real :: vMax
  end type IMAGETYPE

  type PVIMAGETYPE
     real, pointer :: pixel(:,:) ! in arcsec
     real, pointer :: vAxis(:)   ! in km/s
     real, pointer :: pAxis(:)
     integer :: nv
     integer :: np
     type(VECTOR) :: slitDirection
     type(VECTOR) :: slitPosition
     real :: slitWidth, slitLength   ! in arcsec
     
  end type PVIMAGETYPE


  contains


   function initImage(nsize, rSize, vMin, vMax)

     type(IMAGETYPE) :: initImage
     integer :: nsize
     real :: scale, rSize
     real :: vmax, vmin

     scale  = rSize / real(nSize)

     allocate(initImage%pixel(-nsize:nsize,-nsize:nsize))
     allocate(initImage%vel(-nsize:nsize,-nsize:nsize))
     allocate(initImage%totWeight(-nsize:nsize,-nsize:nsize))

     initImage%pixel(-nsize:nsize,-nsize:nsize) = STOKESVECTOR(1.e-30,0.,0.,0.)
     initImage%nsize = nsize
     initImage%scale = scale
     initImage%vMin = vMin
     initImage%vMax = vMax
     initImage%vel(-nsize:nsize,-nsize:nsize) = 0.
     initImage%totWeight(-nsize:nsize,-nsize:nsize) = 0.

   end function initImage

   function initPVimage(nv, vMin, vMax, np, pMin, pMax, slitPosition, slitPA, &
                        slitWidth, slitLength)
     type(PVIMAGETYPE) :: initPVimage
     integer :: nv,np
     real :: vMin, vMax, pMin, pMax
     integer :: i 
     real :: slitPA, slitWidth, slitLength, ang
     type(VECTOR) :: slitPosition

     ang = slitPA + pi/2.
     initPVimage%slitDirection = VECTOR(cos(ang), sin(ang), 0.)
     initPVimage%slitPosition = slitPosition
     initPVimage%slitLength = slitLength
     initPVimage%slitWidth = slitWidth
     initPVimage%np = np
     initPVimage%nv = nv

     allocate(initPVimage%pixel(1:nv, 1: np))
     allocate(initPVimage%pAxis(1:np))
     allocate(initPVimage%vAxis(1:nv))
     do i = 1, np
        initPVimage%pAxis(i) = pMin + (real(i-1)/real(np-1))*(pMax-pMin)
     enddo
     do i = 1, nv
        initPVimage%vAxis(i) = vMin + (real(i-1)/real(nv-1))*(vMax-vMin)
     enddo
    initPVimage%pixel = 1.e-20
   end function initPVimage

   subroutine addPhotontoPVimage(thisImage, thisPhoton, viewVec, rotationAxis, thisVel, &
                                 weight, distance)
     type(PVIMAGETYPE) :: thisImage
     type(PHOTON) :: thisPhoton
     type(VECTOR) :: viewVec, rotationAxis
     real :: thisVel
     real :: weight
     real :: distance
     type(VECTOR) :: xProj, yProj, pVec
     type(VECTOR) :: slitNorm
     real :: velinkms
     real :: xDist, yDist
     integer :: iv, ip
     real :: xSlit, ySlit

     velInkms = thisVel * cSpeed/1.e5

     xProj =  rotationAxis .cross. viewVec
     call normalize(xProj)
     yProj = viewVec .cross. xProj
     call normalize(yProj)
     
     xDist = thisPhoton%position .dot. xProj
     yDist = thisPhoton%position .dot. yProj
     xDist = (xDist / distance) * radiansToArcsec
     yDist = (yDist / distance) * radiansToArcSec

     pVec = VECTOR(xDist,yDist,0.) - thisImage%slitPosition

     slitNorm = thisImage%slitDirection .cross. rotationAxis
     call normalize(slitNorm)

     xSlit = slitNorm .dot. pVec
     ySlit = thisImage%slitDirection .dot. pVec


     if ( (xSlit >= (-thisImage%slitWidth/2.)).and. &
          (xSlit <=  (thisImage%slitWidth/2.)) ) then
        if ( (ySlit >= (-thisImage%slitLength/2.)) .and. &
             (ySlit <= ( thisImage%slitLength/2.)) ) then
           call locate(thisImage%pAxis, thisImage%np, ySlit, ip)
           call locate(thisImage%vAxis, thisImage%nv, velinkms, iv)
           thisImage%pixel(iv, ip) = thisImage%pixel(iv,ip) + thisPhoton%stokes%i * weight
        endif
     endif
   end subroutine addPhotontoPVimage
     
     

   subroutine addPhotonToImage(viewVec, rotationAxis, thisImage, thisPhoton, &
                               thisVel, weight)
     
     type(IMAGETYPE) :: thisImage
     type(PHOTON) :: thisPhoton
     type(VECTOR) :: viewVec,  xProj, yProj, rotationAxis
     real :: xDist, yDist
     integer :: xPix, yPix
     real :: thisVel, velincgs
     real :: weight

     velIncgs = thisVel! * cSpeed/1.e5

!     write(*,*) thisvel,thisimage%vMax,thisimage%vMin

     if ((velincgs < thisImage%vMax) .and. (velincgs > thisImage%vMin)) then


        xProj =  rotationAxis .cross. viewVec
        call normalize(xProj)
        yProj = viewVec .cross. xProj
        call normalize(yProj)

        xDist = thisPhoton%position .dot. xProj
        yDist = thisPhoton%position .dot. yProj

        xPix = nint(xDist / thisImage%scale)
        yPix = nint(yDist / thisImage%scale)


        if ((xPix >= -thisImage%nSize) .and. &
             (yPix >= -thisImage%nSize) .and. &
             (xPix <= thisImage%nSize) .and. &
             (yPix <= thisImage%nSize)) then
           thisImage%pixel(xPix, yPix) = thisImage%pixel(xPix, yPix) + thisPhoton%stokes * weight
           thisImage%vel(xPix,yPix) = thisImage%vel(xPix, yPix) + velincgs * (thisPhoton%stokes%i*weight)
           thisImage%totWeight(xPix,yPix) = thisImage%totWeight(xPix,yPix) + (thisPhoton%stokes%i*weight)
        endif
     endif

     end subroutine addPhotonToImage

     subroutine freeImage(thisImage)
       type(IMAGETYPE) :: thisImage

       deallocate(thisImage%pixel)
       deallocate(thisImage%vel)
       deallocate(thisImage%totWeight)

     end subroutine freeImage

     subroutine freePVImage(thispvImage)
       type(PVIMAGETYPE) :: thispvImage

       deallocate(thisPVImage%pixel)
       deallocate(thisPVImage%paxis)
       deallocate(thisPVImage%vaxis)
     end subroutine freePVImage


     subroutine writeImage(thisImage, filename)

       type(IMAGETYPE) :: thisImage
       character(len=*) :: filename
       real, allocatable :: piximageI(:,:)
       real, allocatable :: piximageQ(:,:)
       real, allocatable :: piximageU(:,:)
       real, allocatable :: piximageV(:,:)
       real, allocatable :: piximageVel(:,:)


       integer :: nSize, nPix, i, j
       nSize = thisImage%nSize
       nPix = 2*thisImage%nSize + 1

       allocate(pixImageI(1:nPix,1:nPix))
       allocate(pixImageQ(1:nPix,1:nPix))
       allocate(pixImageU(1:nPix,1:nPix))
       allocate(pixImageV(1:nPix,1:nPix))
       allocate(pixImageVel(1:nPix,1:nPix))
       pixImageVel = 0.

       do i = 1 , nPix
          do j = 1 , nPix
             pixImageI(i,j) = thisImage%pixel(i-1-nSize,j-1-nSize)%i
             pixImageQ(i,j) = thisImage%pixel(i-1-nSize,j-1-nSize)%q
             pixImageU(i,j) = thisImage%pixel(i-1-nSize,j-1-nSize)%u
             pixImageV(i,j) = thisImage%pixel(i-1-nSize,j-1-nSize)%v
             if (thisImage%totWeight(i-1-nSize,j-1-nSize) /= 0.) then
                pixImageVel(i,j) = thisImage%vel(i-1-nSize,j-1-nSize)/ &
                     thisImage%totWeight(i-1-nSize,j-1-nSize)
             endif
          enddo
       enddo

       call writeImagef77(pixImageI, pixImageQ, pixImageU, &
            pixImageV, pixImageVel, thisImage%scale, &
            thisImage%nSize, nPix, filename)

       deallocate(pixImageI)
       deallocate(pixImageQ)
       deallocate(pixImageU)
       deallocate(pixImageV)
       deallocate(pixImageVel)
     
     end subroutine writeImage
     
     subroutine writePVimage(thisImage, filename, vSys)
       type(PVIMAGETYPE) :: thisImage
       character(len=*) :: filename
       real :: vSys
       real, allocatable :: array(:,:), xAxis(:), yAxis(:)
       integer :: nx, ny

       allocate(array(1:thisImage%nv, 1: thisImage%np))
       allocate(xAxis(1:thisImage%nv))
       allocate(yAxis(1:thisImage%np))
       
       array = log10(thisImage%pixel)
       xAxis = thisImage%vAxis + vSys
       yAxis = thisImage%pAxis
       nx = thisImage%nv
       ny = thisImage%np
       write(*,*) nx,ny,thisimage%nv,thisimage%np
       call writePVimagef77(nx, ny, array, xaxis, yaxis, filename)

       deallocate(array)
       deallocate(xAxis)
       deallocate(yAxis)

     end subroutine writePVimage

     subroutine plotSlitOnImage(thisImage, thisPVimage, device, thisDistance)
       character(len=*) :: device
       type(IMAGETYPE) :: thisImage
       type(PVIMAGETYPE) :: thisPVimage
       real, allocatable :: pixImage(:,:), axis1(:), axis2(:)
       integer :: nPix, i, j, nx, ny, nSize
       real :: fg, bg, tr(6), dt
       real :: thisDistance
       type(VECTOR) :: zAxis, slitNorm, rect(4)

       nSize = thisImage%nSize
       nPix = 2*thisImage%nSize + 1
       allocate(pixImage(1:nPix,1:nPix))
       allocate(Axis1(1:nPix))
       allocate(Axis2(1:nPix))

       zAxis = VECTOR(0.,0.,1.)
       slitNorm = thispvImage%slitDirection .cross. zAxis
       call normalize(slitNorm)


       fg = -1.e30
       bg = 1.e30
       do i = 1 , nPix
          do j = 1 , nPix
             pixImage(i,j) = thisImage%pixel(i-1-nSize,j-1-nSize)%i
             bg = min(bg, pixImage(i,j))
             fg = max(fg, pixImage(i,j))
          enddo
       enddo
       nx = nPix
       ny = nPix

      do i = -nSize, nSize
         axis1(i+nSize+1) = real(i) * thisImage%scale
         axis2(i+nSize+1) = real(i) * thisImage%scale
      enddo
      axis1(1:nPix) = axis1(1:nPix) / thisDistance * radiansToArcsec
      axis2(1:nPix) = axis2(1:nPix) / thisDistance * radiansToArcsec

      dt = (axis1(nx) - axis1(1))/real(nx-1)
      tr(1) = axis1(1) - dt
      tr(2) = dt
      tr(3) = 0.
      dt = (axis2(ny) - axis2(1))/real(ny-1)
      tr(4) = axis2(1) - dt
      tr(5) = 0.
      tr(6) = dt

      call pgbegin(0,device,1,1)
      if (device /= "/xs") then
         call pgpaper(4.,1.)
      endif
      call pgvport(0.1,0.9,0.1,0.9)
      call palette(3)
      call pgsitf(1)
      call pgwnad(axis1(1), axis1(nx), axis2(1), axis2(ny))
      call pgimag(piximage, nx, ny, 1, nx, 1, ny, bg, fg, tr)
      call pgbox('bcnst',0,0,'bcnst',0,0)
      call pglabel("Distance (arcsec)", "Distance (arcsec)", "Slit position")

      rect(1) = thisPVimage%slitPosition - (thisPVimage%slitLength/2.)* &
thisPVimage%slitDirection - (thisPVimage%slitWidth/2.)*slitnorm
      rect(2) = thisPVimage%slitPosition - (thisPVimage%slitLength/2.)* &
thisPVimage%slitDirection + (thisPVimage%slitWidth/2.)*slitnorm
      rect(3) = thisPVimage%slitPosition + (thisPVimage%slitLength/2.)* &
thisPVimage%slitDirection + (thisPVimage%slitWidth/2.)*slitnorm
      rect(4) = thisPVimage%slitPosition + (thisPVimage%slitLength/2.)* &
thisPVimage%slitDirection - (thisPVimage%slitWidth/2.)*slitnorm

      call pgsci(2)
      call pgmove(rect(1)%x, rect(1)%y)
      call pgdraw(rect(2)%x, rect(2)%y)
      call pgdraw(rect(3)%x, rect(3)%y)
      call pgdraw(rect(4)%x, rect(4)%y)
      call pgdraw(rect(1)%x, rect(1)%y)
      call pgsci(1)
      call pgend
    end subroutine plotSlitOnImage
      SUBROUTINE PALETTE(TYPE)
      INTEGER TYPE

      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
      REAL TL(4), TR(4), TG(4), TB(4)

      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/

      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/

      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/

      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/

      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
              0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
              0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
              0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      DATA TL /0.0, 0.5, 0.5, 1.0/
      DATA TR /0.2, 0.6, 0.6, 1.0/
      DATA TG /0.0, 0.0, 0.5, 1.0/
      DATA TB /1.0, 0.0, 0.0, 0.0/

      IF (TYPE.EQ.1) THEN
         CALL PGCTAB(GL, GR, GG, GB, 2, 1.0, 0.5)
      ELSE IF (TYPE.EQ.2) THEN
         CALL PGCTAB(RL, RR, RG, RB, 9, 1.0, 0.5)
      ELSE IF (TYPE.EQ.3) THEN
         CALL PGCTAB(HL, HR, HG, HB, 5, 1.0, 0.5)
      ELSE IF (TYPE.EQ.4) THEN
         CALL PGCTAB(WL, WR, WG, WB, 10, 1.0, 0.5)
      ELSE IF (TYPE.EQ.5) THEN
         CALL PGCTAB(AL, AR, AG, AB, 20, 1.0, 0.5)
      ELSE IF (TYPE.EQ.6) THEN
         CALL PGCTAB(TL, TR, TG, TB, 4, 1.0, 0.5)
      END IF
    END SUBROUTINE PALETTE

    subroutine smoothPVimage(thisImage, vSigma, pSigma)
      type(PVIMAGETYPE) :: thisImage
      real :: vSigma, pSigma
      integer :: i, j, k, ip1, ip2, ip
      integer :: iv1, iv2, iv
      real :: dv
      real, allocatable :: smoothed(:,:)
      real :: tot, dp

      allocate(smoothed(1:thisImage%nv, 1:thisImage%np))
      
      do j = 1, thisImage%nv
         ip = (5.*pSigma)/(thisImage%pAxis(2)-thisImage%pAxis(1))
         do i = 1, thisImage%np
            tot = 0.
            ip1 = max(1,i-ip)
            ip2 = min(thisImage%np, i+ip)
            do k = ip1, ip2
               dp = thisImage%pAxis(i) - thisImage%pAxis(k)
               tot = tot + thisImage%pixel(k,j) * gauss(pSigma, dp)
            enddo
            smoothed(i,j) = tot
         enddo
      enddo
      thisImage%pixel(1:thisImage%nv, 1:thisImage%np) = &
           smoothed(1:thisImage%nv, 1:thisImage%np)

      iv = (5.*vSigma)/(thisImage%vAxis(2)-thisImage%vAxis(1))
      do i = 1, thisImage%np
         do j = 1, thisImage%nv
            tot = 0.
            iv1 = max(1,j-iv)
            iv2 = min(thisImage%nv, j+iv)
            do k = iv1, iv2
               dv = thisImage%vAxis(j) - thisImage%vAxis(k)
               tot = tot + thisImage%pixel(i,k) * gauss(vSigma, dv)
            enddo
            smoothed(i,j) = tot
         enddo
      enddo
      thisImage%pixel(1:thisImage%nv, 1:thisImage%np) = &
           smoothed(1:thisImage%nv, 1:thisImage%np)



      deallocate(smoothed)
    end subroutine smoothPVimage
      
end module image_mod



