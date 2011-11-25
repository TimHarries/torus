  type PVIMAGETYPE
     real, pointer :: pixel(:,:) => null()! in arcsec
     real, pointer :: vAxis(:) => null()  ! in km/s
     real, pointer :: pAxis(:) => null()
     integer :: nv
     integer :: np
     type(VECTOR) :: slitDirection
     type(VECTOR) :: slitPosition
     real :: slitWidth, slitLength   ! in arcsec
     
  end type PVIMAGETYPE

   function initPVimage(nv, vMin, vMax, np, pMin, pMax, slitPosition, slitPA, &
                        slitWidth, slitLength)
     type(PVIMAGETYPE) :: initPVimage
     integer :: nv,np
     real :: vMin, vMax, pMin, pMax
     integer :: i 
     real :: slitPA, slitWidth, slitLength, ang
     type(VECTOR) :: slitPosition

     ang = real(slitPA + pi/2.)
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
     use utils_mod, only: locate
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

     velInkms = real(thisVel * cSpeed/1.e5)

     xProj =  rotationAxis .cross. viewVec
     call normalize(xProj)
     yProj = viewVec .cross. xProj
     call normalize(yProj)
     
     xDist = real(thisPhoton%position .dot. xProj)
     yDist = real(thisPhoton%position .dot. yProj)
     xDist = real((xDist / distance) * radiansToArcsec)
     yDist = real((yDist / distance) * radiansToArcSec)

     pVec = VECTOR(xDist,yDist,0.) - thisImage%slitPosition

     slitNorm = thisImage%slitDirection .cross. rotationAxis
     call normalize(slitNorm)

     xSlit = real(slitNorm .dot. pVec)
     ySlit = real(thisImage%slitDirection .dot. pVec)


     if ( (xSlit >= (-thisImage%slitWidth/2.)).and. &
          (xSlit <=  (thisImage%slitWidth/2.)) ) then
        if ( (ySlit >= (-thisImage%slitLength/2.)) .and. &
             (ySlit <= ( thisImage%slitLength/2.)) ) then
           call locate(thisImage%pAxis, thisImage%np, ySlit, ip)
           call locate(thisImage%vAxis, thisImage%nv, velinkms, iv)
           thisImage%pixel(iv, ip) = real(thisImage%pixel(iv,ip) + thisPhoton%stokes%i * weight)
        endif
     endif
   end subroutine addPhotontoPVimage

     subroutine freePVImage(thispvImage)
       type(PVIMAGETYPE) :: thispvImage

       deallocate(thisPVImage%pixel)
       deallocate(thisPVImage%paxis)
       deallocate(thisPVImage%vaxis)
     end subroutine freePVImage

    subroutine smoothPVimage(thisImage, vSigma, pSigma)
      use utils_mod, only: gauss
      type(PVIMAGETYPE) :: thisImage
      real :: vSigma, pSigma
      integer :: i, j, k, ip1, ip2, ip
      integer :: iv1, iv2, iv
      real :: dv
      real, allocatable :: smoothed(:,:)
      real :: tot, dp

      allocate(smoothed(1:thisImage%nv, 1:thisImage%np))
      
      do j = 1, thisImage%nv
         ip = int((5.*pSigma)/(thisImage%pAxis(2)-thisImage%pAxis(1)))
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

      iv = int((5.*vSigma)/(thisImage%vAxis(2)-thisImage%vAxis(1)))
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

     subroutine ConvertImageToMJanskiesPerStr(image, lambda, distance)
       use utils_mod, only: convertToJanskies
       type(IMAGETYPE) :: image
       real(double) :: lambda
       real(double) :: distance, scale
       real(double) :: angularScale, strad
       integer :: i, j

       scale = 1.d20 / distance**2 ! to per cm^2
       
       angularScale =  (image%xAxisCentre(2) - image%xAxisCentre(1))*1.d10/distance
       strad = angularScale**2 ! str per pix

       if ( strad == 0.0 ) then 
          call writewarning( "ConvertImageToMJanskiesPerStr: strad = 0.0 no conversion performed")
          return
       end if

       do i = 1, image%nx
          do j = 1, image%ny 

! factor 1.d-6 to convert janskies to MJanskies
             image%pixel(i,j)%i = 1.d-6*convertToJanskies(dble(scale*image%pixel(i,j)%i),lambda)/strad
             image%pixel(i,j)%q = 1.d-6*convertToJanskies(dble(scale*image%pixel(i,j)%q),lambda)/strad
             image%pixel(i,j)%u = 1.d-6*convertToJanskies(dble(scale*image%pixel(i,j)%u),lambda)/strad
             image%pixel(i,j)%v = 1.d-6*convertToJanskies(dble(scale*Image%pixel(i,j)%v),lambda)/strad
          enddo
       enddo
     end subroutine ConvertImageToMJanskiesPerStr

