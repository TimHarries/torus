module image_mod


  use phasematrix_mod
  use vector_mod
  use constants_mod
  use photon_mod
  use filter_set_class

  implicit none

  public
  
  type IMAGETYPE
    type(STOKESVECTOR), pointer :: pixel(:,:) => null()
    real, pointer :: vel(:,:) => null()
    real, pointer :: totWeight(:,:) => null()
    real, pointer :: xAxisCentre(:) => null()
    real, pointer :: yAxisCentre(:) => null()
    real, pointer :: xAxisLH(:) => null()
    real, pointer :: yAxisBottom(:) => null()
    integer :: nx
    integer :: ny
    real :: dx, dy
    real :: vMin
    real :: vMax
  end type IMAGETYPE

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


  contains


   function initImage(nx, ny, imageSizeX, imageSizeY, vMin, vMax)

     type(IMAGETYPE) :: initImage
     integer :: nx, ny
     real :: imagesizeX, imageSizeY
     real :: vmax, vmin
     integer :: i

     allocate(initImage%pixel(1:nx,1:ny))
     allocate(initImage%vel(1:nx,1:ny))
     allocate(initImage%totWeight(1:nx,1:ny))

     allocate(initImage%xAxisCentre(1:nx))
     allocate(initImage%yAxisCentre(1:ny))
     allocate(initImage%xAxisLH(1:nx))
     allocate(initImage%yAxisBottom(1:ny))

     initimage%dx = imageSizeX / real(nx)
     initimage%dy = imageSizeY / real(ny)
     
     do i = 1, nx
        initImage%xAxisCentre(i) = -imageSizeX/2. + initimage%dx/2. + initimage%dx*real(i-1)
        initImage%xAxisLH(i) = -imageSizeX/2. + initimage%dx*real(i-1)
     enddo
     do i = 1, ny
        initImage%yAxisCentre(i) = -imageSizeY/2. + initimage%dy/2. + initimage%dy*real(i-1)
        initImage%yAxisBottom(i) = -imageSizeY/2. + initimage%dy*real(i-1)
     enddo

     initImage%pixel = STOKESVECTOR(0.,0.,0.,0.)
     initImage%nx = nx
     initImage%ny = ny

     initImage%vMin = vMin
     initImage%vMax = vMax
     initImage%vel = 0.
     initImage%totWeight = 0.



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
     
     

   subroutine addPhotonToImage(viewVec, rotationAxis, thisImageSet, nImage, thisPhoton, &
                               thisVel, weight, filters, center, lambda0_cont)
     
     integer, intent(in) :: nImage  ! number of images in a set
     type(IMAGETYPE), intent(inout) :: thisImageSet(nImage)
     type(PHOTON) :: thisPhoton
     type(VECTOR) :: viewVec,  xProj, yProj, rotationAxis
     real :: xDist, yDist
     integer :: xPix, yPix
     real :: thisVel, velincgs
     real :: weight
     type(filter_set), intent(in) :: filters     
     type(VECTOR), intent(in) :: center  ! the center of the model space. [10^10cm]
     real, intent(in), optional :: lambda0_cont  ! rest wavelength of contiuum photon
     !
     integer :: i
     real :: filter_response
     real(double) :: lambda_obs

     velIncgs = thisVel! * cSpeed/1.e5 


     ! observed wavelength should be Doppler shifted by local gas velocity
     if (thisPhoton%contPhoton) then
        if (PRESENT(lambda0_cont)) then
           lambda_obs = dble( lambda0_cont*(1.0d0+thisVel) )
        else
           lambda_obs = dble( thisPhoton%lambda*(1.0d0+thisVel) )
        end if
     else
        lambda_obs = dble( thisPhoton%lambda*(1.0d0+thisVel) )  
     end if

!     write(*,*) thisvel,thisimage%vMax,thisimage%vMin


     do i = 1, nImage
!        if ((velincgs < thisImageSet(i)%vMax) .and. (velincgs > thisImageSet(i)%vMin)) then
           xProj =  rotationAxis .cross. viewVec
           call normalize(xProj)
           yProj = viewVec .cross. xProj
           call normalize(yProj)
           
!           xDist = (thisPhoton%position - center) .dot. xProj
!           yDist = (thisPhoton%position - center) .dot. yProj
           xDist = (thisPhoton%position) .dot. xProj
           yDist = (thisPhoton%position) .dot. yProj
           

           call pixelLocate(thisImageSet(i), xDist, yDist, xPix, yPix)

           if ((xPix >= 1) .and. &
                (yPix >= 1) .and. &
                (xPix <= thisImageSet(i)%nx) .and. &
                (yPix <= thisImageSet(i)%ny)) then
              ! using a filter response function in filter_set_class here
              filter_response = real( pass_a_filter(filters, i, lambda_obs ))
              
              
              thisImageSet(i)%pixel(xPix, yPix) = thisImageSet(i)%pixel(xPix, yPix)  &
                   + thisPhoton%stokes * weight * filter_response
              thisImageSet(i)%vel(xPix,yPix) = thisImageSet(i)%vel(xPix, yPix)  &
                   + velincgs * (thisPhoton%stokes%i*weight*filter_response)
              thisImageSet(i)%totWeight(xPix,yPix) = thisImageSet(i)%totWeight(xPix,yPix) &
                   + (thisPhoton%stokes%i*weight*filter_response)

           endif
!        endif
     end do

     end subroutine addPhotonToImage

     subroutine freeImage(thisImage)
       type(IMAGETYPE) :: thisImage

       deallocate(thisImage%pixel)
       deallocate(thisImage%vel)
       deallocate(thisImage%totWeight)
       NULLIFY(thisImage%pixel)
       NULLIFY(thisImage%vel)
       NULLIFY(thisImage%totWeight)

     end subroutine freeImage

     subroutine freePVImage(thispvImage)
       type(PVIMAGETYPE) :: thispvImage

       deallocate(thisPVImage%pixel)
       deallocate(thisPVImage%paxis)
       deallocate(thisPVImage%vaxis)
     end subroutine freePVImage


     subroutine writeImage(thisImage, filename, objectDistance, inArcsec, lambda_center, bandwidth)
       
       type(IMAGETYPE),intent(in) :: thisImage
       character(len=*), intent(in) :: filename
       real(double), intent(in) :: objectDistance ! in [cm]
       logical, intent(in)      :: inArcsec       ! if T, the dimension of images are in arcsec
       real(double), intent(in) :: lambda_center ! of the filter in [A]
       real(double), intent(in) :: bandwidth     ! of the filter in [A]
       !
       real, allocatable :: piximageI(:,:)
       real, allocatable :: piximageQ(:,:)
       real, allocatable :: piximageU(:,:)
       real, allocatable :: piximageV(:,:)
       real, allocatable :: piximageVel(:,:)
       real(double) :: area
       real :: imScale
       real :: physicalPixelSize
       real :: angularPixelSize
!       real :: pixelArea
       real(double) :: scale, convert, c_in_m_per_sec

       integer :: nPix, i, j !, nSize

!       nSize = thisImage%nSize
!       nPix = 2*thisImage%nSize + 1

       npix = thisImage%nx ! assumes square image!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       physicalPixelSize = 1.e10 * (thisImage%dx * thisImage%dy)
       angularPixelSize = physicalPixelSize / objectDistance
!       pixelArea = (angularPixelSize * radtoDeg * 3600.)**2 ! square arcsec

       area = objectDistance**2

       allocate(pixImageI(1:nPix,1:nPix))
       allocate(pixImageQ(1:nPix,1:nPix))
       allocate(pixImageU(1:nPix,1:nPix))
       allocate(pixImageV(1:nPix,1:nPix))
       allocate(pixImageVel(1:nPix,1:nPix))
       pixImageVel = 0.


       ! Note: stokes parameters are in erg/s/pix^2/sr, and we now converting 
       ! them to Stokes Flux using a right conversion factor.

       physicalPixelSize = thisImage%dx * 1.e10
       angularPixelSize = physicalPixelSize / objectDistance
!       pixelArea = (angularPixelSize * radtoDeg * 3600.)**2 ! square arcsec
       
       ! 
       area =(thisImage%dx)**2 * 1.0d20  ! cm^2

       scale = ( (thisImage%dx*1.e10)/objectDistance )**2  
       ! Conversion factor (1 Jy = 9.57e-13 erg/s/cm^2/A at 5600 A)
       ! and the fact that 
       !             F_lambda[erg/s/cm^2/A] = (c/lambda^2) * F_nu [erg/s/cm^2/Hz] 
       !                                    = (c/lambda^2) * 1.0e-23 F_nu[Jy]
       c_in_m_per_sec = 2.998d8 ! [m/s]
       convert = c_in_m_per_sec/((lambda_center*1.d-10)**2)  * 1.0d-23 * 1.0d-10 
       ! --- should be in per angstrome now.

       convert = 1.0e20/convert/area/bandwidth ! [Jy/(erg/s/cm^2/A)/A] = [Jy/(erg/s/cm^2)]
       ! -- NB: 1.0e20 factor is the offset made in energyPerPhoton in main routine!

       scale = scale * convert

       scale = 1.d0 ! energyperphoton scaling factor

       write(*,*) "OBJECT DISTANCE ",objectDistance
       scale = scale * (1.d0 / objectDistance**2) ! erg/s/pix/cm^2

       scale = scale * 1.d-4 ! erg/s/pix/m^2

       scale = scale * 1.d-7 ! W/m^2/pix

       write(*,*) "Scaling factor ", scale
!       !
!       write(*,*) " "
!       write(*,*) "================================================================="
!       write(*,*) "Image Scale Factor =", scale, " [Jy/(erg/s)*sr*(cm^2*A)/(cm^2*A)]"
!       write(*,*) "================================================================="
!       write(*,*) " "
!       !

       do i = 1 , nPix
          do j = 1 , nPix
             pixImageI(i,j) = thisImage%pixel(i,j)%i * scale ![Jy]
             pixImageQ(i,j) = thisImage%pixel(i,j)%q * scale ![Jy]
             pixImageU(i,j) = thisImage%pixel(i,j)%u * scale ![Jy]
             pixImageV(i,j) = thisImage%pixel(i,j)%v * scale ![Jy]
             ! Now the images are in Janskies!

             if (thisImage%totWeight(i,j) /= 0.) then
                pixImageVel(i,j) = thisImage%vel(i,j)/ &
                     thisImage%totWeight(i,j)
             endif
          enddo
       enddo

       if (.not.inArcsec) then
          imscale = thisImage%dx
       else
          imscale = ((thisImage%dx * 1.e10) / objectDistance) * radtodeg * 3600.
       endif

       call writeImagef77(pixImageI, pixImageQ, pixImageU, &
            pixImageV, pixImageVel, imscale, nPix, filename)

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


    subroutine simply_writeImage(thisImage, filename)
       
       type(IMAGETYPE),intent(in) :: thisImage
       character(len=*), intent(in) :: filename
       !
       real, allocatable :: piximageI(:,:)
       real, allocatable :: piximageQ(:,:)
       real, allocatable :: piximageU(:,:)
       real, allocatable :: piximageV(:,:)
       real, allocatable :: piximageVel(:,:)
       real :: imScale

       integer :: nPix, i, j

       nPix = thisImage%nx


       allocate(pixImageI(1:nPix,1:nPix))
       allocate(pixImageQ(1:nPix,1:nPix))
       allocate(pixImageU(1:nPix,1:nPix))
       allocate(pixImageV(1:nPix,1:nPix))
       allocate(pixImageVel(1:nPix,1:nPix))
       pixImageVel = 0.


       

       do i = 1 , nPix
          do j = 1 , nPix
             pixImageI(i,j) = thisImage%pixel(i,j)%i
             pixImageQ(i,j) = thisImage%pixel(i,j)%q
             pixImageU(i,j) = thisImage%pixel(i,j)%u
             pixImageV(i,j) = thisImage%pixel(i,j)%v
          enddo
       enddo

       imscale = thisImage%dx
       
       call writeImagef77(pixImageI, pixImageQ, pixImageU, &
            pixImageV, pixImageVel, imscale,  nPix, filename)

       deallocate(pixImageI)
       deallocate(pixImageQ)
       deallocate(pixImageU)
       deallocate(pixImageV)
       deallocate(pixImageVel)
     
     end subroutine simply_writeImage

     subroutine writePPMAimage(tfile, rimage, gimage, bimage, nx, ny)
       character(len=*) :: tfile
       integer :: nx, ny
       integer :: rimage(nx,ny), gimage(nx,ny), bimage(nx, ny)
       integer :: i, j

       open(33, file=tfile, form="formatted", status="unknown")
       write(33,'(a,i3.3,1x,i3.3,1x,i3.3)') "P3 ", nx, ny, 255
       do j = 1, ny
          do i = 1, nx
             write(33, '(i3.3,1x,i3.3,1x,i3.3)') rimage(i,j), gimage(i,j), bimage(i,j)
          enddo
       enddo
       close(33)
     end subroutine writePPMAimage

     subroutine writeFalseColourPPM(tfile, image, nImage)
       character(len=*) :: tfile
       integer :: nImage
       type(IMAGETYPE) :: image(:)
       integer, allocatable :: rImage(:,:), gImage(:,:), bImage(:,:)
       integer :: nx, ny
       real(double) :: f5, f95
       real :: t
       integer :: i,j

       nx = image(1)%nx
       ny = image(1)%ny
       allocate(rImage(nx, ny))
       allocate(gImage(nx, ny))
       allocate(bImage(nx, ny))
       

       f5 = imagePercentile(image(1), 5.)
       f95 = imagePercentile(image(1), 95.)
       write(*,*) f5, f95
       do i = 1, nx
          do j = 1, ny
             t = (min(f95, max(f5, image(1)%pixel(i,j)%i))-f5)/(f95-f5)
             rImage(i, j) = int(255.*t)
          enddo
       enddo

       f5 = imagePercentile(image(2), 5.)
       f95 = imagePercentile(image(2), 95.)
       write(*,*) f5, f95
       do i = 1, nx
          do j = 1, ny
             t = (min(f95, max(f5, image(2)%pixel(i,j)%i))-f5)/(f95-f5)
             gImage(i, j) = int(255.*t)
          enddo
       enddo

       f5 = imagePercentile(image(3), 5.)
       f95 = imagePercentile(image(3), 95.)
       write(*,*) f5, f95
       do i = 1, nx
          do j = 1, ny
             t = (min(f95, max(f5, image(3)%pixel(i,j)%i))-f5)/(f95-f5)
             bImage(i, j) = int(255.*t)
          enddo
       enddo


       write(*,*) "calling writepmmaimage"
       call writePPMAimage(tfile, rImage, gImage, bImage, nx, ny)
       write(*,*) "done"

       deallocate(rimage, gimage, bimage)

     end subroutine writeFalseColourPPM

     real function imagePercentile(image, limit) result (out)
       type(IMAGETYPE) :: image
       real :: limit
       real, allocatable :: values(:)
       integer :: nValues
       integer :: i, j, k, n
       

       nValues = (image%nx*image%ny)

       n = 0
       allocate(values(nValues))
       do i = 1, image%nx
          do j = 1, image%ny
             if (image%pixel(i,j)%i > 1.e-29) then
                n = n + 1
                values(n) = image%pixel(i,j)%i
             endif
          enddo
       enddo
       call sort(n, values)
       k = nint(real(n)*limit/100.)
       out = values(k)
       deallocate(values)
     end function imagePercentile

!
!*******************************************************************************
!
!! WRITE_IMAGE creates a FITS primary array containing a 2-D image.
!

     subroutine writeFitsImage(image, filename, objectDistance, type)

! Arguments

       use input_variables, only: lamStart, lamEnd
       type(IMAGETYPE), intent(in)   :: image
       character (len=*), intent(in) :: filename, type
       real(double) :: objectDistance

#ifdef USECFITSIO
! Local variables
       integer :: status,unit,blocksize,bitpix,naxis,naxes(2)
       integer :: group,fpixel,nelements
       real, allocatable :: array(:,:)
       real(double) :: scale, dlam, lamCen

       logical :: simple,extend

       allocate(array(1:image%nx, 1:image%ny))
       call writeInfo("Writing fits image",TRIVIAL)
       status=0
       !
       !  Delete the file if it already exists, so we can then recreate it.
       !
       call deleteFitsFile ( filename, status )
       !
       !  Get an unused Logical Unit Number to use to open the FITS file.
       !
       call ftgiou ( unit, status )
       !
       !  Create the new empty FITS file.
       !
       blocksize=1
       call ftinit(unit,filename,blocksize,status)
       !
       !  Initialize parameters about the FITS image (300 x 200 16-bit integers).
       !
       simple=.true.
       bitpix=-32
       naxis=2
       naxes(1)=image%nx
       naxes(2)=image%ny
       extend=.true.
       !
       !  Write the required header keywords.
       !
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       !
       !  Write the array to the FITS file.
       !
       group=1
       fpixel=1
       nelements=naxes(1)*naxes(2)

       scale = 1.d20
       scale = scale / (objectDistance**2) 
       scale = scale * 1.d4
       scale = scale * 1.d-7
       dlam = (lamEnd - lamStart) * 2.d0 * angstoMicrons
       lamCen = (lamStart + lamEnd) / 2.d0 * angstoMicrons

       scale = scale / dlam * lamCen

       

       write(*,*) "objectdist, scale ",objectDistance,scale
       select case(type)
          case("intensity")
             array = image%pixel%i * scale
          case("stokesq")
             array = image%pixel%q * scale
          case("stokesu")
             array = image%pixel%u * scale
          case("pol")
             array = 100.*sqrt(image%pixel%q**2 + image%pixel%u**2)/image%pixel%i
          case DEFAULT
             write(*,*) "Unknown type in writefitsimage ",type
       end select

       call ftppre(unit,group,fpixel,nelements,array,status)
       !
       !  Write another optional keyword to the header.
       !
!       call ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)
       call ftpkye(unit,'CDELT1',4.848137e-10,6,' ',status)
       call ftpkye(unit,'CDELT2',4.848137e-10,6,' ',status)
       call ftpkye(unit,'CRVAL1',4.721411605,10,' ',status)
       call ftpkye(unit,'CRVAL2',-0.412388335,10,' ',status)
       call ftpkyj(unit,'CRPIX1',image%nx/2,' ',status)
       call ftpkyj(unit,'CRPIX2',image%ny/2,' ',status)
       !
       !  Close the file and free the unit number.
       !
       call ftclos(unit, status)
       call ftfiou(unit, status)
       !
       !  Check for any error, and if so print out error messages
       !
       if (status > 0) then
          call printFitserror(status)
       end if
#endif

     end subroutine writeFitsImage

     subroutine deleteFitsFile(filename,status)
       
! Arguments
       character ( len = * ) filename
       integer :: status

#ifdef USECFITSIO
! Local variables
       integer unit,blocksize

       !
       !  Simply return if status is greater than zero.
       !
       if (status > 0) then
          return
       end if
       !
       !  Get an unused Logical Unit Number to use to open the FITS file
       !
       call ftgiou ( unit, status )
       !
       !  Try to open the file, to see if it exists
       !
       call ftopen ( unit, filename, 1, blocksize, status )

       if ( status == 0 ) then
          !
          !  File was opened;  so now delete it 
          !
          call ftdelt(unit,status)

       else if (status == 103)then
          !
          !  File doesn't exist, so just reset status to zero and clear errors
          !
          status=0
          call ftcmsg

       else
          !
          !  There was some other error opening the file; delete the file anyway
          !
          status=0
          call ftcmsg
          call ftdelt(unit,status)
       end if
       !
       !  Free the unit number for later reuse.
       !
       call ftfiou(unit, status)
#endif

     end subroutine deleteFitsFile

     subroutine printFitsError(status)

       !
       !*******************************************************************************
       !
       !! PRINT_ERROR prints out the FITSIO error messages to the user.
       !
       ! Arguments
       integer :: status
#ifdef USECFITSIO

       ! Local variables
       character ( len = 30 ) errtext
       character ( len = 80 ) errmessage
       !
       !  Check if status is OK (no error); if so, simply return.
       !
       if (status <= 0) then
          return
       end if
       !
       !  Get the text string which describes the error
       !
       call ftgerr(status,errtext)
       print *,'FITSIO Error Status =',status,': ',errtext
       !
       !  Read and print out all the error messages on the FITSIO stack
       !
       call ftgmsg(errmessage)
       do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
       end do
#endif
     end subroutine printFitsError

     subroutine pixelLocate(image, xDist, yDist, ix, iy)
       type(IMAGETYPE) :: image
       real :: xDist, yDist
       integer :: ix, iy

       ix = 0
       iy = 0
       if ( (xDist >= (image%xAxisCentre(1)-image%dx/2.)).and. &
            (xDist <= (image%xAxisCentre(image%nx)+image%dx/2.)).and. &
            (yDist >= (image%yAxisCentre(1)-image%dy/2.)).and. &
            (yDist <= (image%yAxisCentre(image%ny)+image%dy/2.)) ) then
          call locate(image%xAxisLH, image%nx, xDist, ix)
          call locate(image%yAxisBottom, image%ny, yDist, iy)
       endif
     end subroutine pixelLocate
            

      
end module image_mod




