module image_mod

  use constants_mod
  use vector_mod
  use messages_mod
  use phasematrix_mod
  use photon_mod

  implicit none

  public

  private :: pixelLocate
  
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

     
   function initImage(nx, ny, imageSizeX, imageSizeY, vMin, vMax, xOffset, yOffset)

     type(IMAGETYPE) :: initImage
     integer :: nx, ny
     real :: imagesizeX, imageSizeY
     real :: vmax, vmin
     real, optional :: xOffset, yOffset

     integer :: i
     real :: thisxOffset, thisyOffset

! Check the requested image size is sensible
! If not then set size to 1 to prevent problems with unallocated arrays
     if ( nx < 1 ) then 
        call writewarning("initImage: nx < 1")
        nx = 1
     end if

     if ( ny < 1 ) then 
        call writewarning("initImage: ny < 1")
        ny = 1
     end if

! Check that the image has a non-zero size otherwise axes will all be zero. 
     if (imageSizeX == 0.0 ) then 
        call writeWarning ("initImage: imageSizeX = 0")
     end if
     if (imageSizeY == 0.0 ) then 
        call writeWarning ("initImage: imageSizeY = 0")
     end if

! Set offsets if required
     if ( present(xOffset) ) then 
        thisxOffset=xOffset
     else
        thisxOffset = 0.0
     end if

     if ( present(yOffset) ) then 
        thisyOffset=yOffset
     else
        thisyOffset = 0.0
     end if

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
        initImage%xAxisCentre(i) = -imageSizeX/2. + initimage%dx/2. + initimage%dx*real(i-1) + thisxOffset
        initImage%xAxisLH(i) = -imageSizeX/2. + initimage%dx*real(i-1) + thisxOffset
     enddo
     do i = 1, ny
        initImage%yAxisCentre(i) = -imageSizeY/2. + initimage%dy/2. + initimage%dy*real(i-1) + thisyOffset
        initImage%yAxisBottom(i) = -imageSizeY/2. + initimage%dy*real(i-1) + thisyOffset
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
     
#ifdef PHOTOION
   subroutine addPhotonToPhotoionImage(observerDirection, thisImage, thisPhoton, totalFlux)

     type(VECTOR), intent(in) :: observerDirection
     type(IMAGETYPE), intent(inout) :: thisImage
!     type(PHOTON), intent(in) :: thisPhoton
     type(PHOTON) :: thisPhoton
     real(double), intent(inout) :: totalFlux
     type(VECTOR) :: xProj, yProj
     real :: xDist, yDist
     integer :: xPix, yPix

     type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.d0)

     xPix = 0; yPix = 0

     xProj =  zAxis .cross. observerDirection
     call normalize(xProj)
     yProj = observerDirection .cross. xProj
     call normalize(yProj)
     xDist = (thisPhoton%position) .dot. xProj
     yDist = (thisPhoton%position) .dot. yProj
           

     call pixelLocate(thisImage, xDist, yDist, xPix, yPix)

     if ( (xPix >= 1)            .and.(yPix >= 1) .and. &
          (xPix <= thisImage%nx) .and.(yPix <= thisImage%ny)) then

        if(thisPhoton%weight == 0.d0) thisPhoton%weight = 1.d0

        thisImage%pixel(xPix, yPix) = thisImage%pixel(xPix, yPix)  &
             + thisPhoton%stokes * oneOnFourPi * exp(-thisPhoton%tau) * thisPhoton%weight
     endif

     totalFlux = totalFlux + thisPhoton%stokes%i * oneOnFourPi * exp(-thisPhoton%tau) * thisPhoton%weight

   end subroutine addPhotonToPhotoionImage
#endif

   subroutine addPhotonToImage(viewVec, rotationAxis, thisImageSet, nImage, thisPhoton, &
                               thisVel, weight, filters, center, lambda0_cont)
     use filter_set_class
     use phasematrix_mod

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
     type(vector) :: test
     test = center

     velIncgs = thisVel! * cSpeed/1.e5 
     xPix = 0; yPix = 0


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

     subroutine writeFalseColourPPM(tfile, image)

       character(len=*) :: tfile
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
       use utils_mod, only: sort
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

! Convert from ergs/s/A to MJy/sr (10^6 ergs/s/cm^2/Hz/sr)
! Note there is no distance dependance as this is per cm^2 AND per sr.
     subroutine ConvertArrayToMJanskiesPerStr(array, lambda, dx, distance)
       real, intent(inout)      :: array(:,:)
!       real, intent(in)         :: lambda
       real         :: lambda
       real(double), intent(in) :: dx, distance
       real(double), parameter :: FluxToJanskies     = 1.e23_db ! ergs s^-1 cm^2 Hz^1
       real(double), parameter :: FluxToMegaJanskies = FluxToJanskies * 1.e-6_db
       real(double), parameter :: PerAngstromToPerCm = 1.e8_db
       real(double) :: nu, PerAngstromToPerHz, strad, scale
       logical :: pointTest = .true.
       integer i, j
       
       strad = (dx*1.d10/distance)**2
       scale = 1.d20/distance**2
       write(*,*) "dx ", dx
       write(*,*) "distance ",distance
       write(*,*) "ang (arcsec) ", sqrt(strad)*radtodeg*3600.d0

!       if(lambda == 0.d0) lambda = 6.e8

       nu = cspeed / ( real(lambda,db) * angstromtocm)
       PerAngstromToPerHz = PerAngstromToPerCm * (cSpeed / nu**2)

       if(pointTest) then
          do i = 1, 128
             do j = 1, 128
                if((array(i,j)*scale) /= 0.0) then
                   print *, "non zero flux at pixel (",i,",",j,")"
                   print *, "F = ", array(i,j)*scale
                   print *, "F/inPixel ", array(i,j)*scale/(128.**2)
!                   print *, "F/cm ", array(i,j)*scale/(128.**2)
                   print *, "F/strad = ", array(i,j)*scale/strad
                   write(*,*) "flux in mjy ",array(i,j)*fluxtomegajanskies * perAngstromtoperhz * &
                        scale,lambda
                   print *, " "
                end if
             end do
          end do
       end if
       
       ! Factor of 1.0e20 converts dx to cm from Torus units
       array = FluxToMegaJanskies * PerAngstromToPerHz  * array * scale / strad

       print *, "FluxToMegaJanskies ", FluxToMegaJanskies
       print *, "PerAngstromToPerHz ", PerAngstromToPerHz
       print *, "scale ", scale
       print *, "strad ", strad
       print *, "PerAngstromToPerCm ", PerAngstromToPerCm
       print *, "cSpeed ", cSpeed
       print *, "nu ", nu
       print *, "angstromtocm ", angstromtocm
!       print *, "array ", array

     end subroutine ConvertArrayToMJanskiesPerStr
!
!*******************************************************************************
!
!! WRITE_IMAGE creates a FITS primary array containing a 2-D image.
!

#ifdef USECFITSIO
     subroutine writeFitsImage(image, filename, objectDistance, type)

! Arguments
       
       use input_variables, only: lamStart, ImageinArcSec
       type(IMAGETYPE),intent(in) :: image
       character (len=*), intent(in) :: filename, type
       real(double) :: objectDistance
! Local variables
       integer :: status,unit,blocksize,bitpix,naxis,naxes(2)
       integer :: group,fpixel,nelements
       real, allocatable :: array(:,:)
       real(double) :: scale,  dx, dy

       logical :: simple,extend

       dx = image%xAxisCentre(2) - image%xAxisCentre(1)
       dy = image%yAxisCentre(2) - image%yAxisCentre(1)


       allocate(array(1:image%nx, 1:image%ny))
       call writeInfo("Writing fits image to: "//trim(filename),TRIVIAL)
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

       scale = 1.
       array = 1.e-30
       select case(type)
          case("intensity")
             array = image%pixel%i * scale
!             print *, "image%pixel%i", image%pixel%i
          case("stokesq")
             where (image%pixel%i /= 0.d0) 
                array = image%pixel%q * scale
             end where
          case("stokesu")
             where (image%pixel%i /= 0.d0) 
                array = image%pixel%u * scale 
             end where
          case("pol")
             where (image%pixel%i /= 0.d0) 
                array = sqrt(image%pixel%q**2 + image%pixel%u**2)
             end where
          
          case DEFAULT
             write(*,*) "Unknown type in writefitsimage ",type
       end select

!       if(lamStart = 0.d0) 
      
       call ConvertArrayToMJanskiesPerStr(array, lamstart, dx, objectDistance)
       call ftppre(unit,group,fpixel,nelements,array,status)
       !
       !  Write another optional keyword to the header.
       !
!       call ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)



       if (imageinarcsec) then
          dx = ((dx * 1.d10)/objectDistance)*radiansToArcsec
          dy = ((dy * 1.d10)/objectDistance)*radiansToArcsec
          call ftpkys(unit,'CTYPE1'," X","x axis", status)
          call ftpkys(unit,'CUNIT1', "arcsec", "x axis unit", status)
          call ftpkys(unit,'CTYPE2'," Y","y axis", status)
          call ftpkyd(unit,'PIXEL', dx*1000.d0, 10, " ", status)
       else
          dx = (dx * 1.d10)/autocm
          dy = (dy * 1.d10)/autocm
          call ftpkys(unit,'CTYPE1'," X","x axis", status)
          call ftpkys(unit,'CUNIT1', "AU", "x axis unit", status)
          call ftpkys(unit,'CTYPE2'," Y","y axis", status)
          call ftpkys(unit,'CUNIT2', "AU", "y axis unit", status)
       endif



       call ftpkyd(unit,'CDELT1',dx,10,' ',status)
       call ftpkyd(unit,'CDELT2',dy,10 ,' ',status)
       call ftpkyd(unit,'CRVAL1',0.d0,10,' ',status)
       call ftpkyd(unit,'CRVAL2',0.d0,10,' ',status)
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

     End subroutine writeFitsImage
#endif

#ifdef USECFITSIO
     subroutine deleteFitsFile(filename,status)
       
! Arguments
       character ( len = * ) filename
       integer :: status

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

     end subroutine deleteFitsFile
#endif

#ifdef USECFITSIO
     subroutine printFitsError(status)

       !
       !*******************************************************************************
       !
       !! PRINT_ERROR prints out the FITSIO error messages to the user.
       !
       ! Arguments
       integer :: status

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
     end subroutine printFitsError
#endif

     subroutine pixelLocate(image, xDist, yDist, ix, iy)
       use utils_mod, only: locate
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
          if (xDist >= image%xAxisCentre(image%nx)-image%dx/2.) ix = image%nx
          if (yDist >= image%yAxisCentre(image%ny)-image%dy/2.) iy = image%ny
       endif
     end subroutine pixelLocate
            
  subroutine createLucyImage(grid, viewVec, lambda, xArray, nLambda, source, nSource)
    use input_variables, only : npix, setimageSize, vmin, vmax
    use source_mod, only: SOURCETYPE, getElement, I_nu, distanceToSource
#ifdef USECFITSIO
    use input_variables, only : griddistance
#endif
    use amr_mod, only: distanceToGridFromOutside, distanceToCellBoundary, findSubcellLocal, returnKappa, inOctal, &
         returnScatteredIntensity
    use atom_mod, only: bnu
    use utils_mod, only: findILambda
    use gridtype_mod, only: GRIDTYPE
    use octal_mod, only: OCTAL
    type(OCTAL), pointer :: thisOctal
    type(SOURCETYPE) :: source(:)
    integer :: nSource
    integer :: subcell
    real :: lambda, xArray(:)
    integer :: ilambda, nLambda, sourceNumber
    real(double) :: tVal
    type(GRIDTYPE) :: grid
    type(IMAGETYPE) :: image
    type(VECTOR) :: viewVec, xVec, yVec, position, photoDirection, thisVec
    integer :: iElement
    type(VECTOR) :: pVec
    real(double) :: cosTheta
    real(double) :: distToGrid, distToSource, currentDistance
    integer :: ix, iy
    real(double) :: i0, tau, dtau, kappaAbs, jnu, kappaSca, iScattered
    logical :: ok, hitsource
    integer, parameter :: nTheta = 11, nPhi = 10
    integer :: iTheta, iPhi
    real(double) :: thisTheta, thisPhi
    real(double) :: scale, dlam, lamCen, lamStart, lamEnd
    real(double) :: objectDistance


    lamStart = 1.d4
    lamEnd = 1.01e4
    objectDistance = 2.25558e-8 * pctocm

    scale = 1.d20
    scale = scale / (objectDistance**2) 
    scale = scale * 1.d4
    scale = scale * 1.d-7
    dlam = (lamEnd - lamStart) * 2.d0 * angstoMicrons
    lamCen = (lamStart + lamEnd) / 2.d0 * angstoMicrons


    ilambda = findIlambda(lambda, xArray, nLambda, ok)

    image = initImage(npix, npix, setimageSize/1.e10, setimageSize/1.e10, vmin, vmax) 

    xVec = zHat .cross. viewVec
    call normalize(xVec)
    yVec = xVec .cross. viewVec
    call normalize(yVec)

    thisVec = (-1.d0)*viewVec
    thisTheta = acos(thisvec%z)
    thisPhi = atan2(thisVec%y,thisVec%x)

    if (thisPhi < 0.d0) thisPhi = thisPhi + twoPi 
    iTheta = nint((thisTheta / pi) * dble(nTheta-1))+1
    iphi = nint((thisPhi / twoPi) * dble(nPhi-1))+1

    write(*,*) "Scattered Light indices " , itheta, iphi

    do ix = 1, image%nx
       do iy = 1, image%ny
          position = (-10.d0*grid%octreeRoot%subcellSize)*viewVec + &
               (dble(image%xAxisCentre(ix)) * xVec) + (dble(image%yAxisCentre(iy)) * yVec)
          distToGrid = distanceToGridFromOutside(grid, position, viewVec) 
          position = position + (1.d-1*grid%halfSmallestSubcell+distToGrid) * viewVec


          call distanceToSource(source, nSource, position, viewVec, hitSource, disttosource, sourcenumber)
          if (.not.hitSource) distToSource = 1.d30

          if (hitSource) then
             pVec = (position + (viewVec * distToSource) - source(sourceNumber)%position)
             call normalize(pVec)
             cosTheta = -1.d0*(pVec.dot.viewVec)
             photoDirection = pVec
             call normalize(PhotoDirection)
             iElement = getElement(source(sourcenumber)%surface, photoDirection)
          endif

          thisOctal => grid%octreeRoot
          subcell = 1
          i0 = 0.d0
          tau = 0.d0
          currentDistance  = 0.d0
          do while(inOctal(grid%octreeRoot, position).and.(currentDistance < disttoSource))
             call findSubcelllocal(position, thisOctal, subcell)
             call distanceToCellBoundary(grid, position, viewVec, tVal, sOctal=thisOctal)
             call returnKappa(grid, thisOctal, subcell, ilambda = ilambda, kappaAbs = kappaAbs, kappaSca = kappaSca)
             dTau = (kappaAbs + kappaSca) * tVal

             jnu = kappaAbs * bnu(cspeed/(lambda*angstromTocm), dble(thisOctal%temperature(subcell)))


             iScattered = returnScatteredIntensity(position,thisOctal, subcell, (-1.d0)*viewVec)
!             write(*,*) iScattered
!             write(*,*) "subcell ",thisOctal%scatteredIntensity(subcell,:,:)
             i0 = i0 + exp(-tau) * (1.d0 - exp(-dtau)) * iScattered

!             if (kappaAbs .ne. 0.d0) then
!                
!                snu = jnu/kappaAbs
!                i0 = i0 +  exp(-tau) * (1.d0-exp(-dtau))*snu
!             else
!                snu = tiny(snu)
!                i0 = i0 + tiny(i0)
!             endif
             tau = tau + dtau
             position = position + (tval+1.d-3*grid%halfSmallestSubcell) * viewVec
          end do
          if (hitSource) then
             i0 = i0 + i_nu(source(sourceNumber), cSpeed/(lambda*angstromtocm), iElement, cosTheta)*exp(-tau)
             write(*,*) "hit core ", &
            i_nu(source(sourceNumber), cSpeed/(lambda*angstromtocm), iElement, cosTheta)*exp(-tau), tau, lambda
          endif
          image%pixel(ix, iy)%i = i0 * scale
       end do
    end do

#ifdef USECFITSIO
    call writeFitsimage(image, "test.fits", griddistance*pctocm, "intensity")
#endif
  end subroutine createLucyImage

#ifdef MPI
! Gather an MPI distributed image. 
  subroutine collateImages(thisImage)
    implicit none
    include 'mpif.h'
    type(IMAGETYPE) :: thisImage
    real, allocatable :: tempRealArray(:), tempRealArray2(:)
    real(double), allocatable :: tempDoubleArray(:), tempDoubleArray2(:)
    integer :: ierr

     allocate(tempRealArray(SIZE(thisImage%pixel)))
     allocate(tempRealArray2(SIZE(thisImage%pixel)))
     allocate(tempDoubleArray(SIZE(thisImage%pixel)))
     allocate(tempDoubleArray2(SIZE(thisImage%pixel)))
     tempRealArray = 0.0
     tempRealArray2 = 0.0
     tempDoubleArray = 0.0_db
     tempDoubleArray2 = 0.0_db

     call writeInfo ("Collating images...", FORINFO)
     tempDoubleArray = reshape(thisImage%pixel%i,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%i = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%i))

     tempDoubleArray = reshape(thisImage%pixel%q,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%q = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%q))

     tempDoubleArray = reshape(thisImage%pixel%u,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%u = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%u))

     tempDoubleArray = reshape(thisImage%pixel%v,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%pixel%v = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%v))


     tempRealArray = reshape(thisImage%vel,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%vel = reshape(tempRealArray2,SHAPE(thisImage%vel))

     tempRealArray = reshape(thisImage%totWeight,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,0,MPI_COMM_WORLD,ierr)
     thisImage%totWeight = reshape(tempRealArray2,SHAPE(thisImage%totWeight))

     call writeInfo ("Done.", FORINFO)
     deallocate(tempRealArray)
     deallocate(tempRealArray2)
     deallocate(tempDoubleArray)
     deallocate(tempDoubleArray2)

  end subroutine collateImages

#endif

end module image_mod


