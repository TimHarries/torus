module image_mod

  use constants_mod
  use vector_mod
  use messages_mod
  use phasematrix_mod
  use photon_mod
#ifdef FITSCUBE
  use datacube_mod
#endif

  implicit none

  public :: initImage, freeImage, addPhotonToImage, createLucyImage, writeFalseColourPPM,  ConvertArrayToMJanskiesPerStr
#ifdef PHOTOION
  public :: addPhotonToPhotoionImage
#endif
#ifdef USECFITSIO
  public writeFitsImage
#endif

  private :: pixelLocate, writePPMAimage, imagePercentile, &
      dumpLine, dumpPointTestData
  
  type IMAGETYPE
    type(STOKESVECTOR), pointer :: pixel(:,:) => null()
    real, pointer :: totWeight(:,:) => null()
    real, pointer :: xAxisCentre(:) => null()
    real, pointer :: yAxisCentre(:) => null()
    real, pointer :: xAxisLH(:) => null()
    real, pointer :: yAxisBottom(:) => null()
    integer, pointer :: nSamples(:,:) => null()
    integer :: nx
    integer :: ny
    real :: dx, dy
  end type IMAGETYPE

  contains

    subroutine smoothImage(image, fwhmPixels)
      type(IMAGETYPE) :: image
      real(double) :: fwhmPixels
      real(double), allocatable :: unsmoothedImage(:,:), smoothedImage(:,:)
      allocate(unsmoothedImage(1:image%nx, 1:image%ny))
      allocate(smoothedImage(1:image%nx, 1:image%ny))      

      unsmoothedImage(1:image%nx,1:image%ny) = image%pixel(1:image%nx,1:image%ny)%i
      call applySmooth(unSmoothedImage, smoothedImage, fwhmPixels)
      image%pixel(1:image%nx,1:image%ny)%i = smoothedImage(1:image%nx,1:image%ny)

      unsmoothedImage(1:image%nx,1:image%ny) = image%pixel(1:image%nx,1:image%ny)%q
      call applySmooth(unSmoothedImage, smoothedImage, fwhmPixels)
      image%pixel(1:image%nx,1:image%ny)%q = smoothedImage(1:image%nx,1:image%ny)

      unsmoothedImage(1:image%nx,1:image%ny) = image%pixel(1:image%nx,1:image%ny)%u
      call applySmooth(unSmoothedImage, smoothedImage, fwhmPixels)
      image%pixel(1:image%nx,1:image%ny)%u = smoothedImage(1:image%nx,1:image%ny)

      unsmoothedImage(1:image%nx,1:image%ny) = image%pixel(1:image%nx,1:image%ny)%v
      call applySmooth(unSmoothedImage, smoothedImage, fwhmPixels)
      image%pixel(1:image%nx,1:image%ny)%v = smoothedImage(1:image%nx,1:image%ny)
     
      deallocate(unsmoothedImage, smoothedImage)

    end subroutine smoothImage

    subroutine applySmooth(inputImage, outputImage, fwhmPixels)

      real(double) :: inputImage(:,:), outputImage(:,:), fwhmPixels, sigma, tot, f
      integer :: nx, ny, i , j, m, i1, j1

      sigma = fwhmPixels/2.355d0
      m  = nint(6.d0 * sigma)
      nx = size(inputimage,1)
      ny = size(inputimage,2)
      outputImage = 0.d0
      do i = 1, nx
         do j = 1, ny
            tot = 0.d0
            do  i1= i - m, i + m
               do j1 = j- m, j + m
                  f = exp(- ( (dble(i1-i)**2 / (2.d0*sigma**2)) + (dble(j1-j)**2 / (2.d0*sigma**2))) )

                  outputImage(i,j) = outputImage(i,j) + f * inputImage(min(max(i1,1),nx),min(max(j1,1),ny))
                  tot = tot + f
               enddo
            enddo
            outputImage(i,j) = outputImage(i,j) / tot
         enddo
      enddo
    
    end subroutine applySmooth
      


   function initImage(imNum)
     use image_utils_mod

     type(IMAGETYPE) :: initImage
     integer, intent(in) :: imNum

     integer :: i
     integer :: nx, ny
     real    :: imagesizeX, imageSizeY
     real    :: thisxOffset, thisyOffset

     ! Number of pixels
     nx = getImagenPixelsX(imNum)
     ny = getImagenPixelsY(imNum)

     ! Size of image
     imagesizeX = getImageSizeX(imNum)/1.0e10
     imageSizeY = getImageSizeY(imNum)/1.0e10

     ! Position of image centre
     call getImageOffsets(imNum, thisxOffset, thisyOffset)

     allocate(initImage%pixel(1:nx,1:ny))
     allocate(initImage%totWeight(1:nx,1:ny))
     allocate(initImage%nSamples(1:nx,1:ny))
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
     initImage%nSamples = 0
     initImage%totWeight = 0.

   end function initImage

#ifdef FITSCUBE
  subroutine addImageSliceToCube(cube, imageSlice, iLambda)
    type(DATACUBE) :: cube
    type(IMAGETYPE) :: imageSlice
    integer :: iLambda

    cube%intensity(1:cube%nx, 1:cube%ny, iLambda) = real(imageSlice%pixel(1:cube%nx, 1:cube%ny)%i)

  end subroutine addImageSliceToCube
#endif

#ifdef PHOTOION
   subroutine addPhotonToPhotoionImage(observerDirection, thisImage, thisPhoton, totalFlux)

     type(VECTOR), intent(in) :: observerDirection
     type(IMAGETYPE), intent(inout) :: thisImage
     type(PHOTON), intent(in) :: thisPhoton
     real(double), intent(inout) :: totalFlux
     type(VECTOR) :: xProj, yProj
     real :: xDist, yDist
     integer :: xPix, yPix

     type(VECTOR), parameter :: zAxis = VECTOR(0.d0, 0.d0, 1.d0)

     xPix = 0; yPix = 0

     xProj =  observerDirection .cross. zAxis
     call normalize(xProj)
     yProj = observerDirection .cross. xProj
     call normalize(yProj)
     xDist = real((thisPhoton%position) .dot. xProj)
     yDist = real((thisPhoton%position) .dot. yProj)
           

     call pixelLocate(thisImage, xDist, yDist, xPix, yPix)

     if ( (xPix >= 1)            .and.(yPix >= 1) .and. &
          (xPix <= thisImage%nx) .and.(yPix <= thisImage%ny)) then

        thisImage%pixel(xPix, yPix) = thisImage%pixel(xPix, yPix)  &
             + thisPhoton%stokes * oneOnFourPi * exp(-thisPhoton%tau) * thisPhoton%weight

        thisImage%nSamples(xPix, yPix) = thisImage%nSamples(xPix, yPix) + 1

     endif

     totalFlux = totalFlux + thisPhoton%stokes%i * oneOnFourPi * exp(-thisPhoton%tau) * thisPhoton%weight

   end subroutine addPhotonToPhotoionImage
#endif

   subroutine addPhotonToImage(viewVec, rotationAxis, thisImageSet, nImage, thisPhoton, &
                               thisVel, weight, filters, positionAngle, lambda0_cont)
     use inputs_mod, only : imageOrigin
     use filter_set_class
     use phasematrix_mod

     integer, intent(in) :: nImage  ! number of images in a set
     type(IMAGETYPE), intent(inout) :: thisImageSet(nImage)
     type(PHOTON) :: thisPhoton
     type(VECTOR) :: viewVec,  xProj, yProj, rotationAxis
     real :: xDist, yDist
     real(double) :: positionAngle!, ang
!     real :: r
     integer :: xPix, yPix
     real  :: thisVel
     real :: weight
     type(filter_set), intent(in) :: filters     
     real, intent(in), optional :: lambda0_cont  ! rest wavelength of contiuum photon
     !
     integer :: i
     real :: filter_response
     real(double) :: lambda_obs, r, ang

     xPix = 0; yPix = 0

     positionAngle=positionAngle

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


     do i = 1, nImage

        xProj = (-1.d0)*(rotationAxis.cross.viewVec)

        call normalize(xProj)
        yProj = (-1.d0)*(viewVec .cross. xProj)

        call normalize(yProj)
           
        xDist = real((thisPhoton%position-imageOrigin) .dot. xProj)
        yDist = real((thisPhoton%position-imageOrigin) .dot. yProj)

        r = sqrt(xDist**2 + yDist**2)
        ang = atan2(yDist, xDist)
        ang = ang + positionAngle
        xDist = real(r * cos(ang))
        yDist = real(r * sin(ang))
           
        call pixelLocate(thisImageSet(i), xDist, yDist, xPix, yPix)

        if ((xPix >= 1) .and. &
             (yPix >= 1) .and. &
             (xPix <= thisImageSet(i)%nx) .and. &
             (yPix <= thisImageSet(i)%ny)) then
           ! using a filter response function in filter_set_class here
           filter_response = real( pass_a_filter(filters, i, lambda_obs ))
!           write(*,*) xpix,ypix,thisPhoton%stokes%i*weight*filter_response, thisPhoton%stokes%i, weight, filter_response
                 
           thisImageSet(i)%pixel(xPix, yPix) = thisImageSet(i)%pixel(xPix, yPix)  &
                + thisPhoton%stokes * weight * filter_response
           thisImageSet(i)%totWeight(xPix,yPix) = real(thisImageSet(i)%totWeight(xPix,yPix) &
                + (thisPhoton%stokes%i*weight*filter_response))

        endif

     end do

     end subroutine addPhotonToImage

     subroutine freeImage(thisImage)
       type(IMAGETYPE) :: thisImage

       deallocate(thisImage%pixel)
       deallocate(thisImage%totWeight)
       deallocate(thisImage%nSamples)
       NULLIFY(thisImage%pixel)
       NULLIFY(thisImage%totWeight)
       NULLIFY(thisImage%nSamples)
     end subroutine freeImage

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
             t = real((min(f95, max(f5, image(1)%pixel(i,j)%i))-f5)/(f95-f5))
             rImage(i, j) = int(255.*t)
          enddo
       enddo

       f5 = imagePercentile(image(2), 5.)
       f95 = imagePercentile(image(2), 95.)
       write(*,*) f5, f95
       do i = 1, nx
          do j = 1, ny
             t = real((min(f95, max(f5, image(2)%pixel(i,j)%i))-f5)/(f95-f5))
             gImage(i, j) = int(255.*t)
          enddo
       enddo

       f5 = imagePercentile(image(3), 5.)
       f95 = imagePercentile(image(3), 95.)
       write(*,*) f5, f95
       do i = 1, nx
          do j = 1, ny
             t = real((min(f95, max(f5, image(3)%pixel(i,j)%i))-f5)/(f95-f5))
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
                values(n) = real(image%pixel(i,j)%i)
             endif
          enddo
       enddo
       call sort(n, values)
       k = nint(real(n)*limit/100.)
       out = values(k)
       deallocate(values)
     end function imagePercentile


! Convert from ergs/s/A to MJy/sr (10^6 ergs/s/cm^2/Hz/sr)
! Note there is no distance dependance as this is per cm^2 AND per sr.
     subroutine ConvertArrayToMJanskiesPerStr(array, lambda, dx, distance, samplings, pointTest, cylinderTest)
       use inputs_mod, only : dumpCut, cutType, sliceIndex
       real, intent(inout)      :: array(:,:)
       integer, intent(inout), optional      :: samplings(:,:)
!       real, intent(in)         :: lambda
       real         :: lambda
       real(double), intent(in) :: dx, distance
       real(double), parameter :: FluxToJanskies     = 1.e23_db ! ergs s^-1 cm^2 Hz^1
       real(double), parameter :: FluxToMegaJanskies = FluxToJanskies * 1.e-6_db
       real(double), parameter :: PerAngstromToPerCm = 1.e8_db
       real(double) :: nu, PerAngstromToPerHz, strad, scale
       logical, optional :: pointTest, cylinderTest


       strad = (dx*1.d10/distance)**2
       scale = 1.d20/distance**2

       nu = cspeed / ( real(lambda,db) * angstromtocm)
       PerAngstromToPerHz = PerAngstromToPerCm * (cSpeed / nu**2)

       ! Factor of 1.0e20 converts dx to cm from Torus units
       array = real(FluxToMegaJanskies * PerAngstromToPerHz  * array * scale / strad)

       if(present(pointTest)) then
          call dumpPointTestData(array, strad, perAngstromToPerHz)
       end if


       if(present(cylinderTesT) .or. dumpCut) then
          call dumpLine(array, strad, perAngstromToPerHz, dx, scale, samplings, sliceIndex, &
               cutType)         
       end if
     end subroutine ConvertArrayToMJanskiesPerStr

! Convert from ergs/s/A to MJy/sr (10^6 ergs/s/cm^2/Hz/sr)
! Note there is no distance dependance as this is per cm^2 AND per sr.
     subroutine ConvertArrayToJanskiesPerBeam(array, lambda, dx, distance, beamArea)
       real, intent(inout)      :: array(:,:)
!       real, intent(in)         :: lambda
       real         :: lambda
       real(double), intent(in) :: dx, distance, beamArea
       real(double), parameter :: FluxToJanskies     = 1.e23_db ! ergs s^-1 cm^2 Hz^1
       real(double), parameter :: PerAngstromToPerCm = 1.e8_db
       real(double) :: nu, PerAngstromToPerHz, strad, scale, beamAreaInPixels

       strad = (dx*1.d10/distance)**2
       scale = 1.d20/distance**2

       beamAreaInPixels = beamArea * distance**2 / (dx*1.d10)**2
       nu = cspeed / ( real(lambda,db) * angstromtocm)
       PerAngstromToPerHz = PerAngstromToPerCm * (cSpeed / nu**2)

       ! Factor of 1.0e20 converts dx to cm from Torus units
       array = real(FluxToJanskies * PerAngstromToPerHz  * array * scale * beamAreaInPixels)
     end subroutine ConvertArrayToJanskiesPerBeam

! Convert from ergs/s/A to magnitudes per square arcsec
     subroutine ConvertArrayToMagPerSqArcsec(array, lambda, dx, distance)
       use utils_mod, only : returnMagnitude
       real, intent(inout)      :: array(:,:)
       character(len=1) :: mag
       real(double), intent(in) :: distance, dx
       real(double), parameter :: FluxToJanskies     = 1.e23_db ! ergs s^-1 cm^2 Hz^1
       real(double), parameter :: PerAngstromToPerCm = 1.e8_db
       real(double) :: strad, scale
       real :: lambda
       integer :: i, j

       if (abs(lambda-1.22e4)/1.22e4 < 1.d-4) then
          mag = "J"
       else if (abs(lambda-1.63e4)/1.63e4 < 1.d-4) then
          mag = "H"
       else if (abs(lambda-2.19e4)/2.19e4 < 1.d-4) then
          mag = "K"
       else
          write(*,*) "unknown magnitude for lambda ",lambda
          stop
       endif
       strad = (dx*1.d10/distance)**2
       scale = 1.d20/distance**2
       array = array * real(scale / strad)
       array = array / 4.254517e10   ! from per str to per arcsec^2
       do i = 1, SIZE(array,1)
          do j = 1, SIZE(array,2)
             array(i,j) = real(returnMagnitude(dble(max(1.e-30,array(i,j))), mag))
          enddo
       enddo

     end subroutine ConvertArrayToMagPerSqArcsec

! Convert from ergs/s/A to Jy/Pix
     subroutine ConvertArrayToJanskysPerPix(array, lambda, distance)
       real, intent(inout)      :: array(:,:)
       real         :: lambda
       real(double), parameter :: FluxToJanskies     = 1.e23_db ! ergs s^-1 cm^2 Hz^1
       real(double), parameter :: FluxToMegaJanskies = FluxToJanskies * 1.e-6_db
       real(double), parameter :: PerAngstromToPerCm = 1.e8_db
       real(double) :: nu, PerAngstromToPerHz, scale, distance


       nu = cspeed / ( real(lambda,db) * angstromtocm)
       PerAngstromToPerHz = PerAngstromToPerCm * (cSpeed / nu**2)
       scale = 1.d20/distance**2

       array = real(FluxToJanskies * PerAngstromToPerHz  * array * scale)

     end subroutine ConvertArrayToJanskysPerPix


     subroutine dumpLine(array, strad, perAngstromToPerHz, dx, scale, samplings,sliceIndex, &
          cutType)

       real, intent(inout)      :: array(:,:)       
       integer, intent(inout)      :: samplings(:,:)       
       real(double) :: inImage, strad, perAngstromToPerHz, scale
       real(double), parameter :: FluxToJanskies     = 1.e23_db ! ergs s^-1 cm^$
       real(double), parameter :: FluxToMegaJanskies = FluxToJanskies * 1.e-6_db
       real(double) :: dx, r
       integer :: i, j
       character(len=*) :: cutType
       integer :: sliceIndex
       logical, save :: firstTime = .true.
       integer :: nx, ny

       open (123, file="pixelFile.dat", status="unknown")

       nx = size(array(:,1))
       ny = size(array(1,:))

       r = -((dx*nx)/2.d0) + (dx/2.d0)
       
       if (cutType == "vertical") then
          if(firstTime) then
             print *, "Dumping horizontal cut through pixel ", sliceIndex
             firstTime = .false.
          end if
          do i = 1, nx
             do j = 1, ny
                inImage = (array(i, j)*strad)/(FluxToMegaJanskies*PerAngstromToPerHz * scale)
                if(j == sliceIndex) then
                   write(123, *) r, inImage, samplings(i,j)
                   r = r + dx
                end if
             end do
          end do
       else if (cutType == "horizontal") then
          if(firstTime) then
             print *, "Dumping horizontal cut through pixel ", sliceIndex
             firstTime = .false.
          end if
          do i = 1, nx
             do j = 1, ny
                inImage = (array(i, j)*strad)/(FluxToMegaJanskies*PerAngstromToPerHz * scale)
                if(i == sliceIndex) then
                   write(123, *) r, inImage, samplings(i,j)
                   r = r + dx
                end if
             end do
          end do
       else
          print *, "Unrecognized cut type in dumpCut", cutType
          stop
       end if

       close(123)

     end subroutine dumpLine


     subroutine dumpPointTestData(array, strad, perAngstromToPerHz)

       real, intent(in)      :: array(:,:)
       real(double) :: inImage, strad, perAngstromToPerHz
       real(double), parameter :: FluxToJanskies     = 1.e23_db ! ergs s^-1 cm^2 Hz^1
       real(double), parameter :: FluxToMegaJanskies = FluxToJanskies * 1.e-6_db
       integer :: i, j
       logical :: found=.false.


          do i = 1, 201
             do j = 1, 201
                if((array(i,j)) /= 0.0) then
                   if(found) then
                   open (123, file="image_flux.dat", status="old")
                      print *, "More than one pixel with non-zero flux in point source test."
                      print *, "Halting test..."
                      write(123,*) "FAIL: More than one pixel has non zero flux"
                      stop
                      close(123)
                   else
                      open (123, file="image_flux.dat", status="old", position="append")
                      found = .true.
                      print *, "non zero flux at pixel (",i,",",j,")"
                      inImage = (array(i, j)*strad)/(FluxToMegaJanskies*PerAngstromToPerHz * 1.d20)
                      print *, "Flux at image: ", inImage
                      write(123, *) inImage
                      print *, " "
                      close(123)
                   end if
                end if
             end do
          end do
     end subroutine

!
!*******************************************************************************
!
!! WRITE_IMAGE creates a FITS primary array containing a 2-D image.
!

#ifdef USECFITSIO
     subroutine writeFitsImage(image, filename, objectDistance, type, lambdaImage, &
          pointTest, cylinderTest)

       use inputs_mod, only : fwhmPixels, beamArea
       use fits_utils_mod
       use utils_mod, only : returnFlux
       use image_utils_mod

! Arguments
       type(IMAGETYPE),intent(in) :: image
       character (len=*), intent(in) :: filename, type
       character(len=80) :: rFile
       real(double) :: objectDistance
       real, intent(in), optional :: lambdaImage
       logical, optional :: pointTest, cylinderTest
! Local variables
       integer :: status,unit,blocksize,bitpix,naxis,naxes(2)
       integer :: group,fpixel,nelements
       real, allocatable :: array(:,:)
       integer, allocatable :: samplings(:,:)
       integer :: i, j,k, n
       real(double) :: rMin, rMax, r, tot, starFlux

       real(double) :: scale,  dx, dy, phi
       logical :: simple,extend
       logical :: oldFilePresent
       character(len=80) :: message
       real(double) :: refValX, refValY ! co-ordinate values in reference pixel
       real(double) :: xt, yt

       dx = image%xAxisCentre(2) - image%xAxisCentre(1)
       dy = image%yAxisCentre(2) - image%yAxisCentre(1)


       allocate(array(1:image%nx, 1:image%ny))
       allocate(samplings(1:image%nx, 1:image%ny))
       call writeInfo("Writing fits image to: "//trim(filename),TRIVIAL)

       call checkBitpix(FitsBitpix)

       status=0
       !
       !  Delete the file if it already exists, so we can then recreate it.
       !
       inquire(file=filename, exist=oldFilePresent)
       if (oldFilePresent) then 
          call writeInfo("Removing old file", FORINFO)
          call deleteFile(filename, status)
       end if

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
       bitpix=fitsbitpix
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

       if (fwhmPixels > 0.d0) call smoothImage(image, fwhmpixels)

       scale = 1.
       array = 1.e-30
       select case(type)
          case("intensity")
             array = real(image%pixel%i * scale)
!             print *, "image%pixel%i", image%pixel%i
          case("stokesq")
             where (image%pixel%i /= 0.d0) 
                array = real(image%pixel%q/image%pixel%i)
             end where
          case("stokesu")
             where (image%pixel%i /= 0.d0) 
                array = real(image%pixel%u/image%pixel%i)
             end where
          case("pol")
             array = real(sqrt(image%pixel%q**2 + image%pixel%u**2))

          case("polr2")
             array = real(sqrt(image%pixel%q**2 + image%pixel%u**2))
             if (abs(lambdaImage - 1.22d4)/1.22d4 < 1.d-4) then
                starFlux = returnFlux(7.3d0, "J")
                starFlux = 10.d0**(-0.4d0 * 7.3d0) 
             else if (abs(lambdaImage - 1.63d4)/1.63d4 < 1.d-4) then
                starFlux = returnFlux(6.9d0, "H")
                starFlux = 10.d0**(-0.4d0 * 6.9d0)
             else
                starFlux = 1.d0
             endif
             do i = 1, image%nx
                do j = 1, image%ny
                   xt = (( image%xAxisCentre(i) * 1.d10)/objectDistance)*radiansToArcsec
                   yt = (( image%yAxisCentre(j) * 1.d10)/objectDistance)*radiansToArcsec
                   array(i,j) = array(i,j) * real(fourPi * (xt**2 + yt**2) /starFlux )
                enddo
             enddo

          case("qr")
             do i = 1, image%nx
                do j = 1, image%ny
                   phi = atan2(image%xAxisCentre(i), image%yAxisCentre(j))
                   array(i,j) = real(cos(2.d0*phi) * image%pixel(i,j)%q + sin(2.d0*phi)*image%pixel(i,j)%u)
                enddo
             enddo

          case("ur")
             do i = 1, image%nx
                do j = 1, image%ny
                   phi = atan2(image%xAxisCentre(i), image%yAxisCentre(j))
                   array(i,j) = real(-sin(2.d0*phi) * image%pixel(i,j)%q + cos(2.d0*phi)*image%pixel(i,j)%u)
                enddo
             enddo

          case("pa")
             array = real(-0.5*atan2(image%pixel%u,image%pixel%q)*radtodeg)
             where (array < 0.e0) 
                array = array + 180.e0
             end where
             where (array > 180.e0) 
                array = array - 180.e0
             end where
          case DEFAULT
             write(*,*) "Unknown type in writefitsimage ",type
       end select



       samplings = 0

! Convert pixel units if a wavelength has been provided
       if (present(lambdaImage)) then 
          write(message,"(a,f14.2,a)") "Converting flux units using lambda= ", &
               lambdaImage, " Angstrom"
          call writeInfo(message,FORINFO)
          if(present(pointTest)) then
             call ConvertArrayToMJanskiesPerStr(array, lambdaImage, dx, objectDistance, &
                  samplings, pointTest=.true.)
          else if(present(cylinderTest)) then
             samplings = image%nSamples
             call ConvertArrayToMJanskiesPerStr(array, lambdaImage, dx, objectDistance, &
                  samplings, cylinderTest=.true.)
          else
             select case (trim(getFluxUnits()))
             case("Jy/beam")
                call ConvertArrayToJanskiesPerBeam(array, lambdaImage, dx, objectDistance, beamArea)
             case("MJy/str")
                call ConvertArrayToMJanskiesPerStr(array, lambdaImage, dx, objectDistance, samplings)
             case("Jy/pix")
                call ConvertArrayToJanskysPerPix(array, lambdaImage, objectDistance)
             case("mag/arcsec2")
                call ConvertArrayToMagPerSqArcsec(array, lambdaImage, dx, objectDistance)
             case DEFAULT
                call writeFatal("Flux unit not recognised: "//trim(getFluxUnits()))
             end select
          end if
       else
          call writeInfo("No wavelength provided, not converting axis units")
       endif

       ! Add keywords for bitpix=16 and bitpix=8 
       call addScalingKeywords(maxval(array), minval(array), unit, bitpix)

       call ftppre(unit,group,fpixel,nelements,array,status)

       !
       !  Write another optional keyword to the header.
       !
!       call ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)

       select case (getFluxUnits())
          case("MJy/str")
             call ftpkys(unit,'BUNIT', "MJY/STR", "units of image values", status)
          case("Jy/pix")
             call ftpkys(unit,'BUNIT', "JY/PIXEL", "units of image values", status)
          case("Jy/beam")
             call ftpkys(unit,'BUNIT', "JY/BEAM", "units of image values", status)
          case("mag/arcsec2")
             call ftpkys(unit,'BUNIT', "MAG/ARCSEC2", "units of image values", status)
          case DEFAULT
             call writeFatal("Flux unit not recognised: "//trim(getFluxUnits()))
       end select




       ! write keywords and set values which depend on the axis units
       select case (getAxisUnits())
       case ("arcsec")
          dx = ((dx * 1.d10)/objectDistance)*radtodeg
          dy = ((dy * 1.d10)/objectDistance)*radtodeg
          refValX = (( image%xAxisCentre(1) * 1.d10)/objectDistance)*radiansToArcsec
          refValY = (( image%yAxisCentre(1) * 1.d10)/objectDistance)*radiansToArcsec
          refValX = 0.
          refValY = 0.
          call ftpkys(unit,'CUNIT1', "deg", "x axis unit", status)
          call ftpkys(unit,'CUNIT2', "deg", "y axis unit", status)
       case ("au", "AU")
          dx = (dx * 1.d10)/autocm
          dy = (dy * 1.d10)/autocm
          refValX = image%xAxisCentre(1) * 1.d10 / autocm
          refValY = image%yAxisCentre(1) * 1.d10 / autocm
          call ftpkys(unit,'CUNIT1', "AU", "x axis unit", status)
          call ftpkys(unit,'CUNIT2', "AU", "y axis unit", status)
       case ("pc","PC")
          dx = (dx * 1.d10)/pctocm
          dy = (dy * 1.d10)/pctocm
          refValX = image%xAxisCentre(1) * 1.d10 / pctocm
          refValY = image%yAxisCentre(1) * 1.d10 / pctocm
          call ftpkys(unit,'CUNIT1', "PC", "x axis unit", status)
          call ftpkys(unit,'CUNIT2', "PC", "y axis unit", status)
       case ("cm")
          dx = dx * 1.d10
          dy = dy * 1.d10
          refValX = image%xAxisCentre(1) * 1.d10
          refValY = image%yAxisCentre(1) * 1.d10
          call ftpkys(unit,'CUNIT1', "cm", "x axis unit", status)
          call ftpkys(unit,'CUNIT2', "cm", "y axis unit", status)
       case default
          call writeFatal("Unrecognised units for image axis")
       end select

       if (getAxisUnits() == "arcsec") then 
          ! write x-axis keywords
          call ftpkys(unit,'CTYPE1',"RA---SIN","x axis", status)
          call ftpkyd(unit,'CRPIX1',dble(image%nx/2.d0),-3,'reference pixel',status)
          call ftpkyd(unit,'CDELT1',dx,10,' ',status)
          call ftpkyd(unit,'CROTA1',0.d0,10,' ',status)
          call ftpkyd(unit,'CRVAL1',refValX,-5,'coordinate value at reference point',status)

          ! write y-axis keywords
          call ftpkys(unit,'CTYPE2',"DEC--SIN","y axis", status)
          call ftpkyd(unit,'CRPIX2',dble(image%ny/2.d0),-3,'reference pixel',status)
          call ftpkyd(unit,'CDELT2',dy,10 ,' ',status)
          call ftpkyd(unit,'CROTA2',0.d0,10,' ',status)
          call ftpkyd(unit,'CRVAL2',refValY,-5,'coordinate value at reference point',status)

          call ftpkyd(unit,'CD1_1',dx,10,' ',status)
          call ftpkyd(unit,'CD1_2',0.d0,10,' ',status)
          call ftpkyd(unit,'CD2_1',0.d0,10,' ',status)
          call ftpkyd(unit,'CD2_2',dy,10,' ',status)

       else
          ! write x-axis keywords
          call ftpkys(unit,'CTYPE1'," X","x axis", status)
          call ftpkyd(unit,'CRPIX1',0.5_db,-3,'reference pixel',status)
          call ftpkyd(unit,'CDELT1',dx,10,' ',status)
          call ftpkyd(unit,'CRVAL1',refValX,-3,'coordinate value at reference point',status)

          ! write y-axis keywords
          call ftpkys(unit,'CTYPE2'," Y","y axis", status)
          call ftpkyd(unit,'CRPIX2',0.5_db,-3,'reference pixel',status)
          call ftpkyd(unit,'CDELT2',dy,10 ,' ',status)
          call ftpkyd(unit,'CRVAL2',refValY,-3,'coordinate value at reference point',status)

       endif

       !
       !  Close the file and free the unit number.
       !
       call ftclos(unit, status)
       call ftfiou(unit, status)

       rfile = filename(1:(len(trim(filename))-5))//".dat"
       open(22,file=rfile,status="unknown",form="formatted")
       do k = 1,100
          rMin = dble(k-1)/100.d0 * image%xAxisCentre(image%nx)
          rMax = dble(k)/100.d0 *  image%xAxisCentre(image%nx)
          tot = 0.d0
          n = 0
          do  i = 1, image%nx
             do j = 1, image%ny
                r = sqrt(image%xAxisCentre(i)**2 + image%yAxisCentre(j)**2)
                if ((r >= rMin).and.(r < rMax)) then
                   n  = n + 1
                   tot = tot + array(i,j)
                endif
             enddo
          enddo
          write(22,*) 0.5d0*(rMin+rMax)*1.d10/objectDistance*radianstoarcsec,tot/dble(n)
       enddo
       close(22)
                
          
       !
       !  Check for any error, and if so print out error messages
       !
       if (status > 0) then
          call printFitserror(status)
       end if

     End subroutine writeFitsImage
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
            
  subroutine createLucyImage(grid, lambda, xArray, nLambda, source, nSource, imNum)
    use image_utils_mod
    use source_mod, only: SOURCETYPE, getElement, I_nu, distanceToSource
#ifdef USECFITSIO
    use inputs_mod, only : griddistance
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
    real(double) :: scale
    real(double) :: objectDistance
    integer, intent(in) :: imNum

    objectDistance = 2.25558e-8 * pctocm
    scale = 1.d20
    scale = scale / (objectDistance**2) 
    scale = scale * 1.d4
    scale = scale * 1.d-7

    ilambda = findIlambda(lambda, xArray, nLambda, ok)

    image   = initImage(imNum) 
    viewVec = getImageViewVec(imNum)

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
    call writeFitsimage(image, "test.fits", griddistance*pctocm, "intensity", lambda)
#endif
  end subroutine createLucyImage

#ifdef MPI
! Gather an MPI distributed image. 
  subroutine collateImages(thisImage, dest)
    use mpi
    implicit none

    type(IMAGETYPE) :: thisImage
    integer, optional, intent(in) :: dest
    integer :: thisDest ! destination of the reduce operation
    real, allocatable :: tempRealArray(:), tempRealArray2(:)
    real(double), allocatable :: tempDoubleArray(:), tempDoubleArray2(:)
    integer, allocatable :: tempIntArray(:), tempIntArray2(:)
    integer :: ierr

! By default reduce to rank zero process but allow other destinations
    if ( present(dest) ) then 
       thisDest = dest
    else
       thisDest = 0
    end if

     allocate(tempRealArray(SIZE(thisImage%pixel)))
     allocate(tempRealArray2(SIZE(thisImage%pixel)))
     allocate(tempDoubleArray(SIZE(thisImage%pixel)))
     allocate(tempDoubleArray2(SIZE(thisImage%pixel)))
     allocate(tempIntArray(SIZE(thisImage%pixel)))
     allocate(tempIntArray2(SIZE(thisImage%pixel)))
     tempRealArray = 0.0
     tempRealArray2 = 0.0
     tempDoubleArray = 0.0_db
     tempDoubleArray2 = 0.0_db
     tempIntArray = 0
     tempIntArray = 2

     call writeInfo ("Collating images...", FORINFO)
     tempDoubleArray = reshape(thisImage%pixel%i,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,thisDest,MPI_COMM_WORLD,ierr)
     thisImage%pixel%i = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%i))

     tempDoubleArray = reshape(thisImage%pixel%q,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,thisDest,MPI_COMM_WORLD,ierr)
     thisImage%pixel%q = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%q))

     tempDoubleArray = reshape(thisImage%pixel%u,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,thisDest,MPI_COMM_WORLD,ierr)
     thisImage%pixel%u = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%u))

     tempDoubleArray = reshape(thisImage%pixel%v,(/SIZE(tempDoubleArray)/))
     call MPI_REDUCE(tempDoubleArray,tempDoubleArray2,SIZE(tempDoubleArray),MPI_DOUBLE_PRECISION,&
                     MPI_SUM,thisDest,MPI_COMM_WORLD,ierr)
     thisImage%pixel%v = reshape(tempDoubleArray2,SHAPE(thisImage%pixel%v))

     tempIntArray = reshape(thisImage%nSamples,(/SIZE(tempIntArray)/))
     call MPI_REDUCE(tempIntArray,tempIntArray2,SIZE(tempIntArray),MPI_INTEGER,&
                     MPI_SUM,thisDest,MPI_COMM_WORLD,ierr)
     thisImage%nSamples = reshape(tempIntArray2,SHAPE(thisImage%nSamples))


     tempRealArray = reshape(thisImage%totWeight,(/SIZE(tempRealArray)/))
     call MPI_REDUCE(tempRealArray,tempRealArray2,SIZE(tempRealArray),MPI_REAL,&
                     MPI_SUM,thisDest,MPI_COMM_WORLD,ierr)
     thisImage%totWeight = reshape(tempRealArray2,SHAPE(thisImage%totWeight))

     call writeInfo ("Done.", FORINFO)
     deallocate(tempRealArray)
     deallocate(tempRealArray2)
     deallocate(tempDoubleArray)
     deallocate(tempDoubleArray2)
     deallocate(tempIntArray)
     deallocate(tempIntArray2)

  end subroutine collateImages




#endif

#ifdef USECFITSIO
  subroutine deletefile(filename,status)

    !  A simple little routine to delete a FITS file

    integer :: status,unit,blocksize
    character(len=*) filename

    !  Simply return if status is greater than zero
    if (status .gt. 0)return

    !  Get an unused Logical Unit Number to use to open the FITS file
    call ftgiou(unit,status)

    !  Try to open the file, to see if it exists
    call ftopen(unit,filename,1,blocksize,status)

    if (status .eq. 0)then
       !         file was opened;  so now delete it 
       call ftdelt(unit,status)
    else if (status .eq. 103)then
       !         file doesn't exist, so just reset status to zero and clear errors
       status=0
       call ftcmsg
    else
       !         there was some other error opening the file; delete the file anyway
       status=0
       call ftcmsg
       call ftdelt(unit,status)
    end if

!  Free the unit number for later reuse
      call ftfiou(unit, status)
    end subroutine deletefile
#endif
end module image_mod

