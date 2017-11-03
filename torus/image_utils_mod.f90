! This module holds image parameters in an array of the 
! imageParameters derived type. 
! Parameters should be set using a call to setImageParams and 
! accessed with the appropriate accessor function. 
!
! D. Acreman, September 2011

module image_utils_mod

  use vector_mod
  implicit none
  public :: setNumImages, setImageParams, getImageWavelength, getImageType, getImagenPixelsX, &
       getImagenPixelsY, getImageFilename, getAxisUnits, getImageInc, getImagePA, &
       getImageViewVec, getnImage, getImageOffsets, getImageSizeX, getImageSizeY, getFluxUnits

  private

!
! Derived type holding image parameters
!
  TYPE imageParameters
     real :: lambda
     character(len=80) :: type
     character(len=80) :: filename
     integer :: nPixelsX,   nPixelsY
     real    :: ImageSizeX, ImageSizeY
     real    :: inclination
     real    :: positionAngle
     real    :: offsetX, offsetY
     real :: time
     type(VECTOR) :: viewVec
  END TYPE imageParameters

!
! Array of image parameters to use for this run
!
  TYPE(imageParameters), private, allocatable, save :: myImages(:)
  integer, save :: numImages = 0 

! Parameters shared by all images
  character(len=10), save :: myAxisUnits(20)
  character(len=12), save :: myFluxUnits(20)

contains 

!
! Allocate the array containing the image parameters
!
  subroutine setNumImages(n)
    integer, intent(in) :: n 

    numImages = n 
    if (allocated(myImages)) deallocate(myImages)
    allocate(myImages(n))

  end subroutine setNumImages

!
! Populate ith element of the image parameter array.
! This routine also does some basic validation of the 
! values it has been given. 
!
  subroutine setImageParams(i, lambda, type, filename, nPixels, axisUnits, fluxUnits, &
       imageSize, aspectRatio, inclination, positionAngle, offsetX, offsetY, gridDistance, viewVec, imagetime)
    use messages_mod
    use constants_mod, only: autocm, pctocm, pi
    implicit none 

! Arguments
    integer, intent(in) :: i 
    real, intent(in)    :: lambda
    character(len=80), intent(in) :: type
    character(len=80), intent(in) :: filename
    integer, intent(in) :: nPixels
    character(len=10), intent(in) :: axisUnits
    character(len=12), intent(in) :: fluxUnits
    real, intent(in) :: imageSize
    real, intent(in) :: aspectRatio
    real, intent(in) :: inclination, positionAngle
    real, intent(in) :: offsetX, offsetY
    real, intent(in) :: gridDistance
    type(VECTOR), optional :: viewVec
    real, intent(in), optional :: imageTime

! Local variables
    character(len=80) :: message

    if ( i > numImages ) then 
       write(message,*) "Trying to set parameters for image ", i, &
            "but only ", numImages, " are allocated."
       call WriteFatal(message)
    end if

! Image wavelength
    if (lambda <= 0.0) then 
       write(message,*) "Image wavelength must be > 0. It is ", lambda
       call writeWarning(message)
    end if
    myImages(i)%lambda = lambda

! Type of image
    myImages(i)%type = type

! Filename for FITS file
    myImages(i)%filename = filename

! Aspect ratio of the image
    if (aspectRatio <= 0.0 ) then 
       write(message,*) "Aspect ratio must be positive ", aspectRatio
       call writeFatal(message)
    end if

! Number of pixels 
    if ( nPixels < 1) then 
       write(message,*) "Number of pixels should be > 0. It is ", nPixels
       call writeFatal(message)
    end if

    myImages(i)%nPixelsY = nPixels
    myImages(i)%nPixelsX = int(nPixels*aspectRatio)

! Axis units (e.g. au, pc, arcsec)
    myAxisUnits(i) = axisUnits

! flux units
    myFluxUnits(i) = fluxUnits

! Size of image
!
    if (imageSize <= 0.0 ) then 
       write(message,*) "Image size must be positive. It is ", imageSize
       call writeWarning(message)
    end if
!
! Set size of image y axis in cm
    select case (myAxisunits(i))
    case ("au","AU")
       myImages(i)%ImageSizeY = imageSize * real(autocm)
    case ("pc","PC")
       myImages(i)%ImageSizeY = imageSize * real(pctocm)
    case ("cm")
       myImages(i)%ImageSizeY = imageSize
    case ("arcsec")
       myImages(i)%ImageSizeY = & 
            real( (imageSize / 3600.0) * (pi/180.0) * gridDistance * pcToCm)
    case default
       myImages(i)%ImageSizeY = imageSize
       write(message,*) "Unrecognised units for image axis: ", trim(axisunits)
       call writeWarning(message)
    end select

! Set size of image x axis, taking the aspect ratio into account. If the number of 
! pixels is the same for both axes don't scale to avoid rounding errors. 
    if (myImages(i)%nPixelsX == myImages(i)%nPixelsY) then 
       myImages(i)%ImageSizeX = myImages(i)%ImageSizeY
    else
! Otherwise scale by the ratio of the number of pixels (not aspect ratio as 
! this could be different due to rounding to a whole number of pixels). 
        myImages(i)%ImageSizeX = myImages(i)%ImageSizeY &
             * real(myImages(i)%nPixelsX / myImages(i)%nPixelsY)
     endif

! Set image inclination and postion angle
    myImages(i)%inclination   = inclination
    myImages(i)%positionAngle = positionAngle
    myImages(i)%viewVec = setImageViewVec(i)
    if (PRESENT(viewVec)) myImages(i)%viewVec = viewVec
! Set the position of the image centre
    myImages(i)%offsetX = offsetX
    myImages(i)%offsetY = offsetY

    myImages(i)%time = 0.d0
    if (PRESENT(imageTime)) myImages(i)%time = imageTime

  end subroutine setImageParams

!
! Return the wavelength of the ith image
!
  real function getImageWavelength(i)
    
    implicit none
    integer, intent(in) :: i 

    getImageWavelength=myImages(i)%lambda

  end function getImageWavelength

!
! Return the type of the ith image
!
  function getImageType(i)

    character(len=80) :: getImageType
    integer, intent(in) :: i 

    if (i <= numImages) then 
       getImageType=myImages(i)%type
    else
       getImageType(1:4)="none"
    endif

  end function getImageType

!
! Return the number of x-axis pixels in the ith image
!
  integer function getImagenPixelsX(i)

     integer, intent(in) :: i 
     getImagenPixelsX=myImages(i)%nPixelsX

  end function getImagenPixelsX

!
! Return the number of y-axis pixels in the ith image
!
  integer function getImagenPixelsY(i)

     integer, intent(in) :: i 
     getImagenPixelsY=myImages(i)%nPixelsY

  end function getImagenPixelsY

!
! Return the filename of the ith image
!
  function getImageFilename(i)

    character(len=80) :: getImageFilename
    integer, intent(in) :: i 
    getImageFilename=myImages(i)%filename

  end function getImageFilename

!
! Return axis units
!
  function getAxisUnits(i)

    integer :: i
    character(len=10) :: getAxisUnits
    getAxisUnits=MyAxisUnits(i)

  end function getAxisUnits

  function getFluxUnits(i)

    integer :: i
    character(len=12) :: getFluxUnits
    getFluxUnits=MyFluxUnits(i)

  end function getFluxUnits

! Return size (length or angular) of ith image
! If an invalid image number is requested then -1 is returned. 
  function getImageSizeX(i)
    integer, intent(in) :: i 
    real :: getImageSizeX

    if (i <= numImages) then 
       getImageSizeX = myImages(i)%ImageSizeX
    else
       getImageSizeX = -1
    end if

  end function getImageSizeX

  function getImageSizeY(i)
    integer, intent(in) :: i 
    real :: getImageSizeY

    if (i <= numImages) then 
       getImageSizeY = myImages(i)%ImageSizeY
    else
       getImageSizeY = -1
    end if

  end function getImageSizeY

! Return inclination of ith image in radians
! If an invalid image number is requested then -1 is returned. 
  function getImageInc(i)
    integer, intent(in) :: i 
    real :: getImageInc

    if (i <= numImages) then 
       getImageInc = myImages(i)%inclination
    else
       getImageInc = -1
    end if

  end function getImageInc

! Return position angle of ith image in radians
! If an invalid image number is requested then -1 is returned. 
  function getImagePA(i)
    integer, intent(in) :: i 
    real :: getImagePA

    if (i <= numImages) then 
       getImagePA = myImages(i)%positionAngle
    else
       getImagePA = -1
    end if

  end function getImagePA

! Return the viewing vector for the ith image
  function setImageViewVec(i)
    use vector_mod, only: vector, rotateZ
    integer, intent(in) :: i
    TYPE(vector) :: setImageViewVec, tempVec

! Make a unit vector with the required inclination
! Inclination is spherical polar theta
    tempVec%x = -sin(myImages(i)%inclination)
    tempVec%y = 0.0
    tempVec%z = -cos(myImages(i)%inclination)
    
! Now apply the required position angle
! Position angle is spherical polar phi
    setImageViewVec = rotateZ(tempVec,dble(myImages(i)%positionAngle))

    tempVec%x = cos(myImages(i)%positionAngle) * sin(myImages(i)%inclination)
    tempVec%y = sin(myImages(i)%positionAngle) * sin(myImages(i)%inclination)
    tempVec%z = cos(myImages(i)%inclination)

    setImageViewVec = (-1.d0) * tempVec


  end function setImageViewVec

! Return the viewing vector for the ith image
  type(VECTOR) function getImageViewVec(i)
    integer, intent(in) :: i

    getImageViewVec = myImages(i)%viewVec

  end function getImageViewVec

! Return the number of images
  function getnImage()
    integer ::  getnImage
    getnImage = numImages
  end function getnImage

! Return x,y position of the image centre
  subroutine getImageOffsets(i, offsetX, offsetY)
    integer, intent(in) :: i
    real, intent(out) :: offsetX, offsetY

    offsetX = myImages(i)%offsetX
    offsetY = myImages(i)%offsetY

  end subroutine getImageOffsets

end module image_utils_mod
