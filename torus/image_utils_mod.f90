! This module holds image parameters in an array of the 
! imageParameters derived type. 
! Parameters should be set using a call to setImageParams and 
! accessed with the appropriate accessor function. 
!
! D. Acreman, September 2011

module image_utils_mod

  implicit none

  public :: setNumImages, setImageParams, getImageWavelength, getImageType, getImagenPixels, &
       getImageFilename, getAxisUnits, getImageSize, getImageInc, getImagePA, getImageViewVec

  private

!
! Derived type holding image parameters
!
  TYPE imageParameters
     real :: lambda
     character(len=80) :: type
     character(len=80) :: filename
     integer :: nPixels
     real    :: ImageSize
     real    :: inclination
     real    :: positionAngle
  END TYPE imageParameters

!
! Array of image parameters to use for this run
!
  TYPE(imageParameters), private, allocatable, save :: myImages(:)
  integer, save :: numImages = 0 

! Parameters shared by all images
  character(len=10), save :: myAxisUnits

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
  subroutine setImageParams(i, lambda, type, filename, nPixels, axisUnits, &
       imageSize, inclination, positionAngle)
    use messages_mod
    use constants_mod, only: autocm, pctocm, radToDeg
    implicit none 

! Arguments
    integer, intent(in) :: i 
    real, intent(in)    :: lambda
    character(len=80), intent(in) :: type
    character(len=80), intent(in) :: filename
    integer, intent(in) :: nPixels
    character(len=10), intent(in) :: axisUnits
    real, intent(in) :: imageSize
    real, intent(in) :: inclination, positionAngle

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

! Number of pixels 
    if ( nPixels < 1) then 
       write(message,*) "Number of pixels should be > 0. It is ", nPixels
       call writeWarning(message)
    end if
    myImages(i)%nPixels = nPixels

! Axis units (e.g. au, pc, arcsec)
    myAxisUnits = axisUnits

! Size of image
!
    if (imageSize <= 0.0 ) then 
       write(message,*) "Image size must be positive. It is ", imageSize
       call writeWarning(message)
    end if
!
! Set image size in cm or arcsec
    select case (myAxisunits)
    case ("au","AU")
       myImages(i)%ImageSize = imageSize * real(autocm)
    case ("pc","PC")
       myImages(i)%ImageSize = imageSize * real(pctocm)
    case ("cm")
       myImages(i)%ImageSize = imageSize
    case ("arcsec")
       myImages(i)%ImageSize = imageSize
    case default
       myImages(i)%ImageSize = imageSize
       write(message,*) "Unrecognised units for image axis: ", trim(axisunits)
       call writeWarning(message)
    end select

! Set image inclination and postion angle
    myImages(i)%inclination   = inclination
    myImages(i)%positionAngle = positionAngle

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
! Return the number of pixels in the ith image
!
  integer function getImagenPixels(i)

     integer, intent(in) :: i 
     getImagenPixels=myImages(i)%nPixels

  end function getImagenPixels

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
  function getAxisUnits()

    character(len=10) :: getAxisUnits
    getAxisUnits=MyAxisUnits

  end function getAxisUnits

! Return size (length or angular) of ith image
! If an invalid image number is requested then -1 is returned. 
  function getImageSize(i)
    integer, intent(in) :: i 
    real :: getImageSize

    if (i <= numImages) then 
       getImageSize = myImages(i)%ImageSize
    else
       getImageSize = -1
    end if

  end function getImageSize

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
  function getImageViewVec(i)
    use vector_mod, only: vector, rotateZ
    integer, intent(in) :: i
    TYPE(vector) :: getImageViewVec, tempVec

! Make a unit vector with the required inclination
! Inclination is spherical polar theta
    tempVec%x = -sin(myImages(i)%inclination)
    tempVec%y = 0.0
    tempVec%z = -cos(myImages(i)%inclination)
    
! Now apply the required position angle
! Position angle is spherical polar phi
    getImageViewVec = rotateZ(tempVec,dble(myImages(i)%positionAngle))

  end function getImageViewVec

end module image_utils_mod
