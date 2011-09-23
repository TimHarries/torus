! This module holds image parameters in an array of the 
! imageParameters derived type. 
! Parameters should be set using a call to setImageParams and 
! accessed with the appropriate accessor function. 
!
! D. Acreman, September 2011

module image_utils_mod

  implicit none

  public :: setNumImages, setImageParams, getImageWavelength, getImageType, getImagenPixels, &
       getImageFilename

  private

!
! Derived type holding image parameters
!
  TYPE imageParameters
     real :: lambda
     character(len=80) :: type
     character(len=80) :: filename
     integer :: nPixels
  END TYPE imageParameters

!
! Array of image parameters to use for this run
!
  TYPE(imageParameters), private, allocatable, save :: myImages(:)
  integer, save :: numImages = 0 

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
  subroutine setImageParams(i, lambda, type, filename, nPixels)
    use messages_mod
    implicit none 

! Arguments
    integer, intent(in) :: i 
    real, intent(in)    :: lambda
    character(len=80), intent(in) :: type
    character(len=80), intent(in) :: filename
    integer, intent(in) :: nPixels

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
    getImageType=myImages(i)%type

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

end module image_utils_mod
