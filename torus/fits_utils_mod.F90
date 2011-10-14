! Module containing general FITS file utility subroutines. 
! Anything specific to fits images should go in image_mod and 
! anything specific to data cubes should go in datacube_mod. 
! 
! D. Acreman, July 2011

module fits_utils_mod

  implicit none

  contains 

    ! Report whether the value of bitpix is valid
    subroutine checkBitpix (thisBitpix)
      use messages_mod

      integer, intent(in)  :: thisBitPix
      character(len=40)    :: message

      if ( thisBitPix ==   8   .or. thisBitPix ==  16 .or. thisBitPix ==  32 &
           .or. thisBitPix == -32   .or. thisBitPix == -64 ) then 
    
         write(message,'(a,i3)') "Writing FITS files using BITPIX= ", thisBitPix
         call writeInfo(message)

      else

         write(message,'(a,i3)') "Invalid value of  BITPIX= ", thisBitPix
         call writeWarning(message)

      end if

    end subroutine checkBitpix

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

!
! Set scaling values used when packing data into an integer representation
! These need to be set correctly in order for the packing to work
!
    subroutine addScalingKeywords(max, min, unit, bitpix)
      use kind_mod
      real, intent(in) :: max, min        ! data max and min values
      integer, intent(in) :: unit         ! fits file unit number
      integer, intent(in) :: bitpix
      real(double) :: bscale, bzero       ! values of keywords
      integer :: blank
      integer :: status

      if ( bitpix == 8 ) then 
         ! unsigned integer
         bscale = (max - min) / 255.0
         bzero  = min
         blank  = 0
      elseif (  bitpix == 16) then
         ! twos complement integer
         bscale = (max - min) / 65534.0
         bzero = (max + min) / 2.0
         blank = -32768
      elseif ( bitpix == 32 ) then
         ! twos complement integer
         ! subtract 1000 to provide some leeway so the conversion doesn't overflow
         bscale = (max - min) / real ( (2_bigInt**bitpix)-1000_bigInt )
         bzero = (max + min) / 2.0
         blank = int((2_bigInt**bitpix)/2_bigInt * (-1))
      else
         ! float so no scaling required
         return
      endif

      status = 0
      call ftpkyd(unit,'BSCALE',   bscale,-9,'Pixel scaling',status)
      if (status > 0) call printFitsError(status)
      call ftpkyd(unit,'BZERO',    bzero, -9,'Pixel offset',status)
      if (status > 0) call printFitsError(status)
      call ftpkyj(unit,'BLANK', blank, 6, ' ', status)
      if (status > 0) call printFitsError(status)

    end subroutine addScalingKeywords

#endif

end module fits_utils_mod
