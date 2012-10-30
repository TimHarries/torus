#ifdef USECFITSIO

! Module for setting up a Torus grid from data held in a FITS file

! To do: 
! 1: Set filename as an input parameter
! 2. Interpolation for when cell centres don't match axes in FITS file

module gridFromFitsFile
  use kind_mod 
  use messages_mod
  implicit none

  public :: read_fits_file_for_grid, assign_from_fitsfile
  private :: setup_axes, read_lutable

! Filename
  character(len=*), parameter, private :: filename = "vh1.fits"

! Arrays to hold data from fits files
! Number of pixels
  integer, private, save :: axis_size(3)
! Axes
  real(double), private, save, allocatable :: xaxis(:), yaxis(:), zaxis(:)
  real(double), private, save :: xMin, xMax, dx
  real(double), private, save :: yMin, yMax, dy
  real(double), private, save :: zMin, zMax, dz
! Values
  real(single), private, save, allocatable :: density(:,:,:), temperature(:,:,:)

! If true use the ideal gas EOS to calculate temperature
  logical, parameter, private :: idealGas=.false.

! File for look-up table for temperature conversion. 
  character(len=*), parameter, private :: lu_table="tempdependent_mu_and_gamma.dat"
  integer, parameter, private :: npd = 91
  real, private, save :: pdlogt(npd), pd(npd)

  contains

!-------------------------------------------------------------------------------
! Reads in a FITS file which will be used to set up the Torus AMR grid
! Based on the FITS cookbook
! D. Acreman, October 2012

    subroutine read_fits_file_for_grid

      use fits_utils_mod, only: printFitsError
      use constants_mod, only: Rgas

      character(len=80) :: message

! FITS file variables and parameters
      integer :: status, unit, blocksize, bitpix
      integer, parameter :: readwrite=0 ! Open file read only 
      integer, parameter :: group=0
      integer :: naxis, nfound, npixels, hdutype
      integer :: i, j, k, n
      real :: grad1
      character(len=80) :: record, comment
      logical :: anynull
      real(single) :: maxdata, mindata

      if (.not.idealGas) call read_lutable

      call writeInfo ("Initialising grid from a FITS file", IMPORTANT)
      call writeInfo ("Reading FITS file "//filename, FORINFO)
      status = 0

! Get a free LUN then open the file
      call ftgiou(unit, status)
      call ftopen(unit, filename,readwrite,blocksize,status)

! Find out how many axes there are
      call ftgkyj(unit,"NAXIS",naxis,comment,status)
      write(message,*) "There are ", naxis, "axes"
      call writeInfo(message,TRIVIAL)

! Find out the size of the axis
      axis_size(:)=1
      call ftgknj(unit,"NAXIS",1,naxis,axis_size(1:naxis),nfound,status)
      write(message,*) "Axis sizes: ", axis_size
      call writeInfo(message,TRIVIAL)
      npixels = axis_size(1) * axis_size(2) * axis_size(3)

      allocate(density     (axis_size(1), axis_size(2), axis_size(3)) )
      allocate(temperature (axis_size(1), axis_size(2), axis_size(3)) )

! Find out what this field is called. VH_DATA will be set for files written by the VH-1 hydro code
      call ftgkys(unit,"VH_DATA",record,comment,status)
      if (status==0) then
         write(message,*) "Setting density from field with VH_DATA= ", trim(record)
         call writeInfo(message,TRIVIAL)
      else
         call writeWarning("Did not find VH_DATA keyword. Will assume this is density")
      endif

! Check the bitpix value (-32 is single precision)
      call ftgidt(unit,bitpix,status)
      if (bitpix == -32) then
         call ftgpve(unit,group,1,npixels,1e-33,density(:,:,:),anynull,status)
      else
         write(*,*) "bitpix=", bitpix
         call writeFatal("FIX ME: cannot handle bitpix != -32")
      endif

! Write out some information about the data we have just read as a basic check 
      mindata = minval(density(:,:,:))
      maxdata = maxval(density(:,:,:))
      write(message,*) "Minimum data value= ", mindata
      call writeInfo(message,FORINFO)
      write(message,*) "Maximum data value= ", maxdata
      call writeInfo(message,FORINFO)

! Move to the next extension
      call ftmrhd(unit,1,hdutype,status)
      ! Find out what this field is called. VH_DATA will be set for files written by the VH-1 hydro code
      call ftgkys(unit,"VH_DATA",record,comment,status)
      if (status==0) then
         write(message,*) "Reading pressure from field with VH_DATA= ", trim(record)
         call writeInfo(message,TRIVIAL)
      else
         call writeWarning("Did not find VH_DATA keyword. Will assume this is pressure")
      endif

! Read in the data. Assume bitpix=-32 as before
      call ftgpve(unit,group,1,npixels,1e-33,temperature(:,:,:),anynull,status)

! Write out some information about the data we have just read as a basic check 
      mindata = minval(temperature(:,:,:))
      maxdata = maxval(temperature(:,:,:))
      write(message,*) "Minimum data value= ", mindata
      call writeInfo(message,FORINFO)
      write(message,*) "Maximum data value= ", maxdata
      call writeInfo(message,FORINFO)

! Now convert pressure to temperature 

      if (idealGas) then 
         call writeInfo("Converting pressure to temperature using ideal gas EOS",TRIVIAL)
         temperature = temperature / (density * real(Rgas,si))
         mindata = minval(temperature(:,:,:))
         maxdata = maxval(temperature(:,:,:))
         write(message,*) "Minimum data value= ", mindata, " K"
         call writeInfo(message,FORINFO)
         write(message,*) "Maximum data value= ", maxdata, " K"
         call writeInfo(message,FORINFO)
      else
         call writeInfo("Converting pressure to temperature using look-up table",TRIVIAL)

! Calculate temperature from p/rho (temperature array holds p at this point).
         temperature(:,:,:) = temperature(:,:,:)/density(:,:,:)
         do k=1, axis_size(3)
            do j=1, axis_size(2)
               do i=1, axis_size(1)
                  if (temperature(i,j,k).lt.pd(1)) then 
                     temperature(i,j,k) = 10**pdlogt(1)
                  else if (temperature(i,j,k).gt.pd(npd)) then
                     temperature(i,j,k) = 10**pdlogt(npd)
                  else 
                     do n=1,npd
                        if(temperature(i,j,k).ge.pd(n) .AND. temperature(i,j,k).lt.pd(n+1)) then
                           grad1 = (pdlogt(n+1)-pdlogt(n)) / (log10(pd(n+1))-log10(pd(n)))
                           temperature(i,j,k)=10**(pdlogt(n)+grad1*(log10(temperature(i,j,k))-log10(pd(n))))
                        end if
                     end do
                  endif
               end do
            end do
         end do
      end if

! Close the file and free the LUN
      call ftclos(unit, status)
      call ftfiou(unit, status)

! If there were any errors then print them out
      call printFitsError(status)

! Now set up the co-ordinate axes. For now we're assuming they are not held in the FITS file. 
      call setup_axes()

    end subroutine read_fits_file_for_grid

!-------------------------------------------------------------------------------
    subroutine read_lutable
      real :: dummy 
      integer :: line

      call writeInfo ("Reading look-up table from "//lu_table, FORINFO)
      open(36,file=lu_table,status='old')
      do line=1,npd
         read(36,*) pdlogt(line),dummy, dummy, dummy, dummy, dummy, pd(line), dummy
      end do
      close(36) 
      
    end subroutine read_lutable

!-------------------------------------------------------------------------------

    subroutine setup_axes
      use constants_mod

      character(len=80) :: message
      integer :: i

      call writeInfo("Setting up axes for a 30pc box", IMPORTANT)

! Minimum and maximum extent of the domain
      xMin= -15.0 * pcToCm
      xMax = 15.0 * pcToCm ! assume a 30pc domain

      write (message,*) "Minimum x= ", xMin, " cm"
      call writeInfo(message, TRIVIAL)
      write (message,*) "Maximum x= ", xMax, " cm"
      call writeInfo(message, TRIVIAL)

! Cell size
      dx = abs(xMax-xMin) / real(axis_size(1),db)
      write (message,*) "Cell size dx= ", dx, " cm"
      call writeInfo(message, TRIVIAL)

      allocate(xaxis(axis_size(1)), yaxis(axis_size(2)), zaxis(axis_size(3)) )

! Set up cell centred positions in the xaxis array
! Values are in cm at this point
      do i=1, axis_size(1)
         xaxis(i) = ( ( real(i)-0.5 ) * dx ) + xMin
      enddo

! Convert to Torus units
      xaxis(:) = xaxis(:) * 1.0d-10
      dx = dx  * 1.0d-10
      xMin = xMin * 1.0d-10
      xMax = xMax * 1.0d-10

! Assume the domain is a cube
      yaxis(:) = xaxis(:)
      dy = dx
      yMin = xMin
      yMax = xMax

      zaxis(:) = xaxis(:)
      dz = dx
      zMin = xMin
      zMax = xMax

    end subroutine setup_axes

!-------------------------------------------------------------------------------

    subroutine assign_from_fitsfile(thisOctal, subcell)
      use octal_mod

      TYPE(OCTAL), intent(inout) :: thisOctal
      integer, intent(in) :: subcell

      TYPE(vector) :: thisCentre
      integer :: thisI, thisJ, thisK

      thisCentre = subcellcentre(thisOctal, subcell)

      thisI = int( ( (thisCentre%x-xMin) / dx) + 0.5 )
      thisJ = int( ( (thisCentre%y-yMin) / dy) + 0.5 )
      thisK = int( ( (thisCentre%z-zMin) / dz) + 0.5 )

      thisI = max( min(thisI,axis_size(1)), 1)
      thisJ = max( min(thisJ,axis_size(2)), 1)
      thisK = max( min(thisK,axis_size(3)), 1)

      thisOctal%rho(subcell) = density(thisI, thisJ, thisK)
      thisOctal%temperature(subcell) = temperature(thisI, thisJ, thisK)
      thisOctal%etaCont(subcell) = 0.
      thisOctal%inFlow(subcell) = .true.
      thisOctal%velocity = VECTOR(0.,0.,0.)
      thisOctal%biasCont3D = 1.
      thisOctal%etaLine = 1.e-30

    end subroutine assign_from_fitsfile

!-------------------------------------------------------------------------------

    subroutine deallocate_gridFromFitsFile

      if (allocated (density))     deallocate(density)
      if (allocated (temperature)) deallocate(temperature)
      if (allocated (xaxis))       deallocate(xaxis)
      if (allocated (yaxis))       deallocate(yaxis)
      if (allocated (zaxis))       deallocate(zaxis)

    end subroutine deallocate_gridFromFitsFile

end module gridFromFitsFile
#endif
