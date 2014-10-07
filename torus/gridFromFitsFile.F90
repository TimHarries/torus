! Module for setting up a Torus grid from data held in a FITS file

module gridFromFitsFile

#ifdef USECFITSIO
  use kind_mod 
  use messages_mod
  use constants_mod
  implicit none

  public :: read_fits_file_for_grid, assign_from_fitsfile, setGridFromFitsParameters, checkFitsSplit
  private :: setup_axes, read_lutable

  logical, private, save :: isRequired = .false.

! Filename and type (VH1, pion)
  character(len=80), private, save :: filename="none"
  character(len=80), private, save :: filetype

! Arrays to hold data from fits files
! Number of pixels
  integer, private, save :: axis_size(3)
! Axes
  real(double), private, save, allocatable :: xaxis(:), yaxis(:), zaxis(:)
  real(double), private, save :: xMin, xMax, dx
  real(double), private, save :: yMin, yMax, dy
  real(double), private, save :: zMin, zMax, dz

! Size of the x-axis
  real(double), parameter :: xSize   = 30.0*pctocm
  real(double), parameter :: xCentre =  0.0*pctocm

! Values
  real(single), private, save, allocatable :: density(:,:,:), temperature(:,:,:)
  real(double), private, save, allocatable :: density_double(:,:,:), temperature_double(:,:,:), ionfrac_double(:,:,:), &
       tr1_double(:,:,:)

! If true use the ideal gas EOS to calculate temperature
  logical, parameter, private :: idealGas=.false.

! File for look-up table for temperature conversion. 
  character(len=*), parameter, private :: lu_table="tempdependent_mu_and_gamma.dat"
  integer, parameter, private :: npd = 91
  real, private, save :: pdlogt(npd), pd(npd)

! Grid geometry
  logical, private, save :: amr2d, amr3d

  contains

! Called from inputs_mod to set parameters for this module
    subroutine setGridFromFitsParameters(infile,twod,threed)
      character(len=80) :: infile
      logical, intent(in) :: twod, threed

      isRequired = .true.
      filename=infile
      amr2d = twod
      amr3d = threed

      if (amr2d.and.amr3d) call writeWarning("Both 2D and 3D are set")
      if (( .not.amr2d) .and. (.not.amr3d) ) call writeWarning("Neither 2D nor 3D is set")

    end subroutine setGridFromFitsParameters

    subroutine read_fits_file_for_grid

      use fits_utils_mod, only: printFitsError

! FITS file variables and parameters
      integer :: status, unit, blocksize
      integer, parameter :: readwrite=0 ! Open file read only 
      integer, parameter :: group=0
      character(len=80) :: record, comment
      logical :: found_file


      call writeInfo ("Initialising grid from a FITS file", IMPORTANT)

      inquire (file=filename, exist=found_file)
      if (.not. found_file) then
         call writeFatal("The file called "//trim(filename)//" does not exist")
      endif
      call writeInfo ("Reading FITS file "//filename, FORINFO)
      status = 0

! Get a free LUN then open the file
      call ftgiou(unit, status)
      call ftopen(unit, filename,readwrite,blocksize,status)

      call ftgkys(unit,"VH_DATA",record,comment,status)
      if (status==0) then
         fileType = "VH1"
      else
         fileType = "pion"
      endif

! Close the file and free the LUN
      call ftclos(unit, status)
      call ftfiou(unit, status)

! If there were any errors then print them out
      call printFitsError(status)

      select case(fileType)
         case("VH1")
            call read_fits_file_for_grid_vh1
         case("pion")
            call read_fits_file_for_grid_pion
         case DEFAULT
            call writeFatal("Fits file type not recognised")
            stop
      end select


    end subroutine read_fits_file_for_grid

!-------------------------------------------------------------------------------
! Reads in a FITS file which will be used to set up the Torus AMR grid
! Based on the FITS cookbook
! D. Acreman, October 2012

    subroutine read_fits_file_for_grid_vh1

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
      logical :: anynull, found_file
      real(single) :: maxdata, mindata

      if (.not.idealGas) call read_lutable

      call writeInfo ("Initialising grid from a FITS file", IMPORTANT)

      inquire (file=filename, exist=found_file)
      if (.not. found_file) then
         call writeFatal("The file called "//trim(filename)//" does not exist")
      endif
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
npd_loop:            do n=1,npd
                        if(temperature(i,j,k).ge.pd(n) .AND. temperature(i,j,k).lt.pd(n+1)) then
                           grad1 = (pdlogt(n+1)-pdlogt(n)) / (log10(pd(n+1))-log10(pd(n)))
                           temperature(i,j,k)=10**(pdlogt(n)+grad1*(log10(temperature(i,j,k))-log10(pd(n))))
                           exit npd_loop
                        end if
                     end do npd_loop
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

    end subroutine read_fits_file_for_grid_vh1

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Reads in a FITS file which will be used to set up the Torus AMR grid
! Based on the FITS cookbook
! D. Acreman, October 2012

    subroutine read_fits_file_for_grid_pion

      use fits_utils_mod, only: printFitsError

      character(len=80) :: message

! FITS file variables and parameters
      integer :: status, unit, blocksize
      integer, parameter :: readwrite=0 ! Open file read only 
      integer, parameter :: group=0
      integer :: naxis, npixels, hdutype, i
      character(len=80) :: record, comment
      logical :: anynull
      real(single) :: xmin0, xmin1, xmin2
      real(single) :: xmax0, xmax1, xmax2
      real(double), parameter :: nullvall = 1.d-33


      call writeInfo ("Initialising grid from a pion FITS file", IMPORTANT)
      call writeInfo ("Reading FITS file "//filename, FORINFO)
      status = 0

! Get a free LUN then open the file
      call ftgiou(unit, status)
      call ftopen(unit, filename,readwrite,blocksize,status)

! Find out how many axes there are
      call ftgkyj(unit,"GRIDNDIM",naxis,comment,status)
      write(message,*) "There are ", naxis, " axes"
      call writeInfo(message,TRIVIAL)

! Find out the size of the axis
      axis_size(:)=1
      call ftgkyj(unit,"NGRID0",axis_size(1),comment,status)
      call ftgkyj(unit,"NGRID1",axis_size(2),comment,status)
      call ftgkyj(unit,"NGRID2",axis_size(3),comment,status)
      write(message,*) "Axis sizes: ", axis_size
      call writeInfo(message,TRIVIAL)
      npixels = axis_size(1) * axis_size(2) * axis_size(3)

      allocate(density_double     (axis_size(1), axis_size(2), axis_size(3)) )
      allocate(temperature_double (axis_size(1), axis_size(2), axis_size(3)) )
      allocate(ionfrac_double (axis_size(1), axis_size(2), axis_size(3)) )
      allocate(tr1_double (axis_size(1), axis_size(2), axis_size(3)) )

      call ftgkye(unit,"XMIN0",xmin0,comment,status)
      call ftgkye(unit,"XMIN1",xmin1,comment,status)
      call ftgkye(unit,"XMIN2",xmin2,comment,status)

      call ftgkye(unit,"XMAX0",xmax0,comment,status)
      call ftgkye(unit,"XMAX1",xmax1,comment,status)
      call ftgkye(unit,"XMAX2",xmax2,comment,status)

      allocate(xAxis(1:axis_size(1)))
      allocate(yAxis(1:axis_size(2)))
      allocate(zAxis(1:axis_size(3)))

      dx = (xMax0 - xMin0)/real(axis_size(1))
      dy = (xMax1 - xMin1)/real(axis_size(2))
      dz = (xMax2 - xMin2)/real(axis_size(3))
      do i = 1, axis_size(1)
         xAxis(i) = xMin0 + dx * real(i-1)
      enddo
      do i = 1, axis_size(2)
         yAxis(i) = xMin1 + dy * real(i-1)
      enddo
      do i = 1, axis_size(3)
         zAxis(i) = xMin2 + dz * real(i-1)
      enddo

      xAxis = xAxis / 1.e10
      yAxis = yAxis / 1.e10
      zAxis = zAxis / 1.e10
      dx = dx / 1.e10
      dy = dy / 1.e10
      dz = dz / 1.e10

      if (writeoutput) then
         write(*,*) "xAxis runs from ",xAxis(1), xAxis(axis_size(1))
         write(*,*) "yAxis runs from ",yAxis(1), xAxis(axis_size(2))
         write(*,*) "zAxis runs from ",zAxis(1), xAxis(axis_size(3))
      endif
! Move to the next extension
      call ftmrhd(unit,1,hdutype,status)
      call ftgkys(unit,"EXTNAME",record,comment,status)
      if (status==0) then
         write(message,*) "Reading density from field with EXTNAME= ", trim(record)
         call writeInfo(message,TRIVIAL)
      else
         call writeWarning("Did not find EXTNAME keyword. Will assume this is density")
      endif

! Read in the data.
      call ftgpvd(unit,group,1,npixels,nullvall,density_double(:,:,:),anynull,status)


      call FTMNHD(unit, 0, "TR0", 0, status)
      call ftgkys(unit,"EXTNAME",record,comment,status)
      if (status==0) then
         write(message,*) "Reading ionization fraction from field with EXTNAME= ", trim(record)
         call writeInfo(message,TRIVIAL)
      else
         call writeFatal("Could not find ionization fraction HDU")
         stop
      endif

! Read in the data.
      call ftgpvd(unit,group,1,npixels,nullvall,ionfrac_double(:,:,:),anynull,status)


      call FTMNHD(unit, 0, "TR1", 0, status)
      call ftgkys(unit,"EXTNAME",record,comment,status)
      if (status==0) then
         write(message,*) "Reading TR1 from field with EXTNAME= ", trim(record)
         call writeInfo(message,TRIVIAL)
      else
         call writeFatal("Could not find TR1 in HDU")
         stop
      endif

! Read in the data.
      call ftgpvd(unit,group,1,npixels,nullvall,tr1_double(:,:,:),anynull,status)


      call FTMNHD(unit, 0, "Eint", 0, status)
      call ftgkys(unit,"EXTNAME",record,comment,status)
      if (status==0) then
         write(message,*) "Reading temperature from field with EXTNAME= ", trim(record)
         call writeInfo(message,TRIVIAL)
      else
         call writeFatal("Could not find temperature HDU")
         stop
      endif

! Read in the data.
      call ftgpvd(unit,group,1,npixels,nullvall,temperature_double(:,:,:),anynull,status)


! Close the file and free the LUN
      call ftclos(unit, status)
      call ftfiou(unit, status)

! If there were any errors then print them out
      call printFitsError(status)
    end subroutine read_fits_file_for_grid_pion

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

      call writeInfo("Setting up axes for a cubic box with identical x,y and z axes", IMPORTANT)

! Minimum and maximum extent of the domain
      xMin=  xCentre - (0.5 * xSize)
      xMax = xCentre + (0.5 * xSize)

      write (message,*) "Minimum x= ", xMin/pctocm, " pc"
      call writeInfo(message, TRIVIAL)
      write (message,*) "Maximum x= ", xMax/pctocm, " pc"
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

      if (allocated(density_double)) then
         call assign_from_fitsfile_interp(thisOctal, subcell)
      else
         call assign_from_fitsfile_nointerp(thisOctal, subcell)
      endif
    end subroutine assign_from_fitsfile

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------

    subroutine assign_from_fitsfile_nointerp(thisOctal, subcell)
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

    end subroutine assign_from_fitsfile_nointerp

!-------------------------------------------------------------------------------
    logical function checkFitsSplit(thisOctal, subcell)
      use octal_mod

      TYPE(OCTAL), intent(inout) :: thisOctal
      integer, intent(in) :: subcell
      type(VECTOR) :: rVec

      select case(fileType)
         case("VH1")
            call checkFitsSplit_vh1
         case("pion")
            call checkFitsSplit_pion
         case DEFAULT
            call writeFatal("Fits file type not recognised")
            stop
      end select

    contains 

      subroutine checkFitsSplit_vh1
        integer :: thisI, thisJ, thisK
        real(double) :: cell_temperature

        rvec = subcellCentre(thisOctal, subcell)
        thisI = int( ( (rVec%x-xMin) / dx) + 0.5 )
        thisJ = int( ( (rVec%y-yMin) / dy) + 0.5 )
        thisK = int( ( (rVec%z-zMin) / dz) + 0.5 )

        thisI = max( min(thisI,axis_size(1)), 1)
        thisJ = max( min(thisJ,axis_size(2)), 1)
        thisK = max( min(thisK,axis_size(3)), 1)

! Split where temperature is below 1500K i.e where dust will be present
        cell_temperature = temperature(thisI, thisJ, thisK)
        checkFitsSplit = (cell_temperature < 1500.0) 

      end subroutine checkFitsSplit_vh1

      subroutine checkFitsSplit_pion
        logical :: inPionDomain
        rVec = thisOctal%centre
        if (amr2d) then
           inPionDomain = (rVec%x >= yAxis(1)).and.(rVec%x <= yAxis(axis_size(2))).and. &
                (rVec%z >= xAxis(1)).and.(rVec%z <= xAxis(axis_size(1)))
        else if (amr3d) then
           inPionDomain = (rVec%x >= xAxis(1)).and.(rVec%x <= xAxis(axis_size(1))).and. &
                (rVec%y >= xAxis(2)).and.(rVec%y <= yAxis(axis_size(2))).and. &
                (rVec%z >= xAxis(3)).and.(rVec%z <= zAxis(axis_size(3)))
        endif
        checkFitsSplit = .false.
        if (inPionDomain) then
           if (amr2d) then
              if (thisOctal%subcellSize > 0.5d0*min(dx,dy)) checkFitsSplit = .true.
           endif
        endif
      end subroutine checkFitsSplit_pion

    end function checkFitsSplit

!-------------------------------------------------------------------------------

    subroutine assign_from_fitsfile_interp(thisOctal, subcell)
      use octal_mod
      use utils_mod, only : locate


      TYPE(OCTAL), intent(inout) :: thisOctal
      integer, intent(in) :: subcell

      TYPE(vector) :: rVec
      integer :: thisI, thisJ, thisK
      real :: u, v, w 
      logical :: inPionDomain
      rVec = subcellcentre(thisOctal, subcell)
      thisOctal%rho(subcell) = 1.d-30
      inPionDomain = .false.
      if (amr2d) then
         inPionDomain = (rVec%z >= xAxis(1)).and.(rVec%z <= xAxis(axis_size(1))).and. &
              (rVec%x >= yAxis(1)).and.(rVec%x <= yAxis(axis_size(2)))
      else if (amr3d) then
         inPionDomain = (rVec%x >= xAxis(1)).and.(rVec%x <= xAxis(axis_size(1))).and. &
              (rVec%y >= xAxis(2)).and.(rVec%y <= yAxis(axis_size(2))).and. &
              (rVec%z >= xAxis(3)).and.(rVec%z <= zAxis(axis_size(3)))
      endif


      if (inPionDomain) then
         if (amr3d) then
            call locate(xAxis, axis_size(1), rVec%x, thisI)
            call locate(yAxis, axis_size(2), rVec%y, thisJ)
            call locate(zAxis, axis_size(3), rVec%z, thisK)
            u = real((rVec%x - xAxis(thisI))/(xAxis(thisI+1)-xAxis(thisI)))
            v = real((rVec%y - yAxis(thisJ))/(yAxis(thisJ+1)-yAxis(thisJ)))
            w = real((rVec%z - zAxis(thisK))/(zAxis(thisK+1)-zAxis(thisK)))
            thisOctal%rho(subcell) =   density_double(thisI  ,  thisJ  ,thisK  ) * (1.d0-u)*(1.d0-v)*(1.d0-w) &
                                    +  density_double(thisI+1,  thisJ  ,thisK  ) * (     u)*(1.d0-v)*(1.d0-w) &
                                    +  density_double(thisI  ,  thisJ+1,thisK  ) * (1.d0-u)*(     v)*(1.d0-w) &
                                    +  density_double(thisI+1,  thisJ+1,thisK  ) * (     u)*(     v)*(1.d0-w) &
                                    +  density_double(thisI  ,  thisJ  ,thisK+1) * (1.d0-u)*(1.d0-v)*(     w) &
                                    +  density_double(thisI+1,  thisJ  ,thisK+1) * (     u)*(1.d0-v)*(     w) &
                                    +  density_double(thisI  ,  thisJ+1,thisK+1) * (1.d0-u)*(     v)*(     w) &
                                    +  density_double(thisI+1,  thisJ+1,thisK+1) * (     u)*(     v)*(     w)  
            thisOctal%temperature(subcell) =   real(temperature_double(thisI  ,  thisJ  ,thisK  ) * (1.d0-u)*(1.d0-v)*(1.d0-w) &
                                    +  temperature_double(thisI+1,  thisJ  ,thisK  ) * (     u)*(1.d0-v)*(1.d0-w) &
                                    +  temperature_double(thisI  ,  thisJ+1,thisK  ) * (1.d0-u)*(     v)*(1.d0-w) &
                                    +  temperature_double(thisI+1,  thisJ+1,thisK  ) * (     u)*(     v)*(1.d0-w) &
                                    +  temperature_double(thisI  ,  thisJ  ,thisK+1) * (1.d0-u)*(1.d0-v)*(     w) &
                                    +  temperature_double(thisI+1,  thisJ  ,thisK+1) * (     u)*(1.d0-v)*(     w) &
                                    +  temperature_double(thisI  ,  thisJ+1,thisK+1) * (1.d0-u)*(     v)*(     w) &
                                    +  temperature_double(thisI+1,  thisJ+1,thisK+1) * (     u)*(     v)*(     w)  )
!            thisOctal%ionfrac(subcell,2) =   ionfrac_double(thisI  ,  thisJ  ,thisK  ) * (1.d0-u)*(1.d0-v)*(1.d0-w) &
!                                    +  ionfrac_double(thisI+1,  thisJ  ,thisK  ) * (     u)*(1.d0-v)*(1.d0-w) &
!                                    +  ionfrac_double(thisI  ,  thisJ+1,thisK  ) * (1.d0-u)*(     v)*(1.d0-w) &
!                                    +  ionfrac_double(thisI+1,  thisJ+1,thisK  ) * (     u)*(     v)*(1.d0-w) &
!                                    +  ionfrac_double(thisI  ,  thisJ  ,thisK+1) * (1.d0-u)*(1.d0-v)*(     w) &
!                                    +  ionfrac_double(thisI+1,  thisJ  ,thisK+1) * (     u)*(1.d0-v)*(     w) &
!                                    +  ionfrac_double(thisI  ,  thisJ+1,thisK+1) * (1.d0-u)*(     v)*(     w) &
!                                    +  ionfrac_double(thisI+1,  thisJ+1,thisK+1) * (     u)*(     v)*(     w)  
!            thisOctal%ionFrac(subcell,1) = 1.d0 - thisOctal%ionFrac(subcell,2)
         else if (amr2d) then
            call locate(xAxis, axis_size(1), rVec%z, thisI)
            call locate(yAxis, axis_size(2), rVec%x, thisJ)
            
            u = real((rVec%z - xAxis(thisI))/(xAxis(thisI+1)-xAxis(thisI)))
            v = real((rVec%x - yAxis(thisJ))/(yAxis(thisJ+1)-yAxis(thisJ)))
            thisOctal%rho(subcell) =     density_double(thisI,thisJ,1) * (1.d0-u)*(1.d0-v) &
                                    +  density_double(thisI+1,thisJ,1) * (     u)*(1.d0-v) &
                                   +  density_double(thisI, thisJ+1,1) * (1.d0-u)*(     v) &
                                   +  density_double(thisI+1,thisJ+1,1)* (     u)*(     v) 
            thisOctal%temperature(subcell) = real( temperature_double(thisI,thisJ,1) * (1.d0-u)*(1.d0-v) &
                                    +  temperature_double(thisI+1,thisJ,1) * (     u)*(1.d0-v) &
                                   +  temperature_double(thisI, thisJ+1,1) * (1.d0-u)*(     v) &
                                   +  temperature_double(thisI+1,thisJ+1,1)* (     u)*(     v) )
            if (.not.associated(thisOctal%dustTypeFraction)) then
               allocate(thisOctal%dustTypeFraction(1:thisOctal%maxChildren,1))
            endif
            thisOctal%dustTypeFraction(subcell,1) =  0.01d0 * (1.d0-(   tr1_double(thisI,thisJ,1) * (1.d0-u)*(1.d0-v) &
                                    +  tr1_double(thisI+1,thisJ,1) * (     u)*(1.d0-v) &
                                   +  tr1_double(thisI, thisJ+1,1) * (1.d0-u)*(     v) &
                                   +  tr1_double(thisI+1,thisJ+1,1)* (     u)*(     v) ))
            if (thisOctal%temperature(subcell) > 1.d6) thisOctal%dustTypeFraction(subcell,1) = 0.d0


         endif
      endif

      thisOctal%etaCont(subcell) = 0.
      thisOctal%inFlow(subcell) = .true.
      thisOctal%velocity = VECTOR(0.,0.,0.)
      thisOctal%biasCont3D = 1.
      thisOctal%etaLine = 1.e-30

    end subroutine assign_from_fitsfile_interp

!-------------------------------------------------------------------------------

      
    subroutine deallocate_gridFromFitsFile

      if (allocated (density))     deallocate(density)
      if (allocated (temperature)) deallocate(temperature)
      if (allocated (xaxis))       deallocate(xaxis)
      if (allocated (yaxis))       deallocate(yaxis)
      if (allocated (zaxis))       deallocate(zaxis)


      if (allocated (density_double))     deallocate(density_double)
      if (allocated (temperature_double))     deallocate(temperature_double)
      if (allocated (ionfrac_double))     deallocate(ionfrac_double)
      if (allocated (tr1_double))     deallocate(tr1_double)

    end subroutine deallocate_gridFromFitsFile
#endif
end module gridFromFitsFile

