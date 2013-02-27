! Read a Flash HDF file, based on Ross Parkin's raytrace.c
! D. Acreman February 2013

! To do: read in additional data for calculating cell positions
!        populate Torus grid with these values
!        de-allocate arrays when they are finished with

module gridFromFlash

  use kind_mod
  use messages_mod

! Public module variables
  logical :: flashFileRequired=.true.

! Private module variables
  character(len=*), parameter, private :: infile="bow_hdf5_chk_0048"
  integer, parameter, private :: maxblocks=11056
  integer, parameter, private :: NXB = 8
  integer, parameter, private :: NYB = 8
  integer, parameter, private :: NZB = 8
  character(len=80), private :: message

  real(kind=db), private, allocatable, save :: density(:,:,:,:)
  real(kind=db), private, allocatable, save :: temperature(:,:,:,:)
  real(kind=db), private, allocatable, save :: pressure(:,:,:,:)

contains

#ifdef USEHDF
  subroutine read_flash_hdf

  USE HDF5

  implicit none

  integer        :: error
  integer(HID_T) :: file_id

  call writeInfo("")
  call writeInfo("Reading a Flash HDF file",FORINFO)

! Initialise HDF interface
  CALL h5open_f (error)

! Open existing file read only
  CALL h5fopen_f (infile, H5F_ACC_RDONLY_F, file_id, error)
  if ( error == 0 ) then 
     call writeInfo("Opened file "//infile,FORINFO)
  else
     call writeFatal("Error opening "//infile)
  endif

  allocate (density(maxblocks, NXB, NYB, NZB))
  allocate (temperature(maxblocks, NXB, NYB, NZB))
  allocate (pressure(maxblocks, NXB, NYB, NZB))

! Read in the data
  call read_gridvar("dens",density)
  call read_gridvar("temp",temperature)
  call read_gridvar("pres", pressure)

! Will need to read some more things here to work out the cell co-ordintes

! Close the file
  CALL h5fclose_f(file_id, error)

! Terminate HDF interface
  CALL h5close_f(error)

  call writeInfo("Finished reading Flash file",FORINFO)
  call writeInfo("")

contains

!----------------------------------------------------------------------------

  subroutine read_gridvar(varName,buf)

    character(len=*), intent(in) :: varName
    real(kind=db), intent(out) :: buf(maxblocks, NXB, NYB, NZB)

    integer(kind=HID_T) :: dataSet, dataSpace
    integer(kind=HSIZE_T) :: count(4), offset(4)
    integer(kind=HSIZE_T) :: dims(4)

    call writeInfo("Reading grid variable "//varName,TRIVIAL)

    ! Open data set
    call h5Dopen_f(file_id, varName, dataSet, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening data set "//varName)
    end if

! Get an identifier for the data space
    call h5dget_space_f(dataSet, dataSpace, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening data space")
    end if

! Define hyperslab
    count = (/maxblocks, NXB, NYB, NZB/)
    offset(:) = 0
    call h5sselect_hyperslab_f(dataSpace, H5S_SELECT_SET_F, offset, count, error)
    if ( error /= 0 ) then
       call writeFatal("Error selecting hyperslab")
    end if

! Read the data
    call h5dread_f(dataSet, H5T_NATIVE_DOUBLE, buf, dims, error)
    if ( error == 0 ) then
       write(message,*) "Maximum data value= ", maxval(buf)
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum data value= ", minval(buf)
       call writeInfo(message,TRIVIAL)
    else 
       call writeFatal("Error reading data")
    end if

! Close data space
    call h5Sclose_f(dataSpace, error)
    if ( error /= 0 ) then
       call writeFatal("Error closing data space")
    end if

! Close data set
    call h5dclose_f(dataSet, error)
    if ( error /= 0 ) then
       call writeFatal("Error closing data set")
    end if

  end subroutine read_gridvar

end subroutine read_flash_hdf

subroutine assign_from_flash(thisOctal, subcell)

  use octal_mod
  use vector_mod
  implicit none

  real(db), parameter :: rho_bg=1.0e-23_db

  TYPE(OCTAL), intent(inout) :: thisOctal
  integer, intent(in) :: subcell
  integer, save :: i=1

  thisOctal%rho(subcell) = density(i,1,1,1) !rho_bg
  i = i+1
  if (i>maxblocks) i=1

  thisOctal%dustTypeFraction(subcell,:) = 0.0

  thisOctal%temperature(subcell) = 10000.
  thisOctal%etaCont(subcell) = 0.
  thisOctal%nh(subcell) = thisOctal%rho(subcell) / mHydrogen
  thisOctal%ne(subcell) = thisOctal%nh(subcell)
  thisOctal%nhi(subcell) = 1.e-8
  thisOctal%nhii(subcell) = thisOctal%ne(subcell)
  thisOctal%inFlow(subcell) = .true.
  thisOctal%velocity = VECTOR(0.,0.,0.)
  thisOctal%biasCont3D = 1.
  thisOctal%etaLine = 1.e-30

end subroutine assign_from_flash

#else

! Stub subroutine in case Torus has been built without HDF5 support
subroutine read_flash_hdf
  write(*,*) "Torus was built without HDF support so cannot read a Flash HDF file"
  write(*,*) "Rebuild Torus with hdf=yes"
  STOP
end subroutine read_flash_hdf

#endif


end module gridFromFlash
