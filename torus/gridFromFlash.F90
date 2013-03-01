! Read a Flash HDF file, based on Ross Parkin's raytrace.c
! D. Acreman February 2013

! To do: read in additional data for calculating cell positions
!        populate Torus grid with these values

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
  integer, parameter, private :: MDIM =3 ! number of dimensions
  character(len=80), private :: message

  real(kind=db), private, allocatable, save :: density(:,:,:,:)
  real(kind=db), private, allocatable, save :: temperature(:,:,:,:)
! Block centre co-ordinates
  real(kind=db), private, allocatable, save :: bc_coords(:,:)
! Block refinement level
  integer, private, allocatable, save :: lrefine(:)
! Block bounding boxes
  real(kind=db), private, allocatable, save :: boundBox(:,:,:)

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
  allocate (bc_coords(maxblocks,3))
  allocate (lrefine(maxblocks))
  allocate (boundBox(maxblocks,3,2))

! Read in the data
  call read_gridvar("dens",density)
  call read_gridvar("temp",temperature)

! Read block information
  call read_lrefine
  call read_coords
  call read_bndbox

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

!--------------------------------------------------------------------------------------

  subroutine read_lrefine

    integer(kind=HID_T) :: dataSet, dataSpace
    integer(kind=HSIZE_T) :: count(1), offset(1)
    integer(kind=HSIZE_T) :: dims(1)

    call writeInfo("Reading refinement levels",TRIVIAL)

! Open data set
    call h5Dopen_f(file_id, "refine level", dataSet, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening refine level data set")
    end if

! Get an identifier for the data space
    call h5dget_space_f(dataSet, dataSpace, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening refine level data space")
    end if

! Define hyperslab
    count(:) = (/maxblocks/)
    offset(:) = (/0/)
    call h5sselect_hyperslab_f(dataSpace, H5S_SELECT_SET_F, offset, count, error)
    if ( error /= 0 ) then
       call writeFatal("Error selecting hyperslab")
    end if

! Read the data
    call h5dread_f(dataSet, H5T_NATIVE_INTEGER, lrefine, dims, error)
    if ( error == 0 ) then
       write(message,*) "Maximum refinement level= ", maxval(lrefine(:))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum refinement level= ", minval(lrefine(:))
       call writeInfo(message,TRIVIAL)
    else 
       call writeFatal("Error reading refinement data")
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


  end subroutine read_lrefine

!--------------------------------------------------------------------------------------

  subroutine read_bndbox

! First dimension is the cell number, size=num blocks
! Second dimension is the co-ordinate (x,y,z), size=3
! Third dimension is lower/upper bound, size=2

    integer(kind=HID_T) :: dataSet, dataSpace
    integer(kind=HSIZE_T) :: count(3), offset(3)
    integer(kind=HSIZE_T) :: dims(3)

    call writeInfo("Reading block bounding boxes",TRIVIAL)

    ! Open data set
    call h5Dopen_f(file_id, "bounding box", dataSet, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening bounding box data set")
    end if

! Get an identifier for the data space
    call h5dget_space_f(dataSet, dataSpace, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening bounding box data space")
    end if

! Define hyperslab
    count(:) = (/maxblocks,3,2/)
    offset(:) = (/0,0,0/)
    call h5sselect_hyperslab_f(dataSpace, H5S_SELECT_SET_F, offset, count, error)
    if ( error /= 0 ) then
       call writeFatal("Error selecting hyperslab")
    end if

! Read the data
    call h5dread_f(dataSet, H5T_NATIVE_DOUBLE, boundBox, dims, error)
    if ( error == 0 ) then
       write(message,*) "Maximum x bound= ", maxval(boundBox(:,1,2))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum x bound= ", minval(boundBox(:,1,1))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Maximum y bound= ", maxval(boundBox(:,2,2))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum y bound= ", minval(boundBox(:,2,1))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Maximum z bound= ", maxval(boundBox(:,3,2))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum z bound= ", minval(boundBox(:,3,1))
       call writeInfo(message,TRIVIAL)
    else 
       call writeFatal("Error reading refinement data")
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

! Convert to Torus units
    boundBox(:,:,:) = boundBox(:,:,:) * 1.0e-10_db

  end subroutine read_bndbox

!--------------------------------------------------------------------------------------

  subroutine read_coords

    integer(kind=HID_T) :: dataSet, dataSpace
    integer(kind=HSIZE_T) :: count(2), offset(2)
    integer(kind=HSIZE_T) :: dims(2)

    call writeInfo("Reading coordinates",TRIVIAL)

! Open data set
    call h5Dopen_f(file_id, "coordinates", dataSet, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening coordinates data set")
    end if

! Get an identifier for the data space
    call h5dget_space_f(dataSet, dataSpace, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening coordinates data space")
    end if

! Define hyperslab
    count(:) = (/maxblocks, 2/)
    offset(:) = (/0,0/)
    call h5sselect_hyperslab_f(dataSpace, H5S_SELECT_SET_F, offset, count, error)
    if ( error /= 0 ) then
       call writeFatal("Error selecting hyperslab")
    end if

! Read the data
    call h5dread_f(dataSet, H5T_NATIVE_DOUBLE, bc_coords, dims, error)
    if ( error == 0 ) then
       write(message,*) "Maximum x value= ", maxval(bc_coords(:,1))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum x value= ", minval(bc_coords(:,1))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Maximum y value= ", maxval(bc_coords(:,2))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum y value= ", minval(bc_coords(:,2))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Maximum z value= ", maxval(bc_coords(:,3))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum z value= ", minval(bc_coords(:,3))
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

! Convert to Torus units
    bc_coords(:,:) = bc_coords(:,:) * 1.0e-10

  end subroutine read_coords

end subroutine read_flash_hdf

!--------------------------------------------------------------------------------------

subroutine assign_from_flash(thisOctal, subcell)

  use octal_mod
  use vector_mod
  implicit none

  real(db), parameter :: rho_bg=1.0e-30_db

  TYPE(OCTAL), intent(inout) :: thisOctal
  integer, intent(in) :: subcell
  integer :: i
  TYPE(VECTOR) :: thisTorusCellCentre
  integer :: lrefineClosest, iClosest
  real(double) :: xdist, ydist, zdist
  real(double) :: dx, dy, dz
  integer :: icell, jcell, kcell

  thisTorusCellCentre = subCellCentre(thisOctal, subcell)
  lrefineClosest = -1
  iClosest = -1

  do i=1, maxblocks

! Is this subcell in this block?

     if ( (boundBox(i,1,1) < thisTorusCellCentre%x).and. &
          (boundBox(i,1,2) > thisTorusCellCentre%x).and. &
          (boundBox(i,2,1) < thisTorusCellCentre%y).and. &
          (boundBox(i,2,2) > thisTorusCellCentre%y).and. &
          (boundBox(i,3,1) < thisTorusCellCentre%z).and. &
          (boundBox(i,3,2) > thisTorusCellCentre%z) ) then 

! If the block is higher refinement than any previous block containing 
! the point then use the current block
        if (lrefine(i) > lrefineClosest) then
           iClosest       = i
           lrefineClosest = lrefine(i)
        end if

     end if
  end do

  if (iClosest /= -1 ) then

     xdist = thisTorusCellCentre%x - boundBox(iClosest,1,1)
     dx    = (boundBox(iClosest,1,2) - boundBox(iClosest,1,1)) / real(NXB,db)
     icell = int(xdist / dx) + 1 

     ydist = thisTorusCellCentre%y - boundBox(iClosest,2,1)
     dy    = (boundBox(iClosest,2,2) - boundBox(iClosest,2,1)) / real(NYB,db)
     jcell = int(ydist / dy) + 1 

     zdist = thisTorusCellCentre%z - boundBox(iClosest,3,1)
     dz    = (boundBox(iClosest,3,2) - boundBox(iClosest,3,1)) / real(NZB,db)
     kcell = int(zdist / dz) + 1 

     thisOctal%rho(subcell) = density(iClosest,icell,jcell,kcell) 

  else
     thisOctal%rho(subcell) = rho_bg
  endif

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

! Stub subroutines in case Torus has been built without HDF5 support
subroutine read_flash_hdf
  write(*,*) "Torus was built without HDF support so cannot read a Flash HDF file"
  write(*,*) "Rebuild Torus with hdf=yes"
  STOP
end subroutine read_flash_hdf

subroutine assign_from_flash(thisOctal, subcell)

  use octal_mod
  use vector_mod
  implicit none

  real(db), parameter :: rho_bg=1.0e-23_db

  TYPE(OCTAL) :: thisOctal
  integer     :: subcell

  thisOctal%rho(subcell) = 1.0e-33_db
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

#endif

subroutine deallocate_gridFromFlash

  if (allocated(density))     deallocate (density)
  if (allocated(temperature)) deallocate (temperature)
  if (allocated(bc_coords))   deallocate (bc_coords)
  if (allocated(lrefine))     deallocate (lrefine)
  if (allocated(boundBox))    deallocate (boundBox)

end subroutine deallocate_gridFromFlash

end module gridFromFlash
