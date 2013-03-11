! Read a Flash HDF file, written with help from Ross Parkin's raytrace.c
! D. Acreman February 2013

module gridFromFlash

  use kind_mod
  use messages_mod

  public :: setGridFromFlashParameters, assign_from_flash, flashFileRequired, &
       read_flash_hdf, deallocate_gridfromflash

  private

  logical :: isRequired=.false.

  character(len=80), private  :: infile
  real(double), private       :: ySlice
  logical, private            :: reflectYaxis=.false.

  integer, private            :: maxblocks
  integer, parameter, private :: NXB = 8
  integer, parameter, private :: NYB = 8
  integer, parameter, private :: NZB = 8
  integer, parameter, private :: MDIM =3 ! number of dimensions
  character(len=80), private  :: message

  real(kind=db), private, allocatable, save :: density(:,:,:,:)
  real(kind=db), private, allocatable, save :: temperature(:,:,:,:)
  real(kind=db), private, save :: minRho, minTem

! Block refinement level
  integer, private, allocatable, save :: lrefine(:)
  integer, private, save :: maxRefine ! maximum level of refinement
! Block bounding boxes
  real(kind=db), private, allocatable, save :: boundBox(:,:,:)

contains

! Accessor function
  logical function flashFileRequired()
    flashFileRequired = isRequired
  end function flashFileRequired

! Set module variables from values the parameters file. Called from inputs_mod.
  subroutine setGridFromFlashParameters(flashfilename, numblocks, slice, doReflectY)

    implicit none 

    integer, intent(in) :: numblocks
    character(len=80), intent(in) :: flashfilename
    real(double), intent(in) :: slice
    logical, intent(in) :: doReflectY

    if (numblocks > 0) then
       maxBlocks = numBlocks
    else
       write(message,*) "Number of blocks must be positive.", numBlocks
       call writeFatal(message)
    endif

    infile = flashfilename

! It is currently assumed that the slice is taken in the y-direction
    yslice = slice

! Currently only the y-axis can be reflected
    reflectYaxis = doReflectY

! The read subroutine is actvated when setGridFromFlashParameters is called
    isRequired=.true.

  end subroutine setGridFromFlashParameters

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

  allocate (density(NXB, NYB, NZB,maxblocks))
  allocate (temperature(NXB, NYB, NZB,maxblocks))
  allocate (lrefine(maxblocks))
  allocate (boundBox(2,3,maxblocks))

! Read in the data
  call read_gridvar("dens",density,minRho)
  call read_gridvar("temp",temperature,minTem)

! Read block information
  call read_lrefine
  call read_bndbox

! Close the file
  CALL h5fclose_f(file_id, error)

! Terminate HDF interface
  CALL h5close_f(error)

  call writeInfo("Finished reading Flash file",FORINFO)
  call writeInfo("")

contains

!----------------------------------------------------------------------------

  subroutine read_gridvar(varName,buf,minValue)

    character(len=*), intent(in) :: varName
    real(kind=db), intent(out) :: buf(maxblocks, NXB, NYB, NZB)
    real(kind=db), intent(out) :: minValue

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
    count = (/NXB, NYB, NZB, maxblocks/)
    offset(:) = 0
    call h5sselect_hyperslab_f(dataSpace, H5S_SELECT_SET_F, offset, count, error)
    if ( error /= 0 ) then
       call writeFatal("Error selecting hyperslab")
    end if

! Read the data
    call h5dread_f(dataSet, H5T_NATIVE_DOUBLE, buf, dims, error)
    minValue = minval(buf)
    if ( error == 0 ) then
       write(message,*) "Maximum data value= ", maxval(buf)
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum data value= ", minvalue
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
       maxRefine =  maxval(lrefine(:))
       write(message,*) "Maximum refinement level= ", maxRefine
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
    count(:) = (/2,3,maxblocks/)
    offset(:) = (/0,0,0/)
    call h5sselect_hyperslab_f(dataSpace, H5S_SELECT_SET_F, offset, count, error)
    if ( error /= 0 ) then
       call writeFatal("Error selecting hyperslab")
    end if

! Read the data
    call h5dread_f(dataSet, H5T_NATIVE_DOUBLE, boundBox, dims, error)
    if ( error == 0 ) then
       write(message,*) "Maximum x bound= ", maxval(boundBox(2,1,:))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum x bound= ", minval(boundBox(1,1,:))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Maximum y bound= ", maxval(boundBox(2,2,:))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum y bound= ", minval(boundBox(1,2,:))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Maximum z bound= ", maxval(boundBox(2,3,:))
       call writeInfo(message,TRIVIAL)
       write(message,*) "Minimum z bound= ", minval(boundBox(1,3,:))
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

end subroutine read_flash_hdf

!--------------------------------------------------------------------------------------

subroutine assign_from_flash(thisOctal, subcell)

  use octal_mod
  use vector_mod
  implicit none

  TYPE(OCTAL), intent(inout) :: thisOctal
  integer, intent(in) :: subcell
  integer :: i
  TYPE(VECTOR) :: thisTorusCellCentre
  integer :: lrefineClosest, iClosest
  real(double) :: xdist, ydist, zdist
  real(double) :: dx, dy, dz
  real(double) :: thisY ! For treating reflective y-axis
  integer :: icell, jcell, kcell

  thisTorusCellCentre = subCellCentre(thisOctal, subcell)
  lrefineClosest = -1
  iClosest = -1

  if (reflectYaxis) then 
     thisY = abs(thisTorusCellCentre%y)
  else
     thisY = thisTorusCellCentre%y
  endif

blocks:  do i=1, maxblocks

! Is this subcell in this block?
     if ( thisOctal%twoD) then

        if ( (boundBox(1,1,i) <= thisTorusCellCentre%x).and. &
             (boundBox(2,1,i) >  thisTorusCellCentre%x).and. &
             (boundBox(1,2,i) <= ySlice).and. &
             (boundBox(2,2,i) >  ySlice).and. &
             (boundBox(1,3,i) <= thisTorusCellCentre%z).and. &
             (boundBox(2,3,i) >  thisTorusCellCentre%z) ) then 

! If the block is higher refinement than any previous block containing 
! the point then use the current block
           if (lrefine(i) > lrefineClosest) then
              iClosest       = i
              lrefineClosest = lrefine(i)
! If we've found a block at the highest level of refinement then we need look 
! no further so exit the loop
              if (lrefineClosest == maxRefine) exit blocks
           end if

        end if

     elseif (thisOctal%threeD) then

        if ( (boundBox(1,1,i) <= thisTorusCellCentre%x).and. &
             (boundBox(2,1,i) >  thisTorusCellCentre%x).and. &
             (boundBox(1,2,i) <= thisy).and. &
             (boundBox(2,2,i) >  thisy).and. &
             (boundBox(1,3,i) <= thisTorusCellCentre%z).and. &
             (boundBox(2,3,i) >  thisTorusCellCentre%z) ) then 

! If the block is higher refinement than any previous block containing 
! the point then use the current block
           if (lrefine(i) > lrefineClosest) then
              iClosest       = i
              lrefineClosest = lrefine(i)
! If we've found a block at the highest level of refinement then we need look 
! no further so exit the loop
              if (lrefineClosest == maxRefine) exit blocks
           end if

        end if

     else
        call writeFatal("Octal is not 2D or 3D")
     end if
  end do blocks

! Find the closest cell within the block
  if (iClosest /= -1 ) then

     xdist = thisTorusCellCentre%x - boundBox(1,1,iClosest)
     dx    = (boundBox(2,1,iClosest) - boundBox(1,1,iClosest)) / real(NXB,db)
     icell = int(xdist / dx) + 1 

     if ( thisOctal%twoD) then
        ydist = ySlice - boundBox(1,2,iClosest)
     elseif (thisOctal%threeD) then
        ydist = thisy - boundBox(1,2,iClosest)
     else
        call writeFatal("Octal is not 2D or 3D")
     endif
     dy    = (boundBox(2,2,iClosest) - boundBox(1,2,iClosest)) / real(NYB,db)
     jcell = int(ydist / dy) + 1

     zdist = thisTorusCellCentre%z - boundBox(1,3,iClosest)
     dz    = (boundBox(2,3,iClosest) - boundBox(1,3,iClosest)) / real(NZB,db)
     kcell = int(zdist / dz) + 1 

     thisOctal%rho(subcell) = density(icell,jcell,kcell,iClosest) 
     thisOctal%temperature(subcell) = temperature(icell,jcell,kcell,iClosest)
  else
     thisOctal%rho(subcell) = minRho
     thisOctal%temperature(subcell) = minTem
  endif

! Choose where to put dust
  if (thisOctal%temperature(subcell) > 1.5e5) then
     thisOctal%dustTypeFraction(subcell,:) = 0.0
  else
     thisOctal%dustTypeFraction(subcell,:) = 1.0
  end if

! Velocity is set to zero for now, add here if required
  thisOctal%velocity = VECTOR(0.,0.,0.)

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
  implicit none

  TYPE(OCTAL) :: thisOctal
  integer     :: subcell

  thisOctal%rho(subcell) = 1.0e-33_db
  call writeFatal("Called assign_from_flash in build without HDF support")

end subroutine assign_from_flash

#endif

subroutine deallocate_gridFromFlash

  if (allocated(density))     deallocate (density)
  if (allocated(temperature)) deallocate (temperature)
  if (allocated(lrefine))     deallocate (lrefine)
  if (allocated(boundBox))    deallocate (boundBox)

end subroutine deallocate_gridFromFlash

end module gridFromFlash
