! Read a Flash HDF file, written with help from Ross Parkin's raytrace.c
! D. Acreman February 2013

module gridFromFlash

  use kind_mod
  use messages_mod

  public :: setGridFromFlashParameters, assign_from_flash, flashFileRequired, &
       read_flash_hdf, deallocate_gridfromflash, splitFlash, copy_flash_bounds, &
       assign_from_flash_silcc

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

!  ! SILCC Hermite sinks
!  integer, parameter, private :: nSinkProp=106 ! no. of sink particle properties
!  integer, parameter, private :: maxNsinks=100 ! max no. of sink particles

  real(kind=db), private, allocatable, save :: density(:,:,:,:)
  real(kind=db), private, allocatable, save :: temperature(:,:,:,:)
  real(kind=db), private, allocatable, save :: tdust(:,:,:,:)
  real(kind=db), private, allocatable, save :: ihp(:,:,:,:)
  real(kind=db), private, allocatable, save :: vx(:,:,:,:)
  real(kind=db), private, allocatable, save :: vy(:,:,:,:)
  real(kind=db), private, allocatable, save :: vz(:,:,:,:)
!  real(kind=db), private, allocatable, save :: sinkList(:,:)
#ifdef USEHDF
  real(kind=db), private, save :: minRho, minTem, minVx, minVy, minVz, minTdust, minihp
  real(kind=db), private, save :: maxRho, maxTem, maxVx, maxVy, maxVz, maxTdust, maxihp
#endif
  integer, private, save :: iBlockLoc=0

! Block refinement level
  integer, private, allocatable, save :: lrefine(:)
#ifdef USEHDF
  integer, private, save :: maxRefine ! maximum level of refinement
#endif
! Block bounding boxes
  real(kind=db), private, allocatable, save :: boundBox(:,:,:)
! block is a leaf (its children cells are not split further)
  logical, private, allocatable, save :: isLeaf(:)

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

! The read subroutine is activated when setGridFromFlashParameters is called
    isRequired=.true.

  end subroutine setGridFromFlashParameters

#ifdef USEHDF
  subroutine read_flash_hdf(silcc)

  USE HDF5

  implicit none

  integer        :: error
  integer(HID_T) :: file_id
  logical, optional :: silcc


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
  if (present(silcc)) then
     if (silcc) then
        allocate (tdust(NXB, NYB, NZB,maxblocks))
        allocate (ihp(NXB, NYB, NZB,maxblocks))
!        allocate (sinkList(nSinkProp,maxNsinks))
     endif
  endif
  allocate (vx(NXB, NYB, NZB,maxblocks))
  allocate (vy(NXB, NYB, NZB,maxblocks))
  allocate (vz(NXB, NYB, NZB,maxblocks))
  allocate (lrefine(maxblocks))
  allocate (boundBox(2,3,maxblocks))
  allocate (isLeaf(maxblocks))

! Read in the data
  call read_gridvar("dens", density,     minRho, maxRho)
  call read_gridvar("temp", temperature, minTem, maxTem)
  if (present(silcc)) then
     if (silcc) then
        call read_gridvar("tdus", tdust, minTdust, maxTdust)
        call read_gridvar("ihp /", ihp, minihp, maxihp) ! "ihp " has trailing space
!        call read_silcc_sinks(sinkList)
     endif
  endif
  call read_gridvar("velx", vx,          minVx,  maxVx)
  call read_gridvar("vely", vy,          minVy,  maxVy)
  call read_gridvar("velz", vz,          minVz,  maxVz)

! Read block information
  call read_nleaves
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

  subroutine read_gridvar(varName,buf,minValue,maxValue)

    character(len=*), intent(in) :: varName
    real(kind=db), intent(out) :: buf(maxblocks, NXB, NYB, NZB)
    real(kind=db), intent(out) :: minValue, maxValue

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
    maxValue = maxval(buf)
    if ( error == 0 ) then
       write(message,*) "Maximum data value= ", maxvalue
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

  subroutine read_nleaves

    integer(kind=HID_T) :: dataSet, dataSpace
    integer(kind=HSIZE_T) :: count(2), offset(2)
    integer(kind=HSIZE_T) :: dims(2)
    integer :: buf(15,maxblocks)
    integer :: nleaves,iblock

    call writeInfo("Reading gid for no. of leaf blocks",TRIVIAL)

! Open data set
    call h5Dopen_f(file_id, "gid", dataSet, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening gid data set")
    end if

! Get an identifier for the data space
    call h5dget_space_f(dataSet, dataSpace, error)
    if ( error /= 0 ) then
       call writeFatal("Error opening gid data space")
    end if

! Define hyperslab
    count = (/maxblocks, 15/)
    offset(:) = 0
    call h5sselect_hyperslab_f(dataSpace, H5S_SELECT_SET_F, offset, count, error)
    if ( error /= 0 ) then
       call writeFatal("Error selecting hyperslab")
    end if

! Read the data
    nleaves = 0
    call h5dread_f(dataSet, H5T_NATIVE_INTEGER, buf, dims, error)
    if ( error == 0 ) then
      ! block is a leaf if it doesn't have children, i.e. all child gid == -1
       do iblock = 1, maxblocks
          if (all(buf(8:15, iblock) == -1)) then
             nleaves = nleaves + 1
             isLeaf(iblock) = .true.
          else
             isLeaf(iblock) = .false.
          endif
       enddo
    else
       call writeFatal("Error reading gid data")
    end if

    write(message,*) "nleaves = ", nleaves
    call writeInfo(message,TRIVIAL)
    if (nleaves < 1) then
       call writeFatal("Error counting leaves")
    endif

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


  end subroutine read_nleaves

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

!--------------------------------------------------------------------------------------

! finished, should work. Check nSinkProp and maxNsinks are correct
!  subroutine read_silcc_sinks(list)
!
!    character(len=*), intent(in) :: varName
!    real(kind=db), intent(out) :: buf(maxblocks, NXB, NYB, NZB)
!    real(kind=db), intent(out) :: minValue, maxValue
!
!    integer(kind=HID_T) :: dataSet, dataSpace
!    integer(kind=HSIZE_T) :: count(4), offset(4)
!    integer(kind=HSIZE_T) :: dims(4)
!
!    call writeInfo("Reading sink list",TRIVIAL)
!
!    ! Open data set
!    call h5Dopen_f(file_id, "sinkList", dataSet, error)
!    if ( error /= 0 ) then
!       call writeFatal("Error opening data set sinkList")
!    end if
!
!! Get an identifier for the data space
!    call h5dget_space_f(dataSet, dataSpace, error)
!    if ( error /= 0 ) then
!       call writeFatal("Error opening data space")
!    end if
!
!! Define hyperslab
!    count = (/nSinkProp*maxNsinks/)
!    offset(:) = (/0/)
!    call h5sselect_hyperslab_f(dataSpace, H5S_SELECT_SET_F, offset, count, error)
!    if ( error /= 0 ) then
!       call writeFatal("Error selecting hyperslab")
!    end if
!
!! Read the data
!    call h5dread_f(dataSet, H5T_NATIVE_DOUBLE, buf, dims, error)
!    if ( error /= 0 ) then
!       call writeFatal("Error reading sink data")
!    end if
!    ! buf is 1d array with nSinkProp properties for each sink, one after the other
!    do i = 1, maxNsinks
!       j1 = 1 + (i-1)*nSinkProp
!       j2 = nSinkProp + (i-1)*nSinkProp
!       sinkList(1:nSinkProp,i) = buf(j1:j2)
!    enddo
!
!! Close data space
!    call h5Sclose_f(dataSpace, error)
!    if ( error /= 0 ) then
!       call writeFatal("Error closing data space")
!    end if
!
!! Close data set
!    call h5dclose_f(dataSet, error)
!    if ( error /= 0 ) then
!       call writeFatal("Error closing data set")
!    end if
!
!  end subroutine read_silcc_sinks

end subroutine read_flash_hdf

!--------------------------------------------------------------------------------------
! Set subcell density, temperature and velocity

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

! if this block is not a leaf block, go to the next block
     if (.not.isLeaf(i)) cycle

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

! Store velocity as fraction of c.
     thisOctal%velocity(subcell) = VECTOR( (vx(icell,jcell,kcell,iClosest) / cspeed) , &
                                           (vy(icell,jcell,kcell,iClosest) / cspeed) , &
                                           (vz(icell,jcell,kcell,iClosest) / cspeed) )

  else
     thisOctal%rho(subcell) = minRho
     thisOctal%temperature(subcell) = minTem
     thisOctal%velocity(subcell) = VECTOR(0.0_db, 0.0_db, 0.0_db)
  endif

end subroutine assign_from_flash

!--------------------------------------------------------------------------------------
! Set source array properties from SILCC sink properties

! AA unfinished...
!subroutine assign_silcc_sinks(source, nSource)
!  use source_mod
!  implicit none
!
!  do i = 1, maxNsinks
!     if (sinkexists) then
!        nSource = nSource + 1
!
!        ! sinkList is in cgs
!        source(nSource)%position%x = sinkList(1,i) / 1.d10
!        source(nSource)%position%y = sinkList(2,i) / 1.d10
!        source(nSource)%position%z = sinkList(3,i) / 1.d10
!
!        source(nSource)%mass = sinkList(60,i)
!        source(nSource)%teff = sinkList(70,i)
!
!        uvlum = sinkList(62,i)
!        source(nSource)%luminosity = uvlum
!     endif
!  enddo
!
!end subroutine assign_from_flash

!--------------------------------------------------------------------------------------

subroutine assign_from_flash_silcc(thisOctal, subcell)

  use octal_mod
  use vector_mod
  implicit none

  TYPE(OCTAL), intent(inout) :: thisOctal
  integer, intent(in) :: subcell
  integer :: i
  TYPE(VECTOR) :: thisTorusCellCentre, blockCentre
  integer :: iClosest
  real(double) :: xdist, ydist, zdist
  real(double) :: dx, dy, dz
  real(double) :: dx_cell, dy_cell, dz_cell
  integer :: icell, jcell, kcell
  logical :: subcellInsideBlock


  thisTorusCellCentre = subCellCentre(thisOctal, subcell)

  if (.not. thisOctal%threeD) then
     call writeFatal("Octal is not 3D")
  endif

  iClosest = -1

! find the leaf block containing this subcell
blocks:  do i=1, maxblocks
     if (isleaf(i)) then
        blockCentre = VECTOR((boundbox(1,1,i)+boundbox(2,1,i))/2.d0, &
                             (boundbox(1,2,i)+boundbox(2,2,i))/2.d0, &
                             (boundbox(1,3,i)+boundbox(2,3,i))/2.d0)
        dx = boundbox(2,1,i)-boundbox(1,1,i)
        dy = boundbox(2,2,i)-boundbox(1,2,i)
        dz = boundbox(2,3,i)-boundbox(1,3,i)
        subcellInsideBlock = inCell(blockCentre, thisTorusCellCentre, dx, dy, dz)
        if (subcellInsideBlock) then
           iClosest = i
           exit blocks
        endif
     endif
  end do blocks

! Find the closest Flash cell within the block
  if (iClosest /= -1 ) then
     ! TODO
     ! loop through blocks
     ! < 1 or > 8 is outside the block
     ! if outside, skip
     ! if in side, we have the cell

     xdist = thisTorusCellCentre%x - boundBox(1,1,iClosest)
     dx_cell    = dx / real(NXB,db)
     icell = int(xdist / dx_cell) + 1

     ydist = thisTorusCellCentre%y - boundBox(1,2,iClosest)
     dy_cell    = dy / real(NYB,db)
     jcell = int(ydist / dy_cell) + 1

     zdist = thisTorusCellCentre%z - boundBox(1,3,iClosest)
     dz_cell    = dz / real(NZB,db)
     kcell = int(zdist / dz_cell) + 1

     thisOctal%rho(subcell) = density(icell,jcell,kcell,iClosest)
     thisOctal%temperature(subcell) = temperature(icell,jcell,kcell,iClosest)
     thisOctal%tdust(subcell) = tdust(icell,jcell,kcell,iClosest)
     thisOctal%ionFrac(subcell,2) = ihp(icell,jcell,kcell,iClosest)

! Store velocity as fraction of c.
     thisOctal%velocity(subcell) = VECTOR( (vx(icell,jcell,kcell,iClosest) / cspeed) , &
                                           (vy(icell,jcell,kcell,iClosest) / cspeed) , &
                                           (vz(icell,jcell,kcell,iClosest) / cspeed) )

  else
     write(*,*) "torus cell centre ", thisTorusCellCentre
     write(*,*) "cell size", thisOctal%subcellsize
     call writeFatal("subcell didn't find a corresponding leaf block")
!     thisOctal%rho(subcell) = minRho
!     thisOctal%temperature(subcell) = minTem
!     thisOctal%tdust(subcell) = minTdust
!     thisOctal%velocity(subcell) = VECTOR(0.0_db, 0.0_db, 0.0_db)
  endif

end subroutine assign_from_flash_silcc

! determine whether to split this torus subcell
function splitFlash(thisOctal, subcell) result(split)
  use octal_mod
  use vector_mod
  implicit none

  TYPE(OCTAL), intent(inout) :: thisOctal
  integer, intent(in) :: subcell
  logical :: split
  integer :: i
  TYPE(VECTOR) :: thisTorusCellCentre, blockCentre
  real(double) :: dx, dy, dz, cellsize
  logical :: blockInsideSubcell
!  logical :: subcellInsideBlock

  ! this returns split=.true. until subcell size matches block size
  ! by comparing Torus ndepth with Flash refine level

  split = .false.

  thisTorusCellCentre = subCellCentre(thisOctal, subcell)

blocks:  do i=1, maxblocks
     blockCentre = VECTOR((boundbox(1,1,i)+boundbox(2,1,i))/2.d0, &
                          (boundbox(1,2,i)+boundbox(2,2,i))/2.d0, &
                          (boundbox(1,3,i)+boundbox(2,3,i))/2.d0)
     dx = boundbox(2,1,i)-boundbox(1,1,i) ! block size
     dy = boundbox(2,2,i)-boundbox(1,2,i)
     dz = boundbox(2,3,i)-boundbox(1,3,i)
     cellSize = thisOctal%subcellSize
     blockInsideSubcell = inCell(blockCentre, thisTorusCellCentre, cellSize, cellSize, cellSize)
!     subcellInsideBlock = inCell(blockCentre, thisTorusCellCentre, dx, dy, dz)

     if (blockInsideSubcell) then
!        if (thisOctal%ndepth < lrefine(i)-1) then
!           split = .true.
!           exit blocks
!        endif
        if (cellsize > min(dx,dy,dz)*1.3) then
           split = .true.
           exit blocks
        endif
     endif
  end do blocks
end function splitFlash

 logical function inCell(a, b, dx, dy, dz)
   use vector_mod
   implicit none
   type(VECTOR), intent(in) :: a, b
   real(double), intent(in) :: dx, dy, dz
   real(double) :: diffx, diffy, diffz
   real (double) :: halfdx, halfdy, halfdz
   ! calculate if a and b are within 0.5 dx of each other
   ! e.g. a is cell centre, b is block centre, dx is cell size

   halfdx = dx/2.0
   halfdy = dy/2.0
   halfdz = dz/2.0

   diffx = b%x - a%x
   diffy = b%y - a%y
   diffz = b%z - a%z

   if ( abs(diffx) <= halfdx .and. abs(diffy) <= halfdy .and. abs(diffz) <= halfdz) then
      inCell = .true.
   else
      inCell = .false.
   endif

 end function inCell

#else

! Stub subroutines in case Torus has been built without HDF5 support
subroutine read_flash_hdf(silcc)
  implicit none
  logical, optional :: silcc
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

function splitFlash(thisOctal, subcell) result(split)
  use octal_mod
  implicit none
  TYPE(OCTAL) :: thisOctal
  integer     :: subcell
  logical :: split
  split = .false.
end function splitFlash

subroutine assign_from_flash_silcc(thisOctal, subcell)
  use octal_mod
  implicit none

  TYPE(OCTAL) :: thisOctal
  integer     :: subcell

  thisOctal%rho(subcell) = 1.0e-33_db
end subroutine assign_from_flash_silcc

#endif

subroutine deallocate_gridFromFlash

  if (allocated(density))     deallocate (density)
  if (allocated(temperature)) deallocate (temperature)
  if (allocated(tdust))       deallocate (tdust)
  if (allocated(ihp))         deallocate (ihp)
  if (allocated(vx))          deallocate (vx)
  if (allocated(vy))          deallocate (vy)
  if (allocated(vz))          deallocate (vz)
  if (allocated(lrefine))     deallocate (lrefine)
  if (allocated(boundBox))    deallocate (boundBox)

end subroutine deallocate_gridFromFlash

end module gridFromFlash
