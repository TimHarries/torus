module vtk_mod
!
!

! written by tjh
! module to write VTK format files - see
! http://public.kitware.com/VTK/pdf/file-formats.pdf

  use kind_mod
  use constants_mod
  use utils_mod
  use amr_mod
  use mpi_amr_mod
  use messages_mod

  implicit none

  public :: writeVtkfile

  private :: writePoints, writeIndices, writeValue

  logical :: writeHeader

contains



  subroutine writePoints(grid, vtkFilename, nPoints)
    type(GRIDTYPE) :: grid
    integer :: lunit = 69
    integer :: nPoints
    character(len=*) :: vtkFilename

    open(lunit, file=vtkFilename, form="formatted", status="old", position="append")
    if (writeHeader) then
       write(69,'(a,i10,a)') "POINTS ",nPoints, " float"
    endif


    call recursiveWritePoints(grid%octreeRoot, lunit, grid)
    close(lunit)

  contains

    recursive subroutine recursiveWritePoints(thisOctal,lunit, grid)
      type(OCTAL), pointer :: thisOctal, child

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE) :: grid
      integer :: lunit
      integer :: subcell, i
      real :: xp, xm, yp, ym, zm, zp, d, r1, r2, phi, dphi
      type(OCTALVECTOR) :: rVec
      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call recursiveWritePoints(child,lunit, grid)
                  exit
               end if
            end do
         else

#ifdef MPI
            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal) .and. grid%splitOverMPI) cycle
#endif
            if (thisOctal%threed) then
               if (.not.thisOctal%cylindrical) then
                  rVec = subcellCentre(thisOctal,subcell)
                  d = thisOctal%subcellSize/2.d0
                  xp = REAL(rVec%x + d)
                  xm = REAL(rVec%x - d)
                  yp = REAL(rVec%y + d)
                  ym = REAL(rVec%y - d)
                  zp = REAL(rVec%z + d)
                  zm = REAL(rVec%z - d)
                  
                  write(lunit,*) xm, ym, zm
                  
                  write(lunit,*) xp, ym, zm
                  
                  write(lunit,*) xm, yp, zm
                  
                  write(lunit,*) xp, yp, zm
                  
                  write(lunit,*) xm, ym, zp
                  
                  write(lunit,*) xp, ym, zp
                  
                  write(lunit,*) xm, yp, zp
                  
                  write(lunit,*) xp, yp, zp
               else
                  rVec = subcellCentre(thisOctal, subcell)
                  d = thisOctal%subcellSize/2.d0
                  zp = REAL(rVec%z + d)
                  zm = REAL(rVec%z - d)
                  r1 = sqrt(rVec%x**2 + rVec%y**2) - d
                  r2 = sqrt(rVec%x**2 + rVec%y**2) + d
                  phi = atan2(rVec%y, rVec%x)
                  dphi = returndPhi(thisOctal)
                  write(lunit,*) r1*cos(phi-dphi), r1*sin(phi-dphi), zm

                  write(lunit,*) r1*cos(phi+dphi), r1*sin(phi+dphi), zm

                  write(lunit,*) r2*cos(phi+dphi), r2*sin(phi+dphi), zm

                  write(lunit,*) r2*cos(phi-dphi), r2*sin(phi-dphi), zm

                  write(lunit,*) r1*cos(phi-dphi), r1*sin(phi-dphi), zp

                  write(lunit,*) r1*cos(phi+dphi), r1*sin(phi+dphi), zp

                  write(lunit,*) r2*cos(phi+dphi), r2*sin(phi+dphi), zp

                  write(lunit,*) r2*cos(phi-dphi), r2*sin(phi-dphi), zp
               endif
            else
               rVec = subcellCentre(thisOctal,subcell)
               d = thisOctal%subcellSize/2.d0
               xp = REAL(rVec%x + d)
               xm = REAL(rVec%x - d)
               zp = REAL(rVec%z + d)
               zm = REAL(rVec%z - d)

               write(lunit,*) xm, zm, 0.

               write(lunit,*) xp, zm, 0.

               write(lunit,*) xm, zp, 0.

               write(lunit,*) xp, zp, 0.
            endif


         endif
      enddo


    end subroutine recursiveWritePoints
  end subroutine writePoints

  subroutine writeIndices(grid, vtkFilename, nPoints, nCells, iOffset)
    type(GRIDTYPE) :: grid
    integer :: lunit = 69
    integer :: nCells, nPoints, iOffset, i, nCount
    character(len=*) :: vtkFilename

    open(lunit, file=vtkFilename, form="formatted", status="old", position="append")
    if (writeHeader) then
       write(lunit, '(a, i10, i10)') "CELLS ",nCells, nPoints+nCells
    endif

    nCount = 0
    call recursiveWriteIndices(grid%octreeRoot, lunit, nCount, iOffset, grid)
    close(lunit)

  contains

    recursive subroutine recursiveWriteIndices(thisOctal,lunit, nCount, iOffset, grid)
      type(OCTAL), pointer :: thisOctal, child

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE) :: grid
      integer :: lunit
      integer :: subcell, i
      integer :: iOffset, nCount

      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call recursiveWriteIndices(child,lunit, nCount, iOffset, grid)
                  exit
               end if
            end do
         else

#ifdef MPI
            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal) .and. grid%splitOverMPI) cycle
#endif

            if (thisOctal%threed) then
               write(lunit, '(9i10)') 8, nCount + iOffset,&
                    nCount + iOffset + 1, &
                    nCount + iOffset + 2, &
                    nCount + iOffset + 3, &
                    nCount + iOffset + 4, &
                    nCount + iOffset + 5, &
                    nCount + iOffset + 6, &
                    nCount + iOffset + 7
               nCount = nCount + 8
            else

               write(lunit, '(5i10)') 4, nCount + iOffset,&
                    nCount + iOffset + 1, &
                    nCount + iOffset + 2, &
                    nCount + iOffset + 3
               nCount = nCount + 4


            endif
         endif
      enddo


    end subroutine recursiveWriteIndices
  end subroutine writeIndices


  subroutine writeValue(grid, vtkFilename, valueType, nCells)
    type(GRIDTYPE) :: grid
    integer :: lunit = 69
    integer :: nCells
    character(len=*) :: valueType
    character(len=*) :: vtkFilename
    logical :: vector, scalar


    select case (valueType)
       case("velocity","hydrovelocity")
          scalar = .false.
          vector = .true.
       case DEFAULT
          scalar = .true.
          vector = .false.
    end select

    open(lunit, file=vtkFilename, form="formatted", status="old", position="append")
    if (writeHeader) then
       if (scalar) then
          write(69,'(a,a,a)') "SCALARS ",trim(valueType)," float"
          write(69, '(a)') "LOOKUP_TABLE default"
       endif
       if (vector) then
          write(69, '(a,a,a)') "VECTORS ", trim(valueType), " float"
       endif

    endif

    call recursiveWriteValue(grid%octreeRoot, valueType, grid)
    close(lunit)

  contains

    recursive subroutine recursiveWriteValue(thisOctal, valueType, grid)
      type(OCTAL), pointer :: thisOctal, child

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE), intent(in) :: grid
      integer :: lunit = 69
      integer :: subcell, i
      character(len=*) :: valueType

      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call recursiveWriteValue(child, valueType, grid)
                  exit
               end if
            end do
         else

#ifdef MPI
            if ( (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) .and. grid%splitOverMPI ) cycle
#endif

            select case (valueType)
               case("rho")
                  write(lunit, *) real(thisOctal%rho(subcell))
               case("hydrovelocity")
                  write(lunit, *) real(thisOctal%rhou(subcell)/thisOctal%rho(subcell)), &
                       real(thisOctal%rhov(subcell)/thisOctal%rho(subcell)), &
                       real(thisOctal%rhow(subcell)/thisOctal%rho(subcell))
               case("velocity")
                  if (thisOctal%threed) then
                     write(lunit, *) thisOctal%velocity(subcell)%x*cspeed/1.e5, &
                          thisOctal%velocity(subcell)%y*cspeed/1.e5, thisOctal%velocity(subcell)%z*cspeed/1.e5
                  else
                     write(lunit, *) thisOctal%velocity(subcell)%x*cspeed/1.e5, &
                           thisOctal%velocity(subcell)%z*cspeed/1.e5, 0.
                     endif
               case("temperature")
                  write(lunit, *) real(thisOctal%temperature(subcell))
               case("phi")
                  write(lunit, *) real(thisOctal%phi_i(subcell))
               case DEFAULT
                  write(*,*) "Cannot write vtk type ",trim(valueType)
             end select


         endif
      enddo


    end subroutine recursiveWriteValue
  end subroutine writeValue

  subroutine writeVtkFile(grid, vtkFilename, valueTypeFilename)
    use input_variables, only : cylindrical
#ifdef MPI
    include 'mpif.h'
#endif
    type(GRIDTYPE) :: grid
    character(len=*) :: vtkFilename
    integer :: nValueType
    character(len=20) :: valueType(50)
    character(len=*), optional ::  valueTypeFilename
    integer :: nCells, nPoints
    integer :: lunit = 69
    integer :: nOctals, nVoxels, ierr, iOffset, i, iType
    integer :: nPointOffset
#ifdef MPI
    integer :: nPointsArray(64)
    integer, allocatable :: iOffsetArray(:)
    integer :: myRank, nThreads, iThread
#endif


#ifdef MPI
! just return if the grid is not decomposed and MPI job and not zero rank thread
    if ((.not.grid%splitOverMpi).and.(myRankGlobal /= 0)) goto 666
#endif

#ifdef MPI
! just return if the grid is decomposed and MPI job and this is rank zero thread
    if (grid%splitOverMpi.and.(myRankGlobal == 0)) goto 666
#endif

    

    if (PRESENT(valueTypeFilename)) then
       nValueType = 1
       open(29, file=valueTypeFilename, status="old", form="formatted")
10     continue
       read(29,'(a)',end=20) valueType(nValueType)
       nValueType = nValueType + 1
       goto 10
20     continue
       nValueType = nValueType - 1
       close(29)
    else
       nValueType = 2
       valueType(1) = "rho"
       valueType(2) = "velocity"
    endif


    if (grid%octreeRoot%threed) then
       nPointOffset = 8
    else if (grid%octreeRoot%twoD) then
       nPointOffset = 4
    else if (grid%octreeRoot%oneD) then
       nPointOffset = 2
    endif

    call countVoxels(grid%octreeRoot, nOctals, nVoxels)  


#ifdef MPI
    if (grid%splitOverMpi) then
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)
       allocate(iOffsetArray(1:nThreads-1))
       call countSubcellsMPI(grid, nVoxels, nSubcellArray = iOffsetArray)
       iOffsetArray(2:nThreads-1) = iOffsetArray(1:nThreads-2)
       iOffsetArray(1) = 0
       do i = 2, nThreads-1
          iOffsetArray(i) = iOffsetArray(i) + iOffsetArray(i-1)
       enddo
       iOffsetArray = iOffsetArray * nPointOffset 
    endif
#endif


    nCells = nVoxels
    nPoints = nCells * nPointOffset
    writeHeader = .true.
#ifdef MPI
    if (myRankGlobal /= 1 .and. grid%splitOverMpi) then
       writeHeader = .false.
    endif
#endif
    
    if (writeHeader) then
       open(lunit,file=vtkFilename, form="formatted", status="unknown")
       write(lunit,'(a)') "# vtk DataFile Version 2.0"
       write(lunit,'(a,a)') "TORUS AMR data"
       write(lunit,'(a)') "ASCII"
       write(lunit,'(a)') "DATASET UNSTRUCTURED_GRID"
       close(lunit)
    endif

    if (.not.grid%splitOverMPI) call writePoints(grid, vtkFilename, nPoints)

#ifdef MPI
    if (grid%splitOverMpi) then
    do iThread = 1, nThreadsGlobal
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       if (iThread == myRankGlobal) then
          call writePoints(grid, vtkFilename, nPoints)
       endif
    enddo
 endif
#endif

    if (.not.grid%splitOverMPI) call writeIndices(grid, vtkFilename, nPoints, nCells, 0)

#ifdef MPI
    if (grid%splitOverMpi) then

    do iThread = 1, nThreadsGlobal
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       if (iThread == myRankGlobal) then
           call writeIndices(grid, vtkFilename, nPoints, nCells, iOffsetArray(myRankGlobal))
       endif
    enddo
 endif
#endif

    if (writeHeader) then
       open(lunit,file=vtkFilename, form="formatted", status="old",position="append")
       write(lunit, '(a, i10)') "CELL_TYPES ",nCells
       do i = 1, nCells
          if ((nPointOffset == 8).and.(.not.cylindrical)) write(lunit, '(a)') "11"
          if ((nPointOffset == 8).and.(cylindrical)) write(lunit, '(a)') "12"
          if (nPointOffset == 4) write(lunit, '(a)') "8"
       enddo
       write(69, '(a,  i10)') "CELL_DATA ",nCells
       close(lunit)
    endif

    do iType = 1, nValueType
    
       if (.not.grid%splitOverMPI) call writeValue(grid, vtkFilename, valueType(iType), nCells)

#ifdef MPI
    if (grid%splitOverMpi) then

       do iThread = 1, nThreadsGlobal
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          if (iThread == myRankGlobal) then
             call writeValue(grid, vtkFilename, valueType(iType), nCells)
          endif
       enddo
    endif
#endif

    enddo

666 continue

  end subroutine writeVtkFile




end module vtk_mod

