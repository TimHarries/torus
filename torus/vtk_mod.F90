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
  use messages_mod

  implicit none


contains

  subroutine writeVtkFile(grid, vtkFilename, dataType)
    type(GRIDTYPE) :: grid
    character(len=*) :: vtkFilename, dataType
    real, allocatable :: pointx(:), pointy(:), pointz(:), value(:)
    integer :: nOctals, nVoxels, nCell, nPoint, i
    integer, allocatable :: cellPointIndex(:,:)

    open(69,file=vtkFilename, form="formatted", status="unknown")

    write(69,'(a)') "# vtk DataFile Version 2.0"
    write(69,'(a,a)') "TORUS AMR data: ",trim(dataType)
    write(69,'(a)') "ASCII"
    write(69,'(a)') "DATASET UNSTRUCTURED_GRID"

    call countVoxels(grid%octreeRoot, nOctals, nVoxels)  

    nCell = nVoxels
    nPoint = nCell * 8
    allocate(pointx(1:nPoint), pointy(1:nPoint), pointz(1:nPoint), value(1:nCell))
    allocate(cellPointIndex(nCell,8))

    ncell = 0 
    nPoint = 0 
    call getPointandCellData(grid%octreeRoot, pointx, pointy, pointz, npoint, cellPointIndex, value, nCell)
    
    write(69,'(a,i10,a)') "POINTS ",nPoint, " float"
    do i = 1, nPoint
       write(69, *) pointx(i), pointy(i), pointz(i)
    enddo

    write(69, '(a, i10, i10)') "CELLS ",nCell, nPoint+nCell
    do i = 1, nCell
       write(69, '(a,8i10)') " 8 ", cellPointIndex(i,1:8)-1
    enddo

    write(69, '(a, i10)') "CELL_TYPES ",nCell
    do i = 1, nCell
       write(69, '(a)') "11"
    enddo

    write(69, '(a,  i10)') "CELL_DATA ",nCell
    write(69,'(a)') "SCALARS rho float"
    write(69, '(a)') "LOOKUP_TABLE default"
    do i = 1, nCell
       write(69, *) value(i)
    enddo
    close(69)
  end subroutine writeVtkFile


  recursive subroutine getPointandCellData(thisOctal, pointx, pointy, pointz, npoint, cellPointIndex, value, nCell)
    type(OCTAL), pointer :: thisOctal, child
    real :: pointx(:), pointy(:), pointz(:), value(:)
    integer :: nPoint, cellPointIndex(:,:), ncell
    integer :: subcell, i
    type(OCTALVECTOR) :: rVec
    real :: d, xp, xm, yp, ym, zp, zm

    do subcell = 1, thisOctal%maxChildren
       if (thisOctal%hasChild(subcell)) then
          ! find the child
          do i = 1, thisOctal%nChildren, 1
             if (thisOctal%indexChild(i) == subcell) then
                child => thisOctal%child(i)
                call  getPointandCellData(child, pointx, pointy, pointz, npoint, cellPointIndex, value, nCell)
                exit
             end if
          end do
       else
          nCell = nCell + 1
          rVec = subcellCentre(thisOctal,subcell)
          d = thisOctal%subcellSize/2.d0
          xp = REAL(rVec%x + d)
          xm = REAL(rVec%x - d)
          yp = REAL(rVec%y + d)
          ym = REAL(rVec%y - d)
          zp = REAL(rVec%z + d)
          zm = REAL(rVec%z - d)
          
          value(nCell) = thisOctal%rho(subcell)
          nPoint = nPoint + 1
          pointx(nPoint) = xm
          pointy(nPoint) = ym
          pointz(nPoint) = zm
          cellPointIndex(nCell, 1) = nPoint

          nPoint = nPoint + 1
          pointx(nPoint) = xp
          pointy(nPoint) = ym
          pointz(nPoint) = zm
          cellPointIndex(nCell, 2) = nPoint

          nPoint = nPoint + 1
          pointx(nPoint) = xm
          pointy(nPoint) = yp
          pointz(nPoint) = zm
          cellPointIndex(nCell, 3) = nPoint


          nPoint = nPoint + 1
          pointx(nPoint) = xp
          pointy(nPoint) = yp
          pointz(nPoint) = zm
          cellPointIndex(nCell, 4) = nPoint

          nPoint = nPoint + 1
          pointx(nPoint) = xm
          pointy(nPoint) = ym
          pointz(nPoint) = zp
          cellPointIndex(nCell, 5) = nPoint

          nPoint = nPoint + 1
          pointx(nPoint) = xp
          pointy(nPoint) = ym
          pointz(nPoint) = zp
          cellPointIndex(nCell, 6) = nPoint

          nPoint = nPoint + 1
          pointx(nPoint) = xm
          pointy(nPoint) = yp
          pointz(nPoint) = zp
          cellPointIndex(nCell, 7) = nPoint

          nPoint = nPoint + 1
          pointx(nPoint) = xp
          pointy(nPoint) = yp
          pointz(nPoint) = zp
          cellPointIndex(nCell, 8) = nPoint

       end if
    enddo

  end subroutine getPointandCellData


end module vtk_mod

