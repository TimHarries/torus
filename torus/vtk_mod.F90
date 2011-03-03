module vtk_mod
!
!

! written by tjh
! module to write VTK format files - see
! http://public.kitware.com/VTK/pdf/file-formats.pdf

  use kind_mod
  use ion_mod
  use dimensionality_mod
  use constants_mod
  use utils_mod
  use amr_mod
  use mpi_amr_mod
  use messages_mod
  use vector_mod
  use mpi_global_mod

  implicit none

  interface writeVtkFile
     module procedure writeVtkFileAMR
     module procedure writeVtkFileSource
  end interface

  public :: writeVtkfile !, writeIfritfile

  private :: writePoints, writeIndices, writeValue

  logical :: writeHeader

contains



  subroutine writePoints(grid, vtkFilename, nPoints, ncells,  xml)
    type(GRIDTYPE) :: grid
    integer :: lunit = 69
    integer :: nPoints, nCells
    logical :: xml
    character(len=*) :: vtkFilename

    open(lunit, file=vtkFilename, form="formatted", status="old", position="append")
    if (writeHeader) then
       if (.not.xml) then
          write(69,'(a,i12,a)') "POINTS ",nPoints, " float"
       else
          write(69,'(a,i12,a,i12,a)') "<Piece NumberOfPoints='",nPoints,"' NumberOfCells='",nCells,"'>"
          write(69,'(a)') "<Points>"
          write(69,'(a)') "<DataArray type='Float32' Name='Position' NumberOfComponents='3' format='ascii'>"
       endif
    endif


    call recursiveWritePoints(grid%octreeRoot, lunit, grid)
    if (writeheader.and.xml) then
       write(69,'(a)') "</DataArray>"
       write(69,'(a)') "</Points>"
    endif
    close(lunit)
       

  contains

    recursive subroutine recursiveWritePoints(thisOctal,lunit, grid)
      use octal_mod, only: returndPhi

      type(OCTAL), pointer :: thisOctal, child

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE) :: grid
      integer :: lunit
      integer :: subcell, i
      real :: xp, xm, yp, ym, zm, zp, d, r1, r2, phi, dphi
      real :: phi1, phi2, phiStart, phiEnd
      integer :: nPhi, iPhi
      type(VECTOR) :: rVec
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
                  phi1 = phi - dphi
                  phi2 = phi + dphi
                  nPhi = max(1, nint(dphi/(10.d0*degtorad)))
                  do iPhi = 1, nPhi
                     phiStart = phi1 + (phi2 - phi1) * dble(iphi-1)/dble(nphi)
                     phiEnd   = phi1 + (phi2 - phi1) * dble(iphi)/dble(nphi)

                     write(lunit,*) r1*cos(phiStart), r1*sin(phiStart), zm
                     
                     write(lunit,*) r1*cos(phiEnd), r1*sin(phiEnd), zm
                     
                     write(lunit,*) r2*cos(phiEnd), r2*sin(phiEnd), zm
                     
                     write(lunit,*) r2*cos(phiStart), r2*sin(phiStart), zm
                     
                     write(lunit,*) r1*cos(phiStart), r1*sin(phiStart), zp
                     
                     write(lunit,*) r1*cos(phiEnd), r1*sin(phiEnd), zp
                     
                     write(lunit,*) r2*cos(phiEnd), r2*sin(phiEnd), zp
                     
                     write(lunit,*) r2*cos(phiStart), r2*sin(phiStart), zp
                  enddo
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

  subroutine getPoints(grid, nPoints, points)
    type(GRIDTYPE) :: grid
    real :: points(:,:)
    integer :: nPoints

    nPoints = 0

    call recursiveGetPoints(grid, grid%octreeRoot, nPoints, Points)

  contains

    recursive subroutine recursiveGetPoints(grid, thisOctal,nPoints, points)
      use octal_mod, only: returndPhi

      type(OCTAL), pointer :: thisOctal, child

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE) :: grid
      real :: points(:,:)
      integer :: nPoints
      integer :: subcell, i
      real :: xp, xm, yp, ym, zm, zp, d, r1, r2, phi, dphi
      real :: phi1, phi2, phiStart, phiEnd
      integer :: nPhi, iPhi
      type(VECTOR) :: rVec
      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call recursiveGetPoints(grid, child, nPoints, points)
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
                  
                  nPoints = nPoints + 1
                  points(1, nPoints) = xm
                  points(2, nPoints) = ym
                  points(3, nPoints) = zm

                  nPoints = nPoints + 1
                  points(1, nPoints) = xp
                  points(2, nPoints) = ym
                  points(3, nPoints) = zm


                  nPoints = nPoints + 1
                  points(1, nPoints) = xm
                  points(2, nPoints) = yp
                  points(3, nPoints) = zm

                  nPoints = nPoints + 1
                  points(1, nPoints) = xp
                  points(2, nPoints) = yp
                  points(3, nPoints) = zm

                  nPoints = nPoints + 1
                  points(1, nPoints) = xm
                  points(2, nPoints) = ym
                  points(3, nPoints) = zp

                  nPoints = nPoints + 1
                  points(1, nPoints) = xp
                  points(2, nPoints) = ym
                  points(3, nPoints) = zp

                  nPoints = nPoints + 1
                  points(1, nPoints) = xm
                  points(2, nPoints) = yp
                  points(3, nPoints) = zp

                  nPoints = nPoints + 1
                  points(1, nPoints) = xp
                  points(2, nPoints) = yp
                  points(3, nPoints) = zp

               else
                  rVec = subcellCentre(thisOctal, subcell)
                  d = thisOctal%subcellSize/2.d0
                  zp = REAL(rVec%z + d)
                  zm = REAL(rVec%z - d)
                  r1 = sqrt(rVec%x**2 + rVec%y**2) - d
                  r2 = sqrt(rVec%x**2 + rVec%y**2) + d
                  phi = atan2(rVec%y, rVec%x)
                  dphi = returndPhi(thisOctal)
                  phi1 = phi - dphi
                  phi2 = phi + dphi
                  nPhi = max(1, nint(dphi/(10.d0*degtorad)))
                  do iPhi = 1, nPhi
                     phiStart = phi1 + (phi2 - phi1) * dble(iphi-1)/dble(nphi)
                     phiEnd   = phi1 + (phi2 - phi1) * dble(iphi)/dble(nphi)


                     nPoints = nPoints + 1
                     points(1, nPoints) = r1*cos(phiStart)
                     points(2, nPoints) = r1*sin(phiStart)
                     points(3, nPoints) = zm


                     nPoints = nPoints + 1
                     points(1, nPoints) = r1*cos(phiEnd)
                     points(2, nPoints) = r1*sin(phiEnd)
                     points(3, nPoints) = zm

                     nPoints = nPoints + 1
                     points(1, nPoints) = r2*cos(phiEnd)
                     points(2, nPoints) = r2*sin(phiEnd)
                     points(3, nPoints) = zm

                     nPoints = nPoints + 1
                     points(1, nPoints) = r2*cos(phiStart)
                     points(2, nPoints) = r2*sin(phiStart)
                     points(3, nPoints) = zm

                     nPoints = nPoints + 1
                     points(1, nPoints) = r1*cos(phiStart)
                     points(2, nPoints) = r1*sin(phiStart)
                     points(3, nPoints) = zp

                     nPoints = nPoints + 1
                     points(1, nPoints) = r1*cos(phiEnd)
                     points(2, nPoints) = r1*sin(phiEnd)
                     points(3, nPoints) = zp


                     nPoints = nPoints + 1
                     points(1, nPoints) = r2*cos(phiEnd)
                     points(2, nPoints) = r2*sin(phiEnd)
                     points(3, nPoints) = zp

                     nPoints = nPoints + 1
                     points(1, nPoints) = r2*cos(phiStart)
                     points(2, nPoints) = r2*sin(phiStart)
                     points(3, nPoints) = zp
                  enddo
               endif
            else
               rVec = subcellCentre(thisOctal,subcell)
               d = thisOctal%subcellSize/2.d0
               xp = REAL(rVec%x + d)
               xm = REAL(rVec%x - d)
               zp = REAL(rVec%z + d)
               zm = REAL(rVec%z - d)

               nPoints = nPoints + 1
               points(1, nPoints) = xm
               points(2, nPoints) = zm
               points(3, nPoints) = 0.

               nPoints = nPoints + 1
               points(1, nPoints) = xp
               points(2, nPoints) = zm
               points(3, nPoints) = 0.

               nPoints = nPoints + 1
               points(1, nPoints) = xm
               points(2, nPoints) = zp
               points(3, nPoints) = 0.


               nPoints = nPoints + 1
               points(1, nPoints) = xp
               points(2, nPoints) = zp
               points(3, nPoints) = 0.

            endif

         endif
      enddo


    end subroutine recursiveGetPoints
  end subroutine getPoints

  subroutine writeIndices(grid, vtkFilename, nPoints, nCells, iOffset, xml)
    type(GRIDTYPE) :: grid
    integer :: lunit = 69
    integer :: nCells, nPoints, iOffset, nCount
    character(len=*) :: vtkFilename
    logical :: xml

    open(lunit, file=vtkFilename, form="formatted", status="old", position="append")
    if (writeHeader) then
       if (.not.xml) then
          write(lunit, '(a, i10, i10)') "CELLS ",nCells, nPoints+nCells
       else
          write(lunit,'(a)') "<Cells>"
          write(lunit,'(a)') "<DataArray type='Int32' Name='connectivity' format='ascii' >"
       endif
    endif

    nCount = 0
    call recursiveWriteIndices(grid%octreeRoot, lunit, nCount, iOffset, grid, xml)
    if (writeHeader.and.xml) then
       write(lunit,'(a)') "</DataArray>"
    endif
    close(lunit)

  contains

    recursive subroutine recursiveWriteIndices(thisOctal,lunit, nCount, iOffset, grid, xml)
      use input_variables, only : cylindrical
      type(OCTAL), pointer :: thisOctal, child
      logical :: xml

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE) :: grid
      integer :: lunit
      integer :: subcell, i,j,n
      integer :: iOffset, nCount

      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call recursiveWriteIndices(child,lunit, nCount, iOffset, grid, xml)
                  exit
               end if
            end do
         else

#ifdef MPI
            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal) .and. grid%splitOverMPI) cycle
#endif

            if (.not.xml) then

               if (thisOctal%threed) then
                  if (.not.cylindrical) then
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
                     n = max(1, nint(returnDphi(thisOctal)/(10.d0*degtorad)))
                     do j = 1 , n
                        write(lunit, '(9i10)') 8, nCount + iOffset,&
                             nCount + iOffset + 1, &
                             nCount + iOffset + 2, &
                             nCount + iOffset + 3, &
                             nCount + iOffset + 4, &
                             nCount + iOffset + 5, &
                             nCount + iOffset + 6, &
                             nCount + iOffset + 7
                        nCount = nCount + 8
                     enddo
                  endif
               else
                  
                  write(lunit, '(5i10)') 4, nCount + iOffset,&
                       nCount + iOffset + 1, &
                       nCount + iOffset + 2, &
                       nCount + iOffset + 3
                  nCount = nCount + 4
                  

               endif
            else

               if (thisOctal%threed) then
                  if (.not.cylindrical) then
                     write(lunit, '(8i10)') nCount + iOffset,&
                          nCount + iOffset + 1, &
                          nCount + iOffset + 2, &
                          nCount + iOffset + 3, &
                          nCount + iOffset + 4, &
                          nCount + iOffset + 5, &
                          nCount + iOffset + 6, &
                          nCount + iOffset + 7
                     nCount = nCount + 8
                  else
                     n = max(1, nint(returnDphi(thisOctal)/(10.d0*degtorad)))
                     do j = 1 , n
                        write(lunit, '(8i10)') nCount + iOffset,&
                             nCount + iOffset + 1, &
                             nCount + iOffset + 2, &
                             nCount + iOffset + 3, &
                             nCount + iOffset + 4, &
                             nCount + iOffset + 5, &
                             nCount + iOffset + 6, &
                             nCount + iOffset + 7
                        nCount = nCount + 8
                     enddo
                  endif
               else
                  
                  write(lunit, '(4i10)') nCount + iOffset,&
                       nCount + iOffset + 1, &
                       nCount + iOffset + 2, &
                       nCount + iOffset + 3
                  nCount = nCount + 4
                  

               endif


            endif
         endif
      enddo


    end subroutine recursiveWriteIndices
  end subroutine writeIndices


  subroutine writeOffsetsXML(grid, vtkFilename, iOffset)
    type(GRIDTYPE) :: grid
    integer :: lunit = 69
    integer :: iOffset, nCount
    character(len=*) :: vtkFilename

    open(lunit, file=vtkFilename, form="formatted", status="old", position="append")
    if (writeHeader) then
       write(lunit,'(a)') "<DataArray type='Int32' Name='offsets' format='ascii' >"
    endif

    nCount = 0
    call recursiveWriteOffsetsXML(grid%octreeRoot, lunit, nCount, iOffset, grid)
    if (writeHeader) then
       write(lunit,'(a)') "</DataArray>"
    endif
    close(lunit)

  contains

    recursive subroutine recursivewriteOffsetsXML(thisOctal,lunit, nCount, iOffset, grid)
      use input_variables, only : cylindrical
      type(OCTAL), pointer :: thisOctal, child

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE) :: grid
      integer :: lunit
      integer :: subcell, i,j,n
      integer :: iOffset, nCount

      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call recursiveWriteOffsetsXML(child,lunit, nCount, iOffset, grid)
                  exit
               end if
            end do
         else

#ifdef MPI
            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal) .and. grid%splitOverMPI) cycle
#endif


            if (thisOctal%threed) then
               if (.not.cylindrical) then
                  nCount = nCount + 8
                  write(lunit, '(i10)') nCount
               else
                  n = max(1, nint(returnDphi(thisOctal)/(10.d0*degtorad)))
                  do j = 1 , n
                     nCount = nCount + 8
                     write(lunit, '(i10)') nCount
                  enddo
               endif
            else

               nCount = nCount + 4
               write(lunit, '(i10)') nCount

            endif
         endif
      enddo

    end subroutine recursivewriteOffsetsXML
  end subroutine writeOffsetsXML

  subroutine getOffsets(grid, offsets, iOffset)
    type(GRIDTYPE) :: grid
    integer ::   iOffset, nCount
    integer :: offsets(:)


    ioffSet = 0 
    nCount = 0
    call getOffsetsRecursive(grid%octreeRoot, nCount, offsets, iOffset, grid)

  contains

    recursive subroutine getOffsetsRecursive(thisOctal, nCount, offsets, iOffset, grid)
      use input_variables, only : cylindrical
      type(OCTAL), pointer :: thisOctal, child

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE) :: grid
      integer :: subcell, i,j,n
      integer :: iOffset, nCount
      integer :: offsets(:)

      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call getOffsetsRecursive(child, nCount, offsets, iOffset, grid)
                  exit
               end if
            end do
         else

#ifdef MPI
            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal) .and. grid%splitOverMPI) cycle
#endif


            if (thisOctal%threed) then
               if (.not.cylindrical) then
                  nCount = nCount + 8
                  iOffset = iOffset + 1
                  offsets(ioffset) =  nCount
               else
                  n = max(1, nint(returnDphi(thisOctal)/(10.d0*degtorad)))
                  do j = 1 , n
                     nCount = nCount + 8
                     iOffset = iOffset + 1
                     offsets(ioffset) =  nCount
                  enddo
               endif
            else

               nCount = nCount + 4
               iOffset = iOffset + 1
               offsets(ioffset) =  nCount

            endif
         endif
      enddo

    end subroutine getOffsetsRecursive
  end subroutine getOffsets


  subroutine writeValue(grid, vtkFilename, valueType, xml)
    type(GRIDTYPE) :: grid
    integer :: lunit = 69
    character(len=*) :: valueType
    character(len=*) :: vtkFilename
    logical :: vectorvalue, scalarvalue
    logical :: xml


    select case (valueType)
       case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom")
          scalarvalue = .false.
          vectorvalue = .true.
       case DEFAULT
          scalarvalue = .true.
          vectorvalue = .false.
    end select

    open(lunit, file=vtkFilename, form="formatted", status="old", position="append")
    if (writeHeader) then
       if (.not.xml) then
          if (scalarvalue) then
             write(69,'(a,a,a)') "SCALARS ",trim(valueType)," float"
             write(69, '(a)') "LOOKUP_TABLE default"
          endif
          if (vectorvalue) then
             write(69, '(a,a,a)') "VECTORS ", trim(valueType), " float"
          endif
       else
          if (scalarvalue) then
             write(69,'(a,a,a)') "<DataArray type='Float32' Name='", &
                  trim(valuetype),"' format='ascii' >"
          endif
          if (vectorvalue) then
             write(69, '(a,a,a)') "<DataArray type='Float32' NumberOfComponents='3' Name='", &
                  trim(valueType), "' format='ascii' >"
          endif
             
       endif
    endif

    call recursiveWriteValue(grid%octreeRoot, valueType, grid, xml)

    if (writeheader.and.xml) then
       write(69,'(a)') "</DataArray>"
    endif

    close(lunit)

  contains

    recursive subroutine recursiveWriteValue(thisOctal, valueType, grid, xml)
      use input_variables, only : lambdasmooth, cylindrical
      type(OCTAL), pointer :: thisOctal, child
      logical :: xml
#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE), intent(in) :: grid
      type(VECTOR) :: rVec, vel
      integer :: lunit = 69
      integer :: subcell, i, j, nVal
      integer, save :: iLambda
      real :: value
      real(double) :: kAbs, kSca
      character(len=*) :: valueType
      real, parameter :: min_single_prec = 1.0e-37
      logical, save :: firstTime = .true.

      kabs = 0.d0
      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call recursiveWriteValue(child, valueType, grid, xml)
                  exit
               end if
            end do
         else

#ifdef MPI
            if ( (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) .and. grid%splitOverMPI ) cycle
#endif

            nVal = 1
            if (thisOctal%threed.and.cylindrical) then
               nVal = max(1, nint(returndPhi(thisOctal)/(10.d0*degtorad)))
            endif


            do j = 1, nVal
            select case (valueType)
               case("rho")

!                  if(thisOctal%rho(subcell) .lt. 1.d-37) then
!                     write(lunit, *) 1e-37 ! floating underflow warning
!                  else
                     write(lunit, *) thisOctal%rho(subcell)
!                  endif

! used for accuracy testing in molecular_mod
               case("interprho")
!                  write(lunit, *) thisOctal%interprho(subcell)

               case("J=0")
!                  if(thisOctal%molecularlevel(subcell,1) .lt. 1.d-37) then
!                     write(lunit, *) 1e-37
!                  else
                     write(lunit, *) thisOctal%molecularlevel(1,subcell)
!                  endif

               case("J=1")
!                  if(thisOctal%molecularlevel(subcell,2) .lt. 1.d-37) then
!                     write(lunit, *) 1e-37
!                  else
                     write(lunit, *) thisOctal%molecularlevel(2,subcell)
!                  endif

               case("J=2")
                  write(lunit, *) thisOctal%molecularlevel(3,subcell)

               case("J=3")
                  write(lunit, *) thisOctal%molecularlevel(4,subcell)

               case("J=4")
                  write(lunit, *) thisOctal%molecularlevel(5,subcell)

               case("J=5")
                  write(lunit, *) thisOctal%molecularlevel(6,subcell)

               case("J=10")
                  write(lunit, *) thisOctal%molecularlevel(11,subcell)

               case("J=16")
                  write(lunit, *) thisOctal%molecularlevel(17,subcell)

               case("dI")
                  write(lunit, *) thisOctal%newmolecularlevel(1,subcell)

               case("dIattenuated")
                  write(lunit, *) thisOctal%newmolecularlevel(2,subcell)

               case("i0")
                  write(lunit, *) thisOctal%newmolecularlevel(3,subcell)

! galLon and galLat re-use storage used for i0 and dIattenuated 
               case("galLon")
                  write(lunit, *) thisOctal%newmolecularlevel(2,subcell)

               case("galLat")
                  write(lunit, *) thisOctal%newmolecularlevel(3,subcell)

               case("crossing")
                  write(lunit, *) thisOctal%newmolecularlevel(4,subcell)

               case("tauacrosscell")
                  write(lunit, *) thisOctal%newmolecularlevel(5,subcell)

               case("tau10")
                  write(lunit, *) thisOctal%tau(1,subcell)

               case("tau21")
                  write(lunit, *) thisOctal%tau(2,subcell)

               case("tau32")
                  write(lunit, *) thisOctal%tau(3,subcell)

               case("tau43")
                  write(lunit, *) thisOctal%tau(4,subcell)

               case("tau54")
                  write(lunit, *) thisOctal%tau(5,subcell)

               case("tau65")
                  write(lunit, *) thisOctal%tau(6,subcell)

               case("tau76")
                  write(lunit, *) thisOctal%tau(7,subcell)

               case("tau87")
                  write(lunit, *) thisOctal%tau(8,subcell)

               case("level0error")
                  write(lunit, *) 10**((real(thisOctal%levelconvergence(1,subcell)) / 6553.6) - 4)

               case("level1error")
                  write(lunit, *) 10**((real(thisOctal%levelconvergence(2,subcell)) / 6553.6) - 4)

               case("level2error")
                  write(lunit, *) 10**((real(thisOctal%levelconvergence(3,subcell)) / 6553.6) - 4)

               case("level3error")
                  write(lunit, *) 10**((real(thisOctal%levelconvergence(4,subcell)) / 6553.6) - 4)

               case("level4error")
                  write(lunit, *) 10**((real(thisOctal%levelconvergence(5,subcell)) / 6553.6) - 4)

               case("niter")
                  write(lunit, *) int(thisOctal%convergence(subcell)/100)

               case("nh2")
                  write(lunit, *) real(thisOctal%nh2(subcell))

               case("convergence")
                  write(lunit, *) mod(thisOctal%convergence(subcell),1.0)

               case("adot")
                  write(lunit, *) real(thisOctal%adot(subcell))

               case("slowestlevel")
                  write(lunit, *) floor(mod(thisoctal%convergence(subcell),100.0))

               case("molabundance")
                  write(lunit, *) thisOctal%molabundance(subcell)

               case("bnu")
                  write(*,*) thisOctal%bnu(1,subcell)
                  write(lunit, *) real(thisOctal%bnu(1,subcell))

               case("dc0")
                  write(lunit, *) thisOctal%molecularlevel(1,subcell) * thisOctal%departcoeff(1,subcell)

               case("dc1")
                  write(lunit, *) thisOctal%molecularlevel(2,subcell) * thisOctal%departcoeff(2,subcell)

               case("dc2")
                  write(lunit, *) thisOctal%molecularlevel(3,subcell) * thisOctal%departcoeff(3,subcell)

               case("dc3")
                  write(lunit, *) thisOctal%molecularlevel(4,subcell) * thisOctal%departcoeff(4,subcell)

               case("dc4")
                  write(lunit, *) thisOctal%molecularlevel(5,subcell) * thisOctal%departcoeff(5,subcell)

               case("jnu10")
                  write(lunit, *) thisOctal%jnu(1,subcell)

               case("dust1")
                  write(lunit, *) real(thisOctal%dustTypeFraction(subcell,1))

               case("dust2")
                  write(lunit, *) real(thisOctal%dustTypeFraction(subcell,2))

               case("bias")
                  write(lunit, *) real(thisOctal%biasCont3d(subcell))

               case("mpithread")
                  write(lunit, *) real(thisOctal%mpithread(subcell))

               case("bcond")
                  write(lunit, *) real(thisOctal%boundaryCondition(subcell))

               case("deltaT")
                  write(lunit, *) real(thisOctal%temperature(subcell)-thisOctal%oldtemperature(subcell))


               case("scattered")
                  write(lunit, *) real(thisOctal%scatteredIntensity(subcell,5,3))

               case("hydrovelocity")
                  if (thisOctal%threeD) then
                     write(lunit, *) real(thisOctal%rhou(subcell)/thisOctal%rho(subcell)), &
                          real(thisOctal%rhov(subcell)/thisOctal%rho(subcell)), &
                          real(thisOctal%rhow(subcell)/thisOctal%rho(subcell))
                  else
                     write(lunit, *) real(thisOctal%rhou(subcell)/thisOctal%rho(subcell)), &
                          real(thisOctal%rhow(subcell)/thisOctal%rho(subcell)), &
                          real(thisOctal%rhov(subcell)/thisOctal%rho(subcell))
                  endif

               case("radmom")
                     write(lunit, *) real(thisOctal%radiationMomentum(subcell)%x/1.d20), &
                          real(thisOctal%radiationMomentum(subcell)%y/1.d20), &
                          real(thisOctal%radiationMomentum(subcell)%z/1.d20)
               case("velocity")
                     ! stop vectors from showing up in visit if too big
                     if(thisoctal%velocity(subcell)%x .ge. 1.d0) then
                        thisoctal%velocity(subcell)%x = 0.d0
                        thisoctal%velocity(subcell)%y = 0.d0
                        thisoctal%velocity(subcell)%z = 0.d0
                     endif
                     if (thisOctal%threed) then
                        write(lunit, *) real(thisOctal%velocity(subcell)%x*cspeed/1.e5), &
                             real(thisOctal%velocity(subcell)%y*cspeed/1.e5), &
                             real(thisOctal%velocity(subcell)%z*cspeed/1.e5)
                     else
                        write(lunit, *) real(thisOctal%velocity(subcell)%x*cspeed/1.e5), &
                             real(thisOctal%velocity(subcell)%z*cspeed/1.e5), &
                             real(thisOctal%velocity(subcell)%y*Cspeed/1.e5)
                     endif
              case("cornervel")
                 rVec = subcellCentre(thisOctal, subcell)
                 vel = amrGridVelocity(grid%octreeRoot,rvec,startOctal=thisOctal,&
                      actualSubcell=subcell) 
                 write(lunit, *) real(vel%x*cspeed/1.e5), real(vel%y*cspeed/1.e5), &
                      real(vel%z*cspeed/1.e5)

!               case("quadvelocity")
!                     write(lunit, *) thisOctal%quadvelocity(subcell)%x*cspeed/1.e5, &
!                          thisOctal%quadvelocity(subcell)%y*cspeed/1.e5, thisOctal%quadvelocity(subcell)%z*cspeed/1.e5

!               case("linearvelocity")
!                     write(lunit, *) thisOctal%linearvelocity(subcell)%x*cspeed/1.e5, &
!                          thisOctal%linearvelocity(subcell)%y*cspeed/1.e5, thisOctal%linearvelocity(subcell)%z*cspeed/1.e5


               case("ne")
                  write(lunit, *) real(thisOctal%ne(subcell))

               case("pressure")
                  write(lunit, *) real(thisOctal%pressure_i(subcell))



               case("inflow")
                  if (thisOctal%inflow(subcell)) then
                     write(lunit, *) 1.
                  else
                     write(lunit, *) 0.
                  endif

               case("NH2")
                  write(lunit, *) real( thisOctal%NH2(subcell) )

               case("fixedtemp")
                  if (thisOctal%fixedTemperature(subcell)) then
                     write(lunit, *) 1.
                  else
                     write(lunit, *) 0.
                  endif

               case("HI")
                  write(lunit, *) real(thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon)))

               case("HeI")
                  write(lunit, *) real(thisOctal%ionfrac(subcell,returnIonNumber("He I", grid%ion, grid%nIon)))

               case("HeII")
                  write(lunit, *) real(thisOctal%ionfrac(subcell,returnIonNumber("He II", grid%ion, grid%nIon)))

               case("OI")
                  write(lunit, *) real(thisOctal%ionfrac(subcell,returnIonNumber("O I", grid%ion, grid%nIon)))

               case("OII")
                  write(lunit, *) real(thisOctal%ionfrac(subcell,returnIonNumber("O II", grid%ion, grid%nIon)))

               case("OIII")
                  write(lunit, *) real(thisOctal%ionfrac(subcell,returnIonNumber("O III", grid%ion, grid%nIon)))

               case("temperature")
                  write(lunit, *) real(thisOctal%temperature(subcell))

               case("chiline")
                  write(lunit, *) max( real(thisOctal%chiline(subcell)), min_single_prec )


               case("microturb")
                  write(lunit, *) max( real(thisOctal%microturb(subcell)*cspeed/1.e5), min_single_prec )

               case("etaline")
                  write(lunit, *) max ( real(thisOctal%etaline(subcell)), min_single_prec )

               case("sourceline")
                  write(lunit, *) max ( real(thisOctal%etaline(subcell)), min_single_prec )/ &
                       max( real(thisOctal%chiline(subcell)), min_single_prec )

               case("tau")
                  if (firstTime) then
                     call locate(grid%lamArray, grid%nLambda, lambdaSmooth, ilambda)
                     firstTime = .false.
                  endif
                  call returnKappa(grid, thisOctal, subcell, ilambda=ilambda,&
                       kappaSca=ksca, kappaAbs=kabs)
                  value = thisOctal%subcellSize * (ksca + kabs)
!                  write(*,*) "sca ",thisOctal%subcellSize * ksca,ilambda,lambdasmooth
                  write(lunit, *) real(value)

               case("ross")

                  call returnKappa(grid, thisOctal, subcell, rosselandKappa=kabs)
                  value = thisOctal%subcellsize * kabs * thisOctal%rho(subcell) * 1.e10

                  write(lunit, *) real(value)


               case("dusttype")
                  write(lunit,*) real(thisOctal%dustTypeFraction(subcell,1))

               case("etacont")
                  write(lunit, *) real(thisOctal%etaCont(subcell))

               case("crossings")
                  value = thisOctal%ncrossings(subcell)
                  if (thisOctal%diffusionApprox(subcell)) value = 1.e6
                  write(lunit, *) real(value)
                  
               case("phi")
                  write(lunit, *) real(thisOctal%phi_i(subcell))

               case("q_i")
                  write(lunit, *) real(thisOctal%q_i(subcell))

               case("u_i")
                  write(lunit, *) real(thisOctal%u_interface(subcell))

               case("rhoe")
                  write(lunit, *) real(thisOctal%rhoe(subcell))

               case("edens_s")
                  write(lunit, *) real(thisOctal%photonEnergyDensityFromSource(subcell))

               case("edens_g")
                  write(lunit, *) real(thisOctal%photonEnergyDensityFromGas(subcell))

               case("rhou")
                  write(lunit, *) real(thisOctal%rhou(subcell))

               case("rhov")
                  write(lunit, *) real(thisOctal%rhov(subcell))

               case("rhow")
                  write(lunit, *) real(thisOctal%rhow(subcell))

               case("q_i-1")
                  write(lunit, *) real(thisOctal%q_i_minus_1(subcell))

               case("ghosts")
                  if (thisOctal%ghostCell(subcell)) then
                     write(lunit, *) 1.
                  else
                     write(lunit,*) 0.
                  endif

               case("edges")
                  if (thisOctal%edgeCell(subcell)) then
                     write(lunit, *) 1.
                  else
                     write(lunit,*) 0.
                  endif

               case DEFAULT
                  write(*,*) "Cannot write vtk type ",trim(valueType)
             end select
          enddo


       endif
    enddo


    end subroutine recursiveWriteValue
  end subroutine writeValue

  subroutine writeVTKfileSource(nSource, source, vtkFilename)
    use source_mod
    use mpi_global_mod
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    character(len=*) :: vtkFilename
    integer, parameter :: lunit = 37
    integer :: nPoints, nElements, iSource
    integer :: i
    type(VECTOR) :: cVec, aVec, v1, v2, v3, v4, zAxis
    real(double) :: dphi, dtheta
    integer :: nCount
    

    if (myrankGlobal /=0 ) goto 666
    zAxis = VECTOR(0.d0, 0.d0, 1.d0)

    open(lunit,file=vtkFilename, form="formatted", status="unknown")
    write(lunit,'(a)') "# vtk DataFile Version 2.0"
    write(lunit,'(a,a)') "TORUS AMR data"
    write(lunit,'(a)') "ASCII"
    write(lunit,'(a)') "DATASET UNSTRUCTURED_GRID"

    nPoints = 0 
    do iSource = 1 , nSource
       nPoints = nPoints + source(iSource)%surface%nElements
    enddo
    nElements = nPoints
    nPoints = nPoints * 4
    write(lunit,'(a,i10,a)') "POINTS ",nPoints, " float"

    do iSource = 1, nSource
       do i = 1, source(iSource)%surface%nElements
          cVec = source(iSource)%surface%element(i)%position
          dphi = source(iSource)%surface%element(i)%dphi
          dtheta = source(iSource)%surface%element(i)%dtheta
          aVec = cVec.cross.zAxis
          call normalize(aVec)

          v1 = arbitraryRotate(cVec, -dtheta/2.d0, aVec)
          v1 = rotateZ(v1, dphi/2.d0)

          v1 = v1 + source(isource)%position

          v2 = arbitraryRotate(cVec, -dtheta/2.d0, aVec)
          v2 = rotateZ(v2, -dphi/2.d0)

          v2 = v2 + source(isource)%position

          v3 = arbitraryRotate(cVec, dtheta/2.d0, aVec)
          v3 = rotateZ(v3, dphi/2.d0)

          v3 = v3 + source(isource)%position

          v4 = arbitraryRotate(cVec, dtheta/2.d0, aVec)
          v4 = rotateZ(v4, -dphi/2.d0)

          v4 = v4 + source(isource)%position

          write(lunit,*) v1%x, v1%y, v1%z
          write(lunit,*) v2%x, v2%y, v2%z
          write(lunit,*) v3%x, v3%y, v3%z
          write(lunit,*) v4%x, v4%y, v4%z
       enddo
    enddo

    write(lunit, '(a, i10, i10)') "CELLS ",nElements, nPoints+nElements
    ncount = 0
    do iSource = 1, nSource
       do i = 1, source(iSource)%surface%nElements
          write(lunit, '(5i10)') 4, nCount,&
                    nCount + 1, &
                    nCount + 2, &
                    nCount + 3
               nCount = nCount + 4
       enddo
    enddo
    write(lunit, '(a, i10)') "CELL_TYPES ",nElements
    do i = 1, nElements
       write(lunit, '(a)') "8"
    enddo
    write(lunit, '(a,  i10)') "CELL_DATA ",nElements

    write(lunit,'(a,a,a)') "SCALARS ","temperature"," float"
    write(lunit, '(a)') "LOOKUP_TABLE default"

    do iSource = 1, nSource
       do i = 1, source(iSource)%surface%nElements
          write(lunit,*) source(isource)%surface%element(i)%temperature
       enddo
    enddo

    write(lunit,'(a,a,a)') "SCALARS ","rho"," float"
    write(lunit, '(a)') "LOOKUP_TABLE default"

    do iSource = 1, nSource
       do i = 1, source(iSource)%surface%nElements
          write(lunit,*) 0.5
       enddo
    enddo
    close(lunit)
666 continue
  end subroutine writeVTKfileSource

  subroutine writeVtkFileAMR(grid, vtkFilename, valueTypeFilename, valueTypeString, xml)
    use input_variables, only : cylindrical, usebinaryxmlvtkfiles
#ifdef MPI
    include 'mpif.h'
#endif
    type(GRIDTYPE) :: grid
    character(len=*) :: vtkFilename
    character(len=100) :: xmlFilename
    integer :: nValueType
    character(len=20) :: valueType(50)
    character(len=*), optional ::  valueTypeFilename
    character(len=*), optional ::  valueTypeString(:)
    character(len=1), pointer :: string(:), headstring(:)
    logical, optional :: xml
    integer :: nCells, nPoints
    integer :: lunit = 69
    integer :: nOctals, nVoxels, i, iType
    integer :: nPointOffset
    logical :: ascii, writeXml
    integer :: iSize(1)
    integer, allocatable :: icell(:)
#ifdef MPI
    integer :: ierr
    integer, allocatable :: iOffsetArray(:)
    integer :: myRank, nThreads, iThread
#endif

    if (useBinaryXMLVtkFiles) then
       xmlFilename = adjustl(vtkFilename)
       xmlFilename = xmlFilename(1:index(xmlFilename," ")-5)//".vtu"
       call writeXMLVtkFileAMR(grid, xmlFilename, valueTypeFilename, valueTypeString)
       goto 666
    endif

    ascii = .true.
    writeXML = .false.
    if (PRESENT(xml)) then
       writeXml = .true.
       ascii = .false.
    endif


       

!
#ifdef MPI
! just return if the grid is not decomposed and MPI job and not zero rank thread
    if ((.not.grid%splitOverMpi).and.(myRankGlobal /= 0)) goto 666
#endif
!
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
       if (.not.grid%splitOverMpi) then
          nValueType = 4
          valueType(1) = "rho"
          valueType(2) = "velocity"
          valueType(3) = "temperature"
          valueType(4) = "inflow"
       else
          nValueType = 4
          valueType(1) = "rho"
          valueType(2) = "velocity"
          valueType(3) = "temperature"
          valueType(4) = "mpithread"
       endif
    endif

    if (PRESENT(valueTypeString)) then
       nValueType = SIZE(valueTypeString)
       valueType(1:nValueType) = valueTypeString(1:nValueType)
    endif


    if (grid%octreeRoot%threed) then
       nPointOffset = 8
    else if (grid%octreeRoot%twoD) then
       nPointOffset = 4
    else if (grid%octreeRoot%oneD) then
       nPointOffset = 2
    endif

    call countVoxels(grid%octreeRoot, nOctals, nVoxels)  
    if (cylindrical) then
       nVoxels = 0
       call countVoxelsCylindrical(grid%octreeRoot, nVoxels)
    endif

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
       if (ascii) then
          open(lunit,file=vtkFilename, form="formatted", status="unknown")
          write(lunit,'(a)') "# vtk DataFile Version 2.0"
          write(lunit,'(a,a)') "TORUS AMR data"
          write(lunit,'(a)') "ASCII"
          write(lunit,'(a)') "DATASET UNSTRUCTURED_GRID"
          write(lunit,'(a)') "FIELD FieldData 2"
          write(lunit,'(a)') "TIME 1 1 double"
          write(lunit,*) grid%currentTime
          write(lunit,'(a)') "CYCLE 1 1 int"
          write(lunit,*) grid%idump
          close(lunit)
       else
          open(lunit,file=vtkFilename, form="formatted", status="unknown")
          write(lunit,'(a)') "<?xml version='1.0'?>"
          write(lunit,'(a)') "<VTKFile xmlns='VTK' type='UnstructuredGrid' version='0.1' byte_order='BigEndian' >"
          write(lunit,'(a)') "<UnstructuredGrid>"
          close(lunit)
       endif

    endif

    if (.not.grid%splitOverMPI) call writePoints(grid, vtkFilename, nPoints, nCells, xml=writexml)

#ifdef MPI
    if (grid%splitOverMpi) then
    do iThread = 1, nThreadsGlobal
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       if (iThread == myRankGlobal) then
          call writePoints(grid, vtkFilename, nPoints, nCells, xml=writexml)
       endif
    enddo
 endif
#endif

    if (.not.grid%splitOverMPI) call writeIndices(grid, vtkFilename, nPoints, nCells, 0, xml=writexml)

#ifdef MPI
    if (grid%splitOverMpi) then

    do iThread = 1, nThreadsGlobal
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       if (iThread == myRankGlobal) then
           call writeIndices(grid, vtkFilename, nPoints, nCells, iOffsetArray(myRankGlobal), xml=writexml)
       endif
    enddo
 endif
#endif

 if (writexml) then
    if (.not.grid%splitOverMPI) call writeOffsetsXML(grid, vtkFilename, 0)

#ifdef MPI
    if (grid%splitOverMpi) then

    do iThread = 1, nThreadsGlobal
       call MPI_BARRIER(amrCOMMUNICATOR, ierr)
       if (iThread == myRankGlobal) then
           call writeOffsetsXML(grid, vtkFilename,  iOffsetArray(myRankGlobal))
       endif
    enddo
 endif
#endif
endif
    if (writeHeader) then
       if (.not.writexml) then
          open(lunit,file=vtkFilename, form="formatted", status="old",position="append")
          write(lunit, '(a, i10)') "CELL_TYPES ",nCells
          do i = 1, nCells
             if ((nPointOffset == 8).and.(.not.cylindrical)) write(lunit, '(a)') "11"
             if ((nPointOffset == 8).and.(cylindrical)) write(lunit, '(a)') "12"
             if (nPointOffset == 4) write(lunit, '(a)') "8"
          enddo
          write(69, '(a,  i10)') "CELL_DATA ",nCells
          close(lunit)
       else
          allocate(iCell(1:nCells))
          do i = 1, nCells
             if ((nPointOffset == 8).and.(.not.cylindrical)) icell(i) = 11
             if ((nPointOffset == 8).and.(cylindrical)) icell(i) = 12
             if (nPointOffset == 4) icell(i) = 8
          enddo

          open(lunit,file=vtkFilename, form="formatted", status="old",position="append")
          write(lunit, '(a)') "<DataArray type='Int32' Name='types' format='binary' >"
          isize(1) = nCells*4
          call base64encode(isize, headstring)
          call base64encode(icell, string)
          write(lunit,'(10000a1)') headstring(1:8),string(1:SIZE(string))
          write(lunit, '(a)') "</DataArray>"
          write(lunit, '(a)') "</Cells>"
          close(lunit)
          deallocate(icell)
       endif
    endif


    if (writeheader.and.writexml) then
       open(lunit,file=vtkFilename, form="formatted", status="old",position="append")
       write(lunit,'(a)') "<CellData>"
       close(lunit)
    endif

    do iType = 1, nValueType
    
       if (.not.grid%splitOverMPI) &
            call writeValue(grid, vtkFilename, valueType(iType), writexml)

#ifdef MPI
    if (grid%splitOverMpi) then

       do iThread = 1, nThreadsGlobal
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          if (iThread == myRankGlobal) then
             call writeValue(grid, vtkFilename, valueType(iType), writexml)
          endif
       enddo
    endif
#endif
    
 enddo

    if (writeheader.and.writexml) then
       open(lunit,file=vtkFilename, form="formatted", status="old",position="append")
       write(lunit,'(a)') "</CellData>"
       close(lunit)
    endif


 if (writeOutput.and.writeXml) then
    open(lunit,file=vtkFilename, form="formatted", status="old", position="append")
    write(lunit,'(a)') "</Piece>"
    write(lunit,'(a)') "</UnstructuredGrid>"
    write(lunit,'(a)') "</VTKFile>"
    close(lunit)
 endif


666 continue

  end subroutine writeVtkFileAMR

!   subroutine writeIfitFile(grid, ifritFilename)
!      use input_variables, only : minDepthAMR, maxDepthAMR
!      type(GRIDTYPE) :: grid
!      character(len=*) :: ifritFilename
!      integer :: fileID, nValueType, nPointOffset
!      integer :: i, iType
!      character(len=20) :: valueType(50)
!      integer :: nCells, nPoints, nOctals, nVoxels
!#ifdef MPI
!      include 'mpif.h'
!      integer :: ierr
!      integer, allocatable :: iOffsetArray(:)
!      integer :: myRank, nThreads, iThread
!#endif
!
!#ifdef MPI
!   if ((.not.grid%splitOverMpi).and.(myRankGlobal /= 0)) goto 101
!#endif
!#ifdef MPI
!   if((grid%splitOverMpi) .and. (myRankGlobal == 0)) goto 101
!#endif
!
!     if(minDepthAMR == maxDepthAMR) then
!        open(fileID, file=ifritFilename, status="new", form="formatted")
!
!        !This is certainly subject to change
!        if(grid%splitOverMpi) then
!           nValueType = 5
!           valueType(1) = "rho"
!           valueType(2) = "velocity"
!           valueType(3) = "temperature"
!           valueType(4) = "HI"
!           valueType(5) = "mpithread"
!        end if
!
!        if(grid%octreeRoot%threeD) then
!           nPointOffset = 8
!           write(fileID, *) 2.**maxDepthAMR, 2.**maxDepthAMR, 2.**maxDepthAMR
!        else if(grid%octreeRoot%twoD) then
!           nPointOffset = 4
!           write(fileID, *) 2.**maxDepthAMR, 2.**maxDepthAMR
!        else if(grid%octreeRoot%oneD) then
!           nPointOffset = 2
!           write(fileID, *) 2.**maxDepthAMR
!        end if
!
!
!        call countVoxels(grid%octreeRoot, nOctals, nVoxels)
!#ifdef MPI
!      if(grid%splitOverMpi) then
!         call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!         call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
!         call MPI_COMM_SIZE(MPI_COMM_WORLD, nThreads, ierr)
!         allocate(iOffsetArray(1:nThreads-1))
!       call countSubcellsMPI(grid, nVoxels, nSubcellArray = iOffsetArray)
!       iOffsetArray(2:nThreads-1) = iOffsetArray(1:nThreads-2)
!       iOffsetArray(1) = 0
!       do i = 2, nThreads-1
!          iOffsetArray(i) = iOffsetArray(i) + iOffsetArray(i-1)
!       enddo
!       iOffsetArray = iOffsetArray * nPointOffset
!    endif
!#endif
!
!
!       nCells = nVoxels
!       nPoints = nCells * nPointOffset
!
!!#ifdef MPI
!!    if (grid%splitOverMpi) then
!!       do iThread = 1, nThreadsGlobal
!!          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!!          if (iThread == myRankGlobal) then
!!             call writePoints(grid, ifritFilename, nPoints)
!!          endif
!!       enddo
!!    endif
!!#endif
!
!!#ifdef MPI
!!       if (grid%splitOverMpi) then
!!!
!          do iThread = 1, nThreadsGlobal
!             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!             if (iThread == myRankGlobal) then
!                call writeIndices(grid, ifritFilename, nPoints, nCells, iOffsetArray(myRankGlobal))
!             endif
!          enddo
!       endif
!#endif
!
!#ifdef MPI
!       if (grid%splitOverMpi) then
!
!          do iThread = 1, nThreadsGlobal
!             call MPI_BARRIER(amrCOMMUNICATOR, ierr)
!             if (iThread == myRankGlobal) then
!                call writeValue(grid, ifritFilename, valueType(iType))
!             endif
!	  end do
!       endif
!#endif
!      else
!      !Write data for test 5 from Cosmological Radiative Transfer Comparison T!est
!         print *, "Not a fixed grid, IFRIT file write failed"
!      end if
!
!#ifdef MPI
!     101 continue
!#endif
!
!   end subroutine writeIfritFile


    recursive subroutine countVoxelsCylindrical(thisOctal, nVoxels)
      type(OCTAL), pointer :: thisOctal, child
      integer :: subcell, i, nVoxels

      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call countVoxelsCylindrical(child, nVoxels)
                  exit
               end if
            end do
         else
            nVoxels = nVoxels + max(1, nint(returndphi(thisOctal)/(10.d0*degtorad)))
         endif
      enddo
    end subroutine countVoxelsCylindrical


  subroutine base64encode(i32Array, string)
    implicit none
    integer :: j 
    integer :: i32array(:)
    integer(kind=1), allocatable :: iArray(:)
    integer(kind=1) :: i6
    integer :: iCount, lineCount
    integer :: i32temp
    character, pointer :: string(:)
    integer :: nBytes, k, npad

    nBytes = SIZE(i32Array)*4
    nPad = 0
    if (mod(nBytes,3) /= 0) then
       nPad = 3-mod(nBytes,3)
    endif
    allocate(iArray(1:nBytes))
    iArray = transfer(i32Array, iArray)
    allocate(string(1:(nBytes*8/6)+10))
!    write(*,*) iArray
    icount = 0
    linecount = 0
    do k = 1, nBytes,3

       if ((k+2) <= nBytes) then
          i32temp = iarray(k)*256**2 +iarray(k+1)*256 + iarray(k+2)
       else
          if ((k+1) == nBytes) then
             i32temp = iarray(k)*256**2 +iarray(k+1)*256 
          else if ((k) == nBytes) then
              i32temp = iarray(k)*256**2
           endif
        endif

       i32temp = ishft(i32temp, 8)
!       write(*,'(a,b32.32)') "i32temp ",i32temp
       do j = 1, 4
          i6 = ibits(i32temp,26,6)
!          write(*,'(a,b6.6)') "i ",i6
          icount = icount + 1
          lineCount = lineCount + 1
          string(icount) = returnbase64char(i6)
          if (lineCount == 76) then
             lineCount = 0
             iCount = iCount + 1
             string(iCount) = char(13)
             iCount = iCount + 1
             string(iCount) = char(10)
          endif
          i32temp = ishft(i32temp, 6)
       enddo
    enddo
    if (nPad == 2) then
!       string((iCount-1):iCount) = "=="
    endif
    if (nPad == 1) then
       string(iCount:iCount) = "="
    endif
  end subroutine base64encode

function returnBase64Char(i) result(c)
  character(len=1) :: c
  integer(kind=1) :: i 

  if ((i >= 0).and.(i<26)) then
     c = char(i + ichar("A")) 
  else if ((i >= 26).and.(i < 52)) then
     c = char(i - 26 + ichar("a"))
  else if ((i >= 52).and.(i < 62)) then
     c = char(i-52+ichar("0"))
  else if (i==62) then
     c = "+"
  else if (i==63) then
     c = "/"
  endif
end function returnBase64Char

  subroutine writeXMLVtkFileAMR(grid, vtkFilename, valueTypeFilename, valueTypeString)
    use input_variables, only : cylindrical
#ifdef MPI
    include 'mpif.h'
#endif
    type(GRIDTYPE) :: grid
    character(len=*) :: vtkFilename
    integer :: nValueType
    character(len=20) :: valueType(50)
    character(len=*), optional ::  valueTypeFilename
    character(len=*), optional ::  valueTypeString(:)
    integer :: nCells, nPoints
    integer :: lunit = 69
    integer :: nOctals, nVoxels, i
    integer :: nPointOffset
    integer(kind=1), allocatable :: cellTypes(:)
    integer, allocatable :: offsets(:)
    integer, allocatable :: connectivity(:)
    real, pointer :: points(:,:)
    integer :: ioff(20), nCellsGlobal
    integer :: nBytesPoints, nBytesCellTypes, nBytesConnect, nBytesOffsets, nBytesData
    real, pointer :: rArray(:,:)
    real :: float
    integer :: int, j, iValues, n
    integer(kind=1) :: int1
    Character(len=200) :: buffer
    character(len=1) :: lf
    character(len=10) :: offset, str1, str2
    logical :: vectorValue, scalarValue
    integer :: sizeOfFloat, sizeOfInt, sizeOfInt1
#ifdef MPI
    integer, allocatable :: nSubcellArray(:), ierr, ithread
#endif
    float = 0.
    int = 0
    int1 = 0

    sizeOfFloat = 4
    sizeofInt = 4
    sizeofint1 = 1
#ifdef MEMCHECK
    sizeOfFloat = sizeof(float)
    sizeOfInt = sizeof(int)
    sizeOfInt1 = sizeof(int1)
#endif

!
#ifdef MPI
! just return if the grid is not decomposed and MPI job and not zero rank thread
    if ((.not.grid%splitOverMpi).and.(myRankGlobal /= 1)) goto 666
#endif
!
#ifdef MPI
! just return if the grid is decomposed and MPI job and this is rank zero thread
    if (grid%splitOverMpi.and.(myRankGlobal == 0)) goto 666
#endif

    lf = char(10)

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
       if (.not.grid%splitOverMpi) then
          nValueType = 4
          valueType(1) = "rho"
          valueType(2) = "velocity"
          valueType(3) = "temperature"
          valueType(4) = "inflow"
       else
          nValueType = 4
          valueType(1) = "rho"
          valueType(2) = "velocity"
          valueType(3) = "temperature"
          valueType(4) = "mpithread"
       endif
    endif

    if (PRESENT(valueTypeString)) then
       nValueType = SIZE(valueTypeString)
       valueType(1:nValueType) = valueTypeString(1:nValueType)
    endif


    if (grid%octreeRoot%threed) then
       nPointOffset = 8
    else if (grid%octreeRoot%twoD) then
       nPointOffset = 4
    else if (grid%octreeRoot%oneD) then
       nPointOffset = 2
    endif


    call countVoxels(grid%octreeRoot, nOctals, nVoxels)  
    if (cylindrical) then
       nVoxels = 0
       call countVoxelsCylindrical(grid%octreeRoot, nVoxels)
    endif
    ncells = nVoxels
    nCellsGlobal = nVoxels

#ifdef MPI
    if (grid%splitOverMpi) then
       allocate(nSubcellArray(1:nHydroThreadsGlobal))
       call countSubcellsMPI(grid, nVoxels, nSubcellArray=nSubcellArray)
       ncells = nSubcellArray(myRankGlobal)
       deallocate(nSubcellArray)
    endif
#endif
    nCellsGlobal = nVoxels
    nPoints = nCellsGlobal * nPointOffset

    writeHeader = .true.
#ifdef MPI
    if (myRankGlobal /= 1) then
       writeHeader = .false.
    endif
#endif

    allocate(points(1:3,1:nPoints))
    call getPoints(grid, nPoints, points)

#ifdef MPI
    if (grid%splitOverMpi)  &
    call concatArrayOverAllThreads(points, nPoints)
#endif


    allocate(offsets(1:nCellsGlobal))
    allocate(connectivity(1:nPoints))
    do i = 1, nPoints
       connectivity(i) = i-1
    enddo

    allocate(cellTypes(1:nCellsGlobal))


    do i = 1, nCellsGlobal
       offsets(i) = i * nPointOffset
    enddo


    if ((nPointOffset == 8).and.(.not.cylindrical)) cellTypes = 11
    if ((nPointOffset == 8).and.(cylindrical)) cellTypes = 12
    if (nPointOffset == 4) cellTypes = 8




    nBytesPoints = sizeofFloat * 3 * nPoints
    nBytesConnect =  sizeofInt * nPoints
    nBytesOffsets = sizeofint * nCellsGlobal
    nBytesCellTypes = sizeofint1 * nCellsGlobal
    

    ioff(1) = 0
    ioff(2) = ioff(1) + sizeofint + nBytesPoints
    ioff(3) = ioff(2) + sizeofint + nBytesConnect
    ioff(4) = ioff(3) + sizeofint + nBytesOffsets
    ioff(5) = ioff(4) + sizeofint + nBytesCellTypes
    
    do i = 1, nValueType
       select case (valueType(i))
       case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom")
          scalarvalue = .false.
          vectorvalue = .true.
       case DEFAULT
          scalarvalue = .true.
          vectorvalue = .false.
       end select
       if (scalarvalue) then
          j =  nCellsGlobal * sizeoffloat
       else
          j =  nCellsGlobal * sizeoffloat * 3
       endif
       ioff(i+5) = ioff(i+4) + sizeofint + j
    enddo


    if (writeheader) then
       open(lunit, file=vtkFilename, form="unformatted",access="stream",status="replace")
       buffer = '<?xml version="1.0"?>'//lf
       write(lunit) trim(buffer)
       buffer = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
       write(lunit) trim(buffer)
       buffer = '  <UnstructuredGrid>'//lf 
       write(lunit) trim(buffer)
       write(str1(1:10),'(i10)') nPoints
       write(str2(1:10),'(i10)') nCellsGlobal
       buffer = '    <Piece NumberOfPoints="'//str1//'" NumberOfCells="'//str2//'">'//lf
       write(lunit) trim(buffer)
       buffer = '      <Points>'//lf
       write(lunit) trim(buffer)
       write(offset(1:10),'(i10)') ioff(1)
       buffer = '       <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" format="appended" offset="'&
            //offset//'" />'//lf
       write(lunit) trim(buffer)
       buffer = '      </Points>'//lf
       write(lunit) trim(buffer)
       buffer = '      <Cells>'//lf
       write(lunit) trim(buffer)
       write(offset(1:10),'(i10)') ioff(2)
       buffer = '        <DataArray type="Int32" Name="connectivity" format="appended" offset="'//offset//'" />'//lf
       write(lunit) trim(buffer)
       write(offset(1:10),'(i10)') ioff(3)
       buffer = '        <DataArray type="Int32" Name="offsets" format="appended" offset="'//offset//'" />'//lf
       write(lunit) trim(buffer)
       write(offset(1:10),'(i10)') ioff(4)
       buffer = '        <DataArray type="UInt8" Name="types" format="appended" offset="'//offset//'" />'//lf
       write(lunit) trim(buffer)
       buffer = '      </Cells>'//lf
       write(lunit) trim(buffer)
       buffer = '      <CellData>'//lf
       write(lunit) trim(buffer)


       do i = 1, nValueType
          select case (valueType(i))
          case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom")
             scalarvalue = .false.
             vectorvalue = .true.
          case DEFAULT
             scalarvalue = .true.
             vectorvalue = .false.
          end select

          write(offset(1:10),'(i10)') ioff(i+4)

          if (scalarvalue) then
             buffer = '        <DataArray type="Float32" Name="'//trim(valueType(i))//&
                  '" format="appended" offset="'//offset//'" />'//lf
             write(lunit) trim(buffer)
          else
             buffer = '        <DataArray type="Float32" NumberOfComponents="3" Name="'//trim(valueType(i))//&
                  '" format="appended" offset="'//offset//'" />'//lf
             write(lunit) trim(buffer)
          endif
       enddo


       buffer = '     </CellData>'//lf
       write(lunit) trim(buffer)
       buffer = '    </Piece>'//lf
       write(lunit) trim(buffer)
       buffer = '  </UnstructuredGrid>'//lf
       write(lunit) trim(buffer)
       buffer = '  <AppendedData encoding="raw">'//lf
       write(lunit) trim(buffer)
       close(lunit)

       open(lunit, file=vtkFilename, form="unformatted", position="append", access="stream")
       buffer = '_'
!       write(lunit) iachar("_")
       write(lunit) buffer(1:1)
       write(lunit) nbytesPoints  , ((points(i,j),i=1,3),j=1,nPoints)
       write(lunit) nbytesConnect , (connectivity(i),i=1,nPoints)
       write(lunit) nbytesOffsets  , (offsets(i),i=1,nCellsGlobal)
       write(lunit) nbytesCellTypes    , (cellTypes(i),i=1,nCellsGlobal)
       close(lunit)
    end if

#ifdef MPI

    do ithread = 1, nHydroThreadsGlobal
       if (iThread == myRankGlobal) then
#endif          


    do iValues = 1, nValueType
       select case (valueType(iValues))
       case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom")
             scalarvalue = .false.
             vectorvalue = .true.
          case DEFAULT
             scalarvalue = .true.
             vectorvalue = .false.
          end select
          if (vectorvalue) then
             allocate(rArray(1:3, 1:nCells))
          else
             allocate(rArray(1:1, 1:nCells))
          endif
          call getValues(grid, valueType(iValues), rarray)
          n = size(rarray,2)
          if (writeheader) then
             open(lunit, file=vtkFilename, form="unformatted", position="append", access="stream")
             if (vectorvalue) then
                nBytesData = ncellsGlobal*3*sizeoffloat
                write(lunit) nBytesData, &
                     ((rArray(i,j),i=1,3),j=1,nCells)
             else
                nBytesData = nCellsGlobal * sizeoffloat
                write(lunit) nBytesData, &
                     (rArray(1,i),i=1,nCells)
             endif
             deallocate(rArray)
             close(lunit)
          else
             open(lunit, file=vtkFilename, form="unformatted", position="append", access="stream")
             if (vectorvalue) then
                nBytesData = ncellsGlobal*3*sizeoffloat
                write(lunit) &
                     ((rArray(i,j),i=1,3),j=1,nCells)
             else
                nBytesData = nCellsGlobal * sizeoffloat
                write(lunit) &
                     (rArray(1,i),i=1,nCells)
             endif
             deallocate(rArray)
             close(lunit)
          endif
       enddo

#ifdef MPI
    endif
    call mpi_barrier(amrCommunicator, ierr)
 enddo
#endif          



       if (writeheader) then
          open(lunit, file=vtkFilename, form="unformatted", position="append", access="stream")
          buffer = lf//'  </AppendedData>'//lf
          write(lunit) trim(buffer)
          buffer = '</VTKFile>'//lf
          write(lunit) trim(buffer)
          endfile(lunit)
          close(lunit)
       endif
       goto 666

666 continue
  end subroutine writeXMLVtkFileAMR


  subroutine getValues(grid, valuetype, rarray)
    type(GRIDTYPE) :: grid
    character(len=*) :: valueType
    real :: rArray(:,:)
    integer :: n

    n = 0 

    call getValuesRecursive(grid, grid%octreeRoot, valueType, rArray, n)

  end subroutine getValues

    recursive subroutine getValuesRecursive(grid, thisOctal, valueType, rArray, n)
      use input_variables, only : lambdaSmooth
      type(OCTAL), pointer :: thisOctal, child
      type(GRIDTYPE) :: grid
      character(len=*) :: valueType
      real :: rArray(:,:)
      integer :: i, subcell, n
      integer, save ::  iLambda
      real(double) :: ksca, kabs, value
      type(VECTOR) :: rVec, vel
      real, parameter :: min_single_prec = 1.0e-37
      logical, save :: firstTime=.true.

!$OMP THREADPRIVATE (firstTime)

      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call getValuesRecursive(grid, child, valueType, rArray, n)
                  exit
               end if
            end do
         else

#ifdef MPI
            if ( (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) .and. grid%splitOverMPI ) cycle
#endif

            n = n + 1

            select case (valueType)
               case("rho")

                  rArray(1, n) = returnPhysicalUnitDensity(thisOctal%rho(subcell))

               case("J=0")
                  rArray(1, n) = thisOctal%molecularlevel(1,subcell)

               case("J=1")
                  rArray(1, n) = thisOctal%molecularlevel(2,subcell)

               case("J=2")
                  rArray(1, n) = thisOctal%molecularlevel(3,subcell)

               case("J=3")
                  rArray(1, n) = thisOctal%molecularlevel(4,subcell)

               case("J=4")
                  rArray(1, n) = thisOctal%molecularlevel(5,subcell)

               case("J=5")
                  rArray(1, n) = thisOctal%molecularlevel(6,subcell)

               case("J=10")
                  rArray(1, n) = thisOctal%molecularlevel(11,subcell)

               case("J=16")
                  rArray(1, n) = thisOctal%molecularlevel(17,subcell)

               case("dI")
                  rArray(1, n) = thisOctal%newmolecularlevel(1,subcell)

               case("dIattenuated")
                  rArray(1, n) = thisOctal%newmolecularlevel(2,subcell)

               case("i0")
                  rArray(1, n) = thisOctal%newmolecularlevel(3,subcell)

! galLon and galLat re-use storage used for i0 and dIattenuated 
               case("galLon")
                  rArray(1, n) = thisOctal%newmolecularlevel(2,subcell)

               case("galLat")
                  rArray(1, n) = thisOctal%newmolecularlevel(3,subcell)

               case("crossing")
                  rArray(1, n) = thisOctal%newmolecularlevel(4,subcell)

               case("tauacrosscell")
                  rArray(1, n) = thisOctal%newmolecularlevel(5,subcell)

               case("tau10")
                  rArray(1, n) = thisOctal%tau(1,subcell)

               case("tau21")
                  rArray(1, n) = thisOctal%tau(2,subcell)

               case("tau32")
                  rArray(1, n) = thisOctal%tau(3,subcell)

               case("tau43")
                  rArray(1, n) = thisOctal%tau(4,subcell)

               case("tau54")
                  rArray(1, n) = thisOctal%tau(5,subcell)

               case("tau65")
                  rArray(1, n) = thisOctal%tau(6,subcell)

               case("tau76")
                  rArray(1, n) = thisOctal%tau(7,subcell)

               case("tau87")
                  rArray(1, n) = thisOctal%tau(8,subcell)

               case("level0error")
                  rArray(1, n) = 10**((real(thisOctal%levelconvergence(1,subcell)) / 6553.6) - 4)

               case("level1error")
                  rArray(1, n) = 10**((real(thisOctal%levelconvergence(2,subcell)) / 6553.6) - 4)

               case("level2error")
                  rArray(1, n) = 10**((real(thisOctal%levelconvergence(3,subcell)) / 6553.6) - 4)

               case("level3error")
                  rArray(1, n) = 10**((real(thisOctal%levelconvergence(4,subcell)) / 6553.6) - 4)

               case("level4error")
                  rArray(1, n) = 10**((real(thisOctal%levelconvergence(5,subcell)) / 6553.6) - 4)

               case("niter")
                  rArray(1, n) = int(thisOctal%convergence(subcell)/100)

               case("nh2")
                  rArray(1, n) = real(thisOctal%nh2(subcell))

               case("convergence")
                  rArray(1, n) = mod(thisOctal%convergence(subcell),1.0)

               case("adot")
                  rArray(1, n) = real(thisOctal%adot(subcell))

               case("slowestlevel")
                  rArray(1, n) = floor(mod(thisoctal%convergence(subcell),100.0))

               case("molabundance")
                  rArray(1, n) = thisOctal%molabundance(subcell)

               case("bnu")
                  rArray(1, n) = real(thisOctal%bnu(1,subcell))

               case("dc0")
                  rArray(1, n) = thisOctal%molecularlevel(1,subcell) * thisOctal%departcoeff(1,subcell)

               case("dc1")
                  rArray(1, n) = thisOctal%molecularlevel(2,subcell) * thisOctal%departcoeff(2,subcell)

               case("dc2")
                  rArray(1, n) = thisOctal%molecularlevel(3,subcell) * thisOctal%departcoeff(3,subcell)

               case("dc3")
                  rArray(1, n) = thisOctal%molecularlevel(4,subcell) * thisOctal%departcoeff(4,subcell)

               case("dc4")
                  rArray(1, n) = thisOctal%molecularlevel(5,subcell) * thisOctal%departcoeff(5,subcell)

               case("jnu10")
                  rArray(1, n) = thisOctal%jnu(1,subcell)

               case("dust1")
                  rArray(1, n) = real(thisOctal%dustTypeFraction(subcell,1))

               case("dust2")
                  rArray(1, n) = real(thisOctal%dustTypeFraction(subcell,2))

               case("bias")
                  rArray(1, n) = real(thisOctal%biasCont3d(subcell))

               case("mpithread")
                  rArray(1, n) = real(thisOctal%mpithread(subcell))

               case("bcond")
                  rArray(1, n) = real(thisOctal%boundaryCondition(subcell))

               case("deltaT")
                  rArray(1, n) = real(thisOctal%temperature(subcell)-thisOctal%oldtemperature(subcell))


               case("scattered")
                  rArray(1, n) = real(thisOctal%scatteredIntensity(subcell,5,3))

               case("hydrovelocity")
                  if (thisOctal%threeD) then
                     rArray(1, n) = real(returnPhysicalUnitSpeed(thisOctal%rhou(subcell)/thisOctal%rho(subcell))/1.e5)
                     rArray(2, n) = real(returnPhysicalUnitSpeed(thisOctal%rhov(subcell)/thisOctal%rho(subcell))/1.e5)
                     rArray(3, n) = real(returnPhysicalUnitSpeed(thisOctal%rhow(subcell)/thisOctal%rho(subcell))/1.e5)
                  else
                     rArray(1, n) = real(returnPhysicalUnitSpeed(thisOctal%rhou(subcell)/thisOctal%rho(subcell))/1.e5)
                     rArray(2, n) = real(returnPhysicalUnitSpeed(thisOctal%rhow(subcell)/thisOctal%rho(subcell))/1.e5)
                     rArray(3, n) = real(returnPhysicalUnitSpeed(thisOctal%rhov(subcell)/thisOctal%rho(subcell))/1.e5)
                  endif

               case("radmom")
                     rArray(1, n) = real(thisOctal%radiationMomentum(subcell)%x/1.d20)
                     rArray(2, n) = real(thisOctal%radiationMomentum(subcell)%y/1.d20)
                     rArray(3, n) = real(thisOctal%radiationMomentum(subcell)%z/1.d20)
               case("velocity")
                     ! stop vectors from showing up in visit if too big
                     if(thisoctal%velocity(subcell)%x .ge. 1.d0) then
                        thisoctal%velocity(subcell)%x = 0.d0
                        thisoctal%velocity(subcell)%y = 0.d0
                        thisoctal%velocity(subcell)%z = 0.d0
                     endif
                     if (thisOctal%threed) then
                        rArray(1, n) = real(thisOctal%velocity(subcell)%x*cspeed/1.e5)
                        rArray(2, n) = real(thisOctal%velocity(subcell)%y*cspeed/1.e5)
                        rArray(3, n) = real(thisOctal%velocity(subcell)%z*cspeed/1.e5)
                     else
                        rArray(1, n) = real(thisOctal%velocity(subcell)%x*cspeed/1.e5)
                        rArray(2, n) = real(thisOctal%velocity(subcell)%z*cspeed/1.e5)
                        rArray(3, n) = real(thisOctal%velocity(subcell)%y*Cspeed/1.e5)
                     endif
              case("cornervel")
                 rVec = subcellCentre(thisOctal, subcell)
                 vel = amrGridVelocity(grid%octreeRoot,rvec,startOctal=thisOctal,&
                      actualSubcell=subcell) 
                 rArray(1, n) = real(vel%x*cspeed/1.e5)
                 rArray(2, n) = real(vel%y*cspeed/1.e5)
                 rArray(3, n) = real(vel%z*cspeed/1.e5)

               case("ne")
                  rArray(1, n) = real(thisOctal%ne(subcell))

               case("pressure")
                  rArray(1, n) = real(thisOctal%pressure_i(subcell))

               case("inflow")
                  if (thisOctal%inflow(subcell)) then
                     rArray(1, n) = 1.
                  else
                     rArray(1, n) = 0.
                  endif

               case("NH2")
                  rArray(1, n) = real( thisOctal%NH2(subcell) )

               case("fixedtemp")
                  if (thisOctal%fixedTemperature(subcell)) then
                     rArray(1, n) = 1.
                  else
                     rArray(1, n) = 0.
                  endif

               case("HI")
                  rArray(1, n) = real(thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon)))

               case("HeI")
                  rArray(1, n) = real(thisOctal%ionfrac(subcell,returnIonNumber("He I", grid%ion, grid%nIon)))

               case("HeII")
                  rArray(1, n) = real(thisOctal%ionfrac(subcell,returnIonNumber("He II", grid%ion, grid%nIon)))

               case("OI")
                  rArray(1, n) = real(thisOctal%ionfrac(subcell,returnIonNumber("O I", grid%ion, grid%nIon)))

               case("OII")
                  rArray(1, n) = real(thisOctal%ionfrac(subcell,returnIonNumber("O II", grid%ion, grid%nIon)))

               case("OIII")
                  rArray(1, n) = real(thisOctal%ionfrac(subcell,returnIonNumber("O III", grid%ion, grid%nIon)))

               case("temperature")
                  rArray(1, n) = real(thisOctal%temperature(subcell))

               case("chiline")
                  rArray(1, n) = max( real(thisOctal%chiline(subcell)), min_single_prec )


               case("microturb")
                  rArray(1, n) = max( real(thisOctal%microturb(subcell)*cspeed/1.e5), min_single_prec )

               case("etaline")
                  rArray(1, n) = max ( real(thisOctal%etaline(subcell)), min_single_prec )

               case("sourceline")
                  rArray(1, n) = max ( real(thisOctal%etaline(subcell)), min_single_prec )/ &
                       max( real(thisOctal%chiline(subcell)), min_single_prec )

               case("tau")
                  if (firstTime) then
                     call locate(grid%lamArray, grid%nLambda, lambdaSmooth, ilambda)
                     firstTime = .false.
                  endif
                  call returnKappa(grid, thisOctal, subcell, ilambda=ilambda,&
                       kappaSca=ksca, kappaAbs=kabs)
                  value = thisOctal%subcellSize * (ksca + kabs)
!                  write(*,*) "sca ",thisOctal%subcellSize * ksca,ilambda,lambdasmooth
                  rArray(1, n) = real(value)

               case("ross")

                  call returnKappa(grid, thisOctal, subcell, rosselandKappa=kabs)
                  value = thisOctal%subcellsize * kabs * thisOctal%rho(subcell) * 1.e10

                  rArray(1, n) = real(value)


               case("dusttype")
                  rArray(1, n) =  real(thisOctal%dustTypeFraction(subcell,1))

               case("etacont")
                  rArray(1, n) = real(thisOctal%etaCont(subcell))

               case("crossings")
                  value = thisOctal%ncrossings(subcell)
                  if (thisOctal%diffusionApprox(subcell)) value = 1.e6
                  rArray(1, n) = real(value)
                  
               case("phi")
                  rArray(1, n) = real(thisOctal%phi_i(subcell))

               case("q_i")
                  rArray(1, n) = real(thisOctal%q_i(subcell))

               case("u_i")
                  rArray(1, n) = real(thisOctal%u_interface(subcell))

               case("rhoe")
                  rArray(1, n) = real(thisOctal%rhoe(subcell))

               case("edens_s")
                  rArray(1, n) = real(thisOctal%photonEnergyDensityFromSource(subcell))

               case("edens_g")
                  rArray(1, n) = real(thisOctal%photonEnergyDensityFromGas(subcell))

               case("rhou")
                  rArray(1, n) = real(thisOctal%rhou(subcell))

               case("rhov")
                  rArray(1, n) = real(thisOctal%rhov(subcell))

               case("rhow")
                  rArray(1, n) = real(thisOctal%rhow(subcell))

               case("q_i-1")
                  rArray(1, n) = real(thisOctal%q_i_minus_1(subcell))

               case("ghosts")
                  if (thisOctal%ghostCell(subcell)) then
                     rArray(1, n) = 1.
                  else
                     rArray(1, n) =  0.
                  endif

               case("edges")
                  if (thisOctal%edgeCell(subcell)) then
                     rArray(1, n) = 1.
                  else
                     rArray(1, n) =  0.
                  endif

               case DEFAULT
                  write(*,*) "Cannot write vtk type ",trim(valueType)
             end select


         endif
      enddo
    end subroutine getValuesRecursive

#ifdef MPI

  subroutine concatArrayOverAllThreads(array, n)
    include 'mpif.h'
    real, pointer :: array(:,:)
    integer :: n, tempInt(1), oldN, nStart, nEnd, i
    real, pointer :: bigArray(:,:)
    real, allocatable :: tempStorage(:)
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 51

    oldn = n
    allocate(tempStorage(1:(size(array,1)*size(array,2))))
    call MPI_ALLREDUCE(n, tempInt, 1, MPI_INTEGER, MPI_SUM,  amrCommunicator, ierr)
    if (myrankGlobal == 1) n = tempInt(1)
    if (myRankGlobal == 1) then
       allocate(bigArray(1:SIZE(array,1), 1:n))
       bigArray(:,1:oldN) = array(:,1:oldN)
       nStart = oldN+1
       nEnd = nStart + oldN - 1
       do i = 2, nHydroThreadsGlobal
          call MPI_RECV(tempStorage, size(array,1)*size(array,2), MPI_REAL, i, tag, MPI_COMM_WORLD, status, ierr)
          array = reshape(tempStorage, shape(array))
          bigArray(1:size(array,1),nStart:nEnd) = array
          nStart = nStart + oldN
          nEnd = nEnd + oldN 
       enddo
       deallocate(array)
       array => bigArray
    else
       tempStorage = reshape(array, (/SIZE(array,1)*size(Array,2)/))
       call MPI_SEND(tempStorage, SIZE(tempStorage), MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
    endif
    call mpi_barrier(amrCommunicator, ierr)
  end subroutine concatArrayOverAllThreads

#endif


end module vtk_mod

