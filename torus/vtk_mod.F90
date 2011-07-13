module vtk_mod
  
!
!

! written by tjh
! module to write VTK format files - see
! http://public.kitware.com/VTK/pdf/file-formats.pdf

  use kind_mod
  use dimensionality_mod
  use constants_mod
  use utils_mod
  use amr_mod
  use mpi_amr_mod
  use messages_mod
  use vector_mod
  use mpi_global_mod
#ifdef PHOTOION
  use ion_mod, only: returnIonNumber
#endif

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
                  if (.not.thisOctal%cylindrical) then
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
                  if (.not.thisOctal%cylindrical) then
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
               if (.not.thisOctal%cylindrical) then
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
               if (.not.thisOctal%cylindrical) then
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
      use inputs_mod, only : lambdasmooth
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
            if (thisOctal%threed.and.thisOctal%cylindrical) then
               nVal = max(1, nint(returndPhi(thisOctal)/(10.d0*degtorad)))
            endif


            do j = 1, nVal
            select case (valueType)
               
               case("sourceCont")
                  write(lunit, *) thisOctal%normSourceContribution(subcell, 1)
               
               case("logRho")
                  write(lunit, *) log10(thisOctal%rho(subcell))
               case("rho")

!                  if(thisOctal%rho(subcell) .lt. 1.d-37) then
!                     write(lunit, *) 1e-37 ! floating underflow warning
!                  else
                     write(lunit, *) thisOctal%rho(subcell)
!                  endif

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

#ifdef PHOTOION
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
#endif

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
                  write(*,*) "Cannot write vtk type A",trim(valueType)
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

  subroutine writeVTKfileNbody(nSource, source, vtkFilename, grid)
    use source_mod
    use mpi_global_mod
    type(GRIDTYPE) :: grid
    integer :: nSource
    type(SOURCETYPE) :: source(:)
    character(len=*) :: vtkFilename
    integer, parameter :: lunit = 37
    integer :: nPoints, nElements, iSource
    integer :: i
    type(VECTOR) :: cVec, aVec, v1, v2, v3, v4, zAxis
    real(double) :: dphi, dtheta
    integer :: nCount
    

#ifdef MPI
    if (myrankGlobal /=1 ) goto 666
#endif
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
          call normalize(cVec)
          cVec = cVec * grid%halfSmallestSubcell/4.d0
          cVec = cVec * source(isource)%radius/(rsol/1.d10)
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

    write(lunit,'(a,a,a)') "SCALARS ","mass"," float"
    write(lunit, '(a)') "LOOKUP_TABLE default"

    do iSource = 1, nSource
       do i = 1, source(iSource)%surface%nElements
          write(lunit,*) source(isource)%mass/msol
       enddo
    enddo

    write(lunit,'(a,a,a)') "SCALARS ","teff"," float"
    write(lunit, '(a)') "LOOKUP_TABLE default"

    do iSource = 1, nSource
       do i = 1, source(iSource)%surface%nElements
          write(lunit,*) source(isource)%teff
       enddo
    enddo

    close(lunit)
#ifdef MPI
666 continue
#endif
  end subroutine writeVTKfileNbody

  subroutine writeVtkFileAMR(grid, vtkFilename, valueTypeFilename, valueTypeString, xml)
    use inputs_mod, only : cylindrical, usebinaryxmlvtkfiles
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
    integer :: iThread
#endif

    if (useBinaryXMLVtkFiles) then
       xmlFilename = adjustl(vtkFilename)
       xmlFilename = xmlFilename(1:index(xmlFilename," ")-5)//".vtu"
       if (.not.grid%octreeRoot%oned) &
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
!    if (grid%splitOverMpi.and.(myRankGlobal == 0)) goto 666
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
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       allocate(iOffsetArray(1:nHydroThreadsGlobal))
       if (myrankGlobal /= 0) call countSubcellsMPI(grid, nVoxels, nSubcellArray = iOffsetArray)
       iOffsetArray(1) = 0
       do i = 2, nHydrothreadsGlobal-1
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
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
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
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
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
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
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
!          call base64encode(isize, headstring)
!          call base64encode(icell, string)
!          write(lunit,'(10000a1)') headstring(1:8),string(1:SIZE(string))
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
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
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
!      use inputs_mod, only : minDepthAMR, maxDepthAMR
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


  subroutine base64encode(writeheader, string, iCount, iArray32, iArray8, iArray64, float32, &
       inputnPad, inputiPad, outnpad, outipad, padend, nBytesHeader)
    implicit none
    logical :: writeheader
    integer, optional :: inputnPad, outnPad
    integer, optional :: nBytesHeader
    logical, optional :: padEnd
    integer :: npad
    integer,optional :: inputipad(:), outiPad(:)
    integer :: j 
    real, optional :: float32(:)
    integer(bigInt), optional :: iarray64(:)
    integer, optional :: iarray32(:)
    integer(kind=1), optional :: iArray8(:)
    integer(kind=1), allocatable :: iArray(:)
    integer(kind=1) :: i6
    integer :: iCount
    integer :: ik, ik1, ik2, i32temp
    character, pointer :: string(:)
    integer :: nBytes, k
    logical :: write32, dopadend

    if (PRESENT(outnpad)) outnpad = 0
    doPadEnd = .true.
    if (present(padend)) dopadend = padend

    if (PRESENT(iArray32)) nBytes = SIZE(iArray32)*4 
    if (PRESENT(float32))  nBytes = SIZE(float32)*4 
    if (PRESENT(iArray64)) nBytes = SIZE(iArray64)*8 
    if (PRESENT(iArray8))  nBytes = SIZE(iArray8)

    if (writeheader) then
       allocate(iArray(1:nBytes+4))
       if (PRESENT(nBytesHeader)) then
          iarray(1:4) = transfer(nBytesHeader, iarray(1:4))
       else
          iarray(1:4) = transfer(nBytes, iarray(1:4))
       endif
       if (PRESENT(iArray64))   iArray(5:) = transfer(iArray64, iArray(5:))
       if (PRESENT(iArray32))   iArray(5:) = transfer(iArray32, iArray(5:))
       if (PRESENT(float32))    iArray(5:) = transfer(float32, iArray(5:))
       if (PRESENT(iArray8))    iArray(5:) = transfer(iArray8, iArray(5:))
       nBytes = nBytes + 4
    endif
    if ((.not.writeheader).and.present(inputnPad)) then
       nBytes = nBytes + inputnPad
       allocate(iArray(1:nBytes))
       if (inputnPad > 0) then
          iArray(1:inputNPad) = inputipad(1:inputnpad)
       endif
       if (PRESENT(iArray64))   iArray((inputnpad+1):) = transfer(iArray64, iArray((inputnpad+1):))
       if (PRESENT(iArray32))   iArray((inputnpad+1):) = transfer(iArray32, iArray((inputnpad+1):))
       if (PRESENT(float32))    iArray((inputnpad+1):) = transfer(float32, iArray((inputnpad+1):))
       if (PRESENT(iArray8))    iArray((inputnpad+1):) = transfer(iArray8, iArray((inputnpad+1):))
    endif
    allocate(string(1:((nBytes+4)*8/6)+20))


    nPad = 0
    if (mod(nBytes,3) /= 0) then
       nPad = 3-mod(nBytes,3)
    endif

!    write(*,*) iArray
    icount = 0
    do k = 1, nBytes,3

       write32 = .true.
       if ((k+2) <= nBytes) then
          i32temp = iarray(k)*256**3 +iarray(k+1)*256**2 + iarray(k+2)*256
          i32temp = 0
          ik = iArray(k)
          ik1 = iArray(k+1)
          ik2 = iarray(k+2)
          call mvbits(ik,0,8,i32temp,24)
          call mvbits(ik1,0,8,i32temp,16)
          call mvbits(ik2,0,8,i32temp,8)
       else
          if ((k+1) == nBytes) then
             i32temp = iarray(k)*256**3 +iarray(k+1)*256**2
             i32temp = 0
             ik = iArray(k)
             ik1 = iArray(k+1)
             call mvbits(ik,0,8,i32temp,24)
             call mvbits(ik1,0,8,i32temp,16)
             if (PRESENT(outiPad)) then
                outnPad = 2
                outipad(1) = iArray(k)
                outipad(2) = iArray(k+1)
                write32 = .false.
             endif
          else if ((k) == nBytes) then
             i32temp = iarray(k)*256**3
             i32temp = 0
             ik = iArray(k)
             call mvbits(ik,0,8,i32temp,24)
             if (PRESENT(outiPad)) then
                outnPad = 1
                outipad(1) = iArray(k)
                write32 = .false.
             endif
           endif
        endif
        if (present(padend)) then
           if (padend) write32 = .true.
        endif
        if (write32) then
           do j = 1, 4
              i6 = ibits(i32temp,26,6)
              icount = icount + 1
              string(icount) = returnbase64char(i6)
              i32temp = ishft(i32temp, 6)
           enddo
        endif
    enddo
    if (dopadEnd) then
       if (nPad == 2) then
          string(iCount:iCount) = "="
          string((iCount-1):(iCount-1)) = "="
       endif
       if (nPad == 1) then
          string(iCount:iCount) = "="
       endif
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
  else
     write(*,*) "encoder called with wrong number ",i
  endif
end function returnBase64Char

subroutine writeXMLVtkFileAMR(grid, vtkFilename, valueTypeFilename, valueTypeString)
#ifdef MPI
  include 'mpif.h'
#endif
  type(GRIDTYPE) :: grid
  character(len=*) :: vtkFilename
  integer :: nValueType
  character(len=20) :: valueType(50)
  character(len=*), optional ::  valueTypeFilename
  character(len=*), optional ::  valueTypeString(:)
  integer :: nCells, nPoints, nPointsGlobal
  integer :: lunit = 69
  integer :: nOctals, nVoxels, i
  integer :: nPointOffset
  integer(kind=1), allocatable :: cellTypes(:)
  integer, allocatable :: offsets(:)
  integer, allocatable :: connectivity(:)
  real, pointer :: points(:,:)
  integer :: nCellsGlobal
  integer :: nBytesPoints, nBytesCellTypes, nBytesConnect, nBytesOffsets
  real, pointer :: rArray(:,:)
  real :: float
  integer :: int,  iValues
  integer(kind=1) :: int1
  Character(len=200) :: buffer
  character(len=1) :: lf
  character, pointer :: pstring(:)
  character(len=12) :: str1, str2
  logical :: vectorValue, scalarValue
  integer :: sizeOfFloat, sizeOfInt, sizeOfInt1
  integer :: nString
  integer(kind=1), allocatable :: itest(:)
  real, allocatable :: float32(:)
  integer(kind=bigInt), allocatable :: itest64(:)
#ifdef MPI
  integer, allocatable :: nSubcellArray(:)
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
10   continue
     read(29,'(a)',end=20) valueType(nValueType)
     nValueType = nValueType + 1
     goto 10
20   continue
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
  if (grid%octreeRoot%cylindrical) then
     nVoxels = 0
     call countVoxelsCylindrical(grid%octreeRoot, nVoxels)
  endif
  ncells = nVoxels
  nCellsGlobal = nVoxels

#ifdef MPI
  if (grid%splitOverMpi) then
     allocate(nSubcellArray(1:nHydroThreadsGlobal))
     if (myrankGlobal /= 0) then
        call countSubcellsMPI(grid, nVoxels, nSubcellArray=nSubcellArray)
        ncells = nSubcellArray(myRankGlobal)
        deallocate(nSubcellArray)
     endif
  endif
#endif
  nCellsGlobal = nVoxels
  nPointsGlobal = nCellsGlobal * nPointOffset
  nPoints = ncells * nPointOffset

  writeHeader = .true.
#ifdef MPI
  if ((myRankGlobal /= 1).and.(grid%splitOverMPI)) then
     writeHeader = .false.
  endif
#endif

  allocate(points(1:3,1:nPoints))
  call getPoints(grid, nPoints, points)


  allocate(offsets(1:nCellsGlobal))
  allocate(connectivity(1:nPointsGlobal))
  do i = 1, nPointsGlobal
     connectivity(i) = i-1
  enddo

  allocate(cellTypes(1:nCellsGlobal))

  if ((nPointOffset == 8).and.(.not.grid%octreeRoot%cylindrical)) cellTypes = 11
  if ((nPointOffset == 8).and.(grid%octreeRoot%cylindrical)) cellTypes = 12
  if (nPointOffset == 4) cellTypes = 8


  do i = 1, nCellsGlobal
     offsets(i) = i * nPointOffset
  enddo


  nBytesPoints = sizeofFloat * 3 * nPointsGlobal
  nBytesConnect =  sizeofInt * nPointsGlobal
  nBytesOffsets = sizeofint * nCellsGlobal
  nBytesCellTypes = sizeofint1 * nCellsGlobal

  if (writeheader) then
     open(lunit, file=vtkFilename, form="unformatted",access="stream",status="replace")
     buffer = '<?xml version="1.0"?>'//lf
     write(lunit) trim(buffer)
     buffer = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
     write(lunit) trim(buffer)
     buffer = '  <UnstructuredGrid>'//lf 
     write(lunit) trim(buffer)
     write(str1(1:12),'(i10)') nPointsGlobal
     write(str2(1:12),'(i10)') nCellsGlobal
     buffer = '    <Piece NumberOfPoints="'//str1//'" NumberOfCells="'//str2//'">'//lf
     write(lunit) trim(buffer)
     buffer = '      <Points>'//lf
     write(lunit) trim(buffer)
     buffer = '       <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" format="binary">'//lf
     write(lunit) trim(buffer)
     close(lunit)
  endif
     

  if (writeheader.and.(.not.grid%splitoverMPI)) then
     allocate(float32(1:(nPointsGlobal*3)))
     float32 = RESHAPE(points, (/SIZE(float32)/))
     call base64encode(writeheader, pstring, nString, float32=float32)
     open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
     write(lunit) pstring(1:nString)
     deallocate(pString)
     deallocate(float32)
     close(lunit)
  else if (grid%splitOverMPI) then
#ifdef MPI
     call writePointsDecomposed(grid, vtkFilename, lunit, nPoints, nPointsGlobal)
#endif
  endif

  if (writeheader) then
     open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)
     buffer = '      </Points>'//lf
     write(lunit) trim(buffer)
     buffer = '      <Cells>'//lf
     write(lunit) trim(buffer)

     buffer = '        <DataArray type="Int64" Name="connectivity" format="binary">'//lf
     write(lunit) trim(buffer)
     allocate(itest64(1:SIZE(connectivity)))
     itest64 = connectivity
     call base64encode(writeheader, pstring, nString, iArray64=iTest64)
     write(lunit) pstring(1:nString)
     deallocate(pString)
     deallocate(itest64)
     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)


     buffer = '        <DataArray type="Int64" Name="offsets" format="binary">'//lf
     write(lunit) trim(buffer)
     allocate(itest64(1:SIZE(offsets)))
     itest64 = offsets
     call base64encode(writeheader, pstring, nString, iArray64=itest64)
     write(lunit) pstring(1:nString)
     deallocate(pString)
     deallocate(itest64)
     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)


     buffer = '        <DataArray type="UInt8" Name="types" format="binary">'//lf
     write(lunit) trim(buffer)

     allocate(itest(1:SIZE(cellTypes)))
     itest = cellTypes
     call base64encode(writeheader, pstring, nString, iArray8=itest)
     write(lunit) pstring(1:nString)
     deallocate(pstring)
     deallocate(itest)
     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)

     buffer = '      </Cells>'//lf
     write(lunit) trim(buffer)
     buffer = '      <CellData>'//lf
     write(lunit) trim(buffer)
     close(lunit)

  endif


  deallocate(points, celltypes, connectivity, offsets)

  do ivalues = 1, nValueType
     select case (valueType(iValues))
     case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom")
        scalarvalue = .false.
        vectorvalue = .true.
     case DEFAULT
        scalarvalue = .true.
        vectorvalue = .false.
     end select

     if (scalarvalue) then
        if (writeheader) then
           open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
           buffer = '        <DataArray type="Float32" Name="'//trim(valueType(iValues))//&
                '" format="binary">'//lf
           write(lunit) trim(buffer)
           close(lunit)
        endif

        if (writeheader.and.(.not.grid%splitOverMPI)) then
           allocate(rArray(1:1, 1:nCellsGlobal))
           call getValues(grid, valueType(iValues), rarray)
           allocate(float32(1:nCellsGlobal))
           float32 = rarray(1,1:nCellsGlobal)
           call base64encode(writeheader, pstring, nString, float32=float32)
           open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
           write(lunit) pstring(1:nString)
           close(lunit)
           deallocate(pString)
           deallocate(float32)
           deallocate(rArray)
        else if (grid%splitOverMPI) then
#ifdef MPI
           call writeDomainDecomposed(valueType(iValues), grid, vtkFilename, lunit, nCells, nCellsGlobal)
#endif
        endif

        if (writeheader) then
           open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
           buffer = lf//'        </DataArray>'//lf
           write(lunit) trim(buffer)
           close(lunit)
        endif

     else
        if (writeheader) then
           open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
           buffer = '        <DataArray type="Float32" NumberOfComponents="3" Name="'//trim(valueType(iValues))//&
                '" format="binary">'//lf
           write(lunit) trim(buffer)
           close(lunit)
        endif

        if (writeheader.and.(.not.grid%splitOverMPI)) then
           allocate(rArray(1:3, 1:nCells))
           call getValues(grid, valueType(iValues), rarray)
           allocate(float32(1:nCellsGlobal*3))
           float32 = RESHAPE(rarray, (/SIZE(float32)/))
           call base64encode(writeheader,pstring, nString, float32=float32)
           open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
           write(lunit) pstring(1:nString)
           close(lunit)
           deallocate(pString)
           deallocate(float32)
           deallocate(rArray)
        else if (grid%splitOverMPI) then
#ifdef MPI
           call writeDomainDecomposed(valueType(iValues), grid, vtkFilename, lunit, nCells, nCellsGlobal)
#endif
        endif
           


        if (writeheader) then
           open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
           buffer = lf//'        </DataArray>'//lf
           write(lunit) trim(buffer)
           close(lunit)
        endif
     endif
  enddo
  if (writeheader) then
     open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
     buffer = '     </CellData>'//lf
     write(lunit) trim(buffer)
     buffer = '    </Piece>'//lf
     write(lunit) trim(buffer)
     buffer = '  </UnstructuredGrid>'//lf
     write(lunit) trim(buffer)
     buffer = '</VTKFile>'//lf
     write(lunit) trim(buffer)
     endfile(lunit)
     close(lunit)
  endif

#ifdef MPI
666 continue
#endif
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
      use inputs_mod, only : lambdaSmooth
      type(OCTAL), pointer :: thisOctal, child
      type(GRIDTYPE) :: grid
      character(len=*) :: valueType
      real :: rArray(:,:)
      integer :: i, subcell, n, iVal, nVal
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



            nVal = 1
            if (thisOctal%threed.and.thisOctal%cylindrical) then
               nVal = max(1, nint(returndPhi(thisOctal)/(10.d0*degtorad)))
            endif

            do iVal = 1, nVal
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

               case("jnu")
                  rArray(1, n) = thisOctal%biasLine3d(subcell)

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

               case("haschild")
                  if (thisOctal%haschild(subcell)) then
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

#ifdef PHOTOION
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
#endif

               case("sourceCont")
                  rArray(1, n) = real(thisOctal%normSourceContribution(subcell, 1))

               case("logRho")
                  rArray(1, n) = log10(returnPhysicalUnitDensity(thisOctal%rho(subcell)))

               case("temperature")
                  rArray(1, n) = real(thisOctal%temperature(subcell))

               case("chiline")
                  rArray(1, n) = max( real(thisOctal%chiline(subcell)), min_single_prec )


               case("microturb")
                  rArray(1, n) = max( real(thisOctal%microturb(subcell)*cspeed/1.e5), min_single_prec )

               case("etaline")
                  rArray(1, n) = max ( real(thisOctal%etaline(subcell)), min_single_prec )

               case("n4")
                  rArray(1, n) = max ( real(thisOctal%atomLevel(subcell,1,4)), min_single_prec )

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
                  write(*,*) "Cannot write vtk type B",trim(valueType)
             end select
          enddo


         endif
      enddo
    end subroutine getValuesRecursive


#ifdef MPI
  subroutine writeDomainDecomposed(valueType, grid, vtkFilename, lunit, nCells, nCellsGlobal)
    include 'mpif.h'
    character(len=*) :: valueType, vtkFilename
    type(GRIDTYPE) :: grid
    integer :: nCells, nCellsGlobal
    integer :: lunit
    integer :: iThread
    logical :: padEnd
    integer :: ierr
    integer :: nString
    integer :: outnPad, outipad(2)
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 51
    integer :: tempInt(3), npad, ipad(2)
    real, pointer :: rArray(:,:)
    real, allocatable :: float32(:)
    character, pointer :: pstring(:)

    logical :: scalarvalue, vectorvalue

     select case (valueType)
     case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom")
        scalarvalue = .false.
        vectorvalue = .true.
     case DEFAULT
        scalarvalue = .true.
        vectorvalue = .false.
     end select
     
     do iThread = 1, nHydroThreadsGlobal
        if (myrankGlobal /= 0) then
        if (iThread == myRankGlobal) then

           nPad = 0
           if (myRankGlobal /= 1) then
              call MPI_RECV(tempInt, 3, MPI_INTEGER, iThread-1, tag, MPI_COMM_WORLD, status, ierr)
              npad = tempInt(1) 
              ipad(1:2) = tempInt(2:3)
           endif
           
           outnpad = 0
           padEnd = .false.
           if (myRankGlobal == nHydroThreadsGlobal) padEnd = .true.
           if (scalarValue) then
              allocate(rArray(1:1, 1:nCells))
              call getValues(grid, valueType, rarray)
              allocate(float32(1:nCells))
              float32 = rarray(1,1:nCells)
              call base64encode(writeheader, pstring, nString, float32=float32, inputNpad=npad, inputipad=ipad, &
                   outnpad=outnpad, outipad=outipad, padEnd=padEnd,nBytesHeader=nCellsGlobal*4)
                 
              open(lunit, file=vtkfilename, form="unformatted",access="stream",status="old",position="append")
              write(lunit) pstring(1:nString)
              close(lunit)
              deallocate(pString)
              deallocate(float32)
              deallocate(rArray)
           else
              allocate(rArray(1:3, 1:nCells))
              call getValues(grid, valueType, rarray)
              allocate(float32(1:nCells*3))
              float32 = RESHAPE(rarray, (/SIZE(float32)/))
              call base64encode(writeheader, pstring, nString, float32=float32, inputNpad=npad, inputipad=ipad, &
                   outnpad=outnpad, outipad=outipad, padEnd=padEnd, nBytesHeader=nCellsGlobal*4*3)

              open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
              write(lunit) pstring(1:nString)
              close(lunit)
              deallocate(pString)
              deallocate(float32)
              deallocate(rArray)
           endif

           if (myrankGlobal < nHydroThreadsGlobal) then
              tempInt(1) = outnpad
              tempInt(2:(1+outnPad)) = outipad(1:outnPad)
              call MPI_SEND(tempInt, 3, MPI_INTEGER, iThread+1, tag, MPI_COMM_WORLD, ierr)
           endif
        endif
        endif
     end do
        call MPI_BARRIER(amrCommunicator, ierr)
   end subroutine writeDomainDecomposed

  subroutine writePointsDecomposed(grid, vtkFilename, lunit, nPoints, nPointsGlobal)
    include 'mpif.h'
    character(len=*) ::vtkFilename
    type(GRIDTYPE) :: grid
    integer :: nPoints, nPointsGlobal
    integer :: lunit
    integer :: iThread
    logical :: padEnd
    integer :: ierr
    integer :: nString
    integer :: outnPad, outipad(2)
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 51
    integer :: tempInt(3), npad, ipad(2)
    real, pointer :: rArray(:,:)
    real, allocatable :: float32(:)
    character, pointer :: pstring(:)

     
     do iThread = 1, nHydroThreadsGlobal
        if (myrankGlobal /= 0) then
           if (iThread == myRankGlobal) then
              
              nPad = 0
              if (myRankGlobal /= 1) then
                 call MPI_RECV(tempInt, 3, MPI_INTEGER, iThread-1, tag, MPI_COMM_WORLD, status, ierr)
                 npad = tempInt(1) 
                 ipad(1:2) = tempInt(2:3)
              endif
              
              outnpad = 0
              padEnd = .false.
              if (myRankGlobal == nHydroThreadsGlobal) padEnd = .true.
              allocate(rArray(1:3, 1:nPoints))
              call getPoints(grid, nPoints, rArray)
              allocate(float32(1:nPoints*3))
              float32 = RESHAPE(rarray, (/SIZE(float32)/))
              call base64encode(writeheader, pstring, nString, float32=float32, inputNpad=npad, inputipad=ipad, &
                   outnpad=outnpad, outipad=outipad, padEnd=padEnd, nBytesHeader=nPointsGlobal*4*3)
              
              open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
              write(lunit) pstring(1:nString)
              close(lunit)
              deallocate(pString)
              deallocate(float32)
              
              if (myrankGlobal < nHydroThreadsGlobal) then
                 tempInt(1) = outnpad
                 tempInt(2:(1+outnPad)) = outipad(1:outnPad)
                 call MPI_SEND(tempInt, 3, MPI_INTEGER, iThread+1, tag, MPI_COMM_WORLD, ierr)
              endif
           endif
        endif
     end do
     call MPI_BARRIER(amrCommunicator, ierr)
   end subroutine writePointsDecomposed
#endif
end module vtk_mod

