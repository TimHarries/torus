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
  use zlib_mod
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
          write(lunit,'(a,i12,a)') "POINTS ",nPoints, " float"
       else
          write(lunit,'(a,i12,a,i12,a)') "<Piece NumberOfPoints='",nPoints,"' NumberOfCells='",nCells,"'>"
          write(lunit,'(a)') "<Points>"
          write(lunit,'(a)') "<DataArray type='Float32' Name='Position' NumberOfComponents='3' format='ascii'>"
       endif
    endif


    call recursiveWritePoints(grid%octreeRoot, lunit, grid)
    if (writeheader.and.xml) then
       write(lunit,'(a)') "</DataArray>"
       write(lunit,'(a)') "</Points>"
    endif
    close(lunit)
       

  contains

    recursive subroutine recursiveWritePoints(thisOctal,lunit, grid)
      use octal_mod, only: returndPhi
#ifdef MPI
      use mpi
#endif
      use inputs_mod, only : vtkIncludeGhosts
      type(OCTAL), pointer :: thisOctal, child
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
            if (associated(thisOctal%ghostCell)) then
               if (thisOctal%ghostCell(subcell).and.(.not.vtkincludeGhosts)) cycle
            endif
#endif
            if (thisOctal%threed) then
               if (.not.thisOctal%cylindrical) then
                  rVec = subcellCentre(thisOctal,subcell)
                  d = real(thisOctal%subcellSize/2.d0)
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
                  d = real(thisOctal%subcellSize/2.d0)
                  zp = REAL(rVec%z + d)
                  zm = REAL(rVec%z - d)
                  r1 = real(sqrt(rVec%x**2 + rVec%y**2) - d)
                  r2 = real(sqrt(rVec%x**2 + rVec%y**2) + d)
                  phi = real(atan2(rVec%y, rVec%x))
                  dphi = real(returndPhi(thisOctal))
                  phi1 = real(phi - dphi)
                  phi2 = real(phi + dphi)
                  nPhi = max(1, nint(dphi/(10.d0*degtorad)))
                  do iPhi = 1, nPhi
                     phiStart = real(phi1 + (phi2 - phi1) * dble(iphi-1)/dble(nphi))
                     phiEnd   = real(phi1 + (phi2 - phi1) * dble(iphi)/dble(nphi))

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
               d = real(thisOctal%subcellSize/2.d0)
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

  subroutine getPoints(grid, nPoints, points, includeGhosts)
    type(GRIDTYPE) :: grid
    real :: points(:,:)
    integer(bigint) :: nPoints
    logical :: includeGhosts

    nPoints = 0

    call recursiveGetPoints(grid, grid%octreeRoot, nPoints, Points, includeGhosts)

  contains

    recursive subroutine recursiveGetPoints(grid, thisOctal,nPoints, points, includeGhosts)
      use octal_mod, only: returndPhi
      use inputs_mod, only : hydrodynamics
#ifdef MPI
      use mpi
#endif

      type(OCTAL), pointer :: thisOctal, child
      logical :: includeGhosts
      type(GRIDTYPE) :: grid
      real :: points(:,:)
      integer(bigint) :: nPoints
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
                  call recursiveGetPoints(grid, child, nPoints, points, includeGhosts)
                  exit
               end if
            end do
         else

#ifdef MPI
            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal) .and. grid%splitOverMPI) cycle
            if (hydrodynamics) then
               if (checkGhost(thisOctal, subcell).and.(.not.includeGhosts))cycle
            endif
#endif
            if (thisOctal%threed) then
               if (.not.thisOctal%cylindrical) then
                  rVec = subcellCentre(thisOctal,subcell)
                  d = real(thisOctal%subcellSize/2.d0)
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
                  d = real(thisOctal%subcellSize/2.d0)
                  zp = REAL(rVec%z + d)
                  zm = REAL(rVec%z - d)
                  r1 = real(sqrt(rVec%x**2 + rVec%y**2) - d)
                  r2 = real(sqrt(rVec%x**2 + rVec%y**2) + d)
                  phi = real(atan2(rVec%y, rVec%x))
                  dphi = real(returndPhi(thisOctal))
                  phi1 = phi - dphi
                  phi2 = phi + dphi
                  nPhi = max(1, nint(dphi/(10.d0*degtorad)))
                  do iPhi = 1, nPhi
                     phiStart = real(phi1 + (phi2 - phi1) * dble(iphi-1)/dble(nphi))
                     phiEnd   = real(phi1 + (phi2 - phi1) * dble(iphi)/dble(nphi))


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
               d = real(thisOctal%subcellSize/2.d0)
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
#ifdef MPI
      use mpi
#endif
      use inputs_mod, only : vtkIncludeGhosts, hydrodynamics
      type(OCTAL), pointer :: thisOctal, child
      logical :: xml
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

            if (hydrodynamics) then
               if (thisOctal%ghostCell(subcell).and.(.not.vtkincludeGhosts))cycle
            endif
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
#ifdef MPI
      use mpi
#endif
      use inputs_mod, only : vtkIncludeGhosts, hydrodynamics
      type(OCTAL), pointer :: thisOctal, child

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
            if (hydrodynamics) then
               if (checkGhost(thisOctal, subcell).and.(.not.vtkincludeGhosts))cycle
            endif

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
#ifdef MPI
      use mpi
#endif
      type(OCTAL), pointer :: thisOctal, child

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
       case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom","radforce","fvisc")
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
#ifdef MPI
      use mpi
#endif
      type(OCTAL), pointer :: thisOctal, child
      logical :: xml
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
!      real :: u2, cs, eThermal

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
                          0.
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
                  value = real(thisOctal%subcellSize * (ksca + kabs))
                  write(lunit, *) real(value)

               case("ross")

                  call returnKappa(grid, thisOctal, subcell, rosselandKappa=kabs)
                  value = real(thisOctal%subcellsize * kabs * thisOctal%rho(subcell) * 1.e10)

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

               case("phigas")
                  write(lunit, *) real(thisOctal%phi_gas(subcell))

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
               case("diff")
                  if (thisOctal%diffusionApprox(subcell)) then
                     write(lunit,*)  1.
                  else
                     write(lunit,*) 0.
                  endif


               case DEFAULT
                  write(*,*) "Cannot write ascii vtk type ",trim(valueType)
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
    

    if (myrankWorldGlobal /=0 ) goto 666
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

  subroutine writeVTKfileNbody(nSource, source, vtkFilename)
    use inputs_mod, only : donBodyOnly
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
    

#ifdef MPI
    if (myrankWorldGlobal /=1 ) goto 666
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
          if (donBodyOnly) then
             cVec = cvec * 0.01d0
          else
             cVec = cVec * source(isource)%accretionRadius/1.d10
          endif
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

  subroutine writeVtkFileAMR(grid, RootvtkFilename, valueTypeFilename, valueTypeString, xml)
    use inputs_mod, only : cylindrical, usebinaryxmlvtkfiles, parallelVTUfiles, noVtkGrid
    use inputs_mod, only : iModel
    use utils_mod, only : findMultiFilename

#ifdef MPI
    use mpi
#endif
    type(GRIDTYPE) :: grid
    character(len=*) :: rootVTKfilename
    character(len=80) :: vtkFilename
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

    if(noVtkGrid) return

    call findMultiFilename(rootVTKFilename, iModel, vtkfilename)


    if (useBinaryXMLVtkFiles) then
       xmlFilename = adjustl(vtkFilename)
       xmlFilename = xmlFilename(1:index(xmlFilename," ")-5)//".vtu"
       if (.not.grid%octreeRoot%oned) then
          if (.not.parallelVTUfiles) then
             call writeXMLVtkFileAMR(grid, xmlFilename, valueTypeFilename, valueTypeString)
          else
             call writeParallelXMLVtkFileAMR(grid, vtkFilename, valueTypeFilename, valueTypeString)
          endif
       endif
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
    if (grid%splitOverMpi.and.(myHydroSetGlobal /=0)) goto 666
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
       call MPI_BARRIER(amrCommunicator, ierr)
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
       call MPI_BARRIER(amrCommunicator, ierr)
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
       call MPI_BARRIER(amrCommunicator, ierr)
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
       call MPI_BARRIER(amrCommunicator, ierr)
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
          call MPI_BARRIER(amrCommunicator, ierr)
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
!      use mpi
!      type(GRIDTYPE) :: grid
!      character(len=*) :: ifritFilename
!      integer :: fileID, nValueType, nPointOffset
!      integer :: i, iType
!      character(len=20) :: valueType(50)
!      integer :: nCells, nPoints, nOctals, nVoxels
!#ifdef MPI
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

   subroutine convertAndCompress(iBytes, iHeader, iArray64, iArray32, fArray32, iArray8)
     integer(kind=bigInt), optional :: iArray64(:)
     integer, optional :: iArray32(:)
     integer(kind=1), optional :: iArray8(:)
     real, optional :: fArray32(:)
     integer(kind=1), pointer :: iBytes(:), iHeader(:),iBytesUncompressed(:)
     integer (bigint) :: nBytes, nBytesUncompressed
     integer(kind=4) :: i
     integer(bigint), allocatable :: sizeCompressedBlock(:)
     integer(kind=1), pointer :: thisBlock(:), iTemp(:)
     integer(bigint) :: iCurrent
     integer :: blockSize, nBlocks, lastBlockSize
     integer(bigint) :: iStart, iEnd
     integer(kind=1), pointer :: compressedBlock(:) => null()
     integer :: i4
     integer(kind=bigint) :: nBig

     if (associated(iBytes)) deallocate(iBytes)
     if (associated(iHeader)) deallocate(iHeader)

     if (PRESENT(iArray64)) then
        nBytes = 8 * SIZE(iarray64)
        allocate(iBytesUncompressed(1:nBytes))
        iBytesUncompressed = transfer(iArray64, iBytesUnCompressed)
     endif

     if (PRESENT(iArray32)) then
        nBytes = 4 * SIZE(iarray32)
        allocate(iBytesUncompressed(1:nBytes))
        iBytesUncompressed = transfer(iArray32, iBytesUnCompressed)
     endif

     if (PRESENT(iArray8)) then
        nBytes = SIZE(iarray8)
        allocate(iBytesUncompressed(1:nBytes))
        iBytesUncompressed = transfer(iArray8, iBytesUnCompressed)
     endif

     if (PRESENT(fArray32)) then
        nBytes = 4 * SIZE(farray32)
        allocate(iBytesUncompressed(1:nBytes))
        iBytesUncompressed = transfer(fArray32, iBytesUnCompressed)
     endif

     nBytesUncompressed = SIZE(iBytesUncompressed, kind=bigint)
     nbig = 2**16
!     blockSize = max(nBytesUncompressed/4, nbig)
     blockSize = 2**15
     nBlocks = int(nBytesUncompressed / blockSize)
     lastBlockSize = int(nBytesUncompressed - nBlocks * blockSize )
     if (lastBlockSize > 0) nBlocks = nBlocks + 1
     if (lastBlockSize == 0) lastBlockSize = blockSize


     allocate(sizeCompressedBlock(1:nBlocks))

     allocate(iTemp(1:nBytesUncompressed))
     iCurrent = 1
!     allocate(compressedBlock(1:blockSize))
     do i = 1, nBlocks
        if (i < nBlocks) then
           iStart = 1 + (i-1)*blockSize
           iEnd = iStart + blockSize - 1
           allocate(thisBlock(blockSize))
        else
           iStart = 1 + (i-1)*blockSize
           iEnd = iStart + lastBlockSize - 1
           allocate(thisBlock(lastBlockSize))
        endif
        thisBlock = iBytesUncompressed(iStart:iEnd)

#ifdef USEZLIB
        call compressBytes(thisBlock, SIZE(thisBlock, kind=bigint), compressedBlock, &
             sizeCompressedBlock(i))
#endif
        iTemp(iCurrent:iCurrent+SizeCompressedBlock(i)-1) = compressedBlock(1:sizeCompressedBlock(i))
        deallocate(thisBlock)
        iCurrent = iCurrent + sizeCompressedblock(i)
        deallocate(compressedBlock)
     enddo
     allocate(iBytes(1:SUM(sizeCompressedBlock(1:nBlocks))))
     iBytes(1:SIZE(ibytes)) = itemp(1:SIZE(iBytes))
     deallocate(iTemp)

     allocate(iHeader(1:(12 + nBlocks*4)))
     iHeader(1 : 4) = transfer(nBlocks, iHeader(1 : 4))
     iHeader(5 : 8) = transfer(blockSize, iHeader(5 : 8))
     iHeader(9 :12) = transfer(lastBlockSize, iHeader(9:12))
     do i = 1, nBlocks
        i4 = int(sizeCompressedBlock(i))
        iHeader(13 + (i-1)* 4: 13 + i*4 - 1) = transfer(i4, &
             iHeader(13 + (i-1)* 4: 13 + i*4 - 1))
     enddo
     deallocate(sizeCompressedBlock)
     deallocate(iBytesUncompressed)

   end subroutine convertAndCompress

   subroutine convertAndCompressInSections(iBytes, iHeader, iArray64, iArray32, fArray32, iArray8, firstTime, lastTime)
     integer(kind=bigInt), optional :: iArray64(:)
     integer, optional :: iArray32(:)
     integer(kind=1), optional :: iArray8(:)
     logical, optional :: firstTime, lastTime
     logical :: doFirstTime, doLastTime
     real, optional :: fArray32(:)
     integer(kind=1), pointer :: iBytes(:), iHeader(:),iBytesUncompressed(:)
     integer (bigint) :: nBytes, nBytesUncompressed
     
     integer(kind=4) :: i
     integer(bigint), allocatable :: sizeCompressedBlock(:)
     integer(kind=1), pointer :: thisBlock(:), iTemp(:)
     integer(bigint) :: iCurrent
     integer, parameter :: blockSize = 2**15
     integer(kind=1), save :: leftOverBytes(blockSize)
     integer, save :: nLeftOverBytes
     integer :: nBlocks, lastBlockSize
     integer(bigint) :: iStart, iEnd
     integer(kind=1), pointer :: compressedBlock(:) => null()
     integer :: i4
     integer, save :: totalBlocks
     integer, allocatable, save :: savedBlockSize(:)
     integer, allocatable :: tempBlocks(:)
     integer(kind=1), pointer :: newIBytes(:) => null()
     
     doFirstTime = .false.
     doLastTime = .false.
     if (PRESENT(firstTime)) doFirstTime = firstTime
     if (PRESENT(lastTime)) doLastTime = lastTime
     
     if (doFirstTime) then
        if (associated(iBytes)) deallocate(iBytes)
        nLeftOverBytes = 0
     endif
     
     if (associated(iHeader)) deallocate(iHeader)

     if (PRESENT(iArray64)) then
        nBytes = 8 * SIZE(iarray64) + nLeftOverBytes
        allocate(iBytesUncompressed(1:nBytes))
        if (nLeftOverBytes > 0) iBytesUncompressed(1:nLeftOverBytes) = leftOverBytes(1:nLeftOverBytes)
        iBytesUncompressed(nLeftOverBytes+1:nBytes) = transfer(iArray64, iBytesUnCompressed(nLeftOverBytes+1:nBytes))
     endif

     if (PRESENT(iArray32)) then
        nBytes = 4 * SIZE(iarray32) + nLeftOverBytes
        allocate(iBytesUncompressed(1:nBytes))
        if (nLeftOverBytes > 0) iBytesUncompressed(1:nLeftOverBytes) = leftOverBytes(1:nLeftOverBytes)
        iBytesUncompressed(nLeftOverBytes+1:nBytes) = transfer(iArray32, iBytesUnCompressed(nLeftOverBytes+1:nBytes))
     endif

     if (PRESENT(iArray8)) then
        nBytes = SIZE(iarray8) + nLeftOverBytes
        allocate(iBytesUncompressed(1:nBytes))
        if (nLeftOverBytes > 0) iBytesUncompressed(1:nLeftOverBytes) = leftOverBytes(1:nLeftOverBytes)
        iBytesUncompressed(nLeftOverBytes+1:nBytes) = transfer(iArray8, iBytesUnCompressed(nLeftOverBytes+1:nBytes))
     endif

     if (PRESENT(fArray32)) then
        nBytes = 4 * SIZE(farray32) + nLeftOverBytes
        allocate(iBytesUncompressed(1:nBytes))
        if (nLeftOverBytes > 0) iBytesUncompressed(1:nLeftOverBytes) = leftOverBytes(1:nLeftOverBytes)
        iBytesUncompressed(nLeftOverBytes+1:nBytes) = transfer(fArray32, iBytesUnCompressed(nLeftOverBytes+1:nBytes))
     endif

     nBytesUncompressed = SIZE(iBytesUncompressed, kind=bigint)
     nBlocks = int(nBytesUncompressed / blockSize)
     lastBlockSize = int(nBytesUncompressed - nBlocks * blockSize )

     if (doLastTime) then
        if (lastBlockSize > 0) nBlocks = nBlocks + 1
     else
        if (lastBlockSize > 0) then
           nLeftOverBytes = lastBlockSize
           leftOverBytes(1:nLeftOverBytes) = iBytesUncompressed(nBytesUnCompressed - lastBlockSize + 1:nBytesUncompressed)
        else
           nLeftOverBytes = 0
        endif
        lastBlockSize = blocksize
     endif

     if (lastBlockSize == 0) lastBlockSize = blockSize


     allocate(sizeCompressedBlock(1:nBlocks))

     allocate(iTemp(1:nBytesUncompressed))
     iCurrent = 1

     do i = 1, nBlocks
        if (i < nBlocks) then
           iStart = 1 + (i-1)*blockSize
           iEnd = iStart + blockSize - 1
           allocate(thisBlock(blockSize))
        else
           iStart = 1 + (i-1)*blockSize
           iEnd = iStart + lastBlockSize - 1
           allocate(thisBlock(lastBlockSize))
        endif
        thisBlock = iBytesUncompressed(iStart:iEnd)

#ifdef USEZLIB
        call compressBytes(thisBlock, SIZE(thisBlock, kind=bigint), compressedBlock, &
             sizeCompressedBlock(i))
#endif
        iTemp(iCurrent:iCurrent+SizeCompressedBlock(i)-1) = compressedBlock(1:sizeCompressedBlock(i))
        deallocate(thisBlock)
        iCurrent = iCurrent + sizeCompressedblock(i)
        deallocate(compressedBlock)
     enddo
     if (doFirstTime) then
        nBytes = SUM(sizeCompressedBlock(1:nBlocks))
        allocate(iBytes(1:nBytes))
        iBytes(1:nBytes) = itemp(1:nBytes)        
     else
        nBytes = SUM(sizeCompressedBlock(1:nBlocks))
        allocate(newIBytes(1:size(iBytes)+nBytes))
        newIbytes(1:size(iBytes)) = iBytes(1:size(iBytes))
        newiBytes(size(iBytes)+1:size(iBytes)+nBytes) = &
             iTemp(1:nBytes)
        deallocate(iBytes)
        iBytes => newIbytes
     endif
     deallocate(iTemp)

     if (doFirstTime) then
        if (allocated(savedBlockSize)) deallocate(savedBlockSize)
        totalBlocks = nBlocks
        allocate(savedBlockSize(1:totalBlocks))
        savedBlockSize(1:totalBlocks) = int(sizeCompressedBlock(1:nBlocks))
     else
        allocate(tempBlocks(1:totalBlocks+nBlocks))
        tempBlocks(1:totalBlocks) = savedBlockSize(1:totalBlocks)
        deallocate(savedBlockSize)
        allocate(savedBlockSize(1:totalBlocks+nBlocks))
        savedBlockSize(1:totalBlocks) = tempBlocks(1:totalBlocks)
        savedBlockSize(totalBlocks+1:totalBlocks+nBlocks) = int(sizeCompressedBlock(1:nBlocks))
        totalBlocks = totalBlocks + nBlocks
        deallocate(tempblocks)
     endif


     allocate(iHeader(1:(12 + totalBlocks*4)))
     iHeader(1 : 4) = transfer(totalBlocks, iHeader(1 : 4))
     iHeader(5 : 8) = transfer(blockSize, iHeader(5 : 8))
     iHeader(9 :12) = transfer(lastBlockSize, iHeader(9:12))
     do i = 1, totalBlocks
        i4 = savedBlockSize(i)
        iHeader(13 + (i-1)* 4: 13 + i*4 - 1) = transfer(i4, &
             iHeader(13 + (i-1)* 4: 13 + i*4 - 1))
     enddo
     deallocate(sizeCompressedBlock)
     deallocate(iBytesUncompressed)
     if (PRESENT(lastTime)) deallocate(savedBlockSize)

   end subroutine convertAndCompressInSections

  subroutine base64encode(writeheader, string, iCount, iArray32, iArray8, iArray64, float32, &
       inputnPad, inputiPad, outnpad, outipad, padend, nBytesHeader, debug)
    implicit none
    logical :: writeheader
    integer, optional :: inputnPad, outnPad
    integer, optional :: nBytesHeader
    logical, optional :: debug
    logical :: writeDebug
    logical, optional :: padEnd
    integer :: npad
    integer,optional :: inputipad(:), outiPad(:)
    integer(bigint) :: j 
    real, optional :: float32(:)
    integer(bigInt), optional :: iarray64(:)
    integer, optional :: iarray32(:)
    integer(kind=1), optional :: iArray8(:)
    integer(kind=1), allocatable :: iArray(:)
    integer(kind=1) :: i6
    integer(bigint) :: iCount
    integer :: ik, ik1, ik2, i32temp, i4
    character, pointer :: string(:)
    integer(bigint) :: nBytes, k
    logical :: write32, dopadend

    writeDebug = .false.
    if (PRESENT(debug)) writeDebug=debug

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
          i4 = int(nBytes)
          iarray(1:4) = transfer(i4, iarray(1:4))
       endif
       if (PRESENT(iArray64))   iArray(5:) = transfer(iArray64, iArray(5:))
       if (PRESENT(iArray32))   iArray(5:) = transfer(iArray32, iArray(5:))
       if (PRESENT(float32))    iArray(5:) = transfer(float32, iArray(5:))
       if (PRESENT(iArray8))    iArray(5:) = transfer(iArray8, iArray(5:))
       nBytes = nBytes + 4
    else if ((.not.writeheader).and.present(inputnPad)) then
       nBytes = nBytes + inputnPad
       allocate(iArray(1:nBytes))
       if (inputnPad > 0) then
          iArray(1:inputNPad) = int(inputipad(1:inputnpad),kind=1)
       endif
       if (PRESENT(iArray64))   iArray((inputnpad+1):) = transfer(iArray64, iArray((inputnpad+1):))
       if (PRESENT(iArray32))   iArray((inputnpad+1):) = transfer(iArray32, iArray((inputnpad+1):))
       if (PRESENT(float32))    iArray((inputnpad+1):) = transfer(float32, iArray((inputnpad+1):))
       if (PRESENT(iArray8))    iArray((inputnpad+1):) = iarray8
    else
       allocate(iArray(1:nBytes))
       if (PRESENT(iArray64))   iArray = transfer(iArray64, iArray)
       if (PRESENT(iArray32))   iArray = transfer(iArray32, iArray)
       if (PRESENT(float32))    iArray = transfer(float32, iArray)
       if (PRESENT(iArray8))    iArray = iarray8
    endif
    if (writedebug) then
       write(*,*) "last byte of iarray8 ",iarray8(SIZE(iarray8))
       write(*,*) "last byte of iArray ",iArray(SIZE(iarray))
    endif

    allocate(string(1:((nBytes+4)*8/6)+20))

    if (writeDebug) write(*,*) "nBytes ",nBytes

    nPad = 0
    if (mod(nBytes,3_bigint) /= 0) then
       nPad = 3-int(mod(nBytes,3_bigint))
    endif

    if (writeDebug) write(*,*) "nPad ",nPad

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
          if (writeDebug.and.((k+2)==nBytes)) then
             write(*,*) "ik, ik1, ik2 ",ik,ik1,ik2
             write(*,*) "last (k+2 = nbytes)",i32temp
          endif

       else
          if ((k+1) == nBytes) then
             i32temp = iarray(k)*256**3 +iarray(k+1)*256**2
             i32temp = 0
             ik = iArray(k)
             ik1 = iArray(k+1)
             call mvbits(ik,0,8,i32temp,24)
             call mvbits(ik1,0,8,i32temp,16)
             if (writeDebug) write(*,*) "last (k+1 = nbytes)",i32temp
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
             if (writeDebug) write(*,*) "last (k=nbytes)",i32temp
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
              i6 = int(ibits(i32temp,26,6),kind=1)
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
    deallocate(iArray)
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
  use mpi
#endif
  use inputs_mod, only : vtkIncludeGhosts
  type(GRIDTYPE) :: grid
  character(len=*) :: vtkFilename
  integer :: nValueType
  character(len=20) :: valueType(50)
  character(len=*), optional ::  valueTypeFilename
  character(len=*), optional ::  valueTypeString(:)
  integer :: nCells
  integer(bigint) :: nPointsGlobal, nPoints
  integer :: lunit = 69
  integer :: nOctals, nVoxels
  integer(bigint) :: i
  integer :: nPointOffset
  integer(kind=1), allocatable :: cellTypes(:)
  integer(kind=bigint), allocatable :: offsets(:)
  integer(kind=bigint), allocatable :: connectivity(:)
  real, pointer :: points(:,:)
  integer :: nCellsGlobal
  integer(kind=bigint) :: nBytesPoints, nBytesCellTypes, nBytesConnect, nBytesOffsets
  real, pointer :: rArray(:,:)
  real :: float
  integer :: int,  iValues
  integer(kind=1) :: int1
  Character(len=200) :: buffer
  character(len=1) :: lf
  character, pointer :: pstring(:), pstring2(:)
  character(len=12) :: str1, str2
  logical :: vectorValue, scalarValue
  integer :: sizeOfFloat, sizeOfInt, sizeOfInt1
  integer(kind=bigint) :: blocksize, iBig, iStart, iEnd
  integer :: outnPad, outipad(2), oldoutipad(2), oldoutnpad
  integer(bigint) :: nString, nString2
  real, allocatable :: float32(:)
  integer(kind=1), pointer :: iBytes(:) => null(), iHeader(:) => null()
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
!  sizeOfFloat = sizeof(float)
!  sizeOfInt = sizeof(int)
!  sizeOfInt1 = sizeof(int1)
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

  ! just return if not the first set of domain decomposed threads
  if (grid%splitOverMpi.and.(myHydroSetGlobal /= 0)) goto 666 
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
        call countSubcellsMPI(grid, nVoxels, nSubcellArray=nSubcellArray, includeGhosts=vtkIncludeGhosts)
        ncells = nSubcellArray(myRankGlobal)
     endif
     deallocate(nSubcellArray)
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


  if (Writeoutput) then


!     nBytesPoints = sizeofFloat * 3 * nPointsGlobal
!     nBytesConnect =  sizeofInt * nPointsGlobal
!     nBytesOffsets = sizeofint * nCellsGlobal
!     nBytesCellTypes = sizeofint1 * nCellsGlobal

     nBytesPoints = 4 * 3 * nPointsGlobal
     nBytesConnect = 4 * nPointsGlobal
     nBytesOffsets = 4 * nCellsGlobal
     nBytesCellTypes = 1 * nCellsGlobal
  endif

  if (writeheader) then
     open(lunit, file=vtkFilename, form="unformatted",access="stream",status="replace")
     buffer = '<?xml version="1.0"?>'//lf
     write(lunit) trim(buffer)
     buffer = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'//lf
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
     allocate(points(1:3,1:nPoints))
     call getPoints(grid, nPoints, points, vtkincludeGhosts)
     allocate(float32(1:(nPointsGlobal*3)))
     float32 = RESHAPE(points, (/SIZE(float32)/))
!     call base64encode(writeheader, pstring, nString, float32=float32)



     open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
!     write(lunit) pstring(1:nString)

     call convertandcompress(iBytes, iHeader, farray32=float32)
     call base64encode(.false., pstring, nString, iArray8=iHeader)
     call base64encode(.false., pstring2, nString2, iArray8=iBytes)
     write(lunit) pstring(1:nString), pstring2(1:nString2)
     deallocate(pString, pstring2)
     deallocate(iBytes, iHeader)
     deallocate(float32, points)
     close(lunit)
  else if (grid%splitOverMPI) then
#ifdef MPI
     call writePointsDecomposed(grid, vtkFilename, lunit, nPoints, vtkincludeGhosts)
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

     buffer = '        <DataArray type="Int64" Name="connectivity" format="binary"'// &
          ' >'//lf
     write(lunit) trim(buffer)

     blockSize = nPointsGlobal / 10
     do i = 1, 10

        iStart = (i-1) * blocksize + 1
        iEnd = iStart + blockSize - 1
        if (i == 10) iEnd = nPointsGlobal

        allocate(connectivity(iStart:iEnd))
        do iBig = iStart,iEnd
           connectivity(iBig) = iBig-1
        enddo
        if (i == 1) then
           call convertandcompressInSections(iBytes, iHeader, iarray64=connectivity, firstTime=.true.)
        else if (i == 10) then
           call convertandcompressInSections(iBytes, iHeader, iarray64=connectivity, lastTime=.true.)
        else
           call convertandcompressInSections(iBytes, iHeader, iarray64=connectivity)
        endif
        deallocate(connectivity)
     enddo


     call base64encode(.false., pstring, nString, iArray8=iHeader)
     write(lunit) pstring(1:nString)
     blockSize = size(iBytes,kind=bigint)/10
     do i = 1, 10
        iStart = (i-1) * blocksize + 1
        iEnd = iStart + blockSize - 1
        if (i == 10) iEnd = size(iBytes,kind=bigInt)
        if (i == 1) then
           call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                outnpad = outnpad, outipad=outipad, padEnd=.false.)
        else if (i == 10) then
           call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                inputnpad = oldoutnpad, inputipad = oldoutipad, &
                padend = .true.)
        else
           call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                inputnpad = oldoutnpad, inputipad = oldoutipad, &
                outnpad = outnpad, outipad = outipad, padend = .false.)
        endif
        oldoutnpad = outnpad
        oldoutipad = Outipad
        write(lunit) pstring2(1:nString2)
        deallocate(pstring2)
     end do

     deallocate(pString)
     deallocate(iBytes, iHeader)

  
     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)


     buffer = '        <DataArray type="Int64" Name="offsets" format="binary">'//lf
     write(lunit) trim(buffer)


     allocate(offsets(1:nCellsGlobal))
     do i = 1, nCellsGlobal
        offsets(i) = i * nPointOffset
     enddo



     call convertandcompress(iBytes, iHeader, iarray64=offsets)
     call base64encode(.false., pstring, nString, iArray8=iHeader)
     call base64encode(.false., pstring2, nString2, iArray8=iBytes)
     write(lunit) pstring(1:nString)
     write(lunit) pstring2(1:nString2)
     deallocate(pString, pstring2)
     deallocate(iBytes, iHeader, offsets)
 


    buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)


     buffer = '        <DataArray type="UInt8" Name="types" format="binary">'//lf
     write(lunit) trim(buffer)

     allocate(cellTypes(1:nCellsGlobal))

     if ((nPointOffset == 8).and.(.not.grid%octreeRoot%cylindrical)) cellTypes = 11
     if ((nPointOffset == 8).and.(grid%octreeRoot%cylindrical)) cellTypes = 12
     if (nPointOffset == 4) cellTypes = 8



     call convertandcompress(iBytes, iHeader, iarray8=cellTypes)
     call base64encode(.false., pstring, nString, iArray8=iHeader)
     call base64encode(.false., pstring2, nString2, iArray8=iBytes)
     write(lunit) pstring(1:nString), pstring2(1:nString2)
     deallocate(pString, pstring2)
     deallocate(iBytes, cellTypes)

     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)

     buffer = '      </Cells>'//lf
     write(lunit) trim(buffer)
     buffer = '      <CellData>'//lf
     write(lunit) trim(buffer)
     close(lunit)

  endif

  do ivalues = 1, nValueType
     select case (valueType(iValues))
     case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom","radforce","fvisc")
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
           call getValues(grid, valueType(iValues), rarray, includeGhosts=vtkIncludeGhosts)
           allocate(float32(1:nCellsGlobal))
           float32 = rarray(1,1:nCellsGlobal)
!           call base64encode(writeheader, pstring, nString, float32=float32)
           open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")
!           write(lunit) pstring(1:nString)

           call convertandcompress(iBytes, iHeader, farray32=float32)
           call base64encode(.false., pstring, nString, iArray8=iHeader)
           call base64encode(.false., pstring2, nString2, iArray8=iBytes)
           write(lunit) pstring(1:nString), pstring2(1:nString2)
           deallocate(pString, pstring2)
           deallocate(iBytes)


           close(lunit)
!           deallocate(pString)
           deallocate(float32)
           deallocate(rArray)
        else if (grid%splitOverMPI) then
#ifdef MPI
           call writeDomainDecomposed(valueType(iValues), grid, vtkFilename, lunit, nCells)
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
           call getValues(grid, valueType(iValues), rarray, includeGhosts=vtkIncludeGhosts)
           allocate(float32(1:nCellsGlobal*3))
           float32 = RESHAPE(rarray, (/SIZE(float32)/))
!           call base64encode(writeheader,pstring, nString, float32=float32)
           open(lunit, file=vtkFilename, form="unformatted",access="stream",status="old",position="append")

           call convertandcompress(iBytes, iHeader, farray32=float32)
           call base64encode(.false., pstring, nString, iArray8=iHeader)
           call base64encode(.false., pstring2, nString2, iArray8=iBytes)
           write(lunit) pstring(1:nString), pstring2(1:nString2)
           deallocate(pString, pstring2)
           deallocate(iBytes)

!           write(lunit) pstring(1:nString)
           close(lunit)
!           deallocate(pString)
           deallocate(float32)
           deallocate(rArray)
        else if (grid%splitOverMPI) then
#ifdef MPI
           call writeDomainDecomposed(valueType(iValues), grid, vtkFilename, lunit, nCells)
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
     close(lunit)
  endif

#ifdef MPI
666 continue
#endif
  if (associated(iHeader)) deallocate(iHeader)
end subroutine writeXMLVtkFileAMR


  subroutine getValues(grid, valuetype, rarray, includeGhosts)
    type(GRIDTYPE) :: grid
    character(len=*) :: valueType
    real :: rArray(:,:)
    integer :: n
    logical :: includeGhosts

    n = 0 

    call getValuesRecursive(grid, grid%octreeRoot, valueType, rArray, n, includeGhosts)

  end subroutine getValues

    recursive subroutine getValuesRecursive(grid, thisOctal, valueType, rArray, n, includeGhosts)
      use inputs_mod, only : lambdaSmooth, hydrodynamics
      type(OCTAL), pointer :: thisOctal, child
      logical :: includeGhosts
      type(GRIDTYPE) :: grid
      character(len=*) :: valueType
      real :: rArray(:,:)
      integer :: i, subcell, n, iVal, nVal
      integer, save ::  iLambda
      real(double) :: ksca, kabs, value, v
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
                  call getValuesRecursive(grid, child, valueType, rArray, n, includeGhosts)
                  exit
               end if
            end do
         else

#ifdef MPI
            if ( (.not.octalOnThread(thisOctal, subcell, myRankGlobal)) .and. grid%splitOverMPI ) cycle
            if (hydrodynamics) then
               if (checkGhost(thisOctal, subcell).and.(.not.includeGhosts)) cycle
            endif
#endif



            nVal = 1
            if (thisOctal%threed.and.thisOctal%cylindrical) then
               nVal = max(1, nint(returndPhi(thisOctal)/(10.d0*degtorad)))
            endif

            do iVal = 1, nVal
            n = n + 1

            select case (valueType)
               case("rho")

                  rArray(1, n) = real(returnPhysicalUnitDensity(thisOctal%rho(subcell)))

               case("J=0")
                  rArray(1, n) = real(thisOctal%molecularlevel(1,subcell))

               case("J=1")
                  rArray(1, n) = real(thisOctal%molecularlevel(2,subcell))

               case("J=2")
                  rArray(1, n) = real(thisOctal%molecularlevel(3,subcell))

               case("J=3")
                  rArray(1, n) = real(thisOctal%molecularlevel(4,subcell))

               case("J=4")
                  rArray(1, n) = real(thisOctal%molecularlevel(5,subcell))

               case("J=5")
                  rArray(1, n) = real(thisOctal%molecularlevel(6,subcell))

               case("J=10")
                  rArray(1, n) = real(thisOctal%molecularlevel(11,subcell))

               case("J=16")
                  rArray(1, n) = real(thisOctal%molecularlevel(17,subcell))

               case("dI")
                  rArray(1, n) = real(thisOctal%newmolecularlevel(1,subcell))

               case("dIattenuated")
                  rArray(1, n) = real(thisOctal%newmolecularlevel(2,subcell))

               case("i0")
                  rArray(1, n) = real(thisOctal%newmolecularlevel(3,subcell))

! galLon and galLat re-use storage used for i0 and dIattenuated 
               case("galLon")
                  rArray(1, n) = real(thisOctal%newmolecularlevel(2,subcell))

               case("galLat")
                  rArray(1, n) = real(thisOctal%newmolecularlevel(3,subcell))

               case("crossing")
                  rArray(1, n) = real(thisOctal%newmolecularlevel(4,subcell))

               case("tauacrosscell")
                  rArray(1, n) = real(thisOctal%newmolecularlevel(5,subcell))

               case("tau10")
                  rArray(1, n) = real(thisOctal%tau(1,subcell))

               case("tau21")
                  rArray(1, n) = real(thisOctal%tau(2,subcell))

               case("tau32")
                  rArray(1, n) = real(thisOctal%tau(3,subcell))

               case("tau43")
                  rArray(1, n) = real(thisOctal%tau(4,subcell))

               case("tau54")
                  rArray(1, n) = real(thisOctal%tau(5,subcell))

               case("tau65")
                  rArray(1, n) = real(thisOctal%tau(6,subcell))

               case("tau76")
                  rArray(1, n) = real(thisOctal%tau(7,subcell))

               case("tau87")
                  rArray(1, n) = real(thisOctal%tau(8,subcell))

               case("level0error")
                  rArray(1, n) = real(10**((real(thisOctal%levelconvergence(1,subcell)) / 6553.6) - 4))

               case("level1error")
                  rArray(1, n) = real(10**((real(thisOctal%levelconvergence(2,subcell)) / 6553.6) - 4))

               case("level2error")
                  rArray(1, n) = real(10**((real(thisOctal%levelconvergence(3,subcell)) / 6553.6) - 4))

               case("level3error")
                  rArray(1, n) = real(10**((real(thisOctal%levelconvergence(4,subcell)) / 6553.6) - 4))

               case("level4error")
                  rArray(1, n) = real(10**((real(thisOctal%levelconvergence(5,subcell)) / 6553.6) - 4))

               case("niter")
                  rArray(1, n) = real(int(thisOctal%convergence(subcell)/100))

               case("nh2")
                  rArray(1, n) = real(real(thisOctal%nh2(subcell)))

               case("convergence")
                  rArray(1, n) = real(mod(thisOctal%convergence(subcell),1.0))

               case("adot")
                  rArray(1, n) = real(real(thisOctal%adot(subcell)))

               case("slowestlevel")
                  rArray(1, n) = real(floor(mod(thisoctal%convergence(subcell),100.0)))

               case("molabundance")
                  rArray(1, n) = real(thisOctal%molabundance(subcell))

               case("bnu")
                  rArray(1, n) = real(real(thisOctal%bnu(1,subcell)))

               case("dc0")
                  rArray(1, n) = real(thisOctal%molecularlevel(1,subcell) * thisOctal%departcoeff(1,subcell))

               case("dc1")
                  rArray(1, n) = real(thisOctal%molecularlevel(2,subcell) * thisOctal%departcoeff(2,subcell))

               case("dc2")
                  rArray(1, n) = real(thisOctal%molecularlevel(3,subcell) * thisOctal%departcoeff(3,subcell))

               case("dc3")
                  rArray(1, n) = real(thisOctal%molecularlevel(4,subcell) * thisOctal%departcoeff(4,subcell))

               case("dc4")
                  rArray(1, n) = real(thisOctal%molecularlevel(5,subcell) * thisOctal%departcoeff(5,subcell))

               case("jnu10")
                  rArray(1, n) = real(thisOctal%jnu(1,subcell))

               case("jnu")
                  rArray(1, n) = real(thisOctal%biasLine3d(subcell))

               case("biasline")
                  rArray(1, n) = real(thisOctal%biasLine3d(subcell))

               case("dust1")
                  rArray(1, n) = real(real(thisOctal%dustTypeFraction(subcell,1)))

               case("dust2")
                  rArray(1, n) = real(real(thisOctal%dustTypeFraction(subcell,2)))

               case("dust3")
                  rArray(1, n) = real(real(thisOctal%dustTypeFraction(subcell,3)))

               case("dust4")
                  rArray(1, n) = real(real(thisOctal%dustTypeFraction(subcell,4)))

               case("dust5")
                  rArray(1, n) = real(real(thisOctal%dustTypeFraction(subcell,5)))


               case("q11")
                  rArray(1, n) = real(thisOctal%qViscosity(subcell,1,1))
               case("q22")
                  rArray(1, n) = real(thisOctal%qViscosity(subcell,2,2))
               case("q33")
                  rArray(1, n) = real(thisOctal%qViscosity(subcell,3,3))

               case("bias")
                  rArray(1, n) = real(real(thisOctal%biasCont3d(subcell)))

               case("mpithread")
                  rArray(1, n) = real(real(thisOctal%mpithread(subcell)))

               case("bcond")
                  rArray(1, n) = real(real(thisOctal%boundaryCondition(subcell)))


               case("deltaT")
                  rArray(1, n) = real(real(thisOctal%temperature(subcell)-thisOctal%oldtemperature(subcell)))


               case("scattered")
                  rArray(1, n) = real(real(thisOctal%scatteredIntensity(subcell,5,3)))

               case("hydrovelocity")
                  if (thisOctal%threeD) then
                     rArray(1, n) = real(real(returnPhysicalUnitSpeed(thisOctal%rhou(subcell)/thisOctal%rho(subcell))/1.e5))
                     rArray(2, n) = real(returnPhysicalUnitSpeed(thisOctal%rhov(subcell)/thisOctal%rho(subcell))/1.e5)
                     rArray(3, n) = real(returnPhysicalUnitSpeed(thisOctal%rhow(subcell)/thisOctal%rho(subcell))/1.e5)
                  else
                     rArray(1, n) = real(real(returnPhysicalUnitSpeed(thisOctal%rhou(subcell)/thisOctal%rho(subcell))/1.e5))
                     rArray(2, n) = real(returnPhysicalUnitSpeed(thisOctal%rhow(subcell)/thisOctal%rho(subcell))/1.e5)
                     rArray(3, n) = 0.
                  endif

               case("radmom")
                     rArray(1, n) = real(thisOctal%radiationMomentum(subcell)%x)
                     rArray(2, n) = real(thisOctal%radiationMomentum(subcell)%y)
                     rArray(3, n) = real(thisOctal%radiationMomentum(subcell)%z)

               case("radforce")
                     v = cellVolume(thisOctal, subcell)*1.d30
                     if (thisOctal%threed) then
                        rArray(1, n) = real(thisOctal%kappaTimesFlux(subcell)%x/cSpeed)
                        rArray(2, n) = real(thisOctal%kappaTimesFlux(subcell)%y/cSpeed)
                        rArray(3, n) = real(thisOctal%kappaTimesFlux(subcell)%z/cSpeed)
                     else
                        rArray(1, n) = real(thisOctal%kappaTimesFlux(subcell)%x/cSpeed)
                        rArray(2, n) = real(thisOctal%kappaTimesFlux(subcell)%z/cSpeed)
                        rArray(3, n) = 0.
                     endif

               case("fvisc")
                        rArray(1, n) = real(thisOctal%fViscosity(subcell)%x)
                        rArray(2, n) = real(thisOctal%fViscosity(subcell)%z)
                        rArray(3, n) = 0.


               case("velocity")
                     ! stop vectors from showing up in visit if too big
                     if(thisoctal%velocity(subcell)%x .ge. 1.d0) then
                        thisoctal%velocity(subcell)%x = 0.d0
                        thisoctal%velocity(subcell)%y = 0.d0
                        thisoctal%velocity(subcell)%z = 0.d0
                     endif
                     if (thisOctal%threed) then
                        rArray(1, n) = real(real(thisOctal%velocity(subcell)%x*cspeed/1.e5))
                        rArray(2, n) = real(thisOctal%velocity(subcell)%y*cspeed/1.e5)
                        rArray(3, n) = real(thisOctal%velocity(subcell)%z*cspeed/1.e5)
                     else
                        rArray(1, n) = real(real(thisOctal%velocity(subcell)%x*cspeed/1.e5))
                        rArray(2, n) = real(thisOctal%velocity(subcell)%z*cspeed/1.e5)
                        rArray(3, n) = real(thisOctal%velocity(subcell)%y*Cspeed/1.e5)
                     endif
              case("cornervel")
                 rVec = subcellCentre(thisOctal, subcell)
                 vel = amrGridVelocity(grid%octreeRoot,rvec,startOctal=thisOctal,&
                      actualSubcell=subcell) 
                 rArray(1, n) = real(real(vel%x*cspeed/1.e5))
                 rArray(2, n) = real(vel%y*cspeed/1.e5)
                 rArray(3, n) = real(vel%z*cspeed/1.e5)

               case("ne")
                  rArray(1, n) = real(real(thisOctal%ne(subcell)))

               case("pressure")
                  rArray(1, n) = real(real(thisOctal%pressure_i(subcell)))

               case("inflow")
                  if (thisOctal%inflow(subcell)) then
                     rArray(1, n) = real(1.)
                  else
                     rArray(1, n) = real(0.)
                  endif

               case("haschild")
                  if (thisOctal%haschild(subcell)) then
                     rArray(1, n) = real(1.)
                  else
                     rArray(1, n) = real(0.)
                  endif

               case("NH2")
                  rArray(1, n) = real(real( thisOctal%NH2(subcell) ))

               case("fixedtemp")
                  if (thisOctal%fixedTemperature(subcell)) then
                     rArray(1, n) = real(1.)
                  else
                     rArray(1, n) = real(0.)
                  endif

#ifdef PHOTOION
               case("HI")
                  rArray(1, n) = real(real(thisOctal%ionfrac(subcell,returnIonNumber("H I", grid%ion, grid%nIon))))

               case("HeI")
                  rArray(1, n) = real(real(thisOctal%ionfrac(subcell,returnIonNumber("He I", grid%ion, grid%nIon))))

               case("HeII")
                  rArray(1, n) = real(real(thisOctal%ionfrac(subcell,returnIonNumber("He II", grid%ion, grid%nIon))))

               case("OI")
                  rArray(1, n) = real(real(thisOctal%ionfrac(subcell,returnIonNumber("O I", grid%ion, grid%nIon))))

               case("OII")
                  rArray(1, n) = real(real(thisOctal%ionfrac(subcell,returnIonNumber("O II", grid%ion, grid%nIon))))

               case("OIII")
                  rArray(1, n) = real(real(thisOctal%ionfrac(subcell,returnIonNumber("O III", grid%ion, grid%nIon))))
#endif

               case("sourceCont")
                  rArray(1, n) = real(real(thisOctal%normSourceContribution(subcell, 1)))

               case("logRho")
                  rArray(1, n) = real(log10(returnPhysicalUnitDensity(thisOctal%rho(subcell))))

               case("temperature")
                  rArray(1, n) = real(real(thisOctal%temperature(subcell)))

               case("chiline")
                  if (.not.associated(thisOctal%chiline)) then
                     rArray(1, n) = 0.d0
                  else
                     rArray(1, n) = real(max( real(thisOctal%chiline(subcell)), min_single_prec ))
                  endif


               case("microturb")
                  rArray(1, n) = real(max( real(thisOctal%microturb(subcell)*cspeed/1.e5), min_single_prec ))

               case("etaline")
                  rArray(1, n) = real(max ( real(thisOctal%etaline(subcell)), min_single_prec ))

               case("n4")
                  rArray(1, n) = real(max ( real(thisOctal%atomLevel(subcell,1,4)), min_single_prec ))

               case("sourceline")
                  rArray(1, n) = real(max ( real(thisOctal%etaline(subcell)), min_single_prec )/ &
                       max( real(thisOctal%chiline(subcell)), min_single_prec ))

               case("tau")
                  if (firstTime) then
                     call locate(grid%lamArray, grid%nLambda, lambdaSmooth, ilambda)
                     firstTime = .false.
                  endif
                  call returnKappa(grid, thisOctal, subcell, ilambda=ilambda,&
                       kappaSca=ksca, kappaAbs=kabs)
                  value = thisOctal%subcellSize * (ksca + kabs)
                  rArray(1, n) = real(real(value))

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

               case("diff")
                  if (thisOctal%diffusionApprox(subcell)) then
                     value = 1.d0
                  else
                     value = 0.d0
                  endif
                  rArray(1, n) = real(value)

               case("mpistore")
                  if (associated(thisOctal%mpiBoundaryStorage)) then
                     value = 1.d0
                  else
                     value = 0.d0
                  endif
                  rArray(1, n) = real(value)
                  
               case("phi")
                  rArray(1, n) = real(thisOctal%phi_i(subcell))

               case("phigas")
                  rArray(1, n) = real(thisOctal%phi_gas(subcell))

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
  subroutine writeDomainDecomposed(valueType, grid, vtkFilename, lunit, nCells)
    use mpi
    use inputs_mod, only : vtkIncludeGhosts
    character(len=*) :: valueType, vtkFilename
    type(GRIDTYPE) :: grid
    integer :: nCells
    integer :: lunit
    integer :: iThread
    integer :: ierr
    integer(bigint) :: nString, nString2
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 51
    real, pointer :: rArray(:,:)
    real, allocatable :: float32(:)
    character, pointer :: pstring(:)
    character, pointer :: pstring2(:)
    integer(kind=1), pointer :: iBytes(:) => null(), iHeader(:) => null()
    logical :: scalarvalue, vectorvalue
    integer :: j, k
    integer(kind=bigint) :: blockSize, i, istart, iend
    integer :: outnpad, outipad(2), oldoutnpad, oldoutipad(2)



     select case (valueType)
     case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom","radforce","fvisc")
        scalarvalue = .false.
        vectorvalue = .true.
     case DEFAULT
        scalarvalue = .true.
        vectorvalue = .false.
     end select

     
     if (myrankGlobal == 1) then
        if (scalarValue) then
           allocate(float32(1:nCells))
           allocate(rArray(1:1, 1:nCells))
           call getValues(grid, valueType, rarray(1:1,1:nCells), includeGhosts=vtkIncludeGhosts)
           float32(1:nCells) = RESHAPE(rArray(1:1, 1:nCells),(/ncells/))
        else
           allocate(float32(1:nCells*3))
           allocate(rArray(1:3, 1:nCells))
           call getValues(grid, valueType, rarray, includeGhosts=vtkIncludeGhosts)
           float32(1:nCells*3) = RESHAPE(rarray, (/SIZE(float32(1:ncells*3))/))
        endif
        call convertandcompressInSections(iBytes, iHeader, farray32=float32, firstTime=.true.)
        deallocate(float32, rArray)
        do iThread = 2, nHydroThreadsGlobal
           call MPI_RECV(k, 1, MPI_INTEGER, iThread, tag, MPI_COMM_WORLD, status, ierr)
           allocate(float32(1:k))
           call MPI_RECV(float32, k, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
           if (iThread /= nHydroThreadsglobal) then
              call convertandcompressInSections(iBytes, iHeader, farray32=float32)
           else
              call convertandcompressInSections(iBytes, iHeader, farray32=float32, lastTime=.true.)
           endif
           deallocate(float32)
        enddo

        call base64encode(.false., pstring, nString, iArray8=iHeader)
        open(lunit, file=vtkFilename, form="unformatted", status="old", access="stream", position="append")
        write(lunit) pstring(1:nString)
        blockSize = size(iBytes,kind=bigInt) / 10
        do i = 1, 10
           iStart = (i-1) * blocksize + 1
           iEnd = iStart + blockSize - 1
           if (i == 10) iEnd = size(iBytes,kind=bigInt)
           if (i == 1) then
              call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                   outnpad = outnpad, outipad=outipad, padEnd=.false.)
           else if (i == 10) then
              call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                   inputnpad = oldoutnpad, inputipad = oldoutipad, &
                   padend = .true.)
           else
              call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                   inputnpad = oldoutnpad, inputipad = oldoutipad, &
                   outnpad = outnpad, outipad = outipad, padend = .false.)
           endif
           oldoutnpad = outnpad
           oldoutipad = Outipad
           write(lunit) pstring2(1:nString2)
           deallocate(pstring2)
        end do
 

        close(lunit)


        deallocate(pString)
        deallocate(iBytes, iHeader)

     else
     
        if (scalarValue) then
           allocate(float32(1:nCells))
           allocate(rArray(1:1, 1:nCells))
           call getValues(grid, valueType, rarray(1:1,1:nCells), includeGhosts=vtkIncludeghosts)
           float32(1:nCells) = reshape(rArray(1:1, 1:nCells),(/ncells/))
           j = nCells
        else
           allocate(float32(1:nCells*3))
           allocate(rArray(1:3, 1:nCells))
           call getValues(grid, valueType, rarray, vtkIncludeGhosts)
           float32 = RESHAPE(rarray, (/SIZE(float32)/))
           j = ncells*3
        endif
        call MPI_SEND(j, 1, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, ierr)
        call MPI_SEND(float32, j, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
        deallocate(float32, rArray)
     endif
        
     call MPI_BARRIER(amrCommunicator, ierr)
   end subroutine writeDomainDecomposed

  subroutine writePointsDecomposed(grid, vtkFilename, lunit, nPoints, includeGhosts)
    use mpi
    use timing
    character(len=*) ::vtkFilename
    logical :: includeGhosts
    type(GRIDTYPE) :: grid
    integer(bigint) :: nPoints
    integer :: lunit
    integer :: iThread
    integer :: ierr
    integer(bigint) :: nString, nstring2
    integer :: status(MPI_STATUS_SIZE)
    integer, parameter :: tag = 55
    real, pointer :: rArray(:,:)
    real, allocatable :: float32(:)
    character, pointer :: pstring(:)
    character, pointer :: pstring2(:)
    integer(kind=1), pointer :: iBytes(:) => null(), iHeader(:) => null()
    integer(bigint) :: j, k
    integer :: outnpad, oldoutnpad
    integer :: outipad(2),  oldoutipad(2)
    integer(kind=bigint) :: blocksize, iStart, iEnd, i

     if (myrankGlobal == 1) then

        allocate(float32(1:nPoints*3))
        allocate(rArray(1:3, 1:nPoints))
        call getPoints(grid, nPoints, rarray, includeGhosts)
        float32(1:nPoints*3) = RESHAPE(rarray(1:3,1:nPoints), (/SIZE(float32(1:nPoints*3))/))
        call convertandcompressInSections(iBytes, iHeader, farray32=float32,firstTime=.true.)
        deallocate(float32, rArray)

        do iThread = 2, nHydroThreadsGlobal
           call MPI_RECV(k, 1, MPI_INTEGER8, iThread, tag, MPI_COMM_WORLD, status, ierr)
           allocate(float32(1:k))
           call MPI_RECV(float32, k, MPI_REAL, iThread, tag, MPI_COMM_WORLD, status, ierr)
           if (iThread /= nHydroThreadsGlobal) then
              call convertandcompressInSections(iBytes, iHeader, farray32=float32)
           else
              call convertandcompressInSections(iBytes, iHeader, farray32=float32, lastTime = .true.)
           endif
           deallocate(float32)
        enddo

        call base64encode(.false., pstring, nString, iArray8=iHeader)
        open(lunit, file=vtkFilename, form="unformatted", status="old", access="stream", position="append")
        write(lunit) pstring(1:nString)
        blockSize = size(iBytes,kind=bigInt) / 10
        do i = 1, 10
           iStart = (i-1) * blocksize + 1
           iEnd = iStart + blockSize - 1
           if (i == 10) iEnd = size(iBytes,kind=bigInt)
           if (i == 1) then
              call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                   outnpad = outnpad, outipad=outipad, padEnd=.false.)
           else if (i == 10) then
              call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                   inputnpad = oldoutnpad, inputipad = oldoutipad, &
                   padend = .true.)
           else
              call base64encode(.false., pstring2, nString2, iArray8=iBytes(iStart:iEnd), &
                   inputnpad = oldoutnpad, inputipad = oldoutipad, &
                   outnpad = outnpad, outipad = outipad, padend = .false.)
           endif
           oldoutnpad = outnpad
           oldoutipad = Outipad
           write(lunit) pstring2(1:nString2)
           deallocate(pstring2)
        end do
 

        close(lunit)


        deallocate(pString)
        deallocate(iBytes, iHeader)

     else
     
        allocate(float32(1:int(nPoints,kind=bigint)*3))
        allocate(rArray(1:3, 1:int(nPoints,kind=bigint)))

        call getPoints(grid, nPoints, rarray, includeGhosts)
        float32 = RESHAPE(rarray, (/SIZE(float32)/))
        j = int(nPoints,kind=bigint)*3
        call MPI_SEND(j, 1, MPI_INTEGER8, 1, tag, MPI_COMM_WORLD, ierr)
        call MPI_SEND(float32, j, MPI_REAL, 1, tag, MPI_COMM_WORLD, ierr)
        deallocate(float32, rArray)
     endif
     

     call MPI_BARRIER(amrCommunicator, ierr)
   end subroutine writePointsDecomposed
#endif


subroutine writeParallelXMLVtkFileAMR(grid, vtkFilename, valueTypeFilename, valueTypeString)
#ifdef MPI
  use mpi
#endif
  use inputs_mod, only : vtkIncludeGhosts
  type(GRIDTYPE) :: grid
  character(len=*) :: vtkFilename
  integer :: nValueType
  character(len=20) :: valueType(50)
  character(len=*), optional ::  valueTypeFilename
  character(len=*), optional ::  valueTypeString(:)
  integer :: nCells
  integer(bigint) :: nPointsGlobal, nPoints
  integer :: lunit = 69
  integer :: nOctals, nVoxels
  integer(bigint) :: i
  integer :: nPointOffset
  integer(kind=1), allocatable :: cellTypes(:)
  integer(kind=bigint), allocatable :: offsets(:)
  integer(kind=bigint), allocatable :: connectivity(:)
  real, pointer :: points(:,:)
  integer :: nCellsGlobal
  integer(kind=bigint) :: nBytesPoints, nBytesCellTypes, nBytesConnect, nBytesOffsets
  real, pointer :: rArray(:,:)
  real :: float
  integer :: int,  iValues
  integer(kind=1) :: int1
  Character(len=200) :: buffer
  character(len=1) :: lf
  character, pointer :: pstring(:), pstring2(:)
  character(len=12) :: str1, str2
  logical :: vectorValue, scalarValue
  integer :: sizeOfFloat, sizeOfInt, sizeOfInt1
  integer(bigint) :: nString, nString2, nUniquePoints
  integer :: j
  character(len=80) :: pVtkfilename
  real, allocatable :: float32(:)
  integer(kind=1), pointer :: iBytes(:) => null(), iHeader(:) => null()
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

!  sizeofFloat = sizeof(float)
!  sizeofInt = sizeof(int)
!  sizeofInt1 = sizeof(int1)

#endif

#ifdef MPI
  ! just return if the grid is decomposed and MPI job and this is rank zero thread
  if (grid%splitOverMpi.and.(myRankGlobal == 0)) goto 666

  ! just return if not the first set of domain decomposed threads
  if (grid%splitOverMpi.and.(myHydroSetGlobal /= 0)) goto 666 

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
        call countSubcellsMPI(grid, nVoxels, nSubcellArray=nSubcellArray, includeGhosts=vtkIncludeGhosts)
        ncells = nSubcellArray(myRankGlobal)
        deallocate(nSubcellArray)
     endif
  endif
#endif
  nCellsGlobal = nVoxels
  nPointsGlobal = nCellsGlobal * nPointOffset
  nPoints = ncells * nPointOffset

  writeHeader = .true.

  if (myrankGlobal == 1) then
     j  = index(vtkFilename, ".") -1
     open(lunit, file=vtkFilename(1:j)//".pvtu", status="unknown", access="stream")
     write(lunit) '<?xml version="1.0"?>'//lf
     write(lunit) '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
     write(lunit) '  <PUnstructuredGrid GhostLevel="0">'//lf
     write(lunit) '    <PCellData Scalars="scalars">'//lf
     do i = 1, nValueType
        select case (valueType(i))
        case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom","radforce","fvisc")
           scalarvalue = .false.
           vectorvalue = .true.
        case DEFAULT
           scalarvalue = .true.
           vectorvalue = .false.
        end select
        
        if (scalarvalue) then
           write(lunit) '      <PDataArray type="Float32" Name="',trim(valueType(i)),'"/>'//lf
        else
           write(lunit) '      <PDataArray type="Float32" Name="',trim(valueType(i)),'" NumberOfComponents="3"/>'//lf
        endif
     enddo
     write(lunit) '    </PCellData>'//lf

     write(lunit) '  <PPoints>'//lf
     write(lunit) '    <PDataArray type="Float32" NumberOfComponents="3"/>'//lf
     write(lunit) '  </PPoints>'//lf

     do i = 1, nHydroThreadsGlobal
        j  = index(vtkFilename, ".") -1
        write(pVtkFilename, '(a,a,i3.3,a)') vtkFilename(1:j),"_",i,".vtu"
        write(lunit) '  <Piece Source="'//trim(pVtkFilename)//'"/>'//lf
     enddo
     write(lunit) '  </PUnstructuredGrid>'//lf
     write(lunit) '</VTKFile>'//lf
     close(lunit)
  endif


  j  = index(vtkFilename, ".") -1
  write(pVtkFilename, '(a,a,i3.3,a)') vtkFilename(1:j),"_",myrankGlobal,".vtu"
     

  if (Writeoutput) then
!     nBytesPoints = sizeofFloat * 3 * nPointsGlobal
!     nBytesConnect =  sizeofInt * nPointsGlobal
!     nBytesOffsets = sizeofint * nCellsGlobal
!     nBytesCellTypes = sizeofint1 * nCellsGlobal

     nBytesPoints = 4 * 3 * nPointsGlobal
     nBytesConnect =  4 * nPointsGlobal
     nBytesOffsets = 4 * nCellsGlobal
     nBytesCellTypes = 1 * nCellsGlobal
  endif

  allocate(points(1:3, 1:nPoints))
  allocate(connectivity(1:nPoints))

  call getUniquePoints(grid, nPoints, points, connectivity, nUniquePoints, includeGhosts=vtkIncludeGhosts)


  if (writeheader) then
     open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="replace")
     buffer = '<?xml version="1.0"?>'//lf
     write(lunit) trim(buffer)
     buffer = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'//lf
     write(lunit) trim(buffer)
     buffer = '  <UnstructuredGrid>'//lf 
     write(lunit) trim(buffer)
     write(str1(1:12),'(i10)') nUniquePoints
     write(str2(1:12),'(i10)') nCells
     buffer = '    <Piece NumberOfPoints="'//str1//'" NumberOfCells="'//str2//'">'//lf
     write(lunit) trim(buffer)
     buffer = '      <Points>'//lf
     write(lunit) trim(buffer)
     buffer = '       <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" format="binary">'//lf
     write(lunit) trim(buffer)
     close(lunit)
  endif


  if (writeheader) then 
     allocate(float32(1:(nUniquePoints*3)))
     float32 = RESHAPE(points(1:3,1:nUniquePoints), (/SIZE(float32)/))
     open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")

     call convertandcompress(iBytes, iHeader, farray32=float32)
     call base64encode(.false., pstring, nString, iArray8=iHeader)
     call base64encode(.false., pstring2, nString2, iArray8=iBytes)
     write(lunit) pstring(1:nString), pstring2(1:nString2)
     deallocate(pString, pstring2)
     deallocate(iBytes, iHeader)
     deallocate(float32, points)
     close(lunit)
  endif

  if (writeheader) then
     open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")
     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)
     buffer = '      </Points>'//lf
     write(lunit) trim(buffer)
     buffer = '      <Cells>'//lf
     write(lunit) trim(buffer)

     buffer = '        <DataArray type="Int64" Name="connectivity" format="binary"'// &
          ' >'//lf
     write(lunit) trim(buffer)

     call convertandcompress(iBytes, iHeader, iArray64=connectivity)
     call base64encode(.false., pstring, nString, iArray8=iHeader)
     call base64encode(.false., pstring2, nString2, iArray8=iBytes)
     write(lunit) pstring(1:nString), pstring2(1:nString2)
     deallocate(pString, pstring2)
     deallocate(iBytes, iHeader)
     deallocate(connectivity)


     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)


     buffer = '        <DataArray type="Int64" Name="offsets" format="binary">'//lf
     write(lunit) trim(buffer)


     allocate(offsets(1:nCells))
     do i = 1, nCells
        offsets(i) = i * nPointOffset
     enddo



     call convertandcompress(iBytes, iHeader, iarray64=offsets)
     call base64encode(.false., pstring, nString, iArray8=iHeader)
     call base64encode(.false., pstring2, nString2, iArray8=iBytes)
     write(lunit) pstring(1:nString)
     write(lunit) pstring2(1:nString2)
     deallocate(pString, pstring2)
     deallocate(iBytes, iHeader, offsets)
 




    buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)


     buffer = '        <DataArray type="UInt8" Name="types" format="binary">'//lf
     write(lunit) trim(buffer)

     allocate(cellTypes(1:nCells))

     if ((nPointOffset == 8).and.(.not.grid%octreeRoot%cylindrical)) cellTypes = 11
     if ((nPointOffset == 8).and.(grid%octreeRoot%cylindrical)) cellTypes = 12
     if (nPointOffset == 4) cellTypes = 8



     call convertandcompress(iBytes, iHeader, iarray8=cellTypes)
     call base64encode(.false., pstring, nString, iArray8=iHeader)
     call base64encode(.false., pstring2, nString2, iArray8=iBytes)
     write(lunit) pstring(1:nString), pstring2(1:nString2)
     deallocate(pString, pstring2)
     deallocate(iBytes, cellTypes)

     buffer = lf//'        </DataArray>'//lf
     write(lunit) trim(buffer)

     buffer = '      </Cells>'//lf
     write(lunit) trim(buffer)
     buffer = '      <CellData>'//lf
     write(lunit) trim(buffer)
     close(lunit)

  endif

  do ivalues = 1, nValueType
     select case (valueType(iValues))
     case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel","radmom","radforce","fvisc")
        scalarvalue = .false.
        vectorvalue = .true.
     case DEFAULT
        scalarvalue = .true.
        vectorvalue = .false.
     end select

     if (scalarvalue) then
        if (writeheader) then
           open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")
           buffer = '        <DataArray type="Float32" Name="'//trim(valueType(iValues))//&
                '" format="binary">'//lf
           write(lunit) trim(buffer)
           close(lunit)
        endif

        if (writeheader) then
           allocate(rArray(1:1, 1:nCells))
           call getValues(grid, valueType(iValues), rarray, vtkIncludeGhosts)
           allocate(float32(1:nCells))
           float32 = rarray(1,1:nCells)
!           call base64encode(writeheader, pstring, nString, float32=float32)
           open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")
!           write(lunit) pstring(1:nString)

           call convertandcompress(iBytes, iHeader, farray32=float32)
           call base64encode(.false., pstring, nString, iArray8=iHeader)
           call base64encode(.false., pstring2, nString2, iArray8=iBytes)
           write(lunit) pstring(1:nString), pstring2(1:nString2)
           deallocate(pString, pstring2)
           deallocate(iBytes)


           close(lunit)
!           deallocate(pString)
           deallocate(float32)
           deallocate(rArray)
        endif

        if (writeheader) then
           open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")
           buffer = lf//'        </DataArray>'//lf
           write(lunit) trim(buffer)
           close(lunit)
        endif

     else
        if (writeheader) then
           open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")
           buffer = '        <DataArray type="Float32" NumberOfComponents="3" Name="'//trim(valueType(iValues))//&
                '" format="binary">'//lf
           write(lunit) trim(buffer)
           close(lunit)
        endif

        if (writeheader) then
           allocate(rArray(1:3, 1:nCells))
           call getValues(grid, valueType(iValues), rarray, includeGhosts=vtkIncludeGhosts)
           allocate(float32(1:nCells*3))
           float32 = RESHAPE(rarray, (/SIZE(float32)/))
!           call base64encode(writeheader,pstring, nString, float32=float32)
           open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")

           call convertandcompress(iBytes, iHeader, farray32=float32)
           call base64encode(.false., pstring, nString, iArray8=iHeader)
           call base64encode(.false., pstring2, nString2, iArray8=iBytes)
           write(lunit) pstring(1:nString), pstring2(1:nString2)
           deallocate(pString, pstring2)
           deallocate(iBytes)

!           write(lunit) pstring(1:nString)
           close(lunit)
!           deallocate(pString)
           deallocate(float32)
           deallocate(rArray)
        endif
           


        if (writeheader) then
           open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")
           buffer = lf//'        </DataArray>'//lf
           write(lunit) trim(buffer)
           close(lunit)
        endif
     endif
  enddo
  if (writeheader) then
     open(lunit, file=pvtkFilename, form="unformatted",access="stream",status="old",position="append")
     buffer = '     </CellData>'//lf
     write(lunit) trim(buffer)
     buffer = '    </Piece>'//lf
     write(lunit) trim(buffer)
     buffer = '  </UnstructuredGrid>'//lf
     write(lunit) trim(buffer)
     buffer = '</VTKFile>'//lf
     write(lunit) trim(buffer)
     close(lunit)
  endif

#ifdef MPI
666 continue
#endif
  if (associated(iHeader)) deallocate(iHeader)
end subroutine writeParallelXMLVtkFileAMR


  subroutine getUniquePoints(grid, nPoints, points, connectivity, nUniquePoints,includeGhosts)
    type(GRIDTYPE) :: grid
    real :: points(:,:)
    integer(bigint) :: connectivity(:)
    integer(bigint) :: nPoints, nUniquePoints
    logical :: includeGhosts

    nPoints = 0
    nUniquePoints = 0

    call recursiveGetUniquePoints(grid, grid%octreeRoot, nPoints, Points, connectivity, nUniquePoints, includeGhosts)

  contains

    recursive subroutine recursiveGetUniquePoints(grid, thisOctal, nPoints, points, connectivity, nUniquePoints, includeGhosts)
      use octal_mod, only: returndPhi
#ifdef MPI
      use mpi
#endif
      use inputs_mod, only : hydrodynamics

      type(OCTAL), pointer :: thisOctal, child
      integer(bigint) :: connectivity(:), nUniquePoints
      type(GRIDTYPE) :: grid
      real :: points(:,:)
      integer(bigint) :: nPoints
      integer :: subcell, i
      real :: xp, xm, yp, ym, zm, zp, d, r1, r2, phi, dphi
      real :: phi1, phi2, phiStart, phiEnd
      integer :: npointsInCell
      real :: pointCell(3, 8)
      integer :: nPhi, iPhi
      logical :: includeGhosts

      type(VECTOR) :: rVec
      do subcell = 1, thisOctal%maxChildren
         if (thisOctal%hasChild(subcell)) then
            ! find the child
            do i = 1, thisOctal%nChildren, 1
               if (thisOctal%indexChild(i) == subcell) then
                  child => thisOctal%child(i)
                  call recursiveGetUniquePoints(grid, child, nPoints, points, connectivity, nUniquePoints, includeGhosts)
                  exit
               end if
            end do
         else

#ifdef MPI
            if (.not.octalOnThread(thisOctal, subcell, myRankGlobal) .and. grid%splitOverMPI) cycle
            if (hydrodynamics) then
               if (checkGhost(thisOctal, subcell).and.(.not.includeGhosts)) cycle
            endif
#endif
            if (thisOctal%threed) then
               if (.not.thisOctal%cylindrical) then
                  rVec = subcellCentre(thisOctal,subcell)
                  d = real(thisOctal%subcellSize/2.d0)
                  
                  nPointsInCell = 8

                  xp = REAL(rVec%x + d)
                  xm = REAL(rVec%x - d)
                  yp = REAL(rVec%y + d)
                  ym = REAL(rVec%y - d)
                  zp = REAL(rVec%z + d)
                  zm = REAL(rVec%z - d)
                  
                  pointCell(1, 1) = xm
                  pointCell(2, 1) = ym
                  pointCell(3, 1) = zm

                  pointCell(1, 2) = xp
                  pointCell(2, 2) = ym
                  pointCell(3, 2) = zm

                  pointCell(1, 3) = xm
                  pointCell(2, 3) = yp
                  pointCell(3, 3) = zm

                  pointCell(1, 4) = xp
                  pointCell(2, 4) = yp
                  pointCell(3, 4) = zm

                  pointCell(1, 5) = xm
                  pointCell(2, 5) = ym
                  pointCell(3, 5) = zp

                  pointCell(1, 6) = xp
                  pointCell(2, 6) = ym
                  pointCell(3, 6) = zp

                  pointCell(1, 7) = xm
                  pointCell(2, 7) = yp
                  pointCell(3, 7) = zp

                  pointCell(1, 8) = xp
                  pointCell(2, 8) = yp
                  pointCell(3, 8) = zp

                  call addPoints(points, nUniquePoints, connectivity, nPoints, nPointsInCell, pointCell)
!                  if (myrankGlobal==1)write(*,*) "nPoints ",nPoints, " Unique Points ", nUniquePoints

               else
                  rVec = subcellCentre(thisOctal, subcell)
                  d = real(thisOctal%subcellSize/2.d0)
                  zp = REAL(rVec%z + d)
                  zm = REAL(rVec%z - d)
                  r1 = real(sqrt(rVec%x**2 + rVec%y**2) - d)
                  r2 = real(sqrt(rVec%x**2 + rVec%y**2) + d)
                  phi = real(atan2(rVec%y, rVec%x))
                  dphi = real(returndPhi(thisOctal))
                  phi1 = phi - dphi
                  phi2 = phi + dphi
                  nPhi = max(1, nint(dphi/(10.d0*degtorad)))
                  do iPhi = 1, nPhi
                     phiStart = real(phi1 + (phi2 - phi1) * dble(iphi-1)/dble(nphi))
                     phiEnd   = real(phi1 + (phi2 - phi1) * dble(iphi)/dble(nphi))


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
               d = real(thisOctal%subcellSize/2.d0)
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


    end subroutine recursiveGetUniquePoints
  end subroutine getUniquePoints

  subroutine addPoints(points, nUniquePoints, connectivity, nPoints, nPointsInCell, pointCell)
    real :: points(:,:)
    integer(bigint) :: nUniquePoints
    integer(bigint) :: connectivity(:)
    integer(bigint) :: nPoints, i, j, one
    real :: pointCell(:,:)
    integer :: nPointsInCell
    logical :: found
    one = 1
    do i = 1, nPointsInCell
       found = .false.
       do j = max(one,nUniquePoints-100), nUniquePoints
          if (theSamePoint(points(1:3,j), pointCell(1:3,i))) then
             connectivity(nPoints+i) = j-1
             found = .true.
             exit
          endif
       enddo
       if (.not.found) then
          nUniquePoints = nUniquePoints + 1
          points(1:3,nUniquePoints) = pointCell(1:3,i)
          connectivity(nPoints+i) = nUniquePoints-1
       endif
    enddo
    npoints = npoints + nPointsincell
  end subroutine addPoints

  logical function theSamePoint(p1, p2)
    real :: p1(3), p2(3)
    real, parameter :: tol = epsilon(p1)

    theSamePoint = .true.
    if (abs(p1(1)-p2(1)) > tol) then
       theSamePoint = .false.
    else if (abs(p1(2)-p2(2)) > tol) then
       theSamePoint = .false.
    else if (abs(p1(3)-p2(3)) > tol) then
       theSamePoint = .false.
    endif
  end function theSamePoint
    
  logical function checkGhost(thisOctal, subcell)
    use inputs_mod, only : cylindricalHydro
    type(VECTOR) :: rVec
    type(OCTAL), pointer :: thisOctal
    integer :: subcell

    checkGhost = .false.
    
    if (thisOctal%ghostCell(subcell)) checkGhost = .true.


    if (cylindricalHydro) then
       rVec = subcellCentre(thisOctal, subcell)
       if (rVec%x > 0.d0) checkGhost = .false.
    endif
  end function checkGhost

end module vtk_mod

