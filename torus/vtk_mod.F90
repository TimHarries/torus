module vtk_mod
!
!

! written by tjh
! module to write VTK format files - see
! http://public.kitware.com/VTK/pdf/file-formats.pdf

  use kind_mod
  use ion_mod
  use constants_mod
  use utils_mod
  use amr_mod
  use mpi_amr_mod
  use messages_mod
  use vector_mod

  implicit none

  interface writeVtkFile
     module procedure writeVtkFileAMR
     module procedure writeVtkFileSource
  end interface

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
      use octal_mod, only: returndPhi

      type(OCTAL), pointer :: thisOctal, child

#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE) :: grid
      integer :: lunit
      integer :: subcell, i
      real :: xp, xm, yp, ym, zm, zp, d, r1, r2, phi, dphi
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
    integer :: nCells, nPoints, iOffset, nCount
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


  subroutine writeValue(grid, vtkFilename, valueType)
    type(GRIDTYPE) :: grid
    integer :: lunit = 69
    character(len=*) :: valueType
    character(len=*) :: vtkFilename
    logical :: vectorvalue, scalarvalue


    select case (valueType)
       case("velocity","hydrovelocity","linearvelocity","quadvelocity", "cornervel")
          scalarvalue = .false.
          vectorvalue = .true.
       case DEFAULT
          scalarvalue = .true.
          vectorvalue = .false.
    end select

    open(lunit, file=vtkFilename, form="formatted", status="old", position="append")
    if (writeHeader) then
       if (scalarvalue) then
          write(69,'(a,a,a)') "SCALARS ",trim(valueType)," float"
          write(69, '(a)') "LOOKUP_TABLE default"
       endif
       if (vectorvalue) then
          write(69, '(a,a,a)') "VECTORS ", trim(valueType), " float"
       endif

    endif

    call recursiveWriteValue(grid%octreeRoot, valueType, grid)
    close(lunit)

  contains

    recursive subroutine recursiveWriteValue(thisOctal, valueType, grid)
      use input_variables, only : lambdasmooth
      type(OCTAL), pointer :: thisOctal, child
#ifdef MPI
      include 'mpif.h'  
#endif
      type(GRIDTYPE), intent(in) :: grid
      type(VECTOR) :: rVec, vel
      integer :: lunit = 69
      integer :: subcell, i
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
               case("velocity")
                     ! stop vectors from showing up in visit if too big
                     if(thisoctal%velocity(subcell)%x .ge. 1.d0) then
                        thisoctal%velocity(subcell)%x = 0.d0
                        thisoctal%velocity(subcell)%y = 0.d0
                        thisoctal%velocity(subcell)%z = 0.d0
                     endif
                     write(lunit, *) real(thisOctal%velocity(subcell)%x*cspeed/1.e5), &
                          real(thisOctal%velocity(subcell)%y*cspeed/1.e5), &
                          real(thisOctal%velocity(subcell)%z*cspeed/1.e5)
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



               case("inflow")
                  if (thisOctal%inflow(subcell)) then
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
                  write(lunit, *) max( real(thisOctal%microturb(subcell)), min_single_prec )

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
       nPoints = source(iSource)%surface%nElements
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

          v2 = arbitraryRotate(cVec, -dtheta/2.d0, aVec)
          v2 = rotateZ(v2, -dphi/2.d0)

          v3 = arbitraryRotate(cVec, dtheta/2.d0, aVec)
          v3 = rotateZ(v3, dphi/2.d0)

          v4 = arbitraryRotate(cVec, dtheta/2.d0, aVec)
          v4 = rotateZ(v4, -dphi/2.d0)


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

  subroutine writeVtkFileAMR(grid, vtkFilename, valueTypeFilename, valueTypeString)
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
    integer :: nOctals, nVoxels, i, iType
    integer :: nPointOffset
#ifdef MPI
    integer :: ierr
    integer, allocatable :: iOffsetArray(:)
    integer :: myRank, nThreads, iThread
#endif

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
       write(lunit,'(a)') "FIELD FieldData 2"
       write(lunit,'(a)') "TIME 1 1 double"
       write(lunit,*) grid%currentTime
       write(lunit,'(a)') "CYCLE 1 1 int"
       write(lunit,*) grid%idump
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
    
       if (.not.grid%splitOverMPI) call writeValue(grid, vtkFilename, valueType(iType))

#ifdef MPI
    if (grid%splitOverMpi) then

       do iThread = 1, nThreadsGlobal
          call MPI_BARRIER(amrCOMMUNICATOR, ierr)
          if (iThread == myRankGlobal) then
             call writeValue(grid, vtkFilename, valueType(iType))
          endif
       enddo
    endif
#endif

    enddo

#ifdef MPI
666 continue
#endif

  end subroutine writeVtkFileAMR




end module vtk_mod

