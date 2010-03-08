module setupamr_mod

  use amr_mod
  use vtk_mod
  use vector_mod
  use messages_mod
  USE constants_mod
  USE octal_mod, only: OCTAL, wrapperArray, octalWrapper, subcellCentre, cellVolume, &
       allocateattribute, copyattribute, deallocateattribute
  use gridtype_mod, only:   gridtype
  USE parallel_mod, ONLY:   torus_abort


  implicit none

contains

  subroutine setupFogel(grid, filename, speciesName)
    use input_variables, only : rinner, rOuter
    type(GRIDTYPE) :: grid
    character(len=*) :: filename, speciesName
    integer, parameter :: maxR = 200, maxZ = 200
    integer :: nr
    integer, allocatable :: nz(:)
    real(double), allocatable :: r(:), z(:,:), rho(:,:), t(:,:)
    real(double), allocatable :: abundance(:,:)
    integer :: iSpecies
    character(len=80) :: cJunk
    real(double) :: junk(5)
    real(double), allocatable :: rev(:)
    integer :: i, j
    logical :: gridConverged

    rinner = 5.*autocm/1.d10
    router = 400.*autocm/1.d10
    allocate(r(maxR), nz(maxR), z(maxR, maxZ), rho(maxR, maxZ), &
         t(maxR,maxZ), abundance(maxR,maxZ))

    select case (speciesName)
       case("CN")
          iSpecies = 1
       case("HCN")
          iSpecies = 2
       case("C")
          iSpecies = 3
       case("C+")
          iSpecies = 4
       case("H2CO")
          iSpecies = 5
       case DEFAULT
          call writeFatal("setupFogel: unknown species")
    end select
    

    open(20, file=filename, status="old", form="formatted")
    read(20,'(A)') cJunk
    read(20,'(A)') cJunk
    
    nr = 1
    nz(1) = 0

10  continue
    read(20,'(a)',end=30) cJunk
    if (len(trim(cJunk)) < 10) then ! blank line
       read(20,'(a)',end=30) cJunk
       read(20,'(a)',end=30) cJunk
       nr = nr + 1
       nz(nr) = 0
       goto 10
    endif

    nz(nr) = nz(nr) + 1
    read(cJunk,*) r(nr), z(nr,nz(nr)), rho(nr,nz(nr)), t(nr,nz(nr)), junk(1:5)
    abundance(nr, nz(nr)) = junk(iSpecies)
    goto 10
30  continue
    close(20)
    nr = nr - 1
    r = r * autocm/1.d10
    z = z * autocm/1.d10

    do i = 1, nr
       allocate(rev(1:nz(i)))
       do j = 1, nz(i)
          rev(j) = z(i,nz(i) - j + 1)
       enddo
       z(i,1:nz(i)) = rev(1:nz(i))

       do j = 1, nz(i)
          rev(j) = rho(i,nz(i) - j + 1)
       enddo
       rho(i,1:nz(i)) = rev(1:nz(i))

       do j = 1, nz(i)
          rev(j) = t(i,nz(i) - j + 1)
       enddo
       t(i,1:nz(i)) = rev(1:nz(i))

       do j = 1, nz(i)
          rev(j) = abundance(i,nz(i) - j + 1)
       enddo
       abundance(i,1:nz(i)) = rev(1:nz(i))
       deallocate(rev)
    enddo



    call splitGridFogel(grid%octreeRoot, grid, r, z, nr, nz, rho, t, abundance)
    do
       gridConverged = .true.
       call myScaleSmooth(3., grid, &
            gridConverged,  inheritProps = .true., interpProps = .false.)
       if (gridConverged) exit
    end do
    call writeWarning("Need to check whether abundances are really relative to N(H) rather than N(H_2)")
    call fillGridFogel(grid%octreeRoot, grid, r, z, nr, nz, rho, t, abundance)

    call writeVtkFile(grid, "fogel.vtk",  valueTypeString=(/"rho         ",&
         "temperature ","molabundance"/))

  end subroutine setupFogel

  recursive subroutine splitGridFogel(thisOctal, grid, r, z, nr, nz, rho, t, abundance)
    use input_variables, only : minDepthAMR, maxDepthAMR
    type(GRIDTYPE) :: grid
    type(OCTAL) :: thisOctal
    integer :: iIndex
    real(double) :: r(:), z(:,:), rho(:,:), t(:,:), abundance(:,:)
    integer :: nr, nz(:)
    logical :: splitInAzimuth
    logical :: split
    integer :: iSubcell, i,j,k
    type(VECTOR) :: rVec
    real(double) :: thisR, thisZ,s
    logical :: outsideGrid

    splitInAzimuth = .false.


    DO iSubcell = 1, thisOctal%maxChildren

       split = .false.
       rVec = subcellCentre(thisOctal, iSubcell)
       s = thisOctal%subcellSize/2.d0
       thisR = sqrt(rVec%x**2 + rVec%y**2)
       thisZ = abs(rVec%z)
       outsideGrid = .false.
       if ((thisR+s < r(1)).or.(thisR-s > r(nr))) outsideGrid = .true.
       if (.not.outsideGrid) then
          call locate(r, nr, thisR, i)
          if (thisZ-s > z(i,nz(i))) outsidegrid = .true.
       endif
       if (.not.OutsideGrid) then
          call locate(z(i,1:nz(i)), nz(i), thisZ, j)
          if (thisOctal%subcellSize > (r(i+1)-r(i))) split = .true.
          if (thisOctal%subcellSize > (z(i,j+1)-z(i,j))) then
             split = .true.
          endif
       endif
       if (thisOctal%nDepth >= maxDepthAMR) split = .false.
       if (thisOctal%nDepth <= minDepthAMR) split = .true.
       IF (split) then

          CALL addNewChild(thisOctal, iSubcell, grid, adjustGridInfo=.TRUE., &
               inherit = .true., &
               splitAzimuthally=splitInAzimuth)

          if (.not.thisOctal%hasChild(isubcell)) then
             write(*,*) "add child failed in splitGrid"
             do
             enddo
          endif



       END IF

    END DO

    do i = 1, thisOctal%nChildren
       if (.not.thisOctal%hasChild(thisOctal%indexchild(i))) then
          write(*,*) "octal children messed up"
          do ; enddo
          endif
       enddo

       do i = 1, thisOctal%maxChildren
          k = -99
          if (thisOctal%hasChild(i)) then
             do j = 1, thisOctal%nChildren
                if (thisOctal%indexChild(j) == i) then
                   k = j
                   exit
                endif
             enddo
             if (k==-99) then
                write(*,*) "subcell screwup"
                do
                enddo
             endif
          endif
       enddo

       if (any(thisOctal%haschild(1:thisOctal%maxChildren)).and.(thisOctal%nChildren==0)) then
          write(*,*) "nchildren screw up"
          do;enddo
          endif

          DO iIndex = 1, thisOctal%nChildren

             call  splitGridFogel(thisOctal%child(iIndex), grid, r, z, nr, nz, rho, t, abundance)

          END DO

        end subroutine splitGridFogel

        recursive subroutine fillGridFogel(thisOctal, grid, r, z, nr, nz, rho, t, abundance)
          use input_variables, only : mcore
          type(GRIDTYPE) :: grid
          type(octal), pointer   :: thisOctal
          type(octal), pointer  :: child 
          type(VECTOR) :: rVec
          real(double) :: thisR, thisZ,fac1,fac2,fac3
          integer :: subcell, i, j, k1,k2
          logical :: outsideGrid
          real(double) :: r(:), z(:,:), rho(:,:), t(:,:), abundance(:,:)
          integer :: nr, nz(:)

          do subcell = 1, thisOctal%maxChildren
             if (thisOctal%hasChild(subcell)) then
                ! find the child
                do i = 1, thisOctal%nChildren, 1
                   if (thisOctal%indexChild(i) == subcell) then
                      child => thisOctal%child(i)
                      call fillGridFogel(child, grid, r, z, nr, nz, rho, t, abundance)
                      exit
                   end if
                end do
             else

                rVec = subcellCentre(thisOctal, subcell)
                thisR = sqrt(rVec%x**2 + rVec%y**2)
                thisZ = abs(rVec%z)

                outsideGrid = .false.

                if ((thisR < r(1)).or.(thisR > r(nr))) outsideGrid = .true.
                if (.not.outsideGrid) then
                   call locate(r, nr, thisR, j)
                   if ((thisZ < z(j,1)).or.(thisZ > z(j,nz(j)))) outsideGrid = .true.
                endif

                thisOctal%rho(subcell) = 1.d-25
                thisOctal%temperature(subcell) = 3.d0
                thisOctal%molAbundance(subcell) = 1.d-20
                thisOctal%velocity(subcell) = VECTOR(0.d0, 0.d0, 0.d0)
                if (.not.outsideGrid) then
                   call locate(z(j,1:nz(j)), nz(j), thisZ, k1)
                   call locate(z(j+1,1:nz(j+1)), nz(j+1), thisZ, k2)
                   fac1 = (thisR - r(j))/(r(j+1)-r(j))
                   fac2 = (thisZ - z(j,k1))/(z(j,k1+1)-z(j,k1))
                   fac3 = (thisZ - z(j+1,k2))/(z(j+1,k2+1)-z(j+1,k2))

                   thisOctal%rho(subcell) = (1.d0-fac1)*(1.d0-fac2)* rho(j,k1) + &
                        (     fac1)*(1.d0-fac3)* rho(j+1,k2) + &
                        (1.d0-fac1)*(     fac2)* rho(j,k1+1) + &
                        (     fac1)*(     fac3)* rho(j+1,k2+1)  


                   thisOctal%temperature(subcell) = (1.d0-fac1)*(1.d0-fac2)* t(j,k1) + &
                        (     fac1)*(1.d0-fac3)* t(j+1,k2) + &
                        (1.d0-fac1)*(     fac2)* t(j,k1+1) + &
                        (     fac1)*(     fac3)* t(j+1,k2+1)  


                   thisOctal%molAbundance(subcell) = (1.d0-fac1)*(1.d0-fac2)* abundance(j,k1) + &
                        (     fac1)*(1.d0-fac3)* abundance(j+1,k2) + &
                        (1.d0-fac1)*(     fac2)* abundance(j,k1+1) + &
                        (     fac1)*(     fac3)* abundance(j+1,k2+1)  

                   thisOctal%microturb(subcell) = max(1d-8,sqrt((2.d-10 * kerg * thisOctal%temperature(subcell) &
                        / (28.0 * amu)) + 0.3**2) &
                        / (cspeed * 1e-5)) ! mu is 0.3km/s subsonic turbulence
                   mcore = 0.5d0 * mSol
                   thisOctal%velocity(subcell) = keplerianVelocity(rvec, grid)
                   thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.d0 * mHydrogen) 

                   thisOctal%molAbundance(subcell) = thisOctal%molAbundance(subcell) * 2.d0
                   if (subcell == 1) CALL fillVelocityCorners(thisOctal,grid,keplerianVelocity,thisOctal%threed)
                endif
                thisOctal%nh2(subcell) = thisOctal%rho(subcell)/(2.d0 * mHydrogen) 

             endif
          enddo
        end subroutine fillGridFogel



end module setupamr_mod
